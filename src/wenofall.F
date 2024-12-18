c  
c
c
c  Fallout with 5th-order WENO scheme (Flux form)
c
c  fw(ix,kz) = -4.0*dtdz4*vt(ix,kz+1) 
c  (flux at w-point through the bottom of the box uses the Vt at the 
c   scalar point)
c
      subroutine wenofall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,vt,
     :  a,db1,imapz,mzdist,ia,id,xfall,fs,loop,rks)
      
c
c Performs 6th-order Crowley, monotonically limited fallout on an x-z slice of
c  data (or a column if nx=1 on input)
c
c  nx,ny,nz : model dimensions
c  nor      : buffer zone around grid (usually nor=3 in SAM)
c  gt       : grid transformation array (stretch factors)
c  ngt      : number of variables in the gt array (ngt=10 in SAM)
c  dtp      : model time step
c  jgs      : y-index of the slice
c  vt       : x-z array containing the terminal speeds of the thing being advected
c  a        : 3-D array of the thing being advected (only do work on the one slice)
c  db1      : 2-D density array
c  imapz    : index of gt array for z-direction map (stretch) factors
c  ia       : scalar index of advected variable ( a(:::ia) )
c  id       : flag: id=1 use density (e.g. for mixing ratio), 
c                   id=0 no density (charge density or number concentration)
c  xfall    : array for precip/charge fallout values at the ground
c

      implicit none
      
      integer nx,ny,nz,nor,ngt,jgs,imapz,mzdist,na,ia
      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)
c      real gz(nz),z1d(nz,4)
      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
c      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
      real fs(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
c      real db(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
      real db1(nx,nz)
      real vt(nx,nz)
      real dtp,dz,dtz1
      real qtmp1(-nor:nx+nor,-nor:nz+nor)
      real qtmp2(-nor:nx+nor,-nor:nz+nor)
      real cmax
      real xfall(nx,ny,na)
      real t0(-nor:nx+nor,-nor:nz+nor)
      integer loop,rks
      

      real  qim2, qim1, qi, qip1, qip2, vv, dir
      real*8  beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
      real, parameter    :: f30 =  7./12., f31 = 1./12
      real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-15
      integer, parameter :: pw = 2

      integer  i, k, im1, jm1, km1
      
      
      integer    np
      parameter (np=6)
      real pz(nz,np,np)

      integer ix,jy,kz,ndfall,n
      integer kd1,id1,kstag,istag
      parameter (kd1=1,id1=1,kstag=1,istag=1)
      
      real dtdz4,dt,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,dtfac
      integer kz1
      
      dtfac = 1.0
      
      jy = jgs
      dt = dtp
      
      DO kz = -nor,nz+nor
        DO ix = -nor,nx+nor
         t0(ix,kz) = 0.0
         qtmp1(ix,kz) = 0.0
        ENDDO
      ENDDO

c      write(6,*) 'dz,ngt,imapz,jgs,ia,id',dz,ngt,imapz,jgs,ia,id
c      write(6,*) 'ny,nor,na ',ny,nor,na 
      
      cmax = 0.0
      DO kz = 1,nz-1
        DO ix = 1,nx-1
         cmax = Max(cmax, vt(ix,kz)*gz(kz))
         qtmp1(ix,kz) = a(ix,jgs,kz,ia)
         qtmp2(ix,kz) = 0.0
        ENDDO
      ENDDO
      
c first check if fallout is worth doing
      IF ( cmax .eq. 0.0 ) THEN
        RETURN
      ENDIF

c      write(6,*) 'jgs,ia,id,cmax',jgs,ia,id,cmax

      
c
c Check if we need to subdivide the time step for stability purposes
c      
      
      IF ( dtp*cmax .lt. 0.98 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*cmax) + 1)
c
c DONT SUBDIVIDE STEP HERE!  Must do it in a loop that calls this 
c  subroutine so that fall speeds get updated!
c        dt = dtp/Float(ndfall)
c        write(6,*) 'jgs, cmax*dtp,ndfall = ',jgs,cmax*dtp,ndfall
      ENDIF
      

      DO n = 1,ndfall


      DO k = 1,     nz-kstag-kd1
       km1 = max(k-1,1)
        DO i = 1,     nx-istag

         vv  = Min(0.0, -vt(i,k)) ! 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))

         IF(  k .ge. 4 .and. k .le. nz-3 ) THEN

 ! assume dir = -1.0
 
            qip2 = qtmp1(i,k-2)
            qip1 = qtmp1(i,k-1)
            qi   = qtmp1(i,k  )
            qim1 = qtmp1(i,k+1)
            qim2 = qtmp1(i,k+2)
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw
         wi1 = gi1 / (eps + beta1)**pw
         wi2 = gi2 / (eps + beta2)**pw
    
         sumwk = wi0 + wi1 + wi2
    
         t0(i,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

c         ELSEIF( k .ge. 3 .and. k .le. nz-2 ) THEN

c         t0(i,k) = vv * ( f30 * (qtmp1(i,k  ) + qtmp1(i,k-1))  
c     :                  - f31 * (qtmp1(i,k+1) + qtmp1(i,k-2))  
c     :                  - f31 * (qtmp1(i,k+1) - qtmp1(i,k-2) - 3.*(qtmp1(i,k)-qtmp1(i,k-1))) )

c         ELSEIF ( k .le. 2 .or. k .eq. nz-1 ) THEN

c         t0(i,k) = vv * 0.5 * (qtmp1(i,k) + qtmp1(i,km1))

         ELSE  !( k .eq. 1  ) THEN

         t0(i,k) = -vt(i,k) * (qtmp1(i,k) )

         ENDIF
     
        ENDDO
      ENDDO




      

      do 4033 ix = 1,     nx-id1
c      t0(ix, 0) = 0.5*(vt(ix,1)+vt(ix,2))*qtmp1(ix,1)*dt/dz
      t0(ix, 0) =  -vt(ix,1)*qtmp1(ix,1)
      t0(ix,nz-kstag) = 0.0
      t0(ix,nz) = 0.0
 4033 continue
 
      DO kz = 1, nz-kstag-kd1
       DO ix = 1,     nx-id1

c without zero divergence

       IF ( id .eq. 1 ) THEN
       qtmp1(ix,kz) =  qtmp1(ix,kz)
     :   + (t0(ix,kz)*db1(ix,kz) 
     :    - t0(ix,kz+kd1)*db1(ix,Min(nz-1,kz+1)))
     >  *gz(kz)/db1(ix,kz)    
     
c      >  dtz1*(vt(ix,kz+1)*db(ix,jy,kz+1)*qtmp1(ix,kz+1)
c     >       -vt(ix,kz)*db(ix,jy,kz)*qtmp1(ix,kz))
c     < /db(ix,jy,kz)
     
       
       ELSE
       qtmp1(ix,kz) =  qtmp1(ix,kz)
     :   + (t0(ix,kz)
     :   - t0(ix,kz+kd1))
     >  *gz(kz)
       
       ENDIF

       IF ( kz .eq. 1 .and. loop .eq. rks ) THEN
         xfall(ix,jy,ia) = xfall(ix,jy,ia)  
     :       - t0(ix, 0)*dt*gz(kz)
       ENDIF

      ENDDO
      ENDDO
c
c
c
 
      ENDDO ! n

      do kz = 1,nz-1-1 
      do ix = 1,nx-1
        fs(ix,jgs,kz,ia) =   qtmp1(ix,kz) - a(ix,jgs,kz,ia)
c     :   + (t0(ix,kz)*db1(ix,kz) 
c     :    - t0(ix,kz+kd1)*db1(ix,Min(nz-1,kz+1)))
c     >  *gz(kz)/db1(ix,kz)    
      enddo
      enddo


      
      RETURN
      END

c
c2345678901234567890123456789012345678901234567890123456789012345678912
c
c
