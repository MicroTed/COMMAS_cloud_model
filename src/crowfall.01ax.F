#include "sam.def.h"
c  
c
c
c  Fallout with 6th-order crowley scheme (Flux form)
c
c  version 01a: fw(ix,kz) = -4.0*dtdz4*vt(ix,kz+1) 
c  (flux at w-point through the bottom of the box uses the Vt at the 
c   scalar point)
c
#ifdef SWM
      subroutine crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,vt,
     :  a,db1,imapz,mzdist,ia,id,xfall)
#else
      subroutine crowfall(nx,ny,nz,nor,na,gt,ngt,dtp,dz,jgs,vt,
     :  a,db1,imapz,mzdist,ia,id,xfall)
#endif

#ifdef MPI
      USE COMMASMPI_MODULE
#endif
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
#ifdef SWM
      parameter(ng1 = 1)
      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
c      real a(nx,ny,nz,na)
#else
      real gt(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,ngt) ! grid transformation
      parameter(ng1 = 0)
#endif
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
c      real db(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
      real db1(nx,nz+1)
      real vt(nx,nz+1)
      real dtp,dz,dtz1
      real qtmp1(-nor:nx+nor,-nor:nz+nor)
      real qtmp2(-nor:nx+nor,-nor:nz+nor)
      real cmax
      real xfall(nx,ny,na)
      real t0(-nor:nx+nor,-nor:nz+nor)
      real t1(-nor:nx+nor,-nor:nz+nor)
      real t2(-nor:nx+nor,-nor:nz+nor)
      real t3(-nor:nx+nor,-nor:nz+nor)
      real t4(-nor:nx+nor,-nor:nz+nor)
      real t5(-nor:nx+nor,-nor:nz+nor)
      real t6(-nor:nx+nor,-nor:nz+nor)
      real t7(-nor:nx+nor,-nor:nz+nor)
      
      real fw(-nor:nx+nor,-nor:nz+nor)
      
      
      integer    np
      parameter (np=6)
c      real pz(nz,np,np)
      real pz(-nor+1:nz+nor,np,np)
      
      integer ix,jy,kz,ndfall,n
      integer kd1,id1,kstag,istag
      parameter (kd1=1,id1=1,kstag=1,istag=1)
      
      real dtdz4,dt,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,dtfac
      integer kz1

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.
#endif

      dtfac = 1.0
      
      jy = jgs
      dt = dtp

#ifdef MPI
      kzb = -nor
      kze = ktile+nor

      ixb = -nor
      ixe = itile+nor

      DO kz = kzb,kze
        DO ix = ixb,ixe
#else
      DO kz = -nor,nz+nor
        DO ix = -nor,nx+nor
#endif
         t0(ix,kz) = 0.0
         qtmp1(ix,kz) = 0.0
        ENDDO
      ENDDO

c      write(6,*) 'dz,ngt,imapz,jgs,ia,id',dz,ngt,imapz,jgs,ia,id
c      write(6,*) 'ny,nor,na ',ny,nor,na 
      
      cmax = 0.0
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO kz = kzb,kze
        DO ix = ixb,ixe
#else
      DO kz = 1,nz-1
        DO ix = 1,nx-1
#endif
#ifdef SWM
         cmax = Max(cmax, vt(ix,kz)*gz(kz))
#else
         cmax = Max(cmax,
     :     vt(ix,kz)/dz*gt(ix,jgs,kz,imapz))
#endif
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
c  create z-direction polynomial coefficients
c
c      call crwcff(np,nz,nz,kstag,kd1,pz)
      call crwcff(np,nz,nz,nor,kstag,kd1,nz,1,nz,1,nz,pz)
c
c
      
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


      dtdz4 =  (0.25)*dt/dz

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do  kz = 1,     nz-kstag-kd1
      do  ix = 1,     nx-istag
#endif
        fw(ix,kz) = -4.0*dtdz4*vt(ix,kz+1)
      ENDDO
      ENDDO

c
c  compute polynomials

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
c C$DOACROSS LOCAL(kz,ix)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,ix)
      do  kz = 1,     nz-kstag-kd1
      do  ix = 1,     nx-istag
#endif
      t1(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,1,1) +
     >  fw(ix,kz)*(pz(kz,1,2) +
     >  fw(ix,kz)*(pz(kz,1,3) +
     >  fw(ix,kz)*(pz(kz,1,4) +
     >  fw(ix,kz)*(pz(kz,1,5) +
     >  fw(ix,kz)*(pz(kz,1,6) ))))))
      t2(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,2,1) +
     >  fw(ix,kz)*(pz(kz,2,2) +
     >  fw(ix,kz)*(pz(kz,2,3) +
     >  fw(ix,kz)*(pz(kz,2,4) +
     >  fw(ix,kz)*(pz(kz,2,5) +
     >  fw(ix,kz)*(pz(kz,2,6) ))))))

      t3(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,3,1) +
     >  fw(ix,kz)*(pz(kz,3,2) +
     >  fw(ix,kz)*(pz(kz,3,3) +
     >  fw(ix,kz)*(pz(kz,3,4) +
     >  fw(ix,kz)*(pz(kz,3,5) +
     >  fw(ix,kz)*(pz(kz,3,6) ))))))

      t4(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,4,1) +
     >  fw(ix,kz)*(pz(kz,4,2) +
     >  fw(ix,kz)*(pz(kz,4,3) +
     >  fw(ix,kz)*(pz(kz,4,4) +
     >  fw(ix,kz)*(pz(kz,4,5) +
     >  fw(ix,kz)*(pz(kz,4,6) ))))))

      t5(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,5,1) +
     >  fw(ix,kz)*(pz(kz,5,2) +
     >  fw(ix,kz)*(pz(kz,5,3) +
     >  fw(ix,kz)*(pz(kz,5,4) +
     >  fw(ix,kz)*(pz(kz,5,5) +
     >  fw(ix,kz)*(pz(kz,5,6) ))))))

      t6(ix,kz) = 
     >  fw(ix,kz)*(pz(kz,6,1) +
     >  fw(ix,kz)*(pz(kz,6,2) +
     >  fw(ix,kz)*(pz(kz,6,3) +
     >  fw(ix,kz)*(pz(kz,6,4) +
     >  fw(ix,kz)*(pz(kz,6,5) +
     >  fw(ix,kz)*(pz(kz,6,6) ))))))


      ENDDO
      ENDDO

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do  kz = 1,     nz-kstag-kd1
      do  ix = 1,     nx-istag
#endif
      t0(ix,kz) = 
     >  t1(ix,kz)*qtmp1(ix,kz-2)+
     >  t2(ix,kz)*qtmp1(ix,kz-1)+
     >  t3(ix,kz)*qtmp1(ix,kz  )+
     >  t4(ix,kz)*qtmp1(ix,kz+1)+
     >  t5(ix,kz)*qtmp1(ix,kz+2)+
     >  t6(ix,kz)*qtmp1(ix,kz+3)
      ENDDO
      ENDDO

c
c  monotonic face values
c
c      if ( icrwmn .eq. 1 ) then 
      
#ifdef MPI
      kzb = 0
      kze = nor

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do  kz = 0,     nor
      do  ix = 1,     nx-istag
#endif
      qtmp1(ix,-kz)   = qtmp1(ix,kz+1)
      qtmp1(ix,nz+kz) = qtmp1(ix,nz-kstag-kz)
      ENDDO
      ENDDO
c
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag-kd1
      do ix = 1,     nx-istag ! -id1
#endif
      t0(ix,kz) = -t0(ix,kz) / (fw(ix,kz) + 1.0e-20)
      qsr  = max(sign((1.0),fw(ix,kz)), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*qtmp1(ix,kz-kd1*1)
     >     + qsl*qtmp1(ix,kz+kd1*2)
      qd   = qsr*qtmp1(ix,kz+kd1*1)
     >     + qsl*qtmp1(ix,kz      )
      qc   = qsr*qtmp1(ix,kz      )
     >     + qsl*qtmp1(ix,kz+kd1*1)
      qr   = qu + (qc - qu) / (abs(fw(ix,kz)) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),
     >           abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
      t0(ix,kz) =
     > -fw(ix,kz) *
     > (
     >   ((1.0)-qxp) * qc
     > +        qxp
     > * (      qxs  * min(max(t0(ix,kz), qc), min(qr, qd))
     > + ((1.0)-qxs) * max(min(t0(ix,kz), qc), max(qr, qd)) )
     < )
      enddo !! ix
      enddo !! kz
      

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do ix = ixb,ixe
#else
      do ix = 1,     nx-id1
#endif
c      t0(ix, 0) = 0.5*(vt(ix,1)+vt(ix,2))*qtmp1(ix,1)*dt/dz
      t0(ix, 0) =  vt(ix,1)*qtmp1(ix,1)*dt/dz
      t0(ix,nz-kstag) = 0.0
      fw(ix, 0) = 0.0
      fw(ix,nz-kstag) = 0.0
      enddo

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      DO kz = kzb,kze
       DO ix = ixb,ixe
#else
      DO kz = 1, nz-kstag-kd1
       DO ix = 1,     nx-id1
#endif
c without zero divergence

       IF ( id .eq. 1 ) THEN
       qtmp1(ix,kz) =  qtmp1(ix,kz)
     :   + (t0(ix,kz)*db1(ix,kz) 
     :   - t0(ix,kz-kd1)*db1(ix,Max(1,kz-kd1)))
#ifdef SWM
     >  *dz*gz(kz)/db1(ix,kz)    
#else
     >  * (gt(ix,jy,kz,imapz))/db1(ix,kz)    
#endif
     
c      >  dtz1*(vt(ix,kz+1)*db(ix,jy,kz+1)*qtmp1(ix,kz+1)
c     >       -vt(ix,kz)*db(ix,jy,kz)*qtmp1(ix,kz))
c     < /db(ix,jy,kz)
     
       
       ELSE
       qtmp1(ix,kz) =  qtmp1(ix,kz)
     :   + (t0(ix,kz)
     :   - t0(ix,kz-kd1))
#ifdef SWM
     >  *dz*gz(kz)
#else
     >  * (gt(ix,jy,kz,imapz))
#endif
       
       ENDIF

#ifdef MPI
       IF ( kzbeg .eq. 1) THEN
#endif
       IF ( kz .eq. 1 ) THEN
         xfall(ix,jy,ia) = xfall(ix,jy,ia)  
#ifdef SWM
     :       + t0(ix, 0)*dz*gz(kz)
#else
     :       + t0(ix, 0)*gt(ix,jy,kz,imapz)
#endif
       ENDIF
#ifdef MPI
       ENDIF
#endif
      ENDDO
      ENDDO
c
c
c
 
      ENDDO ! n

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1,nz-1-1 
      do ix = 1,nx-1
#endif
        a(ix,jgs,kz,ia) =  qtmp1(ix,kz) 
      enddo
      enddo


      
      RETURN
      END

c
c2345678901234567890123456789012345678901234567890123456789012345678912
c
c
