c
c     ##################################################################
c     ######                                                      ######
c     ######       SUBROUTINE STRZSOR  (localsor2 w/ vert. strz)  ######
c     ######                                                      ######
c     ##################################################################
      SUBROUTINE STRZSOR(mask,rhs,phi,nx,ny,nz,nzt,i0,j0,k0,mdel,
     :    cofx,cofy,cofz, bcx, bcy,
     :     ix1,jy1,kz1,ix2,jy2,kz2,itype,rjac1,omega1,ibc,eps,ndebug)
      
c#######################################################################
c  Variable declarations      
c  Will set ibc = 0 for now, regardless of input value
c
c 11/9/2004: Fixed a bug in the sweep over the even z planes:
c            Had k=jst,kfn,2 instead of k=kst,kfn,2.  OOPS!  This bug
c            has probably been there since going to vertical stretching.
c            (Also had a 'jst' instead of 'ist' in another spot)
c
c  2/6/2003 Fixed a couple bugs: the IF statement setting 'ibw' had
c            'ibe' in the ELSE part by mistake (should have been 'ibw')
c           Also, since 'nz' is now the top of the extended domain, added
c           the input variable 'nzt', which is the top of the dynamics
c           domain and is now used in the check for closeness to the top
c           boundary.
c
c  spring 2002 vertically stretched version
c
c  8/11/99  When a local domain gets close to a side or the top, 
c           move the local boundary out to the full domain boundary
c           and use Neumann conditions. 
c  8/15/99  Also when the local domain gets close to the ground, 
c           automatically extend domain.
c#######################################################################
      implicit none
      
      integer ndebug
c      parameter (ndebug=1)
      integer  maxits,nx,ny,nz,nzt  ! nzt is top of dynamics domain

      real cofx(nx,3),cofy(ny,3),cofz(0:nz,3)
      
      integer      ::  bcx, bcy

      integer ibc  ! 1 for Neumann, 0 for Dirichlet on top and sides,
                         !  (D on bottom)
                   ! 2 for periodic on sides with Dirichlet top and bottom
      real rjac,rhs(nx,ny,0:nz),phi(-1:nx+2,-1:ny+2,0:nz+1), eps
      real rjac1
      real mask(-1:nx+2,-1:ny+2,nz)
c      real tem(nx,ny,nz)
      parameter (maxits=200)
      integer minits
      integer nout
      parameter (nout = 20)
      integer i0,j0,k0         ! center of local box (the lower-left 
                               !         endpoint of the new bond)
      integer mdel             ! 'radius' of local box
      integer idel             ! temporary box 'radius'
      integer i1,j1,k1,i2,j2,k2  ! limits of relaxation
      integer ix1,jy1,kz1,ix2,jy2,kz2  ! limits of relaxation or initial box
      integer itype              ! flag to 
                                 !   0: start with small box around i0,j0,k0
                                 !   1: start with box around flash region:
                                 !       i.e. (i1,j1,k1,i2,j2,k2)
      
      integer ist,ifn,jst,jfn,kst,kfn

      integer nzmax
      parameter(nzmax=1000)
      integer ipass,k,ksw,i,isw,n,j,jsw(nzmax,2)
      real residmax(nzmax)
      real anormc,anormf,resid,err1 ! ,omega
      real omega,omega1
      real err(nzmax),anorm(nzmax)
      real anormh1,anormh2,anorm1
      real frac
      parameter(frac = 1.0/6.0)
      integer accel
      parameter (accel = 0)
c      parameter (omega1 = 1.9)  ! optimum for 169x169x83 grid, initial relaxation
c      parameter (omega1 = 1.45)
c      parameter (rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 ) )
      real remax1,phim,phiminv
      
      integer ibe,ibw,jbs,jbn,kbt ! type of boundary on east/west/north/south/top
                                  !  boundary of box.  
      integer mh   !  threshhold distance (in gridpoints) from a domain
                   !   boundary for extending local box directly to the boundary 
      parameter (mh = 10)            
      
      integer   :: nxa, nyc, nze, nper
c #################################################      
c Begin Executable Code      
c #################################################      
      IF (nz.gt.nzmax) THEN
       write(6,*) '(nz.gt.nzmax) !!!! make nzmax bigger!',nz,nzmax
       STOP
      END IF
      rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 )
      ksw = 1
      DO i=1,2
        isw = ksw     
         DO k=1,nz
           jsw(k,i) = isw
           isw = 3 - isw
         END DO ! k
        ksw = 3 - ksw
      END DO ! i
      
c Boundary condition flags
      nxa = 1
      nyc = 1
      nze = 1
      nper = 1
      IF ( bcx .eq. 2 ) THEN 
         nxa = 0
         nper = 0
      ENDIF
      IF ( bcy .eq. 2 ) THEN
         nyc = 0
         nper = 0
      ENDIF
      
c      minits = (mdel + 1)/2
       minits = 1
c      ibc = 0
c
c determine boundaries of max box so we can find phim
c
      IF (itype .eq. 0) THEN
       i1 = Max(2, i0 - mdel)
       i2 = Min(nx-1, i0 + mdel + 1)
       j1 = Max(2, j0 - mdel)
       j2 = Min(ny-1, j0 + mdel + 1)
       k1 = Max(1, k0 - mdel)
       k2 = Min(nz-1, k0 + mdel + 1)
      ELSE
       i1 = Max(2, ix1 - mdel)
       i2 = Min(nx-1, ix2 + mdel + 1)
       j1 = Max(2, jy1 - mdel) 
       j2 = Min(ny-1, jy2 + mdel + 1)
       k1 = Max(1, kz1 - mdel)
       k2 = Min(nz-1, kz2 + mdel + 1)
      END IF
      
         
c setup initial norm.  Stop when anorm < eps*anormf
      anormf = 0.00
      phim = 0.0
      DO k=k1,k2
      DO j=j1,j2
      DO i=i1,i2
        phim = Max( phim,Abs(phi(i,j,k)) )
        anormf = anormf + Abs(rhs(i,j,k))
      END DO
      END DO
      END DO
      IF (ndebug .ge. 1) write(6,'(a,e13.6)') 'anormf = ',anormf
      IF (ndebug .ge. 1) write(6,'(a,e13.6)') 'phim = ',phim
      phiminv = 1.0
      IF (phim .gt. 0.0) phiminv = 1.0/phim
      
c start SOR iterations
      IF (accel .eq. 1) omega = 1.00
      
      IF (ndebug .ge. 1) write(6,*) 'start relax'
c        IF (ndebug .ge. 1 ) 
c     :     write(6,'(a)')
c     :   'n,anormc,anormc/anormf,remax1,remax1*phiminv'
c     :   'iteration,anormc/anormf,err1'

c
c SOR LOOP
c
      DO n=1,maxits
      
c current box boundaries:
      IF (itype .eq. 0) THEN
       idel = Min(2*n + mdel/2,mdel)
       i1 = Max(2, i0 - idel)
       i2 = Min(nx-1, i0 + idel + 1)
       j1 = Max(2, j0 - idel)
       j2 = Min(ny-1, j0 + idel + 1)
       k1 = Max(2, k0 - idel)
       k2 = Min(nz-1, k0 + idel + 1)
      ELSE
       idel = Min(2*n,nz)
       i1 = Max(2, ix1 - idel)
       i2 = Min(nx-1, ix2 + idel + 1)
       j1 = Max(2, jy1 - idel) 
       j2 = Min(ny-1, jy2 + idel + 1)
       k1 = Max(2, kz1 - idel)
       k2 = Min(nz-1, kz2 + idel + 1)
      END IF
      

c check if close to a boundary (other than ground) and if so,
c  then extend the local box to the boundary and set flag to
c  use Neumann B.C.s
      IF ( i1 .le. mh ) THEN 
         ibw = 1
c         i1 = 2
         i1 = 1
         IF ( bcx .eq. 2 ) i2 = nx - 1
      ELSE 
         ibw = 0
      ENDIF
      IF ( i2 .ge. nx-mh ) THEN 
         ibe = 1
         i2 = nx-1
         IF ( bcx .eq. 2 ) i1 = 1
      ELSE 
         ibe = 0
      ENDIF
      IF ( j1 .le. mh ) THEN 
         jbs = 1
c         j1 = 2
         j1 = 1
         IF ( bcy .eq. 2 ) j2 = ny - 1
      ELSE 
         jbs = 0
      ENDIF
      IF ( j2 .ge. ny-mh ) THEN 
         jbn = 1
         j2 = ny-1
         IF ( bcy .eq. 2 ) j1 = 1
      ELSE 
         jbn = 0
      ENDIF
      IF ( k2 .ge. nzt-10 ) THEN 
         kbt = 1
         k2 = nz-1
      ELSE 
         kbt = 0
      ENDIF
      IF ( k1 .le. 5 ) THEN 
         k1 = 1
      ENDIF
        
c set up Neumann boundary conditions for du/dn = 0.0
c on top and side edges.  Bottom edge assumed to be Dirichlet
c and set before calling SOR.
c  top/bottom
       IF (ibc .eq. 1 .or. kbt .eq. 1 ) THEN
        DO j=j1,j2
        DO i=i1,i2
          phi(i,j,nz) = phi(i,j,nz-1)
c          phi(i,j,1) = 0.0
        END DO
        END DO
       END IF
c east/west
       IF (ibc .eq. 1 .or. ibw .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
c          phi(1,j,k) = phi(2,j,k)
          phi(0,j,k) = phi(1,j,k)
        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. ibe .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
          phi(nx,j,k) = phi(nx-1,j,k)
        END DO
        END DO
       END IF
c north/south
       IF (ibc .eq. 1 .or. jbs .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
c          phi(i,1,k) = phi(i,2,k)
          phi(i,0,k) = phi(i,1,k)
        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. jbn .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
          phi(i,ny,k) = phi(i,ny-1,k)
        END DO
        END DO
       END IF
        
c       ELSEIF (ibc .eq. 2) THEN ! periodic on sides
c east/west
c        DO k=k1,k2
c        DO j=j1,j2
c          phi(1,j,k) = phi(nx-1,j,k)
c          phi(nx,j,k) = phi(2,j,k)
c        END DO
c        END DO
c north/south
c        DO k=k1,k2
c        DO i=i1,i2
c          phi(i,1,k) = phi(i,ny-1,k)
c          phi(i,ny,k) = phi(i,2,k)
c        END DO
c        END DO
c        
c       END IF

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
c        DO k=Max(1,k1-1),k2+1
cc         anorm(k) = 0.0
cc         err(k)   = 0.0
c         residmax(k) = 0.0
c        END DO
        remax1 = 0.0
        
         anormc = 0.00
         err1 = 0.00
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = i1
      ifn = i2
      jst = j1
      jfn = j2
c      if (nyd.eq.1) jfn = ny-1
      kst = k1
c      if (nze.eq.1) kst = 3
      kfn = k2

      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'n,i1,j1,k1,i2,j2,k2 = ',n,i1,j1,k1,i2,j2,k2,nz,nzt
      ENDIF
        
c      if (nzf.eq.1) kfn = nz-1
c
c     set periodic boundaries if necessary
c
c
c   red (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) or even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   red(x,y) on even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      remax1 = 0.0
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask),
!$OMP+PRIVATE(i,j,k,resid),
C$OMP+REDUCTION (max : remax1)
      do k=kst,kfn
      do j=jst,jfn
        do i=ist,ifn
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid =  rhs(i,j,k)-(
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1) +
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))*phi(i,j,k)  )
          remax1 = Max(remax1,Abs(resid))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
        IF (ndebug .ge. 1 ) ! .and. (n/nout)*nout.eq.n) 
     :     write(0,'(a,i5,4(3x,e13.6))') 'n,anormc...',
     :     n,anormc,anormf,remax1,remax1*phiminv ! ,
!     :     write(6,'(i5,4(3x,e13.6))') 
!     :     n,anormc,anormc/anormf,remax1,remax1*phiminv ! ,
c     :     anormh1,anormh1/anormf
       
        IF (remax1*phiminv .lt. eps .and. n .ge. minits ) GOTO 999

      END DO ! (n=1,maxits)
       
 999  CONTINUE
c       write(6,'(a,i7)') 'number of iterations = ',n
c       write(6,'(i5,4(3x,e13.6))') 
c     :     n,anormc,anormc/anormf,remax1,remax1*phiminv 
      RETURN
      END


c#######################################################################
c#######################################################################
c#######################################################################


c     ##################################################################
c     ######                                                      ######
c     ######     SUBROUTINE STRZSOR2D  (localsor2 w/ vert. strz)  ######
c     ######                                                      ######
c     ##################################################################
      SUBROUTINE STRZSOR2D(mask,rhs,phi,nx,ny,nz,nzt,i0,j0,k0,mdel,
     :    cofx,cofy,cofz, bcx, bcy,
     :     ix1,jy1,kz1,ix2,jy2,kz2,itype,rjac1,omega1,ibc,eps,ndebug)
      
c#######################################################################
c  Variable declarations      
c  Will set ibc = 0 for now, regardless of input value
c
c 11/9/2004: Fixed a bug in the sweep over the even z planes:
c            Had k=jst,kfn,2 instead of k=kst,kfn,2.  OOPS!  This bug
c            has probably been there since going to vertical stretching.
c            (Also had a 'jst' instead of 'ist' in another spot)
c
c  2/6/2003 Fixed a couple bugs: the IF statement setting 'ibw' had
c            'ibe' in the ELSE part by mistake (should have been 'ibw')
c           Also, since 'nz' is now the top of the extended domain, added
c           the input variable 'nzt', which is the top of the dynamics
c           domain and is now used in the check for closeness to the top
c           boundary.
c
c  spring 2002 vertically stretched version
c
c  8/11/99  When a local domain gets close to a side or the top, 
c           move the local boundary out to the full domain boundary
c           and use Neumann conditions. 
c  8/15/99  Also when the local domain gets close to the ground, 
c           automatically extend domain.
c#######################################################################
      implicit none
      
      integer ndebug
c      parameter (ndebug=1)
      integer  maxits,nx,ny,nz,nzt  ! nzt is top of dynamics domain

      real cofx(nx,3),cofy(ny,3),cofz(0:nz,3)
      
      integer      ::  bcx, bcy

      integer ibc  ! 1 for Neumann, 0 for Dirichlet on top and sides,
                         !  (D on bottom)
                   ! 2 for periodic on sides with Dirichlet top and bottom
      real rjac,rhs(nx,ny,0:nz),phi(-1:nx+2,-1:ny+2,0:nz+1), eps
      real rjac1
      real mask(-1:nx+2,-1:ny+2,nz)
c      real tem(nx,ny,nz)
      parameter (maxits=200)
      integer minits
      integer nout
      parameter (nout = 20)
      integer i0,j0,k0         ! center of local box (the lower-left 
                               !         endpoint of the new bond)
      integer mdel             ! 'radius' of local box
      integer idel             ! temporary box 'radius'
      integer i1,j1,k1,i2,j2,k2  ! limits of relaxation
      integer ix1,jy1,kz1,ix2,jy2,kz2  ! limits of relaxation or initial box
      integer itype              ! flag to 
                                 !   0: start with small box around i0,j0,k0
                                 !   1: start with box around flash region:
                                 !       i.e. (i1,j1,k1,i2,j2,k2)
      
      integer ist,ifn,jst,jfn,kst,kfn

      integer nzmax
      parameter(nzmax=4000)
      integer ipass,k,ksw,i,isw,n,j,jsw(nzmax,2)
      real residmax(nzmax)
      real anormc,anormf,err1 ! ,omega
      real omega,omega1
      real err(nzmax),anorm(nzmax)
      real anormh1,anormh2,anorm1
      real frac
      parameter(frac = 1.0/6.0)
      integer accel
      parameter (accel = 0)
c      parameter (omega1 = 1.9)  ! optimum for 169x169x83 grid, initial relaxation
c      parameter (omega1 = 1.45)
c      parameter (rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 ) )
      real phim,phiminv
      double precision :: resid,remax1
      
      integer ibe,ibw,jbs,jbn,kbt ! type of boundary on east/west/north/south/top
                                  !  boundary of box.  
      integer mh   !  threshhold distance (in gridpoints) from a domain
                   !   boundary for extending local box directly to the boundary 
      parameter (mh = 10)            
      
      integer   :: nxa, nyc, nze, nper
      
      logical redblack
      
c #################################################      
c Begin Executable Code      
c #################################################      
      IF (nz.gt.nzmax) THEN
       write(6,*) '(nz.gt.nzmax) make nzmax bigger!',nz,nzmax
       STOP
      END IF
      rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 )
      ksw = 1
      DO i=1,2
        isw = ksw     
         DO k=1,nz
           jsw(k,i) = isw
           isw = 3 - isw
         END DO ! k
        ksw = 3 - ksw
      END DO ! i
      
c Boundary condition flags
      nxa = 1
      nyc = 1
      nze = 1
      nper = 1
      IF ( bcx .eq. 2 ) THEN 
         nxa = 0
         nper = 0
      ENDIF
      IF ( bcy .eq. 2 ) THEN
         nyc = 0
         nper = 0
      ENDIF
      
c      minits = (mdel + 1)/2
       minits = 2
c      ibc = 0
c
c determine boundaries of max box so we can find phim
c
c      IF (itype .eq. 0) THEN
      IF ( .false. ) THEN
       i1 = Max(2, i0 - mdel)
       i2 = Min(nx-1, i0 + mdel + 1)
       j1 = Max(2, j0 - mdel)
       j2 = Min(ny-1, j0 + mdel + 1)
       k1 = Max(1, k0 - mdel)
       k2 = Min(nz-1, k0 + mdel + 1)
      ELSE
       i1 = Max(2, ix1 - mdel)
       i2 = Min(nx-1, ix2 + mdel + 1)
       j1 = Max(2, jy1 - mdel) 
       j2 = Min(ny-1, jy2 + mdel + 1)
       k1 = Max(1, kz1 - mdel)
       k2 = Min(nz-1, kz2 + mdel + 1)
      END IF
      
      IF ( ny .le. 2 ) THEN ! 2D case
        j1 = 1
        j2 = 1
      ENDIF
      
         
c setup initial norm.  Stop when anorm < eps*anormf
      anormf = 0.00
      phim = 0.0
      DO k=k1,k2
      DO j=j1,j2
      DO i=i1,i2
        phim = Max( phim,Abs(phi(i,j,k)) )
        anormf = anormf + Abs(rhs(i,j,k))
      END DO
      END DO
      END DO
      IF (ndebug .ge. 1) write(0,'(a,e13.6)') 'anormf = ',anormf
      IF (ndebug .ge. 1) write(0,'(a,e13.6)') 'phim = ',phim
      IF ( anormf .le. 0.0 ) THEN
        write(0,*) 'strzsor2d: problem with anormf = ',anormf, 'phim = ',phim
        write(0,*) 'strzsor2d: i1,i2,j1,j2,k1,k2 = ',i1,i2,j1,j2,k1,k2
        write(0,*) 'ny = ',ny
        stop
      ENDIF
      phiminv = 1.0
      IF (phim .gt. 0.0) phiminv = 1.0/phim
      
c start SOR iterations
      IF (accel .eq. 1) omega = 1.00
      
      IF (ndebug .ge. 1) write(6,*) 'start relax'
c        IF (ndebug .ge. 1 ) 
c     :     write(6,'(a)')
c     :   'n,anormc,anormc/anormf,remax1,remax1*phiminv'
c     :   'iteration,anormc/anormf,err1'

c
c SOR LOOP
c
      DO n=1,maxits
      
c current box boundaries:
c      IF (itype .eq. 0) THEN
      IF ( .false. ) THEN
       idel = Min(2*n + mdel/2,mdel)
       i1 = Max(2, i0 - idel)
       i2 = Min(nx-1, i0 + idel + 1)
       j1 = Max(2, j0 - idel)
       j2 = Min(ny-1, j0 + idel + 1)
       k1 = Max(2, k0 - idel)
       k2 = Min(nz-1, k0 + idel + 1)
      ELSE
       idel = Min(2*n,nz)
       i1 = Max(2, ix1 - idel)
       i2 = Min(nx-1, ix2 + idel + 1)
       j1 = Max(2, jy1 - idel) 
       j2 = Min(ny-1, jy2 + idel + 1)
       k1 = Max(2, kz1 - idel)
       k2 = Min(nz-1, kz2 + idel + 1)
      END IF
      
      IF ( ny .le. 2 ) THEN ! 2D case
        j1 = 1
        j2 = 1
      ENDIF

c check if close to a boundary (other than ground) and if so,
c  then extend the local box to the boundary and set flag to
c  use Neumann B.C.s
      IF ( i1 .le. mh ) THEN 
         ibw = 1
c         i1 = 2
         i1 = 1
         IF ( bcx .eq. 2 ) i2 = nx - 1
      ELSE 
         ibw = 0
      ENDIF
      IF ( i2 .ge. nx-mh ) THEN 
         ibe = 1
         i2 = nx-1
         IF ( bcx .eq. 2 ) i1 = 1
      ELSE 
         ibe = 0
      ENDIF
      IF ( j1 .le. mh ) THEN 
         jbs = 1
c         j1 = 2
         j1 = 1
         IF ( bcy .eq. 2 ) j2 = ny - 1
      ELSE 
         jbs = 0
      ENDIF
      IF ( j2 .ge. ny-mh ) THEN 
         jbn = 1
         j2 = 1 ! set to 1 for 2D case ! ny-1
         IF ( bcy .eq. 2 ) j1 = 1
      ELSE 
         jbn = 0
      ENDIF
      IF ( k2 .ge. nzt-10 ) THEN 
         kbt = 1
         k2 = nz-1
      ELSE 
         kbt = 0
      ENDIF
      IF ( k1 .le. 5 ) THEN 
         k1 = 1
      ENDIF
        
c set up Neumann boundary conditions for du/dn = 0.0
c on top and side edges.  Bottom edge assumed to be Dirichlet
c and set before calling SOR.
c  top/bottom
       IF (ibc .eq. 1 .or. kbt .eq. 1 ) THEN
!        DO j=j1,j2
        j = 1
        DO i=i1,i2
          phi(i,j,nz) = phi(i,j,nz-1)
c          phi(i,j,1) = 0.0
        END DO
!        END DO
       END IF
c east/west
       IF (ibc .eq. 1 .or. ibw .eq. 1 ) THEN
        DO k=k1,k2
!        DO j=j1,Max(1,j2)
c          phi(1,j,k) = phi(2,j,k)
          j=1
          phi(0,j,k) = phi(1,j,k)
!        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. ibe .eq. 1 ) THEN
        DO k=k1,k2
!        DO j=j1,Max(1,j2)
          j = 1
          phi(nx,j,k) = phi(nx-1,j,k)
!        END DO
        END DO
       END IF
c north/south
       IF (ibc .eq. 1 .or. jbs .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
c          phi(i,1,k) = phi(i,2,k)
          phi(i,0,k) = phi(i,1,k)
        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. jbn .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
          phi(i,ny,k) = phi(i,ny-1,k)
        END DO
        END DO
       END IF
        
c       ELSEIF (ibc .eq. 2) THEN ! periodic on sides
c east/west
c        DO k=k1,k2
c        DO j=j1,j2
c          phi(1,j,k) = phi(nx-1,j,k)
c          phi(nx,j,k) = phi(2,j,k)
c        END DO
c        END DO
c north/south
c        DO k=k1,k2
c        DO i=i1,i2
c          phi(i,1,k) = phi(i,ny-1,k)
c          phi(i,ny,k) = phi(i,2,k)
c        END DO
c        END DO
c        
c       END IF

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
c        DO k=Max(1,k1-1),k2+1
cc         anorm(k) = 0.0
cc         err(k)   = 0.0
c         residmax(k) = 0.0
c        END DO
        remax1 = 0.0
        
         anormc = 0.00
         err1 = 0.00
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = i1
      ifn = i2
      jst = j1
      jfn = j2
c      if (nyd.eq.1) jfn = ny-1
      kst = k1
c      if (nze.eq.1) kst = 3
      kfn = k2

      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'strzsor2d: n,i1,j1,k1,i2,j2,k2 = ',n,i1,j1,k1,i2,j2,k2,nz,nzt
       write(6,*) 'strzsor2d: ibc,ibw = ',ibc,ibw

      do k=kst,kfn
      do i=ist,ifn
        j=jst
          IF ( .not. (-1.e30 < phi(i,j,k) .and. phi(i,j,k) < 1.e30 ) ) THEN
            write(0,*) 'Problem in strzsor2d o! phi = ',phi(i,j,k),i,k
!            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
!            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k)
            STOP
          ENDIF
      end do
      end do

      ENDIF

      
      
      redblack = .true.
      
      
      IF ( .not. redblack ) THEN ! gauss-seidel
      
      do k=kst,kfn
      do i=ist,ifn
        j=jst
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
c     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofz(k,3))
         ENDIF
      end do
      end do
      
      ELSE
        
c      if (nzf.eq.1) kfn = nz-1
c
c     set periodic boundaries if necessary
c
c
c   red (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist,ifn,2
        j=jst
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid = phi(i,j,k)
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
c     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofz(k,3))
          IF ( ndebug >= 1 ) THEN
          IF ( .not. (-1.e30 < phi(i,j,k) .and. phi(i,j,k) < 1.e30 ) ) THEN
            write(0,*) 'Problem in strzsor2d 1! phi = ',phi(i,j,k),i,k
            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k), resid
            write(0,*) cofx(i,1),cofx(i,2),cofx(i,3)
            write(0,*) cofz(k,1),cofz(k,2),cofz(k,3)
            STOP
          ENDIF
          ENDIF
         ENDIF
      end do
      end do

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) or even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist+1,ifn,2
        j=jst
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid = phi(i,j,k)
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
c     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofz(k,3))
          IF ( ndebug >= 1 ) THEN
          IF ( .not. (-1.e30 < phi(i,j,k) .and. phi(i,j,k) < 1.e30 ) ) THEN
            write(0,*) 'Problem in strzsor2d 2! phi = ',phi(i,j,k),i,k
            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k), resid
            write(0,*) cofx(i,1),cofx(i,2),cofx(i,3)
            write(0,*) cofz(k,1),cofz(k,2),cofz(k,3)
            STOP
          ENDIF
          ENDIF
         ENDIF
      end do
      end do
c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist+1,ifn,2
        j=jst
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid = phi(i,j,k)
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
c     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofz(k,3))
          IF ( ndebug >= 1 ) THEN
          IF ( .not. (-1.e30 < phi(i,j,k) .and. phi(i,j,k) < 1.e30 ) ) THEN
            write(0,*) 'Problem in strzsor2d 3! phi = ',phi(i,j,k),i,k
            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k), resid
            write(0,*) cofx(i,1),cofx(i,2),cofx(i,3)
            write(0,*) cofz(k,1),cofz(k,2),cofz(k,3)
            STOP
          ENDIF
          ENDIF
         ENDIF
      end do
      end do

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   red(x,y) on even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist,ifn,2
         j=jst
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid = phi(i,j,k)
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
c     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofz(k,3))
          IF ( ndebug >= 1 ) THEN
          IF ( .not. (-1.e30 < phi(i,j,k) .and. phi(i,j,k) < 1.e30 ) ) THEN
            write(0,*) 'Problem in strzsor2d 4! phi = ',phi(i,j,k),i,k
            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k), resid
            write(0,*) cofx(i,1),cofx(i,2),cofx(i,3)
            write(0,*) cofz(k,1),cofz(k,2),cofz(k,3)
            STOP
          ENDIF
          ENDIF
         ENDIF
      end do
      end do
      
      ENDIF ! redblack

        IF (ndebug .ge. 1 ) THEN ! .and. (n/nout)*nout.eq.n) 
          write(0,*) 'strzsor2d: calculate residual'
        ENDIF
      remax1 = 0.0
! !$OMP PARALLEL DO 
! !$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask),
! !$OMP+PRIVATE(i,j,k,resid),
! C$OMP+REDUCTION (max : remax1)
      do k=kst,kfn
      j=jst
        do i=ist,ifn
         IF ( mask(i,j,k) .eq. 1.0 ) THEN
          resid =  rhs(i,j,k)-(
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1) +
     +      (cofx(i,3)+cofz(k,3))*phi(i,j,k)  )
          IF ( .not. (-1.e30 < resid .and. resid < 1.e30 ) ) THEN
            write(0,*) 'Problem ins strzsor2d! resid = ',resid,i,k
            write(0,*) phi(i-1,j,k  ), phi(i+1,j,k  )
            write(0,*) phi(i  ,j,k-1), phi(i  ,j,k+1)
            write(0,*) rhs(i,j,k), phi(i,j,k)
            write(0,*) cofx(i,1),cofx(i,2),cofx(i,3)
            write(0,*) cofz(k,1),cofz(k,2),cofz(k,3)
            STOP
          ENDIF
          remax1 = Max(remax1,Abs(resid))
         ENDIF
        end do
      end do

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
        IF (ndebug .ge. 1 ) THEN ! .and. (n/nout)*nout.eq.n) 
          write(0,*) 'strzsor2d: end of pass ',n
          write(6,'(i5,5(3x,e13.6))') 
     :     n,anormc,anormf,anormc/anormf,remax1,remax1*phiminv ! ,
     
        ENDIF
c     :     anormh1,anormh1/anormf
       
        
        IF (remax1*phiminv .lt. eps .and. n .ge. minits ) EXIT

      END DO ! (n=1,maxits)
       
 999  CONTINUE
c       write(6,'(a,i7)') 'number of iterations = ',n
c       write(6,'(i5,4(3x,e13.6))') 
c     :     n,anormc,anormc/anormf,remax1,remax1*phiminv 
      RETURN
      END

c#######################################################################
c#######################################################################
c#######################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######       SUBROUTINE STRZSORtst  (localsor2 w/ vert. strz)  ######
c     ######                                                      ######
c     ##################################################################
      SUBROUTINE STRZSORtst(mask,rhs,phi,nx,ny,nz,nzt,i0,j0,k0,mdel,
     :    cofx,cofy,cofz, bcx, bcy,
     :     ix1,jy1,kz1,ix2,jy2,kz2,itype,rjac1,omega1,ibc,eps,ndebug)
      
c#######################################################################
c  Variable declarations      
c  Will set ibc = 0 for now, regardless of input value
c
c 11/9/2004: Fixed a bug in the sweep over the even z planes:
c            Had k=jst,kfn,2 instead of k=kst,kfn,2.  OOPS!  This bug
c            has probably been there since going to vertical stretching.
c            (Also had a 'jst' instead of 'ist' in another spot)
c
c  2/6/2003 Fixed a couple bugs: the IF statement setting 'ibw' had
c            'ibe' in the ELSE part by mistake (should have been 'ibw')
c           Also, since 'nz' is now the top of the extended domain, added
c           the input variable 'nzt', which is the top of the dynamics
c           domain and is now used in the check for closeness to the top
c           boundary.
c
c  spring 2002 vertically stretched version
c
c  8/11/99  When a local domain gets close to a side or the top, 
c           move the local boundary out to the full domain boundary
c           and use Neumann conditions. 
c  8/15/99  Also when the local domain gets close to the ground, 
c           automatically extend domain.
c#######################################################################
      implicit none
      
      integer ndebug
c      parameter (ndebug=1)
      integer  maxits,nx,ny,nz,nzt  ! nzt is top of dynamics domain

      real cofx(nx,3),cofy(ny,3),cofz(0:nz,3)
      
      integer      ::  bcx, bcy

      integer ibc  ! 1 for Neumann, 0 for Dirichlet on top and sides,
                         !  (D on bottom)
                   ! 2 for periodic on sides with Dirichlet top and bottom
      real rjac,rhs(0:nx+1,0:ny+1,0:nz+1),phi(0:nx+1,0:ny+1,0:nz+1), eps
      real rjac1
      real mask(0:nx+1,0:ny+1,0:nz+1)
c      real tem(nx,ny,nz)
      parameter (maxits=200)
      integer minits
      integer nout
      parameter (nout = 20)
      integer i0,j0,k0         ! center of local box (the lower-left 
                               !         endpoint of the new bond)
      integer mdel             ! 'radius' of local box
      integer idel             ! temporary box 'radius'
      integer i1,j1,k1,i2,j2,k2  ! limits of relaxation
      integer ix1,jy1,kz1,ix2,jy2,kz2  ! limits of relaxation or initial box
      integer itype              ! flag to 
                                 !   0: start with small box around i0,j0,k0
                                 !   1: start with box around flash region:
                                 !       i.e. (i1,j1,k1,i2,j2,k2)
      
      integer ist,ifn,jst,jfn,kst,kfn

      integer nzmax
      parameter(nzmax=4000)
      integer ipass,k,ksw,i,isw,n,j,jsw(nzmax,2)
      real residmax(nzmax)
      real anormc,anormf,resid,err1 ! ,omega
      real omega,omega1
      real err(nzmax),anorm(nzmax)
      real anormh1,anormh2,anorm1
      real frac
      parameter(frac = 1.0/6.0)
      integer accel
      parameter (accel = 0)
c      parameter (omega1 = 1.9)  ! optimum for 169x169x83 grid, initial relaxation
c      parameter (omega1 = 1.45)
c      parameter (rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 ) )
      real remax1,phim,phiminv
      
      integer ibe,ibw,jbs,jbn,kbt ! type of boundary on east/west/north/south/top
                                  !  boundary of box.  
      integer mh   !  threshhold distance (in gridpoints) from a domain
                   !   boundary for extending local box directly to the boundary 
      parameter (mh = 10)            
      
      integer   :: nxa, nyc, nze, nper
c #################################################      
c Begin Executable Code      
c #################################################      
      IF (nz.gt.nzmax) THEN
       write(6,*) '(nz.gt.nzmax) !!!! make nzmax bigger!',nz,nzmax
       STOP
      END IF
      rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 )
      ksw = 1
      DO i=1,2
        isw = ksw     
         DO k=1,nz
           jsw(k,i) = isw
           isw = 3 - isw
         END DO ! k
        ksw = 3 - ksw
      END DO ! i
      
c Boundary condition flags
      nxa = 1
      nyc = 1
      nze = 1
      nper = 1
      IF ( bcx .eq. 2 ) THEN 
         nxa = 0
         nper = 0
      ENDIF
      IF ( bcy .eq. 2 ) THEN
         nyc = 0
         nper = 0
      ENDIF
      
c      minits = (mdel + 1)/2
       minits = 1
c      ibc = 0
c
c determine boundaries of max box so we can find phim
c
      IF (itype .eq. 0) THEN
       i1 = Max(2, i0 - mdel)
       i2 = Min(nx-1, i0 + mdel + 1)
       j1 = Max(2, j0 - mdel)
       j2 = Min(ny-1, j0 + mdel + 1)
       k1 = Max(1, k0 - mdel)
       k2 = Min(nz-1, k0 + mdel + 1)
      ELSE
       i1 = Max(2, ix1 - mdel)
       i2 = Min(nx-1, ix2 + mdel + 1)
       j1 = Max(2, jy1 - mdel) 
       j2 = Min(ny-1, jy2 + mdel + 1)
       k1 = Max(1, kz1 - mdel)
       k2 = Min(nz-1, kz2 + mdel + 1)
      END IF
      
         
c setup initial norm.  Stop when anorm < eps*anormf
      anormf = 0.00
      phim = 0.0
      DO k=k1,k2
      DO j=j1,j2
      DO i=i1,i2
        phim = Max( phim,Abs(phi(i,j,k)) )
        anormf = anormf + Abs(rhs(i,j,k))
      END DO
      END DO
      END DO
      IF (ndebug .ge. 1) write(6,'(a,e13.6)') 'anormf = ',anormf
      IF (ndebug .ge. 1) write(6,'(a,e13.6)') 'phim = ',phim
      phiminv = 1.0
      IF (phim .gt. 0.0) phiminv = 1.0/phim
      
c start SOR iterations
      IF (accel .eq. 1) omega = 1.00
      
      IF (ndebug .ge. 1) write(6,*) 'start relax'
c        IF (ndebug .ge. 1 ) 
c     :     write(6,'(a)')
c     :   'n,anormc,anormc/anormf,remax1,remax1*phiminv'
c     :   'iteration,anormc/anormf,err1'

c
c SOR LOOP
c
      DO n=1,maxits
      
c current box boundaries:
      IF (itype .eq. 0) THEN
       idel = Min(2*n + mdel/2,mdel)
       i1 = Max(2, i0 - idel)
       i2 = Min(nx-1, i0 + idel + 1)
       j1 = Max(2, j0 - idel)
       j2 = Min(ny-1, j0 + idel + 1)
       k1 = Max(2, k0 - idel)
       k2 = Min(nz-1, k0 + idel + 1)
      ELSE
       idel = Min(2*n,nz)
       i1 = Max(2, ix1 - idel)
       i2 = Min(nx-1, ix2 + idel + 1)
       j1 = Max(2, jy1 - idel) 
       j2 = Min(ny-1, jy2 + idel + 1)
       k1 = Max(2, kz1 - idel)
       k2 = Min(nz-1, kz2 + idel + 1)
      END IF
      

c check if close to a boundary (other than ground) and if so,
c  then extend the local box to the boundary and set flag to
c  use Neumann B.C.s
      IF ( i1 .le. mh ) THEN 
         ibw = 1
c         i1 = 2
         i1 = 1
         IF ( bcx .eq. 2 ) i2 = nx - 1
      ELSE 
         ibw = 0
      ENDIF
      IF ( i2 .ge. nx-mh ) THEN 
         ibe = 1
         i2 = nx-1
         IF ( bcx .eq. 2 ) i1 = 1
      ELSE 
         ibe = 0
      ENDIF
      IF ( j1 .le. mh ) THEN 
         jbs = 1
c         j1 = 2
         j1 = 1
         IF ( bcy .eq. 2 ) j2 = ny - 1
      ELSE 
         jbs = 0
      ENDIF
      IF ( j2 .ge. ny-mh ) THEN 
         jbn = 1
         j2 = ny-1
         IF ( bcy .eq. 2 ) j1 = 1
      ELSE 
         jbn = 0
      ENDIF
      IF ( k2 .ge. nzt-10 ) THEN 
         kbt = 1
         k2 = nz-1
      ELSE 
         kbt = 0
      ENDIF
      IF ( k1 .le. 5 ) THEN 
         k1 = 1
      ENDIF
        
c set up Neumann boundary conditions for du/dn = 0.0
c on top and side edges.  Bottom edge assumed to be Dirichlet
c and set before calling SOR.
c  top/bottom
       IF (ibc .eq. 1 .or. kbt .eq. 1 ) THEN
        DO j=j1,j2
        DO i=i1,i2
          phi(i,j,nz) = phi(i,j,nz-1)
c          phi(i,j,1) = 0.0
        END DO
        END DO
       END IF
c east/west
       IF (ibc .eq. 1 .or. ibw .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
c          phi(1,j,k) = phi(2,j,k)
          phi(0,j,k) = phi(1,j,k)
        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. ibe .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
          phi(nx,j,k) = phi(nx-1,j,k)
        END DO
        END DO
       END IF
c north/south
       IF (ibc .eq. 1 .or. jbs .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
c          phi(i,1,k) = phi(i,2,k)
          phi(i,0,k) = phi(i,1,k)
        END DO
        END DO
       END IF
       IF (ibc .eq. 1 .or. jbn .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
          phi(i,ny,k) = phi(i,ny-1,k)
        END DO
        END DO
       END IF
        
c       ELSEIF (ibc .eq. 2) THEN ! periodic on sides
c east/west
c        DO k=k1,k2
c        DO j=j1,j2
c          phi(1,j,k) = phi(nx-1,j,k)
c          phi(nx,j,k) = phi(2,j,k)
c        END DO
c        END DO
c north/south
c        DO k=k1,k2
c        DO i=i1,i2
c          phi(i,1,k) = phi(i,ny-1,k)
c          phi(i,ny,k) = phi(i,2,k)
c        END DO
c        END DO
c        
c       END IF

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
c        DO k=Max(1,k1-1),k2+1
cc         anorm(k) = 0.0
cc         err(k)   = 0.0
c         residmax(k) = 0.0
c        END DO
        remax1 = 0.0
        
         anormc = 0.00
         err1 = 0.00
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = i1
      ifn = i2
      jst = j1
      jfn = j2
c      if (nyd.eq.1) jfn = ny-1
      kst = k1
c      if (nze.eq.1) kst = 3
      kfn = k2

      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'n,i1,j1,k1,i2,j2,k2 = ',n,i1,j1,k1,i2,j2,k2,nz,nzt
      ENDIF
        
c      if (nzf.eq.1) kfn = nz-1
c
c     set periodic boundaries if necessary
c
c
c   red (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) or even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) on odd z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
      do i=ist,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   red(x,y) on even z planes
c
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask)
!$OMP+PRIVATE(i,j,k)
      do k=kst+1,kfn,2
      do i=ist,ifn,2
        do j=jst,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      do i=ist+1,ifn,2
        do j=jst+1,jfn,2
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
         ENDIF
        end do
      end do
      end do

      remax1 = 0.0
!$OMP PARALLEL DO 
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask),
!$OMP+PRIVATE(i,j,k,resid),
C$OMP+REDUCTION (max : remax1)
      do k=kst,kfn
      do j=jst,jfn
        do i=ist,ifn
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          resid =  rhs(i,j,k)-(
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1) +
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))*phi(i,j,k)  )
          remax1 = Max(remax1,Abs(resid))
         ENDIF
        end do
      end do
      end do

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
        IF (ndebug .ge. 1 ) ! .and. (n/nout)*nout.eq.n) 
     :     write(6,'(i5,4(3x,e13.6))') 
     :     n,anormc,anormc/anormf,remax1,remax1*phiminv ! ,
c     :     anormh1,anormh1/anormf
       
        
        IF (remax1*phiminv .lt. eps .and. n .ge. minits ) GOTO 999

      END DO ! (n=1,maxits)
       
 999  CONTINUE
c       write(6,'(a,i7)') 'number of iterations = ',n
c       write(6,'(i5,4(3x,e13.6))') 
c     :     n,anormc,anormc/anormf,remax1,remax1*phiminv 
      RETURN
      END


c#######################################################################
c#######################################################################
c#######################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######       SUBROUTINE STRZSORtst  (localsor2 w/ vert. strz)  ######
c     ######                                                      ######
c     ##################################################################
      SUBROUTINE STRZSORmpi(mask,rhs,phi,nx,ny,nz,nzt1,i0,j0,k0,mdel0,
     :    cofx,cofy,cofz, bcx, bcy,ibcn,ibcs,ibce,ibcw,
     :     ix1,jy1,kz1,ix2,jy2,kz2,itype,rjac1,omega1,ibc,eps,ndebug,
     :     numits)
      
       USE ELEC_MODULE
       USE COMMASMPI_MODULE


c#######################################################################
c  Variable declarations      
c  Will set ibc = 0 for now, regardless of input value
c
c 11/9/2004: Fixed a bug in the sweep over the even z planes:
c            Had k=jst,kfn,2 instead of k=kst,kfn,2.  OOPS!  This bug
c            has probably been there since going to vertical stretching.
c            (Also had a 'jst' instead of 'ist' in another spot)
c
c  2/6/2003 Fixed a couple bugs: the IF statement setting 'ibw' had
c            'ibe' in the ELSE part by mistake (should have been 'ibw')
c           Also, since 'nz' is now the top of the extended domain, added
c           the input variable 'nzt1', which is the top of the dynamics
c           domain and is now used in the check for closeness to the top
c           boundary.
c
c  spring 2002 vertically stretched version
c
c  8/11/99  When a local domain gets close to a side or the top, 
c           move the local boundary out to the full domain boundary
c           and use Neumann conditions. 
c  8/15/99  Also when the local domain gets close to the ground, 
c           automatically extend domain.
c#######################################################################
      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
#endif
      
      integer ndebug,numits
c      parameter (ndebug=1)
      integer  maxits,nx,ny,nz,nzt1  ! nzt1 is top of dynamics domain

      real cofx(nx,3),cofy(ny,3),cofz(0:nz,3)
      
      integer      ::  bcx, bcy

      integer ibc  ! 1 for Neumann, 0 for Dirichlet on top and sides,
                         !  (D on bottom)
                   ! 2 for periodic on sides with Dirichlet top and bottom
      integer ibcn,ibcs,ibce,ibcw
      real rjac,rhs(0:nx+1,0:ny+1,0:nz+1),phi(0:nx+1,0:ny+1,0:nz+1), eps
      real rjac1
      real mask(0:nx+1,0:ny+1,0:nz+1)
c      real tem(nx,ny,nz)
      parameter (maxits=200)
      integer minits
      integer nout
      parameter (nout = 20)
      integer i0,j0,k0         ! center of local box (the lower-left 
                               !         endpoint of the new bond)
      integer mdel0             ! 'radius' of local box
      integer idel             ! temporary box 'radius'
      integer i1,j1,k1,i2,j2,k2  ! limits of relaxation
      integer ix1,jy1,kz1,ix2,jy2,kz2  ! limits of relaxation or initial box
      integer itype              ! flag to 
                                 !   0: start with small box around i0,j0,k0
                                 !   1: start with box around flash region:
                                 !       i.e. (i1,j1,k1,i2,j2,k2)
      
      integer ist,ifn,jst,jfn,kst,kfn

      integer nzmax
      parameter(nzmax=1000)
      integer ipass,k,ksw,i,isw,n,j,jsw(nzmax,2)
      real residmax(nzmax)
      real anormc,anormf,resid,err1 ! ,omega
      real anormfred,anormcred
      real omega,omega1
      real err(nzmax),anorm(nzmax)
      real anormh1,anormh2,anorm1
      real frac
      parameter(frac = 1.0/6.0)
      integer accel
      parameter (accel = 0)
c      parameter (omega1 = 1.9)  ! optimum for 169x169x83 grid, initial relaxation
c      parameter (omega1 = 1.45)
c      parameter (rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 ) )
      real remax1,phim,phiminv
      real remax1red,phimred
      
      integer ibe,ibw,jbs,jbn,kbt ! type of boundary on east/west/north/south/top
                                  !  boundary of box.  
      integer mh   !  threshhold distance (in gridpoints) from a domain
                   !   boundary for extending local box directly to the boundary 
      parameter (mh = 10)            
      
      integer   :: nxa, nyc, nze, nper
      
      integer :: pts,ptstart, ptend, ptstride
      integer :: ibeg, iend
      
      integer ng

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      
c #################################################      
c Begin Executable Code      
c #################################################      
      IF (ndebug .ge. 1) write(0,*) 'starting STRZSORmpi',my_rank
      
      numits = 0
      
      IF (nz.gt.nzmax) THEN
       write(6,*) '(nz.gt.nzmax) !!!! make nzmax bigger!',nz,nzmax
       STOP
      END IF
      rjac = Sqrt( 1.0 - (2.0/omega1 - 1)**2 )
      ksw = 1
      DO i=1,2
        isw = ksw     
         DO k=1,nz
           jsw(k,i) = isw
           isw = 3 - isw
         END DO ! k
        ksw = 3 - ksw
      END DO ! i
      
c Boundary condition flags
      nxa = 1
      nyc = 1
      nze = 1
      nper = 1
      IF ( bcx .eq. 2 ) THEN 
         nxa = 0
         nper = 0
      ENDIF
      IF ( bcy .eq. 2 ) THEN
         nyc = 0
         nper = 0
      ENDIF
      
c      minits = (mdel + 1)/2
       minits = 1
c      ibc = 0
c
c determine boundaries of max box so we can find phim
c
      IF (itype .eq. 0) THEN
       i1 = Max(2, i0 - mdel)
       i2 = Min(nx-1, i0 + mdel + 1)
       j1 = Max(2, j0 - mdel)
       j2 = Min(ny-1, j0 + mdel + 1)
       k1 = Max(1, k0 - mdel)
       k2 = Min(nz-1, k0 + mdel + 1)
      ELSE
       i1 = Max(2, ix1 - mdel)
       i2 = Min(nx-1, ix2 + mdel + 1)
       j1 = Max(2, jy1 - mdel) 
       j2 = Min(ny-1, jy2 + mdel + 1)
       k1 = Max(1, kz1 - mdel)
       k2 = Min(nz-1, kz2 + mdel + 1)
      END IF
      
         
c setup initial norm.  Stop when anorm < eps*anormf
      anormf = 0.00
      phim = 0.0
c      DO k=k1,k2
c      DO j=j1,j2
c      DO i=i1,i2
      
      DO k=1, nz-1
      DO j=1, ny
      DO i=1, nx
        phim = Max( phim,Abs(phi(i,j,k)) )
        anormf = anormf + Abs(rhs(i,j,k))
      END DO
      END DO
      END DO

#ifdef MPI


       phimred = phim

       CALL MPI_AllReduce(phim, phimred, 1, MPI_REAL, 
     :                   MPI_MAX, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       phim = phimred


       anormfred = anormf

       CALL MPI_AllReduce(anormf, anormfred, 1, MPI_REAL, 
     :                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       anormf = anormfred

#endif

      IF (ndebug .ge. 1) write(0,'(a,e13.6)') 'anormf = ',anormf
      IF (ndebug .ge. 1) write(0,'(a,e13.6)') 'phim = ',phim
      phiminv = 1.0
      IF (phim .gt. 0.0) phiminv = 1.0/phim
      
c start SOR iterations
      IF (accel .eq. 1) omega = 1.00
      
      IF (ndebug .ge. 1) write(6,*) 'start relax'
c        IF (ndebug .ge. 1 ) 
c     :     write(6,'(a)')
c     :   'n,anormc,anormc/anormf,remax1,remax1*phiminv'
c     :   'iteration,anormc/anormf,err1'

c
c SOR LOOP
c
      DO n=1,maxits
      
c current box boundaries:
      IF (itype .eq. 0) THEN
       idel = Min(2*n + mdel/2,mdel)
       i1 = Max(2, i0 - idel)
       i2 = Min(nx-1, i0 + idel + 1)
       j1 = Max(2, j0 - idel)
       j2 = Min(ny-1, j0 + idel + 1)
       k1 = Max(2, k0 - idel)
       k2 = Min(nz-1, k0 + idel + 1)
      ELSE
       idel = Min(2*n,nz)
       i1 = Max(2, ix1 - idel)
       i2 = Min(nx-1, ix2 + idel + 1)
       j1 = Max(2, jy1 - idel) 
       j2 = Min(ny-1, jy2 + idel + 1)
       k1 = Max(2, kz1 - idel)
       k2 = Min(nz-1, kz2 + idel + 1)
      END IF
      

c check if close to a boundary (other than ground) and if so,
c  then extend the local box to the boundary and set flag to
c  use Neumann B.C.s
      IF ( i1 .le. mh ) THEN 
         ibw = 1
c         i1 = 2
         i1 = 1
         IF ( bcx .eq. 2 ) i2 = nx - 1
      ELSE 
         ibw = 0
      ENDIF
      IF ( i2 .ge. nx-mh ) THEN 
         ibe = 1
         i2 = nx-1
         IF ( bcx .eq. 2 ) i1 = 1
      ELSE 
         ibe = 0
      ENDIF
      IF ( j1 .le. mh ) THEN 
         jbs = 1
c         j1 = 2
         j1 = 1
         IF ( bcy .eq. 2 ) j2 = ny - 1
      ELSE 
         jbs = 0
      ENDIF
      IF ( j2 .ge. ny-mh ) THEN 
         jbn = 1
         j2 = ny-1
         IF ( bcy .eq. 2 ) j1 = 1
      ELSE 
         jbn = 0
      ENDIF
      IF ( k2 .ge. nzt1-10 ) THEN 
         kbt = 1
         k2 = nz-1
      ELSE 
         kbt = 0
      ENDIF
      IF ( k1 .le. 5 ) THEN 
         k1 = 1
      ENDIF
        
c set up Neumann boundary conditions for du/dn = 0.0
c on top and side edges.  Bottom edge assumed to be Dirichlet
c and set before calling SOR.
c  top/bottom
       IF (ibc .eq. 1 .or. kbt .eq. 1 ) THEN
        DO j=j1,j2
        DO i=i1,i2
          phi(i,j,nz) = phi(i,j,nz-1)
c          phi(i,j,1) = 0.0
        END DO
        END DO
       END IF
c east/west
!       IF (ibc .eq. 1 .or. ibw .eq. 1 ) THEN
       IF (ibcw .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
c          phi(1,j,k) = phi(2,j,k)
          phi(0,j,k) = phi(1,j,k)
        END DO
        END DO
       END IF
!       IF (ibc .eq. 1 .or. ibe .eq. 1 ) THEN
       IF (ibce .eq. 1 ) THEN
        DO k=k1,k2
        DO j=j1,j2
          phi(nx,j,k) = phi(nx-1,j,k)
        END DO
        END DO
       END IF
c north/south
!       IF (ibc .eq. 1 .or. jbs .eq. 1 ) THEN
       IF (ibcs .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
c          phi(i,1,k) = phi(i,2,k)
          phi(i,0,k) = phi(i,1,k)
        END DO
        END DO
       END IF
!       IF (ibc .eq. 1 .or. jbn .eq. 1 ) THEN
       IF (ibcn .eq. 1 ) THEN
        DO k=k1,k2
        DO i=i1,i2
          phi(i,ny,k) = phi(i,ny-1,k)
        END DO
        END DO
       END IF
        
c       ELSEIF (ibc .eq. 2) THEN ! periodic on sides
c east/west
c        DO k=k1,k2
c        DO j=j1,j2
c          phi(1,j,k) = phi(nx-1,j,k)
c          phi(nx,j,k) = phi(2,j,k)
c        END DO
c        END DO
c north/south
c        DO k=k1,k2
c        DO i=i1,i2
c          phi(i,1,k) = phi(i,ny-1,k)
c          phi(i,ny,k) = phi(i,2,k)
c        END DO
c        END DO
c        
c       END IF

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
        
c        DO k=Max(1,k1-1),k2+1
cc         anorm(k) = 0.0
cc         err(k)   = 0.0
c         residmax(k) = 0.0
c        END DO
        remax1 = 0.0
        
         anormc = 0.00
         err1 = 0.00
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = i1
      ifn = i2
      jst = j1
      jfn = j2
c      if (nyd.eq.1) jfn = ny-1
      kst = k1
c      if (nze.eq.1) kst = 3
      kfn = k2

      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'n,i1,j1,k1,i2,j2,k2 = ',n,i1,j1,k1,i2,j2,k2,nz,nzt1
      ENDIF
        
c      if (nzf.eq.1) kfn = nz-1
c
c     set periodic boundaries if necessary
c
c
c   red (x,y) on odd z planes
c
      j1 = Ny
      i1 = Nx
      k1 = Nz
            ptstart  = 0
            ptend    = 1
            ptstride = 1
      
         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

            DO  k=1, k1
               DO j=1, j1
                  !
                  ibeg=mod((iGslg-1)+(j+jGslg-1)+(k+kGslg-1)+pts,2) + 1
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
                 phi(i,j,k) = (rhs(i,j,k) - (
     +                 cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +                 cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +                 cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +                 (cofx(i,3)+cofy(j,3)+cofz(k,3))
                ENDIF
                  
                  END DO
                  !
               END DO
               !
            END DO

c      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)

        ng = 1

#ifdef MPI        
        IF ( pts == ptend ) THEN
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,1,
     :        w_proc(my_rank),e_proc(my_rank),westward_tag,phi)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,1,
     :        w_proc(my_rank),e_proc(my_rank),eastward_tag,phi)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,1,
     :        n_proc(my_rank),s_proc(my_rank),southward_tag,phi)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,1,
     :        n_proc(my_rank),s_proc(my_rank),northward_tag,phi)

        ENDIF
#endif

      ENDDO


      remax1 = 0.0
c !$OMP PARALLEL DO 
c !$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz,mask),
c !$OMP+PRIVATE(i,j,k,resid),
c C$OMP+REDUCTION (max : remax1)
c      do k=kst,kfn
c      do j=jst,jfn
c        do i=ist,ifn
      DO k=1, k1
      DO j=1, j1
      DO i=1, i1
         IF ( mask(i,j,k) .gt. 0.5e30 ) THEN
          resid =  rhs(i,j,k)-(
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1) +
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))*phi(i,j,k)  )
          remax1 = Max(remax1,Abs(resid))
         ENDIF
        end do
      end do
      end do


#ifdef MPI
       
       remax1red = remax1

       CALL MPI_AllReduce(remax1, remax1red, 1, MPI_REAL, 
     :                   MPI_MAX, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       remax1 = remax1red
#endif
        
        IF (ndebug .ge. 1 ) ! .and. (n/nout)*nout.eq.n) 
     :     write(6,'(i5,4(3x,e13.6))') 
     :     n,anormc,anormc/anormf,remax1,remax1*phiminv ! ,
c     :     anormh1,anormh1/anormf
       
        
        IF (remax1*phiminv .lt. eps .and. n .ge. minits ) GOTO 999

      END DO ! (n=1,maxits)
       
 999  CONTINUE
        numits = n
c       write(6,'(a,i7)') 'number of iterations = ',n
c       write(6,'(i5,4(3x,e13.6))') 
c     :     n,anormc,anormc/anormf,remax1,remax1*phiminv 
      RETURN
      END
