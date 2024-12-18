!#include "sam.def.h"
#define SWM
!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
      subroutine boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,vt,   &
     &  a,db1,imapz,mzdist,ia,id,xfall)
!
! Box-lagrangian fallout scheme of Kato (1995, J. Met. Soc. Japan)
!
!  Written by ERM 2/28/2002
!
! 7.6.2005: Note: could also pass inverse of density to save some divides
!
! 4.24.2005 Changed input density to a pre-filled slab that has either the
!           air density or a constant value of 1.  This saves doing for each
!           variable.  Also added imn,imx,kmn,kmx to see if it makes any difference
!
!

      USE COMMASMPI_MODULE
      USE MICRO_MODULE

      implicit none
      
      integer nx,ny,nz,nor,ngt,jgs,imapz,na,ia,mzdist
      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)

!#ifdef CM1
!      integer, parameter :: norz = 1
!      real gz(-norz+ng1:nz+norz),z1d(-norz+ng1:nz+norz,4)
!!      real a(nx,ny,nz,na)
!      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na) ! quantity to be 'advected'
!#else
      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
!      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
!#endif
!      real db(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor) ! air density
      real vt(nx,nz+1)  ! terminal speed for a
      real dtp,dz,dtz1
      real qtmp1(nx,nz+1),qtmp2(nx,nz+1)
      real cmax
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      real zw(0:nz+1),zs(0:nz+1),dzw(nz+1)
      real db1(nx,nz+1)
      real ZL2,ZL1,zt,zb,dbrat
      integer k0, L2,L1
      integer ndebug1
            
      integer ix,jy,kz,ndfall,n,k
      integer iv1,iv2
      real tmp
      integer imn,imx,kmn,kmx

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.

! ###################################################################

!      write(0,*) 'boxfall : ia = ',ia

      jy = jgs

      iv1 = 0
      iv2 = 0

      imn = nx-1
      imx = 1
      kmn = nz-1
      kmx = 1

      cmax = 0.0

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO kz = kzb,kze
        DO ix = ixb,ixe
         cmax = Max(cmax, vt(ix,kz)*z1d(kz,3)) !*gz(kz))
         
         IF ( nprock > 1 .or. ipconc > 5 ) THEN ! limit fall speeds to prevent needing lagrange switch
           vt(ix,kz) = Min( vt(ix,kz), 0.95/(z1d(kz,3)*dtp) )
         ENDIF
!         print*, 'boxfall',ix,jy,kz,ia,a(ix,jy,kz,ia)

         qtmp1(ix,kz) = a(ix,jgs,kz,ia)
         qtmp2(ix,kz) = 0.0
         
         IF ( a(ix,jgs,kz,ia) .ne. 0.0 ) THEN
           imn = Min(ix,imn)
           imx = Max(ix,imx)
           kmn = Min(kz,kmn)
           kmx = Max(kz,kmx)
         ENDIF
        ENDDO
      ENDDO
      

      kzb = 0
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO kz = kzb,kze
       zs(kz) = z1d(kz,1)
       zw(kz) = z1d(kz+1,2) ! height of _top_ of the scalar box kz

!      IF ( myprock == 1 ) THEN
!        write(0,*) my_rank, 'boxfall: k,zs,zw,kz = ',kz,zs(kz),zw(kz),z1d(kz,1),z1d(kz,2)
!      ENDIF

      ENDDO

!       zw(0) = Max( 0.0, z1d(1,2) )

!#ifdef SWM
!      DO kz = 1,nz-1
!       dzw(kz) = 1./z1d(kz,4) !  zw(kz) - zw(kz-1)
!      ENDDO
!#else
!      DO kz = 1,nz-1
!       dzw(kz) = zw(kz) - zw(kz-1)
!      ENDDO
!#endif
      
! first check if fallout is worth doing
      IF ( cmax .eq. 0.0 .or. imn .gt. imx ) THEN
        RETURN
      ENDIF
      

!
!  Only want to divide/multiply by air density for mixing ratio, so
!   for others (charge density, rime time, etc.) just set temp dens = 1.0
!      
!      DO kz=1,nz-1
!        DO ix=1,nx-1
!         IF ( id .eq. 1 ) THEN
!           db1(ix,kz) = db(ix,jy,kz)
!         ELSE
!           db1(ix,kz) = 1.0
!         ENDIF
!        ENDDO
!      ENDDO
      

#ifdef MPI
      kzb = 1
      kze = ktile+1
!      if (kzbeg .le. kmn) kzb = kmn
!      if (kzend .ge. kmx) kze = kmx-kzbeg+1
      kze = kmx

      ixb = imn
      ixe = itile
!      if (ixbeg .le. imn) ixb = imn
!      if (ixend .ge. nxend) ixe = ixend-ixbeg
      ixe = imx

      do kz = kzb,kze
       k0 = kz
      do ix = ixb,ixe
#else
      do kz = kmn,kmx ! 1,nz-1
       k0 = kz
      do ix = imn,imx ! 1,nx-1
#endif
        zt = zw(kz)   - dtp*vt(ix,kz)
        zb = zw(kz-1) - dtp*vt(ix,kz)
        
        IF ( vt(ix,kz) .gt. 400.0 ) THEN
          ndebug1 = 1
        ELSE
          ndebug1 = 0
        ENDIF
!        
! find w-levels just above zt and zb
!
#ifdef MPI
        kzb = ktile+1
        kze = 0
        if (kzbeg .eq. nzbeg) kze = 0
        if (kzend .ge. kzbeg-1+k0)    kzb = k0

        DO k=kzb,kze,-1
#else
        DO k=k0,0,-1
#endif
          IF ( zw(k) .ge. zt ) THEN
            L2 = k
            ZL2 = zw(k)
          ENDIF
          IF ( zw(k) .ge. zb ) THEN
            L1 = k
            ZL1 = zw(k)
          ENDIF
        ENDDO
        IF ( ndebug1 .eq. 1 ) THEN
          write(6,*) 'k0,ix,jy,L2,ZL2,L1,ZL1,zt,zb,vt,ia,a',   &
     &      k0,ix,jy,L2,ZL2,L1,ZL1,zt,zb,vt(ix,kz),ia,a(ix,jgs,kz,ia)
           write(6,*) 'dtp = ',dtp
         ENDIF
       
      IF ( L2 .eq. k0 .and. L1 .eq. k0 - 1 ) THEN ! Eulerian, CFL < 1.0
       iv1 = iv1 + 1

#ifdef SWM
       qtmp2(ix,kz) =  qtmp2(ix,kz) -   &
     &         qtmp1(ix,kz)*(ZL2-zt)*z1d(kz,3) ! *gz(kz) ! /dz*gt(ix,jy,kz,imapz)
       IF ( kz .gt. 1 ) THEN
       qtmp2(ix,kz-1) = qtmp2(ix,kz-1) +   &
     &         qtmp1(ix,kz)*(ZL2-zt) *z1d(kz-1,3)* & ! *gz(kz-1)* ! /dz*gt(ix,jy,kz-1,imapz)*
     &         db1(ix,kz)/db1(ix,kz-1)
       ELSE
         xfall(ix,jy,ia) = xfall(ix,jy,ia) +     &
     &         qtmp1(ix,kz)*(ZL2-zt)*z1d(kz,3)*  & ! *gz(kz)* ! /dz*gt(ix,jy,1,imapz)*
     &         db1(ix,kz)/db1(ix,1)
       ENDIF
#else
       qtmp2(ix,kz) =  qtmp2(ix,kz) -   &
     &         qtmp1(ix,kz)*(ZL2-zt)/dz*gt(ix,jy,kz,imapz)
       IF ( kz .gt. 1 ) THEN
       qtmp2(ix,kz-1) = qtmp2(ix,kz-1) + &
     &         qtmp1(ix,kz)*(ZL2-zt)/dz*gt(ix,jy,kz-1,imapz)* &
     &         db1(ix,kz)/db1(ix,kz-1)
       ELSE
         xfall(ix,jy,ia) = xfall(ix,jy,ia) + &
     &         qtmp1(ix,kz)*(ZL2-zt)/dz*gt(ix,jy,1,imapz)*&
     &         db1(ix,kz)/db1(ix,1)
       ENDIF
#endif
 
!c the following 4 lines also work
!       tmp = qtmp1(ix,kz)*dtp*vt(ix,kz)/dz*gt(ix,jy,kz,imapz)
!       qtmp2(ix,kz) =  qtmp2(ix,kz) - tmp
!       tmp = qtmp1(ix,kz)*dtp*vt(ix,kz)/dz*gt(ix,jy,kz-1,imapz)
!       qtmp2(ix,kz-1) = qtmp2(ix,kz-1) + tmp*db1(ix,kz)/db1(ix,kz-1)

      ELSE

! here we do the box-lagrangian scheme if v*dt/dz > 1.0
       iv2 = iv2 + 1
!       write(0,*) my_rank, 'boxfall: kz,ia,vt,l1,l2,k0 = ',kz, ia,vt(ix,kz),qtmp1(ix,kz),l1,l2,k0
!
! First assume that all stuff leaves the box, then put back in later
!
        qtmp2(ix,k0) = qtmp2(ix,k0) - qtmp1(ix,k0)
      
      IF ( zt .gt. 0.0 .and. k0 .gt. 0 ) THEN
       DO k=L2,L1,-1
        dbrat = db1(ix,k0)/db1(ix,Max(1,k))
#ifdef SWM
       IF ( k .eq. 0 ) THEN
         xfall(ix,jy,ia) = xfall(ix,jy,ia) + &
     &         qtmp1(ix,kz)*(-zb)*dbrat*gz(kz-1)  ! /dz*gt(ix,jy,kz-1,imapz)*dbrat
        ELSEIF ( zw(k) .gt. zt .and. zw(k-1) .lt. zt ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) + &
     &        qtmp1(ix,k0)*(zt-zw(k-1))*gz(k)*dbrat  !  /dz*gt(ix,jy,k,imapz)*dbrat
        ELSEIF ( zw(k) .lt. zt .and. zw(k-1) .gt. zb ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) +    &
     &                  qtmp1(ix,k0)*dbrat
        ELSEIF ( zw(k) .gt. zb .and. zw(k-1) .lt. zb ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) +  &
     &         qtmp1(ix,k0)*(zw(k)-zb)*gz(k)*dbrat !  /dz*gt(ix,jy,k,imapz)*dbrat
        ENDIF
#else
        IF ( k .eq. 0 ) THEN
         xfall(ix,jy,ia) = xfall(ix,jy,ia) +  &
     &         qtmp1(ix,kz)*(-zb)/dz*gt(ix,jy,kz-1,imapz)*dbrat
        ELSEIF ( zw(k) .gt. zt .and. zw(k-1) .lt. zt ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) +    &
     &        qtmp1(ix,k0)*(zt-zw(k-1))/dz*gt(ix,jy,k,imapz)*dbrat
        ELSEIF ( zw(k) .lt. zt .and. zw(k-1) .gt. zb ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) +  &
     &                  qtmp1(ix,k0)*dbrat
        ELSEIF ( zw(k) .gt. zb .and. zw(k-1) .lt. zb ) THEN
          qtmp2(ix,k) = qtmp2(ix,k) +   &
     &         qtmp1(ix,k0)*(zw(k)-zb)/dz*gt(ix,jy,k,imapz)*dbrat
        ENDIF
#endif
       ENDDO
      ELSEIF ( zt .lt. 0.0 ) THEN
         xfall(ix,jy,ia) = xfall(ix,jy,ia) +   &
     &                     qtmp1(ix,k0)*db1(ix,k0)/db1(ix,1)
      ENDIF
        
!        IF ( zt .gt. 0.0 .and. zb .lt. 0.0 ) THEN
!          xfall(ix,jy,ia) = xfall(ix,jy,ia) +  &
!     &             qtmp1(ix,k0)*(db1(ix,k0)*(-zb))/(db1(ix,1)*dzw(1))
!        ENDIF
        
      ENDIF 
      
      enddo
      enddo
! 
! print out when the box scheme was needed and for what variable
!
!      IF ( iv2 .gt. 0 ) THEN
!        write(0,*) my_rank, 'boxfall: jy,ia,iv1,iv2 = ',jy,ia,iv1,iv2
!      ENDIF  

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
      do ix = 1,nx-1
#endif
        a(ix,jgs,kz,ia) =  a(ix,jgs,kz,ia) + qtmp2(ix,kz) 
      enddo
      enddo


      
      RETURN
      END


!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
      subroutine fallout(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,vt,   &
     &  a,db1,imapz,mzdist,ia,id,xfall,dtz1)
!
! First-order, upwind fallout scheme
!
!  Written by ERM 6/10/2011
!
!
!

      USE COMMASMPI_MODULE
      USE MICRO_MODULE

      implicit none
      
      integer nx,ny,nz,nor,ngt,jgs,imapz,na,ia,mzdist
      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)

!#ifdef CM1
!      integer, parameter :: norz = 1
!      real gz(-norz+ng1:nz+norz),z1d(-norz+ng1:nz+norz,4)
!!      real a(nx,ny,nz,na)
!      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na) ! quantity to be 'advected'
!#else
      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
!      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
!#endif
!      real db(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor) ! air density
      real vt(nx,nz+1)  ! terminal speed for a
      real dtp,dz
      real qtmp1(nx,nz+1),qtmp2(nx,nz+1)
      real cmax
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      real zw(0:nz+1),zs(0:nz+1),dzw(nz+1)
      real db1(nx,nz+1),dtz1(nx,nz+1)
      real ZL2,ZL1,zt,zb,dbrat
      integer k0, L2,L1
      integer ndebug1
            
      integer ix,jy,kz,ndfall,n,k
      integer iv1,iv2
      real tmp
      integer imn,imx,kmn,kmx

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.

! ###################################################################

      jy = jgs

      iv1 = 0
      iv2 = 0

      imn = nx-1
      imx = 1
      kmn = nz-1
      kmx = 1

      cmax = 0.0

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      qtmp1(:,:) = 0.0
      
!      write(0,*) 'fallout : ia = ',ia
      
      DO kz = kzb,kze
        DO ix = ixb,ixe
         cmax = Max(cmax, vt(ix,kz)*z1d(kz,3)) !*gz(kz))
         
         qtmp1(ix,kz) = a(ix,jgs,kz,ia)*vt(ix,kz)*db1(ix,kz)
         qtmp2(ix,kz) = 0.0
         
         IF ( a(ix,jgs,kz,ia) .ne. 0.0 ) THEN
           imn = Min(ix,imn)
           imx = Max(ix,imx)
           kmn = Min(kz,kmn)
           kmx = Max(kz,kmx)
         ENDIF
        ENDDO
      ENDDO
      
      kmn = Max(1,kmn-1)

      kzb = 0
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
! first check if fallout is worth doing
      IF ( cmax .eq. 0.0 .or. imn .gt. imx ) THEN
        RETURN
      ENDIF
      

      
      IF ( kmn == 1 .and. myprock == 1 ) THEN
      
      kz = 1
#ifdef MPI
      kzb = kmn
      kze = ktile+1
!      if (kzbeg .le. kmn) kzb = kmn
!      if (kzend .ge. kmx) kze = kmx-kzbeg+1
      kze = kmx

      ixb = imn
      ixe = itile
!      if (ixbeg .le. imn) ixb = imn
!      if (ixend .ge. nxend) ixe = ixend-ixbeg
      ixe = imx

!      do kz = kzb,kze
      do ix = ixb,ixe
#else
      kzb = kmn
!      do kz = kmn,kmx ! 1,nz-1
      do ix = imn,imx ! 1,nx-1
#endif        

         xfall(ix,jy,ia) = xfall(ix,jy,ia) + dtp*a(ix,jgs,kz,ia)*vt(ix,kz)*z1d(kz,3)  ! *gz(kz)* ! /dz*gt(ix,jy,1,imapz)*
 
!c the following 4 lines also work
!       tmp = qtmp1(ix,kz)*dtp*vt(ix,kz)/dz*gt(ix,jy,kz,imapz)
!       qtmp2(ix,kz) =  qtmp2(ix,kz) - tmp
!       tmp = qtmp1(ix,kz)*dtp*vt(ix,kz)/dz*gt(ix,jy,kz-1,imapz)
!       qtmp2(ix,kz-1) = qtmp2(ix,kz-1) + tmp*db1(ix,kz)/db1(ix,kz-1)

      
      enddo
!      enddo
      
      ENDIF
! 
! print out when the box scheme was needed and for what variable
!
!      IF ( iv2 .gt. 0 ) THEN
!        write(0,*) my_rank, 'boxfall: jy,ia,iv1,iv2 = ',jy,ia,iv1,iv2
!      ENDIF  

!#ifdef MPI
      kzb = kmn
      kze = kmx
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = imn
      ixe = imx
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      do ix = ixb,ixe
!#else
!      do kz = 1,nz-1
!      do ix = 1,nx-1
!#endif
        a(ix,jgs,kz,ia) =  a(ix,jgs,kz,ia) + dtp*dtz1(ix,kz)*(qtmp1(ix,kz+1) - qtmp1(ix,kz) )
      enddo
      enddo


      
      RETURN
      END

!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
      subroutine fallouttak(nx,ny,nz,nor,na,z1d,dtp,dz,jgs,vt,   &
     &  a,ia,xfall,rhovt,sx,lfall,lnfall)
!
! First-order, upwind fallout scheme for use with bin physics
!
!  Written by ERM 6/10/2011
!
!
!

      USE COMMASMPI_MODULE
      USE MICRO_MODULE

      implicit none
      
      integer nx,ny,nz,nor,ngt,jgs,na,ia,lfall,lnfall
!      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)
      double precision, intent(in) :: sx ! particle mass (kg)

!#ifdef CM1
!      integer, parameter :: norz = 1
!      real gz(-norz+ng1:nz+norz),z1d(-norz+ng1:nz+norz,4)
!!      real a(nx,ny,nz,na)
!      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na) ! quantity to be 'advected'
!#else
!      real gz(-nor+ng1:nz+nor)
       real z1d(-nor+ng1:nz+nor,4)
!      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
      real rhovt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
!#endif
!      real db(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor) ! air density
      real*8 ::  vt ! terminal speed for a, which is a constant for a size bin
      real dtp,dz
      real qtmp1(nx,nz+1),qtmp2(nx,nz+1)
      real cmax
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      real zw(0:nz+1),zs(0:nz+1),dzw(nz+1)
!      real dtz1(nx,nz+1)
      real ZL2,ZL1,zt,zb,dbrat
      integer k0, L2,L1
      integer ndebug1
            
      integer ix,jy,kz,n,k
      integer iv1,iv2
      real tmp
      integer imn,imx,kmn,kmx
      
      real    :: dtptmp
      integer :: ndfall

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.

! ###################################################################

      jy = jgs

      iv1 = 0
      iv2 = 0

      imn = nx-1
      imx = 1
      kmn = nz-1
      kmx = 1

      cmax = 0.0

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      qtmp1(:,:) = 0.0
      
      DO kz = kzb,kze
        DO ix = ixb,ixe
         cmax = Max(cmax, vt*rhovt(ix,jy,kz)*z1d(kz,3)) !*gz(kz))
         
         qtmp1(ix,kz) = a(ix,jgs,kz,ia)*vt*rhovt(ix,jy,kz)
         qtmp2(ix,kz) = 0.0
         
         IF ( a(ix,jgs,kz,ia) .ne. 0.0 ) THEN
           imn = Min(ix,imn)
           imx = Max(ix,imx)
           kmn = Min(kz-1,kmn)
           kmx = Max(kz,kmx)
         ENDIF
        ENDDO
      ENDDO
      
      kmn = Max(1,kmn-1)

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
! first check if fallout is worth doing
      IF ( cmax .eq. 0.0 .or. imn .gt. imx ) THEN
        RETURN
      ENDIF
      
      IF ( dtp*cmax .lt. 0.7 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*cmax/0.7) + 1)
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp/Real(ndfall)
!        write(0,*) 'subdivide fallout on my_rank = ',my_rank
!        write(0,*) 'for il,jyc = ',il,jy,dtp*vtmax
      ELSE
        dtptmp = dtp
      ENDIF

      DO n = 1,ndfall

      IF ( n > 1 ) THEN
      
      kmn = Max(1,kmn-1)
      DO kz = kzb,kze
        DO ix = ixb,ixe
         
         qtmp1(ix,kz) = a(ix,jgs,kz,ia)*vt*rhovt(ix,jy,kz)
         qtmp2(ix,kz) = 0.0
 

        ENDDO
      ENDDO
      ENDIF
      
      IF ( kmn == 1 .and. myprock == 1 .and. lfall > 0 ) THEN
      
      kz = 1
#ifdef MPI
      kzb = kmn
      kze = ktile+1
!      if (kzbeg .le. kmn) kzb = kmn
!      if (kzend .ge. kmx) kze = kmx-kzbeg+1
      kze = kmx

      ixb = imn
      ixe = itile
!      if (ixbeg .le. imn) ixb = imn
!      if (ixend .ge. nxend) ixe = ixend-ixbeg
      ixe = imx

!      do kz = kzb,kze
      do ix = ixb,ixe
#else
      kzb = kmn
!      do kz = kmn,kmx ! 1,nz-1
      do ix = imn,imx ! 1,nx-1
#endif        

         xfall(ix,jy,lfall) = xfall(ix,jy,lfall) +          &
     &          1.e6*dtptmp*a(ix,jgs,kz,ia)*vt*rhovt(ix,jy,kz)*z1d(kz,3)*sx  ! 1.e6 to convert cm^-3 to m^-3 in num. conc.

        IF ( lnfall > 0 ) THEN
         xfall(ix,jy,lnfall) = xfall(ix,jy,lnfall) +          &
     &          1.e6*dtptmp*a(ix,jgs,kz,ia)*vt*rhovt(ix,jy,kz)*z1d(kz,3)  ! 1.e6 to convert cm^-3 to m^-3 in num. conc.
        ENDIF
      
      enddo
!      enddo
      
      ENDIF
! 
! print out when the box scheme was needed and for what variable
!
!      IF ( iv2 .gt. 0 ) THEN
!        write(0,*) my_rank, 'boxfall: jy,ia,iv1,iv2 = ',jy,ia,iv1,iv2
!      ENDIF  

#ifdef MPI
      kzb = kmn
      kze = kmx
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = imn
      ixe = imx
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = kmn,kmx
      do ix = imn,imx
#endif
        a(ix,jgs,kz,ia) =  a(ix,jgs,kz,ia) + dtptmp*z1d(kz,3)*(qtmp1(ix,kz+1) - qtmp1(ix,kz) )
      enddo
      enddo
      
      ENDDO ! ndfall


      
      RETURN
      END




