!-----------------------------------------------------------------------------
!
! /////////////////////          BEGIN          \\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE  BUOY     ////////////////////
!
! BUOY computes the buoyancy term for the vertical equation of motion
!-----------------------------------------------------------------------------
! Created by Louis Wicker, May 2004
!
! Updated for new data structures
!-----------------------------------------------------------------------------
 SUBROUTINE BUOY(fw,s,st,sinit,nx,ny,nz,ns)

  USE GRID_MODULE
  USE PARAM_MODULE
!#ifdef MPI
  USE COMMASMPI_MODULE
!#endif

  implicit none

  integer, INTENT(IN)    :: nx,ny,nz,ns
  real,    INTENT(INOUT) :: fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
  real,    INTENT(IN)    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
  TYPE(VARIABLE)         :: s(ns)
  real,  INTENT(IN)  :: sinit(-ng+1:nz+ng,ns)
  
   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
! Local variables

  integer i, j, k, n
  integer kt,kb

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

  integer :: ixb, jyb, kzb
  integer :: ixe, jye, kze


!-----------------------------------------------------------------------------
! BUOYANCY term


  DO n = 1,ns 
       
!   print *, n, s(n)%name, s(n)%buotype
       
   SELECT CASE( s(n)%buotype )

! THETA

    CASE ( 1 )
    
      kb = 2
      kt = 1

!#ifdef MPI
     kzb = -ng+2
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 2
     if (kzend .eq. nzend) kze = kzend-kzbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

! erm added temporary array fz2 to thwart issues with vectorization
! on Intel SSE processors.  Inconsistent floating point results from
! vectorization produce noise, so can either reduce optimization or force
! high-precision temporaries.  fx2 should stay in cache and be fast.
! Cuts the number of divisions in half, anyway

     do j = jyb,jye ; do i = ixb,ixe
     fz2(i,j,kb) = st(i,j,kzb-1,n)/sinit(kzb-1,n) - 1.0
     ENDDO; ENDDO

     do k = kzb,kze 

     do j = jyb,jye ; do i = ixb,ixe
     fz2(i,j,kt) = st(i,j,k,n)/sinit(k  ,n) - 1.0
     ENDDO; ENDDO
     
     do j = jyb,jye ; do i = ixb,ixe
!#else
!     
!     fz(ixb:ixe,jyb:jye,kb) = st(ixb:ixe,jyb:jye,1,n)/sinit(1  ,n) - 1.0
!     DO k = 2,nz-1; DO j = 1,ny-1; DO i = 1,nx-1
!#endif
      fw(i,j,k) = fw(i,j,k) + 0.5*g*(fz2(i,j,kt) + fz2(i,j,kb))
!                * (st(i,j,k,n)/sinit(k  ,n) - 1.0                    &
!                  +st(i,j,k-1,n)/sinit(k-1,n) - 1.0)
     ENDDO ; ENDDO
     
      kb = 3 - kb
      kt = 3 - kt

     ENDDO

! QV MIXING RATIO

    CASE ( 2 )

      kb = 2
      kt = 1

!#ifdef MPI
     kzb = -ng+2
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 2
     if (kzend .eq. nzend) kze = kzend-kzbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     do j = jyb,jye ; do i = ixb,ixe
     fz2(i,j,kb) = st(i,j,kzb-1,n) - sinit(kzb-1,n)
     ENDDO; ENDDO


     do k = kzb,kze
     
     do j = jyb,jye ; do i = ixb,ixe
     fz2(i,j,kt) = st(i,j,k,n) - sinit(k  ,n)
     ENDDO; ENDDO
     
     do j = jyb,jye ; do i = ixb,ixe
!#else
!     DO k = 2,nz-1; DO j = 1,ny-1; DO i = 1,nx-1
!#endif
      fw(i,j,k) = fw(i,j,k) + 0.5*g*0.61*(fz2(i,j,kt) + fz2(i,j,kb))

!      fw(i,j,k) = fw(i,j,k) + 0.5*g                                  &
!                * 0.61*(st(i,j,k ,n ) - sinit(k  ,n)                 &
!                       +st(i,j,k-1,n) - sinit(k-1,n))
     ENDDO ; ENDDO
     
      kb = 3 - kb
      kt = 3 - kt
     
     ENDDO

! MIXING RATIO FOR HYDROMETEORS

    CASE ( 3 )

!#ifdef MPI
     kzb = -ng+2
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 2
     if (kzend .eq. nzend) kze = kzend-kzbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
!#else
!     DO k = 2,nz-1; DO j = 1,ny-1; DO i = 1,nx-1
!#endif
      fw(i,j,k) = fw(i,j,k) - 0.5*g*(st(i,j,k,n)+st(i,j,k-1,n))
     ENDDO ; ENDDO ; ENDDO

    END SELECT
          
   ENDDO

 RETURN
 END SUBROUTINE BUOY
