
!------------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE MFILTER6  ////////////////////
!
! FILTER6 performs a monotonic sixth-order filter on the passed variable (Xue 2000, MWR)
! Filtering only performed in the horiztonal.
!------------------------------------------------------------------------------
! Created by Louis Wicker and Ted Mansell 6/2006
! Latest update: 06-28-2006
! Update notes:
!------------------------------------------------------------------------------
      SUBROUTINE MFILTER6(s,sb,fs,svar,dt,nx,ny,nz,gx,gy)

   USE GRID_MODULE
   USE PARAM_MODULE
   USE COMMASMPI_MODULE

      implicit none

! Passed variables

   integer                :: nx, ny, nz
   real                   :: dt
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real                   :: sb(-ng+1:nz+ng)

   TYPE(VARIABLE)         :: gx(4) , gy(4)
   TYPE(VARIABLE)         :: svar

! Local variables

      integer i, j, k, is, js, ks

      real            :: f6dmph, f4dmph, f2dmph
      real, parameter :: f6damp = 1.0/64.
      real            :: f1, sgrad
      real            :: qid
      real, parameter :: fdamp0 = 0.25
      real            :: fdamp
      real            :: xmask(-ng+1:nx+ng), ymask(-ng+1:ny+ng)

      integer         :: imn,imx,jmn,jmx,kmn,kmx,imx1,jmx1,kmx1
      integer         :: im1, im2, ip1, jm1, jm2, jp1, km1, km2, kp1

      integer imxb, jmxb
      integer i1,i2,i3
      integer j1,j2,j3
      integer k1,k2,k3

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
   
      logical :: debug_mpi = .false.

      real,  allocatable :: fx(:,:), fy(:,:)


!------------------------------------------------------------------------------
! Set consts:  f2dmp should have same damping at 2 dx as f4dmp

      f6dmph = f6damp * damph
      f4dmph = f6dmph * 4.0
      f2dmph = f4dmph * 4.0
      fdamp  = fdamp0 * damph
      
      is = svar%istag
      js = svar%jstag
      ks = svar%kstag


!------------------------------------------------------------------------------
! Special horizontal filter

!       f4dmph = f4damp * damph * (1.0 + 2.0*(float(k)/float(nz-1))**4)
!       f2dmph = f4dmph * 4.0

#ifdef MPI
    imn = 1+is
    imx = nxend-1+is
    jmn = 1+js
    jmx = nyend-1+js
!    kmn = 1+ks
    kmn = 1
    kmx = nzend-1 
    imx1 = nxend-1
    jmx1 = nyend-1
    kmx1 = nzend-1
    IF ( bcy == 2 ) THEN
      jmn = 1
!      jmx1 = nyend
      jmx = nyend+js
    ENDIF
    IF ( bcx == 2 ) THEN
      imn = 1
!      imx1 = nxend
      imx = nxend+is
    ENDIF
#else

    imn = nx-1
    imx = 1
    jmn = ny-1
    jmx = 1
    kmn = nz-1
    kmx = 1

    DO k = 1,nz-1+ks
     DO j = 1,ny-1+js
      DO i = 1,nx-1+is

       IF ( s(i,j,k) .ne. 0.0 ) THEN
       
         imn = Min(imn,Max(1+is,i-2))
         imx = Max(imx,Min(i+2,nx-1+is))
         jmn = Min(jmn,Max(1+js,j-2))
         jmx = Max(jmx,Min(j+2,ny-1+js))
         kmn = Min(kmn,Max(1+ks,k-2))
         kmx = Max(kmx,Min(k+2,nz-1+ks))
            
       ENDIF

      ENDDO
     ENDDO
    ENDDO
#endif

    IF ( imn .gt. imx ) RETURN

    imx1 = Min(imx,nx-1)
    jmx1 = Min(jmx,ny-1)
    kmx1 = Min(kmx,nz-1)

#ifdef MPI
    i3 = nxend-3+is
    i2 = nxend-2+is
    i1 = nxend-1+is

    j3 = nyend-3+js
    j2 = nyend-2+js
    j1 = nyend-1+js

    k3 = nzend-3+ks
    k2 = nzend-2+ks
    k1 = nzend-1
    
    imxb = imx
    jmxb = jmx
#else
    i3 = nx-3+is
    i2 = nx-2+is
    i1 = nx-1+is

    j3 = ny-3+js
    j2 = ny-2+js
    j1 = ny-1+js

    k3 = nz-3+ks
    k2 = nz-2+ks
    k1 = nz-1

    imxb = imx
    jmxb = jmx
#endif

    
    allocate (  fx(-ng+1:nx+ng,-ng+1:ny+ng) )
    allocate (  fy(-ng+1:nx+ng,-ng+1:ny+ng) )
    fx(:,:) = 0.0
    fy(:,:) = 0.0

!  xmask(:) = 1.0 ; xmask(1) = 0.0 ; xmask(nx-1+is:nx) = 0.0
!  ymask(:) = 1.0 ; ymask(1) = 0.0 ; ymask(ny-1+js:ny) = 0.0

  xmask(:) = 1.0 
  IF ( bcx .ne. 2 ) THEN
#ifdef MPI
    if (ixbeg .eq. nxbeg) xmask(-ng+1:1) = 0.0
    if (ixend .eq. nxend) then
      ixb = imx-ixbeg+1
      ixe = ixend-ixbeg+1+ng

      xmask(ixb:ixe) = 0.0
    endif
#else
    xmask(-ng+1:1) = 0.0 ; xmask(imx:nx+ng) = 0.0
#endif
  ELSE
#ifndef MPI
    imx1 = imx1 + is
!    fx(nx+is,1:ny-1+js,1:nz-1+ks) = fx(1,1:ny-1+js,1:nz-1+ks)
#endif
  ENDIF

  ymask(:) = 1.0 
  IF ( bcy .ne. 2 ) THEN
#ifdef MPI
!    if (jybeg .eq. nybeg)  ymask(-ng+1:1) = 0.0
    if (jyend .eq. nyend) then
      jyb = jmx-jybeg+1
      jye = jyend-jybeg+1+ng

      ymask(jyb:jye) = 0.0
    endif
#else
    ymask(-ng+1:1) = 0.0 ; ymask(jmx:ny+ng) = 0.0
#endif
  ENDIF

#ifdef MPI
  IF ( bcy .ne. 2 )  THEN
   if (jybeg .eq. nybeg)  ymask(-ng+1:1) = 0.0
  ENDIF
#endif


        
#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn)  ixb = imn-ixbeg+1
    if (ixend .ge. imxb) ixe = imxb-ixbeg+1
#else
    ixb = imn
    ixe = imx
    jyb = jmn
    jye = jmx
    kzb = kmn
    kze = kmx 
#endif

    DO k = kzb,kze

       IF ( nx .gt. 2 ) THEN
         DO j = jyb,jye
          DO i = ixb,ixe 

            im1 = max(i-1,1)
            im2 = max(i-2,1)
            ip1 = min(i+1,nx-1+is)
            
            qid = 0.0

      
!            IF( i .ge. 4 .and. i .le. nx-3+is ) THEN
            IF( (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN

             qid = (-s(i-3,j,k)+5.*s(i-2,j,k)-10.*s(i-1,j,k)+10.*s(i,j,k)-5.*s(i+1,j,k)+s(i+2,j,k))/64.
             qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i-1,j,k))) )

!            ELSEIF( i .eq. 3 .or. i .eq. nx-2+is ) THEN
           ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN

              qid = -(s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)) )/ 16.
              qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i-1,j,k))) )

!            ELSE
            ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
!            ELSEIF( i .eq. 2 .or. i .eq. nx-1+is ) THEN

              qid = (s(i,j,k)-s(i-1,j,k))/4.0
              qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i-1,j,k))) )

!            ELSE
!              CYCLE
            ENDIF


            fx(i,j)  =  - fdamp * qid/(dt*gx(3+is)%flt1d(i))

          ENDDO
         ENDDO
!        ENDDO
            
      ENDIF

!--------------------------------------------------------------------------
! COMPUTE Y-INTERFACE FLUXES

        IF( ny .gt. 2 ) THEN

!          DO k = kmn,kmx 

           DO j = jyb,jye

            jm1 = max(j-1,1)
            jm2 = max(j-2,1)
            jp1 = min(j+1,ny-1+is)

            DO i = ixb,ixe

              qid = 0.0

!              IF( j .ge. 4 .and. j .le. ny-3+js ) THEN
              IF( (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN

                qid = (-s(i,j-3,k)+5.*s(i,j-2,k)-10.*s(i,j-1,k)+10.*s(i,j,k)-5.*s(i,j+1,k)+s(i,j+2,k))/64.
                qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i,j-1,k))) )

!              ELSEIF( j .eq. 3 .or. j .eq. ny-2+js ) THEN
              ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN

                qid = -(s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)) ) / 16.
                qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i,j-1,k))) )

!              ELSE
              ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
!              ELSEIF( j .eq. 2 .or. j .eq. ny-1+js ) THEN

                qid = (s(i,j,k)-s(i,j-1,k))/4.0
                qid = qid*amax1(0.0, sign(1.0,qid*(s(i,j,k)-s(i,j-1,k))) )

!              ELSE
!               CYCLE
              ENDIF
      
              fy(i,j)  = - fdamp*qid/(dt*gy(3+js)%flt1d(j))

            ENDDO
           ENDDO

        ENDIF
 
   DO j = jyb,jmx1 
    DO i = ixb,imx1 


      fs(i,j,k) =  fs(i,j,k)                                                                &
            &      - ymask(j)*(fy(i  ,j+1)        - fy(i,j)       )*gy(3+js)%flt1d(j) &
            &      - xmask(i)*(fx(i+1,j  )        - fx(i,j)       )*gx(3+is)%flt1d(i)


     ENDDO
    ENDDO

        ENDDO ! k

      deallocate( fx )
      deallocate( fy )
      RETURN
      END SUBROUTINE MFILTER6
!------------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE FILTER2  ////////////////////
!
! FILTER2 performs a second order filter on the passed variable
!------------------------------------------------------------------------------
! Created by Louis Wicker, April 6, 1988
! Latest update: 01-22-92
! Update notes:
!------------------------------------------------------------------------------
! .01           ! F4damp (f4damp * .0625 is actual coefficient per time step)
!  set damph = fdamph0/dt
      SUBROUTINE FILTER2(s,sb,fs,r,is,js,ks,nx,ny,nz)

   USE GRID_MODULE
   USE PARAM_MODULE

      implicit none

      integer nx, ny, nz, is, js, ks
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real sb(nz)

! Local variables

      integer i, j, k
      real r (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

!------------------------------------------------------------------------------
! Set consts

      real :: f2dmph
      real :: f2damp = .25 
      f2dmph = f2damp * damph

!------------------------------------------------------------------------------
! Create perturbation fields

      DO k = 1,nz
       DO j = 1,ny
        DO i = 1,nx

         r(k,j,i) = s(k,j,i) - sb(k)

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! X pass

      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1
        DO i = 2,nx-2+is

         fs(i,j,k) = fs(i,j,k)    &
     &             + f2dmph*(r(i+1,j,k) + r(i-1,j,k) - 2.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! Y pass

      DO k = 1+ks,nz-1
       DO j = 2,ny-2+js
        DO i = 1+is,nx-1

         fs(i,j,k) = fs(i,j,k)   &
     &             + f2dmph*(r(i,j+1,k) + r(i,j-1,k) - 2.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! Z pass

      DO k = 2,nz-2+ks
       DO j = 1+js,ny-1
        DO i = 1+is,nx-1

         fs(i,j,k) = fs(i,j,k)   &
     &             + f2dmph * (r(i,j,k+1) + r(i,j,k-1) - 2.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FILTER2
!------------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE FILTER4  ////////////////////
!
! FILTER4 performs a fourth order filter on the passed variable
!------------------------------------------------------------------------------
! Created by Louis Wicker, April 6, 1988
! Latest update: 01-22-92
! Update notes:
!------------------------------------------------------------------------------
      SUBROUTINE FILTER4(s,sb,fs,r,svar,nx,ny,nz)

   USE GRID_MODULE
   USE PARAM_MODULE

      implicit none

! Passed variables

      integer nx, ny, nz
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real sb(nz)
      real r (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

   TYPE(VARIABLE) :: svar

! Local variables

      integer i, j, k, is, js, ks

      real :: f4dmph, f2dmph
      real :: f4damp = .0625 

!------------------------------------------------------------------------------
! Set consts:  f2dmp should have same damping at 2 dx as f4dmp

      f4dmph = f4damp * damph
      f2dmph = f4dmph * 4.0

      is = svar%istag
      js = svar%jstag
      ks = svar%kstag

!------------------------------------------------------------------------------
! Special horizontal filter

!       f4dmph = f4damp * damph * (1.0 + 2.0*(float(k)/float(nz-1))**4)
!       f2dmph = f4dmph * 4.0

!------------------------------------------------------------------------------
! Create perturbation fields
!------------------------------------------------------------------------------
      DO k = 1,nz
       DO j = 1,ny
        DO i = 1,nx

         r(i,j,k) = s(i,j,k) - sb(k)

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! X pass

      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1
        DO i = 3,nx-3+is

         fs(i,j,k) = fs(i,j,k)    &
     &             - f4dmph *   ((r(i+2,j,k) + r(i-2,j,k))  &
     &                      - 4.*(r(i+1,j,k) + r(i-1,j,k))  &
     &                      + 6.* r(i,j,k) )

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 2nd order filter at the boundaries

      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1

        fs(2,j,k) = fs(2,j,k)  &
     &   + f2dmph * (r(3,j,k)  + r(1,j,k) - 2.* r(2,j,k))

        fs(nx-2+is,j,k) = fs(nx-2+is,j,k)           &
     &   + f2dmph * (r(nx-1+is,j,k)+r(nx-3+is,j,k)  &
     &           -2.*r(nx-2+is,j,k))

       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! Y pass

      DO k = 1+ks,nz-1
       DO j = 3,ny-3+js
        DO i = 1+is,nx-1

         fs(i,j,k) = fs(i,j,k)                             &
     &             - f4dmph *   ((r(i,j+2,k) + r(i,j-2,k)) &
     &                      - 4.*(r(i,j+1,k) + r(i,j-1,k)) &
     &                      + 6.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 2nd order filter at the boundaries

      DO k = 1+ks,nz-1
       DO i = 1+is,nx-1

        fs(i,2,k) = fs(i,2,k)          &
     &   + f2dmph * ( r(i,3,k)  + r(i,1,k)    - 2.* r(i,2,k) )
        fs(i,ny-2+js,k) = fs(i,ny-2+js,k)               &
     &   + f2dmph * ( r(i,ny-1+js,k) + r(i,ny-3+js,k)   &
     &                - 2.* r(i,ny-2+js,k) )

       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! Z pass

      DO k = 3,nz-3+ks
       DO j = 1+js,ny-1
        DO i = 1+is,nx-1

         fs(i,j,k) = fs(i,j,k)                             &
     &             - f4dmph *   ((r(i,j,k+2) + r(i,j,k-2)) &
     &                      - 4.*(r(i,j,k+1) + r(i,j,k-1)) &
     &                      + 6.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 2nd order filter at the boundaries

      DO j = 1+js,ny-1
       DO i = 1+is,nx-1

        fs(i,j,2) = fs(i,j,2) &
     &   + f2dmph * (r(i,j,3) + r(i,j,1) - 2.* r(i,j,2))
        fs(i,j,nz-2+ks) = fs(i,j,nz-2+ks)               &
     &   + f2dmph * ( r(i,j,nz-1+ks) + r(i,j,nz-3+ks)  &
     &          - 2.* r(i,j,nz-2+ks) )

       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FILTER4

!------------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE FILTER6  ////////////////////
!
! FILTER6 performs a sixth-order filter on the passed variable
!------------------------------------------------------------------------------
! Created by Louis Wicker, April 6, 1988
! Latest update: 01-22-92
! Update notes:
!------------------------------------------------------------------------------
      SUBROUTINE FILTER6(s,sb,fs,r,svar,nx,ny,nz)

   USE GRID_MODULE
   USE PARAM_MODULE

      implicit none

! Passed variables

      integer nx, ny, nz
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real sb(nz)

   TYPE(VARIABLE) :: svar

   real :: r (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

! Local variables

      integer i, j, k, is, js, ks

      real            :: f6dmph, f4dmph, f2dmph
      real, parameter :: f6damp = 1.0/64.
      real            :: f1, sgrad
      

!------------------------------------------------------------------------------
! Set consts:  f2dmp should have same damping at 2 dx as f4dmp

      f6dmph = f6damp * damph
      f4dmph = f6dmph * 4.0
      f2dmph = f4dmph * 4.0
      
      is = svar%istag
      js = svar%jstag
      ks = svar%kstag

!------------------------------------------------------------------------------
! Special horizontal filter

!       f4dmph = f4damp * damph * (1.0 + 2.0*(float(k)/float(nz-1))**4)
!       f2dmph = f4dmph * 4.0

!------------------------------------------------------------------------------
! Create perturbation fields

      DO k = 1,nz
       DO j = 1,ny
        DO i = 1,nx

         r(i,j,k) = s(i,j,k) - sb(k)

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! X pass

      IF ( nx .gt. 2 ) THEN
      
      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1
        DO i = 4,nx-4+is

         f1 =        f6dmph *  ((r(i+3,j,k) + r(i-3,j,k))  &
     &                     - 6.*(r(i+2,j,k) + r(i-2,j,k))  &
     &                     +15.*(r(i+1,j,k) + r(i-1,j,k))  &
     &                     -20.* r(i,j,k))

        sgrad = r(i+1,j,k) - r(i-1,j,k)
        fs(i,j,k) = fs(i,j,k) + f1  ! *Max( 0.0, Sign(1.0,f1*sgrad) )
        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 4th/2nd order filter at the boundaries

      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1

        f1 =   f2dmph * (r(3,j,k) + r(1,j,k) - 2.* r(2,j,k))
        
        sgrad = r(3,j,k) - r(1,j,k)
        fs(2,j,k) = fs(2,j,k) + f1  ! *Max( 0.0, Sign(1.0,f1*sgrad) )

        f1 =       - f4dmph *   ((r(5,j,k) + r(1,j,k))  &
     &                      - 4.*(r(4,j,k) + r(2,j,k))  &
     &                      + 6.* r(3,j,k) )

        sgrad = r(4,j,k) - r(2,j,k)
        fs(3,j,k) = fs(3,j,k) + f1  ! *Max( 0.0, Sign(1.0,f1*sgrad) )

        f1 =            - f4dmph *   ((r(nx-1+is,j,k) + r(nx-5+is,j,k))  &
     &                           - 4.*(r(nx-2+is,j,k) + r(nx-4+is,j,k))  &
     &                           + 6.* r(nx-3+is,j,k) )

        sgrad = r(nx-2,j,k) - r(nx-4,j,k)
        fs(nx-3+is,j,k) = fs(nx-3+is,j,k) + f1  ! *Max( 0.0, Sign(1.0,f1*sgrad) )

        f1 =            + f2dmph * (r(nx-1+is,j,k)+r(nx-3+is,j,k)  &
     &                          -2.*r(nx-2+is,j,k))

        sgrad = r(nx-1,j,k) - r(nx-3,j,k)
        fs(nx-2+is,j,k) = fs(nx-2+is,j,k) + f1 ! *Max( 0.0, Sign(1.0,f1*sgrad) )

       ENDDO
      ENDDO
      
      ENDIF

!------------------------------------------------------------------------------
! Y pass

      IF ( ny .gt. 2 ) THEN
      
      DO k = 1+ks,nz-1
       DO j = 4,ny-4+js
        DO i = 1+is,nx-1

          fs(i,j,k) = fs(i,j,k)                            &
     &             + f6dmph *  ((r(i,j+3,k) + r(i,j-3,k))  &
     &                     - 6.*(r(i,j+2,k) + r(i,j-2,k))  &
     &                     +15.*(r(i,j+1,k) + r(i,j-1,k))  &
     &                     -20.* r(i,j,k))

        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 4th/2nd order filter at the boundaries

      DO k = 1+ks,nz-1
       DO i = 1+is,nx-1

        fs(i,2,k) = fs(i,2,k)     &
     &   + f2dmph * ( r(i,3,k)  + r(i,1,k)    - 2.* r(i,2,k) )

        fs(i,3,k) = fs(i,3,k)     &
     &             - f4dmph *   ((r(i,5,k) + r(i,1,k)) &
     &                      - 4.*(r(i,4,k) + r(i,2,k)) &
     &                      + 6.* r(i,3,k) )

        fs(i,ny-3+js,k) = fs(i,ny-3+js,k) &
     &                  - f4dmph *   ((r(i,ny-1+js,k) + r(i,ny-5+js,k)) &
     &                           - 4.*(r(i,ny-2+js,k) + r(i,ny-4+js,k)) &
     &                           + 6.* r(i,ny-3+js,k) )

        fs(i,ny-2+js,k) = fs(i,ny-2+js,k) &
     &   + f2dmph * ( r(i,ny-1+js,k) + r(i,ny-3+js,k)  &
     &                - 2.* r(i,ny-2+js,k) )

       ENDDO
      ENDDO
      
      ENDIF

!------------------------------------------------------------------------------
! Z pass

    IF ( .false. ) THEN
      
      DO k = 4,nz-4+ks
       DO j = 1+js,ny-1
        DO i = 1+is,nx-1

         fs(i,j,k) = fs(i,j,k)                            &
     &             + f6dmph *  ((r(i,j,k+3) + r(i,j,k-3)) &
     &                     - 6.*(r(i,j,k+2) + r(i,j,k-2)) &
     &                     +15.*(r(i,j,k+1) + r(i,j,k-1)) &
     &                     -20.* r(i,j,k))


        ENDDO
       ENDDO
      ENDDO

!------------------------------------------------------------------------------
! 4th/2nd order filter at the boundaries

      DO j = 1+js,ny-1
       DO i = 1+is,nx-1

        fs(i,j,2) = fs(i,j,2)     &
     &            + f2dmph * (r(i,j,3) + r(i,j,1) - 2.* r(i,j,2))

        fs(i,j,3) = fs(i,j,3)                          &
     &            - f4dmph *   ((r(i,j,5) + r(i,j,1))  &
     &                     - 4.*(r(i,j,4) + r(i,j,2))  &
     &                     + 6.* r(i,j,3) )

        fs(i,j,nz-3+ks) = fs(i,j,nz-3+ks)                                &
     &                  - f4dmph *   ((r(i,j,nz-1+ks) + r(i,j,nz-5+ks))  &
     &                           - 4.*(r(i,j,nz-2+ks) + r(i,j,nz-4+ks))  &
     &                           + 6.* r(i,j,nz-3+ks) )

        fs(i,j,nz-2+ks) = fs(i,j,nz-2+ks)                           &
     &                  + f2dmph * (r(i,j,nz-1+ks) + r(i,j,nz-3+ks) &
     &                        - 2.* r(i,j,nz-2+ks))

       ENDDO
      ENDDO
      
     ENDIF 

      RETURN
      END SUBROUTINE FILTER6
