!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\     MODULE FORCE      ////////////////////
!     
!-----------------------------------------------------------------------------

  MODULE FORCE_MODULE
  
!-----------------------------------------------------------------------------

  implicit none

  integer, parameter :: nfmax = 10 ! Max number of forcing areas

  integer ::  iwshap(nfmax) = 1 ! for W forcing
  integer ::  ivshap(nfmax) = 1
  integer ::  iushap(nfmax) = 1
  integer ::  iqshap(nfmax) = 1
  integer ::  ichgshap(nfmax) = 1
  integer ::  itshap(nfmax) = 1
  integer ::  isshap(nfmax) = 1 ! for QV/SS forcing
  integer ::  isslcl(nfmax) = 1 ! for QV/SS forcing (use LCL as lower boundary of forcing)
  
  integer ::  iforce = 0
  
  integer ::  iforcetyp = 1
  integer ::  nwfor     = 0  ! number of W forcing regions (max of nfmax)
  integer ::  nssfor    = 0  ! number of QV/SS forcing regions (max of nfmax)
  integer ::  nccnfor   = 0  ! number of CCN forcing regions (max of nfmax)
  integer ::  nchgfor   = 0  ! number of charge forcing regions (max of nfmax)
  integer ::  iwforce   = 0  ! forcing flag for W
  integer ::  issforce  = 0  ! forcing for QV (supersaturation)
  integer ::  iqrforce  = 0  ! forcing for QR
  integer ::  iqcforce  = 0  ! forcing for QC
  integer ::  iccnufforce  = 0  ! forcing for CCNUF
  integer ::  ichgforce  = 0  ! forcing for charge
  integer ::  itopforce = 0
  integer ::  itopforcetype = 1 ! 1 = rain, 2 = graupel, etc., 3 = graupel as 100% rain
  integer ::  itopforcerainopt = 1 ! rain options (binforce) 1 = default, 2 = MP using topforcen0 and topforcerainrate
  real    ::  forcerandfac = 0.0 ! factor for random perturbations to forcing (using ranarray2d). Must use bubbletype 2 or 3
  real    ::  wnaylor = 10.0 ! target W for Naylor/Gilmore forcing
  real    ::  topforceq = 0.0
  real    ::  topforceN = 0.0
  real    ::  topforceX = 0.0 ! xlocation
  real    ::  topforceY = 0.0 ! ylocation
  real    ::  topforceR = 4000.0 ! radius of forcing region (meters)
  real    ::  topforcedia = 1.5e-3 ! mean particle diameter (meters)
  real    ::  topforcen0  = 8.0e6 ! Initial intercept parameter of distribution (-1 to use topforcedia instead)
  real    ::  topforcerainrate  = 54.0 ! rainrate to use with MP rain distribution
  real    ::  topforcealpha = 0    ! Initial shape parameter

  real    ::  qrfmeso(nfmax) = 0.0 ! value of QR for forcing
  real    ::  qcfmeso(nfmax) = 0.0 ! value of QC for forcing
  real    ::  qhfmeso(nfmax) = 0.0 ! value of QH for forcing
  real    ::  ssfmeso(nfmax) = 0.0 ! value of QV/SS for forcing
  real    ::  xssfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  yssfcen(nfmax) = 0.0 ! y-center of forcing

  real    ::  rssfrad(nfmax) = 0.0 ! radius of cylindrical forcing region
  
  real    ::  wfmeso(nfmax) = 0.0 ! value of dw/dt for accel forcing or 'alpha' for Naylor forcing
  real    ::  ufmeso(nfmax) = 0.1 ! value of U for forcing 
  real    ::  umaxmeso(nfmax) = 60.0 ! max value of U when forcing is turned off
  real    ::  xwfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  ywfcen(nfmax) = 0.0 ! y-center of forcing
  real    ::  zwfcen(nfmax) = 0.0 ! z-center of forcing

  real    ::  xccnfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  yccnfcen(nfmax) = 0.0 ! y-center of forcing
  real    ::  zccnfcen(nfmax) = 0.0 ! z-center of forcing

  real    ::  chgratefor(nfmax) = 0.0 ! value of d(charge)/dt
  real    ::  xchgfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  ychgfcen(nfmax) = 0.0 ! y-center of forcing
  real    ::  zchgfcen(nfmax) = 0.0 ! z-center of forcing

  real    ::  xchgcldfcen(nfmax) = 0.0 ! x-center of 'cloud' around the charge
  real    ::  ychgcldfcen(nfmax) = 0.0 ! y-center
  real    ::  zchgcldfcen(nfmax) = 0.0 ! z-center

  real    ::  xwfcennew(nfmax) = 0.0 ! x-center of forcing
  real    ::  ywfcennew(nfmax) = 0.0 ! y-center of forcing
  real    ::  zwfcennew(nfmax) = 0.0 ! z-center of forcing

  real    ::  xwfrad(nfmax) = 0.0 ! radii of ellipsoidal forcing region
  real    ::  ywfrad(nfmax) = 0.0 
  real    ::  zwfrad(nfmax) = 0.0 

  real    ::  xwfmov(nfmax) = 0.0 ! movement rate of forcing region
  real    ::  ywfmov(nfmax) = 0.0 
  real    ::  zwfmov(nfmax) = 0.0 

  real    ::  xccnfrad(nfmax) = 0.0 ! radii of ellipsoidal forcing region
  real    ::  yccnfrad(nfmax) = 0.0
  real    ::  zccnfrad(nfmax) = 0.0

  real    ::  xchgfrad(nfmax) = 0.0 ! radii of cylindrical forcing region
  real    ::  ychgfrad(nfmax) = -1.0
  real    ::  zchgfrad(nfmax) = 0.0

  real    ::  xchgcldfrad(nfmax) = 0.0 ! radii of cylindrical 'cloud' around charge forcing region
  real    ::  ychgcldfrad(nfmax) = -1.0
  real    ::  zchgcldfrad(nfmax) = 0.0


  integer ::  twstrt(nfmax) = 0   ! starting time of W forcing
  integer ::  twstop(nfmax) = 0   ! ending time of W forcing

  integer ::  tccnstrt(nfmax) = 0   ! starting time of CCNuf forcing
  integer ::  tccnstop(nfmax) = 0   ! ending time of CCNuf forcing

  integer ::  tsstrt(nfmax) = 0   ! starting time of QV/SS forcing
  integer ::  tsstop(nfmax) = 0   ! ending time of QV/SS forcing

  real ::  tsslow(nfmax) = 253.15   ! low (i.e., top) temperature of QV/SS forcing
  real ::  tsshigh(nfmax) = 285.15   ! high (i.e., bottom) temperature of QV/SS forcing
  
  integer :: nphys = 1   ! number of substeps for microphysics
  integer :: igamrain = 1  ! 1 = gamma of volume, 2 = gamma of diameter
  integer :: igamsnow = 1  ! 1 = gamma of volume, 2 = gamma of diameter
  
  real :: chargeperparticle = 100.e-15
  


! Parker (2008) low level cooling
  integer              :: low_level_cooling_flag = 0 ! 0 = off; 1 = on
  real                 :: low_level_cooling_depth = 1000. ! meters
  real                 :: low_level_cooling_rate  = 3.0 ! degrees per hour
  integer              :: low_level_cooling_start = 10800
  integer              :: low_level_cooling_end   = 21600

 
  real, allocatable :: tt13sl(:,:) ! U flux
  real, allocatable :: tt23sl(:,:) ! V flux
  real, allocatable :: td13sl(:,:) ! deformation term
  real, allocatable :: td23sl(:,:) ! deformation term
  real, allocatable :: twt3sl(:,:) ! TKE buoyancy term
  real, allocatable :: thf3sl(:,:) ! TH flux
  real, allocatable :: tmf3sl(:,:) ! QV flux
  real, allocatable :: tdh3sl(:,:) ! not used


  CONTAINS
  
!-----------------------------------------------------------------------------
   SUBROUTINE FORCE_INIT
   RETURN
   END SUBROUTINE FORCE_INIT
   
!-----------------------------------------------------------------------------
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE FORCE      ////////////////////
!     general forcing functions
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE FORCE()
   
!-----------------------------------------------------------------------------
   implicit none
   
   
   RETURN
   
   END SUBROUTINE FORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE WFORCE      ////////////////////
!     general forcing functions
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE WFORCE(nx,ny,nz,dt,w,wt,fw,u,ut,fu,t0,             &
                     gxt,gyt,gzt,loop,time,time_real,    &
                     uinit,vinit,ugrid,vgrid,ranarr)

   USE PARAM_MODULE, only: RKSCHEME,ng,pii,luno,bcy
#ifdef MPI
   USE COMMASMPI_MODULE
#endif

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz
   real    :: dt
   real    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: wt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: ut(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng)
   real, intent(in) :: ranarr(-ng+1:nx+ng,-ng+1:ny+ng)
   real    :: ugrid, vgrid
   integer :: loop
   integer :: time
   real    :: time_real
   
   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius, tmp
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: gamma

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
#endif

!-----------------------------------------------------------------------------

   work_to_do = .false.
   idofor(:)  = 0
   
   
   DO n = 1,nwfor


    IF ( time_real .ge. twstrt(n) .and. time_real .le. twstop(n) ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ENDIF

! set movement to mean layer wind
      IF ( xwfmov(n) .le. -999. ) THEN
       wtot = 0.0
       u1 = 0.0
       v1 = 0.0

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
    
       do k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        z = gzt(k,1)

         radius = sqrt(((z-zwfcen(n))/zwfrad(n))**2)
      
       IF ( radius .lt. 1.0 ) THEN
       
        IF ( xwfmov(n) .eq. -999. ) THEN       ! weighted mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         wtot = wtot + (cos(pii*radius)+1.)/2./gzt(k,3)
        ELSEIF ( xwfmov(n) .eq. -1000. ) THEN  ! mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
        ENDIF
       ELSEIF ( xwfmov(n) .eq. -1001. ) THEN ! mean wind through depth of domain
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 10.e3 ) THEN ! mean wind through 10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 6.e3 ) THEN ! mean wind through 6 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1003. .and. z .gt. 5.e3 .and. z .lt. 10.e3 ) THEN ! mean wind 5-10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ENDIF
       
      ENDDO

       xwfmov(n) = u1/wtot
       ywfmov(n) = v1/wtot
       
       write(luno,*) 'wforce: n, xwfmov, ywfmov = ',n,xwfmov(n),ywfmov(n)
      
      ENDIF

       IF ( Abs(xwfmov(n)) .lt. 900. ) THEN
        xwfcennew(n) = xwfcen(n) + Max(0.0, time_real - twstrt(n))*xwfmov(n)
!        print*, 'Update xwfcen(',n,') = ',xwfcen(n),dt,xwfmov(n)
       ELSE
         xwfcennew(n) = xwfcen(n) 
       ENDIF
       IF ( Abs(ywfmov(n)) .lt. 900. ) THEN
        ywfcennew(n) = ywfcen(n) + Max(0.0, time_real - twstrt(n))*ywfmov(n)
       ELSE
        ywfcennew(n) = ywfcen(n) 
       ENDIF
       IF ( Abs(zwfmov(n)) .lt. 900. ) THEN
        zwfcennew(n) = zwfcen(n) ! + Max(0.0, time_real - twstrt(n))*zwfmov(n)
       ELSE
        zwfcennew(n) = zwfcen(n) 
       ENDIF
      
   ENDDO
   

!    IF ( .not. work_to_do )  RETURN
    IF ( work_to_do ) THEN ! RETURN
    
    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
   DO n = 1,nwfor
     IF ( abs(iwshap(n)) .eq. 1 ) THEN
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((y-ywfcennew(n))/ywfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ELSE
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
!       IF ( iwshap(n) .lt. 0 .and. dist .le. 1.0 ) THEN
        tmp = 0.0
        IF ( forcerandfac /= 0.0 ) tmp = forcerandfac*ranarr(i,j)
        
        IF ( iwshap(n) .eq. 1 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
         IF ( wfmeso(n) .gt. 0.0 ) THEN
           t0(i,j,k) = Max(t0(i,j,k), (wfmeso(n)+tmp)*(cos(pii*radius)+1.)/2.)
         ELSE
           t0(i,j,k) = (wfmeso(n)+tmp)*(cos(pii*radius)+1.)/2.
         ENDIF
        ELSEIF ( iwshap(n) .eq. -1 .and. w(i,j,k) .ge. -0.1 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
          t0(i,j,k) = Max(t0(i,j,k), (wfmeso(n)+tmp)*(cos(pii*radius)+1.)/2.)
        ENDIF
!       ENDIF
       ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

     ELSEIF ( abs(iwshap(n)) .eq. 2 ) THEN ! nudge toward 10 m/s (Naylor et al. 2012, MWR)
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((y-ywfcennew(n))/ywfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ELSE
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
!       IF ( iwshap(n) .lt. 0 .and. dist .le. 1.0 ) THEN
        IF ( iwshap(n) .eq. 2 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
           t0(i,j,k) = wfmeso(n)*Max(0.0, wnaylor*(cos(pii*radius)+1.)/2. - wt(i,j,k) )
        ELSEIF ( iwshap(n) .eq. -2  ) THEN
           t0(i,j,k) =(wfmeso(n))*Max(0.0, wnaylor*(cos(pii*radius)+1.)/2. - w(i,j,k) )
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
!          t0(i,j,k) = Max(t0(i,j,k), wfmeso(n)*(cos(pii*radius)+1.)/2.)
        ENDIF
!       ENDIF
       ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF
 
      ELSEIF ( abs(iwshap(n)) .eq. 3 ) THEN ! nudge U (for workshop squall line case)
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

      gamma = 1.0
      IF (time_real >= twstop(n)-300.) THEN
        gamma = 1.0+(0.0-1.0)*( time_real - (twstop(n)-300.) )/(300.)
      ENDIF

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
!         IF ( ny .gt. 2 ) THEN
!         radius = sqrt(((x-xwfcen(n))/xwfrad(n))**2  &
!                      +((y-ywfcen(n))/ywfrad(n))**2  &
!                      +((z-zwfcen(n))/zwfrad(n))**2)
!         ELSE
! 2D forcing -- assumes periodic domain or channel simulation
          radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
!         ENDIF

!        IF ( Abs(x-xwfcennew(n)) < xwfrad(n) .and. Abs(z - zwfcen(n)) < zwfrad(n) ) THEN
        IF ( ( bcy == 2 .and. Abs(x-xwfcennew(n)) < xwfrad(n) .and.  z < zwfrad(n) ) .or. &
     &       ( Abs(y-ywfcennew(n))/ywfrad(n) <= 1. .and. Abs(x-xwfcennew(n)) < xwfrad(n) .and.  z < zwfrad(n) )   ) THEN
!        IF ( radius < 1. ) THEN
!          uten1(i,j,k)=uten1(i,j,k)+0.10*gamma           &
!                       *cos(0.5*pi*(xh(i)-0.0)/10000.0)  &
!                       *((cosh(2.5*(zh(i,j,k)-0.0)/10000.0))**(-2))

         IF ( abs( u(i,j,k) ) <= umaxmeso(n) ) THEN
         fu(i,j,k) = fu(i,j,k) + ufmeso(n)*gamma           &
                       *cos(0.5*pii*(x-xwfcennew(n))/xwfrad(n))  &
                       *((cosh(2.5*z/zwfrad(n)))**(-2))
!                       *((cosh(2.5*(z-zwfcen(n))/zwfrad(n)))**(-2))
         ENDIF
        ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

     ENDIF ! iwshap
   
   ENDDO

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
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
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
        fw(i,j,k) = fw(i,j,k) + t0(i,j,k)
        
      ENDDO
     ENDDO
    ENDDO
    
    ENDIF ! work_to_do

!   DO n = 1,nwfor
!    IF ( loop .eq. RKSCHEME ) THEN
!       IF ( Abs(xwfmov(n)) .lt. 900. ) THEN
!        xwfcen(n) = xwfcen(n) + dt*xwfmov(n)
!!        print*, 'Update xwfcen(',n,') = ',xwfcen(n),dt,xwfmov(n)
!       ENDIF
!       IF ( Abs(ywfmov(n)) .lt. 900. ) THEN
!        ywfcen(n) = ywfcen(n) + dt*ywfmov(n)
!       ENDIF
!       IF ( Abs(zwfmov(n)) .lt. 900. ) THEN
!        zwfcen(n) = zwfcen(n) + dt*zwfmov(n)
!       ENDIF
!     ENDIF
!   ENDDO
   
   RETURN
   
   END SUBROUTINE WFORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE SSFORCE      ////////////////////
!     Water vapor forcing based on Fierro et al. 2012.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE SSFORCE(nx,ny,nz,ns,dt,w,st,pt,t0,piinit1d,sbase,        &
                     gxt,gyt,gzt,time,time_real,loop,den,   &
                     uinit,vinit,ugrid,vgrid)

   USE PARAM_MODULE, only: ng, pii, RKSCHEME, luno
   USE COMMASMPI_MODULE
   USE INDEX_MODULE

!-----------------------------------------------------------------------------
   implicit none

#ifdef MPI
   include 'mpif.h'
#endif
   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   integer :: time
   real    :: time_real
   integer :: loop
   real    :: uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng)
   real    :: ugrid, vgrid
   real    :: den(-ng+1:nz+ng,2)

   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: temlo,temhi
   real    :: temp, qvs
   real    :: fac, diff,xmas,xmas0
   double precision :: qvadded,qcadded,qradded,qhadded,dv
   double precision  mpitotindp(15), mpitotoutdp(15)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   real :: temg1d(-ng+1:nz+ng)
   real :: hgtlcl(nx,ny)
!-----------------------------------------------------------------------------


!      write(*,*) 'In SSforce'

   work_to_do = .false.
   idofor(:)  = 0
   
   temlo = 253.15
   temhi = 273.15
   
   qvadded = 0.d0
   qradded = 0.d0
   qcadded = 0.d0
   qhadded = 0.d0
   
   DO n = 1,nssfor

! set movement to mean layer wind
      IF ( xwfmov(n) .le. -999. ) THEN
       wtot = 0.0
       u1 = 0.0
       v1 = 0.0

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
    
       do k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        z = gzt(k,1)

         radius = sqrt(((z-zwfcen(n))/zwfrad(n))**2)
      
       IF ( radius .lt. 1.0 ) THEN
       
        IF ( xwfmov(n) .eq. -999. ) THEN       ! weighted mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         wtot = wtot + (cos(pii*radius)+1.)/2./gzt(k,3)
        ELSEIF ( xwfmov(n) .eq. -1000. ) THEN  ! mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
        ENDIF
       ELSEIF ( xwfmov(n) .eq. -1001. ) THEN ! mean wind through depth of domain
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 10.e3 ) THEN ! mean wind through 10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 6.e3 ) THEN ! mean wind through 6 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1003. .and. z .gt. 5.e3 .and. z .lt. 10.e3 ) THEN ! mean wind 5-10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ENDIF
       
      ENDDO

       xwfmov(n) = u1/wtot
       ywfmov(n) = v1/wtot
       
       write(luno,*) 'ssforce: n, xwfmov, ywfmov = ',n,xwfmov(n),ywfmov(n)
      
      ENDIF
      

    IF ( time_real .ge. tsstrt(n) .and. time_real .le. tsstop(n) ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ENDIF

   ENDDO

    IF ( work_to_do ) THEN
    
    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg-1

!       jyb = -ng+1
!       jye = jtile+ng
       jyb = 1
       jye = jtile
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

!       ixb = -ng+1
!       ixe = itile+ng
       ixb = 1
       ixe = itile
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

! compute base state temperature profile
   DO k = kzb,kze
    temg1d(k) = sbase(k,1)*piinit1d(k)
   ENDDO

   CALL adaslcl(nx,ny,nz,gzt,pt,piinit1d,st(-ng+1,-ng+1,-ng+1,lt),st(-ng+1,-ng+1,-ng+1,lv),hgtlcl)

   DO n = 1,nssfor
     IF ( abs(isshap(n)) .eq. 1 ) THEN
      
      IF ( idofor(n) .eq. 1 ) THEN

!       write(*,*) 'SS forcing at time = ',time_real
       do k = kzb,kze
        z = gzt(k,1)
        
!        IF ( temg1d(k) >= tsslow(n) .and. temg1d(k) <= tsshigh(n) ) THEN
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
         
          dv = 1.d0/(gxt(i,3)*gyt(j,3)*gzt(k,3))

         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xssfcen(n))/rssfrad(n))**2  &
                      +((y-yssfcen(n))/rssfrad(n))**2)
         ELSE
         radius = sqrt(((x-xssfcen(n))/rssfrad(n))**2 )
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
         IF ( ssfmeso(n) .gt. 0.0 ) THEN
! supersat calculation...
        IF ( Abs(isshap(n)) .eq. 1 ) THEN
         fac = (cos(pii*radius)+1.)/2.
          IF ( isshap(n) == -1 .and.   w(i,j,k) .le. -0.3 ) THEN
            ! turn off forcing in downdrafts (note that w is at w points, but not worrying about that just yet)
            fac = 0.0
          ENDIF
        ELSE
         fac = 1.0
        ENDIF

         IF ( fac > 0.0 ) THEN
!           IF ( temg1d(k) >= tsslow(n) .and. temg1d(k) <= tsshigh(n) ) THEN
           IF ( temg1d(k) >= tsslow(n) .and. ( ( isslcl(n) > 0 .and. gzt(k,1) >= hgtlcl(i,j) ) .or. &
     &            ( isslcl(n) <= 0 .and. temg1d(k) <= tsshigh(n)) )  ) THEN
             temp  = (piinit1d(k) + pt(i,j,k) )*st(i,j,k,lt)
             qvs   = 380.*exp(17.27*(temp - 273.)/(temp - 36.))/(1.e5*(piinit1d(k) + pt(i,j,k))**3.509)
         
!         write(*,*) 'SSForce: i,k,qv,qvs,qvforce = ',i,k,st(i,j,k,lv),qvs,  &
!     &             st(i,j,k,lv) + (ssfmeso(n)*qvs - st(i,j,k,lv))*(cos(pii*radius)+1.)/2.
!         write(*,*) 'temp, pi1d, pt, lt = ',temp,piinit1d(k),pt(i,j,k),st(i,j,k,lt)
         
             qvadded = qvadded + Max(0.0, (ssfmeso(n)*qvs - st(i,j,k,lv)))*fac*dv*den(k,2)
             st(i,j,k,lv) = st(i,j,k,lv) + Max(0.0, (ssfmeso(n)*qvs - st(i,j,k,lv)))*fac  ! *Min(1., dt/30.)
           ENDIF
           
!         IF ( qrfmeso(n) > 0. .and. iqrforce > 0 .and. temg1d(k) >= 263.0 .and. temg1d(k) <= 283.0 ) THEN
         IF ( qrfmeso(n) > 0. .and. iqrforce > 0 .and. temg1d(k) >= 263.0 .and.   &
           (  gzt(k+1,1) >= hgtlcl(i,j) .or. ( isslcl(n) <= 0 .and. temg1d(k) <= tsshigh(n)) ) ) THEN
           diff = Max( 0., (qrfmeso(n)*fac - st(i,j,k,lr))/60.*dt )
           IF ( diff > 0. ) THEN
             st(i,j,k,lr) = st(i,j,k,lr) + diff
             qradded = qradded + diff*dv*den(k,2)
             IF ( lnr > 0 ) THEN  ! note that number conc. is still a "mixing ratio" at this point
               xmas0 = 1000.*(pii/6.)*(300.e-6)**3 ! set minimum size for new drops
               IF ( st(i,j,k,lnr) > 0.1 ) THEN
                 xmas = st(i,j,k,lr)/st(i,j,k,lnr)
               ELSE
                 xmas = 0.0
               ENDIF
               xmas = Max( xmas, xmas0 )
               st(i,j,k,lnr) = st(i,j,k,lnr) + diff/xmas
             ENDIF
           ENDIF
         ENDIF

         IF ( qcfmeso(n) > 0. .and. iqcforce > 0 .and. temg1d(k) >= 233.0 .and.   &
           (  gzt(k+1,1) >= hgtlcl(i,j) .or. ( isslcl(n) <= 0 .and. temg1d(k) <= tsshigh(n)) ) ) THEN
           diff = Max( 0., (qcfmeso(n)*fac - st(i,j,k,lc))/60.*dt )
           IF ( diff > 0. ) THEN
             st(i,j,k,lc) = st(i,j,k,lc) + diff
             qcadded = qcadded + diff*dv*den(k,2)
             IF ( lnc > 0 ) THEN  ! note that number conc. is still a "mixing ratio" at this point
               xmas0 = 1000.*(pii/6.)*(10.e-6)**3 ! set minimum size for new drops
               IF ( st(i,j,k,lnc) > 0.1 ) THEN
                 xmas = st(i,j,k,lc)/st(i,j,k,lnc)
               ELSE
                 xmas = 0.0
               ENDIF
               xmas = Max( xmas, xmas0 )
               st(i,j,k,lnc) = st(i,j,k,lnc) + diff/xmas
             ENDIF
           ENDIF
         ENDIF

         IF ( qhfmeso(n) > 0. .and. iqrforce > 0 .and. temg1d(k) >= 253.0 .and. temg1d(k) <= 273.0 ) THEN
           diff = Max( 0., (qhfmeso(n)*fac - st(i,j,k,lh))/60.*dt )
           IF ( diff > 0. ) THEN
             st(i,j,k,lh) = st(i,j,k,lh) + diff
             qhadded = qhadded + diff*dv*den(k,2)
             IF ( lnh > 0 ) THEN  ! note that number conc. is still a "mixing ratio" at this point
               xmas0 = 800.*(pii/6.)*(300.e-6)**3 ! set minimum size for new drops
               IF ( st(i,j,k,lnr) > 0.1 ) THEN
                 xmas = st(i,j,k,lh)/st(i,j,k,lnh)
               ELSE
                 xmas = 0.0
               ENDIF
               xmas = Max( xmas, xmas0 )
               st(i,j,k,lnh) = st(i,j,k,lnh) + diff/xmas
               st(i,j,k,lvh) = st(i,j,k,lvh) + diff/800.
             ENDIF
           ENDIF
         ENDIF
         
         ENDIF ! fac > 0.0
         
!         st(i,j,k,lv) = Max( st(i,j,k,lv), ssfmeso(n)*qvs*(cos(pii*radius)+1.)/2.)
! assumes A,B, and C constant, no graupel
!           supsat = max(q(i,k),qmin)-qs(i,k,1)
!           qdiff = max(q(i,k),qmin)-qs(i,k,1)! -0.001*qs(i,k,1)
!
!
! !          IF (1.eq.0) THEN ! CTRL
!
!
!            angle(i,k)=gridlight2(i,k)/10. !
!            angle2(i,k)=((1000*qrs(i,k,3))**1.5)/1 !
!
!           if (gridlight2(i,k).ge.1) then
!           if (qrs(i,k,3).lt.0.005) then ! if Qgraupel greater than 5 g/kg no action. ! based on MATLAB plot
!           if (q(i,k).lt.0.95*qs(i,k,1)) then ! if ambient Qvapor already  >= 95% of env qsat no action
!               if (t(i,k).lt.273.15.and.t(i,k).gt.253.15) then  ! increase Qv between 0 and -20C only
!               q(i,k)= 0.95*qs(i,k,1)+0.1*qs(i,k,1)*tanh(0.25*angle(i,k))*(1.-tanh(0.25*angle2(i,k))) !
!               endif
!           endif
!           endif
!           endif
!         

!           t0(i,j,k) = Max(t0(i,j,k), wfmeso(n)*(cos(pii*radius)+1.)/2.)
         ENDIF
!        ENDIF
!       ENDIF
       ENDIF

         ENDDO
        ENDDO
!        ENDIF
       ENDDO
      
      ENDIF
      
     ENDIF ! iwshap

   
   ENDDO

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = qvadded
       mpitotindp(2)  = qradded

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
#endif

      IF ( my_rank == 0 ) THEN
#ifdef MPI
       qvadded = mpitotoutdp(1)
       qradded = mpitotoutdp(2)
#endif
        write(luno,*) 'qvadded,qradded (Mg)= ',qvadded*1.e-3,qradded*1.e-3
       ENDIF

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
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
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
!        fw(i,j,k) = fw(i,j,k) + t0(i,j,k)
        
      ENDDO
     ENDDO
    ENDDO
    
    ENDIF ! work_to_do

   DO n = 1,nwfor
    IF ( loop .eq. RKSCHEME ) THEN
       IF ( Abs(xwfmov(n)) .lt. 900. ) THEN
        xssfcen(n) = xssfcen(n) + dt*xwfmov(n)
!        print*, 'Update xwfcen(',n,') = ',xwfcen(n),dt,xwfmov(n)
       ENDIF
       IF ( Abs(ywfmov(n)) .lt. 900. ) THEN
        yssfcen(n) = yssfcen(n) + dt*ywfmov(n)
       ENDIF
!       IF ( Abs(zwfmov(n)) .lt. 900. ) THEN
!        zssfcen(n) = zssfcen(n) + dt*zwfmov(n)
!       ENDIF
     ENDIF
   ENDDO


   RETURN
   
   END SUBROUTINE SSFORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE QFORCE      ////////////////////
!     Hydrometeor forcing at level nz-1 for sedimentation tests
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE QFORCE(nx,ny,nz,ns,dt,st,pt,t0,piinit1d,sbase,den,    &
                     gxt,gyt,gzt,time,time_real)

   USE PARAM_MODULE, only: ng, pii
   USE COMMASMPI_MODULE
   USE INDEX_MODULE
   USE MICRO_MODULE, only: rho_qh,rho_qhl,imurain

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   integer :: time
   real    :: time_real
   real    :: den(-ng+1:nz+ng,2)
   
   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: temlo,temhi
   real    :: temp, tmpq, tmpn, tmpz
   real    :: fac, diff,xmas,xmas0
   integer :: lx, lnx, lzx, lxw, lvx
   real    :: xden, alp, g1
   real    :: xpi = 3.14159
   real    :: xv,dtmp
   integer :: ifirst = 0 ! setting value here saves the value going forward

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   real :: temg1d(-ng+1:nz+ng)
!-----------------------------------------------------------------------------


!      write(0,*) 'In Qforce'

   work_to_do = .false.
   idofor(:)  = 0
   
   temlo = 253.15
   temhi = 273.15
   
   n = 1
   nssfor = 1
! only one region for now, and do it for all time.
   
   DO n = 1,nssfor

!
!    IF ( time_real .ge. tsstrt(n) .and. time_real .le. tsstop(n) ) THEN
      IF ( ifirst == 0 ) THEN
        work_to_do = .true.
      ENDIF
      idofor(n) = 1
!    ENDIF

   ENDDO

    IF ( work_to_do ) THEN
    
!    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
!       write(0,*) 'QFORCE: 1'
       
       IF ( itopforce <= 1 ) THEN
         kzb = nz-2
       ELSE
         ifirst = 1
         kzb = nz-2-itopforce
       ENDIF
       kze = nz-2
!       if (kzbeg .eq. nzbeg) kzb = 1
!       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

   lxw = 0
   lvx = 0
   
   IF ( itopforcetype == 1 .or. itopforcetype == 10 ) THEN ! rain
     lx = lr
     lnx = lnr
     lzx = lzr
     xden = 1000.
     alp = rnu
   ELSEIF ( itopforcetype == 2 .or. itopforcetype == 20 ) THEN ! graupel
     lx = lh
     lnx = lnh
     lzx = lzh
     lvx = lvh
     xden = rho_qh
     alp = alphah
   !  write(0,*) 'QFORCE Graupel: xden,alp, q, diam',xden,alp,topforceq,topforcedia
     IF ( itopforcetype == 20 .and. time_real <= dt ) THEN
       kzb = 1
     ENDIF
   ELSEIF ( itopforcetype == 3 ) THEN ! mixed-phase graupel
     lx = lh
     lxw = lhw
     lnx = lnh
     lzx = lzh
     xden = 1000.
     lvx = lvh
     alp = alphah
   ELSEIF ( itopforcetype == 4 ) THEN ! hail
     lx = lhl
     lnx = lnhl
     lzx = lzhl
     lvx = lvhl
     lxw = 0 ! lhlw
     xden = rho_qhl
     alp = alphahl
   ELSEIF ( itopforcetype == 5 ) THEN ! snow
     lx = ls
     lnx = lns
     lzx = 0
     xden = 100.
     alp = snu
   ELSEIF ( itopforcetype == 6 ) THEN ! frozen drops
     lx = lf
     lnx = lnf
     lzx = lzf
     lvx = lvf
     lxw = 0 ! lhlw
     xden = rho_qh
     alp = alphah
   ELSEIF ( itopforcetype == 7 ) THEN ! cloud ice
     lx = li
     lnx = lni
     lzx = 0
     xden = 900.
     alp = cinu
   ELSE
     write(0,*) 'Bad value for itopforcetype = ',itopforcetype
     call commasmpi_abort()
   ENDIF

   DO n = 1,nssfor
     IF ( abs(iqshap(n)) >= 1 ) THEN
      
      IF ( idofor(n) .eq. 1 ) THEN

!       write(*,*) 'SS forcing at time = ',time_real
       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)

         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-topforceX)/topforceR)**2  &
                      +((y-topforceY)/topforceR)**2)
         ELSE
         radius = sqrt(((x-topforceX)/topforceR)**2 )
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
        IF ( topforceq .gt. 0.0 ) THEN
! supersat calculation...
        IF ( iqshap(n) .eq. 1 .or. iqshap(n) .eq. 2) THEN
         fac = (cos(xpi*radius)+1.)/2.
        ELSE
         fac = 1.0
        ENDIF

         tmpq = topforceq*fac
         
         
         st(i,j,k,lx) =  tmpq
         
         IF ( lxw > 0 ) st(i,j,k,lxw) = st(i,j,k,lx)
         
         IF ( lvx > 0 ) THEN
           st(i,j,k,lvx) = tmpq/xden ! leave as mixing ratio here.
         ENDIF
         
         IF ( lnx > 0 ) THEN
           IF ( tmpq > 0. ) THEN
            IF ( topforceN > 0 ) THEN
              tmpn = topforceN*(fac + 1.)*0.5
              st(i,j,k,lnx) =  tmpn/den(k,1)
            ELSE
              IF ( iqshap(n) .eq. 1 ) THEN
               temp = topforcedia*(fac + 1.)*0.5
              ELSE
               temp = topforcedia
              ENDIF
              dtmp = temp
              xv = xpi/6.*temp**3
              tmpn = tmpq*den(k,1)/(xden*xv)
              st(i,j,k,lnx) =  tmpn/den(k,1)
            ENDIF
           ELSE
           st(i,j,k,lnx) = 0.0
           ENDIF
         ENDIF
         
         IF ( lzx > 0 ) THEN ! reflectivity moment
            IF ( itopforcetype == 1 .and. imurain == 3 ) THEN ! rain
             temp = den(k,1)*tmpq/(1000.*Max(1.0e-9,tmpn))
             g1 = 36.*(alp+2.0)/((alp+1.0)*xpi**2)
             tmpz = g1*den(k,1)**2*(tmpq)*tmpq/(xden**2*tmpn)
            ELSEIF ( itopforcetype >= 2 .or. imurain == 1 ) THEN
            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/ &
     &            ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            tmpz = g1*den(k,1)**2*(6*tmpq)**2/(tmpn*(xpi*xden)**2)
            ENDIF
            st(i,j,k,lzx) = tmpz/den(k,1) ! convert to mixing ratio
         
         ENDIF
         
         IF ( tmpq > 0.0 .and. i == 3 ) THEN
 !        write(0,*) 'QFORCE Graupel: q,n,z', st(i,j,k,lx), st(i,j,k,lnx), st(i,j,k,lzx) 


!         write(0,*) 'i,xden = ',i,lvx,xden,st(i,j,k,lvx), den(k,1)*st(i,j,k,lx)/(st(i,j,k,lvx)*den(k,1))
!         xden = den(k,1)*st(i,j,k,lx)/(st(i,j,k,lvx)*den(k,1))
!         xv = den(k,1)*st(i,j,k,lx)/(xden*st(i,j,k,lnx)*den(k,1))
         
!         write(0,*) 'xv,dia = ',xv,dtmp,(6.0*xv/xpi)**(1./3.)
         ENDIF
         
         ENDIF ! topforceq
       ENDIF ! radius

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF
      
     ENDIF ! iwshap
   
   ENDDO

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
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
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
!        fw(i,j,k) = fw(i,j,k) + t0(i,j,k)
        
      ENDDO
     ENDDO
    ENDDO
    
    ENDIF ! work_to_do


   RETURN
   
   END SUBROUTINE QFORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE BINFORCE   ////////////////////
!     Hydrometeor forcing at level nz-1 for sedimentation tests -- Takahashi bin micro
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE BINFORCE(nx,ny,nz,ns,dt,st,pt,t0,piinit1d,sbase,den,    &
                     gxt,gyt,gzt,time,time_real, onedoutput)

   USE PARAM_MODULE, only: ng, pii
   USE COMMASMPI_MODULE
   USE INDEX_MODULE, only: lr,li,lh,lhl,lnr,lni,lnh,lnhl,lhw,   &
                           lscr,lsci,lsch,lschl,                &
                           lvh,lvhl,                            &
                           alphar,alphas,alphah,alphahl,        &
                           ntakid,ntakpd,ntakrd,ntakit,ntakid,rnu,snu
   USE MICRO_MODULE, only: rho_qh,rho_qhl,ipelec,takrhoi
   use takcommon,       only: lmax,lfmax,kfmax,iimax,kimax,lr100,ls250,hjo=>dj,hjos=>dji
   use comm1,           only: rhof
   use comm2,           only: sr,sx,sxi,srf,sri,fw,vw,bew,cdw,szr,szf,szi      &
     &                            ,sri,sfi,di,sxf,sxi,sh,ff,vf,fi,vi       

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   integer :: time
   real    :: time_real
   real    :: den(-ng+1:nz+ng,2)
   integer :: onedoutput
   
   integer :: n,l,il,nh
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: temlo,temhi
   real    :: temp, tmpq, tmpn, tmpz
   real    :: fac, diff,xmas,xmas0
   integer :: lx, lnx, lzx, lxw, lvx, lscx
   real,save    :: xden, alp, g1
   real    :: xpi = 3.14159
   real    :: xv,dtmp,rdia,totn,totq,totq2
   real    :: gamma
   real    :: xvr0, xvr1, xvri, vrnu, rwdia0, gfrnu, d2v, xmass0, lam
   real,save    :: aa,bb,mu,nu,tmp
   
   double precision :: masstot,ztmpi
   
   real, save :: rdrda(ntakpd),rda(ntakpd),rma(ntakpd)=0.0,rva(ntakpd),rvna(ntakpd),rdmda(ntakpd)
   real, save :: rmna(ntakpd)

      real    :: xrtmp(ntakpd),xhtmp(ntakpd),xhltmp(ntakpd)
      real    :: zrtmp(ntakpd),zhtmp(ntakpd),zhltmp(ntakpd)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   real :: temg1d(-ng+1:nz+ng)
!-----------------------------------------------------------------------------


!      write(0,*) 'In BINFORCE'

    IF ( sr(1) == 0.d0 ) THEN
      write(0,*) 'BINFORCE call takinitbin 1'
      call takinitbin
      write(0,*) 'BINFORCE call takinitbin, sx(1) = ',sx(1),sxf(1,1),sxf(1,2),lmax,lfmax

!        il = 2 ! 1 for graupel, 2 for frozen drops
!        DO l = 1,ntakpd
!         rma(l)  = 1.e-3*sxf(l,il) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
!         rva(l)  = rma(l)/(1.e3*rhof(il))  ! (m^3) volume (mass/1000.)
!         rda(l)  = 1.e-2*sx(1) ! (cm to meters) (6.*rma(l)/(pi*xdn))**(1./3.)
!!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
!         rvna(l) = 1.
!         rdrda(l) = rda(l)/hjo
!        ENDDO
    ENDIF

   work_to_do = .false.
   idofor(:)  = 0
   
   temlo = 253.15
   temhi = 273.15
   
   n = 1
   nssfor = 1
! only one region for now, and do it for all time.
   
   DO n = 1,nssfor

!
    IF ( tsstrt(n) == 0 .and. tsstop(n) == 0 ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ELSEIF ( time_real .ge. tsstrt(n) .and. time_real .le. tsstop(n) ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ENDIF

   ENDDO


!      write(0,*) 'BINFORCE work_to_do = ',work_to_do
      
    IF ( work_to_do ) THEN
    
!    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
!       write(0,*) 'BINFORCE: 1'
       
       kzb = nz-2
       kze = nz-2
!       if (kzbeg .eq. nzbeg) kzb = 1
!       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

   lxw = 0
   lvx = 0
   
   IF ( itopforcetype == 1 .or. itopforcetype == 10) THEN ! rain
     nh = ntakrd
     lx = lr
     lnx = lnr
     lscx = lscr
     lzx = 0 ! lzr
     xden = 1000.
     alp = alphar
     
          IF ( itopforcerainopt == 2 ) THEN
             igamrain = 2
             alp = 0
          ENDIF
          
         IF ( igamrain .eq. 1 ) THEN ! gamma of volume
         alp = rnu
         ENDIF
     IF ( itopforcetype == 10 .and. time_real <= dt ) THEN
       kzb = 1
     ENDIF
          nu = rnu
          mu = 1.
          
          bb = ( gamma( (nu + 1.)/mu ) / gamma( (nu+2.)/mu ) )**(-mu)
          
          aa = mu*bb**((nu+1.)/mu)/gamma( (nu+1.)/mu )
          
          
        IF ( rma(1) == 0.0 ) THEN
       write(0,*) 'binforce: set up rain'
        DO l = 1,nh
         rma(l)  = 1.e-3*sx(l) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
         rva(l)  = rma(l)/(xden)  ! (m^3) volume (mass/1000.)
!         rda(l)  = 1.e-2*sr(l) ! (cm to meters) (6.*rma(l)/(pi*xdn))**(1./3.)
         rda(l)  = 2.*1.e-2*sr(l) ! (cm to meters and radius to diameter) (6.*rma(l)/(pi*xdn))**(1./3.)
         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
         
         IF ( igamrain .eq. 1 ) THEN ! gamma of volume
         alp = rnu
         rvna(l) = (alp + 1.)**(alp + 1.) * rva(l)**alp
!         rdrda(l) = rda(l)/hjo
         rdrda(l) = 3.*rva(l) /hjo


!          rvna(l) = 1.
!          rdrda(l) = rda(l)/hjo
          
          write(0,*) 'l,rma,rma,rva,rvna,rda,rdrda = ',l,rma(l),rma(1)*exp(3.*(l-1)/hjo),rva(l),rvna(l),rda(l),rdrda(l)

         write(0,*) 'l,rma,rva,rvna,rda,rdrda = ',l,rma(l),rva(l),rvna(l),rda(l),rdrda(l)
         ELSE      ! gamma of diameter   
         rvna(l) = 1.
         rdrda(l) = rda(l)/hjo
         write(0,*) 'l,rma,rda,rdrda = ',l,rma(l),rda(l),rdrda(l)
         ENDIF
        
        ENDDO
        ENDIF
        
   ELSEIF ( itopforcetype == 2 .or. itopforcetype == 20 ) THEN ! graupel
     nh = ntakpd
     lx = lh
     lnx = lnh
     lscx = lsch
     lzx = 0 ! lzh
     lvx = lvh
     xden =  1.e3*rhof(1) ! rho_qh
     alp = alphah
     IF ( itopforcetype == 20 .and. time_real <= dt ) THEN
       kzb = 1
     ENDIF
        il = 1 ! 1 for graupel, 2 for frozen drops
        IF ( rma(1) == 0.0 ) THEN
       write(0,*) 'binforce: set up graupel'
         DO l = 1,ntakpd
         rma(l)  = 1.e-3*sxf(l,il) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
         rva(l)  = rma(l)/(xden)  ! (m^3) volume (mass/1000.)
         rda(l)  = 2.*1.e-2*sr(l) ! (cm to meters and radius to diameter) (6.*rma(l)/(pi*xdn))**(1./3.)
         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
         rvna(l) = 1.
         rdrda(l) = rda(l)/hjo
         write(0,*) 'l,rma,rda,sr2,rdrda = ',l,rma(l),rda(l),2.*1.e-2*sr(l),rdrda(l)
        ENDDO
        ENDIF
   ELSEIF ( .false. .and. itopforcetype == 3 ) THEN ! mixed-phase graupel
     nh = ntakpd
     lx = lh
     lxw = lhw
     lnx = lnh
     lscx = lsch
     lzx = 0 !lzh
     xden = 1000.
     lvx = 0 ! lvh
     alp = alphah
   ELSEIF ( itopforcetype == 4 ) THEN ! hail
     nh = ntakpd
     lx = lhl
     lnx = lnhl
     lscx = lschl
     lzx = 0 ! lzhl
     lvx = lvhl
     lxw = 0 ! lhlw
     xden =  1.e3*rhof(2) ! rho_qhl
     alp = alphahl
        il = 2 ! 1 for graupel, 2 for frozen drops
        IF ( rma(1) == 0.0 ) THEN
        write(0,*) 'binforce: set up hail'
        DO l = 1,ntakpd
         rma(l)  = 1.e-3*sxf(l,il) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
         rva(l)  = rma(l)/(xden)  ! (m^3) volume (mass/1000.)
         rda(l)  = 2.*1.e-2*sr(l) ! (cm to meters and radius to diameter) (6.*rma(l)/(pi*xdn))**(1./3.)
         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
         rvna(l) = 1.
         rdrda(l) = rda(l)/hjo
          write(0,*) 'l,rma,rda,sr2,rdrda = ',l,rma(l),rda(l),2.*1.e-2*sr(l),rdrda(l)
        ENDDO
        ENDIF
   ELSEIF ( itopforcetype == 5 .or. itopforcetype/10 == 5) THEN ! ice crystals
!     write(0,*) 'Binforce, ice crystals'
     nh = ntakid
     lx = li
     lnx = lni
     lscx = lsci
     lzx = 0 ! lzh
     lvx = 0 ! lvh
!     xden =  1.e3*0.9 !
     xden =  1.e3*takrhoi !

          nu = snu
          mu = 1.
          
          bb = ( gamma( (nu + 1.)/mu ) / gamma( (nu+2.)/mu ) )**(-mu)
          
          aa = mu*bb**((nu+1.)/mu)/gamma( (nu+1.)/mu )
          
        il = 1 ! use il for thickness parameter
        IF ( itopforcetype >= 50 ) il = Min(ntakit, Max(1, itopforcetype - 50) )
        lnx = lni + (il-1)*ntakid
        lscx = lsci  + (il-1)*ntakid
         d2v = 2.
        
        IF ( rma(1) == 0.0 ) THEN
        write(0,*) 'binforce: set up ice crystals, lnx,lscx,lsci = ',lnx,lscx,lsci
        DO l = 1,ntakid
         rma(l)  = 1.e-3*sxi(l,il) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
         rva(l)  = rma(l)/(xden)  ! (m^3) volume (mass/1000.)
        ! need to do distribution on mass instead of volume
        ! rva(l)  = pii/6.*(2.*1.e-2*sri(l))**2 ! rma(l)/(xden)  ! (m^3) volume (mass/1000.)
         rda(l)  = 2.*1.e-2*sri(l) ! (cm to meters and radius to diameter) (6.*rma(l)/(pi*xdn))**(1./3.)
!         rva(l)  = pii/6. * rda(l)**3
!         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
!         d2v = 2.474
         IF ( igamsnow == 1 ) THEN ! gamma of volume
           alp = snu
           rvna(l) = (alp + 1.)**(alp + 1.) * rva(l)**alp
           rmna(l) = (alp + 1.)**(alp + 1.) * rma(l)**alp
!           rdrda(l) = 3.*rva(l) /hjos
           rdrda(l) = d2v*rva(l) /hjos
           rdmda(l) = 2.*rma(l)/hjos
!          rvna(l) = 1.
!          rdrda(l) = rda(l)/hjos
          
          IF ( topforceN > 0.0 ) THEN
            write(0,*) 'l,rma,rma,rda,rdmda,rmna = ',l,rma(l),rma(1)*exp(d2v*(l-1)/hjos),rda(l),rdmda(l),rmna(l)
          ELSE
            write(0,*) 'l,rma,rma,rda,rdrda = ',l,rma(l),rma(1)*exp(d2v*(l-1)/hjos),rda(l),rdrda(l)
          ENDIF
         ELSE
           alp = alphas
           rvna(l) = 1.
           rdrda(l) = rda(l)/hjos
          write(0,*) 'l,rma,rda,sr2,rdrda = ',l,rma(l),rda(l),2.*1.e-2*sri(l),rdrda(l)
         ENDIF
        ENDDO
        ENDIF
   
   ELSE
     write(0,*) 'Bad value for itopforcetype = ',itopforcetype
     call commasmpi_abort()
   ENDIF

!   write(0,*) 'Binforce, nssfor loop'
   DO n = 1,nssfor
!       write(0,*) 'bin forcing: n,iqshap,idfor = ',n,iqshap(n),idofor(n)
     IF ( abs(iqshap(n)) >= 1 ) THEN
      
      IF ( idofor(n) .eq. 1 ) THEN

!       write(0,*) 'bin forcing at time = ',time_real
       xmass0 = 0.0
       do k = kzb,kze
        z = gzt(k,1)

          IF ( topforceN > 0 ) THEN
            xmass0 = den(k,1)*topforceq/topforceN
          ENDIF
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)

         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-topforceX)/topforceR)**2  &
                      +((y-topforceY)/topforceR)**2)
         ELSE
         radius = sqrt(((x-topforceX)/topforceR)**2 )
         ENDIF

!       write(0,*) 'i,j,k,radius = ',i,j,k,radius

       IF ( radius .lt. 1.0 ) THEN
        IF ( topforceq .gt. 0.0 ) THEN
! supersat calculation...
        IF ( iqshap(n) .eq. 1 .or. iqshap(n) == 2 ) THEN
         fac = (cos(xpi*radius)+1.)/2.
        ELSE
         fac = 1.0
        ENDIF

         tmpq = topforceq*fac

         
         
         st(i,j,k,lx) =  tmpq
         
!         IF ( lxw > 0 ) st(i,j,k,lxw) = st(i,j,k,lx)
         
!         IF ( lvx > 0 ) THEN
!           st(i,j,k,lvx) = tmpq/xden ! leave as mixing ratio here.
!         ENDIF

!            write(0,*) 'i,k,tmpq = ',i,k,tmpq
         
         tmpn = 0.0
!         IF ( lnx > 0 ) THEN
           IF ( tmpq > 0. ) THEN


   !         IF ( topforceN > 0 ) THEN
   !           tmpn = topforceN*(fac + 1.)*0.5
   !           st(i,j,k,lnx) =  tmpn/den(k,1)
   !           
   !         ELSE
              IF ( iqshap(n) == 2 ) THEN
                temp = topforcedia
              ELSE
                temp = topforcedia*(fac + 1.)*0.5
              ENDIF
             dtmp = temp
             rdia = temp*((3. + alp)*(2. + alp)*(1.0 + alp))**(-1./3.)
   !         ENDIF
           IF ( itopforcetype == 5 ) THEN
            xv =  1.e-6*(xpi*(100.*dtmp/2.)**2)*0.0141*(100.*dtmp)**0.474
           ELSE
             xv = xpi/6.*temp**3
           ENDIF
           
            IF ( topforceN > 0 ) THEN
              tmpn = topforceN*(fac + 1.)*0.5
!              st(i,j,k,lnx) =  tmpn/den(k,1)
              
            ELSE
             tmpn = tmpq*den(k,1)/(xden*xv)
            ENDIF

!           st(i,j,k,lnx) =  tmpn/den(k,1)
!           ELSE
!           st(i,j,k,lnx) = 0.0
        
        IF ( igamrain == 1 .and. ( itopforcetype == 1 .or. itopforcetype == 10 ) ) THEN
          rwdia0 = temp
          xvr0 = pii/6.0*(rwdia0)**3  
          rdia = rwdia0/6**(1./3.)

          xvr1 = xvr0
          xv   = xvr0
          xvri = 1./xvr1
          gfrnu = gamma(1. + alp)
          vrnu = (xvri)**alp/(xvr1 * gfrnu )
          tmpn = tmpq*den(k,1)/(xden*xvr0)
          
          
        ENDIF
          IF ( itopforcerainopt == 2 .and. ( itopforcetype == 1 .or. itopforcetype == 10 ) ) THEN
            lam = 4.1e3*topforcerainrate**(-0.21)
            tmpn =   topforcen0/lam
          ENDIF

        IF ( igamsnow == 1 .and. itopforcetype == 5) THEN
          rwdia0 = temp
          xvr0 = pii/6.0*(rwdia0)**3  
!            xvr0 =  1.e-6*(xpi*(100.*dtmp/2.)**2)*0.0141*(100.*dtmp)**0.474
          rdia = rwdia0/6**(1./3.)

          xvr1 = xvr0
          xv   = xvr0
          IF ( xmass0 > 0.0 ) THEN
            xvri = 1./xmass0
            gfrnu = gamma(1. + alp)
            vrnu = (xvri)**(alp+1.)/( gfrnu )
           ! tmpn = tmpq*den(k,1)/(xden*xvr0)
          ELSE
            xvri = 1./xvr1
            gfrnu = gamma(1. + alp)
            vrnu = (xvri)**alp/(xvr1 * gfrnu )
            tmpn = tmpq*den(k,1)/(xden*xvr0)
          ENDIF
        ENDIF

          IF ( i == nx/2 .and. time_real < 21.  ) THEN
            write(0,*) 'binforce: time_real = ',time_real
            write(0,*) 'binforce: i,k = ',i,k
            write(0,*) 'rdia,xv,tmpn,tmpq,den = ',rdia,xv,tmpn,tmpq,den(k,1)
            write(0,*) 'xmass0 = ',xmass0*1.e6
            write(0,*) 'ls250,xden,vrnu,alp,xvri = ',ls250,xden,vrnu,alp,xvri
            write(0,*) 'dtmp,gfrnu,rwdia0,A,B = ',dtmp,gfrnu,rwdia0,aa,bb
          ENDIF

         totn = 0.0
         totq = 0.0
         totq2 = 0.0
         tmp = 0.0
         IF ( my_rank == 0 .and. time_real == dt .and. i == onedoutput  ) THEN
           IF ( itopforcetype == 1 .or. itopforcetype == 10 ) THEN
           write(0,*) 'init spectra at i,k = ',onedoutput,k
!           write(0,*) 'il, sr, nr, nh, nhl, Nr(D), Nh(D), Nhl(D)'
           write(0,*) 'il, sr, nr, Nr(D)'
           ELSEIF ( itopforcetype == 5 .or. itopforcetype == 50 ) THEN
           write(0,*) 'init spectra at i,k = ',onedoutput,k
           write(0,*) 'il, sri, ni, Ni(D)'
           ENDIF
         ENDIF
         masstot = 0.d0
         ztmpi = 0.d0
         DO l = 1,nh
         ! NOTE: this is called from a point where num. conc. has been converted to
         !       number mixing ratio, so have to divide by air density.
         IF ( igamrain == 1 .and. ( itopforcetype == 1 .or. itopforcetype == 10 ) ) THEN ! gamma volume for rain
         
          IF ( xmass0 > 0.0 ) THEN
            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*rmna(l)*vrnu*Exp(  -(alp+1.)*rma(l)/xmass0 )*rdmda(l)/den(k,1)
          ELSE
!            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*xvri**alp*xv*rvna(l)*vrnu*Exp(  -(alp+1)*rva(l)*xvri )*rdrda(l)
            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*rvna(l)*vrnu*Exp(  -(alp+1.)*rva(l)*xvri )*rdrda(l)/den(k,1)
          ENDIF

            masstot = masstot + st(i,j,k,lnx+l-1)*sx(l)

!           st(i,j,k,lnx+l-1) = 1.e-6*3.*tmpn*aa*rda(l)**(3.*nu+2.)/(rwdia0**(3.*nu+3.) )*Exp( -bb*(rda(l)/rwdia0)**(3.*mu) )*rdrda(l)
!           tmp = 1.e-6*3.*tmpn*aa*rda(l)**(3.*nu+2.)/(rwdia0**(3.*nu+3.) )*Exp( -bb*(rda(l)/rwdia0)**(3.*mu) )

          xrtmp(l) = st(i,j,k,lnx+l-1)
          

         ELSEIF ( igamsnow == 1 .and. ( itopforcetype == 5 .or. itopforcetype == 50 ) ) THEN ! gamma volume for snow
         
!            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*xvri**alp*xv*rvna(l)*vrnu*Exp(  -(alp+1)*rva(l)*xvri )*rdrda(l)
          IF ( xmass0 > 0.0 ) THEN
            ! rmna(l) = (alp + 1.)**(alp + 1.) * rma(l)**alp
            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*rmna(l)*vrnu*Exp(  -(alp+1.)*rma(l)/xmass0 )*rdmda(l)/den(k,1)
            masstot = masstot + st(i,j,k,lnx+l-1)*sxi(l,1)
            ztmpi = ztmpi+st(i,j,k,lnx+l-1)*(((6.0/pii))*sxi(l,1))**2
            
            IF ( i == onedoutput .and. k == kze ) THEN
              write(0,*) 'l,rma,rdmda,st = ',l,rma(l),rdmda(l),st(i,j,k,lnx+l-1)*den(k,1)
            ENDIF
          ELSE
            st(i,j,k,lnx+l-1) = 1.e-6*tmpn*rvna(l)*vrnu*Exp(  -(alp+1.)*rva(l)*xvri )*rdrda(l)/den(k,1)
            masstot = masstot + st(i,j,k,lnx+l-1)*sxi(l,1)
          ENDIF

            IF ( lscr > 1 .and. ipelec > 0 ) THEN
              st(i,j,k,lscx+l-1) = 1.e6*st(i,j,k,lnx+l-1)*chargeperparticle
            ENDIF

!           st(i,j,k,lnx+l-1) = 1.e-6*3.*tmpn*aa*rda(l)**(3.*nu+2.)/(rwdia0**(3.*nu+3.) )*Exp( -bb*(rda(l)/rwdia0)**(3.*mu) )*rdrda(l)
!           tmp = 1.e-6*3.*tmpn*aa*rda(l)**(3.*nu+2.)/(rwdia0**(3.*nu+3.) )*Exp( -bb*(rda(l)/rwdia0)**(3.*mu) )
          xrtmp(l) = st(i,j,k,lnx+l-1)
         ELSE
           IF ( alp < 14. ) THEN ! {
            st(i,j,k,lnx+l-1) = 1.e-6*tmpn/(den(k,1)*gamma(1.+alp)*rdia**(1.+alp))*rda(l)**alp*Exp( -rda(l)/rdia )*rdrda(l)
            IF ( lscr > 1 .and. ipelec > 0 ) THEN
              st(i,j,k,lscx+l-1) = 1.e6*st(i,j,k,lnx+l-1)*chargeperparticle
            ENDIF
           IF ( lvx > 0 ) THEN
             st(i,j,k,lvx+l-1) = st(i,j,k,lnx+l-1)*rva(l) ! rma(l)/xden
           ENDIF
           
           ELSE ! {
            st(i,j,k,lnx+l-1) = 0.0
            IF ( l > 1 .and. l < nh ) THEN
              IF ( 0.5*(rda(l-1)+rda(l)) < dtmp .and.  dtmp < 0.5*(rda(l)+rda(l+1)) ) THEN
                 st(i,j,k,lnx+l-1) = 1.e-6*tmpn/(den(k,1))
                IF ( lscr > 1 .and. ipelec > 0 ) THEN
                  st(i,j,k,lscx+l-1) = 1.e6*st(i,j,k,lnx+l-1)*chargeperparticle
                ENDIF
              ENDIF
            ENDIF
           ENDIF ! } }
          ENDIF
          totn = totn + 1.e6*st(i,j,k,lnx+l-1) ! *den(k,1)
          totq = totq + 1.e6*st(i,j,k,lnx+l-1)*rma(l) ! *den(k,1)

            masstot = masstot + st(i,j,k,lnx+l-1)*rma(l)

!          IF ( l > ls250 ) totq2 = totq2 + 1.e6*st(i,j,k,lnx+l-1)*rma(l)
!          totq2 = totq2 + 1.e3*st(i,j,k,lnx+l-1)*sxf(l,2)*den(k,1)
          
          xrtmp(l) = st(i,j,k,lnx+l-1)
          
          IF ( ny == 2 .and. i == nx/2 .and. time_real < 21. ) THEN
            write(0,*) 'l,lnx+l-1,rda,rdrda,rvnast,totn,totq,tmp = ',l,lnx+l-1,rda(l),rdrda(l),rvna(l),  &
     &          st(i,j,k,lnx+l-1),totn,1000.*totq/den(k,1),tmp ! ,1000.*totq2/den(k,1)
          ENDIF
         ENDDO
         
         !  write(0,*) 'masstot,q(g/kg),ztmpi = ',1.e6*masstot*den(k,1),1.e6*masstot, 10.0*Log10(Max(1.d-40,1.d12*0.224*ztmpi)) 

           ENDIF
         
         IF ( time_real == dt .and. i == onedoutput ) THEN
           IF ( my_rank == 0 .and. ( itopforcetype == 1 .or. itopforcetype == 10 )  ) THEN
             DO il = 1,lfmax
!         write(6,'(i2,7(2x,1pe12.5))') il, sr(il),xrtmp(il),xhtmp(il),xhltmp(il),   &
!     &          xrtmp(il)/(2.*sr(il)/hjo),xhtmp(il)/(2.*sr(il)/hjo),xhltmp(il)/(2.*sr(il)/hjo)
            write(0,'(i2,7(2x,1pe12.5))') il, sr(il),xrtmp(il), xrtmp(il)/(2.*sr(il)/hjo)
             ENDDO
           ELSEIF ( my_rank == 0 .and. ( itopforcetype == 5 .or. itopforcetype == 50 )  ) THEN
             DO il = 1,iimax
!         write(6,'(i2,7(2x,1pe12.5))') il, sr(il),xrtmp(il),xhtmp(il),xhltmp(il),   &
!     &          xrtmp(il)/(2.*sr(il)/hjo),xhtmp(il)/(2.*sr(il)/hjo),xhltmp(il)/(2.*sr(il)/hjo)
             write(0,'(i2,7(2x,1pe12.5))') il, sri(il),xrtmp(il), xrtmp(il)/(2.*sri(il)/hjos)
             ENDDO
            ENDIF
          ENDIF
         
         IF ( tmpq > 0.0 ) THEN
!         write(0,*) 'i,xden = ',i,lvx,xden,st(i,j,k,lvx), den(k,1)*st(i,j,k,lx)/(st(i,j,k,lvx)*den(k,1))
!         xden = den(k,1)*st(i,j,k,lx)/(st(i,j,k,lvx)*den(k,1))
!         xv = den(k,1)*st(i,j,k,lx)/(xden*st(i,j,k,lnx)*den(k,1))
         
!         write(0,*) 'xv,dia = ',xv,dtmp,(6.0*xv/xpi)**(1./3.)
         ENDIF
         
         ENDIF ! topforceq
       ENDIF ! radius

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF
      
     ENDIF ! iwshap
   
   ENDDO

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
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
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
!        fw(i,j,k) = fw(i,j,k) + t0(i,j,k)
        
      ENDDO
     ENDDO
    ENDDO
    
    ENDIF ! work_to_do


   RETURN
   
   END SUBROUTINE BINFORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE CCNUFFORCE      ////////////////////
!     Forcing of ultrafine CCN
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE CCNUFFORCE(nx,ny,nz,ns,dt,st,pt,t0,piinit1d,sbase,den,    &
                     gxt,gyt,gzt,time,time_real)

   USE PARAM_MODULE, only: ng, pii
   USE COMMASMPI_MODULE
   USE INDEX_MODULE
   USE MICRO_MODULE, only: rho_qh,rho_qhl,imurain

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   integer :: time
   real    :: time_real
   real    :: den(-ng+1:nz+ng,2)

   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: temlo,temhi
   real    :: temp, tmpq, tmpn, tmpz
   real    :: fac, diff,xmas,xmas0
   integer :: lx, lnx, lzx, lxw, lvx
   real    :: xden, alp, g1
   real    :: xpi = 3.14159
   real    :: xv,dtmp
   integer :: ifirst = 0

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   real :: temg1d(-ng+1:nz+ng)
!-----------------------------------------------------------------------------


!      write(0,*) 'In Qforce'

   work_to_do = .false.
   idofor(:)  = 0

   temlo = 253.15
   temhi = 273.15

   n = 1
   nssfor = 1
! only one region for now, and do it for all time.

   DO n = 1,nccnfor

    IF ( time_real .ge. tccnstrt(n) .and. time_real .le. tccnstop(n) ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ENDIF

   ENDDO

    IF ( work_to_do ) THEN

!    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
!       write(0,*) 'CCNFORCE: 1'


       kzb = 1
       kze = nz-2
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

   lxw = 0
   lvx = 0

   lx = lccnuf
   lnx = lccnuf
   
   IF ( lccnuf < 1 ) THEN
     write(0,*) 'cannot do CCNUF forcing because lccnuf = ',lccnuf
     return
   ENDIF


   DO n = 1,nccnfor
     IF ( abs(iqshap(n)) >= 1 ) THEN

      IF ( idofor(n) .eq. 1 ) THEN

!       write(*,*) 'SS forcing at time = ',time_real
       do k = kzb,kze
        z = gzt(k,1)

        do j = jyb,jye
        y = gyt(j,1)

         do i = ixb,ixe
         x = gxt(i,1)

         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xccnfcen(n))/xccnfrad(n))**2  &
                      +((y-yccnfcen(n))/yccnfrad(n))**2  &
                      +((z-zccnfcen(n))/zccnfrad(n))**2)
         ELSE
         radius = sqrt(((x-xccnfcen(n))/xccnfrad(n))**2  &
                      +((z-zccnfcen(n))/zccnfrad(n))**2)
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
        IF ( topforceN .gt. 0.0 ) THEN
! supersat calculation...
!         IF ( iqshap(n) .eq. 1 .or. iqshap(n) .eq. 2) THEN
!          fac = (cos(xpi*radius)+1.)/2.
!         ELSE
!          fac = 1.0
!         ENDIF
         fac = 1.0

         tmpq = topforceN*fac

         st(i,j,k,lccnuf) =  tmpq

         
         ENDIF ! topforceq
       ENDIF ! radius

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF
      
     ENDIF ! iwshap
   
   ENDDO

   
    ENDIF ! work_to_do


   RETURN
   
   END SUBROUTINE CCNUFFORCE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE CHARGEFORCE  ////////////////////
!     Introduces charge into cylindrical regions at a linear rate dQ/dt
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE CHARGEFORCE(nx,ny,nz,ns,dt,st,pt,t0,piinit1d,sbase,den,    &
                     gxt,gyt,gzt,time,time_real)

   USE PARAM_MODULE, only: ng, pii
   USE COMMASMPI_MODULE
   USE INDEX_MODULE
   USE MICRO_MODULE, only: rho_qh,rho_qhl,imurain

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   integer :: time
   real    :: time_real
   real    :: den(-ng+1:nz+ng,2)

   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius,radiusz
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: temlo,temhi
   real    :: temp, tmpq, tmpn, tmpz
   real    :: fac, diff,xmas,xmas0
   integer :: lx, lnx, lzx, lxw, lvx, lscx
   real    :: xden, alp, g1
   real    :: xpi = 3.14159
   real    :: xv,dtmp
   integer :: ifirst = 0 ! setting value here saves the value going forward

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   real :: temg1d(-ng+1:nz+ng)
!-----------------------------------------------------------------------------


!      write(0,*) 'In Qforce'

   work_to_do = .false.
   idofor(:)  = 0


   DO n = 1,nchgfor

        work_to_do = .true.
      idofor(n) = 1

   ENDDO

    IF ( work_to_do ) THEN

       kzb = 1
       kze = ktile
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

     lx = li
     lnx = lni
     lscx = lsci
     IF ( lscx < 1 ) THEN
       write(0,*) 'No charge variable! Must turn on elec!'
       STOP
     ENDIF

   IF ( ifirst == 0 .and. xchgcldfrad(1) > 0. ) THEN
   ! set up cloud
       n = 1
       do k = kzb,kze
        z = gzt(k,1)

        do j = jyb,jye
        y = gyt(j,1)

         do i = ixb,ixe
         x = gxt(i,1)

         ! check if point is within a cylinder
         IF ( ny .gt. 2 ) THEN

         radius = sqrt(((x-xchgcldfcen(n))/xchgcldfrad(n))**2  &
                      +((y-ychgcldfcen(n))/xchgcldfrad(n))**2)
         ELSE
         radius = sqrt(((x-xchgcldfcen(n))/xchgcldfrad(n))**2 )
         ENDIF
         
         radiusz = Abs(z-zchgcldfcen(n))/zchgcldfrad(n) ! cylinder depth


        IF ( radius .lt. 1.0 .and. radiusz <= 1.0 ) THEN

          IF ( ifirst == 0 ) st(i,j,k,lc) =  0.5e-3 ! set nominal cloud droplet mixing ratio

        ENDIF ! radius

         ENDDO
        ENDDO
       ENDDO
   
   ENDIF ! cloud setup

   DO n = 1,nchgfor
     IF ( abs(ichgshap(n)) >= 1 ) THEN

      IF ( idofor(n) .eq. 1 ) THEN

!       write(*,*) 'CHG forcing at time = ',time_real
       do k = kzb,kze
        z = gzt(k,1)

        do j = jyb,jye
        y = gyt(j,1)

         do i = ixb,ixe
         x = gxt(i,1)

         ! check if point is within a cylinder
         IF ( ny .gt. 2 ) THEN

         radius = sqrt(((x-xchgfcen(n))/xchgfrad(n))**2  &
                      +((y-ychgfcen(n))/xchgfrad(n))**2)
         ELSE
         radius = sqrt(((x-xchgfcen(n))/xchgfrad(n))**2 )
         ENDIF
         
         radiusz = Abs(z-zchgfcen(n))/zchgfrad(n) ! cylinder depth


        IF ( radius .lt. 1.0 .and. radiusz <= 1.0 ) THEN

          IF ( iqshap(n) .eq. 1 ) THEN
            fac = 1.0
          ELSE ! (iqshap(n) .eq. 2) THEN
            fac = (cos(xpi*radius)+1.)/2.
          ENDIF

          tmpq = chgratefor(n)*fac/den(k,1) ! set as charge 'mixing ratio' (converted to C/m**3 in microphysics)

          st(i,j,k,lscx) = st(i,j,k,lscx) + tmpq*dt

          IF ( ifirst == 0 ) st(i,j,k,lx) =  1.0e-3 ! set nominal cloud ice mixing ratio as placeholder

        ENDIF ! radius

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF ! idofor
      
     ENDIF ! ichgshap
   
   ENDDO ! n


      ifirst = 1 ! set flag that hydrometeor mixing ratio has been set

    ENDIF ! work_to_do

   RETURN
   
   END SUBROUTINE CHARGEFORCE

!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE BFOR      ////////////////////
!     Forces boundaries
!-----------------------------------------------------------------------------
   SUBROUTINE BFOR()
   
   implicit none
   
   
   RETURN
   
   END SUBROUTINE BFOR

!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE BLCOOL      ////////////////////
!     Forces boundaries
!-----------------------------------------------------------------------------
!   SUBROUTINE BLCOOL(nx,ny,nz,dt,time_real,ns,den,pi,piinit,ut,vt,gxt,gyt,gzt,st,sbase) 
   SUBROUTINE BLCOOL(nx,ny,nz,ns,dt,w,st,pt,t0,piinit1d,sbase,        &
                     gxt,gyt,gzt,time_real,loop,den,   &
                     uinit,vinit,ugrid,vgrid,ut,vt)

   USE PARAM_MODULE, only: ng, luno
   USE COMMASMPI_MODULE
!   USE INDEX_MODULE
    USE INDEX_MODULE, only : lt, lv

!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,ns
   real    :: dt
   real    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   real    :: pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pert. exner pi
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: ut(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: vt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
      real :: pi(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real :: piinit(-ng+1:nz+ng)
   real    :: piinit1d(-ng+1:nz+ng)
   real    :: sbase(-ng+1:nz+ng,ns)
   real    :: time_real
   integer :: loop
   real    :: uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng)
   real    :: ugrid, vgrid
   real    :: den(-ng+1:nz+ng,2)

! local vars

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: temp, qvs, tref
   real    :: fac, diff,xmas,xmas0
   
   
   IF ( low_level_cooling_flag < 1 ) RETURN

! Adjust boundary layer temperature as in Parker (2008)

      IF ( time_real > low_level_cooling_start .and. time_real <= low_level_cooling_end ) THEN


!
!  compute forcing based on shape I (cosine squared function)
!
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg-1

!       jyb = -ng+1
!       jye = jtile+ng
       jyb = 1
       jye = jtile
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

!       ixb = -ng+1
!       ixe = itile+ng
       ixb = 1
       ixe = itile
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

         tref = piinit1d(1)*sbase(1,lt) - (time_real - low_level_cooling_start)*low_level_cooling_rate/3600.
         
         write(luno,*) 'BLCOOL: time, tref ', time_real,tref


       DO k = kzb,kze
        z = gzt(k,1)
        
         IF ( z < low_level_cooling_depth ) THEN
         
!         tref = piinit1d(1)*sbase(1,lt) - (time_real - low_level_cooling_start)*low_level_cooling_rate
         
         DO j = jyb,jye
         y = gyt(j,1)
         
           DO i = ixb,ixe
           x = gxt(i,1)
         
           temp = st(i,j,k,lt)*( piinit1d(k) + pt(i,j,k) )
           IF ( temp > tref ) THEN
             st(i,j,k,lt) = tref/( piinit1d(k) + pt(i,j,k) )
           ENDIF
           
           ENDDO ! i
         
         ENDDO ! j
         
         ENDIF
         
        ENDDO ! k


      ENDIF
   
   
   RETURN
   
   END SUBROUTINE BLCOOL


!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! STRAKAS ATMOSPHERIC MODEL  (SAM)
!   Designed by Jerry M. Straka
!   Version 10.1 (01 JANUARY 2002) 
!
! SUBROUTINE SFCFLX
!    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      subroutine sfcflx(ugrid,vgrid,nx,ny,nz,dx,dy,dt,time,ns,den,pi,piinit               &
       ,ut,vt,gxt,gyt,gzt,st)                 
       
      USE PARAM_MODULE, only : g,rd,rcp,ng
      USE COMMASMPI_MODULE
      USE INDEX_MODULE, only : lt, lv
      USE INIT_MODULE, only: prg0,thg0 ,uspeed,timlndst,thg1,isstbeg,isstend,isfcl

!
! declarations for surface layer formulation
!
       implicit none

!
      integer        id1, jd1, istag, jstag, ig, ithg, iprg, itemg,itema,ivmean
      integer        nx,     ny,     nz,   ns, iuwind, ivwind,time
      integer        ix,     jy,     kz, nsphys, nssoil, iqvg, itflux, iqflux
      real            dzsl,    ugrid, vgrid
      real           dx, dy, dz, dt
      real           cdval
      real         Pres
      parameter (id1=1,jd1=1,istag=1,jstag=1,iuwind = 19, ivwind = 20, nsphys=29, nssoil=10,  &
                ithg = 1,iqvg =2, iprg = 3,itema=28,itemg = 29, itflux=12,iqflux=13,ivmean=18)
!      parameter (id1=1,jd1=1,istag=1,jstag=1,iuwind = 19, ivwind = 20, nsphys=29, nssoil=10, &
!      ithg = 1,iqvg =2, iprg = 3,itema=28,itemg = 29, itflux=12,iqflux=13,ivmean=18)
      real :: den(-ng+1:nz+ng,2)
!      real :: gtz(nz)
      real :: pi(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: piinit(-ng+1:nz+ng)
      real :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      real :: ut(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: vt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real, allocatable :: frsea(:,:)
      real, allocatable :: sphys(:,:,:)
      real, allocatable :: ssoil(:,:,:,:)

       real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)

      real :: xlandb = 0.0
      real :: ylandb = 0.0
      real :: xlandf = 0.0
      real :: ylandf = 0.0
      real :: fland = 0.0
      
      real :: fac, thg

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: im1, jm1

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag


      if (.not. allocated(tt13sl)) THEN
      allocate(tt13sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(tt23sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(td13sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(td23sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(twt3sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(thf3sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(tmf3sl(-ng+1:nx+ng,-ng+1:ny+ng))
      allocate(tdh3sl(-ng+1:nx+ng,-ng+1:ny+ng))
      ENDIF
      
      allocate ( frsea(-ng+1:nx+ng,-ng+1:ny+ng) )
      allocate(sphys(-ng+1:nx+ng,-ng+1:ny+ng,nsphys))
      allocate(ssoil(-ng+1:nx+ng,-ng+1:ny+ng,1,nssoil))
      sphys(:,:,:) = 0.0
!
!
!..constants in surface-layer and surface energy budget formulations
!
!


!      dzsl = 0.5*dz
!
! zero tt13sl, tt23sl, thf3sl, tmf3sl
!
!     write(0,*) 'Start sfcflx ',my_rank
     
    kzb = 1
    kze = ktile
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    jyb = 1
    jye = jtile
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg

    ixb = 1
    ixe = itile
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    
    IF ( .true. ) THEN
    DO jy = jyb-1,jye+1
     DO ix = ixb-1,ixe+1
      tt13sl(ix,jy) = 0.0 !0.0
      tt23sl(ix,jy) = 0.0
      td13sl(ix,jy) = 0.0
      td23sl(ix,jy) = 0.0
      twt3sl(ix,jy) = 0.0
      thf3sl(ix,jy) = 0.0
      tmf3sl(ix,jy) = 0.0
     ENDDO
    ENDDO
    
    ELSE

    DO jy = jyb-1,jye+1
     DO ix = ixb-1,ixe+1
      tt13sl(ix,jy) = -1.e32 !0.0
      tt23sl(ix,jy) = -1.e32
      td13sl(ix,jy) = -1.e32
      td23sl(ix,jy) = -1.e32
      twt3sl(ix,jy) = -1.e32
      thf3sl(ix,jy) = -1.e32
      tmf3sl(ix,jy) = -1.e32
     ENDDO
    ENDDO
    
    ENDIF


! vmean, u, v in the sfclayer ( first grid point in vertical ) 


! set kz = 1 for entire routine - DO NOT change

      kz = 1

    DO jy = jyb-1,jye+1
     DO ix = ixb-1,ixe+1

      sphys(ix,jy,iuwind) =            &
       (0.5)                           &
      *(ut(ix      ,jy,kz) + ugrid    &
      + ut(ix+istag,jy,kz) + ugrid)

      sphys(ix,jy,ivwind) =            &
       (0.5)                           &
      *(vt(ix,jy,kz      ) + vgrid    &
      + vt(ix,jy+jstag,kz) + vgrid)

     ENDDO
    ENDDO

!     write(0,*) 'sfcflx A ',my_rank


    DO jy = jyb,jye 
     DO ix = ixb,ixe
!c
!c  mean wind speed in sfc layer (first model layer)
!c
      sphys(ix,jy,ivmean) =                    &
       max(                                    &
          (sphys(ix,jy,iuwind)**2              &
          +sphys(ix,jy,ivwind)**2)**(0.5)      &
!          , 0.00)
          , 1.0e-02)

! convective wind adjustment
!
!      vcon =
!     >  (2.0)*(max((ssoil(ix,jy,01,ithg)-ac(ix,jy,01,01)),0.0))**(0.5)
!
!  update mean wind speed with vcon
!
!      sphys(ix,jy,ivmean) =
!     >  (sphys(ix,jy,ivmean)**2 + vcon**2)**(0.50)

     ENDDO
    ENDDO

! read in data at first time in sfcl


! Adjust sea surface temp (SST) if needed:
      thg = thg0
      IF ( isstend > isstbeg .and. thg0 /= thg1 ) THEN
          IF ( time >= isstend ) THEN
            thg = thg1
          ELSEIF ( time <= isstbeg ) THEN
            thg = thg0
          ELSE
             thg = thg0 + (time-isstbeg)*(thg1 - thg0)/(isstend - isstbeg)
          ENDIF
      ENDIF

     ig = 1
     DO jy = jyb,jye 
      DO ix = ixb,ixe
!       ssoil(ix,jy,ig,ithg)  = thg0 ! redefined later
!       ssoil(ix,jy,ig,iprg)  = prg0 ! redefined later
!       sphys(ix,jy,itemg) =  ssoil(ix,jy,01,ithg)*(((ssoil(ix,jy,01,iprg))/(1.0e05))**(rcp))
       ! this could be replaced by a scalar since it is constant 
       sphys(ix,jy,itemg) =  thg*(((prg0)/(1.0e05))**(rcp))
      ENDDO
     ENDDO

!     write(0,*) 'sfcflx B ',my_rank

! The following is for simple hurricane simulations
!   if you set 
!   icdsfs = 2

! straight coast - moving land with simulation
! forward speed of land (uspeed) in m/s
! timlndst is time to start land
! fland = land-sea barrier straight coast

!      uspeed = 8.0
!      timlndst = 37000.0


      IF ( time .gt. timlndst .and. uspeed > 0. ) THEN
        fland = Max(time-timlndst,0.0)*uspeed
        DO jy = jyb,jye 
          ylandb = 0.0
          ylandf = 0.0
          DO ix = ixb,ixe
            xlandb =  gxt(ix,1) ! (ix)*dx
            xlandf =  gxt(ix+1,1) ! (ix+1)*dx
            if ( fland .gt. xlandb ) then
              if ( fland .gt. xlandf ) then  ! grid box totally over land
                frsea(ix,jy) = 0.0
              end if
            end if
            if ( fland .lt. xlandb ) then
              if ( fland .lt. xlandf ) then  ! grid box totally over sea
                frsea(ix,jy) = 1.0
              end if
            end if
            if ( fland .gt. xlandb ) then
            if ( fland .lt. xlandf ) then  ! grid box over coastline
            frsea(ix,jy) = (fland-xlandb)/dx
            end if
            end if
          end do
        end do

      ELSE
        do jy = 1,ny
          do ix = 1,nx
            frsea(ix,jy) = 1.0
          enddo
        enddo
      ENDIF


!     write(0,*) 'sfcflx C ',my_rank

! end of frsea set up
!
!
! SEMI-SLIP ON HORIZONTAL VELOCITY (ISFCL = -1)
!
!  set kz=1 for the first theta point
!
! Surface dissipative heating parameterization following Bister and 
!   Emanuel 1998? (KE to IE conversion)
!
! start moisture and heat flux
!
      kz = 1
    DO jy = jyb,jye 
     DO ix = ixb,ixe
!
      if ( frsea(ix,jy) .lt. 0.50 ) then
        cdval = 0.0
      ELSE
!        cdval = (1.1e-03) + (4.0e-5)*sphys(ix,jy,ivmean)
        cdval = (0.5e-03) + (4.0e-5)*Min(40.0,sphys(ix,jy,ivmean)) ! ~Powell et al. (2003, Nature, Fig. 3)
      end if
!
      sphys(ix,jy,itema) = (st(ix,jy,01,lt))*(pi(ix,jy,01)+piinit(01))
      Pres=100000.*(piinit(01)+pi(ix,jy,01))**3.509 
      

!  
! NOTE: unless needed in another subroutine, can replace array (iprg) with a scalar
!       write(0,*) 'iprg parts = ', pres,g,1./gzt(kz,3),rd,sphys(ix,jy,itema)
      ssoil(ix,jy,01,iprg) =                        &
       (Pres)                                      &
!      * ( exp(g*(0.5)*(dz/gtz(kz))                  &
      * ( exp(g*(0.5)*(1./gzt(kz,3))                  &
       /(rd*(0.5)*(sphys(ix,jy,itema)))) )

! NOTE: unless needed in another subroutine, can replace array (ithg) with a scalar
      ssoil(ix,jy,01,ithg) =                        &
       sphys(ix,jy,itemg)                           &
       *(((1.0e+05)/(ssoil(ix,jy,01,iprg)))**(rcp))

! NOTE: unless needed in another subroutine, can replace array (iqvg) with a scalar
      ssoil(ix,jy,01,iqvg) =                                    &
       (380.0/(ssoil(ix,jy,01,iprg)))*exp(17.2693882*           &
       (sphys(ix,jy,itemg)-273.16)/(sphys(ix,jy,itemg)-35.86))

! NOTE: this is not used here before being deallocated
      sphys(ix,jy,itflux) = cdval*sphys(ix,jy,ivmean)     &
       * (ssoil(ix,jy,01,ithg)-st(ix,jy,kz,lt))           &
       * gzt(kz,3)
!       / (dz/gtz(kz))

! NOTE: this is not used here before being deallocated
      sphys(ix,jy,iqflux) = cdval*sphys(ix,jy,ivmean)     &
       * (ssoil(ix,jy,01,iqvg)-st(ix,jy,kz,lv))           &
       * gzt(kz,3)
!       / (dz/gtz(kz))

      fac = Min( cdval*sphys(ix,jy,ivmean), 0.75/dt) ! limit flux to less than the difference to target value
     ! Theta flux
      thf3sl(ix,jy) =                                     &
         -den(kz,1)*fac        &  ! cdval*sphys(ix,jy,ivmean)
       * (ssoil(ix,jy,01,ithg)-st(ix,jy,kz,lt))


     ! QV flux
      tmf3sl(ix,jy) =                                    &
         -den(kz,1)*fac         &  ! cdval*sphys(ix,jy,ivmean)
       * (ssoil(ix,jy,01,iqvg)-st(ix,jy,kz,lv))

!       write(91,*) 'sfc fluxes: i,j,TH,QV = ',ix,jy, thf3sl(ix,jy), tmf3sl(ix,jy)
!       write(91,*) 'thf3sl parts = ', -den(kz,1),cdval,sphys(ix,jy,ivmean),ssoil(ix,jy,01,iqvg),st(ix,jy,kz,lv)

      ! buoyancy flux term in TKE generation
      twt3sl(ix,jy) =                                    &
         fac       &  ! cdval*sphys(ix,jy,ivmean)                       &
       * (ssoil(ix,jy,01,ithg)-st(ix,jy,kz,lt))
     
!      IF ( (tmf3sl(ix,jy) .gt. 0.0) .or. (thf3sl(ix,jy) .gt. 0.0)) THEN 
!      write(luno,*) 'tmf3sl(ix,jy) ix  jy =' , tmf3sl(ix,jy),ix, jy
!      write(luno,*) 'thf3sl(ix,jy) ix  jy =' , thf3sl(ix,jy),ix, jy
!      ENDIF

!      write(luno,*) 'cdval, den (kz,1) =' ,cdval, den(kz,1),sphys(ix,jy,ivmean),ssoil(ix,jy,kz,ithg),st(ix,jy,kz,lv),st(ix,jy,kz,lt),ssoil(ix,jy,kz,ithg)
 
      end do
      end do

! end moisture and heat flux

!     write(0,*) 'sfcflx D (comms) ',my_rank

! MPI

#ifdef MPI
  ! need to communicate sphys

    IF ( number_of_processes .gt. 1 ) THEN
    
        
!        IF ( my_rank == 0 ) THEN
!          DO ix = 1,nsphys
!          write(0,*) 'i, sphys = ',ix,sphys(nx+1,ny/2,ix)
!          ENDDO
!        ENDIF
        
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,nsphys,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,sphys)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,nsphys,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,sphys)

        ENDIF

!        IF ( my_rank == 0 ) THEn
!          DO ix = 1,nsphys
!          write(0,*) 'i, sphys2 = ',ix,sphys(nx+1,ny/2,ix)
!          ENDDO
!        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,nsphys,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,sphys)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,nsphys,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,sphys)

        ENDIF


    
    ENDIF


#endif

!     write(0,*) 'sfcflx E ',my_rank

! u-velocity

!      do 30101 jy = 1,     ny-jstag
!      do 30102 ix = 1+id1, nx-id1
      DO jy = jyb,jye 
       DO ix = ixb,ixe
         im1 = ix - 1
         IF ( nxbeg == ixbeg ) im1 = Max(1,ix-1)

         cdval = (1.1e-03)+(4.0e-5)*( (0.50)*(sphys(ix,jy,ivmean)+sphys(im1,jy,ivmean)) )

         if ( frsea(ix,jy) .lt. 0.5 ) then
          cdval = 0.02
         end if

         td13sl(ix,jy) = 2.0*(ut(ix,jy,kz)+ugrid)*gzt(kz,3)

         tt13sl(ix,jy) =                                                &
           cdval*(ut(ix,jy,kz)+ugrid)                                  &
           *( (0.50)*(sphys(ix,jy,ivmean)+sphys(im1,jy,ivmean)) )    &
           *( (0.50)*(den(1,1)+den(1,1)) )
      
      ENDDO
      ENDDO
      
!30102 continue
!30101 continue

! v-velocity

!     write(0,*) 'sfcflx F ',my_rank

!      do 30111 jy = 1+jd1, ny-jd1
!      do 30112 ix = 1,     nx-istag
      DO jy = jyb,jye 
       DO ix = ixb,ixe
        jm1 = jy - 1
        IF ( nybeg == jybeg ) jm1 = Max(1,jy-1)

      cdval = (1.1e-03) + (4.0e-5)*( (0.50)*(sphys(ix,jy,ivmean)+sphys(ix,jm1,ivmean)) )

      if ( frsea(ix,jy) .lt. 0.5 ) then
      cdval = 0.02
      end if

      td23sl(ix,jy) =                                        &
        2.0*(vt(ix,jy,kz)+vgrid)*gzt(kz,3)

      tt23sl(ix,jy) =                                                 & 
        cdval*(vt(ix,jy,kz)+vgrid)                                  &
        *( (0.50)*(sphys(ix,jy,ivmean)+sphys(ix,jm1,ivmean)) )    &
        *( (0.50)*(den(1,1)+den(1,1)) )
    ENDDO
    ENDDO
! 30112 continue
! 30111 continue

!     write(0,*) 'sfcflx G ',my_rank

#ifdef MPI
  ! need to communicate td13sl,td23sl

    IF ( number_of_processes .gt. 1 ) THEN
    
        
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,td13sl)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,td13sl)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,td13sl)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,td13sl)

        ENDIF


        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,td23sl)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,td23sl)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,td23sl)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,td23sl)

        ENDIF

        
        ! communicate tmf3sl (QV flux)
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,tmf3sl)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,tmf3sl)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,tmf3sl)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,tmf3sl)

        ENDIF


        ! communicate thf3sl (TH flux)
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,thf3sl)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,1,ng,ng,0,ng,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,thf3sl)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,thf3sl)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,1,ng,ng,0,ng,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,thf3sl)

        ENDIF

    
    ENDIF


#endif

!     write(0,*) 'sfcflx deallocate ',my_rank
      
      
      deallocate ( frsea )
      deallocate(ssoil)
      deallocate(sphys)
      
!      write(0,*) 'sfcflx done'

      return

      end SUBROUTINE sfcflx



   
END MODULE FORCE_MODULE
   