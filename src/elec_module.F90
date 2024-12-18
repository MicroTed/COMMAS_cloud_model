!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        ELEC_MODULE        ////////////////////
!
!
!
!===========================================================================
 MODULE ELEC_MODULE
 
 implicit none


      
! Electrification parameters

  integer :: elec_on_time = -1     ! time (seconds) to turn on charge separation.
  integer :: elec_ramp_time = 0   ! time (interval) for linear ramp after elec_on_time 
                                   ! (i.e., linear factor on chg sep to smoothly turn on elec)
                                   ! full charging rate is achieved at time = elec_on_time + elec_ramp_time
  
  logical :: fairweather = .true.  ! true: set up background with fairweather ions/e-field
                                   ! false: set up background with zero ions.
  logical :: largeion = .false.   ! true: use large ion category 
                                   ! false: does not use large ion category
  integer :: icorona = 1
  integer :: iondriftorder = 6  ! 6th order Crowley or 1st order upstream
  real    :: emaxc   = 5.0e3
  
  integer :: nic_noliq = 0  ! switch for NIC without riming (MST, 2006)
  real    :: qc_noliq = 0.1 ! charge per collision (in fC) for non-riming NIC.
  
  integer :: ieopt = 0       ! = 0 no additional electrification options
                             ! = 1 allow charge separation during melting based on Drake (1968)
                             ! = 2 allow charge separation during depositional growth based on Dong and Hallett (1992)
                             ! = 3 Canosa and List (1993)
  logical :: drake = .false. 
  integer :: ihallettsaund = 0 ! Hallett and Saunders (1979) ice mult. charging
  real    :: xhallettsaund = 5.e-14 ! charge separation per ejected crystal
  
  integer :: ibdcnd = 0, jbdcnd = 0  ! lateral BCs; for future use
  integer :: kbdcnd = 1 ! set potential at top & bottom
  integer :: kdbcdn = 0 ! not used?
  
  integer :: maxbzone = 20 ! maximum number of lateral grid points in extended potential domain.

  integer :: nliter = 200, ipotslv = 1, ilghtest = 0
  integer :: ilight = 2 ! 1 = bulk cylinder scheme; 2 = DBM branched scheme; 3 = MSZ
  integer :: nhelec = 999999, iscreen = 0
  
  integer :: nonigrd = 0          ! nonigrd = 0 : UMIST schemes (S91, SP98, etc)
                                  ! nonigrd = -1 : Takahashi variants
                                  ! nonigrd = 2  : Gardiner/Ziegler
  integer :: ftauopt = 1         ! GZ charging method for shifting Trev in the ftau function: 
                                 !  1 = Ziegler/Straka scaled 
                                 !  2 = Left-right shift of the function (doesn't change magnitude)
  integer :: isaund = -999
  real    :: rgard   = 9
  real    :: rarfac  = 1.0
  integer :: iraropt = 1  ! 1 = mass-wgt; 2 = num-wgt; 3 = area-wgt; 4 = diam-wgt; 5 = avg. of mass and num-wgt Vt
  real    :: rgard1
  real    :: trever = -15.0
  integer :: isctemopt = 2  ! option for roll-off factor for nonigrd=2: 1=old quadratic; 2 = cosine
  real    :: nic_min_temp = -33. ! Low temperature cut-off for cosine roll-off (nonigrd=2)
  
  real    :: cidiamin = 0.0  ! minimum crystal size for charging
  
  real    :: erbnd  = 0.01
  real    :: fdgt   = 1.0
  real    :: costhe = 0.35

  real    :: delqnsa = -50.0e-14
  real    :: delqxsa =  50.0e-14
  real    :: delqnsb = -2.0e-14
  real    :: delqxsb =  2.0e-14
  real    :: delqnia = -50.0e-14
  real    :: delqxia =  50.0e-14
  
  real    :: delqnra = 0.0
  real    :: delqxra = 0.0
  
  real    :: delqxw =  1.0e-10
  real    :: delqnw = -1.0e-10
  
  real    :: scxacymax = 3000.e-12 ! Maximum charging rate magnitude (C s-1 m-3)
      
  real    :: scippmx = 20.0e-12
  real    :: scwppmx = 20.0e-12

  real    :: esctot = 1.0e-13               ! threshold for charge budget reporting
  
  real    :: tindmn = 233, tindmx = 298.0  ! min/max temperatures to allow inductive collisional charging
  
  integer :: energymethod = 1  ! 1 = pot*rho; 2 = E^2
  
  real    :: crgfac = 1.0
  
  logical :: do_ionatt_conduction = .true. ! Turn on/off the conduction part of ion attachment
  logical :: do_ionatt_diffusion = .true. ! Turn on/off the diffusion part of ion attachment

  real    :: ionatt_conduction_factor = 1.0 ! For testing only!
  real    :: ionatt_diffusion_factor = 1.0 ! For testing only!

!
!   BoxMG control parameters (for MPI)
!


  integer ::  bmg_cycletype = 1 ! 0=FMG, 1=n-cycle
  integer ::  bmg_maxiter = 2
  integer ::  bmg_sep     = 0  ! whether to use separable mod or not: 0=BMG default, 1=simplified 1-D coeffs for fine grid
  integer ::  bmg_cgsolve   = 1  ! 0 = LU direct solver; 1 = boxmg_serial solve
  real    ::  bmg_tol     = 1.d-3
  
  integer :: nlxdg, nlydg, nlzdg, ngxdg, ngydg, ngzdg
  integer :: igsdg, jgsdg, kgsdg
  integer :: igsdg0=0, jgsdg0=0, kgsdg0=0
  
  integer :: NOGdg
  integer :: NOGcdg
  integer :: NFdg
  integer :: NCbmgdg
  integer :: NCIdg
  integer :: NSOdg
  integer :: NSORdg
  integer :: NCBWdg
  integer :: NCUdg
  integer :: NMSGidg
  integer :: NMSGrdg

  integer :: nlxlg, nlylg, nlzlg, ngxlg, ngylg, ngzlg
  integer :: igslg, jgslg, kgslg
  integer :: igslg0=0, jgslg0=0, kgslg0=0
  
  integer :: NOGlg
  integer :: NOGclg
  integer :: NFlg
  integer :: NCbmglg
  integer :: NCIlg
  integer :: NSOlg
  integer :: NSORlg
  integer :: NCBWlg
  integer :: NCUlg
  integer :: NMSGilg
  integer :: NMSGrlg

! BoxMG arrays for dynamics grid (dg):

  integer, parameter :: ip_RESdg = 3
  integer, parameter :: ip_SORdg = 5
  integer, parameter :: ip_CIdg  = 6
  integer, parameter :: ip_iGdg  = 10
  integer, parameter :: ip_MSGdg = 14
  integer, parameter :: ip_MSG_BUFdg = 16
  integer, parameter :: id_BMG3_SETUPdg = 52
  integer, parameter :: id_BMG3_MAX_ITERSdg = 61


      integer, parameter :: NBMG_pWORKdg = 18
      integer, parameter :: NBMG_InWORKdg = 7
      
      INTEGER   BMG_pWORKdg(NBMG_pWORKdg)
      LOGICAL   BMG_InWORKdg(NBMG_InWORKdg)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSIdg, pSRdg 


      integer, parameter :: NBMG_iPARMSdg = 72
      integer, parameter :: NBMG_rPARMSdg = 35
      integer, parameter :: NBMG_IOFLAGdg = 44
      
      INTEGER   BMG_iPARMSdg(NBMG_iPARMSdg)
      REAL      BMG_rPARMSdg(NBMG_rPARMSdg)
      LOGICAL   BMG_IOFLAGdg(NBMG_IOFLAGdg)

      INTEGER   NFmdg, NOGmdg, NSOmdg
      INTEGER   NBMG_iWORKdg, NBMG_rWORKdg
      INTEGER   NBMG_iWORK_PLdg, NBMG_rWORK_PLdg
      INTEGER   NBMG_iWORK_CSdg, NBMG_rWORK_CSdg



      REAL, allocatable, dimension(:) :: Qdg, QFdg, SOdg

      INTEGER, allocatable, dimension(:)  :: BMG_iWORKdg
      REAL, allocatable, dimension(:)  :: BMG_rWORKdg
  
      INTEGER, allocatable, dimension(:) ::  BMG_iWORK_PLdg
      REAL,  allocatable, dimension(:) ::  BMG_rWORK_PLdg
  
      INTEGER, allocatable, dimension(:)  ::   BMG_iWORK_CSdg
      REAL,  allocatable, dimension(:)  ::   BMG_rWORK_CSdg

      INTEGER, allocatable, dimension(:) ::  pMSGdg, pMSGSOdg


      INTEGER   NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg

      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRIDdg

      INTEGER, parameter :: NBMG_MSG_pGRIDdg = 8
      INTEGER  BMG_MSG_pGRIDdg(NBMG_MSG_pGRIDdg)

! BoxMG arrays for lightning grid (lg):

      INTEGER   BMG_pWORKlg(NBMG_pWORKdg)
      LOGICAL   BMG_InWORKlg(NBMG_InWORKdg)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSIlg, pSRlg 

      INTEGER   BMG_iPARMSlg(NBMG_iPARMSdg)
      REAL      BMG_rPARMSlg(NBMG_rPARMSdg)
      LOGICAL   BMG_IOFLAGlg(NBMG_IOFLAGdg)

      INTEGER   NFmlg, NOGmlg, NSOmlg
      INTEGER   NBMG_iWORKlg, NBMG_rWORKlg
      INTEGER   NBMG_iWORK_PLlg, NBMG_rWORK_PLlg
      INTEGER   NBMG_iWORK_CSlg, NBMG_rWORK_CSlg

      REAL, allocatable, dimension(:) :: Qlg, QFlg, SOlg

      INTEGER, allocatable, dimension(:)  :: BMG_iWORKlg
      REAL, allocatable, dimension(:)  :: BMG_rWORKlg
  
      INTEGER, allocatable, dimension(:) ::  BMG_iWORK_PLlg
      REAL,  allocatable, dimension(:) ::  BMG_rWORK_PLlg
  
      INTEGER, allocatable, dimension(:)  ::   BMG_iWORK_CSlg
      REAL,  allocatable, dimension(:)  ::   BMG_rWORK_CSlg

      INTEGER, allocatable, dimension(:) ::  pMSGlg, pMSGSOlg

      INTEGER   NBMG_MSG_iGRIDmlg, NBMG_MSG_iGRIDlg

      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRIDlg

      INTEGER  BMG_MSG_pGRIDlg(NBMG_MSG_pGRIDdg)

!
! Lightning params:
!
! for ilight = 1
  real    :: chgthr   = 0.1e-9  ! threshold charge density for ilight = 1
  real    :: elgt1    = 0.7     ! fractional reduction in charge in excess of chgthr
  real    :: elght1   = 150000. ! threshold e-field to init lightning
  real    :: lightrad = 8000.0  ! radius (meters) of influence (half-width of square) for bulk lightning (ilight=1).  
                                ! Set < 0 for whole domain (not recommended, though).
                                ! lightrad is turned into a number of gridpoints as nr = Nint(lightrad/dx)
  real    :: lightextendmsz = 6000. ! when comparing footprints of positive and negative discharge region, remove points 
                                    ! that are beyond this distance of overlap. This is to prevent one region from having a
                                    ! much larger area than the other
! for ilight = 2
  real    :: elgtthx  = 100.0e3
  real    :: elgtthn  =  90.0e3
  real    :: elgtdel  =  10.0e3
  real    :: elgtfdel = 0.9
  real    :: einitmax = 180.0e3   ! imposed max for height-dependent initiation threshold (ibrkd > 1)
  real    :: overvolt = 1.0
  integer :: ibrkd    = 4         ! (=1 for constant elgtthx, =2 for breakeven field, =4 for new breakdown profile (about 1.5*old))
  integer :: iusecl   = 0         ! (=0 for old init threshold, =1 to check in vertical for minimum characteristic length)
  integer :: ieint    = 1         ! (1=use const. eint; 2=use feint)
  real    :: eint     = 100.0     ! channel internal field for propagation
  real    :: feint    = 0.01      ! fraction of breakdown threshold for internal field (if ieint=2)
  real    :: fecrit   = 0.47      !
  real    :: fefac(2) = 1.0       !
  integer :: iseed    = 789
  integer :: intfg    = 1        ! how often to run 'full grid' relaxation (every intfg steps)
  integer :: intfgmpi = 25       ! how often to call BOXMG solver during lightning (MPI only)
  integer :: mdel     = 15
  integer :: ibal     = 2
  real    :: cgfr     = 1.0
  real    :: cgthres  = 2500.
  integer :: ibranch  = 0
  integer :: maxpicks = 1 ! 10
  integer :: pickmaskwidth = 2
  integer :: multipickmod = 20 ! add an extra pick for each number of added steps
  integer :: nohold   = 0  ! if > 0, then always adjust ref. pot. after each step with nohold number of iterations
                           !    = 0, then normal ibal2 with limit on reference shift
                           ! if < 0, then normal ibal2 but no limit on reference shift
  real    :: dslight  = 500.
  real    :: dslightz = 500.
  integer :: lightintx = -1
  integer :: lightintz = -1
  integer :: l2nde   = 2       ! l2nde=1 for 1st order, =2 for 2nd order electric field next to channel
  integer :: intsormpi = 2     ! interval for calling strzsormpi right after new point
  real    :: soreps  = 5.0e-3  ! convergence epsilon for strzsor
  real    :: soreps1 = 5.0e-3  ! convergence epsilon for strzsor on interval 'intsormpi'
  integer :: ltgeta  = 1       ! power law for propagation probability
  integer :: igrid   = 0       ! 0 = interpolate to points between grid centers
                               ! 1 = subdivide gridvolumes so that lightnign volumes
                               !     are within original volume and inherit the same value
                               !     i.e., values are injected, not interpolated, which can
                               !     enhance gradients but preserves total charge.

  integer :: initwire = 0   ! 0 = no lightning initialization by a wire
                            ! 1 = initiate a wire at a certain time and location
                            ! 2 = initiate a wire at a certain E threshold (fraction of E-init)
  integer :: ixinitwire,jyinitwire,kzinitwire,itinitwire ! x,y,z,t for wire initiation
  real    :: efracinitwire = 1.0 ! for initwire=2
  real    :: wirelength = 1000.   ! length of wire (meters) extended vertically
  
  integer ::             ieswi = 1, ieswir = 1, ieswip = 1, ieswc = 1, ieswr = 0
  integer :: ieglsw = 1, iegli = 1, ieglir = 1, ieglip = 1, ieglc = 1, ieglr = 0
  integer :: iegmsw = 1, iegmi = 1, iegmir = 1, iegmip = 1, iegmc = 1, iegmr = 0
  integer :: ieghsw = 1, ieghi = 1, ieghir = 1, ieghip = 1, ieghc = 1, ieghr = 0
  integer :: iefwsw = 1, iefwi = 1, iefwir = 1, iefwip = 1, iefwc = 1, iefwr = 0
  integer :: iehwsw = 1, iehwi = 1, iehwir = 1, iehwip = 1, iehwc = 1, iehwr = 0
  integer :: iehlsw = 1, iehli = 1, iehlir = 1, iehlip = 1, iehlc = 1, iehlr = 0
  integer :: iehlis = 1, iehis = 1, iesis = 1, iefwis = 1
  
  integer :: irand = 0
  integer :: isa   = 0
  integer :: ixst0  = -1, jyst0 = -1, kzst0 = -1, itst0 = 0 ! x,y,z,t for choosing an initiation point
  
! for ilight = 3 (MSZ parameterization)
      real    :: elgtfestopcg = 0.15         ! fraction of Emag for stopping channel for CG
      real    :: zgrnd = -1. ! 1500.0              ! threshold height to declare channel end at ground
      real    :: tgrnd = 268.15              ! threshold temperature to declare channel end at ground

   logical :: hstretch = .false. ! flag for horizontal stretching

!
! fair weather field 
!
  integer, parameter :: nzfair = 1000
  real               :: ezfair  (nzfair)
  real               :: potfair (nzfair)
  real               :: rhofair (nzfair)
  real               :: ezfairw (nzfair)

! other shared variables

!  real            :: eztop
  
  integer, parameter :: nzmax1 = 1000, nzmax1lgt = 6*nzmax1
  real               :: z1d2(nzmax1,4)
  real, allocatable  :: x1d2(:) ! dx on boxmg solver grid
  real, allocatable  :: y1d2(:) ! dy on boxmg solver grid
  real, allocatable  :: x1d2lgt(:) ! dx on boxmg solver lightning grid
  real, allocatable  :: y1d2lgt (:) ! dy on boxmg solver lightning grid
  real               :: dlz ! average dz for extended potential domain
  integer            :: nz1d

  real               :: z1d2lgt(nzmax1lgt,4)
  real               :: dlzlgt
  integer            :: nz1dlgt
  real               :: ezfairlgt  (nzmax1lgt)
  real               :: potfairlgt (nzmax1lgt)
  real               :: rhofairlgt (nzmax1lgt)
  real               :: ezfairwlgt (nzmax1lgt)



  CONTAINS
  
  SUBROUTINE ELEC_MODULE_INIT
  RETURN
  END SUBROUTINE ELEC_MODULE_INIT
!
 END MODULE ELEC_MODULE
