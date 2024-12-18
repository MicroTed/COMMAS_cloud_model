!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        INDEX_MODULE        ////////////////////
!
!
!
!===========================================================================
 MODULE INDEX_MODULE
 
 implicit none

  integer :: imixice = 1
  integer :: icespheres = 0 ! turn ice spheres (frozen droplets) on (1) or off (0)
  integer :: frozendrops = 0 ! turn frozen drops on (1) or off (0)
  integer :: iraintypes = 0  ! turn on rain sources 
  integer :: nraintypes = 3  ! number of rain source types
  integer :: inetchargetend = 0 ! turn on net charge tendency output arrays
  integer :: ioutput_xtrachgsep = 0
  integer :: iashtypes = 0
! Define indices for microphysics.  Set default values for ZIEG scheme

  integer, parameter :: lqmx = 30
  integer, parameter :: lt = 1
  integer, parameter :: le = 0

! Constants for bin microphysics.

  integer, parameter :: nch = 22 ! Straka HCM bins
  real               :: hm(nch) = 0.0, hmvol(nch) = 0.0, hmdn(nch) = 0.0
  real, parameter    :: hjo = 3.40,hmmin = 4.7e-10 
  
  integer, parameter :: ntakrd = 34 ! 33 ! number of precip liquid bins in Takahashi
  integer, parameter :: ntakpd = 45 ! number of precip ice bins in Takahashi
  integer, parameter :: ntakid =  21 ! 10 ! 21 !  5 number of ice crystal bins in Takahashi (diameter)
  integer, parameter :: ntakit =   5 !   2 number of ice crystal bins in Takahashi (thickness)
  integer :: ntakht = 2   !  number of rimed ice types (2=graupel + hail; 1=graupel only)
  integer, parameter :: nmore  =   3 ! number of extra bin characteristics (charge, volume, water fraction)
  

! Mixing Ratio
  integer, parameter :: lv = 2
  integer, parameter :: lc = 3
  integer :: lr = 4
  integer :: lm = 0
  integer :: li = 5
  integer :: lir = 0
  integer :: lis = 0
  integer :: ls = 6
  integer :: lgl = 0
  integer :: lgm = 0
  integer :: lgh = 0
  integer :: lf = 0
  integer :: lh = 7
  integer :: lip = 0
  integer :: lhl = 0
  integer :: lhw = 0
  integer :: lfw = 0
  integer :: lsw = 0
  integer :: lhlw = 0
  integer :: lhlwlg = 0
  integer :: lhwlg = 0
  integer :: lhab = 7
  integer :: lqb = lc
  integer :: lqe = 7
  integer :: lg = 0
  integer :: lrauto = 0
  integer :: lrshed = 0
  integer :: lrmelt = 0
  integer :: lash1 = 0
  integer :: lash2 = 0
  integer :: lash3 = 0
  !  parameter ( lqb=lc, lqe=lhl )

   integer, parameter :: lqi = lc !  index of first mixing ratio (qc)

! Number concentration

  integer :: lccna = 0 ! activated ccn (default accumulation mode)
  integer :: lccnaco = 0 ! activated ccn (co = coarse)
  integer :: lccnanu = 0 ! activated ccn (nu = 'nuclei' (aitken))
  integer :: lcina = 0 ! activated ice nuclei
  integer :: lccn = 8
  integer :: lccnuf = 0
  integer :: lcin = 0  ! ice nuclei (Takahashi)
  integer :: lcng = 0  ! crystals from HM process (takahashi)

  integer :: lcn_nu = 0 ! 27 ! cloud nuclei ! need to check no conflict with other variables
  integer :: lcn_ac = 0 ! 28 ! accumulation mode
  integer :: lcn_co = 0 ! 29 ! coarse mode aerosol
  integer :: lcinp  = 0 ! 30 ! ice nuclei


  integer :: lnc = 9
  integer :: lnr = 10
  integer :: lnm = 0
  integer :: lni = 11
  integer :: lnis = 0 ! ice spheres
  integer :: lnir = 0
  integer :: lnip = 0
  integer :: lns = 12
  integer :: lngl = 0
  integer :: lngm = 0
  integer :: lngh = 0
  integer :: lnf = 0
  integer :: lnh = 13
  integer :: lnhl = 0
  integer :: lnhf = 0 ! graupel number from frozen drops (only if lf = 0)
  integer :: lnhlf = 0 ! hail number from frozen drops
  
  integer :: lnchaff
  integer :: lnox
  integer :: lnox_a
  integer :: lnox_b
  integer :: lnox_c
  integer :: lnox_d
  integer :: lnox_e
  integer :: lnox_f
  integer :: lanox
  integer :: lco
!  integer :: lnshi

! Max supersaturation
  integer :: lsat = 15
  integer :: lsati = 16
  integer :: lss = 17
  integer :: lppert = 0

! microphysics rates (special output)
  integer :: lqaut   = 0
  integer :: lqacc   = 0
  integer :: lqrevap = 0
  integer :: lqcevap = 0
  integer :: lqcond  = 0
  integer :: lqdep   = 0
  integer :: lqmelt  = 0
  integer :: lqsub   = 0
  integer :: lchlcnh = 0 ! graupel-hail conversion rate
  integer :: lchlcnf = 0 ! FD-hail conversion rate
  integer :: ldhlcnh = 0 ! graupel-hail conversion diameter
  integer :: ldhlcnf = 0 ! FD-hail conversion diameter

  integer :: lqh2hl  = 0 ! conversion from graupel to hail
  integer :: lqhl2h  = 0 ! conversion from hail to graupel (Tak bin)
  
  integer :: ld0shdrate ! production of rain drops from melting in range 0 to D1 (4.5 mm drops)
  integer :: ld1shdrate ! production of rain drops from melting in range D1 to D2 (4.5 mm drops)
  integer :: ld2shdrate ! production of rain drops from melting in range D2 to D3 (3 mm drops)
  integer :: ld3shdrate ! production of rain drops from melting in range D3 to inf (1.5 mm drops)

! fall speeds (special output)
  integer :: lmvr  = 0 ! mass-weighted rain, etc.
  integer :: lmvh  = 0
  integer :: lmvf  = 0
  integer :: lmvhl = 0
  integer :: lmvi  = 0
  integer :: lmvs  = 0
  integer :: lnvhl  = 0 ! number-weighted hail fall speed
  integer :: lzvhl  = 0 ! Z-weighted hail fall speed
  integer :: lmnvhl = 0 ! ratio of mass to number-weighted hail fall speed
  
! hail parameters
  integer :: idfw  = 0 ! mean mass frozen drop diameter
  integer :: idhw  = 0 ! mean mass graupel diameter
  integer :: idhl  = 0 ! mean mass hail diameter
  integer :: idmhl = 0 ! mass-weighted hail diameter
  integer :: idnhl  = 0 ! hail characteristic diameter
  integer :: ialphahl = 0 ! hail shape parameter
  integer :: ialphah = 0 ! graupel shape parameter
  integer :: ialphaf = 0 ! frozen drop shape parameter
  integer :: ialphar = 0 ! rain shape parameter
!  integer :: i
  
! Other output variables
  integer :: ld0   = 0 ! rain median volume diameter
  integer :: lairpress = 0 ! Air pressure
  integer :: lairtem   = 0 ! Air temperature
  integer :: lthetav   = 0 ! Virtual pot. temp.
  integer :: lssw3d = 0   ! SSw
  integer :: lssi3d = 0   ! SSi
  integer :: lrh3d  = 0   ! RHw
  integer :: lcwdia = 0   ! droplet diameter
  integer :: lcwnu  = 0   ! droplet shape parameter
  integer :: lswdia = 0   ! snow diameter
  integer :: lswmass = 0   ! snow mean mass
  integer :: lswdn = 0    ! snow density
  integer :: lswagg = 0    ! snow aggregation
  integer :: lnctak = 0
  integer :: lnrtak = 0
  integer :: lnhtak = 0
  integer :: lnhltak = 0
  integer :: lnitak = 0
  integer :: lhwdn = 0    ! graupel density
  integer :: lfwdn = 0    ! frozen drop density
  integer :: lhldn = 0    ! hail density
  integer :: lkmt  = 0
  integer :: lkht  = 0
  integer :: ltkediss = 0

  integer :: ldbzr = 0   ! rain reflectivity (dbz)
  integer :: ldbzh = 0   ! rain reflectivity (dbz)
  integer :: ldbzhl = 0   ! rain reflectivity (dbz)
  integer :: ldbzi = 0   ! rain reflectivity (dbz)
  integer :: lwvzf = 0
  integer :: lvdbz = 0
  integer :: ldbzchangeh = 0 ! increases in Z from sedimentation for graupel
  integer :: ldbzchanger = 0 ! increases in Z from sedimentation for rain

  integer :: lxrarh   = 0
  integer :: lxcrgis  = 0
  integer :: lxehw    = 0
  integer :: lxrarhl  = 0
  integer :: lxcrhis  = 0
  integer :: lxehlw   = 0
  
  integer :: lhwind = 0 ! Horizontal wind speed
  integer :: ioutput_hwind3d = 0 ! flag to output hwind
  integer :: ioutput_temC3d = 1 ! flag to output temperature in Celsius
  integer :: ioutput_pres3d = 1 ! flag to output pressure in Pa
  integer :: ioutput_thv3d  = 1 ! flag to output theta-v in Pa
  integer :: ioutput_workshop = 0 ! flag to output special workshop diagnostics
  integer :: ioutput_flshr8km = 0 ! special GLM output. Must have dx=1000 and nx and ny divisible by 8
  integer :: ioutput_snowstuff = 0 ! diagnostic snow arrays (density, mean diam)
  integer :: ioutput_takrates = 0 ! terms for Tak bin
  integer :: ioutput_mltshedsizerates = 0 ! output rates of different drop sizes from melting
  integer :: ioutput_sedmeltstuff = 0 ! output for melting and sedimentation 
  integer :: ioutput_vzf = 0 ! output refl-wgt fall speed and W-VZF
  integer :: ioutput_icedensity = 1 ! output density of graupel, hail, frozen drops
  integer :: ioutput_tkediss = 0 ! tke dissipation rate; must have iturbenhance > 0
  integer :: ioutput_dbzsedchange = 0 ! increases in Z from sedimentation
  
  integer :: ioutput_ssw = 0 ! flag to output supersaturation wrt liquid  
  integer :: ioutput_ssi = 0 ! flag to output supersaturation wrt ice  
  integer :: ioutput_rh = 0  ! flag to output relative humidity
  integer :: ioutput_ssmx = 1 ! flag to output maximum supersaturation wrt liquid  
  integer :: ioutput_cwdia = 0 ! flag to output droplet diameter
  integer :: ioutput_cnu = 0 ! flag to output droplet diameter

! Rime history

  integer :: lrtc = 0
  integer :: lrtr = 0
  integer :: lrti = 0
  integer :: lrtir = 0
  integer :: lrts = 0
  integer :: lrtgl = 0
  integer :: lrtgm = 0
  integer :: lrtgh = 0
  integer :: lrtf = 0
  integer :: lrth = 0
  integer :: lrtip = 0
  integer :: lrthl = 0
  integer :: lrshi = 0


! Space charge

  integer :: lscw = 15
  integer :: lscr = 16
  integer :: lscm = 0
  integer :: lsci = 17
  integer :: lscis = 0
  integer :: lscir = 0
  integer :: lscs = 18
  integer :: lscgl = 0
  integer :: lscgm = 0
  integer :: lscgh = 0
  integer :: lscf = 0
  integer :: lsch = 19
  integer :: lscip = 0
  integer :: lschl = 0
  integer :: lscwi = 0
  integer :: lscpi = 0
  integer :: lscni = 0
  integer :: lscpli = 0
  integer :: lscnli = 0  
  integer :: lschab = 19

  integer :: lscb = 15
  integer :: lsce = 21
  integer :: lsceq = 19

  integer, parameter :: lscmx = 100

! Particle volume

  integer :: lvi = 0
  integer :: lvs = 0
  integer :: lvgl = 0
  integer :: lvgm = 0
  integer :: lvgh = 0
  integer :: lvf = 0
  integer :: lvh = 0
  integer :: lvhl = 0

! reflectivity (6th moment)

  integer :: lzr = 0
  integer :: lzm = 0
  integer :: lzi = 0
  integer :: lzs = 0
  integer :: lzgl = 0
  integer :: lzgm = 0
  integer :: lzgh = 0
  integer :: lzf = 0
  integer :: lzh = 0
  integer :: lzhl = 0
  
  real :: alphamax = 15.
  real :: alphamin = 0.
  real :: rnumin = -0.8
  real :: rnumax = 15.0

! elec array indices
  integer, parameter :: iex      = 1
  integer, parameter :: iey      = 2
  integer, parameter :: iez      = 3
  integer, parameter :: iemag    = 4
  integer, parameter :: ipot     = 5
  integer, parameter :: ieflshn  = 6
  integer, parameter :: ieflshp  = 7
  integer, parameter :: ieinit   = 8
  integer, parameter :: icghis   = 9
  integer, parameter :: icghw    = 10
  integer, parameter :: iscnet   = 11
  integer, parameter :: iemagstag= 12

! net charge tendency arrays
  integer :: ichgtndadv = 0
  integer :: ichgtndadvsn = 0
  integer :: ichgtndmix = 0
  integer :: ichgtndnic = 0
  integer :: ichgtndsed = 0
  integer :: ichgtndion = 0
  integer :: ichgtndlgt = 0
  integer :: ichgtndave = 0
  
  integer :: iex2      = 0
  integer :: iey2      = 0
  integer :: iez2      = 0
  integer :: icgaddl  = 0
  integer :: icgnoliq = 0
  integer :: nprecip = 4
  integer :: neelec   = 12
  integer :: neelec2d = 0
  
  integer :: iprecip = 0
  integer :: iprainacc = 0
  integer :: iprainacc2 = 0
  integer :: ipwfrat = 0
  integer :: ipwfacc = 0
  integer :: iphailacc = 0
  integer :: iphailacc2 = 0
  integer :: iphailfacc = 0
  integer :: iphailfacc2 = 0
  integer :: iphailnumacc = 0
  integer :: iphailnumacc2 = 0
  integer :: iphailfnumacc = 0
  integer :: iphailfnumacc2 = 0
  integer :: iphaildiam = 0
  
  integer :: iprainaccauto = 0 ! rain from autoconversion
  integer :: iprainaccmelt = 0 ! rain from melting
  integer :: iprainaccshed = 0 ! rain from shedding
  integer :: iprainaccauto2 = 0 ! rain from autoconversion
  integer :: iprainaccmelt2 = 0 ! rain from melting
  integer :: iprainaccshed2 = 0 ! rain from shedding

  integer :: ipfdauto = 0 ! frozen drops from warm rain
  integer :: ipfdmelt = 0 ! FD from meltwater rain
  integer :: ipfdshed = 0 ! FD from shed rain
  integer :: ipfdauto2 = 0 ! frozen drops from warm rain
  integer :: ipfdmelt2 = 0 ! FD from meltwater rain
  integer :: ipfdshed2 = 0 ! FD from shed rain

  integer :: ipgracc = 0
  integer :: ipgracc2 = 0
  integer :: ipgrnumacc = 0
  integer :: ipgrnumacc2 = 0
  integer :: ipgrdiam = 0

  integer :: ipfdacc = 0
  integer :: ipfdnumacc = 0
  integer :: ipfdacc2 = 0
  integer :: ipfdnumacc2 = 0
  integer :: ipfddiam = 0

  integer :: idbzcomp = 0
  integer :: ihwindsfcmx = 0
  integer :: iwzsfcmax = 0
  integer :: iwzsfcmin = 0
  integer :: iflshr = 0
  integer :: iflshsrc = 0
  integer :: iflshr8km = 0
  integer :: iflshfed = 0
  integer :: iflshfod = 0
  integer :: iflshfodic = 0
  integer :: iflshfodcgn = 0
  integer :: iflshfodcgp = 0
  integer :: iflshfedic = 0
  integer :: iflshfedicp = 0
  integer :: iflshfedicn = 0
  integer :: iflshfedcgn = 0
  integer :: iflshfedcgp = 0
  integer :: ipnic2d = 0
  integer :: innic2d = 0
  integer :: ihailmax2d = 0
  integer :: ihailmaxk1 = 0

  integer :: nxtra   = 1

  real :: cnoh0 = 4.0e+5 
  real :: hwdn1 = 700.0

  real    :: alphar  = 0.0 ! shape parameter for ZIEG rain (only used for imurain == 1)
  real    :: alphai  = 0.0 ! shape parameter for ZIEG ice crystals
  real    :: alphas  = 0.0 ! shape parameter for ZIEG snow
  real    :: alphah  = 0.0 ! shape parameter for ZIEG graupel
  real    :: alphahl = 1.0 ! shape parameter for ZIEG hail

  real    :: dmuh    = 1.0  ! power in exponential part (graupel)
  real    :: dmuhl   = 1.0  ! power in exponential part (hail)
  
  real            :: cnu = 0.0
  real, parameter :: rnu = -0.8, snu = -0.8, cinu = 0.0
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
  
  real xnu(lc:lqmx) ! 1st shape parameter (mass)
  real xmu(lc:lqmx) ! 2nd shape parameter (mass)
  real dnu(lc:lqmx) ! 1st shape parameter (diameter)
  real dmu(lc:lqmx) ! 2nd shape parameter (diameter)
  
  real ax(lc:lqmx)
  real bx(lc:lqmx)
  real :: fx(lc:lqmx) = 1.

!
! max and min mean volumes
!
!      real :: xvcmn, xvcmx = 2.89e-13  ! min, max droplet volumes
      real xvrmn, xvrmx0  ! min, max rain volumes
      real xvsmn, xvsmx0  ! min, max snow volumes
      real xvfmn, xvfmx  ! min, max frozen drop volumes
      real xvgmn, xvgmx  ! min, max graupel volumes
      real :: xvhmn = 0., xvhmx = 0.  ! min, max graupel volumes
      real xvhmn0, xvhmx0  ! min, max hail volumes
      real xvhlmn0, xvhlmx0  ! default min, max lg hail volumes
      real :: xvhlmn = 0., xvhlmx = 0. ! min, max lg hail volumes

!      real, parameter :: cwradn = 2.5e-7, xcradmn = cwradn    ! minimum radius
      real, parameter :: cwradn = 2.0e-6, xcradmn = cwradn    ! minimum radius
      real, parameter :: cwradx = 60.e-6, xcradmx = cwradx    ! maximum radius

      real, parameter :: dhlmn0 = 0.3e-3, dhlmx0 = 40.e-3
      real, parameter :: dhmn0 = 0.3e-3
      real :: dhmn = dhmn0, dhmx = -1.
      real :: dhlmn = dhlmn0, dhlmx = -1.
      real, parameter :: dsmx0 = 10.e-3
      real :: dsmx = -1.0 ! dsmx0
      real :: xvsmx

!      parameter( xvcmn=4.188e-18 )   ! mks  min volume = 1 micron radius
      real, parameter :: xvcmn=0.523599*(2.*cwradn)**3    ! mks  min volume = 2.5 micron radius
      real, parameter :: xvcmx=0.523599*(2.*xcradmx)**3    ! mks  min volume = 2.5 micron radius

      real, parameter :: cwmasn = 1000.*xvcmn   ! minimum mass, defined by radius of 5.0e-6
      real, parameter :: cwmasx = 1000.*xvcmx   ! maximum mass, defined by radius of 50.0e-6
      real, parameter :: cwmasn5 = 1000.*0.523599*(2.*5.0e-6)**3 !  5.23e-13

      real, parameter :: xvimn=0.523599*(2.*5.e-6)**3    ! mks  min volume = 2.5 micron radius
      real, parameter :: xvimx=0.523599*(2.*1.e-3)**3    ! mks  min volume = 2.5 micron radius

      real     :: xvdmx = -1.0 ! use to set maximum rain diameter (meters); default 6.0e-3 (xvrmx0)
      real     :: xvrmx
      parameter( xvrmn=0.523599*(80.e-6)**3, xvrmx0=0.523599*(6.e-3)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvsmn=0.523599*(0.01e-3)**3, xvsmx0=0.523599*(dsmx0)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvfmn=0.523599*(0.1e-3)**3, xvfmx=0.523599*(15.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvgmn=0.523599*(0.1e-3)**3, xvgmx=0.523599*(10.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhmn0=0.523599*(dhmn0)**3, xvhmx0=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhlmn0=0.523599*(dhlmn0)**3, xvhlmx0=0.523599*(dhlmx0)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3


!      parameter ( cnoh0=4.0e+5 )
!      parameter ( hwdn1 = 700.0 )

 CONTAINS
! ---------------------------------------
 SUBROUTINE INDEX_MODULE_INIT
 ! This dummy subroutine is here because gfortran on mac would not load nonparameter 
 !  data values without a contained subroutine.
   RETURN
 END SUBROUTINE INDEX_MODULE_INIT

 END MODULE INDEX_MODULE
