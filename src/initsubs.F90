MODULE INIT_MODULE

!-----------------------------------------------------------------------------
!
! SOUNDING PARAMETERS
!
!-----------------------------------------------------------------------------
  implicit none
  
  character(LEN = 100) :: sndfile   = ' '
  real                 :: xsnd
  real                 :: ysnd
  integer              :: sndtype   = 1
!  integer              :: inittype  


! BBLE parameters

  integer, parameter :: maxbub   = 100
  integer            :: nbble    = 1
  integer            :: ibbleseed = -1
  integer            :: bbletype = 1
  integer            :: bblsp(maxbub) = 0  ! 'special' bubble parameters
  real               :: tbble(maxbub) = 0.0
  real               :: xrad (maxbub) = 0.0
  real               :: yrad (maxbub) = 0.0
  real               :: zrad (maxbub) = 0.0
  real               :: xcntr(maxbub) = 0.0
  real               :: ycntr(maxbub) = 0.0
  real               :: zcntr(maxbub) = 0.0
  real               :: dthranpert = 0.1  ! default maximum temperature pert for bbletype = 2


! BBLE parameters


! WK PARAMS

  integer              :: wtype  = 1
  real                 :: Us     = 10.
  real                 :: Usmv   = 0.     ! constant subtracted from U
  real                 :: Uz     = 2500.
  real                 :: Usl    = 0.     ! magnitude of lower shear layer (below Ubase)
  real                 :: Ubase  = 0.     ! Base of upper shear layer
  real                 :: psfc   = 100000.
  real                 :: tsfc   = 300.
  real                 :: qsfc   = 14.
  real                 :: wkexponent = 1.25
  real                 :: chgtwrf = 2250.
  real                 :: ztr     = 12000.
  real                 :: zbl     = -1.
  real                 :: p_bl    = -1.
  real                 :: thbl    = -1. ! 308
  real                 :: dthdzbl = -1.
  real                 :: dzcaplayer =  0. ! depth of inhibition layer
  real                 :: dthcaplayer = 0.  ! delta-theta over the inhibition layer
  real                 :: thtr    = 343.
  real                 :: ttr     = 213.
  real                 :: rhmaxin = 0.98 ! max RH for input sounding
  real                 :: rhmax  = 0.90
  real                 :: rhmax2 = 0.25
  real                 :: zrhmax2 = -100.
  real                 :: zrhdel  = 700. ! distance for transition from rhmax to rhmax2
  integer              :: ntr    = 1      ! number of iterations for WK profile to get right qv-pi balance
  logical              :: rebalance = .false. ! whether to update PI after bubble init
  logical              :: moisten   = .false. ! whether to moisten bubble to keep same RH as environment
  real                 :: wkumax1 = 7.0, wkumax2 = 31.0 ! Dennis and Kumjian (DK2017) params
  real                 :: wkvmax = 7.0 ! DK2017
  real                 :: wkh1 = 2000., wkh2 = 6000. ! hodo height params (DK2017)
  real                 :: ccn_height1 = 1.e10, ccn_upperair = -1. ! height to start low ccn value

! LOW LEVEL PARAMS

  integer              :: shape = 0
  real                 :: dudz0 = 0.0
  real                 :: dudz1 = 0.0
  real                 :: dudz2 = 0.0
  real                 :: z0    = 99000.
  real                 :: z1    = 98000.
  real                 :: z2    = 97000.

! Extra wind control
  integer              :: imodw1 = 0, imodw2 = 0, imodw3 = 0, ich = 0
  real                 :: ushrl1=0.0,ushrlb1=0.0,ushrlt1=0.0,vshrl1=0.0,vshrlb1=0.0,vshrlt1=0.0
  real                 :: ushrl2=0.0,ushrlb2=0.0,ushrlt2=0.0,vshrl2=0.0,vshrlb2=0.0,vshrlt2=0.0
  real                 :: ushrl3=0.0,ushrlb3=0.0,ushrlt3=0.0,vshrl3=0.0,vshrlb3=0.0,vshrlt3=0.0

  double precision, allocatable :: z1dinit(:,:)
  

! Tropical Cyclone initialization
  integer              :: bogusvortex = 0  ! 1 = Old Kanak vortex; 2 = WRF Rotunno-Emanuel balanced vortex; 3 = Nolan init
  real                 :: vorthgt = 15000.
  real                 :: zetamax = 13.e-4
  real                 :: radzero = 280000.
  real                 :: timint = 0.
  real                 :: vtmax = 15.0   ! tang. wind at radius of max winds (WRF TC init)
  real                 :: vtrad = 82500.  ! radius of max winds (WRF TC init)
  real                 :: r0  = 412500.   ! outer radius (WRF TC init)
  
  real                 :: vortexp = 0.5 ! Stern and Nolan modified Rankine vortex
  real                 :: vort_rcut = 600000 ! cutoff radius
  real                 :: vort_zmax = 1500.
  real                 :: vort_beta = 2.0
  real                 :: vort_lz   = 3175.
  real                 :: vort_lz2  = -1.
  integer              :: iadd_wind_perts = 0, ipertseed = 344
  real                 :: pert_base = 0.1

  real    :: prg0     = 101100.
  real    :: thg0     = 301.7 ! initial SST
  real    :: thg1     = 301.7 ! target SST
  real    :: isstbeg  = 0     ! time to start transitioning from thg0 to thg1
  real    :: isstend  = 0     ! time to end transitioning from thg0 to thg1
  real    :: uspeed   = 0.
  real    :: timlndst = 1.e9
  integer :: isfcl = 0

! EXTERNAL MESOSCALE DOMAIN PARAMS

  integer              :: inhom = 0
  integer              :: numinhom = 0
  integer, parameter   :: maxsnd  = 10
  integer              :: nxmeso = 100
  integer              :: nymeso = 100
  integer              :: nzmeso           ! nzmeso = nz - 1 as set in commas.F90
  integer              :: nmtimes = 1
  integer              :: nmint   = 1000
  integer              :: nvbls   = 4
  real                 :: dxmeso = 5000.0
  real                 :: dymeso = 5000.0
  real                 :: dzmeso = 200.0
  character(LEN = 100) :: sndinhom(maxsnd)  = ' '
  real                 :: xsndinhom(maxsnd)
  real                 :: ysndinhom(maxsnd)
  integer              :: ivelh1d = 0
  integer              :: wnditp = 0
  integer              :: irot = 0
  real                 :: rotdeg = 0.0
  integer              :: khomog           ! khomog either read in or set in commas.F90
  real                 :: sf_north = 100000.0
  real                 :: sf_south = 60000.0
  integer              :: dencorr = 0
  real                 :: cden = 0.0
  real                 :: hden = 0.0
  real                 :: yshiftmag = 0.0
  real                 :: yshiftdepth = 700.
 

!
!.... mesoscale domain array abmeso(nxmeso,nymeso,nzmeso,nvbls,nmtimes)
!

    real, allocatable  :: abmeso(:,:,:,:,:)

!
!.... CLZ (3-8-13): arrays now dynamically allocated to avoid stack problem and save values
!

    real, allocatable  :: athet(:,:,:,:)
    real, allocatable  :: aqrat(:,:,:,:)
    real, allocatable  :: au(:,:,:,:)
    real, allocatable  :: av(:,:,:,:)

!     mesoscale work arrays lag_thet, lag_qrat, lag_u, and lag_v

    real, allocatable    :: lag_thet(:,:,:)
    real, allocatable    :: lag_qrat(:,:,:) 
    real, allocatable    :: lag_u(:,:,:)
    real, allocatable    :: lag_v(:,:,:)

!
!.... CLZ (3-12-13): mesoscale heterogeneous environmental sounding arrays
!

    real, allocatable    :: tz_meso(:,:)
    real, allocatable    :: qz_meso(:,:)
    real, allocatable    :: u1d_meso(:,:)
    real, allocatable    :: v1d_meso(:,:)




  NAMELIST /homog_init/ nbble,      &
                        ibbleseed,  &
                        bbletype,   &
                        dthranpert, &
                        bblsp,      &
                        tbble,      &
                        xrad,       &
                        yrad,       &
                        zrad,       &
                        xcntr,      &
                        ycntr,      &
                        zcntr,      &
                        sndfile,    &
                        xsnd,       &
                        ysnd,       &
                        sndtype,    &
                        wtype,      &
                        ntr,        &
                        rebalance,  &
                        moisten,    &
                        Us,         &
                        Usmv,       &
                        Uz,         &
                        Usl,        &
                        Ubase,      &
                        psfc,       &
                        tsfc,       &
                        wkexponent, &
                        chgtwrf,    &
                        ztr, thtr, ttr, &
                        zbl, p_bl, thbl, dthdzbl, &
                        dzcaplayer, dthcaplayer, &
                        qsfc,       &
                        rhmaxin,    &
                        rhmax,      &
                        rhmax2,     &
                        zrhmax2,    &
                        zrhdel,    &
                        shape,      &
                        dudz0,      &
                        dudz1,      &
                        dudz2,      &
                        z0,         &
                        z1,         &
                        z2,         &
                        inhom,      &
                        numinhom,   &
                        sndinhom,   &
                        xsndinhom,  &
                        ysndinhom,  &
                        ivelh1d,    &
                        wnditp,     &
                           irot,    &
                         rotdeg,    &
                         khomog,    &
                        nxmeso,     &
                        nymeso,     &
                        nzmeso,     &
                        nmint,      &
                        nmtimes,    &
                        nvbls,      &
                        sf_north,   &
                        sf_south,   &
                        dencorr,    &
                        yshiftmag, yshiftdepth, &
                        cden,       &
                        hden,       &
                        dxmeso,     &
                        dymeso,     &
                        dzmeso,     &
                        bogusvortex, vorthgt, zetamax, radzero,timint, &
                        r0, vtmax, vtrad,    &
                        vortexp,vort_rcut,vort_zmax,vort_beta,vort_lz,vort_lz2, &
                        iadd_wind_perts, ipertseed, pert_base, &
                        wkumax1, wkumax2, wkvmax, wkh1, wkh2,ccn_height1, ccn_upperair, &
                        imodw1, imodw2, imodw3, ich,                                     &
                        ushrl1,ushrlb1,ushrlt1,vshrl1,vshrlb1,vshrlt1,                   &
                        ushrl2,ushrlb2,ushrlt2,vshrl2,vshrlb2,vshrlt2,                   &
                        ushrl3,ushrlb3,ushrlt3,vshrl3,vshrlb3,vshrlt3






END MODULE INIT_MODULE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


MODULE GRIDPARAM_MODULE

!-----------------------------------------------------------------------------
!
! GRID PARAMETERS      
! 
!-----------------------------------------------------------------------------
  
  USE COMMASMPI_MODULE, only: nxend,nyend,nzend
  
  implicit none
  
  integer              :: nx         = 51
  integer              :: ny         = 51
  integer              :: nz         = 43
  real                 :: dt         = 6.0
  integer              :: nsmall     = 6
  integer              :: n_middle   = 101
  integer              :: nbndlyr    = 0
  real                 :: dx_stretch = -1.0
  real                 :: dy_stretch = -1.0
  real                 :: dz_stretch = -1.0
  real                 :: max_stretch = 1.1    ! maximum stretch factor
  real                 :: dzmax      = 700.
  real                 :: dxmax      = 10000.
  real                 :: dymax      = 10000.
  real                 :: dzmaxtop   = 700.    ! max upper level dz
  real                 :: rtop       = 1.0     ! upper level stretch factor
  real                 :: ztopstr    = 100000. ! height to start upper level stretching
  real                 :: xdomain    = 50000.
  real                 :: ydomain    = 50000.
  real                 :: zdomain    = 18000.
  real                 :: x_sw_loc   = 0.0
  real                 :: y_sw_loc   = 0.0
  real                 :: lat        =  35.0
  real                 :: lon        = -95.0
  real                 :: hgt        = 0.0

! non-namelist variables
  real                 :: dx
  real                 :: dy
  real                 :: dz
                        
  NAMELIST /gridn/nx,         &     
                  ny,         &     
                  nz,         &      
                  dt,         &     
                  nsmall,     &      
                  xdomain,    &     
                  ydomain,    &     
                  zdomain,    &     
                  x_sw_loc,   &
                  y_sw_loc,   &
                  lat,        &
                  lon,        &
                  hgt,        &
                  dx_stretch, &
                  dy_stretch, &
                  dz_stretch, &
                  max_stretch, &
                  dzmax,      &
                  dxmax,dymax,  &
                  dzmaxtop,   &
                  rtop,       &
                  ztopstr,    &
                  n_middle,   &
                  nbndlyr,    &
                  nxend,      &
                  nyend,      &
                  nzend

END MODULE GRIDPARAM_MODULE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! SUBROUTINE INIT_GRID initializes the grid constants and parameters
!
!-----------------------------------------------------------------------
! APS JUN07: Note that there are some crude vertical MPI statements in here
!            that work as is only for ktile=nzend (i.e., no MPI in the vert).  
!            First need to build a vertical tile communication system.  
!            Then come back and fix the vertical integrations, particularly 
!            for sounding subs.
!
 SUBROUTINE INIT_GRID(gd) !,                                           &
!                      dt, nsmall, dx, dy, dz,                       &
!                      n_middle, nbndlyr,                           &
!                      dx_stretch,  dy_stretch, dz_stretch, dzmax,   &
!                      x_sw_loc, y_sw_loc, lat, lon, hgt, dzmaxtop, rtop, ztopstr)

  USE GRID_MODULE
  USE FILE_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE INDEX_MODULE
  USE GRIDPARAM_MODULE
  USE INIT_MODULE
  USE COMMASMPI_MODULE
  USE ELEC_MODULE, only : iseed

  implicit none

!  include 'param.h'

! Passed variables

  TYPE(GRID) :: gd

!  real    :: dx, dy, dz
!  real    :: dx_stretch, dy_stretch, dz_stretch, dzmax
!  real    :: x_sw_loc, y_sw_loc, lat, lon, hgt
!  integer :: dt
!  integer :: n_middle, nbndlyr, nsmall
!  real    ::  dzmaxtop, rtop, ztopstr

! Local vars

!  integer, pointer :: nx, ny, nz
  integer          :: i, j, k
  integer          :: dum

  real, pointer :: gc(:)
  real, pointer :: ge(:)
  real, pointer :: dgc(:)
  real, pointer :: dge(:)
  real, allocatable :: gctemp(:)
  real, allocatable :: getemp(:)

#ifdef OPENDX
  real, pointer :: xcdx(:,:)
  real, pointer :: xedx(:,:)

  real, pointer :: ycdx(:,:)
  real, pointer :: yedx(:,:)

  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)
  
  real xc1, xe1, yc1, ye1, zc1, ze1
#endif
  real zheight
  
  integer  :: gcindex, ngs
!  external zheight

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_GRID BEGIN')

!-----------------------------------------------------------------------------
! GET GRID ATTRIBUTES SPECIFYING DIMENSIONS

!  CALL GET_ATTRIBUTE(gd,'NX',nx)
!  CALL GET_ATTRIBUTE(gd,'NY',ny)
!  CALL GET_ATTRIBUTE(gd,'NZ',nz)

!--------------------------------------------------------------------------------------
! DEFER TESTING OF HORIZONTAL STRETCHING FOR MPI

!#ifdef MPI
!    if (dx_stretch .ge. 0. .OR. dy_stretch .ge. 0.) then
!     write(luno,*) 'INIT_GRID: HORIZ STRETCHING NOT TESTED FOR MPI YET'
!     write(luno,*) 'INIT_GRID: SET DX_STRETCH AND DY_STRETCH TO -1.0'
!     STOP
!    endif
!#endif

!--------------------------------------------------------------------------------------
  IF( dx_stretch .lt. 0.0 ) THEN
   dx_stretch = dx
   dxmax = dx
  ENDIF
  IF( dy_stretch .lt. 0.0 ) THEN
   dy_stretch = dy
      dymax = dy ! Max(dx,dy)
!   n_middle   = nx 
  ENDIF

  IF( dz_stretch .lt. 0.0 ) THEN
     dz_stretch = dz
     dzmax      = dz
  ENDIF

#ifdef MPI
  write(luno,*) 'XDOMAIN(RECALC FOR SUBDOM) = ', dx*(nx-1)
  write(luno,*) 'NX / DX / DX_MID           = ', nx, dx, dx_stretch
  write(luno,*)
  write(luno,*) 'YDOMAIN(RECALC FOR SUBDOM) = ', dy*(ny-1)
  write(luno,*) 'NY / DY / DY_MID           = ', ny, dy, dy_stretch
  write(luno,*)
  write(luno,*) 'ZDOMAIN(RECALC FOR SUBDOM) = ', dz*(nz-1)
  write(luno,*) 'NZ / DZ / DZ_BOT           = ', nz, dz, dz_stretch
  write(luno,*)
  write(luno,*) 'NO. OF PTS FOR INNER GRID  = ', n_middle
  write(luno,*)
  write(luno,*) 'TIME STEP (SEC)         DT = ', dt
  write(luno,*)
#else
  write(luno,*) 'XDOMAIN                    = ', dx*(nx-1)
  write(luno,*) 'NX / DX / DX_MID           = ', nx, dx, dx_stretch
  write(luno,*)
  write(luno,*) 'YDOMAIN                    = ', dy*(ny-1)
  write(luno,*) 'NY / DY / DY_MID           = ', ny, dy, dy_stretch
  write(luno,*)
  write(luno,*) 'ZDOMAIN                    = ', dz*(nz-1)
  write(luno,*) 'NZ / DZ / DZ_BOT           = ', nz, dz, dz_stretch
  write(luno,*)
  write(luno,*) 'NO. OF PTS FOR INNER GRID  = ', n_middle
  write(luno,*)
  write(luno,*) 'TIME STEP (SEC)         DT = ', dt
  write(luno,*)
#endif
!--------------------------------------------------------------------------
! Fill in constants
 
  IF( .not. SET_VARIABLE(gd,'NX',        nx)            ) write(6,*) 'INIT_GRID:  Problem setting NX'
  IF( .not. SET_VARIABLE(gd,'NY',        ny)            ) write(6,*) 'INIT_GRID:  Problem setting NY'
  IF( .not. SET_VARIABLE(gd,'NZ',        nz)            ) write(6,*) 'INIT_GRID:  Problem setting NZ'
  IF( .not. SET_VARIABLE(gd,'NG',        ng)            ) write(6,*) 'INIT_GRID:  Problem setting NG'
#ifdef MPI
  IF( .not. SET_VARIABLE(gd,'NXEND',     nxend)         ) write(6,*) 'INIT_GRID:  Problem setting NXEND'
  IF( .not. SET_VARIABLE(gd,'NYEND',     nyend)         ) write(6,*) 'INIT_GRID:  Problem setting NYEND'
  IF( .not. SET_VARIABLE(gd,'NZEND',     nzend)         ) write(6,*) 'INIT_GRID:  Problem setting NZEND'
#else
  IF( .not. SET_VARIABLE(gd,'NXEND',     nx)            ) write(6,*) 'INIT_GRID:  Problem setting NXEND'
  IF( .not. SET_VARIABLE(gd,'NYEND',     ny)            ) write(6,*) 'INIT_GRID:  Problem setting NYEND'
  IF( .not. SET_VARIABLE(gd,'NZEND',     nz)            ) write(6,*) 'INIT_GRID:  Problem setting NZEND'
#endif
  IF( .not. SET_VARIABLE(gd,'NSMALL',    nsmall)        ) write(6,*) 'INIT_GRID:  Problem setting NSMALL'
  IF( .not. SET_VARIABLE(gd,'BCX',       bcx)           ) write(6,*) 'INIT_GRID:  Problem setting BCX'
  IF( .not. SET_VARIABLE(gd,'BCY',       bcy)           ) write(6,*) 'INIT_GRID:  Problem setting BCY'
  IF( .not. SET_VARIABLE(gd,'BCZ',       bcz)           ) write(6,*) 'INIT_GRID:  Problem setting BCZ'
  IF( .not. SET_VARIABLE(gd,'N_INNER',   n_middle)      ) write(6,*) 'INIT_GRID:  Problem setting N_INNER'
  IF( .not. SET_VARIABLE(gd,'N_BNDLYR',  nbndlyr)       ) write(6,*) 'INIT_GRID:  Problem setting N_BNDLYR'
  IF( .not. SET_VARIABLE(gd,'HOLE_FILL', hole_fill)     ) write(6,*) 'INIT_GRID:  Problem setting HOLE_FILL'
  IF( .not. SET_VARIABLE(gd,'AUTOCONV',  autoconversion)) write(6,*) 'INIT_GRID:  Problem setting AUTOCONV'
  IF( .not. SET_VARIABLE(gd,'IPCONC',   ipconc)         ) write(6,*) 'INIT_GRID:  Problem setting IPCONC'
  IF( .not. SET_VARIABLE(gd,'ICHAFF',   ichaff)         ) write(6,*) 'INIT_GRID:  Problem setting ICHAFF'
  IF( .not. SET_VARIABLE(gd,'INUCOPT',  inucopt)        ) write(6,*) 'INIT_GRID:  Problem setting INUCOPT'
  IF( .not. SET_VARIABLE(gd,'IPELEC',   ipelec)         ) write(6,*) 'INIT_GRID:  Problem setting IPELEC'
  IF( .not. SET_VARIABLE(gd,'LGTSEED',  iseed)          ) write(6,*) 'INIT_GRID:  Problem setting LGTSEED'
  IF( .not. SET_VARIABLE(gd,'SST',      thg0)           ) write(6,*) 'INIT_GRID:  Problem setting SST'
  IF( .not. SET_VARIABLE(gd,'IENKF',    ienkf)          ) write(6,*) 'INIT_GRID:  Problem setting IENKF'
  IF( .not. SET_VARIABLE(gd,'ISFCPHYS', isfcphys)       ) write(6,*) 'INIT_GRID:  Problem setting ISFCPHYS'
  IF( .not. SET_VARIABLE(gd,'IMURAIN',  imurain)        ) write(6,*) 'INIT_GRID:  Problem setting IMURAIN'

  IF( .not. SET_VARIABLE(gd,'DT',         dt)        ) write(6,*) 'INIT_GRID:  Problem setting DT'
  IF( .not. SET_VARIABLE(gd,'DX',         dx)        ) write(6,*) 'INIT_GRID:  Problem setting DX'
  IF( .not. SET_VARIABLE(gd,'DY',         dy)        ) write(6,*) 'INIT_GRID:  Problem setting DY'
  IF( .not. SET_VARIABLE(gd,'DZ',         dz)        ) write(6,*) 'INIT_GRID:  Problem setting DZ'
  IF( .not. SET_VARIABLE(gd,'XDOMAIN',  xdomain)     ) write(6,*) 'INIT_GRID:  Problem setting XDOMAIN'
  IF( .not. SET_VARIABLE(gd,'YDOMAIN',  ydomain)     ) write(6,*) 'INIT_GRID:  Problem setting YDOMAIN'
  IF( .not. SET_VARIABLE(gd,'ZDOMAIN',  zdomain)     ) write(6,*) 'INIT_GRID:  Problem setting ZDOMAIN'
  IF( .not. SET_VARIABLE(gd,'DX_STRETCH', dx_stretch)) write(6,*) 'INIT_GRID:  Problem setting DX_STRETCH'
  IF( .not. SET_VARIABLE(gd,'DY_STRETCH', dy_stretch)) write(6,*) 'INIT_GRID:  Problem setting DY_STRETCH'
  IF( .not. SET_VARIABLE(gd,'DZ_STRETCH', dz_stretch)) write(6,*) 'INIT_GRID:  Problem setting DZ_STRETCH'
  IF( .not. SET_VARIABLE(gd,'CORIOLIS',   coriol)    ) write(6,*) 'INIT_GRID:  Problem setting CORIOLIS'
  IF( .not. SET_VARIABLE(gd,'CDM',        drag)      ) write(6,*) 'INIT_GRID:  Problem setting CDM'
! IF( .not. SET_VARIABLE(gd,'PSFC',       psfc)      ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
! IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc)      ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
! IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc)      ) write(6,*) 'INIT_GRID:  Problem setting QSFC'
!  IF( .not. SET_VARIABLE(gd,'UGRID',      0.0)       ) write(6,*) 'INIT_GRID:  Problem setting UGRID'
!  IF( .not. SET_VARIABLE(gd,'VGRID',      0.0)       ) write(6,*) 'INIT_GRID:  Problem setting VGRID'
  IF( .not. SET_VARIABLE(gd,'XG_POS',     x_sw_loc)  ) write(6,*) 'INIT_GRID:  Problem setting XG_POS'
  IF( .not. SET_VARIABLE(gd,'YG_POS',     y_sw_loc)  ) write(6,*) 'INIT_GRID:  Problem setting YG_POS'
  IF( .not. SET_VARIABLE(gd,'DX_CORNER',  0.0)       ) write(6,*) 'INIT_GRID:  Problem setting DX_CORNER'
  IF( .not. SET_VARIABLE(gd,'DY_CORNER',  0.0)       ) write(6,*) 'INIT_GRID:  Problem setting DY_CORNER'
  IF( .not. SET_VARIABLE(gd,'RHO_QR',     rho_qr)    ) write(6,*) 'INIT_GRID:  Problem setting RHO_QR'
  IF( .not. SET_VARIABLE(gd,'RHO_QS',     rho_qs)    ) write(6,*) 'INIT_GRID:  Problem setting RHO_QS'
  IF( .not. SET_VARIABLE(gd,'RHO_QH',     rho_qh)    ) write(6,*) 'INIT_GRID:  Problem setting RHO_QH'
  IF( .not. SET_VARIABLE(gd,'RHO_QHL',    rho_qhl)   ) write(6,*) 'INIT_GRID:  Problem setting RHO_QHL'
  IF( .not. SET_VARIABLE(gd,'CNOR',       cnor)      ) write(6,*) 'INIT_GRID:  Problem setting CNOR'
  IF( .not. SET_VARIABLE(gd,'CNOS',       cnos)      ) write(6,*) 'INIT_GRID:  Problem setting CNOS'
  IF( .not. SET_VARIABLE(gd,'CNOH',       cnoh)      ) write(6,*) 'INIT_GRID:  Problem setting CNOH'
  IF( .not. SET_VARIABLE(gd,'CNOHL',      cnohl)     ) write(6,*) 'INIT_GRID:  Problem setting CNOHL'
  IF( .not. SET_VARIABLE(gd,'RHO_QR_10',  rho_qr_i10) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QR_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QS_10',  rho_qs_i10) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QS_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QF_10',  rho_qf_i10) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QF_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QGL_10', rho_qgl_i10)) write(6,*) 'INIT_GRID:  Problem setting RHO_QGL_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QGM_10', rho_qgm_i10)) write(6,*) 'INIT_GRID:  Problem setting RHO_QGM_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QGH_10', rho_qgh_i10)) write(6,*) 'INIT_GRID:  Problem setting RHO_QGH_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QH_10',  rho_qh_i10) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QH_10'
  IF( .not. SET_VARIABLE(gd,'RHO_QHL_10', rho_qhl_i10)) write(6,*) 'INIT_GRID:  Problem setting RHO_QHL_10'
  IF( .not. SET_VARIABLE(gd,'CNOR_10',    cnor_i10)   ) write(6,*) 'INIT_GRID:  Problem setting CNOR_10'
  IF( .not. SET_VARIABLE(gd,'CNOS_10',    cnos_i10)   ) write(6,*) 'INIT_GRID:  Problem setting CNOS_10'
  IF( .not. SET_VARIABLE(gd,'CNOF_10',    cnof_i10)   ) write(6,*) 'INIT_GRID:  Problem setting CNOF_10'
  IF( .not. SET_VARIABLE(gd,'CNOGL_10',   cnogl_i10)  ) write(6,*) 'INIT_GRID:  Problem setting CNOGL_10'
  IF( .not. SET_VARIABLE(gd,'CNOGM_10',   cnogm_i10)  ) write(6,*) 'INIT_GRID:  Problem setting CNOGM_10'
  IF( .not. SET_VARIABLE(gd,'CNOGH_10',   cnogh_i10)  ) write(6,*) 'INIT_GRID:  Problem setting CNOGH_10'
  IF( .not. SET_VARIABLE(gd,'CNOH_10',    cnoh_i10)   ) write(6,*) 'INIT_GRID:  Problem setting CNOH_10'
  IF( .not. SET_VARIABLE(gd,'CNOHL_10',   cnohl_i10)  ) write(6,*) 'INIT_GRID:  Problem setting CNOHL_10'
  ! DTD: added MY variables
  IF( .not. SET_VARIABLE(gd,'NTC_MY',     ntc_my) )    write(6,*) 'INIT_GRID:  Problem setting NTC_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QR_MY',  rho_qr_my) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QR_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QI_MY',  rho_qi_my) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QI_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QS_MY',  rho_qs_my) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QS_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QG_MY',  rho_qg_my) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QG_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QH_MY',  rho_qh_my) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QH_MY'
  IF( .not. SET_VARIABLE(gd,'CNOR_MY',    cnor_my) )   write(6,*) 'INIT_GRID:  Problem setting CNOR_MY'
  IF( .not. SET_VARIABLE(gd,'CNOS_MY',    cnos_my) )   write(6,*) 'INIT_GRID:  Problem setting CNOS_MY'
  IF( .not. SET_VARIABLE(gd,'CNOG_MY',    cnog_my) )   write(6,*) 'INIT_GRID:  Problem setting CNOG_MY'
  IF( .not. SET_VARIABLE(gd,'CNOH_MY',    cnoh_my) )   write(6,*) 'INIT_GRID:  Problem setting CNOH_MY'
  IF( .not. SET_VARIABLE(gd,'ALPHAR_MY',  alphar_my) )   write(6,*) 'INIT_GRID:  Problem setting ALPHAR_MY'
  IF( .not. SET_VARIABLE(gd,'ALPHAI_MY',  alphai_my) )   write(6,*) 'INIT_GRID:  Problem setting ALPHAI_MY'
  IF( .not. SET_VARIABLE(gd,'ALPHAS_MY',  alphas_my) )   write(6,*) 'INIT_GRID:  Problem setting ALPHAS_MY'
  IF( .not. SET_VARIABLE(gd,'ALPHAG_MY',  alphag_my) )   write(6,*) 'INIT_GRID:  Problem setting ALPHAG_MY'
  IF( .not. SET_VARIABLE(gd,'ALPHAH_MY',  alphah_my) )   write(6,*) 'INIT_GRID:  Problem setting ALPHAH_MY'
  ! DTD: end MY variables
  IF( .not. SET_VARIABLE(gd,'RHO_QG_TAK',  rho_qh_tak) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QG_MY'
  IF( .not. SET_VARIABLE(gd,'RHO_QH_TAK',  rho_qhl_tak) ) write(6,*) 'INIT_GRID:  Problem setting RHO_QH_MY'
  IF( .not. SET_VARIABLE(gd,'QCMINCON',   qcmincwrn) ) write(6,*) 'INIT_GRID:  Problem setting QCMINCON'
  IF( .not. SET_VARIABLE(gd,'CWDIAP',     cwdiap)    ) write(6,*) 'INIT_GRID:  Problem setting CWDIAP'
  IF( .not. SET_VARIABLE(gd,'CWDISP',     cwdisp)    ) write(6,*) 'INIT_GRID:  Problem setting CWDISP'
  IF( .not. SET_VARIABLE(gd,'CCN',        ccn)       ) write(6,*) 'INIT_GRID:  Problem setting CCN'
  IF( .not. SET_VARIABLE(gd,'CCNUF',      ccnuf)     ) write(6,*) 'INIT_GRID:  Problem setting CCNUF'
  IF( .not. SET_VARIABLE(gd,'CCNAC',      ccnac)     ) write(6,*) 'INIT_GRID:  Problem setting CCNAC'
  IF( .not. SET_VARIABLE(gd,'LAT',        lat)       ) write(6,*) 'INIT_GRID:  Problem setting LAT'
  IF( .not. SET_VARIABLE(gd,'LON',        lon)       ) write(6,*) 'INIT_GRID:  Problem setting LON'
  IF( .not. SET_VARIABLE(gd,'HGT',        hgt)       ) write(6,*) 'INIT_GRID:  Problem setting HGT'
  IF( .not. SET_VARIABLE(gd,'TRESTART',   60 )       ) write(6,*) 'INIT_GRID:  Problem setting TRESTART'
  IF( .not. SET_VARIABLE(gd,'THISTORY',   30 )       ) write(6,*) 'INIT_GRID:  Problem setting THISTORY'
  IF( .not. SET_VARIABLE(gd,'TVIS5D',     999999)    ) write(6,*) 'INIT_GRID:  Problem setting TVIS5D'
  IF( .not. SET_VARIABLE(gd,'TPRINT',     60 )       ) write(6,*) 'INIT_GRID:  Problem setting TPRINT'
  IF( .not. SET_VARIABLE(gd,'TIME_STOP',  60 )       ) write(6,*) 'INIT_GRID:  Problem setting STOP'
  IF( .not. SET_VARIABLE(gd,'ALPHAH',     alphah )   ) write(6,*) 'INIT_GRID:  Problem setting ALPHAH'
  IF( .not. SET_VARIABLE(gd,'ALPHAHL',    alphahl )  ) write(6,*) 'INIT_GRID:  Problem setting ALPHAHL'
  IF( .not. SET_VARIABLE(gd,'ALPHAI',     alphai )   ) write(6,*) 'INIT_GRID:  Problem setting ALPHAI'
  IF( .not. SET_VARIABLE(gd,'ALPHAS',     alphas )   ) write(6,*) 'INIT_GRID:  Problem setting ALPHAS'
  IF( .not. SET_VARIABLE(gd,'EHSLFO0',    ehslfo0 )  ) write(6,*) 'INIT_GRID:  Problem setting EHSLFO0'
!--------------------------------------------------------------------------
! ERROR CHECK HORIZONTAL GRID STRETCH

!#ifdef MPI
  IF((dx .gt. dx_stretch .and. nxend-20 .lt. n_middle ) .or. ( dy .gt. dy_stretch .and. nyend-20 .lt. n_middle ) ) THEN
 ! (nxend-20 .lt. n_middle .or. nyend-20 .lt. n_middle)) THEN
!#else
!  IF((dx .gt. dx_stretch .or. dy .gt. dy_stretch) .and. (nx-20 .lt. n_middle .or. ny-20 .lt. n_middle)) THEN
!#endif
   write(luno,*) 'INIT_GRID:  INSUFFICIENT POINTS TO CREATE STRETCHED GRID'
   write(luno,*) 'INIT_GRID:  STOPPING RUN!!!!'
   STOP

  ENDIF 


!--------------------------------------------------------------------------
! Generate east-west grid
! Retrieve the pointers to the x-grid variables

  CALL GET_VARIABLE(gd, 'XC',   gc)
  CALL GET_VARIABLE(gd, 'XE',   ge)
  CALL GET_VARIABLE(gd, 'DXC',  dgc)
  CALL GET_VARIABLE(gd, 'DXE',  dge)

#ifdef OPENDX
  CALL GET_VARIABLE(gd, 'XCDX',   xcdx)
  CALL GET_VARIABLE(gd, 'XEDX',   xedx)
#endif

  gcindex  = GET_VARIABLE_INDEX(gd, 'XC')
  ngs = gd%var(gcindex)%ng

#ifdef MPI
  ixb = -ng+1
  ixe = itile+ng
  if (ixbeg .eq. nxbeg) ixb = 1
  if (ixend .eq. nxend) ixe = ixend-ixbeg+1

#ifdef OPENDX
  DO dum = 1,3
  DO i = ixb,ixe
   xcdx (dum,i) = 0.0
   xedx (dum,i) = 0.0
  ENDDO ; ENDDO
#endif
#else
#ifdef OPENDX
  xcdx (1:3,1:nx) = 0.0
  xedx (1:3,1:nx) = 0.0
#endif 
#endif

  IF( dx .le. dx_stretch*1.001 ) THEN

#ifdef MPI
   ixb = -ng+1
   ixe = itile+ng
   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg+1

   DO i = ixb,ixe
    gc(i) = (float(ixbeg+i-1) - 0.5) * dx    ! SGX
    ge(i) = (float(ixbeg+i-1) - 1.0) * dx    ! UGX
   END DO
#else
   DO i = 1,nx
    gc(i) = (float(i) - 0.5) * dx    ! SGX
    ge(i) = (float(i) - 1.0) * dx    ! UGX
   ENDDO
#endif

   dgc(:) = 1.0 / dx
   dge(:) = 1.0 / dx

#ifdef MPI
   write(luno,*) 'INIT_GRID:  X-SUBDOMAIN  LEN  = ', ge(nx)
#else
   write(luno,*) 'INIT_GRID:  X-DOMAIN  LEN  = ', ge(nx)
#endif

  ELSE

#ifdef MPI
   ALLOCATE(gctemp(-ngs+1:nxend+ngs))
   ALLOCATE(getemp(-ngs+1:nxend+ngs))
   CALL HGRID(dx, dx_stretch, nxend, n_middle, gctemp, getemp, ngs, dxmax)
   ixb = -ng+1
   ixe = itile+ng

!   if (ixbeg .eq. nxbeg) ixb = 1
!   if (ixend .eq. nxend) ixe = ixend-ixbeg+1

   DO i = ixb,ixe
     gc(i)=gctemp(i+ixbeg-1)
     ge(i)=getemp(i+ixbeg-1)
!     write(luno,*) 'rank,i,gc,ge = ',my_rank,i,gc(i),ge(i)
   ENDDO
!   write(luno,*) ixb,ixe,itile
   DO i=1,nxend
!   write(luno,*) 'i,getempX:  ',i,getemp(i)
   ENDDO
!   write(luno,*) getemp
   write(luno,*)'gcx:',gctemp(1),gc(1),':',gc(itile-1),  gctemp(nxend-1)
   write(luno,*)'gex:',getemp(1),ge(1),':',ge(itile),    getemp(nxend)
   DEALLOCATE(gctemp)
   DEALLOCATE(getemp)
#else
   CALL HGRID(dx, dx_stretch, nx, n_middle, gc, ge, ngs, dxmax)
   DO i = 1,nx
!     write(luno,*) 'i,gc,ge = ',i,gc(i),ge(i)
   ENDDO
#endif

! Compute 1 / DX factors for stretched grid

#ifdef MPI
   if (ixbeg .eq. nxbeg)    &
#endif
   dgc(1) = 0.5 / (gc(1) - ge(1))

#ifdef MPI
   ixb = -ng+1
   ixe = itile+ng-1
   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg+1

   DO i = ixb,ixe
                             dgc(i) = 1.0 / ( ge(i+1) - ge(i)   ) ! DXC(i)
    IF( (ixbeg+i-1) .ne. 1 ) dge(i) = 1.0 / ( gc(i)   - gc(i-1) ) ! DXE(i)
    IF( (ixbeg+i-1) .eq. 1 ) dge(i) = 1.0 / ( 2.0*gc(i) ) ! DXE(i)
!     write(luno,*) 'i,dgc,dge = ',i,dgc(i),dge(i)
   ENDDO
#else

   DO i = 1,nx-1
                   dgc(i) = 1.0 / ( ge(i+1) - ge(i)   ) ! DXC(i)
    IF( i .ne. 1 ) dge(i) = 1.0 / ( gc(i)   - gc(i-1) ) ! DXE(i)
    IF( i .eq. 1 ) dge(i) = 1.0 / ( 2.0*gc(i) ) ! DXE(i)
   ENDDO
#endif

   if (ixend .eq. nxend ) dge(nx) = 0.5 / (ge(nx) - gc(nx-1))

   if (ixend .eq. nxend) THEN
     dge(nx) = 0.5 / (ge(nx) - gc(nx-1))
     dge(nx+1) = dge(nx)
     dgc(nx:nx+1) = dgc(nx-1)
    endif

  ENDIF

#ifdef MPI
   ixb = -ng+1
   ixe = itile+ng
   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg
   
#ifdef OPENDX
   DO i = ixb,ixe
    xcdx(3,i) = 0.001*gc(i)
   ENDDO
#endif

   ixb = -ng+1
   ixe = itile+ng
   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg+1
   
#ifdef OPENDX
   DO i = ixb,ixe
    xedx(3,i)   = 0.001*ge(i)
   ENDDO
#endif
#else
#ifdef OPENDX
    xcdx(3,1:nx-1) = 0.001*gc(1:nx-1)
    xedx(3,1:nx)   = 0.001*ge(1:nx)
#endif
#endif

#ifdef OPENDX
    xc1 = xcdx(3,1)
    xe1 = xedx(3,1)
#endif  

!--------------------------------------------------------------------------
! Generate south-north grid
! Retrieve the pointers to the x-grid variables

  CALL GET_VARIABLE(gd,'YC',   gc)
  CALL GET_VARIABLE(gd,'YE',   ge)
  CALL GET_VARIABLE(gd,'DYC',  dgc)
  CALL GET_VARIABLE(gd,'DYE',  dge)

#ifdef OPENDX
  CALL GET_VARIABLE(gd, 'YCDX',   ycdx)
  CALL GET_VARIABLE(gd, 'YEDX',   yedx)
#endif

  gcindex  = GET_VARIABLE_INDEX(gd, 'YC')
  ngs = gd%var(gcindex)%ng

#ifdef MPI
  jyb = -ng+1
  jye = jtile+ng
  if (jybeg .eq. nybeg) jyb = 1
  if (jyend .eq. nyend) jye = jyend-jybeg+1
  
#ifdef OPENDX
  DO dum=1,3
  DO j=jyb,jye
   ycdx (dum,j) = 0.0
   yedx (dum,j) = 0.0  
  END DO
  END DO
#endif

#else
#ifdef OPENDX
  ycdx (1:3,1:ny) = 0.0
  yedx (1:3,1:ny) = 0.0
#endif
#endif

  IF( dy .le. dy_stretch*1.001 ) THEN

#ifdef MPI
   jyb = -ng+1
   jye = jtile+ng
   if (jybeg .eq. nybeg) jyb = 1
   if (jyend .eq. nyend) jye = jyend-jybeg+1

   DO j = jyb,jye
    gc(j) = (float(jybeg+j-1) - 0.5) * dy    ! SGX
    ge(j) = (float(jybeg+j-1) - 1.0) * dy    ! UGX
   ENDDO
#else   
   DO j = 1,ny
    gc(j) = (float(j) - 0.5) * dy    ! SGX
    ge(j) = (float(j) - 1.0) * dy    ! UGX
   ENDDO
#endif

   dgc(:) = 1.0 / dy
   dge(:) = 1.0 / dy

#ifdef MPI
   write(luno,*) 'INIT_GRID:  Y-SUBDOMAIN  LEN  = ', ge(ny)
#else
   write(luno,*) 'INIT_GRID:  Y-DOMAIN  LEN  = ', ge(ny)
#endif

  ELSE

#ifdef MPI
   ALLOCATE(gctemp(-ngs+1:nyend+ngs))
   ALLOCATE(getemp(-ngs+1:nyend+ngs))
   CALL HGRID(dy, dy_stretch, nyend, n_middle, gctemp, getemp, ngs ,dymax)
   jyb = -ngs+1
   jye = jtile+ngs
!   if (jybeg .eq. nybeg) jyb = 1
!   if (jyend .eq. nyend) jye = jyend-jybeg+1
   DO j = jyb, jye
     ge(j)=getemp(j+jybeg-1)
     gc(j)=gctemp(j+jybeg-1)
   ENDDO
!   write(luno,*) 'getempY:  ', getemp
!   write(luno,*)jyb,jye
   write(luno,*)'gcy:',gctemp(1),gc(1),':',gc(jtile-1),gctemp(nyend-1)
   write(luno,*)'gey:',getemp(1),ge(1),':',ge(jtile),  getemp(nyend)
   DEALLOCATE(gctemp)
   DEALLOCATE(getemp)
#else
   CALL HGRID(dy, dy_stretch, ny, n_middle, gc, ge, ngs ,dymax)
#endif

! Compute 1 / DY factors for stretched grid
#ifdef MPI
   if (jybeg .eq. nybeg)   &
#endif
   dgc(1) = 0.5 / (gc(1) - ge(1))

#ifdef MPI
   jyb = -ng+2
   jye = jtile+ng-1
   if (jybeg .eq. nybeg) jyb = 1
   if (jyend .eq. nyend) jye = jyend-jybeg+1

   DO j = jyb,jye
                             dgc(j) = 1.0 / ( ge(j+1) - ge(j)   ) ! DYC(j)
    IF( (jybeg+j-1) .ne. 1 ) dge(j) = 1.0 / ( gc(j)   - gc(j-1) ) ! DYE(j)
    IF( (jybeg+j-1) .eq. 1 ) dge(j) = 1.0 / ( 2.0*gc(j) ) ! DYE(j)
   ENDDO
#else

   DO j = 1,ny-1
                   dgc(j) = 1.0 / ( ge(j+1) - ge(j)   ) ! DYC(j)
    IF( j .ne. 1 ) dge(j) = 1.0 / ( gc(j)   - gc(j-1) ) ! DYE(j)
    IF( j .eq. 1 ) dge(j) = 1.0 / ( 2.0*gc(j) ) ! DYE(j)
   ENDDO
#endif

   if (jyend .eq. nyend) THEN
     dge(ny) = 0.5 / (ge(ny) - gc(ny-1))
     dge(ny+1) = dge(ny)
     dgc(ny:ny+1) = dgc(ny-1)
    endif

  ENDIF

#ifdef MPI
   jyb = -ng+1
   jye = jtile+ng
   if (jybeg .eq. nybeg) jyb = 1
   if (jyend .eq. nyend) jye = jyend-jybeg
  
#ifdef OPENDX
   DO j = jyb,jye
     ycdx(2,j) = 0.001*gc(j)
   END DO
#endif
   jyb = -ng+1
   jye = jtile+ng
   if (jybeg .eq. nybeg) jyb = 1
   if (jyend .eq. nyend) jye = jyend-jybeg+1
  
#ifdef OPENDX
   DO j = jyb,jye
     yedx(2,j)   = 0.001*ge(j)
   END DO   
#endif   
#else
#ifdef OPENDX
  ycdx(2,1:ny-1) = 0.001*gc(1:ny-1)
  yedx(2,1:ny)   = 0.001*ge(1:ny)
#endif
#endif


#ifdef OPENDX
  yc1 = ycdx(2,1)
  ye1 = yedx(2,1)
#endif
!--------------------------------------------------------------------------
! Generate vertical grid
! Retrieve the pointers to the z-grid variables

  CALL GET_VARIABLE(gd,'ZC',   gc)
  CALL GET_VARIABLE(gd,'ZE',   ge)
  CALL GET_VARIABLE(gd,'DZC',  dgc)
  CALL GET_VARIABLE(gd,'DZE',  dge)

#ifdef OPENDX
  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
#endif

  gcindex  = GET_VARIABLE_INDEX(gd, 'ZC')
  ngs = gd%var(gcindex)%ng

#ifdef MPI
  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .eq. nzbeg) kzb = 1
  if (kzend .eq. nzend) kze = kzend-kzbeg+1

#ifdef OPENDX
  DO dum=1,3  
  DO k=kzb,kze
   zcdx (dum,k) = 0.0
   zedx (dum,k) = 0.0
  ENDDO
  ENDDO
#endif

#else  
#ifdef OPENDX
  zcdx (1:3,1:nz) = 0.0
  zedx (1:3,1:nz) = 0.0
#endif
#endif
   
   IF ( .not. allocated(z1dinit) ) THEN
     allocate( z1dinit(-ngs+1:nzend+ngs,4) ) ! zc,ze,dzc,dze
     z1dinit(:,:) = 0.0d0
   ENDIF
   
   CALL ZGRID(dz, dz_stretch, nzend, nbndlyr, z1dinit(-ngs+1,1), z1dinit(-ngs+1,2), ngs)
   
   DO k = 1,nzend-1
     write(luno,*) 'k,z1dinit = ',k,z1dinit(k,1),z1dinit(k,2)
   ENDDO
   
   kzb = -ngs+1
   kze = ktile+ngs
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg+1
   DO k = kzb, kze
     gc(k)=z1dinit(k+kzbeg-1,1)
     ge(k)=z1dinit(k+kzbeg-1,2)
   ENDDO
!   write(luno,*) 'getempY:  ', getemp
!   write(luno,*)jyb,jye
   write(luno,*)'gcz:',z1dinit(1,1),gc(1),':',gc(ktile-1),z1dinit(nzend-1,1)
   write(luno,*)'gez:',z1dinit(1,2),ge(1),':',ge(ktile),  z1dinit(nzend,2)

  z1dinit(1,3) = 0.5 / (z1dinit(1,1) - z1dinit(1,2))
  
  DO k = 1,nzend-1
                            z1dinit(k,3) = 1.0 / ( z1dinit(k+1,2) - z1dinit(k,2)   ) ! DYC(k)
   IF ( k .ne. 1 ) z1dinit(k,4) = 1.0 / ( z1dinit(k,1)   - z1dinit(k-1,1) ) ! DYE(k)

  ENDDO
  
   z1dinit(nzend,4) = 0.5 / (z1dinit(nzend,2) - z1dinit(nzend-1,1))
   z1dinit(nzend,3) = z1dinit(nzend-1,3)
   z1dinit(1,4)  = z1dinit(1,3)

   kzb = -ngs+1
   kze = ktile+ngs
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg+1
   DO k = kzb, kze
     dgc(k)=z1dinit(k+kzbeg-1,3)
     dge(k)=z1dinit(k+kzbeg-1,4)
   ENDDO
  
! Compute 1 / DZ factors for stretched grid

!#ifdef MPI
!  if (kzbeg .eq. nzbeg)    &
!#endif
!  dgc(1) = 0.5 / (gc(1) - ge(1))
!
!#ifdef MPI
!  kzb = -ng+1
!  kze = ktile+ng
!  if (kzbeg .eq. nzbeg) kzb = 1
!  if (kzend .eq. nzend) kze = kzend-kzbeg
!  
!  DO k = kzb,kze
!                            dgc(k) = 1.0 / ( ge(k+1) - ge(k)   ) ! DYC(k)
!   IF( (kzbeg+k-1) .ne. 1 ) dge(k) = 1.0 / ( gc(k)   - gc(k-1) ) ! DYE(k)
!  ENDDO
!#else
!
!  DO k = 1,nz-1
!                  dgc(k) = 1.0 / ( ge(k+1) - ge(k)   ) ! DYC(k)
!   IF( k .ne. 1 ) dge(k) = 1.0 / ( gc(k)   - gc(k-1) ) ! DYE(k)
!  ENDDO
!#endif
!
!#ifdef MPI
!  if (kzend .eq. nzend)   dge(nzend) = 0.5 / (ge(nzend) - gc(nzend-1))
!#else
!  dge(nz) = 0.5 / (ge(nz) - gc(nz-1))
!#endif
!
!  dge(1)  = dgc(1)


#ifdef MPI
  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .eq. nzbeg) kzb = 1
  if (kzend .eq. nzend) kze = kzend-kzbeg

#ifdef OPENDX
  DO k = kzb,kze
   zcdx(1,k) = 0.001*gc(k)
  ENDDO
#endif

  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .eq. nzbeg) kzb = 1
  if (kzend .eq. nzend) kze = kzend-kzbeg+1

#ifdef OPENDX
  DO k = kzb,kze  
   zedx(1,k)   = 0.001*ge(k)
  ENDDO
#endif
#else  
#ifdef OPENDX
  zcdx(1,1:nz-1) = 0.001*gc(1:nz-1)
  zedx(1,1:nz)   = 0.001*ge(1:nz)
#endif
#endif

#ifdef OPENDX
  zc1 = zcdx(1,1)
  ze1 = zedx(1,1)
#endif
!  xcdx(2,1:nx) = yc1
!  xedx(2,1:nx) = yc1
!  xcdx(1,1:nx) = zc1
!  xedx(1,1:nx) = zc1
!
!  ycdx(3,1:ny) = xc1
!  yedx(3,1:ny) = xc1
!  ycdx(1,1:ny) = zc1
!  yedx(1,1:ny) = zc1
!
!  zcdx(3,1:nz) = xc1
!  zedx(3,1:nz) = xc1
!  zcdx(2,1:nz) = yc1
!  zedx(2,1:nz) = yc1

  luno = FILE_CLOSE('INIT_GRID END')

  RETURN

CONTAINS

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE HGRID  <<<<<<<<<<<<<<<<<<<<<<<<<<< !
!
!-------------------------------------------------------------------------------
! Creates the horizontally stretched grid mapping for the model
!
! N_INNER is the number of points in the middle of the grid where dx_middle 
! is constant
!-------------------------------------------------------------------------------

 SUBROUTINE HGRID(dh, dh_middle, n_total, n_middle, gxc, gxe, ng, dxmax)

   USE GRID_MODULE
!   USE PARAM_MODULE
   USE COMMASMPI_MODULE
   implicit none
   
!   include 'param.h'

   real dh, dh_middle
   integer n_total, n_middle, ng
   real gxc(-ng+1:n_total+ng)
   real gxe(-ng+1:n_total+ng)
   real dxmax
!   parameter( dxmax = 10000. )

   integer n1, n
   double precision stretch, xx, xmid, fmid, hlen

   real zheight
!   external zheight

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

! Use Newton interation to find grid coefficients

    n1      = (n_total) / 2

!    hlen    = dh * n1 
    hlen    = dh * ( n_total - 1 ) * 0.5 
    xx      = 1.0d0
    stretch = 1.0d0

!   print *, 'HGRID: N1   = ', n1
!   print *, 'HGRID: HLEN = ', hlen

    IF( dh .gt. dh_middle ) THEN

     DO n = 1,100

      IF( abs(xx) .gt. 1.0e-10 ) THEN
       xx   = xx * 0.5
       xmid = stretch + xx
       IF ( n_total /= (n_total/2)*2 ) THEN
       fmid = ZHEIGHT(dh_middle,xmid,n1,dxmax,n_middle/2,1.e30,1.0,dxmax) - hlen 
       ELSE
       fmid = ZHEIGHT(dh_middle,xmid,n1,dxmax,n_middle/2,1.e30,1.0,dxmax) - hlen - dh_middle/2
       ENDIF
       IF( fmid .le. 0.0 ) stretch = xmid
      ENDIF

     ENDDO

    ENDIF

    write(luno,*)
    IF( stretch .gt. 1.1 ) THEN
     write(luno,*) 'STRETCH FAC TOO BIG! - NUMERICAL ERRORS WILL BE LARGE'
     write(luno,*) 'STRETCH FAC  = ',stretch
     write(luno,*) 'INCREASE NUMBER OF HORIZONTAL POINTS in 3d.run'
     IF ( my_rank == 0 ) THEN
     write(0,*) 'STRETCH FAC TOO BIG! - NUMERICAL ERRORS WILL BE LARGE'
     write(0,*) 'STRETCH FAC  = ',stretch
     write(0,*) 'INCREASE NUMBER OF HORIZONTAL POINTS in 3d.run'
     ENDIF
#ifdef MPI
      CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
     CALL commasmpi_abort()
    ELSE
     write(luno,*) 'HGRID:  STRETCH FAC  = ',stretch
     write(luno,*) 'HGRID:  DOMAIN  LEN  = ',2*ZHEIGHT(dh_middle,stretch,n1,dxmax,n_middle/2,1.e30,1.0,dxmax)
    ENDIF

!----------------------------------------------------------------------
! Create grid


!    write(lune,*) 'HGRID: my_rank, n_total,n1,hlen = ',my_rank,n_total,n1,hlen

    gxe(1)    = 0.0
!    gxe(n1+1) = hlen

    IF ( n_total /= (n_total/2)*2 ) THEN
    gxe(n1+1) = hlen
    
    DO n = 1,n1
!     write(lune,*) 'HGRID: myrank, ',my_rank,n1+n+1,n1-n+1
!     IF ( n1+n+1 .le. n_total ) THEN
     gxe(n1+n+1) = hlen + ZHEIGHT(dh_middle,stretch,n,dxmax,n_middle/2,1.e30,1.0,dxmax)
!     ENDIF
     gxe(n1-n+1) = hlen - ZHEIGHT(dh_middle,stretch,n,dxmax,n_middle/2,1.e30,1.0,dxmax)
    ENDDO

    ELSE

    DO n = 1,n1
     gxe(n1+n) = hlen - dh_middle/2 + ZHEIGHT(dh_middle,stretch,n,dxmax,n_middle/2,1.e30,1.0,dxmax)
    ENDDO

    DO n = 1,n1
     gxe(n1-n+1) = hlen + dh_middle/2 - ZHEIGHT(dh_middle,stretch,n,dxmax,n_middle/2,1.e30,1.0,dxmax)
    ENDDO
 
    ENDIF

    DO n = 1,n_total-1
      gxc(n) = 0.5*(gxe(n) + gxe(n+1))
!     write(lune,*) 'HGRID: myrank, n,gxc,gxe, ',my_rank,n,gxc(n),gxe(n)
    ENDDO

    gxc(n_total) = 2*gxc(n_total-1) - gxc(n_total-2) 

  END SUBROUTINE HGRID

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE ZGRID  <<<<<<<<<<<<<<<<<<<<<<<<<<< !
!
!-------------------------------------------------------------------------------
! Creates the vertical grid using a geometric stretch 
!
! Option - set nbndlyr > 0 input deck.  A suggested value is ~ nz/4 to create a 
! layer having constant resolution near the surface.
!
!-------------------------------------------------------------------------------

 SUBROUTINE ZGRID(dz, dz_stretch, nz, nbndlyr, gzc, gze, ng)

!   USE PARAM_MODULE
#ifdef MPI
   USE COMMASMPI_MODULE
#endif
   implicit none

!   include 'param.h'

   real dz, dz_stretch
   integer nbndlyr, nz, ng
   double precision gzc(-ng+1:nz+ng), gze(-ng+1:nz+ng)
#ifdef MPI
   real gzct(-ng+1:nz+ng), gzet(-ng+1:nz+ng)
#endif

!   real dzmax
!   parameter( dzmax = 700. )
   integer n, k
   double precision stretch, zx, xmid, fmid, ztop

   real zheight
!   external zheight

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   integer :: ktmp
#endif

! Use Newton interation to find grid coefficients

   ztop    = dz * (nz-1)
   zx      = 1.0d0
   stretch = 1.0d0

   IF( dz .gt. dz_stretch ) THEN

    DO n = 1,50

     IF( abs(zx) .gt. 1.0e-12 ) THEN
      zx   = zx * 0.5
      xmid = stretch + zx
      fmid = ZHEIGHT(dz_stretch,xmid,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop) - ztop
!      write(6,*) 'Stretch: ',n,zx,xmid,fmid,fmid + ztop
      IF( fmid .le. 0.0 ) stretch = xmid
      IF ( fmid .eq. 0.0d0 ) EXIT
     ENDIF

    ENDDO

   ENDIF

   write(luno,*)
   IF( stretch .gt. max_stretch ) THEN
    write(6,*) 'STRETCH FAC TOO BIG! - NUMERICAL ERRORS WILL BE LARGE'
    write(6,*) 'STRETCH FAC  = ',stretch
    write(6,*) 'INCREASE NZ in 3d.run'
    STOP
   ELSE
    write(luno,*) 'ZGRID:  STRETCH FAC  = ',stretch
    write(luno,*) 'ZGRID:  DOMAIN  HGT  = ',ZHEIGHT(dz_stretch,stretch,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
   ENDIF

!----------------------------------------------------------------------
! Create grid

!#ifdef MPI
!   gzct(:) = 0.0
!   gzet(:) = 0.0
!
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!
!   DO k = kzb, kze
!    ktmp = kzbeg+k-1
!    gzet(k+1) = ZHEIGHT(dz_stretch,stretch,ktmp,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
!    gzct(k)   = 0.5 * ( gzet(k) + gzet(k+1) )
!   ENDDO
!
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!   do k = kzb, kze
!    gze(k) = gzet(k) !gzet(kzbeg + k-1)
!    gzc(k) = gzct(k) !gzct(kzbeg + k-1)
!   enddo
!
!   if (kzend.eq.nzend) gzc(nz) = 2*gzc(nz-1) - gzc(nz-2) 
!
!#else
   gze(1) = 0.0
   DO k = 1,nz-1
    gze(k+1) = ZHEIGHT(dz_stretch,stretch,k,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
    gzc(k)   = 0.5 * ( gze(k) + gze(k+1) )
   ENDDO

   gzc(nz) = 2.*gzc(nz-1) - gzc(nz-2) 

!#endif

  END SUBROUTINE ZGRID

END SUBROUTINE INIT_GRID
 
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_PERT_BUBBLE  <<<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain with thermal perturubations
!
!-----------------------------------------------------------------------------
 SUBROUTINE INIT_PERT_BBLE(gd) !, tbble, xrad, yrad, zrad, xcntr, ycntr, zcntr, bbletype, nbble, bblsp)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE INIT_MODULE
  USE INDEX_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif
  
! Passed variables

  TYPE(GRID) :: gd

!  integer nbble, bbletype, ibbleseed
!  integer bblsp  (nbble)
!  real    tbble  (nbble)
!  real    xrad   (nbble)
!  real    yrad   (nbble)
!  real    zrad   (nbble)
!  real    xcntr  (nbble)
!  real    ycntr  (nbble)
!  real    zcntr  (nbble)

! Local vars

  integer i, j, k, i1, j1
  integer n

  real, allocatable :: ranarray(:,:),ranarraysmth(:,:)
  real, allocatable :: thpert(:,:,:)
  integer           :: iranseed
  integer           :: ibsrc
  real              :: dfac, dth, rndnum, radius
  character(LEN=20) :: strtmp,strtmp1
  double precision  :: dthtot, dthrantot, dthrantotv, dthran
  double precision  :: dthrantotn, dthrantotp
  double precision  :: dthrantotvn, dthrantotvp
  real              :: dv
  double precision  :: fn, fp
  real              :: vr, g1
  
  real              :: pres, qvs, rh, qvs1
  real              :: fac
  
  integer, pointer :: nx, ny, nz
  real             :: ugrid1, vgrid1
!  integer, pointer :: nxend, nyend, nzend

  real, allocatable :: pztmp(:)
  real              :: tv0, tv1

  TYPE(ATTRIBUTE), pointer   :: attr, microphys
  
  real, pointer :: xc(:)
  real, pointer :: yc(:)
  real, pointer :: zc(:)
  real, pointer :: dxc(:)
  real, pointer :: dyc(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: tz(:)
  real, pointer :: pz(:)
  real, pointer :: qz(:)
  real, pointer :: th(:,:,:)
  real, pointer :: qv3(:,:,:)
  real, pointer :: cn(:,:,:)
  real, pointer :: cn2(:,:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: qr(:,:,:)
  real, pointer :: zr(:,:,:)
  real, pointer :: pi3(:,:,:)
  real, pointer :: p3(:,:,:)
  real, pointer :: tv3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: r3(:,:,:)
  real, pointer :: hwind3(:,:,:)
  real, pointer :: ranarr(:,:)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
#ifdef MPI

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)


#endif

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)
!  CALL GET_ATTRIBUTE(gd, 'NXEND', nxend)
!  CALL GET_ATTRIBUTE(gd, 'NYEND', nyend)
!  CALL GET_ATTRIBUTE(gd, 'NZEND', nzend)

  CALL GET_VARIABLE(gd,'XC',xc)
  CALL GET_VARIABLE(gd,'YC',yc)
  CALL GET_VARIABLE(gd,'ZC',zc)
  CALL GET_VARIABLE(gd,'DXC',dxc)
  CALL GET_VARIABLE(gd,'DYC',dyc)
  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'TH',th)
  CALL GET_VARIABLE(gd,'QV',qv3)
  CALL GET_VARIABLE(gd,'U',u3)
  CALL GET_VARIABLE(gd,'V',v3)
  CALL GET_VARIABLE(gd,'UGRID',ugrid1)
  CALL GET_VARIABLE(gd,'VGRID',vgrid1)
  CALL GET_VARIABLE(gd,'RHO',r3)
  CALL GET_VARIABLE(gd,'RANARRAY2D',ranarr)

  CALL GET_ATTRIBUTE (gd, 'MICROPHYS', microphys)
  
#ifdef MPI
  allocate ( ranarray(-ng+1:nxend+ng,-ng+1:nyend+ng) )
  allocate ( ranarraysmth(-ng+1:nxend+ng,-ng+1:nyend+ng) )
#else
  allocate ( ranarray(-ng+1:nx+ng,-ng+1:ny+ng) )
  allocate ( ranarraysmth(-ng+1:nx+ng,-ng+1:ny+ng) )
#endif
  allocate ( thpert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  
  iranseed = -Abs(ibbleseed)
  dfac = 0.0
  ranarray(:,:) = 0.0
  ranarraysmth(:,:) = 0.0
  thpert(:,:,:) = 0.0

  IF ( bbletype .ge. 2 ) THEN
     IF ( my_rank == 0 .and. bbletype .ge. 2 ) THEN
      write(luno,*)
      write(luno,*) 'Initializing Random Perturbations in initial bubble(s)'
     ENDIF

    IF ( bbletype .gt. 9 ) THEN  ! initial hack for creating a different seed -- use ibbleseed now instead
     iranseed = -bbletype
     write(strtmp,*) bbletype

#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .gt. len(strtmp)) ixe = len(strtmp)-ixbeg+1
     
     do i = ixb,ixe
#else
     DO i = 1,len(strtmp)
#endif
       IF ( strtmp(i:i) .ne. ' ' ) THEN
       write(strtmp1,'(a)') strtmp(i:i)
       read(strtmp1,'(i1)') bbletype
       write(strtmp1,'(a)') strtmp(i+1:len(strtmp))
       read(strtmp1,*) iranseed
       iranseed = -iranseed
       exit
       ENDIF
     ENDDO
     IF ( my_rank == 0 ) THEN
       write(luno,*) 'iranseed,bbletype,strtmp,strtmp1,dthranpert = ',iranseed,bbletype,strtmp,strtmp1,dthranpert
     ENDIF
    ENDIF

    dfac = 1./3.
    
#ifdef MPI
!    ixb = -ng+1
!    ixe = itile+ng
!    if (ixbeg .eq. nxbeg) ixb = 3
!    if (ixend .eq. nxend) ixe = ixend-ixbeg-2
!    
!    jyb = -ng+1
!    jye = jtile+ng
!    if (jybeg .eq. nybeg) jyb = 3
!    if (jyend .eq. nyend) jye = jyend-jybeg-2
!    
    DO j = 3,nyend-2
     DO i = 3,nxend-2
#else    
    DO j = 3,ny-3
     DO i = 3,nx-3
#endif
       ranarray(i,j) = 2.0*(rndnum(iranseed) - 0.5)  ! values of -1 to 1
     ENDDO
    ENDDO
    
    ixb = -ng+1
    ixe = itile+ng
    ixb = 1
    ixe = ixend-ixbeg+1
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = -ng+1
    jye = jtile+ng
    jyb = 1
    jye = jyend-jybeg+1
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg

     DO j = jyb, jye
      DO i = ixb, ixe
        ranarr(i,j) = ranarray(ixbeg-1+i,jybeg-1+j)
      ENDDO
     ENDDO

#ifdef MPI
    DO j = 3,nyend-2,4
     DO i = 3,nxend-2,4
#else    
    DO j = 3,ny-3,4
     DO i = 3,nx-3,4
#endif
       DO  j1 = 0,3
       DO  i1 = 0,3
         ranarraysmth(i+i1,j+j1) = ranarray(i,j)
       ENDDO
       ENDDO
     ENDDO
    ENDDO
  ENDIF ! ( bbletype .ge. 2 )

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_PERT_BBLE BEGIN')

!-----------------------------------------------------------------------
! PRINT OUT SOME BASIC STUFF

  write(luno,*)
  write(luno,*) 'SUBROUTINE SET3D1:  BUBBLE INITIALIZATION'
  write(luno,*)
  write(luno,*)

  IF( bbletype .lt. 2 ) THEN
    write(luno,*) 'Initializing bubbles with simple thermal'
    write(luno,*)
  ENDIF
  IF( bbletype .eq. 2 ) THEN
    write(luno,*) 'Initializing bubbles with 0.1 K random perturbations'
    write(luno,*)
  ENDIF
  IF( bbletype .eq. 3 ) THEN
    write(luno,*) 'Initializing bubbles using Kogan random perturbations'
    write(luno,*)
  ENDIF
    CALL GET_VARIABLE(gd,'PI',pi3)
    
  DO n = 1,nbble
    write(luno,*) 'BUBBLE #      ', n
    write(luno,*) 'BUBBLE MAG  = ', tbble(n)
    write(luno,*) 'X-CENTER    = ', xcntr(n)
    write(luno,*) 'Y-CENTER    = ', ycntr(n)
    write(luno,*) 'Z-CENTER    = ', zcntr(n)
    write(luno,*) 'X-RADIUS    = ', xrad(n)
    write(luno,*) 'Y-RADIUS    = ', yrad(n)
    write(luno,*) 'Z-RADIUS    = ', zrad(n)
    write(luno,*)


    IF ( bblsp(n) .eq. 3 .or. bblsp(n) == 9 ) THEN
!      CALL GET_VARIABLE(gd,'U',u3)
!      CALL GET_VARIABLE(gd,'V',v3)
    ENDIF
    IF ( bblsp(n) .eq. 5 ) THEN
      IF ( microphys%str(1:4) == 'DRYC' ) THEN
        print*,'get qc array'
        CALL GET_VARIABLE(gd,'QC',cn)
      ELSEIF ( microphys%str(1:4) == 'DRYT' ) THEN
        print*,'get tracer array'
        CALL GET_VARIABLE(gd,'TRACER1',cn)
      ELSE
        print*,'get cccn array'
        CALL GET_VARIABLE(gd,'CCCN',cn)
      ENDIF
    ENDIF
    IF ( bblsp(n) .eq. 6 ) THEN
      print*,'get rain q and n, check for crw'
!      i = -1 ! GET_VARIABLE_INDEX(gd, 'CRW',nofail=1)
      IF ( lnr > 1 ) THEN
      CALL GET_VARIABLE(gd,'CRW',cn)
      ENDIF
      IF ( lzr > 1 ) THEN
      CALL GET_VARIABLE(gd,'ZRW',zr)
      ENDIF
      CALL GET_VARIABLE(gd,'QR' ,qr)
    ENDIF
    IF ( bblsp(n) .eq. 7  ) THEN
      IF ( lnchaff .gt. 1 ) THEN
      print*,'get chaff concentration array'
      CALL GET_VARIABLE(gd,'CCHAFF',cn)
      ELSE
       write(0,*) 'CHAFF ARRAY not here!! STOP!'
       CALL commasmpi_abort()
      ENDIF
    ENDIF
    IF ( bblsp(n) .eq. 8 ) THEN
      print*,'get graupel q and n'

      IF ( microphys%str(1:3) == 'TAK' ) THEN
      
       write(0,*) 'get CHL30'
       CALL GET_VARIABLE(gd,'CHL30',cn)
       CALL GET_VARIABLE(gd,'CHL40',cn2)
       write(0,*) 'got CHL30'

      ELSE

!      i = -1 ! GET_VARIABLE_INDEX(gd, 'CRW',nofail=1)
      IF ( lnh > 1 ) THEN
      CALL GET_VARIABLE(gd,'CHW',cn)
      ENDIF
      IF ( lzh > 1 ) THEN
      CALL GET_VARIABLE(gd,'ZHW',zr)
      ENDIF
      CALL GET_VARIABLE(gd,'QH' ,qr)
      
      ENDIF
      
    ENDIF

    IF ( bblsp(n) .eq. 20  ) THEN
      IF ( lash1 .gt. 1 ) THEN
      print*,'get ash mass array'
      CALL GET_VARIABLE(gd,'QASH1',cn)
      ELSE
       write(0,*) 'ASH1 ARRAY not here!! STOP!'
       CALL commasmpi_abort()
      ENDIF
    ENDIF

  ENDDO

!-----------------------------------------------------------------------------
! TEMPERATURE perturbations

    dthtot      = 0.0d0
    dthrantot   = 0.0d0
    dthrantotv  = 0.0d0
    dthrantotn  = 0.0d0
    dthrantotp  = 0.0d0
    dthrantotvn = 0.0d0
    dthrantotvp = 0.0d0

#ifdef MPI
    ixb = -ng+1
    ixe = itile+ng
    ixb = 1
    ixe = ixend-ixbeg+1
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = -ng+1
    jye = jtile+ng
    jyb = 1
    jye = jyend-jybeg+1
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg


    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg



    DO k = kzb, kze
     DO j = jyb, jye
      DO i = ixb, ixe
       dth = 0.0

! in the MPI case, dv is used as a flag for the interior of the 
! domain.  Dont want to integrate into ghostzones for normalization
! because that would double-count.
       dv = 0.0
       IF ( i .ge. 1 .and. j .ge. 1 .and. k .ge. 1 ) THEN
       IF (  i .le. itile .and. i .le. ixe .and. &
             j .le. jtile .and. j .le. jye .and. &
             k .le. ktile .and. k .le. kze ) THEN
       dv = 1./(dxc(i)*dyc(j)*dzc(k))
       ENDIF
       ENDIF
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
       dth = 0.0
       dv = 1./(dxc(Min(i,nx-1))*dyc(Min(j,ny-1))*dzc(Min(k,nz-1)))
#endif

       dthran = 0.0d0
       
       DO n = 1,nbble

! "Regular" thermal perturbation bubble

       IF ( bblsp(n) .eq. 0 ) THEN
       
         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

   ! IF radius (X/Y) is large, create a 2D bubble
 
         IF( xrad(n) .gt. 1.0e6 ) &
           radius = sqrt(((yc(j)-ycntr(n))/yrad(n))**2 +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( yrad(n) .gt. 1.0e6 ) &
           radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2 +((zc(k)-zcntr(n))/zrad(n))**2)

   ! Two types of bubbles:  one with a cos^2(x) distribution, the other cos(x)

         IF( radius .le. 1.0 ) THEN
           
           fac = cos(0.5*pii*radius)**2
           
           IF  ( tbble(n) .gt. 0.0 ) THEN
            
             dth = Max( dth, tbble(n)*fac )
           
           ELSE
         
             dth = tbble(n)*fac
           
           ENDIF
#ifdef MPI
           IF( bbletype .eq. 2 ) dthran = dthran + dthranpert*ranarray(ixbeg-1+i,jybeg-1+j)
           IF( bbletype .eq. 3 ) dthran = dthran + dthranpert*ranarray(ixbeg-1+i,jybeg-1+j)*fac
#else
           IF( bbletype .eq. 2 ) dthran = dthran + dthranpert*ranarray(i,j)
           IF( bbletype .eq. 3 ) dthran = dthran + dthranpert*ranarray(i,j)*fac
#endif

          IF ( moisten ) THEN
          
            pres  = psfc*pz(k)**(cp/rd)
            qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
            rh    = qz(k) / qvs
            qvs1  =  380.*exp(17.27*(pz(k)*(tz(k)+dth)-273.16) / (pz(k)*(tz(k)+dth)- 36.)) / pres
            qv3(i,j,k) = rh*qvs1
          
          ENDIF

        ENDIF
 
! Cold pool - dam break style of bryan's cm1

        ELSEIF ( bblsp(n) .eq. 3 ) THEN 
        
         zcntr(n) = 0.   !keep cold pool at sfc
        
         IF ( xc(i).le.xcntr(n) .and. zc(k).le.zrad(n) ) THEN
          dth = tbble(n)*(zrad(n)-zc(k))/zrad(n)
          IF ( xrad(n) <= 0.0 ) u3(i,j,k) = Abs( xrad(n) ) ! set u-component of wind
          IF ( yrad(n) <= 0.0 ) v3(i,j,k) = Abs( yrad(n) ) ! set v-component of wind
         ELSE
          dth = 0.
         END IF
 
! Cold pool - bubble style
 
       ELSEIF ( bblsp(n) .eq. 4 ) THEN 
       
        zcntr(n) = 0.   !keep cold pool at sfc    
        
        IF ( xc(i) .gt. xcntr(n) ) THEN  !2D cylinder
          radius = sqrt( ((xc(i)-xcntr(n))/xrad(n))**2    &
                        +((zc(k)-zcntr(n))/zrad(n))**2 )
        ELSE
          radius = sqrt( ((zc(k)-zcntr(n))/zrad(n))**2 )
        END IF
         
        IF( radius .le. 1.0 ) THEN
          dth = tbble(n)*cos(0.5*pii*radius)**2
        ENDIF
        
!         
        
       ELSEIF ( bblsp(n) .eq. 5 ) THEN  ! hack something here
       
         
         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( radius .le. 1.0 ) THEN

           IF ( microphys%str(1:4) == 'DRYC' ) THEN
              cn(i,j,k) = 1.e-3*(cos(0.5*pii*radius)**2)
           ELSEIF ( microphys%str(1:4) == 'DRYT' ) THEN
              cn(i,j,k) = 1.0 !  1.e-3*(1. - radius) ! *(cos(0.5*pii*radius)**2)
           ELSE
              cn(i,j,k) = cn(i,j,k)*(1.0 - 0.99*cos(0.5*pii*radius)**2)
           ENDIF
           
         ENDIF
       
       ELSEIF ( bblsp(n) .eq. 6 ) THEN  ! hack something here


         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( radius .le. 1.0 ) THEN
           
            qr(i,j,k) = tbble(n)*1.0e-3
            IF ( lnr > 1 ) THEN
              cn(i,j,k) = qr(i,j,k)*1.0e5*pz(k)**2.509/(287.04*tz(k))/(1000.*0.523599*(1.5e-3)**3)
            ENDIF
            IF ( lzr > 1 .and. cn(i,j,k) > 1.e-6 ) THEN
              vr = 1.0e5*pz(k)**2.509/(287.04*tz(k))*qr(i,j,k)/(1000.*cn(i,j,k))
              zr(i,j,k) = 36.*(rnumin+2.0)*cn(i,j,k)*vr**2/((rnumin+1.0)*pii**2)
            ENDIF
!            IF ( GET_VARIABLE_INDEX(gd, 'CRW',1) > 1 ) THEN
!               cn(i,j,k) = qr(i,j,k)*1.0e5*pz(k)**2.509/(287.04*tz(k))/(1000.*0.523599*(0.5e-3)**3)
!            ENDIF
           
         ENDIF

       ELSEIF ( bblsp(n) .eq. 7 ) THEN  ! initialize chaff density
       
         
         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( radius .le. 1.0 ) THEN
           ! assume 30 m/s airspeed, 10**7 fiber/minute release rate
           ! time in box = dx/airspeed
           ! number of fibers = (time in box)*(release rate)
            cn(i,j,k) = tbble(n)/dxc(i)/30. * (1.e7/60.) &
                  *dxc(i)*dyc(j)*dxc(k)*cos(0.5*pii*radius)**2
           
         ENDIF
       
       ELSEIF ( bblsp(n) .eq. 8 ) THEN  ! graupel


         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( radius .le. 1.0 ) THEN
           
          IF ( microphys%str(1:3) == 'TAK' ) THEN
!            cn(i,j,k) = tbble(n)
            cn2(i,j,k) = 0.01*tbble(n)
          ELSE
            qr(i,j,k) = tbble(n)*1.0e-3
            IF ( lnh > 1 ) THEN
              cn(i,j,k) = qr(i,j,k)*1.0e5*pz(k)**2.509/(287.04*tz(k))/(1000.*0.523599*(1.0e-3)**3)
            ENDIF
            IF ( lzh > 1 .and. cn(i,j,k) > 1.e-6 ) THEN
              vr = 1.0e5*pz(k)**2.509/(287.04*tz(k))*qr(i,j,k)/(600.*cn(i,j,k))
              g1 = (6.0)*(5.0)*(4.0)/((3.0)*(2.0)*(1.0))
              zr(i,j,k) = g1*vr**2*cn(i,j,k)
            ENDIF
!            IF ( GET_VARIABLE_INDEX(gd, 'CRW',1) > 1 ) THEN
!               cn(i,j,k) = qr(i,j,k)*1.0e5*pz(k)**2.509/(287.04*tz(k))/(1000.*0.523599*(0.5e-3)**3)
!            ENDIF
          ENDIF ! microphys
         ENDIF
       
! Cold pool - dam break style of bryan's cm1, but reduce QV to keep the same RH (avoids saturating the layer)

        ELSEIF ( bblsp(n) .eq. 9 .or. bblsp(n) .eq. 10 .or. bblsp(n) .eq. 11 ) THEN 
        
         zcntr(n) = 0.   !keep cold pool at sfc
        
         IF ( xc(i) .le. ( xcntr(n) ) .and. zc(k).le.zrad(n) ) THEN
          
          fac = 1.0
!          IF ( xc(i) .gt.  xcntr(n) ) THEN
!           fac = 1. - (xc(i) - xcntr(n))/xrad(n) 
!          ENDIF
          IF ( bblsp(n) .eq. 10 ) THEN
            dth = fac*(tbble(n) + xrad(n) )*(zrad(n)-zc(k)+zc(1))/zrad(n) - xrad(n)
          ELSEIF ( bblsp(n) .eq. 11 ) THEN
            dth = fac*(tbble(n) )*(zrad(n) + xrad(n) - zc(k) + zc(1))/(zrad(n) + xrad(n))
          ELSE
            dth = fac*tbble(n)*(zrad(n)-zc(k)+zc(1))/zrad(n)
          ENDIF
          
          pres  = psfc*pz(k)**(cp/rd)
          qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
          rh    = qz(k) / qvs
          qvs1  =  380.*exp(17.27*(pz(k)*(tz(k)+dth)-273.16) / (pz(k)*(tz(k)+dth)- 36.)) / pres
          qv3(i,j,k) = rh*qvs1

          IF (  xc(i) .ge. xcntr(n) - 25000.) THEN ! add dthranpert degree perts
            dthran = dthran + dthranpert*ranarraysmth(ixbeg-1+i,jybeg-1+j)
           ENDIF

          IF ( xrad(n) <= 0.0 ) u3(i,j,k) = u3(i,j,k) + Abs( xrad(n) ) ! set u-component of wind
          IF ( yrad(n) <= 0.0 ) v3(i,j,k) = Abs( yrad(n) ) ! set v-component of wind
         ELSE
          dth = 0.
         END IF

       ELSEIF ( bblsp(n) .eq. 20 ) THEN  ! initialize ash cloud
       
         
         radius = sqrt(((xc(i)-xcntr(n))/xrad(n))**2  &
                      +((yc(j)-ycntr(n))/yrad(n))**2  &
                      +((zc(k)-zcntr(n))/zrad(n))**2)

         IF( radius .le. 1.0 ) THEN
            cn(i,j,k) = tbble(n)
           
         ENDIF
       
 
        ENDIF ! bblsp(n)

       ENDDO ! n = 1,nbble
           
         IF ( dth .ne. 0.0 ) THEN   ! dth = 0 if outside of radius
           
           th(i,j,k) = tz(k) + dth 


           thpert(i,j,k) = thpert(i,j,k) + dthran

           IF ( dv .gt. 0.0 ) THEN
           dthtot     = dthtot + dth
           dthrantot  = dthrantot + dthran
           dthrantotn = dthrantotn + Min( 0.0d0, dthran )
           dthrantotp = dthrantotp + Max( 0.0d0, dthran )
           
           dthrantotv = dthrantotv + dthran*dv
           dthrantotvn = dthrantotvn + Min( 0.0d0, dthran*dv )
           dthrantotvp = dthrantotvp + Max( 0.0d0, dthran*dv )
           ENDIF

         ENDIF

      ENDDO
     ENDDO
    ENDDO 

#ifdef MPI
!    write(luno,*) 'pre-sum dthrantotv,dthrantot = ',dthrantotv,dthrantot
!    write(luno,*) 'pre-sum dthrantotvp/n = ',dthrantotvp,dthrantotvn
#endif

#ifdef MPI
       mpitotin(1)  =  dthtot
       mpitotin(2)  =  dthrantot
       mpitotin(3)  =  dthrantotn
       mpitotin(4)  =  dthrantotp
       mpitotin(5)  =  dthrantotv
       mpitotin(6)  =  dthrantotvn
       mpitotin(7)  =  dthrantotvp
       
       n = 7

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       dthtot = mpitotout(1)  
       dthrantot = mpitotout(2)  
       dthrantotn = mpitotout(3)  
       dthrantotp = mpitotout(4)  
       dthrantotv = mpitotout(5)  
       dthrantotvn = mpitotout(6)  
       dthrantotvp = mpitotout(7)  
     
#endif

     IF ( my_rank == 0 ) THEN
      write(luno,*) 'Total dth perturbation is ',dthtot
      write(luno,*) 'Total dthran perturbation is ',dthrantot
      write(luno,*) 'Pos/neg dthran perturbation are ',dthrantotp,dthrantotn
     ENDIF

!
! Normalize random perts to zero
!
    IF ( dthrantot .ne. 0.0 ) THEN
    fn = 0.5*(dthrantotvp - dthrantotvn)/Abs(dthrantotvn)
    fp = 0.5*(dthrantotvp - dthrantotvn)/Abs(dthrantotvp)

#ifdef MPI
!    write(lune,*) 'my_rank,fn,fp = ',my_rank,fn,fp,dthrantotvp,dthrantotvn
#endif

    dthrantot   = 0.0d0
    dthrantotv  = 0.0d0
    dthrantotn  = 0.0d0
    dthrantotp  = 0.0d0
    dthrantotvn = 0.0d0
    dthrantotvp = 0.0d0

#ifdef MPI
    ixb = -ng+1
    ixe = itile+ng
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = -ng+1
    jye = jtile+ng
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg


    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    DO k = kzb, kze
     DO j = jyb, jye
      DO i = ixb, ixe
       dv = 0.0
       IF ( i .ge. 1 .and. j .ge. 1 .and. k .ge. 1 ) THEN
       IF (  i .le. itile .and. i .le. ixe .and. &
             j .le. jtile .and. j .le. jye .and. &
             k .le. ktile .and. k .le. kze ) THEN
       dv = 1./(dxc(i)*dyc(j)*dzc(k))
       ENDIF
       ENDIF
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
       dv = 1./(dxc(i)*dyc(j)*dzc(k))
#endif
      
       IF ( thpert(i,j,k) .ge. 0.0 ) THEN
         dthran = fp*thpert(i,j,k)
       ELSE
         dthran = fn*thpert(i,j,k)
       ENDIF
         th(i,j,k) = th(i,j,k) + dthran

           IF ( dv .gt. 0.0 ) THEN
           dthrantot  = dthrantot + dthran
           dthrantotn = dthrantotn + Min( 0.0d0, dthran )
           dthrantotp = dthrantotp + Max( 0.0d0, dthran )
           dthrantotv = dthrantotv + dthran*dv
           dthrantotvn = dthrantotvn + Min( 0.0d0, dthran*dv )
           dthrantotvp = dthrantotvp + Max( 0.0d0, dthran*dv )
           ENDIF
      
      ENDDO
     ENDDO
    ENDDO 

#ifdef MPI
       mpitotin(1)  =  0.0 ! dthtot
       mpitotin(2)  =  dthrantot
       mpitotin(3)  =  dthrantotn
       mpitotin(4)  =  dthrantotp
       mpitotin(5)  =  dthrantotv
       mpitotin(6)  =  dthrantotvn
       mpitotin(7)  =  dthrantotvp
       
       n = 7

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
!       dthtot = mpitotout(1)  
       dthrantot = mpitotout(2)  
       dthrantotn = mpitotout(3)  
       dthrantotp = mpitotout(4)  
       dthrantotv = mpitotout(5)  
       dthrantotvn = mpitotout(6)  
       dthrantotvp = mpitotout(7)  
     
#endif

     IF ( my_rank == 0 ) THEN
!      write(luno,*) 'Norm. Total dth perturbation is ',dthtot
      write(luno,*) 'Norm. Total dthran perturbation is ',dthrantot
      write(luno,*) 'Norm. Pos/neg dthran perturbation are ',dthrantotp,dthrantotn
     ENDIF
    
    ENDIF
!    th(1:nx-1,1:ny-1,1:nz-1) = th(1:nx-1,1:ny-1,1:nz-1) + thpert(1:nx-1,1:ny-1,1:nz-1)

! rebalance hydrostatically for cold pool init
        IF ( Any(bblsp(:) .eq. 3) .or. Any(bblsp(:) .eq. 9) .or. Any(bblsp(:) .eq. 10) &
             .or. Any(bblsp(:) .eq. 11) .or. rebalance ) THEN

#ifdef MPI
    ixb = -ng+1
    ixe = itile+ng
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = -ng+1
    jye = jtile+ng
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg


    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    allocate( pztmp(kze) )

     DO j = jyb, jye
      DO i = ixb, ixe

#else
    kze = nz-1
    kzb = 1
    allocate( pztmp(kze) )
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif

       DO k = kze,kzb,-1
        
        IF ( k .eq. kze ) THEN
          tv0   = th(i,j,kze) * (1.0 + 0.61*qv3(i,j,kze))
          pztmp(kze) = pz(k) ! (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
        ELSE
          tv1   =  th(i,j,k) * (1.0 + 0.61*qv3(i,j,k))
          pztmp(k) = pztmp(k+1) + 2.0*g / ((tv0+tv1)*cp*dze(k+1))
          tv0   = tv1
        ENDIF
        
       ENDDO
        
       DO k = kzb,kze
         pi3(i,j,k) = pztmp(k) - pz(k)
         r3(i,j,k)   = 1.0e5*pztmp(k)**cvr/(rd*th(i,j,k) * (1.0 + 0.61*qv3(i,j,k)))
       ENDDO
        
       ENDDO
      ENDDO
 
      deallocate( pztmp )
      ENDIF

! initialize diagnostic air temperature and pressure, if they are in the grid
     IF ( lairtem >= 1 ) THEN
      CALL GET_VARIABLE(gd,'T',t3)
      IF ( lairpress >= 1 ) CALL GET_VARIABLE(gd,'P',p3)
      IF ( lthetav >= 1 ) CALL GET_VARIABLE(gd,'THETA_V',tv3)
#ifdef MPI
    ixb = -ng+1
    ixe = itile+ng
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = -ng+1
    jye = jtile+ng
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg


    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    DO k = kze,kzb,-1
     DO j = jyb, jye
      DO i = ixb, ixe

#else
    kze = nz-1
    kzb = 1
    DO k = kze,kzb,-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif

         t3(i,j,k) = th(i,j,k)*(pz(k) + pi3(i,j,k))
         IF ( lairpress >= 1 ) THEN
          p3(i,j,k) = 1.0e5*(pz(k) + pi3(i,j,k))**3.509
         ENDIF
         IF ( lthetav >= 1 ) THEN
          tv3(i,j,k) = th(i,j,k)*(1. + 0.61*qv3(i,j,k))
         ENDIF
        
       ENDDO
      ENDDO
      ENDDO
 
     ENDIF
 
 

     IF ( lhwind >= 1 ) THEN
!      write(0,*) 'get hwind, my_rank = ',my_rank
      CALL GET_VARIABLE(gd,'HWIND',hwind3)
    ixb = 1
    ixe = itile
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = 1
    jye = jtile
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg


    kzb = 1
    kze = ktile
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    DO k = kze,kzb,-1
     DO j = jyb, jye
      DO i = ixb, ixe


         ! not adding grid motion here, because u and v are ground-relative at this point
         hwind3(i,j,k) = Sqrt( (0.5*(u3(i,j,k) + u3(i+1,j,k)) )**2 +   &
                               (0.5*(v3(i,j,k) + v3(i,j+1,k)) )**2 )
        
        
       ENDDO
      ENDDO
      ENDDO

     ENDIF


!-----------------------------------------------------------------------------
! CLOSE THE OUTPUT FILE

   luno = FILE_CLOSE('INIT_PERT_BBLE END')
   deallocate ( ranarray )
   deallocate ( thpert )

 END SUBROUTINE INIT_PERT_BBLE


 MODULE TROPICAL_CYCLONE
 
 CONTAINS
 
 
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_TC  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes a bogus vortex for idealized tropical cyclone simulation.
! Adapted from WRF code (ERM)
!
!-----------------------------------------------------------------------
 SUBROUTINE INIT_TC(gd, rtime,nx1,ny1,nz1,ng1,ut,vt)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE
  USE GRIDPARAM_MODULE, only: xdomain,ydomain,lat,dx_stretch

  USE COMMASMPI_MODULE

#ifdef MPI
   use mpi
#endif

   implicit none

! Passed variables

  TYPE(GRID) :: gd
  real, intent(in) :: rtime  ! current model time (is zero for base vortex initialization)
  integer, intent(in) :: nx1,ny1,nz1,ng1
  real, optional, intent(inout) :: ut (-ng1+1:nx1+ng1,-ng1+1:ny1+ng1,-ng1+1:nz1+ng1)
  real, optional, intent(inout) :: vt (-ng1+1:nx1+ng1,-ng1+1:ny1+ng1,-ng1+1:nz1+ng1)
!  integer, intent(in) :: istep


! Local vars

  integer i, j, k, m, l, ii
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer nxin,nyin,nzin,xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes, nltimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)

  real, pointer :: xc(:)
  real, pointer :: xe(:)
  real, pointer :: dxc(:)
  real, pointer :: dxe(:)

  real, pointer :: yc(:)
  real, pointer :: ye(:)
  real, pointer :: dyc(:)
  real, pointer :: dye(:)


  real, pointer :: ccn1d(:)
  real, pointer :: ccnuf1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: p3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)
  real, pointer :: wz(:,:,:)
  real, pointer :: dbz(:,:,:)
  
!  real  :: u3t(nxend,nyend,nzend)
!  real  :: v3t(nxend,nyend,nzend)
!  real  :: t3t(nxend,nyend,nzend)
!  real  :: q3t(nxend,nyend,nzend)

      double precision, allocatable :: zeta(:,:,:) !  real zeta(nxl,nyl,nzl)
      double precision, allocatable :: zetaf(:,:,:) ! real zetaf(nxl,nyl,nzl)
      double precision, allocatable :: psi(:,:,:) !    real psi(-nor:nxl+nor,-nor:nyl+nor,-nor:nzl+nor)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable, save :: ub(:,:,:)
  real, allocatable, save :: vb(:,:,:)
  integer,allocatable :: mask(:,:,:)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
!  double precision, allocatable :: t1(:,:,:)
!  double precision, allocatable :: t2(:,:,:)
!  double precision, allocatable :: t3(:,:,:)
!  real, allocatable, dimension(:) ::  den

      double precision, allocatable ::  vamp(:) ! (nz)
      double precision, allocatable ::  zetan(:) ! (nz)
      double precision, allocatable ::  zetaave(:) ! (nz)
      double precision, allocatable ::  zetafave(:) ! (nz)
      double precision, allocatable ::  zetatmp(:) ! (nz)
      
  real dx,dy,dz,dt
  real xcent,ycent,rad
  integer ix,jy,kz
  integer, parameter  :: irelaxord = 2


      real alfa
      real calfa1
      real calfa2
      real eps, epscheck
      real omega, rjac
      real sum0
      real sum1
      real xmax
      integer :: ipass,isw,ksw,jsw1
      integer :: icent,jcent
      integer, allocatable :: jsw(:,:)
      
      integer :: kk = 1
      double precision, parameter :: pi1 = 3.14159265359
!      double precision, parameter :: pi1 = 3.151592653
      double precision, parameter :: pi2 = 3.14159265359

!      real           vmax,vmax1
      real           vmin,vmin1
      real           umax,umax1
      real           umin,umin1
      integer       itimint
      real  :: fac,tfrac,tmp
      real :: tv0,tv1
      real :: lz
      
      real pztmp(nz1)

!-----------------------------------------------------------------------
!               USER SETTINGS

!  Parameters for analytic vortex:
!  Reference:  Rotunno and Emanuel, 1987, JAS, p. 549

!     real :: r0   !  =   412500.0     ! outer radius (m)
!     real :: rmax  ! =    82500.0     ! approximate radius of max winds (m)
!     real :: vmax  ! =       15.0     ! approximate value of max wind speed (m/s)
     real :: zdd   ! =    20000.0     ! depth of vortex (m)
     
!     real :: dd2
!     integer :: nref, kref

! variables/arrays for analytic vortex:
    integer :: nref,kref,nloop,i1,i2
    real :: dd1,dd2,xref,vr,e1,tx,px,qx,ric,rjc,rr,diff
    real*8 :: rmax,vmax,frac,angle
    real, dimension(:), allocatable :: rref,zref,th0,qv0,thv0,prs0,pi0,rh0
    real, dimension(:,:), allocatable :: vref,piref,pref,thref,thvref,qvref
    
    real :: rannum, pert
    real :: ran0

!  real                 :: zetamax = 13.e-4
!  real                 :: radzero = 280000.
!  real                 :: timint = 0.

! other settings:

    real :: fcor  ! Coriolis parameter (1/s), which is set according to the value of "lat"
!    real :: sst    =  28.0          ! sea-surface temperature (Celsius)

   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: t0           = 300.
   REAL    , PARAMETER :: p0           = p1000mb
   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: r_v          = 461.6
   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL    , PARAMETER ::  EP_2=R_d/R_v

!-----------------------------------------------------------------------
       real fcor1,fcor2
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   integer :: ixb1,jyb1,ixe1,jye1

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag


  IF ( rtime > timint ) THEN
    IF ( allocated( ub ) ) THEN 
      deallocate( ub, vb )
    ENDIF
    RETURN
    
  ENDIF

!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)
  CALL GET_VARIABLE(gd, 'DT', dt)

  CALL GET_VARIABLE(gd, 'XC',   xc)
  CALL GET_VARIABLE(gd, 'XE',   xe)
  CALL GET_VARIABLE(gd, 'DXC',  dxc)
  CALL GET_VARIABLE(gd, 'DXE',  dxe)

  CALL GET_VARIABLE(gd, 'YC',   yc)
  CALL GET_VARIABLE(gd, 'YE',   ye)
  CALL GET_VARIABLE(gd, 'DYC',  dyc)
  CALL GET_VARIABLE(gd, 'DYE',  dye)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
!  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
!  CALL GET_VARIABLE(gd,'CCNIUFNIT',ccnuf1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
  CALL GET_VARIABLE(gd,'WZ', wz)
  CALL GET_VARIABLE(gd,'DBZ', dbz)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'PI',p3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

   ixb = -ng+1
   ixe = itile+ng

   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg+1

!   DO i = ixb,ixe
!     write(luno,*) 'rank,i,xc,dxc = ',my_rank,i,xc(i),dxc(i)
!   ENDDO

       ! assume vortex center at center of domain
       xcent = xdomain/2.
       ycent = ydomain/2.
       
       ! put center squarely on a grid point
       ! need to fix for horizontal stretch!
       ! Actually is OK for horizontal stretch, because the stretch is symmetric.
       IF ( nxend/2 >= ixbeg .and. nxend/2 <= ixend .and. &
            nyend/2 >= jybeg .and. nyend/2 <= jyend  ) THEN
         icent = nxend/2 - ixbeg + 1
         jcent = nyend/2 - jybeg + 1
         xcent = xc(icent)
         ycent = yc(jcent)
         IF ( rtime == 0.0 ) THEN
           write(0,*) 'xcent,ycent = ',xcent,ycent
           write(0,*) 'icent,jcent = ',icent,jcent
         ENDIF
       ELSE
         xcent = -1.
         ycent = -1.
       ENDIF

#ifdef MPI
     ! mpimax for xcent and ycent
       
       mpitotin(1) = xcent
       mpitotin(2) = ycent

       CALL MPI_AllReduce(mpitotin, mpitotout, 2, MPI_DOUBLE_PRECISION,           &
     &                   MPI_MAX, my_comm, mpi_error_code)

       xcent = mpitotout(1)
       ycent = mpitotout(2)
       
#endif

   fcor = 2.0*(7.292e-5)*Sin(degtorad*lat)

     zdd = vorthgt
     rmax = vtrad
     vmax = vtmax
!     r0 = radzero

!-----------------------------------------------------------------------
!  Analytic vortex.
!  Reference:  Rotunno and Emanuel, 1987, JAS, p. 549

    dd2 = 2.0 * rmax / ( r0 + rmax )

   ! nref = 1 + int( float(nxend-nxbeg+1)/2.0 )
    nref = 1 + Int( xdomain/dx_stretch )
    
    kref = nzend-1

!    print *,'  ids,ide,kds,kds  = ',nxbeg,nxend,nzbeg,nzend
!    print *,'  its,ite,kts,kts  = ',ixbeg,ixend,kzbeg,kzend
    IF ( my_rank == 0 ) print *,'  nref,fcor        = ',nref,fcor
    IF ( my_rank == 0 ) print *,'  r0,rmax,vmax,zdd = ',r0,rmax,vmax,zdd

    allocate(  rref(nref)         )
    allocate(  zref(0:kref+1)     )
    allocate(   th0(0:kref+1)     )
    allocate(   qv0(0:kref+1)     )
    allocate(  thv0(0:kref+1)     )
    allocate(  prs0(0:kref+1)     )
    allocate(   pi0(0:kref+1)     )
    allocate(   rh0(0:kref+1)     )
    allocate(  vref(nref,0:kref+1))
    allocate( piref(nref,0:kref+1))
    allocate(  pref(nref,0:kref+1))
    allocate( thref(nref,0:kref+1))
    allocate(thvref(nref,0:kref+1))
    allocate( qvref(nref,0:kref+1))

    ! get base state:
    IF ( my_rank == 0 ) print *,'  zref,th0,qv0,thv0:'
    do k=1,kref
      th0(k) = tz(k)
      qv0(k) = qz(k)
      thv0(k) = th0(k)*(1.0+qv0(k)/epsilon)/(1.0+qv0(k))
      zref(k) =  zc(k) ! 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
      IF ( my_rank == 0 ) print *,k,zref(k),th0(k),qv0(k),thv0(k)
    enddo

    IF ( my_rank == 0 ) print *,'  pi0,rh0:'
    do k=1,kref
!      prs0(k) = grid%p(1,k,1)+grid%pb(1,k,1)
      pi0(k) = pz(k) !  (prs0(k)/p0)**(r_d/cp)
      E1=1000.0*SVP1*EXP(SVP2*(th0(k)*pi0(k)-SVPT0)/(th0(k)*pi0(k)-SVP3))
          pres  = psfc*pz(k)**(cp/rd)
!          qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
       prs0(k) = pres
      qvs = EP_2*E1/(prs0(k)-E1)
      rh0(k) = qv0(k)/qvs
      IF ( my_rank == 0 ) print *,k,pi0(k),rh0(k)
    enddo

    
    zref(0) = -zref(1)
    zref(kref+1) = zref(kref)+(zref(kref)-zref(kref-1))

      rref(:)   = 0.0
      vref(:,:) = 0.0
     piref(:,:) = 0.0
      pref(:,:) = 0.0
     thref(:,:) = 0.0
    thvref(:,:) = 0.0


    do i=1,nref
      rref(i) = dx_stretch*(float(i-1)+0.5)
    enddo

    IF ( my_rank == 0 ) THEN
    print *,'  zref,dz:'
    do k=0,kref+1
      if( k.ge.2 .and. k.le.kref )then
        print *,k,zref(k),zref(k)-zref(k-1)
      else
        print *,k,zref(k)
      endif
    enddo
    ENDIF


    IF ( my_rank == 0 ) print *,'  vref:'
    
    IF ( bogusvortex == 2 ) THEN
    do k=1,kref
      do i=1,nref
        if(rref(i).lt.r0)then
          dd1 = 2.0 * rmax / ( rref(i) + rmax )
          vr = sqrt( vmax**2 * (rref(i)/rmax)**2     &
          * ( dd1 ** 3 - dd2 ** 3 ) + 0.25*fcor*fcor*rref(i)*rref(i) )   &
                  - 0.5 * fcor * rref(i)
        else 
          vr = 0.0
        endif
        if(zref(k).lt.zdd)then
          vref(i,k) = vr * (zdd-zref(k))/(zdd-0.0)
        else
          vref(i,k) = 0.0
        endif
        if(k.eq.1 .and.  my_rank == 0 )  print *,i,rref(i),vref(i,k)
      enddo
    enddo
    
    ELSEIF ( bogusvortex == 3 ) THEN
    
    IF ( vort_lz2 < 0 ) vort_lz2 = vort_lz
    
    do k=1,kref
    
      IF ( zref(k) <= vort_zmax ) THEN
        lz = vort_lz
      ELSE
        lz = vort_lz2
      ENDIF
    
      do i=1,nref
 !       if (rref(i).lt.r0) then
           
           ! modified Rankine (Stern and Nolan 2011)
           IF ( rref(i) <= rmax ) THEN
             
             vr = vmax*rref(i)/rmax
           
           ELSE
           
             vr = vmax*( rmax/rref(i) )**vortexp
           
           ENDIF
           
            vr = vr*exp( -(rref(i)/vort_rcut)**4 )
        
!        else 
!          vr = 0.0
!        endif
        
        dd1 = exp( -( (zref(k) - vort_zmax)/lz)**vort_beta / vort_beta)
        vref(i,k) = vr * dd1
        
        if(k.eq.1 .and.  my_rank == 0 )  print *,i,rref(i),vref(i,k)
      enddo
    enddo

    ENDIF


    IF ( my_rank == 0 ) print *,'  Iterate:'
    DO nloop=1,20

      ! thref is perturbation theta
      ! thvref is "full" theta-v (although divided by (1. + qvref))
      ! piref is pert. Pi
      
      ! get qv and thv from rh and th:
      do k=1,kref
      do i=1,nref
        tx = (pi0(k)+piref(i,k))*(th0(k)+thref(i,k))
        px = p0*((pi0(k)+piref(i,k))**(cp/r_d))
        E1 = 1000.0*SVP1*EXP(SVP2*(tx-SVPT0)/(tx-SVP3))
        qvs = EP_2*E1/(px-E1)
        qvref(i,k) = rh0(k)*qvs ! saturation mixing ratio * RH
        thvref(i,k)=(th0(k)+thref(i,k))*(1.0+qvref(i,k)/epsilon)/(1.0+qvref(i,k))  
                                       
      enddo
      enddo

      ! get nondimensional pressure perturbation (piref):
      do k=1,kref
        piref(nref,k)=0.0
        do i=nref,2,-1 
          piref(i-1,k) = piref(i,k)                                       &
       + (rref(i-1)-rref(i))/(cp*0.5*(thvref(i-1,k)+thvref(i,k))) * 0.5 * &
           ( vref(i  ,k)*vref(i  ,k)/rref(i)                              &
            +vref(i-1,k)*vref(i-1,k)/rref(i-1)                            &
             + fcor * ( vref(i,k) + vref(i-1,k) ) )
        enddo
      enddo

      do i=1,nref
        piref(i,   0) = piref(i, 1)
        piref(i,kref+1) = piref(i,kref)
      enddo

      ! get potential temperature perturbation (thref):
      do k=2,kref
      do i=1,nref
        thref(i,k) = 0.5*( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))/(zref(k+1)-zref(k))     &
                          +cp*0.5*(thvref(i,k)+thvref(i,k-1))*(piref(i,k)-piref(i,k-1))/(zref(k)-zref(k-1)) )   &
                        *thv0(k)/g
        thref(i,k)=(thv0(k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+qvref(i,k)/epsilon)-th0(k)
      enddo
      enddo

      k=1
      do i=1,nref
        thref(i,k) = ( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))/(zref(k+1)-zref(k)) )   &
                        *thv0(k)/g
        thref(i,k)=(thv0(k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+qvref(i,k)/epsilon)-th0(k)
      enddo

      IF ( my_rank == 0 ) print *,'  th,qv,pi = ',nloop,thref(1,1),qvref(1,1),piref(1,1)
      IF ( my_rank == 0 ) print *,'  th,qv,pi = ',nloop,thref(nref-1,1),qvref(nref-1,1),piref(nref-1,1)

    ENDDO   ! enddo for iteration

      IF ( my_rank == 0 ) THEN
        
        print *,' i, thref,qvref,piref,thvref = '
        DO i = 1,nref
         write(*,'(i4,4(2x,1pe13.5))') i,thref(i,1),qvref(i,1),piref(i,1),thvref(i,1)
        ENDDO
      ENDIF


    ! reference (total) pressure:
    do k=1,kref
    do i=1,nref
      pref(i,k) = p0*( ( pi0(k)+piref(i,k) )**(cp/r_d) )
    enddo
    enddo

    ! analytic axisymmetric vortex is ready ... now interpolate to 3D grid:
    ! (note:  vortex is placed in center of domain)

       ixb = -ng+1
       ixe = itile+ng
       ixe1 = itile
       IF ( ixbeg == nxbeg ) ixb = 1
       
       IF (ixend .eq. nxend) THEN
        ixe = ixend-ixbeg+1
        ixe1 = ixe - 1
       ENDIF

       jyb = -ng+1
       IF ( jybeg == nybeg ) jyb = 1
       
       jye = jtile+ng
       jye1 = jtile
       IF (jyend .eq. nyend) THEN 
         jye = jyend-jybeg+1
         jye1 = jye - 1
       ENDIF

!      ixb = 1
!      ixe = itile+1
!      jyb = 1
!      jye = itile+1
      
      DO j = jyb,jye
      DO i = ixb,ixe

      rr = sqrt( (xc(i)-xcent)**2 + (yc(j)-ycent)**2 )
      rr = min( rr , rref(nref) )
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              diff = rref(ii)-rr
            enddo
            i2 = max( ii , 2 )
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))

            do k=1,nz-1
              p3(i,j,k) = piref(i1,k)+(piref(i2,k)-piref(i1,k))*frac ! - pz(k)
              t3(i,j,k) = thref(i1,k)+(thref(i2,k)-thref(i1,k))*frac + tz(k)
              q3(i,j,k) = qvref(i1,k)+(qvref(i2,k)-qvref(i1,k))*frac
            enddo

       IF ( .false. ) THEN
       kzb = 1
       kze = nz-1
       
       DO k = kze,kzb,-1
        
        IF ( k .eq. kze ) THEN
          tv0   = t3(i,j,kze) * (1.0 + 0.601*q3(i,j,kze))
          pztmp(kze) = pz(k) ! (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
        ELSE
          tv1   =  t3(i,j,k) * (1.0 + 0.601*q3(i,j,k))
          pztmp(k) = pztmp(k+1) + 2.0*g / ((tv0+tv1)*cp*dze(k+1))
          tv0   = tv1
        ENDIF
        
       ENDDO
        
       DO k = kzb,kze
         p3(i,j,k) = pztmp(k) - pz(k)
       ENDDO
       
       ENDIF

     ! U-component
      rr = sqrt( (xe(i)-xcent)**2 + (yc(j)-ycent)**2 )
      rr = min( rr , rref(nref) )
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              diff = rref(ii)-rr
            enddo
            i2 = max( ii , 2 )
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))

              angle = datan2(dble(yc(j)-ycent),dble(xe(i)-xcent))

            do k=1,nz

             pert = 0.0
             IF ( iadd_wind_perts /= 0 ) THEN
               rannum = 2.0*( ran0(ipertseed) - 0.5 )
              IF ( rr < vtrad*1.5 ) THEN
                pert = pert_base*rannum
              ENDIF
             ENDIF

              u3(i,j,k) = u1d(k) -( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*sin(angle) + pert
            enddo

     ! V-component
      rr = sqrt( (xc(i)-xcent)**2 + (ye(j)-ycent)**2 )
      rr = min( rr , rref(nref) )
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              diff = rref(ii)-rr
            enddo
            i2 = max( ii , 2 )
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))

              angle = datan2(dble(ye(j)-ycent),dble(xc(i)-xcent))

            do k=1,nz

             pert = 0.0
             IF ( iadd_wind_perts /= 0 ) THEN
               rannum = 2.0*( ran0(ipertseed) - 0.5 )
              IF ( rr < vtrad*1.5 ) THEN
                pert = pert_base*rannum
              ENDIF
             ENDIF

              v3(i,j,k) =  v1d(k) + ( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*cos(angle) + pert
            enddo


      ENDDO ! i
      ENDDO ! j

     deallocate(  rref )
     deallocate(  zref )
     deallocate(   th0 )
     deallocate(   qv0 )
     deallocate(  thv0 )
     deallocate(  prs0 )
     deallocate(   pi0 )
     deallocate(   rh0 )
     deallocate(  vref )
     deallocate( piref )
     deallocate(  pref )
     deallocate( thref )
     deallocate(thvref )
     deallocate( qvref )


 END SUBROUTINE INIT_TC
 
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE HURRFORCE  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes a bogus vortex for idealized tropical cyclone simulation.
! Adapted from SAM code used by A. Fierro
!
!-----------------------------------------------------------------------
 SUBROUTINE HURRFORCE(gd, rtime,nx1,ny1,nz1,ng1,ut,vt)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE
  USE GRIDPARAM_MODULE, only: xdomain,ydomain

  USE COMMASMPI_MODULE

#ifdef MPI
   use mpi
#endif

   implicit none

! Passed variables

  TYPE(GRID) :: gd
  real, intent(in) :: rtime  ! current model time (is zero for base vortex initialization)
  integer, intent(in) :: nx1,ny1,nz1,ng1
  real, optional, intent(inout) :: ut (-ng1+1:nx1+ng1,-ng1+1:ny1+ng1,-ng1+1:nz1+ng1)
  real, optional, intent(inout) :: vt (-ng1+1:nx1+ng1,-ng1+1:ny1+ng1,-ng1+1:nz1+ng1)
!  integer, intent(in) :: istep


! Local vars

  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer nxin,nyin,nzin,xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes, nltimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)

  real, pointer :: xc(:)
  real, pointer :: xe(:)
  real, pointer :: dxc(:)
  real, pointer :: dxe(:)

  real, pointer :: yc(:)
  real, pointer :: ye(:)
  real, pointer :: dyc(:)
  real, pointer :: dye(:)


  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
!  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)
  real, pointer :: wz(:,:,:)
  real, pointer :: dbz(:,:,:)
  
!  real  :: u3t(nxend,nyend,nzend)
!  real  :: v3t(nxend,nyend,nzend)
!  real  :: t3t(nxend,nyend,nzend)
!  real  :: q3t(nxend,nyend,nzend)

      double precision, allocatable :: zeta(:,:,:) !  real zeta(nxl,nyl,nzl)
      double precision, allocatable :: zetaf(:,:,:) ! real zetaf(nxl,nyl,nzl)
      double precision, allocatable :: psi(:,:,:) !    real psi(-nor:nxl+nor,-nor:nyl+nor,-nor:nzl+nor)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable, save :: ub(:,:,:)
  real, allocatable, save :: vb(:,:,:)
  integer,allocatable :: mask(:,:,:)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  double precision, allocatable :: t1(:,:,:)
  double precision, allocatable :: t2(:,:,:)
  double precision, allocatable :: t3(:,:,:)
!  real, allocatable, dimension(:) ::  den

      double precision, allocatable ::  vamp(:) ! (nz)
      double precision, allocatable ::  zetan(:) ! (nz)
      double precision, allocatable ::  zetaave(:) ! (nz)
      double precision, allocatable ::  zetafave(:) ! (nz)
      double precision, allocatable ::  zetatmp(:) ! (nz)
      
  real dx,dy,dz,dt
  real xcent,ycent,rad
  integer ix,jy,kz
  integer, parameter  :: irelaxord = 2


      real alfa
      real calfa1
      real calfa2
      real eps, epscheck
      real omega, rjac
      real sum0
      real sum1
      real xmax
      integer :: ipass,isw,ksw,jsw1
      integer :: icent,jcent
      integer, allocatable :: jsw(:,:)
      
      integer :: kk = 1
      double precision, parameter :: pi1 = 3.14159265359
!      double precision, parameter :: pi1 = 3.151592653
      double precision, parameter :: pi2 = 3.14159265359

      real           vmax,vmax1
      real           vmin,vmin1
      real           umax,umax1
      real           umin,umin1
      integer       itimint
      real  :: fac,tfrac,tmp
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   integer :: ixb1,jyb1,ixe1,jye1

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag


  IF ( rtime > timint ) THEN
    IF ( allocated( ub ) ) THEN 
      deallocate( ub, vb )
    ENDIF
    RETURN
    
  ENDIF

!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)
  CALL GET_VARIABLE(gd, 'DT', dt)

  CALL GET_VARIABLE(gd, 'XC',   xc)
  CALL GET_VARIABLE(gd, 'XE',   xe)
  CALL GET_VARIABLE(gd, 'DXC',  dxc)
  CALL GET_VARIABLE(gd, 'DXE',  dxe)

  CALL GET_VARIABLE(gd, 'YC',   yc)
  CALL GET_VARIABLE(gd, 'YE',   ye)
  CALL GET_VARIABLE(gd, 'DYC',  dyc)
  CALL GET_VARIABLE(gd, 'DYE',  dye)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
!  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
  CALL GET_VARIABLE(gd,'WZ', wz)
  CALL GET_VARIABLE(gd,'DBZ', dbz)
!  CALL GET_VARIABLE(gd,'KM',k3)
!  CALL GET_VARIABLE(gd,'TH',t3)
!  CALL GET_VARIABLE(gd,'QV',q3)

   ixb = -ng+1
   ixe = itile+ng

   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg+1

!   DO i = ixb,ixe
!     write(luno,*) 'rank,i,xc,dxc = ',my_rank,i,xc(i),dxc(i)
!   ENDDO

       ! assume vortex center at center of domain
       xcent = xdomain/2.
       ycent = ydomain/2.
       
       ! put center squarely on a grid point
       ! need to fix for horizontal stretch!
       ! Actually is OK for horizontal stretch, because the stretch is symmetric.
       IF ( nxend/2 >= ixbeg .and. nxend/2 <= ixend .and. &
            nyend/2 >= jybeg .and. nyend/2 <= jyend  ) THEN
         icent = nxend/2 - ixbeg + 1
         jcent = nyend/2 - jybeg + 1
         xcent = xc(icent)
         ycent = yc(jcent)
         IF ( rtime == 0.0 ) THEN
           write(0,*) 'xcent,ycent = ',xcent,ycent
           write(0,*) 'icent,jcent = ',icent,jcent
         ENDIF
       ELSE
         xcent = -1.
         ycent = -1.
       ENDIF

#ifdef MPI
     ! mpimax for xcent and ycent
       
       mpitotin(1) = xcent
       mpitotin(2) = ycent

       CALL MPI_AllReduce(mpitotin, mpitotout, 2, MPI_DOUBLE_PRECISION,           &
     &                   MPI_MAX, my_comm, mpi_error_code)

       xcent = mpitotout(1)
       ycent = mpitotout(2)
       
#endif

      itimint = Int(timint/dt)

      calfa1 = 0.1
      calfa2 = radzero*0.9
      alfa = -alog(calfa1)/(calfa2**4)

       ixe = itile
       ixe1 = itile
       IF (ixend .eq. nxend) THEN
        ixe = ixend-ixbeg+1
        ixe1 = ixe - 1
       ENDIF

       ixb1 = 1
       IF (ixbeg == nxbeg ) ixb1 = 2
       
       jye = jtile
       jye1 = jtile
       IF (jyend .eq. nyend) THEN 
         jye = jyend-jybeg+1
         jye1 = jye - 1
       ENDIF

       jyb1 = 1
       IF (jybeg == nybeg ) jyb1 = 2
       
       IF ( rtime <= 0.0 .or. (rtime < timint .and. .not. allocated(ub) ) ) THEN ! second check is in case of restart while nudging should be happening
       
       IF ( my_rank == 0 ) THEN
         write(0,*) 'Vortex init:'
         write(0,*) 'ixe,jye = ',ixe,jye
         write(0,*) 'ixe1,jye1 = ',ixe1,jye1
         write(0,*) 'ixb1,jyb1 = ',ixb1,jyb1
         write(0,*) 'ixend,nxend,jyend,nyend =',ixend,nxend,jyend,nyend
       ENDIF
       
       allocate (zeta(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (zetaf(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (psi(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (t1(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (t2(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (t3(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (u(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (v(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (ub(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate (vb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))
       allocate( vamp(nz), zetan(nz), zetaave(nz) , zetafave(nz),jsw(nz,2), zetatmp(nz) )
       allocate (mask(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng))

      ksw = 1
      DO i=1,2
        isw = ksw     
         DO k=1,nz
           jsw(k,i) = isw
           isw = 3 - isw
         END DO ! k
        ksw = 3 - ksw
      END DO ! i

        DO kz=1,nz
 !         jsw(kz,1) = 1
 !         jsw(kz,2) = 2
 !         write(0,*) 'kz,jsw1,2 = ',kz,jsw(kz,1),jsw(kz,2)
        ENDDO
         

      zeta(:,:,:) = 0.0
      zetaf(:,:,:) = 0.0
      t2(:,:,:) = 0.0
      psi(:,:,:) = 0.0
      u(:,:,:) = 0.0
      v(:,:,:) = 0.0
      mask(:,:,:) = 0
      
      do kz = 1,nz
      zetan(kz) = 0.0
      zetaave(kz) = 0.0
      zetafave(kz) = 0.0
      do jy = 1,jye
      do ix = 1,ixe

      rad = sqrt( (xc(ix)-xcent)**2 + (yc(jy)-ycent)**2 )

      if ( rad .eq. 0.   ) then
      zeta(ix,jy,kz) = zetamax
      zetaave(kz) = zetaave(kz) + zeta(ix,jy,kz)
      zetaf(ix,jy,kz) = sin(pi1*rad/radzero)
      zetafave(kz) = zetafave(kz) + sin(pi1*rad/radzero)
      end if

      if ( rad .lt. radzero ) then
      zeta(ix,jy,kz) = zetamax*exp( -alfa*(rad**4) )
      zetaave(kz) = zetaave(kz) + zeta(ix,jy,kz)
      zetaf(ix,jy,kz) = sin(pi1*rad/radzero)
      zetafave(kz) = zetafave(kz) + sin(pi1*rad/radzero)
      end if

      if ( rad .ge. radzero ) then
      zeta(ix,jy,kz) = 0.0
      end if


      end do
      end do
      end do
      
! MPI to do: sum up zetafave and zetaave

#ifdef MPI



     ! mpisum zetaave
       
       CALL MPI_AllReduce(zetaave, zetatmp, nz, MPI_DOUBLE_PRECISION,           &
     &                   MPI_SUM, my_comm, mpi_error_code)

       DO k = 1,nz
         zetaave(k) = zetatmp(k)
       ENDDO

     ! mpisum zetafave
       
       CALL MPI_AllReduce(zetafave, zetatmp, nz, MPI_DOUBLE_PRECISION,           &
     &                   MPI_SUM, my_comm, mpi_error_code)

       DO k = 1,nz
         zetafave(k) = zetatmp(k)
       ENDDO
       
#endif

!
!  compute zetan
!
      do kz = 1,nz
      zetaave(kz) = zetaave(kz)/(nxend*nyend)
      zetafave(kz) = zetafave(kz)/(nxend*nyend)
      zetan(kz) = zetaave(kz) / zetafave(kz)
      IF ( my_rank == 0 ) write(0,*) 'z,z,z',zetaave(kz),zetafave(kz),zetan(kz)
      end do
!c
!c  net zero vorticity
!c
      sum0 = 0.0
      sum1 = 0.0
      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe
      sum0 = sum0 + zeta(ix,jy,kz)
      rad = sqrt( (xc(ix)-xcent)**2 + (yc(jy)-ycent)**2 )
      if ( rad .lt. radzero ) then
      zeta(ix,jy,kz) = zeta(ix,jy,kz) - zetaf(ix,jy,kz)*zetan(kz)
      zetaf(ix,jy,kz) = zetaf(ix,jy,kz) * zetan(kz)
      mask(ix,jy,kz) = 1
      else if ( rad .ge. radzero ) then
        zeta(ix,jy,kz) = 0.0
      end if
!      IF ( rad <= xcent - 3.*dx ) mask(ix,jy,kz) = 1
      sum1 = sum1 + zeta(ix,jy,kz)
      end do
      end do
      end do


! MPI to do: sum up sum0 and sum1
! MPI TO DO: communicate edges of mask
#ifdef MPI
      

     ! mpisum
       
       mpitotin(1) = sum0
       mpitotin(2) = sum1

       CALL MPI_AllReduce(mpitotin, mpitotout, 2, MPI_DOUBLE_PRECISION,           &
     &                   MPI_SUM, my_comm, mpi_error_code)

       sum0 = mpitotout(1)
       sum1 = mpitotout(2)
       
#endif

      IF ( my_rank == 0 ) THEN
        write(0,*) 'sum0 zeta', sum0
        write(0,*) 'sum1 zeta', sum1
      ENDIF

#ifdef MPI        
!        IF ( pts == ptend ) THEN
        westward_tag = 201
        CALL sendrecv_westward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,zeta)

        eastward_tag = 202
        CALL sendrecv_eastward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,zeta)

        southward_tag = 203
        CALL sendrecv_southward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,zeta)

        northward_tag = 204
        CALL sendrecv_northward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,zeta)


        westward_tag = 201
        CALL sendrecv_westward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,zetaf)

        eastward_tag = 202
        CALL sendrecv_eastward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,zetaf)

        southward_tag = 203
        CALL sendrecv_southward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,zetaf)

        northward_tag = 204
        CALL sendrecv_northward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,zetaf)

! do not have a comms routine for integer array, so copy mask into t3 temporarily
      t3(:,:,:) = mask(:,:,:)
      
        westward_tag = 201
        CALL sendrecv_westward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,t3)

        eastward_tag = 202
        CALL sendrecv_eastward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,t3)

        southward_tag = 203
        CALL sendrecv_southward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,t3)

        northward_tag = 204
        CALL sendrecv_northward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,t3)
      
      mask(:,:,:) = Nint(t3(:,:,:))
      t3(:,:,:) = 0.0
      
      IF ( my_rank == 0 ) THEN
!       write(0,*) 'mask at nx/2,ny ',mask(nx/2,ny-1,1),mask(nx/2,ny,1),mask(nx/2,ny+1,1)
!       write(0,*) 'mask at nx,ny/2 ',mask(nx-1,ny/2,1),mask(nx,ny/2,1),mask(nx+1,ny/2,1)
      
      ENDIF
      

#endif


!
!  solve del^2 psi = - zeta where psi is the steram function
!
      omega = cos(pi2/nxend)+cos(pi2/nyend)
      omega = 4./(2.+sqrt(4-omega**2))
      rjac = Sqrt( 1.0 - (2.0/omega - 1)**2 )

      IF ( my_rank == 0 ) write(6,*) 'optimal omega is ', omega
!
!      do kz = kk,nz

      omega = cos(pi2/nxend)+cos(pi2/nyend)
      omega = 4./(2.+sqrt(4-omega**2))
!      rjac = Sqrt( 1.0 - (2.0/omega - 1)**2 )
      rjac = (cos(pi2/nxend)+cos(pi2/nyend))/2.0 ! Jacobi interation spectral radius

      do kz = kk,nz-1
      do jy = 0,jye+1
      do ix = 0,ixe+1
      psi(ix,jy,kz) = 0.0
      end do
      end do
      end do
!
!  itertation for sor
!
      do kz = kk,nz-1
      do jy = 1,jye
      do ix = 1,ixe
      t3(ix,jy,kz) = -zeta(ix,jy,kz)/(dxc(ix)*dyc(jy))
      end do
      end do
      end do
!
      if ( irelaxord .eq. 2 ) then
      do m=1,2000
!
!  set t1=t2 (psi_old) and t3=zeta
!
      do kz = kk,nz-1
      Where( mask(:,:,kz) == 0 ) t2(:,:,kz) = 0.0
      end do

      do kz = kk,nz-1
      do jy = 0,jye+1
      do ix = 0,ixe+1
      t1(ix,jy,kz) = t2(ix,jy,kz)
      end do
      end do
      end do
!
!  derivative boudaries
!
     IF ( .false. ) THEN
      IF ( nybeg == jybeg ) THEN
      do ix = 1,ixe
!      t2(ix,1,kz) = 0.0
      t2(ix,0,kz) = 0.0
      end do
      ENDIF

      IF ( nyend == jyend ) THEN
      do ix = 1,ixe
!      t2(ix,ny,kz) = 0.0
      t2(ix,ny+1,kz) = 0.0
      end do
      ENDIF

      IF ( nxbeg == ixbeg ) THEN
      do jy = 1,jye
!      t2(1,jy,kz) = 0.0
      t2(0,jy,kz) = 0.0
      end do
      ENDIF

      IF ( nxend == ixend ) THEN
      do jy = 1,jye
!      t2(nx,jy,kz) = 0.0
      t2(nx+1,jy,kz) = 0.0
      end do
      ENDIF
     ENDIF
      

!
!
!  middle
!
      eps = 1.e-5

      IF ( .true. ) THEN
      ! not sure why, but acceleration does not seem to help (actually slows convergence a little), so turned off for now.
        IF (m.eq.1) THEN
      !        omega = 1.00/(1.00 - 0.50*rjac**2)
          ELSE
      !        omega = 1.00/(1.00 - 0.250*omega*rjac**2 )
          END IF

      ! red-black alternating relaxation
     !  omega = 1.0
      jsw1 = 2
      DO ipass=1,2

           isw = jsw1 ! jsw(kz,ipass)
           jsw1 = 3-jsw1
           ! isw = 3-jsw1
!            DO  jy=jyb1,jye1

        DO kz = kk,nz-1

#ifdef MPI
            DO  jy=jye1,jyb1,-1
                  ibeg=mod((ixbeg-1)+(jy+jybeg-1)+(kz+kzbeg-1)+ipass-1,2) + 1 ! code to ensure that the tile starts
                  iend=2*((nx-ibeg)/2)+ibeg                                   ! processing on a red/black point
                                                                              ! regarless of nx,ny being even/odd.
              DO ix=ibeg,iend,2 ! isw,ixe1,2
#else
            DO  jy=jye1,jyb1,-1
              DO ix=isw+1,ixe1,2
#endif
       t2(ix,jy,kz)=t2(ix,jy,kz) +0.25*omega       &
                 *(t2(ix+1,jy,kz)+t2(ix-1,jy,kz)  &
                  +t2(ix,jy+1,kz)+t2(ix,jy-1,kz)  &
               -4.*t2(ix,jy,kz)                   &
                  +t3(ix,jy,kz))


              ENDDO
             isw = 3 - isw
           ENDDO
          ENDDO
        
! MPI to do: communicate t1 and t2

#ifdef MPI        
!        IF ( ipass == 2 ) THEN
        westward_tag = 201
        CALL sendrecv_westward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,t2)

        eastward_tag = 202
        CALL sendrecv_eastward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,t2)

        southward_tag = 203
        CALL sendrecv_southward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,t2)

        northward_tag = 204
        CALL sendrecv_northward_dp(nx,ny,nz,ng,ng,ng,ng,1,          &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,t2)

!        ENDIF
#endif

      ENDDO
      
      ELSE
       ! serial only: Cannot use this in MPI
      do kz = kk,nz-1
      do jy=jyb1,jye1
      do ix=ixb1,ixe1
      t2(ix,jy,kz)=t2(ix,jy,kz) +0.25*omega       &
                 *(t2(ix+1,jy,kz)+t2(ix-1,jy,kz)  &
                  +t2(ix,jy+1,kz)+t2(ix,jy-1,kz)  &
               -4.*t2(ix,jy,kz)                   &
                  +t3(ix,jy,kz))
      end do
      end do
      end do
      
      ENDIF



!


      IF ( m .gt. 1 ) THEN
      xmax = -1000000.00
      epscheck = 0.0
      OUTER:       do kz = kk,nz-1
      DO jy=jyb1,jye1
      DO ix=ixb1,ixe1
        IF ( mask(ix,jy,kz) == 1 .and. t2(ix,jy,kz) /= 0.0d0 ) THEN
        xmax = max(abs((t2(ix,jy,kz)-t1(ix,jy,kz))/t2(ix,jy,kz)),xmax)
        IF (abs((t2(ix,jy,kz)-t1(ix,jy,kz))/t2(ix,jy,kz)).gt.eps)  THEN
          epscheck = 1.0
!          exit OUTER
        ENDIF
        ENDIF
      end do
      end do
      end do OUTER
      
      ELSE
        epscheck = 1.0
      ENDIF

!      go to 990
!
! MPI to do: reduce all max on epscheck
#ifdef MPI
       
       mpitotin(1) = epscheck
       mpitotin(2) = xmax

       CALL MPI_AllReduce(mpitotin, mpitotout, 2, MPI_DOUBLE_PRECISION,           &
     &                   MPI_MAX, my_comm, mpi_error_code)

       epscheck = mpitotout(1)
       xmax = mpitotout(2)
       
#endif

      IF ( my_rank == 0 .and. ( m < 10 .or.  (m/100)*100 == m ) ) THEN
        write(0,*) 'm,xmax = ',m,xmax
      ENDIF
      IF ( epscheck == 0.0 ) EXIT

      ! test
!        IF ( .false. .and. (m/100)*100 == m .and. m < 1999 ) THEN
!         do jy=jyb1,jye1
!         do ix=ixb1,icent
!           t2(ix,jy,kz)=t2(icent+(icent-ix),jy,kz) 
!         end do
!         end do
!        ENDIF
      
      
      end do ! mloop
      
! 998   continue

! MPI to do: reduce max on xmax
      IF ( my_rank == 0 ) print*,'xmax = ',xmax,'m = ',m,'kz = ',kz
      IF ( .false. .and. my_rank == 0 .and. kz == 1 ) THEN
    !   DO ix = 1,ixe
    !    write(0,*) 'ix,t2 = ',ix,t2(ix,ny/2,kz)
    !   ENDDO
        
        DO jy = jcent-10,jcent+10
          DO ix = icent-10,icent+10
             IF ( ix-icent == jy-jcent ) THEN
                write(0,*) 'ix,jy,psi = ',ix-icent,jy-jcent,t2(ix,jy,kz),t3(ix,jy,kz),zeta(ix,jy,kz)
              ENDIF
           ENDDO
        ENDDO

        DO jy = jcent-10,jcent+10
          DO ix = icent-10,icent+10
             IF ( ix-icent == -jy+jcent ) THEN
                write(0,*) 'ix,jy,psi = ',ix-icent,jy-jcent,t2(ix,jy,kz),t3(ix,jy,kz),zeta(ix,jy,kz)
              ENDIF
           ENDDO
        ENDDO

          jy = jcent
          DO ix = icent-10,icent+10
             write(0,*) 'ix,jy,psi = ',ix-icent,jy-jcent,t2(ix,jy,kz),t3(ix,jy,kz),zeta(ix,jy,kz)
          ENDDO
          DO ix = 1,10
            write(0,*) 'ix,nx-ix+1 ',ix,t2(ix,jy,kz),t2(nx+1-ix,jy,kz)
          ENDDO
        
      ENDIF

      do kz = kk,nz-1
      do jy = 0,jye+1
        do ix = 0,ixe+1
          psi(ix,jy,kz) = t2(ix,jy,kz)
        end do
      end do
      end do

!
!  end do kz=kk,nz
!  end if irelaxord = 2
!
      end if ! irelaxord = 2
      
!      end do ! kz loop

!
!  compute u and v from vamp and psi
!
      vorthgt = Min( vorthgt, 15000.)
      
      do kz = nz-1,kk,-1
      if ( zc(kz) <= vorthgt ) then
      vamp(kz) = 1.0
      elseif ( zc(kz) > vorthgt .and. zc(kz) < 15000. ) then
  !    vamp(kz) = exp(-( ((zc(kz)-vorthgt)/2000.0)**2))
  
!       vamp(kz) = (15001.0-vorthgt-( zc(kz)-vorthgt ) )/ (15001.0-vorthgt)
       vamp(kz) = Max(0.0, 15000.0-vorthgt-( zc(kz)-vorthgt ) )/ (15000.0-vorthgt)
      
      else
        vamp(kz) = 0.0
      endif
   !    IF ( my_rank == 0 ) write(0,*) 'vamp',kz,vamp(kz)
      end do

      vmax = 0.0
      vmin = 0.0
      umax = 0.0
      umin = 0.0
!
! u vel
!
      fac = 1.0 !  1./Float(itimint)
      do kz = kk,nz-1
      do jy = 1,jye
      do ix = 1,itile
      ub(ix,jy,kz) = -vamp(kz)*(psi(ix,jy,kz)-psi(ix,jy-1,kz))*dyc(jy) !  /(dy)
!      u3(ix,jy,kz) = -vamp(kz)*(-1./24.*psi(ix,jy+1,kz)+9./8.*psi(ix,jy,kz)-9./8.*psi(ix,jy-1,kz)+1./24.*psi(ix,jy-2,kz))*dyc(jy) !  /(dy)
      u(ix,jy,kz) = ub(ix,jy,kz)
!      ub(ix,jy,kz) = fac*u3(ix,jy,kz)
      umax = max(umax,u(ix,jy,kz))
      umin = min(umin,u(ix,jy,kz))
      end do
      end do
      end do
!
! v vel
!
      do kz = kk,nz-1
      do jy = 1,jtile
      do ix = 1,ixe
      tmp = vamp(kz)*(psi(ix,jy,kz)-psi(ix-1,jy,kz))*dxc(ix) ! /(dx)
      IF ( .not. ( tmp > -1.e3 .and. tmp < 1.e3 ) ) THEN
        write(0,*) 'problem with v on rank ',my_rank,ix,jy,kz
        write(0,*) 'psi,psi(ix-1),dxc,vamp = ',psi(ix,jy,kz),psi(ix-1,jy,kz),dxc(ix) ,vamp(kz)
        write(0,*) 
        call commasmpi_abort()
      ENDIF
      vb(ix,jy,kz) = tmp
!      v3(ix,jy,kz) = vamp(kz)*(-1./24.*psi(ix+1,jy,kz)+9./8.*psi(ix,jy,kz)-9./8.*psi(ix-1,jy,kz)+1./24.*psi(ix-2,jy,kz))*dyc(jy) !  /(dy)
       v(ix,jy,kz) = vb(ix,jy,kz)
!      vb(ix,jy,kz)= fac*v(ix,jy,kz)
      vmax = max(vmax,v(ix,jy,kz))
      vmin = min(vmin,v(ix,jy,kz))
      end do
      end do
      end do

!      write(0,*) 'my_rank', my_rank,'vmax,vmin = ',vmax,vmin,'umax,umin = ',umax,umin


      ! compute divergence and put into WZ array
      do k = kk,nz-1
      kz=k

      do jy = 1,jye
      do ix = 1,ixe
           ! A-grid velocities
!           t1(ix,jy,k) = -vamp(k)*(psi(ix+1,jy,kz)-psi(ix-1,jy,kz))*0.5*dxc(ix) 
!           t2(ix,jy,k) =  vamp(k)*(psi(ix,jy+1,kz)-psi(ix,jy-1,kz))*0.5*dyc(jy) 
      end do
      end do

      do j = 1,jye
      do i = 1,ixe
!          wz(i,j,k)  =   (u(i+1,j,  k  ) - u(i,j,k))*dxc(i)   &
!                       + (v(i,  j+1,k  ) - v(i,j,k))*dyc(j)

       !   A-grid
!          dbz(i,j,k)  =   ( (t1(i+1,j,  k  ) - t1(i-1,j,k)))*0.5*dxc(i)   &
!                       + (t2(i,  j+1,k  ) - t2(i,j-1,k))*0.5*dyc(j)
                       
!          wz(i,j,k) = zeta(i,j,k)
!          dbz(i,j,k) = psi(i,j,k)
      end do
      end do
      end do

#ifdef MPI
       
       mpitotin(1) = umax
       mpitotin(2) = vmax

       CALL MPI_AllReduce(mpitotin, mpitotout, 4, MPI_DOUBLE_PRECISION,           &
     &                   MPI_MAX, my_comm, mpi_error_code)

       umax = mpitotout(1)
       vmax = mpitotout(2)

       mpitotin(1) = umin
       mpitotin(2) = vmin

       CALL MPI_AllReduce(mpitotin, mpitotout, 4, MPI_DOUBLE_PRECISION,           &
     &                   MPI_MIN, my_comm, mpi_error_code)

       umin = mpitotout(1)
       vmin = mpitotout(2)
       
#endif

      IF ( my_rank == 0 ) print*,'vmax,vmin = ',vmax,vmin,'umax,umin = ',umax,umin



       deallocate (zeta)
       deallocate (zetaf)
       deallocate (psi)
       deallocate (t1)
       deallocate (t2)
       deallocate (t3)
       deallocate( vamp, zetan, zetaave , zetafave, zetatmp )
       deallocate ( u, v )

       ENDIF ! rtime <= 0.0


!      print*,'time',rtime,'forcing u and v to balance p'
      
!      tfrac = float(nstep)/Float(itimint) ! *dtscl/timint
!      tfrac = float(nstep)/Float(Max(1,nstep-1)) ! *dtscl/timint
      IF ( timint < dt ) THEN
         tfrac = 1.0 ! rtime/timint
! v-vel
      do kz = kk,nz-1
      do jy = 1,jtile
      do ix = 1,ixe
        v3(ix,jy,kz) =  vb(ix,jy,kz)
      end do
      end do
      end do

! u vel
!
      do kz = kk,nz-1
      do jy = 1,jye
      do ix = 1,itile
        u3(ix,jy,kz) = ub(ix,jy,kz)
      end do
      end do
      end do

      ELSE
         tfrac = dt/600. !  dt/(0.25 * timint) !  rtime/timint
         IF ( rtime == 0.0 ) tfrac = 0.0
      IF ( my_rank == 0 ) print*, 'time,tfrac,timint = ',rtime,tfrac,timint

      vmax = 0.0
      vmin = 0.0
      umax = 0.0
      umin = 0.0
      vmax1 = 0.0
      vmin1 = 0.0
      umax1 = 0.0
      umin1 = 0.0
      
!
! v vel
!


      IF ( present( vt ) ) THEN
      do kz = kk,nz-1
      do jy = 1,jtile
      do ix = 1,ixe
!        vb(ix,jy,kz) = tfrac*vb(ix,jy,kz)
!        v3(ix,jy,kz) =  vb(ix,jy,kz)
        vt(ix,jy,kz) = vt(ix,jy,kz)  +  tfrac*(vb(ix,jy,kz) - vt(ix,jy,kz) )
!!        v(ix,jy,kz) = tfrac*vb(ix,jy,kz,mv)
!        vc(ix,jy,kz,mv)= v(ix,jy,kz)
!        vd(ix,jy,kz,mv)= v(ix,jy,kz)
!        vn(ix,jy,kz,mv)= v(ix,jy,kz)
!      vmax = max(vmax,v(ix,jy,kz))
!      vmin = min(vmin,v(ix,jy,kz))
!      vmax1 = max(vmax1,vb(ix,jy,kz,mv))
!      vmin1 = min(vmin1,vb(ix,jy,kz,mv))
      end do
      end do
      end do
      ENDIF
!      print*,'vmax1,vmin1 = ',vmax1,vmin1
!      print*,'vmax,vmin = ',vmax,vmin
!
! u vel
!
      IF ( present( ut ) ) THEN
      do kz = kk,nz-1
      do jy = 1,jye
      do ix = 1,itile
        ut(ix,jy,kz) = ut(ix,jy,kz) + tfrac*(ub(ix,jy,kz) - ut(ix,jy,kz))
      end do
      end do
      end do
      ENDIF
!      print*,'umax1,umin1 = ',umax1,umin1
!      print*,'umax,umin = ',umax,umin

      ENDIF




  IF ( rtime >= timint ) THEN
    IF ( allocated( ub ) ) THEN 
      deallocate( ub, vb )
    ENDIF
  ENDIF

      
 END SUBROUTINE HURRFORCE

 END MODULE TROPICAL_CYCLONE

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_IDEAL_DEF  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to user defined background state
!
!-----------------------------------------------------------------------
 SUBROUTINE INIT_IDEAL_DEF(gd) !,                                       &
!                               sndtype, sndfile,                         &
!                               wtype, Us, Uz, psfc, tsfc, qsfc, rhmax, &
!                               shape, dudz0, dudz1, dudz2, z0, z1, z2)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

!  integer sndtype
!  character(LEN = 100) :: sndfile   
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars

  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer nxin,nyin,nzin,xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes, nltimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


!  real lag_thet(62,62,26),lag_qrat(62,62,26)   
!  real lag_u(62,62,26),lag_v(62,62,26)


  real dxin,dyin,dzin,dx,dy,dz
  real grid_xin(62),grid_yin(62),grid_zin(26)  
  real f1,f2,f3,f4,f5,f6,f7,f8
  real grid_x(1000),grid_y(1000),grid_z(1000)  
  real tv(25),ufl(11),vfl(11)
  real gxin,gyin,gx,gy
  real nrint,nlint

  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)


  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)

  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)
  
  real  :: u3t(nxend,nyend,nzend)
  real  :: v3t(nxend,nyend,nzend)
  real :: t3t(nxend,nyend,nzend)
  real  :: q3t(nxend,nyend,nzend)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  


!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    outsoundfile = attr%str(ibeg:iend)//'.output.sound'



  allocate ( lag_thet(62,62,26) )
  allocate ( lag_qrat(62,62,26) )
  allocate ( lag_u(62,62,26) )
  allocate ( lag_v(62,62,26) )

!  allocate( cz(nz) )
  allocate ( den(nz) )

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_IDEAL_DEF BEGIN')

      open(14, file='ideal_def.int',status='unknown')

#ifdef MPI
      read(14,4022)(((t3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(14,4022)(((q3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(14,4022)(((u3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(14,4022)(((v3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      
    do i=ixbeg,ixend
     do j=jybeg,jyend
      do k=kzbeg,kzend  
      
      t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3t(i,j,k)
      q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3t(i,j,k)
      u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3t(i,j,k)                
      v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3t(i,j,k)    
      end do
     end do
    end do
      
#else
      read(14,4022)(((t3(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(14,4022)(((q3(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(14,4022)(((u3(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(14,4022)(((v3(m,l,k),m=1,nx),l=1,ny),k=1,nz)
#endif
      close(14)

!....msb...make format equal to nx

 4022 format (150(f10.4,1x))

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_IDEAL_DEF END')
!-----------------------------------------------------------------------------

 END SUBROUTINE INIT_IDEAL_DEF



!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_IDEAL_DEF_MOV  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to user defined background state
!
!-----------------------------------------------------------------------

 SUBROUTINE INIT_IDEAL_DEF_MOV(gd, x_sw_loc, y_sw_loc)  

!                               sndtype, sndfile,                         &
!                               wtype, Us, Uz, psfc, tsfc, qsfc, rhmax, &
!                               shape, dudz0, dudz1, dudz2, z0, z1, z2)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE
  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

!  integer sndtype
!  character(LEN = 100) :: sndfile   
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars

  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


!  real lag_thet(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_qrat(nxmeso+1,nymeso+1,nzmeso+1) 
!  real lag_u(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_v(nxmeso+1,nymeso+1,nzmeso+1)


  real dx,dy,dz
  real grid_xin(nxmeso+1),grid_yin(nymeso+1)
  real grid_zin(nzmeso+1)  
  real f1,f2,f3,f4,f5,f6,f7,f8
  real grid_x(1000),grid_y(1000),grid_z(1000)  
  real tv(25),ufl(11),vfl(11)
  real gxin,gyin,gx,gy
  real nrint,height(nzmeso),dbl,H,invh
  real x_sw_loc,y_sw_loc
  
  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)


  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)

  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  


!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    outsoundfile = attr%str(ibeg:iend)//'.output.sound'
!  allocate( cz(nz) )
  allocate(den(nz))

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_USER_DEF_MOV BEGIN')


!
!.... allocate the mesoscale domain array abmeso and work arrays
!

      allocate ( abmeso(nxmeso, nymeso, nzmeso, nvbls, nmtimes) )

      allocate  ( athet(nxmeso, nymeso, nzmeso,2) )
      allocate  ( aqrat(nxmeso, nymeso, nzmeso,2) )
      allocate  ( au(nxmeso, nymeso, nzmeso,2) )
      allocate  ( av(nxmeso, nymeso, nzmeso,2) )


      allocate ( lag_thet(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_qrat(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_u(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_v(nxmeso + 1, nymeso + 1, nzmeso + 1) )




     H = 1500.0     
     dbl = 3000.0
     invh = 12000.0



    do n = 1, nmtimes


    do i = 1, nxmeso
    do j = 1, nymeso
    do k = 1, nzmeso

     height(k) = (k - 1) * dzmeso

      if (height(k) .lt. H) then
      abmeso(i,j,k,3,n) = (4.5 - 0.30 * (k - 1))
      else
      abmeso(i,j,k,3,n) = 0.0
      endif
      
      abmeso(i,j,k,4,n) = 1.0 + 4.0 * tanh((i - 110)/1.0)   
     
      if(k .eq. 1) then
      abmeso(i,j,k,1,n) = 311.0 - 1.0 * tanh((i - 110)/1.0)
      abmeso(i,j,k,2,n) = 0.011 + 0.004 * tanh((i - 110)/1.0)

      else
      
      if(height(k) .lt. H) then
      abmeso(i,j,k,1,n) = abmeso(i,j,k-1,1,n)
      abmeso(i,j,k,2,n) = abmeso(i,j,k-1,2,n)
       
      elseif (height(k) .ge. H .and. height(k) .lt. dbl) then
      abmeso(i,j,k,1,n) = 312.00
      abmeso(i,j,k,2,n) = 0.007
      
      elseif (height(k) .ge. dbl .and. height(k) .lt. invh) then     
      abmeso(i,j,k,1,n) = 313.35 + 0.4 * (k - 15)
      abmeso(i,j,k,2,n) = 0.00005 - 0.000001 * (k - 15)

      else   
      abmeso(i,j,k,1,n) = abmeso(i,j,k-1,1,n) + 3.0
      abmeso(i,j,k,2,n) = 0.000001
      endif
      
      endif ! k=1

      abmeso(i,j,k,2,n) = 0.0

    
!       if(i.le.114) then
!       abmeso(i,j,k,1,n)=310.0
!       abmeso(i,j,k,2,n)=0.0002
!       abmeso(i,j,k,3,n)=0.0
!       abmeso(i,j,k,4,n)=-6.0
!             
!       elseif(i.ge.117) then
!       abmeso(i,j,k,1,n)=310.0
!       abmeso(i,j,k,2,n)=0.0002
!       abmeso(i,j,k,3,n)=0.0
!       abmeso(i,j,k,4,n)=12.0
! 
!       else            
! 
!       abmeso(i,j,k,1,n)=310.0
!       abmeso(i,j,k,2,n)=0.0002
!       abmeso(i,j,k,3,n)=0.0
!       abmeso(i,j,k,4,n)=-6.0 + 6*(i-114)
!             
!       end if
          
      enddo
      enddo
      enddo
 
 
      enddo




    do i = 1, nxmeso
    do j = 1, nymeso
    do k = 1, nzmeso
    
    lag_thet(i,j,k) = abmeso(i,j,k,1,1)
    lag_qrat(i,j,k) = abmeso(i,j,k,2,1)
    lag_u(i,j,k) = abmeso(i,j,k,3,1)
    lag_v(i,j,k) = abmeso(i,j,k,4,1)

	  if(i.eq.2.and.j.eq.2) then
      print*,lag_qrat(i,j,k),lag_thet(i,j,k),k,j,'thet'
      endif

    enddo
    enddo
    enddo
    
    
      do i=1,nxmeso +1
      do j=1,nymeso +1
      do k=1,nzmeso +1
      lag_thet(nxmeso+1,j,k)=lag_thet(nxmeso,j,k)
      lag_thet(i,nymeso+1,k)=lag_thet(i,nymeso,k)
      lag_thet(i,j,nzmeso+1)=lag_thet(i,j,nzmeso)
      lag_qrat(nxmeso+1,j,k)=lag_qrat(nxmeso,j,k)
      lag_qrat(i,nymeso+1,k)=lag_qrat(i,nymeso,k)
      lag_qrat(i,j,nzmeso+1)=lag_qrat(i,j,nzmeso)
      lag_u(nxmeso+1,j,k)=lag_u(nxmeso,j,k)
      lag_u(i,nymeso+1,k)=lag_u(i,nymeso,k)
      lag_u(i,j,nzmeso+1)=lag_u(i,j,nzmeso)
      lag_v(nxmeso+1,j,k)=lag_v(nxmeso,j,k)
      lag_v(i,nymeso+1,k)=lag_v(i,nymeso,k)
      lag_v(i,j,nzmeso+1)=lag_v(i,j,nzmeso)
      enddo
      enddo
      enddo

!--------------------------------------------------
!
! Tri-linear interpolation of input data to model grid
! msb 10/30/07
!
! dxmeso,dymeso,dzmeso => input grid spacing
! nxmeso,nymeso,nzmeso => input grid dimensions
!


  do i=1,nxmeso+1
  grid_xin(i)=(real(i)-1)*dxmeso
  enddo

  do i=1,nymeso+1
  grid_yin(i)=(real(i)-1)*dymeso
  enddo

  do i=1,nzmeso+1
  grid_zin(i)=(real(i)-1)*dzmeso
  enddo

#ifdef MPI
  do i=nxbeg,nxend
#else  
  do i=1,nx
#endif 
  grid_x(i)=(real(i)-1)*dx + x_sw_loc
  enddo

#ifdef MPI
  do i=nybeg,nyend
#else  
  do i=1,ny
#endif 
  grid_y(i)=(real(i)-1)*dy + y_sw_loc
  enddo

!  do i=1,ny
!  grid_y(i)=(real(i)-1)*dy
!  enddo


!...msb 5/28/08 change for stretched vertical grid

#ifdef MPI
  grid_z(nzbeg)=0.0  
  do i=nzbeg+1,nzend
#else  
  grid_z(1)=0.0
  do i=2,nz
#endif 
  grid_z(i)=grid_z(i-1)+(1/dzc(i-1))
  enddo



   gxin=(real(nxmeso)-1)*dxmeso
   gyin=(real(nymeso)-1)*dymeso

#ifdef MPI   
   gx=(real(nxend)-1)*dx
   gy=(real(nyend)-1)*dy
   xmx=(gxin/gx)*(nxend-1) + 1
   ymx=(gyin/gy)*(nyend-1) + 1
#else   
   gx=(real(nx)-1)*dx
   gy=(real(ny)-1)*dy
   xmx=(gxin/gx)*(nx-1) + 1
   ymx=(gyin/gy)*(ny-1) + 1   
#endif



#ifdef MPI
   do i=ixbeg,ixend
    do j=jybeg,jyend
     do k=kzbeg,kzend
#else     
   do i=1,nx
    do j=1,ny
     do k=1,nz
#endif

#ifdef MPI
!  ix1= (real(nxmeso)-1)/(real(nxend)-1)*(gx/gxin)*(i-1) + 1

   ix1=int(grid_x(i)/dxmeso) +1
  
  if(ix1.gt.nxmeso) then
  ix1=nxmeso
  end if
  ix2=ix1+1

!  iy1= (real(nymeso)-1)/(real(nyend)-1)*(gy/gyin)*(j-1) + 1

  iy1 = int(grid_y(i)/dymeso) +1

  if(iy1.gt.nymeso) then
  iy1=nymeso
  end if
  iy2= iy1+1

#else 
    
!  ix1= (real(nxmeso)-1)/(real(nx)-1)*(gx/gxin)*(i-1) + 1

   ix1=int(grid_x(i)/dxmeso) +1
  
  if(ix1.gt.nxmeso) then
  ix1=nxmeso
  end if
  ix2=ix1+1
  
!  iy1= (real(nymeso)-1)/(real(ny)-1)*(gy/gyin)*(j-1) + 1

   iy1 = int(grid_y(i)/dymeso) +1
  
  if(iy1.gt.nymeso) then
  iy1=nymeso
  end if
  iy2= iy1+1
  
#endif

    

!.....msb 5/28/08 change for stretched vertical grid
   iz1=grid_z(k)/dzmeso + 1
  iz2= iz1 +1
!  print*,iz1,iz2,grid_z(k),dzmeso
!  iz1= (real(nzmeso)-1)/(real(nz)-1)*(k-1) + 1  
!  iz2= iz1 +1
!    print*,dzmeso,1/dzc(k),grid_z(k)/dzmeso,iz1,iz2,k

!     print*,i,j,k,ix1,ix2,iy1,iy2,iz1,iz2,'n' 

     f1=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
     f2=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
     f3=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
     f4=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))

     f5=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
     f6=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
     f7=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
     f8=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))

!     print*,f1,f2,f3,f4,f5,f6,f7,f8,'h'

#ifdef MPI

     u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_u(ix1,iy1,iz1) + f2*lag_u(ix2,iy1,iz1) + f3*lag_u(ix1,iy2,iz1)+ &
      f4*lag_u(ix2,iy2,iz1) +f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

     v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
      f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

     t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+f3*lag_thet(ix1,iy2,iz1)+ &
       f4*lag_thet(ix2,iy2,iz1)+f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+        &
       f8*lag_thet(ix2,iy2,iz2)

     q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+f3*lag_qrat(ix1,iy2,iz1)+ &
      f4*lag_qrat(ix2,iy2,iz1)+f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#else

     u3(i,j,k)=f1*lag_u(ix1,iy1,iz1)+f2*lag_u(ix2,iy1,iz1)+f3*lag_u(ix1,iy2,iz1)+f4*lag_u(ix2,iy2,iz1)+ &
     f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

     v3(i,j,k)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
     f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

     t3(i,j,k)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+f3*lag_thet(ix1,iy2,iz1)+f4*lag_thet(ix2,iy2,iz1)+ &
     f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+f8*lag_thet(ix2,iy2,iz2)

     q3(i,j,k)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+f3*lag_qrat(ix1,iy2,iz1)+f4*lag_qrat(ix2,iy2,iz1)+ &
     f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#endif


!       if(k.eq.1.and.my_rank.eq.1) then
!       print*, v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1), t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
!        print*,lag_v(ix1,iy1,iz1),lag_v(ix2,iy1,iz1),lag_v(ix1,iy2,iz1),lag_v(ix2,iy2,iz1), &
!        lag_v(ix1,iy1,iz2),lag_v(ix2,iy1,iz2),lag_v(ix1,iy2,iz2),lag_v(ix2,iy2,iz2),  &
!        ix1,ix2,iy1,iy2,iz1,iz2
!       print*,i-ixbeg+1,j-jybeg+1,k-kzbeg+1,i,j,k,ixbeg,jybeg,kzbeg
!
!       end if

    enddo
   enddo
  enddo

!---------------------------------------------------
!
! End tri-linear interpolation


! msb 02/26/10
! Extend grid

#ifdef MPI
  do i=ixbeg,ixend
  do j=jybeg,jyend
  do k=kzbeg,kzend
#else
  do i=1,nx
  do j=1,ny
  do k=1,nz
#endif

#ifdef MPI

  if(i.gt.xmx.and.j.le.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,j,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,j,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,j,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,1,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(i,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(i,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(i,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(i,1,k)
  end if
 
#else
 
  if(i.gt.xmx.and.j.le.ymx) then
  u3(i,j,k)=u3(xmx,j,k)
  v3(i,j,k)=v3(xmx,j,k)
  t3(i,j,k)=t3(xmx,j,k)
  q3(i,j,k)=q3(xmx,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(xmx,ymx,k)
  v3(i,j,k)=v3(xmx,ymx,k)
  t3(i,j,k)=t3(xmx,ymx,k)
  q3(i,j,k)=q3(xmx,ymx,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(i,ymx,k)
  v3(i,j,k)=v3(i,ymx,k)
  t3(i,j,k)=t3(i,ymx,k)
  q3(i,j,k)=q3(i,ymx,k)
  end if

#endif

  enddo
  enddo
  enddo



!---------------------------
!
! Option to eliminate superadiabatic layer at the boundaries
!

    iadiabat=0

    if(iadiabat.eq.1) then

   do j=1,ny
    do i=1,nx
    tv(25)=t3(i,j,25)*(1 + 0.61*q3(i,j,25))
     do k=nz-1,1,-1

      tv(k)=t3(i,j,k)*(1 + 0.61*q3(i,j,k))

      if((i.eq.1).or.(i.eq.nx).or.(j.eq.1).or.(j.eq.ny)) then
      if(tv(k).gt.tv(k+1)) then

      tv(k)=tv(k+1)
      t3(i,j,k)=tv(k)/(1 + 0.61*q3(i,j,k))
!      print*,t3(i,j,k),tv(k),q3(i,j,k),i,j,k

      end if
      end if

     enddo
    enddo
   enddo

   endif

 

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_USER_DEF END')
!-----------------------------------------------------------------------------

 END SUBROUTINE INIT_IDEAL_DEF_MOV




!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_DATA3D_DEF_MOV  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to user-defined MESOSCALE background state prescribed from real data (e.g., soundings)
!
!-----------------------------------------------------------------------

 SUBROUTINE INIT_DATA3D_DEF_MOV(gd, x_sw_loc, y_sw_loc)  


  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

!  integer sndtype
!  character(LEN = 100) :: sndfile   
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars


  integer jybl0
  real xblw0
  integer ixblw
  real xble01
  real xble02
  real xble03
  real xble04
  real xble05
  real xble06

  real rotang
  real xmeso
  real xwgt
  real xblw
  real xble1
  real xble2
  real xble3
  real xble4
  real xble5
  real xble6



  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile
  character(LEN=100)         :: insndfile
  integer                    :: ibeg,iend







!  real lag_thet(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_qrat(nxmeso+1,nymeso+1,nzmeso+1) 
!  real lag_u(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_v(nxmeso+1,nymeso+1,nzmeso+1)


  real dx
  real dy
  real dz
  real grid_xin(nxmeso+1)
  real grid_yin(nymeso+1)
  real grid_zin(nzmeso+1)  
  real f1
  real f2
  real f3
  real f4
  real grid_x(1000)
  real grid_y(1000)
  real grid_z(1000)  
  real tv(25)
  real ufl(11)
  real vfl(11)
  real gxin,gyin,gx,gy
  real nrint,height(nzmeso),dbl,H,invh
  real x_sw_loc,y_sw_loc
  
  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)


  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)

  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den

!.... arrays now defined in init_module
!    abmeso(i,j,k,1,n) = tz_meso(k)
!    abmeso(i,j,k,2,n) = qz_meso(k)
!    abmeso(i,j,k,3,n) = u1d_meso(k)
!    abmeso(i,j,k,4,n) = v1d_meso(k)

  real           :: uzwork(nzmeso)
  real           :: vzwork(nzmeso)
  real           :: tzwork(nzmeso)
  real           :: qzwork(nzmeso)



!--------------------------------------------------------------------------
! local vars and code cannibalized from internal subprogram SND2

! Define arrays for reading in the sounding


    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax)
    real tv0, tv1
!    real z0, tv0, tv1
    real z, an, ugnd

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac



!--------------------------------------------------------------------------




!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  

!-----------------------------------------------------------------------








!
!.... begin executable code
!









! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)

!
!.... retrieve base-state sounding data arrays (see grid_create.F90 for definitions)
!

  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)



  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
!    outsoundfile = attr%str(ibeg:iend)//'.output.sound'


!  allocate( cz(nz) )
  allocate(den(nz))



!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_DATA3D_DEF_MOV BEGIN')
 


!
!.... allocate the mesoscale domain array abmeso and work arrays
!

      allocate ( abmeso(nxmeso, nymeso, nzmeso, nvbls, nmtimes) )

      allocate  ( athet(nxmeso, nymeso, nzmeso,2) )
      allocate  ( aqrat(nxmeso, nymeso, nzmeso,2) )
      allocate  ( au(nxmeso, nymeso, nzmeso,2) )
      allocate  ( av(nxmeso, nymeso, nzmeso,2) )


      allocate ( lag_thet(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_qrat(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_u(nxmeso + 1, nymeso + 1, nzmeso + 1) )
      allocate ( lag_v(nxmeso + 1, nymeso + 1, nzmeso + 1) )

!
!.... CLZ (3-12-13): mesoscale heterogeneous environmental sounding arrays
!

      allocate ( tz_meso(nzmeso, maxsnd) )
      allocate ( qz_meso(nzmeso, maxsnd) )
      allocate ( u1d_meso(nzmeso, maxsnd) )
      allocate ( v1d_meso(nzmeso, maxsnd) )


      H = 1500.0     
      dbl = 3000.0
      invh = 12000.0


!
!.... CLZ (3-12-13):  CAVEAT EMPTOR! jybl0 is same for all soundings in an east-west array
!

      jybl0 = ((ysnd + 0.01)/dymeso) + 1.0

!
!.... compute mesoscale x-coordinates of ixblw
!

      xblw0 = xsnd
      ixblw = ((xsnd + 0.01)/dxmeso) + 1.0




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------




!--------------------------------------------------------------------------
! INITIALIZE HETEROGENEOUS SOUNDING (ASSUMES SNDTYPE = 2)

      abmeso(:,:,:,:,:) = 0.0
      u1d_meso(:,:) = 0.0
      v1d_meso(:,:) = 0.0
      tz_meso(:,:) = 0.0
      qz_meso(:,:) = 0.0

!
!.... CLZ (3-12-13):  read in the (# = numinhom) mesoscale environmental soundings
!

      do j = 1, numinhom

      insndfile = sndinhom(j)

      write(luno,*)
      write(luno,*)
      write(luno,*)'j=',j,' reading in sndinhom=',insndfile

      uzwork(:) = 0.0
      vzwork(:) = 0.0
      tzwork(:) = 0.0
      qzwork(:) = 0.0

!      if (numinhom .ge. 1) xble01 = xsndinhom(1)
!      if (numinhom .ge. 2) xble02 = xsndinhom(2)
!      if (numinhom .ge. 3) xble03 = xsndinhom(3)
!      if (numinhom .ge. 4) xble04 = xsndinhom(4)
!      if (numinhom .ge. 5) xble05 = xsndinhom(5)
!      if (numinhom .ge. 6) xble06 = xsndinhom(6)

      if (j .eq. 1) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom1.snd'
      xble01 = xsndinhom(j)

      elseif (j .eq. 2) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom2.snd'
      xble02 = xsndinhom(j)

      elseif (j .eq. 3) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom3.snd'
      xble03 = xsndinhom(j)

      elseif (j .eq. 4) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom4.snd'
      xble04 = xsndinhom(j)

      elseif (j .eq. 5) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom5.snd'
      xble05 = xsndinhom(j)

      elseif (j .eq. 6) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom6.snd'
      xble06 = xsndinhom(j)

      endif

      CALL SNDREAD_INHOM3D(insndfile,outsoundfile,nzmeso,uzwork,vzwork,tzwork,qzwork)



!.... CLZ (3/4/13): assume nzmeso = nz - 1

      do k = 1, nzmeso

!     if (k .lt. nzmeso) then

     u1d_meso(k,j)   = uzwork(k)      ! U
     v1d_meso(k,j)   = vzwork(k)      ! V
     tz_meso(k,j)   = tzwork(k)      ! TH
     qz_meso(k,j)   = qzwork(k)      ! QV

!     write(luno,*) 'k=',k,' u1d_meso=',u1d_meso(k)
!     write(luno,*) 'k=',k,' v1d_meso=',v1d_meso(k)
!     write(luno,*) 'k=',k,' tz_meso=',tz_meso(k)
!     write(luno,*) 'k=',k,' qz_meso=',qz_meso(k)

!
!.... enddo k = 1, nzmeso
!

     enddo

!
!.... enddo j = 1, numinhom
!

     enddo



!--------------------------------------------------------------------------
! DONE INITIALIZING HETEROGENEOUS SOUNDINGS





!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------




      write(luno,*)
      write(luno,*)'ivelh1d=',ivelh1d,' irot=',irot
      write(luno,*)


!
!.... rotdeg is input rotation of isolines (deg from North)
!

    if (irot .eq. 1) then

    rotang = degtorad * rotdeg

!
!.... endif (irot .eq. 0) then
!

    endif



!.... now assign abmeso field values from series of interpolated 1-D soundings

    do n = 1, nmtimes


    do k = 1, nzmeso

!
!.... CLZ (7/30/09): rotate wind vector to the ground-relative flow streamline
!

    if (ivelh1d .eq. 0 .and. irot .eq. 0) then

!.... CLZ (3-11-13):  horizontally homogeneous winds from base-state sounding

    rotang = atan(u1d(k)/v1d(k))
!    rotang = atan(vbz(kz,mu)/vbz(kz,mv))    

    elseif (ivelh1d .gt. 0 .and. irot .eq. 0) then

!.... CLZ (3-14-13):  horizontally homogeneous winds from selected eastern sounding

    rotang = atan(u1d_meso(k,ivelh1d)/v1d_meso(k,ivelh1d))

    endif

      write(luno,*)'k=',k,' rotang=',rotang


!.... Looping through j, compute x-coordinates of gradient edges (xblw and xble)

    do j = 1, nymeso

!.... rotate west/east edge coordinates of each gradient with wind

      xblw = xblw0 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 1) xble1 = xble01 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 2) xble2 = xble02 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 3) xble3 = xble03 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 4) xble4 = xble04 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 5) xble5 = xble05 + float(j - jybl0) * dymeso * tan(rotang)
      if (numinhom .ge. 6) xble6 = xble06 + float(j - jybl0) * dymeso * tan(rotang)

!    if(j .eq. 1) then
!    write(luno,*)'xblw=',xblw,' xble1=',xble1,' xble2=',xble2,' xble3=',xble3
!    endif


    do i = 1, nxmeso


!.... compute x-coordinate in mesoscale domain, xmeso

      xmeso = (i - 1) * dxmeso

!.... abmeso (x, y, z, 1, n) - theta
!.... abmeso (x, y, z, 2, n) - qv


!.... West of gradient zone, use western (base-state) sounding.

      if (xmeso .le. xblw) then

      abmeso(i, j, k, 1, n) = tz(k)
      abmeso(i, j, k, 2, n) = qz(k)


!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*) 'k=',k,'< xblw=',xblw,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif

!.... In 1st gradient zone, weight linearly from western to eastern (1) sounding.

      elseif (xmeso .gt. xblw .and. xmeso .le. xble1) then

      xwgt = (xmeso - xblw) / (xble1 - xblw)
      abmeso(i, j, k, 1, n) = (1.0 - xwgt) * tz(k) + xwgt * tz_meso(k, 1)
      abmeso(i, j, k, 2, n) = (1.0 - xwgt) * qz(k) + xwgt * qz_meso(k, 1)

!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*)'xwgt=',xwgt
!    write(luno,*)'k=',k,'> xblw=',xblw,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif

!.... In 2nd gradient zone, weight linearly from eastern (1) to eastern (2) sounding.

      elseif (xmeso .gt. xble1 .and. xmeso .le. xble2) then

      xwgt = (xmeso - xble1) / (xble2 - xble1)
      abmeso(i, j, k, 1, n) = (1.0 - xwgt) * tz_meso(k, 1) + xwgt * tz_meso(k, 2)
      abmeso(i, j, k, 2, n) = (1.0 - xwgt) * qz_meso(k, 1) + xwgt * qz_meso(k, 2)

!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*)'xwgt=',xwgt
!    write(luno,*)'k=',k,'> xble1=',xble1,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif

!.... In 3rd gradient zone, weight linearly from eastern (2) to eastern (3) sounding.

      elseif (xmeso .gt. xble2 .and. xmeso .le. xble3) then

      xwgt = (xmeso - xble2) / (xble3 - xble2)
      abmeso(i, j, k, 1, n) = (1.0 - xwgt) * tz_meso(k, 2) + xwgt * tz_meso(k, 3)
      abmeso(i, j, k, 2, n) = (1.0 - xwgt) * qz_meso(k, 2) + xwgt * qz_meso(k, 3)

!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*)'xwgt=',xwgt
!    write(luno,*)'k=',k,'> xble2=',xble2,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif

!.... In 4th gradient zone, weight linearly from eastern (3) to eastern (4) sounding.

      elseif (xmeso .gt. xble3 .and. xmeso .le. xble4) then

      xwgt = (xmeso - xble3) / (xble4 - xble3)
      abmeso(i, j, k, 1, n) = (1.0 - xwgt) * tz_meso(k, 3) + xwgt * tz_meso(k, 4)
      abmeso(i, j, k, 2, n) = (1.0 - xwgt) * qz_meso(k, 3) + xwgt * qz_meso(k, 4)

!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*)'xwgt=',xwgt
!    write(luno,*)'k=',k,'> xble2=',xble2,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif


!.... East of easternmost gradient zone, use eastern (4) sounding (i.e., 5th in array).

      elseif (xmeso .gt. xble4) then

      abmeso(i, j, k, 1, n) = tz_meso(k, 4)
      abmeso(i, j, k, 2, n) = qz_meso(k, 4)

!    if(i .eq. 2 .and. j .eq. 2) then
!    write(luno,*)'k=',k,'> xble3=',xble3,' abmeso(i,j,k,1,1)=',abmeso(i,j,k,1,1) 
!    endif

      endif



!
!.... CLZ (2/5/07): debug printout of abmeso theta and qv.
!
!      if(jy.eq.20.and.kz.eq.10) then
!      write(6,*) 'ix = ',ix,', abmeso(pt)/abmeso(qv) = ', abmeso(ix,jy,kz,1),abmeso(ix,jy,kz,2)
!      endif


      if (ivelh1d .eq. 0) then

!.... CLZ (3-11-13):  horizontally homogeneous winds from base-state sounding

      abmeso(i, j, k, 3, n) = u1d(k)
      abmeso(i, j, k, 4, n) = v1d(k)
    
      elseif (ivelh1d .gt. 0) then

!.... CLZ (3-14-13):  horizontally homogeneous winds from ivelh1d-th selected sounding

      abmeso(i, j, k, 3, n) = u1d_meso(k,ivelh1d)
      abmeso(i, j, k, 4, n) = v1d_meso(k,ivelh1d)

!.... CLZ (3-27-14):  reset base-state wind profile from ivelh1d-th selected sounding

      u1d(k) = u1d_meso(k,ivelh1d)
      v1d(k) = v1d_meso(k,ivelh1d)

      endif

!
!.... enddo k = 1, nzmeso loop through mesoscale domain prescribing theta and qv
!
    enddo
    enddo
    enddo

!
!.... CLZ (3/20/13): test mandatory elimination of heterogeneities above k = khomog
!

      do k = khomog, nzmeso
      do j = 1, nymeso
      do i = 1, nxmeso 

      abmeso(i, j, k, 1, n) = tz(k)
      abmeso(i, j, k, 2, n) = qz(k)

      enddo
      enddo
      enddo



    enddo

!.... CLZ (11/10/11): done with prescribing sounding-based environment (abmeso)




    do i = 1, nxmeso
    do j = 1, nymeso
    do k = 1, nzmeso
    
    lag_thet(i,j,k) = abmeso(i,j,k,1,1)
    lag_qrat(i,j,k) = abmeso(i,j,k,2,1)
    lag_u(i,j,k) = abmeso(i,j,k,3,1)
    lag_v(i,j,k) = abmeso(i,j,k,4,1)

	  if(i.eq.2.and.j.eq.2) then
	  write(luno,*) k,'qrat=',lag_qrat(i,j,k),' theta=',lag_thet(i,j,k)
      print*,lag_qrat(i,j,k),lag_thet(i,j,k),k,j,'thet'
      endif

    enddo
    enddo
    enddo
    
    
      do i=1,nxmeso +1
      do j=1,nymeso +1
      do k = 1, nzmeso +1
      lag_thet(nxmeso + 1, j, k) = lag_thet(nxmeso, j, k)
      lag_thet(i,nymeso+1,k)=lag_thet(i,nymeso,k)
      lag_thet(i,j,nzmeso+1)=lag_thet(i,j,nzmeso)
      lag_qrat(nxmeso+1,j,k)=lag_qrat(nxmeso,j,k)
      lag_qrat(i,nymeso+1,k)=lag_qrat(i,nymeso,k)
      lag_qrat(i,j,nzmeso+1)=lag_qrat(i,j,nzmeso)
      lag_u(nxmeso+1,j,k)=lag_u(nxmeso,j,k)
      lag_u(i,nymeso+1,k)=lag_u(i,nymeso,k)
      lag_u(i,j,nzmeso+1)=lag_u(i,j,nzmeso)
      lag_v(nxmeso+1,j,k)=lag_v(nxmeso,j,k)
      lag_v(i,nymeso+1,k)=lag_v(i,nymeso,k)
      lag_v(i,j,nzmeso+1)=lag_v(i,j,nzmeso)
      enddo
      enddo
      enddo
    






!--------------------------------------------------
!
! Tri-linear interpolation of input data to model grid
! msb 10/30/07
!

#ifdef MPI
  do i = nxbeg, nxend
#else  
  do i = 1, nx
#endif 
  grid_x(i) = x_sw_loc + (real(i) - 1.0) * dx
  enddo

#ifdef MPI
  do i = nybeg, nyend
#else  
  do i = 1, ny
#endif 
  grid_y(i) = y_sw_loc + (real(i) - 1.0) * dy
  enddo


!.....CLZ (3/26/13):  mesoscale and model grid levels coincide vertically

#ifdef MPI
  do i = nzbeg, nzend - 1
#else  
  do i = 1, nz - 1
#endif
  grid_z(i) = zc(i)
  enddo


! dxmeso,dymeso,dzmeso => input grid spacing
! nxmeso,nymeso,nzmeso => input grid dimensions
!


  do i = 1, nxmeso + 1
  grid_xin(i) = (real(i) - 1.0) * dxmeso
  enddo

  do i = 1, nymeso + 1
  grid_yin(i) = (real(i) - 1.0) * dymeso
  enddo

! CLZ (3/1/13): check value of nzmeso

    write(luno,*)
    write(luno,*) 'nzmeso=',nzmeso

!
!  soundings, abmeso, and model grid now share same vertical levels
!

!    write(luno,*)
!    write(luno,*) 'print common levels of mesoscale and model grids'


  do i = 1, nzmeso + 1

   grid_zin(i) = grid_z(i)
!    write(luno,*)'at level=',i,'grid_zin(i)=',grid_zin(i)

  enddo




   gxin = (real(nxmeso) - 1.0) * dxmeso
   gyin = (real(nymeso) - 1.0) * dymeso

#ifdef MPI   
   gx=(real(nxend)-1)*dx
   gy=(real(nyend)-1)*dy
   xmx=(gxin/gx)*(nxend-1) + 1
   ymx=(gyin/gy)*(nyend-1) + 1
#else   
   gx=(real(nx)-1)*dx
   gy=(real(ny)-1)*dy
   xmx=(gxin/gx)*(nx-1) + 1
   ymx=(gyin/gy)*(ny-1) + 1   
#endif



#ifdef MPI
   do i=ixbeg,ixend
    do j=jybeg,jyend

     do k=kzbeg,kzend
#else     
   do i=1,nx
    do j=1,ny

     do k=1,nz
#endif

#ifdef MPI
!  ix1= (real(nxmeso)-1)/(real(nxend)-1)*(gx/gxin)*(i-1) + 1

   ix1 = int(grid_x(i)/dxmeso) + 1
  
  if (ix1 .gt. nxmeso) then
  ix1 = nxmeso
  endif

  ix2 = ix1 + 1

!  iy1= (real(nymeso)-1)/(real(nyend)-1)*(gy/gyin)*(j-1) + 1

  iy1 = int(grid_y(j)/dymeso) + 1

  if (iy1 .gt. nymeso) then
  iy1 = nymeso
  endif

  iy2 = iy1 + 1

#else 
    
!  ix1= (real(nxmeso)-1)/(real(nx)-1)*(gx/gxin)*(i-1) + 1

   ix1 = int(grid_x(i)/dxmeso) + 1
  
  if (ix1 .gt. nxmeso) then
  ix1 = nxmeso
  endif

  ix2=ix1+1
  
!  iy1= (real(nymeso)-1)/(real(ny)-1)*(gy/gyin)*(j-1) + 1

   iy1 = int(grid_y(j)/dymeso) + 1
  
  if (iy1 .gt. nymeso) then
  iy1 = nymeso
  endif

  iy2 = iy1 + 1
  
#endif

    

!.....msb 5/28/08 change for stretched vertical grid

!
! CLZ (3/1/20): grids coincident in vertical, so use index k
!

   iz1 = k
   iz2 = iz1 + 1
!   iz1 = (grid_z(k) / dzmeso) + 1
!   iz2= iz1 + 1

! CLZ (3/1/13): new indices iz1 and iz2?????

!    if (i .eq. ixbeg .and. j.eq. jybeg) then
!    write(luno,*)
!    write(luno,*) 'grid_z(k)=',grid_z(k),'grid_zin(k)=',grid_zin(k)

! CLZ (3/1/13): back-check for bad iz1 or iz2 > nzmeso + 1 (61)

!    write(luno,*) 'iz1=',iz1,' iz2=',iz2
!    endif


!  print*,iz1,iz2,grid_z(k),dzmeso
!  iz1= (real(nzmeso)-1)/(real(nz)-1)*(k-1) + 1  
!  iz2= iz1 +1
!    print*,dzmeso,1/dzc(k),grid_z(k)/dzmeso,iz1,iz2,k

!     print*,i,j,k,ix1,ix2,iy1,iy2,iz1,iz2,'n' 




! CLZ (3/1/13): horizontal bilinear interpolation

      f1 = ( (grid_xin(ix2) - grid_x(i)) * (grid_yin(iy2) - grid_y(j)) )/(dxmeso * dymeso)
      f2 = ( (grid_x(i) - grid_xin(ix1)) * (grid_yin(iy2) - grid_y(j)) )/(dxmeso * dymeso)
      f3 = ( (grid_xin(ix2) - grid_x(i)) * (grid_y(j) - grid_yin(iy1)) )/(dxmeso * dymeso)
      f4 = ( (grid_x(i) - grid_xin(ix1)) * (grid_y(j) - grid_yin(iy1)) )/(dxmeso * dymeso)

!     f1=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
!     f2=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
!     f3=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))
!     f4=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxmeso*dymeso*dzmeso))

!     f5=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
!     f6=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
!     f7=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))
!     f8=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxmeso*dymeso*dzmeso))

!     print*,f1,f2,f3,f4,f5,f6,f7,f8,'h'

#ifdef MPI


     u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_u(ix1,iy1,iz1) + f2*lag_u(ix2,iy1,iz1) + &
           f3*lag_u(ix1,iy2,iz1) + f4*lag_u(ix2,iy2,iz1) 
     v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_v(ix1,iy1,iz1) + f2*lag_v(ix2,iy1,iz1) + &
           f3*lag_v(ix1,iy2,iz1) + f4*lag_v(ix2,iy2,iz1)
     t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_thet(ix1,iy1,iz1) + f2*lag_thet(ix2,iy1,iz1) + &
           f3*lag_thet(ix1,iy2,iz1) + f4*lag_thet(ix2,iy2,iz1)
     q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_qrat(ix1,iy1,iz1) + f2*lag_qrat(ix2,iy1,iz1) + &
           f3*lag_qrat(ix1,iy2,iz1) + f4*lag_qrat(ix2,iy2,iz1)

!     u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = f1*lag_u(ix1,iy1,iz1) + f2*lag_u(ix2,iy1,iz1) + f3*lag_u(ix1,iy2,iz1)+f4*lag_u(ix2,iy2,iz1) + &
!     f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

!     v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
!     f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

!     t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+f3*lag_thet(ix1,iy2,iz1)+f4*lag_thet(ix2,iy2,iz1)+ &
!     f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+f8*lag_thet(ix2,iy2,iz2)

!     q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+f3*lag_qrat(ix1,iy2,iz1)+f4*lag_qrat(ix2,iy2,iz1)+ &
!     f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#else

     u3(i,j,k) = f1*lag_u(ix1,iy1,iz1) + f2*lag_u(ix2,iy1,iz1) + f3*lag_u(ix1,iy2,iz1) + f4*lag_u(ix2,iy2,iz1)
     v3(i,j,k) = f1*lag_v(ix1,iy1,iz1) + f2*lag_v(ix2,iy1,iz1) + f3*lag_v(ix1,iy2,iz1) + f4*lag_v(ix2,iy2,iz1)
     t3(i,j,k) = f1*lag_thet(ix1,iy1,iz1) + f2*lag_thet(ix2,iy1,iz1) + f3*lag_thet(ix1,iy2,iz1) + f4*lag_thet(ix2,iy2,iz1)
     q3(i,j,k) = f1*lag_qrat(ix1,iy1,iz1) + f2*lag_qrat(ix2,iy1,iz1) + f3*lag_qrat(ix1,iy2,iz1) + f4*lag_qrat(ix2,iy2,iz1)


!     u3(i,j,k)=f1*lag_u(ix1,iy1,iz1)+f2*lag_u(ix2,iy1,iz1)+f3*lag_u(ix1,iy2,iz1)+f4*lag_u(ix2,iy2,iz1)+ &
!     f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

!     v3(i,j,k)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
!     f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

!     t3(i,j,k)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+f3*lag_thet(ix1,iy2,iz1)+f4*lag_thet(ix2,iy2,iz1)+ &
!     f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+f8*lag_thet(ix2,iy2,iz2)

!     q3(i,j,k)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+f3*lag_qrat(ix1,iy2,iz1)+f4*lag_qrat(ix2,iy2,iz1)+ &
!     f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#endif


!       if(k.eq.1.and.my_rank.eq.1) then
!       print*, v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1), t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
!        print*,lag_v(ix1,iy1,iz1),lag_v(ix2,iy1,iz1),lag_v(ix1,iy2,iz1),lag_v(ix2,iy2,iz1), &
!        lag_v(ix1,iy1,iz2),lag_v(ix2,iy1,iz2),lag_v(ix1,iy2,iz2),lag_v(ix2,iy2,iz2),  &
!        ix1,ix2,iy1,iy2,iz1,iz2
!       print*,i-ixbeg+1,j-jybeg+1,k-kzbeg+1,i,j,k,ixbeg,jybeg,kzbeg
!
!       end if

    enddo
   enddo
  enddo

!---------------------------------------------------
!
! End tri-linear interpolation


! msb 02/26/10
! Extend grid

#ifdef MPI
  do i=ixbeg,ixend
  do j=jybeg,jyend
  do k=kzbeg,kzend
#else
  do i=1,nx
  do j=1,ny
  do k=1,nz
#endif

#ifdef MPI

  if(i.gt.xmx.and.j.le.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,j,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,j,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,j,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,1,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(i,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(i,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(i,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(i,1,k)
  end if
 
#else
 
  if(i.gt.xmx.and.j.le.ymx) then
  u3(i,j,k)=u3(xmx,j,k)
  v3(i,j,k)=v3(xmx,j,k)
  t3(i,j,k)=t3(xmx,j,k)
  q3(i,j,k)=q3(xmx,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(xmx,ymx,k)
  v3(i,j,k)=v3(xmx,ymx,k)
  t3(i,j,k)=t3(xmx,ymx,k)
  q3(i,j,k)=q3(xmx,ymx,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(i,ymx,k)
  v3(i,j,k)=v3(i,ymx,k)
  t3(i,j,k)=t3(i,ymx,k)
  q3(i,j,k)=q3(i,ymx,k)
  end if

#endif

  enddo
  enddo
  enddo

 

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_DATA3D_DEF_MOV END')
!-----------------------------------------------------------------------------



!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM DEFINITIONS (Fortran90 stuff....)

  CONTAINS



!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SNDREAD_INHOM3D:  Reads in a sounding from an input file and interpolates
!                            the sounding to the grid.

   SUBROUTINE SNDREAD_INHOM3D(insndfile,outsoundfile,nzmeso,uzwork,vzwork,tzwork,qzwork)

! Define arrays for reading in the sounding

    implicit none


    integer nzmeso

    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax)
    real uzwork(nzmeso)
    real vzwork(nzmeso)
    real tzwork(nzmeso)
    real qzwork(nzmeso)

!
!.... internal work array
!

    real pzwork(nzmeso)
    
    real z0, tv0, tv1
    real z, an, ugnd

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac
    character(LEN=255)         :: outsoundfile
    character(LEN=100)         :: insndfile


! Error check


!------------------------------------------------------------------------------
! First read in the surface pressure, temperature, and mixing ratio
! Format is z,th,qv,u,v because of stretched grid
!     Sounding input file format: ascii
!     z (meters), theta (K), qv (g/kg),  u (m/s), v (m/s)
!------------------------------------------------------------------------------

    CALL FILE_MESSAGE('READING SOUNDING FROM: '//insndfile)

!
!.... open the input sounding unit
!

    open(unit= 17, file= insndfile, status= 'old', form= 'formatted')
!    open(unit= 17, file= sndfile, status= 'old', form= 'formatted')
    rewind(17)

! Read in sfc pressure, theta, qv

    read(17,*) p000, t000, q000

    write(luno,*) p000, t000, q000

   psfc = p000
   tsfc = t000
   qsfc = q000
   qfac = 1.0
   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qfac = 0.001
   
   qsfc = qfac*qsfc


! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! Read in sounding values

    nsnd = 0
    DO k = 1,nmax
	 read(17,*,iostat=ios) zsnd(k),tsnd(k),qvsnd(k),usnd(k),vsnd(k)
	 qvsnd(k) = qfac*qvsnd(k)
       IF( ios < 0 ) EXIT
	 nsnd = nsnd + 1
    ENDDO

!
!.... close the input sounding unit
!

    close(17)




#ifdef MPI
    if (kzbeg .eq. nzbeg)  zsfc = zsnd(1)
#else
    zsfc = zsnd(1)
#endif

   IF( .not. SET_VARIABLE(gd,'HGT',       zsfc) ) write(6,*) 'INIT_GRID:  Problem setting HGT'

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .le. nsnd ) kze = nsnd-kzbeg
    DO k = kzb,kze
#else
    DO k = 1,nsnd
#endif
      zsnd(k) = zsnd(k) - zsfc
    ENDDO

    zsnd(nsnd) = zsnd(nsnd) + 0.1

! Check to make sure the sounding levels span the computational grid

#ifdef MPI
    if (kzbeg .eq. nzbeg .AND. zsnd(1) .gt. zc(1) ) THEN
#else
    IF( zsnd(1) .gt. zc(1) ) THEN
#endif
	write(luno,*) 'ZMIN of SOUNDING   = ',zsnd(1)
	write(luno,*) 'ZMIN of MODEL GRID = ',zc(1)
	write(luno,*) 'ZMIN OF INPUT SOUNDING IS > ZMIN OF GRID'
	write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
	stop
    ENDIF

#ifdef MPI
    IF( zsnd(nsnd) .lt. zc(nzend-1) ) THEN
#else
    IF( zsnd(nsnd) .lt. zc(nz-1) ) THEN
#endif
	write(luno,*) 'ZMAX of SOUND      = ',zsnd(nsnd)
	write(luno,*) 'ZMAX of MODEL GRID = ',zc(nz-1)
	write(luno,*) 'ZMAX OF INPUT SOUNDING IS < ZMAX OF MODEL GRID'
	write(luno,*) 'Extrapolating sounding with isothermal layer....'
!      write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
!      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
!      stop
    ENDIF

! Now interpolate sounding to grid

    write(luno,*)
    write(luno,*) 'INTERPOLATING SOUNDING TO GRID '

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg
    DO k = kzb,kze
#else
    DO k = 1,nz-1
#endif
      z0 = zc(k)

      IF ( z0 .lt. zsnd(nsnd) ) THEN

!    real uzwork(nzmeso)
!    real vzwork(nzmeso)
!    real tzwork(nzmeso)
!    real qzwork(nzmeso)

        uzwork(k) = zfind(zsnd, usnd,  nsnd, z0)
        vzwork(k) = zfind(zsnd, vsnd,  nsnd, z0)
        tzwork(k) = zfind(zsnd, tsnd,  nsnd, z0)
        qzwork(k) = zfind(zsnd, qvsnd, nsnd, z0)

!        u1d(k) = zfind(zsnd, usnd,  nsnd, z0)
!        v1d(k) = zfind(zsnd, vsnd,  nsnd, z0)
!        tz(k) = zfind(zsnd, tsnd,  nsnd, z0)
!        qz(k) = zfind(zsnd, qvsnd, nsnd, z0)

        ktop = k

        IF ( k .eq. 1 ) THEN
          tv0   = tzwork(1) * (1.0 + 0.601*qzwork(1))
!          tv0   = tz(1) * (1.0 + 0.601*qz(1))
          pzwork(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
        ELSE
          tv1   = tzwork(k) * (1.0 + 0.601*qzwork(k))
!          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pzwork(k) = pzwork(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
          tmptr = pzwork(k)*tzwork(k)
!          tmptr = pz(k)*tz(k)
        ENDIF

      ELSE

        uzwork(k) = uzwork(ktop)
        vzwork(k) = vzwork(ktop)
        qzwork(k) = qzwork(ktop)
        tzwork(k) = tzwork(ktop) * (Exp((z0 - zc(ktop) - zsnd(1)) * g/cp/tmptr ) )

!        u1d(k) = u1d(ktop)
!        v1d(k) = v1d(ktop)
!        qz(k) = qz(ktop)
!        tz(k) = tz(ktop)*(Exp((z0-zc(ktop) - zsnd(1))*g/cp/tmptr ) )

          tv1   = tzwork(k) * (1.0 + 0.601*qzwork(k))
!          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pzwork(k) = pzwork(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1

      ENDIF

!      IF( qz(k) .gt. 0.05 ) qz(k) = qz(k) / 1000.

    ENDDO


! Write sounding out to file

    IF ( ny .le. 2 ) THEN
      v1d(:) = 0.0
    ELSEIF ( nx .le. 2 ) THEN
      u1d(:) = 0.0
    ENDIF
    

#ifdef MPI
  if(my_rank == 0) then
#endif
    open(21,file=outsoundfile, status='unknown')
    rewind(21)
!    write(21,*) 'MODEL SND'
    write(21,'(1x,2(f9.2,2x),f9.4)') psfc, t000, qsfc*1000.

    DO k = 1,nz-1
!     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), tzwork(k),1000.*qzwork(k),uzwork(k),vzwork(k),pzwork(k)
!     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),psfc
!     write(21,'(1x, 6(f10.2,2x))') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),p000

    ENDDO
    close(21)

#ifdef MPI
  endif
#endif
   END SUBROUTINE SNDREAD_INHOM3D



 END SUBROUTINE INIT_DATA3D_DEF_MOV




!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_DATA3D_DEF_FRONT  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain with 2 soundings to create an E-W oriented stationary front
! Added 9/19/2019 by RLM
!
!-----------------------------------------------------------------------

 SUBROUTINE INIT_DATA3D_DEF_FRONT(gd, x_sw_loc, y_sw_loc)  


  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

!  integer sndtype
!  character(LEN = 100) :: sndfile   
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars


  integer jybl0
  real xblw0
  integer ixblw
  real xble01
  real xble02
  real xble03
  real xble04
  real xble05
  real xble06

  real xmeso
  real ymeso,yemeso
  real yble01
  real xwgt
  real xblw
  real xble1
  real xble2
  real xble3
  real xble4
  real xble5
  real xble6
!  real sf_north
!  real sf_south
!  real cden
!  real hden
!  real dencorr
  real linwgt
  real ywgt,yewgt
  real ycoord_wgt,yecoord_wgt
  real pi



  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile
  character(LEN=100)         :: insndfile
  integer                    :: ibeg,iend



  real dx
  real dy
  real dz
  real grid_xin(nxmeso+1)
  real grid_yin(nymeso+1)
  real grid_zin(nzmeso+1)  
  real f1
  real f2
  real f3
  real f4
  real grid_x(1000)
  real grid_y(1000)
  real grid_z(1000)  
  real tv(25)
  real ufl(11)
  real vfl(11)
  real gxin,gyin,gx,gy
  real nrint,height(nzmeso),dbl,H,invh
  real x_sw_loc,y_sw_loc
  real :: yshift
  
  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
  real, pointer :: xc(:)
  real, pointer :: xe(:)
  real, pointer :: yc(:)
  real, pointer :: ye(:)
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)


  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)

  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)

  real, allocatable :: cz(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den


  real :: uzwork(nzmeso + 1)
  real :: vzwork(nzmeso + 1)
  real :: tzwork(nzmeso + 1)
  real :: qzwork(nzmeso + 1)
  real :: pzwork(nzmeso + 1)
  
!  real, pointer  :: u1d_meso(nz,numinhom)
!  real, pointer  :: v1d_meso(nz,numinhom)
!  real, pointer  :: tz_meso(nz,numinhom)
!  real, pointer  :: qz_meso(nz,numinhom)



!--------------------------------------------------------------------------
! local vars and code cannibalized from internal subprogram SND2

! Define arrays for reading in the sounding


    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax)
    real tv0, tv1
!    real z0, tv0, tv1
    real z, an, ugnd

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac



!--------------------------------------------------------------------------




!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  

!-----------------------------------------------------------------------

!
!.... begin executable code
!

! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'XC',    xc)
  CALL GET_VARIABLE(gd,'XE',    xe)
  CALL GET_VARIABLE(gd,'YC',    yc)
  CALL GET_VARIABLE(gd,'YE',    ye)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)

!
!.... retrieve base-state sounding data arrays (see grid_create.F90 for definitions)
!

  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)



  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)
  
  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)



!  allocate( cz(nz) )
  allocate(den(nz))



!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_DATA3D_DEF_FRONT BEGIN')
  
     allocate ( tz_meso(nz, maxsnd) )
     allocate ( qz_meso(nz, maxsnd) )
     allocate ( u1d_meso(nz, maxsnd) )
     allocate ( v1d_meso(nz, maxsnd) )
 

      H = 1500.0     
      dbl = 3000.0
      invh = 12000.0
      pi = 3.1415927

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!--------------------------------------------------------------------------
! INITIALIZE HETEROGENEOUS SOUNDING (ASSUMES SNDTYPE = 2)

!      abmeso(:,:,:,:,:) = 0.0
      u1d_meso(:,:) = 0.0
      v1d_meso(:,:) = 0.0
      tz_meso(:,:) = 0.0
      qz_meso(:,:) = 0.0

!
!.... CLZ (3-12-13):  read in the (# = numinhom) mesoscale environmental soundings
!

      do j = 1, numinhom

      insndfile = sndinhom(j)

      write(luno,*)
      write(luno,*)
      write(luno,*)'j=',j,' reading in sndinhom=',insndfile

      uzwork(:) = 0.0
      vzwork(:) = 0.0
      tzwork(:) = 0.0
      qzwork(:) = 0.0


      if (j .eq. 1) then
      outsoundfile = attr%str(ibeg:iend)//'.mdl_inhom1.snd'
      yble01 = ysndinhom(j)

      elseif (j .gt. 1) then
      write(luno,*) 'oops, for some reason it thinks it wants another sounding'

      
      endif

     ! CALL SNDREAD_INHOM3D(insndfile,outsoundfile,nzmeso,uzwork,vzwork,tzwork,qzwork)

      CALL SND2X(gd,tzwork,pzwork,qzwork,dzc(1:nzmeso+1),dze(1:nzmeso+1),zc(1:nzmeso+1),nzmeso+1,0,psfc,tsfc,qsfc,rhmaxin,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,uzwork,vzwork,z0,z1,dudz0,rotdeg, &
                   outsoundfile,insndfile)   ! IF SND_TYPE = 1, initialize with WK sounding + wind profile

     open(21,file=outsoundfile,form='formatted', status='unknown', position='append')
!   rewind(21)
!   write(21,*) psfc, tsfc, qsfc*1000.
     DO k = 1,nzend-1
      write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), tzwork(k), qzwork(k)*1000., uzwork(k), vzwork(k), &
                  100.*1000.*pz(k)**3.508
     ENDDO

   close(21)



!.... CLZ (3/4/13): assume nzmeso = nz - 1

      do k = 1, nz

!     if (k .lt. nzmeso) then

     u1d_meso(k,j)   = uzwork(k)      ! U
     v1d_meso(k,j)   = vzwork(k)      ! V
     tz_meso(k,j)   = tzwork(k)      ! TH
     qz_meso(k,j)   = qzwork(k)      ! QV


!
!.... enddo k = 1, nzmeso
!

     enddo

!
!.... enddo j = 1, numinhom
!

     enddo



!--------------------------------------------------------------------------
! DONE INITIALIZING HETEROGENEOUS SOUNDINGS





!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------



!.... now put in the stationary front using the base state and input sounding


    do k = 1, nz
    
    yshift = 0.0
    IF ( yshiftmag /= 0.0 .and. zc(k) < yshiftdepth ) THEN
      yshift = yshiftmag*(1.0 - zc(k)/yshiftdepth)
  !    write(0,*) 'k,z,yshift = ',k,zc(k),yshift
    ENDIF
    

    do j = 1, ny 


    do i = 1, nx


!... RLM define north and south side of domain with cool and warm sounding

       ymeso = yc(j) ! (j - 1) * dy
       yemeso = ye(j) ! = yc(j) - dy/2.0
       
!   calculate average u value that way inflow component is the same -> no horizontal shearing instability

      u3(i, j, k) = (u1d(k) + u1d_meso(k,1)) / 2.0
       
!     Cool side of sounding (northern 0215 sounding)       

      if (ymeso .ge. sf_north - yshift) then

        t3(i, j, k) = tz(k)
        q3(i, j, k) = qz(k)
        v3(i, j, k) = v1d(k)

         ! catch the V point just south of sf_north
         IF (yemeso < sf_north - yshift) THEN
          yecoord_wgt = (yemeso - (sf_south - yshift)) / (sf_north - sf_south)
          yewgt = (sin((pi * ycoord_wgt) + (pi/2.0)) + 1.0) / 2.0
          if (yewgt .gt. 1.0) yewgt = 1.0
          if (yewgt .lt. 0.0) yewgt = 0.0
          v3(i, j, k) = (1.0 - yewgt) * v1d(k) + yewgt * v3(i, j, k)
         
         ENDIF
      
!     Warm side of sounding (splice 0300 and 0215 sounding) 
      
      elseif (ymeso .le. sf_south - yshift) then
        t3(i, j, k) = tz_meso(k, 1)
        q3(i, j, k) = qz_meso(k, 1)
        
  !     Density speed correction to v      
          IF ( dencorr == 1 ) THEN
            linwgt = 1.0 - (zc(k) / hden)
            if (linwgt .lt. 0.0) linwgt = 0.0
          ELSE
            linwgt = 0.0
          ENDIF
        v3(i, j, k) = (linwgt * cden) + v1d_meso(k, 1)
   
      
!.... In 1st gradient zone, weight linearly from western to eastern (1) sounding.

      elseif (ymeso .gt. sf_south-yshift .and. ymeso .le. sf_north-yshift ) then
      
      !.... adjust v component that gets put in the sin wave to match environment
        IF ( dencorr == 1 ) THEN
          linwgt = 1.0 - (zc(k) / hden)
          if (linwgt .lt. 0.0) linwgt = 0.0
        ELSE
          linwgt = 0.0
        ENDIF
      v3(i, j, k) = (linwgt * cden) + v1d_meso(k, 1)
      
      !.... calculate weights using sine function      

      ycoord_wgt = (ymeso - (sf_south-yshift) ) / (sf_north - sf_south)
      ywgt = (sin((pi * ycoord_wgt) + (pi/2.0)) + 1.0) / 2.0
      if (ywgt .gt. 1.0) ywgt = 1.0
      if (ywgt .lt. 0.0) ywgt = 0.0

      yecoord_wgt = (yemeso - (sf_south-yshift)) / (sf_north - sf_south)
      yewgt = (sin((pi * yecoord_wgt) + (pi/2.0)) + 1.0) / 2.0
      if (yewgt .gt. 1.0) yewgt = 1.0
      if (yewgt .lt. 0.0) yewgt = 0.0
      
!.... Want weight to start out higher for south side and decrease as you go north      
      t3(i, j, k) = (1.0 - ywgt) * tz(k) + ywgt * tz_meso(k, 1)
      q3(i, j, k) = (1.0 - ywgt) * qz(k) + ywgt * qz_meso(k, 1)
      v3(i, j, k) = (1.0 - yewgt) * v1d(k) + yewgt * v3(i, j, k)
      
      endif 
      
!
!.... enddo k = 1, nzmeso loop through mesoscale domain prescribing theta and qv
!
    enddo
    enddo
    enddo

 

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_DATA3D_STAT_FRONT END')
!-----------------------------------------------------------------------------



!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM DEFINITIONS (Fortran90 stuff....)

  CONTAINS



!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SNDREAD_INHOM3D:  Reads in a sounding from an input file and interpolates
!                            the sounding to the grid.

   SUBROUTINE SNDREAD_INHOM3D(insndfile,outsoundfile,nzmeso,uzwork,vzwork,tzwork,qzwork)

! Define arrays for reading in the sounding

    implicit none


    integer nzmeso

    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax)
    real uzwork(nzmeso)
    real vzwork(nzmeso)
    real tzwork(nzmeso)
    real qzwork(nzmeso)

!
!.... internal work array
!

    real pzwork(nzmeso)
    
    real z0, tv0, tv1
    real z, an, ugnd

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac
    character(LEN=255)         :: outsoundfile
    character(LEN=100)         :: insndfile


! Error check


!------------------------------------------------------------------------------
! First read in the surface pressure, temperature, and mixing ratio
! Format is z,th,qv,u,v because of stretched grid
!     Sounding input file format: ascii
!     z (meters), theta (K), qv (g/kg),  u (m/s), v (m/s)
!------------------------------------------------------------------------------

    CALL FILE_MESSAGE('READING SOUNDING FROM: '//insndfile)

!
!.... open the input sounding unit
!

    open(unit= 17, file= insndfile, status= 'old', form= 'formatted')
!    open(unit= 17, file= sndfile, status= 'old', form= 'formatted')
    rewind(17)

! Read in sfc pressure, theta, qv

    read(17,*) p000, t000, q000

    write(luno,*) p000, t000, q000

   psfc = p000
   tsfc = t000
   qsfc = q000
   qfac = 1.0
   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qfac = 0.001
   
   qsfc = qfac*qsfc


! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! Read in sounding values

    nsnd = 0
    DO k = 1,nmax
	 read(17,*,iostat=ios) zsnd(k),tsnd(k),qvsnd(k),usnd(k),vsnd(k)
	 qvsnd(k) = qfac*qvsnd(k)
       IF( ios < 0 ) EXIT
	 nsnd = nsnd + 1
    ENDDO

!
!.... close the input sounding unit
!

    close(17)




#ifdef MPI
    if (kzbeg .eq. nzbeg)  zsfc = zsnd(1)
#else
    zsfc = zsnd(1)
#endif

   IF( .not. SET_VARIABLE(gd,'HGT',       zsfc) ) write(6,*) 'INIT_GRID:  Problem setting HGT'

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .le. nsnd ) kze = nsnd-kzbeg
    DO k = kzb,kze
#else
    DO k = 1,nsnd
#endif
      zsnd(k) = zsnd(k) - zsfc
    ENDDO

    zsnd(nsnd) = zsnd(nsnd) + 0.1

! Check to make sure the sounding levels span the computational grid

#ifdef MPI
    if (kzbeg .eq. nzbeg .AND. zsnd(1) .gt. zc(1) ) THEN
#else
    IF( zsnd(1) .gt. zc(1) ) THEN
#endif
	write(luno,*) 'ZMIN of SOUNDING   = ',zsnd(1)
	write(luno,*) 'ZMIN of MODEL GRID = ',zc(1)
	write(luno,*) 'ZMIN OF INPUT SOUNDING IS > ZMIN OF GRID'
	write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
	stop
    ENDIF

#ifdef MPI
    IF( zsnd(nsnd) .lt. zc(nzend-1) ) THEN
#else
    IF( zsnd(nsnd) .lt. zc(nz-1) ) THEN
#endif
	write(luno,*) 'ZMAX of SOUND      = ',zsnd(nsnd)
	write(luno,*) 'ZMAX of MODEL GRID = ',zc(nz-1)
	write(luno,*) 'ZMAX OF INPUT SOUNDING IS < ZMAX OF MODEL GRID'
	write(luno,*) 'Extrapolating sounding with isothermal layer....'
!      write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
!      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
!      stop
    ENDIF

! Now interpolate sounding to grid

    write(luno,*)
    write(luno,*) 'INTERPOLATING SOUNDING TO GRID '

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg
    DO k = kzb,kze
#else
    DO k = 1,nz-1
#endif
      z0 = zc(k)

      IF ( z0 .lt. zsnd(nsnd) ) THEN

!    real uzwork(nzmeso)
!    real vzwork(nzmeso)
!    real tzwork(nzmeso)
!    real qzwork(nzmeso)

        uzwork(k) = zfind(zsnd, usnd,  nsnd, z0)
        vzwork(k) = zfind(zsnd, vsnd,  nsnd, z0)
        tzwork(k) = zfind(zsnd, tsnd,  nsnd, z0)
        qzwork(k) = zfind(zsnd, qvsnd, nsnd, z0)

!        u1d(k) = zfind(zsnd, usnd,  nsnd, z0)
!        v1d(k) = zfind(zsnd, vsnd,  nsnd, z0)
!        tz(k) = zfind(zsnd, tsnd,  nsnd, z0)
!        qz(k) = zfind(zsnd, qvsnd, nsnd, z0)

        ktop = k

        IF ( k .eq. 1 ) THEN
          tv0   = tzwork(1) * (1.0 + 0.601*qzwork(1))
!          tv0   = tz(1) * (1.0 + 0.601*qz(1))
          pzwork(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
        ELSE
          tv1   = tzwork(k) * (1.0 + 0.601*qzwork(k))
!          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pzwork(k) = pzwork(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
          tmptr = pzwork(k)*tzwork(k)
!          tmptr = pz(k)*tz(k)
        ENDIF

      ELSE

        uzwork(k) = uzwork(ktop)
        vzwork(k) = vzwork(ktop)
        qzwork(k) = qzwork(ktop)
        tzwork(k) = tzwork(ktop) * (Exp((z0 - zc(ktop) - zsnd(1)) * g/cp/tmptr ) )

!        u1d(k) = u1d(ktop)
!        v1d(k) = v1d(ktop)
!        qz(k) = qz(ktop)
!        tz(k) = tz(ktop)*(Exp((z0-zc(ktop) - zsnd(1))*g/cp/tmptr ) )

          tv1   = tzwork(k) * (1.0 + 0.601*qzwork(k))
!          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pzwork(k) = pzwork(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1

      ENDIF

!      IF( qz(k) .gt. 0.05 ) qz(k) = qz(k) / 1000.

    ENDDO


! Write sounding out to file

    IF ( ny .le. 2 ) THEN
      v1d(:) = 0.0
    ELSEIF ( nx .le. 2 ) THEN
      u1d(:) = 0.0
    ENDIF
    

#ifdef MPI
  if(my_rank == 0) then
#endif
    open(21,file=outsoundfile, status='unknown')
    rewind(21)
!    write(21,*) 'MODEL SND'
    write(21,'(1x,2(f9.2,2x),f9.4)') psfc, t000, qsfc*1000.

    DO k = 1,nz-1
!     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), tzwork(k),1000.*qzwork(k),uzwork(k),vzwork(k),pzwork(k)
!     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),psfc
!     write(21,'(1x, 6(f10.2,2x))') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),p000

    ENDDO
    close(21)

#ifdef MPI
  endif
#endif
   END SUBROUTINE SNDREAD_INHOM3D



 END SUBROUTINE INIT_DATA3D_DEF_FRONT








!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE HYDRO_BALANCE  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to user defined background state
!
!-----------------------------------------------------------------------

 SUBROUTINE HYDRO_BALANCE(gd) 

! SUBROUTINE HYDRO_BALANCE(gd) !,                                       &
!                               sndtype, sndfile,                         &
!                               wtype, Us, Uz, psfc, tsfc, qsfc, rhmax, &
!                               shape, dudz0, dudz1, dudz2, z0, z1, z2)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

! Local vars


  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer nxin,nyin,nzin,xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes, nltimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


!  real lag_thet(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_qrat(nxmeso+1,nymeso+1,nzmeso+1) 
!  real lag_u(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_v(nxmeso+1,nymeso+1,nzmeso+1)


  real dxin,dyin,dzin,dx,dy,dz
  real grid_xin(nxmeso+1)
  real grid_yin(nymeso+1)
  real grid_zin(nzmeso+1)
  real f1,f2,f3,f4,f5,f6,f7,f8
  real grid_x(1000),grid_y(1000),grid_z(1000)  
  real tv(25),ufl(11),vfl(11)
  real gxin,gyin,gx,gy
  real nrint,nlint
  
  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)

  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)
  real, pointer :: p3(:,:,:)

  real  :: pp
    
  real  :: u3t(nxend,nyend,nzend)
  real  :: v3t(nxend,nyend,nzend)
  real  :: t3t(nxend,nyend,nzend)
  real  :: q3t(nxend,nyend,nzend)
  real  :: den(nzend)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
!  real, allocatable, dimension(:) ::  den
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  






!
!.... begin executable code
!




!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)
  CALL GET_VARIABLE(gd,'PI',p3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    outsoundfile = attr%str(ibeg:iend)//'.output.sound'


!
!.... CLZ(3-6-13): mesoscale work arrays were defined in INIT_MODULE,
!                    and allocated in INIT_*_DEF_*
!


!  allocate( cz(nz) )
!  allocate(den(nz))


!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'HYDRO_BALANCE BEGIN')
  write(luno,*) 


!#ifdef MPI

    write(luno,*) 'kzbeg=',kzbeg,' kzend=',kzend

    do k = kzbeg, nz - 1
!.... CLZ (3/4/13): force upper limit of nz - 1
!    do k = kzbeg, kzend

    if (k .eq. kzend) then
    write(luno,*) 'k=',k,' pz=',pz(k),' tz=',tz(k),' qz=',qz(k)
    endif

    pp = (pz(k)**(cp/rd))*p00
    den(k) = pp/(rd*pz(k)*tz(k)*(1 + 0.61*qz(k)))
    write(luno,*) 'k=',k,' pp=',pp,' den(k)=',den(k)

    enddo

      write(luno,*)
      write(luno,*)

    do i = ixbeg, ixend
     do j = jybeg, jyend



      do k = kzend - 2, kzbeg, -1

!      if(i .eq. ixbeg .and. j .eq. jybeg) then
!      write(luno,*)'about to compute p3 at i,j,k=',i,j,k
!      endif
        
!     p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)  &
!      - (g*0.5*(den(k-kzbeg+1)+den(k-kzbeg+2)))*   &             
!     ((0.5*(t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)+t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)) &
!       -0.5*(tz(k-kzbeg+1)+tz(k-kzbeg+2)))/(0.5*(tz(k-kzbeg+1)+tz(k-kzbeg+2))))*(1/dzc(k))            

     p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1) = p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)  &
      - (g*0.5*(den(k-kzbeg+1)+den(k-kzbeg+2)))*     &             
     ((0.5*(t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)*       &
     (1 + 0.61*q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1))    &
     +t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)*             &
     (1 + 0.61*q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)))   &
       -0.5*(tz(k-kzbeg+1)*(1 + 0.61*qz(k-kzbeg+1))+ &
       tz(k-kzbeg+2)*(1 + 0.61*qz(k-kzbeg+2))))/     &
       (0.5*(tz(k-kzbeg+1)*(1 + 0.61*qz(k-kzbeg+1))+ &
       tz(k-kzbeg+2)*(1 + 0.61*qz(k-kzbeg+2)))))*(1/dzc(k)) 

      if(i .eq. 2 .and. j .eq. 2) then
      write(luno,*) 'k=',k,' zc=',zc(k),' p3=',p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
      write(luno,*) 'tzk=',tz(k-kzbeg+1),' tzkp1=',tz(k-kzbeg+2)
      write(luno,*) 't3k=',t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
      write(luno,*)
      endif


      enddo
     enddo
    enddo


    do i = ixbeg, ixend
     do j = jybeg, jyend
      do k = kzbeg, kzend - 1

      if(i.eq.2.and.j.eq.2) then
      print*,p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1),'pres',k
      print*,tz(k-kzbeg+1),tz(k-kzbeg+2),t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
      endif


      p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=     &
      p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)/     &
      (den(k)*cp*tz(k-kzbeg+1)*(1 + 0.61*qz(k-kzbeg+1)))   
      

!      if(i.eq.100.and.j.eq.100) then
!      print*,p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1),'pi',k
!      end if



!      p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=     &
!      p3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)/     &
!      (den(k)*rd*0.5*(t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)+    &
!      t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+2)))   
      
      enddo
     enddo
    enddo

      
!#else
!
!      do k=1,nz
!
!      pp=(pz(k)**(cp/rd))*p00
!      den(k)=pp/(rd*pz(k)*tz(k)*(1 + .61*qz(k)))
!    
!      end do
!   
!      do i=1,nx
!      do j=1,ny
!
!      p3(i,j,nz)=0.00
!      p3(i,j,nz-1)=0.00
!
!      do k=nz-2,1
!
!      p3(i,j,k)=p3(i,j,k+1) - ((g*0.5*(den(k)+den(k+1))*   &
!      (0.5*(t3(i,j,k)*(1 + 0.61*q3(i,j,k))+                &
!      t3(i,j,k+1)*(1 + 0.61*q3(i,j,k+1)))                  &
!       - 0.5*(tz(k)*(1 + 0.61*qz(k))-tz(k+1)*              &
!       (1 + 0.61*qz(k+1))))*                               & 
!      (1/dzc(k)))/(den(k)*rd*0.5*(t3(i,j,k)*               &
!      (1 + 0.61*q3(i,j,k))+t3(i,j,k+1)*(1 + 0.61*q3(i,j,k+1)))))
!      
!      enddo
!      enddo
!      enddo
! 
!#endif
      close(14)

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('HYDRO_BALANCE END')
!-----------------------------------------------------------------------------

 END SUBROUTINE HYDRO_BALANCE












!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_USER_DEF  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to user defined background state
!
!-----------------------------------------------------------------------
 SUBROUTINE INIT_USER_DEF(gd) !,                                       &
!                               sndtype, sndfile,                         &
!                               wtype, Us, Uz, psfc, tsfc, qsfc, rhmax, &
!                               shape, dudz0, dudz1, dudz2, z0, z1, z2)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE INIT_MODULE

  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
   include "mpif.h"
#endif

! Passed variables

  TYPE(GRID) :: gd

!  integer sndtype
!  character(LEN = 100) :: sndfile   
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars

  integer i, j, k, m, l
  integer n, ns, s, iadiabat
  integer lccn, lscpi, lscni
  integer lmupzc, lmunzc, lmupze, lmunze
  integer nxin,nyin,nzin,xmx,ymx
  integer ix1,ix2,iy1,iy2,iz1,iz2
  integer nrtimes, nltimes
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile   
  integer                    :: ibeg,iend


!  real lag_thet(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_qrat(nxmeso+1,nymeso+1,nzmeso+1) 
!  real lag_u(nxmeso+1,nymeso+1,nzmeso+1)
!  real lag_v(nxmeso+1,nymeso+1,nzmeso+1)


  real dxin,dyin,dzin,dx,dy,dz
  real grid_xin(nxmeso+1),grid_yin(nymeso+1)
  real grid_zin(nzmeso+1)
  real f1,f2,f3,f4,f5,f6,f7,f8
  real grid_x(1000),grid_y(1000),grid_z(1000)  
  real tv(25),ufl(11),vfl(11)
  real gxin,gyin,gx,gy
  real nrint,nlint

  character(LEN=70)  :: infile(100)
  character(LEN=70)  :: inrad(100)


  integer, pointer :: nx, ny, nz
!  real, pointer :: dx, dy, dz
  real, pointer :: zc(:)
  real, pointer :: ze(:)
  real, pointer :: dzc(:)
  real, pointer :: dze(:)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)
  real, pointer :: kz(:)
  real, pointer :: zcdx(:,:)
  real, pointer :: zedx(:,:)


  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)

  real, pointer :: ccn1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
!  real, pointer :: k3(:,:,:)

  real, allocatable :: cz(:)
!  real, allocatable :: grid_x(:)
!  real, allocatable :: grid_y(:)
!  real, allocatable :: grid_z(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den
  
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  






!
!.... begin executable code
!









!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd, 'DX', dx)
  CALL GET_VARIABLE(gd, 'DY', dy)
  CALL GET_VARIABLE(gd, 'DZ', dz)

!  CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!  CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
  CALL GET_VARIABLE(gd,'ZC',    zc)
  CALL GET_VARIABLE(gd,'ZE',    ze)
  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'PIINIT',pz)
  CALL GET_VARIABLE(gd,'KMINIT',kz)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'DZC',dzc)
  CALL GET_VARIABLE(gd,'DZE',dze)

!  CALL GET_VARIABLE(gd,'VS',v2s)
!  CALL GET_VARIABLE(gd,'VN',v2n)
!  CALL GET_VARIABLE(gd,'UW',u2w)
!  CALL GET_VARIABLE(gd,'UE',u2e)

  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    outsoundfile = attr%str(ibeg:iend)//'.output.sound'






!
!.... allocate the mesoscale domain array abmeso and work arrays
!

   allocate ( lag_thet(62, 62, 26) )
   allocate ( lag_qrat(62, 62, 26) )
   allocate ( lag_u(62, 62, 26) )
   allocate ( lag_v(62, 62, 26) )




!  allocate( cz(nz) )
  allocate(den(nz))



!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_USER_DEF BEGIN')

      open(14, file='user_def.int',status='unknown')
      read(14,4022)(((lag_thet(m,l,k),m=1,61),l=1,61),k=1,25)
      read(14,4022)(((lag_qrat(m,l,k),m=1,61),l=1,61),k=1,25)
      read(14,4022)(((lag_u(m,l,k),m=1,61),l=1,61),k=1,25)
      read(14,4022)(((lag_v(m,l,k),m=1,61),l=1,61),k=1,25)
      close(14)

      do i=1,62
      do j=1,62
      do k=1,26
      lag_thet(62,j,k)=lag_thet(61,j,k)
      lag_thet(i,62,k)=lag_thet(i,61,k)
      lag_thet(i,j,26)=lag_thet(i,j,25)
      lag_qrat(62,j,k)=lag_qrat(61,j,k)
      lag_qrat(i,62,k)=lag_qrat(i,61,k)
      lag_qrat(i,j,26)=lag_qrat(i,j,25)
      lag_u(62,j,k)=lag_u(61,j,k)
      lag_u(i,62,k)=lag_u(i,61,k)
      lag_u(i,j,26)=lag_u(i,j,25)
      lag_v(62,j,k)=lag_v(61,j,k)
      lag_v(i,62,k)=lag_v(i,61,k)
      lag_v(i,j,26)=lag_v(i,j,25)
      end do
      end do
      end do

 4022 format (61(f10.4,1x))

!      do i=1,25
!      print*,lag_thet(2,2,i),lag_u(2,2,i),lag_qrat(2,2,i),'lags'
!      end do


!  DO k = 1,nz
!   DO j = 1,ny
!    DO i = 1,nx
!
!     u3(i,j,k)   = lag_u(i,j,k)      ! U
!     v3(i,j,k)   = lag_v(i,j,k)      ! V
!     t3(i,j,k)   = lag_thet(i,j,k)      ! TH
!     q3(i,j,k)   = lag_qrat(i,j,k)      ! QV
!
!    ENDDO
!   ENDDO
!  ENDDO


!......fill in around xpol
!
!      do k=1,11
!      ufl(k)=(lag_u(43,45,k)+lag_u(51,45,k))/2.0
!      vfl(k)=(lag_v(47,41,k)+lag_v(47,49,k))/2.0
!      end do
!
!      do i=45,49
!      do j=43,47
!      do k=1,11
!      lag_u(i,j,k)=ufl(k)
!      lag_v(i,j,k)=vfl(k)
!      end do
!      end do
!      end do
!
!......fill in around Dow 3
!
!      do k=1,11
!      ufl(k)=(lag_u(53,19,k)+lag_u(61,19,k))/2.0
!      vfl(k)=(lag_v(57,15,k)+lag_v(57,23,k))/2.0
!      end do
!
!      do i=55,59
!      do j=17,21
!      do k=1,11
!      lag_u(i,j,k)=ufl(k)
!      lag_v(i,j,k)=vfl(k)
!      end do
!      end do
!      end do
!
!......fill in around SR1
!
!      do k=1,11
!      ufl(k)=(lag_u(21,5,k)+lag_u(29,5,k))/2.0
!      vfl(k)=(lag_v(25,1,k)+lag_v(25,9,k))/2.0
!      end do
!
!      do i=23,27
!      do j=3,7
!      do k=1,11
!      lag_u(i,j,k)=ufl(k)
!      lag_v(i,j,k)=vfl(k)
!      end do
!      end do
!      end do
!
!
!......additional fill
!
!      do k=1,11
!      ufl(k)=(lag_u(46,45,k)+lag_u(54,45,k))/2.0
!      vfl(k)=(lag_v(50,41,k)+lag_v(50,49,k))/2.0
!      end do
!
!      do i=48,52
!      do j=43,47
!      do k=1,11
!      lag_u(i,j,k)=ufl(k)
!      lag_v(i,j,k)=vfl(k)
!      end do
!      end do
!      end do



!--------------------------------------------------
!
! Tri-linear interpolation of input data to model grid
! msb 10/30/07
!
! dxin,dyin,dzin => input grid spacing
! nxin,nyin,nzin => input grid dimensions
!

!  dxin=500.0
!  dyin=500.0
!  dzin=250.0
!  nxin=61
!  nyin=61
!  nzin=25

      open(13,file='tdlbc.in',status='old')
      read(13,3101) nrtimes, nrint
 3101 format(1x,i2,1x,f5.0)
      do i=1,nrtimes
      read(13,3012) inrad(i)
      end do
      read(13,3101) nltimes, nlint
 3012 format(1x,a70)
      do i=1,nltimes
      read(13,3012) infile(i)
      end do
      read(13,3013) dxin,dyin,dzin
      read(13,3014) nxin,nyin,nzin
 3013 format(3(1x,f5.0))
 3014 format(3(1x,i3))
      close(13)


!  dx=((real(nxin)-1)/(real(nx)-1))*dxin
!  dy=((real(nyin)-1)/(real(ny)-1))*dyin
!  dz=((real(nzin)-1)/(real(nz)-1))*dzin

  do i=1,nxin+1
  grid_xin(i)=(real(i)-1)*dxin
  enddo

  do i=1,nyin+1
  grid_yin(i)=(real(i)-1)*dyin
  enddo

  do i=1,nzin+1
  grid_zin(i)=(real(i)-1)*dzin
  enddo

#ifdef MPI
  do i=nxbeg,nxend
#else  
  do i=1,nx
#endif 
  grid_x(i)=(real(i)-1)*dx
  enddo

#ifdef MPI
  do i=nybeg,nyend
#else  
  do i=1,ny
#endif 
  grid_y(i)=(real(i)-1)*dy
  enddo

!  do i=1,ny
!  grid_y(i)=(real(i)-1)*dy
!  enddo



!...msb 5/28/08 change for stretched vertical grid




#ifdef MPI
  grid_z(nzbeg)=0.0  
  do i=nzbeg+1,nzend
#else  
  grid_z(1)=0.0
  do i=2,nz
#endif 
  grid_z(i)=grid_z(i-1)+(1/dzc(i-1))
  enddo



   gxin=(real(nxin)-1)*dxin
   gyin=(real(nyin)-1)*dyin

#ifdef MPI   
   gx=(real(nxend)-1)*dx
   gy=(real(nyend)-1)*dy
   xmx=(gxin/gx)*(nxend-1) + 1
   ymx=(gyin/gy)*(nyend-1) + 1
#else   
   gx=(real(nx)-1)*dx
   gy=(real(ny)-1)*dy
   xmx=(gxin/gx)*(nx-1) + 1
   ymx=(gyin/gy)*(ny-1) + 1   
#endif



#ifdef MPI
   do j=jybeg,jyend
    do i=ixbeg,ixend
     do k=kzbeg,kzend
#else     
   do j=1,ny
    do i=1,nx
     do k=1,nz
#endif

#ifdef MPI
  ix1= (real(nxin)-1)/(real(nxend)-1)*(gx/gxin)*(i-1) + 1
  
  if(ix1.gt.nxin) then
  ix1=nxin
  end if
  ix2=ix1+1

  iy1= (real(nyin)-1)/(real(nyend)-1)*(gy/gyin)*(j-1) + 1
  if(iy1.gt.nyin) then
  iy1=nyin
  end if
  iy2= iy1+1

#else 
    
  ix1= (real(nxin)-1)/(real(nx)-1)*(gx/gxin)*(i-1) + 1
  
  if(ix1.gt.nxin) then
  ix1=nxin
  end if
  ix2=ix1+1
  
  iy1= (real(nyin)-1)/(real(ny)-1)*(gy/gyin)*(j-1) + 1
  if(iy1.gt.nyin) then
  iy1=nyin
  end if
  iy2= iy1+1
  
#endif

    

!.....msb 5/28/08 change for stretched vertical grid
   iz1=grid_z(k)/dzin + 1
  iz2= iz1 +1
!  print*,iz1,iz2,grid_z(k),dzin
!  iz1= (real(nzin)-1)/(real(nz)-1)*(k-1) + 1  
!  iz2= iz1 +1
!    print*,dzin,1/dzc(k),grid_z(k)/dzin,iz1,iz2,k

!     print*,i,j,k,ix1,ix2,iy1,iy2,iz1,iz2,'n' 

     f1=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f2=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f3=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f4=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))

     f5=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f6=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f7=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f8=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))

!     print*,f1,f2,f3,f4,f5,f6,f7,f8,'h'

#ifdef MPI

     u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_u(ix1,iy1,iz1)+f2*lag_u(ix2,iy1,iz1)+ &
         f3*lag_u(ix1,iy2,iz1)+f4*lag_u(ix2,iy2,iz1)+ &
     f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

     v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+&
         f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
     f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

     t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+&
         f3*lag_thet(ix1,iy2,iz1)+f4*lag_thet(ix2,iy2,iz1)+ &
     f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+f8*lag_thet(ix2,iy2,iz2)

     q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+&
         f3*lag_qrat(ix1,iy2,iz1)+f4*lag_qrat(ix2,iy2,iz1)+ &
     f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#else

     u3(i,j,k)=f1*lag_u(ix1,iy1,iz1)+f2*lag_u(ix2,iy1,iz1)+f3*lag_u(ix1,iy2,iz1)+f4*lag_u(ix2,iy2,iz1)+ &
     f5*lag_u(ix1,iy1,iz2)+f6*lag_u(ix2,iy1,iz2)+f7*lag_u(ix1,iy2,iz2)+f8*lag_u(ix2,iy2,iz2)

     v3(i,j,k)=f1*lag_v(ix1,iy1,iz1)+f2*lag_v(ix2,iy1,iz1)+f3*lag_v(ix1,iy2,iz1)+f4*lag_v(ix2,iy2,iz1)+ &
     f5*lag_v(ix1,iy1,iz2)+f6*lag_v(ix2,iy1,iz2)+f7*lag_v(ix1,iy2,iz2)+f8*lag_v(ix2,iy2,iz2)

     t3(i,j,k)=f1*lag_thet(ix1,iy1,iz1)+f2*lag_thet(ix2,iy1,iz1)+f3*lag_thet(ix1,iy2,iz1)+f4*lag_thet(ix2,iy2,iz1)+ &
     f5*lag_thet(ix1,iy1,iz2)+f6*lag_thet(ix2,iy1,iz2)+f7*lag_thet(ix1,iy2,iz2)+f8*lag_thet(ix2,iy2,iz2)

     q3(i,j,k)=f1*lag_qrat(ix1,iy1,iz1)+f2*lag_qrat(ix2,iy1,iz1)+f3*lag_qrat(ix1,iy2,iz1)+f4*lag_qrat(ix2,iy2,iz1)+ &
     f5*lag_qrat(ix1,iy1,iz2)+f6*lag_qrat(ix2,iy1,iz2)+f7*lag_qrat(ix1,iy2,iz2)+f8*lag_qrat(ix2,iy2,iz2)

#endif


!       if(k.eq.1.and.my_rank.eq.1) then
!       print*, q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1), t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)
!       print*,i-ixbeg+1,j-jybeg+1,k-kzbeg+1,i,j,k,ixbeg,jybeg,kzbeg
!       end if

    enddo
   enddo
  enddo

!---------------------------------------------------
!
! End tri-linear interpolation


! msb 02/26/10
! Extend grid

#ifdef MPI
  do i=ixbeg,ixend
  do j=jybeg,jyend
  do k=kzbeg,kzend
#else
  do i=1,nx
  do j=1,ny
  do k=1,nz
#endif

#ifdef MPI

  if(i.gt.xmx.and.j.le.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,j,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,j,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,j,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(1,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(1,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(1,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(1,1,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=u3(i,1,k)
  v3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=v3(i,1,k)
  t3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=t3(i,1,k)
  q3(i-ixbeg+1,j-jybeg+1,k-kzbeg+1)=q3(i,1,k)
  end if
 
#else
 
  if(i.gt.xmx.and.j.le.ymx) then
  u3(i,j,k)=u3(xmx,j,k)
  v3(i,j,k)=v3(xmx,j,k)
  t3(i,j,k)=t3(xmx,j,k)
  q3(i,j,k)=q3(xmx,j,k)
!  print*,xmx,ymx,i,j,k,u3(xmx,j,k)
  else if(i.gt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(xmx,ymx,k)
  v3(i,j,k)=v3(xmx,ymx,k)
  t3(i,j,k)=t3(xmx,ymx,k)
  q3(i,j,k)=q3(xmx,ymx,k)
  else if(i.lt.xmx.and.j.gt.ymx) then
  u3(i,j,k)=u3(i,ymx,k)
  v3(i,j,k)=v3(i,ymx,k)
  t3(i,j,k)=t3(i,ymx,k)
  q3(i,j,k)=q3(i,ymx,k)
  end if

#endif

  end do
  end do
  end do



!---------------------------
!
! Option to eliminate superadiabatic layer at the boundaries
!

    iadiabat=0

    if(iadiabat.eq.1) then

   do j=1,ny
    do i=1,nx
    tv(25)=t3(i,j,25)*(1 + 0.61*q3(i,j,25))
     do k=nz-1,1,-1

      tv(k)=t3(i,j,k)*(1 + 0.61*q3(i,j,k))

      if((i.eq.1).or.(i.eq.nx).or.(j.eq.1).or.(j.eq.ny)) then
      if(tv(k).gt.tv(k+1)) then

      tv(k)=tv(k+1)
      t3(i,j,k)=tv(k)/(1 + 0.61*q3(i,j,k))
!      print*,t3(i,j,k),tv(k),q3(i,j,k),i,j,k

      end if
      end if

     end do
    end do
   end do

   end if



!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_USER_DEF END')
!-----------------------------------------------------------------------------

 END SUBROUTINE INIT_USER_DEF

















!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE INIT_BACKGROUND_1D  <<<<<<<<<<<<<<<<<<<<<
!
!-------------------------------------------------------------------------------
!
! Initializes the domain to XY homogeneous background state
!
!-----------------------------------------------------------------------
 SUBROUTINE INIT_BACKGROUND_1D(gd,prtbase) !,                                        &
!                               sndtype, sndfile,                          &
!                               wtype, Us, Uz, psfc, tsfc, qsfc, rhmax,  &
!                               shape, dudz0, dudz1, dudz2, z0, z1, z2)


  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE MICRO_MODULE
  USE ELEC_MODULE
  USE INIT_MODULE
  USE COMMASMPI_MODULE

  implicit none

! Passed variables

  TYPE(GRID) :: gd
  logical :: prtbase ! whether to print out the base state or not

!  integer sndtype
!  character(LEN = 100) :: sndfile
!  integer shape
!  real    dudz0, dudz1, dudz2
!  real    z0,    z1,    z2
!  integer wtype
!  real    Us
!  real    Uz
!  real    psfc
!  real    tsfc
!  real    qsfc
!  real    rhmax

! Local vars

  integer i, j, k
  integer n, ns, s
  integer lcin, lccn, lscpi, lscni, lscpli
  integer lccnuf, lcn_ac, lcn_nu, lcn_co
  integer lmupzc, lmunzc, lmupze, lmunze
  character(LEN=name_length)   names
  real qvs, pres, zfac
  TYPE(ATTRIBUTE), pointer   :: attr, microphys
  character(LEN=255)         :: outsoundfile
  integer                    :: ibeg,iend
  
  integer, pointer :: nx, ny, nz
  real, pointer :: zc1(:)
  real, pointer :: ze1(:)  
  real, pointer :: dzc1(:)
  real, pointer :: dze1(:)
  real, pointer :: u1d1(:)
  real, pointer :: v1d1(:)
  real, pointer :: tz1(:)
  real, pointer :: qz1(:)
  real, pointer :: kz1(:)

  real, pointer :: u1d1b(:)
  real, pointer :: v1d1b(:)
  real, pointer :: tz1b(:)
  real, pointer :: qz1b(:)
  real, pointer :: pz1b(:)

  real, pointer :: mupzc(:)
  real, pointer :: munzc(:)
  real, pointer :: mupze(:)
  real, pointer :: munze(:)
  
  real, pointer :: ccn1d(:)
  real, pointer :: ccnuf1d(:)
  real, pointer :: nion1d(:)
  real, pointer :: pion1d(:)
  real, pointer :: pz1(:)
  real, pointer :: v2s(:,:)
  real, pointer :: v2n(:,:)
  real, pointer :: u2e(:,:)
  real, pointer :: u2w(:,:)
  real, pointer :: u3(:,:,:)
  real, pointer :: v3(:,:,:)
  real, pointer :: w3(:,:,:)
  real, pointer :: t3(:,:,:)
  real, pointer :: q3(:,:,:)
  real, pointer :: r3(:,:,:)
!  real, pointer :: k3(:,:,:)

  real, allocatable :: cz(:)
  real, pointer :: c3(:,:,:)
  real, pointer :: c4(:,:,:)
  real, allocatable, dimension(:) ::  den
  
  real, allocatable, dimension(:) :: u1d,v1d,tz,qz,pz,kz,  &
!                   nion1d,pion1d,         &
                   zc,ze,dzc,dze,  &
                   temp,thetav,pinit,rinit,rhw,rhi,dewpt,riinit,rivinit


  real :: ps,td,bb,bbv,bsh
  real :: dz
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

  integer :: ixb, jyb, kzb
  integer :: ixe, jye, kze

!-----------------------------------------------------------------------
! Reference the variables

  CALL GET_ATTRIBUTE(gd, 'NX', nx)
  CALL GET_ATTRIBUTE(gd, 'NY', ny)
  CALL GET_ATTRIBUTE(gd, 'NZ', nz)

  CALL GET_VARIABLE(gd,'ZC',    zc1)
  CALL GET_VARIABLE(gd,'ZE',    ze1)
  CALL GET_VARIABLE(gd,'UINIT', u1d1)
  CALL GET_VARIABLE(gd,'VINIT', v1d1)
  CALL GET_VARIABLE(gd,'THINIT',tz1)
  CALL GET_VARIABLE(gd,'QVINIT',qz1)
  CALL GET_VARIABLE(gd,'PIINIT',pz1)
  CALL GET_VARIABLE(gd,'KMINIT',kz1)
  CALL GET_VARIABLE(gd,'CCNINIT',ccn1d)
  CALL GET_VARIABLE(gd,'CCNUFINIT',ccnuf1d)
  CALL GET_VARIABLE(gd,'CNIONINIT',nion1d)
  CALL GET_VARIABLE(gd,'CPIONINIT',pion1d)

  CALL GET_VARIABLE(gd,'UINIT0', u1d1b)
  CALL GET_VARIABLE(gd,'VINIT0', v1d1b)
  CALL GET_VARIABLE(gd,'THINIT0',tz1b)
  CALL GET_VARIABLE(gd,'QVINIT0',qz1b)
  CALL GET_VARIABLE(gd,'PIINIT0',pz1b)

  CALL GET_VARIABLE(gd,'DZC',dzc1)
  CALL GET_VARIABLE(gd,'DZE',dze1)

!mpi!  CALL GET_VARIABLE(gd,'VS',v2s)
!mpi!  CALL GET_VARIABLE(gd,'VN',v2n)
!mpi!  CALL GET_VARIABLE(gd,'UW',u2w)
!mpi!  CALL GET_VARIABLE(gd,'UE',u2e)
 
  CALL GET_VARIABLE(gd,'U', u3)
  CALL GET_VARIABLE(gd,'V', v3)
  CALL GET_VARIABLE(gd,'W', w3)
!  CALL GET_VARIABLE(gd,'KM',k3)
  CALL GET_VARIABLE(gd,'TH',t3)
  CALL GET_VARIABLE(gd,'QV',q3)
  CALL GET_VARIABLE(gd,'RHO',r3)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)

!
!.... CLZ (3/16/13): customize output sounding file name based on inhom
!

    if (inhom .ne. 4) then
     outsoundfile = attr%str(ibeg:iend)//'.output.sound'
    elseif (inhom .eq. 4) then
     outsoundfile = attr%str(ibeg:iend)//'.mdl_base.snd'
    endif

!  allocate( cz(nz) )
  allocate(den(-ng+1:nzend+ng))


!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INIT_BACKGROUND_1D BEGIN')

!-----------------------------------------------------------------------

  write(luno,*)
  write(luno,*) 'SOUNDING TYPE        = ',sndtype
  write(luno,*) 'WK SND WIND          = ',wtype
  write(luno,*) 'WK SND US            = ',Us
  write(luno,*) 'WK SND Zmax          = ',Uz
  write(luno,*) 'OUTPUT SOUNDING FILE = ',outsoundfile
  write(luno,*)

!--------------------------------------------------------------------------
! INITIALIZE SOUNDING

  allocate( u1d(-ng+1:nzend+ng),v1d(-ng+1:nzend+ng),tz(-ng+1:nzend+ng),  &
             qz(-ng+1:nzend+ng),pz(-ng+1:nzend+ng),kz(-ng+1:nzend+ng) )
  
  allocate( zc(-ng+1:nzend+ng),ze(-ng+1:nzend+ng),dzc(-ng+1:nzend+ng),dze(-ng+1:nzend+ng) )

!  allocate( mupzc(-ng+1:nzend+ng),munzc(-ng+1:nzend+ng),mupze(-ng+1:nzend+ng),munze(-ng+1:nzend+ng) )
!  allocate( nion1d(-ng+1:nzend+ng), pion1d(-ng+1:nzend+ng) )
  
  
  DO k = 1,nzend
    zc(k) = z1dinit(k,1)
    ze(k) = z1dinit(k,2)
    dzc(k) = z1dinit(k,3)
    dze(k) = z1dinit(k,4)
  ENDDO
  
  u1d(:) = 0.0
  v1d(:) = 0.0
  kz (:) = 0.0
  qz (:) = 0.0
  
  IF( sndtype .eq. -1 ) CALL SNDM1()  ! IF SND_TYPE = 1, initialize with constant temperature and density environment with near-saturation (liquid or ice)


  IF( sndtype .eq. 0 ) CALL SND0X(gd,tz,pz,qz,dzc,dze,nzend,ng,psfc,tsfc,qsfc)   ! IF SND_TYPE = 0, initialize with adiabatic, dry, no-flow environment


  IF( sndtype .eq. 1 ) CALL SND1X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   outsoundfile,wkumax1,wkumax2,wkvmax,wkh1,wkh2)   ! IF SND_TYPE = 1, initialize with WK sounding + wind profile


!  IF( sndtype .eq. 2 ) CALL SND2()   ! IF SND_TYPE = 2, initialize with sounding in sndfile
  IF( sndtype .eq. 2 ) CALL SND2X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmaxin,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   outsoundfile,sndfile)   ! IF SND_TYPE = 1, initialize with WK sounding + wind profile


  IF( sndtype .eq. 3 ) CALL SND3()   ! IF SND_TYPE = 3, user supplied sounding


  IF( sndtype .eq. 4 ) CALL SND4()   ! IF SND_TYPE = 4, user supplied sounding

!
! eventually have only rank0 call snd and then broadcast the sounding arrays (only need one to read the input
! sounding that way).  (Or have rank0 read the sounding in snd2 and broadcast contents.)
!
   kzb = -ng+1
   kze = ktile+ng
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg+1
   DO k = kzb, kze
     u1d1(k) = u1d(k+kzbeg-1)
     v1d1(k) = v1d(k+kzbeg-1)
     tz1(k)  = tz(k+kzbeg-1)
     qz1(k)  = qz(k+kzbeg-1)
     pz1(k)  = pz(k+kzbeg-1)
     kz1(k)  = kz(k+kzbeg-1)
     u1d1b(k) = u1d(k+kzbeg-1)
     v1d1b(k) = v1d(k+kzbeg-1)
     tz1b(k)  = tz(k+kzbeg-1)
     qz1b(k)  = qz(k+kzbeg-1)
     pz1b(k)  = pz(k+kzbeg-1)
   ENDDO

!--------------------------------------------------------------------------
! 3D arrays

! Initialize to base states (not valid for all variables)

!#ifdef MPI
  ixb = -ng+1
  ixe = itile+ng
  if (ixbeg .eq. nxbeg) ixb = 1
  if (ixend .eq. nxend) ixe = ixend-ixbeg+1

  jyb = -ng+1
  jye = jtile+ng
  if (jybeg .eq. nybeg) jyb = 1
  if (jyend .eq. nyend) jye = jyend-jybeg+1

  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .eq. nzbeg) kzb = 1
  if (kzend .eq. nzend) kze = kzend-kzbeg+1
  
  DO k = kzb,kze
!    write(0,*) 'k,u,v,t,q = ',k,u1d1(k),v1d1(k),tz1(k),qz1(k),pz(k)
   DO j = jyb,jye
    DO i = ixb,ixe
!#else  
!  DO k = 1,nz
!   DO j = 1,ny
!    DO i = 1,nx
!#endif

     u3(i,j,k)   = u1d1(k)      ! U
     v3(i,j,k)   = v1d1(k)      ! V
!     k3(i,j,k)   = kz1(k)      ! KM
     t3(i,j,k)   = tz1(k)      ! TH
     q3(i,j,k)   = qz1(k)      ! QV
     w3(i,j,k)   = 0.0
     IF ( k < kze ) THEN
       r3(i,j,k)   = 1.0e5*pz(k)**cvr/(rd*tz1(k) * (1.0 + 0.61*qz1(k)))
     ELSE
       r3(i,j,k) = r3(i,j,k-1)
     ENDIF

    ENDDO
   ENDDO

   gd%var(GET_VARIABLE_INDEX(gd,  'U') )%base1d(k) = u1d1(k)
   gd%var(GET_VARIABLE_INDEX(gd,  'V') )%base1d(k) = v1d1(k)
   gd%var(GET_VARIABLE_INDEX(gd, 'QV') )%base1d(k) = qz1(k)
   gd%var(GET_VARIABLE_INDEX(gd, 'TH') )%base1d(k) = tz1(k)
  IF ( tke_type == 1 ) THEN
   gd%var(GET_VARIABLE_INDEX(gd, 'KM') )%base1d(k) = kz1(k)
  ELSE
   gd%var(GET_VARIABLE_INDEX(gd, 'TKE') )%base1d(k) = kz1(k)
  ENDIF
  ENDDO

! Init outflow boundary velocities

#ifdef MPI
  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .eq. nzbeg) kzb = 1
  if (kzend .eq. nzend) kze = kzend-kzbeg+1

  DO k = kzb,kze
#else
  DO k = 1,nz
#endif
!mpi!   u2e(:,k) = u1d(k)
!mpi!   u2w(:,k) = u1d(k)
!mpi!   v2s(:,k) = v1d(k)
!mpi!   v2n(:,k) = v1d(k)
  ENDDO

! set up air density:

!#ifdef MPI
!  kzb = -ng+1
!  kze = ktile+ng
!  if (kzbeg .eq. nzbeg) kzb = 1
!  if (kzend .eq. nzend) kze = kzend-kzbeg
!
!  DO k = kzb,kze
!#else
  DO k = 1,nzend-1
!#endif
    den(k) = pz(k)**cvr * p00 / (rd * tz(k))
  ENDDO

!
! Initialize CCN (and ion concentrations later....)
!
  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  s      = GET_VARIABLE_INDEX(gd, 'TH')

  lcn_ac = 0
  lcn_nu = 0
  lcn_co = 0
  lccn = 0
  lccnuf = 0
  lcin = 0
  lscpi = 0
  lscni = 0
  lmupzc = 0
  lmunzc = 0
  lmupze = 0
  lmunze = 0

  CALL GET_ATTRIBUTE (gd, 'MICROPHYS', microphys)

  DO n = s,s + ns - 1
!    print*, 'var name = ', gd%var(n)%name
    names = gd%var(n)%name
    IF ( names .eq. 'CCCN' ) THEN
      lccn = n
      CALL GET_VARIABLE(gd,'CCCN',c3)
      CALL SETCCN()
    ENDIF

    IF ( names .eq. 'CCCNAC' ) THEN
      lcn_ac = n
      CALL GET_VARIABLE(gd,'CCCNAC',c3)
      CALL SETCCNAC(lcn_ac,ccnac)
    ENDIF

    IF ( names .eq. 'CCCNNU' ) THEN
      lcn_nu = n
      CALL GET_VARIABLE(gd,'CCCNNU',c3)
      CALL SETCCNAC(lcn_nu,ccnnu)
    ENDIF

    IF ( names .eq. 'CCCNCO' ) THEN
      lcn_co = n
      CALL GET_VARIABLE(gd,'CCCNCO',c3)
      CALL SETCCNAC(lcn_co,ccnco)
    ENDIF

    IF ( names .eq. 'CCCNUF' ) THEN
      lccnuf = n
      CALL GET_VARIABLE(gd,'CCCNUF',c3)
      CALL SETCCNUF()
    ENDIF

    IF ( names .eq. 'CCIN' ) THEN
      lcin = n
      CALL GET_VARIABLE(gd,'CCIN',c3)
      CALL SETCIN()
    ENDIF

    IF ( names .eq. 'CPION' ) THEN
! get indices for CPION and CNION and mobility arrays
! then set up 1d arrays and fill 3d arrays with values
      lscpi = n
      lscni = GET_VARIABLE_INDEX(gd, 'CNION')
   
      lmupzc = GET_VARIABLE_INDEX(gd,'MUPOSZC')
      lmunzc = lmupzc + 1
      lmupze = lmunzc + 1
      lmunze = lmupze + 1

      CALL GET_VARIABLE(gd,'MUPOSZC',mupzc)
      CALL GET_VARIABLE(gd,'MUNEGZC',munzc)
      CALL GET_VARIABLE(gd,'MUPOSZE',mupze)
      CALL GET_VARIABLE(gd,'MUNEGZE',munze)
  
      CALL GET_VARIABLE(gd,'CPION',c3)
      CALL GET_VARIABLE(gd,'CNION',c4)
      CALL SETION()
    ENDIF

! turn on large ion processes if category is set (to unset, comment lines out in grid_create)
    IF ( names .eq. 'CPLION' ) THEN
     largeion = .true.
    END IF

  ENDDO
  
  IF ( prtbase .and. my_rank == 0 ) THEN
   allocate(    temp(-ng+1:nzend+ng),thetav(-ng+1:nzend+ng),pinit(-ng+1:nzend+ng),  &
                rinit(-ng+1:nzend+ng),rhw(-ng+1:nzend+ng),rhi(-ng+1:nzend+ng),dewpt(-ng+1:nzend+ng), &
                riinit(-ng+1:nzend+ng), rivinit(-ng+1:nzend+ng) )

   rivinit(:) = 0.0
   DO k = 1,nzend-1
     temp(k)   = tz(k)*pz(k)
     thetav(k) = tz(k)*(1.0 + 0.61*qz(k))
     pinit(k)  = 1000.*pz(k)**3.508
     rinit(k)  = pinit(k)*100./(287.04*tz(k)*pz(k))
     rhw(k)    = 3.8*exp(17.27*(tz(k)*pz(k)-273.16)/(tz(k)*pz(k)-36.))/pinit(k)
     rhw(k)    = qz(k) / rhw(k)
     rhi(k)    = 3.8*exp(21.87*(tz(k)*pz(k)-273.16)/(tz(k)*pz(k)-7.66))/pinit(k)
     rhi(k)    = qz(k) / rhi(k)
       
     ps        = 1.0e3*pz(k)**3.509
     td        = qz(k)*ps/(0.622+0.001*qz(k))                ! vapor pressure
     td        = max(td,0.001)                                     ! avoid problems near zero
     dewpt(k)  = 273.16 + (243.5/( (17.67/alog(td/6.112)) - 1.0))  ! Bolton's approximation
   ENDDO
   
   thetav(nzend) = tz(nzend)*(1.0 + 0.61*qz(nzend))

    DO k = 2,nzend-1
     bb          = 2.0*g*z1dinit(k,4)*(tz(k)-tz(k-1))/(tz(k)+tz(k-1))
     bbv         = 2.0*g*z1dinit(k,4)*(thetav(k)-thetav(k-1))/(thetav(k)+thetav(k-1))
     bsh         = (u1d(k)-u1d(k-1))**2 + (v1d(k)-v1d(k-1))**2
     bsh         = bsh * z1dinit(k,4)**2
     riinit (k-1) = min(bb/(bsh+1.0e-10),50.)
     rivinit(k-1) = min(bbv/(bsh+1.0e-10),50.)
    ENDDO

    riinit (nzend-1) = riinit(nzend-2)
    rivinit(nzend-1) = rivinit(nzend-2)

    write(luno,*)
    write(luno,"(1x,4x,'HEIGHT',5x,'PRESS',5x,'DEN',5x,'THETA',5x,'TEMP',4x,'THETAV',3x,'DEWPT',6x, &
   &  'QV',10x,'U',10x,'V',7x,'RHw',6x,'RHi',6x,'RI',6x,'RIV',6x,'PI0',6x,'1/PI0')")

    DO k = nzend-1,1,-1
     write(luno,"(1x,f11.2,3x,f7.2,3x,f7.4,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f7.4,3x,2(f8.3,3x),6(f6.2,3x))")  &
              z1dinit(k,1),pinit(k),rinit(k),tz(k),temp(k),thetav(k),dewpt(k),                           &
              qz(k)*1000.,u1d(k),v1d(k), rhw(k)*100., rhi(k)*100.,riinit(k),rivinit(k),pz(k),1./pz(k)
    ENDDO

    CALL GET_VARIABLE(gd,'DZ',dz)

    write(luno,*)
    write(luno,*)
    write(luno,*) 'z1d: level, zc, ze, dzc, dze'
    DO k = 1,nzend
      write(luno,'(i3,2x,5(1pe13.5,2x))') k, z1dinit(k,1),z1dinit(k,2),z1dinit(k,3),z1dinit(k,4), 1./z1dinit(k,3)
    ENDDO

     open(21,file=outsoundfile,form='formatted', status='unknown', position='append')
!   rewind(21)
!   write(21,*) psfc, tsfc, qsfc*1000.
     DO k = 1,nzend-1
      write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), tz(k), qz(k)*1000., u1d(k), v1d(k), 100.*pinit(k)
     ENDDO

   close(21)

   deallocate(temp,thetav,pinit,rinit,rhw,rhi,dewpt,riinit,rivinit)
  ENDIF ! prtbase
  
  
  
  
  
  deallocate ( den )
  deallocate( u1d,v1d,tz,qz,pz,kz )
!  deallocate( mupzc,munzc,mupze,munze )
!  deallocate( nion1d, pion1d )

!-----------------------------------------------------------------------------
! OPEN THE OUTPUT FILE

  luno = FILE_CLOSE('INIT_BACKGROUND_1D END')

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM DEFINITIONS (Fortran90 stuff....)

  CONTAINS

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SNDM1:  initialize with constant temperature and density environment with near-saturation (liquid or ice)
! Allows 1-D shear
!
  SUBROUTINE SNDM1()
  implicit none
    real cai,caw,cbi,cbw
    real :: chgt,shgt
    real z, an, ugnd, tv0, tv1

!   real psfc
!   parameter( psfc = 1.0e5)

      cai = 21.87455
      caw = 17.2693882
      cbi = 7.66
      cbw = 35.86

   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.


! SNDM1 does a moist, stationary, adiabatic base state with shear

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! PI INIT via vertical integration of hydrostatic equation

!#ifdef MPI
!   if (kzbeg .eq. nzbeg) tz(1) = tsfc 
!   if (kzbeg .eq. nzbeg) pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tsfc*cp)
!   
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 2
!   if (kzend .eq. nzend) kze = kzend-kzbeg+ng
!   DO k = kzb,kze
!    tz(k) = tsfc
!    pz(k) = pz(k-1) - g / (tsfc*cp*dze(k))
!   ENDDO
!#else
   tz(1) = tsfc 
   pz(1) = (psfc/1.e5)**rcp  ! - 0.5*g/(dzc(1)*tsfc*cp)

   DO k = 2,nzend-1
    tz(k) = tsfc
    pz(k) = pz(k-1) ! - g / (tsfc*cp*dze(k))
    pres  = psfc ! *pz(k)**(cp/rd)
    IF ( tsfc >= 273.15 ) THEN ! liquid water sat.
      qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    ELSEIF ( tsfc < 273.15 ) THEN ! set to ice sat
      qvs = 380*exp(cai*(pz(k)*tz(k)-273.15)/(pz(k)*tz(k)-cbi))/pres
    ENDIF
    qz(k) = min( qvs*0.998, qsfc )
   ENDDO
   qz(1) = qz(2)
   
!    qz(:) = 0.0
   
!#endif

   IF( wtype .eq. 1 ) THEN
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .lt. Ubase) u1d(k) = Usl*(zc(k)/Ubase)-Usl - Usmv
     IF( zc(k) .lt. (Uz+Ubase) .and. zc(k) .ge. Ubase ) u1d(k) = Us * ( (zc(k)-Ubase)/Uz ) - Usmv
     IF( zc(k) .ge. (Uz+Ubase) ) u1d(k) = Us - Usmv
    ENDDO
   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 2 ) THEN

	 ugnd = Us / pii
!#ifdef MPI
!     kzb = -ng+1
!     kze = ktile+ng
!     if (kzbeg .eq. nzbeg) kzb = 1
!     if (kzend .eq. nzend) kze = kzend-kzbeg
!     DO k = kzb,kze
!#else
     DO k = 1,nzend-1
!#endif
	  z = zc(k)
	 IF( z .le. Uz ) THEN
	  an = degtorad*(180.0 * (1.0 - z / Uz) )
	  u1d(k) = ugnd * (1.0 + cos(an))
	  v1d(k) = ugnd * sin(an)
	 ENDIF
	 IF( z .gt. Uz ) THEN
	  u1d(k) = 2*ugnd
	  v1d(k) = 0.0
	 ENDIF
	ENDDO
	
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 3 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF


! CREATE 2D Weisman Profile

   IF( wtype .eq. 4 ) THEN

      v1d(:) = 0.0

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
	 z = zc(k)
	 IF( z .le. Uz ) THEN
        u1d(k) = Us * ( zc(k)/Uz ) 
	 ENDIF
	 IF( z .gt. Uz .and. z .le. z0 ) THEN
	  u1d(k) = Us
	 ELSEIF ( z .gt. z0 .and. z .le. z1 ) THEN
        u1d(k) = Us + dudz0*(z - z0)/(z1 - z0 )
	 ELSEIF ( z .gt. z1 ) THEN
        u1d(k) = Us + dudz0
	 ENDIF
	ENDDO
	
   ENDIF


! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 5 ) THEN

      ugnd = (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 7000. ) THEN
        u1d(k) = ugnd + Us*(5./7.)*(z - 2000.)/5000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 7000.  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2.5km, then linear shear extended to 7km
! recreation of shear in WRF quarter circle supercell
   IF( wtype .eq. 6 ) THEN
   
    chgt = 2500.0
    shgt = 7000.0

      ugnd = (chgt*Us/shgt) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. chgt ) THEN
        an = degtorad*(90.0 * (2.0 - z / chgt) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. chgt .and. z .le. shgt ) THEN
        u1d(k) = ugnd + Us*((shgt-chgt)/shgt)*(z - chgt)/(shgt - chgt)
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. shgt  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

 END SUBROUTINE SNDM1


!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SND0:  Create dry, adiabatic no flow base state

  SUBROUTINE SND0()

   real psfc, tsfc
   parameter( psfc = 1.0e5, tsfc = 300. )


! SND0 does a dry, stationary, adiabatic base state

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! PI INIT via vertical integration of hydrostatic equation

!#ifdef MPI
!   if (kzbeg .eq. nzbeg) tz(1) = tsfc 
!   if (kzbeg .eq. nzbeg) pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tsfc*cp)
!   
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 2
!   if (kzend .eq. nzend) kze = kzend-kzbeg+ng
!   DO k = kzb,kze
!    tz(k) = tsfc
!    pz(k) = pz(k-1) - g / (tsfc*cp*dze(k))
!   ENDDO
!#else
   tz(1) = tsfc 
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tsfc*cp)

   DO k = 2,nzend-1
    tz(k) = tsfc
    pz(k) = pz(k-1) - g / (tsfc*cp*dze(k))
   ENDDO
   
    qz(:) = 0.0
   
!#endif

 END SUBROUTINE SND0

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SND1:  Creates WK82 sounding and wind profiles

  SUBROUTINE SND1()
  
  implicit none

   real :: chgt,shgt

! WK thermodynamic sounding parameters

!   real ztr    ;  parameter( ztr   = 12000. )
!   real thtr   ;  parameter( thtr  = 343.   )
!   real ttr    ;  parameter( ttr   = 213.   )

! Local arrays

   real rh(-2:1000)   !!-ng+1:1000
   real z, an, ugnd, ugnd2, tv0, tv1,v1,thet
   integer k0

   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.

! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! BASIC WK82 Sounding

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     rh(k) = 1. - 0.75 * zfac
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
     rh(k) = 0.25
    ENDIF
    rh(k) = Min( rh(k), rhmax )
    
    IF ( zrhmax2 > 0.0 ) THEN
      IF ( z > zrhmax2 ) rh(k) = Min( rh(k), rhmax2 )
    ENDIF

   ENDDO 
   
   qz(:) = 0.0
   
   DO i = 1,ntr

! PI INIT via vertical integration of hydrostatic equation
! Here we dont know Qv yet, and I am too lazy to iterate to find correct pressure

!#ifdef MPI
!   if (kzbeg .eq. nzbeg)   tv0   = tz(1) * (1.0 + 0.0*qz(1))
!   if (kzbeg .eq. nzbeg)   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#else
   tv0   = tz(1) * (1.0 + 0.61*qz(1))
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#endif

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 2
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 2,nzend-1
!#endif
    tv1   = tz(k) * (1.0 + 0.61*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1
   ENDDO

! CREATE Qv profile
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 1
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    pres  = psfc*pz(k)**(cp/rd)
    qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    qz(k) = min( rh(k) * qvs, qsfc )

   ENDDO
   
   ENDDO ! i

! CREATE 1D initial wind profile

   IF( wtype .eq. 0 ) THEN
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
      u1d(k) = Us * Tanh( zc(k)/Uz )
 
    ENDDO

   ENDIF
   
   IF( wtype .eq. 1 ) THEN
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .lt. Ubase) u1d(k) = Usl*(zc(k)/Ubase)-Usl - Usmv
     IF( zc(k) .lt. (Uz+Ubase) .and. zc(k) .ge. Ubase ) u1d(k) = Us * ( (zc(k)-Ubase)/Uz ) - Usmv
     IF( zc(k) .ge. (Uz+Ubase) ) u1d(k) = Us - Usmv
!     print*,'Here! k,zc(k),u1d(k)',k,zc(k),u1d(k)
 
    ENDDO

   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 2 ) THEN

	 ugnd = Us / pii
!#ifdef MPI
!     kzb = -ng+1
!     kze = ktile+ng
!     if (kzbeg .eq. nzbeg) kzb = 1
!     if (kzend .eq. nzend) kze = kzend-kzbeg
!     DO k = kzb,kze
!#else
     DO k = 1,nzend-1
!#endif
	  z = zc(k)
	 IF( z .le. Uz ) THEN
	  an = degtorad*(180.0 * (1.0 - z / Uz) )
	  u1d(k) = ugnd * (1.0 + cos(an))
	  v1d(k) = ugnd * sin(an)
	 ENDIF
	 IF( z .gt. Uz ) THEN
	  u1d(k) = 2*ugnd
	  v1d(k) = 0.0
	 ENDIF
	ENDDO
	
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 3 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF


! CREATE 2D Weisman Profile

   IF( wtype .eq. 4 ) THEN

      v1d(:) = 0.0

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
	 z = zc(k)
	 IF( z .le. Uz ) THEN
        u1d(k) = Us * ( zc(k)/Uz ) 
	 ENDIF
	 IF( z .gt. Uz .and. z .le. z0 ) THEN
	  u1d(k) = Us
	 ELSEIF ( z .gt. z0 .and. z .le. z1 ) THEN
        u1d(k) = Us + dudz0*(z - z0)/(z1 - z0 )
	 ELSEIF ( z .gt. z1 ) THEN
        u1d(k) = Us + dudz0
	 ENDIF
	ENDDO
	
   ENDIF


! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 5 ) THEN

      ugnd = (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 7000. ) THEN
        u1d(k) = ugnd + Us*(5./7.)*(z - 2000.)/5000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 7000.  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2.5km, then linear shear extended to 7km
! recreation of shear in WRF quarter circle supercell
   IF( wtype .eq. 6 ) THEN
   
    chgt = 2500.0
    shgt = 7000.0

      ugnd = (chgt*Us/shgt) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. chgt ) THEN
        an = degtorad*(90.0 * (2.0 - z / chgt) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. chgt .and. z .le. shgt ) THEN
        u1d(k) = ugnd + Us*((shgt-chgt)/shgt)*(z - chgt)/(shgt - chgt)
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. shgt  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC). Code to match CM1
! Hodograph length is about 35 m/s
   IF( wtype .eq. 7 ) THEN

      ugnd = 7.0 ! (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      ugnd2 = 31.0
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = 90.*degtorad*( z / 2000.) 
        u1d(k) = ugnd * (1.0 - cos(an))
        v1d(k) = ugnd * sin(an)
       ELSEIF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + (ugnd2-ugnd)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ELSE
        u1d(k) = ugnd2
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! rotate winds if needed

    IF ( rotdeg /= 0.0 ) THEN
      
      
      DO k = 1,nzend
        v1 = Sqrt( u1d(k)**2 + v1d(k)**2 )
        IF ( v1d(k) == 0.0 ) THEN
          IF ( u1d(k) >= 0.0 ) THEN
            thet = 0.0
          ELSE 
            thet = pii
          ENDIF
         ELSEIF ( u1d(k) == 0.0 ) THEN
          IF ( v1d(k) >= 0.0 ) THEN
            thet = pii/2
          ELSE 
            thet = -pii/2
          ENDIF
         ELSE
           thet = Atan( v1d(k)/u1d(k) )
         ENDIF
         
         thet = thet + degtorad * rotdeg
         v1d(k) = v1*Sin(thet)
         u1d(k) = v1*Cos(thet)
         
      ENDDO
    ENDIF

! Write sounding out to file

  if(my_rank == 0) then
   open(21,file=outsoundfile, status='unknown')
   rewind(21)
   write(21,'(1x,2(f9.2,2x),f9.4)') psfc, tsfc, qsfc*1000.
     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsfc,1000.*qsfc,0.0,0.0,psfc
!   DO k = 1,nz-1
!    write(21,'(1x, 5(f10.2,2x))') zc(k), tz(k), qz(k)*1000., u1d(k), v1d(k)
!   ENDDO

   close(21)

   end if

  END SUBROUTINE SND1

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SND2:  Reads in a sounding from an input file and interpolates
!                            the sounding to the grid.

   SUBROUTINE SND2()

! Define arrays for reading in the sounding

    implicit none

    real :: chgt,shgt

    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax), temk(nmax)
    real z0, tv0, tv1
    real z, an, ugnd, ugnd2

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac
    real newp000
    logical :: rescalep
      real :: pzsfc0,pzsfc, tksfc

! Error check


!------------------------------------------------------------------------------
! First read in the surface pressure, temperature, and mixing ratio
! Format is z,th,qv,u,v because of stretched grid
!     Sounding input file format: ascii
!     z (meters), theta (K), qv (g/kg),  u (m/s), v (m/s)
!------------------------------------------------------------------------------

    CALL FILE_MESSAGE('READING SOUNDING FROM: '//sndfile)

    open(unit=17,file=sndfile,status='old',form='formatted')
    rewind(17)

! Read in sfc pressure, theta, qv

    read(17,*) p000, t000, q000

    IF ( p000 < 0.0 ) THEN
    
    read(17,*) newp000
    p000 = Abs(p000)
    rescalep = .true.
    write(luno,*) 'rescaling for new sfc pressure ',newp000
    write(luno,*) 'WARNING: This does not work correctly'
    STOP
    IF( newp000 .lt. 1.0e4 ) newp000 = newp000*100.
    
    ELSE
    
    rescalep = .false.
    
    ENDIF
    
    
    write(luno,*) p000, t000, q000

   psfc = p000
   tsfc = t000
   qsfc = q000
   qfac = 1.0
   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qfac = 0.001
   
   qsfc = qfac*qsfc


! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! Read in sounding values

    nsnd = 0
    DO k = 1,nmax
	 read(17,*,iostat=ios) zsnd(k),tsnd(k),qvsnd(k),usnd(k),vsnd(k)
	 qvsnd(k) = qfac*qvsnd(k)
       IF( ios < 0 ) EXIT
	 nsnd = nsnd + 1
    ENDDO
    
    close(17)

!#ifdef MPI
!    if (kzbeg .eq. nzbeg)  zsfc = zsnd(1)
!#else
    zsfc = zsnd(1)
!#endif

!   IF( .not. SET_VARIABLE(gd,'HGT',       zsfc) ) write(6,*) 'INIT_GRID:  Problem setting HGT'

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .le. nsnd ) kze = nsnd-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nsnd
!#endif
      zsnd(k) = zsnd(k) - zsfc
    ENDDO

    zsnd(nsnd) = zsnd(nsnd) + 0.1

! Check to make sure the sounding levels span the computational grid

!#ifdef MPI
!    if (kzbeg .eq. nzbeg .AND. zsnd(1) .gt. zc(1) ) THEN
!#else
    IF( zsnd(1) .gt. zc(1) ) THEN
!#endif
	write(luno,*) 'ZMIN of SOUNDING   = ',zsnd(1)
	write(luno,*) 'ZMIN of MODEL GRID = ',zc(1)
	write(luno,*) 'ZMIN OF INPUT SOUNDING IS > ZMIN OF GRID'
	write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
	stop
    ENDIF

!#ifdef MPI
    IF( zsnd(nsnd) .lt. zc(nzend-1) ) THEN
!#else
!    IF( zsnd(nsnd) .lt. zc(nz-1) ) THEN
!#endif
	write(luno,*) 'ZMAX of SOUND      = ',zsnd(nsnd)
	write(luno,*) 'ZMAX of MODEL GRID = ',zc(nz-1)
	write(luno,*) 'ZMAX OF INPUT SOUNDING IS < ZMAX OF MODEL GRID'
	write(luno,*) 'Extrapolating sounding with isothermal layer....'
!      write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
!      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
!      stop
    ENDIF

! Now interpolate sounding to grid

    write(luno,*)
    write(luno,*) 'INTERPOLATING SOUNDING TO GRID '

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
      z0 = zc(k)
      IF ( z0 .lt. zsnd(nsnd) ) THEN
        u1d(k) = zfind(zsnd, usnd,  nsnd, z0)
        v1d(k) = zfind(zsnd, vsnd,  nsnd, z0)
        tz(k) = zfind(zsnd, tsnd,  nsnd, z0)
        qz(k) = zfind(zsnd, qvsnd, nsnd, z0)
        ktop = k
        IF ( k .eq. 1 ) THEN
          tv0   = tz(1) * (1.0 + 0.601*qz(1))
          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
          temk(k) = pz(k)*tz(k)
        ELSE
          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
          tmptr = pz(k)*tz(k)
          temk(k) = pz(k)*tz(k)
        ENDIF
      ELSE
        u1d(k) = u1d(ktop)
        v1d(k) = v1d(ktop)
        qz(k) = qz(ktop)
        tz(k) = tz(ktop)*(Exp((z0-zc(ktop) - zsnd(1))*g/cp/tmptr ) )

          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
          temk(k) = pz(k)*tz(k)

      ENDIF

!      IF( qz(k) .gt. 0.05 ) qz(k) = qz(k) / 1000.

    ENDDO
    
    IF ( rescalep ) THEN
      
      pzsfc0 = (psfc/1.e5)**rcp ! original surface Exner Pi
      tksfc = pzsfc0 * tsfc
      write(luno,*) 'old tsfc = ',tsfc
      ! get new surface Ex Pi and new theta using (real) surface temperature
      pzsfc = (newp000/1.e5)**rcp
      tsfc = tksfc/pzsfc
      write(luno,*) 'new tsfc = ',tsfc
      psfc = newp000
      t000 = tsfc
      
      tsnd(1) = tsfc
      
      DO k = 1,nzend-1
       write(luno,*) 'k, old pz = ',k,pz(k),tz(k)
       pz(k) = pz(k) - pzsfc0 + pzsfc
       tz(k) = temk(k)/pz(k)
       write(luno,*) 'k, temp pz = ',k,pz(k),tz(k)
      ENDDO
      
      DO k = 1,nzend-1
        
        IF ( k .eq. 1 ) THEN
           
           tv0   = tz(1) * (1.0 + 0.601*qz(1))
           pz(1) = (newp000/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
           tz(k) = temk(k)/pz(k)
!          tv0   = tz(1) * (1.0 + 0.601*qz(1))
!          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!          temk(k) = pz(k)*tz(k)
        ELSE
          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
!          tmptr = pz(k)*tz(k)
!          temk(k) = pz(k)*tz(k)
        ENDIF
        
        
        
      ENDDO

   IF( .not. SET_VARIABLE(gd,'PSFC',     newp000) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

    ENDIF

! PI INIT via vertical integration of hydrostatic equation
! Here we know Qv, so we can integrate using virtual temperature

!    tv0   = tz(1) * (1.0 + 0.601*qz(1))
!    pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)

    DO k = 2,nzend-1

!     tv1   = tz(k) * (1.0 + 0.601*qz(k))
!     pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!     tv0   = tv1

    ENDDO


! CREATE 1D initial wind profile

   IF( wtype .eq. 10 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .le. Uz ) u1d(k) = Us * ( zc(k)/Uz )
     IF( zc(k) .gt. Uz ) u1d(k) = Us 
     v1d(k) = 0.0

    ENDDO

   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 20 ) THEN

	ugnd = Us / pii

!#ifdef MPI
!    kzb = 1
!    kze = ktile
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
	 z = zc(k)
	 IF( z .le. Uz ) THEN
	  an = degtorad*(180.0 * (1.0 - z / Uz) )
	  u1d(k) = ugnd * (1.0 + cos(an))
	  v1d(k) = ugnd * sin(an)
	 ENDIF
	 IF( z .gt. Uz ) THEN
	  u1d(k) = 2*ugnd
	  v1d(k) = 0.0
	 ENDIF
	ENDDO
	
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 30 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF



! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 50 ) THEN

      ugnd = (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 7000. ) THEN
        u1d(k) = ugnd + Us*(5./7.)*(z - 2000.)/5000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 7000.  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2.5km, then linear shear extended to 7km
! recreation of shear in WRF quarter circle supercell
   IF( wtype .eq. 60 ) THEN
   
    chgt = 2500.0
    shgt = 7000.0

      ugnd = (chgt*Us/shgt) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. chgt ) THEN
        an = degtorad*(90.0 * (2.0 - z / chgt) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. chgt .and. z .le. shgt ) THEN
        u1d(k) = ugnd + Us*((shgt-chgt)/shgt)*(z - chgt)/(shgt - chgt)
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. shgt  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC). Code to match CM1
! Hodograph length is about 35 m/s
   IF( wtype .eq. 70 ) THEN

      ugnd = 7.0 ! (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      ugnd2 = 31.0
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = 90.*degtorad*( z / 2000.) 
        u1d(k) = ugnd * (1.0 - cos(an))
        v1d(k) = ugnd * sin(an)
       ELSEIF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + (ugnd2-ugnd)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ELSE
        u1d(k) = ugnd2
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! Flip winds:
   IF( wtype .eq. 1000 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        u1d(k) = -u1d(k)
        v1d(k) = -v1d(k)
      ENDDO

   ENDIF

! Swap winds:
   IF( wtype .eq. 1001 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        z = u1d(k)
        u1d(k) = v1d(k)
        v1d(k) = z
      ENDDO

   ENDIF

! Rotate winds:
   IF( wtype .eq. 1002 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        z = u1d(k)
        u1d(k) = -v1d(k)
        v1d(k) = z
      ENDDO

   ENDIF

! Write sounding out to file

    IF ( ny .le. 2 ) THEN
      v1d(:) = 0.0
    ELSEIF ( nx .le. 2 ) THEN
      u1d(:) = 0.0
    ENDIF
    

  if(my_rank == 0) then
! write out first two lines of sounding: surface values and first input sounding level (taken as surface)
    open(21,file=outsoundfile, status='unknown')
    rewind(21)
!    write(21,*) 'MODEL SND'
    write(21,'(1x,2(f9.2,2x),f9.4)') psfc, t000, qsfc*1000.
!    DO k = 1,nz-1
     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),psfc
!     write(21,'(1x, 6(f10.2,2x))') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),p000
!    ENDDO
    close(21)

  endif
   END SUBROUTINE SND2

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SND3:  Special USER soundings

   SUBROUTINE SND3()

! WK thermodynamic sounding parameters

!   real psfc   ;  parameter( psfc  = 1.0e5  ) 
!   real tsfc   ;  parameter( tsfc  = 300.   )
!   real qsfc   ;  parameter( qsfc  = .014   )
!    real ztr    ;  parameter( ztr   = 12000. )
!    real thtr   ;  parameter( thtr  = 343.   )
!    real ttr    ;  parameter( ttr   = 213.   )

! Curving wind shear parameters

    real shear0 ; parameter( shear0 = 0.024 )    ! Shear in layer between 0-z0
    real shear1 ; parameter( shear1 = 0.012 )    ! Shear between z0-z1
    real shear2 ; parameter( shear2 = 0.006 )    ! Shear between z1-z2
!   real z0     ; parameter( z0     = 0333. )    ! Height of turning layer
!   real z1     ; parameter( z1     = 1000. )    ! Height of 1st linear shear layer
!   real z2     ; parameter( z2     = 7000. )    ! Height of 2nd linear shear layer
!   real angle  ; parameter( angle  = 112.5 )    ! angle of low-level flow

    real degtorad ; parameter( degtorad = pii / 180.)
    real z, tv0, tv1, U0

! Namelist declarations

!    real psfc
!    real tsfc
!    real qsfc

!    NAMELIST /lowlevel_params/ shape, dudz0, dudz1, dudz2, z0, z1, z2, psfc, tsfc, qsfc

! Local arrays

    real rh(1000)

    write(luno,'(1x,71("-"))')
    write(luno,*) ' SND3 ROUTINE CALLED'

! Error check

   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.

! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! BASIC WK82 Sounding

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     rh(k) = 1. - 0.75 * zfac
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
     rh(k) = 0.25
    ENDIF 

    rh(k) = Min( rh(k), rhmax )
    
   ENDDO

! PI INIT via vertical integration of hydrostatic equation
! Here we dont know Qv yet, and I am too lazy to interate to find correct pressure

!#ifdef MPI
!   
!   if(kzbeg .eq. nzbeg) tv0   = tz(1) * (1.0 + 0.0*qz(1))
!   if(kzbeg .eq. nzbeg) pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 2
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!    tv1   = tz(k) * (1.0 + 0.0*qz(k))
!    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!    tv0   = tv1
!   ENDDO
!#else

   tv0   = tz(1) * (1.0 + 0.0*qz(1))
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
   
   DO k = 2,nzend-1
    tv1   = tz(k) * (1.0 + 0.0*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1
   ENDDO
!#endif

! CREATE Qv profile

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    pres  = psfc*pz(k)**(cp/rd)
    qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    qz(k) = min( rh(k) * qvs, qsfc )

   ENDDO
    
! Hack for experimental runs

    U0 = shear0 * z0 / sqrt(2.0)

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
	 z = zc(k)

	 IF( z .le. z0 ) then

	  IF( shape .eq. 0 ) THEN
	    u1d(k) = - U0 * z / z0
	  ELSE
	    u1d(k) =   U0 * (-2.0 + (z / z0))
	  ENDIF

	  v1d(k)   =   U0 * (z / z0)

	 ENDIF

	 IF( z .gt. z0 .and. z .le. z1 ) then

	  u1d(k) =   U0 * ( (z-z0) / (z1-z0) - 1.0 ) 
	  v1d(k) =   U0 * ( (z-z0) / (z1-z0) + 1.0 )

	 ENDIF

	 IF( z .gt. z1 .and. z .le. z2 ) then
	  u1d(k) = shear2 * (z-z1) 
	  v1d(k) = 2.0 * U0
	 ENDIF

	 IF( z .gt. z2 ) THEN

	  u1d(k) = shear2 * (z2-z1) 
	  v1d(k) = 2.0 * U0

	 ENDIF

	ENDDO


!        DO k = 1,nz
!
!          z = z1d(k,1)
!
!         IF (z .le. z0) THEN
!          an = degtorad*(270. - 180.*z/z0)
!          z1d(k,5)  = shear0*cos(an)
!          z1d(k,6)  = shear0*sin(an)
!         ENDIF
!
!         IF (z .gt. z0 .and. z .le. z1) THEN
!          z1d(k,5) = shear1 * (z-z0)
!          z1d(k,6) = shear0
!         ENDIF
!
!         IF (z .gt. z1 .and. z .le. z2) THEN
!          z1d(k,6) = shear0
!          z1d(k,5) = shear1 * (z1-z0) + shear2 * (z-z1)
!         ENDIF
!
!         IF( z .gt. z2 ) THEN
!          z1d(k,6) = shear0
!          z1d(k,5) = shear1 * (z1-z0) + shear2 * (z2-z1)
!         ENDIF
!
!        ENDDO                      
!
! Hack for experimental runs
!
!        DO k = 1,nz
!
!         z = z1d(k,1)
!
!         IF (z .le. z0) THEN
!          z1d(k,5)  = shear0*cos(angle*degtorad)*z/z0
!          z1d(k,6)  = shear0*sin(angle*degtorad)*z/z0
!         ENDIF
!
!         IF (z .gt. z0 .and. z .le. z1) THEN
!          z1d(k,5) = shear0*cos(angle*degtorad) + shear1 * (z-z0)
!          z1d(k,6) = shear0*sin(angle*degtorad)
!         ENDIF
!
!         IF (z .gt. z1 .and. z .le. z2) THEN
!          z1d(k,5) = shear0*cos(angle*degtorad) + shear1*(z1-z0) + shear2*(z-z1)
!          z1d(k,6) = shear0*sin(angle*degtorad)
!         ENDIF
!
!         IF( z .gt. z2 ) THEN
!          z1d(k,5) = shear0*cos(angle*degtorad) + shear1*(z1-z0) + shear2*(z2-z1)
!          z1d(k,6) = shear0*sin(angle*degtorad)
!         ENDIF
!
!        ENDDO                      

! Write sounding out to file

  if(my_rank == 0) then
    open(21,file=outsoundfile, status='unknown')
    rewind(21)
!    write(21,*) 'MODEL SND'
    write(21,'(1x,2(f9.2,2x),f9.4)') psfc, tsfc, qsfc*1000.
!    DO k = 1,nz-1
!      write(21,'(1x, 5(f10.2,2x))') zc(k), tz(k), qz(k)*1000., u1d(k), v1d(k)
!    ENDDO
    close(21)

  endif

   END SUBROUTINE SND3
   
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM SND4:  Creates WK82 sounding and wind profiles

  SUBROUTINE SND4()

! WK thermodynamic sounding parameters

!   real ztr    ;  parameter( ztr   = 12000. )
!   real thtr   ;  parameter( thtr  = 343.   )
!   real ttr    ;  parameter( ttr   = 213.   )
   real zz1, pp1, pold, tt1
   real deltbl
   integer ibl

! Local arrays

   real rh(1000)
   real z, an, ugnd, tv0, tv1

   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.

! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! BASIC WK82 Sounding

   pp1 = 81000.0
   deltbl = 0.2
   
!   pp1 = 83000.0
!   deltbl = 0.5

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     rh(k) = 1. - 0.75 * zfac
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
     rh(k) = 0.25
    ENDIF
    rh(k) = Min( rh(k), rhmax )

   ENDDO 

! PI INIT via vertical integration of hydrostatic equation
! Here we dont know Qv yet, and I am too lazy to interate to find correct pressure

   tv0   = tz(1) * (1.0 + 0.0*qz(1))
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)

   DO k = 2,nzend-1

    tv1   = tz(k) * (1.0 + 0.0*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1

   ENDDO

! CREATE Qv profile
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    pres  = psfc*pz(k)**(cp/rd)
     IF ( pres .lt. pp1 .and. pold .ge. pp1 ) THEN
       ibl = k - 1
       zz1 = zc(k-1) + (zc(k) - zc(k-1))*(pold - pp1)/(pold - pres)
       tt1 = tz(k-1) + (tz(k) - tz(k-1))*(pold - pp1)/(pold - pres)
!       print*, 'pp1,zz1,tt1 = ',pp1,zz1,tt1
     ENDIF
    qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    qz(k) = min( rh(k) * qvs, qsfc )
    pold = pres

   ENDDO


! create linear lapse rate in BL
   
   DO k = ibl,1,-1
     tz(k) = tt1 - deltbl*(zz1 - zc(k))/zz1
   ENDDO
   tsfc = tt1 - deltbl


! recalculate exner pressure
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 2
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
   DO k = 2,nzend-1
!#endif
    tv1   = tz(k) * (1.0 + 0.0*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1

   ENDDO
   

! CREATE 1D initial wind profile

   IF( wtype .eq. 0 ) THEN
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
      u1d(k) = Us * Tanh( zc(k)/Uz )
 
    ENDDO

   ENDIF


   IF( wtype .eq. 1 ) THEN
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .lt. Ubase) u1d(k) = Usl*(zc(k)/Ubase)-Usl - Usmv
     IF( zc(k) .lt. (Uz+Ubase) .and. zc(k) .ge. Ubase ) u1d(k) = Us * ( (zc(k)-Ubase)/Uz ) - Usmv
     IF( zc(k) .ge. (Uz+Ubase) ) u1d(k) = Us - Usmv
    ENDDO

   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 2 ) THEN

	ugnd = Us / pii

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
	 z = zc(k)
	 IF( z .le. Uz ) THEN
	  an = degtorad*(180.0 * (1.0 - z / Uz) )
	  u1d(k) = ugnd * (1.0 + cos(an))
	  v1d(k) = ugnd * sin(an)
	 ENDIF
	 IF( z .gt. Uz ) THEN
	  u1d(k) = 2*ugnd
	  v1d(k) = 0.0
	 ENDIF
	ENDDO
	
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 3 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   DO k = kzb,kze
!#else
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! Write sounding out to file

   if (my_rank ==0) then
   open(21,file=outsoundfile, status='unknown')
   rewind(21)
   write(21,'(1x,2(f9.2,2x),f9.4)') psfc, tsfc, qsfc*1000.
     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsfc,1000.*qsfc,0.0,0.0,psfc
!   DO k = 1,nz-1
!    write(21,'(1x, 5(f10.2,2x))') zc(k), tz(k), qz(k)*1000., u1d(k), v1d(k)
!   ENDDO

   close(21)

   endif

  END SUBROUTINE SND4
   
!--------------------------------------------------------------------------
! Subroutine  SETCCN:  initializes CCN concentrations
!
   SUBROUTINE SETCCN (  )
   
   implicit none
   real cwccn
   real, parameter :: rho00 = 1.225
   
   cwccn = Abs(ccn)
   IF ( microphys%str(1:3) == 'TAK' ) THEN
!     convert mks to cgs
     cwccn = 1.e-6*cwccn
   ENDIF
   
     IF ( ccn .lt. 0.0 ) THEN
!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
!#else     
!     DO k = 1,nz
!#endif
       gd%var(lccn)%base1d(k) = cwccn
       
       IF ( zc(k) > ccn_height1 .and. ccn_upperair > 0.0 ) THEN
         gd%var(lccn)%base1d(k) = ccn_upperair
       ENDIF
     ENDDO
     
     ELSE

!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     DO k = kzb,kze
!#else     
!     DO k = 1,nz-1
!#endif
       gd%var(lccn)%base1d(k) = Min(cwccn, cwccn*den(kzbeg+k-1)/rho00) ! den(1)
!       write(luno,*) 'CCN = ',k,gd%var(lccn)%base1d(k)

       IF ( zc(k) > ccn_height1 .and. ccn_upperair > 0.0 ) THEN
         gd%var(lccn)%base1d(k) = ccn_upperair*den(kzbeg+k-1)/rho00
       ENDIF

     ENDDO
     
     IF ( myprock == nprock ) THEN
       gd%var(lccn)%base1d(kzend-kzbeg+1) = gd%var(lccn)%base1d(nz-1)*den(nzend-1)/den(nzend-2)
     ENDIF
!       gd%var(lccn)%base1d(nz) = gd%var(lccn)%base1d(nz-1)*den(nz-1)/den(nz-2)
!       write(luno,*) 'CCN = ',k,gd%var(lccn)%base1d(nz)
!#endif
     ENDIF

! set base state init array:
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
      ccn1d(k) = gd%var(lccn)%base1d(k)
     ENDDO
#else  
     ccn1d(1:nz) = gd%var(lccn)%base1d(1:nz)
#endif

!#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     
     DO k = kzb,kze
      DO j = jyb,jye
       DO i = ixb, ixe
!#else
!      DO k = 1,nz-1
!       DO j = 1,ny-1
!        DO i = 1,nx-1
!#endif
          c3(i,j,k)   = gd%var(lccn)%base1d(k)      ! QV
        ENDDO
       ENDDO

!       gd%var(lccn)%base1d(k) = cz(k)

      ENDDO

   END SUBROUTINE SETCCN 

!--------------------------------------------------------------------------
! Subroutine  SETCCNAC:  initializes CCNAC concentrations
!
   SUBROUTINE SETCCNAC(lcn, ccnac)
   
   implicit none
   integer lcn
   real cwccn, ccnac
   real, parameter :: rho00 = 1.225
   
   cwccn = Abs(ccnac)
   IF ( microphys%str(1:3) == 'TAK' ) THEN
!     convert mks to cgs
     cwccn = 1.e-6*cwccn
   ENDIF
   
     IF ( ccnac .lt. 0.0 ) THEN
!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
!#else     
!     DO k = 1,nz
!#endif
       gd%var(lcn)%base1d(k) = cwccn
       
       IF ( zc(k) > ccn_height1 .and. ccn_upperair > 0.0 ) THEN
         gd%var(lcn)%base1d(k) = ccn_upperair
       ENDIF
     ENDDO
     
     ELSE

!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     DO k = kzb,kze
!#else     
!     DO k = 1,nz-1
!#endif
       gd%var(lcn)%base1d(k) = Min(cwccn, cwccn*den(kzbeg+k-1)/rho00) ! den(1)
!       write(luno,*) 'CCN = ',k,gd%var(lcn)%base1d(k)

       IF ( zc(k) > ccn_height1 .and. ccn_upperair > 0.0 ) THEN
         gd%var(lccn)%base1d(k) = ccn_upperair*den(kzbeg+k-1)/rho00
       ENDIF

     ENDDO
     
     IF ( myprock == nprock ) THEN
       gd%var(lcn)%base1d(kzend-kzbeg+1) = gd%var(lcn)%base1d(nz-1)*den(nzend-1)/den(nzend-2)
     ENDIF
!       gd%var(lcn)%base1d(nz) = gd%var(lcn)%base1d(nz-1)*den(nz-1)/den(nz-2)
!       write(luno,*) 'CCN = ',k,gd%var(lcn)%base1d(nz)
!#endif
     ENDIF

! set base state init array:
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
      ccn1d(k) = gd%var(lcn)%base1d(k)
     ENDDO
#else  
     ccn1d(1:nz) = gd%var(lcn)%base1d(1:nz)
#endif

!#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     
     DO k = kzb,kze
      DO j = jyb,jye
       DO i = ixb, ixe
!#else
!      DO k = 1,nz-1
!       DO j = 1,ny-1
!        DO i = 1,nx-1
!#endif
          c3(i,j,k)   = gd%var(lcn)%base1d(k)      ! QV
        ENDDO
       ENDDO

!       gd%var(lcn_ac)%base1d(k) = cz(k)

      ENDDO

   END SUBROUTINE SETCCNAC

!--------------------------------------------------------------------------
! Subroutine  SETCCNUF:  initializes UF CCN concentrations
!
   SUBROUTINE SETCCNUF (  )
   
   implicit none
   real cwccn
   real, parameter :: rho00 = 1.225
   
   ! write(0,*) 'setccnuf'
   
   cwccn = Max(0.0,ccnuf)
   IF ( microphys%str(1:3) == 'TAK' ) THEN
!     convert mks to cgs
     cwccn = 1.e-6*cwccn
   ENDIF
   
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     DO k = kzb,kze
       gd%var(lccnuf)%base1d(k) = Min(cwccn, cwccn*den(kzbeg+k-1)/rho00) ! den(1)
       write(luno,*) 'CCNUF = ',k,gd%var(lccnuf)%base1d(k)
     ENDDO
     
     IF ( myprock == nprock ) THEN
       gd%var(lccnuf)%base1d(kzend-kzbeg+1) = gd%var(lccnuf)%base1d(nz-1)*den(nzend-1)/den(nzend-2)
     ENDIF

! set base state init array:
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
      ccnuf1d(k) = gd%var(lccnuf)%base1d(k)
     ENDDO
#else  
     ccnuf1d(1:nz) = gd%var(lccnuf)%base1d(1:nz)
#endif

!#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     
     DO k = kzb,kze
      DO j = jyb,jye
       DO i = ixb, ixe
!#else
!      DO k = 1,nz-1
!       DO j = 1,ny-1
!        DO i = 1,nx-1
!#endif
          c3(i,j,k)   = gd%var(lccnuf)%base1d(k)      ! QV
        ENDDO
       ENDDO

!       gd%var(lccn)%base1d(k) = cz(k)

      ENDDO

   END SUBROUTINE SETCCNUF

   
!--------------------------------------------------------------------------
! Subroutine  SETCIN:  initializes CIN concentrations
!
   SUBROUTINE SETCIN (  )
   
   implicit none
   real ccin
   real, parameter :: rho00 = 1.225
   
   ccin = Abs(cin)
   IF ( microphys%str(1:3) == 'TAK' ) THEN
!     convert mks to cgs
     ccin = 1.e-6*ccin
   ENDIF
   
     IF ( ccn .lt. 0.0 ) THEN
!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
!#else     
!     DO k = 1,nz
!#endif
       gd%var(lcin)%base1d(k) = ccin
     ENDDO
     
     ELSE

!#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     DO k = kzb,kze
!#else     
!     DO k = 1,nz-1
!#endif
       gd%var(lcin)%base1d(k) = Min(ccin, ccin*den(kzbeg+k-1)/rho00) ! den(1)
     ENDDO
     
     IF ( myprock == nprock ) THEN
       gd%var(lcin)%base1d(kzend-kzbeg+1) = gd%var(lcin)%base1d(nz-1)*den(nzend-1)/den(nzend-2)
     ENDIF
!#endif
     ENDIF

! set base state init array:
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
      ccn1d(k) = gd%var(lcin)%base1d(k)
     ENDDO
#else  
     ccn1d(1:nz) = gd%var(lcin)%base1d(1:nz)
#endif

!#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg

     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     
     DO k = kzb,kze
      DO j = jyb,jye
       DO i = ixb, ixe
!#else
!      DO k = 1,nz-1
!       DO j = 1,ny-1
!        DO i = 1,nx-1
!#endif
          c3(i,j,k)   = gd%var(lcin)%base1d(k)      ! QV
        ENDDO
       ENDDO

!       gd%var(lccn)%base1d(k) = cz(k)

      ENDDO

   END SUBROUTINE SETCIN 


!--------------------------------------------------------------------------
! Subroutine  SETION:  initializes ion concentrations and mobilities
!
   SUBROUTINE SETION (  )
   
   implicit none
   
   real :: cion(nz,2)
   
   real, parameter :: ec = 1.602e-19 ! fundamental unit of charge
   real, parameter :: jc = -2.0e-12  ! fair weather current density (Amp m**-2) (Helsdon & Farley 1987; Chiu (1978) has -2.7e-12, quoting Gish (1944))
   real, parameter :: eperao = 8.8592e-12 ! permittivity of air

   double precision :: ez,dezdz
   
   double precision, parameter :: ezfairo = -80.00d0
   double precision, parameter :: efa1 = 4.5d-3, efa2 = 3.8d-4, efa3 = 1.0d-4
   double precision, parameter :: efb1 = 0.50, efb2 = 0.65, efb3 = 0.10
!   double precision, parameter :: ezfairs = -100.0d0        ! fair field at surface

   


! ion mobilities from Shreve 1970 / Mansell et al 2005
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
#else
     DO k = 1,nzend
#endif
      mupzc(k) = 1.4e-4*Exp(1.4e-4*zc(k))
      munzc(k) = 1.9e-4*Exp(1.4e-4*zc(k))
      mupze(k) = 1.4e-4*Exp(1.4e-4*ze(k))
      munze(k) = 1.9e-4*Exp(1.4e-4*ze(k))
     ENDDO
! fair weather field from Gish 1944 / Mansell et al 2005     
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
#else       
     DO k = 1,nz
#endif
      IF ( fairweather ) THEN
      ez = ezfairo *                             &
                ( efb1*exp(-efa1*zc(k))          &
                 +efb2*exp(-efa2*zc(k))          &
                 +efb3*exp(-efa3*zc(k)) )
      dezdz = ezfairo *                            &
                ( -efa1*efb1*exp(-efa1*zc(k))      &
                  -efa2*efb2*exp(-efa2*zc(k))      & 
                  -efa3*efb3*exp(-efa3*zc(k)) )

! positive ion concentration (e.g., Chiu 1978; Takahashi 1979; Helsdon and Farley 1987)
      cion(k,1) = ( Sqrt(crgfac)*jc/ez + munzc(k)*eperao*dezdz )/  &
                     ( ec*(mupzc(k)+munzc(k)) )
! negative ion concentration
      cion(k,2) = cion(k,1) - eperao*dezdz/ec 
      
      ELSE
       cion(k,1) = 0.0
       cion(k,2) = 0.0
      ENDIF
      
     ENDDO

! set base state ion density arrays
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
#else 
     DO k = 1,nz
#endif
      gd%var(lscpi)%base1d(k) = cion(k,1)
      gd%var(lscni)%base1d(k) = cion(k,2)
!      write(luno,*) 'CPION = ',k,gd%var(lscpi)%base1d(k)
!      write(luno,*) 'CNION = ',k,gd%var(lscni)%base1d(k)
!       write(6,*) 'CNION,CPION,NETION = ',k,gd%var(lscni)%base1d(k),gd%var(lscpi)%base1d(k),-gd%var(lscni)%base1d(k)+gd%var(lscpi)%base1d(k)
     ENDDO     

! set ion denisty init arrays CPIONINIT, CNIONINIT
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1
     DO k = kzb,kze
      pion1d(k) = gd%var(lscpi)%base1d(k)
      nion1d(k) = gd%var(lscni)%base1d(k)
     ENDDO
#else 
     pion1d(1:nz) = gd%var(lscpi)%base1d(1:nz)
     nion1d(1:nz) = gd%var(lscni)%base1d(1:nz)
#endif
! fill CPION and CNION 3d arrays    
#ifdef MPI
     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg
     
     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg
     
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg
     
     DO k = kzb,kze
      DO j = jyb,jye
       DO i = ixb, ixe
#else
     DO k = 1,nz-1
      DO j = 1,ny-1
       DO i = 1,nx-1
#endif
         c3(i,j,k)   = gd%var(lscpi)%base1d(k)      
         c4(i,j,k)   = gd%var(lscni)%base1d(k)      
       ENDDO
      ENDDO
     ENDDO   
   
   
   END SUBROUTINE SETION 

 END SUBROUTINE INIT_BACKGROUND_1D

!--------------------------------------------------------------------------
! FUNCTION ZHEIGHT:  Computes the height of a geometrically stretched grid
!                    with a few wrinkles:  It can have a layer of constant
!                    dz at the bottom 'n1' layers thick, it also limits
!                    the size of dz at the top of the model to be 'dzmax'.

 REAL FUNCTION ZHEIGHT(dzbot,r,nz,dzmax,n1,zctop,ztopr,dzmax2)

  implicit none
  integer nz, n1, k, k2
  integer n2
  real dzbot
  double precision r
  double precision sum
  double precision dznew, dzmaxdp, dzmax2dp
  real dzmax
  real zctop  ! height for upper level stretch
  real ztopr  ! upper level stretch factor
  real dzmax2 ! maximum upper dz
  real dzm

  sum = 0.0d0
  dzmaxdp = dzmax
  dzmax2dp = dzmax2
  
!  zctop = 10000.
!  ztopr = 1.09
!  dzmax2 = 1000. ! 2*dzmax
  
  n2 = 0

  DO k = 1,nz

   IF( k .le. n1 ) THEN
    dznew=dzbot
   ELSE
    k2=k-n1
    dznew = Min(dzbot * r**(k2-1),dzmaxdp)
      IF ( sum .ge. zctop ) THEN
       IF ( n2 .eq. 0 ) dzm = Min(dznew,dzmaxdp)
       n2 = n2 + 1
       dznew = Min(dzm * ztopr**n2, dzmax2)
      ENDIF
   ENDIF
   sum = sum + dznew

  ENDDO

  ZHEIGHT = sum

 RETURN
 END FUNCTION ZHEIGHT

!------------------------------------------------------------------------------
!
!   /////////////////////        BEGIN       \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   FUNCTION ZFIND   ////////////////////
!
! Zfind computes the value of a variable for a given hgt
!------------------------------------------------------------------------------
 REAL FUNCTION ZFIND(z,var,np,hgt)
 
  implicit none

  integer np, n, n0
  real z(np), var(np), hgt, linterp
  external linterp

!------------------------------------------------------------------------------
! Find the two hgt points which straddle hgt
! Check to see whether z is an increasing or decreasing function
!------------------------------------------------------------------------------

  IF( z(np) .gt. z(1) ) THEN
   DO n = 1,np
    IF( hgt .ge. z(n) .and. hgt .lt. z(n+1) ) then
	n0 = n
	EXIT
    ENDIF
   ENDDO
  ELSE
   DO n = 1,np
    IF( hgt .le. z(n) .and. hgt .gt. z(n+1) ) then
     n0 = n
     EXIT
    ENDIF
   ENDDO
  ENDIF
!------------------------------------------------------------------------------
! If no point found, then tell you and leave town
!------------------------------------------------------------------------------
!     if( n .eq. np ) then
!           print *, 'FUNCTION ZFIND -- NO POINT FOUND!!!!!!'
!           zfind = -999.
!           return
!     endif
!------------------------------------------------------------------------------
! Do a linear interpolation to get value
! Check for special cases
!------------------------------------------------------------------------------
  IF( n0 .eq. 1 .and. hgt .eq. z(1) ) THEN
    zfind = var(1)
    RETURN
  ENDIF

  IF( n0 .eq. np ) THEN
    zfind = var(np)
    RETURN
  ENDIF

  zfind = linterp(var(n0),var(n0+1),z(n0),z(n0+1),hgt)

 RETURN
 END FUNCTION ZFIND
!------------------------------------------------------------------------------
! Linear interpolation function
!------------------------------------------------------------------------------
 REAL FUNCTION LINTERP(q0, q1, z0, z1, z)

  implicit none
  
  real q0, q1, z0, z1, z, dz
  
  dz = z1 - z0
  linterp = ( (z1 - z) / dz ) * q0 + ( ( z - z0 ) / dz ) * q1

 RETURN
 END FUNCTION LINTERP


!--------------------------------------------------------------------------
! EXTERNAL SUBPROGRAM SND0:  Create dry, adiabatic no flow base state

  SUBROUTINE SND0X(gd,tz,pz,qz,dzc,dze,nzend,ng,psfc,tsfc,qsfc)

   USE GRID_MODULE
   USE PARAM_MODULE, only: g,rcp,cp
   
   implicit none
   
   integer :: nzend,ng
   real    :: tz(-ng+1:nzend+ng),pz(-ng+1:nzend+ng),qz(-ng+1:nzend+ng),dzc(-ng+1:nzend+ng),dze(-ng+1:nzend+ng)
   
   type(grid) :: gd
   
   real psfc, tsfc, qsfc
!   parameter( psfc = 1.0e5, tsfc = 300. )
   integer :: k

    qsfc = 0.0
! SND0 does a dry, stationary, adiabatic base state

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! PI INIT via vertical integration of hydrostatic equation

!#ifdef MPI
!   if (kzbeg .eq. nzbeg) tz(1) = tsfc 
!   if (kzbeg .eq. nzbeg) pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tsfc*cp)
!   
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 2
!   if (kzend .eq. nzend) kze = kzend-kzbeg+ng
!   DO k = kzb,kze
!    tz(k) = tsfc
!    pz(k) = pz(k-1) - g / (tsfc*cp*dze(k))
!   ENDDO
!#else
   tz(1) = tsfc 
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tsfc*cp)

   DO k = 2,nzend-1
    tz(k) = tsfc
    pz(k) = pz(k-1) - g / (tsfc*cp*dze(k))
   ENDDO
   
    qz(:) = 0.0
   
!#endif

 END SUBROUTINE SND0X
 
!--------------------------------------------------------------------------
! EXTERNAL SUBPROGRAM SND1:  Creates WK82 sounding and wind profiles
! CALL       SND1X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
!                  zrhmax2,zrhdel,wtype,Us,Uz,Usl,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
!                  outsoundfile) 
  SUBROUTINE SNDWKTHERMO(tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   iqflag)
  
!   USE GRID_MODULE
   USE PARAM_MODULE, only: g,rcp,cp,rd,degtorad,pii
   USE commasmpi_module, only: my_rank
   USE init_module, only : wkexponent,ztr,thtr,ttr,zbl,thbl,dthdzbl,chgtwrf,dzcaplayer,dthcaplayer
   
   implicit none

   integer :: nzend,ng
   real    :: tz(-ng+1:nzend+ng),pz(-ng+1:nzend+ng),qz(-ng+1:nzend+ng),dzc(-ng+1:nzend+ng),dze(-ng+1:nzend+ng),zc(-ng+1:nzend+ng)
   real    :: u1d(-ng+1:nzend+ng),v1d(-ng+1:nzend+ng)
   
!   type(grid) :: gd
   real psfc, tsfc, qsfc, rhmax, rhmax2, zrhmax2,zrhdel
   integer :: wtype,ntr
   real    :: Us,Uz,Usl,Usmv, Ubase, z0, z1, dudz0, rotdeg
   integer, intent(in) :: iqflag ! whether to set qv or not
      

   real :: chgt,shgt,pres
   real :: qvs
   integer :: i,k

! WK thermodynamic sounding parameters

!   real ztr    ;  parameter( ztr   = 12000. )
!   real thtr   ;  parameter( thtr  = 343.   )
!   real ttr    ;  parameter( ttr   = 213.   )

! Local arrays

   real rh(-2:1000)   !!-ng+1:1000
   real z, an, ugnd, ugnd2, tv0, tv1,v1,thet
   real zfac, zfac2, tmp, zdiff
   integer k0

!    IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
!    IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.
! 
! ! SET some defaults
! 
!    IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
!    IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
!    IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! BASIC WK82 Sounding
!                         zbl, p_bl, thbl
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     IF ( zbl > 0. .and. ( thbl > tsfc .or. dthdzbl > 0.0 ) ) THEN
       IF ( dthdzbl > 0.0 ) THEN
        ! let dthdzbl supercede thbl
        thbl = tsfc + zbl*dthdzbl
       ENDIF
       IF ( z <= zbl ) THEN
         zfac  = ( z / zbl ) ! linear increase from surface to zbl
         tz(k) = tsfc + ( thbl - tsfc ) * zfac
         rh(k) = rhmax ! 1. - 0.75 * zfac
       ELSEIF ( dzcaplayer > 0. .and. z <= zbl + dzcaplayer ) THEN
         zfac  = ( (z - zbl) / dzcaplayer ) ! linear increase from surface to zbl
         tz(k) =  thbl + ( dthcaplayer ) * zfac
         zfac2  = ( (z - zbl) / (ztr - zbl) ) ** wkexponent
         rh(k) = Min( rhmax, 1. - 0.75 * zfac2)
       ELSE
         zfac  = ( (z - zbl - dzcaplayer) / (ztr - zbl - dzcaplayer) ) ** wkexponent
         tz(k) = thbl + dthcaplayer + ( thtr - thbl - dthcaplayer ) * zfac
         zfac2  = ( (z - zbl) / (ztr - zbl) ) ** wkexponent
         rh(k) = 1. - 0.75 * zfac2
       ENDIF
!       write(0,*) 'k,zc rh1 = ',k,zc(k),zfac,rh(k)
     ELSE
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     rh(k) = 1. - 0.75 * zfac
!       write(0,*) 'k,zc rh2 = ',k,zc(k),zfac,rh(k)
     ENDIF
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
     rh(k) = 0.25
    ENDIF
    rh(k) = Min( rh(k), rhmax )
    
    IF ( zrhmax2 > 0.0 .and. z > zrhmax2 ) THEN
    
      IF ( z > zrhmax2 .and. z <  zrhdel + zrhmax2 ) THEN
        tmp =  rhmax - (rhmax - rhmax2)*(z - zrhmax2)/zrhdel
      ELSE
        tmp = rhmax2
      ENDIF
        rh(k) = Min( rh(k), tmp )
      
    ENDIF

!    write(0,*) 'k,zc = ',k,zc(k),zfac,tz(k),rh(k),zbl,thbl

   ENDDO 
   
   IF ( iqflag > 0 ) THEN
   qz(:) = 0.0
!   v1d(:) = 0.0
!   u1d(:) = 0.0
   
   DO i = 1,ntr

! PI INIT via vertical integration of hydrostatic equation
! Here we dont know Qv yet, and I am too lazy to interate to find correct pressure

!#ifdef MPI
!   if (kzbeg .eq. nzbeg)   tv0   = tz(1) * (1.0 + 0.0*qz(1))
!   if (kzbeg .eq. nzbeg)   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#else
   tv0   = tz(1) * (1.0 + 0.61*qz(1))
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#endif

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 2
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 2,nzend-1
!#endif
    tv1   = tz(k) * (1.0 + 0.61*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1
   ENDDO

! CREATE Qv profile
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 1
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    pres  = 100000.*pz(k)**(cp/rd) ! needs to be 1e5 here, not psfc, since pz is already scaled to psfc
    qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    qz(k) = min( rh(k) * qvs, qsfc )

   ENDDO
   
   ENDDO ! i
   
   ENDIF ! iqflag

   END SUBROUTINE SNDWKTHERMO
!--------------------------------------------------------------------------
! EXTERNAL SUBPROGRAM SND1:  Creates WK82 sounding and wind profiles
! CALL       SND1X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
!                  zrhmax2,zrhdel,wtype,Us,Uz,Usl,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
!                  outsoundfile) 
  SUBROUTINE SND1X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   outsoundfile,wkumax1,wkumax2,wkvmax,wkh1,wkh2)
  
   USE GRID_MODULE
   USE PARAM_MODULE, only: g,rcp,cp,rd,degtorad,pii
   USE commasmpi_module, only: my_rank
   USE init_module, only : wkexponent,ztr,thtr,ttr,zbl,thbl,dthdzbl,chgtwrf,dzcaplayer,dthcaplayer, &
                           imodw1, imodw2, imodw3, ich,                                     &
                           ushrl1,ushrlb1,ushrlt1,vshrl1,vshrlb1,vshrlt1,                   &
                           ushrl2,ushrlb2,ushrlt2,vshrl2,vshrlb2,vshrlt2,                   &
                           ushrl3,ushrlb3,ushrlt3,vshrl3,vshrlb3,vshrlt3


   
   implicit none

   integer :: nzend,ng
   real    :: tz(-ng+1:nzend+ng),pz(-ng+1:nzend+ng),qz(-ng+1:nzend+ng),dzc(-ng+1:nzend+ng),dze(-ng+1:nzend+ng),zc(-ng+1:nzend+ng)
   real    :: u1d(-ng+1:nzend+ng),v1d(-ng+1:nzend+ng)
   
   type(grid) :: gd
   real psfc, tsfc, qsfc, rhmax, rhmax2, zrhmax2,zrhdel
   integer :: wtype,ntr
   real    :: Us,Uz,Usl,Usmv, Ubase, z0, z1, dudz0, rotdeg
   real    :: wkumax1,wkumax2,wkvmax,wkh1,wkh2,zshr
   
   character(*) :: outsoundfile
   

   real :: chgt,shgt,pres
   real :: qvs
   integer :: i,k

! WK thermodynamic sounding parameters

!   real ztr    ;  parameter( ztr   = 12000. )
!   real thtr   ;  parameter( thtr  = 343.   )
!   real ttr    ;  parameter( ttr   = 213.   )

! Local arrays

   real rh(-2:1000)   !!-ng+1:1000
   real z, an, ugnd, ugnd2, tv0, tv1,v1,thet
   real zfac, zfac2, tmp, zdiff
   real u0,u1,u2,u3,v0,v2,v3
   integer k0

   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qsfc = qsfc/1000.

! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! BASIC WK82 Sounding
!                         zbl, p_bl, thbl
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     IF ( zbl > 0. .and. ( thbl > tsfc .or. dthdzbl > 0.0 ) ) THEN
       IF ( dthdzbl > 0.0 ) THEN
        ! let dthdzbl supercede thbl
        thbl = tsfc + zbl*dthdzbl
       ENDIF
       IF ( z <= zbl ) THEN
         zfac  = ( z / zbl ) ! linear increase from surface to zbl
         tz(k) = tsfc + ( thbl - tsfc ) * zfac
         rh(k) = rhmax ! 1. - 0.75 * zfac
       ELSEIF ( dzcaplayer > 0. .and. z <= zbl + dzcaplayer ) THEN
         zfac  = ( (z - zbl) / dzcaplayer ) ! linear increase from surface to zbl
         tz(k) =  thbl + ( dthcaplayer ) * zfac
         zfac2  = ( (z - zbl) / (ztr - zbl) ) ** wkexponent
         rh(k) = Min( rhmax, 1. - 0.75 * zfac2)
       ELSE
         zfac  = ( (z - zbl - dzcaplayer) / (ztr - zbl - dzcaplayer) ) ** wkexponent
         tz(k) = thbl + dthcaplayer + ( thtr - thbl - dthcaplayer ) * zfac
         zfac2  = ( (z - zbl) / (ztr - zbl) ) ** wkexponent
         rh(k) = 1. - 0.75 * zfac2
       ENDIF
!       write(0,*) 'k,zc rh1 = ',k,zc(k),zfac,rh(k)
     ELSE
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     rh(k) = 1. - 0.75 * zfac
!       write(0,*) 'k,zc rh2 = ',k,zc(k),zfac,rh(k)
     ENDIF
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
     rh(k) = 0.25
    ENDIF
    rh(k) = Min( rh(k), rhmax )
    
    IF ( zrhmax2 > 0.0 .and. z > zrhmax2 ) THEN
    
      IF ( z > zrhmax2 .and. z <  zrhdel + zrhmax2 ) THEN
        tmp =  rhmax - (rhmax - rhmax2)*(z - zrhmax2)/zrhdel
      ELSE
        tmp = rhmax2
      ENDIF
        rh(k) = Min( rh(k), tmp )
      
    ENDIF

!    write(0,*) 'k,zc = ',k,zc(k),zfac,tz(k),rh(k),zbl,thbl

   ENDDO 
   
   qz(:) = 0.0
   v1d(:) = 0.0
   u1d(:) = 0.0
   
   DO i = 1,ntr

! PI INIT via vertical integration of hydrostatic equation
! Here we dont know Qv yet, and I am too lazy to interate to find correct pressure

!#ifdef MPI
!   if (kzbeg .eq. nzbeg)   tv0   = tz(1) * (1.0 + 0.0*qz(1))
!   if (kzbeg .eq. nzbeg)   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#else
   tv0   = tz(1) * (1.0 + 0.61*qz(1))
   pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!#endif

!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 2
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 2,nzend-1
!#endif
    tv1   = tz(k) * (1.0 + 0.61*qz(k))
    pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
    tv0   = tv1
   ENDDO

! CREATE Qv profile
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if(kzbeg .eq. nzbeg) kzb = 1
!   if(kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    pres  = 100000.*pz(k)**(cp/rd) ! needs to be 1e5 here, not psfc, since pz is already scaled to psfc
    qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
    qz(k) = min( rh(k) * qvs, qsfc )

   ENDDO
   
   ENDDO ! i

! CREATE 1D initial wind profile

   IF( wtype .eq. 0 ) THEN
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
      u1d(k) = Us * Tanh( zc(k)/Uz )
 
    ENDDO

   ENDIF
   
   IF( wtype .eq. 1 .or. wtype .eq. -1 ) THEN
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .lt. Ubase) u1d(k) = Usl*(zc(k)/Ubase)-Usl - Usmv
     IF( zc(k) .lt. (Uz+Ubase) .and. zc(k) .ge. Ubase ) u1d(k) = Us * ( (zc(k)-Ubase)/Uz ) - Usmv
     IF( zc(k) .ge. (Uz+Ubase) ) u1d(k) = Us - Usmv
     IF ( wtype == -1 ) THEN
       v1d(k) = u1d(k)
       u1d(k) = 0.0
     ENDIF
!     print*,'Here! k,zc(k),u1d(k)',k,zc(k),u1d(k)
 
    ENDDO

   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 2 ) THEN

     ugnd = Us / pii
!#ifdef MPI
!     kzb = -ng+1
!     kze = ktile+ng
!     if (kzbeg .eq. nzbeg) kzb = 1
!     if (kzend .eq. nzend) kze = kzend-kzbeg
!     DO k = kzb,kze
!#else
     DO k = 1,nzend-1
!#endif
      z = zc(k)
     IF( z .le. Uz ) THEN
      an = degtorad*(180.0 * (1.0 - z / Uz) )
      u1d(k) = ugnd * (1.0 + cos(an))
      v1d(k) = ugnd * sin(an)
     ENDIF
     IF( z .gt. Uz ) THEN
      u1d(k) = 2*ugnd
      v1d(k) = 0.0
     ENDIF
    ENDDO
    
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick/candy cane (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 3 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF


! CREATE 2D Weisman Profile

   IF( wtype .eq. 4 ) THEN

      v1d(:) = 0.0

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
	 z = zc(k)
	 IF( z .le. Uz ) THEN
        u1d(k) = Us * ( zc(k)/Uz ) 
	 ENDIF
	 IF( z .gt. Uz .and. z .le. z0 ) THEN
	  u1d(k) = Us
	 ELSEIF ( z .gt. z0 .and. z .le. z1 ) THEN
        u1d(k) = Us + dudz0*(z - z0)/(z1 - z0 )
	 ELSEIF ( z .gt. z1 ) THEN
        u1d(k) = Us + dudz0
	 ENDIF
	ENDDO
	
   ENDIF


! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 5 ) THEN

      ugnd = (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 7000. ) THEN
        u1d(k) = ugnd + Us*(5./7.)*(z - 2000.)/5000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 7000.  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2.5km, then linear shear extended to 7km
! recreation of shear in WRF quarter circle supercell
   IF( wtype .eq. 6 ) THEN
   
    chgt = chgtwrf ! 2300.0
    shgt = 7000.0

!      ugnd = (chgt*Us/shgt) / (pii/2.)  ! radius of quarter circle
      ugnd = (Us/3.0)/(pii/2.0)
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. chgt ) THEN
        an = degtorad*(90.0 * (2.0 - z / chgt) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. chgt .and. z .le. shgt ) THEN
!        u1d(k) = ugnd + Us*((shgt-chgt)/shgt)*(z - chgt)/(shgt - chgt)
        u1d(k) = ugnd + Us*(2./3.)*(z - chgt)/(shgt - chgt)
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. shgt  ) THEN
!        u1d(k) = ugnd + Us*(5./7.)
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC). Code to match CM1
! Hodograph length is about 35 m/s
   IF( wtype .eq. 7 ) THEN

      ugnd = 7.0 ! (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      ugnd2 = 31.0
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = 90.*degtorad*( z / 2000.) 
        u1d(k) = ugnd * (1.0 - cos(an))
        v1d(k) = ugnd * sin(an)
       ELSEIF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + (ugnd2-ugnd)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ELSE
        u1d(k) = ugnd2
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC). 
! Dennis and Kumjian (2017) adjustable parameters

   IF( wtype .eq. 8 ) THEN

      ugnd = wkumax1 ! (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      ugnd2 = wkumax2 ! 31.0
            
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. wkh1 ) THEN
        an = 90.*degtorad*( z / wkh1) 
        u1d(k) = ugnd * (1.0 - cos(an))
        v1d(k) = wkvmax * sin(an)
       ELSEIF( z .gt. wkh1 .and. z .le. wkh2 ) THEN
        u1d(k) = ugnd + (ugnd2-ugnd)*(z - wkh1)/(wkh2 - wkh1)
        v1d(k) = wkvmax
       ELSE
        u1d(k) = ugnd2
        v1d(k) = wkvmax
       ENDIF
      ENDDO
      
   ENDIF

! 
!  for ich .eq. 1, set up LP storm hodograph.  
!   uset is now an overall rotation angle (in degrees) 
!     (applied later below)
!     
! 
      IF ( ich .eq. 1 ) THEN

       ugnd = Us / (pii)
       u0 = ugnd
       v0 = ugnd
       zshr = Uz
       k0 = 1
       DO k = 1,nzend-1
        z = zc(k)
        IF( z .le. Uz ) THEN
          an = degtorad*(90.0 * (2.0 - z / Uz) )
          u1d(k) = ugnd * (1.0 + cos(an))
          v1d(k) = ugnd * sin(an)
          k0 = k
        ENDIF
        IF( z .gt. Uz ) THEN
          u1d(k) = ugnd ! u1d(k0)
          v1d(k) = ugnd ! v1d(k0)
        ENDIF
       ENDDO

       if ( imodw1 >= 1 ) then
        ! Find U and V at top of respective shear layers
         u1 = u0 + ushrl1*(ushrlt1 - ushrlb1)
         v1 = v0 + vshrl1*(vshrlt1 - vshrlb1)
        DO k = 1,nzend-1
         z = zc(k)
         if ( z .ge. vshrlb1 .and. z .le. vshrlt1 ) then
!c       v(kz) = v(kzcurv) + (z(kz)-vshrlb1)*vshrl1
           v1d(k) = v0 + vshrl1*(z - vshrlb1) !  v(k-1) + dz*vshrl1
!        kztopv = kz
         elseif ( z > vshrlt1 ) then
           v1d(k) = v1
         end if
         if ( z .ge. ushrlb1 .and. z .le. ushrlt1 ) then
!c       u(kz) = u(kzcurv) + (z(kz)-ushrlb1)*ushrl1
!          u(kz) = u(kz-1) + dz*ushrl1
           u1d(k) = u0 + ushrl1*(z - ushrlb1) !  v(k-1) + dz*vshrl1
!          kztopu = kz
         elseif ( z > ushrlt1 ) then
           u1d(k) = u1
         end if
        ENDDO
        u0 = u1
        v0 = v1
       end if

       if ( imodw2 >= 1 ) then
        ! Find U and V at top of respective shear layers
         u1 = u0 + ushrl2*(ushrlt2 - ushrlb2)
         v1 = v0 + vshrl2*(vshrlt2 - vshrlb2)
        DO k = 1,nzend-1
         z = zc(k)
         if ( z .ge. vshrlb2 .and. z .le. vshrlt2 ) then
!c       v(kz) = v(kzcurv) + (z(kz)-vshrlb1)*vshrl1
           v1d(k) = v0 + vshrl2*(z - vshrlb2) !  v(k-1) + dz*vshrl1
!        kztopv = kz
         elseif ( z > vshrlt2 ) then
           v1d(k) = v1
         end if
         if ( z .ge. ushrlb2 .and. z .le. ushrlt2 ) then
!c       u(kz) = u(kzcurv) + (z(kz)-ushrlb1)*ushrl1
!          u(kz) = u(kz-1) + dz*ushrl1
           u1d(k) = u0 + ushrl2*(z - ushrlb2) !  v(k-1) + dz*vshrl1
!          kztopu = kz
         elseif ( z > ushrlt2 ) then
           u1d(k) = u1
         end if
        ENDDO
        u0 = u1
        v0 = v1
       end if

       if ( imodw3 >= 1 ) then
        ! Find U and V at top of respective shear layers
         u1 = u0 + ushrl3*(ushrlt3 - ushrlb3)
         v1 = v0 + vshrl3*(vshrlt3 - vshrlb3)
        DO k = 1,nzend-1
         z = zc(k)
         if ( z .ge. vshrlb3 .and. z .le. vshrlt3 ) then
!c       v(kz) = v(kzcurv) + (z(kz)-vshrlb1)*vshrl1
           v1d(k) = v0 + vshrl3*(z - vshrlb3) !  v(k-1) + dz*vshrl1
!        kztopv = kz
         elseif ( z > vshrlt3 ) then
           v1d(k) = v1
         end if
         if ( z .ge. ushrlb3 .and. z .le. ushrlt3 ) then
!c       u(kz) = u(kzcurv) + (z(kz)-ushrlb1)*ushrl1
!          u(kz) = u(kz-1) + dz*ushrl1
           u1d(k) = u0 + ushrl3*(z - ushrlb3) !  v(k-1) + dz*vshrl1
!          kztopu = kz
         elseif ( z > ushrlt3 ) then
           u1d(k) = u1
         end if
        ENDDO
        u0 = u1
        v0 = v1
       end if

!        DO k = 1, nzend
!        z = zc(k)
!        kzcurv = (zshr-zsfc)/dz + 1
!        if ( z(kz) .le. zshr ) then
!         dd = (kz-1)*degturn/(kzcurv-1)
!      >       + deginit
!         vshr = ushr/2.
!           u(kz) = vshr*cos(4.*atan(1.)/180.*(270.-dd))
!           v(kz) = vshr*sin(4.*atan(1.)/180.*(270.-dd))
!           u(kz) = u(kz) + vshr
!         else
!         u(kz) = u(kzcurv)
!         v(kz) = v(kzcurv)
!         end if
!       end if
      ENDIF ! ich
! rotate winds if needed

    IF ( rotdeg /= 0.0 ) THEN
      
      
      DO k = 1,nzend
        v1 = Sqrt( u1d(k)**2 + v1d(k)**2 )
        IF ( v1d(k) == 0.0 ) THEN
          IF ( u1d(k) >= 0.0 ) THEN
            thet = 0.0
          ELSE 
            thet = pii
          ENDIF
         ELSEIF ( u1d(k) == 0.0 ) THEN
          IF ( v1d(k) >= 0.0 ) THEN
            thet = pii/2
          ELSE 
            thet = -pii/2
          ENDIF
         ELSE
           thet = Atan( v1d(k)/u1d(k) )
         ENDIF
         
         thet = thet + degtorad * rotdeg
         v1d(k) = v1*Sin(thet)
         u1d(k) = v1*Cos(thet)
         
      ENDDO
    ENDIF

! Write sounding out to file

  if(my_rank == 0) then
   open(21,file=outsoundfile, status='unknown')
   rewind(21)
   write(21,'(1x,2(f9.2,2x),f9.4)') psfc, tsfc, qsfc*1000.
     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsfc,1000.*qsfc,0.0,0.0,psfc
!   DO k = 1,nz-1
!    write(21,'(1x, 5(f10.2,2x))') zc(k), tz(k), qz(k)*1000., u1d(k), v1d(k)
!   ENDDO

   close(21)

   end if

  END SUBROUTINE SND1X
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! EXTERNAL SUBPROGRAM SND2X:  Reads in a sounding from an input file and interpolates
!                            the sounding to the grid.
!--------------------------------------------------------------------------

  SUBROUTINE SND2X(gd,tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmaxin,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   outsoundfile,sndfile)
  
   USE GRID_MODULE
   USE PARAM_MODULE, only: g,rcp,cp,rd,degtorad,pii,luno
   USE commasmpi_module, only: my_rank
   USE FILE_MODULE, only: file_message,file_close
   USE init_module, only : wkexponent,ztr,thtr,ttr,zbl,thbl,dthdzbl
   
   implicit none

   integer :: nzend,ng
   real    :: tz(-ng+1:nzend+ng),pz(-ng+1:nzend+ng),qz(-ng+1:nzend+ng),dzc(-ng+1:nzend+ng),dze(-ng+1:nzend+ng),zc(-ng+1:nzend+ng)
   real    :: u1d(-ng+1:nzend+ng),v1d(-ng+1:nzend+ng)
   
   type(grid) :: gd
   real psfc, tsfc, qsfc, rhmaxin, rhmax, rhmax2, zrhmax2, zrhdel
   integer :: qvorrh = 0
   integer :: wtype,ntr
   real    :: Us,Uz,Usl,Usmv, Ubase, z0, z1, dudz0, rotdeg
   
   character(*) :: outsoundfile
   character(*) :: sndfile
   
   real :: chgt,shgt,pres,temp,qsfcwk
   real :: qvs
   integer :: i,k
   integer :: nx,ny
   real zfac

!--------------------------------------------------------------------------
! Define arrays for reading in the sounding

    integer nmax, nsnd, ios
    parameter ( nmax = 5000 )
    real rh(-2:nmax)   !!-ng+1:1000
    real zsnd(nmax), tsnd(nmax), qvsnd(nmax)
    real usnd(nmax), vsnd(nmax), temk(nmax)
    real psnd(nmax), rhsnd(nmax)
    real tv0, tv1
    real z, an, ugnd, ugnd2

    real zfind, p000, t000, q000
    external zfind
    integer ktop
    real tmptr,zsfc, qfac
    real newp000
    logical :: rescalep
      real :: pzsfc0,pzsfc, tksfc

! Error check


!------------------------------------------------------------------------------
! First read in the surface pressure, temperature, and mixing ratio
! Format is z,th,qv,u,v because of stretched grid
!     Sounding input file format: ascii
!     z (meters), theta (K), qv (g/kg),  u (m/s), v (m/s)
!------------------------------------------------------------------------------

    CALL GET_VARIABLE(gd, 'NX', nx)
    CALL GET_VARIABLE(gd, 'NY', ny)

    CALL FILE_MESSAGE('READING SOUNDING FROM: '//sndfile)

    open(unit=17,file=sndfile,status='old',form='formatted')
    rewind(17)

! Read in sfc pressure, theta, qv

    read(17,*) p000, t000, q000

    IF ( p000 < 0.0 ) THEN
    
    read(17,*) newp000
    p000 = Abs(p000)
    rescalep = .true.
!    write(luno,*) 'rescaling for new sfc pressure ',newp000
!    write(luno,*) 'WARNING: This does not work correctly'
    STOP
    IF( newp000 .lt. 1.0e4 ) newp000 = newp000*100.
    
    ELSE
    
    rescalep = .false.
    
    ENDIF
    
    
!    write(luno,*) p000, t000, q000

   qsfcwk = qsfc ! save input value from namelist
   IF ( wtype == 102 .or. wtype == 101 ) THEN
     t000 = tsfc
     p000 = psfc
   ELSE
     tsfc = t000
     psfc = p000
   ENDIF

   IF ( wtype == 102 ) THEN
     q000 = qsfc
   ELSE
     qsfc = q000
   ENDIF
   
   qfac = 1.0
   IF ( qsfc < 0 ) THEN
      qvorrh = 1
   ELSE
      qvorrh = 0
   ENDIF
   IF( psfc .lt. 1.0e4 ) psfc = psfc*100.
   IF( qsfc .gt. 0.05  ) qfac = 0.001
   
   qsfc = qfac*qsfc


! SET some defaults

   IF( .not. SET_VARIABLE(gd,'PSFC',       psfc) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

! Read in sounding values

    nsnd = 0
    DO k = 1,nmax
       IF ( qvorrh == 0 ) THEN
       read(17,*,iostat=ios) zsnd(k),tsnd(k),qvsnd(k),usnd(k),vsnd(k)
       qvsnd(k) = qfac*qvsnd(k)
       ELSE
       read(17,*,iostat=ios) zsnd(k),tsnd(k),rhsnd(k),usnd(k),vsnd(k), psnd(k)
       temp = tsnd(k) + 273.16
       tsnd(k) = tsnd(k)*(psfc/psnd(k))**rcp
        
        pres  = psnd(k) ! psfc*pz(k)**(cp/rd)
     !   pztmp = (psnd(k)/psfc)**rcp
        qvs   = 380.*exp(17.27*(temp-273.16) / (temp- 36.)) / pres
        qvsnd(k) = 0.01*rhsnd(k) * qvs
       
       write(0,*) zsnd(k),tsnd(k),qvsnd(k),usnd(k),vsnd(k), psnd(k)
       ENDIF
       IF( ios < 0 ) EXIT
       nsnd = nsnd + 1
    ENDDO
    
    close(17)

!#ifdef MPI
!    if (kzbeg .eq. nzbeg)  zsfc = zsnd(1)
!#else
    zsfc = zsnd(1)
!#endif

!   IF( .not. SET_VARIABLE(gd,'HGT',       zsfc) ) write(6,*) 'INIT_GRID:  Problem setting HGT'

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .le. nsnd ) kze = nsnd-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nsnd
!#endif
      zsnd(k) = zsnd(k) - zsfc
    ENDDO

    zsnd(nsnd) = zsnd(nsnd) + 0.1

! Check to make sure the sounding levels span the computational grid

!#ifdef MPI
!    if (kzbeg .eq. nzbeg .AND. zsnd(1) .gt. zc(1) ) THEN
!#else
    IF( zsnd(1) .gt. zc(1) ) THEN
!#endif
      write(luno,*) 'ZMIN of SOUNDING   = ',zsnd(1)
      write(luno,*) 'ZMIN of MODEL GRID = ',zc(1)
      write(luno,*) 'ZMIN OF INPUT SOUNDING IS > ZMIN OF GRID'
      write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
      stop
    ENDIF

!#ifdef MPI
    IF( zsnd(nsnd) .lt. zc(nzend-1) ) THEN
!#else
!    IF( zsnd(nsnd) .lt. zc(nz-1) ) THEN
!#endif
      write(luno,*) 'ZMAX of SOUND      = ',zsnd(nsnd)
      write(luno,*) 'ZMAX of MODEL GRID = ',zc(nzend-1)
      write(luno,*) 'ZMAX OF INPUT SOUNDING IS < ZMAX OF MODEL GRID'
      write(luno,*) 'Extrapolating sounding with isothermal layer....'
!      write(luno,*) 'STOPPING RUN....PLEASE CORRECT INPUT FILE'
!      luno = FILE_CLOSE('ERROR IN SND2 !!!!!!!!! - STOPING EXECUTION!')
!      stop
    ENDIF

! Now interpolate sounding to grid

    write(luno,*)
    write(luno,*) 'INTERPOLATING SOUNDING TO GRID '

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
      z0 = zc(k)
      IF ( z0 .lt. zsnd(nsnd) ) THEN
        u1d(k) = zfind(zsnd, usnd,  nsnd, z0)
        v1d(k) = zfind(zsnd, vsnd,  nsnd, z0)
        tz(k) = zfind(zsnd, tsnd,  nsnd, z0)
        qz(k) = zfind(zsnd, qvsnd, nsnd, z0)
        ktop = k


        IF ( k .eq. 1 ) THEN
          tv0   = tz(1) * (1.0 + 0.601*qz(1))
          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
          temk(k) = pz(k)*tz(k)
        ELSE
          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k)) ! dze = 1/dz

         pres  = 100000.*pz(k)**3.508 ! psfc*pz(k)**(cp/rd)
         qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
         rh(k) = 0
         IF ( qvs > 0.0 .or. qz(k) > 0.0 ) THEN
           rh(k) = qz(k)/qvs
          ! IF ( my_rank == 0 ) write(0,*) 'k,rh = ',k,rh(k),qz(k),qvs,pres,psfc*pz(k)**(cp/rd),pz(k)
           IF ( rh(k) > rhmaxin ) THEN
             DO i = 1,3
               qz(k) = rhmaxin * qvs
               tv1   = tz(k) * (1.0 + 0.601*qz(k))
               pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
               pres  = 100000.*pz(k)**3.508 ! psfc*pz(k)**(cp/rd)
               qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
             ENDDO
           ENDIF
           
           
         ELSE
           rh(k) = 0
         ENDIF

          tv0   = tv1
          tmptr = pz(k)*tz(k)
          temk(k) = pz(k)*tz(k)
        ENDIF
      ELSE
        u1d(k) = u1d(ktop)
        v1d(k) = v1d(ktop)
        qz(k) = qz(ktop)
        tz(k) = tz(ktop)*(Exp((z0-zc(ktop) - zsnd(1))*g/cp/tmptr ) )

          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
          temk(k) = pz(k)*tz(k)

      ENDIF

  

!      IF( qz(k) .gt. 0.05 ) qz(k) = qz(k) / 1000.

    ENDDO
    
    IF ( rescalep ) THEN
      
      pzsfc0 = (psfc/1.e5)**rcp ! original surface Exner Pi
      tksfc = pzsfc0 * tsfc
      write(luno,*) 'old tsfc = ',tsfc
      ! get new surface Ex Pi and new theta using (real) surface temperature
      pzsfc = (newp000/1.e5)**rcp
      tsfc = tksfc/pzsfc
      write(luno,*) 'new tsfc = ',tsfc
      psfc = newp000
      t000 = tsfc
      
      tsnd(1) = tsfc
      
      DO k = 1,nzend-1
       write(luno,*) 'k, old pz = ',k,pz(k),tz(k)
       pz(k) = pz(k) - pzsfc0 + pzsfc
       tz(k) = temk(k)/pz(k)
       write(luno,*) 'k, temp pz = ',k,pz(k),tz(k)
      ENDDO
      
      DO k = 1,nzend-1
        
        IF ( k .eq. 1 ) THEN
           
           tv0   = tz(1) * (1.0 + 0.601*qz(1))
           pz(1) = (newp000/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
           tz(k) = temk(k)/pz(k)
!          tv0   = tz(1) * (1.0 + 0.601*qz(1))
!          pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)
!          temk(k) = pz(k)*tz(k)
        ELSE
          tv1   = tz(k) * (1.0 + 0.601*qz(k))
          pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
          tv0   = tv1
!          tmptr = pz(k)*tz(k)
!          temk(k) = pz(k)*tz(k)
        ENDIF
        
        
        
      ENDDO

   IF( .not. SET_VARIABLE(gd,'PSFC',     newp000) ) write(6,*) 'INIT_GRID:  Problem setting PSFC'
   IF( .not. SET_VARIABLE(gd,'TSFC',       tsfc) ) write(6,*) 'INIT_GRID:  Problem setting TSFC'
   IF( .not. SET_VARIABLE(gd,'QSFC',       qsfc) ) write(6,*) 'INIT_GRID:  Problem setting QSFC'

    ENDIF

! overwrite theta
   IF ( wtype == 100 ) THEN

! BASIC WK82 Sounding
!                         zbl, p_bl, thbl
!#ifdef MPI
!   kzb = -ng+1
!   kze = ktile+ng
!   if (kzbeg .eq. nzbeg) kzb = 1
!   if (kzend .eq. nzend) kze = kzend-kzbeg
!   
!   DO k = kzb,kze
!#else
   DO k = 1,nzend-1
!#endif
    z = zc(k)
    IF( z .le. ztr ) THEN
     IF ( zbl > 0. .and. ( thbl > tsfc .or. dthdzbl > 0.0 ) ) THEN
       IF ( dthdzbl > 0.0 ) THEN
        ! let dthdzbl supercede thbl
        thbl = tsfc + zbl*dthdzbl
       ENDIF
       IF ( z <= zbl ) THEN
         zfac  = ( z / zbl ) ! linear increase from surface to zbl
         tz(k) = tsfc + ( thbl - tsfc ) * zfac
       ELSE
         zfac  = ( (z - zbl) / (ztr - zbl) ) ** wkexponent
         tz(k) = thbl + ( thtr - thbl ) * zfac
       ENDIF
     ELSE
     zfac  = ( z / ztr ) ** wkexponent
     tz(k) = tsfc + ( thtr - tsfc ) * zfac
     ENDIF
    ELSE
     tz(k) = thtr * exp (g * (z - ztr) / (cp * ttr) )
    ENDIF

!    write(0,*) 'k,zc = ',k,zc(k),zfac,tz(k),rh(k)

   ENDDO 
   
   
   ENDIF

   IF ( wtype == 101 ) THEN
     ! use subroutine to set thermo
     tsnd(1) = tsfc
     i = 0
     call SNDWKTHERMO(tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   i)

   ELSEIF ( wtype == 102 ) THEN
     ! use subroutine to set thermo and moisture
     i = 1
     qvsnd(1) = qsfc
     tsnd(1) = tsfc
     call SNDWKTHERMO(tz,pz,qz,dzc,dze,zc,nzend,ng,psfc,tsfc,qsfc,rhmax,rhmax2,  &
                   zrhmax2,zrhdel,wtype,Us,Uz,Usl,Usmv,Ubase,ntr,u1d,v1d,z0,z1,dudz0,rotdeg, &
                   i)


   ENDIF

! PI INIT via vertical integration of hydrostatic equation
! Here we know Qv, so we can integrate using virtual temperature

!    tv0   = tz(1) * (1.0 + 0.601*qz(1))
!    pz(1) = (psfc/1.e5)**rcp - 0.5*g/(dzc(1)*tv0*cp)

    DO k = 2,nzend-1

!     tv1   = tz(k) * (1.0 + 0.601*qz(k))
!     pz(k) = pz(k-1) - 2.0*g / ((tv0+tv1)*cp*dze(k))
!     tv0   = tv1

    ENDDO


! CREATE 1D initial wind profile

   IF( wtype .eq. 10 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
     IF( zc(k) .le. Uz ) u1d(k) = Us * ( zc(k)/Uz )
     IF( zc(k) .gt. Uz ) u1d(k) = Us 
     v1d(k) = 0.0

    ENDDO

   ENDIF

! CREATE 2D Weisman Profile

   IF( wtype .eq. 20 ) THEN

      ugnd = Us / pii

!#ifdef MPI
!    kzb = 1
!    kze = ktile
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
    DO k = 1,nzend-1
!#endif
       z = zc(k)
       IF( z .le. Uz ) THEN
        an = degtorad*(180.0 * (1.0 - z / Uz) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. Uz ) THEN
        u1d(k) = 2*ugnd
        v1d(k) = 0.0
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- hockey stick (1/4 circle to 2km, then linear shear
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 30 ) THEN

      ugnd = (Us/3.0) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + Us*(2./3.)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 6000.  ) THEN
        u1d(k) = ugnd + Us*(2./3.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF



! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC)
   IF( wtype .eq. 50 ) THEN

      ugnd = (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = degtorad*(90.0 * (2.0 - z / 2000.) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. 2000. .and. z .le. 7000. ) THEN
        u1d(k) = ugnd + Us*(5./7.)*(z - 2000.)/5000.0
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. 7000.  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2.5km, then linear shear extended to 7km
! recreation of shear in WRF quarter circle supercell
   IF( wtype .eq. 60 ) THEN
   
    chgt = 2500.0
    shgt = 7000.0

      ugnd = (chgt*Us/shgt) / (pii/2)  ! radius of quarter circle
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. chgt ) THEN
        an = degtorad*(90.0 * (2.0 - z / chgt) )
        u1d(k) = ugnd * (1.0 + cos(an))
        v1d(k) = ugnd * sin(an)
       ENDIF
       IF( z .gt. chgt .and. z .le. shgt ) THEN
        u1d(k) = ugnd + Us*((shgt-chgt)/shgt)*(z - chgt)/(shgt - chgt)
        v1d(k) = ugnd
       ENDIF
       IF( z .gt. shgt  ) THEN
        u1d(k) = ugnd + Us*(5./7.)
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! CREATE 2D Weisman Profile -- candy cane (1/4 circle to 2km, then linear shear extended to 7km
! recreation of shear in Weisman, Gilmore, and Wicker (1998 SLSC). Code to match CM1
! Hodograph length is about 35 m/s
   IF( wtype .eq. 70 ) THEN

      ugnd = 7.0 ! (2.0*Us/7.0) / (pii/2)  ! radius of quarter circle
      ugnd2 = 31.0
      
      
!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else      
      DO k = 1,nzend
!#endif
       z = zc(k)
       IF( z .le. 2000.0 ) THEN
        an = 90.*degtorad*( z / 2000.) 
        u1d(k) = ugnd * (1.0 - cos(an))
        v1d(k) = ugnd * sin(an)
       ELSEIF( z .gt. 2000. .and. z .le. 6000. ) THEN
        u1d(k) = ugnd + (ugnd2-ugnd)*(z - 2000.)/4000.0
        v1d(k) = ugnd
       ELSE
        u1d(k) = ugnd2
        v1d(k) = ugnd
       ENDIF
      ENDDO
      
   ENDIF

! Flip winds:
   IF( wtype .eq. 1000 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        u1d(k) = -u1d(k)
        v1d(k) = -v1d(k)
      ENDDO

   ENDIF

! Swap winds:
   IF( wtype .eq. 1001 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        z = u1d(k)
        u1d(k) = v1d(k)
        v1d(k) = z
      ENDDO

   ENDIF

! Rotate winds:
   IF( wtype .eq. 1002 ) THEN

!#ifdef MPI
!    kzb = -ng+1
!    kze = ktile+ng
!    if (kzbeg .eq. nzbeg) kzb = 1
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!    DO k = kzb,kze
!#else
      DO k = 1,nzend-1
!#endif
        z = u1d(k)
        u1d(k) = -v1d(k)
        v1d(k) = z
      ENDDO

   ENDIF

! Write sounding out to file

    IF ( ny .le. 2 ) THEN
      v1d(:) = 0.0
    ELSEIF ( nx .le. 2 ) THEN
      u1d(:) = 0.0
    ENDIF
    

  if(my_rank == 0) then
! write out first two lines of sounding: surface values and first input sounding level (taken as surface)
    open(21,file=trim(outsoundfile), status='unknown')
    rewind(21)
!    write(21,*) 'MODEL SND'
    write(21,'(1x,2(f9.2,2x),f9.4)') psfc, t000, qsfc*1000.
!    DO k = 1,nz-1
     k = 1
     write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),psfc
!     write(21,'(1x, 6(f10.2,2x))') 0.0, tsnd(k),1000.*qvsnd(k),usnd(k),vsnd(k),p000
!    ENDDO
    close(21)

  endif
   END SUBROUTINE SND2X

