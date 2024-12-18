!-------------------------------------------------------------------------------
! 
!
!         <<<<<<<<<<<<<<<  GRID_DEFINE_FROM_LIST  >>>>>>>>>>>>>>>>>>>>  
!  
!
! Written by Lou Wicker, Nov/Dec 2005
!  
! I tend to hate input files, so here is a way to define your model variables
!   internally for initialization.  Files get corrupted, versions change, etc.
!
! Its stupid, but the easiest way is to write the information to a file
!    and let the GRID_DEFINE_FROM_FILE read the file....
! --- UPDATE (at some point...) Info is written to a string array instead of a file to avoid the
!       I/O, which could be significant on MPI. The ugly part is having to tick the array index by
!       1 each time, thus all the necessary lines with i = i + 1
!
! GRID_DEFINE_FROM_FILE is called at the bottom 
!-------------------------------------------------------------------------------

 SUBROUTINE GRID_DEFINE_FROM_LIST(gd, microphys, ipconc, ichaff, inucopt, ienkf, isfcphys)

  USE STRING_MODULE
  USE INDEX_MODULE, only: neelec2d,nch,ntakpd,ntakid,ntakit,ntakrd,imixice,icespheres,   &
                          frozendrops,nprecip,iraintypes,nraintypes,iashtypes, &
                          inetchargetend,ioutput_hwind3d,ioutput_temC3d,ioutput_pres3d,        &
                          ioutput_thv3d,ioutput_workshop,ioutput_xtrachgsep,ioutput_flshr8km,  &
                          ioutput_ssw,ioutput_ssi,ioutput_rh,ioutput_ssmx,ioutput_cwdia,ioutput_cnu,       &
                          ioutput_snowstuff, ioutput_takrates, ioutput_mltshedsizerates, &
                          ioutput_sedmeltstuff, ioutput_vzf, ioutput_icedensity, ioutput_tkediss,  &
                          ioutput_dbzsedchange
  USE MICRO_MODULE, only: idoniconly,irenuc,icenucopt,takcwnucopt,ibinhlmlr,ibinhmlr,ccnuf,takkfmax, &
                          in_freeze_rain_first,ac_opt,iturbenhance
  USE ELEC_MODULE, only: ieopt
  USE FORCE_MODULE, only: iccnufforce
  USE PARAM_MODULE, only: tke_type

  implicit none
  
  type(GRID)                      :: gd
  character(LEN=*), INTENT(INOUT) :: microphys
  integer                         :: ipconc
  integer                         :: ichaff
  integer                         :: inucopt
  integer, optional               :: ienkf
  integer, optional               :: isfcphys
  logical                         :: elec, pwf_flag, hl_flag
  
  integer                         :: isfc

  integer, parameter              :: nattr = 25 ! when Conventions is turned on
!  integer, parameter              :: nattr = 24
  integer                         :: nvar  = 86
  integer                         :: ikf,i,it,llen,n,m

  character(LEN=10), parameter :: dummy_char = '----------'
  character(LEN=14), parameter :: internal_filename = 'modelfield.tmp'
 
  character(len=2) :: tmpnum
  character(len=1) :: tmpdigit
  character(len=name_length) :: tmpname
  character(len=40) :: tmpdescript
  
  ! flags for turning on appropriage fields
  
  integer :: hasqv = 1, hasqc = 0, hasqr = 0, hasqi = 0, hasqs = 0, hasqh = 0, hasqhl = 0
  integer :: hasqsw = 0, hasqhw = 0, hasqhlw = 0
  integer :: hasccn = 0, hascin = 0, hasccna = 0, hascina = 0
  integer :: hasnc = 0, hasnr = 0, hasni = 0, hasns = 0, hasnh = 0, hasnhl = 0
  integer :: haszc = 0, haszr = 0, haszi = 0, haszs = 0, haszh = 0, haszhl = 0
  
  microphys = UCASE(microphys)
  elec      = .false.
  isfc  = 0
  IF ( present(isfcphys) ) isfc = isfcphys

  llen = index(microphys,' ')-1
  
  IF ( present(ienkf) ) THEN
    ikf = ienkf
  ELSE
    ikf = 0
  ENDIF
    
! hack to force ikf  
!  ikf = 1

  IF ( microphys(1:7) .eq. 'KESSLER' ) THEN
  
    hasqc = 1
    hasqr = 1

  ELSEIF ( microphys(1:3) .eq. 'SCG' ) THEN

    hasqc = 1
    hasqr = 1

  ELSEIF ( microphys(1:7) .eq. 'WARMLFO' ) THEN

    hasqc = 1
    hasqr = 1

  ELSEIF ( microphys(1:3) .eq. 'LFO' ) THEN

    hasqc = 1
    hasqr = 1
    hasqi = 1
    hasqs = 1
    hasqh = 1

  ELSEIF ( microphys(1:4) .eq. 'THOM' ) THEN

    hasqc = 1
    hasqr = 1
    hasqi = 1
    hasqs = 1
    hasqh = 1
    
    hasnr = 1
    hasni = 1
    
  ELSEIF ( microphys(1:4) .eq. 'MORR' ) THEN

    hasqc = 1
    hasqr = 1
    hasqi = 1
    hasqs = 1
    hasqh = 1
    
    hasnr = 1
    hasni = 1
    hasns = 1
    hasnh = 1
    

  ELSEIF ( microphys(1:8) .eq. 'ZIEGELEC' .or. microphys .eq. 'ZIEGHELEC' .or.   &
           microphys(1:5) .eq. 'ZIEGE'    .or. microphys(1:6) .eq. 'ZIEGHE' ) THEN

    elec = .true.

  ELSEIF ( microphys(1:4) .eq. 'ZIEG' ) THEN

  ELSEIF ( microphys(1:8) .eq. 'WARMZIEG' ) THEN

  ELSEIF ( microphys(1:4) .eq. 'ZVDE' .or. microphys(1:5) .eq. 'ZVDHE' .or.   &
           microphys(1:5) .eq. 'ZVDME' .or. microphys(1:6) .eq. 'ZVDHVE' .or. &
           microphys(1:6) .eq. 'ZVDHME' .or. microphys(1:6) .eq. 'ZVDMHE' ) THEN
  
   elec = .true.

  ELSEIF ( microphys(1:4) .eq. 'ZMRE' .or. microphys(1:5) .eq. 'ZMRHE' .or.   &
           microphys(1:5) .eq. 'ZMRME' .or. microphys(1:6) .eq. 'ZMRHVE' .or. &
           microphys(1:6) .eq. 'ZMRHME' .or. microphys(1:6) .eq. 'ZMRMHE' ) THEN
  
   elec = .true.

  ELSEIF ( microphys(1:3) .eq. 'ZVD' ) THEN

  ELSEIF ( microphys(1:3) .eq. 'ZMR' ) THEN

  ELSEIF ( microphys(1:4) .eq. 'ICE3' ) THEN


       
  ELSEIF ( microphys(1:8) .eq. 'ICE10DME' ) THEN

    elec = .true.

  ELSEIF ( microphys(1:7) .eq. 'ICE10DM' ) THEN


  ELSEIF ( microphys(1:6) .eq. 'ICE10E' ) THEN

    elec = .true.

  ELSEIF ( microphys(1:5) .eq. 'ICE10' ) THEN

  ELSEIF ( microphys(1:3) .eq. 'DRY' ) THEN
  
  ELSEIF ( microphys(1:2) .eq. 'MY') THEN

  ELSEIF ( microphys(1:3) .eq. 'HCM') THEN

  ELSEIF ( microphys(1:3) .eq. 'TAK') THEN
     IF ( microphys(llen:llen) == 'E' .or.  microphys(llen-1:llen-1) == 'E' ) elec = .true.
  ELSE    

    write(6,*) 'GRID_DEFINE_FROM_LIST:  UNRECOGNIZED OPTION: ', microphys, 'FOR MICROPHYSICS SCHEME  STOP!'
    write(0,*) 'GRID_DEFINE_FROM_LIST:  UNRECOGNIZED OPTION: ', microphys, 'FOR MICROPHYSICS SCHEME  STOP!'
    STOP

  ENDIF

  pwf_flag = .false.
  IF ( microphys(1:4) .eq. 'ZVDM' .or. microphys(1:5) .eq. 'ZVDHM' .or. microphys(1:5) .eq. 'ZIEGM' &
        .or. microphys(1:6) .eq. 'ZIEGHM' .or. microphys(1:4) .eq. 'ZMRM' .or. microphys(1:5) .eq. 'ZMRHM') THEN
      pwf_flag = .true.
  ENDIF
  
  hl_flag = .false. 
  
  IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:4) .eq. 'ZIEGH' .or. microphys(1:4) .eq. 'ZMRH') THEN
      hl_flag = .true.
  ENDIF

!----------------------------------------------------------------------------------
! Open the file for writing...

!  open(51,file=internal_filename,status='unknown',form = 'formatted' )
!  rewind(51)

!  open(51,status='scratch',form = 'formatted' )
  numlines = 0
  i = 0
!----------------------------------------------------------------------------------------------------
! Write out the number of attrs and variables

 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
   
 i = i + 1
 write(f(i),*) nattr
  IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_LIST:  NATTR = ', nattr
   
 i = i + 1
 write(f(i),*) nvar
  IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_LIST:  NVAR  = ', nvar
   
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file

!----------------------------------------------------------------------------------------------------
! Write out the grid attributes

101 format(1x,a,2x,a,2x,i3,2x,g15.5)

 i = i + 1
 write(f(i),101) 'NX                  ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NY                  ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NZ                  ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NG                  ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NXEND               ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NYEND               ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'NZEND               ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'PREFIX_NAME         ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'MEMBER              ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'SIZE_OF_ENSEMBLE    ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'CREATION_DATE       ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'CREATION_TIME       ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'MEMBER_NAME         ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'OUTPUT_FILE_NAME    ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'COARDS              ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'Conventions         ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'MICROPHYS           ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'V5DFIELDS           ', 'str', -1, -1.0
 i = i + 1
 write(f(i),101) 'YEAR                ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'MONTH               ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'DAY                 ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'HOUR                ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'MINUTE              ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'SECOND              ', 'int', -1, -1.0
 i = i + 1
 write(f(i),101) 'TILE_INDEX          ', 'int', -1, -1.0
! i = i + 1
! write(f(i),101) 'CREATOR             ', 'str', -1, -1.0
! i = i + 1
! write(f(i),101) 'NVAR                ', 'int', -1, -1.0
 
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file
 i = i + 1
 write(f(i),*) dummy_char                        ! this is just to skip a comment line in the file

!----------------------------------------------------------------------------------
! Write out variable information

2 format(1x,a10,2x,i1,2x,a5,2x,8(i2,2x),a15,2x,a)
4 format(1x,a,2x,i3,2x,a,2x,8(i2,2x),a,2x,a)

! pd (posdef):  0 = not positive definite
!               1 = is positive definite
!
! dyntype :  1 = do not scale by density (i.e., convert quantity to 'mixing ratio')
!            0 =  convert quantity to 'mixing ratio' for advection/mixing
!
! phytype : 0 = apply turb. mixing
!           1 or anything nonzero : do not apply turb mixing
!
! buotype :  0 = not used for buoyancy calc
!            2 = is used for buoy calc (vapor mixing ratio)
!            3 = is used for buoy calc (condensed water mixing ratio)
!               name         dim  type   td  stag   pd  dyn phy buo  unit              description

! INTEGER CONSTANTS

 i = i + 1
 write(f(i),2) 'NX        ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON U-GRID IN X         '
 i = i + 1
 write(f(i),2) 'NY        ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON V-GRID IN Y         '
 i = i + 1
 write(f(i),2) 'NZ        ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON W-GRID IN Z         '
 i = i + 1
 write(f(i),2) 'NG        ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF COMPUTATIONAL BUFFER ZONES    '
 i = i + 1
 write(f(i),2) 'NXEND     ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON FULL U-GRID IN X    '
 i = i + 1
 write(f(i),2) 'NYEND     ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON FULL V-GRID IN Y    '
 i = i + 1
 write(f(i),2) 'NZEND     ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF POINTS ON FULL W-GRID IN Z    '
 i = i + 1
 write(f(i),2) 'NSCALAR   ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF SCALAR FIELDS                 '
 i = i + 1
 write(f(i),2) 'NSMALL    ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'NUMBER OF SMALL TIME STEPS              '
 i = i + 1
 write(f(i),2) 'BCX       ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'BOUNDARY CONDITION FOR X                '
 i = i + 1
 write(f(i),2) 'BCY       ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'BOUNDARY CONDITION FOR Y                '
 i = i + 1
 write(f(i),2) 'BCZ       ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'BOUNDARY CONDITION FOR Z                '
 i = i + 1
 write(f(i),2) 'N_INNER   ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', '# POINTS FOR HOR STRETCHED GRID         '
 i = i + 1
 write(f(i),2) 'N_BNDLYR  ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', '# BLR POINTS GRID WHERE DZBOT CONSTANT  '
 i = i + 1
 write(f(i),2) 'HOLE_FILL ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'HOLE FILLING PARAMETER FOR MICROPHYSICS '
 i = i + 1
 write(f(i),2) 'AUTOCONV  ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'AUTOCONVERSION FLAG                     '
 i = i + 1
 write(f(i),2) 'TIME      ', 0, 'icnst', 1, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'LOCAL MODEL TIME IN SECONDS FROM START  '
 i = i + 1
 write(f(i),2) 'TRESTART  ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'DT IN SECONDS BETWEEN RESTART DUMPS     '
 i = i + 1
 write(f(i),2) 'THISTORY  ', 0, 'icnst', 1, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'DT IN SECONDS BETWEEN HISTORY DUMPS     '
 i = i + 1
 write(f(i),2) 'TVIS5D    ', 0, 'icnst', 1, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'DT IN SECONDS BETWEEN VIS5D DUMPS       '
 i = i + 1
 write(f(i),2) 'TPRINT    ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'DT IN SECONDS BETWEEN PRINT OUTPUT      '
 i = i + 1
 write(f(i),2) 'TSTAT     ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'DT IN SECONDS BETWEEN STAT OUTPUT       '
 i = i + 1
 write(f(i),2) 'TIME_STOP ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'LOCAL TIME IN SECS MODEL INTEGRATED TO  '
 i = i + 1
 write(f(i),2) 'IPCONC    ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CONCENTRATION  FLAG                     '
  i = i + 1
 write(f(i),2) 'INUCOPT   ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'FLAG FOR ACTIVATED NUCLEI PREDICTION    '
 i = i + 1
 write(f(i),2) 'IPELEC    ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'ELECTRIFICATION FLAG                    '
 i = i + 1
 write(f(i),2) 'LGTSEED   ', 0, 'icnst', 1, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'Random number seed for lightning        '
! i = i + 1
! write(f(i),2) 'track     ', 0, 'icnst', 1, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'Y-TRANSLATION SPEED OF GRID             '
 i = i + 1
 write(f(i),2) 'SST       ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'kelvin         ', 'Sea Surface Temperature                 '
 IF ( ichaff .ge. 0 ) THEN 
! allow ichaff < 0 as flag for backwards compatibility (avoids error mssgs in reading/restarting)
 i = i + 1
 write(f(i),2) 'ICHAFF    ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CHAFF FLAG                              '
 ENDIF
 i = i + 1
 write(f(i),2) 'IENKF     ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'ENKF FLAG                               '
 i = i + 1
 write(f(i),2) 'ISFCPHYS  ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SURFACE PHYSICS FLAG                    '
 i = i + 1
 write(f(i),2) 'IMURAIN   ', 0, 'icnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'RAIN SECOND SHAPE PARAMETER             '

! REAL CONSTANTS
 i = i + 1
 write(f(i),2) 'DT        ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'seconds        ', 'TIME STEP                               '
 i = i + 1
 write(f(i),2) 'DX        ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'MEAN X-GRID SPACING                     '
 i = i + 1
 write(f(i),2) 'DY        ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'MEAN Y-GRID SPACING                     '
 i = i + 1
 write(f(i),2) 'DZ        ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'MEAN Z-GRID SPACING                     '
 i = i + 1
 write(f(i),2) 'DX_STRETCH', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'STRETCH PARAMETER FOR X-GRID            '
 i = i + 1
 write(f(i),2) 'DY_STRETCH', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'STRETCH PARAMETER FOR Y-GRID            '
 i = i + 1
 write(f(i),2) 'DZ_STRETCH', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'STRETCH PARAMETER FOR Z-GRID            '
 i = i + 1
 write(f(i),2) 'DX_CORNER ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'DX OF SW CORNER DUE TO BOX MOTION       '
 i = i + 1
 write(f(i),2) 'DY_CORNER ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'DY OF SW CORNER DUE TO BOX MOTION       '
 i = i + 1
 write(f(i),2) 'XDOMAIN   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'E-W Domain extent                       '
 i = i + 1
 write(f(i),2) 'YDOMAIN   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'S-N Domain extent                       '
 i = i + 1
 write(f(i),2) 'ZDOMAIN   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'Vertical Domain extent                  '
 i = i + 1
 write(f(i),2) 'CORIOLIS  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'second-1       ', 'CORIOLIS VALUE                          '
 i = i + 1
 write(f(i),2) 'CDM       ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'DRAG COEFFICENT IF BCZ == 1             '
 i = i + 1
 write(f(i),2) 'PSFC      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'hPa            ', 'SURFACE PRESSURE                        '
 i = i + 1
 write(f(i),2) 'TSFC      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kelvin         ', 'SURFACE TEMPERATURE AT GROUND           '
 i = i + 1
 write(f(i),2) 'QSFC      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg kg-1        ', 'SURFACE MIXING RATIO AT GROUND          '
 i = i + 1
 write(f(i),2) 'XBAR      ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters         ', 'USED TO TRACK STORM MOTION              '
 i = i + 1
 write(f(i),2) 'YBAR      ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters         ', 'USED TO TRACK STORM MOTION              '
 i = i + 1
 write(f(i),2) 'UGRID0    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'X-TRANSLATION SPEED OF GRID AT T=0      '
 i = i + 1
 write(f(i),2) 'VGRID0    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'Y-TRANSLATION SPEED OF GRID AT T=0      '
 i = i + 1
 write(f(i),2) 'UGRID     ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'X-TRANSLATION SPEED OF GRID             '
 i = i + 1
 write(f(i),2) 'VGRID     ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'Y-TRANSLATION SPEED OF GRID             '
 i = i + 1
 write(f(i),2) 'XG_POS    ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters         ', 'LOCATION OF SW BOX CORNER               '
 i = i + 1
 write(f(i),2) 'YG_POS    ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'meters         ', 'LOCATION OF SW BOX CORNER               '
 i = i + 1
 write(f(i),2) 'RHO_QR    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF RAINWATER                    '
 i = i + 1
 write(f(i),2) 'CNOR      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF RAINWATER            '
 i = i + 1
 write(f(i),2) 'RHO_QS    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF SNOW                         '
 i = i + 1
 write(f(i),2) 'CNOS      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF SNOW                 '
 i = i + 1
 write(f(i),2) 'RHO_QH    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF HAIL                         '
 i = i + 1
 write(f(i),2) 'RHO_QHL   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF LARGE HAIL                   '
 i = i + 1
 write(f(i),2) 'CNOH      ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF HAIL                 '
 i = i + 1
 write(f(i),2) 'CNOHL     ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF LARGE HAIL           '
 i = i + 1
 write(f(i),2) 'RHO_QR_10 ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF RAINWATER                    '
 i = i + 1
 write(f(i),2) 'RHO_QS_10 ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF SNOW                         '
 i = i + 1
 write(f(i),2) 'RHO_QGL_10', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF GRAUPEL-L                    '
 i = i + 1
 write(f(i),2) 'RHO_QGM_10', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF GRAUPEL-M                    '
 i = i + 1
 write(f(i),2) 'RHO_QGH_10', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF GRAUPEL-H                    '
 i = i + 1
 write(f(i),2) 'RHO_QF_10 ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF FROZEN DROPS                 '
 i = i + 1
 write(f(i),2) 'RHO_QH_10 ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF HAIL                         '
 i = i + 1
 write(f(i),2) 'RHO_QHL_10', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF LARGE HAIL                   '
 i = i + 1
 write(f(i),2) 'CNOR_10   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF RAINWATER            '
 i = i + 1
 write(f(i),2) 'CNOS_10   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF SNOW                 '
 i = i + 1
 write(f(i),2) 'CNOGL_10  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF GRAUPEL-L            '
 i = i + 1
 write(f(i),2) 'CNOGM_10  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF GRAUPEL-M            '
 i = i + 1
 write(f(i),2) 'CNOGH_10  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF GRAUPEL-H            '
 i = i + 1
 write(f(i),2) 'CNOF_10   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF FROZEN DROPS         '
 i = i + 1
 write(f(i),2) 'CNOH_10   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF HAIL                 '
 i = i + 1
 write(f(i),2) 'CNOHL_10  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SLOPE INTERCEPT OF LARGE HAIL           '
 i = i + 1
 write(f(i),2) 'QCMINCON  ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg kg-1        ', 'KESSLER QC THRESHOLD FOR AUTOCONVERSION '
 i = i + 1
 write(f(i),2) 'CWDIAP    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'MEAN DIAMETER FOR Ferrier/Berry AC      '
 i = i + 1
 write(f(i),2) 'CWDISP    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'MEAN DISP FOR Ferrier/Berry AC          '
 i = i + 1
 write(f(i),2) 'CCN       ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CCN CONCENTRATION FOR Ferrier/Berry AC  '
 i = i + 1
 write(f(i),2) 'CCNUF     ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CCN CONCENTRATION ULTRA FINE            '
 i = i + 1
 write(f(i),2) 'CCNAC     ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CN_AC CONCENTRATION                     '
 i = i + 1
 write(f(i),2) 'CCNNU     ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CN_NU CONCENTRATION                     '
 i = i + 1
 write(f(i),2) 'CCNCO     ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'CN_CO CONCENTRATION                     '
 i = i + 1
 write(f(i),2) 'LAT       ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'degrees        ', 'STARTING LATITUDE OF SW CORNER OF GRID  '
 i = i + 1
 write(f(i),2) 'LON       ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'degrees        ', 'STARTING LONGITUDE OF SW CORNER OF GRID '
 i = i + 1
 write(f(i),2) 'HGT       ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'HEIGHT OF THE GRID                      '
 i = i + 1
 write(f(i),2) 'ALPHAH    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF GRAUPEL              '
 i = i + 1
 write(f(i),2) 'ALPHAHL   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF HAIL                 '
 i = i + 1
 write(f(i),2) 'ALPHAI    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF ICE CRYSTALS         '
 i = i + 1
 write(f(i),2) 'ALPHAS    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF SNOW                 '
 i = i + 1
 write(f(i),2) 'EHSLFO0   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'FACTOR in graupel-snow collection       '
! i = i + 1
! write(f(i),2) 'RAIN_TOT  ', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'kg m-2         ', 'TOTAL ACCUMULATED RAIN                  '
! i = i + 1
! write(f(i),2) 'PRECIP_TOT', 0, 'rcnst', 1, 0,0,0, -1, -1, -1, 0, 'kg m-2         ', 'TOTAL ACCUMULATED PRECIPITATION         '

! DTD: Milbrandt-Yau scheme constants
 i = i + 1
 write(f(i),2) 'NTC_MY    ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'meters-3       ', 'FIXED NT OF CLOUD FOR MY SCHEME         '
 i = i + 1
 write(f(i),2) 'RHO_QR_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF RAINWATER FOR MY SCHEME      '
 i = i + 1
 write(f(i),2) 'RHO_QI_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF ICE CRYSTALS FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'RHO_QS_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF SNOW FOR MY SCHEME           '
 i = i + 1
 write(f(i),2) 'RHO_QG_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF GRAUPEL FOR MY SCHEME        '
 i = i + 1
 write(f(i),2) 'RHO_QH_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF HAIL FOR MY SCHEME           '
 i = i + 1
 write(f(i),2) 'CNOR_MY   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'm4             ', 'INTERCEPT PARAM OF RAIN FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'CNOS_MY   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'm4             ', 'INTERCEPT PARAM OF SNOW FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'CNOG_MY   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'm4             ', 'INTERCEPT PARAM OF GRAUPEL FOR MY SCHEME'
 i = i + 1
 write(f(i),2) 'CNOH_MY   ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'm4             ', 'INTERCEPT PARAM OF HAIL FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'ALPHAR_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF RAIN FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'ALPHAI_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF ICE FOR MY SCHEME    '
 i = i + 1
 write(f(i),2) 'ALPHAS_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF SNOW FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'ALPHAG_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF GRAUPEL FOR MY SCHEME'
 i = i + 1
 write(f(i),2) 'ALPHAH_MY ', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'SHAPE PARAMETER OF HAIL FOR MY SCHEME   '
 i = i + 1
 write(f(i),2) 'RHO_QG_TAK', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF Takahashi graupel            '
 i = i + 1
 write(f(i),2) 'RHO_QH_TAK', 0, 'rcnst', 0, 0,0,0, -1, -1, -1, 0, 'kg meters-3    ', 'DENSITY OF Takahashi hail/frozen drops  '

! X1D
 i = i + 1
 write(f(i),2) 'XC        ', 1, 'x1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN X               '
 i = i + 1
 write(f(i),2) 'XE        ', 1, 'x1d  ', 0, 1,0,0, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN X            '
 i = i + 1
 write(f(i),2) 'DXC       ', 1, 'x1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR X-GRID CENTERS     '
 i = i + 1
 write(f(i),2) 'DXE       ', 1, 'x1d  ', 0, 1,0,0, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR X-GRID EDGE        '
#ifdef OPENDX
 i = i + 1
 write(f(i),2) 'XCDX      ', 2, 'x2d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN X for DX        '
 i = i + 1
 write(f(i),2) 'XEDX      ', 2, 'x2d  ', 0, 1,0,0, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN X for DX     '
#endif
! Y1D
 i = i + 1
 write(f(i),2) 'YC        ', 1, 'y1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN Y               '
 i = i + 1
 write(f(i),2) 'YE        ', 1, 'y1d  ', 0, 0,1,0, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN Y            '
 i = i + 1
 write(f(i),2) 'DYC       ', 1, 'y1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR Y-GRID CENTERS     '
 i = i + 1
 write(f(i),2) 'DYE       ', 1, 'y1d  ', 0, 0,1,0, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR Y-GRID EDGE        '
#ifdef OPENDX
 i = i + 1
 write(f(i),2) 'YCDX      ', 2, 'y2d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN Y for DX        '
 i = i + 1
 write(f(i),2) 'YEDX      ', 2, 'y2d  ', 0, 0,1,0, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN Y for DX     '
#endif
! Z1D
 i = i + 1
 write(f(i),2) 'ZC        ', 1, 'z1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN Z               '
 i = i + 1
 write(f(i),2) 'ZE        ', 1, 'z1d  ', 0, 0,0,1, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN Z            '
 i = i + 1
 write(f(i),2) 'DZC       ', 1, 'z1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR Z-GRID CENTERS     '
 i = i + 1
 write(f(i),2) 'DZE       ', 1, 'z1d  ', 0, 0,0,1, -1, -1, -1, 0, 'meters-1       ', '1 / GRID SPACING FOR Z-GRID EDGE        '
 i = i + 1
 write(f(i),2) 'UINIT     ', 1, 'z1d  ', 1, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'BASE STATE U                            '
 i = i + 1
 write(f(i),2) 'VINIT     ', 1, 'z1d  ', 1, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'BASE STATE V                            '
 i = i + 1
 write(f(i),2) 'WINIT     ', 1, 'z1d  ', 0, 0,0,1, -1, -1, -1, 0, 'meters sec-1   ', 'BASE STATE W                            '
 i = i + 1
 write(f(i),2) 'PIINIT    ', 1, 'z1d  ', 1, 0,0,0, -1, -1, -1, 0, 'count          ', 'BASE STATE PI                           '
 i = i + 1
 write(f(i),2) 'KMINIT    ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'meters2 sec-1  ', 'BASE STATE KM                           '
 i = i + 1
 write(f(i),2) 'KMBASERM  ', 1, 'z1d  ', 1, 0,0,0,  1, -1, -1, 0, 'meters2 sec-1  ', 'BASE STATE KM SUBTRACTION               '
 i = i + 1
 write(f(i),2) 'THINIT    ', 1, 'z1d  ', 1, 0,0,0,  1, -1, -1, 0, 'kelvin         ', 'BASE STATE THETA                        '
 i = i + 1
 write(f(i),2) 'QVINIT    ', 1, 'z1d  ', 1, 0,0,0,  1, -1, -1, 0, 'kg kg-1        ', 'BASE STATE QV                           '
 i = i + 1
 write(f(i),2) 'CCNUFINIT ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'meters-3       ', 'BASE STATE CCCNUF                       '
 i = i + 1
 write(f(i),2) 'CCNINIT   ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'meters-3       ', 'BASE STATE CCCN                         '
 i = i + 1
 write(f(i),2) 'CPIONINIT ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'BASE STATE POSITIVE ION CONCENTRATION   '
 i = i + 1
 write(f(i),2) 'CNIONINIT ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'BASE STATE NEGATIVE ION CONCENTRATION   '
 i = i + 1
 write(f(i),2) 'MUPOSZC   ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'POSITIVE ION MOBILITY AT SCALAR POINTS  '
 i = i + 1
 write(f(i),2) 'MUNEGZC   ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'NEGATIVE ION MOBILITY AT SCALAR POINTS  '
 i = i + 1
 write(f(i),2) 'MUPOSZE   ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'POSITIVE ION MOBILITY AT EDGE POINTS    '
 i = i + 1
 write(f(i),2) 'MUNEGZE   ', 1, 'z1d  ', 0, 0,0,1,  1, -1, -1, 0, 'meters-3       ', 'NEGATIVE ION MOBILITY AT EDGE POINTS    '
 i = i + 1
 write(f(i),2) 'UINIT0    ', 1, 'z1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'BASE STATE U AT TIME 0                  '
 i = i + 1
 write(f(i),2) 'VINIT0    ', 1, 'z1d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters sec-1   ', 'BASE STATE V AT TIME 0                  '
 i = i + 1
 write(f(i),2) 'PIINIT0   ', 1, 'z1d  ', 0, 0,0,0, -1, -1, -1, 0, 'count          ', 'BASE STATE PI AT TIME 0                 '
 i = i + 1
 write(f(i),2) 'THINIT0   ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'kelvin         ', 'BASE STATE THETA AT TIME 0              '
 i = i + 1
 write(f(i),2) 'QVINIT0   ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'kg kg-1        ', 'BASE STATE QV AT TIME 0                 '

! IF ( (.not. elec) .and. (.not. idoniconly )) THEN
IF ( (.not. elec) ) THEN
 i = i + 1
 write(f(i),2) 'EX        ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD X-COMPONENT              '
 i = i + 1
 write(f(i),2) 'EY        ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD Y-COMPONENT              '
 i = i + 1
 write(f(i),2) 'EZ        ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD Z-COMPONENT              '
 i = i + 1
 write(f(i),2) 'EMAG      ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD MAGNITUDE                '
 i = i + 1
 write(f(i),2) 'POTENTIAL ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt           ', 'ELECTRIC POTENTIAL                      '
 i = i + 1
 write(f(i),2) 'FLSHN     ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'NEGATIVE CHANNEL SEGMENTS               '
 i = i + 1
 write(f(i),2) 'FLSHP     ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'POSITIVE CHANNEL SEGMENTS               '
 i = i + 1
 write(f(i),2) 'FLSHI     ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH INITIATION POINTS                 '
 i = i + 1
 write(f(i),2) 'RSCGHIS   ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'C m-3 s-1      ', 'NONINDUCTIVE GRAUPEL CHARGING RATE      '
 i = i + 1
 write(f(i),2) 'RSCGHW    ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'C m-3 s-1      ', 'INDUCTIVE GRAUPEL CHARGING RATE         '
 i = i + 1
 write(f(i),2) 'SCNET     ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'C m-3          ', 'NET CHARGE DENSITY                      '
! i = i + 1
! write(f(i),2) 'FLSHR     ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH RATE                              '
! i = i + 1
! write(f(i),2) 'FLSHSRC    ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH RATE                              '
! i = i + 1
! write(f(i),2) 'FLSHR8KM  ', 1, 'z1d  ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH RATE                              '
! i = i + 1
! write(f(i),2) 'D-EX      ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD X-COMPONENT CHANGE       '
! i = i + 1
! write(f(i),2) 'D-EY      ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD Y-COMPONENT CHANGE       '
! i = i + 1
! write(f(i),2) 'D-EZ      ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD Z-COMPONENT CHANGE       '
! i = i + 1
! write(f(i),2) 'D-EMAG    ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt meters-1  ', 'ELECTRIC FIELD MAGNITUDE CHANGE         '
! i = i + 1
! write(f(i),2) 'DPOTENTIAL', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'volt           ', 'ELECTRIC POTENTIAL CHANGE               '
ENDIF

IF ( ikf == 0 .and. ioutput_vzf == 0) THEN
! i = i + 1
! write(f(i),2) 'SSFORCE   ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'count          ', 'SATURATION FORCING                      '
 i = i + 1
 write(f(i),2) 'VZF       ', 1, 'z1d  ', 0, 0,0,0,  0, -1, -1, 0, 'count          ', 'POWER-WEIGHTED FALL VELOCITY            '
ENDIF

#ifdef OPENDX
 i = i + 1
 write(f(i),2) 'ZCDX      ', 2, 'z2d  ', 0, 0,0,0, -1, -1, -1, 0, 'meters         ', 'SCALAR GRID POSITION IN Z for DX        '
 i = i + 1
 write(f(i),2) 'ZEDX      ', 2, 'z2d  ', 0, 0,0,1, -1, -1, -1, 0, 'meters         ', 'STAGGERED GRID POSITION IN Z for DX     '
#endif

 i = i + 1
 write(f(i),2) 'RANARRAY2D', 2, 'xy2d ', 0, 0,0,0, -1, -1, -1, 0, 'count          ', '2D Random Pert array                    '

! XY2D
IF ( ikf > 0 ) THEN
! i = i + 1
! write(f(i),2) 'SSFORCE   ', 2, 'xy2d ', 0, 0,0,0,  1, -1, -1, 0, 'count          ', 'SATURATION FORCING                      '
ENDIF

 i = i + 1
 write(f(i),2) 'RAIN_RAT  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2 s-1     ', 'SFC PRECIP RATE FOR LIQUID              '
 i = i + 1
 write(f(i),2) 'HAIL_RAT  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2 s-1     ', 'SFC PRECIP RATE FOR ICE                 '
 i = i + 1
 write(f(i),2) 'RAIN_ACC  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP               '
 i = i + 1
 write(f(i),2) 'HAIL_ACC  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED FROZEN PRECIP               '

 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACC2 ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP (grid motion) '

 IF ( microphys(1:4) == 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or. microphys(1:3) == 'TAK' ) THEN
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_ACCH ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED HAIL                        '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_ACCH2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED HAIL (grid motion)          '
 IF ( frozendrops >= 1 ) THEN
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_ACCHF', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED HAIL FROM FD                '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_NACCF', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED HAIL NUMBER from FD         '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_ACHF2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED HAIL FROM FD (gm)           '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_NACF2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED HAIL NUMBER from FD (gm)    '
 ENDIF
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_NACC ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED HAIL NUMBER                 '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_NACC2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED HAIL NUMBER (grid motion)   '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'HAIL_DIAM ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'mm             ', 'MEAN MASS HAIL Diameter (mm)            '
ENDIF

 IF ( microphys(1:3) == 'ZVD' .or. microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) == 'TAK' ) THEN
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'GR_DIAM   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'mm             ', 'MEAN MASS Graupel Diameter (mm)         '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'GR_NACC   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED Graupel NUMBER              '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'GR_NACC2  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED Graupel NUMBER (grid motion)'
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'GR_ACC    ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED Graupel Mass                '
  nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'GR_ACC2   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED Graupel Mass (grid motion)  '
 
IF ( frozendrops >= 1 ) THEN
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_DIAM   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'mm             ', 'MEAN MASS FD Diameter (mm)              '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_NACC   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED FD NUMBER                   '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_ACC    ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED FROZEN DROPS                '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_NACC2  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2            ', 'ACCUMULATED FD NUMBER (gm)              '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_ACC2   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED FROZEN DROPS (gm)           '
ENDIF 
ENDIF

  IF ( iraintypes > 0 .and. ( microphys(1:3) .eq. 'ZVD' .or. microphys(1:4) .eq. 'ZIEG' ) ) THEN

 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACC_A', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM WARM RAIN'
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACC_M', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM MELTING  '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACC_S', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM SHEDDING '

 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACCA2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM WARM RAIN'
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACCM2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM MELTING  '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'RAIN_ACCS2', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP FROM SHEDDING '

 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_AUTO   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from warm rain             '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_SHED   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from shedding              '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_MELT   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from melting               '

 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_AUTO2  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from warm rain (gm)        '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_SHED2  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from shedding (gm)         '
 nprecip = nprecip + 1
 i = i + 1
 write(f(i),2) 'FD_MELT2  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'Frozen Drops from melting (gm)          '

  ENDIF
! nprecip = 4 : Note value already set to 4 in index_module
! turn this on later when there is a variable for the size of the precip array
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'COMP_DBZ  ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'dBZ            ', 'COMPOSITE REFLECTIVITY                  '
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'WZ_SFC_MAX', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2 s-2        ', 'MAX SURFACE VORTICITY OVER HISTORY TIME '
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'WZ_SFC_MIN', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm-2 s-2        ', 'MIN SURFACE VORTICITY OVER HISTORY TIME '
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'UV_SFC_MAX', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm s-1          ', 'MAX SURFACE WIND SPEED OVER HISTORY TIME'

 IF ( microphys(1:3) == 'ZVD' .or. microphys(1:4) .eq. 'ZIEG'  ) THEN
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'HAIL_MAX2D', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm              ', 'MAX HAIL DIAMETER IN COLUMN             '

 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'HAIL_MAXK1', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'm              ', 'MAX HAIL DIAMETER AT SURFACE (k=1)      '
 
 ENDIF
 
 IF ( pwf_flag ) THEN
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'PWF_RAT   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'SFC LIQUID FRACTION RATE                '
 i = i + 1
 nprecip = nprecip + 1
 write(f(i),2) 'PWF_ACC   ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID FRACTION             '
 ENDIF

! i = i + 1
! nprecip = nprecip + 1
! write(f(i),2) 'RAIN_ACC2 ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED LIQUID PRECIP GRND REL      '
! i = i + 1
! nprecip = nprecip + 1
! write(f(i),2) 'HAIL_ACC2 ', 2, 'xy2d ', 1, 0,0,0,  1, -1, -1, 0, 'kg m-2         ', 'ACCUMULATED FROZEN PRECIP GRND REL      '

 neelec2d = 0
  it = 0
! IF ( elec .or. ikf > 0 .or. idoniconly ) THEN 
 IF ( elec .or. ikf > 0 ) THEN 
 it = 1
 neelec2d = 0
 
 IF ( elec ) THEN
! neelec2d = neelec2d + 1
! i = i + 1
! write(f(i),2) 'FLSHR     ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH RATE IN A COLUMN                  '

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHSRC   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH SOURCES IN A COLUMN               '

 IF ( ioutput_flshr8km > 0 ) THEN
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSH8KM   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH RATE IN A COLUMN 8-KM FED         '
 ENDIF
 
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFED   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH EXTENT DENSITY (Column)            '

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFOD   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'FLASH ORIGIN RATE IN A COLUMN            '

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFODIC ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'IC FLASH ORIGIN RATE IN A COLUMN         '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFODCGN', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'CGN FLASH ORIGIN RATE IN A COLUMN        '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFODCGP', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'CGP FLASH ORIGIN RATE IN A COLUMN        '

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFEDIC ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'IC FLASH FED                             '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFEDICP', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'IC FLASH FED POS CHAN                    '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFEDICN', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'IC FLASH FED NEG CHAN                    '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFEDCGN', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'CGN FLASH FED                            '
 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'FLSHFEDCGP', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'count          ', 'CGP FLASH FED                            '

 ENDIF

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'PNIC_2D   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'C m-2          ', 'POS NONINDUCTIVE CHARGING IN A COLUMN   '

 neelec2d = neelec2d + 1
 i = i + 1
 write(f(i),2) 'NNIC_2D   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'C m-2          ', 'NEG NONINDUCTIVE CHARGING IN A COLUMN   '


! neelec2d = neelec2d + 1
! i = i + 1
! write(f(i),2) 'TNIC_2D   ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'C m-2          ', 'TOTAL NONINDUCTIVE CHARGING IN A COLUMN  '

 ENDIF ! elec
 
!.....msb add sfc phyics variables
  it = 0 ; 
  
  IF ( isfc >=1 ) THEN 
    it = 1
 i = i + 1
 write(f(i),2) 'TCANP     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'K              ', 'CANOPY POTENTIAL TEMPERATURE            '
 i = i + 1
 write(f(i),2) 'WCANP     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'kg kg-1        ', 'CANOPY MIXING RATIO                     '
 i = i + 1
 write(f(i),2) 'QAV       ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'kg kg-1        ', 'MIXING RATIO WITHIN FOLIAGE             '
 i = i + 1
 write(f(i),2) 'VEG       ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'fraction       ', 'VEGETATION DENSITY (FRACTION)           '
 i = i + 1
 write(f(i),2) 'STYPE     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'number         ', 'SOIL TYPE                               '
 i = i + 1
 write(f(i),2) 'TSRFC     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'K              ', 'SURFACE POTENTIAL TEMPERATURE           '
 i = i + 1
 write(f(i),2) 'WSFC      ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'kg kg-1        ', 'SURFACE MIXING RATIO                    '
 i = i + 1
 write(f(i),2) 'TSOIL     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'K              ', 'SOIL POTENTIAL TEMPERATURE              '
 i = i + 1
 write(f(i),2) 'WSOIL     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'kg m-3         ', 'SOIL MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'VLAI      ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'number         ', 'LEAF AREA INDEX                         '
 i = i + 1
 write(f(i),2) 'ALBEDO    ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'fraction       ', 'ALBEDO                                  '
 i = i + 1
 write(f(i),2) 'ROUGH     ', 2, 'xy2d ',it, 0,0,0,  1,  1, -1, 0, 'm              ', 'ROUGHNESS LENGTH                        '
 i = i + 1
 write(f(i),2) 'EFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'm-2 s-1        ', 'VEGETATION MOISTURE FLUX                '
 i = i + 1
 write(f(i),2) 'FFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'm-2 s-1        ', 'EVAPORATION MOISTURE FLUX               '
 i = i + 1
 write(f(i),2) 'UFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'm-1 s-2        ', 'U MOMENTUM FLUX                         '
 i = i + 1
 write(f(i),2) 'VFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'm-1 s-2        ', 'V MOMENTUM FLUX                         '
 i = i + 1
 write(f(i),2) 'TFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'W m-2          ', 'SFC SENSIBLE HEAT FLUX                  '
 i = i + 1
 write(f(i),2) 'QFLX      ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'W m-2          ', 'SFC MOISTURE FLUX                       '
 i = i + 1
 write(f(i),2) 'RADSW     ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'W m-2          ', 'NET SW RADIATION                        '
 i = i + 1
 write(f(i),2) 'RADLW     ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'W m-2          ', 'NET LW RADIATION                        '

!......sfc physics

 i = i + 1
 write(f(i),2)  'KH        ', 2, 'xy2d ',it, 0,0,0,  1, -1, -1, 0, 'm2 s-1         ', 'HEAT EDDY MIXING COEFFICIENT            '

  ENDIF

! XZ2D
!mpi!  i = i + 1
! write(f(i),2) 'VS        ', 2, 'xz2d ', 1, 0,1,0,  0, -1, -1, 0, 'meters sec-1   ', 'NORMAL VELOCITY ON SOUTH BOUNDARY       '
!mpi!  i = i + 1
! write(f(i),2) 'VN        ', 2, 'xz2d ', 1, 0,1,0,  0, -1, -1, 0, 'meters sec-1   ', 'NORMAL VELOCITY ON NORTH BOUNDARY       '

! YZ2D
!mpi!  i = i + 1
! write(f(i),2) 'UW        ', 2, 'yz2d ', 1, 1,0,0,  0, -1, -1, 0, 'meters sec-1   ', 'NORMAL VELOCITY ON WEST BOUNDARY        '
!mpi!  i = i + 1
! write(f(i),2) 'UE        ', 2, 'yz2d ', 1, 1,0,0,  0, -1, -1, 0, 'meters sec-1   ', 'NORMAL VELOCITY ON EAST BOUNDARY        '

! XYZ3D
IF ( ikf > 0  .or. ioutput_vzf > 0 ) THEN
 i = i + 1
 write(f(i),2) 'VZF       ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'count          ', 'POWER-WEIGHTED FALL VELOCITY            '
ENDIF
 i = i + 1
 write(f(i),2) 'DBZ       ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'RADAR REFLECTIVITY                      '
 i = i + 1
 write(f(i),2) 'WZ        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'sec-1          ', 'VERTICAL VORTICITY                      '

! KEEP WZ as the last "regular" diagnostic variable.  Put "Extra" diag. variables between WZ and EX.

 IF ( tke_type == 2 ) THEN
  i = i + 1
  write(f(i),2) 'KMT       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters2 sec-1  ', 'SCAL MIXING COEFFICIENT                 '
  i = i + 1
  write(f(i),2) 'KHT       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters2 sec-1  ', 'MOM MIXING COEFFICIENT                  '
 ENDIF

 IF ( ioutput_tkediss > 0 .and. iturbenhance > 0 ) THEN
  i = i + 1
  write(f(i),2) 'TKE_DISS  ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters2 sec-3  ', 'TKE Dissipation Rate                    '
 ENDIF
 
 IF ( ioutput_vzf > 0 ) THEN
  i = i + 1
  write(f(i),2) 'WVZF      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'count          ', 'W-POWER-WEIGHTED FALL VELOCITY          '
 ENDIF

 IF ( ioutput_pres3d > 0 ) THEN
  i = i + 1
  write(f(i),2) 'P         ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'Pa             ', 'AIR PRESSURE                            '
 ENDIF
 IF ( ioutput_temC3d > 0 ) THEN
  i = i + 1
  write(f(i),2) 'T         ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kelvin         ', 'AIR TEMPERATURE                         '
 ENDIF
 IF ( ioutput_thv3d > 0 ) THEN
  i = i + 1
  write(f(i),2) 'THETA_V   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kelvin         ', 'Virtual Pot. Temp.                      '
 ENDIF

  IF ( ioutput_hwind3d > 0 ) THEN
  i = i + 1
  write(f(i),2) 'HWIND     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'm s-1          ', 'Horizontal Wind Speed                   '
  ENDIF

  IF ( ioutput_ssw > 0 ) THEN
  i = i + 1
  write(f(i),2) 'SSAT      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT LIQUID              '
  ENDIF

  IF ( ioutput_ssi > 0 ) THEN
  i = i + 1
  write(f(i),2) 'SSATI     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT ICE                 '
  ENDIF


  IF ( ioutput_cwdia > 0 ) THEN
  i = i + 1
  write(f(i),2) 'CWDIA     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'microns        ', 'Droplet Mean Diameter                   '
  ENDIF

  IF ( ioutput_cwdia > 0 ) THEN
  i = i + 1
  write(f(i),2) 'CWNU      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'microns        ', 'Droplet shape parameter                 '
  ENDIF

  IF ( ioutput_snowstuff > 0 ) THEN
  i = i + 1
  write(f(i),2) 'SWDIA     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Snow Mean Diameter                      '

  i = i + 1
  write(f(i),2) 'SWMASS    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mg             ', 'Snow Mean Mass                          '

    i = i + 1
  write(f(i),2) 'SWAGG     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 's-1            ', 'Snow Aggregation rate                   '

   IF ( ioutput_snowstuff > 1 ) THEN
    i = i + 1
  write(f(i),2) 'SWDN      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg m-3         ', 'Snow Density                            '
   
   ENDIF

  ENDIF

 IF ( ioutput_icedensity > 0 ) THEN
   i = i + 1
  write(f(i),2) 'HWDN      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg m-3         ', 'Graupel Density                     '

  IF ( hl_flag ) THEN
  i = i + 1
  write(f(i),2) 'HLDN      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg m-3         ', 'Hail Density                        '
  ENDIF

  IF ( frozendrops == 1 ) THEN
    i = i + 1
  write(f(i),2) 'FWDN      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg m-3         ', 'Frozen Drop Density                     '
  ENDIF
 ENDIF

  IF ( ioutput_rh > 0 ) THEN
  i = i + 1
  write(f(i),2) 'RH        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'Relative Humidity                       '
  ENDIF


  IF ( microphys(llen:llen) == 'X' .or.  microphys(llen-1:llen-1) == 'X' ) THEN

!  i = i + 1
!  write(f(i),2) 'HWIND     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'm s-1          ', 'Horizontal Wind Speed                   '

!  i = i + 1
!  write(f(i),2) 'SSAT      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT WATER               '
!  i = i + 1
!  write(f(i),2) 'SSATI     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT ICE                 '
!  i = i + 1
!  write(f(i),2) 'PPERT     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'Pa             ', 'PERTURBATION PRESSURE                   '

!  ENDIF

! i = i + 1
! write(f(i),2) 'CRCNW     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'NUMBER AUTOCONVERSION RATE              '
!
! i = i + 1
! write(f(i),2) 'CRFRZ     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'BIGG FREEZING RATE                      '
  
 ENDIF

  IF ( microphys(1:3) == 'TAK' ) THEN
  i = i + 1
  write(f(i),2) 'CCWTAK    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'cm-3           ', 'DROPLET NUMBER                          '
  i = i + 1
  write(f(i),2) 'CCITAK    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'cm-3           ', 'Ice Crystal NUMBER                      '
  i = i + 1
  write(f(i),2) 'CRWTAK    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'cm-3           ', 'Rain NUMBER                          '
  i = i + 1
  write(f(i),2) 'CHWTAK    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'cm-3           ', 'Graupel NUMBER                          '
  i = i + 1
  write(f(i),2) 'CHLTAK    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'cm-3           ', 'Hail NUMBER                          '
  
  ENDIF
! special fields for squall line workshop output
  IF ( ( microphys(1:1) == 'Z' .or. microphys(1:3) == 'TAK' ) .and. (microphys(llen:llen) == 'W' .or. ioutput_workshop > 0 ) ) THEN

  i = i + 1
  write(f(i),2) 'VR        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'RAIN FALL SPEED                         '
  i = i + 1
  write(f(i),2) 'VG        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'GRAUPEL FALL SPEED                      '
  IF ( frozendrops == 1 ) THEN
  i = i + 1
  write(f(i),2) 'VF        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'Frozen drop FALL SPEED                  '
  ENDIF
  IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or. microphys(1:5) .eq. 'ZVDMH' &
       .or. microphys(1:4) .eq. 'ZMRH' .or. ( microphys(1:3) == 'TAK' .and. takkfmax >= 2 ) ) THEN
  i = i + 1
  write(f(i),2) 'VH        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'HAIL FALL SPEED                         '
  ENDIF
  i = i + 1
  write(f(i),2) 'VI        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'Cloud ICE FALL SPEED                    '
  i = i + 1
  write(f(i),2) 'VS        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'SNOW FALL SPEED                         '
  i = i + 1
  write(f(i),2) 'D0        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'RAIN MEDIAN DIAMETER                    '

  i = i + 1
  write(f(i),2) 'QAUT      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'CLOUD-RAIN AUTOCONVERION RATE           '
  i = i + 1
  write(f(i),2) 'QACC      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QRACW RAIN AC OF DROPLETS               '
  i = i + 1
  write(f(i),2) 'QREVAP    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QRCEV RAIN EVAPORATION RATE             '
  i = i + 1
  write(f(i),2) 'QCEVAP    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QCEV CLOUD EVAPORATION RATE             '
  i = i + 1
  write(f(i),2) 'QCOND     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QCOND CLOUD CONDENSATION RATE           '
  i = i + 1
  write(f(i),2) 'QDEP      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PDEP TOTAL DEPOSITION RATE              '
  i = i + 1
  write(f(i),2) 'QMELT     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PMLT TOTAL MELTING RATE                 '
  i = i + 1
  write(f(i),2) 'QSUB      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PSUB TOTAL SUBLIMATION RATE             '

  IF ( ioutput_ssw <= 0 ) THEN ! do not duplicate if already added
  i = i + 1
  write(f(i),2) 'SSAT      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT WATER               '
  ENDIF
  IF ( ioutput_ssi <= 0 ) THEN ! do not duplicate if already added
  i = i + 1
  write(f(i),2) 'SSATI     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'SUPERSATURATION WRT ICE                 '
  ENDIF
    
  ENDIF

  IF ( ( microphys(1:1) == 'Z' .or. microphys(1:3) == 'TAK' ) .and. ( ioutput_dbzsedchange > 0 ) ) THEN
  i = i + 1
  write(f(i),2) 'DBZCHANGEH', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'Change in graupel dBZ from sedimentation'
  i = i + 1
  write(f(i),2) 'DBZCHANGER', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'Change in rain dBZ from sedimentation   '

  ENDIF
  
  IF ( ( microphys(1:1) == 'Z' .or. microphys(1:3) == 'TAK' ) .and. ( ioutput_sedmeltstuff > 0 ) ) THEN
  
  IF ( hl_flag ) THEN
  IF ( ioutput_workshop == 0 ) THEN
  i = i + 1
  write(f(i),2) 'VH        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'HAIL Mass-wgt FALL SPEED                '
  ENDIF
  
!   i = i + 1
!   write(f(i),2) 'VNH       ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'HAIL NUM-WGT FALL SPEED                 '
! 
!   i = i + 1
!   write(f(i),2) 'VZH       ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'HAIL Z-WGT FALL SPEED                   '
! 
!   i = i + 1
!   write(f(i),2) 'VMNHL     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'ratio          ', 'HAIL M/N FALL SPEED Ratio               '

  i = i + 1
  write(f(i),2) 'ALPHAHL3D ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'ratio          ', 'HAIL Shape Parameter                    '

  ENDIF

  IF ( frozendrops == 1 ) THEN
  i = i + 1
  write(f(i),2) 'FWDIA     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Frozen drops Mean Mass Diameter         '
  ENDIF

  IF ( hl_flag ) THEN
! if hail is on, because otherwise hldia becomes hwdia
  i = i + 1
  write(f(i),2) 'HWDIA     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Graupel Mean Mass Diameter              '
  ENDIF
  
  i = i + 1
  write(f(i),2) 'HLDIA     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Hail Mean Mass Diameter                 '
  i = i + 1
  write(f(i),2) 'HLMWDIA   ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Hail Mass-wgt Diameter                  '
  i = i + 1
  write(f(i),2) 'HLCHARDIA ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Hail Characteristic Diameter            '

  IF ( hl_flag ) THEN
  i = i + 1
  write(f(i),2) 'CHLCNH    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-3 s-1        ', 'Graupel->Hail # conversion rate         '

  i = i + 1
  write(f(i),2) 'DHLCNH    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'Graupel->Hail conversion diam           '
  ENDIF
  
  i = i + 1
  write(f(i),2) 'ALPHAR3D  ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'ratio          ', 'Rain Shape Parameter                    '

  i = i + 1
  write(f(i),2) 'ALPHAH3D  ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'ratio          ', 'Graupel Shape Parameter                 '


  IF ( frozendrops >= 1 ) THEN
  i = i + 1
  write(f(i),2) 'ALPHAF3D  ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'ratio          ', 'Frozen drops Shape Parameter            '

  IF ( hl_flag ) THEN
  i = i + 1
  write(f(i),2) 'CHLCNF    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-3 s-1        ', 'FrozenDrops->Hail # conversion rate     '
  i = i + 1
  write(f(i),2) 'DHLCNF    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'FrozenDrops->Hail conversion diam       '
  ENDIF ! hl_flag
  ENDIF
  
  ENDIF

! special fields for squall line workshop output
  IF ( ( microphys(1:4) == 'MORR' .or. microphys(1:4) == 'THOM' ) .and.  ioutput_workshop > 0 ) THEN

!  i = i + 1
!  write(f(i),2) 'VR        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'RAIN FALL SPEED                         '
!  i = i + 1
!  write(f(i),2) 'VG        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'GRAUPEL FALL SPEED                      '
!  i = i + 1
!  write(f(i),2) 'VI        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'Cloud ICE FALL SPEED                    '
!  i = i + 1
!  write(f(i),2) 'VS        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'm-1 s-1        ', 'SNOW FALL SPEED                         '
  i = i + 1
  write(f(i),2) 'D0        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'mm             ', 'RAIN MEDIAN DIAMETER                    '
!
!  i = i + 1
!  write(f(i),2) 'QAUT      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'CLOUD-RAIN AUTOCONVERION RATE           '
!  i = i + 1
!  write(f(i),2) 'QACC      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QRACW RAIN AC OF DROPLETS               '
!  i = i + 1
!  write(f(i),2) 'QREVAP    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QRCEV RAIN EVAPORATION RATE             '
!  i = i + 1
!  write(f(i),2) 'QCEVAP    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QCEV CLOUD EVAPORATION RATE             '
!  i = i + 1
!  write(f(i),2) 'QCOND     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'QCOND CLOUD CONDENSATION RATE           '
!  i = i + 1
!  write(f(i),2) 'QDEP      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PDEP TOTAL DEPOSITION RATE              '
!  i = i + 1
!  write(f(i),2) 'QMELT     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PMLT TOTAL MELTING RATE                 '
!  i = i + 1
!  write(f(i),2) 'QSUB      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg kg-1 s-1    ', 'PSUB TOTAL SUBLIMATION RATE             '

  ENDIF


! special fields for Takahashi bin
  IF ( (  microphys(1:3) == 'TAK' )  ) THEN
 i = i + 1
 write(f(i),2) 'DBZR      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'RAIN RADAR REFLECTIVITY                 '
 i = i + 1
 write(f(i),2) 'DBZH      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'GRAUPEL RADAR REFLECTIVITY              '
 IF ( takkfmax >= 2 ) THEN
 i = i + 1
 write(f(i),2) 'DBZHL     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'HAIL RADAR REFLECTIVITY                 '
 ENDIF
 i = i + 1
 write(f(i),2) 'DBZI      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'dBZ            ', 'ICE RADAR REFLECTIVITY                  '

  IF ( ioutput_takrates > 0 ) THEN ! do not duplicate if already added
  i = i + 1
  write(f(i),2) 'QH2HL     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg s-1  ', 'GRAUPEL TO HAIL CONVERSION               '
  i = i + 1
  write(f(i),2) 'QHL2H     ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg s-1  ', 'HAIL TO GRAUPEL CONVERSION               '
  ENDIF



  ENDIF
  
  IF ( ioutput_mltshedsizerates > 0 .and. (  microphys(1:3) == 'TAK' .or.  microphys(1:3) == 'ZVD' ) ) THEN
  i = i + 1
  write(f(i),2) 'D0SHED    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 's-1            ', 'Drop shedding rate  0-D1               '
  i = i + 1
  write(f(i),2) 'D1SHED    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 's-1            ', 'Drop shedding rate D1-D2               '
  i = i + 1
  write(f(i),2) 'D2SHED    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 's-1            ', 'Drop shedding rate D2-D3               '
  i = i + 1
  write(f(i),2) 'D3SHED    ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 's-1            ', 'Drop shedding rate D3-Inf              '
  ENDIF
  

! IF ( elec .or. idoniconly ) THEN
IF ( elec ) THEN

! unstaggered: (i.e., at scalar points)
! i = i + 1
! write(f(i),2) 'EX        ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD X-COMPONENT              '
! i = i + 1
! write(f(i),2) 'EY        ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD Y-COMPONENT              '
! i = i + 1
! write(f(i),2) 'EZ        ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD Z-COMPONENT              '

! staggered:
 i = i + 1
 write(f(i),2) 'EX        ', 3, 'xyz3d', 1, 1,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD X-COMPONENT              '
 i = i + 1
 write(f(i),2) 'EY        ', 3, 'xyz3d', 1, 0,1,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD Y-COMPONENT              '
 i = i + 1
 write(f(i),2) 'EZ        ', 3, 'xyz3d', 1, 0,0,1,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD Z-COMPONENT              '

 i = i + 1
 write(f(i),2) 'EMAG      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'ELECTRIC FIELD MAGNITUDE                '
 i = i + 1
 write(f(i),2) 'POTENTIAL ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt           ', 'ELECTRIC POTENTIAL                      '
 i = i + 1
 write(f(i),2) 'FLSHN     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'count          ', 'NEGATIVE CHANNEL SEGMENTS               '
 i = i + 1
 write(f(i),2) 'FLSHP     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'count          ', 'POSITIVE CHANNEL SEGMENTS               '
 i = i + 1
 write(f(i),2) 'FLSHI     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'count          ', 'FLASH INITIATION POINTS                 '
 i = i + 1
 write(f(i),2) 'RSCGHIS   ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'NONINDUCTIVE GRAUPEL CHARGING RATE      '
 i = i + 1
 write(f(i),2) 'RSCGHW    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'INDUCTIVE GRAUPEL CHARGING RATE         '
 i = i + 1
 write(f(i),2) 'SCNET     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3          ', 'NET CHARGE DENSITY                      '

 i = i + 1
 write(f(i),2) 'EMAGSTAG  ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'MAX ELECTRIC FIELD MAGNITUDE            '

! unstaggered:
! IF ( microphys(llen:llen) == '2' ) THEN
 IF ( .false. ) THEN
 i = i + 1
 write(f(i),2) 'EX2       ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'UNSTAG ELECTRIC FIELD X-COMPONENT       '
 i = i + 1
 write(f(i),2) 'EY2       ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'UNSTAG ELECTRIC FIELD Y-COMPONENT       '
 i = i + 1
 write(f(i),2) 'EZ2       ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'volt meters-1  ', 'UNSTAG ELECTRIC FIELD Z-COMPONENT       '
 ENDIF
 
 IF ( microphys(llen:llen) == 'X' .or. ioutput_xtrachgsep > 0 .or. ieopt > 0 ) THEN

 i = i + 1
 write(f(i),2) 'RSCADDL   ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'CHARGING RATE VIA ADDITIONAL MECHANISMS '
 i = i + 1
 write(f(i),2) 'RSCNOLIQ  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 s-1      ', 'LIQUID-FREE GRAUPEL CHARGING RATE       '

 ENDIF
 
 IF ( inetchargetend > 0 ) THEN

 i = i + 1
 write(f(i),2) 'SCTNDADV  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE ADVECTION TENDENCY           '

 i = i + 1
 write(f(i),2) 'SCTNDADVSN', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE ADVECTION TENDENCY SIGN      '
 
 i = i + 1
 write(f(i),2) 'SCTNDMIX  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE MIXING TENDENCY              '
 
 i = i + 1
 write(f(i),2) 'SCTNDNIC  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'GRAUPEL/HAIL CHARGING (INDUCT + NIC)    '
 
 i = i + 1
 write(f(i),2) 'SCTNDSED  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE SEDIMENTATION TENDENCY       '
 
 i = i + 1
 write(f(i),2) 'SCTNDION  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE ION DRIFT TENDENCY           '
 
 i = i + 1
 write(f(i),2) 'SCTNDLGT  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3 min-1    ', 'NET CHARGE LIGHTNING TENDENCY           '
 
 i = i + 1
 write(f(i),2) 'SCTNDAVG  ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'C m-3          ', 'AVERAGE NET CHARGE                      '
 
 ENDIF
 
ENDIF ! elec


 i = i + 1
 write(f(i),2) 'U         ', 3, 'xyz3d', 1, 1,0,0,  0,  1,  0, 0, 'meters sec-1   ', 'X-WIND COMPONENT                        '
 i = i + 1
 write(f(i),2) 'V         ', 3, 'xyz3d', 1, 0,1,0,  0,  1,  0, 0, 'meters sec-1   ', 'Y-WIND COMPONENT                        '
 i = i + 1
 write(f(i),2) 'W         ', 3, 'xyz3d', 1, 0,0,1,  0,  1,  0, 0, 'meters sec-1   ', 'Z-WIND COMPONENT                        '
 i = i + 1
 write(f(i),2) 'PI        ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'count          ', 'PERT. EXNER                             '

 i = i + 1
 write(f(i),2) 'RHO       ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'kg m-1         ', 'MOIST AIR DENSITY                       '

 IF ( tke_type == 1 ) THEN
 i = i + 1
 write(f(i),2) 'KM        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters2 sec-1  ', 'MIXING COEFFICIENT                      '
 ELSE
 i = i + 1
! write(f(i),2) 'KM        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters sec-1   ', 'TURBULENCE ENERGY (SQRT)                '
 write(f(i),2) 'TKE       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'meters sec-1   ', 'TURBULENCE ENERGY (SQRT)                '
 ENDIF
! Start scalars

 i = i + 1
 write(f(i),2) 'TH        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 1, 'kelvin         ', 'POTENTIAL TEMPERATURE                   '

! THIS IS IMPORTANT:  ALL VARIABLES FROM TH onward are treated as the scalar array - so dont put any variable below here you
! dont want treated as part of a 4D scalar array.

 IF ( microphys(1:3) .eq. 'DRY' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '

  IF ( microphys .eq. 'DRYC' ) THEN
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
  ENDIF

  IF ( microphys .eq. 'DRYT' ) THEN
 i = i + 1
 write(f(i),2) 'TRACER1   ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg kg-1        ', 'SCALAR TRACER 1                         '
  ENDIF

  IF ( microphys .eq. 'DRYR' ) THEN
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN CONCENTRATION                       '
! i = i + 1
! write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '
  ENDIF
  
 ENDIF
 
 
 IF ( microphys(1:7) .eq. 'KESSLER' .or.  microphys(1:7) .eq. 'WARMLFO' .or.  &
      microphys(1:8) .eq. 'WARMZIEG' .or. microphys(1:3) .eq. 'SCG' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '

 ENDIF
  
 IF ( microphys(1:3) .eq. 'LFO' .or. microphys(1:4) .eq. 'ZIEG' .or. microphys(1:1) .eq. 'Z'  &
      .or. microphys(1:4) .eq. 'THOM' .or. microphys(1:4) .eq. 'MORR' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '

 IF ( microphys(1:4) .eq. 'ZMRM' .or. microphys(1:5) .eq. 'ZMRHM') THEN
 i = i + 1
 write(f(i),2) 'QM        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'MELT/SHED RAIN MIXING RATIO             '
 ENDIF

 i = i + 1
 IF ( imixice == 1 ) THEN ! turb mixing for ice is on
 write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ICE MIXING RATIO                        '
 ELSE ! turb mixing for ice is off
 write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1, -1, 3, 'kg kg-1        ', 'ICE MIXING RATIO                        '
 ENDIF

 IF ( icespheres >= 1 ) THEN ! turn on ice spheres (frozen droplets)
 i = i + 1
 write(f(i),2) 'QIS       ', 3, 'xyz3d', 1, 0,0,0,  1,  1, -1, 3, 'kg kg-1        ', 'ICE SPHERE MIXING RATIO                 '
 ENDIF
 
 i = i + 1
 write(f(i),2) 'QS        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'SNOW MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QH        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GRAUPEL MIXING RATIO                    '
 IF ( frozendrops == 1 ) THEN ! turn on frozen drops
 i = i + 1
 write(f(i),2) 'QF        ', 3, 'xyz3d', 1, 0,0,0,  1,  1, -1, 3, 'kg kg-1        ', 'FROZEN DROP MIXING RATIO                '
 ENDIF


  IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or. microphys(1:5) .eq. 'ZVDMH' &
       .or. microphys(1:4) .eq. 'ZMRH' ) THEN
 i = i + 1
 write(f(i),2) 'QHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'HAIL MIXING RATIO                       '
  ENDIF


  IF ( microphys(1:4) .eq. 'ZVDM' .or. microphys(1:5) .eq. 'ZVDHM' .or. microphys(1:5) .eq. 'ZIEGM' &
        .or. microphys(1:6) .eq. 'ZIEGHM' .or. microphys(1:4) .eq. 'ZMRM' .or. microphys(1:5) .eq. 'ZMRHM') THEN
 i = i + 1
 write(f(i),2) 'QSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON SNOW MIXING RATIO       '
 i = i + 1
 write(f(i),2) 'QHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON GRAUPEL MIXING RATIO    '
 IF ( frozendrops == 1 ) THEN ! turn on frozen drops
 i = i + 1
 write(f(i),2) 'QFW       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON FRODROP MIXING RATIO    '
 ENDIF
  IF ( ipconc == 8 .and. ibinhmlr >= 1 ) THEN
 i = i + 1
 write(f(i),2) 'QHWLG     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON LARGE GRAUPEL           '
  ENDIF


  ENDIF
  
  IF ( microphys(1:5) .eq. 'ZVDHM' .or. microphys(1:5) .eq. 'ZVDMH' .or.  &
       microphys(1:6) .eq. 'ZIEGHM' .or. microphys(1:5) .eq. 'ZMRHM' ) THEN
 i = i + 1
 write(f(i),2) 'QHLW      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON HAIL MIXING RATIO       '

  IF ( ipconc == 8 .and. ibinhlmlr >= 1 ) THEN
! i = i + 1
! write(f(i),2) 'QHLWLG    ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'LIQUID WATER ON LARGE HAIL              '
  ENDIF
  ENDIF
  

 ENDIF


IF ( microphys(1:8) .eq. 'WARMZIEG'  ) THEN

 IF ( ipconc .ge. 2 ) THEN
 i = i + 1
 write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD CONCENTRATION                     '
 ENDIF
 IF ( ipconc .ge. 3 ) THEN
 i = i + 1
 write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '
 ENDIF
 
 IF ( ipconc .ge. 2 .and. ioutput_ssmx > 0 ) THEN
 i = i + 1
 write(f(i),2) 'SSMX      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'MAX SUPERSATURATION                     '
 ENDIF
! i = i + 1
! write(f(i),2) 'VHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'coulomb m-3    ', 'HAIL PARTICLE VOLUME                    '

  IF ( ipconc .ge. 8 ) THEN
 i = i + 1
 write(f(i),2) 'ZRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z meters-3     ', 'RAIN REFLECTIVITY                       '
  ENDIF

 ENDIF

IF ( microphys(1:4) .eq. 'THOM' ) THEN

 i = i + 1
 write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '

 i = i + 1
 write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'ICE CRYSTAL CONCENTRATION               '

ENDIF

IF ( microphys(1:4) .eq. 'MORR' ) THEN

 i = i + 1
 write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '

 i = i + 1
 write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'ICE CRYSTAL CONCENTRATION               '

 i = i + 1
 write(f(i),2) 'CSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'SNOW CONCENTRATION                      '

 i = i + 1
 write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GRAUPEL CONCENTRATION                   '

ENDIF

IF ( microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) .eq. 'ZVD' .or.  microphys(1:3) .eq. 'ZMR' ) THEN

 IF ( ipconc .ge. 2 ) THEN
 IF ( inucopt == 1 .or. inucopt == 3 .or. irenuc >= 3 .or. ac_opt >= 11 ) THEN ! ( .or. takcwnucopt >= 2 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NUMBER OF ACTIVATED CCN                 '
 IF ( ac_opt == 22 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNACO   ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NUMBER OF ACTIVATED CO-CCN              '
 i = i + 1
 write(f(i),2) 'CCCNANU   ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NUMBER OF ACTIVATED NU-CCN              '
 ENDIF
 ENDIF

 IF ( inucopt == 2 .or. inucopt == 3 .or. icenucopt == 4 .or. in_freeze_rain_first /= 0 ) THEN
 i = i + 1
  IF ( imixice == 1 ) THEN ! turb mixing for ice is on
 write(f(i),2) 'CCINA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NUMBER OF ACTIVATED ICE NUCLEI           '
  ELSE  ! turb mixing for ice is off
 write(f(i),2) 'CCINA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0, -1, 0, 'meters-3       ', 'NUMBER OF ACTIVATED ICE NUCLEI           '
  ENDIF
 ENDIF

 IF ( inucopt == -1  ) THEN
 i = i + 1
 write(f(i),2) 'CCIN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD ICE NUCLEI CONCENTRATION          '
 ENDIF

 IF ( Abs(ccnuf) > 0.01 .or. iccnufforce > 0 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNUF    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'UF CCN CONCENTRATION                    '
 ENDIF

! IF ( ac_opt == 0 ) THEN
 IF ( ac_opt > 0 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNAC    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN_AC CONCENTRATION                     '
 IF ( ac_opt == 2 .or. ac_opt == 22 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNNU    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN_NU CONCENTRATION                     '
 i = i + 1
 write(f(i),2) 'CCCNCO    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN_CO CONCENTRATION                     '
 ENDIF

 ENDIF

 i = i + 1
 write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN CONCENTRATION                       '
 
 i = i + 1
 write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD CONCENTRATION                     '
 ENDIF
 
 IF ( ipconc .le. -2 ) THEN
 i = i + 1
 write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD CONCENTRATION                     '
 ENDIF
 
 IF ( ipconc .ge. 3 ) THEN
 i = i + 1
 write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '

 IF ( microphys(1:4) .eq. 'ZMRM' .or. microphys(1:5) .eq. 'ZMRHM') THEN
 i = i + 1
 write(f(i),2) 'CMW       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'MELT/SHED RAIN CONCENTRATION            '
 ENDIF

 ENDIF
 IF ( Abs(ipconc) .ge. 1 ) THEN
 i = i + 1
  IF ( imixice == 1 ) THEN ! turb mixing for ice is on
 write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'COLM ICE CONCENTRATION                  '
  ELSE  ! turb mixing for ice is off
 write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0, -1, 0, 'meters-3       ', 'COLM ICE CONCENTRATION                  '
  ENDIF

 IF ( icespheres >= 1 ) THEN ! turn on ice spheres (frozen droplets)
 i = i + 1
  IF ( imixice == 1 ) THEN ! turb mixing for ice is on
 write(f(i),2) 'CIS       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'SPHERE ICE CONCENTRATION                '
  ELSE  ! turb mixing for ice is off
 write(f(i),2) 'CIS       ', 3, 'xyz3d', 1, 0,0,0,  1,  0, -1, 0, 'meters-3       ', 'SPHERE ICE CONCENTRATION                '
  ENDIF
 ENDIF
 ENDIF

!
!  ipconc = -1  :  single moment ice and diagnose CCI (ipconc = 0)
!  ipconc = -2  :  single moment ice and diagnose CCI and CCW (ipconc = 0)
!  ipconc = -3  :  two moment ice and diagnose CCW (ipconc = 1)
!
 
 IF ( ipconc .lt. 0 ) THEN
  IF ( ipconc .eq. -3 ) THEN
    ipconc = 1
  ELSE
    ipconc = 0
  ENDIF
 ENDIF
! i = i + 1
! write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN CONCENTRATION                       '
! i = i + 1
! write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD CONCENTRATION                     '
! i = i + 1
! write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '
! i = i + 1
! write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'COLM ICE CONCENTRATION                  '
  
  IF ( ipconc .ge. 5 ) THEN
 i = i + 1
 write(f(i),2) 'CSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'SNOW CONCENTRATION                      '
 
 i = i + 1
 write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GRAUPEL CONCENTRATION                   '

 IF ( frozendrops == 1 ) THEN ! turn on frozen drops
 i = i + 1
 write(f(i),2) 'CFW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'FROZEN DROP CONCENTRATION              '
 ENDIF

  IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or. microphys(1:5) .eq. 'ZVDMH'  ) THEN
 i = i + 1
 write(f(i),2) 'CHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'HAIL CONCENTRATION                      '

 IF ( frozendrops == 2 .and. ipconc >= 8 ) THEN ! require 3-moment so that sedimentation does not need to include number correction
 i = i + 1
 write(f(i),2) 'CHF       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GRAUPEL CONCENTRATION FROM FD           '
 ENDIF


 IF ( frozendrops >= 1 .and. ipconc >= 8 ) THEN ! require 3-moment so that sedimentation does not need to include number correction
 i = i + 1
 write(f(i),2) 'CHLF      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'HAIL CONCENTRATION FROM FD              '
 ENDIF

  ENDIF
  ENDIF

  IF ( ipconc .ge. 6 ) THEN

 i = i + 1
 write(f(i),2) 'ZHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z meters-3     ', 'GRAUPEL REFLECTIVITY                    '

 IF ( frozendrops == 1 ) THEN ! turn on frozen drops
 i = i + 1
 write(f(i),2) 'ZFW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z meters-3     ', 'FROZEN DROP REFLECTIVITY                '
 ENDIF

  ENDIF

  IF ( ipconc .ge. 7 .and. ipconc < 9 ) THEN
  IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or. microphys(1:5) .eq. 'ZVDMH'  ) THEN
 i = i + 1
 write(f(i),2) 'ZHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z meters-3     ', 'HAIL REFLECTIVITY                       '
  ENDIF
  ENDIF

  IF ( ipconc .ge. 8 ) THEN
 i = i + 1
 write(f(i),2) 'ZRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z meters-3     ', 'RAIN REFLECTIVITY                       '
 IF ( microphys(1:4) .eq. 'ZMRM' .or. microphys(1:5) .eq. 'ZMRHM') THEN
 i = i + 1
 write(f(i),2) 'ZMW       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'MELT/SHED RAIN REFLECTIVITY             '
 ENDIF
  ENDIF

  IF ( ipconc .ge. 2 .and. ioutput_ssmx > 0 ) THEN
  
 i = i + 1
 write(f(i),2) 'SSMX      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'MAX SUPERSATURATION                     '
  ENDIF

! i = i + 1
! write(f(i),2) 'VHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'coulomb m-3    ', 'HAIL PARTICLE VOLUME                    '

 ENDIF
 
 IF ( microphys(1:3) .eq. 'ZVD' ) THEN

! i = i + 1
! write(f(i),2) 'VSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'SNOW PARTICLE VOLUME                    '

 i = i + 1
 write(f(i),2) 'VHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'GRAUPEL PARTICLE VOLUME                 '

 IF ( frozendrops == 1 ) THEN ! turn on frozen drops
 i = i + 1
 write(f(i),2) 'VFW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'FROZEN DROP PARTICLE VOLUME             '
 ENDIF

 ENDIF

! IF ( ( microphys(1:4) .eq. 'ZVDH' .and. .not. microphys(1:5) .eq. 'ZVDHF') .or. microphys(1:5) .eq. 'ZVDMH' ) THEN ! this line for testing
 IF ( microphys(1:4) .eq. 'ZVDH'  .or. microphys(1:5) .eq. 'ZVDMH' ) THEN
 i = i + 1
 write(f(i),2) 'VHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'HAIL PARTICLE VOLUME                    '
 ENDIF

  IF ( iraintypes > 0 .and. ( microphys(1:3) .eq. 'ZVD' .or. microphys(1:4) .eq. 'ZIEG' ) ) THEN
 i = i + 1
 write(f(i),2) 'QRAUTO    ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'RAIN FROM AUTOCONVERSION                '
 i = i + 1
 write(f(i),2) 'QRSHED    ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'RAIN FROM SHEDDING (WET GROWTH)         '
 i = i + 1
 write(f(i),2) 'QRMELT    ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'kg kg-1        ', 'RAIN FROM MELTING                       '
  nraintypes = 3
  ENDIF

 IF ( microphys(1:5) .eq. 'ICE10' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'COLM ICE MIXING RATIO                   '
 i = i + 1
 write(f(i),2) 'QIP       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'PLATE ICE MIXING RATIO                  '
 i = i + 1
 write(f(i),2) 'QIR       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RIME ICE MIXING RATIO                   '
 i = i + 1
 write(f(i),2) 'QS        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'SNOW MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QGL       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GL  MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QGM       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GM  MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QGH       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GH  MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QF        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'FD  MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QH        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'HAIL MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'LARGE HAIL MIXING RATIO                 '

 i = i + 1
 write(f(i),2) 'RTC       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'CLOUD TIME                              '
! i = i + 1
! write(f(i),2) 'RTI       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'COLM ICE TIME                           '
! i = i + 1
! write(f(i),2) 'RTIR      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'RIME ICE TIME                           '
! i = i + 1
! write(f(i),2) 'RTIP      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'PLATE ICE TIME                          '
! i = i + 1
! write(f(i),2) 'RTS       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'SNOW TIME                               '
! i = i + 1
! write(f(i),2) 'RTGL      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'GL  TIME                                '
! i = i + 1
! write(f(i),2) 'RTGM      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'GM  TIME                                '
! i = i + 1
! write(f(i),2) 'RTGH      ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'GH  TIME                                '
! i = i + 1
! write(f(i),2) 'RTF       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  1, 0, 'seconds        ', 'FD  TIME                                '

 ENDIF

 IF ( microphys(1:7) .eq. 'ICE10DM' ) THEN

 IF ( ipconc .ge. 2 ) THEN
 i = i + 1
 write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CCN CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CLOUD CONCENTRATION                     '
 ENDIF

 IF ( ipconc .ge. 3 ) THEN
 i = i + 1
 write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RAIN CONCENTRATION                      '
 ENDIF
 IF ( ipconc .ge. 1 ) THEN
 i = i + 1
 write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'COLM ICE CONCENTRATION                  '
 i = i + 1
 write(f(i),2) 'CIP       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'RIME ICE CONCENTRATION                  '
 i = i + 1
 write(f(i),2) 'CIR       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'PLATE ICE CONCENTRATION                 '
 ENDIF
 IF ( ipconc .ge. 4 ) THEN
 i = i + 1
 write(f(i),2) 'CSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'SNOW CONCENTRATION                      '
 ENDIF

 IF ( ipconc .ge. 5 ) THEN
 i = i + 1
 write(f(i),2) 'CGL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GL  CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CGM       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GM  CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CGH       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'GH  CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CFW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'FD  CONCENTRATION                       '
 i = i + 1
 write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'HAIL CONCENTRATION                      '
 i = i + 1
 write(f(i),2) 'CHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'LARGE HAIL CONCENTRATION                '
 ENDIF
  
 IF ( ioutput_ssmx > 0 ) THEN
 i = i + 1
 write(f(i),2) 'SSMX      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'count          ', 'MAX SUPERSATURATION                     '
 ENDIF
 
 ENDIF

 IF ( microphys(1:5) .eq. 'ZIEGE' .or. microphys(1:4) .eq. 'ZVDE' .or. microphys(1:5) .eq. 'ZVDHE' .or. &
      microphys(1:5) .eq. 'ZVDME' .or. microphys(1:6) .eq. 'ZVDHVE' .or. microphys(1:6) .eq. 'ZIEGHE'   &
      .or. microphys(1:6) .eq. 'ZVDHME' .or. microphys(1:6) .eq. 'ZVDMHE') THEN
  
 i = i + 1
 write(f(i),2) 'SCCW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'CLOUD CHARGE DENSITY                    '
 i = i + 1
 write(f(i),2) 'SCRW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'RAIN CHARGE DENSITY                     '
 i = i + 1
 write(f(i),2) 'SCCI      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'COLM ICE CHARGE DENSITY                 '
 IF ( icespheres >= 1 ) THEN
 i = i + 1
 write(f(i),2) 'SCIS      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'SPHERE ICE CHARGE DENSITY               '
 ENDIF
 i = i + 1
 write(f(i),2) 'SCSW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'SNOW CHARGE DENSITY                     '
 i = i + 1
 write(f(i),2) 'SCHW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'GRAUPEL CHARGE DENSITY                  '

 IF ( frozendrops == 1 ) THEN
 i = i + 1
 write(f(i),2) 'SCFW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'FD  CHARGE DENSITY                      '
 ENDIF
 
 IF ( microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH' .or.  microphys(1:5) .eq. 'ZVDMH' ) THEN
 i = i + 1
 write(f(i),2) 'SCHL      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'HAIL CHARGE DENSITY                     '
 ENDIF
 i = i + 1
 write(f(i),2) 'CPION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'POS ION CONCENTRATION                   '
 i = i + 1
 write(f(i),2) 'CNION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NEG ION CONCENTRATION                   '

 IF ( .false. ) THEN
 i = i + 1
 write(f(i),2) 'CPLION    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'POS LARGE ION CONCENTRATION             '
 i = i + 1
 write(f(i),2) 'CNLION    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NEG LARGE ION CONCENTRATION             '
 ENDIF
 
  IF ( ichaff == 1 ) THEN

 i = i + 1
 write(f(i),2) 'CCHAFF    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'CHAFF CONCENTRATION                     '

  ELSEIF ( ichaff == 2 ) THEN

 i = i + 1
 write(f(i),2) 'CLNOX     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO             '

  ELSEIF ( ichaff == 3 ) THEN

 i = i + 1
 write(f(i),2) 'CLNOX     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO             '
 i = i + 1
 write(f(i),2) 'CLNOX_A   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB A       '
 i = i + 1
 write(f(i),2) 'CLNOX_B   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB B       '
 i = i + 1
 write(f(i),2) 'CLNOX_C   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB C       '
 i = i + 1
 write(f(i),2) 'CLNOX_D   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB D       '
 i = i + 1
 write(f(i),2) 'CLNOX_E   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB E       '
 i = i + 1
 write(f(i),2) 'CLNOX_F   ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'LNOX MOLECULAR MIXING RATIO SLAB F       '
 i = i + 1
 write(f(i),2) 'CANOX     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'AMBIENT NOX MIXING RATIO                 '
 i = i + 1
 write(f(i),2) 'CCO       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 0, 'mol mol-1      ', 'AMBIENT CO MIXING RATIO                  '

  ENDIF
 ENDIF

 IF ( microphys(1:8) .eq. 'ICE10DME' .or. microphys(1:6) .eq. 'ICE10E' ) THEN
  
 i = i + 1
 write(f(i),2) 'SCCW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'CLOUD CHARGE DENSITY                    '
 i = i + 1
 write(f(i),2) 'SCRW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'RAIN CHARGE DENSITY                     '
 i = i + 1
 write(f(i),2) 'SCCI      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'COLM ICE CHARGE DENSITY                 '
 i = i + 1
 write(f(i),2) 'SCIR      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'RIME ICE CHARGE DENSITY                 '
 i = i + 1
 write(f(i),2) 'SCIP      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'PLATE ICE CHARGE DENSITY                '
 i = i + 1
 write(f(i),2) 'SCSW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'SNOW CHARGE DENSITY                     '
 i = i + 1
 write(f(i),2) 'SCGL      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'GL  CHARGE DENSITY                      '
 i = i + 1
 write(f(i),2) 'SCGM      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'GM  CHARGE DENSITY                      '
 i = i + 1
 write(f(i),2) 'SCGH      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'GH  CHARGE DENSITY                      '
 i = i + 1
 write(f(i),2) 'SCFW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'FD  CHARGE DENSITY                      '
 i = i + 1
 write(f(i),2) 'SCHW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'HAIL CHARGE DENSITY                     '
 i = i + 1
 write(f(i),2) 'SCHL      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'LARGE HAIL CHARGE DENSITY               '
 i = i + 1
 write(f(i),2) 'CPION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'POS ION CONCENTRATION                   '
 i = i + 1
 write(f(i),2) 'CNION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NEG IONCONCENTRATION                   '

 ENDIF
 
 ! DTD: Added variables for MY scheme; mostly re-used existing ZVD variable names where appropriate
 ! Also, concentration and reflectivity variables are scaled by density within scheme, so they
 ! are passed back to the model dynamics already in "per mass" units
 ! UPDATE 09/27/2011: MY concentration and reflectivity now scaled in dynamics just like ZVD
 
 IF ( microphys(1:2) .eq. 'MY' ) THEN
  i = i + 1
  write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
  i = i + 1
  write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
  i = i + 1
  write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '
  i = i + 1
  write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ICE MIXING RATIO                        '
  i = i + 1
  write(f(i),2) 'QS        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'SNOW MIXING RATIO                       '
  i = i + 1
  write(f(i),2) 'QH        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GRAUPEL MIXING RATIO                    '
  i = i + 1
  write(f(i),2) 'QHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'HAIL MIXING RATIO                       '
  
  IF ( microphys(1:3) .eq. 'MY2' .or. microphys(1:3) .eq. 'MY3' ) THEN
    i = i + 1
    write(f(i),2) 'CCW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'CLOUD CONCENTRATION                     '
    i = i + 1
    write(f(i),2) 'CRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'RAIN CONCENTRATION                      '
    i = i + 1
    write(f(i),2) 'CCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'ICE CONCENTRATION                       '
    i = i + 1
    write(f(i),2) 'CSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'SNOW CONCENTRATION                      '
    i = i + 1
    write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'GRAUPEL CONCENTRATION                   '
    i = i + 1
    write(f(i),2) 'CHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'HAIL CONCENTRATION                      '
  END IF
  
  IF ( microphys(1:3) .eq. 'MY3' ) THEN
    i = i + 1
    write(f(i),2) 'ZRW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z kg-1        ', 'RAIN REFLECTIVITY                       '
    i = i + 1
    write(f(i),2) 'ZCI       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z kg-1        ', 'ICE REFLECTIVITY                        '
    i = i + 1
    write(f(i),2) 'ZSW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z kg-1        ', 'SNOW REFLECTIVITY                       '
    i = i + 1
    write(f(i),2) 'ZHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z kg-1        ', 'GRAUPEL REFLECTIVITY                    '
    i = i + 1
    write(f(i),2) 'ZHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'Z kg-1        ', 'HAIL REFLECTIVITY                       '
  END IF
  
 ENDIF

 IF ( microphys(1:3) .eq. 'HCM' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ICE MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QS        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'SNOW MIXING RATIO                       '
 
 DO n = 1,nch
 write(tmpnum,'(i2.2)') n
!  write(0,*) 'tmpnum = ',tmpnum
  IF ( n == 1 ) THEN
   tmpname = 'QH'
  ELSE
   tmpname = 'QH'//trim(tmpnum)
  ENDIF
!  write(0,*) 'tmpname = ',tmpname
 tmpdescript = 'HCM GRAUPEL SIZE '//trim(tmpnum)
 i = i + 1
 write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
 ENDDO

  IF ( microphys(1:4) .eq. 'HCMV' ) THEN

  DO n = 1,nch
   write(tmpnum,'(i2.2)') n
    IF ( n == 1 ) THEN
     tmpname = 'VHW'
    ELSE
     tmpname = 'VHW'//tmpnum
    ENDIF
    tmpdescript = 'HCM GRAUPEL VOLUME '//tmpnum
    i = i + 1
!    write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  tmpdescript ! 'GRAUPEL MIXING RATIO                    '
   write(f(i),2) tmpname, 3, 'xyz3n', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ',  trim(tmpdescript) !'GRAUPEL PARTICLE VOLUME                 '
   ENDDO

  ENDIF

  IF ( microphys(1:4) .eq. 'HCMD' ) THEN

  DO n = 1,nch
   write(tmpnum,'(i2.2)') n
    IF ( n == 1 ) THEN
     tmpname = 'VHW'
    ELSE
     tmpname = 'VHW'//tmpnum
    ENDIF
    tmpdescript = 'HCM GRAUPEL DENSITY '//tmpnum
    i = i + 1
!    write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  tmpdescript ! 'GRAUPEL MIXING RATIO                    '
   write(f(i),2) tmpname, 3, 'xyz3n', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ',  trim(tmpdescript) !'GRAUPEL PARTICLE DENSITY                 '
   ENDDO

  ENDIF

     IF ( .false. ) THEN ! testing for how many variables before parallel I/O fails
      ! DROPLETS/RAIN
      DO n = 1,ntakrd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CRW'
        ELSE
         tmpname = 'CRW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK RAIN BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
! write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'GRAUPEL CONCENTRATION                   '
       ENDDO

      ! ICE CRYSTALS
      DO m = 1,ntakit
       write(tmpdigit,'(i1)') m
      DO n = 1,ntakid
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 .and. m == 1 ) THEN
         tmpname = 'CCI'
        ELSE
         tmpname = 'CCI'//trim(tmpdigit)//'_'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK CRYSTAL BIN '//trim(tmpdigit)//'_'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       ENDDO

      IF ( .false. ) THEN
      ! GRAUPEL
      DO n = 1,ntakpd
!      DO n = 1,5
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CHW'
        ELSE
         tmpname = 'CHW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK GRAUPEL BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!      write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
      ENDIF ! .false.
      
      IF ( .true. ) THEN
      ! HAIL
!      DO n = 1,ntakpd
      DO n = 1,10 ! ntakpd/2
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CHL'
        ELSE
         tmpname = 'CHL'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK HAIL BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!      write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       ENDIF ! .false
       
      ENDIF ! true/false

 
 ENDIF

 IF ( iashtypes > 0 ) THEN
 i = i + 1
 write(f(i),2) 'QASH1     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ASH TYPE 1                              '
   IF ( iashtypes > 1 ) THEN
 i = i + 1
 write(f(i),2) 'QASH2     ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ASH TYPE 2                              '
   ENDIF
 ENDIF

 IF ( microphys(1:3) .eq. 'TAK' ) THEN
   
 i = i + 1
 write(f(i),2) 'QV        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 2, 'kg kg-1        ', 'VAPOR MIXING RATIO                      '
! QC -- QHL are actually diagnostic, but are here for buoyancy calculations. What's a few more scalars, eh?
 i = i + 1
 write(f(i),2) 'QC        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'CLOUD MIXING RATIO                      '
 i = i + 1
 write(f(i),2) 'QR        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'RAIN MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QI        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'ICE MIXING RATIO                        '
 i = i + 1
 write(f(i),2) 'QS        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'SNOW MIXING RATIO                       '
 i = i + 1
 write(f(i),2) 'QH        ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'GRAUPEL MIXING RATIO                    '
 IF ( takkfmax >= 2 ) THEN
 i = i + 1
 write(f(i),2) 'QHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ', 'HAIL MIXING RATIO                       '
 ENDIF
 
 IF ( microphys(1:5) .eq. 'TAKHV' ) THEN
 ! do not increment 'i' because we want QHL to be overwritten
 ! But for now add volume for both graupel and hail
! i = i + 1
! write(f(i),2) 'VHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'GRAUPEL PARTICLE VOLUME                 '
! i = i + 1
! write(f(i),2) 'VHL       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ', 'HAIL PARTICLE VOLUME                    '
 ENDIF

!   IF ( microphys(1:4) .eq. 'TAKR' .or. microphys(1:4) .eq. 'TAKH' ) THEN

      ! DROPLETS/RAIN
      DO n = 1,ntakrd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CRW'
        ELSE
         tmpname = 'CRW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK RAIN BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
! write(f(i),2) 'CHW       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'kg-1          ', 'GRAUPEL CONCENTRATION                   '
       ENDDO

!   ENDIF


!  IF ( microphys(1:4) .eq. 'TAKH' ) THEN

      ! ICE CRYSTALS
      DO m = 1,ntakit
       write(tmpdigit,'(i1)') m
      DO n = 1,ntakid
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 .and. m == 1 ) THEN
         tmpname = 'CCI'
        ELSE
         tmpname = 'CCI'//trim(tmpdigit)//'_'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK CRYSTAL BIN '//trim(tmpdigit)//'_'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       ENDDO
      
      
      ! GRAUPEL
      DO n = 1,ntakpd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CHW'
        ELSE
         tmpname = 'CHW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK GRAUPEL BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!      write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
      
      
      ! HAIL
      IF ( takkfmax >= 2 ) THEN
      DO n = 1,ntakpd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'CHL'
        ELSE
         tmpname = 'CHL'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK HAIL BIN '//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
!      write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       ENDIF
      
!       write(0,*) 'microphys = ',microphys
        IF ( microphys(1:5) .eq. 'TAKHV' ) THEN
!       write(0,*) 'TAKHV: set up volumes i = ',i
        DO n = 1,ntakpd
         write(tmpnum,'(i2.2)') n
          IF ( n == 1 ) THEN
           tmpname = 'VHW'
          ELSE
           tmpname = 'VHW'//tmpnum
          ENDIF
          tmpdescript = 'TAKHV GRAUPEL VOLUME '//tmpnum
          i = i + 1
      !    write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  tmpdescript ! 'GRAUPEL MIXING RATIO                    '
         write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ',  trim(tmpdescript) !'GRAUPEL PARTICLE VOLUME                 '
         ENDDO

        IF ( takkfmax >= 2 ) THEN
        DO n = 1,ntakpd
         write(tmpnum,'(i2.2)') n
          IF ( n == 1 ) THEN
           tmpname = 'VHL'
          ELSE
           tmpname = 'VHL'//tmpnum
          ENDIF
          tmpdescript = 'TAKHV HAIL VOLUME '//tmpnum
          i = i + 1
      !    write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  1,  0, 3, 'kg kg-1        ',  tmpdescript ! 'GRAUPEL MIXING RATIO                    '
         write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'm3 m-3         ',  trim(tmpdescript) !'GRAUPEL PARTICLE VOLUME                 '
         ENDDO
         ENDIF
      
!       write(0,*) 'TAKHV: set up volumes i = ',i

        ENDIF  ! TAKHV


      IF ( elec ) THEN
! CHARGE
      
      ! DROPLETS/RAIN
      DO n = 1,ntakrd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'SCRW'
        ELSE
         tmpname = 'SCRW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK RAIN BIN CHARGE'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
! write(f(i),2) 'SCHW      ', 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ', 'GRAUPEL CHARGE DENSITY                  '
       ENDDO

!   ENDIF


!  IF ( microphys(1:4) .eq. 'TAKH' ) THEN

      ! ICE CRYSTALS
      DO m = 1,ntakit
       write(tmpdigit,'(i1)') m
      DO n = 1,ntakid
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 .and. m == 1 ) THEN
         tmpname = 'SCCI'
        ELSE
         tmpname = 'SCCI'//trim(tmpdigit)//'_'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK CRYSTAL BIN CHARGE'//trim(tmpdigit)//'_'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       ENDDO
      
      
      ! GRAUPEL
      DO n = 1,ntakpd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'SCHW'
        ELSE
         tmpname = 'SCHW'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK GRAUPEL BIN CHARGE'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
      
      
      ! HAIL
      IF ( takkfmax >= 2 ) THEN
      DO n = 1,ntakpd
       write(tmpnum,'(i2.2)') n
      !  write(0,*) 'tmpnum = ',tmpnum
        IF ( n == 1 ) THEN
         tmpname = 'SCHL'
        ELSE
         tmpname = 'SCHL'//trim(tmpnum)
        ENDIF
      !  write(0,*) 'tmpname = ',tmpname
       tmpdescript = 'TAK HAIL BIN CHARGE'//trim(tmpnum)
       i = i + 1
       write(f(i),2) tmpname, 3, 'xyz3d', 1, 0,0,0,  0,  0,  0, 0, 'coulomb m-3    ',  trim(tmpdescript) ! 'GRAUPEL MIXING RATIO                    '
       ENDDO
       
       ENDIF
      
 i = i + 1
 write(f(i),2) 'CPION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'POS ION CONCENTRATION                   '
 i = i + 1
 write(f(i),2) 'CNION     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'NEG ION CONCENTRATION                   '

      ENDIF ! elec

       
!     ENDIF ! TAKH

 IF ( ccnuf > 1.0 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNUF    ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'meters-3       ', 'UF CCN CONCENTRATION                    '
 ENDIF
 
  i = i + 1
 write(f(i),2) 'CCCN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'CCN CONCENTRATION                       '

 i = i + 1
 write(f(i),2) 'CCIN      ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'CLOUD ICE NUCLEI CONCENTRATION          '

 i = i + 1
 write(f(i),2) 'CCINA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'NUMBER OF ACTIVATED ICE NUCLEI          '

 IF ( takcwnucopt >= 3 ) THEN
 i = i + 1
 write(f(i),2) 'CCCNA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'NUMBER OF ACTIVATED CCN                 '
 ENDIF
! i = i + 1
! write(f(i),2) 'CCCNA     ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'NUMBER OF ACTIVATED CCN                 '

 i = i + 1
 write(f(i),2) 'CNG       ', 3, 'xyz3d', 1, 0,0,0,  1,  0,  0, 0, 'cm-3           ', 'HALLETT-MOSSOP ICE CONCENTRATION        '

 IF ( ioutput_ssmx > 0  ) THEN
 i = i + 1
 write(f(i),2) 'SSMX      ', 3, 'xyz3d', 1, 0,0,0,  0,  1,  0, 0, 'percent        ', 'MAX SUPERSATURATION                     '
 ENDIF
 
 ENDIF ! if TAK


 
 i = i + 1 ! need one more for EOF
 write(f(i),'(a)') 'EOF'
 
 numlines = i
 
! close(51)
! rewind(51)
  
! CALL GRID_DEFINE_FROM_FILE(gd,internal_filename)
 CALL GRID_DEFINE_FROM_FILE(gd)

!   close(51)
 
 CALL INDEX_INIT(gd,microphys)

! Write out grid registry for conviences

 IF ( my_rank == 0 ) THEN
 open(61, file='ncommas.registry', form = 'formatted', status='unknown')
 DO i = 1, numlines
   write(61,FMT='(132a)') f(i)
 enddo
 close(61)
 ENDIF
    

 END SUBROUTINE GRID_DEFINE_FROM_LIST

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
! 

 SUBROUTINE INDEX_INIT(gd,microphys)

   USE INDEX_MODULE
   
   implicit none
   
   type(GRID)                      :: gd
   character(LEN=*), INTENT(INOUT) :: microphys
   
!  LOCAL VARS
   integer   :: s, ns, ex, wz
   
! ##########################################################
! 
     ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
     s  = GET_VARIABLE_INDEX(gd, 'TH') - 1
     ex  = GET_VARIABLE_INDEX(gd, 'EX') - 1
     wz  = GET_VARIABLE_INDEX(gd, 'WZ')

! check for extra diagnostic electrification vars
     icgaddl  = Max( 0, GET_VARIABLE_INDEX(gd, 'RSCADDL', 1) - ex )
     icgnoliq = Max( 0, GET_VARIABLE_INDEX(gd, 'RSCNOLIQ', 1) - ex )

     iex2  = Max( 0, GET_VARIABLE_INDEX(gd, 'EX2', 1) - ex )
     iey2  = Max( 0, GET_VARIABLE_INDEX(gd, 'EY2', 1) - ex )
     iez2  = Max( 0, GET_VARIABLE_INDEX(gd, 'EZ2', 1) - ex )

!     icghis  = Max( 0, GET_VARIABLE_INDEX(gd, 'RSCGHIS', 1) - ex ) 

     ichgtndadv = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDADV', 1) - ex )
     ichgtndadvsn= Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDADVSN', 1) - ex )
     ichgtndmix = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDMIX', 1) - ex )
     ichgtndnic = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDNIC', 1) - ex )
     ichgtndsed = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDSED', 1) - ex )
     ichgtndion = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDION', 1) - ex )
     ichgtndlgt = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDLGT', 1) - ex )
     ichgtndave = Max( 0, GET_VARIABLE_INDEX(gd, 'SCTNDAVG', 1) - ex )

     IF ( iez2 > 1      ) neelec = Max( neelec, iez2 )
     IF ( icgnoliq > 1  ) neelec = Max( neelec, icgnoliq )
     IF ( icgaddl > 1   ) neelec = Max( neelec, icgaddl )
     IF ( ichgtndave > 1) neelec = Max( neelec, ichgtndave )
!     write(0,*) 'icgaddl,icgnoliq, neelec = ',icgaddl,icgnoliq,neelec

! check for vars for extra diagnostic array
     lsat  = Max( 0, GET_VARIABLE_INDEX(gd, 'SSAT', 1) - wz )
     lsati = Max( 0, GET_VARIABLE_INDEX(gd, 'SSATI', 1) - wz )
     lppert = Max( 0, GET_VARIABLE_INDEX(gd, 'PPERT', 1) - wz )

     lqaut    = Max( 0, GET_VARIABLE_INDEX(gd, 'QAUT', 1) - wz )
     lqacc    = Max( 0, GET_VARIABLE_INDEX(gd, 'QACC', 1) - wz )
     lqrevap  = Max( 0, GET_VARIABLE_INDEX(gd, 'QREVAP', 1) - wz )
     lqcevap  = Max( 0, GET_VARIABLE_INDEX(gd, 'QCEVAP', 1) - wz )
     lqcond   = Max( 0, GET_VARIABLE_INDEX(gd, 'QCOND', 1) - wz )
     lqdep    = Max( 0, GET_VARIABLE_INDEX(gd, 'QDEP', 1) - wz )
     lqmelt   = Max( 0, GET_VARIABLE_INDEX(gd, 'QMELT', 1) - wz )
     lqsub    = Max( 0, GET_VARIABLE_INDEX(gd, 'QSUB', 1) - wz )

     lqh2hl    = Max( 0, GET_VARIABLE_INDEX(gd, 'QH2HL', 1) - wz )
     lqhl2h    = Max( 0, GET_VARIABLE_INDEX(gd, 'QHL2H', 1) - wz )

     ld0shdrate    = Max( 0, GET_VARIABLE_INDEX(gd, 'D0SHED', 1) - wz )
     ld1shdrate    = Max( 0, GET_VARIABLE_INDEX(gd, 'D1SHED', 1) - wz )
     ld2shdrate    = Max( 0, GET_VARIABLE_INDEX(gd, 'D2SHED', 1) - wz )
     ld3shdrate    = Max( 0, GET_VARIABLE_INDEX(gd, 'D3SHED', 1) - wz )

     lmvr  = Max( 0, GET_VARIABLE_INDEX(gd, 'VR', 1) - wz )
     lmvh  = Max( 0, GET_VARIABLE_INDEX(gd, 'VG', 1) - wz )
     lmvhl = Max( 0, GET_VARIABLE_INDEX(gd, 'VH', 1) - wz )
     lmvf  = Max( 0, GET_VARIABLE_INDEX(gd, 'VF', 1) - wz )
     lmvi  = Max( 0, GET_VARIABLE_INDEX(gd, 'VI', 1) - wz )
     lmvs  = Max( 0, GET_VARIABLE_INDEX(gd, 'VS', 1) - wz )
     ld0   = Max( 0, GET_VARIABLE_INDEX(gd, 'D0', 1) - wz )
     
     lnctak = Max( 0, GET_VARIABLE_INDEX(gd, 'CCWTAK', 1) - wz )
     lnitak = Max( 0, GET_VARIABLE_INDEX(gd, 'CCITAK', 1) - wz )
     lnrtak = Max( 0, GET_VARIABLE_INDEX(gd, 'CRWTAK', 1) - wz )
     lnhtak = Max( 0, GET_VARIABLE_INDEX(gd, 'CHWTAK', 1) - wz )
     lnhltak = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLTAK', 1) - wz )

     lash1 = Max( 0, GET_VARIABLE_INDEX(gd, 'QASH1', 1) - wz )
     lash2 = Max( 0, GET_VARIABLE_INDEX(gd, 'QASH2', 1) - wz )
     lash3 = Max( 0, GET_VARIABLE_INDEX(gd, 'QASH3', 1) - wz )


     lnvhl = Max( 0, GET_VARIABLE_INDEX(gd, 'VNH', 1) - wz )
     lzvhl = Max( 0, GET_VARIABLE_INDEX(gd, 'VZH', 1) - wz )
     lmnvhl = Max( 0, GET_VARIABLE_INDEX(gd, 'VMNHL', 1) - wz )
     ialphahl = Max( 0, GET_VARIABLE_INDEX(gd, 'ALPHAHL3D', 1) - wz )
     idfw = Max( 0, GET_VARIABLE_INDEX(gd, 'FWDIA', 1) - wz )
     idhw = Max( 0, GET_VARIABLE_INDEX(gd, 'HWDIA', 1) - wz )
     idhl = Max( 0, GET_VARIABLE_INDEX(gd, 'HLDIA', 1) - wz )
     idmhl = Max( 0, GET_VARIABLE_INDEX(gd, 'HLMWDIA', 1) - wz )
     idnhl = Max( 0, GET_VARIABLE_INDEX(gd, 'HLCHARDIA', 1) - wz )
     ialphah = Max( 0, GET_VARIABLE_INDEX(gd, 'ALPHAH3D', 1) - wz )
     ialphaf = Max( 0, GET_VARIABLE_INDEX(gd, 'ALPHAF3D', 1) - wz )
     ialphar = Max( 0, GET_VARIABLE_INDEX(gd, 'ALPHAR3D', 1) - wz )
     ldbzchangeh = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGEH', 1) - wz )
     ldbzchanger = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGER', 1) - wz )

     lchlcnh = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNH', 1) - wz )
     lchlcnf = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNF', 1) - wz )
     ldhlcnh = Max( 0, GET_VARIABLE_INDEX(gd, 'DHLCNH', 1) - wz )
     ldhlcnf = Max( 0, GET_VARIABLE_INDEX(gd, 'DHLCNF', 1) - wz )

     ldbzr  = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZR', 1) - wz )
     ldbzh  = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZH', 1) - wz )
     ldbzhl = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZHL', 1) - wz )
     ldbzi  = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZI', 1) - wz )
     lwvzf  = Max( 0, GET_VARIABLE_INDEX(gd, 'WVZF', 1) - wz )

     lkmt  = Max( 0, GET_VARIABLE_INDEX(gd, 'KMT', 1) - wz )
     lkht  = Max( 0, GET_VARIABLE_INDEX(gd, 'KHT', 1) - wz )
     ltkediss = Max( 0, GET_VARIABLE_INDEX(gd, 'TKE_DISS', 1) - wz )

     lairpress = Max( 0, GET_VARIABLE_INDEX(gd, 'P', 1) - wz )
     lairtem   = Max( 0, GET_VARIABLE_INDEX(gd, 'T', 1) - wz )
     lthetav   = Max( 0, GET_VARIABLE_INDEX(gd, 'THETA_V', 1) - wz )
     lhwind    = Max( 0, GET_VARIABLE_INDEX(gd, 'HWIND', 1) - wz )
!     lssw3d    = Max( 0, GET_VARIABLE_INDEX(gd, 'SSAT', 1) - wz )
!     lssi3d    = Max( 0, GET_VARIABLE_INDEX(gd, 'SSATI', 1) - wz )
     lrh3d    = Max( 0, GET_VARIABLE_INDEX(gd, 'RH', 1) - wz )
     lcwdia    = Max( 0, GET_VARIABLE_INDEX(gd, 'CWDIA', 1) - wz )
     lcwnu    = Max( 0, GET_VARIABLE_INDEX(gd, 'CWNU', 1) - wz )
     lswdia    = Max( 0, GET_VARIABLE_INDEX(gd, 'SWDIA', 1) - wz )
     lswmass   = Max( 0, GET_VARIABLE_INDEX(gd, 'SWMASS', 1) - wz )
     lswdn     = Max( 0, GET_VARIABLE_INDEX(gd, 'SWDN',  1) - wz )
     lhwdn     = Max( 0, GET_VARIABLE_INDEX(gd, 'HWDN',  1) - wz )
     lhldn     = Max( 0, GET_VARIABLE_INDEX(gd, 'HLDN',  1) - wz )
     lfwdn     = Max( 0, GET_VARIABLE_INDEX(gd, 'FWDN',  1) - wz )
     lswagg    = Max( 0, GET_VARIABLE_INDEX(gd, 'SWAGG', 1) - wz )
!     write(0,*) 'lswagg,lswdn,lswdia = ',lswagg,lswdn,lswdia

     iprecip = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_RAT', 1)  )

     iprainacc = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACC', 1) - iprecip + 1  )
     iprainacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACC2', 1) - iprecip + 1  )

     iphailacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACCH', 1) - iprecip + 1  )
     iphailacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACCH2', 1) - iprecip + 1  )
     iphailfacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACCHF', 1) - iprecip + 1  )
     iphailfacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACHF2', 1) - iprecip + 1  )
     iphailnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACC', 1) - iprecip + 1  )
     iphailnumacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACC2', 1) - iprecip + 1  )
     iphailfnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACCF', 1) - iprecip + 1  )
     iphailfnumacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACF2', 1) - iprecip + 1  )
     iphaildiam = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_DIAM', 1) - iprecip + 1  )

     iprainaccauto = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACC_A', 1) - iprecip + 1  )
     iprainaccmelt = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACC_M', 1) - iprecip + 1  )
     iprainaccshed = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACC_S', 1) - iprecip + 1  )
     iprainaccauto2 = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACCA2', 1) - iprecip + 1  )
     iprainaccmelt2 = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACCM2', 1) - iprecip + 1  )
     iprainaccshed2 = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_ACCS2', 1) - iprecip + 1  )

     ipfdauto = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_AUTO', 1) - iprecip + 1  )
     ipfdmelt = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_MELT', 1) - iprecip + 1  )
     ipfdshed = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_SHED', 1) - iprecip + 1  )
     ipfdauto2 = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_AUTO2', 1) - iprecip + 1  )
     ipfdmelt2 = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_MELT2', 1) - iprecip + 1  )
     ipfdshed2 = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_SHED2', 1) - iprecip + 1  )

     ipgracc = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_ACC', 1) - iprecip + 1  )
     ipgracc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_ACC2', 1) - iprecip + 1  )
     ipgrnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_NACC', 1) - iprecip + 1  )
     ipgrnumacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_NACC2', 1) - iprecip + 1  )
     ipgrdiam = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_DIAM', 1) - iprecip + 1  )

     ipfdacc = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_ACC', 1) - iprecip + 1  )
     ipfdnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_NACC', 1) - iprecip + 1  )
     ipfdacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_ACC2', 1) - iprecip + 1  )
     ipfdnumacc2 = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_NACC2', 1) - iprecip + 1  )
     ipfddiam = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_DIAM', 1) - iprecip + 1  )


     ipwfrat = Max( 0, GET_VARIABLE_INDEX(gd, 'PWF_RAT', 1) - iprecip + 1  )
     ipwfacc = Max( 0, GET_VARIABLE_INDEX(gd, 'PWF_ACC', 1) - iprecip + 1  )
     idbzcomp = Max( 0, GET_VARIABLE_INDEX(gd, 'COMP_DBZ', 1)   - iprecip + 1 )
     ihwindsfcmx = Max( 0, GET_VARIABLE_INDEX(gd, 'UV_SFC_MAX', 1)   - iprecip + 1 )
     ihailmax2d = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_MAX2D', 1)   - iprecip + 1 )
     ihailmaxk1 = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_MAXK1', 1)   - iprecip + 1 )
     iwzsfcmax = Max( 0, GET_VARIABLE_INDEX(gd, 'WZ_SFC_MAX', 1)   - iprecip + 1 )
     iwzsfcmin = Max( 0, GET_VARIABLE_INDEX(gd, 'WZ_SFC_MIN', 1)   - iprecip + 1 )
     iflshr  = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHR', 1)   - iprecip + 1 )
     iflshsrc  = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHSRC', 1)   - iprecip + 1 )
     iflshr8km = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSH8KM', 1)   - iprecip + 1 )
     iflshfed = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFED', 1)   - iprecip + 1 )
     iflshfod = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFOD', 1)   - iprecip + 1 )
     iflshfodic = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFODIC', 1)   - iprecip + 1 )
     iflshfodcgn = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFODCGN', 1)   - iprecip + 1 )
     iflshfodcgp = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFODCGP', 1)   - iprecip + 1 )
     iflshfedic = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFEDIC', 1)   - iprecip + 1 )
     iflshfedicp = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFEDICP', 1)   - iprecip + 1 )
     iflshfedicn = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFEDICN', 1)   - iprecip + 1 )
     iflshfedcgn = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFEDCGN', 1)   - iprecip + 1 )
     iflshfedcgp = Max( 0, GET_VARIABLE_INDEX(gd, 'FLSHFEDCGP', 1)   - iprecip + 1 )
     ipnic2d = Max( 0, GET_VARIABLE_INDEX(gd, 'PNIC_2D', 1) - iprecip + 1 )
     innic2d = Max( 0, GET_VARIABLE_INDEX(gd, 'NNIC_2D', 1) - iprecip + 1 )
     
!     write(0,*) 'ipnic2d,innic2d = ',ipnic2d,innic2d
!     write(0,*) 'iflshfedcgp = ',iflshfedcgp
     
! Find vars in scalar array

! Theta
!     lt = Max( 0, GET_VARIABLE_INDEX(gd, 'TH', 1) - s )

! Mixing Ratios
!     lv = Max( 0, GET_VARIABLE_INDEX(gd, 'QV', 1) - s )
!     lc = Max( 0, GET_VARIABLE_INDEX(gd, 'QC', 1) - s )
     lr = Max( 0, GET_VARIABLE_INDEX(gd, 'QR', 1) - s )
     lm = Max( 0, GET_VARIABLE_INDEX(gd, 'QM', 1) - s )
     li = Max( 0, GET_VARIABLE_INDEX(gd, 'QI', 1) - s )
     lis = Max( 0, GET_VARIABLE_INDEX(gd, 'QIS', 1) - s )
     ls = Max( 0, GET_VARIABLE_INDEX(gd, 'QS', 1) - s )
     lh = Max( 0, GET_VARIABLE_INDEX(gd, 'QH', 1) - s )

     lir = Max( 0, GET_VARIABLE_INDEX(gd, 'QIR', 1) - s )
     lip = Max( 0, GET_VARIABLE_INDEX(gd, 'QIP', 1) - s )
     lgl = Max( 0, GET_VARIABLE_INDEX(gd, 'QGL', 1) - s )
     lgm = Max( 0, GET_VARIABLE_INDEX(gd, 'QGM', 1) - s )
     lgh = Max( 0, GET_VARIABLE_INDEX(gd, 'QGH', 1) - s )
     lf  = Max( 0, GET_VARIABLE_INDEX(gd, 'QF',  1) - s )
     lhl = Max( 0, GET_VARIABLE_INDEX(gd, 'QHL', 1) - s )

     lrtc = Max( 0, GET_VARIABLE_INDEX(gd, 'RTC', 1) - s )
     lrti = Max( 0, GET_VARIABLE_INDEX(gd, 'RTI', 1) - s )
     lrts = Max( 0, GET_VARIABLE_INDEX(gd, 'RTS', 1) - s )
     lrth = Max( 0, GET_VARIABLE_INDEX(gd, 'RTH', 1) - s )

     lrtir = Max( 0, GET_VARIABLE_INDEX(gd, 'RTIR', 1) - s )
     lrtip = Max( 0, GET_VARIABLE_INDEX(gd, 'RTIP', 1) - s )
     lrtgl = Max( 0, GET_VARIABLE_INDEX(gd, 'RTGL', 1) - s )
     lrtgm = Max( 0, GET_VARIABLE_INDEX(gd, 'RTGM', 1) - s )
     lrtgh = Max( 0, GET_VARIABLE_INDEX(gd, 'RTGH', 1) - s )
     lrtf  = Max( 0, GET_VARIABLE_INDEX(gd, 'RTF',  1) - s )
     lrthl = Max( 0, GET_VARIABLE_INDEX(gd, 'RTHL', 1) - s )

! Mixed phase hydrometeors
     lsw  = Max( 0, GET_VARIABLE_INDEX(gd, 'QSW', 1) - s )
     lhw  = Max( 0, GET_VARIABLE_INDEX(gd, 'QHW', 1) - s )
     lfw  = Max( 0, GET_VARIABLE_INDEX(gd, 'QFW', 1) - s )
     lhlw = Max( 0, GET_VARIABLE_INDEX(gd, 'QHLW', 1) - s )
     lhwlg  = Max( 0, GET_VARIABLE_INDEX(gd, 'QHWLG', 1) - s )
     lhlwlg = Max( 0, GET_VARIABLE_INDEX(gd, 'QHLWLG', 1) - s )
     lrauto = Max( 0, GET_VARIABLE_INDEX(gd, 'QRAUTO', 1) - s )
     lrshed = Max( 0, GET_VARIABLE_INDEX(gd, 'QRSHED', 1) - s )
     lrmelt = Max( 0, GET_VARIABLE_INDEX(gd, 'QRMELT', 1) - s )

! Concentrations
     lccna = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNA',1) - s )
     lccnaco = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNACO',1) - s )
     lccnanu = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNANU',1) - s )
     lccn = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCN',1) - s )
     lcn_ac = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNAC',1) - s )
     lcn_nu = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNNU',1) - s )
     lcn_co = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNCO',1) - s )
     lccnuf = Max( 0, GET_VARIABLE_INDEX(gd, 'CCCNUF',1) - s )
     lcin = Max( 0, GET_VARIABLE_INDEX(gd, 'CCIN',1) - s )
     lcina = Max( 0, GET_VARIABLE_INDEX(gd, 'CCINA',1) - s )
     lcng = Max( 0, GET_VARIABLE_INDEX(gd, 'CNG',1) - s )
     lnc  = Max( 0, GET_VARIABLE_INDEX(gd, 'CCW', 1) - s )
     lnr  = Max( 0, GET_VARIABLE_INDEX(gd, 'CRW', 1) - s )
     lnm  = Max( 0, GET_VARIABLE_INDEX(gd, 'CMW', 1) - s )
     lni  = Max( 0, GET_VARIABLE_INDEX(gd, 'CCI', 1) - s )
     lnis = Max( 0, GET_VARIABLE_INDEX(gd, 'CIS', 1) - s )
     lns  = Max( 0, GET_VARIABLE_INDEX(gd, 'CSW', 1) - s )
     lnh  = Max( 0, GET_VARIABLE_INDEX(gd, 'CHW', 1) - s )

     lnir = Max( 0, GET_VARIABLE_INDEX(gd, 'CIR', 1) - s )
     lnip = Max( 0, GET_VARIABLE_INDEX(gd, 'CIP', 1) - s )
     lngl = Max( 0, GET_VARIABLE_INDEX(gd, 'CGL', 1) - s )
     lngm = Max( 0, GET_VARIABLE_INDEX(gd, 'CGM', 1) - s )
     lngh = Max( 0, GET_VARIABLE_INDEX(gd, 'CGH', 1) - s )
     lnf  = Max( 0, GET_VARIABLE_INDEX(gd, 'CFW', 1) - s )
     lnhl = Max( 0, GET_VARIABLE_INDEX(gd, 'CHL', 1) - s )
     lnhf = Max( 0, GET_VARIABLE_INDEX(gd, 'CHF', 1) - s )
     lnhlf= Max( 0, GET_VARIABLE_INDEX(gd, 'CHLF', 1) - s )

     lnchaff = Max( 0, GET_VARIABLE_INDEX(gd, 'CCHAFF', 1) - s )
     lnox    = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX', 1) - s )
     lnox_a  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_A', 1) - s )
     lnox_b  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_B', 1) - s )
     lnox_c  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_C', 1) - s )
     lnox_d  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_D', 1) - s )
     lnox_e  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_E', 1) - s )
     lnox_f  = Max( 0, GET_VARIABLE_INDEX(gd, 'CLNOX_F', 1) - s )
     lanox    = Max( 0, GET_VARIABLE_INDEX(gd, 'CANOX', 1) - s )
     lco      = Max( 0, GET_VARIABLE_INDEX(gd, 'CCO', 1) - s )
     
     lss   = Max( 0, GET_VARIABLE_INDEX(gd, 'SSMX', 1) - s )

! Space charge variables
     lscw = Max( 0, GET_VARIABLE_INDEX(gd, 'SCCW', 1) - s )
     lscr = Max( 0, GET_VARIABLE_INDEX(gd, 'SCRW', 1) - s )
     lscm = Max( 0, GET_VARIABLE_INDEX(gd, 'SCMW', 1) - s )
     lsci = Max( 0, GET_VARIABLE_INDEX(gd, 'SCCI', 1) - s )
     lscis= Max( 0, GET_VARIABLE_INDEX(gd, 'SCIS', 1) - s )
     lscs = Max( 0, GET_VARIABLE_INDEX(gd, 'SCSW', 1) - s )
     lsch = Max( 0, GET_VARIABLE_INDEX(gd, 'SCHW', 1) - s )

     lscir = Max( 0, GET_VARIABLE_INDEX(gd, 'SCIR', 1) - s )
     lscip = Max( 0, GET_VARIABLE_INDEX(gd, 'SCIP', 1) - s )
     lscgl = Max( 0, GET_VARIABLE_INDEX(gd, 'SCGL', 1) - s )
     lscgm = Max( 0, GET_VARIABLE_INDEX(gd, 'SCGM', 1) - s )
     lscgh = Max( 0, GET_VARIABLE_INDEX(gd, 'SCGH', 1) - s )
     lscf  = Max( 0, GET_VARIABLE_INDEX(gd, 'SCFW', 1) - s )
     lschl = Max( 0, GET_VARIABLE_INDEX(gd, 'SCHL', 1) - s )

     lscpi = Max( 0, GET_VARIABLE_INDEX(gd, 'CPION', 1) - s )
     lscni = Max( 0, GET_VARIABLE_INDEX(gd, 'CNION', 1) - s )

     lscpli = Max( 0, GET_VARIABLE_INDEX(gd, 'CPLION', 1) - s )
     lscnli = Max( 0, GET_VARIABLE_INDEX(gd, 'CNLION', 1) - s )

! Particle volume:
     lvi  = Max( 0, GET_VARIABLE_INDEX(gd, 'VCI', 1) - s )
     lvs  = Max( 0, GET_VARIABLE_INDEX(gd, 'VSW', 1) - s )
     lvgl = Max( 0, GET_VARIABLE_INDEX(gd, 'VGL', 1) - s )
     lvgm = Max( 0, GET_VARIABLE_INDEX(gd, 'VGM', 1) - s )
     lvgh = Max( 0, GET_VARIABLE_INDEX(gd, 'VGH', 1) - s )
     lvf  = Max( 0, GET_VARIABLE_INDEX(gd, 'VFW', 1) - s )
     lvh  = Max( 0, GET_VARIABLE_INDEX(gd, 'VHW', 1) - s )
     lvhl = Max( 0, GET_VARIABLE_INDEX(gd, 'VHL', 1) - s )

! reflectivity (6th moment)

     lzr  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZRW', 1) - s )
     lzm  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZMW', 1) - s )
     lzi  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZCI', 1) - s )
     lzs  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZSW', 1) - s )
     lzgl = Max( 0, GET_VARIABLE_INDEX(gd, 'ZGL', 1) - s )
     lzgm = Max( 0, GET_VARIABLE_INDEX(gd, 'ZGM', 1) - s )
     lzgh = Max( 0, GET_VARIABLE_INDEX(gd, 'ZGH', 1) - s )
     lzf  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZFW', 1) - s )
     lzh  = Max( 0, GET_VARIABLE_INDEX(gd, 'ZHW', 1) - s )
     lzhl = Max( 0, GET_VARIABLE_INDEX(gd, 'ZHL', 1) - s )

     
     IF (  microphys(1:4) .eq. 'ZVDH' .or. microphys(1:5) .eq. 'ZIEGH'  .or.  microphys(1:5) .eq. 'ZVDMH' ) THEN
       lhab = lhl
       lg   = lh
       lschab = lschl
     ELSEIF ( microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) .eq. 'ZVD' ) THEN
       lhab = Max(lh,lf)
       lg   = lh
       lschab = lsch
     ELSEIF ( microphys(1:3) .eq. 'LFO' ) THEN
       lhab = lh
       lg   = lh
       lschab = lsch
     ELSEIF ( microphys(1:5) .eq. 'ICE10' ) THEN
       lhab = lhl
       lg   = lgl
       lschab = lschl
     ELSEIF ( microphys(1:8) .eq. 'WARMZIEG' ) THEN
       lhab = lr
       lg   = lr
     ELSEIF ( microphys(1:3) .eq. 'DRY' ) THEN
       lhab = lv
     ELSE
       IF ( lhl .gt. 1 ) THEN
         lhab = lhl
         lg   = lh
         lschab = lschl
       ELSEIF ( lh .gt. 1 ) THEN
         lhab = Max(lh,lf)
         lg   = lh
         lschab = lsch
       ELSEIF ( lr .gt. 1 ) THEN
         lhab = lr
         lg   = lr
       ENDIF
     ENDIF
     
     IF ( lhab .gt. lqmx ) THEN
       write(0,*) 'STOP! lhab exceeds lqmx! Increase lqmx in index_module!!',lhab,lqmx
     ENDIF

!     lqi = lc
!     lqb  = lc
     lqe  = lhab
     lscb = lscw
     IF ( lscnli > 1 ) THEN
       lsce = lscnli
     ELSE
       lsce = lscni
     ENDIF
     lsceq= lschab
     
     IF ( microphys(1:4) .eq. 'TAKE' ) THEN
       lscb = lscr
       lsceq = lscpi-1
       lqb = lnr
       lqe = lnhl + ntakpd - 1
     ENDIF
     
!     IF ( lsceq .gt. lscmx ) THEN
!       write(0,*) 'ERROR! lsceq .gt. lscmx! Increase lscmx!!!'
!       STOP
!     ENDIF
     
!     print*, 'lt,lv,lc,lr,li,ls,lh = ',lt,lv,lc,lr,li,ls,lh
!     print*, 'lccn,lnc,lnr,lni,lns,lnh,lss = ',lccn,lnc,lnr,lni,lns,lnh,lss
!     print*, 'lscw = ',lscw

 
 END SUBROUTINE INDEX_INIT

!-------------------------------------------------------------------------------
! 
!
!         <<<<<<<<<<<<<<<  GRID_DEFINE_FROM_FILE  >>>>>>>>>>>>>>>>>>>>  
!  
!
! Written by Lou Wicker, Nov/Dec 2005
!  
!-------------------------------------------------------------------------------

!  SUBROUTINE GRID_DEFINE_FROM_FILE(gd,filename)
  
  SUBROUTINE GRID_DEFINE_FROM_FILE(gd)
  
   USE PARAM_MODULE

   implicit none
   
!   character(LEN=*) filename
   integer, parameter :: maxvars = 1000
   character(LEN=255) lines(maxvars),line
   type(GRID) :: gd

   integer nvar, nattr, n, count, fend
   integer ibeg, iend
   character(LEN=10) dummy_char
   character(LEN=type_length) type
   integer i,iadd
   
!   open(unit=51, file=filename, form='formatted')
!   rewind(51)
   i = 0
   
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file

 i = i + 1
   read(f(i),*) nattr
   IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_FILE:  NATTR = ', nattr
   allocate(gd%attr(nattr))

 i = i + 1
   read(f(i),*) nvar
!   IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_FILE:  NVAR  = ', nvar

 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file

! Read in grid attributes

101 format(1x,a,2x,a,2x,i3,2x,g15.5)

   DO n = 1,nattr

  i = i + 1
   read(f(i),101) gd%attr(n)%name, gd%attr(n)%type
    gd%attr(n)%str(1:desc_length) = ' '            ! Did this because non-initialized strings get garbage put in them...

    IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_FILE:  ATTR  = ', gd%attr(n)%name

   ENDDO
   
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file
 i = i + 1
   read(f(i),*) dummy_char 				! this is just to skip a comment line in the file

! Read in grid variables (we have omitted the "index" parameter - that is set later...)

102 format(1x,a,2x,i1,2x,a,2x,8(i2,2x),a,2x,a)

    n = 0
    fend = 0

    DO n = 1,maxvars
   i = i + 1
   read(f(i),'(a)') line
     IF ( line(1:3) .ne. 'EOF' ) THEN
!       write(6,'(i3,a)') n, line(1:50)
       lines(n) = line
     ELSE
       nvar = n - 1
       fend = 1
       EXIT
     ENDIF
    
    ENDDO
    
    IF ( fend .ne. 1 ) THEN
      write(0,*) 'DEFINE_GRID_FROM_FILE: Too many variables! Current max is ',maxvars
      STOP
    ENDIF

   IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_FILE:  NVAR  = ', nvar
   allocate(gd%var(nvar))

   DO n = 1,nvar

!  i = i + 1
!   read(f(i),102) gd%var(n)%name,    &
    read(lines(n),102)              &
                 gd%var(n)%name,    &
                 gd%var(n)%dim,     &
!                 gd%var(n)%ng,      &
                 gd%var(n)%type,    &
                 gd%var(n)%tdepend, &
                 gd%var(n)%istag,   &
                 gd%var(n)%jstag,   &
                 gd%var(n)%kstag,   &
                 gd%var(n)%pdef,    &
                 gd%var(n)%dyntype, &
                 gd%var(n)%phytype, &
                 gd%var(n)%buotype, &
                 gd%var(n)%unit,    &
                 gd%var(n)%description
                 
!                 IF ( gd%var(n)%type == 'nxyz4d' ) THEN
!                   gd%var(n)%nbins = gd%var(n)%dim
!                   gd%var(n)%dim  = 4
!                   iadd = gd%var(n)%nbins
!                 ELSE
!                   gd%var(n)%nbins = 0
                   iadd = 1
!                 ENDIF

    CALL STRING_LIMITS(gd%var(n)%name, ibeg, iend)
    
    IF ( gd%var(n)%dim .eq. 3 ) THEN ! .or. gd%var(n)%type(3:4) .eq. '2d' ) THEN
      IF ( gd%var(n)%tdepend .ge. 1 ) THEN
        gd%var(n)%field = gd%var(n)%name(ibeg:iend)//', scalar, series'
      ELSE
        gd%var(n)%field = gd%var(n)%name(ibeg:iend)//', scalar'
      ENDIF
    ELSE
      gd%var(n)%field = CHAR(0)
      gd%var(n)%positions = CHAR(0)
    ENDIF
    
    IF( n .eq. 1 ) THEN
      type = gd%var(1)%type
      count = 1 
    ENDIF

    IF( type .ne. gd%var(n)%type ) THEN
      type = gd%var(n)%type
      count = 1
    ENDIF
    
    IF ( gd%var(n)%dim .eq. 3 .or. gd%var(n)%dim .eq. 4 ) THEN
      gd%var(n)%basedim = 1
      
    ELSE
      gd%var(n)%basedim = 0
    ENDIF

    if (gd%var(n)%dim .eq. 0) then
     gd%var(n)%ng = 0 
    else
     gd%var(n)%ng = ng
    endif
    
    IF ( gd%var(n)%dim .ge. 3 ) THEN

     gd%var(n)%positions = ' '

#ifdef OPENDX
     IF ( gd%var(n)%dim .eq. 3 .or. gd%var(n)%type .eq. 'yz2d' .or. gd%var(n)%type .eq. 'xz2d' ) THEN
      CALL STRING_LIMITS(gd%var(n)%positions, ibeg, iend)
      
      IF ( gd%var(n)%kstag .eq. 0 ) THEN
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'ZCDX, product;'
      ELSE
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'ZEDX, product;'
      ENDIF
    ENDIF

     IF ( gd%var(n)%dim .eq. 3 .or. gd%var(n)%type .eq. 'xy2d' .or. gd%var(n)%type .eq. 'yz2d' ) THEN
      CALL STRING_LIMITS(gd%var(n)%positions, ibeg, iend)
      
      IF ( gd%var(n)%jstag .eq. 0 ) THEN
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'YCDX, product;'
      ELSE
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'YEDX, product;'
      ENDIF
     ENDIF

     IF ( gd%var(n)%dim .eq. 3 .or. gd%var(n)%type .eq. 'xy2d' .or. gd%var(n)%type .eq. 'xz2d' ) THEN
     
      CALL STRING_LIMITS(gd%var(n)%positions, ibeg, iend)
     IF ( gd%var(n)%istag .eq. 0 ) THEN
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'XCDX, product;'
      ELSE
        gd%var(n)%positions = gd%var(n)%positions(ibeg:iend)//'XEDX, product;'
      ENDIF
     ENDIF
#endif

    ENDIF

    gd%var(n)%index = count
    count = count + iadd

    IF( DEBUG .or. DEBUG_IO ) write(6,*) 'DEFINE_GRID_FROM_FILE:  VAR  = ', gd%var(n)%name, &
                                                                            gd%var(n)%type, &
                                                                            gd%var(n)%index
   ENDDO


  END SUBROUTINE GRID_DEFINE_FROM_FILE

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GRID_SET_ATTRIBUTES(gd, prefix, ne, member, coards  &
                                ,nx, ny, nz                      &
                                ,microphys, v5dflds              &
#ifdef MPI                                
                                ,nxe,nye,nze,tileindex)
#else 
                                )
#endif

   USE PARAM_MODULE

   implicit none 
   
   type(GRID)         :: gd
   character(LEN = *) :: prefix
   integer            :: ne
   integer            :: member
   integer            :: nx, ny, nz
   integer            :: coards(6)
   character(LEN = *) :: microphys
   character(LEN = *) :: v5dflds
#ifdef MPI
   integer            :: nxe, nye, nze
   integer            :: tileindex
#else
   integer            :: tileindex = 1
#endif
   integer n
   logical return
   integer ibeg, iend
   character(LEN=8) date
   character(LEN=10) time
   character(LEN=5) zone
   character(LEN=40) dummy
   character(LEN=8)  num
   character(LEN=4)  num1,num2
   character         yr*4, mth*2, day*2, hr*2, min*2, sec*2 
   character(LEN=100) tmpname

   IF( .not. SET_ATTRIBUTE(gd, 'NX', nx) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NX ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NY', ny) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NY ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NZ', nz) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NZ ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NG', ng) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NG ATTRIBUTE'

#ifdef MPI
   IF( .not. SET_ATTRIBUTE(gd, 'NXEND', nxe) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NXEND ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NYEND', nye) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NYEND ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NZEND', nze) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NZEND ATTRIBUTE'
#else
   IF( .not. SET_ATTRIBUTE(gd, 'NXEND', nx) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NXEND ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NYEND', ny) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NYEND ATTRIBUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'NZEND', nz) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING NZEND ATTRIBUTE'
#endif

   IF( .not. SET_ATTRIBUTE(gd, 'SIZE_OF_ENSEMBLE', ne) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING SIZE_OF_ENSEMBLE'
   IF( .not. SET_ATTRIBUTE(gd, 'MEMBER',       member) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING MEMBER NUMBER'
   IF( .not. SET_ATTRIBUTE(gd, 'TILE_INDEX',tileindex) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING TILE INDEX'

   CALL DATE_AND_TIME(date,time,zone)
   dummy(1:10) = date(1:4) // '-' //date(5:6) // '-' //date(7:8)

   IF( .not. SET_ATTRIBUTE(gd, 'CREATION_DATE', dummy(1:10)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING CREATION DATE'

   dummy(1:14) = time(1:4) // ' ' // zone(1:5) // ' GMT'
   IF( .not. SET_ATTRIBUTE(gd, 'CREATION_TIME', dummy(1:14)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING CREATION TIME'

! Parse the coards time information

   write(yr, '(i4.4)') coards(1)
   write(mth,'(i2.2)') coards(2)
   write(day,'(i2.2)') coards(3)
   write(hr, '(i2.2)') coards(4)
   write(min,'(i2.2)') coards(5)
   write(sec,'(i2.2)') coards(6)

   dummy(1:33) = 'seconds since ' // yr // '-' // mth // '-' // day // ' ' // hr // ':' // min // ':' // sec 

   IF( .not. SET_ATTRIBUTE(gd, 'COARDS', dummy(1:33)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING COARDS'

   IF( .not. SET_ATTRIBUTE(gd, 'Conventions', 'CF-1.7') ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING Conventions'

   CALL STRING_LIMITS(prefix, ibeg, iend)

   write(num1, '(a,i3.3)') '.', member
#ifdef MPI
   write(num2, '(a,i3.3)') '.', tileindex

   num = num1//num2
#else 
   num = num1
#endif

   
!   IF( .not. SET_ATTRIBUTE(gd,'PREFIX_NAME',prefix(ibeg:iend)) )          write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING RUNNAME'
   IF( .not. SET_ATTRIBUTE(gd,'PREFIX_NAME',prefix) )  write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING RUNNAME'

    tmpname = prefix(ibeg:iend)//trim(num)//'.nc'
    
   IF( .not. SET_ATTRIBUTE(gd,'MEMBER_NAME',tmpname) )       &
              write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING MEMBERNAME'
    tmpname = prefix(ibeg:iend)//trim(num)//'.out'

   IF( .not. SET_ATTRIBUTE(gd,'OUTPUT_FILE_NAME',tmpname) ) &
              write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING OUTPUTFILE'
   CALL STRING_LIMITS(microphys, ibeg, iend)

   IF( .not. SET_ATTRIBUTE(gd,'MICROPHYS',microphys(ibeg:iend)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING MICROPHYS'

!   CALL STRING_LIMITS(v5dflds, ibeg, iend)

   IF( .not. SET_ATTRIBUTE(gd,'V5DFIELDS',v5dflds) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING VIS5D FIELD NAMES'
!   IF( .not. SET_ATTRIBUTE(gd,'V5DFIELDS',v5dflds(ibeg:iend)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING VIS5D FIELD NAMES'

! Set the time of the model run

   IF( .not. SET_ATTRIBUTE(gd, 'YEAR',   coards(1)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING YEAR'
   IF( .not. SET_ATTRIBUTE(gd, 'MONTH',  coards(2)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING MONTH'
   IF( .not. SET_ATTRIBUTE(gd, 'DAY',    coards(3)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING DAY'
   IF( .not. SET_ATTRIBUTE(gd, 'HOUR',   coards(4)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING HOUR'
   IF( .not. SET_ATTRIBUTE(gd, 'MINUTE', coards(5)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING MINUTE'
   IF( .not. SET_ATTRIBUTE(gd, 'SECOND', coards(6)) ) write(6,*) 'GRID_SET_ATTRIBUTES:  PROBLEM SETTING SECOND'

  END SUBROUTINE GRID_SET_ATTRIBUTES
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GRID_ALLOCATE(gd, nx, ny, nz)

   USE PARAM_MODULE

   implicit none 

   type(GRID), target        :: gd
   integer            :: nx, ny, nz

   integer            :: n, nxyz3d, ngs
   
   real, parameter :: spval = 0.0  ! 1.e32

! count the number of 3d arrays
   
   ngs = ngv
   
   nxyz3d = 0

   DO n = 1,size(gd%var)
    SELECT CASE( gd%var(n)%dim )
    CASE DEFAULT                                         ! Do nothing
    CASE ( 3 )        
     IF( gd%var(n)%type(1:4) .eq. 'xyz3' ) THEN 
       nxyz3d = nxyz3d + 1
     ENDIF
    CASE ( 4 )        
     IF( gd%var(n)%type .eq. 'xyz4d' ) THEN 
       IF ( gd%var(n)%name == 'QHCM' ) THEN
!       nxyz3d = nxyz3d + gd%var(n)%nbins
       ELSEIF ( gd%var(n)%name == 'VHCM' ) THEN
!       nxyz3d = nxyz3d + gd%var(n)%nbins
       ELSEIF ( gd%var(n)%name == 'QGTAK' ) THEN
       ELSEIF ( gd%var(n)%name == 'QHTAK' ) THEN
       ELSEIF ( gd%var(n)%name == 'QRTAK' ) THEN
       ELSEIF ( gd%var(n)%name == 'QITAK' ) THEN
       ENDIF
     ENDIF
    END SELECT
   ENDDO

    
    IF ( nxyz3d .ge. 1 ) THEN
     allocate( gd%xyz3d%flt4d(-ngs+1:nx+ngs,-ngs+1:ny+ngs,-ngs+1:nz+ngs,nxyz3d) )
    ENDIF


! Allocate memory for arrays

   nxyz3d = 0
  
   DO n = 1,size(gd%var)
    
    nullify ( gd%var(n)%flt3d ) ! explicitly disassociate flt3d
    
    SELECT CASE( gd%var(n)%dim )
    
    CASE DEFAULT                                         ! Do nothing
    
    CASE ( 1 )        
     IF( gd%var(n)%type .eq. 'x1d' ) allocate(gd%var(n)%flt1d(-ngs+1:nx+ngs))
     IF( gd%var(n)%type .eq. 'y1d' ) allocate(gd%var(n)%flt1d(-ngs+1:ny+ngs))
     IF( gd%var(n)%type .eq. 'z1d' ) allocate(gd%var(n)%flt1d(-ngs+1:nz+ngs))
     gd%var(n)%flt1d(:) = 0.0
     gd%var(n)%ng = ngs
    
    CASE ( 2 )        
     IF( gd%var(n)%type .eq. 'xy2d' ) allocate(gd%var(n)%flt2d(-ngs+1:nx+ngs,-ngs+1:ny+ngs))
     IF( gd%var(n)%type .eq. 'xz2d' ) allocate(gd%var(n)%flt2d(-ngs+1:nx+ngs,-ngs+1:nz+ngs))
     IF( gd%var(n)%type .eq. 'yz2d' ) allocate(gd%var(n)%flt2d(-ngs+1:ny+ngs,-ngs+1:nz+ngs))
     IF( gd%var(n)%type .eq. 'x2d' ) allocate(gd%var(n)%flt2d(3,-ngs+1:nx+ngs))
     IF( gd%var(n)%type .eq. 'y2d' ) allocate(gd%var(n)%flt2d(3,-ngs+1:ny+ngs))
     IF( gd%var(n)%type .eq. 'z2d' ) allocate(gd%var(n)%flt2d(3,-ngs+1:nz+ngs))
     gd%var(n)%flt2d(:,:) = 0.0
     gd%var(n)%ng = ngs
    
    CASE ( 3 )        
     IF( gd%var(n)%type(1:4) .eq. 'xyz3' ) THEN 
       nxyz3d = nxyz3d + 1
!       allocate( gd%var(n)%flt3d(-ngs+1:nx+ngs,-ngs+1:ny+ngs,-ngs+1:nz+ngs) )
!       gd%var(n)%flt3d(-ngs+1,-ngs+1,-ngs+1) => gd%xyz3d%flt4d(-ngs+1:nx+ngs,-ngs+1:ny+ngs,-ngs+1:nz+ngs,nxyz3d)
!       gd%var(n)%flt3d => gd%xyz3d%flt4d(-ngs+1:nx+ngs,-ngs+1:ny+ngs,-ngs+1:nz+ngs,nxyz3d)
!
! is there a better way to assign a pointer to a section of a larger array
! and set the lower bounds to be something besides 1,1,1 (i.e., not become aliased)?  This trick uses
! a subroutine to set the pointer
!       CALL GRID_POINT(gd%xyz3d%flt4d(-ngs+1:nx+ngs,-ngs+1:ny+ngs,-ngs+1:nz+ngs,nxyz3d), nx, ny, nz, ngs, gd%var(n)%flt3d)
       CALL GRID_POINT(gd%xyz3d%flt4d(-ngs+1,-ngs+1,-ngs+1,nxyz3d), nx, ny, nz, ngs, gd%var(n)%flt3d)

       allocate( gd%var(n)%base1d(-ngs+1:nz+ngs) )

     ENDIF
     gd%var(n)%flt3d(:,:,:) = 0.0

     gd%var(n)%flt3d(:,:,nz:nz+ngs) = spval
      IF ( gd%var(n)%name == 'W' ) gd%var(n)%flt3d(:,:,nz) = 0.     
     gd%var(n)%flt3d(:,:,-ngs+1:0) = spval

!     IF ( gd%var(n)%name == 'DBZ' ) gd%var(n)%flt3d(:,:,:) = 0./0.
     gd%var(n)%base1d(:) = 0.0
     gd%var(n)%ng = ngs
    
    END SELECT
   
   ENDDO

  END SUBROUTINE GRID_ALLOCATE
  
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GRID_POINT(a, nx, ny, nz, ng, p)
  
  implicit none
  
  integer :: nx,ny,nz,ng
  real, target :: a(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
  real, pointer :: p(:,:,:)
  
 
  p => a
  

  END SUBROUTINE GRID_POINT
  
