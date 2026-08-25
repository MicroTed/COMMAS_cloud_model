!
!  parameters for ncar graphics plots
!
MODULE PLOT_MODULE

   implicit none
   
   integer :: ifile2 = 1
   
   integer :: thistorysp = 0 ! fake thistory to use for hack to generate stat file
   
   logical :: timeheight = .false.
   
   character(len=60)  :: stitle = 'COMMAS RUN'   ! name of simulation for plot title
   
   character(len=120)  ::  inname,outname,inname2,flstname,cdum
   character(len=120)  ::  thfile  ! name of file with time-height ascii data
   
   integer, parameter :: imsz = 0
   
   real               :: dzplot = 500.
   
   logical            :: gridrelative = .true.    ! whether to adjust U and V by ugrid,vgrid
   logical            :: precipgroundrelative = .true.  ! whether to use the ground-relative "2" precip accums
   
   logical            :: newmotion = .false.      ! options to change relative motion to umvnew, vmvnew
   real               :: umvnew = 0.0, vmvnew = 0.0
   real               :: umove  = 0.0,  vmove  = 0.0
   
!   real               :: dtp
   
   integer            :: ijump = 0
   
   integer            :: itemplt = 1 ! flag to plot temperature on right; might not be used, use itlab instead
   
   integer            :: ithetae = 0 ! flag to compute theta-e, 1 = yes
   
   integer            :: icolor = 0  ! for preset color contour lines; may not work anymore
   
   integer            :: icdash = 2  ! set =1 for dash pattern, 
                                     ! =2 for solid with color inegclr (easier to modify in Illustrator)
                                     ! =3 for dash with same color as positive contours
   integer            :: idash  = 21845
   integer            :: inegclr = 58 ! color for negative value contours
   
   integer            :: ipmore = 1 ! =1 makes plots with pre-chosen information (labels, tickmarks, etc.).  Setting 
                                    ! this value to 0 makes plots with minimal information (quicker plots).  
   
   integer            :: ixlab = 1, iylab = 1
   
   integer            :: itlab = 2  ! whether to plot temperature ticks at 0, -10, -20, ets.
   
   integer            :: itckni = 5, itckji = 10
   
   integer            :: itckno = 2, itckjo = 6
   
   integer            :: nhic = 0   ! old conrec param
   
   integer            :: iclu = 3 
   
   integer            :: lets = 1
   
   real               :: cwm  = 1.0
   
   integer            :: ipqual = 2   ! not used
    
   real               :: gksdot = 2.0 ! not used
   
   integer            :: isforz = 1
   real               :: sforz  = 2.3
   
   real               :: xfacb = 0.001, yfacb = 0.001, zfacb = 0.001
   real               :: xdivb = 10000, ydivb = 10000, zdivb = 2000
   
   integer            :: iotyp = 0 ! 0 = ncgm; 20 = color ps; 21 = color eps (1 frame only)
   
   integer            :: isum  = 0

   integer            :: ixave  = 0, jyave = 0
   
   integer            :: ixmdiv = 0, jymdiv = 0, kzmdiv = 4
   
   integer            :: radarmask = 0  ! flag to plot wind vectors only where dbz exceeds threshold of lowdbz
   
   real               :: lowdbz = 5.   ! cutoff for wind vectors when radarmask = 1
   
   logical            :: statout = .false.

   logical            :: lchgratevolume = .false.
   real               :: chgratevolume = 1.e-12
   
   real               :: labelsize = 0.0 ! used to set LLS if .ne. 0.0
   
   integer            :: iaddxy = 0  ! whether to add origin offset
   
   integer            :: iplotcint = 1  ! whether to print 'cntr interval' info line
   integer            :: iplottitle = 1  ! whether to print title info line
   integer            :: iplottime = 1  ! how to print time info line
   integer            :: iplotidx = 1  ! whether to print field info line
   integer            :: iplotcolorbar = 1
   integer            :: iplotrefvector = 1 ! whether to print the reference vector
   
   real               :: vectoffsetx = 0.0 ! X offset position for reference vector
   real               :: vectoffsety = -0.2 ! Y offset position for reference vector
   real               :: vectcharsize = 0.015 ! vector text block character size
   integer            :: itimeoffset = 0 ! integer time to add/subtract from model time (for forecast)
   
   integer            :: ilaborhv = -1 ! -1 for default; 0 for horiz; 1 for vertical
   integer            :: ihack = 1, jhack = 1, khack = 1
   integer            :: ibase_thproc = 30 ! old is 27
   integer            :: ibase_micth = 24 ! old is 24
   logical            :: lraintyp2dswap = .false. ! for swapping 2D FD src arrays (shed and melt)

   real               :: thlengthscale = -1.0 ! normalization length in meters for time-height plots to convert from 
                                              ! 'per model level' to 'per vertical length scale'

   integer  :: recalcdbz = 0  ! =1 sets iusewetgraupel=1 and calls radardbz
   
   integer :: dbz_color_option = 1 ! 1 = old default, 2 = Ziegler/Fierro
   
   ! Dual-pol options
   logical :: do_dualpol = .false.
!   character(len=100) :: rsafndir = '/Users/ted.mansell/develop/pyCAPS-PRS_1.0/S-band' ! directory with scattering tables
   integer :: MPflg = 2
   real :: wavelen = 107.0             ! Radar wavelength (mm) e.g.: 107.0 = S-band, 32.0 = X-band
   real :: beamwidth = 1.0             ! Radar beamwidth (deg)
   character(len=256) :: dirscatt = '/Users/ted.mansell/develop/pyCAPS-PRS_1.0/S-band/MFflg1/' ! Directory containing scattering amplitude
   integer ::MFflg = 1                   ! 0 - Based on Jung et al. (2008) diagnostic water fraction. Some modifications were not thoroughly inspected.
                            ! 1 - NEW (06/27/2012) diagnostic method by DTD
                            ! 2 - Melted fraction from microphysics (If the microphysics doesn't predict water fraction on
                            ! ice, then this option effectively treats all ice as dry in the polarimetric calculations.  Useful
                            ! for comparing with options 0 or 1)
                            ! 3 - Linear melting fraction based on air temperature
                            ! 4 - Melting fraction based on air temperature. But dry and wet scattering amplitutes are averaged.
   logical :: addtorain = .True.  ! Add excess water on graupel and hail to rain? This has no effect for MFflg = 0,2
   real :: wgfac_set = 1.0            ! Reduction factor for water shell during wet growth.  Set to 1 (original) for no reduction
   logical :: dualpol_fulldomain = .false.

   integer ::    &
       ldbz = 1, &
       ldbzh = 2, &
       ldbzv = 3, &
       lzdr = 4,  &
       ldbzr = 5, &
       lzhv  = 6, &
       lkdp  = 7, &
       lahh  = 8, &
       lavv  = 9, &
       lrhv  = 10
       
   integer :: ndpvars = 10

  REAL, allocatable :: qr3d(:,:,:),qs3d(:,:,:),qg3d(:,:,:),qh3d(:,:,:)
  REAL, allocatable :: Ntr(:,:,:),Nts(:,:,:),Ntg(:,:,:),Nth(:,:,:)
  REAL, allocatable :: alphar3d(:,:,:),alphas3d(:,:,:),alphag3d(:,:,:),alphah3d(:,:,:)
  REAL, allocatable :: qsw(:,:,:),qgw(:,:,:),qhw(:,:,:)
  REAL, allocatable :: logZ(:,:,:),sumZh(:,:,:),sumZv(:,:,:),logZdr(:,:,:)
  REAL, allocatable :: logZdrrain(:,:,:)
  REAL, allocatable :: sumZhv(:,:,:),Kdp(:,:,:),Ahh(:,:,:),Avv(:,:,:),rhv(:,:,:)
  REAL, allocatable :: rhoms(:,:,:),rhomg(:,:,:),rhomh(:,:,:),rhograup(:,:,:),rhohail(:,:,:)
  REAL, allocatable :: tair(:,:,:),rhoair(:,:,:)
  integer, allocatable :: dpmask(:,:,:)

   
   integer, parameter :: maxclr = 20
   integer :: customcolors(maxclr) = -1
   real    :: customlevels(maxclr) = -999.
   
   real    :: cxmin_plot = 0.0 ! to override default cxmin
   real    :: zxmin_plot = 0.0 ! to override default zxmin

   real               :: rfac1 = 1.0, rfac2 = 0.0, rfac3 = 1.0, rfac4 = 1.0
   NAMELIST /plot/                 &
                     ifile2,       &
                     thistorysp,   &
                     outname,      &
                     timeheight,   &
                     stitle,       &
                     dzplot,       &
                     gridrelative, &
                     precipgroundrelative, &
                     umvnew,       &
                     vmvnew,       &
                     ithetae,      &
                     icolor,       &
                     icdash,       &
                     idash,        &
                     inegclr,      &
                     ipmore,       &
                     ixlab,        &
                     iylab,        &
                     itlab,        &
                     itckni,       &
                     itckji,       &
                     itckno,       &
                     itckjo,       &
                     nhic,         &
                     iclu,         &
                     lets,         &
                     cwm,          &
                     ipqual,       &
                     gksdot,       &
                     isforz,       &
                     sforz,        &
                     xfacb, yfacb, zfacb, &
                     xdivb, ydivb, zdivb, &
                     iotyp,        &
                     isum,         &
                     ixave, jyave, &
                     thfile,       &
                     ixmdiv, jymdiv, kzmdiv, &
                     radarmask,    &
                     lowdbz,       &
                     statout,      &
                     labelsize,    &
                     lchgratevolume, &
                     chgratevolume, &
                     iaddxy,        &
                     iplotcint,     &
                     itimeoffset,   &
                     ilaborhv,      &
                     recalcdbz,     &
                     rfac1, rfac2,  &
                     rfac3, rfac4,  &
                     customcolors,  &
                     customlevels,  &
                     ihack,jhack,khack, &
                     dbz_color_option,  &
                     cxmin_plot,        &
                     zxmin_plot,        &
                     iplottitle,        &
                     iplottime,        &
                     iplotidx,         &
                     iplotcolorbar,    &
                     iplotrefvector,   &
                     ibase_thproc,     &
                     vectoffsetx, vectoffsety, vectcharsize, &
                     thlengthscale
                     
!  NAMELIST /dualpol/ do_dualpol, &
!                     dirscatt,   &
!                     MPflg, MFflg, wgfac_set
                     

END MODULE PLOT_MODULE

