#define PUBPRIV public

!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        MICRO_MODULE        ////////////////////
!
!
!
!===========================================================================
 MODULE MICRO_MODULE
 
 implicit none
!
! Microphysics Parameters
!
  
  character(LEN=15) :: microp

  logical, parameter :: lwsm6 = .false. ! act like wsm6 for some single moment interactions

! some constants from WSM6
  real, parameter  :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
  real, parameter  :: roqimax = 2.08e22*dimax**8
  
! Set up default values of constants
  
  real :: rho_qr = 1000., cnor = 8.0e6  ! LFO rain params
  real :: rho_qs =  100., cnos = 3.0e6  ! LFO snow params
  real :: rho_qh =  900., cnoh = 4.0e4  ! LFO/ZIEG graupel params
  real :: rho_qf =  900., cnof = 4.0e5  ! LFO/ZIEG graupel params
  real :: rho_qhl=  900., cnohl = 4.0e3 ! ZIEG hail params
  
  real :: rho_qh_max = 900.
  real :: rho_qf_max = 900.
  real :: rho_qhl_max = 900. 
  
  real :: cnohmn  = 1.e-2 ! minimum intercept for 2-moment graupel (alphah < 0.5)
  real :: cnohlmn = 1.e-2 ! minimum intercept for 2-moment hail (alphahl < 0.5)

  real :: rho_qr_i10  = 1000., cnor_i10  = 8.0e6  ! ICE10 rain params
  real :: rho_qs_i10  =  100., cnos_i10  = 8.0e6  ! ICE10 snow params
  real :: rho_qf_i10  =  800., cnof_i10  = 4.0e5  ! ICE10 graupel/hail params
  real :: rho_qgl_i10 =  300., cnogl_i10 = 4.0e5  ! ICE10 graupel/hail params
  real :: rho_qgm_i10 =  500., cnogm_i10 = 2.0e5  ! ICE10 graupel/hail params
  real :: rho_qgh_i10 =  700., cnogh_i10 = 1.0e5  ! ICE10 graupel/hail params
  real :: rho_qh_i10  =  800., cnoh_i10  = 4.0e4  ! ICE10 graupel/hail params
  real :: rho_qhl_i10 =  900., cnohl_i10 = 1.0e3  ! ICE10 graupel/hail params
  real :: hdnmn  = 170.0  ! minimum graupel density (for variable density graupel)
  real :: fdnmn  = 170.0  ! minimum frozen drop density (for variable density graupel)
  real :: hldnmn = 500.0  ! minimum hail density (for variable density hail)
  
  ! DTD: Added parameters for Milbrandt and Yau microphysics scheme
  ! the cno values are constant intercept parameters specified for the single-moment scheme
  ! and have no effect on the double- or triple-moment versions.
  ! the alphas are constant shape parameters specified for the single- and double-moment scheme
  ! and have no effect for the "2.5-moment" diagnostic-alpha or triple-moment versions.
 
  real :: ntc_my = 1.0e8  ! Fixed NC for cloud for single-moment MY 
  real :: rho_qr_my = 1000., cnor_my = 8.0e6, alphar_my = 0.0   ! MY rain params
  real :: rho_qi_my =  500.,                  alphai_my = 0.0   ! MY ice params
  real :: rho_qs_my =  100., cnos_my = 1.0e7, alphas_my = 0.0   ! MY snow params
  real :: rho_qg_my =  400., cnog_my = 4.0e5, alphag_my = 0.0   ! MY graupel params
  real :: rho_qh_my =  900., cnoh_my = 1.0e5, alphah_my = 0.0   ! MY hail params
  
  real :: rho_qh_tak  =  300.  ! TAKAHASHI graupel density
  real :: rho_qhl_tak =  900.  ! TAKAHASHI hail/frozen drops density

! Other values for reference

! real :: rho_qr = 1000., cnor = 1.0e7  ! RAINDROPS       / MOST PARTICLES
! real :: rho_qr = 1000., cnor = 1.0e6  ! RAINDROPS       / FEWER PARTICLES THAN KESSLER
! real :: rho_qr = 1000., cnor = 1.0e5  ! RAINDROPS       / FEWEST PARTICLES PROBABLY REASONABLE

! real :: rho_qs =  450., cnos = 1.0e5  ! SNOW            / VALUE USED BY WAKIMOTO & FOVELL

! real :: rho_qh =  900., cnoh = 4.0e3  ! VERY LARGE HAIL / FEW PARTICLES
! real :: rho_qh =  900., cnoh = 4.0e8  ! LARGE HAIL      / LFO VALUES
! real :: rho_qh =  400., cnoh = 4.0e6  ! GRAUPEL & HAIL  / MORE PARTICLES
! real :: rho_qh =  400., cnoh = 4.0e8  ! GRAUPEL         / THE MOST PARTICLES
  
! Autoconversion parameters
      
  integer :: autoconversion = 1         ! 0 = classic, 1 = Berry autoconversion
                                        ! this is for LFO; for 10ICE use ircnw
  real    :: qcmincwrn      = 2.0e-3    ! qc threshold for autonconversion (LFO; for 10ICE use qminrncw for ircnw != 5)
  real    :: cwdiap         = 20.0e-6   ! threshold diameter of cloud drops (Ferrier 1994 autoconversion)
  real    :: cwdisp         = 0.15      ! assume droplet dispersion parameter (can be 0.3 for maritime)
  real    :: ccn            = 0.6e+09   ! Central plains CCN value
  real    :: ccnuf          = 0          ! ultra fine CCN
  real    :: ccnac          = 0          ! Accum. mode CCN
  real    :: ccnnu          = 0          ! Nuclei mode CCN
  real    :: ccnco          = 0          ! Coarse mode CCN
! real    :: ccn            = 1.0e+09   ! Western plains CCN value
! real    :: ccn            = 0.3e+09   ! Maritime CCN value
  integer :: iauttim        = 1         ! 10-ice rain delay flag
  real    :: auttim         = 300.      ! 10-ice rain delay time
  real    :: qcwmntim       = 1.0e-5    ! 10-ice rain delay min qc for time accrual

  real    :: cin            = 0.001e+09   ! Ice nuclei concentration (also used as background value for DeMott (2010) icenucopt=4
  real    :: naer           = 1.0e+06 ! background value for DeMott (2010) icenucopt=4
! Hole-filling parameters

  integer :: hole_fill      = 0
    


! sedimentation flags  
! itfall -> 0 = 1st order 'lagrangian box' fallout; 2 = 6th-order monotonic (not recommended)
! iscfall, infall -> fallout options for charge and number concentration, respectively
!                    1 = mass-weighted fall speed; 2 = number-weighted fallspeed.
  integer :: itfall = 0
  integer :: iscfall = 1
  integer :: irfall = -1
  integer :: isfall = -1
  integer :: iifall = 2
  integer :: infall = 4   ! 0 -> uses number-wgt for N; NO correction applied
                          ! 1 -> uses mass-weighted fallspeed for N ALWAYS
                          ! 2 -> uses number-wgt for N and mass-weighted correction for N (Method II in Mansell, 2010 JAS)
                          ! 3 -> uses number-wgt for N and Z-weighted correction for N (Method I in Mansell, 2010 JAS)
                          ! 4 -> Hybrid of 2 and 3: Uses minimum N from each method (z-wgt and m-wgt corrections) (Method I+II in Mansell, 2010 JAS)
                          ! 5 -> uses number-wgt for N and uses average of N-wgt and q-wgt instead of Max.
  integer :: infalln = 0 ! =1 to apply factor to Vn
  integer :: imydiagalpha = 0 ! apply MY diagnostic shape parameter for fall speeds (1=for fall speed only; 2=also for microphysics rates);
                              ! 3=Milbrandt (2010) for sedimentation (adjust alpha for num-wgt Vt)
  real    :: rainfallfac = 1.0 ! factor to adjust rain fall speed (single moment only)
  real    :: snowfallfac = 1.25 ! factor to adjust rain fall speed
  integer :: ifallsedonly = 0 ! 0 = apply gr/hl Vt factor always; 1 = apply Vt factor only for sedimentation; 2 = apply only for microphysics rates
  real    :: graupelfallfac = 1.0 ! factor to adjust graupel fall speed
  real    :: fdfallfac = -1.0 ! factor to adjust graupel fall speed
  real    :: hailfallfac = 1.0 ! factor to adjust hail fall speed
  real    :: icefallfac = 1.5 ! factor to adjust ice fall speed
  integer :: icefallopt = 3 ! 1= default, 2 = Ferrier ice fall speed, 3 = Ferrier adjusted for slightly higher speeds
  integer :: icdx   = 6 ! (graupel) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
  integer :: icdxhl = 6 ! (hail) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
  real    :: axh = 106.95, bxh = 0.65
  real    :: axf = 75.7149, bxf = 0.5
  real    :: axhl = 206.984, bxhl = 0.6384
  real    :: cdhmin = 0.45, cdhmax = 0.8        ! defaults for graupel (icdx=4)
  real    :: cdhdnmin = 500., cdhdnmax = 800.0  ! defaults for graupel (icdx=4)
  real    :: cdhlmin = 0.45, cdhlmax = 0.6      ! defaults for hail (icdx=4)
  real    :: cdhldnmin = 500., cdhldnmax = 800.0  ! defaults for hail (icdx=4)
  real    :: vtmaxsed = 70. ! Limit on fall speed (m/s, all moments) for sedimentation calculations. Not applied to fall speeds for microphysical rates
  
  ! DTD note 05/25/2012: I realize the following are somewhat redundant with the logic of infall above,
  ! but infall only affects the double-moment case.  The following flags are brute force flags to turn size-sorting off 
  ! for 2 or 3-moment by setting the N-weighted and Z-weighted fall speeds for rain,graupel, and hail equal
  ! to the mass-weighted.
  
  integer :: rssflg = 1   ! Rain size-sorting allowed (1, default), or disallowed (0).  If 0, sets N and Z-weighted fall speeds to q-weighted value
  integer :: sssflg = 1   ! As above but for snow
  integer :: hssflg = 1   ! As above but for graupel; hssflg = 2 allows size sorting for alphah < alphamax; -1 allows sorting for T < -0.5C
  integer :: hlssflg = 1  ! As above but for hail; hlssflg = 2 allows size sorting for alphahl < alphamax; -1 allows sorting for T < -0.5C

! 10-ice input flags

  integer :: ndebug = -1, ncdebug = 0
  integer :: iptime = 2
  integer :: ipconc = 0
  integer :: inucopt = 0
  integer :: ichaff = 0
  real    :: chaffconc = 0 ! for setting a constant chaff density everywhere (only works when ichaff /= 1)
  integer :: ilimit = 0
  
  real :: cimn = 1.0e3, cimx = 1.0e6
  
! Params for qtodbz:
  integer  :: iuseferrier = 1  ! =1: use dry graupel only from Ferrier 1994; = 0: Use Smith (wet graupel)
  integer  :: idbzci      = 1  ! =0 off; =1 treat as ice spheres; =2 parameterization from Heymsfield (1977, JAS) using qi
  integer  :: iusewetgraupel = 1 ! =1 to turn on use of QHW for graupel reflectivity (only for ZVDM -- mixedphase)
                                 ! =2 turn on for graupel density less than 300. only 
  integer  :: iusewethail = 1 ! =1 to turn on use of QHLW for hail reflectivity (only for ZVDM -- mixedphase)
                                 ! =2 turn on for graupel density less than 300. only 
  integer  :: iusewetsnow = 0 ! =1 to turn on use of QSW for snow reflectivity (only for ZVDM -- mixedphase)
                                 ! =2 turn on for snow density less than 300. only 
  integer  :: icorrecthaildbz = 1 ! =1 to adjust hail number conc. from gr->hl conversion to keep correct Z
  integer  :: icorrectfddbz = 1 ! =1 to adjust graupel/FD number conc. from rain freezing to keep correct Z
  real     :: zxmincorr = 1.e-15 ! minimum Z to run correction to C
  real     :: cxmincorr = 1.e-3 ! minimum C to run correction to C
  
  real    :: rhofrz = 900 ! density of freezing drops
  real    :: ifrzg = 1.0 ! fraction of frozen drops (Bigg freezing) going to graupel. 1=freeze all rain to graupel, 0=freeze all to hail
  real    :: ifiacrg = 1.0 ! fraction of frozen drops (3-component freezing qiacr) going to graupel. 1=freeze all rain to graupel, 0=freeze all to hail
  real    :: ifrzs = 1.0 ! fraction of small frozen drops going to snow. 1=freeze rain to snow, 0=freeze to cloud ice
  real    :: ffrzs = 0.0 ! fraction of other initiated cloud ice going to snow (vapor nuc., contact frz, HR enh.. NOT HM mult.). 1=send all to snow, 0=send to cloud ice
  real    :: f2h = 1.0 ! fraction of cloud ice conversion going to graupel (vs. frozen drops). For testing
  integer :: irwfrz = 1 ! compute total rain that can freeze (checks heat budget)
  integer :: irimtim = 0 ! future use
!  integer :: infdo = 1   ! 1 = calculate number-weighted fall speeds
  
  integer, PUBPRIV :: irimdenopt = 1 ! = 1 for default Heymsfeld and Pflaum (1985) (as used in Mansell et al. 2010); 
                                     ! = 2 for experimental Cober and List (1993)
                                     ! = 3 for experimental Macklin as used in Farley (1987, JAM)
                                     ! = 4 Saunders and Hosseini 2001
  real    :: rimc1 = 300.0, rimc2 = 0.44  ! rime density coeff. and power (Default Heymsfield and Pflaum, 1985)
  real    :: rimc3 = 170.0                ! minimum rime density
  real    :: rimc4 = 900.0                ! maximum rime density
  real    :: rimtim = 120.0               ! cut-off rime time (10ICE)
  real    :: eqtot = 1.0e-9               ! threshold for mass budget reporting
  real, PUBPRIV :: rimdenvwgt = 0.0 ! weight (0-1) given to number-weighted fall speed when calculating rime density
  
  integer :: ireadmic = 0
  real    :: alfarim = 0.0                ! 10ICE
  real    :: qrimmnc = 1.e-12, qrimmnp = 1.e-12  ! 10ICE
  
  integer :: idiagnosecnu = 0 ! =1 to diagnose cnu based on Chandrakar et al. 2016 data; =2 for Geoffroy et al. (2010, ACP)
  integer :: iccwflg = 1     ! sets max size of first droplets in parcel to 4 micron radius (in two-moment liquid)
                             ! (first nucleation is done with a KW sat. adj. step)
  integer :: issfilt = 0     ! flag to turn on filtering of supersaturation field
  integer :: icnuclimit = 0  ! limit droplet nucleation based on Konwar et al. (2012) and Chandrakar et al. (2016)
  integer :: irenuc = 5      ! options for renucleation of droplets within the cloud
                             ! =2 renucleation following Twomey/Cohard&Pinty
                             ! =5 Similar to 7 but can produce extra activated nuclei from the 'smaller' CCN at higher SS
                             ! =7 New renucleation that requires prediction of the number of activated nuclei
  integer :: irenuc3d = 0      ! =1 to include horizontal gradient in renucleation of droplets within the cloud 
  real    :: renucfrac = 0.0 ! = 0 : cnuc = cwccn
                             ! = 1 : cnuc = actual available CCN
                             ! otherwise cnuc = cwccn*(1. - renufrac) + ccnc(1:ngscnt)*renucfrac
  real    :: ssf2kmax = 10.0 ! max value for ssf**cck in irenuc=4 or 5 (was default 1.05, but usually want "unlimited" for irenuc=5
  integer :: i_uf_or_ccn = 0      ! 0 = regular UF; 1 = treat UF as regular ccn (add to qccn)
  logical :: luseccn = .true.  ! = true to use predicted CCN
                             ! = false to use constant background
  logical, parameter :: invertccn = .false. ! =true for base state of ccn=0, =false for ccn initialized in the base state
  logical :: restoreccn = .true.  ! whether or not to restore CCN when droplets evaporate and nudge CCN back to base state (qccn) (only applies if CCNA is NOT predicted)
                                  ! Refs: Takahashi (1973, JGR), Feingold et al. (1996, JGR), Lebo and Seinfeld (2011, ACP), Xue et al. (2010, JAS)
  real    :: ccntimeconst = 3600.  ! time constant for CCN restore (either for CCNA or when restoreccn = true)
  real    :: restoreccnfrac = 1.0  ! fraction of evaporated droplets that restore CCN
  real    :: cck = 0.6       ! exponent in Twomey expression 
  real    :: ciintmx = 1.0e6
 !real    :: xcradmx=40e-6 ! this is defined now in index_module
  
  real    :: cwccn, qccn ! , cwmasn,cwmasx
  real    :: ccwmx
  
  integer :: idocw = 1, idorw = 1, idoci = 1, idoir = 1, idoip = 1, idosw = 1
  integer :: idogl = 1, idogm = 1, idogh = 1, idofw = 1, idohw = 1, idohl = 1
!  integer :: ido(3:14) = / 12*1 /


! 0,2, 5.00e-10, 1, 0, 0, 0      : itype1,itype2,cimas0,icfn,ihrn,ibfc,iacr
  integer :: itype1 = 0, itype2 = 2  ! controls Hallett-Mossop process
  integer :: in_freeze_rain_first = 0 ! =1 use IN to freezed rain drops (if none, then freeze droplets)
  integer :: icenucopt = 1       ! =1 Meyers/Ferrier primary ice nucleation; =2 Thompson/Cooper, =3 Phillips (Meyers/Demott)
                                 ! = 4 DeMott et al. 2010 (uses input of 'cin' for background n_aer)
  integer :: icfn = 2                ! contact freezing: 0 = off; 1 = hack (ok for single moment); 2 = full Cotton/Meyers version
  integer :: ihrn = 0            ! Hobbs-Rangno ice multiplication (Ferrier, 1994; use in 10-ice only)
  integer :: ibfc = 1            ! Flag to use Bigg freezing on droplets (recommend default of 1 = on)
  real    :: cwfrz2snowfrac = 0.0 ! fraction of freezing droplet mass to send to snow
  real    :: cwfrz2snowratio = 5. ! Assumed number of frozen droplets in a cluster
  integer :: iremoveqwfrz = 1    ! Whether to remove (=1) or not (=0) the newly-frozen cloud droplets (ibfc=1) from the CWC used for charge separation
  integer :: iacr = 2            ! Flag for drop contact freezing with crytals 
                                 ! (0=off; 1=drops > 500micron diameter; 2 = > 300micron)
  integer :: icrcev = 1          ! 1 = old crcev; 2 = crcev scaled by vtrain ratio (num/mass); 3 = set to zero
  integer :: icracr = 1          ! Flag to turn off rain self-collection (=0 to turn off)
  integer :: icracrthresh = 1    ! For rain self-coll. thresh. use: 1 = mean diam of 2mm; 2 = rain median volume diam of 1.9mm
  integer :: ibfr = 2            ! Flag for Bigg freezing conversion of freezing drops to graupel 
                                 ! (1=min graupel size is vr1mm; 2=use min size of dfrz, 5= as for 2 and apply dbz conservation)
  integer :: ibiggsmallrain = 0  ! 1 = When rain is too small, freeze none to graupel and send all to snow (experimental)
  integer :: ibiggopt = 2        ! 1 = old Bigg; 2 = experimental Bigg (only for imurain = 1, however)
  integer :: iacrsize = 5        ! assumed min size of drops freezing by capture
                                 !  1: > 500 micron diam
                                 !  2: > 300 micron
                                 !  3: > 40 micron
                                 !  4: all sizes
                                 !  5: > 150 micron (only for imurain = 1)
  real    :: cimas0 = 6.62e-11   ! default mass (kg) of Hallett-Mossop crystals
                                 ! 6.62e-11kg results in half the diam. (60 microns) of old default value of 5.0e-10
  real    :: cimas1 = 6.88e-13   ! default mass (kg) of new ice crystals
  real    :: splintermass = 6.88e-13 ! kg
  real    :: cfnfac = 0.1        ! Hack factor that goes with icfn=1
  integer :: iscni = 4           ! default option for ice crystal aggregation/conversion to snow
  real    :: fscni = 1.0         ! factor for calculating cscni
  logical :: imeyers5 = .false.  ! .false.=off, true=on for Meyers ice nucleation for temp > -5 C
  real    :: dmincw = 15.0e-6    ! minimum droplet diameter for collection for iehw=3
  integer :: iehw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data
  integer :: iefw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data
  integer :: iehlw = 1           ! 0 -> ehlw=ehlw0; 1 -> old ehlw; 2 -> test ehlw with Mason table data
                                 ! For ehw/ehlw = 1, ehw0/ehlw0 act as maximum limit on collection efficiency (defaults are 1.0)
  integer :: ierw = 1            ! for single-moment rain (LFO/Z) 
  integer :: iehr0c = 0          ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  integer :: iehlr0c = 0         ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  real    :: eiw0 = 0.5          ! constant or max assumed ice-crystal-droplet collection efficiency
  real    :: esw0 = 0.5          ! constant or max assumed snow-droplet collection efficiency
  real    :: ehw0 = 0.9 ! 0.5    ! constant or max assumed graupel-droplet collection efficiency
  real    :: efw0 = 0.9 ! 0.5    ! constant or max assumed frozen-drop-droplet collection efficiency
  real    :: erw0 = 1.0          ! constant assumed rain-droplet collection efficiency
  real    :: ehlw0 = 0.9 ! 0.75  ! constant or max assumed hail-droplet collection efficiency
  real    :: ehr0 = 1.0          ! constant or max assumed graupel-rain collection efficiency
  real    :: efr0 = 1.0          ! constant or max assumed graupel-rain collection efficiency
  real    :: ehlr0 = 1.0         ! constant or max assumed hail-rain collection efficiency
  real   , PUBPRIV :: exwmindiam = 0.0    ! minimum diameter of droplets for riming. If set > 0, will exclude that fraction of mass/number from accretion (idea from Furtado and Field 2017 JAS but also Fierro and Mansell 2017)
  
  real    :: esilfo0 = 1.0       ! factor for LFO collection efficiency of snow for cloud ice.
  real    :: ehslfo0 = 1.0       ! factor for LFO collection efficiency of hail/graupel for snow.
  
  integer :: ircnw    = 5        ! single-moment warm-rain autoconversion option.  5= Ferrier 1994.
  real    :: qminrncw = 2.0e-3   ! qc threshold for rain autoconversion (NA for ircnw=5)
  
  integer :: iqcinit = 2         ! For ZVDxx schemes, flag to choose which way to initialize droplets
                                 ! 1 = Soong-Ogura adjustment
                                 ! 2 = Saturation adjustment to value of ssmxinit
                                 ! 3 = KW adjustment
  
  real    :: ssmxinit = 0.4      ! saturation percentage to adjust down to for initial cloud
                                 ! formation (ZVDxx scheme only)
  
  real    :: ewfac = 1.0         ! hack factor applied to graupel and hail collection eff. for droplets
  real    :: eii0 = 0.1 ,eii1 = 0.1  ! graupel-crystal sticking eff. parameters: eii0*exp(eii1*min(temcg(mgs),0.0))
                                     ! set eii1 = 0 to get a constant value of eii0
  real    :: eii0hl = 0.2 ,eii1hl = 0.0  ! hail-crystal sticking eff. parameters: eii0hl*exp(eii1hl*min(temcg(mgs),0.0))
                                     ! set eii1hl = 0 to get a constant value of eii0hl
  real, PUBPRIV :: ewi_dcmin = 15.0e-06 ! minimum droplet diameter for nonzero ewi
  real, PUBPRIV :: ewi_dimin = 30.0e-06 ! minimum ice crystal diameter for nonzero ewi
  real    :: eri0 = 0.1              ! rain efficiency to collect ice crystals
  real    :: eri_cimin = 10.e-6      ! minimum ice crystal diameter for collection by rain
  real    :: esi0 = 0.1              ! linear factor in snow-ice collection efficiency
  real    :: ehs0 = 0.1 ,ehs1 = 0.1  ! graupel/hail-snow sticking eff. parameters: ehs0*exp(ehs1*min(temcg(mgs),0.0))
                                     ! set ehs1 = 0 to get a constant value of ehs0
  integer :: iessopt = 1  ! 1 = Original (no factor); 2 = factor based on wvel; 3 = factor based on SSI (fac=0 for undersat);
                          ! 4 = as 3 but sets min factor of 0.1 and goes to full value at 0.5% SSI
  real    :: ess0 = 0.5 ,ess1 = 0.05 ! snow aggregation coefficients: ess0*exp(ess1*min(temcg(mgs),0.0))
                                     ! set ess1 = 0 to get a constant value of ess0
      ! Note: for Takahashi, must set takessopt=2 to use esstem1/esstem2
  real    :: esstem1 = -15.  ! lower temperature where snow aggregation turns on
  real    :: esstem2 = -10.  ! higher temperature for linear ramp of ess from zero at esstem1 to formula value at esstem2
  real    :: essrmax = 0.02  ! maximum snow radius (meters) for csacs
  real    :: essfrac1 = 0.5  ! snow mass fraction 1 for aggregation roll-off
  real    :: essfrac2 = 0.75 ! snow mass fraction 2 for aggregation roll-off
  integer :: iessec0flag = 0 ! flag to activate aggregation roll-off
  real    :: ehsfrac = 1.0   ! multiplier for graupel collection efficiency in wet growth
  real    :: ehimin = 0.0 ! Minimum collection efficiency (graupel/hail - ice crystal)
  real    :: ehimax = 1.0 ! Maximum collection efficiency (graupel/hail - ice crystal)
  real    :: ehsmax = 0.5 ! Maximum collection efficiency (graupel - snow)
  real    :: ecollmx = 0.5 ! Maximum collision efficiency for graup/hail with ice; used only for charging rates
  integer :: iglcnvi = 1  ! flag for riming conversion from cloud ice to rimed ice/graupel
  integer :: iglcnvs = 2  ! flag for conversion from snow to rimed ice/graupel
  
  real    :: rz          ! reflectivity conservation factor for graupel/rain
                         ! now calculated in icezvd_dr.F from alphah and rnu
                         ! currently only used for graupel melting to rain
  real    :: rzhl        ! reflectivity conservation factor for hail/rain
                         ! now calculated in icezvd_dr.F from alphahl and rnu

  real    :: rzs     ! reflectivity conservation factor for snow(imusnow=3) with rain (imurain=1)

  real    :: alphahacx = 0.0 ! assumed minimum shape parameter for zhacw and zhacr

  real    :: fconv = 1.0  ! factor to boost max graupel depletion by riming conversions in 10ICE
  
  real    :: rg0 = 400.0  ! reference graupel density for graupel fall speed
  
  integer :: rcond = 2    ! (Z only) rcond = 2 includes rain condensation in loop with droplet condensation
                                    ! 0 = no condensation on rain; 1 = bulk condensation on rain
  integer :: icond = 1    ! (Z only) icond = 1 calculates ice deposition (crystals and snow) BEFORE droplet condensation
                          ! icond = 2 does not work (intended to calc. dep in loop with droplet cond.)
  integer, PUBPRIV :: iqis0 = 2    ! = 1 for normal qis; = 2 to set qis to use T = 0C when T > 0C  
  real, PUBPRIV    :: dfrz = 0.15e-3  ! minimum diameter of frozen drops from Bigg freezing (used for vfrz) for iacr > 1
                            ! and for ciacrf for iacr=4
  real, PUBPRIV    :: dfrzs= 0.10e-3  !  diameter for snow from Bigg freezing 
  real, PUBPRIV    :: dmlt = 3.0e-3  ! maximum (mean volume) diameter for rain melting from graupel and hail
  real, PUBPRIV    :: dshd = 1.0e-3  ! nominal diameter for drops shed from graupel/hail
  integer, PUBPRIV :: ivshdgs   = 1  ! 0 = use dshd for all shedding (non-mixedphase); 1 = use vshdgs with sheddiam
  real, PUBPRIV    :: fshed2cld = 1.0 ! fraction reduction with ished2cld=2
  integer, PUBPRIV :: ished2cld = 0  ! 1: Send shed liquid (from wet growth) to cloud droplets
                                     ! 2: Reduce collection to offset shedding
                                     ! 3: As for 2 but only for T < 0
                                     ! 4: As for 2 but only for T > 0
  integer :: ihmlt = 2      ! 1=old melting with vmlt; 2=new melting using mean volume diam of graupel/hail
  integer :: imltshddmr = 2 ! 0 (default)=mean diameter of drops produced during melting+shedding as before (using mean diameter of graupel/hail
                            ! and max mean diameter of rain)
                            ! 1=new method where mean diameter of rain during melting is adjusted linearly downward 
                            ! toward 3 mm for large (> sheddiam) graupel and hail, to take into account shedding of 
                            ! smaller drops.  sheddiam0 controls the size of graupel/hail above which the assumed 
                            ! mean diameter of rain is set to 3 mm
                            ! Only valid for ihmlt = 2 for ZVD(H) but also applies to ZVD(H)M
                            ! 2 = method that sets the resulting rain size ( vshdgs ) according to the mass-weighted diameter of the ice

   real  :: mltdiam1 = 9.0e-3, mltdiam2 = 16.0e-3, mltdiam3 = 19.0e-3, mltdiam4 = 200.0e-3, mltdiam05 = 4.5e-3

  integer :: nsplinter = 0  ! number of ice splinters per freezing drop, if negative, then per resulting graupel particle
                            ! Set nslpinter >= 1000 to turn on Lawson 2015 splintering option
                            ! nslpinter = 1001 (NSSL) applies temperature-based factor from Sullivan et al. 2018
                            ! nslpinter = 1000 (TAK) applies temperature-based factor from Sullivan et al. 2018
  real    :: lawson_splinter_fac = 2.5e-11  ! constant in Lawson et al. (2015, JAS) for ice particle production from freezing drops
  integer :: isnwfrac = 0   ! 0= no snow fragmentation; 1 = turn on snow fragmentation (Schuur, 2000)

  integer :: denscale = 1  ! 1=scale num. conc. and charge by air density for advection, 0=turn off for comparison
  
  real, PUBPRIV  :: qhdpvdn = -1.
  real, PUBPRIV  :: qhacidn = -1.
  
  logical :: mixedphase = .false.   ! .false.=off, true=on to include mixed phase graupel
  logical :: qsdenmod = .false.     ! true = modify snow density by linear interpolation of snow and rain density
  logical :: qhdenmod = .false.     ! true = modify graupel density by linear interpolation of graupel and rain density
  logical :: qsvtmod = .false.      ! true = modify snow fall speed by linear interpolation of snow and rain vt
  real    :: sheddiam   = 8.0e-03  ! minimum diameter of graupel before shedding occurs
  real    :: sheddiamlg = 10.0e-03  ! diameter of hail to use fwmlarge
  real    :: sheddiam0  = 20.0e-03  ! diameter of hail at which all water is shed
  
  real    :: fwmhtmptem = -15. ! temperature at which fwmhtmp fully switches to liquid water only being on large particles
  integer :: ifwmhtmptemopt = 1 ! option to use fwmhtmptem (1) or dwet (2) for max liquid at T < 0.
  integer :: ifwmhopt = 2 ! option for calculating maximum liquid fraction when fwmh and/or fwmhl is set to -1
                          ! 1 = maximum based on size of maximum mass diameter (deprecated)
                          ! 2 = integrate over spectrum for maximum liquid 

  integer :: ihxw2rain = 0 ! = 0 no transfer
                           ! = 1 transfer completely melted (99.5%) graupel/hail to rain when fwmh/fwmhl is set to -1.
  
  real    :: fwms = 0.5 ! maximum liquid water fraction on snow
  real    :: fwmh = 0.5 ! maximum liquid water fraction on graupel (also for frozen drops species)
  real    :: fwmhl = 0.5 ! maximum liquid water fraction on hail
  real    :: fwmlarge = 0.2 ! maximum liquid water fraction on hail larger than sheddiam
  integer :: ifwmfall = 0   ! whether to interpolate toward rain fall speed for graupel and hail
                            ! when diam < sheddiam and liquid fraction is predicted (0=no, 1=yes)
  
  logical :: rescale_high_alpha = .false.  ! whether to rescale number. conc. when alpha = alphamax (3-moment only)
  logical :: rescale_low_alpha = .true.    ! whether to rescale Z (graupel/hail) when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphar = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphah = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphahl = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)

  real, parameter :: alpharmax = 8. ! limited for rwvent calculation
  
  integer ::  ihlcnh = -1  ! which graupel -> hail conversion to use
                           ! 1 = Milbrandt and Yau (2005) using Ziegler 1985 wet growth diameter
                           ! 2 = Straka and Mansell (2005) conversion using size threshold
                           ! 3 = Conversion using wet growth diameter
  real    :: hlcnhdia = 1.e-3 ! threshold diameter for graupel -> hail conversion for ihlcnh = 1 option.
  real    :: hlcnhqmin = 0.1e-3 ! minimum graupel mass content for graupel -> hail conversion (ihlcnh = 1)
  real    :: hldia1 = 10.0e-3  ! threshold diameter for graupel -> hail conversion for ihlcnh = 2 option.
  integer :: incwet = 0    ! flag to do wet growth only on D > D_wet
  integer :: iusedw = 0    ! flag to use experimental wet growth ice diameter for gr -> hl conversion (=1 turns on)
  real    :: dwmin  = 5.0e-3 ! Minimum diameter with iusedw (can stay at 0 or be set to something larger)
  real    :: dwetmin  = 2.0e-3 ! Minimum diameter for wet growth (and shedding) with incwet (can stay at 5mm or be set to something larger)
  real    :: dwmax  = 15.e-3 ! for ihlcnh, always convert this size and larger whether or not there is wet growth
  real    :: dwtempmin = 242. ! lowest temperature to allow wet growth conversion to hail
  real    :: dwehwmin = 0.   ! Minimum ehw to use to find wet growth diameter (if > ehw0, then wet growth diam becomes smaller)
  real    :: dg0thresh = 0.15 ! graupel wet growth diameter above which we say do not bother
  real    :: wetgrthtoffset = -1. ! maximum temperature (Celcius) for wet growth
  real    :: hailcnvtoffset = -2. ! maximum temperature (Celcius) for hail conversion
  integer :: ifddenfac = 0  ! = 1 to use density threshold to count FD as GR when converting to HL
  real    :: fddenthresh = 500. ! if ifddenfac > 0, then hail from FD with lower density are considered to come from graupel
  integer :: icvhl2h = 0   ! allow conversion of hail back to graupel when hail density gets close to minimum allowed
  real    :: hlcnvtimeconst = 1.0 ! time constant (seconds) for graupel -> hail conversion (based on Ferrier 1994) for ihlcnh=3

  integer :: imurain      = 1 ! 3 for gamma-volume, 1 for gamma-diameter DSD for rain.
  integer :: imusnow      = 3 ! 3 for gamma-volume, 1 for gamma-diameter DSD for snow (=1 NOT IMPLEMENTED!!).
  integer :: iturbenhance = 0 ! warm-rain collision enhancement
                              ! 1 = enhance autoconversion only
                              ! 2 = add rain collection of cloud
                              ! 3 = add rain self-collection
  integer :: isedonly     = 0 ! 1 = only do sedimentation and skip other microphysics
                              ! 2 = allow sedimentation and normal microphysics but not droplet nucleation/condensation
                              ! 3 = (same as 0: does sedimentation and all microphysics)
                              ! 4 = allow sedimentation and only run the microphysics for axtra fields; no droplet nucleation/condensation (otherwise like 1)
  integer :: iferwisventr = 2 ! =1 for Ferrier rwvent, =2 for Wisner rwvent (imurain=1)
  integer :: izwisventr   = 2 ! =1 for old Ziegler rwvent, =2 for Wisner-style rwvent (imurain=3)
  integer :: iresetmoments = 0 ! if =1 , then set all moments to zero when one of them is zero (3-moment only)
                               ! if = -1, then reset moments only if q or Z is zero, but calculate new N if q and Z are valid and N = 0
                               ! = 0, then N or Z is recalculated if either, but not both, is missing (below threshold) (and must have q > qxmin)
  integer :: imaxdiaopt    = 3 ! = 1 use mean diameter for breakup (for rain)
                               ! = 2 use maximum mass diameter for breakup (for rain)
                               ! = 3 use mass-weighted diameter for breakup (for rain)
  integer :: irainbreak    = -1 ! 1 = on (no diameter dependence) (recommend using option 2)
                                ! 2 = (recommended) as for 1, but apply factor of 1-ec0 to turn off a smaller diameter (ec0 is rain self-coll factor)
                                ! 10 = as for 1, but sets ec0=1 for rain self-collection (i.e., no passive breakup); set higher rainbreakfac for this option
                                ! 11 = breakup for DSD tail only; uses draintail etc.
  real    :: rainbreakfac   = 2.5e6 ! 2.0e6 for lower hand fit; 2.542e6 for 'best' fit
  real    :: draintail      = 10.e-3 ! starting size for rain breakup (irainbreak = 11)
  real    :: drsmall        = 1.e-3 ! size of small drops from breakup (irainbreak = 11)
  real    :: qrbrthresh1    = 0.1e-3 ! lower threshold rain content (kg/m^3) for large drop breakup (irainbreak=11)
  real    :: qrbrthresh2    = 1.0e-3 ! upper threshold rain content (kg/m^3) for large drop breakup (irainbreak=11)
  integer :: dmrauto       = 0 ! = -1 no limiter on crcnw
                               ! =  0 limit crcnw when qr > 1.2*L (Cohard-Pinty 2002)
                               ! =  1 DTD version based on MY code
                               ! =  2 DTD mass-weighted test
                               ! =  3 Milbrandt version (from Cohard and Pinty's code)
  integer :: dmropt = 0 ! extra option for crcnw
  integer :: dmhlopt = 0 ! options for graupel -> hail conversion
  integer :: irescalerainopt = 3 ! 0 = default option
                                 ! 1 = qx(mgs,lc) > qxmin(lc) 
                                 ! 2 = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
                                 ! 3 = temcg(mgs) > rescale_tempthresh .and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
  real    :: rescale_wthresh = 3.0
  real    :: rescale_tempthresh = 0.0
  real    :: alpharaut = 0.0 ! MY2005 for autoconversion
  real    :: cxmin = 1.e-8  ! threshold cutoff for number concentration
  real    :: zxmin = 1.e-28 ! threshold cutoff for reflectivity moment
  
  integer :: ithompsoncnoh = 0 ! For single moment graupel only
                               ! 0 = fixed intercept
                               ! 1 = intercept based on graupel mass

  integer :: ivhmltsoak = 1   ! 0=off, 1=on : flag to simulate soaking (graupel/hail) during melting 
                             ! when liquid fraction is not predicted
  logical :: iwetsoak = .true. ! soak and freeze during wet growth or not
  integer :: ioldlimiter = 0 ! test switch for old(=0) or new(=2) size limiter at the end of GS for 3-moment categories (1 and 3 are for testing, see code for details)
  integer :: igrplfall = 0   ! DTD: temporary test switch for choosing between graupel fall speed parameters to make them more or less like hail.  Added because I got tired
                             ! of switching between versions of icezvd_dr.F90 and recompiling
                             ! 0 = original, 1 = Ferrier hail fall speed parameters, 2 = set cdxh to 0.45 (hail-like) 
  integer :: isnowfall = 2   ! Option for choosing between snow fall speed parameters
                             ! 1 = original Zrnic et al. (Mansell et al. 2010)
                             ! 2 = Ferrier 1994 (results in slower fall speeds)
                             ! 3 = Test using Cox fall speed as a function of mean particle *mass*
  integer :: isnowdens = 1   ! Option for choosing between snow density options
                             ! 1 = constant of 100 kg m^-3
                             ! 2 = Option based on Cox 
  integer :: ibiggsnow   = 3 ! 1 = switch conversion over to snow for small frozen drops from Bigg freezing
                             ! 2 = switch conversion over to snow for small frozen drops from rain-ice interaction
                             ! 3 = switch conversion over to snow for small frozen drops from both
                             ! 0 = no conversion, all frozen drops go to graupel
  real    :: biggsnowdiam = -1.0 ! If >0, use for ibiggsnow threshold
  
  integer :: ixtaltype = 1 ! =1 column, =2 disk (similar to Takahashi)
  
  real    :: evapfac     = 1.0 ! Multiplier on rain evaporation rate (also for liquid on snow, graupel, and hail)
  real    :: depfac      = 1.0 ! Multiplier on graupel/hail sublimation/deposition rate
  real    :: meltfac     = 1.0 ! Multiplier on graupel/hail melting rate
  
  integer :: ibinhmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of graupel, and appropriate shed drop sizes 
                           ! =2 to test melting by temporary bins
  integer :: ibinhlmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of hail, and appropriate shed drop sizes 
                            ! =2 to test melting by temporary bins
  integer :: ibincracr = 0 ! =1 use bin code to compute rain self-collection and breakup
  integer :: ibinnum   = 2  ! number of bins for melting of smaller ice (for ibinhmlr = 1)
  integer, PUBPRIV :: iqhacrmlr = 1  ! turn on/off qhacrmlr (rain collection for melting)
  integer, PUBPRIV :: iqhlacrmlr = 1  ! turn on/off qhlacrmlr
  integer, PUBPRIV :: iqhacwshr = 1  ! turn on/off qhacw for T > 0
  integer, PUBPRIV :: iqhlacwshr = 1  ! turn on/off qhlacw for T > 0
  real, PUBPRIV :: binmlrmxdia = 40.e-3 ! threshold diameter (graupel/hail) to switch bin-bulk melting to use standard chmlr
  real, PUBPRIV :: binmlrzrrfac = 1.0 ! factor for reflectivity change ice that sheds while melting
  real, PUBPRIV :: snowmeltdia = 0 ! If nonzero, sets the size of rain drops from melting snow.
  real, PUBPRIV :: alphasmlr0 = 14.0 ! shape parameter for drops formed from melting/shedding snow
  
  real, PUBPRIV :: delta_alphamlr = 0.5 ! offset from alphamax at which melting does not further collapse the shape parameter
  integer :: ibinhacw = 0  ! NOT USED YET =1 to test riming by temporary bins
  real    :: chmlrmult = 0.0 ! used to blend bin and bulk values of chmlr:  chlmlr(mgs) = (chlmlr(mgs) + chmlrmult*chlmlrtmp)/(1.+ chmlrmult)
  
  integer :: iqvsopt = 1 ! =0 use old default for tabqvs; =1 use Bolton formulation (Rogers and Yau) and es/(p-es) instead of approximated es/p
  
  logical :: wrfccn = .true. ! true to pass ccn to WRF version, false to use assumed ccn (as in option 17 for WRF)
  logical :: wrfdefaults = .false. ! If true, then do not set microphysics values for WRF microphysics version in solver
                                   ! Note that wrf init routine will still try to read nssl_mp_params namelist if present
  
  integer :: imaxsupopt = 4 ! how to treat saturation adjustment in two-moment droplets
                            ! 1 = add droplets with same mean mass as current droplets
                            ! 2 = add droplets with minimum radius of 30 microns
                            ! 3 = only add 1.5*cxmin to number concentration (allow max size to apply)
                            ! 4 = add droplets with minimum radius of 20 microns
  real    :: maxsupersat = 1.9 ! maximum supersaturation ratio, above which a saturation adustment is done
  real    :: maxlowtempss = 1.08 ! Sat. ratio threshold for allowing droplet nucleation at T < tfrh
  real    :: ssmxuf = 4.0 ! supersaturation at which to start using "ultrafine" CCN (if ccnuf > 0.)
  
!  Not using ihailrain. Use iehr0c and iehlr0c to turn on graupel/hail collection of rain at T > tfr (sheds all)
  integer :: ihailrain = 0 ! NOT USED =0 default turns off graupel/hail rain collection for T > tfr
                           ! =1  turns on graupel/hail rain collection for T < tfr (old behavior)
  
! put ipelec here for now....
  integer :: ipelec = 0
  logical :: idoniconly = .false.
  integer :: jchgs = 3  ! number of points near boundary where charging is turned off (to keep lightning from getting wonky)
  integer :: jchgn = 2
  integer :: ichge = 3
  integer :: ichgw = 2
  real    :: charging_border = 4000. ! width of no-charging zone from boundary
  

! stuff for HCM

  integer :: iotype = 1
  integer :: nttsav = 0
  integer :: iosave,is1b,is1e,is2b,is2e

! for Takahashi

  integer :: takccntype = 4 ! 1 = maritime, 2 = continental, 3 = single mass with index of taknucsize
  integer :: taknucsize = 3
  integer :: takcoagopt = 1 ! 1 = Berry and Reinhart (original); 2 = Bott 1998 (for liquid drops)
  integer :: takaggopt = 1 ! 1 = Berry and Reinhart (original); 2 = Bott 1998 (for ice crystals)
  integer :: takbottopt = 2 ! drops: 1 = Bott 1998 LFM; 2 = Bott 2000 EFM
  integer :: takbottoptice = 2 ! ice crystals: 1 = Bott 1998 LFM; 2 = Bott 2000 EFM
  double precision :: takgminflux = 1.d-10 ! flux limit in Bott scheme
  logical :: takbottrateonly = .true. ! true for independent calc., false for dependent (original) rates
  integer :: takshedsmall = 1 ! 1 = shed smaller drops from large melting hail; 0 = melt at largest drop size
  real    :: takshedsize1 = 0.15 ! diameter (cm) of drop shed from ice with D > 1.9 cm
  real    :: takshedsize2 = 0.3 ! diameter (cm) of drop shed from ice with D < 1.9 cm and D > 0.8 cm
  real    :: takshedsize3 = 0.45 ! diameter (cm) of drop shed from ice with D < 1.6 cm and D > 0.8 cm
  integer :: numshedregimes = 3
  integer :: takbreakup = 4 ! 0 = all breakup off; 
                            ! 1 = spontaneous breakup off (do not use!); 
                            ! 2 = spontan. breakup on only,
                            ! 3 = collisional breakup on only, 
                            ! 4 = both spontan. and collisional
  real    :: taksponbreaksize = 0.3  ! minimum radius (cm) of drops experiencing spontaneous breakup; old value was 0.05, but use 0.33 based on Young 1975 using only r>3mm
  real    :: taksponbreakmax  = 0.34 ! maximum radius (cm) to use in Srivastava spontaneous breakup exponential term
  integer :: takikcollmax = 2 ! max thickness to consider ice aggregation. 1 = only unrimed crystals; 2 = unrimed+lightly rimed
  integer :: takessopt = 2 ! =1 use full range of temperature; =2 to set ess=0 for T < -25
  real    :: takess0 = 1.0 ,takess1 = 1./7. ! snow aggregation coefficients: ess0*exp(ess1*min(temcg(mgs),0.0))
                                     ! set ess1 = 0 to get a constant value of ess0
                                     ! 1/7 = 0.143
  logical :: takrestorcn = .false. ! whether to restore CCN when small droplets/ice crystals completely evaporate.
  integer :: takbiggopt = 2 ! =1 for Wisner (original) =2 for alternative (Bigg fit to Fig. 1)
  real    :: takifrzg = 0.0 ! fraction of frozen drops going to graupel, (1-takifrzg) goes to frozen drops
  integer :: takgrmassopt = 1 ! =1 for graupel/hail bins have same radii as rain (recalc. mass)
                           ! =2 for graupel/hail bins have same MASS as rain (recalc. radii)

  logical :: takdodmwdt = .true.
  logical :: takdocrim = .true.
  logical :: takdocondsn = .true.
  logical :: takdocmelt = .true.
  logical :: takdocollsn = .true.
  logical :: takdosed = .true.
  real    :: takcondsndt  = 1.0 ! sub time step for condsn
  integer :: takrventopt  = 1 ! 1=original ventillation for liquid drops; 2=new ventillation (Beard and Pruppacher)
  integer :: takhventopt  = 1 ! 1=original ventillation for ice spheres; 2=new ventillation (Beard and Pruppacher)
  integer :: takiventopt  = 1 ! 1=original ventillation for ice spheres; 2=new ventillation (Beard and Pruppacher)
  
  integer :: takkfmax = 2 ! used to set value of kfmax
  integer :: takvf = 2  ! =0 for original graupel/hail fall speeds
                        ! =1 use icdx options (Wisner formula)
                        ! =2 Use the smaller of the options 0 and 1
  integer :: takvfw = 0  ! =0 for original droplet/rain fall speeds
                        ! =1 bulk fit to rain (only for takvf = 1)
  integer :: takviopt = 1 ! = 1 for original ice crystal fall speed; =2 Original but no Best number check; = 3 for Cox V(m)
  real    :: takrw1 = 6.35e-3 ! old value was 4.e-3  ! switch-over radius for drop fall speed (Tak 76, what radius to switch from eq. 10 to eq. 12
  real    :: takcxmin = 1.e-14
  integer :: takrhovtopt = 1 ! = 0 no adjustment to collision rates for air density on fall speeds
                             ! = 1 rhovt = Sqrt(rho00/dn(i,j,k))
                             ! = 2 rhovt = Sqrt( Sqrt(rho00/dn(i,j,k)) )
  integer :: takicenucopt = 1 ! = 0 uses original Fletcher
                              ! = 1 uses Meyers/Ferrier
                              ! = 2 uses Cooper
                              ! = 3 uses Philips-Demott
  integer :: takicenucsize = 2 ! = 1 uses original crystal size distribution
                              ! = 2 puts new crystals only into smallest size bin 1,1
  integer :: takfdiopt    = 2 ! frozen drop to ice crystal conversion
                              ! = 0 off; =1 original; =2 only matching mass
  real    :: takfdirmax   = 2.2e-3 ! maximum radius (cm) of frozen drops for conversion to ice crystals
  integer :: takicethickopt = 0 ! =1 thick plates (as Khain and Sednev), =2 hexagonal plates (thinner); =0 original Tak76 of 10microns
  real*8  :: takicethickness = 1.d-3 ! Ice crystal thickness (cm) for takicethickopt = 0 
  integer :: takdieff = 1  ! option for effective diam. used for ice fall speed; 1: dieff=di, 2: dieff calc. from sxi for vi
  integer :: takcwnucopt = 7  ! = 0 uses original
                              ! = 1 uses Twomey for number of nucleated drops
                              ! = 7 uses same as irenuc = 7
  integer :: takhallettopt = 1  ! = 0 turn off Hallett-Mossop
                           ! = 1 Hallett Mossop with Takahashi temperature function (top-hat function)
  integer :: takhallettsize = 1 ! size bin (diameter) for HM crystals (NOTE: Use takhallettmass instead)
  real    :: takhallettmass = 6.88e-10 ! units of grams
  integer :: takgrhlcv = 1 ! = 0 do not allow graupel->hail conversion 
                           ! = 1 turn on graupel-> hail conversion for graupel in wet growth
  integer :: takhlgrcv = 1 ! hail/frozen drop -> graupel conversion
  integer :: taknicfd  = 0 ! =0 off; =1 do charge separation with small frozen drops
  integer :: takehiopt = 0 ! =0 used dkfi; =1 use dkfi1 (ehi = 1)
  integer :: takfrzopt = 1 ! =1 to put smaller frozen drops into ice crystals intead of frozen drops
  real    :: takrhoi   = 0.5 ! old default was 0.9 ! ice crystal density (cgs units)
  double precision  :: taksri0   = 2.0d-3 ! radius of smallest ice crystal (cm), so 2.0e-3 cm = 20e-3 mm = 20 microns
  integer :: itakgrowth = 1 ! 1 = Kovetz and Olund; 2 = Khain et al. 2008 3-moment remapping
  integer :: itakmeltgrowth = 1 ! 1 = Kovetz and Olund; 2 = Khain et al. 2008 3-moment remapping
  integer :: idichgopt  = 1 ! 0 = dichg = di; 1 = dichg is for column with same mass
  logical :: taknewshed = .true.
  integer :: takshedopt = 1 ! 0 = no wet growth (no shedding); 1 = wet growth sheds 1mm drops; 
                            ! 2 = shed drops as in melting for larger stones; 1mm drops from D < 9mm
                            ! 3 = As for 2, but shed same-mass drops for D < 9mm
  integer :: takrimeffopt = 0 ! for graupel-liquid collection: 0 = use table values (dkwf), 1 = use ehw0 for all (dkwf1)
  
  logical :: usenucond = .true.
  logical :: usensslgs = .true.
  logical :: useoldgs  = .false.

! T.Iguchi Y2021 Update
      integer :: ac_opt = 0 ! option flag for: (1 and 2 currently for NUWRF only)
                                       ! 0 : normal NSSL CCN physics
                                       ! 1 : accumulation mode CN following Fridland et al. (2012, 2017),
                                       !     where CCN number is sum of unactivated CCN and droplet concentrations
                                       ! 2 : As for 1 but have three modes
!      logical :: ac_only = .true.  ! flag for considering ac_mode of CN only, or all nu,ac,co modes (still under construction)

      logical :: arg_para = .false.  ! flag for Abdul-Razzak_and_Ghan parameterization works similarly to flag_qndrop, and neglects irenuc, ccna(mgs), and cnuc(mgs)
      real, parameter :: nu_pmr = 7.5 * 1.e-3 * 1.e-6  ! meter; these parameter values follow Cheng et al. (2007QJ)
      real, parameter :: nu_pgw = 0.53  ! Unlike original Abdul-Razzak_and_Ghan, this value is used without log (Cheng et al. 2007QJ)
      real :: nu_kappa = 0.07  ! 0.61 activates too easily; ammonium sulfate as CCN (Petters and Kreidenweis, 2007ACP)
      real, parameter :: ac_pmr = 3.8 * 1.e-2 * 1.e-6  ! meter
      real, parameter :: ac_pgw = 0.69
      real :: ac_kappa = 0.61  ! ammonium sulfate as CCN (Petters and Kreidenweis, 2007ACP)
      real, parameter :: co_pmr = 0.51 * 1.e-6  ! meter
      real, parameter :: co_pgw = 0.77
      real :: co_kappa = 0.61  ! ammonium sulfate as CCN (Petters and Kreidenweis, 2007ACP)

      logical :: dm15_para = .false.  ! flag for DeMott et al. (2015) parameterization for heterogenous freezing, regardless of "ibfc"

      real :: ac_wthresh = 10.0 ! for W < ac_wthresh, use max of sswater and diagnosed SS; otherwise use sswater

! Morrison scheme
  integer :: morr_rimed_ice = 0 ! 0 = graupel-like; 1 = hail-like
  real    :: morr_droplet_conc = 250. ! per cm**3 constant droplet number concentration
  real    :: morr_dnr_max = 2800.e-6 ! Max. characteristic diam for Morrison scheme: lammaxr = 1/morr_dnr_max
  logical :: morr_dbzsoak = .false.
  integer :: imorrgdnglimit = 0
  real    :: morrdnglimit = 2000.E-6
! Thompson
  real    :: thom_autoconv_fac = 1.0 ! factor in denominator of autoconversion formula. For thom381 only. It did not exist in 381 and in 4.1 it is 200
!
!  gamma function lookup table
!
      integer ngm0,ngm1,ngm2
      parameter (ngm0=3001,ngm1=500,ngm2=500)
      double precision, parameter :: dgam = 0.01, dgami = 100.
      double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

      integer, parameter :: nqiacralpha =  300 ! 240 !480 ! 240 ! 120 ! number of points to discretize the shape parameter for 0 to maxalphalu
      integer, parameter :: nqiacrratio =  400 ! 100 ! 500 !50  ! 25
      real,    parameter :: maxratiolu = 100. ! 25.
      real,    parameter :: maxalphalu = 15. ! maximum shape parameter in lookup table
      real,    parameter :: minalphalu = -0.95 ! minimum shape parameter in lookup table
      real,    parameter :: dqiacralpha = maxalphalu/Float(nqiacralpha), dqiacrratio = maxratiolu/Float(nqiacrratio) 
      real,    parameter :: dqiacrratioinv = 1./dqiacrratio, dqiacralphainv = 1./dqiacralpha
      integer, parameter :: ialpstart = minalphalu*dqiacralphainv
      real :: ciacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      real :: qiacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      real :: ziacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      double precision :: gamxinflu(0:nqiacrratio,ialpstart:nqiacralpha,13,2) ! last index for graupel (1) or hail (2)

    integer, parameter :: ngdnmm = 9
    real :: mmgraupvt(ngdnmm,3)  ! Milbrandt and Morrison (2013) fall speed coefficients for graupel/hail

! for 3-moment collection coefficients
!      real dab0lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
!      real dab1lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005


    DATA mmgraupvt(:,1) / 50., 150., 250., 350., 450., 550., 650., 750., 850./
    DATA mmgraupvt(:,2) / 62.923, 94.122, 114.74, 131.21, 145.26, 157.71, 168.98, 179.36, 189.02 /
    DATA mmgraupvt(:,3) / 0.67819, 0.63789, 0.62197, 0.61240, 0.60572, 0.60066, 0.59663, 0.59330, 0.59048 /

      integer, parameter :: nqsat=1000001  ! (nqsat=20001)
      real, parameter :: fqsat=0.002,fqsati=1./fqsat
      real :: tabqvs(nqsat)=0,tabqis(nqsat)=0,dtabqvs(nqsat)=0,dtabqis(nqsat)=0

      real, parameter :: cbi = 7.66
      real, parameter :: cbw = 35.86
      real, parameter :: cai = 21.87455
      real, parameter :: caw = 17.2693882

! More formulae at https://www.eas.ualberta.ca/jdwilson/EAS372_13/Vomel_CIRES_satvpformulae.html
      real, parameter :: cbwbolton = 29.65 ! constants for Bolton formulation (Note that 273.15 - 243.5 = 29.65)
      real, parameter :: cawbolton = 17.67
      real, parameter :: esbolton = 6.112e2 ! factor of 1.e2 to convert mb to Pa


 CONTAINS

! #######################################################################
   SUBROUTINE LFO_INIT()
    
   RETURN
   END SUBROUTINE LFO_INIT

! #######################################################################

   SUBROUTINE MAKEGMOI()
   ! build lookup table for the complete gamma function (argument values from 0.01 to 30)
     implicit none
     double precision gamma_dp,arg
     integer igam
     
     gmoi(0) = 1.d32
     do igam = 1,ngm0
      arg = dgam*igam
      gmoi(igam) = gamma_dp(arg)
!      write(97,*) igam,gmoi(igam)
      end do
     
     
   END SUBROUTINE MAKEGMOI

! #######################################################################

   SUBROUTINE MAKEqiacrratio(bxh,bxhl)
     ! build lookup table to compute the number and mass fractions of rain drops 
     ! (imurain=1) greater than a given diameter. Used for qiacr and ciacr
     ! Uses incomplete gamma functions
     !
     ! The terms with bxh or bxhl will be off if the actual bxh or bxhl is different from the base value (icdx=6 option)
     !
     implicit none
      real :: bxh, bxhl ! note that this will not work correctly with icdx=6 option (Milbrandt-Morrison)
     
      integer :: i,j
      real    :: alp,ratio
      double precision  :: x,y,y2,y3,y4,y5h, y5hl,y6h, y6hl,ym1,y7
      real gamxinf
      real gamma
      double precision  gamxinfdp
      double precision  gamma_dpr
      
      DO j = ialpstart,nqiacralpha
      alp = float(j)*dqiacralpha
      y = gamma_dpr(real(1.+alp))
      y2 = gamma_dpr(real(2.+alp))
      y3 = gamma_dpr(real(3.+alp))
  !    y4 = gamma_dpr(real(4.+alp))
  !    y5h = gamma_dpr(real(5.+alp))
  !    y5hl = gamma_dpr(real(5.+alp))
  !    y6h = gamma_dpr(real(5.5+alp+0.5*bxh))
  !    y6hl = gamma_dpr(real(5.5+alp+0.5*bxhl))
      DO i = 0,nqiacrratio
        ratio = float(i)*dqiacrratio
        x = gamxinfdp( 1.+alp, ratio )
!        write(0,*) 'i, x/y = ',i, x/y
        ciacrratio(i,j) = x/y
        ! graupel (.,.,.,1)
        gamxinflu(i,j,1,1) = x/y
        gamxinflu(i,j,2,1) = gamxinfdp( 2.0+alp, ratio )/y
        gamxinflu(i,j,3,1) = gamxinfdp( 2.5+alp+0.5*bxh, ratio )/y
!        gamxinflu(i,j,5,1) = gamxinfdp( 5.0+alp, ratio )/y
!        gamxinflu(i,j,6,1) = gamxinfdp( 5.5+alp+0.5*bxh, ratio )/y
        gamxinflu(i,j,5,1) = (gamma_dpr(real(5.0+alp)) - gamxinfdp( 5.0+alp, ratio ))/y
        gamxinflu(i,j,6,1) = (gamma_dpr(real(5.5+alp+0.5*bxh)) - gamxinfdp( 5.5+alp+0.5*bxh, ratio ))/y
        gamxinflu(i,j,9,1) = gamxinfdp( 1.0+alp, ratio )/y
        gamxinflu(i,j,10,1)= gamxinfdp( 4.0+alp, ratio )/y

        gamxinflu(i,j,12,1) = gamxinfdp( 2.0+alp, ratio )/y2
        gamxinflu(i,j,13,1) = gamxinfdp( 3.0+alp, ratio )/y3
        
        ! hail (.,.,.,2)
        gamxinflu(i,j,1,2) = gamxinflu(i,j,1,1)
        gamxinflu(i,j,2,2) = gamxinflu(i,j,2,1)
        gamxinflu(i,j,3,2) = gamxinfdp( 2.5+alp+0.5*bxhl, ratio )/y
        gamxinflu(i,j,5,2) = gamxinflu(i,j,5,1)
!        gamxinflu(i,j,5,2) = (gamma_dpr(5.0+alp) - gamxinfdp( 5.0+alp, ratio ))/y
        gamxinflu(i,j,6,2) = (gamma_dpr(real(5.5+alp+0.5*bxhl)) - gamxinfdp( 5.5+alp+0.5*bxhl, ratio ))/y
        gamxinflu(i,j,9,2) = gamxinflu(i,j,9,1)
        gamxinflu(i,j,10,2)= gamxinflu(i,j,10,1)

        gamxinflu(i,j,12,2) = gamxinflu(i,j,12,1)
        gamxinflu(i,j,13,2) = gamxinflu(i,j,13,1)

      IF ( alp > 1.1 ) THEN
!       gamxinflu(i,j,7,1) = gamxinfdp( alp - 1., ratio )/y
       gamxinflu(i,j,7,1) = (gamma_dpr(real(alp - 1.)) - gamxinfdp( alp - 1., ratio ))/y 
!       gamxinflu(i,j,8,1) = gamxinfdp( alp - 0.5 + 0.5*bxh, ratio )/y
       gamxinflu(i,j,8,1) = (gamma_dpr(real(alp - 0.5 + 0.5*bxh)) - gamxinfdp( alp - 0.5 + 0.5*bxh, ratio ))/y
!       gamxinflu(i,j,8,2) = gamxinfdp( alp - 0.5 + 0.5*bxhl, ratio )/y
       gamxinflu(i,j,8,2) = (gamma_dpr(real(alp - 0.5 + 0.5*bxhl)) - gamxinfdp( alp - 0.5 + 0.5*bxhl, ratio ))/y
      ELSE
!       gamxinflu(i,j,7,1) = gamxinfdp( .1, ratio )/y
       gamxinflu(i,j,7,1) = (gamma_dpr(real(0.1)) - gamxinfdp( 0.1, ratio ) )/y
!       gamxinflu(i,j,8,1) = gamxinfdp( 1.1 - 0.5 + 0.5*bxh, ratio )/y
!       gamxinflu(i,j,8,2) = gamxinfdp( 1.1 - 0.5 + 0.5*bxhl, ratio )/y
       gamxinflu(i,j,8,1) = (gamma_dpr(real(1.1 - 0.5 + 0.5*bxh)) - gamxinfdp( 1.1 - 0.5 + 0.5*bxh, ratio ) )/y
       gamxinflu(i,j,8,2) = (gamma_dpr(real(1.1 - 0.5 + 0.5*bxhl)) - gamxinfdp( 1.1 - 0.5 + 0.5*bxhl, ratio ) )/y
      ENDIF
        
        gamxinflu(i,j,7,2) = gamxinflu(i,j,7,1)
        
        
      ENDDO
      ENDDO
      ciacrratio(0,:) = 1.0

      DO j = ialpstart,nqiacralpha
      alp = float(j)*dqiacralpha
      y = gamma_dpr(real(4.+alp))
      y7 = gamma_dpr(real(7.+alp))
      DO i = 0,nqiacrratio
        ratio = float(i)*dqiacrratio

        ! mass fraction
        x = gamxinfdp( 4.+alp, ratio )
!        write(0,*) 'i, x/y = ',i, x/y
        qiacrratio(i,j) = x/y

        gamxinflu(i,j,4,1) = x/y
        gamxinflu(i,j,4,2) = x/y

        ! reflectivity fraction
        x = gamxinfdp( 7.+alp, ratio )
        ziacrratio(i,j) = x/y7
        gamxinflu(i,j,11,1) = x/y7
        gamxinflu(i,j,11,2) = x/y7
      ENDDO
      ENDDO
      qiacrratio(0,:) = 1.0

   END SUBROUTINE MAKEqiacrratio
      

! #######################################################################

 END MODULE MICRO_MODULE

