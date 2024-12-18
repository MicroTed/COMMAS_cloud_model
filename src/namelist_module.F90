!-----------------------------------------------------------------------------------------------------------------------------------
!  
!  TRMM MODULE:  parameters for subdomain analysis
!

MODULE TRMM_MODULE
   
  implicit none
  
  integer :: gridtimes(100) = 0
   
  integer :: gridx1(100) = 0
  integer :: gridx2(100) = 0
  integer :: gridy1(100) = 0
  integer :: gridy2(100) = 0
  integer :: ibsd,iesd,jbsd,jesd  ! b_egin, e_nd of x (i) and y (j) subdomain
  logical :: lpredict = .true.
  character(LEN = 255) :: trmmoutfile = 'none', trmmstatfile = 'none'
   
  NAMELIST /trmm/                 &
                   gridtimes,     &
                   gridx1,        &
                   gridx2,        &
                   gridy1,        &
                   gridy2,        &
                   trmmoutfile,   &
                   trmmstatfile
                   
  CONTAINS
  
   SUBROUTINE TRMM_MODULE_INIT()
   END SUBROUTINE TRMM_MODULE_INIT
                    
END MODULE TRMM_MODULE

!-----------------------------------------------------------------------------------------------------------------------------------
!                 
! ENS_INIT NAMELIST
!

MODULE ENS_INIT_NML

  implicit none
  
  integer              :: ens_init_type                 = -1
  integer              :: seed                          = 2147483562
  character(LEN = 100) :: ens_init_file                 = 'none'
  real                 :: ens_init_file_time            = 0.0
  real                 :: mean_u_profile_pert           = 0.0
  real                 :: mean_v_profile_pert           = 0.0
  real                 :: mean_t_profile_pert           = 0.0
  real                 :: mean_td_profile_pert          = 0.0
  real                 :: member_u_profile_pert         = 0.0
  real                 :: member_v_profile_pert         = 0.0
  real                 :: member_t_profile_pert         = 0.0
  real                 :: member_td_profile_pert        = 0.0

  real                 :: mean_u_profile_pert_bottom    = 0.0
  real                 :: mean_u_profile_pert_top       = 0.0
  real                 :: mean_v_profile_pert_bottom    = 0.0
  real                 :: mean_v_profile_pert_top       = 0.0
  real                 :: mean_t_profile_pert_bottom    = 0.0
  real                 :: mean_t_profile_pert_top       = 0.0
  real                 :: mean_td_profile_pert_bottom   = 0.0
  real                 :: mean_td_profile_pert_top      = 0.0

  real                 :: mean_qv_profile_pert_bottom   = 0.0
  real                 :: mean_qv_profile_pert_top      = 0.0

  real                 :: member_u_profile_pert_bottom  = 0.0
  real                 :: member_u_profile_pert_top     = 0.0
  real                 :: member_v_profile_pert_bottom  = 0.0
  real                 :: member_v_profile_pert_top     = 0.0
  real                 :: member_t_profile_pert_bottom  = 0.0
  real                 :: member_t_profile_pert_top     = 0.0
  real                 :: member_td_profile_pert_bottom = 0.0
  real                 :: member_td_profile_pert_top    = 0.0

  real                 :: member_qv_profile_pert_bottom = 0.0
  real                 :: member_qv_profile_pert_top    = 0.0

  integer              :: pertmethod1d                  = 3     ! Type of perts added to initial sounding of ensemble members 
                                                                !    using 'member_x_profile_pert':
                                                                !    1 = uncorrelated random noise (eigenpert)
                                                                !    2 = cosine perturbations 
                                                                !    3 = uncorrelated random noise (old version)
                                                                !    4 = random noise
  integer              :: n_levels_to_perturb           = -1    ! How many levels (starting at ground) to perturb in sounding 
                                                                ! used with pertmethod1d, value of -1 causes reset to nz-1
  integer              :: fix_adiabatic_layers          = 1
  integer              :: fix_saturated_layers          = 0
  integer              :: pertmethod3d                  = 2     ! type of perturbations added to 3D fields of ensemble members 
                                                                !    1 = random noise
                                                                !    2 = Gaussian bubbles

  real                 :: rbubh                  = 7500.
  real                 :: rbubv                  = 2000.
  integer              :: nb                     = 2
  real                 :: upert                  = 0.0
  real                 :: vpert                  = 0.0
  real                 :: wpert                  = 0.0
  real                 :: tpert                  = 5.0
  real                 :: qpert                  = 0.0
  real                 :: tdpert                 = 0.0  
  real                 :: xbmin                  = -105000.0      ! Region to add bubbles, coordinates are geo-reference
  real                 :: xbmax                  = -55000.0       ! Region to add bubbles, coordinates are geo-reference
  real                 :: ybmin                  = -60000.0       ! Region to add bubbles, coordinates are geo-reference
  real                 :: ybmax                  = 0.0                  ! Region to add bubbles, coordinates are geo-reference
  real                 :: zbmin                  = 250.0          ! Region to add bubbles, coordinates are geo-reference
  real                 :: zbmax                  = 2250.0         ! Region to add bubbles, coordinates are geo-reference

  integer              :: isetcnox               = 0
  integer              :: isetccn                = 0
  integer              :: isethl                 = 0
  integer              :: isetehs                = 0
  integer              :: write_profiles         = 0

! old deprecated vars
  integer              :: perttype               = -1
  integer              :: pertmethod             = -1
  integer              :: nzob                   = 0


  NAMELIST /ens_init/ ens_init_type,                  &
                      seed,                           &
                      ens_init_file,                  &
                      ens_init_file_time,             &
                      mean_u_profile_pert_bottom   ,  &
                      mean_u_profile_pert_top      ,  &
                      mean_v_profile_pert_bottom   ,  &
                      mean_v_profile_pert_top      ,  &
                      mean_t_profile_pert_bottom   ,  &
                      mean_t_profile_pert_top      ,  &
                      mean_td_profile_pert_bottom  ,  &
                      mean_qv_profile_pert_bottom,    &
                      mean_qv_profile_pert_top,       &
                      mean_td_profile_pert_top     ,  &
                      member_u_profile_pert_bottom ,  &
                      member_u_profile_pert_top    ,  &
                      member_v_profile_pert_bottom ,  &
                      member_v_profile_pert_top    ,  &
                      member_t_profile_pert_bottom ,  &
                      member_t_profile_pert_top    ,  &
                      member_td_profile_pert_bottom,  &
                      member_td_profile_pert_top   ,  &
                      member_qv_profile_pert_bottom,  &
                      member_qv_profile_pert_top   ,  &
                      pertmethod1d,                   &
                      n_levels_to_perturb,            &
                      fix_adiabatic_layers,           &
                      fix_saturated_layers,           &
                      pertmethod3d,                   &
                      perttype,                       &
                      rbubh,                          &
                      rbubv,                          &
                      nb,                             &
                      nzob,                           &
                      upert,                          &
                      vpert,                          &
                      wpert,                          &
                      tpert,                          &
                      tdpert,                         &
                      qpert,                          &
                      xbmin,                          &
                      xbmax,                          &
                      ybmin,                          &
                      ybmax,                          &
                      zbmin,                          &
                      zbmax,                          &
                      isetcnox,                       &
                      isetccn, isethl,                &
                      isetehs,                        &
                      write_profiles 
                      
  CONTAINS
  
   SUBROUTINE ENS_INIT_INIT()
   END SUBROUTINE ENS_INIT_INIT

END MODULE ENS_INIT_NML

!-----------------------------------------------------------------------------------------------------------------------------------
!
!  Correct_ensemble NAMELIST:
!

MODULE CORRECT_ENS_NML

  implicit none
  
  integer :: add_bubbles           = 1           ! Should thermal bubbles be added in/near where ensemble-mean
                                                 !    reflectivity is too low? (1=yes, 0=no)

  real    :: min_refl_diff         = 30.0        ! Minimum difference between observed and model reflectivity
                                                 !    for deciding where to add bubbles

  real    :: cref_depth            = 18000.      ! Depth in (meters) from which compute composite reflectivity

  integer :: add_pert_in_high_refl = 0           ! Should perturbations be added to high-reflectivity
                                                 !    regions? (1=yes, 0=no)
  integer :: pert_type             = 2           ! Perturbation type:  1=noise, 2=smooth

  integer :: which_reflectivity    = 2           ! 1 = add perturbations where ensemble-mean reflectivity is high
                                                 ! 2 = add perturbations where observed reflectivity is high
                                                 ! 3 = add perturbations where either observed or mean reflect is high
  real :: refl_thresh_for_pert     = 20.0        ! minimum reflectivity threshold (dBZ) for adding perturbations
  real :: lhsmth                   = 4000.       ! horizontal length scale (m) for smooth perturbations
  real :: lvsmth                   = 2000.       ! vertical length scale (m) for smooth perturbations

  integer, parameter :: maxbins    = 20          ! Number of time bins for varying the perturbations

  real, dimension(maxbins+1) :: pert_start_times   =  1.e15
  real, dimension(maxbins+1) :: bubble_start_times =  -1
 
  real :: u_noise(maxbins)         = 0.0         ! std. dev. of u noise (m/s), before smoothing
  real :: v_noise(maxbins)         = 0.0         ! std. dev. of v noise (m/s), before smoothing
  real :: w_noise(maxbins)         = 0.0         ! std. dev. of w noise (m/s), before smoothing
  real :: t_noise(maxbins)         = 0.0         ! std. dev. of t noise (K), before smoothing
  real :: td_noise(maxbins)        = 0.0         ! std. dev. of dewpoint noise (K), before smoothing
  real :: qv_noise(maxbins)        = 0.0         ! std. dev. of qv noise (g/kg), before smoothing
  real :: qc_noise(maxbins)        = 0.0         ! std. dev. of qc noise (g/kg), before smoothing
  real :: qr_noise(maxbins)        = 0.0         ! std. dev. of qr noise (g/kg), before smoothing
  real :: qi_noise(maxbins)        = 0.0         ! std. dev. of qi noise (g/kg), before smoothing
  real :: qs_noise(maxbins)        = 0.0         ! std. dev. of qs noise (g/kg), before smoothing
  real :: qh_noise(maxbins)        = 0.0         ! std. dev. of qh noise (g/kg), before smoothing
  
  real :: cor_xbmin                = -500000.0   ! region to add noise or bubbles (default all domain)
  real :: cor_xbmax                =  500000.0   ! region to add noise or bubbles (default all domain)
  real :: cor_ybmin                = -500000.0   ! region to add noise or bubbles (default all domain)
  real :: cor_ybmax                =  500000.    ! region to add noise or bubbles (default all domain)
  real :: cor_zbmin                =  250.       ! region to add noise or bubbles (default boundary layer) 
  real :: cor_zbmax                =  2500.      ! region to add noise or bubbles (default boundary layer)

  real :: noise_zmin                =  0.        ! region to add noise (default boundary layer) 
  real :: noise_zmax                =  12000.    ! region to add noise (default boundary layer)
  
  real :: cor_tpert                = 2.0         ! TH-perturbation for correct_ensemble (not often used)
  real :: cor_rbubh                = 7500.       ! Horiz. radius of bubble perturbations
  real :: cor_rbubv                = 2000.       ! Vertical radius of bubble perturbations

  real :: t_correct             = -100.  ! magnitude of theta perturbation (K) used in correct_ensemble
  real :: rh_correct            = -100.  ! horizontal radius (m) of bubbles
  real :: rv_correct            = -100.  ! vertical radius (m) of bubbles

  namelist /correct_ens/ add_bubbles,             &
                         min_refl_diff,           &
                         add_pert_in_high_refl,   &
                         pert_type,               &
                         cref_depth,              &
                         which_reflectivity,      &
                         refl_thresh_for_pert,    &
                         pert_start_times,        &
                         bubble_start_times,      &
                         lhsmth,                  &
                         lvsmth,                  &
                         cor_tpert,               &
                         cor_rbubh,               &
                         cor_rbubv,               &
                         cor_xbmin,               &
                         cor_xbmax,               &
                         cor_ybmin,               &
                         cor_ybmax,               &
                         cor_zbmin,               &
                         cor_zbmax,               &
                         noise_zmin,              &
                         noise_zmax,              &
                         u_noise,                 &
                         v_noise,                 &
                         w_noise,                 &
                         t_noise,                 &
                         td_noise,                &
                         qv_noise,                &
                         qc_noise,                &
                         qr_noise,                &
                         qi_noise,                &
                         qs_noise,                &
                         qh_noise

END MODULE CORRECT_ENS_NML

MODULE ENKF_NML

  implicit none
  

!-----------------------------------------------------------------------------------------------------------------------------------
!
!
!  ENKF NAMELIST
!
! write_anal     = 0               # Flag to determine whether the forcast and analysis fields are written to the file
! unfold_vr = 1                    # should Doppler velocities be unfolded?  (1=yes, 0=no)
! vr_diff_reject = 999.9           # Doppler velocity rejection criterion for maximum
!                                  #    difference between ensemble mean and observation
! map_proj = 0                     # map projection (for relating lat, lon to x, y):
!                                  #    0=flat earth
! cutoff = 2                       # localization cutoff type:  1=sharp, 2=smooth
! zcutoff = 19000.                  # Altitude cut-off above which obs are not assimilated
! rhoriz = 6000.0                  # horizontal cutoff radius (m) for covariance estimation
! rvert = 6000.0                   # vertical cutoff radius (m) for covariance estimation
! rhoriz2 = 2000.0                 # horizontal cutoff radius (m) for low-reflectivity observations
! rvert2 = 2000.0                  # vertical cutoff radius (m) for low-reflectivity observations
! reflflag1 = 1                    # reflectivity processing method:  1=standard EnKF,
!                                  #    2=standard EnKF for high values, special processing for low values
! reflflag2 = 0                    # replace bad/missing reflectivity values with mindbz? 0=no, 1=yes
! lowrefonly = 0                   # assimilate only low reflectivity values? 0=no, 1=yes
! inflateb = 1.0                   # inflation factor, before filter
! inflated = 1.0                   # inflation factor, during filter
! inflatea = 1.0                   # inflation factor, after filter

  integer :: unfold_vr       = 1       !
  real    :: vr_diff_reject  = 999.9   !
  integer :: map_proj        = 0       !
  real    :: cor_thresh      = 0.0     !  Minimum abs. value of correlation coeff. to allow update
  integer :: cutoff          = 2       !
  real    :: zcutoff         = 19000.0 !
  real    :: rzcutoff        = 1000.0   ! distance for 1/e exponential damping of 'af' above zcutoff
  real    :: rhoriz          = 6000.0  !
  real    :: rvert           = 6000.0  !
  real    :: rhoriz2         = 2000.0  !
  real    :: rvert2          = 2000.0  !
  integer :: reflflag1       = 1       !
  integer :: reflflag2       = 0       !
  integer :: lowrefonly      = 0       !
  real    :: startlowrefonly = -1.     !
  integer :: write_fcst      = 0       !
  integer :: write_anal      = 0       !
  real    :: inflateb        = 1.0     !
  real    :: inflated        = 1.0     !
  real    :: inflatea        = 1.0     !

  real    :: lowdbz          = 10.0
  
  real    :: q_pwr = 0.0, c_pwr = 0.0, z_pwr = 0.0  ! exponents for scaling q,c,z in the filter

  logical :: write_ensemble  = .true.
  logical :: parallel_algorithm = .false.
  
  integer :: vr_option = 1             ! Option for estimating vertical component of
                                       ! Vr from reflectivity
                                       ! = 0 : old method of wt(dBZ)
                                       ! = 1 : Use (or calculate) power-weight 
                                       !       fall velocity (VZF)

  integer :: glm_obs_op_type = 1      ! 1 = graupel volume
                                      ! 2 = ice water path
  integer :: glm_obs_op_option = 1    ! options within type (see obs_op subroutine code)
  real    :: glm_obs_op_dx = 8000.    ! half-width of box for GLM obs_op

  integer :: update_dbzvzf = 0        ! whether to update the DBZ and VZF fields or not
                                      ! 0 = do not update DBZ/VZF fields
                                      ! 1 = update DBZ/VZF fields and set ob_parallel_algorithm = .true.
                                      ! -1 = update DBZ/VZF fields and NOT set ob_parallel_algorithm flag
                                      ! 2 = ONLY FOR TESTING: set ob_parallel_algorithm = .true. but DO NOT update DBZ/VZF

  integer :: dbz_shutoff_time = -1    ! time (model time) to stop assimilating DBZ. Switches dbz obs to verification
  
  namelist /enkf_param/ unfold_vr,       &
                        vr_diff_reject,  &
                        map_proj,        &
                        cor_thresh,      &
                        zcutoff,         &
                        rzcutoff,        &
                        cutoff,          &
                        rhoriz,          &
                        rvert,           &
                        rhoriz2,         &
                        rvert2,          &
                        reflflag1,       &
                        reflflag2,       &
                        lowrefonly,      &
                        startlowrefonly, &
                        inflateb,        &
                        inflated,        &
                        inflatea,        &
                        lowdbz,          &
                        write_fcst,      &
                        write_anal,      &
                        write_ensemble,  &
                        parallel_algorithm, &
                        vr_option,       &
                        q_pwr,           &
                        c_pwr,           &
                        z_pwr,           &
                        glm_obs_op_type, &
                        glm_obs_op_option,&
                        glm_obs_op_dx,    &
                        update_dbzvzf,    &
                        dbz_shutoff_time


END MODULE ENKF_NML

!-----------------------------------------------------------------------------------------------------------------------------------
!
! Synthob NAMELIST

MODULE SYNTHOB_NML

  implicit none
  
  integer, parameter :: n_sfc_ob_types = 4  ! # of different surface ob types in namelist
  
  integer nxyz3dtruth                       ! number of fields in true model state
  real    refl_threshold_for_vr
  integer radar_loc_flag                    ! 1=location specified by radar_lat and radar_lon, 2=U observed, 3=V observed
  real    radar_lat                         ! radar latitude (deg)
  real    radar_lon                         ! radar longitude (deg)
  real    rand_error_refl, rand_error_vr
  real    bias_error_refl, bias_error_vr
  integer :: sfc_ob_types(n_sfc_ob_types) = 0        ! 1 if the following observation type should be created; otherwise 0:
                                                     !   (1) u10m, (2) v10m, (3) temp2m, and (4) qv2m
  integer :: sfc_ob_distribution(2)       = 0        ! sfc_ob_distribution(1).eq.0 indicates gridded locations,
                                                     !   on a grid sfc_ob_distribution(2) X sfc_ob_distribution(2) points
                                                     ! sfc_ob_distribution(1).eq.1 indicates random locations,
                                                     !   with sfc_ob_distribution(2) total points
  real ::    sfc_ob_xloc(2)              = 0.0       ! surface obs are created for sfc_ob_xloc(1) <= x <= sfc_ob_xloc(2)
  real ::    sfc_ob_yloc(2)              = 0.0       ! surface obs are created for sfc_ob_yloc(1) <= y <= sfc_ob_yloc(2)
  real ::    sfc_ob_rand_error(n_sfc_ob_types) = 1.0 ! random errors for the following observation types:
                                                     !   (1) u10m, (2) v10m, (3) temp2m, and (4) qv2m
  real ::    sfc_ob_bias_error(n_sfc_ob_types) = 1.0 ! bias errors for the following observation types:
                                                     !   (1) u10m, (2) v10m, (3) temp2m, and (4) qv2m

  real ::    glm_lightning_dx = 8000.       !  horizontal resolution (m) of synthetic GLM lightning flash rate
  real ::    glm_lightning_z  = 6500.       !  assumed altitude (m) of synthetic GLM lightning flash rate
  real ::    glm_rand_error   = 0.
  real ::    glm_bias_error   = 0.
  integer :: glm_averagingwindow   = 1           ! how many time levels to average over. Assumes centered for odd, "upwind in time" for even
  
  character(len=100) :: flashfile = 'none'

  real dbz_thres
  integer radar_sample_flag                 ! 1=Point sampling, 2=Volumetric sampling, 3=Simple
   
  namelist /synthob_param/ nxyz3dtruth,            &
                          refl_threshold_for_vr,   &
                          radar_loc_flag,          &
                          radar_lat,               &
                          radar_lon,               &
                          rand_error_refl,         &
                          rand_error_vr,           &
                          bias_error_refl,         &
                          bias_error_vr,           &
                          sfc_ob_types,            &
                          sfc_ob_distribution,     &
                          sfc_ob_xloc,             &
                          sfc_ob_yloc,             &
                          sfc_ob_rand_error,       &
                          sfc_ob_bias_error,       &
                          glm_lightning_dx,        &
                          glm_rand_error,          &
                          glm_bias_error,          &
                          glm_lightning_z,         &
                          dbz_thres,               &
                          radar_sample_flag,       &
                          flashfile,               &
                          glm_averagingwindow

   real    beamwidth                         ! half-power beamwidth
   real    eff_bw                            ! effective beamwidth
   real    azim_interval                     ! Azimuth sampling interval
   real    range_interval                    ! gate lenth (range sampling interval )
   integer vcp_num                           ! 11: Storm mode 14 tilts, 12: Precip mode 9 tilts
   integer pts_az, pts_el, pts_rg
   integer samp_az, samp_el, samp_rg
   integer sweeps

   namelist /rad_vol_param/ beamwidth,               &
                            eff_bw,                  &
                            azim_interval,           &
                            range_interval,          &
                            vcp_num,                 &
                            sweeps,                  &
                            samp_az,                 &
                            samp_rg,                 &
                            samp_el,                 &
                            pts_az,                  &
                            pts_el,                  &
                            pts_rg

END MODULE SYNTHOB_NML
!-----------------------------------------------------------------------------------------------------------------------------------
!  
!  RUN_ATTRIBUTE NAMELIST
!

MODULE RUN_ATT_NML

  implicit none
  
  character(LEN = 80) :: prefix
  character(LEN = 80) :: rstprefix = ' '  ! name of run to use for restart, but output new files using 'prefix'
  integer             :: member   = 0
  integer             :: start    = -1
  integer             :: stop     = 1800
  integer             :: new_dt   = 0
  integer             :: thistory = 900
  integer             :: trestart = 3600
  integer             :: tprint   = 60
  integer             :: tstat    = 60
  integer             :: tvis5d   = 999999
  integer             :: tplot    = 999999 ! only used by complot
  real                :: ugrid    = 0.0
  real                :: vgrid    = 0.0
  character(LEN = 15) :: microphys = 'LFO'   ! Default setting is LFO microphysics
  character(LEN = 255):: v5dflds   = '/THETA/THETAP/T/TD/U/V/W/DBZ/QC/QR/'
  integer             :: ne        = 1
  integer             :: year      = 1982
  integer             :: month     = 1
  integer             :: day       = 1
  integer             :: hour      = 00
  integer             :: minute    = 00
  integer             :: second    = 00
  integer             :: v5dstridex = 1  ! used to coarsify vis5d output to allow more time levels of data
  integer             :: v5dstridez = 1  ! used to coarsify vis5d output to allow more time levels of data
  integer             :: itrack = -1
  real                :: track_xmid = 55000., track_ymid = 50000. ! x,y of target location for meso tracking
  real                :: track_xmid_delta = -1., track_ymid_delta = -1.
  real                :: track_tau = 200., track_alpha = 0.75, track_zeta0 = 0.001 ! meso tracking params

  character(LEN=10),dimension(400) :: v5dvarlist = ' ' ! CHAR(0)

  NAMELIST /run/              &      ! prefix,      &
                 member,      &
                 start,       &
                 stop,        &
                 new_dt,      &
                 tprint,      &
                 tstat,       &
                 thistory,    &
                 trestart,    &
                 tvis5d,      &
                 tplot,       &
                 ugrid,       &
                 vgrid,       &
                 itrack, track_xmid, track_ymid, track_tau, track_alpha, track_zeta0, &
                 track_xmid_delta, track_ymid_delta
 
 CONTAINS
 
   SUBROUTINE RUN_ATT_MODULE_INIT()
   END SUBROUTINE RUN_ATT_MODULE_INIT


END MODULE RUN_ATT_NML


!===================================================================================================================================
!
!
! 3D.RUN NAMELIST PARAMETERS
!
!-------------------------------------------------------------------------------
! 
! NAMELIST MODULE:  manages systems namelists, including I/O and info
!

MODULE NAMELIST_MODULE
 
 USE MICRO_MODULE      ! from micro_module.f90
 USE ELEC_MODULE       ! from elec_module.f90
 USE STRING_MODULE     ! from string_module.f90
 USE PARAM_MODULE      ! from param_module.f90
 USE FORCE_MODULE, only: nphys,iforce,iforcetyp,nwfor,iwforce,ivshap,iwshap,iushap,  &    ! from force_module.F90
                        iqshap,   itshap,   isshap, isslcl,   wfmeso, ufmeso, umaxmeso,  xwfcen,   ywfcen,         &
                        zwfcen,   xwfrad,   ywfrad,   zwfrad,   xwfmov,   ywfmov,         &
                        zwfmov,   nssfor,   issforce, ssfmeso,  iqrforce, qrfmeso,        &
                        qhfmeso,  xssfcen,  yssfcen,  rssfrad,  tsstrt,   tsstop,         &
                        tsslow,   tsshigh,  twstrt,   twstop,   itopforce,itopforcetype,  &
                        topforceq,topforceN,topforceX,topforceY,topforceR,topforcedia,    &
                        topforcen0,  topforcealpha,  igamrain,chargeperparticle,          &
                        itopforcerainopt, topforcerainrate, wnaylor, forcerandfac,        &
                   !     prg0,thg0 ,uspeed,timlndst,thg1,isstbeg,isstend,                  &
                        force_init,low_level_cooling_flag ,low_level_cooling_depth,       &
                        low_level_cooling_rate ,low_level_cooling_start,                  &
                        low_level_cooling_end,                                            &
                        iccnufforce,nccnfor,tccnstrt,tccnstop,xccnfcen,yccnfcen,zccnfcen, &
                        xccnfrad,yccnfrad,zccnfrad,                                       &
                        ichgshap, nchgfor, chgratefor, xchgfcen, ychgfcen, zchgfcen,      &
                        xchgfrad, ychgfrad, zchgfrad, ichgforce,                          &
                        xchgcldfrad, ychgcldfrad, zchgcldfrad, xchgcldfcen, ychgcldfcen, zchgcldfcen



 USE FILE_MODULE       ! from file_module.f90
 USE BALLOON_MODULE    ! from balloon.f90
 USE INDEX_MODULE      ! from index_module.f90
 USE INIT_MODULE       ! from initsubs.F90
 USE GRIDPARAM_MODULE  ! from initsubs.F90
 USE COMMASMPI_MODULE  ! from commasmpi_module.f90
 USE TRMM_MODULE       ! from above
 USE ENKF_NML          ! from above
 USE ENS_INIT_NML      ! from above
 USE CORRECT_ENS_NML   ! from above
 USE SYNTHOB_NML       ! from above
 USE RUN_ATT_NML       ! from above
 USE TRAJ_MODULE       ! from traj_module.F90
 USE BSS_NML           ! from box.o
 USE VIS5D_MODULE, only: alphamax_v5d, alphamin_v5d

 implicit none
               
!-----------------------------------------------------------------------------
!                 
! LFO PARAMETERS 
!                 
!-----------------------------------------------------------------------------

!  real                  :: rho_qr         = 1000.     
!  real                  :: cnor           = 8.0e6 
!  real                  :: rho_qs         = 100.
!  real                  :: cnos           = 3.0e6
!  real                  :: rho_qh         = 700.
!  real                  :: cnoh           = 4.0e5
!  integer               :: autoconversion = 1
!
!  NAMELIST /lfo_params/ rho_qr,         &
!                        cnor,           &
!                        rho_qs,         &
!                        cnos,           &
!                        rho_qh,         &
!                        cnoh,           &
!                        autoconversion

!-----------------------------------------------------------------------------
!                 
! PARAM PARAMETERS 
!                 
!-----------------------------------------------------------------------------

  NAMELIST /paramh/ &
            coriol,     &
            drag,       &
            rayd_hgt,   &
            rayd_mag,   &
            hsponge,    &
            hrayd_mag,  &
            lfwds,      &
            detrendpi,  &
            icrwmp,     &
            icrwmn, icrwmn1st, &
            mix_type,   &
            len_type,   &
            tke_type,   &
            iusebvfreq, &
            rmbasekm,   &
            Lmax,       &
            ahighk,     &
            Cm,         &
            Ce,         &
            Pr,         &
            km_thresh,  &
            kmgen_thresh, &
            kmshear_thresh, &
            kmbasefac,  &
            vert_adv_scheme, &
            MONOTONIC,  &
            RKSCHEME,   &
            rkstepping, &
            damph,      &
            bcx,        &
            bcy,        &
            bcz,        &
            atypes1,    &
            atypes2,    &
            atypem1,    &
            atypem2,    &
            atype2qv,   &
            atype2th,   &
            atype2tke,  &
            iwensmooth, &
            dxt,        &
            kdiv0,      &
            nphys,      &
            nocollapse, &
            ainflo,     &
            ainfloqv,   &
            ainflom,    &
            zinflo_u_s, &
            zinflo_u_n, &
            zinflo_v_w, &
            zinflo_v_e, &
            zinflo_w_e, &
            zinflo_w_w, &
            zinflo_w_s, &
            zinflo_w_n, &
            ifilt,      &
            izero_opt,  &
            iuvwadv,    &
            isfcphys,   &
            isfcl,prg0,thg0 ,uspeed,timlndst,   &
            thg1,isstbeg,isstend, &
            vert_implicit, &
            Cspd,          &
            rotunno_radiation

  NAMELIST /mpi_params/ nxt,           &
                        nyt,           &
                        nzt,           &
                        nxprocs,         &
                        nyprocs,         &
                        nzprocs,         &
                        verbose_mpi

  NAMELIST /lfo_params/ rho_qr,       &
                        cnor,         &
                        rho_qs,       &
                        cnos,         &
                        rho_qh,       &
                        cnoh,         &
                        rho_qhl,      &
                        cnohl,        &
                        qcmincwrn,    &
                        cwdiap,       &
                        cwdisp,       &
                        ccn,          &
                        hole_fill,    &
                        autoconversion, &
                        iuseferrier

  NAMELIST /attributes/ prefix,      &
                        rstprefix,   &
                        microphys,   &
                        v5dflds,     &
                        ne,          &
                        member,      &
                        year,        &
                        month,       &
                        day,         &
                        hour,        &
                        minute,      &
                        second,      &
                        v5dvarlist,  &
                        v5dstridex,  &
                        v5dstridez,  &
                        historyoutput, &
                        historyinput,  &
                        restartformat, &
                        parallelio,    &
                        parallel_compress_on, &
                        parallelio_in,    &
                        parallelio_type,  &
                        netcdfversion, &  ! defined in file_module, default value is now 4 (if NC4 is defined, then netcdf4 will be the default format)
                        deflate_level, &
                        deflate,       &
                        restart_separate, &
                        onedoutput

NAMELIST /output_options/               &
                        inetchargetend, &
                        ioutput_xtrachgsep, &
                        ioutput_hwind3d, &
                        ioutput_temC3d,  &
                        ioutput_pres3d,  &
                        ioutput_thv3d,   &
                        ioutput_workshop, &
                        ioutput_flshr8km, &
                        ioutput_ssw,      &
                        ioutput_ssi,      &
                        ioutput_rh,       &
                        ioutput_cwdia,    &
                        ioutput_cnu,      &
                        ioutput_ssmx,     &
                        ioutput_snowstuff,&
                        ioutput_takrates, &
                        ioutput_mltshedsizerates, &
                        ioutput_sedmeltstuff, &
                        ioutput_vzf,      &
                        ioutput_icedensity, &
                        ioutput_tkediss,    &
                        ioutput_dbzsedchange, &
                        alphamax_v5d, alphamin_v5d
  
  NAMELIST /micro_params/               &
                        ndebug, ncdebug,&
                        iptime,         &
                        ipconc, ilimit, &
                        ichaff, chaffconc, inucopt,&
                        imixice,        &
                        icespheres,     &
                        frozendrops, iraintypes,    &
                        iashtypes,      &
                        rho_qr,         &
                        cnor,           &
                        rho_qs,         &
                        cnos,           &
                        rho_qh,         &
                        cnoh,           &
                        rho_qhl,        &
                        cnohl,          &
                        rho_qh_max,     &
                        rho_qhl_max,    &
                        rho_qr_i10,     &
                        cnor_i10,       &
                        rho_qs_i10,     &
                        cnos_i10,       &
                        rho_qf_i10,     &
                        cnof_i10,       &
                        rho_qgl_i10,    &
                        cnogl_i10,      &
                        rho_qgm_i10,    &
                        cnogm_i10,      &
                        rho_qgh_i10,    &
                        cnogh_i10,      &
                        rho_qh_i10,     &
                        cnoh_i10,       &
                        rho_qhl_i10,    &
                        cnohl_i10,      &
                        rho_qh_tak,     &
                        rho_qhl_tak,    &
                        qcmincwrn,      &
                        cwdiap,         &
                        cwdisp,         &
                        ccn,ccnuf,ccnac,&
                        ccnnu, ccnco,   &
                        cin,naer,       &
                        hole_fill,      &
                        autoconversion, &
                        iuseferrier,    &
                        iusewetgraupel, &
                        iusewethail,    &
                        iusewetsnow,    &
                        idbzci,         &
                        cimn, cimx,     &
                        vtmaxsed,       &
                        itfall,iscfall, &
                        infall,infalln, imydiagalpha, &
                        irfall,isfall,iifall,  &
                        rssflg,         &
                        sssflg,         &
                        hssflg,         &
                        hlssflg,        &
                        ireadmic,       &
                        irimdenopt,rimdenvwgt,     &
                        rimc1, rimc2, rimc3, rimc4,   &
                        alfarim,        &
                        qrimmnc,qrimmnp,&
                        idiagnosecnu,   &
                        iccwflg,        &
                        issfilt,        &
                        icnuclimit,     &
                        irenuc, i_uf_or_ccn,  &
                        irenuc3d,       &
                        restoreccn, ccntimeconst, restoreccnfrac, cck, &
                        renucfrac, ssf2kmax,      &
                        luseccn,        &
                        ac_opt, arg_para, nu_kappa, ac_kappa, co_kappa, ac_wthresh, &
!                        xcradmx,        &
                        ciintmx,        &
                        itype1, itype2, &
                        icenucopt,in_freeze_rain_first,      &
                        icfn, ihrn,     &
                        ibfc, iacr, icracr, &
                        icrcev, icracrthresh, &
                        cwfrz2snowfrac, cwfrz2snowratio, &
                        iremoveqwfrz,   &
                        ibfr,           &
                        ibiggopt,       &
                        ibiggsmallrain, &
                        rhofrz, ifrzg,ifiacrg,  &
                        ifrzs,ffrzs,f2h,    &
                        iacrsize,       &
                        cimas0, cimas1, cfnfac, &
                        splintermass,   &
                        ewfac,          &
                        eii0, eii1,     &
                        eri0, esi0,     &
                        ewi_dcmin, ewi_dimin, &
                        eri_cimin, dmincw,     &
                        eii0hl, eii1hl, &
                        ehs0, ehs1,     &
                        ess0, ess1, iessopt,    &
                        esstem1,esstem2, &
                        essrmax, essfrac1, essfrac2, iessec0flag, &
                        ehslfo0,        &
                        ehsfrac,        &
                        ircnw, qminrncw,&
                        iauttim,        &
                        auttim,         &
                        iglcnvi,        &
                        iglcnvs,        &
                        rz,             &  ! rz is deprecated as namelist var -- recalculated in microphysics now
                        alphahacx,      &
                        fconv,          &
                        eqtot,          &
                        rcond, icond, iqis0,   &
                        imeyers5,       &
                        iehw,iefw,iehlw, &
                        ierw,           &
                        iehr0c,iehlr0c, &
                        alphai,         &
                        alphar,         &
                        alphas,         &
                        alphah,         &
                        alphahl,        &
                        dmuh,           &
                        dmuhl,          &
                        cnu,            &
                        iscni,fscni,    &
                        dfrz,           &
                        dmlt, dshd, ished2cld, ivshdgs, &
                        rainfallfac,    &
                        icefallfac,     &
                        snowfallfac,    &
                        ifallsedonly,   &
                        graupelfallfac, &
                        fdfallfac,      &
                        hailfallfac,    &
                        icefallopt,     &
                        icdx,icdxhl,    &
                        axh,bxh,axf,bxf,axhl,bxhl, &
                        cdhmin, cdhmax,       &
                        cdhdnmin, cdhdnmax,   &
                        cdhlmin, cdhlmax,     &
                        cdhldnmin, cdhldnmax, &
                        ihmlt,          &
                        ehimin,         &
                        ehimax,         &
                        ehsmax,         &
                        ecollmx,        &
                        eiw0, esw0,     &
                        ehw0, ehlw0, efw0,   &
                        ehr0, ehlr0, efr0,   &
                        erw0,           &
                        exwmindiam,     &
                        nsplinter,      &
                        lawson_splinter_fac, &
                        iqcinit,        &
                        ssmxinit,      &
                        denscale,       &
                        xvdmx,          &
                        dhmn, dhmx,     &
                        dhlmn, dhlmx,     &
                        dsmx,           &
                        fwms,fwmh,fwmhl,  &
                        ifwmhopt,         &
                        ihxw2rain,        &
                        fwmlarge,         &
                        ifwmfall,         &
                        iturbenhance,     &
                        qsdenmod,qhdenmod, &
                        qsvtmod,          &
                        alphamin,alphamax, &
                        isnwfrac,          &
                        rescale_low_alpha, &
                        rescale_low_alphar, &
                        rescale_low_alphah, &
                        rescale_low_alphahl, &
                        rescale_high_alpha, &
                        ihlcnh, hldia1,iusedw, dwehwmin, dwmin, dwmax, dwtempmin, dg0thresh, &
                        hlcnvtimeconst, &
                        wetgrthtoffset, hailcnvtoffset, &
                        incwet,dwetmin,    &
                        ifddenfac, fddenthresh, &
                        icvhl2h, hldnmn,hdnmn,    &
                        hlcnhdia, hlcnhqmin, &
                        isedonly,           &
                        iresetmoments,      &
                        cxmin, zxmin,       &
                        imurain,            &
                        iferwisventr,       &
                        izwisventr,         &
                        qhdpvdn,            &
                        qhacidn,            &
                        sheddiam,sheddiamlg, &
                        sheddiam0,           &
                        mltdiam1,mltdiam2,mltdiam3,mltdiam4,mltdiam05, &
                        fwmhtmptem, ifwmhtmptemopt, &
                        imaxdiaopt,          &
                        irainbreak, rainbreakfac, draintail, drsmall, qrbrthresh1, qrbrthresh2, &
                        ithompsoncnoh,       &
                        cnohmn,             &
                        ivhmltsoak,         &
                        ioldlimiter,        &
                        igrplfall,          &
                        isnowfall,          &
                        isnowdens,          &
                        ibiggsnow,biggsnowdiam,   &
                        ixtaltype,          &
                        evapfac,meltfac,    &
                        depfac,             &
                        dmrauto,irescalerainopt,dmropt, dmhlopt,          &
                        rescale_tempthresh, rescale_wthresh, &
                        alpharaut,          &
                        takccntype, taknucsize,        &
                        takcoagopt, takaggopt, takbottopt, takbottrateonly, &
                        takbottoptice, takgminflux, &
                        takshedsmall, takbreakup, taksponbreaksize, taksponbreakmax, &
                        takshedsize1, takshedsize2,takshedsize3,numshedregimes, &
                        takikcollmax,takess0,takess1, &
                        takrestorcn,    &
                        takbiggopt,     &
                        takifrzg,       &
                        takgrmassopt,   &
                        takvf,takkfmax, &
                        takvfw,         &
                        takrw1,         &
                        takviopt,       &
                        takcxmin,       &
                        takrhovtopt,    &
                        takcondsndt,    &
                        takdodmwdt,     &
                        takdocrim, takdocondsn, takdocmelt, takdocollsn, takdosed, &
                        takessopt,      &
                        takrventopt,    &
                        takhventopt,    &
                        takiventopt,    &
                        takicenucopt,   &
                        takicenucsize,  &
                        takfdiopt,      &
                        takfdirmax,     &
                        takicethickopt, &
                        takicethickness,&
                        takdieff,       &
                        takcwnucopt,    &
                        takhallettopt,  &
                        takgrhlcv,takhlgrcv,      &
                        taknicfd,                 &
                        takehiopt,                &
                        takfrzopt,                &
                        takrhoi,                  &
                        taksri0,                  &
                        itakgrowth,itakmeltgrowth,&
                        taknewshed,takshedopt,takrimeffopt,    &
                        idichgopt,                &
                        takhallettsize,           &
                        ibinhmlr,ibinhlmlr,imltshddmr, binmlrmxdia, binmlrzrrfac,ibinnum,   &
                        ibincracr,      &
                        iqhacrmlr, iqhlacrmlr, iqhacwshr, iqhlacwshr, &
                        snowmeltdia, alphasmlr0,   &
                        delta_alphamlr, &
                        ibinhacw, chmlrmult,      &
                        iqvsopt,     &
                        wrfccn, wrfdefaults,        &
                        usenucond,      &
                        usensslgs,      &
                        useoldgs ,      &
                        imaxsupopt, maxsupersat,ssmxuf,    &
                        ihailrain,      &
                        morr_rimed_ice, morr_droplet_conc, morr_dnr_max, &
                        imorrgdnglimit,morrdnglimit, &
                        thom_autoconv_fac, &
                        iwetsoak

  NAMELIST /nssl_mp_params/               &
                        ndebug, ncdebug,&
                        iusewetgraupel, &
                        iusewethail,    &
                        iusewetsnow,    &
                        idbzci,         &
                        vtmaxsed,       &
                        itfall,iscfall, &
                        infall,infalln, imydiagalpha, &
                        irfall,isfall,iifall, &
                        rssflg,         &
                        sssflg,         &
                        hssflg,         &
                        hlssflg,        &
                        irimdenopt,rimdenvwgt,     &
                        rimc1, rimc2, rimc3, rimc4,   &
                        idiagnosecnu,   &
                        icnuclimit,     &
                        irenuc, i_uf_or_ccn, &
                        restoreccn, ccntimeconst, restoreccnfrac, cck, &
                        ciintmx,        &
                        itype1, itype2, &
                        icenucopt,in_freeze_rain_first,      &
                        naer,           &
                        icfn,           &
                        ibfc, iacr, icracr, &
                        icrcev, icracrthresh, &
                        cwfrz2snowfrac, cwfrz2snowratio, &
                        ibfr,           &
                        ibiggopt,       &
                        ibiggsmallrain, &
                        ifrzg,ifiacrg,  &
                        ifrzs,ffrzs,    &
                        iacrsize,       &
                        cimas0, cimas1, cfnfac, &
                        splintermass,   &
                        ewfac,          &
                        eii0, eii1,     &
                        eri0, esi0,     &
                        ewi_dcmin, ewi_dimin, &
                        eri_cimin, dmincw,     &
                        eii0hl, eii1hl, &
                        ehs0, ehs1,     &
                        ess0, ess1, iessopt,    &
                        esstem1,esstem2, &
                        ircnw, qminrncw,& ! single-moment only
                        iglcnvi,        &
                        iglcnvs,        &
                        alphahacx,      &
                        fconv,          &
                        eqtot,          &
                        imeyers5,       &
                        iehw,iefw,iehlw, &
                        ierw,           &
                        iehr0c,iehlr0c, &
                        alphai,         &
                        alphar,         &
                        alphas,         & ! note that alphah and alphahl come through physics namelist
                        cnu,            &
                        iscni,fscni,    &
                        dfrz,           &
                        dmlt,           &
                        rainfallfac,    &
                        icefallfac,     &
                        snowfallfac,    &
                        graupelfallfac, &
                        fdfallfac,      &
                        hailfallfac,    &
                        icefallopt,     &
                        icdx,icdxhl,    &
                        axh,bxh,axf,bxf,axhl,bxhl, &
                        cdhmin, cdhmax,       &
                        cdhdnmin, cdhdnmax,   &
                        cdhlmin, cdhlmax,     &
                        cdhldnmin, cdhldnmax, &
                        ihmlt,          &
                        ehimin,         &
                        ehimax,         &
                        ehsmax,         &
                        ecollmx,        &
                        eiw0, esw0,     &
                        ehw0, ehlw0, erw0,   &
                        ehr0, ehlr0, efr0,   &
                        erw0,           &
                        exwmindiam,     &
                        nsplinter,      &
                        lawson_splinter_fac, &
                        iqcinit,        &
                        ssmxinit,      &
                        xvdmx,          &
                        dhmn, dhmx,     &
                        dhlmn, dhlmx,     &
                        fwms,fwmh,fwmhl,  &
                        ifwmhopt,         &
                        ihxw2rain,        &
                        fwmlarge,         &
                        ifwmfall,         &
                        iturbenhance,     &
                        qsdenmod,qhdenmod, &
                        qsvtmod,          &
                        alphamin,alphamax, &
                        isnwfrac,          &
                        rescale_low_alpha, &
                        rescale_low_alphar, &
                        rescale_low_alphah, &
                        rescale_low_alphahl, &
                        rescale_high_alpha, &
                        ihlcnh, hldia1,iusedw, dwehwmin, dwmin, dwmax, dwtempmin, &
                        hlcnvtimeconst, &
                        incwet,dwetmin,    &
                        ifddenfac, fddenthresh, &
                        icvhl2h, hldnmn,hdnmn,    &
                        hlcnhdia, hlcnhqmin, &
                        isedonly,           &
                        iresetmoments,      &
                        cxmin, zxmin,       &
                        imurain,            &
                        iferwisventr,       &
                        izwisventr,         &
                        qhdpvdn,            &
                        qhacidn,            &
                        sheddiam,sheddiamlg, &
                        sheddiam0,           &
                        mltdiam1,mltdiam2,mltdiam3,mltdiam4,mltdiam05, &
                        fwmhtmptem, ifwmhtmptemopt, &
                        imaxdiaopt,          &
                        irainbreak, rainbreakfac, draintail, drsmall, qrbrthresh1, qrbrthresh2, &
                        ithompsoncnoh,       &
                        cnohmn,             &
                        ivhmltsoak,         &
                        ioldlimiter,        &
                        isnowfall,          &
                        isnowdens,          &
                        ibiggsnow,          &
                        ixtaltype,          &
                        evapfac,    &
                        depfac,             &
                        dmrauto,irescalerainopt, dmropt,dmhlopt,     &
                        rescale_tempthresh, rescale_wthresh, &
                        ibinhmlr,ibinhlmlr,imltshddmr, binmlrmxdia, binmlrzrrfac,ibinnum,   &
                        ibincracr,      &
                        iqhacrmlr, iqhlacrmlr, iqhacwshr, iqhlacwshr, &
                        snowmeltdia, alphasmlr0

  
  NAMELIST /ice10_params/               &
                        ndebug, ncdebug,&
                        iptime,         &
                        ipconc, ilimit, &
                        ichaff, chaffconc, inucopt,&
                        imixice,        &
                        icespheres,     &
                        frozendrops, iraintypes,    &
                        iashtypes,      &
                        rho_qr,         &
                        cnor,           &
                        rho_qs,         &
                        cnos,           &
                        rho_qh,         &
                        cnoh,           &
                        rho_qhl,        &
                        cnohl,          &
                        rho_qh_max,     &
                        rho_qhl_max,    &
                        rho_qr_i10,     &
                        cnor_i10,       &
                        rho_qs_i10,     &
                        cnos_i10,       &
                        rho_qf_i10,     &
                        cnof_i10,       &
                        rho_qgl_i10,    &
                        cnogl_i10,      &
                        rho_qgm_i10,    &
                        cnogm_i10,      &
                        rho_qgh_i10,    &
                        cnogh_i10,      &
                        rho_qh_i10,     &
                        cnoh_i10,       &
                        rho_qhl_i10,    &
                        cnohl_i10,      &
                        rho_qh_tak,     &
                        rho_qhl_tak,    &
                        qcmincwrn,      &
                        cwdiap,         &
                        cwdisp,         &
                        ccn,ccnuf,ccnac,&
                        ccnnu, ccnco,   &
                        cin,naer,       &
                        hole_fill,      &
                        autoconversion, &
                        iuseferrier,    &
                        iusewetgraupel, &
                        iusewethail,    &
                        iusewetsnow,    &
                        idbzci,         &
                        cimn, cimx,     &
                        vtmaxsed,       &
                        itfall,iscfall, &
                        infall,infalln,imydiagalpha, &
                        irfall,isfall,iifall,  &
                        rssflg,         &
                        sssflg,         &
                        hssflg,         &
                        hlssflg,        &
                        ireadmic,       &
                        irimdenopt,rimdenvwgt,     &
                        rimc1, rimc2, rimc3, rimc4,   &
                        alfarim,        &
                        qrimmnc,qrimmnp,&
                        idiagnosecnu,   &
                        iccwflg,        &
                        issfilt,        &
                        icnuclimit,     &
                        irenuc, i_uf_or_ccn, &
                        irenuc3d,       &
                        restoreccn, ccntimeconst, restoreccnfrac, cck, &
                        renucfrac, ssf2kmax,     &
                        luseccn,        &
                        ac_opt, arg_para, nu_kappa, ac_kappa, co_kappa, ac_wthresh, &
!                        xcradmx,        &
                        ciintmx,        &
                        itype1, itype2, &
                        icenucopt,in_freeze_rain_first,      &
                        icfn, ihrn,     &
                        ibfc, iacr, icracr, &
                        cwfrz2snowfrac, cwfrz2snowratio, &
                        iremoveqwfrz,   &
                        ibfr,           &
                        ibiggopt,       &
                        ibiggsmallrain, &
                        rhofrz, ifrzg,ifiacrg,  &
                        ifrzs,ffrzs,f2h,    &
                        iacrsize,       &
                        cimas0, cimas1, cfnfac, &
                        splintermass,   &
                        ewfac,          &
                        eii0, eii1,     &
                        eri0, esi0,     &
                        ewi_dcmin, ewi_dimin, &
                        eri_cimin, dmincw,      &
                        eii0hl, eii1hl, &
                        ehs0, ehs1,     &
                        ess0, ess1, iessopt,    &
                        esstem1,esstem2, &
                        essrmax, essfrac1, essfrac2, iessec0flag, &
                        ehslfo0,        &
                        ehsfrac,        &
                        ircnw, qminrncw,&
                        iauttim,        &
                        auttim,         &
                        iglcnvi,        &
                        iglcnvs,        &
                        rz,             &  ! rz is deprecated as namelist var -- recalculated in microphysics now
                        alphahacx,      &
                        fconv,          &
                        eqtot,          &
                        rcond, icond, iqis0,   &
                        imeyers5,       &
                        iehw,iefw,iehlw, &
                        ierw,           &
                        iehr0c,iehlr0c, &
                        alphai,         &
                        alphar,         &
                        alphas,         &
                        alphah,         &
                        alphahl,        &
                        dmuh,           &
                        dmuhl,          &
                        cnu,            &
                        iscni,fscni,    &
                        dfrz,           &
                        dmlt, dshd, ished2cld, ivshdgs, &
                        rainfallfac,    &
                        icefallfac,     &
                        snowfallfac,    &
                        ifallsedonly,   &
                        graupelfallfac, &
                        fdfallfac,      &
                        hailfallfac,    &
                        icefallopt,     &
                        icdx,icdxhl,    &
                        axh,bxh,axf,bxf,axhl,bxhl, &
                        cdhmin, cdhmax,       &
                        cdhdnmin, cdhdnmax,   &
                        cdhlmin, cdhlmax,     &
                        cdhldnmin, cdhldnmax, &
                        ihmlt,          &
                        ehimin,         &
                        ehimax,         &
                        ehsmax,         &
                        ecollmx,        &
                        eiw0, esw0,     &
                        ehw0, ehlw0, efw0,   &
                        ehr0, ehlr0, efr0,   &
                        erw0,           &
                        exwmindiam,     &
                        nsplinter,      &
                        lawson_splinter_fac, &
                        iqcinit,        &
                        ssmxinit,       &
                        denscale,       &
                        xvdmx,          &
                        dhmn, dhmx,     &
                        dhlmn, dhlmx,     &
                        dsmx,           &
                        fwms,fwmh,fwmhl, &
                        ifwmhopt,         &
                        ihxw2rain,        &
                        fwmlarge,        &
                        ifwmfall,         &
                        iturbenhance,    &
                        qsdenmod,qhdenmod, &
                        qsvtmod,          &
                        alphamin,alphamax, &
                        isnwfrac,          &
                        rescale_low_alpha, &
                        rescale_low_alphar, &
                        rescale_low_alphah, &
                        rescale_low_alphahl, &
                        rescale_high_alpha, &
                        ihlcnh, hldia1, iusedw, dwehwmin, dwmin, dwmax, dwtempmin, dg0thresh, &
                        hlcnvtimeconst, &
                        wetgrthtoffset, hailcnvtoffset, &
                        incwet,dwetmin,      &
                        ifddenfac, fddenthresh, &
                        icvhl2h, hldnmn,hdnmn,    &
                        hlcnhdia, hlcnhqmin, &
                        isedonly,           &
                        iresetmoments,      &
                        cxmin, zxmin,       &
                        imurain,            &
                        iferwisventr,       &
                        izwisventr,         &
                        qhdpvdn,            &
                        qhacidn,            &
                        sheddiam,sheddiamlg, &
                        sheddiam0,           &
                        mltdiam1,mltdiam2,mltdiam3,mltdiam4,mltdiam05, &
                        fwmhtmptem, ifwmhtmptemopt, &
                        imaxdiaopt,          &
                        irainbreak, rainbreakfac, draintail, drsmall, qrbrthresh1, qrbrthresh2, &
                        ithompsoncnoh,       &
                        cnohmn,             &
                        ivhmltsoak,         &
                        ioldlimiter,        &
                        igrplfall,          &
                        isnowfall,          &
                        isnowdens,          &
                        ibiggsnow,biggsnowdiam,   &
                        ixtaltype,          &
                        evapfac,meltfac,    &
                        depfac,             &
                        dmrauto,irescalerainopt, dmropt, dmhlopt,          &
                        rescale_tempthresh, rescale_wthresh, &
                        alpharaut,          &
                        takccntype,taknucsize,         &
                        takcoagopt, takaggopt, takbottopt, takbottrateonly, &
                        takbottoptice, takgminflux, &
                        takshedsmall, takbreakup, taksponbreaksize, taksponbreakmax, &
                        takshedsize1, takshedsize2,takshedsize3,numshedregimes, &
                        takikcollmax,takess0,takess1, &
                        takrestorcn,    &
                        takbiggopt,     &
                        takifrzg,       &
                        takgrmassopt,   &
                        takvf,takkfmax, &
                        takvfw,         &
                        takrw1,         &
                        takviopt,       &
                        takcxmin,       &
                        takrhovtopt,    &
                        takcondsndt,    &
                        takdodmwdt,     &
                        takdocrim, takdocondsn, takdocmelt, takdocollsn, takdosed, &
                        takessopt,      &
                        takrventopt,    &
                        takhventopt,    &
                        takiventopt,    &
                        takicenucopt,   &
                        takicenucsize,  &
                        takfdiopt,      &
                        takfdirmax,     &
                        takicethickopt, &
                        takicethickness,&
                        takdieff,       &
                        takcwnucopt,    &
                        takhallettopt,  &
                        takgrhlcv,takhlgrcv,      &
                        taknicfd,                 &
                        takehiopt,                &
                        takfrzopt,                &
                        takrhoi,                  &
                        taksri0,                  &
                        itakgrowth,itakmeltgrowth,&
                        taknewshed,takshedopt,takrimeffopt,    &
                        idichgopt,                &
                        takhallettsize,           &
                        ibinhmlr,ibinhlmlr,imltshddmr, binmlrmxdia, binmlrzrrfac, ibinnum,   &
                        ibincracr,      &
                        iqhacrmlr, iqhlacrmlr, iqhacwshr, iqhlacwshr, &
                        snowmeltdia, alphasmlr0,    &
                        delta_alphamlr, &
                        ibinhacw, chmlrmult,      &
                        iqvsopt,     &
                        wrfccn, wrfdefaults,        &
                        usenucond,      &
                        usensslgs,      &
                        useoldgs,       &
                        imaxsupopt, maxsupersat,ssmxuf,    &
                        ihailrain,      &
                        morr_rimed_ice, morr_droplet_conc, morr_dnr_max, morr_dbzsoak, &
                        imorrgdnglimit,morrdnglimit, &
                        thom_autoconv_fac, &
                        iwetsoak

  
  NAMELIST /my_params/                  &
                        ntc_my,         &
                        cnor_my,        &
                        cnos_my,        &
                        cnog_my,        &
                        cnoh_my,        &
                        rho_qr_my,      &
                        rho_qi_my,      &
                        rho_qs_my,      &
                        rho_qg_my,      &
                        rho_qh_my,      &
                        alphar_my,      &
                        alphai_my,      &
                        alphas_my,      &
                        alphag_my,      &
                        alphah_my

  NAMELIST /elec_params/                &
                        ipelec,         &
                        charging_border,&
                        icorona,iondriftorder,do_ionatt_conduction, do_ionatt_diffusion, &
                        ionatt_conduction_factor, ionatt_diffusion_factor, &
                        emaxc,          &
                        nliter,         &
                        ilight,         &
                        nonigrd,        &
                        ftauopt,        &
                        nic_noliq,      &
                        qc_noliq,       &
                        ieopt,          &
                        drake,          &
                        ihallettsaund,   &
                        xhallettsaund,   &
                        rgard,          &
                        isaund,         &
                        rarfac,         &
                        iraropt,        &
                        trever,         &
                        cidiamin,       &
                        erbnd,          &
                        costhe,         &
                        delqnsa,delqxsa,&
                        delqnsb,delqxsb,&
                        delqnia,delqxia,&
                        delqxw,delqnw,  &
                        scxacymax,      &
                        elgtthx,        &
                        elgtthn,        &
                        elgtdel,        &
                        einitmax,       &
                        overvolt,       &
                        elgtfdel,       &
                        ibrkd,          &
                        elgtfestopcg,   &
                        zgrnd, tgrnd,   &
                        ieint,          &
                        eint,           &
                        feint,          &
                        fecrit,         &
                        fefac,          &
                        iseed,          &
                        intfg,          &
                        intfgmpi,       &
                        mdel,           &
                        ibal,           &
                        cgthres,        &
                        dslight,        &
                        l2nde,          &
                        nohold,         &
                        maxpicks,       &
                        multipickmod,   &
                        pickmaskwidth,  &
                        cgfr,           &
                        soreps,         &
                        intsormpi,      &
                        lightintx,      &
                        lightintz,      &
                        ltgeta,         &
                        elec_on_time,   &
                        elec_ramp_time, &
                        fairweather,    &
                        igrid,          &
                        isa,            &
                        irand,          &
                        ixst0,          &
                        jyst0,          &
                        kzst0,          &
                        initwire,efracinitwire,wirelength, &
                        esctot,         &
                        bmg_cycletype,  &
                        bmg_maxiter,    &
                        bmg_sep,        &
                        bmg_cgsolve,    &
                        bmg_tol,        &
                        maxbzone,       &
                        chgthr,         &
                        elgt1,          &
                        elght1,         &
                        lightrad,lightextendmsz,       &
                        idoniconly,         &
                        energymethod,   &
                        crgfac,         &
                        iehlis, iehis, iesis, &
                        isctemopt,      &
                        nic_min_temp
                        
   NAMELIST /forcing/                   &
                        iforce,         &
                        iforcetyp,      &
                        nwfor,          &
                        iwforce,        &
                        iwshap,         &
                        ivshap,         &
                        iushap,         &
                        iqshap,         &
                        itshap,         &
                        isshap,isslcl,  &
                        wfmeso,         &
                        ufmeso,umaxmeso,&
                        xwfcen,         &
                        ywfcen,         &
                        zwfcen,         &
                        xwfrad,         &
                        ywfrad,         &
                        zwfrad,         &
                        xwfmov,         &
                        ywfmov,         &
                        zwfmov,         &
                        nssfor,         &
                        issforce,       &
                        ssfmeso,        &
                        iqrforce,       &
                        qrfmeso,        &
                        qhfmeso,        &
                        xssfcen,        &
                        yssfcen,        &
                        rssfrad,        &
                        tsstrt,         &
                        tsstop,         &
                        tsslow,         &
                        tsshigh,        &
                        twstrt,         &
                        twstop,         &
                        itopforce,      &
                        itopforcetype,  &
                        topforceq,      &
                        topforceN,      &
                        topforceX,      &
                        topforceY,      &
                        topforceR,      &
                        topforcedia,    &
                        topforcen0,     &
                        topforcealpha,  &
                        itopforcerainopt, topforcerainrate, &
                        igamrain,       &
                        chargeperparticle, &
                        low_level_cooling_flag, &
                        low_level_cooling_depth, &
                        low_level_cooling_rate, & 
                        low_level_cooling_start, &
                        low_level_cooling_end, &
                        iccnufforce,nccnfor,tccnstrt,tccnstop,xccnfcen,yccnfcen,zccnfcen, &
                        xccnfrad,yccnfrad,zccnfrad, &
                        ichgshap, nchgfor, chgratefor, xchgfcen, ychgfcen, zchgfcen,      &
                        xchgfrad, ychgfrad, zchgfrad, ichgforce, &
                        xchgcldfrad, ychgcldfrad, zchgcldfrad, xchgcldfcen, ychgcldfcen, zchgcldfcen, &
                        wnaylor, forcerandfac
                        
                        

   NAMELIST /balloon/                   &
                        numsnd,         &
                        brise,          &
                        istpos,         &
                        ijkst,          &
                        xyzst0,          &
                        tstsnd,         &
                        icont,          &
                        sndname        

   NAMELIST /trajectories/              &
                        itraj,          &
                        ntrajtype,      &
                        ixtrj1, ixtrj2, &
                        jytrj1, jytrj2, &
                        kztrj1, kztrj2, &
                        dxtraj,         &
                        dytraj,         &
                        dztraj,         &
                        time_traj1,     &
                        time_traj2,     &
                        riserate,       &
                        chgavex,        &
                        chgavez,        &
                        iverttraj

CONTAINS

!-------------------------------------------------------------------------------

 LOGICAL FUNCTION READ_NAMELIST(filename,namelist0)
 
  USE COMMASMPI_MODULE
  
  implicit none

  character(LEN=*) filename
  character(LEN=*) namelist0
  character(LEN=50) namelist

  integer ibeg, iend, ib, ie, lenv5dflds
  logical if_exist
  integer k
  integer istat,iunit
  character(len=1000) :: line
  integer, parameter :: error_unit = 0

  CALL STRING_LIMITS(filename, ibeg, iend)

  INQUIRE(file=filename(ibeg:iend), exist=if_exist)

  IF( .NOT. if_exist ) THEN
    IF ( my_rank == 0 ) THEN
      write(0,*) 'READ_NAMELIST:  ERROR - INPUT FILE:  ', filename(ibeg:iend), ' DOES NOT EXIST!!!'
      write(0,*) 'READ_NAMELIST:  ERROR - DOES NOT HAVE A FILE TO READ, EXITING'
    ENDIF
    STOP
  ENDIF
  
  iunit = 15

  open(15,file=filename(ibeg:iend),status='old',form='formatted')
  rewind(15)

  namelist = namelist0
  namelist = UCASE(namelist)
  CALL STRING_LIMITS(namelist, ibeg, iend)

  SELECT CASE( namelist(ibeg:iend) )

    CASE( 'ATTRIBUTE')
      read(15,NML=attributes)
      READ_NAMELIST = .true.
      lenv5dflds = 1
      
#ifndef MPI
      parallelio = .false.
#endif

      DO k = 1,size(v5dvarlist)
       IF ( v5dvarlist(k) .eq. ' ' ) EXIT
       IF ( k .eq. 1 ) v5dflds = '/'
         CALL STRING_LIMITS(v5dvarlist(k), ib, ie)
         write(v5dflds(lenv5dflds+1:lenv5dflds + ie - ib + 3),'(a)') v5dvarlist(k)(ib:ie)//'/'
         lenv5dflds = lenv5dflds + ie - ib + 2
      ENDDO
      
      IF ( rstprefix .eq. ' ' ) THEN
         filehead = prefix
      ELSE
         filehead = rstprefix
      ENDIF
         CALL STRING_LIMITS(filehead, ib, ie)
         lfilehead = ie - ib + 1
         membernumber = member

    CASE( 'GRIDN' )
      read(15,NML=gridn,iostat=istat)
      READ_NAMELIST = .true.
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading namelist: ',namelist, '-- not found or bad token'
      ENDIF

    CASE( 'MPI_PARAMS' )
      read(15,NML=mpi_params,iostat=istat)
      IF ( istat .ne. 0 ) THEN
        READ_NAMELIST = .false.
        IF ( number_of_processes > 1 ) write(0,*) my_rank,' : Problem reading mpi_params!'
      ELSE
        READ_NAMELIST = .true.
      ENDIF

    CASE( 'HOMOG_INIT' )
      read(15,NML=homog_init,iostat=istat)
      IF( nbble .gt. maxbub ) THEN
       IF ( my_rank == 0 ) THEN
         write(0,*) 'READ_NAMELIST:  Too many bubbles in HOMG_INIT!!!'
         write(0,*) 'READ_NAMELIST:  NBBLE > MAXBUB, increase maxbub to more than:  ', maxbub
       ENDIF
       STOP
      ELSE
        READ_NAMELIST = .true.
      ENDIF
      IF ( iforce .ne. 0  .and. my_rank == 0 ) write(0,*) 'iforce in homog_init is deprecated.  Please define in forcing'
      IF ( istat .ne. 0 .and. my_rank == 0  ) THEN
        write(0,*) 'Problem reading namelist: ',namelist,' -- not found or bad token'
      ENDIF

    CASE( 'FORCING' )
!        write(0,*) 'READ_NAMELIST:  Looking for Forcing namelist'
       CALL FORCE_INIT
      read(15,NML=forcing,iostat=istat)
      IF ( istat .eq. 0 ) THEN
!        write(0,*) 'READ_NAMELIST: Forcing namelist found'
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: Problem reading Forcing namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'BALLOON' )
!        write(0,*) 'READ_NAMELIST:  Looking for Forcing namelist'
      read(15,NML=balloon,iostat=istat)
      IF ( istat .eq. 0 ) THEN
!        write(0,*) 'READ_NAMELIST: balloon namelist found'
        READ_NAMELIST = .true.
      ELSE
!        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: Problem reading Balloon namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'TRAJECTORIES' )
      read(15,NML=trajectories,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSEIF ( istat == -1 ) THEN ! status of -1 means end of file was reached
         IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: Problem reading TRAJECTORIES namelist: Namelist not found'
         READ_NAMELIST = .false.
      ELSE
          IF ( my_rank == 0 ) THEN
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(error_unit,'(A)') 'Invalid line in namelist TRAJECTORIES: '//trim(line)
          ENDIF
          CALL COMMASMPI_ABORT()
          READ_NAMELIST = .false.
      ENDIF
      
    CASE( 'OUTPUT_OPTIONS' )
      read(15,NML=output_options,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSEIF ( istat == -1 ) THEN ! status of -1 means end of file was reached
         IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: Problem reading OUTPUT_OPTIONS namelist: Namelist not found'
         READ_NAMELIST = .false.
      ELSE
          IF ( my_rank == 0 ) THEN
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(error_unit,'(A)') 'Invalid line in namelist OUTPUT_OPTIONS: '//trim(line)
          ENDIF
          CALL COMMASMPI_ABORT()
          READ_NAMELIST = .false.
      ENDIF

    CASE( 'ENS_INIT' )
      read(15,NML=ens_init,iostat=istat)
      READ_NAMELIST = .true.
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading namelist: ',namelist
      ENDIF

    CASE( 'RUN' )
!     CALL RUN_ATT_MODULE_INIT
      read(15,NML=run,iostat=istat)
      READ_NAMELIST = .true.
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading namelist: ',namelist
      ENDIF

    CASE( 'LFO_PARAMS' )
      CALL LFO_INIT
      read(15,NML=lfo_params,iostat=istat)
      READ_NAMELIST = .true.
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading namelist: ',namelist
      ENDIF

    CASE( 'ICE10_PARAMS', 'MICRO_PARAMS' )
      CALL INDEX_MODULE_INIT
      read(15,NML=micro_params,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH MICRO_PARAMS namelist: not found or bad token'
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: Checking for ICE10_PARAMS, istat =', istat
        READ_NAMELIST = .false.

          IF ( istat /= -1 ) THEN ! istat = -1 means end of file; istat > 0 is some kind of bad read
          backspace(iunit)
          read(iunit,fmt='(A)') line
           write(error_unit,'(A)') &
          'Invalid line in namelist: '//trim(line)
          ENDIF

        rewind(15)
        read(15,NML=ice10_params,iostat=istat)
      
         IF ( istat /= 0 .and. istat /= -1 ) THEN
           write(error_unit,*) 'Error reading ice10_params. istat = ',istat
           backspace(iunit)
           read(iunit,fmt='(A)') line
           write(error_unit,'(A)') &
          'Invalid line in namelist: '//trim(line)
          CALL COMMASMPI_ABORT()
         ENDIF

      ENDIF
      rewind(15)
      CALL ELEC_MODULE_INIT
      IF ( my_rank == 0 ) write(0,*) 'read elec_module'
      read(15,NML=elec_params,iostat=istat)
      
      IF ( my_rank == 0 ) write(0,*) 'icorona = ',icorona
      IF ( istat .eq. 0 ) THEN
        IF ( ieopt == 0 .and. drake ) ieopt = 1 ! backwards compatibility with namelist that have 'drake'
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH ELEC_PARAMS namelist: not found or bad token'
        ipelec = 0
      ENDIF
      IF ( (microphys .eq. 'ICE10' .or. microphys .eq. 'ICE10E' .or. microphys(1:3) .eq. 'LFO'   &
            .or. microphys(1:7) .eq. 'WARMLFO' ) .and. &
            ipconc .gt. 0 ) ipconc = 0
!      IF ( microphys(1:1) .eq. 'Z' .or. microphys .eq. 'WARMZIEG' ) ipconc = Max(ipconc, 0)
      IF ( .not. microphys(1:6) .eq. 'ICE10E'   .and. .not. microphys(1:8) .eq. 'ICE10DME' .and.    & 
     &     .not. microphys(1:5) .eq. 'ZIEGE'    .and. .not. microphys(1:6) .eq. 'ZIEGHE' .and.  &
     &     .not. microphys(1:4) .eq. 'ZVDE'     .and. .not. microphys(1:5) .eq. 'ZVDME' .and.   &
     &     .not. microphys(1:5) .eq. 'ZVDHE'    .and. .not. microphys(1:6) .eq. 'ZVDHVE' .and.  &
     &     .not. microphys(1:6) .eq. 'ZVDMHE'   .and. .not. microphys(1:6) .eq. 'ZVDHME'.and.  &
     &     .not. microphys(1:4) .eq. 'TAKE' ) ipelec = Min(ipelec,0)

      IF ( ipconc < 5 .and. ( microphys(1:4) .eq. 'ZVDM' .or. microphys(1:5) .eq. 'ZVDHM' ) ) THEN
        write(0,*) 'Incompatible options: Must have ipconc >= 5 for ',microphys
        call commasmpi_abort()
      ENDIF
      IF ( ipconc <= 5 ) THEN
        ibinhmlr  = 0
        ibinhlmlr = 0
      ENDIF
    CASE( 'MY_PARAMS' )
      CALL INDEX_MODULE_INIT
      read(15,NML=my_params,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
!        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH MY_PARAMS namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'PARAMH' )
      CALL PARAM_MODULE_INIT
      read(15,NML=paramh,iostat=istat)
      READ_NAMELIST = .true.
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading namelist: ',namelist
      ENDIF

    CASE( 'CORRECT_ENS' )
      read(15,NML=correct_ens,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM with CORRECT_ENSEMBLE namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'TRMM' )
      CALL TRMM_MODULE_INIT
      read(15,NML=trmm,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH TRMM namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'ENKF_PARAM' )
      read(15,NML=enkf_param,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH ENKF_PARAM namelist: not found or bad token'
        READ_NAMELIST = .false.
      ENDIF

    CASE( 'BSS_INIT' )
      read(15,NML=bss_init,iostat=istat)
      IF ( istat .eq. 0 ) THEN
        READ_NAMELIST = .true.
      ELSE
        IF ( my_rank == 0 ) THEN
          write(0,*) 'READ_NAMELIST: PROBLEM WITH BSS_INIT namelist: not found or bad token'
        
          IF ( istat == -1 ) THEN
            write(0,*) 'Namelist not found'
          ELSE
             backspace(iunit)
             read(iunit,fmt='(A)') line
             write(error_unit,'(A)') &
            'Invalid line in namelist: '//trim(line)
          ENDIF
        
        ENDIF
        READ_NAMELIST = .false.
      ENDIF

    CASE DEFAULT
      write(0,*) 'READ_NAMELIST:  ERROR - UNKNOWN NAMELIST REQUESTED:  ', namelist(ibeg:iend)
      write(0,*) 'READ_NAMELIST:  ERROR - NAMELIST WAS NOT READ!!'
      READ_NAMELIST = .false.

   END SELECT

  close(15)

 END FUNCTION READ_NAMELIST

!-------------------------------------------------------------------------------

 SUBROUTINE WRITE_NAMELISTS(iunit)
 
  implicit none

  integer iunit
  integer istat

  write(iunit,'(a)') '! Namelist values for the SIMULATION'
  
  write(iunit,NML=attributes)
  
  write(iunit,'(a)')
  
  write(iunit,NML=gridn)
  
  write(iunit,NML=homog_init)
  
  write(iunit,NML=forcing,iostat=istat)
  
  write(iunit,NML=ens_init)
  
  write(iunit,NML=correct_ens)

  write(iunit,NML=enkf_param)

  write(iunit,NML=run)
  
  write(iunit,NML=output_options)
  
  write(iunit,NML=synthob_param)
  
  write(iunit,NML=lfo_params)
  
  write(iunit,NML=ice10_params,iostat=istat)

!  write(iunit,NML=micro_params,iostat=istat)
  
  write(iunit,NML=elec_params,iostat=istat)
  
  write(iunit,NML=paramh)

  write(iunit,NML=balloon)

      write(iunit,NML=enkf_param)
      
      write(iunit,NML=correct_ens)

      write(iunit,NML=bss_init)

      write(iunit,NML=trajectories)
      
      write(iunit,NML=nssl_mp_params)
      

 END SUBROUTINE WRITE_NAMELISTS

END MODULE NAMELIST_MODULE
