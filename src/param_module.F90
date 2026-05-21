 
!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        PARAM_MODULE       ////////////////////
!
!
!
!
!===========================================================================

MODULE PARAM_MODULE

!----------------------------------------------------------------------------------------------------------------------------------
!
!                                        ////////////////////  CONSTANTS   \\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!----------------------------------------------------------------------------------------------------------------------------------

      real,    parameter :: t00      = 300.0 
      real,    parameter :: p00      = 1.0e5 
      real,    parameter :: rd       = 287.04
      real,    parameter :: rw       = 461.50
      real,    parameter :: cp       = 1004.
      real,    parameter :: cpv      = 1870.
      real,    parameter :: Lav      = 2.5e6
      real,    parameter :: Cl       = 4190.
      real,    parameter :: epsilon  = 0.622
      real,    parameter :: epsilon2 = rd/rw
      real,    parameter :: cthres   = 1.0e-4
      real,    parameter :: g        = 9.806
      real,    parameter :: Cv       = 717.
      real               :: Cspd     = 300.      ! sound wave speed
      real,    parameter :: pii      = 3.1415926 
      real,    parameter :: gcp      = g / cp
      real,    parameter :: rcp      = rd / cp
      real,    parameter :: cvr      = cv / rd
      real,    parameter :: degtorad = pii/180.
      real,    parameter :: ae       = 1.21*6371.23   ! 6/5 earth radius in km
      real,    parameter :: adjrng   = 0.5
      real,    parameter :: adjazm   = 0.5
      real,    parameter :: rtd      = 1./degtorad

!----------------------------------------------------------------------------------------------------------------------------------
!
!                                        ////////////////////  PARAMS   \\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!----------------------------------------------------------------------------------------------------------------------------------

      integer            :: luno                 ! unit number for standard out 
      integer            :: lune                 ! unit number for standard error
      integer            :: ng       = 3         ! No. of ghost zones
      integer            :: ngv      = 3         ! No. of ghost zones in variable arrays
      integer            :: hsponge  = 0         ! No. grid points from boundary for lateral sponge/0=no sponge 
      real,    parameter :: alpha    = 0.1       ! parameter for vertically implicit small-step solver
      real               :: dxt      = 25.0      ! gravity wave speed for radiative boundary conditions
      real,    parameter :: lapse    = 0.0    
      real               :: kdiv0    = 0.05      ! divergence damping coefficient
      logical            :: vert_implicit = .true. ! flag for implicit or explicit vertical solver in small step
      real               :: zPBLhgt  = 1000.     ! Hgt of bnd layer in O'brien mixing profile
      real               :: hrayd_mag= 0.0025    ! Horizontal rayleigh damping magnitude
      real               :: coriol   = 0.0e-04   ! Coriolis parameter for f-plane approximation
      real               :: drag     = 1.0e-2    ! In original run, was semi-slip/rigid with drag=.0005
      real               :: rayd_hgt = 15000.    ! Height to begin rayleigh damping
      real               :: rayd_mag = 0.001     ! Vertical rayleigh damper magnitude
      real               :: damph    = 0.01      ! Filter coefficient

!  RK advection options
!  atypes1 (scalars) and atypem1 (momentum) are for the first steps
!  atypes2 (scalars) and atypem2 (momentum) are for the last step

      integer :: atypes1 = 13
      integer :: atypes2 = 1

      integer :: atype2qv = 1
      integer :: atype2th = 1
      integer :: atype2tke = 1

      integer :: atypem1 = 0
      integer :: atypem2 = 0
      
      integer :: iuvwadv = 1

      integer :: iwensmooth = 1 ! WENO smoothness weights option

      logical :: nocollapse = .true.  ! true to turn off stencil collapse at boundaries

!----------------------------------------------------------------------------------------------------------------------------------
! Ice parameterization particle densities and intercept parameters in MKS units
!  ---> moved to MICRO_MODULE namelist
!----------------------------------------------------------------------------------------------------------------------------------
! Boundary condition parameters ==> periodic and rigid laternal boundaries not implemented yet!

      integer :: bcx = 1         ! rigid = 0, open = 1 
      integer :: bcy = 1         ! rigid = 0, open = 1
      integer :: bcz = 0         ! free slip=  0, KW drag = 1

!----------------------------------------------------------------------------------------------------------------------------------
! Parameters controling mixing parameterization

      integer            :: mix_type = 1        ! [-1,0,1]=> [ahighk, SMAG, TKE]
      integer            :: len_type = 0        ! [0,1,2,3] => [volume,dz,BLR,constant (Lmax)]
      integer            :: rmbasekm = 0        ! [0,1]   => [dont rm, subtract kmbase method]
      integer            :: tke_type = 1        ! [1,2,3] => [km, sqrt(E), E]
      integer            :: iusebvfreq = 0      ! use Brunt-Vaisalla freq for TKE length (Deardorff 1980)
      real               :: Cm       = 0.21  
      real               :: Ce       = 0.70
      real               :: Pr       = 2.5
      real               :: Lmax     = 500. 
      real               :: ahighk   = 1000.
      real               :: km_thresh = 0.0   ! 0.01 set small values of KM to zero (TKE option only)
      real               :: kmgen_thresh = 0.0   ! 0.001 set small values of KM source term to zero (not sink, though)
      real               :: kmshear_thresh = 0.0  ! 0.001 set small values of shear term (shear_tke) to zero.
      real               :: kmbasefac = 1.0     ! factor to boost kmbase array

!----------------------------------------------------------------------------------------------------------------------------------
! Flags for model configuration options

      integer            :: RKSCHEME        = 3
      integer            :: rkstepping      = 1
      logical            :: MONOTONIC       = .true.
      integer            :: vert_adv_scheme = 5
      logical            :: lfwds           = .false.       ! .true.  => Use RK advection u,v,w,theta,tke, Crowley for the rest
                                                            ! .false. => Use RK advection for all advection
      logical            :: lfwdth          = .false.       ! .true. => use Crowley for theta and qv, also
      integer            :: icrwmp          = 1             ! Crowley: 1=1pass, 2=2-pass
      integer            :: icrwmn          = 1             ! Crowley monotonic filter: 0=off, 1=on
      integer            :: icrwmn1st       = 1             ! 1st order Crowley monotonic filter: 0=off, 1=on
      real               :: ainflo          = 0.5           ! set to >0 (e.g. 0.5) to turn on inflow nudging
                                                            !   for scalars that have base state of 0 (e.g., hydrometeor mixing ratios)
      real               :: ainfloqv        = 0.0           ! set to >0 (e.g. 0.5) to turn on inflow nudging
      real               :: ainflom         = 0.0           ! set to >0 (e.g. 0.5) to turn on inflow nudging for momentum (parallel component)
      real               :: zinflo_u_s      = 0             ! height to nudge U on south boundary
      real               :: zinflo_u_n      = 0             ! height to nudge U on north boundary
      real               :: zinflo_v_w      = 0             ! height to nudge V on west boundary
      real               :: zinflo_v_e      = 0             ! height to nudge V on east boundary
      real               :: zinflo_w_e      = 0             ! height to nudge W on ease boundary
      real               :: zinflo_w_w      = 0             ! height to nudge W on west boundary
      real               :: zinflo_w_s      = 0             ! height to nudge W on south boundary
      real               :: zinflo_w_n      = 0             ! height to nudge W on north boundary
      integer            :: ifilt           = 0             ! 1 = turn on 6th-order monotonic filter for scalars in conjunction with Crowley only
      integer            :: izero_opt       = 1             ! 1 = allow izero=1 in advection
      
      integer            :: isfcphys        = 0             ! flag for turning on surface physics
      logical            :: detrendpi       = .true.        ! flag for removing average pi perturbation   
      
      integer            :: rotunno_radiation = 0            ! Pseudo-radiation from Rotunno and Emanuel 1987 (JAS)
      
   CONTAINS
   
     SUBROUTINE PARAM_MODULE_INIT
     RETURN
     END SUBROUTINE PARAM_MODULE_INIT
     
     
END MODULE PARAM_MODULE

