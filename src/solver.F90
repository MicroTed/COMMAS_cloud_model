#define CHGELEC
#ifdef NOELEC
#undef CHGELEC
#endif
!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE SOLVER   ////////////////////
!     
!-----------------------------------------------------------------------------
  SUBROUTINE SOLVER(gd,                 &
                    u,  uinit ,         &            ! U,  UINIT
                    v,  vinit ,         &            ! V,  VINIT
                    w,  winit ,         &            ! W,  WINIT
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    rho,                &
                    kmbaserm,           &
                    s,  sinit,          &            ! S,  SINIT
                    uinit0, vinit0,     &
                    precip,             &            ! PRECIP
                    gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                    gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                    gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                 &            ! DT
                    ugrid, vgrid,       &            ! GRID MOTION
                    ugrid0, vgrid0,     &            ! GRID MOTION at t=0
                    lat, lon,           &            ! Latitude and Longitude
                    nsmall,             &            ! NSMALL
                    microphys,          &            ! microphysical scheme
                    nx, ny, nz, ns,     &            ! NX,NY,NZ,NS
                    io_flag,            &            ! IO_FLAG
                    st, ut, vt, wt, pt, &
                    fu, fv, fw, fp,     &
                    fs, kt, t0, ft,divv,kmt,kht,&
                    precip_old,         &
                    dbz, vzf, wz, xtra, &
                    elec, cion,         &
                    muz,                &
                    x_sw_loc,y_sw_loc,  &                                   
                    time, time_real, tstat, lstt, tstop, z1d4, onedoutput, this, &
                    run_file)

   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE MICRO_MODULE
   USE PARAM_MODULE
   USE INDEX_MODULE
   USE FORCE_MODULE
   USE TRAJ_MODULE
   USE INIT_MODULE, only: inhom,bogusvortex,timint,isfcl ! ,HURRFORCE
   USE TROPICAL_CYCLONE
#ifdef USEWRF
   use module_mp_nssl_2mom, only : nssl_2mom_driver, nssl_2mom_init, &
                   imurain_wrf => imurain, &
                   infall_wrf  => infall,  &
                   ccn_wrf     => ccn,     &
                  ! ccnuf_wrf   => ccnuf,    &
                   isedonly_wrf=> isedonly, &
                   dfrz_wrf    => dfrz,     &
                   iferwisventr_wrf => iferwisventr, &
                   izwisventr_wrf => izwisventr,  &
                   dmrauto_wrf    => dmrauto,     &
                   imaxdiaopt_wrf => imaxdiaopt,  &
                   ioldlimiter_wrf => ioldlimiter, &
                   ehw0_wrf => ehw0, &
                   ehlw0_wrf => ehlw0, &
                   eri0_wrf => eri0, &
                   ibfc_wrf => ibfc,  &
                   itype2_wrf => itype2, &
                   eii0hl_wrf => eii0hl, &
                   eii1hl_wrf => eii1hl, &
                   icenucopt_wrf => icenucopt, &
                   hlcnhdia_wrf => hlcnhdia, &
                   hlcnhqmin_wrf => hlcnhqmin, &
                   iehw_wrf => iehw, &
                   iehlw_wrf => iehlw, &
                   ibiggopt_wrf => ibiggopt, &
                   ibiggsnow_wrf => ibiggsnow, &
                   dhmn_wrf => dhmn, &
                   dhmx_wrf => dhmx, &
                   iacr_wrf => iacr, &
                   icracr_wrf => icracr, &
                   iacrsize_wrf => iacrsize, &
                   irenuc_wrf => irenuc,     &
                   alphah_wrf => alphah,     &
                   alphahl_wrf => alphahl,   &
                   ibinhmlr_wrf => ibinhmlr, &
                   ibinhlmlr_wrf => ibinhlmlr, &
                   iturbenhance_wrf => iturbenhance, &
                   ifrzs_wrf => ifrzs, &
                   ffrzs_wrf => ffrzs, &
                   hdnmn_wrf => hdnmn, &
                   hldnmn_wrf => hldnmn, &
                   isnowdens_wrf => isnowdens, &
                   nsplinter_wrf => nsplinter, &
                   xvdmx_wrf => xvdmx,   &
                   evapfac_wrf => evapfac, &
                   rainfallfac_wrf => rainfallfac, &
                   icefallfac_wrf => icefallfac, &
                   snowfallfac_wrf => snowfallfac, &
                   icefallopt_wrf => icefallopt, &
                   ihlcnh_wrf => ihlcnh, &
                   ndebug_wrf => ndebug, &
                   ac_opt_wrf => ac_opt, &
                   ac_wthresh_wrf => ac_wthresh, &
                   iqvsopt_wrf => iqvsopt
#endif
   USE module_mp_thompson, only: mp_gt_driver, thompson_init
   USE module_mp_thompson411, only: mp_gt_driver411, thompson_init411
   USE module_mp_thompson381, only: mp_gt_driver381, thompson_init381, thom_autoconv_fac_381 => thom_autoconv_fac
   USE MODULE_MP_MORR_TWO_MOMENT, only: MORR_TWO_MOMENT_INIT, MP_MORR_TWO_MOMENT, &
                    evapfac_morr => evapfac, &
                    icracr_morr => icracr, &
                    rssflg_morr => rssflg , &
                    sssflg_morr => sssflg , &
                    hssflg_morr => hssflg, &
                    ndcnst
                    
   USE COMMASMPI_MODULE
#ifdef TAKON
       use takcommon, only: capth,capw,accel,rdry,ratio,sheat !, &
!     &                      nqsat_tak=>nqsat,fqsat_tak=>fqsat,fqsati_tak=>fqsati, &
!     &                      tabqvs_tak=>tabqvs,tabqis_tak=>tabqis,dtabqvs_tak=>dtabqvs,dtabqis_tak=>dtabqis
#endif

#ifdef MPI
    USE mpi
#endif

!-----------------------------------------------------------------------------
   
   implicit none
 
#ifdef MPI
!  INCLUDE "mpif.h"
#endif

!-----------------------------------------------------------------------------
! GRID DEFINITIONS

   TYPE(GRID)         :: gd 

   character(LEN = *) :: microphys
   integer            :: mscheme,nscalar ! Needed for Milbrandt-Yau scheme
   integer            :: nx, ny, nz, ns
   integer            :: nsmall    
   real               :: dt
   real               :: ugrid, vgrid, ugrid0, vgrid0
   real               :: lat, lon
   logical            :: io_flag

   real unorm
   real vnorm

   TYPE(VARIABLE)     :: u, uinit
   TYPE(VARIABLE)     :: v, vinit
   TYPE(VARIABLE)     :: w, winit
   TYPE(VARIABLE)     :: pi, piinit
   TYPE(VARIABLE)     :: km, kminit, kmbaserm
   TYPE(VARIABLE)     :: rho
   TYPE(VARIABLE)     :: s(ns), sinit(2)     ! We are going to try and create arrays of the scalar variables needed here
   TYPE(VARIABLE)     :: uinit0, vinit0
   TYPE(VARIABLE)     :: precip(nprecip+neelec2d)
   TYPE(VARIABLE)     :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE)     :: dbz, vzf, wz
   TYPE(VARIABLE)     :: xtra(nxtra)
   TYPE(VARIABLE)     :: elec(neelec)
   TYPE(VARIABLE)     :: cion(2)
   TYPE(VARIABLE)     :: muz(4)


   real :: preciptmp(-ng+1:nx+ng,-ng+1:ny+ng,nprecip+neelec2d)
   real :: precip_old(-ng+1:nx+ng,-ng+1:ny+ng,nprecip)

   double precision :: preciptot(nprecip+neelec2d), preciptotall(nprecip+neelec2d)

   real :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)

   real :: ut (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: vt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: wt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: pt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: kt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: kmt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: kht(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

   real :: fu (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fv (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fw (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fp (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fs (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: ft (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: divv (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real, allocatable, save :: dpdt(:,:,:)
   real, allocatable, save :: tke_diss(:,:,:)
   real, allocatable, save :: khh(:,:,:),khv(:,:,:),kmh(:,:,:),kmv(:,:,:)
   
   real :: z1d4(nzend,4)
   
   integer :: onedoutput
   
   integer :: tstat, this
   logical :: lstt
   real    :: x_sw_loc,y_sw_loc
   character(LEN = 120)      :: run_file
      
! Local variables
      
   real, pointer :: ranarr(:,:)
   integer :: infileunit = -1
   logical, parameter :: debugsolver = .false.
   logical, parameter :: truetime = .true. ! turns on barriers to get timing of each section for the slowest process
   integer :: loop, ge, i, j, k, n, m, nstep, atype
   integer :: ntimestep = 0
   real    :: dx, dy, dz
   real    :: dx1, dy1, dz1, dyl, mlen
   real    :: dts, sdt
   real    :: den(-ng+1:nz+ng,2), rrp(-ng+1:nz+ng,2), rrm(-ng+1:nz+ng,2)
   double precision :: deninv(-ng+1:nz+ng,2)
   integer :: nsrk,nsfwd
      
   real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng), gtz(-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   integer :: nst = 1                                ! dummy timestep; used to flip xy|yx ordering in Crowley advection
   save nst
   integer :: iadiv
   real    :: dt1
   integer :: nsub
   real    :: tmp
   real    :: gamma
   real    :: cons1 = 0.0

   integer :: ix,jy,kz,ia, num2d
   integer :: i1,i2, j1,j2, k1,k2
   integer :: is, js, ks
   integer :: im1, jm1, ip1, jp1
   real    :: wmax, wmin, vmax, a
   integer :: time, tstop
   real    :: time_real

   real    :: q(2*lqmx), pii_total, t_total
   real    :: qtodbz

   real, allocatable, save :: sbase(:,:)

   real wzz

   logical, parameter :: filter = .false.
   integer, parameter :: ihole  = 2

! Flag to dump out dpdt stats for noise estimates

   logical, parameter :: print_noise = .false.

! microphysics flags for number of sub-steps and minimum dbz

      
   real, parameter    :: mindbz = 0.0
!   integer, parameter :: nphys = 1
!   iuvwadv is defd in param_module, read in namelist ! 1 = normal; 2 = box scheme; 3 = turns off wind update for pure scalar advection
   
   logical, parameter :: ltkeall = .true. ! .false. 
   logical, parameter :: ltkeon = .true. ! .false. 

   integer, save  :: ifirst = 0
   
   integer, save :: llen
   
   integer, save :: ib,ie,jb,je,kb,ke,ni,nj,nk
   
   integer :: nrain
   
   real, allocatable, save :: dz3d(:,:,:)
   
   real, allocatable, save :: pb(:), db(:)
   
   real, allocatable, save :: chgadvtemp(:,:,:,:)
   
   integer :: tflag ! flag to store turb mixing tendency (for supersaturation tendency)

      real, parameter :: ec= 1.602e-19, ecinv=1./ec  ! fundamental unit of charge

 ! RK-5 coefficients from Hu et al. (1996, J. Comp. Phys.)
    REAL(kind=8), PARAMETER :: ssp_coef(5) = (/0.197707993,0.237179241,0.3333116,0.5,1.0/)
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   integer :: westward_tag, eastward_tag
   integer :: northward_tag, southward_tag
   integer :: downward_tag, upward_tag

   logical :: debug_mpi = .false.

   integer       :: nampi, nb

! Takahashi microphysics
          real   (kind=8) :: a1,a2,a3,parta,pertrb                            &
     &                       ,term,term3,term4,term5,terma,work,worka          &
     &                       ,cap,farm,tempe,dummy,fbuf,pbuf

      real :: temq
      integer :: l

   
! Set advection options for rksteps:  
!  atypes1 (scalars) and atypem1 (momentum) are for the first steps
!  atypes2 (scalars) and atypem2 (momentum) are for the last step
!
!  IF( ATYPE .eq. 0 ) CALL ADVECT5()
!  IF( ATYPE .eq. 1 ) CALL ADVECT_WENO()
!  IF( ATYPE .eq. 2 ) CALL ADVECT_MONO(0.0,0.0,.false.)
!  IF( ATYPE .eq. 3 ) CALL ADVECT_MONO(0.0,0.0,.true.)
!  IF( ATYPE .eq. 4 ) CALL ADVECT_MONO(damph,0.0,.false.)
!  IF( ATYPE .eq. 5 ) CALL ADVECT_MONO(damph,0.0,.true.)
!  IF( ATYPE .eq. 6 ) CALL ADVECT_MONO(0.0,1.0,.true.)
!  IF( ATYPE .eq. 7 ) CALL ADVECT_MONO(damph,1.0,.false.)
!  IF( ATYPE .eq. 8 ) CALL ADVECT_MONO(damph,1.0,.true.)

!   integer, parameter :: atypes1 = 0
!   integer, parameter :: atypes2 = 5

!   integer, parameter :: atypem1 = 0
!   integer, parameter :: atypem2 = 1

! Vorticity function statements

   wzz(i,j,k,im1,jm1)=(v%flt3d(i,j,k)-v%flt3d(im1,j,k))*gx(4)%flt1d(i)-(u%flt3d(i,j,k)-u%flt3d(i,jm1,k))*gy(4)%flt1d(j)
      
!==================================================================================================================================

  IF ( debugsolver ) write(0,*) 'SOLVER3D:  BEGIN , my_rank = ',my_rank
  IF ( debugsolver ) write(0,*) ''
  IF ( debugsolver ) write(0,*) 'SOLVER3D:  NX/NY/NZ/NS/NG ', nx, ny, nz, ns, ng
  IF ( debugsolver ) write(0,*) 'RKSCHEME = ',RKSCHEME
       IF ( .not. ( RKSCHEME == 2 .or. RKSCHEME == 3 .or. RKSCHEME == 5 ) ) THEN
         write(0,*) 'Unsupported value of RKSCHEME = ',RKSCHEME
         call commasmpi_abort()
       ENDIF

!==================================================================================================================================
! COMPUTE locally need variables

  IF ( ntimestep == 0 ) THEN
    ntimestep = 1
  ELSE
    ntimestep = int(time/dt)
  ENDIF
#ifndef MPI
  IF ( print_noise .and. .not. allocated( dpdt ) ) THEN
    allocate( dpdt(nx-1,ny-1,nz-1) )
  ENDIF
#endif
  
  IF ( .not. allocated( tke_diss ) ) THEN
    IF ( iturbenhance > 0 ) THEN
      allocate( tke_diss(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
    ELSE
      allocate( tke_diss(1,1,1) )
    ENDIF
  ENDIF

  IF ( .not. allocated( khh ) ) THEN
     allocate( khh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
     allocate( khv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
     allocate( kmh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
     allocate( kmv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  ENDIF
  
  CALL cld_cpu('TIMESTEP')  

   CALL GET_VARIABLE (gd, 'DX', dx)
   CALL GET_VARIABLE (gd, 'DY', dy)
   CALL GET_VARIABLE (gd, 'DZ', dz)
   CALL GET_VARIABLE (gd,'RANARRAY2D',ranarr)

   nrain = nprecip+neelec2d
   
   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  DX/DY/DZ ', dx, dy, dz

   IF ( .not. allocated( sbase ) ) THEN

   allocate ( sbase(-ng+1:nz+ng,ns) )
   sbase(:,:) = 0.0

   DO n = 1,2
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
#else
     DO k = 1,nz
#endif
       sbase(k,n) = sinit(n)%flt1d(k)
     ENDDO
   ENDDO

   DO n = 3,ns
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
#else
     DO k = 1,nz
#endif
      sbase(k,n) = s(n)%base1d(k)
     ENDDO
   ENDDO

   ELSE

   sbase(:,:) = 0.0
   DO n = 1,2
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
       sbase(k,n) = sinit(n)%flt1d(k)
     ENDDO
   ENDDO

   DO n = 3,ns
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
      sbase(k,n) = s(n)%base1d(k)
     ENDDO
   ENDDO


   ENDIF


   preciptmp(:,:,:) = 0.0
   DO k = 1,nprecip+neelec2d
#ifdef MPI
   ixb = 1
   ixe = itile
   if(ixend .eq. nxend) ixe = ixend-ixbeg
   
   jyb = 1
   jye = jtile
   if(jyend .eq. nyend) jye = jyend-jybeg
   
   do j = jyb,jye
    do i = ixb, ixe
#else
   DO j = 1,ny-1
     DO i = 1,nx-1
#endif
       preciptmp(i,j,k) = precip(k)%flt2d(i,j)
     ENDDO
   ENDDO 
   ENDDO

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  LFWDS '

! FOR RK/Forward:

   IF ( lfwds ) THEN
!     IF ( lfwdth ) THEN
!       nsrk  = 1  ! leave TH on RK3
!       nsfwd = ns
!     ELSE
       nsrk  = 2
       nsfwd = ns
!     ENDIF
   ELSE
     nsrk = ns
     nsfwd = 0
   ENDIF

     IF ( .not. allocated( pb ) ) THEN
       allocate( pb(-ng+1:nz+ng) )
       allocate( db(-ng+1:nz+ng) )
     ENDIF
      kzb = -1
      kze = ktile+2
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO kz=kzb,kze

       db(kz)= 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sbase(kz,lt))
       pb(kz)= 1.0e5*piinit%flt1d(kz)**3.509
        
      ENDDO

!! Constants needed
   IF ( ifirst == 0 ) THEN ! {
 
 
#ifdef CHGELEC
   IF ( inetchargetend > 0 ) THEN
     allocate ( chgadvtemp(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,2) )
   ELSE
     allocate ( chgadvtemp(1,1,1,1) )
   ENDIF
#else
   allocate ( chgadvtemp(1,1,1,1) )
#endif
      
      IF ( microphys(1:4) == 'THOM' ) THEN 
!        write(0,*) 'solver: Calling thomson_init'
        IF ( microphys(1:7) == 'THOM411' ) THEN
          CALL thompson_init411
        ELSEIF ( microphys(1:7) == 'THOM381' ) THEN
          thom_autoconv_fac_381 = thom_autoconv_fac
          CALL thompson_init381
        ELSE
          CALL thompson_init
        ENDIF
!        write(0,*) 'solver: done with thomson_init'

        IF ( .not. allocated( dz3d ) ) THEN
          allocate( dz3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
          
          DO kz=kzb,kze
            dz3d(:,:,kz) = 1./gz(3)%flt1d(kz)
          ENDDO
          
        ENDIF
        
          ib = -ng+1
          jb = -ng+1
          kb = -ng+1
          
          ie = nx+ng
          je = ny+ng
          ke = nz+ng
          
          ni = ixend-ixbeg+1
          IF ( ixend == nxend ) ni = ixend-ixbeg
          nj = jyend-jybeg+1
          IF ( jyend == nyend ) nj = jyend-jybeg
          IF ( ny == 2 ) nj = 1
          nk = nzend-1
        
      ENDIF

      IF ( microphys(1:4) == 'MORR' ) THEN 
!        write(0,*) 'solver: Calling MORR_TWO_MOMENT_INIT'

        evapfac_morr = evapfac
        icracr_morr = icracr
        rssflg_morr = rssflg
        sssflg_morr = sssflg
        hssflg_morr = hssflg
        
        CALL MORR_TWO_MOMENT_INIT(morr_rimed_ice,morr_dnr_max)
        
        ndcnst = morr_droplet_conc
        
!        write(0,*) 'solver: done with MORR_TWO_MOMENT_INIT'

        IF ( .not. allocated( dz3d ) ) THEN
          allocate( dz3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
          
          DO kz=kzb,kze
            dz3d(:,:,kz) = 1./gz(3)%flt1d(kz)
          ENDDO
          
        ENDIF
        
              !         CS = RHOSN*PI/6.
              !         DS = 3.
              !          CONS1=GAMMA(1.+DS)*CS
          cons1 = gamma(1. + 3.0)*100.*pii/6.
        
          ib = -ng+1
          jb = -ng+1
          kb = -ng+1
          
          ie = nx+ng
          je = ny+ng
          ke = nz+ng
          
          ni = ixend-ixbeg+1
          IF ( ixend == nxend ) ni = ixend-ixbeg
          nj = jyend-jybeg+1
          IF ( jyend == nyend ) nj = jyend-jybeg
          IF ( ny == 2 ) nj = 1
          nk = nzend-1
        
      ENDIF

     llen = index(microphys,' ')-1
     llen = len( trim(microphys) )
     
#ifdef USEWRF
     
     infileunit = -1
     
     IF ( iuvwadv >= 3 ) ntimestep = 1  ! set for first time step even if a restart
     
     ccn_wrf     = ccn
    ! IF ( ccnuf > 0.0 ) ccnuf_wrf  = ccnuf
     imurain_wrf = imurain
     
     IF ( .not. wrfdefaults ) THEN
     infall_wrf  = infall
     isedonly_wrf = isedonly
     dfrz_wrf    = dfrz
     iferwisventr_wrf = iferwisventr
     izwisventr_wrf = izwisventr
     dmrauto_wrf    = dmrauto
     imaxdiaopt_wrf = imaxdiaopt
     ioldlimiter_wrf = ioldlimiter
     ehw0_wrf = ehw0
     ehlw0_wrf = ehlw0
     eri0_wrf = eri0
     ibfc_wrf = ibfc
     itype2_wrf = itype2
     eii0hl_wrf = eii0hl
     eii1hl_wrf = eii1hl
     icenucopt_wrf = icenucopt
     hlcnhdia_wrf = hlcnhdia
     hlcnhqmin_wrf = hlcnhqmin
     iehw_wrf = iehw
     iehlw_wrf = iehlw
     ibiggopt_wrf = ibiggopt
     ibiggsnow_wrf = ibiggsnow
     dhmn_wrf = dhmn
     dhmx_wrf = dhmx
     iacr_wrf = iacr
     icracr_wrf = icracr
     iacrsize_wrf = iacrsize
     irenuc_wrf = irenuc
     alphah_wrf = alphah
     alphahl_wrf = alphahl
     ibinhmlr_wrf = ibinhmlr
     ibinhlmlr_wrf = ibinhlmlr
     iturbenhance_wrf = iturbenhance
     ifrzs_wrf = ifrzs
     ffrzs_wrf = ffrzs
     hdnmn_wrf = hdnmn
     hldnmn_wrf = hldnmn
     isnowdens_wrf = isnowdens
     nsplinter_wrf = nsplinter
     xvdmx_wrf = xvdmx
     evapfac_wrf = evapfac
     icefallopt_wrf = icefallopt
     rainfallfac_wrf = rainfallfac
     icefallfac_wrf = icefallfac
     snowfallfac_wrf = snowfallfac
     ihlcnh_wrf = ihlcnh
     ndebug_wrf = ndebug
     ac_opt_wrf = ac_opt
     ac_wthresh_wrf = ac_wthresh
     iqvsopt_wrf = iqvsopt
     
     ELSE

      infileunit = 29
      open(infileunit,file=run_file,status='old',form='formatted',action='read')
 
     ENDIF

        IF ( microphys(llen:llen) == 'F' ) THEN
!         write(0,*) 'SOLVER: set up for WRF micro'
          IF ( microphys(1:4) == 'ZVDH' ) THEN
            IF ( ipconc .eq. 5 .or. ipconc .eq. 8 ) THEN
!            write(0,*) 'lccnuf = ',lccnuf, (lccnuf/Max(1,lccnuf))
              IF ( microphys(1:5) == 'ZVDHM' ) THEN
              CALL nssl_2mom_init(ipctmp=ipconc,mixphase=1,nssl_hail_on=.true., nssl_ccn_on=.true.,infileunit=infileunit, &
                                  nssl_ufccn=(lccnuf/Max(1,lccnuf)))
              ELSE
!              CALL nssl_2mom_init(ipctmp=ipconc,mixphase=0,ihvol=1,infileunit=infileunit, &
!                                  nssl_ufccn=(lccnuf/Max(1,lccnuf)))
!              write(0,*) 'solver: ipconc = ',ipconc
              CALL nssl_2mom_init(ipctmp=ipconc,mixphase=0,nssl_hail_on=.true.,infileunit=infileunit, &
                                  nssl_ufccn=(lccnuf/Max(1,lccnuf)))
              ENDIF
            ELSE
             write(0,*) 'Must have ipconc = 5 or 8 for WRF version of ZVDH, ipconc =',ipconc
             call commasmpi_abort()
            ENDIF
          
          ELSEIF ( microphys(1:3) == 'ZVD' ) THEN
            IF ( ipconc .eq. 5 .or. ipconc .eq. 8 ) THEN
             IF ( ffrzs >= 0.999 ) THEN ! need new option for ihvol=-2 to combine ice and snow
             CALL nssl_2mom_init(ipctmp=ipconc,mixphase=0,nssl_hail_on=.false.,nssl_icecrystals_on=.false., &
                                 infileunit=infileunit, nssl_ufccn=(lccnuf/Max(1,lccnuf)))
             ELSE
             CALL nssl_2mom_init(ipctmp=ipconc,mixphase=0,nssl_hail_on=.false.,infileunit=infileunit, &
                                  nssl_ufccn=(lccnuf/Max(1,lccnuf)))
             ENDIF
            ELSE
             write(0,*) 'Must have ipconc = 5 or 8 for WRF version of ZVD'
             call commasmpi_abort()
            ENDIF
          
          ELSE
             write(0,*) 'Bad options for WRF version of ZVD'
             call commasmpi_abort()
          
          ENDIF
          
          IF ( wrfdefaults ) THEN
            IF (  irenuc /= irenuc_wrf ) THEN
              IF ( my_rank == 0 ) write(0,*) 'For wrfdefaults=true, must set same irenuc in both namelists!'
               call commasmpi_abort()
            ENDIF
          ENDIF

        IF ( .not. allocated( dz3d ) ) THEN
          allocate( dz3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
          
          DO kz=kzb,kze
            dz3d(:,:,kz) = 1./gz(3)%flt1d(kz)
          ENDDO
          
        ENDIF
        
          ib = -ng+1
          jb = -ng+1
          kb = -ng+1
          
          ie = nx+ng
          je = ny+ng
          ke = nz+ng
          
          ni = ixend-ixbeg+1
          IF ( ixend == nxend ) ni = ixend-ixbeg
          nj = jyend-jybeg+1
          IF ( jyend == nyend ) nj = jyend-jybeg
          IF ( ny == 2 ) nj = 1
          nk = nzend-1
        
        ENDIF
#endif /* USEWRF */

#ifdef MPI

    IF ( number_of_processes .gt. 1 ) THEN

!
! Here we communicate the 1-d grid distances and intervals (x,dx,y,dy) on the first step
!

      nampi = 4

        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF

      CALL cld_cpu('MPICOM-SOLVER')

      IF ( nproci .gt. 1 ) THEN

     ! copy gx(1-4) into gxt to save separate communications on each gx/gy variable
      ! gxt(:,1) = xc; 2=xe; 3=dxc (1/dx); 4 = dxe (1/dx)
       DO i = 1,4
        gxt(:,i) = gx(i)%flt1d(:)
       ENDDO

       westward_tag = 10001
       CALL sendrecv_westward(nx,1,1,ng,0,0,ng,nampi,  &
         w_proc(my_rank),e_proc(my_rank),westward_tag, gxt(-ng+1,1) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,1,1,ng,0,0,ng,nampi,  &
         w_proc(my_rank),e_proc(my_rank),eastward_tag, gxt(-ng+1,1) )
   
       IF ( myproci < nproci .or. bcx /= 2 ) THEN
         DO i = 1,4
           gx(i)%flt1d(:) = gxt(:,i)
         ENDDO
       ELSEIF ( myproci == nproci .and. bcx == 2 ) THEN

       ! for periodic case, do not overwrite the high end of the eastern domain. Otherwise xe(nxend) becomes 0.
         DO i = 1,2
           gx(i)%flt1d(-ng+1:1) = gxt(-ng+1:1,i)
         ENDDO
         DO i = 3,4
           gx(i)%flt1d(:) = gxt(:,i)
         ENDDO
       
       ENDIF

      ENDIF
    
      IF ( nprocj .gt. 1 ) THEN
       
     ! copy gy(1-4) into gyt to save separte communications on each gx/gy variable
      
       DO i = 1,4
        gyt(:,i) = gy(i)%flt1d(:)
       ENDDO

       southward_tag = 10003
       CALL sendrecv_southward(1,ny,1,0,ng,0,ng,nampi,  &
         n_proc(my_rank),s_proc(my_rank),southward_tag, gyt(-ng+1,1) )

       northward_tag = 10004
       CALL sendrecv_northward(1,ny,1,0,ng,0,ng,nampi,  &
         n_proc(my_rank),s_proc(my_rank),northward_tag, gyt(-ng+1,1) )
      
       DO i = 1,4
         gy(i)%flt1d(:) = gyt(:,i)
       ENDDO
   
      ENDIF

      IF ( nprock .gt. 1 ) THEN

     ! copy gx(1-4) into gxt to save separate communications on each gx/gy variable
      
       DO i = 1,4
        gzt(:,i) = gz(i)%flt1d(:)
       ENDDO
       
       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,nampi,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, gzt(-ng+1,1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,nampi,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, gzt(-ng+1,1) )
   
       DO i = 1,4
         gz(i)%flt1d(:) = gzt(:,i)
       ENDDO

       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,ns,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, sbase(-ng+1,1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,ns,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, sbase(-ng+1,1) )

       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, piinit%flt1d(-ng+1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, piinit%flt1d(-ng+1) )

       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, uinit%flt1d(-ng+1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, uinit%flt1d(-ng+1) )

       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, vinit%flt1d(-ng+1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,1,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, vinit%flt1d(-ng+1) )

      ENDIF

        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF

      CALL cld_cpu('MPICOM-SOLVER')
    
    ENDIF ! number_of_processes .gt. 1 
#endif

    IF ( myprock == nprock ) THEN

    DO i = 3,3
       gz(i)%flt1d(nz:nz+ng) = gz(i)%flt1d(nz-1)
    ENDDO

    DO i = 4,4
       gz(i)%flt1d(nz+1:nz+ng) = gz(i)%flt1d(nz)
    ENDDO

    DO i = 1,2
    DO n=1,ng
       gz(i)%flt1d(nz+n) = gz(i)%flt1d(nz) + n/gz(3)%flt1d(nz-1)
    ENDDO
    ENDDO

    ENDIF

#ifdef TAKON
! set lookup table for ice/water saturation

!
! Build lookup table for saturation mixing ratio (Soong and Ogura 73)
!
!      cai = 21.87455
!      caw = 17.2693882
!      cbi = 7.66
!      cbw = 35.86

      do l = 1,nqsat
      temq = 163.15 + (l-1)*fqsat
      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
      dtabqvs(l) = ((-caw*(-273.15 + temq))/(temq - cbw)**2 + & 
     &                 caw/(temq - cbw))*tabqvs(l)
      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
      dtabqis(l) = ((-cai*(-273.15 + temq))/(temq - cbi)**2 + & 
     &                 cai/(temq - cbi))*tabqis(l)
      end do

#endif
    
   ENDIF ! } ifirst
      
   gtx(:) = gx(3)%flt1d(:)*dx
   gty(:) = gy(3)%flt1d(:)*dy
   gtz(:) = gz(3)%flt1d(:)*dz
       
   DO i = 1,4
     gxt(:,i) = gx(i)%flt1d(:)
     gyt(:,i) = gy(i)%flt1d(:)
     gzt(:,i) = gz(i)%flt1d(:)
   ENDDO


   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  set den'

   ge  = -ng + 1

! Density variables (den(:,1) = scalar-pt, den(:,2) = w-pt)


#ifdef MPI
   kzb = -ng+1
   kze = ktile+ng
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg
     
   do k = kzb,kze
    den(k,1) = 1.0e5*piinit%flt1d(k)**2.509/(rd*sbase(k,1) * (1.0 + 0.61*sbase(k,2)) )
!   IF ( my_rank == 0 ) write(luno,*) my_rank, 'k,den1,ub,vb = ',k,den(k,1),uinit%flt1d(k),vinit%flt1d(k)
   enddo

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  set den2'

   if (kzbeg .eq. nzbeg) den(-ng+1:0,1)      = den(1,1)
   if (kzend .eq. nzend) den(nz:nz+ng,1)     = den(nz-1,1)
   if (kzbeg .eq. nzbeg) den(-ng+1:1,2)      = den(-ng+1:1,1)
   if (kzend .eq. nzend) den(nz:nz+ng,2)     = den(nz:nz+ng,1)

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  set den3'

   kzb = -ng+2
   kze = ktile+ng
   if (kzbeg .eq. nzbeg) kzb = 2
   if (kzend .eq. nzend) kze = kzend-kzbeg

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  set den4: kzb,kze =',kzb,kze
     
   do k = kzb,kze
    den(k,2) = 0.5*(den(k,1) + den(k-1,1))
   IF ( debugsolver ) write(0,*) my_rank, 'k,den1 = ',k,den(k,1)
   enddo
   
#else
   den(1:nz-1,1) = 1.0e5*piinit%flt1d(1:nz-1)**2.509/(rd*sbase(1:nz-1,1) * (1.0 + 0.61*sbase(1:nz-1,2)))
   den(0,1)      = den(1,1)
   den(nz,1)     = den(nz-1,1)
   den(0:1,2)    = den(0:1,1)
   den(nz,2)     = den(nz,1)
   den(2:nz-1,2) = 0.5*(den(2:nz-1,1) + den(1:nz-2,1))
#endif

! Reciprocal density variables for flux calculations

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  set rrp'

   rrp(:,:)      = 1.0  ! init to 1.0 for all undefined reciprocals
   rrm(:,:)      = 1.0  ! init to 1.0 for all undefined reciprocals

#ifdef MPI
   kzb = -ng+2
   kze = ktile+ng-1
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg

   do k = kzb, kze
    rrp(k,1) = den(k+1,2) / den(k,1)    ! scalar correction top
    rrm(k,1) = den(k,2)   / den(k,1)    ! scalar correction bottom
    deninv(k,1) = 1.d0/den(k,1)
!    write(luno,*) my_rank ,'k,rrp,rrm1 = ',k,rrp(k,1),rrm(k,1)
   enddo
   
   kzb = -ng+2
   kze = ktile+ng
   if (kzbeg .eq. nzbeg) kzb = -ng+2
   if (kzend .eq. nzend) kze = kzend-kzbeg+ng

   do k = kzb, kze
    rrp(k,2) = den(k,1)   / den(k,2)   ! w      correction top
    rrm(k,2) = den(k-1,1) / den(k,2)   ! w      correction bottom
   enddo
#else
   rrp(1:nz-1,1) = den(2:nz,  2) / den(1:nz-1,1)   ! scalar correction top
   rrm(1:nz-1,1) = den(1:nz-1,2) / den(1:nz-1,1)   ! scalar correction bottom
   rrp(2:nz-1,2) = den(2:nz-1,1) / den(2:nz-1,2)   ! w      correction top
   rrm(2:nz-1,2) = den(1:nz-2,1) / den(2:nz-1,2)   ! w      correction bottom
   deninv(1:nz-1,1) = 1.d0/den(1:nz-1,1)
#endif

  CALL cld_cpu('TIMESTEP')


!=================================================================================================================================
! MPI: FILL GHOST ZONES ON FIRST STEP OR RESTART



   IF ( ifirst == 0 ) THEN
     ifirst = 1
#ifdef MPI
!   IF ( time .gt. dt ) THEN ! assume this is a restart and communicate ghostzones

    IF ( number_of_processes .gt. 1 ) THEN

    CALL cld_cpu('MPICOM-1st')

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  first 3D comms'

    nampi = s(1)%index - u%index + ns
    
   IF ( nproci > 1 ) THEN

   westward_tag = 101
   CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),westward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   eastward_tag = 102
   CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),eastward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))
   
   ENDIF

   IF ( nprocj > 1 ) THEN

   southward_tag = 103
   CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),southward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   northward_tag = 104
   CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),northward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   ENDIF
   
   IF ( nprock > 1 ) THEN

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  first 3D comms start k'

   downward_tag = 101
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  first 3D comms halfway k'

   upward_tag = 102
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  first 3D comms done k'
   
   ENDIF

    CALL cld_cpu('MPICOM-1st')
    
    ENDIF
   
!    ENDIF
#endif
   ENDIF ! ifirst


!=================================================================================================================================

!==================================================================================================================================
! COPY data into tmp arrays for RK time stepping

  CALL cld_cpu('TIMESTEP-CP1') 
  CALL cld_cpu('TIMESTEP')

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  copy arrays'
   
#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1

    jyb = -ng+1
    jye = jtile+ng

    ixb = -ng+1
    ixe = itile+ng

    do k = kzb,kze
     do j = jyb,jye
      do i = ixb,ixe
       ut(i,j,k) = u%flt3d(i,j,k)
       vt(i,j,k) = v%flt3d(i,j,k)
       wt(i,j,k) = w%flt3d(i,j,k)
       pt(i,j,k) = pi%flt3d(i,j,k)
       kt(i,j,k) = km%flt3d(i,j,k)
      enddo
     enddo
    enddo


#else
    DO k = 1,nz
    ut(1:nx,1:ny,k) = u%flt3d(1:nx,1:ny,k) 
    vt(1:nx,1:ny,k) = v%flt3d(1:nx,1:ny,k) 
    wt(1:nx,1:ny,k) = w%flt3d(1:nx,1:ny,k) 
    pt(1:nx,1:ny,k) = pi%flt3d(1:nx,1:ny,k) 
    kt(1:nx,1:ny,k) = km%flt3d(1:nx,1:ny,k) 
   ENDDO
#endif

! check for variables to convert to mixing ratio for advection/mixing
   DO n = 1,ns

   IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  check scaling for scalar ',n

     IF ( s(n)%dyntype .eq. 0 .and. denscale >= 1 ) THEN
#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg+ng

    jyb = -ng+1
    jye = jtile+ng
    if (jyend .eq. nyend) jye = jyend-jybeg+ng

    ixb = -ng+1
    ixe = itile+ng
    if (ixend .eq. nxend) ixe = ixend-ixbeg+ng

    do k = kzb,kze
     do j = jyb,jye
      do i = ixb,ixe
       s(n)%flt3d(i,j,k) = s(n)%flt3d(i,j,k)/den(k,1)
       st(i,j,k,n) = s(n)%flt3d(i,j,k)
      enddo
     enddo
    enddo
#else
      DO k = 1,nz-1
      s(n)%flt3d(1:nx-1,1:ny-1,k) = s(n)%flt3d(1:nx-1,1:ny-1,k)/den(k,1)
      st(1:nx-1,1:ny-1,k,n) = s(n)%flt3d(1:nx-1,1:ny-1,k)
      ENDDO
#endif

    ELSE
    
#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng

    jyb = -ng+1
    jye = jtile+ng
    if (jyend .eq. nyend) jye = jyend-jybeg+1+ng

    ixb = -ng+1
    ixe = itile+ng
    if (ixend .eq. nxend) ixe = ixend-ixbeg+1+ng

    do k = kzb,kze
     do j = jyb,jye
      do i = ixb,ixe
       st(i,j,k,n) = s(n)%flt3d(i,j,k)
      enddo
     enddo
    enddo
#else
    st(1:nx-1,1:ny-1,1:nz-1,n) = s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1)
#endif

   ENDIF
   
   ENDDO  !! n

  IF ( allocated(chgadvtemp) ) chgadvtemp(:,:,:,:) = 0.0
  
  CALL cld_cpu('TIMESTEP')
  CALL cld_cpu('TIMESTEP-CP1') 

!==================================================================================================================================
! RK LOOP

   DO loop = 1,RKSCHEME

#ifdef MPI
    if (debug_mpi) write(0,*) my_rank,"SOLVER: START OF RK LOOP=",loop
#endif

     IF ( rkstepping == 1 ) THEN ! 1/3, 1/2, 1
       IF ( RKSCHEME == 2 .or. RKSCHEME == 3 ) THEN
         nstep = NInt( nsmall / float((RKSCHEME+1)-loop) + 0.1 )
         nstep = Max(1,nstep)
         sdt   = dt     / float((RKSCHEME+1)-loop)
       ELSEIF ( RKSCHEME == 5 ) THEN
         sdt  = dt * ssp_coef(loop)
         nstep = NInt( ssp_coef(loop)*nsmall + 0.1 )
         nstep = Max(1,nstep)
         
       ELSE
!         write(0,*) 'Improper value of RKSCHEME!'
       ENDIF
     ELSE ! IF ( rkstepping == 2 ) THEN ! 1/2, 1/2, 1
       IF( loop < RKSCHEME ) THEN
         nstep = nsmall / 2
         sdt   = dt / 2.0
       ELSE
         nstep = nsmall
         sdt   = dt 
       ENDIF
     ENDIF
     dts   = sdt / float(nstep)

!   DO k = 1,nz
!    st(1:nx,1:ny,k,1) = sinit(1)%flt1d(k) 
!   ENDDO

  IF ( bcx .eq. 2 .and. nproci == 1 ) THEN    ! w/e-periodic

#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg

    do k = kzb,kze
     do j = jyb,jye
       ut(nx,j,k) = ut(1,j,k)
     enddo
    enddo
#else
    ut(nx,1:ny-1,1:nz-1) = ut(1,1:ny-1,1:nz-1)
#endif

   DO n = 1,ng

#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg+1

    do k = kzb,kze
     do j = jyb,jye

       ut( nx+n,j,k) = ut(n+1 ,j,k)
       ut( 1-n ,j,k) = ut(nx-n,j,k) 

       wt( 1-n   ,j,k) = wt(nx-n,j,k) 
       wt( nx+n-1,j,k) = wt(n   ,j,k)
       vt( 1-n   ,j,k) = vt(nx-n,j,k) 
       vt( nx+n-1,j,k) = vt(n   ,j,k)
     enddo
    enddo
#else
    ut( nx+n,1:ny,1:nz) = ut(n+1 ,1:ny,1:nz)
    ut( 1-n ,1:ny,1:nz) = ut(nx-n,1:ny,1:nz) 

    wt( 1-n   ,1:ny,1:nz) = wt(nx-n,1:ny,1:nz) 
    wt( nx+n-1,1:ny,1:nz) = wt(n   ,1:ny,1:nz)
    vt( 1-n   ,1:ny,1:nz) = vt(nx-n,1:ny,1:nz) 
    vt( nx+n-1,1:ny,1:nz) = vt(n   ,1:ny,1:nz)
#endif
   ENDDO
   
#ifndef MPI
   gxt(nx,3) = gxt(1,3)
   gxt(nx,1) = gxt(1,1)
   
   gx(3)%flt1d(nx) = gx(3)%flt1d(1)
   gx(1)%flt1d(nx) = gx(1)%flt1d(1)
#endif
   
  ENDIF
  
!  IF ( bcx /= 2 ) THEN
  ! fill in outer ghost zones. Used for electricity extended boundaries
  ! Use separate IF statements since it can happen that nproci=1
    IF ( myproci == nproci ) THEN
     gxt(nx:nx+ng,3) = gxt(nx-1,3)
     gtx(nx:nx+ng) = gtx(nx-1)
     gx(3)%flt1d(nx:nx+ng) = gx(3)%flt1d(nx-1)
    ENDIF
    IF ( myproci == 1 ) THEN
     gxt(-ng+1:0,3) = gxt(1,3)
     gtx(-ng+1:0) = gtx(1)
     gx(3)%flt1d(-ng+1:0) = gx(3)%flt1d(1)
    
    ENDIF
!  ENDIF


  IF ( bcy .eq. 2 .and. nprocj == 1 ) THEN    ! n/s-periodic
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    do k = kzb,kze
     do i = ixb,ixe
      vt(i,ny,k) = vt(i,1,k)
     enddo
    enddo
#else
    vt(1:nx-1,ny,1:nz-1) = vt(1:nx-1,1,1:nz-1)
#endif

   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg+1

    do k = kzb,kze
     do i = ixb,ixe
      vt( i, ny+n ,k)   = vt(i,n+1 ,k) 
      vt( i,  1-n ,k)   = vt(i,ny-n,k) 
      wt( i,  1-n ,k)   = wt(i,ny-n,k) 
      wt( i,ny+n-1,k)   = wt(i,  n ,k)   
      ut( i,  1-n ,k)   = ut(i,ny-n,k) 
      ut( i,ny+n-1,k)   = ut(i,  n ,k)
     enddo
    enddo
#else
      vt(1:nx,ny+n ,1:nz)   = vt(1:nx,n+1 ,1:nz) 
      vt(1:nx, 1-n ,1:nz)   = vt(1:nx,ny-n,1:nz) 
      wt( 1:nx,1-n ,1:nz)   = wt(1:nx,ny-n,1:nz) 
      wt( 1:nx,ny+n-1,1:nz) = wt(1:nx,n ,1:nz)   
      ut( 1:nx,1-n ,1:nz)   = ut(1:nx,ny-n,1:nz) 
      ut( 1:nx,ny+n-1,1:nz) = ut(1:nx,n ,1:nz)
#endif
   ENDDO

#ifndef MPI
   gyt(ny,3) = gyt(1,3)
   gyt(ny,1) = gyt(1,1)
   
   gy(3)%flt1d(ny) = gy(3)%flt1d(1)
   gy(1)%flt1d(ny) = gy(1)%flt1d(1)
#endif

  ENDIF

!  IF ( bcy /= 2 ) THEN
  ! fill in outer ghost zones. Used for electricity extended boundaries
  ! Use separate IF statements since it can happen that nprocj=1
    IF ( myprocj == nprocj ) THEN
      gyt(ny:ny+ng,3) = gyt(ny-1,3)
      gty(ny:ny+ng) = gty(ny-1)
      gy(3)%flt1d(ny:ny+ng) = gy(3)%flt1d(ny-1)
    ENDIF
    IF ( myprocj == 1 ) THEN
      gyt(-ng+1:0,3) = gyt(1,3)
      gty(-ng+1:0) = gty(1)
      gy(3)%flt1d(-ng+1:0) = gy(3)%flt1d(1)
    ENDIF
!  ENDIF


  IF( bcy .eq. 0 ) THEN    ! mirror plane on north boundary

    vt(1:nx-1,ny,1:nz-1) = vt(1:nx-1,ny-1,1:nz-1)
   DO n = 1,ng
    vt( 1:nx,ny+n-1,1:nz) = vt(1:nx,ny-n ,1:nz) 
    wt( 1:nx,ny+n-1,1:nz) = wt(1:nx,ny-n ,1:nz)   
    ut( 1:nx,ny+n-1,1:nz) = ut(1:nx,ny-n ,1:nz)   
   ENDDO

!   gyt(ny,3) = gyt(ny-1,3)
!   gyt(ny,1) = gyt(ny-1,1)
   
!   gy(3)%flt1d(ny) = gy(3)%flt1d(1)
!   gy(1)%flt1d(ny) = gy(1)%flt1d(1)

  ENDIF


     IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  RK-DYNAMICS => LOOP/NSTEP/SDT:',loop,nstep,sdt

!==================================================================================================================================
! PART I A: PROGNOSTIC TKE

      IF ( debugsolver ) write(luno,*) 'SOLVER3D:  STARTING TKE'

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

      IF ( mix_type .eq. 1 ) THEN

       atype = atypes1 
       IF( loop .eq. RKSCHEME ) atype = atype2tke

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif
       CALL cld_cpu('ADVECT SCALAR')

       ft(:,:,:) = 0.0

       i = 1
!       IF ( atype == 12 ) i = 0
       
       CALL ADVECT(gd,kt,km%flt3d,ft,ut,vt,wt,                   &
                   gxt,gyt,gzt,rrp(-ng+1,1),rrm(-ng+1,1),sdt,dt, &
                   nx,ny,nz,                                  &
                   km,atype,i,time,time_real,ugrid,vgrid,x_sw_loc,y_sw_loc)

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif
       CALL cld_cpu('ADVECT SCALAR')

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  TKE ADVECTED'

       CALL cld_cpu('TKE')

       IF ( loop .eq. RKSCHEME .or. ltkeall ) THEN

        fu(:,:,:) = 0.0
        fv(:,:,:) = 0.0

          
        CALL TKE_SHEAR(fu,ut,vt,wt,gx,gy,gz,nx,ny,nz,ng)
        CALL TKE_BUOY (fv,s,st,piinit%flt1d,sbase(-ng+1,1),gzt,nx,ny,nz,ns)

        IF ( tke_type == 1 ) THEN ! km is km
        
        tflag = 0
        CALL MIX_SCAL (kt,ft,kt,kminit%flt1d,gxt,gyt,gzt,2.0,nx,ny,nz,t0,den,divv,'KM',tflag)  ! Note: Prtl is set to 2.0 on purpose here
        
        
         IF ( loop .eq. RKSCHEME .and. ifilt >= 1 ) THEN
!         IF ( loop .eq. RKSCHEME ) THEN
           CALL cld_cpu('FILTER') 

           CALL MFILTER6(kt,sbase(-ng+1,3),ft,km,sdt,nx,ny,nz,gx,gy)

           CALL cld_cpu('FILTER') 
         ENDIF

        CALL RAYDAMP  (kt,ft,   kminit%flt1d,gzt,0.0,nx,ny,nz,km)

! UPDATE TKE with rhs forcing (result is stored in kt), except for mix_type = 1, in which case the tendency is added to ft, and km is limited by ahighk
! shear=fu, buoy=fv, km%flt3d not actually used
        CALL COMPUTEKM(kt,km%flt3d,fu,fv,ft,ut,vt,wt,tke_diss,                               &
                        uinit%flt1d,vinit%flt1d,piinit%flt1d,sbase(-ng+1,1),kmbaserm%flt1d,  &
                        gxt,gyt,gzt,sdt,ugrid,vgrid,nx,ny,nz,ns)

        kmt(ixb:ixe,jyb:jye,kzb:kze) = kt(ixb:ixe,jyb:jye,kzb:kze)
        kht(ixb:ixe,jyb:jye,kzb:kze) = kt(ixb:ixe,jyb:jye,kzb:kze)

           IF ( ltkediss >= 1 ) THEN
             xtra(ltkediss)%flt3d(ixb:ixe,jyb:jye,kzb:kze) = tke_diss(ixb:ixe,jyb:jye,kzb:kze) ! fv(ixb:ixe,jyb:jye,kzb:kze) !
           ENDIF
        
        ELSEIF ( tke_type == 2 ) THEN ! km is SqrtE
        
        tflag = 0
        
        
        
         IF ( loop .eq. RKSCHEME .and. ifilt >= 1 ) THEN
!         IF ( loop .eq. RKSCHEME ) THEN
           CALL cld_cpu('FILTER') 

           CALL MFILTER6(kt,sbase(-ng+1,3),ft,km,sdt,nx,ny,nz,gx,gy)

           CALL cld_cpu('FILTER') 
         ENDIF
!        ENDIF

        CALL RAYDAMP  (kt,ft,   kminit%flt1d,gzt,0.0,nx,ny,nz,km)

! UPDATE TKE with rhs forcing (result is stored in kt), except for mix_type = 1, in which case the tendency is added to ft, and km is limited by ahighk
! shear=fu, buoy=fv, km%flt3d not actually used
!      SUBROUTINE COMPUTETKE(km,kt,shear,buoy,ft,u,v,w,kmt,kht,              &
!                            uinit,vinit,pinit,tinit,kmbaserm,               &
!                            gx,gy,gz,dt,ugrid,vgrid,nx,ny,nz,ns)
        CALL COMPUTETKE(kt,km%flt3d,fu,fv,ft,ut,vt,wt,kmt,kht,tke_diss,khh,khv,kmh,kmv,      &
                        uinit%flt1d,vinit%flt1d,piinit%flt1d,sbase(-ng+1,1),kmbaserm%flt1d,  &
                        gxt,gyt,gzt,sdt,ugrid,vgrid,nx,ny,nz,ns)

           IF ( lkmt >= 1 )  THEN
           
             !write(luno,*) 'set kmt/kht'
             xtra(lkmt)%flt3d(ixb:ixe,jyb:jye,kzb:kze) = kmt(ixb:ixe,jyb:jye,kzb:kze) ! fu(ixb:ixe,jyb:jye,kzb:kze) !
             xtra(lkht)%flt3d(ixb:ixe,jyb:jye,kzb:kze) = kht(ixb:ixe,jyb:jye,kzb:kze) ! fv(ixb:ixe,jyb:jye,kzb:kze) !
           
           ENDIF
           
           IF ( ltkediss >= 1 ) THEN
             xtra(ltkediss)%flt3d(ixb:ixe,jyb:jye,kzb:kze) = tke_diss(ixb:ixe,jyb:jye,kzb:kze) ! fv(ixb:ixe,jyb:jye,kzb:kze) !
           ENDIF

        ENDIF ! tke_type
       ENDIF
        ! calculate Km in kmt
        
        IF ( tke_type == 1 ) THEN
         CALL MIX_SCAL (kt,ft,kht,kminit%flt1d,gxt,gyt,gzt,2.0,nx,ny,nz,t0,den,divv,'TKE',tflag)  ! Note: Prtl is set to 2.0 on purpose here
        ELSEIF ( tke_type == 2 ) THEN
         CALL MIX_SCAL_HV(kt,ft,khh,khv,kminit%flt1d,gxt,gyt,gzt,2.0,nx,ny,nz,t0,den,divv,'TKE',tflag)  ! Note: Prtl is set to 2.0 on purpose here
        ENDIF

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  AFTER RHS TKE'

! UPDATE TKE advective tendencies
! Not doing this here anymore, rather at the end of the RK loop.  Otherwise have to communicate for MPI
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

!         do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!         DO k =1,nz-1 ; DO j = 1,ny-1 ; DO i = 1,nx-1
#endif
!          kt(i,j,k) = Min( ahighk, Max( kt(i,j,k) + sdt*ft(i,j,k), 0.0 ) )
!         ENDDO ; ENDDO ; ENDDO

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif
       CALL cld_cpu('TKE')

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  TKE UPDATED, LOOP = ', loop

      ENDIF

!==================================================================================================================================
! PART I B: SMAGORINSKY FORMULATION

      IF( mix_type .eq. 0 .and. loop .eq. 1 ) THEN

       CALL cld_cpu('TKE')

!        CALL TKE_SHEAR(fu,u%flt3d,v%flt3d,w%flt3d,gx,gy,gz,nx,ny,nz,ng)
        CALL TKE_SHEAR(fu,ut,vt,wt,gx,gy,gz,nx,ny,nz,ng)
        CALL TKE_BUOY (fv,s,st,piinit%flt1d,sbase(-ng+1,1),gzt,nx,ny,nz,ns)
        CALL COMPUTEKM(kt,km%flt3d,fu,fv,ft,ut,vt,wt,tke_diss,                               &
                        uinit%flt1d,vinit%flt1d,piinit%flt1d,sbase(-ng+1,1),kmbaserm%flt1d,  &
                        gxt,gyt,gzt,sdt,ugrid,vgrid,                          &
                        nx,ny,nz,ns)
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

        kmt(ixb:ixe,jyb:jye,kzb:kze) = kt(ixb:ixe,jyb:jye,kzb:kze)

       CALL cld_cpu('TKE')

       IF ( debugsolver ) write(0,*) 'SOLVER3D:  SMAGORINSKY KM COMPUTED...LOOP = ',loop

      ENDIF

!==================================================================================================================================
! PART I C: Constant value

      IF( mix_type .eq. -1 ) THEN

       CALL cld_cpu('TKE')
        
        IF ( debugsolver ) write(luno,*) 'SOLVER3D:  CONSTANT KM COMPUTED...LOOP = ',loop, ' ahighk = ', ahighk
        
!        kt(:,:,:) = 0.0
        
        CALL COMPUTEKM(kt,km%flt3d,fu,fv,ft,ut,vt,wt,tke_diss,                               &
                        uinit%flt1d,vinit%flt1d,piinit%flt1d,sbase(-ng+1,1), &
                        kmbaserm%flt1d,gxt,gyt,gzt,sdt,ugrid,vgrid,nx,ny,nz,ns)
       
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif
       CALL cld_cpu('TKE')

        kmt(ixb:ixe,jyb:jye,kzb:kze) = kt(ixb:ixe,jyb:jye,kzb:kze)
       
       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  CONSTANT KM COMPUTED...LOOP = ',loop
      
      ENDIF

!==================================================================================================================================
! PART II: Compute buoyancy forcing before scalar update


     CALL cld_cpu('SMLSTEP')  

#ifdef MPI
     CALL cld_cpu('SOLVER-WAIT-1')  
!     CALL MPI_BARRIER(my_comm, mpi_error_code)
     CALL cld_cpu('SOLVER-WAIT-1')  
#endif
      fp(:,:,:) = 0.0

      CALL BUOY(fp, s, st, sbase, nx, ny, nz, ns)

      IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST BUOYANCY'
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

     CALL cld_cpu('SMLSTEP')  

!==================================================================================================================================
! PART III:  SCALAR INTEGRATION

! SCALAR LOOP

     IF ( debugsolver ) write(luno,*) 'BEGIN SCALAR INTEGRATION, RK-LOOP = ', loop

        IF ( bcx .eq. 2 .and. loop .eq. RKSCHEME) THEN
         DO n = 1,ng
         
#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!          jyb = 1
!          jye = jtile
!          if (jyend .eq. nyend) jye = jyend-jybeg+1
!
!          do k = kzb,kze
!           do j = jyb,jye
!            kt( nx-1+n,j,k) = kt(n   ,j,k)
!            kt( 1-n   ,j,k) = kt(nx-n,j,k)
!           enddo
!          enddo

#else
          kt( nx-1+n,1:ny,1:nz) = kt(n   ,1:ny,1:nz)
          kt( 1-n   ,1:ny,1:nz) = kt(nx-n,1:ny,1:nz)
          kmt( nx-1+n,1:ny,1:nz) = kmt(n   ,1:ny,1:nz)
          kmt( 1-n   ,1:ny,1:nz) = kmt(nx-n,1:ny,1:nz)
#endif
         ENDDO
!        ELSE
!         DO n = 1,ng
!          kt( nx-1+n,1:ny,1:nz) = kt(nx-1,1:ny,1:nz)
!          kt( 1-n   ,1:ny,1:nz) = kt(1   ,1:ny,1:nz)
!         ENDDO
        ENDIF

        IF ( bcy .eq. 2 .and. loop .eq. RKSCHEME) THEN
         DO n = 1,ng
#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!          ixb = 1
!          ixe = itile
!          if (ixend .eq. nxend) ixe = ixend-ixbeg+1
!
!          do k = kzb,kze
!           do i = ixb,ixe
!            kt(i,ny-1+n,k) = kt(i,n,   k) 
!            kt(i, 1-n  ,k) = kt(i,ny-n,k)
!           enddo
!          enddo
#else
          kt(1:nx,ny-1+n,1:nz) = kt(1:nx,n,   1:nz) 
          kt(1:nx, 1-n  ,1:nz) = kt(1:nx,ny-n,1:nz)
          kmt(1:nx,ny-1+n,1:nz) = kmt(1:nx,n,   1:nz) 
          kmt(1:nx, 1-n  ,1:nz) = kmt(1:nx,ny-n,1:nz)
#endif
         ENDDO
        ENDIF


     IF ( nsrk > 0 ) THEN
     DO n = 1,nsrk

       IF ( debugsolver ) write(luno,*) 'SCALAR INTEGRATION, N = ',n, s(n)%name 
       IF ( debugsolver .and. debug_mpi ) write(0,*) my_rank, 'SCALAR INTEGRATION, N = ',n, s(n)%name 


! ADVECT SCALAR

       CALL cld_cpu('ADVECT SCALAR')

        IF( loop .eq. RKSCHEME ) THEN
          atype = atypes2
          IF (  s(n)%name .eq. 'TH' .and. atype2th /= 1 ) atype = atypem2
          IF (  s(n)%name .eq. 'QV' .and. atype2qv /= 1 ) atype = atypem2
        ELSE
          atype = atypes1
          IF (  s(n)%name .eq. 'TH' .and. atype2th /= 1 ) atype = atypem1
          IF (  s(n)%name .eq. 'QV' .and. atype2qv /= 1 ) atype = atypem1
        ENDIF
!        IF( loop .eq. RKSCHEME .and. MONOTONIC .and. s(n)%name .eq. 'TH' ) atype = atype2th
!        IF( loop .eq. RKSCHEME .and. MONOTONIC .and. s(n)%name .eq. 'QV' ) atype = atype2qv

! Here, i is the flag to find the minimum box around the non-zero values of a variable and do advection only within that box.

        i = 1  ! default to find box

!        IF ( atype == 12 ) i = 0

        IF ( s(n)%name .eq. 'TH' )    i = 0
        IF ( s(n)%name .eq. 'QV' )    i = 0
        IF ( s(n)%name .eq. 'CCCN' )  i = 0
        IF ( s(n)%name .eq. 'CPION' ) i = 0
        IF ( s(n)%name .eq. 'CNION' ) i = 0
        IF ( s(n)%name .eq. 'CCIN' )  i = 0

! Set fs = 0
#ifndef MPI
       IF ( i /= 0 ) THEN
       CALL cld_cpu('TIMESTEP') 
         fs(:,:,:) = 0.0
       CALL cld_cpu('TIMESTEP') 
       ENDIF
#else
       fs(:,:,:) = 0.0
#endif

        CALL ADVECT(gd,st(ge,ge,ge,n),s(n)%flt3d,fs,ut,vt,wt,           &
                    gxt,gyt,gzt,rrp(-ng+1,1),rrm(-ng+1,1),sdt,dt,       &
                    nx,ny,nz,                                        &
                    s(n),atype,i,time,time_real,ugrid,vgrid,x_sw_loc,y_sw_loc)
                    
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif
       
        CALL cld_cpu('ADVECT SCALAR')

        IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST SCALAR ADVECTION VARIABLE = ',n, s(n)%name 
        IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  PAST SCALAR ADVECTION VARIABLE = ',n, s(n)%name 

! FIT terms on last step... (mixing, rayleigh damping, etc)

        CALL cld_cpu('MIX')

         IF ( loop .eq. RKSCHEME .and. ltkeon ) THEN
         
#ifdef CHGELEC

         IF ( loop .eq. RKSCHEME ) THEN

          IF ( ichgtndadv > 1 ) THEN
           IF (  s(n)%name(1:2) == 'SC' ) THEN
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndadv)%flt3d(i,j,k) = elec(ichgtndadv)%flt3d(i,j,k) + sdt*den(k,1)*fs(i,j,k)
                   chgadvtemp(i,j,k,1) = chgadvtemp(i,j,k,1) + sdt*den(k,1)*fs(i,j,k)
               enddo; enddo; enddo
           ELSEIF ( s(n)%name(1:5) == 'CPION' ) THEN
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndadv)%flt3d(i,j,k) = elec(ichgtndadv)%flt3d(i,j,k) + ec*sdt*den(k,1)*fs(i,j,k)
                   chgadvtemp(i,j,k,1) = chgadvtemp(i,j,k,1) + ec*sdt*den(k,1)*fs(i,j,k)
               enddo; enddo; enddo
           ELSEIF ( s(n)%name(1:5) == 'CNION' ) THEN
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndadv)%flt3d(i,j,k) = elec(ichgtndadv)%flt3d(i,j,k) - ec*sdt*den(k,1)*fs(i,j,k)
                   chgadvtemp(i,j,k,1) = chgadvtemp(i,j,k,1) - ec*sdt*den(k,1)*fs(i,j,k)
               enddo; enddo; enddo
           ENDIF
          ENDIF

          IF ( ichgtndadvsn > 1 ) THEN
          
!           dbz%flt3d(:,:,:) = 0.0
!           
!           IF (  s(n)%name(1:2) == 'SC' ) THEN
!               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
!                   elec(ichgtndadvsn)%flt3d(i,j,k) = elec(ichgtndadvsn)%flt3d(i,j,k) + sdt*den(k,1)*fs(i,j,k)
!               enddo; enddo; enddo
!           ELSEIF ( s(n)%name(1:5) == 'CPION' ) THEN
!               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
!                   elec(ichgtndadvsn)%flt3d(i,j,k) = elec(ichgtndadvsn)%flt3d(i,j,k) + ec*sdt*den(k,1)*fs(i,j,k)
!               enddo; enddo; enddo
!           ELSEIF ( s(n)%name(1:5) == 'CNION' ) THEN
!               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
!                   elec(ichgtndadvsn)%flt3d(i,j,k) = elec(ichgtndadvsn)%flt3d(i,j,k) - ec*sdt*den(k,1)*fs(i,j,k)
!               enddo; enddo; enddo
!           ENDIF
          ENDIF

         ENDIF

           IF ( ichgtndmix > 1 .and. ( s(n)%name(1:2) == 'SC'  .or. s(n)%name(1:5) == 'CNION'  &
                                       .or. s(n)%name(1:5) == 'CPION' ) ) THEN
             dbz%flt3d(:,:,:) = fs(:,:,:)
           ENDIF

#endif

           IF ( s(n)%phytype .eq. 0 ) THEN
             IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  CALL MIX_SCAL = ',n, s(n)%name 
             IF ( tke_type == 1 ) THEN
               IF ( microphys(1:3) == 'TAK' .and. lss > 0 .and. ( s(n)%name .eq. 'TH' .or. s(n)%name .eq. 'QV') ) THEN
                 tflag = 1
                 IF ( s(n)%name .eq. 'TH' ) THEN ! store TH mixing tendency in dbz array
                   CALL MIX_SCAL(st(ge,ge,ge,n),fs,kht,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,dbz%flt3d,s(n)%name,tflag)
                 ELSEIF (s(n)%name .eq. 'QV' ) THEN ! store QV mixing tendency in wz array
                   CALL MIX_SCAL(st(ge,ge,ge,n),fs,kht,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,wz%flt3d,s(n)%name,tflag)
                 ENDIF
               ELSE
                 tflag = 0
                 CALL MIX_SCAL(st(ge,ge,ge,n),fs,kht,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,divv,s(n)%name,tflag)
               ENDIF
             ELSEIF ( tke_type == 2 ) THEN
               IF ( microphys(1:3) == 'TAK' .and. lss > 0 .and. ( s(n)%name .eq. 'TH' .or. s(n)%name .eq. 'QV') ) THEN
                 tflag = 1
                 IF ( s(n)%name .eq. 'TH' ) THEN ! store TH mixing tendency in dbz array
                   CALL MIX_SCAL_HV(st(ge,ge,ge,n),fs,khh,khv,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,dbz%flt3d,s(n)%name,tflag)
                 ELSEIF (s(n)%name .eq. 'QV' ) THEN ! store QV mixing tendency in wz array
                   CALL MIX_SCAL_HV(st(ge,ge,ge,n),fs,khh,khv,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,wz%flt3d,s(n)%name,tflag)
                 ENDIF
               ELSE
                 tflag = 0
                 CALL MIX_SCAL_HV(st(ge,ge,ge,n),fs,khh,khv,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,divv,s(n)%name,tflag)
               ENDIF
             ENDIF
           ENDIF

#ifdef CHGELEC
           IF ( ichgtndmix > 1 .and. s(n)%name(1:2) == 'SC' ) THEN
!             IF ( my_rank == 0 ) write(0,*) 'solver: ichgtndmix, name = ',s(n)%name
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndmix)%flt3d(i,j,k) = elec(ichgtndmix)%flt3d(i,j,k) + sdt*den(k,1)*( fs(i,j,k) - dbz%flt3d(i,j,k) )
               enddo; enddo; enddo
           ENDIF
           IF ( ichgtndmix > 1 .and. s(n)%name(1:5) == 'CPION' ) THEN
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndmix)%flt3d(i,j,k) = elec(ichgtndmix)%flt3d(i,j,k) + ec*sdt*den(k,1)*( fs(i,j,k) - dbz%flt3d(i,j,k) )
               enddo; enddo; enddo
           ENDIF
           IF ( ichgtndmix > 1 .and. s(n)%name(1:5) == 'CNION' ) THEN
               do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
                   elec(ichgtndmix)%flt3d(i,j,k) = elec(ichgtndmix)%flt3d(i,j,k) - ec*sdt*den(k,1)*( fs(i,j,k) - dbz%flt3d(i,j,k) )
               enddo; enddo; enddo
           ENDIF
#endif

           IF ( s(n)%name .eq. 'TH' .or. s(n)%name .eq. 'QV' ) THEN
             IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  CALL RAYDAMP = ',n, s(n)%name 
             CALL RAYDAMP(st(ge,ge,ge,n),fs,sbase(-ng+1,n),gzt,0.0,nx,ny,nz,s(n))
           ENDIF

           IF ( s(n)%name .eq. 'TH' .and. rotunno_radiation > 0 ) THEN
             ! Pseudo-radiation from Rotunno and Emanuel 1987 (JAS)
             ! limited to rate of 2 deg per day (1 deg per 12 hr)
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

             DO k = kzb,kze
               DO j = jyb,jye
                 DO i = ixb,ixe
                  tmp = Min( 1.0, Max( -1.0, s(n)%flt3d(i,j,k) - sbase(k,n) ) )
                  fs(i,j,k) = fs(i,j,k) - tmp/(12.*3600.)
                 ENDDO 
               ENDDO 
             ENDDO
             
           ENDIF

           IF ( s(n)%name .eq. 'QV' .and. rotunno_radiation >= 2 ) THEN
             ! hack to nudge QV
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

             DO k = kzb,kze
               DO j = jyb,jye
                 DO i = ixb,ixe
                  tmp = Min( 0.0, Max( -1.0, s(n)%flt3d(i,j,k) - sbase(k,n) ) )
                  fs(i,j,k) = fs(i,j,k) - tmp/(12.*3600.)
                 ENDDO 
               ENDDO 
             ENDDO
             
           ENDIF


         ENDIF



         IF ( loop .eq. RKSCHEME .and. ifilt >= 1 ) THEN
           CALL cld_cpu('FILTER') 

           CALL MFILTER6(st(ge,ge,ge,n),sbase(-ng+1,n),fs,s(n),sdt,nx,ny,nz,gx,gy)

           CALL cld_cpu('FILTER') 
         ENDIF

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif


        CALL cld_cpu('MIX')

        IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST SCALAR MIXING/RAYD, N = ', n, s(n)%name 
#ifdef MPI
        IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  PAST SCALAR MIXING/RAYD, N = ', n, s(n)%name 
        IF ( debugsolver ) CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
 
        CALL cld_cpu('TIMESTEP') 

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
           do j = jyb,jye ; do i = ixb,ixe
#else
         DO k =1,nz-1 ; DO j = 1,ny-1 ; DO i = 1,nx-1 
#endif
           st(i,j,k,n) = s(n)%flt3d(i,j,k) + sdt*fs(i,j,k)
         ENDDO ; ENDDO ; ENDDO

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

        CALL cld_cpu('TIMESTEP') 

        IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST SCALAR UPDATE, nscalar = ',n
#ifdef MPI
        IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:  PAST SCALAR UPDATE, nscalar = ',n
        IF ( debugsolver ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

        ENDDO
        
        ENDIF ! n = 1,nsrk
 
!        IF ( ihole .eq. 1 .or. ( atypes1 .eq. 0 .and. loop .lt. RKSCHEME ) ) THEN
        IF ( ihole .eq. 1 ) THEN
          CALL cld_cpu('HOLE-FILLING') 

          CALL HOLEFILL()

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

          CALL cld_cpu('HOLE-FILLING') 
        ENDIF

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  COMPLETED SCALAR UPDATE, RK-LOOP = ',loop
       IF ( debugsolver .and. debug_mpi ) write(0,*) my_rank, 'SOLVER3D:  COMPLETED SCALAR UPDATE, RK-LOOP = ',loop

!==================================================================================================================================
! PART IV: VELOCITY & PRESSURE

       CALL cld_cpu('TIMESTEP')  
#ifdef MPI
        
       fv(:,:,:) = 0.0
       fu(:,:,:) = 0.0
       fw(:,:,:) = 0.0

#endif

       CALL cld_cpu('TIMESTEP')  

       IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  FU/FV/FW = 0'

! COMPUTE THE ADVECTION OF U/V/W

       CALL cld_cpu('ADVECT UVW')

        atype = atypem1
        IF( loop .eq. RKSCHEME ) atype = atypem2

        IF ( iuvwadv .eq. 1 ) THEN
        
          IF ( nx .gt. 2 ) THEN

            
            IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  FU'
            CALL ADVECT(gd,ut,u%flt3d,fu,ut,vt,wt,                     &
                        gxt,gyt,gzt,rrp(-ng+1,1),rrm(-ng+1,1),sdt,dt,  &
                        nx,ny,nz,                                   &
                        u,atype,0,time,time_real,ugrid,vgrid,x_sw_loc,y_sw_loc)
          ENDIF
           
          IF ( ny .gt. 2 ) THEN

            IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  FV'
            CALL ADVECT(gd,vt,v%flt3d,fv,ut,vt,wt,                     &
                        gxt,gyt,gzt,rrp(-ng+1,1),rrm(-ng+1,1),sdt,dt,  &
                        nx,ny,nz,                                   &
                        v,atype,0,time,time_real,ugrid,vgrid,x_sw_loc,y_sw_loc)
          ENDIF

            IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  FW'
            CALL ADVECT(gd,wt,w%flt3d,fw,ut,vt,wt,                     &
                        gxt,gyt,gzt,rrp(-ng+1,2),rrm(-ng+1,2),sdt,dt,  &
                        nx,ny,nz,                                   &
                        w,atype,0,time,time_real,ugrid,vgrid,x_sw_loc,y_sw_loc)


            IF ( debugsolver ) write(0,*) my_rank,'SOLVER3D:  DONE FW'

#ifdef MPI
           IF ( debugsolver ) CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
        
        ELSEIF ( iuvwadv .eq. 2 ) THEN
        
          fu(:,:,:) = 0.0
          fv(:,:,:) = 0.0
          fw(:,:,:) = 0.0
         CALL boxadvect( nx,ny,nz,ng,         &
                         sdt,dx,dy,dz,gtz,    &
                         ut,u%flt3d,          &  ! note ut and u%flt3d are received as 4-D arrays (u,v,w)
                         fu, fv, fw   )
        
        ELSEIF ( iuvwadv >= 3 ) THEN ! to keep winds constant -- kinematic model
          fu(:,:,:) = 0.0
          fv(:,:,:) = 0.0
          fw(:,:,:) = 0.0
          fp(:,:,:) = 0.0
          nstep = 1
          dts = sdt

        ELSEIF ( iuvwadv == 0 ) THEN ! to keep winds constant -- kinematic model
          fu(:,:,:) = 0.0
          fv(:,:,:) = 0.0
          fw(:,:,:) = 0.0
!          fp(:,:,:) = 0.0
          
        ENDIF

#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

       CALL cld_cpu('ADVECT UVW')

       CALL cld_cpu('W_DAMP')

            CALL W_DAMP(fw)
       CALL cld_cpu('W_DAMP')

         IF ( iforce .ne. 0 ) THEN
          CALL cld_cpu('FORCE')

          IF ( iwforce .ne. 0 ) THEN
          CALL WFORCE (nx,ny,nz,sdt,w%flt3d,wt,fw,u%flt3d,ut,fu,t0,         &
                       gxt,gyt,gzt,loop,time,time_real,       &
                       uinit0%flt1d,vinit0%flt1d,ugrid0,vgrid0,ranarr)
          ENDIF
          
          CALL cld_cpu('FORCE')
         ENDIF

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST UVW ADVECT'

! FIT terms on last step...

       CALL cld_cpu('MIX')

        IF ( loop .eq. RKSCHEME .and. ltkeon ) THEN

         IF ( tke_type == 1 ) THEN
           CALL MIX_VELO(ut,vt,wt,fu,fv,fw,kmt,kmh,kmv,kt,fs,uinit%flt1d,vinit%flt1d,gxt,gyt,gzt,ugrid,vgrid,nx,ny,nz,t0,den)
         ELSE
           CALL MIX_VELO(ut,vt,wt,fu,fv,fw,kmt,kmh,kmv,kt,fs,uinit%flt1d,vinit%flt1d,gxt,gyt,gzt,ugrid,vgrid,nx,ny,nz,t0,den)
         ENDIF
         CALL RAYDAMP(ut,fu,uinit%flt1d,gzt,ugrid,nx,ny,nz,u)
         CALL RAYDAMP(vt,fv,vinit%flt1d,gzt,vgrid,nx,ny,nz,v)
         CALL RAYDAMP(wt,fw,winit%flt1d,gzt,0.0,  nx,ny,nz,w)

        ENDIF
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

       CALL cld_cpu('MIX')

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST MOMENTUM MIXING/RAYD'

! CORIOLIS ACCELERATION

       IF ( coriol .ne. 0.0 ) THEN
       CALL cld_cpu('TIMESTEP')  

        CALL CORIOLIS(ut,vt,wt,fu,fv,fw,uinit%flt1d,vinit%flt1d,ugrid,vgrid,lat,nx,ny,nz)

       CALL cld_cpu('TIMESTEP')  
       ENDIF

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST CORIOLIS'


! SMALL TIME STEP
! MPI: if mpi is used, u/v/w/p ghost zones are filled in small step

       IF ( iuvwadv < 3 ) THEN !{
       CALL cld_cpu('SMLSTEP')



#ifdef MPI

!     CALL cld_cpu('SOLVER-WAIT-2')  
!     CALL MPI_BARRIER(my_comm, mpi_error_code)
!     CALL cld_cpu('SOLVER-WAIT-2')  

       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       do k = kzb,kze
#else
        DO k = 1,nz-1
#endif
          fw(:,:,k) = fw(:,:,k) + fp(:,:,k)
        ENDDO

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  BEFORE SMLSTEP'

#ifndef MPI
        IF ( loop .eq. RKSCHEME .and. print_noise ) THEN
         dpdt(:,:,:) = pi%flt3d(1:nx-1,1:ny-1,1:nz-1)
        ENDIF
#endif

        CALL SMLSTEP(u%flt3d, ut, fu,                                 &
                     v%flt3d, vt, fv,                                 &
                     w%flt3d, wt, fw,                                 &
                     pi%flt3d, pt, fp,                                &
                     piinit%flt1d, sbase(-ng+1,1), sbase(-ng+1,2),    &
                     uinit%flt1d, vinit%flt1d,                        &
                     gxt, gyt, gzt,                                   &
                     dts, nx, ny, nz, nstep)


#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

       CALL cld_cpu('SMLSTEP')

#ifndef MPI
       IF ( loop .eq. RKSCHEME .and. print_noise ) THEN
          dpdt(:,:,:) = abs((pt(1:nx-1,1:ny-1,1:nz-1) - dpdt(:,:,:))) / dt
       ENDIF
#endif

       
       IF ( detrendpi ) THEN
       CALL cld_cpu('DTREND')  

        CALL DTREND()
#ifdef MPI
        IF ( truetime ) THEN
          CALL MPI_BARRIER(my_comm, mpi_error_code)
        ENDIF
#endif

       CALL cld_cpu('DTREND')  
       
       ENDIF ! detrendpi

       
       ENDIF !} iuvwadv < 3
       
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
         rho%flt3d(i,j,k)   =  1.0e5*(piinit%flt1d(k) + pi%flt3d(i,j,k))**cvr/ &
              (rd*s(lt)%flt3d(i,j,k) * (1.0 + 0.61*s(lv)%flt3d(i,j,k)))
         
         ENDDO ; ENDDO ; ENDDO

       
#ifdef MPI
        IF ( debugsolver ) CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif

       IF ( debugsolver ) write(luno,*) 'SOLVER3D:  PAST SMLSTEP'

      IF ( mix_type .eq. 1 ) THEN
! UPDATE TKE advective tendencies
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
         DO k =1,nz-1 ; DO j = 1,ny-1 ; DO i = 1,nx-1
#endif
!          kt(i,j,k) = Min( ahighk, Max( kt(i,j,k) + sdt*ft(i,j,k), 0.0 ) )
          kt(i,j,k) = kt(i,j,k) + sdt*ft(i,j,k)
          
          IF ( iturbenhance == 0 .and. tke_type == 1 ) THEN
            kt(i,j,k) = max(min(kt(i,j,k), ahighk),0.0) 
          ELSE
            kt(i,j,k) = Max(0.0, kt(i,j,k) )
          ENDIF

          IF ( tke_type == 1 ) THEN
            IF ( kt(i,j,k) < km_thresh ) kt(i,j,k) = 0.0
            kmt(i,j,k) = kt(i,j,k)
          ENDIF

         ENDDO ; ENDDO ; ENDDO
      
      ENDIF

! MPI: Update scalars at end of RK LOOP

#ifdef MPI

    IF ( number_of_processes .gt. 1  ) THEN

    CALL cld_cpu('MPICOM-SOLVER')

    IF ( loop .lt. RKSCHEME ) THEN

      nampi = 1 + nsrk

! Note: kt is the 3d array at index 0 of the st arrary (allocated as one block)

   IF ( nproci > 1 ) THEN

   westward_tag = 1001
   CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),westward_tag,kt(-ng+1,-ng+1,-ng+1) )

   eastward_tag = 1002
   CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),eastward_tag,kt(-ng+1,-ng+1,-ng+1) )

   ENDIF

   IF ( nprocj > 1 ) THEN
   
   southward_tag = 1003
   CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),southward_tag,kt(-ng+1,-ng+1,-ng+1) )

   northward_tag = 1004
   CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),northward_tag,kt(-ng+1,-ng+1,-ng+1) )

   ENDIF
   
   IF ( nprock > 1 ) THEN

   downward_tag = 1005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,kt(-ng+1,-ng+1,-ng+1))

   upward_tag = 1006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,kt(-ng+1,-ng+1,-ng+1))
   
   ENDIF

    
    ENDIF

    IF ( truetime ) THEN
    
     CALL cld_cpu('SOLVER-WAIT-3')  
     CALL MPI_BARRIER(my_comm, mpi_error_code)
     CALL cld_cpu('SOLVER-WAIT-3')  
     
    ENDIF

    CALL cld_cpu('MPICOM-SOLVER')

    
    ENDIF

#endif


!==================================================================================================================================
! END MAIN RK-LOOP 

   IF ( debugsolver ) write(luno,*) 'SOLVER3D:  END OF RK-LOOP',loop

   ENDDO

!==================================================================================================================================
! If print_noise == True, dump out stats for dpdt here.  Note - does not work in MPI mode

#ifndef MPI
  IF ( print_noise ) THEN
    write(luno,"(a6,2x,i6,2x,4(g14.6,2x))") "MDPDT:", time, SUM(dpdt(:,:,1))/(float(nx-1)*float(ny-1)), &
                                                            SUM(dpdt)/(float(nx-1)*float(ny-1)*float(nz-1))
  ENDIF
#endif


!==================================================================================================================================
! FORWARD IN TIME ADVECTION/RHS for SCALARS

   IF ( nsfwd .eq. ns ) THEN

  IF( bcx .eq. 2 .and. nproci == 1 ) THEN    ! periodic
#ifdef MPI
    IF ( nproci == 1 ) THEN
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg

    do k = kzb,kze ; do j = jyb,jye
     ut(nx,j,k) = ut(1,j,k)
     u%flt3d(nx,j,k) = u%flt3d(1,j,k)
    enddo ; enddo
    ENDIF
#else
    ut(nx,1:ny-1,1:nz-1) = ut(1,1:ny-1,1:nz-1)
    u%flt3d(nx,1:ny-1,1:nz-1) = u%flt3d(1,1:ny-1,1:nz-1)
#endif
   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg+1

    do k = kzb,kze ; do j = jyb,jye
     ut( nx+n,j,k) = ut(n+1 ,j,k)
     ut( 1-n ,j,k) = ut(nx-n,j,k) 

     wt( 1-n   ,j,k) = wt(nx-n,j,k) 
     wt( nx+n-1,j,k) = wt(n   ,j,k)
     vt( 1-n   ,j,k) = vt(nx-n,j,k) 
     vt( nx+n-1,j,k) = vt(n   ,j,k)
    enddo ; enddo
#else
    ut( nx+n,1:ny,1:nz) = ut(n+1   ,1:ny,1:nz)
    ut( 1-n ,1:ny,1:nz) = ut(nx-n,1:ny,1:nz) 

    wt( 1-n   ,1:ny,1:nz) = wt(nx-n,1:ny,1:nz) 
    wt( nx+n-1,1:ny,1:nz) = wt(n     ,1:ny,1:nz)
    vt( 1-n   ,1:ny,1:nz) = vt(nx-n,1:ny,1:nz) 
    vt( nx+n-1,1:ny,1:nz) = vt(n     ,1:ny,1:nz)
#endif
   ENDDO
      
  ENDIF

  IF( bcy .eq. 2 .and. nprocj == 1  ) THEN    ! periodic
#ifdef MPI
!    kzb = 1
!    kze = ktile
!    if (kzend .eq. nzend) kze = kzend-kzbeg
!
!    ixb = 1
!    ixe = itile
!    if (ixend .eq. nxend) ixe = ixend-ixbeg
!
!    do k = kzb,kze ; do i = ixb,ixe
!     vt(i,ny,k) = vt(i,1,k)
!     v%flt3d(i,ny,k) = v%flt3d(i,1,k)
!    enddo ; enddo
#else
    vt(1:nx-1,ny,1:nz-1) = vt(1:nx-1,1,1:nz-1)
    v%flt3d(1:nx-1,ny,1:nz-1) = v%flt3d(1:nx-1,1,1:nz-1)
#endif
   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg+1

    do k = kzb,kze ; do i = ixb,ixe
     vt(i, ny+n ,k)   = vt(i,n+1 ,k) 
     vt(i,  1-n ,k)   = vt(i,ny-n,k) 
     wt(i,  1-n ,k)   = wt(i,ny-n,k) 
     wt(i,ny+n-1,k)   = wt(i,  n ,k)   
     ut(i,  1-n ,k)   = ut(i,ny-n,k) 
     ut(i,ny+n-1,k)   = ut(i,  n ,k)   
    enddo ; enddo
#else
    vt(1:nx,ny+n ,1:nz)   = vt(1:nx,n+1 ,1:nz) 
    vt(1:nx, 1-n ,1:nz)   = vt(1:nx,ny-n,1:nz) 
    wt( 1:nx,1-n ,1:nz)   = wt(1:nx,ny-n,1:nz) 
    wt( 1:nx,ny+n-1,1:nz) = wt(1:nx,n ,1:nz)   
    ut( 1:nx,1-n ,1:nz)   = ut(1:nx,ny-n,1:nz) 
    ut( 1:nx,ny+n-1,1:nz) = ut(1:nx,n ,1:nz)   
#endif
   ENDDO

  ENDIF


      nst = Int(time_real/dt) + 1
!      nst = nst + 1
             
      CALL SETFUVW (nx,ny,nz,dx,dy,dz,dt,  &
                    u%flt3d,ut,v%flt3d,vt,w%flt3d,wt,gtx,gty,gtz,fu,fv,fw)


#ifdef MPI
! communicate vertically if needed

   IF ( nprock > 1 ) THEN

      nampi = 1 
      nb = 1

   downward_tag = 1005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,kt(-ng+1,-ng+1,-ng+1))

   upward_tag = 1006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,kt(-ng+1,-ng+1,-ng+1))
   
   ENDIF

#endif

        IF ( bcx .eq. 2 .and. nproci == 1 ) THEN
         DO n = 1,ng
#ifdef MPI
         kzb = 1
         kze = ktile
         if (kzend .eq. nzend) kze = kzend-kzbeg+1

         jyb = 1
         jye = jtile
         if (jyend .eq. nyend) jye = jyend-jybeg+1

         do k = kzb,kze ; do j = jyb,jye
          kmt( nx-1+n,j,k) = kmt(n   ,j,k)
          kmt( 1-n   ,j,k) = kmt(nx-n,j,k)
         enddo; enddo
#else
          kmt( nx-1+n,1:ny,1:nz) = kmt(n   ,1:ny,1:nz)
          kmt( 1-n   ,1:ny,1:nz) = kmt(nx-n,1:ny,1:nz)
#endif
         ENDDO
!        ELSE
!         DO n = 1,ng
!          kmt( nx-1+n,1:ny,1:nz) = kmt(nx-1,1:ny,1:nz)
!          kmt( 1-n   ,1:ny,1:nz) = kmt(1   ,1:ny,1:nz)
!         ENDDO
        ENDIF

        IF ( bcy .eq. 2 .and. nprocj == 1 ) THEN
         DO n = 1,ng
#ifdef MPI
!         kzb = 1
!         kze = ktile
!         if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!         ixb = 1
!         ixe = itile
!         if (ixend .eq. nxend) ixe = ixend-ixbeg+1
!
!         do k = kzb,kze ; do i = ixb,ixe
!          kmt(i,ny-1+n,k) = kmt(i,n,   k) 
!          kmt(i, 1-n  ,k) = kmt(i,ny-n,k)
!         enddo; enddo
#else
          kmt(1:nx,ny-1+n,1:nz) = kmt(1:nx,n,   1:nz) 
          kmt(1:nx, 1-n  ,1:nz) = kmt(1:nx,ny-n,1:nz)
#endif
         ENDDO
        ENDIF

! ADVECT SCALAR

                   
     CALL cld_cpu('ADVECT SCALAR')
     CALL cld_cpu('ADVECT CROWLEY')

      DO j = 1,icrwmp
      DO n = nsrk+1,nsfwd

! Here, i is the flag to find the minimum box around the non-zero values of a variable and do advection only within that box.

        i     = 1
        iadiv = 1 ! s(n)%dyntype

        IF ( s(n)%name .eq. 'TH' )  i = 0
        IF ( s(n)%name .eq. 'QV' )  i = 0
        IF ( s(n)%name .eq. 'CCCN' )  i = 0
        IF ( s(n)%name .eq. 'CPION' ) i = 0
        IF ( s(n)%name .eq. 'CNION' ) i = 0

#ifdef MPI
!         i = 0
#endif

!  IF ( bcx .eq. 2 ) THEN  ! periodic
!
!   DO ia = 1,ng
!    st( nx-1+ia,1:ny,1:nz,n) = st(ia   ,1:ny,1:nz,n)
!    st( 1-ia   ,1:ny,1:nz,n) = st(nx-ia,1:ny,1:nz,n)
!   ENDDO
!   
!  ENDIF
!
!  IF( bcy .eq. 2 ) THEN    ! periodic
!
!   DO ia = 1,ng
!    st(1:nx-1,ny-1+ia,1:nz-1,n) = st(1:nx-1,ia   ,1:nz-1,n) 
!    st(1:nx-1,1-ia   ,1:nz-1,n) = st(1:nx-1,ny-ia,1:nz-1,n)
!   ENDDO
!   
!  ENDIF

! #ifndef MPI
       CALL ADVECTCRW(nx, ny, nz, ns, nst,             &
                      dt, dx, dy, dz,  gtx, gty, gtz,  &
                      sbase(-ng+1,n), s(n)%flt3d, st, den(-ng+1,1), u%flt3d, v%flt3d, w%flt3d,  &
                      ft,  fu, fv, fw, iadiv, n, i, t0, s(n),j) 
! #endif
      ENDDO !! n

#ifdef MPI
! do communication for 6-pass option if needed
      IF ( icrwmp == 2 .and. j == 1 ) THEN
      nampi = nsfwd - nsrk

   IF ( nproci > 1 ) THEN

   westward_tag = 1001
   CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),westward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1) )

   eastward_tag = 1002
   CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),eastward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1) )

   ENDIF

   IF ( nprocj > 1 ) THEN
   
   southward_tag = 1003
   CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),southward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1) )

   northward_tag = 1004
   CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),northward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1) )

   ENDIF
   
   IF ( nprock > 1 ) THEN

   downward_tag = 1005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1))

   upward_tag = 1006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,st(-ng+1,-ng+1,-ng+1,nsrk+1))
   
   ENDIF
      
      
      ENDIF
#endif
      ENDDO !! j

     CALL cld_cpu('ADVECT CROWLEY')
     CALL cld_cpu('ADVECT SCALAR')
     


     DO n = nsrk+1,nsfwd

       IF ( debugsolver ) write(0,*) 'SCALAR INTEGRATION, N = ',n

! Copy scalar arrays into temporary storage that have computational zones

       CALL cld_cpu('TIMESTEP') 

#ifdef MPI
!        kzb = 1
!        kze = ktile
!        if (kzend .eq. nzend) kze = kzend-kzbeg

!        do k = kzb,kze
#else
!        DO k = 1,nz-1
#endif
!          fs(:,:,k) = 0.0
!        ENDDO

          fs(:,:,:) = 0.0

       CALL cld_cpu('TIMESTEP') 

       IF ( debugsolver ) write(0,*) 'SOLVER3D:  PAST CROWLEY SCALAR ADVECTION VARIABLE = ',n
      
! FIT terms on last step... (mixing, rayleigh damping, etc)

       IF ( iuvwadv < 3 ) THEN
       CALL cld_cpu('MIX')
       CALL cld_cpu('MIX-CRW')

           IF ( s(n)%phytype .eq. 0 ) THEN
!            CALL MIX_SCAL(st(ge,ge,ge,n),fs,kmt,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den)
            tflag = 0
            IF ( tke_type == 1 ) THEN
              CALL MIX_SCAL(s(n)%flt3d,fs,kht,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,divv,s(n)%name,tflag)
            ELSE
              CALL MIX_SCAL_HV(s(n)%flt3d,fs,khh,khv,sbase(-ng+1,n),gxt,gyt,gzt,Pr,nx,ny,nz,t0,den,divv,s(n)%name,tflag)
            ENDIF
           ENDIF
         
          IF ( s(n)%name .eq. 'TH' .or. s(n)%name .eq. 'QV' ) THEN
!            CALL RAYDAMP(st(ge,ge,ge,n),fs,sbase(-ng+1,n),gzt,0.0,nx,ny,nz,s(n))
            CALL RAYDAMP(s(n)%flt3d,fs,sbase(-ng+1,n),gzt,0.0,nx,ny,nz,s(n))
          ENDIF

       CALL cld_cpu('MIX-CRW')
       CALL cld_cpu('MIX')
       ENDIF

       IF ( debugsolver ) write(0,*) 'SOLVER3D:  PAST CROWLEY SCALAR MIXING/RAYD, N = ', n
      
! UPDATE SCALAR

      IF ( ifilt >= 1 ) THEN
       CALL cld_cpu('FILTER') 
       
       CALL MFILTER6(st(ge,ge,ge,n),sbase(-ng+1,n),fs,s(n),dt,nx,ny,nz,gx,gy)
       
       CALL cld_cpu('FILTER') 
      ENDIF
 
       CALL cld_cpu('TIMESTEP') 

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
        DO k =1,nz-1 ; DO j = 1,ny-1 ; DO i = 1,nx-1 
#endif
          st(i,j,k,n) = st(i,j,k,n) + dt*fs(i,j,k)
        ENDDO ; ENDDO ; ENDDO

       CALL cld_cpu('TIMESTEP') 

       IF ( debugsolver ) write(0,*) 'SOLVER3D:  PAST FIT SCALAR UPDATE, nscalar = ',n

     ENDDO !! n

   ENDIF !! ( nsfwd .eq. ns )


         IF ( iforce .ne. 0 ) THEN
          CALL cld_cpu('FORCE')
          
          IF ( issforce .ne. 0 ) THEN
          CALL SSFORCE (nx,ny,nz,ns,sdt,w%flt3d,st,pt,t0,piinit%flt1d(-ng+1),sbase,       &
                       gxt,gyt,gzt,time,time_real,loop,den,      &
                       uinit%flt1d,vinit%flt1d,ugrid,vgrid)
          ENDIF

          IF ( itopforce .ne. 0 ) THEN
           IF ( microphys(1:3) == 'TAK' ) THEN
             CALL BINFORCE (nx,ny,nz,ns,sdt,st,pt,t0,piinit%flt1d(-ng+1),sbase,den,    &
                         gxt,gyt,gzt,time,time_real, onedoutput)
           ELSE
             CALL QFORCE (nx,ny,nz,ns,sdt,st,pt,t0,piinit%flt1d(-ng+1),sbase,den,    &
                         gxt,gyt,gzt,time,time_real)
           ENDIF
          ENDIF

          IF ( iccnufforce .ne. 0 ) THEN
             CALL CCNUFFORCE (nx,ny,nz,ns,sdt,st,pt,t0,piinit%flt1d(-ng+1),sbase,den,    &
                         gxt,gyt,gzt,time,time_real)
          ENDIF

          IF ( ichgforce .ne. 0 ) THEN
             CALL CHARGEFORCE (nx,ny,nz,ns,sdt,st,pt,t0,piinit%flt1d(-ng+1),sbase,den,    &
                         gxt,gyt,gzt,time,time_real)
          ENDIF
          
          CALL cld_cpu('FORCE')
         ENDIF
         
   IF ( bogusvortex == 1 .and. ( time_real <= timint + dt*2 )  ) THEN
          CALL cld_cpu('FORCE')
            CALL HURRFORCE(gd, time_real,nx,ny,nz,ng,ut=ut,vt=vt)
          CALL cld_cpu('FORCE')
   ENDIF
   
   IF ( low_level_cooling_flag >= 1 ) THEN
             CALL BLCOOL (nx,ny,nz,ns,sdt,w%flt3d,st,pt,t0,piinit%flt1d(-ng+1),sbase,       &
                       gxt,gyt,gzt,time_real,loop,den,      &
                       uinit%flt1d,vinit%flt1d,ugrid,vgrid,ut,vt)
    ENDIF
         

#ifdef TAKON
! Estimate supersaturation. For the most part follows Clark 1973... adapted from Takahashi's code in march.F

! NEED TO INITIALIZE VALUES FIRST!!!!
!          write(0,*) 'solver: ratio,cap,accel,rdry,capw,sheat = ',ratio,cap,accel,rdry,capw,sheat
!          cap   = capw
!          term  = ratio*cap*accel/rdry ! factor for H1 in eq. 22 of Clark 1973
!          work  = cap*cap*ratio/(sheat*rdry) ! Lv^2 / (Rv Cp) where sheat = Cp, ratio/rdry = 0.622/Rd = 1/Rv; cap = Lv


#endif


!==================================================================================================================================
! HOLE filling pos-def variables

   IF ( ihole .eq. 2 ) THEN
   CALL cld_cpu('HOLE-FILLING') 

    CALL HOLEFILL()

   CALL cld_cpu('HOLE-FILLING') 
   ENDIF

!==================================================================================================================================
! PART V:  Copy arrays back

#ifdef MPI
! communicate vertically if needed (for microphysics sedimentation and vertical gradients)

   IF ( nprock > 1 ) THEN

      nampi = 1 + ns
      nb = ng

   downward_tag = 1005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,kmt(-ng+1,-ng+1,-ng+1))

   upward_tag = 1006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,kmt(-ng+1,-ng+1,-ng+1))

   nampi = 1
   nb = ng
   
   downward_tag = 1005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,pt(-ng+1,-ng+1,-ng+1))

   upward_tag = 1006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,pt(-ng+1,-ng+1,-ng+1))

   ENDIF

#endif

   CALL cld_cpu('TIMESTEP-CP2') 
   CALL cld_cpu('TIMESTEP') 

#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

    do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
      u%flt3d(i,j,k)  = ut(i,j,k) 
      v%flt3d(i,j,k)  = vt(i,j,k) 
      w%flt3d(i,j,k)  = wt(i,j,k)
      ! DTD: In the case of the MY scheme, both past and current (updated) time levels for 
      ! pressure (pi) are needed, so postpone this until after the microphysics is called.
      if ( microphys(1:2) /= 'MY' ) then
        pi%flt3d(i,j,k) = pt(i,j,k) 
      endif
      km%flt3d(i,j,k) = kt(i,j,k) 
    enddo; enddo; enddo
#else
    DO k = 1,nz
      u%flt3d(1:nx,1:ny,k)  = ut(1:nx,1:ny,k) 
      v%flt3d(1:nx,1:ny,k)  = vt(1:nx,1:ny,k) 
      w%flt3d(1:nx,1:ny,k)  = wt(1:nx,1:ny,k) 
      ! DTD: In the case of the MY scheme, both past and current (updated) time levels for 
      ! pressure (pi) are needed, so postpone this until after the microphysics is called.
      if ( microphys(1:2) /= 'MY' ) then
        pi%flt3d(1:nx,1:ny,k) = pt(1:nx,1:ny,k) 
      endif
      km%flt3d(1:nx,1:ny,k) = kt(1:nx,1:ny,k) 
    ENDDO 
#endif

! check for variables to convert back from mixing ratio
    DO n = 1,ns

     IF ( s(n)%dyntype .eq. 0 .and. denscale >= 1 .and. microphys(1:4) /= 'THOM' .and. microphys(1:4) /= 'MORR' ) THEN !{
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
!       st(i,j,k,n) = den(k,1)*st(i,j,k,n)
      
      ! DTD: In the case of the MY scheme both past and current (updated) time levels for
      ! scalars are needed so scale both st and sn here instead of scaling and copying back
      ! to sn in one step.
      if ( microphys(1:2) == 'MY' ) then
        s(n)%flt3d(i,j,k) = den(k,1)*s(n)%flt3d(i,j,k)
        st(i,j,k,n) = den(k,1)*st(i,j,k,n)
      else
        s(n)%flt3d(i,j,k) = den(k,1)*st(i,j,k,n)
      endif
      enddo; enddo; enddo
#else
      DO k = 1,nz-1
!      st(1:nx-1,1:ny-1,k,n) = den(k,1)*st(1:nx-1,1:ny-1,k,n)

      ! DTD: In the case of the MY scheme both past and current (updated) time levels for
      ! scalars are needed so scale both st and sn here instead of scaling and copying back
      ! to sn in one step.
      if ( microphys(1:2) == 'MY' ) then
        do j=1,ny-1
          do i=1,nx-1
            s(n)%flt3d(i,j,k) = den(k,1)*s(n)%flt3d(i,j,k)
            st(i,j,k,n) = den(k,1)*st(i,j,k,n)
          end do
        end do
      else ! Carry on as normal for other microphysics schemes.
        s(n)%flt3d(1:nx-1,1:ny-1,k) = den(k,1)*st(1:nx-1,1:ny-1,k,n)
      endif
      ENDDO
#endif

    ELSE ! } {

#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      ! DTD: In the case of the MY scheme both past and current (updated) time levels for
      ! scalars are needed so postpone this copying until later.
      if( microphys(1:2) /= 'MY' ) then
        do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
          s(n)%flt3d(i,j,k) = st(i,j,k,n)
        enddo; enddo; enddo
      endif
#else
      ! DTD: In the case of the MY scheme both past and current (updated) time levels for
      ! scalars are needed so postpone this copying until later.
      if( microphys(1:2) /= 'MY' ) then
        s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1) = st(1:nx-1,1:ny-1,1:nz-1,n)
      endif
#endif
    
    ENDIF ! }
    
    ENDDO !! n

   CALL cld_cpu('TIMESTEP') 
   CALL cld_cpu('TIMESTEP-CP2') 

   IF ( debugsolver ) write(0,*) 'SOLVER3D:  COMPLETED UPDATING DEPDENDENT VARIABLES'

#ifdef MPI
   kzb = 1
   kze = ktile
   if (kzend .eq. nzend) kze = kzend-kzbeg

   jyb = 1
   jye = jtile
   if (jyend .eq. nyend) jye = jyend-jybeg

   ixb = 1
   ixe = itile
   if (ixend .eq. nxend) ixe = ixend-ixbeg

   DO k = kzb,kze
    DO j = jyb,jye
     DO i = ixb,ixe
#else
   DO k = 1,nz-1
    DO j = 1,ny-1
     DO i = 1,nx-1
#endif
       IF ( .not. ( u%flt3d(i,j,k) .lt. 500. .and. u%flt3d(i,j,k) .gt. -500. ) .or.   &
            .not. ( v%flt3d(i,j,k) .lt. 500. .and. v%flt3d(i,j,k) .gt. -500. ) .or.   &
            .not. ( w%flt3d(i,j,k) .lt. 500. .and. w%flt3d(i,j,k) .gt. -500. ) ) THEN
#ifdef MPI
         write(luno,*) 'Winds blowing up! STOP! my_rank=,',my_rank,'i,j,k = ',i,j,k
#else
         write(luno,*) 'Winds blowing up! STOP! i,j,k = ',i,j,k
#endif
         write(0,*) 'Winds blowing up! STOP! my_rank=,',my_rank,'i,j,k = ',i,j,k
         write(luno,*) 'u,v,w = ',u%flt3d(i,j,k),v%flt3d(i,j,k),w%flt3d(i,j,k)
         CALL commasmpi_abort()
         STOP
       ENDIF


     ENDDO
    ENDDO
   ENDDO

!==================================================================================================================================
! Surface physics, if needed.

         IF ( isfcl .eq. -1 ) THEN

         CALL sfcflx(ugrid,vgrid,nx,ny,nz,dx,dy,dt,time,ns,den,pi%flt3d,piinit%flt1d,ut,vt,gxt,gyt,gzt,st)

         ENDIF
    
    IF ( isfcphys >= 1 ) THEN

     IF ( debugsolver ) write(0,*) 'SOLVER: calling sfc_physics'
     CALL cld_cpu('SFC_PHYSICS') 

      call  sfc_interface(gd,               &
                    pi, piinit,             &            ! PI, PIINIT
                    km, kminit,             &            ! KM, KINIT
                    u, v, w,                &
                    s, ns,                  &
                    gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1), & 
                    preciptmp,                 &
                    lat,lon, gtx, gty,      &
                    time,nx,ny,nz,dt,ugrid,vgrid,nst,dx,dy)

     IF ( debugsolver ) write(0,*) 'SOLVER: done with sfc_physics'
   
     CALL cld_cpu('SFC_PHYSICS') 
    
    ENDIF


!=================================================================================================================================
! MPI: FILL GHOST ZONES BEFORE MP

#ifdef MPI

!       write(luno,*) 'Max Min before Microphysics'
!       CALL PRINT_MPI( gd )


#endif

!=================================================================================================================================
! PART VI:  MICROPHYSICS

!mpidebug
  do n = 1,ns
#ifdef MPI
    if( debug_mpi .and. s(n)%name(1:1) == 'Q' .and. my_rank==0) then
!    if(my_rank==0 .and. debug_mpi) then
     write(0,"('--------------SOLVER: before microphysics --------------')")
     write(0,"(4x,'k',8x,'ut',12x,'vt',12x,'wt',12x,'pt',9x,A,8x,A,8x,'fs')") 'temp',trim(s(n)%name)
     i = nx/2 ! 21
     j = ny/2 ! 22
     do k = -ng+1,nz+ng
      if (s(n)%name(1:1) == 'Q') then
      write(0,"(i2,1x,i2,7(2x,es12.5))") &
        my_rank,k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        st(i,j,k,n)*1000,s(n)%flt3d(i,j,k)*1000,fs(i,j,k)
      else
      write(0,"(i2,1x,i2,7(2x,es12.5))") &
        my_rank,k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        st(i,j,k,n),s(n)%flt3d(i,j,k),fs(i,j,k)
      endif
     enddo
    endif
#else
    if( .false. .and. debugsolver .and. s(n)%name(1:1) == 'Q' ) then
     write(0,"('--------------SOLVER: before microphysics --------------')")
     write(0,"(4x,'k',8x,'ut',12x,'vt',12x,'wt',12x,'pt',9x,A,8x,A,8x,'fs')") 'temp',trim(s(n)%name)
     i = nx/2 ! 21
     j = ny/2 ! 22
     do k = -ng+1,nz+ng
      if (s(n)%name(1:1) == 'Q') then
      write(0,"(i2,7(2x,es12.5))") &
        k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        st(i,j,k,n)*1000,s(n)%flt3d(i,j,k)*1000,fs(i,j,k)
      else
      write(0,"(i2,7(2x,es12.5))") &
        k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        st(i,j,k,n),s(n)%flt3d(i,j,k),fs(i,j,k)
      endif
     enddo
    endif
#endif
  end do
!end mpidebug

      IF ( microphys .eq. 'KESSLER'  ) THEN

       CALL cld_cpu('MICROPHYSICS')

        CALL KESSLER(s(lt)%flt3d,                         & ! potential temperature
                     s(lv)%flt3d,                         & ! water vapor mixing ratio
                     s(lc)%flt3d,                         & ! cloud water mixing ratio at t
                     s(lr)%flt3d,                         & ! rain    water mixing ratio at t
                     preciptmp,                          & ! precip array
                     sbase(-ng+1,1), piinit%flt1d,       & ! base state theta and pi
                     gz(3)%flt1d,                        & ! vertical grid spacing
                     dt,                                 & ! time step
                     fu, fv, fw, fp,                     & ! scratch arrays
                     nx, ny, nz)                           ! grid dimensions
        IF ( io_flag ) THEN
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          jyb = 1
          jye = jtile
          if (jyend .eq. nyend) jye = jyend-jybeg

          ixb = 1
          ixe = itile
          if (ixend .eq. nxend) ixe = ixend-ixbeg

          do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
#endif
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(lt)%flt3d(i,j,k)
 
              q(lr)            = s(lr)%flt3d(i,j,k)
             
             IF ( associated (vzf%flt3d) ) THEN
             
             CALL qtovzf(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                               mindbz, microphys, 0, dbz%flt3d(i,j,k), vzf%flt3d(i,j,k))
             ELSE
              dbz%flt3d(i,j,k) = qtodbz(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                                 mindbz, microphys, 0)
             ENDIF
            ENDDO
           ENDDO
          ENDDO
      
        ENDIF
  
        IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST KESSLER'

        CALL cld_cpu('MICROPHYSICS')

      ELSEIF ( microphys .eq. 'DRYC' .or. microphys .eq. 'DRYR' .or. &
               microphys .eq. 'KESSLERQC'  ) THEN

       CALL cld_cpu('MICROPHYSICS')

        IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:   BEFORE DRYC'

      IF ( iuvwadv < 3 ) THEN
      CALL KESSLERQC(s(lt)%flt3d,                         & ! potential temperature
                     s(lv)%flt3d,                         & ! water vapor mixing ratio
                     s(lc)%flt3d,                         & ! cloud water mixing ratio at t
                     preciptmp,                          & ! precip array
                     sbase(-ng+1,1), piinit%flt1d,       & ! base state theta and pi
                     gz(3)%flt1d,                        & ! vertical grid spacing
                     dt,                                 & ! time step
                     fu, fv, fw, fp,                     & ! scratch arrays
                     nx, ny, nz)                           ! grid dimensions
  
        IF ( io_flag ) THEN
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          jyb = 1
          jye = jtile
          if (jyend .eq. nyend) jye = jyend-jybeg

          ixb = 1
          ixe = itile
          if (ixend .eq. nxend) ixe = ixend-ixbeg

          do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
#endif
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(lt)%flt3d(i,j,k)
 
              dbz%flt3d(i,j,k) = mindbz

              IF ( lthetav >= 1 ) xtra(lthetav)%flt3d(i,j,k)  = &
                    gd%xyz3d%flt4d(i,j,k,km%index+lt)*(1. + 0.61*gd%xyz3d%flt4d(i,j,k,km%index+lv) )
              IF ( lairtem >= 1 ) xtra(lairtem)%flt3d(i,j,k)  = &
                   gd%xyz3d%flt4d(i,j,k,km%index+lt)*(pi%flt3d(i,j,k) + piinit%flt1d(k))
              IF ( lairpress >= 1 ) xtra(lairpress)%flt3d(i,j,k)  = &
                    1.e5*(pi%flt3d(i,j,k) + piinit%flt1d(k))**3.509
  
            ENDDO
           ENDDO
          ENDDO
      
        ENDIF

        
        ENDIF
  
        IF ( debugsolver ) write(0,*) my_rank, 'SOLVER3D:   PAST DRYC'

        CALL cld_cpu('MICROPHYSICS')


      ELSEIF ( microphys .eq. 'SCG'  ) THEN
  
       CALL cld_cpu('MICROPHYSICS')
  
        CALL ZHANG  (s(lt)%flt3d,                         & ! potential temperature
                     s(lv)%flt3d,                         & ! water vapor mixing ratio
                     s(lc)%flt3d,                         & ! cloud water mixing ratio at t
                     s(lr)%flt3d,                         & ! rain    water mixing ratio at t
                     preciptmp,                          & ! precip array
                     sbase(-ng+1,1), piinit%flt1d,           & ! base state theta and pi
                     gz(3)%flt1d,                        & ! vertical grid spacing
                     dt,                                 & ! time step
                     fu, fv, fw, fp,                     & ! scratch arrays
                     nx, ny, nz)                           ! grid dimensions
  
        IF ( io_flag ) THEN
  
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
  
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(lt)%flt3d(i,j,k)
 
              q(lr)             = s(lr)%flt3d(i,j,k)
              dbz%flt3d(i,j,k) = qtodbz(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                                 mindbz, microphys, 0)
  
            ENDDO
           ENDDO
          ENDDO
      
        ENDIF
  
        IF ( debugsolver ) print *,'SOLVER3D:   PAST SCG Microphysics' 

        CALL cld_cpu('MICROPHYSICS')

      ELSEIF ( microphys .eq. 'WARMLFO' ) THEN
  
       CALL cld_cpu('MICROPHYSICS')

         CALL WARMLFO(s(lt)%flt3d,                         & ! potential temperature
                      s(lv)%flt3d,                         & ! water vapor mixing ratio
                      s(lc)%flt3d,                         & ! cloud water mixing ratio at t
                      s(lr)%flt3d,                         & ! rain    water mixing ratio at t
                      preciptmp,                          & ! precip array
                      sbase(-ng+1,1), piinit%flt1d,       & ! base state theta and pi
                      gz(3)%flt1d,                        & ! vertical grid spacing
                      dt,                                 & ! time step
                      nx, ny, nz)                           ! grid dimensions

        IF ( io_flag ) THEN
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          jyb = 1
          jye = jtile
          if (jyend .eq. nyend) jye = jyend-jybeg

          ixb = 1
          ixe = itile
          if (ixend .eq. nxend) ixe = ixend-ixbeg

          do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
#endif
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(lt)%flt3d(i,j,k)

              q(lr)             = s(lr)%flt3d(i,j,k)
              dbz%flt3d(i,j,k) = qtodbz(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                                 mindbz, microphys, 0)

            ENDDO
           ENDDO
          ENDDO

        ENDIF

        IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST WARMLFO'

       CALL cld_cpu('MICROPHYSICS')


      ELSEIF ( microphys(1:3) .eq. 'LFO' ) THEN

       CALL cld_cpu('MICROPHYSICS')

       DO n = 1,nphys


        CALL LFO_ICE_DRIVE(piinit%flt1d, sbase(-ng+1,1), sbase(-ng+1,2),         &   ! base state
                           pi%flt3d,                                             &   ! pi'
                           gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),         &     ! scalars
                           preciptmp,                                            &   ! precip
                           dt/float(nphys), gz(3)%flt1d,                         &   ! time step, vertical grid spacing
                           nx, ny, nz, ns, ng )

       ENDDO



       IF ( io_flag ) THEN
#ifdef MPI
       kzb = 1
       kze = ktile
       if (kzend .eq. nzend) kze = kzend-kzbeg
 
       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg
 
       ixb = 1
       ixe = itile
       if (ixend .eq. nxend) ixe = ixend-ixbeg
 
       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else  
         DO k = 1,nz-1
          DO j = 1,ny-1
           DO i = 1,nx-1
#endif
             pii_total = piinit%flt1d(k) + pi%flt3d(i,j,k)
             t_total     = s(lt)%flt3d(i,j,k)

             q(lr) = s(lr)%flt3d(i,j,k)
             q(li) = s(li)%flt3d(i,j,k)
             q(ls) = s(ls)%flt3d(i,j,k)
             q(lh) = s(lh)%flt3d(i,j,k)
             IF ( associated (vzf%flt3d) ) THEN
             
             CALL qtovzf(4, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                               mindbz, microphys, 0, dbz%flt3d(i,j,k), vzf%flt3d(i,j,k))
             ELSE
             dbz%flt3d(i,j,k) = qtodbz(4, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                               mindbz, microphys, 0)
             ENDIF

           IF ( lairtem >= 1 )  THEN
           
             xtra(lairtem)%flt3d(i,j,k)  = s(lt)%flt3d(i,j,k)*(pi%flt3d(i,j,k) + piinit%flt1d(k))
           
           ENDIF
           
           ENDDO
          ENDDO
         ENDDO

       ENDIF

       IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST LFO_ICE_DRIVE'

       CALL cld_cpu('MICROPHYSICS')

      ELSEIF( microphys(1:3) .eq. 'TAK' ) THEN

       IF ( debugsolver ) write(0,*) 'SOLVER3D:   CALL TAK_DRIVE'

       dt1 = dt/float(nphys)


#ifdef TAKON
!       CALL cld_cpu('MICROPHYSICS')
      DO n = 1,nphys
      CALL TAK_DRIVE(gd, st,                                                &
                      piinit%flt1d, sbase,                                  &     ! base state
                      pi%flt3d,                                             &     ! pi'
                       gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),        &     ! scalars
                      preciptmp, precip_old,                                &     ! precip
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,z1d4,                                             &     ! vertical grid height
                       nx, ny, nz, ng, ng, ns,                              &
                      u%flt3d,v%flt3d,w%flt3d,km%flt3d,                     &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                      t0,ut,vt,wt,divv,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d, io_flag, time, time_real, dt,     &
                      elec, cion, muz, vzf,chgadvtemp, gxt, gyt, gzt,                  &
                      tstat, this, lstt, bcx, bcy, tstop,xtra(1)%flt3d)
!                      tstat, lstt, 1, 1)

        ENDDO
!       CALL cld_cpu('MICROPHYSICS')
#else
      IF ( my_rank == 0 ) write(0,*) 'Must run commastak to use TAK microphysics!'
      call commasmpi_abort()
#endif

      IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST HCM_DRIVE'

      ELSEIF( microphys(1:3) .eq. 'HCM' ) THEN

       IF ( debugsolver ) write(0,*) 'SOLVER3D:   CALL HCM_DRIVE'

       dt1 = dt/float(nphys)

#ifdef HCMON
       CALL cld_cpu('MICROPHYSICS')
      DO n = 1,nphys
      CALL HCM_DRIVE(gd,                                                   &
                      piinit%flt1d, sbase,                                  &     ! base state
                      pi%flt3d,                                             &     ! pi'
                       gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),        &     ! scalars
                      preciptmp, precip_old,                                &     ! precip
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,z1d4,                                             &     ! vertical grid height
                       nx, ny, nz, ng, ng, ns,                              &
                      u%flt3d,v%flt3d,w%flt3d,km%flt3d,                     &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                      t0,ut,vt,wt,divv,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d, io_flag, time, time_real, dt,     &
                      elec, cion, muz, vzf, gxt, gyt, gzt,                  &
                      tstat, this, lstt, bcx, bcy, tstop,xtra(1)%flt3d)
!                      tstat, lstt, 1, 1)
       ENDDO
       CALL cld_cpu('MICROPHYSICS')
#else
      IF ( my_rank == 0 ) write(0,*) 'Must run commashcm to use HCM microphysics!'
      call commasmpi_abort()
#endif

      IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST HCM_DRIVE'

      ELSEIF( microphys(1:5) .eq. 'ICE10' .or.        &
              microphys(1:4) .eq. 'ZIEG'  .or.        &
              microphys(1:3) .eq. 'ZVD'   .or.        &
              microphys(1:3) .eq. 'ZMR'   .or.        &
              microphys(1:8) .eq. 'WARMZIEG' ) THEN

       IF ( debugsolver ) write(0,*) 'SOLVER3D:   CALL ICE10_DRIVE'

       dt1 = dt/float(nphys)


       IF ( microphys(1:5) .eq. 'ICE10' ) THEN
       DO n = 1,nphys
#ifdef OLDSTUFF
      CALL ICE10_DRIVE(gd,                                                  &
                       piinit%flt1d, sbase,                                 &     ! base state
                       pi%flt3d,                                            &     ! pi'
                       gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),        &     ! scalars
                       preciptmp, precip_old,                                &     ! precip
                       dt1,                                                 &     ! time step
                       gx(3)%flt1d,                                         &     ! horizontal grid spacing
                       gy(3)%flt1d,                                         &     ! horizontal grid spacing
                       gz(3)%flt1d,                                         &     ! vertical grid spacing
                       gzt,z1d4,                                            &     ! vertical grid height
                       nx, ny, nz, ng, ng, ns,                              &
                       u%flt3d,v%flt3d,w%flt3d,km%flt3d,                    &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                       t0,ut,vt,wt,divv,fu,fv,fw,fp,fs,                       &     ! temporary arrays
                       dx,dy,dz,                                            &     ! average dx,dy,dz
                       nphys,luno,1,1,1,                                    &     ! istag,jstag,kstag
                       dbz%flt3d, io_flag, time, time_real, dt1,    &
                       elec, cion, muz, vzf,chgadvtemp,  gxt, gyt, gzt,                 &
                       tstat, this, lstt, bcx, bcy, tstop,xtra(1)%flt3d)
#endif
         ENDDO

       ELSEIF ( microphys(1:1) .eq. 'Z' ) THEN

        IF ( .false. .and. isedonly /= 0 ) THEN
         write(luno,*) 'Before micro: k,q,c'
         i = km%index
         jy = 1
         ix = nx/2 - 1
         DO kz = nz-1,1,-1
           write(luno,*) kz,gd%xyz3d%flt4d(ix,jy,kz,i+lh),gd%xyz3d%flt4d(ix,jy,kz,i+lnh), &
             gd%xyz3d%flt4d(ix,jy,kz,i+lnh)*den(kz,1) ! &
           !,gd%xyz3d%flt4d(ix,jy,kz,i+lzh),gd%xyz3d%flt4d(ix,jy,kz,i+lzh)/den(kz,1)
         ENDDO
        ENDIF

        IF ( microphys(llen:llen) == 'F' ) THEN

#ifdef USEWRF
        
         CALL cld_cpu('MICROPHYSICS')
          i = km%index

        kzb = -1
        kze = ktile+2
        IF ( kzbeg == nzbeg ) kzb = 1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO kz=kzb,kze
          DO ix=ixb,ixe
            DO jy=jyb,jye
              ut(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/  &
     &         (287.04* gd%xyz3d%flt4d(ix,jy,kz, i+lt))
              fp(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509
              divv(ix,jy,kz) = piinit%flt1d(kz)+pi%flt3d(ix,jy,kz)

              IF ( iturbenhance > 0 ) THEN

                 fs(ix,jy,kz) = tke_diss(ix,jy,kz)

!                 IF ( len_type .eq. 0 .or. len_type .eq. 3 ) dyl = 1./(gxt(ix,3)*gyt(jy,3)*gzt(kz,3))**(1./3.)
!         
!                 IF( len_type .eq. 0 ) mlen = Cm*dyl ! (dx*dy*dz)**(1./3.)
!                 IF( len_type .eq. 1 ) mlen = Cm/gzt(kz,3) ! * dzz(kz)
!                 IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gzt(kz,1)) + 1./Lmax)
!                 IF( len_type .eq. 3 ) mlen = Cm * Lmax
!                 IF( len_type .eq. 4 ) THEN
!                   IF ( gzt(kz,1) .gt. 1400. ) THEN
!                    mlen = Cm*dyl
!                   ELSE
!                    mlen = Cm / (1./(0.4*gzt(kz,1)) + 1./dyl)
!                   ENDIF
!                 ENDIF
!                 fs(ix,jy,kz) =  0.5*Ce*Cm*(km%flt3d(ix,jy,kz)/mlen)**2 
                      
              ENDIF

            ENDDO
          ENDDO
!          write(0,*) 'k,p,dn = ',kz,fp(1,1,kz),ut(1,1,kz)
        ENDDO
        
!        write(0,*) 'solver: my_rank, nxend, ixbeg,ixend,ib,ie,ni = ',my_rank, nxend, ixbeg,ixend,ib,ie,ni 
        

          IF ( microphys(1:5) == 'ZVDHM' ) THEN ! {
            IF ( ipconc .eq. 5 ) THEN !{
            
            DO n = 1,nphys
            call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               qhl = gd%xyz3d%flt4d(ib,jb,kb, i+lhl),       &
                               qsw  = gd%xyz3d%flt4d(ib,jb,kb, i+lsw),       &
                               qhw  = gd%xyz3d%flt4d(ib,jb,kb, i+lhw),       &
                               qhlw  = gd%xyz3d%flt4d(ib,jb,kb, i+lhlw),       &
                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),       &
                               cna  = gd%xyz3d%flt4d(ib,jb,kb, i+lccna),       &
                               f_cna = ( lccna > 1 ),                          &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),       &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               chl = gd%xyz3d%flt4d(ib,jb,kb, i+lnhl),      &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               vhl = gd%xyz3d%flt4d(ib,jb,kb, i+lvhl),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                              !nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
             
             ENDDO
             ENDIF !}

          ELSEIF ( microphys(1:4) == 'ZVDH' ) THEN
            
           DO n = 1,nphys
            IF ( ipconc .eq. 5 ) THEN
             IF ( wrfccn ) THEN
             call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               qhl = gd%xyz3d%flt4d(ib,jb,kb, i+lhl),      &
                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),     &
                               cna  = gd%xyz3d%flt4d(ib,jb,kb, i+lccna),   &
                               f_cna = ( lccna > 1 ),                      &
                               cna_co  = gd%xyz3d%flt4d(ib,jb,kb, i+lccnaco), &
                               f_cnaco = ( lccnaco > 1 ),                     &
                               cna_nu  = gd%xyz3d%flt4d(ib,jb,kb, i+lccnanu), &
                               f_cnanu = ( lccnanu > 1 ),                     &
                               cnuf= gd%xyz3d%flt4d(ib,jb,kb, i+lccnuf),   &
                               f_cnuf = ( lccnuf > 1 ),                    &
                               cn_ac  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_ac),&
                               f_cnac = ( lcn_ac > 1 ),                    &
                               cn_nu  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_nu),&
                               f_cnnu = ( lcn_nu > 1 ),                    &
                               cn_co  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_co),&
                               f_cnco = ( lcn_co > 1 ),                    &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),      &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),      &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               chl = gd%xyz3d%flt4d(ib,jb,kb, i+lnhl),     &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               vhl = gd%xyz3d%flt4d(ib,jb,kb, i+lvhl),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                              !nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
               ELSE
              call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               qhl = gd%xyz3d%flt4d(ib,jb,kb, i+lhl),       &
!                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),       &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),       &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               chl = gd%xyz3d%flt4d(ib,jb,kb, i+lnhl),      &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               vhl = gd%xyz3d%flt4d(ib,jb,kb, i+lvhl),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                              !nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
                
               ENDIF
            ELSEIF ( ipconc .eq. 8 ) THEN
             call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               qhl = gd%xyz3d%flt4d(ib,jb,kb, i+lhl),      &
                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),     &
                               cna  = gd%xyz3d%flt4d(ib,jb,kb, i+lccna),   &
                               f_cna = ( lccna > 1 ),                      &
                               cni  = gd%xyz3d%flt4d(ib,jb,kb, i+lcina),   &
                               f_cina = ( lcina > 1 ),                     &
                               cnuf= gd%xyz3d%flt4d(ib,jb,kb, i+lccnuf),   &
                               f_cnuf = ( lccnuf > 1 ),                    &
                               cn_ac  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_ac),&
                               f_cnac = ( lcn_ac > 1 ),                    &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),      &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),      &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               chl = gd%xyz3d%flt4d(ib,jb,kb, i+lnhl),     &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               vhl = gd%xyz3d%flt4d(ib,jb,kb, i+lvhl),     &
                               zrw = gd%xyz3d%flt4d(ib,jb,kb, i+lzr),      &
                               zhw = gd%xyz3d%flt4d(ib,jb,kb, i+lzh),      &
                               zhl = gd%xyz3d%flt4d(ib,jb,kb, i+lzhl),     &
                               pii = divv,                     & ! full pi
                               p   =  fp,                    & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                    & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                    &
                               itimestep = ntimestep,        &
                              RAIN = preciptmp,              &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                              !nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
            ENDIF
            
            ENDDO ! n

          ELSEIF ( microphys(1:3) == 'ZVD' ) THEN
            
           DO n = 1,nphys !{
            IF ( ipconc .eq. 5 ) THEN
             IF ( wrfccn ) THEN
             !  write(0,*) 'solver: calling wrf zvd'
             call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),     &
                               cna  = gd%xyz3d%flt4d(ib,jb,kb, i+lccna),   &
                               f_cna = ( lccna > 1 ),                      &
                               cnuf= gd%xyz3d%flt4d(ib,jb,kb, i+lccnuf),   &
                               f_cnuf = ( lccnuf > 1 ),                    &
                               cn_ac  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_ac),&
                               f_cnac = ( lcn_ac > 1 ),                    &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),      &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),      &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                              !nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
               ELSE
              call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
!                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),       &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),       &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
                               !nproc = 1, f_thproc=.false.,   &
                             ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
                
               ENDIF
            ELSEIF ( ipconc .eq. 8 ) THEN
             call nssl_2mom_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qh  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               cn  = gd%xyz3d%flt4d(ib,jb,kb, i+lccn),     &
                               cna  = gd%xyz3d%flt4d(ib,jb,kb, i+lccna),   &
                               f_cna = ( lccna > 1 ),                      &
                               cnuf= gd%xyz3d%flt4d(ib,jb,kb, i+lccnuf),   &
                               f_cnuf = ( lccnuf > 1 ),                    &
                               cn_ac  = gd%xyz3d%flt4d(ib,jb,kb, i+lcn_ac),&
                               f_cnac = ( lcn_ac > 1 ),                    &
                               ccw = gd%xyz3d%flt4d(ib,jb,kb, i+lnc),      &
                               crw = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),      &
                               cci = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               csw = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               chw = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
                               vhw = gd%xyz3d%flt4d(ib,jb,kb, i+lvh),      &
                               zrw = gd%xyz3d%flt4d(ib,jb,kb, i+lzr),      &
                               zhw = gd%xyz3d%flt4d(ib,jb,kb, i+lzh),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               w   =  w%flt3d,               &
                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dtp = dt1,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              axtra = xtra(1)%flt3d, nxtra3d = nxtra, &
                              tkediss = fs,                       &
                              dbz = dbz%flt3d,                    &
                              vzf = wt,                           &
!                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
!                              tcond = qbudget(1),            &
!                              tevac = qbudget(2),            &
!                              tevar = qbudget(5),            &
!                              train = qbudget(6),            &
!                              rr    = vt,                  & ! must calculate inside or before
                              diagflag = io_flag,             &
!                              nproc = 1, f_thproc=.false.,   &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
            ENDIF
            
            ENDDO ! n }

!          ELSEIF ( microphys(1:3) == 'ZVD' ) THEN
!            IF ( ipconc .eq. 5 ) THEN
!             CALL nssl_2mom_init(ipctmp=5,mixphase=0,ihvol=-1,eqtset_tmp=0)
!            ELSE
!             write(0,*) 'Must have ipconc = 5 for WRF version of ZVD'
!             call commasmpi_abort()
!            ENDIF
            
          
          ELSE
             write(0,*) 'Bad options for WRF version of ZVD'
             call commasmpi_abort()
          
          ENDIF ! }
          
        
 
       IF (  associated( vzf%flt3d ) .and. io_flag ) THEN
        
        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
              vzf%flt3d(ix,jy,kz)= wt(ix,jy,kz)
            ENDDO
          ENDDO
        ENDDO
        
       ENDIF
        
       IF ( io_flag .and. ( lairpress >= 1 .or. lairtem >= 1 .or. lthetav >= 1 .or. ld0 >= 1 ) ) THEN

        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
            
              IF ( lthetav >= 1 ) xtra(lthetav)%flt3d(ix,jy,kz)  = &
                    gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(1. + 0.61*gd%xyz3d%flt4d(ix,jy,kz,i+lv) )
              IF ( lairtem >= 1 ) xtra(lairtem)%flt3d(ix,jy,kz)  = &
                   gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(pi%flt3d(ix,jy,kz) + piinit%flt1d(kz))
              IF ( lairpress >= 1 ) xtra(lairpress)%flt3d(ix,jy,kz)  = &
                    1.e5*(pi%flt3d(ix,jy,kz) + piinit%flt1d(kz))**3.509

              IF ( ld0 >= 1 .and. isedonly == 1 ) THEN
              !       LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

                IF ( gd%xyz3d%flt4d(ix,jy,kz, i+lr) > 1.e-14 .and. gd%xyz3d%flt4d(ix,jy,kz, i+lnr) > 1.e-8 ) THEN
                  tmp = (pii*1000.*gd%xyz3d%flt4d(ix,jy,kz, i+lnr)/gd%xyz3d%flt4d(ix,jy,kz, i+lr)  )**(1./3.)
!        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
!          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
                  xtra(ld0)%flt3d(ix,jy,kz) = 1000.0*(3.67+alphar)/tmp ! multiply by 1000 to get mm
                ELSE
                  xtra(ld0)%flt3d(ix,jy,kz) = 0.0
                ENDIF
              ENDIF


            ENDDO
          ENDDO
        ENDDO
        
        ENDIF
          
#ifdef MPI
         CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
         CALL cld_cpu('MICROPHYSICS')
#else

          write(0,*) 'Solver not compiled to use WRF version microphysics! Set flag -DUSEWRF'
             call commasmpi_abort()
#endif
       
        ELSE ! microphys(llen:llen) == 'F'


          DO n = 1,nphys

          CALL ZIEG_DRIVE(gd,                                                   &
                      piinit%flt1d, sbase,                                  &     ! base state (pinit, ab)
                      pi%flt3d,                                             &     ! pi' (p2)
                       gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),        &     ! scalars (an array)
                      preciptmp, precip_old,                                &     ! precip (xfalltot, precip_old)
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,z1d4,                                             &     ! vertical grid height
                       nx, ny, nz, ng, ng, ns,                              &
                      u%flt3d,v%flt3d,w%flt3d,km%flt3d,                     &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                      t0,ut,vt,wt,divv,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d, io_flag, time, time_real, dt,              &
                      elec, cion, muz, vzf, chgadvtemp,                     &
                      gxt, gyt, gzt,                                        &
                      tstat, this, lstt, bcx, bcy, tstop,xtra(1)%flt3d)
!                      tstat, lstt, 1, 1)
          ENDDO


       IF ( io_flag .and. ld0 >= 1 ) THEN

        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
            
 
              IF ( ld0 >= 1 .and. isedonly == 1 ) THEN
              !       LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

                IF ( gd%xyz3d%flt4d(ix,jy,kz, km%index+lr) > 1.e-14 .and. gd%xyz3d%flt4d(ix,jy,kz, km%index+lnr) > 1.e-8 ) THEN
                  tmp = (pii*1000.*gd%xyz3d%flt4d(ix,jy,kz, km%index+lnr)/gd%xyz3d%flt4d(ix,jy,kz, km%index+lr)  )**(1./3.)
!        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
!          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
                  xtra(ld0)%flt3d(ix,jy,kz) = 1000.0*(3.67+alphar)/tmp ! multiply by 1000 to get mm
                ELSE
                  xtra(ld0)%flt3d(ix,jy,kz) = 0.0
                ENDIF
              ENDIF


            ENDDO
          ENDDO
        ENDDO
        
        ENDIF

          
         ENDIF ! microphys(llen:llen) == 'F'

        IF ( .false. .and. isedonly /= 0 ) THEN
         write(luno,*) 'After micro: k,q,c'
         i = km%index
         jy = 1
         ix = nx/2 - 1
         DO kz = nz-1,1,-1
           write(luno,*) kz,gd%xyz3d%flt4d(ix,jy,kz,i+lh),gd%xyz3d%flt4d(ix,jy,kz,i+lnh), &
             gd%xyz3d%flt4d(ix,jy,kz,i+lnh)*den(kz,1) ! &
!           ,gd%xyz3d%flt4d(ix,jy,kz,i+lzh),gd%xyz3d%flt4d(ix,jy,kz,i+lzh)/den(kz,1)
         ENDDO
        ENDIF
       
       ELSEIF  ( microphys(1:8) .eq. 'WARMZIEG' ) THEN
      
#ifdef OLDSTUFF
      CALL WARMZIEG_DRIVE(gd,                                               &
                      piinit%flt1d, sbase,                                  &     ! base state
                      pi%flt3d,                                             &     ! pi'
                      gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),        &     ! scalars
                      preciptmp,                                            &     ! precip
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,z1d4,                                             &     ! vertical grid height
                       nx, ny, nz, ng, ng, ns,                              &
                      u%flt3d,v%flt3d,w%flt3d,km%flt3d,                     &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                      t0,ut,vt,wt,divv,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d, io_flag, time, time_real, dt,     &
                      elec, cion, muz, vzf, gxt, gyt, gzt,                  &
                      tstat, this, lstt, bcx, bcy, tstop,xtra(1)%flt3d)
#endif
       ENDIF  !! ( microphys .eq. ICE10 )

      IF ( debugsolver ) write(0,*) 'SOLVER3D:   PAST ICE10_DRIVE'



       ELSEIF ( microphys(1:4) .eq. 'THOM' ) THEN

         CALL cld_cpu('MICROPHYSICS')
          i = km%index

        kzb = -1
        kze = ktile+2
        IF ( kzbeg == nzbeg ) kzb = 1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO kz=kzb,kze
          DO ix=ixb,ixe
            DO jy=jyb,jye
!              ut(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/  &
!     &         (287.04* gd%xyz3d%flt4d(ix,jy,kz, i+lt))
              fp(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509
              divv(ix,jy,kz) = piinit%flt1d(kz)+pi%flt3d(ix,jy,kz)
            ENDDO
          ENDDO
!          write(0,*) 'k,p,dn = ',kz,fp(1,1,kz),ut(1,1,kz)
        ENDDO

        IF ( microphys(1:7) .eq. 'THOM411' ) THEN
             call mp_gt_driver411(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qg  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               nr = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               ni = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
!                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               w  =  w%flt3d,                &
                               dt_in = dt,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              refl_10cm = dbz%flt3d,                    &
                              diagflag = io_flag,             &
                              do_radar_ref = 1,               &
                              has_reqc=0, has_reqi=0, has_reqs=0, &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
        ELSEIF ( microphys(1:7) .eq. 'THOM381' ) THEN
             call mp_gt_driver381(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qg  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               nr = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               ni = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
!                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               w  =  w%flt3d,                &
                               dt_in = dt,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              refl_10cm = dbz%flt3d,                    &
                              diagflag = io_flag,             &
                              do_radar_ref = 1,               &
                              has_reqc=0, has_reqi=0, has_reqs=0, &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
        ELSE
        IF (  associated( vzf%flt3d ) ) THEN
             call mp_gt_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qg  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               nr = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               ni = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
!                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dt_in = dt,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              refl_10cm = dbz%flt3d,                    &
                              VT_DBZ_WT = vzf%flt3d,                    &
                              diagflag = io_flag,             &
                              do_radar_ref = 1,               &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)

        ELSE
        
             call mp_gt_driver(                          &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qg  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               nr = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               ni = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
!                               dn  =  ut,                   & ! must calculate inside or before
                               dz  =  dz3d,                  &
                               dt_in = dt,                     &
                               itimestep = ntimestep,            &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              refl_10cm = dbz%flt3d,                    &
                              diagflag = io_flag,             &
                              do_radar_ref = 1,               &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
        ENDIF
        ENDIF
        
       IF ( io_flag .and. ( lairpress >= 1 .or. lairtem >= 1 .or. lthetav >= 1 ) ) THEN

        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
            
              IF ( lthetav >= 1 ) xtra(lthetav)%flt3d(ix,jy,kz)  = &
                    gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(1. + 0.61*gd%xyz3d%flt4d(ix,jy,kz,i+lv) )
              IF ( lairtem >= 1 ) xtra(lairtem)%flt3d(ix,jy,kz)  = &
                   gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(pi%flt3d(ix,jy,kz) + piinit%flt1d(kz))

              IF ( ld0 >= 1 ) THEN
              !       LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

                IF ( gd%xyz3d%flt4d(ix,jy,kz, i+lr) > 1.e-14 .and. gd%xyz3d%flt4d(ix,jy,kz, i+lnr) > 1.e-8 ) THEN
                  tmp = (pii*997.*gd%xyz3d%flt4d(ix,jy,kz, i+lnr)/gd%xyz3d%flt4d(ix,jy,kz, i+lr)  )**(1./3.)
                  xtra(ld0)%flt3d(ix,jy,kz) = 1000.0*3.67/tmp ! multiply by 1000 to get mm
                ELSE
                  xtra(ld0)%flt3d(ix,jy,kz) = 0.0
                ENDIF
              ENDIF


            ENDDO
          ENDDO
        ENDDO
        
        ENDIF
      
    DO n = 1,ns

     IF ( s(n)%dyntype .eq. 0 .and. denscale >= 1  ) THEN !{
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
        s(n)%flt3d(i,j,k) = den(k,1)*s(n)%flt3d(i,j,k)
      enddo; enddo; enddo
#else
      DO k = 1,nz-1
      ! st(1:nx-1,1:ny-1,k,n) = den(k,1)*st(1:nx-1,1:ny-1,k,n)
       s(n)%flt3d(1:nx-1,1:ny-1,k) = den(k,1)*s(n)%flt3d(1:nx-1,1:ny-1,k)

      ENDDO
#endif
     ENDIF
     ENDDO
     
     
         CALL cld_cpu('MICROPHYSICS')



       ELSEIF ( microphys(1:4) .eq. 'MORR' ) THEN

         CALL cld_cpu('MICROPHYSICS')
          i = km%index

       
        kzb = -1
        kze = ktile+2
        IF ( kzbeg == nzbeg ) kzb = 1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO kz=kzb,kze
          DO ix=ixb,ixe
            DO jy=jyb,jye
!              ut(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/  &
!     &         (287.04* gd%xyz3d%flt4d(ix,jy,kz, i+lt))
              fp(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509
              divv(ix,jy,kz) = piinit%flt1d(kz)+pi%flt3d(ix,jy,kz)
            ENDDO
          ENDDO
!          write(0,*) 'k,p,dn = ',kz,fp(1,1,kz),ut(1,1,kz)
        ENDDO

       
             call MP_MORR_TWO_MOMENT(                          &
                               itimestep = ntimestep,            &
                               th  = gd%xyz3d%flt4d(ib,jb,kb, i+lt),       &
                               qv  = gd%xyz3d%flt4d(ib,jb,kb, i+lv),       &
                               qc  = gd%xyz3d%flt4d(ib,jb,kb, i+lc),       &
                               qr  = gd%xyz3d%flt4d(ib,jb,kb, i+lr),       &
                               qi  = gd%xyz3d%flt4d(ib,jb,kb, i+li),       &
                               qs  = gd%xyz3d%flt4d(ib,jb,kb, i+ls),       &
                               qg  = gd%xyz3d%flt4d(ib,jb,kb, i+lh),       &
                               nr = gd%xyz3d%flt4d(ib,jb,kb, i+lnr),       &
                               ni = gd%xyz3d%flt4d(ib,jb,kb, i+lni),      &
                               ns = gd%xyz3d%flt4d(ib,jb,kb, i+lns),      &
                               ng = gd%xyz3d%flt4d(ib,jb,kb, i+lnh),      &
!                               rho  =  ut,                   & ! must calculate inside or before
                               pii = divv,                    & ! full pi
                               p   =  fp,                   & ! must calculate inside or before
                               dt_in = dt,                     &
                               dz =  dz3d,                  &
                               w  =  w%flt3d,                &
                              RAIN = preciptmp,                   &
                              nrain = nrain,                 &
                              refl_10cm = dbz%flt3d,                    &
                              diagflag = io_flag,             &
                              do_radar_ref = 1,               &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &  
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)

        
       IF ( io_flag .and. ( lairpress >= 1 .or. lairtem >= 1 .or. lthetav >= 1 ) ) THEN

        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
            
              IF ( lthetav >= 1 ) xtra(lthetav)%flt3d(ix,jy,kz)  = &
                    gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(1. + 0.61*gd%xyz3d%flt4d(ix,jy,kz,i+lv) )
              IF ( lairtem >= 1 ) xtra(lairtem)%flt3d(ix,jy,kz)  = &
                   gd%xyz3d%flt4d(ix,jy,kz,i+lt)*(pi%flt3d(ix,jy,kz) + piinit%flt1d(kz))
              
              IF ( ld0 >= 1 ) THEN
              !       LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

                IF ( gd%xyz3d%flt4d(ix,jy,kz, i+lr) > 1.e-14 .and. gd%xyz3d%flt4d(ix,jy,kz, i+lnr) > 1.e-8 ) THEN
                  tmp = (pii*997.*gd%xyz3d%flt4d(ix,jy,kz, i+lnr)/gd%xyz3d%flt4d(ix,jy,kz, i+lr)  )**(1./3.)
                  xtra(ld0)%flt3d(ix,jy,kz) = 1000.0*3.67/tmp ! multiply by 1000 to get mm
                ELSE
                  xtra(ld0)%flt3d(ix,jy,kz) = 0.0
                ENDIF
              ENDIF

              IF ( lswmass >= 1 .or. lswdia >= 1 ) THEN
              !     
              !       LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
              !         RHOSN = 100.
              !         CS = RHOSN*PI/6.
              !         DS = 3.
              !          CONS1=GAMMA(1.+DS)*CS




                IF ( gd%xyz3d%flt4d(ix,jy,kz, i+ls) > 1.e-14 .and. gd%xyz3d%flt4d(ix,jy,kz, i+lns) > 1.e-8 ) THEN
                  tmp = gd%xyz3d%flt4d(ix,jy,kz, i+ls)/gd%xyz3d%flt4d(ix,jy,kz, i+lns) 
                ELSE
                  tmp = 0.0
                ENDIF
                
                IF ( lswmass >= 1 ) THEN
                  xtra(lswmass)%flt3d(ix,jy,kz) = 1.e6*tmp ! multiply by 1000 to get mg
                ENDIF
                
                IF ( lswdia >= 1 ) THEN ! mean mass diameter (units of mm)
                  xtra(lswdia)%flt3d(ix,jy,kz) = 1000.*( tmp/100.*6./pii)**(1./3.)
                ENDIF
                
              ENDIF

            ENDDO
          ENDDO
        ENDDO
        
        ENDIF
      
    DO n = 1,ns

     IF ( s(n)%dyntype .eq. 0 .and. denscale >= 1  ) THEN !{
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
        s(n)%flt3d(i,j,k) = den(k,1)*s(n)%flt3d(i,j,k)
      enddo; enddo; enddo
#else
      DO k = 1,nz-1
      ! st(1:nx-1,1:ny-1,k,n) = den(k,1)*st(1:nx-1,1:ny-1,k,n)
       s(n)%flt3d(1:nx-1,1:ny-1,k) = den(k,1)*s(n)%flt3d(1:nx-1,1:ny-1,k)

      ENDDO
#endif
     ENDIF
     ENDDO
     
     
     
         CALL cld_cpu('MICROPHYSICS')



    
   ELSEIF (microphys(1:2) .eq. 'MY') THEN ! Milbrandt and Yau 1,2,2.5,3-moment microphysics
    
    ! Determine which version of the scheme we are calling and the number of microphysics scalars
    
    !WRITE(0,*) 'SOLVER3D: About to call MY_DRIVE'
    
    IF(microphys(1:3) == 'MY1') THEN
      mscheme = 1
      nscalar = 6
    ELSE IF (microphys(1:5) == 'MY2DA') THEN
      mscheme = 3
      nscalar = 12
    ELSE IF (microphys(1:3) == 'MY2') THEN
      mscheme = 2
      nscalar = 12
    ELSE IF (microphys(1:3) == 'MY3') THEN
      mscheme = 4
      nscalar = 17
    END IF
    
    !WRITE(0,*) 'nx,ny,nz',nx,ny,nz
    
    CALL MY_DRIVE(gd,mscheme,nscalar,                                   &
                  piinit%flt1d,                                         &     ! base state pi
                  pi%flt3d,                                             &     ! pi at t
                  gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),         &     ! scalars at t
                  pt,                                                   &     ! pi at t+1
                  st,                                                   &     ! scalars at t+1
                  preciptmp,                                            &     ! precip
                  dt,gzt,                                               &     ! time step, height                                                                              
                  nx, ny, nz, ng, ng, ns, luno,                         &
                  w%flt3d,                                              &     ! vn (u=1,v=2,w=3)  !!mpidebug: removed(1,1,1) after flt3d
                  t0,ut,                                                &     ! temporary arrays
                  dbz%flt3d,io_flag)                                           ! reflectivity
                  

   ENDIF !! ( microphys )

   ! DTD: swap time levels for pi and scalars here in the case of the MY scheme

   IF (microphys(1:2) == 'MY') THEN
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

    do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
        pi%flt3d(i,j,k) = pt(i,j,k) 
    enddo; enddo; enddo
#else
    DO k = 1,nz
        pi%flt3d(1:nx,1:ny,k) = pt(1:nx,1:ny,k) 
    ENDDO 
#endif

    DO n = 1,ns

#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

     do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
      s(n)%flt3d(i,j,k) = st(i,j,k,n)
     enddo; enddo; enddo

#else
     s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1) = st(1:nx-1,1:ny-1,1:nz-1,n)

#endif

    ENDDO !! n

  END IF ! microphys == MY

     IF ( iuvwadv == 4 .and. microphys(1:2) /= 'MY') THEN ! reset theta and qv

     kzb = -ng+1
     kze = ktile+ng
     if (kzbeg .eq. nzbeg) kzb = 1
     if (kzend .eq. nzend) kze = kzend-kzbeg+1

     jyb = -ng+1
     jye = jtile+ng
     if (jybeg .eq. nybeg) jyb = 1
     if (jyend .eq. nyend) jye = jyend-jybeg+1

     ixb = -ng+1
     ixe = itile+ng
     if (ixbeg .eq. nxbeg) ixb = 1
     if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
        s(lt)%flt3d(i,j,k) = st(i,j,k,lt)
        s(lv)%flt3d(i,j,k) = st(i,j,k,lv)
      enddo; enddo; enddo
     
     ENDIF


!mpidebug
  do n = 1,ns
#ifdef MPI
    if( debug_mpi .and. s(n)%name(1:1) == 'Q' .and. my_rank==0) then
!    if(my_rank==0 .and. debug_mpi) then
     write(0,"('--------------SOLVER: after microphysics --------------')")
     write(0,"(4x,'k',8x,'ut',12x,'vt',12x,'wt',12x,'pt',9x,A,8x,A,8x,'fs')") 'temp',trim(s(n)%name)
     i = nx/2 ! 21
     j = ny/2 ! 22
     do k = -ng+1,nz+ng
      if (s(n)%name(1:1) == 'Q') then
      write(0,"(i2,1x,i2,6(2x,es12.5))") &
        my_rank,k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        s(n)%flt3d(i,j,k)*1000,fs(i,j,k)
      else
      write(0,"(i2,1x,i2,6(2x,es12.5))") &
        my_rank,k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        s(n)%flt3d(i,j,k),fs(i,j,k)
      endif
     enddo
    endif
#else
    if( debugsolver .and. s(n)%name(1:1) == 'Q' ) then
     write(0,"('--------------SOLVER: after microphysics --------------')")
     write(0,"(4x,'k',8x,'ut',12x,'vt',12x,'wt',12x,'pt',9x,A,8x,A,8x,'fs')") 'temp',trim(s(n)%name)
     i = nx/2 ! 21
     j = ny/2 ! 22
     do k = -ng+1,nz+ng
      if (s(n)%name(1:1) == 'Q') then
      write(0,"(i2,6(2x,es12.5))") &
        k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        s(n)%flt3d(i,j,k)*1000,fs(i,j,k)
      else
      write(0,"(i2,6(2x,es12.5))") &
        k,ut(i,j,k),vt(i,j,k),wt(i,j,k),pt(i,j,k),                                &
        s(n)%flt3d(i,j,k),fs(i,j,k)
      endif
     enddo
    endif
#endif
  end do
!end mpidebug


!=================================================================================================================================
! Clean up some stuff...

!   deallocate ( sbase )

       CALL cld_cpu('TIMESTEP')

! MPI: Send 2D data for advection
#ifdef MPI
    IF ( number_of_processes .gt. 1 ) THEN
    
    CALL cld_cpu('MPI-COMM-SFC')
         
        num2d = nprecip+neelec2d
        nb = ng ! number of ghost zones needed to comm. 3 zones for 6th-order 2D Crowley advection
        
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,num2d,ng,ng,0,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,preciptmp)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,num2d,ng,ng,0,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,preciptmp)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,num2d,ng,ng,0,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,preciptmp)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,num2d,ng,ng,0,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,preciptmp)

        ENDIF

    CALL cld_cpu('MPI-COMM-SFC')
    
    ENDIF

#endif

! 2D advection of precip duplicate fields for grid motion
    IF ( ugrid /= 0.0 .or. vgrid /= 0.0 ) THEN
    nst = Nint(time_real/dt) + 1
      DO n = 1,num2d
        IF ( n == iprainacc2 .or. &
             n == iphailacc2 .or. &
             n == iphailfacc2 .or. &
             n == iphailnumacc2 .or. &
             n == iphailfnumacc2 .or. &
             n == iprainaccauto2  .or. &
             n == iprainaccmelt2  .or. &
             n == iprainaccshed2  .or. &
             n == ipfdauto2 .or. &
             n == ipfdmelt2 .or. &
             n == ipfdshed2 .or. &
             n == ipgracc2  .or. &
             n == ipgrnumacc2  .or. &
             n == ipfdacc2 .or. &
             n == ipfdnumacc2 ) THEN
           call ADVECTCRW2D          &
                (nx,ny,nst,             &
                 dt,dx,dy, gtx,gty,   &
                 preciptmp(-ng+1,-ng+1,n), ugrid,vgrid,           &
                 fu)
        ENDIF
      ENDDO
    ENDIF

   preciptot(:) = 0.0d0
   
   IF ( myprock == 1 ) THEN

#ifdef MPI
   kzb = 1
   kze = nprecip+neelec2d

   jyb = 1
   jye = jtile
   if (jyend .eq. nyend) jye = jyend-jybeg

   ixb = 1
   ixe = itile
   if (ixend .eq. nxend) ixe = ixend-ixbeg

   DO j = jyb,jye
    DO i = ixb,ixe
      a = 1./(gx(3)%flt1d(i) * gy(3)%flt1d(j))
     DO k = kzb,kze 
#else
   DO j = 1,ny-1
    DO i = 1,nx-1
      a = 1./(gx(3)%flt1d(i) * gy(3)%flt1d(j))
     DO k = 1,nprecip+neelec2d
#endif
       precip(k)%flt2d(i,j) = preciptmp(i,j,k)
       IF ( k <= 4 ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipwfrat > 0 .and. k == ipwfrat ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipwfacc > 0 .and. k == ipwfacc ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipgracc > 0 .and. k == ipgracc ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipfdacc > 0 .and. k == ipfdacc ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( iphailacc > 0 .and. k == iphailacc ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( iphailfacc > 0 .and. k == iphailfacc ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( iprainaccauto > 0 .and. k == iprainaccauto ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( iprainaccmelt > 0 .and. k == iprainaccmelt ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( iprainaccshed > 0 .and. k == iprainaccshed ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipfdauto > 0 .and. k == ipfdauto ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipfdshed > 0 .and. k == ipfdshed ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       IF ( ipfdmelt > 0 .and. k == ipfdmelt ) preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
     ENDDO
    ENDDO 
   ENDDO
   
   ENDIF

   
#ifdef MPI
  CALL MPI_Allreduce(preciptot, preciptotall, nprecip+neelec2d, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)
  
  preciptot(:) = preciptotall(:)
#endif 
   
   

   IF ( io_flag .and. my_rank == 0 ) THEN
     write(luno,*) 'Total rainfall = ',preciptot(3)
     write(luno,*) 'Total hailfall = ',preciptot(4)
     IF ( iphailacc > 0 .and. lhl > 1 ) write(luno,*) 'Total hail = ',preciptot(iphailacc)
     IF ( iphailfacc > 0 .and. lhl > 1 ) write(luno,*) 'Total hail from FD = ',preciptot(iphailfacc)
     IF ( ipwfacc > 0 ) write(luno,*) 'Total PWF rainfall = ',preciptot(ipwfacc)
     IF ( ipgracc > 0 .and. lh > 1 ) write(luno,*) 'Total graupel = ',preciptot(ipgracc)
     IF ( ipfdacc > 0 .and. lf > 1 ) write(luno,*) 'Total frozen drops = ',preciptot(ipfdacc)
     IF ( iprainaccauto > 0 ) write(luno,*) 'Total AUTO rainfall = ',preciptot(iprainaccauto)
     IF ( iprainaccshed > 0 ) write(luno,*) 'Total SHED rainfall = ',preciptot(iprainaccshed)
     IF ( iprainaccmelt > 0 ) write(luno,*) 'Total MELT rainfall = ',preciptot(iprainaccmelt)
     IF ( ipfdauto > 0 ) write(luno,'(a,1pe12.5)') 'Total FD_AUTO = ',preciptot(ipfdauto)
     IF ( ipfdshed > 0 ) write(luno,'(a,1pe12.5)') 'Total FD_SHED = ',preciptot(ipfdshed)
     IF ( ipfdmelt > 0 ) write(luno,'(a,1pe12.5)') 'Total FD_MELT = ',preciptot(ipfdmelt)
   ENDIF

   IF ( debugsolver ) THEN
#ifdef MPI
   kzb = 1
   kze = ktile
   if (kzend .eq. nzend) kze = kzend-kzbeg

   jyb = 1
   jye = jtile
   if (jyend .eq. nyend) jye = jyend-jybeg

   ixb = 1
   ixe = itile
   if (ixend .eq. nxend) ixe = ixend-ixbeg

   do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
   DO k = 1,nz-1
    DO j = 1,ny-1
     DO i = 1,nx-1
#endif
       IF ( .not. ( u%flt3d(i,j,k) .lt. 500. .and. u%flt3d(i,j,k) .gt. -500. ) .or.   &
            .not. ( v%flt3d(i,j,k) .lt. 500. .and. v%flt3d(i,j,k) .gt. -500. ) .or.   &
            .not. ( w%flt3d(i,j,k) .lt. 500. .and. w%flt3d(i,j,k) .gt. -500. ) .or.   &
            .not. ( ut(i,j,k) .lt. 500. .and. ut(i,j,k) .gt. -500. ) .or.   &
            .not. ( vt(i,j,k) .lt. 500. .and. vt(i,j,k) .gt. -500. ) .or.   &
            .not. ( wt(i,j,k) .lt. 500. .and. wt(i,j,k) .gt. -500. )) THEN
#ifdef MPI
!         write(0,*) 'SOLVER: little after microphysics calls'
         write(luno,*) 'Winds blowing up! STOP! my_rank=',my_rank,'i,j,k = ',i,j,k
#else
         write(luno,*) 'Winds blowing up! STOP! i,j,k = ',i,j,k
#endif
         write(0,*) 'Winds blowing up! STOP! my_rank=',my_rank,'i,j,k = ',i,j,k
         write(luno,*) 'u,v,w = ',u%flt3d(i,j,k),v%flt3d(i,j,k),w%flt3d(i,j,k)
         write(luno,*) 'ut,vt,wt = ',ut(i,j,k),vt(i,j,k),wt(i,j,k)
         STOP
       ENDIF
       
!       IF (  s(lv)%flt3d(i,j,k) .gt. 50.e-3 .or. s(lr)%flt3d(i,j,k) .gt. 50.e-3 ) THEN
!       
!#ifdef MPI
!         write(luno,*) 'Scalar blowing up! STOP! my_rank=',my_rank,'i,j,k = ',i,j,k
!#else
!         write(luno,*) 'Scalar blowing up! STOP! i,j,k = ',i,j,k
!#endif
!         write(luno,*) 'qv,qr = ', s(lv)%flt3d(i,j,k), s(lr)%flt3d(i,j,k)
!         write(luno,*) 'u,v,w,km = ',u%flt3d(i,j,k),v%flt3d(i,j,k),w%flt3d(i,j,k),km%flt3d(i,j,k)
!         STOP
!
!       ENDIF

     ENDDO
    ENDDO
   ENDDO
   
   ENDIF

       CALL cld_cpu('TIMESTEP')

!=================================================================================================================================
! MPI: FILL GHOST ZONES

#ifdef MPI

    IF ( number_of_processes .gt. 1 ) THEN

    CALL cld_cpu('MPICOM-SOLVER')

    nampi = s(1)%index - u%index + ns

   IF ( nproci > 1 ) THEN
   westward_tag = 10001
   CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),westward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index) )

   eastward_tag = 10002
   CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        w_proc(my_rank),e_proc(my_rank),eastward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))
   ENDIF

   IF ( nprocj > 1 ) THEN
   southward_tag = 10003
   CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),southward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   northward_tag = 10004
   CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        n_proc(my_rank),s_proc(my_rank),northward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))
   ENDIF

   IF ( nprock > 1 ) THEN

   downward_tag = 10005
   CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),downward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))

   upward_tag = 10006
   CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,ng,nampi,  &
        d_proc(my_rank),u_proc(my_rank),upward_tag,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,u%index))
   
   ENDIF

    CALL cld_cpu('MPICOM-SOLVER')

    ENDIF
#endif

!
! Do trajectories if needed
!
  IF ( itraj > 0 ) THEN
    CALL cld_cpu('TRAJECTORIES')
!    write(0,*) 'solver1: qh,nh = ',st(20,20,30,lh),st(20,20,30,lnh)
    CALL TRAJ(nx,ny,nz,ns,dt,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),   &
     &        sbase,u%flt3d,v%flt3d,w%flt3d,t0,fu,fv,piinit%flt1d, &
     &        pi%flt3d, km%flt3d, dbz%flt3d, elec, fw,fp,fs,  &
!     &        pi%flt3d, km%flt3d, &
     &        gxt,gyt,gzt,time,time_real,    &
     &        uinit%flt1d,vinit%flt1d,ugrid,vgrid,microphys,dx,dy)
    
!    write(0,*) 'solver2: qh,nh = ',st(20,20,30,lh),st(20,20,30,lnh)
    
    CALL cld_cpu('TRAJECTORIES') 
  ENDIF


!=================================================================================================================================
! Calculate vertical vorticity if it is a write or print time.
 
    IF ( io_flag ) THEN

       CALL cld_cpu('TIMESTEP')

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
!      IF ( jybeg == nybeg ) jyb = 3
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
!      IF ( ixbeg == nxbeg ) ixb = 3
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do k = kzb,kze

       do j = jyb,jye
        jm1 = j-1
        jp1 = j+1
        if (jybeg .eq. nybeg) jm1 = max(j-1,1)
        if (jyend .eq. nyend) jp1 = min(j+1,ny-1)

        do i = ixb,ixe
          im1 = i-1
          ip1 = i+1
          if (ixbeg .eq. nxbeg) im1 = max(i-1,1)
          if (ixend .eq. nxend) ip1 = min(i+1,nx-1)
#else
      DO k = 1,nz-1
       DO j = 1,ny-1

        jm1 = max(j-1,1)
        jp1 = min(j+1,ny-1)

        DO i = 1,nx-1
          im1 = max(i-1,1)
          ip1 = min(i+1,nx-1)
#endif

          wz%flt3d(i,j,k)  = .25 * (wzz(i,j,k,im1,jm1) + wzz(i,jp1,k,im1,j)+wzz(ip1,j,k,i,jm1) + wzz(ip1,jp1,k,i,j))
          
          IF ( nxtra > 0 .and. lhwind > 0 ) THEN
          xtra(lhwind)%flt3d(i,j,k) = Sqrt( (ugrid + 0.5*(u%flt3d(i,j,k) + u%flt3d(i+1,j,k)) )**2 +   &
                                             (vgrid + 0.5*(v%flt3d(i,j,k) + v%flt3d(i,j+1,k)) )**2 )
          ENDIF
          
          IF ( k == 1 .and. idbzcomp > 1 ) THEN
            precip(idbzcomp)%flt2d(i,j) = Maxval( dbz%flt3d(i,j,kzb:kze) )
          ENDIF
          
        ENDDO
       ENDDO
     ENDDO
     
     

       CALL cld_cpu('TIMESTEP')

    ENDIF

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

       IF ( ihwindsfcmx > 0 ) THEN
!        write(0,*) 'SOLVER: ihwindsfcmx,nprecip = ',ihwindsfcmx,nprecip
        DO j = jyb,jye
        DO i = ixb,ixe
          tmp = Sqrt( (ugrid + u%flt3d(i,j,1))**2 + (vgrid + v%flt3d(i,j,1))**2 )
!          write(0,*) 'i,j,tmp = ',i,j,tmp
          precip(ihwindsfcmx)%flt2d(i,j) = Max(tmp, precip(ihwindsfcmx)%flt2d(i,j) )
        ENDDO
        ENDDO
       ENDIF

       
       IF ( iwzsfcmax > 0 .or. iwzsfcmin > 0 ) THEN

      jyb = 1
      IF ( jybeg == nybeg ) jyb = 3
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg-2

      ixb = 1
      IF ( ixbeg == nxbeg ) ixb = 3
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-2

!        write(0,*) 'SOLVER: iwzsfcmax,iwzsfcmin = ',iwzsfcmax,iwzsfcmin

        DO j = jyb,jye
         jm1 = j-1
         if (jybeg .eq. nybeg) jm1 = max(j-1,1)

        DO i = ixb,ixe
          im1 = i-1
          if (ixbeg .eq. nxbeg) im1 = max(i-1,1)

             tmp = wzz(i,j,1,im1,jm1)

         IF ( iwzsfcmax > 0 ) precip(iwzsfcmax)%flt2d(i,j) = Max(tmp, precip(iwzsfcmax)%flt2d(i,j) )
         
         IF ( iwzsfcmin > 0 ) precip(iwzsfcmin)%flt2d(i,j) = Min(tmp, precip(iwzsfcmin)%flt2d(i,j) )
        
        ENDDO
        ENDDO
       
       ENDIF



     IF ( inhom == 4 ) THEN ! hacks on boundaries for nonhomogeneous setup
!
!.... CLZ (3/31/14): check to over-write u-component with base state value at south boundary given normal inflow
!

      IF (jybeg .eq. nybeg ) THEN ! check for south boundary

      kzb = 1
!      kze = 12  ! CAVEAT EMPTOR: hardwired kze value for over-writing u-component in 1,kze
      kze = 25  ! CAVEAT EMPTOR: hardwired kze value for over-writing u-component in 1,kze

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do k = kzb,kze

!.... debug
!      write(6,*) 'k=',k,' u=',u%flt3d(30,1,k),' ub=',uinit%flt1d(k),' ug=',ugrid

!       do j = jyb,jye
        do i = ixb,ixe

        vnorm = 0.5 * (v%flt3d(i,1,k) + v%flt3d(i+1,1,k))
        
        if (vnorm .gt. 0.0) then
        u%flt3d(i,1,k) = uinit%flt1d(k) - ugrid  ! subtract grid motion from base-state u-component
        endif


        ENDDO
!       ENDDO
      ENDDO

      ENDIF




!
!.... CLZ (3/31/14): check to over-write v-component with base state value at east boundary given normal inflow
!

      IF (ixend .eq. nxend ) THEN ! check for east boundary

      kzb = 1
!      kze = 12  ! CAVEAT EMPTOR: hardwired kze value for over-writing u-component in 1,kze
      kze = 25  ! CAVEAT EMPTOR: hardwired kze value for over-writing u-component in 1,kze

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      do k = kzb,kze

!.... debug
!      write(6,*) 'k=',k,' u=',u%flt3d(30,1,k),' ub=',uinit%flt1d(k),' ug=',ugrid

        do j = jyb,jye
!        do i = ixb,ixe

        unorm = 0.5 * (u%flt3d(nx,j,k) + u%flt3d(nx,j+1,k))

        if (unorm .lt. 0.0) then
        v%flt3d(nx,j,k) = vinit%flt1d(k) - vgrid  ! subtract grid motion from base-state v-component
        endif


        ENDDO
!       ENDDO
      ENDDO

      ENDIF
      
      ENDIF ! inhom == 4





   IF ( debugsolver ) write(0,*) 'SOLVER3D:    EXITING...'

  RETURN
    
!-----------------------------------------------------------------------------
!
!   /////////////////////          END          \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE SOLVER   ////////////////////
!     
!-----------------------------------------------------------------------------

  CONTAINS
  
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE HOLEFILL    ////////////////////
!     
!-----------------------------------------------------------------------------
  SUBROUTINE HOLEFILL()

   USE COMMASMPI_MODULE
   
   implicit none

#ifdef MPI
  INCLUDE "mpif.h"
#endif
   
   double precision :: tmp1, tmp2, tmpsum, tmpfrac, tmp1all, tmp2all
   integer count

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      real dentmp(nz)
      integer istop
    
!  fv = grid box volume * air density

   istop = 0
#ifdef MPI
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


   do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
!       fv(i,j,k) = den(k,1)/(gx(3)%flt1d(i) * gy(3)%flt1d(j) * gz(3)%flt1d(k))
       fv(i,j,k) = 1./(gx(3)%flt1d(i) * gy(3)%flt1d(j) * gz(3)%flt1d(k))

            IF ( debugsolver .and. .not. ( st(i,j,k,lt) .lt. 900. .and. st(i,j,k,lt) .gt. 100.) ) THEN
              write(luno,*) 'THETA blowing up! STOP! i,j,k = ',i,j,k
              write(luno,*) 'TH = ', st(i,j,k,lt), s(lt)%flt3d(i,j,k), fs(i,j,k)
              write(luno,*) 'THm1= ', st(i,j,k-1,lt), s(lt)%flt3d(i,j,k-1), fs(i,j,k-1)

              istop = istop + 1
              IF( istop > 10 ) call commasmpi_abort()
            ENDIF
      ENDDO
     ENDDO
    ENDDO
    
    


    DO n = 1,ns

     IF ( debugsolver ) THEN
      IF ( s(n)%name(1:1) .eq. 'Q'  ) THEN
       DO k = 1,nz-1
         DO j = 1,ny-1
          DO i = 1,nx-1

       
            IF ( .not. ( st(i,j,k,n) .lt. 50.e-3 .and. st(i,j,k,n) .gt. -1.e-3) ) THEN
       
              write(luno,*) 'Scalar blowing up! STOP! n,i,j,k = ',n,i,j,k
              write(luno,*) s(n)%name,' = ', st(i,j,k,n)
              write(luno,*) 'qv,qr = ', st(i,j,k,lv), st(i,j,k,lr)
!              write(luno,*) 'u,v,w,km = ',u%flt3d(i,j,k),v%flt3d(i,j,k),w%flt3d(i,j,k),km%flt3d(i,j,k)

              STOP
       
            ENDIF

          ENDDO                   
         ENDDO
        ENDDO
   
       ENDIF ! Q
      ENDIF ! debugsolver
     
     IF ( (s(n)%name(1:1) .eq. 'Q' .or. s(n)%name(1:1) .eq. 'V'  .or.   &
           s(n)%name(1:6) .eq. 'CCHAFF' .or.                                 &
           s(n)%name(1:1) .eq. 'C' .or.                                 &
           s(n)%name(1:1) .eq. 'Z') .and.                               &
            s(n)%name(1:2) .ne. 'QV' .and. s(n)%pdef .ge. 1 ) THEN


        IF ( s(n)%dyntype == 1 .or. ( s(n)%dyntype == 0 .and. denscale >= 1 ) ) THEN
           dentmp(1:nz) = den(1:nz,1)
        ELSE
           dentmp(1:nz) = 1.0
        ENDIF
        
        tmp1 = 0.0d0
        tmp2 = 0.0d0
#ifdef MPI
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

        do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
        DO k = 1,nz-1
         DO j = 1,ny-1
          DO i = 1,nx-1
#endif
            tmp1 = tmp1 + Min( 0.0, st(i,j,k,n) * fv(i,j,k) * dentmp(k) )
            tmp2 = tmp2 + Max( 0.0, st(i,j,k,n) * fv(i,j,k) * dentmp(k) )
          ENDDO
         ENDDO
        ENDDO

#ifdef MPI


    mpitotindp(1)  = tmp1
    mpitotindp(2)  = tmp2

!   write(luno,*) 'tile tmp1,tmp2 = ',tmp1,tmp2,my_rank
  CALL MPI_Allreduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)
  
  tmp1 = mpitotoutdp(1)
  tmp2 = mpitotoutdp(2)

  IF ( debugsolver )    write(luno,*) 'var, tmp1,tmp2 = ',s(n)%name,tmp1,tmp2,my_rank

#endif 
      
        IF ( tmp1 .lt. 0.0d0 .and. tmp2 .gt. Abs(tmp1) ) THEN
          tmpfrac = (tmp2 + tmp1)/tmp2
          IF ( tmpfrac .lt. 0.7 ) THEN
           IF ( my_rank == 0 ) THEN
           write(0,*) 'PROBLEM with ',s(n)%name ,' tmpfrac,tmp1,tmp2 = ',tmpfrac,tmp1,tmp2
           write(0,*) 'my_rank = ',my_rank
           write(luno,*) 'PROBLEM with ',s(n)%name ,' tmpfrac,tmp1,tmp2 = ',tmpfrac,tmp1,tmp2
           write(luno,*) 'my_rank = ',my_rank
           ENDIF
           IF ( tmpfrac .lt. 0.55 .and. tmp2 .gt. 0.1 ) THEN
#ifdef MPI
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

           do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
           DO k = 1,nz-1
            DO j = 1,ny-1
             DO i = 1,nx-1
#endif
              IF ( st(i,j,k,n) .ne. 0.0 ) write(luno,*) 'i,j,k,st = ',i,j,k,st(i,j,k,n)
             ENDDO
            ENDDO
           ENDDO
          
            CALL commasmpi_abort()
            STOP
           ENDIF
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
          st(i,j,k,n)  = Max( 0.0d0, st(i,j,k,n)*tmpfrac )
        enddo ; enddo ; enddo
#else
          st(1:nx-1,1:ny-1,1:nz-1,n)  = Max( 0.0d0, st(1:nx-1,1:ny-1,1:nz-1,n)*tmpfrac )
#endif
!         write(*,*) time,': deficit in ', s(n)%name,' is ',0.001*tmp1,0.001*tmp2, Abs(tmp1)/tmp2*100
        ENDIF
      
      ENDIF


    ENDDO
    
    RETURN

    END SUBROUTINE HOLEFILL
!-----------------------------------------------------------------------------
!
!   /////////////////////           END            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE HOLEFILL    ////////////////////
!     
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE DTREND      ////////////////////
!     
!-----------------------------------------------------------------------------
!
! DTREND removes the mean of the perturbation pressure field
!------------------------------------------------------------------------------
!
! Created By Louis Wicker: 8/21/89
! Latest update: 01/02/01
!
!------------------------------------------------------------------------------
  SUBROUTINE DTREND()


      IMPLICIT NONE

#ifdef MPI
  INCLUDE "mpif.h"
#endif

      integer i, j, k
      double precision :: pavg, denom

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
    

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

      pavg = 0.0

      DO k = kzb,kze
       DO j = jyb,jye
        DO i = ixb,ixe
         pavg = pavg + pi%flt3d(i,j,k)
        ENDDO
       ENDDO
      ENDDO

#ifdef MPI

    mpitotindp(1)  = pavg

!   write(luno,*) 'tile pavg = ',my_rank,pavg
  CALL MPI_Allreduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)
  
  pavg = mpitotoutdp(1)

#endif 

      denom = (nxend-1)*(nyend-1)*(nzend-1)
      pavg = pavg/denom

   kzb = 1-ng
   kze = ktile+ng
   if (kzbeg .eq. nzbeg) kzb = 1
   if (kzend .eq. nzend) kze = kzend-kzbeg

   jyb = 1-ng
   jye = jtile+ng
   if (jybeg .eq. nybeg) jyb = 1
   if (jyend .eq. nyend) jye = jyend-jybeg

   ixb = 1-ng
   ixe = itile+ng
   if (ixbeg .eq. nxbeg) ixb = 1
   if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO k = kzb,kze
       DO j = jyb,jye
        DO i = ixb,ixe
         pi%flt3d(i,j,k) = pi%flt3d(i,j,k) - pavg
        ENDDO
       ENDDO
      ENDDO

    RETURN

    END SUBROUTINE DTREND
!-----------------------------------------------------------------------------
!
!   /////////////////////           END            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE DTREND      ////////////////////
!     
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE W_DAMP      ////////////////////
!     
!-----------------------------------------------------------------------------
!
! W_DAMP calculates vertical CFL conditions and in future will optionally apply damping to w
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
  SUBROUTINE W_DAMP(fw)


      IMPLICIT NONE

#ifdef MPI
  INCLUDE "mpif.h"
#endif
      integer i, j, k, maxi, maxj, maxk, km1
      real :: max_vert_cfl,max_vert_cflw,vert_cfl,vert_cflw,maxdub,maxdubw, maxdz
      INTEGER :: some,some1

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      real :: fw (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

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

      max_vert_cfl = 0.0
      max_vert_cflw = 0.0
      some = 0
      some1 = 0

      DO k = kzb,kze
       km1 = Max(1,k-1)
       DO j = jyb,jye
        DO i = ixb,ixe
          vert_cfl = dt*wt(i,j,k)*gzt(k,3)
          
          IF ( vert_cfl > max_vert_cfl ) THEN
            max_vert_cfl = vert_cfl
          ENDIF
          
          IF ( vert_cfl > 1.0 ) some1 = some1 + 1
          IF ( vert_cfl > 2.0 ) some  = some  + 1

          vert_cflw = dt*0.5*(wt(i,j,k) + wt(i,j,km1))*gzt(k,4)
          IF ( vert_cflw > max_vert_cflw ) THEN
            max_vert_cflw = vert_cflw
          ENDIF
          
        ENDDO
       ENDDO
      ENDDO


#ifdef MPI

    mpitotindp(1)  = max_vert_cfl
    mpitotindp(2)  = max_vert_cflw

!   write(luno,*) 'tile pavg = ',my_rank,pavg
  CALL MPI_Allreduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)
  
  max_vert_cfl = mpitotoutdp(1)
  max_vert_cflw = mpitotoutdp(2)


    mpitotindp(1)  = some
    mpitotindp(2)  = some1

!   write(luno,*) 'tile pavg = ',my_rank,pavg
  CALL MPI_Allreduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)
  
  some = mpitotoutdp(1)
  some1 = mpitotoutdp(2)

#endif 

    IF ( my_rank == 0 ) THEN
      IF ( max_vert_cfl > 1. .or. max_vert_cflw > 1. ) THEN
        write(luno,'(a,2(f9.2,1x),3(i8,1x))') 'max_vert_cfl,max_vert_cflw,some1,some,loop = ', &
              max_vert_cfl,max_vert_cflw,some1,some,loop
      ELSEIF ( loop == RKSCHEME ) THEN
        write(luno,'(a,2(f9.2,1x),3(i8,1x))') 'max_vert_cfl,max_vert_cflw,some1,some,loop = ', &
              max_vert_cfl,max_vert_cflw,some1,some,loop
      ENDIF
    ENDIF

    END SUBROUTINE W_DAMP
!-----------------------------------------------------------------------------
!
!   /////////////////////           END            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE W_DAMP      ////////////////////
!     
!-----------------------------------------------------------------------------



   END SUBROUTINE SOLVER
