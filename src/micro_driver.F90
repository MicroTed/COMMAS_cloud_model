#ifdef CM1 
!/* special setup for CM1, but is obsolete as it now uses the WRF module */

#undef CHGELEC
#define ICE3
#define SWM

#else

#include "sam.def.h"
! #define CHGELEC

#endif

#ifdef LFO3
#define ICE3
#undef CHGELEC
#endif

#ifdef NOELEC
#undef CHGELEC
#endif

#ifndef RKIND
#define RKIND 4
#endif

! only use this for 10-ice with SWM for now
!#ifdef SWM
!#undef ICE3
!#endif
! a little driver for calling electrified microphysics.
!
! 'driver.sgi.bg'  'bg' = 'below ground':  set up so that the ground
!   potential is at -dz/2 instead of +dz/2.  To do this, some
!   array equivalencing is implemented.
!  set value of ibg = +1 for above-ground boundary,
!  set value of ibg = -1 for below-ground boundary
!
! 03.2007    (APS) Added print statements for chosen microphysical process 
!
! 9.2005     Added a bunch of vars to the v5d output loop. 
!            Moved frozen drops from hail family to graupel family for 
!            both vis5d and bulk reporting
!
!
! 2.14.2004  Changed energy calculation to the rho*phi formula, replacing
!            the E^2 formula.
!
! 11.17.2003 ion.03.v5d.f: provisions for on-the-fly vis5d output
!
!
! 10.23.2002 version ion.2a.f:
!            set potential buffer zone to be at least 4km deep or 8 gridpoints,
!            whichever is larger instead of just 8 gridpoints
!
!  5.10.2002  New version (.mud.2.f) adding option to use mudpack to solve 
!             for electric potential if istrz .ne. 0
!
!  1/8/2002  .2.  Convert to all local or passed in arrays (i.e. get rid of sam.dims.h)
!
! 1/4/2002  .1f.  Set up array bounds for SLM lightning so they'll be automatically 
!            generated
!
! 12/14/2001 version sam.08.driver.lm.omp.1e.f
!            This revision has the 'elec' array passed in so that it doesn't require
!            changing the domain values here.  (Also pass in arrays needed by
!            the microphysics)
!        
!
!
! 11/26/2001 New version (sam.08.driver.lm.omp.1d.f) to allow compatability
!            with the combined regular and stand-alone versions of the lightning
!            driver *sll* and the lightning solver.
!
!
! 9/12/2001
!       Increase height of Poisson solver domain. (already are extending
!          the lateral boundaries)
!
! 8/17/2001  Add reporting of updraft mass flux, total graupel/hail/ice crystal mass,
!            graupel volume (> 0.5 g/kg)
!            Maybe also upward ice crystal mass flux?
!
! 7/23/2001  Test dividing charge variables by air density before going to
!            advection, then multiply back again before microphysics. 
!            Result:  Much improved conservation of net charge by advection/diffusion!
! 7/21/2001  version sam.08.driver.lowmem.1a.f
!            Start keeping track of Wilson ion damping on rain (reported in new 
!                array xfall)
!            Added report of charge falling to ground (which includes the wilson
!                damping charge).  
!            Added report of total charge after lightning cycle as a check
!                before going to advection/diffusion.
!            
! 5/99     Added monolayer corona charge.  This just limits Ez at surface
!          and doesn't do much else (not advected).  Decay time constant
!          is about 1.5 minutes.
!
! 7/11/99  Extend domain for FFT solver to get a better solution
!          (buffer region set to zero charge)
! 7/14/99  setup so that fields are always dumped when specified
!          by nhelec without changing the delay for lightning dumps
!
! 9/2/99  Stop saving the elec array (in history dumps).  From now on just use the net
!         charge density to recalculate potential and e-field components. 
!
! 10/?/99 Added calculation of total electric field energy and the change
!         after each flash. 
!
! 11/06/99 Added calculations of total net positive and negative charge
!          in the storm and in the magic corona layer
!
! 8/2000  Fixed a bug for when ilight > 0 but nliter = 0 (i.e. charging is 
!         active but lightning turned off)
!


#ifdef WARM
#undef ICE3
#undef ICE10
       subroutine WARMZIEG_DRIVE(gd,                        &
     &                       pinit, ab,                     & ! base state
     &                       p2,                            & ! pi_prime
     &                       an,                            &  ! scalars
     &                       xfalltot,                      &  ! precip
     &                       dtp1, gx,gy,gz,z1d,z1d4,       &  ! time step, vertical grid spacing
     &                       nx, ny, nz, nor, norz, na,     &
     &                       u,v,w,km,                      &  ! vn (u=1,v=2,w=3)
     &                       t0,t1,t2,t3,t4,t5,t6,t7,t8,t9, & ! temporary arrays
     &                       dx,dy,dz,                      & ! nstep, start,tstat,
     &                       ntmul,iunit,istag,jstag,kstag, &
     &                       dbz, io_flag, time, time_real, dt, &
     &                       elec, cion, muz, vzf,             &
     &                       gxt, gyt, gzt,                    & ! electricity arrays
     &                       tstat, this, lstt, bcx1, bcy1, tstop, axtra ) ! axtra=extra scalars
#else
#ifdef ICE3
#undef ICE10
#ifdef LFO3
       subroutine LFO3_DRIVE(                              &
#else
       subroutine ZIEG_DRIVE(                              &
#endif
#ifndef CM1
     &                    gd,                              &
#endif
     &                    pinit, ab,                       & ! base state
#endif
#if defined (HCMON)
#undef ICE10
       subroutine HCM_DRIVE(                               &
     &                    gd,                              &
     &                    pinit, ab,                       & ! base state
#endif
#if defined (TAKON)
#undef ICE10
       subroutine TAK_DRIVE(                               &
     &                    gd, st,                          &
     &                    pinit, ab,                       & ! base state
#endif
#ifdef ICE10
       subroutine ICE10_DRIVE(gd,                          &
     &                        pinit, ab,                   &  ! base state
#endif
     &                        p2,                          & ! pi_prime
     &                        an,                          &  ! scalars
     &                        xfalltot, precip_old,        &  ! precip
     &                       dtp1, gx,gy,gz,z1d,z1d4,     &  ! time step, vertical grid spacing
     &                       nx, ny, nz, nor, norz, na,    &
     &                       u,v,w,km,                      &  ! vn (u=1,v=2,w=3)
     &                       t0,t1,t2,t3,t4,t5,t6,t7,t8,t9, & ! temporary arrays
     &                       dx,dy,dz,                      & ! nstep, start,tstat,
     &                       ntmul,iunit,istag,jstag,kstag, &
     &                       dbz, io_flag, time, time_real, dt, &
#ifndef CM1
     &                       elec, cion, muz, vzf,chgadvtemp,    &
#else
     &                       pn, pb, dn, db,                   &
#endif
     &                       gxt, gyt, gzt,                    & ! electricity arrays
     &                       tstat, this, lstt, bcx1, bcy1, tstop, axtra ) ! axtra=extra scalars
#endif

      USE FILE_MODULE
#ifndef CM1
      USE VIS5D_MODULE
      USE GRID_MODULE
#endif
      USE MICRO_MODULE
#ifdef CHGELEC
      USE ELEC_MODULE
      USE ELEC_STUFF
#endif
      USE INDEX_MODULE
#ifndef CM1
      USE CPUTIME_MODULE
      USE TRMM_MODULE, only : lpredict,ibsd,jbsd,iesd,jesd
#endif
      USE COMMASMPI_MODULE
      USE PARAM_MODULE, only : g, pi=>pii

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"

! #ifdef BOXMG
! #include        "BMG_constants_f77.h"
!      INCLUDE   'BMG_workspace_f77.h'
!      INCLUDE   'BMG_parameters_f77.h'
! #endif

#endif

      integer nx,ny,nz,na,nor,norz,nba,nv,ia,il

#ifndef CM1
      TYPE(GRID)         :: gd 
      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      TYPE(VARIABLE)     :: vzf
#endif

      integer, parameter    :: ng1 = 1

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real               :: km(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-norz+ng1:nz+norz,4)

      integer :: tstat, this
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

!      include 'sam.dims.h'

      
      integer iunit

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real chgadvtemp(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,2) ! This is only a full size array if CHGELEC > 0, otherwise is (1,1,1,1)

      real pinit(-norz+ng1:nz+norz)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi
      real pb(-norz+ng1:nz+norz)
      real db(-norz+ng1:nz+norz)
      real ab(-norz+ng1:nz+norz,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
#ifdef TAKON
      real st(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
#endif
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)
!      real an(nx,ny,nz,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      logical :: io_flag, io_flag_tmp

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)
      real precip_old(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip)


      integer ntmul
      real dtp1

#ifdef CM1
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

#else
      real, allocatable :: dn(:,:,:), pn(:,:,:)
#endif
!
! Maximum supersaturation for a parcel
!
      
      
      integer imicro
      logical iflag
      
!      integer lqi  !  index of first mixing ratio (qc)
!      parameter (lqi = 4)
!      integer lscwi
!      parameter ( lscwi = 17 )

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt


      save nstep
      data nstep/0/
      parameter ( nv = 3 )

      real dzc(nz)                     ! 1/dz(k)
      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-norz+ng1:nz+norz)
      real z1d(-norz+ng1:nz+norz,4)
      real z1d4(nzend,4)

      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
!     integer nxl,nyl,nzl
      integer nx1,ny1,nz1

      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n,k0
      real x,tx
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

!
!  Vis5D output stuff
!

      integer v5dwriteappend
      integer myv5dupdate

      integer idatime,hh,mm,ss
      
      real ccimx
      parameter ( ccimx = 200.0 )
      
      integer istat1

      logical ifv5d
!      integer nzdbz
#ifndef CM1
      data ifv5d/.true./
#else
      data ifv5d/.false./
#endif

#ifdef CM1
      logical, parameter :: lpredict = .true.
      integer            :: ibsd,iesd,jbsd,jesd
      logical            :: lv5dwrite,lncwrite,lrst,lprt
      
#endif
!      data idoradar /.false./
#ifndef CM1
      integer iv5dinterval,iv5dstart
      data iv5dinterval / 2000 /
      
      save ifv5d,iv5dinterval
!      save iv5dcnt
      integer nr, nr0, nr1, nr2, nc !, nl(MAXVARS)

      integer varnum, it, iv
      character*10 varname(MAXVARS)
!      integer idates(MAXTIMES)
!      integer itimes(MAXTIMES)
      integer itimes0,idates0
!      integer compressmode
!      integer projection
      real proj_args(100)
!      integer vertical
      real vert_args(MAXLEVELS)


!     initialize the variables to missing values
      data nr,nc / IMISSING, IMISSING /
!      data (nl(i),i=1,MAXVARS) / MAXVARS*IMISSING /
      data (varname(i),i=1,MAXVARS) / MAXVARS*"          " /
!      data (idates(i),i=1,MAXTIMES) / MAXTIMES*IMISSING /
!      data (itimes(i),i=1,MAXTIMES) / MAXTIMES*IMISSING /
      data (proj_args(i),i=1,100) / 100*MISSING /
      data (vert_args(i),i=1,100) / 100*MISSING /
      
!      data compressmode / 1 /
!      data projection / 0 /
!      data vertical / 2 /
      data iv5dstart / 1 /
      
      save varname,varnum,it,iv5dstart
      save proj_args,vert_args
#endif
      real :: smin, smax, fac

!
!
!      
      logical lnstat
!      real dslightz   ! lightning resolution
!      save dslightz
!      integer nxslm,nyslm,nzslm    !  dimensions for the SLM lightning
!      integer nxb2,nyb2,nzbt2,nzb2,nb2
!      integer nxb2s(2),nyb2s(2),nzbt2s(2),nzb2s(2),nb2s(2)

!      save nxslm,nyslm,nzslm    !  dimensions for the SLM lightning
!      save nxb2s,nyb2s,nzbt2s,nzb2s,nb2s
!      save nxb2,nyb2,nzbt2,nzb2,nb2
      
      integer icldtop
      real    zcldtop
      integer, allocatable, save :: kcldtop(:,:)
      real,    allocatable, save :: thproc(:,:),thproctot(:,:)

#define NUMPROCTMAX 101
! 41 and add 3 for time-height shed rates from graupel, FD, and hail
#define NUMPROCT 44
      
      integer, save :: numproc = NUMPROCT
      integer, parameter :: ntq1    = 4
      
      real, allocatable, save :: proctot(:) ! (numproc)
      character(len=40), allocatable, save :: procnames(:)
      
      real tmp, tmpg, tmph, tmps, tmpi, tmphnum, tmphfnum
      double precision :: tmppwfs, tmppwf, tmppwfh, tmppwff, tmppwfhl
      double precision :: tmpdp, cwch
!      real, save :: cwmasn,cwmasx


      

!      real, allocatable, save :: tq1(:,:,:)  !mp process temporary
      real, allocatable :: tq1(:,:,:,:)  !mp process temporary



!
! Local temp array
!

!      real scrnchg(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)

       real,allocatable :: tt0(:,:,:), tt7(:,:,:)
       real,allocatable :: tem1(:,:,:),tem2(:,:,:),tem3(:,:,:)
       real,allocatable :: ssat(:,:,:), ssfilt(:,:,:), ssati(:,:,:)
!      real tt0(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real tt7(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!
!
      integer iwrite,iwritecg

!
!
!  electricity
!

      
      
!
!      integer ntem
!      parameter(ntem=200)
      real, allocatable,save :: temc(:)
      real, allocatable :: xfall(:,:,:) ! (nx,ny,na). xfall is q*vt*dt/dz(1), 
                         ! which is q*(fraction of cell landing on the ground)
                         ! Note that the fraction could be > 1 if substepped because vt*dt > dz(1)
                         ! Mass rate (per m^2) would be xfall*rho0*dz/dt
                         ! Accumulation (per m^2) is xfall*dz*rho0 (per time step)
                         ! Total mass landing in grid cell is xfall*dz*dx*dy*rho0 
                         ! Total number is xfall*dz*dx*dy

      integer :: nrate2d = 3
      real, allocatable :: rate2d(:,:,:) ! (nx,ny,na) For 2D (column) microphysics rates 
      real, allocatable, save :: rate2dtot(:)

     
      integer initcond
!      integer :: firstsolve = 0
      data initcond /0/
      save initcond
      real zlayer
!
! Miscellaneous Variables      
!
      integer ifile
      integer ix,jy,kz,kk,i1,j1,k1,i2,j2,k2
      integer iliter,lgtstp,numlgt,loccur
      character(LEN=80) fnelec,fncld,flight,flgtadj,fethr,fnscd,fnpre
      character(LEN=6)  rstime
!      character*5  rstime5
      integer irstime
      character(LEN=3)  nmliter

      real times
      integer ierr,nyy,nxx
      integer idelay ,idelaycg ! count min timesteps until data dump can happen
      data idelay /60/
      data idelaycg /60/
      save idelay,idelaycg
      logical lhelec, lnsave, ldumpsp
      
      real wmax,wmin,vmax,tmin,ssmax,ssatmax,condmax, depmax
      real :: condmaxmax = 0.0, depmaxmax = 0.0
      
      double precision udv5n, udv5, udv10, udv20
      real dbzmx, dbzi


      real timetd1,timetd2
      real tarray(2),etime
!
! MP Process Max/Min/Average
!
      integer iavgcnt,iqproc
      
      double precision qprocmax, qprocmin, qprocavg
      real,allocatable,save,dimension(:) :: qproczmx, qproczmn, qproczav
      real qprocthr
      parameter (qprocthr = 1.e-06)

!
!  other test quantities
!
      double precision,allocatable,save :: udmf(:)  ! updraft mass flux (w > 3.0)
      double precision,allocatable,save :: udmf1(:)  ! updraft mass flux (w > 0.5)
      double precision,allocatable,save :: udmf10(:)  ! updraft mass flux (w > 10.0)
      double precision,allocatable,save :: uicmf(:) ! upward ice crystal mass flux
      double precision,allocatable,save :: grvol(:) ! graupel volume
      double precision,allocatable,save :: grvol1(:) ! graupel volume
      double precision,allocatable,save :: hlvol(:) ! hail volume
      double precision,allocatable,save :: hlvol1(:) ! hail volume
      double precision,allocatable,save :: grms(:) ! graupel mass
      double precision,allocatable,save :: hlms(:) ! hail mass

      double precision,allocatable,save :: hlnum(:) ! hail number
      double precision,allocatable,save :: hlnumfd(:) ! hail number from frozen drops

      double precision,allocatable,save :: fdvol(:) ! frozen drop volume
      double precision,allocatable,save :: fdvol1(:) ! frozen drop volume
      double precision,allocatable,save :: fdms(:) ! frozen drop mass
      double precision,allocatable,save :: fdv(:) ! frozen drop particle volume

      double precision,allocatable,save :: grv(:) ! graupel particle volume
      double precision,allocatable,save :: hlv(:) ! hail particle volume
      double precision,allocatable,save :: ticms(:) ! total ice crystal mass
      double precision,allocatable,save :: rnms(:)  ! rain mass
      double precision,allocatable,save :: liqms(:)  ! mass of liquid water on ice
      double precision,allocatable,save :: swms(:)  ! snow mass
      double precision,allocatable,save :: qcms(:)  ! cloud droplet mass
      double precision,allocatable,save :: qvms(:)  ! vapor mass
      double precision,allocatable,save :: iwcmx(:) ! max precip ice content
      double precision,allocatable,save :: grpotc(:)  ! gravitational potential energy of cloud
      double precision,allocatable,save :: grpotp(:)  ! gravitational potential energy of precip
      double precision,allocatable,save :: udke(:)  ! updraft kinetic energy (w > 0)
      double precision,allocatable,save :: ddke(:)  ! downdraft kinetic energy (w < 0)
      
      double precision grvoltot, hlvoltot, grmstot, hlmstot, ticmstot,swmstot
      double precision hlnumtot, hlnumfdtot, hlnumfalltot, hlnumfdfalltot, hlfdmasstot
      double precision fdvoltot, fdvoltot1, fdvoltota, fdvoltotb, fdmstot, fdmstota, fdmstotb
      double precision grvoltot1, hlvoltot1
      double precision  grmstota, hlmstota, grmstotb, hlmstotb
      double precision  grvtota, hlvtota
      double precision  grvoltota, hlvoltota, grvoltotb, hlvoltotb
      double precision  rnmstot, rnmstota, rnmstotb
      double precision  liqtot
      double precision  qctot, qvtot, qctota, qctotb
      double precision :: udketot, ddketot, grpotctot, grpotptot

      integer, parameter :: ntot = 100
      double precision  mpitotin(ntot), mpitotout(ntot)
      double precision, allocatable :: mpitotinth(:,:),mpitotoutth(:,:)
      real, allocatable :: mpitotinthr(:,:),mpitotoutthr(:,:)
      
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      real poo, rd, cp, cap
      integer lt273,lt263,lt253,lt243,lt233
      
!      save temc,
      save lt273,lt263,lt253,lt243,lt233
!      save dinv
!      common/strz/z1d2,dlz,nz1d
      
      real eztop, pottop
      common/bndcz/eztop, pottop
      
      real dbt

      real qcmin,qimin,qrmin,qsmin,qhmin
      
      real gt1

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

! #ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.
!      logical :: debug_mpi = .true.

      if (debug_mpi) write(0,*) "DRIVER: ENTERED SUBROUTINE",my_rank


! #endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc Begin Execute
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      allocate (    tt0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate (    tt7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate (   ssat(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      ssat(:,:,:) = 0.0
      allocate ( ssfilt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate (  ssati(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )

#ifndef CM1
      allocate ( dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz) )
      allocate ( pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz) )
#else
      ibsd = 1
      iesd = nx 
      jbsd = 1 
      jesd = ny
#endif
      
      allocate ( xfall(nx,ny,na) )
      nrate2d = nraintypes
      allocate ( rate2d(nx,ny,nrate2d) )
      IF ( .not. allocated( rate2dtot ) ) allocate ( rate2dtot(nrate2d) )
      
!      IF ( .not. allocated( tq1) ) allocate ( tq1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate ( tq1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,ntq1) )
      tq1(:,:,:,:) = 0.0
      
      IF ( .not. allocated( kcldtop ) ) THEN 
        allocate( kcldtop(nx,ny) )
        IF ( iraintypes >= 1 ) THEN
          numproc = numproc + nraintypes 
        ENDIF
        IF ( ldbzchanger >= 1 ) THEN
          numproc = numproc + 1 
        ENDIF
        IF ( ldbzchangeh >= 1 ) THEN
          numproc = numproc + 1 
        ENDIF
        IF ( numproc > NUMPROCTMAX ) THEN
          IF ( my_rank == 0 ) write(0,*) 'Error: increase max value of numproc for format statements'
          call commasmpi_abort()
        ENDIF
        
        allocate( thproc(nzend,numproc) )
        allocate( thproctot(nzend,numproc) )
        allocate( proctot(numproc) )
        allocate( procnames(numproc) )
        thproctot(:,:) = 0.0
      ENDIF

      IF ( .not. allocated( udmf ) ) THEN
         allocate ( udmf(nzend),  &
                    udmf1(nzend),  &
                    udmf10(nzend),  &
                    uicmf(nzend),  &
                    grvol(nzend),  &
                    grvol1(nzend),  &
                    hlvol(nzend),  &
                    hlvol1(nzend),  &
                    hlnum(nzend), hlnumfd(nzend), &
                    grms(nzend),  &
                    hlms(nzend),  &
                    grv(nzend),  &
                    hlv(nzend),  &
                    ticms(nzend),  &
                    rnms(nzend),  &
                    liqms(nzend),  &
                    swms(nzend),  &
                    qcms(nzend),  &
                    qvms(nzend),  &
                    iwcmx(nzend), &
                    grpotc(nzend), &
                    grpotp(nzend), &
                    udke(nzend), &
                    ddke(nzend) &
                    )
          
          allocate( qproczmx(nzend), qproczmn(nzend), qproczav(nzend) )
          allocate( temc(nzend) )
          temc(:) = -1000.
         
!         IF ( lf > 1 ) THEN
         allocate ( fdms(nzend),  &
                    fdvol(nzend),  &
                    fdvol1(nzend),  &
                    fdv(nzend) )
         
!         ENDIF
          
      ENDIF
      
      bcy = Min( 1, bcy1 ) ! this turns off periodic mainly for electrification concerns -- lightning not completely working in periodic
     ! bcy = bcy1  ! use this line for testing periodic lightning; microphysics doesn't otherwise care about periodic BC
      bcx = Min( 1, bcx1 )


!      cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
!      cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
      
#ifdef CHGELEC
      
      IF ( .not. allocated( lgtth ) ) THEN 
        allocate( lgtth(nzend,14) )
        lgtth(:,:) = 0.0
        allocate( uz(nz,4) )
        allocate( cioncp(nz,2) )
      IF ( ipelec .ge. 1 ) THEN

        allocate( sc(lscb:lsceq) )
        allocate( sc2(lscb:lsceq) )

      iestag = elec(iex)%istag
      CALL cld_cpu('ELECTRICITY')

!      CALL ELEC_STUFF_INIT()
      
      DO k = 1,4
        uz(1:nz,k)  = muz(k)%flt1d(1:nz)  ! copy of ion mobilities
      ENDDO
      
      DO k = 1,2
        cioncp(1:nz,k)  = cion(k)%flt1d(1:nz) ! copy of fairweather ion concentrations
      ENDDO
      
      CALL cld_cpu('ELECTRICITY')
      ENDIF

      ENDIF
#endif

      CALL cld_cpu('MICROPHYSICS')


      
#ifdef ICE10
      IF ( ndebug .ge. 2 ) THEN
      write(6,*) 'Check mixing ratio:'
      qcmin = 1.0e-09
      qimin = 1.0e-09
      qrmin = 1.0e-07
      qsmin = 1.0e-07
      qhmin = 1.0e-07
        DO kz=1,nz-1
        ix = nx/2
        jy = ny/2
!        DO jy=1,ny-1
!        DO ix=1,nx-1
          if (                                &
     &     an(ix,jy,kz,lc)  .gt. qcmin   .or. &
     &     an(ix,jy,kz,li)  .gt. qimin   .or. &
     &     an(ix,jy,kz,lip) .gt. qimin   .or. &
     &     an(ix,jy,kz,lir) .gt. qimin   .or. &
     &     an(ix,jy,kz,lr)  .gt. qrmin   .or. &
     &     an(ix,jy,kz,ls)  .gt. qsmin   .or. &
     &     an(ix,jy,kz,lgl) .gt. qsmin   .or. &
     &     an(ix,jy,kz,lgm) .gt. qsmin   .or. &
     &     an(ix,jy,kz,lgh) .gt. qsmin   .or. &
     &     an(ix,jy,kz,lf)  .gt. qsmin   .or. &
     &     an(ix,jy,kz,lhl) .gt. qhmin   .or. &
     &     an(ix,jy,kz,lh)  .gt. qhmin ) then
      
      write(6,*) 'check1 :',ix,jy,kz,    &
     & an(ix,jy,kz,lv),an(ix,jy,kz,lir), &
     & an(ix,jy,kz,ls),an(ix,jy,kz,lgl), &
     & an(ix,jy,kz,lv)  +                &
     &     an(ix,jy,kz,lc)  +            &
     &     an(ix,jy,kz,li)  +            &
     &     an(ix,jy,kz,lip) +            &
     &     an(ix,jy,kz,lir) +            &
     &     an(ix,jy,kz,lr)  +            &
     &     an(ix,jy,kz,ls)  +            &
     &     an(ix,jy,kz,lgl) +            &
     &     an(ix,jy,kz,lgm) +            &
     &     an(ix,jy,kz,lgh) +            &
     &     an(ix,jy,kz,lf)  +            &
     &     an(ix,jy,kz,lhl) +            &
     &     an(ix,jy,kz,lh) ,             &
     &     w(ix,jy,kz)
        ENDIF
!        ENDDO
!        ENDDO
        ENDDO
       ENDIF
#endif

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG start"
!      ifv5d = .false. !mpidebug
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG start"
!      ifv5d = .false.  !mpidebug
#endif

#ifndef CM1
      IF ( ifv5d ) THEN  !mpidebug: note these vals for vis5d
        allocate ( tem1(Max(2,ny-1),nx-1,nz-1) )
        allocate ( tem2(nx-1,ny-1,nz-1) )
        IF ( vis5dstridex*vis5dstridez .ne. 1 ) THEN
          allocate ( tem3(Max(2,(ny-1)/vis5dstridex),(nx-1)/vis5dstridex,(nz-1)/vis5dstridez) )
        ENDIF
      ENDIF
#endif
      
      dtp = ntmul*dtp1

#ifdef SWM
      
!      nstep = nstep + 1
      nstep = int(time_real/dt)

!
! dxx = dx(x), etc.
!

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO ix=ixb,ixe
#else
      ixb = 1
      ixe = nx-1
      DO ix=1,nx-1
#endif
        dxx(ix) = 1.0/gx(ix)
      ENDDO

#ifdef MPI
      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      DO jy=jyb,jye
#else
      jyb = 1
      jye = ny-1
      DO jy=1,ny-1
#endif
        dyy(jy) = 1.0/gy(jy)
      ENDDO

#ifdef MPI
      kzb = 1
      kze = ktile
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO kz=kzb,kze
#else
      DO kz=1,nz-1
#endif
        dzz(kz) = 1.0/gz(kz)
      ENDDO

! Set up 1D vertical arrays

!      write(6,*) 'density and pressure base arrays'

#ifdef MPI
      kzb = -1
      kze = ktile+2
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO kz=kzb,kze
#else
      DO kz = 1,nz-1
#endif
       db(kz)= 1.0e5*pinit(kz)**2.509/(287.04*ab(kz,lt))
       pb(kz)= 0.0 ! 1.0e5*pinit(kz)**3.509
       tmp = 1.0e5*pinit(kz)**3.509
!        write(6,*) kz,db(kz),pb(kz),ab(kz,lt)

#ifdef MPI
        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO ix=ixb,ixe
         DO jy=jyb,jye
#else
        DO ix=1,nx-1
          DO jy=1,ny-1
#endif
!          IF ( Abs(p2(ix,jy,kz)) .gt. 0.2 ) THEN
!           print*,'ice driver: pinit,p2,an(lt) = ',pinit(kz),p2(ix,jy,kz),an(ix,jy,kz,lt)
!     &        ,u(ix,jy,kz),w(ix,jy,kz)
!          ENDIF
            dn(ix,jy,kz)= 1.0e5*(pinit(kz)+p2(ix,jy,kz))**2.509/(287.04*an(ix,jy,kz,lt) )
            pn(ix,jy,kz)= 1.0e5*(pinit(kz)+p2(ix,jy,kz))**3.509 ! - pb(kz)
            IF ( lppert > 1 ) THEN
             axtra(ix,jy,kz,lppert) = pn(ix,jy,kz) - tmp
            ENDIF
          ENDDO
        ENDDO
!        write(6,*) 'pn(nx/2,ny/2,kz),pb = ',p2(nx/2,ny/2,kz),pb(kz),
!     &     db(kz),dn(nx/2,ny/2,kz),an(nx/2,ny/2,kz,lt)
      ENDDO
#endif

      ldumpsp = .false.
      
      poo = 1.0e+05
      rd = 287.04
      cp = 1004.0
      cap = rd/cp 

      
!      IF ( neelec .le. 1 .or. na .lt. lscb) THEN
!        ipelec = 0
!      ENDIF

#ifdef CHGELEC
      IF ( ipelec .gt. 0 .and. .not. associated(elec(1)%flt3d) ) THEN
        write(0,*) 'Problem! Must choose microphys option with electrification!!'
        STOP
      ENDIF
#endif

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 0"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 0"
#endif

! ####################################################
      IF (initcond .eq. 0) THEN ! {
         
         IF ( xvdmx .gt. 0.0 ) THEN
         xvrmx = 0.523599*(xvdmx)**3
         ELSE
         xvrmx = xvrmx0
         ENDIF
         
         IF ( dhmn < 0.0 ) THEN
           IF ( my_rank == 0 ) THEN
             write(0,*) 'cannot have negative dhmn! Set to a positive value or use default.', dhmn
           ENDIF
           call commasmpi_abort()
         ENDIF
         
          xvhmn = 0.523599*(dhmn)**3

         IF ( dhmx <= 0.0 ) THEN
           xvhmx = xvhmx0
         ELSE
           xvhmx = 0.523599*(dhmx)**3
         ENDIF
         
         IF ( irainbreak == -1 ) THEN
           IF ( ipconc >= 6 ) THEN
             irainbreak = 2
           ELSE
             irainbreak = 0
           ENDIF
         ENDIF

         IF ( ihlcnh <= 0 ) THEN
           IF ( ipconc == 5 ) THEN
            ihlcnh = 3
           ELSEIF ( ipconc >= 6 ) THEN
            ihlcnh = 3
           ENDIF
         ENDIF

         IF ( ipconc == 5 .and. imorrgdnglimit >= 2 ) THEN
         ! convert morrdnglimit to xvhmx equivalent
           cwch = ((3. + alphah)*(2. + alphah)*(1.0 + alphah))**(-1./3.)
           xvhmx = pi/6.0*(morrdnglimit/cwch)**3
           dhmx = morrdnglimit/cwch
         ENDIF

         IF ( dhlmx <= 0.0 ) THEN
           xvhlmx = xvhlmx0
         ELSE
           xvhlmx = 0.523599*(dhlmx)**3
         ENDIF
         xvhlmn = 0.523599*(dhlmn)**3

         IF ( dsmx <= 0.0 ) THEN
           xvsmx = xvsmx0
         ELSE
           xvsmx = 0.523599*(dsmx)**3
         ENDIF
         
!         varnum = numvars
!         nzdbz = 0
               
! ###
!        initcond = 1
        
#ifdef CHGELEC
      CALL ELEC_STUFF_INIT(nzend)
      IF ( ipelec .ge. 1 ) THEN ! {

      CALL cld_cpu('MICROPHYSICS')
      CALL cld_cpu('ELECTRICITY')

        call ELEC_DRIVE1(gd,                                   &
     &                     pinit, ab,                          &        ! base state
     &                         p2,                             &       ! pi'
     &                         an,                             &        ! scalars
     &                         xfalltot,                       &        ! precip
     &                        dtp1, gx,gy,gz,z1d,              &        ! time step, vertical grid spacing
     &                        nx, ny, nz, nor, na,             &      
     &                        u,v,w,                           &        ! vn 
     &                        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,   &       ! temporary arrays
     &                        dx,dy,dz,                        &       ! nstep, start,tstat,
     &                        ntmul,iunit,istag,jstag,kstag,   &      
     &                        dbz, io_flag, time, time_real, dt, &
     &                        elec, cion, muz, gxt, gyt, gzt,        & ! electricity arrays
     &                        tstat, lstt, bcx, bcy, tstop,        &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi )
        CALL cld_cpu('ELECTRICITY')
        CALL cld_cpu('MICROPHYSICS')

        ENDIF ! } ipelec .ge. 1
#endif


!  calculate temperature profile
         write(iunit,*) 'Temperature profile: kz, temp (K), kzbeg = ',kzbeg
        call flush(iunit)

      END IF ! (initcond .eq. 0) }

#ifdef MPI

        kzb = ktile+2
        kze = -2
        IF ( kzbeg == nzbeg ) kze = 1
        if (kzend .eq. nzend) kzb = kzend-kzbeg

        DO kz=kzb,kze,-1
#else
        DO kz=nz-1,1,-1
#endif
          temc(kzbeg-1+kz) = pinit(kz)*ab(kz,lt)

!         write(iunit,'(1x,i3,2x,f8.3)') kz,temc(kz)
!         write(iunit,*) pn(1,1,kz),pb(1,1,kz),cion(kz,1),cion(kz,2),
!     &    cion(kz,1)-cion(kz,2)
! end of temporary cut for SWM
        ENDDO

#ifdef MPI
      IF ( nprock > 1 ) THEN
 ! communicate temc array using max from each level
      allocate( mpitotinth((nzend),1))
      allocate(mpitotoutth((nzend),1))
      
      mpitotinth(:,:) = 0.0d0
      
      DO k = kzbeg,kzend
       mpitotinth(k, 1) = temc(k)
      ENDDO

      n = 1*(nzend)

      CALL MPI_AllReduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_Max, my_comm, mpi_error_code)
      
       DO k = 1,nzend
        temc(k)   =  mpitotoutth(k, 1)
       ENDDO
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )
      
      ENDIF
#endif
!
!  Find levels closest to 0, -10, and -20 Celsius for reporting of UDMF
!
      lt273 = nzend - 1
      lt263 = nzend - 1
      lt253 = nzend - 1
      lt243 = nzend - 1
      lt233 = nzend - 1

#ifdef MPI
      kzb = nzend-1
      kze = 1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg

      do kz=kzb,kze,-1
#else
      DO kz=nz-1,1,-1
#endif
        IF ( temc(kz)  .lt. 273.0 ) THEN
          lt273 = Max(1,kz-1)
        ENDIF
        IF ( temc(kz) .lt. 263.0 ) THEN
          lt263 = Max(1,kz-1)
        ENDIF
        IF ( temc(kz)  .lt. 253.0 ) THEN
          lt253 = Max(1,kz-1)
        ENDIF
        IF ( temc(kz)  .lt. 243.0 ) THEN
          lt243 = Max(1,kz-1)
        ENDIF
        IF ( temc(kz)  .lt. 233.0 ) THEN
          lt233 = Max(1,kz-1)
        ENDIF
      ENDDO
      
        IF ( Abs(temc(lt273) - 273.15 ) .gt.    &
     &       Abs(temc(lt273+1) - 273.15 ) ) THEN
          lt273 = lt273 + 1
        ENDIF
        IF ( Abs(temc(lt263) - 263.15 ) .gt.    &
     &       Abs(temc(lt263+1) - 263.15 ) ) THEN
          lt263 = lt263 + 1
        ENDIF
        IF ( Abs(temc(lt253) - 253.15 ) .gt.    &
     &       Abs(temc(lt253+1) - 253.15 ) ) THEN
          lt253 = lt253 + 1
        ENDIF
        IF ( Abs(temc(lt243) - 243.15 ) .gt.    &
     &       Abs(temc(lt243+1) - 243.15 ) ) THEN
          lt243 = lt243 + 1
        ENDIF
        IF ( Abs(temc(lt233) - 233.15 ) .gt.    &
     &       Abs(temc(lt233+1) - 233.15 ) ) THEN
          lt233 = lt233 + 1
        ENDIF

        IF ( initcond .eq. 0 ) THEN !{
        initcond = 1
        
        IF ( my_rank == 0 ) THEN
        write(iunit,*) 'lt273 = ',lt273
        write(iunit,*) 'lt263 = ',lt263
        write(iunit,*) 'lt253 = ',lt253
        write(iunit,*) 'lt243 = ',lt243
        write(iunit,*) 'lt233 = ',lt233
        ENDIF

        END IF ! (initcond .eq. 0) }

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 2"
#else
      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 2"
#endif


! ###


#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 3"
#else
      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 3"
#endif

! ####################################################
      
       lnstat = .false.
       

      if (idelay .ge. 10) THEN
        iwrite = 1
      else 
        iwrite = 0
        idelay = idelay + 1
      end if

      if (idelaycg .ge. 10) THEN
        iwritecg = 1
      else 
        iwritecg = 0
        idelaycg = idelaycg + 1
      end if

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          t1(ix,jy,kz) = 0.5*Sqrt(                    &
     &       ((u(ix,jy,kz)+u(ix+istag,jy,kz)))**2 +   &
     &       ((v(ix,jy,kz)+v(ix,jy+jstag,kz)))**2 +   &
     &       ((w(ix,jy,kz)+w(ix,jy,kz+kstag)))**2 )
        END DO
        END DO
        END DO

        wmin = 0.0
        wmax = 0.0
        vmax = 0.0
        
        i  = 1; j  = 1; k  = 1
        i1 = 1; j1 = 1; k1 = 1
        i2 = 1; j2 = 1; k2 = 1
        

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          IF ( wmax .lt. w(ix,jy,kz) ) THEN
            wmax = w(ix,jy,kz) 
            i = ix
            j = jy
            k = kz
          END IF
          IF ( wmin .gt. w(ix,jy,kz) ) THEN
            wmin = w(ix,jy,kz)
            i2 = ix
            j2 = jy
            k2 = kz
          END IF
          IF ( vmax .lt. t1(ix,jy,kz) ) THEN
            vmax = t1(ix,jy,kz) 
            i1 = ix
            j1 = jy
            k1 = kz
          END IF
        END DO
        END DO
        END DO

#ifdef MPI
! find global integrated rate max
       mpitotin(1)  =  wmax
       mpitotin(2)  =  vmax
       mpitotin(3)  = -wmin
       n = 3

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       wmax = mpitotout(1)  
       vmax = mpitotout(2)  
       wmin = -mpitotout(3)
      ENDIF

! find global integrated rate max
!       mpitotin(1)  = wmin
!       n = 1

!      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, 0, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
!       wmin = mpitotout(1)  
!      ENDIF
        IF ( my_rank == 0 ) THEN
        write(iunit,'(a,f12.5,3i6)') 'wmax = ',wmax ! ,i,j,k
        write(iunit,'(a,f12.5,3i6)') 'vmax = ',vmax ! ,i1,j1,k1
        write(iunit,'(a,f12.5,3i6)') 'wmin = ',wmin ! ,i2,j2,k2
        ENDIF
#else
        write(iunit,'(a,f12.5,3i6)') 'wmax = ',wmax,i,j,k
        write(iunit,'(a,f12.5,3i6)') 'vmax = ',vmax,i1,j1,k1
        write(iunit,'(a,f12.5,3i6)') 'wmin = ',wmin,i2,j2,k2

#endif
        

!        DO ia = 1,na
!        DO kz=1,nz-kd1
!        DO jy=1,ny-jd1
!        DO ix=1,nx-id1
!          IF ( .not. ( an(ix,jy,kz,ia)  .lt. 1.0e10  .and.   &
!     &              an(ix,jy,kz,ia)  .gt. -1.0e10  )   ) THEN
!           write(iunit,'(a,4i5,1x,1pe13.5)')   &
!     &      'NaN at ',ix,jy,kz,ia,an(ix,jy,kz,ia)
!          END IF
!        END DO
!        END DO
!        END DO
!        END DO
      
!      timetd1 = etime(tarray)
!      timetd1 = tarray(1)

!  
!  Count up total net positive (chgpos) and net negative charge (chgneg)
!


#ifdef CHGELEC

      IF ( ipelec .ge. 1 .and. lsce .gt. 1) THEN

      CALL cld_cpu('MICROPHYSICS')
      CALL cld_cpu('ELECTRICITY')
      
      CALL ELEC_DRIVE2     (gd,                                      &
     &                     pinit, ab,                                &  ! base state
     &                         p2,                                   & ! pi'
     &                         an,                                   &  ! scalars
     &                         xfalltot,                             &  ! precip
     &                        dtp1, gx,gy,gz,z1d,                    &  ! time step, vertical grid spacing
     &                        nx, ny, nz, nor, na,                   &
     &                        u,v,w,                                 &  ! vn 
     &                        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,         & ! temporary arrays
     &                        dx,dy,dz,                              & ! nstep, start,tstat,
     &                        ntmul,iunit,istag,jstag,kstag,         &
     &                        dbz, io_flag, time, time_real, dt, &
     &                        elec, cion, muz, gxt, gyt, gzt,        & ! electricity arrays
     &                        tstat, lstt, bcx, bcy, tstop,        &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,chgadvtemp )

      CALL cld_cpu('ELECTRICITY')
 
      CALL cld_cpu('MICROPHYSICS')

      ENDIF ! ipelec

#endif

        xfall(:,:,:) = 0.0
        rate2d(:,:,:) = 0.0
        rate2dtot(:) = 0.0

! Microphysics calls

       iflag = lv5dwrite .or. lncwrite .or. lrst .or. lprt .or. lstt

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: ABOUT TO MAKE MICROPHYSICS CALL"
#else
      if (ndebug .gt. 0) write(iunit,*) "DRIVER: ABOUT TO MAKE MICROPHYSICS CALL "
#endif

      IF ( lpredict ) THEN
#if 0
#ifndef MPI
      call ice10q01                       &
     &  (ntmul,nx,ny,nz,na,nba,nv,nstep   &
     &  ,nor,istag,jstag,kstag,itopo   &
     &  ,iwrite                          &
     &  ,cwmasn,cwmasx                       &
     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d   &
     &  ,xfall                                        &
     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9,tq1   &
     &  ,ab,an,db,dn,p2,pinit   &
     &  ,pb,pn,u,v,w,iunit   &
     &  ,ssat,ssfilt,t0,t7,ssati   &
     &  ,elec,bcx,bcy )
#endif
#endif
      thproc(:,:) = 0.0
#ifdef ICE3
#ifdef LFO3
!      write(0,*) 'call ice_zieg, nstep = ',nstep
!      call lfoice   &
!     &  (ntmul,nx,ny,nz,na,nba,nv,nstep,time_real   &
!     &  ,nor,norz,istag,jstag,kstag,itopo   &
!     &  ,iwrite   &
!     &  ,cwmasn,cwmasx   &
!     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d   &
!     &  ,xfall   &
!     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9,tq1,ntq1   &
!     &  ,ab,an,db,dn,p2,pinit   &
!     &  ,pb,pn,u,v,w,iunit   &
!     &  ,ssat,ssfilt,t0,t7,ssati,thproc,numproc )

#else
!      write(0,*) 'call ice_zieg, nstep = ',nstep
      io_flag_tmp = .true.
      call ice_zieg   &
     &  (ntmul,nx,ny,nz,na,nba,nv,nstep,time_real   &
     &  ,nor,norz,istag,jstag,kstag,itopo   &
     &  ,iwrite,io_flag_tmp   &
     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d   &
     &  ,xfall,rate2d,nrate2d, xfalltot   &
     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9,tq1,ntq1   &
     &  ,ab,an,db,dn,p2,pinit   &
     &  ,pb,pn,u,v,w,km,iunit   &
     &  ,ssat,ssfilt,t0,t7,ssati,thproc,numproc   &
#ifndef CM1
     &  ,elec, axtra  &
#endif
     &  ,bcx,bcy)
!      write(0,*) 'done ice_zieg, nstep = ',nstep

#endif
     
!      write(0,*) 'done ice_zieg'

!      call ice3elec
!     &  (ntmul,nx,ny,nz,na,nba,nv,nstep
!     &  ,nor,istag,jstag,kstag,itopo
!     &  ,iwrite,ipelec,nliter,ipotslv,ilghtest,ilight,nhelec,iscreen
!     &  ,cwccn,cwmasn,cwmasx,ipconc,icorona,cimn,cimx
!     &  ,kdbcdn
!     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d
!     &  ,xfall
!     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9
!     &  ,ab,an,db,dn,p2,pinit
!     &  ,pb,pn,u,v,w,iunit
!     &  ,ssat,ssfilt,t0,t7)
!cc     &  ,elec,neelec)
#endif

#ifdef HCMON
! call hcm
      CALL cld_cpu('MICROPHYSICS')
       call hcm_gs    &
     &  (nstep,nx,ny,nz,na    &
     &  ,nor,norz,istag,jstag,kstag    &
     &  ,dtp1,dz    &
     &  ,ab,an,db,dn, w    &
     &  ,pb,pn,p2,pinit    &
     &  ,dbz,iflag)
      CALL cld_cpu('MICROPHYSICS')
#endif

#ifdef TAKON
!      CALL cld_cpu('MICROPHYSICS')
! call TAKAHASHI
       call takinterface    &
     &  (nstep,nx,ny,nz,na    &
     &  ,nor,norz,istag,jstag,kstag    &
     &  ,dtp1,dx,dy,dz,dxx,dyy,dzz,z1d    &
     &  ,ab,an,db,dn, w , st   &
     &  ,st(-nor+1,-nor+1,-nor+1,lnr)  &
     &  ,st(-nor+1,-nor+1,-nor+1,lni)  &
     &  ,st(-nor+1,-nor+1,-nor+1,lnh)  &
     &  ,an(-nor+1,-nor+1,-nor+1,lv)  &
     &  ,an(-nor+1,-nor+1,-nor+1,lccn)  &
     &  ,an(-nor+1,-nor+1,-nor+1,lcin)  &
     &  ,an(-nor+1,-nor+1,-nor+1,lcng)  &
     &  ,pb,pn,p2,pinit                &
     &  ,xfall                         &
     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9  &
     &  ,elec, axtra, vzf  &
     &  ,dbz,iflag,iunit,thproc,numproc)

!      CALL cld_cpu('MICROPHYSICS')
#endif

#ifdef WARM
      call warmzieg   &
     &  (ntmul,nx,ny,nz,na,nba,nv,nstep   &
     &  ,nor,istag,jstag,kstag,itopo   &
     &  ,iwrite   &
     &  ,cwmasn,cwmasx   &
     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d   &
     &  ,xfall   &
     &  ,tt0,t1,t2,t3,t4,t5,t6,tt7,t8,t9   &
     &  ,ab,an,db,dn,p2,pinit   &
     &  ,pb,pn,u,v,w,iunit   &
     &  ,ssat,ssfilt,t0,t7,ssati,thproc,numproc,axtra)
#endif

      DO i = 1,numproc
       DO kz=kzbeg,kzend
         thproctot(kz,i) = thproctot(kz,i) + thproc(kz,i)
       ENDDO
      ENDDO
!      write(0,*) 'thproctot n-1,n: ', my_rank, thproctot(1,numproc-1),thproctot(1,numproc)

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DONE WITH MICROPHYSICS CALL"
#else
      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DONE WITH MICROPHYSICS CALL"
#endif

!
! Update the fallout totals
!
      IF ( ifv5d ) THEN
        tem2(1:nx-1,1:ny-1,1:nz-1) = t9(1:nx-1,1:ny-1,1:nz-1)
      ENDIF
      
      dvt = 1./dtp1

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      DO ix = ixb,ixe
        DO jy = jyb,jye
#else
      DO ix = 1,nx-1
        DO jy = 1,ny-1
#endif
          dv = dn(ix,jy,1)*dzz(1) ! *dxx(ix)*dyy(jy)
          
          ! rain rate = xfall*airdens*(dz/dt)*(1/raindens)*(1000. mm/m) = mm per m**2 (or kg/m**2 since raindens = 1000.)
          tmp = xfall(ix,jy,lr)
          tmppwfs = 0.0
          tmppwfh = 0.0
          tmppwff = 0.0
          tmppwfhl = 0.0
          IF ( lsw .gt. 1  ) tmppwfs = xfall(ix,jy,lsw)
          IF ( lhw .gt. 1  ) tmppwfh = xfall(ix,jy,lhw)
          IF ( lfw .gt. 1  ) tmppwff = xfall(ix,jy,lfw)
          IF ( lhlw .gt. 1 ) tmppwfhl = xfall(ix,jy,lhlw)
          tmppwf = tmppwfh + tmppwff + tmppwfhl
          tmp = tmp + tmppwf + tmppwfs ! all liquid (rain + liquid fractions)
          
          IF ( myprock == 1 ) THEN
          xfalltot(ix,jy,1) = tmp*dv*dvt ! total precip rate
          xfalltot(ix,jy,3) = xfalltot(ix,jy,3) + tmp*dv ! tot precip accum

            IF ( iprainacc2 > 0 ) THEN
              xfalltot(ix,jy,iprainacc2) = xfalltot(ix,jy,iprainacc2) + dv*tmp
            ENDIF

            IF ( ipwfrat > 0 .and. ipwfacc > 0 ) THEN
             xfalltot(ix,jy,ipwfrat) = (tmppwfs+tmppwf)*dv*dvt
             xfalltot(ix,jy,ipwfacc) = xfalltot(ix,jy,ipwfacc) + (tmppwfs+tmppwf)*dv
            ENDIF
          ENDIF
          
          tmp = 0.0
          IF ( lrauto .gt. 1 .and. iprainaccauto > 1 ) THEN
             xfalltot(ix,jy,iprainaccauto) = xfalltot(ix,jy,iprainaccauto) + dv*xfall(ix,jy,lrauto)
             IF ( iprainaccauto2 > 0 ) THEN
             xfalltot(ix,jy,iprainaccauto2) = xfalltot(ix,jy,iprainaccauto2) + dv*xfall(ix,jy,lrauto)
             ENDIF
          ENDIF
          IF ( lrmelt .gt. 1 .and. iprainaccmelt > 1 ) THEN
             xfalltot(ix,jy,iprainaccmelt) = xfalltot(ix,jy,iprainaccmelt) +  &
                       dv*(tmppwf +tmppwfs + xfall(ix,jy,lrmelt))
             IF ( iprainaccmelt2 > 0 ) THEN
              xfalltot(ix,jy,iprainaccmelt2) = xfalltot(ix,jy,iprainaccmelt2) +  &
                       dv*(tmppwf +tmppwfs + xfall(ix,jy,lrmelt))
             ENDIF
          ENDIF
          IF ( lrshed .gt. 1 .and. iprainaccshed > 1 ) THEN
             xfalltot(ix,jy,iprainaccshed) = xfalltot(ix,jy,iprainaccshed) + dv*xfall(ix,jy,lrshed)
             IF ( iprainaccshed2 > 0 ) THEN
               xfalltot(ix,jy,iprainaccshed2) = xfalltot(ix,jy,iprainaccshed2) + dv*xfall(ix,jy,lrshed)
             ENDIF
          ENDIF


          IF ( ipfdauto > 1 ) THEN
             xfalltot(ix,jy,ipfdauto) = xfalltot(ix,jy,ipfdauto) + rate2d(ix,jy,1)
             IF ( ipfdauto2 > 0 ) THEN
               xfalltot(ix,jy,ipfdauto2) = xfalltot(ix,jy,ipfdauto2) + rate2d(ix,jy,1)
             ENDIF
             rate2dtot(1) = rate2dtot(1) + rate2d(ix,jy,1)
          ENDIF
          IF ( ipfdshed > 1 ) THEN
             xfalltot(ix,jy,ipfdshed) = xfalltot(ix,jy,ipfdshed) + rate2d(ix,jy,2)
             IF ( ipfdshed2 > 0 ) THEN
               xfalltot(ix,jy,ipfdshed2) = xfalltot(ix,jy,ipfdshed2) + rate2d(ix,jy,2)
             ENDIF
             rate2dtot(2) = rate2dtot(2) + rate2d(ix,jy,2)
          ENDIF
          IF ( ipfdmelt > 1 ) THEN
             xfalltot(ix,jy,ipfdmelt) = xfalltot(ix,jy,ipfdmelt) + rate2d(ix,jy,3)
             IF ( ipfdmelt2 > 0 ) THEN
               xfalltot(ix,jy,ipfdmelt2) = xfalltot(ix,jy,ipfdmelt2) + rate2d(ix,jy,3)
             ENDIF
             rate2dtot(3) = rate2dtot(3) + rate2d(ix,jy,3)
          ENDIF


          ! hail rate = xfall*airdens*(dz/dt) = kg per m**2 per second
          tmpdp = 0
          IF ( lg > lr ) THEN
          DO il = lg,lhab
           tmpdp = tmpdp + xfall(ix,jy,il)
          ENDDO

           tmpdp = tmpdp - tmppwf
!           IF ( lhw  .gt. 1 ) tmpdp = Max( 0.0, tmpdp - xfall(ix,jy,lhw))
!           IF ( lfw  .gt. 1 ) tmpdp = Max( 0.0, tmpdp - xfall(ix,jy,lfw))
!           IF ( lhlw .gt. 1 ) tmpdp = Max( 0.0, tmpdp - xfall(ix,jy,lhlw) )

          ENDIF
          
          IF ( myprock == 1 ) THEN
          xfalltot(ix,jy,2) = tmpdp*dv*dvt ! total ice rate
          xfalltot(ix,jy,4) = xfalltot(ix,jy,4) + tmpdp*dv
          ENDIF
          
          IF ( lh > 1 .and. ipgracc > 0 ) THEN
             tmpdp = 0
!             IF ( lhw > 1 ) THEN
!              ! subtract liquid 
!               tmpdp = xfall(ix,jy,lhw)
!             ENDIF
             ! graupel accumulation
             xfalltot(ix,jy,ipgracc) = xfalltot(ix,jy,ipgracc) + dv*Max( 0.0, xfall(ix,jy,lh) - tmppwfh )

             IF ( ipgracc2 > 0 ) THEN
               xfalltot(ix,jy,ipgracc2) = xfalltot(ix,jy,ipgracc2) + dv*Max( 0.0, xfall(ix,jy,lh) - tmppwfh )
             ENDIF

             IF ( lnh > 1 ) THEN
             ! graupel number accumulation (total)
               xfalltot(ix,jy,ipgrnumacc) = xfalltot(ix,jy,ipgrnumacc) + dv*xfall(ix,jy,lnh)
               IF ( ipgrnumacc2 > 0 ) THEN
                 xfalltot(ix,jy,ipgrnumacc2) = xfalltot(ix,jy,ipgrnumacc2) + dv*xfall(ix,jy,lnh)
               ENDIF
             ENDIF
             IF ( ipgrdiam > 0 .and. io_flag) THEN
                IF ( xfalltot(ix,jy,ipgracc) - precip_old(ix,jy,ipgracc) > 1.e-8 .and. &
                     (xfalltot(ix,jy,ipgrnumacc) - precip_old(ix,jy,ipgrnumacc))*(60./this)*50.**2 > 10. ) THEN ! at least 10 stones per 50^2 m^2 per minute
                   
                   tmpdp = (xfalltot(ix,jy,ipgracc) - precip_old(ix,jy,ipgracc))/  &
                         (xfalltot(ix,jy,ipgrnumacc) - precip_old(ix,jy,ipgrnumacc)) ! mean gr mass in grid square
                   xfalltot(ix,jy,ipgrdiam) = 1000.*(tmpdp*6.0/(3.14159*900.0) )**(1.0/3.0)
                ELSE
                   xfalltot(ix,jy,ipgrdiam) = 0.0
                ENDIF
             ENDIF
          
          ENDIF

          IF ( lf > 1 .and. ipfdacc > 0 ) THEN
             tmpdp = 0
!             IF ( lfw > 1 ) THEN
!              ! subtract liquid 
!               tmpdp = xfall(ix,jy,lfw)
!             ENDIF
             ! graupel accumulation
             xfalltot(ix,jy,ipfdacc) = xfalltot(ix,jy,ipfdacc) + dv*Max( 0.0, xfall(ix,jy,lf) - tmppwff )
             IF ( ipfdacc2 > 0 ) THEN
               xfalltot(ix,jy,ipfdacc2) = xfalltot(ix,jy,ipfdacc2) + dv*Max( 0.0, xfall(ix,jy,lf) - tmppwff )
             ENDIF
             IF ( lnf > 1 ) THEN
             ! graupel number accumulation (total)
               xfalltot(ix,jy,ipfdnumacc) = xfalltot(ix,jy,ipfdnumacc) + dv*xfall(ix,jy,lnf)
               IF ( ipfdnumacc2 > 0 ) THEN
                 xfalltot(ix,jy,ipfdnumacc2) = xfalltot(ix,jy,ipfdnumacc2) + dv*xfall(ix,jy,lnf)
               ENDIF
             ENDIF
             IF ( ipfddiam > 0 .and. io_flag) THEN
                IF ( xfalltot(ix,jy,ipfdacc) - precip_old(ix,jy,ipfdacc) > 1.e-8 .and. &
                     (xfalltot(ix,jy,ipfdnumacc) - precip_old(ix,jy,ipfdnumacc))*(60./this)*50.**2 > 10. ) THEN ! at least 10 stones per 50^2 m^2 per minute
                   
                   tmpdp = (xfalltot(ix,jy,ipfdacc) - precip_old(ix,jy,ipfdacc))/  &
                         (xfalltot(ix,jy,ipfdnumacc) - precip_old(ix,jy,ipfdnumacc)) ! mean gr mass in grid square
                   xfalltot(ix,jy,ipfddiam) = 1000.*(tmpdp*6.0/(3.14159*900.0) )**(1.0/3.0)
                ELSE
                   xfalltot(ix,jy,ipfddiam) = 0.0
                ENDIF
             ENDIF
          
          ENDIF

          IF ( lhl > 1 .and. iphailacc > 0 ) THEN
!             tmpdp = 0
!             IF ( lhlw > 1 ) THEN
!              ! subtract liquid 
!               tmpdp = xfall(ix,jy,lhlw)
!             ENDIF
             ! hail accumulation
             xfalltot(ix,jy,iphailacc) = xfalltot(ix,jy,iphailacc) + dv*Max( 0.0, xfall(ix,jy,lhl) - tmppwfhl )

             IF ( iphailacc2 > 0 ) THEN
              xfalltot(ix,jy,iphailacc2) = xfalltot(ix,jy,iphailacc2) + dv*Max( 0.0, xfall(ix,jy,lhl) - tmppwfhl )
             ENDIF

             IF ( iphailfacc > 0 .and. lnhlf > 1 .and. lnhl > 1 ) THEN
             IF ( xfall(ix,jy,lnhl) > 0.0 ) THEN
             ! hail mass accumulation from frozen drops
             xfalltot(ix,jy,iphailfacc) = xfalltot(ix,jy,iphailfacc) +  &
                              dv*Max( 0.0, xfall(ix,jy,lhl) - tmppwfhl )*xfall(ix,jy,lnhlf)/xfall(ix,jy,lnhl)
               IF ( iphailfacc2 > 0 ) THEN
                 xfalltot(ix,jy,iphailfacc2) = xfalltot(ix,jy,iphailfacc2) +  &
                              dv*Max( 0.0, xfall(ix,jy,lhl) - tmppwfhl )*xfall(ix,jy,lnhlf)/xfall(ix,jy,lnhl)
               ENDIF
             ENDIF
             ENDIF
             IF ( lnhl > 1 ) THEN
             ! hail number accumulation (total)
               xfalltot(ix,jy,iphailnumacc) = xfalltot(ix,jy,iphailnumacc) + dv*xfall(ix,jy,lnhl)
               IF ( iphailnumacc2 > 0 ) THEN
                xfalltot(ix,jy,iphailnumacc2) = xfalltot(ix,jy,iphailnumacc2) + dv*xfall(ix,jy,lnhl)
               ENDIF
               IF ( lnhlf > 1 .and. iphailfnumacc > 1 ) THEN
                ! hail number accumulation from frozen drops
                 xfalltot(ix,jy,iphailfnumacc) = xfalltot(ix,jy,iphailfnumacc) +  dv*xfall(ix,jy,lnhlf)
                 IF ( iphailfnumacc2 > 0 ) THEN
                 xfalltot(ix,jy,iphailfnumacc2) = xfalltot(ix,jy,iphailfnumacc2) +  dv*xfall(ix,jy,lnhlf)
                 ENDIF
               ENDIF
             ENDIF
             IF ( iphaildiam > 0 .and. io_flag) THEN
                IF ( xfalltot(ix,jy,iphailacc) - precip_old(ix,jy,iphailacc) > 1.e-8 .and. &
                     (xfalltot(ix,jy,iphailnumacc) - precip_old(ix,jy,iphailnumacc))*(60./this)*50.**2 > 10. ) THEN ! at least 10 stones per 50^2 m^2 per minute [1 per (15.8m)^2]
                   
                   tmpdp = (xfalltot(ix,jy,iphailacc) - precip_old(ix,jy,iphailacc))/  &
                         (xfalltot(ix,jy,iphailnumacc) - precip_old(ix,jy,iphailnumacc)) ! mean hail mass in grid square
                   xfalltot(ix,jy,iphaildiam) = 1000.*(tmpdp*6.0/(3.14159*900.0) )**(1.0/3.0) ! assume ice density of 900; factor of 1000 to convert m to mm
                ELSE
                  xfalltot(ix,jy,iphaildiam) = 0.0
                ENDIF
             ENDIF
          ENDIF
          
        ENDDO
      ENDDO
      
      
      
!      IF ( na .ge. lccn ) THEN
!      write(iunit,*) 'check CCN (dr2)'
!      i1 = 0
!        DO kz=1,nz-1
!        DO jy=1,ny-1
!        DO ix=1,nx-1
!         IF ( i1 .gt. 10 ) exit
!         IF ( .not. (an(ix,jy,kz,lccn) .gt. -100. .and.
!     :               an(ix,jy,kz,lccn) .lt. 1.0e10 ) ) THEN
!          write(iunit,*) 'CCN: ', ix,jy,kz,an(ix,jy,kz,lccn)
!          i1 = i1 + 1
!         ENDIF
!        ENDDO
!        ENDDO
!        ENDDO
!      ENDIF


#ifdef CHGELEC
        IF ( imudstop .eq. 1 .and. ipelec .ge. 1 ) THEN
          write(iunit,*) 'STOP: imudstop .eq. 1 .and. ipelec .ge. 1 '
          IF ( my_rank == 0 ) THEN
            write(3,*) 'STOP: imudstop .eq. 1 .and. ipelec .ge. 1 '
          ENDIF
          write(0,*) 'STOP: imudstop .eq. 1 .and. ipelec .ge. 1 '
          STOP
        ENDIF
#endif /* chgelec */
        
        tmin = 4000.0
#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
!          tt0(ix,jy,kz) = an(ix,jy,kz,lt)
!     :        *((pn(ix,jy,kz)+pb(ix,jy,kz))/1.0e5)**(287.04/1004.)
          IF ( tmin .gt. tt0(ix,jy,kz) ) THEN
            tmin = tt0(ix,jy,kz) 
            i = ix
            j = jy
            k = kz
          END IF
          
        END DO
        END DO
        END DO

#ifdef MPI
! find global integrated rate max
       mpitotin(1)  = tmin
       n = 1

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       tmin = mpitotout(1)  
      ENDIF

        IF ( my_rank == 0 ) THEN
        write(iunit,'(a,f12.5,3i6)') 'tmin = ',tmin ! ,i,j,k
        ENDIF
#else
        write(iunit,'(a,f12.5,3i6)') 'tmin = ',tmin,i,j,k
#endif
        


        IF ( lss .gt. 1 ) THEN
        ssmax = 0.0
        i  = 1; j  = 1; k  = 1
        i2 = 1; j2 = 1; k2 = 1

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          IF ( ssmax  .lt.  an(ix,jy,kz,lss) ) THEN
            ssmax =  an(ix,jy,kz,lss)
            i = ix
            j = jy
            k = kz
          END IF
        END DO
        END DO
        END DO

        ssatmax = 0.0

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          IF ( ssatmax  .lt.  ssat(ix,jy,kz) ) THEN
            ssatmax =  ssat(ix,jy,kz)
            i2 = ix
            j2 = jy
            k2 = kz
          END IF
        END DO
        END DO
        END DO

#ifdef MPI
! find global integrated rate max
       mpitotin(1)  =  ssmax
       mpitotin(2)  =  ssatmax
       n = 2

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       ssmax = mpitotout(1)  
       ssatmax = mpitotout(2)  
      ENDIF

#endif

        IF ( my_rank == 0 ) THEN
          write(iunit,*) 'ssmax = ',ssmax ! ,' at ',i,j,k
          write(iunit,*) 'ssatmax = ',ssatmax ! ,' at ',i2,j2,k2
        ENDIF

        ENDIF

        IF ( lqcond >= 1 ) THEN
        condmax = 0.0
        depmax = 0.0

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          condmax = Max( condmax, axtra(ix,jy,kz,lqcond) )
          IF ( lqdep >= 1 ) THEN
          depmax = Max( depmax, axtra(ix,jy,kz,lqdep) )
          END IF
        END DO
        END DO
        END DO

#ifdef MPI
! find global integrated rate max
       mpitotin(1)  =  condmax
       mpitotin(2)  =  depmax
       n = 2

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       condmax = mpitotout(1)  
       depmax = mpitotout(2)  
      ENDIF

#endif

        IF ( my_rank == 0 ) THEN
          condmaxmax = Max( condmaxmax, condmax )
          depmaxmax = Max( depmaxmax, depmax )
          write(iunit,*) 'condmax = ',condmax, condmaxmax
          write(iunit,*) 'depmax = ',depmax, depmaxmax
        ENDIF

        ENDIF


#ifdef CHGELEC

!        IF ( ifv5d .and. lstt .and.
!     :       ( lchgratn .or. lchgrati) ) THEN

! #ifndef TAKON

        IF ( ipelec .ge. 2 .or. idoniconly) THEN

        CALL cld_cpu('MICROPHYSICS')
        CALL cld_cpu('ELECTRICITY')
        
!        write(0,*) 'maxval t3,t8,t2,t6 =',maxval(t3),maxval(t8),maxval(t2),maxval(t6)

#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          dxdy = dxx(ix)*dyy(jy)
          dv = dxdy*dzz(kz)

          elec(icghis)%flt3d(ix,jy,kz) = t1(ix,jy,kz) + t3(ix,jy,kz) + t8(ix,jy,kz) + t2(ix,jy,kz) + t6(ix,jy,kz)
          elec(icghw )%flt3d(ix,jy,kz) = t5(ix,jy,kz)
          IF ( ichgtndnic > 1 ) elec(ichgtndnic )%flt3d(ix,jy,kz) = elec(ichgtndnic )%flt3d(ix,jy,kz) +    &
     &           dtp1*(elec(icghis)%flt3d(ix,jy,kz) +  t5(ix,jy,kz) )
          IF ( icgaddl > 1 ) elec(icgaddl )%flt3d(ix,jy,kz) = tq1(ix,jy,kz,3)
          IF ( icgnoliq > 1 ) elec(icgnoliq)%flt3d(ix,jy,kz) = tq1(ix,jy,kz,4)
          
          IF ( ipnic2d > 0 .and. innic2d > 0 ) THEN
           xfalltot(ix,jy,ipnic2d) = xfalltot(ix,jy,ipnic2d) + Max( 0.0, elec(icghis)%flt3d(ix,jy,kz) )*dv
           xfalltot(ix,jy,innic2d) = xfalltot(ix,jy,innic2d) + Min( 0.0, elec(icghis)%flt3d(ix,jy,kz) )*dv
          ENDIF
        END DO
        END DO
        END DO
        
        CALL cld_cpu('ELECTRICITY')
        CALL cld_cpu('MICROPHYSICS')
        
        ENDIF !! (ipelec .ge. 2)
!#endif /* takon */
#endif /* chgelec */
        
      ELSE !! lpredict = false
    
#ifdef CHGELEC
        IF ( ipelec .ge. 2 ) THEN

        t1(:,:,:) = 0.0
        t2(:,:,:) = 0.0
        t4(:,:,:) = 0.0
        t5(:,:,:) = 0.0
        t6(:,:,:) = 0.0
        t8(:,:,:) = 0.0

!    NOTE: This is only done for recreating time-height data from history data. (lpredict = false)
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
          t3(ix,jy,kz) = elec(icghis)%flt3d(ix,jy,kz) 
          t5(ix,jy,kz) = elec(icghw )%flt3d(ix,jy,kz)
        END DO
        END DO
        END DO
        
        
        ENDIF
#endif /* chgelec */


      ENDIF ! lpredict

!#ifdef TAKON
!        IF ( ipelec .ge. 2 ) THEN
!        t1(:,:,:) = 0.0
!        t2(:,:,:) = 0.0
!        t4(:,:,:) = 0.0
!        t5(:,:,:) = 0.0
!        t6(:,:,:) = 0.0
!        t8(:,:,:) = 0.0
!
!!    NOTE: This is only done for recreating time-height data from history data. (lpredict = false)
!        DO kz=1,nz-kd1
!        DO jy=1,ny-jd1
!        DO ix=1,nx-id1
!          t3(ix,jy,kz) = elec(icghis)%flt3d(ix,jy,kz) 
!          t5(ix,jy,kz) = elec(icghw )%flt3d(ix,jy,kz)
!        END DO
!        END DO
!        END DO
!
!        ENDIF
!
!#endif /* takon */

#ifdef CHGELEC

      IF ( ipelec .ge. 1 ) THEN ! {

      CALL cld_cpu('MICROPHYSICS')
      CALL cld_cpu('ELECTRICITY')

      CALL ELEC_DRIVE3     (gd,                                   &
     &                     pinit, ab,                              &    ! base state
     &                         p2,                                 &   ! pi'
     &                         an,                                 &    ! scalars
     &                         xfalltot,                           &    ! precip
     &                        dtp1, gx,gy,gz,z1d,                  &    ! time step, vertical grid spacing
     &                        nx, ny, nz, nor, na,                 &  
     &                        u,v,w,                               &    ! vn 
     &                        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,       &   ! temporary arrays
     &                        dx,dy,dz,                            &   ! nstep, start,tstat,
     &                        ntmul,iunit,istag,jstag,kstag,       &  
     &                        dbz, io_flag, time, time_real, dt, &
     &                        elec, cion, muz, gxt, gyt, gzt,        & ! electricity arrays
     &                        tstat, lstt, bcx, bcy, tstop,       & 
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi )

      CALL cld_cpu('ELECTRICITY')
      CALL cld_cpu('MICROPHYSICS')

      ENDIF ! } ipelec

#endif /* chgelec */

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: start time-height"
#endif

!      tq1(1:nx-1,1:ny-1,1:nz-1) = t1(1:nx-1,1:ny-1,1:nz-1)
      

!      timetd2 = etime(tarray)
!      timetd2 = tarray(1)
!      write (iunit,*) 'tdmphy1 call time is', timetd2 - timetd1

!
! At this point, the temporary arrays have various rates in them,
! so now find and print out the max and mins
!

      qprocmax = 0.d0
      qprocmin = 1.d0
      qprocavg = 0.d0

#ifdef CHGELEC
      cswimn = 0.0
      cswimx = 0.0
      cswwmn = 0.0
      cswwmx = 0.0
      cghsmn = 0.0
      cghsmx = 0.0
      cghimn = 0.0
      cghimx = 0.0
      cghwmn = 0.0
      cghwmx = 0.0
      cgsmn = 0.0
      cgsmx = 0.0
      cgimn = 0.0
      cgimx = 0.0
      chsmn = 0.0
      chsmx = 0.0
      chimn = 0.0
      chimx = 0.0

      ctswin = 0.0
      ctswip = 0.0
      ctswwn = 0.0
      ctswwp = 0.0
      ctghsn = 0.0
      ctghsp = 0.0
      ctghin = 0.0
      ctghip = 0.0
      ctghwn = 0.0
      ctghwp = 0.0
      ctgsn = 0.0
      ctgsp = 0.0
      ctgin = 0.0
      ctgip = 0.0
      cthsn = 0.0
      cthsp = 0.0
      cthin = 0.0
      cthip = 0.0
      
      cthsdn = 0.0 ! drake and other mechanisms negative to snow/graupel
      cthsdp = 0.0 
      
#endif /* chgelec */

      udv5n = 0.0
      udv5 = 0.0
      udv10 = 0.0
      udv20 = 0.0
      hlnumfalltot = 0.0
      hlnumfdfalltot = 0.0
      hlfdmasstot = 0.0

        icldtop = 0
        zcldtop = 0.0
        kcldtop(:,:) = 0

!      IF ( nstep/nstat*nstat .eq. nstep .or. nstep .eq. 1 ) THEN
!      cgszmn2(kz) = 0.0
!      chszmn2(kz) = 0.0
!      cgizmn2(kz) = 0.0
!      chizmn2(kz) = 0.0
!      cghwzmn2(kz) = 0.0
!      cgszmx2(kz) = 0.0
!      chszmx2(kz) = 0.0
!      cgizmx2(kz) = 0.0
!      chizmx2(kz) = 0.0
!      cghwzmx2(kz) = 0.0
!      cswizmn2(kz) = 0.0
!      cswizmx2(kz) = 0.0
!      cghszmn2(kz) = 0.0
!      cghszmx2(kz) = 0.0
!      cghizmn2(kz) = 0.0
!      cghizmx2(kz) = 0.0
!      
!      ENDIF

#ifndef MPI
! C$DOACROSS LOCAL(kz,jy,ix,tmp), 
! C$&   REDUCTION(ctswin,ctswip,ctghsn,ctghsp,ctghin,ctghip,
! C$&    ctghwn,ctghwp,ctgsn,ctgsp,ctgin,ctgip,cthsn,cthsp,cthin,cthip)
#ifdef CHGELEC
!$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,tmp,dxdy,dv,dv9,tmpg,tmph,gt1,dbt)  &
!$OMP ,REDUCTION (+ : ctswin,ctswip,ctghsn,ctghsp,ctghin,ctghip,cthsdn,cthsdp,    &
!$omp    ctghwn,ctghwp,ctgsn,ctgsp,ctgin,ctgip,cthsn,cthsp,cthin,cthip, &
!$omp    udv5n, udv5, udv10, udv20), REDUCTION(Max : icldtop,zcldtop)
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,tmp,dxdy,dv,dv9,tmpg,tmph,gt1,dbt)  &
!$OMP  ,REDUCTION (+ : udv5n, udv5, udv10, udv20), REDUCTION(Max : icldtop,zcldtop)
#endif
#endif

      DO kz = 1,nzend
      
      qproczmx(kz) = 0.
      qproczmn(kz) = 1.
      qproczav(kz) = 0.
      iavgcnt = 0
      
#ifdef CHGELEC
      cgszmn(kz) = 0.0
      chszmn(kz) = 0.0
      cgizmn(kz) = 0.0
      chizmn(kz) = 0.0
      cghwzmn(kz) = 0.0
      cgszmx(kz) = 0.0
      chszmx(kz) = 0.0
      cgizmx(kz) = 0.0
      chizmx(kz) = 0.0
      cghwzmx(kz) = 0.0
      cswizmn(kz) = 0.0
      cswizmx(kz) = 0.0
      cswwzmn(kz) = 0.0
      cswwzmx(kz) = 0.0
      cghszmn(kz) = 0.0
      cghszmx(kz) = 0.0
      cghizmn(kz) = 0.0
      cghizmx(kz) = 0.0
#endif

      udmf(kz) = 0.0  ! updraft mass flux
      udmf1(kz) = 0.0  ! updraft mass flux
      udmf10(kz) = 0.0  ! updraft mass flux
      uicmf(kz) = 0.0 ! upward ice crystal mass flux
      grvol(kz) = 0.0 ! graupel 'echo' volume
      grvol1(kz) = 0.0 ! graupel 'echo' volume
      hlvol(kz) = 0.0 ! hail 'echo' volume
      hlnum(kz) = 0.0 ! hail number
      hlnumfd(kz) = 0.0 ! number of hail particles from frozen drops
      hlvol1(kz) = 0.0 ! hail 'echo' volume
      grms(kz) = 0.0  ! graupel mass
      hlms(kz) = 0.0  ! hail mass
      grv(kz) = 0.0  ! graupel particle volume
      hlv(kz) = 0.0  ! hail particle volume
      ticms(kz) = 0.0 ! total ice crystal mass
      rnms(kz) = 0.0  ! total rain mass
      liqms(kz) = 0.0  ! total mass of liquid on ice particles
      swms(kz) = 0.0  ! snow mass
      qcms(kz) = 0.0  ! cloud droplets
      qvms(kz) = 0.0  ! vapor
      iwcmx(kz) = 0.0
      grpotc(kz) = 0.0
      grpotp(kz) = 0.0
      udke(kz)  = 0.0
      ddke(kz)  = 0.0

!      IF ( lf > 1 ) THEN
      fdms(kz) = 0.0
      fdvol(kz) = 0.0
      fdvol1(kz) = 0.0
      fdv(kz)  = 0.0
!      ENDIF
      
      ENDDO

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: start 3D time-height"
#endif

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1
#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: start 3D time-height",ixe,jye,kze
#endif

      DO kz = kzb,kze
#else
      DO kz = 1, nz-1*kd1
#endif

#ifdef MPI
!       jyb = 1
!       jye = jtile
! !      if (jybeg .le. jbsd) jyb = jbsd-jybeg+1
! !      if (jyend .ge. jesd) jye = jesd-jybeg+1
! 
!       ixb = 1
!       ixe = itile
! !      if (ixbeg .le. ibsd) ixb = ibsd-ixbeg+1
! !      if (ixend .ge. iesd) ixe = iesd-ixbeg+1
! 
! 
!       if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag
!       if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag


      DO jy = jyb,jye
      DO ix = ixb,ixe
#else
!      DO jy = 1, ny-1*jd1
!      DO ix = 1, nx-1*id1
      DO jy = jbsd,jesd ! 1,ny-1
      DO ix = ibsd,iesd ! 1,nx-1
#endif
      dxdy = dxx(ix)*dyy(jy)
      dv = dxdy*dzz(kz)

      dv9 = 1.0e-9*dv
      
      IF ( w(ix,jy,kz) > 0.0 ) THEN
        udke(kzbeg-1+kz) = udke(kzbeg-1+kz) + 0.5*dn(ix,jy,kz)*w(ix,jy,kz)**2*dv
      ELSE
        ddke(kzbeg-1+kz) = ddke(kzbeg-1+kz) + 0.5*dn(ix,jy,kz)*w(ix,jy,kz)**2*dv
      ENDIF
      
      
      
      IF ( w(ix,jy,kz) .gt. 3.0 ) THEN
        udmf(kzbeg-1+kz) = udmf(kzbeg-1+kz) + w(ix,jy,kz)*dxdy
      ENDIF
      IF ( w(ix,jy,kz) .gt. 0.5 ) THEN
        udmf1(kzbeg-1+kz) = udmf1(kzbeg-1+kz) + w(ix,jy,kz)*dxdy
      ENDIF
      IF ( w(ix,jy,kz) .gt. 10.0 ) THEN
        udmf10(kzbeg-1+kz) = udmf10(kzbeg-1+kz) + w(ix,jy,kz)*dxdy
      ENDIF

      IF ( w(ix,jy,kz) .le. -5.0 ) THEN
        udv5n = udv5n + dv9 !  1.0/gt(ix,jy,kz,imapz)
      ENDIF
      IF ( w(ix,jy,kz) .ge.  5.0 ) THEN
        udv5 = udv5 + dv9 !  1.0/gt(ix,jy,kz,imapz)
      ENDIF
      IF ( w(ix,jy,kz) .ge.  10.0 ) THEN
        udv10 = udv10 + dv9 !  1.0/gt(ix,jy,kz,imapz)
      ENDIF
      IF ( w(ix,jy,kz) .ge.  20.0 ) THEN
        udv20 = udv20 + dv9 !  1.0/gt(ix,jy,kz,imapz)
      ENDIF

      tmp = 0.0
#ifndef WARM
#ifdef ICE10
      tmp =( an(ix,jy,kz,li)  + an(ix,jy,kz,lir) + an(ix,jy,kz,lip) )
#else
! ELSEIF ( na .ge. 8 ) THEN ! LFO
       tmp = an(ix,jy,kz,li) + an(ix,jy,kz,ls)
!      ELSE
!       tmp = 0.0
!      ENDIF
#endif
#endif
      IF ( lis > 1 ) tmp = tmp + an(ix,jy,kz,lis)
      uicmf(kzbeg-1+kz) = uicmf(kzbeg-1+kz) + Max( 0.0,w(ix,jy,kz) )*tmp*dxdy
      IF ( tmp .ge. 1.0e-5 .or. an(ix,jy,kz,lc) .ge. 1.0e-5 ) THEN
        icldtop = Max(icldtop,kz)
        kcldtop(ix,jy) = Min(icldtop+2,nz-1)

        zcldtop = z1d(icldtop,1)

      ENDIF
      
      tmpi = tmp
      ticms(kzbeg-1+kz) = ticms(kzbeg-1+kz) + tmp*dv
      
      tmpg = 0.0
      IF ( lh .gt. 1 ) THEN

#ifndef WARM
#ifdef ICE10
      tmpg = an(ix,jy,kz,lgl) + an(ix,jy,kz,lgm) + an(ix,jy,kz,lgh) + &
     &        an(ix,jy,kz,lf) 
#elif defined ICE3
!      ELSEIF ( na .ge. 8 ) THEN ! LFO
!      IF ( na .ge. 8 ) THEN
       tmpg = an(ix,jy,kz,lh) 
       IF ( lf > 1 ) tmpg = tmpg + an(ix,jy,kz,lf) 
!      ENDIF
!      ELSE ! Kessler
!      tmpg = 0.0
!      ENDIF
#elif defined HCMON
       tmpg = an(ix,jy,kz,lh) 
       DO i = lh+1,lh+nch-1
        tmpg = tmpg + an(ix,jy,kz,i)
       ENDDO
#elif defined TAKON
      tmpg = an(ix,jy,kz,lh) 
#endif
#endif
      ENDIF
      
      IF ( tmpg .gt. 0.1e-3 ) THEN
        grvol1(kzbeg-1+kz) = grvol1(kzbeg-1+kz) + dv9
        IF ( tmpg .gt. 0.5e-3 ) THEN
         grvol(kzbeg-1+kz) = grvol(kzbeg-1+kz) + dv9
        ENDIF
      ENDIF
      
      grms(kzbeg-1+kz) = grms(kzbeg-1+kz) + tmpg*dv
      
      IF ( lf > 1 ) THEN ! get frozen drop mass
        fdms(kzbeg-1+kz) = fdms(kzbeg-1+kz) + an(ix,jy,kz,lf)*dv
      ELSEIF ( lnhf > 1 ) THEN ! used tracked FD fraction in graupel for FD mass
        IF ( an(ix,jy,kz,lnh) > 1.e-8 ) THEN
          fdms(kzbeg-1+kz) = fdms(kzbeg-1+kz) + an(ix,jy,kz,lh)*dv*Min( 1.0, an(ix,jy,kz,lnhf)/an(ix,jy,kz,lnh) ) 
        ENDIF
      ENDIF
      IF ( lvh .gt. 1 ) grv(kzbeg-1+kz) = grv(kzbeg-1+kz) + an(ix,jy,kz,lvh)*dv
      IF ( lvf .gt. 1 ) grv(kzbeg-1+kz) = grv(kzbeg-1+kz) + an(ix,jy,kz,lvf)*dv

      tmp = 0.0
      tmph = 0.0
      tmphnum = 0.0
      tmphfnum = 0.0
#ifdef ICE10
      IF ( lh .gt. 1 ) THEN ! this is for 10-ice
      tmph = an(ix,jy,kz,lh) + an(ix,jy,kz,lhl) ! + an(ix,jy,kz,lf)
      tmp = tmph
#elif defined ICE3 || defined TAKON
      IF ( lhl .gt. 1 ) THEN ! this is hail in the ZIEG scheme
      tmph = an(ix,jy,kz,lhl)
      IF ( lnhl > 1 ) THEN
        tmphnum = an(ix,jy,kz,lnhl)
        IF ( lnhlf > 1 ) tmphfnum = an(ix,jy,kz,lnhlf)
      ENDIF
      tmp = 0.0
#else
      IF ( .false. ) THEN
#endif
      ELSE 
      tmph = 0.0
      ENDIF
      
      IF ( tmph .gt. 0.1e-3 ) THEN
        hlvol1(kzbeg-1+kz) = hlvol1(kzbeg-1+kz) + dv9 ! 1.0
        IF ( tmph .gt. 0.5e-3 ) THEN
         hlvol(kzbeg-1+kz) = hlvol(kzbeg-1+kz) + dv9 ! 1.0
        ENDIF
      ENDIF
      
      IF ( lnhlf > 1 .and. lnhl > 1 ) THEN
        IF ( an(ix,jy,kz,lnhl) > 1.e-8 ) THEN
        hlfdmasstot = hlfdmasstot + dn(ix,jy,kz)*dv*an(ix,jy,kz,lhl)*Min(1.0, an(ix,jy,kz,lnhlf)/an(ix,jy,kz,lnhl))
        ENDIF

        IF ( kzbeg-1+kz == 1 ) THEN
         hlnumfalltot = hlnumfalltot + xfall(ix,jy,lnhl)*dzz(1)*dxx(ix)*dyy(jy)
         hlnumfdfalltot = hlnumfdfalltot + xfall(ix,jy,lnhlf)*dzz(1)*dxx(ix)*dyy(jy)
        ENDIF
      ENDIF

      IF ( lvhl .gt. 1 ) hlv(kzbeg-1+kz) = hlv(kzbeg-1+kz) + an(ix,jy,kz,lvhl)*dv
      
      hlms(kzbeg-1+kz) = hlms(kzbeg-1+kz) + tmph*dv
      
      iwcmx(kzbeg-1+kz) = Max( iwcmx(kzbeg-1+kz), tmpg + tmph )

      qvms(kzbeg-1+kz) = qvms(kzbeg-1+kz) + an(ix,jy,kz,lv)*dv
      
      qcms(kzbeg-1+kz) = qcms(kzbeg-1+kz) + an(ix,jy,kz,lc)*dv

      rnms(kzbeg-1+kz) = rnms(kzbeg-1+kz) + an(ix,jy,kz,lr)*dv
      
      hlnum(kzbeg-1+kz)   = hlnum(kzbeg-1+kz) + tmphnum*dv
      hlnumfd(kzbeg-1+kz) = hlnumfd(kzbeg-1+kz) + tmphfnum*dv

      IF ( lsw > 1 )  liqms(kzbeg-1+kz) = liqms(kzbeg-1+kz) + an(ix,jy,kz,lsw)*dv
      IF ( lhw > 1 )  liqms(kzbeg-1+kz) = liqms(kzbeg-1+kz) + an(ix,jy,kz,lhw)*dv
      IF ( lhlw > 1 ) liqms(kzbeg-1+kz) = liqms(kzbeg-1+kz) + an(ix,jy,kz,lhlw)*dv
      
      
      IF ( ls .gt. 1 ) THEN
        swms(kzbeg-1+kz) = swms(kzbeg-1+kz) + an(ix,jy,kz,ls)*dv
        tmps = an(ix,jy,kz,ls)
      ELSE
        swms(kzbeg-1+kz) = 0.0
        tmps = 0.0
      ENDIF
      
      grpotc(kzbeg-1+kz) = grpotc(kzbeg-1+kz) + g*z1d(kzbeg-1+kz,1)*dv*(an(ix,jy,kz,lc) + tmpi)
      grpotp(kzbeg-1+kz) = grpotp(kzbeg-1+kz) + g*z1d(kzbeg-1+kz,1)*dv*(tmps + tmpg + tmph + an(ix,jy,kz,lr) )

#ifdef MPI
!      if (debug_mpi) write(0,*) "DRIVER: DEBUG 5"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 5"
#endif

!  Microphysics processes domain max/min/avg calculations
      
       if (tq1(ix,jy,kz,1) .ne. 0.) then
        qproczmx(kzbeg-1+kz) = Max( tq1(ix,jy,kz,1), qproczmx(kzbeg-1+kz) )  !aps
        qproczmn(kzbeg-1+kz) = Min( tq1(ix,jy,kz,1), qproczmn(kzbeg-1+kz) )  !aps
      
       if (qproczmx(kzbeg-1+kz) .ne. 0. .or. qproczmn(kzbeg-1+kz) .ne. 1.) then
        iavgcnt = iavgcnt + 1
        qproczav(kzbeg-1+kz) = qproczav(kzbeg-1+kz) + tq1(ix,jy,kz,1)
       end if 
          
! compute the average at the last point in the loop 
#ifdef MPI
        if (ixbeg-1+ix .eq. nxend-1) then
#else
        if (ix .eq. nx-1*id1) then
#endif
!         qproczav(kzbeg-1+kz) = qproczav(kzbeg-1+kz) / iavgcnt !aps
        endif
       
      end if !if tq1 > 0
      
#ifdef CHGELEC
! Charging totals      
      
      IF ( ipelec .gt. 0 ) THEN
      gt1 = dxx(ix)*dyy(jy)*dzz(kz)

      cswizmn(kzbeg-1+kz) = Min( cswizmn(kzbeg-1+kz), t1(ix,jy,kz) )
      cswizmx(kzbeg-1+kz) = Max( cswizmx(kzbeg-1+kz), t1(ix,jy,kz) )

      cswwzmn(kzbeg-1+kz) = Min( cswwzmn(kzbeg-1+kz), t4(ix,jy,kz) )
      cswwzmx(kzbeg-1+kz) = Max( cswwzmx(kzbeg-1+kz), t4(ix,jy,kz) )

      ctswin = ctswin + Min( 0.0, t1(ix,jy,kz) )*gt1
      ctswip = ctswip + Max( 0.0, t1(ix,jy,kz) )*gt1

      ctswwn = ctswwn + Min( 0.0, t4(ix,jy,kz) )*gt1
      ctswwp = ctswwp + Max( 0.0, t4(ix,jy,kz) )*gt1

      cthsdn = cthsdn + Min( 0.0, tq1(ix,jy,kz,3) )*gt1
      cthsdp = cthsdp + Max( 0.0, tq1(ix,jy,kz,3) )*gt1


!      cthsdn = cthsdn + Min( 0.0, t4(ix,jy,kz) )*gt1
!      cthsdp = cthsdp + Max( 0.0, t4(ix,jy,kz) )*gt1
      cghszmn(kzbeg-1+kz) = Min( cghszmn(kzbeg-1+kz), t2(ix,jy,kz) + t6(ix,jy,kz) )
      cghszmx(kzbeg-1+kz) = Max( cghszmx(kzbeg-1+kz), t2(ix,jy,kz) + t6(ix,jy,kz) )
      ctghsn = ctghsn + Min( 0.0, t2(ix,jy,kz) + t6(ix,jy,kz) )*gt1
      ctghsp = ctghsp + Max( 0.0, t2(ix,jy,kz) + t6(ix,jy,kz) )*gt1
      cghizmn(kzbeg-1+kz) = Min( cghizmn(kzbeg-1+kz), t3(ix,jy,kz) + t8(ix,jy,kz) )
      cghizmx(kzbeg-1+kz) = Max( cghizmx(kzbeg-1+kz), t3(ix,jy,kz) + t8(ix,jy,kz) )
      ctghin = ctghin + Min( 0.0, t3(ix,jy,kz) + t8(ix,jy,kz) )*gt1
      ctghip = ctghip + Max( 0.0, t3(ix,jy,kz) + t8(ix,jy,kz) )*gt1
      cghwzmn(kzbeg-1+kz) = Min( cghwzmn(kzbeg-1+kz), t5(ix,jy,kz) )
      cghwzmx(kzbeg-1+kz) = Max( cghwzmx(kzbeg-1+kz), t5(ix,jy,kz) )
      ctghwn = ctghwn + Min( 0.0, t5(ix,jy,kz) )*gt1
      ctghwp = ctghwp + Max( 0.0, t5(ix,jy,kz) )*gt1
      cgszmn(kzbeg-1+kz) = Min( cgszmn(kzbeg-1+kz), t2(ix,jy,kz) )
      cgszmx(kzbeg-1+kz) = Max( cgszmx(kzbeg-1+kz), t2(ix,jy,kz) )
      ctgsn = ctgsn + Min( 0.0, t2(ix,jy,kz) )*gt1
      ctgsp = ctgsp + Max( 0.0, t2(ix,jy,kz) )*gt1
      cgizmn(kzbeg-1+kz) = Min( cgizmn(kzbeg-1+kz), t3(ix,jy,kz) )
      cgizmx(kzbeg-1+kz) = Max( cgizmx(kzbeg-1+kz), t3(ix,jy,kz) )
      ctgin = ctgin + Min( 0.0, t3(ix,jy,kz) )*gt1
      ctgip = ctgip + Max( 0.0, t3(ix,jy,kz) )*gt1
      chszmn(kzbeg-1+kz) = Min( chszmn(kzbeg-1+kz), t6(ix,jy,kz) )
      chszmx(kzbeg-1+kz) = Max( chszmx(kzbeg-1+kz), t6(ix,jy,kz) )
      cthsn = cthsn + Min( 0.0, t6(ix,jy,kz) )*gt1
      cthsp = cthsp + Max( 0.0, t6(ix,jy,kz) )*gt1
      chizmn(kzbeg-1+kz) = Min( chizmn(kzbeg-1+kz), t8(ix,jy,kz) )
      chizmx(kzbeg-1+kz) = Max( chizmx(kzbeg-1+kz), t8(ix,jy,kz) )
      cthin = cthin + Min( 0.0, t8(ix,jy,kz) )*gt1
      cthip = cthip + Max( 0.0, t8(ix,jy,kz) )*gt1
      
      ctswinz(kzbeg-1+kz) = ctswinz(kzbeg-1+kz) +  Min( 0.0, t1(ix,jy,kz) )*gt1
      ctswipz(kzbeg-1+kz) = ctswipz(kzbeg-1+kz) +  Max( 0.0, t1(ix,jy,kz) )*gt1
      ctghsnz(kzbeg-1+kz) = ctghsnz(kzbeg-1+kz) + Min( 0.0, t2(ix,jy,kz) + t6(ix,jy,kz) )*gt1
      ctghspz(kzbeg-1+kz) = ctghspz(kzbeg-1+kz) + Max( 0.0, t2(ix,jy,kz) + t6(ix,jy,kz) )*gt1
      ctghinz(kzbeg-1+kz) = ctghinz(kzbeg-1+kz) + Min( 0.0, t3(ix,jy,kz) + t8(ix,jy,kz) )*gt1
      ctghipz(kzbeg-1+kz) = ctghipz(kzbeg-1+kz) + Max( 0.0, t3(ix,jy,kz) + t8(ix,jy,kz) )*gt1
      ctghwnz(kzbeg-1+kz) = ctghwnz(kzbeg-1+kz) + Min( 0.0, t5(ix,jy,kz) )*gt1
      ctghwpz(kzbeg-1+kz) = ctghwpz(kzbeg-1+kz) + Max( 0.0, t5(ix,jy,kz) )*gt1
      ctgsnz(kzbeg-1+kz) = ctgsnz(kzbeg-1+kz) +  Min( 0.0, t2(ix,jy,kz) )*gt1
      ctgspz(kzbeg-1+kz) = ctgspz(kzbeg-1+kz) +  Max( 0.0, t2(ix,jy,kz) )*gt1
      ctginz(kzbeg-1+kz) = ctginz(kzbeg-1+kz) +  Min( 0.0, t3(ix,jy,kz) )*gt1
      ctgipz(kzbeg-1+kz) = ctgipz(kzbeg-1+kz) +  Max( 0.0, t3(ix,jy,kz) )*gt1
      cthsnz(kzbeg-1+kz) = cthsnz(kzbeg-1+kz) +  Min( 0.0, t6(ix,jy,kz) )*gt1
      cthspz(kzbeg-1+kz) = cthspz(kzbeg-1+kz) +  Max( 0.0, t6(ix,jy,kz) )*gt1
      cthinz(kzbeg-1+kz) = cthinz(kzbeg-1+kz) +  Min( 0.0, t8(ix,jy,kz) )*gt1
      cthipz(kzbeg-1+kz) = cthipz(kzbeg-1+kz) +  Max( 0.0, t8(ix,jy,kz) )*gt1
      ENDIF ! ipelec
#endif
      END DO
      END DO

      dbt = db(kz)
      
      udmf(kzbeg-1+kz) = udmf(kzbeg-1+kz)*dbt! *dx*dy   ! updraft mass flux
      udmf1(kzbeg-1+kz) = udmf1(kzbeg-1+kz)*dbt! *dx*dy   ! updraft mass flux1
      udmf10(kzbeg-1+kz) = udmf10(kzbeg-1+kz)*dbt! *dx*dy   ! updraft mass flux1
      uicmf(kzbeg-1+kz) = uicmf(kzbeg-1+kz)*dbt! *dx*dy ! upward ice crystal mass flux
!      grvol(kzbeg-1+kz) = grvol(kzbeg-1+kz) !*dx*dy*dz/gt(1,1,kz,imapz)*1.0e-9 ! graupel volume (km**3)
!      grvol1(kzbeg-1+kz) = grvol1(kzbeg-1+kz) !*dx*dy*dz/gt(1,1,kz,imapz)*1.0e-9 ! graupel volume (km**3)
!      hlvol(kzbeg-1+kz) = hlvol(kzbeg-1+kz) !*dx*dy*dz/gt(1,1,kz,imapz)*1.0e-9 ! hail volume (km**3)

      IF ( lnhlf > 1 ) THEN
        hlnum(kzbeg-1+kz) = hlnum(kzbeg-1+kz) !*dbt 
      ENDIF

      IF ( lf > 1 .or. lnhf > 1 ) fdms(kzbeg-1+kz) = fdms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! graupel mass
      grms(kzbeg-1+kz) = grms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! graupel mass
      hlms(kzbeg-1+kz) = hlms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! hail mass
      ticms(kzbeg-1+kz) = ticms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz) ! total ice crystal mass
      rnms(kzbeg-1+kz) = rnms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! rain mass
      liqms(kzbeg-1+kz) = liqms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! rain mass
      swms(kzbeg-1+kz) = swms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! snow mass
      qcms(kzbeg-1+kz) = qcms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! droplet mass
      qvms(kzbeg-1+kz) = qvms(kzbeg-1+kz)*dbt ! *dx*dy*dz/gt(1,1,kz,imapz)  ! vapor mass
      iwcmx(kzbeg-1+kz) = 1.e3*iwcmx(kzbeg-1+kz)*dbt ! convert to g/m**3
      grpotc(kzbeg-1+kz) = grpotc(kzbeg-1+kz)*dbt
      grpotp(kzbeg-1+kz) = grpotp(kzbeg-1+kz)*dbt

      END DO

#ifndef WARM

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: cloud-top all_reduce"
! find global integrated rate max
       mpitotin(1)  =  zcldtop
       IF ( icldtop > 0 ) THEN
       mpitotin(2)  =  icldtop+kzbeg-1
       ELSE
       mpitotin(2)  =  0
       ENDIF
       
       n = 2

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)

       
       zcldtop = mpitotout(1)  
       icldtop = Nint( mpitotout(2) )


#endif

      IF ( my_rank == 0 ) THEN
        write(iunit,*) 'icldtop,zcldtop = ',icldtop, zcldtop
      ENDIF
#endif
      
       udketot = 0
       ddketot = 0
       grpotctot = 0
       grpotptot = 0
       grvoltot = 0.0
       grvoltot1 = 0.0
       hlvoltot = 0.0
       hlvoltot1 = 0.0
       hlnumtot = 0.0
       hlnumfdtot = 0.0
       grvoltota = 0.0
       hlvoltota = 0.0
       grvoltotb = 0.0
       hlvoltotb = 0.0
       grmstot = 0.0
       hlmstot = 0.0
       hlnumtot = 0.0
       hlnumfdtot = 0.0
       grmstota = 0.0
       hlmstota = 0.0
       grvtota = 0.0
       hlvtota = 0.0
       grmstotb = 0.0
       hlmstotb = 0.0
       ticmstot = 0.0
       rnmstot = 0.0
       rnmstota = 0.0
       rnmstotb = 0.0
       liqtot = 0.0d0
       swmstot = 0.0
       qctot = 0.0
       qctota = 0.0
       qctotb = 0.0
       qvtot = 0.0
       iqproc = 0

       fdvoltot = 0.0
       fdvoltot1 = 0.0
       fdmstot = 0.0
       fdmstota = 0.0
       fdmstotb = 0.0

       proctot(:) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      DO kz = kzb,kze
#else
      DO kz = 1, nz-1*kd1
#endif
       grvoltot = grvoltot + grvol(kzbeg-1+kz)
       grvoltot1 = grvoltot1 + grvol1(kzbeg-1+kz)
       hlvoltot = hlvoltot + hlvol(kzbeg-1+kz)
       hlvoltot1 = hlvoltot1 + hlvol1(kzbeg-1+kz)
       grmstot = grmstot + grms(kzbeg-1+kz)
       hlmstot = hlmstot + hlms(kzbeg-1+kz)
       hlnumtot = hlnumtot + hlnum(kzbeg-1+kz)
       hlnumfdtot = hlnumfdtot + hlnumfd(kzbeg-1+kz)
       IF ( lf > 1 .or. lnhf > 1 ) fdmstot = fdmstot + fdms(kzbeg-1+kz)
       ticmstot = ticmstot + ticms(kzbeg-1+kz)
       rnmstot = rnmstot + rnms(kzbeg-1+kz)
       swmstot = swmstot + swms(kzbeg-1+kz)
       liqtot = liqtot + liqms(kzbeg-1+kz)
       qctot = qctot + qcms(kzbeg-1+kz)
       qvtot = qvtot + qvms(kzbeg-1+kz)
       udketot = udketot + udke(kzbeg-1+kz)
       ddketot = ddketot + ddke(kzbeg-1+kz)
       grpotctot = grpotctot + grpotc(kzbeg-1+kz)
       grpotptot = grpotptot + grpotp(kzbeg-1+kz)
       
       DO i = 1,numproc
         proctot(i) = proctot(i) + thproc(kzbeg-1+kz,i)
       ENDDO
       
       IF ( kzbeg-1+kz .gt. lt273 ) THEN
         grvoltota = grvoltota + grvol(kzbeg-1+kz)
         hlvoltota = hlvoltota + hlvol(kzbeg-1+kz)
         grmstota = grmstota + grms(kzbeg-1+kz)
         hlmstota = hlmstota + hlms(kzbeg-1+kz)
         IF ( lf > 1 .or. lnhf > 1 ) fdmstota = fdmstota + fdms(kzbeg-1+kz)
         grvtota = grvtota + grv(kzbeg-1+kz)
         hlvtota = hlvtota + hlv(kzbeg-1+kz)
         rnmstota = rnmstota + rnms(kzbeg-1+kz)
         qctota = qctota + qcms(kzbeg-1+kz)
       ELSE
         grvoltotb = grvoltotb + grvol(kzbeg-1+kz)
         hlvoltotb = hlvoltotb + hlvol(kzbeg-1+kz)
         grmstotb = grmstotb + grms(kzbeg-1+kz)
         hlmstotb = hlmstotb + hlms(kzbeg-1+kz)
         IF ( lf > 1 .or. lnhf > 1 ) fdmstotb = fdmstotb + fdms(kzbeg-1+kz)
         rnmstotb = rnmstotb + rnms(kzbeg-1+kz)
         qctotb = qctotb + qcms(kzbeg-1+kz)
       ENDIF
       
      qprocmax = Max( qprocmax, qproczmx(kzbeg-1+kz) )
      qprocmin = Min( qprocmin, qproczmn(kzbeg-1+kz) )

      if (qproczav(kzbeg-1+kz) .ne. 0.) then
       iavgcnt = iavgcnt + 1
       qprocavg = qprocavg + qproczav(kzbeg-1+kz)
      end if
      
      if (qproczav(kzbeg-1+kz) .ge. qprocthr) then
       iqproc = iqproc + 1
      end if
      
#ifdef CHGELEC
      cswimn = Min( cswimn, cswizmn(kzbeg-1+kz) )
      cswimx = Max( cswimx, cswizmx(kzbeg-1+kz) )
      cswwmn = Min( cswwmn, cswwzmn(kzbeg-1+kz) )
      cswwmx = Max( cswwmx, cswwzmx(kzbeg-1+kz) )
      cghsmn = Min( cghsmn, cghszmn(kzbeg-1+kz) )
      cghsmx = Max( cghsmx, cghszmx(kzbeg-1+kz) )
      cghimn = Min( cghimn, cghizmn(kzbeg-1+kz) )
      cghimx = Max( cghimx, cghizmx(kzbeg-1+kz) )
      cghwmn = Min( cghwmn, cghwzmn(kzbeg-1+kz) )
      cghwmx = Max( cghwmx, cghwzmx(kzbeg-1+kz) )
      cgsmn = Min( cgsmn, cgszmn(kzbeg-1+kz) )
      cgsmx = Max( cgsmx, cgszmx(kzbeg-1+kz) )
      cgimn = Min( cgimn, cgizmn(kzbeg-1+kz) )
      cgimx = Max( cgimx, cgizmx(kzbeg-1+kz) )
      chsmn = Min( chsmn, chszmn(kzbeg-1+kz) )
      chsmx = Max( chsmx, chszmx(kzbeg-1+kz) )
      chimn = Min( chimn, chizmn(kzbeg-1+kz) )
      chimx = Max( chimx, chizmx(kzbeg-1+kz) )
#endif
      END DO
      
!Print domain max/min/avg for microphysics processes

#ifdef MPI
     

! find global charging rate max
       mpitotin(1)  = qprocavg
       mpitotin(2)  = dble(iavgcnt)

      CALL MPI_Reduce(mpitotin, mpitotout, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
       qprocavg = mpitotout(1)  
       iavgcnt = Nint(mpitotout(2))
      ENDIF

       mpitotin(1)  = qprocmax

      CALL MPI_Reduce(mpitotin, mpitotout, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
       qprocmax = mpitotout(1)  
      ENDIF

       mpitotin(1)  = qprocmin

      CALL MPI_Reduce(mpitotin, mpitotout, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
       qprocmin = mpitotout(1)  
      ENDIF


      
#endif
      
      if (iavgcnt .ne. 0) qprocavg = qprocavg / iavgcnt
      
      IF ( my_rank == 0 ) THEN
        write(iunit, *) 'qproc max =', qprocmax
!         write(iunit, *) 'qproc max,min,avg=', qprocmax,qprocmin,qprocavg
!         write(iunit, *) 'iavgcnt = ',iavgcnt
!         write(iunit, *) 'qavg exceeds thresh',iqproc,'number of times'
      ENDIF
!      DO kz = 1, nz-1
!      DO jy = 1, ny-1
!      DO ix = 1, nx-1
!       if (qprocmax .eq. tq1(ix,jy,kz) .and. qprocmax .ne. 0.) write(iunit, *) 'qmax at i j k ',ix,jy,kz
!       if (qprocmin .eq. tq1(ix,jy,kz) .and. qprocmin .ne. 1.) write(iunit, *) 'qmin at i j k ',ix,jy,kz
!      END DO
!      END DO
!      END DO

#ifdef MPI
#ifdef CHGELEC
     
      IF ( ipelec .ge. 1 ) THEN

! find global charging rate max
       mpitotin(1)  = cswimx
       mpitotin(2)  = cghsmx
       mpitotin(3)  = cghimx
       mpitotin(4)  = cghwmx
       mpitotin(5)  = cgsmx
       mpitotin(6)  = cgimx
       mpitotin(7)  = chsmx
       mpitotin(8)  = chimx
       mpitotin(9)  = cswwmx

      CALL MPI_Reduce(mpitotin, mpitotout, 9, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       cswimx = mpitotout(1)  
       cghsmx = mpitotout(2)  
       cghimx = mpitotout(3)  
       cghwmx = mpitotout(4)  
       cgsmx = mpitotout(5)  
       cgimx = mpitotout(6)  
       chsmx = mpitotout(7)  
       chimx = mpitotout(8) 
       cswwmx = mpitotout(9)  
      ENDIF

! find global charging rate min:

       mpitotin(1)  = cswimn
       mpitotin(2)  = cghsmn
       mpitotin(3)  = cghimn
       mpitotin(4)  = cghwmn
       mpitotin(5)  = cgsmn
       mpitotin(6)  = cgimn
       mpitotin(7)  = chsmn
       mpitotin(8)  = chimn
       mpitotin(9)  = cswwmn

      CALL MPI_Reduce(mpitotin, mpitotout, 9, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       cswimn = mpitotout(1)  
       cghsmn = mpitotout(2)  
       cghimn = mpitotout(3)  
       cghwmn = mpitotout(4)  
       cgsmn = mpitotout(5)  
       cgimn = mpitotout(6)  
       chsmn = mpitotout(7)  
       chimn = mpitotout(8) 
       cswwmn = mpitotout(9)  
      ENDIF

! find global integrated rate totals
       mpitotin(1)  = ctswin
       mpitotin(2)  = ctghsn
       mpitotin(3)  = ctghin
       mpitotin(4)  = ctghwn
       mpitotin(5)  = ctgsn
       mpitotin(6)  = ctgin
       mpitotin(7)  = cthsn
       mpitotin(8)  = cthin
       
       mpitotin(9)  = ctswip
       mpitotin(10)  = ctghsp
       mpitotin(11)  = ctghip
       mpitotin(12)  = ctghwp
       mpitotin(13)  = ctgsp
       mpitotin(14)  = ctgip
       mpitotin(15)  = cthsp
       mpitotin(16)  = cthip

       mpitotin(17)  = cthsdn
       mpitotin(18)  = cthsdp

       mpitotin(19)  = ctswwn
       mpitotin(20)  = ctswwp

      CALL MPI_Reduce(mpitotin, mpitotout, 20, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       ctswin = mpitotout(1)  
       ctghsn = mpitotout(2)  
       ctghin = mpitotout(3)  
       ctghwn = mpitotout(4)  
       ctgsn = mpitotout(5)  
       ctgin = mpitotout(6)  
       cthsn = mpitotout(7)  
       cthin = mpitotout(8) 

       ctswip = mpitotout(9)  
       ctghsp = mpitotout(10)  
       ctghip = mpitotout(11)  
       ctghwp = mpitotout(12)  
       ctgsp = mpitotout(13)  
       ctgip = mpitotout(14)  
       cthsp = mpitotout(15)  
       cthip = mpitotout(16) 
       cthsdn = mpitotout(17)
       cthsdp = mpitotout(18)
       ctswwn = mpitotout(19)
       ctswwp = mpitotout(20)
      ENDIF

      ENDIF ! ipelec
! endif for chgelec:
#endif /* chgelec */
      
      allocate( mpitotinth((nzend),6))
      allocate(mpitotoutth((nzend),6))
      
      mpitotinth(:,:) = 0.0d0
      
      DO k = kzbeg,kzend
       mpitotinth(k, 1) = udmf(k)
       mpitotinth(k, 2) = udmf1(k)
       mpitotinth(k, 3) = udmf10(k)
       mpitotinth(k, 4) = uicmf(k)
       mpitotinth(k, 5) = udke(k)
       mpitotinth(k, 6) = ddke(k)
      ENDDO

      n = 6*(nzend)

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
       DO k = 1,nzend
        udmf(k)   =  mpitotoutth(k, 1)
        udmf1(k)  =  mpitotoutth(k, 2)
        udmf10(k) =  mpitotoutth(k, 3)
        uicmf(k)  =  mpitotoutth(k, 4)
        udke(k)   =  mpitotoutth(k, 5)
        ddke(k)   =  mpitotoutth(k, 6)
       ENDDO
      ENDIF
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )
      

#endif
      
#ifdef CHGELEC
      IF ( ipelec .ge. 1 ) THEN

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a)') 'gridpoint min/max charging rates'
      write(iunit,'(2(a,1pe12.5))') 'cswimn,cswimx = ',cswimn,', ',cswimx
      write(iunit,'(2(a,1pe12.5))') 'cghsmn,cghsmx = ',cghsmn,', ',cghsmx
      write(iunit,'(2(a,1pe12.5))') 'cghimn,cghimx = ',cghimn,', ',cghimx
      write(iunit,'(2(a,1pe12.5))') 'cghwmn,cghwmx = ',cghwmn,', ',cghwmx
      write(iunit,'(2(a,1pe12.5))') 'cgsmn,cgsmx = ',cgsmn,', ',cgsmx
      write(iunit,'(2(a,1pe12.5))') 'cgimn,cgimx = ',cgimn,', ',cgimx
      write(iunit,'(2(a,1pe12.5))') 'chsmn,chsmx = ',chsmn,', ',chsmx
      write(iunit,'(2(a,1pe12.5))') 'chimn,chimx = ',chimn,', ',chimx
      write(iunit,'(2(a,1pe12.5))') 'cswwmn,cswwmx = ',cswwmn,', ',cswwmx

! write out the integrated charge rates
      write(iunit,'(a)') 'Integrated pos/neg charging rates:'
      write(iunit,'(2(a,1pe12.5))') 'ctswin,ctswip = ',ctswin,', ',ctswip
      write(iunit,'(2(a,1pe12.5))') 'ctghsn,ctghsp = ',ctghsn,', ',ctghsp
      write(iunit,'(2(a,1pe12.5))') 'ctghin,ctghip = ',ctghin,', ',ctghip
      write(iunit,'(2(a,1pe12.5))') 'ctghwn,ctghwp = ',ctghwn,', ',ctghwp
      write(iunit,'(2(a,1pe12.5))') 'ctgsn,ctgsp = ',ctgsn,', ',ctgsp
      write(iunit,'(2(a,1pe12.5))') 'ctgin,ctgip = ',ctgin,', ',ctgip
      write(iunit,'(2(a,1pe12.5))') 'cthsn,cthsp = ',cthsn,', ',cthsp
      write(iunit,'(2(a,1pe12.5))') 'cthin,cthip = ',cthin,', ',cthip
      write(iunit,'(2(a,1pe12.5))') 'cthsdn,cthsdp = ',cthsdn,', ',cthsdp
      write(iunit,'(2(a,1pe12.5))') 'ctswwn,ctswwp = ',ctswwn,', ',ctswwp

      ENDIF
      
      ENDIF ! ipelec

! endif for chgelec:
#endif /* chgelec */

#ifdef MPI

       mpitotin(1)  = qctot
       mpitotin(2)  = qvtot
       mpitotin(3)  = grvoltot
       mpitotin(4)  = hlvoltot
       mpitotin(5)  = grvoltot1
       mpitotin(6)  = hlvoltot1
       mpitotin(7)  = grvoltota
       mpitotin(8)  = hlvoltota
       mpitotin(9)  = grvoltotb
       mpitotin(10) = hlvoltotb
       mpitotin(11) = grmstot
       mpitotin(12) = hlmstot
       mpitotin(13) = ticmstot
       mpitotin(14) = grmstota
       mpitotin(15) = hlmstota
       mpitotin(16) = grmstotb
       mpitotin(17) = hlmstotb
       mpitotin(18) = swmstot
       mpitotin(19) = rnmstot
       mpitotin(20) = rnmstota
       mpitotin(21) = rnmstotb
       mpitotin(22) = udv5n
       mpitotin(23) = udv5
       mpitotin(24) = udv10
       mpitotin(25) = udv20
       mpitotin(26) = grvtota
       mpitotin(27) = hlvtota
       mpitotin(28) = qctota
       mpitotin(29) = qctotb
       mpitotin(30) = liqtot
       mpitotin(31) = udketot
       mpitotin(32) = ddketot
       mpitotin(33) = grpotctot
       mpitotin(34) = grpotptot
       mpitotin(35) = fdmstot
       mpitotin(36) = fdmstota
       mpitotin(37) = fdmstotb
       mpitotin(38) = hlnumtot
       mpitotin(39) = hlnumfdtot
       mpitotin(40) = hlnumfalltot
       mpitotin(41) = hlnumfdfalltot
       mpitotin(42) = hlfdmasstot
!mpitotin(ntot), mpitotout(ntot)

       k = 42

       DO i = 1,numproc ! k+1,k+numproc
         mpitotin(k+i) = proctot(i)
       ENDDO

       n = k+numproc
       IF ( n > ntot ) THEN
         IF( my_rank == 0 ) write(0,*) 'micro_driver: increase size of ntot! ntot,n = ',ntot,n
         call commasmpi_abort()
       ENDIF

!      CALL MPI_Allreduce(mpitotin, mpitotout, ntot, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       IF ( my_rank == 0 ) THEN
        qctot     =   mpitotout(1) 
        qvtot     =   mpitotout(2) 
        grvoltot  =   mpitotout(3) 
        hlvoltot  =   mpitotout(4) 
        grvoltot1 =   mpitotout(5) 
        hlvoltot1 =   mpitotout(6) 
        grvoltota =   mpitotout(7) 
        hlvoltota =   mpitotout(8) 
        grvoltotb =   mpitotout(9) 
        hlvoltotb =   mpitotout(10)
        grmstot   =   mpitotout(11)
        hlmstot   =   mpitotout(12)
        ticmstot  =   mpitotout(13)
        grmstota  =   mpitotout(14)
        hlmstota  =   mpitotout(15)
        grmstotb  =   mpitotout(16)
        hlmstotb  =   mpitotout(17)
        swmstot   =   mpitotout(18)
        rnmstot   =   mpitotout(19)
        rnmstota  =   mpitotout(20)
        rnmstotb  =   mpitotout(21)
        udv5n     =   mpitotout(22)
        udv5      =   mpitotout(23)
        udv10     =   mpitotout(24)
        udv20     =   mpitotout(25)
        grvtota   =   mpitotout(26)
        hlvtota   =   mpitotout(27)
        qctota    =   mpitotout(28)
        qctotb    =   mpitotout(29)
        liqtot    =   mpitotout(30)
        udketot   =   mpitotout(31)
        ddketot   =   mpitotout(32)
        grpotctot =   mpitotout(33)
        grpotptot =   mpitotout(34)
        fdmstot   =   mpitotout(35)
        fdmstota  =   mpitotout(36)
        fdmstotb  =   mpitotout(37)
        hlnumtot   =  mpitotout(38)
        hlnumfdtot =  mpitotout(39)
        hlnumfalltot   =  mpitotout(40)
        hlnumfdfalltot =  mpitotout(41)
        hlfdmasstot = mpitotout(42)
        DO i = 1,numproc ! k+1,k+numproc
          proctot(i) = mpitotout(k+i)
        ENDDO
       ENDIF

      IF ( iraintypes > 0 ) THEN
        DO i = 1,nraintypes
          mpitotin(i) = rate2dtot(i)
        ENDDO

      n = nraintypes
      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
        DO i = 1,nraintypes
          rate2dtot(i) = mpitotout(i)
        ENDDO
      
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
#ifndef WARM
     IF ( lf > 0 .or. lnhf > 1 ) THEN
!      write(iunit,'(a,2(1x,1pe12.5))') 'frozen drop volumes',    &
!     &  fdvoltot, fdvoltot1
      write(iunit,'(a,3(1x,1pe12.5))')    &
     & 'frozen drop masses tot/above/below 0C',   & 
     &  fdmstot, fdmstota, fdmstotb
     
     ENDIF
     
     IF ( lnhlf > 1 ) THEN
      write(iunit,'(a,4(1x,1pe12.5))') 'Hail/HailFD number, fallout',    &
     &  hlnumtot, hlnumfdtot, hlnumfalltot, hlnumfdfalltot
      write(iunit,'(a,2(1x,1pe12.5))') 'HailFD mass ',hlfdmasstot,hlmstot
     ENDIF

      write(iunit,'(a,4(1x,1pe12.5))') 'graupel, hail volumes',    &
     &  grvoltot, hlvoltot, grvoltot1, hlvoltot1
      write(iunit,'(a,4(1x,1pe12.5))')    &
     & 'graupel/hail volumes above/below 0C',   &
     &  grvoltota, hlvoltota, grvoltotb, hlvoltotb
      write(iunit,'(a,3(1x,1pe12.5))') 'graupel, hail, ice masses',    &
     &  grmstot, hlmstot, ticmstot
      write(iunit,'(a,4(1x,1pe12.5))')    &
     & 'graupel/hail masses above/below 0C',   & 
     &  grmstota, hlmstota, grmstotb, hlmstotb
      tmp = 0.0
      tmpg = 0.0
      IF ( grvtota > 0.0 ) tmp = grmstota/grvtota
      IF ( hlvtota > 0.0 ) tmpg = hlmstota/hlvtota
      write(iunit,'(a,4(1x,1pe12.5))')    &
     & 'graupel/hail particle volumes/densities above 0C',    &
     &  grvtota, hlvtota, tmp, tmpg ! grmstota/(grvtota + 1.e-30), hlmstota/(hlvtota + 1.e-30)
      write(iunit,'(a,1(1x,1pe12.5))') 'snow mass', swmstot
#endif
      write(iunit,'(a,3(1x,1pe12.5))')    &
     & 'rain masses tot,above/below 0C',    &
     &  rnmstot, rnmstota, rnmstotb
      write(iunit,'(a,(1x,1pe12.5))')    &
     & 'total liquid on ice ',    &
     &  liqtot
     
      ENDIF
     
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,5(1x,1pe12.5))') 'udmf at 0,-10,-20,and -30:',   &
     &   udmf(lt273),udmf(lt263),udmf(lt253),udmf(lt243),udmf(lt233)
      write(iunit,'(a,5(1x,1pe12.5))') 'udmf1 at 0,-10,-20,and -30:',   &
     &   udmf1(lt273),udmf1(lt263),udmf1(lt253),udmf1(lt243),udmf1(lt233)
      write(iunit,'(a,5(1x,1pe12.5))') 'udmf10 at 0,-10,-20,and -30:',   &
     &   udmf10(lt273),udmf10(lt263),udmf10(lt253),udmf10(lt243),udmf10(lt233)
#ifndef WARM
      write(iunit,'(a,3(1x,1pe12.5))') 'uicmf at -10, -20, -30:',   &
     &   uicmf(lt263),uicmf(lt253),uicmf(lt243)
#endif

      write(iunit,'(a,4(1x,1pe12.5))') 'udv  -5, 5, 10, 20:',   &
     &   udv5n, udv5, udv10, udv20
!     &   udv5n*dx*dy*dz, udv5*dx*dy*dz, udv10*dx*dy*dz, udv20*dx*dy*dz
      write(iunit,'(a,2(1x,1pe12.5))') 'qctot, qvtot', qctot, qvtot
      write(iunit,'(a,2(1x,1pe12.5))') 'qc mass above/below 0C', qctota, qctotb

      write(iunit,'(a,4(1x,1pe12.5))') 'udketot, ddketot :',udketot, ddketot
      write(iunit,'(a,4(1x,1pe12.5))') 'grpotctot, grpotptot :',grpotctot, grpotptot
!      write(iunit,'(a,4(1x,1pe12.5))') 'grpotc, grpotp :',grpotc, grpotp


#ifdef WARM
      write(iunit,'(a,1x,1pe12.5)') 'Total mass: ', qvtot+qctot+rnmstot
#endif

#if defined ( ICE3 ) || defined ( HCMON )  || defined(TAKON)
      write(iunit,'(a,1x,1pe12.5)') 'Total mass: ',   &
     & qvtot+qctot+rnmstot+ticmstot+swmstot+hlmstot
#endif

#ifdef ICE10
      write(iunit,'(a,1x,1pe12.5)') 'Total mass: ',   &
     & qvtot+qctot+rnmstot+ticmstot+swmstot+grmstot+hlmstot
#endif

#if defined( ICE3) || defined(TAKON)
 11   format('processes: ',NUMPROCTMAX(1x,1pe12.5))
      write(iunit,11) (proctot(i), i=1,numproc)
!      write(iunit,'(a,2(1x,1pe12.5))') 'Process16: ',proctot(16),proctot(17)
#endif

       IF ( iraintypes > 0 ) THEN
        write(iunit,'(a,5(1x,1pe12.5))') 'FD source rates(A/S/M): ', (rate2dtot(i), i=1,nrate2d)
       ENDIF

      
      ENDIF ! my_rank == 0
      
! write layer max/mins to .stat file
      IF ( lstt ) THEN

#ifdef MPI
#ifdef CHGELEC
      IF ( ipelec .ge. 1 ) THEN
      allocate( mpitotinth((nzend),10))
      allocate(mpitotoutth((nzend),10))
      
      mpitotinth(:,:) = 0.0d0
      
      DO k = kzbeg,kzend

       mpitotinth(k, 1) =  cgizmx(k)
       mpitotinth(k, 2) = -cgizmn(k)
       mpitotinth(k, 3) =  cgszmx(k)
       mpitotinth(k, 4) = -cgszmn(k)
       mpitotinth(k, 5) =  chizmx(k)
       mpitotinth(k, 6) = -chizmn(k)
       mpitotinth(k, 7) =  chszmx(k)
       mpitotinth(k, 8) = -chszmn(k)
       mpitotinth(k, 9) =  cghwzmx(k)
       mpitotinth(k,10) = -cghwzmn(k)
       
      ENDDO

      n = 10*(nzend)

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
       DO k = 1,nzend
       cgizmx(k)  =  mpitotoutth(k, 1) 
       cgizmn(k)  = -mpitotoutth(k, 2) 
       cgszmx(k)  =  mpitotoutth(k, 3) 
       cgszmn(k)  = -mpitotoutth(k, 4) 
       chizmx(k)  =  mpitotoutth(k, 5) 
       chizmn(k)  = -mpitotoutth(k, 6) 
       chszmx(k)  =  mpitotoutth(k, 7) 
       chszmn(k)  = -mpitotoutth(k, 8) 
       cghwzmx(k) =  mpitotoutth(k, 9) 
       cghwzmn(k) = -mpitotoutth(k,10) 
       ENDDO
      ENDIF
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )

      allocate( mpitotinth((nzend),6))
      allocate(mpitotoutth((nzend),6))
      
      mpitotinth(:,:) = 0.0d0
      
      DO k = kzbeg,kzend
       mpitotinth(k, 1) =  ctghsnz(k)
       mpitotinth(k, 2) =  ctghspz(k)
       mpitotinth(k, 3) =  ctghinz(k)
       mpitotinth(k, 4) =  ctghipz(k)
       mpitotinth(k, 5) =  ctghwnz(k)
       mpitotinth(k, 6) =  ctghwpz(k)
      ENDDO

       n = 6*(nzend)

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
       ctghsnz(1:nzend)  =  mpitotoutth(1:nzend, 1) 
       ctghspz(1:nzend)  =  mpitotoutth(1:nzend, 2) 
       ctghinz(1:nzend)  =  mpitotoutth(1:nzend, 3) 
       ctghipz(1:nzend)  =  mpitotoutth(1:nzend, 4) 
       ctghwnz(1:nzend)  =  mpitotoutth(1:nzend, 5) 
       ctghwpz(1:nzend)  =  mpitotoutth(1:nzend, 6) 
      ENDIF
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )
       
       ENDIF ! ipelec .ge. 1
! endif for chgelec
#endif

       k0 = 15
       k = k0+numproc

       allocate( mpitotinth((nzend),k))
       allocate(mpitotoutth((nzend),k))
       
      mpitotinth(:,:) = 0.0d0
      
      DO kz = kzbeg,kzend
       mpitotinth(kz, 1) =  grvol(kz)
       mpitotinth(kz, 2) =  hlvol(kz)
       mpitotinth(kz, 3) =  grms(kz)
       mpitotinth(kz, 4) =  hlms(kz)
       mpitotinth(kz, 5) =  ticms(kz)
       mpitotinth(kz, 6) =  rnms(kz)
       mpitotinth(kz, 7) =  swms(kz)
       mpitotinth(kz, 8) =  grvol1(kz)
       mpitotinth(kz, 9) =  hlvol1(kz)
       mpitotinth(kz,10) =  grv(kz)
       mpitotinth(kz,11) =  hlv(kz)
       mpitotinth(kz,12) =  iwcmx(kz)
       mpitotinth(kz,13) =  grpotc(kz)
       mpitotinth(kz,14) =  grpotp(kz)
       IF ( lf > 1 .or. lnhf > 1 ) THEN
       mpitotinth(kz, 15) =  fdms(kz)
       ENDIF
!       mpitotinth(kz,13) =  liqms(kz)
      ENDDO
       
       DO i = 1,numproc
        DO kz = kzbeg,kzend
         mpitotinth(kz,k0+i) =  thproctot(kz,i)
        ENDDO
       ENDDO

!      write(0,*) 'mpitotinth k-1,k: ', my_rank, mpitotinth(1,k-1),mpitotinth(1,k)

      ! first sum columns 1-11
      n = (k0)*(nzend)

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

      ! column 12 is a max function -- overwrites outth
      n = nzend
      
      CALL MPI_Reduce(mpitotinth(1,12), mpitotoutth(1,12), n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      ! then sum the rest of the columns
      n = (numproc)*(nzend)

      CALL MPI_Reduce(mpitotinth(1,k0+1), mpitotoutth(1,k0+1), n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
       grvol(1:nzend)  =  mpitotoutth(1:nzend, 1) 
       hlvol(1:nzend)  =  mpitotoutth(1:nzend, 2) 
       grms(1:nzend)   =  mpitotoutth(1:nzend, 3) 
       hlms(1:nzend)   =  mpitotoutth(1:nzend, 4) 
       ticms(1:nzend)  =  mpitotoutth(1:nzend, 5) 
       rnms(1:nzend)   =  mpitotoutth(1:nzend, 6) 
       swms(1:nzend)   =  mpitotoutth(1:nzend, 7) 
       grvol1(1:nzend) =  mpitotoutth(1:nzend, 8) 
       hlvol1(1:nzend) =  mpitotoutth(1:nzend, 9) 
       grv(1:nzend)   =   mpitotoutth(1:nzend,10) 
       hlv(1:nzend)   =   mpitotoutth(1:nzend,11) 
       iwcmx(1:nzend)  =  mpitotoutth(1:nzend,12) 
       grpotc(1:nzend) = mpitotoutth(1:nzend,13)
       grpotp(1:nzend) = mpitotoutth(1:nzend,14)
       IF ( lf > 1 .or. lnhf > 1 ) fdms(1:nzend)   =  mpitotoutth(1:nzend,15) 
!       liqms(1:nzend)  =  mpitotoutth(1:nzend,13) 
       DO i = 1,numproc
         thproctot(1:nzend,i) =  mpitotoutth(1:nzend,k-numproc+i)
       ENDDO
!      write(0,*) 'mpitotoutth k-1,k: ', my_rank, mpitotoutth(1,k-1),mpitotoutth(1,k)
      ENDIF

       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )


#endif

#ifdef CHGELEC
      IF ( my_rank == 0 ) THEN
       write(3,'(a,i6)') 'Layer charging rates at time = ',time
      ENDIF

      IF ( ipelec .ge. 1 ) THEN

      IF ( my_rank == 0 ) THEN
      write(3,'(2a)') 'kz,cgizmx,cgizmn,cgszmx,cgszmn,',  &
     & 'chizmx,chizmn,chszmx,chszmn,cghwzmx,cghwzmn'

#ifdef MPI
      kzb = 1
      kze = nzend-1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg+1-kd1

      DO kz=kze,kzb,-1
#else
      DO kz=nz-kd1,1,-1
#endif
      write(3,'(1x,i3,10(2x,1pe13.5))' )  &
     & kz,cgizmx(kz),cgizmn(kz),cgszmx(kz),cgszmn(kz),  &
     & chizmx(kz),chizmn(kz),chszmx(kz),chszmn(kz),  &
     & cghwzmx(kz),cghwzmn(kz)
      END DO

       write(3,'(2a)') 'kz,ctghsnz,ctghspz,ctghinz,',  &
     & 'ctghipz,ctghwnz,ctghwpz,ctswinz,ctswipz'

#ifdef MPI
      kzb = 1
      kze = nzend-1
!       if (kzend .eq. nzend) kzb = kzend-kzbeg+1-kd1

       DO kz=kze,kzb,-1
#else
       DO kz=nz-kd1,1,-1
#endif
        write(3,'(1x,i3,8(2x,1pe13.5))' )  &
     & kz,  &
     &   ctghsnz(kz)*dtp*60./(tstat),  &
     &   ctghspz(kz)*dtp*60./(tstat),  &
     &   ctghinz(kz)*dtp*60./(tstat),  &
     &   ctghipz(kz)*dtp*60./(tstat),  &
     &   ctghwnz(kz)*dtp*60./(tstat),  &
     &   ctghwpz(kz)*dtp*60./(tstat),  &
     &   ctswinz(kz)*dtp*60./(tstat),  &
     &   ctswipz(kz)*dtp*60./(tstat)
       END DO

       ENDIF ! my_rank == 0

#ifdef MPI
      kzb = 1
      kze = nzend-1
!       if (kzend .eq. nzend) kze = kzend-kzbeg+1

       DO kz = kzb,kze
#else
       DO kz=1,nz-1
#endif
        ctswinz(kz) = 0.0
        ctswipz(kz) = 0.0
        ctghsnz(kz) = 0.0
        ctghspz(kz) = 0.0
        ctghinz(kz) = 0.0
        ctghipz(kz) = 0.0
        ctghwnz(kz) = 0.0
        ctghwpz(kz) = 0.0
        ctgsnz(kz) = 0.0
        ctgspz(kz) = 0.0
        ctginz(kz) = 0.0
        ctgipz(kz) = 0.0
        cthsnz(kz) = 0.0
        cthspz(kz) = 0.0
        cthinz(kz) = 0.0
        cthipz(kz) = 0.0
       ENDDO


      ENDIF ! ipelec .ge. 1
! endif for chgelec
#endif


      IF ( my_rank == 0 ) THEN

      write(3,'(a,i6)') 'Layer volumes/masses/rates at time = ',time
      write(3,'(a)')   &
     &      'kz,udmf,uicmf,grvol,hlvol,grms,hlms,ticms,rnms,snms,iwcmx,grvol1,hlvol1,grv,hlv,udke,ddke,grpotc,grpotp,fdms'
      kzb = 1
      kze = nzend-1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg+1-kd1

      DO kz=kze,kzb,-1
      tmp = 0.0
      IF ( grv(kz) .gt. 1.e-20 ) tmp = grms(kz)/grv(kz)
      write(3,'(1x,i3,20(2x,1pe13.5))' )  &
     & kz, udmf(kz), uicmf(kz), grvol(kz), hlvol(kz),  &
     & grms(kz), hlms(kz), ticms(kz), rnms(kz),swms(kz),iwcmx(kz),grvol1(kz),hlvol1(kz),  &
     & grv(kz),hlv(kz),tmp,udke(kz),ddke(kz),grpotc(kz),grpotp(kz),fdms(kz)
      END DO


      write(3,'(a,i6)') 'Layer microphysics processes at time = ',time,' nzend = ',nzend
      write(3,'(a)')   &
     &      'kz,crfrz,ciacrf,chcnsh,chcnih'
      
       DO kz=kze,kzb,-1
         write(3,111 ) kz,0.001*z1d4(kz,1),temc(kz),(thproctot(kz,i), i=1,numproc)
  111    format (1x,i3,1x,f6.2,1x,f6.2,NUMPROCTMAX(2x,1pe13.5))
       ENDDO

       
      ENDIF  ! my_rank == 0
      
      thproctot(:,:) = 0.0
      
      ENDIF  ! lstt


! save mixing ratio and space charge in compact gather arrays      
!#ifdef ICE3
!      IF ( na .ge. 8 ) THEN
!      CALL mic2savlfo(nx,ny,nz,nor,na,
!     :  an,rstime,irstime,
!     :  istag,jstag,kstag,ipelec )
!      ENDIF
!#else
!      CALL mic2sav(nx,ny,nz,nor,na,
!     :  an,rstime,irstime,
!     :  istag,jstag,kstag,ipelec,ipconc )
!#endif
!      
!      END IF ! lhelec
!

#ifdef CHGELEC
!
!  DOIN' THE ELECTRICITY THING!!!
!
!

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 6"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 6"
#endif

      if ( ipelec .ge. 2 ) then !{

      CALL cld_cpu('MICROPHYSICS')
      CALL cld_cpu('ELECTRICITY')

       IF ( .not. lpredict ) THEN

#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          jyb = 1
          jye = jtile
!          if (jybeg .le. jbsd) jyb = jbsd-jybeg+1
!          if (jyend .ge. jesd) jye = jesd-jybeg+1

          ixb = 1
          ixe = itile
!          if (ixbeg .le. ibsd) ixb = ibsd-ixbeg+1
!          if (ixend .ge. iesd) ixe = iesd-ixbeg+1

      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

          DO kz = kzb,kze
           DO jy = jyb,jye
            DO ix = ixb,ixe
#else
          DO kz = 1,nzend-kzbeg
           DO jy = jbsd,jesd ! 1,ny-1
            DO ix = ibsd,iesd ! 1,nx-1
#endif
              lgtth(kzbeg-1+kz,chans) = lgtth(kzbeg-1+kz,chans) + Abs(elec(ieflshn)%flt3d(ix,jy,kz))  &
     &                           + Abs(elec(ieflshp)%flt3d(ix,jy,kz))
              lgtth(kzbeg-1+kz,chansp) = lgtth(kzbeg-1+kz,chansp) + elec(ieflshp)%flt3d(ix,jy,kz)
              lgtth(kzbeg-1+kz,chansn) = lgtth(kzbeg-1+kz,chansn) + elec(ieflshn)%flt3d(ix,jy,kz)
              lgtth(kzbeg-1+kz,cinit) =  lgtth(kzbeg-1+kz,cinit) + elec(ieinit)%flt3d(ix,jy,kz)
            ENDDO
           ENDDO
          ENDDO


         GOTO 998  ! skip to deallocate statements at end of routine
       ENDIF 

      CALL ELEC_DRIVE_LTGLOOP (gd,                                    &
     &                     pinit, ab,db,                              & ! base state
     &                         p2, tt0,                               &! pi'
     &                         an,                                    & ! scalars
     &                         xfalltot,                              & ! precip
     &                        dtp1, gx,gy,gz,z1d,                     & ! time step, vertical grid spacing
     &                        nx, ny, nz, nor, na,                    & 
     &                        u,v,w,                                  & ! vn 
     &                        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,          & ! temporary arrays
     &                        dx,dy,dz,                               & ! nstep, start,tstat,
     &                        ntmul,iunit,istag,jstag,kstag,          &
     &                        dbz, io_flag, time, time_real, dt,  &
     &                        elec, cion, muz, gxt, gyt, gzt,         & ! electricity arrays
     &                        tstat, lstt, bcx, bcy, bcx1, bcy1, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,    &
     &                        xfall,kcldtop,icldtop,zcldtop,          &
     &                        idelay,idelaycg)



!#ifndef SWM


#ifndef BOXMG
      IF ( ion .eq. 1 .and. number_of_processes .eq. 1 ) THEN
#else
      IF ( ion .eq. 1 ) THEN !{
#endif

      CALL cld_cpu('IONS')
#ifdef BOXMG

         nbw = nbwmg(1)
         nbe = nbemg(1)
         nbs = nbsmg(1)
         nbn = nbnmg(1)

#else

         nbw = (nnxs(1) - nx - 1)/2
         nbe = nnxs(1) - nx - 1 - nbw
         nbs = (nnys(1) - ny - 1)/2
         nbn = nnys(1) - ny - 1 - nbs
#endif

         IF ( ny .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         IF ( bcx .eq. 2 ) THEN ; nbw = -1 ; nbe = 0 ; ENDIF
         IF ( bcy .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

      IF ( .true. ) THEN
      call ionstep  &
     & (nx,ny,nz,na,nba,nor,nstep,  &
     &  dtp1,dx,dy,dz,gxt,gyt,gzt,dxx,dyy,dzz,  &
     &  an,dn,pn,  &
     &  db,pb,  &
     &  elec,uz,cioncp,  &
     &  cwccn,ipconc,cimn,cimx, &
     &  t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,  &
     &  tt0,tt7,  &
     &   nnxs,nnys,nnzs,nztop,  &
     &   iixps,jjyqs,kkzrs,iiexs,jjeys,kkezs,  &
     &   nbw,nbe,nbs,nbn,  &
     &   llworks,istretch,  &
     &   ibg,iunit,  &
     &   id1,jd1,kd1,istag,jstag,kstag,  &
     &   microp,iestag,bcx,bcy,  &
     &   gd%xyz3d%flt4d(-nor+1,-nor+1,-nor+1,elec(iex)%index))
 
      ENDIF

      IF ( ndebug .ge. 2 ) THEN
       write(iunit,*) 'post ion values'

!#ifdef MPI
      kzb = ktile
      kze = 1
      if (kzend .eq. nzend) kzb = kzend-kzbeg+1

      do kz = kzb,kze,-1
!#else
!      DO kz=nz,1,-1
!#endif
       write(iunit,'(a,i3,5(2x,1pe12.5))')   &
     &  'kz,test: ',kz,gzt(kz,1),  &
     &   an(nx/2,ny/2,kz,lscpi)-cioncp(kz,1), cioncp(kz,1),  &
     &   an(nx/2,ny/2,kz,lscni)-cioncp(kz,2), cioncp(kz,2)
      ENDDO
      ENDIF

      CALL cld_cpu('IONS')

      ENDIF !}
!#endif
      CALL cld_cpu('ELECTRICITY')
      CALL cld_cpu('MICROPHYSICS')
      end if ! } ipelec .gt. 2

      IF ( ipelec .ge. 1 ) THEN

      CALL cld_cpu('MICROPHYSICS')
      CALL cld_cpu('ELECTRICITY')

      CALL ELEC_DRIVE4     (gd,                                     &
     &                     pinit, ab,                               &   ! base state
     &                         p2,                                  &  ! pi'
     &                         an,                                  &   ! scalars
     &                         xfalltot,                            &   ! precip
     &                        dtp1, gx,gy,gz,z1d,                   &   ! time step, vertical grid spacing
     &                        nx, ny, nz, nor, na,                  & 
     &                        u,v,w,                                &   ! vn 
     &                        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,        &  ! temporary arrays
     &                        dx,dy,dz,                             &  ! nstep, start,tstat,
     &                        ntmul,iunit,istag,jstag,kstag,        & 
     &                        dbz, io_flag, time, time_real, dt, &
     &                        elec, cion, muz, gxt, gyt, gzt,        & ! electricity arrays
     &                        tstat, lstt, bcx, bcy, tstop,        &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,   &
     &                        kcldtop)

!#ifndef SWM
      

      CALL cld_cpu('ELECTRICITY')
      CALL cld_cpu('MICROPHYSICS')
      
      ENDIF ! ipelec .ge. 1

! endif for CHGELEC
#endif

!#endif
! cut for SWM

!#endif
! cut for SWM


! #ifndef SWM
#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 9"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 9"
#endif

!
! Vis5D output
!


       IF ( lv5dwrite .or. lncwrite .or. lrst .or. lprt .or. lstt .or. io_flag ) THEN
!         IF ( .not. lstt ) THEN
#ifndef LFO3
#ifndef TAKON
#ifndef HCMON
         
#ifdef CM1
         call radardd02(nx,ny,nz,nor,na,an,tt0,  &
     &    dbz,dn, nzend ,cnoh0,rho_qh, ipconc, iunit, microp, 1, 0)
#else
         IF ( associated( vzf%flt3d ) ) THEN
         call radardd02(nx,ny,nz,nor,na,an,tt0,  &
     &    dbz,dn, nzend ,cnoh0,rho_qh, ipconc, iunit, microp, 1, 1, vzf%flt3d)

           IF ( lwvzf > 0 ) THEN
             DO kz = kzb,kze
              DO jy = jyb,jye
               DO ix = ixb,ixe
                 axtra(ix,jy,kz,lwvzf) = 0.5*(w(ix,jy,kz)+w(ix,jy,kz+kstag)) - vzf%flt3d(ix,jy,kz)
               ENDDO
              ENDDO
             ENDDO
           ENDIF
         
         ELSE

         call radardd02(nx,ny,nz,nor,na,an,tt0,  &
     &    dbz,dn, nzend ,cnoh0,rho_qh, ipconc, iunit, microp, 1, 0, t0) ! t0 being passed but should not be used

         ENDIF
#endif

#endif
#endif
#endif
         
!         ENDIF

!#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!          jyb = 1
!          jye = jtile
!          if (jyend .eq. nyend) jye = jyend-jybeg+1
!
!          ixb = 1
!          ixe = itile
!          if (ixend .eq. nxend) ixe = ixend-ixbeg+1
!
!          DO k=kzb,kze
!            DO j=jyb,jye
!              DO i=ixb,ixe
!#else
!          DO k=1,nz
!            DO j=1,ny
!              DO i=1,nx
!#endif
!                dbz(i,j,k) = tt7(i,j,k)
!              ENDDO
!            ENDDO
!          ENDDO
         
        ENDIF

 998  CONTINUE

#ifndef CM1

      IF ( lstt ) THEN
!
! Print out lightning time-height info
!

#ifdef CHGELEC

      IF ( ipelec .ge. 1 ) THEN

#ifdef MPI

!        allocate( lgtth(nz,14) )

      allocate(mpitotoutthr((nzend), 14))
      

      n = 14*(nzend)

      CALL MPI_Reduce(lgtth, mpitotoutthr, n, MPI_REAL, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        DO k = 1,14
          lgtth(1:nzend,k) = mpitotoutthr(1:nzend,k)
        ENDDO
      ENDIF
       
       deallocate( mpitotoutthr )

#endif

      IF ( my_rank == 0 ) THEN
      write(3,'(a,i6)') 'Lightning sums for time = ',time
       write(3,'(3a)') ' Time Altitude Channels chanp chann ',  &
     &  'cgchan cgchanp cgchann icchan icchnp icchann ',  &
     &  'cinit cginit icinit cgninit cgpinit'


!#ifdef MPI
      kzb = nzend-kd1
      kze = 1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg+1-kd1

      DO kz=kzb,kze,-1
!#else
!      DO kz=nz-kd1,1,-1
!#endif
          write(3,'(1x,f8.1,1x,f6.2,14(1x,f7.0))')   &
     &    float(time)/60.0,0.001*z1d4(kz,1), (lgtth(kz,i), i=1,14)
!     :     chans(kz,it), chansp(kz,it), chansn(kz,it), 
!     :     cgchans(kz,it), cgchansp(kz,it), cgchansn(kz,it), 
!     :     icchans(kz,it), icchansp(kz,it), icchansn(kz,it), 
!     :     cinit(kz,it), 
!     :     cginit(kz,it), Float(icinit(kz,it)),
!     :     cgninit(kz,it),cgpinit(kz,it)
      END DO
      
      ENDIF

      lgtth(:,:) = 0.0
      
      ENDIF
! endif for CHGELEC
#endif

!
! print out time-hieght maximum reflectivity
!

      IF ( my_rank == 0 ) THEN
      write(3,'(a,i6)') 'Max dBZ and integrated UDV at time = ',time
      write(3,'(a,a)') 'Time, altitude, dBZ-max, dBZ-I ',  &
     &        'udvm5 udv5 udv10 udv20'
      ENDIF

      allocate(mpitotinthr((nzend), 6))
      allocate(mpitotoutthr((nzend), 6))

      mpitotinthr(:,:) = 0.0
      mpitotinthr(1:nzend,1) = -10
      
!#ifdef MPI
         kzb = 1
         kze = ktile
         if (kzend .eq. nzend) kze = kzend-kzbeg

         DO kz = kzb,kze
!#else
!         DO kz = 1,nzend
!#endif
           dbzmx = -100.0
           dbzi = 0.0
           udv5n = 0.
           udv5 = 0.
           udv10 = 0.
           udv20 = 0.
           
!           dv = dx*dy*dz/gt(1,1,kz,imapz)

#ifdef MPI
           jyb = 1
           jye = jtile
!           if (jybeg .le. jbsd) jyb = jbsd-jybeg+1
!           if (jyend .ge. jesd) jye = jesd-jybeg+1

           ixb = 1
           ixe = itile
!           if (ixbeg .le. ibsd) ixb = ibsd-ixbeg+1
!           if (ixend .ge. iesd) ixe = iesd-ixbeg+1

      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag
           
           DO jy = jyb,jye
            DO ix = ixb,ixe
#else
!           DO jy= 1,ny-1
!            DO ix = 1,nx-1
           DO jy = jbsd,jesd
            DO ix = ibsd,iesd
#endif
              dxdy = dxx(ix)*dyy(jy)
              dv = dxdy*dzz(kz)
!              dbzi = dbzi + Max( 0.0, tt7(ix,jy,kz) )
!              dbzmx = Max( dbzmx, tt7(ix,jy,kz) )
              mpitotinthr(kzbeg-1+kz,2) = mpitotinthr(kzbeg-1+kz,2)  + Max( 0.0, dbz(ix,jy,kz) )
!              dbzi = dbzi + Max( 0.0, dbz(ix,jy,kz) )
              mpitotinthr(kzbeg-1+kz,1) = Max( mpitotinthr(kzbeg-1+kz,1), dbz(ix,jy,kz) )
!              dbzmx = Max( dbzmx, dbz(ix,jy,kz) )
              IF ( w(ix,jy,kz) .le. -5.0 ) THEN
               mpitotinthr(kzbeg-1+kz,3) = mpitotinthr(kzbeg-1+kz,3) + dv
!               udv5n = udv5n + dv
              ENDIF
              IF ( w(ix,jy,kz) .ge.  5.0 ) THEN
                mpitotinthr(kzbeg-1+kz,4) = mpitotinthr(kzbeg-1+kz,4) + dv
!                udv5 = udv5 + dv
              ENDIF

              IF ( w(ix,jy,kz) .ge. 10.0 ) THEN
                mpitotinthr(kzbeg-1+kz,5) = mpitotinthr(kzbeg-1+kz,5) + dv
!                udv10 = udv10 + dv
              ENDIF

              IF ( w(ix,jy,kz) .ge. 20.0 ) THEN
                mpitotinthr(kzbeg-1+kz,6) = mpitotinthr(kzbeg-1+kz,6) + dv
!                udv20 = udv20 + dv
              ENDIF
            ENDDO
           ENDDO
         ENDDO

#ifdef MPI
      n = 5*(nzend)

      CALL MPI_Reduce(mpitotinthr(1,2), mpitotoutthr(1,2), n, MPI_REAL, MPI_SUM, 0, my_comm, mpi_error_code)

      n = (nzend)

      CALL MPI_Reduce(mpitotinthr(1,1), mpitotoutthr(1,1), n, MPI_REAL, MPI_MAX, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        DO k = 1,6
          mpitotinthr(1:nzend,k) = mpitotoutthr(1:nzend,k)
        ENDDO
      ENDIF

#endif

        IF ( my_rank == 0 ) THEN
!#ifdef MPI
         kzb = 1
         kze = nzend
!         if (kzend .eq. nzend) kze = kzend-kzbeg

         DO kz = kzb,kze
!#else
!         DO kz = 1,nzend
!#endif
           write(3,'(1x,f8.1,1x,f6.2,2(1x,f9.2),4(1x,1pe13.5))')   &
     &     float(time)/60.0,0.001*z1d4(kz,1),mpitotinthr(kz,1),mpitotinthr(kz,2),  &
     &      mpitotinthr(kz,3),mpitotinthr(kz,4),mpitotinthr(kz,5),mpitotinthr(kz,6)
!     :     float(time)/60.0,0.001*gzt(kz,1),dbzmx,dbzi,
!     :      udv5n,udv5,udv10,udv20
         ENDDO
        ENDIF


      deallocate(mpitotinthr)
      deallocate(mpitotoutthr)

      ENDIF  ! lstt

      IF ( ifv5d ) THEN
        IF ( lv5dwrite ) THEN
      
        nr=ny
        nr0 = Max(2,nr-1)
        nr1 = Min(2,ny-1)
        nr2 = Max(1,ny-2)
        nc=nx
!        varnum = 1
        nl(1)=nz-1
      
          times = (nstep)*(dtp)
        idatime = time ! Int(times)
        hh = idatime/3600
        idatime = idatime - 3600*hh
        mm = idatime/60
        ss = idatime - 60*mm
     
        it = itv5d
#ifdef MPI
!          write(0,*) 'DRIVER: REMEMBER TO TURN V5D FILE BACK ON IN DRIVER AFTER IMPLEMENTATION'
#else
          write(iunit,'(a,a)') 'Writing Vis5D file: ',v5dfilename
#endif
!        itimes(it+1) = 10000*hh+100*mm+ss
!        write(0,*) 'itimes(it+1) = ',it+1, itimes(it+1)
!        idates(it+1) = 010101
        
        IF ( it .gt. 1 ) THEN
        itimes(it) = 10000*hh+100*mm+ss
        write(iunit,*) 'itimes(it) = ',it, itimes(it)
        idates(it) = 2006001  ! YYYYDDD
        ENDIF

        
!        write(0,*) 'DRIVER: varnum = ',varnum
        
        DO iv = 1, numvars
        

        IF ( iv5dwritten(iv) .eq. 0 ) THEN
          iv5dwritten(iv) = 1
        ELSE
          CYCLE
        ENDIF

         fac = 1./Float(vis5dstridex*vis5dstridex*vis5dstridez)

         tem1(:,:,:) = 0.0

!        write(0,*) 'DRIVER: Writing V5D variable ',iv,', ',varname(iv)
        dunit = '                   '
        write(dunit,'(a)') units(idx(iv))
! #ifndef MPI
        call v5dsetunits(iv, dunit )
! #endif

        
        IF ( idx(iv) .eq. 410 .and. ipelec .gt. 0 ) THEN ! scnet

!#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg
!
!          ixb = 1
!          ixe = itile
!          if (ixend .ge. nc) ixe = nc-ixbeg+1
!
!          jyb = 1
!          jye = jtile
!          if (jyend .ge. nr) jye = nr-jybeg+1
!
!          DO k=kzb,kze
!            DO i=ixb,ixe
!              DO j=jyb,jye
!#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
!#endif

#ifdef MPI
                tem1(j,i,k) = 1.0e9*t0(i,nr-j,k)  !mpidebug: need to fix this! 
#else
                tem1(j,i,k) = 1.0e9*t0(i,nr-j,k)
#endif
              ENDDO
            ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 421  .and. ipconc .eq. 0 ) THEN ! cci

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
               IF ( lni .gt. 1 ) THEN
                tem1(j,i,k) =  1.0e-3*an(i,nr-j,k,lni)   ! per liter
               ELSE
                tem1(j,i,k) =   1.0e-3*tt7(i,nr-j,k)   ! per liter
               ENDIF
            ENDDO
           ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 432 .and. ipconc .eq. 0 ) THEN

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
               IF ( lni .gt. 1 ) THEN
                tem1(j,i,k) =  Min( ccimx, 1.0e-3*an(i,nr-j,k,lni))  ! per liter
               ELSE
                tem1(j,i,k) =  Min( ccimx, 1.0e-3*tt7(i,nr-j,k))  ! per liter
               ENDIF
            ENDDO
           ENDDO
          ENDDO
         
         
         ELSEIF ( idx(iv) .eq. 434 .and. ipconc .ge. 2 ) THEN ! SS

#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                IF ( lsat > 1 ) THEN
                tem1(j,i,k) = axtra(i,nr-j,k,lsat)
                ELSE
                tem1(j,i,k) = ssat(i,nr-j,k)
                ENDIF
#else
                IF ( lsat > 1 ) THEN
                tem1(j,i,k) = axtra(i,nr-j,k,lsat)
                ELSE
                tem1(j,i,k) = ssat(i,nr-j,k)
                ENDIF
#endif
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 435 .and. ipconc .ge. 2 ) THEN ! SSI
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                IF ( lsati > 1 ) THEN
                tem1(j,i,k) = axtra(i,nr-j,k,lsati)
                ELSE
                tem1(j,i,k) = ssati(i,nr-j,k)
                ENDIF
#else
                IF ( lsati > 1 ) THEN
                tem1(j,i,k) = axtra(i,nr-j,k,lsati)
                ELSE
                tem1(j,i,k) = ssati(i,nr-j,k)
                ENDIF
#endif
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 437 .and. ipconc .ge. 2 ) THEN ! ssmax
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                tem1(j,i,k) = an(i,nr-j,k,lss)  !mpidebug: need to fix this!
#else
                tem1(j,i,k) = an(i,nr-j,k,lss)
#endif
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 439 .and. ipconc .ge. 2 ) THEN ! SSNEW
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                tem1(j,i,k) = ssfilt(i,nr-j,k)  !mpidebug: need to fix this!
#else
                tem1(j,i,k) = ssfilt(i,nr-j,k)
#endif
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 464 .and. ipconc .eq. 0 ) THEN ! cidia (columns)

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( lni .gt. 1 ) THEN
               IF ( an(i,nr-j,k,li) .gt. 1.e-12 .and. an(i,nr-j,k,lni) .gt. 1.e-1 ) THEN
                 tmp = an(i,nr-j,k,li)*dn(i,nr-j,k)/an(i,nr-j,k,lni)  ! mass per crystal
                 tem1(j,i,k) =  1.e6*0.1871*(tmp**(0.3429))  ! factor of 1.e6 to convert to microns
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                IF ( an(i,nr-j,k,li) .gt. 1.e-12 .and. tt7(i,nr-j,k) .gt. 1.e-1 ) THEN
                  tmp = an(i,nr-j,k,li)*dn(i,nr-j,k)/tt7(i,nr-j,k)  ! mass per crystal
                  tem1(j,i,k) =  1.e6*0.1871*(tmp**(0.3429))  ! factor of 1.e6 to convert to microns
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 790  ) THEN ! crfrz
         tmp = 0.0
!#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg
!
!          ixb = 1
!          ixe = itile
!          if (ixend .ge. nc) ixe = nc-ixbeg+1
!
!          jyb = 1
!          jye = jtile
!          if (jyend .ge. nr) jye = nr-jybeg+1
!
!          DO k=kzb,kze
!            DO i=ixb,ixe
!              DO j=jyb,jye
!#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
!#endif
                tmp = tmp + tq1(i,nr-j,k,1)
#ifdef MPI
                tem1(j,i,k) = tq1(i,nr-j,k,1)  !mpidebug: need to fix this!
#else
                tem1(j,i,k) = tq1(i,nr-j,k,1)
#endif

            ENDDO
           ENDDO
          ENDDO
          
          write(iunit,*) 'crfrz tot = ',tmp

         ELSEIF ( idx(iv) .eq. 791  ) THEN ! chcnsi
         
         tmp = 0.0
!#ifdef MPI
!          kzb = 1
!          kze = ktile
!          if (kzend .eq. nzend) kze = kzend-kzbeg
!
!          ixb = 1
!          ixe = itile
!          if (ixend .ge. nxend) ixe = nc-1
!
!          jyb = 1
!          jye = jtile
!          if (jyend .ge. nr) jye = nr-jybeg+1
!
!          DO k=kzb,kze
!            DO i=ixb,ixe
!              DO j=jyb,jye
!#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
!#endif
                tmp = tmp + tq1(i,nr-j,k,2)

#ifdef MPI
                tem1(j,i,k) = tq1(i,nr-j,k,2)  !mpidebug: need to fix this!
#else
                tem1(j,i,k) = tq1(i,nr-j,k,2)
#endif


            ENDDO
           ENDDO
          ENDDO

          write(iunit,*) 'chcnsi tot = ',tmp

         ELSEIF ( idx(iv) .eq. 800  ) THEN ! qproc
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                tem1(j,i,k) = 1.e12*tq1(i,nr-j,k,1)
#else
                tem1(j,i,k) = 1.e12*tq1(i,nr-j,k,1)
#endif

            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 860  ) THEN ! qhshr
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
                tem1(j,i,k) = 1.e6*tem2(i,nr-j,k)
#else
                tem1(j,i,k) = 1.e6*tem2(i,nr-j,k)
#endif

            ENDDO
           ENDDO
          ENDDO



#ifdef ICE10
         ELSEIF ( idx(iv) .eq. 481 .and. ipconc .ge. 5 ) THEN ! N0GL
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
              IF ( an(i,nr-j,k,lgl) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngl) .gt. 1.e-3 ) THEN
                tmp = (300.*an(i,nr-j,k,lngl))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgl))
                tmpg = an(i,nr-j,k,lngl)*(tmp*(3.14159))**(1./3.)  !mpidebug: need to fix this!
#else
              IF ( an(i,nr-j,k,lgl) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngl) .gt. 1.e-3 ) THEN
                tmp = (300.*an(i,nr-j,k,lngl))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgl))
                tmpg = an(i,nr-j,k,lngl)*(tmp*(3.14159))**(1./3.)
#endif
                IF ( tmpg .gt. 1. ) THEN
                  tem1(j,i,k) = Log10( tmpg )
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 482 .and. ipconc .ge. 5 ) THEN ! N0gm
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
              IF ( an(i,nr-j,k,lgm) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngm) .gt. 1.e-3 ) THEN
                tmp = (500.*an(i,nr-j,k,lngm))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgm))
                tmpg = an(i,nr-j,k,lngm)*(tmp*(3.14159))**(1./3.)  !mpidebug: need to fix this!
#else
              IF ( an(i,nr-j,k,lgm) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngm) .gt. 1.e-3 ) THEN
                tmp = (500.*an(i,nr-j,k,lngm))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgm))
                tmpg = an(i,nr-j,k,lngm)*(tmp*(3.14159))**(1./3.)
#endif
                IF ( tmpg .gt. 1. ) THEN
                  tem1(j,i,k) = Log10( tmpg )
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 483 .and. ipconc .ge. 5 ) THEN ! N0gh
#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          jyb = 1
          jye = jtile
          if (jyend .ge. nr) jye = nr-jybeg+1

          DO k=kzb,kze
            DO i=ixb,ixe
              DO j=jyb,jye
#else
          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
#endif

#ifdef MPI
              IF ( an(i,nr-j,k,lgh) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngh) .gt. 1.e-3 ) THEN
                tmp = (700.*an(i,nr-j,k,lngh))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgh))
                tmpg = an(i,nr-j,k,lngh)*(tmp*(3.14159))**(1./3.)  !mpidebug: need to fix this!
#else
              IF ( an(i,nr-j,k,lgh) .gt. 1.e-6 .and.   &
     &             an(i,nr-j,k,lngh) .gt. 1.e-3 ) THEN
                tmp = (700.*an(i,nr-j,k,lngh))/  &
     &                 (dn(i,nr-j,k)*an(i,nr-j,k,lgh))
                tmpg = an(i,nr-j,k,lngh)*(tmp*(3.14159))**(1./3.)
#endif
                IF ( tmpg .gt. 1. ) THEN
                  tem1(j,i,k) = Log10( tmpg )
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

#endif


         
         ELSE
           iv5dwritten(iv) = 0
           CYCLE
         ENDIF ! idx check
         
!
! Copy slab to j=2 for 2D simulation because vis5d requires a minium of 2 rows.
!
#ifdef MPI
         IF ( nyend .eq. 2 ) THEN
#else
         IF ( ny .eq. 2 ) THEN
#endif

#ifdef MPI
          kzb = 1
          kze = ktile
          if (kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if (ixend .ge. nc) ixe = nc-ixbeg+1

          DO k = kzb,kze
            DO i = ixb,ixe
#else
          DO k=1,nz-1
            DO i=1,nc-1
#endif
              tem1(2,i,k) =  tem1(1,i,k)
            ENDDO
          ENDDO
          ENDIF

           itimes0 = itimes(it)
           idates0 = idates(it)
           
!#ifndef MPI
         n = v5dwriteappend(it,iv,tem1,itimes0,idates0)

         IF ( vis5dstridex == 1 .and. vis5dstridez == 1 ) THEN
         n = v5dwriteappend(it,iv,tem1,itimes0,idates0)
         
         ELSE
          tem3(:,:,:) = 0.0
 !         fac = 1./Float(vis5dstridex*vis5dstridex*vis5dstridez)
          DO k=1,(nz-1)/vis5dstridez
             DO i=1,Max(1,(nc-1)/vis5dstridex)
               DO j=1,Max(1, (nr-1)/vis5dstridex )
               
               smax = 0.0
               smin = 0.0
               DO k1 = 1,vis5dstridez
               DO i1 = 1,Min(nx-1,vis5dstridex)
               DO j1 = 1,Min(ny-1,vis5dstridex)
                ix = (i-1)*vis5dstridex + i1
                jy = (j-1)*vis5dstridex + j1
                kz = (k-1)*vis5dstridez + k1
                smin = Min( smin, tem1(jy,ix,kz) )
                smax = Max( smax, tem1(jy,ix,kz) )
                tem3(j,i,k) = tem3(j,i,k) + tem1(jy,ix,kz)
               ENDDO
               ENDDO
               ENDDO
               tem3(j,i,k) = Max( smin, Min( smax, fac * tem3(j,i,k) ) )
             ENDDO
            ENDDO
           ENDDO
         
          IF ( ny .eq. 2 ) THEN
           DO k=1,nz-1
             DO i=1,(nc-1)/vis5dstridex
               tem3(2,i,k) =  tem3(1,i,k)
             ENDDO
           ENDDO
          ENDIF
 
         n = v5dwriteappend(it,iv,tem3,itimes0,idates0)
         
 !        write(0,*) 'writeappend: n,iv,idx(iv) = ',
 !     :      n,iv,idx(iv),itimes(it),idates(it)
         ENDIF
!#endif
        ENDDO ! iv

!#ifndef MPI
          n = myv5dupdate()
!#endif
!
! reset the lightning and initiation arrays
!        
        
        ENDIF
      

      deallocate ( tem1 )
      deallocate ( tem2 )
      IF ( allocated ( tem3 ) ) deallocate ( tem3 )
      
      ENDIF ! ifv5d

! endif for CM1
#endif

! #endif
! cut for SWM
!
      
 999  CONTINUE
 
      deallocate ( tt0 )
      deallocate ( tt7 )
      deallocate ( ssat )
      deallocate ( ssati )
      deallocate ( ssfilt )

#ifndef CM1
      deallocate ( dn )
      deallocate ( pn )
#endif
      deallocate ( xfall )
      deallocate ( rate2d )
      
      deallocate ( tq1 )

 
      CALL cld_cpu('MICROPHYSICS')

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: EXITING SUBROUTINE"
#endif

      RETURN
      END


#ifdef CM1
      subroutine cld_cpu(instring)
      implicit none
      character(*) instring
      return
      end
#endif

#if defined(SWM) & defined(ICE10)

#if defined(__ppc__)
! #####################################################################
      subroutine flush(i)
      implicit none
      integer i
       call flush_(i)
      RETURN
      END

#endif
!
! ###############################################################      
!  SUBROUTINE LDTEMPS
! ###############################################################      
      SUBROUTINE LDTEMPS (nxl,nyl,nzl,nor,ngs,ngscnt,nstep,nsave,  &
     &  igs,jy,kgs,  &
     &  scsaci,scsacip,scsacir,scglacs,  &
     &  scgmacs,scghacs,  &
     &  scfacs,schacs,schlacs, scglaci,   &
     &  scglacip, scglacir, scgmaci,  &
     &  scgmacip, scgmacir, scghacir,scghacip,  &
     &  scghaci, scfacir,  &
     &  scfacip, scfaci, schacir, schacip,   &
     &  schaci, schlacir,  &
     &  schlacip, schlaci, scsacw, scglacw,   &
     &  scgmacw, scghacw,  &
     &  scfacw, schacw, schlacw, dezcomp,   &
     &  vtxbar, &
!     :  vthbar, vthlbar,vtglbar,vtgmbar,vtghbar,vtfbar,  &
     &  sxsaci,sxsacip,sxsacir,sxglacs,sxglaci,sxglacip,sxglacir,  &
     &  sxgmacs,sxgmaci,sxgmacip,sxgmacir,sxghacs,sxghacir,sxghacip,  &
     &  sxghaci,sxfacs,sxfaci,sxfacir,sxfacip,sxhacs,sxhaci,sxhacir,  &
     &  sxhacip,sxhlacs,sxhlaci,sxhlacir,sxhlacip,  &
     &  t1,t2,t3,t4,t5,t6,t8,t9)
      
      
       USE INDEX_MODULE

      implicit none
!      include 'sam.index.ion.h'
      integer nxl,nyl,nzl,nor,ngs,ngscnt,mgs
      integer nstep,nsave
      integer jy,igs(ngs),kgs(ngs)
      real scsaci(ngs),scsacip(ngs),scsacir(ngs),scglacs(ngs),  &
     &  scgmacs(ngs),scghacs(ngs),  &
     &  scfacs(ngs),schacs(ngs),schlacs(ngs), scglaci(ngs),   &
     &  scglacip(ngs), scglacir(ngs), scgmaci(ngs),  &
     &  scgmacip(ngs), scgmacir(ngs), scghacir(ngs),scghacip(ngs),  &
     &  scghaci(ngs), scfacir(ngs),  &
     &  scfacip(ngs), scfaci(ngs), schacir(ngs), schacip(ngs),   &
     &  schaci(ngs), schlacir(ngs),  &
     &  schlacip(ngs), schlaci(ngs), scsacw(ngs), scglacw(ngs),   &
     &  scgmacw(ngs), scghacw(ngs),  &
     &  scfacw(ngs), schacw(ngs), schlacw(ngs), dezcomp(ngs)
!     :  vthbar(ngs), vthlbar(ngs),
!     :  vtglbar(ngs),vtgmbar(ngs),vtghbar(ngs),vtfbar(ngs)

      real vtxbar(ngscnt,lc:lhab,2)
      
      real sxsaci(ngs),  &
     & sxsacip(ngs),  &
     & sxsacir(ngs),  &
     & sxglacs(ngs),  &
     & sxglaci(ngs),  &
     & sxglacip(ngs),  &
     & sxglacir(ngs),  &
     & sxgmacs(ngs),  &
     & sxgmaci(ngs),  &
     & sxgmacip(ngs),  &
     & sxgmacir(ngs),  &
     & sxghacs(ngs),  &
     & sxghacir(ngs),  &
     & sxghacip(ngs),  &
     & sxghaci(ngs),  &
     & sxfacs(ngs),  &
     & sxfaci(ngs),  &
     & sxfacir(ngs),  &
     & sxfacip(ngs),  &
     & sxhacs(ngs),  &
     & sxhaci(ngs),  &
     & sxhacir(ngs),  &
     & sxhacip(ngs),  &
     & sxhlacs(ngs),  &
     & sxhlaci(ngs),  &
     & sxhlacir(ngs),  &
     & sxhlacip(ngs) 
!
! external temporary arrays
!
      real t1(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t2(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t3(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t4(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t5(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t6(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t8(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)
      real t9(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor)

! ###############################################################      

!DIR$ IVDEP
      do mgs = 1,ngscnt
      t1(igs(mgs),jy,kgs(mgs)) = scsaci(mgs) +   &
     &                           scsacip(mgs) +   &
     &                           scsacir(mgs)
      t2(igs(mgs),jy,kgs(mgs)) = scglacs(mgs) +  &
     &                           scgmacs(mgs) +  &
     &                           scghacs(mgs) +  &
     &                           scfacs(mgs)
      t3(igs(mgs),jy,kgs(mgs)) =   &
     &                           scglaci(mgs) +  &
     &                           scglacip(mgs) +  &
     &                           scglacir(mgs) +  &
     &                           scgmaci(mgs) +  &
     &                           scgmacip(mgs) +  &
     &                           scgmacir(mgs) +  &
     &                           scghacir(mgs) +  &
     &                           scghacip(mgs) +  &
     &                           scghaci(mgs) +  &
     &                           scfacir(mgs) +  &
     &                           scfacip(mgs) +  &
     &                           scfaci(mgs)
      t4(igs(mgs),jy,kgs(mgs)) = scsacw(mgs)
      t5(igs(mgs),jy,kgs(mgs)) = scglacw(mgs) +  &
     &                           scgmacw(mgs) +  &
     &                           scghacw(mgs) +  &
     &                           scfacw(mgs) +  &
     &                           schacw(mgs) +  &
     &                           schlacw(mgs)
      t6(igs(mgs),jy,kgs(mgs)) =   &
     &                           schacs(mgs) +  &
     &                           schlacs(mgs)
      
      t8(igs(mgs),jy,kgs(mgs)) =   &
     &                           schacir(mgs) +  &
     &                           schacip(mgs) +  &
     &                           schaci(mgs) +  &
     &                           schlacir(mgs) +  &
     &                           schlacip(mgs) +  &
     &                           schlaci(mgs)
      t9(igs(mgs),jy,kgs(mgs)) =  vtxbar(mgs,lhl,1)
      end do
      
      RETURN
      END

#endif


