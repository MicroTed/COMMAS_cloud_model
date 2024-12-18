#ifndef RKIND
#define RKIND 4
#endif


 MODULE ELEC_STUFF
 
 implicit none

#ifdef MPI

#ifdef BOXMG
#include        "BMG_constants.h"
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'
#endif

#endif

      integer ion
      parameter (ion=1)
      
!
! Mudpack variables
!
!      integer imudstop
!      save imudstop
!      integer nzmax,nztop
!      parameter (nzmax1 = 400)
!      real z1d2(nzmax1,4)
!      real dlz
!      real    dzmax
!      save    dzmax
      real  ::  delzd, delzu, r
!      logical istretch
!      save    istretch
!      integer iixps(2),jjyqs(2),kkzrs(2)
!      integer iiexs(2),jjeys(2),kkezs(2)
!      integer*8 llworks(2), n1,n2,n3
!      integer nnxs(2),nnys(2),nnzs(2),nzs
      
      logical   ::  ltmp

      
!      save imud,nbz2,nbw,nbe,nbs,nbn,nztop
!c      save z1d2
!      save iixps,jjyqs,kkzrs,iiexs,jjeys,kkezs,llworks,nnxs,nnys,nnzs
      

!
!    Ion stuff
!
!      real cwmasn,cwmasx
!      save cwmasn,cwmasx
!      real cinccn(nzmax1) !, cwc1(nzmax1)
      
      
!      real uz(nz,4)  ! copy of ion mobilities
!      real cioncp(nz,2) ! copy of fairweather ion concentrations

      real sctot, sctot1n, sctot1p, qtot, sctot2n, sctot2p
      real sctot3n, sctot3p
      real scredisn, scredisp
      real scqmin,sc1
      parameter ( scqmin = 1.0e-10 )

      double precision chgneg,chgpos,chgnet ! total net negative/positive charge in storm (not including corona)
      double precision chgnegcld,chgposcld,chgnetcld ! total net negative/positive charge in storm (not including corona)
      double precision, allocatable :: chgposz(:),chgnegz(:),chgnetz(:)
      double precision, allocatable :: chgpionz(:),chgnionz(:),chgnetionz(:)
      double precision chgneg0,chgpos0
      real cnaneg,cnapos ! total net negative/positive corona charge

!
!  dummy vars for stand-alone lightning
!
      integer ixst,jyst,kzst
      integer numic,numcgn,numcgp,nkzmn(500)
!      save numic,numcgn,numcgp
      character (LEN=80) :: outname

      double precision xfallp, xfalln, xfalltp, xfalltn
      
!
!  arrays used in FISHPAK for solver for eleectric potential field
!
      integer nxb,nyb,nb,nzbt,nzb


! vars for magical corona layer
      
!      real sccna(nx,ny)
!      save sccna
      real ez1,ezmin,ezmax
      real cnamx,cnamn
      real cnmx
      parameter ( cnmx = 1.5e-9 )
      real delc
!      real emaxc
!      parameter (emaxc = 5.0e3)   

!
!  electrical permitivity of air C / (N m**2) -  check the units
!
      real eperao,dslay
      parameter (eperao  = 8.8592e-12 )
      
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      
!
!  screening layer depth
!
      parameter( dslay = 500.0 )

!
!  constants for boundary conditions for poisson solver
!

! other stuff for poisson solver
!      integer ierror
!      real xs, xf, ys, yf, zs, zf  ! domain limits in meters
!      real pertrb,elmbda
      
!
!  WILSON CORONA
!
      real wc1,zctop
      integer wc2
      parameter(wc1 = 0.1, wc2 = 5)
      integer nzctop
      real czcfac 
      
!      save potfair,ezfair,rhofair,ezfairw
      
      double precision ezfairo
      parameter(ezfairo = -80.00d0)
      double precision efa1, efa2, efa3, efb1, efb2, efb3
      parameter(efa1 = 4.5d-3, efa2 = 3.8d-4, efa3 = 1.0d-4, &
                efb1 = 0.50,   efb2 = 0.65,   efb3 = 0.10    )
     
      double precision ezfairs        ! fair field at surface
      parameter(  ezfairs = -100.0d0 )

      double precision dezdz,zgrd

!
! data: 1D
!
      real  bcondc(21)
      real  bconda(21)
      real  bzcond(21)
!
      integer ncond
      parameter (ncond=400)
      real  conda(ncond), condc(ncond)
!
!  polar conductivity
!
      data bconda/5.5e-14,6.e-14,7.e-14,7.8e-14,8.8e-14,  &
       1.e-13,1.1e-13,1.3e-13,1.5e-13,1.7e-13,1.9e-13,    &
       2.1e-13,2.3e-13,2.6e-13,2.8e-13,3.1e-13,3.4e-13,   &
       3.8e-13,4.2e-13,4.5e-13,4.9e-13/


! Energy density

      double precision energynew, energyold
      save energynew, energyold  
      data energynew /0.0d0/    
      
      real energychange

      integer  icgyn  ! =0 no cg, =1 cg happened

      real :: emax = 0.0

  real, allocatable :: uz(:,:) ! uz(nz,4)  ! copy of ion mobilities
  real, allocatable :: cioncp(:,:) ! cioncp(nz,2) ! copy of fairweather ion concentrations

  integer :: iestag = 1 ! NOTE: Setting here is moot.  Must change stag setup in grid_create to use unstaggered
  
      logical :: istretch

      integer :: nzmax,nztop
      real    :: dzmax

      integer :: imud, imudstop
      integer :: nbz2,nbw,nbe,nbs,nbn
      integer :: nbwmg(2),nbemg(2),nbsmg(2),nbnmg(2)

      integer :: iixps(2),jjyqs(2),kkzrs(2)
      integer :: iiexs(2),jjeys(2),kkezs(2)
      integer*8 :: llworks(2)
      integer :: nnxs(2),nnys(2),nnzs(2),nzs

      integer :: nxslm,nyslm,nzslm    !  dimensions for the SLM lightning
      integer :: nxb2,nyb2,nzbt2,nzb2,nb2
      integer :: nxb2s(2),nyb2s(2),nzbt2s(2),nzb2s(2),nb2s(2)

      double precision, allocatable :: sc(:),sc2(:) ! sc(lscmx), sc2(lscmx)
      double precision scion(6),scion2(6)

      real, allocatable ::  ctswinz(:),ctswipz(:)
      real, allocatable ::  ctswwnz(:),ctswwpz(:)
      real, allocatable ::  ctghsnz(:),ctghspz(:)
      real, allocatable ::  ctghinz(:),ctghipz(:)
      real, allocatable ::  ctghwnz(:),ctghwpz(:)
      real, allocatable ::  ctgsnz(:),ctgspz(:)
      real, allocatable ::  ctginz(:),ctgipz(:)
      real, allocatable ::  cthsnz(:),cthspz(:)
      real, allocatable ::  cthinz(:),cthipz(:)

      real, allocatable ::  cgszmn(:),chszmn(:)
      real, allocatable ::  cgizmn(:),chizmn(:)
      real, allocatable ::  cghwzmn(:)
      real, allocatable ::  cgszmx(:),chszmx(:)
      real, allocatable ::  cgizmx(:),chizmx(:)
      real, allocatable ::  cghwzmx(:)
      real, allocatable ::  cswizmn(:),cswizmx(:)
      real, allocatable ::  cswwzmn(:),cswwzmx(:)
      real, allocatable ::  cghszmn(:),cghszmx(:)
      real, allocatable ::  cghizmn(:),cghizmx(:)

!       save imud,nbz2,nbw,nbe,nbs,nbn,nztop
!      save iixps,jjyqs,kkzrs,iiexs,jjeys,kkezs,llworks,nnxs,nnys,nnzs

!c
!c charging rate max/mins
!c
      real cswimn,cswimx
      real cswwmn,cswwmx
      real cghsmn,cghsmx
      real cghimn,cghimx
      real cghwmn,cghwmx
      real cgsmn,cgsmx
      real cgimn,cgimx
      real chsmn,chsmx
      real chimn,chimx

      
      real, allocatable :: lgtth(:,:)
     
!c
!c integrated charging rate pos/neg
!c
      real ctswin,ctswip
      real ctswwn,ctswwp
      real ctghsn,ctghsp
      real ctghin,ctghip
      real ctghwn,ctghwp
      real ctgsn,ctgsp
      real ctgin,ctgip
      real cthsn,cthsp
      real cthin,cthip
      real cthsdn,cthsdp

      integer :: firstsolve = 0

!
!  indices time-height arrays
!
      
      integer, parameter :: chans    = 1
      integer, parameter :: chansp   = 2
      integer, parameter :: chansn   = 3
      integer, parameter :: cgchans  = 4
      integer, parameter :: cgchansp = 5
      integer, parameter :: cgchansn = 6
      integer, parameter :: icchans  = 7
      integer, parameter :: icchansp = 8
      integer, parameter :: icchansn = 9
      integer, parameter :: cinit    = 10
      integer, parameter :: cginit   = 11
      integer, parameter :: icinit   = 12
      integer, parameter :: cgninit  = 13
      integer, parameter :: cgpinit  = 14

#ifdef BOXMG
!
!  Stuff for BOXMG solver (used for MPI)
!
! ------------------------------------------------
!      Multigrid/Workspace Memory Allocation: 
! ------------------------------------------------
 
      !
      ! Workspace pointers and logicals
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)
      LOGICAL   BMG_InWORK(NBMG_InWORK)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSI, pSR 

! -------------------------------------------------
!      Multigrid:  Variables
! -------------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Plane Relaxation
      !
      INTEGER, allocatable, dimension(:) ::  BMG_iWORK_PL
      REAL*RKIND,  allocatable, dimension(:) ::  BMG_rWORK_PL

      !
      ! Workspace: Generic
      !
      INTEGER, allocatable, dimension(:)  :: BMG_iWORK
      REAL*RKIND,  allocatable, dimension(:)  :: BMG_rWORK

      !
      ! Workspace: Coarse-grid Solve
      !
      INTEGER, allocatable, dimension(:)  ::   BMG_iWORK_CS
      REAL*RKIND,  allocatable, dimension(:)  ::   BMG_rWORK_CS

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*RKIND, allocatable, dimension(:) :: Q, QF, SO

      !
      ! Miscellaneous
      !
      INTEGER   NCbmg, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*RKIND    TOL, TOL_SAVE

!      INTEGER NGx, NGy, NGz
      INTEGER   iGs, jGs, kGs, NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER, allocatable, dimension(:) ::  pMSG, pMSGSO

! -------------------------------------------------
!      MPI/MSG:  Variables
! -------------------------------------------------

      INTEGER   NBMG_MSG_iGRIDm, NBMG_MSG_iGRID

      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRID
      
      INTEGER  BMG_MSG_iGRIDdum
      INTEGER  BMG_MSG_pGRID(NBMG_MSG_pGRID)
      
      INTEGER   NMSGi, NMSGr
      

! --------------------------------------------------
!      Local Variables:
! --------------------------------------------------

      INTEGER   NFm, NOGm, NSOm
      INTEGER   NBMG_iWORK, NBMG_rWORK
      INTEGER   NBMG_iWORK_PL, NBMG_rWORK_PL
      INTEGER   NBMG_iWORK_CS, NBMG_rWORK_CS
      
      real*8  dt1, dt2

#endif

 CONTAINS
 
 SUBROUTINE ELEC_STUFF_INIT(nzend)
  implicit none
  integer nzend
        allocate( ctswinz(nzend), &
        ctswipz(nzend), &
        ctswwnz(nzend), ctswwpz(nzend), &
        ctghsnz(nzend), &
        ctghspz(nzend), &
        ctghinz(nzend), &
        ctghipz(nzend), &
        ctghwnz(nzend), &
        ctghwpz(nzend), &
        ctgsnz(nzend), &
        ctgspz(nzend), &
        ctginz(nzend), &
        ctgipz(nzend), &
        cthsnz(nzend), &
        cthspz(nzend), &
        cthinz(nzend), &
        cthipz(nzend), &
        cgszmn(nzend),chszmn(nzend), &
        cgizmn(nzend),chizmn(nzend), &
        cghwzmn(nzend), &
        cgszmx(nzend),chszmx(nzend), &
        cgizmx(nzend),chizmx(nzend), &
        cghwzmx(nzend), &
        cswizmn(nzend),cswizmx(nzend), &
        cswwzmn(nzend),cswwzmx(nzend), &
        cghszmn(nzend),cghszmx(nzend), &
        cghizmn(nzend),cghizmx(nzend), &
        chgposz(nzend),chgnegz(nzend),chgnetz(nzend), &
        chgpionz(nzend),chgnionz(nzend),chgnetionz(nzend) &
         )
  
  
  RETURN
 END SUBROUTINE ELEC_STUFF_INIT
!
 END MODULE ELEC_STUFF

! ####################################################################
! ####################################################################

       subroutine ELEC_DRIVE1(gd,                                   &
     &                     pinit, ab,                                &  ! base state
     &                         p2,                                   &  ! pi'
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
     &                        tstat, lstt, bcx, bcy, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi )

      USE FILE_MODULE
      USE GRID_MODULE
      USE GRIDPARAM_MODULE, only: dx_stretch,dy_stretch
      USE MICRO_MODULE
      USE ELEC_MODULE ! , only: hstretch
      USE ELEC_STUFF
      USE INDEX_MODULE
      USE CPUTIME_MODULE
      USE TRMM_MODULE, only : lpredict
      USE COMMASMPI_MODULE

      implicit none


#ifdef MPI
      INCLUDE "mpif.h"
#endif

      integer nx,ny,nz,na,nor,nba,nv,ia,il

      TYPE(GRID)         :: gd 

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer :: tstat
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

      
      integer iunit
      integer, parameter    :: ng1 = 1

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)

      logical io_flag

      integer ntmul
      real dtp1

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt
      parameter ( nv = 3 )

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)
      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
      integer nx1,ny1,nz1

      real eztop, pottop
      common/bndcz/eztop, pottop

! local variables
      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n
      integer ix,jy,kz
      real x,tx
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

      double precision :: piontot, niontot ! , chgpos, chgneg
      double precision :: potmx, potmn

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

! #ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

!      if (debug_mpi) write(0,*) "DRIVER: ENTERED SUBROUTINE",my_rank


! #endif
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      scion2(:) = 0.0d0
      sc2(:) = 0.0d0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1

      DO kz=kzb,kze
#else
      DO kz=1,nz
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

        imudstop = 0

          imud = 1

! #ifdef MPI

#if defined (MPI) && defined (BOXMG)

       nxslm = itile
       nyslm = jtile
       nzslm = ktile
       
       nnxs(1) = nxslm
       nnys(1) = nyslm
       nnzs(1) = nzslm
       
       write(iunit,*) 'ELEC_DRIVE1: nnxs,nnys,nnzs = ',nnxs(1),nnys(1),nnzs(1)
       
       IF ( .true. ) THEN

         IF ( Abs( dx_stretch - dy_stretch ) > 1. ) THEN
           IF ( my_rank == 0 ) write(0,*) 'dx_stretch /= dy_stretch! These must be equal for lightning.'
           CALL commasmpi_abort()
         ENDIF
         IF ( lightintx .ge. 1 ) THEN
           dslight = dx_stretch/Float(lightintx)
         ELSE
           lightintx = Nint(dx_stretch/dslight)
         ENDIF

       write(iunit,*) 'ELEC_DRIVE1: dslight, dx/lightintx = ',dslight, dx/lightintx
       
       IF ( (dx_stretch > 0.0 .and. dx_stretch /= dx ) .or.  &
            (dy_stretch > 0.0 .and. dy_stretch /= dy )  )  THEN
         write(iunit,*) 'dx_str,dx,dy_str,dy = ',dx_stretch, dx, dy_stretch, dy 
         hstretch = .true.
       ENDIF
       
       nbwmg(:) = 0
       nbemg(:) = 0
       nbsmg(:) = 0
       nbnmg(:) = 0
       
       IF ( ixbeg == nxbeg ) nbwmg(1) = Min( maxbzone, Max(NInt(4000./dx) , 12) )
       nnxs(1) = nnxs(1) + nbwmg(1)
       
       IF ( jybeg == nybeg ) nbsmg(1) = Min( maxbzone, Max(NInt(4000./dy) , 12) )
       nnys(1) = nnys(1) + nbsmg(1)
       

       IF ( ixend == nxend ) nbemg(1) = Min( maxbzone, Max(NInt(4000./dx) , 12) )
       nnxs(1) = nnxs(1) + nbemg(1)
       
       IF ( jyend == nyend ) nbnmg(1) = Min( maxbzone, Max(NInt(4000./dy) , 12) )
       nnys(1) = nnys(1) + nbnmg(1)
       
       write(iunit,*) 'ELEC_DRIVE1: nbwmg,s,e,n = ',nbwmg(1),nbemg(1),nbsmg(1),nbnmg(1)
       
       nxslm = nxslm*lightintx
       nyslm = nyslm*lightintx
!       nzslm = (nz-1)*lightintz + 1

       nnxs(2) = nnxs(1)*lightintx
       nnys(2) = nnys(1)*lightintx

       nbwmg(2) = nbwmg(1)*lightintx
       nbemg(2) = nbemg(1)*lightintx
       nbsmg(2) = nbsmg(1)*lightintx
       nbnmg(2) = nbnmg(1)*lightintx

       
       ENDIF
       
       
       IF ( kzend == nzend ) nnzs(1) = nnzs(1) + Max(NInt(12000./dz) , 16)
       
       
       
!       nnxs(2) = nnxs(1)
!       nnys(2) = nnys(1)
!       nnzs(2) = nnzs(1)

        dzmax = dzz(nz-1)

         IF ( dzz(1) .eq. dzz(nz-1) ) THEN ! dz = constant (unstretched)

           dslightz = dz/Max( 1, NInt(dz/dslight) )

           IF ( lightintz .eq. -1 ) THEN 
             lightintz = Max( 1, NInt( dz/dslight ) )
             write(iunit,*) 'Unstretched vertical grid. Setting lightintz = ',lightintz
           ENDIF
           IF ( lightintz .ge. 1 ) dslightz = dz/Float(lightintz)
           
           nzslm = (nz-1)*lightintz + 1

           istretch = .false.

         ELSE ! vertically stretched

           istretch = .true.

           IF ( lightintz .eq. 2 ) THEN
             nzslm = lightintz*nz
             dslightz = dzmax/lightintz
           
           ELSEIF ( lightintz .gt. 2 ) THEN
             write(0,*) 'lightintz .gt. 2 not supported yet for stretched vertical grid!'
             STOP
           ELSE
             lightintz = 1
             dslightz = dzmax
             nzslm = nz
           ENDIF
         ENDIF

         nnzs(2) = nnzs(1)*lightintz

        write(iunit,*) 'dzmax = ',dzmax

#else

       call mudpacksetup(nx,ny,nz,dx,dy,dz,     &
     &           nyslm,nxslm,nzslm,             &
     &           dxx,dyy,dzz,dzmax,             &
     &           nxb2s,nyb2s,nzbt2s,nzb2s,nb2s, &
     &           iixps,jjyqs,kkzrs,             &
     &           iiexs,jjeys,kkezs,             &
     &           llworks,                       &
     &           nnxs,nnys,nnzs,                &
     &           bcx,bcy,                       &
     &           istretch,imudstop,iunit)

#endif
         
         nztop = nnzs(1)
         
         nz1d = nnzs(1)
         nz1dlgt = nnzs(2)
         
!c
!c  set up extended transformation array
!c

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1

      DO k=kzb,kze
#else
      DO k=1,nz
#endif
        z1d2(k,3) = (dz*gzt(k,3)) ! gt(1,1,k,imapz)    ! MFC
        z1d2(k,4) = (dz*gzt(k,4)) ! 0.5*(gt(1,1,k-1,imapz) + gt(1,1,k,imapz))  ! MFE
        z1d2(k,1) = gzt(k,1) ! gt(1,1,k,mzdist)   ! scalar height
        z1d2(k,2) = gzt(k,2) ! gt(1,1,k,mzdist)-0.5*dz/gt(1,1,k,imapz)  ! W-height
      ENDDO

!c      dzmax = gt(1,1,nz,mzdist) - gt(1,1,nz-1,mzdist)
!c      write(iunit,*) 'dzmax = ',dzmax
      
      IF ( nz .lt. nnzs(1) ) THEN

#ifdef MPI
       kzb = nzend
       kze = nnzs(1)
       if (kzend .ge. nnzs(1)) kze = nnzs(1)-kzbeg+1

       DO k=kzb,kze
#else
       DO k=nz,nnzs(1)
#endif

          z1d2(k,3) = z1d2(nz-1,3)
          z1d2(k,4) = z1d2(nz-1,4)
          z1d2(k,2) = z1d2(k-1,2) + dzmax
          z1d2(k,1) =  z1d2(k,2)+ 0.5*dzmax ! 0.5 * ( z1d2(k,2) + z1d2(k+1,2) ) ! z1d2(k-1,1) + dzmax
!c          z1d2(k,3) = z1d2(nz,3)
!c          z1d2(k,4) = z1d2(nz,4)
!c          z1d2(k,1) = z1d2(k-1,1) + dzmax
!c          z1d2(k,2) = z1d2(k-1,2) + dzmax
       ENDDO

#ifdef MPI
       IF ( my_rank .eq. 0 ) THEN
       kzb = nzend
       kze = nnzs(1)
       if (kzend .ge. nnzs(1)) kze = nnzs(1)-kzbeg+1

       DO k=kzb,kze
#else
       DO k=nz,nnzs(1)
#endif
        write(iunit,                                            &
     &   '(1x,i3,1x,f7.0,1x,f7.0,1x,f7.4,1x,f7.4,2(1x,f7.4))')  &
     &    k,                                                    &
     &   z1d2(k,1),z1d2(k,2),z1d2(k,3),z1d2(k,4),               &
     &   (z1d2(k,3)+z1d2(max(1,k-1),3))/2.,                     &
     &   (z1d2(k,4)+z1d2(min(nz,k+1),4))/2.
       ENDDO 
#ifdef MPI
       ENDIF ! my_rank .eq. 0
#endif
      ENDIF !( nz .lt. nnzs(1) )
      
      nzs = nnzs(1)

#ifdef MPI
      kzb = nzend-2
      kze = nzs-1
      if (kzend .ge. nzs-1) kze = nzs-kzbeg

      DO k = kzb,kze
#else
      DO k = nz-2,nzs-1
#endif
                      z1d2(k,3) = dz / ( z1d2(k+1,2) - z1d2(k,2) )   ! MFC(k)
       IF( k .ne. 1 ) z1d2(k,4) = dz / ( z1d2(k,1)   - z1d2(k-1,1) ) ! MFE(k)

      ENDDO
      z1d2(nzs,3) = dz / (2.0*(z1d2(nzs,2) - z1d2(nzs-1,1)))
      z1d2(nzs,4) = dz / (2.0*(z1d2(nzs,2) - z1d2(nzs-1,1)))

      dlz = (z1d2(nnzs(1)-1,1) + z1d2(1,1))/(nnzs(1)-1)
      write(iunit,*) 'dlz = ',dlz
      
      
      write(iunit,*) 'recompute stretch factors for extended grid'

#ifdef MPI
      kzb = 1
      kze = nnzs(1)
      if (kzend .ge. nnzs(1)) kze = nnzs(1)-kzbeg+1

      DO k=kzb,kze
#else
      DO k=1,nnzs(1)
#endif
        DO j=3,4
          z1d2(k,j) = z1d2(k,j)*dlz/dz
        ENDDO
        
        ! now dz = dlz/z1d2(k,3)
        write(iunit,    &
     &  '(1x,i3,1x,f10.3,1x,f10.3,1x,f10.7,1x,f10.7,2(1x,f10.7),1x,f10.3)') &
     &    k,                        &
     &   z1d2(k,1),z1d2(k,2),z1d2(k,3),z1d2(k,4), &
     &   (z1d2(k,3)+z1d2(max(1,k-1),3))/2.,       &
     &   (z1d2(k,4)+z1d2(min(nnzs(1),k+1),4)*dlz/dz)/2., dlz/z1d2(k,3)
      ENDDO
!        gt(ix,jy,kz,imapz)
        
        
!        ENDIF ! imud .eq. 1
        
        
        call flush(iunit)
        write(iunit,*) 'Max Corona density set at ',cnmx

!c 
!c  Set up Fair weather field
!c
      
      IF (  lsce .gt. 1 .and. lpredict ) THEN ! {

!      CALL cld_cpu('MICROPHYSICS')
!      CALL cld_cpu('ELECTRICITY')

      IF ( nzfair .lt. nztop ) THEN !{
        write(iunit,*)                                       &
     &  'STOP in sam.08.driver.x, nzfair must be > nztop!',  &
     &    ' nzfair,nz=',nzfair,nztop
        write(0,*) 'STOP in sam.08.driver.x, nzfair must be > nztop!',  &
     &    ' nzfair,nz=',nzfair,nztop
      STOP
      ENDIF ! }
      
      IF ( imud .eq. 1 ) THEN ! {
      write(iunit,*) 'kz, z, potfair, ezfairw,ezfairs,muz(3)%flt1d(kz): '
      do kz = 1,nztop
       IF ( fairweather ) THEN
!c      zgrd = .5*dz + (kz-1)*dz
      zgrd = z1d2(kz,1)
      ezfair(kz) = ezfairo *            &
     &           ( efb1*exp(-efa1*zgrd) &
     &            +efb2*exp(-efa2*zgrd) &
     &            +efb3*exp(-efa3*zgrd) )
      
       dezdz = ezfairo *                      &
     &           ( -efa1*efb1*exp(-efa1*zgrd) &
     &             -efa2*efb2*exp(-efa2*zgrd) &
     &             -efa3*efb3*exp(-efa3*zgrd) )
       
       rhofair(kz) = dezdz*eperao

      zgrd = z1d2(kz,2)
      ezfairw(kz) = ezfairo *           &
     &           ( efb1*exp(-efa1*zgrd) &
     &            +efb2*exp(-efa2*zgrd) &
     &            +efb3*exp(-efa3*zgrd) )
       ELSE
         ezfair(kz) = 0.0
         rhofair(kz) = 0.0
         ezfairw(kz) = 0.0
       ENDIF
      end do

#ifdef MPI
      kzb = 1
      kze = nztop-kd1+1
!      if (kzend .ge. nztop-kd1) kze = nztop-kd1-kzbeg+1

      do kz = kzb,kze
#else
      do kz = 1,nztop-kd1
#endif
      IF ( fairweather ) THEN
      if ( kz.eq.1 ) then
! erm 7/24/2008 changed factor from 0.5 to 1.0 to get consistent result as from 
!  Poisson solver, which "thinks" the ground is really at k=-1
      potfair(kz) = -(1.0)*(ezfair(1)+ezfairs)*dlz*(0.5)/ z1d2(kz,3)
      end if
      if ( kz.gt.1 ) then
      potfair(kz) = potfair(kz-1) - (0.5)*(ezfair(kz)+ezfair(Max(1,kz-1)))*dlz/z1d2(kz,3)
      end if
      
      ELSE
        potfair(kz) = 0.0
      ENDIF ! fairweather
      
       IF ( kz .le. nz) THEN
        write(iunit,*) kz,z1d2(kz,1),potfair(kz),ezfairw(kz),  &
     &   ezfair(kz),rhofair(kz),muz(1)%flt1d(kz) ! ,uz(kz,3), uz(kz,2),uz(kz,4)
       ELSE
        write(iunit,*) kz,z1d2(kz,1),potfair(kz),ezfairw(kz), &
     &   ezfair(kz),rhofair(kz)
       ENDIF
      end do  
       
      ELSE ! go here if imud .ne. 1 

#ifdef MPI
      kzb = 1
      kze = nztop-kd1
!      if (kzend .ge. nztop) kze = kzend-nztop+1

      do kz = kzb,kze
#else
      do kz = 1,nztop
#endif
      zgrd = .5*dz + (kz-1)*dz
      IF ( fairweather ) THEN
      ezfair(kz) = ezfairo *            &
     &           ( efb1*exp(-efa1*zgrd) &
     &            +efb2*exp(-efa2*zgrd) &
     &            +efb3*exp(-efa3*zgrd) )
      ELSE
       ezfair(kz) = 0.0
      ENDIF
      end do

#ifdef MPI
      kzb = 1
      kze = nztop-kd1
!      if (kzend .ge. nztop-kd1) kze = kzend-(nztop-kd1)+1

      do kz = kzb,kze
#else
      do kz = 1,nztop-kd1
#endif
      IF ( fairweather ) THEN
      if ( kz.eq.1 ) then
      potfair(kz) = -(0.5)*(ezfair(1)+ezfairs)*dz*(0.5)
      end if
      if ( kz.gt.1 ) then
      potfair(kz) = potfair(kz-1) - (0.5)*(ezfair(kz)+ezfair(kz-1))*dz
      end if 

      ELSE
        potfair(kz) = 0.0
      ENDIF ! fairweather

      end do  

      ENDIF ! }imud .eq. 1
      
      eztop = ezfair(nztop)
      pottop = potfair(nztop)

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

         nbwmg(1) = nbw
         nbemg(1) = nbe
         nbsmg(1) = nbs
         nbnmg(1) = nbn

         nbwmg(2) = nbwmg(1)*lightintx
         nbemg(2) = nbemg(1)*lightintx
         nbsmg(2) = nbsmg(1)*lightintx
         nbnmg(2) = nbnmg(1)*lightintx

!         nbe = (nnxs(1) - nx - 1)/2
!         nbw = nnxs(1) - nx - 1 - nbe
!         nbn = (nnys(1) - ny - 1)/2
!         nbs = nnys(1) - ny - 1 - nbn
#endif
         
         IF ( ny .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF
!         IF ( bcx .eq. 2 ) THEN ; nbw = -1 ; nbe = 0 ; ENDIF
!         IF ( bcy .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

      write(iunit,*) 'nztop,nbw,nbe,nbs,nbn,eztop,bcx,bcy = ', nztop,nbw,nbe,nbs,nbn,eztop,bcx,bcy
      write(iunit,*) 'gx,gy: ',gx(-nor+1),gy(1)

        IF ( .not. allocated( x1d2 ) ) THEN
!          allocate( x1d2(-nbw:nx+2+nbe) )
!          allocate( y1d2(-nbs:ny+2+nbn) )
          allocate( x1d2(0:nx+2+nbe+nbw) )
          allocate( y1d2(0:ny+2+nbn+nbs) )
!      real gx(-nor+1:nx+nor)
!      real gy(-nor+1:ny+nor)
          DO i = -nbw,nx+2+nbe
            IF ( i >= -nor+1 .and. i <= nx+nor ) THEN
              x1d2(i+nbw) = 1.0/gx(i)
            ELSEIF ( i < -nor+1 ) THEN
              x1d2(i+nbw) = 1.0/gx(-nor+1)
            ELSE
              x1d2(i+nbw) = 1.0/gx(nx+nor)
            ENDIF
            write(iunit,*) 'i,x1d2 = ',i+nbw,x1d2(i+nbw)
          ENDDO
          
          IF ( ny > 2 ) THEN
          DO j = -nbs,ny+2+nbn
            IF ( j >= -nor+1 .and. j <= ny+nor ) THEN
              y1d2(j+nbs) = 1.0/gy(j)
            ELSEIF ( j < -nor+1 ) THEN
              y1d2(j+nbs) = 1.0/gy(-nor+1)
            ELSE
              y1d2(j+nbs) = 1.0/gy(ny+nor)
            ENDIF
            write(iunit,*) 'j,y1d2 = ',j+nbs,y1d2(j+nbs)
          
          ENDDO
          ELSE
             y1d2(:) = 1.0/gy(1)
          ENDIF
        
        ENDIF
        
        IF ( .not. allocated( x1d2lgt ) ) THEN

       ! subdivide each dynamics grid box into lightintx pieces of 
       ! of equal size dx/lightnintx. Not worrying about stretch within each
       ! large dx.

!          allocate( x1d2lgt(-nbwmg(2):nxslm+1+nbemg(2) ) )
!          allocate( y1d2lgt(-nbsmg(2):nyslm+1+nbnmg(2) ) )
!          allocate( x1d2lgt(-nbwmg(2):lightintx*(nx+2+nbe) ) )
!          allocate( y1d2lgt(-nbsmg(2):lightintx*(ny+2+nbn) ) )
          allocate( x1d2lgt(0:lightintx*(nx+2+nbe)+nbwmg(2) ) )
          allocate( y1d2lgt(0:lightintx*(ny+2+nbn)+nbsmg(2) ) )
          
!      real gx(-nor+1:nx+nor)
!      real gy(-nor+1:ny+nor)
          DO i = -nbw,nx+1+nbe
            x1d2lgt(i*lightintx+nbwmg(2)) = x1d2(i+nbw)/lightintx
            
            DO n = 1,lightintx
              x1d2lgt(i*lightintx+n+nbwmg(2)) = x1d2(i+nbw)/lightintx
            ENDDO
          ENDDO
          
          DO i = -nbwmg(2),nxslm+1+nbemg(2)
            write(iunit,*) 'i,x1d2lgt = ',i+nbwmg(2),x1d2lgt(i+nbwmg(2))
          ENDDO
          
          IF ( ny > 2 ) THEN
          DO j = -nbs,ny+1+nbn
            y1d2lgt(j*lightintx+nbsmg(2)) = y1d2(j+nbs)/lightintx
            DO n = 1,lightintx
              y1d2lgt(j*lightintx+n+nbsmg(2)) = y1d2(j+nbs)/lightintx
            ENDDO
          ENDDO

          DO j = -nbsmg(2),nyslm+1+nbnmg(2)
            write(iunit,*) 'j,y1d2lgt = ',j+nbsmg(2),y1d2lgt(j+nbsmg(2))
          ENDDO
          ENDIF
        
        ENDIF

        dlzlgt = dlz/lightintz
        
        IF ( nnzs(1) .eq. nnzs(2) ) THEN
         potfairlgt(1:nnzs(1)) = potfair(1:nnzs(1))
         ezfairlgt(1:nnzs(1))  = ezfair(1:nnzs(1))
         rhofairlgt(1:nnzs(1)) = rhofair(1:nnzs(1))
         ezfairwlgt(1:nnzs(1)) = ezfairw(1:nnzs(1))
         
         DO k = 1,4
           z1d2lgt(1:nztop,k) = z1d2(1:nztop,k)
         ENDDO
         
        ELSEIF ( lightintz .ge. 2 ) THEN ! subdivide grid and interpolate to midpoints

            write(iunit,*) 'set up subdivided vertical grid'
        
        IF ( lightintz .eq. 2 ) THEN  
!       First set up heights.  Subdivide between w-points (except first)
#ifdef MPI
        kzb = 1
        kze = nztop
!        if (kzend .eq. nzend) kze = nztop

        DO k=kzb,kze
#else
        DO k=1,nztop
#endif
!           z1d2lgt(2*k,1) = z1d2(k,1)
           z1d2lgt(2*k-1,2) = z1d2(k,2)
!           z1d2lgt(2*k,1) = gzt(k,1) ! gt(1,1,k,mzdist)   ! scalar height
!           z1d2lgt(2*k,2) = gzt(k,2) ! gt(1,1,k,mzdist)-0.5*dz/gt(1,1,k,imapz)  ! W-height
        
           IF ( k .eq. 1 ) THEN
!            z1d2lgt(2*k-1,1) =  z1d2(k,1)
            z1d2lgt(2*k-1,2) =  z1d2(k,2)
            z1d2lgt(2*k  ,2) =  z1d2(k+1,2)/lightintz
            z1d2lgt(2*k-1,1) =  z1d2lgt(2*k  ,2)/lightintz
           ELSE
!             IF ( k .eq. 2 ) THEN  ! check for stretching
!              delzu = z1d2(2,1) - z1d2(1,1)
!              delzd = 2*z1d2(1,1)
!             ELSE
!              delzu = z1d2(k  ,1) - z1d2(k-1,1)
!              delzd = z1d2(k-1,1) - z1d2(k-2,1)
!             ENDIF

              delzu = z1d2(k+1,2) - z1d2(k  ,2)
              delzd = z1d2(k  ,2) - z1d2(k-1,2)

! stretch factor for ilightintz = 2, where x is the subdivided stretch factor and r is 
! the full stretch factor:
!  x = (-1.0 + Sqrt(  1.0 + 8.0*r ) ) / 2.0
!  reduces to x = r when r = 1.0 (unstretched region)
            
              r = delzu/delzd
              IF ( k .eq. nztop ) r = 1.0
              x = (-1.0 + Sqrt(  1.0 + 8.0*r ) ) / 2.0
              write(iunit,*) 'k,r,x = ',k,r,x, delzu,delzd
              z1d2lgt(2*k,2) = z1d2(k,2) + x*( delzd )/lightintz
!              z1d2lgt(2*k-1,1) = z1d2(k-1,1) + x*( delzd )/lightintz
!              z1d2lgt(2*k-1,2) = z1d2(k-1,2) + x*( z1d2(k-1,2) - z1d2(k-2,2) )/lightintz
 
              z1d2lgt(2*k-1,1) =  0.5*( z1d2lgt(2*k  ,2) + z1d2lgt(2*k-1,2) )
              z1d2lgt(2*k-2,1) =  0.5*( z1d2lgt(2*k-1,2) + z1d2lgt(2*k-2,2) )

           ENDIF

            write(iunit,*) k, z1d2(k  ,1), z1d2lgt(Max(1,2*k-2)  ,1),  z1d2(k  ,2), z1d2lgt(2*k-1  ,2)
            write(iunit,*) k, z1d2(k  ,1), z1d2lgt(2*k-1    ,1),  z1d2(k  ,2), z1d2lgt(2*k    ,2)

!            write(iunit,*) k,   z1d2(k  ,2), z1d2lgt(2*k-1  ,2)
!            write(iunit,*) k,   z1d2(k  ,2), z1d2lgt(2*k    ,2)
            
        ENDDO
        
        ELSE !  ( lightintz .gt. 2 ) --> must be unstretched grid!

      n=1
      DO k = 1,nz1dlgt ! lightintz*nztop-1
             
        if (Mod(k,lightintz) .eq. 0) n=n+1
        z1d2lgt(k,2) = dslightz*(k-1)
        z1d2lgt(k,1) = z1d2lgt(k,2) + dslightz/2

      ENDDO
        
        
        
        ENDIF
        


!       Now set up dzc and dze
#ifdef MPI
        kzb = 1
        kze = nz1dlgt-1 ! lightintz*nztop - 1
        if (kzend .eq. nzend) kze = nz1dlgt - 1 ! lightintz*nztop - 1

        DO k=kzb,kze
#else
        DO k=1,nz1dlgt-1 ! lightintz*nztop - 1
#endif
        
        
              z1d2lgt(k,3) = dlzlgt / ( z1d2lgt(k+1,2) - z1d2lgt(k,2) )   ! MFC(k)
              IF( k .ne. 1 ) THEN
                z1d2lgt(k,4) = dlzlgt / ( z1d2lgt(k,1)   - z1d2lgt(k-1,1) ) ! MFE(k)
              ELSE
                z1d2lgt(k,4) = z1d2lgt(k,3) ! MFE(k)
              ENDIF
              
              IF ( k .eq. nz1dlgt-1 ) THEN ! lightintz*nztop - 1 ) THEN
                z1d2lgt(k+1,4) = z1d2lgt(k,4)
                z1d2lgt(k+1,3) = z1d2lgt(k,3) 
              ENDIF

        
        ENDDO


!       Now print up dzc and dze
#ifdef MPI
        kzb = 1
!        kze = kzend
        kze = nztop
        if (kzend .eq. nzend) kze = nztop

        DO k=kzb,kze
#else
        DO k=1,nztop
#endif
        
            write(iunit,*) k, z1d2(k  ,3), z1d2lgt(2*k-1  ,3),  z1d2(k  ,4), z1d2lgt(2*k-1  ,4)
            write(iunit,*) k, z1d2(k  ,3), z1d2lgt(2*k    ,3),  z1d2(k  ,4), z1d2lgt(2*k    ,4)
       
        ENDDO


#ifdef MPI
        kzb = 1
!        kze = kzend
        kze = nz1dlgt ! lightintz*nztop
        if (kzend .eq. nzend) kze = nz1dlgt ! lightintz*nztop

        DO kz=kzb,kze
#else
        DO kz=1,nz1dlgt ! lightintz*nztop
#endif
       IF ( fairweather ) THEN
!c      zgrd = .5*dz + (kz-1)*dz
      zgrd = z1d2lgt(kz,1)
      ezfairlgt(kz) = ezfairo *          &
     &           ( efb1*exp(-efa1*zgrd)  &
     &            +efb2*exp(-efa2*zgrd)  &
     &            +efb3*exp(-efa3*zgrd) )
      
       dezdz = ezfairo *                       &
     &           ( -efa1*efb1*exp(-efa1*zgrd)  &
     &             -efa2*efb2*exp(-efa2*zgrd)  &
     &             -efa3*efb3*exp(-efa3*zgrd) )
       
       rhofairlgt(kz) = dezdz*eperao

      zgrd = z1d2lgt(kz,2)
      ezfairwlgt(kz) = ezfairo *          &
     &           ( efb1*exp(-efa1*zgrd)   &
     &            +efb2*exp(-efa2*zgrd)   &
     &            +efb3*exp(-efa3*zgrd) )
       ELSE
         ezfairlgt(kz) = 0.0
         rhofairlgt(kz) = 0.0
         ezfairwlgt(kz) = 0.0
       ENDIF
      end do

      write(iunit,*) 'Fairweather on lightning grid, dlzlgt = ',dlzlgt
      write(iunit,*) 'kz, potfair, ezfair, rhofair'
#ifdef MPI
        kzb = 1
!        kze = kzend
        kze = nz1dlgt-1 ! lightintz*nztop - 1
        if (kzend .eq. nzend) kze = nz1dlgt-1 ! lightintz*nztop - 1

        DO kz=kzb,kze
#else
        DO kz=1,nz1dlgt-1 ! lightintz*nztop - kd1
#endif
      
      IF ( fairweather ) THEN
      if ( kz.eq.1 ) then
! erm 7/24/2008 changed factor from 0.5 to 1.0 to get consistent result as from 
!  Poisson solver, which "thinks" the ground is really at k=-1
      potfairlgt(kz) = -(1.0)*(ezfairlgt(1)+ezfairs)*dlzlgt*(0.5)/   &
     &         z1d2lgt(kz,3)
      end if
      if ( kz.gt.1 ) then
      potfairlgt(kz) = potfairlgt(kz-1)                            &
     &             -(0.5)*(ezfairlgt(kz)+ezfairlgt(kz-1))*dlzlgt/  &
     &         z1d2lgt(kz,3)
      end if 
      
      ELSE
        potfairlgt(kz) = 0.0
      ENDIF ! fairweather
      
        write(iunit,*) kz,z1d2lgt(kz,1),potfairlgt(kz),ezfairwlgt(kz),  &
     &   ezfairlgt(kz),rhofairlgt(kz)
      end do  
        
!        STOP
        
        ENDIF


!      CALL cld_cpu('ELECTRICITY')
       
       ENDIF ! (lsce .gt. 1)
!c         
!c      DO j=1,ny
!c      DO i=1,nx
!c        sccna(i,j) = 0.0
!c      END DO
!c      END DO
!
!c      DO k=lscw,lschl
!c
!c
!c  set xfalltot to zero in main3 now
!c      DO k=1,na
!c      DO j=1,ny
!c      DO i=1,nx
!c        xfalltot(i,j,k) = 0.0
!c      ENDDO
!c      ENDDO
!c      ENDDO
!      
!c      DO kz=-nor,nz+nor
!c      DO jy=-nor,ny+nor
!c      DO ix=-nor,nx+nor
!c        IF ( dc(ix,jy,kz) .gt. 0.0 ) THEN
!c          dinv(ix,jy,kz) = 1.0/dc(ix,jy,kz)
!c        ELSE
!c          dinv(ix,jy,kz) = 1.0
!c        ENDIF
!c      ENDDO
!c      ENDDO
!c      ENDDO
      
      
      
!c write settings for magical corona:
      
      IF ( ipelec .ge. 1 .and. lpredict ) THEN !{
!       CALL cld_cpu('ELECTRICITY')
     
      write(iunit,*) 'ibg = ',ibg
      IF (ibg.eq.1) THEN
        write(iunit,*) 'setting potential = 0.0 at k=1 (z = +dz/2)'
      ELSEIF (ibg .eq. -1) THEN
        write(iunit,*) 'setting potential = 0.0 at k=0 (z = -dz/2)'
        write(iunit,*) 'corona onset field, emaxc = ',emaxc
      END IF
!c regular grid:

#ifdef BOXMG

         piontot = 0.0d0
         niontot = 0.0d0
         chgpos = 0.0d0
         chgneg = 0.0d0
         potmx = -1.e32
         potmn =  1.e32

         nbw = nbwmg(1)
         nbe = nbemg(1)
         nbs = nbsmg(1)
         nbn = nbnmg(1)


      IF ( time_real .gt. dt ) THEN !{

      t0(:,:,:) = 0.0

!c C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
!c C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!c C$OMP+ REDUCTION (+ : chgpos,chgneg)
      
      
      
      DO ia = lscb,lsceq

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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO


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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
       dv = dxx(ix)*dyy(jy)*dzz(kz)

!        IF ( t0(ix,jy,kz) .ge. 0.0 ) THEN
!          chgpos = chgpos + t0(ix,jy,kz)*dv
!        ELSE
!          chgneg = chgneg + t0(ix,jy,kz)*dv
!        ENDIF

       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) &
     &                             + ec*(an(ix,jy,kz,lscpli)-an(ix,jy,kz,lscnli))
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))
       ENDIF
!        piontot = piontot + an(ix,jy,kz,lscpi)
!        niontot = niontot + an(ix,jy,kz,lscni)
        
!        potmx = Max( potmx, elec(ipot)%flt3d(ix,jy,kz) )
!        potmn = Min( potmn, elec(ipot)%flt3d(ix,jy,kz) )
      
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
!       mpitotindp(1)  = chgpos
!       mpitotindp(2)  = chgneg
!       mpitotindp(3)  = piontot 
!       mpitotindp(4)  = niontot
      
!      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      

      IF ( my_rank == 0 ) THEN
!       chgpos = mpitotoutdp(1)
!       chgneg = mpitotoutdp(2)   
      ENDIF
#endif
      
!       IF ( my_rank == 0 ) THEN
!         write(iunit,*) 'domain chgpos,chgneg,net = ',mpitotoutdp(1),mpitotoutdp(2),mpitotoutdp(1)-mpitotoutdp(2)
!         write(iunit,*) 'domain piontot,niontot,net = ',mpitotoutdp(3),mpitotoutdp(4),mpitotoutdp(3)-mpitotoutdp(4)
!       ENDIF

       IF ( my_rank == 0 ) THEN
!        write(iunit,*) 'chgpos,chgneg = ',chgpos,chgneg,chgpos-chgneg
!        write(iunit,*) 'piontot,niontot = ',piontot, niontot,piontot-niontot
        write(iunit,*) 'call mgsolve 1: nbw,nbe,nbs,nbn = ',nbw,nbe,nbs,nbn
        write(iunit,*) 'number_of_processes = ',number_of_processes
       ENDIF

#ifdef MPI
! find global integrated rate max
!       mpitotindp(1)  = potmx
!       mpitotindp(2)  = -potmn
      
!      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)
      

#endif
!      IF ( my_rank == 0 ) THEN
!        write(iunit,*) 'domain potential max/min = ',mpitotoutdp(1),-mpitotoutdp(2)
!      ENDIF
!      write(iunit,*) 'local potential max/min = ',potmx,potmn

      CALL cld_cpu('BOXMG')
      CALL cld_cpu('BOXMG-INITDG')
       
      CALL mgsetupdg(nx,ny,nz,nor,    &
     &   dx,dy,dz,gxt,gyt,gzt,        &
     &   dxx,dyy,dzz,                 &
     &   nnxs,nnys,nnzs,              &
     &   nbw,nbe,nbs,nbn,             &
     &   istretch,                    &
     &   iunit, bcx, bcy )

!c         write(iunit,*) 'driver0: call putf'
         
!c         IF ( my_rank == 0 ) THEN
!c           DO k = 1,NLzdg
!c             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!c           ENDDO
!c         ENDIF
         
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,   &
     &               nbw,nbe,nbs,nbn,                                 &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,        &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 1               &
     &              )
!c         write(iunit,*) 'done putf'

         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()
         

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!c         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg 0'

       i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = 1.e-5
       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 12
       firstsolve = 1
         IF (my_rank .eq. 0) THEN
            write(iunit,*) ' boxmg isetup = ', BMG_iPARMSdg(id_BMG3_SETUPdg)
         ENDIF
         
         
         ! Do not do a solve here, but leave BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) 
         ! set to 0 so that setup will get done when solve is called at the regular time
         IF ( .true. ) THEN

         CALL commas_SymStd_SOLVE_boxmgdg(                                     &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg, &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                      &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,  &
     &             SOdg, NSOdg,                                                   &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                    &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                      &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,               &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                          &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                          &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                          &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                          &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                  &
     &             pMSGdg, pMSGSOdg,                                         &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,         &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg, &
     &             my_comm                                                               &
     &             )
!c
!c Turn off future setup:
         BMG_iPARMSdg(id_BMG3_SETUPdg) = 2

         ENDIF

       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = bmg_tol

!c         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         
         
         dt2 = MPI_Wtime()

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for firstsolve = ', dT2 - dT1
         ENDIF

!c      IF ( my_rank .eq. 0 ) THEN
!c          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!c     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!c     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!c     :              )
!c      ENDIF

         
         IF ( .false. ) THEN ! leave the elec array alone when restarting
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,     &
     &               gxt,gyt,gzt,                   &
     &               nbw,nbe,nbs,nbn,iunit,         &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz &
     &              )
         ENDIF

      CALL cld_cpu('BOXMG-INITDG')

      CALL cld_cpu('BOXMG')
     
      ENDIF !}

#else \* ifdef BOXMG *\

         nbw = (nnxs(1) - nx - 1)/2
         nbe = nnxs(1) - nx - 1 - nbw
         nbs = (nnys(1) - ny - 1)/2
         nbn = nnys(1) - ny - 1 - nbs

!         nbe = (nnxs(1) - nx - 1)/2
!         nbw = nnxs(1) - nx - 1 - nbe
!         nbn = (nnys(1) - ny - 1)/2
!         nbs = nnys(1) - ny - 1 - nbn
#endif


         IF ( ny .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF
         IF ( bcx .eq. 2 ) THEN ; nbw = -1 ; nbe = 0 ; ENDIF
         IF ( bcy .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         nbz2 = nnzs(1) - nz
         nzbt2 = nz + nbz2
         nzb2 = nzbt2 + 1
         nxb = nnxs(1)
         nyb = nnys(1)
!c
!c
      IF ( time_real .eq. dt  ) THEN !{
      
      t0(:,:,:) = 0.0
      
      IF ( fairweather ) THEN
!c
!c  Set up what the model thinks is the fair-weather Ez
!c
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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
      t0(ix,jy,kz) =  ( ec*(cion(1)%flt1d(kz)-cion(2)%flt1d(kz))  )
      elec(ipot)%flt3d(ix,jy,kz) = potfair(kz)
      end do
      end do
      end do
      
      ENDIF ! fairweather
      
!      print*, 'INITIAL EFIELD SOLUTION'


#ifdef BOXMG

       
      CALL cld_cpu('BOXMG-INITDG')
!       IF ( my_rank == 0 ) THEN
        write(iunit,*) 'call mgsolve 2: nbw,nbe,nbs,nbn = ',nbw,nbe,nbs,nbn
        write(iunit,*) 'number_of_processes = ',number_of_processes
!       ENDIF

       
      CALL mgsetupdg(nx,ny,nz,nor, &
     &   dx,dy,dz,gxt,gyt,gzt,     &
     &   dxx,dyy,dzz,              &
     &   nnxs,nnys,nnzs,           &
     &   nbw,nbe,nbs,nbn,          &
     &   istretch,                 &
     &   iunit, bcx, bcy )



         
!c         write(iunit,*) 'driver1: call putf'
         
!c         IF ( my_rank == 0 ) THEN
!c           DO k = 1,NLzdg
!c             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!c           ENDDO
!c         ENDIF
         
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,   &
     &               nbw,nbe,nbs,nbn, &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 1    &
     &              )
!c         write(iunit,*) 'done putf'


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!c         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg 1'

!c       BMG_IOFLAGdg(iBMG3_BUG_PARAMETERS) = .true.
       i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = 1.e-5
       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 12
       firstsolve = 1

#ifdef __APPLE__
!c        BMG_iPARMSdg(id_BMG3_CG_COMM) = BMG_CG_ALLGATHER
#endif


         CALL commas_SymStd_SOLVE_boxmgdg(                                        &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg, &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                      &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,  &
     &             SOdg, NSOdg,                                                   &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                    &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                      &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,               &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                          &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                          &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                          &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                          &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                  &
     &             pMSGdg, pMSGSOdg,                                         &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,         &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg, &
     &             my_comm                                        &
     &             )

       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = bmg_tol

!c         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dt2 = MPI_Wtime()

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for boxmg first = ', dT2 - dT1
         ENDIF

!c      IF ( my_rank .eq. 0 ) THEN
!c          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!c     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!c     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!c     :              )
!c      ENDIF

!c
!c Turn off future setup:
         BMG_iPARMSdg(id_BMG3_SETUPdg) = 2
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,     &
     &               gxt,gyt,gzt,                  &
     &               nbw,nbe,nbs,nbn,iunit,        &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz             &
     &              )

      CALL cld_cpu('BOXMG-INITDG')

#else


      IF ( fairweather ) THEN

      IF ( number_of_processes .eq. 1 ) THEN
      CALL cld_cpu('EFIELD')
      DO kz=1,12

      ltmp = .false.
      IF ( kz .eq. 12 ) ltmp = .true.  ! calculate efield on the last step if iestag=1

      CALL mud(nx,ny,nz,nor,   &
     &   dx,dy,dz,gxt,gyt,gzt,   &
     &   nnxs(1),nnys(1),nnzs(1), &
     &   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1), &
     &   nbw,nbe,nbs,nbn,          &
     &   llworks(1),istretch,      &
     &   t0,ibg,                   &
     &   elec, iestag, iunit, ltmp, bcx, bcy )

      ENDDO

      IF ( iestag .eq. 0 ) THEN

      call efield                      &
     &  (nx,ny,nz                      &
     &  ,id1,jd1,kd1,istag,jstag,kstag &
     &  ,dx,dy,dz,elec                 &
     &  ,gxt(1,3),gyt(1,3),gzt(1,3)    &
     &  ,iunit)
      ENDIF

      CALL cld_cpu('EFIELD')
      ENDIF ! number_of_processes .eq. 1
      
      ENDIF


#endif


      
      ENDIF !}
!      CALL cld_cpu('ELECTRICITY')
!      CALL cld_cpu('MICROPHYSICS')
      
      IF ( ilight == 3 ) THEN
      ! since MSZ lightning is on the dynamics grid, need to set lightning for same grid
        iGslg0 = ixbeg - 1
        jGslg0 = jybeg - 1
        kGslg0 = kzbeg - 1
      
      ENDIF
      
      ENDIF !} ipelec .ge. 1 .and. lpredict



       RETURN
       END

! ####################################################################


       subroutine ELEC_DRIVE2(gd,                                   &
     &                     pinit, ab,                                &  ! base state
     &                         p2,                                   &  ! pi'
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
     &                        tstat, lstt, bcx, bcy, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,chgadvtemp )

      USE FILE_MODULE
      USE VIS5D_MODULE
      USE GRID_MODULE
      USE MICRO_MODULE
      USE ELEC_MODULE
      USE ELEC_STUFF
      USE INDEX_MODULE
      USE CPUTIME_MODULE
      USE TRMM_MODULE, only : ibsd,jbsd,iesd,jesd
      USE COMMASMPI_MODULE

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
#endif

      integer nx,ny,nz,na,nor,nba,nv,ia,il

      TYPE(GRID)         :: gd 

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer :: tstat
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

      
      integer iunit
      integer, parameter    :: ng1 = 1

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real chgadvtemp(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,2)

      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)

      logical io_flag

      integer ntmul
      real dtp1

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt
      parameter ( nv = 3 )

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)
      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
      integer nx1,ny1,nz1

! local variables
      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n
      integer ix,jy,kz
      real x,tx,tmp
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      logical :: debug_mpi

      integer, parameter :: ntot = 500
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
!
! ###########################################################################
!
      t0(:,:,:) = 0.0

      chgneg = 0.0d0
      chgpos = 0.0d0
      chgnet = 0.0d0
! C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
! C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
! C$OMP+ REDUCTION (+ : chgpos,chgneg,chgnet)
      
      
      DO ia = lscb,lsceq
      sc(ia) = 0.0
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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
        sc(ia) = sc(ia) + an(ix,jy,kz,ia) ! /gt(ix,jy,kz,imapz)
      end do
      end do
      end do
      ENDDO
      
      scion(1) = 0.0
      scion(2) = 0.0
      scion(3) = 0.0
      scion(4) = 0.0
      scion(5) = 0.0
      scion(6) = 0.0
      
      x = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd-jybeg+1
!      if (jyend .ge. jesd) jye = jesd-jybeg+1

      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd-ixbeg+1
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1

      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
!      do jy = 1, ny-1*jd1
!      do ix = 1, nx-1*id1
      DO jy = jbsd,jesd ! 1,ny-1
      DO ix = ibsd,iesd ! 1,nx-1
#endif
       dv = dxx(ix)*dyy(jy)*dzz(kz)
        scion(1) = scion(1) + Max(0.0,an(ix,jy,kz,lscpi))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        scion(2) = scion(2) + Min(0.0,an(ix,jy,kz,lscpi))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        scion(3) = scion(3) + Max(0.0,an(ix,jy,kz,lscni))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        scion(4) = scion(4) + Min(0.0,an(ix,jy,kz,lscni))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        scion(5) = scion(5) + (an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))

        scion(6) = scion(6) + (an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv ! *dz*gz(kz) ! /gt(ix,jy,kz,imapz)
        
       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  &
     &                             + ec*(an(ix,jy,kz,lscpli)-an(ix,jy,kz,lscnli))
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) 
       ENDIF
      
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      END IF
      
      chgnet = chgnet + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      
      ! temporary advection+mixing tendency, assuming iscnet has value from start of time step
      
!      IF ( Abs(chgadvtemp(ix,jy,kz,1)) > Abs(x) ) x = chgadvtemp(ix,jy,kz,1)
      
      IF ( ichgtndadvsn > 1 ) THEN
        tmp = chgadvtemp(ix,jy,kz,1) ! (t0(ix,jy,kz) - elec(iscnet)%flt3d(ix,jy,kz) )
        IF ( (tmp > 0. .and. elec(iscnet)%flt3d(ix,jy,kz) >= 0. )  .or.  &
             (tmp < 0. .and. elec(iscnet)%flt3d(ix,jy,kz) <= 0. )  .or.  &
             (tmp > 0. .and. elec(iscnet)%flt3d(ix,jy,kz) + tmp >= 0. )  .or.  &
             (tmp < 0. .and. elec(iscnet)%flt3d(ix,jy,kz) + tmp <= 0. )) THEN
          elec(ichgtndadvsn)%flt3d(ix,jy,kz) = elec(ichgtndadvsn)%flt3d(ix,jy,kz) + chgadvtemp(ix,jy,kz,1)
        ENDIF
     !   elec(ichgtndadv)%flt3d(ix,jy,kz) = elec(ichgtndadv)%flt3d(ix,jy,kz) + (t0(ix,jy,kz) - elec(iscnet)%flt3d(ix,jy,kz) )
      ENDIF
      
      elec(iscnet)%flt3d(ix,jy,kz) = t0(ix,jy,kz)
      
      end do
      end do
      end do
      

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 4"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 4"
#endif

!c      dxdy = dxx(ix)*dyy(jy)
!c      dv = dxdy*dzz(kz)
#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg
       mpitotindp(3)  = chgnet
       n = 3
       DO ia = 1,6
         n = n + 1
!         mpitotindp(2+ia) = scion(ia)
         mpitotindp(n) = scion(ia)
       ENDDO
       
!       IF ( n + lsceq - lscb + 1 <= ntot ) THEN
       IF ( microp(1:3) /= 'TAK' ) THEN
       DO ia = lscb,lsceq
         n = n + 1
!         mpitotindp(9+ia-lscb) = sc(ia)
         mpitotindp(n) = sc(ia)
       ENDDO
       ENDIF
       
       IF ( n .gt. ntot ) THEN
         write(0,*) 'Driver: PROBLEM WITH mpitotindp n,ntot = ',n,ntot
         call commasmpi_abort()
       ENDIF

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
       chgnet = mpitotoutdp(3)  
       n = 3
       DO ia = 1,6
         n = n + 1
         scion(ia) = mpitotoutdp(n)
       ENDDO
!       IF ( n + lsceq - lscb + 1 <= ntot ) THEN
       IF ( microp(1:3) /= 'TAK' ) THEN
       DO ia = lscb,lsceq
         n = n + 1
         sc(ia) = mpitotoutdp(n)
       ENDDO
       ENDIF
      ENDIF
#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')     &
     &  'premic pos/neg/net charge (C):',  &
     &      chgpos,chgneg,(chgpos+chgneg)
!     :      chgpos*dx*dy*dz,chgneg*dx*dy*dz,(chgpos+chgneg)*dx*dy*dz
!      write(iunit,'(a,1(2x,1pe15.8))') 'chgnet = ',chgnet*dx*dy*dz
      write(iunit,'(a,1(2x,1pe15.8))') 'chgnet = ',chgnet

      dv = ec ! dx*dy*dz*ec

!c      write(iunit,'(a,6(2x,1pe15.8))') 
!c     : 'scion1,2,3,4,5,6:',dv*scion(1),
!c     :  dv*scion(2),dv*scion(3),dv*scion(4),dv*scion(5),dv*scion(6)

      DO ia = 1,6
        write(iunit,'(a,i1,3(2x,1pe15.8))') 'scion',ia,dv*scion(ia),  &
     &     dv*(scion(ia) - scion2(ia)), dv*scion(ia) - dv*scion2(ia)
      ENDDO
     
      dv = 1.0 ! dx*dy*dz
      
!      IF ( lscb < 100 ) THEN
      IF ( microp(1:3) /= 'TAK' ) THEN
      DO ia = lscb,lsceq
        write(iunit,'(i3,3(2x,1pe15.8))') ia, sc(ia)*dv,  &
     &   (sc(ia) - sc2(ia))*dv, dv*sc(ia) - dv*sc2(ia)
      ENDDO
      ENDIF
      
      ENDIF ! my_rank


       RETURN
       END
       
! ####################################################################


       subroutine ELEC_DRIVE3(gd,                                   &
     &                     pinit, ab,                                &  ! base state
     &                         p2,                                   &  ! pi'
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
     &                        tstat, lstt, bcx, bcy, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi )

      USE FILE_MODULE
      USE VIS5D_MODULE
      USE GRID_MODULE
      USE MICRO_MODULE
      USE ELEC_MODULE
      USE ELEC_STUFF
      USE INDEX_MODULE
      USE CPUTIME_MODULE
!      USE TRMM_MODULE
      USE COMMASMPI_MODULE

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
#endif

      integer nx,ny,nz,na,nor,nba,nv,ia,il

      TYPE(GRID)         :: gd 

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer :: tstat
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

      
      integer iunit
      integer, parameter    :: ng1 = 1

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)

      logical io_flag

      integer ntmul
      real dtp1

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt
      parameter ( nv = 3 )

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)
      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
      integer nx1,ny1,nz1

! local variables
      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n
      integer ix,jy,kz
      real x,tx
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      logical :: debug_mpi

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
!
! ###########################################################################
!


!c  
!c  Count up total net positive (chgpos) and net negative charge (chgneg)
!c


      chgneg = 0.0d0
      chgpos = 0.0d0
      
      t0(:,:,:) = 0.0
      
! C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
! C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
! C$OMP+ REDUCTION (+ : chgpos,chgneg)
      
      DO ia = lscb,lsceq
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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO

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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))

       dv = dxx(ix)*dyy(jy)*dzz(kz)
        
       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  &
     &                             + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) 
       ENDIF
      
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      END IF
      
      
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
      ENDIF
#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))') &
     &    'postmic pos/neg/net charge (C):', chgpos, chgneg, (chgpos+chgneg)
      ENDIF



       RETURN
       END

! ####################################################################

#define ICE10
#define SWM

       subroutine ELEC_DRIVE_LTGLOOP(gd,                             &
     &                     pinit, ab,db,                             &  ! base state
     &                         p2, tt0,                              &  ! pi'
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
!     &                        tstat, lstt, bcx, bcy, tstop,         &
     &                        tstat, lstt, bcx, bcy, bcx1, bcy1, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,   &
     &                        xfall,kcldtop,icldtop,zcldtop,         &
     &                        idelay, idelaycg)

      USE FILE_MODULE
      USE VIS5D_MODULE
      USE GRID_MODULE
      USE GRIDPARAM_MODULE, only: dx_stretch,dy_stretch
      USE MICRO_MODULE
      USE ELEC_MODULE
      USE ELEC_STUFF
      USE INDEX_MODULE
      USE CPUTIME_MODULE
      USE TRMM_MODULE, only : ibsd,jbsd,iesd,jesd
      USE COMMASMPI_MODULE
      USE TRAJ_MODULE, only: elecrates,numflashtraj

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
#endif

      integer nx,ny,nz,na,nor,nba,nv,ia,il

      TYPE(GRID)         :: gd 

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer :: tstat
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

      
      integer iunit
      integer, parameter    :: ng1 = 1

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real tt0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! air temperature (kelvin)
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)
      real xfall(nx,ny,na)
      integer kcldtop(nx,ny)
      integer icldtop
      real zcldtop
      integer idelay, idelaycg

      logical io_flag

      integer ntmul
      real dtp1

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt
      parameter ( nv = 3 )

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)
      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
      integer nx1,ny1,nz1

! local variables
      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n,count
      integer ix,jy,kz
      real x,tx, tmp
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)
      
      integer :: numlgt, iliter, lgtstp
      integer :: loccur

      logical lhelec, lnsave, ldumpsp
      integer iwrite,iwritecg

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      logical :: debug_mpi

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      real, parameter ::  alnox = 0.34e21, blnox = 1.30e16  ! parameters from Wang et al. 1998 (JGR)
      real, parameter ::  avogadro = 6.0221415e23           ! Avogadro's number
      real            ::  air_molmass  ! molecular mass of air
      double precision :: totnox,nox,totlen
!
! ###########################################################################
!

       dtp = dtp1
!
!  Add up the charge that fell to ground 
!

        xfallp = 0.0
        xfalln = 0.0
        xfalltp = 0.0
        xfalltn = 0.0

! CVD$L SKIP
      DO ia=lscb,lsceq


#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd-jybeg+1
!      if (jyend .ge. jesd) jye = jesd-jybeg+1

      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd-ixbeg+1
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1

      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      DO jy = jyb,jye
      DO ix = ixb,ixe
#else
!CVD$L SKIP
!      DO jy=1,ny-1
! CVD$L SKIP
!      DO ix=1,nx-1
      DO jy = jbsd,jesd ! 1,ny-1
      DO ix = ibsd,iesd ! 1,nx-1
#endif
        dxdy = dxx(ix)*dyy(jy)
        dv = dxdy*dzz(1)
!#ifndef SWM
!        xfalltot(ix,jy,ia) = xfalltot(ix,jy,ia) + xfall(ix,jy,ia)
!        xfalltp = xfalltp + Max(xfalltot(ix,jy,ia),0.0)*dv
!        xfalltn = xfalltn + Min(xfalltot(ix,jy,ia),0.0)*dv
!#endif
        xfallp = xfallp + Max(xfall(ix,jy,ia),0.0)*dv
        xfalln = xfalln + Min(xfall(ix,jy,ia),0.0)*dv
      ENDDO
      ENDDO
      ENDDO


#ifdef MPI

       mpitotindp(1) = xfallp
       mpitotindp(2) = xfalln

       n = 2

!      CALL MPI_Allreduce(mpitotin, mpitotout, ntot, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
      CALL MPI_Reduce(mpitotindp, mpitotoutdp, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       IF ( my_rank == 0 ) THEN
        xfallp = mpitotoutdp(1) 
        xfalln = mpitotoutdp(2) 
       ENDIF


#endif

       IF ( my_rank == 0 ) THEN

      write(iunit,'(a,3(2x,1pe15.8))')'Charge fallout this step (pos,neg,tot): ', &
     &   xfallp,xfalln,(xfallp+xfalln)

       ENDIF


!c
!c  lightning iterations per time step
!c
      numlgt = 0
!      t1(:,:,:) = 0.0


       IF ( allocated( elecrates ) ) THEN

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
          elecrates(ix,jy,kz,1) = 0.0
          elecrates(ix,jy,kz,2) = 0.0
          elecrates(ix,jy,kz,3) = 0.0
          elecrates(ix,jy,kz,4) = 0.0
         ENDDO
         ENDDO
         ENDDO
       
       ENDIF


!#ifndef SWM
      do iliter = 1,nliter+1

      IF ( my_rank == 0 ) THEN
        IF ( iliter .le. nliter ) write(iunit,*) 'ILITER',iliter
      ENDIF

!
!      if ( ipotslv .ge. 1 ) then 
!c
!c  set boundary conditions
!c
!c
!c  ezfairo (V/m)
!c
!c
!c constants:  efa1 1/m
!c             efa2 1/m
!c
!c
!c
!c
!c
!c  set the forcing term:  t0 (net charge density)
!c  Also count up total net positive (chgpos) and net negative charge (chgneg)
!c
      chgneg = 0.0d0
      chgpos = 0.0d0
      chgnegcld = 0.0d0
      chgposcld = 0.0d0
      cnapos = 0.0
      cnaneg = 0.0

      CALL cld_cpu('EFIELD')

      t0(:,:,:) = 0.0

!c !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!c !$OMP+ REDUCTION (+ : chgpos,chgneg)
      
      DO ia = lscb,lsceq

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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO


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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
       
       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  &
     &                             + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) 
       ENDIF
        dv = dxx(ix)*dyy(jy)*dzz(kz)

#ifdef MPI
      IF ( jybeg-1+jy .ge. jbsd .and. jybeg-1+jy .le. jesd .and.   &
     &     ixbeg-1+ix .ge. ibsd .and. ixbeg-1+ix .le. iesd ) THEN 
#else
      IF ( jy .ge. jbsd .and. jy .le. jesd .and.       &
     &     ix .ge. ibsd .and. ix .le. iesd ) THEN
#endif
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv
        IF ( kz .le. kcldtop(ix,jy) ) chgposcld = chgposcld + t0(ix,jy,kz)*dv
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv
        IF ( kz .le. kcldtop(ix,jy) ) chgnegcld = chgnegcld + t0(ix,jy,kz)*dv
      END IF
      ENDIF
      
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg
       mpitotindp(3)  = chgposcld
       mpitotindp(4)  = chgnegcld

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
       chgposcld = mpitotoutdp(3)
       chgnegcld = mpitotoutdp(4)
      ENDIF
#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))') 'net1 pos/neg/net charge (C):',  &
     &      chgpos,chgneg,(chgpos+chgneg)
      write(iunit,'(a,3(2x,1pe15.8))') 'storm net1 pos/neg/net charge (C):',  &
     &      chgposcld,chgnegcld,(chgposcld+chgnegcld)
      ENDIF
!      write(iunit,'(a,2(2x,1pe15.8))') 
!     :     'net1 positive/negative corona (C):',
!     :      cnapos,cnaneg
      


       IF ( imud .eq. 1 ) THEN
       
        
!c regular grid:
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

!         nbe = (nnxs(1) - nx - 1)/2
!         nbw = nnxs(1) - nx - 1 - nbe
!         nbn = (nnys(1) - ny - 1)/2
!         nbs = nnys(1) - ny - 1 - nbn

#endif


         IF ( ny .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         IF ( bcx .eq. 2 ) THEN ; nbw = -1 ; nbe = 0 ; ENDIF
         IF ( bcy .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         nbz2 = nnzs(1) - nz
         nzbt2 = nz + nbz2
         nzb2 = nzbt2 + 1
         nxb = nnxs(1)
         nyb = nnys(1)
      
      IF ( time_real .eq. dt  .and. iliter .eq. 1) THEN
        k = 10
      ELSE
        k = 1
      ENDIF  
      

#if defined ( BOXMG ) && defined ( MPI )

      CALL cld_cpu('BOXMG')

!         write(iunit,*) 'driver2: call putf'
         
!c         IF ( my_rank == 0 ) THEN
!c           DO k = 1,NLzdg
!c             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!c           ENDDO
!c         ENDIF
         
!         write(iunit,*) 'call putf, my_rank,nbw, etc. = ',my_rank,nbw,nbe,nbs,nbn,NLxdg,NLydg,nlzdg,ngxdg,ngydg,ngzdg
!c         CALL PUTF( SOdg, QFdg, QFdg, Qdg, t0, nx,ny,nz,nor,     
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,  &
     &               nbw,nbe,nbs,nbn,                            &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 0          &
     &              )
!         write(iunit,*) 'done putf, my_rank,nbw, etc. = ',my_rank,nbw,nbe,nbs,nbn,NLxdg,NLydg,nlzdg,ngxdg,ngydg,ngzdg


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg 2'

       
       i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)
       
       IF ( firstsolve == 0 ) THEN
         firstsolve = 1
         BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 10
         BMG_rPARMSdg(id_BMG3_STOP_TOL) = 1.e-5
       ELSEIF ( iliter .gt. 1 ) THEN
!         BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 4
!         BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 2
!         BMG_iPARMSdg(id_BMG3_NRELAX_DOWN )   = 2
       ENDIF


         CALL commas_SymStd_SOLVE_boxmgdg(                                        &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg, &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                      &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RES)), NFdg, NCbmgdg,    &
     &             SOdg, NSOdg,                                                   &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SOR)), NSORdg,                      &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CI)), NCIdg,                        &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iG)), NOGdg, NOGcdg,                &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                               &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                               &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                               &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                          &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSG)), NMSGidg,                  &
     &             pMSGdg, pMSGSOdg,                                         &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,         &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUF)), NMSGrdg, &
     &             my_comm                                        &
     &             )

       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = bmg_tol

         BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 1
         BMG_iPARMSdg(id_BMG3_NRELAX_DOWN ) = 1

!         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dt2 = MPI_Wtime()

!       For restart, this is the first time setup is turned off
         IF ( BMG_iPARMSdg(id_BMG3_SETUPdg) < 2 ) BMG_iPARMSdg(id_BMG3_SETUPdg) = 2

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for boxmg = ', dT2 - dT1
         ENDIF

!c      IF ( my_rank .eq. 0 ) THEN
!c          CALL PRINTQ( SOdg, QFdg, Qdg,                       &
!c     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, &  
!c     &               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    &
!c     &              )
!c      ENDIF
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,                &
     &               gxt,gyt,gzt,                             &
     &               nbw,nbe,nbs,nbn,iunit,                   &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz          &
     &              )

      CALL cld_cpu('BOXMG')


#else

      IF ( number_of_processes .eq. 1 ) THEN
!c      CALL cld_cpu('EFIELD')
#ifdef MPI
      kzb = 1
      kze = k
      if (kzbeg .gt. k) kze = -1  !!mpidebug: check this ... don't need if k>10

      DO kz=kzb,kze
#else
      DO kz=1,k 
#endif
      ltmp = .false.
      IF ( kz .eq. k ) ltmp = .true.
      CALL mud(nx,ny,nz,nor,          &
     &   dx,dy,dz,gxt,gyt,gzt,        &
     &   nnxs(1),nnys(1),nnzs(1),     &
     &   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1), &
     &   nbw,nbe,nbs,nbn,             &
     &   llworks(1),istretch,         &
     &   t0,ibg,                      &
     &   elec, iestag, iunit, ltmp, bcx, bcy)
      ENDDO
!c      CALL cld_cpu('EFIELD')
      
      ENDIF  ! number_of_processes

#endif  ! boxmg
      
      ENDIF ! (imud .eq. 1)
       


      IF ( number_of_processes .eq. 1 ) THEN
      IF ( iestag .eq. 0 ) THEN
!c      CALL cld_cpu('EFIELD')
      call efield                      &
     &  (nx,ny,nz                      &
     &  ,id1,jd1,kd1,istag,jstag,kstag &
     &  ,dx,dy,dz,elec                 &
     &  ,gxt(1,3),gyt(1,3),gzt(1,3)    &
     &  ,iunit)
      ENDIF
      ENDIF

      CALL cld_cpu('EFIELD')


      emax = 0.0

      energyold = energynew
      energynew = 0.0

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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
! !$DOACROSS LOCAL(kz,jy,ix), REDUCTION(emax,energynew)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : energynew), reduction(max : emax )
      do kz = 1,nz-kd1
      do jy = 1,ny-jd1
      do ix = 1,nx-id1
#endif
      dv = dxx(ix)*dyy(jy)*dzz(kz)

!c      emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
!c Old version of calculating energy:
!c      energynew = energynew + (( elec(ix,jy,kz,iemag) )**2)/
!c     :                         gt(ix,jy,kz,imapz) 
       
         emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
       IF (energymethod == 1 ) THEN
         energynew = energynew + 0.5*elec(ipot)%flt3d(ix,jy,kz)*t0(ix,jy,kz)*dv
       ELSEIF (energymethod == 2 ) THEN
         energynew = energynew + 0.5*elec(iemag)%flt3d(ix,jy,kz)**2*dv
       ENDIF

      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = energynew

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       energynew = mpitotoutdp(1)  
!      ENDIF

! find global max E
       mpitotindp(1)  = emax

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)

       emax = mpitotoutdp(1)  

#endif
      

       IF (energymethod == 2 ) THEN
         energynew = 0.5*eperao*energynew
       ENDIF
      
      IF ( time .gt. ntstart .or. iliter .gt. 1  ) THEN
       IF ( iliter .gt. 1 .and. energynew - energyold .gt. 0.0 ) THEN
         write(iunit,*) 'WARNING: Lightning increased the total energy!',time,iliter
        IF ( my_rank == 0 ) THEN
          write(16,*) '999999 rejected for energy violation'
        ENDIF
!         write(0,*) 'WARNING: Lightning increased the total energy!',
!     :     time,iliter
#ifdef MPI
!c        IF ( number_of_processes .gt. 1 ) THEN
!c          CALL commasmpi_abort()
!c          STOP
!c        ENDIF
#endif
      IF ( Abs(icgyn) .eq. 1 ) then

        IF ( loccur .eq. 3 ) then
          numcgp = numcgp - 1
        ELSEIF ( loccur .eq. 2 ) then
          numcgn = numcgn - 1
        ENDIF
      
      ELSE
        numic = numic - 1
      ENDIF
      
        numlgt = numlgt - 1


       IF ( my_rank == 0 ) THEN
       IF ( energyold .gt. 0.0 ) THEN
       energychange = 100.0*(energynew - energyold)/energyold
       ELSE
       energychange = 0.0
       ENDIF
       write(iunit,'(a,3(1x,1pe13.5),1x,0pf7.2,a)')   &
     &       'Energy info: old,new,difference: ',     &
     &                 energyold, energynew, energynew - energyold,energychange,' %'
       
       ENDIF


!c Take the channels out
!c
          IF ( ixst .gt. 0 .and. jyst .gt. 0 .and. kzst .gt. 0 ) THEN
          elec(ieinit)%flt3d(ixst,jyst,kzst) = elec(ieinit)%flt3d(ixst,jyst,kzst) - 1

          IF ( allocated(flshi) ) THEN
            flshi(ixst,jyst,kzst) = flshi(ixst,jyst,kzst) - 1
          ENDIF
            IF ( iflshfod > 1 ) THEN
               xfalltot(ixst,jyst,iflshfod) = xfalltot(ixst,jyst,iflshfod) - 1  ! remove 1 from flshfod if removed flash started in the column.
            ENDIF
            IF ( iflshfodic > 1 .and. loccur == 1 ) THEN
               xfalltot(ixst,jyst,iflshfodic) = xfalltot(ixst,jyst,iflshfodic) - 1  ! adds 1 to flshfod if IC flash starts in the column.
            ENDIF
            IF ( iflshfodcgn > 1 .and. loccur == 2 ) THEN
               xfalltot(ixst,jyst,iflshfodcgn) = xfalltot(ixst,jyst,iflshfodcgn) - 1  ! adds 1 to flshfod if CGN flash starts in the column.
            ENDIF
            IF ( iflshfodcgp > 1 .and. loccur == 3 ) THEN
               xfalltot(ixst,jyst,iflshfodcgp) = xfalltot(ixst,jyst,iflshfodcgp) - 1  ! adds 1 to flshfod if CGP flash starts in the column.
            ENDIF
          
          ENDIF


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

           DO kz = kzb,kze
            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
           DO kz=1,nz-1
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
               elec(ieflshn)%flt3d(ix,jy,kz) = elec(ieflshn)%flt3d(ix,jy,kz) - Min(0.0,t2(ix,jy,kz))
               elec(ieflshp)%flt3d(ix,jy,kz) = elec(ieflshp)%flt3d(ix,jy,kz) - Max(0.0,t2(ix,jy,kz))

               t2(ix,jy,-1) =  t2(ix,jy,-1) + Abs(t2(ix,jy,kz) )

             IF ( allocated(flshp) ) THEN
               flshp(ix,jy,kz) = flshp(ix,jy,kz) - Max(0.0,t2(ix,jy,kz))
             ENDIF

             IF ( allocated(flshn) ) THEN
               flshn(ix,jy,kz) = flshn(ix,jy,kz) - Min(0.0,t2(ix,jy,kz))
             ENDIF

              lgtth(kzbeg-1+kz,chans) = lgtth(kzbeg-1+kz,chans) - Abs( t2(ix,jy,kz) )
              lgtth(kzbeg-1+kz,chansp) = lgtth(kzbeg-1+kz,chansp) - Max( 0.0, t2(ix,jy,kz) )
              lgtth(kzbeg-1+kz,chansn) = lgtth(kzbeg-1+kz,chansn) - Abs( Min( 0.0, t2(ix,jy,kz) ) )
             
             IF ( loccur .eq. 1 ) THEN
               lgtth(kzbeg-1+kz,icchans)  = lgtth(kzbeg-1+kz,icchans) - Abs( t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,icchansp) = lgtth(kzbeg-1+kz,icchansp) - Max( 0.0, t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,icchansn) = lgtth(kzbeg-1+kz,icchansn) - Abs( Min( 0.0, t2(ix,jy,kz) ) )
             ELSEIF ( loccur .ge. 2 ) THEN
               lgtth(kzbeg-1+kz,cgchans)  = lgtth(kzbeg-1+kz,cgchans) - Abs( t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,cgchansp) = lgtth(kzbeg-1+kz,cgchansp) - Max( 0.0, t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,cgchansn) = lgtth(kzbeg-1+kz,cgchansn) - Abs( Min( 0.0, t2(ix,jy,kz) ) )
             ENDIF

             ENDDO
            ENDDO
           ENDDO

           IF ( iflshfed > 0 .or. iflshsrc > 0 ) THEN
#ifdef MPI
            jyb = 1
            jye = jtile
            if (jyend .eq. nyend) jye = jyend-jybeg

            ixb = 1
            ixe = itile
            if (ixend .eq. nxend) ixe = ixend-ixbeg

            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
              IF ( iflshfed > 0 ) THEN 
                xfalltot(ix,jy,iflshfed) = xfalltot(ix,jy,iflshfed) - Min( 1.0, t2(ix,jy,-1) )
              ENDIF
              
              IF ( iflshsrc > 0  ) THEN
                xfalltot(ix,jy,iflshsrc) = xfalltot(ix,jy,iflshsrc) - t2(ix,jy,-1)
              ENDIF
              t2(ix,jy,-1) = 0.0
             ENDDO
            ENDDO
            ENDIF



       IF ( kzst .ne. 0 ) THEN ! for mpi, only the tile where lightning starts has kzst .ne. 0

         lgtth(kzst,cinit) = lgtth(kzst,cinit) - 1
         IF ( loccur .eq. 1 ) lgtth(kzst,icinit)  = lgtth(kzst,icinit) - 1
         IF ( loccur .eq. 2 ) lgtth(kzst,cgninit) = lgtth(kzst,cgninit) - 1
         IF ( loccur .eq. 3 ) lgtth(kzst,cgpinit) = lgtth(kzst,cgpinit) - 1
         IF ( loccur .eq. 2 .or. loccur .eq. 3 ) lgtth(kzst,cginit) = lgtth(kzst,cginit) - 1
         
       ENDIF
!c
!c take back the lightning charge
!c
        chgneg = 0.0d0
        chgpos = 0.0d0


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

        do kz = kzb,kze
        do jy = jyb,jye
        do ix = ixb,ixe
#else
        do kz = 1,nz-kd1
        do jy = 1,ny-jd1
        do ix = 1,nx-id1
#endif
        dv = dxx(ix)*dyy(jy)*dzz(kz)
        IF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
          an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni)+t3(ix,jy,kz)/ec
        ELSE 
          an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi)-t3(ix,jy,kz)/ec
        ENDIF

         t0(ix,jy,kz) = t0(ix,jy,kz) - t3(ix,jy,kz)

        IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
          chgpos = chgpos + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        ELSE
          chgneg = chgneg + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        END IF
        end do
        end do
        end do

         write(iunit,*) 'Remove lightning charge and recalculate potential and energy'

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
     &                0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
      ENDIF
#endif

      write(iunit,'(a,3(2x,1pe15.8))') 'net1 pos/neg/net charge (C):',chgpos,chgneg,(chgpos+chgneg)
!      write(iunit,'(a,2(2x,1pe15.8))') 
!     :  'net1 positive/negative corona (C):',
!     :      cnapos,cnaneg

      CALL cld_cpu('EFIELD')

#if defined (BOXMG) && defined (MPI)

      CALL cld_cpu('BOXMG')

!         write(iunit,*) 'driver3: call putf'
         
!c         IF ( my_rank == 0 ) THEN
!c           DO k = 1,NLzdg
!c             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!c           ENDDO
!c         ENDIF
         
!c         CALL PUTF( SOdg, QFdg, QFdg, Qdg, t0, nx,ny,nz,nor,     
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,     &
     &               nbw,nbe,nbs,nbn,                                                   &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,                          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 0                                 &
     &              )
!         write(iunit,*) 'done putf'


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!c         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg 3'

!         i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)
!       
!         BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 4
!         BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 2
!         BMG_iPARMSdg(id_BMG3_NRELAX_DOWN )   = 2

         CALL commas_SymStd_SOLVE_boxmgdg(                                   &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg, &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                   &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg, &
     &             SOdg, NSOdg,                                              &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                &   
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                   &  
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,            &   
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                          &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                          &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                          &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                          &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,               &   
     &             pMSGdg, pMSGSOdg,                                         &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,        & 
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg,&
     &             my_comm                                        &
     &             )



!         BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i
!         BMG_rPARMSdg(id_BMG3_STOP_TOL) = bmg_tol

!         BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 1
!         BMG_iPARMSdg(id_BMG3_NRELAX_DOWN ) = 1

!c         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'


         CALL MPI_Barrier(my_comm, mpi_error_code)
         dt2 = MPI_Wtime()

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for boxmg = ', dT2 - dT1
         ENDIF

!c      IF ( my_rank .eq. 0 ) THEN
!c          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!c     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!c     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!c     :              )
!c      ENDIF

!c
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,     &
     &               gxt,gyt,gzt,                  &
     &               nbw,nbe,nbs,nbn,iunit, &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz  &
     &              )

      CALL cld_cpu('BOXMG')

#else

      IF ( number_of_processes .eq. 1 ) THEN

      CALL mud(nx,ny,nz,nor,          &
     &   dx,dy,dz,gxt,gyt,gzt,         &
     &   nnxs(1),nnys(1),nnzs(1),      &
     &   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1), &
     &   nbw,nbe,nbs,nbn,               &
     &   llworks(1),istretch,           &
     &   t0,ibg,                        &
     &   elec, iestag, iunit, .true., bcx, bcy)
     
       ENDIF
#endif

      CALL cld_cpu('EFIELD')

      energynew = 0.0
      emax = 0.0


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

        do kz = kzb,kze
        do jy = jyb,jye
        do ix = ixb,ixe
#else
! C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
! C$OMP+ REDUCTION (+ : energynew), reduction(max : emax )
        do kz = 1,nz-kd1
        do jy = 1,ny-jd1
        do ix = 1,nx-id1
#endif
        dv = dxx(ix)*dyy(jy)*dzz(kz)

         emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
       IF (energymethod == 1 ) THEN
         energynew = energynew + 0.5*elec(ipot)%flt3d(ix,jy,kz)*t0(ix,jy,kz)*dv
       ELSEIF (energymethod == 2 ) THEN
         energynew = energynew + 0.5*elec(iemag)%flt3d(ix,jy,kz)**2*dv
       ENDIF

        end do
        end do
        end do

#ifdef MPI
! find global integrated electrical energy
       mpitotindp(1)  = energynew

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       energynew = mpitotoutdp(1)  

! find global max E
       mpitotindp(1)  = emax

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)

       emax = mpitotoutdp(1)  

#endif

       ELSEIF ( iliter > 1 ) THEN ! flash was fine.  Add LNOX if turned on
         IF  ( lnox > 1 .and. ilight == 2 ) THEN
          air_molmass = 28.96e-3
          totnox = 0.0d0
          totlen = 0.0d0
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

           DO kz = kzb,kze
            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
           DO kz=1,nz-1
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
              IF ( t1(ix,jy,kz) > 0.0 ) THEN
              dv = dxx(ix)*dyy(jy)*dzz(kz)

! Formula from Wang et al 1998 (JGR) as in Barthe et al. 2007 (JGR)
! t1 is the channel length in the grid box
! NO mass 30.0061e-3 kg/mol; 
! mol mass of dry air is 28.96 10^(-3) kg
! This is not quite correct yet, since air density should be for moist air rather than dry
! but lower density for air and lower molecular mass should compensate errors somewhat.
              
      !        an(ix,jy,kz,lnox) = an(ix,jy,kz,lnox) + &
      !           t1(ix,jy,kz)*air_molmass*(alnox + blnox*pn(ix,jy,kz))/(avogadro*dv*dn(ix,jy,kz))
              nox = t1(ix,jy,kz)*(alnox + blnox*(pn(ix,jy,kz)+pb(kz)))/(dn(ix,jy,kz)*avogadro)
              an(ix,jy,kz,lnox) = an(ix,jy,kz,lnox) + air_molmass*nox/(dv)
              totnox = totnox + nox
              totlen = totlen + t1(ix,jy,kz)
              ENDIF
             ENDDO
            ENDDO
           ENDDO

#ifdef MPI
! find global integrated electrical energy
       mpitotindp(1)  = totnox
       mpitotindp(2)  = totlen

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       IF ( my_rank == 0 ) THEN
         totnox = mpitotoutdp(1)
         totlen = mpitotoutdp(2)
       ENDIF
#endif
       IF ( my_rank == 0 ) THEN
       write(iunit,'(a,2(1x,1pe13.5))')    &
     &       'Moles LNOX and length for this flash: ', totnox,totlen
       ENDIF
           
         ENDIF
       ENDIF !( iliter .gt. 1 .and. energynew - energyold .gt. 0.0 ) 


       IF ( energyold .gt. 0.0 ) THEN
       energychange = 100.0*(energynew - energyold)/energyold
       ELSE
       energychange = 0.0
       ENDIF
       
       IF ( my_rank == 0 ) THEN
       write(iunit,'(a,3(1x,1pe13.5),1x,0pf7.2,a)')    &
     &       'Energy info: old,new,difference: ',      &
     &                 energyold, energynew, energynew - energyold,energychange, ' %'
       ENDIF

      END IF ! ( time .gt. ntstart .or. iliter .gt. 1  )

!      IF ( my_rank == 0 ) THEN
!      write(3,*) 'E-MAX=',time,emax
!      ENDIF

#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 7"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 7"
#endif
!c
!c balloon sounding
!c
!c

! set starting values for field change arrays
       IF ( allocated( elecrates ) .and. iliter .eq. 1 ) THEN

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
          elecrates(ix,jy,kz,1) = elec(iex)%flt3d(ix,jy,kz)
          elecrates(ix,jy,kz,2) = elec(iey)%flt3d(ix,jy,kz)
          elecrates(ix,jy,kz,3) = elec(iez)%flt3d(ix,jy,kz)
          elecrates(ix,jy,kz,4) = elec(ipot)%flt3d(ix,jy,kz)
         ENDDO
         ENDDO
         ENDDO
       
       ENDIF


      IF ( number_of_processes .eq. 1 ) THEN
      IF ( iliter .eq. 1 ) THEN

      CALL balloonf                                 &
     &  (ntmul,time,tstop,nstep,nx,ny,nz,na,iunit,  &
     &   nor,istag,jstag,kstag,                     &
     &   elec,                                      &
     &   dtp1,dx,dy,dz,                             &
     &   dxx,dyy,dzz,                               &
     &   gxt,gyt,gzt,                               &
     &   t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,             &
     &   an,pb,pn,                                  &
     &   u,v,w)                                     

      
      END IF !  ( iliter .eq. 1 ) 
      ENDIF


#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 7a"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 7a"
#endif

!c check if this is the extra loop after max number of flashes
      IF ( iliter .gt. nliter ) EXIT ! GOTO 3999
!c
!c  if ilight > 0 
!c
      lgtstp = 0
      loccur = 0


      IF ( ilight .gt. 0 ) THEN !{
      
      icgyn = 0
!      isa = 0

!c lightning grid:
#ifdef BOXMG

         nbw = nbwmg(2)
         nbe = nbemg(2)
         nbs = nbsmg(2)
         nbn = nbnmg(2)

#else

         nbw = (nnxs(2) - nxslm - 1)/2
         nbe = nnxs(2) - nxslm - 1 - nbw
         nbs = (nnys(2) - nyslm - 1)/2
         nbn = nnys(2) - nyslm - 1 - nbs

!         nbe = (nnxs(2) - nxslm - 1)/2
!         nbw = nnxs(2) - nxslm - 1 - nbe
!         nbn = (nnys(2) - nyslm - 1)/2
!         nbs = nnys(2) - nyslm - 1 - nbn
#endif

         IF ( ny .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         IF ( bcx .eq. 2 ) THEN ; nbw = -1 ; nbe = 0 ; ENDIF
         IF ( bcy .eq. 2 ) THEN ; nbs = -1 ; nbn = 0 ; ENDIF

         nbz2 = nnzs(2) - nzslm
         nzbt2 = nzslm + nbz2
         nzb2 = nzbt2 + 1
        
!c        IF ( dslightz .lt. dzmax ) THEN
!c          icldtop = icldtop * Int(dzmax/dslightz + 0.001) 
!c        ENDIF


!#ifdef SWM
#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 7b"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 7b"
#endif
        
#ifndef BOXMG
      IF ( number_of_processes .eq. 1 ) THEN
#else
      IF ( .true. ) THEN
#endif
      CALL cld_cpu('LIGHTNING')

      IF ( irand == 1 ) THEN
        ixst = ixst0
        jyst = jyst0
        kzst = kzst0
      ENDIF
      
#ifdef USE_SLM
      IF ( ilight == 1 .or. ilight == 2 ) THEN
      call mudlight                                 &
     &  (time_real,time,nstep,iliter,lgtstp,numlgt,loccur     &
     &  ,nor,nx,ny,nz                                &
     &  ,na                                           &
     &  ,id1,jd1,kd1,istag,jstag,kstag,icldtop,zcldtop &
     &  ,icgyn                                         &
     &  ,dtp1,dx_stretch,dy_stretch,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9    &
     &  ,elec,an,gx,gy,gz,gxt,gyt,gzt                  &
     &  ,dxx,dyy,dzz                                   &
     &  ,ixst,jyst,kzst                      &
     &  ,numic,numcgn,numcgp,outname,nkzmn,            &
     &   nnxs(2),nnys(2),nnzs(2),                      &
     &   nnxs,nnys,nnzs,                                &
     &   iixps(2),jjyqs(2),kkzrs(2),iiexs(2),jjeys(2),kkezs(2), &
     &   nbw,nbe,nbs,nbn,nbz2,                                 &
     &   llworks(2), istretch, dzmax,                           &
     &   nxslm,nyslm,nzslm,  &  ! nb2,nxb2,nyb2,nzbt2,nzb2,      
     &   iunit, bcx, bcy, db) 
      ELSEIF ( ilight == 3 ) THEN
#else
      IF ( ilight == 1 .or. ilight == 2 .or. ilight == 3) THEN
#endif
       call mszlight                                 &
     &  (time_real,time,nstep,iliter,lgtstp,numlgt,loccur     &
     &  ,nor,nx,ny,nz                                &
     &  ,na                                           &
     &  ,id1,jd1,kd1,istag,jstag,kstag,icldtop,zcldtop &
     &  ,icgyn                                         &
     &  ,dtp1,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9    &
     &  ,tt0,elec,an,gx,gy,gz,gxt,gyt,gzt                  &
     &  ,dxx,dyy,dzz                                   &
     &  ,ixst,jyst,kzst                      &
     &  ,numic,numcgn,numcgp,outname,nkzmn,            &
     &   nnxs(2),nnys(2),nnzs(2),                      &
     &   nnxs,nnys,nnzs,                                &
     &   iixps(2),jjyqs(2),kkzrs(2),iiexs(2),jjeys(2),kkezs(2), &
     &   nbw,nbe,nbs,nbn,nbz2,                                 &
     &   llworks(2), istretch, dzmax,                           &
     &   nxslm,nyslm,nzslm,  &  ! nb2,nxb2,nyb2,nzbt2,nzb2,      
!     &   iunit, dn, pn, bcx, bcy, db) 
     &   iunit, dn, pn, bcx, bcy, bcx1, bcy1, db, pb) 
      ELSE
       write(0,*) 'Bad value for ilight! Must be 1,2, or 3. ilight = ',ilight
       call commasmpi_abort()
      ENDIF

      CALL cld_cpu('LIGHTNING')
      
      ELSE
        lgtstp = 1
      ENDIF

   ! save current seed value
    IF( .not. SET_VARIABLE(gd,'LGTSEED',  iseed)          ) write(6,*) 'ELEC_DRIVER:  Problem setting LGTSEED'

       
!#else

!       loccur = 0
!#endif


       if ( loccur .ge. 1 ) then ! {

!         write(3,*) 'Lightning occurred at step',nstep
!c
!c  add lightning densities to total if lightning happened and vis5d output is on
!c         
!         IF ( ifv5d .and. nstep .ge. iv5dstart - iv5dinterval ) THEN 
         
          IF ( ilight > 1 .and. ixst .gt. 0 .and. jyst .gt. 0 .and. kzst .gt. 0 ) THEN
            elec(ieinit)%flt3d(ixst,jyst,kzst) = elec(ieinit)%flt3d(ixst,jyst,kzst) + 1
            IF ( allocated(flshi) ) THEN
              flshi(ixst,jyst,kzst) = flshi(ixst,jyst,kzst) + 1
            ENDIF
            IF ( iflshfod > 1 ) THEN
               xfalltot(ixst,jyst,iflshfod) = xfalltot(ixst,jyst,iflshfod) + 1  ! adds 1 to flshfod if flash starts in the column.
            ENDIF
            IF ( iflshfodic > 1 .and. loccur == 1 ) THEN
               xfalltot(ixst,jyst,iflshfodic) = xfalltot(ixst,jyst,iflshfodic) + 1  ! adds 1 to flshfod if IC flash starts in the column.
            ENDIF
            IF ( iflshfodcgn > 1 .and. loccur == 2 ) THEN
               xfalltot(ixst,jyst,iflshfodcgn) = xfalltot(ixst,jyst,iflshfodcgn) + 1  ! adds 1 to flshfod if CGN flash starts in the column.
            ENDIF
            IF ( iflshfodcgp > 1 .and. loccur == 3 ) THEN
               xfalltot(ixst,jyst,iflshfodcgp) = xfalltot(ixst,jyst,iflshfodcgp) + 1  ! adds 1 to flshfod if CGP flash starts in the column.
            ENDIF
          ENDIF


#ifdef MPI
            jyb = 1
            jye = jtile
            if (jyend .eq. nyend) jye = jyend-jybeg

            ixb = 1
            ixe = itile
            if (ixend .eq. nxend) ixe = ixend-ixbeg

            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
              t2(ix,jy,-2:0) = 0.0
             ENDDO
            ENDDO


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

           DO kz = kzb,kze
            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
           DO kz=1,nz-1
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
               elec(ieflshn)%flt3d(ix,jy,kz) =   &
     &              elec(ieflshn)%flt3d(ix,jy,kz) + Min(0.0,t2(ix,jy,kz))
               elec(ieflshp)%flt3d(ix,jy,kz) =   &
     &              elec(ieflshp)%flt3d(ix,jy,kz) + Max(0.0,t2(ix,jy,kz))
             
!               t2(ix,jy,-1) = Max( t2(ix,jy,-1), Min(1.0, Abs(t2(ix,jy,kz) ) ) )
               t2(ix,jy,-1) = t2(ix,jy,-1) + Abs(t2(ix,jy,kz) )
               t2(ix,jy,-2) = t2(ix,jy,-2) + Max( 0.0,t2(ix,jy,kz) )
               t2(ix,jy, 0) = t2(ix,jy, 0) + Min( 0.0,t2(ix,jy,kz) )
             IF ( allocated(flshp) ) THEN
               flshp(ix,jy,kz) = flshp(ix,jy,kz) + Max(0.0,t2(ix,jy,kz))
             ENDIF

             IF ( allocated(flshn) ) THEN
               flshn(ix,jy,kz) = flshn(ix,jy,kz) + Min(0.0,t2(ix,jy,kz))
             ENDIF

              lgtth(kzbeg-1+kz,chans) = lgtth(kzbeg-1+kz,chans) + Abs( t2(ix,jy,kz) )
              lgtth(kzbeg-1+kz,chansp) = lgtth(kzbeg-1+kz,chansp) + Max( 0.0, t2(ix,jy,kz) )
              lgtth(kzbeg-1+kz,chansn) = lgtth(kzbeg-1+kz,chansn) + Abs( Min( 0.0, t2(ix,jy,kz) ) )
             
             IF ( loccur .eq. 1 ) THEN
               lgtth(kzbeg-1+kz,icchans)  = lgtth(kzbeg-1+kz,icchans) + Abs( t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,icchansp) = lgtth(kzbeg-1+kz,icchansp) + Max( 0.0, t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,icchansn) = lgtth(kzbeg-1+kz,icchansn) + Abs( Min( 0.0, t2(ix,jy,kz) ) )
             ELSEIF ( loccur .ge. 2 ) THEN
               lgtth(kzbeg-1+kz,cgchans)  = lgtth(kzbeg-1+kz,cgchans) + Abs( t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,cgchansp) = lgtth(kzbeg-1+kz,cgchansp) + Max( 0.0, t2(ix,jy,kz) )
               lgtth(kzbeg-1+kz,cgchansn) = lgtth(kzbeg-1+kz,cgchansn) + Abs( Min( 0.0, t2(ix,jy,kz) ) )
             ENDIF
             
             ENDDO
            ENDDO
           ENDDO
           
!            tmp = 0.0
           IF ( iflshsrc > 0 .or. iflshfed > 0 ) THEN
#ifdef MPI
            jyb = 1
            jye = jtile
            if (jyend .eq. nyend) jye = jyend-jybeg

            ixb = 1
            ixe = itile
            if (ixend .eq. nxend) ixe = ixend-ixbeg

            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
              IF ( iflshfed > 0 ) THEN
                xfalltot(ix,jy,iflshfed) = xfalltot(ix,jy,iflshfed) + Min( 1.0, t2(ix,jy,-1) ) ! adds 1 to flshr if flash goes through the column.
              ENDIF

              IF ( iflshsrc > 0 ) THEN
                xfalltot(ix,jy,iflshsrc) = xfalltot(ix,jy,iflshsrc) + t2(ix,jy,-1)  ! adds all "sources" (channel segments) in the column.
              ENDIF
!              tmp = tmp + xfalltot(ix,jy,5)

              IF ( iflshfedic > 0 .and. loccur == 1 ) THEN
                xfalltot(ix,jy,iflshfedic) = xfalltot(ix,jy,iflshfedic)  + Min( 1.0, t2(ix,jy,-1) )
              ENDIF

              IF ( iflshfedicp > 0 .and. loccur == 1 ) THEN
                xfalltot(ix,jy,iflshfedicp) = xfalltot(ix,jy,iflshfedicp)  + Min( 1.0, t2(ix,jy,-2) )
              ENDIF
              IF ( iflshfedicn > 0 .and. loccur == 1 ) THEN
                xfalltot(ix,jy,iflshfedicn) = xfalltot(ix,jy,iflshfedicn)  + Max( -1.0, t2(ix,jy,0) )
              ENDIF

              IF ( iflshfedcgn > 0 .and. loccur == 2 ) THEN
                xfalltot(ix,jy,iflshfedcgn) = xfalltot(ix,jy,iflshfedcgn)  + Min( 1.0, t2(ix,jy,-1) )
              ENDIF

              IF ( iflshfedcgp > 0 .and. loccur == 3 ) THEN
                xfalltot(ix,jy,iflshfedcgp) = xfalltot(ix,jy,iflshfedcgp)  + Min( 1.0, t2(ix,jy,-1) )
              ENDIF

             ENDDO
            ENDDO
            
           ENDIF
            
! hack for 8-km FED. Assumes that dx=1000m and that nxt and nyt are multiples of 8 !
           IF ( iflshr8km > 0 .and. Abs(dxx(1) - 1000.) < 1.0 .and. Mod(ny,8) == 0 .and. Mod(nx,8) == 0) THEN
!           IF ( iflshr8km > 0  ) THEN

!            write(iunit,*) 'check xfalltot(8km) = '
           
           tmp = 0.0
#ifdef MPI
            jyb = 4
            jye = jtile
            if (jyend .eq. nyend) jye = jyend-jybeg

            ixb = 4
            ixe = itile
            if (ixend .eq. nxend) ixe = ixend-ixbeg

            DO jy = jyb,jye,8
             DO ix = ixb,ixe,8
#else
            DO jy=4,ny-1,8
             DO ix=4,nx-1,8
#endif
              count = 0
              DO j = jy-3,jy+4
               DO i = ix-3,ix+4
                 IF ( t2(i,j,-1) > 0 ) count = 1
               ENDDO
              ENDDO
              
              IF ( count == 1 ) THEN
!              write(iunit,*) 'add to 8km FED',ix,jy
              DO j = jy-3,jy+4
               DO i = ix-3,ix+4
                 xfalltot(i,j,iflshr8km) = xfalltot(i,j,iflshr8km) + count ! adds 1 to flshr if flash goes through the column.
               ENDDO
              ENDDO
                 tmp = tmp + xfalltot(ix,jy,iflshr8km)
              ENDIF
             ENDDO
            ENDDO

!            write(iunit,*) 'total in xfalltot(8km) = ',tmp
         

           ELSEIF ( iflshr8km > 0 ) THEN

#ifdef MPI
            jyb = 1
            jye = jtile
            if (jyend .eq. nyend) jye = jyend-jybeg

            ixb = 1
            ixe = itile
            if (ixend .eq. nxend) ixe = ixend-ixbeg

            DO jy = jyb,jye
             DO ix = ixb,ixe
#else
            DO jy=1,ny-1
             DO ix=1,nx-1
#endif
              xfalltot(ix,jy,iflshr8km) = xfalltot(ix,jy,iflshr8km) + t2(ix,jy,-1)  ! adds total sources
!              tmp = tmp + xfalltot(ix,jy,5)
             ENDDO
            ENDDO
           
           ENDIF
            
            t2(:,:,-1) = 0.0


!         ENDIF ! ( ifv5d )
       
       IF ( kzst .ne. 0 ) THEN ! for mpi, only the tile where lightning starts has kzst .ne. 0

         lgtth(kzst,cinit) = lgtth(kzst,cinit) + 1
         IF ( loccur .eq. 1 ) lgtth(kzst,icinit)  = lgtth(kzst,icinit) + 1
         IF ( loccur .eq. 2 ) lgtth(kzst,cgninit) = lgtth(kzst,cgninit) + 1
         IF ( loccur .eq. 3 ) lgtth(kzst,cgpinit) = lgtth(kzst,cgpinit) + 1
         IF ( loccur .eq. 2 .or. loccur .eq. 3 ) lgtth(kzst,cginit) = lgtth(kzst,cginit) + 1

       ENDIF
       
       end if !}
       
!       if ( loccur .eq. 0 ) then
!       write(3,*) 'Lightning did not occur at step',nstep
!       end if
       
       ENDIF !}

#ifndef SWM 
! MIGHT NEVER USE THIS IN COMMAS
!
!  write some electricity stuff 
!
      lhelec = (nstep/nhelec*nhelec .eq. nstep) .and. (iliter .eq. 1)
      if ( lhelec .or. (loccur .ge. 1 )  ) then ! {

      if ( ( iwrite .ne. 0 .and. loccur .ge. 1 ) .or.               &
     &     ( iwritecg .ne. 0 .and. Abs(icgyn) .eq. 1 ) .or. lhelec  &
     &     ) THEN
       
       IF (loccur .ge. 1) THEN
       idelay = 0
       iwrite = 0
       END IF
       IF (Abs(icgyn) .ge. 1) THEN
       idelaycg = 0
       iwritecg = 0
       END IF
       
       
       icgyn = 0

      times = (nstep-1)*(dtp1)
      write(nmliter,915) iliter
 915  format(i3.3)

        write(rstime,timfmt) int(dtp1*(nstep-1))
        irstime = ihttima

!c       IF ( int(times) .gt. 99999 ) THEN
!c        write(rstime,'(i6.6)') int(times)
!c        irstime = 6
!c      ELSE
!c        write(rstime,'(i5.5)') int(times)
!c        irstime = 5
!c      ENDIF

!c      write(rstime,912) int(times)
!c 912  format(i5.5)

!c      DO jy=1,ny-jd1
!c      DO ix=1,nx-id1
!c        t0(ix,jy,1) = t0(ix,jy,1) + sccna(ix,jy)
!c      END DO
!c      END DO

      fnelec = hdrcfs(1:lhdrcfs)//filehead(1:lfilehead)  &
     & //'.nelec'//'.'//nmliter//'.'//rstime(1:irstime)
        CALL asnctl ('NEWLOCAL', 1, ierr)
        CALL asnfile(fnelec, '-F f77 -N ieee', ierr)
      
      open(unit=9,file=fnelec,form='unformatted')
     
!c       write(9) sccna
       
       
      
      IF ( lhelec .and. loccur .le. 0   ) THEN
        write(9) 0  !  nothing more to save
      ELSEIF ( lhelec .and. loccur .ge. 1  ) THEN
        write(9) 1  ! just save gather of flash locations and charge deposited
        CALL ELEC1SAV(nx,ny,nz,nor,na,9,an,t4,t3,t0,istag,jstag,kstag )
      ELSE
        write(9) 2  ! stand-alone mini-dump
        CALL ELEC2SAV(nx,ny,nz,nor,na,lhelec,9,an,t4,t3,t0,istag,jstag,kstag )
      
      END IF
      

      close(9)

      
      
      end if ! (iwrite .ne. 0)

      end if ! } lhelec 
!c
!c  end ipelec 
!c
!c  end ipot
!c
#endif

      if ( lgtstp .eq. 1 ) EXIT ! go to 3999
      if ( ilight .eq. 0 ) EXIT ! go to 3999
!c
!      end if
!c
      end do ! iliter


! set starting values for field change arrays
       IF ( allocated( elecrates ) ) THEN

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
          elecrates(ix,jy,kz,5) = elec(iex)%flt3d(ix,jy,kz) - elecrates(ix,jy,kz,1) 
          elecrates(ix,jy,kz,6) = elec(iey)%flt3d(ix,jy,kz) - elecrates(ix,jy,kz,2) 
          elecrates(ix,jy,kz,7) = elec(iez)%flt3d(ix,jy,kz) - elecrates(ix,jy,kz,3) 
          elecrates(ix,jy,kz,8) = elec(ipot)%flt3d(ix,jy,kz) - elecrates(ix,jy,kz,4) 
         ENDDO
         ENDDO
         ENDDO
       
       ENDIF


!      chgneg = 0.0d0
!      chgpos = 0.0d0
!      chgnegcld = 0.0d0
!      chgposcld = 0.0d0
!      cnapos = 0.0
!      cnaneg = 0.0
!
!!      CALL cld_cpu('EFIELD')
!
!      t0(:,:,:) = 0.0
!
!!c !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!!c !$OMP+ REDUCTION (+ : chgpos,chgneg)
!      
!      DO ia = lscb,lsceq
!
!#ifdef MPI
!      kzb = 1
!      kze = ktile
!      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1
!
!      jyb = 1
!      jye = jtile
!      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1
!
!      ixb = 1
!      ixe = itile
!      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1
!
!      do kz = kzb,kze
!      do jy = jyb,jye
!      do ix = ixb,ixe
!#else
!      do kz = 1, nz-1*kd1
!      do jy = 1, ny-1*jd1
!      do ix = 1, nx-1*id1
!#endif
!        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
!      end do
!      end do
!      end do
!      ENDDO
!
!
!#ifdef MPI
!      kzb = 1
!      kze = ktile
!      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1
!
!      jyb = 1
!      jye = jtile
!      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1
!
!      ixb = 1
!      ixe = itile
!      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1
!
!      do kz = kzb,kze
!      do jy = jyb,jye
!      do ix = ixb,ixe
!#else
!      do kz = 1, nz-1*kd1
!      do jy = 1, ny-1*jd1
!      do ix = 1, nx-1*id1
!#endif
!        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
!        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
!       
!       IF ( largeion ) THEN
!       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  &
!     &                             + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
!       ELSE
!       t0(ix,jy,kz) = t0(ix,jy,kz) + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) 
!       ENDIF
!        dv = dxx(ix)*dyy(jy)*dzz(kz)
!
!#ifdef MPI
!      IF ( jybeg-1+jy .ge. jbsd .and. jybeg-1+jy .le. jesd .and.   &
!     &     ixbeg-1+ix .ge. ibsd .and. ixbeg-1+ix .le. iesd ) THEN 
!#else
!      IF ( jy .ge. jbsd .and. jy .le. jesd .and.       &
!     &     ix .ge. ibsd .and. ix .le. iesd ) THEN
!#endif
!      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
!        chgpos = chgpos + t0(ix,jy,kz)*dv
!      ELSE
!        chgneg = chgneg + t0(ix,jy,kz)*dv
!      END IF
!      ENDIF
!
!     ! lightning charge tendency, assuming iscnet has value from start of time step
!      IF ( ichgtndlgt > 1 ) THEN
!        elec(ichgtndlgt)%flt3d(ix,jy,kz) = elec(ichgtndlgt)%flt3d(ix,jy,kz) + (t0(ix,jy,kz) - elec(iscnet)%flt3d(ix,jy,kz) )
!      ENDIF
!
!      elec(iscnet)%flt3d(ix,jy,kz) = t0(ix,jy,kz)
!      
!      end do
!      end do
!      end do
!
!#ifdef MPI
!! find global integrated rate max
!       mpitotindp(1)  = chgpos
!       mpitotindp(2)  = chgneg
!
!      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
!
!       
!      IF ( my_rank == 0 ) THEN
!       chgpos = mpitotoutdp(1)  
!       chgneg = mpitotoutdp(2)  
!      ENDIF
!#endif
!      
!      IF ( my_rank == 0 ) THEN
!      write(iunit,'(a,3(2x,1pe15.8))') 'postlight pos/neg/net charge (C):',  &
!     &      chgpos,chgneg,(chgpos+chgneg)
!      ENDIF
!!      write(iunit,'(a,2(2x,1pe15.8))') 
!!     :     'net1 positive/negative corona (C):',
!!     :      cnapos,cnaneg
!      


#ifdef MPI
      if (debug_mpi) write(0,*) "DRIVER: DEBUG 8"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 8"
#endif

!3999  continue

      IF ( my_rank == 0 ) THEN
        write(iunit,'(a,f9.2,a,i3)') 'NUMBER OF FLASHES THIS TIME STEP',time_real,' : ',numlgt
      ENDIF
      
      IF ( allocated(elecrates) ) numflashtraj = numlgt


       RETURN
       END

! ####################################################################


       subroutine ELEC_DRIVE4(gd,                                   &
     &                     pinit, ab,                                &  ! base state
     &                         p2,                                   &  ! pi'
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
     &                        tstat, lstt, bcx, bcy, tstop,         &
     &                        nstep, dxx,dyy,dzz, dn,pn,debug_mpi,    &
     &                        kcldtop)

      USE FILE_MODULE
      USE VIS5D_MODULE
      USE GRID_MODULE
      USE MICRO_MODULE
      USE ELEC_MODULE
      USE ELEC_STUFF
      USE INDEX_MODULE
      USE CPUTIME_MODULE
!      USE TRMM_MODULE
      USE COMMASMPI_MODULE
      USE RUN_ATT_NML, only: thistory

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
#endif

      integer nx,ny,nz,na,nor,nba,nv,ia,il

      TYPE(GRID)         :: gd 

      real               :: u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      TYPE(VARIABLE)     :: elec(neelec)
      TYPE(VARIABLE)     :: cion(2)
      TYPE(VARIABLE)     :: muz(4)
      
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer :: tstat
      logical :: lstt
      
      integer  :: bcx, bcy, bcx1, bcy1
      

      
      integer iunit
      integer, parameter    :: ng1 = 1

!
! external temporary arrays
!
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
!      real vn(nx,ny,nz,nv)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)
      integer kcldtop(nx,ny)

      logical io_flag

      integer ntmul
      real dtp1

      integer nstart, nstop, nstep, nstep0
      integer time, tstop
      real    time_real
      real dt
      parameter ( nv = 3 )

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)
      integer nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      
      real dtp,dx,dy,dz,dv,dvt,dxdy,dv9
      integer nx1,ny1,nz1

! local variables
      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer i,j,k,n
      integer ix,jy,kz
      real x,tx
      double precision :: scx(lc:lhab)
      
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      logical :: debug_mpi
      integer, parameter :: ndebugd = 0

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      double precision, allocatable :: mpitotinth(:,:),mpitotoutth(:,:)
!
! ###########################################################################
!


!c  
!c  Also count up total net positive (chgpos) and net negative charge (chgneg)
!c

      IF ( ndebugd >= 1 ) write(0,*) 'start of elec_drive4'
      chgneg = 0.0d0
      chgpos = 0.0d0
      chgnegcld = 0.0d0
      chgposcld = 0.0d0
      cnapos = 0.0
      cnaneg = 0.0

      t0(:,:,:) = 0.0

!c C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
!c C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!c C$OMP+ REDUCTION (+ : chgpos,chgneg,chgposcld,chgnegcld)
      
      
      DO ia = lscb,lsceq
      sc2(ia) = 0.0d0


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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        dv = dxx(ix)*dyy(jy)*dzz(kz)
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
        sc2(ia) = sc2(ia) + an(ix,jy,kz,ia)*dv ! /gt(ix,jy,kz,imapz)
      end do
      end do
      end do
      ENDDO


#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1

      DO kz = kzb,kze
#else
      DO kz = 1, nz
#endif
        chgposz(kz) = 0.0d0
        chgnegz(kz) = 0.0d0
        chgnetz(kz) = 0.0d0
        chgpionz(kz) = 0.0d0
        chgnionz(kz) = 0.0d0
        chgnetionz(kz) = 0.0d0
      ENDDO
        
      scion2(1) = 0.0
      scion2(2) = 0.0
      scion2(3) = 0.0
      scion2(4) = 0.0
      scion2(5) = 0.0
      scion2(6) = 0.0
      
      scx(:) = 0.0


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

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
        dv = dxx(ix)*dyy(jy)*dzz(kz)
        scion2(1) = scion2(1) + &
     &       Max(0.0,an(ix,jy,kz,lscpi))*dv ! /gt(ix,jy,kz,imapz)
        scion2(2) = scion2(2) + &
     &       Min(0.0,an(ix,jy,kz,lscpi))*dv ! /gt(ix,jy,kz,imapz)
        scion2(3) = scion2(3) + &
     &       Max(0.0,an(ix,jy,kz,lscni))*dv ! /gt(ix,jy,kz,imapz)
        scion2(4) = scion2(4) + &
     &       Min(0.0,an(ix,jy,kz,lscni))*dv ! /gt(ix,jy,kz,imapz)
        scion2(5) = scion2(5) + &
     &    (an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv ! /gt(ix,jy,kz,imapz)
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))

        scion2(6) = scion2(6) +  &
     &    (an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv ! /gt(ix,jy,kz,imapz)

       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz) &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  &
     &  + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz) &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni)) 
       ENDIF
      
        chgpionz(kz) = chgpionz(kz) + ec*(an(ix,jy,kz,lscpi))*dv
        chgnionz(kz) = chgnionz(kz) + ec*( - an(ix,jy,kz,lscni))*dv
        chgnetionz(kz) = chgnetionz(kz) +   &
     &         ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))*dv
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        chgposz(kz) = chgposz(kz) + t0(ix,jy,kz)*dv
        IF ( kz .le. kcldtop(ix,jy) ) chgposcld = chgposcld + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        IF ( kz .le. kcldtop(ix,jy) ) chgnegcld = chgnegcld + t0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        chgnegz(kz) = chgnegz(kz) + t0(ix,jy,kz)*dv
      END IF

     ! ion drift+corona tendency
      IF ( ichgtndion > 1 ) THEN
       IF ( ndebugd >= 1 ) write(0,*) 'write to ichgtndion = ',ichgtndion
        elec(ichgtndion)%flt3d(ix,jy,kz) = elec(ichgtndion)%flt3d(ix,jy,kz) + (t0(ix,jy,kz) - elec(iscnet)%flt3d(ix,jy,kz) )
      ENDIF

      elec(iscnet)%flt3d(ix,jy,kz) = t0(ix,jy,kz)
      
      IF ( ichgtndave > 1 ) THEN
       IF ( ndebugd >= 1 ) write(0,*) 'write to ichgtndave = ',ichgtndave
        elec(ichgtndave)%flt3d(ix,jy,kz) = elec(ichgtndave)%flt3d(ix,jy,kz) + t0(ix,jy,kz)*dt/thistory
      ENDIF
      
      end do
      end do
      end do
#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg
       mpitotindp(3)  = chgposcld
       mpitotindp(4)  = chgnegcld
       n = 4
       DO ia = 1,6
         n = n + 1
         mpitotindp(n) = scion2(ia)
       ENDDO
       
       IF ( microp(1:3) /= 'TAK' ) THEN
       DO ia = lscb,lsceq
         n = n + 1
         mpitotindp(n) = sc2(ia)
       ENDDO
       ENDIF
       

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
       chgposcld = mpitotoutdp(3)
       chgnegcld = mpitotoutdp(4)
       
       n = 4
       
       DO ia = 1,6
         n = n + 1
         scion2(ia) = mpitotoutdp(n)
       ENDDO

       IF ( microp(1:3) /= 'TAK' ) THEN
       DO ia = lscb,lsceq
         n = n + 1
         sc2(ia) = mpitotoutdp(n)
       ENDDO
       ENDIF
      ENDIF
#endif

      IF ( lstt .or. time .eq. ntstart ) THEN

#ifdef MPI

      allocate( mpitotinth((nz-1), 5))
      allocate(mpitotoutth((nz-1), 5))
      
       mpitotinth(1:nz-1, 1) =  chgposz(1:nz-1)
       mpitotinth(1:nz-1, 2) =  chgnegz(1:nz-1)
       mpitotinth(1:nz-1, 3) =  chgpionz(1:nz-1)
       mpitotinth(1:nz-1, 4) =  chgnionz(1:nz-1)
       mpitotinth(1:nz-1, 5) =  chgnetionz(1:nz-1)

      n = 5*(nz-1)

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
       chgposz(1:nz-1)  =  mpitotoutth(1:nz-1, 1) 
       chgnegz(1:nz-1)  =  mpitotoutth(1:nz-1, 2) 
       chgpionz(1:nz-1)  =  mpitotoutth(1:nz-1, 3) 
       chgnionz(1:nz-1)  =  mpitotoutth(1:nz-1, 4) 
       chgnetionz(1:nz-1)  =  mpitotoutth(1:nz-1, 5) 
      ENDIF
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )

#endif

      IF ( my_rank == 0 ) THEN
        write(3,'(a,i6)') 'Layer charge totals= ',time
        write(3,'(a)') 'kz,chgpos,chgneg,chgnet,posion,negion,netion'

#ifdef MPI
      kzb = ktile
      kze = 1
      if (kzend .eq. nzend) kzb = kzend-kzbeg

      do kz = kzb,kze,-1
#else
      DO kz = nz-1,1,-1
#endif
        chgpionz(kz) = chgpionz(kz) ! *dx*dy*dz/gt(1,1,kz,imapz)
        chgnionz(kz) = chgnionz(kz) ! *dx*dy*dz/gt(1,1,kz,imapz)
        chgnetionz(kz) = chgnetionz(kz) ! *dx*dy*dz/gt(1,1,kz,imapz)

        chgposz(kz) = chgposz(kz) ! *dx*dy*dz/gt(1,1,kz,imapz)
        chgnegz(kz) = chgnegz(kz) ! *dx*dy*dz/gt(1,1,kz,imapz)
          write(3,'(1x,i3,6(2x,1pe15.8))' )                       &
     &      kz,chgposz(kz),chgnegz(kz),chgposz(kz)+chgnegz(kz),   &
     &      chgpionz(kz),chgnionz(kz),chgnetionz(kz)
      ENDDO
      
      ENDIF ! my_rank
      ENDIF ! lstt .or. time .eq. ntstart

!C#ifdef MPI
!C! find global integrated rate max
!C       mpitotindp(1)  = chgpos
!C       mpitotindp(2)  = chgneg
!C
!C      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
!C
!C       
!C      IF ( my_rank == 0 ) THEN
!C       chgpos = mpitotoutdp(1)  
!C       chgneg = mpitotoutdp(2)  
!C      ENDIF
!C#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))') 'postlight pos/neg/net charge (C):',  &
     &      chgpos,chgneg,(chgpos+chgneg)
      write(iunit,'(a,3(2x,1pe15.8))') 'storm postlight pos/neg/net charge (C):', &
     &      chgposcld,chgnegcld,(chgposcld+chgnegcld)
      
      DO ia = 1,6
        write(iunit,'(a,i1,1(2x,1pe15.8))') 'scion2:',ia,scion2(ia)*ec
      ENDDO
     
!      IF ( lscb < 100 ) THEN
      IF ( microp(1:3) /= 'TAK' ) THEN
!      DO ia = lscb,lsceq
!        write(iunit,'(i3,1(2x,1pe15.8))') ia, sc2(ia)
        IF ( lscw > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'scwtot = ', sc2(lscw)
        IF ( lscr > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'scrtot = ', sc2(lscr)
        IF ( lsci > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'scitot = ', sc2(lsci)
        IF ( lscs > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'scstot = ', sc2(lscs)
        IF ( lsch > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'schtot = ', sc2(lsch)
        IF ( lschl > 1 ) write(iunit,'(a,1(2x,1pe15.8))') 'schltot = ', sc2(lschl)
!      ENDDO
      ENDIF
      ENDIF


!c
!c  END OF THE ELECTRICITY STUFF
!      IF ( my_rank == 0 ) THEN
!        write(iunit,*) 'NUMLGT',numlgt
!      ENDIF

      IF ( ndebugd >= 1 ) write(0,*) 'end of elec_drive4'


       RETURN
       END
