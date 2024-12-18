! CVD$F SKIP
! NOTE: fracn fix is enabled!!!!
!
! 11/2019: (erm) Swapped meanings of t2 and t4 arrays to match branched scheme and 
!                  to correct the "drho" centers, which were backwards (dipole was still correct)
!                Added some communication from ichoose rank to print out more info about
!                  the starting point and trajectory end points
!                Now start from the neighboring w-point with the largest |Ez|
!                  instead of at the scalar point, so that Emag should decrease in both 
!                  directions of the channel.
!
! 1.5.2004  treat ibal=2 same as ibal = 5 for leaving a fraction of charge on CG channels
!
! 05.13.2003  Added option ibrkd > 2 to use elgtthx as the max value for the 
!             breakeven field option.
!
! 10.27.2002 For ibal=5, allow some charge on downward part of CG (up to 25% of 
!            the charge on upward leaders)
!
! 8/25/2002  Fixed a minor bug.  For ibrkd=1 didn't set the values of zlev, so 
!            now it is fixed and the model levels print out correctly.
!
!
! 1/17/2002  Now have more variables that are from inlightdbm,
!             including 
!                   ieint (selects constant or fractional internal field)
!                   eint1 (input value for constant internal field)
!                   ibal  (selects channel neutrality balancing method)
!                   fefac (individual propagation settings for pos/neg ends)
!
! 1/4/2001  .1e.  Now use automatic arrays in SLPM... This removes the need to set
!           the slpm arrays manually.
!
!  11/26/2001 Merged the regular and stand-alone versions so that both are always
!             up to date.  Stand-alone behavior is controlled by 
!             variable 'isa' (>=1 for stand-alone, otherwise normal)
!
!  9/8/2000  sam.sll.dbm.omp.1a.f:
!          Stopped saving the charge from branches outside the cloud.  Now
!          changed to drop that charge and rescale the other end's charge
!          so that the IC remains neutral.
!
!          Also reversed sign of t4 array (use with sam.spark.2.3b.f or later)
!          so that t4 is positive where lightning charge (t3) is also positive.
!
! 6/27/99 new version to use stochastic lightning model for flashes
!         For now, compute IC flash discharge by using the smaller
!         magnitude of the two branching regions.
!         Not sure how to handle CG's yet.
!
!  CG's taken care of by slpm subroutine.
!
! 8/28/99 Add reporting of total positive and negative charge deposited
!         onto each hydrometeor category by a flash
!
! 11/06/99 pass in net charge in array t0, so don't need to recalculate it here.
!
!  START OF LIGHT
!
!
      subroutine mszlight    &
     &  (time_real,time,nstep,iliter,lgtstp,numlgt,loccur    &
     &  ,nor,nx,ny,nz    &
     &  ,na    &
     &  ,id1,jd1,kd1,istag,jstag,kstag,icldtop,zcldtop    &
     &  ,icgyn    &
     &  ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9    &
     &  ,tt0,elec,an,gx,gy,gz,gxt,gyt,gzt    &
     &  ,dxx,dyy,dzz    &
     &  ,ixst,jyst,kzst    &
     &  ,numic,numcgn,numcgp,outname,nkzmn,    &
     &   nnx,nny,nnz,    &
     &   nnxs,nnys,nnzs,    &
     &   iixp,jjyq,kkzr,iiex,jjey,kkez,    &
     &   nbw,nbe,nbs,nbn,nbz2,    &
     &   llwork, istretch, dzmax,    &
     &   nxslm,nyslm,nzslm,  & ! nb,nxb,nyb,nzbt,nzb,    &
     &   iunit,dn,pn, bcx, bcy, bcx1, bcy1, db, pb)
     
       USE GRID_MODULE
       USE INDEX_MODULE
       USE ELEC_MODULE, only: zgrnd,tgrnd,elgtfestopcg,ibal, cgfr, ibrkd,z1d2, &
                              elgtthx,elgtthn,elgtfdel,elgtdel,lightrad,ilight,  &
                              efracinitwire, overvolt, iseed, isa, elgt1,   &
                              igslg0, jgslg0, kgslg0, iGsdg, jGsdg, kGsdg, &
                              iGsdg0, jGsdg0, kGsdg0, dlz, eint, lightextendmsz
       USE COMMASMPI_MODULE
!
!
!      subroutine light
!     >  (nstep,ilight,iliter,nliter,lgtstp,numlgt,loccur,ilghtest
!     >  ,nor,nx,ny,nz,nx,ny,nz
!     >  ,ngt,na,id1,jd1,kd1,istag,jstag,kstag
!     >  ,iex,iey,iez,iemag,ipot,imapx,imapy,imapz,ig13,ig23,igsr
!     >  ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
!     >  ,elec,gt,ab,ac,ad,an)
!
!  include file containes indices in arrays a* and some of elec
!
!      include 'sam.jmshab.h'
!
#ifdef MPI
      use mpi
#endif
      implicit none

#ifdef MPI
!      INCLUDE "mpif.h"
       INTEGER :: my_mpi_status(MPI_Status_size)
#endif
      
      TYPE(VARIABLE)     :: elec(neelec)
      
      integer iunit
      
      real dzz(nz),dxx(nx),dyy(ny)         ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)

!      real gx(nx)
!      real gy(ny)
!      real gz(nz)

!      real    :: gxt(nx,4), gyt(ny,4), gzt(nz,4)
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)
      real    :: db(-nor+1:nz+nor)
      real    :: pb(-nor+1:nz+nor)
!      include 'sam.files.h'
      
      integer, parameter    :: ng1 = 1
      integer iixp,jjyq,kkzr,iiex,jjey,kkez
      integer*8 llwork
      integer nnx,nny,nnz,nnzt
      integer nnxs(2),nnys(2),nnzs(2)
      integer nnytmp
      integer nbw,nbe,nbs,nbn,nbz2
      integer nbntmp
      
      integer bcx, bcy, bcx1, bcy1  ! boundary condition flags
      
!      real,allocatable :: pot2(:,:,:) 
!      real pot2(nnx,nny,nnz)
!      real,allocatable :: rhs(:,:,:)
!      real,allocatable :: rho(:,:,:)
            
!      real work(llwork)

      real  dzmax
      logical istretch, ldumpsp
      integer nxslm,nyslm,nzslm
      integer nyslmtmp
      
      
      integer nb
!      parameter(nb=12)
      
      integer nxb,nyb,nzbt,nzb
!      parameter(nxb=nx+2*nb+1, nyb=ny+2*nb+1, nzbt=nz+nb)
!    nzb = nz+nb+1 
      
!      real f(nxb,nyb,nzbt)
!      real rho(nxb,nyb,nzbt)
!      real pot2(nxb,nyb,nzb) 
!      
!      real bdxs2(nyb,nzb),bdxf2(nyb,nzb)
!      real bdys2(nxb,nzb),bdyf2(nxb,nzb)

      character*100 line
      character*4 nmliterx
      character*6  rstime
      integer istat1
      integer icheck   !  flag to make slpm check that the dimensions are ok
      integer icldtop
      real    zcldtop
      integer na2
      parameter (na2=20)
      real chgrt(2,na2)
      real chgrtn,chgrtp
      real scal, scaln, scalp
      integer idrop
      real scdrop,scdropnet
      data scdropnet/0.0/
      double precision chg,sckeepp,sckeepn

      double precision chgneg,chgpos,chgtot,chgpnts ! total net negative/positive charge in storm (not including corona)
      
      integer iredo
      integer kzbrkd,iter,iread
      integer ix,jy,kz,kz1,loccur,i,j,k,n,m
      integer time,nstep,nor,nx,ny,nz,ngt,na,nba
      real    time_real
      integer id1,jd1,kd1,istag,jstag,kstag
      

!
!  vars for stand-alone lightning
!
      integer ifile
      integer ixst,jyst,kzst
      integer nkzmn(nz)
      character(len=80) fnelec,outname
      integer nlit,nliter2

      real frachl
      real frach
      real fracr
      real fracw
      real fraci
      real fracs
      real fractot


      integer pchnl,nchnl   !  number of segments in neg/pos leaders
      integer nttim,nenum,ninfo
      integer llx,lly,llz,lex,ley,lez
      integer lemag,lpot,lstop,lemag2,lcharge
      integer ndebug,itest
      integer lgtinit
      integer ieout,ndiscg
      integer ndisic,ndisca
      real pi,eperao,zsfc
      integer iliter,lgtstp,numlgt
      integer iopen,ifail
      data iopen /0/
      integer ilgt,ietotx,ietoty,ietotz
      integer :: nlgt1 = 0
      save nlgt1
      integer ixmx,jymx,kzmx,ixmn,jymn,kzmn
      integer iemaxx,iemaxy,iemaxz
      integer lgtgo
      integer islgt,jslgt,kslgt
      integer islgt_local,jslgt_local,kslgt_local
      real    t0mx,t0mn,scramt,cndamt,times
      integer in,it,ic,jc,kc,iep,je,ke
      integer :: icdom,jcdom,iepdom,jedom  ! domain-relative indices
      integer  idebug,idtsign, isign1,isign2
      real    eps,facx,facy,facz,facxe,facye,facze
      parameter (eps = 1.0e-7)
      real    telec,tcharge,dtcomp,tem
      integer i1,j1,k1,ichneu,idist
      real :: rdist
      real :: e1,e2
      integer ichpos,ichneg,sizetest
      real    scthend1,scthend2,scthendp,scthendn
      real    scthp,scthn,potthrp, potthrn
      integer isendn,jsendn,ksendn,isendp,jsendp,ksendp
      integer iic,ica,icg,icgp,icgn,i2,j2,k2,icgpflg,icgnflg
      integer iredis,iicp
      integer :: numic,numca=0,numcg=0
      integer numcgp,iicn,numcgn, icountp,icountn
      real    rndnum,cons1,cons2
      double precision scnetpi,scnetni,volzone
      integer irepeat,nscnetp,nscnetn,nscnet,nscnetnp
      double precision    voln,volp,scnetp,scnetn,scnet,scnetnp,scnetnn
      integer nscnetnn,ia
      
      character*40 flgtchn
      character*6 rstim,htail
      character*2 rnl
      character*11 yesno
      
!      real zgrnd   !  threshold height to declare channel end at ground
!      parameter (zgrnd = 1500.0)
      integer ichend(2)   ! tells whether channel end at
                         ! 1 = ground
                         ! 2 = cloud
                         ! 3 = air
!
      integer icgyn,icgpn
      integer ipend,inend
      real qci  ! interpolated value of cloud content at end of channel
      real chgmin
!      
!  vars for dipole calculation
!
      double precision qpx,qpy,qpz    ! dipole moment vector components
      double precision qpxn,qpyn,qpzn ! negative drho centroid
      double precision qpxp, qpyp, qpzp  ! positive drho centroid
      double precision qpxm, qpym, qpzm  ! dipole midpoint location
      double precision qpxcn,qpycn,qpzcn
      double precision qpxcp,qpycp,qpzcp  ! neg and pos centroid locations
      double precision qnetn,qnetp
      integer npntn,npntp
      double precision fracn,fracp,cgfrac
      double precision multp,multn
      double precision sctarg
      integer npsurfp,npsurfn  ! number points inside pot. surface
      integer nppntp,nppntn   ! number of potential flash points
      double precision :: total_nox
!      
!  vars for lightning channel 'cone'
!
!    rbase,theta define initial radius and angle of cone (read from inlight2)
!
!        channel
!           ^        edge of cone
!  \    |   |   |    /
!   \   |   |   |   /
!    \  |   |   |  /
!     \ |   |   |t/    (t=theta)
!      \|___|___|/
!           |<->|
!           rbase
      
      integer imx,jmy,kmz
      real rtest,chanlen,rchan,dl
      real,save :: tanth
      !  rtest is the distance from a gridpoint to the channel test point
      !  chanlen is the length along the channel of the test point
      !  rchan = rbase + chanlen*Tan(theta)
      !  tanth = Tan(theta)
      integer nchadd(2)  !  number of points added along channel
      integer isignt
      integer ichcnt
      real chcnt
      real delx,dely,delz  !  relative coords of channel point from corner
      
! variables from 'inlight2'
!      real    :: elgtthx  = 150.0e3          ! elgtthx (e-max before discharge)
!      real    :: elgtthn  = 115.0e3          ! elgtthn (end of lightning 'bolt' at e<elgtthn)
!      real    :: elgtdel  = 15.0e3           ! elgtdel (discharge possible [random] where emax-elgtthn<e) 
!      real    :: elgtfdel =  0.9             ! elgtfdel e magnitude fraction of delta for [f(z)] breakdown

      real    :: scth     = 0.11e-9          ! scth (space charge threshold for discharge region +/- scth)
      real    :: fscth    = 0.10             ! fscth (fraction of scth removed in iteration part of scheme)
      real    :: fprcnt   = 0.7             ! fprcnt (percent of charge removed in discharge volume)

! Alex test:
!      real    :: scth     = 0.01e-9          ! scth (space charge threshold for discharge region +/- scth)
!      real    :: fscth    = 0.90             ! fscth (fraction of scth removed in iteration part of scheme)
!      real    :: fprcnt   = 0.5             ! fprcnt (percent of charge removed in discharge volume)

      real    :: scthmn   = -1.0e-12         ! scthmn (minimum value [negative space charge] for work)
      real    :: scthmx   = 1.0e-12          ! scthmx (maximum value [positive space charge] for work)
      integer :: niter    = 2                ! niter (number of iterations on space charge threshold)
      integer :: icountth = 04               ! icountth (number of points required in discharge region)
      integer :: nretry   = 06               ! nretry (number of channel retries)
      real    :: rbase    = 250.0            ! rbase (meters)
      real    :: theta    = 22.0             ! theta (degrees)
!      real    :: elgtfestopcg = 0.15         ! fraction of Emag for stopping channel for CG
!
      real, parameter :: delq0 = 1.e-9       ! reference charge density change to scale NOx production
      
      integer :: idownward
      integer :: ifoundic
      
      save iread
      
       integer ninit  ! number of possible initiation points
      integer, save :: itermax ! maximum iterations for flash search
      integer :: itern, iterp ! which iteration of neg/pos channel volume to stop at
      real, allocatable, save :: chgavailiter(:),chgavailitersum(:) ! array to sum up available lightning charge in each wildfire iteration
      real, allocatable, save :: voliter(:),volitersum(:)
      real :: chanzn(nz),chanzp(nz), channtot,chanptot
      integer iretry   ! counter for nretry
       integer npick
       parameter(npick = 50000)
       integer ipick(npick,3)
       integer ichoose
       integer idum

      integer :: ichooserank = 0
      integer :: ichoosenumber
      integer :: ichoosenextrank

      integer ioffset, joffset, koffset
      integer ishift, jshift, kshift
      integer NLxchoose, NLychoose, NLzchoose
      
!      integer, parameter :: isa = 0
      
! variable of lnox calculation
      integer itt
      real amin,bmin,cmin
      real alnox,blnox  ! parameters from Wang et al. 1998 (JGR)
      real avogadro     ! Avogadro's number
      real air_molmass  ! molecular mass of air
      parameter(alnox=0.34e21,blnox=1.30e16,avogadro=6.0221415e23,air_molmass=28.96e-3)
      double precision :: tmp

!
!  nttim:  maximum number of intervals in a lightning channel
!  nenum:  number of ends to a channel from a common location
!          (a positive and a negative end)
!  ninfo:  number of parameters in the trje array 
!
      parameter  (nttim=4000,nenum=2,ninfo=12)
!
!  llx:    x location of trajectory of lightning channel in trje
!  lly:    y location of trajectory of lightning channel in trje
!  llz:    z location of trajectory of lightning channel in trje
!  lex:    x component of e-field along lightning channel in trje
!  ley:    y component of e-field along lightning channel in trje
!  lez:    z component of e-field along lightning channel in trje
!  lemag:  magnitude of e along lightning channel in traj
!  lpot:   potential along lightning channel in traj
!  lstop:  stop indictator along lightning channel in traj
!  lemag2: magnitude from lex,y,z of e along lightning channel in traj
!  lcharge:  charge along lightning channel in trje
!  llic,lljc,llkc : i,j,k of lower-left corner of x,y,z point
!  ndebug:  ndebug parameter (1=on)
!  itest:   itest    parameter (1=on)
!  ilgtinit: iligtinit is a parameter that must be = 1 for lightning
!            discharges to be calculated in the code
!  ibrkd:    ibrkd is they of break down; 1=based on emax; 2=based
!            on emax(z)
!  ieout:    ieout is interval to output lightning information
!
      parameter  (llx=1,lly=2,llz=3,lex=4,ley=5,lez=6)
      parameter  (lemag=7,lpot=8,lstop=9,lemag2=10)
      parameter  (lcharge=11)
      integer, parameter :: ltemg = 12
!      integer, parameter :: llic = 12, lljc = 13, llkc = 14
      parameter  (ndebug=1,itest=0)
      parameter  (lgtinit=1)
!      ,ibrkd=2)
      parameter  (ieout=1)
      integer :: find_indexlgt
!
!  ndis** = 1 (redistribute charge based on sfc area)
!  ndis** = 2 (redistribute charge based just removing a fraction)
!    cg = cloud-ground
!    ic = intracloud
!    ca = cloud-air
!
      parameter  (ndiscg=3)
      parameter  (ndisic=3)
      parameter  (ndisca=1)
!
      integer    nttimend(nenum),idxlgt(nenum,nttim,3),nttimendcg(nenum)
!
      parameter  (eperao=8.8592e-12)
      parameter  (zsfc=0.0)

      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
!
      real elgt
      double precision etot
      real dtelec(nenum),dt,dtp,dx,dy,dz
      real dxi,dyi,dzi
      real dxi2,dyi2,dzi2
!
      integer nemaxx(1000),nemaxy(1000),nemaxz(1000)
      integer kelec(1000)
      real emax, zlev(1000), ebrkd(1000), ebrkdp(1000)
      real emaxkz(1000),emaxlightkz(1000),exmaxkz(1000),eymaxkz(1000),ezmaxkz(1000)
      real emag, ez, zinit, pot
      save zlev,ebrkd,ebrkdp
      real,save :: estop(1000)
      real,save :: estopfrac
      real estopcg(1000)
      real trje(nenum,nttim,ninfo)

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

      real tt0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)

      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      save scdropnet
      
      integer nxyz,nwslv
      
      double precision dv
      real, parameter ::  rho00 = 1.225
      real :: dax, day, daz, phiw, phie, phis, phin, phid, phiu


      real,allocatable :: pot2(:,:,:) 
      real,allocatable :: rhs(:,:,:)
      real,allocatable :: rho(:,:,:)
      real,allocatable :: work(:)
      
      real,allocatable :: Qmask(:,:,:,:)
      
      real,    allocatable, save :: xy_init(:,:)
      integer, allocatable, save :: ij_reduce(:,:)
      real, allocatable, save :: flsh_map(:,:,:)
      integer, allocatable, save :: ijmask(:,:)
      integer, save :: nxm, nym
      integer              :: nr,nrflsh

      real, allocatable, save :: nox(:,:,:) ! (-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real, allocatable, save :: lnox1(:,:,:) ! (-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real noxtmp
      
      real  :: x,y,r
      
      real :: ran0mpi
      
      integer :: ie
      integer, save :: klgtmin = 3 ! lowest level to extend lightning: to avoid corona layers
      real, parameter :: zlgtmin = 600.


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

      integer, parameter :: ntot = 50
      integer :: ichoose_send_buffer(ntot)
      integer :: ichoose_recv_buffer(ntot)
      real    :: rchoose_send_buffer(ntot)
      integer :: mpitotinint(ntot), mpitotoutint(ntot)
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: nz1
      integer :: itype,imxmn
      real  mpitotin(ntot), mpitotout(ntot)
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      integer :: irank1,irank2
   
      integer :: listsize

      integer :: mxmn_size,max_rank,min_rank
      integer :: loc_send_size,loc_recv_size
      integer :: exyz_send_size,exyz_recv_size

      integer,allocatable,save :: mpiloc(:,:), mpiloc0(:,:,:)
      integer,allocatable,save :: loc_send_buffer(:,:),loc_recv_buffer(:,:)

      real,allocatable,save :: mpiexyz(:,:), mpiexyz0(:,:,:)
      real,allocatable,save :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      real,allocatable,save :: max_send_buffer(:,:),max_recv_buffer(:,:)
      real,allocatable,save :: min_send_buffer(:,:),min_recv_buffer(:,:)

      real,allocatable,save :: mpimxmn(:,:)
      integer :: init_send_size, init_recv_size, ninittmp, ninittot

      integer,allocatable,save :: init_send_buffer(:)
      integer,allocatable,save :: init_recv_buffer(:,:)
      integer,allocatable,save :: mpi_init0(:,:)
      
      
      integer, save :: boxmginit = 0

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      
      integer :: ntasks,mytask
      
! #####################################################################
!   Begin Executable Code
! #####################################################################


      ldumpsp = .false.
      ntasks = number_of_processes
      mytask = my_rank

      IF ( ndebug >= 2 ) THEN
       write(iunit,*) 'MSZ Start: ',time_real,time,nstep,iliter,lgtstp,numlgt,loccur    &
     &  ,nor,nx,ny,nz    &
     &  ,na    &
     &  ,id1,jd1,kd1,istag,jstag,kstag,icldtop,zcldtop    &
     &  ,icgyn    &
     &  ,dtp,dx,dy,dz
       ENDIF
!
!  read in variables
!
      if ( iread .eq. 0 ) then
!
      iread = 1
      
      icountth = Min(icountth, Max(2, Nint(icountth*1000./dx)))
!      
!      write(iunit,*) 'theta = ',theta 
      tanth = Tan(theta*4.*atan(1.)/180.0)

       itermax = Max( nxend, nyend, nz)
       itermax = Min( itermax-1, Int(100000./dx) ) ! limit lightning half-width to 100 km
       write(91,*) 'itermax = ',itermax
       IF ( .not. allocated( chgavailiter ) ) THEN
         allocate( chgavailiter(-itermax:itermax) )
         allocate( chgavailitersum(-itermax:itermax) )
         allocate( voliter(-itermax:itermax) )
         allocate( volitersum(-itermax:itermax) )
       ENDIF

!      write(6,*) 'mdel = ',mdel
!      IF ( mdel .lt. 10  ) THEN
!        write (iunit,*) 'Must have mdel at least = 10, resetting mdel=10'
!        mdel = 10
!      ENDIF
!      IF ( mdel .lt. 15 .and. dslight .lt. 251.00 ) THEN
!        write(iunit,*) 'Must have mdel at least =15 for dslight .le. 250'
!        mdel = 15
!      ENDIF
!
!      write(iunit,*) 'mdel = ',mdel
      
!      read(ifile,*) ibal
!      write(6,*) 'ibal = ',ibal
!      read(ifile,'(a)',end=101) line
!        read(line,*,iostat=istat1)  ibal, cgfr
!        IF ( istat1 .ne. 0 ) THEN
!          write(6,*) 'read error, iostat = ',istat1
          IF ( ibal .eq. 2 ) THEN
            IF(  cgfr .lt. 1.0 ) THEN
              write(iunit,*) 'setting cgfr ',1.0
              cgfr = 1.0
            ENDIF
          ELSE
          write(iunit,*) 'setting cgfr ',0.25
          cgfr = 0.25
          ENDIF
!        ENDIF
      write(iunit,*) 'ibal,cgfr = ',ibal,cgfr

!      read(ifile,'(a)',end=101) line
!        read(line,*,iostat=istat1)  cgthres
!        IF ( istat1 .ne. 0 ) THEN
!          write(iunit,*) 'read error, iostat = ',istat1
!          write(iunit,*) 'setting cgthres to default value of ',2001.
!          cgthres = 2001.
!        ENDIF
!       write(iunit,*) 'cgthres = ',cgthres
!      
!      GOTO 102
!      
! 101  CONTINUE
!      
!      write(iunit,*) 'No value given for cgthres. Setting to 2001. m'
!      cgthres = 2001.0
!      
! 102  CONTINUE
      
!
!      write(iunit,*) 'done reading inlight3'
!
!      IF ( isa .lt. 1 ) THEN
!        close(ifile)
!      ENDIF
!
!
!  compute ebrkz (breakdown field as a function of z)
!  function computed by Straka using Marshall and Rust plot
!  provided by MacGorman
!
      DO kz = 1,nz-1
       IF ( z1d2(kz,1) < zlgtmin ) THEN
         klgtmin = Max(kz,klgtmin)
       ENDIF
      ENDDO

       IF (ibrkd .ge. 2) THEN
       
      write(iunit,'(a)') 'Compute break-even field:'
      do kz = 1,nz
!      zlev(kz) = .5*dz + (kz-1)*dz
      zlev(kz) = z1d2(kz,1)
!      ebrkd(kz) = 197.81
!     >          - 1.9606e-2*zlev(kz)
!     >          + 5.7664e-7*zlev(kz)**2
      ebrkd(kz) = 1.208*167.0*Exp(-zlev(kz)/8.4e3)
      ebrkd(kz) = ebrkd(kz) * 1.e3
      kz1 = Min(kz,nz-1)
      ebrkdp(kz) = 284.e3*db(kz1)/rho00

      IF (ibrkd .eq. 2) THEN
        ebrkd(kz)= Min( ebrkd(kz), 125.0e3 )
      ELSEIF ( ibrkd .eq. 3 ) THEN
        ebrkd(kz)= Min( ebrkd(kz), elgtthx )
      ELSEIF ( ibrkd .eq. 4 ) THEN
        ebrkd(kz)= Min( ebrkdp(kz), 180.0e3 )
! Alex test:
!        ebrkd(kz)= Max( ebrkdp(kz), 65.0e3 )

      ELSEIF ( ibrkd .eq. 5 ) THEN
        ebrkd(kz)= Min( ebrkdp(kz), elgtthx )
      ENDIF
      
      IF ( ebrkd(kz) .ne. 0.0 ) THEN
        dv = ebrkdp(kz)/ebrkd(kz)
      ELSE
        dv = 0.0
      ENDIF
      
      write(iunit,'(i3,2x,f9.4,2x,1pe13.5,2x,1pe13.5,0pf9.3)')     &
     &          kz,1.0e-3*zlev(kz),ebrkd(kz), ebrkdp(kz), dv ! ebrkdp(kz)/ebrkd(kz)
      end do
!
       ELSE 

        write(iunit,*) 'Using uniform initiation threshold of ',elgtthx
        DO kz=1,nz
          ebrkd(kz) = elgtthx
         ! elgtfdel = (elgtthx - elgtdel)/elgtthx
          zlev(kz) = z1d2(kz,1)
        END DO
       END IF 
       
       IF (ibrkd .gt. 2) ibrkd = 2
       
       icheck = 1

      nxm = Int(abs(lightrad)/dx) + 1
      nym = nxm
      allocate( ijmask(-nxm:nxm,-nym:nym) )
      
      ijmask(:,:) = 0
      
       DO j = -nym,nym
        DO i = -nxm,nxm
          r = dx*Sqrt( Float(i**2 + j**2) )
          IF ( r <= Abs(lightrad) ) ijmask(i,j) = 1
        ENDDO
       ENDDO

       estopfrac = elgtthn/elgtthx
       IF (ibrkd .eq. 2 .or. ibrkd == 4 ) THEN
        DO kz=1,nz
!         estop(kz) = Min( (elgtthn)*ebrkd(kz)/elgtthx,elgtthn)
         estop(kz) =  estopfrac*ebrkd(kz)
         estopcg(kz) = elgtfestopcg*ebrkd(kz)
         IF ( my_rank == 0 ) write(0,*) 'kz,estop: ',kz,estop(kz),ebrkd(kz),elgtthn,elgtthx
        END DO
        
       ELSE 
        DO kz=1,nz
          estop(kz) = elgtthn
        END DO
       END IF 
      
      end if  ! (iread .eq. 0)

      IF ( lnox > 1 .and. .not. allocated( nox ) ) THEN
        allocate( nox(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) )
        allocate( lnox1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) )
      ENDIF

!
!
!
!  output from electrification model
!
!
!
!      htail = 'hsll'
!      if (iopen .ne. 1 .and. isa .lt. 1 ) then
!      iopen=1
!      open(24,file=hdrcfs(1:lhdrcfs)//filehead(1:lfilehead)//'.hlgttyp'
!     :   //htail,form='formatted',status='unknown')
!      end if
!
!  constants/flags
!
!  ifail:  flag to indicate to calling routine that 
!  there was a disasterous failure
!
      ifail = 0
!
!  ilgt:  flag to indicate lightning
!
      ilgt = 0
!
!  loccur:  another flag to indicate lightning occurs
!
      loccur = 0
!
!  irrational number
!
      pi = 4.*atan(1.)
! 
!  zero all temporary 3-d arrays but t0
!  how each is used is described later (I may equivalence these to
!  other names later). 
!
!      call zerond(nx,ny,nz,1,nor,t0) ! net charge -- is passed in now!
      t1(:,:,:) = 0.0
      t3(:,:,:) = 0.0
      t4(:,:,:) = 0.0
!      t5(:,:,:) = 0.0
      t6(:,:,:) = 0.0
!      t7(:,:,:) = 0.0
!      t8(:,:,:) = 0.0
!      t9(:,:,:) = 0.0


!
!  LIGHTNING PARAMETERIZATION 3: MacGorman, Straka, and Ziegler (2001)
!
!  Definitions of temporary arrays (I might equivalence these to more
!  meaningful names later.)
!  t0  RHOtot -- space charge density
!  t1  flag where brk down occurs
!  t2  +/- flag for + pot/rho or - pot/rho
!  t3  rho adj
!  t4  flag where lightning is
!  t5  pot
!  t6  total hydrometeor surface area
!  t7  flag where q>th, rho>th
!  t8  qc + qi
!  t9  condensate total
!
!  lightning based on intercloud, cloud-to-ground, and
!  cloud air discharges
!
      if ( ilight .eq. 3 ) then
!
      if ( iliter .eq. 1 ) then
      lgtstp = 0
      end if
!
!  Array holding flag (t7) showing co-position of cloud/condensate&charge
!  if scramt>1.e-14 (t0) and cndamt>1.e-7 (t9).  If both conditions are
!  met we consider ic discarge initiations.  If these are met and one 
!  end of the channel is within 500m of the ground then the discarge is a cg.
!  Otherwise, if one end of the channel is in the condensate/cloud and 
!  the other is not we have a cloud to air discharge.
!  Remember, total space charge is in t0
!
!      t0mx=0
!      t0mn=0
!      do kz = 1,nz-kd1
!      do jy = 1,ny-jd1
!      do ix = 1,nx-id1
!
!  be careful and zero some variables
!
!      t7(ix,jy,kz)  =  0.0
!      t8(ix,jy,kz)  =  0.0
!      t9(ix,jy,kz)  =  0.0
!
!  cloud amount total
!
!      DO ia=lqb,lqe
!      t9(ix,jy,kz)  =  t9(ix,jy,kz) 
!     >   + an(ix,jy,kz,ia) 
!      ENDDO

!      if ( t0(ix,jy,kz) .gt. t0mx ) then
!      t0mx = t0(ix,jy,kz) 
!      ixmx = ix
!      jymx = jy
!      kzmx = kz
!      end if
!      if ( t0(ix,jy,kz) .lt. t0mn ) then
!      t0mn = t0(ix,jy,kz) 
!      ixmn = ix
!      jymn = jy
!      kzmn = kz
!      end if
!
!  put potential into temporary variable t5
!
!      t5(ix,jy,kz) = elec(ipot)%flt3d(ix,jy,kz)
!
!  find cndamt (condensate amount) and put into t9
!
      
!      cndamt =
!     >  an(ix,jy,kz,lc) + an(ix,jy,kz,li) +
!     >  an(ix,jy,kz,lr) + an(ix,jy,kz,ls) +
!     >  an(ix,jy,kz,lf) + an(ix,jy,kz,lh) +
!     >  an(ix,jy,kz,lgl) + an(ix,jy,kz,lgm) +
!     >  an(ix,jy,kz,lgh) + an(ix,jy,kz,lir) +
!     >  an(ix,jy,kz,lip) + an(ix,jy,kz,lhl)
!      t9(ix,jy,kz)  =  cndamt
!
!  if abs(scramt).gt.1.e-14 and abs(cndamt).gt.1.e-7 then flag t7 as 1.01
!
!      t7(ix,jy,kz)  =  0.0
!      if ( abs( t0(ix,jy,kz) ).gt.1.e-14 .and.
!     >     abs(t9(ix,jy,kz)).gt.1.e-7 ) then
!      t7(ix,jy,kz)  = 1.01
!      end if
!      end do
!      end do
!      end do

#ifdef MPI
      nz1 = nz-kd1
      
      IF ( .not. allocated( mpimxmn ) ) THEN
      allocate(mpimxmn(nz1,4))
      allocate(max_send_buffer(2,nz1))
      allocate(max_recv_buffer(2,nz1))
      allocate(min_send_buffer(2,nz1))
      allocate(min_recv_buffer(2,nz1))

      allocate(mpiloc(nz1,2))
      allocate(mpiloc0(nz1,2,0:number_of_processes-1))
      allocate(loc_send_buffer(nz1,2))
      allocate(loc_recv_buffer(nz1*number_of_processes,2))

!      real,allocatable :: mpiexyz(:,:), mpiexyz0(:,:,:)
!      real,allocatable :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      allocate(mpiexyz(nz1,3))
      allocate(mpiexyz0(nz1,3,0:number_of_processes-1))
      allocate(exyz_send_buffer(nz1,3))
      allocate(exyz_recv_buffer(nz1*number_of_processes,3))

      allocate(init_send_buffer(2))
      allocate(init_recv_buffer(number_of_processes,2))

      allocate(mpi_init0(2,0:number_of_processes-1))
      ENDIF

#endif

     nr = NInt( lightrad/dx )
     nrflsh = Max(2, NInt(lightextendmsz/dx))
     IF ( .not. allocated( ij_reduce ) ) THEN
       allocate( ij_reduce(-nr+1:nx+nr,-nr+1:ny+nr))
       allocate( xy_init(-nr+1:nx+nr,-nr+1:ny+nr))
       allocate( flsh_map(-nrflsh+1:nx+nrflsh,-nrflsh+1:ny+nrflsh,2) )
     ENDIF


!
!  find largest emag, but first zero emag
!
      if ( ndebug .ge. 1 ) write(iunit,*) 'emax'
      emax = 0.0
      iemaxx = 1
      iemaxy = 1
      iemaxz = 1
!
!  find largest emag as a f(z) and locations, 
!  but first zero emag and locations
!
      ixb = 1
      ixe = itile
      
!      IF ( myproci == 1 .and. bcx1 /= 2 ) ixb = 2
      IF ( myproci == nproci .and. bcx1 /= 2 ) ixe = ixend - ixbeg
      
      jyb = 1
      jye = jtile
      
!      IF ( myprocj == 1 .and. bcy1 /= 2 ) jyb = 2
      IF ( myprocj == nprocj .and. bcy1 /= 2  ) jye = jyend - jybeg
      IF ( ny .le. 2 ) THEN
        jyb = 1
        jye = 1
      ENDIF

      do kz = 1,nz-kd1
!
      emaxkz(kz) = 0
      emaxlightkz(kz) = -1.
      nemaxx(kz) = 1
      nemaxy(kz) = 1
      nemaxz(kz) = kz
!
      do jy = jyb,jye ! Min(ny-1,2),Max(1,ny-2*jd1)
      do ix = ixb,ixe ! 2,nx-2*id1
!
      if ( elec(iemag)%flt3d(ix,jy,kz) .gt. emaxkz(kz) ) then
      nemaxx(kz) = ix
      nemaxy(kz) = jy
      nemaxz(kz) = kz
      emaxkz(kz) = elec(iemag)%flt3d(ix,jy,kz)
      exmaxkz(kz) = elec(iex)%flt3d(ix,jy,kz)
      eymaxkz(kz) = elec(iey)%flt3d(ix,jy,kz)
      ezmaxkz(kz) = elec(iez)%flt3d(ix,jy,kz)
      IF ( emaxkz(kz) >= emaxkz(kz) ) THEN
        emaxlightkz(kz) = emaxkz(kz)
      ENDIF
      end if
!
      if ( elec(iemag)%flt3d(ix,jy,kz) .gt. emax ) then
      iemaxx = ix
      iemaxy = jy
      iemaxz = kz
      emax = elec(iemag)%flt3d(ix,jy,kz)
      end if
!
      end do
      end do


#ifdef MPI
      mpimxmn(kz,1) = exmaxkz(kz) ! elec(iex)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,2) = eymaxkz(kz) ! elec(iey)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,3) = ezmaxkz(kz) ! elec(iez)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,4) = emaxkz(kz)
      mpiloc(kz,1) = nemaxx(kz)
      mpiloc(kz,2) = nemaxy(kz)
!      mpiloc(kz,3) = nemaxz(kz)

#endif
      end do
!

#ifdef MPI

      IF ( ntasks > 1 ) THEN
! This section only needed for terrain (i.e., WRF version)
!      DO kz = 1,nz1
!
!       exyz_send_buffer(kz,1) = emaxlightkz(kz)
!
!      ENDDO
!      mxmn_size = nz1
!
!! find out the max emag at each level, and which processor it is on
!      CALL MPI_Allreduce(exyz_send_buffer, exyz_recv_buffer, mxmn_size, MPI_REAL,  &
!     &                  MPI_MAX, my_comm, mpi_error_code)
!
!      IF (mpi_error_code /= MPI_SUCCESS) THEN
!
!        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Allreduce attempt"
!
!        CALL MPI_Abort(my_comm, mpi_error_code,mpi_abort_error_code)
!      ENDIF
!
!! put global maximum emax into mpimxmn
!      DO kz = 1,nz1
!
!       emaxlightkz(kz) = exyz_recv_buffer(kz,1)
!!       write(iunit,*) 'kz, emaxlightkz,ezmaxkz = ',kz,emaxlightkz(kz),ezmaxkz(kz)
!
!      ENDDO



      DO kz = 1,nz1

       max_send_buffer(1,kz) = mpimxmn(kz,4)
       max_send_buffer(2,kz) = real(my_rank)

      ENDDO
    
      mxmn_size = nz1

! find out the max emag at each level, and which processor it is on
      CALL MPI_Allreduce(max_send_buffer, max_recv_buffer, mxmn_size, MPI_2REAL,      &
     &                  MPI_MAXLOC, my_comm, mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Allreduce attempt"

        CALL MPI_Abort(my_comm, mpi_error_code,mpi_abort_error_code)
      ENDIF

! put global maximum emax into mpimxmn
      emax = 0.0
      DO kz = 1,nz1

       mpimxmn(kz,4) = max_recv_buffer(1,kz)
       emaxkz(kz) = max_recv_buffer(1,kz)

      ENDDO

! send e-field components of emag max
      DO kz = 1,nz1
       exyz_send_buffer(kz,1) = mpimxmn(kz,1)
       exyz_send_buffer(kz,2) = mpimxmn(kz,2)
       exyz_send_buffer(kz,3) = mpimxmn(kz,3)
      ENDDO

      exyz_send_size = nz1*3
      exyz_recv_size = nz1*3

      CALL MPI_Gather(exyz_send_buffer,exyz_send_size,MPI_REAL,    &
     &               exyz_recv_buffer,exyz_recv_size,MPI_REAL,     &
     &               0,my_comm,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiexyz0 = reshape(exyz_recv_buffer,(/nz1,3,number_of_processes/))

           DO kz = 1,nz1
       
            irank1 = int(max_recv_buffer(2,kz))

            exmaxkz(kz) = mpiexyz0(kz, 1, irank1)
            eymaxkz(kz) = mpiexyz0(kz, 2, irank1)
            ezmaxkz(kz) = mpiexyz0(kz, 3, irank1)

          ENDDO
        ENDIF


! send locations (x,y) of emag max
      DO kz = 1,nz1
       loc_send_buffer(kz,1) = ixbeg-1+mpiloc(kz,1)
       loc_send_buffer(kz,2) = jybeg-1+mpiloc(kz,2)
      ENDDO

      loc_send_size = nz1*2
      loc_recv_size = nz1*2

      CALL MPI_Gather(loc_send_buffer,loc_send_size,MPI_INTEGER,    &
     &               loc_recv_buffer,loc_recv_size,MPI_INTEGER,     &
     &               0,my_comm,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiloc0 = reshape(loc_recv_buffer,(/nz1,2,number_of_processes/))

!          DO itype = 1,2 ! type is really the ix or jy index
           DO kz = 1,nz1
       
            irank1 = int(max_recv_buffer(2,kz))

            mpiloc(kz,1) = mpiloc0(kz,1,irank1)
            mpiloc(kz,2) = mpiloc0(kz,2,irank1)
            nemaxx(kz) = mpiloc0(kz,1,irank1)
            nemaxy(kz) = mpiloc0(kz,2,irank1)

          ENDDO
!         ENDDO

        DO kz = nz1,1,-1
!         write(iunit,'(a,3(i4,1x),2(g15.3,1x))') 
!     :  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',mpimxmn(kz,4),
!     :  iemaxy, iemaxz, emax, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)

!         write(iunit,'(a,3(i4,1x),3(g15.5,1x))') 
!     :  'gIEMAXX, gIEMAXY, gIEMAXZ, gEMAX =: ',
!     :  mpiloc(kz,1) , mpiloc(kz,2) , kz, exmaxkz(kz),eymaxkz(kz),ezmaxkz(kz) !, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)
          
        ENDDO


      ENDIF !! my_rank==0
      
      ENDIF ! ntasks

#endif

      DO kz = 1,nz-1
      if ( emaxkz(kz) .gt. emax ) then
      iemaxx = nemaxx(kz)
      iemaxy = nemaxy(kz)
      iemaxz = kz
      emax = emaxkz(kz) 
      end if
      ENDDO

      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'Location and value of max(|E|)'
!      write(iunit,'(a,3(i4,1x),2(g15.3,1x))') 
!     :  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',iemaxx,
!     :  iemaxy, iemaxz, emax, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)

      write(iunit,'(a,3(i4,1x),2(g15.3,1x))')     &
     &  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',iemaxx,    &
     &  iemaxy, iemaxz, emax, emaxkz(iemaxz) ! elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)
!
      write(iunit,'(a)')     &
     & 'Z, NEMAXX, NEMAXY, NEMAXZ, EMAXKZ, EZ, EH, yesno =: '
     
       ENDIF
!      do kz = nz-kd1,1,-1
!      end do
!
!  check if lightning should occur (if emax .lt. elgtthx )
!  there are two methods....One based on a typical value
!  for the whole atmosphere....another based on a breakdown
!  efield that is a function of z
!
!  test if maximum breakdown approach is going to allow lightning
!  note that here lgtstp is the important parameter
!
      if ( ibrkd .eq. 0 ) then
      lgtstp = 0
      if ( emax .lt. elgtthx ) then
      IF ( my_rank == 0 ) write(iunit,*) 'NO LIGHTNING AT STEP',nstep
      yesno = 'NO'
      lgtstp = 1
      go to 2999
      else
      IF ( my_rank == 0 ) write(iunit,*) 'START LIGHTNING AT STEP',nstep
      lgtstp = 0
      end if
      end if
!
!  test if variable breakdown approach is going to allow lightning
!  note that here lgtgo is the important parameter
!
      if ( ibrkd .ge. 2 ) then
      
      lgtstp = 0
      lgtgo = 0
      do kz = nz-kd1,1,-1
      if ( emaxkz(kz) .lt. efracinitwire*overvolt*ebrkd(kz) ) then
!      write(iunit,'(a,i4,1x,f7.0)') 
!     :       'NO LIGHTNING AT STEP AND Z OF',nstep,zlev(kz)
      yesno = 'NO LIGHT'
      lgtstp = 1
      else
      kzbrkd = kz
!      write(iunit,'(a,i4,1x,f7.0,1x,i4)') 
!     >   'START LIGHTNING AT STEP AND Z OF', nstep,zlev(kz),kzbrkd
      yesno = 'START LIGHT'
      lgtgo = 1
      end if
      IF ( my_rank == 0 ) THEN
      write(iunit,'(3(i4,1x),1x,3(1x,1pe12.4),1x,a,1x,0pf7.0)')     &
     &   nemaxx(kz),nemaxy(kz), nemaxz(kz),    &
     &   emaxkz(kz), ezmaxkz(kz),    &
     &  Sqrt(exmaxkz(kz)**2 + eymaxkz(kz)**2 ),yesno,    &
     &   zlev(kz)
!      write(iunit,'(3i4,2x,3(1x,1pe12.4),1x,a,1x,0pf7.0)') 
!     >   nemaxx(kz),nemaxy(kz), nemaxz(kz),
!     >   emaxkz(kz), elec(iez)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz)),
!     >  Sqrt(elec(iex)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 +
!     >       elec(iey)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 ),yesno,
!     >   zlev(kz)
      
      ENDIF
      end do
      if ( lgtgo .eq. 0 ) go to 2999 
      IF (lgtgo .eq. 1) lgtstp = 0
      
      end if
      
!      write(iunit,*) 'Lightning not ready! Return!'
!      lgtstp = 1
!      GOTO 2999

!      IF ( isa .eq. 1 ) THEN
!        nliter2 = nliter
!      ELSE
        nliter2 = 1
!      ENDIF
      
      IF ( ny .le. 2 ) THEN
!        lgtstp = 1
!        RETURN
      ENDIF

#ifdef MPI
!      IF ( number_of_processes .ne. 1 ) THEN
!        write(0,*) 'lightning flag ON in MPI'
!        lgtstp = 1
!        RETURN
!      ENDIF

#endif
       t2(:,:,:) = 0.
       t4(:,:,:) = 0.
     
      DO nlit=1,nliter2

!#######################################################################
!  AT THIS POINT, IT HAS BEEN DECIDED THAT LIGHTNING HAPPENS . . .
!  NOW WE CAN CALL THE SLPM LIGHTNING MODULE . . . .     
!#######################################################################
!
      
! make list of possible initiation points
      ninit = 0
      DO kz = 2,nz-1
      DO jy = jyb,jye
      DO ix = ixb,ixe
       IF (  elec(iemag)%flt3d(ix,jy,kz)  >= ebrkd(kz)*elgtfdel ) THEN
        ninit = ninit + 1
        IF (ninit .gt. npick) THEN
          write(0,'(a)') 'Ee-gad! Over-ran the ipick array! ABORT!'
          call commasmpi_abort()
        END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
      ENDDO  ! k loop
      
    ! write(0,*) 'myrank, ninit = ',my_rank,ninit

#ifdef MPI
! send e-field components of emag max
      nz1 = 2
        init_send_buffer(1) = my_rank
        init_send_buffer(2) = ninit

      init_send_size = nz1
      init_recv_size = nz1

      CALL MPI_Gather(init_send_buffer,init_send_size,MPI_INTEGER,  &
     &               init_recv_buffer,init_recv_size,MPI_INTEGER,   &
     &               0,my_comm,mpi_error_code)

       IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "discharge_msz: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpi_init0 = reshape(init_recv_buffer,(/2,number_of_processes/))

           ninittot = 0
           DO kz = 0,number_of_processes-1
       
!        write(iunit,'(1x,i6,a,i6)')  mpi_init0(2, kz), ' possible starting points on rank =',kz
            ninittot = ninittot +  mpi_init0(2, kz)

          ENDDO
          
        write(iunit,'(1x,i6,a,i6)') ninittot, ' possible starting points' ! on rank =',my_rank
          

       ENDIF

       mpitotinint(1)  = ninit

      CALL MPI_AllReduce(mpitotinint, mpitotoutint, 1, MPI_INTEGER, MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       ninittot = mpitotoutint(1)  
          

#else

       ninittot = ninit

#endif

      IF ( my_rank == 0 ) write(iunit,'(1x,i6,a)') ninittot, ' possible starting points'
!
         IF ( ninittot == 0 ) THEN
           lgtstp = 1
           loccur = 0
           GOTO 2999
         ENDIF
         
         nretry = Max(1,ninittot/2)
!
! put continue here to make a new channel if previous one turned out 
! to be a dud.  Maybe implement a list of starting points....
!
      iretry = 0
      chgmin = (1.0-(fscth)*(niter-1))*scth

 1101 CONTINUE
       t2(:,:,:) = 0.
       t4(:,:,:) = 0.
       
       IF ( lnox > 1 ) THEN
       lnox1(:,:,:) = 0. 
       nox(:,:,:) = 0. 
       ENDIF
        
      iretry = iretry+1
      IF (iretry .gt. nretry) THEN
         write(iunit,*) 'Too many channel tries...abort. iretry,nretry = ',iretry,nretry
           lgtstp = 1
           loccur = 0
!         call commasmpi_abort()
         go to 2999
      END IF
            
      if ( ndebug .ge. 1 ) write(iunit,*)  'find lightning initiation pt. Try # ',iretry 
      IF (ninittot .eq. 0) THEN
        write(iunit, '(a)') ' No more starting points, try next timestep'
        lgtstp = 1
        loccur = 0
        GOTO 2999
      END IF

      ichoose = int(real(ninittot)*ran0mpi(iseed)) + 1
      IF ( my_rank == 0 ) THEN
      write (iunit,'(a,i6)') 'ichoose = ',ichoose
!      write(iunit,'(a,3(1x,i3))') 'Picked point ',islgt,jslgt,kslgt
      ENDIF

      ichooserank = 0
#ifdef MPI
! find which processor has the initiation point
      IF ( my_rank == 0 ) THEN
       ninittmp = 0
           ichooserank = -1
           DO m = 0,number_of_processes-1
       
!            write(iunit,*) 'm = ',m
            
            IF ( ichoose .le. ninittmp + mpi_init0(2, m) ) THEN
               write(iunit,*) 'ichooserank = ',mpi_init0(1,m)
               ichooserank = mpi_init0(1,m)
               ichoosenumber = ichoose - ninittmp
               EXIT
!            ELSE
!              write(iunit,*) 'init point not on processor ',mpi_init0(1,m)
            ENDIF
            
            ninittmp = ninittmp + mpi_init0(2, m)

          ENDDO
          
          IF ( ichooserank == -1 ) THEN
            write(0,*) 'Problem choosing init point! ichooserank = ',ichooserank
          ENDIF

          ichoose_send_buffer(1) = ichooserank
          ichoose_send_buffer(2) = ichoosenumber
          
      ELSE
          ichoose_send_buffer(1) = -1
          ichoose_send_buffer(2) = -1
      ENDIF ! (my_rank == 0 )

       
       CALL MPI_Bcast(ichoose_send_buffer, 2, MPI_INTEGER, 0, my_comm, mpi_error_code)
       ichooserank = ichoose_send_buffer(1)
       ichoose = ichoose_send_buffer(2)
       
       IF ( ndebug .ge. 1 )  write(iunit,*) 'my_rank: ichooserank,ichoose = ',my_rank,ichooserank,ichoose
#endif

       
      IF ( my_rank == ichooserank ) THEN
      
      i1 = ipick(ichoose,1) ! + Igslg - 1
      j1 = ipick(ichoose,2) ! + Jgslg - 1
      k1 = ipick(ichoose,3)
      ixst = i1
      jyst = j1
      kzst = k1
!      write(iunit,'(a,3(1x,i3))') 'Picked starting point ',i1,j1,k1
! update the points array (move up the points lower on the list)     
      DO idum=ichoose,ninit-1
        ipick(idum,1) = ipick(idum+1,1)
        ipick(idum,2) = ipick(idum+1,2)
        ipick(idum,3) = ipick(idum+1,3)
      END DO

          ichoose_send_buffer(1) = i1
          ichoose_send_buffer(2) = j1
          ichoose_send_buffer(3) = k1
          ichoose_send_buffer(4) = ixbeg-1 ! iGsdg - 1 !igslg0
          ichoose_send_buffer(5) = jybeg-1 ! jGsdg - 1 !jgslg0
          ichoose_send_buffer(6) = kzbeg-1 ! kGsdg - 1 !kgslg0
          ichoose_send_buffer(7) = Nx
          ichoose_send_buffer(8) = Ny
          ichoose_send_buffer(9) = Nz
          
          mpitotin(1) = elec(iemag)%flt3d(i1,j1,k1) ! emag
          mpitotin(2) = elec(iez)%flt3d(i1,j1,k1) ! ez
          mpitotin(3) = elec(ipot)%flt3d(i1,j1,k1) ! pot
          
          write(iunit,*) 'EMAG,k-1,k,k+1: ',elec(iemag)%flt3d(i1,j1,k1-1),elec(iemag)%flt3d(i1,j1,k1),elec(iemag)%flt3d(i1,j1,k1+1)
          write(iunit,*) 'EZ,k-1,k,k+1: ',elec(iez)%flt3d(i1,j1,k1-1),elec(iez)%flt3d(i1,j1,k1),elec(iez)%flt3d(i1,j1,k1+1)
          write(iunit,*) 'x,y,z = ',gxt(i1,1),gyt(j1,1),gzt(k1,1)
      
       ninit = ninit - 1
      ENDIF ! my_rank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
#endif
       i1 = ichoose_send_buffer(1)
       j1 = ichoose_send_buffer(2)
       k1 = ichoose_send_buffer(3)
       ixst = i1 + ichoose_send_buffer(4)
       jyst = j1 + ichoose_send_buffer(5)
       kzst = k1 + ichoose_send_buffer(6)
       ioffset = ichoose_send_buffer(4)
       joffset = ichoose_send_buffer(5)
       koffset = ichoose_send_buffer(6)
       NLxchoose = ichoose_send_buffer(7)
       NLychoose = ichoose_send_buffer(8)
       NLzchoose = ichoose_send_buffer(9)

#ifdef MPI
       CALL MPI_Bcast(mpitotin, 3, MPI_REAL, ichooserank, my_comm, mpi_error_code)
#endif
        emag = mpitotin(1)
        ez   = mpitotin(2)
        pot  = mpitotin(3)

      islgt = ixst
      jslgt = jyst
      kslgt = kzst
      
      ! reset ixst etc. to local index for the tile that has the init point and to zero for everyone else
      IF ( my_rank == ichooserank ) THEN
        ixst = i1
        jyst = j1
        kzst = k1
      ELSE
        ixst = 0 
        jyst = 0 
        kzst = 0
        kslgt = 1
      ENDIF

      islgt_local = i1
      jslgt_local = j1
      kslgt_local = k1

! one fewer init point now      
      ninittot = ninittot - 1
! find neighboring point with largest E and make it the first bond.

    IF ( my_rank == 0 ) THEN
      write(iunit,*) 'LIGHTNING INITIATION LOCATION AND |E|'
      write(iunit,*) 'NSTEP,ILITER =: ',nstep,iliter
      write(iunit,*) 'ISLGT, JSLGT, KSLGT = ',islgt_local+ioffset,jslgt_local+joffset, kslgt_local
      write(iunit, '(a,4(1x,1pe12.5))') 'emag,pot,ez,zc = ',emag,pot,ez,z1d2(kslgt_local,1)
    ENDIF
    
!    write(0,*) 'rank,ioffset,joffset = ',my_rank,ioffset,joffset,igsdg0,jgsdg0

!
! For now, the lightning channel will stay vertical to avoid making a streamline across an MPI patch boundary
!
!
!  draw a streamline along the efield lines until the efield threshold
!  is found and store the potential value for the + and - regions
!
      if ( ndebug .ge. 1 )write(iunit,*) 'generate +/- channels'
      dt = 6.0 ! dtp
      times = time_real ! dt*(nstep-1)
!
!  Initialize lstop--indicates trajectory end.
!
      trje(1,1,lstop) = 0.
      trje(2,1,lstop) = 0.
!
!  initiate position of trajectories
!
      trje(1,1,llx) = (islgt-1)*dx
      trje(2,1,llx) = (islgt-1)*dx
      trje(1,1,lly) = (jslgt-1)*dy
      trje(2,1,lly) = (jslgt-1)*dy
      trje(1,1,llz) =  z1d2(kslgt,1) ! (kslgt-1)*dz
      trje(2,1,llz) =  z1d2(kslgt,1) ! (kslgt-1)*dz
!
!  Set nttimend
!
      nttimend(1) = 0
      nttimend(2) = 0
      nttimendcg(1) = 0
      nttimendcg(2) = 0
!
! set end indicators to zero for safety      
      ichend(1) = 0
      ichend(2) = 0
!      
! zero out the number of points added in 'channel cone'
!
      nchadd(1) = 0
      nchadd(2) = 0

! chosen tile with starting point now draws a trajectory (vertical for now in MPI)
      IF ( my_rank == ichooserank ) THEN !{

       ! shift starting height to w/Ez point with largest efield:
       e1 = elec(iez)%flt3d(i1,j1,k1)   ! w-point below the scalar point
       e2 = elec(iez)%flt3d(i1,j1,k1+1)
       IF ( Abs(e1) > Abs(e2) ) THEN
         trje(1,1,llz) =  z1d2(kslgt,2) ! (kslgt-1)*dz
         trje(2,1,llz) =  z1d2(kslgt,2) ! (kslgt-1)*dz
       ELSE
         trje(1,1,llz) =  z1d2(kslgt+1,2) ! (kslgt-1)*dz
         trje(2,1,llz) =  z1d2(kslgt+1,2) ! (kslgt-1)*dz
       ENDIF
         
!
!
!  Trajectory number counter
!
      do in = 1,nenum
      
      ifoundic = 0
      idownward = 0
!
!  Iteration counter on drawing channel
!  (Number of iterations to find decenet channel and volume)
!
      do it = 1,nttim-1
!
!
!  actual interpolation (mlint1 does tri-linear interpolation)
!
!
!     subroutine mlint1
!    > (ndebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,
!    >  nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
!
!
!  find lower left corner of grid cell that trajectory
!  end point is currently in for interpolation of scalars
!
!      if ( ndebug .ge. 2 ) write(0,*) 'get ic, jc, kc, '

! subtract off igslg0/jgslg0 to translate to local coordinates

      icdom = ifix((trje(in,it,llx)+eps)/dx)  + 1
      ic = icdom - igslg0
      jcdom = ifix((trje(in,it,lly)+eps)/dy) + 1
      jc = jcdom  - jgslg0
      iepdom = ifix((trje(in,it,llx)+eps)/dx + 0.5) + 1  
      iep = iepdom - igslg0
      jedom = ifix((trje(in,it,lly)+eps)/dy + 0.5) + 1
      je = jedom - jgslg0
!      kc = ifix((trje(in,it,llz)+eps)/dz) + 1
      kc = find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)
      ke = find_indexlgt(trje(in,it,llz),z1d2(1,2),nz)
      facx = (trje(in,it,llx)-(icdom-1)*dx)/dx
      facy = (trje(in,it,lly)-(jcdom-1)*dy)/dy
!      facz = (trje(in,it,llz)-(kc-1)*dz)/dz
      facz = (trje(in,it,llz)- z1d2(kc,1))/dzz(kc) ! *z1d2(kc,3)
      facxe = (trje(in,it,llx)-(iepdom-1-0.5)*dx)/dx
      facye = (trje(in,it,lly)-(jedom-1-0.5)*dy)/dy
      facze = (trje(in,it,llz)- z1d2(ke,2))*gzt(ke,4) ! /(dlz*z1d2(ke,4)) ! dzz(kc) ! *z1d2(kc,3)
!
      idxlgt(in,it,1) = icdom
      idxlgt(in,it,2) = jcdom
      idxlgt(in,it,3) = kc
!     
!   print if in test mode.
!
      if ( ndebug .ge. 2 ) then
      write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
      write(iunit,*) 'FINDING SCALARS'
      write(iunit,*) 'IN=:',in,'IT=:',it,' FACX,Y,Z:',  facx, facy, facz
      write(iunit,*) ' FACXE,YE,ZE:',  facxe, facye, facze
      write(iunit,*) 'IN=:',in,' IC,JC,KC:',  IC,JC,KC
      write(iunit,*) 'iep,je,ke =',iep,je,ke
      write(iunit,*) trje(in,it,llz), z1d2(kc,1),z1d2(kc,3)
      write(iunit,*)  'ze(ke) = ',gzt(ke,2),z1d2(ke,2)
      write(iunit,*)  'dze = ',gzt(ke,4),z1d2(ke,4),dlz,z1d2(ke,4)*dlz
      write(iunit,*)  'dzz = ',dzz(ke),dzz(kc)
      end if
!
      if ( ic .ge. 1 .and. ic .le. ixe .and.    & 
     &     jc .ge. 1 .and. jc .le. jye .and.    & 
     &     kc .ge. 1 .and. kc .le. nz-2*kd1 ) then !{ {
!
!  Interpolate
!  efield stuff
!
      idebug = 0
      if ( ndebug .gt. 0 ) idebug = 1
      ie = iex
      call mlint2    & 
     & (idebug, elec(ie)%flt3d, -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, 1, 1,    & 
     &  nx, ny, nz, telec, facxe, facy, facz,iep,jc,kc, 1)
      trje(in,it,(lex-1+ie)) = telec              

      ie = iey
      call mlint2    & 
     & (idebug, elec(ie)%flt3d, -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, 1, 1,    & 
     &  nx, ny, nz, telec, facx, facye, facz,ic,je,kc, 1)
      trje(in,it,(lex-1+ie)) = telec              

      ie = iez
      call mlint2    & 
     & (idebug, elec(ie)%flt3d, -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, 1, 1,    & 
     &  nx, ny, nz, telec, facx, facy, facze,ic,jc,ke, 1)
      trje(in,it,(lex-1+ie)) = telec              

      do ie = iemag,lpot
      call mlint2    & 
     & (idebug, elec(ie)%flt3d, -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, 1, 1,    & 
     &  nx, ny, nz, telec, facx, facy, facz,ic,jc,kc, 1)
      trje(in,it,(lex-1+ie)) = telec              
      end do
!
!  Interpolate
!  total charge sctot(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
!
      call mlint2  &
     & (idebug, t0,  -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, &
     &  1, 1, nx, ny, nz, tcharge, facx, facy, facz,ic,jc,kc, 1)
      trje(in,it,lcharge) = tcharge              

      call mlint2  &
     & (idebug, tt0,  -nor+ng1, nx+nor, -nor+ng1, ny+nor, -nor+ng1, nz+nor, &
     &  1, 1, nx, ny, nz, tem, facx, facy, facz,ic,jc,kc, 1)
      trje(in,it,ltemg) = tem
!
!  Compute
!  lemag from lex,ley,lez
!
      trje(in,it,lemag2) =     & 
     &  ((trje(in,it,lex)**2)+(trje(in,it,ley)**2)    & 
     &  +(trje(in,it,lez)**2))**(0.5)
!
!  Set "pseudo" time step:  They are very small increments...With
!  a breakdown of 100 kV/m they are < 0.1 of a grid interval
!
      if ( in .eq. 1 ) idtsign = 1  ! toward neg. charge (pos. leader)
      if ( in .eq. 2 ) idtsign = -1 ! toward pos. charge (neg. leader)
!
      dtelec(in) = idtsign*dt*1.0e-7
!
      dtcomp = 0.0
      if ( it .gt. 1 ) then  
      dtcomp = (.05)*min(dx,dzz(kc))    & 
     &       / max(trje(in,it,lemag2),elgtthn)
      dtelec(in) = idtsign*dtcomp
      if ( ndebug .ge. 2 ) write(iunit,*) 'dtelec',dtelec(in),min(dx,dzz(kc)),    & 
     &  trje(in,it,lemag2),dtelec(in)*trje(in,it,lemag2)
      end if
!
!  Integrate forward in time (along streamline)
!
      IF ( .false. .and. ntasks == 1 ) THEN ! only do horizontal if single processor -- for now.
!      IF ( ntasks == 1 ) THEN ! only do horizontal if single processor -- for now.
      trje(in,it+1,llx) =     & 
     &  trje(in,it,llx) + dtelec(in) * trje(in,it,lex)
      trje(in,it+1,lly) =     & 
     &  trje(in,it,lly) + dtelec(in) * trje(in,it,ley)
      ELSE
      trje(in,it+1,llx) = trje(in,it,llx)
      trje(in,it+1,lly) = trje(in,it,lly)
      ENDIF
      trje(in,it+1,llz) =     & 
     &  trje(in,it,llz) + dtelec(in) * trje(in,it,lez)
     
     IF ( it == 1 .and. dtelec(in) * trje(in,it,lez) < 0.0 ) idownward = 1
!
!  Print debug stuff
!
      if ( ndebug .ge. 2 ) then
      write(iunit,*) 'TRAJ'
      write(iunit,*) 'in, it, times'
      write(iunit,*) 'trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)'
      write(iunit,*) 'trje(in,it+1,llx),trje(in,it+1,lly),trje(in,it+1,llz)'
      write(iunit,*) ' dtelec(in) * trje(in,it,lez), dtelec(in), trje(in,it,lez)'
      write(iunit,*) 'dt,dtp,idtsign,dtcomp'
      write(iunit,*) 'trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)'
      write(iunit,*) 'trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)'
      write(iunit,*) 'trje(in,it,lcharge),estop,kc'
      write(iunit,*) 'elec(ipot)%flt3d(ic+1,jc,kc,ipot),elec(ipot)%flt3d(ic+1,jc+1,kc,ipot)'
      write(iunit,*) 'elec()%flt3d(ic+1,jc,kc+1,ipot),elec()%flt3d(ic+1,jc+1,kc+1,ipot)'
      write(iunit,*) 'elec()%flt3d(ic,jc,kc,ipot),elec()%flt3d(ic,jc,kc,ipot)'
      write(iunit,*) 'Ex, Ey, Ez'
      write(iunit,*) in, it, times
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it+1,llx),trje(in,it+1,lly),trje(in,it+1,llz)
      write(iunit,*)  dtelec(in) * trje(in,it,lez), dtelec(in), trje(in,it,lez)
      write(iunit,*) dt,dtp,idtsign,dtcomp
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge),estop(kc),kc
      write(iunit,*) elec(ipot)%flt3d(ic+1,jc,kc),elec(ipot)%flt3d(ic+1,jc+1,kc)
      write(iunit,*) elec(ipot)%flt3d(ic+1,jc,kc+1),elec(ipot)%flt3d(ic+1,jc+1,kc+1)
      write(iunit,*) elec(ipot)%flt3d(ic,jc,kc),elec(ipot)%flt3d(ic,jc,kc)
      write(iunit,*) elec(iex)%flt3d(ic,jc,kc),elec(iey)%flt3d(ic,jc,kc),elec(iez)%flt3d(ic,jc,kc)
      write(iunit,*) elec(iex)%flt3d(iep,jc,kc),elec(iey)%flt3d(ic,je,kc),elec(iez)%flt3d(ic,jc,ke)
      write(iunit,*) elec(iex)%flt3d(iep+1,jc,kc),elec(iey)%flt3d(ic,je+1,kc),elec(iez)%flt3d(ic,jc,ke+1)
      write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
      end if

      nttimendcg(in) = it
!
! check for ground flash .... 
      if (  idownward == 1 .and. &
          (trje(in,it,llz) .lt. zgrnd .or. (zgrnd < 0. .and. trje(in,it,ltemg) > tgrnd)  ) &
         .and. &
       ( ( trje(in,1,lez) >  10.e3 .and. trje(in,1,lpot) < -50.e6 .and. trje(in,it,lcharge) > chgmin ) .or.   &
         ( trje(in,1,lez) < -10.e3 .and. trje(in,1,lpot) >  50.e6 .and. trje(in,it,lcharge) < -chgmin) ) ) then
      nttimend(in) = it
      nttimendcg(in) = it
      ichend(in) = 1 
      go to 1999
      end if


      if ( (trje(in,it,lemag2) .lt. Min(estop(kc),estopfrac*emag))  .and. ifoundic == 0 ) then
!      trje(in,it,lstop) = 10.
        IF ( ( in .eq. 1 .and. trje(in,it,lcharge) .lt. -chgmin ) .or.    & 
     &     ( in .eq. 2 .and. trje(in,it,lcharge) .gt.  chgmin ) .or.    & 
     &     trje(in,it,lemag2) .lt. 0.1*estop(kc) ) THEN
       IF ( ndebug .ge. 2 ) THEN
        write(iunit,*) 'END OF TRAJ for in = ',in
        write(iunit,*) 'kc,trje,estop = ',kc,trje(in,it,lemag2),estop(kc)
       ENDIF
        ifoundic = 1
        nttimend(in) = Max(1,it-1)
        ichend(in) = 2  ! set end type as 'cloud' for now -- check for air later
          IF ( idownward == 0 ) THEN
            go to 1999
          ENDIF
        END IF
      ELSEIF ( trje(in,it,lemag2) .lt. estopcg(kc) ) THEN
    ! get here if we hit the estopcg going downward, but didn't reach zgrnd
    ! leave nttimend and ichend at the values when estop was reached
          go to 1999
      end if
 
      
!
!  Else ic,jc, kc at boundary
!  set lstop > 1 and nttimend(in) = it-1
!
      else ! } {
!
!      trje(in,it,lstop) = 10.
      nttimend(in) = max(it-1,1)
      go to 1999
!
!  end if on computing trajectory
!
      end if ! } }
!
!  end do on it
!
      end do  ! it
!
1999  continue

      
      it = nttimend(in)

!
!  Bound computed values....set it+1 values = it (for plotting purposes)
!
      trje(in,it+1,llz) = max(0.0,trje(in,it,llz))
      trje(in,it+1,llz) = min( z1d2(nz-1,1),trje(in,it,llz))
      trje(in,it+1,lly) = max(0.0,trje(in,it,lly))
      trje(in,it+1,lly) = min((ny-2*jstag)*dy,trje(in,it,lly))
      trje(in,it+1,llx) = max(0.0,trje(in,it,llx))
      trje(in,it+1,llx) = min((nx-2*istag)*dx,trje(in,it,llx))
      trje(in,it+1,lex) = trje(in,it,lex)
      trje(in,it+1,ley) = trje(in,it,ley)
      trje(in,it+1,lez) = trje(in,it,lez)
      trje(in,it+1,lemag) = trje(in,it,lemag)
      trje(in,it+1,lemag2) = trje(in,it,lemag2)
      trje(in,it+1,lpot) = trje(in,it,lpot)
      trje(in,it+1,lcharge) = trje(in,it,lcharge)
      idxlgt(in,it+1,1) = idxlgt(in,it,1)
      idxlgt(in,it+1,2) = idxlgt(in,it,2)
      idxlgt(in,it+1,3) = idxlgt(in,it,3)
!     
      write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
      write(iunit,*) 'TRAJ END'
      write(iunit,*) in, it, times
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge)
      write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
!
      if ( ndebug.ge.1 ) write(iunit,*) 'nttimend',nttimend(1),nttimend(2)
!
!  dump lightning x, i coordinates
!
      if ( ilgt .gt. 0 ) then
!     write(rnl,911) int(iliter)
!911  format(i2.2)
!     write(rstime,912) int(times)
!912  format(i6.6)
!     flgtchn = 'hlgtchn'//rnl//rstime
!     open(unit=9,file=flgtchn,form='formatted')
!     do in = 1,nenum
!     do it = 1,nttimend(in)
!     write(9,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
!     end do
!     end do
!     close(9)
      end if
!
      end do  ! in = 1, nenum

      ENDIF ! } my_rank == ichooserank
!
#ifdef MPI
        CALL MPI_Bcast(nttimend, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
        CALL MPI_Bcast(nttimendcg, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
        CALL MPI_Bcast(ichend, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
     !   call mpi_bcast( nttimend

       ! send trje info to all processors
        n = nenum*nttim*ninfo
        CALL MPI_Bcast(trje, n, MPI_REAL, ichooserank, my_comm, mpi_error_code)
      
!      IF ( ichooserank /= 0 ) THEN
!       ! communicate trje from ichooserank to zero
!       IF ( my_rank == 0 .or. my_rank == ichooserank ) THEN
!        ! trje(nenum,nttim,ninfo)
!         n = nenum*nttim*ninfo
        
!         westward_tag = 10
!         IF ( my_rank == ichooserank ) THEN
!           Call MPI_Send(trje, n, MPI_REAL, 0, westward_tag, my_comm, mpi_error_code)
!         ENDIF

         IF ( my_rank == 0 ) THEN
!           Call MPI_Recv(trje, n, MPI_REAL, ichooserank, westward_tag, my_comm, my_mpi_status, mpi_error_code)
      DO in = 1,nenum
       it = nttimend(in)
   !   write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
      write(iunit,*) 'TRAJ END'
      write(iunit,*) in, it, times
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge),trje(in,it,ltemg)
       tmp = trje(in,it,llz)
       it = nttimendcg(in)
       IF ( trje(in,it,llz) < tmp ) THEN
      write(iunit,*) 'last downward point: ', it
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge),trje(in,it,ltemg)
      write(iunit,*) 'starting potential = ',trje(in,1,lpot)
      
      write(iunit,*) 'zgrnd,ez,pot = ',zgrnd,trje(in,1,lez),trje(in,1,lpot)

!      if ( trje(in,it,llz) .lt. zgrnd .and. idownward == 1 .and. &
!       ( ( trje(in,1,lez) >  10.e3 .and. trje(in,1,ipot) < -50.e6 ) .or.   &
!         ( trje(in,1,lez) < -10.e3 .and. trje(in,1,ipot) >  50.e6 ) ) ) then
        
        ENDIF

!          DO it = 2,nttimend(in)
!            write(iunit,*) 'it,dl = ',it,trje(in,it,llz)-trje(in,it-1,llz)
!          ENDDO
      write(iunit,*) '==================================================='
  !    write(iunit,*) '==================================================='
        
        ENDDO
      if ( ndebug.ge.1 ) write(iunit,*) 'nttimend',nttimend(1),nttimend(2)
        ENDIF
        
!       ENDIF
      
!      ENDIF
#endif
      

!  Determine which end of channel is at pos and neg charge.
!  Also, permit no air-to-air end point discharges or other un-physical
!  discharges.  Initiate isign1 & isign2 = 999
!
      isign1 = 999
      isign2 = 999
      
      chgmin = (1.0-(fscth)*(niter-1))*scth
!
!  Set flags for isign1:  if end of channel is in very small charge
!  then isign1 = 999
!
      IF (ichend(1) .gt. 1) THEN
!      if ( trje(1,nttimend(1),lcharge) .gt. (1.0e-14) ) then
      if ( trje(1,nttimend(1),lcharge) .gt. chgmin*0.01 ) then
      isign1 =  1
      end if
!      if ( trje(1,nttimend(1),lcharge) .lt. (-1.0e-14) ) then
      if ( trje(1,nttimend(1),lcharge) .lt. -chgmin*0.01 ) then
      isign1 = -1
      end if
      if ( abs(trje(1,nttimend(1),lcharge)) .lt. chgmin*0.01 ) then
      isign1 = 999
      end if
      
        IF ( trje(1,2,lpot) <  trje(1,1,lpot) .and. isign1 /= -1 ) THEN
          isign1 = -1
          trje(1,nttimend(1),lcharge) = -chgmin*0.1
        ELSEIF ( trje(1,2,lpot) >  trje(1,1,lpot) .and. isign1 /= 1 ) THEN
          isign1 = 1
          trje(1,nttimend(1),lcharge) = chgmin*0.1
        ENDIF
      
      END IF  !  (ichend(1) .ne. 1) 
!
!  Set flags for isign2:  if end of channel is in very small charge
!  then isign2 = 999
!
      IF (ichend(2) .gt. 1) THEN
      if ( trje(2,nttimend(2),lcharge) .gt. chgmin*0.01 ) then
      isign2 =  1
      end if
      if ( trje(2,nttimend(2),lcharge) .lt. -chgmin*0.01 ) then
      isign2 = -1
      end if
      if ( abs(trje(2,nttimend(2),lcharge)) .lt. chgmin*0.01 ) then
      isign2 = 999
      end if

        IF ( trje(2,2,lpot) >  trje(2,1,lpot)  .and. isign2 /= 1  ) THEN
          isign2 = 1
          trje(2,nttimend(2),lcharge) = chgmin*0.1
        ELSEIF ( trje(2,2,lpot) <   trje(2,1,lpot)  .and. isign2 /= -1  ) THEN
          isign2 = -1
          trje(2,nttimend(2),lcharge) = -chgmin*0.1
        ENDIF

      END IF  !  (ichend(2) .ne. 1) 
!

      IF ( my_rank == ichooserank ) THEN !{

!  Now calculate lnox using by chanlen and dv
        
     IF ( lnox > 1 ) THEN
     DO in = 1,2 
      chanlen = 0.0
      itt=2
        DO it=itt,nttimend(in)
            ic = ifix((trje(in,it,llx)+eps)/dx) + 1 - igslg0
            jc = ifix((trje(in,it,lly)+eps)/dy) + 1 - jgslg0  ! set lower-left corner
            kc = idxlgt(in,it,3) ! find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)
         dl = Sqrt(    &
     &      (trje(in,it,llx) - trje(in,it-1,llx))**2 +    &
     &      (trje(in,it,lly) - trje(in,it-1,lly))**2 +    &
     &      (trje(in,it,llz) - trje(in,it-1,llz))**2 )
         chanlen = chanlen + dl

!     IF (chanlen .ge. dzz(ic,kc,jc)) THEN
     IF (chanlen .ge. dl .and. dl > 0.0) THEN
!        amin=4000.
!          DO ix=ic-1,ic+1
!            IF (Abs(trje(in,it,llx)-(ix-1)*dx).lt.amin) THEN
!              amin=Abs(trje(in,it,llx)-(ix-1)*dx)
!                IF (amin.lt.1.) THEN
!                  ic=ix
!                END IF
!             END IF
!           END DO
           
        cmin=500.
!         DO kz=kc-1,kc+1
!            IF (Abs(trje(in,it,llz)-(kz-1)*dzz(kz)).lt.cmin) THEN
!             cmin=Abs(trje(in,it,llz)-(kz-1)*dzz(kz))
!!                IF (cmin.lt.1.) THEN
!                   kc=kz
!                END IF
!             END IF
!           END DO
           
!        bmin=4000.
!          DO jy=jc-1,jc+1
!            IF (Abs(trje(in,it,lly)-(jy-1)*dy).lt.bmin) THEN
!              bmin=Abs(trje(in,it,lly)-(jy-1)*dy)
!                IF (bmin.lt.1.) THEN
!                   jc=jy
!                END IF
!             END IF
!           END DO
      ! Here, chanlen is the channel length, pn = pressure, dn=air density, dv=dx*dy*dzz (cell volume)
              dv=dx*dy*dzz(kc)
              noxtmp = dl*(alnox+blnox*(pn(ic,jc,kc)+pb(kc)))/avogadro
              nox(ic,jc,kc) = nox(ic,jc,kc) + noxtmp
              lnox1(ic,jc,kc)=lnox1(ic,jc,kc)+air_molmass*noxtmp/(dn(ic,jc,kc)*dv)    
      
      END IF
          
        END DO
        
          itt = it  
          chanlen = 0.0
        
      END DO
      
      ENDIF ! lnox > 1
            
! Now look for points within a widening cone along the channel in 
!  each direction that ends above ground
!
      DO in = 1,2
       IF (in .eq. 1) isignt = isign1 ! sign of ambient charge at end of channel
       IF (in .eq. 2) isignt = isign2
        nchadd(in) = 0
      IF (ichend(in) .gt. 1 .and. isignt .ne. 999) THEN
        chanlen = 0.0
        chcnt = 1.0
        DO it=2,nttimend(in)
         chanlen = chanlen + Sqrt(    &
     &      (trje(in,it,llx) - trje(in,it-1,llx))**2 +    &
     &      (trje(in,it,lly) - trje(in,it-1,lly))**2 +    &
     &      (trje(in,it,llz) - trje(in,it-1,llz) )**2 )
          IF (chanlen .gt. dzz(kslgt)*chcnt) THEN
            chcnt = chcnt + 1.0  ! increase to next check length
            rchan = (rbase + chanlen*tanth)  ! radius of sphere to look into         
            ic = ifix((trje(in,it,llx)+eps)/dx) + 1 - igslg0
            jc = ifix((trje(in,it,lly)+eps)/dy) + 1 - jgslg0  ! set lower-left corner
!            kc = ifix((trje(in,it,llz)+eps)* z1d2(kz,3)) + 1
            kc = idxlgt(in,it,3) ! find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)
            imx = Int(rchan/dx) 
            jmy = Int(rchan/dy)   ! how many grid lengths can radius reach
            kmz = Int(rchan/(dlz*z1d2(kc,3)))
            delx = trje(in,it,llx) - Real(ic - 1)*dx ! relative coords of
            dely = trje(in,it,lly) - Real(jc - 1)*dy ! channel point from
            delz = trje(in,it,llz) - z1d2(kc,1) ! Real(kc - 1)*dz ! lowerleft corner
            rchan = (rchan)**2  ! use r^2 to check and save doing square roots
            
            DO kz = Max(1,kc-kmz),Min(nz,kc+kmz+1)
            DO jy = Max(jyb,jc-jmy),Min(jye,jc+jmy+1)
            DO ix = Max(ixb,ic-imx),Min(ixe,ic+imx+1)
               rtest = (dx*Real(ix-ic) - delx)**2 +    &
     &                 (dy*Real(jy-jc) - dely)**2 +    &
     &                 (dzz(kc) - delz)**2
               IF (rtest .le. rchan .and.     &
     &             Abs(t2(ix,jy,kz)) .lt. 0.5)  THEN
               IF ( Int(Sign(1.1,t0(ix,jy,kz))) .eq. isignt .and.    &
     &               Abs(t0(ix,jy,kz)) .gt. scth ) THEN
                    t2(ix,jy,kz) = -2.1*Real(isignt) 
                    t4(ix,jy,kz) = -1.01*Real(isignt)
                    nchadd(in) = nchadd(in) + 1
                   IF ( lnox > 1 ) THEN
                    dv=dx*dy*dzz(kz)
                    noxtmp = Sqrt(rtest)*(alnox+blnox*pn(ix,jy,kz))/avogadro   !here rtest is chanlen
                    nox(ix,jy,kz) = nox(ix,jy,kz) + noxtmp
                    lnox1(ix,jy,kz) = lnox1(ix,jy,kz)+air_molmass*noxtmp/(dn(ix,jy,kz)*dv)    
                   ENDIF
               END IF
               END IF        
            END DO
            END DO
            END DO
          END IF ! (chanlen .gt. 500.0*chcnt)
        END DO  ! it
        
        write(iunit,'(2(a,i3))') 'added ',nchadd(in),    &
     &       ' points along channel end ',in
        
      END IF    ! (ichend(in) .gt. 1 .and. isignt .ne. 999) 
      END DO    ! in
      


      IF (ichend(1) .gt. 1) THEN
! decide between cloud/precip vs. CG
      in = 1
!      it = nttimend(in)
!      i1 = idxlgt(in,it,1)
!      j1 = idxlgt(in,it,2)
!      k1 = idxlgt(in,it,3)
!      facx = (trje(in,it,llx)-(i1-1)*dx)/dx
!      facy = (trje(in,it,lly)-(j1-1)*dy)/dy
!      facz = (trje(in,it,llz)-z1d2(k1,1))/dzz(i1,k1,j1)

! leaving out the check for cloudy/precip vs. clear air for now
!       IF ( k1 > 5 ) THEN 
!        ichend(in) = 2 ! IC
!       ENDIF
        
      ENDIF

! leaving out the check for cloudy/precip vs. clear air for now for end 2


      ENDIF ! } my_rank == ichooserank


!
!
!  so, now have isign = 999 for a flash at ground or in clear air
!  otherwise, isign shows sign of charge in the cloud/precip at the
!  ends of the channel.  Also have isign = 999 if space charge is
!  below the absolute minimum we can go to.
!
!  Channel will fail for ground-ground or air-air or ground-air
!
      ifail = 0
      
   !  write(iunit,'(a,2(i3))') 'ichend(1),ichend(2)',ichend(1),ichend(2)
      IF ( ichend(1) .eq. 0 .or. ichend(2) .eq. 0 ) THEN
        IF ( my_rank == ichooserank .or. my_rank == 0 ) THEN
        write(iunit,'(a)') 'FAILED LIGHTNING: one end not determined'
        ENDIF
        ifail = 1 !GOTO 1101
      END IF 
     
     IF ( my_rank == ichooserank .or. my_rank == 0 ) THEN
      if ( ndebug .ge. 1 ) write(iunit,*) 'isign1,isign2',isign1,isign2
      ENDIF
!
      if ( isign1 .eq. isign2 .and. ifail == 0 ) then
        IF ( my_rank == ichooserank .or. my_rank == 0 ) THEN
      write(iunit,*) 'Same sign of charge for both ends of lightning'
      write(iunit,*) 'FAILED LIGHTNING END SIGN PROBLEM-try again'
      write(iunit,*) 'isign1,trje(1,nttimend(1),lcharge)'
      write(iunit,*) isign1,trje(1,nttimend(1),lcharge)
      write(iunit,*) 'isign2,trje(2,nttimend(2),lcharge)'
      write(iunit,*) isign2,trje(2,nttimend(2),lcharge)
      ENDIF
!      go to 2999
      lgtstp = 1
      ifail = 1 ! GOTO 1101
      end if
!      
! check for ends in cloud but not enough charge
!  
! Change this to check instead for added points along channel....
! if at least two points were added, then allow flash.
      IF (ichend(1).eq.2 .and. isign1.eq.999 .and.  (nchadd(1) .lt. icountth .and. ifail == 0)  ) THEN
       write(iunit,'(a)') 'FAILED LIGHTNING at end1: in cloud but < chgmin'
      ifail = 1 ! GOTO 1101
      end if
      IF (ichend(2).eq.2 .and. isign2.eq.999 .and. (nchadd(2) .lt. icountth .and. ifail == 0) ) THEN
       write(iunit,'(a)') 'FAILED LIGHTNING at end2: in cloud but < chgmin'
      ifail = 1 ! GOTO 1101
      end if


      
!      IF ( my_rank /= ichooserank ) ifail = 0
      
!#ifdef MPI
!       ichoose_send_buffer(1) = ifail
!       ichoose_send_buffer(2) = lgtstp
!       CALL MPI_Bcast(ichoose_send_buffer, 2, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
!       ifail = ichoose_send_buffer(1)
!       lgtstp = ichoose_send_buffer(2)
!#endif
       
       IF ( ifail > 0 ) GO TO 1101

!
!  Set indices to be used in rest of code to indicate which channel 
!  ends in positive charge and which in negative charge.
!
      ichneu = 0
      ichpos = 0
      ichneg = 0
!      if ( idxlgt(1,nttimend(1),3) .gt. 1 ) then
      IF ( ichend(1) .eq. 2 ) THEN
      if ( trje(1,nttimend(1),lcharge) .eq. 0.0 ) then
      ichneu=1
      end if
      if ( trje(1,nttimend(1),lcharge) .gt. 0.0 ) then
      ichpos=1
      ichneg=2
      end if
      if ( trje(1,nttimend(1),lcharge) .lt. 0.0 ) then
      ichneg=1
      ichpos=2
      end if
      end if
!
!      if ( idxlgt(2,nttimend(2),3) .gt. 1 ) then
      IF ( ichend(2) .eq. 2 ) THEN
      if ( trje(2,nttimend(2),lcharge) .eq. 0.0 ) then
      ichneu=2
      end if
      if ( trje(2,nttimend(2),lcharge) .gt. 0.0 ) then
      ichpos=2
      ichneg=1
      end if
      if ( trje(2,nttimend(2),lcharge) .lt. 0.0 ) then
      ichpos=1
      ichneg=2
      end if
      end if
!
      if ( ndebug .ge. 1 ) then
      write(iunit,*) 'ichneu,ichpos,ichneg',ichneu,ichpos,ichneg
      end if
!
!  Setting end points of lighting and other assorted goodies...
!
!      iretry = 0
!      iterst = 1
!
! 88887 continue

!
!  Iteration to ensure that space charge and potential values at end of
!  each channel is adequate and within tolerance chosen by user.
!  After niter iterations an attempt to discharge is aborted and is 
!  retried later.
!
!      do 89999 iter = iterst,niter
!
!  scthend# are temporary values for the charge at the end of 
!  channel 1 and 2.
!
      scthend1 =  1.0e10
      scthend2 = -1.0e10
      IF (ichend(ichpos).eq.2) scthend1 = trje(1,nttimend(1),lcharge)
      IF (ichend(ichneg).eq.2) scthend2 = trje(2,nttimend(2),lcharge)
!
!  scthendp and scthendn are the values of charge at the ends of 
!  the channels. 
!  
      scthendp =  1.0e10
      scthendn = -1.0e10
      IF (ichend(ichpos).eq.2) scthendp=trje(ichpos,nttimend(ichpos),lcharge)
      IF (ichend(ichneg).eq.2) scthendn=trje(ichneg,nttimend(ichneg),lcharge)
!
!  scthp and scthn are working values (minimums) of the space charge 
!  at the end of the channel required to outline the area of the 
!  discharge.
!
!      IF (ichend(ichpos).eq.2) scthp = scth ! Min(scth,scthendp,Abs(scthendn))
!      IF (ichend(ichneg).eq.2) scthn = -scth !-Min(scth,scthendp,Abs(scthendn))
      IF (ichend(ichpos).eq.2) scthp = Min(scth,scthendp,Abs(scthendn))
      IF (ichend(ichneg).eq.2) scthn = -Min(scth,scthendp,Abs(scthendn))
!
!  potthrp and potthrn are the values of potential at the ends of 
!  the channels.
!
      potthrp = 1.0e20
      potthrn = -1.0e20
      IF (ichend(ichpos).eq.2) potthrp=max(trje(1,nttimend(1),lpot),trje(2,nttimend(2),lpot))
      IF (ichend(ichneg).eq.2) potthrn=min(trje(1,nttimend(1),lpot),trje(2,nttimend(2),lpot))
!
!  isendn, jsendn, ksendn are the indices of the lower left hand corner
!  of the grid volume where the end of the negative channel is found.
!
      isendn = idxlgt(ichneg,nttimend(ichneg),1)
      jsendn = idxlgt(ichneg,nttimend(ichneg),2)
      ksendn = idxlgt(ichneg,nttimend(ichneg),3)
!
!  isendp, jsendp, ksendp are the indices of the lower left hand corner
!  of the grid volume where the end of the positive channel is found.
!
      isendp = idxlgt(ichpos,nttimend(ichpos),1)
      jsendp = idxlgt(ichpos,nttimend(ichpos),2)
      ksendp = idxlgt(ichpos,nttimend(ichpos),3)
!
      if ( ndebug .ge. 1 ) then
      write(iunit,*) 'scthend1',trje(1,nttimend(1),lcharge)
      write(iunit,*) 'scthend2',trje(2,nttimend(2),lcharge)
      write(iunit,*) 'sendn',isendn,jsendn,ksendn
      write(iunit,*) 'sendp',isendp,jsendp,ksendp
      end if
      if ( ndebug.ge.1 ) then
      write(iunit,*) 'scthp',scthp
      write(iunit,*) 'scthn',scthn
      write(iunit,*) 'scthendp',scthendp
      write(iunit,*) 'scthendn',scthendn
      write(iunit,*) 'potthrp',potthrp
      write(iunit,*) 'potthrn',potthrn
      write(iunit,*) 'ichend p/n = ',ichend(ichpos),ichend(ichneg)
      end if

      IF ( my_rank == ichooserank ) THEN !{
      ENDIF ! } my_rank == ichooserank
      

#ifdef MPI
    IF ( ntasks > 1 ) THEN
      IF ( my_rank == ichooserank ) THEN

          rchoose_send_buffer(1) = isendn
          rchoose_send_buffer(2) = jsendn
          rchoose_send_buffer(3) = ksendn
          rchoose_send_buffer(4) = isendp
          rchoose_send_buffer(5) = jsendp
          rchoose_send_buffer(6) = ksendp
          rchoose_send_buffer(7) = scthp
          rchoose_send_buffer(8) = scthn
          rchoose_send_buffer(9) = potthrp
          rchoose_send_buffer(10) = potthrn
          rchoose_send_buffer(11) = ichpos
          rchoose_send_buffer(12) = ichneg
          rchoose_send_buffer(13) = ichend(1)
          rchoose_send_buffer(14) = ichend(2)

      ENDIF ! } my_rank == ichooserank


       k = 14
       
       CALL MPI_Bcast(rchoose_send_buffer, k, MPI_REAL, ichooserank, my_comm, mpi_error_code)

       IF ( my_rank /= ichooserank ) THEN
       isendn = Nint( rchoose_send_buffer(1) )
       jsendn = Nint( rchoose_send_buffer(2) )
       ksendn = Nint( rchoose_send_buffer(3) )
       isendp = Nint( rchoose_send_buffer(4) )
       jsendp = Nint( rchoose_send_buffer(5) )
       ksendp = Nint( rchoose_send_buffer(6) )
       scthp   = rchoose_send_buffer(7)
       scthn   = rchoose_send_buffer(8)
       potthrp = rchoose_send_buffer(9)
       potthrn = rchoose_send_buffer(10)
       ichpos    = Nint( rchoose_send_buffer(11) )
       ichneg    = Nint( rchoose_send_buffer(12) )
       ichend(1) = Nint( rchoose_send_buffer(13) )
       ichend(2) = Nint( rchoose_send_buffer(14) )
       
       ENDIF

    ENDIF
#endif
       
!       IF ( lgtstp .gt. 0 ) GOTO 2999

!  Count number of ICs and CGs and CGPs and CGNs and CAs
!
      ilgt = 0
      iic  = 0
      ica  = 0
      icg  = 0
      icgp = 0
      icgn = 0
      icgpn = 0
      icgyn = 0

!      i1=idxlgt(1,nttimend(1),1)
!      j1=idxlgt(1,nttimend(1),2)
!      k1=idxlgt(1,nttimend(1),3)
!      i2=idxlgt(2,nttimend(2),1)
!      j2=idxlgt(2,nttimend(2),2)
!      k2=idxlgt(2,nttimend(2),3)
      
      IF ( ichend(1) == 1 .or. ichend(2) == 1 ) THEN
        icgyn = 1
        IF ( ichend(ichpos) .eq. 2 ) icgpn = 1
        IF ( ichend(ichneg) .eq. 2 ) icgpn = -1
      ENDIF

!
! Is discharge an ic
!
      if ( icgyn .eq. 0 ) then

      iic = 1
      ilgt = 1
      iredis = ndisic
      numic = numic + 1
      numlgt = numlgt + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      IC DISCHARGE'
      write(iunit,*) '         TOTAL NUMBER OF ICs',numic
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!      write(24,'(i5,a,f13.3,3x,i4)') nstep,' IC ', dtp*(nstep - 1),
!     : iliter
      ENDIF

      end if
!

      if ( icgyn .eq. 1 ) then
!
      icg = 1
      ilgt = 1
      iredis = ndiscg
      numcg = numcg + 1
      numlgt = numlgt + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      CG DISCHARGE'
      write(iunit,*) '         TOTAL NUMBER OF CGs',numcg
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
!      icgyn = 1
!
!
      if ( icgpn .eq. 1 ) then
      iicp = 1
      ilgt = 1
      loccur = 3
      numcgp = numcgp + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      DISCHARGE IS POSITIVE'
      write(iunit,*) '         TOTAL NUMBER OF CGPs',numcgp
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!        write(24,'(i5,a,f13.3,3x,i4)') nstep,' CG POS ', 
!     :  dtp*(nstep - 1), iliter
      ENDIF
!      end if
!
      elseif ( icgpn .eq. -1 ) then
      iicn = 1
      ilgt = 1
      loccur = 2
      numcgn = numcgn + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      DISCHARGE IS NEGATIVE'
      write(iunit,*) '         TOTAL NUMBER OF CGNs',numcgn
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!       write(24,'(i5,a,f13.3,3x,i4)') nstep,' CG NEG ', dtp*(nstep - 1),
!     : iliter
      ENDIF
      end if
!
      end if  ! ( icgyn.eq.1 ) 


!
!
!  flag the regions using t2 to state where there is potential for
!  lightning discharge:  sc>scth and pot>potthr 
!  also find scdel
!
!  Set t2 arrays...the first guess array...
!
!    number of points in each potential surface
      npsurfp = 0
      npsurfn = 0
      nppntp = 0
      nppntn = 0
      
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe

! calculate electrical energy in t5 array (units of MJ)
                dv = dxx(ix)*dyy(jy)*dzz(kz)
                dax = dyy(jy)*dzz(kz)
                day = dxx(ix)*dzz(kz)
                daz = dxx(ix)*dyy(jy)
! Have tried interpolating phi to the faces, but it doesn't seem to work as well as
! using the scalar point phi for all faces.
                phiw =  elec(ipot)%flt3d(ix,jy,kz)
                phie =  elec(ipot)%flt3d(ix,jy,kz)
                phis =  elec(ipot)%flt3d(ix,jy,kz)
                phin =  elec(ipot)%flt3d(ix,jy,kz)
                phid =  elec(ipot)%flt3d(ix,jy,kz)
                phiu =  elec(ipot)%flt3d(ix,jy,kz)

               t5(ix,jy,kz) = 1.e-6*0.5*eperao*(elec(iemag)%flt3d(ix,jy,kz)**2*dv + &
                  ( -elec(iex)%flt3d(ix  ,jy  ,kz  )*phiw*dax    &
                    +elec(iex)%flt3d(ix+1,jy  ,kz  )*phie*dax    &
                    -elec(iey)%flt3d(ix  ,jy  ,kz  )*phis*day    &
                    +elec(iey)%flt3d(ix  ,jy+1,kz  )*phin*day    &
                    -elec(iez)%flt3d(ix  ,jy  ,kz  )*phid*dax    &
                    +elec(iez)%flt3d(ix  ,jy  ,kz+1)*phiu*daz    &
                  ) )

! already zeroed and may have added points already...
!       t2(ix,jy,kz) = 0.0
      ! use rectilinear distance for potential change along channel
!      idist = Abs((ixbeg + ix) - islgt) + Abs( (jybeg + jy ) - jslgt)
      rdist = Sqrt(float( ((ixbeg + ix) - islgt)**2 + ( (jybeg + jy ) - jslgt)**2 ))
      if ( ichend(ichpos) .eq. 2 ) then
!        if ( elec(ipot)%flt3d(ix,jy,kz) .gt. potthrp + eint*(dx*idist + Abs(zlev(kslgt) - zlev(kz))) ) then
        if ( elec(ipot)%flt3d(ix,jy,kz) .gt. potthrp + eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt) - zlev(kz))**2) ) then
!         if ( elec(ipot)%flt3d(ix,jy,kz) .gt. potthrp ) then ! + eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt) - zlev(kz))**2) ) then
        npsurfp = npsurfp + 1
          if ( t0(ix,jy,kz) .gt. scthp ) then
!            t2(ix,jy,kz) = Max( 1.01, t2(ix,jy,kz) )
            t2(ix,jy,kz) = Min( -1.01, t2(ix,jy,kz) )
            nppntp = nppntp + 1
          end if
        end if
      end if
      if ( ichend(ichneg) .eq. 2) then
!        if ( elec(ipot)%flt3d(ix,jy,kz) .lt. potthrn - eint*(dx*idist + Abs(zlev(kslgt) - zlev(kz)))  ) then
        if ( elec(ipot)%flt3d(ix,jy,kz) .lt. potthrn - eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt) - zlev(kz))**2)  ) then
!        if ( elec(ipot)%flt3d(ix,jy,kz) .lt. potthrn ) then ! - eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt) - zlev(kz))**2)  ) then
         npsurfn = npsurfn + 1
          if ( t0(ix,jy,kz) .lt. scthn ) then
!            t2(ix,jy,kz) = Min( -1.01, t2(ix,jy,kz) )
            t2(ix,jy,kz) = Max( 1.01, t2(ix,jy,kz) )
            nppntn = nppntn + 1
          end if
        end if
      end if

      end do
      end do
      end do
      
#ifdef MPI
! find global integrated rate max
         mpitotinint(1)  = npsurfp
         mpitotinint(2)  = npsurfn
         mpitotinint(3)  = nppntp
         mpitotinint(4)  = nppntn


        CALL MPI_AllReduce(mpitotinint, mpitotoutint, 4, MPI_INTEGER, MPI_SUM, my_comm, mpi_error_code)

        npsurfp = mpitotoutint(1)
        npsurfn = mpitotoutint(2)
        nppntp = mpitotoutint(3)
        nppntn = mpitotoutint(4)
       
#endif

      write(iunit,*) 'potthrp,potthrn = ',potthrp,potthrn
      write(iunit,*) 'pnts in pos/neg pot. surfs ',npsurfp,npsurfn
      write(iunit,*) 'pot pnts in pos/neg regions',nppntp,nppntn
!
!  search for positive continuous discharge region
!
      itype = -1
      cons1 = 0.0
      cons2 = 0.0
      icountp = 0
      icgnflg = 0
      if ( ichend(ichpos) .eq. 2 ) then
!      subroutine searchl    &
!     &  (itype,icount,is,js,ks,id1,jd1,kd1    &
!     &  ,nx,ny,nz,nxl,nyl,nzl,nor    &
!     &  ,cons1,cons2,tthr,tflg)

      IF ( mytask == ichooserank ) THEN
! set starting point
        t4(isendp-igslg0,jsendp-jgslg0,ksendp) = itype*itermax
      ENDIF
      
      call searchl   &
     &  (itype,itermax,icountp,isendp,jsendp,ksendp,id1,jd1,kd1   &
     &  ,nx,ny,nz,nx,ny,nz,nor   &
     &  ,cons1,cons2,t2,t4,icgnflg,klgtmin)
      end if
!
!  search for negative continuous discharge region = positive leader breakdown region
!
      itype = 1
      cons1 = 0
      cons2 = 0
      icountn = 0
      icgpflg = 0
      if ( ichend(ichneg) .eq. 2 ) then
      IF ( mytask == ichooserank ) THEN 
! set starting point
        t4(isendn-igslg0,jsendn-jgslg0,ksendn) = itype*itermax
      ENDIF
      call searchl   &
     &  (itype,itermax,icountn,isendn,jsendn,ksendn,id1,jd1,kd1   &
     &  ,nx,ny,nz,nx,ny,nz,nor   &
     &  ,cons1,cons2,t2,t4,icgpflg,klgtmin)
      end if

      IF ( lnox > 1 ) THEN
! calculate lnox in positive and negative continuous discharge regions
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      IF (abs(t4(ix,jy,kz)).lt.abs(itype)*itermax.and.abs(t4(ix,jy,kz)).gt.1.01) THEN
       dv=dx*dy*dzz(kz)
!       nox(ix,jy,kz)=min(dx,dzz(kz))*(alnox+blnox*pn(ix,jy,kz))/(dn(ix,jy,kz)*avogadro)
!        nox(ix,jy,kz)=dzz(kz)*(alnox+blnox*pn(ix,jy,kz))/(db(ix,jy,kz)*avogadro)
         noxtmp = dx*(alnox+blnox*pn(ix,jy,kz))/avogadro
         nox(ix,jy,kz) = nox(ix,jy,kz) + noxtmp
       lnox1(ix,jy,kz) = lnox1(ix,jy,kz)+air_molmass*noxtmp/(dn(ix,jy,kz)*dv)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      IF ( icgpflg /= 0 .or. icgnflg /= 0 ) THEN
!        write(iunit,*) 'CG? icgpflg,icgnflg = ',icgpflg,icgnflg
        IF ( icgpflg == 0 .or. icgnflg == 0 ) THEN
           icgpflg = Min( 1, icgpflg )
           icgnflg = Min( 1, icgnflg )
        ELSE ! both ends at ground, so treat as IC
         icgpflg = 0
         icgnflg = 0
        ENDIF
      ENDIF


!
      if ( ndebug .ge. 1 ) then
      write(iunit,*) 'ICNP: icountn,icountp=',icountn,icountp
      end if
!
      sizetest = 1

      if ( icountn .lt. icountth .and. ichend(ichneg) .eq. 2 )  then
      write(iunit,*) 'icountn = ',icountn
      write(iunit,*) 'RETRY:IRETRY:ICOUNTN:STEP ',iretry,icountn,nstep
!      iterst = iterst+1
!      iretry = iretry+1
!      if ( iterst .gt. niter ) then
      write(iunit,*) 'Cant find enough points..try a new starting point'
      
       sizetest = 0
      
!      go to 2999
!      end if
!      go to 1101
      end if
!
      if ( icountp .lt. icountth .and. ichend(ichpos) .eq. 2 ) then
      write(iunit,*) 'icountp = ',icountp
!      iterst = iterst+1
      write(iunit,*) 'RETRY:IRETRY:ICOUNTP:STEP ',iretry,icountp,nstep
!      iretry = iretry+1
!      if ( iterst .gt. niter ) then
      write(iunit,*) 'Cant find enough points..try a new starting point'
      
       sizetest = 0
      
!      go to 2999
!      end if
!      go to 1101
      end if

      IF ( sizetest == 0 ) THEN
        numlgt = numlgt - 1
        IF ( icgyn .eq. 0 ) THEN
          numic = numic - 1
        ELSE
          numcg = numcg - 1
          IF ( icgpn .eq. 1 ) THEN
            numcgp = numcgp - 1
          ELSEIF ( icgpn .eq. -1 ) THEN
            numcgn = numcgn - 1
          ENDIF
        ENDIF
       
        go to 1101
        
      ENDIF


!
!  define t3 using t4 to define regions where actual discharge 
!  will occur...
!
      irepeat = 0
59999 continue
      irepeat = irepeat +1
      scnetpi = 0.0
      scnetni = 0.0
!
      xy_init(:,:) = 0
      flsh_map(:,:,:) = 0.0
      chgavailiter(:) = 0.0
     
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      
      t3(ix,jy,kz) = 0.0
       dv=dx*dy*dzz(kz)

       IF ( t4(ix,jy,kz) /= 0.0 ) THEN
         xy_init(ix,jy) = 1. ! using xy_init to designate the lightning footprint
         IF ( t4(ix,jy,kz) > 0.5 ) THEN
           flsh_map(ix,jy,1) = 1
         ELSE
           flsh_map(ix,jy,2) = 1
         ENDIF
       ENDIF

      if ( ichend(ichneg) .eq. 2) then
       if ( t4(ix,jy,kz) .gt. 0.5 ) then
 !       if ( t0(ix,jy,kz) .gt. scthp ) then
        if ( t0(ix,jy,kz) .lt. scthn ) then
!         t3(ix,jy,kz) = -fprcnt*Min(0.0, (t0(ix,jy,kz)+scth) )
         t3(ix,jy,kz) = -fprcnt*Min(0.0, (t0(ix,jy,kz)) )
          n = Nint( t4(ix,jy,kz) )
         chgavailiter(n) = chgavailiter( n ) + t3(ix,jy,kz)*dv
!         write(0,*) 'myrank, t4,t3 = ',my_rank,t4(ix,jy,kz),t3(ix,jy,kz)
!         write(0,*) 'myrank, n,chgavil= ',my_rank,n,chgavailiter( n )
!         IF ( t2(ix,jy,kz) .gt. 1.5 )  t3(ix,jy,kz) = -fprcnt*t0(ix,jy,kz)
         scnetpi = scnetpi + t0(ix,jy,kz)
        end if
       end if
      end if
      
      if (ichend(ichpos) .eq. 2 ) then
       if ( t4(ix,jy,kz) .lt. -0.5 ) then
!        if ( t0(ix,jy,kz) .lt. scthn ) then
        if ( t0(ix,jy,kz) .gt. scthp ) then
!         t3(ix,jy,kz) = -fprcnt*Max(0.0, (t0(ix,jy,kz)-scth) )
         t3(ix,jy,kz) = -fprcnt*Max(0.0, (t0(ix,jy,kz)) )
          n = Nint( t4(ix,jy,kz) )
    !     write(0,*) 'myrank, t4,t3 = ',my_rank,t4(ix,jy,kz),t3(ix,jy,kz)
         chgavailiter(n) = chgavailiter(n) + t3(ix,jy,kz)*dv
!         IF ( t2(ix,jy,kz) .lt. -1.5 ) t3(ix,jy,kz) = -fprcnt*t0(ix,jy,kz)
         scnetni = scnetni + t0(ix,jy,kz)
        end if
       end if
      end if
      
      end do
      end do
      end do

!
#ifdef MPI
       mpitotindp(1) = scnetpi
       mpitotindp(2) = scnetni

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetpi = mpitotoutdp(1)
       scnetni = mpitotoutdp(2)
       
       n = 2*itermax+1
       CALL MPI_AllReduce(chgavailiter, chgavailitersum, n, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
       
       chgavailiter(-itermax:itermax) = chgavailitersum(-itermax:itermax)
       chgavailitersum(:) = 0.0
       
#endif
        
        IF ( my_rank == 0 ) THEN
          write(iunit,*) 'scnetpi,scnetni = ',scnetpi,scnetni
        ENDIF

        chgavailitersum(itermax) = chgavailiter(itermax)
        chgavailitersum(-itermax) = chgavailiter(-itermax)
!          IF ( my_rank == 0 ) THEN
!            write(iunit,*) 'n, chgavailitersum(n), chgavailitersum(-n)'
!          ENDIF
        DO n = itermax-1,1,-1
          chgavailitersum( n) = chgavailitersum( n+1) + chgavailiter( n)
          chgavailitersum(-n) = chgavailitersum(-n-1) + chgavailiter(-n)
!          IF ( my_rank == 0 ) THEN
!            write(iunit,*) n, chgavailitersum(n), chgavailitersum(-n),chgavailiter(n)
!          ENDIF
        ENDDO


#ifdef MPI


      IF ( number_of_processes .gt. 1 .and. nr > 0 .and. ninittot > 0 ) THEN
        
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,1,nr,nr,0,nr,1,                &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,xy_init)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,xy_init)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,xy_init)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,xy_init)


        westward_tag = 201
        CALL sendrecv_westward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,flsh_map)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,flsh_map)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,flsh_map)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,flsh_map)
    
      ENDIF
#endif


      ij_reduce(:,:) = 0
      
      IF ( nr > 0 ) THEN ! create circular region around each flash point 
       DO jy = -nr+1,ny+nr
        DO ix = -nr+1,nx+nr
          IF ( xy_init(ix,jy) > 0.5 ) THEN
            
            DO j = Max(-nr+1,jy-nr),Min(ny+nr, jy+nr)
             DO i = Max(-nr+1,ix-nr),Min(nx+nr, ix+nr)
               ij_reduce(i,j) = Max( ij_reduce(i,j) , ijmask( i-ix, j-jy ) )
             ENDDO
            ENDDO
            
          ENDIF
        ENDDO
       ENDDO
      ELSE
        ij_reduce(:,:) = 1
      ENDIF

       ! for IC flash, compare footprints of positive and negative discharge regions
       ! to limit extent of one beyond the other
       IF ( ichend(ichneg) .eq. 2 .and. ichend(ichpos) .eq. 2 ) THEN
       DO jy = 1,ny
        DO ix = 1,nx
          IF ( flsh_map(ix,jy,1) > 0.5 .or. flsh_map(ix,jy,2) > 0.5 ) THEN
            
            ! search local box for points in footprints of pos and neg drho
            i1 = 0
            i2 = 0
            DO j = Max(-nrflsh+1,jy-nrflsh),Min(ny+nrflsh, jy+nrflsh)
             DO i = Max(-nrflsh+1,ix-nrflsh),Min(nx+nrflsh, ix+nrflsh)
               i1 = i1 + Nint(flsh_map(i,j,1))
               i2 = i2 + Nint(flsh_map(i,j,2))
             ENDDO
            ENDDO
            
            IF ( flsh_map(ix,jy,1) > 0.5 .and. i2 == 0 ) THEN
              flsh_map(ix,jy,1) = 0.0
              t4(ix,jy,1:nz-1) = 0.0
              t3(ix,jy,1:nz-1) = 0.0
            ENDIF
            IF ( flsh_map(ix,jy,2) > 0.5 .and. i1 == 0 ) THEN
              flsh_map(ix,jy,2) = 0.0
              t4(ix,jy,1:nz-1) = 0.0
              t3(ix,jy,1:nz-1) = 0.0
            ENDIF
            
          ENDIF
        ENDDO
       ENDDO
       ENDIF


!
!      if ( ndebug .ge. 1 ) then
!      write(iunit,*) 'scnetni,scnetpi',scnetni,scnetpi 
!      end if
!
!  accumululate scnet and nscnet
!
      if ( ndebug .ge. 1 ) write(iunit,*) 'set scnet,nscnet'
!      volzone = dx*dy*dz
      voln = 0.0 ! icountn*volzone
      volp = 0.0 ! icountp*volzone
      scnetp = 0.0
      nscnetp = 0
      scnetn = 0.0
      nscnetn = 0
      scnet = 0.0
      nscnet = 0
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       volzone = dv ! /gt(1,1,kz,imapz)
      if ( ichend(ichpos) .eq. 2 ) then
      if ( t4(ix,jy,kz) .gt. 0.5 ) then
!      if ( t0(ix,jy,kz) .gt. scthp ) then
      if ( t0(ix,jy,kz) .lt. scthn ) then
      scnetp = scnetp + t3(ix,jy,kz)*dv
      nscnetp = nscnetp + 1
      scnet = scnet + t3(ix,jy,kz)*dv
      nscnet = nscnet + 1
      volp = volp + dv
      end if
      end if
      end if
      if ( ichend(ichneg) .eq. 2 ) then
      if ( t4(ix,jy,kz) .lt. -0.5 ) then
!      if ( t0(ix,jy,kz) .lt. scthn ) then
      if ( t0(ix,jy,kz) .gt. scthp ) then
      scnetn = scnetn + t3(ix,jy,kz)*dv
      nscnetn = nscnetn + 1
      scnet = scnet + t3(ix,jy,kz)*dv
      nscnet = nscnet + 1
      voln = voln + dv
      end if
      end if
      end if
      end do
      end do
      end do
!
!      if ( ndebug .ge. 1 ) write(iunit,*) 'comp scnetp',scnetp,nscnetp
!      if ( ndebug .ge. 1 ) write(iunit,*) 'comp scnetn',scnetn,nscnetn
!      if ( ndebug .ge. 1 ) write(iunit,*) 'comp scnet',scnet,nscnet
!      if ( ndebug .ge. 1 ) write(iunit,*) 'comp volp',volp
!      if ( ndebug .ge. 1 ) write(iunit,*) 'comp voln',voln
!      if ( ndebug .ge. 1 ) write(iunit,*) 'COULOMB SCNETP',nstep,scnetp ! *volzone
!      if ( ndebug .ge. 1 ) write(iunit,*) 'COULOMB SCNETN',nstep,scnetn ! *volzone
!      if ( ndebug .ge. 1 ) write(iunit,*) 'COULOMB SCNET',nstep,scnet ! *volzone
!
!      lgtstp = 1
!      return

!
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  If scnetn or scnetp is too small to accomodate the charge in
!  voln or volp, respectively, then make a volume sufficient
!  to accomodate the charge so breakdown will not occur.
!  after that go to 59999 and repeat this  procedure above.
!  only one attempt will be required.
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      if ( irepeat .gt. 1 ) then
      write(iunit,*) 'failed lightning owing to failure to redistribute charge'
      lgtstp =1
      GOTO 2999 ! return
      end if

       idrop = 0
       scnetpi = 0.0
       scnetni = 0.0
       chgrtp = 0.0
       chgrtn = 0.0
       scdrop = 0.0
       sckeepp = 0.0
       sckeepn = 0.0
      IF ( isa .lt. 1 .and. iredis .eq. 1 ) THEN

       write(0,*) 'sll: defunct option; check code!'
       call commasmpi_abort()
      
!
!  accumululate scnet and nscnet
!

#ifndef MPI
! C$DOACROSS LOCAL(kz,jy,ix)
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(kz,jy,ix)
#endif
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
      t6(ix,jy,kz) = 0.0
      t7(ix,jy,kz) = 0.0
      if ( abs(t4(ix,jy,kz)) .gt. 0.5 ) then
       
     
!
!  If no mass is around, then zero the lightning charge,
!    except when the ambient charge (t0) is not zero, in which
!    case we apply t3, but only up to the amount in t0.
!
      IF ( t6(ix,jy,kz) .gt. 1.0e-12 ) THEN
        t7(ix,jy,kz) = 1/t6(ix,jy,kz)
!        chgrtn = chgrtn + Min( 0.0, t3(ix,jy,kz) )
!        chgrtp = chgrtp + Max( 0.0, t3(ix,jy,kz) )
      ELSE
       IF ( t0(ix,jy,kz) .lt. 0.0 .and. t3(ix,jy,kz) .gt. 0.0) THEN
        chg = Min ( Abs(t0(ix,jy,kz)), t3(ix,jy,kz) )
        scdrop = scdrop + (t3(ix,jy,kz) - chg)*dv ! /gt(ix,jy,kz,imapz)
        t3(ix,jy,kz) = chg
        sckeepp = sckeepp + chg*dv ! /gt(ix,jy,kz,imapz)
        t7(ix,jy,kz) = 1.0
       ELSEIF ( t0(ix,jy,kz) .gt. 0.0 .and. t3(ix,jy,kz) .lt. 0.0) THEN
         chg = Max ( -t0(ix,jy,kz), t3(ix,jy,kz) )
        scdrop = scdrop + (t3(ix,jy,kz) - chg)*dv ! /gt(ix,jy,kz,imapz)
        t3(ix,jy,kz) = chg
        sckeepn = sckeepn + chg*dv ! /gt(ix,jy,kz,imapz)
        t7(ix,jy,kz) = 1.0
       ELSE
        scdrop = scdrop + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
!        t3(ix,jy,kz) = 0.0
!        t4(ix,jy,kz) = 0.0
       ENDIF
        idrop = idrop + 1
        scnetpi = scnetpi + Max(0.0, t0(ix,jy,kz) )*dv ! /gt(ix,jy,kz,imapz)
        scnetni = scnetni + Min(0.0, t0(ix,jy,kz) )*dv ! /gt(ix,jy,kz,imapz)
      ENDIF
      end if
      end do
      end do
      end do
      
       scdropnet = scdropnet + scdrop
!      IF ( my_rank == 0 ) THEN
!      write(iunit,'(a,i6,2x,1pe12.5)')    &
!     &  'Points, net charge dropped outside cloud: ',    &
!     &   idrop,scdrop
!      write(iunit,'(a,2(1x,1pe12.5))')    &
!     &  'Net kept pos/neg lightning charge = ',    &
!     &   sckeepp,sckeepn
!      write(iunit,'(a,1pe12.5)')    &
!     &  'Cumulative net charge dropped outside cloud: ',    &
!     &   scdropnet
!      write(iunit,'(a,2(1x,1pe12.5))')    &
!     &  'Ambient net pos/neg charge outside cloud: ',    &
!     &   scnetpi, scnetni
!      ENDIF
       ENDIF !   ( isa .lt. 1 ) 

      
!      scalp = 1.0
!      scaln = 1.0
!      IF ( iic .eq. 1 ) THEN
!        IF ( Abs(chgrtn) .gt. chgrtp ) THEN
!          scaln = chgrtp/Abs(chgrtn)
!        ELESEIF ( Abs(chgrtn) .lt. chgrtp ) THEN
!          scalp = Abs(chgrtn)/chgrtp
!        ENDIF
!      ENDIF



      scnetpi = 0.0d0
      scnetni = 0.0d0
      if ( ndebug .ge. 2 ) write(iunit,*) 'set scnet,nscnet'
      scnetp = 0.0d0
      nscnetp = 0
      scnetn = 0.0d0
      nscnetn = 0
      scnet = 0.0d0
      nscnet = 0
      voln = 0.0d0
      volp = 0.0d0
      chgtot = 0.0d0
      chgpnts = 0.0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      chgavailiter(:) = 0.0
      voliter(:) = 0.0
      
      do kz = 1,nz-kd1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       volzone = dv ! /gt(1,1,kz,imapz)
        if ( t4(ix,jy,kz) .gt. 0.5 ) then
          scnetpi = scnetpi + t0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          scnetp = scnetp + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetp = nscnetp + 1
          scnet = scnet + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnet = nscnet + 1
          volp = volp + volzone
          n = Nint( t4(ix,jy,kz) )
          chgavailiter(n) = chgavailiter(n) + t3(ix,jy,kz)*dv
          voliter(n) = voliter(n) + volzone
        end if
        if ( t4(ix,jy,kz) .lt. -0.5 ) then
          scnetni = scnetni + t0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
          scnet = scnet + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnet = nscnet + 1
          voln = voln + volzone
          n = Nint( t4(ix,jy,kz) )
          chgavailiter(n) = chgavailiter(n) + t3(ix,jy,kz)*dv
          voliter(n) = voliter(n) + volzone
        end if
        
        IF ( ij_reduce(ix,jy) == 1 ) THEN ! count up charge in all points in the flash columns
          chgtot = chgtot + t0(ix,jy,kz)*dv
          IF ( kz == 1 ) chgpnts = chgpnts + 1
        ENDIF
        
      end do
      end do
      end do

#ifdef MPI
       n = 2*itermax+1
       CALL MPI_AllReduce(chgavailiter, chgavailitersum, n, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
       
       chgavailiter(-itermax:itermax) = chgavailitersum(-itermax:itermax)
       chgavailitersum(:) = 0.0

       CALL MPI_AllReduce(voliter, chgavailitersum, n, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
       voliter(-itermax:itermax) = chgavailitersum(-itermax:itermax)
       chgavailitersum(:) = 0.0


! find global integrated rate max
       mpitotindp(1)  = scnetpi
       mpitotindp(2)  = scnetp
       mpitotindp(3)  = nscnetp
       mpitotindp(4)  = volp

       mpitotindp(5)  = scnetni
       mpitotindp(6)  = scnetn
       mpitotindp(7)  = nscnetn
       mpitotindp(8)  = voln

       mpitotindp(9)  = scnet
       mpitotindp(10) = nscnet

       mpitotindp(11)  = idrop
       mpitotindp(12)  = scdrop
       mpitotindp(13)  = sckeepp
       mpitotindp(14)  = sckeepn
       mpitotindp(15)  = scdropnet
       mpitotindp(16) = chgtot
       mpitotindp(17) = chgpnts

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 17, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetpi = mpitotoutdp(1)
       scnetp  = mpitotoutdp(2)
       nscnetp = Nint( mpitotoutdp(3) )
       volp    = mpitotoutdp(4)

       scnetni = mpitotoutdp(5)
       scnetn  = mpitotoutdp(6)
       nscnetn = Nint( mpitotoutdp(7) )
       voln    = mpitotoutdp(8)

       scnet   = mpitotoutdp(9)
       nscnet  = Nint( mpitotoutdp(10) )

       idrop     = Nint( mpitotoutdp(11) )
       scdrop    = mpitotoutdp(12)
       sckeepp   = mpitotoutdp(13)
       sckeepn   = mpitotoutdp(14)
       scdropnet = mpitotoutdp(15)
       chgtot    = mpitotoutdp(16)
       chgpnts   = mpitotoutdp(17)

#endif

        chgavailitersum(itermax) = chgavailiter(itermax)
        chgavailitersum(-itermax) = chgavailiter(-itermax)
        volitersum( itermax) = voliter( itermax)
        volitersum(-itermax) = voliter(-itermax)
          IF ( my_rank == 0 ) THEN
!            write(iunit,*) 'n, chgavailitersum(n), chgavailitersum(-n)'
          ENDIF
        DO n = itermax-1,1,-1
          chgavailitersum( n) = chgavailitersum( n+1) + chgavailiter( n)
          chgavailitersum(-n) = chgavailitersum(-n-1) + chgavailiter(-n)

          volitersum( n) = volitersum( n+1) + voliter( n)
          volitersum(-n) = volitersum(-n-1) + voliter(-n)

          IF ( my_rank == 0 ) THEN
!            write(iunit,*) n, chgavailitersum(n), chgavailitersum(-n),volitersum(n),voliter(n)
          ENDIF
        ENDDO

!
!      voln = nscnetn*volzone
!      volp = nscnetp*volzone
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,i6,2x,1pe12.5)')    &
     &  'Points, net charge dropped outside cloud: ',    &
     &   idrop,scdrop
      write(iunit,'(a,2(1x,1pe12.5))')    &
     &  'Net kept pos/neg lightning charge = ',    &
     &   sckeepp,sckeepn
      write(iunit,'(a,1pe12.5)')    &
     &  'Cumulative net charge dropped outside cloud: ',    &
     &   scdropnet
!      write(iunit,'(a,2(1x,1pe12.5))')    &
!     &  'Ambient net pos/neg charge outside cloud: ',    &
!     &   scnetpi, scnetni
      ENDIF
      
      IF ( my_rank == 0 ) THEN
      if ( ndebug .ge. 1 ) then
      write(iunit,'(a,2(1x,1pe12.5))') 'Ambient charge scnetni,scnetpi ',    &
     &    scnetni,scnetpi
      end if
!
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &      'comp scnetp',scnetp,nscnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &     'comp scnetn',scnetn,nscnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &     'comp scnet',scnet,nscnet
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')     &
     &     'comp volp',volp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')     &
     &     'comp voln',voln
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')     &
     &     'COULOMB SCNETP',nstep,scnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')     &
     &     'COULOMB SCNETN',nstep,scnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3(3x,g20.9))')     &
     &     'COULOMB SCNET',nstep, scnet, chgtot, chgpnts
      ENDIF
!      IF ( isa .lt. 1 ) THEN

!      write(24,'(a,i6,3x,g20.9,2(1x,i6))') 'COULOMB SCNETP',nstep,
!     >                       scnetp*volzone, nscnetp,pchnl
!      write(24,'(a,i6,3x,g20.9,2(1x,i6))') 'COULOMB SCNETN',nstep,
!     >                       scnetn*volzone, nscnetn,nchnl
!      write(24,'(a,i6,3x,g20.9)') 'COULOMB SCNET',nstep,
!     >                            scnet*volzone
      
!      ENDIF
!
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  If scnetn or scnetp is to0 small to accomodate the charge in
!  voln or volp, respectively, then make a volume sufficient
!  to accomodate the charge so breakdown will not occur.
!  after that go to 59999 and repeat this  procedure above.
!  only one attempt will be required.
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      irepeat = 0
      if ( irepeat .gt. 1 ) then
      write(iunit,*) 'failed lightning owing to failure to redistribute charge'
      lgtstp =1
      GOTO 9898
      end if
!
!
!  compute scdelp (t3) [modified from scdel or t3]
!
!
      
      fracn = 1.0d0
      fracp = 1.0d0
      multp = 0.0d0
      multn = 0.0d0

!
! if ibal = 0 or 3 then want to use average maybe...
!
      if ( iic .eq. 1 ) then  
      if ( nscnet .gt. 0 ) then
!        sctarg = Abs( scnetn - scnet/Float(nscnet)*Float(nscnetn) )
        sctarg = Min( Abs(scnetn),Abs(scnetp) )
!        
! for ibal = 0 or 3 set sctarg to minimum of average or 2.0*the smaller.
!
        IF ( ibal .eq. 0 .or. ibal .ge. 3 ) THEN

         IF ( ibal .eq. 0 ) THEN
           sctarg = Min( 0.5d0*(Abs(scnetn) + Abs(scnetp)), 1.5d0*sctarg )
         ELSE
           sctarg = Min( 0.5d0*(Abs(scnetn) + Abs(scnetp)), 1.1d0*sctarg )
         ENDIF
         
         
          fracn = sctarg/Abs(scnetn)
          multn = 0.0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0

        ELSEIF ( ibal == 1 ) THEN
          
        if ( Abs(scnetn) .gt. Abs(scnetp) ) then
          IF ( chgtot < 0 ) THEN
          ! allow slightly unbalanced flash to account for net charge, but not more than the net charge
          sctarg = Min( Abs(scnetn),Min( 1.25*Abs(scnetp), Abs(scnetp) - chgtot)  )
          ELSE
          sctarg = Min( Abs(scnetn),Abs(scnetp) )
          ENDIF
          fracn = sctarg/Abs(scnetn)
          multn = 0.0d0
          fracp = 1.0d0
          multp = 0.0d0
        elseif ( Abs(scnetn) .lt. Abs(scnetp) ) then
          IF ( chgtot >  0 ) THEN
          sctarg = Min( Min( 1.25*Abs(scnetn), Abs(scnetn) + chgtot ) ,Abs(scnetp) )
          ELSE
          sctarg = Min( Abs(scnetn),Abs(scnetp) )
          ENDIF
          fracn = 1.0d0
          multn = 0.0d0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0d0
        end if
          
        ELSE ! ibal = 2
!        write(iunit,*)  'sctarg,scnetn,scnetp = ',sctarg,scnetn,scnetp
         
         IF ( .true. ) THEN
         tmp = Min(  Abs(scnetn), Abs(scnetp) )
         iterp = 1
         itern = -1
          IF ( Abs(scnetn) .gt. Abs(scnetp) ) THEN
            ! find index itern that just exceeds scnetp
            DO n = -itermax,-1
             ! IF ( my_rank == 0 ) write(iunit,*) 'n:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                itern = n
                EXIT
              ENDIF
            ENDDO
            scnetn = chgavailitersum(itern)
            voln = volitersum(itern)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetn = ',scnetn
              write(iunit,*) 'New voln = ',voln
            ENDIF
            WHERE ( (t4 > itern + 0.5) .and. t4 < -0.5 )
              t3 = 0.0
            END WHERE
          ELSE
            ! find index iterp that just exceeds scnetn
            DO n = itermax,1,-1
             ! IF ( my_rank == 0 ) write(iunit,*) 'p:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                iterp = n
                EXIT
              ENDIF
            ENDDO

            scnetp = chgavailitersum(iterp)
            volp = volitersum(iterp)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetp = ',scnetp
              write(iunit,*) 'New volp = ',volp
            ENDIF
            WHERE ( (t4 < iterp - 0.5) .and. t4 > 0.5 )
              t3 = 0.0
            END WHERE
          
          ENDIF
          
          IF ( my_rank == 0 ) THEN
            write(iunit,*) 'itern,iterp = ',itern,iterp
          ENDIF
          
         ENDIF
         
         ! find in
        
        if ( Abs(scnetn) .gt. Abs(scnetp) ) then
          fracn = sctarg/Abs(scnetn)
          multn = 0.0d0
          fracp = 1.0d0
          multp = 0.0d0
        elseif ( Abs(scnetn) .lt. Abs(scnetp) ) then
          fracn = 1.0d0
          multn = 0.0d0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0d0
        end if
        ENDIF
      end if
       IF ( my_rank == 0 ) THEN
         write(iunit,'(a,2(1x,f12.5))') 'fracn,multn',fracn,multn
         write(iunit,'(a,2(1x,f12.5))') 'fracp,multp',fracp,multp
       ENDIF
      end if

      if ( iic .eq. 1 ) then
      scnetnp = 0.0d0
      nscnetnp = 0
      scnetnn = 0.0d0
      nscnetnn = 0
      nscnetn = 0
      scnetn = 0.0d0
        if ( nscnet .gt. 0 ) then
          scnet = scnet/nscnet
        else 
          scnet = 0. 
        end if
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')    &
     &      'scnet=',scnet,nscnet
      if ( ndebug .ge. 2 ) write(iunit,'(a)') 'set scdelp--t3'

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      channtot = 0.0
      chanptot = 0.0
      chanzn(:) = 0.0
      chanzp(:) = 0.0
      
      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
      if ( t4(ix,jy,kz) .gt. 0.5 .and. t3(ix,jy,kz) /= 0.0) then
      t3(ix,jy,kz) = fracp*t3(ix,jy,kz) - multp*scnet
      scnetnp = scnetnp + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetnp = nscnetnp + 1
      scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetn = nscnetn + 1
      chanzp(kz) = chanzp(kz) + dx
      chanptot = chanptot + dx
      end if
      if ( t4(ix,jy,kz) .lt. -0.5 .and. t3(ix,jy,kz) /= 0.0) then
      t3(ix,jy,kz) = fracn*t3(ix,jy,kz) - multn*scnet
      scnetnn = scnetnn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetnn = nscnetnn + 1
      scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetn = nscnetn + 1
      chanzn(kz) = chanzn(kz) + dx
      channtot = channtot + dx
      end if
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = scnetnp
       mpitotindp(2)  = nscnetnp
       mpitotindp(3)  = scnetnn
       mpitotindp(4)  = nscnetnn

       mpitotindp(5)  = scnetn
       mpitotindp(6)  = nscnetn
       mpitotindp(7)  = channtot
       mpitotindp(8)  = chanptot

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 8, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetnp  = mpitotoutdp(1)
       nscnetnp = Nint( mpitotoutdp(2) )
       scnetnn  = mpitotoutdp(3) 
       nscnetnn = Nint( mpitotoutdp(4) )

       scnetn = mpitotoutdp(5)
       nscnetn = Nint( mpitotoutdp(6) )
       channtot = mpitotoutdp(7)
       chanptot = mpitotoutdp(8)
#endif


      IF ( my_rank == 0 ) THEN

      if ( nscnetn .gt. 0 ) then
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',    &
     &  scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',    &
     &  scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET',    &
     &  scnetn, nscnetn
      write(iunit,'(a,3(1x,g20.9))') 'CHANNEL LENGTH N,P,TOT',    &
     &  channtot, chanptot, channtot + chanptot
     
      ENDIF

!        IF ( isa .lt. 1 ) THEN

!          write(24,'(a,2(1x,g20.9))') 'scnetn, scnetp',
!     :    scnetnn*(volzone),scnetnp*(volzone)
      
!        ENDIF
      
      end if
      
      
!
! calculate dipole moment of the neutralizing charge
!  (change in dipole moment)
       qpx = 0.0d0
       qpy = 0.0d0
       qpz = 0.0d0
       qpxn = 0.0d0
       qpyn = 0.0d0
       qpzn = 0.0d0
       qpxp = 0.0d0
       qpyp = 0.0d0
       qpzp = 0.0d0
       qpxm = 0.0d0
       qpym = 0.0d0
       qpzm = 0.0d0
       qpxcn = 0.0d0
       qpycn = 0.0d0
       qpzcn = 0.0d0
       qpxcp = 0.0d0
       qpycp = 0.0d0
       qpzcp = 0.0d0
       qnetn = 0.0d0
       qnetp = 0.0d0
       npntn = 0
       npntp = 0
      do kz = 1,nz-1
      do jy = 1,jye ! ny-1
      do ix = 1,ixe ! nx-1
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       IF (t4(ix,jy,kz) .lt. -1.0e-18) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
        qpzn = qpzn +     &
     &         z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        qnetn = qnetn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF (t4(ix,jy,kz) .gt. 1.0e-18) THEN
        qpxp = qpxp + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyp = qpyp + gyt(jy,1)*t3(ix,jy,kz)*dv
        qpzp = qpzp +     &
     &         z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        qnetp = qnetp + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntp = npntp + 1
       END IF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = qpxn
       mpitotindp(2)  = qpyn
       mpitotindp(3)  = qpzn

       mpitotindp(4)  = qpxp
       mpitotindp(5)  = qpyp
       mpitotindp(6)  = qpzp

       mpitotindp(7)  = qnetn
       mpitotindp(8)  = npntn

       mpitotindp(9)  = qnetp
       mpitotindp(10) = npntp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 10, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       qpxn = mpitotoutdp(1)
       qpyn = mpitotoutdp(2)
       qpzn = mpitotoutdp(3)

       qpxp = mpitotoutdp(4)
       qpyp = mpitotoutdp(5)
       qpzp = mpitotoutdp(6)

       qnetn = mpitotoutdp(7)
       npntn = Nint( mpitotoutdp(8) )

       qnetp = mpitotoutdp(9)
       npntp = Nint( mpitotoutdp(10) )
#endif

! dipole moment components (C km)
        qpx = (qpxn + qpxp)*0.001d0
        qpy = (qpyn + qpyp)*0.001d0
        qpz = (qpzn + qpzp)*0.001d0
! location of midpoint of dipole (x,y,z in km)        
       IF ( qnetp .ne. 0.0 ) THEN
       qpxcn = -qpxn*0.001/(qnetp)
       qpycn = -qpyn*0.001/(qnetp)
       qpzcn = -qpzn*0.001/(qnetp)
       qpxcp = qpxp*0.001/(qnetp)
       qpycp = qpyp*0.001/(qnetp)
       qpzcp = qpzp*0.001/(qnetp)
        qpxm = (-qpxn + qpxp)*0.001/(2.0*qnetp)
        qpym = (-qpyn + qpyp)*0.001/(2.0*qnetp)
        qpzm = (-qpzn + qpzp)*0.001/(2.0*qnetp)
        ENDIF


      IF ( my_rank == 0 ) THEN
      
      write(iunit,'(a,2(1x,i7))') 'npntn, npntp ',npntn,npntp
      
      write(iunit,'(a,2(1x,1pe12.5))')     &
     &    'qnetn, qnetp',qnetn*(volzone),qnetp*(volzone)

      IF ( isa .lt. 1 ) THEN

!      write(24,'(a,2(1x,1pe12.5))') 
!     :    'qnetn, qnetp',qnetn*(volzone),qnetp*(volzone)
      
      ENDIF
      
        write(iunit,'(a,4(1x,1pe12.5))') 'dipole moment (px,py,pz) ',    &
     &       qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )

      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,1pe12.5))') 'dipole moment (px,py,pz) ',
!     :       qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',    &
     &    qpxm,qpym,qpzm   
      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',
!     :    qpxm,qpym,qpzm   
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'pos. drho center (x,y,z) ',    &
     &    qpxcp,qpycp,qpzcp   
        write(iunit,'(a,3(1x,f12.5))') 'neg. drho center (x,y,z) ',    &
     &    qpxcn,qpycn,qpzcn   

      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,f12.5))') 'pos. drho center (x,y,z) ',
!     :    qpxcp,qpycp,qpzcp   
!        write(24,'(a,3(1x,f12.5))') 'neg. drho center (x,y,z) ',
!     :    qpxcn,qpycn,qpzcn   
      ENDIF
      
      ENDIF ! my_rank == 0
        
      ELSE  ! cg or ca flash: calculate charge center
      
      IF ( ibal .eq. 2 .and. cgfr .ge. 1.0 ) THEN
        cgfrac = 1.0
        fracp = 1.0
        fracn = 1.0
        multp = 0.0
        multn = 0.0
      ELSE
        cgfrac = cgfr
      ENDIF
      
      IF ( (ibal .eq. 5)  .or.     &
     &    ( ibal .eq. 2 .and. cgfr .lt. 1.0 ) ) THEN
       IF ( my_rank == 0 ) THEN
        write(iunit,*) 'scnetp,scnetn',scnetp,scnetn
       ENDIF
      IF ( icgpn .eq. 1 ) THEN
! DISCHARGE IS POSITIVE, set scnetp <= cgfrac*scnetn
      fracn = 1.0
      multp = 0.0
      multn = 0.0
       IF ( Abs(scnetp) .gt. cgfrac*(Abs(scnetn))) THEN
          fracp = cgfrac*Abs(scnetn/scnetp)
       ELSE
          fracp = 1.0
       ENDIF
        
      ELSEIF ( icgpn .eq. -1 ) THEN
! DISCHARGE IS NEGATIVE, set scnetn <= cgfrac*scnetp
      fracp = 1.0
      multp = 0.0
      multn = 0.0
       IF ( Abs(scnetn) .gt. cgfrac*(Abs(scnetp))) THEN
          fracn = cgfrac*Abs(scnetp/scnetn)
       ELSE
          fracn = 1.0
       ENDIF
      
      ENDIF ! icgpn
        
        IF ( my_rank == 0 ) THEN
          write(iunit,'(a,(1x,f12.5))') 'fracn',fracn
          write(iunit,'(a,(1x,f12.5))') 'fracp',fracp
        ENDIF
      
      ENDIF ! (ibal .eq. 5 .or. (ibal .eq. 2 .and. cgfr .lt. 1.0) )

      scnetnp = 0.0
      nscnetnp = 0
      scnetnn = 0.0
      nscnetnn = 0
      scnetn = 0.0
      nscnetn = 0
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')    &
     &      'scnet=',scnet,nscnet
      if ( ndebug .ge. 2 ) write(iunit,'(a)') 'set scdelp--t3'

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
       do jy = 1,jye
        do ix = 1,ixe
        dv = dxx(ix)*dyy(jy)*dzz(kz)
        if ( t4(ix,jy,kz) .gt. 0.5 ) then
          t3(ix,jy,kz) = fracp*t3(ix,jy,kz)
          scnetnp = scnetnp + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetnp = nscnetnp + 1
          scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
        end if
        if ( t4(ix,jy,kz) .lt. -0.5 ) then
          t3(ix,jy,kz) = fracn*t3(ix,jy,kz)
          scnetnn = scnetnn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetnn = nscnetnn + 1
          scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
        end if
        end do
       end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = scnetnp
       mpitotindp(2)  = nscnetnp
       mpitotindp(3)  = scnetnn
       mpitotindp(4)  = nscnetnn

       mpitotindp(5)  = scnetn
       mpitotindp(6)  = nscnetn

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 6, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetnp  = mpitotoutdp(1)
       nscnetnp = Nint( mpitotoutdp(2) )
       scnetnn  = mpitotoutdp(3) 
       nscnetnn = Nint( mpitotoutdp(4) )

       scnetn = mpitotoutdp(5)
       nscnetn = Nint( mpitotoutdp(6) )
#endif

      if ( nscnetn .gt. 0 .and. my_rank == 0 ) then
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',    &
     &  scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',    &
     &  scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET ',    &
     &  scnetn, nscnetn
      ENDIF ! nscnetn
        
      
       qpx = 0.0
       qpy = 0.0
       qpz = 0.0
       qpxn = 0.0
       qpyn = 0.0
       qpzn = 0.0
       qpxcn = 0.0
       qpycn = 0.0
       qpzcn = 0.0
       qpxcp = 0.0
       qpycp = 0.0
       qpzcp = 0.0
       qnetn = 0.0
       qnetp = 0.0
       npntn = 0
       npntp = 0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
       do jy = 1,jye
        do ix = 1,ixe

       IF ((t4(ix,jy,kz)) .ne. 0.0) THEN

         dv = dxx(ix)*dyy(jy)*dzz(kz)

       IF ( (icgpn .eq. -1 .or. ibal .eq. 5 .or. ibal .eq. 2)    &
     &      .and. t3(ix,jy,kz) .gt. 0.0 ) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*t3(ix,jy,kz)

        qpzn = qpzn +     &
     &         z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

        qnetn = qnetn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF ( (icgpn .eq. 1 .or. ibal .eq. 5 .or. ibal .eq. 2)    &
     &      .and. t3(ix,jy,kz) .lt. 0.0 ) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*t3(ix,jy,kz)

        qpzn = qpzn +     &
     &         z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

        qnetn = qnetn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ENDIF

       ENDIF ! t4
      end do
      end do
      end do


#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = qpxn
       mpitotindp(2)  = qpyn
       mpitotindp(3)  = qpzn

!       mpitotindp(4)  = qpxp
!       mpitotindp(5)  = qpyp
!       mpitotindp(6)  = qpzp

       mpitotindp(7)  = qnetn
       mpitotindp(8)  = npntn

!       mpitotindp(9)  = qnetp
!       mpitotindp(10) = npntp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 8, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       qpxn = mpitotoutdp(1)
       qpyn = mpitotoutdp(2)
       qpzn = mpitotoutdp(3) 

!       qpxp = mpitotoutdp(4)
!       qpyp = mpitotoutdp(5)
!       qpzp = mpitotoutdp(6)

       qnetn = mpitotoutdp(7)
       npntn = Nint( mpitotoutdp(8) )

!       qnetp = mpitotoutdp(9)
!       npntp = Nint( mpitotoutdp(10) )
#endif



! dipole moment components (C km)
!        qpx = (qpxn + qpxp)*volzone*dx*0.001
!        qpy = (qpyn + qpyp)*volzone*dy*0.001
!        qpz = (qpzn + qpzp)*volzone*dz*0.001
! location of midpoint of dipole (x,y,z in km)
       IF ( qnetn .ne. 0.0 ) THEN
         qpxcn = qpxn*0.001/(qnetn)
         qpycn = qpyn*0.001/(qnetn)
         qpzcn = qpzn*0.001/(qnetn)
       ELSE
         qpxcn = 0.0
         qpycn = 0.0
         qpzcn = 0.0
       ENDIF
       
       IF ( my_rank == 0 ) THEN
       write(iunit,'(a,1pe12.5)')     &
     &   'Charge center of ',qnetn ! *volzone
       write(iunit,'(a,3(2x,f11.4))')     &
     &   'at x,y,z (km)', qpxcn,qpycn,qpzcn
       ENDIF

      IF ( isa .lt. 1 ) THEN
!       write(24,'(a,1pe12.5)') 
!     :   'Charge center of ',qnetn*volzone
!       write(24,'(a,3(2x,f11.4))') 
!     :   'at x,y,z (km)', qpxcn,qpycn,qpzcn
      ENDIF

!       write(24,'(7i6)') nstep, ilgt,iic,ica,icg,icgp,icgn
!       write(24,'(a,1pe12.5,a,3(2x,0pf11.4))') 
!     :   'Charge center of ',qnetn*volzone,' C is at x,y,z (km)',
!     :    qpxcn,qpycn,qpzcn
      
      end if  ! ( iic .eq. 1 ) 
!
!  sum total surface area exposed by hydrometeor species
!
      if ( ndebug .ge. 10 ) write(iunit,*) 'sum total particle area'
!
!  now redistribute the charge....
!
      if ( ndebug .ge. 2 ) write(iunit,*)'redistribute charge',lscb,lsceq
      if ( nscnet .gt. 0 ) then

      IF ( isa .lt. 1 ) THEN
!
!  redistribute charge based on surface area  NO LONGER AVAILABLE!
!  IF YOU WANT TO DO THIS, YOU WILL HAVE TO RECALCULATE THE SURFACE AREAS
!  HERE.
!
      if ( iredis .eq. 1 ) then
       
       write(0,*) 'iredis=1 is not available anymore!  Stopping!'
       call commasmpi_abort()
! C$DOACROSS LOCAL(ia,kz,jy,ix,chgrtp,chgrtn )
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(ia,kz,jy,ix,chgrtp,chgrtn)
      do ia = lscb,lsceq
       chgrtp = 0.0
       chgrtn = 0.0
      do kz = 1,nz-kd1
      do jy = 1,ny-jd1
      do ix = 1,nx-id1
       dv = dxx(ix)*dyy(jy)*dzz(kz)
      if ( abs(t4(ix,jy,kz)) .gt. 0.5 ) then
      if ( t6(ix,jy,kz) .gt. 1.e-12 ) then
!c      IF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
!        scal = t3(ix,jy,kz)*ab(ix,jy,kz,ia)*t7(ix,jy,kz)
!c      ELSE 
!c        scal = t3(ix,jy,kz)*ab(ix,jy,kz,ia)*t7(ix,jy,kz)
!c      ENDIF
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia) + scal
      chgrtn = chgrtn + Min( 0.0, scal )*dv ! /gt(ix,jy,kz,imapz)
      chgrtp = chgrtp + Max( 0.0, scal )*dv ! /gt(ix,jy,kz,imapz)
      end if
!      ELSEIF ( ia.eq.lsci .and. t7(ix,jy,kz) .gt. 0.5 ) THEN
!        scal = t3(ix,jy,kz)
!        an(ix,jy,kz,ia) = an(ix,jy,kz,ia) + scal
!        chgrtn = chgrtn + Min( 0.0, scal )/gt(ix,jy,kz,imapz)
!        chgrtp = chgrtp + Max( 0.0, scal )/gt(ix,jy,kz,imapz)
      end if
      end do ! ix
      end do ! jy
      end do ! kz
       chgrt(1,ia-lscb+1) = chgrtp
       chgrt(2,ia-lscb+1) = chgrtn
      end do ! ia
      
       chgrtp = 0.0
       chgrtn = 0.0

      write(iunit,'(a)') 'lightning charge breakdown by category: pos,neg'
      DO ia = lscb,lsceq
       chgrtp = chgrtp + chgrt(1,ia-lscb+1)
       chgrtn = chgrtn + chgrt(2,ia-lscb+1)
       write(iunit,'(i3,2(2x,1pe12.5))' ) ia,    &
     &       chgrt(1,ia-lscb+1),chgrt(2,ia-lscb+1)
      END DO ! ia
      
      write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',chgrtp,chgrtn
      
      end if
!
!  redistribute based on fraction
!
      if ( iredis .eq. 2 ) then
      do ia = lscb,lsceq
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      if ( abs(t4(ix,jy,kz)) .gt. 0.5 ) then
      an(ix,jy,kz,ia) = elgt1*an(ix,jy,kz,ia)
      end if
      end do
      end do
      end do
      end do
      end if
!
!  redistribute to ions
!
      IF ( iredis .eq. 3 ) THEN

      chgrtp = 0.0
      chgrtn = 0.0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-kd1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
      IF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni)-t3(ix,jy,kz)/ec
        chgrtn = chgrtn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSEIF ( t3(ix,jy,kz) .gt. 0.0 ) THEN
        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi)+t3(ix,jy,kz)/ec
        chgrtp = chgrtp + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ENDIF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgrtn
       mpitotindp(2)  = chgrtp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       chgrtn = mpitotoutdp(1)
       chgrtp  = mpitotoutdp(2)
#endif

      IF ( my_rank == 0 ) THEN
       write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',    &
     &  chgrtp,chgrtn
       ENDIF
      END IF
!
!  end nscnet
!
      if ( ndebug .ge. 1 ) write(iunit,*) 'redistribute charge done'
      end if
      
      ENDIF  ! ( isa .lt. 1 ) 
      ENDDO ! ilit

!
2999  continue
!
!  zero approprioate parts of the ab array
!
!      IF ( isa .ne. 1 ) THEN
!      if ( iliter .eq. nliter .or. lgtstp .eq. 1 ) then
! C$DOACROSS LOCAL(ia,kz,jy,ix)
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(ia,kz,jy,ix)
!      do ia = lscb,lsceq
!      do kz = 1,nz
!      do jy = 1,ny
!      do ix = 1,nx
!      ab(ix,jy,kz,ia) = 0.0
!      end do
!      end do
!      end do
!      end do
!      end if
!      
!      ENDIF !  ( isa .ne. 1 )
      

!    determine the clnox      
     IF ( lnox > 1 ) THEN
      total_nox = 0.0d0
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe     
!         tmp = 0.1*Min(1.0, Max(0.0,Abs(t3(ix,jy,kz)) - scth )/delq0) ! only make lnox for delta-q exceeding scth
         tmp = 0.1*Min(1.0, Max(0.0,Abs(t3(ix,jy,kz)) )/delq0) ! only make lnox for delta-q /= 0 and scale by delq0
         an(ix,jy,kz,lnox)=an(ix,jy,kz,lnox)+ tmp*lnox1(ix,jy,kz)
         dv=dx*dy*dzz(kz)
         total_nox = total_nox + tmp*nox(ix,jy,kz)
!         IF ( kz == 2 .and. tmp > 0.0 ) THEN
!           write(0,*) 'scdisch1: ',ix,jy,kz, scdisch0(ix,jy,kz)
!           write(0,*) 'scdisch2: ',scdisch(ix,jy,kz),tmp,lnox(ix,jy,kz)
!         ENDIF
      end do
      end do
      end do

       write(0,*) 'total_nox1 = ',total_nox

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = total_nox

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       total_nox = mpitotoutdp(1)
#endif

!       cnoxtot = total_nox
       write(0,*) 'total_nox = ',total_nox
       
      ENDIF


!
!
!  put the lightning position in the t4 array (holding the lightning
!  location)...
!
!
49999  continue
!
!
!
      if ( ilgt .eq. 1 )  then

! flush lightning summary file buffer
      IF ( isa .lt. 1 ) THEN
        CALL flush(24)  
      ENDIF
      
!       do in = 1,2
!       do it=1,nttimend(in) 
!        ic = ifix((trje(in,it,llx)+0.1)/dx+0.5) + 1
!        jc = ifix((trje(in,it,lly)+0.1)/dy+0.5) + 1
!        kc = ifix((trje(in,it,llz)+0.1)/dz+0.5) + 1
!        t4(ic,jc,kc) = 2.01
!       end do
!       end do 


!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

       DO kz = 1,nz-kd1
       DO jy = 1,jye
       DO ix = 1,ixe
       
!        IF (Abs(t2(ix,jy,kz)) .lt. 0.5 .and.     &
!     &      Abs(t4(ix,jy,kz)) .gt. 0.5 ) THEN
!            t4(ix,jy,kz) = Sign( 0.6,t4(ix,jy,kz) )
!        ELSE
!            t4(ix,jy,kz) = t2(ix,jy,kz)
!        END IF

         IF ( Abs(t4(ix,jy,kz)) > 0.5 .and. t3(ix,jy,kz) /= 0 ) THEN
           t2(ix,jy,kz) =  Sign(1.0,t4(ix,jy,kz)) 
         ELSE
           t2(ix,jy,kz)  = 0.0
         ENDIF
         IF ( kz == 1 .and.  t2(ix,jy,kz) /= 0.0 ) THEN
!          write(iunit,*) 'got to kz=1? STOP. ix,jy,kz,t2 = ',ix,jy,kz,t2(ix,jy,kz),t3(ix,jy,kz),t4(ix,jy,kz) 
!          call commasmpi_abort()
         ENDIF
         
       END DO
       END DO
       END DO
!      do it = 1,nttimend(1)
!      do is = idxlgt(1,it,1),idxlgt(1,it,1)
!      do js = idxlgt(1,it,2),idxlgt(1,it,2)
!      do ks = idxlgt(1,it,3),idxlgt(1,it,3)
!      t4(is,js,ks) = 2.01
!      end do
!      end do
!      end do
!      end do
!
!      do it = 1,nttimend(2)
!      do is = idxlgt(2,it,1),idxlgt(2,it,1)
!      do js = idxlgt(2,it,2),idxlgt(2,it,2)
!      do ks = idxlgt(2,it,3),idxlgt(2,it,3)
!      t4(is,js,ks) = 2.01
!      end do
!      end do
!      end do
!      end do
!
      end if
!
!
      end if  !   ( ilight .eq. 3 ) 
!
      if ( ndebug .ge. 1 ) write(iunit,*) 'done'
      if ( ilgt .ge. 1 ) loccur = Max(1, loccur)
!
!
!  END OF LIGHTNING PARAMETERIZATION VERSION 2
!
!
!
 9898 continue

!       deallocate ( pot2 )
!       deallocate ( rhs )
!       deallocate ( rho )

      return
      end 
!
!
! #####################################################################
! #####################################################################

!
      subroutine searchl    &
     &  (itype,itermax,icount,is,js,ks,id1,jd1,kd1    &
     &  ,nx,ny,nz,nxl,nyl,nzl,nor    &
     &  ,cons1,cons2,tthr,tflg,icgflg,klgtmin)
! 
      USE COMMASMPI_MODULE
      USE ELEC_MODULE, only: igslg0, jgslg0, kgslg0
      
      implicit none
!
#ifdef MPI
      INCLUDE 'mpif.h'
#else
      integer :: mytask = 0
#endif
      integer :: klgtmin
      integer ndebug
      parameter (ndebug=1)
      integer n,icheck,itype,icount,icount2,numflg
      integer :: ichecksum, icountsum
      integer, intent(IN) :: itermax
      integer icgflg
      integer i0,j0,k0
      integer i1,j1,k1
      integer i2,j2,k2
      integer is,js,ks
      integer ix,jy,kz
      integer id1,jd1,kd1
      integer nx,ny,nz,nmax
      integer nxl,nyl,nzl,nor
      integer nxs0,nys0,nzs0
      integer nxs1,nys1,nzs1
      integer nxs2,nys2,nzs2
      integer nxe0,nye0,nze0
      integer nxe1,nye1,nze1
      integer nxe2,nye2,nze2
      integer :: iunit = 91
 
      integer mynxs1,mynys1
      integer mynxe1,mynye1
      integer mynxe2,mynye2
      integer mynxs2,mynys2
!
      real cons1,cons2,xsval
      real tthr(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor) ! t2
      real tflg(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor) ! t4
!
      integer westward_tag,eastward_tag,southward_tag,northward_tag

!  -----------------------------------

      nmax = itermax ! max(nx,nz,ny) + 2
#ifdef MPI
! communication for tflg
      IF ( number_of_processes > 1 ) THEN 
        IF ( nproci > 1 ) THEN
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,tflg)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,tflg)
        ENDIF

        IF ( nprocj > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,tflg)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,tflg)
        ENDIF


        IF ( nproci > 1 ) THEN
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,tthr)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,tthr)
        ENDIF

        IF ( nprocj > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,tthr)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,tthr)
        ENDIF
      
      ENDIF
#endif

      icgflg = 0
!
      nzs0=ks
      nze0=ks
      nys0=js
      nye0=js
      nxs0=is
      nxe0=is
!      do k0 = nzs0,nze0
!      do j0 = nys0,nye0
!      do i0 = nxs0,nxe0
!      if ( itype .eq. 1 ) then
!      tflg(i0,j0,k0) =nmax
!      end if
!      if ( itype .eq. -1 ) then
!      tflg(i0,j0,k0) =-nmax
!      end if
!      end do
!      end do
!      end do
!
      do n = 1,nmax
      if ( ndebug .ge. 2 ) write(iunit,*) 'n,is,js,ks',n,is,js,ks
      icheck = 0
      numflg = 0
      xsval = (nmax-n)
      nxs1 = max(1,is-n)
      nxe1 = min(nxend,is+n)
      nys1 = max(1,js-n)
      nye1 = min(nyend,js+n)
      nzs1 = max(klgtmin,ks-n)
      nze1 = min(nz,ks+n)
      
! get limits relative to this tile:
      mynxs1 = nxs1 - igslg0
      mynxe1 = nxe1 - igslg0
      mynys1 = nys1 - jgslg0
      mynye1 = nye1 - jgslg0
!      mynxs1 = nxs1 - ixbeg + 1 ! igslg0
!      mynxe1 = nxe1 - ixbeg + 1 ! igslg0
!      mynys1 = nys1 - jybeg + 1 ! jgslg0
!      mynye1 = nye1 - jybeg + 1 ! jgslg0
      
!      write(6,*) 'mynxs1 etc = ',mynxs1,mynxe1,mynys1,mynye1
      
      IF ( ( mynxs1 <= nx .and. mynxe1 >= 1 ) .and. &
     &     ( mynys1 <= ny .and. mynye1 >= 1 ) ) THEN
      
      nxs1 = Max(mynxs1,1)
      nxe1 = Min(mynxe1,nx)
      nys1 = Max(mynys1,1)
      nye1 = Min(mynye1,ny)

!      write(6,*) 'nxs1 etc = ',nxs1,nxe1,nys1,nye1
      
      do k1 = nzs1,nze1
      do j1 = nys1,nye1
      do i1 = nxs1,nxe1
      nxs2 = max(0,i1-1)
      nxe2 = min(nx+1,i1+1)
      nys2 = max(0,j1-1)
      nye2 = min(ny+1,j1+1)
      nzs2 = max(klgtmin,k1-1)
      nze2 = min(nz,k1+1)
      do k2 = nzs2,nze2
      do j2 = nys2,nye2
      do i2 = nxs2,nxe2
!
!  itype .eq. 1
!
      if ( itype .eq. 1 ) then
      if ( tthr(i1,j1,k1) .gt. 0.5 ) then
      numflg = numflg + 1
      if ( tflg(i2,j2,k2) .gt. 0.5 ) then 
      if(tflg(i1,j1,k1).lt.0.5) then
      icheck = icheck+1
      tflg(i1,j1,k1) = xsval
     
      end if
      end if
      end if
      end if
!
!  itype .eq. -1
!
      if ( itype .eq. -1 ) then
      if ( tthr(i1,j1,k1) .lt. -0.5 ) then
      numflg = numflg + 1
      if ( tflg(i2,j2,k2) .lt. -0.5 ) then 
      if(tflg(i1,j1,k1).gt.-0.5) then
      icheck = icheck+1
      tflg(i1,j1,k1) = -xsval
      
      end if
      end if
      end if
      end if
!
      end do
      end do
      end do
!
      end do
      end do
      end do
      
      ENDIF ! limits within tile
!

#ifdef MPI
! communication for tflg
      IF ( number_of_processes > 1 ) THEN 
        IF ( nproci > 1 ) THEN
!        write(0,*) 'send_west, ntask_x = ',ntasks_x
!        CALL sendrecv_westward(nx,ny,nz,1,1,0,1,1,  &
!        CALL sendrecv_westward(nxslm,nyslm,nzslm,ghosti1,0,0,ghosti1,1,  &
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,tflg)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,tflg)
        ENDIF

        IF ( nprocj > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,tflg)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,tflg)
        ENDIF
      
      ENDIF
#endif
#ifdef MPI
!
!
! MPI: reduce icheck as SUM
!
!      myicheck = icheck
      IF ( number_of_processes > 1 ) THEN
        CALL MPI_AllReduce(icheck, ichecksum, 1, MPI_INTEGER, MPI_SUM, my_comm, mpi_error_code)
        icheck = ichecksum
      ENDIF
#endif

      if ( ndebug .ge. 2 ) then
      write(iunit,*) 'search iteration at n,icheck',n,icheck
      write(iunit,*) 'icount,',icount
      write(iunit,*)  nzs1,nze1
      write(iunit,*)  nys1,nye1
      write(iunit,*)  nxs1,nxe1

      IF ( nzs1 == 1 ) THEN
!       write(iunit,*) 'got to kz=1? STOP. is,js,ks = ',is,js,ks
!       call commasmpi_abort()
      ENDIF
      end if
!
      if ( icheck .eq. 0 ) EXIT ! go to 777
!
      end do
777   continue
!
      if ( ndebug .ge. 3 ) then
      write(iunit,*) 'done'
      do jy = 1,ny
      write(iunit,*) 'jy',jy
      do kz = nz,1,-1
      write(iunit,101) (tflg(ix,jy,kz),ix=is-10,is+10)
101   format(1x,22f3.0)
      end do
      end do
      end if
!
      icount = 0
!
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
      if ( itype .eq. 1 ) then
      if ( tflg(ix,jy,kz) .gt. 0.5 ) then
!      tflg(ix,jy,kz) = 1.01
      icount = icount + 1
      end if
      end if
      if ( itype .eq. -1 ) then
      if ( tflg(ix,jy,kz) .lt. -0.5 ) then
!      tflg(ix,jy,kz) = -1.01
      icount = icount + 1
      end if
      end if
      IF ( kz == 1 .and. tflg(ix,jy,kz) /= 0.0 .and. Int(Sign(1., tflg(ix,jy,kz))) == itype ) icgflg = icgflg + 1
      end do
      end do
      end do

!        IF ( Any(tflg(:,:,1)  /= 0 ) ) icgflg = 1
      icount2 = icount
#ifdef MPI
!
! MPI: reduce icount as SUM
!
      IF ( number_of_processes > 1 ) THEN
        CALL MPI_AllReduce(icount, icountsum, 1, MPI_INTEGER, MPI_SUM, my_comm, mpi_error_code)
        icount = icountsum
      ENDIF
#endif

      if ( ndebug .ge. 1 ) then
      write(iunit,*) 'icount',icount,icount2,icheck
      end if
!
      return
      end
!
!
!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######             INTEGER FUNCTION FIND_INDEX              ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     This function returns the array index (here, the value returned by
!     find_index is designated as i) such that x is between xa(i) and xa(i+1).
!     If x is less than xa(1), then i=0 is returned.  If x is greater than
!     xa(n), then i=n is returned.  It is assumed that the values of
!     xa increase monotonically with increasing i.
!
!############################################################################
!
!     Author:  David Dowell (based on "locate" algorithm in Numerical Recipes)
!
!     Creation Date:  17 November 2004
!
!############################################################################

      integer function find_indexlgt(x, xa, n)

      implicit none

      integer n                          ! array size
      real xa(n)                         ! array of locations
      real x                             ! location of interest
      integer il, im, iu                 ! lower and upper limits, and midpoint


      il = 0
      iu = n+1
      
      IF ( x .gt. xa(n) ) THEN
        il = n+2
      ELSEIF ( x .lt. xa(1) ) THEN
        il = -1
      ELSE

 10   if ((iu-il).gt.1) then
        im=(il+iu)/2
        if (x.ge.xa(im)) then
          il=im
        else
          iu=im
        endif
        go to 10
      endif
      
      ENDIF

!      if ( (il.eq.0) .and. (n.gt.1) .and.
!           ( (xa(1)-x) .lt. 0.001*(xa(2)-xa(1)) ) ) then
!        il = 1
!      endif
!      if ( (il.eq.n) .and. (n.gt.1) .and.
!           ( (x-xa(n)) .lt. 0.001*(xa(n)-xa(n-1)) ) ) then
!        il = n-1
!      endif

      find_indexlgt = il

      return
      end
