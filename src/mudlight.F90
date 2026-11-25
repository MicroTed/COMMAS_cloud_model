!CVD$F SKIP
! NOTE: fracn fix is enabled!!!!
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
      subroutine mudlight      &
         (time_real,time,nstep,iliter,lgtstp,numlgt,loccur      &
         ,nor,nx,ny,nz      &
         ,na      &
         ,id1,jd1,kd1,istag,jstag,kstag,icldtop,zcldtop      &
         ,icgyn      &
         ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
         ,elec,an,gx,gy,gz,gxt,gyt,gzt      &
         ,dxx,dyy,dzz      &
         ,ixst,jyst,kzst      &
         ,numic,numcgn,numcgp,outname,nkzmn,      &
          nnx,nny,nnz,      &
          nnxs,nnys,nnzs,      &
          iixp,jjyq,kkzr,iiex,jjey,kkez,      &
          nbw,nbe,nbs,nbn,nbz2,      &
          llwork, istretch, dzmax,      &
          nxslm,nyslm,nzslm,      &  ! nb,nxb,nyb,nzbt,nzb, 
          iunit, bcx, bcy, db)
     
       USE GRID_MODULE
       USE INDEX_MODULE
       USE ELEC_MODULE
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
      implicit none

#ifdef MPI
      INCLUDE "mpif.h"
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
!      include 'sam.files.h'
      
      integer, parameter    :: ng1 = 1
      integer iixp,jjyq,kkzr,iiex,jjey,kkez
      integer*8 llwork
      integer nnx,nny,nnz,nnzt
      integer nnxs(2),nnys(2),nnzs(2)
      integer nnytmp
      integer nbw,nbe,nbs,nbn,nbz2
      integer nbntmp
      
      integer bcx, bcy  ! boundary condition flags
      
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

      double precision chgneg,chgpos ! total net negative/positive charge in storm (not including corona)
      
      integer iretry   ! counter for nretry
      integer iredo
      integer kzbrkd,iter,iread
      integer ix,jy,kz,kz1,loccur,i,j
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
      integer lgtgo,imkseed,ichoose
      data imkseed /0/
      save imkseed
      integer islgt,jslgt,kslgt
      real    t0mx,t0mn,scramt,cndamt,times
      integer in,it,ic,jc,kc
      integer  idebug,idtsign, isign1,isign2
      real    eps,facx,facy,facz
      parameter (eps = 1.0e-7)
      real    telec,tcharge,dtcomp
      integer i1,j1,k1,ichneu
      integer ichpos,ichneg
      real    scthend1,scthend2,scthendp,scthendn
      real    scthp,scthn,potthrp, potthrn
      integer isendn,jsendn,ksendn,isendp,jsendp,ksendp
      integer iic,ica,icg,icgp,icgn,i2,j2,k2
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
      real rtest,chanlen,rchan,tanth
      !  rtest is the distance from a gridpoint to the channel test point
      !  chanlen is the length along the channel of the test point
      !  rchan = rbase + chanlen*Tan(theta)
      !  tanth = Tan(theta)
      integer nchadd(2)  !  number of points added along channel
      integer isignt
      integer ichcnt
      real chcnt
      real delx,dely,delz  !  relative coords of channel point from corner
      
! variables from 'inlightdbm'
      real scth, fscth
      real  scthmn, scthmx, fprcnt
      integer  niter, icountth, nretry
      real rbase,theta
      real ecrit
!
      
      save iread
      
       integer ninit,ninitlocal  ! number of possible initiation points
       integer npick
       parameter(npick = 10000)
       integer ipick(npick,3)
       integer idum

!
!  nttim:  maximum number of intervals in a lightning channel
!  nenum:  number of ends to a channel from a common location
!          (a positive and a negative end)
!  ninfo:  number of parameters in the trje array 
!
      parameter  (nttim=4000,nenum=2,ninfo=11)
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
      parameter  (ndebug=0,itest=0)
      parameter  (lgtinit=1)
!      ,ibrkd=2)
      parameter  (ieout=1)
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
      integer    nttimend(nenum),idxlgt(nenum,nttim,3)
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
      real emax, emag, zlev(1000), ebrkd(1000), ebrkdp(1000)
      real emaxkz(1000),exmaxkz(1000),eymaxkz(1000),ezmaxkz(1000)
      save zlev,ebrkd,ebrkdp
      real estop(1000)
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

      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      
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
      integer,allocatable :: imask(:,:,:,:)
      
      real,    allocatable :: xy_init(:,:)
      integer, allocatable :: ij_reduce(:,:)
      integer, allocatable, save :: ijmask(:,:)
      integer, save :: nxm, nym
      integer              :: nr
      
      real  :: x,y,r


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: nz1
      integer :: itype,imxmn
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
      
      integer, save :: boxmginit = 0

      integer, parameter :: ntot = 50
      real  mpitotin(ntot), mpitotout(ntot)
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      
! #####################################################################
!   Begin Executable Code
! #####################################################################


      ldumpsp = .false.
      
!
!  read in variables
!
      if ( iread .eq. 0 ) then
!
      iread = 1
!      
!      numic = 0
!      numcgp = 0
!      numcgn = 0
!      
!      IF ( isa .ge. 1 ) THEN
!        ifile = 5
!      ELSE
!        ifile = 47
!        open(ifile,file='inlightdbm',form='formatted')
!      ENDIF
!c
!      write(iunit,*) 'reading inlightdbm'
!c
!c  read elgtthx:  e magnitude maximum threshold for breakdown
!c
!      read(ifile,*) elgtthx
!      write(iunit,*) 'elgtthx = ',elgtthx
!c
!c  read elgtthn:  e magnitude mminimum threshold to stop
!c  lightning channels
!c
!      read(ifile,*) elgtthn
!      write(iunit,*) 'elgtthn = ',elgtthn
!c
!c  read elgtdel:  e magnitude delta to determine area for breakdown
!c
!      read(ifile,*) elgtdel
!      write(iunit,*) 'elgtdel = ',elgtdel
!c
!c  read elgtfdel:  e magnitude fraction of delta to 
!c  determine variable e magnitude [f(z)] for breakdown
!c
!      read(ifile,*) elgtfdel
!      write(iunit,*) 'elgtfdel = ',elgtfdel
!c
!c  read scth:  space charge magnitude threshold maximum for breakdown
!c
!      read(ifile,*) scth
!      write(iunit,*) 'scth = ',scth
!c
!c  read scth:  fraction multiplier for space charge threshold maximum
!c  maximum for lightning volume determination
!c
!      read(ifile,*) fscth
!      write(iunit,*) 'fscth = ',fscth
!c
!c  read  scthmn:  space charge magnitude minimum:  NOT USED
!c
!      read(ifile,*) scthmn
!      write(iunit,*) 'scthmn = ',scthmn
!c
!c  read  scthmp:  space charge magnitude maximum:  NOT USED
!c
!      read(ifile,*) scthmx
!      write(iunit,*) 'scthmx = ',scthmx
!c
!c  fractional percen used in space charge redistribution in lightning
!c  volume
!c
!      read(ifile,*) fprcnt
!      write(iunit,*) 'fprcnt = ',fprcnt
!c
!c  niter:  iteration to ensure scthp, scthn, potthrp, and potthrn
!c  are reasonable for determining lightning volume...also used
!c  to ensure there are enough points in each volume (+/-) where
!c  lightning occurs.
!c
!      read(ifile,*) niter
!      write(iunit,*) 'niter = ',niter
!c
!c  icountth:  minimum number of points required in either positive
!c  or negative charge are where lightning occurs...sort of a wag
!c  type parameter...if the number icountth of points is not met
!c  another attempt is made at determining the lightning volume
!c  reducing scth.  After niter times the attempt is aborted.
!c
!      read(ifile,*) icountth
!      write(iunit,*) 'icountth = ',icountth
!c
!c
!      read(ifile,*) nretry
!      write(iunit,*) 'nretry = ',nretry 
!c
!c
!      read(ifile,*) rbase
!      write(iunit,*) 'rbase = ',rbase 
!c
!c
!      read(ifile,*) theta
!      write(iunit,*) 'theta = ',theta 
!      tanth = Tan(theta*4.*atan(1.)/180.0)
!c
!c
!      read(ifile,*) ibrkd
!      write(iunit,*) 'ibrkd = ',ibrkd
!c
!      read(ifile,*) ieint,eint
!      write(6,*) 'ieint,eint = ',ieint,eint
!      read(ifile,*) feint
!      write(6,*) 'fraction of breakdown threshold for eint: ',feint
!      read(ifile,*) fecrit, fefac(1), fefac(2)
!      write(6,*) 'fecrit, fefac(1), fefac(2) = ',
!     :            fecrit, fefac(1), fefac(2)
!      read(ifile,*) iseed
!      write(6,*) 'iseed = ',iseed
!c      read(ifile,*)  intfg
!      read(ifile,'(a)',end=101) line
!        read(line,*,iostat=istat1)  intfg, ibranch
!        IF ( istat1 .ne. 0 ) THEN
!          write(6,*) 'read error, iostat = ',istat1
!          write(6,*) 'setting ibranch ',0
!          ibranch = 0
!        ENDIF
!
!      write(6,*) 'intfg,ibranch = ',intfg,ibranch
!      read(ifile,*) mdel
!      write(6,*) 'mdel = ',mdel
      IF ( mdel .lt. 10  ) THEN
        write (iunit,*) 'Must have mdel at least = 10, resetting mdel=10'
        mdel = 10
      ENDIF
      IF ( mdel .lt. 15 .and. dslight .lt. 251.00 ) THEN
        write(iunit,*) 'Must have mdel at least =15 for dslight .le. 250'
        mdel = 15
      ENDIF

      write(iunit,*) 'mdel = ',mdel
      
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
      ebrkdp(kz) = 284.e3*db(kz1)/rho00 ! based on Dwyer 2003, although probably not the actual initiation mechanism

      IF (ibrkd .eq. 2) THEN
        ebrkd(kz)= Min( ebrkd(kz), einitmax )
      ELSEIF ( ibrkd .eq. 3 ) THEN
        ebrkd(kz)= Min( ebrkd(kz), einitmax )
      ELSEIF ( ibrkd .eq. 4 ) THEN
        ebrkd(kz)= Min( ebrkdp(kz), einitmax )
      ELSEIF ( ibrkd .eq. 5 ) THEN
        ebrkd(kz)= Min( ebrkdp(kz), einitmax )
      ENDIF
      
      IF ( ebrkd(kz) .ne. 0.0 ) THEN
        dv = ebrkdp(kz)/ebrkd(kz)
      ELSE
        dv = 0.0
      ENDIF
      
      write(iunit,'(i3,2x,f9.4,2x,1pe13.5,2x,1pe13.5,0pf9.3)')       &
                 kz,1.0e-3*zlev(kz),ebrkd(kz), ebrkdp(kz), dv ! ebrkdp(kz)/ebrkd(kz)
      end do
!
       ELSE 

        write(iunit,*) 'Using uniform initiation threshold of ',elgtthx
        DO kz=1,nz
          ebrkd(kz) = elgtthx
          elgtfdel = (elgtthx - elgtdel)/elgtthx
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
      
      end if  ! (iread .eq. 0)
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
!
!  LIGHTNING DISCHARGE:  SCHEME 1.....
!
!
!
!  This scheme (ilight=1) is the simplist of all.  If |E| => elgtt1 
!  (400 kV/m; REF) then the space charge is decreased by some percent
!  elgtp1 (30%:  Rawlins, 1981, QJRMS).
!  
      if ( ilight .eq. 1 ) then

      t2(:,:,:) = 0.0 ! for number of channel points
      t4(:,:,:) = 0.0 ! for lightning-deposited charge

      IF ( lightrad < 0 ) THEN 
        nr = 0
      ELSE
        nr = NInt( lightrad/dx )
      ENDIF
      
      allocate( xy_init(-nr+1:nx+nr,-nr+1:ny+nr))
      allocate( ij_reduce(-nr+1:nx+nr,-nr+1:ny+nr))
      
      xy_init(:,:) = 0.0

      IF ( lightrad < 0 ) THEN 
        ij_reduce(:,:) = 1
      ELSE
        ij_reduce(:,:) = 0
      ENDIF
!
!  find largest etot and location, but first zero etot and location
!  indices first
!
      etot = 0.0
      ietotx = 0
      ietoty = 0
      ietotz = 0
      ninit = 0
!
      ixe = nx
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      
      jye = ny
      IF ( myprocj == nprocj ) jye = jyend - jybeg

      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe
!
      if ( elec(iemag)%flt3d(ix,jy,kz) .gt. etot ) then
      ietotx = ix
      ietoty = jy
      ietotz = kz
      ixst = ix
      jyst = jy
      kzst = kz
      etot = elec(iemag)%flt3d(ix,jy,kz)
      end if

      
      if ( (ibrkd == 1 .and. elec(iemag)%flt3d(ix,jy,kz) .gt. elght1) .or.       &
            (ibrkd > 1 .and. elec(iemag)%flt3d(ix,jy,kz) .gt. ebrkd(kz) ) ) then
        xy_init(ix,jy) = 1.0
        ninit = ninit + 1
        elec(ieinit)%flt3d(ix,jy,kz) = elec(ieinit)%flt3d(ix,jy,kz) + 1
      ENDIF

!
      end do
      end do
      end do
      
      ninitlocal = ninit
!
#ifdef MPI

! find global integrated rate max
       mpitotindp(1)  = etot
       mpitotindp(2)  = ninit

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)

       etot = mpitotoutdp(1)
       ninit = Nint( mpitotoutdp(2) )



! find global integrated rate max
       mpitotindp(1)  = etot

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, my_comm, mpi_error_code)

       etot = mpitotoutdp(1)  

      IF ( number_of_processes .gt. 1 .and. nr > 0 .and. ninit > 0 ) THEN
    

        
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,1,nr,nr,0,nr,1,  &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,xy_init)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,1,nr,nr,0,nr,1,  &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,xy_init)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,xy_init)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,xy_init)
    
      ENDIF
#endif


      IF ( my_rank == 0 ) THEN
        write(iunit,*) 'ETOT,ninit', ETOT,ninit
        write(iunit,*) 'Location and value of max(|E|)'
        write(iunit,*) 'IETOTX, IETOTY, IETOTZ, ETOT =: ',ietotx,ietoty, ietotz, etot
      ENDIF
!
!  decrease space charge by percentage elgtp1
!
      chgpos = 0
      chgneg = 0
      if ( ninit > 0 ) then
      nlgt1 = nlgt1 + 1

      IF ( nr > 0 ) THEN ! create square region around each init point
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
      
      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1
      
      do jy = 1,jye
      do ix = 1,ixe
      
      chg = 0.0d0

       IF ( ij_reduce(ix,jy) > 0 ) THEN

        do kz = 1,nz-1

!      DO ia = lscb,lsceq
!        chg = chg +  an(ix,jy,kz,ia) 
!        an(ix,jy,kz,ia) = an(ix,jy,kz,ia)*elgt1
!      ENDDO

       dv = dxx(ix)*dyy(jy)*dzz(kz)

!
! t1 is the non-ion charge
!
!       IF ( chg .ne. 0.0 ) THEN
!       t1(ix,jy,kz) = t0(ix,jy,kz) 
!       ELSE
!       t1(ix,jy,kz) = ! t0(ix,jy,kz) 
!     :  chg + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  
!       ENDIF
       
! add up available positive and negative charge excess (greater density than chgthr) within the cylider volumes
      IF ( t0(ix,jy,kz) .gt. chgthr ) THEN
        t1(ix,jy,kz) = t0(ix,jy,kz) - chgthr
        chgpos = chgpos + t1(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        t2(ix,jy,kz) = -1
      ELSEIF ( t0(ix,jy,kz) .lt. -chgthr ) THEN
        t1(ix,jy,kz) = t0(ix,jy,kz) + chgthr
        chgneg = chgneg + t1(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        t2(ix,jy,kz) = 1
      END IF
      
      end do

      ENDIF

      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
#endif

      
       chgrtp = 0.0
       chgrtn = 0.0
       x = 0.0
       y = 0.0
      
      chg = (1.0 - elgt1)*Max(chgpos, Abs(chgneg))
      IF ( chg .gt. chgpos ) chg = chgpos
      IF ( chg .gt. Abs(chgneg) ) chg = Abs(chgneg)
      
      IF ( Abs(chgneg) > 1.e-3 ) x = chg/(ec*Abs(chgneg))
      IF ( Abs(chgpos) > 1.e-3 ) y = chg/(ec*Abs(chgpos))
      
      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe

      IF ( ij_reduce(ix,jy) > 0 ) THEN

       dv = dxx(ix)*dyy(jy)*dzz(kz)
      
!      chg = chg*(1.0 - elgt1)
      IF ( t1(ix,jy,kz) .lt. 0.0 ) THEN ! local chg is neg. so add pos. ions
        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) - t1(ix,jy,kz)*x
        chgrtp = chgrtp + t1(ix,jy,kz)*x*ec*dv
      ELSEIF ( t1(ix,jy,kz) .gt. 0.0 ) THEN ! local chg is pos. so add neg. ions
        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) + t1(ix,jy,kz)*y
        chgrtn = chgrtn + t1(ix,jy,kz)*y*ec*dv
      ENDIF
      
      ENDIF

      end do
      end do
      end do

      loccur = 1

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgrtp
       mpitotindp(2)  = chgrtn

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       chgrtp = mpitotoutdp(1)  
       chgrtn = mpitotoutdp(2)  
#endif
      
!
!  write out that there is  a lightning discharge
!
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '______________________________'
      write(iunit,*) 'LIGHTNING DISCHARGE:  SCHEME 1'
      write(iunit,*) 'TIME STEP IS: ',time_real
      write(iunit,*) 'TOTAL NUMBER OF DISCHARGES IS: ',nlgt1
      write(iunit,*) 'chgpos,chgneg,chg : ',chgpos,chgneg,chg
      write(iunit,*) 'chgrtp,chgrtn,x,y = ',chgrtp,chgrtn,x,y
      write(iunit,*) '______________________________'
      ENDIF
!
      ELSE
       lgtstp = 1
      loccur = 0
      ENDIF
!
!  end of case 1
!
      deallocate( xy_init   )
      deallocate( ij_reduce )

      end if
!
!
!
!  END OF LIGHTNING PARAMETERIZATION VERSION 1
!
!
!
!
!  LIGHTNING PARAMETERIZATION 2
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
      if ( ilight .eq. 2 ) then
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
      ENDIF

#endif

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
      
!      IF ( myproci == 1 ) ixb = 2
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      
      jyb = 1
      jye = jtile
      
!      IF ( myprocj == 1 ) jyb = 2
      IF ( myprocj == nprocj ) jye = jyend - jybeg
      IF ( ny .le. 2 ) THEN
        jyb = 1
        jye = 1
      ENDIF

      do kz = 1,nz-kd1
!
      emaxkz(kz) = 0
      nemaxx(kz) = 1
      nemaxy(kz) = 1
      nemaxz(kz) = kz
!
      do jy = jyb,jye ! Min(ny-1,2),Max(1,ny-2*jd1)
      do ix = ixb,ixe ! 2,nx-2*id1
!
      emag = elec(iemag)%flt3d(ix,jy,kz)
!      IF ( iemagstag > 1 ) THEN 
!        emag = Max( emag, elec(iemagstag)%flt3d(ix,jy,kz) )
!      ENDIF
      if ( emag .gt. emaxkz(kz) ) then
      nemaxx(kz) = ix
      nemaxy(kz) = jy
      nemaxz(kz) = kz
      emaxkz(kz) = emag ! elec(iemag)%flt3d(ix,jy,kz)
      exmaxkz(kz) = elec(iex)%flt3d(ix,jy,kz)
      eymaxkz(kz) = elec(iey)%flt3d(ix,jy,kz)
      ezmaxkz(kz) = elec(iez)%flt3d(ix,jy,kz)
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
      DO kz = 1,nz1

       max_send_buffer(1,kz) = mpimxmn(kz,4)
       max_send_buffer(2,kz) = real(my_rank)

      ENDDO
    
      mxmn_size = nz1

! find out the max emag at each level, and which processor it is on
      CALL MPI_Allreduce(max_send_buffer, max_recv_buffer, mxmn_size, MPI_2REAL,        &
                         MPI_MAXLOC, my_comm, mpi_error_code)

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

      CALL MPI_Gather(exyz_send_buffer,exyz_send_size,MPI_REAL,      &
                      exyz_recv_buffer,exyz_recv_size,MPI_REAL,       &
                      0,my_comm,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiexyz0 = reshape(exyz_recv_buffer,(/nz1,3,number_of_processes/))

           DO kz = 1,nz1
       
            irank1 = nint(max_recv_buffer(2,kz))

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

      CALL MPI_Gather(loc_send_buffer,loc_send_size,MPI_INTEGER,      &
                      loc_recv_buffer,loc_recv_size,MPI_INTEGER,       &
                      0,my_comm,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiloc0 = reshape(loc_recv_buffer,(/nz1,2,number_of_processes/))

!          DO itype = 1,2 ! type is really the ix or jy index
           DO kz = 1,nz1
       
            irank1 = nint(max_recv_buffer(2,kz))

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

      write(iunit,'(a,3(i4,1x),2(g15.3,1x))')       &
         'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',iemaxx,      &
         iemaxy, iemaxz, emax, emaxkz(iemaxz) ! elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)
!
      write(iunit,'(a)')       &
        'Z, NEMAXX, NEMAXY, NEMAXZ, EMAXKZ, EZ, EH, yesno =: '
     
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
      if ( ibrkd .ge. 1 ) then
      
      lgtstp = 0
      lgtgo = 0
      do kz = nz-kd1,1,-1
      if ( emaxkz(kz) .lt. efracinitwire*overvolt*ebrkd(kz) .or. kz >= icldtop ) then ! don't count points above cloud top
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
      write(iunit,'(3(i4,1x),1x,3(1x,1pe12.4),1x,a,1x,0pf7.0)')       &
          nemaxx(kz),nemaxy(kz), nemaxz(kz),      &
          emaxkz(kz), ezmaxkz(kz),      &
         Sqrt(exmaxkz(kz)**2 + eymaxkz(kz)**2 ),yesno,      &
          zlev(kz)
!      write(iunit,'(3i4,2x,3(1x,1pe12.4),1x,a,1x,0pf7.0)') 
!     >   nemaxx(kz),nemaxy(kz), nemaxz(kz),
!     >   emaxkz(kz), elec(iez)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz)),
!     >  Sqrt(elec(iex)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 +
!     >       elec(iey)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 ),yesno,
!     >   zlev(kz)
      
      ENDIF
      end do
      
      IF ( ndebug >= 1 )  write(0,*) 'lgtgo,my_rank, nstep = ',lgtgo,my_rank,nstep
      
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
      
      DO nlit=1,nliter2

!#######################################################################
!  AT THIS POINT, IT HAS BEEN DECIDED THAT LIGHTNING HAPPENS . . .
!  NOW WE CAN CALL THE SLPM LIGHTNING MODULE . . . .     
!#######################################################################
!
   
!
!  zero t4 (tflg) and t2
!
      t1(:,:,:) = 0.0 ! for channel length
      t2(:,:,:) = 0.0 ! for number of channel points
      t4(:,:,:) = 0.0 ! for lightning-deposited charge


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

               t5(ix,jy,kz) = 1.e-6*0.5*eperao*(elec(iemag)%flt3d(ix,jy,kz)**2*dv +       &
                    ( -elec(iex)%flt3d(ix  ,jy  ,kz  )*phiw*dax          &
                      +elec(iex)%flt3d(ix+1,jy  ,kz  )*phie*dax          &
                      -elec(iey)%flt3d(ix  ,jy  ,kz  )*phis*day          &
                      +elec(iey)%flt3d(ix  ,jy+1,kz  )*phin*day          &
                      -elec(iez)%flt3d(ix  ,jy  ,kz  )*phid*dax          &
                      +elec(iez)%flt3d(ix  ,jy  ,kz+1)*phiu*daz          &
                    ) )

      end do
      end do
      end do

#ifdef MPI
      IF ( number_of_processes .gt. 1  ) THEN
            
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,t5)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,t5)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,t5)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,t5)
    
      ENDIF
#endif

!
!#######################################################################
! CALL THE SLPM     
!#######################################################################
         pchnl = 0
         nchnl = 0
         icheck = 0
         nnzt = nnz - 1
         
!         IF ( isa .eq. 1 ) iliter = nlit
      
!      allocate ( pot2(nnx,nny,nnz) )
!      allocate ( rhs(nnx,nny,nnz) )
!      allocate ( rho(nnx,nny,nnz) )


      iredo = 0

 17   CONTINUE   
      
      IF ( ny .gt. 2 ) THEN
       nnytmp = nny
       nyslmtmp = nyslm
!       write(iunit,*) 'mudlight: nny,nyslm = ',nny,nyslm
      ELSE
       nnytmp = 1
       nyslmtmp = 1
      ENDIF

#if defined (MPI) && defined (BOXMG)

      IF ( boxmginit == 0 ) THEN
        boxmginit = 1

      CALL mgsetuplg(nxslm,nyslm,nzslm,nor,      &
          dx,dy,dz,gxt,gyt,gzt,      &
          dxx,dyy,dzz,      &
          nnxs,nnys,nnzs,      &
          nbw,nbe,nbs,nbn,      &
          istretch,      &
          iunit, bcx, bcy )
        
         
         CALL PUTFLGT( SOlg, QFlg, QFlg, Qlg, Qlg, t0, elec(ipot)%flt3d,nxslm,nyslm,nzslm,nor,           &
                      nbw,nbe,nbs,nbn,      &
                      NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg,         &
                      iGslg, jGslg, kGslg, dslight, dslight, dz, 1          &
                     )
      

!        write(iunit,*) 'mudlight: nbw,nbe,nbs,nbn = ',nbw,nbe,nbs,nbn
!        write(iunit,*) 'mgsetuplg: my_rank,nbw,e,s,n: ',my_rank,nbw,nbe,nbs,nbn
!        write(iunit,*) 'mgsetuplg: my_rank,nnx,nnytmp: ',my_rank,nnx,nnytmp,nnz
!        write(iunit,*) 'mgsetuplg: my_rank,nlxdg: ',my_rank,nlxdg,nlydg,nlzdg
!        write(iunit,*) 'mgsetuplg: my_rank,nlxlg: ',my_rank,nlxlg,nlylg,nlzlg
!        write(iunit,*) 'mgsetuplg: nnz,nnzt,nlzlg: ',nnz,nnzt,nlzlg
      ENDIF
        
!        write(0,*) 'lightning flag ON in MPI'
!        lgtstp = 1
!        RETURN

      allocate ( rho(NLxlg+2,NLylg+2,NLzlg+2) )
      allocate ( work(2) )
      allocate ( Qmask(NLxlg+2,NLylg+2,NLzlg+2,3) )
      allocate ( imask(NLxlg+2,NLylg+2,NLzlg+2,1) )

       BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 1

      IF ( ndebug >= 1 ) write(0,*) 'call spark ,my_rank, nstep = ',my_rank,nstep

      call mudslpmmpi      &
         (time_real,time,nstep,iliter,lgtstp,pchnl,nchnl,ldumpsp      &
         ,nor,nx,ny,nz,NLxlg,NLylg,NLzlg,elec(ipot)%flt3d      &
         ,icldtop,zcldtop      &
         ,icgyn,icgpn,ebrkd      &
         ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
         ,dzmax,istretch      &
         ,gxt,gyt,gzt,dxx,dyy,dzz      &
         ,ixst,jyst,kzst,outname,nkzmn,icheck,      &
          nnx, nnytmp, nnz, nnzt,      &
          iixp,jjyq,kkzr,iiex,jjey,kkez,      &
          nbw,nbe,nbs,nbn,nbz2,               &
!     :   llwork,      &
          nxslm,nyslmtmp,nzslm,nb,nnx,nnytmp,nnzt,nnz,nxyz,nwslv,      &
          Qlg,Qlg,qmask(1,1,1,2),QFlg,QFlg,rho,rho,qmask(1,1,1,1),Qlg,      &
          numic,numcgn,numcgp,qmask(1,1,1,3),qmask(1,1,1,3),imask(1,1,1,1),      &
           iunit,bcx,bcy)

       deallocate ( rho )
       deallocate ( work )
       deallocate ( Qmask )

       BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = bmg_maxiter

#else
      
      allocate ( pot2(nnx,nnytmp,nnz) )
      allocate ( rhs(nnx,nnytmp,nnz) )
      allocate ( rho(nnx,nnytmp,nnz) )
      allocate ( work(llwork) )

      call mudslpm      &
         (time_real,time,nstep,iliter,lgtstp,pchnl,nchnl,ldumpsp      &
         ,nor,nx,ny,nz,elec(ipot)%flt3d,elec(iex)%flt3d      &
         ,elec(iey)%flt3d,elec(iez)%flt3d,elec(iemag)%flt3d      &
         ,icldtop,zcldtop      &
         ,icgyn,icgpn,ebrkd      &
         ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
         ,dzmax,istretch      &
         ,gxt,gyt,gzt,dxx,dyy,dzz      &
         ,ixst,jyst,kzst,outname,nkzmn,icheck,      &
          nnx, nnytmp, nnz, nnzt,      &
          iixp,jjyq,kkzr,iiex,jjey,kkez,      &
          nbw,nbe,nbs,nbn,nbz2,      &
          llwork,      &
          nxslm,nyslmtmp,nzslm,nb,nnx,nnytmp,nnzt,nnz,nxyz,nwslv,      &
          pot2,pot2,pot2(1,1,2),rhs,rhs,rho,rho,work,      &
          numic,numcgn,numcgp,      &
           iunit,bcx,bcy)


      
       deallocate ( pot2 )
       deallocate ( rhs )
       deallocate ( rho )
       deallocate ( work )

#endif
!
! check for flag to dump info when lightning goes awry (turned off for now)
!
!      IF ( .false. ) THEN
!      
!       IF ( lgtstp .eq. 2 .and. iredo .le. 4 ) THEN
!         iredo = iredo + 1
!         ldumpsp = .false.
!         lgtstp = 0
!         GOTO 17
!       ENDIF
!       
!       IF ( ldumpsp .and. isa .ne. 1 ) THEN
! 
!      write(rstime,timfmt) int(dtp*(nstep-1))
!         write(nmliterx,'(i4.4)') iliter+1000
!        fnelec = hdrcfs(1:lhdrcfs)//filehead(1:lfilehead)
!     : //'.nelecx'//'.'//nmliterx//'.'//rstime(1:ihttima)
!      
!      open(unit=97,file=fnelec,form='unformatted')
!     
!        write(97) 2  ! stand-alone mini-dump
!        CALL ELEC2SAV(nx,ny,nz,nor,na,ldumpsp,97,
!     :  an,
!     :  t4,t3,t0,
!     :  istag,jstag,kstag )
!      
!       
!       write(97) t0
!
!       chgpos = 0
!       chgneg = 0
!       DO kz=1,nz-1
!       DO jy=1,ny-1
!       DO ix=1,nx-1
!        IF ( t0(ix,jy,kz) .gt. 0 ) THEN
!          chgpos = chgpos + t0(ix,jy,kz)/gt(ix,jy,kz,imapz)
!        ELSE
!          chgneg = chgneg + t0(ix,jy,kz)/gt(ix,jy,kz,imapz)
!        END IF
!       END DO
!       END DO
!       END DO
!       
!      
!      write(iunit,'(a,3(2x,1pe12.5))') 't0 pos/neg/net charge (C):',
!     :      chgpos*dx*dy*dz,chgneg*dx*dy*dz,(chgpos+chgneg)*dx*dy*dz
!
!       
!       ENDIF
!       
!      ENDIF  ! ( .false. )
       
       IF ( lgtstp .gt. 0 ) GOTO 2999

!  Count number of ICs and CGs and CGPs and CGNs and CAs
!
      ilgt = 0
      iic  = 0
      ica  = 0
      icg  = 0
      icgp = 0
      icgn = 0
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

      IF ( isa .lt. 1 .and. iredis .eq. 1 ) THEN

       write(0,*) 'sll: defunct option; check code!'
       STOP
      
!
!  accumululate scnet and nscnet
!
      scnetpi = 0.0
      scnetni = 0.0
       chgrtp = 0.0
       chgrtn = 0.0
       scdrop = 0.0
       sckeepp = 0.0
       sckeepn = 0.0
       idrop = 0

#ifndef MPI
! C$DOACROSS LOCAL(kz,jy,ix)
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(kz,jy,ix)
#endif
      do kz = 1,nz-kd1
      do jy = 1,ny-jd1
      do ix = 1,nx-id1
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
        t3(ix,jy,kz) = 0.0
        t4(ix,jy,kz) = 0.0
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
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,i6,2x,1pe12.5)')      &
         'Points, net charge dropped outside cloud: ',      &
          idrop,scdrop
      write(iunit,'(a,2(1x,1pe12.5))')      &
         'Net kept pos/neg lightning charge = ',      &
          sckeepp,sckeepn
      write(iunit,'(a,1pe12.5)')      &
         'Cumulative net charge dropped outside cloud: ',      &
          scdropnet
      write(iunit,'(a,2(1x,1pe12.5))')      &
         'Ambient net pos/neg charge outside cloud: ',      &
          scnetpi, scnetni
      ENDIF
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

      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

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
        end if
        if ( t4(ix,jy,kz) .lt. -0.5 ) then
          scnetni = scnetni + t0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
          scnet = scnet + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnet = nscnet + 1
          voln = voln + volzone
        end if
      end do
      end do
      end do

#ifdef MPI
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

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 10, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
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
#endif

!
!      voln = nscnetn*volzone
!      volp = nscnetp*volzone
      
      IF ( my_rank == 0 ) THEN
      if ( ndebug .ge. 1 ) then
      write(iunit,'(a,2(1x,1pe12.5))') 'Ambient charge scnetni,scnetpi ',      &
           scnetni,scnetpi
      end if
!
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')       &
             'comp scnetp',scnetp,nscnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')       &
            'comp scnetn',scnetn,nscnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')       &
            'comp scnet',scnet,nscnet
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')       &
            'comp volp',volp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')       &
            'comp voln',voln
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')       &
            'COULOMB SCNETP',nstep,scnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')       &
            'COULOMB SCNETN',nstep,scnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')       &
            'COULOMB SCNET',nstep, scnet
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
!  If scnetn or scnetp is to small to accomodate the charge in
!  voln or volp, respectively, then make a volume sufficient
!  to accomodate the charge so breakdown will not occur.
!  after that go to 59999 and repeat this  procedure above.
!  only one attempt will be required.
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      irepeat = 0
      if ( irepeat .gt. 1 ) then
      print*,'failed lightning owing to failure to redistribute charge'
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

        ELSE ! ibal = 2
!        print*, 'sctarg,scnetn,scnetp = ',sctarg,scnetn,scnetp
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
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')      &
             'scnet=',scnet,nscnet
      if ( ndebug .ge. 2 ) write(iunit,'(a)') 'set scdelp--t3'

      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
      if ( t4(ix,jy,kz) .gt. 0.5 ) then
      t3(ix,jy,kz) = fracp*t3(ix,jy,kz) - multp*scnet
      scnetnp = scnetnp + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetnp = nscnetnp + 1
      scnetn = scnetn + t3(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetn = nscnetn + 1
      end if
      if ( t4(ix,jy,kz) .lt. -0.5 ) then
      t3(ix,jy,kz) = fracn*t3(ix,jy,kz) - multn*scnet
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


      IF ( my_rank == 0 ) THEN

      if ( nscnetn .gt. 0 ) then
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',      &
         scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',      &
         scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET',      &
         scnetn, nscnetn
     
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
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       IF (t3(ix,jy,kz) .lt. 0) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
        qpzn = qpzn +       &
                z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        qnetn = qnetn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF (t3(ix,jy,kz) .gt. 0) THEN
        qpxp = qpxp + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyp = qpyp + gyt(jy,1)*t3(ix,jy,kz)*dv
        qpzp = qpzp +       &
                z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
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
       qpxcn = -qpxn*0.001/(-qnetn)
       qpycn = -qpyn*0.001/(-qnetn)
       qpzcn = -qpzn*0.001/(-qnetn)
       qpxcp = qpxp*0.001/(qnetp)
       qpycp = qpyp*0.001/(qnetp)
       qpzcp = qpzp*0.001/(qnetp)
        qpxm = (-qpxn + qpxp)*0.001/(qnetp-qnetn)
        qpym = (-qpyn + qpyp)*0.001/(qnetp-qnetn)
        qpzm = (-qpzn + qpzp)*0.001/(qnetp-qnetn)
        ENDIF


      IF ( my_rank == 0 ) THEN
      
      write(iunit,'(a,2(1x,i7))') 'npntn, npntp ',npntn,npntp
      
      write(iunit,'(a,2(1x,1pe12.5))')       &
           'qnetn, qnetp',qnetn*(1.0),qnetp*(1.0)

      IF ( isa .lt. 1 ) THEN

!      write(24,'(a,2(1x,1pe12.5))') 
!     :    'qnetn, qnetp',qnetn*(volzone),qnetp*(volzone)
      
      ENDIF
      
        write(iunit,'(a,4(1x,1pe12.5))') 'dipole moment (px,py,pz) ',      &
              qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )

      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,1pe12.5))') 'dipole moment (px,py,pz) ',
!     :       qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',      &
           qpxm,qpym,qpzm   
      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',
!     :    qpxm,qpym,qpzm   
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'pos. drho center (x,y,z) ',      &
           qpxcp,qpycp,qpzcp   
        write(iunit,'(a,3(1x,f12.5))') 'neg. drho center (x,y,z) ',      &
           qpxcn,qpycn,qpzcn   

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
      
      IF ( (ibal .eq. 5)  .or.       &
           ( ibal .eq. 2 .and. cgfr .lt. 1.0 ) ) THEN
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
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')      &
             'scnet=',scnet,nscnet
      if ( ndebug .ge. 2 ) write(iunit,'(a)') 'set scdelp--t3'

      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

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
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',      &
         scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',      &
         scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET ',      &
         scnetn, nscnetn
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

      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
       do jy = 1,jye
        do ix = 1,ixe

       IF ((t4(ix,jy,kz)) .ne. 0.0) THEN

         dv = dxx(ix)*dyy(jy)*dzz(kz)

       IF ( (icgpn .eq. -1 .or. ibal .eq. 5 .or. ibal .eq. 2)      &
             .and. t3(ix,jy,kz) .gt. 0.0 ) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*t3(ix,jy,kz)

        qpzn = qpzn +       &
                z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

        qnetn = qnetn + t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF ( (icgpn .eq. 1 .or. ibal .eq. 5 .or. ibal .eq. 2)      &
             .and. t3(ix,jy,kz) .lt. 0.0 ) THEN
        qpxn = qpxn + gxt(ix,1)*t3(ix,jy,kz)*dv
        qpyn = qpyn + gyt(jy,1)*t3(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*t3(ix,jy,kz)

        qpzn = qpzn +       &
                z1d2(kz,1)*t3(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

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
       write(iunit,'(a,1pe12.5)')       &
          'Charge center of ',qnetn ! *volzone
       write(iunit,'(a,3(2x,f11.4))')       &
          'at x,y,z (km)', qpxcn,qpycn,qpzcn
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
      if ( ndebug .ge. 10 ) print*,'sum total particle area'
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
       STOP
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
       write(iunit,'(i3,2(2x,1pe12.5))' ) ia,      &
              chgrt(1,ia-lscb+1),chgrt(2,ia-lscb+1)
      END DO ! ia
      
      write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',chgrtp,chgrtn
      
      end if
!
!  redistribute based on fraction
!
      if ( iredis .eq. 2 ) then
      do ia = lscb,lsceq
      do kz = 1,nz-kd1
      do jy = 1,ny-jd1
      do ix = 1,nx-id1
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

      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

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
       write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',      &
         chgrtp,chgrtn
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


      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      jye = ny
      IF ( myprocj == nprocj ) jye = ny-1

       DO kz = 1,nz-kd1
       DO jy = 1,jye
       DO ix = 1,ixe
       
        IF (Abs(t2(ix,jy,kz)) .lt. 0.5 .and.       &
             Abs(t4(ix,jy,kz)) .gt. 0.5 ) THEN
            t4(ix,jy,kz) = Sign( 0.6,t4(ix,jy,kz) )
        ELSE
            t4(ix,jy,kz) = t2(ix,jy,kz)
        END IF
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
      end if  !   ( ilight .eq. 2 ) 
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
!

!
!
!  END OF LIGHT
!
!
!      subroutine mudslpm1
!     >  (nstep,iliter,lgtstp,pchnl,nchnl,ldumpsp
!     >  ,nor,nx,ny,nz,elec,neelec
!     >  ,ngt,icldtop,zcldtop,cgthres
!     >  ,icgyn,ibrkd,icgpn,ebrkd,fefac,ibal
!     >  ,ieint,eint,fecrit,feint,elgtfdel,iseed,intfg,mdel
!     >  ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
!     >  ,gt,imapz,mzdist,dzmax,istretch
!     >  ,ixst,jyst,kzst,outname,nkzmn,icheck,
!     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
!     :   nbw,nbe,nbs,nbn,nbz2,
!     :   llwork,
!     :  nxslm,nyslm,nzslm,nb,nnx,nny,nnzt,nnz,nxyz,nwslv
!     :   ,nliter,numic,numcgn,numcgp,
!     :    potfair,rhofair,cgfr)
!      
!      integer nnx,nny,nnz
!
!      real,allocatable :: pot2(:,:,:) 
!      real,allocatable :: rhs(:,:,:)
!      real,allocatable :: rho(:,:,:)
!      real,allocatable :: work(:)
!
!      allocate ( pot2(nnx,nny,nnz) )
!      allocate ( rhs(nnx,nny,nnz) )
!      allocate ( rho(nnx,nny,nnz) )
!      allocate ( work(llwork) )
!      
!            call mudslpm
!     >  (nstep,iliter,lgtstp,pchnl,nchnl,ldumpsp
!     >  ,nor,nx,ny,nz,elec,neelec
!     >  ,ngt,icldtop,zcldtop,cgthres
!     >  ,icgyn,ibrkd,icgpn,ebrkd,fefac,ibal
!     >  ,ieint,eint,fecrit,feint,elgtfdel,iseed,intfg,mdel
!     >  ,dtp,dx,dy,dz,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
!     >  ,gt,imapz,mzdist,dzmax,istretch
!     >  ,ixst,jyst,kzst,outname,nkzmn,icheck,
!     :   nnx,nny,nnz,nnzt,
!     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
!     :   nbw,nbe,nbs,nbn,nbz2,
!     :   llwork,
!     :  nxslm,nyslm,nzslm,nb,nnx,nny,nnzt,nnz,nxyz,nwslv,
!     :   pot2,pot2,pot2(1,1,2),rhs,rhs,rho,rho,work
!     :   ,nliter,numic,numcgn,numcgp,
!     :   potfair,rhofair,cgfr)
!      
!       deallocate ( pot2 )
!       deallocate ( rhs )
!       deallocate ( rho )
!       deallocate ( work )
!      
!      return
!      end
!
!
