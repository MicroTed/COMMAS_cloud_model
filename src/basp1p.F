c
c c23428789012342878901234287890123428789012342878901234287890123428789012
c 11/2005 Update for double moment (set ipconc based on na)
c
c
c 5/18/2000 Add capability to read flash file (flash composites) for contouring
c           (no change to input file; uses name of 'ht1' file to create filename
c            of the flash composite file at the same dump time.)
c
c 8/18/99 new version allows use of compact gather array dump files
c         instead of full history dump.
c         use ihdmpt = 3 to read 
c            filnm1: q and sc (post-microphysics) (puts into a array)
c            filnm2: read the history.base file
c         use ihdmpt = 4 to read filnm1 and filnm2: 
c              post- (a) and pre- (ab) microphysics
c 11/99 adapted to use arbitrary oblique vertical cutting plane.  Spacing along
c       the slice is assumed to be dx.  Cues for cutting through some maximum value still
c       work but go in through the second value in the ib ie pair.
c
      subroutine cplotba(gd,                 
     :   u,  uinit ,                     ! U,  UINIT
     :   v,  vinit ,                     ! V,  VINIT
     :   w,  winit ,                     ! W,  WINIT
     :   expi, piinit,                     ! PI, PIINIT
     :   km, kminit,                     ! KM, KINIT
     :   s,  sinit0,                     ! S,  SINIT
     :   precip,                         ! PRECIP
     :   gx,                             ! XCNTR, XEDGE, DXC, DXE
     :   gy,                             ! YCNTR, YEDGE, DYC, DYE
     :   gz,                             ! ZCNTR, ZEDGE, DZC, DZE
     :   ugrid, vgrid,                   ! GRID MOTION
     :   nx, ny, nz, na,                 ! NX,NY,NZ,NS
     :   st,                             ! scalar array
     :   dn, pn, tt0, rho,
     :   dbz, wz, elec, xtra,
     :   sinit,
     :   iunit,itime,member,dx,dy,dz,dtp,iframe,precip_old)
c                    (dx,dy,dz,nx,ny,nz,na,nor,nv)

      USE PLOT_MODULE
      USE GRID_MODULE
      USE PARAM_MODULE
      USE MICRO_MODULE
      USE INDEX_MODULE
      use dualpara,  only: setgrplhl, model_dsd
      use rsa_table, only: set_nbins

      implicit none

      TYPE(GRID) :: gd

      TYPE(VARIABLE) :: u, uinit
      TYPE(VARIABLE) :: v, vinit
      TYPE(VARIABLE) :: w, winit
      TYPE(VARIABLE) :: expi, piinit
      TYPE(VARIABLE) :: km, kminit
      TYPE(VARIABLE) :: s(na), sinit0(2)     ! We are going to try and create arrays of the scalar variables needed here
      TYPE(VARIABLE) :: precip(nprecip+neelec2d)
      TYPE(VARIABLE) :: gx(4), gy(4), gz(4)
      TYPE(VARIABLE) :: dbz, wz

      TYPE(VARIABLE) :: xtra(nxtra)
      TYPE(VARIABLE) :: elec(neelec)
      
      real           :: sinit(nz,na)

      real :: dn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! air density
      real :: pn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! pressure
      real :: tt0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! temperature
      
      real :: precip_old(-ng+1:nx+ng,-ng+1:ny+ng,nprecip)
      
      real ugrid, vgrid     


!      real    :: gxt(nx,4), gyt(ny,4), gzt(nz,4)
      real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
      
      integer iunit
      integer member
      
      real maxval, newval
      real VMX,VMN,dmx,VL,VR,VB1,VT,UL,UR,UB,UT,LL,vrl
      integer inx,iny
      logical lskyb,lskyt
      integer iskyb,iskyt
      
      integer icgsymbol,icgxy ! whether to plot symbols for +/-CG initiation points
      data icgsymbol/0/
      integer igmname   ! whether to set the name of the ncgm file
      data igmname/0/
      
      integer icontent
      
      real rdm  ! dummy value
      integer idm ! dummy value
      
      real cross1(2,2),cross2(2,2),dx1,dy1
      real triang(4,2)
      real diamond(5,2)
      
      real scale
      
      integer ieps, inameps, epsflag,ifirsteps
      character*120 outnameps
      save inameps

      
      integer numtime   ! number of time dumps to cycle through.
      integer ntime ! counter to keep track of time cycle
      real dtdump   ! time interval between dumps in cycle
      integer :: ifirst = 0
      save ifirst
      integer indx,iargc,iargc0
      character*120 inputfile
      
      integer iwdt   ! flag whether wdt files have been opened already
      data iwdt/0/
      integer iwdtas,iwdtbin  ! units for ascii and binary output
      character*40 fwdt  ! base file name for wdt output
      real umv,vmv,xoffset,yoffset
      real umv1,vmv1,dz1
      integer ioffset
      data ioffset /0/
      data xoffset /0.0/
      data yoffset /0.0/
      
      integer istat1
      integer iexp,iclr,nzset,n2set
      integer ifll
      real hgt,hgt0,ycn,tmp,tmp2,hgtset,hgtlow
      real ZCD(2),XCD(2)
      logical fshortht
      integer nliter
      integer istrz  ! flag =0 for unstretched vertical, =1 for stretched
      
      character*100 line
      
c      integer ig13,ig23
      integer iscrn  ! flag for existence of screening charge array
      integer ion    ! flag for existence of ion number densities
      integer mzdist,kdbcdn
      parameter (mzdist = 4 )
      integer nxb2,nyb2,nzbt2,nzb2,nb2,n
      integer nxb2s(2),nyb2s(2),nzbt2s(2),nzb2s(2),nb2s(2)
      real chgneg,chgpos ! total net negative/positive charge in storm (not including corona)
      real cnaneg,cnapos ! total net negative/positive corona charge
c
c Mudpack variables
c
      integer imud,imudstop
      integer nzmax,nzmax1,nz1d,nztop
      parameter (nzmax1 = 400)
      real z1d2(nzmax1,4)
      real dlz,dzmax
      integer iixps(2),jjyqs(2),kkzrs(2)
      integer iiexs(2),jjeys(2),kkezs(2),llworks(2)
      integer nnxs(2),nnys(2),nnzs(2),nzs

      integer nbz2,nbw,nbe,nbs,nbn
      
      save imud,nbz2,nbw,nbe,nbs,nbn,nztop
c      save z1d2
      save iixps,jjyqs,kkzrs,iiexs,jjeys,kkezs,llworks,nnxs,nnys,nnzs
      
      real eztop

      real templt(nzmax1),temkz(nzmax1)
      real  felv(nzmax1), felf(nzmax1), fels(nzmax1)
      real temg,temcg

      
      common/bndcz/eztop
      
      real rhomx,rhomn
      integer ngt,nmax,itime,iframe
      character*6 stime,stimeplot
      character*7 smin, shr
      character*4 smin2
!      real gr
      
!      integer imsz  ! =1 for msz flashes, =0 for SLM flashes
      character*120 xfname(5),bname,segname
      integer iusef,ixtraf

      logical fenele,feht1,femic,femica,febase,feflst,fesll

!      character*80 inname,outname,inname2,flstname,cdum
      character*80 inname3,inname4,basename,insll
      character*10 titlelabel  ! name of the field
      real xib,xib0,xie,yjb,yje,zkb,zke,dd,dd0,phi,pi,dd1,dd2
      integer ipltyp,nl
      
      real theta1
      
      integer jl,ju,ji,ix,ix2,kz,jy,no,np,nnc,ia,na !,nh,
      integer id1,jd1,kd1,istag,jstag,kstag
      integer itype,irm,im,imore,mu,mv,mw,nsmth
      integer itck1,itck2,nplot,ntimep(100)
      real cstag
      
c      include 'sam.plot.dims'
      
      real tem1(nx-1,ny-1,3)  ! Ezsfc, reflectivity, lgt
c      real tem2(nx,ny)
c      real tem3(nx,ny)
      
      real hdist(nx+ny),zk(nz)
      real zdb(nz),scalez(nz)
      integer nhdist
      integer nx,ny,nz
c      parameter(nx=nx,ny=ny,nz=nz)
      integer nxelec,nyelec,nzelec
      parameter (np=1,nnc=1,ngt=10)
      parameter (id1=1,jd1=1,kd1=1,istag=1,jstag=1,kstag=1)
      integer nxstag,nystag,nzstag

c      parameter (neelec=6) ! ,nxelec=nx,nyelec=ny,nzelec=nz)
c      integer iex,iey,iez,iemag,ipot
c      parameter (iex=1,iey=2,iez=3,iemag=4,ipot=5)
      
      integer neletyp  !  indicator of what is in the nelec file
      data neletyp /1/

       integer
     :        iqi,iqc,iqr,iqs,iqg,iqh,
     :        isci,isccw,iscr,iscs,iscg,isch
c       parameter( iqi=4,iqc=5,iqr=6,iqs=7,iqg=8,iqh=9,
c     :     isci=10,isccw=11,iscr=12,iscs=13,iscg=14,isch=15 )

c
c  arrays used in FISHPAK for solver for electric potential (not used anymore)
c
      integer    imapx
      parameter (imapx=07)
      integer    imapy
      parameter (imapy=06)
      integer    imapz
      parameter (imapz=05)
      integer    ig13
      parameter (ig13=09)
      integer    ig23
      parameter (ig23=09)
      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)

      integer nxb,nyb,nb
      parameter(nb = 4)
c      parameter(nxb = nxelec + 2*nb + 1, nyb = nyelec + 2*nb + 1)
c      integer nwslv
c      parameter (nwslv=nxb+nyb+5*nzelec+30+nxb*2
c     >          +7*nxb/2+7*nyb/2)
c      real wslv(nwslv)
c      real bdxs(nyb,nzelec),bdxf(nyb,nzelec)
c      real bdys(nxb,nzelec),bdyf(nxb,nzelec)
c      real bdzs(nxb,nyb),bdzf(nxb,nyb)

c      real pot(nxb,nyb,nzelec) 
c      real potw(-nb:nxm+nb,-nb:nym+nb,nzelec)
c      equivalence( pot(1,1,1), potw(-nb,-nb,1) )
c      real work2(nyb,nxb),work2a(nyb,nxb)
c      real work3(nxb,nyb),work3a(nxb,nyb)
     
c
c
c  electricity
c
c      integer  neelec

c      integer  iex,iey,iez,iemag,ipot,iezb
c      parameter (iex=1,iey=2,iez=3,iemag=4,ipot=5,iezb=6)
      
      logical istretch
      data istretch/.true./
      
c arrays for phi=0 below ground
      integer nzb
c      parameter (nzb=nzelec+1)
c      real bdxs2(nyb,nzb),bdxf2(nyb,nzb)
c      real bdys2(nxb,nzb),bdyf2(nxb,nzb)
c      real bdzs2(nxb,nyb),bdzf2(nxb,nyb)

c      real pot2(nxb,nyb,nzb) 
c      real pot2w(-nb:nxm+nb,-nb:nym+nb,nzb)
c      equivalence( pot2(1,1,1),pot2w(-nb,-nb,1) )
c      equivalence(pot(1,1,1), pot2(1,1,2))
c      equivalence(bdxs(1,1), bdxs2(1,2))
c      equivalence(bdxf(1,1), bdxf2(1,2))
c      equivalence(bdys(1,1), bdys2(1,2))
c      equivalence(bdyf(1,1), bdyf2(1,2))

c
c fair weather field 
c
      integer nzfair
      parameter (nzfair=200)
      real ezfair(nzfair),potfair(nzfair),rhofair(nzfair)
      real ezfairw(nzfair)
      save potfair,ezfair,rhofair,ezfairw
      
      real ezfairo
      parameter(ezfairo = -80.00)
      real efa1, efa2, efa3, efb1, efb2, efb3
      parameter(efa1 = 4.5e-3, efa2 = 3.8e-4, efa3 = 1.0e-4,
     :          efb1 = 0.50,   efb2 = 0.65,   efb3 = 0.10    )
     
      real ezfairs        ! fair field at surface
      parameter(  ezfairs = -100.0 )
      real dezdz
      real zgrd
c vars for magical corona layer
      
c      real sccna(nx,ny)

c
c  electrical permitivity of air C / (N m**2) -  check the units
c
      real eperao,dslay
      parameter (eperao  = 8.8592e-12 )
      
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      
c
c  constants for boundary conditions for poisson solver
c
      integer ibdcnd
      integer jbdcnd
      integer kbdcnd
c other stuff for poisson solver
      integer ierror
      real xs, xf, ys, yf, zs, zf  ! domain limits in meters
      real pertrb,elmbda
      
c      include 'sam.index.f'
c
c      integer lt,le
c      parameter (lt=1,le=2)
c
c   mixing ratio
c
c      integer lv,lc,lr,li,lir,ls
c      integer lgl,lgm,lgh,lf,lh
c      integer lip,lhl,lhab

c      integer lqb,lqe

c      parameter (lv=3,lc=4,lr=5,li=6,lir=7,ls=8)
c      parameter (lgl=9,lgm=10,lgh=11,lf=12,lh=13)
c      parameter (lip=14,lhl=15,lhab=15)
c
c  number concentration (was 16 to 27)
c
c      integer lnc,lnr,lni,lnir,lns
c      integer lngl,lngm,lngh,lnf
c      integer lnh,lnip,lnhl
      integer lnshi
c      parameter (lnshi = 64)
c      parameter (lnc=0+lnshi,lnr=1+lnshi,lni=2+lnshi,lnir=3+lnshi)
c      parameter (lns=4+lnshi,lngl=5+lnshi,lngm=6+lnshi,lngh=7+lnshi)
c      parameter (lnf=8+lnshi,lnh=9+lnshi,lnip=10+lnshi,lnhl=11+lnshi)
c
c  space charge (28 to 39)
c
c      integer lscw,lscr,lsci,lscir,lscs
c      integer lscgl,lscgm,lscgh,lscf,lsch
c      integer lscip,lschl

c      integer lscb,lsce,lsceq,lscwi

c      parameter (lscw=28,lscr=29,lsci=30,lscir=31,lscs=32)
c      parameter (lscgl=33,lscgm=34,lscgh=35,lscf=36,lsch=37)
c      parameter (lscip=38,lschl=39)
c
c  temperature (was 40 to 51)
c
c      integer ltpc,ltpr,ltpi,ltpir,ltps
c      integer ltpgl,ltpgm,ltpgh,ltpf,ltph
c      integer ltpip,ltphl
c      integer ltshi
c      parameter (ltshi = 52)
c      parameter (ltpc=0+ltshi,ltpr=1+ltshi,ltpi=2+ltshi)
c      parameter (ltpir=3+ltshi,ltps=4+ltshi,ltpgl=5+ltshi)
c      parameter (ltpgm=6+ltshi,ltpgh=7+ltshi,ltpf=8+ltshi)
c      parameter (ltph=9+ltshi,ltpip=10+ltshi,ltphl=11+ltshi)
c
c  mean cloud content exposure (was 52 to 63)
c
c      integer lamcw,lamcr,lamci,lamcir,lamcs
c      integer lamcgl,lamcgm,lamcgh,lamcf,lamch
c      integer lamcip,lamchl
c      integer lashi
c      parameter (lashi = 16)
c      parameter (lamcw=0+lashi,lamcr=1+lashi,lamci=2+lashi)
c      parameter (lamcir=3+lashi,lamcs=4+lashi,lamcgl=5+lashi)
c      parameter (lamcgm=6+lashi,lamcgh=7+lashi,lamcf=8+lashi)
c      parameter (lamch=9+lashi,lamcip=10+lashi,lamchl=11+lashi)
c
c  tau  (was 64 to 75)
c
c      integer lrtc,lrtr,lrti,lrtir,lrts
c      integer lrtgl,lrtgm,lrtgh,lrtf,lrth
c      integer lrtip,lrthl
c      integer lrshi
c      parameter (lrshi = 40)
c      parameter (lrtc=0+lrshi,lrtr=1+lrshi,lrti=2+lrshi)
c      parameter (lrtir=3+lrshi,lrts=4+lrshi,lrtgl=5+lrshi)
c      parameter (lrtgm=6+lrshi,lrtgh=7+lrshi,lrtf=8+lrshi)
c      parameter (lrth=9+lrshi,lrtip=10+lrshi,lrthl=11+lrshi)
      
c

c      real flashm(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) !  flash initiation points
c      real flashp(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! number of 'positive' segments (negative leaders)
c      real flashn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! number of 'negative' segments (positive leaders)
c      real drhom(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      integer nht,nhsg,nhg,nhsh,nhnh,nhl,nhgh,lhcm,nhcm,nxyz2
      parameter (nht=30,nhsg=2,nhg=3,nhsh=4,nhnh=5,nhl=6,nhgh=7)
      parameter (lhcm=7,nhcm=14)
      parameter (mu=1,mv=2,mw=3)
      parameter (no=0,nmax=1000)
      parameter (cstag=0.5)
      parameter (nxyz2=nmax+2)
!      integer ifield(nx,ny,nz)
      integer ifile
      integer :: ibase, ibase_mic, ibase_mic_add
      character*6 labdx,labdy,labdz
      character*6 lprn1,lprn2,lprn3
      integer ntit
      parameter(ntit=1200)
      character*40 title(ntit)
      character*80 fildum,filnm1,filnm2,filnm3,filnm4,filnm5
      character*60 label
      integer ilabel
      character*40 header, tlab, xlab, ylab
      character*9 lmaxvl,lminvl,lmin,lcnt,ltim,lloc
      character*28 label1,label2,label3,label4,label5,label6
      character*50 blab1,blab2,filnmtraj

      character*71 xtitle
      integer titlen,xtitlen
      character*1 labcrv
      integer ncbar
      parameter (ncbar=129)
      character*6 labcbar(129)
      integer ifill(129)
      integer nh
      parameter (nh = 1 )
      integer lszhdn(nz,nh)
      integer irast2d(nmax,nmax)
      real hd(nh,nz),hdn(nz,nh)
      real spv(2),dumhd(nh),dumhdn(nh)
c      real gt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ngt)
c      real ht(-ng+1:nx+ng,-ng+1:ny+ng,nht)
      real st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na+1)
      integer nv
c      real v(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nv)
c      real p(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real a2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2,na2)
c      real v2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2,nv2)
c      real p2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2)
c      real d2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2)
c      real ab(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na+1)
c      real vb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nv)
c      real pb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real db(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    :: pb(nz)
      real    :: db(nz)
 
      real thdnm(nh),zhdnm(nz,nh) ! ,hdnx(nxm,nym,nz,nh)
      real, allocatable, save, dimension(:,:,:) :: rh, the, temc
      
!      real rh(nx,ny,nz),the(nx,ny,nz)

c      real dbz(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real temc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      

c      real qc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real qi(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real rho(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real scrnchg(-ng:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real algt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real abz(nz,na),vbz(nz,3),zb(nz),pbz(nz),eye(3)
      real slab3(nxyz2,nxyz2),piz(nz),dnz(nz)
      real qhd(nx,na),qhcn(nx,na),qhdn(nx,na),sdum(nx,na)
      real qhnum(nx,na)
      real time,time0,xmx,xmn,tiso
      real plt2da(nmax,nmax),cntz(nz)
c      integer LZDT
c      parameter(LZDT = nmax*nmax)
c      real ZREG(LZDT)
c      real zdat(nmax,nmax)
c      equivalence (zreg(1),zdat(1,1))
      real plt2dp(nmax,nmax),plt2dn(nmax,nmax)
      real tmp2d(nmax,nmax)
      real, allocatable, save, dimension(:,:,:) :: plt3dp, plt3dn
!      real plt3dp(nx,ny,nz),plt3dn(nx,ny,nz)
      real xvec(nmax,nmax),yvec(nmax,nmax)
      real zvec(nmax,nmax)
      real qq(np),chi(1000),clo(1000) ! ,cin(1000)
      real xplt(5000),yplt(5000),zplt(5000)
      integer iflag

c      real t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real t1(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real t2(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real t3(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

c      real tema(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real tke(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real tkeb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real thetap(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real theta(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real thetab(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real qvb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

c      real scionp(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real scionn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      
c     equivalence(t0,dum)
C
C Declare the required real and integer workspaces.
C
        integer LRWK,LIWK
        Parameter(LRWK=20000000,LIWK=20000000)
        real RWRK(LRWK)
        integer IWRK(LIWK)
C
C Declare an array to hold an area map.
C
        integer LIAMA
        parameter (LIAMA=20000000)
        integer IAMA(LIAMA)
C
C Declare the arrays needed by ARSCAM for x/y coordinates.
C
        integer MCS
        parameter(MCS=200000000)
        real XCRA(MCS),YCRA(MCS)
C
C Declare the arrays needed by ARSCAM for area and group identifiers.
C
        integer MAI
        parameter (MAI = 40)
        integer IARA(MAI),IGRA(MAI)
C
C Declare an array to hold the GKS "aspect source flags".
C
        integer IASF(13)
C
C Declare the routine which will draw contour lines, avoiding labels.
C
        EXTERNAL DRAWCL
C
C Declare the routine which does the shading.
C
        EXTERNAL SHDHR
        EXTERNAL SHDHL
        EXTERNAL SHDHH
        EXTERNAL SHDHV
        EXTERNAL SHDHX
        EXTERNAL SHDD2
        EXTERNAL SHDD4
        EXTERNAL SHDD6
        EXTERNAL SHDD8
        EXTERNAL SHDFIL
        EXTERNAL SHDGRAY1
        EXTERNAL SHDBLUE1
        EXTERNAL SHDRED1
        EXTERNAL SHDBLUE2
        EXTERNAL SHDRED2
        EXTERNAL SHDGREEN1
        EXTERNAL SHDCLR
c
        EXTERNAL SHDH1
        EXTERNAL SHDH2
        EXTERNAL SHDDOT
        EXTERNAL COLRAM
c        external SETDBZCLR
        external COLDBZ
        external REDPOS
        external BLUNEG
        external BLURED
        external CSTCLR
C
C Initialize the values in the aspect-source-flag array.
C
        DATA IASF / 13*1 /
C
      integer i,j,k
      integer icomp
      real dzadd,dxi,dyi,dzi,flov,hiv
      integer ichanflag

c      real dt0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! temperature
c      real chan(nx,ny,nz)
c      real dat(nx,ny,nz)
c      real drho(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real flash(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      equivalence( algt(-ng,-ng,-ng), flash(-ng,-ng,-ng) )

c      real chan2(nx,ny,nz)
c      real drho2(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
c      real flash2(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      integer nxt,nyt,nzt
      integer nend
      parameter(nend=2)
      integer nxyz3
      parameter ( nxyz3 = 10000/nend )
      integer numchan(2)
      real signchan(2)
c      integer ibonds(nxyz3,3,2,nend) ! list of connected points
c      real chanrho(nxyz3,nend)  !  channel charge densities
      real dum1,dum2
      integer jend
      real ds,eps
      real chanq(nend),frac(nend)
      integer icg,icgn,icgp
      character*10 flshtitle
      character*1 flshnum
      integer idum,i1,j1,k1,i2,j2,k2,iend
      real minq,volinv,x,y,z
      
c      real   rwtot2d(nx,ny),rwrt2d(nx,ny)
c     >      ,swtot2d(nx,ny),swrt2d(nx,ny),gltot2d(nx,ny)
c     >  ,glrt2d(nx,ny)
c     >      ,gmtot2d(nx,ny),gmrt2d(nx,ny),ghtot2d(nx,ny)
c     >   ,ghrt2d(nx,ny)
c     >      ,fwtot2d(nx,ny),fwrt2d(nx,ny),hwtot2d(nx,ny)
c     >   ,hwrt2d(nx,ny)
c     >      ,prt2d(nx,ny),ptot2d(nx,ny),prti2d(nx,ny)
c     >   ,ptoti2d(nx,ny)
c     >      ,prtl2d(nx,ny),ptotl2d(nx,ny)
c      
      real dtp,dx,dy,dz,dzsave,dzin,dxi2,dyi2,dzi2
      real dt,dt0
      integer ispv, ismthx, ismthy, ismthz
      real facyv,facxv
      integer jypu,jymu,ixpv,ixmv,ioverlay,icnt,ioverpnt,ioverpnt2
      integer idone,iforceol,idx,idx0,ifread,ivint,ifvint
      real fhiv,mcl,cll
      integer ido,ipware,ishade,icrfil
      integer ib,ie,jb,je,kb,kstrzb,kb1,ke,icrtyp,idolab
      integer irev  ! flag to reverse xlabs
      integer ilabor,ilablab,ishft,ishadcnt0
      integer itshft
      real df
      integer isoflg,imx,imn,jmx,jmn,kmx,kmn,imax
      real dxp,dyp,dzp,d1,d1p,d1pobl,fac1,d02,d2p,fac2
      real dzt
      integer ixdiv,jydiv, kzdiv
      integer i3,itdiv
      real axb,ayb,axe,aye,axep,awrk
      integer iform,isign,ilabcrv,ilabcrvtyp,it,ii,ICLV
      real clev,disfac,timew
      integer NOCL
      integer lunr,lunw
      integer lunr1,lunw1
      real xfacb0,xdivb0
      real yfac,zfac,zdiv
      integer nsetc,nsetv,istat
      real lengthv
      real qwv,pres,thetae,tmr,qsat
      integer kzlow,kzhigh
      real rclon,rclop,rchin,rchip,rcinn,rcinp,xmult
      real rlo,rhi,rin,sf
      real sdx,sdy,xfac,xdiv,ydiv,sdz,xmax,xmin,uava,vava
      real xmax2,xmax3,xmax4
      integer nclr,nxp,nyp,nzp,n1,n2
      real rx,ry,rxy,xbarr,xbarl,ybart,ybarb
      integer n3,m,mjrx,mnrx,mjry,mnry,ndot,numtraj,num,nt
      real width,rclo
      integer nclv,mxyz,mxyz2,mxyz2s
      integer posdef
      
      integer ishdclr
c
c  common 'plot'
c
c      common /ppppp/ ab,db,pb,vb,a,d,p,v,elec,algt,
c     > umove,vmove,dx,dy,dz,dxi2,dyi2,dzi2,temc,flashp,flashn,
c     > drhom,flashm,rho,dbz,the,rh

      common /oblique/xib,xie,yjb,yje,zkb,zke,ipltyp,phi,dd,nl,d1pobl
      common /bastuff/femic,inname3
      common /clrstuff/ ishdclr
      common /linestuff/ cll

        integer ifll0
        common/brparam/ ifll0
        data ifll0/0/

c      common /indices/lt,le,lv,lc,lr,li,lir,ls,lgl,lgm,lgh,
c     :    lf,lh,lip,lhl,lhab,lqb,lqe,
c     :    lscw,lscr,lsci,lscir,lscs,
c     :    lscgl,lscgm,lscgh,lscf,lsch,
c     :    lscip,lschl,lscb,lsce

c
c  conrec common block 'conre4'
c
      integer IOFFP, SPVAL
      integer isizel,  isizem,  isizep,   nrep,
     >                ncrt,   ilab,   nulbll,  ioffd,
     >                ext,    ioffm,  isolid,  nla, nlm
      real            xlt,    ybt,     side
      COMMON /CONRE1/ IOFFP      ,SPVAL
      common /conre4/ isizel,  isizem,  isizep,   nrep,
     >                ncrt,   ilab,   nulbll,  ioffd,
     >                ext,    ioffm,  isolid,  nla,
     >                nlm,    xlt,    ybt,     side
c
c  color for conrec common block 'recint'
c
      integer irecmj,irecmn,irectx
      common /recint/ irecmj,irecmn,irectx
c
      real chsval, chsang
      common /hsaw/  chsval, chsang

      common/strzb/kstrzb,hgtlow
      common/strz/z1d2,dlz,nz1d

      integer maxt,itcount
      parameter(maxt=600)
      real timearr(maxt),times
      real heights(nz)
      
      character*200 dumstr
      
      integer maxc
      parameter(maxc=100)
      real charr(nz,maxt,maxc)
      real chtotarr(nz,maxc)
      
      real cinit(nz,maxt),cinittot(nz) ! whole domain inits
      real cvinit(nz,maxt),cvinittot(nz) ! convective region inits
      real aninit(nz,maxt),aninittot(nz) ! anvil region inits
      
      real chans(nz,maxt),chanstot(nz) ! whole domain
c      real chansn(nz,maxt),chansp(nz,maxt) ! whole domain
      real chanstotn(nz),chanstotp(nz) ! whole domain

c      real cvchans(nz,maxt),cvchanstot(nz)  ! convective region
c      real cvchansn(nz,maxt),cvchansp(nz,maxt) ! convective region
c      real cvchanstotn(nz),cvchanstotp(nz) ! convective region

c      real cgchans(nz,maxt),cgchanstot(nz)  ! cg only 
c      real cgchansn(nz,maxt),cgchansp(nz,maxt) ! cg only 
c      real cgchanstotn(nz),cgchanstotp(nz) ! cg only 
c      real cginit(nz,maxt),cginittot(nz) ! convective region inits

c      real cgpchans(nz,maxt),cgpchanstot(nz)  ! cgp only 
c      real cgpchansn(nz,maxt),cgpchansp(nz,maxt) ! cgp only 
c      real cgpchanstotn(nz),cgpchanstotp(nz) ! cgp only 
c      real cgpinit(nz,maxt),cgpinittot(nz) ! convective region inits

c      real cgnchans(nz,maxt),cgnchanstot(nz)  ! cgn only 
c      real cgnchansn(nz,maxt),cgnchansp(nz,maxt) ! cgn only 
c      real cgnchanstotn(nz),cgnchanstotp(nz) ! cgn only 
c      real cgninit(nz,maxt),cgninittot(nz) ! convective region inits

c      real icchans(nz,maxt),icchanstot(nz)  ! ic only 
c      real icchansn(nz,maxt),icchansp(nz,maxt) ! ic only 
c      real icchanstotn(nz),icchanstotp(nz) ! ic only 
c      real icinit(nz,maxt),icinittot(nz) ! convective region inits

c      real anchans(nz,maxt),anchanstot(nz) ! anvil region (non-convective)
c      real anchansn(nz,maxt),anchansp(nz,maxt) ! anvil region (non-convective)
c      real anchanstotn(nz),anchanstotp(nz) ! anvil region (non-convective)

c      real drhop(nz,maxt), drhon(nz,maxt), drhomag(nz,maxt)
c      real drhototp(nz), drhototn(nz), drhotot(nz)
      
c      real cvdrhop(nz,maxt), cvdrhon(nz,maxt), cvdrhomag(nz,maxt)
c      real cvdrhototp(nz), cvdrhototn(nz), cvdrhotot(nz)
      
c      real andrhop(nz,maxt), andrhon(nz,maxt), andrhomag(nz,maxt)
c      real andrhototp(nz), andrhototn(nz), andrhotot(nz)

      real segheight(nz)
      real skyline(2*(nz+2),2)
        
         integer cginit
         parameter ( cginit = 11 ) 
         integer cgninit
         parameter ( cgninit = 13 )
         integer cgpinit
         parameter ( cgpinit = 14 )

c        Equivalence(chtotarr(1,1), chanstot(1))
c        Equivalence(chtotarr(1,2), chanstotp(1))
c        Equivalence(chtotarr(1,3), chanstotn(1))
c        Equivalence(chtotarr(1,4), cgchanstot(1))
c        Equivalence(chtotarr(1,5), cgchanstotp(1))
c        Equivalence(chtotarr(1,6), cgchanstotn(1))
c        Equivalence(chtotarr(1,7), icchanstot(1))
c        Equivalence(chtotarr(1,8), icchanstotp(1))
c        Equivalence(chtotarr(1,9), icchanstotn(1))
c        Equivalence(chtotarr(1,10), cinittot(1))
c        Equivalence(chtotarr(1,11), cginittot(1))
c        Equivalence(chtotarr(1,12), icinittot(1))
      
       real xgpos, ygpos

       INTEGER :: ibgn,jbgn,kbgn(nx,ny),kend(nx,ny)
c ###################################################################
c #################################################      
c Begin Executable Code      
c #################################################   

      IF ( .not. allocated( rh ) ) THEN
       allocate( rh(nx,ny,nz) )
       allocate( the(nx,ny,nz) )
       allocate( temc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
       allocate( plt3dp(nx,ny,nz) )
       allocate( plt3dn(nx,ny,nz) )
      ENDIF
      
      IF ( do_dualpol .and. .not. allocated( logZ ) ) THEN

      allocate( qr3d(nx,ny,nz),qs3d(nx,ny,nz),qg3d(nx,ny,nz),qh3d(nx,ny,nz) )
      allocate( Ntr(nx,ny,nz),Nts(nx,ny,nz),Ntg(nx,ny,nz),Nth(nx,ny,nz) )
      allocate( alphar3d(nx,ny,nz),alphas3d(nx,ny,nz),alphag3d(nx,ny,nz),alphah3d(nx,ny,nz) )
      allocate( qsw(nx,ny,nz),qgw(nx,ny,nz),qhw(nx,ny,nz) )
      allocate( logZ(nx,ny,nz),sumZh(nx,ny,nz),sumZv(nx,ny,nz),logZdr(nx,ny,nz) )
      allocate( logZdrrain(nx,ny,nz) )
      allocate( sumZhv(nx,ny,nz),Kdp(nx,ny,nz),Ahh(nx,ny,nz),Avv(nx,ny,nz),rhv(nx,ny,nz) )
      allocate( rhoms(nx,ny,nz),rhomg(nx,ny,nz),rhomh(nx,ny,nz) )
      allocate( rhograup(nx,ny,nz),rhohail(nx,ny,nz),tair(nx,ny,nz),rhoair(nx,ny,nz) )
      allocate( dpmask(nx,ny,nz) )
       
      ENDIF
      
      IF ( do_dualpol ) THEN
      qsw(:,:,:) = 0.0
      qgw(:,:,:) = 0.0
      qhw(:,:,:) = 0.0 
      ENDIF

       IF ( gridrelative ) THEN
         uava = ugrid
         vava = vgrid
         write(6,*) 'uava,vava = ',uava,vava
       ELSE
         uava = 0.0
         vava = 0.0
       ENDIF
         write(6,*) 'uava,vava = ',uava,vava
       
!       CALL GET_VARIABLE (gd, 'DX', dx)
!       CALL GET_VARIABLE (gd, 'DY', dy)
!       CALL GET_VARIABLE (gd, 'DZ', dz)
       
       IF ( iaddxy .eq. 1 ) THEN
         CALL GET_VARIABLE (gd, 'XG_POS', xgpos)
         CALL GET_VARIABLE (gd, 'YG_POS', ygpos)
       ELSE
         xgpos = 0.0
         ygpos = 0.0
       ENDIF

!       write(*,*) 'copy gx'
       DO i = 1,4
         gxt(:,i) = gx(i)%flt1d(:)
         gyt(:,i) = gy(i)%flt1d(:)
         gzt(:,i) = gz(i)%flt1d(:)
       ENDDO

       nxstag=nx-istag
       nystag=ny-jstag
       nzstag=nz-kstag

!       write(*,*) 'copy sinit'

      DO n = 1,2
        DO k = 1,nz
          sinit(k,n) = sinit0(n)%flt1d(k)
        ENDDO
      ENDDO

      ! write(0,*) 'copy sinit'

      DO kz = 1,nz-1
        db(kz) = 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sinit(kz,lt))
        pb(kz) = 1.0e5*piinit%flt1d(kz)**3.509

        temkz(kz) = sinit(kz,lt)*piinit%flt1d(kz) - 273.16
c       write(0,*) 'kz, temkz =',kz,temkz(kz),piinit%flt1d(kz),sinit(kz,lt),db(kz),pb(kz)
        DO jy=1,ny-1
          DO ix=1,nx-1
            pn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+
     :                      expi%flt3d(ix,jy,kz))**3.509 - pb(kz)
            dn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+
     :                       expi%flt3d(ix,jy,kz))**2.509/
     :                      (rd*sinit(kz,1)*(1.0+0.61*sinit(kz,2)))
            tt0(ix,jy,kz) = s(lt)%flt3d(ix,jy,kz)*
     :                      (piinit%flt1d(kz)+expi%flt3d(ix,jy,kz))
          ENDDO
        ENDDO
      ENDDO

       IF ( recalcdbz .ge. 1 ) THEN
         CALL RADARDD02(nx,ny,nz,ng,na,st    
     :       ,tt0,dbz%flt3d,dn, nz-1 ,cnoh,rho_qh, ipconc, 6, microp, 0, 0)
       ENDIF


      IF ( do_dualpol ) THEN
      rhoair(1:nx,1:ny,1:nz) = dn(1:nx,1:ny,1:nz)
        dpmask(:,:,:) = 0
        DO k = 1,nz
        DO j = 1,ny-1
        DO i = 1,nx-1
          IF ( dbz%flt3d(i,j,k) > 0 ) dpmask(i,j,k) = 1
        ENDDO
        ENDDO
        ENDDO

      IF ( microp(1:3) == 'TAK' ) THEN

           call setgrplhl(1,1,iusewetsnow)
            rhograup(i,j,k) = rho_qh_tak
            rhohail(i,j,k) = rho_qhl_tak
           qr3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lr)
           qs3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,ls)
           qg3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lh)
           qh3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lhl)


      ELSEIF ( microp(1:1) == 'Z' ) THEN

        qs3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,ls)
        IF ( lns > 1 ) THEN
           Nts(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lns)
        ENDIF

        IF ( lzr > 0 ) THEN
           sumZh(:,:,:) = rho_qr
           qr3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lr)
           Ntr(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnr)
           Ahh(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lzr)
           tmp = 1./3.
          write(0,*) 'call alpha iter for rain'
          call solve_alpha_iter(nx,ny,nz,rhoair,tmp,qr3d,
     :     ntr,ahh,sumZh,alphar3d)
          write(0,*) 'done call alpha iter for rain'
           sumZh(:,:,:) = 0.0
        ELSE
           qr3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lr)
           IF ( lnr > 1 ) Ntr(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnr)
          alphar3d(:,:,:) = alphar
        ENDIF

        alphas3d(:,:,:) = -0.4 ! alphas

        i = 0
        j = 0
        IF ( lh > 1 ) i = 1
        IF ( lhl > 1 ) j = 1
        call setgrplhl(i,j,iusewetsnow)
        call model_dsd(MPflg,cnor,cnos,cnohl,cnoh,rho_qs,rho_qhl,rho_qh)
        call set_nbins(0)
        
        IF ( lh > 1 .and. lnh > 1 ) THEN
           qg3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lh)
           Ntg(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnh)
        ENDIF

        IF ( lhl > 1 .and. lnhl > 1 ) THEN
           qh3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lhl)
           Nth(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnhl)
        ENDIF

         IF ( lsw > 1 ) THEN
           qsw(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lsw)
         ENDIF

         IF ( lhw > 1 ) THEN
           qgw(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lhw)
         ENDIF

         IF ( lhlw > 1 ) THEN
           qhw(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lhlw)
         ENDIF


       IF ( lh > 1 ) THEN
           qg3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lh)
           Ntg(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnh)
        DO k = 1,nz
        DO j = 1,ny-1
        DO i = 1,nx-1
         qg3d(i,j,k) = st(i,j,k,lh)
         rhograup(i,j,k) = 0.0
         IF ( st(i,j,k,lh) .gt. 1.e-12 ) THEN
            rhograup(i,j,k) = rho_qh
          IF ( lvh > 1 ) THEN
           IF ( st(i,j,k,lvh) .gt. 1.e-28  ) THEN
             rhograup(i,j,k) = dn(i,j,k)*st(i,j,k,lh)/st(i,j,k,lvh)
             rhograup(i,j,k) = Min( 1000., Max( 50., rhograup(i,j,k) ) )
           ENDIF
          ENDIF
         ENDIF
          
        rhomg(i,j,k)  = rhograup(i,j,k)
        
        ENDDO
        ENDDO
        ENDDO
        ENDIF
 
        IF ( lzh > 0 ) THEN
           Ahh(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lzh)
          write(0,*) 'call alpha iter for graupel'
          call solve_alpha_iter(nx,ny,nz,rhoair,1./3.,qg3d,
     :     ntg,ahh,rhograup,alphag3d)
          write(0,*) 'done call alpha iter for graupel'
        ELSE
          alphag3d(:,:,:) = alphah
        ENDIF

        qh3d(:,:,:) = 0.0

        
        IF ( lhl > 1 .and. lnhl > 1 ) THEN
           qh3d(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lhl)
           Nth(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lnhl)
        ENDIF
        rhomh(:,:,:) = rho_qhl

        IF ( lhl > 1 ) THEN


        DO k = 1,nz
        DO j = 1,ny
        DO i = 1,nx

         qh3d(i,j,k) = st(i,j,k,lhl)
         rhohail(i,j,k) = 0.0
         IF ( st(i,j,k,lhl) .gt. 1.e-12 ) THEN
            rhohail(i,j,k) = rho_qhl
          IF ( lvhl > 1 ) THEN
           IF ( st(i,j,k,lvhl) .gt. 1.e-28  ) THEN
             rhohail(i,j,k) = dn(i,j,k)*st(i,j,k,lhl)/st(i,j,k,lvhl)
             rhohail(i,j,k) = Min( 1000., Max( 50., rhohail(i,j,k) ) )
           ENDIF
          ENDIF
         ENDIF
          
        rhomh(i,j,k)  = rhohail(i,j,k)

        ENDDO
        ENDDO
        ENDDO
        ENDIF
 
        IF ( lhl > 1 .and. lzhl > 0 ) THEN
           Ahh(1:nx,1:ny,1:nz) = st(1:nx,1:ny,1:nz,lzhl)
           
          call solve_alpha_iter(nx,ny,nz,rhoair,1./3.,qh3d,
     :     nth,ahh,rhohail,alphah3d)
        ELSEIF ( lhl > 1 ) THEN
          alphah3d(:,:,:) = alphahl        
        ENDIF

      
      sumZh(:,:,:) = 0.0
      sumZv(:,:,:) = 0.0
      logZdr(:,:,:) = 0.0
      logZdrrain(:,:,:) = 0.0
      sumZhv(:,:,:) = 0.0
      Kdp(:,:,:) = 0.0
      Ahh(:,:,:) = 0.0
      Avv(:,:,:) = 0.0
      rhv(:,:,:) = 0.0
      rhoms(:,:,:) = rho_qs
      tair(1:nx,1:ny,1:nz) = tt0(1:nx,1:ny,1:nz)
      rhoair(1:nx,1:ny,1:nz) = dn(1:nx,1:ny,1:nz)
      
      DO kz = nz-1,1,-1
        write(0,*) 'kz,tair = ',kz,tair(nx/2,ny/2,kz)
      ENDDO
      
!      logZ(1:nx,1:ny,1:nz) = dbz%flt3d(1:nx,1:ny,1:nz)
!      logZ(1,1,nz) = 1.e11
      
      IF ( dualpol_fulldomain ) THEN ! calculate full domain
       ibgn = 1
       iend = nx-1
       jbgn = 1
       jend = ny-1
       kbgn(:,:) = 1
       kend(:,:) = nz-1
      
       call refl_rsa_array(addtorain,MPflg,MFflg,nx,ny,nz,dirscatt,wavelen,ibgn,iend,jbgn,jend,kbgn,kend,  
     :      rhoair,qr3d,qs3d,qg3d,qh3d,Ntr,Nts,Ntg,Nth,alphar3d,alphas3d,alphag3d,alphah3d,
     :      qsw,qgw,qhw,rhograup,rhohail, 
     :      logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh,tair,
     :       dpmask, dbz%flt3d(1:nx,1:ny,1:nz) )
        
!        write(0,*) 'k,logz,logzdr,rhv,alpr,alph,alphl,qgw,qhw'
!        i = 110
!        j = 92
!        DO k = nz-1,1,-1
!          write(0,*) k,logz(i,j,k),logzdr(i,j,k),rhv(i,j,k),alphar3d(i,j,k),alphag3d(i,j,k),
!     :           alphah3d(i,j,k),qg3d(i,j,k),st(i,j,k,lh),qgw(i,j,k),qhw(i,j,k),qr3d(i,j,k),st(i,j,k,lr)
!        ENDDO
       
       ENDIF
       
      ENDIF ! microp = Z
      ENDIF ! do_dualpol

      ! write(0,*) 'finished setup'

c       lt=1
c       le=2
c       lv=3
c       lc=4
c       lr=5
c       li=6
       
c       IF ( na .le. 42 ) THEN
c         ipconc = 0
c       ELSEIF ( na .eq. 46 ) THEN
c         ipconc = 3
c       ELSEIF ( na .eq. 47 ) THEN
c         ipconc = 4
c       ELSEIF ( na .eq. 53 ) THEN
c         ipconc = 5
c       ENDIF

c       IF ( na .ge. 40 ) THEN
c       iqi=4
c       iqc=5
c       iqr=6
c       iqs=7
c       iqg=8
c       iqh=9
c       isci=10
c       isccw=11
c       iscr=12
c       iscs=13
c       iscg=14
c       isch=15 
c        
c        lscwi = 17 
c       lscw=0+lscwi
c       lscr=1+lscwi
c       lsci=2+lscwi
c       lscip=3+lscwi
c       lscir=4+lscwi
c       lscs=5+lscwi
c       lscgl=6+lscwi
c       lscgm=7+lscwi
c       lscgh=8+lscwi
c       lscf=9+lscwi
c       lsch=10+lscwi
c       lschl=11+lscwi
c
c         lv=3
c         lc=4
c         lr=5
c         li=6
c         lip=7
c         lir=8
c         ls=9
c         lgl=10
c         lgm=11
c         lgh=12
c         lf=13
c         lh=14
c         lhl=15
c         lhab=15
c         
c         lnshi = 42
c         lnip=2+lnshi
c         lnr=4+lnshi
c         lns=5+lnshi
c         lngl=6+lnshi
c         lngm=7+lnshi
c         lngh=8+lnshi
c         lnf=9+lnshi
c         lnh=10+lnshi
c         lnhl=11+lnshi
c        
c        lqb=lc
c        lqe=lhl 
c
c       lscb=lscwi
c       lsce=lschl
c       lsceq = lsce
c       
c       ELSEIF ( na .ge. 8 ) THEN
c       iqi=4
c       iqc=5
c       iqr=6
c       iqs=7
c       iqg=8
c       iqh=9
c
c       isci=10
c       isccw=11
c       iscr=12
c       iscs=13
c       iscg=14
c       isch=15 
c
c       lir=na+1
c       ls=7
c       lgl=na+1
c       lgm=na+1
c       lgh=na+1
c       lf=na+1
c       lh=8
c       lip=na+1
c       lhl=na+1
c       lhab=lh
c        
c        lscwi = na+1 ! 17 
c       lscw=0+lscwi
c       lscr=0+lscwi
c       lsci=0+lscwi
c       lscir=0+lscwi
c       lscs=0+lscwi
c       lscgl=0+lscwi
c       lscgm=0+lscwi
c       lscgh=0+lscwi
c       lscf=0+lscwi
c       lsch=0+lscwi
c       lscip=00+lscwi
c       lschl=0+lscwi
c        
c        lqb=lc
c        lqe=lh 
c
c       lscb=lscwi
c       lsce=lschl
c        IF ( na .ge. 15 ) THEN
c        lscwi = 9
c       lscw=0+lscwi
c       lscr=1+lscwi
c       lsci=2+lscwi
c       lscs=3+lscwi
c       lsch=4+lscwi
c
c       lscb=lscwi
c       lsce=lsch
c        
c        ENDIF
c       ENDIF
c
c       lsceq = lsce
c
      
      ion = 0
      
      IF ( lscpi >= 1 ) ion = 1
      IF ( lscpli >= 1 ) ion = 2

      
      imud = 1
      xmax4 = 0.0
c      ifirst = 0
      ieps = 1
      epsflag = 0
      TIME = 0
      NUMTIME = 0
      ntime = 0
        
        write(6,*) 'Start! ifirst = ', ifirst
        
        imudstop = 0
        feflst = .false.

c
c
c  open input data file
c  this file has parametres and file names
c
c      open(22,file='inpltba',status='old')
c read from stdin
      lunr = iunit
      lunw = 23
      
c      indx = iargc( )
      
c      IF ( indx .ge. 1 ) THEN
c        call getarg ( 1, inputfile )
c        write(0,'(a,a)') 'opening file: ',inputfile
c        lunr = 55
c        open(unit=lunr,file=inputfile,form='formatted',status='old')
c        IF ( indx .gt. 1 ) THEN
c         write(0,*) 'Only reads first file! One argument, please!'
c        ENDIF
c      ELSE
c      
c      ENDIF


 909  CONTINUE
       ifirsteps = 0

      IF ( lunr .ne. 5 ) THEN
        rewind(lunr)
      ENDIF
      
      IF ( ifirst .eq. 0 ) THEN
c        read(lunr,*) lunr1,lunw1
c        write(*,*) lunr1,lunw1
c
        open(23,file='pltinfo',form='formatted',status='unknown')
      ELSE
       read(lunr,'(a72)') line
      ENDIF
c
c  read file name of data, and ploting and scaling parameters
c
c 9999 continue
      write(*,*) 
     > 'INPUT: imore, imore>   0 (plot more data after this one)'
      write(*,*) '             , imore=  0 (plot 1 data set)'
      write(*,*) '             , imore<  0 (stop)'
c      read(lunr,*) imore
c
      write(*,*) 'INPUT: simulation title'
c      read(lunr,*) stitle
      write(*,*) 'title = ', stitle
      titlen = len(stitle)
      CALL strlnth(stitle,titlen)
      
c      write(*,*) 'INPUT: in file name with data to be plotted'
c      read(lunr,*) filnm1,ihdmpt,icomp,imicro
c      read(lunr,*) filnm2
c      read(lunr,*) filnm3,filnm4
c      read(lunr,*) filnm5
        icomp = 0
        
c        read(lunr,*) outname
c        write(6,*) 'outname= ',  outname     
        

c        read(lunr,*) imsz  ! flag for an msz flash
c        write(6,*) 'imsz = ',imsz 
c        read(lunr,*) ifile2
c        write(6,*) 'ifile2 = ',ifile2
c          read(lunr,*) basename
c          write(6,*) 'basename = ',basename   ! *base file
c        read(lunr,*) bname
c        write(6,*) 'bname = ',bname  ! base filename *nelec* file
        
c          read(lunr,*) inname2
c          write(6,*) 'inname2 = ',inname2   ! *ht1pa* file
c          read(lunr,*) inname4
c          write(6,*) 'inname4 = ',inname4  ! *mica* file 
c          read(lunr,*) inname3
c          write(6,*) 'inname3 = ',inname3   ! *mic* file

c      read(lunr,'(a)') line
      IF ( ifirst .eq. 0 ) THEN
        inameps = len(outname)
        CALL strlnth(outname,inameps)
c        inameps = inameps+3+4

c        read(line,*) time,nliter
c        ntime = 0
c        IF ( time .lt. 0 ) THEN
c          read(line,*) time,nliter,dtdump
c          time0 = Abs(time)
c          time = time0
c          numtime = nliter
c          nliter = 1
c        ELSE
c          dtdump = 0.0
c          numtime = 1
c        ENDIF
c      ELSE
       nliter = 1
      ENDIF
     
      dzin = dzplot
      nulbll = 0
      ioffp  = 0
      spval  = 0.0
      dzadd  = 0.0
c      read(lunr,*) umove,vmove,dzadd  ! ,timew
c      write(*,*) umove,vmove,dzadd ! ,timew

c
      xfacb0 = xfacb
      xdivb0 = xdivb
      
      
c      IF ( iotyp .eq. 20 ) THEN
c      Open(unit=9,file=outname,form='formatted',status='unknown')
c      ELSE
c      Open(unit=9,file=outname,form='unformatted',status='unknown')
c      ENDIF

      
      pi = 4.0*Atan(1.0)

c
c  open gks
c
c      call opngks
C
C  Open "workstation" for NCAR GRAPHICS 
C
        IF ( ifirst .eq. 0 ) THEN
        ifirst =  1
        ifirsteps = 1
        !write(0,*) 'call gopks'
        CALL GOPKS(6,0)
        IF ( iotyp .eq. 0 ) THEN
          iotyp = 1
          igmname = 1
           write(outnameps,'(a,a)') 
     :             outname(1:inameps),'.ncgm'
            CALL ngsetc('ME',outnameps)
        ENDIF
        IF ( iotyp .ne. 1 ) THEN
          IF ( iotyp .eq. 20 ) THEN
           write(outnameps,'(a,a)') 
     :             outname(1:inameps),'.ps'
            CALL ngsetc('ME',outnameps)
          ELSE
           write(outnameps,'(a,a,i2.2,a,i3.3,a)') 
     :             outname(1:inameps),'.',ieps,'.',iframe,'.eps'
            CALL ngsetc('ME',outnameps)
c            write(0,*) 'outnameps: ',outnameps

          ENDIF
        ENDIF
        CALL GOPWK(1,9,iotyp)  ! ncgm output
        CALL GACWK(1)
        call agseti('FRAME.',2)  !  set no automatic frame advance; allows overlays
         IF ( iotyp .ne. 1 ) THEN
           call ngseti('WO',1)
           call ngseti('PA',10000)
           call ngsetr('NO',2.)
         ENDIF

c        ELSEIF (ifirsteps .eq. 0 .and. iotyp.eq.21 ) THEN
        
        ELSEIF ( epsflag .eq. 1 ) THEN
          
          epsflag = 0
          
        CALL GOPKS(6,0)
        write(outnameps,'(a,a,i2.2,a,i3.3,a)') 
     :             outname(1:inameps),'.',ieps,'.',iframe,'.eps'
        CALL ngsetc('ME',outnameps)
        write(0,*) ieps,'outnameps: ',outnameps
        CALL GOPWK(1,9,iotyp)  ! eps output
        CALL GACWK(1)
        call agseti('FRAME.',2)  !  set no automatic frame advance; allows overlays
         call ngseti('WO',1)
           call ngseti('PA',10000)
         call ngsetr('NO',2.)

c
c        set some parameters
c
         do i = 1,129
         ifill(i) = i
         end do

         CALL setclr
         
         CALL GSTXFP (-29,0)  ! courier font, string precision

        ELSEIF ( ( ifirst > 0 .and. iotyp == 21 ) ) THEN
        CALL GOPKS(6,0)
           write(outnameps,'(a,a,i2.2,a,i3.3,a)') 
     :             outname(1:inameps),'.',ieps,'.',iframe,'.eps'
        write(0,*) ieps,'outnameps: ',outnameps

            CALL ngsetc('ME',outnameps)
        CALL GOPWK(1,9,iotyp)  ! ncgm output
        CALL GACWK(1)
        call agseti('FRAME.',2)  !  set no automatic frame advance; allows overlays
         IF ( iotyp .ne. 1 ) THEN
           call ngseti('WO',1)
           call ngseti('PA',10000)
           call ngsetr('NO',2.)
         ENDIF
        
        ENDIF
C
C Set all the GKS aspect source flags to "individual".
C
c      CALL GSASF (IASF)
c
c
c  set some parameters
c
      do i = 1,129
      ifill(i) = i
      end do
c
      CALL setclr
        IF ( iotyp .ne. 1 ) THEN
          CALL GSTXFP (-29,0)  ! times-roman font, string precision
        ENDIF
c
c  over ride plots (not in use yet)
c
c     read(*,*) ifvint, fhiv, iflengthv
c     write(*,*) 'ifvint, fhiv, iflengthv'
c     write(*,*) ifvint, fhiv, iflengthv
c
      ibase = ibase_thproc
      ibase_mic_add = ibase_micth-24
      CALL setplottitles(ntit,title,ibase,ibase_micth,ibase_mic_add)
c  constants
c
      dxi = 1.0/dx
      dyi = 1.0/dy
      dzi = 1.0/dzin
      dxi2 = 0.5*dxi
      dyi2 = 0.5*dyi
      dzi2 = 0.5*dzi
      
 
 9091  CONTINUE
     
c      DO kz=-ng+1,nz+ng
c      DO jy=-ng+1,ny+ng
c      DO ix=-ng+1,nx+ng
        rho (:,:,:) = 0.0
c        flashm (ix,jy,kz) = 0.0
c        flashp (ix,jy,kz) = 0.0
c        flashn (ix,jy,kz) = 0.0
c        drhom (ix,jy,kz) = 0.0
c        scionp (ix,jy,kz) = 0.0
c        scionn (ix,jy,kz) = 0.0
c      ENDDO
c      ENDDO
c      ENDDO
      
      write (6,*) 'zero!'

      DO kz=1,nz
      DO jy=1,ny
      DO ix=1,nx
c        dat (ix,jy,kz) = 0.0
c        chan (ix,jy,kz) = 0.0
      ENDDO
      ENDDO
      ENDDO
      
c
c  conrec and velvct parameters
c    set ioffm=1 to turn off conrec message
c
      ioffm = 1
      nsetc = -1 
      flov = 0.
      hiv = 20.
      nsetv = -1
      lengthv = 60
      ispv = 4
      spv(1) = 999.
      spv(2) = 999.
      nplot = 0
c
c 1999 continue
c
c
      if ( ijump .ge. 1 ) go to 777
      
c      itime = Int(time)
      write(stime,'(i6.6)') itime
      write(stimeplot,'(i6)') itime + itimeoffset
      write(smin,'(f7.2)') float(itime + itimeoffset)/60.0
      write(smin2,'(i4)') Nint((itime + itimeoffset)/60.0)
      write(shr,'(f7.2)') float(itime + itimeoffset)/3600.0
      xtitle = stitle(1:titlen) ! //' ( '//stime//' s)'
      xtitlen = titlen  ! +5+6

      IF ( .false. ) THEN !{
c
c  read raw data
c
c check for flash file.  Use name of ht1 file to create the filename
c

      idum = len(bname)
      CALL strlnth(bname,idum)
      flstname = bname(1:idum)//'.flst.'//stime
      write(6,*) flstname

      INQUIRE(file=flstname, exist = feflst )
       IF ( .not. feflst ) THEN
         write(6,*) 'file not found: ',flstname
       END IF      

        idum = len(bname)
        CALL strlnth(bname,idum)
c nelec name        
        write(inname,'(a,i2.2,a)') 
     :   bname(1:idum)//'.nelec.',nliter,'.'//stime
        write(6,*) inname
c ht1pa name
        inname2 = bname(1:idum)//'.ht1pa.'//stime
        write(6,*) inname2
c mica file        
        inname4 = bname(1:idum)//'.mica'//stime
        write(6,*) inname4
c mic file        
        inname3 = bname(1:idum)//'.mic'//stime
        write(6,*) inname3
c sll file
        insll = bname(1:idum)//'.hlgttyphsll'
        write(6,*) insll

c
c read nelec file
        ifile = 10   
      iscrn = 0
        
      INQUIRE(file=inname, exist = fenele )
      IF( .not. fenele ) THEN
        write(6,*) 'Could not open file ',inname
        
c        IF ( dtdump .gt. 0.0 ) THEN
c          GOTO 9010 ! Stop
c        ELSE
          GOTO 1010
c        ENDIF
      END IF

      write(6,*) 'reading file, inname = ',inname
         
        open(unit=ifile,file=inname,status='old',
     :       form='unformatted',iostat= istat )
        
c        read(ifile) sccna
        read(ifile) neletyp
        
        write(6,*) 'neletyp = ',neletyp

        IF ( neletyp .eq. -1 ) THEN
          iscrn = 1
c          read(ifile) scrnchg
c          read(ifile) neletyp
        ENDIF
        
        IF ( neletyp .eq. 0 ) THEN 
           ichanflag = 0
        ELSEIF ( neletyp .eq. 1 ) THEN  !  flash info for regular dump
          ichanflag = 1
        write(6,*) 'got here 0'
c         CALL read1ele(nx,ny,nz,ng,ifile,
c     :        flash,drho )
         
        write(6,*) 'got here 1'
        ELSEIF ( neletyp .eq. 2 ) THEN ! stand-alone flash dump
c
c NOTE that 'rho' is the preflash net space charge density and
c  that the other charge densities are _post-flash_
c         
         ichanflag = 1
         write(6,*) 'neletyp=',neletyp,'  Read charge etc.'
         
c         CALL read2ele(nx,ny,nz,ng,na,a,ifile,
c     :        flash,drho,rho,
c     :        iqi,iqc,iqr,iqs,iqg,iqh,
c     :        isci,isccw,iscr,iscs,iscg,isch )
        
c        read(ifile) rho ! rho (space charge)
c        read(ifile) drho ! scadj (drho)
c        read(ifile) flash ! flash
c        read(ifile) qi ! cld ice tot (qci + qip + qir)
c        read(ifile) pre ! pre
c        read(ifile) w  ! w
c        read(ifile) qc
c        read(ifile) chan
        
        END IF
        close (ifile)
        

 1010  CONTINUE
      
        fshortht = .false.      
c      INQUIRE(file=inname2, exist = feht1 )
c      IF( .not. feht1 .and. ifile2 .ge. 1 ) THEN
c        write(6,*) 'file not found: ',inname2
c        GOTO 1020
c      END IF
c      IF (ifile2 .ge. 1) THEN
cc read *ht1pa* file
c        
c        write(6,*) 'reading file ',inname2
c        ifile = 13   
c        open(unit=ifile,file=inname2,status='old',
c     :       form='unformatted',iostat= istat )
c        
c        read(ifile) temc  ! temperature (kelvin)
c        read(ifile) qv  ! vapor
c        read(ifile) v   ! 3-D winds (grid-relative)
c        read(ifile) p   ! pressure perturbation
c        
c
c        read(ifile,end=1015) 
c     >                   rwtot2d,rwrt2d
c     >                  ,swtot2d,swrt2d,gltot2d,glrt2d
c     >                  ,gmtot2d,gmrt2d,ghtot2d,ghrt2d
c     >                  ,fwtot2d,fwrt2d,hwtot2d,hwrt2d
c     >                  ,prt2d,ptot2d,prti2d,ptoti2d
c     >                  ,prtl2d,ptotl2d
c        read(ifile,end=1015) theta
c        read(ifile,end=1015) tke
c
c        read(ifile,end=1016) scionp
c        read(ifile,end=1016) scionn
c        
c        ion = 1
c        
c        GOTO 1018
c 1015   CONTINUE
c         fshortht = .true.
c       write(6,*) 'premature EOF in ht1 file'
c 1016   CONTINUE
c         write(6,*) 'No ions in ht1 file!'
c          
c         ion = 0
c         
c 1018   CONTINUE       
c        DO kz=-ng+1,nz+ng
c        DO jy=-ng+1,ny+ng
c        DO ix=-ng+1,nx+ng
c          st(ix,jy,kz,lv) = qv(ix,jy,kz)
c          st(ix,jy,kz,lt) = theta(ix,jy,kz)
c          st(ix,jy,kz,le) = tke(ix,jy,kz)
c          IF (ix .eq. nx ) THEN
c           temc(ix,jy,kz) = temc(ix-1,jy,kz)
c          ELSEIF ( jy .eq. ny ) THEN
c           temc(ix,jy,kz) = temc(ix,jy-1,kz)
c          ENDIF
c        END DO
c        END DO
c        END DO
c        
c        close (ifile)
c        write(6,*) 'done reading file ',inname2
c        
c      END IF
c 1020  CONTINUE
        
c      INQUIRE(file=inname4, exist = femica )  
c       IF ( .not. femica .and. ifile2 .ge. 1) THEN
c         write(6,*) 'file not found: ',inname4
c       END IF      
c      INQUIRE(file=inname3, exist = femic )
c       IF ( .not. femic .and. ifile2 .ge. 1 ) THEN
c         write(6,*) 'file not found: ',inname3
c       END IF      
c
c      INQUIRE(file=insll, exist = fesll )
c       IF ( .not. fesll .and. ifile2 .lt. 0 ) THEN
c         write(6,*) 'file not found: ',insll
c       END IF      
c
c      INQUIRE(file=basename, exist = febase )
c       IF ( .not. febase  ) THEN
c         write(6,*) 'file not found: ',basename
c       END IF    
c       
c       IF ( .not. (femica .or. femic .or. feht1 .or. fenele)) THEN
c        IF ( dtdump .gt. 0.0 ) THEN
c          GOTO 9010 ! Stop
c        ENDIF
c       ENDIF
c         
c
c      IF (ifile2 .ge. 2 .and. femica) THEN
c        write(6,*) 'reading file ',inname4
c       IF ( na .ge. 40 ) THEN
c       CALL ba2read(nx,ny,nz,ng,na,a,inname4,ipconc)
c       ELSE
c         CALL lfo2read(nx,ny,nz,ng,na,a,inname4)
c       ENDIF
c        IF ( ifile2 .eq. 2 .or. .not. femic) THEN
c          IF ( na .ge. 40 ) THEN
c            CALL ba2read(nx,ny,nz,ng,na,ab,inname4,ipconc)
c          ELSE
c            CALL lfo2read(nx,ny,nz,ng,na,ab,inname4)
c          ENDIF
cc          CALL ba2read(nx,ny,nz,ng,na,ab,inname4)
c        END IF
c       write(6,*) 'done reading file ',inname4
c      END IF
c
c      IF (ifile2 .ge. 3 .and. femic) THEN
c        write(6,*) 'reading file ',inname3
c       CALL baread(nx,ny,nz,ng,na,ab,tema,-10,inname3,ipconc)
c       write(6,*) 'done reading file ',inname3
c      END IF
c      
c        uava = 0.0
c        vava = 0.0
c
c      IF ( febase ) THEN
c        write(6,*) 'reading file ',basename
c        ifile = 15  
c        open(unit=ifile,file=basename,status='old',
c     :       form='unformatted',iostat= istat )
c        
c        rewind ifile
cc        read(ifile) thetab,tkeb,qvb,db,pb,vb,gt,ht
cc       read(ifile) umv,vmv
c        read(ifile) umv1,vmv1,dx1,dy1,dz1
c        
c        
c        IF ( dz1 .eq. 0.0 ) THEN
c          write(0,*) 'dx,dy,dz missing from base file:',
c     :         ' using input values for dx and dy'
c         dzmax = gt(1,1,nz,mzdist) - gt(1,1,nz-1,mzdist)
c         write(6,*) 'dzmax = ',dzmax
c
c         dz = dzmax*gt(1,1,nz,imapz)
c         write(0,*) 'setting dz = ',dz
c       ELSE
c          write(0,*) 'dx,dy,dz read from base file:',dx1,dy1,dz1
c          dz = dz1
c          dx = dx1
c          dy = dy1
c          dzmax = dz/gt(1,1,nz-1,imapz)
c        ENDIF
c        
c        DO kz=-ng+1,nz+ng
c        DO jy=-ng+1,ny+ng
c        DO ix=-ng+1,nx+ng
c          ab(ix,jy,kz,lv) = qvb(ix,jy,kz)
c          ab(ix,jy,kz,lt) = thetab(ix,jy,kz)
c          ab(ix,jy,kz,le) = tkeb(ix,jy,kz)
c        END DO
c        END DO
c        END DO
c        
c      dzsave = dz
c        DO kz=1,nz
c         z1d2(kz,1) = heights(kz)
c        ENDDO
c      DO k=1,nz
c        z1d2(k,3) = gt(1,1,k,imapz)
c        z1d2(k,4) = 0.5*(gt(1,1,k-1,imapz) + gt(1,1,k,imapz))
c        z1d2(k,1) = gt(1,1,k,mzdist)
c        z1d2(k,2) = gt(1,1,k,mzdist)-0.5*dz/gt(1,1,k,imapz)
c      ENDDO

      ENDIF ! .false. !}
      
      dzsave = dz
      
      DO k=1,nz-1
        z1d2(k,3) = (dz*gzt(k,3)) ! gt(1,1,k,imapz)    ! MFC
        z1d2(k,4) = (dz*gzt(k,4)) ! 0.5*(gt(1,1,k-1,imapz) + gt(1,1,k,imapz))  ! MFE
        z1d2(k,1) = gzt(k,1) ! gt(1,1,k,mzdist)   ! scalar height
        z1d2(k,2) = gzt(k,2) ! gt(1,1,k,mzdist)-0.5*dz/gt(1,1,k,imapz)  ! W-height
      ENDDO
        dzmax = 1.0/gzt(nz-1,3)
        z1d2(nz,2) = z1d2(nz-1,2) + dzmax
c
c  create templt 
c


      templt(1) = temkz(1)
      do kz = 1,nz
        z = (kz - 1)*dzin
        zk(kz) = z
       DO k = 1,nz
         IF ( z1d2(k,1) .le. z .and. z1d2(k+1,1) .gt. z ) THEN
          tmp = (z-z1d2(k,1))*
     :      (temkz(k+1) - temkz(k))/(z1d2(k+1,1)-z1d2(k,1))
          
          templt(kz) = temkz(k) + tmp
          ENDIF
        ENDDO
        
c       write(0,*) kz,z,temkz(kz),templt(kz)
      end do

        IF ( gzt(1,3) .eq. gzt(nz,3) ) THEN
         istrz = 0
        ELSE
         istrz = 1
        ENDIF
        
c        close(ifile)
c      END IF

        IF ( ifile2 .le. -1 ) THEN
         segname = thfile !  bname

         DO jy=1,maxc
         DO ix=1,maxt
         DO kz=1,nz
           charr(kz,ix,jy) = 0.0
         ENDDO
         ENDDO
         ENDDO
         
         IF ( ifile2 .eq. -1 ) THEN
         idx0 = 900
         ELSEIF ( ifile2 .eq. -2 ) THEN
         idx0 = 1000
         ELSEIF ( ifile2 .eq. -3 ) THEN
         idx0 = 950
         ELSEIF ( ifile2 .eq. -4 ) THEN  ! radar time-height (and updraft vol. time-height)
         idx0 = 990
         ELSEIF ( ifile2 .eq. -5 ) THEN
         idx0 = 980
         ELSEIF ( ifile2 .eq. -6 ) THEN ! micro (hydrometeor) time-height
         idx0 = 1100
         ELSEIF ( ifile2 .eq. -7 ) THEN ! micro process time-height
         idx0 = 1020
         ELSE
          write(0,*) 'ifile2 .lt. -7 not allowed! Stopping'
          STOP
         ENDIF

      Open(unit=25,file=segname,form='formatted',status='old')
      
      read(25,'(a)') dumstr
      itcount = 0

      icgp = 0
      icgn = 0
      icg = 0

      IF ( ifile2 .eq. -1 ) THEN
      DO it=1,maxt

          DO kz=1,nz-1
          read(25,'(a)',end=600) dumstr
c          
      istat1 = 0
!      read(dumstr,'(f8.1,1x,f6.2,14(1x,f7.0))',iostat=istat1)
!     :       timearr(it), heights(kz), 
!     :     (charr(kz,it,i), i=1,14)

      read(dumstr,*,iostat=istat1)
     :       timearr(it), heights(kz), 
     :     (charr(kz,it,i), i=1,14)
c       write(6,*) 'ltg: ',timearr(it), heights(kz), charr(kz,it,1),charr(kz,it,10)

c     :       chans(kz,it), chansp(kz,it), chansn(kz,it), 
c     :       cgchans(kz,it), cgchansp(kz,it), cgchansn(kz,it), 
c     :       icchans(kz,it), icchansp(kz,it), icchansn(kz,it), 
c     :       cinit(kz,it), 
c     :       cginit(kz,it), icinit(kz,it),
c     :       cgninit(kz,it),cgpinit(kz,it)
         icg = icg + charr(kz,it,cginit)
         icgn = icgn + charr(kz,it,cgninit)
         icgp = icgp + charr(kz,it,cgpinit)

c         IF ( istat1 .ne. 0 ) THEN
c      read(dumstr,'(f8.1,1x,f6.2,12(1x,f7.0))',iostat=istat1)
c     :     timearr(it), heights(kz), 
c     :     (charr(kz,it,i), i=1,12)

c     :     chans(kz,it), chansp(kz,it), chansn(kz,it), 
c     :     cgchans(kz,it), cgchansp(kz,it), cgchansn(kz,it), 
c     :     icchans(kz,it), icchansp(kz,it), icchansn(kz,it), 
c     :     cinit(kz,it), 
c     :     cginit(kz,it), icinit(kz,it)
         
c         ENDIF
         END DO

!       DO kz=50,1,-1
!        write(6,'(i3,50(1x,i2))') it+39, (Int(charr(i,it,10)), i=1,50)
!       ENDDO
       
c       write(6,*) 'it,timearr= ',it,timearr(it)

       itcount = itcount + 1

      END DO
      ELSEIF ( ifile2 .eq. -2 ) THEN
! Time, Altitude, udmf, uicmf, grvol, hlvol, grms, hlms, ticms, rnms, snms, iwcmx, wmax, wmin

       DO it=1,maxt
      
         IF ( dumstr(89:92) .eq. 'wmin' ) THEN
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     (charr(kz,it,i), i=1,12)
          ENDDO
         ELSE
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     (charr(kz,it,i), i=1,9)
c     charr(kz,it,2), charr(kz,it,3), 
c     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6), 
c     :     charr(kz,it,7), charr(kz,it,8), charr(kz,it,9) 
c     :     ,charr(kz,it,10), 
c     :     charr(kz,it,11), charr(kz,it,12)
cc     :     drhop(kz,it), drhon(kz,it), drhomag(kz,it)
cc     :     ,cvdrhop(kz,it), cvdrhon(kz,it), cvdrhomag(kz,it),
cc     :     andrhop(kz,it), andrhon(kz,it), andrhomag(kz,it) 
        
           ENDDO
         ENDIF
        itcount = itcount + 1
       END DO

      ELSEIF ( ifile2 .eq. -3 ) THEN

      DO it=1,maxt
      
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     charr(kz,it,1), charr(kz,it,2), charr(kz,it,3), 
     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6), 
     :     charr(kz,it,7), charr(kz,it,8) ! , charr(kz,it,9) 
c     :     ,charr(kz,it,10), 
c     :     charr(kz,it,11), charr(kz,it,12)
cc     :     drhop(kz,it), drhon(kz,it), drhomag(kz,it)
cc     :     ,cvdrhop(kz,it), cvdrhon(kz,it), cvdrhomag(kz,it),
cc     :     andrhop(kz,it), andrhon(kz,it), andrhomag(kz,it) 
        
        ENDDO
       itcount = itcount + 1
      END DO

      ELSEIF ( ifile2 .eq. -4 ) THEN

          read(25,'(a)') dumstr
          it = 1
          kz = 1
         istat1 = 0
          read(dumstr,*,iostat=istat1)  timearr(it), heights(kz),
     :     charr(kz,it,1), charr(kz,it,2), charr(kz,it,3), 
     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6) 

      IF ( istat1 .eq. 0 ) THEN
      write(6,*) '*** Reading updraft volumes ***'
      rewind (25)
      read(25,'(a)') dumstr
      DO it=1,maxt
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     charr(kz,it,1), charr(kz,it,2) , charr(kz,it,3), 
     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6) 
        
        ENDDO
       itcount = itcount + 1
      END DO
      
      ELSE  ! old file with reflectivity only
      
      write(6,*) '*** NOT Reading updraft volumes ***'
      rewind (25)
      read(25,'(a)') dumstr
      DO it=1,maxt
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     charr(kz,it,1), charr(kz,it,2) 
        
        ENDDO
       itcount = itcount + 1
      END DO
      ENDIF

      ELSEIF ( ifile2 .eq. -5 ) THEN

      DO it=1,maxt
      
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     charr(kz,it,1), charr(kz,it,2), charr(kz,it,3), 
     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6) !, 
c     :     charr(kz,it,7), charr(kz,it,8) ! , charr(kz,it,9) 
c     :     ,charr(kz,it,10), 
c     :     charr(kz,it,11), charr(kz,it,12)
cc     :     drhop(kz,it), drhon(kz,it), drhomag(kz,it)
cc     :     ,cvdrhop(kz,it), cvdrhon(kz,it), cvdrhomag(kz,it),
cc     :     andrhop(kz,it), andrhon(kz,it), andrhomag(kz,it) 
        
        ENDDO
       itcount = itcount + 1
      END DO

      ELSEIF ( ifile2 .eq. -6 ) THEN

       DO it=1,maxt
      
          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     (charr(kz,it,i), i=1,ibase_micth)

c     :     charr(kz,it,1), charr(kz,it,2), charr(kz,it,3), 
c     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6), 
c     :     charr(kz,it,7), charr(kz,it,8), charr(kz,it,9) 
c     :     ,charr(kz,it,10), 
c     :     charr(kz,it,11), charr(kz,it,12)
cc     :     drhop(kz,it), drhon(kz,it), drhomag(kz,it)
cc     :     ,cvdrhop(kz,it), cvdrhon(kz,it), cvdrhomag(kz,it),
cc     :     andrhop(kz,it), andrhon(kz,it), andrhomag(kz,it) 
        
           
           charr(kz,it,25) = sum( charr(kz,it,12:24) )  ! total hydrometeor mass
           charr(kz,it,26) = sum( charr(kz,it,16:24) ) + charr(kz,it,14)  ! total precip mass
           charr(kz,it,27) =  charr(kz,it,13)  + charr(kz,it,15)  ! qc + qi
           charr(kz,it,28) =  charr(kz,it,22)  + charr(kz,it,24)  ! qg + qhl
!           charr(kz,it,ibase_micth+1) = sum( charr(kz,it,12:24) )  ! total hydrometeor mass
!           charr(kz,it,ibase_micth+2) = sum( charr(kz,it,16:24) ) + charr(kz,it,14)  ! total precip mass
!           charr(kz,it,ibase_micth+3) =  charr(kz,it,13)  + charr(kz,it,15)  ! qc + qi
           
           ENDDO
        itcount = itcount + 1
c        write(0,*) timearr(it),(charr(20,it,i), i=1,24)
       END DO

      ELSEIF ( ifile2 .eq. -7 ) THEN

       DO it=1,maxt
      
        DO kz=1,nz-1
         temg = Max( Min( sinit(kz,lt)*piinit%flt1d(kz), 313.15 ), 233.15 )
         temcg = Max( Min( sinit(kz,lt)*piinit%flt1d(kz), 273.15 ), 223.15 ) - 273.15
         felv(kz) = 2500837.367 * (273.15/temg)**((0.167)+(3.67e-4)*temg)
         felf(kz) = 333690.6098 + (2030.61425)*temcg - (10.46708312)*temcg**2
         fels(kz) = felv(kz) + felf(kz)
        ENDDO


          DO kz=1,nz-1
          read(25,*,end=600) 
     :     timearr(it), heights(kz), 
     :     (charr(kz,it,i), i=1,ibase)

c     :     charr(kz,it,1), charr(kz,it,2), charr(kz,it,3), 
c     :     charr(kz,it,4), charr(kz,it,5), charr(kz,it,6), 
c     :     charr(kz,it,7), charr(kz,it,8), charr(kz,it,9) 
c     :     ,charr(kz,it,10), 
c     :     charr(kz,it,11), charr(kz,it,12)
cc     :     drhop(kz,it), drhon(kz,it), drhomag(kz,it)
cc     :     ,cvdrhop(kz,it), cvdrhon(kz,it), cvdrhomag(kz,it),
cc     :     andrhop(kz,it), andrhon(kz,it), andrhomag(kz,it) 
        
           charr(kz,it,ibase+1) = charr(kz,it,16) + charr(kz,it,9) ! for total latent heating (but missing homogeneous qc freezing in old code)


           charr(kz,it,ibase+2) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (felf(kz)*(charr(kz,it,22) + charr(kz,it,19))
     &  +fels(kz)*(charr(kz,it,20)+charr(kz,it,21))  ! psub(mgs)    
     &  +felv(kz)*(charr(kz,it,18) + charr(kz,it,17)) ) !  pvap(mgs))

           charr(kz,it,ibase+3) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (felf(kz)*(charr(kz,it,22) + charr(kz,it,19)) ) 

           charr(kz,it,ibase+4) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (fels(kz)*(charr(kz,it,20)+charr(kz,it,21)) ) 

           charr(kz,it,ibase+5) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (felv(kz)*(charr(kz,it,18) + charr(kz,it,17)) ) !  pvap(mgs))

           charr(kz,it,ibase+6) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (felv(kz)*(charr(kz,it,18) ) ) !  evaporation only

           charr(kz,it,ibase+7) =
     &      1.0/(1004.0*piinit%flt1d(kz)*db(kz))*   
     &  (felf(kz)*(charr(kz,it,19))
     &  +fels(kz)*(charr(kz,it,21))  ! psub(mgs)    
     &  +felv(kz)*(charr(kz,it,18) ) ) !  pvap(mgs))

           ENDDO
        itcount = itcount + 1
c        write(0,*) timearr(it),(charr(20,it,i), i=1,24)
       END DO
      
      ENDIF
      
 600  CONTINUE 
 

       write(*,*) 'icg,icgn,icgp = ',icg,icgn,icgp
        
        Close(unit=25)
        
        
         GOTO 777
        ENDIF


      IF ( feflst ) THEN
        write(6,*) 'reading flst file ',flstname
        ifile = 15  
        open(unit=ifile,file=flstname,status='old',
     :       form='unformatted',iostat= istat )
        
c        read(ifile) flashm
c        read(ifile) flashp
c        read(ifile) flashn
c        read(ifile) drhom
        
        close(ifile)
      END IF
      

       
c
c Calculate potential and electric field.  If the 'nelec' file
c    contains space charge, then use that, otherwise add up from
c    'mica' file.
c
c
        imud = 0

c
c Calculate potential and electric field.  If the 'nelec' file
c    contains space charge, then use that, otherwise add up from
c    'mica' file.
c
c
      
      
      
c      dzmax = dz/gt(1,1,nz-1,imapz)
        dzmax = 1.0/gzt(nz-1,3)
      write(6,*) 'dzmax,dz = ',dzmax,dz
c      dslightz = dzmax
c      write(6,*) 'setting dslight = dzmax =',dslight

c         nxslm = NInt((nx-1)*dx/dslight) + 1
c         nyslm = NInt((ny-1)*dy/dslight) + 1
c
c IF vertical spacing is constant, then can have dslight < dz
c
         IF ( gzt(1,3) .eq. gzt(nz-1,3) ) THEN
c           nzslm = NInt((nz-1)*dz/dslight) + 1
           istretch = .false.
c           dslightz = dslight
         ELSE
c           dslightz = dzmax
c           nzslm = nz
           istretch = .true.
         ENDIF

           istretch = .true.
         
         write(6,*) 'istretch = ',istretch

      IF ( imud .eq. 1 ) THEN
         
         nb2s(1) = Max(NInt(4000./dx) , 8)
c         nb2s(2) = Max(NInt(4000./dslight) , 8)

         nxb2s(1) = nx + 2*nb2s(1) + 1
         nyb2s(1) = ny + 2*nb2s(1) + 1
         nzbt2s(1) = nz + Max(8,NInt(4000./dzmax + 0.5) ) ! nb2s(2)
         nzb2s(1) = nzbt2s(1) + 1

c         nxb2s(2) = nxslm + 2*nb2s(2) + 1
c         nyb2s(2) = nyslm + 2*nb2s(2) + 1
c         nzbt2s(2) = nzslm + Max(8,NInt(4000./dslightz + 0.5) ) ! nb2s(2)
c         nzb2s(2) = nzbt2s(2) + 1
         
         DO n=1,1
         IF ( nxb2s(n) .le. 41 ) THEN
           iixps(n) = 5
           iiexs(n) = 4
         ELSEIF ( nxb2s(n) .le. 49 ) THEN
           iixps(n) = 3
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 65 ) THEN
           iixps(n) = 2
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 81 ) THEN
           iixps(n) = 5
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 97 ) THEN
           iixps(n) = 3
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 113 ) THEN
           iixps(n) = 7
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 129 ) THEN
           iixps(n) = 2
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 161 ) THEN
           iixps(n) = 5
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 193 ) THEN
           iixps(n) = 3
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 257 ) THEN
           iixps(n) = 2
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 321 ) THEN
           iixps(n) = 5
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 385 ) THEN
           iixps(n) = 3
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 449 ) THEN
           iixps(n) = 7
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 513 ) THEN
           iixps(n) = 2
           iiexs(n) = 9
         ELSEIF ( nxb2s(n) .le. 641 ) THEN
           iixps(n) = 5
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 768 ) THEN
           iixps(n) = 3
           iiexs(n) = 9
         ELSE
           imudstop = 1
           write(6,*) 'WARNING! nxb2 too large for mudpack, will stop'
           STOP
         ENDIF
         
         IF ( nyb2s(n) .eq. nxb2 ) THEN
           jjyqs(n) = iixps(n)
           jjeys(n) = iiexs(n)
         ELSEIF ( nyb2s(n) .le. 41 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 4
         ELSEIF ( nyb2s(n) .le. 49 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 65 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 81 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 97 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 113 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 129 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 161 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 193 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 257 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 321 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 385 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 449 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 513 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 9
         ELSEIF ( nyb2s(n) .le. 641 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 768 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 9
         ELSE
           imudstop = 1
           write(6,*) 'WARNING! nyb2 too large for mudpack, will stop'
           STOP
         ENDIF
         
         IF ( nzb2s(n) .le. 41 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 4
         ELSEIF ( nzb2s(n) .le. 49 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 65 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 81 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 97 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 113 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 129 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 7
         ELSEIF ( nzb2s(n) .le. 161 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 193 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 7
         ELSE
           imudstop = 1
           write(6,*) 'WARNING! nzb2 too large for mudpack, will stop'
           STOP
         ENDIF
         
         nnxs(n) = iixps(n)*2**(iiexs(n)-1)+1
         nnys(n) = jjyqs(n)*2**(jjeys(n)-1)+1
         nnzs(n) = kkzrs(n)*2**(kkezs(n)-1)+1
         llworks(n) = 7*(nnxs(n)+2)*(nnys(n)+2)*(nnzs(n)+2)/2
         nz1d = nnzs(n)
         IF ( nz1d .gt. nzmax1 ) THEN
           write(6,*) 'Increase size of nzmax1 in driver to at least',
     :      nnzs(n)
           write(0,*) 'Increase size of nzmax1 in driver to at least',
     :      nnzs(n)
           STOP
         ENDIF
         nxb2s(n) = nnxs(n)
         nyb2s(n) = nnys(n)
         
         write(6,*) 'mudpack dimensions for grid ',n,' are ',
     :      nnxs(n),nnys(n),nnzs(n)

         ENDDO ! n
         
         nztop = nnzs(1)
c
c  set up extended transformation array
c

      DO k=1,nz
c        z1d2(k,3) = gt(1,1,k,imapz)    ! MFC
c        z1d2(k,4) = 0.5*(gt(1,1,k-1,imapz) + gt(1,1,k,imapz))  ! MFE
c        z1d2(k,1) = gt(1,1,k,mzdist)   ! scalar height
c        z1d2(k,2) = gt(1,1,k,mzdist)-0.5*dz/gt(1,1,k,imapz)  ! W-height
      ENDDO

c      dzmax = gt(1,1,nz,mzdist) - gt(1,1,nz-1,mzdist)
c      write(6,*) 'dzmax = ',dzmax
      
      IF ( nz .lt. nnzs(1) ) THEN
       DO k=nz,nnzs(1)
          z1d2(k,3) = z1d2(nz-1,3)
          z1d2(k,4) = z1d2(nz-1,4)
          z1d2(k,2) = z1d2(k-1,2) + dzmax
          z1d2(k,1) =  z1d2(k,2)+ 0.5*dzmax ! 0.5 * ( z1d2(k,2) + z1d2(k+1,2) ) ! z1d2(k-1,1) + dzmax
c          z1d2(k,3) = z1d2(nz,3)
c          z1d2(k,4) = z1d2(nz,4)
c          z1d2(k,1) = z1d2(k-1,1) + dzmax
c          z1d2(k,2) = z1d2(k-1,2) + dzmax
       ENDDO
       
       DO k=nz,nnzs(1)
        write(6,'(1x,i3,1x,f7.0,1x,f7.0,1x,f7.4,1x,f7.4,2(1x,f7.4))') 
     :    k,
     :   z1d2(k,1),z1d2(k,2),z1d2(k,3),z1d2(k,4),
     :   (z1d2(k,3)+z1d2(max(1,k-1),3))/2.,
     :   (z1d2(k,4)+z1d2(min(nz,k+1),4))/2.
       ENDDO 
      ENDIF
      
      nzs = nnzs(1)
      DO k = nz-2,nzs-1
                      z1d2(k,3) = dz / ( z1d2(k+1,2) - z1d2(k,2) )   ! MFC(k)
       IF( k .ne. 1 ) z1d2(k,4) = dz / ( z1d2(k,1)   - z1d2(k-1,1) ) ! MFE(k)

      ENDDO
      z1d2(nzs,3) = dz / (2.0*(z1d2(nzs,2) - z1d2(nzs-1,1)))
      z1d2(nzs,4) = dz / (2.0*(z1d2(nzs,2) - z1d2(nzs-1,1)))

      dlz = (z1d2(nnzs(1)-1,1) + z1d2(1,1))/(nnzs(1)-1)
      write(6,*) 'dlz = ',dlz
      
      
      write(6,*) 'recompute stretch factors for extended grid'
      DO k=1,nnzs(1)
        DO j=3,4
          z1d2(k,j) = z1d2(k,j)*dlz/dz
        ENDDO
        write(6,'(1x,i3,1x,f7.0,1x,f7.0,1x,f7.4,1x,f7.4,2(1x,f7.4))') 
     :    k,
     :   z1d2(k,1),z1d2(k,2),z1d2(k,3),z1d2(k,4),
     :   (z1d2(k,3)+z1d2(max(1,k-1),3))/2.,
     :   (z1d2(k,4)+z1d2(min(nz,k+1),4))/2.
      ENDDO

c 
c  Set up Fair weather field
c
      
      IF ( nzfair .lt. nz ) THEN
        write(6,*) 'STOP in sam.08.driver.x, nzfair must be > nz!',
     :    ' nzfair,nz=',nzfair,nz
        write(0,*) 'STOP in sam.08.driver.x, nzfair must be > nz!',
     :    ' nzfair,nz=',nzfair,nz
      STOP
      ENDIF
      
      IF ( imud .eq. 1 ) THEN
      do kz = 1,nztop
c      zgrd = .5*dz + (kz-1)*dz
      zgrd = z1d2(kz,1)
      ezfair(kz) = ezfairo *
     >           ( efb1*exp(-efa1*zgrd)
     >            +efb2*exp(-efa2*zgrd)
     >            +efb3*exp(-efa3*zgrd) )
      
       dezdz = ezfairo *
     >           ( -efa1*efb1*exp(-efa1*zgrd)
     >             -efa2*efb2*exp(-efa2*zgrd)
     >             -efa3*efb3*exp(-efa3*zgrd) )
       
       rhofair(kz) = dezdz*eperao

      zgrd = z1d2(kz,2)
      ezfairw(kz) = ezfairo *
     >           ( efb1*exp(-efa1*zgrd)
     >            +efb2*exp(-efa2*zgrd)
     >            +efb3*exp(-efa3*zgrd) )
       IF ( kz .le. nz) THEN
        write(6,*) 'kz,ezfairw,ezfairs: ',kz,ezfairw(kz),
     :   ezfair(kz)
cc     ,uz(kz,1),uz(kz,3), uz(kz,2),uz(kz,4)
       ENDIF
      end do
c
      do kz = 1,nztop-kd1
      if ( kz.eq.1 ) then
      potfair(kz) = -(0.5)*(ezfair(1)+ezfairs)*dz*(0.5)/
     :         z1d2(kz,3)
      end if
      if ( kz.gt.1 ) then
      potfair(kz) = potfair(kz-1)
     >             -(0.5)*(ezfair(kz)+ezfair(kz-1))*dz/
     :         z1d2(kz,3)
      end if 
      end do  
       
      ELSE
       
      do kz = 1,nztop
      zgrd = .5*dz + (kz-1)*dz
      ezfair(kz) = ezfairo *
     >           ( efb1*exp(-efa1*zgrd)
     >            +efb2*exp(-efa2*zgrd)
     >            +efb3*exp(-efa3*zgrd) )
      end do
c
      do kz = 1,nztop-kd1
      if ( kz.eq.1 ) then
      potfair(kz) = -(0.5)*(ezfair(1)+ezfairs)*dz*(0.5)
      end if
      if ( kz.gt.1 ) then
      potfair(kz) = potfair(kz-1)
     >             -(0.5)*(ezfair(kz)+ezfair(kz-1))*dz
      end if 
      end do  

      ENDIF ! imud .eq. 1
      
      eztop = ezfair(nztop)
         nbw = (nnxs(1) - nx - 1)/2
         nbe = nnxs(1) - nx - 1 - nbw
         nbs = (nnys(1) - ny - 1)/2
         nbn = nnys(1) - ny - 1 - nbs
      write(6,*) 'nztop,nbw,nbe,nbs,nbn,eztop = ',
     :          nztop,nbw,nbe,nbs,nbn,eztop


      
      ENDIF ! imud.eq.1
      
      write(6,*) 'here 1, neletyp=',neletyp

      chgneg = 0.0
      chgpos = 0.0
      cnapos = 0.0
      cnaneg = 0.0
      IF (  neletyp .ne. 2 .and. ifile2 .ge. 1 .and. lscb > 0 .and. lsceq > 0 ) THEN !{ {
      write(6,*) 'here 2'
      DO ia = lscb,lsceq
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
        rho(ix,jy,kz) = rho(ix,jy,kz) + st(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO

      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
       IF ( ion .eq. 1 ) THEN
       rho(ix,jy,kz) = rho(ix,jy,kz) 
     :     + ec*(st(ix,jy,kz,lscpi) - st(ix,jy,kz,lscni) )
       ELSEIF ( ion .eq. 2 ) THEN
       rho(ix,jy,kz) = rho(ix,jy,kz) 
     :     + ec*(st(ix,jy,kz,lscpi) - st(ix,jy,kz,lscni) )
     :     + ec*(st(ix,jy,kz,lscpli) - st(ix,jy,kz,lscnli) )
       ENDIF

      IF ( rho(ix,jy,kz) .gt. 0 ) THEN
        chgpos = chgpos + rho(ix,jy,kz)/gzt(kz,1)
      ELSE
        chgneg = chgneg + rho(ix,jy,kz)/gzt(kz,1)
      END IF
      
     
c      potw(ix,jy,kz) = -rho(ix,jy,kz) /eperao
          
        END DO
        END DO
        END DO

        
        write(6,*) 'chgpos,chgneg=',chgpos*dx*dy,chgneg*dx*dy
        
        IF ( chgpos .eq. 0.0 .and. chgneg .eq. 0.0 ) GOTO 111
        GOTO 111 
      
      ELSEIF (  neletyp .eq. 2 ) THEN !} {

c      write(0,*) 'here 3'

       DO kz=1,nz
       DO jy=1,ny
       DO ix=1,nx
c         t0(ix,jy,kz) =   rho(ix,jy,kz) 
         IF (kz .eq. 1) THEN
c           t0(ix,jy,kz) = t0(ix,jy,kz) - sccna(ix,jy)/eperao
c           rho(ix,jy,kz) = rho(ix,jy,kz) + sccna(ix,jy)
         END IF
     
c         potw(ix,jy,kz) = -t0(ix,jy,kz) /eperao
       END DO
       END DO
       END DO
      
      ELSE 
        write(6,*) 'WARNING: CANT CALCULATE POTENTIAL!'
        GOTO 111  
      END IF !}}
      
      write(6,*) 'before mud call, dz,dlz = ',dz,dlz

c regular grid:
         nbw = (nnxs(1) - nx - 1)/2
         nbe = nnxs(1) - nx - 1 - nbw
         nbs = (nnys(1) - ny - 1)/2
         nbn = nnys(1) - ny - 1 - nbs
         write(6,*) 'nbw,nbe,nbs,nbn=',nbw,nbe,nbs,nbn
c         nbz2 = nnzs(1) - nz
c         nzbt2 = nz + nbz2
c         nzb2 = nzbt2 + 1
c         nxb = nnxs(1)
c         nyb = nnys(1)
       
c
C      DO kz=1,4
C      CALL mud(nx,ny,nz,ng,
C     :   dx,dy,dz,
C     :   nnxs(1),nnys(1),nnzs(1),
C     :   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1),
C     :   nbw,nbe,nbs,nbn,
C     :   llworks(1),istretch,
C     :   rho,ibg,ibdcnd,jbdcnd,kbdcnd,kdbcdn,
C     :   neelec,elec,potfair,rhofair)
C      ENDDO
      
c
      write(6,*) 'after mud call'
C
C      call efield
C     >  (4,ng,nx,ny,nz,ngt
C     >  ,neelec
C     >  ,id1,jd1,kd1,istag,jstag,kstag
C     >  ,iex,iey,iez,iemag,ipot
C     >  ,imapx,imapy,imapz,ig13,ig23,1
C     >  ,dx,dy,dz,elec,gzt
C     >  ,ibg)    
C
c
c calculate Ez at w-points
c
c      DO kz=1,nz-1
c      DO jy=1,ny-1
c      DO ix=1,nx-1
c      
c       elec(ix,jy,kz,iezb) = 
c     > -(elec(ix,jy,kz,ipot)-elec(ix,jy,kz-1,ipot))/dz*
c     :   gt(ix,jy,kz,imapz)
c
c      ENDDO
c      ENDDO
c      ENDDO

c      END IF

 111  CONTINUE
 
      write(6,*) 'here 3a'
 
c      if ( nsmth .ne. 0 ) then
cc
c      ismthx = 1
c      ismthy = 1
c      ismthz = 1
cc
c      call smth
c     >  (nsmth,ismthx,ismthy,ismthz,1,
c     >   nx,ny,nz,na,ng,nx,ny,nz,nx,ny,nz,
c     >   1,1,1,01,
c     >   nx-1,ny-1,nz-1,na,
c     >   t0,a,ab)
cc
c      call smth
c     >  (nsmth,ismthx,ismthy,ismthz,1,
c     >   nx,ny,nz,03,ng,nx,ny,nz,nx,ny,nz,
c     >   1,1,1,01,
c     >   nx,ny-1,nz-1,01,
c     >   t0,v,vb)
cc
c      call smth
c     >  (nsmth,ismthx,ismthy,ismthz,1,
c     >   nx,ny,nz,03,ng,nx,ny,nz,nx,ny,nz,
c     >   1,1,1,02,
c     >   nx-1,ny,nz-1,02,
c     >   t0,v,vb)
cc
c      call smth
c     >  (nsmth,ismthx,ismthy,ismthz,1,
c     >   nx,ny,nz,03,ng,nx,ny,nz,nx,ny,nz,
c     >   1,1,1,03,
c     >   nx-1,ny-1,nz,03,
c     >   t0,v,vb)
cc
c      call smth
c     >  (nsmth,ismthx,ismthy,ismthz,0,
c     >   nx,ny,nz,01,ng,nx,ny,nz,nx,ny,nz,
c     >   1,1,1,01,
c     >   nx-1,ny-1,nz-1,01,
c     >   t0,p,pb)
cc
c      end if
      
      write(6,*) 'here 3b'
c                
c  comparison arrays---from differences x2-x
c
c      if ( icomp .eq. 1 ) then
c      print*,'icomp',icomp
c      open(1,file=filnm2,form='unformatted')
c      read(1) a2,d2,p2,v2
c      close(unit=1)
c      do 96301 kz = 1,nz
c      do 96301 jy = 1,ny
c      do 96301 ix = 1,nx
c      a2(ix,jy,kz,01) = a2(ix,jy,kz,01)-st(ix,jy,kz,01)
c      a2(ix,jy,kz,02) = a2(ix,jy,kz,02)-st(ix,jy,kz,02)
c      a2(ix,jy,kz,03) = a2(ix,jy,kz,03)-st(ix,jy,kz,03)
c      a2(ix,jy,kz,04) = a2(ix,jy,kz,04)-st(ix,jy,kz,04)
c      a2(ix,jy,kz,05) = a2(ix,jy,kz,05)-st(ix,jy,kz,05)
c      if ( na .gt. 5 ) then
c      a2(ix,jy,kz,06) = a2(ix,jy,kz,06)-st(ix,jy,kz,06)
c      a2(ix,jy,kz,07) = a2(ix,jy,kz,07)-st(ix,jy,kz,07)
c      a2(ix,jy,kz,08) = a2(ix,jy,kz,08)-st(ix,jy,kz,08)
c      end if
c      v2(ix,jy,kz,01) = v2(ix,jy,kz,01)-v(ix,jy,kz,01)
c      v2(ix,jy,kz,02) = v2(ix,jy,kz,02)-v(ix,jy,kz,02)
c      v2(ix,jy,kz,03) = v2(ix,jy,kz,03)-v(ix,jy,kz,03)
c      p2(ix,jy,kz) = p2(ix,jy,kz)-p(ix,jy,kz)
c96301 continue
c      end if
c
c  create zb
c
      do kz = 1,nz
      zb(kz) = z1d2(kz,1) ! (kz-1)*dz + dzadd
      end do 
c
c  prepare derived and reorganized data
c
c
c.. theta-e and RH calculation
c
      print*,'THETA-E START'
      do kz = 1,nz-1
      do jy = 1,ny-1
      do ix = 1,nx-1
c      tc(ix,jy,kz) = st(ix,jy,kz,lt)
      pres = pn(ix,jy,kz)+pb(kz)
      if ( ithetae .ge. 1 ) then
        qwv = max(st(ix,jy,kz,lv) , 0.0 )
        theta1 = st(ix,jy,kz,lt)
        the(ix,jy,kz) = thetae(qwv,theta1,pres)
      endif
      tmr = st(ix,jy,kz,lt)*(pres/100000.0)**(287./1004.)
      qsat = (380./pres)*exp(17.2693882*(tmr-273.16)/(tmr-35.86))
      rh(ix,jy,kz) = 100.*(st(ix,jy,kz,lv))/qsat
      end do
      end do
      end do
      print*,'THETA-E DONE'
c
c
c
777   continue
c
c
c
      write(*,*) 'DONE WITH INITIALIZATION: START PLOTING'
c
c
c
      write(*,*) 'START WITH PLOTTING'
c
c
c
      ioverlay = 0
      iclr = 0
      ICNT = 0
 2210 continue
      write(*,*) 'IOVERLAY = ',ioverlay
      if ( ioverlay .gt. 0 ) icnt = 0
      if ( ioverlay .lt. 0 ) write(*,*) '*** OVERLAY PLOT ***'
      if ( ioverlay .lt. 0 ) then
      ioverpnt = 0
      else
      ioverpnt = 1
      end if
c
c
      write(*,*) 'idebg=1'
      idone = 0
      IFORCEOL = 0
c
      if ( ioverlay .eq. -17 .or. ioverlay .eq. -2
     >  .or. ioverlay .eq. -19 ) then
      iforceol = 1
      end if
c
      if ( ioverlay .eq. -17 .or. ioverlay .eq. -23 
     > .or. ioverlay .eq. -32 .or. ioverlay .eq. -40 ) then
      iforceol = 1
      end if
c
      if ( iforceol .eq. 1 ) then
c
      if ( ioverlay .eq. -2 ) then
      write(*,*) 'DEFAULT TEMPERATURE AND CLOUD  OVERLAY'
      rclon = -40.
      rclop = -40.
      rchin = 00.
      rchip = 00.
      rcinn = -10.
      rcinp = -10.
      idx = 2
      xmult = 1.0
      idx = -ioverlay
      idone = 1
      end if
c
      if ( ioverlay .eq. -17 ) then
      write(*,*) 'DEFAULT CLOUD OVERLAY'
      ifread = 0
      rclon = 0.0001
      rclop = 0.0001
      rchin = 0.0001
      rchip = 0.0001
      rcinn = 0.0000001
      rcinp = 0.0000001
      idx = -ioverlay
      xmult = 1.0
      idone = 1
      end if
c
      if ( ioverlay .eq. -23 .or. ioverlay .eq. -32) then
      write(*,*) 'DEFAULT VECTOR OVERLAY'
      ifread = 0
c      ivint = ifvint
c      hiv = fhiv
c     lengthv = iflengthv
      ivint = 2
      hiv = 60.
      lengthv = 60
c      idx = -ioverlay
      idone = 1
      xmult = 1.0
      end if
c
      if ( ioverlay .eq. -19 ) then
      write(*,*) 'DEFAULT TEMPERATURE AND CLOUD  OVERLAY'
      rclon = -40.
      rclop = -40.
      rchin = 00.
      rchip = 00.
      rcinn = -10.
      rcinp = -10.
      idx = 2
      ioverlay = -17
      xmult = 1.0
      end if
c
      if ( ioverlay .eq. -40 ) then
      write(*,*) 'DEFAULT VECTOR AND CLOUD  OVERLAY'
c     ivint = ifvint
c     hiv = fiv
c     lengthv = iflengthv
      ivint = 3
      hiv = 30.
      lengthv = 60
      idx = 23
      ioverlay = -17
      xmult = 1.0
      end if
c
      end if
c
      write(*,*) 'idebg=2'
c
      icnt = icnt + 1
c
      if ( icnt .eq. 1 ) then
      nplot = nplot + 1
      end if
c
      write(*,*) 'idebg=3'
c
c
c  INPUT PLOTTING INFORMATION
c
c
      write(*,*) 'INPUT: ido, ipwr, ipltyp, iover, idx, ib, ie, jb, je'
      write(*,*) '(cont), kb, ke, xmult,rlo, rhi rin,tiso,sf,is,sl,mcl'
      write(*,*) '(cont), cll'
      write(*,*) 'NOTE: rlo=xeye=alwv, rhi=yeye=ahhv, rin=zeye=length'
      write(*,*) 'NOTE: tiso=ivint'
c      read(lunr,*) ido, ipware, ipltyp, ioverlay,
c     >  idx, ib, ie, jb, je, kb, ke,
c     >  xmult, rlo, rhi, rin, tiso, sf, ishade,icrfil,mcl,cll
      
 10    Continue
      read(lunr,'(a)') line
      write(6,*) 'line= ',line
      IF (line(1:1) .eq. ';' .or. line(1:1) .eq. ' ' .or.
     :    line(1:1) .eq. '#' ) GOTO 10

      IF (line(1:1) .eq. '&' ) THEN
       read(lunr,NML=plot)
       GOTO 10
      ENDIF
      
      read(line,*) ido
      IF ( ido .eq. 2 ) THEN
      read(line,*) ido, ipware, ipltyp, xib0, xie, yjb, yje, zkb, zke
      write(6,*) 'ido, ipware, ipltyp, xib0, xie, yjb, yje, zkb, zke',
     : ido, ipware, ipltyp, xib0, xie, yjb, yje, zkb, zke
        IF ( ifile2 .eq. -1 .and. 
     :    .not. (ipltyp .eq. 10 .or. ipltyp .eq. 11) ) THEN
          write(6,*) 'For ifile2 = -1 Must have ipltyp=10 or 11!'
          write(6,*) 'Stopping now!  ipltyp= ', ipltyp
          STOP
        ELSEIF ( ifile2 .eq. -2 .and. 
     :    .not. (ipltyp .eq. 10 ) ) THEN
          write(6,*) 'For ifile2 = -2 Must have ipltyp=10!'
          write(6,*) 'Stopping now!  ipltyp= ', ipltyp
          STOP
        ENDIF
        
        IF ( do_dualpol .and. .not. dualpol_fulldomain ) THEN
        
         ibgn = Max(1,  Nint(xib0)-1)
         iend = Min(nx-1, Nint(xie)+1)
         jbgn = Max(1,  Nint(yjb)-1)
         jend = Min(ny-1, Nint(yje)+1)
         IF ( Nint(zkb) == Nint(zke) ) THEN
         ! XY slab
           kbgn(:,:) = Max(1, Nint(zkb))
           kend(:,:) = Min(nz-1, Nint(zke))
         ELSE
           ! z is km, so just do the whole slice
          kbgn(:,:) = 1
          kend(:,:) = nz-1
         ENDIF
         write(0,*) 'ibgn,etc = ',ibgn,iend,jbgn,jend

!     qsw(:,:,:) = 0.0
!     qgw(:,:,:) = 0.0
!     qhw(:,:,:) = 0.0 
!     
!     sumZh(:,:,:) = 0.0
!     sumZv(:,:,:) = 0.0
!     logZdr(:,:,:) = 0.0
!     logZdrrain(:,:,:) = 0.0
!     sumZhv(:,:,:) = 0.0
!     Kdp(:,:,:) = 0.0
!     Ahh(:,:,:) = 0.0
!     Avv(:,:,:) = 0.0
!     rhv(:,:,:) = 0.0
!     rhoms(:,:,:) = rho_qs
      
!         logZ(ibgn:iend,jbgn:jend,1:nz) = dbz%flt3d(ibgn:iend,jbgn:jend,1:nz)
!         logZ(1,1,nz) = 1.e11 ! set to use dbz as mask -- DEPRECATED for dpmask

         tair(1:nx,1:ny,1:nz) = tt0(1:nx,1:ny,1:nz)
         rhoair(1:nx,1:ny,1:nz) = dn(1:nx,1:ny,1:nz)

       IF ( microp(1:1) == 'Z' ) THEN
       call refl_rsa_array(addtorain,MPflg,MFflg,nx,ny,nz,dirscatt,wavelen,ibgn,iend,jbgn,jend,kbgn,kend,  
     :      rhoair,qr3d,qs3d,qg3d,qh3d,Ntr,Nts,Ntg,Nth,alphar3d,alphas3d,alphag3d,alphah3d,
     :      qsw,qgw,qhw,rhograup,rhohail, 
     :      logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh,tair, ! note that rhoms,rhomg,rhomh are overwritten (for wet particles)
     :       dpmask, dbz%flt3d(1:nx,1:ny,1:nz) )
        ELSEIF (microp(1:3) == 'TAK' ) THEN
        
           call refl_rsa_array_takahashi(MFflg,nx,ny,nz,ng,ntakrd,ntakit,ntakid,ntakpd,ntakpd,
     :      ibgn,iend,jbgn,jend,kbgn,kend,dirscatt,wavelen,rhoair,qr3d,qs3d,qg3d,qh3d,  
     :      st(-ng+1,-ng+1,-ng+1,lnr),st(-ng+1,-ng+1,-ng+1,lni),st(-ng+1,-ng+1,-ng+1,lnh),
     :      st(-ng+1,-ng+1,-ng+1,lnhl),
     :      qsw,qgw,qhw,rhograup,rhohail, 
     :      logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh, dpmask)        
        
        ENDIF
        
        
!       write(0,*) 'k,dbz,logz,logzdr,rhv,alpr,alph,alphl,qgw,qhw,rhograup,rhomg'
!       i = (ibgn + iend)/2
!       j = (jbgn + jend)/2
!       DO k = kbgn(i,j),kend(i,j)
!         write(0,*) k,dbz%flt3d(i,j,k)-logz(i,j,k),logz(i,j,k),logzdr(i,j,k),rhv(i,j,k),alphar3d(i,j,k),alphag3d(i,j,k),
!    :           alphah3d(i,j,k),qg3d(i,j,k),st(i,j,k,lh),qgw(i,j,k),qhw(i,j,k),qr3d(i,j,k),st(i,j,k,lr),
!    :           rhograup(i,j,k),rhomg(i,j,k)
!       ENDDO
        
        ENDIF
        
       GOTO 10
      ELSEIF ( ido .lt. 0 ) THEN
        write(6,*) 'Stop signal encountered.  Buh-bye!'
        IF ( numtime .eq. 1 .or. ntime .eq. numtime) THEN
          GOTO 9010
        ELSE
         ntime = ntime + 1
         time = time0 + dtdump*ntime
         RETURN ! GOTO 909
        ENDIF
        
      ELSEIF ( ido .eq. 3 ) THEN ! change to a new time
        IF ( numtime .gt. 1 ) THEN
         write(0,*) 'Cannot use ido=3 with time cycle!!!'
         GOTO 9010
        ELSE
         read(line,*) ido,time
         GOTO 9091
        ENDIF
      ELSEIF ( ido .eq. 4 ) THEN ! change line scaling
         read(line,*) ido,scale
         call ngsetr('NO',scale)
         GOTO 10
      ELSEIF ( ido .eq. 5 ) THEN ! vector label controls
c  read in your control values here.
c         read(line,*) ido
c         
         GOTO 10
c      ELSEIF ( ido .eq. 98 .and. ioffset .eq. 0 ) THEN ! wdt output
c          ioffset = 1
c          read(line,*) ido,xoffset,yoffset
c          xoffset = xoffset + umv*dtdump
c          yoffset = yoffset + vmv*dtdump
c      ELSEIF ( ido .eq. 99 ) THEN ! wdt output
c         IF ( iwdt .eq. 0 ) THEN  ! open binary and ascii file
c            iwdt = 1
c            read(line,*) ido,umv,vmv,fwdt,xoffset,yoffset   ! base file name for wdt output
c            idum = len(fwdt)
c            CALL strlnth(fwdt,idum)  ! now idum is actual length of fwdt
c            
c            IF ( ioffset .eq. 0 .and. 
c     :        (xoffset.ne.0.0 .or. yoffset.ne.0.0) ) THEN
c                 ioffset = 1
c                 xoffset = xoffset + umv*dtdump
c                 yoffset = yoffset + vmv*dtdump
c            ENDIF
c            
c            iwdtas = 98
c            iwdtbin = 99
c            
c            OPEN(unit=iwdtas,form='formatted',status='unknown',
c     :           file= fwdt(1:idum)//'.txt' )
cc            OPEN(unit=iwdtbin,form='unformatted',status='unknown',
cc     :           file= fwdt(1:idum)//'.bin' )
c          write(iwdtas,502) 
c     :       nx-1,ny-1,dx,dy
c 502  Format(2(i3,','), f7.2,',',f7.2)
c
c         ENDIF
c          
c          
cc          write(iwdtas,'(f8.2,a,1x,2f10.2)') 
c          write(iwdtas,501) 
c     :       time,smin,xoffset,yoffset
c          
c 501  Format(f8.2,',',a,',',1x,f10.2,',',f10.2)
c          xoffset = xoffset + umv*dtdump
c          yoffset = yoffset + vmv*dtdump
c          
c          DO jy=1,ny-1
cc            ix = 50
c            DO ix=1,nx-1
c               tem1(ix,jy,1) = elec(iez)%flt3d(ix,jy,1)
c               tem1(ix,jy,2) = dbz(ix,jy,1)
c               tem1(ix,jy,3) = flashn(ix,jy,1) + flashp(ix,jy,1)
c               write(iwdtas,'(2(i3,","),1x,f12.3,",",f8.3,",",f4.0)') 
cc               write(iwdtas,503), ix,jy,
c     :           ix,jy, tem1(ix,jy,1),tem1(ix,jy,2),tem1(ix,jy,3)
c            ENDDO
c          ENDDO
c
c 503  Format(2(i3,','),1x,f12.3,',',f8.3,',',f4.0)
c
cc          Call writewdt(nx-1,ny-1,iwdtas,tem1,
c         
c         GOTO 10
      ELSEIF ( ido .eq. 0 ) THEN
        GOTO 10
      ELSE
        IF ( epsflag .eq. 1 ) THEN
          epsflag = 0
          
        CALL GOPKS(6,0)
        write(outnameps,'(a,a,i2.2,a)') 
     :             outname(1:inameps),'.',ieps,'.eps'
        CALL ngsetc('ME',outnameps)
c        write(0,*) ieps,'outnameps: ',outnameps
        CALL GOPWK(1,9,iotyp)  ! eps output
        CALL GACWK(1)
        call agseti('FRAME.',2)  !  set no automatic frame advance; allows overlays
         call ngseti('WO',1)
         call ngsetr('NO',2.)
         call ngseti('PA',10000)


c
c        set some parameters
c
         do i = 1,129
         ifill(i) = i
         end do

         CALL setclr
         
         CALL GSTXFP (-29,0)  ! courier font, string precision
 
        ENDIF
      read(line,*) ido, ioverlay,
     >  idx, 
     >  xmult, rlo, rhi, rin, tiso, sf, ishade,icrfil,mcl,cll,ilab

      write(6,*) 'idx= ',idx,' = ',title(Abs(idx))

          ioverpnt2 = 0
        IF ( ioverlay .gt. 0 .and. ipltyp .ne. 11 ) THEN 
          IF ( ioverlay .ge. 1 ) THEN
          ioverpnt2 = 1
c          ELSEIF ( ioverlay .eq. 2 ) THEN
c           ioverpnt2 = 1
          ENDIF
c        ELSE
c          ioverpnt = 0
        ENDIF
      ENDIF
      xib = xib0
      
      ib = Int(xib)
      ie = Int(xie)
      jb = Int(yjb)
      je = Int(yje)
      kb = Int(zkb)
      ke = Int(zke)
      kb1 = 0
      
      irev = 0
      IF ( jb .gt. je .and. ipltyp .ne. 5 ) THEN
        irev = 1
        i = jb
        jb = je
        je = i
      ENDIF
      
       d1pobl = 0.0
      IF ( ipltyp .eq. 5 .and. xib0 .lt. 0.0 ) THEN
       ib = 0
       d1pobl = -xib0*1000.0
       IF ( yje .le. 0.0  .or. 
     :   ( yje.ge.0.0 .and. ( yje .le. Abs(xib0) ) ) ) THEN
          IF ( yje .ge. 0.0 ) THEN
             write(0,*) 'WARNING: resetting yje from ',yje,' to ',
     :        Abs(yje) + Abs(xib0)
             write(0,*) 'Maybe you shoud make yje negative?'
          ENDIF
         yje = Abs(yje) + Abs(xib0)
       ENDIF
      ELSE
       d1pobl = 0.0
      ENDIF
      
      nzset = 0
      n2set = 0
      hgtset = 0.0
      kstrzb = 0
      IF ( kb .le. 0 ) THEN
       ke = nz-1
       write(*,*) 'ke, nz = ',ke, nz
       IF ( zke*1000. .gt. z1d2(nz-1,1) ) THEN
        write(0,*) 'resetting zke down to heighest allowed'
        zke = Int(z1d2(nz-1,1)/1000)
c        nzset = nz-1
        write(0,*) 'zke = ',zke
       ENDIF

c
c find lower limit
c
       hgtlow = 0.0 ! -1.0e9
       IF ( kb .eq. 0 ) THEN
         !write(*,*) 'resetting zkb to zero'
         kb = 1
         kb1 = 0
         kstrzb = 0
         hgtlow = 0.0
       ELSE
       DO kz=1,nz-1
!         write(0,*) 'kb loop: kz,z1d2,zkb,-zkb*1000' , kz,z1d2(kz,1),zkb,-zkb*1000
         IF ( z1d2(kz,1) .le. -zkb*1000. ) THEN
           kb1 = kz
           kb = kz
           kstrzb = kz - 1
           hgtlow = -zkb
!           write(0,*) 'kb loop IF : kz,z1d2,zkb,-zkb*1000' , kz,z1d2,zkb,-zkb*1000
         ENDIF
       ENDDO
!         write(0,*) 'resetting kb to ',kb,' kb1,kstrzb,hgtlow = ',kb1,kstrzb,hgtlow
       
       ENDIF

c find upper limit
       DO kz=1,nz-1
         IF ( z1d2(kz,1) .le. zke*1000. ) THEN
           ke = Min(nz-1,kz+2)   ! 11/7/2016 erm added "+2" to load data one more level so that contours fill the plot
           hgtset = zke
           n2set = Int((hgtset-hgtlow)*1000./dzin + 0.01) + 1
           nzset = kz-kb1+2   ! 11/7/2016 erm added "+2" to load data one more level so that contours fill the plot
           write(*,*) 'kz,zke,z1d2(kz,1),ke = ',kz,zke,z1d2(kz,1),ke
         ENDIF
       ENDDO
       dz = dzin
      ELSE
       dz = dzsave
c       hgtset = 0.001*(gt(1,1,n2,mzdist) + 
c     :   0.5*(gt(1,1,n2+1,mzdist)-gt(1,1,n2,mzdist)))
      ENDIF
 
      do kz = 1,nz
        z = (kz - 1)*dz
        DO k = 1,nz
         IF ( z1d2(k,1) .le. z .and. z1d2(k+1,1) .gt. z ) THEN
          tmp = (z-z1d2(k,1))*
     :      (temkz(k+1) - temkz(k))/(z1d2(k+1,1)-z1d2(k,1))
          
          templt(kz) = temkz(k) + tmp
          ENDIF
        ENDDO
        
      end do
      
      
!      xfacb = xfacb0
!      xdivb = xdivb0
      
      IF ( ipltyp .eq. 10 ) THEN
c        ie = Min(ie,itcount)
      ENDIF

      

c compute slice length 
      IF (ipltyp .eq. 5 ) THEN
        xib = Max( 1.0, xib )
        yjb = Max( 1.0, yjb )
        dd0 = yje*1000
        je = 0
        write(6,*) 'dd0 =',dd0
        phi = xie
        ie = 0
c        kb = 1

c check for a plot at some extreme value
      if (  ie .lt. -001  ) then
c      print*,'EXTREMES',ib,jb,kb
c
c  plots at extremes of w+
c
      if ( (ie) .eq. -99 ) then
      call extrms(1,ng,1,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,w%flt3d)
      end if
c
c  plots at extremes of qh
c
      if ( (ie) .eq. -98 ) then
      call extrms(lh,ng,na,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,st)
      end if
c
c  plots at extremes of qr
c
      if ( (ie) .eq. -97  ) then
      call extrms(lr,ng,na,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,st)
      end if
c
c now compute location of 'left' edge of the slice       
c
        IF ( phi .gt. 0.0 .and. phi .lt. 90.0 ) THEN
           dd1 = Float(imx - 1)/Cos(phi)
           dd2 = Float(jmx - 1)/Sin(phi)
           dd = Min( dd1, dd2 )
           IF ( dd0 .gt. 0.0 ) THEN
            dd = Min(dd,dd0)
           ENDIF
           write(6,*) 'dd=',dd
           xib = Float(imx) - dd*Cos(phi)
           yjb = Float(imx) - dd*Sin(phi)
        ELSEIF ( phi .lt. 0.0 .and. phi .gt. -90 ) THEN
           dd1 = Float(imx - 1)/Cos(phi)
           dd2 = Float(ny-1-jmx)/Abs(Sin(phi))
           dd = Min( dd1, dd2 )
           IF ( dd0 .gt. 0.0 ) THEN
            dd = Min(dd,dd0)
           ENDIF
           xib = Float(imx) - dd*Cos(phi)
           yjb = Float(imx) + dd*Abs( Sin(phi) )
        END IF
        
        xib = Max ( xib, 1.0 )
        yjb = Min ( Float(ny-1), Max(1.0,yjb) )
            
      END IF ! ( ie .lt. -1 )
c
c now find length of slice and how many horizontal points are needed      
c              
        IF ( phi .gt. 0.0 .and. phi .lt. 90.0 ) THEN
          dd1 = dy*(ny - yjb)/Sin(phi*pi/180.0)
          dd2 = dx*(nx - xib)/Cos(phi*pi/180.0)
          dd = Min(dd1,dd2)
           IF ( dd0 .gt. 0.0 ) THEN
            dd = Min(dd,dd0)
           ENDIF
          write(6,*) 'dd,dd1,dd2,dd0',dd,dd1,dd2,dd0
          nl = Int( (dd-d1pobl)/dx +1 )
        ELSEIF ( phi .lt. 0.0 .and. phi .gt. -90 ) THEN
          dd1 = dy*(1.0 - yjb)/Sin(phi*pi/180.0)
          dd2 = dx*(nx - xib)/Cos(phi*pi/180.0)
          dd = Min(dd1,dd2)
           IF ( dd0 .gt. 0.0 ) THEN
            dd = Min(dd,dd0)
           ENDIF
          write(6,*) 'dd,dd1,dd2,dd0',dd,dd1,dd2,dd0
          nl = Int( (dd-d1pobl)/dx +1)
        
        ELSEIF ( phi .eq. 0.0 ) THEN  ! xz slice
          ipltyp = 2
          IF ( dd .eq. 0.0 ) THEN
            ib = Int(xib + 0.001)
            ie = nx-1
          ELSE
            ib = Int(xib + 0.001)
            ie = Min(nx-1, ib + Int(dd/dx))
          END IF
          jb = Int(yjb + 0.001)
          je = jb
c          kb = 0
c          ke = 0
        ELSEIF ( Abs(phi) .eq. 90.0 ) THEN  ! yz slice
          ipltyp = 3
          ib = Int(xib + 0.001)
          ie = ib
          jb = 0
          je = 0
c          kb = 0
c          ke = 0
        ELSE
          write(6,*) 'value for phi is out of range!'
          STOP
        END IF
      END IF
      
      ivint = tiso+1.e-20
      flov = rlo
      hiv = rhi
      lengthv = rin
c
      write(*,*) 'idebg=4'
c
      ifll = 0
      nclr=15
      icrtyp = 0
      ilabor = 0
      idolab = 0
        ishadcnt0 = 0
      IF ( ishade == 1 ) ishade = 59 ! set to black and avoid SHDH1
      
      if ( ishade+icrfil .ne. 0 ) then
c      read(lunr,*) chsval,chsang,icrtyp,nclr,rfac1,rfac2,rfac3,rfac4
c     >  ,idolab,ilabor,ilablab,ishft,ishadcnt0
        chsval = 0.1
        chsang = 0
        IF ( icrfil .le. 1 ) THEN
        icrtyp = 0
        nclr = 15
        ENDIF
!        rfac1 = 0
!        rfac2 = 0
!        rfac3 = 0
!        rfac4 = 0
        IF ( (idx .eq. 22 .or. idx .eq. 35 .or. idx .eq. 991 .or. idx == 601) 
     :          .and. icrfil .eq. 1 ) THEN
          ilabor = 0  ! 0 for horiz (bottom), 1 for vertical (right)
          idolab = 1
       ELSEIF ( icrfil .eq. 2 .or. icrfil .eq. 3 .or. ipware .eq. 2 ) THEN
          ilabor = 1
          idolab = 1
       ELSE
          idolab = 0
          ilabor = 0
        ENDIF
        
        IF ( iplotcolorbar == 0 ) idolab = 0
        
        IF ( ilaborhv .eq. 0 .or. ilaborhv .eq. 1 ) ilabor = ilaborhv
        
        ilablab = 0
        ishft = 0
        ishadcnt0 = 0
        IF ( Abs( ishade ) .ge. 130 ) THEN
         ifll = 1
        ENDIF
        IF ( Abs( ishade ) .gt. 337 ) THEN
         write(6,*) 'ishade is OUT OF RANGE! = ',ishade
         STOP
        ENDIF

      write(*,*) 'chsval etc.', chsval,chsang,icrtyp,nclr,rfac1,rfac2,rfac3,rfac4
     >  ,idolab,ilabor,ilablab,ishft,ishadcnt0
      
       ELSEIF ( ipware .eq. 2 ) THEN
          ilabor = 0
          idolab = 1
          ilablab = 1
      end if
c
      write(*,*) 'idebg=5'
c
      WRITE(*,*) 'ido, ipwr, iptyp, iov, idx, ib, ie, jb, je, kb, ke,'
      WRITE(*,*) 'xmult, rlo, rhi, rin, tiso, sf, ishad,icolf,mcl,cll'
      write(*,*) ido, ipware, ipltyp, ioverlay, 
     > idx, ib, ie, jb, je, kb, ke,
     > xmult, rlo, rhi, rin, tiso, sf, ishade,icrfil,mcl,cll

      write(6,*) 'idx= ',idx,' = ',title(Abs(idx))
      write(*,*) 'ishade,icrfil',ishade,icrfil
c
      write(*,*) 'idebg=6'
c
CH  read in DF
      df = 1.
c
      write(lprn2,9111) nplot
9111  format(i6.5)
      IF ( idx .gt. 0 ) THEN
      if ( ipltyp .eq. 1 ) then
      write(lprn1,9111) kb
      write(lunw,*) 'XY PLT:',title(Abs(idx)),':K=',lprn1,': #= ',lprn2
      write(23,*) 'XY PLT:',title(Abs(idx)),' :K=',lprn1,': #= ',lprn2
      end if
      if ( ipltyp .eq. 2 ) then
      write(lprn1,9111) jb
      write(lunw,*) 'XZ PLT:',title(Abs(idx)),' :J=',lprn1,': #= ',lprn2
      write(23,*) 'XZ PLT:',title(Abs(idx)),' :J=',lprn1,': #= ',lprn2
      end if
      if ( ipltyp .eq. 3 ) then
      write(lprn1,9111) ib
      write(lunw,*) 'YZ PLT:',title(Abs(idx)),' :I=',lprn1,': #= ',lprn2
      write(23,*) 'YZ PLT:',title(Abs(idx)),' :I=',lprn1,': #= ',lprn2
      end if
      if ( ipltyp .eq. 4 ) then
      write(lunw,*) 'XYZ PLT OF ',title(Abs(idx)),': # ',nplot
      write(23,*) 'XYZ PLT OF ',title(Abs(idx)),': # ',nplot
      end if
      if ( ipltyp .eq. 5 ) then
      write(lprn1,9111) jb
      write(lunw,*) 'OBLIQ PLT:',title(Abs(idx)),' :J=',lprn1,': #= ',lprn2
      write(23,*) 'OBLIQ PLT:',title(Abs(idx)),' :J=',lprn1,': #= ',lprn2
      end if
      ENDIF ! idx.gt.0
c
      write(*,*) 'idebg=7'
c
      if ( ipltyp .lt. 1 .or. ido .lt. 0 ) go to 9010
c
      if ( ido .eq. 0 ) go to  2210
c
      write(*,*) 'idebg=8'
c
      sdx = 1.
      sdz = sf 
      xfac = xfacb
      yfac = yfacb
      zfac = zfacb
      xdiv = xdivb
      ydiv = ydivb
      zdiv = zdivb

      if ( ipltyp .eq. 1 ) then
      sdx = 1.
      sdy = sf 
      xfac = xfacb
      yfac = yfacb
      zfac = zfacb
      xdiv = xdivb
      ydiv = ydivb
      zdiv = zdivb
      end if
      if ( ipltyp .eq. 2 .or. ipltyp .eq. 5) then
      sdx = 1.
      sdz = sf 
      xfac = xfacb
      yfac = yfacb
      zfac = zfacb
      xdiv = xdivb
      ydiv = ydivb
      zdiv = zdivb
      end if
      if ( ipltyp .eq. 3 ) then
      sdy = 1.
      sdz = sf 
      xfac = xfacb
      yfac = yfacb
      zfac = zfacb
      xdiv = xdivb
      ydiv = ydivb
      zdiv = zdivb
      end if
c
      if ( isforz .eq. 1 .and. sf .eq. 1.0 ) then
      sdz = sforz
      end if
c
      if ( ipltyp .eq. 4 ) then
      eye(1) = rlo
      eye(2) = rhi
      eye(3) = rin
c      isoflg = sf+0.001
      isoflg = 7
      end if
c
      if ( ipltyp .eq. 23 .or. ipltyp .eq. 32 ) then
      ivint = tiso+0.001 
      hiv = rhi
      lengthv = rin+0.001
      end if
c
      rclop = rlo
      rchip = rhi
      rcinp = rin
      rclon = rlo
      rchin = rhi
      rcinn = rin
c
c     end if
c
      write(*,*) 'idebg=9'
c
      if ( (ib .eq. 0 .or. ie .eq. 0) .and. ifile2.ge.0 ) then
      ib = 1
      ie = nxstag
      nxp = nxstag
      end if
c
      if ( jb .eq. 0 .or. je .eq. 0 ) then
      jb = 1
      je = nystag
      nyp = nystag
      end if
c
      if ( kb .eq. 0 .or. ke .eq. 0 ) then
      kb = 1
      ke = nzstag
      nzp = nzstag
      end if
c
      if ( ipltyp .ne. 5 .and. 
     >    ((ib) .lt. -001 
     >             .or. 
     >     (jb) .lt. -001 
     >             .or. 
     <     (kb) .lt. -001 ) ) then
      print*,'EXTREMES',ib,jb,kb
c
c  plots at extremes of w+
c
      if ( (ib) .eq. -99 .or. 
     > (jb) .eq. -99 .or. (kb) .eq. -99 ) then
      call extrms(1,ng,1,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,w%flt3d)
      end if
c
c  plots at extremes of qh
c
      if ( (ib) .eq. -98 .or. 
     > (jb) .eq. -98 .or. (kb) .eq. -98 ) then
      call extrms(lh,ng,na,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,st)
      end if
c
c  plots at extremes of qr
c
      if ( (ib) .eq. -97 .or. 
     > (jb) .eq. -97 .or. (kb) .eq. -97 ) then
      call extrms(lr,ng,na,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,st)
      end if
c
c  plots at center of domain
c
      if ( (ib) .eq. -91 ) then
      imx = nx/2
      end if
      if ( (jb) .eq. -91 ) then
      jmx = ny/2
      end if
      if ( (kb) .eq. -91 ) then
      kmx = nz/2
      end if
c
c  plots at level above melting level
c
      if ( (kb) .eq. -94 ) then
      do 6795 kz = 1 ,nz
      if ( templt(kz) .lt. 0.0 ) then 
      kmx = kz 
      go to 6796
      end if
6795  continue
6796  continue
      end if
c
      if ( ib .lt. -90 ) then
      ib = imx 
      ie = imx  
      end if
      if ( jb .lt. -90 ) then
      jb = jmx 
      je = jmx  
      end if
      if ( kb .lt. -90 ) then
      kb = kmx 
      ke = kmx  
      end if
      print*,'EXTREMES',ib,ie,jb,je,kb,ke,xmax,xmin
      end if
      nxp = ie-ib+1
      nyp = je-jb+1
      nzp = ke-kb+1
!      write(*,*) 'DEBUG: nzp = ',nzp,ke,kb
      dzt = z1d2(ke+1,2) - z1d2(kb,2)
      IF ( nzset .gt. 0 ) THEN
       dzt = (hgtset-hgtlow)*1000.
       dzp = hgtlow*1000.
      ELSE
       dzp = z1d2(kb,1) - z1d2(1,1)
      ENDIF
      dxp = (ib-1)*dx
      dyp = (jb-1)*dy
c      dzp = (kb-1)*dz
      ixdiv = ifix(xdiv/dx + .01/dx)
      jydiv = ifix(ydiv/dy + .01/dy)
      kzdiv = ifix(zdiv/dz + .01/dz)
      
      dxp = gxt(ib,1) - gxt(1,1)
      dyp = gyt(jb,1) - gyt(1,1)
      
      ! hack to round to nearest 10km
      IF ( dxp > 250000. ) dxp = 10000.*Nint(dxp/10000.)
      IF ( dyp > 250000. ) dyp = 10000.*Nint(dyp/10000.)
      
      write(6,*) 'dx,dy,dz,dzt,dxp,dyp,dzp = ',dx,dy,dz,dzt,dxp,dyp,dzp
      
c      write(0,*) 'nzset,n2set,mod(n2set-1,kzdiv) = ',
c     :    nzset,n2set,mod(n2set-1,kzdiv)
      write(6,*) 'xdiv,ydiv,zdiv',xdiv,ydiv,zdiv
      write(6,*) 'ixdiv,jydiv,kzdiv',ixdiv,jydiv,kzdiv
c
      write(*,*) 'idebg=10'
c
c  XY PLOT
c
      if ( ipltyp .eq. 1 ) then
       IF ( (idx .eq. 327 .or. idx .eq. 328 .or. idx == 329 .or.
     :     idx .eq. 324 .or. idx .eq. 325 .or. idx == 326 ) .and.
     :     ishade .lt. 0 .and. ( ke .eq. 1 .or. idx == 329 ) )  THEN
         icgsymbol = 2
         ishade = Abs(ishade)
       ENDIF
      d1 = dx
      d1p = dxp + xgpos
      fac1 = xfac
      d02 = dy
      d2p = dyp + ygpos
      fac2 = yfac
      n1 = nxp
      n2 = nyp
      i1 = ixdiv
      i2 = jydiv
      idm = 1
      call getslb(gd,idx,ib,jb,kb,idm,nxp,nyp,nzp
     >  ,nmax,nmax,1,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     >  ,xmx,xmn,cstag,plt2dp,plt2dn,xvec,yvec,zvec,uava,vava,
     :   nx,ny,nz,na,ng,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)
      rx = min(1.0,(sdx*dx*float(nxp))/(sdy*dy*float(nyp)))
      ry = min(1.0,(sdy*dy*float(nyp))/(sdx*dx*float(nxp)))
      write(lloc,911) zb(kb)*zfacb + (1-kstag)*dzadd
      label1 = '                    X-Y SLAB'
      label6 = 'z(km) ='//lloc
      end if 
c
      write(*,*) 'idebg=11'
c
c  XZ PLOT
c
      if ( ipltyp .eq. 2 ) then
       IF ( ishade .lt. 0 .and. ( ke .eq. 1 .or. idx == 329 ) )  THEN
         icgsymbol = 1
         ishade = Abs(ishade)
       ENDIF
      d1 = dx
      d1p = dxp + xgpos
      fac1 = xfac
      d02 = dzt/(nzp-1) ! dz
      d2p = dzp
      fac2 = zfac
      n1 = nxp
      n2 = nzp
      i1 = ixdiv
      i2 = kzdiv

      IF ( nzset .gt. 0 ) THEN
       d02 = dz
       n2 = nzset
      ENDIF
       idm=1
      call getslb(gd,idx,ib,jb,kb,idm,nxp,nyp,nzp
     >  ,nmax,1,nmax,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     >  ,xmx,xmn,cstag,plt2dp,plt2dn,xvec,zvec,yvec,uava,vava,
     :   nx,ny,nz,na,ng,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)
c      rx = min(1.0,(sdx*dx*float(nxp))/(sdz*dz*float(nzp)))
      rx = min(1.0,(sdx*dx*float(nxp))/( sdz*dzt ))
c      ry = min(1.0,(sdz*dz*float(nzp))/(sdx*dx*float(nxp)))
      ry = min(1.0,(sdz*dzt)/(sdx*dx*float(nxp)))
      write(lloc,911) dy*yfac*(float(jb-1) + 0.5)
      label1 = '                    X-Z SLAB'
      label6 = 'y(km) ='//lloc
      print*,'WIERD',label6,lloc,yfac,jb,dy,sf,sdz,sforz
      write(*,*) 'CVC',xfacb, yfacb, zfacb, xdivb, ydivb, zdivb
      end if 
c
      write(*,*) 'idebg=13'
c
c  YZ PLOT
c
      if ( ipltyp .eq. 3 ) then
      d1 = dy
      d1p = dyp + ygpos
      fac1 = yfac
      d02 = dzt/(nzp-1) ! dz
      d2p = dzp
      fac2 = zfac
      n1 = nyp
      n2 = nzp
      i1 = jydiv
      i2 = kzdiv

      IF ( nzset .gt. 0 ) THEN
       d02 = dz
       n2 = nzset
      ENDIF
      idm = 1
      call getslb(gd,idx,ib,jb,kb,idm,nxp,nyp,nzp
     >  ,1,nmax,nmax,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     >  ,xmx,xmn,cstag,plt2dp,plt2dn,zvec,xvec,yvec,uava,vava,
     :   nx,ny,nz,na,ng,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)
c      rx = min(1.0,(sdy*dy*float(nyp))/(sdz*dz*float(nzp)))
      rx = min(1.0,(sdy*dy*float(nyp))/(sdz*dzt))
c      ry = min(1.0,(sdz*dz*float(nzp))/(sdy*dy*float(nyp)))
      ry = min(1.0,(sdz*dzt)/(sdy*dy*float(nyp)))
c      write(lloc,'(f6.2)') dx*xfac*float(ib-1)
      write(lloc,911) dx*xfac*(ib-0.5)
      label1 = '                    Y-Z SLAB'
      label6 = 'x(km) ='//lloc

      IF ( irev .eq. 1 ) THEN
       DO kz = 1,n2
        DO ix = 1,n1
         tmp2d(ix,kz) = plt2dp(ix,kz)
        ENDDO
        DO ix = 1,n1
         plt2dp(ix,kz) = tmp2d(n1-ix+1,kz)
        ENDDO
        DO ix = 1,n1
         tmp2d(ix,kz) = plt2dn(ix,kz)
        ENDDO
        DO ix = 1,n1
         plt2dn(ix,kz) = tmp2d(n1-ix+1,kz)
        ENDDO
        DO ix = 1,n1
         tmp2d(ix,kz) = yvec(ix,kz)
        ENDDO
        DO ix = 1,n1
         yvec(ix,kz) =  tmp2d(n1-ix+1,kz)
        ENDDO
        DO ix = 1,n1
         tmp2d(ix,kz) = xvec(ix,kz)
        ENDDO
        DO ix = 1,n1
         xvec(ix,kz) = -tmp2d(n1-ix+1,kz)
        ENDDO
       ENDDO
      
      ENDIF

      end if 
c
c  OBLIQUE PLOT
c
      if ( ipltyp .eq. 5 ) then
      d1 = dx   ! ok
c      d1p = (xib-1.0)*Cos(phi*pi/180.0)*dx + 
c     :      (yjb-1.0)*Sin(phi*pi/180.0)*dy
      d1p = 0.0
      d1p=d1pobl
      fac1 = xfac
      d02 = dzt/(nzp-1) ! dz
      d2p = dzp
      fac2 = zfac
      n1 =  nl    ! ok
      n2 = nzp   ! ok
      i1 = ixdiv ! ok
      i2 = kzdiv ! ok

      IF ( nzset .gt. 0 ) THEN
       d02 = dz
       n2 = nzset
      ENDIF
      
      DO j=1,nmax
      DO i=1,nmax
      plt2dp(i,j) = 0.0
      plt2dn(i,j) = 0.0
      ENDDO
      ENDDO

c      write(0,*) '1: idx,xmx,xmn',idx,xmx,xmn,ivint
      idm = 1
      call getslb(gd,idx,1,1,kb,idm,nx,ny,nzp
     >  ,nmax,1,nmax,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     >  ,xmx,xmn,cstag,plt2dp,plt2dn,xvec,zvec,yvec,uava,vava,
     :   nx,ny,nz,na,ng,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)
     
c      write(0,*) '2: idx,xmx,xmn',idx,xmx,xmn,ivint
      
c      rx = min(1.0,(sdx*dx*float(nxp))/(sdz*dz*float(nzp)))
c      ry = min(1.0,(sdz*dz*float(nzp))/(sdx*dx*float(nxp)))
c      rx = min(1.0,(sdx*dx*float(nxp))/(sdz*dzt))
      rx = min(1.0,(sdx*dx*float(nl))/(sdz*dzt))
      ry = min(1.0,(sdz*dzt)/(sdx*dx*float(nl)))
      write(6,*) 'nxp = ',nxp
c      write(lloc,'(f5.1)') dy*yfac*float(jb-1)
      lloc = '         '
      label1 = '                   XY-Z SLAB'
c      label6 = 'x,y(km),t:'//lloc
      write(label6,'(a,f5.1,1x,f5.1,1x,f5.1)') 
     :  'x,y,t:',dx*xfac*xib ,dy*yfac*(yjb-1.0),phi
      print*,'WIERD',label6,lloc,yfac,jb,dy
      write(*,*) 'CVC',xfacb, yfacb, zfacb, xdivb, ydivb, zdivb
      end if 
c
c  ISOSFC
c
      if ( ipltyp .eq. 4 ) then
      d1 = dx
      d1p = dxp
      fac1 = xfac
      d1 = dy
      d1p = dyp
      fac1 = yfac
      d02 = dz
      d2p = dzp
      fac2 = zfac
      n1 = nxp
      n2 = nyp
      n3 = nzp
      i1 = ixdiv
      i2 = jydiv
      i3 = kzdiv
      idm = 1
      call getslb(gd,idx,ib,jb,kb,idm,nxp,nyp,nzp
     >  ,nx,ny,nz,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     >  ,xmx,xmn,cstag,plt3dp,plt3dn,plt3dn,plt3dn,plt3dn,uava,vava,
     :   nx,ny,nz,na,ng,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)
      if ( dx*nxp .gt. dy*nyp ) then
c      rx = min(1.0,(sdx*dx*float(nxp))/(sdz*dz*float(nzp)))
c      ry = min(1.0,(sdz*dz*float(nzp))/(sdx*dx*float(nxp)))
      rx = min(1.0,(sdx*dx*float(nxp))/(sdz*dzt))
      ry = min(1.0,(sdz*dzt)/(sdx*dx*float(nxp)))
      else
c      rx = min(1.0,(sdy*dy*float(nyp))/(sdz*dz*float(nzp)))
c      ry = min(1.0,(sdz*dz*float(nzp))/(sdy*dy*float(nyp)))
      rx = min(1.0,(sdy*dy*float(nyp))/(sdz*dzt))
      ry = min(1.0,(sdz*dzt)/(sdy*dy*float(nyp)))
      end if
      write(lloc,911) tiso
      label1 = '                     SURFACE'
      label6 = 'surface value='//lloc
      end if
c
      write(*,*) 'idebg=14'
c
c
c  TZ PLOT (time-height)
c
      if ( ipltyp .eq. 10 ) then
         IF ( ishade .lt. 0 .and. 
     :       ( Abs(idx).eq.912 .or. Abs(idx).eq.913 .or.
     :         Abs(idx).eq.914 )) THEN
             icgsymbol = 1
             ishade = Abs(ishade)
         ENDIF
      IF ( Abs(idx) .lt. 901 ) THEN
        write(6,*) 'Bad idx for ipltyp=10! idx = ',idx
        STOP
      ENDIF
      dt = timearr(2) - timearr(1)
      write(6,*) 'dt = ',dt
c      dx = dt
      IF ( dt .eq. 1.0 ) THEN
        dxp = (ib-1)*dt 
      ELSE
        ib = Nint( (ib - 1)/dt ) + 1
        dxp = (ib-1)*dt 
        nxp = Nint( (nxp - 1)/dt )
        ie  = Nint( ie/dt )
      ENDIF
      d1 = dt
      d1p = dxp
      fac1 = xfac
      d02 = dzt/(nzp-1) ! dz
      d2p = dzp
      fac2 = zfac
      n1 = nxp+1
      n2 = nzp

      IF ( nzset .gt. 0 ) THEN
       d02 = dz
       n2 = nzset
      ENDIF
      write(6,*) 'dxp,etc = ',dxp,ib,fac1,d02,n1,n2

      itdiv = ifix(xdiv/dt + .01/dt)
      i1 = itdiv
      i2 = kzdiv
      
      write(6,*) 'itdiv,i1,i2 = ',itdiv,i1,i2,xdiv,dt
      
c      write(lunw,200) 'xmax=',xmx,' imx=',imx,' jmx=',jmx, 'kmx=',kmx
c      write(lunw,200) 'xmin=',xmn,' imn=',imn,' jmn=',jmn, 'kmn=',kmn

      xmax = 0.0
      xmax2 = 0.0
      xmax3 = 0.0
      xmin = 0.0
      jmx = 1
      jmn = 0
      imx = 1
      imn = 0
      kmx = 1
      kmn = 0
      
      DO j=1,nmax
      DO i=1,nmax
      plt2dp(i,j) = 0.0
      plt2dn(i,j) = 0.0
      ENDDO
      ENDDO
      
      icontent = 0
      IF ( (idx .ge. 1131 .and. idx .le. 1154) .or. idx == 1158 ) THEN
        icontent = 1
        idx0 = 1130
        
        DO kz = 1,nz-1
         zdb(kz) = db(kz)
        ENDDO
      ELSE
        DO kz = 1,nz-1
          zdb(kz) = 1.0
        ENDDO
      ENDIF
      IF ( thlengthscale < 0.0 ) THEN
        scalez(:) = 1.0
      ELSE
        DO kz = 1,nz-1
          scalez(kz) = thlengthscale*gzt(kz,3) ! normalize from per grid level to per length 
        !  write(0,*) 'scalez = ',kz,scalez(kz)
        ENDDO
        scalez(nz) = scalez(nz-1)
      ENDIF
      itshft = Nint( timearr(1)/dt ) - 1
      write(6,*) ' itshft,ib,ie = ',itshft,ib,ie
      xmax3 = 0
      DO k=1,nzp
       kz = (k - 1) + kb
       DO ix=Max(ib,1),ie+1
       i = ix - ib + 1
       IF ( ix .gt. itshft ) THEN
        ix2 = Min(ie-itshft,ix-itshft-1)
        ix2 = Max(1,ix2)
        xmax2 = Max(xmax2, charr(kz,ix2,idx-idx0)*zdb(kz)*scalez(kz) )
        IF ( charr(kz,ix2,idx-idx0) .ge. 0.0 ) THEN
          plt2dp(i,k) = charr(kz,ix2,idx-idx0)*zdb(kz)*scalez(kz)
c
c The following code keeps the contour lines from going off the bottom
c  of the plot (keeps a closed contour around CGs)
c
           IF ( (idx .eq. 902 .or. idx .eq. 903) .and. ishade .lt. 100
     :      .and. kz .eq. 1 ) THEN
c            write(0,*) idx,ishade,kz,charr(kz,ix2,idx-idx0)
             plt2dp(i,k) = Min(0.7, charr(kz,ix2,idx-idx0) )
           ENDIF
       !   xmax3 = Max(xmax3, plt2dp(i,k)  )
          IF ( plt2dp(i,k) .gt. xmax .and. k <= n2 ) THEN
           xmax = plt2dp(i,k)
           IF ( k < nzp ) THEN
            xmax3 = xmax
           ENDIF
c           write(6,*) 'xmax = ',xmax
           imx = ix
           kmx = kz
          ENDIF
        ELSE
          plt2dn(i,k) = charr(kz,ix2,idx-idx0)*zdb(kz)
           IF ( (idx .eq. 902 .or. idx .eq. 903) .and. ishade .lt. 100
     :      .and. kz .eq. 1 ) THEN
           plt2dp(i,k) = Max(-0.7, charr(kz,ix2,idx-idx0) )
           ENDIF
          IF ( plt2dn(i,k) .lt. xmin ) THEN
           xmin = plt2dn(i,k)
           imn = ix
           kmn = kz
          ENDIF
        ENDIF
       ENDIF
       
!       write(6,*) 'i,k,plt2dp = ',i,k,kz,plt2dp(i,k)
       ENDDO
      ENDDO
      
      IF ( icontent .eq. 1 ) THEN
        icontent = 0
        idx0 = 1100
      ENDIF

      xmx = xmax
      xmn = xmin
      write(6,200) 'xmax=',xmx,' imx=',imx,' jmx=',jmx, 'kmx=',kmx
      write(6,200) 'xmin=',xmn,' imn=',imn,' jmn=',jmn, 'kmn=',kmn
      write(6,*) 'xmax2,xmax3 = ',xmax2,xmax3,plt2dp(imx,kmx)

      IF ( ioverlay .ge. 1 ) THEN
       xmax4 = 0.0
       xmax3 = 0.0
       xmax2 = 0.0
      ENDIF

c      call getslb(idx,ib,jb,kb,ivint,nxp,nyp,nzp
c     >  ,nmax,1,nmax,istag,jstag,kstag
c     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
c     >  ,xmx,xmn,cstag,plt2dp,plt2dn,xvec,zvec,yvec,uava,vava)
c      rx = min(1.0,(sdx*dx*float(nxp))/( sdz*dzt ))
      rx = min(1.0,(sdx*dt*500.*float(nxp))/( sdz*dzt ))
c      ry = min(1.0,(sdz*dzt)/(sdx*dx*float(nxp)))
      ry = min(1.0,(sdz*dzt)/(sdx*dt*500.*float(nxp)))
      write(lloc,911) dy*yfac*float(jb-1)
      label1 = '                   Time-Height'
      label6 = 'Time-Height'
      print*,'WIERD',label6,lloc,yfac,jb,dy,sf,sdz,sforz
      write(*,*) 'CVC',xfacb, yfacb, zfacb, xdivb, ydivb, zdivb
      end if 
c
      write(*,*) 'idebg=14a'

c
c
c  segments-height PLOT (vertical skyline)
c
      if ( ipltyp .eq. 11 ) then
      IF ( Abs(idx) .lt. 901 ) THEN
        write(6,*) 'Bad idx for ipltyp=10! idx = ',idx
        STOP
      ENDIF
      dt = timearr(2) - timearr(1)
      write(6,*) 'dt = ',dt
c      dx = dt
      dxp = 0.0
      d1 = dt
      d1p = 0.0
      d02 = dzt/(nzp-1) ! dz
c      write(0,*) 'd02,dz= ',d02,dz
      d2p = dzp
      fac2 = zfac
      n1 = nxp+1
      n2 = nzp

      IF ( nzset .gt. 0 ) THEN
       d02 = dz
       n2 = nzset
      ENDIF

      i2 = kzdiv
      

      xmax = 0.0
      xmax2 = 0.0
      xmax3 = 0.0
      xmin = 0.0
      jmx = 1
      jmn = 0
      imx = 1
      imn = 0
      kmx = 1
      kmn = 0
      
      DO kz=1,nz
       segheight(kz) = 0.0
      ENDDO
      
      DO j=1,nmax
      DO i=1,nmax
       plt2dp(i,j) = 0.0
       plt2dn(i,j) = 0.0
      ENDDO
      ENDDO
      
      itshft = Nint( timearr(1) ) - 1
      write(6,*) ' itshft,ib,ie = ',itshft,ib,ie
      DO kz=1,nz
       DO ix=Max(ib,1),ie
       i = ix - ib + 1
       IF ( ix .gt. itshft ) THEN
        ix2 = Min(ie-itshft,ix-itshft)
        IF ( charr(kz,ix2,idx-idx0) .ge. 0.0 ) THEN
          segheight(kz) = segheight(kz) + charr(kz,ix2,idx-idx0)
          IF ( segheight(kz) .gt. xmax ) THEN
           xmax = segheight(kz)
           imx = ix
           kmx = kz
          ENDIF
        ELSE
          IF ( plt2dn(i,kz) .lt. xmin ) THEN
           xmin = plt2dn(i,kz)
           imn = ix
           kmn = kz
          ENDIF
        ENDIF
       ENDIF
       ENDDO
      ENDDO

      xmx = xmax
      xmn = xmin
      write(6,*) 'xmax=',xmx,' kmx=',kmx,' height= ',heights(kmx)
      write(6,200) 'xmin=',xmn,' imn=',imn,' jmn=',jmn, 'kmn=',kmn
      

c
c Note: 
c  xdiv = width of major division
c  d1 = width of minor division
c  xfac = scaling factor (usually 0.001 or 1.0 ) (dylab can show scale factor)
c
c      IF ( (Int(xmax/4000) + 1)*4000 .ge. 100000. ) THEN
c        xdiv = 2000.
c        d1 = 10000.
c        xfac = .001
c        n1 = Int(xmax/xdiv) + 1
c        xmax2 = xdiv*(n1)
c        n1 = Nint(xmax2/d1) + 1
c      ELSEIF ( xmax .ge. 10000. ) THEN
      IF ( xmax .ge. 100000. ) THEN
        n1 = 2
        xdiv = (Int(xmax/20000) + 1)*20000/2
        d1 = 10000.0 ! xdiv/10.
        xfac = 1.
c        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1

      ELSEIF ( xmax .ge. 10000. ) THEN
        n1 = 2
        xdiv = (Int(xmax/4000) + 1)*4000/2
        d1 = 2000.0 ! xdiv/10.
        xfac = 1.
c        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 5000. ) THEN
        xdiv = 2000.
        d1 = 1000.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 2000. ) THEN
        xdiv = 1000.
        d1 = 200.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 1400. ) THEN
        xdiv = 500.
        d1 = 100.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 800. ) THEN
        xdiv = 200.
        d1 = 50.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 400. ) THEN
        xdiv = 100.
        d1 = 20.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 150. ) THEN
        xdiv = 50.
        d1 = 10.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 30. ) THEN
        xdiv = 10.
        d1 = 5.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSEIF ( xmax .ge. 15. ) THEN
        xdiv = 5.
        d1 = 1.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ELSE
        xdiv = 2.
        d1 = 1.
        xfac = 1.
        n1 = Int(xmax/xdiv) + 1
        xmax2 = xdiv*(n1)
        n1 = Nint(xmax2/d1) + 1
      ENDIF
      write(6,*) 'xmax2,xmax3 = ',xmax2,xmax3,plt2dp(imx,kmx)

      fac1 = xfac
      
c      write(0,*) 'n1,xfac,d1,xdiv,xmax2= ',n1,xfac,d1,xdiv,xmax2
      
      write(6,*) '1:xmax4= ',xmax4
      IF ( xmax2 .gt. xmax4 ) THEN
        IF ( xmax4 .eq. 0.0 ) THEN
          xmax4 = xmax2
      dt = 1.0
c      d1 = 1.0/x
      itdiv = ifix(xdiv/d1 + .01/d1)
      i1 = itdiv
        ELSE
          write(0,*) 'segheight overlay has larger xmax2! Reorder!'
          STOP
        ENDIF
      ENDIF

      write(6,*) '2:xmax4= ',xmax4
      

       hgt = 1000.0 

      write(6,*) 'n1,i1,hgt = ',n1,i1,hgt

      skyline(1,1) = 0.0
      skyline(1,2) = 0.0
c      skyline(2,1) = 0.0
c      skyline(2,2) = segheight(1)
      DO kz=2, 2*nz,2
        skyline(kz,1) = skyline(kz-1,1)
        skyline(kz,2) = segheight(kz/2)/xmax4
        
        skyline(kz+1,1) = z1d2(kz/2,2)/hgt
        skyline(kz+1,2) = skyline(kz,2)
      ENDDO

      lskyb = .true.
      iskyb = 1
      DO kz=2, 2*nz,2
       lskyb = ( lskyb .and. skyline(kz,2) .eq. 0.0 )
       IF ( lskyb ) iskyb = kz
      ENDDO

      lskyt = .true.
      iskyt = 2*nz
      
      DO kz=2*nz, 2, -2
       lskyt = ( lskyt .and. skyline(kz,2) .eq. 0.0 )
       IF ( lskyt ) iskyt = kz
      ENDDO
      
      IF ( iskyt .lt. iskyb ) THEN
       iskyb = iskyt
      ENDIF
      
      IF ( ioverlay .ge. 1 ) THEN
       xmax4 = 0.0
       xmax2 = 0.0
       iclr = 1
      ENDIF
      
      ZCD(1) = 0.0
      ZCD(2) = 1.
      XCD(1) = 0.0
      XCD(2) = 1.
      
      IF ( ioverpnt .eq. 1 ) ioverpnt2 = 3

c      rx = min(1.0,(sdx*dt*500.*float(nxp))/( sdz*dzt ))
c      ry = min(1.0,(sdz*dzt)/(sdx*dt*500.*float(nxp)))
      rx = 1.0
      ry = 1.0
      write(lloc,911) dy*yfac*float(jb-1)
      label1 = '                   Time-Height'
      label6 = 'Total-Height'
      print*,'WIERD',label6,lloc,yfac,jb,dy,sf,sdz,sforz
      write(*,*) 'CVC',xfacb, yfacb, zfacb, xdivb, ydivb, zdivb
      
c      STOP
      
      end if 
c
      write(*,*) 'idebg=14b'
      
      
      if ( .not. (idx .eq. 23 .or. idx .eq. 32)) then
      write(lunw,*) 'MAXIMUM and MINIMUM of PLOT, idx=',idx
      write(lunw,200) 'xmax=',xmx,' imx=',imx,' jmx=',jmx, 'kmx=',kmx
      write(lunw,200) 'xmin=',xmn,' imn=',imn,' jmn=',jmn, 'kmn=',kmn
  200 format(1x,a5,1pe13.5,3(a5,i5.5))
      end if
c
      axb = 0.15
      ayb = 0.20
      rxy = max(rx,ry)
      axe = axb + rx*(0.55)/rxy
      aye = ayb + ry*(0.55)/rxy
      axep =  0.85-axb
      if ( ilabor .eq. 0 ) then
      xbarr = axe
      xbarl = axb
      ybart = 0.10
      ybarb = 0.05
      end if
      if ( ilabor .eq. 1 ) then
        IF ( itlab .ge. 1 .and. itlab .le. 3 ) THEN
          xbarr = Min(axe + 0.025 + 0.1 , .95)
          xbarl = Min(axe + 0.025, .85)
        ELSE
          xbarr = Min(axe + 0.001 + 0.1 , .95)
          xbarl = Min(axe + 0.001, .85)
        ENDIF
      ybart = aye
      ybarb = ayb
      end if
      label = title(Abs(idx))
c      call gsclip(0)
      call set(0.,1.0,0.,1.0,0.,1.0,0.,1.0,1)
c
      write(*,*) 'ioverlay',ioverlay
c     if (ioverlay .eq. -2 ) then
      write(*,*) 'idebg=14a'
      write(*,*) 'ioverlay ,ioverlayplt map stuff'
c
      iexp = 0
        IF ( xmult .eq. 0.0 ) THEN
          tmp = 0.5*( Abs(xmx)  + Abs(xmn) )
          IF ( tmp .gt. 1.e+5 .or. tmp.lt.0.1 ) THEN
           iexp =  -3*(Int( Log10(tmp) )/3) 
c           write(0,*) 'idx,xmx,xmn,tmp,iexp=',idx,xmx,xmn,tmp,iexp
           IF ( iexp .gt. 2 ) iexp = iexp + 3
           xmult = 10.0**( iexp )
           rchip = rchip*xmult
           rclon = rclon*xmult
           rcinp = rcinp*xmult
          ENDIF
        ENDIF
      IF ( rcinp .eq. 0.0 .and. rchip .eq. 0.0 .and. rclon.eq.0.0) THEN
        rchip = xmx*Abs(xmult)
        rclon = xmn*Abs(xmult)
        nclv = 8
         tmp = (rchip - rclon)/(nclv-1)
        rcinp = 0.9*tmp
          rchip = rchip - 0.05*tmp
          rclon = rclon + 0.05*tmp
        IF ( xmult .eq. 1.0 ) THEN
          tmp = 0.5*( Abs(xmx)  + Abs(xmn) )
          IF ( (tmp .gt. 1.e+4 .or. tmp.lt.0.1) .and.
     :         Abs(tmp).gt.1e-15   ) THEN
           iexp =  -3*(Int( Log10(tmp) )/3) 
           IF ( iexp .gt. 2 ) iexp = iexp + 3
           IF ( Log10(tmp) .gt. -3.0 .and. Log10(tmp) .lt. -1. ) THEN
            iexp = 3
           ENDIF
           write(6,*) 'idx,xmx,xmn,tmp,iexp=',idx,xmx,xmn,tmp,iexp
           xmult = 10.0**( iexp )
           rchip = rchip*xmult
           rclon = rclon*xmult
           rcinp = rcinp*xmult
          ENDIF
        ENDIF
        rclop = rclon
        write(6,*) 'rchip,rclop,rcinp,xmult= ',rchip,rclop,rcinp,xmult
      ENDIF
      write(lcnt,911) rcinp
      write(ltim,911) time
      write(lmaxvl,912) xmx*xmult
      IF ( xmx*xmult .ge. 10000. ) THEN
       write(lmaxvl,'(f9.0)') xmx*xmult
      ENDIF

      IF ( Abs(xmult) .ne. 1.0 ) THEN
       idum = len(label)
       CALL strlnth(label,idum)
       write(titlelabel, '(a,1pe7.1)') 'x ',1./xmult
       write(6,*) 'LABEL = ',titlelabel
       label = label(1:idum+1)//titlelabel
       write(6,*) 'LABEL = ',label
      ENDIF

       ilabel = len(label)
       CALL strlnth(label,ilabel)

      write(lminvl,912) xmn*xmult
c
      do 876 m = 0,100
      awrk = rclon + m*rcinn
      if ( awrk .gt. xmn*xmult ) then
      xmn = awrk
      go to 877
      end if
 876  continue
 877  continue
c      if ( rclon .gt. -0.00025 ) awrk = rclon
      write(lmin,911) max( rclon, awrk )
c
 912  format(f9.2)
 911  format(f8.3)
c
      label2 = 'TIME(MIN) = '//ltim
      label3 = 'cntr interval ='//lcnt
      label4 = 'min cntr ='//lmin
      label5 = 'mx/mn='//lmaxvl//'/'//lminvl

c
      if ( ipltyp .ne. 4 ) then
      if ( ioverpnt .eq. 1 ) then
      PRINT*,'LABEL6',label6
      IF ( iplotidx > 0) THEN
      if ( ipware .ge. 2) then
       if ( ilabor .eq. 0 ) then
         ! print field name over plot
         call pwrit(axb,(ayb-.18),label,ilabel,lets,0,-1)
       end if
       if ( ilabor .eq. 1 ) then
       call pwrit((.98),ayb,label,ilabel,lets,90,-1)
       end if
      end if
c      call pwrit(axb,(aye+.10),label,40,lets,0.0,-1)
      call pwrit(axb,(aye+.055),label,ilabel,lets,0,-1)
c      call pwrit(axep,(aye+0.13),label1,28,lets,0,1)
c      call pwrit(axep,(aye+0.07),label5,28,lets,0,1)
c      call pwrit(axb ,(aye+0.07),label6,28,lets,0,-1)
      ENDIF
      
      
c erm strln
      titlen = len(label5) ! z(km)  string (on line with mx/mm)
      CALL strlnth(label5,titlen)
      IF ( iplotidx > 0 ) THEN
        call pwrit(axep,(aye+0.035),label5,titlen,lets,0,1)
      ENDIF
      titlen = len(label6) ! plot field index name
      CALL strlnth(label6,titlen)
c      write(0,*)
      IF ( iplotidx > 0 ) THEN
        call pwrit(axb ,(aye+0.035),label6,titlen,lets,0,-1)
      ENDIF
      if ( ipmore .gt. 0 ) then
c      call pwrit(axb,(aye+0.10),stitle,60,min(3,2*lets),0,-1)
        IF ( iplottitle == 1 ) THEN
          call pwrit(axb,(aye+0.10),xtitle,xtitlen,min(3,2*lets),0,-1)
        ELSEIF ( iplottitle > 1 ) THEN
          call pwrit(axb,(aye+0.018),xtitle,xtitlen,min(3,2*lets),0,-1)
         ENDIF
         
      IF ( ifile2 .ge. 0 ) THEN
        IF ( itimeoffset .eq. 0 ) THEN
          IF ( iplottime == 1 ) THEN
           call pwrit(axb,(aye+0.075),
     :     stimeplot//' s ='//smin//' min ='//shr//' hr',33,
     :      1,0,-1)
           ELSEIF ( iplottime > 1 ) THEN

! yo + 2.0*yvalp
           call pwrit(axb,(ayb-0.1+0.06),
     :     trim(smin2)//' min',8,
     :      lets*2,0,-1)
           
           ENDIF
         ELSE
        call pwrit(axb,(aye+0.075),'Forecast time ='//smin//' min',26,
     :   1,0,-1)
         ENDIF
      ENDIF
      
      
c      call pwrit(axep,(aye+0.10),
c     >  '               3-D',18,lets,0.0,1)
c      call pwrit(axb ,(aye+0.13),label2,28,lets,0.0,-1)
c      call pwrit(axb ,(aye+0.04),label3,28,lets,0.0,-1)
c      call pwrit(axep,(aye+0.04),label4,28,lets,0.0,1)
      IF ( ipltyp .ne. 11 .and. iplotcint .ge. 1 ) THEN
      titlen = len(label3)
      CALL strlnth(label3,titlen)
      call pwrit(axb ,(aye+0.015),label3,titlen,lets,0,-1)
      titlen = len(label4)
      CALL strlnth(label4,titlen)
      call pwrit(axep,(aye+0.015),label4,titlen,lets,0,1)
      ENDIF

      end if
      end if
c
      if ( ipltyp .ne. 4 ) then
      if ( ioverpnt .eq. 1 ) then
      if ( ipmore .gt. 0 ) then
      if ( ixlab .gt. 0 ) then
c      subroutine dxlabs(it,xo,xn,yo,yn,nx,yvalp,ixdiv,dx,dxp,xfac)
c      write(6,*) 'dxlabs: ',ipltyp,axb,axe,(ayb-.1),aye,n1
c     >           ,i1,d1,d1p,fac1
      
      write(6,*) 'call dxlabs: ',ipltyp,axb,axe,(ayb-.1),aye,n1,0.03
     >           ,i1,d1,d1p,fac1,irev
      call dxlabs(ipltyp,axb,axe,(ayb-.1),aye,n1,0.03
     >           ,i1,d1,d1p,fac1,irev)
      
      end if
c     if ( ipltyp .eq. 2 .or. ipltyp .eq. 3 ) d2p = dzadd
      if ( iylab .gt. 0 ) then
c      subroutine dylabs(it,yo,yn,xo,xn,ny,xvalp
c     >                 ,jydiv,dy,dyp,yfac,itlab,templt,nz)
      IF ( ipltyp .ne. 1 ) THEN
c      write(0,*) 'axb,axe,ayb,aye= ',axb,axe,ayb,aye
c      write(0,*) 'nzset,d02,d2p,nz= ',nzset,d02,d2p,nz
        IF ( nzset .ne. 0 ) THEN
c        write(0,*) 'call dylabs1, ipltyp = ',ipltyp,itlab
          call dylabs(ipltyp,ayb,aye,0.0,axe,n2set,0.03,axb
     >           ,i2,d02,d2p,fac2,itlab,templt,zk,nz,kb1)
        ELSE
c        write(0,*) 'call dylabs2, ipltyp = ',ipltyp,itlab
          call dylabs(ipltyp,ayb,aye,0.0,axe,n2,0.03,axb
     >           ,i2,d02,d2p,fac2,itlab,templt,zk,nz,kb1)
        ENDIF
      ELSE
c        write(0,*) 'call dylabs2, ipltyp = ',ipltyp,itlab
            call dylabs(ipltyp,ayb,aye,0.0,axe,n2,0.03,axb
     >           ,i2,d02,d2p,fac2,itlab,templt,zk,ny,kb1)
      ENDIF
      end if
      end if
      end if
      end if
      write(*,*) 'idebg=14b'
cc
c     end if
c
      end if
c
      write(*,*) 'idebg=15'
c
      call set(axb,axe,ayb,aye,0.,1.,0.,1.,1)
      if ( ipltyp .ne. 4 ) then
      if ( ioverpnt2 .eq. 3 ) then
      if ( ipmore .gt. 0 ) then
      if ( mod(n1-1,i1) .eq. 0 
     >          .and. 
     <    ( mod(n2-1,i2) .eq. 0 .or. 
     :      (mod(n2set-1,i2).eq.0 .and. n2set.ne.0) ) ) then
      mjrx = n1/i1
      mnrx = i1
      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
      IF ( n2set .ne. 0 ) THEN
      mjry = n2set/i2
      ELSE
      mjry = n2/i2
      ENDIF
      mnry = i2
      IF ( jymdiv .ne. 0 ) mnry = jymdiv
      IF ( (ipltyp .eq. 2 .or. ipltyp .eq. 3 .or. ipltyp == 5) 
     :       .and. kzmdiv .ne. 0 ) mnry = kzmdiv
      write(6,*) '1 mnrx,mnry,i1,i2 = ',mnrx,mnry,i1,i2
c      write(0,*) 'here1'
      call tick4(itckji,itckni,itckji,itckni)
      call perim(mjrx,mnrx,mjry,mnry)
      call tick4(-itckjo,-itckno,-itckjo,-itckno)
      call perim(mjrx,mnrx,mjry,mnry)
      ELSEIF ( mod(n1-1,i1) .eq. 0 
     >          .and. 
     <    ( mod(n2-1,2) .eq. 0.or. 
     :      (mod(n2set-1,2).eq.0 .and. n2set.ne.0) ) ) then
      mjrx = n1/i1
      mnrx = i1 
      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
      IF ( n2set .ne. 0 ) THEN
      mjry = n2set/i2
      ELSE
      mjry = n2/i2
      ENDIF
c      mjry = n2/2
      mnry = 2
c      write(0,*) 'here2'
      write(6,*) '2 mnrx,mnry,i1,i2 = ',mnrx,mnry,i1,i2
      call tick4(itckji,itckni,itckji,itckni)
      call perim(mjrx,mnrx,mjry,mnry)
      call tick4(-itckjo,-itckno,-itckjo,-itckno)
      call perim(mjrx,mnrx,mjry,mnry)
      ELSEIF( mod(n1-1,i1) .eq. 0  ) THEN
      mjrx = n1/i1
      mnrx = i1 
      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
c      write(0,*) 'here3'
      call tick4(itckji,itckni,itckni,itckji)
      write(6,*) '3 mnrx,mnry,i1,i2 = ',mnrx,mnry,i1,i2
        IF ( n2set .ne. 0 ) THEN
         call perim(mjrx,mnrx,(n2set-1),0)
        ELSE
         call perim(mjrx,mnrx,(n2-1),0)
        ENDIF
      else
c      write(0,*) 'here4'
      call tick4(itckji,itckni,itckji,itckni)
c      call perim(n1/i1,3,n2/i2,3)
c      call perim((n1-1),0,(n2-1),0)
      write(6,*) '4 mnrx,mnry,i1,i2 = ',mnrx,mnry,i1,i2
        IF ( n2set .ne. 0 ) THEN
         call perim(2,(n1-3),2,(n2set-3))
        ELSE
         call perim(2,(n1-3),2,(n2-3))
        ENDIF
      end if 
      end if
      if ( ipmore .eq. 0 ) then
      call perim((1),0,(1),0)
      end if
      end if
      end if
c
      if ( ipltyp .ne. 4 ) then
      IF ( ipltyp .eq. 2 .or. ipltyp .eq. 3 .or.
     :     ipltyp .eq. 5 .or. ipltyp .eq. 10
     :      ) THEN
       write(6,*) 'set map transformation'
       dd = nl*dx
       
       IF ( hgtset .eq. 0.0 ) THEN
       hgt = 0.001*(gzt(n2,1) + 
     :   0.5*(gzt(n2+1,1)-gzt(n2,1)))
       hgt0 = 0.0 ! gt(1,1,1,mzdist)
       ELSE
        hgt = hgtset
        hgt0 = hgtlow
       ENDIF
       write(6,*) 'axb,axe,ayb,aye,hgt= ',axb,axe,ayb,aye,hgt
cc       call set(axb,axe,ayb,aye,0.,dd,0.,hgt,1)
c       call set(axb,axe,ayb,aye,0.,1.,0.,hgt,1)
c test version with hgt0:
       call set(axb,axe,ayb,aye,0.,1.,hgt0,hgt,1)
       CALL CPSETI('SET - DO SET CALL FLAG',0)
       CALL CPSETR('XC1 - X COORDINATE AT INDEX 1',0.)
c       CALL CPSETR('XCM - X COORDINATE AT INDEX M',dd)
       CALL CPSETR('XCM - X COORDINATE AT INDEX M',1.0)

c
c set YC1 = YCN and MAP=3 so that yinp is in range 1-N
c  for CPMPXY routine
c
c       istretch = .false.
       IF ( istretch ) THEN
        CALL CPSETR('YC1 - Y COORDINATE AT INDEX 1',1.)
        CALL CPSETR('YCN - Y COORDINATE AT INDEX N',1.)
        CALL CPSETI ('MAP - MAPPING FLAG',3)
       ELSE
        CALL CPSETI ('MAP - MAPPING FLAG',0)
       ENDIF
      
      ELSEIF ( ipltyp .eq. 11 ) THEN
       IF ( hgtset .eq. 0.0 ) THEN
        hgt = 0.001*(gzt(n2,1) + 
     :    0.5*(gzt(n2+1,1)-gzt(n2,1)))
        hgt0 = 0.0
       ELSE
        hgt = hgtset
        hgt0 = hgtlow
       ENDIF
c       write(6,*) 'axb,axe,ayb,aye,hgt= ',axb,axe,ayb,aye,hgt
c       call set(axb,axe,ayb,aye,0.,1.,0.,hgt,1)
       call set(axb,axe,ayb,aye,0.,1.,hgt0,hgt,1)
      ELSE
       call set(axb,axe,ayb,aye,0.,1.,0.,1.,1)
       CALL CPSETI('SET - DO SET CALL FLAG',0)
       CALL CPSETR('XC1 - X COORDINATE AT INDEX 1',0.)
       CALL CPSETR('XCM - X COORDINATE AT INDEX M',1.)
       CALL CPSETR('YC1 - Y COORDINATE AT INDEX 1',0.)
       CALL CPSETR('YCN - Y COORDINATE AT INDEX N',1.)
       CALL CPSETI ('MAP - MAPPING FLAG',0)
      
      ENDIF
      
      IF ( ipltyp .eq. 11 ) THEN
c
c set color of line
c       
       IF ( Abs( ishade ) .ge. 20 .and. Abs( ishade ) .lt. 100 ) THEN
         ishdclr = Abs( ishade )
         call GSPLCI(ishdclr)
        ELSEIF ( ishdclr .eq. 0 .and. iotyp .ge. 20 ) THEN
         call GSPLCI(99)
       ENDIF
c
c set width of line
c
       CALL GSLWSC(cll)

c
c draw curve
c       
c       CALL CURVE(skyline(1,2),skyline(1,1),2*(nzp+1)+1)
c       CALL CURVE(skyline(iskyb,2),skyline(iskyb,1),iskyt+1)
       CALL CURVE(skyline(iskyb+1,2),skyline(iskyb+1,1),iskyt-iskyb)

c
c reset color to black (default line color) and width to 1.
c
       call GSPLCI(1)
       CALL GSLWSC(1.0)
       
!      ELSEIF ( icgsymbol .eq. 1 .or. (ipltyp == 2 .and. idx == 329  )) THEN ! {
      ELSEIF ( icgsymbol .eq. 1 ) THEN 
      
      icgsymbol = 0
        IF ( idx .eq. 913 .or. idx == 329 ) THEN  ! neg CGs
c          write(0,*) 'set cross marker'
          CALL GSMK( 5 )
        ELSEIF ( idx .eq. 914 ) THEN ! pos CGs
c          write(0,*) 'set plus marker'
          CALL GSMK( 2 )
        ELSEIF ( idx .eq. 912 ) THEN ! IC
c          write(0,*) 'set circle marker'
          CALL GSMK( 4 )
        ENDIF
c
c set color of polymarker
c
        IF ( Abs( ishade ) .ge. 30 .and. Abs( ishade ) .lt. 100 ) THEN
         ishdclr = Abs( ishade )
         CALL GSPMCI(ishdclr)
         call GSPLCI(ishdclr)
        ELSEIF (ishdclr .eq. 0 .and. iotyp .ge. 20 ) THEN
         call GSPLCI(99)
        ENDIF
        
       CALL GSMKSC(1.0)   ! scale size of marker  
c
c set width of line
c
       CALL GSLWSC(cll)
              
           IF ( iotyp .ge. 20 ) THEN
             CALL NGSETI('CA',2)  ! sets square caps on end of line
           ENDIF

c       write(0,*) 'loop around points'
       DO ix = 1,n1
        DO kz = 1,n2
          IF ( plt2dp(ix,kz) .gt. 0 ) THEN
          DO i = 1, Min(5,Nint(plt2dp(ix,kz)))
           jy = 1
c           x = Float(ix-1)/Float(n1)  ! too far left
c           x = Float(ix)/Float(n1)  ! slightly off to right
           dx1 = 0.75*0.5/n1
!           dy1 = 0.5*hgt/n2
           dy1 = 0.75*0.5*hgt/float(n1)*(aye - ayb)/(axe - axb)
c If there are multiple inits at a point, successive ones will be offset to the
c left and right to keep them visible.
          x = Float(ix-1)/Float(n1-1) + Int((i)/2)*((-1)**i)*0.2/(n1) ! Right on!
          y = 0.001*gzt(kz,1) + hgtlow ! fixme kz
c           write(0,*) 'ix,kz= ',ix,kz
c            CALL GPM(jy,x,y)
            
           IF ( idx .eq. 913 .or. idx == 329 ) THEN  ! neg CG, make a cross

! make a triangle
C            triang(1,1) = x
C            triang(1,2) = y+dy1
C            triang(2,1) = x-dx1
C            triang(2,2) = y-dy1
C
C            triang(3,1) = x+dx1
C            triang(3,2) = y-dy1
C            triang(4,1) = x
C            triang(4,2) = y+dy1
C             
C            CALL CURVE(triang(1,1),triang(1,2),4)


            cross1(1,1) = x-dx1
            cross1(1,2) = y-dy1
            cross1(2,1) = x+dx1
            cross1(2,2) = y+dy1

            cross2(1,1) = x-dx1
            cross2(1,2) = y+dy1
            cross2(2,1) = x+dx1
            cross2(2,2) = y-dy1
             
            CALL CURVE(cross1(1,1),cross1(1,2),2)
            CALL CURVE(cross2(1,1),cross2(1,2),2)
           ELSEIF ( idx .eq. 914 ) THEN   ! pos CG, make a plus

! make a cross:
C            cross1(1,1) = x-dx1
C            cross1(1,2) = y-dy1
C            cross1(2,1) = x+dx1
C            cross1(2,2) = y+dy1
C
C            cross2(1,1) = x-dx1
C            cross2(1,2) = y+dy1
C            cross2(2,1) = x+dx1
C            cross2(2,2) = y-dy1
C
C            CALL CURVE(cross1(1,1),cross1(1,2),2)
C            CALL CURVE(cross2(1,1),cross2(1,2),2)

            cross1(1,1) = x-dx1
            cross1(2,1) = x+dx1
            cross1(1,2) = y
            cross1(2,2) = y

            cross2(1,1) = x
            cross2(2,1) = x
            cross2(1,2) = y+dy1
            cross2(2,2) = y-dy1
             
            CALL CURVE(cross1(1,1),cross1(1,2),2)
            CALL CURVE(cross2(1,1),cross2(1,2),2)
           ELSEIF ( idx .eq. 912 ) THEN ! IC -> diamond
!            dx1 = 0.5*(dx1+dy1)
!            dy1 = 2.*dx1*n2/float(n1)

! make a plus:
C            cross1(1,1) = x-dx1
C            cross1(2,1) = x+dx1
C            cross1(1,2) = y
C            cross1(2,2) = y
C
C            cross2(1,1) = x
C            cross2(2,1) = x
C            cross2(1,2) = y+dy1
C            cross2(2,2) = y-dy1
C
C            CALL CURVE(cross1(1,1),cross1(1,2),2)
C            CALL CURVE(cross2(1,1),cross2(1,2),2)

c            triang(1,1) = x
c            triang(1,2) = y+dy1
c            triang(2,1) = x-dx1
c            triang(2,2) = y-dy1

c            triang(3,1) = x+dx1
c            triang(3,2) = y-dy1
c            triang(4,1) = x
c            triang(4,2) = y+dy1
             
c            CALL CURVE(triang(1,1),triang(1,2),4)

            diamond(1,1) = x
            diamond(1,2) = y+dy1
            diamond(2,1) = x+dx1
            diamond(2,2) = y

            diamond(3,1) = x
            diamond(3,2) = y-dy1
            diamond(4,1) = x-dx1
            diamond(4,2) = y

            diamond(5,1) = diamond(1,1)
            diamond(5,2) = diamond(1,2)
            CALL CURVE(diamond(1,1),diamond(1,2),5)

           ENDIF
            
          ENDDO
          ENDIF
        ENDDO
       ENDDO
c       write(0,*) 'end loop around points'
           IF ( iotyp .ge. 20 ) THEN
             CALL NGSETI('CA',1)  ! sets round caps on end of line
           ENDIF

c
c reset color to black (default line color) and width to 1.
c
       call GSPLCI(1)
       CALL GSLWSC(1.0)


      ELSEIF ( icgsymbol .eq. 2 ) THEN ! } { horizontal overlay
      icgsymbol = 0
        j = 5
        IF ( idx .eq. 325 ) THEN  ! neg CGs
c          write(0,*) 'set cross marker'
          CALL GSMK( 5 )
        ELSEIF ( idx .eq. 326 .or. idx .eq. 329 ) THEN  ! init points
          write(0,*) 'set flash INIT cross marker'
          CALL GSMK( 5 )
          j = 1
        ELSEIF ( idx .eq. 324 ) THEN ! pos CGs
c          write(0,*) 'set plus marker'
          CALL GSMK( 2 )
        ENDIF

         ishdclr = Abs( ishade )
         CALL GSPMCI(ishdclr)
         call GSPLCI(ishdclr)

       CALL GSMKSC(10.0)   ! scale size of marker  
c
c set width of line
c
       CALL GSLWSC(cll)
              
           IF ( iotyp .ge. 20 ) THEN
             CALL NGSETI('CA',2)  ! sets square caps on end of line
           ENDIF

c       write(0,*) 'loop around points'
       DO ix = 1,n1
        DO jy = 1,n2
          IF ( plt2dp(ix,jy) .gt. 0 ) THEN
          DO i = 1, Min(j,Nint(plt2dp(ix,jy)))
c           x = Float(ix-1)/Float(n1)  ! too far left
c           x = Float(ix)/Float(n1)  ! slightly off to right
           dx1 = 1.0/n1
           dy1 = 1.0/n2
          x = Float(ix-1)/Float(n1-1) + 
     :         Int((i+2)/4)*((-1)**i)*0.2/(n1)
          y = Float(jy-1)/Float(n2-1) + 
     :         Int((i+2)/4)*((-1)**Int((i-1)/2) )*0.2/(n2) 
c           x = ix-1
c           y = jy-1
c           x = 0.5
c           y = 0.5

c           write(0,*) 'idx,ix,jy,n1,n2,num= ',
c     :                 idx,ix,jy,n1,n2,plt2dp(ix,jy)
c            CALL GPM(jy,x,y)
            
           IF ( idx .eq. 325 .or. idx .eq. 328 ) THEN ! cross for negative
            cross1(1,1) = x-dx1
            cross1(2,1) = x+dx1
            cross1(1,2) = y-dy1
            cross1(2,2) = y+dy1

            cross2(1,1) = x-dx1
            cross2(2,1) = x+dx1
            cross2(1,2) = y+dy1
            cross2(2,2) = y-dy1
             
            CALL CURVE(cross1(1,1),cross1(1,2),2)
            CALL CURVE(cross2(1,1),cross2(1,2),2)
           ELSEIF ( idx .eq. 326 .or. idx .eq. 329 ) THEN ! cross for INIT
            cross1(1,1) = x-dx1
            cross1(2,1) = x+dx1
            cross1(1,2) = y-dy1
            cross1(2,2) = y+dy1

            cross2(1,1) = x-dx1
            cross2(2,1) = x+dx1
            cross2(1,2) = y+dy1
            cross2(2,2) = y-dy1
             
            CALL CURVE(cross1(1,1),cross1(1,2),2)
            CALL CURVE(cross2(1,1),cross2(1,2),2)
           ELSEIF ( idx .eq. 324 .or. idx .eq. 327 ) THEN ! plus for positive
            cross1(1,1) = x-dx1
            cross1(2,1) = x+dx1
            cross1(1,2) = y
            cross1(2,2) = y

            cross2(1,1) = x
            cross2(2,1) = x
            cross2(1,2) = y+dy1
            cross2(2,2) = y-dy1
             
            CALL CURVE(cross1(1,1),cross1(1,2),2)
            CALL CURVE(cross2(1,1),cross2(1,2),2)
           ENDIF
            
          ENDDO
          ENDIF
        ENDDO
       ENDDO

           IF ( iotyp .ge. 20 ) THEN
             CALL NGSETI('CA',1)   ! sets round caps on end of line
           ENDIF
c       write(0,*) 'end loop around points'

c
c reset color to black (default line color) and width to 1.
c
       call GSPLCI(1)
       CALL GSLWSC(1.0)
        
        ! }
      elseif ( .not. (idx .eq. 23 .or. idx .eq. 32) .and.
     :     ipltyp .ne. 11  ) then
c
c  mult plt2d by xmult
c
      tmp = 0.
      tmp2 = 0.
      do 20101 ix = 1,n1
      do 20102 kz = 1,n2
      
      tmp = Max(tmp,plt2dp(ix,kz))
      tmp2 = Min(tmp2,plt2dn(ix,kz))
      
      
      plt2dp(ix,kz) = plt2dp(ix,kz)*xmult
      plt2dn(ix,kz) = plt2dn(ix,kz)*xmult
      plt2da(ix,kz) = plt2dn(ix,kz)+plt2dp(ix,kz)
20102 continue
20101 continue
c
      write(*,*) 'idebg=16, xmult,min,max=',xmult,tmp2,tmp
c
c  GCA TYPE PLOT (Grid cell array: pseudo raster)
c
      if ( ipware .eq. 2) then
      write(*,*) 'idebg=16a, GCA Plot',idx,icrtyp,nclr,iform
      iform = 1
      ishft = 0
      icrtyp = 0
      call rastset(idx,icrtyp,nclr,iform,rfac1,rfac2,rfac3,rfac4
     >  ,plt2da,nmax,nmax,irast2d,ncbar,labcbar,ifill
     >  ,rlo,rhi,rin,idolab,ilablab,n1,n2) 
      write(*,*) 'idebg=16b, GCA Plot'
      write(*,*) 'idebg=16c, GCA Plot',nmax,n1,n2
       imax = 0
!       imn = 0
       DO j = 1,n2
        DO i = 1,n1
         IF ( irast2d(i,j) .ne. 0 ) THEN
!         write(*,*) 'i,j,irast2d = ',i,j,irast2d(i,j)
         imax = Max(imax,irast2d(i,j))
         ENDIF
        ENDDO
       ENDDO
       write(6,*) 'Max value in irast2d is ',imax
   !   nclr = imax
      call dfclrs2(nclr,icrtyp,ishft)
       
!      call gca(0.,0.,1.,1.,nmax,nmax,1,1,n1,n2,irast2d)
!      call gca(0.,0.,1.,1.,nmax,nmax,2,2,n1,n2,irast2d)
!      call gca(0.,0.,1.,1.,nmax,nmax,1,1,n1-1,n2-1,irast2d)
! GCA fills grid cells, but normal corner point is at the center point of contour plots, so they
! are essentially 1/2 cell offset. To compensate start the GCA plot at (-1/2,-1/2) offset to match
! an contour plot overlay 
      call gca(-0.5/n1,-0.5/n2,(float(n1)+0.5)/n1,(float(n2)+0.5)/n2,nmax,nmax,1,1,n1,n2,irast2d)
      write(*,*) 'idebg=16d, GCA Plot'
      call dfclrs2(nclr,icrtyp,ishft)
      write(*,*) 'idebg=16e, GCA Plot'
      call gsfais(1)
      call sfseti('TY',0)
!      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,BLURED)
      if ( idolab .eq. 1  .and. nclr > 2 ) then
      ifll0 = 1 !  228 - Nint(rfac4) + 1
      if ( ilabor .eq. 1 ) then
      call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >  nclr,0.3,1.0,ifill,1,labcbar,nclr+1,2)  
      end if
      if ( ilabor .eq. 0 ) then
      call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >  nclr,1.0,0.3,ifill,1,labcbar,nclr+1,1)  
      end if
      ifll0 = 0
      end if
      if ( ipltyp .ne. 4 ) then
      if ( ioverpnt .gt. 0 ) then
      if ( ipmore .gt. 0 ) then
      if ( mod(n1-1,i1) .eq. 0 
     >          .and. 
     <     mod(n1-1,i2) .eq. 0 ) then
      mjrx = n1/i1
      mnrx = i1 
      mjry = n2/i2
      mnry = i2
      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
      IF ( jymdiv .ne. 0 ) mnry = jymdiv
      call tick4(itckji,itckni,itckji,itckni)
      call perim(mjrx,mnrx,mjry,mnry)
      call tick4(-itckjo,-itckno,-itckjo,-itckno)
      call perim(mjrx,mnrx,mjry,mnry)
      else
      call tick4(itckji,itckni,itckji,itckni)
      call perim((n1-1),0,(n2-1),0)
      end if 
      end if
      if ( ipmore .eq. 0 ) then
      call perim((1),0,(1),0)
      end if
      end if
      end if
      end if
c
c  OLD CONREC TYPE PLOTS
c
      write(*,*) 'idebg=17'
c
      if ( ipware .eq. 1) then
      if ( icolor .eq. 1 ) then
c
c  plot positive numbers red
c
      isign = 1
      call dfclr(idx,01,irecmj,irecmn,irectx,ndot)
      call conrec(plt2dp,nmax,n1,n2,
     >              rclop,rchip,rcinp,nsetc,nhic,ndot)
c
c  plot negative numbers blue
c
      isign = -1
      call dfclr(idx,-1,irecmj,irecmn,irectx,ndot)
      call conrec(plt2dn,nmax,n1,n2,
     >              rclon,rchin,rcinn,nsetc,nhic,ndot)
      else
c
c  no color  
c 
      ndot = 1
      if ( icdash .eq. 1 ) ndot = -682
      call conrec(plt2da,nmax,n1,n2,
     >              rclon,rchip,rcinp,nsetc,nhic,ndot)
      endif
c
      endif
c
c  CONPAC SCHEME FOR PLOTS 
c
c
      write(*,*) 'idebg=18'
c
      if ( ipware .eq. 0 ) then
      print*,'USING CONPAC'
c
c
c  SET SOME PARAMETERS
C
C
      nclr = max(nclr,4)
      call DFCLRS(nclr,icrtyp,ishft)
      call gsmksc(GKSDOT)
c
c
c
c  If drawing 'trajectory' then call curve
c
c
      if (idx.eq.997.or.idx.eq.998.or.idx.eq.999) then
      read(22,*) numtraj
      do num = 1,numtraj
      read(22,*) filnmtraj,width,ilabcrv,ilabcrvtyp
      open(29,file=filnmtraj,form='formatted')
      read(29,*) nt
      do it = 1,nt
      read(29,*) ii,xplt(it),yplt(it),zplt(it)
      zplt(it) = zplt(it)
      xplt(it) = xplt(it)/(max(1,(nx-2))*dx)
      yplt(it) = yplt(it)/(max(1,(ny-2))*dy)
      zplt(it) = zplt(it)/(max(1,(nz-2))*dz)
      if ( ilabcrv .ge. 1 ) then
      if ( it/ilabcrv*ilabcrv .eq. it .or. it .eq. 1 ) then
      if ( ilabcrvtyp.eq.9 ) labcrv = "#"
      if ( ilabcrvtyp.eq.9 ) labcrv = "*"
      if ( ilabcrvtyp.eq.9 ) labcrv = "+"
      if ( ilabcrvtyp.eq.0 ) labcrv = "o"
      if ( ilabcrvtyp.eq.1 ) labcrv = "+"
      if ( ilabcrvtyp.eq.-1 ) labcrv = "-"
      if ( idx .eq. 997) then
      call pwrit(xplt(it),(zplt(it)+.01),labcrv,1,2,0,0)
      end if
      if ( idx .eq. 998) then
      call pwrit(yplt(it),(zplt(it)+.01),labcrv,1,2,0,0)
      end if
      if ( idx .eq. 999) then
      call pwrit(xplt(it),yplt(it),labcrv,1,2,0,0)
      end if
      end if
      end if
      end do
      close(29)
      call gslwsc(width)
      if (idx.eq.997) call curve(xplt,zplt,nt)
      if (idx.eq.998) call curve(yplt,zplt,nt)
      if (idx.eq.999) call curve(xplt,yplt,nt)
      end do
      end if
c
c
      if ( rcinp .ne. 0.0 ) then
      nclv = Int((rchip-rclon)/rcinp +0.001) + 1
      end if
      if ( rchip .eq. rclon ) nclv = 1
      
      IF ( icrfil == 3 ) THEN ! { custom levels and colors
        nclv = 0
        DO i = 1,maxclr
         IF ( customlevels(i) /= -999.0 ) THEN
          nclv = nclv + 1
          write(0,*) 'i, customlevel = ',i,customlevels(i),customcolors(i)
         ELSE
          exit
         ENDIF
        ENDDO
        write(0,*) 'number of custom levels = ', nclv
      
      ENDIF ! }
      
      IF ( xmult .lt. 0.0 ) THEN
        tmp = xmx
        xmx = xmult*xmn
        xmn = xmult*tmp
        write(*,*) 'idebg: xmult<0; xmn,xmx = ',xmn,xmx
      ENDIF
      
      IF ( icrfil .eq. 2 .or. icrfil .eq. 4 ) THEN !{
       
       IF ( rcinp .eq. 0  ) THEN ! { auto number of contour levels
        IF ( rclon .eq. 0.0 .or. xmn .ge. 0.0 ) THEN 
          nclv = 10
          posdef = 1
        ELSE
          nclv = 19
          posdef = 0
        ENDIF
       CALL CPSETI ('CLS - CONTOUR LEVEL SELECTION FLAG', nclv/2)
       CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER', 0.)
       CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
       CALL cppkcl(plt2da,rwrk,iwrk)
       CALL CPGETR ('CMN - CONTOUR MINIMUM', rclon)
       CALL CPGETR ('CMX - CONTOUR MAXIMUM', newval)
       CALL CPGETR ('CIS - CONTOUR INTERVAL SPECIFIER', rcinp)
c       write(*,*) 'CMN,CMX,CIS = ',rclon,newval,rcinp
       CALL CPGETI ('NCL - NUMBER OF CONTOUR LINES',NOCL)
c       write(*,*) 'NOCL = ',nocl
       rclon = 0
       rchip = 0
       DO i = 1,nocl
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
          CALL CPGETR ('CLV - CONTOUR LEVEL VALUE',CLEV)
          IF ( i .eq. 1 ) rcinp = clev
          IF ( i .eq. 2 ) rcinp = clev - rcinp
          rclon = Min( rclon, clev )
          rchip = Max( rchip, clev )
c          write(*,*) 'Contour number, clev = ',i,clev
       ENDDO
       newval = Max( Abs(rclon), Abs(rchip) )
c       rclon = -newval
       IF ( posdef .eq. 1 ) rclon = 0
       rclop = rclon
       rchip = newval
c       write(*,*) 'rclon,rchip,rcinp = ',rclon,rchip,rcinp
       nclv = nocl
       
       ELSE ! } {
         maxval = Max( Abs(xmn), Abs(xmx) )
         IF ( rchip .ne. 0.0 ) maxval = rchip
         newval = rcinp*Int(maxval/rcinp)
        IF ( xmn .ge. 0.0 .or. rlo .eq. 0.0 ) THEN 
          rclon = 0.0
          posdef = 1
        ELSE
          rclon = -newval
          nclv = 19
          posdef = 0
        ENDIF
          nclv = newval*(posdef + 1)/rcinp
       CALL CPSETI ('CLS - CONTOUR LEVEL SELECTION FLAG', nclv)
       CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER', rcinp)
       CALL CPSETI ('LIS - LABEL INTERVAL SPECIFIER', 2)
       CALL CPSETR ('CMN - CONTOUR MINIMUM', rclon)
       CALL CPSETR ('CMX - CONTOUR MAXIMUM', newval)
       CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
       CALL cppkcl(plt2da,rwrk,iwrk)
c       CALL CPGETR ('CMN - CONTOUR MINIMUM', rclon)
c       CALL CPGETR ('CMX - CONTOUR MAXIMUM', newval)
c       CALL CPGETR ('CIS - CONTOUR INTERVAL SPECIFIER', rcinp)
c       write(*,*) 'CMN,CMX,CIS = ',rclon,newval,rcinp
       CALL CPGETI ('NCL - NUMBER OF CONTOUR LINES',NOCL)
       write(*,*) 'NOCL = ',nocl
       DO i = 1,nocl
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
          CALL CPGETR ('CLV - CONTOUR LEVEL VALUE',CLEV)
c          write(*,*) 'Contour number, clev = ',i,clev
       ENDDO
       nclv = nocl
       ENDIF ! }
       
      
      
      
      ELSEIF ( icrfil == 3 ) THEN ! } {
C
       rclon = customlevels(1)
      CALL CPSETC ('CFT - CONSTANT-FIELD LABEL TEXT STRING',' ')
      CALL CPSETI ('CLS - CONTOUR LEVEL SELECTION FLAG',0)
      CALL CPSETI ('NCL - NUMBER OF CONTOUR LEVELS',nclv)

C       CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER', 0.)
c       CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
c       CALL cppkcl(plt2da,rwrk,iwrk)
c       CALL CPGETR ('CMN - CONTOUR MINIMUM', rclon)
C       CALL CPGETR ('CMX - CONTOUR MAXIMUM', newval)
C       CALL CPGETR ('CIS - CONTOUR INTERVAL SPECIFIER', rcinp)
Cc       write(*,*) 'CMN,CMX,CIS = ',rclon,newval,rcinp
C       CALL CPGETI ('NCL - NUMBER OF CONTOUR LINES',NOCL)
Cc       write(*,*) 'NOCL = ',nocl
C       rclon = 0
C       rchip = 0
C       DO i = 1,nocl
C          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',I)
C          CALL CPGETR ('CLV - CONTOUR LEVEL VALUE',CLEV)
C          IF ( i .eq. 1 ) rcinp = clev
C          IF ( i .eq. 2 ) rcinp = clev - rcinp
C          rclon = Min( rclon, clev )
C          rchip = Max( rchip, clev )
Cc          write(*,*) 'Contour number, clev = ',i,clev
C       ENDDO
C       newval = Max( Abs(rclon), Abs(rchip) )
Cc       rclon = -newval
C       IF ( posdef .eq. 1 ) rclon = 0
C       rclop = rclon
C       rchip = newval
Cc       write(*,*) 'rclon,rchip,rcinp = ',rclon,rchip,rcinp
C       nclv = nocl

      
      ELSE ! }
      
      print*,'nclv',nclv
      CALL CPSETC ('CFT - CONSTANT-FIELD LABEL TEXT STRING',' ')
      CALL CPSETI ('CLS - CONTOUR LEVEL SELECTION FLAG',0)
      CALL CPSETI ('NCL - NUMBER OF CONTOUR LEVELS',nclv)

      ENDIF

      CALL CPSETR ('CWM - CHARACTER-WIDTH MULTIPLIER', CWM)

C
C Make PLOTCHAR use high-quality characters.
C
      CALL PCSETI ('QU - QUALITY FLAG',IPQUAL)
c
c  Informational Label text
c
      if ( ioffm .eq. 1 ) then
      call cpsetc('ILT - INFORMATIONAL LABEL TEXT',' ')
      end if
C
C Turn on the positioning of labels by the penalty scheme. (LLP)
C
      if ( ilab .ge. 1 ) then
c      CALL CPSETI ('LLP - LINE LABEL POSITIONING',2)
c       CALL CPSETI ('LLP - LINE LABEL POSITIONING', Max(2,ilab))
       CALL CPSETI ('LLP - LINE LABEL POSITIONING', ilab)
       CALL CPSETI ('LLO - LINE LABEL ORIENTATION', 1)
      else
      CALL CPSETI ('LLP - LINE LABEL POSITIONING',0)
      endif
      
      IF ( labelsize .ne. 0.0 ) THEN
        CALL CPSETR('LLS - LINE LABEL SIZE',labelsize)
       ENDIF
C
C Turn on the drawing of the high and low label boxes. (HLB)
C Tell CONPACK to delete high/low labels which overlap the informational
C label or another high/low label, but to move those which overlap the
C edge inward a little. (HLO)
C
      if ( nhic .ge. 1 ) then
        if ( icolor .eq. 1 ) then
           call CPSETI('HLC - HIGH/LOW LABEL COLOR INDEX', 1)  
        end if
        if ( nhic .ge. 2 ) then
          CALL CPSETI ('HLB - HIGH/LOW LABEL BOX FLAG',1)
         end if
      else
        CALL CPSETI ('HLO - HIGH/LOW LABEL OVERLAP FLAG',7)
        CALL CPSETC ('HLT - HIGH/LOW LABEL TEXT',' ')
      end if  !  ( nhic .ge. 1 )
c
c  SET STUFF
c
      write(*,*) 'idebg=19; NCLv = ',NCLV
      
      labcbar(:) = '   '
c
      IF (  icrfil .le. 20 ) THEN
       DO ICLV=1,NCLV ! {
c
c  Set parameter array indedex (PIA)
c
       CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',ICLV)
c
c  Set countour levels
c
      clev = (iclv-1)*rcinp + rclop
      IF ( icrfil == 3 ) THEN
        clev = customlevels(iclv)
        ishade = customcolors(iclv) + 100
      ENDIF
      
      CALL CPSETR ('CLV - CONTOUR LEVEL VALUE',CLEV)
c
c  set labeling of line labels
c
      if ( ilab .ge. 1 ) then
c        IF(MOD(ICLV,ILAB).EQ.1.AND.ILAB.NE.1) THEN
        IF ( mcl .le. -1. ) THEN
          IF(MOD(ICLV,Abs(Nint(mcl))).EQ.0) THEN
           IF ( cll .gt. 0.0 ) THEN
             CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',iclu)
             CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',cll)
           ENDIF
          if ( cll .le. 0.0 ) then
            CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',2) ! labels w/o contour lines
          end if
          ENDIF
        ELSEIF(MOD(ICLV,ILAB).EQ.0) THEN
           IF ( cll .gt. 0.0 ) THEN
             CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',iclu)
             CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',cll)
           ENDIF
          if ( cll .le. 0.0 ) then
            CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',2) ! labels w/o contour lines
          end if
        ELSE
           IF ( cll .gt. 0.0 ) THEN
              CALL CPSETI ('CLU - CONTOUR LEVEL USE',1)
              CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',cll)
           ELSE
              CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',1)
           ENDIF
       END IF
c      CALL CPGETI ('CLU - CONTOUR LEVEL USE FLAG',iclu)
c        IF( CLEV .EQ. 0) THEN
c          CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',1)
c        END IF
        if ( NCLV .EQ. 1 ) THEN
           CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',cll)
c           CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',1)
        END IF
      elseif ( cll .le. 0.0 ) THEN
          CALL CPSETI ('CLU - CONTOUR LEVEL USE',0)  ! no contour line if ilab < 1 and cll=0
      ELSE ! ( ilab .lt. 1 )
          CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',cll)
          CALL CPSETI ('CLU - CONTOUR LEVEL USE',1)  ! no labels if ilab < 1
      end if
c
c  make contours
c
      if ( mcl .eq. 0 ) then
      CALL CPSETI ('CLU - CONTOUR LEVEL USE FLAG',0)
      end if
c
c  color
c
      if ( ioverlay .eq. -17 ) then 
      CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',2.0)
      end if
      if ( icolor .eq. 1 ) then
      if ( clev .eq. 0 ) then
      CALL CPSETI ('CLC - CONTOUR LINE COLOR',1)
      CALL CPSETI ('LLC - LINE LABEL COLOR',1)
      CALL CPSETR ('CLL - CONTOUR-LINE LINE WIDTH',1.5)
      end if
      if ( clev .le. 0 ) then
      call dfclr(idx,-1,irecmj,irecmn,irectx,ndot)
      CALL CPSETI ('CLC - CONTOUR LINE COLOR',irecmj)
      CALL CPSETI ('LLC - LINE LABEL COLOR',irecmj)
      end if
      if ( clev .ge. 0 ) then
      call dfclr(idx,01,irecmj,irecmn,irectx,ndot)
      CALL CPSETI ('CLC - CONTOUR LINE COLOR',irecmj)
      CALL CPSETI ('LLC - LINE LABEL COLOR',irecmj)
      end if
      end if
c
c Dash
c
      if ( clev .lt. 0. .and. Abs( ishade ) .lt. 100) then
      if ( icdash .eq. 1 ) then
       CALL CPSETI ('CLD - CONTOUR LINE DASH PATTERN', idash)
        CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', inegclr)
      ELSEIF ( icdash .eq. 2 ) THEN
        CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', inegclr)
      ELSEIF ( icdash .eq. 3 ) THEN
       CALL CPSETI ('CLD - CONTOUR LINE DASH PATTERN', idash)
      end if
      
      end if
c
c  Shading
c
      write(*,*) 'idebg=20, clev,ishade = ',clev,ishade
c
      if ( ishade .eq. 0 ) then
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      end if
c
c      write(0,*) 'icrfil,idx= ',icrfil,idx
      if ( icrfil .eq. 1 .and.  
     :  (idx .eq. 22 .or. idx .eq. 35 .or. idx .eq. 991 .or. idx == 601) ) then
      IF (CLEV.GE.RLO) THEN ! {
      
      idum = Int(clev/5.0)
      ifill(iclv) = idum+200
      IF ( clev - Int(clev) .eq. 0.0 ) THEN
        IF ( clev .lt. 1000. .and. clev .gt. -100. ) THEN
        write(labcbar(iclv),'(i3)') Int(clev) 
        ELSE
        write(labcbar(iclv),'(i5)') Int(clev)         
        ENDIF
      ELSE
        write(labcbar(iclv),'(f5.1)') clev 
      ENDIF
c      write(0,*) 'iclv,idum,clev,rlo = ',iclv,idum,clev,rlo
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',idum)
       IF ( iclv .eq. 1 ) THEN
         idum = 0
       ELSE
         idum = ifill(iclv-1) - 200
       ENDIF
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',idum)
      
      ELSE 
      
      ifill(iclv) = idum
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      
      END IF ! }
      
      elseif ( ishade .ge. 1 ) then
      
      IF (CLEV.GT.RLO) THEN
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      ELSE 
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',1)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',2)
      END IF
      
      end if
      
      if ( ishade .le. -1 ) then
      IF (CLEV.LT.RHI) THEN
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      ELSE
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',2)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',1)
      END IF
      end if
      
c
      if ( icrfil .eq. 2 ) then ! {
      IF (CLEV.GE.rclon) THEN ! {
      
      idum = Min(8, Max(-8, Int(clev/rchip*8.0) ) )
      IF ( clev .lt. 0.0 ) ifill(iclv+1) = Min(228,idum+229)
      IF ( clev .ge. 0.0 ) ifill(iclv+1) = idum+229
c      write(*,*) 'clev, idum = ',clev,idum,ifill(iclv+1)
      idum = idum + 10
      IF ( clev - Int(clev) .eq. 0.0 ) THEN
        IF ( clev .lt. 1000. .and. clev .gt. -100. ) THEN
        write(labcbar(iclv+1),'(i3)') Int(clev) 
        write(labcbar(iclv+2),'(a)') '   '
        ELSE
        write(labcbar(iclv+1),'(i6)') Int(clev) 
        write(labcbar(iclv+2),'(a)') '   '
        ENDIF
      ELSE
        write(labcbar(iclv+1),'(f5.1)') clev 
      ENDIF
c      write(0,*) 'iclv,idum,clev,ifill = ',iclv,idum,clev,ifill(iclv+1),Max(0,ifill(iclv)),
c     :       Max(0,ifill(iclv) - 220)
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',idum)
       IF ( iclv .eq. 1 ) THEN
c         idum = 1
         idum = Max( 1, idum - 1 )
         ifill(iclv) = ifill(iclv+1)-1
       ELSE
         idum = ifill(iclv) - 219
       ENDIF
c      write(0,*) 'aib = ',idum
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',idum)
      
      ELSE 
      
      ifill(iclv) = 229 ! idum
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      
      END IF !}
      ENDIF !}
      
      if ( icrfil .eq. 3 ) then !{
      IF (CLEV.GE.rclon) THEN ! {
      
      idum =  iclv ! Min(8, Max(-8, Int(clev/rchip*8.0) ) )
      ifill(iclv) = customcolors(iclv)
      ifill(iclv+1) = customcolors(iclv+1)
!      IF ( clev .lt. 0.0 ) ifill(iclv+1) = Min(228,idum+229)
!      IF ( clev .ge. 0.0 ) ifill(iclv+1) = idum+229
      write(*,*) 'clev, idum = ',clev,idum,ifill(iclv+1)
!      idum = idum + 10
      IF ( clev - Int(clev) .eq. 0.0 ) THEN !{
        IF ( clev .lt. 1000. .and. clev .gt. -100. ) THEN
        write(labcbar(iclv),'(i3)') Int(clev) 
        write(labcbar(iclv+1),'(a)') '   '
        ELSE
        write(labcbar(iclv),'(i6)') Int(clev) 
        write(labcbar(iclv+1),'(a)') '   '
        ENDIF
      ELSE
        IF ( idx == 607 ) THEN
        write(labcbar(iclv),'(f4.2)') clev 
        ELSE
        write(labcbar(iclv),'(f5.1)') clev 
        ENDIF
        write(labcbar(iclv+1),'(a)') '   '
      ENDIF ! }
c      write(0,*) 'iclv,idum,clev,ifill = ',iclv,idum,clev,ifill(iclv+1),Max(0,ifill(iclv)),
c     :       Max(0,ifill(iclv) - 220)
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',idum)
       IF ( iclv .eq. 1 ) THEN
!         idum = 0
!         idum = Max( 1, idum - 1 )
!         ifill(iclv) = ifill(iclv+1)-1
       ELSE
!         idum = ifill(iclv) - 219
       ENDIF
       idum = idum - 1
c      write(0,*) 'aib = ',idum
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',idum)
      
      ELSE !}
      
      ifill(iclv) = 0 ! 229 ! idum
      CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
      CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
      
      END IF 
      ENDIF !} icrfil .eq. 3
      
            
      IF ( Abs( ishade ) .ge. 20 .and. Abs( ishade ) .lt. 100 ) THEN
         ishdclr = Abs( ishade )
        CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', ishdclr)
      ELSEIF ( ishdclr .eq. 0 .and. iotyp .ge. 20
     :         .and. clev .ge. 0.) THEN
        ishdclr = 99
        CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', ishdclr)
        ishdclr = 0
c         call GSPLCI(99)
      END IF

        IF ( ishade .eq. 200 .and. 
     :  (idx .eq. 22 .or. idx .eq. 35 .or. idx .eq. 991 .or. idx == 601) ) THEN
        ishdclr = 200+Int(clev/5.0)
        CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', ishdclr)
        ENDIF

      end do ! } iclv
      
      ENDIF ! icrfil .le. 20


        IF ( icrfil .eq. 1 .and.  
     :  (idx .eq. 22 .or. idx .eq. 35 .or. idx .eq. 991 .or. idx == 601) ) THEN

c gsfais may or may not be needed....
         CALL GSFAIS (1)
         CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
         CALL ARINAM (IAMA,LIAMA)

         CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
         CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,COLDBZ)
c
c Add contour lines if cll > 0.0 
c
         IF ( cll .gt. 0.0 ) THEN
           CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
         ENDIF

         nclr = nclv
         if ( idolab .eq. 1  .and. nclr > 2 ) then
          if ( ilabor .eq. 1 ) then
           call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >       nclr,0.3,1.0,ifill,1,labcbar,nclr+1,2)
          end if
          if ( ilabor .eq. 0 ) then
            call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >          nclr,1.0,0.3,ifill,1,labcbar,nclr+1,1)
          end if
         end if

        ENDIF
        
        IF ( icrfil .eq. 2 ) THEN

c gsfais may or may not be needed....
         CALL GSFAIS (1)
!         CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
         CALL ARINAM (IAMA,LIAMA)

         CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
c         IF ( posdef .eq. 0 ) THEN
           CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,BLURED)
c         ELSE
c           CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,REDPOS)
c         ENDIF

c
c Add contour lines if cll > 0.0 
c
         IF ( cll .gt. 0.0 ) THEN
           CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
         ENDIF

         nclr = nclv+1-posdef
c         IF ( posdef .eq. 0 ) THEN
c           DO i = 1,Min(18,nclr)
c           ifill(i) = 220 + i
c           ENDDO
c         ELSE
c           DO i = 1,Min(9,nclr)
c           ifill(i) = 229 + i
c           ENDDO
c         ENDIF
         if ( idolab .eq. 1  .and. nclr > 2 ) then
          if ( ilabor .eq. 1 ) then
           call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >       nclr,0.3,1.0,ifill,1,labcbar,nclr+1,1)
           write(*,*) 'DEBUG20a: nclr = ',nclr
          end if
          if ( ilabor .eq. 0 ) then
            call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >          nclr,1.0,0.3,ifill,1,labcbar,nclr+1,1)
            write(*,*) 'DEBUG20b: nclr = ',nclr
          end if
         end if

        ENDIF
       

        IF ( icrfil .eq. 3 ) THEN

c gsfais may or may not be needed....
         CALL GSFAIS (1)


         CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
         CALL ARINAM (IAMA,LIAMA)

         CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
         CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,CSTCLR)

c
c Add contour lines if cll > 0.0 
c
         IF ( cll .gt. 0.0 ) THEN
           CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
         ENDIF

         nclr = nclv
!         nclr = nclv+1-posdef
c         IF ( posdef .eq. 0 ) THEN
c           DO i = 1,Min(18,nclr)
c           ifill(i) = 220 + i
c           ENDDO
c         ELSE
c           DO i = 1,Min(9,nclr)
c           ifill(i) = 229 + i
c           ENDDO
c         ENDIF
         if ( idolab .eq. 1 .and. nclr > 2 ) then
          if ( ilabor .eq. 1 ) then
           write(*,*) 'DEBUG20a3a: nclr = ',nclr
           call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >       nclr,0.3,1.0,ifill,1,labcbar,nclr+1,1)
           write(*,*) 'DEBUG20a3: nclr = ',nclr
          end if
          if ( ilabor .eq. 0 ) then
            call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >          nclr,1.0,0.3,ifill,1,labcbar,nclr+1,1)
            write(*,*) 'DEBUG20b3: nclr = ',nclr
          end if
         end if

        ENDIF

        IF ( icrfil .eq. 4 ) THEN

c gsfais may or may not be needed....
         CALL GSFAIS (1)
!         CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)
         CALL ARINAM (IAMA,LIAMA)

         CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
c         IF ( posdef .eq. 0 ) THEN
           CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,BLUNEG)
c         ELSE
c           CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,REDPOS)
c         ENDIF

c
c Add contour lines if cll > 0.0 
c
         IF ( cll .gt. 0.0 ) THEN
           CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
         ENDIF

         nclr = nclv+1-posdef
c         IF ( posdef .eq. 0 ) THEN
c           DO i = 1,Min(18,nclr)
c           ifill(i) = 220 + i
c           ENDDO
c         ELSE
c           DO i = 1,Min(9,nclr)
c           ifill(i) = 229 + i
c           ENDDO
c         ENDIF
         if ( idolab .eq. 1 ) then
          if ( ilabor .eq. 1 ) then
           call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >       nclr,0.3,1.0,ifill,1,labcbar,nclr+1,1)
           write(*,*) 'DEBUG20a: nclr = ',nclr
          end if
          if ( ilabor .eq. 0 ) then
            call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >          nclr,1.0,0.3,ifill,1,labcbar,nclr+1,1)
            write(*,*) 'DEBUG20b: nclr = ',nclr
          end if
         end if

        ENDIF
       

c
c
      write(*,*) 'idebg=21'
c
c
      if ( icrfil .eq. 0  )  then

c ARINAM - Initializes an area map
      CALL ARINAM (IAMA,LIAMA)
      
      tmp=0.0
      tmp2=0.0
      DO j=1,n2
      DO i=1,n1
        tmp = max(tmp, plt2da(i,j))
        tmp2 = min(tmp2, plt2da(i,j))
      ENDDO
      ENDDO
      
      write(*,*) 'idebug=21a, min,max = ',tmp2,tmp
      
      CALL CPRECT (plt2da,nmax, n1,   n2,   RWRK,LRWK,IWRK,LIWK)



c      CALL GSLWSC (2.)
c     CALL CPBACK (PLT2DA,RWRK,IWRK)
c GSLWSC (Set linewidth scale factor) 
      CALL GSLWSC (1.)
c ARINAM - Initializes an area map
c      CALL ARINAM (IAMA,LIAMA)
c CPLBAM - Adds label boxes
      CALL CPLBAM (plt2da,RWRK,IWRK,IAMA)
      
      IF ( ifll .eq. 0  .or. ishade .eq. 200 ) THEN
c CPCLDM - draw contour lines masked by an existing area map
      CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
c CPLBDR - Draws labels 
      CALL CPLBDR (plt2da,RWRK,IWRK)
      ENDIF
C
C Now, add the zero contour line to the area map.
C
      IF ( ISHADCNT0 .eq. 1 ) THEN
      CALL GSLWSC (2.)
      CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
      CALL GSLWSC (1.)
      ELSE
c      CALL GSLWSC (4.)
      CALL CPCLAM (plt2da,RWRK,IWRK,IAMA)
      CALL GSLWSC (1.)
      END IF
C
C Scan the area map.  The routine SHDRH1 will be called to shade the
C areas with hatching
C
      write(*,*) 'ISHADE IS: ', ishade
      if ( ishade .eq. 01 .or. ishade .eq. -01 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDH1)
      end if
      if ( ishade .eq. 02 .or. ishade .eq. -02 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDH2)
      end if
      if ( ishade .eq. 03 .or. ishade .eq. -03 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDDOT)
      end if
      if ( ishade .eq. 04 .or. ishade .eq. -04 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDFIL)
      end if
C
      if ( ishade .eq. 11 .or. ishade .eq. -11 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDHV)
      end if
      if ( ishade .eq. 12 .or. ishade .eq. -12 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDHH)
      end if
      if ( ishade .eq. 13 .or. ishade .eq. -13 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDHR)
      end if
      if ( ishade .eq. 14 .or. ishade .eq. -14 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDHL)
      end if
      if ( ishade .eq. 15 .or. ishade .eq. -15 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDHX)
      end if
      if ( ishade .eq. 16 .or. ishade .eq. -16 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDD2)
      end if
      if ( ishade .eq. 17 .or. ishade .eq. -17 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDD4)
      end if
      if ( ishade .eq. 18 .or. ishade .eq. -18 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDD6)
      end if
      if ( ishade .eq. 19 .or. ishade .eq. -19 ) then
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDD8)
      end if
!      if ( ishade .eq. 20 .or. ishade .eq. -20 ) then
!      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDFIL)
!      end if

!      IF ( ishade .eq. 21 .or. ishade .eq. -21 ) THEN
!         ishdclr = 51
!         CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDCLR)
!      END IF
      
        IF ( ishade .eq. 201 .and.  
     :   (idx .eq. 22 .or. idx .eq. 35 .or. idx .eq. 991 .or. idx == 601) ) THEN
        write(0,*) 'call arscam with coldbz'
        CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,COLDBZ)
        ENDIF

      IF ( ( Abs( ishade ) .ge. 130 .and. Abs( ishade ) .lt. 200) .or.
     :       Abs(ishade) .gt. 300 ) THEN
         
         IF ( Abs(ishade) .gt. 300 ) THEN 
           ishdclr = Abs( ishade ) - 100
         ELSE
           ishdclr = Abs( ishade ) - 100
         ENDIF
c        IF ( ishade .eq. 200 .and. idx .eq. 22 ) THEN
c         CALL setdbzclr(ishdclr, clev)
c        ENDIF
        
      write(*,*) 'idebg=21z, ishdclr= ',ishdclr
        CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,SHDCLR)
       
       IF ( cll .gt. 0.0 ) THEN
        CALL CPCLDM (plt2da,RWRK,IWRK,IAMA,DRAWCL)
c CPLBDR - Draws labels 
        CALL CPLBDR (plt2da,RWRK,IWRK)
       ENDIF
      END IF
         

      end if ! ( icrfil .eq. 0 )
c
c
c   if icrfil .eq. 1
c
      write(*,*) 'idebg=22'
c
      if ( icrfil .eq. 22 ) then
      iform = -1
      call rastset(idx,icrtyp,nclr,iform,rfac1,rfac2,rfac3,rfac4
     >  ,plt2da,nmax,nmax,irast2d,ncbar,labcbar,ifill
     >  ,rlo,rhi,rin,idolab,ilablab,n1,n2) 
      CALL GSFAIS (1)
      CALL DFCLRS(nclr,icrtyp,ishft)
c     CALL CPSETI('NOF - NUMERIC OMMISSION FLAGS',0)
      CALL CPRECT (plt2da,nmax,n1,n2,RWRK,LRWK,IWRK,LIWK)
      CALL ARINAM (IAMA,20000)
      CALL CPLBAM (plt2da,RWRK,IWRK,IAMA)
      CALL ARSCAM (IAMA,XCRA,YCRA,MCS,IARA,IGRA,MAI,COLRAM)
      CALL GSPLCI (0) 
      CALL CPCLDR (plt2da,RWRK,IWRK)
      CALL GSPLCI (1) 
      CALL LBSETI ('CBL - COLOR OF BOX LINES',0)
      call dfclrs(nclr,icrtyp,ishft)
      call gsfais(1)
      call sfseti('TY',0)
      if ( idolab .eq. 1 ) then
      if ( ilabor .eq. 1 ) then
      call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >  nclr,0.3,1.0,ifill,1,labcbar,nclr,2)
      end if
      if ( ilabor .eq. 0 ) then
      call lblbar(ilabor,xbarl,xbarr,ybarb,ybart,
     >  nclr,1.0,0.3,ifill,1,labcbar,nclr,1)
      end if
      end if
      end if
C
C
      write(*,*) 'idebg=23'
c
C
      end if ! ( icrfil .eq. 1 )
c
c
      write(*,*) 'idebg=24'
c
c
c
      endif ! (idx .ne. 23)?
c
c
      MNRX = i1
      MNRY = i2
      if ( ipltyp .ne. 4 .and. ipltyp .ne. 11 ) then
      if ( ioverpnt2 .eq. 1 ) then
      if ( ipmore .gt. 0 ) then
      axb = 0.15
      ayb = 0.20
      rxy = max(rx,ry)
      axe = axb + rx*(0.55)/rxy
      aye = ayb + ry*(0.55)/rxy
      call set(axb,axe,ayb,aye,0.,1.,0.,1.,1)      
      if ( mod(n1-1,i1) .eq. 0 
     >          .and. 
     <     ( mod(n2-1,i2) .eq. 0 .or. 
     :      (mod(n2set-1,i2).eq.0 .and. n2set.ne.0) ) ) then
      mjrx = n1/i1
      mnrx = i1 
      mjry = n2/i2
      mnry = i2

      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
      IF ( jymdiv .ne. 0 ) mnry = jymdiv
      IF ( (ipltyp .eq. 2 .or. ipltyp .eq. 3 .or. ipltyp == 5) 
     :       .and. kzmdiv .ne. 0 ) mnry = kzmdiv

      IF ( n2set .ne. 0 ) THEN
      mjry = n2set/i2
      ELSE
      mjry = n2/i2
      ENDIF
c      mnry = i2
      write(6,*) 'a: mnrx,mnry,i1,i2 = ', mnrx,mnry,i1,i2
      call tick4(itckji,itckni,itckji,itckni)
      call perim(mjrx,mnrx,mjry,mnry)
      call tick4(-itckjo,-itckno,-itckjo,-itckno)
      call perim(mjrx,mnrx,mjry,mnry)
      ELSEIF ( mod(n1-1,i1) .eq. 0 
     >          .and. 
     <     mod(n2-1,2) .eq. 0 ) then
      mjrx = n1/i1
      mnrx = i1 
      mjry = n2/2
      mnry = 2

      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
c      IF ( jymdiv .ne. 0 ) mnry = jymdiv

      write(6,*) 'b: mnrx,mnry,i1,i2 = ', mnrx,mnry,i1,i2
      call tick4(itckji,itckni,itckji,itckni)
      call perim(mjrx,mnrx,mjry,mnry)
      call tick4(-itckjo,-itckno,-itckjo,-itckno)
      call perim(mjrx,mnrx,mjry,mnry)
      ELSEIF( mod(n1-1,i1) .eq. 0  ) THEN
      mjrx = n1/i1
      mnrx = i1 
      IF ( ixmdiv .ne. 0 ) mnrx = ixmdiv
c      IF ( jymdiv .ne. 0 ) mnry = jymdiv
      write(6,*) 'c: mnrx,mnry,i1,i2 = ', mnrx,mnry,i1,i2
      call tick4(itckji,itckni,itckni,itckji)
      call perim(mjrx,mnrx,(n2-1),0)
      else
      write(6,*) 'd: mnrx,mnry,i1,i2 = ', mnrx,mnry,i1,i2
      call tick4(itckji,itckni,itckji,itckni)
c      call perim(n1/i1,3,n2/i2,3)
c      call perim((n1-1),0,(n2-1),0)
      call perim(2,(n1-3),2,(n2-3))
      end if 
      end if
      if ( ipmore .eq. 0 ) then
      write(6,*) 'e: mnrx,mnry,i1,i2 = ', mnrx,mnry,i1,i2
      call perim((1),0,(1),0)
      end if
      end if
      end if

      if ( idx .eq. 23 .or. idx .eq. 32 ) then
      print*,'ISPV,SPV',ispv,spv
      IF ( ipltyp .eq. 2 .or. ipltyp .eq. 3 .or.
     :     ipltyp .eq. 5 ) THEN
       write(6,*) 'set map transformation'
       dd = nl*dx
c       hgt = 0.001*(gt(1,1,nz,mzdist) + 0.5*dzmax)
c       hgt0 = 0.0 ! gt(1,1,1,mzdist)
       IF ( hgtset .eq. 0.0 ) THEN
        hgt = 0.001*(gzt(n2,1) + 
     :    0.5*(gzt(n2+1,1)-gzt(n2,1)))
        hgt0 = 0.0
       ELSE
        hgt = hgtset
        hgt0 = hgtlow
       ENDIF
       call set(axb,axe,ayb,aye,0.,1.,hgt0,hgt,1)
c       write(0,*) 'axb,axe,ayb,aye',axb,axe,ayb,aye
c       write(0,*) 'hgt0,hgt,nzset,n1,n2 ',hgt0,hgt,nzset,n1,n2
       CALL VVSETI('SET - DO SET CALL FLAG',0)
c       CALL VVSETR('XC1 - X COORDINATE AT INDEX 1',0.)
c       CALL VVSETR('XCM - X COORDINATE AT INDEX M',1.0)

       CALL VVSETR('XC1 - X COORDINATE AT INDEX 1',axb)
       CALL VVSETR('XCM - X COORDINATE AT INDEX M',axe)

       CALL VVSETR('YC1 - Y COORDINATE AT INDEX 1',1.)
c       CALL VVSETR('YCN - Y COORDINATE AT INDEX N',1.)
       CALL VVSETR('YCN - Y COORDINATE AT INDEX N',Float(nzset))
c        CALL VVSETI ('MAP - MAPPING FLAG',3)
c       IF ( istretch ) THEN
        CALL VVSETI ('MAP - MAPPING FLAG',3)
c       ELSE
c        CALL VVSETI ('MAP - MAPPING FLAG',0)
c       ENDIF
       
       CALL VVSETI('CPM', -1)

      ELSEIF ( ipltyp .eq. 1 ) THEN
       call set(axb,axe,ayb,aye,0.,1.,0.,1.,1)
       write(0,*) 'axb,axe,ayb,aye',axb,axe,ayb,aye
c       CALL VVSETI('CPM', -1)
c        CALL VVSETI ('MAP - MAPPING FLAG',0)
 
       CALL VVSETI('SET - DO SET CALL FLAG',0)
       CALL VVSETR('XC1 - X COORDINATE AT INDEX 1',0.)
       CALL VVSETR('XCM - X COORDINATE AT INDEX M',1.)
       CALL VVSETR('YC1 - Y COORDINATE AT INDEX 1',0.)
       CALL VVSETR('YCN - Y COORDINATE AT INDEX N',1.)
       CALL VVSETI ('MAP - MAPPING FLAG',0)
      
      ENDIF

      
c      call velvct(xvec,nmax,yvec,nmax,n1,n2,flov,hiv,
c     >            nsetv,lengthv,ispv,spv)

       write(0,*) 'lengthv,flov',lengthv,flov
       
       rdm = 0.0
       idm = 0
       
       IF ( ivint .lt. 10 ) THEN
       inx = ivint
       iny = ivint
       ELSE
       inx = ivint/10            ! use tens digit for x interval
       iny = ivint - 10*inx      ! use ones digit for y interval
       ENDIF
       CALL VVSETI('XIN',inx)
       CALL VVSETI('YIN',iny)
       
c       ix = 2
c       write(6,*) 'ix,kstrzb',ix,kstrzb
c       DO jy = 1,n2
cc       DO ix = 2,2
c        write(6,*) 'i,j,u,v',ix,jy,xvec(ix,jy),yvec(ix,jy),
c     :    0.5*(v(ib+1,jb,jy+kstrzb,1)+v(ib+2,jb,jy+kstrzb,1))
cc       ENDDO
c       ENDDO
       
       call vvinit(xvec,nmax,yvec,nmax,rdm,idm,n1,n2,rdm,idm)
       IF ( lengthv .gt. 0 ) THEN
         call vvsetr('VRM', lengthv)
         CALL VVSETR('AMN - Arrow Head Minimum Size',0.001)
       ENDIF
       IF ( cll .gt. 0.0 .and. cll .lt. 1.9999 ) THEN
          write(6,*) 'setting LWD to ',cll
          call vvsetr('LWD', cll)
       ENDIF
       call vvsetr('VLC',flov)
c       vmx = hiv
c       vmn = flov
       call vvsetr('VHC',hiv)
      CALL VVGETR('VMX - Vector Maximum Magnitude', VMX)
      CALL VVGETR('VMN - Vector Minimum Magnitude', VMN)
      
      
      CALL VVSETC('MNT - Minimum Vector Text String', ' ')  ! turns off min vector text block


      IF ( iplotrefvector >= 1 ) THEN
        CALL VVSETR('MXX', vectoffsetx)  ! - Maximum Vector Text Block X Coordinate
        CALL VVSETR('MXY', vectoffsety) ! - Maximum Vector Text Block Y Coordinate
        CALL VVSETC('MXT - Maximum Vector Text String', 'Ref. Vector')  

        CALL VVSETR('MXS', vectcharsize) ! - Maximum Vector Text Block Character Size
      ELSE
        CALL VVSETR('MXY', -0.2) ! - Maximum Vector Text Block Y Coordinate
        CALL VVSETC('MXT - Maximum Vector Text String', '')  

        CALL VVSETR('MXS', 0.015) ! - Maximum Vector Text Block Character Size
      ENDIF
      
c      CALL VVSETR('VLC - Vector Low Cutoff Value', VMN+0.1*(VMX-VMN))
       CALL VVGETR('DMX - NDC Maximum Vector Size', DMX)

      CALL GETSET(VL,VR,VB1,VT,UL,UR,UB,UT,LL)
      IF ( lengthv .gt. 0.0 ) THEN
       VRL =  DMX / (VR - VL) *vmx/lengthv
       CALL VVSETR('VRL - Vector Reference Length', vrl)

c       CALL VVSETR('VRL - Vector Reference Length', dmx*vmx/lengthv)
      ENDIF 
       write(0,*) 'vmx,vmn,dmx,vrl,vrl0 = ',
     :    vmx,vmn,dmx,vrl,DMX / (VR - VL)
       
         ishdclr = Abs( ishade )
        IF ( ishdclr .gt. 0 ) THEN
         call GSPLCI(ishdclr)
        ELSEIF ( iotyp .ge. 20 ) THEN
         call GSPLCI(99)
        ENDIF
       call vvectr(xvec,yvec,rdm,idm,idm,rdm)

         call GSPLCI(1)
      
      end if !  
      end if ! ( ipltyp .ne. 4 )
c
      write(*,*) 'idebg=25'
c
      if ( ipltyp .eq. 4 ) then
      mxyz = max(nxp,nyp)
      mxyz2 = max(mxyz,nzp) +2
      isign = 1
      call dfclr(idx,01,irecmj,irecmn,irectx,ndot)
      call gsplci(irecmj)
      print*,nx,nxp,ny,nyp,nz,nzp
      call isosrf(plt3dp,nx,nxp,ny,nyp,nzp,eye,mxyz2,
     >            slab3,tiso,isoflg)
      isign = -1
      call dfclr(idx,-1,irecmj,irecmn,irectx,ndot)
      call gsplci(irecmj)
      call isosrf(plt3dn,nx,nxp,ny,nyp,nzp,eye,mxyz2,
     >            slab3,tiso,isoflg)
      do  k = 1,nz
      do  j = 1,ny
      do  i = 1,nx
       plt3dn(i,j,k) = 0.0
      enddo
      enddo
      enddo
      call gsplci(1)
      call isosrf(plt3dn,nx,nxp,ny,nyp,nzp,eye,mxyz2,
     >            slab3,1000.,1)
      end if
c
      write(*,*) 'idebg=26'
c
      IF ( ioverlay .gt. 0 .or. idone .eq. 1 ) THEN
        IF ( iotyp .eq. 21 ) then
         write(0,*) 'done eps frame ',ieps
         call clsgks

        ieps = ieps + 1
        epsflag = 1
        
        ELSE
        call frame
        ENDIF
      ENDIF
      write(*,*) 'idebg=27'
c
      if ( idone .eq. 1 ) ioverlay = Max(1,ioverlay)
      write(*,*) 'idebg=28'
c
c
      go to 2210
c      write(*,*) 'idebg=29'
c
 9010 continue
      write(*,*) 'idebg=30'
c
c
c      write(lmin,911) max( rclon, awrk )
c      label1 = '                    X-Z SLAB'
c      label6 = 'y(km) ='//lloc
c      write(lloc,911) dx*xfac*float(ix-1)
c      label1 = '                    Y-Z SLAB'
c      label6 = 'x(km) ='//lloc
c
c  
c
c    END OF PROGRAM
c
c
c
c      IF ( time .eq. lasttime )
      
c      call clsgks
c      close (unit=28)
c
c
c
      return
      end
c
c
c  end of program
c
c
       subroutine setplottitles(ntit,title,ibase,ibase_micth,ibase_mic_add)
       
       USE MICRO_MODULE
       USE INDEX_MODULE, only: lis,lnis
       
       implicit none
       integer ntit,ibase,ibase_micth,ibase_mic_add
      
       character*40 title(ntit)
       
       integer k
c
c  titles of plots
c
       
c
c noninductive charging rates: snow-ice
      title(800) = 'SCSACI (pC/(m**3 s))'
      title(801) = 'SCSACIP (pC/(m**3 s))'
      title(802) = 'SCSACIR (pC/(m**3 s))'
      
c noninductive charging rates: low den graup-snow,ice
      title(803) = 'SCGLACS (pC/(m**3 s))'
      title(804) = 'SCGLACI (pC/(m**3 s))'
      title(805) = 'SCGLACIR (pC/(m**3 s))'
      title(806) = 'SCGLACIP (pC/(m**3 s))'
      
c noninductive charging rates: med den graup-snow,ice
      title(807) = 'SCGMACS (pC/(m**3 s))'
      title(808) = 'SCGMACI (pC/(m**3 s))'
      title(809) = 'SCGMACIR (pC/(m**3 s))'
      title(810) = 'SCGMACIP (pC/(m**3 s))'
      
c noninductive charging rates: high den graup-snow,ice
      title(811) = 'SCGHACS (pC/(m**3 s))'
      title(812) = 'SCGHACI (pC/(m**3 s))'
      title(813) = 'SCGHACIR (pC/(m**3 s))'
      title(814) = 'SCGHACIP (pC/(m**3 s))'
      
c noninductive charging rates: hail-snow,ice
      title(815) = 'SCHACS (pC/(m**3 s))'
      title(816) = 'SCHACI (pC/(m**3 s))'
      title(817) = 'SCHACIR (pC/(m**3 s))'
      title(818) = 'SCHACIP (pC/(m**3 s))'
      
c noninductive charging rates: large hail-snow,ice
      title(819) = 'SCHLACS (pC/(m**3 s))'
      title(820) = 'SCHLACI (pC/(m**3 s))'
      title(821) = 'SCHLACIR (pC/(m**3 s))'
      title(822) = 'SCHLACIP (pC/(m**3 s))'

c noninductive charging rates: frozen drops-snow,ice
      title(823) = 'SCFACS (pC/(m**3 s))'
      title(824) = 'SCFACI (pC/(m**3 s))'
      title(825) = 'SCFACIR (pC/(m**3 s))'
      title(826) = 'SCFACIP (pC/(m**3 s))'

c grouped noninductive charging rates: 
      title(827) = 'Graupel-ice net (pC/(m**3 s))'
      title(828) = 'Graupel-snow net (pC/(m**3 s))'
      title(829) = 'Hail/FD-ice net (pC/(m**3 s))'
      title(830) = 'Hail/FD-snow net (pC/(m**3 s))'

c inductive charging rates:
      title(831) = 'SCSACW (pC/(m**3 s))'
      title(832) = 'SCGLACW (pC/(m**3 s))'
      title(833) = 'SCGMACW (pC/(m**3 s))'
      title(834) = 'SCGHACW (pC/(m**3 s))'
      title(835) = 'SCHACW (pC/(m**3 s))'
      title(836) = 'SCHLACW (pC/(m**3 s))'
      title(837) = 'SCFACW (pC/(m**3 s))'

c grouped inductive charging rates:
      title(838) = 'Graupel tot ind. chg rate (pC/(m**3 s))'
      title(839) = 'Hail/FD tot ind. chg rate (pC/(m**3 s))'

c rime accretion rates:
      title(840) = 'RAR Snow (g/(m**2 s))'
      title(841) = 'RAR GL (g/(m**2 s))'
      title(842) = 'RAR GM (g/(m**2 s))'
      title(843) = 'RAR GH (g/(m**2 s))'
      title(844) = 'RAR H (g/(m**2 s))'
      title(845) = 'RAR HL (g/(m**2 s))'
      title(846) = 'RAR FD (g/(m**2 s))'

c more group noninductive rates 
      title(847) = 'SCGL-ICE NET (pC/(m**3 s))'
      title(848) = 'SCGM-ICE NET (pC/(m**3 s))'
      title(849) = 'SCGH-ICE NET (pC/(m**3 s))'
      title(850) = 'SCH-ICE NET (pC/(m**3 s))'
      title(851) = 'SCHL-ICE NET (pC/(m**3 s))'
      title(852) = 'SCFD-ICE NET (pC/(m**3 s))'

      title(853) = 'SCS_I'
      
      title(854) = 'CGHI_NET'
      title(855) = 'CGHS_NET'
      title(856) = 'CGHIS_NET'
      
      title(857) = 'CGHW_NET'
      title(858) = 'CGH_NET'

      title(1) = 'POT. TEMP. PERTURBATION (K) '
      title(2) = 'TEMPERATURE (C)'
      title(3) = 'PERTURBATION U -WIND (M/S)'
      title(4) = 'PERTURBATION V -WIND (M/S)'
      title(5) = 'W (M/S)'
      title(6) = 'TOTAL U -WIND (M/S)  '
      title(7) = 'TOTAL V -WIND (M/S)  '
      title(8) = 'PERTURBATION PRESSURE (MB) '
      title(9) = 'MIXING COEF.(M2/S2) FROM TKE EQ.'
      title(10) = 'WATER VAPOR MIXING RATIO (G/KG) '
      title(11) = 'CLOUD WATER MIXING RATIO (G/KG)'
      title(12) = 'RAIN WATER MIXING RATIO (G/KG) '
      title(13) = 'CLOUD ICE MIXING RATIO (G/KG)  '
      title(14) = 'SNOW AGGREGATE MIXING RATIO (G/KG)'
      title(15) = 'TOTAL HAIL MIXING RATIO (G/KG)'
      title(16) = 'CLOUD WATER-ICE (G/KG)'
      title(17) = 'VISIBLE CLOUD (G/KG)'
      title(18) = 'TOTAL ICE (G/KG)'
      title(19) = 'Total Liquid Water Content (g/m^3)'
      title(20) = 'TOTAL PRECIPITATION WATER (G/KG)'
      title(21) = 'TOTAL CONDENSATE (G/KG)'
      title(22) = 'RADAR (DBZ)'
      title(23) = 'WIND VECTOR'  ! wind vectors?
      title(24) = 'RELATIVE HUMIDITY (*100)'
      title(25) = 'HORIZONTAL DIV/CONVERGENCE ( *1000/S)'
      title(26) = 'VERTICAL VORTICITY ( *1000/S)'
      title(27) = 'WINDSPEED (M/S)'
      title(28) = 'hack point'
      title(29) = 'THETA-E (K)'
      title(30) = 'DELh[W] (*1000/s)'
      title(32) = 'ELECTRIC FIELD VECTOR'
      title(33) = 'Potential Temperature (K) '
      title(34) = 'TEMPERATURE (K)'
      title(35) = 'COMPOSITE RADAR (DBZ)'
      title(36) = 'THETA-V (K)'
      title(37) = 'THETA-V Perturbation (K)'
      title(38) = 'Ice Water Path kg/m^2'
      title(39) = 'Net Buoyancy 10e-2 m/s^2'
      title(40) = 'Therm. Buoyancy 10e-2 m/s^2'
      title(41) = 'Load Buoyancy 10e-2 m/s^2'
c 
      title(42) = 'QSW (G/KG)'
      title(43) = 'QSI (G/KG)'
      title(44) = 'QSS (G/KG)'
      title(45) = 'SSI (%)'
      title(46) = 'SSW (%)'
      title(47) = 'CWDIA (micron)'
      title(48) = 'CIDIA (micron)'
      title(49) = 'RAIN Mean Vol. Diam (micron)'
      title(50) = 'SNOW Mean Vol. Diam (micron)'
      title(51) = 'Graupel Mean Vol. Diam (micron)'
      title(52) = 'Hail Mean Vol. Diam (micron)'
!      title(52) = 'CWDN (KG/M3)'
      title(53) = 'CIDN (KG/M3)'
      title(54) = 'RWDN (KG/M3)'
      title(55) = 'SWDN (KG/M3)'
      title(56) = 'Pressure (mb)'
      title(57) = 'Pressure (mb)' ! with mark at approx. lowest pressure
      title(28) = 'HWDN (KG/M3)'
      title(57) = 'VTWBAR (M/S)'
      title(58) = 'VTIBAR (M/S)'
      title(59) = 'VTHBAR (M/S)'
      title(60) = 'VTRBAR (M/S)'
      title(61) = 'VTSBAR (M/S)'
      title(62) = 'QRACW (*1E6 1/S)'
      title(63) = 'QSACW-T (*1E6 1/S)'
      title(64) = 'QSACW-S (*1E6 1/S)'
      title(65) = 'QHACW (*1E6 1/S)'
      title(66) = 'QRACI-T (*1E6 1/S)'
      title(67) = 'QRACI-S (*1E6 1/S)'
      title(68) = 'QRACI-H (*1E6 1/S)'
      title(69) = 'QIACR-T (*1E6 1/S)'
      title(70) = 'QIACR-S (*1E6 1/S)'
      title(71) = 'QIACR-H (*1E6 1/S)'
      title(72) = 'QSACI (*1E6 1/S)'
      title(73) = 'QHACI (*1E6 1/S)'
      title(74) = 'QHACR (*1E6 1/S)'
      title(75) = 'QHACS (*1E6 1/S)'
      title(76) = 'QRACS-T (*1E6 1/S)'
      title(77) = 'QRACS-1 (*1E6 1/S)'
      title(78) = 'QRACS-2 (*1E6 1/S)'
      title(79) = 'QSACR-T (*1E6 1/S)'
      title(80) = 'QSACR-2 (*1E6 1/S)'
      title(81) = 'QSFW (*1E6 1/S)'
      title(82) = 'QSFI (*1E6 1/S)'
      title(83) = 'QRCNW (*1E6 1/S)'
      title(84) = 'QSCNI (*1E6 1/S)'
      title(85) = 'QHCNS (*1E6 1/S)'
      title(86) = 'QRFRZ-T (*1E6 1/S)'
      title(87) = 'QRFRZ-1 (*1E6 1/S)'
      title(88) = 'QRFRZ-2 (*1E6 1/S)'
      title(89) = 'QSMLR (*1E6 1/S)'
      title(90) = 'QHMLR (*1E6 1/S)'
      title(91) = 'QSDSV (*1E6 1/S)'
      title(92) = 'QHDSV (*1E6 1/S)'
      title(93) = 'QHSHR (*1E6 1/S)'
      title(94) = 'QHWET (*1E6 1/S)'
      title(95) = 'QHDRY (*1E6 1/S)'
      title(96) = 'QRCEV (*1E6 1/S)'
      title(97) = 'QSCEV (*1E6 1/S)'
      title(98) = 'QHCEV (*1E6 1/S)'
      title(99) = 'QIINT (*1E6 1/S)'
      title(100) = 'CWFLX (*1E6 G/GS)'
      title(101) = 'CIFLX (*1E6 G/GS)'
      title(102) = 'RWFLX (*1E6 G/GS)'
      title(103) = 'SWFLX (*1E6 G/GS)'
      title(104) = 'HWFLX (*1E6 G/GS)'
      title(105) = 'PQWVI (*1E6 G/GS)'
      title(106) = 'PQWVD (*1E6 G/GS)'
      title(107) = 'PQCWI (*1E6 G/GS)'
      title(108) = 'PQCWD (*1E6 G/GS)'
      title(109) = 'PQIWI (*1E6 G/GS)'
      title(110) = 'PQIWD (*1E6 G/GS)'
      title(111) = 'PQRWI (*1E6 G/GS)'
      title(112) = 'PQRWD (*1E6 G/GS)'
      title(113) = 'PQSWI (*1E6 G/GS)'
      title(114) = 'PQSWD (*1E6 G/GS)'
      title(115) = 'PQHWI (*1E6 G/GS)'
      title(116) = 'PQHWD (*1E6 G/GS)'
      title(117) = 'QCWCE (*1E6 G/GS)'
      title(118) = 'QCIDS (*1E6 G/GS)'
      title(119) = 'QWFZI (*1E6 G/GS)'
      title(120) = 'QIMLW (*1E6 G/GS)'
      title(121) = 'LOADING: G/G'
      title(122) = 'TOTAL HEATING: C/S'
      title(123) = 'DEP/SUB HEATING: C/S'
      title(124) = 'FRZ/MLT HEATING: C/S'
      title(125) = 'COND/EVAP HEATING: C/S'
      title(126) = 'ERI' 
      title(127) = 'ESI'
      title(128) = 'EHI'
      title(129) = 'ERW'
      title(130) = 'ESW'
      title(131) = 'EHW'
      title(132) = 'ERS'
      title(133) = 'EHR'
      title(134) = 'EHS'
      title(135) = 'EIC'
      title(136) = 'ESCI'
      title(137) = 'EHSCNV'
c
c     title(140) = 'SMALL HAIL EMBRYOS: 0.1<D (CM) (G/G)'
c     title(141) = 'HAIL EMBRYOS AND GRAUPEL: 0.1<D<0.5(CM) (G/G)'
c     title(142) = 'SMALL HAIL: 0.5<D<1.0(CM) (G/G)'
c     title(143) = 'MEDIUM HAIL: 1.0<D<2.5(CM) (G/G)'
c     title(144) = 'LARGE HAIL: 2.5<D<4.0(CM) (G/G)'
c     title(145) = 'GIANT HAIL: 4.0<D(CM) (G/G)'
c
      title(140) = 'SMALL GRAUPEL: 0.05<D (CM) (KG/KG)'
      title(141) = 'GRAUPEL: 0.05<D<0.5 (CM) (KG/KG)'
      title(142) = 'SMALL HAIL: 0.5<D<2.0 (CM) (KG/KG)'
      title(143) = 'LARGE HAIL: 2.0>D (CM) (KG/KG)'
c
      title(146) = 'HDN: TOTAL HAIL: 4.0<D(CM) (G/G)'
      title(147) = 'HDN: SMALL HAIL EMBRYOS: 0.1<D (CM) (G/G)'
      title(148) = 'HDN: HAIL EMBRYOS AND GRAUPEL: 0.1<D<0.5(CM) (G/G)'
      title(149) = 'HDN: SMALL HAIL: 0.5<D<1.0(CM) (G/G)'
      title(150) = 'HDN: MEDIUM HAIL: 1.0<D<2.5(CM) (G/G)'
      title(151) = 'HDN: LARGE HAIL: 2.5<D<4.0(CM) (G/G)'
      title(152) = 'HDN: GIANT HAIL: 4.0<D(CM) (G/G)'
c
      title(180) = 'SHEAR HAZARD REGION'
      title(181) = 'HORIZ. SHEAR X (1.0E+03) >2.5/S'
      title(182) = 'VERT. WIND SHEAR X (1.0E+02) >8.3/S'
      title(183) = 'VERT. WIND HAZARD REGION >2.0 M/S'

      title(221) = 'dBZ Net Oversorting'
      title(222) = 'dBZ Undersorting'
      title(223) = 'dBZ Over/undersorting'
      title(224) = 'dBZ Gross Oversorting'
      title(225) = 'dBZ Size Limiter Change'
      title(226) = 'Hail number rate from graupel' ! # m-3 s-1 averaged over thistory
      title(227) = 'Hail number rate from frozen drops' ! # m-3 s-1 averaged over thistory
      title(228) = 'Total Hail number rate' ! # m-3 s-1 averaged over thistory

      title(239) = 'TOTAL HAIL (N/M**3)'
      title(240) = 'SMALL GRAUPEL: 0.05<D (CM) (N/M**3)'
      title(241) = 'GRAUPEL: 0.05<D<0.5 (CM) (N/M**3)'
      title(242) = 'SMALL HAIL: 0.5<D<2.0 (CM) (N/M**3)'
      title(243) = 'LARGE HAIL: 2.0>D (CM) (N/M**3)'
c
c
      title(251) = 'HORIZ. W GRAD. (*1.0E+3 M/S**2)'
      title(252) = 'VERT. PRESSURE GRAD. (*1.0E+2 N/KG)'
c
c These are for SAM binary output:
      title(261)  = 'PRECIP PRECIP RATE (MM/HR)'
      title(262)  = 'PRECIP-ICE PRECIP RATE (MM/HR)'
      title(263)  = 'PRECIP-LIQ PRECIP RATE (MM/HR)'
      title(264)  = 'CLOUD WATER PRECIP RATE (MM/HR)'
      title(265)  = 'RAIN WATER PRECIP RATE (MM/HR)'
      title(266)  = 'CLOUD ICE-COLUMN PRECIP RATE (MM/HR)'
      title(267)  = 'SNOW AGGREGATE PRECIP RATE (MM/HR)'
      title(268)  = 'LOW DEN GRAUPEL PRECIP RATE (MM/HR)'
      title(269)  = 'MED DEN GRAUPEL PRECIP RATE (MM/HR)'
      title(270)  = 'HI DEN GRAUPEL PRECIP RATE (MM/HR)'
      title(271)  = 'FROZEN DROP PRECIP RATE (MM/HR)'
      title(272)  = 'HAIL PRECIP RATE (MM/HR)'
      title(273)  = 'LARGE HAIL PRECIP RATE (MM/HR)'
      title(274)  = 'CLOUD ICE-RIMED PRECIP RATE (MM/HR)'
      title(275)  = 'CLOUD ICE-PLATE PRECIP RATE (MM/HR)'

! for COMMAS 2D precip fields
      title(281)  = 'PRECIP PRECIP RATE (MM/HR)'
      title(282)  = 'PRECIP-ICE PRECIP RATE (MM/HR)'
      title(283)  = 'PRECIP-LIQ PRECIP RATE (MM/HR)'
      title(284)  = 'Graupel Precip Rate (MM/HR)'
      title(285)  = 'Hail Precip Rate (MM/HR)'
      title(286)  = 'Frozen drop Precip Rate (MM/HR)'
      title(287)  = 'Hail % from Frozen drops mass fraction'

      title(288)  = 'Rain-auto Accum. (mm)'
      title(289)  = 'Rain-melt Accum. (mm)'
      title(290)  = 'Rain-shed Accum. (mm)'

      title(291)  = 'Rain Accum. (mm)'
      title(292)  = 'PRECIP-ICE PRECIP RATE (mm)'
      title(293)  = 'Rain2 Accum. (mm) (grid mot.)'
!      title(293)  = 'PRECIP-LIQ PRECIP RATE (mm)'
      title(294)  = 'Graupel Precip Accum. (mm)'
      title(295)  = 'Hail Precip Accum. (mm)'
      title(296)  = 'Hail from Frozen drop Precip Accum. (mm)'
      title(297)  = 'Hail % from Frozen drops num. fraction'
      title(298)  = 'Hail mean diameter at sfc. (mm) (thistory)'
      title(299)  = 'Hail mean diameter at sfc. (mm) (total)'
      title(300)  = 'Hail % from Frozen drops num. fraction (dthis)'
      title(301)  = 'Hail Precip Accum. (mm) (dthis)'

      title(302)  = 'FD from Rain-Auto (number)'
      title(303)  = 'FD from Rain-Shed (number)'
      title(304)  = 'FD from Rain-Melt (number)'
      title(305)  = 'FD from Rain-Auto2 (number)' ! grid motion
      title(306)  = 'FD from Rain-Shed2 (number)' ! grid motion
      title(307)  = 'FD from Rain-Melt2 (number)' ! grid motion

      title(308)  = 'Graupel Precip Accum.2 (mm)' ! grid motion
      title(309)  = 'Hail Precip Accum.2 (mm)' ! grid motion

      title(310)  = 'Graupel mean diameter at sfc. (mm) (thistory)'
      title(311)  = 'Graupel mean diameter at sfc. (mm) (total)'
      title(312)  = 'FD Precip Accum. (mm)'
      title(313)  = 'FD Precip Accum.2 (mm)' ! grid motion
      title(314)  = 'Hail mean diameter2 at sfc. (mm) (total)'
      title(315)  = 'Graupel mean diameter2 at sfc. (mm) (total)'
     
c
!      title(332) = 'DELTA POT. TEMP. PERURBATION (K)'
!      title(333) = 'DELTA U -WIND (M/S)'
!      title(334) = 'DELTA V -WIND (M/S)'
!      title(335) = 'DELTA W -WIND (M/S)'
!      title(336) = 'DELTA PRESSURE (MB)'
!      title(337) = 'DELTA MIXING COEF.(M2/S2) FROM TKE EQ.'
!      title(338) = 'DELTA WATER VAPOR MIXING RATIO (G/KG)'
!      title(339) = 'DELTA CLOUD WATER MIXING RATIO (G/KG)'
!      title(340) = 'DELTA RAIN WATER MIXING RATIO (G/KG)'
!      title(341) = 'DELTA CLOUD ICE MIXING RATIO (G/KG)'
!      title(342) = 'DELTA SNOW AGGREGATE MIXING RATIO (G/KG)'
!      title(343) = 'DELTA TOTAL HAIL MIXING RATIO (G/KG)'
c
      title(320) = 'LIGHTNING MAP'
      title(321) = 'XZ LIGHTNING FOOTPRINT MAP'
      title(322) = 'YZ LIGHTNING FOOTPRINT MAP'
      title(323) = 'XZ LIGHTNING FOOTPRINT MAP'

      title(324) = 'Positive Leader Composite'  ! no isum integration
      title(325) = 'Negative Leader Composite'  ! no isum integration
      title(326) = 'Flash Init Points'          ! no isum integration

      title(327) = 'Positive Leader Composite'
      title(328) = 'Negative Leader Composite'
      title(329) = 'Flash Init Points'
      title(330) = 'Total Leader Composite'
      title(331) = 'Flash Extent Density'
      title(332) = 'Flash Extent Density (GLM 8km)'
      title(333) = 'Flash Extent Density (GLM 8km)' ! hack version using iflshr max within box
c
      title(393) = 'Droplet + Rain charge (nC m-3)'
      title(394) = 'TOTAL CLOUD ICE CHARGE (nC/m**3)'
      title(395) = 'TOTAL GRAUPEL + HAIL CHARGE (nC/m**3)'
      IF ( microp(1:5) == 'ICE10') THEN
      title(396) = 'LARGE HAIL CHARGE (nC/m**3)'
      ELSE
      title(396) = 'HAIL CHARGE (nC/m**3)'
      ENDIF
      title(399) = 'CLOUD WATER CHARGE (nC/m**3)'
      title(400) = 'CLOUD ICE-COLUMN CHARGE (nC/m**3)'
      title(401) = 'CLOUD ICE-PLATE CHARGE (nC/m**3)'
      title(402) = 'CLOUD ICE-RIMED CHARGE (nC/m**3)'
      title(403) = 'RAIN WATER CHARGE (nC/m**3)'
      title(404) = 'SNOW AGGREGATE CHARGE (nC/m**3)'
      title(405)  = 'LOW DEN GRAUPEL CHARGE (nC/m**3)' 
      title(406)  = 'MED DEN GRAUPEL CHARGE (nC/m**3)'
      title(407)  = 'HI DEN GRAUPEL CHARGE (nC/m**3)'
      title(408)  = 'FROZEN DROP CHARGE (nC/m**3)'
      IF ( microp(1:5) == 'ICE10') THEN
      title(409) = 'HAIL CHARGE (nC/m**3)'
      ELSE
      title(409) = 'GRAUPEL CHARGE (nC/m**3)'
      ENDIF
      title(410) = 'Net Charge (nC/m**3)'
      title(411) = 'ELECTRIC FIELD: X-DIR (kV/m)'
      title(412) = 'ELECTRIC FIELD: Y-DIR (kV/m)'
      title(413) = 'ELECTRIC FIELD: Z-DIR (kV/m)'
      title(414) = 'ELECTRIC FIELD: |E| (kV/m)'
      title(415) = 'ELECTRIC POTENTIAL'
      title(416) = 'TOTAL GRAUPEL CHARGE (nC/m**3)'
      title(417) = 'TOTAL ICE+Snow CHARGE (nC/m**3)'
      title(418) = 'TOTAL CHARGE X (nC/m**3)'
      title(419) = 'E-Field: Ez at w-pnts (kV/m)'
c
      title(420) = 'CLOUD WATER CONCENTRATION (#/cm**3)'
      title(421) = 'CLOUD ICE-COLUMN CONCENTRATION (#/L)'
      title(422) = 'CLOUD ICE-PLATE CONCENTRATION (#/L)'
      IF ( lnis >= 1 ) THEN
      title(423) = 'CLOUD ICE-Nonint CONCENTRATION (#/L)'
      ELSE
      title(423) = 'CLOUD ICE-RIMED CONCENTRATION (#/L)'
      ENDIF
      title(424) = 'RAIN WATER CONCENTRATION (#/L)'
      title(425) = 'SNOW AGGREGATE CONCENTRATION (#/L)'
      title(426)  = 'LOW DEN GRAUPEL CONCENTRATION (#/m**3)' 
      title(427)  = 'MED DEN GRAUPEL CONCENTRATION (#/m**3)'
      title(428)  = 'HI DEN GRAUPEL CONCENTRATION (#/m**3)'
      title(429)  = 'FROZEN DROP CONCENTRATION (#/m**3)'
      title(430) = 'TOTAL HAIL CONCENTRATION (#/m**3)'
      title(431) = 'TOTAL GRAUPEL CONCENTRATION (#/m**3)'
      title(432) = 'TOTAL ICE CONCENTRATION (#/m**3)'
c
      title(434) = 'SSW - Supersat. wrt liquid' ! supersaturation w.r.t. liquid water before microphysics
      title(435) = 'SSI - Supersat. wrt ice' ! supersaturation w.r.t. ice
c
      title(440) = 'Chaff Fiber Concentration (#/1000 m**3)' ! number of chaff fibers
      title(441) = 'LNO Concentration (ppbv)' ! lightning-produced NO/NOx
      title(442) = 'Ambient NO Concentration (ppbv)' ! Ambient (background) NO/NOx
      title(443) = 'Ambient CO Concentration (ppbv)' ! Ambient (background) CO

c      title(440) = 'CLOUD WATER TEMPERATURE (K)'
c      title(441) = 'CLOUD ICE-COLUMN TEMPERATURE (K)'
c      title(442) = 'CLOUD ICE-PLATE TEMPERATURE (K)'
c      title(443) = 'CLOUD ICE-RIMED TEMPERATURE (K)'
c      title(444) = 'RAIN WATER TEMPERATURE (K)'
c      title(445) = 'SNOW AGGREGATE TEMPERATURE (K)'
      title(446)  = 'LOW DEN GRAUPEL TEMPERATURE (K)' 
      title(447)  = 'MED DEN GRAUPEL TEMPERATURE (K)'
      title(448)  = 'HI DEN GRAUPEL TEMPERATURE (K)'
      title(449)  = 'FROZEN DROP TEMPERATURE (K)'
!      title(450) = 'TOTAL HAIL TEMPERATURE (K)'
      title(450) = 'Cloud Top Temperature (C)'
      title(451) = 'CCN Concentration (#/cm**3)'
      title(452) = 'Activated CCN (#/cm**3)'
      title(453) = 'Activated CIN  (#/cm**3)'
c
      title(460) = 'CLOUD WATER MEAN CLOUD AMOUNT'
      title(461) = 'CLOUD ICE-COLUMN MEAN CLOUD AMOUNT'
      title(462) = 'CLOUD ICE-PLATE MEAN CLOUD AMOUNT'
      title(463) = 'CLOUD ICE-RIMED MEAN CLOUD AMOUNT'
      title(464) = 'RAIN WATER MEAN CLOUD AMOUNT'
      title(465) = 'SNOW AGGREGATE MEAN CLOUD AMOUNT'
      title(466)  = 'LOW DEN GRAUPEL MEAN CLOUD AMOUNT' 
      title(467)  = 'MED DEN GRAUPEL MEAN CLOUD AMOUNT'
      title(468)  = 'HI DEN GRAUPEL MEAN CLOUD AMOUNT'
      title(469)  = 'FROZEN DROP MEAN CLOUD AMOUNT'
      title(470) = 'TOTAL HAIL MEAN CLOUD AMOUNT'
c
      title(480) = 'CLOUD WATER TAU (SEC)'
      title(481) = 'CLOUD ICE-COLUMN TAU (SEC)'
      title(482) = 'CLOUD ICE-PLATE TAU (SEC)'
      title(483) = 'CLOUD ICE-RIMED TAU (SEC)'
      title(484) = 'RAIN WATER TAU (SEC)'
      title(485) = 'SNOW AGGREGATE TAU (SEC)'
      title(486)  = 'LOW DEN GRAUPEL TAU (SEC)' 
      title(487)  = 'MED DEN GRAUPEL TAU (SEC)'
      title(488)  = 'HI DEN GRAUPEL TAU (SEC)'
      title(489)  = 'FROZEN DROP TAU (SEC)'
      title(490) = 'TOTAL HAIL TAU (SEC)'
c
c
      title(503)  = 'Cloud Water Content (g/m**3)'
      title(504)  = 'CLOUD WATER Mixing Ratio (g/kg)'
      title(505)  = 'RAIN WATER Mixing Ratio (g/kg)'
      title(506)  = 'CLOUD ICE-COLUMN Mixing Ratio (g/kg)'
      title(507)  = 'SNOW AGGREGATE Mixing Ratio (g/kg)'
      title(508)  = 'LOW DEN GRAUPEL Mixing Ratio (g/kg)'
      title(509)  = 'MED DEN GRAUPEL Mixing Ratio (g/kg)'
      title(510)  = 'HI DEN GRAUPEL Mixing Ratio (g/kg)'
      title(511)  = 'Frozen Drop Mixing Ratio (g/kg)'
      title(512)  = 'HAIL Mixing Ratio (g/kg)'
      title(513)  = 'LARGE HAIL Mixing Ratio (g/kg)'
      
      IF ( microp(1:5) .eq. 'ICE10' ) THEN
      title(512)  = 'HAIL Mixing Ratio (g/kg)'
      title(513)  = 'LARGE HAIL Mixing Ratio (g/kg)'
      ELSEIF ( microp(1:3) .eq. 'LFO' ) THEN
      title(512)  = 'Graupel/Hail Mixing Ratio (g/kg)'
      title(513)  = 'LARGE HAIL Mixing Ratio (g/kg)'
      ELSEIF ( microp(1:1) .eq. 'Z' ) THEN
      title(512)  = 'Graupel Mixing Ratio (G/KG)'
      title(513)  = 'Hail Mixing Ratio (G/KG)'
      ELSEIF ( microp(1:3) .eq. 'TAK' ) THEN
      title(512)  = 'Graupel Mixing Ratio (G/KG)'
      title(513)  = 'Hail Mixing Ratio (G/KG)'
      ENDIF

      IF ( lis >= 1 ) THEN
      title(514)  = 'CLOUD ICE-Nonint. Mixing Ratio (G/KG)'
      ELSE
      title(514)  = 'CLOUD ICE-RIMED Mixing Ratio (G/KG)'
      ENDIF
      title(515)  = 'CLOUD ICE-PLATE MIXING RATIO (G/KG)'
      
c
      title(516)  = 'PLATE-WARM MIXING RATIO (G/KG)'
      title(517)  = 'NEEDLE-WARM MIXING RATIO (G/KG)'
      title(518)  = 'COLUMN-WARM MIXING RATIO (G/KG)'
      title(519)  = 'SECTOR-WARM MIXING RATIO (G/KG)'
      title(520)  = 'DENDRITE MIXING RATIO (G/KG)'
      title(521)  = 'SECTOR-COLD MIXING RATIO (G/KG)'
      title(522)  = 'COLUMN-COLD MIXING RATIO (G/KG)'
      title(523)  = 'BULLET-COLD MIXING RATIO (G/KG)'
      title(524)  = 'EQUILIBRIUM MIXING RATIO (G/KG)'
      title(525)  = 'FROZEN DROP MIXING RATIO (G/KG)'

c
      title(528) = 'CLOUD ICE Columns + Plates (G/KG)'
      title(529) = 'TOTAL CLOUD ICE MIXING RATIO (G/KG)'
      title(530) = 'TOTAL GRAUPEL MIXING RATIO (G/KG)'
      title(531) = 'TOTAL CLOUD MIXING RATIO (G/KG)'
      title(532) = 'TOTAL CLOUD + SNOW MIXING RATIO (G/KG)'
      title(533) = 'TOTAL LIQUID WATER MIXING RATIO (G/KG)'
      title(534) = 'TOTAL CLOUD ICE + SNOW MIXING RATIO (G/KG)'
      title(535) = 'TOTAL PRECIP ICE MIXING RATIO (G/KG)'
      title(536) = 'TOTAL ICE WATER MIXING RATIO (G/KG)'
      title(537) = 'TOTAL CONDENSATE MIXING RATIO (G/KG)'
      title(538) = 'TOTAL PRECIPITATION MIXING RATIO (G/KG)'
      title(539) = 'TOTAL GRAUPEL + HAIL MIXING RATIO (G/KG)'
      title(540) = 'TOTAL HAIL MIXING RATIO (G/KG)'
c
      title(541)  = 'Graupel Particle Density (kg/m**3)'
      title(542)  = 'Liquid on Snow Mixing Ratio (G/KG)'
      title(543)  = 'Liquid on Graupel Mixing Ratio (G/KG)'
      title(544)  = 'Liquid on Hail Mixing Ratio (G/KG)'
      title(545)  = 'Frozen Drop Particle Density (kg/m**3)'
      title(546)  = 'Hail Particle Density (kg/m**3)'
      title(547)  = 'Graupel Particle Density (dualpol) (kg/m**3)'
      title(548)  = 'Hail Particle Density (dualpol) (kg/m**3)'
      
c
c reflectivities (6th moments):
c
      title(570)  = 'ZRW (dBZ)' 
      title(571)  = 'ZSW (dBZ)' 
      title(572)  = 'ZGL (dBZ)' 
      title(573)  = 'ZGM (dBZ)' 
      title(574)  = 'ZGH (dBZ)' 
      title(575)  = 'ZF (dBZ)' 
      title(576)  = 'ZHW (dBZ)' 
      title(577)  = 'ZHL (dBZ)' 
      title(578)  = 'ZRW2 (dBZ)' 
      title(579)  = 'ZCI (dBZ)' 
c
      title(581)  = 'PRECIP TOTAL (MM)'
      title(582)  = 'PRECIP-ICE TOTAL (MM)'
      title(583)  = 'PRECIP-LIQ TOTAL (MM)'
      title(584)  = 'CLOUD WATER TOTAL (MM)'
      title(585)  = 'RAIN WATER TOTAL (MM)'
      title(586)  = 'CLOUD ICE-COLUMN TOTAL (MM)'
      title(587)  = 'SNOW AGGREGATE TOTAL (MM)'
      title(588)  = 'LOW DEN GRAUPEL TOTAL (MM)'
      title(589)  = 'MED DEN GRAUPEL TOTAL (MM)'
      title(590)  = 'HI DEN GRAUPEL TOTAL (MM)'
      title(591)  = 'FROZEN DROP TOTAL (MM)'
      title(592)  = 'HAIL TOTAL (MM)'
      title(593)  = 'LARGE HAIL TOTAL (MM)'
      title(594)  = 'CLOUD ICE-RIMED TOTAL (MM)'
      title(595)  = 'CLOUD ICE-PLATE TOTAL (MM)'
c
      title(601) = 'Zh (dBZ) (dualpara)'
      title(602) = 'Zv (dB)'
      title(603) = 'Zdr (dB)'
      title(604) = 'Zdp (dB)'
      title(605) = 'Kdp (deg/km)'
      title(606) = 'Delta'
      title(607) = 'RHOhv'
      title(608) = 'LDR (dB)'
      title(609) = 'Hdr'
      title(610) = 'Hdw'
      title(611) = 'Zdr-rain (dB)'
      title(621) = 'Fi (Kdp)'
      title(622) = 'Fl (Kdp)'
      title(623) = 'Fi (Zdp)'
      title(624) = 'Fl (Zdp)'
c
      title(701) = 'HYDROMETEOR SUMMARY A'
      title(702) = 'HYDROMETEOR SUMMARY B'
      title(703) = 'HYDROMETEOR SUMMARY A: SCORE'
      title(704) = 'HYDROMETEOR SUMMARY B: SCORE'
      title(705) = 'HYDROMETEOR DETAILED A'
      title(706) = 'HYDROMETEOR DETAILED B'
      title(707) = 'HYDROMETEOR DETAILED A: SCORE'
      title(708) = 'HYDROMETEOR DETAILED B: SCORE'
c
      title(710) = 'HAIL TYPES'
      title(711) = 'Hail shape parameter' ! alphahl
      title(712) = 'Hail Content (g/m**3)'
      title(713) = 'Hail number concentration'
      title(714) = 'Hail Mean Volume Diameter'
      title(715) = 'HAIL TERMINAL VELCOITY'
      title(716) = 'Hail Maximum Mass Diameter' !'HAIL-SMALL: SCORE'
      title(717) = 'Hail Mass-weighted Diameter' ! 'HAIL-LARGE: SCORE'
      title(718) = 'HAIL-GIANT: SCORE'
      title(719) = 'HAIL-WET: SCORE'
      title(720) = 'HAIL-DRY: SCORE'
c
      title(721) = 'Frozen Drop Shape parameter'
      title(722) = 'Frozen Drop Content (g/m**3)'
      title(723) = 'Frozen Drop Number'
      title(724) = 'Frozen Drop Mean Volume Diameter'
      title(725) = 'Frozen Drop TERMINAL VELCOITY'
      title(726) = 'Frozen Drop Maximum Mass Diameter'
      title(727) = 'Frozen Drop Mass-weighted Diameter'
c
      title(730) = 'GRAUPEL TYPES'
      title(731) = 'Graupel Shape parameter'
      title(732) = 'Graupel Content (g/m**3)'
      title(733) = 'Graupel Number'
      title(734) = 'Graupel Mean Volume Diameter'
      title(735) = 'GRAUPEL TERMINAL VELCOITY'
      title(736) = 'Graupel Maximum Mass Diameter'
      title(737) = 'Graupel Mass-weighted Diameter'
c
      title(750) = 'RAIN TYPES'
      title(751) = 'Rain Shape parameter'
      title(752) = 'RAIN CONTENT (g/m**3)'
      title(753) = 'RAIN NUMBER'
      title(754) = 'RAIN Mean Volume Diameter'
      title(755) = 'RAIN TERMINAL VELCOITY'
      title(756) = 'RAIN Maximum Mass Diameter'
      title(757) = 'RAIN Mass-weighted Diameter'
      title(758) = 'RAIN Median Volume Diameter'
      title(759) = 'RAIN-LARGE-LOW CONC: SCORE'
c
      title(770) = 'SNOW TYPES'
      title(771) = 'SNOW TYPES: SCORE'
      title(772) = 'SNOW CONTENT'
      title(773) = 'SNOW NUMBER'
      title(774) = 'SNOW DIAMETER'
      title(775) = 'SNOW TERMINAL VELCOITY'
      title(776) = 'SNOW-WET AGGREGATE: SCORE'
      title(777) = 'SNOW-DRY AGGREGATE: SCORE'
      title(778) = 'SNOW-COLUMN: SCORE'
      title(779) = 'SNOW-NEEDLE: SCORE'
      title(780) = 'SNOW-PLATE: SCORE'
c
      title(790) = 'RAIN/HAIL TYPES'
      title(791) = 'RAIN/HAIL TYPES: SCORE'
      title(792) = 'RAIN/HAIL-SMALL: SCORE'
      title(793) = 'RAIN/HAIL-LARGE: SCORE'
      title(794) = 'RAIN/HAIL-GIANT: SCORE'

      
      title(901) = 'Total lightning segments'
      title(902) = 'All Positive segments'
      title(903) = 'All Negative segments'

      title(904) = 'CG: Total lightning segments'
      title(905) = 'CG: All Positive segments'
      title(906) = 'CG: All Negative segments'
      
      title(907) = 'IC: Total lightning segments'
      title(908) = 'IC: All Positive segments'
      title(909) = 'IC: All Negative segments'

      title(910) = 'Total lightning initiations'
      title(911) = 'All CG:  initiations'
      title(912) = 'IC:  initiations'
      title(913) = '-CG:  initiations'
      title(914) = '+CG:  initiations'
      
      title(919) = 'Lightning dipole/monopole'

c
c  Min/Max charge separation rates
c
      title(951) = 'cgizmn'
      title(952) = 'cgizmx'
      title(953) = 'cghwzmn'
      title(954) = 'cghwzmx'

c
c  Neg/pos integrated charge separation rates (C per minute)
c
      title(955) = 'ctghinz'
      title(956) = 'ctghipz'
      title(957) = 'ctghwnz'
      title(958) = 'ctghwpz'

      title(981) = 'Positive Charge'
      title(982) = 'Negative Charge'
      title(983) = 'Net Charge'
      title(984) = 'Positive Ion Charge'
      title(985) = 'Negative Ion Charge'
      title(986) = 'Net Ion Charge'

      title(991) = 'Max Reflectivity (dBZ)'
      title(992) = 'Integrated Reflectivity'
      title(993) = 'Updraft vol < -5 m/s'
      title(994) = 'Updraft vol >  5 m/s'
      title(995) = 'Updraft vol > 10 m/s'
      title(996) = 'Updraft vol > 20 m/s'
      
      title(1001) = 'UDMF kg/s'
      title(1002) = 'UICMF'
      title(1003) = 'Graupel Volume (km**3)'

      title(1004) = 'Hail Volume (km**3)'
      title(1005) = 'Graupel Mass (kg)'
      title(1006) = 'Hail Mass (kg)'
      
      title(1007) = 'Ice crystal Mass (kg)'
      title(1008) = 'Rain Mass (kg)'
      title(1009) = 'Snow Mass (kg)'

      title(1010) = 'IWC Max (g/m**3)'  ! Ice Water content
      title(1011) = 'W-Max (m/s)'
      title(1012) = 'W-Min (m/s)'

!
! Microphysics process time-height: the order is determined by output from getstat.procth.pl
!
      title(1021) = 'crfrzf'
      title(1022) = 'ciacrf'
      title(1023) = 'chcnsh'
      title(1024) = 'chcnih'
      title(1025) = 'qhacw+qhacr+qhshr'
      title(1026) = 'qracw'
      title(1027) = 'qrcnw'
      title(1028) = 'vhacw+vhacr+vhshr'
      title(1029) = 'total latent heating (excl. condensation)'
      title(1030) = 'Rime splinter'
      title(1031) = 'Drop Freeze splinters'
      title(1032) = 'qrfrzf'
      title(1033) = 'qiacrf'
      title(1034) = 'crcnw'
      title(1035) = 'rain drop prod by melting/shedding'
      title(1036) = 'Condensation/Evaporation Latent Heating (cloud+rain)'
      title(1037) = 'Condensation mass rate'
      title(1038) = 'Rain+cloud evaporation (pevap)'
      title(1039) = 'pmlt'
      title(1040) = 'pdep'
      title(1041) = 'psub'
      title(1042) = 'pfrz'
      title(1043) = 'Frozen drop production (crfrz+ciarcr)'
      title(1044) = 'qrfrzf+qiacrf'
      title(1045) = 'Rime Density'
      title(1046) = 'Graupel number conversion (chcnsh+chcnih)'
      title(1047) = 'All ice splinters (chmul1+csplinter)'
      title(1020+ibase+1) = 'Total latent heating'
      title(1020+ibase+2) = 'Total latent heating rate (recalc)'
      title(1020+ibase+3) = 'Frz/mlt latent heating rate (recalc)'
      title(1020+ibase+4) = 'Sub/dep latent heating rate (recalc)'
      title(1020+ibase+5) = 'Cond/evap latent heating rate (recalc)'
      title(1020+ibase+6) = 'Evap latent heating rate (recalc)'
      title(1020+ibase+7) = 'Evap+sub+melt latent cooling rate (recalc)'
      IF ( ibase > 27 ) THEN
      title(1048) = 'Frozen drops from rain-auto'
      title(1049) = 'Frozen drops from rain-shed'
      title(1050) = 'Frozen drops from rain-melt'
      title(1051) = 'Hail number from graupel' ! chlcnhhl
      title(1052) = 'Hail number from frozen drops' ! chlcnfhl
      title(1053) = 'Droplet contact freezing (#)' ! chlcnfhl
      title(1054) = 'Ice crystal nucleation (#)' ! chlcnfhl
      ENDIF


c
c individual hydrometeor time-height max value and integrated
c
c  Maximum mixing ratio:
      title(1101) = 'Cloud water (kg/kg) (max)'
      title(1102) = 'Rain (kg/kg) (max)'
      title(1103) = 'Ice Columns (kg/kg) (max)'

      title(1104) = 'Rimed ice (kg/kg) (max)'
      title(1105) = 'Snow (kg/kg) (max)'

      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1106) = 'Rain-Auto max (kg/kg)'
      title(1107) = 'Rain-Shed max (kg/kg)'
      title(1108) = 'Rain-Melt max (kg/kg)'
      ELSE
      title(1106) = 'Low Dens. Graupel (kg/kg) (max)'
      title(1107) = 'Med Dens. Graupel (kg/kg) (max)'
      title(1108) = 'High Dens. Graupel (kg/kg) (max)'
      ENDIF

      title(1109) = 'Frozen Drops (kg/kg) (max)'
      title(1110) = 'Small Hail (kg/kg) (max)'
      title(1111) = 'Ice Plates (kg/kg) (max)'
      title(1112) = 'Large Hail (kg/kg) (max)'
c
c Integrated mass (mixing ratio) per model level:
c
      title(1113) = 'Cloud water mass (kg)'
      title(1114) = 'Rain mass (kg)'
      title(1115) = 'Ice Columns mass (kg)'

      title(1116) = 'Rimed ice mass (kg)'
      title(1117) = 'Snow mass (kg)'

      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1118) = 'Rain-Auto mass (kg)'
      title(1119) = 'Rain-Shed mass (kg)'
      title(1120) = 'Rain-Melt mass (kg)'
      ELSE
      title(1118) = 'Low Dens. Graupel mass (kg)'
      title(1119) = 'Med Dens. Graupel mass (kg)'
      title(1120) = 'High Dens. Graupel mass (kg)'
      ENDIF

      title(1121) = 'Frozen Drops mass (kg)'

      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1122) = 'Graupel mass (kg)'
      ELSE
      title(1122) = 'Small Hail mass (kg)'
      ENDIF
      title(1123) = 'Ice Plates mass (kg)'
      
      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1124) = 'Hail mass (kg)'
      ELSE
      title(1124) = 'Large Hail mass (kg)'
      ENDIF
      title(1125) = 'Total hydrometeor mass (kg)'
      title(1126) = 'Total precipitation mass (kg)'
      title(1127) = 'Total non-precipitation mass (kg)'
      title(1128) = 'Graupel+FD/Hail (kg)'

c  Maximum content:

      title(1131) = 'Cloud water (kg/m**3) (max)'
      title(1132) = 'Rain (kg/m**3) (max)'
      title(1133) = 'Ice Columns (kg/m**3) (max)'

      title(1134) = 'Rimed ice (kg/m**3) (max)'
      title(1135) = 'Snow (kg/m**3) (max)'

      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1136) = 'Rain-Auto mass (kg/m**3) (max)'
      title(1137) = 'Rain-Shed mass (kg/m**3) (max)'
      title(1138) = 'Rain-Melt mass (kg/m**3) (max)'
      ELSE
      title(1136) = 'Low Dens. Graupel (kg/m**3) (max)'
      title(1137) = 'Med Dens. Graupel (kg/m**3) (max)'
      title(1138) = 'High Dens. Graupel (kg/m**3) (max)'
      ENDIF

      title(1139) = 'Frozen Drops (kg/m**3) (max)'

      title(1140) = 'Small Hail (kg/m**3) (max)'
      title(1141) = 'Ice Plates (kg/m**3) (max)'
      title(1142) = 'Large Hail (kg/m**3) (max)'

! Mass per level (as 1113-1127 but multiplied by air density)

      title(1143) = 'Cloud water (kg/level)'
      title(1144) = 'Rain (kg/level)'
      title(1145) = 'Ice Columns (kg/level)'

      title(1146) = 'Rimed ice (kg/level)'
      title(1147) = 'Snow (kg/level)'

      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1148) = 'Rain-Auto mass (kg/level)'
      title(1149) = 'Rain-Shed mass (kg/level)'
      title(1150) = 'Rain-Melt mass (kg/level)'
      ELSE
      title(1148) = 'Low Dens. Graupel (kg/level)'
      title(1149) = 'Med Dens. Graupel (kg/level)'
      title(1150) = 'High Dens. Graupel (kg/level)'
      ENDIF

      title(1151) = 'Frozen Drops (kg/level)'


      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1152) = 'Graupel (kg/level)'
      ELSE
      title(1152) = 'Small Hail (kg/level)'
      ENDIF


      title(1153) = 'Ice Plates (kg/level)'

      title(1154) = 'Large Hail (kg/level)'
      IF ( microp(1:1) .eq. 'Z' ) THEN
      title(1154) = 'Hail (kg/level)'
      ELSE
      title(1154) = 'Large Hail (kg/level)'
      ENDIF

      title(1158) = 'Graupel+FD/Hail (kg/level)'

      RETURN
      END

c
c #################################################################
c
      SUBROUTINE CPMPXY(IMAP,XINP,YINP,XOTP,YOTP)
C 
C Transform contours to overlay various mapping transformations:
C IMAP= 0 - Cartesian data: no transformation necessary
C IMAP= 1 - Lat/Lon transformation
C IMAP=-1 - inverse Lat/Lon transformation
C IMAP= 2 - Rho/Theta transformation
C IMAP=-2 - inverse Rho/Theta transformation
C IMAP= 3 - stretched height 
C 
C 
      implicit none
      integer imap
      real xinp,yinp,xotp,yotp
      integer nzmax,nzmax1,nz1d
      parameter(nzmax1=400)
      real z1d2(nzmax1,4)
      real dlz
      integer iy
      real z1,z2,zr,dify,y
      integer kstrzb
      real hgtlow
      common/strzb/kstrzb,hgtlow
      common/strz/z1d2,dlz,nz1d
      
C 
C Handle the EZMAP case ...
C 
      IF (ABS(IMAP).EQ.1) THEN
         IF (IMAP.GT.0) THEN
            CALL MAPTRA (YINP,XINP,XOTP,YOTP)
         ELSE
            CALL MAPTRI (XINP,YINP,YOTP,XOTP)
         END IF
C 
C ... the polar coordinate case ...
C 
      ELSE IF (ABS(IMAP).EQ.2) THEN
         IF (IMAP.GT.0) THEN
            XOTP=XINP*COS(.017453292519943*YINP)
            YOTP=XINP*SIN(.017453292519943*YINP)
         ELSE
            XOTP=SQRT(XINP*XINP+YINP*YINP)
            YOTP=57.2957795130823*ATAN2(YINP,XINP)
         END IF
C 
C grid stretched in z-direction
C 
C Pressure transformation in the X direction
      ELSEIF(IMAP.EQ.3) THEN
C
C The height transformation in X direction is linear
C
         XOTP = XINP
C
C Find next lowest y data point
C
         Y = YINP + kstrzb
C
C Distance between next lowest data point and contour point
C
         IY=INT(Y)
c      IF ( kstrzb .gt. 0 .and. (xinp.gt.0.5 .and. xinp.lt.0.6 )) THEN
c         write(6,*) 'y,iy,xinp',yinp,iy,xinp! ,z1d2(iy,1)
c         ENDIF
         IF ( y .lt. 1 .and. (xinp.gt.0.5 .and. xinp.lt.0.51 )) THEN
         write(6,*) 'y,iy,xinp',yinp,iy,xinp,y*41. ! ,z1d2(iy,1)
         ENDIF
c         IF ( iy .eq. 0 ) THEN
c           y = yinp*41 + 1
c           iy = int(y)
c         ENDIF
c         IF ( IY .eq. 0 ) IY = 1
         
         DIFY=Y-FLOAT(IY)
C
C Linear interpolation between points to give height at contour point
C
         Z1 = z1d2(iy,1)
         Z2 = z1d2(iy+1,1)
         IF ( kstrzb .eq. 0 ) THEN
         ZR = Z1 + DIFY*(Z2-Z1)
         ELSE
         ZR = Z1 + DIFY*(Z2-Z1) - hgtlow
         ENDIF 
         YOTP=0.001*zr
c         write(6,*) 'yinp,iy,yotp',yinp,iy,yotp,z1d2(iy,3)
C 
C If IMAP isnt specified as above, then do an identity transformation.
C 
      ELSE
         XOTP = XINP
         YOTP = YINP
      ENDIF
C 
      RETURN
      END

c #####################################################################
      SUBROUTINE VVUMXY (X,Y,U,V,UVM,XB,YB,XE,YE,IST)
C
C This is a user modifiable routine that allows custom projections of
C the vector space. X and Y give the vector position within the domain
C of the data space. By default, this space is coincident with the
C grid space (i.e. 1 through dimension lengths of the U and V arrays).
C The vector endpoints are output in fractional coordinates (NDC space).
C Note that this is different from the old MXF and MYF routines, which
C output in 'plotter coordinate' space. It also differs from the 
C Conpack routine CPMPXY, which returns values in user space. 
C 
C VVUMXY (Velocity Vector -- User Map X,Y) is called whenever 
C the internal parameter MAP is set to a value other than 0, 1, or 2.
C
C Based on the magnitude and direction of the vector the start and 
C ending points of the vector are returned in NDC space.
C
C Input parameters:
C
C X,Y   -- vector position in the user coordinate system
C U,V   -- vector components from the U,V arrays for this position
C UVM   -- magnitude of the U,V components (supplied for convenience
C          and efficiency - but note that many mappings do not need 
C          this value)
C
C Output parameters:
C
C XB,YB -- starting point of the vector in fractional coordinates
C          (NDC space), before offset based on the value of the
C          vector positioning parameter, VPO.
C XE,YE -- ending point of the vector in fractional coordinates
C          (NDC space), before offset based on the value of the
C          vector positioning parameter, VPO.
C IST   -- status results of the mapping: 0 indicates success -- any
C          non-zero value causes VVECTR to discard the vector at this
C          location.
C
C
C The mapping common block: made available to user mapping routines
C
      integer nzmax,nzmax1,nz1d
      parameter(nzmax1=400)
      real z1d2(nzmax1,4)
      real dlz
      integer iy
      real z1,z2,zr
      integer icall
      data icall/0/
      common/strzb/kstrzb,hgtlow
      common/strz/z1d2,dlz,nz1d

      COMMON /VVMAP/
     +                IMAP       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                SXDC       ,SYDC       ,NXCT       ,NYCT       ,
     +                RLEN       ,LNLG       ,INVX       ,INVY       ,
     +                ITRT       ,IWCT       ,FW2W       ,FH2H       ,
     +                DVMN       ,DVMX       ,RBIG       ,IBIG
C
      SAVE /VVMAP/
C
C Math constants
C
      PARAMETER (PDTOR  = 0.017453292519943,
     +           PRTOD  = 57.2957795130823,
     +           P1XPI  = 3.14159265358979,
     +           P2XPI  = 6.28318530717959,
     +           P1D2PI = 1.57079632679489,
     +           P5D2PI = 7.85398163397448) 
C
C --------------------------------------------------------------------
C
C Description of VVMAP contents:
C
C IMAP                - value of the internal parameter 'MAP'
C XVPL,XVPR,YVPB,YVPT - the currently set viewport values. (GETSET
C                       arguments 1, 2, 3, and 4)
C WXMN,WXMX,WYMN,WYMX - the min and max boundaries of user coordinate
C                       space, (usually but not always equivalent to
C                       window coordinates). WXMN and WYMN are true
C                       minimum values even one or both axes is 
C                       inverted. (i.e. WXMN = MIN of GETSET args 5,6;
C                       WYMN = MIN of GETSET args 7,8)
C XLOV,XHIV,YLOV,YHIV - min and max boundaries of the data space, by
C                       default equivalent to the array grid space.
C                       XLOV and YLOV are not necessarily less than 
C                       XHIV and YHIV respectively.
C SXDC,SYDC           - Scaling factors for converting vector component
C                       values into lengths in NDC space.
C NXCT,NYCT           - Length of each dimension of the U and V 
C                       component arrays.
C RLEN                - Length of the maximum vector in user 
C                       coordinates.
C LNLG                - The linear/log mode (GETSET argument 9)
C INVX,INVY           - User coordinates inversion flags: 
C                       0 - not inverted, 1 - inverted
C ITRT                - value of the internal parameter TRT
C IWCT                - not currently used
C FW2W,FH2H           - scale factors for converting from fraction of
C                       viewport width/height to NDC width/height 
C DVMN,DVMX           - min/max vector lengths in NDC
C RBIG,IBIG           - machine dependent maximum REAL/INTEGER values
C
C --------------------------------------------------------------------
C User-defined common block...
C
C --------------------------------------------------------------------
C
C     IF (IMAP .EQ. =user defined mapping code=) THEN
C
C User defined mapping code....
C
C     ELSE
      
!      write(6,*) 'IMAP = ',imap
      
      IF (IMAP .EQ. 3) THEN
        dxvpt = XVPR-XVPL  ! ndc width
        dyvpt = YVPT-YVPB  ! ndc height
        
        dxndc = WXMX-WXMN  ! user width
        dyndc = WYMX-WYMN  ! user height
        
        IF ( icall .eq. 0 ) THEN
        icall = 1
        write(6,*) 'dxvpt,dyvpt,dxndc,dyndc',dxvpt,dyvpt,dxndc,dyndc
        write(6,*) 'XVPL,XVPR,YVPB,YVPT',XVPL,XVPR,YVPB,YVPT
        write(6,*) 'WXMN,WXMX,WYMN,WYMX',WXMN,WXMX,WYMN,WYMX
        write(6,*) 'X,Y,U,V,UVM,XB,YB,XE,YE'
        
        ENDIF
        
        xb = x
C
C Distance between next lowest data point and contour point
C
c         YINP = Y + kstrzb
         YINP = Y
         IY=INT(YINP)
         DIFY=YINP-FLOAT(IY)
C
C Linear interpolation between points to give height at contour point
C
         IF ( kstrzb .gt. 0 ) THEN
           z0 = z1d2(kstrzb,1)
         ELSE
           z0 = 0.0
         ENDIF
         IF ( iy+kstrzb .ge. 1 ) THEN
           Z1 = z1d2(iy+kstrzb,1)
           Z2 = z1d2(iy+kstrzb+1,1)
         ELSE 
           z1 = 0.0
           z2 = 0.0
         ENDIF
c         write(6,*) 'yinp,z1,z2,kstrzb,u,v',yinp,z1,z2,kstrzb,u,v
         ZR = Z1 + DIFY*(Z2-Z1)
c         IF ( kstrzb .eq. 0 ) THEN
c         ZR = Z1 + DIFY*(Z2-Z1)
c         ELSE
c         ZR = Z1 + DIFY*(Z2-Z1) !  - hgtlow*1000.0
c         ENDIF 
c         zt = z1d2(Int(WYMX-1),1)  
cc     :    +  0.5*(z1d2(Int(WYMX-1),1) - z1d2(Int(WYMX-2),1))

c         yb = YVPB + dyvpt*zr/zt
         
         zr = 0.001*zr
         yb = cufy(zr)
         
c         xc = 0.0
c         yc = 0.0
         
         
         xe = xb + SXDC*U
         ye = yb + SYDC*V

!        write(6,*) X,Y,U,V,UVM,XB,YB,XE,YE
        
      ELSE
C
C Default mapping:
C
C Note that the default version of VVUMXY included with
C NCAR Graphics 3.2 incorrectly fails to convert the input X and Y
C location into fractional coordinates (NDC space).
C
         IF (X.LT.WXMN .OR. X.GT.WXMX .OR. 
     +        Y.LT.WYMN .OR. Y.GT.WYMX) THEN
            IST = -1
            RETURN
         END IF
         XB=CUFX(X)
         YB=CUFY(Y)
         XE=XB+U*SXDC
         YE=YB+V*SYDC

      END IF
C
C Done.
C
      RETURN
C
      END



