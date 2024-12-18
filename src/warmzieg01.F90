!#include "sam.def.h"
!#define ICE10
#define SWM
#undef ELEC
!#define SAM
!
! Things to do:
!
!  Fix Rain evaporation for gamma function (ipconc >= 3)
!
!  convert cloud ice to snow as in Ferrier 1994 (change only mass in cloud ice),
!    then can try turning off direct conversion from cloud ice to graupel and rimed ice
!
!  look at an iterative check on overdepletion;  need to be careful with two-moment
!
!  check ice supersaturation in two-moment.  Getting enough deposition, or need 
!      to do sat adj. when cloud droplets are all gone?
!
!  
!
! new comment
!
! Fix use of gt for SWM IN FALLOUT ROUTINES
!
!  How to remove hl for ipconc=5?  Need to preprocess?
!
!   When the charging rates are moved to a subroutine, need to move the
!   call to be after the wet growth calculations -- or at least the 
!   splashing stuff.  Think about this....
!
!  Think about what to do with cracif
!
!    FIX BIGSAVE! (done)
!
!    Replace qv0 with qx(mgs,lv)? No. qv0 is base val
!
! Need to look at limiting supersaturation to 1 or so by nucleation/condensation
!
!  put in temperature-dependent function for homogeneous freezing
!
!c--------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------
!
      subroutine warmzieg &
     &  (ntmul,nx,ny,nz,na,nba,nv,nstep &
     &  ,nor,istag,jstag,kstag,itopo &
     &  ,iwrite &
     &  ,cwmasn,cwmasx &
     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d &
     &  ,xfall &
     &  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9 &
     &  ,ab,an,db,dn,p2,pinit &
     &  ,pb,pn,u,v,w,iunit &
     &  ,ssat,ssfilt,t00,t77,ssati,thproc,numproc,axtra)
!     >  ,elec,neelec)

       USE MICRO_MODULE
       USE INDEX_MODULE, only: lt,lc,lr,li,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lhw,lhlw, &
                               lnc,lnr,lni,lns,lnh,lnhl,cinu,dmuh,dnu,dmu,xnu,xmu,dmuhl,rnu,cnu,snu, &
                               lss,lsat,lsati,xvcmx,xvcmn,xvrmn,xvrmx,lccn,rnumin,rnumax, &
                               lqmx,nxtra,xcradmx
       USE CPUTIME_MODULE
       USE COMMASMPI_MODULE
       USE TRAJ_MODULE
!
!--------------------------------------------------------------------------
!                                
!     Ziegler 1985 parameterized microphysics (also Zrnic et al. 1993)
!     1)  cloud water
!     2)  rain
!     3)  column ice 
!     6)  snow
!     11) graupel/hail
!
!--------------------------------------------------------------------------
!
! Notes:
!
!  5/12/2006:  Converted qsacw/csacw and qsaci/csaci to Z93
!
!  5/12/2006:  Put a threshold on Bigg rain freezing.  If the frozen drops
!              have an average volume less than xvhmn, then the drops are put
!              into snow instead of graupel/hail.
!
!              Fixed bug when vapor deposition was limited.
!
!  5/13/2006:  Note that qhacr has a large effect, but Z85 did not include it.
!              Turned off qsacr (set to zero).
!
!--------------------------------------------------------------------------
!
!  general declarations
!
!--------------------------------------------------------------------------
!
!
!
      implicit none
!
      logical, parameter :: fallonly = .false.
      
      integer ng1
      integer iunit
      parameter(ng1 = 1)
      
      real,allocatable :: db1(:,:),db0(:,:),ztmp(:,:,:),dtz1(:,:),dtz0(:,:)
      real qvex
      
      integer iraincv, icgxconv
      parameter ( iraincv = 1, icgxconv = 1)
      real ffrz

      real ccwtmp
      
      double precision dp1
      
      real frac
      
      real dtp1
      integer ntmul
      
!      real rar  ! rime accretion rate as calculated from qxacw


! a few vars for time-split fallout      
      real vtmax
      integer n,ndfall
      
      integer nstep
      integer nx,ny,nz,na,nba,nv
      integer ng
      integer nor,imapz,mzdist,istag,jstag,kstag,itopo
      integer iwrite
      real dtp,dx,dy,dz
      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)

      integer numproc
      real thproc(nzend,numproc)

      real dv

      real dtptmp
!      integer nxl,nyl,nzl

      integer itest,nidx,id1,jd1,kd1
      parameter (itest=1)
      parameter (nidx=10)
      parameter (id1=1,jd1=1,kd1=1)
      integer ierr
      integer iend,imake,iread
      save imake
      integer ix,jy,kz, il,in, ic, ir, icp1, irp1
      integer i,j
      real slope1, slope2
      real x, x1, x2, x3, del, y
!      integer nxm,nym,nzm
      real eps,eps2
      parameter (eps=1.e-20,eps2=1.e-5)
!
!  electrical permitivity of air C / (N m**2) -  check the units
!
      real eperao
      parameter (eperao  = 8.8592e-12 )
      
      real ec,eci  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      parameter (eci = 1.0/ec)
!
!      include 'swm.index.warmzieg.h'
!
      
      

      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      integer lz(lc:lhab)
      integer lvol(lc:lhab)
      
!
! temporary arrays-self contained-sizes
!
!
!  conversion parameters
!
      integer mqcw,mqxw,mtem,mrho,mtim
      parameter (mqcw=21,mqxw=21,mtem=21,mrho=5,mtim=6)

      real xftim,xftimi,yftim, xftem,yftem, xfqcw,yfqcw, xfqxw,yfqxw
      parameter (xftim=0.05,xftimi = 1./xftim,yftim=1.)
      parameter (xftem=0.5,yftem=1.)
      parameter (xfqcw=2000.,yfqcw=1.)
      parameter (xfqxw=2000.,yfqxw=1.)
      
!
!  charge fallout arrays
!
      real xfall(nx,ny,na)
      real xfall0(nx,ny)
!      real gt0(-nor+ng1:nx+nor,-nor+ng1:1+nor,-nor+ng1:nz+nor,ngt)

!
! params read in from inmicro
!

      real dtfac
      parameter ( dtfac = 1.0 )

      

!      real dtrim

      


      real cckm,ccne,ccnefac,cnexp
      save cckm,ccne,ccnefac,cnexp
      real beta
!
!
      

!
! external temporary arrays
!
      real t00(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t77(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

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
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,nxtra)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)

      real u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      real ssati(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real ssat(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ssfilt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
! 
!  declarations microphyscs and for gather/scatter
!
      integer nxmpb,nzmpb,nxz
      integer jgs,mgs,ngs,numgs
      parameter (ngs=500)
      integer ntt
      parameter (ntt=300)

      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs)
      integer nsvcnt
!      integer isave(ntt)
      integer ncuse
      parameter (ncuse=0)
      integer il0(ngs),il5(ngs)
!      integer il1m(ngs),il2m(ngs),il3m(ngs),il4m(ngs),il5m(ngs)
!      integer nqsat
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati
!      parameter (fqsat=0.002,fqsati=1./fqsat)
!      parameter (fqsat=0.01,fqsati=1./fqsat)
!
!
!
      real :: axx(ngs,lh:lhab), bxx(ngs,lh:lhab)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)

!      real cai,caw,cbi,cbw
      real tdtol,tfrcbw,tfrcbi,thnuc
      
      real tfr,tfrh
      parameter ( tfr = 273.15, tfrh = 233.15)
      
      real cp, rd
      parameter ( cp = 1004.0, rd = 287.04 )
      
      real cpi
      parameter ( cpi = 1./cp )
      
      real poo,cap
      parameter ( cap = rd/cp, poo = 1.0e+05 )

!      real tmxs(ntt),xmxs(ntt),xmns(ntt)
!
!
!  gamma function
!
!      integer ngm0,ngm1,ngm2
!      parameter (ngm0=1000,ngm1=500,ngm2=500)
!      real gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)
!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs)
      real cnuc(ngs)
      real sscb  ! 'cloud base' SS threshold
      parameter ( sscb = 2.0 )
      integer idecss  ! flag to turn on (=1) decay of ssmax when no cloud or ice crystals
      parameter ( idecss = 1 )
      integer iba ! flag to do condensation/nucleation in 1st or 2nd loop
                  ! =0 to use ad to calculate SS
                  ! =1 to use an at end of main jy loop to calculate SS
      parameter (iba = 1)
      integer ifilt   ! =1 to filter ssat, =0 to set ssfilt=ssat
      parameter ( ifilt = 0 ) 
      real temp1,temp2
      real ssmax(ngs)       ! maximum SS experienced by a parcel
      real ssmx
      real bfnu
      parameter ( bfnu = (rnu + 2.0)/(rnu + 1.0)  )
      real ventr, ventc
      real ventrx(ngs)
      save ventr, ventc
      real volb, aa1, aa2
      double precision t2s, xdp
      double precision xl2p(ngs),rb(ngs)
      parameter ( aa1 = 9.44e15, aa2 = 5.78e3 ) ! a1 in Ziegler
! snow parameters:
      real cexs, cecs
      parameter ( cexs = 0.1, cecs = 0.5 )
      real rvt      ! ratio of collection kernels (Zrnic et al, 1993)
      parameter ( rvt = 0.104 )
      real cautn(ngs), rh(ngs), nh(ngs)
      real rhoinv(ngs)
      double precision ec0(ngs)
      
      integer kbound
      real ac1,bc, taus, c1,d1,e1,f1,p380
      real d1r
      real c1sw   ! integration factor for snow melting with snu = -0.8
      save c1sw
      real, parameter :: vr1mm = 5.23599e-10 ! volume of 1mm diameter sphere (m**3)
      real rhosm
      parameter ( rhosm = 500. )
      integer nc ! condensation step
      real dtcon,dtcon1,dtcon2 ! condensation time step (dtcon*nc = dtp)
      real delta
      integer ltemq1,ltemq1m ! ,ltemq1m2
      real dqv,qv1,ss1,ss2,qvs1,dqvs,dtemp,dt1   ! temporaries for condensation
!      real  dtemp2,ss1m2
      real ssi1, ssi2, dqvis, dqvii,qis1
      real dqvr, dqc, dqr
      real qv1m,qvs1m,ss1m
      real cwmastmp 
      real dcloud,dcloud2
      real cn(ngs) 
!      real xvc(ngs), xvr(ngs)
!      real xvs(ngs),xvgl(ngs),xvgm(ngs),xvgh(ngs),xvf(ngs)
!      real xvh(ngs),xvhl(ngs)
      real,allocatable :: xv(:,:)
      real mwfac
!      parameter ( mwfac = 6.0**(1./3.) ) 
!      ! factor for mass-weighted rain volume diameter
!      real wijk ! wvel
      real  es(ngs) ! ss(ngs),
      real ssbar(ngs) !, sqsat(ngs)
      real ssf(ngs),ssfkp1(ngs),ssfkm1(ngs),ssat0(ngs)
      real ssfjp1(ngs),ssfjm1(ngs)
      real ssfip1(ngs),ssfim1(ngs)
      real supcb, supmx
      parameter (supcb=0.5,supmx=238.0)
      real r2dxm, r2dym, r2dzm
      real dssdz, dssdy, dssdx
!      real cck
!      real xvcmn, xvcmx  ! min, max droplet volumes
!      real xvrmn, xvrmx  ! min, max rain volumes
!      real xvsmn, xvsmx  ! min, max snow volumes
!      real xvfmn, xvfmx  ! min, max frozen drop volumes
!      real xvgmn, xvgmx  ! min, max graupel volumes
!      real xvhmn, xvhmx  ! min, max hail volumes
!      real xvhlmn, xvhlmx  ! min, max lg hail volumes
      real xvmn(lc:lhab), xvmx(lc:lhab)
!      PARAMETER(XVCMN=4.188E-12,XVCMX=6.54E-8)   ! CGS
!      PARAMETER(XVRMN=2.8866E-7,XVRMX=4.1887E-3) ! CGS
!      parameter( xvcmn=4.188e-18 )   ! mks  min volume = 1 micron radius
!      parameter( xvcmx=6.54e-14 )    ! mks, 25 micron max radius
!      parameter( xvcmx=2.89e-13 )    ! mks, 41 micron max radius

!      data xvcmx/2.89e-13/
!      save xvcmx
!      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(3.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvsmn=0.523599*(0.1e-3)**3, xvsmx=0.523599*(6.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvfmn=0.523599*(0.1e-3)**3, xvfmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvgmn=0.523599*(0.1e-3)**3, xvgmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhmn=0.523599*(0.3e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhlmn=0.523599*(1.e-3)**3, 
!     :           xvhlmx=0.523599*(40.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      real rwmasn,rwmasx

      real vgra
      parameter ( vgra = 0.523599*(1.0e-3)**3 )
     
      real epsi,d
      parameter (epsi = 0.622, d = 0.266)
      real r1,qevap ! ,slv
      
      real vr,nrx,chw,g1,qr,z,z1,rdi,alp,xnutmp,xnuc
      real tmp, ctmp, qtmp

      integer infdo
      
!      real svc(ngs)  !  droplet volume
!
!  misc
!
      
!      real tabqvs(nqsat),tabqis(nqsat),dtabqvs(nqsat)
!      save tabqvs,tabqis,dtabqvs
      real temp(ngs),tempc(ngs),tempbz(ngs),temp0(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) ! ,tembzg(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real elv(ngs),elf(ngs),els(ngs),elvs(ngs),elss(ngs)
      real gamw(ngs),gams(ngs)   !   qciavl(ngs),
      real tsqr(ngs),ssi(ngs),ssw(ngs)
      real cc3(ngs),cqv1(ngs),cqv2(ngs)
      real f5, prod0, qvs0  ! Kessler condensation factor
      real cwmcnd(ngs)
      real cwc1,cwc0
!      real qcwdif(ngs) ! ,dcwnc
      real cwncmx(ngs),cwncmn(ngs),cwnccn(ngs)
      real cwmasn,cwmasx
      real cwmasn5
      real cwradn
      real cimcnd(ngs) ! ,qcidif(ngs)  ! ,dcinc
      real cincmx(ngs),cincmn(ngs)
!      real cinccn(nz)
!      real cinc(ngs)    !  ,qcitmp(ngs)
      real cimasn,cimasx,ccimx
      real pi,pid4
      real ar,br,cs,ds,gf7,gf6,gf5,gf4,gf3,gf2,gf1
      real gf73rds, gf83rds
      real aradcw,bradcw,cradcw,dradcw,cwrad,rwrad,rwradmn
      parameter ( rwradmn = 50.e-6 )
      
!
!  other arrays
!
      
      
      real fvds(ngs),fvce(ngs)
      real fvent(ngs)
!
      real fai(ngs),fav(ngs),fbi(ngs),fbv(ngs)
      real felv(ngs),fels(ngs),felf(ngs)
      real felvs(ngs),felss(ngs)      !   ,felfs(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs)
      real fschm(ngs),fpndl(ngs)
      real fgamw(ngs),fgams(ngs)
      real fcqv1(ngs),fcqv2(ngs),fcc3(ngs)  
!
!
      
      real qxmin(lc:lqmx)
      save qxmin
      
      real,allocatable :: qx(:,:) ! qx(ngs,lv:lhab)
      real,allocatable :: qxw(:,:) ! qx(ngs,lv:lhab)
      real,allocatable :: cx(:,:) ! cx(ngs,lc:lhab)
      real,allocatable :: cxmxd(:,:) ! cx(ngs,lc:lhab)
      real,allocatable :: qxmxd(:,:) ! cx(ngs,lc:lhab)
      real,allocatable :: alpha(:,:)  ! alpha(ngs,lc:lhab)
      real,allocatable :: g1x(:,:)  ! alpha(ngs,lc:lhab)
      real,allocatable :: zx(:,:) ! qx(ngs,lv:lhab)
      real,allocatable :: zxmxd(:,:) ! cx(ngs,lc:lhab)

      real g1shr, alphashr
      real g1mlr, alphamlr
!
      real rwvent(ngs)
!
!      real xivent(ngs,nhab),xix(ngs,nhab),xirey(ngs,nhab)
!      real xpvent(ngs)
!
!
!      real amccw(ngs) !,amcci(ngs),amcrw(ngs),amcsw(ngs),amchw(ngs)
!
!
      
      real, allocatable :: vtxbar(:,:,:)
      real, allocatable :: xmas(:,:)
      real, allocatable :: xdn(:,:)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real, allocatable :: xdia(:,:,:)
!      real vtxbar(ngs,nhab),ximas(ngs,nhab),xidn(ngs,nhab)
!
!
!
      real rwcap(ngs)

!
      real qcmxd(ngs),qrmxd(ngs)
!
      real ccmxd(ngs),crmxd(ngs)

!
!
!
!
!
!  concentration arrays...
!
!
!      real cwacigl(ngs) !, cwacigm(ngs), cwacigh(ngs)
!      real cwacipgl(ngs) ! , cwacipgm(ngs), cwacipgh(ngs)

!      real chcngh(ngs), chcngm(ngs), chcngl(ngs), chcnf(ngs)
!      real chlcnh(ngs)
!      real cracif(ngs), ciacrf(ngs)
!      real cracipf(ngs), cipacrf(ngs)
!      real  cipacip(ngs) !  
      real cracr(ngs)

!
!      real ciint(ngs), crfrz(ngs), crfrzf(ngs), crfrzs(ngs)
!      real cicint(ngs)  !  , ciracir(ngs), ciaci(ngs)
!      real ciacir(ngs), ciraci(ngs)
!      real cipint(ngs) !, cipacwi(ngs)
!      real ciracip(ngs)  !  , ciihrp(ngs)
!
!      real ciacw(ngs), cwaci(ngs), cwacii(ngs)
!      real ciacr(ngs), craci(ngs)
!      real ciracw(ngs), cwacir(ngs)
!      real ciracr(ngs)
!      real cracir(ngs)
!      real cipacw(ngs), cwacip(ngs)
!      real cipacr(ngs), 
!      real cracip(ngs)
!      real cipaci(ngs)
!
!      real csacw(ngs),   cwacs(ngs)
!      real csacwgl(ngs)
!      real cwacsgl(ngs)
!      real csacwgm(ngs)
!      real cwacsgm(ngs)
!      real csacwgh(ngs)
!      real cwacsgh(ngs)
!      real csacwf(ngs)
!      real cwacsf(ngs)
!      real csacr(ngs),   cracs(ngs)
!      real csacrgl(ngs)
!      real cracsgl(ngs)
!      real csacrgm(ngs)
!      real cracsgm(ngs)
!      real csacrgh(ngs)
!      real cracsgh(ngs)
!      real csacrf(ngs),  cracsf(ngs)
!      real csaci(ngs),   csacs(ngs)
!      real csacir(ngs),  ciracs(ngs)
!      real csacip(ngs),  cipacs(ngs)
!      real ciacs(ngs)
!
!
      real cracw(ngs)

      real crcnw(ngs) ! ,ciacwi(ngs)

      real crcev(ngs) ! ,crmlr(ngs)
!      real crshr(ngs)
!
!
! arrays for w-ac-x ;  x-ac-w
!
!
!
      real qrcnw(ngs)
      real zrcnw(ngs),zracr(ngs),zracw(ngs),zrcev(ngs)
!

      real qracw(ngs) ! qwacr(ngs),
!
!
!
!
!
!  conversions
!

      real cninm(ngs),cnina(ngs),cninp(ngs),wvel(ngs),wvelkm1(ngs)
      real uvel(ngs),vvel(ngs)
!
!
!
!
!      real qrztot(ngs),qrzmax(ngs),qrzfac(ngs)
      real qrcev(ngs)
!      real qrshr(ngs)
!
!
!      real exwidia(nhab),exwwdia(nhab)

!      real eiw(ngs),eii(ngs),eiri(ngs),eipir(ngs) ! eww(ngs),
      real erw(ngs) ! ,esw(ngs),eglw(ngs),eghw(ngs),efw(ngs)
!      real ehxw(ngs),ehlw(ngs),egmw(ngs),ehw(ngs) ! eaw(ngs),
      real err(ngs) ! ,esr(ngs),eglr(ngs),eghr(ngs),efr(ngs)
!      real ehxr(ngs),ehlr(ngs),egmr(ngs) ! ,eipr(ngs),ear(ngs)
!      real eri(ngs),esi(ngs),egli(ngs),eghi(ngs),efi(ngs)
!      real ehxi(ngs),ehli(ngs),egmi(ngs),ehi(ngs) ! eai(ngs),
!      real ers(ngs),ess(ngs),egls(ngs),eghs(ngs),efs(ngs),ehs(ngs)
!      real ehscnv(ngs)
!      real ehxs(ngs),ehls(ngs),egms(ngs),egmip(ngs) ! eas(ngs),

      real ew(8,6)
      real cwr(8,2)  ! radius and inverse of interval
      data cwr / 2.0, 3.0, 4.0, 6.0,  8.0,  10.0, 15.0,  20.0 , & ! radius
     &           1.0, 1.0, 0.5, 0.5,  0.5,   0.2,  0.2,  1.  /   ! inverse of interval
      integer icwr(ngs), igwr(ngs), irwr(ngs), ihlr(ngs)
      real grad(6,2) ! graupel radius and inverse of interval
      data grad / 100., 200., 300., 400., 600., 1000.,   &
     &            1.e-2,1.e-2,1.e-2,5.e-3,2.5e-3, 1.    /
!droplet radius: 2     3     4     6     8    10    15    20
      data ew /0.03, 0.07, 0.17, 0.41, 0.58, 0.69, 0.82, 0.88,  & ! 100
!     :         0.07, 0.13, 0.27, 0.48, 0.65, 0.73, 0.84, 0.91,  ! 150
     &         0.10, 0.20, 0.34, 0.58, 0.70, 0.78, 0.88, 0.92,  & ! 200
     &         0.15, 0.31, 0.44, 0.65, 0.75, 0.83, 0.96, 0.91,  & ! 300
     &         0.17, 0.37, 0.50, 0.70, 0.81, 0.87, 0.93, 0.96,  & ! 400
     &         0.17, 0.40, 0.54, 0.71, 0.83, 0.88, 0.94, 0.98,  & ! 600
     &         0.15, 0.37, 0.52, 0.74, 0.82, 0.88, 0.94, 0.98 / ! 1000
!     :         0.11, 0.34, 0.49, 0.71, 0.83, 0.88, 0.94, 0.95 / ! 1400

      

!      real ehip(ngs),ehlip(ngs),ehlir(ngs)
!      real erir(ngs),esir(ngs),eglir(ngs),egmir(ngs),eghir(ngs)
!      real efir(ngs),ehir(ngs),eirw(ngs),eirir(ngs),ehr(ngs)
!      real erip(ngs),esip(ngs),eglip(ngs),eghip(ngs)
!      real efip(ngs),eipi(ngs),eipw(ngs),eipip(ngs)
!
!  arrays for production terms
!
      real ptotal(ngs) ! , pqtot(ngs) 
     
!
      real pqcwi(ngs),pqrwi(ngs),pqwvi(ngs)
!      real pqswi(ngs),pqhwi(ngs),pqcii(ngs)
!      real pqgli(ngs),pqghi(ngs),pqfwi(ngs)
!      real pqgmi(ngs),pqhli(ngs) ! ,pqhxi(ngs)
!      real pqiri(ngs),pqipi(ngs) ! pqwai(ngs),
!
      real pqcwd(ngs),pqrwd(ngs),pqwvd(ngs)
!      real pqswd(ngs),pqhwd(ngs),pqcid(ngs)
!      real pqgld(ngs),pqghd(ngs),pqfwd(ngs)
!      real pqgmd(ngs),pqhld(ngs) ! ,pqhxd(ngs)
!      real pqird(ngs),pqipd(ngs) ! pqwad(ngs),
!
!      real pqxii(ngs,nhab),pqxid(ngs,nhab)
!
      real  pctot(ngs)
!      real  pcipi(ngs), pcipd(ngs)
!      real  pciri(ngs), pcird(ngs)
      real  pccwi(ngs), pccwd(ngs)
!      real  pccii(ngs), pccid(ngs)
      real  pcrwi(ngs), pcrwd(ngs)

      real  pzrwi(ngs), pzrwd(ngs)

!      real  pcswi(ngs), pcswd(ngs)
!      real  pchwi(ngs), pchwd(ngs)
!      real  pchli(ngs), pchld(ngs)
!      real  pcfwi(ngs), pcfwd(ngs)
!      real  pcgli(ngs), pcgld(ngs)
!      real  pcgmi(ngs), pcgmd(ngs)
!      real  pcghi(ngs), pcghd(ngs)
!
!      real  psctot(ngs) 
!      real  psccwi(ngs), psccwd(ngs)
!      real  psccii(ngs), psccid(ngs)
!      real  pscipi(ngs), pscipd(ngs)
!      real  psciri(ngs), pscird(ngs)
!      real  pscrwi(ngs), pscrwd(ngs)
!      real  pscswi(ngs), pscswd(ngs)
!      real  pschwi(ngs), pschwd(ngs)
!      real  pschli(ngs), pschld(ngs)
!      real  pscfwi(ngs), pscfwd(ngs)
!      real  pscgli(ngs), pscgld(ngs)
!      real  pscgmi(ngs), pscgmd(ngs)
!      real  pscghi(ngs), pscghd(ngs)

!      real  psccwmi(ngs), psccwmd(ngs)
!      real  psccimi(ngs), psccimd(ngs)
!      real  pscipmi(ngs), pscipmd(ngs)
!      real  pscirmi(ngs), pscirmd(ngs)
!      real  pscrwmi(ngs), pscrwmd(ngs)
!      real  pscswmi(ngs), pscswmd(ngs)
!      real  pschwmi(ngs), pschwmd(ngs)
!      real  pschlmi(ngs), pschlmd(ngs)
!      real  pscfwmi(ngs), pscfwmd(ngs)
!      real  pscglmi(ngs), pscglmd(ngs)
!      real  pscgmmi(ngs), pscgmmd(ngs)
!      real  pscghmi(ngs), pscghmd(ngs)
!
!  other arrays
!
!
      real wvdf(ngs),tka(ngs) !,akvisc(ngs),ci(ngs),cw(ngs),thdf(ngs)
!      real dqisdt(ngs)
      real advisc(ngs) !dqwsdt(ngs), ,schm(ngs),pndl(ngs)

      real qss0(ngs)

      real advisc0,advisc1,tka0

!      real qsacip(ngs)
      real pbz(ngs),pres(ngs),presp(ngs),pres0(ngs)
      real pk(ngs)
      real rho0(ngs),dnz(ngs),pi0(ngs),piz(ngs)
      real rhovt(ngs)
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
!      real thsave(ngs)
!      real pceds(ngs) ! ,ppceds(ngs),pmceds(ngs)
      real times
!      real qwfzi(ngs) ! ,qimlw(ngs)
!      real ptwfzi(ngs),ptimlw(ngs)
      real pvap(ngs),ptem(ngs) ! ,psub(ngs),pfrz(ngs)
!      real fload(ngs)
!      character*80 filnam
!      character*15 rrshcm
!      character*2  headr1
!      character*5  rstime
!      character*2  nmliter
!
!  iholef = 1 to do hole filling technique version 1
!  which uses all hydrometerors to do hole filling of all hydrometeors
!  iholef = 2 to do hole filling technique version 2
!  which uses an individual hydrometeror species to do hole 
!  filling of a species of a hydrometeor
!
!  iholen = interval that hole filling is done
!
      integer  iholef
      integer  iholen
      parameter (iholef = 9)
      parameter (iholen = 10000000)
      real  cqtotn !,cqtotn1
      real  cctotn
!      real  citotn
      real  crtotn
!      real  cstotn
      real  cvtotn
!      real  cftotn
!      real  cgltotn
!      real  cghtotn
!      real  chtotn
      real  cqtotp ! ,cqtotp1
      real  cctotp
!      real  citotp
!      real  ciptotp
      real  crtotp
!      real  cstotp
      real  cvtotp
!      real  cftotp
!      real  chltotp
!      real  chxtotp
!      real  cgltotp
!      real  cgmtotp
!      real  cghtotp
!      real  chtotp
      real  cqfac
      real  ccfac
!      real  cifac
!      real  cipfac
      real  crfac
!      real  csfac
      real  cvfac
!      real  cffac
!      real  cglfac
!      real  cghfac
!      real  chfac
!
      
      real,allocatable :: xvt(:,:,:,:) ! (nx,nz,2,lc:lhab) ! 1=mass-weighted, 2=number-weighted
      real,allocatable :: tmpn(:,:,:)
      real,allocatable :: tmpn2(:,:,:)
!     
!
!
!   Miscellaneous variables
!
!      real alpha(ngs)
      integer ireadqf ! ,lrho,lqsw,lqgl,lqgm ,lqgh ! ,ltim,ltem,lqcw,lqfw
      integer imkgam ! , lqrw
!      real vt
      data imkgam /0/
      integer igam        
      real vt
      real arg, gamma  ! gamma is a function  
!      real gaml02
!      real erbnd1, fdgt1, costhe1
      real qeps,qctol ! ,erbnd ,fdgt,costhe ,dslay,dxi,dyi,dzi,dxi2
      real cp608,cv,bta1,cnit,dragh,dnz00,rho00,pii
      real qccrit,gf4br,gf4ds,gf4p5, gf3ds, gf1ds,gr

!      real cnoi,cnoip,cnoir,cnor,cnos,cnogl,cnogm,cnogh,cnof,cnoh
!      real cnohl,
!      real rwdnmx,cwdnmx,cidnmx,xidnmx,swdnmx,gldnmx,gmdnmx
!      real cirdn0, cwdn0, rwdn0, swdn0, gldn0
!      real gmdn0, ghdn0, fwdn0, hwdn0, hldn0
      
      real xdn0(lc:lhab)

!      real ghdnmx,fwdnmx,hwdnmx,hldnmx,rwdnmn,cwdnmn,xidnmn,cidnmn
!      real swdnmn,gldnmn,gmdnmn,ghdnmn,fwdnmn
      integer l,itermax ,ltemq,inumgs !, idelq ! , ib
!      real hwdnmn,hldnmn,
      real c1f3,brz,arz,rw,temq ! ,cmn,cmi40,cmi50
!      real ri50,vti50,bsfw,cm50a,a,cm40b,cm50b
      real cnin20,cnin10,cnin1a
      real cnin2a,cnin2b,ssival,tqvcon
!      real cd(5)
      real cdx(lc:lhab)
      real cno(lc:lhab) ! , cnox
!      real cval,aval,eval,fval,gval ,qsign,ftelwc,qconkq
!      real qconm,qconn,cfce15,gf8,gf4i,gf1a,
      real gf1p5,qdiff,argrcnw,gf3p5
      real c4,bradp,bl2,bt2 ! ,dtrh,hrifac, hdia0,hdia1,civenta,civentb
!      real civentc,civentd,civente,civentf,civentg,cireyn,xcivent
!      real cipventa,cipventb,cipventc,cipventd,cipreyn,cirventa
!      real cirventb
      integer igmrwa,igmrwb ! ,igmswa, igmswb,igmfwa,igmfwb,igmhwa,igmhwb
      real rwventa ,rwventb ! ,swventa,swventb,fwventa,fwventb,fwventc
!      real hwventa,hwventb 
      integer iptotal ! counts number of times that eqtol is exceeded
      integer ilock1,ilock2,ilock3,ilock4,ilockc1,ilockc2
      real ptotalmx,ptotalmn
      real psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p
      integer nscmax,nscmin,iwetg
      
!      real    hwventc, hlventa, hlventb,  hlventc
!      real  glventa, glventb, glventc 
!      real   gmventa, gmventb,  gmventc, ghventa, ghventb, ghventc 

!      real  dzfacp,  dzfacm,  cmassin,  cwdiar ! , cwmasr
!      real  rimmas, rhobar
!      real   argtim !, argqcw, argqxw, argtem
!      real   frcswsw, frcswgl, frcswgm, frcswgh, frcswfw, frcswsw1
!      real   frcglgl, frcglgm, frcglgh,  frcglfw, frcglgl1
!      real   frcgmgl, frcgmgm, frcgmgh,  frcgmfw, frcgmgm1
!      real   frcghgl, frcghgm, frcghgh,  frcghfw,  frcghgh1
!      real   frcfwgl, frcfwgm, frcfwgh, frcfwfw,  frcfwfw1
!      real   frcswrsw, frcswrgl,  frcswrgm,  frcswrgh, frcswrfw
!      real   frcswrsw1
!      real   frcrswsw, frcrswgl, frcrswgm, frcrswgh, frcrswfw
!      real  frcrswsw1 
!      real  frcglrgl, frcglrgm, frcglrgh,  frcglrfw, frcglrgl1
!      real  frcrglgl  
!      real  frcrglgm,  frcrglgh, frcrglfw, frcrglgl1  
!      real  frcgmrgl, frcgmrgm, frcgmrgh, frcgmrfw,  frcgmrgm1
!      real  frcrgmgl, frcrgmgm,  frcrgmgh, frcrgmfw, frcrgmgm1
!      real  sum,  qweps,  gf2a, gf4a, dqldt, dqidt, dqdt
!      real frcghrgl, frcghrgm, frcghrgh, frcghrfw, frcghrgh1, frcrghgl
!      real frcrghgm, frcrghgh,  frcrghfw, frcrghgh1
      real    a1 ! ,a2,a3,a4,a5,a6
!      real   gamss
!      real cdw, cdi, denom1, denom2, delqci1, delqip1 ! , dtz1, dtz2
!      real cirtotn,  ciptotn, cgmtotn, chltotn,  cirtotp
!      real  cgmfac, chlfac,  cirfac
!      integer igmhla, igmhlb, igmgla, igmglb, igmgma,  igmgmb
!      integer igmgha, igmghb
!      integer idqis, item, itim0 ! ,  itim
!      integer  iqgl, iqgm, iqgh, iqrw, iqsw ! ,iqcw, iqfw
      integer  ia ! ,itertd

      save iread
      
!     volatile iptotal,ptotalmx,ptotalmn,nscmax,nscmin,iwetg,
!    &   psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      real    :: t0p3, t0p1


#ifdef MPI
      INCLUDE "mpif.h"
      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      integer :: downward_tag, upward_tag

      logical :: debug_mpi = .false.

      integer       :: nampi,nb
#endif
!
! ####################################################################
!
!  Start routine
!
! ####################################################################
!
!
!  slope intercepts
!
      IF ( ipconc > 3) ipconc = 3
      infdo  = 1

      IF ( infall .ge. 3 .or. lzr > 1 ) THEN
         infdo = 2
      ENDIF
      
      IF ( ngs .lt. nz ) THEN
       write(0,*) 'Error in WARMZIEG: Must have ngs .ge. nz!'
       STOP
      ENDIF
      
      allocate ( db1(nx,nz+1) )
      allocate ( db0(nx,nz+1) )
      allocate ( dtz1(nx,nz+1) )
      allocate ( dtz0(nx,nz+1) )
      allocate ( xvt(nx,nz+1,3,lc:lhab) )
      allocate ( tmpn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) )
      allocate ( tmpn2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) )
      allocate ( ztmp(-nor+ng1:nx+nor,-nor+ng1:nz+nor,lr:lhab) )

      
      ng = nor

!      write(0,*) my_rank, 'WZ: lnc,lnr,lc,lhab = ',lnc,lnr,lc,lhab,lr
      ln(lc) = lnc
      ln(lr) = lnr

      ipc(lc) = 2
      ipc(lr) = 3
      
      lz(:) = 0
      lz(lr) = lzr
      
      lvol(:) = 0

      cno(:) = 1.


      xnu(lc) = 0.0
      xmu(lc) = 1.
      
      xnu(lr) = -0.8
      xmu(lr) = 1.

      dnu(lc) = 3.*xnu(lc) + 2. ! alphac
      dmu(lc) = 3.*xmu(lc)

      dnu(lr) = 3.*xnu(lr) + 2. ! alphar
      dmu(lr) = 3.*xmu(lr)
      
!
!  density maximums and minimums
!
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
!
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0

      xdn0(lc) = 1000.0
      xdn0(lr) = 1000.0

!
!  Set terminal velocities...
!    also set drag coefficients
!
      cdx(lr) = 0.60
      cdx(lc) = 2.0

      
      dtp = ntmul*dtp1
      IF (my_rank == 0 ) write(iunit,*) 'IN TDMPHY'
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix)

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: TOP OF CODE - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('kz',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',10x,'temp',10x,'lss',10x,'lv',10x,'ss')")
        do kz = 1,nz
         write(0,"(i3,8(1x,es12.5))") kz,an(21,1,kz,lc), an(21,1,kz,lr), &
     &                                    an(21,1,kz,lnc),an(21,1,kz,lnr), &
     &                                    an(21,1,kz,lt),an(21,1,kz,lss), &
     &                                    an(21,1,kz,lv),ssat(21,1,kz)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: TOP OF CODE - - - - - - - - - "
        write(0,"('kz',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',10x,'temp',10x,'lss',10x,'lv',10x,'ss')")
        do kz = 1,nz
         write(0,"(i3,8(1x,es12.5))") kz,an(21,22,kz,lc), an(21,22,kz,lr), &
     &                                    an(21,22,kz,lnc),an(21,22,kz,lnr), &
     &                                    an(21,22,kz,lt),an(21,22,kz,lss), &
     &                                    an(21,22,kz,lv),ssat(21,22,kz)
        end do
#endif
!end mpidebug
      ENDIF

       call sett09s(nx,ny,nz,nor,nor,na,lt, &
     &             t77,pinit,p2,an,pn,pb, &
     &          t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00 )

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
       if (my_rank==0) then
       write(0,*) " - - - - - - - - - WARMZIEG: after temp set - - - - - - - - - "
       write(0,*) "my_rank=",my_rank
       write(0,"('k',7x,'t0',12x,'t00',11x,'t77',11x,'lt')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,t0(21,1,kz),t00(21,1,kz),t77(21,1,kz),an(21,1,kz,lt)
       end do
       end if
#else
       write(0,*) " - - - - - - - - - WARMZIEG: after temp set - - - - - - - - - "
       write(0,"('k',7x,'t0',12x,'t00',11x,'t77',11x,'lt')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,t0(21,22,kz),t00(21,22,kz),t77(21,22,kz),an(21,22,kz,lt)
       end do
#endif
!end mpidebug
       ENDIF


!
!  open input data file for the microphysics code
!
      if ( iread .eq. 0 ) then
!
      iread = 1
!

      write(iunit,*) 'INPUT:  ndebug (-1=no debug)'
      write(iunit,*) 'ndebug= ',ndebug
      
      write(iunit,*) 'ipconc = ',ipconc



      ffrz = fconv
      
      write(iunit,*) 'ircnw, qminrncw=', ircnw, qminrncw

      write(iunit,*) 'cwdiap=', cwdiap

      write(iunit,*) 'cwdisp= ', cwdisp

      write(iunit,*) 'iauttim= ',iauttim

      write(iunit,*) 'auttim= ',auttim

      write(iunit,*) 'qcwmntim= ', qcwmntim

      cwccn = ccn
      alfarim = Min (0.0,alfarim)
      write(iunit,*) 'alfarim = ', alfarim

         pi = 4.0*atan(1.0)
!         xvcmx = (4./3.)*pi*xcradmx**3

         write(iunit,*) 'iccwflg, irenuc = ', iccwflg, irenuc
         
         write(iunit,*) 'xcradmx, xvcmx = ', xcradmx, xvcmx
         
         write(iunit,*) 'ciintmx = ',ciintmx
         
      cck = 0.6
      cckm = cck-1.
!      ccn = cwccn
      ccnefac =  (1.63/(cck * beta(3./2., cck/2.)))**(cck/(cck + 2.0))
      cnexp   = (3./2.)*cck/(cck+2.0)
!      ccne = 0.9893*1.e6*(1.e-6*Abs(cwccn))**0.77
      write(iunit,*) 'cwccn, cck, ccne = ',cwccn,cck,ccne

      IF ( .false. .and. nstep .eq. 1 .and. na .ge. lccn ) THEN
       write(iunit,*) 'Initializing CCN field, nx,ny,nz,na,lccn=', &
     &          nx,ny,nz,na,lccn
       kbound = 1

#ifdef MPI
       kzb = 1
       kze = ktile+1
       if(kzend .eq. nzend) kze = kzend-kzbeg
       
       DO kz=kzb,kze
#else
       DO kz=1,nz-1
#endif
        IF ( cwccn .gt. 0.0 ) THEN

        cwccn = ccn
        ab(kz,lccn) = cwccn*db(kz)/db(kbound)
        write(iunit,*) 'CCN = ',kz,ab(kz,lccn)
#ifdef MPI
        jyb = -nor+1
        jye = jtile+nor

        ixb = -nor+1
        ixe = itile+nor
       
        DO jy = jyb,jye ; DO ix = ixb,ixe
#else
        DO jy = -nor+ng1,ny+nor
         DO ix = -nor+ng1,nx+nor
#endif
           an(ix,jy,kz,lccn) = ab(kz,lccn) ! cwccn*dn(ix,jy,kz)/dn(ix,jy,kbound)
!           ab(kz,lccn) = an(ix,jy,kz,lccn)
#ifdef MPI
          IF ( kzbeg-1+kz .eq. 1 ) THEN
#else
          IF ( kz .eq. 1 ) THEN
#endif
           an(ix,jy,0,lccn) = ab(kz,lccn) ! cwccn*dn(ix,jy,kz)/dn(ix,jy,kbound)
          ENDIF
#ifdef MPI
          IF ( kzbeg-1+kz .eq. nzend-1 ) THEN
#else
          IF ( kz .eq. nz - 1 ) THEN
#endif
           an(ix,jy,nz,lccn) = ab(kz,lccn) ! cwccn*dn(ix,jy,kz)/dn(ix,jy,kbound)
          ENDIF
         ENDDO
        ENDDO
        
        ELSE
        
        ab(kz,lccn) = Abs(cwccn)
        write(iunit,*) 'CCN = ',kz,ab(kz,lccn)
#ifdef MPI
        jyb = -nor+1
        jye = jtile+nor

        ixb = -nor+1
        ixe = itile+nor

        DO jy = jyb,jye ; DO ix = ixb,ixe
#else
        DO jy = -nor+ng1,ny+nor
         DO ix = -nor+ng1,nx+nor
#endif
           an(ix,jy,kz,lccn) = Abs(cwccn)
!           ab(kz,lccn) = an(ix,jy,kz,lccn)
#ifdef MPI
          IF ( kzbeg-1+kz .eq. 1 ) THEN
#else
          IF ( kz .eq. 1 ) THEN
#endif
           an(ix,jy,0,lccn) = Abs(cwccn) 
          ENDIF
#ifdef MPI
          IF ( kzbeg-1+kz .eq. nzend-1 ) THEN
#else
          IF ( kz .eq. nz - 1 ) THEN
#endif
           an(ix,jy,nz,lccn) = Abs(cwccn)
          ENDIF
         ENDDO
        ENDDO

        ENDIF
       ENDDO
      ENDIF

      IF ( cwccn .lt. 0.0 ) THEN
      cwccn = Abs(cwccn)
      ccwmx = cwccn
      ELSE
      ccwmx = cwccn*1.4
      ENDIF

      write(iunit,'(a,i2,2x,i2,1x,1pe12.5,3i3)')  &
     &  'itype1, itype2,cimas0, icfn, ihrn, ibfc = ', &
     &    itype1, itype2, cimas0, icfn, ihrn, ibfc


      IF ( irfall .lt. 0 ) irfall = infall
      IF ( lzr > 0 ) irfall = 0
      IF ( my_rank == 0 ) write(iunit,*) 'itfall, irfall,infall = ',itfall, irfall,infall

!      dtrim = rimtim
!      argtim = ((rimtim-0.5/xftim)*xftim + yftim)
!      itim0 = min(max(ifix(argtim),1), mtim)
!      write(iunit,*) 'itim0 = ',itim0


      
      ireadqf = 0

!

      f5 = 237.3 * 17.27 * 2.5e6 / cp 

      do l = 1,nqsat
      temq = 163.15 + (l-1)*fqsat
      tabqvs(l) = exp(caw*(temq-273.16)/(temq-cbw))
      dtabqvs(l) = ((-caw*(-273.16 + temq))/(temq - cbw)**2 + &
     &                 caw/(temq - cbw))*tabqvs(l)
      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
      end do
!
      end if    ! ( iread .eq. 0 ) 
!
!

!  Gamma function
!
!
!     compute fgoi(arg*100)
!
      if ( imkgam .eq. 0 ) then
      imkgam = 1
      call makegmoi()
!      do igam = 1,1000
!      arg = 0.01*igam
!      gmoi(igam) = gamma(arg)
!c      write(97,*) igam,gmoi(igam)
!      end do
      end if
!
!
!  electricity constants
!
!  mixing ratio epsilon
!
      qeps  = 1.0e-20
!
!  cloud mixing ratio tolerance for screening layer parameterization
!
      qctol = 1.0e-5
!
!  rebound efficiency (erbnd)
!
!
!  screening layer depth
!
!      dslay = 500.0
!
!  grid intervals
!
!      dxi   = (1.0)/max(dx, qeps)
!      dyi   = (1.0)/max(dy, qeps)
!      dzi   = (1.0)/max(dz, qeps)  ! not used
!      dxi2  = (0.5)/max(dx, qeps)
!      dyi2  = (0.5)/max(dy, qeps)
!      dzi2  = (0.5)/max(dz, qeps)  ! not used
!
!      times = (nstep-1)*(dtp1*dtfac)
      times = (nstep)*(dtp1*dtfac)
!      headr1 = 'pa'
!
!  constants
!
      cp608 = 0.608
      cv = 717.0
      ar = 841.99666  
      br = 0.8
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
      rho00 = 1.225
!      cs = 4.83607122
!      ds = 0.25
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pi = 4.0*atan(1.0)
      pii = 1./pi
      pid4 = pi/4.0 
!      qscrit = 6.0e-04
      gf1 = gamma(1.0)
      gf1p5 = 0.8862269255
      gf2 = gamma(2.0)
      gf3 = gamma(3.0)
      gf3p5 = 3.32335097
      gf4 = gamma(4.0)
      gf5 = gamma(5.0)
      gf6 = gamma(6.0)
      gf7 = gamma(7.0)
      gf4br = gamma(4.0+br)
      gf4ds = gamma(4.0+ds)
      gf4p5 = gamma(4.0+0.5)
      gf3ds = gamma(3.0+ds)
      gf1ds = gamma(1.0+ds)
      gr = 9.8
      gf73rds = gamma(7./3.)
      gf83rds = gamma(8./3.)
!
!  constants
!
      c1f3 = 1.0/3.0
!
!  general constants for microphysics
!
      brz = 100.0
      arz = 0.66

      call setqxminz(qxmin)
      
      
      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx


      tdtol = 1.0e-05
      thnuc = 233.15
      rw = 461.5
      advisc0 = 1.832e-05
      advisc1 = 1.718e-05
      tka0 = 2.43e-02
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
!
!
      if ( imake .eq. 0 ) then
      imake = 1

      ventr   = Gamma(rnu + 4./3.)/(rnu + 1.)**(1./3.)/Gamma(rnu + 1.)
      ventc   = Gamma(cnu + 4./3.)/(cnu + 1.)**(1./3.)/Gamma(cnu + 1.)

      xnu(lc) = 0.0
      xmu(lc) = 1.
      
      xnu(lr) = -0.8
      xmu(lr) = 1.

!
!  Set up look up tables for supersaturation w.r.t. liq and ice
!
!VD$L SKIP
!      do l = 1,nqsat
!      temq = 163.15 + (l-1)*fqsat
!      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
!      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
!      end do
            
!
      itermax = 0
      DO mgs = 1,ngs
        il0(mgs) = 1
      ENDDO
      
      c1sw = Gamma(snu + 4./3.)*(snu + 1.0)**(-1./3.)/gamma(snu + 1.0) 

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      DO kz = kzb,kze
#else
      DO kz = 1,nz-kstag
#endif
      piz(kz) = pinit(kz) ! (pb(kz)/poo)**cap ! pinit(kz)+p2(igs(mgs),jy,kgs(mgs)) 
      dnz(kz) = db(kz)
      pbz(kz) = pb(kz)
      temp(kz) = piz(kz)*ab(kz,lt)

      tempc(kz) = temp(kz) - tfr
      tempbz(kz) = temp(kz) 
      elv(kz) = 2500300. - 2369.276*tempc(kz) 
      elf(kz) = 335717. - 2369.276*tempc(kz) 
      els(kz) = elv(kz) + elf(kz)
      elvs(kz) = elv(kz)*elv(kz)
      elss(kz) = els(kz)*els(kz)
      gamw(kz) = elv(kz)*cpi/piz(kz)
      gams(kz) = els(kz)*cpi/piz(kz)
      cqv1(kz) = 4098.0258*piz(kz)*gamw(kz)
      cqv2(kz) = 5807.6953*piz(kz)*gams(kz)
      cc3(kz) = cpi*elf(kz)/piz(kz)
      
      ENDDO

!
!  cw constants in mks units
!
!      cwmasn = 4.25e-15  ! radius of 1.0e-6
      cwmasn = 5.23e-13  ! radius of 5.0e-6
      cwmasn5 =  5.23e-13
      cwradn = 5.0e-6
      cwmasx = 5.25e-10  ! radius of 50.0e-6
      mwfac = 6.0**(1./3.)
      IF ( ipconc .ge. 2 ) THEN
        cwmasn = xvmn(lc)*1000.
        cwradn = 1.0e-6
        cwmasx = xvmx(lc)*1000.
      ENDIF
        rwmasn = xvmn(lr)*1000.
        rwmasx = xvmx(lr)*1000.

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
#else
      do kz = 1,nz-kstag 
#endif
!      cwdn(kz) = 1000.0
      cwmcnd(kz) = 5.0e-12  
      cwncmn(kz) = 1.0e+01
      cwncmx(kz) = 1.0e+10
      cwnccn(kz) = cwccn*dnz(kz)/dnz(1)  
      cwc1 = 6.0/(pi*1000.)
 
      end do

      cwc0 = pii ! 6.0*pii
!
!  ci constants in mks units
!
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter
#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
#else
      do kz = 1,nz-kstag 
#endif
!      cidn(kz) = 900.0
      cimcnd(kz) = 1.0e-12 
      cincmn(kz) = 1.0e+01
      cincmx(kz) = 1.0e+10
!      cinccn(kz) = 1.0e+09
      end do

! 
!  constants for paramerization
!
      if ( ncuse .eq. 0 ) then

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kstag
       
      do kz = kzb,kze
#else
      do kz = 1,nz-kstag
#endif
      wvdf(kz) = (2.11e-05)*((temp(kz)/tfr)**1.94)* &
     &  (101325.0/(pbz(kz)))
      advisc(kz) = advisc0*(416.16/(temp(kz)+120.0))* &
     &  (temp(kz)/296.0)**(1.5)
      tka(kz) = tka0*advisc(kz)/advisc1
      end do
      end if
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank.eq.0) write(0,*) 'WARMZIEG: dbg = 1'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 1'
#endif
!
!
!
!  end if for constants (imake)
!
      end if
!
!  zero the precip flux arrays (2d)
!
!      if ( nstep .eq. nstart ) then
!      do jy = 1,ny-jstag
!      do ix = 1,nx-istag
!      ciptot2d(ix,jy) = 0.0
!      cirtot2d(ix,jy) = 0.0
!      cwtot2d(ix,jy) = 0.0
!      rwtot2d(ix,jy) = 0.0
!      swtot2d(ix,jy) = 0.0
!      fwtot2d(ix,jy) = 0.0
!      gltot2d(ix,jy) = 0.0
!      gmtot2d(ix,jy) = 0.0
!      ghtot2d(ix,jy) = 0.0
!      hwtot2d(ix,jy) = 0.0
!      hltot2d(ix,jy) = 0.0
!      end do
!      end do
!      if ( ihabdo .lt. 1 ) then
!      do jy = 1,ny-jstag
!      do ix = 1,nx-istag
!      citot2d(ix,jy) = 0.0
!      end do
!      end do
!      end if     
!      end if     

!
!  set 3-d temperature and saturation mixing ratio variables 
!
      ierr = 0

      IF ( ndebug .ge. 0 ) write(iunit,*) 'WARMZIEG: ltemq set'

#ifdef MPI
      kzb = 0
      kze = ktile+1
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ltemq)
!! c$omp PARALLEL DO  IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix,ltemq)
      do kz = 1,nz-1
      do jy = 1,ny-jstag
      do ix = 1,nx-istag
#endif

!  temperature

!      t0(ix,jy,kz) = an(ix,jy,kz,lt)
!     >  *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
      
      
!      t0(ix,jy,kz) = 
!     :     theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)))
!
!  look-up index
! 
      ltemq = Int( (t0(ix,jy,kz)-163.15)/fqsat+1.5 )
      IF ( ltemq .ge. nqsat .or. ltemq .le. 0 ) THEN
        ierr = ierr + 1
        write(iunit,*) 'WARNING 1: ltemq .ge. nqsat! ierr = ',ierr
        write(iunit,*) 'ix,jy,kz,t0,an,ltemq,t77 = ',ix,jy,kz,t0(ix,jy,kz), &
     &    an(ix,jy,kz,lt), ltemq,t77(ix,jy,kz)
        write(0,*) 'WARNING 1: ltemq .ge. nqsat! ierr = ',ierr
        write(0,*) 'ix,jy,kz,t0,an,ltemq = ',ix,jy,kz,t0(ix,jy,kz), &
     &    an(ix,jy,kz,lt), ltemq
        ltemq = Max(1, Min (nqsat,ltemq))
       
      ENDIF
      IF ( ierr .gt. 4 ) THEN
       write(0,*) 'too many warnings! STOP'
       write(iunit,*) 'too many warnings! STOP'
       STOP
      ENDIF
!
! saturation mixing ratio for liquid
!
      a1 = t00(ix,jy,kz)

      t8(ix,jy,kz) = a1*tabqvs(ltemq)

!
! saturation mixing ratio for ice
!
      t9(ix,jy,kz) = a1*tabqis(ltemq)
!
      end do
      end do
      end do

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
       if (my_rank==0) then
       write(0,*) " - - - - - - - - - WARMZIEG: ltemq set - - - - - - - - - "
       write(0,*) "my_rank=",my_rank
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,1,kz,lc),an(21,1,kz,lr),an(21,1,kz,lnc),an(21,1,kz,lnr)
       end do
       end if
#else
       write(0,*) " - - - - - - - - - WARMZIEG: ltemq set - - - - - - - - - "
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
#endif
      ENDIF
!end mpidebug
      
!      IF ( ierr .gt. 4 ) THEN
!        call flush(iunit)
!        STOP
!      ENDIF
!
!
!  Stuff for Meyers et al. (1992) and Ferrier (1994) ice nucleation
!
!  constants
!
      cnin20 = 1.0e3
      cnin10 = 5.0e1
      cnin1a = 4.5
      cnin2a = 12.96
      cnin2b = 0.639
!
!  calculate rate of nucleation
!

      IF ( ndebug .ge. 0 ) write(iunit,*) 'WARMZIEG: t7 set'

#ifdef MPI
      kzb = 0
      kze = ktile+1
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ssival)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix,ssival)
      do kz = 1,nz-1
      do jy = 1,ny-jstag
      do ix = 1,nx-istag
#endif

! filter supersaturation field (if ifilt = 1)
!
      ssival = Min(t8(ix,jy,kz),max(an(ix,jy,kz,lv),0.0))/t9(ix,jy,kz)
!
!      t7(ix,jy,kz) = 0.0
!
      if ( ssival .gt. 1.0 ) then
!
      if ( t0(ix,jy,kz).le.268.15 ) then
!     if ( t0(ix,jy,kz).le.273.15 ) then
        
       dp1 = cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,jy,kz) = Min(dp1, 1.0d30)
      end if
      
!
!   have turned off nucleation by Meyer at warm temperatures
!
!     if ( t0(ix,jy,kz).lt.273.55.and.t0(ix,jy,kz).gt.268.15 ) then
!     qvapor = max(an(ix,jy,kz,lv),0.0) 
!     ssifac(mgs) = 0.0
!     if ( (qvapor-t9(ix,jy,kz)) .gt. 1.0e-5 ) then
!     if ( (t8(ix,jy,kz)-t9(ix,jy,kz)) .gt. 1.0e-5 ) then
!     ssifac(mgs) = (qvapor-t9(ix,jy,kz))
!    >             /(t8(ix,jy,kz)-t9(ix,jy,kz))
!     ssifac(mgs) = ssifac(mgs)**cnin1a   
!     end if
!     end if
!     t7(ix,jy,kz) = cnin10*ssifac(mgs)
!    >  *exp(-(t0(ix,jy,kz)-tfr)*bta1)
!     end if
!
      end if
!
      end do
      end do
      end do

#ifdef MPI
      kzb = -nor+1
      kze = ktile+nor

      jyb = -nor+1
      jye = jtile+nor

      ixb = -nor+1
      ixe = itile+nor

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix)
      DO kz=-nor+ng1,nz+nor
      DO jy=-nor+ng1,ny+nor
      DO ix=-nor+ng1,nx+nor
#endif
        t8(ix,jy,kz) = 0.0
        t9(ix,jy,kz) = 0.0
      END DO
      END DO
      END DO
!
!  set save counter (number of saves):  nsvcnt
!
      nsvcnt = 0
      iend = 0
      iptotal = 0
      nscmax = 0
      nscmin = 0
      ptotalmx = 0.0
      ptotalmn = 0.0
      psctotmx = 0.0
      psctotmn = 0.0     
      psctot1 = 0.0
      psctot2 = 0.0 
      sctot1n = 0.0
      sctot1p = 0.0
      ilock1 = 0
      ilock2 = 0
      ilock3 = 0
      iwetg = 0
      ilock4 = 0
      ilockc1 = 1
      ilockc2 = 2
     

!      timetd1 = etime(tarray)
!      timetd1 = tarray(1)

! 
!$     ndebug = -1
! cmic$  cncall
!***********************************************************
!  start jy loop
!***********************************************************
!
      IF ( ndebug .ge. 0 ) write(iunit,*) 'WARMZIEG: start jy loop'
!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
       if (my_rank==0) then
       write(0,*) " - - - - - - - - - WARMZIEG: about to start jy loop - - - - - - - - - "
       write(0,*) "my_rank=",my_rank
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,1,kz,lc),an(21,1,kz,lr),an(21,1,kz,lnc),an(21,1,kz,lnr)
       end do
       end if
#else
       write(0,*) " - - - - - - - - - WARMZIEG: about to start jy loop - - - - - - - - - "
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
#endif
      ENDIF
!end mpidebug

      db0(:,:) = 1.0
      DO kz = 1,nz+1
       dtz0(:,kz) = z1d(kz,3)
      ENDDO

#ifdef MPI
      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      do 9999 jy = jyb,jye
#else
      do 9999 jy = 1,ny-jstag
#endif
!  VERY IMPORTANT:  SET jgs = jy
!
      jgs = jy
      
!
!  Do fallout stuff now if itfall=0 (boxfall) or 2 (crowfall)
!
      IF ( itfall .eq. 0 .or. itfall .eq. 1 ) THEN

!
!  zero the precip flux arrays (2d)
!

      xvt(:,:,:,:) = 0.0

      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze 
      DO ix = ixb,ixe
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo


      DO kz = kzb,kze
      DO ix = ixb,ixe
       dtz1(ix,kz) = z1d(kz,3)/db1(ix,kz)
      ENDDO
      ENDDO

! set up vt arrays:
!      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
!     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)

      IF ( ndebug .ge. 0 ) write(iunit,*) 'WARMZIEG: call ziegfall, jy =',jy ! ,jyb,jye
      IF ( ndebug .ge. 0 ) write(0,*) my_rank,'WARMZIEG: call ziegfall, jy =',jy ! ,jyb,jye

      call ziegfall(nx,ny,nz,nor,nor,na,dtp,dz,jgs,ngs, &
     &  xvt, &
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, &
     &  rwmasn,rwmasx,cwradn, &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn, &
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

      DO il = lc,lhab

      vtmax = 0.0

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
      
      do kz = kzb,kze
      do ix = ixb,ixe
      
! fix gt for SWM
!      DO il = lr,lhab
      vtmax = Max(vtmax,xvt(ix,kz,1,il)*gz(kz))
      vtmax = Max(vtmax,xvt(ix,kz,2,il)*gz(kz))
      vtmax = Max(vtmax,xvt(ix,kz,3,il)*gz(kz))
!      ENDDO

      ENDDO
      ENDDO


      
      IF ( dtp*vtmax .lt. 0.7 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*vtmax/0.7) + 1)
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp/Real(ndfall)
!        write(0,*) 'subdivide fallout on my_rank = ',my_rank
!        write(0,*) 'for il,jyc = ',il,jy,dtp*vtmax
      ELSE
        dtptmp = dtp
      ENDIF

      DO n = 1,ndfall

      IF ( n .ge. 2 ) THEN
!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0
      
      kzb=1
      kze=ktile+1
      if(kzend.eq.nzend) kze=kzend-kzbeg

      ixb=1
      ixe=itile
      if(ixend.eq.nxend) ixe=ixend-ixbeg
      
      do kz = kzb,kze
      DO ix = ixb,ixe
       db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

      call ziegfall(nx,ny,nz,nor,nor,na,dtptmp,dz,jgs,ngs, &
     &  xvt, &
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, &
     &  rwmasn,rwmasx,cwradn, &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn, &
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

      ENDIF ! (n .ge. 2)

        IF ( il >= lr .and. ( infall .eq. 3 .or. infall .eq. 4 ) ) THEN
!         DO il = lr,lhab
           IF ( (il .eq. lr .and. irfall .eq. infall .and. lzr < 1)  ) THEN
            call calczgr(nx,ny,nz,nor,na,an,ixe,kze,ztmp, &
     &         db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), &
               xvmn(il), xvmx(il), lvol(il), rho_qh , -1)
           ENDIF
!         ENDDO
        ENDIF
!     :  rwmasn,rwmasx)

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3b"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 1b'
#endif

! mixing ratio

!      DO il = lc,lhab
!      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), 
!     &             an,db1,imapz,mzdist,il,1,xfall)
      call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), &
     &             an,db1,imapz,mzdist,il,1,xfall,dtz1)
!      ENDDO



#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3c"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3c'
#endif


! reflectivity

        IF ( lz(il) .gt. 1 ) THEN
!         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), 
!     &              an,db0,imapz,mzdist,lz(il),0,xfall)
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), &
     &              an,db0,imapz,mzdist,lz(il),0,xfall,dtz0)
        ENDIF

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3d"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3d'
#endif

! mean amount

      
      IF ( ipconc .gt. 0 ) THEN !{
!      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN

      IF ( ( infall .ge. 2  ) .and. lz(il) .lt. 1) THEN !{
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!


          kzb=1
          kze=ktile+1
          if (kzend.eq.nzend) kze = kzend-kzbeg

          ixb=1
          ixe=itile
          if (ixend.eq.nxend) ixe = ixend-ixbeg

!        IF ( infall .eq. 3 .and. ( il .eq. lr .or. il .eq. lh .or. il .eq. lhl ) ) THEN
        IF ( ( infall .eq. 3 .or. infall .eq. 4 ) .and. ( &
     &      ( il .eq. lr .and. irfall .eq. infall) ) ) THEN
!          call calczgr(nx,ny,nz,nor,na,an,ixe,kze,
!     :       z,db1,jgs,ipconc, microp, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il) )

          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn2(ix,jy,kz) = ztmp(ix,kz,il)
            ENDDO
          ENDDO
          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
        
        ELSE
          
          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO

        ENDIF

      ENDIF !}


#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3f"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3f'
#endif

       in = 2
       IF ( infall .eq. 1 ) in = 1

!         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), 
!     &        an,db0,imapz,mzdist,ln(il),0,xfall)
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall,dtz0)


         IF ( lz(il) .lt. 1 ) THEN ! if not 3-moment, run one of the correction schemes
         IF ( (infall .ge. 2 .or. infall .eq. 3) .and. ( il .eq. lr  )) THEN
!     :        .or. il .eq. lhl )) THEN
           
           xfall0(:,jgs) = 0.0

           IF ( ( infall .eq. 3 .or. infall .eq. 4 ) .and.  &
     &        ( il .ge. lh .or. (il .eq. lr .and. irfall .eq. infall) ) ) THEN
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), & 
!     &         tmpn2,db0,imapz,mzdist,1,0,xfall0)           
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &         tmpn,db0,imapz,mzdist,1,0,xfall0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), &
     &         tmpn2,db0,imapz,mzdist,1,0,xfall0,dtz0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), &
     &         tmpn,db0,imapz,mzdist,1,0,xfall0,dtz0)
           ELSE
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), 
!     &         tmpn,db0,imapz,mzdist,1,0,xfall0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), &
     &         tmpn,db0,imapz,mzdist,1,0,xfall0,dtz0)
           ENDIF

           IF ( ( infall .eq. 3 .or. infall .eq. 4 ) .and. ( (il .eq. lr .and. irfall .eq. infall) &
     &           ) ) THEN
! "Method I" - dbz correction
             kze = ktile+1
             if (kzend .eq. nzend) kze = kzend-kzbeg

             ixe = itile
             if (ixend .eq. nxend) ixe = ixend-ixbeg

             call calcnfromz(nx,ny,nz,nor,na,an,tmpn2,ixe,kze, &
     &       ztmp,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il),tmpn,  &
     &       lvol(il), rho_qh, infall )

!             IF ( il == lr ) THEN
!             DO kz = kzb,kze
!              DO ix = ixb,ixe
!                IF ( xvt(ix,kz,3,il) .ne. xvt(ix,kz,1,il) ) THEN
!                  write(6,*) 'DR: xvt1,3 = ',xvt(ix,kz,1,il),xvt(ix,kz,3,il)
!                ENDIF
!              ENDDO
!             ENDDO
!             ENDIF

!             DO kz = kzb,kze
!              DO ix = ixb,ixe
!               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), 0.5* ( an(ix,jgs,kz,ln(il)) + tmpn(ix,jy,kz) ))
!              
!              ENDDO
!             ENDDO           
           ELSEIF ( infall .eq. 5 .and. il .ge. lh .or. ( il == lr .and. irfall == 5 ) ) THEN

             kze = ktile+1
             if (kzend .eq. nzend) kze = kzend-kzbeg

             ixe = itile
             if (ixend .eq. nxend) ixe = ixend-ixbeg


             DO kz = kzb,kze
              DO ix = ixb,ixe
               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), 0.5* ( an(ix,jgs,kz,ln(il)) + tmpn(ix,jy,kz) ))
              
              ENDDO
             ENDDO           
           ELSEIF ( .not. (il .eq. lr .and. irfall .eq. 0) ) THEN
! "Method II" M-wgt N-fallout correction
             kzb = 1
             kze = ktile+1
             if (kzend .eq. nzend) kze = kzend-kzbeg

             ixb = 1
             ixe = itile
             if (ixend .eq. nxend) ixe = ixend-ixbeg

             DO kz = kzb,kze
              DO ix = ixb,ixe

               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
              
              ENDDO
             ENDDO
           ENDIF 
           ENDIF ! lz(il) .lt. 1
           
!          ENDDO

         ENDIF
        ENDIF
!      ENDDO


      ENDIF !}

      ENDDO ! n=1,ndfall
      ENDDO ! il

      
      ELSEIF ( itfall .eq. 2 .or. itfall .eq. 4 ) THEN
      
!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0
      
#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze ; DO ix = ixb,ixe
#else
      do kz = 1,nz-1
      DO ix = 1,nx-1
#endif
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

      call ziegfall(nx,ny,nz,nor,nor,na,dtp,dz,jgs,ngs, &
     &  xvt, &
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, &
     &  rwmasn,rwmasx,cwradn, &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn, &
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

      vtmax = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
      do ix = 1,nx-1
#endif
      
! fix gt for SWM
      DO il = lr,lhab
      vtmax = Max(vtmax,xvt(ix,kz,1,il)*gz(kz))
      ENDDO

      ENDDO
      ENDDO


      
      IF ( dtp*vtmax .lt. 0.9 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*vtmax) + 1)
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp
        dtp = dtp/Real(ndfall)
      ENDIF

      DO n = 1,ndfall

      IF ( n .ge. 2 ) THEN
!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0
      
#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze ; DO ix = ixb,ixe
#else
      do kz = 1,nz-1
      DO ix = 1,nx-1
#endif
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

      call ziegfall(nx,ny,nz,nor,nor,na,dtp,dz,jgs,ngs, &
     &  xvt, &
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, &
     &  rwmasn,rwmasx,cwradn, &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn, &
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

      ENDIF ! (n .ge. 2)

! mixing ratio

      DO il = lc,lhab
      IF ( itfall .eq. 2 ) THEN
      call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), &
     &  an,db1,imapz,mzdist,il,1,xfall)
      ELSEIF ( itfall .eq. 4 ) THEN
      call wenofall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), &
     &  an,db1,imapz,mzdist,il,1,xfall)
      ENDIF
      ENDDO
      


!      infall = 2
!      DO il = lc,lhab
!        IF ( ipconc .ge. ipc(il) ) THEN
!         call crowfall(nx,ny,nz,nor,na,gt,ngt,dtp,dz,jgs,xvt(1,1,infall,il),
!     :        an,db0,imapz,mzdist,ln(il),0,xfall)
!        ENDIF
!      ENDDO
      IF ( ipconc .gt. 0 ) THEN
      DO il = lc,lhab

      IF ( ipconc .ge. ipc(il) ) THEN
      IF ( infall .eq. 2 ) THEN
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!
#ifdef MPI
          kzb = 1
          kze = ktile+1
          if(kzend .eq. nzend) kze = kzend-kzbeg

          ixb = 1
          ixe = itile
          if(ixend .eq. nxend) ixe = ixend-ixbeg

          do kz = kzb,kze ; DO ix = ixb,ixe
#else
          DO kz = 1,nz
            DO ix = 1,nx
#endif
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
      ENDIF

        IF ( itfall .eq. 2 ) THEN
         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
        ELSEIF ( itfall .eq. 4 ) THEN
         call wenofall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
        ENDIF

!         call boxfall(nx,ny,nz,nor,na,gt,ngt,dtp,dz,jgs,xvt(1,1,infall,il),
!     :        an,db0,imapz,mzdist,ln(il),0,xfall)

         IF ( infall .eq. 2 .and. ( il .eq. lr  )) THEN

        IF ( itfall .eq. 2 ) THEN
         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
        ELSEIF ( itfall .eq. 4 ) THEN
         call wenofall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
        ENDIF
!           call boxfall(nx,ny,nz,nor,1,gt,ngt,dtp,dz,jgs,xvt(1,1,1,il),
!     :       tmpn,db0,imapz,mzdist,1,0,xfall0)

#ifdef MPI
           kzb = 1
           kze = ktile
           if(kzend .eq. nzend) kze = kzend-kzbeg

           ixb = 1
           ixe = itile
           if(ixend .eq. nxend) ixe = ixend-ixbeg

           do kz = kzb,kze ; DO ix = ixb,ixe
#else
           DO kz=1,nz-1
           DO ix=1,nx-1
#endif
             an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
             
           ENDDO
           ENDDO
           
!          ENDDO

         ENDIF
        ENDIF
      ENDDO
      
      ENDIF


      ENDDO ! n=1,ndfall
      
      IF ( ndfall .gt. 1 ) THEN
        dtp = dtptmp
      ENDIF

      ENDIF ! itfall.eq.2


!
!  zero the precip rate arrays (2d)
!
!      do ix = 1,nx-istag
!      ciprt2d(ix,jy) = 0.0
!      cirrt2d(ix,jy) = 0.0
!      cwrt2d(ix,jy) = 0.0
!      rwrt2d(ix,jy) = 0.0
!      swrt2d(ix,jy) = 0.0
!      fwrt2d(ix,jy) = 0.0
!      glrt2d(ix,jy) = 0.0
!      gmrt2d(ix,jy) = 0.0
!      ghrt2d(ix,jy) = 0.0
!      hwrt2d(ix,jy) = 0.0
!      hlrt2d(ix,jy) = 0.0
!      prt2d(ix,jy) = 0.0
!      prti2d(ix,jy) = 0.0
!      prtl2d(ix,jy) = 0.0
!      end do
!      do ix = 1,nx-istag
!      cirt2d(ix,jy) = 0.0
!      end do
!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
       if (jy == 22) then
#ifdef MPI
       if (my_rank==0) then
       write(0,*) " - - - - - - - - - WARMZIEG: before gather - - - - - - - - - "
       write(0,*) "my_rank=",my_rank
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
       end if
#else
       write(0,*) " - - - - - - - - - WARMZIEG: before gather - - - - - - - - - "
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
#endif
       end if
       ENDIF
!end mpidebug

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
      DO ix = ixb,ixe
#else
      do kz = 1,nz-1
      DO ix = 1,nx-1
#endif
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo
      
      vtmax = 0.0
      
!      IF ( fallonly ) cycle ! GOTO 3999
      IF ( fallonly ) GOTO 3999
       
!      do kz = 1,nz-kstag
!      do ix = 1,nx-istag
!      civt(ix,kz) = 0.0
!      end do
!      end do
!
!..Gather microphysics  
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank.eq.0) write(0,*) 'WARMZIEG: GATHER STAGE 1'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: GATHER STAGE 1'
#endif

      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
#ifdef MPI
      ixe = itile
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg+1-istag

      nxz = ixe*(kzend-kzbeg-kstag)
      numgs = nxz/ngs + 1
#endif

      do 1000 inumgs = 1,numgs 

      ngscnt = 0

#ifdef MPI
      kzb = nzmpb
      kze = ktile
!      if(kzbeg .le. nzmpb .and. kzend .gt. nzmpb) kzb = nzmpb-kzbeg+1
      if(kzend .eq. nzend) kze = kzend-kzbeg-kstag

      ixb = nxmpb
      ixe = itile
!      if(ixbeg .le. nxmpb .and. ixend .gt. nxmpb) ixb = nxmpb-ixbeg+1
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do 1001 kz = kzb,kze 
      do 1002 ix = nxmpb,ixe
#else
      do 1001 kz = nzmpb,nz-kstag-1 
      do 1002 ix = nxmpb,nx-istag
#endif

       pqs(kz) = t00(ix,jy,kz)

       theta(kz) = an(ix,jy,kz,lt) 
       temg(kz) = t0(ix,jy,kz)
       temcg(kz) = temg(kz) - tfr
       tqvcon = temg(kz)-cbw
       ltemq = (temg(kz)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       qvs(kz) = pqs(kz)*tabqvs(ltemq)
       qis(kz) = pqs(kz)*tabqis(ltemq)

       qss(kz) = qvs(kz)
       if ( temg(kz) .lt. tfr ) then
        qss(kz) = qis(kz)
       end if
!
       if ( an(ix,jy,kz,lv)  .gt. qss(kz) .or. &
     &      an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.  &
     &      an(ix,jy,kz,lr)  .gt. qxmin(lr) ) then
        ngscnt = ngscnt + 1
        igs(ngscnt) = ix
        kgs(ngscnt) = kz
        if ( ngscnt .eq. ngs ) goto 1100
       end if
 1002 continue
      nxmpb = 1
 1001 continue
!      if ( jy .eq. (ny-jstag) ) iend = 1
 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank.eq.0) write(0,*) 'WARMZIEG: DBG=5'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: DBG=5'
#endif
!      write(0,*) 'allocating qc
      allocate (  qx(ngscnt,lv:lhab) )
      allocate (  qxw(ngscnt,ls:lhab) )
      allocate (  cx(ngscnt,lc:lhab) )
      allocate (  cxmxd(ngscnt,lc:lhab) )
      allocate (  qxmxd(ngscnt,lv:lhab) )
      allocate (  xv(ngscnt,lc:lhab) )
      allocate (  vtxbar(ngscnt,lc:lhab,3) )
      allocate (  xmas(ngscnt,lc:lhab) )
      allocate (  xdn(ngscnt,lc:lhab) )
      allocate (  xdia(ngscnt,lc:lhab,3) )
      allocate (  alpha(ngscnt,lc:lhab) )
!      IF ( ipconc .ge. 6 ) THEN
      allocate (  zx(ngscnt,lr:lhab) )
      allocate (  zxmxd(ngscnt,lr:lhab) )
      allocate (  g1x(ngscnt,lr:lhab) )
!      ENDIF
      
      xv(:,:) = 0.0
      xmas(:,:) = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:) = 0.0
      zx(:,:) = 0.0
      alpha(:,:) = 0.0
      ventrx(:) = ventr
!
!  define temporaries for state variables to be used in calculations
!
      do 1010 mgs = 1,ngscnt
        temp0(mgs) = tempbz(kgs(mgs))
        theta0(mgs) = ab(kgs(mgs),lt)
        thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt) &
     &              - ab(kgs(mgs),lt)
        theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
        pres0(mgs) = pb(kgs(mgs))
        qv0(mgs) = ab(kgs(mgs),lv)
        qwvp(mgs) = an(igs(mgs),jy,kgs(mgs),lv)  &
     &            - ab(kgs(mgs),lv)

        presp(mgs) = pn(igs(mgs),jy,kgs(mgs))
        pres(mgs) = presp(mgs) + pres0(mgs)
        rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
        rhoinv(mgs) = 1.0/rho0(mgs)
        rhovt(mgs) = Sqrt(rho00/rho0(mgs))
        pi0(mgs) = piz(kgs(mgs))
        temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
        pk(mgs)   = t77(igs(mgs),jy,kgs(mgs))
        temcg(mgs) = temg(mgs) - tfr
        qss0(mgs) = (380.0)/(pres(mgs))
        pqs(mgs) = (380.0)/(pres(mgs))
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
!
        il5(mgs) = 0
        if ( temg(mgs) .lt. tfr ) then 
         il5(mgs) = 1
        end if
 1010 continue
!
!
!  set temporaries for microphysics variables
!
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

!
!  6th moments
!

      IF ( lzr > 1 ) THEN
       zx(:,:) = 0.0
#ifdef MPI
      if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = load Z'
#endif
         DO mgs = 1,ngscnt
          zx(mgs,lr) = Max( an(igs(mgs),jy,kgs(mgs),lzr), 0.0 ) 
         ENDDO
             
      ENDIF

!
!  set concentrations
!
!      ssmax = 0.0
      
      cx(:,:) = 0.0
      
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        IF ( lccn .gt. 1 ) THEN
         ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        ELSE
         ccnc(mgs) = 0.0
        ENDIF
       end do
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
          cx(mgs,lr) = 0.0
        ELSEIF ( cx(mgs,lr) .le. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) THEN
          qx(mgs,lr) = 0.0
        ELSE
          cx(mgs,lr) = Max( 1.e-9, cx(mgs,lr) )
        ENDIF
       end do

!
!  set factors
!
      IF ( ndebug .ge. 0 ) write(iunit,*) 'WARMZIEG: set factors'

      do mgs = 1,ngscnt
!
      ssi(mgs) = qx(mgs,lv)/qis(mgs)
      ssw(mgs) = qx(mgs,lv)/qvs(mgs)
!
      tsqr(mgs) = temg(mgs)**2
!
      temgx(mgs) = min(temg(mgs),313.15)
      temgx(mgs) = max(temgx(mgs),233.15)
      felv(mgs) = 2500837.367  &
     &  * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
!
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),223.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      felf(mgs) = 333690.6098  &
     &  + (2030.61425)*temcgx(mgs) &
     &  - (10.46708312)*temcgx(mgs)**2
!
      fels(mgs) = felv(mgs) + felf(mgs)
!
      felvs(mgs) = felv(mgs)*felv(mgs)
      felss(mgs) = fels(mgs)*fels(mgs)
!
      fgamw(mgs) = felv(mgs)*cpi/pi0(mgs)
      fgams(mgs) = fels(mgs)*cpi/pi0(mgs)
!
      fcqv1(mgs) = 4098.0258*pi0(mgs)*fgamw(mgs)
      fcqv2(mgs) = 5807.6953*pi0(mgs)*fgams(mgs)
      fcc3(mgs) = cpi*felf(mgs)/pi0(mgs)
!
      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)* &
     &  (101325.0/(pres(mgs)))
!
      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &  (temg(mgs)/296.0)**(1.5)
!
      fakvisc(mgs) = fadvisc(mgs)/rho0(mgs)
!
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fci(mgs) = (2.118636 + 0.007371*(temcgx(mgs)))*(1.0e+03)
!
      if ( temg(mgs) .lt. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4203.1548  + (1.30572e-2)*((temcgx(mgs)-35.)**2) &
     &                 + (1.60056e-5)*((temcgx(mgs)-35.)**4)
      end if
      if ( temg(mgs) .ge. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),308.15)
      temcgx(mgs) = max(temcgx(mgs),273.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4243.1688  + (3.47104e-1)*(temcgx(mgs)**2)
      end if
!
      ftka(mgs) = tka0*fadvisc(mgs)/advisc1
      fthdf(mgs) = ftka(mgs)/(cp*rho0(mgs))
!
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))
      fpndl(mgs) = (fakvisc(mgs)/fthdf(mgs))
!
      fai(mgs) = (fels(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbi(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qis(mgs)))
      fav(mgs) = (felv(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbv(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qvs(mgs)))
!
      end do
!
!
!   ice habit fractions
!
!
!
!  Set density
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank.eq.0) write(0,*) 'WARMZIEG: Set density'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Set density'
#endif

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
      end do

      g1shr = 1.0
      g1mlr = 1.0
      
      ! set base g1x in case rain is not 3-moment
       IF ( ipconc >= 6 ) THEN
         il = lr
         DO mgs = 1,ngscnt
           g1x(mgs,il) = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
         ENDDO
       ENDIF

!  Find shape parameter rain

      IF ( lzr > 1 ) THEN ! { RAIN SHAPE PARAM
      
       IF ( imurain == 3 ) THEN
         alphashr = 0.0
         g1shr = 36.*(alphashr+2.0)/((alphashr+1.0)*pi**2)
         alphamlr = -2.0/3.0
         g1mlr = 36.*(alphamlr+2.0)/((alphamlr+1.0)*pi**2)
       ELSEIF ( imurain == 1 ) THEN
         alphashr = 2.0
         g1shr = (6.0 + alphashr)*(5.0 + alphashr)*(4.0 + alphashr)/ &
     &            ((3.0 + alphashr)*(2.0 + alphashr)*(1.0 + alphashr))
         alphamlr = 0.0
         g1mlr = (6.0 + alphamlr)*(5.0 + alphamlr)*(4.0 + alphamlr)/ &
     &            ((3.0 + alphamlr)*(2.0 + alphamlr)*(1.0 + alphamlr))
       ENDIF
      
      CALL cld_cpu('Z-MOMENT-1r')  
          il = lr
          DO mgs = 1,ngscnt
          

         IF ( iresetmoments == 1 .or. iresetmoments == il  ) THEN
         IF ( zx(mgs,lr) <= zxmin ) THEN
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           qx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
         ELSEIF ( cx(mgs,lr) <= cxmin ) THEN
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         ENDIF
         ENDIF
         
         IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*xdn(mgs,lr)**2)
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
          ENDIF
        
          IF ( zx(mgs,lr) > 0.0 ) THEN
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            vr = xv(mgs,lr)
!            z = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z,xv
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

! check for artificial breakup (rain larger than allowed max size)
        IF (  xv(mgs,il) .gt. xvmx(il) .or. (ioldlimiter == 3 .and. xv(mgs,il) .gt. xvmx(il)/8.)) THEN
          tmp = cx(mgs,il)
          IF ( ioldlimiter == 3 ) THEN ! MY-style active breakup
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = 1.0 ! x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.e4*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.e4*x1**2), x3)
          ELSE ! simple cutoff 
            xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          ENDIF
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))

          IF ( tmp < cx(mgs,il) ) THEN ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           IF ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) THEN

            IF ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alphar .and. alp <= rnumin ) THEN
             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = zx(mgs,il)
            ENDIF
           ENDIF
           
         ! set g1x to use as G factor later. If alpha is in the range ( rnumin < alpha < rnumax ), then 
         ! this will be the same as computing G from alpha.  If alpha = rnumax, however, it probably means that
         ! the moments aren't matched correctly, so we compute G from the moments instead so that the dZ/dt rates
         ! stay consistent with dN/dt and dq/dt.
           IF ( alp >= rnumax - 0.01 ) THEN
             g1x(mgs,il) = zx(mgs,il)/(cx(mgs,il)*xv(mgs,lr)**2)
           ELSE
             g1x(mgs,il) = g1
           ENDIF
           
           tmp = alpha(mgs,lr) + 4./3.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma(alpha(mgs,lr) + 1.)
           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))
           
           ENDIF
          ENDIF
          
          ENDIF
          
          ENDDO
        CALL cld_cpu('Z-MOMENT-1r')  
        ENDIF ! }
!
!  set some values for ice nucleation
!
      do mgs = 1,ngscnt
!      uvel(mgs) = 0.0 
!      vvel(mgs) = 0.0 
      uvel(mgs) = (0.5)*(u(igs(mgs),jgs,kgs(mgs)) &
     &                  +u(igs(mgs)+1,jgs,kgs(mgs)))
      vvel(mgs) = (0.5)*(v(igs(mgs),jgs,kgs(mgs)) &
     &                  +v(igs(mgs),jgs+1,kgs(mgs)))
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)+1) &
     &                  +w(igs(mgs),jgs,kgs(mgs)))

#ifdef MPI
      if (kzbeg-1+kgs(mgs) .eq. 1) then
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                    +w(igs(mgs),jgs,1))
        kgsm(mgs) = 1
        kgsp(mgs) = kgs(mgs)+1
      else
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                    +w(igs(mgs),jgs,kgs(mgs)-1))
        kgsm(mgs) = kgs(mgs)-1
        if (kzend-1+kgs(mgs) .ge. nzend-1) then
          kgsp(mgs) = nz-1
        else
          kgsp(mgs) = kgs(mgs)+1
        endif
      endif
#else
      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))
      kgsm(mgs) = max(kgs(mgs)-1,1)
      kgsp(mgs) = min(kgs(mgs)+1,nz-1)
#endif

      cninm(mgs) = t7(igs(mgs),jgs,kgsm(mgs))
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      cninp(mgs) = t7(igs(mgs),jgs,kgsp(mgs))
      end do

!
!  Set a couple of cloud variables...
!

!      SUBROUTINE setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
!     :                 xmas,xdn,xvmn,xvmx,xv,cdx,
!     :                 ipconc,ndebug)

      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno, &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx, &
     &                 ipconc,ndebug,ngs,nz,igs,kgs,cwnccn,fadvisc, &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx, &
     &                 itype1,itype2,temcg,0,alpha,alpha,axx,bxx,0)

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
       if (jy == 22) then
#ifdef MPI
       if (my_rank==0) then
       write(0,*) " - - - - - - - - - WARMZIEG: after setvtz - - - - - - - - - "
       write(0,*) "my_rank=",my_rank
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
       end if
#else
       write(0,*) " - - - - - - - - - WARMZIEG: after setvtz - - - - - - - - - "
       write(0,"('k',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
       do kz = 1,nz
         write(0,"(i3,4(1x,es12.5))") kz,an(21,22,kz,lc),an(21,22,kz,lr),an(21,22,kz,lnc),an(21,22,kz,lnr)
       end do
#endif
       end if
       ENDIF
!end mpidebug

!
!  Set number concentrations (need xdia from setvt)
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank.eq.0) write(0,*) 'WARMZIEG: Set concentration'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Set concentration'
#endif
      
      if ( ipconc .lt. 5 ) then
      do mgs = 1,ngscnt


      IF ( ipconc .lt. 3 ) THEN
!      cx(mgs,lr) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lh) )  then
!      cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
      end if
      ENDIF


      end do
      end if
      
      IF ( ipconc .ge. 2 ) THEN
      DO mgs = 1,ngscnt
        rb(mgs) = 0.5*xdia(mgs,lc,1)*((1./(1.+cnu)))**(1./6.)
        xl2p(mgs) = Max(0.0d0, 2.7e-2*xdn(mgs,lc)*cx(mgs,lc)*xv(mgs,lc)*   &
     &           ((0.5e20*rb(mgs)**3*xdia(mgs,lc,1))-0.4) )
        IF ( rb(mgs) .gt. 3.51e-6 ) THEN
!          rh(mgs) = Max( 0.5d0*xdia(mgs,lc,1), 6.3d-4/(1.d6*(rb(mgs) - 3.5d-6)) )
          rh(mgs) = Max( 41.d-6, 6.3d-4/(1.d6*(rb(mgs) - 3.5d-6)) )
        ELSE
          rh(mgs) = 41.d-6
        ENDIF
        IF ( xl2p(mgs) .gt. 0.0 ) THEN
          nh(mgs) = 4.2d9*xl2p(mgs)
        ELSE
          nh(mgs) = 1.e30
        ENDIF
      ENDDO
      ENDIF
      
!
!
!
!              
!
!  maximum depletion tendency by any one source
!
!
      if( ndebug .ge. 0 ) THEN
#ifdef MPI
       if (my_rank==0) write(0,*) 'WARMZIEG: Set depletion max/min1'
#else
       write(0,*) 'WARMZIEG: Set depletion max/min1'
#endif
       call flush(iunit)
      endif
      do mgs = 1,ngscnt

      frac = 0.1
      qcmxd(mgs)  = frac*qx(mgs,lc)/dtp
      qrmxd(mgs)  = frac*qx(mgs,lr)/dtp
      end do
!
      if( ndebug .ge. 0 ) THEN
#ifdef MPI
       if (my_rank==0) write(0,*) 'WARMZIEG: Set depletion max/min2'
#else
       write(0,*) 'WARMZIEG: Set depletion max/min2'
#endif
       call flush(iunit)
      endif

      do mgs = 1,ngscnt
!  
      if ( qx(mgs,lc) .le. qxmin(lc) ) then
      ccmxd(mgs)  = 0.20*cx(mgs,lc)/dtp
      else
      IF ( ipconc .ge. 2 ) THEN
        ccmxd(mgs)  = frac*cx(mgs,lc)/dtp
      ELSE
        ccmxd(mgs)  = frac*qx(mgs,lc)/(xmas(mgs,lc)*rho0(mgs)*dtp)
      ENDIF
      end if
!

      crmxd(mgs)  = 0.10*cx(mgs,lr)/dtp

      ccmxd(mgs)  = frac*cx(mgs,lc)/dtp
      crmxd(mgs)  = frac*cx(mgs,lr)/dtp

      qxmxd(mgs,lv) = Max(0.0, frac*(qx(mgs,lv)-qss(mgs))/dtp )

      DO il = lc,lhab
       qxmxd(mgs,il) = frac*qx(mgs,il)/dtp
       cxmxd(mgs,il) = frac*cx(mgs,il)/dtp
      ENDDO

      end do
!
      if( ndebug .ge. 0 ) THEN
#ifdef MPI
       if (my_rank==0) write(0,*) 'WARMZIEG: Setup inductive charging...'
#else
       write(0,*) 'WARMZIEG: Setup inductive charging...'
#endif
       call flush(iunit)
      endif
 
!
!
!
!
!
!  microphysics source terms (1/s) for mixing ratios 
!
!
!
!  Collection efficiencies:
!
! 
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
       if (my_rank==0) write(0,*) 'WARMZIEG: Set collection efficiencies'
#else
       write(0,*) 'WARMZIEG: Set collection efficiencies'
#endif
      ENDIF
      
      do mgs = 1,ngscnt
!
!
!
      erw(mgs) = 0.0
!
      err(mgs) = 0.0
!
      icwr(mgs) = 1
      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN
       cwrad = 0.5*xdia(mgs,lc,1)
      DO il = 1,8
         IF ( cwrad .ge. 1.e-6*cwr(il,1) ) icwr(mgs) = il
      ENDDO
      ENDIF

      irwr(mgs) = 1
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
         rwrad = 0.5*xdia(mgs,lr,3)  ! changed to mean volume diameter (10/6/06)
      DO il = 1,6
         IF ( rwrad .ge. 1.e-6*grad(il,1) ) irwr(mgs) = il
      ENDDO
      ENDIF
!
!
!
!  Rain: Collection (cxc) efficiencies
!
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lc).gt.qxmin(lc) ) then

       IF ( lnr .gt. 1 ) THEN
       erw(mgs) = 1.0
       
       ELSE

!      cwrad = 0.5*xdia(mgs,lc,1)
!      erw(mgs) =
!     >  min((aradcw + cwrad*(bradcw + cwrad*
!     <  (cradcw + cwrad*(dradcw)))), 1.0)
!       IF ( xdia(mgs,lc,1) .lt. 2.4e-06 .or. xdia(mgs,lr,1) .le. 50.0e-6 ) THEN
!          erw(mgs)=0.0
!       ENDIF
!       erw(mgs) = ew(icwr(mgs),igwr(mgs))
! interpolate along droplet radius
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = irwr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lr,1)
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)

!       write(iunit,*) 'slop1: ',slope1,slope2,ew(ic,ir),cwr(ic,2)

       x1 = ew(ic,  ir) + slope1*Max(0.0, (cwrad - cwr(ic,1)) )
       x2 = ew(icp1,ir) + slope2*Max(0.0, (cwrad - cwr(ic,1)) )
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       erw(mgs) = Max(0.0, x1 + slope1*Max(0.0, (rwrad - grad(ir,1)) ))
       
!       write(iunit,*) 'erw: ',erw(mgs),1.e6*cwrad,1.e6*rwrad,ic,ir,x1,x2
!       write(iunit,*)
       
       erw(mgs) = Max(0.0, erw(mgs) )
       IF ( rwrad .lt. 50.e-6 ) THEN
         erw(mgs) = 0.0
       ELSEIF (  rwrad .lt. 100.e-6 ) THEN  ! linear change from zero at 50 to erw at 100 microns
         erw(mgs) = erw(mgs)*(rwrad - 50.e-6)/50.e-6
       ENDIF
       
       ENDIF
      end if
      IF ( cx(mgs,lc) .le. 0.0 ) erw(mgs) = 0.0
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lr).gt.qxmin(lr) ) then
      err(mgs)=1.0
      end if

      
      ENDDO
!
!
!
!
!  Collection growth equations....
!
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: Collection: rain collects xxxxx'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Collection: rain collects xxxxx'
#endif
!
      do mgs = 1,ngscnt
      qracw(mgs) =  0.0
      IF ( qx(mgs,lr) .gt. qxmin(lr) .and. erw(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .lt. 3 ) THEN
       IF ( erw(mgs) .gt. 0.0 .and. qx(mgs,lr) .gt. 1.e-7 ) THEN
       vt = (ar*(xdia(mgs,lc,1)**br))*rhovt(mgs)
       qracw(mgs) =    &
     &   (0.25)*pi*erw(mgs)*qx(mgs,lc)*cx(mgs,lr) &
!     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *Max(0.0, vtxbar(mgs,lr,1)-vt)   &
     &  *(  gf3*xdia(mgs,lr,2)    &
     &    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )  
       ENDIF
      ELSE

       rwrad = 0.5*xdia(mgs,lr,1)
        IF ( rwrad .gt. rh(mgs) ) THEN ! .or. cx(mgs,lr) .gt. nh(mgs) ) THEN
         IF ( rwrad .gt. rwradmn ) THEN
!      DM1CCC=A2*XNC*XNR*XVC*(((CNU+2.)/(CNU+1.))*XVC+XVR)       ! (A12)
           qracw(mgs) = erw(mgs)*aa2*cx(mgs,lr)*cx(mgs,lc)*xmas(mgs,lc)*   &
     &        ((cnu + 2.)*xv(mgs,lc)/(cnu + 1.) + xv(mgs,lr))/rho0(mgs) !*rhoinv(mgs)
         ELSE

!      DM1CCC=A1*XNC*XNR*(((CNU+3.)*(CNU+2.)/(CNU+1.)**2)*XVC**3+ ! (A14)
!     1 ((RNU+2.)/(RNU+1.))*XVC*XVR**2)

           qracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*xdn(mgs,lc)*   &
     &        ((cnu + 3.)*(cnu + 2.)*xv(mgs,lc)**3/(cnu + 1.)**2 +    &
     &         (alpha(mgs,lr) + 2.)*xv(mgs,lc)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))/rho0(mgs) !*rhoinv(mgs)

!           xvc = xv(mgs,lc)*(1.e6)
!           xvr = xv(mgs,lr)*1.e6
           
!           qracw(mgs) = 1.e-18*(aa1*xvc*cx(mgs,lr)*cx(mgs,lc)*xdn(mgs,lc)*
!     :        ((cnu + 3.)*(cnu + 2.)*xvc**2/(cnu + 1.)**2 + 
!     :         (alpha(mgs,lr) + 2.)*xvr**2/(alpha(mgs,lr) + 1.))/rho0(mgs)) !*rhoinv(mgs)
         ENDIF
        ENDIF
       ENDIF
!       qracw(mgs) = Min(qracw(mgs), qx(mgs,lc))
       qracw(mgs) = Min(qracw(mgs), qcmxd(mgs))
       ENDIF
      end do
!
!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: RAIN ACC CLOUD - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'qracw')")
        do mgs = 1,ngscnt
          write(0,"(i3,1(1x,es12.5))") mgs,qracw(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: RAIN ACC CLOUD - - - - - - - - - "
        write(0,"('mgs',5x,'qracw')")
        do mgs = 1,ngscnt
          write(0,"(i3,1(1x,es12.5))") mgs,qracw(mgs)
        end do
#endif
      ENDIF
!end mpidebug
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: conc 18'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: conc 18'
#endif

      if ( ipconc .ge. 2 ) then
      do mgs = 1,ngscnt
       cracw(mgs) = 0.0
       cracr(mgs) = 0.0
       ec0(mgs) = 1.e9
      IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr)    &
     &      .and. qracw(mgs) .gt. 0.0 ) THEN

       IF ( ipconc .lt. 3 ) THEN
        IF ( erw(mgs) .gt. 0.0 ) THEN
        cracw(mgs) =   &
     &   ((0.25)*pi)*erw(mgs)*cx(mgs,lc)*cx(mgs,lr)   &
     &  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf1*xdia(mgs,lc,2)   &
     &    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lr,1)   &
     &    + gf3*xdia(mgs,lr,2) )
        ENDIF
       ELSE ! IF ( ipconc .ge. 3 .and. 
        IF ( 0.5*xdia(mgs,lr,1) .gt. rh(mgs) ) THEN !  .or. cx(mgs,lr) .gt. nh(mgs) ) THEN
!        IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( 0.5*xdia(mgs,lr,1) .gt. rwradmn ) THEN ! 50.e-6 ) THEN
!      DM0CCC=A2*XNC*XNR*(XVC+XVR)                               ! (A11)
          cracw(mgs) = aa2*cx(mgs,lr)*cx(mgs,lc)*(xv(mgs,lc) + xv(mgs,lr))
        ELSE
!      DM0CCC=A1*XNC*XNR*(((CNU+2.)/(CNU+1.))*XVC**2+            ! (A13)
!     1 ((RNU+2.)/(RNU+1.))*XVR**2)
          cracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*   &
     &        ((cnu + 2.)*xv(mgs,lc)**2/(cnu + 1.) +    &
     &         (alpha(mgs,lr) + 2.)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))
        ENDIF
        ENDIF
       ENDIF
      ENDIF

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: RAIN ACC CLOUD CONC - - - - - - - - - "
        if (mgs ==1) write(0,*) "my_rank=",my_rank
        if (mgs ==1) write(0,"('mgs',5x,'cracw')")
        write(0,"(i3,1(1x,es12.5))") mgs,cracw(mgs)
        end if
#else
        if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: RAIN ACC CLOUD CONC - - - - - - - - - "
        if (mgs ==1) write(0,"('mgs',5x,'cracw')")
        write(0,"(i3,1(1x,es12.5))") mgs,cracw(mgs)
#endif
       ENDIF
!end mpidebug

! Rain self collection (cracr) and break-up (factor of ec0)
!       
        ec0(mgs) = 2.e9
        IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        rwrad = 0.5*xdia(mgs,lr,1)
        IF ( xdia(mgs,lr,1) .gt. 2.0e-3 ) THEN
          ec0(mgs) = 0.0
          cracr(mgs) = 0.0
        ELSE
          IF ( xdia(mgs,lr,1) .lt. 6.1e-4 ) THEN
            ec0(mgs) = 1.0
          ELSE
            ec0(mgs) = Exp(-50.0*(50.0*(xdia(mgs,lr,1) - 6.0e-4)))
          ENDIF
          

          IF ( rwrad .ge. 50.e-6 ) THEN
            cracr(mgs) = ec0(mgs)*aa2*cx(mgs,lr)**2*xv(mgs,lr)
          ELSE
            cracr(mgs) = ec0(mgs)*aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 2.)/(alpha(mgs,lr) + 1.)
          ENDIF
!          cracr(mgs) = Min(cracr(mgs),crmxd(mgs))
        ENDIF !! ( rwrad .gt. 1.0e-3 )
        ENDIF !! ( qx(mgs,lr) .gt. qxmin(lr) )
        
!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: RAIN SELFCOLLECTION - - - - - - - - - "
        if (mgs ==1) write(0,*) "my_rank=",my_rank
        if (mgs ==1) write(0,"('mgs',5x,'cracr')")
          write(0,"(i3,1(1x,es12.5))") mgs,cracr(mgs)
        end if
#else
        if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: RAIN SELFCOLLECTION - - - - - - - - - "
        if (mgs ==1) write(0,"('mgs',5x,'cracr')")
          write(0,"(i3,1(1x,es12.5))") mgs,cracr(mgs)
#endif
       ENDIF
!end mpidebug
        
!      cracw(mgs) = min(cracw(mgs),ccmxd(mgs))
      
      IF ( rcond == 0 ) THEN
!        cracr(mgs) = 0.0
      ENDIF
      
      end do
      end if !! ( ipconc .ge. 2  )
        
!

!
! Ziegler (1985) autoconversion
!
!
      IF ( ipconc .ge. 2 ) THEN
      
      DO mgs = 1,ngscnt
        zrcnw(mgs) = 0.0
        qrcnw(mgs) = 0.0
        crcnw(mgs) = 0.0
        cautn(mgs) = 0.0
      ENDDO

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: BEFORE ZIEG AUTOCONVERSION - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: BEFORE ZIEG AUTOCONVERSION - - - - - - - - - "
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
#endif
      ENDIF
!end mpidebug
      
      DO mgs = 1,ngscnt
!      qracw(mgs) = 0.0
!      cracw(mgs) = 0.0
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1000. .and. temg(mgs) .gt. tfrh+4. ) THEN
         volb = xv(mgs,lc)*(1./(1.+CNU))**(1./2.)
         cautn(mgs) = Min(ccmxd(mgs),   &
     &      ((CNU+2.)/(CNU+1.))*aa1*cx(mgs,lc)**2*xv(mgs,lc)**2)
         cautn(mgs) = Max( 0.0d0, cautn(mgs) )
         IF ( rb(mgs) .le. 7.51d-6 ) THEN
           t2s = 1.d30
!           cautn(mgs) = 0.0
         ELSE
!         XL2P=2.7E-2*XNC*XVC*((1.E12*RB**3*RC)-0.4)
         
!        T2S=3.72E-3/(((1.E4*RB)-7.5)*XNC*XVC) 
!           t2s = 3.72E-3/(((1.e6*rb)-7.5)*cx(mgs,lc)*xv(mgs,lc))
!           t2s = 3.72/(((1.e6*rb(mgs))-7.5)*rho0(mgs)*qx(mgs,lc))
           t2s = 3.72/(1.e6*(rb(mgs)-7.500d-6)*rho0(mgs)*qx(mgs,lc))

           qrcnw(mgs) = Max( 0.0d0, xl2p(mgs)/(t2s*rho0(mgs)) )
!          IF ( iturbenhance > 0 ) THEN
!            tmp = Min(1.5, (3.94 + (5.46-3.94)/(0.15)*tke(mgs))/3.94 )
!            cautn(mgs) = cautn(mgs)*tmp
!            qrcnw(mgs) = qrcnw(mgs)*tmp
!          ENDIF
           crcnw(mgs) = Max( 0.0d0, Min(3.5e9*xl2p(mgs)/t2s,0.5*cautn(mgs)) )
           
           
           IF ( crcnw(mgs) < 1.e-30 ) qrcnw(mgs) = 0.0

           IF ( lzr > 1 .and. qrcnw(mgs) > 0.0 ) THEN
!             vr = rho0(mgs)*qrcnw(mgs)/(1000.*crcnw(mgs))
!             zrcnw(mgs) = 36.*(xnu(lr)+2.0)*crcnw(mgs)*vr**2/((xnu(lr)+1.0)*pi**2)
             vr = rho0(mgs)*qrcnw(mgs)/(1000.)
             zrcnw(mgs) = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
           ENDIF

!           IF (  crcnw(mgs) .gt. cautn(mgs) .and. crcnw(mgs) .gt. 1.0 )
!     :          THEN
!             print*, 'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s,qx(mgs,lr)
!             print*, '            ',qx(mgs,lc),cx(mgs,lc),0.5e6*xdia(mgs,lc,1)
!             print*, '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :         1.e6*(( 3/(4.*pi))*rho0(mgs)*qrcnw(mgs)/
!     :       (crcnw(mgs)*xdn(mgs,lr)))**(1./3.),rh(mgs)*1.e6,rwrad(mgs)
!           ELSEIF ( crcnw(mgs) .gt. 1.0 .and. cautn(mgs) .gt. 0.) THEN
!             print*, 'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s
!             print*, '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :  1.e6*(( 3*pi/4.)*rho0(mgs)*qrcnw(mgs)/
!     :   (crcnw(mgs)*xdn(mgs,lr)))**(1./3.)
!           ENDIF
!           crcnw(mgs) = Min(cautn(mgs),3.5e9*xl2p(mgs)/t2s)
           
!           IF ( qrcnw(mgs) .gt. 0.3e-2 ) THEN
!            print*, 'QRCNW'
!            print*, qrcnw(mgs),crcnw(mgs),cautn(mgs)
!            print*, xl2p,t2s,rho0(mgs),xv(mgs,lc),cx(mgs,lc),qx(mgs,lc)
!            print*, rb,0.5*xdia(mgs,lc,1),mgs,igs(mgs),kgs(mgs)
!           ENDIF
!           qrcnw(mgs) = Min(qrcnw(mgs),qcmxd(mgs))
         ENDIF
         
         
       ENDIF
      ENDDO

!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: AUTOCONVERSION - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: AUTOCONVERSION - - - - - - - - - "
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
#endif
      ENDIF
!end mpidebug

      ELSE

!
!  Berry 1968 auto conversion for rain (Orville & Kopp 1977)
!
!
      if ( ircnw .eq. 4 ) then
      do mgs = 1,ngscnt
!      sconvmix(lcw,mgs) = 0.0
      qrcnw(mgs) =  0.0
      qdiff = max((qx(mgs,lc)-qminrncw),0.0)
      if ( qdiff .gt. 0.0 .and. xdia(mgs,lc,1) .gt. 20.0e-6 ) then
      argrcnw =   &
     &  ((1.2e-4)+(1.596e-12)*(cx(mgs,lc)*1.0e-6)   &
     &  /(cwdisp*qdiff*1.0e-3*rho0(mgs)))
      qrcnw(mgs) = (rho0(mgs)*1e-3)*(qdiff**2)/argrcnw
!      sconvmix(lcw,mgs) = max(sconvmix(lcw,mgs),0.0)
      qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      end do
      
      ENDIF
!
!
!
!  Berry 1968 auto conversion for rain (Ferrier 1994)
!
!
      if ( ircnw .eq. 5 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      qrcnw(mgs) =  0.0
      qccrit = (pi/6.)*(cx(mgs,lc)*cwdiap**3)*xdn(mgs,lc)/rho0(mgs)
      qdiff = max((qx(mgs,lc)-qccrit),0.)
      if ( qdiff .gt. 0.0 .and. cx(mgs,lc) .gt. 1.0 ) then
      argrcnw = &
!     >  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)/(cwdisp*rho0(mgs)*qdiff))   &
     &  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)*1.0e-3/(cwdisp*rho0(mgs)*qdiff))
      qrcnw(mgs) = &
!     >  timflg(mgs)*rho0(mgs)*(qdiff**2)/argrcnw   &
     &  1.0e-3*rho0(mgs)*(qdiff**2)/argrcnw
      qrcnw(mgs) = Min(qxmxd(mgs,lc), (max(qrcnw(mgs),0.0)) )
      
!      write(iunit,*) 'qrcnw,cx =',qrcnw(mgs),cx(mgs,lc),mgs,1.e3*qx(mgs,lc),cno(lr)
      end if
      end do
      end if
      
!
!
!  kessler auto conversion for rain.
!
      if ( ircnw .eq. 2 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      qrcnw(mgs) = (0.001)*max((qx(mgs,lc)-qminrncw),0.0)
      end do
      end if
!
!  c4 = pi/6
!  c1 = 0.12-0.32 for colorado storms...typically 0.3-0.4
!  berry reinhart type conversion (proctor 1988)
!
      if ( ircnw .eq. 1 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      c1 = 0.2
      c4 = pi/(6.0)
      bradp =    &
     & (1.e+06) * ((c1/(0.38))**(1./3.)) * (xdia(mgs,lc,1)*(0.5))
      bl2 =   &
     & (0.027) * ((100.0)*(bradp**3)*(xdia(mgs,lc,1)*(0.5)) - (0.4))
      bt2 = (bradp -7.5) / (3.72)
      qrcnw(mgs) = 0.0
      if ( bl2 .gt. 0.0 .and. bt2 .gt. 0.0 ) then
      qrcnw(mgs) = bl2 * bt2 * rho0(mgs)   &
     &  * qx(mgs,lc) * qx(mgs,lc)
      end if
      end do
      end if
      
      
      
      ENDIF  !  ( ipconc .ge. 2 )

!
!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: after all AUTOCONVERSION - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: after all AUTOCONVERSION - - - - - - - - - "
        write(0,"('mgs',5x,'cracr',9x,'qrcnw',9x,'crcnw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,4(1x,es12.5))") mgs,cracr(mgs),qrcnw(mgs),crcnw(mgs),cautn(mgs)
        end do
#endif
      ENDIF
!end mpidebug

!
!  Ventilation coeficients
!
      do mgs = 1,ngscnt
      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))
      end do
!
!
!
!
      igmrwa = 100.0*2.0
      igmrwb = 100.*((5.0+br)/2.0)
      rwventa = (0.78)*gmoi(igmrwa)  ! 0.78
      rwventb = (0.308)*gmoi(igmrwb) ! 0.562825
      do mgs = 1,ngscnt
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( ipconc .ge. 3 ) THEN
         rwvent(mgs) = ventrx(mgs)*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)
        ELSE
         rwvent(mgs) =   &
     &  (rwventa + rwventb*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
        ENDIF
      ELSE
       rwvent(mgs) = 0.0
      ENDIF
      end do
!

!
!  Vapor Deposition constants
!
      do mgs = 1,ngscnt
      fvds(mgs) =    &
     &  (4.0*pi/rho0(mgs))*(ssi(mgs)-1.0)*   &
     &  (1.0/(fai(mgs)+fbi(mgs)))
      end do
      do mgs = 1,ngscnt
      fvce(mgs) =    &
     &  (4.0*pi/rho0(mgs))*(ssw(mgs)-1.0)*   &
     &  (1.0/(fav(mgs)+fbv(mgs)))
      end do
!
!
!
!  evaporation of rain
!
!
!
      do mgs = 1,ngscnt
!
      qrcev(mgs) = 0.0
      crcev(mgs) = 0.0
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
      rwcap(mgs) = (0.5)*xdia(mgs,lr,1)
      qrcev(mgs) = &
     &  fvce(mgs)*cx(mgs,lr)*rwvent(mgs)*rwcap(mgs)

! this line to allow condensation on rain:
      IF ( rcond .eq. 1 ) THEN
        qrcev(mgs) = min(qrcev(mgs), qxmxd(mgs,lv))
      ELSEIF ( rcond == 0 ) THEN
        qrcev(mgs) = 0.0
        crcev(mgs) = 0.0
! this line to have evaporation only:
      ELSE
        qrcev(mgs) = min(qrcev(mgs), 0.0)
      ENDIF
      qrcev(mgs) = max(qrcev(mgs), -qrmxd(mgs))
      
!      IF ( temg(mgs) .gt. 273.15 .and. qx(mgs,lr) .gt. 0.5e-3 ) THEN
!      IF ( kzbeg-1+kgs(mgs) <= 20 .and. qrcev(mgs) /= 0.0 ) THEN
!        write(iunit,*) 'qr,qrcev = ',igs(mgs),jgs,kzbeg-1+kgs(mgs),qx(mgs,lr)*1000.,Abs(qrcev(mgs)*dtp/qx(mgs,lr))*100.
!        write(iunit,*) 'fvce,cx,rwvent,rwcap = ',fvce(mgs),cx(mgs,lr),rwvent(mgs),rwcap(mgs)
!      ENDIF
      IF ( qrcev(mgs) < 0 ) crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)+1.e-20))*qrcev(mgs)
!      if ( temg(mgs) .lt. 273.15 ) qrcev(mgs) = 0.0
!      crcev(mgs) = (rho0(mgs)/(cx(mgs,lr)+1.e-20))*qrcev(mgs)
!      if ( temg(mgs) .lt. 273.15 ) crcev(mgs) = 0.0
!  
      ENDIF
      end do
!
!
!
!
!mpidebug
      IF ( ndebug .gt. 2 ) THEN
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: evap/cond - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'qrcev',9x,'crcev')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,qrcev(mgs),crcev(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: evap/cond - - - - - - - - - "
        write(0,"('mgs',5x,'qrcev',9x,'crcev')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,qrcev(mgs),crcev(mgs)
        end do
#endif
      ENDIF
!end mpidebug
! 
! 

!
!  vapor to cloud droplets   UP
!
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: dbg = 8'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 8'
#endif
!
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Collection: set 3-component'
!
!  time for riming....
!
!     rimtim = 240.0
!     dtrim = rimtim
!     xacrtim  = 120.0
!     tranfr = 0.50
!     tranfw = 0.50
!
!  coefficients for riming
!
!     rimc1 = 300.00
!     rimc2 = 0.44
!
! 
!
!  first sum all of the shed rain
!
!
!
!

!
!
!
!
      IF ( ipconc .ge. 1 ) THEN
!
!
!  concentration production terms
!
!  YYY
!
!
       DO mgs = 1,ngscnt
       pccwi(mgs) = 0.0
       pccwd(mgs) = 0.0
       pcrwi(mgs) = 0.0
       pcrwd(mgs) = 0.0
       ENDDO
!
!  Cloud water
!
      IF ( ipconc .ge. 2 ) THEN
      
      do mgs = 1,ngscnt
      pccwi(mgs) =  (0.0) ! + (1-il5(mgs))*(-cirmlw(mgs))
      pccwd(mgs) =  &
     &  - cautn(mgs) &
     &  - cracw(mgs)

      IF ( -pccwd(mgs)*dtp .gt. cx(mgs,lc) ) THEN
!       write(0,*) 'OUCH! pccwd(mgs)*dtp .gt. ccw(mgs) ',pccwd(mgs),cx(mgs,lc)
!       write(0,*) 'qc = ',qx(mgs,lc)
!       write(0,*) -ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)-cwfrzc(mgs)-cwctfzc(mgs)
!       write(0,*)  -cracw(mgs) -csacw(mgs)  -chacw(mgs)
!       write(0,*) - cautn(mgs)
       
       frac = -cx(mgs,lc)/(pccwd(mgs)*dtp)
       pccwd(mgs) = -cx(mgs,lc)/dtp
       
        cracw(mgs)   = frac*cracw(mgs)
        cautn(mgs)   = frac*cautn(mgs)
       
!       STOP
      ENDIF

      end do

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: cloud conc prod - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'cracw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,cracw(mgs),cautn(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: cloud conc prod - - - - - - - - - "
        write(0,"('mgs',5x,'cracw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,cracw(mgs),cautn(mgs)
        end do
#endif
!end mpidebug
      ENDIF

      ENDIF

!
!  Rain
!
      IF ( ipconc .ge. 3 ) THEN

      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 9'

      do mgs = 1,ngscnt
      pcrwi(mgs) = &
!     >   cracw(mgs) +  &
     &   crcnw(mgs)
      pcrwd(mgs) = &
     &  +crcev(mgs) &
     &  - cracr(mgs)

      IF ( -pcrwd(mgs)*dtp .gt. cx(mgs,lr) ) THEN
!       write(0,*) 'OUCH! pcrwd(mgs)*dtp .gt. crw(mgs) ',pcrwd(mgs)*dtp,cx(mgs,lr),mgs,igs(mgs),kgs(mgs)
!       write(0,*) -ciacr(mgs) 
!       write(0,*) -crfrz(mgs)
!       write(0,*) -chacr(mgs)
!       write(0,*)  crcev(mgs)
!       write(0,*)  -cracr(mgs)
       
       frac =  -cx(mgs,lr)/(pcrwd(mgs)*dtp)
       pcrwd(mgs) = -cx(mgs,lr)/dtp
        
        crcev(mgs) = frac*crcev(mgs)
        cracr(mgs) = frac*cracr(mgs)
       
!       STOP
      ENDIF
      end do

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: rain conc prod - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',5x,'cracw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,cracw(mgs),cautn(mgs)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: rain conc prod - - - - - - - - - "
        write(0,"('mgs',5x,'cracw',9x,'cautn')")
        do mgs = 1,ngscnt
          write(0,"(i3,2(1x,es12.5))") mgs,cracw(mgs),cautn(mgs)
        end do
#endif
!end mpidebug
      ENDIF

      ENDIF
!
!
!  Balance and checks for continuity.....within machine precision...
!
      do mgs = 1,ngscnt
      pctot(mgs)   = pccwi(mgs) +pccwd(mgs) + &
     &               pcrwi(mgs) +pcrwd(mgs) 
      end do
!
!
      ENDIF ! ( ipconc .ge. 1 )

      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 10'
!
!
!
!
!
!
!  GOGO
!  production terms
!
!
!
!  Vapor
!
      do mgs = 1,ngscnt
      pqwvi(mgs) =  &
     &  -Min(0.0, qrcev(mgs))
      pqwvd(mgs) =   &
     &  -Max(0.0, qrcev(mgs))
      end do
!
!  Cloud water
! 
      do mgs = 1,ngscnt
      pqcwi(mgs) =  (0.0)
      pqcwd(mgs) =  &
     &  -qracw(mgs) -qrcnw(mgs)

      IF ( -pqcwd(mgs)*dtp .gt. qx(mgs,lc) ) THEN
!       write(0,*) 'OUCH! pqcwd(mgs)*dtp .gt. qcw(mgs) ',pqcwd(mgs)*dtp,qx(mgs,lc)
!       write(0,*) il5(mgs)*(-qiacw(mgs))
!       write(0,*) il5(mgs)*(-qwfrzc(mgs))
!       write(0,*) il5(mgs)*(-qwctfzc(mgs))
!       write(0,*) -il5(mgs)*(qicichr(mgs))
!       write(0,*) -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs)
!       write(0,*) -il5(mgs)*(qwfrzp(mgs)+qwctfzp(mgs))

       frac = -qx(mgs,lc)/(pqcwd(mgs)*dtp)
       pqcwd(mgs) = -qx(mgs,lc)/dtp
       
        qracw(mgs)   = frac*qracw(mgs)
        qrcnw(mgs)   = frac*qrcnw(mgs)

!       STOP
      ENDIF

      end do
!
!  Rain
!
      do mgs = 1,ngscnt
      pqrwi(mgs) =   &
     &   qracw(mgs) +qrcnw(mgs) + Max(0.0, qrcev(mgs))
      pqrwd(mgs) =   &
     &  + Min(0.0,qrcev(mgs))
      
      IF ( -pqrwd(mgs)*dtp .gt. qx(mgs,lr) ) THEN
       write(0,*) 'OUCH! pqrwd(mgs)*dtp .gt. qrw(mgs) ',pqrwd(mgs),qx(mgs,lr)
       write(0,*)  qrcev(mgs)
       STOP
      ENDIF
      end do

!
!  rain reflectivity
!
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 11'

      IF ( lzr .gt. 1 ) THEN ! 
       
        DO mgs = 1,ngscnt
        
        zracw(mgs) = 0.0
        zracr(mgs) = 0.0
        zrcev(mgs) = 0.0
           
        IF ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) .gt. 0.0 ) THEN

          tmp = qx(mgs,lr)/cx(mgs,lr)
          g1 = g1x(mgs,lr) ! 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)


        IF ( qracw(mgs) > 0.0 .and. cx(mgs,lr) > 0.0 ) THEN
         zracw(mgs) =  g1x(mgs,lr)*(rho0(mgs)/(1000.))**2*( 2.*tmp * qracw(mgs) )
        ENDIF
        
        IF ( cracr(mgs) > 0.0 .and. cx(mgs,lr) > 0.0  ) THEN
         zracr(mgs) =  g1x(mgs,lr)*(rho0(mgs)/(1000.))**2*( tmp**2 * cracr(mgs) )
        ENDIF

        qtmp = qrcev(mgs)
        ctmp = crcev(mgs)
        
        zrcev(mgs) = g1x(mgs,lr)*(rho0(mgs)/(xdn(mgs,lr)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )
        zrcev(mgs) = Max( zrcev(mgs), -zxmxd(mgs,lr) )
        
        ENDIF

         pzrwi(mgs) = zrcnw(mgs) + zracw(mgs) + zracr(mgs) + Max(0.,zrcev(mgs) )

         pzrwd(mgs) = 0.0  + Min(0.,zrcev(mgs) )

      IF (  ny <= 2 .and.                      &
     &         zx(mgs,lr) + dtp*(pzrwi(mgs)+pzrwd(mgs))  <= 0.0  &
     &         .and. qx(mgs,lr) > 0.01e-3 ) THEN
           
!           g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           g1 = g1x(mgs,lr) ! 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)

           
           write(0,*) 'zrw negative: pzrwi,pzrwd = ',dtp*pzrwi(mgs),dtp*pzrwd(mgs),igs(mgs),kgs(mgs)
           write(0,*) 'qr,cr = ',1.e3*qx(mgs,lr),cx(mgs,lr),an(igs(mgs),jgs,kgs(mgs),lr),an(igs(mgs),jgs,kgs(mgs),lnr)
           write(0,*) 'zx,z,z0 = ',zx(mgs,lr),g1*rho0(mgs)**2*qx(mgs,lr)**2/(xdn(mgs,lr)**2*cx(mgs,lr)), &
     &      an(igs(mgs),jgs,kgs(mgs),lzr)
           write(0,*) 'alpha_r,ventrx = ', alpha(mgs,lr),ventrx(mgs)
           write(0,*) 'zrcev, q ',   dtp*Min(0.,zrcev(mgs) ), dtp*Min(0.0,qrcev(mgs)), dtp*Min(0.0,crcev(mgs))
           write(0,*) 'zrcnw,zracr,zracw = ',dtp*zrcnw(mgs),dtp*zracr(mgs),dtp*zracw(mgs)
           write(0,*) 'temp, qc, rho0 = ',temcg(mgs),1.e3*qx(mgs,lc),rho0(mgs)
           write(0,*) 'tmp,qtmp,ctmp,g1 = ',tmp,qtmp,ctmp,g1
           write(0,*) 'xdia = ',xdia(mgs,lr,3)

          IF ( zx(mgs,lr) > 0.0 ) THEN
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            vr = xv(mgs,lr)
!            z = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z,xv
           DO i = 1,10
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           print*,'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
          ENDIF
         ENDIF

         ENDIF ! ny <= 2
      
        ENDDO
      
      ENDIF

      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 12'
      
!
      IF ( .false. .and. ny .eq. 2 ) THEN
      DO mgs=1,ngscnt
      IF (( pqcwi(mgs) .ne. 0.0 .or. pqcwd(mgs) .ne. 0.0 .or. pccwd(mgs) .ne. 0.0) &
     &   .and. igs(mgs)+ixbeg-1 .lt. 25) THEN
      write(iunit,*)   'Cloud water ',igs(mgs)+ixbeg-1,kgs(mgs), qx(mgs,lc), cx(mgs,lc)
!
      write(iunit,*)   'pqcwi =',pqcwi(mgs) 
      write(iunit,*)   'qracw = ',-qracw(mgs), qx(mgs,lr),qx(mgs,lc)
      write(iunit,*)   'qrcnw = ',-qrcnw(mgs) 
      write(iunit,*)   'pqcwd = ',pqcwd(mgs) 


      write(iunit,*)  'Concentration:'
      write(iunit,*)   'cautn = ', -cautn(mgs) 
      write(iunit,*)   'cracw = ',-cracw(mgs), cx(mgs,lr), cx(mgs,lc)
      write(iunit,*)   'pccwd = ',pccwd(mgs) 
!      ENDIF

!      IF ( pqrwi(mgs) .ne. 0.0 .or. pqrwd(mgs) .ne. 0.0 .or.  
!     :     pcrwi(mgs) .ne. 0.0 .or. pcrwd(mgs) .ne. 0.0) THEN
      write(iunit,*)      'Rain ',igs(mgs)+ixbeg-1,kgs(mgs), qx(mgs,lr), cx(mgs,lr)
!
      write(iunit,*)    'qracw = ',  qracw(mgs)
      write(iunit,*)    'qrcnw = ', qrcnw(mgs)
      write(iunit,*)    'pqrwi = ', pqrwi(mgs)    
      write(iunit,*)    'qrcev = ', qrcev(mgs)
      write(iunit,*)    'pqrwd = ', pqrwd(mgs) 
!
      
      write(iunit,*)  'Rain concentration'
      write(iunit,*)  'pcrwi = ', pcrwi(mgs) 
      write(iunit,*)  'crcnw = ', crcnw(mgs)
      write(iunit,*)  'pcrwd = ', pcrwd(mgs) 
      write(iunit,*)  'crcev = ', crcev(mgs)
      write(iunit,*)   'cracr = ', cracr(mgs)
      
      ENDIF
      
      ENDDO
      ENDIF
!
!
      if ( nstep/1*1 .eq. nstep ) then
      do mgs = 1,ngscnt
!
      ptotal(mgs) = 0.
      ptotal(mgs) = ptotal(mgs)  &
     &  + pqwvi(mgs) + pqwvd(mgs) &
     &  + pqcwi(mgs) + pqcwd(mgs) &
     &  + pqrwi(mgs) + pqrwd(mgs)
!      

!      if ( ( ndebug .ge. 1 .and. abs(ptotal(mgs)) .gt. eqtot )
      if ( (  abs(ptotal(mgs)) .gt. eqtot ) &
!     :    .or. pqswi(mgs)*dtp .gt. 1.e-3
!     :    .or. pqhwi(mgs)*dtp .gt. 1.e-3
!     :     .or. dtp*(pqrwi(mgs)+pqrwd(mgs)) .gt. 10.0e-3 
!     :     .or. dtp*(pccii(mgs)+pccid(mgs)) .gt. 1.e7 
!     :     .or. dtp*(pcipi(mgs)+pcipd(mgs)) .gt. 1.e7  &
     &  .or.  .not. (ptotal(mgs) .lt. 1.0 .and. &
     &            ptotal(mgs) .gt. -1.0)    ) then
      write(iunit,*) 'YIKES! ','ptotal1',nstep,mgs,igs(mgs),jgs, &
     &       kgs(mgs),ptotal(mgs)
     
      write(iunit,*) 't7: ', t7(igs(mgs),jgs,kgs(mgs))
      write(iunit,*)  'cci,ccw,crw,rdia: ',cx(mgs,lc),cx(mgs,lr),0.5*xdia(mgs,lr,1)
      write(iunit,*)  'qc,qi,qr : ',qx(mgs,lc),qx(mgs,lr)
      write(iunit,*)  'rmas, qrcalc : ',xmas(mgs,lr),xmas(mgs,lr)*cx(mgs,lr)/rho0(mgs)
      write(iunit,*)  'vti,vtc,vtr: ',vtxbar(mgs,lc,1),vtxbar(mgs,lr,1)
      write(iunit,*)  'cidia,cwdia,qcmxd: ',xdia(mgs,lc,1),qcmxd(mgs)
      
      
      
      
      write(iunit,*)  'rain cx,xv : ',cx(mgs,lr),xv(mgs,lr)
      write(iunit,*)  'temcg = ', temcg(mgs)
      write(iunit,*)  pqwvi(mgs) ,pqwvd(mgs)
      write(iunit,*)  pqcwi(mgs) ,pqcwd(mgs)
      write(iunit,*)  pqrwi(mgs) ,pqrwd(mgs)
      write(iunit,*) 'END OF OUTPUT OF SOURCE AND SINK'

!
!  print production terms
!
      write(iunit,*)   'Vapor'
!
      write(iunit,*)   -qrcev(mgs)
      write(iunit,*)    pqwvi(mgs) 
      write(iunit,*)    pqwvd(mgs)
!
!
!
!
!
      write(iunit,*)   'Cloud water'
!
      write(iunit,*)   pqcwi(mgs) 
      write(iunit,*)   -qracw(mgs)
      write(iunit,*)   -qrcnw(mgs) 
      write(iunit,*)   pqcwd(mgs) 


      write(iunit,*)  'Concentration:'
      write(iunit,*)   -cautn(mgs) 
      write(iunit,*)   -cracw(mgs)
      write(iunit,*)   pccwd(mgs) 
!
      write(iunit,*)      'Rain '
!
      write(iunit,*)      qracw(mgs)
      write(iunit,*)      qrcnw(mgs)
      write(iunit,*)        pqrwi(mgs)    
      write(iunit,*)        qrcev(mgs)
      write(iunit,*)        pqrwd(mgs) 
!
      
      write(iunit,*)  'Rain concentration'
      write(iunit,*)  pcrwi(mgs) 
      write(iunit,*)    crcnw(mgs)
      write(iunit,*)  pcrwd(mgs) 
      write(iunit,*)   +crcev(mgs)
      write(iunit,*)   cracr(mgs)
!      write(iunit,*)   -il5(mgs)*ciracr(mgs)


!
!  Balance and checks for continuity.....within machine precision...
!
!
      write(iunit,*) 'END OF OUTPUT OF SOURCE AND SINK'
      write(iunit,*) 'PTOTAL',ptotal(mgs)
!
      end if
!
      end do
!

      DO mgs = 1,ngscnt
      IF ( .not. (ptotal(mgs) .lt. eqtot .and. &
     &            ptotal(mgs) .gt. -eqtot) ) THEN
! cmic$ guard
! C$PAR CRITICAL SECTION    
!! c$omp  critical
        iptotal = iptotal + 1
!        write(iunit,*) igs(mgs),jy,kgs(mgs),ptotal(mgs)
        ptotalmx = Max( ptotalmx, ptotal(mgs) )
        ptotalmn = Min( ptotalmn, ptotal(mgs) )
!! c$omp  end critical
! C$PAR END CRITICAL SECTION   
! cmic$ endguard
      END IF
      END DO

      end if ! ( nstep/12*12 .eq. nstep )
!
!  latent heating from phase changes (except qcw, qci cond, and evap)
!
      do mgs = 1,ngscnt
      pvap(mgs) =  &
     &   qrcev(mgs)
      ptem(mgs) =  &
     &  (cpi/pi0(mgs))* &
     &  (felv(mgs)*pvap(mgs))
      thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)
      end do
!

!  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
!
!
      do mgs = 1,ngscnt
      qwvp(mgs) = qwvp(mgs) +      &
     &   dtp*(pqwvi(mgs)+pqwvd(mgs)) 
      qx(mgs,lc) = qx(mgs,lc) + &
     &   dtp*(pqcwi(mgs)+pqcwd(mgs)) 
!      IF ( qx(mgs,lr) .gt. 10.0e-3 )  THEN
!       print*, 'RAIN1a: ',igs(mgs),kgs(mgs),qx(mgs,lr)
!      ENDIF
      qx(mgs,lr) = qx(mgs,lr) + &
     &   dtp*(pqrwi(mgs)+pqrwd(mgs)) 
!      IF ( qx(mgs,lr) .gt. 10.0e-3 ) THEN
!        print*, 'RAIN1b: ',igs(mgs),kgs(mgs),qx(mgs,lr)
!        print*, pqrwi(mgs),pqrwd(mgs)
!       ENDIF
      end do
!
!
!
! concentrations
!
      DO mgs = 1,ngscnt
      IF ( ipconc .ge. 2 ) THEN
      cx(mgs,lc) = cx(mgs,lc) + &
     &   dtp*(pccwi(mgs)+pccwd(mgs)) 
      ENDIF
      IF ( ipconc .ge. 3 ) THEN
      cx(mgs,lr) = cx(mgs,lr) + &
     &   dtp*(pcrwi(mgs)+pcrwd(mgs)) 
      ENDIF
      end do

!
! reflectivity
!
       IF ( lzr .gt. 1 ) THEN
       DO mgs = 1,ngscnt
       zx(mgs,lr) = zx(mgs,lr) +    &
     &   dtp*(pzrwi(mgs)+pzrwd(mgs)) 
       ENDDO
       ENDIF
!
!
!
! start saturation adjustment
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 30a'
!      include 'sam.jms.satadj.sgi'
!
!
      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: start of sat adj - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: start of sat adj - - - - - - - - - "
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr)
        end do
#endif
!end mpidebug
      ENDIF



      DO mgs = 1,ngscnt
!       tq1(igs(mgs),jy,kgs(mgs),1) = qsmlr(mgs)
       thproc(kzbeg-1+kgs(mgs),1) = thproc(kzbeg-1+kgs(mgs),1) + pcrwi(mgs)
       thproc(kzbeg-1+kgs(mgs),2) = thproc(kzbeg-1+kgs(mgs),2) + pcrwd(mgs)
       thproc(kzbeg-1+kgs(mgs),3) = thproc(kzbeg-1+kgs(mgs),3) + crcnw(mgs)
       thproc(kzbeg-1+kgs(mgs),4) = thproc(kzbeg-1+kgs(mgs),4) + crcev(mgs)
       thproc(kzbeg-1+kgs(mgs),5) = thproc(kzbeg-1+kgs(mgs),5) + cracr(mgs)
       thproc(kzbeg-1+kgs(mgs),6) = thproc(kzbeg-1+kgs(mgs),6) + qracw(mgs)
       thproc(kzbeg-1+kgs(mgs),7) = thproc(kzbeg-1+kgs(mgs),7) + qrcnw(mgs)
       thproc(kzbeg-1+kgs(mgs),8) = thproc(kzbeg-1+kgs(mgs),8) + qrcev(mgs)
      ENDDO
!     >   qracw(mgs) +qrcnw(mgs) + Max(0.0, qrcev(mgs))

!
!  scatter precipitation fluxes, and thetap, and hydrometeors
!
!
!
!
!
!DIR$ IVDEP
      do mgs = 1,ngscnt
      t0(igs(mgs),jy,kgs(mgs)) =  temg(mgs)
      end do
!

      do mgs = 1,ngscnt
!
      an(igs(mgs),jy,kgs(mgs),lt) =  &
     &  ab(kgs(mgs),lt) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) =  &
     &  ab(kgs(mgs),lv) + qwvp(mgs) 
!
      
      DO il = lc,lhab
         an(igs(mgs),jy,kgs(mgs),il) = qx(mgs,il) + &
     &     min( an(igs(mgs),jy,kgs(mgs),il), 0.0 )
         qx(mgs,il) = an(igs(mgs),jy,kgs(mgs),il)
      ENDDO


!
      end do
!
      IF ( lzr > 1 ) THEN
      if (ndebug .gt. 0 ) write(0,*) 'put back ZRW'
!
!  6th moments
!
       DO mgs = 1,ngscnt
         an(igs(mgs),jy,kgs(mgs),lzr) = zx(mgs,lr) +   &
     &     min( an(igs(mgs),jy,kgs(mgs),lzr), 0.0 )
         zx(mgs,lr) = an(igs(mgs),jy,kgs(mgs),lzr)
       ENDDO
      ENDIF


      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: sat adj an start - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: sat adj an start - - - - - - - - - "
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr)
        end do
#endif
!end mpidebug
      ENDIF


      if ( ipconc .ge. 1 ) then
!DIR$ IVDEP
      DO il = lc,lhab

        IF ( ipconc .ge. ipc(il) ) THEN
          DO mgs = 1,ngscnt
            an(igs(mgs),jy,kgs(mgs),ln(il)) = Max(cx(mgs,il), 0.0)
          ENDDO
        ENDIF
      ENDDO

      IF (  ipconc .ge. 2 ) THEN
      do mgs = 1,ngscnt
        an(igs(mgs),jy,kgs(mgs),lss) = ssmax(mgs) 
        an(igs(mgs),jy,kgs(mgs),lccn) = Min(ccwmx,ccnc(mgs))
      end do
      ENDIF

      end if 
!
!
!
      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: after scatter - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',7x,'lt',12x,'il',11x,'ss',11x,'ccn')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lt), an(igs(mgs),jy,kgs(mgs),il), &
     &                                    an(igs(mgs),jy,kgs(mgs),lss),an(igs(mgs),jy,kgs(mgs),lccn)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: after scatter - - - - - - - - - "
        write(0,"('mgs',7x,'lt',12x,'il',11x,'ss',11x,'ccn')")
        do mgs = 1,ngscnt
         write(0,"(i3,4(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lt), an(igs(mgs),jy,kgs(mgs),il), &
     &                                    an(igs(mgs),jy,kgs(mgs),lss),an(igs(mgs),jy,kgs(mgs),lccn)
        end do
#endif
!end mpidebug
      ENDIF
!
!
!
      deallocate ( alpha )
      deallocate ( qx )
      deallocate ( qxw )
      deallocate ( cx )
      deallocate ( cxmxd )
      deallocate ( qxmxd )
      deallocate ( xv )
      deallocate ( vtxbar )
      deallocate ( xmas )
      deallocate ( xdn )
      deallocate ( xdia )
      deallocate ( zx )
      deallocate ( zxmxd )
      deallocate ( g1x )

      
 9998 continue  !! ( ngscnt .eq. 0 )

#ifdef MPI
      if ( kzbeg-1+kz .gt. nzend-kstag-1 .and. ixbeg-1+ix .gt. nxend-istag ) then
!      if ( kz .gt. nz-kstag-1 .and. ix .ge. nx-istag) then
        if ( ixbeg-1+ix .eq. nxend-1 ) then
         go to 1200
        elseif ( ix .ge. nx ) then
         go to 1200
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz 
      end if
#else
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
        go to 1200
      else
        nzmpb = kz 
      end if
#endif

#ifdef MPI
      if ( ix .ge. nx-1 ) then
       if ( ixbeg-1+ix .eq. nxend-1 ) then
        nxmpb = 1
       elseif ( ix .ge. nx ) then
        nxmpb = 1
       else
        nxmpb = ix+1
       endif
      else
       nxmpb = ix+1
      end if
#else
      if ( ix+1 .gt. nx-1 ) then
       nxmpb = 1
      else
       nxmpb = ix+1
      end if
#endif

 1000 continue  !! inumgs
 1200 continue
!
!  end of gather scatter (for this jy slice)
!
!  precipitation fallout contributions (if itfall .eq. 3 ) then
!
!
!  Do fallout stuff now if itfall=3
!
      IF ( itfall .eq. 3 ) THEN

!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0
      
#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze ; DO ix = ixb,ixe
#else
      do kz = 1,nz-1
      DO ix = 1,nx-1
#endif
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

! set up vt arrays:
!      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
!     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)

      call ziegfall(nx,ny,nz,nor,nor,na,dtp,dz,jgs,ngs, &
     &  xvt, &
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, &
     &  rwmasn,rwmasx,cwradn, &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn, &
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

! mixing ratio      
      DO il = lc,lhab
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), &
     &  an,db1,imapz,mzdist,il,1,xfall)
      ENDDO
      




      
      IF ( ipconc .ge. 1 ) THEN
!      infall = 2
      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN
         IF ( .not. (lnr .ge. 1 .and. ln(il) .eq. lnr ) ) THEN
         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
         ELSE
           IF ( itfall .eq. 1 ) THEN
            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,2,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
           ELSE
            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il), &
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
           ENDIF
         ENDIF
        ENDIF
      ENDDO
      ENDIF

      
      ENDIF ! itfall.eq.3
!
!
!
! Maximum courant number is dt*vtmax.  If this is greater than 1 then
! need to split the fallout term into smaller steps.
!
!   
!
!
!  end of jy loop
!
 3999 continue
!
!
 9999 continue
!
!
!

!       IF ( fallonly ) GOTO 4999
!
!
!  Ziegler nucleation 
!
      IF ( ipconc .ge. 2 .and. iba .eq. 1 ) THEN

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: zieg nucl start - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('kz',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',10x,'temp',10x,'lss',10x,'ss')")
        do kz = 1,nz
         write(0,"(i3,7(1x,es12.5))") kz,an(21,22,kz,lc), an(21,22,kz,lr), &
     &                                    an(21,22,kz,lnc),an(21,22,kz,lnr), &
     &                                    an(21,22,kz,lt),an(21,22,kz,lss),ssat(21,22,kz)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: zieg nucl start - - - - - - - - - "
        write(0,"('kz',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',10x,'temp',10x,'lss',10x,'ss')")
        do kz = 1,nz
         write(0,"(i3,7(1x,es12.5))") kz,an(21,22,kz,lc), an(21,22,kz,lr), &
     &                                    an(21,22,kz,lnc),an(21,22,kz,lnr), &
     &                                    an(21,22,kz,lt),an(21,22,kz,lss),ssat(21,22,kz)
        end do
#endif
!end mpidebug
      ENDIF

#ifdef MPI
      IF ( .false. .and. number_of_processes .gt. 1 ) THEN

       CALL cld_cpu('MPICOM-MICRO')

       nampi = 2
       nb = 1


       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,   &
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,an(-ng+1,-ng+1,-ng+1,lt) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,an(-ng+1,-ng+1,-ng+1,lt))

       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,   &
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,an(-ng+1,-ng+1,-ng+1,lt))

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,   &
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,an(-ng+1,-ng+1,-ng+1,lt))

        CALL cld_cpu('MPICOM-MICRO')

       ENDIF

#endif

        ssat(:,:,:) = 0.0
      ssfilt(:,:,:) = 0.0

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe


!         t7(ix,jy,kz) = ((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
         temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz)
!     >     *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
        
!        temp1 = 
!     :    theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)) )
        
         t0(ix,jy,kz) = temp1
         
         ltemq = Int( (temp1-163.15)/fqsat+1.5 )
      IF ( ltemq .lt. 1 .or. ltemq .gt. nqsat ) THEN
        write(iunit,*) 'out of range ltemq!',temp1, &
     &      an(ix,jy,kz,lt),t77(ix,jy,kz),t00(ix,jy,kz), &
!     :       thetap(mgs),theta0(mgs),pres(mgs),theta(mgs), &
     &      ltemq,ix,jy,kz 
!        write(iunit,*) an(igs(mgs),jy,kgs(mgs),lt),
!     :   ab(igs(mgs),jy,kgs(mgs),lt),
!     :   t0(igs(mgs),jy,kgs(mgs))
!        write(iunit,*) fcc3(mgs),qx(mgs,lc),qitmp(mgs),dtp,ptem(mgs)
        ltemq = Min( nqsat, Max(1,ltemq) )
!        STOP
! C$PAR END CRITICAL SECTION
      END IF

          c1 = t00(ix,jy,kz)*tabqvs(ltemq)

          ssat(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
                    
!          ssfilt(ix,jy,kz) = ssat(ix,jy,kz)

!#ifdef MPI
!         IF ( nyend .eq. 2 .and. ssfilt(ix,jy,kz) .gt. 1.0 &
!#else
!         IF ( ny .eq. 2 .and. ssfilt(ix,jy,kz) .gt. 1.0 &
!#endif
!     &                                .and. ncdebug .ge. 1 ) THEN
!          write(iunit,*) 'ssfilt1 = ',ssfilt(ix,jy,kz), & 
!     &     an(ix,jy,kz,lv)*1.e3,c1*1.e3,an(ix,jy,kz,lt),ix,kz, & 
!     &     an(ix,jy,kz,lnc)*1.e-6,an(ix,jy,kz,lc)*1.e3
!
!         ENDIF

          IF ( lsat >= 1 )  axtra(ix,jy,kz,lsat)  = ssat(ix,jy,kz)
!          IF ( lsati >= 1 ) axtra(ix,jy,kz,lsati) = ssati(ix,jy,kz)


      IF ( ndebug .gt. 2 ) THEN
!mpidebug
       if (ix == 21 .and. jy == 22) then
#ifdef MPI
        if (my_rank==0) then
         if(kz==1) write(0,*) " - - - - - - - - - WARMZIEG: middle of ssfilt fill - - - - - - - - - "
         if(kz==1) write(0,*) "my_rank=",my_rank
         if(kz==1) write(0,"('kz',7x,'lt',11x,'t77',8x,'ltemq',10x,'t00',7x,'tabqvs',8x,'c1',12x,'lv',12x,'ss',8x,'ssfilt')")
         write(0,"(i3,9(1x,es12.5))") kz,an(ix,jy,kz,lt), t77(ix,jy,kz), &
     &                                    ltemq,t00(ix,jy,kz),tabqvs(ltemq),c1, &
     &                                    an(ix,jy,kz,lv),ssat(ix,jy,kz),ssfilt(ix,jy,kz)
        end if
#else
         if(kz==1) write(0,*) " - - - - - - - - - WARMZIEG: middle of ssfilt fill - - - - - - - - - "
         if(kz==1) write(0,"('kz',7x,'lt',11x,'t77',8x,'ltemq',10x,'t00',7x,'tabqvs',8x,'c1',12x,'lv',12x,'ss',8x,'ssfilt')")
         write(0,"(i3,9(1x,es12.5))") kz,an(ix,jy,kz,lt), t77(ix,jy,kz), &
     &                                    ltemq,t00(ix,jy,kz),tabqvs(ltemq),c1, &
     &                                    an(ix,jy,kz,lv),ssat(ix,jy,kz),ssfilt(ix,jy,kz)
#endif
        end if
!end mpidebug
       ENDIF
       
        ENDDO
       ENDDO
      ENDDO


#ifdef MPI

      IF ( irenuc3d > 0 .or. issfilt > 0 .or. nprock > 1 ) THEN
      IF ( number_of_processes .gt. 1 ) THEN

       CALL cld_cpu('MPICOM-MICRO')

       nampi = 1
       nb = 1
       IF ( issfilt == 1 ) nb = 2


       IF ( irenuc3d > 0 .or. issfilt > 0 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,ssat(-ng+1,-ng+1,-ng+1) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,ssat(-ng+1,-ng+1,-ng+1))
     
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,ssat(-ng+1,-ng+1,-ng+1))

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,ssat(-ng+1,-ng+1,-ng+1))

       ENDIF

       IF ( nprock > 1 ) THEN
       
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,ssat(-ng+1,-ng+1,-ng+1))

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,ssat(-ng+1,-ng+1,-ng+1))

       
       nb = 3
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,t0(-ng+1,-ng+1,-ng+1))

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,    &
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,t0(-ng+1,-ng+1,-ng+1))

        ENDIF
        
        CALL cld_cpu('MPICOM-MICRO')

       ENDIF
       ENDIF

#endif

      IF ( issfilt .eq. 1 ) THEN

      kzb = 0
      kze = ktile+1
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe

!
! filter supersaturation field
!
        IF ( ipconc .ge. 2 .and. ixbeg-1+ix .gt. 1 .and. ixbeg-1+ix .lt. nxend-1  & 
     &       .and. kzbeg-1+kz .gt. 1 .and. kzbeg-1+kz .lt. nzend-1  ) THEN
     
#ifdef MPI
!         write(0,*) "ICEZVD_DR: SS FILTER NEEDS TO BE TESTED FOR MPI"
!         write(0,*) 'need to change code to use im1,ip1 etc.  STOP'
!         STOP
#endif
         ssfilt(ix,jy,kz)=(0.25* & 
     &       (ssat(ix+1,jy,kz) + ssat(ix-1,jy,kz) + &
#ifdef MPI 
     &        ssat(ix,jy+1,kz) + ssat(ix,jy-1,kz) + &
#else
     &        ssat(ix,Min(nyend-1,jy+1),kz) + ssat(ix,Max(1,jyend-1),kz) + &
#endif
     &        ssat(ix,jy,kz+1) + ssat(ix,jy,kz-1))  & 
     &        + 0.5*ssat(ix,jy,kz))/2.
        ELSE
         ssfilt(ix,jy,kz) = ssat(ix,jy,kz)
        ENDIF
!
      END DO
      END DO
      END DO
      
      ELSE
      
       ssfilt(:,:,:) = ssat(:,:,:)
      
      ENDIF
      
!
      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      do 29999 jy = jyb,jye
!
!  VERY IMPORTANT:  SET jgs = jy
!
      jgs = jy
      
!
!..Gather microphysics  
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: Gather stage 2'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Gather stage 2'
#endif
      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1

      do 2000 inumgs = 1,numgs 
      ngscnt = 0

#ifdef MPI
      kzb = nzmpb
      kze = ktile
!      if (kzbeg .le. nzmpb .and. kzend .gt. nzmpb) kzb = nzmpb-kzbeg+1
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      ixb = nxmpb
      ixe = itile
!      if (ixbeg .le. nxmpb .and. ixend .gt. nxmpb) ixb = nxmpb-ixbeg+1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do 2001 kz = kzb,kze
      do 2002 ix = nxmpb,ixe
#else
      do 2001 kz = nzmpb,nz-kstag-1 
      do 2002 ix = nxmpb,nx-istag
#endif
      pqs(kz) = 380.0/(pn(ix,jy,kz)+pb(kz))
      theta(kz) = an(ix,jy,kz,lt) 
      temg(kz) = t0(ix,jy,kz) 

      temcg(kz) = temg(kz) - tfr
      tqvcon = temg(kz)-cbw
      ltemq = (temg(kz)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(kz) = pqs(kz)*tabqvs(ltemq)

      qss(kz) = qvs(kz)
      if ( temg(kz) .lt. tfr ) then
      end if
!
      if ( an(ix,jy,kz,lv)  .gt. qss(kz) .or. &
     &     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.  &
     &     an(ix,jy,kz,lr)  .gt. qxmin(lr)        ) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
#ifdef MPI
      IF ( nyend .eq. 2 .and. ncdebug .ge. 1 ) THEN
#else
      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
#endif
       write(iunit,*) 'add point ',ix+ixbeg-1,kz,ssat(ix,jy,kz), &
     &  an(ix,jy,kz,lv)*1.e3,qss(kz)*1.e3, an(ix,jy,kz,lc)*1.e3, &
     &  an(ix,jy,kz,lnc)*1.e-6
      ENDIF
      if ( ngscnt .eq. ngs ) goto 2100
      end if
 2002 continue !! ix
      nxmpb = 1
 2001 continue !! kz
!      if ( jy .eq. (ny-jstag) ) iend = 1
 2100 continue
      if ( ngscnt .eq. 0 ) go to 29998
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: dbg=5'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg=5'
#endif      
      allocate ( qx(ngscnt,lv:lhab) )
      allocate ( cx(ngscnt,lc:lhab) )
      allocate ( xv(ngscnt,lc:lhab) )
      allocate ( vtxbar(ngscnt,lc:lhab,3) )
      allocate ( xmas(ngscnt,lc:lhab) )
      allocate ( xdn(ngscnt,lc:lhab) )
      allocate ( xdia(ngscnt,lc:lhab,3) )
      allocate (  alpha(ngscnt,lr:lhab) )
      allocate (  zx(ngscnt,lr:lhab) )
      
      qx(:,:) = 0.0
      cx(:,:) = 0.0
      zx(:,:) = 0.0

      xv(:,:) = 0.0
      xmas(:,:) = 0.0

      alpha(:,lr) = xnu(lr)

!
!  define temporaries for state variables to be used in calculations
!
      DO mgs = 1,ngscnt
      
      qx(mgs,lv) = an(igs(mgs),jy,kgs(mgs),lv)
      qx(mgs,lc) = max(an(igs(mgs),jy,kgs(mgs),lc), 0.0) 
      qx(mgs,lr) = max(an(igs(mgs),jy,kgs(mgs),lr), 0.0) 

       temp0(mgs) = tempbz(kgs(mgs))
      theta0(mgs) = ab(kgs(mgs),lt)
      thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt) - ab(kgs(mgs),lt)
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      pres0(mgs) = pb(kgs(mgs))
      qv0(mgs) = ab(kgs(mgs),lv)
      qwvp(mgs) = qx(mgs,lv) - qv0(mgs) ! an(igs(mgs),jy,kgs(mgs),lv) 
!     > - ab(kgs(mgs),lv)

!       theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
       presp(mgs) = pn(igs(mgs),jy,kgs(mgs))
       pres(mgs) = presp(mgs) + pres0(mgs)
       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhoinv(mgs) = 1.0/rho0(mgs)
       pi0(mgs) = pinit(kgs(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
       pk(mgs) = t77(igs(mgs),jy,kgs(mgs)) ! ( pres(mgs) / poo ) ** cap
!       temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap  ! updated temperature
!       temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
       temcg(mgs) = temg(mgs) - tfr
       qss0(mgs) = (380.0)/(pres(mgs))
       pqs(mgs) = (380.0)/(pres(mgs))
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
!
        qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
        es(mgs) = 6.1078e2*tabqvs(ltemq)
        qss(mgs) = qvs(mgs)


        temgx(mgs) = min(temg(mgs),313.15)
        temgx(mgs) = max(temgx(mgs),233.15)
        felv(mgs) = 2500837.367  &
     &   * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))

       il5(mgs) = 0
       if ( temg(mgs) .lt. tfr ) then 
        il5(mgs) = 1
       end if

      
      ENDDO

!
! load concentrations
!
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        IF ( lccn .gt. 1 ) THEN
          ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        ENDIF
       end do
      end if
      
        cnuc(1:ngscnt) = cwccn*(1. - renucfrac) + ccnc(1:ngscnt)*renucfrac
      
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
       end do
      end if

!  Find shape parameter rain

      ventrx(:) = ventr
      
      IF ( lzr > 1 .and. rcond == 2 ) THEN ! { RAIN SHAPE PARAM
      CALL cld_cpu('Z-MOMENT-1')  
          il = lr
          DO mgs = 1,ngscnt

         zx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lzr), 0.0)

         IF ( iresetmoments == 1 .or. iresetmoments == il  ) THEN
         IF ( zx(mgs,il) <= zxmin ) THEN
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ELSEIF ( cx(mgs,il) <= cxmin ) THEN
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
         ENDIF
         
         IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

          xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
          IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
!            xv(mgs,lr) = xvmx(lr)
!            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
          ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
            xv(mgs,lr) = xvmn(lr)
            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z1   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(chw*1000.*1000)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z1   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
          ENDIF
        
          IF ( zx(mgs,lr) > 0.0 ) THEN
            vr = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
!            z1 = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z1 = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z1*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z1 .gt. 0.0 ) THEN
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z1,xv
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           print*,'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

! check for artificial breakup (rain larger than allowed max size)
        IF (  xv(mgs,il) .gt. xvmx(il) ) THEN
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          IF ( tmp < cx(mgs,il) ) THEN ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z1 = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
           IF ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) THEN

            IF ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alphar .and. alp <= rnumin ) THEN

             z1  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z1
             ENDIF
           ENDIF

           tmp = alpha(mgs,lr) + 4./3.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma(alpha(mgs,lr) + 1.)
           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))
           
           ENDIF
          ENDIF
          
          ENDIF
          
          ENDDO
        CALL cld_cpu('Z-MOMENT-1')  
        ENDIF ! }
      
      DO mgs = 1,ngscnt
      
      uvel(mgs) = (0.5)*(u(igs(mgs),jgs,kgs(mgs)) &
     &                  +u(igs(mgs)+1,jgs,kgs(mgs)))
      IF ( ny .gt. 2 ) THEN
       vvel(mgs) = (0.5)*(v(igs(mgs),jgs,kgs(mgs)) &
     &                  +v(igs(mgs),jgs+1,kgs(mgs)))
      ELSE
       vvel(mgs) = 0.0
      ENDIF
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)+1) &
     &                  +w(igs(mgs),jgs,kgs(mgs)))

#ifdef MPI
      if (kzbeg-1+kgs(mgs) .eq. 1) then
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                    +w(igs(mgs),jgs,1))
      else
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                    +w(igs(mgs),jgs,kgs(mgs)-1))
      endif
#else
      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) &
     &                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))
#endif
      IF ( ssfilt(igs(mgs),jgs,kgs(mgs)) .gt. 0.0 ) THEN
        ssbar(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs))
      ELSE
        ssbar(mgs) = ssat(igs(mgs),jgs,kgs(mgs))
      ENDIF

#ifdef MPI
      ssat0(mgs)  = ssat(igs(mgs),jgs,kgs(mgs))
      ssf(mgs)    = ssfilt(igs(mgs),jgs,kgs(mgs))
      
      ssfkp1(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs)+1)
      if (kzend.ge.nzend-1) ssfkp1(mgs) = ssfilt(igs(mgs),jgs,Min(nz-1,kgs(mgs)+1))
      ssfkm1(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs)-1)
      IF ( kzbeg == nzbeg )  ssfkm1(mgs) = ssfilt(igs(mgs),jgs,Max(1,kgs(mgs)-1))
      
      ssfjp1(mgs) = ssfilt(igs(mgs),jgs+1,kgs(mgs))
      if (jyend.eq.nyend) ssfjp1(mgs) = ssfilt(igs(mgs),Min(ny-1,jgs+1),kgs(mgs))
      if (jyend.eq.nyend-1) ssfjp1(mgs) = ssfilt(igs(mgs),Min(ny,jgs+1),kgs(mgs))
      ssfjm1(mgs) = ssfilt(igs(mgs),jgs-1,kgs(mgs))
      if (jybeg.eq.nybeg) ssfjm1(mgs) = ssfilt(igs(mgs),Max(1,jgs-1),kgs(mgs))
      
      ssfip1(mgs) = ssfilt(igs(mgs)+1,jgs,kgs(mgs))
      if (ixend.eq.nxend) ssfip1(mgs) = ssfilt(Min(nx-1,igs(mgs)+1),jgs,kgs(mgs))
      if (ixend.eq.nxend-1) ssfip1(mgs) = ssfilt(Min(nx,igs(mgs)+1),jgs,kgs(mgs))
      ssfim1(mgs) = ssfilt(igs(mgs)-1,jgs,kgs(mgs))
      if (ixbeg.eq.nxbeg) ssfim1(mgs) = ssfilt(Max(1,igs(mgs)-1),jgs,kgs(mgs))
#else
      ssat0(mgs)  = ssat(igs(mgs),jgs,kgs(mgs))
      ssf(mgs)    = ssfilt(igs(mgs),jgs,kgs(mgs))
      ssfkp1(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs)+1)
      ssfkm1(mgs) = ssfilt(igs(mgs),jgs,Max(1,kgs(mgs)-1))
      ssfjp1(mgs) = ssfilt(igs(mgs),Min(ny-1,jgs+1),kgs(mgs))
      ssfjm1(mgs) = ssfilt(igs(mgs),Max(1,jgs-1),kgs(mgs))
      ssfip1(mgs) = ssfilt(Min(nx-1,igs(mgs)+1),jgs,kgs(mgs))
      ssfim1(mgs) = ssfilt(Max(1,igs(mgs)-1),jgs,kgs(mgs))
#endif

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
      if (my_rank==0) then
       if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: after ss gradient - - - - - - - - - "
       if (mgs ==1) write(0,*) "my_rank=",my_rank
       if (mgs ==1) then
        write(0,"('mgs',14x,'ssat',9x,'ssfilt',8x,'ssf',8x,'ssfkp1',7x,'ssfkm1',7x,'ssfjp1',7x,'ssfjm1',7x,'ssfip1',7x,'ssfim1')")
       endif
        write(0,"(4(1x,i3),9(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs), &
     &                               ssat(igs(mgs),jgs,kgs(mgs)),ssfilt(igs(mgs),jgs,kgs(mgs)),ssf(mgs),   &
     &                               ssfkp1(mgs),ssfkm1(mgs),ssfjp1(mgs),ssfjm1(mgs),ssfip1(mgs),ssfim1(mgs)
        end if
#else
       if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: after ss gradient - - - - - - - - - "
       if (mgs ==1) then
        write(0,"('mgs',14x,'ssat',9x,'ssfilt',8x,'ssf',8x,'ssfkp1',7x,'ssfkm1',7x,'ssfjp1',7x,'ssfjm1',7x,'ssfip1',7x,'ssfim1')")
       endif
        write(0,"(4(1x,i3),9(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs),ssat(igs(mgs),jgs,kgs(mgs)), &
     &                               ssfilt(igs(mgs),jgs,kgs(mgs)),ssf(mgs),   &
     &                               ssfkp1(mgs),ssfkm1(mgs),ssfjp1(mgs),ssfjm1(mgs),ssfip1(mgs),ssfim1(mgs)
#endif
!end mpidebug
      ENDIF

#ifdef MPI
      IF ( nyend .eq. 2 .and. ncdebug .ge. 1 ) THEN
#else
      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
#endif
       write(iunit,*) 'point1 ',igs(mgs)+ixbeg-1,kgs(mgs),ssat0(mgs),ssbar(mgs), &
     &  qx(mgs,lv)*1.e3,qss(mgs)*1.e3, qx(mgs,lc)*1.e3, &
     &  cx(mgs,lc)*1.e-6
      ENDIF

      ENDDO

!  Set density
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: Set density'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Set density'
#endif

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
      end do

!
!  cloud water variables
!
#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: Set cloud water variables'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Set cloud water variables'
#endif

      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      IF ( ipconc .ge. 2 .and. cx(mgs,lc) .gt. 1.0e6 ) THEN
        xmas(mgs,lc) = &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
      ELSE
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 0.01 ) THEN
        xmas(mgs,lc) = &
     &     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)), &
     &      xdn(mgs,lc)*xvmx(lc) )
        
        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)
        
       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 0.01 ) THEN
        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        
       ELSE
        xmas(mgs,lc) = cwmasn
       ENDIF
      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3

!        rb(mgs) = 0.5*xdia(mgs,lc,1)*(Sqrt(1./(1.+cnu)))**(1./3.)
!        xl2p(mgs) = 2.7e-2*xdn(mgs,lc)*cx(mgs,lc)*xv(mgs,lc)*
!     :           ((0.5e20*rb(mgs)**3*xdia(mgs,lc,1))-0.4)
!        IF ( rb(mgs) .gt. 3.51e-6 ) THEN
!          rh(mgs) = Max( 0.5*xdia(mgs,lc,1), 6.3e-4/(1.e6*(rb(mgs) - 3.5e-6)) )
!        ELSE
!          rh(mgs) = 1.e9
!        ENDIF
!        nh(mgs) = 4.2e9*xl2p(mgs)

      end do

#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: Set rain variables'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: Set rain variables'
#endif
!
! rain
!
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      
      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-9,cx(mgs,lr)))
!      parameter( xvmn(lr)=2.8866e-13, xvmx(lr)=4.1887e-9 )  ! mks
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
          xv(mgs,lr) = xvmx(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF

        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,1) = (xmas(mgs,lr)*cwc1)**(1./3.)
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

! Inverse exponential version:
!        xdia(mgs,lr,1) =
!     >  (qx(mgs,lr)*rho0(mgs)
!     > /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      ELSE
        xdia(mgs,lr,1) = &
     &  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if

      end do
!
!  Ziegler nucleation 
!
!
! cloud evaporation, condensation, and nucleation
!  sqsat -> qss(mgs)

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        write(0,*) " - - - - - - - - - WARMZIEG: begin zieg nucl - - - - - - - - - "
        write(0,*) "my_rank=",my_rank
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
        do mgs = 1,ngscnt
        write(0,"(i3,8(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr), &
     &                                    ssf(mgs),qx(mgs,lc),qxmin(lc),cx(mgs,lc)
        end do
        end if
#else
        write(0,*) " - - - - - - - - - WARMZIEG: begin zieg nucl - - - - - - - - - "
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
        do mgs = 1,ngscnt
        write(0,"(i3,8(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr), &
     &                                    ssf(mgs),qx(mgs,lc),qxmin(lc),cx(mgs,lc)
        end do
#endif
!end mpidebug
      ENDIF

#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: nucleation'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: nucleation, ngscnt = ',ngscnt
#endif

      DO mgs=1,ngscnt
      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: begin cloud cond - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
        write(0,*) " - - - - - - - - - WARMZIEG: start of nuc/cond loop - - - - - - - - - "
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'tem',11x,'qx',10x,'qxmin',11x,'cx')")
         write(0,"(i3,8(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr), &
     &                                    temg(mgs),qx(mgs,lc),qxmin(lc),cx(mgs,lc)
#endif
      ENDIF


        dcloud = 0.0
        IF ( temg(mgs) .le. tfrh ) THEN
        
        
         CYCLE
        ENDIF

        
      IF( ssat0(mgs) .GT. 0. .OR. ssf(mgs) .GT. 0. ) GO TO 620
!6/4      IF( qvap(mgs) .EQ. qss(mgs) ) GO TO 631
!
!.... EVAPORATION. QV IS LESS THAN qss(mgs).
!.... EVAPORATE CLOUD FIRST
!
      IF ( qx(mgs,lc) .LE. 0. ) GO TO 631
!.... CLOUD EVAPORATION.
! convert input 'cp' to cgs
      R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/ &
     &            (cp*(temg(mgs) - cbw)**2))
      QEVAP= Min( qx(mgs,lc), R1*(qss(mgs)-qvap(mgs)) )
      
      
      IF ( qx(mgs,lc) .LT. QEVAP ) THEN ! GO TO 63
        qwvp(mgs) = qwvp(mgs) + qx(mgs,lc)
        thetap(mgs) = thetap(mgs) - felv(mgs)*qx(mgs,lc)/(cp*pi0(mgs))
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felv(mgs)*qx(mgs,lc)/(cp*pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - qx(mgs,lc)*rho0(mgs)*dv/dtp  ! evaporation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) - qx(mgs,lc)/dtp*felv(mgs)/(cp*pi0(mgs))*3600.0
        ENDIF
        qx(mgs,lc) = 0.
        cx(mgs,lc) = 0.
      ELSE
        qwvp(mgs) = qwvp(mgs) + QEVAP
        qx(mgs,lc) = qx(mgs,lc) - QEVAP
        IF ( qx(mgs,lc) .le. 0. ) cx(mgs,lc) = 0.
        thetap(mgs) = thetap(mgs) - felv(mgs)*QEVAP/(CP*pi0(mgs))
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felv(mgs)*QEVAP/(CP*pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - QEVAP*rho0(mgs)*dv/dtp  ! evaporation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) - QEVAP/dtp*felv(mgs)/(cp*pi0(mgs))*3600.0
        ENDIF
      ENDIF

      GO TO 631


  620 CONTINUE

!.... CLOUD CONDENSATION

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: begin cloud cond - - - - - - - - - "
        if(mgs==1) write(0,*) "my_rank=",my_rank,igs(mgs)+ixbeg-1,jgs+jybeg-1,kgs(mgs)
        if(mgs==1) write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
         write(0,"(i3,8(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr), &
     &                                    ssf(mgs),qx(mgs,lc),qxmin(lc),cx(mgs,lc)
        end if
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: begin cloud cond - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
        write(0,*) " - - - - - - - - - WARMZIEG: begin cloud cond - - - - - - - - - "
        write(0,"('mgs',7x,'qc',12x,'qr',11x,'ccw',11x,'crw',11x,'ssf',11x,'qx',10x,'qxmin',11x,'cx')")
         write(0,"(i3,8(1x,es12.5))") mgs,an(igs(mgs),jy,kgs(mgs),lc), an(igs(mgs),jy,kgs(mgs),lr), &
     &                                    an(igs(mgs),jy,kgs(mgs),lnc),an(igs(mgs),jy,kgs(mgs),lnr), &
     &                                    ssf(mgs),qx(mgs,lc),qxmin(lc),cx(mgs,lc)
#endif
!end mpidebug
      ENDIF


        IF ( qx(mgs,lc) .GT. qxmin(lc) .and. cx(mgs,lc) .ge. 1. ) THEN


!       ac1 =  xdn(mgs,lc)*elv(kgs(mgs))**2*epsi/
!     :        (tka(kgs(mgs))*rw*temg(mgs)**2)
! took out xdn factor because it cancels later...
       ac1 =  felv(mgs)**2*epsi/(tka(kgs(mgs))*rw*temg(mgs)**2)
       

!       bc = xdn(mgs,lc)*rw*temg(mgs)/
!     :       (epsi*wvdf(kgs(mgs))*es(mgs))
! took out xdn factor because it cancels later...
       bc =   rw*temg(mgs)/(epsi*wvdf(kgs(mgs))*es(mgs))


!      
      IF ( ssf(mgs) .gt. 0.0 .or. ssat0(mgs) .gt. 0.0 ) THEN
       IF ( ny .le. 2 ) THEN
!        print*, 'undershoot: ',ssf(mgs),
!     :   ( (qx(mgs,lv) - dcloud)/c1 - 1.0)*100.
       ENDIF


       
       IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN

         IF ( xdia(mgs,lc,1) .le. 0.0 ) THEN
          xmas(mgs,lc) = cwmasn
          xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3 
         ENDIF
        d1 = (1./(ac1 + bc))*4.0*pi*ventc &
     &        *0.5*xdia(mgs,lc,1)*cx(mgs,lc)*rhoinv(mgs)
       ELSE
         d1 = 0.0
       ENDIF

       IF ( rcond .eq. 2 .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
       rwvent(mgs) = ventrx(mgs)*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)

       d1r = (1./(ac1 + bc))*4.0*pi*rwvent(mgs) &
     &        *0.5*xdia(mgs,lr,1)*cx(mgs,lr)*rhoinv(mgs)
       ELSE
       d1r = 0.0
       ENDIF


       
       e1  = felv(mgs)/(CP*pi0(mgs))
       f1 = pk(mgs) ! (pres(mgs)/poo)**cap

!
!  fifth trial to see what happens: 
!
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       ltemq1 = ltemq 
       temp1 = temg(mgs)
       p380 = 380.0/pres(mgs)
       
!       nc = NInt(dtp/Min(1.0d0,0.5*taus))
!       dtcon = dtp/float(nc)
       ss1 = qx(mgs,lv)/qvs(mgs)
       ss2 = ss1
       ssi1 = qx(mgs,lv)/qis(mgs)
       ssi2 = ssi1
       temp2 = temp1
       qv1 = qx(mgs,lv)
       qvs1 = qvs(mgs)
       qis1 = qis(mgs)
       dt1 = 0.0d0

          
!          dtcon = Max(dtcon,0.2d0)
!          nc = Nint(dtp/dtcon)

       ltemq1 = ltemq 
! want to start out with a small time step to handle the steep slope
! and fast changes, then can switch to a larger step (dtcon2) for the
! rest of the big time step.
! base the initial time step (dtcon1) on the slope (delta)
       IF ( Abs(ss1 - 1.0) .gt. 1.e-5 ) THEN
         delta = 0.5*(qv1-qvs1)/(d1*(ss1 - 1.0))
       ELSE
         delta = 0.1*dtp
       ENDIF
! delta is the extrapolated time to get halfway from qv1 to qvs1
! want at least 5 time steps to the halfway point, so multiply by 0.2
! for the initial time step
       dtcon1 = Min(0.05,0.2*delta)
       nc = Max(5,2*NInt( (dtp-4.0*dtcon1)/delta))
       dtcon2 = (dtp-4.0*dtcon1)/nc
!       print*,'delt,taus = ',delta,taus,
!     :   delta/taus,dtcon1,dtcon2
!          print*, 'RK2c start',nc,dtcon,qv1*1.e3,qvs1*1.e3,
!     :       ss1*100.,temp1

!
! 2006.5.20: Test adding rain condensation
!
       n = 1
       dt1 = 0.0d0
       nc = 0
       dqc = 0.0d0
       dqr = 0.0d0
       dqvii = 0
       dqvis = 0

      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
       write(iunit,*) 'point3a ',igs(mgs)+ixbeg-1,kgs(mgs)
       write(iunit,*) 'delta = ',delta,dtcon1,nc,dtcon2,d1r
       write(iunit,*) 'ss1,ssi1,temp2,qv1,qvs1 = ',ss1,ssi1,temp2,qv1,qvs1
      ENDIF
       
       RK2c: DO WHILE ( dt1 .lt. dtp ) 
          nc = 0
          IF ( n .le. 4 ) THEN
            dtcon = dtcon1
          ELSE
            dtcon = dtcon2
          ENDIF
 609       dqv  = -(ss1 - 1.)*d1*dtcon
           dqvr = -(ss1 - 1.)*d1r*dtcon
            dtemp = -0.5*e1*f1*(dqv + dqvr)
!          print*,'RK2c dqv1 = ',dqv
! calculate midpoint values:
           ltemq1m = ltemq1 + Nint(dtemp*fqsat + 0.5)
            dqvs = dtemp*p380*dtabqvs(ltemq1m)
            qv1m = qv1 + dqv + dqvr
!          qv1mr = qv1r + dqvr

            qvs1m = qvs1 + dqvs
            ss1m = qv1m/qvs1m

    ! check for undersaturation when no ice is present, if so, then reduce time step
          IF ( ss1m .lt. 1.  .and. (dqvii + dqvis) .eq. 0.0 ) THEN
!            dtcon = Max(dtcon1,0.5*dtcon)
!            IF ( dtcon .gt. dtcon1 ) GOTO 609
            dtcon = (0.5*dtcon)
            IF ( dtcon .ge. dtcon1 ) THEN
             GOTO 609
            ELSE
 !            print*,'dtcon = ',dtcon,n,dt1,temg(mgs),qx(mgs,lc),qx(mgs,lr),qx(mgs,li),qx(mgs,ls)
 !            print*,'dtcon2 : ',dqc,dqr,dqi,dqs,mgs,igs(mgs),kgs(mgs),ssf(mgs)
             EXIT
            ENDIF
          ENDIF
! calculate full step:
          dqv  = -(ss1m - 1.)*d1*dtcon
          dqvr = -(ss1m - 1.)*d1r*dtcon


!          print*,'RK2a dqv1m = ',dqv
          dtemp = -e1*f1*(dqv + dqvr)
          ltemq1 = ltemq1 + Nint(dtemp*fqsat + 0.5)
          dqvs = dtemp*p380*dtabqvs(ltemq1)

          qv1 = qv1 + dqv + dqvr

          dqc = dqc - dqv
          dqr = dqr - dqvr

          qvs1 = qvs1 + dqvs
          ss1 = qv1/qvs1
          temp1 = temp1 + dtemp
!          write(iunit,*) 'RK2c ',n,qv1*1000.,qvs1*1000.,ss1*100.,temp1,
!     :     dt1+dtcon,dtcon, dqvs, dqv, ss1m
!          write(iunit,*) 'RK2c: dtemp ',dtemp,dqvs,qv1m,qvs1m
!          write(iunit,*) 'RK2c: ltemq1m', ltemq1m, dqvr, dqr
          IF ( temp2 .eq. temp1 .or. ss2 .eq. ss1 .or.  &
     &           ss1 .eq. 1.00 .or.  &
     &      ( n .gt. 10 .and. ss1 .lt. 1.0005 ) ) THEN
!           print*,'RK2c break'
           EXIT
          ELSE
           ss2 = ss1
           ssi2 = ssi1
           temp2 = temp1
           dt1 = dt1 + dtcon
!           IF ( n .ge. 5 ) dtcon = 1.25*dtcon
!           dtcon = Min(dtcon, dtp - dt1)
           n = n + 1
          ENDIF
       ENDDO RK2c
       
        
        dcloud = dqc ! qx(mgs,lv) - qv1
        thetap(mgs) = thetap(mgs) + e1*(DCLOUD + dqr)
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + e1*(DCLOUD + dqr)*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + (DCLOUD + dqr)*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + &
     &     (DCLOUD + dqr)/dtp*felv(mgs)/(cp*pi0(mgs))*3600.0
        ENDIF
        qwvp(mgs) = qwvp(mgs) - (DCLOUD + dqr)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        qx(mgs,lr) = qx(mgs,lr) + dqr

        IF ( lzr > 1 .and. rcond == 2 .and. qx(mgs,lr) .gt. qxmin(lr)   &
     &       .and. cx(mgs,lr) .gt. 1.e-9 ) THEN
          tmp = qx(mgs,lr)/cx(mgs,lr)
          g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
          zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(xdn(mgs,lr)))**2*( 2.*( tmp ) * dqr )
        ENDIF


        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*f1
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        es(mgs) = 6.1078e2*tabqvs(ltemq)
        

      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
!      IF ( ny .eq. 2 .and. igs(mgs)+ixbeg-1 .lt. 25 ) THEN
       write(iunit,*) 'point3b ',igs(mgs)+ixbeg-1,kgs(mgs),dcloud*1.e3 ! ,dcloud2*1.e3
       write(iunit,*) 'ac,bc = ',ac1,bc,d1,e1,f1
       write(iunit,*) 'cwdia,ccw,cwdn,c1 = ',xdia(mgs,lc,1),cx(mgs,lc)*1.e-6, &
     &     xdn(mgs,lc),qx(mgs,lc)*1.e3
      ENDIF
!            
      
      ENDIF  ! dcloud .gt. 0.
     
!       write(iunit,*) 'cond:ix,kz,dcloud,qcw,cwdia,ccw,taus = ',
!     :   igs(mgs),kgs(mgs),dcloud,qx(mgs,lc),
!     :   xdia(mgs,lc,1),cx(mgs,lc),taus

      ELSE  ! qc .le. qxmin(lc)

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 1 - - - - - - - - - "
        if(mgs==1) write(0,*) "my_rank=",my_rank
        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
        end if
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 1 - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 1 - - - - - - - - - "
         write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
#endif
!end mpidebug
      ENDIF
      
        IF ( ssf(mgs) .gt. 0.0 ) THEN ! .and.  ssmax(mgs) .lt. sscb ) THEN

          IF ( iqcinit == 1 ) THEN

!         temp  = pb(k) * t(i,j,k)
         
         qvs0   = 380.*exp(17.27*(temg(mgs)-273.)/(temg(mgs)- 36.))/pk(mgs)

         dcloud = Max(0.0, (qx(mgs,lv)-qvs0) / (1.+qvs0*f5/(temg(mgs)-36.)**2) )
          
          
          ELSEIF ( iqcinit == 3 ) THEN
              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/ &
     &             (cp*(temg(mgs) - cbw)**2))
            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
          
          ELSEIF ( iqcinit == 2 ) THEN
!              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/
!     :             (cp*(temg(mgs) - cbw)**2))
!            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
         ssmx = ssmxinit

         CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,dcloud, & 
     &    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felv,ssmx,pk,ngscnt)

         ENDIF
         
         IF ( ncdebug >= 1 .and. ny == 2) THEN
           write(iunit,*) 'init qc: i,k,dcloud = ',igs(mgs),kgs(mgs),dcloud*1.e3,iqcinit
         ENDIF
         
        ELSE
            dcloud = 0.0
        ENDIF

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 1 - - - - - - - - - "
        if(mgs==1) write(0,*) "my_rank=",my_rank
        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
        end if
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 1 - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 1 - - - - - - - - - "
         write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
#endif
!end mpidebug
      ENDIF

!      PT=PT+SLV*DCLOUD/(CP*PIB(K))
        thetap(mgs) = thetap(mgs) + felv(mgs)*DCLOUD/(CP*pi0(mgs))
        IF ( lsati > 1 ) THEN
!          axtra(igs(mgs),jy,kgs(mgs),lsati) = axtra(igs(mgs),jy,kgs(mgs),lsati) + (DCLOUD)/dtp*felv(mgs)/(cp*pi0(mgs)) 
        ENDIF
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felv(mgs)*DCLOUD/(CP*pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felv(mgs)/(cp*pi0(mgs))*3600.0
        ENDIF

        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) !( pres(mgs) / poo ) ** cap
!        temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        es(mgs) = 6.1078e2*tabqvs(ltemq)

        END IF ! qc .gt. 0.

!        ES=EES(PIB(K)*PT)
!        SQSAT=EPSI*ES/(PB(K)*1000.-ES)

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
      if (my_rank==0) then
       if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: before cloud nucl - - - - - - - - - "
       if (mgs ==1) write(0,*) "my_rank=",my_rank
       if (mgs ==1) write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,"(4(1x,i3),6(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs), &
     &                               dssdx,dssdy,dssdz,uvel(mgs),vvel(mgs),wvel(mgs)
        end if
#else
!       if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: before cloud nucl - - - - - - - - - "
!       if (mgs ==1) write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,*) " - - - - - - - - - WARMZIEG: before cloud nucl - - - - - - - - - "
        write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,"(4(1x,i3),6(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs), &
     &                               dssdx,dssdy,dssdz,uvel(mgs),vvel(mgs),wvel(mgs)
#endif
!end mpidebug
      ENDIF

!mpidebug: need to code cloud nucl for MPI if vertical tiling

!.... CLOUD NUCLEATION
!      T=PIB(K)*PT
!      ES=1.E3*PB(K)*QV/EPSI

      IF ( wvel(mgs) .le. 0. ) GO TO 616
      IF ( cx(mgs,lc) .le. 0. )  GO TO 613                             !TWOMEY (1959) Nucleation
      IF ( kzbeg-1+kgs(mgs) .GT. 1 .and. qx(mgs,lc) .le. qxmin(lc)) GO TO 613  !TWOMEY (1959) Nucleation
      IF ( kzbeg-1+kgs(mgs) .eq. 1 .and. wvel(mgs) .gt. 0. ) GO TO 613         !TWOMEY (1959) Nucleation

!.... ATTEMPT ZIEGLER CLOUD NUCLEATION IN CLOUD INTERIOR UNLESS...
  616 IF ( ssf(mgs) .LE. SUPCB .AND. wvel(mgs) .GT. 0. ) GO TO 631 !... weakly saturated updraft
      IF ( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 .AND.  &
     &    (ssfkp1(mgs) .GE. SUPMX .OR. &
     &     ssf(mgs)    .GE. SUPMX .OR. &
     &     ssfkm1(mgs) .GE. SUPMX)) GO TO 631                      !... too much vapour
      IF (ssf(mgs) .LT. 1.E-10 .OR. ssf(mgs) .GE. SUPMX) GO TO 631 !... at the extremes for ss

!
! get here if ( qc > 0 and ss > supcb) or (w < 0)
!

#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: Entered Ziegler Cloud Nucleation" !mpidebug
#else
      if (ndebug .gt. 0) write(0,*) "ICEZVD_DR: Entered Ziegler Cloud Nucleation" !mpidebug
#endif

      DSSDZ=0.
      DSSDY=0.
      DSSDX=0.
      r2dxm=1./(2.*dx)
      r2dym=1./(2.*dy)
      r2dzm=1./(2.)*gz(kgs(mgs))
      IF ( irenuc >= 0 ) THEN

        IF ( kzend == nzend ) THEN
          t0p3 = t0(igs(mgs),jgs,Min(kze,kgs(mgs)+3))
          t0p1 = t0(igs(mgs),jgs,Min(kze,kgs(mgs)+1))
        ELSE
          t0p3 = t0(igs(mgs),jgs,kgs(mgs)+3)
          t0p1 = t0(igs(mgs),jgs,kgs(mgs)+1)
        ENDIF

      IF ( ( ssf(mgs) .gt. ssmax(mgs) .or.  irenuc .eq. 1 ) & 
     &   .and.  ( ( lccn .lt. 1 .and.  & 
     &            cx(mgs,lc) .lt. cwccn*(Min(1.0,rho0(mgs)))) .or. & 
     &    ( lccn .gt. 1 .and. ccnc(mgs) .gt. 0. )   ) & 
     &    ) THEN
      IF( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 & 
     &  .and. ssf(mgs) .gt. 0.0 & 
     &  .and. ssfkp1(mgs) .LT. SUPMX .and. ssfkp1(mgs) .ge. 0.0  &
     &  .AND. ssfkm1(mgs) .LT. SUPMX .AND. ssfkm1(mgs) .ge. 0.0  & 
     &  .AND. ssfkp1(mgs) .gt. ssfkm1(mgs)  & 
     &  .and. t0p3 .gt. 233.2) THEN
          DSSDZ = (ssfkp1(mgs) - ssfkm1(mgs))*R2DZM
!
! otherwise check for cloud base condition with updraft:
!
        ELSEIF( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 &
!        IF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1 & 
     &  .and. ssf(mgs) .gt. 0.0  .and. wvel(mgs) .gt. 0.0 & 
     &  .and. ssfkp1(mgs) .gt. 0.0   &
     &  .AND. ssfkm1(mgs) .le. 0.0 .and. wvelkm1(mgs) .gt. 0.0 & 
     &  .AND. ssf(mgs) .gt. ssfkm1(mgs)  & 
     &  .and. t0p1 .gt. 233.2) THEN
         DSSDZ = 2.*(ssf(mgs) - ssfkm1(mgs))*R2DZM  ! 1-sided difference
        ENDIF

      IF ( irenuc3d > 0 ) THEN !{
 
      IF( jybeg-1+jgs .GT. 1 .AND. jybeg-1+jgs .LT. nyend-1 &
     &    .and. ssf(mgs) .gt. 0.0 & 
     &    .and. ssfjm1(mgs) .LT. SUPMX .and. ssfjm1(mgs) .ge. supcb   &
     &    .AND. ssfjp1(mgs) .LT. SUPMX .and. ssfjp1(mgs) .ge. supcb ) &
!     :    .AND. ssfjp1(mgs) .GT. ssfjm1(mgs)
     &    DSSDY = (ssfjp1(mgs) - ssfjm1(mgs))*R2DYM
      
      IF( ixbeg-1+igs(mgs) .GT. 1 .AND. ixbeg-1+igs(mgs) .LT. nzend-1 &
     &  .and. ssf(mgs) .gt. 0.0 & 
     &  .and. ssfim1(mgs) .LT. SUPMX .and. ssfim1(mgs) .ge. supcb  &
     &  .AND. ssfip1(mgs) .LT. SUPMX .AND. ssfip1(mgs) .ge. supcb ) &
     &    DSSDX = (ssfip1(mgs) - ssfim1(mgs))*R2DXM
!     :  .AND. ssfip1(mgs) .GT. ssfim1(mgs)
     
       ENDIF !}
       ENDIF
!
!CLZ  IF(wijk.LE.0.) CN=CCN*ssfilt(ix,jy,kz)**CCK
! note: CCN -> cwccn, DELT -> dtp
      c1 = Max(0.0, rho0(mgs)*(qx(mgs,lv) - qss(mgs))/ & 
     &        (xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3))
      IF ( lccn .lt. 1 ) THEN ! not predicting CCN
       CN(mgs) = cwccn*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         ( uvel(mgs)*DSSDX +    &  ! these 3 lines are V dot Grad(SS)
     &           vvel(mgs)*DSSDY +    &  ! where SS = supersaturation
     &           wvel(mgs)*DSSDZ) )      ! probably the vertical gradient dominates
      ELSE ! CCN conc. is predicted
      CN(mgs) =  &
!     :   Min(Min(c1,ccnc(mgs)), cwccn*CCK*ssf(mgs)**CCKM*dtp*
!     :   Min(ccnc(mgs), cwccn*CCK*ssf(mgs)**CCKM*dtp*
!     &    ( cwccn*CCK*ssf(mgs)**CCKM*Min(1.0,dtp)*   &
     &    Min(ccnc(mgs), cnuc(mgs)*CCK*ssf(mgs)**CCKM*dtp*   &
!     &    ( cwccn*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         ( uvel(mgs)*DSSDX +    &
     &           vvel(mgs)*DSSDY +    &
     &           wvel(mgs)*DSSDZ) )  )
!      IF ( cn(mgs) .gt. 0 ) ccnc(mgs) = ccnc(mgs) - cn(mgs)
      ENDIF
      
      IF ( cn(mgs) .gt. 0.0 ) THEN
       IF ( ccnc(mgs) .lt. 5.e7 .and. cn(mgs) .ge. 5.e7 ) THEN
          cn(mgs) = 5.e7
          ccnc(mgs) = 0.0
       ELSEIF ( cn(mgs) .gt. ccnc(mgs) ) THEN
         cn(mgs) = ccnc(mgs)
         ccnc(mgs) = 0.0
       ENDIF
      cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      ENDIF

      ENDIF ! irenuc >= 0

      IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .LE. qxmin(lc)) cx(mgs,lc)=0.
      GO TO 631
!.... NUCLEATION ON CLOUD INFLOW BOUNDARY POINT

  613 CONTINUE
!.... S. TWOMEY (1959)
! Note: get here if there is no previous cloud water and w > 0.
      cn(mgs) = 0.0
      
      IF ( ncdebug .ge. 1 ) THEN
        write(iunit,*) 'at 613: ',qx(mgs,lc),cx(mgs,lc),wvel(mgs),ssmax(mgs),kgs(mgs)
      ENDIF
      
!      IF ( ssmax(mgs) .lt. sscb .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
      IF ( qx(mgs,lc) .gt. qxmin(lc)) THEN
       CN(mgs) =   CCNE*wvel(mgs)**cnexp ! *Min(1.0,1./dtp) ! 0.3465
        IF ( ny .le. 2 .and. cn(mgs) .gt. 0.0    &
     &                    .and. ncdebug .ge. 1 ) THEN ! .and. kgs(mgs) <= 6 ) THEN
          write(iunit,*) 'CN: ',cn(mgs)*1.e-6, cx(mgs,lc)*1.e-6, qx(mgs,lc)*1.e3,   &
     &       wvel(mgs), dcloud*1.e3
          IF ( cn(mgs) .gt. 1.0 ) write(iunit,*) 'cwrad = ',   &
     &       1.e6*(rho0(mgs)*qx(mgs,lc)/cn(mgs)*cwc1)**c1f3,   &
     &   igs(mgs),kgs(mgs),temcg(mgs),    &
     &   1.e3*an(igs(mgs),jgs,kgs(mgs)-1,lc)
        ENDIF
        IF ( iccwflg .eq. 1 ) THEN
          cn(mgs) = Min(cwccn, Max(cn(mgs),   &
     &       rho0(mgs)*qx(mgs,lc)/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3)))
        ENDIF
      ELSE
       cn(mgs) = 0.0
      ENDIF

      IF( CN(mgs) .GT. cx(mgs,lc) ) cx(mgs,lc) = CN(mgs)
      IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .le. qxmin(lc) ) THEN
        cx(mgs,lc) = 0.
      ELSE
        cx(mgs,lc) = Min(cx(mgs,lc),rho0(mgs)*Max(0.0,qx(mgs,lc))/cwmasn)
      ENDIF

  631  CONTINUE

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
      if (my_rank==0) then
       if (mgs ==1) write(0,*) " - - - - - - - - -   - - - - - - - - - "
       if (mgs ==1) write(0,*) "my_rank=",my_rank
       if (mgs ==1) write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,"(4(1x,i3),6(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs), &
     &                               dssdx,dssdy,dssdz,uvel(mgs),vvel(mgs),wvel(mgs)
        end if
#else
!       if (mgs ==1) write(0,*) " - - - - - - - - - WARMZIEG: after cloud nucl - - - - - - - - - "
!       if (mgs ==1) write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,*) " - - - - - - - - - WARMZIEG: after cloud nucl - - - - - - - - - "
        write(0,"('mgs',14x,'dssdx',9x,'dssdy',8x,'dssdz',8x,'uvel',7x,'vvel',7x,'wvel')")
        write(0,"(4(1x,i3),6(1x,es12.5))") mgs,igs(mgs),jgs,kgs(mgs), &
     &                               dssdx,dssdy,dssdz,uvel(mgs),vvel(mgs),wvel(mgs)
#endif
!end mpidebug
      ENDIF

!
! Check for supersaturation greater than ssmx and adjust down
!
       ssmx = 1.1
       qv1 = qv0(mgs) + qwvp(mgs)
       qvs1 = qvs(mgs)

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 2 - - - - - - - - - "
        if(mgs==1) write(0,*) "my_rank=",my_rank
        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
        end if
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 2 - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,*) " - - - - - - - - - WARMZIEG: before qvexcess 2 - - - - - - - - - "
         write(0,"('mgs',5x,'qv1',10x,'qvs1',10x,'ss1',9x,'pres',8x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
#endif
!end mpidebug
      ENDIF

       IF ( qv1 .gt. (ssmx*qvs1) ) THEN
        
         ss1 = qv1/qvs1

        ssmx = 100.*(ssmx - 1.0)
        
        qvex = 0.0

        CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,qvex,   &
     &    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felv,ssmx,pk,ngscnt)


       IF ( ny .le. 2 .and. ncdebug .ge. 1 ) THEN
        write(iunit,'(a,3(f8.3,1x),2i4)')    &
     &   'Large SS1: ',(ss1-1.0)*100.,   &
     &   qx(mgs,lv)*1.e3,qvs(mgs)*1.e3,igs(mgs),kgs(mgs)
        write(iunit,*) 'qc,qwvp,qvex',qx(mgs,lc)*1.e3, qwvp(mgs)*1.e3,qvex*1.e3
        write(iunit,*) 'qv1,dcloud ',qv1*1.e3,dcloud*1.e3
        
       ENDIF

        IF ( qvex .gt. 0.0 ) THEN
        thetap(mgs) = thetap(mgs) + felv(mgs)*qvex/(CP*pi0(mgs))
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felv(mgs)*qvex/(CP*pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + qvex*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + qvex/dtp*felv(mgs)/(cp*pi0(mgs))*3600.0
        ENDIF
        qwvp(mgs) = qwvp(mgs) - qvex
        qx(mgs,lc) = qx(mgs,lc) + qvex
        cn(mgs) = Min( cwnccn(kgs(mgs)), qvex/Max( cwmasn5, xmas(mgs,lc) )  )
        ccnc(mgs) = Max( 0.0, ccnc(mgs) - cn(mgs) )
        cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
        
!        write(iunit,*) 'theta = ',theta0(mgs) + thetap(mgs)

!        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap

        ENDIF

       IF ( ny .le. 2 .and. ncdebug .ge. 1 ) THEN
        write(iunit,*) 'theta = ',theta0(mgs) + thetap(mgs)
       ENDIF
       
       ENDIF

      IF ( ndebug .gt. 2 ) THEN
!mpidebug
#ifdef MPI
        if (my_rank==0) then
        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 2 - - - - - - - - - "
        if(mgs==1) write(0,*) "my_rank=",my_rank
        if(mgs==1) write(0,"('mgs',7x,'qv1',11x,'qvs1',10x,'ss1',10x,'pres',9x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
        end if
#else
!        if(mgs==1) write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 2 - - - - - - - - - "
!        if(mgs==1) write(0,"('mgs',7x,'qv1',11x,'qvs1',10x,'ss1',10x,'pres',9x,'tabqvs',10x,'ssf')")
         write(0,*) " - - - - - - - - - WARMZIEG: after qvexcess 2 - - - - - - - - - "
         write(0,"('mgs',7x,'qv1',11x,'qvs1',10x,'ss1',10x,'pres',9x,'tabqvs',10x,'ssf')")
         write(0,"(i3,6(1x,es12.5))") mgs,qv1,qvs1,qv1/qvs1, pres(mgs),tabqvs(mgs),ssf(mgs)
#endif
!end mpidebug
      ENDIF

!
! Calculate droplet volume and check if it is within bounds.
!  Adjust if necessary
!  


      cx(mgs,lc) = Min( cwnccn(kgs(mgs)), cx(mgs,lc) )
      IF( cx(mgs,lc) .GT. 1.0e7 .AND. qx(mgs,lc) .GT. qxmin(lc)) THEN
!        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc)*xdn(mgs,lc))
        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc))
      ENDIF

!      svc(mgs) = Min( svc(mgs), xvmx(lc) )
!      svc(mgs) = Max( svc(mgs), xvmn(lc) )
      xmas(mgs,lc) = Min( xmas(mgs,lc), cwmasx )
      xmas(mgs,lc) = Max( xmas(mgs,lc), cwmasn )
!      IF( SVC .GT. xvmx(lc) ) SVC=xvmx(lc)
!      IF( SVC .LT. xvmn(lc) ) SVC=xvmn(lc)


      IF( cx(mgs,lc) .GT. 10.e6 .AND. qx(mgs,lc) .GT. qxmin(lc) ) GO TO 681
        ccwtmp = cx(mgs,lc)
        cwmastmp = xmas(mgs,lc)
!       svc(mgs) = xvmn(lc)
       xmas(mgs,lc) = Max(xmas(mgs,lc), cwmasn)
!      IF(qx(mgs,lc) .GT. 0. .AND. cx(mgs,lc) .EQ. 0.) 
       IF(qx(mgs,lc) .GT. qxmin(lc) .AND. cx(mgs,lc) .le. 0.) THEN
          cx(mgs,lc) = Min(0.5*cwccn,rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc))
          xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
       ENDIF
!     :        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(SVC(mgs)*xdn(mgs,lc))
!      IF(cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. 0.) 
      IF(cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. qxmin(lc))    &
     &        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!     :        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!      IF(qx(mgs,lc) .GT. 0. .AND. SVC(mgs) .LT. xvmn(lc)) SVC(mgs) = xvmn(lc)
!      IF(qx(mgs,lc) .GT. 0. .AND. SVC(mgs) .GT. xvmx(lc)) SVC(mgs) = xvmx(lc)
!      IF(qx(mgs,lc) .GT. 0.) cx(mgs,lc) = 
!     :     rho0(mgs)*qx(mgs,lc)/(SVC(mgs)*xdn(mgs,lc))
!      IF(qx(mgs,lc) .GT. 0. .AND. xmas(mgs,lc) .LT. cwmasn) 
      IF(qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .LT. cwmasn)    &
     &          xmas(mgs,lc) = cwmasn
!      IF(qx(mgs,lc) .GT. 0. .AND. xmas(mgs,lc) .GT. cwmasx) 
      IF(qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .GT. cwmasx)    &
     &    xmas(mgs,lc) = cwmasx
!      IF(qx(mgs,lc) .GT. 0.) cx(mgs,lc) = 
      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/Max(cwmasn,xmas(mgs,lc))
!        write(iunit,*) 'ccw check3: i,k,ccw,qcw,cwmas,cwmasn = ',
!     :     igs(mgs),kgs(mgs),1.e-6*cx(mgs,lc),1.e3*qx(mgs,lc),
!     :     xmas(mgs,lc),cwmasn
      ENDIF
        
!      IF ( cx(mgs,lc) .gt. 1.5e9 ) THEN
!         write(iunit,*) 
!     : 'ccw check4: i,k,ccw,ccwtmp,qcw,cwmas,cwmastmp,cwmasn = ',
!     :     igs(mgs),kgs(mgs),1.e-6*cx(mgs,lc),1.e-6*ccwtmp,1.e3*qx(mgs,lc),
!     :     xmas(mgs,lc),cwmastmp,cwmasn
!      ENDIF

 681  CONTINUE
        
      IF ( ipconc .ge. 3 ) THEN

        
        IF(cx(mgs,lr) .GT. 0. .AND. qx(mgs,lr) .GT. qxmin(lr))    &
     &       xv(mgs,lr)=rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
        IF(xv(mgs,lr) .GT. xvmx(lr)) xv(mgs,lr) = xvmx(lr)
        IF(xv(mgs,lr) .LT. xvmn(lr)) xv(mgs,lr) = xvmn(lr)
!        IF ( .not. (cx(mgs,lr) .gt. 0. .and. qx(mgs,lr) .gt. qxmin(lr))) THEN
!          xv(mgs,lr) = xvmn(lr)
!          IF(qx(mgs,lr) .gt. 0. .and. cx(mgs,lr) .eq. 0.) 
!     :       cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*xv(mgs,lr))
!          IF(cx(mgs,lr) .gt. 0. .and. qx(mgs,lr) .gt. 0.) 
!     :       xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
!          IF( qx(mgs,lr) .gt. 0. .and. xv(mgs,lr) .lt. xvmn(lr)) xv(mgs,lr)=xvmn(lr)
!          IF( qx(mgs,lr) .gt. 0. .and. xv(mgs,lr) .gt. xvmx(lr)) xv(mgs,lr)=xvmx(lr)
!          IF( qx(mgs,lr) .gt. 0.) 
!     :        cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*xv(mgs,lr))
!        ENDIF
!      IF ( qx(mgs,lr) .gt. 10.0e-3 )
!     :  print*, 'RAIN2b: ',igs(mgs),kgs(mgs),qx(mgs,lr)

      ENDIF
      


      ENDDO ! mgs

#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: done with cond/nuc'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: done with cond/nuc'
#endif

      DO mgs=1,ngscnt
      IF ( ssf(mgs) .gt. ssmax(mgs)  &
     &  .and. ( idecss .eq. 0 .or. qx(mgs,lc) .gt. qxmin(lc)) ) THEN
        ssmax(mgs) = ssf(mgs)
      ENDIF
      ENDDO
!

      do mgs = 1,ngscnt

      an(igs(mgs),jy,kgs(mgs),lt) =  &
     &  ab(kgs(mgs),lt) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) =  &
     &  ab(kgs(mgs),lv) + qwvp(mgs) 
!
        an(igs(mgs),jy,kgs(mgs),lc) = qx(mgs,lc) +  &
     &    min( an(igs(mgs),jy,kgs(mgs),lc), 0.0 )  
!        qx(mgs,lc) = an(igs(mgs),jy,kgs(mgs),lc)
!
        an(igs(mgs),jy,kgs(mgs),lr) = qx(mgs,lr) +  &
     &    min( an(igs(mgs),jy,kgs(mgs),lr), 0.0 )  
!        qx(mgs,lr) = an(igs(mgs),jy,kgs(mgs),lr)

        IF ( lzr > 1 .and. rcond == 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lzr) = zx(mgs,lr) +  &
     &    min( an(igs(mgs),jy,kgs(mgs),lzr), 0.0 )  
        ENDIF
       IF (  ipconc .ge. 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnc) = Max(cx(mgs,lc) , 0.0)
        an(igs(mgs),jy,kgs(mgs),lss) = Max( 0.0, ssmax(mgs) )
        IF ( lccn .gt. 1 ) THEN
          an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0, Min( ccwmx, ccnc(mgs) ) )
        ENDIF
       ENDIF
       IF (  ipconc .ge. 3 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnr) = Max(cx(mgs,lr) , 0.0)
       ENDIF
      end do

#ifdef MPI
      if (ndebug .gt. 0 .and. my_rank .eq. 0) write(0,*) 'WARMZIEG: start deallocate 2'
#else
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: start deallocate 2'
#endif

      deallocate ( qx )
      deallocate ( cx )
      deallocate ( xv )
      deallocate ( vtxbar )
      deallocate ( xmas )
      deallocate ( xdn )
      deallocate ( xdia )
      deallocate ( alpha )
      deallocate ( zx )

29998 continue

#ifdef MPI
      if ( kzbeg-1+kz .gt. nzend-kstag-1 .and. ixbeg-1+ix .gt. nxend-istag ) then
!      if ( kz .gt. nz-kstag-1 .and. ix .ge. nx-istag) then
        if ( ixbeg-1+ix .eq. nxend-1 ) then
         go to 2200
        elseif ( ix .ge. nx ) then
         go to 2200
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz 
      end if
#else
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
        go to 2200
      else
        nzmpb = kz 
      end if
#endif

#ifdef MPI
      if ( ix .ge. nx-1 ) then
       if ( ixbeg-1+ix .eq. nxend-1 ) then
        nxmpb = 1
       elseif ( ix .ge. nx ) then
        nxmpb = 1
       else
        nxmpb = ix+1
       endif
      else
       nxmpb = ix+1
      end if
#else
      if ( ix+1 .gt. nx-1 ) then
       nxmpb = 1
      else
       nxmpb = ix+1
      end if
#endif

 2000 continue
 2200 continue
!
!  end of gather scatter (for this jy slice)


!
!
!  end of jy loop
!
!
!
29999  continue
!
!
      ENDIF ! ( ipconc .ge. 2 .and. iba .eq. 1 )


 4999  continue
!
!

      IF ( iptotal .gt. 0 ) THEN
        write(iunit,*) 'YIKES! Violated eqtot ',iptotal,' times!'
        write(iunit,*) 'ptotalmn,ptotalmx = ',ptotalmn,ptotalmx
      END IF

      call flush(iunit)
      
!
!
!  set save counter:  nsvcnt
!
!      
      nsvcnt = nsvcnt + 1
!
!/////////////////////////////////////////////////////////////////////
!/////  save microphysical data at interval nsvmp                /////
!/////////////////////////////////////////////////////////////////////
!
!
!  end of nsave
!
!

!
!
!
      if ( nstep/iholen*iholen .eq. nstep ) then
!
!  hole filling routine: done at interval iholen 
!
!  zero values
!
      cqtotn = 0.0
      cctotn = 0.0
!      citotn = 0.0
      crtotn = 0.0
!      cstotn = 0.0
      cvtotn = 0.0
!      chtotn = 0.0
      cqtotp = 0.0
      cctotp = 0.0
!      citotp = 0.0
      crtotp = 0.0
!      cstotp = 0.0
      cvtotp = 0.0
!      chtotp = 0.0
      cqfac = 0.0
      ccfac = 0.0
!      cifac = 0.0
      crfac = 0.0
!      csfac = 0.0
      cvfac = 0.0
!      chfac = 0.0
!
!  sum to find totals (negative and positive)
!
      if ( iholef .eq. 1 ) then

#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77611 kz = kzb,kze
      do 77612 jy = jyb,jye
      do 77613 ix = ixb,ixe
#else
! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotn)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : cqtotn)
      do 77611 kz =  1, nz-kd1
      do 77612 jy =  1, ny-jd1
      do 77613 ix =  1, nx-id1
#endif
      dv = dxx(ix)*dyy(jy)*dzz(kz)
      cqtotn = cqtotn + dn(ix,jy,kz)*dv* &
     &               ( an(ix,jy,kz,lc) &
     &                + an(ix,jy,kz,lr) &
     &                + an(ix,jy,kz,lv) )
77613 continue
77612 continue
77611 continue

!  set negatives to zero

#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77631 kz = kzb,kze
      do 77632 jy = jyb,jye
      do 77633 ix = ixb,ixe
#else
! cmic$ do all autoscope
! cmic$1 shared (an)
! cmic$2 private (kz,jy,ix)
! C$DOACROSS  LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do 77631 kz =  1, nz-kd1
      do 77632 jy =  1, ny-jd1
      do 77633 ix =  1, nx-id1
#endif
      an(ix,jy,kz,lc) = max(an(ix,jy,kz,lc),0.0)
      an(ix,jy,kz,lr) = max(an(ix,jy,kz,lr),0.0)
      an(ix,jy,kz,lv) = max(an(ix,jy,kz,lv),0.0)
      IF ( ipconc .ge. 2 ) THEN
      an(ix,jy,kz,lnc) = max(an(ix,jy,kz,lnc),0.0)
      ENDIF
      IF ( ipconc .ge. 3 ) THEN
      an(ix,jy,kz,lnr) = max(an(ix,jy,kz,lnr),0.0)
      ENDIF
77633 continue
77632 continue
77631 continue
!  resum values
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77651 kz = kzb,kze
      do 77652 jy = jyb,jye
      do 77653 ix = ixb,ixe
#else
! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotp)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : cqtotp)
      do 77651 kz =  1, nz-kd1
      do 77652 jy =  1, ny-jd1
      do 77653 ix =  1, nx-id1
#endif
      dv = dxx(ix)*dyy(jy)*dzz(kz)
      cqtotp = cqtotp + dn(ix,jy,kz)*dv* &
     &               ( an(ix,jy,kz,lc) &
     &                + an(ix,jy,kz,lr) &
     &                + an(ix,jy,kz,lv) )
77653 continue
77652 continue
77651 continue
!  set factors based on ratio of (totn+a)/(totp+a) where a =1.e-20
!
      cqfac = &
     & (cqtotn+1.e-20)/(cqtotp+1.e-20)
      write(iunit,*) 'CQFAC = ',cqfac,cqtotn,cqtotp
      
      IF ( .not. ( cqfac .lt. 2.0 .and. cqfac .gt. 0.0 ) ) THEN
        write(iunit,*) 'Blowing up! STOP!'
        write(0,*) 'Blowing up! STOP!'
        
        DO ia = 1,na
#ifdef MPI
        kzb = 1
        kze = ktile
        if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

        jyb = 1
        jye = jtile
        if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

        ixb = 1
        ixe = itile
        if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

        DO kz = kzb,kze
        DO jy = jyb,jye
        DO ix = ixb,ixe
#else
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
#endif
          IF ( .not. ( an(ix,jy,kz,ia)  .lt. 1.0e10  .and.  &
     &              an(ix,jy,kz,ia)  .gt. -1.0e10  )   ) THEN
           write(iunit,'(4i5,1x,1pe13.5)')  &
     &      ix,jy,kz,ia,an(ix,jy,kz,ia)
          END IF
        END DO
        END DO
        END DO
        END DO
        
        STOP
      END IF
!
!  now reset all values based on factors
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77671 kz = kzb,kze
      do 77672 jy = jyb,jye
      do 77673 ix = ixb,ixe
#else
! cmic$ do all autoscope
! cmic$1 shared (an,cqfac)
! cmic$2 private (kz,jy,ix)
! C$DOACROSS  LOCAL(kz,jy,ix), share(cqfac,an)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do 77671 kz =  1, nz-kd1
      do 77672 jy =  1, ny-jd1
      do 77673 ix =  1, nx-id1
#endif
      an(ix,jy,kz,lc) = an(ix,jy,kz,lc)*cqfac
      an(ix,jy,kz,lr) = an(ix,jy,kz,lr)*cqfac
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv)*cqfac
77673 continue
77672 continue
77671 continue
!
      end if
!
!  end of hole filling version 1
!
!
!  sum to find totals (negative and positive)
!
      if ( iholef .eq. 2 ) then
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77711 kz = kzb,kze
      do 77712 jy = jyb,jye
      do 77713 ix = ixb,ixe
#else
      do 77711 kz =  1, nz-kd1
      do 77712 jy =  1, ny-jd1
      do 77713 ix =  1, nx-id1
#endif
      cctotn = cctotn + dn(ix,jy,kz)*an(ix,jy,kz,lc)
      crtotn = crtotn + dn(ix,jy,kz)*an(ix,jy,kz,lr)
      cvtotn = cvtotn + dn(ix,jy,kz)*an(ix,jy,kz,lv)
77713 continue
77712 continue
77711 continue
!
!  set negatives to zero
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77731 kz = kzb,kze
      do 77732 jy = jyb,jye
      do 77733 ix = ixb,ixe
#else
      do 77731 kz =  1, nz-kd1
      do 77732 jy =  1, ny-jd1
      do 77733 ix =  1, nx-id1
#endif
      an(ix,jy,kz,lc) = max(an(ix,jy,kz,lc),0.0)
      an(ix,jy,kz,lr) = max(an(ix,jy,kz,lr),0.0)
      an(ix,jy,kz,lv) = max(an(ix,jy,kz,lv),0.0)
77733 continue
77732 continue
77731 continue
!77740 continue
!
!  resum values
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77751 kz = kzb,kze
      do 77752 jy = jyb,jye
      do 77753 ix = ixb,ixe
#else
      do 77751 kz =  1, nz-kd1
      do 77752 jy =  1, ny-jd1
      do 77753 ix =  1, nx-id1
#endif
      cctotp = cctotp + dn(ix,jy,kz)*an(ix,jy,kz,lc)
      crtotp = crtotp + dn(ix,jy,kz)*an(ix,jy,kz,lr)
      cvtotp = cvtotp + dn(ix,jy,kz)*an(ix,jy,kz,lv)
77753 continue
77752 continue
77751 continue
!
!  set factors based on ratio of (totn+a)/(totp+a) where a =1.e-20
!
!      csfac =
!     > (cstotn+1.e-20)/(cstotp+1.e-20)
      crfac = &
     & (crtotn+1.e-20)/(crtotp+1.e-20)
!      cifac =
!     > (citotn+1.e-20)/(citotp+1.e-20)
      ccfac = &
     & (cctotn+1.e-20)/(cctotp+1.e-20)
      cvfac = &
     & (cvtotn+1.e-20)/(cvtotp+1.e-20)
!      chfac =
!     > (chtotn+1.e-20)/(chtotp+1.e-20)
!
!  now reset all values based on factors
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      do 77771 kz = kzb,kze
      do 77772 jy = jyb,jye
      do 77773 ix = ixb,ixe
#else
      do 77771 kz =  1, nz-kd1
      do 77772 jy =  1, ny-jd1
      do 77773 ix =  1, nx-id1
#endif
      an(ix,jy,kz,lc) = an(ix,jy,kz,lc)*ccfac
      an(ix,jy,kz,lr) = an(ix,jy,kz,lr)*crfac
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv)*cvfac
77773 continue
77772 continue
77771 continue
!
      end if
!
!  end of hole filling version 2
!        
      endif

#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
            t0(ix,jy,kz) = an(ix,jy,kz,lt)*t77(ix,jy,kz)

        ENDDO
       ENDDO
      ENDDO


      IF ( ipconc .ge. 2 ) THEN
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
         temp1 = t0(ix,jy,kz)

         ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

!          c1 = t00(ix,jy,kz)*tabqvs(ltemq)
          c1 = t00(ix,jy,kz)*tabqis(ltemq)

          ssati(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values

          c1 = t00(ix,jy,kz)*tabqvs(ltemq)

          ssfilt(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
          IF ( ny .eq. 2 .and. ssfilt(ix,jy,kz) .gt. 1.0 &
     &                             .and. ncdebug .ge. 1 ) THEN
          write(iunit,*) 'ssfilt2 = ',ssfilt(ix,jy,kz), &
     &     an(ix,jy,kz,lv)*1.e3,c1*1.e3,an(ix,jy,kz,lt),ix,kz, &
     &     an(ix,jy,kz,lnc)*1.e-6,an(ix,jy,kz,lc)*1.e3
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      ENDIF
! 
!
! Redistribute inappreciable cloud particles and charge
!
! Redistribution everywhere in the domain...
!
      sctot1n = 0.0
      sctot1p = 0.0
      frac = 1.0 ! 0.2
! cmic$ do all autoscope
! cmic$1 private(kz,jy,ix),
! cmic$2 shared(an)
! C$DOACROSS LOCAL(kz,jy,ix),
! C$&        SHARE(an)      
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),SHARED(an)  

      IF ( ipconc .lt. 3 ) THEN
!      IF ( .true. ) THEN

#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
      do jy = 1,ny
      do ix = 1,nx
#endif

!      IF ( an(ix,jy,kz,lr)  .gt. 10.0e-3 )
!     :  print*, 'RAIN3: ',ix,kz,an(ix,jy,kz,lr)

      if ( an(ix,jy,kz,lr) .lt. frac*qxmin(lr) &
!     :     .and. ( an(ix,jy,kz,lr) .ge. 0.0 .or.
!     :       an(ix,jy,kz,lscr) .ne. 0.0 )  &
     &  ) then
        an(ix,jy,kz,lc) = an(ix,jy,kz,lc) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        IF ( ipconc .ge. 3 ) THEN
          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
        ENDIF
      
      end if

      
!
!  for qcw
!
      IF ( an(ix,jy,kz,lc) .le. frac*qxmin(lc) &
!     :    .and. ( an(ix,jy,kz,lc) .gt. 0.0 .or.
!     :       an(ix,jy,kz,lscw) .ne. 0.0 )  &
     &       ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0

        IF ( lccn .lt. 1 ) THEN
         an(ix,jy,kz,lccn) =   &
     &      Min( ccwmx, an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc)) )
        ENDIF
         an(ix,jy,kz,lnc) = 0.0
       
       IF ( idecss .eq. 1  &
     &     .and. an(ix,jy,kz,lss) .gt. 0.0 ) THEN
         an(ix,jy,kz,lss) = an(ix,jy,kz,lss)*Exp(-dtp/10.)
         IF ( ssat(ix,jy,kz) .lt. -5.0 ) THEN
          an(ix,jy,kz,lss)  = 0.0 
         ENDIF 
       ENDIF
       
       ENDIF


      end do
      end do
      end do
      
      ELSE
!
!  alternate test version for ipconc .ge. 3
!  just vaporize stuff to prevent noise in the number concentrations
!
#ifdef MPI
      kzb = 1
      kze = ktile
      if(kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if(jyend .eq. nyend) jye = jyend-jybeg+1

      ixb = 1
      ixe = itile
      if(ixend .eq. nxend) ixe = ixend-ixbeg+1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,nz-1
      do jy = 1,ny
      do ix = 1,nx
#endif

!      IF ( an(ix,jy,kz,lr)  .gt. 10.0e-3 )
!     :  print*, 'RAIN3: ',ix,kz,an(ix,jy,kz,lr)

      if ( an(ix,jy,kz,lr) .lt. frac*qxmin(lr) &
     &  ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        IF ( ipconc .ge. 3 ) THEN
!          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
         IF ( lzr > 1 )  an(ix,jy,kz,lzr) = 0.0
        ENDIF
      
      end if

      
!
!  for qcw
!
      IF ( an(ix,jy,kz,lc) .le. frac*qxmin(lc) &
!     :    .and. ( an(ix,jy,kz,lc) .gt. 0.0 .or.
!     :       an(ix,jy,kz,lscw) .ne. 0.0 )  &
     &       ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0
       IF ( ipconc .ge. 2 ) THEN
        IF ( lccn .lt. 1 ) THEN
         an(ix,jy,kz,lccn) =   &
     &      Min( ccwmx, an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc)) )
        ENDIF
         an(ix,jy,kz,lnc) = 0.0
       
       IF ( idecss .eq. 1  &
     &     .and. an(ix,jy,kz,lss) .gt. 0.0 ) THEN
         an(ix,jy,kz,lss) = an(ix,jy,kz,lss)*Exp(-dtp/10.)
         IF ( ssat(ix,jy,kz) .lt. -5.0 ) THEN
          an(ix,jy,kz,lss)  = 0.0 
         ENDIF 
       ENDIF
       
       ENDIF

      ENDIF

      end do
      end do
      end do
      
      ENDIF

!
!
!
!  return   
!     

      deallocate ( db1 )
      deallocate ( db0 )
      deallocate ( dtz1 )
      deallocate ( dtz0 )
      deallocate ( xvt )
      deallocate ( tmpn )
      deallocate ( tmpn2 )
      deallocate ( ztmp )
      
      return
      end
!
!  end of subroutine
!
!--------------------------------------------------------------------------
!
!  Bottom of code
!
