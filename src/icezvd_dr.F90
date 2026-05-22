#define SWM

#ifdef CM1  /* special setup for CM1 */

#undef CHGELEC

#else

#include "sam.def.h"
#define CHGELEC

#endif

#ifdef NOELEC
#undef CHGELEC
#endif

!#include "sam.def.h"
!#define ICE10
!#define CHGELEC
!#define SAM
!
! Things to do:
!
!  Test using exponential formulation for rain fall speed.  If there is little change
!  from the quadratic, it would be less complicated to use.
!
!  Contact nucleation needs to be fixed up to be similar to Cotton et al. 1986 and Meyers et al 1992.
!
! The following are done?
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
!    Replace qv0 with qx(mgs,lv)? No. qv0 is base val
!
! Need to look at limiting supersaturation to 1 or so by nucleation/condensation
!
!  put in temperature-dependent function for homogeneous freezing
!
!cc--------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------
!
      subroutine ice_zieg    & 
     &  (ntmul,nx,ny,nz,na,nba,nv,nstep,time_real & 
     &  ,nor,norz,istag,jstag,kstag,itopo & 
     &  ,iwrite,io_flag & 
     &  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d & 
     &  ,xfall,rate2d,nrate2d,xfalltot & 
     &  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,tq1,ntq1 & 
     &  ,ab,an,db,dn,p2,pinit & 
     &  ,pb,pn,u,v,w,km,iunit & 
     &  ,ssat,ssfilt,t00,t77,ssati,thproc,numproc &
#ifndef CM1
     &  ,elec, axtra &
#endif
     &  ,bcx,bcy)

#ifndef CM1
       USE GRID_MODULE
       USE CPUTIME_MODULE
#endif
       USE MICRO_MODULE
#ifdef CHGELEC
       USE ELEC_MODULE
#endif
       USE INDEX_MODULE !, only: lt,lc,lr,li,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lhw,lhlw, &
                        !       lnc,lnr,lni,lns,lnh,lnhl,cinu,dmuh,dnu,dmu,xnu,xmu,dmuhl,rnu,cnu,snu, &
                        !       lss,lsat,lsati,xvcmx,xvcmn,xvrmn,xvrmx,lccn,rnumin,rnumax, &
                        !       lqmx,nxtra,lvi,lvs,lvh,lvhl,lzi,lzs,lzr,lzh,lzhl,lsw,alphar,alphamin,alphamax, &
                        !       lscw,lscr,lsci,lscs,lsch,lschl,fx,alphas,alphah,alphahl,xvsmn,xvhmn,xvsmx,xvhmx, &
                        !       xvhlmn,xvhlmx,lmvr,lmvh,lmvhl,lmvs,ld0,lqaut,lqacc,lqrevap,lqcevap,lqcond, &
                        !       lqdep,lqmelt,lqsub,lccna
       USE COMMASMPI_MODULE
       USE TRAJ_MODULE
       USE PARAM_MODULE, only: len_type, Cm, Ce, Lmax, rdorv => epsilon
!       USE module_mp_nssl_2mom, only: hailmaxd
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
!  10/17/2006: added flag (iehw) to select how to calculate ehw
!
!  10/5/2006: switched chacr to integrated version rather than assuming that average rain
!             drop mass does not change.  This acts to reduce rain size somewhat via graupel
!             collection.
!             Use Mason data for ehw, with scaling toward ehw=1 as air density decreases.
!
!  10/3/2006: Turned off Meyers nucleation for T > -5 (can turn on with imeyers5 flag)
!             Turned off contact nucleation in updrafts
!
!  7/24/2006:  Turned on Meyers nucleation for -5 < T < 0
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
#ifndef CM1
      TYPE(VARIABLE)     :: elec(neelec)
#endif
      
      integer            :: bcx, bcy
      
      integer ng1
      integer,parameter :: iunit0 = 0 
      integer :: iunit
      parameter(ng1 = 1)
      
      real,allocatable :: db1(:,:),db0(:,:),z(:,:,:),dtz1(:,:),dtz0(:,:)
      real qvex
      
      integer iraincv, icgxconv
      parameter ( iraincv = 1, icgxconv = 1)

      real ccwtmp,ccitmp ! ,ciptmp,cirtmp
      
      double precision dp1
      
      real frac
      
      real dtp1
      integer ntmul
      
!      logical, parameter :: usenucond = .true.
!      logical, parameter :: usensslgs = .false.
!      logical, parameter :: useoldgs = .false.
      logical, parameter :: do1dsedimentation = .true.
! a few vars for time-split fallout      
      real vtmax, vtmaxall(lc:lhab), vtmaxz(lc:lhab), vtmaxq(lc:lhab),vtmaxn(lc:lhab)
      integer n,ndfall,m,ndfalltot
      
      double precision chgneg,chgpos,sctot
      
      real temgtmp
      integer  nstep ! , nstart, nstop
      real     time_real !  model time in seconds.
      integer nx,ny,nz,na,nba,nv
      integer ng
      integer nor,norz,mzdist,imapz,istag,jstag,kstag,itopo ! ,nht,ngt,igsr
      integer iwrite
      logical, intent(in) :: io_flag
      real dtp,dx,dy,dz
!      real dzc(nz)                         ! 1/dz(k)
      real dzz(nz),dxx(nx),dyy(ny)          ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-norz+ng1:nz+norz)
      real z1d(-norz+ng1:nz+norz,4)

      integer numproc
      real thproc(nzend,numproc)

      real dv

      real dtptmp

      integer itest,nidx,id1,jd1,kd1
      parameter (itest=1)
      parameter (nidx=10)
      parameter (id1=1,jd1=1,kd1=1)
      integer ierr
      integer iend,imake
      integer, save :: iread = 0
      save imake
      integer ix,jy,kz, il, ic, ir, icp1, irp1, in, i
      integer j,k
      real slope1, slope2
      real x1, x2

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
! Takahashi lookup table
!
      integer nlwc,ntem
      parameter (nlwc=30, ntem=31)
      real takalu(0:ntem,0:nlwc) 
      save takalu
!
!  Other elec. vars
!
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      integer lvol(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
      integer lrain(nraintypes)
      integer, save :: linfall(lc:lqmx)
      
      logical ldovol
      logical zeroout(lc:lhab)

!
!  fallout arrays
!
      real xfall(nx,ny,na)
      integer nrate2d
      real rate2d(nx,ny,nrate2d)
      real xfall0(nx,ny)
      real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,nprecip+neelec2d)

!
! params read in from inmicro
!
      real dtfac
      parameter ( dtfac = 1.0 )

      
      real cckm,ccne,ccnefac,cnexp,ccne0
      save cckm,ccne,ccnefac,cnexp,ccne0
      real beta
!
!      integer nsave
!
      integer ido(lc:lqmx)
      save ido
      
      integer iexy(lc:lqmx,lc:lqmx)

      save iexy
      
!
! external temporary arrays
!
      real t00(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t77(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

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
      
      integer :: ntq1
      real tq1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,ntq1)  !qproc temporary

      real t8s, t9s


      real pinit(-norz+ng1:nz+norz)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi
      real pb(-norz+ng1:nz+norz)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real db(-norz+ng1:nz+norz)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
#ifndef CM1
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)
#endif
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real ab(-norz+ng1:nz+norz,na)

      real u(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real v(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real km(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

      real ssati(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)


      real ssat(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real ssfilt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
! 
!  declarations microphysics and for gather/scatter
!
      integer nxmpb,nzmpb,nxz
      integer jgs,mgs,ngs,numgs
      parameter (ngs=500)
      integer ntt
      parameter (ntt=300)

      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs)
      integer nsvcnt

      integer ncuse
      parameter (ncuse=0)
      integer il0(ngs),il5(ngs)

!      integer nqsat
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati
!      parameter (fqsat=0.002,fqsati=1./fqsat)
!      parameter (fqsat=0.01,fqsati=1./fqsat)

!      real cai,caw,cbi,cbw
      real tdtol,temsav,tfrcbw,tfrcbi,thnuc
      
      real tfr,tfrh
      parameter ( tfr = 273.15, tfrh = 233.15)
      
      real cp, rd
      parameter ( cp = 1004.0, rd = 287.04 )
      
      real cpi
      parameter ( cpi = 1./cp )
      
      real poo,cap
      parameter ( cap = rd/cp, poo = 1.0e+05 )


!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs), ccna(ngs)
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
      real temp1,temp2 ! ,ssold
      real :: ssmax(ngs) = 0.0      ! maximum SS experienced by a parcel
      real ssmx
      real dnnet,dqnet
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real bfnu
      parameter ( bfnu = (rnu + 2.0)/(rnu + 1.0)  )
      real ventr, ventrn, ventc
      save ventr, ventrn, ventc
      real ventrx(ngs)
      real ventrxn(ngs)
      real volb, t2s, aa1, aa2
      parameter ( aa1 = 9.44e15, aa2 = 5.78e3 ) ! a1 in Ziegler
! snow parameters:
      real cexs, cecs
      parameter ( cexs = 0.1, cecs = 0.5 )
      real rvt      ! ratio of collection kernels (Zrnic et al, 1993)
      parameter ( rvt = 0.104 )

      real ec0, ex1, ft, rhoinv(ngs),rho_qx
      
      real chw, g1, rd1

      
!      integer kbound
      real ac1,bc, taus, c1,d1,e1,f1,p380,tmp,tmp2 ! , sstdy, super
      real tmpmx, fw
      real x,y,del,r,rtmp,alpr
      real :: alpjj, alpii, xnuii, xnujj
      integer :: ii, jj
      double precision :: vent1,vent2
      real g1palp
      real bs
      real v1, v2
      real d1r, d1i, d1s, e1i
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

      real ssi1, ssi2, dqvi, dqvis, dqvii,qis1
      real dqvr, dqc, dqr, dqi, dqs
      real qv1m,qvs1m,ss1m,ssi1m,qis1m
      real cwmastmp 
      real  dcloud,dcloud2 ! ,as, bs
      real cn(ngs) 

      real mwfac

      real  es(ngs) ! ss(ngs),
      real  eis(ngs)
      real ssbar(ngs) !, sqsat(ngs)
      real ssf(ngs),ssfkp1(ngs),ssfkm1(ngs),ssat0(ngs)
      real ssfjp1(ngs),ssfjm1(ngs)
      real ssfip1(ngs),ssfim1(ngs)
      real supcb, supmx
      parameter (supcb=0.5,supmx=238.0)
      real r2dxm, r2dym, r2dzm
      real dssdz, dssdy, dssdx

      real xvmn(lc:lhab), xvmx(lc:lhab)

      real rwmasn,rwmasx

      real epsi,d
      parameter (epsi = 0.622, d = 0.266)
      real r1,qevap ! ,slv
      
      real vr,nrx,qr,z1,z2,rdi,alp,xnutmp,xnuc
      real ctmp

      integer, save :: infdo = 1

!      real svc(ngs)  !  droplet volume
!
!  misc
!
      integer :: ni,nj,nk
      real :: nr,d0
      
!      real tabqvs(nqsat),tabqis(nqsat),dtabqvs(nqsat),dtabqis(nqsat)
!      save tabqvs,tabqis,dtabqvs,dtabqis
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
!      real delqci(ngs) ! ,delqip(ngs)
!      real temp(ngs),tempc(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) ! ,tembzg(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
!      real elv(ngs),elf(ngs),els(ngs),elvs(ngs),elss(ngs)
      real gamw(ngs),gams(ngs)   !   qciavl(ngs),
      real tsqr(ngs),ssi(ngs),ssw(ngs)
      real cc3(ngs),cqv1(ngs),cqv2(ngs)
      real f5, prod0, qvs0  ! Kessler condensation factor
      real cwc1
      save cwc1
!      real qcwdif(ngs) ! ,dcwnc
      real, allocatable, save :: cwnccnold(:)
      real :: cwnccn(ngs)
      real qcwtmp(ngs),qtmp,qtot(ngs) ! ,cwnc(ngs)
!      real cwmasn,cwmasx
!      real cwmasn5
!      save cwmasn5
!      real cwradn
      real cimasn,cimasx,ccimx
      real pi,pid4
      real ar,br,cs,ds !,gf7,gf6,gf5,gf4,gf3,gf2,gf1
!      real gf73rds, gf83rds
!      real gf43rds, gf53rds
!      real aradcw,bradcw,cradcw,dradcw,cwrad,rwrad,rwradmn
!      parameter ( rwradmn = 50.e-6 )
      
      real cionp(ngs),cionn(ngs)
!
!  other arrays
!
      
      
!      real fwet1(ngs),fwet2(ngs)   !   ,fwet3(ngs)
!      real fmlt1(ngs),fmlt2(ngs)   !   ,fmlt3(ngs)
      real fvds(ngs) ! ,fvce(ngs),fiinit(ngs) ! ,fcinit(ngs)
      real fvent(ngs) !,fraci(ngs),fracl(ngs)
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
      real fcci(ngs), fcip(ngs)
!
       real qitmp(ngs)
      
!      real qxmin(lc:lhab)
      real qxmin(lc:lqmx)
      save qxmin

!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)
      real :: axx(ngs,lh:lhab),bxx(ngs,lh:lhab)
      
      real :: qx(ngs,lv:lhab)
      real :: qxw(ngs,ls:lhab)
      real :: cx(ngs,lc:lhab)
      real :: xv(ngs,lc:lhab)
      real :: vtxbar(ngs,lc:lhab,3)
      real :: xmas(ngs,lc:lhab)
      real :: xdn(ngs,lc:lhab)
      real :: cdxgs(ngs,lc:lhab)
      real :: xdia(ngs,lc:lhab,3)
      real :: vx(ngs,li:lhab)
      real :: alpha(ngs,lc:lhab)
      real :: zx(ngs,lr:lhab)

!
      real swvent(ngs),hwvent(ngs),rwvent(ngs)
      real civent(ngs)
!
!      real xivent(ngs,nhab),xix(ngs,nhab),xirey(ngs,nhab)
!      real xpvent(ngs)
!
!
!      real amccw(ngs) !,amcci(ngs),amcrw(ngs),amcsw(ngs),amchw(ngs)
!
!
      real fsczz(ngs)
      real fschw(ngs),fscsw(ngs)
      real fsccw(ngs),fscci(ngs),fscrw(ngs)
!
      
      real, save :: xdnmx(lc:lqmx), xdnmn(lc:lqmx)
!      real vtxbar(ngs,nhab),ximas(ngs,nhab),xidn(ngs,nhab)
!
!
      real cilen(ngs) ! ,ciplen(ngs)
!
!
      real rwcap(ngs),swcap(ngs)
      real hwcap(ngs)
      real cicap(ngs)
!
!  conversions
!
      real cnina(ngs),wvel(ngs),wvelkm1(ngs)
      real uvel(ngs),vvel(ngs)
!
      real qidpv(ngs),qisbv(ngs)
      real qidsv(ngs),qidsvp(ngs)
!
!
      real qsdpv(ngs),qssbv(ngs),qsdsv(ngs)
!

      
      real da0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real dab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real dab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real da1 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real bb  (lc:lqmx)

      save da0, dab0, dab1, da1, bb

! for 3-moment collection coefficients
      real, save :: dab0lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real, save :: dab1lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005

      real va0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real vab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real vab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real va1 (lc:lqmx)          ! collection coefficients from Seifert 2005

      save va0, vab0, vab1, va1
      
!      real alpha(ngs) ! shape parameter
      
!      save alpha
      
      real exy(ngs,ls:lhab,lc:ls)
      real scxacy(ngs,ls:lhab,lc:ls)
      real cxacy(ngs,ls:lhab,lc:ls)
      real sxxacy(ngs,ls:lhab,li:ls)

!
!  arrays for production terms
!
      real ptotal(ngs) ! , pqtot(ngs) 
     
!
!
!  other arrays
!
!
!      real wvdf(ngs),tka(ngs) !,akvisc(ngs),ci(ngs),cw(ngs),thdf(ngs)
      real dqisdt(ngs) !,advisc(ngs) !dqwsdt(ngs), ,schm(ngs),pndl(ngs)
      
!      save wvdf,tka,advisc

      real qss0(ngs)

      real advisc0,advisc1,tka0

      real qsacip(ngs)
      real pres(ngs),presp(ngs),pres0(ngs)
      real pk(ngs)
      real rho0(ngs),pi0(ngs),piz(ngs)
      save piz
      real rhovt(ngs)
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
!      real pceds(ngs) ! ,ppceds(ngs),pmceds(ngs)
      real times
!      real qwfzi(ngs) ! ,qimlw(ngs)
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
      parameter (iholen = 9)
      real  cqtotn,cqtotn1
      real  cctotn
      real  citotn
      real  crtotn
      real  cstotn
      real  cvtotn
      real  cftotn
      real  cgltotn
      real  cghtotn
      real  chtotn
      real  cqtotp,cqtotp1
      real  cctotp
      real  citotp
      real  ciptotp
      real  crtotp
      real  cstotp
      real  cvtotp
      real  cftotp
      real  chltotp
!      real  chxtotp
      real  cgltotp
      real  cgmtotp
      real  cghtotp
      real  chtotp
      real  cqfac
      real  ccfac
      real  cifac
      real  cipfac
      real  crfac
      real  csfac
      real  cvfac
      real  cffac
      real  cglfac
      real  cghfac
      real  chfac
      
      real ssifac, qvapor
!
      
      real,allocatable :: xvt(:,:,:,:) ! (nx,nz,2,lc:lhab) ! 1=mass-weighted, 2=number-weighted
      real,allocatable :: tmpn(:,:,:)
      real,allocatable :: tmpn2(:,:,:)
  !    real, dimension(its:ite, 1, kts:kte, 3) :: alpha2d
       real,allocatable :: alpha2d(:,:,:,:)
       integer :: nalpha2d
!     
!
!
!   Miscellaneous variables
!
      integer ireadqf,lrho,lqsw,lqgl,lqgm ,lqgh ! ,ltim,ltem,lqcw,lqfw
      integer lqrw,imkgam 
      real vt
      data imkgam /0/
      integer igam        
      real arg
      real, external :: gamma  ! gamma is a function  
      real gaml02, gaml02d500
      real erbnd1, fdgt1, costhe1
      real qeps
      real dyi2,dzi2,cp608,cv,bta1,cnit,dragh,dnz00,rho00,pii
      real qccrit,gf4br,gf4ds,gf4p5, gf3ds, gf1ds,gr

      real hwdn, tmpg
      
      real, save :: xdn0(lc:lqmx)

      integer l ,ltemq,inumgs, idelq ! , ib
      
      integer, save :: ids2li=0 !=0 downscale to small ions only, =1 ds all but cloud water to large ion, =2 ds all

      real c1f3,rw,temq ! ,cmn,cmi40,cmi50
!      real ri50,vti50,bsfw,cm50a,a,cm40b,cm50b
      real cnin20,cnin10,cnin1a
      real cnin2a,cnin2b,ssival,tqvcon
      real cd(5)
      real cdx(lc:lhab)
      real cno(lc:lhab), cnox
      real cval,aval,eval,fval,gval ,qsign,ftelwc,qconkq
      real qconm,qconn,cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5,qdiff,argrcnw
      real c4,bradp,bl2,bt2,dtrh,hrifac, hdia0,hdia1,civenta,civentb
      real civentc,civentd,civente,civentf,civentg,cireyn,xcivent
      real cipventa,cipventb,cipventc,cipventd,cipreyn,cirventa
      real cirventb
      integer igmrwa,igmrwb,igmswa, igmswb,igmfwa,igmfwb,igmhwa,igmhwb
      real rwventa ,rwventb,swventa,swventb,fwventa,fwventb,fwventc
      real hwventa,hwventb 
      integer iptotal ! counts number of times that eqtol is exceeded
      integer ilock1,ilock2,ilock3,ilock4,ilockc1,ilockc2
      real ptotalmx,ptotalmn
      real chlmax
      real psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p
      integer nscmax,nscmin,iwetg, iwetg1
      
      real ptotalmxy(ny),ptotalmny(ny)
      real psctotmxy(ny),psctotmny(ny),psctot1y(ny),psctot2y(ny)
      real sctot1ny(ny),sctot1py(ny)
      integer iwetgy(ny),iwetg1y(ny),iptotaly(ny),nscmaxy(ny),nscminy(ny)
      integer, allocatable :: ngscmaxy(:,:,:),ngscminy(:,:,:),ngscxymax(:,:,:)
      real :: dbzchange(2,nz)

      real    hwventc, hlventa, hlventb,  hlventc
      real  glventa, glventb, glventc 
      real   gmventa, gmventb,  gmventc, ghventa, ghventb, ghventc 

      real  dzfacp,  dzfacm,  cmassin,  cwdiar ! , cwmasr
      real  rimmas, rhobar
      real   argtim, argqcw, argqxw, argtem
      real    a1,a2,a3,a4,a5,a6
      real   gamss
      integer  itertd, ia

      real delbk, delabk, delvbk

      real mlen, dyl
      real tkediss(-nor+1:nx+nor,-norz+ng1:nz+norz)

#ifndef MPI
! !$    volatile iptotal,ptotalmx,ptotalmn,nscmax,nscmin,iwetg,  &
! !$   &   psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p
#endif
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      real    :: t0p1, t0p3


#ifdef MPI
      INCLUDE "mpif.h"
      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      integer :: downward_tag, upward_tag

      logical :: debug_mpi = .false.
      integer       :: nampi, nb

      integer, parameter :: ntot = 50
      real  mpitotin(ntot), mpitotout(ntot)
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

#endif

!
! ####################################################################
!
!  Start routine
!
! ####################################################################
!
#ifndef CHGELEC
      ipelec = 0
#endif


          ni = ixend-ixbeg+1
          IF ( ixend == nxend ) ni = ixend-ixbeg
          nj = jyend-jybeg+1
          IF ( jyend == nyend ) nj = jyend-jybeg
          IF ( ny == 2 ) nj = 1
          nk = nzend-1

!
!  slope intercepts
!
!      ipconc = 5

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG START"
#endif
      
      IF ( ngs .lt. nz ) THEN
       write(0,*) 'Error in ICEZVD: Must have ngs .ge. nz!  ngs,nz = ',ngs,nz
       STOP
      ENDIF
      
      allocate ( db1(nx,nz+1) )
      allocate ( db0(nx,nz+1) )
      allocate ( dtz1(nx,nz+1) )
      allocate ( dtz0(nx,nz+1) )
      allocate (   z(-nor+ng1:nx+nor,-norz+ng1:nz+norz,lr:lhab) )
      allocate ( xvt(nx,nz+1,3,lc:lhab) )
      allocate ( tmpn (-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate ( tmpn2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
        tmpn (:,:,:) = 0.0
        tmpn2(:,:,:) = 0.0
      IF ( .not. allocated( alpha2d ) ) THEN
        IF ( lf > 1 ) THEN
          nalpha2d = 4
        ELSE
          nalpha2d = 3
        ENDIF
        allocate ( alpha2d(-nor+1:nx+nor,1,-norz+ng1:nz+norz,nalpha2d) )
      ENDIF
      
      allocate ( ngscminy(lg:lhab,lc:ls,ny) )
      allocate ( ngscmaxy(lg:lhab,lc:ls,ny) )
      allocate ( ngscxymax(lg:lhab,lc:ls,ny) )

      ng = nor
      
      lsc(lc) = lscw
      lsc(lr) = lscr
      lsc(li) = lsci
      lsc(ls) = lscs
      lsc(lh) = lsch
      IF ( lhl .gt. 1 ) lsc(lhl) = lschl
      IF ( lf .gt. 1 ) lsc(lf) = lscf
      IF ( lis > 1 ) lsc(lis) = lscis

      ln(lc) = lnc
      ln(lr) = lnr
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      IF ( lhl .gt. 1 ) ln(lhl) = lnhl
      IF ( lf .gt. 1 ) ln(lf) = lnf
      IF ( lis > 1 .and. lnis > 1 ) ln(lis) = lnis

      ipc(lc) = 2
      ipc(lr) = 3
      ipc(li) = 1
      ipc(ls) = 4
      ipc(lh) = 5
      IF ( lhl .gt. 1 ) ipc(lhl) = 5
      IF ( lf .gt. 1 ) ipc(lf) = 5
      IF ( lis > 1 ) ipc(lis) = 1
      
      ldovol = .false.
      lvol(:) = 0
      lvol(li) = lvi
      lvol(ls) = lvs
      lvol(lh) = lvh
      IF ( lhl .gt. 1 .and. lvhl .gt. 1 ) lvol(lhl) = lvhl
      IF ( lf .gt. 1 .and. lvf .gt. 1 ) lvol(lf) = lvf
      
      DO il = lc,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO
      
      lz(:) = 0
      lz(lr) = lzr
      lz(li) = lzi
      lz(ls) = lzs
      lz(lh) = lzh
      IF ( lhl .gt. 1 .and. lzhl > 1 ) lz(lhl) = lzhl
      IF ( lf .gt. 1 .and. lzf > 1 ) lz(lf) = lzf
      
      lliq(:) = 0
      lliq(ls) = lsw
      lliq(lh) = lhw
      IF ( lhl .gt. 1 ) lliq(lhl) = lhlw
      IF ( lf .gt. 1 ) lliq(lf) = lfw
      
      lrain(:) = 0
      IF ( iraintypes >= 1 ) THEN
        il = 0
        IF ( lrauto > 1 ) THEN
          il = il + 1
          lrain(il) = lrauto
        ENDIF  
        IF ( lrshed > 1 ) THEN
          il = il + 1
          lrain(il) = lrshed
        ENDIF  
        IF ( lrmelt > 1 ) THEN
          il = il + 1
          lrain(il) = lrmelt
        ENDIF
        IF ( il /= nraintypes ) THEN
          write(0,*) 'Problem with rain types (il /= nraintypes)! il, nraintypes = ', il,nraintypes
          call commasmpi_abort()
        ENDIF
        IF ( ndebug > 1 ) THEN
          DO m = 1,il
            write(0,*) 'raintypes: il,lrain = ',m,lrain(m)
          ENDDO
        ENDIF
      ENDIF

      bx(lr) = 0.85
      ax(lr) = 1647.81
      fx(lr) = 135.477
      
      IF ( icdx == 6 ) THEN
        bx(lh) = 0.6 ! Milbrandt and Morrison (2013) for density of 550.
        ax(lh) = 157.71
      ELSEIF ( icdx < 0 ) THEN
        bx(lh) = bxh
        ax(lh) = axh
      ELSEIF ( icdx > 1 ) THEN
        bx(lh) = 0.5
        ax(lh) = 75.7149
      ELSE ! ( icdx == 0 ) THEN
        bx(lh) = 0.37 ! 0.6  ! Ferrier 1994 graupel
        ax(lh) = 19.3
!      ELSE ! icdx == 0
!        ax(lh) = 206.984 ! Ferrier 1994 hail/frozen drops
!        bx(lh) = 0.6384
      ENDIF
!      bx(lh) = 0.6
      
      IF ( lhl .gt. 1 ) THEN
        IF ( icdxhl == 6 ) THEN
          bx(lhl) = 0.593 ! Milbrandt and Morrison (2013) for density of 750.
          ax(lhl) = 179.36
        ELSEIF (icdxhl < 0 ) THEN
         bx(lhl) = bxhl
         ax(lhl) = axhl
         
        ELSEIF (icdxhl > 1 ) THEN
         bx(lhl) = 0.5
         ax(lhl) = 75.7149
        ELSE ! == 0
          ax(lhl) = 206.984 ! Ferrier 1994
          bx(lhl) = 0.6384
        ENDIF
      ENDIF

      IF ( lf .gt. 1 ) THEN
       IF ( icdx == 6 ) THEN
!          bx(lf) = 0.593 ! Milbrandt and Morrison (2013) for density of 750.
!          ax(lf) = 179.36
          bx(lf) = 0.6 ! Milbrandt and Morrison (2013) for density of 550. for testing
          ax(lf) = 157.71
       ELSEIF ( icdx == 1 ) THEN
          bx(lf) = bxf
          ax(lf) = axf
       ELSEIF ( icdx > 1 ) THEN
          bx(lf) = 0.5
          ax(lf) = 75.7149
       ELSE ! icdx == 0
          ax(lf) = 206.984 ! Ferrier 1994 hail/frozen drops
          bx(lf) = 0.6384
        ENDIF
      ENDIF

      xnu(lc) = cnu
      xmu(lc) = 1.
      
      IF ( imurain == 3 ) THEN
        xnu(lr) = -0.8
        xmu(lr) = 1.
      ELSEIF ( imurain == 1 ) THEN
        xnu(lr) = (alphar - 2.0)/3.0
        xmu(lr) = 1./3.
      ELSE
         write(0,*) 'Value of imurain not valid. Must be 1 or 3, but is ',imurain
         call commasmpi_abort()
      ENDIF

      IF ( imusnow == 3 ) THEN
        xnu(ls) = -0.8
        xmu(ls) = 1.
      ELSEIF ( imusnow == 1 ) THEN
        xnu(ls) = (alphas - 2.0)/3.0
        xmu(ls) = 1./3.
      ELSE
         write(0,*) 'Value of imusnow not valid. Must be 1 or 3, but is ',imusnow
         call commasmpi_abort()
      ENDIF

      xnu(li) = cinu
      xmu(li) = 1.

      IF ( lis >= 1 ) THEN
      xnu(lis) = cinu ! 0.0
      xmu(lis) = 1.
      ENDIF

      dnu(lc) = 3.*xnu(lc) + 2. ! alphac
      dmu(lc) = 3.*xmu(lc)

      dnu(lr) = 3.*xnu(lr) + 2. ! alphar
      dmu(lr) = 3.*xmu(lr)
      
      dnu(ls) = 3.*xnu(ls) + 2. ! alphas
      dmu(ls) = 3.*xmu(ls)


      dnu(lh) = alphah
      dmu(lh) = dmuh

      xnu(lh) = (dnu(lh) - 2.)/3.
      xmu(lh) = dmuh/3.
      
      IF ( imurain == 3 ) THEN ! rain is gamma of volume
      rz =  ((4. + alphah)*(5. + alphah)*(6. + alphah)*(1. + xnu(lr)))/ & 
     &  ((1 + alphah)*(2 + alphah)*(3 + alphah)*(2. + xnu(lr)))

!      IF ( ipconc .lt. 5 ) alphahl = alphah
      
      rzhl =  ((4. + alphahl)*(5. + alphahl)*(6. + alphahl)*(1. + xnu(lr)))/ & 
     &  ((1. + alphahl)*(2. + alphahl)*(3. + alphahl)*(2. + xnu(lr)))

      rzs =  1. ! assume rain and snow are both gamma volume

      ELSE ! rain is gamma of diameter
      
      rz =  ((4. + alphah)*(5. + alphah)*(6. + alphah)*(1. + alphar)*(2. + alphar)*(3. + alphar))/ & 
     &  ((1 + alphah)*(2 + alphah)*(3 + alphah)*(4. + alphar)*(5. + alphar)*(6. + alphar))
      
      rzhl =  ((4. + alphahl)*(5. + alphahl)*(6. + alphahl)*(1. + alphar)*(2. + alphar)*(3. + alphar))/ & 
     &  ((1 + alphahl)*(2 + alphahl)*(3 + alphahl)*(4. + alphar)*(5. + alphar)*(6. + alphar))

      
      rzs =   & 
     &  ((1. + alphar)*(2. + alphar)*(3. + alphar)*(2. + xnu(ls)))/  &
     &  ((4. + alphar)*(5. + alphar)*(6. + alphar)*(1. + xnu(ls)))
       

      ENDIF
      
!      write(0,*) 'rz,rzhl = ', rz,rzhl
       
      IF ( ipconc .lt. 4 ) THEN

      dnu(ls) = alphas
      dmu(ls) = 1.

      xnu(ls) = (dnu(ls) - 2.)/3.
      xmu(ls) = 1./3.
      
      
      ENDIF
      
      IF ( lhl .gt. 1 ) THEN

      dnu(lhl) = alphahl
      dmu(lhl) = dmuhl

      xnu(lhl) = (dnu(lhl) - 2.)/3.
      xmu(lhl) = dmuhl/3.

      ENDIF

      IF ( lf .gt. 1 ) THEN

      dnu(lf) = alphah
      dmu(lf) = dmuh

      xnu(lf) = (dnu(lf) - 2.)/3.
      xmu(lf) = dmuh/3.

      ENDIF
      
      
      CALL setcnoz(cno)
!
!  density maximums and minimums
!
      xdnmx(:) = 900.0
      
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(ls) =  300.0
      xdnmx(lh) =  rho_qh_max ! 900.0
      IF ( lhl .gt. 1 ) xdnmx(lhl) = rho_qhl_max ! 900.0
      IF ( lf .gt. 1 ) xdnmx(lf) = rho_qf_max ! 900.0
!
      xdnmn(:) = 900.0
      
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lh) =  hdnmn
      IF ( lhl .gt. 1 ) xdnmn(lhl) = hldnmn
      IF ( lf .gt. 1 ) xdnmn(lf) = fdnmn

      xdn0(:) = 900.0
      
      xdn0(lc) = 1000.0
      xdn0(li) = 900.0
      xdn0(lr) = 1000.0
      xdn0(ls) = rho_qs ! 100.0
      xdn0(lh) = Min(rho_qh_max, rho_qh) ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      IF ( lhl .gt. 1 ) xdn0(lhl) = Min( rho_qhl_max, rho_qhl) ! 800.0
      IF ( lf .gt. 1 ) xdn0(lf) = Min(rho_qf_max, rho_qf) !  Min( rho_qhl_max, rho_qhl) ! 800.0
      
      IF ( lis > 1 ) THEN
        xdnmx(lis) = 917.0
        xdnmn(lis) = 100.0
        xdn0(lis) = 900.0
      ENDIF
      
      
         IF ( qhdpvdn < 0. ) qhdpvdn = xdnmn(lh)
         IF ( qhacidn < 0. ) qhacidn = xdnmn(lh)

!
!  Set terminal velocities...
!    also set drag coefficients
!
      cdx(lr) = 0.60
      cdx(lh) = 0.8 ! 1.0 ! 0.45
      cdx(ls) = 2.00
      cd(1) = cdx(ls)
      IF ( lhl .gt. 1 ) cdx(lhl) = 0.45
      IF ( lf .gt. 1 ) cdx(lf) = 0.8 ! 0.45  -- 0.8 for testing comparison of graupel and frozen drops


      
      dtp = dtp1 ! ntmul*dtp1
      IF ( my_rank == 0 ) THEN
        write(iunit,*) 'IN TDMPHY'
      ENDIF
!        write(0,*) 'IN TDMPHY: rank, time',my_rank,time_real
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), &
!! c$omp DEFAULT(SHARED),PRIVATE(kz,jy,ix)
!
! zero temporaries, fill t77 (full pi), t0 (temperature), t00 (pressure coeff in Teten's formula)
!
       call sett09s(nx,ny,nz,nor,norz,na,lt, & 
     &              t77,pinit,p2,an,pn,pb, & 
     &              t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00 )  


        IF ( ipconc .eq. 0 .and. lni .gt. 1 ) THEN
      
            an(1:nx,1:ny,1:nz,lni) = 0.0
          
        ENDIF

!
!  open input data file for the microphysics code
!
      if ( iread .eq. 0 ) then
!
      iread = 1
!

#ifdef MPI
!      write(0,*) 'INPUT:  ndebug (-1=no debug)'
!      IF ( my_rank == 0 ) write(0,*) 'ndebug= ',ndebug
#else
      IF ( my_rank == 0 ) write(iunit,*) 'INPUT:  ndebug (-1=no debug)'
      IF ( my_rank == 0 ) write(iunit,*) 'ndebug= ',ndebug
#endif
      
      IF ( my_rank == 0 ) write(iunit,*) 'ipconc = ',ipconc
!      IF ( ipconc .ne. 3 .and. ipconc .ne. 5 ) THEN
!        write(iunit,*) 'WARNING: Unsupported value for ipconc! Must be 3 or 5!'
!        write(iunit,*) 'Setting ipconc to 5.'
!        ipconc = 5
!      ENDIF

      IF ( ipconc .eq. 4 ) THEN
        IF ( my_rank == 0 ) write(iunit,*) 'WARNING: ipconc = 4 is unsupported!'
        IF ( my_rank == 0 ) write(iunit,*) 'Setting ipconc to 5.'
        ipconc = 5
      ENDIF

      IF ( my_rank == 0 ) THEN
        write(iunit,*) 'lr,lh,lnr,lnh,lzr,lzh = ',lr,lh,lnr,lnh,lzr,lzh
        write(iunit,*) 'lf,lnf,lzf = ',lf,lnf,lzf
        write(iunit,*) 'lhab = ',lhab
      ENDIF
      
      ido(lc) = idocw ; ido(lr) = idorw ; ido(li) = idoci
      ido(ls) = idosw
      ido(lh)  = idohw
      IF ( lhl .gt. 1 ) ido(lhl) = idohl
      IF ( lf .gt. 1 ) ido(lf) = idofw
      IF ( lis > 1 ) ido(lis) = 1

      IF ( my_rank == 0 ) write(iunit,'(a,6(1x,i2))') 'ido(lc:ls)',  &  
     & ido(lc), ido(lr), ido(li), ido(ls)

!
      IF ( my_rank == 0 ) write(iunit,'(a,6(1x,i2))') 'ido(lgl:lhl)', & 
     & ido(lh)

      iexy(:,:) = 0
      
      
#ifdef CHGELEC
      iexy(ls,li) = ieswi
      iexy(ls,lc) = ieswc ; iexy(ls,lr) = ieswr


      iexy(lh,ls) = iehwsw ; iexy(lh,li) = iehwi
      iexy(lh,lc) = iehwc  ; iexy(lh,lr) = iehwr

      IF ( lhl .gt. 1 ) THEN
      iexy(lhl,ls) = iehlsw ; iexy(lhl,li) = iehli
      iexy(lhl,lc) = iehlc  ; iexy(lhl,lr) = iehlr
      ENDIF

      IF ( lf .gt. 1 ) THEN
      iexy(lf,ls) = iefwsw ; iexy(lf,li) = iefwi
      iexy(lf,lc) = iefwc  ; iexy(lf,lr) = iefwr
      ENDIF
      
      IF ( lis > 1 ) THEN
        iexy(ls,lis) = iesis
        iexy(lh,lis) = iehis
        IF ( lhl > 1 ) iexy(lhl,lis) = iehlis
        IF ( lf > 1 ) iexy(lf,lis) = iefwis
      ENDIF

      IF ( isaund .le. -999 .or. nonigrd /= 0 ) THEN
      
      IF ( nonigrd .eq. 0 ) THEN
        isaund = Nint(rgard)
        IF ( my_rank == 0 ) write (iunit,*) 'RAR: isaund = ',isaund
      ELSEIF ( nonigrd .eq. 1 ) THEN
         isaund = Nint(rgard)
         IF ( isaund .ge. 1 ) THEN
          IF ( my_rank == 0 ) write(iunit,*) 'using qxacw-based RAR calculation'
         ELSE
         nonigrd = 0
         isaund = 0
          IF ( my_rank == 0 ) write (iunit,*) 'Saunders1991: isaund = ',isaund
         ENDIF
      ELSEIF ( nonigrd .eq. -1 ) THEN
        isaund = Nint(rgard)
        IF ( my_rank == 0 ) write (iunit,*) 'Takahashi: isaund = ',isaund
      ELSEIF ( nonigrd .eq. 2 ) THEN  ! Gardiner/Ziegler scheme
        rgard1 = rgard
        !  isaund = 11  ! not used for nonigrd=2 anymore. Use ftauopt instead
      ENDIF     
      
      ENDIF
      
      IF ( nonigrd .eq. 2 ) rgard1 = rgard

        IF ( my_rank == 0 ) write(iunit,*) 'nonigrd= ',nonigrd, '  rgard= ',rgard,  & 
     &   '  rgard1= ',rgard1, ' isaund = ',isaund, & 
     &   ' ewfac= ',ewfac, ' eii0 = ',eii0, ' eii1 = ',eii1

#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'ircnw, qminrncw=', ircnw, qminrncw
      write(iunit,*) 'cwdiap=', cwdiap
      write(iunit,*) 'cwdisp= ', cwdisp
      write(iunit,*) 'iauttim= ',iauttim
      write(iunit,*) 'auttim= ',auttim
      write(iunit,*) 'qcwmntim= ', qcwmntim
      write(iunit,*) 'rcond = ',rcond
      ENDIF
      
      IF ( fdfallfac < 0.0 ) THEN
         fdfallfac = graupelfallfac
      ENDIF
      linfall(:) = infall
      linfall(lc) = 0
      IF ( irfall .lt. 0 ) irfall = infall
      IF ( isfall .lt. 0 ) isfall = infall
      IF ( iifall .lt. 0 ) iifall = infall
      IF ( lzr > 0 ) irfall = 0
      IF ( lzs > 0 ) isfall = 0
      IF ( lzh > 0 ) linfall(lh) = 0
      IF ( lzhl > 0 .and. lhl > 0 ) linfall(lhl) = 0
      IF ( lzf > 0 .and. lf > 0 ) linfall(lf) = 0
      linfall(lr) = irfall
      linfall(ls) = isfall
      linfall(li) = iifall
      IF ( my_rank == 0 ) write(iunit,*) 'itfall, irfall,infall, iscfall,isfall,iifall = ', &
           itfall, irfall,infall,iscfall,isfall,iifall

      IF ( infall .ne. 1 .or. iscfall .ge. 2 ) THEN
         infdo = 1
      ELSE
         infdo = 0
      ENDIF

      IF ( Any(linfall(:) .ge. 3 ) .or. ipconc .ge. 6 .or. ldbzchanger >= 1 .or. ldbzchangeh >= 1 ) THEN
         infdo = 2
      ENDIF
      
      IF ( my_rank == 0 ) write(iunit,*) 'infdo = ',infdo
#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: DEBUG 0.1"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 0.1'
#endif

      cwccn = ccn

         pi = 4.0*atan(1.0)
!         xvcmx = (4./3.)*pi*xcradmx**3

        IF ( my_rank == 0 ) THEN
         write(iunit,*) 'iccwflg, irenuc = ', iccwflg, irenuc, irenuc3d
         write(iunit,*) 'xcradmx, xvcmx = ', xcradmx, xvcmx
         write(iunit,*) 'ciintmx = ',ciintmx
        ENDIF

#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: DEBUG 0.2"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 0.2'
#endif
        
!      cck = 0.6 ! this is now in micro_module
      cckm = cck-1.
      ccnefac =  (1.63/(cck * beta(3./2., cck/2.)))**(cck/(cck + 2.0))
      cnexp   = (3./2.)*cck/(cck+2.0)
!      ccn = cwccn
! ccne is all the factors with w in eq. A7 in Mansell et al. 2010 (JAS).  The constant changes
! if k (cck) is changed!
!      ccne = 0.9893*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
      ccne0 = ccnefac*1.e6*(1.e-6)**(2./(2.+cck))
      IF ( my_rank == 0 ) write(iunit,*) 'cwccn, cck, ccne, ccnefac, cnexp, ccne0 = ', &
              cwccn,cck,ccne,ccnefac,cnexp,ccne0


      IF ( cwccn .lt. 0.0 ) THEN
      cwccn = Abs(cwccn)
      ccwmx = 50.e9 ! cwccn
      ELSE
      ccwmx = 50.e9 ! cwccn ! *1.4
      ENDIF

      IF ( my_rank == 0 ) write(iunit,'(a,i2,2x,i2,1x,1pe12.5,4i3)')  & 
     &  'itype1, itype2,cimas0, icfn, ihrn, ibfc, iacr = ', & 
     &    itype1, itype2, cimas0, icfn, ihrn, ibfc, iacr



#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: DEBUG 1"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 1'
#endif

#ifdef CHGELEC
      IF ( nonigrd .eq. -1 .and. ipelec .gt. 0 ) THEN
       open(unit=90,file='takahashi.txt',form='formatted')
       DO ix=nlwc,0,-1
        read(90,*) (takalu(jy,ix), jy=0,31)
       END DO
       close(90)
      END IF
      
      IF ( ipelec > 0 ) THEN
        jchgs = Max(3, Int(0.5+charging_border/dy))
        jchgn = Max(2, Int(0.5+charging_border/dy))
        ichge = Max(3, Int(0.5+charging_border/dx))
        ichgw = Max(2, Int(0.5+charging_border/dx))
      ENDIF
#endif
      
      ireadqf = 0

!  Gamma function
!
!
!     compute fgoi(arg*100)
!
      if ( imkgam .eq. 0 ) then
      imkgam = 1
      call makegmoi()
      IF ( imurain == 1 ) THEN
      IF ( lhl > 1 ) THEN
        tmp = bx(lhl)
      ELSEIF ( lf > 1 ) THEN
        tmp = bx(lf)
      ELSE
        tmp = bx(lh)
      ENDIF
      CALL cld_cpu('MAKETABLE')  
       call MAKEqiacrratio(bx(lh),tmp)
      CALL cld_cpu('MAKETABLE')  
      ENDIF
      
!      do igam = 1,1000
!      arg = 0.01*igam
!      gmoi(igam) = gamma(arg)
!c      write(97,*) igam,gmoi(igam)
!      end do
      end if
!

!
! Build lookup table for saturation mixing ratio (Soong and Ogura 73)
!
!      cai = 21.87455
!      caw = 17.2693882
!      cbi = 7.66
!      cbw = 35.86

      do l = 1,nqsat
      temq = 163.15 + (l-1)*fqsat
      IF ( iqvsopt == 0 ) THEN
      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
      dtabqvs(l) = ((-caw*(-273.15 + temq))/(temq - cbw)**2 + & 
     &                 caw/(temq - cbw))*tabqvs(l)
      ELSE
      tabqvs(l) = exp(cawbolton*(temq-273.15)/(temq-cbwbolton))
      dtabqvs(l) = ((-cawbolton*(-273.15 + temq))/(temq - cbwbolton)**2 + & 
     &                 cawbolton/(temq - cbwbolton))*tabqvs(l)
      ENDIF
      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
      dtabqis(l) = ((-cai*(-273.15 + temq))/(temq - cbi)**2 + & 
     &                 cai/(temq - cbi))*tabqis(l)
      end do

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!        / T - 9.44523 * ALOG(T) + 0.014025 * T))

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)

!
!
!  Set collection coefficients (Seifert and Beheng 05)
!
      bb(:) = 1.0/3.0
      bb(li) = 0.3429
      IF ( lis >= 1 ) bb(lis) = 0.3429
      DO il = lc,lhab
        da0(il) = delbk(bb(il), xnu(il), xmu(il), 0)
        da1(il) = delbk(bb(il), xnu(il), xmu(il), 1)
        
!        write(0,*) 'il, da0, da1, xnu, xmu = ', il, da0(il), da1(il), xnu(il), xmu(il)
      ENDDO

      dab0(:,:) = 0.0
      dab1(:,:) = 0.0
      
      DO il = lc,lhab
        DO j = lc,lhab
          IF ( il .ne. j ) THEN
          
            dab0(il,j) = delabk(bb(il), bb(j), xnu(il), xnu(j), xmu(il), xmu(j), 0)
            dab1(il,j) = delabk(bb(il), bb(j), xnu(il), xnu(j), xmu(il), xmu(j), 1)
          
!           write(0,*) 'il, j, dab0, dab1 = ',il, j, dab0(il,j), dab1(il,j)
          ENDIF
        ENDDO
      ENDDO

      dab0lu(:,:,:,:) = 0.0
      dab1lu(:,:,:,:) = 0.0
      
      IF ( ipconc >= 6 ) THEN
      call cld_cpu('MAKETABLE')
      DO il = lc,lhab ! collector
        DO j = lc,lhab ! collected
          IF ( il .ne. j ) THEN

            DO jj = ialpstart,nqiacralpha
                alpjj = float(jj)*dqiacralpha
                xnujj = (alpjj - 2.)/3.
            DO ii = ialpstart,nqiacralpha
                alpii = float(ii)*dqiacralpha
                xnuii = (alpii - 2.)/3.
          
           
            dab0lu(ii,jj,il,j) = delabk(bb(il), bb(j), xnuii, xnujj, xmu(il), xmu(j), 0)
            dab1lu(ii,jj,il,j) = delabk(bb(il), bb(j), xnuii, xnujj, xmu(il), xmu(j), 1)
          
!            IF ( il == lr .and. j == lc .and. alpii == 0.0 .and. xnujj == 0.0 ) THEN
!              write(0,*) 'table: il,j,ii,jj = ',il,j,ii,jj
!              write(0,*) 'alpi,alpj,xnui,xnuj = ',alpii,alpjj,xnuii,xnujj
!              write(0,*) 'dab0lu = ',dab0lu(ii,jj,il,j),dab1lu(ii,jj,il,j) 
!            ENDIF
            ENDDO
            ENDDO
!           write(0,*) 'il, j, dab0, dab1 = ',il, j, dab0(il,j), dab1(il,j)
          ENDIF
        ENDDO
      ENDDO
      call cld_cpu('MAKETABLE')
      
      ENDIF
      
!      write(0,*) 'da0(lh): ',da0(lh),dab1(lh,li),da1(li)
!      write(0,*) 'dab0(lh): ',dab0(lh,li),da0(li)
!      write(0,*) 'da0(lr): ',da0(lr),dab1(lr,li),da1(li)
!      write(0,*) 'dab0(lr): ',dab0(lr,li),da0(li)
!      write(0,*) 'da0(lh,ls): ',da0(ls),dab1(lh,ls),da1(ls)
!      write(0,*) 'dab0(ls): ',dab0(lh,ls),da0(ls)
      
      end if    ! ( iread .eq. 0 ) 
!
!

!
!  electricity constants
!
!  mixing ratio epsilon
!
      qeps  = 1.0e-20

!  rebound efficiency (erbnd)
!
!
!      times = (nstep-1)*(dtp1*dtfac)
      times = (nstep)*(dtp1*dtfac)
!      headr1 = 'pa'
!
!  constants
!
      cp608 = 0.608          ! constant used in conversion of T to Tv
      cv = 717.0             ! specific heat at constant volume - air
      ar = 841.99666         ! rain terminal velocity power law coefficient (LFO)
      br = 0.8               ! rain terminal velocity power law coefficient (LFO)
!      aradcw = -0.27544      !
!      bradcw = 0.26249e+06   !
!      cradcw = -1.8896e+10   !
!      dradcw = 4.4626e+14    !
      bta1 = 0.6             ! beta-1 constant used for ice nucleation by deposition (Ferrier 94, among others)
      cnit = 1.0e-02         ! No for ice nucleation by deposition (Cotton et al. 86)
      dragh = 0.60           ! coefficient used to adjust fall speed for hail versus graupel (Pruppacher and Klett 78)
      dnz00 = 1.225          ! reference/MSL air density
      rho00 = 1.225          ! reference/MSL air density
      qccn = ccn/rho00
!      cs = 4.83607122       ! snow terminal velocity power law coefficient (LFO)
!      ds = 0.25             ! snow terminal velocity power law coefficient (LFO)
!  new values for  cs and ds
      cs = 12.42             ! snow terminal velocity power law coefficient 
      ds = 0.42              ! snow terminal velocity power law coefficient 
      pi = 3.141592653589793
      pii = 1./pi
      pid4 = pi/4.0

      gr = 9.8

!
!  constants
!
      c1f3 = 1.0/3.0
!
!  general constants for microphysics
!
!      cai = 21.87455         ! constant used for saturation mixing ratio wrt ice (Soong and Ogura 73)
!      caw = 17.2693882       ! constant used for saturation mixing ratio wrt water (Soong and Ogura 73)
!      cbi = 7.66             ! constant used for saturation mixing ratio wrt ice (Soong and Ogura 73)
!      cbw = 35.86            ! constant used for saturation mixing ratio wrt water (Soong and Ogura 73)

      f5 = 237.3 * 17.27 * 2.5e6 / cp ! combined constants for rain condensation (Soong and Ogura 73)

      call setqxminz(qxmin)
      
      
      xvmn(lc) = xvcmn
      xvmn(li) = xvimn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn ! 0.523599*(0.01e-3)**3
      xvmn(lh) = xvhmn

      xvmx(lc) = xvcmx
      xvmx(li) = xvimx
      xvmx(lr) = xvrmx
      xvmx(ls) = xvsmx
      xvmx(lh) = xvhmx
      
      IF ( lhl .gt. 1 ) THEN
        xvmn(lhl) = xvhlmn
        xvmx(lhl) = xvhlmx
      ENDIF

      IF ( lf .gt. 1 ) THEN
        xvmn(lf) = xvhmn ! xvfmn ! graupel values for comparison test
        xvmx(lf) = xvhmx ! xvfmx
      ENDIF
      
      IF ( lis > 1 ) THEN
        xvmn(lis) = xvimn
        xvmx(lis) = xvimx
      ENDIF


      tdtol = 1.0e-05
      thnuc = 233.15          ! homogeneous nucleation temperature threshold
      rw = 461.5              ! gas const. for water vapor
      advisc0 = 1.832e-05     ! reference dynamic viscosity (SMT; see Beard & Pruppacher 71)
      advisc1 = 1.718e-05     ! dynamic viscosity constant used in thermal conductivity calc
      tka0 = 2.43e-02         ! reference thermal conductivity
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
!
!
      if ( imake .eq. 0 ) then
      imake = 1
      
      IF ( lhl < 1 ) ifrzg = 1
      
      IF ( imltshddmr == -2 .or. imltshddmr == -3 ) THEN ! hack to be able to force imltshddmr = 3 for testing
        imltshddmr = Abs( imltshddmr )
      ELSEIF ( ipconc <= 5 .or. ( ibinhmlr == 0 .or. ibinhlmlr == 0 ) ) THEN 
        imltshddmr = Min(1, imltshddmr)
      ENDIF

      IF ( imurain == 3 ) THEN
!       IF ( izwisventr == 1 ) THEN
        ventr = Gamma(rnu + 4./3.)/((rnu + 1.)**(1./3.)*Gamma(rnu + 1.)) ! Ziegler 1985
!       ELSE
        ventrn =  Gamma(rnu + 1.5 + br/6.)/(Gamma(rnu + 1.)*(rnu + 1.)**((1.+br)/6. + 1./3.) ) ! adapted from Wisner et al. 1972
!        ventr = Gamma(rnu + 4./3.)/((rnu + 1.)**(1./3.)*Gamma(rnu + 1.)) ! Ziegler 1985
!        ventr  = Gamma(rnu + 4./3.)/Gamma(rnu + 1.) 
!       ENDIF
      ELSE ! imurain == 1
!       IF ( iferwisventr == 1 ) THEN
        ventr = Gamma(2. + alphar)  ! Ferrier 1994
!       ELSEIF ( iferwisventr == 2 ) THEN
        ventrn =  Gamma(alphar + 2.5 + br/2.)/Gamma(alphar + 1.) ! adapted from Wisner et al. 1972
!       ENDIF
      ENDIF
      ventc = Gamma(cnu + 4./3.)/((cnu + 1.)**(1./3.)*Gamma(cnu + 1.))
      c1sw = Gamma(snu + 4./3.)/((snu + 1.0)**(1./3.)*Gamma(snu + 1.0) )
!      print*,'ventr,ventc = ',ventr,ventc

!
!  Set up look up tables for supersaturation w.r.t. liq and ice
!
!CVD$L SKIP
!      do l = 1,nqsat
!      temq = 163.15 + (l-1)*fqsat
!      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
!      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
!      end do

#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: DEBUG 2"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 2'
#endif

!
!  cloud water constants in mks units
!
!      cwmasn = 4.25e-15  ! radius of 1.0e-6
!      cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
!      cwmasn5 =  5.23e-13
!      cwradn = 5.0e-6     ! minimum radius
!      cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
      mwfac = 6.0**(1./3.)
      IF ( ipconc .ge. 2 ) THEN
!        cwmasn = xvmn(lc)*1000.  ! minimum mass, defined by minimum droplet volume
!        cwradn = 1.0e-6          ! minimum radius
!        cwmasx = xvmx(lc)*1000.  ! maximum mass, defined by maximum droplet volume
        
        IF ( my_rank == 0 ) write(iunit,*) 'ICEZVD_DR: cwmasn,cwmasx = ',cwmasn,cwmasx
      ENDIF
        rwmasn = xvmn(lr)*1000.  ! minimum mass, defined by minimum rain volume
        rwmasx = xvmx(lr)*1000.  ! maximum mass, defined by maximum rain volume

      cwc1 = 6.0/(pi*1000.)
!
!  cloud ice constants in mks units
!
      cimasn = 6.88e-13  ! minimum mass
      cimasx = 1.0e-8    ! maximum mass
      ccimx = 5000.0e3   ! max of 5000 per liter

!
!  end if for constants (imake)
!
      end if


!
! Set up thermodynamic 1-D arrays
!
      kzb = 1
      kze = ktile+1
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag
      
      DO kz = kzb,kze
!      temp(kz) = pinit(kz)*ab(kz,lt)

!      tempc(kz) = temp(kz) - tfr
      
      ENDDO


      kzb = 1
      kze = ktile+1
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag
      
      IF ( .not. allocated ( cwnccnold ) ) THEN
        allocate( cwnccnold(nz+1) )
      ENDIF
      
      do kz = kzb,kze
!      cwdn(kz) = 1000.0
      
      IF ( ccn .gt. 0.0 ) THEN
!      cwnccn(kz) = cwccn*db(kz)/db(1)  
      cwnccnold(kz) = cwccn*db(kz)/rho00
      ELSE
      cwnccnold(kz) = cwccn
      ENDIF
 
      enddo
      

! 
!  constants for paramerization
!
      if ( ncuse .eq. 0 ) then
      kzb = 1
      kze = ktile+1
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
!      wvdf(kz) = (2.11e-05)*((temp(kz)/tfr)**1.94)* & 
!     &  (101325.0/(pb(kz)))                            ! diffusivity of water vapor, Hall and Pruppacher (76)
!      advisc(kz) = advisc0*(416.16/(temp(kz)+120.0))* & 
!     &  (temp(kz)/296.0)**(1.5)                         ! dynamic viscosity (SMT; see Beard & Pruppacher 71)
!      tka(kz) = tka0*advisc(kz)/advisc1                 ! thermal conductivity
      enddo
      end if
!
#ifdef MPI
      if (debug_mpi .and. my_rank>=0) write(0,*) my_rank, "ICEZVD_DR: DEBUG 3"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3'
#endif
!
!
!  end if for constants (imake)
!
!      end if
!
!
!  set 3-d temperature and saturation mixing ratio variables 
!
!
!
!  Stuff for Meyers et al. (1992) and Ferrier (1994) ice nucleation
!
!  constants (see Ferrier 94 eqns 4.31-4.34)
!
      cnin20 = 1.0e3
      cnin10 = 5.0e1
      cnin1a = 4.5
      cnin2a = 12.96
      cnin2b = 0.639
      ierr = 0

      kzb = 0
      kze = ktile+1
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ltemq)
!$omp PARALLEL DO  IF (nx*ny*nz .gt. 40*40*42), &
!$omp  DEFAULT(SHARED),  &
!$omp PRIVATE(kz,jy,ix,ltemq,a1,ssival,dp1,qvapor,ssifac,t8s,t9s)
      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
!  temperature
!
!      t0(ix,jy,kz) = an(ix,jy,kz,lt)
!     >  *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
      
      
!      t0(ix,jy,kz) = 
!     :     theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)))
!
!  look-up index (find nearest integer index of Soong & Ogura qvsat lookup table)
! 
      ltemq = Int( (t0(ix,jy,kz)-163.15)/fqsat+1.5 )
      IF ( ltemq .ge. nqsat .or. ltemq .le. 0 ) THEN
#ifndef MPI
!$omp critical
#endif
        ierr = ierr + 1
        IF ( ierr < 3 ) THEN
        write(iunit,*) 'WARNING 1: ltemq .ge. nqsat! ierr = ',ierr
        write(iunit,*) 'ix,jy,kz,t0,an,ltemq,t77,w = ',ix,jy,kz,t0(ix,jy,kz), & 
     &    an(ix,jy,kz,lt), ltemq,t77(ix,jy,kz),w(ix,jy,kz)
        write(0,*) 'WARNING 1: ltemq .ge. nqsat! ierr,my_rank = ',ierr,my_rank
        write(0,*) 'ix,jy,kz,t0,an,ltemq,w = ',ix,jy,kz,t0(ix,jy,kz), & 
     &    an(ix,jy,kz,lt), ltemq,w(ix,jy,kz)
        ENDIF
        ltemq = Min( nqsat, Max(1,ltemq) )
!        IF ( ltemq <= 0 ) THEN
!          t0(ix,jy,kz) = float(2)*float(fqsat) + 163.15
!          ltemq = Int( (t0(ix,jy,kz)-163.15)/fqsat+1.5 )
!        ELSE
!          t0(ix,jy,kz) = float(nqsat-2)*float(fqsat) + 163.15
!          ltemq = Int( (t0(ix,jy,kz)-163.15)/fqsat+1.5 )
!        ENDIF

      ! reset theta
!        an(ix,jy,kz,lt) = t0(ix,jy,kz)/t77(ix,jy,kz)

#ifndef MPI
!$omp end critical
#endif
      ENDIF

      IF ( ierr .gt. nx*ny ) THEN
      
       call commasmpi_abort()
       STOP
!$omp critical
       write(0,*) 'too many warnings! STOP, my_rank = ',my_rank
       write(iunit,*) 'too many warnings! STOP'
!$omp end critical
       EXIT
      ENDIF
!
! saturation mixing ratio
!
      a1 = t00(ix,jy,kz)      !pressure coeff in Teten's formula (380./(pn+pb))
      
      IF ( iqvsopt == 0 ) THEN
        t8s = a1*tabqvs(ltemq)  !saturation mixing ratio wrt water
      ELSE
        t8s = rdorv*esbolton*tabqvs(ltemq)/(pn(ix,jy,kz) + pb(kz) - esbolton*tabqvs(ltemq))
      ENDIF
      t9s = a1*tabqis(ltemq)  !saturation mixing ratio wrt ice

!      end do
!      end do
!      end do
      
!
!  calculate rate of nucleation
!

! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ssival)
! c$omp  PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), &
! c$omp  DEFAULT(SHARED),PRIVATE(kz,jy,ix,ssival,dp1,qvapor,ssifac)
!      do kz = 1,nz-1
!      do jy = 1,ny-jstag
!      do ix = 1,nx-istag
!
! filter supersaturation field (if issfilt = 1)
!
      ssival = Min(t8s,max(an(ix,jy,kz,lv),0.0))/t9s  ! qv/qvi
!
!      t7(ix,jy,kz) = 0.0
!
      if ( ssival .gt. 1.0 ) then
      
      IF ( icenucopt == 1 .or. icenucopt == 10 .or. icenucopt == -10 .or. icenucopt == -11 ) THEN ! Meyers/Ferrier
!
      if ( t0(ix,jy,kz).le.268.15 ) then
!     if ( t0(ix,jy,kz).le.273.15 ) then
        
       dp1 = dn(ix,jy,kz)/rho00*cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,jy,kz) = Min(dp1, 1.0d30)
      end if
      
!
!   have turned off nucleation by Meyer at warm temperatures
!  This is really from Ferrier (1994), eq. 4.31 - 4.34
      IF ( imeyers5 ) THEN
      if ( t0(ix,jy,kz).lt.tfr .and. t0(ix,jy,kz).gt.268.15 ) then
      qvapor = max(an(ix,jy,kz,lv),0.0) 
      ssifac = 0.0
      if ( (qvapor-t9s) .gt. 1.0e-5 ) then
      if ( (t8s-t9s) .gt. 1.0e-5 ) then
      ssifac = (qvapor-t9s) /(t8s-t9s)
      ssifac = ssifac**cnin1a   
      end if
      end if
      t7(ix,jy,kz) = dn(ix,jy,kz)/rho00*cnin10*ssifac*exp(-(t0(ix,jy,kz)-tfr)*bta1)
      end if
      ENDIF
      
      ELSEIF ( icenucopt == 2 ) THEN ! Thompson/Cooper; Note Thompson 2004 has constants of
                                     ! 0.005 and 0.304 because the line function was estimated from Cooper's plot
                                     ! Here, the fit line values from Cooper 1986 are converted. Very little difference 
                                     ! in practice
      
        t7(ix,jy,kz) = 1000.*0.00446684*exp(0.3108*(273.16 - Max(233.0, t0(ix,jy,kz) ) ) ) ! factor of 1000 to convert L**-1 to m**-3
      
      ELSEIF ( icenucopt == 3 ) THEN ! Phillips (Meyers/DeMott)

      if ( t0(ix,jy,kz).le.268.15 .and.  t0(ix,jy,kz) > 243.15 ) then ! Meyers with factor of Psi=0.06
        
       dp1 = 0.06*cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,jy,kz) = Min(dp1, 1.0d30)
      elseif ( t0(ix,jy,kz) <= 243.15 ) then ! Phillips estimate of DeMott et al (2003) data
       dp1 = 1000.*( exp( Min( 57.0 ,cnin2a*(ssival-1.1) ) ) )**0.3
       t7(ix,jy,kz) = Min(dp1, 1.0d30)
      
      end if

      ELSEIF ( icenucopt == 4 ) THEN ! DeMott 2010

        IF ( t0(ix,jy,kz) < 268.16 .and.  t0(ix,jy,kz) > 223.15 .and. ssival > 1.001 ) THEN ! 
      
        ! a = 0.0000594, b = 3.33, c = 0.0264, d = 0.0033,
        ! nint = a*(-Tc)**b * naer**(c*(-Tc) + d)
        ! nint has units of per (standard) liter, so mult by 1.e3 and scale by dn/rho00
        ! naer needs units of cm**-3, so mult by 1.e-6
        
        !  dp1 = 1.e3*0.0000594*(273.16 - t0(ix,jy,kz))**3.33 * (1.e-6*cin*dn(ix,jy,kz))**(0.0264*(273.16 - t0(ix,jy,kz)) + 0.0033)
          dp1 = 1.e3*dn(ix,jy,kz)/rho00*0.0000594*(273.16 - t0(ix,jy,kz))**3.33 * (1.e-6*naer)**(0.0264*(273.16 - t0(ix,jy,kz)) + 0.0033)
          t7(ix,jy,kz) = Min(dp1, 1.0d30)
      
        ELSE
          t7(ix,jy,kz) = 0.0
        ENDIF

      ELSEIF ( icenucopt == -1 ) THEN ! Fletcher (unlimited)
      
        t7(ix,jy,kz) = 1.0e-2*exp(0.6*abs(t0(ix,jy,kz)-273.0)) 

      ELSEIF ( icenucopt == -2 ) THEN ! Fletcher (limited at T = 246K as in Thompson et al. 2004)
      
        t7(ix,jy,kz) = 1.0e-2*exp(0.6*abs(Max(246.0,t0(ix,jy,kz))-273.0)) 
      
      ENDIF ! icenucopt
!
      end if ! ssival > 1.
!
      end do
      end do
      
!       write(iunit,*) 'ice_dr: k,t7 = ',kzbeg-1+kz,t7(nx/2,ny/2,kz),t0(nx/2,ny/2,kz), &
!     &         t00(nx/2,ny/2,kz),an(nx/2,ny/2,kz,lv),an(nx/2,ny/2,kz,li),an(nx/2,ny/2,kz,lc)
      
      end do

      IF ( ierr .gt. 4 ) THEN
!        call flush(iunit)
        STOP
      ENDIF
!
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), &
!! c$omp DEFAULT(SHARED),PRIVATE(kz,jy,ix)
!      DO kz=-nor+ng1,nz+nor
!      DO jy=-nor+ng1,ny+nor
!      DO ix=-nor+ng1,nx+nor
!        t8(ix,jy,kz) = 0.0
!        t9(ix,jy,kz) = 0.0
!      END DO
!      END DO
!      END DO
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
      iwetg1 = 0
      ilock4 = 0
      ilockc1 = 1
      ilockc2 = 2
     
      tkediss(:,:) = 0.0

       vtmaxall(:) = 0.0
       vtmaxn(:) = 0.0
       vtmaxz(:) = 0.0
       vtmaxq(:) = 0.0

! 
!***********************************************************
!  start jy loop 1
!***********************************************************
!
      db0(:,:) = 1.0
      DO kz = 1,nz+1
       dtz0(:,kz) = z1d(kz,3)
      ENDDO
      dbzchange(1:2,1:nz) = 0.0

         IF ( io_flag ) THEN
          IF ( lnvhl     >= 1 ) axtra(:,:,:,lnvhl   ) = 0.
          IF ( lzvhl     >= 1 ) axtra(:,:,:,lzvhl   ) = 0.
          IF ( lmnvhl    >= 1 ) axtra(:,:,:,lmnvhl  ) = 0.
          IF ( ialphahl  >= 1 ) axtra(:,:,:,ialphahl) = 0.
          IF ( ialphah   >= 1 ) axtra(:,:,:,ialphah ) = 0.
          IF ( ialphaf   >= 1 ) axtra(:,:,:,ialphaf ) = 0.
          IF ( ialphar   >= 1 ) axtra(:,:,:,ialphar ) = 0.
          IF ( idfw      >= 1 ) axtra(:,:,:,idfw    ) = 0.
          IF ( idhw      >= 1 ) axtra(:,:,:,idhw    ) = 0.
          IF ( idhl      >= 1 ) axtra(:,:,:,idhl    ) = 0.
          IF ( idmhl     >= 1 ) axtra(:,:,:,idmhl   ) = 0.
          IF ( idnhl     >= 1 ) axtra(:,:,:,idnhl   ) = 0.
          IF ( ld0       >= 1 ) axtra(:,:,:,ld0     ) = 0.
          IF ( lchlcnh   >= 1 ) axtra(:,:,:,lchlcnh ) = 0.
          IF ( lchlcnf   >= 1 ) axtra(:,:,:,lchlcnf ) = 0.
         ENDIF

      IF ( itfall .ne. 3 ) THEN

      CALL cld_cpu('SEDIMENTATION')  

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag


     ndfalltot = 0
      

     IF ( do1dsedimentation ) THEN ! {

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
     
! !$omp PARALLEL DO  
! !$omp DEFAULT(SHARED), &
! !$omp PRIVATE(kz,jy,ix,jgs,xvt,il,vtmax,ndfall,dtptmp,dtp,tmpn)
     do jy = jyb,jye

!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

      call sediment1d(dtp,nx,ny,nz,an,na,nor,norz,xfall,dn,z1d, &
     &                    t0,t7,axtra,io_flag,infdo,jy,ixe,kze,                       &
     &                    ldovol,rho00,cdx,cno,ido,ln,lz,lvol,lsc,ipc,lliq,lrain, &
     &                    xdn0,xvmn,xvmx,xdnmn,xdnmx,qxmin,cwnccnold, &
     &                    vtmaxz,vtmaxq,vtmaxn,vtmaxall,dbzchange, linfall ) 

!     subroutine sediment1d(dtp,nx,ny,nz,an,na,nor,norz,xfall,dn,dz3d,dz3dinv, &
!     &                    t0,t7,infdo,jslab,its,jts,    &

      enddo
     ELSE ! } {

!  Do fallout stuff now if itfall=0 (boxfall) or 2 (crowfall)
!
! !$omp PARALLEL DO  
! !$omp DEFAULT(SHARED), &
! !$omp PRIVATE(kz,jy,ix,jgs,xvt,il,vtmax,ndfall,dtptmp,dtp,tmpn)
     do jy = jyb,jye

!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

!      ix = 53 ; kz = 40
!      write(0,*) 'DR_0: qh,zh,ch = ',an(ix,jy,kz,lh),an(ix,jy,kz,lzh),an(ix,jy,kz,lnh)
      IF ( itfall .eq. 0 .or. itfall .eq. 1 ) THEN


!
!  zero the precip flux arrays (2d)
!

      xvt(:,:,:,:) = 0.0

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3a"
#else
      if ( ndebug .gt. 0 ) print*,'dbg = 3a'
#endif

      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO kz = kzb,kze
      DO ix = ixb,ixe
       db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      ENDDO

      DO kz = kzb,kze
      DO ix = ixb,ixe
       dtz1(ix,kz) = z1d(kz,3)/db1(ix,kz)
      ENDDO
      ENDDO

      IF ( lzh .gt. 1 ) THEN
      DO kz = kzb,kze
      DO ix = ixb,ixe
        an(ix,jy,kz,lzh) = Max( 0., an(ix,jy,kz,lzh) )
      ENDDO
      ENDDO
      ENDIF

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3a2"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3a2'
#endif

! set up vt arrays:
!      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
!     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)


      CALL cld_cpu('ZIEGFALL')  


      call ziegfall(nx,ny,nz,nor,norz,na,dtp,dz,jgs,ngs, & 
     &  xvt, & 
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, & 
     &  rwmasn,rwmasx,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccnold, & 
     &  itype1,itype2,infdo)

!      ix = 38
!      DO kz = 1,nz-1
!       write(iunit,*) 'kz,vt1/2/3 = ',kz,an(ix,jy,kz,lhl),an(ix,jy,kz,lnhl),xvt(ix,kz,1,lhl), &
!           xvt(ix,kz,2,lhl),xvt(ix,kz,3,lhl)
!      ENDDO
      
      CALL cld_cpu('ZIEGFALL')  
     

!       IF ( jy == 1 ) THEN
!        write(0,*) 'io_flag,lhl,lzvhl =',io_flag,lhl,lzvhl
!       ENDIF
       IF ( io_flag  ) THEN
      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
        vtmax = 0.0
        DO kz = kzb,kze
        DO ix = ixb,ixe
        
          IF ( lhl > 1 .and. (lzvhl > 1) ) THEN
          
          IF ( lmvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lmvhl) = xvt(ix,kz,1,lhl)
             vtmax = Max(vtmax,xvt(ix,kz,1,lhl))
          ENDIF
          IF ( lzvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lzvhl) = xvt(ix,kz,3,lhl)
          ENDIF
          IF ( lnvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lnvhl) = xvt(ix,kz,2,lhl)
             IF ( xvt(ix,kz,2,lhl) > 0. .and. lmnvhl > 0 ) THEN
               axtra(ix,jgs,kz,lmnvhl) = xvt(ix,kz,1,lhl) - xvt(ix,kz,2,lhl)
             ENDIF
          ENDIF
          
          ELSEIF ( lh > 1 .and. (lzvhl > 1) ) THEN
          ! if hail is turned off, then use graupel values

          IF ( lmvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lmvhl) = xvt(ix,kz,1,lh)
             vtmax = Max(vtmax,xvt(ix,kz,1,lh))
          ENDIF
          IF ( lzvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lzvhl) = xvt(ix,kz,3,lh)
          ENDIF
          IF ( lnvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lnvhl) = xvt(ix,kz,2,lh)
             IF ( xvt(ix,kz,2,lh) > 0. .and. lmnvhl > 0 ) THEN
               axtra(ix,jgs,kz,lmnvhl) = xvt(ix,kz,1,lh) - xvt(ix,kz,2,lh)
             ENDIF
          ENDIF
          
          ENDIF
          
         
         ENDDO
         ENDDO
!         write(0,*) 'vtmax = ',vtmax
       ENDIF

     DO il = lc,lhab

      vtmax = 0.0 ! this is actually Courant number

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
!       IF ( il == lhl .and. lzhl > 1 ) vtmaxzhl = Max( vtmaxzhl, xvt(ix,kz,3,il) ) ! track max reflectivity-wgt hail fall speed
       vtmaxz(il) = Max( vtmaxz(il), xvt(ix,kz,3,il) ) ! track max reflectivity-wgt fall speed
       vtmaxq(il) = Max( vtmaxq(il), xvt(ix,kz,1,il) ) ! track max mass-wgt fall speed
       vtmaxn(il) = Max( vtmaxn(il), xvt(ix,kz,2,il) ) ! track max number-wgt fall speed

      ! apply limit vtmaxsed (08/20/2015)
        xvt(ix,kz,1,il) = Min( vtmaxsed,  xvt(ix,kz,1,il) )
        xvt(ix,kz,2,il) = Min( vtmaxsed,  xvt(ix,kz,2,il) )
        xvt(ix,kz,3,il) = Min( vtmaxsed,  xvt(ix,kz,3,il) )

      vtmax = Max(vtmax,xvt(ix,kz,1,il)*gz(kz)) ! vt/dz
      vtmax = Max(vtmax,xvt(ix,kz,2,il)*gz(kz))
      vtmax = Max(vtmax,xvt(ix,kz,3,il)*gz(kz))
      
!        IF ( ny == 2 .and. an(ix,jy,kz,il) > qxmin(il) .and. ix == 40 ) THEN
!         write(0,*) 'DR: vt(',il,kz,')=',xvt(ix,kz,2,il),xvt(ix,kz,1,il),xvt(ix,kz,3,il)
!        ENDIF
      
!      ENDDO

       vtmaxall(il) = Max( vtmaxall(il), vtmax ) ! 

      ENDDO
      ENDDO


      
      IF ( dtp*vtmax .lt. 0.7 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*vtmax/0.7) + 1)
        ndfalltot = ndfalltot + 1
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

      call ziegfall(nx,ny,nz,nor,norz,na,dtptmp,dz,jgs,ngs, & 
     &  xvt, & 
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, & 
     &  rwmasn,rwmasx,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccnold, & 
     &  itype1,itype2,infdo)

    ! apply limit vtmaxsed (08/20/2015)
      do kz = kzb,kze
      do ix = ixb,ixe
        xvt(ix,kz,1,il) = Min( vtmaxsed,  xvt(ix,kz,1,il) )
        xvt(ix,kz,2,il) = Min( vtmaxsed,  xvt(ix,kz,2,il) )
        xvt(ix,kz,3,il) = Min( vtmaxsed,  xvt(ix,kz,3,il) )
      ENDDO
      ENDDO


      ENDIF ! (n .ge. 2)

!       IF ( il == 8 ) THEN
!         ix = 10
!         DO kz = kze,kze-10,-1
!           write(91,*) 'k,q,c,vt123 = ',kz,an(ix,jy,kz,il),an(ix,jy,kz,ln(il)),an(ix,jy,kz,lz(il))
!           write(91,*) '              ',xvt(ix,kz,1,il),xvt(ix,kz,2,il),xvt(ix,kz,3,il)
!         ENDDO
!       ENDIF

        IF ( il >= lr .and. ( linfall(il) .eq. 3 .or. linfall(il) .eq. 4 ) .and. ln(il) > 0 ) THEN
            call calczgr(nx,ny,nz,nor,na,an,ixe,kze, & 
     &         z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il), lvol(il), &
               xdn0(il), -1)
        ENDIF

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3b"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 1b'
#endif

! mixing ratio

!      DO il = lc,lhab
!      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &             an,db1,imapz,mzdist,il,1,xfall)
      call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &             an,db1,imapz,mzdist,il,1,xfall,dtz1)
!      ENDDO

!      DO il = ls,lhab
       IF ( il >= ls ) THEN
       IF ( lliq(il) .gt. 1 ) THEN
!       call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &             an,db1,imapz,mzdist,lliq(il),1,xfall)
       call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &             an,db1,imapz,mzdist,lliq(il),1,xfall,dtz1)
       ENDIF
       ENDIF
!      ENDDO

       IF ( il == lr ) THEN
       IF ( iraintypes >= 1 ) THEN
       DO m = 1,nraintypes
       call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,lr), & 
     &             an,db1,imapz,mzdist,lrain(m),1,xfall,dtz1)
       ENDDO
       ENDIF
       ENDIF


#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3c"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3c'
#endif

! volume

      IF ( ldovol .and. il >= li ) THEN
!      DO il = li,lhab
        IF ( lvol(il) .gt. 1 ) THEN
!         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &              an,db0,imapz,mzdist,lvol(il),0,xfall)
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &              an,db0,imapz,mzdist,lvol(il),0,xfall,dtz0)
        ENDIF
!      ENDDO
      ENDIF

! reflectivity

      IF ( ipconc .ge. 6 ) THEN
!      DO il = lr,lhab
        IF ( lz(il) .gt. 1 ) THEN
!         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), & 
!     &              an,db0,imapz,mzdist,lz(il),0,xfall)
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), & 
     &              an,db0,imapz,mzdist,lz(il),0,xfall,dtz0)
        ENDIF
!      ENDDO
      ENDIF

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3d"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3d'
#endif

! mean amount



       IF ( ipelec .ge. 2 ) THEN
      
!      DO il = lc,lhab
!        call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,iscfall,il), & 
!     &               an,db0,imapz,mzdist,lsc(il),0,xfall)
        call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,iscfall,il), & 
     &               an,db0,imapz,mzdist,lsc(il),0,xfall,dtz0)
!      ENDDO

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3e"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3e'
#endif
      
      
      ENDIF ! ipelec
      
      IF ( ipconc .gt. 0 ) THEN !{
!      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN

      IF ( linfall(il) .ge. 2  .and. lz(il) .lt. 1 ) THEN !{
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

        IF ( linfall(il) .eq. 3 .or. linfall(il) .eq. 4 .and. il >= lr ) THEN
        ! Method I or I+II
          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn2(ix,jy,kz) = z(ix,kz,il)
            ENDDO
          ENDDO
          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
        
        ELSE
          ! Method II
          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO

        ENDIF

      ELSE
       tmpn(:,jy,:) = 0.0
      ENDIF !}


#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3f"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3f'
#endif

       in = 2
       IF ( linfall(il) .eq. 1 ) in = 1

!         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), & 
!     &        an,db0,imapz,mzdist,ln(il),0,xfall)
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall,dtz0)


         IF ( lz(il) .lt. 1 ) THEN ! if not 3-moment, run one of the correction schemes
         IF ( linfall(il) >= 2 ) THEN
           
           xfall0(:,jgs) = 0.0

           IF ( linfall(il) .eq. 3 .or. linfall(il) .eq. 4 ) THEN
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), & 
!     &         tmpn2,db0,imapz,mzdist,1,0,xfall0)           
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &         tmpn,db0,imapz,mzdist,1,0,xfall0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,3,il), & 
     &         tmpn2,db0,imapz,mzdist,1,0,xfall0,dtz0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &         tmpn,db0,imapz,mzdist,1,0,xfall0,dtz0)
           ELSE
!             call boxfall(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
!     &         tmpn,db0,imapz,mzdist,1,0,xfall0)
             call fallout(nx,ny,nz,nor,1,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &         tmpn,db0,imapz,mzdist,1,0,xfall0,dtz0)
           ENDIF

           IF ( linfall(il) .eq. 3 .or. linfall(il) .eq. 4 ) THEN
! "Method I" - dbz correction
             kze = ktile+1
             if (kzend .eq. nzend) kze = kzend-kzbeg

             ixe = itile
             if (ixend .eq. nxend) ixe = ixend-ixbeg

             call calcnfromz(nx,ny,nz,nor,na,an,tmpn2,ixe,kze, & 
     &       z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il),tmpn,  & 
     &       lvol(il), xdn0(il), infall , -1)

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
           ELSEIF ( linfall(il) .eq. 5 ) THEN

             kze = ktile+1
             if (kzend .eq. nzend) kze = kzend-kzbeg

             ixe = itile
             if (ixend .eq. nxend) ixe = ixend-ixbeg


             DO kz = kzb,kze
              DO ix = ixb,ixe
               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), 0.5* ( an(ix,jgs,kz,ln(il)) + tmpn(ix,jy,kz) ))
              
              ENDDO
             ENDDO           
           ELSEIF ( .not. (il .eq. lr .and. irfall .eq. 0) .and.  &
                    .not. (il .eq. ls .and. isfall .eq. 0) .and.  &
                    .not. (il .eq. li .and. iifall .eq. 0) ) THEN
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
      
      ! If have 3M hail and frozen drops, then do number of hail from frozen drops (if lnhlf > 1)
       IF ( lhl > 1 .and. il == lhl .and. lnhlf > 1 .and. lzhl > 1 ) THEN
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), & 
     &        an,db0,imapz,mzdist,lnhlf,0,xfall,dtz0)
       
       ENDIF

       IF ( il == lh .and. lnhf > 1 .and. lzh > 1 ) THEN
         call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,in,il), & 
     &        an,db0,imapz,mzdist,lnhf,0,xfall,dtz0)
       
       ENDIF

      ENDDO ! n=1,ndfall
      ENDDO ! il
      


#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3g"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3g'
#endif
      
      ELSEIF ( itfall .eq. 2 ) THEN
      
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

      call ziegfall(nx,ny,nz,nor,norz,na,dtp,dz,jgs,ngs, & 
     &  xvt, & 
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, & 
     &  rwmasn,rwmasx,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccnold, & 
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3h"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3h'
#endif

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
      DO il = lr,lhab
      vtmax = Max(vtmax,xvt(ix,kz,1,il)*gz(kz))
      vtmax = Max(vtmax,xvt(ix,kz,2,il)*gz(kz))
      vtmax = Max(vtmax,xvt(ix,kz,3,il)*gz(kz))
      ENDDO

      ENDDO
      ENDDO


      
      IF ( dtp*vtmax .lt. 0.9 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*vtmax) + 1)
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp/Real(ndfall)
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

      call ziegfall(nx,ny,nz,nor,norz,na,dtptmp,dz,jgs,ngs, & 
     &  xvt, & 
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, & 
     &  rwmasn,rwmasx,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccnold, & 
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)

      ENDIF ! (n .ge. 2)

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3i"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3i'
#endif

! mixing ratio

      DO il = lc,lhab
      call crowfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &  an,db1,imapz,mzdist,il,1,xfall)
      ENDDO

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3j"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3j'
#endif

! mean amount


       IF ( ipelec .ge. 2 ) THEN
      
!      iscfall = 1
      DO il = lc,lhab
        call crowfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,iscfall,il), & 
     &    an,db0,imapz,mzdist,lsc(il),0,xfall)
      ENDDO
      
#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3k"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3k'
#endif

      ENDIF ! ipelec
      
!      infall = 2
!      DO il = lc,lhab
!        IF ( ipconc .ge. ipc(il) ) THEN
!         call crowfall(nx,ny,nz,nor,na,gt,ngt,dtptmp,dz,jgs,xvt(1,1,infall,il),
!     :        an,db0,imapz,mzdist,ln(il),0,xfall)
!        ENDIF
!      ENDDO
      IF ( ipconc .gt. 0 ) THEN
      DO il = lc,lhab

      IF ( ipconc .ge. ipc(il) ) THEN
      IF ( linfall(il) .eq. 2 ) THEN
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!

          kzb = 1
          kze = ktile+1
          if (kzend .eq. nzend) kze = kzend-kzbeg+1

          ixb = 1
          ixe = itile
          if (ixend .eq. nxend) ixe = ixend-ixbeg+1

          DO kz = kzb,kze
            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
      ENDIF

         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,linfall(il),il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
!         call boxfall(nx,ny,nz,nor,na,gt,ngt,dtptmp,dz,jgs,xvt(1,1,infall,il),
!     :        an,db0,imapz,mzdist,ln(il),0,xfall)

         IF ( linfall(il) .eq. 2 .and. & 
     &       ( il .eq. lr .or. (il .ge. ls .and. il .le. lhab) )) THEN
!     :        .or. il .eq. lhl )) THEN

         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
           
!           call boxfall(nx,ny,nz,nor,1,gt,ngt,dtptmp,dz,jgs,xvt(1,1,1,il),
!     :       tmpn,db0,imapz,mzdist,1,0,xfall0)

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICE_ZIEG: DEBUG 3l"
#else
      if (ndebug .gt. 0 ) print*,'dbg = 3l'
#endif

           kzb=1
           kze=ktile
           if(kzend.eq.nzend) kze=kzend-kzbeg

           ixb=1
           ixe=itile
           if(ixend.eq.nxend) ixe=ixend-ixbeg
      
           DO kz = kzb,kze
           DO ix = ixb,ixe
             an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
             
           ENDDO
           ENDDO
           
!          ENDDO

         ENDIF
        ENDIF
      ENDDO
      
      ENDIF


      ENDDO ! n=1,ndfall
      

      ENDIF ! itfall.eq.2

      ENDDO ! jy loop over fallout
      
      ENDIF ! } ! do1dsedimentation

#ifdef CHGELEC
      IF ( ipelec .ge. 1 ) THEN
      chgneg = 0.0d0
      chgpos = 0.0d0
! C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix), &
!! c$omp REDUCTION (+ : chgpos,chgneg)

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-1*kd1

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-1*jd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-1*id1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
      sctot =    &
     & ( an(ix,jy,kz,lscw)    &
     &  + an(ix,jy,kz,lscr)    &
     &  + an(ix,jy,kz,lsci)    &
     &  + an(ix,jy,kz,lscs)   &
     &  + an(ix,jy,kz,lsch)   &
     &  + ec*(an(ix,jy,kz,lscpi)  - an(ix,jy,kz,lscni)) )
       IF ( largeion ) THEN
       sctot = sctot     &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) 
       ENDIF
      IF ( lhl .gt. 1 ) sctot = sctot + an(ix,jy,kz,lschl)
      IF ( lf .gt. 1 ) sctot = sctot + an(ix,jy,kz,lscf)
      IF ( lis .gt. 1 ) sctot = sctot + an(ix,jy,kz,lscis)
      
      dv = dxx(ix)*dyy(jy)*dzz(kz)
      IF ( sctot .gt. 0 ) THEN
        chgpos = chgpos + sctot*dv
      ELSE
        chgneg = chgneg + sctot*dv
      END IF

      ! temporary advection+mixing tendency, assuming iscnet has value from start of time step
      IF ( ichgtndsed > 1 ) THEN
        elec(ichgtndsed)%flt3d(ix,jy,kz) = elec(ichgtndsed)%flt3d(ix,jy,kz) + (sctot - elec(iscnet)%flt3d(ix,jy,kz) )
        elec(iscnet)%flt3d(ix,jy,kz) = sctot
      ENDIF

      
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
      write(iunit,'(a,3(2x,1pe15.8))') 'post-sed pos/neg/net charge (C):',   &
     &      chgpos,chgneg,chgpos+chgneg
      
      ENDIF 
      ENDIF
#endif /* chgelec */

      CALL cld_cpu('SEDIMENTATION')  
      
      
      ENDIF ! itfall .ne. 3

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICEZVD_DR: DEBUG 4"
#else
      if (ndebug .gt. 0 ) print*,'ICEZVD_DR: dbg = 4'
#endif
! 
!***********************************************************
!  start jy loop 2
!***********************************************************


       IF ( ldbzchangeh > 1 ) THEN
         DO kz = 1,nz-1
           thproc(kz,numproc) = thproc(kz,numproc) + dbzchange(2,kz)
         ENDDO
       ENDIF
       IF ( ldbzchanger > 1 ) THEN
         DO kz = 1,nz-1
           thproc(kz,numproc-1) = thproc(kz,numproc-1) + dbzchange(1,kz)
         ENDDO
       ENDIF

      IF ( isedonly == 1 ) GOTO 39999 ! skips microphysics and nucleation
      
      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye=jyend-jybeg+1-jstag


          IF ( lmvr  >= 1 ) axtra(:,:,:,lmvr ) = 0.
          IF ( lmvh  >= 1 ) axtra(:,:,:,lmvh ) = 0.
          IF ( lmvhl >= 1 .and. lhl > 1 .and. lzvhl == 0 ) axtra(:,:,:,lmvhl) = 0.
          IF ( lmvs  >= 1 ) axtra(:,:,:,lmvs ) = 0.
          IF ( lmvi  >= 1 ) axtra(:,:,:,lmvi ) = 0.
          IF ( ld0   >= 1 ) axtra(:,:,:,ld0  ) = 0.
          IF ( lqaut   >= 1 ) axtra(:,:,:,lqaut  ) = 0.
          IF ( lqacc   >= 1 ) axtra(:,:,:,lqacc  ) = 0.
          IF ( lqrevap >= 1 ) axtra(:,:,:,lqrevap) = 0.
          IF ( lqcevap >= 1 ) axtra(:,:,:,lqcevap) = 0.
          IF ( lqcond  >= 1 ) axtra(:,:,:,lqcond ) = 0.
          IF ( lqdep   >= 1 ) axtra(:,:,:,lqdep  ) = 0.
          IF ( lqmelt  >= 1 ) axtra(:,:,:,lqmelt ) = 0.
          IF ( lqsub   >= 1 ) axtra(:,:,:,lqsub  ) = 0.
          IF ( lswdia  >= 1 ) axtra(:,:,:,lswdia ) = 0.
          IF ( lswdn   >= 1 ) axtra(:,:,:,lswdn  ) = 0.
          IF ( lfwdn   >= 1 ) axtra(:,:,:,lfwdn  ) = 0.
          IF ( lhwdn   >= 1 ) axtra(:,:,:,lhwdn  ) = 0.
          IF ( lhldn   >= 1 ) axtra(:,:,:,lhldn  ) = 0.
          IF ( lswagg  >= 1 ) axtra(:,:,:,lswagg ) = 0.
          IF ( ldhlcnh >= 1 ) axtra(:,:,:,ldhlcnh) = 0.
          IF ( lchlcnf >= 1 ) axtra(:,:,:,ldhlcnf) = 0.
!          IF ( lchlcnh >= 1 ) axtra(:,:,:,lchlcnh) = 0.
!          IF ( lchlcnf >= 1 ) axtra(:,:,:,lchlcnf) = 0.

!          IF ( lnvhl     >= 1 ) axtra(:,:,:,lnvhl   ) = 0.
!          IF ( lzvhl     >= 1 ) axtra(:,:,:,lzvhl   ) = 0.
!          IF ( lmnvhl    >= 1 ) axtra(:,:,:,lmnvhl  ) = 0.
!          IF ( ialphahl  >= 1 ) axtra(:,:,:,ialphahl) = 0.
!          IF ( idhl      >= 1 ) axtra(:,:,:,idhl    ) = 0.
!          IF ( idmhl     >= 1 ) axtra(:,:,:,idmhl   ) = 0.
!          IF ( idnhl     >= 1 ) axtra(:,:,:,idnhl   ) = 0.



! !$omp PARALLEL DO  &
! !$omp  DEFAULT(SHARED), &
! !$omp PRIVATE(jy,jgs,alpha2d)
      DO jy = jyb,jye

!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

! things to pass:

! an,t00,tabqvs,tabqis, dtabqvs,dtabqis,gmoi,
!  xdnmx,xdnmn,xdn0,dnu,dmu,xnu,xmu,lsc,ln,ipc,lvol,bb,dab0,dab1,da0,da1
!  ventr,ventc,c1sw
!  cdx,cd,dtp
!  takalu,ido
!  tq1
!
       CALL cld_cpu('ICEZVD_GS')  

      IF ( iturbenhance > 0 ) THEN
      
      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
       
       DO kz = kzb,kze
         DO ix = ixb,ixe
        
          IF ( len_type .eq. 0 .or. len_type .eq. 3 ) dyl = (dxx(ix)*dyy(jy)*dzz(kz))**(1./3.)
        
          IF( len_type .eq. 0 ) mlen = Cm*dyl ! (dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dzz(kz)
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(kz)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
          IF( len_type .eq. 4 ) THEN
            IF ( gz(kz) .gt. 1400. ) THEN
             mlen = Cm*dyl
            ELSE
             mlen = Cm / (1./(0.4*gz(kz)) + 1./dyl)
            ENDIF
          ENDIF
          tkediss(ix,kz) =  0.5*Ce*Cm*(km(ix,jy,kz)/mlen)**2 
      
        ENDDO
       ENDDO
      
      ENDIF

!      ix = 53 ; kz = 40
!      write(0,*) 'DR_1: qh,zh,ch = ',an(ix,jy,kz,lh),an(ix,jy,kz,lzh),an(ix,jy,kz,lnh)

      IF ( ihailmax2d > 0 ) THEN
         alpha2d(ixb:ixe,1,kzb:kze,1) = alphar
         alpha2d(ixb:ixe,1,kzb:kze,2) = alphah
         alpha2d(ixb:ixe,1,kzb:kze,3) = alphahl
         IF ( lf > 1 )  alpha2d(ixb:ixe,1,kzb:kze,4) = alphah
      ENDIF


! isedonly=4 skips rates and just does axtra fields

        call nssl_2mom_gs   &
     &  (nx,ny,nz,na,jybeg+jy  &
     &  ,nor,norz          &
     &  ,dtp,gz       &
     &  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
     &  ,an,dn,p2                  &
     &  ,pn,w,iunit                   &
     &  ,t00,t77,                             &
     &   ventr,ventc,c1sw,jgs,ido,    &
     &   xdnmx,xdnmn,               &
!     &   ln,ipc,lvol,lz,lliq,   &
     &   cdx,                              &
     &   xdn0,tmpn,tkediss  &
     &  ,dxx,dyy,dzz,z1d,io_flag   &
     &  ,ventrn,qxmin, xvmn, xvmx, cno  &
     &  ,cckm,ccne,ccnefac,cnexp,ccne0  &
     &  ,lsc,ln,ipc,lvol,lz,lliq,lrain &
     &  ,dab0,dab1,da0,da1,bb &
     &  ,dab0lu,dab1lu &
     &  ,iexy,takalu &
     &  ,ab,pb, pinit, time_real,   &
     &   ptotalmxy(jy),ptotalmny(jy), & 
     &   psctotmxy(jy),psctotmny(jy),iwetgy(jy),iwetg1y(jy),iptotaly(jy), & 
     &   nscmaxy(jy),nscminy(jy),psctot1y(jy),psctot2y(jy), & 
     &   sctot1ny(jy),sctot1py(jy),            & 
     &   ngscmaxy(lg:lhab,lc:ls,jy),ngscminy(lg:lhab,lc:ls,jy),ngscxymax(lg:lhab,lc:ls,jy), &
     &   thproc,numproc,    &
     &   rate2d,nrate2d,alpha2d,nalpha2d,    &
     &  axtra,tq1,ntq1,   &
     &   ni,nj,nk            &
     & ,elec,ixbeg,nxbeg,nxend,nybeg,nyend  & ! its,ids,ide,jds,jde &
     & )
!       write(0,*) 'done jy = ',jy,Maxval( t2 )

!       IF ( lf > 1 ) THEN
!       
!       DO kz = kzb,kze
!         DO ix = ixb,ixe
!           an(ix,jy,kz,lf) = Max(0.0, an(ix,jy,kz,lf) )
!         ENDDO
!       ENDDO
!       
!       ENDIF

         IF ( ihailmax2d > 0 .and. ihailmaxk1 > 0) THEN
           call hailmaxd(dtp,nx,ny,nz,an,na,nor,norz,alpha2d,dn,qxmin,   &
                xfalltot(-nor+ng1,-nor+ng1,ihailmax2d),xfalltot(-nor+ng1,-nor+ng1,ihailmaxk1), &
                jy,nalpha2d,ixb,ixe,kzb,kze )
         ENDIF


       CALL cld_cpu('ICEZVD_GS')  

       IF ( .not. ( isedonly == 2 .or. isedonly == 4 ) ) THEN 
       CALL NUCOND    &
     &  (nx,ny,nz,na,jy & 
     &  ,nor,norz,dtp,ni &
     &  ,dzz & 
     &  ,t0,t9 & 
     &  ,an,dn,p2 & 
     &  ,pn,w & 
#ifdef COMMAS
     &  ,ventr,ventrn,ventc,qxmin, xdn0, xvmn, xvmx, cno  &
     &  ,cckm,ccne,ccnefac,cnexp,ccne0  &
     &  ,ido,lsc,ln,ipc,lvol,lz,lliq,lrain  &
     &  ,xdnmn,xdnmx &
     &  ,pb, pinit, thproc, numproc, dxx,dyy,dzz &
     &  ,axtra, io_flag &
#endif
     &  ,ssfilt,t00,t77,.false.)


!mpi: need ssfilt(i+1,i-1,j+1...) for ss gradient
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg


      do kz = kzb,kze
        do ix = ixb,ixe

         temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz)
        
         ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )
          IF ( iqvsopt == 0 ) THEN
            c1 = t00(ix,jy,kz)*tabqvs(ltemq)
          ELSE
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pn(ix,jy,kz) + pb(kz) - esbolton*tabqvs(ltemq))
          ENDIF

          ssat(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values

        ENDDO
      ENDDO
        
        ENDIF

! Clean up tiny values of mixing ratio
       IF ( .true. ) THEN
       CALL smallvalues    &
     &  (nx,ny,nz,na,jy & 
     &  ,nor,norz,dtp,ni &
     &  ,t0 & 
     &  ,an,dn,w & 
#ifdef COMMAS
     &  ,qxmin, xdn0, xvmn, xvmx, cno  &
     &  ,ido,lsc,ln,ipc,lvol,lz,lliq,lrain  &
     &  ,xdnmn,xdnmx &
#endif
     &  ,t77,.false.)
       ENDIF

      ENDDO ! jy


! For droplet cond/nuc:
!   wvdf,advisc, tka
!
! arrays for psctotmx,psctotmn,iwetg,iwetg1,iptotal,nscmax,nscmin
!      real psctotmxy(ny),psctotmny(ny)
!      integer iwetgy(ny),iwetg1y(ny),iptotaly(ny),nscmaxy(ny),nscminy(ny)
!
!..Gather microphysics  
!
!
!  CALL icezvd_gs
!
!
!  precipitation fallout contributions (if itfall .eq. 3 ) then
!
!
!  Do fallout stuff now if itfall=3
!
      IF ( .false. .and. itfall .eq. 3 ) THEN

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "ICEZVD_DR: DEBUG 5"
#else
      if (ndebug .gt. 0 ) print*,'ICEZVD_DR: dbg = 5'
#endif

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      DO jy = jyb,jye
!
!  zero the precip flux arrays (2d)
!
      xvt(:,:,:,:) = 0.0


      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
      
      do kz = kzb,kze
      do ix = ixb,ixe
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

! set up vt arrays:
!      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
!     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)

      call ziegfall(nx,ny,nz,nor,norz,na,dtp,dz,jgs,ngs, & 
     &  xvt, & 
     &  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx, & 
     &  rwmasn,rwmasx,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccnold, & 
     &  itype1,itype2,infdo)
!     :  rwmasn,rwmasx)
      
! mixing ratio      
      DO il = lc,lhab
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il), & 
     &  an,db1,imapz,mzdist,il,1,xfall)
      ENDDO
      




       IF ( ipelec .ge. 2 ) THEN
      
!      iscfall = 1
      DO il = lc,lhab
        call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,iscfall,il), & 
     &  an,db0,imapz,mzdist,lsc(il),0,xfall)
      ENDDO
      
      
      ENDIF ! ipelec
      
      IF ( ipconc .ge. 1 ) THEN
!      infall = 2
      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN
         IF ( .not. (lnr .ge. 1 .and. ln(il) .eq. lnr ) ) THEN
         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,linfall(il),il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
         ELSE
           IF ( itfall .eq. 1 ) THEN ! special test case to use number-weighted regardless of infall
            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,2,il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
           ELSE
            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,linfall(il),il), & 
     &        an,db0,imapz,mzdist,ln(il),0,xfall)
           ENDIF
         ENDIF
        ENDIF
      ENDDO
      ENDIF

      
      ENDDO
      
      ENDIF ! itfall.eq.3
!


39999  continue


!
!  Intermediate charge totals check
!
!  
!  Count up total net positive (chgpos) and net negative charge (chgneg)
!
!#ifdef ELEC
#ifdef CHGELEC
      IF ( ipelec .ge. 1 ) THEN
      chgneg = 0.0d0
      chgpos = 0.0d0
! C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix), &
!! c$omp REDUCTION (+ : chgpos,chgneg)

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-1*kd1

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-1*jd1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-1*id1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
      sctot =    &
     & ( an(ix,jy,kz,lscw)    &
     &  + an(ix,jy,kz,lscr)    &
     &  + an(ix,jy,kz,lsci)    &
     &  + an(ix,jy,kz,lscs)   &
     &  + an(ix,jy,kz,lsch)   &
     &  + ec*(an(ix,jy,kz,lscpi)  - an(ix,jy,kz,lscni)) )
       IF ( largeion ) THEN
       sctot = sctot     &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) 
       ENDIF
      IF ( lhl .gt. 1 ) sctot = sctot + an(ix,jy,kz,lschl)
      IF ( lf .gt. 1 ) sctot = sctot + an(ix,jy,kz,lscf)
      IF ( lis .gt. 1 ) sctot = sctot + an(ix,jy,kz,lscis)
      
      dv = dxx(ix)*dyy(jy)*dzz(kz)
      IF ( sctot .gt. 0 ) THEN
        chgpos = chgpos + sctot*dv
      ELSE
        chgneg = chgneg + sctot*dv
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
      write(iunit,'(a,3(2x,1pe15.8))') 'post-jy pos/neg/net charge (C):',   &
     &      chgpos,chgneg,chgpos+chgneg
      
      ENDIF 
      ENDIF

!#endif
!
      jyb=1
      jye=jtile
      if(jyend.eq.nyend) jye=jyend-jybeg

      DO jy = jyb,jye
        ptotalmx = Max( ptotalmx, ptotalmxy(jy) )
        ptotalmn = Min( ptotalmn, ptotalmny(jy) )
        psctotmx = Max( psctotmx, psctotmxy(jy) )
        psctotmn = Min( psctotmn, psctotmny(jy) )
        iwetg = iwetg + iwetgy(jy)
        iwetg1 = iwetg1 + iwetg1y(jy)
        iptotal = iptotal + iptotaly(jy)
        nscmax = nscmax + nscmaxy(jy)
        nscmin = nscmin + nscminy(jy)
        psctot1 = psctot1 + psctot1y(jy)
        psctot2 = psctot2 + psctot2y(jy)
        sctot1n = sctot1n + sctot1ny(jy)
        sctot1p = sctot1p + sctot1py(jy)
      ENDDO

#endif /* chgelec */

#ifdef MPI

       mpitotin(1) = ptotalmx
       mpitotin(2) = psctotmx
       
       k = 2
       
       CALL MPI_Reduce(mpitotin, mpitotout, k, MPI_REAL, MPI_MAX, 0, my_comm, mpi_error_code)
       
       IF ( my_rank == 0 ) THEN
         ptotalmx = mpitotout(1)
         psctotmx = mpitotout(2)
       ENDIF

       mpitotin(1) = ptotalmn
       mpitotin(2) = psctotmn
       
       k = 2
       
       CALL MPI_Reduce(mpitotin, mpitotout, k, MPI_REAL, MPI_MIN, 0, my_comm, mpi_error_code)
       
       IF ( my_rank == 0 ) THEN
         ptotalmn = mpitotout(1)
         psctotmn = mpitotout(2)
       ENDIF

! find global integrated rate max
       mpitotindp(1)  = psctot1
       mpitotindp(2)  = psctot2
       mpitotindp(3)  = sctot1p
       mpitotindp(4)  = sctot1n
       mpitotindp(5)  = iptotal
       mpitotindp(6)  = iwetg
       mpitotindp(7)  = iwetg1
       mpitotindp(8)  = ndfalltot
       
       k = 8
       
#ifdef CHGELEC
       IF ( ipelec .gt. 0 ) THEN
       mpitotindp(k+1) = Sum(ngscmaxy(lh,li,1:jye))
       mpitotindp(k+2) = Sum(ngscminy(lh,li,1:jye))
       mpitotindp(k+3) = Sum(ngscmaxy(lh,ls,1:jye))
       mpitotindp(k+4) = Sum(ngscminy(lh,ls,1:jye))
       mpitotindp(k+5) = Sum(ngscxymax(lh,li,1:jye))
       
       k = k + 5
       
       IF ( lhl .gt. 1 ) THEN
       
       mpitotindp(k+1)  = Sum(ngscmaxy(lhl,li,1:jye))
       mpitotindp(k+2)  = Sum(ngscminy(lhl,li,1:jye))
       mpitotindp(k+3)  = Sum(ngscmaxy(lhl,ls,1:jye))
       mpitotindp(k+4)  = Sum(ngscminy(lhl,ls,1:jye))
       mpitotindp(k+5)  = Sum(ngscxymax(lhl,li,1:jye))
       
       k = k + 5
       
       ENDIF

       IF ( lf .gt. 1 ) THEN
       
       mpitotindp(k+1)  = Sum(ngscmaxy(lf,li,1:jye))
       mpitotindp(k+2)  = Sum(ngscminy(lf,li,1:jye))
       mpitotindp(k+3)  = Sum(ngscmaxy(lf,ls,1:jye))
       mpitotindp(k+4)  = Sum(ngscminy(lf,ls,1:jye))
       mpitotindp(k+5)  = Sum(ngscxymax(lf,li,1:jye))
       
       k = k + 5
       
       ENDIF
       
       ENDIF
#endif /* chgelec */

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, k, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)


       
      IF ( my_rank == 0 ) THEN
       psctot1 = mpitotoutdp(1)  
       psctot2 = mpitotoutdp(2)  
       sctot1p = mpitotoutdp(3)  
       sctot1n = mpitotoutdp(4)  
       iptotal = Int( mpitotoutdp(5) )
       iwetg  = Int( mpitotoutdp(6) )
       iwetg1 = Int( mpitotoutdp(7) )
       ndfalltot = Int( mpitotoutdp(8) )
       
      k = 8

      IF ( iptotal .gt. 0 ) THEN
        write(iunit,*) 'YIKES! Violated eqtot ',iptotal,' times!'
      END IF
      write(iunit,*) 'ptotalmn,ptotalmx = ',ptotalmn,ptotalmx
      write(iunit,*) 'wet growth at ', iwetg, ' points.  iwetg1 = ',iwetg1
      write(iunit,*) 'ndfalltot = ',ndfalltot
#ifdef CHGELEC
      IF ( ipelec .gt. 0 ) THEN
      write(iunit,*) 'Hit graupel-ice max ',  Int( mpitotoutdp(k+1) ),' times'
      write(iunit,*) 'Hit graupel-ice min ',  Int( mpitotoutdp(k+2) ),' times'
      write(iunit,*) 'Hit graupel-snow max ', Int( mpitotoutdp(k+3) ),' times'
      write(iunit,*) 'Hit graupel-snow min ', Int( mpitotoutdp(k+4) ),' times'
      write(iunit,*) 'Hit graupel-ice rate ', Int( mpitotoutdp(k+5) ),' times'
      k = k + 5
      IF ( lhl .gt. 1 ) THEN
      write(iunit,*) 'Hit hail-ice max ',   Int( mpitotoutdp(k+1) ),' times'
      write(iunit,*) 'Hit hail-ice min ',   Int( mpitotoutdp(k+2) ),' times'
      write(iunit,*) 'Hit hail-snow max ',  Int( mpitotoutdp(k+3) ),' times'
      write(iunit,*) 'Hit hail-snow min ',  Int( mpitotoutdp(k+4) ),' times'
      write(iunit,*) 'Hit hail-ice rate ',  Int( mpitotoutdp(k+5) ),' times'
      k = k + 5
      ENDIF

      IF ( lf .gt. 1 ) THEN
      write(iunit,*) 'Hit fdrop-ice max ',   Int( mpitotoutdp(k+1) ),' times'
      write(iunit,*) 'Hit fdrop-ice min ',   Int( mpitotoutdp(k+2) ),' times'
      write(iunit,*) 'Hit fdrop-snow max ',  Int( mpitotoutdp(k+3) ),' times'
      write(iunit,*) 'Hit fdrop-snow min ',  Int( mpitotoutdp(k+4) ),' times'
      write(iunit,*) 'Hit fdrop-ice rate ',  Int( mpitotoutdp(k+5) ),' times'
      ! k = k + 5
      ENDIF

      write(iunit,*) 'Max/Min psctot: ',psctotmx,psctotmn
      ENDIF ! ipelec
#endif /* chgelec */
      ENDIF ! my_rank

#else /* non-MPI */
      
      IF ( iptotal .gt. 0 ) THEN
        write(iunit,*) 'YIKES! Violated eqtot ',iptotal,' times!'
      END IF
      write(iunit,*) 'ptotalmn,ptotalmx = ',ptotalmn,ptotalmx
      write(iunit,*) 'wet growth at ', iwetg, ' points.  iwetg1 = ',iwetg1
      write(iunit,*) 'ndfalltot = ',ndfalltot
#ifdef CHGELEC
      IF ( ipelec .gt. 0 ) THEN
      write(iunit,*) 'Hit graupel-ice max ',  Sum(ngscmaxy(lh,li,1:jye)),' times'
      write(iunit,*) 'Hit graupel-ice min ',  Sum(ngscminy(lh,li,1:jye)),' times'
      write(iunit,*) 'Hit graupel-snow max ', Sum(ngscmaxy(lh,ls,1:jye)),' times'
      write(iunit,*) 'Hit graupel-snow min ', Sum(ngscminy(lh,ls,1:jye)),' times'
      write(iunit,*) 'Hit graupel-ice rate ', Sum(ngscxymax(lh,li,1:jye)),' times'
      IF ( lhl .gt. 1 ) THEN
      write(iunit,*) 'Hit hail-ice max ',  Sum(ngscmaxy(lhl,li,1:jye)),' times'
      write(iunit,*) 'Hit hail-ice min ',  Sum(ngscminy(lhl,li,1:jye)),' times'
      write(iunit,*) 'Hit hail-snow max ', Sum(ngscmaxy(lhl,ls,1:jye)),' times'
      write(iunit,*) 'Hit hail-snow min ', Sum(ngscminy(lhl,ls,1:jye)),' times'
      write(iunit,*) 'Hit hail-ice rate ', Sum(ngscxymax(lhl,li,1:jye)),' times'
      
      ENDIF

      IF ( lf .gt. 1 ) THEN
      write(iunit,*) 'Hit fdrop-ice max ',  Sum(ngscmaxy(lf,li,1:jye)),' times'
      write(iunit,*) 'Hit fdrop-ice min ',  Sum(ngscminy(lf,li,1:jye)),' times'
      write(iunit,*) 'Hit fdrop-snow max ', Sum(ngscmaxy(lf,ls,1:jye)),' times'
      write(iunit,*) 'Hit fdrop-snow min ', Sum(ngscminy(lf,ls,1:jye)),' times'
      write(iunit,*) 'Hit fdrop-ice rate ', Sum(ngscxymax(lf,li,1:jye)),' times'
      
      ENDIF


      write(iunit,*) 'Max/Min psctot: ',psctotmx,psctotmn
      ENDIF ! ipelec
#endif /* chgelec */
      
#endif /* MPI */
      
#ifdef CHGELEC
      IF ( ipelec .gt. 0 ) THEN
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,2(2x,1pe12.5))')    &
     &    'Charge total begin/end',   &
     &     psctot1, psctot2
!,(psctot1+psctot2)*dx*dy*dz
      write(iunit,'(a,3(2x,1pe12.5))')    &
     & 'Charge pos/neg/net',sctot1p,   &
     &   sctot1n,(sctot1p+sctot1n)
      ENDIF
      
      ENDIF
#endif /* chgelec */
!      call flush(iunit)
      


      IF ( io_flag .and. ( lairpress >= 1 .or. lairtem >= 1 .or. lsat >= 1 .or.  &
                lsati >= 1 .or. lrh3d >= 1 .or. lcwdia >= 1 ) ) THEN
#ifdef MPI
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

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
        ! temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz) ! t0(ix,jy,kz)
         temp1 = t0(ix,jy,kz)

          IF ( lairtem >= 1 )   axtra(ix,jy,kz,lairtem)  = temp1
          IF ( lairpress >= 1 ) axtra(ix,jy,kz,lairpress)  = pn(ix,jy,kz) + pb(kz)
          IF ( lthetav >= 1 ) axtra(ix,jy,kz,lthetav)  = an(ix,jy,kz,lt)*(1. + 0.61*an(ix,jy,kz,lv) )
          
          IF ( ( lsat >= 1 .or. lsati >= 1 .or. lrh3d >= 1 ) ) THEN

           temp1 = t0(ix,jy,kz)

           ltemq = Int( (temp1-163.15)/fqsat+1.5 )
           ltemq = Min( nqsat, Max(1,ltemq) )

            IF ( lsati >= 1 ) THEN
              c1 = t00(ix,jy,kz)*tabqis(ltemq)
              axtra(ix,jy,kz,lsati) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
            ENDIF

            IF ( lsat >= 1 .or. lrh3d >= 1 ) THEN
              IF ( iqvsopt == 0 ) THEN
                c1 = t00(ix,jy,kz)*tabqvs(ltemq)
              ELSE
                c1 = rdorv*esbolton*tabqvs(ltemq)/(pn(ix,jy,kz) + pb(kz) - esbolton*tabqvs(ltemq))
              ENDIF
              IF ( lsat >= 1 ) axtra(ix,jy,kz,lsat)  = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
              IF ( lrh3d >= 1 ) axtra(ix,jy,kz,lrh3d)  = 100.*(an(ix,jy,kz,lv)/c1)  ! from "new" values
            ENDIF

           ENDIF

            IF ( lcwdia >= 1 ) THEN
              IF ( an(ix,jy,kz,lc) > qxmin(lc) .and. an(ix,jy,kz,lnc) > cxmin ) THEN
                tmp = dn(ix,jy,kz)*an(ix,jy,kz,lc)/(1000.*an(ix,jy,kz,lnc))
                tmp = 1.e6*(3*tmp/(4.0*3.14159))**(1./3.)
                axtra(ix,jy,kz,lcwdia)  = tmp
              ELSE
               axtra(ix,jy,kz,lcwdia)  = 0.0
              ENDIF
            ENDIF
          

        ENDDO
       ENDDO
      ENDDO

      ENDIF

  ! Find max fall speed, max hail z-wgt speed
#ifdef MPI
! find global integrated rate max
       n = 0
       DO il = lc,lhab
         mpitotindp(n+1)  = vtmaxall(il)
         mpitotindp(n+2)  = vtmaxz(il)
         mpitotindp(n+3)  = vtmaxq(il)
         mpitotindp(n+4)  = vtmaxn(il)
         n = n + 4
       ENDDO

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
!       vtmaxall = mpitotoutdp(1)  
!       vtmaxzhl = mpitotoutdp(2)  
!      ENDIF
       n = 0
       DO il = lc,lhab
         vtmaxall(il) = mpitotoutdp(n+1)
         vtmaxz(il)   = mpitotoutdp(n+2)
         vtmaxq(il)   = mpitotoutdp(n+3)
         vtmaxn(il)   = mpitotoutdp(n+4)
         n = n + 4
       ENDDO

#endif

       DO il = lc,lhab
        write(91,'(a,i2,1x,4(f7.2))') 'Max fallspeed, Z-speed = ',il,vtmaxn(il),  vtmaxq(il), vtmaxz(il),vtmaxall(il)
       ENDDO

      deallocate ( db1 )
      deallocate ( db0 )
      deallocate ( dtz1 )
      deallocate ( dtz0 )
      deallocate ( z )
      deallocate ( xvt )
      deallocate ( tmpn )
      deallocate ( tmpn2 )
      deallocate ( alpha2d )
      
      deallocate ( ngscminy )
      deallocate ( ngscmaxy )
      deallocate ( ngscxymax )

#ifdef MPI
      if ( debug_mpi .and. my_rank>=0 ) write(0,*) my_rank, "END OF ICEZVD_DR - deallocated"
#endif

!        write(0,*) 'END TDMPHY: rank, time',my_rank,time_real

      return
      end
!
!  end of subroutine

! #######################################################################
!  HAILMAXD - calculated maximum expected hail size
! #######################################################################
#ifdef CCPPFLAG
!>\ingroup mod_nsslmp
!! Hail max size subroutine. 
#endif
     subroutine hailmaxd(dtp,nx,ny,nz,an,na,nor,norz,alpha2d,dn,qxmin,  &
     &                    hailmax1d,hailmaxk1,jslab,nalpha2d,its,ite,kts,kte ) 
!
! Calculate maximum hail size from the tail of of the distribution. The value
! of thresh_conc sets the minimum concentration in the integral over (Dmax, Inf).
! This uses the lookup tables for incomplete gamma functions and simply search for
! the expected value (and linearly interpolate) on D.
!
!  Written by ERM 7/2023
!
!
!
      USE INDEX_MODULE
      USE MICRO_MODULE
      
      implicit none

      integer nx,ny,nz,nor,norz,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
      integer, intent(in) :: its,ite,kts,kte ! x-z-range to calculate
      
      integer ng1
      parameter(ng1 = 1)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real qxmin(lc:lqmx)

!      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
      real dtp
      integer, intent(in) :: nalpha2d
      real alpha2d(-nor+1:nx+nor,1,-norz+1:nz+norz,nalpha2d)  ! array for PSD shape parameters
      real  :: hailmax1d(-nor+1:nx+nor,-nor+1:ny+nor),hailmaxk1(-nor+1:nx+nor,-nor+1:ny+nor)
      integer infdo
      integer jslab ! which line of xfall to use
            
      integer ix,jy,kz,ndfall,n,k,il,in
      double precision :: tmp, ratio, del, g1palp
      real, parameter :: dz = 200.

      real :: db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1)
      
      real :: rhovtzx(nz,nx)

      real :: alp, diam, diam1, hwdn
      
!      real, parameter :: cmin = 0.001 ! threshold number per m^3 for maximum diamter (threshold from diag_nwp)
!      DOUBLE PRECISION, PARAMETER:: thresh_conc = 0.0005d0                 ! number conc. of graupel/hail per cubic meter
      real, PARAMETER:: thresh_conc = 0.0005d0                 ! number conc. of graupel/hail per cubic meter
      real :: cwchtmp,cwchltmp, maxdia

!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: plo, phi
      integer :: ialp, i, j

      logical :: debug_mpi = .TRUE.

! ###################################################################


      IF ( lh > 1 ) THEN
        cwchtmp  = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
      ENDIF
      IF ( lhl > 1 ) THEN
        cwchltmp = ((3. + dnu(lhl))*(2. + dnu(lhl))*(1.0 + dnu(lhl)))**(-1./3.)
      ENDIF


      kzb = kts
      kze = kte

      ixb = its
      ixe = ite


      jy = jslab
      jgs = jy


!      hailmax1d(:,jy) = 0.0
!      hailmaxk1(:,jy) = 0.0

      if ( ndebug .gt. 0 ) write(0,*) 'dbg = 3a'


! first graupel, even if hail is also predicted, since graupel can sometime be large on its own
      IF ( lh > 1 .and. lnh > 1 ) THEN
      DO kz = kzb,kze
      DO ix = ixb,ixe
        IF ( (an(ix,jy,kz,lh) .gt. qxmin(lh)) .and. ( an(ix,jy,kz,lnh) .gt. thresh_conc) ) THEN
          IF ( lvh .gt. 1 ) THEN
            IF ( an(ix,jy,kz,lvh) > 1.e-20 ) THEN
            hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
            hwdn = Min(900.0,hwdn)
            ELSE
            hwdn = rho_qh
            ENDIF
          ELSE
            hwdn = rho_qh
          ENDIF

          tmp = 1. + alpha2d(ix,1,kz,2)
          i = Int(dgami*(tmp))
          del = tmp - dgam*i
          g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

          tmp = dn(ix,jy,kz)*an(ix,jy,kz,lh)/(hwdn*an(ix,jy,kz,lnh))
          diam = (6.0*tmp/(3.14159))**(1./3.)
          IF ( lzh > 1 ) THEN ! 3moment
            cwchtmp = ((3. + alpha2d(ix,1,kz,2))*(2. + alpha2d(ix,1,kz,2))*(1.0 + alpha2d(ix,1,kz,2)))**(-1./3.)
          ENDIF
          diam1 = diam*cwchtmp ! characteristic diameter, i.e., 1/lambda
         ! want cxd1 = thresh_conc
         !  tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
         ! cxd1 = cx(mgs,lh)*(tmp)/g1palp
         ! tmp = thresh_conc*g1palp/cx
         ! 
         tmp = thresh_conc*g1palp/an(ix,jy,kz,lnh)
         alp = alpha2d(ix,1,kz,2)
         ! gamxinflu(i,j,luindex,ilh)
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           ratio = 0.0
           maxdia = 0.0
           ! eventually could replace with bisection search, but final value of i is usually small
           ! compared to nqiacrratio
           DO i = 0,nqiacrratio-1
              IF ( gamxinflu(i,j,1,1) >= tmp .and. tmp >= gamxinflu(i+1,j,1,1) ) THEN
               !  interpolate here for FWIW
                ratio = i*dqiacrratio
                del = tmp - gamxinflu(i,j,1,1)
                ratio = (float(i) + del/(gamxinflu(i+1,j,1,1) - gamxinflu(i,j,1,1)))*dqiacrratio
                exit
              ENDIF
           ENDDO
           
           IF ( ratio > 0.0 ) THEN
              maxdia = ratio*diam1 ! units of m
           ENDIF

           IF ( kz == kzb ) THEN
             hailmaxk1(ix,jy) = Max( maxdia, hailmaxk1(ix,jy) )
!             IF ( maxdia > 0.1 ) THEN
!            IF ( an(ix,jy,kz,lh) > 1.e-4 ) THEN
!              write(0,*) 'maxdia,tmp,alp,ratio,diam,diam1= ',maxdia,tmp,alp,ratio,diam*100.,diam1*100.
!              write(0,*) 'hwdn, cxhl, qx, g1palp = ',hwdn, an(ix,jy,kz,lnhl), an(ix,jy,kz,lhl), g1palp
!              write(0,*) 'j,gamxinflu(0,2,4) = ',j,gamxinflu(0,j,1,1),gamxinflu(2,j,1,1), &
!                gamxinflu(4,j,1,1)
!            ENDIF
           ENDIF
           
           hailmax1d(ix,jy) = Max(maxdia, hailmax1d(ix,jy) )

        ! 

        ENDIF

      ENDDO
      ENDDO

      ENDIF ! lh

! then frozen drops
      IF ( lf > 1 .and. lnf > 1 ) THEN
      DO kz = kzb,kze
      DO ix = ixb,ixe
        IF ( an(ix,jy,kz,lf) .gt. qxmin(lf) .and. (an(ix,jy,kz,lnf) .gt. thresh_conc) ) THEN
          IF ( lvf .gt. 1 ) THEN
            IF ( an(ix,jy,kz,lvf) > 1.e-20 ) THEN
            hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
            hwdn = Min(900.0,hwdn)
            ELSE
            hwdn = rho_qf
            ENDIF
          ELSE
            hwdn = rho_qf
          ENDIF

          tmp = 1. + alpha2d(ix,1,kz,4)
          i = Int(dgami*(tmp))
          del = tmp - dgam*i
          g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

          tmp = dn(ix,jy,kz)*an(ix,jy,kz,lf)/(hwdn*an(ix,jy,kz,lnf))
          diam = (6.0*tmp/(3.14159))**(1./3.)
          IF ( lzf > 1 ) THEN ! 3moment
            cwchtmp = ((3. + alpha2d(ix,1,kz,4))*(2. + alpha2d(ix,1,kz,4))*(1.0 + alpha2d(ix,1,kz,4)))**(-1./3.)
          ENDIF
          diam1 = diam*cwchtmp ! characteristic diameter, i.e., 1/lambda
         ! want cxd1 = thresh_conc
         !  tmp = gaminterp(ratio,alpha(mgs,lf),1,1)
         ! cxd1 = cx(mgs,lf)*(tmp)/g1palp
         ! tmp = thresh_conc*g1palp/cx
         ! 
         tmp = thresh_conc*g1palp/an(ix,jy,kz,lnf)
         alp = alpha2d(ix,1,kz,4)
         ! gamxinflu(i,j,luindex,ilh)
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           ratio = 0.0
           maxdia = 0.0
           ! eventually could replace with bisection search, but final value of i is usually small
           ! compared to nqiacrratio
           DO i = 0,nqiacrratio-1
              IF ( gamxinflu(i,j,1,1) >= tmp .and. tmp >= gamxinflu(i+1,j,1,1) ) THEN
               !  interpolate here for FWIW
                ratio = i*dqiacrratio
                del = tmp - gamxinflu(i,j,1,1)
                ratio = (float(i) + del/(gamxinflu(i+1,j,1,1) - gamxinflu(i,j,1,1)))*dqiacrratio
                exit
              ENDIF
           ENDDO
           
           IF ( ratio > 0.0 ) THEN
              maxdia = ratio*diam1 ! units of m
           ENDIF

           IF ( kz == kzb ) THEN
             hailmaxk1(ix,jy) = Max( maxdia, hailmaxk1(ix,jy) )
!             IF ( maxdia > 0.1 ) THEN
!            IF ( an(ix,jy,kz,lh) > 1.e-4 ) THEN
!              write(0,*) 'maxdia,tmp,alp,ratio,diam,diam1= ',maxdia,tmp,alp,ratio,diam*100.,diam1*100.
!              write(0,*) 'hwdn, cxhl, qx, g1palp = ',hwdn, an(ix,jy,kz,lnhl), an(ix,jy,kz,lhl), g1palp
!              write(0,*) 'j,gamxinflu(0,2,4) = ',j,gamxinflu(0,j,1,1),gamxinflu(2,j,1,1), &
!                gamxinflu(4,j,1,1)
!            ENDIF
           ENDIF
           
           hailmax1d(ix,jy) = Max(maxdia, hailmax1d(ix,jy) )

        ! 

        ENDIF

      ENDDO
      ENDDO

      ENDIF ! lh

! And diam for hail if present
      IF ( lhl > 1 .and. lnhl > 1 ) THEN
      DO kz = kzb,kze
      DO ix = ixb,ixe
        IF ( an(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. an(ix,jy,kz,lnhl) .gt. thresh_conc ) THEN
          IF ( lvhl .gt. 1 ) THEN
            IF ( an(ix,jy,kz,lvhl) > 1.e-20 ) THEN
            hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
            hwdn = Min(900.0,hwdn)
            ELSE
            hwdn = rho_qhl
            ENDIF
          ELSE
            hwdn = rho_qhl
          ENDIF

          tmp = 1. + alpha2d(ix,1,kz,3)
          i = Int(dgami*(tmp))
          del = tmp - dgam*i
          g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

          tmp = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/(hwdn*an(ix,jy,kz,lnhl))
          diam = (6.0*tmp/(3.14159))**(1./3.)
          IF ( lzhl > 1 ) THEN ! 3moment
            cwchltmp = ((3. + alpha2d(ix,1,kz,3))*(2. + alpha2d(ix,1,kz,3))*(1.0 + alpha2d(ix,1,kz,3)))**(-1./3.)
          ENDIF
          diam1 = diam*cwchltmp ! characteristic diameter, i.e., 1/lambda
         ! want cxd1 = thresh_conc
         !  tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
         ! cxd1 = cx(mgs,lh)*(tmp)/g1palp
         ! tmp = thresh_conc*g1palp/cx
         ! 
         tmp = thresh_conc*g1palp/an(ix,jy,kz,lnhl)
         alp = alpha2d(ix,1,kz,3)
         ! gamxinflu(i,j,luindex,ilh)
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           ratio = 0.0
           maxdia = 0.0
           ! eventually could replace with bisection search, but final value of i is usually small
           ! compared to nqiacrratio
           DO i = 0,nqiacrratio-1
              IF ( gamxinflu(i,j,1,1) >= tmp .and. tmp >= gamxinflu(i+1,j,1,1) ) THEN
               !  interpolate here for FWIW
                ratio = i*dqiacrratio
                del = tmp - gamxinflu(i,j,1,1)
                ratio = (float(i) + del/(gamxinflu(i+1,j,1,1) - gamxinflu(i,j,1,1)))*dqiacrratio
                exit
              ENDIF
           ENDDO
           
           IF ( ratio > 0.0 ) THEN
              maxdia = ratio*diam1 ! units of m
           ENDIF

           IF ( kz == kzb ) THEN
             hailmaxk1(ix,jy) = Max( maxdia, hailmaxk1(ix,jy) )
!             IF ( maxdia > 0.1 ) THEN
!            IF ( an(ix,jy,kz,lhl) > 1.e-4 ) THEN
!              write(0,*) 'maxdia,tmp,alp,ratio,diam,diam1= ',maxdia,tmp,alp,ratio,diam*100.,diam1*100.
!              write(0,*) 'hwdn, cxhl, qx, g1palp = ',hwdn, an(ix,jy,kz,lnhl), an(ix,jy,kz,lhl), g1palp
!              write(0,*) 'j,gamxinflu(0,2,4) = ',j,gamxinflu(0,j,1,1),gamxinflu(2,j,1,1), &
!                gamxinflu(4,j,1,1)
!            ENDIF
           ENDIF
           
           hailmax1d(ix,jy) = Max(maxdia, hailmax1d(ix,jy) )

        ! 

        ENDIF

      ENDDO
      ENDDO

      ENDIF


     END SUBROUTINE HAILMAXD
! #######################################################################


#ifdef CM1
      subroutine crowfall(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt, & 
     &  an,db1,imapz,mzdist,il,idum,xfall)
      
      write(0,*) 'crowfall not supported in CM1!  Use boxfall (itfall=1)'
      STOP
      
      RETURN
      END


#endif

! #######################################################################
#ifdef CCPPFLAG
!>\ingroup mod_nsslmp
!! Sedimentation driver subroutine. Calls fallout column by column
#endif
     subroutine sediment1d(dtp,nx,ny,nz,an,na,nor,norz,xfall,dn,z1d, &
     &                    t0,t7,axtra,io_flag,infdo,jslab,ixe,kze,    &
     &                    ldovol,rho00,cdx,cno,ido,ln,lz,lvol,lsc,ipc,lliq,lrain, &
     &                    xdn0,xvmn,xvmx,xdnmn,xdnmx,qxmin,cwnccnold,  &
     &                    vtmaxz, vtmaxq, vtmaxn, vtmaxall, dbzchange,linfall)
!
! Sedimentation driver -- column by column
!
!  Written by ERM 10/2011
!
!
!
!      USE INDEX_MODULE
      USE COMMASMPI_MODULE, only: my_rank
      USE INDEX_MODULE, only: lt,lc,lr,li,lis,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lhw,lfw,lhlw, &
                               lnc,lnr,lni,lnis,lns,lnh,lnf,lnhl,cinu,dmuh,dnu,dmu,xnu,xmu,dmuhl,rnu,cnu,snu, &
                               lss,lsat,lsati,xvcmx,xvcmn,xvrmn,xvrmx,lccn,rnumin,rnumax, &
                               lqmx,nxtra,lvi,lvs,lvh,lvf,lvhl,lzi,lzs,lzr,lzh,lzf,lzhl,lsw, &
                               alphar,alphamin,alphamax, lnhf,lnhlf, nraintypes, iraintypes, &
                               cwmasn, cwmasx, cwradn, lmvhl, lzvhl, lnvhl, lmnvhl, ldbzchangeh, ldbzchanger
      USE MICRO_MODULE, only: ipconc, ipelec, infall, irfall, isfall, isedonly

      implicit none

      integer nx,ny,nz,nor,norz,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
!      integer :: its,jts ! SW point of local tile
      
      integer ng1
      parameter(ng1 = 1)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real dz3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real z1d(-norz+ng1:nz+norz,4)
!      real dz3dinv(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)
     
      integer, intent(in) :: lsc(lc:lhab)
      integer, intent(in) :: ln(lc:lhab)
      integer, intent(in) :: ipc(lc:lhab)
      integer, intent(in) :: lvol(lc:lhab)
      integer, intent(in) :: lz(lc:lhab)
      integer, intent(in) :: lliq(li:lhab)
      integer, intent(in) :: lrain(nraintypes)
      integer, intent(in) :: ido(lc:lqmx)
      integer, intent(in) :: linfall(lc:lqmx)
      real,    intent(in) :: xdnmx(lc:lqmx), xdnmn(lc:lqmx), qxmin(lc:lqmx)
      real,    intent(in) :: xdn0(lc:lqmx), xvmn(lc:lhab), xvmx(lc:lhab)
      real,    intent(in) :: cwnccnold(nz)
      real,    intent(inout) :: vtmaxz(lc:lhab), vtmaxq(lc:lhab),vtmaxn(lc:lhab),vtmaxall(lc:lhab)
      real,    intent(inout) :: dbzchange(2,nz)

      
      logical, intent(in) :: ldovol, io_flag

      real, intent(in) ::  rho00
      real, intent(in) ::  cdx(lc:lhab)
      real, intent(in) ::  cno(lc:lhab)
      integer, intent(in) :: ixe, kze

!      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
      real dtp
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      integer infdo
      integer jslab ! which line of xfall to use
            
      integer ix,jy,kz,ndfall,n,k,il,in
      real tmp, vtmax, dtptmp, dtfrac
      double precision :: tmpdp
      real, parameter :: dz = 200.

      real, allocatable :: db1(:,:), dtz1(:,:,:),dz2dinv(:,:),db1inv(:,:) ! db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1)
      real, allocatable :: rhovtzx(:,:)
      real, allocatable :: xfall0(:,:), xvt(:,:,:,:),tmpn(:,:,:),tmpn2(:,:,:),z(:,:,:),tmpn3(:,:,:)
      

      integer :: ngs ! = 512
      integer :: ngscnt,mgs,ipconc0
      

      real, allocatable ::  qx(:,:)
      real, allocatable ::  qxw(:,:)
      real, allocatable ::  cx(:,:)
      real, allocatable ::  xv(:,:)
      real, allocatable ::  vtxbar(:,:,:)
      real, allocatable ::  xmas(:,:)
      real, allocatable ::  xdn(:,:)
      real, allocatable ::  xdia(:,:,:)
      real, allocatable ::  vx(:,:)
      real, allocatable ::  alpha(:,:)
      real, allocatable ::  zx(:,:)
      logical, allocatable :: hasmass(:,:)

      integer, allocatable :: igs(:),kgs(:)
      
      real, allocatable :: rho0(:),temcg(:)

      real, allocatable :: temg(:)
      
      real, allocatable :: rhovt(:)
      
      real, allocatable :: cwnc(:),cinc(:)
      real, allocatable :: fadvisc(:),cwdia(:),cipmas(:)
      
      real, allocatable :: cnina(:),cimas(:)
      
      real, allocatable :: cnostmp(:)
      
      real :: cimasn,cimasx
      
      integer, parameter :: interval_sedi_vt = 1 ! interval for recalculating Vt in sedimentation subloop (only when do_accurate_sedimentation = .true.)
      real   , parameter :: vtmaxsed = 70. ! Limit on fall speed (m/s, all moments) for sedimentation calculations. Not applied to fall speeds for microphysical rates
      logical, parameter :: do_accurate_sedimentation = .true.
      integer, parameter :: ndebug = 0
!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb, m, k1
      integer :: jye
      integer :: plo, phi

! ###################################################################


      allocate( db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1),rhovtzx(nz,nx) )
      allocate( xfall0(nx,ny), xvt(nz+1,nx,3,lc:lhab), tmpn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate( tmpn2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz), z(-nor+ng1:nx+nor,-norz+ng1:nz+norz,lr:lhab))
      allocate( tmpn3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz))

        tmpn(:,:,:) = 0.0
       tmpn2(:,:,:) = 0.0
       tmpn3(:,:,:) = 0.0
           z(:,:,:) = 0.0

      ngs = nz+3
      
      allocate( qx(ngs,lv:lhab),  &
                qxw(ngs,ls:lhab),  &
                cx(ngs,lc:lhab),  &
                xv(ngs,lc:lhab),  &
                vtxbar(ngs,lc:lhab,3),  &
                xmas(ngs,lc:lhab),  &
                xdn(ngs,lc:lhab),  &
                xdia(ngs,lc:lhab,3),  &
                vx(ngs,li:lhab),  &
                alpha(ngs,lc:lhab),  &
                zx(ngs,lr:lhab),     &
                hasmass(nx,lc+1:lhab), &
                igs(ngs),kgs(ngs), &
                rho0(ngs),temcg(ngs),temg(ngs), rhovt(ngs), &
                cwnc(ngs),cinc(ngs), &
                fadvisc(ngs),cwdia(ngs),cipmas(ngs), &
                cnina(ngs),cimas(ngs), &
                cnostmp(ngs) )


      kzb = 1
!      kze = nz

      ixb = 1
!      ixe = nx


      jy = jslab
      jgs = jy

!      write(0,*) 'sed1d: ixe,kze = ',ixe,kze

!
!  zero the precip flux arrays (2d)
!

      xvt(:,:,:,:) = 0.0

      if ( ndebug .gt. 0 ) write(0,*) 'dbg = 3a'


      DO kz = kzb,nz
       k1 = Min(kz,nz-1)
      DO ix = ixb,ixe
       IF ( dn(ix,jy,k1) == 0.0 ) THEN
         write(0,*) 'dn zero at i,j,k = ',ix,jy,k1
       ENDIF
       db1(ix,kz) = dn(ix,jy,k1)
       db1inv(ix,kz) = 1./dn(ix,jy,k1)
       rhovtzx(kz,ix) = Sqrt(rho00*Min(1.0/0.05, db1inv(ix,kz))) ! prevent excessive rhovt
      ENDDO
      ENDDO

      DO kz = kzb,nz
       k1 = Min(kz,nz-1)
      DO ix = ixb,ixe
       dtz1(kz,ix,0) = z1d(k1,3) ! dz3dinv(ix,jy,kz)
       dtz1(kz,ix,1) = z1d(k1,3)*db1inv(ix,kz) 
!       dtz1(ix,kz) = z1d(kz,3)/db1(ix,kz)
       dz2dinv(kz,ix) = z1d(k1,3) ! dz3dinv(ix,jy,kz)
      ENDDO
      ENDDO

      IF ( lzh .gt. 1 ) THEN
!       DO kz = kzb,kze
!       DO ix = ixb,ixe
!         an(ix,jy,kz,lzh) = Max( 0., an(ix,jy,kz,lzh) )
!       ENDDO
!       ENDDO
      ENDIF

      
      DO il = lc+1,lhab
       DO ix = ixb,ixe
!        hasmass(ix,il) = Any( an(ix,jy,:,il) > qxmin(il) )
       ENDDO
      ENDDO




      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3a2'

! loop over columns
      DO ix = ixb,ixe
      
#if defined(MPI) & defined(TIMINGS)
         dt1 = MPI_Wtime()
#endif

      call ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ix, & 
     &  xvt, rhovtzx, & 
     &  an,dn,ipconc,t0,t7,cwmasn,cwmasx, & 
     &  cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx, & 
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,ln,lz,lvol,lliq,              &
     &  infdo,0,cwnccnold     &
     & )

       IF ( io_flag  ) THEN
        vtmax = 0.0
        DO kz = kzb,nz-1
        
          IF ( lhl > 1 .and. (lzvhl > 1) ) THEN
          
          IF ( lmvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lmvhl) = xvt(kz,ix,1,lhl)
             vtmax = Max(vtmax,xvt(kz,ix,1,lhl))
          ENDIF
          IF ( lzvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lzvhl) = xvt(kz,ix,3,lhl)
          ENDIF
          IF ( lnvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lnvhl) = xvt(kz,ix,2,lhl)
             IF ( xvt(kz,ix,2,lhl) > 0. .and. lmnvhl > 0 ) THEN
               axtra(ix,jgs,kz,lmnvhl) = xvt(kz,ix,1,lhl) - xvt(kz,ix,2,lhl)
             ENDIF
          ENDIF
          
          ELSEIF ( lh > 1 .and. (lzvhl > 1) ) THEN
          ! if hail is turned off, then use graupel values

          IF ( lmvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lmvhl) = xvt(kz,ix,1,lh)
             vtmax = Max(vtmax,xvt(kz,ix,1,lh))
          ENDIF
          IF ( lzvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lzvhl) = xvt(kz,ix,3,lh)
          ENDIF
          IF ( lnvhl >= 1 ) THEN
             axtra(ix,jgs,kz,lnvhl) = xvt(kz,ix,2,lh)
             IF ( xvt(kz,ix,2,lh) > 0. .and. lmnvhl > 0 ) THEN
               axtra(ix,jgs,kz,lmnvhl) = xvt(kz,ix,1,lh) - xvt(kz,ix,2,lh)
             ENDIF
          ENDIF
          
          ENDIF
          
         
         ENDDO
!         write(0,*) 'vtmax = ',vtmax
       ENDIF

! loop over each species and do sedimentation for all moments
     DO il = lc,lhab
       IF ( ido(il) == 0 ) CYCLE

!       IF ( .not. hasmass(ix,il) ) CYCLE

!      plo = nz
!      phi = 0


      vtmax = 0.0
      
      do kz = kzb,kze

       vtmaxz(il) = Max( vtmaxz(il), xvt(kz,ix,3,il) ) ! track max reflectivity-wgt fall speed
       vtmaxq(il) = Max( vtmaxq(il), xvt(kz,ix,1,il) ) ! track max mass-wgt fall speed
       vtmaxn(il) = Max( vtmaxn(il), xvt(kz,ix,2,il) ) ! track max number-wgt fall speed


      ! apply limit vtmaxsed (08/20/2015)
      xvt(kz,ix,1,il) = Min( vtmaxsed,  xvt(kz,ix,1,il) )
      xvt(kz,ix,2,il) = Min( vtmaxsed,  xvt(kz,ix,2,il) )
      xvt(kz,ix,3,il) = Min( vtmaxsed,  xvt(kz,ix,3,il) )
      
      vtmax = Max(vtmax,xvt(kz,ix,1,il)*dz2dinv(kz,ix))
      vtmax = Max(vtmax,xvt(kz,ix,2,il)*dz2dinv(kz,ix))
      vtmax = Max(vtmax,xvt(kz,ix,3,il)*dz2dinv(kz,ix))

       vtmaxall(il) = Max( vtmaxall(il), vtmax ) ! 
      
      ENDDO
      
      IF ( vtmax == 0.0 ) CYCLE


      
      IF ( dtp*vtmax .lt. 0.7 ) THEN ! check whether multiple steps are needed.
        ndfall = 1
      ELSE
       IF ( dtp > 20.0 ) THEN ! more stringent subdivision for large time steps
         ndfall = Max(2, Int(dtp*vtmax/0.7) + 1)
       ELSE ! more relaxed for small time steps, but might still be a problem for very thin vertical layers near the ground
         ndfall = 1+Int(dtp*vtmax + 0.301)
       ENDIF
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp/Real(ndfall)
!        write(0,*) 'subdivide fallout on its,jts,ix,plo,phi = ',its,jts,ix,plo,phi
!        write(0,*) 'for il,jsblab,c,ndfall = ',il,jslab,dtp*vtmax,ndfall
      ELSE
        dtptmp = dtp
      ENDIF
      
      dtfrac = dtptmp/dtp


      DO n = 1,ndfall

      IF ( do_accurate_sedimentation .and. n .ge. 2 .and. &
           ( n == interval_sedi_vt*(n/interval_sedi_vt) ) ) THEN !{
!
!  For n >= 2, zero the precip flux arrays (2d) and get new fall speeds
!
      

      xvt(1:nz+1,ix,1:3,il) = 0.0 ! reset to zero because routine will only compute points  
                                  ! with q > qmin (i.e., does not always set every point)

      call ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ix, & 
     &  xvt, rhovtzx, & 
     &  an,dn,ipconc,t0,t7,cwmasn,cwmasx, & 
     &  cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx, & 
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,ln,lz,lvol,lliq,              &
     &  infdo,il,cwnccnold)

      DO kz = kzb,kze
      ! apply limit vtmaxsed (08/20/2015)
        xvt(kz,ix,1,il) = Min( vtmaxsed,  xvt(kz,ix,1,il) )
        xvt(kz,ix,2,il) = Min( vtmaxsed,  xvt(kz,ix,2,il) )
        xvt(kz,ix,3,il) = Min( vtmaxsed,  xvt(kz,ix,3,il) )
      ENDDO




      ENDIF ! } (n .ge. 2)

!       IF ( ix == 10 .and. il == 8 ) THEN
!         DO kz = kze,kze-10,-1
!           write(91,*) 'k,q,c,vt123 = ',kz,an(ix,jy,kz,il),an(ix,jy,kz,ln(il)),an(ix,jy,kz,lz(il))
!           write(91,*) '              ',xvt(kz,ix,1,il),xvt(kz,ix,2,il),xvt(kz,ix,3,il)
!         ENDDO
!       ENDIF


        IF ( il >= lr .and. ( linfall(il) .eq. 3 .or. linfall(il) .eq. 4 ) .and. ln(il) > 0 .and. lz(il) <= 0 ) THEN
            call calczgr(nx,ny,nz,nor,na,an,ixe,kze, & 
     &         z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il), &
               lvol(il), xdn0(il), ix )
        ELSEIF ( ( il == lh .and. ldbzchangeh >= 1) .or. ( il == lr .and. ldbzchanger >= 1 ) ) THEN
            IF ( lz(il) <= 0 ) THEN
            call calczgr(nx,ny,nz,nor,na,an,ixe,kze, & 
     &         z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il), &
               lvol(il), xdn0(il), ix )
            ENDIF
        ENDIF

      if (ndebug .gt. 0 ) write(0,*) 'dbg = 1b'

! mixing ratio

      call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &             an,db1,il,1,xfall,dtz1,ix)

       IF ( il == lr ) THEN
       IF ( iraintypes >= 1 ) THEN
       DO m = 1,nraintypes
        IF ( lrain(m) > 1 ) THEN
!       call fallout(nx,ny,nz,nor,na,gz,z1d,dtptmp,dz,jgs,xvt(1,1,1,lr), & 
!     &             an,db1,imapz,mzdist,lrain(m),1,xfall,dtz1)
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &             an,db1,lrain(m),1,xfall,dtz1,ix)
        ENDIF
       ENDDO
       ENDIF
       ENDIF

       IF ( il >= ls ) THEN
       IF ( lliq(il) .gt. 1 ) THEN
       call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &             an,db1,lliq(il),1,xfall,dtz1,ix)
       ENDIF
       ENDIF

      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3c'

! volume

      IF ( ldovol .and. il >= li ) THEN
        IF ( lvol(il) .gt. 1 ) THEN
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &              an,db1,lvol(il),0,xfall,dtz1,ix)
        ENDIF
      ENDIF

#if defined( CHGELEC ) && !defined( CHGNIONLY )
! charge

      IF ( ipelec > 0 ) THEN
        IF ( lsc(il) .gt. 1 ) THEN
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), &
     &              an,db1,lsc(il),0,xfall,dtz1,ix)
        ENDIF
      ENDIF
#endif
! reflectivity

      IF ( ipconc .ge. 6 ) THEN
        IF ( lz(il) .gt. 1 ) THEN
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,3,il), & 
     &              an,db1,lz(il),0,xfall,dtz1,ix)
        
        ENDIF
      ENDIF


      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3d'

      
      IF ( ipconc .gt. 0 ) THEN !{
      IF ( ipconc .ge. ipc(il) ) THEN

      IF ( ( linfall(il) .ge. 2 .or. (infall .eq. 0 .and. il .lt. lh) ) .and. lz(il) .lt. 1) THEN !{
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!

        IF ( linfall(il) == 3 .or. linfall(il) == 4 ) THEN
          ! set up for method I or I+II
          DO kz = kzb,kze
              tmpn2(ix,jy,kz) = z(ix,kz,il)
          ENDDO
          DO kz = kzb,kze
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
          ENDDO

        ELSE
          ! set up for method II only
          DO kz = kzb,kze
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
          ENDDO


        ENDIF

      ENDIF !}

          IF ( ( il == lh .and. ldbzchangeh >= 1) .or. ( il == lr .and. ldbzchanger >= 1 ) ) THEN
           ! save pre-sedimentation z array into tmpn3
            DO kz = kzb,kze
              tmpn3(ix,jy,kz) = z(ix,kz,il)
            ENDDO
          ENDIF

      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3f'

       in = 2
       IF ( linfall(il) .eq. 1 ) in = 1

         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,in,il), & 
     &        an,db1,ln(il),0,xfall,dtz1,ix)


         IF ( lz(il) .lt. 1 ) THEN ! { if not 3-moment, run one of the correction schemes
           IF ( linfall(il) >= 2 ) THEN  !{
           xfall0(:,jgs) = 0.0

           IF ( linfall(il) == 3 .or. linfall(il) == 4 ) THEN
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,3,il), & 
     &         tmpn2,db1,1,0,xfall0,dtz1,ix) ! sediment temporary Z moment with Z-wgt Vt
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &         tmpn,db1,1,0,xfall0,dtz1,ix) ! sediment temporary N moment with Z-wgt Vt
           ELSE
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &         tmpn,db1,1,0,xfall0,dtz1,ix)
           ENDIF

           IF ( linfall(il) == 3 .or. linfall(il) == 4 ) THEN !{
           ! "Method I" - dbz correction
           ! Uses input tmpn2 (temp. Z-moment) to determine if new N and q values in an(:,:,:,ln(il))
           ! cause an increase in reflectivity moment. If so, either use N from mass-wgt Vt (tmpn) to replace 
           ! new N (infall=3; I) or use smaller N from tmpn or calculated from q and temporary Z (infall=4; I+II)
           ! Uses 'z' array to check if new reflectivity is greater than pre-sedimentation reflectivity
             call calcnfromz(nx,ny,nz,nor,na,an,tmpn2,ixe,kze, & 
     &       z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il),tmpn,  & 
     &       lvol(il), xdn0(il), linfall(il), ix)

           ELSEIF ( linfall(il) .eq. 5 .and. il .ge. lh .or. ( il == lr .and. irfall == 5 ) ) THEN

             DO kz = kzb,kze
               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), 0.5* ( an(ix,jgs,kz,ln(il)) + tmpn(ix,jy,kz) ))
             ENDDO           

           ELSEIF ( .not. (il .eq. lr .and. irfall .eq. 0) .and. .not. (il .eq. ls .and. isfall .eq. 0) ) THEN
! "Method II" M-wgt N-fallout correction

             DO kz = kzb,kze
               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
             ENDDO
           ENDIF  !}

           ENDIF !}

             IF ( ( il == lh .and. ldbzchangeh >= 1) .or. ( il == lr .and. ldbzchanger >= 1 ) ) THEN !{
             ! finding reflectivity increases from sedimentation
!               IF ( .not. (infall == 3 .or. infall == 4 ) ) THEN
                call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,3,il), & 
     &                         tmpn3,db1,1,0,xfall0,dtz1,ix)
!               ENDIF

               ! calculate new species reflectivity
               call calczgr(nx,ny,nz,nor,na,an,ixe,kze, & 
     &           z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il), &
                 lvol(il), xdn0(il), ix )

                IF ( il == lr ) THEN
                 DO kz = kzb,kze
                   IF (  z(ix,kz,il) > 1.e-28 .and. tmpn3(ix,jy,kz) > 1.e-28 .and. z(ix,kz,il) > tmpn3(ix,jy,kz) ) THEN
                   IF ( 10.*log10(1.e18*z(ix,kz,il)) > 10.0 .and. 10.*log10(1.e18*tmpn3(ix,jy,kz)) > 10.0 ) THEN
                    tmpdp = 10.*(Log10(z(ix,kz,il)) - Log10(tmpn3(ix,jy,kz)))
                    axtra(ix,jy,kz,ldbzchanger) = axtra(ix,jy,kz,ldbzchanger) + tmpdp
                    dbzchange(1,kz) = dbzchange(1,kz) + tmpdp
                   ENDIF
                   ENDIF
                 ENDDO
                ENDIF

                IF ( il == lh ) THEN
                 DO kz = kzb,kze
                   IF (  z(ix,kz,il) > 1.e-28 .and. tmpn3(ix,jy,kz) > 1.e-28 .and. z(ix,kz,il) > tmpn3(ix,jy,kz) ) THEN
                   IF ( 10.*log10(1.e18*z(ix,kz,il)) > 10.0 .and. 10.*log10(1.e18*tmpn3(ix,jy,kz)) > 10.0 ) THEN
                    tmpdp = 10.*(Log10(z(ix,kz,il)) - Log10(tmpn3(ix,jy,kz)))
                    axtra(ix,jy,kz,ldbzchangeh) = axtra(ix,jy,kz,ldbzchangeh) + tmpdp !  + 10.*(Log10(z(ix,kz,il)) - Log10(tmpn3(ix,jy,kz)))
                    dbzchange(2,kz) = dbzchange(2,kz) + tmpdp
                   ENDIF
                   ENDIF
                 ENDDO
                ENDIF

             ENDIF !} il==lh

           

         ENDIF !} lz(il) .lt. 1
        ENDIF ! ipconc > ipc


      ENDIF !}

      ! If have 3M hail and frozen drops, then do number of hail from frozen drops (if lnhlf > 1)
      ! Need to have 3M so that no number correction is done (no artificial breakup)
       IF ( lhl > 1 .and. il == lhl .and. lnhlf > 1 .and. lzhl > 1 ) THEN
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,in,il), &
     &              an,db1,lnhlf,0,xfall,dtz1,ix)
       ENDIF

       IF ( il == lh .and. lnhf > 1 .and. lzh > 1 ) THEN
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,in,il), &
     &              an,db1,lnhf,0,xfall,dtz1,ix)
       ENDIF


      ENDDO ! n=1,ndfall
      ENDDO ! il
      
      ENDDO ! ix


      deallocate( db1,dtz1,dz2dinv,db1inv,rhovtzx )
      deallocate( xfall0, xvt, tmpn )
      deallocate( tmpn2, tmpn3, z)

      deallocate( qx,  &
                qxw,  &
                cx,  &
                xv,  &
                vtxbar,  &
                xmas,  &
                xdn,  &
                xdia,  &
                vx,  &
                alpha,  &
                zx,     &
                hasmass, &
                igs,kgs, &
                rho0,temcg,temg, rhovt, &
                cwnc,cinc, &
                fadvisc,cwdia,cipmas, &
                cnina,cimas, &
                cnostmp )

      RETURN
      END SUBROUTINE SEDIMENT1D


! #####################################################################

!
! #####################################################################


!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
#ifdef CCPPFLAG
!>\ingroup mod_nsslmp
!! Column sedimentation fallout subroutine
#endif
      subroutine fallout1d(nx,ny,nz,nor,na,dtp,dtfrac,jgs,vt,   &
     &  a,db1,ia,id,xfall,dtz1,ixcol)
!
! First-order, upwind fallout scheme
!
!  Written by ERM 6/10/2011
!
!
!
      implicit none

      integer nx,ny,nz,nor,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)
      integer :: ixcol

!      real dz3dinv(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
!      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
      real vt(nz+1,nx)  ! terminal speed for a
      real dtp,dtfrac
      real cmax
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      real db1(nx,nz+1),dtz1(nz+1,nx,0:1)

! Local
           
      integer ix,jy,kz,n,k,km1
      integer iv1,iv2
      real tmp
      integer imn,imx,kmn,kmx
      real qtmp1(nz+1)

!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

! ###################################################################

      jy = jgs

      iv1 = 0
      iv2 = 0

      imn = nx
      imx = 1
      kmn = nz
      kmx = 1

      cmax = 0.0

      kzb = 1
      kze = nz-1 ! nz-1 is top scalar var level in commas

      ix  = ixcol

      qtmp1(kze+1) = 0.0
     
      DO kz = kzb,kze
         
         IF ( id == 1 ) THEN
           qtmp1(kz) = a(ix,jgs,kz,ia)*vt(kz,ix)*db1(ix,kz)
         ELSE
           qtmp1(kz) = a(ix,jgs,kz,ia)*vt(kz,ix)
         ENDIF
         
         IF ( a(ix,jgs,kz,ia) .ne. 0.0 ) THEN
           kmn = Min(kz,kmn)
           kmx = Max(kz,kmx)
         ENDIF
      ENDDO
            
      kmn = Max(1,kmn-1)
      
! first check if fallout is worth doing
!      IF ( cmax .eq. 0.0 .or. imn .gt. imx ) THEN
!        RETURN
!      ENDIF
      
      IF ( kmn == 1 ) THEN
      
      kz = 1
         xfall(ix,jy,ia) = xfall(ix,jy,ia) + dtp*dtz1(kz,ix,id)*a(ix,jgs,kz,ia)*vt(kz,ix)
      
      ENDIF

      do kz = 1,kze
        a(ix,jgs,kz,ia) =  a(ix,jgs,kz,ia) + dtp*dtz1(kz,ix,id)*(qtmp1(kz+1) - qtmp1(kz) )
      enddo

      
      RETURN
      END SUBROUTINE FALLOUT1D

! ##############################################################################
! ##############################################################################

