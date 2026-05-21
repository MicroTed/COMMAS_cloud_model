


! prepocessed on "Dec  2 2024" at "13:32:34"





!c--------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------
!

      subroutine nssl_2mom_gs   &
     &  (nx,ny,nz,na,jyslab  &
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
     &   xdn0,tmp3d,tkediss  &
     &  ,dxx,dyy,dzz,z1d,io_flag  &
     &  ,ventrn,qxmin, xvmn, xvmx, cno  &
     &  ,cckm,ccne,ccnefac,cnexp,ccne0  &
     &  ,lsc,ln,ipc,lvol,lz,lliq,lrain &
     &  ,dab0,dab1,da0,da1,bb &
     &  ,dab0lu,dab1lu &
     &  ,iexy,takalu &
     &  ,ab,pb, pinit, time_real, &
     &   ptotalmx,ptotalmn,                &
     &   psctotmx,psctotmn,iwetg,iwetg1,iptotal,   &
     &   nscmax,nscmin,psctot1,psctot2,       &
     &   sctot1n,sctot1p,                &
     &   ngscmaxy,ngscminy,ngscxymax,         &
     &   thproc,numproc,    &
     &   rate2d,nrate2d,alpha2d,nalpha2d,    &
     &   axtra,tq1,ntq1,     &
     &   nxi,nyj,nzk            &
     & ,elec,its,ids,ide,jds,jde &
     & )

       USE GRID_MODULE
       USE INDEX_MODULE
       USE MICRO_MODULE
       USE ELEC_MODULE
       USE CPUTIME_MODULE
       USE PARAM_MODULE, only: bcx,bcy, rd, cp, cap => rcp, poo => p00, rdorv => epsilon
       USE COMMASMPI_MODULE, only: ixbeg,ixend,jybeg,jyend,kzbeg,kzend, &
                                   nxbeg,nxend,nybeg,nyend,nzbeg,nzend, &
                                   ktile,my_rank,commasmpi_abort,myprock,nprock ! itile,jtile,
       USE RUN_ATT_NML, only: thistory
       USE TRAJ_MODULE, only: microrates,nmicrorates, elecrates,nelecrates
        use kind_parameters, only: i4=>int4,r8=>real8
   use takcommon,       only: lmax,lfmax,lfmax1,kfmax,iimax,kimax,lr100,ls250,shedsmall,dj, &
     &                        bka,capw,bmu,hlat
   use comm1,           only: rhof,acomm1
   use comm2,           only: sr,sx,srf,vw,bew,szr,szf,szi      &
     &                            ,sri,sfi,di,sxi,sh,ff,vf,fi,vi   &
     &                            ,sxfmlr,isxfmlr,sxfinv,sxfmlrinv,szr,szf &
     &                            ,qv0tak=>qv0,th0,t0tak=>t0,p0,th0pr
        use comm3,           only: workwp,workwm                       &
     &                            ,termkp,termkm,termfp,termfm                 &
     &                            ,xprf,xprfr,sqrtrefnummks
       use comm4, only: acomm4

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
!  4/27/2009: allows for liquid water to be advected on snow and graupel particles using flag "mixedphase"
!
!  3/14/2007: (APS) added qproc temp to make microphysic process timeseries
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
!  9/14/2007: erm: recalculate vx(lh) after setting xdn(lh) in case xdn was out of allowed range.
!             added parameter rimc3 for minimum rime density.  Default value set at 170. kg/m**3
!             instead of previous use of 100.  (Farley, 1987)
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
!      integer icond 
!      parameter ( icond = 2 )

      integer, parameter :: ng1 = 1

      integer nx,ny,nz,na,nba,nv
      integer nor,norz,istag,jstag,kstag ! ,nht,ngt,igsr
      integer iwrite
      real dtp,dx,dy,dz

      logical, intent(in) :: io_flag

      ! ---------
      ! passed in
      ! ----------
      real ventrn
      real qxmin(lc:lqmx)
!      real cxmin(lc:lqmx)
!      real cwmasn,cwmasx,cwmasn5
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real cno(lc:lhab)
      real cckm,ccne,ccnefac,cnexp,ccne0
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      integer lvol(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
      integer lrain(nraintypes)
      real da0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real dab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real dab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real da1 (lc:lqmx)          ! collection coefficients from Seifert 2005

      real bb  (lc:lqmx)

! for 3-moment collection coefficients
      real :: dab0lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real :: dab1lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005

      real ab(-norz+ng1:nz+norz,na)

      integer iexy(lc:lqmx,lc:lqmx)
      
      real :: time_real
      
      integer :: ntq1
      real tq1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,ntq1)  !qproc temporary

      real dzz(nz),dxx(nx),dyy(ny)          ! dz(k),dx(i),dy(j)
      real z1d(-norz+1:nz+norz,4)

      integer, intent(in) :: numproc
      real, intent(inout) :: thproc(nzend,numproc)

      integer, intent(in) :: nxi,nyj,nzk

      ! ----------
      ! functions
      ! -----------
      real, external ::  delbk, delabk, delvbk
      real, external ::  gaml02, gaml02d500, gaml02d300,polysvp1
      double precision :: gamma_dp, gamma_dpr, gamxinfdp
      real :: gamxinf, gaminterp, gaminterpdebug  ! gamma and gamxinf are functions
      
      
      ! -------
      ! local
      ! --------

      integer, parameter :: wrfchem_flag = 0
      logical, parameter :: has_wetscav = .false.
      
      integer :: imixedphase
      
      real, parameter :: warmonly = 0.0 ! testing parameter, set to 1.0 to reduce to warm-rain physics (ice variables stay zero)
      
      real, parameter :: tfr = 273.15, tfrh = 233.15
      
!      real, parameter :: cp = 1004.0, rd = 287.04 
      
      real, parameter :: cpi = 1./cp 
      
!      real, parameter :: cap = rd/cp, poo = 1.0e+05
      real, parameter :: rovcp = rd/cp

      real, parameter :: eperao  = 8.8592e-12 ! permittivity

      real, parameter :: ec = 1.602e-19 ! fundamental unit of charge
      real, parameter :: eci = 1.0/ec

      real, dimension(1, 1) :: rainprod2d, evapprod2d

      integer :: nrate2d
      real rate2d(nx,ny,nrate2d)

!
! Takahashi lookup table
!
      integer, parameter :: nlwc=30, ntem=31
      real takalu(0:ntem,0:nlwc) 

      real, parameter :: rw = 461.5              ! gas const. for water vapor
      real, parameter :: advisc0 = 1.832e-05     ! reference dynamic viscosity (SMT; see Beard & Pruppacher 71)
      real, parameter :: advisc1 = 1.718e-05     ! dynamic viscosity constant used in thermal conductivity calc
      real, parameter :: tka0 = 2.43e-02         ! reference thermal conductivity
      
      integer, parameter :: eqtset = 0

      real, parameter ::      cv = 717.0             ! specific heat at constant volume - air
      REAL, parameter ::      cvv = 1408.5     ! specific heat of water vapor at constant volume
      REAL, parameter ::      cpl = 4190.0     ! specific heat of liquid water  at constant pressure
      REAL, parameter ::      cpigb = 2106.0

      REAL, parameter :: rho00 = 1.225          ! reference/MSL air density

      real, parameter :: pi = 3.141592653589793
      real, parameter :: piinv = 1./pi
      real, parameter :: c1f3 = 1.0/3.0
      real, parameter :: cwc1 = 6.0/(pi*1000.)

      real, parameter :: ar = 841.99666         ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: br = 0.8               ! rain terminal velocity power law coefficient (LFO)

      integer iptotal ! counts number of times that eqtot is exceeded
      real ptotalmx,ptotalmn
      real psctot1,psctot2,sctot1n,sctot1p
      integer nscmax,nscmin,iwetg, iwetg1

      integer ngscmaxy(lg:lhab,lc:ls), ngscminy(lg:lhab,lc:ls), ngscxymax(lg:lhab,lc:ls)

      integer itile,jtile
      
      real, parameter :: gr = 9.8

      
      integer :: nalpha2d
      real :: alpha2d(-nor+1:nx+nor,-norz+ng1:nz+norz,nalpha2d)

      real, parameter :: tfrdry = 243.15

      logical lrescalelow(lc:lhab)
      real tkediss(-nor+1:nx+nor,-norz+ng1:nz+norz)
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)

      real :: galpharaut
      real :: xvbarmax
      
      integer jyslab,its,ids,ide,jds,jde ! domain boundaries
      integer, intent(in) :: iunit !,iunit0
      real qvex
      integer iraincv, icgxconv
      parameter ( iraincv = 1, icgxconv = 1)
      real ffrz
      real :: ffrzh = 1.0

      real qcitmp,cirdiatmp ! ,qiptmp,qirtmp
      real ccwtmp,ccitmp ! ,ciptmp,cirtmp
      real cpqc,cpci ! ,cpip,cpir
      real cpqc0,cpci0 ! ,cpip0,cpir0
      real scfac ! ,cpip1
      
      double precision dp1
      
      double precision frac, frach, xvfrz, xvbiggsnow
      
      double precision :: timevtcalc
      double precision :: dpt1,dpt2
            
      logical, parameter :: gammacheck = .false.
      integer :: luindex
      double precision :: tmpgam
      logical, parameter :: usegamxinfcnu = .false.
      logical, parameter :: usegamxinf = .false.
      logical, parameter :: usegamxinf2 = .false.
      logical, parameter :: usegamxinf3 = .false.
!      real rar  ! rime accretion rate as calculated from qxacw

! a few vars for time-split fallout      
      real vtmax
      integer n,ndfall
      
      double precision chgneg,chgpos,sctot
      
      real temgtmp

      real pb(-norz+ng1:nz+norz)
      real pinit(-norz+ng1:nz+norz)

      real gz(-norz+1:nz+norz) ! 1./dz
      
      real qimax,xni0,roqi0


      real dv

      real dtptmp
      integer itest,nidx,id1,jd1,kd1
      parameter (itest=1)
      parameter (nidx=10)
      parameter (id1=1,jd1=1,kd1=1)
      integer ierr
      integer iend

      integer ix,kz, il, ic, ir, icp1, irp1, ip1,jp1,kp1
      integer :: jy
      integer i,j,k,i1
      integer kzb,kze
      real slope1, slope2
      real x1, x2, x3, y1
      real eps,eps2
      parameter (eps=1.e-20,eps2=1.e-5)
!
!  Other elec. vars
!
      real  temele
      real  trev
      
      logical ldovol, ishail, ltest, wtest
      logical , parameter :: alp0flag = .false.
!
!
!  wind indicies
!
      integer mu,mv,mw
      parameter (mu=1,mv=2,mw=3)
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
      real dtfac
      parameter ( dtfac = 1.0 )
      integer ido(lc:lqmx)

!      integer iexy(lc:lqmx,lc:lqmx)
!      integer ieswi, ieswir, ieswip, ieswc, ieswr
!      integer ieglsw, iegli, ieglir, ieglip, ieglc, ieglr
!      integer iegmsw, iegmi, iegmir, iegmip, iegmc, iegmr
!      integer ieghsw, ieghi, ieghir, ieghip, ieghc, ieghr
!      integer iefwsw, iefwi, iefwir, iefwip, iefwc, iefwr
!      integer iehwsw, iehwi, iehwir, iehwip, iehwc, iehwr
!      integer iehlsw, iehli, iehlir, iehlip, iehlc, iehlr
!      real delqnsa, delqxsa, delqnsb, delqxsb, delqnia, delqxia
!      real delqnra, delqxra

       real delqnxa(lc:lqmx)
       real delqxxa(lc:lqmx)
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

      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

      real tmp3d(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

! 
!  declarations microphyscs and for gather/scatter
!
      integer nxmpb,nzmpb,nxz
      integer jgs,mgs,ngs,numgs
      parameter (ngs=64) !500)
      integer, parameter :: ngsz = 500
      integer ntt
      parameter (ntt=300)

      real dvmgs(ngs)
      
      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs),kgsm2(ngs)
      integer ncuse
      parameter (ncuse=0)
      integer il0(ngs),il5(ngs),il2(ngs),il3(ngs)
!      integer il1m(ngs),il2m(ngs),il3m(ngs),il4m(ngs),il5m(ngs)
!
      real tdtol,temsav,tfrcbw,tfrcbi
      real, parameter :: thnuc = 235.15
!
!  Ice Multiplication Arrays.
!
      real  fimt1(ngs),fimta(ngs),fimt2(ngs) !,qmul1(ngs),qmul2(ngs)
      real xcwmas
!
!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs),ccin(ngs),cina(ngs),ccna(ngs)
      real cwnccn(ngs)
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
      real :: mwat, mice, dice, mwshed, fwmax, fw, mwcrit, massfactor, tmpdiam
      real, parameter :: shedalp = 3.  ! set 3 for maximum mass diameter (same as area-weighted diameter), 4 for mass-weighted diameter
      real ssmax(ngs)       ! maximum SS experienced by a parcel
      real ssmx
      real dnnet,dqnet
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real bfnu, bfnu0, bfnu1
      parameter ( bfnu0 = (rnu + 2.0)/(rnu + 1.0)  )
      real ventr, ventc
      real volb
      double precision t2s, xdp
      double precision xl2p(ngs),rb(ngs)
      real, parameter :: aa1 = 9.44e15, aa2 = 5.78e3  ! a1 in Ziegler
! snow parameters:
      real, parameter :: cexs = 0.1, cecs = 0.5 
      real, parameter :: rvt = 0.104  ! ratio of collection kernels (Zrnic et al, 1993)
      real, parameter :: kfrag = 1.0e-6 ! rate coefficent for collisional splintering (Schuur & Rutledge 00b)
      real, parameter :: mfrag = 1.0e-10 ! assumed ice fragment mass for collisional splintering (Schuur & Rutledge 00b)
      double precision cautn(ngs), rh(ngs), nh(ngs)
      real ex1, ft, rhoinv(ngs)
      real :: ec0(ngs)
      
      real ac1,bc, c1,d1,e1,f1,p380,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,temp3 ! , sstdy, super
      real :: flim, xmass
      real dw,dwr
      double precision :: tmpz, tmpzmlt
      real :: tmpc
      real ratio, delx, dely
      real dbigg,volt
      real chgtmp,fac,mixedphasefac
      real x,y,y2,del,r,rtmp,alpr
      double precision :: vent1,vent2
      double precision :: g1palp,g4palp
      double precision :: g1palpinf,g4palpinf
      real fqt !charge separation as fn of temperature from Dong and Hallett 1992
      real bs
      real v1, v2
      real d1r, d1i, d1s, e1i
      real c1sw   ! integration factor for snow melting with snu = -0.8
      real, parameter :: vr1mm = 5.23599e-10 ! volume of 1mm diameter sphere (m**3)
      real, parameter :: vr3mm = 5.23599e-10*(3.0/1.)**3   ! volume of a 3 mm diameter sphere (m**3) (Rasmussen et al. 1984b, JAS)
      real, parameter :: vr4p5mm = 5.23599e-10*(4.5/1.)**3 ! volume of 4.5mm diameter sphere (m**3) (Rasmussen et al. 1984b, JAS)
      real vmlt,vshd, vshdgs(ngs,lh:lhab), maxmassfac(lc:lhab)
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
      double precision xvc, xvr
      real mwfac
!      real  es(ngs) ! ss(ngs),
!      real  eis(ngs)

      real rwmasn,rwmasx

      real vgra,vfrz
      parameter ( vgra = 0.523599*(1.0e-3)**3 )
     
!      real, parameter :: epsi = 0.622
!      real, parameter :: d = 0.266
      real :: d, dold, denom,denominv,vth
      double precision :: h1, h2, h3, h4,denomdp, denominvdp
      real r1,qevap ! ,slv
      
      real vr,nrx,chw,g1,qr,z,z1,rdi,alp,xnutmp,xnuc,g1r,rd1,rdia,rmas
      real :: snowmeltmass = 0
      
!      real, parameter :: rhofrz = 900.   ! density of graupel from newly-frozen rain
      real, parameter :: rimedens = 500. ! default rime density

!      real svc(ngs)  !  droplet volume
!
!  contact freezing nucleation
!
      real raero,kaero !assumd aerosol radius, thermal conductivity
      parameter ( raero = 3.e-7, kaero = 5.39e-3 )
      real kb   ! Boltzman constant  J K-1
      parameter (kb = 1.3807e-23)
      
      real knud(ngs),knuda(ngs) !knudsen number and correction factor
      real gtp(ngs)  !G(T,p) = 1/(a' + b')  Cotton 72b
      real dfar(ngs) !aerosol diffusivity
      real fn1(ngs),fn2(ngs),fnft(ngs)
      
      real ccia(ngs)
      real ctfzbd(ngs),ctfzth(ngs),ctfzdi(ngs)
!
!  misc
!
      real ni,nis,nr,d0
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs),dqcitmp(ngs),dqwvtmp(ngs)
      real tempc(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) 
      real temgkm1(ngs), temgkm2(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real elv(ngs),elf(ngs),els(ngs)
      real tsqr(ngs),ssi(ngs),ssw(ngs),ssi0(ngs)
      real qcwtmp(ngs),qtmp,qtot(ngs) 
      real qcond(ngs)
      real ctmp, sctmp
      real cimasn,cimasx,ccimx
      real pid4
      real cs,ds,gf7,gf6,gf5,gf4,gf3,gf2,gf1
      real gcnup1,gcnup2
      real gf73rds, gf83rds
      real gamice73fac, gamsnow73fac
      real gf43rds, gf53rds
      real gamma_sp
      real cnudiag
      real aradcw,bradcw,cradcw,dradcw,cwrad,rwrad,rwradmn
      parameter ( rwradmn = 50.e-6 )
      real dh0
      real dg0(ngs),df0(ngs)
      real dhwet(ngs),dhlwet(ngs),dfwet(ngs)
      
      real clionpmx,clionnmx
      parameter (clionpmx=1.e9,clionnmx=1.e9) ! Takahashi 84
!
!  other arrays

      real fwet1(ngs),fwet2(ngs)   
      real fmlt1(ngs),fmlt2(ngs),fmlt1e(ngs)
      real fvds(ngs),fvce(ngs),fiinit(ngs) 
      real fvent(ngs),fraci(ngs),fracl(ngs)
!
      real fai(ngs),fav(ngs),fbi(ngs),fbv(ngs)
      real felv(ngs),fels(ngs),felf(ngs)
      real felvcp(ngs),felscp(ngs),felfcp(ngs)
      real felvpi(ngs),felspi(ngs),felfpi(ngs)
      real felvs(ngs),felss(ngs)      !   ,felfs(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs) ! heat capacities of ice and liquid
      real fschm(ngs),fpndl(ngs)
      real fgamw(ngs),fgams(ngs)
      real fcqv1(ngs),fcqv2(ngs),fcc3(ngs) 
      
      real cvm,cpm,rmm

      real, parameter ::      cpv = 1885.0       ! specific heat of water vapor at constant pressure
!
      real fcci(ngs), fcip(ngs)
!
      real :: sfm1(ngs),sfm2(ngs)
      real :: gfm1(ngs),gfm2(ngs)
      real :: ffm1(ngs),ffm2(ngs)
      real :: hfm1(ngs),hfm2(ngs)

      logical :: wetsfc(ngs),wetsfchl(ngs),wetsfcf(ngs)
      logical :: wetgrowth(ngs), wetgrowthhl(ngs), wetgrowthf(ngs)

      real qitmp(ngs),qistmp(ngs)
       
      real rzxh(ngs), rzxhl(ngs), rzxhlh(ngs), rzxhlf(ngs)
      real rzxs(ngs), rzxf(ngs)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)
      real cdh(ngs),cdhl(ngs)
      real :: axx(ngs,lh:lhab),bxx(ngs,lh:lhab)
      real vt2ave(ngs)

      real :: qcwresv(ngs), ccwresv(ngs) ! "reserved" droplet mass and number that are too small for accretion
      
      real ::  lfsave(ngs,6)
      real ::  qx(ngs,lv:lhab)
      real ::  qxw(ngs,ls:lhab)
      real ::  qxwlg(ngs,lh:lhab)
      real ::  chxf(ngs,lh:lhab)
      real ::  cx(ngs,lc:lhab)
      real ::  cxmxd(ngs,lc:lhab)
      real ::  qxmxd(ngs,lv:lhab)
      real ::  scx(ngs,lc:lhab)
      real ::  xv(ngs,lc:lhab)
      real ::  xsfca(ngs,lc:lhab)
      real ::  qxrain(ngs,1:nraintypes)
      real ::  qxrainold(ngs,0:nraintypes)
      real ::  pqmelt, pqauto, pqshed
      real ::  qxshedfrac(ngs,lr:lhab,0:nraintypes)
      real ::  vtxbar(ngs,lc:lhab,3)
      real ::  xmas(ngs,lc:lhab)
      real ::  xdn(ngs,lc:lhab)
      real ::  xdntmp(ngs,lc:lhab)
      real ::  cdxgs(ngs,lc:lhab)
      real ::  xdia(ngs,lc:lhab,3)
      real ::  vtwtdia(ngs,lr:lhab) ! sweep-out volume weighted diameter
      real ::  rarx(ngs,ls:lhab)
      real ::  vx(ngs,li:lhab)
      real ::  rimdn(ngs,li:lhab)
      real ::  raindn(ngs,li:lhab)
      real ::  alpha(ngs,lc:lhab)
      real ::  dab0lh(ngs,lc:lhab,lc:lhab)
      real ::  dab1lh(ngs,lc:lhab,lc:lhab)
      real ::  zx(ngs,lr:lhab)
      real ::  zxmxd(ngs,lr:lhab)
      real ::  g1x(ngs,lr:lhab)

      real :: g1xmax,g1xmin
      real :: qsimxdep(ngs) ! max sublimation of qi+qs+qis
      real :: qsimxsub(ngs) ! max depositionof qi+qs+qis
      logical,parameter :: DoSublimationFix = .true.
      real :: qrtmp(ngs),qvtmp(ngs),qctmp(ngs)
      real :: felvcptmp,felscptmp,qsstmp
      real :: thetatmp, thetaptmp, temcgtmp,qvaptmp
      real :: qvstmp, qisstmp, qvptmp, qitmp1, qctmp1
      
      real :: galphrout
      
      real ventrx(ngs)
      real ventrxn(ngs)
      real g1shr, alphashr
      real g1mlr, alphamlr
      real g1smlr, alphasmlr
      real massfacshr, massfacmlr
      
      real :: qhgt8mm ! ice mass greater than 8mm
      real :: qhwgt8mm ! ice + max water mass greater than 8mm
      real :: qhgt10mm ! mass greater than 10mm
      real :: qhgt20mm ! mass greater than 20mm
      real :: fwmhtmp
!      real, parameter :: fwmhtmptem = -15. ! temperature at which fwmhtmp fully switches to liquid water only being on large particles
      real, parameter :: d1t = (6.0 * 0.268e-3/(917.* pi))**(1./3.) ! d1t is the diameter of the ice sphere with the mass (0.268e-3 kg) of an 8mm spherical drop
      real, parameter :: srasheym = 0.1389 ! slope fraction from Rasmussen and Heymsfield 
      real :: dtmp
!
      real swvent(ngs),hwvent(ngs),rwvent(ngs),hlvent(ngs),hwventy(ngs),hlventy(ngs),rwventz(ngs)
      real hxventtmp
      real hlventinc(ngs),hwventinc(ngs)
      integer, parameter :: ndiam = 10
      integer :: numdiam
      real hwvent0(ndiam+4),hlvent0 ! 0 to d1
      real hwvent1,hlvent1 ! d1 to infinity
      real hwvent2,hlvent2 ! d2 to infinity
      real gama0,gamb0
      real gama1,gamb1
      real gama2,gamb2
!      real, parameter :: mltdiam1 = 9.0e-3, mltdiam1p5 = 16.0e-3, mltdiam2 = 19.0e-3, mltdiam3 = 200.0e-3, mltdiam05 = 4.5e-3
      real :: mltdiam(ndiam+4)
      real mltmass0inv,mltmass1inv,mltmass2inv, mltmass1cgs, mltmass2cgs,mltmass3inv, mltmass3cgs
      real qhmlr0, qhmlr05, qhmlr1, qhmlr2,qhmlr3, qhmlr12, qhmlr23
      real qhlmlr0, qhlmlr05, qhlmlr1, qhlmlr2,qhlmlr3, qhlmlr12, qhlmlr23
      real qxd1, cxd1, zxd1 ! mass and number up to mltdiam1
      real qxd05, cxd05 ! mass and number up to mltdiam1/2
      real :: qrbreak, crbreaksmall, crbreak, zrbreak, breakbin
      
      real :: qxd(ndiam+4), cxd(ndiam+4), qhml(ndiam+4), qhml0(ndiam+4)
      real :: dqxd(ndiam+4), dcxd(ndiam+4), dqhml(ndiam+4)
      
      
      real civent(ngs)
      real isvent(ngs)
!
      real xmascw(ngs)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real dnmx
      real :: xdiamxmas(ngs,lc:lhab)
!
      real cilen(ngs) ! ,ciplen(ngs)
!
!
      real rwcap(ngs),swcap(ngs)
      real hwcap(ngs)
      real hlcap(ngs)
      real cicap(ngs)
      real iscap(ngs)

      real qvimxd(ngs)
      real qimxd(ngs),qismxd(ngs),qcmxd(ngs),qrmxd(ngs),qsmxd(ngs),qhmxd(ngs),qhlmxd(ngs)
      real cimxd(ngs),ccmxd(ngs),crmxd(ngs),csmxd(ngs),chmxd(ngs)
      real qfmxd(ngs), cfmxd(ngs)
      real cionpmxd(ngs),cionnmxd(ngs)
      real clionpmxd(ngs),clionnmxd(ngs)

!
!
!  Da electricity array work space...
!
!
      real dezcomp(ngs)
      real sctem(ngs)
      real dellwc(ngs)
      real,parameter :: tcc = -30.
!
!  space charge arrays.....
!
      real scsacw(ngs)
      real scsacr(ngs)
      real sciacw(ngs) ! not present...
      real scsaci(ngs)
      real :: scsacis(ngs)
      real sxsaci(ngs)
!
      real schacw(ngs), schacr(ngs)
      real schaci(ngs), schacs(ngs)
      real :: schacis(ngs)
      real sxhaci(ngs), sxhacs(ngs)

      real :: scfmlr(ngs), scfshr(ngs)
      real sxfaci(ngs), sxfacs(ngs)

      real schmlr(ngs), scsmlr(ngs)
      real scsshr(ngs), schshr(ngs)
      real scscev(ngs), scsdep(ngs)
!
      real schlacw(ngs), schlacr(ngs)
      real schlaci(ngs), schlacs(ngs)
      real :: schlacis(ngs)

      real scfacw(ngs), scfacr(ngs)
      real scfaci(ngs), scfacs(ngs)
      real :: scfacis(ngs)

      real exy(ngs,ls:lhab,lc:ls)
      real scxacy(ngs,ls:lhab,lc:ls)
      real cxacy(ngs,ls:lhab,lc:ls)
      real sxxacy(ngs,ls:lhab,li:ls)

! ions  - one size assumed

      real cionp(ngs),cionn(ngs),clionp(ngs),clionn(ngs) 

!    charge concentration

!    production of space charge per collison/per unit mass?

      real fschw(ngs),fscsw(ngs),fschl(ngs),fscfw(ngs)
      real fschw2(ngs),fscsw2(ngs),fschl2(ngs),fscfw2(ngs)
      real fsccw(ngs),fscci(ngs),fscrw(ngs),fscis(ngs)

!    production of space charge term  + / -

      real  psctot(ngs)
      real  psctotmx,psctotmn
      real  psccwi(ngs), psccwd(ngs) ! cloud water
      real  psccii(ngs), psccid(ngs) ! ice crystals
      real  pscisi(ngs), pscisd(ngs) ! ice spheres
      real  pscrwi(ngs), pscrwd(ngs)  ! rain water
      real  pscswi(ngs), pscswd(ngs)  ! snow
      real  pschwi(ngs), pschwd(ngs)  ! graupel
      real  pscfwi(ngs), pscfwd(ngs)  ! frozen drops
      real  pschli(ngs), pschld(ngs)   ! hail

!ions assume mass = 0
      real  pscpii(ngs), pscpid(ngs) !positive ions
      real  pscnii(ngs), pscnid(ngs) !negative ions
      real  pscplii(ngs), pscplid(ngs) !positive (large) ions
      real  pscnlii(ngs), pscnlid(ngs) !negative (large) ions

!      production of space charge from mass transfer

      real  psccwmi(ngs), psccwmd(ngs)
      real  psccimi(ngs), psccimd(ngs)
      real  pscrwmi(ngs), pscrwmd(ngs)
      real  pscswmi(ngs), pscswmd(ngs)
      real  pschwmi(ngs), pschwmd(ngs)
      real  pscfwmi(ngs), pscfwmd(ngs)
      real  pschlmi(ngs), pschlmd(ngs)


      TYPE(VARIABLE)     :: elec(neelec)

!
!
      ! Hallett-Mossop arrays
      real chmul1(ngs),chlmul1(ngs),csmul1(ngs),csmul(ngs)
      real qhmul1(ngs),qhlmul1(ngs),qsmul1(ngs),qsmul(ngs)
      
      ! splinters from drop freezing
      real csplinter(ngs),qsplinter(ngs)
      real csplinter2(ngs),qsplinter2(ngs)
!
!
!  concentration arrays...
!
      real :: chlcnh(ngs), vhlcnh(ngs), vhlcnhl(ngs)
      real :: chlcnhhl(ngs) ! number of new hail particles (may be different from number of lost graupel)
      real cracif(ngs), ciacrf(ngs)
      real cracr(ngs)

!
      real ciint(ngs), crfrz(ngs), crfrzf(ngs), crfrzs(ngs)
      real cicint(ngs)
      real cipint(ngs)
      real ciacw(ngs), cwacii(ngs) 
      real ciacr(ngs), craci(ngs)
      real csacw(ngs)
      real csacr(ngs)
      real csaci(ngs),   csacs(ngs)
      real cracw(ngs) 
      real chacw(ngs), chacr(ngs)
      real :: chlacw(ngs) 
      real chaci(ngs), chacs(ngs)
!
      real :: chlacr(ngs)
      real :: chlaci(ngs), chlacs(ngs)
      real crcnw(ngs) 
      real cidpv(ngs),cisbv(ngs)
      real cisdpv(ngs),cissbv(ngs)
      real cimlr(ngs),cismlr(ngs)

      real chlsbv(ngs), chldpv(ngs)
      real chlmlr(ngs), chlmlrr(ngs)
      real chlfmlr(ngs)
!      real chlmlrsave(ngs),chlsave(ngs),qhlsave(ngs)
      real chlshr(ngs), chlshrr(ngs)

      real dxshdrate(ngs,0:3)

      real chdpv(ngs),chsbv(ngs)
      real chmlr(ngs),chcev(ngs)
      real chmlrr(ngs)
      real chshr(ngs), chshrr(ngs)

      real csdpv(ngs),cssbv(ngs)
      real csmlr(ngs),csmlrr(ngs),cscev(ngs)
      real csshr(ngs), csshrr(ngs)

      real crcev(ngs)
      real crshr(ngs)
      real cwshw(ngs), qwshw(ngs)
!
!
! arrays for w-ac-x ;  x-ac-w
!
!
!
      real qrcnw(ngs), qwcnr(ngs)
      real zrcnw(ngs),zracr(ngs),zracw(ngs),zrcev(ngs)

      real qracw(ngs) ! qwacr(ngs),
      real qiacw(ngs) !, qwaci(ngs)

      real qsacw(ngs) ! ,qwacs(ngs),
      real qhacw(ngs) ! qwach(ngs),
      real :: qhlacw(ngs), qxacwtmp, qxacrtmp, qxacitmp, qxacstmp !
      real :: cxacstmp,cxacitmp
      real vhacw(ngs), vsacw(ngs), vhlacw(ngs), vhlacr(ngs)

      real qwacf(ngs),qfacw(ngs)
      real qfacr(ngs), qfacrmlr(ngs),qfacwmlr(ngs)
!      real qfacrf(ngs) ! qracff(ngs),
      real qfaci(ngs), qfaci0(ngs) ! ,qiacf(ngs)
      real qfacis(ngs), qfacis0(ngs) ! ,qiacf(ngs)
      real qfacip(ngs) ! ,qipacf(ngs)
      real qfacs(ngs), qfacs0(ngs) ! ,qsacf(ngs)
!      real qfdpv(ngs),qfsbv(ngs) ! qfcnv(ngs),qfevv(ngs),
!      real qfmlr(ngs),qfdsv(ngs) ! ,
!      real qfwet(ngs),qfdry(ngs),qfshr(ngs)
!      real pqfwi(ngs),pqfwd(ngs)
      real qhlcnf(ngs)

      real cfacw(ngs)
      real cfacr(ngs) 
      real cfaci(ngs),cfaci0(ngs)
      real cfacis(ngs),cfacis0(ngs)
      real cfacs(ngs), cfacs0(ngs)
      real chlcnf(ngs),chlcnfhl(ngs),fddenfac(ngs)
      real cfcev(ngs)
      
      real vfacw(ngs), vfacr(ngs)
      real zfacr(ngs), zracf(ngs)
      
      real fwvent(ngs),fwventy(ngs), fwcap(ngs)
      real fwventinc(ngs)

      real cfmlr(ngs),cfsbv(ngs)
      real cfmlrr(ngs)
      
      real cfshr(ngs),cfshrr(ngs)

      real qfcev(ngs)
      real qfmul1(ngs),cfmul1(ngs)
!
      real qsacws(ngs)

!
!  arrays for x-ac-r and r-ac-x; 
!
      real qsacr(ngs),qracs(ngs)
      real qhacr(ngs),qhacrmlr(ngs),qhacwmlr(ngs) ! ,qrach(ngs)
      real vhacr(ngs), zhacr(ngs), zhacrf(ngs), zrach(ngs), zrachl(ngs)
      real qiacr(ngs),qraci(ngs)
      
      real ziacr(ngs)

      real qracif(ngs),qiacrf(ngs),qiacrs(ngs),ciacrs(ngs)

      real :: qhlacr(ngs),qhlacrmlr(ngs), qhlacwmlr(ngs)
      real qsacrs(ngs) !,qracss(ngs)
!
!  ice - ice interactions
!
      real qsaci(ngs)
      real qsacis(ngs)
      real csacis(ngs)
      real qhaci(ngs)
      real qhacs(ngs)

      real :: qhacis(ngs) 
      real :: chacis(ngs) 
      real :: chacis0(ngs)

      real :: csaci0(ngs) ! collision rate only
      real :: csacis0(ngs) ! collision rate only
      real :: chaci0(ngs) ! collision rate only
      real :: chacs0(ngs) ! collision rate only
      real :: chlaci0(ngs)
      real :: chlacis(ngs)
      real :: chlacis0(ngs)
      real :: chlacs0(ngs) 

      real :: qsaci0(ngs) ! collision rate only
      real :: qsacis0(ngs) ! collision rate only
      real :: qhaci0(ngs) ! collision rate only
      real :: qhacis0(ngs) ! collision rate only
      real :: qhacs0(ngs) ! collision rate only
      real :: qhlaci0(ngs)  
      real :: qhlacis0(ngs)
      real :: qhlacs0(ngs) 

      real :: qhlaci(ngs)  
      real :: qhlacis(ngs)
      real :: qhlacs(ngs)
!
!  conversions
!
      real qrfrz(ngs) ! , qirirhr(ngs)
      real zrfrz(ngs), zrfrzf(ngs), zrfrzs(ngs)
      real ziacrf(ngs), zhcnsh(ngs), zhcnih(ngs)
      real zhacw(ngs), zhacs(ngs), zhaci(ngs)
      real zhmlr(ngs), zhdsv(ngs), zhsbv(ngs), zhlcnh(ngs), zhshr(ngs)
      real zfacw(ngs), zfacs(ngs), zfaci(ngs)
      real zfmlr(ngs), zfdsv(ngs), zfsbv(ngs), zhlcnf(ngs), zfshr(ngs), zfshrr(ngs)
      real zhmlrtmp,zhmlr0inf,zhlmlr0inf
      real zhmlrr(ngs),zhlmlrr(ngs),zhshrr(ngs),zhlshrr(ngs),zfmlrr(ngs)
!      real zsmlr(ngs)
      real zsmlrr(ngs), zsshr(ngs), zsshrr(ngs)
      real zhcns(ngs), zhcni(ngs)
      real zhwdn(ngs), zfwdn(ngs) ! change in Z due to density changes
      real zhldn(ngs) ! change in Z due to density changes

      real zhlacw(ngs), zhlacs(ngs), zhlacr(ngs)
      real zhlmlr(ngs), zhldsv(ngs), zhlsbv(ngs), zhlshr(ngs)

      
      real vrfrzf(ngs), viacrf(ngs)
      real qrfrzs(ngs), qrfrzf(ngs)
      real qwfrz(ngs), qwctfz(ngs)
      real cwfrz(ngs), cwctfz(ngs)
      real qwfrzis(ngs), qwctfzis(ngs) ! droplet freezing to ice spheres
      real cwfrzis(ngs), cwctfzis(ngs)
      real qwfrzc(ngs), qwctfzc(ngs) ! droplet freezing to columns
      real cwfrzc(ngs), cwctfzc(ngs)
      real qwfrzp(ngs), qwctfzp(ngs) ! droplet freezing to plates
      real cwfrzp(ngs), cwctfzp(ngs)
      real xcolmn(ngs), xplate(ngs)
      real ciihr(ngs), qiihr(ngs)
      real cicichr(ngs), qicichr(ngs)
      real cipiphr(ngs), qipiphr(ngs)
      real qscni(ngs), cscni(ngs), cscnis(ngs)
      real qscnvi(ngs), cscnvi(ngs), cscnvis(ngs)
      real qhcns(ngs), chcns(ngs), chcnsh(ngs), vhcns(ngs)
      real qscnh(ngs), cscnh(ngs), vscnh(ngs)
      real qhcni(ngs), chcni(ngs), chcnih(ngs), vhcni(ngs)
      real qiint(ngs),qipipnt(ngs),qicicnt(ngs)
      real cninm(ngs),cnina(ngs),cninp(ngs),wvel(ngs),wvelkm1(ngs)
      real tke(ngs)
      real uvel(ngs),vvel(ngs)
!
      real qidpv(ngs),qisbv(ngs) ! qicnv(ngs),qievv(ngs),
      real qimlr(ngs),qidsv(ngs),qisdsv(ngs),qidsvp(ngs) ! ,qicev(ngs)
      real qismlr(ngs)

!
      real qfdpv(ngs),qfsbv(ngs) ! qfcnv(ngs),qfevv(ngs),
      real qfmlr(ngs),qfdsv(ngs) ! ,qfcev(ngs)
      real qfwet(ngs),qfdry(ngs),qfshr(ngs)
      real qfshrp(ngs)
      real vfsoak(ngs), vfmlr(ngs), vfshdr(ngs), vhlcnf(ngs), vhlcnfhl(ngs)
      real cfdpv(ngs)
!
      real :: qhldpv(ngs), qhlsbv(ngs) ! qhlcnv(ngs),qhlevv(ngs),
      real :: qhlmlr(ngs), qhldsv(ngs), qhlmlrsave(ngs)
      real :: qhlwet(ngs), qhldry(ngs), qhlshr(ngs), qxwettmp
!
      real :: qrfz(ngs),qsfz(ngs),qhfz(ngs),qhlfz(ngs)
      real :: qffz(ngs)
!
      real qhdpv(ngs),qhsbv(ngs) ! qhcnv(ngs),qhevv(ngs),
      real qhmlr(ngs),qhdsv(ngs),qhcev(ngs),qhcndv(ngs),qhevv(ngs)
      real qhlcev(ngs), chlcev(ngs)
      real qhwet(ngs),qhdry(ngs),qhshr(ngs)
      real qhshrp(ngs)
      real qhshh(ngs) !accreted water that remains on graupel
      real qhmlh(ngs) !melt water that remains on graupel
      real qhfzh(ngs) !water that freezes on mixed-phase graupel
      real qffzf(ngs) !water that freezes on mixed-phase FD
      real qhlfzhl(ngs) !water that freezes on mixed-phase hail
      
      real qhmlrlg(ngs),qhlmlrlg(ngs) ! melting from the larger diameters
      real qhfzhlg(ngs) !water that freezes on mixed-phase graupel (large sizes)
      real qhlfzhllg(ngs) !water that freezes on mixed-phase hail (large sizes)
      real qhlcevlg(ngs), chlcevlg(ngs)
      real qhcevlg(ngs), chcevlg(ngs)

      real vhfzh(ngs), vffzf(ngs) ! change in volume from water that freezes on mixed-phase graupel, frozen drops
      real vhlfzhl(ngs) !  change in volume from water that freezes on mixed-phase hail

      real vhshdr(ngs) !accreted water that leaves on graupel (mixedphase)
      real vhlshdr(ngs) !accreted water that leaves on hail (mixedphase)
      real vhmlr(ngs) !melt water that leaves graupel (single phase)
      real vhlmlr(ngs) !melt water that leaves hail (single phase)
      real vhsoak(ngs) !  aquired water that seeps into graupel.
      real vhlsoak(ngs) !  aquired water that seeps into hail.
      
!
      real qsdpv(ngs),qssbv(ngs) ! qscnv(ngs),qsevv(ngs),
      real qsmlr(ngs),qsdsv(ngs),qscev(ngs),qscndv(ngs),qsevv(ngs)
      real qswet(ngs),qsdry(ngs),qsshr(ngs)
      real qsshrp(ngs)
      real qsfzs(ngs)
!
!
      real qipdpv(ngs),qipsbv(ngs)
      real qipmlr(ngs),qipdsv(ngs)
!
      real qirdpv(ngs),qirsbv(ngs)
      real qirmlr(ngs),qirdsv(ngs),qirmlw(ngs)
!
      real qgldpv(ngs),qglsbv(ngs)
      real qglmlr(ngs),qgldsv(ngs)
      real qglwet(ngs),qgldry(ngs),qglshr(ngs)
      real qglshrp(ngs)
!
      real qgmdpv(ngs),qgmsbv(ngs)
      real qgmmlr(ngs),qgmdsv(ngs)
      real qgmwet(ngs),qgmdry(ngs),qgmshr(ngs)
      real qgmshrp(ngs)
      real qghdpv(ngs),qghsbv(ngs)
      real qghmlr(ngs),qghdsv(ngs) 
      real qghwet(ngs),qghdry(ngs),qghshr(ngs)
      real qghshrp(ngs)
!
      real qrztot(ngs),qrzmax(ngs),qrzfac(ngs)
      real qrcev(ngs)
      real qrshr(ngs)
      real fsw(ngs),fhw(ngs),fhlw(ngs),ffw(ngs) !liquid water fractions
      real fswmax(ngs),fhwmax(ngs),fhlwmax(ngs) !liquid water fractions
      real ffwmax(ngs)
      real qhcnf(ngs) 
      real :: qhlcnh(ngs)
      real qhcngh(ngs),qhcngm(ngs),qhcngl(ngs)
      
      real :: qhcnhl(ngs), chcnhl(ngs), zhcnhl(ngs), vhcnhl(ngs) ! conversion of low-density hail back to graupel

      real eiw(ngs),eii(ngs),eiri(ngs),eipir(ngs),eisw(ngs)
      real erw(ngs),esw(ngs),eglw(ngs),eghw(ngs),efw(ngs)
      real ehxw(ngs),ehlw(ngs),egmw(ngs),ehw(ngs)
      real err(ngs),esr(ngs),eglr(ngs),eghr(ngs),efr(ngs)
      real ehxr(ngs),ehlr(ngs),egmr(ngs) 
      real eri(ngs),esi(ngs),esis(ngs),egli(ngs),eghi(ngs),efi(ngs),efis(ngs)
      real ehxi(ngs),ehli(ngs),egmi(ngs),ehi(ngs),ehis(ngs),ehlis(ngs)
      real ers(ngs),ess(ngs),egls(ngs),eghs(ngs),efs(ngs),ehs(ngs),ehsfac(ngs)
      real ehscnv(ngs)
      real ehxs(ngs),ehls(ngs),egms(ngs),egmip(ngs) 

      real ehsclsn(ngs),ehiclsn(ngs),ehisclsn(ngs)
      real efsclsn(ngs),eficlsn(ngs),efisclsn(ngs)
      real ehlsclsn(ngs),ehliclsn(ngs),ehlisclsn(ngs)
      real esiclsn(ngs),esisclsn(ngs)

      real :: ehs_collsn = 0.5, ehi_collsn = 1.0
      real :: efs_collsn = 0.5, efi_collsn = 1.0
      real :: ehls_collsn = 1.0, ehli_collsn = 1.0
      real :: esi_collsn = 1.0
      
      real ew(8,6)
      real cwr(8,2)  ! radius and inverse of interval
      data cwr / 2.0, 3.0, 4.0, 6.0,  8.0,  10.0, 15.0,  20.0 , & ! radius
     &           1.0, 1.0, 0.5, 0.5,  0.5,   0.2,  0.2,  1.  /   ! inverse of interval
      integer icwr(ngs), igwr(ngs), irwr(ngs), ihlr(ngs), ifwr(ngs)
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


      real da0lr(ngs),da1lr(ngs)
      real da0lc(ngs),da1lc(ngs)
      real da0lh(ngs)
      real da0lhl(ngs)
      real da0lf(ngs)
      real :: da0lx(ngs,lr:lhab)
      
      real va0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real vab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real vab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real va1 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real ehip(ngs),ehlip(ngs),ehlir(ngs)
      real erir(ngs),esir(ngs),eglir(ngs),egmir(ngs),eghir(ngs)
      real efir(ngs),ehir(ngs),eirw(ngs),eirir(ngs),ehr(ngs)
      real erip(ngs),esip(ngs),eglip(ngs),eghip(ngs)
      real efip(ngs),eipi(ngs),eipw(ngs),eipip(ngs)
!
!  arrays for production terms
!
      real ptotal(ngs) ! , pqtot(ngs)
!
      real pqcwi(ngs),pqcii(ngs),pqrwi(ngs),pqisi(ngs)
      real pqswi(ngs),pqhwi(ngs),pqwvi(ngs)
      real pqgli(ngs),pqghi(ngs),pqfwi(ngs)
      real pqgmi(ngs),pqhli(ngs) ! ,pqhxi(ngs)
      real pqiri(ngs),pqipi(ngs) ! pqwai(ngs),
      real pqlwsi(ngs),pqlwhi(ngs),pqlwhli(ngs),pqlwfi(ngs)
      
      real pqlwlghi(ngs),pqlwlghli(ngs)
      real pqlwlghd(ngs),pqlwlghld(ngs)
      
      real pqrauto(ngs), pqrshed(ngs), pqrmelt(ngs), pqrother(ngs)
      
      

      real pvhwi(ngs), pvhwd(ngs)
      real pvfwi(ngs), pvfwd(ngs)
      real pvhli(ngs), pvhld(ngs)
      real pvswi(ngs), pvswd(ngs)
!
      real pqcwd(ngs),pqcid(ngs),pqrwd(ngs),pqisd(ngs), pqcwdacc(ngs)
      real pqswd(ngs),pqhwd(ngs),pqwvd(ngs)
      real pqgld(ngs),pqghd(ngs),pqfwd(ngs)
      real pqgmd(ngs),pqhld(ngs) ! ,pqhxd(ngs)
      real pqird(ngs),pqipd(ngs) ! pqwad(ngs),
      real pqlwsd(ngs),pqlwhd(ngs),pqlwhld(ngs),pqlwfd(ngs)
!
!      real pqxii(ngs,nhab),pqxid(ngs,nhab)
!
      real  pctot(ngs)
      real  pcipi(ngs), pcipd(ngs)
      real  pciri(ngs), pcird(ngs)
      real  pccwi(ngs), pccwd(ngs), pccwdacc(ngs)
      real  pccii(ngs), pccid(ngs)
      real  pcisi(ngs), pcisd(ngs)
      real  pccin(ngs)
      real  pcrwi(ngs), pcrwd(ngs)
      real  pcswi(ngs), pcswd(ngs)
      real  pchwi(ngs), pchwd(ngs)
      real  pchli(ngs), pchld(ngs)
      real  pcfwi(ngs), pcfwd(ngs)
      real  pcgli(ngs), pcgld(ngs)
      real  pcgmi(ngs), pcgmd(ngs)
      real  pcghi(ngs), pcghd(ngs)

      real  pzrwi(ngs), pzrwd(ngs)
      real  pzhwi(ngs), pzhwd(ngs)
      real  pzfwi(ngs), pzfwd(ngs)
      real  pzhli(ngs), pzhld(ngs)
      real  pzswi(ngs), pzswd(ngs)

!
!  other arrays
!
      real dqisdt(ngs) !,advisc(ngs) !dqwsdt(ngs), ,schm(ngs),pndl(ngs)

      real qss0(ngs)

      real qsacip(ngs)
      real pres(ngs),pipert(ngs)
      real pk(ngs)
      real rho0(ngs),pi0(ngs)
      real rhovt(ngs),sqrtrhovt
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
      real ptwfzi(ngs),ptimlw(ngs)
      real psub(ngs),pvap(ngs),pfrz(ngs),ptem(ngs),pmlt(ngs),pevap(ngs),pdep(ngs),ptem2(ngs)
      
      real cnostmp(ngs)   ! for diagnosed snow intercept
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
      parameter (iholef = 1)
      parameter (iholen = 1)
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
!   Miscellaneous variables
!
      real, parameter :: cwmas30 = 1000.*0.523599*(2.*30.e-6)**3 ! mass of 30-micron radius droplet, for sat. adj.
      real, parameter :: cwmas20 = 1000.*0.523599*(2.*20.e-6)**3 ! mass of 20-micron radius droplet, for sat. adj.
      integer ireadqf,lrho,lqsw,lqgl,lqgm ,lqgh 
      integer lqrw
      real vt
      real arg  ! gamma is a function
      real erbnd1, fdgt1, costhe1
      real qeps
      real dyi2,dzi2,bta1,cnit,dragh,dnz00,pii ! ,cp608
      real qccrit,gf4br,gf4ds,gf4p5, gf3ds, gf1ds
      real gf1palp(ngs) ! for storing Gamma[1.0 + alphar]

      
      real xdn0(lc:lhab)
      real xdn_new,drhodt
      
      integer l ,ltemq,inumgs, idelq

      real brz,arz,temq

      real ssival,tqvcon
      real cdx(lc:lhab)
      real cnox
      real cval,aval,eval,fval,gval ,qsign,ftelwc,qconkq,elecfac,altelecfac
      real qconm,qconn,cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5,qdiff,argrcnw
      real c4,bradp,bl2,bt2,dthr,hrifac, hdia0,hdia1,civenta,civentb
      real civentc,civentd,civente,civentf,civentg,cireyn,xcivent
      real cipventa,cipventb,cipventc,cipventd,cipreyn,cirventa
      real cirventb
      integer igmrwa,igmrwb,igmswa, igmswb,igmfwa,igmfwb,igmhwa,igmhwb
      real rwventa ,rwventb,swventa,swventb,fwventa,fwventb,fwventc
      real hwventa,hwventb
      real    hwventc, hlventa, hlventb,  hlventc
      real  glventa, glventb, glventc
      real   gmventa, gmventb,  gmventc, ghventa, ghventb, ghventc
      real  dzfacp,  dzfacm,  cmassin,  cwdiar 
      real  rimmas, rhobar
      real   argtim, argqcw, argqxw, argtem
      real   frcswsw, frcswgl, frcswgm, frcswgh, frcswfw, frcswsw1
      real   frcglgl, frcglgm, frcglgh,  frcglfw, frcglgl1
      real   frcgmgl, frcgmgm, frcgmgh,  frcgmfw, frcgmgm1
      real   frcghgl, frcghgm, frcghgh,  frcghfw,  frcghgh1
      real   frcfwgl, frcfwgm, frcfwgh, frcfwfw,  frcfwfw1
      real   frcswrsw, frcswrgl,  frcswrgm,  frcswrgh, frcswrfw
      real   frcswrsw1
      real   frcrswsw, frcrswgl, frcrswgm, frcrswgh, frcrswfw
      real  frcrswsw1
      real  frcglrgl, frcglrgm, frcglrgh,  frcglrfw, frcglrgl1
      real  frcrglgl
      real  frcrglgm,  frcrglgh, frcrglfw, frcrglgl1
      real  frcgmrgl, frcgmrgm, frcgmrgh, frcgmrfw,  frcgmrgm1
      real  frcrgmgl, frcrgmgm,  frcrgmgh, frcrgmfw, frcrgmgm1
      real  total,  qweps,  gf2a, gf4a, dqldt, dqidt, dqdt
      real frcghrgl, frcghrgm, frcghrgh, frcghrfw, frcghrgh1, frcrghgl
      real frcrghgm, frcrghgh,  frcrghfw, frcrghgh1
      real    a1,a2,a3,a4,a5,a6
      real   gamss
      real cdw, cdi, denom1, denom2, delqci1, delqip1
      real cirtotn,  ciptotn, cgmtotn, chltotn,  cirtotp
      real  cgmfac, chlfac,  cirfac
      integer igmhla, igmhlb, igmgla, igmglb, igmgma,  igmgmb
      integer igmgha, igmghb
      integer idqis, item, itim0 
      integer  iqgl, iqgm, iqgh, iqrw, iqsw 
      integer  itertd, ia
      
      integer :: infdo
      
      real tau, ewtmp
      real ftau
      
      integer cntnic_noliq
      real     q_noliqmn, q_noliqmx
      real     scsacimn, scsacimx
      
      real :: dtpinv
      
!   arrays for temporary bin space

      integer nbin
      parameter (nbin=ntakpd)  ! number of mass bins for bin model
      real rn(nbin) !,rd(nbin),rm(nbin)
      real rq(nbin),vtr(nbin) !,rdrd(nbin)
      
      double precision, save :: rda(nbin,2),rma(nbin,2),rva(nbin,2),rvna(nbin,2),rdrda(nbin,2),sqrtreynolds(nbin)
      double precision, save :: sdzf(nbin,2) ! 
      double precision :: sqrtvf(nbin),termfp1(nbin),termfm1(nbin)
      
      real :: cracrbin,zracrbin

      double precision :: alpdp,rdiadp,rdiadpinv,tmpdp, xdnfrac, tmpqr, tmpqxrain
      real :: xden,xmlt,cmlt,cmlttot,fventm,fventh,am,ah,felfinv,dmwdt
      integer :: kf
      real(kind=r8) :: savec(45),savex(45),save1,finx(45),xprfr1,xmlta(45),cmlta(45),savexf(45),savez(45),savedz(45),savexh(45)
      real(kind=r8), save :: savesx(45),savesxf(45),sum1(45),sum2(45),sum1a(45),sum2a(45)
      double precision :: szftmp(45,2)
      
      real :: xf(ntakpd,3,ngs) ! ,xprf(ntakpd,2),xprfr(ntakpd,2)
      real :: sxf(ntakpd,3)
      real :: xr(ntakpd,ngs) ! rain
      real :: xc(ntakpd,ngs) ! droplets
      double precision :: totn,totq,totn2,totq2,totnr,totqr,totz,totz2,totzr,totn2a,totq2a,totz2a,totzmlt,totzmlt2
      double precision :: totzr1,totzr2,totzr3,totzr4
      double precision :: totzsmall,totzlarge,totzsmall2,totzlarge2
      double precision :: zmax,xmlttmp,cmlttmp,xftmp
      integer :: izmax

!          real rdrda(ntakpd),rda(ntakpd),rma(ntakpd)=0.0,rvna(ntakpd)

       
       real vtra(nbin)
       real hmmings,hjogs
!       parameter ( hjo = 0.8*7.5*nbin/(41.) )
       parameter (hmmings = 1.e-11, hjogs = 0.8*7.5 )
       
       integer, save :: imake = 0
!       real hmmin,hjo
!       parameter ( hjo = 0.8*7.5*nbin/(41.) )
      real :: qhmlrtmp,qhmlrtmp2, chmlrtmp, chmlrtmpd1inf, chlmlrtmp, zhlmlrtmp, zhlmlrrtmp, qvs0,tmpcmlt

      real :: term1,term2,term3,term4
      real :: qaacw ! combined qsacw-qhacw for WSM6 variation
      real :: cwchtmp

      real, parameter ::  c1r=19.0, c2r=0.6, c3r=1.8, c4r=17.0   ! rain
      real, parameter ::  c1h=5.5, c2h=0.7, c3h=4.5, c4h=8.5   ! Graupel
      real, parameter ::  c1hl=3.7, c2hl=0.3, c3hl=9.0, c4hl=6.5, c5hl=1.0, c6hl=6.5 ! Hail


! inline functions for Newton method
       real :: galpha, dgalpha, zraten, zrateq, zrateqn
       real :: a_in
       logical, parameter :: newton = .false.

! inline function for charging formula
      ftau(tau) = (-1.7e-5)*tau**3 - 0.003*tau**2 - 0.05*tau + 0.13


!      galpha(a_in) = ((4. + a_in)*(5. + a_in)*(6. + a_in))/((1. + a_in)*(2. + a_in)*(3. + a_in))
!      dgalpha(a_in) = (876. + 1260.*a_in + 621.*a_in**2 + 126.*a_in**3 + 9.*a_in**4)/            &
!     &  (36. + 132.*a_in + 193.*a_in**2 + 144.*a_in**3 + 58.*a_in**4 + 12.*a_in**5 + a_in**6)
!
! ####################################################################
!
!  Start routine
!
! ####################################################################



!

      imixedphase = 0
      IF ( mixedphase ) imixedphase = 1
      itile = nxi
      jtile = nyj
!      write(0,*) 'GS: itile,jtile,nx,ny = ',nxi,nyj,nx,ny
!      ktile = nz
!      ixbeg = its
!      jybeg = jyslab
!      ixend = nx
!      jyend = ny
!      kzend = nz
!      nxend = ide
!      nyend = jde
!      nzend = nz
!      kzbeg = 1
!      nzbeg = 1

      istag = 0
      jstag = 0
      kstag = 1


      lrescalelow(:) = rescale_low_alpha
      lrescalelow(lr) = rescale_low_alphar .and. rescale_low_alpha
      lrescalelow(lh) = rescale_low_alphah .and. rescale_low_alpha
      IF ( lf > 1 ) lrescalelow(lf) = rescale_low_alphah .and. rescale_low_alpha
      IF ( lhl > 1 ) lrescalelow(lhl) = rescale_low_alphahl .and. rescale_low_alpha


!
!  slope intercepts
!

      IF ( ngs .lt. nz ) THEN
!       write(0,*) 'Error in ICEZVD: Must have ngs .ge. nz!'
!       STOP
      ENDIF

      cntnic_noliq = 0
      q_noliqmn = 0.0
      q_noliqmx = 0.0
      scsacimn = 0.0
      scsacimx = 0.0

      ldovol = .false.

      DO il = lc,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO

      delqnxa(lc:lqmx) = 0.0 ! delqnia
      delqxxa(lc:lqmx) = 0.0 ! delqxia

      delqnxa(li) = delqnia
      delqxxa(li) = delqxia

      IF ( lis > 1 ) THEN
        delqnxa(lis) = delqnia
        delqxxa(lis) = delqxia
      ENDIF

      delqnxa(ls) = delqnsa
      delqxxa(ls) = delqxsa

      ffrzh = 1
      IF ( lf > 1 ) THEN
        ffrzh = 0
      ENDIF
!      DO il = lc,lhab
!        write(iunit,*) 'delqnxa(',il,') = ',delqnxa(il)
!      ENDDO
      
!
!  density maximums and minimums
!

!
!  Set terminal velocities...
!    also set drag coefficients
!

      dtpinv = 1.d0/dtp

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
!
!  constants
!

!      cp608 = 0.608
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
!      cs = 4.83607122
!      ds = 0.25
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pii = piinv ! 1./pi
      pid4 = pi/4.0 
!      qscrit = 6.0e-04
      gf1 = 1.0 ! gamma(1.0)
      gf1p5 = 0.8862269255  ! gamma(1.5)
      gf2 = 1.0 ! gamma(2.0)
      gf3 = 2.0 ! gamma(3.0)
      gf3p5 = 3.32335097 ! gamma(3.5)
      gf4 = 6.00 ! gamma(4.0)
      gf5 = 24.0 ! gamma(5.0)
      gf6 = 120.0 ! gamma(6.0)
      gf7 = 720.0 ! gamma(7.0)
      gf4br = 17.837861981813607 ! gamma(4.0+br)
      gf4ds = 10.41688578110938 ! gamma(4.0+ds)
      gf4p5 = 11.63172839656745 ! gamma(4.0+0.5)
      gf3ds = 3.0458730354120997 ! gamma(3.0+ds)
      gf1ds = 0.8863557896089221 ! gamma(1.0+ds)

      gf43rds = 0.8929795116 ! gamma(4./3.)
      gf53rds = 0.9027452930 ! gamma(5./3.)
      gf73rds = 1.190639349 ! gamma(7./3.)
      gf83rds = 1.504575488 ! gamma(8./3.)
      
      gamice73fac =  (Gamma_sp(7./3. + cinu))**3/ (Gamma_sp(1. + cinu)**3 * (1. + cinu)**4)
      gamsnow73fac =  (Gamma_sp(7./3. + snu))**3/ (Gamma_sp(1. + snu)**3 * (1. + snu)**4)
      
!      gcnup1 = Gamma_sp(cnu + 1.)
!      gcnup2 = Gamma_sp(cnu + 2.)
!
!  constants
!
!
!  general constants for microphysics
!
      brz = 100.0
      arz = 0.66
      
      bfnu1 = (4. + alphar)*(5. + alphar)*(6. + alphar)/ &
     &       ((1. + alphar)*(2. + alphar)*(3. + alphar))

       galpharaut = (6.+alpharaut)*(5.+alpharaut)*(4.+alpharaut)/ &
     &             ((3.+alpharaut)*(2.+alpharaut)*(1.+alpharaut))
      
      vfrz = 0.523599*(dfrz)**3 
      vmlt = Min(xvmx(lr), 0.523599*(dmlt)**3 )
      vshd = Min(xvmx(lr), 0.523599*(dshd)**3 )

      IF ( snowmeltdia > 0.0 ) THEN
        snowmeltmass = pi/6.0 * 1000. * snowmeltdia**3  ! maximum rain particle mass from melting snow (if snowmeltdia > 0)
      ENDIF

      tdtol = 1.0e-05
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
      
      IF ( mixedphase ) THEN
       ibinhmlr = 0
       ibinhlmlr = 0
      ENDIF
!
!
! #ifdef 1
!      print*,'ventr,ventc = ',ventr,ventc

!
!  Set up look up tables for supersaturation w.r.t. liq and ice
!
!VD$L SKIP
!      do l = 1,nqsat
!      temq = 163.15 + (l-1)*fqsat
!      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
!      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
!      end do

      mltmass0inv = 1.0/( 1000.0* xvmx(lr) ) ! for drops melting from ice with diameter > 1.9cm
      mltmass1inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize1)**3) ) ! for drops melting from ice with diameter > 1.9cm; 0.01 converts cm to m, 0.5 conv. diam to radius
      mltmass2inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize2)**3) ) ! for drops melting from ice with 0.9cm < d < 1.9cm (or 1.6cm to 1.9cm)
      mltmass3inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize3)**3) ) ! for drops melting from ice with 0.9cm < d < 1.6cm
      mltmass1cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize1)**3) 
      mltmass2cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize2)**3) 
      mltmass3cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize3)**3) 
      
!      real, parameter :: mltdiam1 = 9.0e-3, mltdiam2 = 19.0e-3, mltdiam05 = 4.5e-3

      IF ( ibinnum == 1 ) THEN
        numdiam = 1 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = 4.5e-3
      ELSEIF ( ibinnum == 2 ) THEN
        numdiam = 2 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = mltdiam1/6. ! 1.5e-3
        mltdiam(2) = mltdiam1/2. ! 4.5e-3
      ELSEIF ( ibinnum > 2 ) THEN
        numdiam = Min(ibinnum, ndiam)
        DO k = 1,numdiam
          mltdiam(k) = (k - 0.5)*mltdiam1/float(numdiam)
        ENDDO
      
      ELSE
        numdiam = 5 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = 0.5e-3
        mltdiam(2) = 1.0e-3
        mltdiam(3) = 2.0e-3
        mltdiam(4) = 4.0e-3
        mltdiam(5) = 6.0e-3
      ENDIF


      IF ( numshedregimes == 2 ) THEN
        mltdiam(ndiam+1) = mltdiam1 !  9.0e-3
        mltdiam(ndiam+2) = mltdiam3 ! 19.0e-3
        mltdiam(ndiam+3) = mltdiam4 !100.0e-3
      ELSEIF ( numshedregimes == 3 ) THEN
        mltdiam(ndiam+1) = mltdiam1 !  9.0e-3
        mltdiam(ndiam+2) = mltdiam2 ! 16.0e-3
        mltdiam(ndiam+3) = mltdiam3 ! 19.0e-3
        mltdiam(ndiam+4) = mltdiam4 !200.0e-3
      ENDIF


      IF ( imake == 0 ) THEN
!      write(iunit,*) 'mltmass1inv,takshedsize1 = ',mltmass1inv,takshedsize1
!      write(iunit,*) 'mltmass2inv,takshedsize2 = ',mltmass2inv,takshedsize2
        imake = 1
        
!        IF ( ibinhmlr < 2 .and. ibinhlmlr < 2) THEN
!        DO l = 1,nbin
!         rma(l,1)  = hmmings*exp(3.0*(l-1)/hjogs)
!         rva(l,1)  = rma(l,1)/xdn0(lh)  ! volume (mass/1000.)
!         rda(l,1)  = (6.*rma(l,1)/(pi*xdn0(lh)))**(1./3.)
!         rvna(l,1) = 1.
!         rdrda(l,1) = rda(l,1) /hjogs
!         
!         x = rda(l,1)
!                  
!!        write(6,*) 'l,rma,rva,rvna,rda,vtra,rdrda = ',l,rma(l),  &
!!     &     rva(l),rvna(l),rda(l),rdrda(l)
!        ENDDO
!        ENDIF


       IF ( ibinhmlr > 1 .or. ibinhlmlr > 1 .or. ibinhacw > 0 ) THEN
        IF ( my_rank == 0 ) write(0,*) 'must use commastak compile for bin hybrid! ibinhmlr,ibinhacw = ',ibinhmlr,ibinhacw
        STOP
       
       ELSEIF ( ( ibinhmlr == 1 .or. ibinhlmlr == 1) ) THEN
!           call takinitbin
!
!        DO k=1,2
!         kf = 2
!        DO l = 1,ntakpd
!         rma(l,k)  = 1.e-3*sxf(l,kf) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
!!         rva(l,kf)  = rma(l,kf)/(xden)  ! (m^3) volume (mass/1000.)
!         rda(l,k)  = 2.*1.e-2*srf(l,kf) ! (cm to meters and radius to diameter) (6.*rma(l,kf)/(pi*xdn))**(1./3.)
!!         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
!         rvna(l,k) = 1.
!         rdrda(l,k) = rda(l,k)/dj
!         sqrtreynolds(l) = (2.0*0.01*sr(l))**(0.5)
!         IF ( my_rank == 0 ) write(0,*) 'l,rma,rda,sr2,rdrda = ',l,rma(l,kf),rda(l,kf),2.*1.e-2*sr(l),rdrda(l,kf)
!        ENDDO
!        ENDDO
!          do l=1,lmax
!            sum1  (l) = workwp(l)
!            sum2  (l) = workwm(l)
!            savesx(l) = sx    (l)
!          enddo
       
       ENDIF

       ENDIF ! imake
!
      
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag


!
!  cw constants in mks units
!
!      cwmasn = 4.25e-15  ! radius of 1.0e-6
      mwfac = 6.0**(1./3.)
      IF ( ipconc .ge. 2 ) THEN
!        cwmasn = xvmn(lc)*1000.
!        cwradn = 1.0e-6
!        cwmasx = xvmx(lc)*1000.
      ENDIF
        rwmasn = xvmn(lr)*1000.
        rwmasx = xvmx(lr)*1000.

      IF ( biggsnowdiam > 0.0 ) THEN
        xvbiggsnow = (pi/6.0)*biggsnowdiam**3
      ELSE
        xvbiggsnow = xvmn(lh)
      ENDIF

!
!  ci constants in mks units
!
      cimasn = Min(cimas0, cimas1) ! 12 microns for  0.1871*(xmas(mgs,li)**(0.3429))
      cimasx = 1.0e-8   ! 338 microns
      ccimx = 5000.0e3   ! max of 5000 per liter

!
!  constants for paramerization
!
!
!  set save counter (number of saves):  nsvcnt
!
!      nsvcnt = 0
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

      iwetg = 0
      iwetg1 = 0
      ngscmaxy(lg:lhab,lc:ls) = 0
      ngscminy(lg:lhab,lc:ls) = 0
      ngscxymax(lg:lhab,lc:ls) = 0

!      timetd1 = etime(tarray)
!      timetd1 = tarray(1)

!
!***********************************************************
!  start jy loop
!***********************************************************
!

!      do 9999 jy = 1,ny-jstag
!
!  VERY IMPORTANT:  SET jy = jgs
!
      jy = jgs
     
     
!      t1(:,:,:) = 0
!      t2(:,:,:) = 0
!      t3(:,:,:) = 0
!      t4(:,:,:) = 0
!      t5(:,:,:) = 0
!      t6(:,:,:) = 0
!      t8(:,:,:) = 0
      
      IF ( ipconc < 2 ) THEN ! Make a copy of cloud droplet mixing ratio to use for homogeneous freezing
        DO kz = 1,kze
         DO ix = 1,itile
           t9(ix,jy,kz) = an(ix,jy,kz,lc)
         ENDDO
        ENDDO
      ENDIF
      
!
!..Gather microphysics  
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: ENTER GATHER STAGE'


      
      nxmpb = 1
      nzmpb = 1
      nxz = itile*nz
      numgs = nxz/ngs + 1
!      write(0,*) 'ICEZVD_GS: ENTER GATHER STAGE: nx,nz,nxz,numgs,ngs = ',nx,nz,nxz,numgs,ngs

      do 1000 inumgs = 1,numgs
      ngscnt = 0
      
      do kz = nzmpb,kze
      do ix = nxmpb,itile

      pqs(1) = t00(ix,jy,kz)
      pres(1) = pn(ix,jy,kz) + pb(kz)

      theta(1) = an(ix,jy,kz,lt)
      temg(1) = t0(ix,jy,kz)
      temcg(1) = temg(1) - tfr
      tqvcon = temg(1)-cbw
      ltemq = (temg(1)-163.15)/fqsat + 1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      IF ( iqvsopt == 0 ) THEN
        qvs(1) = pqs(1)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(1) = rdorv*esbolton*tabqvs(ltemq)/(pres(1) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(1),POLYSVP1(temg(1),0))
        qvs(1) = rdorv*tmp/(pres(1) - tmp)
      ENDIF

      IF ( iqis0 == 1 .or. temg(1) <= tfr+0.5 ) THEN
        qis(1) = pqs(1)*tabqis(ltemq)
      ELSE
        ltemq = (tfr - 163.15)/fqsat + 1.5
        qis(1) = pqs(1)*tabqis(ltemq)
      ENDIF

      qss(1) = qvs(1)

      if ( temg(1) .lt. tfr ) then
        qss(1) = qis(1)
      end if
!
      ishail = .false.
      IF ( lhl > 1 ) THEN
        IF ( an(ix,jy,kz,lhl)  .gt. qxmin(lhl) ) ishail = .true.
      ENDIF
      IF ( lf > 1 ) THEN
        ishail = ishail .or. an(ix,jy,kz,lf)  .gt. qxmin(lf)
      ENDIF


      
      if ( an(ix,jy,kz,lv)  .gt. qss(1) .or.   &
     &     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.    &
     &     an(ix,jy,kz,li)  .gt. qxmin(li)   .or.   &
     &     an(ix,jy,kz,lr)  .gt. qxmin(lr)   .or.   &
     &     an(ix,jy,kz,ls)  .gt. qxmin(ls)   .or.   &
     &     an(ix,jy,kz,lh)  .gt. qxmin(lh)   .or.  ishail ) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 1100
      end if
      enddo !ix
      nxmpb = 1
      enddo !kz
 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998

      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: dbg = 5, ngscnt = ',ngscnt
      
!      write(0,*) 'allocating qc'

      
      xv(:,:) = 0.0
      xsfca(:,:) = 0.0
      xmas(:,:) = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:) = 0.0
      raindn(:,:) = 900.
      cx(:,:) = 0.0
      IF ( lnhf > 1 .or. lnhlf > 1 ) chxf(:,:) = 0.0
      alpha(:,:) = 0.0
      DO il = li,lhab
        DO mgs = 1,ngscnt
          rimdn(mgs,il)  = rimedens ! xdn0(il)
        ENDDO
      ENDDO
!
!  define temporaries for state variables to be used in calculations
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: dbg = def temps'
      do mgs = 1,ngscnt
      kgsm(mgs) = max(kgs(mgs)-1,1)
      kgsp(mgs) = min(kgs(mgs)+1,nz-1)
      kgsm2(mgs) = Max(kgs(mgs)-2,1)
      theta0(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt) - theta0(mgs)
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      qv0(mgs) = an(igs(mgs),jy,kgs(mgs),lv)
      qwvp(mgs) = an(igs(mgs),jy,kgs(mgs),lv)  - qv0(mgs) ! qv0(mgs) is full qv, so qwvp starts as zero!

      pres(mgs) = pn(igs(mgs),jy,kgs(mgs)) + pb(kgs(mgs))
      pipert(mgs) = p2(igs(mgs),jy,kgs(mgs))
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhoinv(mgs) = 1.0/rho0(mgs)
      rhovt(mgs) = Sqrt(rho00/Max(0.05,rho0(mgs))) ! prevent excessive rhovt
      pi0(mgs) = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      temgkm1(mgs) = t0(igs(mgs),jy,kgsm(mgs))
      temgkm2(mgs) = t0(igs(mgs),jy,kgsm2(mgs))
      pk(mgs)   = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs)) ! t77(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
      qss0(mgs) = (380.0)/(pres(mgs))
      pqs(mgs) = (380.0)/(pres(mgs))
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      IF ( iqvsopt == 0 ) THEN
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
        qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF

      IF ( iqis0 == 1 .or. temg(mgs) <= tfr+0.5 ) THEN
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
      ELSE
        ltemq = (tfr - 163.15)/fqsat + 1.5
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
      ENDIF
      qss(mgs) = qvs(mgs)
!      es(mgs)  = 6.1078e2*tabqvs(ltemq)
!      eis(mgs) = 6.1078e2*tabqis(ltemq)
      cnostmp(mgs) = cno(ls)
!

      il5(mgs) = 0
      if ( temg(mgs) .lt. tfr ) then
      il5(mgs) = 1
      end if
      enddo !mgs
      
      IF ( ipconc < 1 .and. lwsm6 ) THEN
        DO mgs = 1,ngscnt
          tmp = Min( 0.0, temcg(mgs) )
          cnostmp(mgs) = Min( 2.e8, 2.e6*exp(0.12*tmp) )
        ENDDO
      ENDIF


!
! zero arrays that are used but not otherwise set (tm)
!
      do mgs = 1,ngscnt
         qhshr(mgs) = 0.0 
       end do
!
!  set temporaries for microphysics variables
!
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

      qxw(:,:) = 0.0
      qxwlg(:,:) = 0.0



       IF ( iraintypes > 0 ) THEN
       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = load raintypes'
         
!         write(0,*) 'il,lrain = ',il,lrain(il)
         DO mgs = 1,ngscnt
           tmpdp = 0
           DO il = 1,nraintypes
            qxrain(mgs,il) = Max(0.0, Min(qx(mgs,lr), an(igs(mgs),jy,kgs(mgs),lrain(il)) ) )
            qxrainold(mgs,il) = Max(0.0, Min(qx(mgs,lr), an(igs(mgs),jy,kgs(mgs),lrain(il)) ) )
            tmpdp = tmpdp + qxrainold(mgs,il)
           ENDDO
           ! make sure that parts add up to whole:
           IF ( tmpdp > 0.0d0 ) THEN
           DO il = 1,nraintypes
             qxrainold(mgs,il) = qxrainold(mgs,il)*qx(mgs,lr)/tmpdp
           ENDDO
           ENDIF
           
         ENDDO
          DO mgs = 1,ngscnt
           qxrainold(mgs,0) = Max(0.0, qx(mgs,lr) )
          ENDDO
          qxshedfrac(:,:,:) = 0.0
          qxshedfrac(:,lr,2) = 1.0 ! default count all as shed rain
          
       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = load raintypes2'
       ENDIF



!
!  set concentrations
!
!      ssmax = 0.0
      
      
      if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*)  'ICEZVD_GS: dbg = 5b'
      
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
          IF ( qx(mgs,li) .le. qxmin(li) ) THEN
            cx(mgs,li) = 0.0
          ENDIF

        IF ( lcina .gt. 1 ) THEN
         cina(mgs) = an(igs(mgs),jy,kgs(mgs),lcina)
        ELSE
         cina(mgs) = cx(mgs,li)
        ENDIF
        IF ( lcin > 1 ) THEN
         ccin(mgs) = an(igs(mgs),jy,kgs(mgs),lcin)
        ENDIF
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
!        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
        IF ( qx(mgs,lc) .le. qxmin(lc) ) THEN
          cx(mgs,lc) = 0.0
        ENDIF
        IF ( lss > 1 ) THEN
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        ENDIF
        IF ( lccn .gt. 1 ) THEN
         ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        ELSE
         ccnc(mgs) = 0.0
        ENDIF
        IF ( lccna .gt. 1 ) THEN
         ccna(mgs) = an(igs(mgs),jy,kgs(mgs),lccna)
        ELSE
         ccna(mgs) = cx(mgs,lc)
        ENDIF
       end do
!       ELSE
!       cx(mgs,lc) = Abs(ccn)
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
!          cx(mgs,lr) = 0.0
        ELSEIF ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) THEN
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lr)
          qx(mgs,lr) = 0.0
        ELSE
          cx(mgs,lr) = Max( 1.e-9, cx(mgs,lr) )
        ENDIF
       end do
      end if
      if ( ipconc .ge. 4 ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
!          cx(mgs,ls) = 0.0
        ELSEIF ( cx(mgs,ls) .eq. 0.0 .and. qx(mgs,ls) .lt. 3.0*qxmin(ls) ) THEN
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,ls)
          qx(mgs,ls) = 0.0
        ELSE
          cx(mgs,ls) = Max( 1.e-9, cx(mgs,ls) )

         IF ( ilimit .ge. ipc(ls) ) THEN
            tmp = (xdn0(ls)*cx(mgs,ls))/(rho0(mgs)*qx(mgs,ls))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,ls)*(tmp2)
         IF ( cnox .gt. 3.0*cno(ls) ) THEN
           cx(mgs,ls) = 3.0*cno(ls)/tmp2
         ENDIF
         ENDIF
        ENDIF
       end do
      end if
      if ( ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
!          cx(mgs,lh) = 0.0
        ELSEIF ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lh) 
          qx(mgs,lh) = 0.0
        ELSE
          cx(mgs,lh) = Max( 1.e-9, cx(mgs,lh) )
         IF ( ilimit .ge. ipc(lh) ) THEN
            tmp = (xdn0(lh)*cx(mgs,lh))/(rho0(mgs)*qx(mgs,lh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lh)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lh) ) THEN
           cx(mgs,lh) = 3.0*cno(lh)/tmp2
         ENDIF
         ENDIF
        ENDIF

        IF ( lf > 1 ) THEN
          cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)
        IF ( qx(mgs,lf) .le. qxmin(lf) ) THEN
         !   cx(mgs,lf) = 0.0
        ELSEIF ( cx(mgs,lf) .eq. 0.0 .and. qx(mgs,lf) .lt. 3.0*qxmin(lf) ) THEN
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lf) 
          qx(mgs,lf) = 0.0
        ELSE
          cx(mgs,lf) = Max( 1.e-9, cx(mgs,lf) )
          IF ( ilimit .ge. ipc(lf) ) THEN
             tmp = (xdn0(lf)*cx(mgs,lf))/(rho0(mgs)*qx(mgs,lf))
             tmp2 = (tmp*(3.14159))**(1./3.)
             cnox = cx(mgs,lf)*(tmp2)
            IF ( cnox .gt. 3.0*cno(lf) ) THEN
              cx(mgs,lf) = 3.0*cno(lf)/tmp2
            ENDIF
          ENDIF
        ENDIF
        
        ENDIF

        IF ( lnhf > 1 ) THEN
           chxf(mgs,lh) = Min(cx(mgs,lh), Max(an(igs(mgs),jy,kgs(mgs),lnhf), 0.0))
        ENDIF


       end do


      end if

      if ( lhl .gt. 1 .and. ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        IF ( qx(mgs,lhl) .le. qxmin(lhl) ) THEN
          cx(mgs,lhl) = 0.0
        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) THEN
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lhl) 
          qx(mgs,lhl) = 0.0
        ELSE
          cx(mgs,lhl) = Max( 1.e-9, cx(mgs,lhl) )
         IF ( ilimit .ge. ipc(lhl) ) THEN
            tmp = (xdn0(lhl)*cx(mgs,lhl))/(rho0(mgs)*qx(mgs,lhl))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lhl)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lhl) ) THEN
           cx(mgs,lhl) = 3.0*cno(lhl)/tmp2
         ENDIF
         ENDIF
        ENDIF

        IF ( lnhlf > 1 ) THEN ! number of hail from frozen drops
           chxf(mgs,lhl) = Min(cx(mgs,lhl), Max(an(igs(mgs),jy,kgs(mgs),lnhlf), 0.0))
!           chxf(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhlf), 0.0)
        ENDIF

       end do
      end if

!
! Set mean particle volume
!
      IF ( ldovol ) THEN
      
      vx(:,:) = 0.0
      
       DO il = li,lhab
        
        IF ( lvol(il) .ge. 1 ) THEN
        
          DO mgs = 1,ngscnt
            vx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lvol(il)), 0.0)
          ENDDO

        ENDIF

       ENDDO

      ENDIF


!
! Set liquid water fraction
!
      fhw(:) = 0.0
      fsw(:) = 0.0
      fhlw(:) = 0.0
      ffw(:) = 0.0



!
!  6th moments
!

      IF ( ipconc .ge. 6 ) THEN
       zx(:,:) = 0.0
       DO il = lr,lhab
        IF ( lz(il) .gt. 1 ) THEN
         DO mgs = 1,ngscnt
          zx(mgs,il) = Max( an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0 )
         ENDDO
        ENDIF
       ENDDO

      ENDIF

      IF ( ipconc .ge. 6 ) THEN

         tmp = alphamax - 1.0
         g1xmax = (6.0 + tmp)*(5.0 + tmp)*(4.0 + tmp)/ &
     &            ((3.0 + tmp)*(2.0 + tmp)*(1.0 + tmp))
         g1xmin = (6.0 + alphamin)*(5.0 + alphamin)*(4.0 + alphamin)/ &
     &            ((3.0 + alphamin)*(2.0 + alphamin)*(1.0 + alphamin))

       IF ( lz(lr) .lt. 1 ) THEN
         g1x(:,lr) = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/ &
     &            ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))

         
         DO mgs = 1,ngscnt
           IF ( cx(mgs,lr) .gt. 0.0 .and. qx(mgs,lr) .gt. qxmin(lr)  ) THEN
            
            vr = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            IF ( lzr < 1 ) THEN
             IF ( imurain == 3 ) THEN
               zx(mgs,lr) = 3.6476*(rnu+2.0)*cx(mgs,lr)*vr**2/(rnu+1.0)
             ELSE ! imurain == 1
               zx(mgs,lr) = 3.6476*g1x(mgs,lr)*cx(mgs,lr)*vr**2
             ENDIF
            ENDIF
             
           ENDIF
         ENDDO
       ENDIF
      
      ENDIF


         IF ( ipconc == 5 ) THEN
         ! set up factors for ihlcnh=3 conversion
         g1x(:,lr) = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/ &
     &              ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
         g1x(:,lh) = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/ &
     &               ((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))
           IF ( lhl > 0 ) THEN
            g1x(:,lhl) = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/ &
     &               ((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))
           ENDIF
         ENDIF

! load charges
      IF ( ipelec > 0 .and. lscw .gt. 1 ) THEN
        ierr = 0
        DO il = lc,lhab
          DO mgs = 1,ngscnt
            scx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lsc(il))
!            IF ( .not. (scx(mgs,il) > -1.e-6 .and. scx(mgs,il) < 1.e-6 ) ) THEN
!            IF ( Abs (scx(mgs,il) ) > 1000.e-9 .or. Abs (scx(mgs,il) ) > 1.e-9 ) THEN ! DEBUGTED
            IF ( Abs (scx(mgs,il) ) > 1000.e-9 ) THEN 
             ierr = ierr + 1
             write(0,*) 'Problem0a with scx il = ',il,lsc(il),scx(mgs,il),qx(mgs,il),cx(mgs,il)
             write(0,*) 'at ix,jy,kz = ',igs(mgs),jyslab,kgs(mgs)
          IF ( ierr > 10 ) STOP
             IF ( il == ls ) THEN
             IF ( Abs (scx(mgs,il) ) > 10.e-9 ) THEN
               write(0,*) 'myrank,scs,mgs = ',my_rank,scx(mgs,il),mgs
             ENDIF
             ENDIF
            ENDIF
          ENDDO
        ENDDO
          IF ( ierr > 0 ) STOP

      DO mgs = 1,ngscnt
       cionp(mgs) = an(igs(mgs),jy,kgs(mgs),lscpi)
       cionn(mgs) = an(igs(mgs),jy,kgs(mgs),lscni)
      if ( largeion ) then
       clionp(mgs) = an(igs(mgs),jy,kgs(mgs),lscpli)
       clionn(mgs) = an(igs(mgs),jy,kgs(mgs),lscnli)
      end if
      ENDDO ! mgs 

!$PAR CRITICAL SECTION
      DO mgs = 1,ngscnt
       dvmgs(mgs) = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
      DO il = lc,lhab
       psctot1 = psctot1 + scx(mgs,il)*dvmgs(mgs)
      ENDDO
      ENDDO
!$PAR END CRITICAL SECTION

      ELSE
        scx(:,:) = 0.0
      ENDIF

!
!  set shape parameters
!
       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set alpha'
      IF ( imurain == 1 ) THEN
        alpha(:,lr) = alphar
      ELSEIF ( imurain == 3 ) THEN
        alpha(:,lr) = xnu(lr)
      ENDIF
      
      alpha(:,li) = xnu(li)
      alpha(:,lc) = xnu(lc)
        IF ( idiagnosecnu > 0 ) THEN
          DO mgs = 1,ngscnt
            IF ( cx(mgs,lc) > cxmin ) THEN
              IF ( idiagnosecnu == 1 ) THEN ! diagnose cloud drop DSD shape parameter (cnu) based on Chandrakar et al. 2016 (PNAS) data
                alpha(mgs,lc) = cnudiag( cx(mgs,lc) )
              ELSEIF ( idiagnosecnu > 1 ) THEN ! Geoffrey et al. (2010, ACP, eq. 12 for alpha = nu - 1)
                alpha(mgs,lc) = 1.58*rho0(mgs)*qx(mgs,lc)*1.e6 + 0.72 - 1.0
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      IF ( imusnow == 1 ) THEN
        alpha(:,ls) = alphas
      ELSEIF ( imusnow == 3 ) THEN
        alpha(:,ls) = xnu(ls)
      ENDIF

       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set dab'
      
      DO il = lr,lhab
      do mgs = 1,ngscnt
        IF ( il .ge. lg ) alpha(mgs,il) = dnu(il)


        DO ic = lc,lhab
        dab0lh(mgs,il,ic) =  dab0(il,ic) ! dab0(ic,il)
        dab1lh(mgs,il,ic) =  dab1(il,ic) ! dab1(ic,il)
        ENDDO
      end do
      ENDDO

      
!      DO mgs = 1,ngscnt
        DO il = lr,lhab
          da0lx(:,il) = da0(il)
        ENDDO
        da0lh(:) = da0(lh)
        da0lr(:) = da0(lr)
        da1lr(:) = da1(lr)
        da0lc(:) = da0(lc)
        da1lc(:) = da1(lc)
        IF ( idiagnosecnu > 0 ) THEN
          DO mgs = 1,ngscnt
            IF ( cx(mgs,lc) > cxmin ) THEN
              da0lc(mgs) =  delbk(bb(lc), alpha(mgs,lc), xmu(lc), 0)
              da1lc(mgs) =  delbk(bb(lc), alpha(mgs,lc), xmu(lc), 1)
            ENDIF
          ENDDO
        ENDIF


       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set rz'

        IF ( lf > 1 ) da0lf(:) = da0(lf)
        IF ( lzh < 1 .or. lzhl < 1 ) THEN
          rzxhlh(:) = rzhl/rz
        ELSEIF ( lzh > 1 .and. lzhl > 1 ) THEN
          rzxhlh(:) = 1.
        ENDIF
        IF ( lf > 1 ) THEN
          IF ( lzf < 1 .or. lzhl < 1 ) THEN
            rzxhlf(1:ngscnt) = rzhl/rz
          ELSEIF ( lzf > 1 .and. lzhl > 1 ) THEN
            rzxhlf(1:ngscnt) = 1.
          ENDIF
        ENDIF
        IF ( lzr > 1 ) THEN
          rzxh(:) = 1.
          rzxhl(:) = 1.
          rzxf(:) = 1.
        ELSE
          rzxh(:) = rz
          rzxhl(:) = rzhl
          rzxf(:) = rz
        ENDIF
        
        IF ( imurain == 1 .and. imusnow == 3 .and. lzr < 1 ) THEN
          rzxs(:) = rzs
        ELSEIF ( imurain == imusnow .or. lzr > 1 ) THEN
          rzxs(:) = 1.
        ENDIF
 !     ENDDO
      
      IF ( lhl .gt. 1 ) THEN
      DO mgs = 1,ngscnt
        da0lhl(mgs) = da0(lhl)
      ENDDO
      ENDIF
      
      ventrx(:) = ventr
      ventrxn(:) = ventrn
      gf1palp(:) = gamma_sp(1.0 + alphar)

!
!  set factors
!
      do mgs = 1,ngscnt
!
      ssi(mgs) = qx(mgs,lv)/qis(mgs)
      ssw(mgs) = qx(mgs,lv)/qvs(mgs)
!
      tsqr(mgs) = temg(mgs)**2
!
      temgx(mgs) = min(temg(mgs),313.15)
      temgx(mgs) = max(temgx(mgs),233.15)
      felv(mgs) = 2500837.367 * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
!
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),223.15)
      temcgx(mgs) = temcgx(mgs)-273.15

! felf = latent heat of fusion, fels = LH of sublimation, felv = LH of vaporization
      felf(mgs) = 333690.6098 + (2030.61425)*temcgx(mgs) - (10.46708312)*temcgx(mgs)**2
!
      fels(mgs) = felv(mgs) + felf(mgs)
!
      felvs(mgs) = felv(mgs)*felv(mgs)
      felss(mgs) = fels(mgs)*fels(mgs)
      
        IF ( eqtset <= 1 ) THEN
          felvcp(mgs) = felv(mgs)*cpi
          felscp(mgs) = fels(mgs)*cpi
          felfcp(mgs) = felf(mgs)*cpi
        ELSE
          
          ! equations from appendix in Bryan and Morrison (2012, MWR)
          ! note that rw is Rv in the paper, and rd is R.
          
          tmp = qx(mgs,li)+qx(mgs,ls)+qx(mgs,lh)
          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          IF ( lf > 1 ) tmp = tmp + qx(mgs,lf)
          cvm = cv+cvv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)

          IF ( eqtset == 2 ) THEN ! compact form from treating dT/dt = theta*d(pi)/dt + pi*d(theta)dt and then applied to theta assuming constant pi
          felvcp(mgs) = (felv(mgs)-rw*temg(mgs))/cvm
          felscp(mgs) = (fels(mgs)-rw*temg(mgs))/cvm
          felfcp(mgs) = felf(mgs)/cvm
          
          ELSE
           ! equivalent version that applies separate updates of latent heating to theta and pi, when both are returned.

          cpm = cp+cpv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          rmm=rd+rw*qx(mgs,lv)
          
          felvcp(mgs) = (felv(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
          felscp(mgs) = (fels(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
          felfcp(mgs) = felf(mgs)*cv/(cp*cvm)

          felvpi(mgs) = pi0(mgs)*rovcp*(felv(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm
          felspi(mgs) = pi0(mgs)*rovcp*(fels(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm 
          felfpi(mgs) = pi0(mgs)*rovcp*(felf(mgs)/(cvm*temg(mgs)))
          
          ENDIF

        ENDIF
!
      fgamw(mgs) = felvcp(mgs)/pi0(mgs)
      fgams(mgs) = felscp(mgs)/pi0(mgs)
!
      fcqv1(mgs) = 4098.0258*pi0(mgs)*fgamw(mgs)
      fcqv2(mgs) = 5807.6953*pi0(mgs)*fgams(mgs)
      fcc3(mgs) = felfcp(mgs)/pi0(mgs)
!
!  fwvdf = water vapor diffusivity
      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)*(101325.0/(pres(mgs)))
!
! fadvisc = 'd' for dynamic viscosity
! fakvisc = 'k' for kinematic viscosity
      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*(temg(mgs)/296.0)**(1.5) ! dynamic visc.
!
      fakvisc(mgs) = fadvisc(mgs)*rhoinv(mgs) ! divide by rho_air to get kinematic visc. (note the 'k' vs. 'd')
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
      fcw(mgs) = 4203.1548  + (1.30572e-2)*((temcgx(mgs)-35.)**2)   &
     &                 + (1.60056e-5)*((temcgx(mgs)-35.)**4)
      end if
      if ( temg(mgs) .ge. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),308.15)
      temcgx(mgs) = max(temcgx(mgs),273.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4243.1688  + (3.47104e-1)*(temcgx(mgs)**2)
      end if
!
      ftka(mgs) = tka0*fadvisc(mgs)/advisc1  ! thermal conductivity: proportional to dynamic viscosity
      fthdf(mgs) = ftka(mgs)*cpi*rhoinv(mgs)
!
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))  ! Schmidt number
      fpndl(mgs) = (fakvisc(mgs)/fthdf(mgs))  ! Prandl number (only used for bin melting)
!
      fai(mgs) = (fels(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbi(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qis(mgs)))
      fav(mgs) = (felv(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbv(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qvs(mgs)))

      kp1 = Min(nz, kgs(mgs)+1 )
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1)   &
     &                  +w(igs(mgs),jgs,kgs(mgs)))

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
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set density'
!

      do mgs = 1,ngscnt
        xdn(mgs,li) = xdn0(li)
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
        xdn(mgs,ls) = xdn0(ls)
        xdn(mgs,lh) = xdn0(lh)
        IF ( lvol(ls) .gt. 1 ) THEN
         IF ( vx(mgs,ls) .gt. 0.0 .and. qx(mgs,ls) .gt. qxmin(ls) ) THEN
           xdn(mgs,ls) = Min( xdnmx(ls), Max( xdnmn(ls), rho0(mgs)*qx(mgs,ls)/vx(mgs,ls) ) )
         ENDIF
        ENDIF

        IF ( lvol(lh) .gt. 1 ) THEN
         IF ( vx(mgs,lh) .gt. 0.0 .and. qx(mgs,lh) .gt. qxmin(lh) ) THEN
           IF ( mixedphase ) THEN
           ELSE
             dnmx = xdnmx(lh)
           ENDIF
           xdn(mgs,lh) = Min( dnmx, Max( xdnmn(lh), rho0(mgs)*qx(mgs,lh)/vx(mgs,lh) ) )
           vx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/xdn(mgs,lh)
         
         ELSEIF ( vx(mgs,lh) == 0.0 .and. qx(mgs,lh) .gt. qxmin(lh) ) THEN ! if volume is zero, need to initialize the default value

           vx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/xdn(mgs,lh)
         
         ENDIF
        ENDIF

        IF ( lf > 1 ) THEN
        xdn(mgs,lf) = xdn0(lf)
        IF ( lvol(lf) .gt. 1 ) THEN
         IF ( vx(mgs,lf) .gt. 0.0 .and. qx(mgs,lf) .gt. qxmin(lf) ) THEN
           IF ( mixedphase ) THEN
           ELSE
             dnmx = xdnmx(lf)
           ENDIF
           xdn(mgs,lf) = Min( dnmx, Max( xdnmn(lf), rho0(mgs)*qx(mgs,lf)/vx(mgs,lf) ) )
           vx(mgs,lf) = rho0(mgs)*qx(mgs,lf)/xdn(mgs,lf)
         
         ELSEIF ( vx(mgs,lf) == 0.0 .and. qx(mgs,lf) .gt. qxmin(lf) ) THEN ! if volume is zero, need to initialize the default value

           vx(mgs,lf) = rho0(mgs)*qx(mgs,lf)/xdn(mgs,lf)
         
         ENDIF
        ENDIF
        ENDIF



        IF ( lhl .gt. 1 ) THEN

          xdn(mgs,lhl) = xdn0(lhl)
          xdntmp(mgs,lhl) = xdn0(lhl)

          IF ( lvol(lhl) .gt. 1 ) THEN
           IF ( vx(mgs,lhl) .gt. 0.0 .and. qx(mgs,lhl) .gt. qxmin(lhl) ) THEN

           IF ( mixedphase .and. lhlw > 1 ) THEN
           ELSE
             dnmx = xdnmx(lhl)
           ENDIF

             xdn(mgs,lhl) = Min( dnmx, Max( xdnmn(lhl), rho0(mgs)*qx(mgs,lhl)/vx(mgs,lhl) ) )
             vx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/xdn(mgs,lhl)
             xdntmp(mgs,lhl) = xdn(mgs,lhl)
         
           ELSEIF ( vx(mgs,lhl) == 0.0 .and. qx(mgs,lhl) .gt. qxmin(lhl) ) THEN ! if volume is zero, need to initialize the default value

             vx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/xdn(mgs,lhl)
         
           ENDIF
          ENDIF

        ENDIF


      end do

      IF ( ipconc == 5 .and. imydiagalpha == 2 ) THEN

        cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
        
        DO mgs = 1,ngscnt
          !IF ( igs(mgs) == 19 ) write(0,*) 'k,qr,qh,cr,ch = ',kgs(mgs),qx(mgs,lr),cx(mgs,lr),qx(mgs,lh),cx(mgs,lh)
          IF ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > cxmin ) THEN
             xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))            ! 
             xdia(mgs,lr,3) = (xv(mgs,lr)*6.0*cwc1)**(1./3.) 
           !  alpha(mgs,lr) = Min(alphamax, c1r*tanh(c2r*(xdia(mgs,lr,3)*1000. - c3r)) + c4r)
           ! IF ( igs(mgs) == 19 ) write(0,*) 'imy: i,k,alpr,xdia = ',igs(mgs),kgs(mgs),alpha(mgs,lr),xdia(mgs,lr,3)*1000.

            ! Milbrandt & M-C 2010:
             tmp = 4. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*xdia(mgs,lr,3)*cwchtmp

             alpha(mgs,lr) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
          ENDIF
          IF ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) > cxmin ) THEN
!      MY 2005:
             xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))            ! 
             xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
!             alpha(mgs,lh) = Min(alphamax, c1h*tanh(c2h*(xdia(mgs,lh,3)*1000. - c3h)) + c4h)

            ! Milbrandt & M-C 2010:
             tmp = 4. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*xdia(mgs,lh,3)*cwchtmp

             alpha(mgs,lh) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
            ! alphan(mgs,lh) = alpha(mgs,lh)
            
           ! IF ( igs(mgs) == 19 ) write(0,*) 'imy: i,k,alph,xdia = ',igs(mgs),kgs(mgs),alpha(mgs,lh),xdia(mgs,lh,3)*1000.
            il = lh
            DO ic = lc,lh-1 ! lhab
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               IF ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) THEN
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               ELSE ! IF ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               ENDIF
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)
             ENDDO
          ENDIF
!        alpha(:,lr) = 0. ! 10.
!        alpha(:,lh) = 0. ! 10.
          IF ( lhl > 0 ) THEN
          IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) > cxmin ) THEN
             xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))            ! 
             xdia(mgs,lhl,3) = (xv(mgs,lhl)*6.*piinv)**(1./3.)
             IF ( xdia(mgs,lhl,3) < 0.008 ) THEN
               alpha(mgs,lhl) = Min(alphamax, c1hl*tanh(c2hl*(xdia(mgs,lhl,3)*1000. - c3hl)) + c4hl)
             ELSE
               alpha(mgs,lhl) = Min(alphamax, c5hl*xdia(mgs,lhl,3)*1000. + c6hl)
             ENDIF

            il = lhl
            DO ic = lc,lh-1 ! lhab
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               IF ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) THEN
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               ELSE ! IF ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               ENDIF
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)
             ENDDO

          ENDIF
          ENDIF



        ENDDO
      ENDIF
      

       IF ( imurain == 3 ) THEN
         IF ( lzr > 1 ) THEN
           alphashr = 0.0
           alphamlr = -2.0/3.0
           alphasmlr = -2.0/3.0
         ELSE
           alphashr = xnu(lr)
           alphamlr = xnu(lr)
           alphasmlr = xnu(lr)
         ENDIF
!         massfacshr = ( (2. + 3.*(1. +alphashr) )/( 3.*(1. + alphashr) ) )**(1./3.) ! this is the diameter factor
!         massfacmlr = ( (2. + 3.*(1. +alphamlr) )/( 3.*(1. + alphamlr) ) )**(1./3.)
         massfacshr = ( (2. + 3.*(1. +alphashr) )**3/( 3.*(1. + alphashr) ) )  ! this is the mass or volume factor
         massfacmlr = ( (2. + 3.*(1. +alphamlr) )**3/( 3.*(1. + alphamlr) ) )
       ELSEIF ( imurain == 1 ) THEN
         IF ( lzr > 1 ) THEN
           alphashr = 4.0
           alphamlr = 4.0
           alphasmlr = alphasmlr0
         ELSE
           alphashr = alphar
           alphamlr = alphar
           alphasmlr = alphar
         ENDIF
!         massfacshr = (3.0 + alphashr)*((3.+alphashr)*(2.+alphashr)*(1. + alphashr) )**(-1./3.) ! this is the diameter factor
!         massfacmlr = (3.0 + alphamlr)*((3.+alphamlr)*(2.+alphamlr)*(1. + alphamlr) )**(-1./3.)
         massfacshr = (3.0 + alphashr)**3/((3.+alphashr)*(2.+alphashr)*(1. + alphashr) ) ! this is the mass or volume factor
         massfacmlr = (3.0 + alphamlr)**3/((3.+alphamlr)*(2.+alphamlr)*(1. + alphamlr) )
       ENDIF
       
!  Find shape parameter rain

      g1shr = 1.0
      g1mlr = 1.0
      g1smlr = 1.0
 
!      CALL cld_cpu('Z-MOMENT-1')  
      
      IF ( ipconc >= 6 ) THEN
      
      ! set base g1x in case rain is not 3-moment
       IF ( ipconc >= 6 .and. imurain == 3 ) THEN
         il = lr
         DO mgs = 1,ngscnt
!           g1x(mgs,il) = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           g1x(mgs,il) = (alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0))
         ENDDO
       ENDIF

      IF (lzr > 1 ) THEN
       IF ( imurain == 3 ) THEN
         g1shr = (alphashr+2.0)/((alphashr+1.0))
         g1mlr = (alphamlr+2.0)/((alphamlr+1.0))
         g1smlr = (alphasmlr+2.0)/((alphasmlr+1.0))
       ELSEIF ( imurain == 1 ) THEN
!         g1shr = 36.*(6.0 + alphashr)*(5.0 + alphashr)*(4.0 + alphashr)/ &
!     &            (pi**2*(3.0 + alphashr)*(2.0 + alphashr)*(1.0 + alphashr))
         g1shr = (6.0 + alphashr)*(5.0 + alphashr)*(4.0 + alphashr)/ &
     &            ((3.0 + alphashr)*(2.0 + alphashr)*(1.0 + alphashr))
!         g1mlr = 36.*(6.0 + alphamlr)*(5.0 + alphamlr)*(4.0 + alphamlr)/ &
!     &            (pi**2*(3.0 + alphamlr)*(2.0 + alphamlr)*(1.0 + alphamlr))
         g1mlr = (6.0 + alphamlr)*(5.0 + alphamlr)*(4.0 + alphamlr)/ &
     &            ((3.0 + alphamlr)*(2.0 + alphamlr)*(1.0 + alphamlr))
         g1smlr = (6.0 + alphasmlr)*(5.0 + alphasmlr)*(4.0 + alphasmlr)/ &
     &            ((3.0 + alphasmlr)*(2.0 + alphasmlr)*(1.0 + alphasmlr))
       ENDIF
      ENDIF

      IF ( lzr > 1 .and. imurain == 3 ) THEN ! { RAIN SHAPE PARAM
      
      
!      CALL cld_cpu('Z-MOMENT-1r')  
          il = lr
          DO mgs = 1,ngscnt
          

         IF ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) THEN ! .or. qx(mgs,il) <= qxmin(il)  THEN
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3  THEN
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ELSEIF ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3   THEN
         
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         ENDIF
         ENDIF

         IF ( .false. .and. zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
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
           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
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
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

! check for artificial breakup (rain larger than allowed max size)
        IF (  (xv(mgs,il) .gt. xvmx(il) .or. (ioldlimiter >= 2 .and. xv(mgs,il) .gt. xvmx(il)/8.) )) THEN
          tmp = cx(mgs,il)
          IF ( ioldlimiter >= 2 ) THEN ! MY-style active breakup
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          ELSE ! simple cutoff 
            xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          ENDIF
            !xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            !cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))

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
         ! the moments are not matched correctly, so we compute G from the moments instead so that the dZ/dt rates
         ! stay consistent with dN/dt and dq/dt.
           IF ( alp >= rnumax - 0.01 ) THEN
!             g1x(mgs,il) = 6**2*zx(mgs,il)/(cx(mgs,il)*(pi*xv(mgs,lr))**2)
!             g1x(mgs,il) = xdn(mgs,il)*zx(mgs,il)*cx(mgs,il)/((rho0(mgs)*qx(mgs,lr))**2)
             g1x(mgs,il) = (pi*xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((6.*rho0(mgs)*qx(mgs,il))**2)
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
           
           gf1palp(mgs) = y

!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma_sp(alpha(mgs,lr) + 1.)
           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))

           IF ( imurain == 3 .and. izwisventr == 2 ) THEN

           tmp = alpha(mgs,lr) + 1.5 + br/6.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 1.5 + br/6.)/Gamma_sp(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/(y*(alpha(mgs,lr) + 1.)**((1.+br)/6. + 1./3.))
           
! This whole section is imurain == 3, so this branch never runs
!           ELSEIF ( imurain == 1 .and.  iferwisventr == 2 ) THEN
!
!           tmp = alpha(mgs,lr) + 2.5 + br/2.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 1.5 + br/6.)/Gamma_sp(alpha(mgs,lr) + 1.)
!           ventrxn(mgs) = x/y
           
           
           ENDIF
           
           ENDIF
          ENDIF
          
          ENDIF
          
          ENDDO
!        CALL cld_cpu('Z-MOMENT-1r')  
        ENDIF ! }
        
      ENDIF ! ipconc >= 6

!  Find shape parameters for graupel and hail
      IF ( ipconc .ge. 6 ) THEN
            
        DO il = lr,lhab
          
        ! set base values of g1x
          IF ( (.not. ( il == lr .and. imurain == 3 )) .and. ( il == lr .or. il == lh .or. il == lhl .or. il == lf ) ) THEN
          DO mgs = 1,ngscnt
            g1x(mgs,il) = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
          ENDDO
          ENDIF
        
        IF ( lz(il) .gt. 1   .and. ( .not. ( il == lr .and. imurain == 3 )) ) THEN
        
        DO mgs = 1,ngscnt


         IF ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) THEN ! .or. qx(mgs,il) <= qxmin(il) ) THEN
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3 ) THEN
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           zx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ELSEIF ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3  ) THEN
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ENDIF
         ENDIF

         IF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
        
        IF ( qx(mgs,il) .gt. qxmin(il) ) THEN

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        IF ( xv(mgs,il) .lt. xvmn(il) ) THEN
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
        ENDIF

          IF ( zx(mgs,il) > zxmin .and. cx(mgs,il) <= cxmin ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
!            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)
            IF ( cx(mgs,il) < cxmin ) THEN
            ! if resulting concentration is too small, then zero out
              cx(mgs,il) = 0.0
              zx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

              qx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
!            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            g1 = (6.0 + alphamax)*(5.0 + alphamax)*(4.0 + alphamax)/ &
     &            ((3.0 + alphamax)*(2.0 + alphamax)*(1.0 + alphamax))
            zx(mgs,il) = Max(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

            IF ( zx(mgs,il) <= zxmin ) THEN
            ! if resulting reflectivity is still too small, then zero out
              cx(mgs,il) = 0.0
              zx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

              qx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
!               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           ELSE
          
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          IF ( zx(mgs,il) .gt. 0. ) THEN
           
!            rdi = z*(pi/6.*1000.)**2*chw/((rho0(mgs)*qr)**2)
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           alp = Max( alphamin, Min( alphamax, alp ) )
           
         IF ( newton ) THEN
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = alp + ( galpha(alp) - rdi )/dgalpha(alp)
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
           
         ELSE
           DO i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,il))/(Abs(alpha(mgs,il))+1.e-5) .lt. 1. ) EXIT
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
          ENDIF


! check for artificial breakup (graupel/hail larger than allowed max size)
        IF ( imaxdiaopt == 1 .or. il /= lr ) THEN
          xvbarmax = xvmx(il) 
        ELSEIF ( imaxdiaopt == 2 ) THEN ! test against maximum mass diameter
          xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ELSEIF ( imaxdiaopt == 3 ) THEN ! test against mass-weighted diameter
          xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ELSE
          xvbarmax = xvmx(il) 
        ENDIF

        IF (  xv(mgs,il) .gt. xvbarmax .or. (il == lr .and. ioldlimiter >= 2 .and. xv(mgs,il) .gt. xvmx(il)/8.)) THEN
          tmp = cx(mgs,il)
          IF( ioldlimiter >= 2 .and. il == lr) THEN ! MY-style drop limiter for rain
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          ELSE
            xv(mgs,il) = Min( xvbarmax, Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          ENDIF
          IF ( tmp < cx(mgs,il) ) THEN ! artificial breakup has happened, so need to adjust reflectivity and find new shape parameter
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
!             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            ! check if incoming zx is consistent
            ! Z from incoming cx, qx, and alpha
            tmpz = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/tmp
            IF ( tmpz > zx(mgs,il) ) THEN
              tmpc = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/zx(mgs,il)
              cx(mgs,il) = Max(cx(mgs,il), tmpc)
              ! find cx that gives zx
            ENDIF
            zx(mgs,il) = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/cx(mgs,il)
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

            qr  = qx(mgs,il)
            chw = cx(mgs,il)
            z   = zx(mgs,il)

            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
 
           IF ( ( lrescalelow(il) .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) THEN



            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( lrescalelow(il) .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) .and. &
                     .not. ( il == lr .and. .not. rescale_low_alphar ) ) THEN ! alpha = alphamin, so reset Z to prevent growth in C
             wtest = .false.
             IF ( irescalerainopt == 0 ) THEN
               wtest = .false.
             ELSEIF ( irescalerainopt == 1 ) THEN
               wtest = qx(mgs,lc) > qxmin(lc) 
             ELSEIF ( irescalerainopt == 2 ) THEN
               wtest = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             ELSEIF ( irescalerainopt == 3 ) THEN
               wtest = temcg(mgs) > rescale_tempthresh .and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             ENDIF
             
             IF ( il == lr .and. ( wtest ) ) THEN
!             IF ( temcg(mgs) > 0.0 .and. il == lr .and. qx(mgs,lc) > qxmin(lc) ) THEN
             ! certain situations where rain number is adjusted instead of Z. Helps avoid rain being 'zapped' by autoconverted 
             ! drops (i.e., favor preserving Z when alpha tries to go negative)
             chw = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2 ! g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1
             cx(mgs,il) = chw
             an(igs(mgs),jy,kgs(mgs),ln(il)) = chw
             ELSE
             
             ! Usual resetting of reflectivity moment to force consisntency between Q, N, Z, and alpha when alpha = alphamin
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             ENDIF
            ENDIF
           ENDIF
          
          
         ! set g1x to use as G factor later. If alpha is in the range ( rnumin < alpha < rnumax ), then 
         ! this will be the same as computing G from alpha.  If alpha = rnumax, however, it probably means that
         ! the moments are not matched correctly, so we compute G from the moments instead so that the dZ/dt rates
         ! stay consistent with dN/dt and dq/dt.
!          g1x(mgs,il) = zx(mgs,il)*chw*(pi*xdn(mgs,il))**2/(6.*qr*dn(igs(mgs),jy,kgs(mgs)))**2
!          g1x(mgs,il) = g1 ! zx(mgs,il)*cx(mgs,il)/(qr)**2
           IF ( alp >= alphamax - 0.5 ) THEN
!             g1x(mgs,il) = 6**2*zx(mgs,il)/(cx(mgs,il)*(pi*xv(mgs,lr))**2)
!             g1x(mgs,il) = (xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((rho0(mgs)*qx(mgs,il))**2)
             g1x(mgs,il) = (pi*xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((6.*rho0(mgs)*qx(mgs,il))**2)
           ELSE
             g1x(mgs,il) = g1
           ENDIF
          
           ENDIF
          
!          IF ( ny .eq. 2 ) THEN
!          IF ( qr .gt. 1.e-3 ) THEN
!           write(0,*) 'alphah at nstep,i,k = ',dtp*(nstep-1),igs(mgs),kgs(mgs),alpha(mgs,il),qr*1000.
!          ENDIF
!          ENDIF
          
           
           ENDIF ! .true.

          IF ( il == lr ) THEN
           
!           tmp = alpha(mgs,lr) + 4./3.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!           tmp = alpha(mgs,lr) + 1.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma_sp(alpha(mgs,lr) + 1.)
!           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))


           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           gf1palp(mgs) = y

           IF (   iferwisventr == 2 ) THEN
!        ventrn =  Gamma(alphar + 2.5 + br/2.)/Gamma(alphar + 1.) ! adapted from Wisner et al. 1972
           tmp = alpha(mgs,lr) + 2.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami


           ventrxn(mgs) = x/y
           
           ENDIF
           
          ENDIF ! il==lr
 
          
          ELSE ! below mass threshold
!             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
!             z1 = g1*rho0(mgs)**2*(qr)*qr/chw
!             z  = 1.e18*z1*(6./(pi*1000.))**2
!             z  = z1*(6./(pi*1000.))**2
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
          ENDIF ! ( qx(mgs,il) .gt. qxmin(il) )
        
        
        
!        ENDIF
        ENDDO ! mgs

!         CALL cld_cpu('Z-DELABK')  
        
!        IF ( il == lr ) THEN
!          xnutmp = (alpha(mgs,il) - 2.)/3.
!           da0lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
!        ENDIF
        
        IF ( .not. ( il == lr .and. imurain == 3 ) ) THEN
!          CALL cld_cpu('Z-DELABK')  
        DO mgs = 1,ngscnt
          IF ( qx(mgs,il) > qxmin(il) ) THEN
          xnutmp = (alpha(mgs,il) - 2.)/3.
          
!          IF ( .true. ) THEN
          DO ic = lc,lh-1 ! lhab
           IF ( il .ne. ic .and.  qx(mgs,ic) .gt. qxmin(ic)) THEN
             xnuc = xnu(ic)
             IF ( ic == lc .and. idiagnosecnu > 0 ) xnuc = alpha(mgs,lc) ! alpha for droplets is actually nu
             IF ( il /= lr .and. ic == lr .and. lzr > 1 ) THEN
               IF ( imurain == 3 ) THEN
                 xnuc = alpha(mgs,lr) ! alpha is nu already
               ELSE
                 xnuc = ( alpha(mgs,lr) - 2. )/3. ! convert alpha to nu
               ENDIF
             ENDIF
                                 ! delabk(ba,bb,nua,nub,mua,mub,k), where a (il)  is collector and b (ic) is collected
             IF ( .false. ) THEN
             dab0lh(mgs,ic,il) =  delabk(bb(ic), bb(il), xnuc, xnutmp, xmu(ic), xmu(il), 0) !dab0(il,ic)
             dab1lh(mgs,ic,il) =  delabk(bb(ic), bb(il), xnuc, xnutmp, xmu(ic), xmu(il), 1) !dab1(il,ic)
             dab0lh(mgs,il,ic) =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 0) !dab0(il,ic)
             dab1lh(mgs,il,ic) =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 1) !dab1(il,ic)
             ELSE ! use lookup table -- not interpolating yet because table resolution of 0.05 is good enough
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               IF ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) THEN
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               ELSE ! IF ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               ENDIF
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)

!               tmp1 = dab0lu(j,i,ic,il)
!               tmp2 = dab1lu(j,i,ic,il)
!               tmp3 = dab0lu(i,j,il,ic)
!               tmp4 = dab1lu(i,j,il,ic)
!               tmp5 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(ic), xmu(il), 0) !dab0(il,ic)
!               tmp6 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(ic), xmu(il), 1) !dab1(il,ic)
!               tmp5 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 0) !dab0(il,ic)
!               tmp6 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 1) !dab1(il,ic)
               
               IF ( .false. .and. ny <= 2 ) THEN
                 write(0,*)
                 write(0,*) 'bb: ', bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic)
                 write(0,*) 'il,ic = ',il,ic,alpha(mgs,il),i,xnuc,alp,j
                 write(0,*) 'dab0lh,tmp1 = ',dab0lh(mgs,ic,il),tmp1
                 write(0,*) 'dab1lh,tmp2 = ',dab1lh(mgs,ic,il),tmp2
                 write(0,*) 'dab0lh,tmp3 = ',dab0lh(mgs,il,ic),tmp3,tmp5
                 write(0,*) 'dab1lh,tmp4 = ',dab1lh(mgs,il,ic),tmp4,tmp6
               
               ENDIF
             
             ENDIF
             
           ENDIF
          ENDDO

!          ENDIF
           
             da0lx(mgs,il) = delbk(bb(il), xnutmp, xmu(il), 0)
           IF ( il .eq. lh ) THEN
             da0lh(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
            IF ( lzr > 1 ) THEN
             rzxh(mgs) = 1.
            ELSE
             rzxh(mgs) = ((4. + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr)))
            ENDIF
            
            IF ( lzhl < 1 ) THEN
              rzxhlh(mgs) = rzxhl(mgs)/(((4. + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr))))
            ENDIF
           ELSEIF ( il .eq. lf ) THEN
             da0lf(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
            IF ( lzr > 1 ) THEN
             rzxf(mgs) = 1.
            ELSE
             rzxf(mgs) = ((4.0 + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr)))
            ENDIF
           ELSEIF ( il .eq. lhl ) THEN
             da0lhl(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
            IF ( lzr > 1 ) THEN
             rzxhl(mgs) = 1.
            ELSE
             rzxhl(mgs) = ((4.0 + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr)))
            ENDIF
           ELSEIF ( il == lr ) THEN
             xnutmp = (alpha(mgs,il) - 2.)/3.
             da0lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
             da1lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 1)
           ENDIF
          
          ENDIF ! ( qx(mgs,il) > qxmin(il) )
        ENDDO ! mgs
!          CALL cld_cpu('Z-DELABK')  
        ENDIF ! il /= lr

!         CALL cld_cpu('Z-DELABK')  
        
        ENDIF ! lz(il) .gt. 1
        
        ENDDO ! il
          
      ENDIF ! ipconc .ge. 6

!      CALL cld_cpu('Z-MOMENT-1')  

!
!  set some values for ice nucleation
!
      do mgs = 1,ngscnt
      kp1 = Min(nz, kgs(mgs)+1 )
!      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1)   &
!     &                  +w(igs(mgs),jgs,kgs(mgs)))

      
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs))   &
     &                    +w(igs(mgs),jgs,kgsm(mgs)))
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
!      SUBROUTINE setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno, &
!     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,            &
!     &                 ipconc1,ndebug1,ngs,nz,kgs,cwnccn,fadvisc,   &
!     &                 cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
!     &                 itype1a,itype2a,temcg,infdo,alpha)


      infdo = 1
      IF ( rimdenvwgt > 0 ) infdo = 1

!      SUBROUTINE setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno, &
!     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,            &
!     &                 ipconc1,ndebug1,ngs,nz,kgs,cwnccn,fadvisc,   &
!     &                 cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
!     &                 itype1a,itype2a,temcg,infdo,alpha,axh,bxh,axhl,bxhl)
      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,   &
     &                 ipconc,ndebug,ngs,nz,igs,kgs,cwnccn,fadvisc,   &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,   &
     &                 itype1,itype2,temcg,infdo,alpha,alpha,axx,bxx,0) ! ,cdh,cdhl)
!     &                 itype1,itype2,temcg,infdo,alpha,axh,bxh,axhl,bxhl) ! ,cdh,cdhl)

       IF ( isedonly == 4 ) THEN
         goto 1234
       ENDIF


       IF ( lwsm6 .and. ipconc == 0 ) THEN
         tmp = Max(qxmin(lh), qxmin(ls))
         DO mgs = 1,ngscnt
           total = qx(mgs,lh) + qx(mgs,ls)
           IF ( total > tmp ) THEN
             vt2ave(mgs) = (qx(mgs,lh)*vtxbar(mgs,lh,1) + qx(mgs,ls)*vtxbar(mgs,ls,1))/total
           ELSE
             vt2ave(mgs) = 0.0
           ENDIF
         ENDDO
       ENDIF


!
!  Set number concentrations (need xdia from setvt)
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set concentration'
      IF ( ipconc .lt. 1 ) THEN
         cina(1:ngscnt) = cx(1:ngscnt,li)
      ENDIF
      if ( ipconc .lt. 5 ) then
      do mgs = 1,ngscnt


      IF ( ipconc .lt. 3 ) THEN
!      cx(mgs,lr) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lh) )  then
!      cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
!      xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
      end if
      ENDIF

      IF ( ipconc .lt. 4 ) THEN
!      tmp = cx(mgs,ls)
!      cx(mgs,ls) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) )  then
!      cx(mgs,ls) = cno(ls)*xdia(mgs,ls,1)
!      xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
      end if
      ENDIF ! ( ipconc .lt. 4 )

      IF ( ipconc .lt. 5 ) THEN


!      cx(mgs,lh) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) )  then
!      cx(mgs,lh) = cno(lh)*xdia(mgs,lh,1)
!      xv(mgs,lh) = Max(xvmn(lh), rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh)) )
!      xdia(mgs,lh,3) = (xv(mgs,lh)*6./pi)**(1./3.) 
      end if

      ENDIF ! ( ipconc .lt. 5 )

      end do
      end if
      
      IF ( ipconc .ge. 2 ) THEN
      DO mgs = 1,ngscnt
        
        rb(mgs) = 0.5*xdia(mgs,lc,1)*(1./(1.+alpha(mgs,lc)))**(1./6.)
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
! Set mean particle surface area
!
      IF ( ipelec .ge. 1 .and. ieopt .eq. 2) THEN
!      DO il = li,lhab
       do mgs = 1,ngscnt
         xsfca(mgs,ls) = 2*pi*cx(mgs,ls)*xdia(mgs,ls,2)
       enddo
!      ENDDO
       ENDIF
!
!
!              
!
!  maximum depletion tendency by any one source
!
!
      if( ndebug .ge. 0 ) THEN
!mpi!        write(0,*) 'Set depletion max/min1'
      endif
      do mgs = 1,ngscnt
      qvimxd(mgs) = 0.70*(qx(mgs,lv)-qis(mgs))*dtpinv ! depletion by all vap. dep to ice.
      
      IF ( qx(mgs,lc) < qxmin(lc) ) qvimxd(mgs) = 0.99*(qx(mgs,lv)-qis(mgs))*dtpinv ! this makes virtually no difference whatsoever, but what the heck
      
      qvimxd(mgs) = max(qvimxd(mgs), 0.0)

      frac = 0.1d0
      qimxd(mgs)  = frac*qx(mgs,li)*dtpinv
      qcmxd(mgs)  = frac*qx(mgs,lc)*dtpinv
      qrmxd(mgs)  = frac*qx(mgs,lr)*dtpinv
      qsmxd(mgs)  = frac*qx(mgs,ls)*dtpinv
      qhmxd(mgs)  = frac*qx(mgs,lh)*dtpinv
      IF ( lf > 1 ) qfmxd(mgs)  = frac*qx(mgs,lf)*dtpinv
      IF ( lhl > 1 ) qhlmxd(mgs)  = frac*qx(mgs,lhl)*dtpinv
      end do
!
      if( ndebug .ge. 0 ) THEN
!mpi!        write(0,*) 'Set depletion max/min2'
      endif

      do mgs = 1,ngscnt
!
      if ( qx(mgs,lc) .le. qxmin(lc) ) then
      ccmxd(mgs)  = 0.20*cx(mgs,lc)*dtpinv
      else
      IF ( ipconc .ge. 2 ) THEN
        ccmxd(mgs)  = frac*cx(mgs,lc)*dtpinv
      ELSE
        ccmxd(mgs)  = frac*qx(mgs,lc)/(xmas(mgs,lc)*rho0(mgs)*dtp)
      ENDIF
      end if
!
      if ( qx(mgs,li) .le. qxmin(li) ) then
      cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      else
      IF ( ipconc .ge. 1 ) THEN
        cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      ELSE
        cimxd(mgs)  = frac*qx(mgs,li)/(xmas(mgs,li)*rho0(mgs)*dtp)
      ENDIF
      end if
!
!
      crmxd(mgs)  = 0.10*cx(mgs,lr)*dtpinv
      csmxd(mgs)  = frac*cx(mgs,ls)*dtpinv
      chmxd(mgs)  = frac*cx(mgs,lh)*dtpinv

      ccmxd(mgs)  = frac*cx(mgs,lc)*dtpinv
      cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      crmxd(mgs)  = frac*cx(mgs,lr)*dtpinv
      csmxd(mgs)  = frac*cx(mgs,ls)*dtpinv
      chmxd(mgs)  = frac*cx(mgs,lh)*dtpinv

      qxmxd(mgs,lv) = Max(0.0, 0.1*(qx(mgs,lv) - qvs(mgs))*dtpinv)

      DO il = lc,lhab
       qxmxd(mgs,il) = frac*qx(mgs,il)*dtpinv
       cxmxd(mgs,il) = frac*cx(mgs,il)*dtpinv
      ENDDO

      end do




      IF ( ipconc >= 6 ) THEN
      frac = 0.4d0
      zxmxd(:,:) = 0.0
      DO il = lr,lhab
       IF ( lz(il) > 0 .or. ( il == lr ) ) THEN
         DO mgs = 1,ngscnt
           zxmxd(mgs,il) = frac*zx(mgs,il)*dtpinv
         ENDDO
       ENDIF
      ENDDO
      ENDIF



      IF ( ipelec .ge. 1 ) THEN
      DO mgs = 1,ngscnt
      cionpmxd(mgs) = frac*cionp(mgs)*dtpinv
      cionnmxd(mgs) = frac*cionn(mgs)*dtpinv
      
      if (largeion) then
       clionpmxd(mgs) = frac*clionp(mgs)*dtpinv
       clionnmxd(mgs) = frac*clionn(mgs)*dtpinv
      endif
      
      ENDDO
      ENDIF

      if ( ipelec .gt. 1 ) then
!
!  vertical component of e
!
      do mgs = 1,ngscnt
      dezcomp(mgs) = Max(-1.00e5,   &
     &   Min(elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs)),1.00e5 )) 
      end do
!       write(0,*) 'dezcomp',MAXVAL(dezcomp),MAXLOC(dezcomp),MINVAL(dezcomp),MINLOC(dezcomp)
!       write(0,*) 'ELEC3',MAXVAL(elec),MAXLOC(elec),MINVAL(elec),MINLOC(elec)



      IF ( .false. ) THEN
        DO mgs = 1,ngscnt
          IF ( Abs(scx(mgs,lh)) > 0.5e-9 .and. qx(mgs,lh) > 0.1e-3 .and. Abs( elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs))) > 30.e3 ) THEN
          tmp = 4.0*elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs))*scx(mgs,lh)*gamma_sp(alpha(mgs,lh)+1.) &
               *gamma_sp(alpha(mgs,lh)+4.-bxx(mgs,lh)) /   &
             ( axx(mgs,lh)*cdx(lh)*cx(mgs,lh)*pi*rho0(mgs)*gamma_sp(4.+alpha(mgs,lh)) * gamma_sp(3.0+alpha(mgs,lh))* &
              xdia(mgs,lh,1)**(2.0+bxx(mgs,lh)) )
           write(0,*) 'vtg, vtg-star, qg, cg, bxh = ',vtxbar(mgs,lh,1), tmp, qx(mgs,lh), cx(mgs,lh),bxx(mgs,lh)
           write(0,*) ' numer: ',4.0*elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs))*scx(mgs,lh)*gamma_sp(alpha(mgs,lh)+4.-bxx(mgs,lh))
           write(0,*) ' numer parts: ',elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs)), scx(mgs,lh),gamma_sp(alpha(mgs,lh)+4.-bxx(mgs,lh))
           write(0,*) ' denom: ', axx(mgs,lh),cdx(lh),cx(mgs,lh),pi,rho0(mgs),gamma_sp(4.+alpha(mgs,lh)) ,&
              gamma_sp(3.0+alpha(mgs,lh)), xdia(mgs,lh,1)**(2.0+bxx(mgs,lh)) 
              
            tmp = xdia(mgs,lh,1)*elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs))*scx(mgs,lh) / &
               (2560.*cdx(lh)* cx(mgs,lh)* sqrt(rho0(mgs)*pi)* (xdia(mgs,lh,1)/2. )**(3.5))
            write(0,*) ' rawlins: ',tmp
        
         ENDIF
        ENDDO
      ENDIF
!
!  others...
!
      if ( nonigrd .ge. 2 ) then
      
      sctem(:) = 0.0
      
      do mgs = 1,ngscnt
!
!
!  temperature factor:
!      if trever = -21: reduces to jayratene
!      if trever = -10: reduces to takahashi
!
      IF (isaund .eq. 1 .or. isaund .eq. 11) THEN
        trev = trever
      ELSEIF (isaund .eq. 2 .or. isaund .eq. 12) THEN  
        ewtmp = 0.5*1.e3*rho0(mgs)*qx(mgs,lc)
        trev = max(-23.8, (-15.06 * ewtmp - 7.38))  ! Saunders reversal temp.
      ELSE
        trev = trever
      ENDIF

       IF (ftauopt == 2 ) THEN
! shift the charge function ftau so that zero crossing matches local
! reversal temperature.
        IF (temcg(mgs) .le. -7.5) THEN
          sctem(mgs) =  ftau( temcg(mgs) - 21.0 - trev )
        ELSE
          sctem(mgs) = ftau(-7.5 - 21.0 - trev)*(-temcg(mgs)*(1./7.5))
        ENDIF
       ELSEIF ( ftauopt == 1 ) THEN
! do the ziegler/straka shift...
        IF (temcg(mgs) .lt. 0.0) THEN
          sctem(mgs) = ftau(max(-25.0,((-21.0)/trev)*temcg(mgs)) )
        ENDIF
       
       ENDIF

!      temele = max(-25.0,((-21.0)/trever)*temcg(mgs))
!
!      sctem(mgs)  =
!     >   (-1.7e-05)*(temele**3)
!     > + (-3.0e-03)*(temele**2)
!     > + (-5.0e-02)*(temele**1)
!     > + ( 1.3e-01)
!
!  critical (threshold) liquid water content factor (1.0e-04 kg/m**3)
!
      if ( temcg(mgs) .ge. trev ) then
      dellwc(mgs) = (1000.0)*(rho0(mgs)*qx(mgs,lc)-(1.0e-04))
      end if
      if ( temcg(mgs) .lt. trev ) then
      dellwc(mgs) = (1000.0)*(rho0(mgs)*qx(mgs,lc))
      end if
      if ( qx(mgs,lc) .lt. 1.e-6 ) dellwc(mgs) = 0.0
!
! parabolic function to reduce charging at low temperature
!
!      IF ( temcg(mgs) .lt. -30 .and. temcg(mgs) .gt. -43.0 ) THEN
!        sctem(mgs) = sctem(mgs)*(1.0 - ((temcg(mgs)+30.0)/(43.0 - 30.0))**2)
      
      IF ( isctemopt == 1 ) THEN
        IF ( temcg(mgs) .lt. -30 .and. temcg(mgs) .gt. -40.0 ) THEN
          sctem(mgs) = sctem(mgs)*(1.0 - ((temcg(mgs)+30.0)/(40.0 - 30.0))**2)
!        ELSEIF ( temcg(mgs) .le. -43.0 ) THEN
        ELSEIF ( temcg(mgs) .le. -40.0 ) THEN
          sctem(mgs) = 0.0
        END IF
      ELSE
        IF ( temcg(mgs) >= tcc ) THEN
          fac = 1.0
        ELSEIF ( temcg(mgs) < tcc .and. temcg(mgs) .gt. nic_min_temp ) THEN
!         fac = 1.0 - ((temcg+30.0)/(-tmin - 30.0))**2
! Cosine funtion roll-off
          fac = 0.5 *(1. + Cos(pi *(((temcg(mgs) - tcc)/(-nic_min_temp + tcc)))))
        ELSE
          fac = 0.0
        END IF
        
        sctem(mgs) = fac*sctem(mgs)
      ENDIF
      

      end do
      end if
!
      end if
!



    ! default factors between mean volume and maximum mass volume
      maxmassfac(lc) = ( (2. + 3.*(1. + xnu(lc)) )**3/( 3.*(1. + xnu(lc)) ) )
      maxmassfac(li) = ( (2. + 3.*(1. + xnu(li)) )**3/( 3.*(1. + xnu(li)) ) )

      IF ( imurain == 3 ) THEN
        maxmassfac(lr) = ( (2. + 3.*(1. + xnu(lr)) )**3/( 3.*(1. + xnu(lr)) ) )
      ELSE
        maxmassfac(lr) =  (3.0 + alphar)**3/    &
     &                  ((3.+alphar)*(2.+alphar)*(1. + alphar) )
      ENDIF

      IF ( imusnow == 3 ) THEN
        maxmassfac(ls) = ( (2. + 3.*(1. + alphas) )**3/( 3.*(1. + alphas) ) )
      ELSE
        maxmassfac(ls) =  (3.0 + alphas)**3/    &
     &                  ((3.+alphas)*(2.+alphas)*(1. + alphas) )
      ENDIF
      
        maxmassfac(lh) =  (3.0 + alphah)**3/    &
     &                  ((3.+alphah)*(2.+alphah)*(1. + alphah) )

       IF ( lhl > 1 ) THEN
        maxmassfac(lhl) =  (3.0 + alphahl)**3/    &
     &                  ((3.+alphahl)*(2.+alphahl)*(1. + alphahl) )
       ENDIF
      

! calculate maximum mass diameters
    DO il = lc,lhab
      
      
      IF ( il == lc .or. il == li .or. (il == lr .and. imurain == 3 ) .or. (il == ls .and. imusnow == 3 ) ) THEN
       DO mgs = 1,ngscnt
        IF ( qx(mgs,il) > qxmin(il) ) THEN
         IF ( lz(il) < 1 ) THEN
           ! 2-moment
!           xdiamxmas(mgs,il) = xdia(mgs,il,3)*maxmassfac(il)
           xdiamxmas(mgs,il) = maxmassfac(il)
         ELSE
           ! 3-moment
           xdiamxmas(mgs,il) = xdia(mgs,il,1)*(3. + alpha(mgs,il))
!           xdiamxmas(mgs,il) = ( (2. + 3.*(1. + alpha(mgs,il)) )**3/( 3.*(1. + alpha(mgs,il)) ) )
         ENDIF
        ELSE
          xdiamxmas(mgs,il) = 1. ! xdia(mgs,il,3) ! temporary value!!!
        ENDIF
       ENDDO
      
      ELSE
       DO mgs = 1,ngscnt
        IF ( qx(mgs,il) > qxmin(il) ) THEN
         IF ( lz(il) < 1 ) THEN
           ! 2-moment
!           xdiamxmas(mgs,il) = xdia(mgs,il,3)*maxmassfac(il)
           xdiamxmas(mgs,il) = maxmassfac(il)
         ELSE
!          xdiamxmas(mgs,il) = xdia(mgs,il,3)*(3.0 + alpha(mgs,il))*    &
!     &     ((3.+alpha(mgs,il))*(2.+alpha(mgs,il))*(1. + alpha(mgs,il)) )**(-1./3.)
          xdiamxmas(mgs,il) = (3.0 + alpha(mgs,il))**3/    &
     &     ((3.+alpha(mgs,il))*(2.+alpha(mgs,il))*(1. + alpha(mgs,il)) )
          ENDIF
        ELSE
          xdiamxmas(mgs,il) = 1.0 ! xdia(mgs,il,3) ! temporary value!!!
        ENDIF
       ENDDO
      ENDIF
    
    ENDDO ! il

       DO mgs = 1,ngscnt
          DO il = lh,lhab ! graupel and hail only (and frozen drops)
            
            vshdgs(mgs,il) = vshd ! base value
            
            IF ( qx(mgs,il) > qxmin(il) .and. ivshdgs > 0 ) THEN
              
              ! tmpdiam is weighted diameter of d^(shedalp-1), so for shedalp=3, this is the area-weighted diameter or maximum mass diameter.
             ! tmpdiam = (shedalp+alpha(mgs,il))*xdia(mgs,il,1) ! *( xdn(mgs,il)/917. )**(1./3.) ! erm added density factor for equiv. solid ice sphere 10.12.2015
              tmpdiam = (shedalp+0.0)*xdia(mgs,il,1) ! *( xdn(mgs,il)/917. )**(1./3.) ! erm added density factor for equiv. solid ice sphere 10.12.2015
              ! imltshddmr
              IF ( tmpdiam > sheddiam0 ) THEN
                vshdgs(mgs,il) = 0.523599*(1.5e-3)**3/massfacshr ! 1.5mm drops from very large ice
              ELSEIF ( tmpdiam > sheddiam ) THEN ! intermediate size
                vshdgs(mgs,il) = 0.523599*(3.0e-3)**3/massfacshr ! 3.0mm drops from medium-large ice
              ELSE
!                vshdgs(mgs,il) = Min( xvmx(lr), xv(mgs,il)*xdn(mgs,il)*0.001 ) ! size of drop from melted mean ice particle
                vshdgs(mgs,il) = Min( xvmx(lr), 6./pi*xdn(mgs,il)*0.001*tmpdiam**3 )/massfacshr ! size of drop from melted mean ice particle; 0.001 is 1/rhow
              ENDIF
            ENDIF
          ENDDO
       ENDDO

!
!
!  microphysics source terms (1/s) for mixing ratios 
!
!
!
!  Collection efficiencies:
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set collection efficiencies'
!
       exy(:,:,:) = 0.0
      do mgs = 1,ngscnt
!
!
!
      qcwresv(mgs) = 0.0
      ccwresv(mgs) = 0.0
      
      erw(mgs) = 0.0
      esw(mgs) = 0.0
      ehw(mgs) = 0.0
      efw(mgs) = 0.0
      ehlw(mgs) = 0.0
!      ehxw(mgs) = 0.0
!
      err(mgs) = 0.0
      esr(mgs) = 0.0
      il2(mgs) = 0
      il3(mgs) = 0
      ehr(mgs) = 0.0
      ehlr(mgs) = 0.0
!      ehxr(mgs) = 0.0
!
      eri(mgs) = 0.0
      esi(mgs) = 0.0
      ehi(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehi*ehiclsn
      ehis(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehi*ehiclsn
      ehli(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehli*ehliclsn
      ehlis(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehli*ehliclsn
!      ehxi(mgs) = 0.0
!
      ers(mgs) = 0.0
      ess(mgs) = 0.0
      ehs(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehs*ehsclsn
      ehsfac(mgs) = 1.0 ! factor based on ice saturation
      ehls(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehls*ehlsclsn
      ehscnv(mgs) = 0.0
!      ehxs(mgs) = 0.0
!
      eiw(mgs) = 0.0
      eii(mgs) = 0.0
      efsclsn(mgs) = 0.0
      eficlsn(mgs) = 0.0
      efi(mgs) = 0.0
      efis(mgs) = 0.0
      efs(mgs) = 0.0
      efr(mgs) = 0.0
      efw(mgs) = 0.0
      ehsclsn(mgs) = 0.0
      ehiclsn(mgs) = 0.0
      ehlsclsn(mgs) = 0.0
      ehliclsn(mgs) = 0.0
      esiclsn(mgs) = 0.0


! reserve droplets
         IF ( exwmindiam > 0 .and. qx(mgs,lc) > qxmin(lc) ) THEN
           tmp = cx(mgs,lc)*Exp(- (exwmindiam/xdia(mgs,lc,1))**3 )
           ccwresv(mgs) =  Min( cx(mgs,lc), Max( 2.e6,  cx(mgs,lc) -  tmp ) )
           
           tmp = cx(mgs,lc) - ccwresv(mgs)

           volt = pi/6.*(exwmindiam)**3
           qcwresv(mgs) = qx(mgs,lc) - tmp*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
           
           
           IF ( .false. .and. qx(mgs,lc) > 0.1e-3 ) THEN
           
             write(0,*) 'cx,qx,crsv,qrsv = ',cx(mgs,lc),qx(mgs,lc),ccwresv(mgs),qcwresv(mgs)
           
           ENDIF

         ENDIF


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


      igwr(mgs) = 1
!      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
!         rwrad = 0.5*xdia(mgs,lr,1)
! setting erw = 1 always, so now use igwr for graupel
      IF ( qx(mgs,lh) .gt. qxmin(lh) ) THEN
         rwrad = 0.5*xdia(mgs,lh,3)  ! changed to mean volume diameter (10/6/06)
      DO il = 1,6
         IF ( rwrad .ge. 1.e-6*grad(il,1) ) igwr(mgs) = il
      ENDDO
      ENDIF

      IF ( lf > 1 ) THEN
      ifwr(mgs) = 1
      IF ( qx(mgs,lf) .gt. qxmin(lf) ) THEN
         rwrad = 0.5*xdia(mgs,lf,3)  ! changed to mean volume diameter (10/6/06)
      DO il = 1,6
         IF ( rwrad .ge. 1.e-6*grad(il,1) ) ifwr(mgs) = il
      ENDDO
      ENDIF
      ENDIF

      IF ( lhl .gt. 1 ) THEN ! hail is turned on
      ihlr(mgs) = 1
      IF ( qx(mgs,lhl) .gt. qxmin(lhl) ) THEN
         rwrad = 0.5*xdia(mgs,lhl,3)  ! changed to mean volume diameter (10/6/06)
      DO il = 1,6
         IF ( rwrad .ge. 1.e-6*grad(il,1) ) ihlr(mgs) = il
      ENDDO
      ENDIF
      ENDIF

!
!
!  Ice-Ice: Collection (cxc) efficiencies
!
!
      if ( qx(mgs,li) .gt. qxmin(li) ) then
!      IF ( ipconc .ge. 14 ) THEN
!       eii(mgs)=0.1*exp(0.1*temcg(mgs))
!       if ( temg(mgs) .lt. 243.15 .and. qx(mgs,lc) .gt. 1.e-6 ) then
!        eii(mgs)=0.1
!       end if
!      
!      ELSE
        eii(mgs) = exp(0.025*Min(temcg(mgs),0.0))  ! alpha1 from LFO83 (21)
!      ENDIF
      if ( temg(mgs) .gt. 273.15 ) eii(mgs) = 1.0
      end if
!
!
!
!  Ice-cloud water: Collection (cxc) efficiencies
!
!
      eiw(mgs) = 0.0
      if ( qx(mgs,li).gt.qxmin(li) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      
      
      if (xdia(mgs,lc,1).gt.ewi_dcmin .and. xdia(mgs,li,1).gt.ewi_dimin) then
! erm 5/10/2007 test following change:
!      if (xdia(mgs,lc,1).gt.12.0e-06 .and. xdia(mgs,li,1).gt.50.0e-06) then
      eiw(mgs) = eiw0
      end if
      if ( temg(mgs) .ge. 273.15 ) eiw(mgs) = 0.0
      end if

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
       cwrad = 0.5*xdia(mgs,lc,3)
       rwrad = 0.5*xdia(mgs,lr,3)
       
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
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      ers(mgs)=1.0
      end if
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,li).gt.qxmin(li) ) then
!        IF ( vtxbar(mgs,lr,1) .gt. vtxbar(mgs,li,1) .and.
!     :       xdia(mgs,lr,3) .gt. 200.e-6 .and. xdia(mgs,li,3) .gt. 100.e-6 ) THEN
         eri(mgs) = eri0
!      cwrad = 0.5*xdia(mgs,li,3)
!      eri(mgs) =
!     >  1.0*min((aradcw + cwrad*(bradcw + cwrad*
!     <  (cradcw + cwrad*(dradcw)))), 1.0)
!         ENDIF
!       if ( xdia(mgs,li,1) .lt. 10.e-6 ) eri(mgs)=0.0
       if ( xdia(mgs,li,3) .lt. eri_cimin ) eri(mgs)=0.0
      end if
!
!
!  Snow aggregates: Collection (cxc) efficiencies
!
! Modified by ERM with a linear function for small droplets and large
! snow agg. based numerical data from Wang and Ji (1992) in P&K 1997 (Fig. 14-13), which
! allows collection of very small droplets, albeit at low efficiency.  But slow
! fall speeds of snow make up for the efficiency.
!
      esw(mgs) = 0.0
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lc).gt.qxmin(lc) ) then
        esw(mgs) = 0.5
        if ( xdia(mgs,lc,1) .gt. 15.e-6 .and. xdia(mgs,ls,1) .gt. 100.e-6) then
          esw(mgs) = 0.5
        ELSEIF ( xdia(mgs,ls,1) .ge. 500.e-6 ) THEN
          esw(mgs) = Min(0.5, 0.05 + (0.8-0.05)/(40.e-6)*xdia(mgs,lc,1) )
        ENDIF
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lr).gt.qxmin(lr)  &
     &     .and. temg(mgs) .lt. tfr - 1.   &
     &                               ) then
      esr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,ls,1))
      IF ( qx(mgs,ls) < 1.e-4 .and. qx(mgs,lr) < 1.e-4 ) il2(mgs) = 1
      end if
      
      IF ( ipconc < 3 .and. temg(mgs) < tfr .and. qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lr) < 1.e-4 ) THEN
        il3(mgs) = 1
      ENDIF
!
!      if ( qx(mgs,ls).gt.qxmin(ls) ) then
      if ( temcg(mgs) < 0.0 ) then
            
      IF ( ipconc .lt. 4 .or. temcg(mgs) < esstem1 ) THEN
        ess(mgs) = 0.0
!        ess(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
!        ess(mgs)=min(0.1,ess(mgs))
      
      ELSE
      
        fac = Abs(ess0)
        IF ( iessopt == 2 ) THEN ! experimental code
!         IF ( wvel(mgs) > 2.0 .or. wvel(mgs) < -0.5 .or. ssi(mgs) < 1.0 ) THEN
         IF ( wvel(mgs) > 2.0 ) THEN
          ! assume convective cell or downdraft
           fac = 0.0
         ELSEIF ( wvel(mgs) > 1.0 ) THEN ! transition to stratiform range of values
           fac = Max(0.0, 2.0 - wvel(mgs))*fac
         ENDIF
        ELSEIF ( iessopt == 3 ) THEN ! factor based on ice supersat
           IF ( ssi(mgs) <= 1.0 ) THEN
             fac = 0.0
             ehsfac(mgs) = 0.0
           ELSEIF ( ssi(mgs) <= 1.02 ) THEN
             fac = fac*(ssi(mgs) - 1.0)/0.02
             ehsfac(mgs) = (ssi(mgs) - 1.0)/0.02
           ENDIF
        ELSEIF ( iessopt == 4 ) THEN ! factor based on ice supersat; very roughly based on Hosler et al. 1957 (J. Met.)
           IF ( ssi(mgs) <= 1.0 ) THEN
             fac = 0.1
             ehsfac(mgs) = 0.1
           ELSEIF ( ssi(mgs) <= 1.005 ) THEN ! ssi in range of 1.0 to 1.005
             fac = 0.1 + (ssi(mgs) - 1.0)*(fac - 0.1)/(1.005 - 1.0) !  Max(0.1, fac*(ssi(mgs) - 1.0)/0.005)
             ehsfac(mgs) = fac !  Max(0.1, (ssi(mgs) - 1.0)/0.005)
           ENDIF
        ELSEIF ( iessopt == 5 ) THEN ! factor based on ice supersat; very roughly based on Hosler et al. 1957 (J. Met.)
           IF ( ssi(mgs) < 0.90 ) THEN
             fac = 0.1
             ehsfac(mgs) = 0.1
           ELSEIF ( ssi(mgs) < 1.0 ) THEN ! ssi in range of 0.9 to 1.0
             fac = 0.1 + (ssi(mgs) - 0.9)*(fac - 0.1)/(1.0 - 0.9) 
             ehsfac(mgs) = fac ! Max(0.1, 0.1*(1.0 - ssi(mgs))/0.1)
           ENDIF
        ENDIF
        
        IF ( temcg(mgs) > esstem1 .and. temcg(mgs) < esstem2 ) THEN  ! only nonzero for T > esstem1
          ess(mgs) = fac*Exp(ess1*(esstem2) )*(temcg(mgs) - esstem1)/(esstem2 - esstem1) ! linear ramp up from zero at esstem1 to value at esstem2
        ELSEIF ( temcg(mgs) >= esstem2 ) THEN
          ess(mgs) = fac*Exp(ess1*Min( temcg(mgs), 0.0 ) )
        ENDIF
        
      ENDIF
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,li).gt.qxmin(li) ) then
       esiclsn(mgs) = esi_collsn
!      IF ( ipconc .lt. 4 ) THEN
      IF ( ipconc < 1 .and. lwsm6 ) THEN
        esi(mgs) = exp(0.7*min(temcg(mgs),0.0))
      ELSE
        esi(mgs) = esi0*exp(0.1*min(temcg(mgs),0.0))
        esi(mgs) = Min(0.1,esi(mgs))
      ENDIF
      IF ( ipconc .le. 3 ) THEN
       esi(mgs) =  exp(0.025*min(temcg(mgs),0.0)) ! LFO
!       esi(mgs) =  Min(0.5, exp(0.025*min(temcg(mgs),0.0)) ) ! LFO
!       esi(mgs)=0.5*exp(0.1*min(temcg(mgs),0.0))  ! 10ice
      ENDIF
!      ELSE ! zrnic/ziegler 1993
!      esi(mgs)= 0.1 ! 0.5*exp(0.1*min(temcg(mgs),0.0))
!      ENDIF
      if ( temg(mgs) .gt. 273.15 ) esi(mgs) = 0.0
      end if

!
!
!
!
!  Graupel: Collection (cxc) efficiencies
!
!
       xmascw(mgs) = xmas(mgs,lc)
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc).gt.qxmin(lc) ) then !{
       ehw(mgs) = 1.0
       IF ( iehw .eq. 0 ) THEN
       ehw(mgs) = ehw0  ! default value is 1.0
       ELSEIF ( iehw .eq. 1 .or. iehw .eq. 10 ) THEN
      cwrad = 0.5*xdia(mgs,lc,1)
      ehw(mgs) = Min( ehw0,    &
     &  ewfac*min((aradcw + cwrad*(bradcw + cwrad*   &
     &  (cradcw + cwrad*(dradcw)))), 1.0) )
      
       ELSEIF ( iehw .eq. 2 .or. iehw .eq. 10 ) THEN
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = igwr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lh,3)  ! changed to mean volume diameter
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)
 
!        write(iunit,*) 'slop1: ',slope1,slope2,ew(ic,ir),cwr(ic,2)

       x1 = ew(ic,  ir) + slope1*Max(0.0, (cwrad - cwr(ic,1)) )
       x2 = ew(icp1,ir) + slope2*Max(0.0, (cwrad - cwr(ic,1)) )
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       tmp = Max( 0.0, Min( 1.0, x1 + slope1*Max(0.0, (rwrad - grad(ir,1)) ) ) )
       ehw(mgs) = Min( ehw(mgs), tmp )

!       write(iunit,*) 'ehw: ',ehw(mgs),1.e6*cwrad,1.e6*rwrad,ic,ir,x1,x2
!       write(iunit,*)

!       ehw(mgs) = Max( 0.2, ehw(mgs) )
!  assume that ehw = 1 for zero air resistance (rho0 = 0.0) and extrapolate toward that
!      ehw(mgs) = ehw(mgs) + (ehw(mgs) - 1.0)*(rho0(mgs) - rho00)/rho00
!      ehw(mgs) = ehw(mgs) + (1.0 - ehw(mgs))*((Max(0.0,rho00 - rho0(mgs)))/rho00)**2

       ELSEIF ( iehw .eq. 3 .or. iehw .eq. 10 ) THEN ! use fraction of droplets greater than dmincw diameter
         tmp = Exp(- (dmincw/xdia(mgs,lc,1))**3)
         xmascw(mgs) = xmas(mgs,lc) + xdn0(lc)*(pi*dmincw**3/6.0) ! this is the average mass of the droplets with d > dmincw
         ehw(mgs) = Min( ehw(mgs), tmp )
       ELSEIF ( iehw .eq. 4 .or. iehw .eq. 10 ) THEN ! Cober and List 1993, eq. 19-20
         tmp =  &
     &   2.0*xdn(mgs,lc)*vtxbar(mgs,lh,1)*(0.5*xdia(mgs,lc,1))**2 &
     &  /(9.0*fadvisc(mgs)*0.5*xdia(mgs,lh,3))
         tmp = Max( 1.5, Min(10.0, tmp) )
         ehw(mgs) = Min( ehw(mgs), 0.55*Log10(2.51*tmp) )
       ENDIF
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehw(mgs)=0.0

       ehw(mgs) = Min( ehw0, ehw(mgs) )
       
       IF ( ibfc == -1 .and. temcg(mgs) < -41.0 ) THEN
        ehw(mgs) = 0.0
       ENDIF 

      end if !}
!
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lr).gt.qxmin(lr)    &
!     &     .and. temg(mgs) .lt. tfr    &
     &                               ) then
!      ehr(mgs) = Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lh,1))
!      ehr(mgs) = 1.0
       ehr(mgs) = Exp(-(40.e-6)/xdia(mgs,lr,3))*Exp(-40.e-6/xdia(mgs,lh,3))
       ehr(mgs) = Min( ehr0, ehr(mgs) )
      end if
!
      IF ( qx(mgs,ls).gt.qxmin(ls) ) THEN
        IF ( ipconc .ge. 4 ) THEN
        ehscnv(mgs) = ehs0*exp(ehs1*min(temcg(mgs),0.0)) ! for 2-moment, used as default for ehs and ehls. Otherwise not used for snow->graupel conversion
        ELSE
        ehscnv(mgs) = exp(0.09*min(temcg(mgs),0.0))
        ENDIF
        
        IF ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc) >= qxmin(lc)  ) THEN
!          ehsclsn(mgs) = ehs_collsn
!          ehs(mgs) = ehscnv(mgs)*ehsfac(mgs)*Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  )
!        ELSEIF ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc) >= qxmin(lc)  ) then
          ehsclsn(mgs) = ehs_collsn
          IF ( xdia(mgs,ls,3) < 40.e-6 ) THEN
            ehsclsn(mgs) = 0.0
          ELSEIF ( xdia(mgs,ls,3) < 150.e-6 ) THEN
            ehsclsn(mgs) =  ehs_collsn*(xdia(mgs,ls,3) - 40.e-6)/(150.e-6 - 40.e-6)
          ELSE
            ehsclsn(mgs) = ehs_collsn
          ENDIF
!          ehs(mgs) = ehscnv(mgs)*Min(1.0, Max(0., xdn(mgs,lh) - xdnmn(lh)*1.2)/xdnmn(lh)  ) ! shut off qhacs as graupel goes to lowest density
          ehs(mgs) = ehscnv(mgs)*Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  ) ! shut off qhacs as graupel goes to low density; limits scavenging of snow in bright band
!          ehs(mgs) = ehscnv(mgs) ! *Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  ) ! shut off qhacs as graupel goes to low density
          ehs(mgs) = Min(ehs(mgs),ehsmax)
        end if
      ENDIF
!
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehiclsn(mgs) = ehi_collsn
      ehi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehi(mgs) = Min( ehimax, Max( ehi(mgs), ehimin ) )
!      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehi(mgs) = 0.0
      end if

      IF ( lis > 1 ) THEN
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lis).gt.qxmin(lis) ) then
      ehisclsn(mgs) = ehi_collsn
      ehis(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehis(mgs) = Min( ehimax, Max( ehis(mgs), ehimin ) )
!      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehis(mgs) = 0.0
      end if
      ENDIF

!
!
!  Frozen drops: Collection (cxc) efficiencies
!
!
      IF ( lf .gt. 1 ) THEN

      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lc).gt.qxmin(lc) ) then
       IF ( iefw == 3 ) iefw = 3
       IF ( iefw == 4 ) iefw = 4
       efw(mgs) = efw0
       IF ( iefw .eq. 0 ) THEN
       efw(mgs) = efw0  ! default value is 1.0
       ELSEIF ( iefw .eq. 1 .or. iefw .eq. 10 ) THEN
      cwrad = 0.5*xdia(mgs,lc,1)
      efw(mgs) = Min( efw0,    &
     &  ewfac*min((aradcw + cwrad*(bradcw + cwrad*   &
     &  (cradcw + cwrad*(dradcw)))), 1.0) )
      
       ELSEIF ( iefw .eq. 2 .or. iefw .eq. 10 ) THEN
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = ifwr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lf,3)  ! changed to mean volume diameter
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)
       
       x1 = ew(ic,  ir) + slope1*(cwrad - cwr(ic,1))
       x2 = ew(icp1,ir) + slope2*(cwrad - cwr(ic,1))
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       tmp = Max( 0.0, Min( 1.0, x1 + slope1*(rwrad - grad(ir,1)) ) )
         efw(mgs) = Min( efw(mgs), tmp )
       efw(mgs) = Min( efw0, efw(mgs) )
!       ehw(mgs) = Max( 0.2, ehw(mgs) )
!  assume that ehw = 1 for zero air resistance (rho0 = 0.0) and extrapolate toward that
!      ehw(mgs) = ehw(mgs) + (ehw(mgs) - 1.0)*(rho0(mgs) - rho00)/rho00
!      efw(mgs) = efw(mgs) + (1.0 - efw(mgs))*((Max(0.0,rho00 - rho0(mgs)))/rho00)**2

       ELSEIF ( iefw .eq. 3 .or. iefw .eq. 10 ) THEN ! use fraction of droplets greater than 15 micron diameter
         tmp = Exp(- (dmincw/xdia(mgs,lc,1))**3)
         efw(mgs) = Min( efw(mgs), tmp )
       ELSEIF ( iefw .eq. 4 .or. iefw .eq. 10 ) THEN ! Cober and List 1993
         tmp =  &
     &   2.0*xdn(mgs,lc)*vtxbar(mgs,lf,1)*(0.5*xdia(mgs,lc,1))**2 &
     &  /(9.0*fadvisc(mgs)*0.5*xdia(mgs,lf,3))
         tmp = Max( 1.5, Min(10.0, tmp) )
         efw(mgs) = Min( efw(mgs), 0.55*Log10(2.51*tmp) )
       ENDIF
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) efw(mgs)=0.0
       efw(mgs) = Min( efw0, efw(mgs) )

       IF ( ibfc == -1 .and. temcg(mgs) < -41.0 ) THEN 
        efw(mgs) = 0.0
       ENDIF 

      end if
!
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lr).gt.qxmin(lr)    &
!     &     .and. temg(mgs) .lt. tfr    &
     &                               ) then
        efr(mgs) = 1.0
        efr(mgs) = Min( efr0, efr(mgs) )
      else
        efr(mgs) = 0.0
      end if
!
      IF ( qx(mgs,ls).gt.qxmin(ls) ) THEN
        if ( qx(mgs,lf).gt.qxmin(lf)  ) then
          efsclsn(mgs) = efs_collsn
          IF ( xdia(mgs,ls,3) < 40.e-6 ) THEN
            efsclsn(mgs) = 0.0
          ELSEIF ( xdia(mgs,ls,3) < 150.e-6 ) THEN
            efsclsn(mgs) =  efs_collsn*(xdia(mgs,ls,3) - 40.e-6)/(150.e-6 - 40.e-6)
          ELSE
            efsclsn(mgs) = efs_collsn
          ENDIF
          IF ( qx(mgs,lc) > qxmin(lc) ) THEN
          efs(mgs) = ehscnv(mgs)*Min(1.0, Max(0.0,xdn(mgs,lf) - 300.)/300.  ) ! shut off qhacs as graupel goes to low density; limits scavenging of snow in bright band
          efs(mgs) = Min(efs(mgs),ehsmax)
!          efs(mgs) = ehscnv(mgs)
!          efs(mgs) = Min(efs(mgs),ehsmax)
          ENDIF
        end if
      ENDIF
!
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,li).gt.qxmin(li) ) then
      eficlsn(mgs) = efi_collsn
      efi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      efi(mgs) = Min( ehimax, Max( efi(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) efi(mgs) = 0.0
      end if

      IF ( lis > 1 ) THEN
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lis).gt.qxmin(lis) ) then
      efisclsn(mgs) = efi_collsn
      efis(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      efis(mgs) = Min( ehimax, Max( efis(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) efis(mgs) = 0.0
      end if
      ENDIF


      ENDIF ! lf .gt. 1
 
!
!
!  Hail: Collection (cxc) efficiencies
!
!
      IF ( lhl .gt. 1 ) THEN

      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lc).gt.qxmin(lc) ) then
       IF ( iehw == 3 ) iehlw = 3
       IF ( iehw == 4 ) iehlw = 4
       ehlw(mgs) = ehlw0
       IF ( iehlw .eq. 0 ) THEN
       ehlw(mgs) = ehlw0  ! default value is 1.0
       ELSEIF ( iehlw .eq. 1 .or. iehlw .eq. 10 ) THEN
      cwrad = 0.5*xdia(mgs,lc,1)
      ehlw(mgs) = Min( ehlw0,    &
     &  ewfac*min((aradcw + cwrad*(bradcw + cwrad*   &
     &  (cradcw + cwrad*(dradcw)))), 1.0) )
      
       ELSEIF ( iehlw .eq. 2 .or. iehlw .eq. 10 ) THEN
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = ihlr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lhl,3)  ! changed to mean volume diameter
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)
       
       x1 = ew(ic,  ir) + slope1*(cwrad - cwr(ic,1))
       x2 = ew(icp1,ir) + slope2*(cwrad - cwr(ic,1))
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       tmp = Max( 0.0, Min( 1.0, x1 + slope1*(rwrad - grad(ir,1)) ) )
         ehlw(mgs) = Min( ehlw(mgs), tmp )
       ehlw(mgs) = Min( ehlw0, ehlw(mgs) )
!       ehw(mgs) = Max( 0.2, ehw(mgs) )
!  assume that ehw = 1 for zero air resistance (rho0 = 0.0) and extrapolate toward that
!      ehw(mgs) = ehw(mgs) + (ehw(mgs) - 1.0)*(rho0(mgs) - rho00)/rho00
!      ehlw(mgs) = ehlw(mgs) + (1.0 - ehlw(mgs))*((Max(0.0,rho00 - rho0(mgs)))/rho00)**2

       ELSEIF ( iehlw .eq. 3 .or. iehlw .eq. 10 ) THEN ! use fraction of droplets greater than 15 micron diameter
         tmp = Exp(- (dmincw/xdia(mgs,lc,1))**3)
         ehlw(mgs) = Min( ehlw(mgs), tmp )
       ELSEIF ( iehlw .eq. 4 .or. iehlw .eq. 10 ) THEN ! Cober and List 1993
         tmp =  &
     &   2.0*xdn(mgs,lc)*vtxbar(mgs,lhl,1)*(0.5*xdia(mgs,lc,1))**2 &
     &  /(9.0*fadvisc(mgs)*0.5*xdia(mgs,lhl,3))
         tmp = Max( 1.5, Min(10.0, tmp) )
         ehlw(mgs) = Min( ehlw(mgs), 0.55*Log10(2.51*tmp) )
       ENDIF
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehlw(mgs)=0.0
       ehlw(mgs) = Min( ehlw0, ehlw(mgs) )

       IF ( ibfc == -1 .and. temcg(mgs) < -41.0 ) THEN 
        ehlw(mgs) = 0.0
       ENDIF 

      end if
!
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lr).gt.qxmin(lr)    &
!     &     .and. temg(mgs) .lt. tfr    &
     &                               ) then
        ehlr(mgs) = 1.0
       ehlr(mgs) = Min( ehlr0, ehlr(mgs) )
      end if
!
      IF ( qx(mgs,ls).gt.qxmin(ls) ) THEN
        if ( qx(mgs,lhl).gt.qxmin(lhl)  ) then
          ehlsclsn(mgs) = ehls_collsn
          ehls(mgs) = ehscnv(mgs)
          ehls(mgs) = Min(ehls(mgs),ehsmax)
        end if
      ENDIF
!
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehliclsn(mgs) = ehli_collsn
      ehli(mgs)=eii0hl*exp(eii1hl*min(temcg(mgs),0.0))
      ehli(mgs) = Min( ehimax, Max( ehli(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehli(mgs) = 0.0
      end if

      IF ( lis > 1 ) THEN
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lis).gt.qxmin(lis) ) then
      ehlisclsn(mgs) = ehli_collsn
      ehlis(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehlis(mgs) = Min( ehimax, Max( ehlis(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehlis(mgs) = 0.0
      end if
      ENDIF


      ENDIF ! lhl .gt. 1

      ENDDO  ! mgs loop for collection efficiencies

!
!
!
!  Set flags for plates vs. columns
!
!
      do mgs = 1,ngscnt
!
      xplate(mgs) = 0.0
      xcolmn(mgs) = 1.0
!
!      if ( temcg(mgs) .lt. 0. .and. temcg(mgs) .ge. -4. ) then
!      xplate(mgs) = 1.0
!      xcolmn(mgs) = 0.0
!      end if
!c
!      if ( temcg(mgs) .lt. -4. .and. temcg(mgs) .ge. -9. ) then
!      xplate(mgs) = 0.0
!      xcolmn(mgs) = 1.0
!      end if
!c
!      if ( temcg(mgs) .lt. -9. .and. temcg(mgs) .ge. -22.5 ) then
!      xplate(mgs) = 1.0
!      xcolmn(mgs) = 0.0
!      end if
!c
!      if ( temcg(mgs) .lt. -22.5 .and. temcg(mgs) .ge. -90. ) then
!      xplate(mgs) = 0.0
!      xcolmn(mgs) = 1.0
!      end if
!
      end do
      
      

!
!  Frozen drops
!
      IF ( lf <= 1 ) THEN
      qfdpv(1:ngscnt) = 0.0
      qfacr(1:ngscnt) = 0.0
      qfacw(1:ngscnt) = 0.0
      qfacs(1:ngscnt) = 0.0
      qfaci(1:ngscnt) = 0.0
      qfmlr(1:ngscnt) = 0.0
      qfmul1(1:ngscnt) = 0.0

      cfmul1(1:ngscnt) = 0.0
      cfmlr(1:ngscnt) = 0.0
      cfmlrr(1:ngscnt) = 0.0
      cfsbv(1:ngscnt) = 0.0
      cfcev(1:ngscnt) = 0.0
      ENDIF
            

!
!
!
!  Collection growth equations....
!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: rain collects xxxxx'
!
      do mgs = 1,ngscnt
      qracw(mgs) =  0.0
      IF ( qx(mgs,lr) .gt. qxmin(lr) .and. erw(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .lt. 3 ) THEN
       IF ( erw(mgs) .gt. 0.0 .and. qx(mgs,lr) .gt. 1.e-7 ) THEN
       vt = (ar*(xdia(mgs,lc,1)**br))*rhovt(mgs)
       qracw(mgs) =    &
     &   (0.25)*pi*erw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lr) &
!     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *Max(0.0, vtxbar(mgs,lr,1)-vt)   &
     &  *(  gf3*xdia(mgs,lr,2)    &
     &    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )
!       qracw(mgs) = 0.0
!      write(iunit,*) 'qracw,cx =',qracw(mgs),1.e6*xdia(mgs,lr,1),erw(mgs)
!      write(iunit,*) 'qracw,cx =',qracw(mgs),cx(mgs,lc),kgs(mgs),cx(mgs,lr),1.e6*xdia(mgs,lr,1),vtxbar(mgs,lr,1),vt
!      write(iunit,*) 'vtr: ',vtxbar(mgs,lr,1), ar*gf4br/6.0*xdia(mgs,lr,1)**br, rhovt(mgs),
!     :         ar*gf4br/6.0*xdia(mgs,lr,1)**br * rhovt(mgs)
       ENDIF
      ELSE

      IF ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) THEN 
       rwrad = 0.5*xdia(mgs,lr,3)
        IF ( rwrad .gt. rh(mgs) ) THEN ! .or. cx(mgs,lr) .gt. nh(mgs) ) THEN
         IF ( rwrad .gt. rwradmn ) THEN
!      DM1CCC=A2*XNC*XNR*XVC*(((CNU+2.)/(CNU+1.))*XVC+XVR)       ! (A12)
!     NOTE: Result is independent of imurain, assumes mucloud = 3
           qracw(mgs) = erw(mgs)*aa2*cx(mgs,lr)*cx(mgs,lc)*xmas(mgs,lc)*   &
     &        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,lr))/rho0(mgs) !*rhoinv(mgs)
         ELSE

          IF ( imurain == 3 ) THEN

!      DM1CCC=A1*XNC*XNR*(((CNU+3.)*(CNU+2.)/(CNU+1.)**2)*XVC**3+ ! (A14)
!     1 ((RNU+2.)/(RNU+1.))*XVC*XVR**2)

!           qracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*xdn(mgs,lc)*   &
!     &        ((cnu + 3.)*(cnu + 2.)*xv(mgs,lc)**3/(cnu + 1.)**2 +    &
!     &         (alpha(mgs,lr) + 2.)*xv(mgs,lc)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))/rho0(mgs) !*rhoinv(mgs)
! save multiplies by converting cx*xdn*xv/rho0 to qx
           qracw(mgs) = aa1*cx(mgs,lr)*(qx(mgs,lc)-qcwresv(mgs))*   &
     &        ((alpha(mgs,lc) + 3.)*(alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.)**2 +    &
     &         (alpha(mgs,lr) + 2.)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.)) 
           
           ELSE ! imurain == 1

           qracw(mgs) = aa1*cx(mgs,lr)*(qx(mgs,lc)-qcwresv(mgs))*   &
     &        ((alpha(mgs,lc) + 3.)*(alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.)**2 +    &
     &         (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)*xv(mgs,lr)**2/ &
     &          ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.))) 
           
           ENDIF
           
         ENDIF
        ENDIF
        ENDIF
       ENDIF
!       qracw(mgs) = Min(qracw(mgs), qx(mgs,lc))
       qracw(mgs) = Min(qracw(mgs), qcmxd(mgs))
       ENDIF
      end do
!
      do mgs = 1,ngscnt
      qraci(mgs) = 0.0
      qracif(mgs) = 0.0
      craci(mgs) = 0.0
      qracs(mgs) = 0.0
      IF ( eri(mgs) .gt. 0.0 .and. iacr .ge. 1 .and. xdia(mgs,lr,3) .gt. 2.*rwradmn ) THEN
        IF ( ipconc .ge. 3 ) THEN

           tmp = eri(mgs)*aa2*cx(mgs,lr)*cx(mgs,li)*   &
     &        ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,lr))

        qraci(mgs) = Min( qxmxd(mgs,li), tmp*xmas(mgs,li)*rhoinv(mgs) )
        craci(mgs) = Min( cxmxd(mgs,li), tmp )

!       vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +
!     :            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )
!
!          qraci(mgs) = 0.25*pi*eri(mgs)*cx(mgs,lr)*qx(mgs,li)*vt*
!     :         (  da0(lr)*xdia(mgs,lr,3)**2 +
!     :            dab1(lr,li)*xdia(mgs,lr,3)*xdia(mgs,li,3) +
!     :            da1(li)*xdia(mgs,li,3)**2 )
!
!
!       vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +
!     :            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )
!
!          craci(mgs) = 0.25*pi*eri(mgs)*cx(mgs,lr)*cx(mgs,li)*vt*
!     :         (  da0(lr)*xdia(mgs,lr,3)**2 +
!     :            dab0(lr,li)*xdia(mgs,lr,3)*xdia(mgs,li,3) +
!     :            da0(li)*xdia(mgs,li,3)**2 )
!
!          qraci(mgs) = Min( qraci(mgs), qxmxd(mgs,li) )
!          craci(mgs) = Min( craci(mgs), cxmxd(mgs,li) )

        ELSE
          qraci(mgs) =    &
     &     min(   &
     &     (0.25)*pi*eri(mgs)*qx(mgs,li)*cx(mgs,lr)   &
     &    *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))   &
     &    *(  gf3*xdia(mgs,lr,2)    &
     &      + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1)    &
     &      + gf1*xdia(mgs,li,2) )     &
     &    , qimxd(mgs))
        ENDIF
      if ( temg(mgs) .gt. 268.15 ) then
      qraci(mgs) = 0.0
      end if
      ENDIF
      end do
!
      IF ( ipconc < 3 ) THEN
      do mgs = 1,ngscnt
      qracs(mgs) = 0.0
      IF ( ers(mgs) .gt. 0.0 .and. ipconc < 3 ) THEN
       IF ( lwsm6 .and. ipconc == 0 ) THEN
         vt = vt2ave(mgs)
       ELSE
         vt = vtxbar(mgs,ls,1)
       ENDIF
      qracs(mgs) =      &
     &   min(     &
     &   ((0.25)*pi/gf4)*ers(mgs)*qx(mgs,ls)*cx(mgs,lr)     &
     &  *abs(vtxbar(mgs,lr,1)-vt)     &
     &  *(  gf6*gf1*xdia(mgs,ls,2)     &
     &    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lr,1)      &
     &    + gf4*gf3*xdia(mgs,lr,2) )      &
     &  , qsmxd(mgs))
      ENDIF
      end do
      ENDIF

!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: snow collects xxxxx'
!
      do mgs = 1,ngscnt
      qsacw(mgs) =  0.0
      csacw(mgs) =  0.0
      vsacw(mgs) =  0.0
      IF ( esw(mgs) .gt. 0.0 ) THEN

       IF ( ipconc .ge. 4 ) THEN
!      QSACC=CECS*RVT*A2*XNC*XNS*XVC*ROS*
!     *    (((CNU+2.)/(CNU+1.))*XVC+XVS)/RO

!        tmp = esw(mgs)*rvt*aa2*cx(mgs,ls)*cx(mgs,lc)*
!     :        ((cnu + 2.)*xv(mgs,lc)/(cnu + 1.) + xv(mgs,ls))
        tmp = 1.0*rvt*aa2*cx(mgs,ls)*cx(mgs,lc)*   &
     &        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,ls))

        qsacw(mgs) = Min( qxmxd(mgs,lc), tmp*xmas(mgs,lc)*rhoinv(mgs) )
        csacw(mgs) = Min( cxmxd(mgs,lc), tmp )

          IF ( lvol(ls) .gt. 1 ) THEN
             IF ( temg(mgs) .lt. 273.15) THEN
             rimdn(mgs,ls) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,ls,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
             rimdn(mgs,ls) = Min( Max( rimc3, rimdn(mgs,ls) ), rimc4 )
             ELSE
             rimdn(mgs,ls) = 1000.
             ENDIF

           vsacw(mgs) = rho0(mgs)*qsacw(mgs)/rimdn(mgs,ls)

          ENDIF


!        qsacw(mgs) = cecs*aa2*cx(mgs,ls)*cx(mgs,lc)*xmas(mgs,lc)*
!     :        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,ls))*rhoinv(mgs)
       ELSE
!      qsacw(mgs) =
!     >   min(
!     >   ((0.25)*pi)*esw(mgs)*qx(mgs,lc)*cx(mgs,ls)
!     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))
!     >  *(  gf3*xdia(mgs,ls,2)
!     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,lc,1)
!     >    + gf1*xdia(mgs,lc,2) )
!     <  , qcmxd(mgs))

            vt = abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))

          qsacw(mgs) = 0.25*pi*esw(mgs)*cx(mgs,ls)*qx(mgs,lc)*vt*   &
     &         (  da0(ls)*xdia(mgs,ls,3)**2 +     &
     &            dab1(ls,lc)*xdia(mgs,ls,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )
        qsacw(mgs) = Min( qsacw(mgs), qxmxd(mgs,ls) )
        csacw(mgs) = rho0(mgs)*qsacw(mgs)/xmas(mgs,lc)
       ENDIF
      ENDIF
      end do
!
!
      do mgs = 1,ngscnt
      qsaci(mgs) = 0.0
      csaci(mgs) = 0.0
      csaci0(mgs) = 0.0
      IF ( ipconc .ge. 4 ) THEN
      IF ( esi(mgs) .gt. 0.0 .or. ( ipelec > 0 .and. esiclsn(mgs) > 0.0 )) THEN
!      QSCOI=CEXS*RVT*A2*XNCI*XNS*XVCI*ROS*
!     *  (((CINU+2.)/(CINU+1.))*VCIP+XVS)/RO

        tmp = esiclsn(mgs)*rvt*aa2*cx(mgs,ls)*cx(mgs,li)*   &
     &        ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,ls))

        qsaci(mgs) = Min( qxmxd(mgs,li), esi(mgs)*tmp*xmas(mgs,li)*rhoinv(mgs) )
        csaci0(mgs) = tmp
        csaci(mgs) = Min(cxmxd(mgs,li), esi(mgs)*tmp )

!      qsaci(mgs) =
!     >   min(
!     >   ((0.25)*pi)*esi(mgs)*qx(mgs,li)*cx(mgs,ls)
!     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))
!     >  *(  gf3*xdia(mgs,ls,2)
!     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,li,1)
!     >    + gf1*xdia(mgs,li,2) )
!     <  , qimxd(mgs))
      ENDIF
      ELSE ! 
      IF ( esi(mgs) .gt. 0.0 ) THEN
         qsaci(mgs) =    &
     &   min(   &
     &   ((0.25)*pi)*esi(mgs)*qx(mgs,li)*cx(mgs,ls)   &
     &  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))   &
     &  *(  gf3*xdia(mgs,ls,2)    &
     &    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,li,1)    &
     &    + gf1*xdia(mgs,li,2) )     &
     &  , qimxd(mgs))
      ENDIF
      ENDIF
      end do
!
!
!
      do mgs = 1,ngscnt
      qsacr(mgs) = 0.0
      qsacrs(mgs) = 0.0
      csacr(mgs) = 0.0
      IF ( esr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
!       vt = Sqrt((vtxbar(mgs,ls,1)-vtxbar(mgs,lr,1))**2 + 
!     :            0.04*vtxbar(mgs,ls,1)*vtxbar(mgs,lr,1) )
!       qsacr(mgs) = esr(mgs)*cx(mgs,ls)*vt*
!     :     qx(mgs,lr)*0.25*pi*
!     :      (3.02787*xdia(mgs,lr,2) + 
!     :       3.30669*xdia(mgs,ls,1)*xdia(mgs,lr,1) + 
!     :       2.*xdia(mgs,ls,2))
!        qsacr(mgs) = Min( qsacr(mgs), qrmxd(mgs) )
!        csacr(mgs) = qsacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!        csacr(mgs) = min(csacr(mgs),crmxd(mgs))
      ELSE
       IF ( lwsm6 .and. ipconc == 0 ) THEN
         vt = vt2ave(mgs)
       ELSE
         vt = vtxbar(mgs,ls,1)
       ENDIF
       
       qsacr(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*esr(mgs)*qx(mgs,lr)*cx(mgs,ls)   &
     &  *abs(vtxbar(mgs,lr,1)-vt)   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,ls,1)    &
     &    + gf4*gf3*xdia(mgs,ls,2) )    &
     &  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do
!
!
!

      IF ( lf > 1 ) THEN !{
      if (ndebug .gt. 0 ) write(0,*) 'Collection: frozen drops collects xxxxx'

      DO mgs = 1,ngscnt
      qfacw(mgs) = 0.0
      qfacwmlr(mgs) = 0.0
      rarx(mgs,lf) = 0.0
      vfacw(mgs) = 0.0
      vfsoak(mgs) = 0.0
      zfacw(mgs) = 0.0
      
      
!      ENDDO
      
      
      IF ( efw(mgs) .gt. 0.0 ) THEN

        IF ( ipconc .ge. 2 ) THEN !{{

         vt = abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1)) 

          qfacw(mgs) = 0.25*pi*efw(mgs)*cx(mgs,lf)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab1lh(mgs,lf,lc)*xdia(mgs,lf,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 ) 

          qfacw(mgs) = Min( qfacw(mgs), 0.5*qx(mgs,lc)*dtpinv )
        
         IF ( lzf .gt. 1 ) THEN
          tmp = qx(mgs,lf)/cx(mgs,lf)
         ENDIF
        
        ELSE !}
         qfacw(mgs) =    &
     &   min(   &
     &   ((0.25)*pi)*efw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lf)   &
     &  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf3*xdia(mgs,lf,2)    &
     &    + 2.0*gf2*xdia(mgs,lf,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )     &
     &    , 0.5*(qx(mgs,lc)-qcwresv(mgs))*dtpinv)
       
         
       ENDIF !}

          qfacwmlr(mgs) = qfacw(mgs)
          IF ( temg(mgs) > tfr .and. iqhacwshr == 0 ) THEN
            qfacw(mgs) = 0.0
          ENDIF

          IF ( lvol(lf) .gt. 1 .or. lhl .gt. 1 ) THEN ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             IF ( temg(mgs) .lt. 273.15) THEN
               IF ( irimdenopt == 1 ) THEN ! Heymsfield and Pflaum (1985)
               vt = ( (1.0-rimdenvwgt)*vtxbar(mgs,lf,1) + rimdenvwgt*vtxbar(mgs,lf,2) )
               
             rimdn(mgs,lf) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vt )   &
     &                /(temg(mgs)-273.15))**(rimc2)
!             rimdn(mgs,lf) = Min( Max( hdnmn, rimc3, rimdn(mgs,lf) ), rimc4 )
             rimdn(mgs,lf) = Min( Max( rimc3, rimdn(mgs,lf) ), rimc4 )


               ELSEIF ( irimdenopt == 2 ) THEN ! Cober and List (1993)

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lf,1) + rimdenvwgt*vtxbar(mgs,lf,2) )   &
     &                /(temg(mgs)-273.15))
                tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) ) ! have to limit range of "R" because quadratic function starts to decrease (unphysically) at higher values
                
                rimdn(mgs,lf) = 1000.*(0.051 + 0.114*tmp - 0.0055*tmp**2)

               ELSEIF ( irimdenopt == 3 .or. irimdenopt == 4) THEN ! Macklin (3) or Saunders and Hosseini 2001

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lf,1) + rimdenvwgt*vtxbar(mgs,lf,2) )   &
     &                /(temg(mgs)-273.15))
              !  tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                IF ( irimdenopt == 3 ) THEN
                  rimdn(mgs,lf) =  Min(900., Max( 170., 110.*tmp**0.76 ) )
                ELSEIF ( irimdenopt == 4 ) THEN ! Saunders and Hosseini
                  rimdn(mgs,lf) =  Min(917., Max( 10.,  900.0*(1.0 - 0.905**tmp ) ) )
                ENDIF
               
               ENDIF
             ELSE
             rimdn(mgs,lf) = 1000.
             ENDIF
             
             IF ( lvol(lf) > 1 ) vfacw(mgs) = rho0(mgs)*qfacw(mgs)/rimdn(mgs,lf)

          ENDIF
      
        IF ( qx(mgs,lf) .gt. qxmin(lf) .and. ipelec .ge. 1 ) THEN
         rarx(mgs,lf) =     &
     &    qfacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lf,2)*cx(mgs,lf))
        ENDIF
      
      ENDIF  
      end do   
      
!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: frozen drops collects ice/snow'

      do mgs = 1,ngscnt
      qfaci(mgs) = 0.0
      qfaci0(mgs) = 0.0
      IF ( efi(mgs) .gt. 0.0 ) THEN
       IF (  ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,li,1) )

          qfaci0(mgs) = 0.25*pi*eficlsn(mgs)*cx(mgs,lf)*qx(mgs,li)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab1lh(mgs,lf,li)*xdia(mgs,lf,3)*xdia(mgs,li,3) +    &
     &            da1(li)*xdia(mgs,li,3)**2 ) 
          qfaci(mgs) = Min( efi(mgs)*qfaci0(mgs), qimxd(mgs) )
       ENDIF
      ENDIF
      end do   

!
!
      do mgs = 1,ngscnt
      qfacs(mgs) = 0.0
      qfacs0(mgs) = 0.0
      IF ( efs(mgs) .gt. 0.0 ) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,ls,1) )

          qfacs0(mgs) = 0.25*pi*efsclsn(mgs)*cx(mgs,lf)*qx(mgs,ls)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab1lh(mgs,lf,ls)*xdia(mgs,lf,3)*xdia(mgs,ls,3) +    &
     &            da1(ls)*xdia(mgs,ls,3)**2 ) 
      
          qfacs(mgs) = Min( efs(mgs)*qfacs0(mgs), qsmxd(mgs) )

        ENDIF
      ENDIF
      end do   
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: frozen drops collects rain'

      do mgs = 1,ngscnt
      qfacr(mgs) = 0.0
      qfacrmlr(mgs) = 0.0
      vfacr(mgs) = 0.0
      cfacr(mgs) = 0.0
      zfacr(mgs) = 0.0
      IF ( temg(mgs) .gt. tfr ) raindn(mgs,lf) = 1000.0

      IF ( efr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,lr,1) )
     
       qfacr(mgs) = 0.25*pi*efr(mgs)*cx(mgs,lf)*qx(mgs,lr)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab1lh(mgs,lf,lr)*xdia(mgs,lf,3)*xdia(mgs,lr,3) +    &
     &            da1lr(mgs)*xdia(mgs,lr,3)**2 )
!     &            da1(lr)*xdia(mgs,lr,3)**2 )

        qfacr(mgs) = Min( qfacr(mgs), qxmxd(mgs,lr) )

            qfacrmlr(mgs) = qfacr(mgs)
        
        IF ( temg(mgs) > tfr .and. iehr0c == 0 ) THEN
          qfacr(mgs) = 0.0

          IF ( iqhacrmlr == 0 ) THEN
              qfacrmlr(mgs) = -qfacw(mgs)
          ENDIF

        ELSE

!          cfacr(mgs) = qfacr(mgs)*cx(mgs,lr)/qx(mgs,lr)

        cfacr(mgs) = 0.25*pi*efr(mgs)*cx(mgs,lf)*cx(mgs,lr)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab0lh(mgs,lf,lr)*xdia(mgs,lf,3)*xdia(mgs,lr,3) +    &
     &            da0lr(mgs)*xdia(mgs,lr,3)**2 )

          cfacr(mgs) = min(cfacr(mgs),crmxd(mgs))

!          IF ( lzf .gt. 1  .and. cx(mgs,lf) > cxmin ) THEN
!            tmp = qx(mgs,lf)/cx(mgs,lf)
!          ENDIF

         ENDIF ! temg > tfr
            
       ENDIF
          
        IF ( lvol(lf) .gt. 1 .or. lhl .gt. 1 ) THEN ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             IF ( temg(mgs) .lt. 273.15) THEN
             raindn(mgs,lf) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lr,3))   &
     &                *((0.60)*vt)   &
     &                /(temg(mgs)-273.15))**(rimc2)

             raindn(mgs,lf) = Min( Max( rimc3, rimdn(mgs,lf) ), rimc4 )
             ELSE
             raindn(mgs,lf) = 1000.
             ENDIF
             
             IF ( raindn(mgs,lf) > 0.0 )  vfacr(mgs) = rho0(mgs)*qfacr(mgs)/raindn(mgs,lf)
        ENDIF
      ENDIF
      end do

      if (ndebug .gt. 0 ) write(0,*) 'Collection: frozen drops collects rain -- done'

      ELSE !}


      DO mgs = 1,ngscnt
      qfacw(mgs) = 0.0
      vfacw(mgs) = 0.0
      vfsoak(mgs) = 0.0
      zfacw(mgs) = 0.0
      qfacr(mgs) = 0.0
      qfacrmlr(mgs) = 0.0
      vfacr(mgs) = 0.0
      cfacr(mgs) = 0.0
      zfacr(mgs) = 0.0
      qfacs(mgs) = 0.0
      qfacs0(mgs) = 0.0
      qfacis(mgs) = 0.0
      qfacis0(mgs) = 0.0
      qfaci(mgs) = 0.0
      qfaci0(mgs) = 0.0
      ENDDO
      
      ENDIF ! lf > 1


      if (ndebug .gt. 0 ) write(0,*) 'Collection: graupel collects xxxxx'
!
      do mgs = 1,ngscnt
      qhacw(mgs) = 0.0
      qhacwmlr(mgs) = 0.0
      rarx(mgs,lh) = 0.0
      vhacw(mgs) = 0.0
      vhsoak(mgs) = 0.0
      zhacw(mgs) = 0.0
      
      IF ( .false. ) THEN
        vtmax = 1./(gz(kgs(mgs))*dtp)
        vtxbar(mgs,lh,1) = Min( vtmax, vtxbar(mgs,lh,1))
        vtxbar(mgs,lh,2) = Min( vtmax, vtxbar(mgs,lh,2))
        vtxbar(mgs,lh,3) = Min( vtmax, vtxbar(mgs,lh,3))
      ENDIF
      IF ( ehw(mgs) .gt. 0.0 ) THEN

        IF ( ipconc .ge. 2 ) THEN

        IF ( .false. ) THEN  
        qhacw(mgs) = (ehw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lh)*pi*   &
     &    abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*   &
     &    (2.0*xdia(mgs,lh,1)*(xdia(mgs,lh,1) +    &
     &         xdia(mgs,lc,1)*gf73rds) +    &
     &      xdia(mgs,lc,2)*gf83rds))/4.     
     
         ELSE  ! using Seifert coefficients
            vt = abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1)) 

          qhacw(mgs) = 0.25*pi*ehw(mgs)*cx(mgs,lh)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,lc)*xdia(mgs,lh,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 ) 
         
         ENDIF
          qhacw(mgs) = Min( qhacw(mgs), 0.5*qx(mgs,lc)*dtpinv )
        
         IF ( lzh .gt. 1 ) THEN
          tmp = qx(mgs,lh)/cx(mgs,lh)
          
!!          g1 = (6.0 + alpha(mgs,lh))*(5.0 + alpha(mgs,lh))*(4.0 + alpha(mgs,lh))/
!!     :         ((3.0 + alpha(mgs,lh))*(2.0 + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))
!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
         ENDIF
        
        ELSE
         qhacw(mgs) =    &
     &   min(   &
     &   ((0.25)*pi)*ehw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf3*xdia(mgs,lh,2)    &
     &    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )     &
     &    , 0.5*(qx(mgs,lc)-qcwresv(mgs))*dtpinv)
!     <  , qxmxd(mgs,lc))
!     <  , qcmxd(mgs))
       
       
         IF ( lwsm6 .and. qsacw(mgs) > 0.0 .and.  qhacw(mgs) > 0.0) THEN
           qaacw = ( qx(mgs,ls)*qsacw(mgs) + qx(mgs,lh)*qhacw(mgs) )/(qx(mgs,ls) + qx(mgs,lh))
!           qaacw = Min( qaacw, 0.5*(qsacw(mgs) + qhacw(mgs) ) )
           qsacw(mgs) = qaacw
           qhacw(mgs) = qaacw
         ENDIF
         
       ENDIF

          qhacwmlr(mgs) = qhacw(mgs)
          IF ( temg(mgs) > tfr .and. iqhacwshr == 0 ) THEN
            qhacw(mgs) = 0.0
          ENDIF
          
          IF ( lvol(lh) .gt. 1 .or. lhl .gt. 1 ) THEN ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             IF ( temg(mgs) .lt. 273.15) THEN
               IF ( irimdenopt == 1 ) THEN ! Heymsfield and Pflaum (1985)
               vt = ( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )
               
             rimdn(mgs,lh) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vt )   &
     &                /(temg(mgs)-273.15))**(rimc2)
!             rimdn(mgs,lh) = Min( Max( hdnmn, rimc3, rimdn(mgs,lh) ), rimc4 )
             rimdn(mgs,lh) = Min( Max( rimc3, rimdn(mgs,lh) ), rimc4 )

!               IF ( igs(mgs) == 30 ) THEN
!                 write(0,*) 'k,vt: ',kgs(mgs),vt, vtxbar(mgs,lh,1),vtxbar(mgs,lh,2), rhovt(mgs)*axx(mgs,lh)*( (alpha(mgs,lh)+3.)*xdia(mgs,lh,1) )**bxx(mgs,lh)
!                 write(0,*) 'diam: char, mean, maxmass = ',xdia(mgs,lh,1),xdia(mgs,lh,3),(alpha(mgs,lh)+3.)*xdia(mgs,lh,1)
!                 write(0,*) 'ax,bx,cd,xdn = ',axx(mgs,lh),bxx(mgs,lh),cdxgs(mgs,lh),xdn(mgs,lh)
!                 write(0,*) 'vt_char,vt_mean = ',rhovt(mgs)*axx(mgs,lh)*( xdia(mgs,lh,1) )**bxx(mgs,lh),rhovt(mgs)*axx(mgs,lh)*( xdia(mgs,lh,3) )**bxx(mgs,lh)
!                 write(0,*) 'rimdn,alpha = ',rimdn(mgs,lh),alpha(mgs,lh)
!               ENDIF

               ELSEIF ( irimdenopt == 2 ) THEN ! Cober and List (1993)

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )   &
     &                /(temg(mgs)-273.15))
                tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) ) ! have to limit range of "R" because quadratic function starts to decrease (unphysically) at higher values
                
                rimdn(mgs,lh) = 1000.*(0.051 + 0.114*tmp - 0.0055*tmp**2)

               ELSEIF ( irimdenopt == 3 .or. irimdenopt == 4) THEN ! Macklin (3) or Saunders and Hosseini 2001

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )   &
     &                /(temg(mgs)-273.15))
              !  tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                IF ( irimdenopt == 3 ) THEN
                  rimdn(mgs,lh) =  Min(900., Max( 170., 110.*tmp**0.76 ) )
                ELSEIF ( irimdenopt == 4 ) THEN ! Saunders and Hosseini
                  rimdn(mgs,lh) =  Min(917., Max( 10.,  900.0*(1.0 - 0.905**tmp ) ) )
                ENDIF
               
               ENDIF
             ELSE
             rimdn(mgs,lh) = 1000.
             ENDIF
             
             IF ( lvol(lh) > 1 ) vhacw(mgs) = rho0(mgs)*qhacw(mgs)/rimdn(mgs,lh)

          ENDIF
      
        IF ( qx(mgs,lh) .gt. qxmin(lh) .and. ipelec .ge. 1 ) THEN
         rarx(mgs,lh) =     &
     &    qhacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lh,2)*cx(mgs,lh))
        ENDIF
      
      ENDIF  
      end do   
!
!
      do mgs = 1,ngscnt
      qhaci(mgs) = 0.0
      qhaci0(mgs) = 0.0
      IF ( ehi(mgs) .gt. 0.0 ) THEN
       IF (  ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,li,1) )

          qhaci0(mgs) = 0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*qx(mgs,li)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da1(li)*xdia(mgs,li,3)**2 ) 
          qhaci(mgs) = Min( ehi(mgs)*qhaci0(mgs), qimxd(mgs) )
       ELSE
        qhaci(mgs) =    &
     &  min(   &
     &  ((0.25)*pi)*ehi(mgs)*ehiclsn(mgs)*qx(mgs,li)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))   &
     &  *(  gf3*xdia(mgs,lh,2)    &
     &    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,li,1)    &
     &    + gf1*xdia(mgs,li,2) )     &
     &  , qimxd(mgs))
       ENDIF
      ENDIF
      end do   


!
!
      do mgs = 1,ngscnt
      qhacs(mgs) = 0.0
      qhacs0(mgs) = 0.0
      IF ( ehs(mgs) .gt. 0.0 ) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,ls,1) )

          qhacs0(mgs) = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*qx(mgs,ls)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            da1(ls)*xdia(mgs,ls,3)**2 ) 
      
          qhacs(mgs) = Min( ehs(mgs)*qhacs0(mgs), qsmxd(mgs) )

       ELSE
         qhacs(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*ehs(mgs)*ehsclsn(mgs)*qx(mgs,ls)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))   &
     &  *(  gf6*gf1*xdia(mgs,ls,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)   &
     &    + gf4*gf3*xdia(mgs,lh,2) )   &
     &  , qsmxd(mgs))
        ENDIF
      ENDIF
      end do   
!
      do mgs = 1,ngscnt
      qhacr(mgs) = 0.0
      qhacrmlr(mgs) = 0.0
      vhacr(mgs) = 0.0
      chacr(mgs) = 0.0
      zhacr(mgs) = 0.0
      IF ( temg(mgs) .gt. tfr ) raindn(mgs,lh) = 1000.0

      IF ( ehr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lr,1) )
!       qhacr(mgs) = ehr(mgs)*cx(mgs,lh)*vt*
!     :     qx(mgs,lr)*0.25*pi*
!     :      (3.02787*xdia(mgs,lr,2) + 
!     :       3.30669*xdia(mgs,lh,1)*xdia(mgs,lr,1) + 
!     :       2.*xdia(mgs,lh,2))
     
       qhacr(mgs) = 0.25*pi*ehr(mgs)*cx(mgs,lh)*qx(mgs,lr)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,lr)*xdia(mgs,lh,3)*xdia(mgs,lr,3) +    &
     &            da1lr(mgs)*xdia(mgs,lr,3)**2 )
!     &            da1(lr)*xdia(mgs,lr,3)**2 )
!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'qhacr= ',qhacr(mgs),tmp
!!        qhacr(mgs) = Min( qhacr(mgs), qrmxd(mgs) )
!!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!!        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

        qhacr(mgs) = Min( qhacr(mgs), qxmxd(mgs,lr) )

            qhacrmlr(mgs) = qhacr(mgs)
        
        IF ( temg(mgs) > tfr .and. iehr0c == 0 ) THEN
          qhacr(mgs) = 0.0

          IF ( iqhacrmlr == 0 ) THEN
              qhacrmlr(mgs) = -qhacw(mgs)
          ENDIF

        ELSE
!        chacr(mgs) = Min( qhacr(mgs)*rho0(mgs)/xmas(mgs,lr), cxmxd(mgs,lr) )

!       chacr(mgs) = ehr(mgs)*cx(mgs,lh)*vt*
!     :     cx(mgs,lr)*0.25*pi*
!     :      (0.69874*xdia(mgs,lr,2) +
!     :       1.24001*xdia(mgs,lh,1)*xdia(mgs,lr,1) +
!     :       2.*xdia(mgs,lh,2))

        chacr(mgs) = 0.25*pi*ehr(mgs)*cx(mgs,lh)*cx(mgs,lr)*vt*      &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +                     &
     &            dab0lh(mgs,lh,lr)*xdia(mgs,lh,3)*xdia(mgs,lr,3) +  &
     &            da0lr(mgs)*xdia(mgs,lr,3)**2 )

!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'chacr= ',chacr(mgs),tmp

!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

      IF ( lzh .gt. 1 ) THEN
          tmp = qx(mgs,lh)/cx(mgs,lh)

!          g1 = (6.0 + alpha(mgs,lh))*(5.0 + alpha(mgs,lh))*(4.0 + alpha(mgs,lh))/
!     :         ((3.0 + alpha(mgs,lh))*(2.0 + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))
!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhacr(mgs) )
      ENDIF
      ENDIF ! temg > tfr
      
      ELSE
       IF ( lwsm6 .and. ipconc == 0 ) THEN
         vt = vt2ave(mgs)
       ELSE
         vt = vtxbar(mgs,lh,1)
       ENDIF

      qhacr(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*ehr(mgs)*qx(mgs,lr)*cx(mgs,lh)   &
     &  *abs(vt-vtxbar(mgs,lr,1))   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lh,1)   &
     &    + gf4*gf3*xdia(mgs,lh,2) )   &
     &  , qrmxd(mgs))
      
        IF ( temg(mgs) > tfr ) THEN
          IF ( iqhacrmlr >= 1 ) qhacrmlr(mgs) = qhacr(mgs)
          qhacr(mgs) = 0.0
        ENDIF
      
      ENDIF
          IF ( lvol(lh) .gt. 1 .or. lhl .gt. 1 ) THEN ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             IF ( temg(mgs) .lt. 273.15) THEN
             raindn(mgs,lh) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lr,3))   &
     &                *((0.60)*vt)   &
     &                /(temg(mgs)-273.15))**(rimc2)

             raindn(mgs,lh) = Min( Max( rimc3, rimdn(mgs,lh) ), rimc4 )
             ELSE
             raindn(mgs,lh) = 1000.
             ENDIF
             
             IF ( lvol(lh) > 1 )  vhacr(mgs) = rho0(mgs)*qhacr(mgs)/raindn(mgs,lh)
        ENDIF
      ENDIF
      end do

!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: hail collects xxxxx'
!

      do mgs = 1,ngscnt
      qhlacw(mgs) = 0.0
      qhlacwmlr(mgs) = 0.0
      vhlacw(mgs) = 0.0
      vhlsoak(mgs) = 0.0
      IF ( lhl > 1 .and. .true.) THEN
        vtmax = 1./(gz(kgs(mgs))*dtp)
        vtxbar(mgs,lhl,1) = Min( vtmax, vtxbar(mgs,lhl,1))
        vtxbar(mgs,lhl,2) = Min( vtmax, vtxbar(mgs,lhl,2))
        vtxbar(mgs,lhl,3) = Min( vtmax, vtxbar(mgs,lhl,3))
      ENDIF

      IF ( lhl > 0 ) THEN
      rarx(mgs,lhl) = 0.0
      ENDIF

      IF ( lhl .gt. 1 .and. ehlw(mgs) .gt. 0.0 ) THEN


!        IF ( ipconc .ge. 2 ) THEN

            vt = abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))

          qhlacw(mgs) = 0.25*pi*ehlw(mgs)*cx(mgs,lhl)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,lc)*xdia(mgs,lhl,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )


          qhlacw(mgs) = Min( qhlacw(mgs), 0.5*qx(mgs,lc)*dtpinv )

          qhlacwmlr(mgs) = qhlacw(mgs)
          IF ( temg(mgs) > tfr .and. iqhlacwshr == 0 ) THEN
            qhlacw(mgs) = 0.0
          ENDIF

          IF ( lvol(lhl) .gt. 1 ) THEN

             IF ( temg(mgs) .lt. 273.15) THEN
               IF ( irimdenopt == 1 ) THEN ! Heymsfeld and Pflaum (1985)
             rimdn(mgs,lhl) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) ))   &
     &                /(temg(mgs)-273.15))**(rimc2)
             rimdn(mgs,lhl) = Min( Max( hldnmn, rimc3, rimdn(mgs,lhl) ), rimc4 )
               
               ELSEIF ( irimdenopt == 2 ) THEN ! Cober and List (1993)
                tmp = -0.5*(1.e+06)*xdia(mgs,lc,1)   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) )   &
     &                /(temg(mgs)-273.15)
                tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                rimdn(mgs,lhl) = 1000.*(0.051 + 0.114*tmp - 0.005*tmp**2)
               
               ELSEIF ( irimdenopt == 3 .or. irimdenopt == 4) THEN ! Macklin (3) or Saunders and Hosseini 2001
                tmp = -0.5*(1.e+06)*xdia(mgs,lc,1)   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) )  &
     &                /(temg(mgs)-273.15)
              !  tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                IF ( irimdenopt == 3 ) THEN
                  rimdn(mgs,lhl) =  Min(900., Max( 170., 110.*tmp**0.76 ) )
                ELSEIF ( irimdenopt == 4 ) THEN ! Saunders and Hosseini
                  rimdn(mgs,lhl) =  Min(917., Max( 10.,  900.0*(1.0 - 0.905**tmp ) ) )
                ENDIF
               
               ENDIF
             ELSE
             rimdn(mgs,lhl) = 1000.
             ENDIF

             vhlacw(mgs) = rho0(mgs)*qhlacw(mgs)/rimdn(mgs,lhl)

          ENDIF


        IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. ipelec .ge. 1 ) THEN
         rarx(mgs,lhl) =     &
     &    qhlacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lhl,2)*cx(mgs,lhl))
        ENDIF

      ENDIF
      end do

      qhlaci(:) = 0.0
      qhlaci0(:) = 0.0
      IF ( lhl .gt. 1  ) THEN
      do mgs = 1,ngscnt
      IF ( ehli(mgs) .gt. 0.0 ) THEN
       IF (  ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,li,1) )

          qhlaci0(mgs) = 0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*qx(mgs,li)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            da1(li)*xdia(mgs,li,3)**2 )
        ! qhlaci(mgs) = Min( qhlaci(mgs), qimxd(mgs) )
          qhlaci(mgs) = Min( ehli(mgs)*qhlaci0(mgs), qimxd(mgs) )
       ENDIF
      ENDIF
      end do
      ENDIF
!
      qhlacis(:) = 0.0
      qhlacis0(:) = 0.0

      qhlacs(:) = 0.0
      qhlacs0(:) = 0.0
      IF ( lhl .gt. 1 ) THEN
      do mgs = 1,ngscnt
      IF ( ehls(mgs) .gt. 0.0) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,ls,1) )

          qhlacs0(mgs) = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*qx(mgs,ls)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            da1(ls)*xdia(mgs,ls,3)**2 )

          qhlacs(mgs) = Min( ehls(mgs)*qhlacs0(mgs), qsmxd(mgs) )
        ENDIF
      ENDIF
      end do
      ENDIF


      do mgs = 1,ngscnt
      qhlacr(mgs) = 0.0
      qhlacrmlr(mgs) = 0.0
      chlacr(mgs) = 0.0
      vhlacr(mgs) = 0.0
      IF ( lhl .gt. 1 .and. temg(mgs) .gt. tfr ) raindn(mgs,lhl) = 1000.0

      IF ( lhl .gt. 1 .and. ehlr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,lr,1) )

       qhlacr(mgs) = 0.25*pi*ehlr(mgs)*cx(mgs,lhl)*qx(mgs,lr)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,lr)*xdia(mgs,lhl,3)*xdia(mgs,lr,3) +    &
     &            da1lr(mgs)*xdia(mgs,lr,3)**2 )
!     &            da1(lr)*xdia(mgs,lr,3)**2 )
!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'qhacr= ',qhacr(mgs),tmp
!!        qhacr(mgs) = Min( qhacr(mgs), qrmxd(mgs) )
!!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!!        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

        qhlacr(mgs) = Min( qhlacr(mgs), qxmxd(mgs,lr) )

     
        IF ( iqhlacrmlr >= 1 ) qhlacrmlr(mgs) = qhlacr(mgs)
        
        IF ( temg(mgs) > tfr .and. iehlr0c == 0) THEN
          qhlacr(mgs) = 0.0
          IF ( iqhlacrmlr == 0 ) THEN
              qhlacrmlr(mgs) = -qhlacw(mgs)
          ENDIF
        ELSE
        chlacr(mgs) = 0.25*pi*ehlr(mgs)*cx(mgs,lhl)*cx(mgs,lr)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0lh(mgs,lhl,lr)*xdia(mgs,lhl,3)*xdia(mgs,lr,3) +    &
     &            da0lr(mgs)*xdia(mgs,lr,3)**2 )

        chlacr(mgs) = min(chlacr(mgs),crmxd(mgs))

        IF ( lvol(lhl) .gt. 1 ) THEN
         vhlacr(mgs) = rho0(mgs)*qhlacr(mgs)/raindn(mgs,lhl)
        ENDIF
        ENDIF
      ENDIF
      ENDIF
      end do



!
!
!
!
!      if (ndebug .gt. 0 ) write(0,*) 'Collection: Cloud collects xxxxx'

      if (ndebug .gt. 0 ) write(0,*) 'Collection: cloud ice collects xxxx2'
!
      do mgs = 1,ngscnt
      qiacw(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 ) THEN

       vt = Sqrt((vtxbar(mgs,li,1)-vtxbar(mgs,lc,1))**2 +    &
     &            0.04*vtxbar(mgs,li,1)*vtxbar(mgs,lc,1) )

          qiacw(mgs) = 0.25*pi*eiw(mgs)*cx(mgs,li)*qx(mgs,lc)*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab1(li,lc)*xdia(mgs,li,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )

       qiacw(mgs) = Min( qiacw(mgs), qxmxd(mgs,lc) )
      ENDIF
      end do


!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: cloud ice collects xxxx8'
!
      do mgs = 1,ngscnt
      qiacr(mgs) = 0.0
      qiacrf(mgs) = 0.0
      qiacrs(mgs) = 0.0
      ciacrs(mgs) = 0.0
      ciacr(mgs) = 0.0
      ciacrf(mgs) = 0.0
      viacrf(mgs) = 0.0
      csplinter(mgs) = 0.0
      qsplinter(mgs) = 0.0
      csplinter2(mgs) = 0.0
      qsplinter2(mgs) = 0.0
      IF ( iacr .ge. 1 .and. eri(mgs) .gt. 0.0    &
     &     .and. temg(mgs) .le. 270.15 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       ni = 0.0
         IF ( xdia(mgs,li,1) .ge. 10.e-6 ) THEN
          ni = ni + cx(mgs,li)*Exp(- (40.e-6/xdia(mgs,li,1))**3 )
         ENDIF
       IF ( imurain == 1 ) THEN ! gamma of diameter
         IF ( iacrsize /= 4 ) THEN
           IF ( iacrsize .eq. 1 ) THEN
             ratio = 500.e-6/xdia(mgs,lr,1)
           ELSEIF ( iacrsize .eq. 2 ) THEN
             ratio = 300.e-6/xdia(mgs,lr,1)
           ELSEIF ( iacrsize .eq. 3 ) THEN
             ratio = 40.e-6/xdia(mgs,lr,1)
           ELSEIF ( iacrsize .eq. 5 ) THEN
             ratio = 150.e-6/xdia(mgs,lr,1)
           ELSEIF ( iacrsize .eq. 6 ) THEN
             ratio = 60.e-6/xdia(mgs,lr,1)
             ni = cx(mgs,li)
           ENDIF
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
!           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))
           
           ! interpolate along alpha
           
           nr = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)
           
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qr = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)
           
         ELSE ! iacrsize == 4 : use all
           nr = cx(mgs,lr)
           qr = qx(mgs,lr)
         ENDIF

          vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +     &
     &            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )

          qiacr(mgs) = 0.25*pi*eri(mgs)*ni*qr*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab1lh(mgs,li,lr)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da1(lr)*xdia(mgs,lr,3)**2 ) 

          qiacr(mgs) = Min( qrmxd(mgs), qiacr(mgs) )


          ciacr(mgs) = 0.25*pi*eri(mgs)*ni*nr*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab0lh(mgs,li,lr)*xdia(mgs,lr,3)*xdia(mgs,li,3) +    &
     &            da0(lr)*xdia(mgs,lr,3)**2 ) 

          ciacr(mgs) = Min( crmxd(mgs), ciacr(mgs) )
          
!          write(iunit,*) 'qiacr: ',cx(mgs,lr),nr,qx(mgs,lr),qr,qiacr(mgs),ciacr(mgs)
!          write(iunit,*) 'xdia r li = ',xdia(mgs,lr,3),xdia(mgs,li,3),xdia(mgs,lr,1),xdia(mgs,li,1)
!          write(iunit,*) 'i,j,ratio = ',i,j,ciacrratio(i,j),qiacrratio(i,j)
!          write(iunit,*) 'ni,ci = ',ni,cx(mgs,li),qx(mgs,li)

       ELSEIF ( imurain == 3 ) THEN ! gamma of volume
!   Set nr to the number of drops greater than 40 microns.
         arg = 1000.*xdia(mgs,lr,3)
!         nr = cx(mgs,lr)*gaml02( arg )
!        IF ( iacr .eq. 1 ) THEN
         IF ( ipconc .ge. 3 ) THEN
           IF ( iacrsize .eq. 1 ) THEN
            nr = cx(mgs,lr)*gaml02d500( arg )  ! number greater than 500 microns in diameter
           ELSEIF ( iacrsize .eq. 2 .or. iacrsize .eq. 5 ) THEN
            nr = cx(mgs,lr)*gaml02d300( arg )  ! number greater than 300 microns in diameter
           ELSEIF ( iacrsize .eq. 3 ) THEN
            nr = cx(mgs,lr)*gaml02( arg ) ! number greater than 40 microns in diameter
           ELSEIF ( iacrsize .eq. 4 ) THEN
            nr = cx(mgs,lr) ! all raindrops
           ENDIF
         ELSE
         nr = cx(mgs,lr)*gaml02( arg )
         ENDIF
!        ELSEIF ( iacr .eq. 2 ) THEN
!         nr = cx(mgs,lr)*gaml02d300( arg )  ! number greater than 300 microns in diameter
!        ENDIF
       IF ( ni .gt. 0.0 .and. nr .gt. 0.0 ) THEN
       d0 = xdia(mgs,lr,3)
       qiacr(mgs) = xdn(mgs,lr)*rhoinv(mgs)*   &
     &     (0.217239*(0.522295*(d0**5) +    &
     &      49711.81*(d0**6) -    &
     &      1.673016e7*(d0**7)+    &
     &      2.404471e9*(d0**8) -    &
     &      1.22872e11*(d0**9))*ni*nr)
      qiacr(mgs) = Min( qrmxd(mgs), qiacr(mgs) )
      ciacr(mgs) =   &
     &   (0.217239*(0.2301947*(d0**2) +    &
     &      15823.76*(d0**3) -    &
     &      4.167685e6*(d0**4) +    &
     &      4.920215e8*(d0**5) -    &
     &      2.133344e10*(d0**6))*ni*nr)
      ciacr(mgs) = Min( crmxd(mgs), ciacr(mgs) )
!      ciacr(mgs) = qiacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
      ENDIF
      ENDIF
       IF ( iacr .eq. 1 .or. iacr .eq. 3 ) THEN
         ciacrf(mgs) = Min(ciacr(mgs), qiacr(mgs)/(1.0*vr1mm*1000.0)*rho0(mgs) ) ! *rzxh(mgs)
       ELSEIF ( iacr .eq. 2 ) THEN
         ciacrf(mgs) = ciacr(mgs) ! *rzxh(mgs)
       ELSEIF ( iacr .eq. 4 ) THEN
         ciacrf(mgs) = Min(ciacr(mgs), qiacr(mgs)/(1.0*vfrz*1000.0)*rho0(mgs) ) ! *rzxh(mgs)
       ELSEIF ( iacr .eq. 5 ) THEN
         ciacrf(mgs) = ciacr(mgs)*rzxh(mgs)
       ENDIF 
!      crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*27.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
       ENDIF
      
      
      ELSE ! single-moment rain
      qiacr(mgs) =    &
     &  min(        &
     &   ((0.25/gf4)*pi)*eri(mgs)*cx(mgs,li)*qx(mgs,lr)   &
     &  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)    &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1)    &
     &    + gf4*gf3*xdia(mgs,li,2) )     &
     &  , qrmxd(mgs))
      ENDIF
!      if ( temg(mgs) .gt. 268.15 ) then
!      qiacr(mgs) = 0.0
!      ciacr(mgs) = 0.0
!      end if

      IF ( ipconc .ge. 1 ) THEN
        IF ( nsplinter .ge. 1000 ) THEN
        ! Lawson et al. 2015 JAS
         ! ave. diam of freezing drops in microns
           IF ( qiacr(mgs)*dtp > qxmin(lh) .and. ciacr(mgs) > 1.e-3 ) THEN
             tmpdiam = 1.e6*( 6.*qiacr(mgs)/(1000.*pi*ciacr(mgs) ) )**(1./3.) ! avg. diameter of newly frozen drops in microns
              fac = 1.0
              IF ( nsplinter .eq. 1001 ) THEN
             !   fac = 0.2/sqrt(2.0*pi*10.**2)*Exp(-0.5*((258.-temg(mgs))/10.)**2 ) ! temperature dependence from Sullivan et al. 2018 ACP
             ! ELSE
                fac = 0.2*Exp(-0.5*((258.-temg(mgs))/10.)**2 ) ! temperature dependence from Sullivan et al. 2018 ACP
              ENDIF
             csplinter(mgs) = fac*lawson_splinter_fac*tmpdiam**4*ciacr(mgs)
           ENDIF
        ELSEIF ( nsplinter .ge. 0 ) THEN
          csplinter(mgs) = nsplinter*ciacr(mgs)
        ELSE
          csplinter(mgs) = -nsplinter*ciacrf(mgs)
        ENDIF
        qsplinter(mgs) = Min(0.1*qiacr(mgs), csplinter(mgs)*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel
      ENDIF
      
      frach = 1.0
           IF ( ibiggsnow == 2 .or. ibiggsnow == 3 ) THEN
           IF ( ciacr(mgs) > qxmin(lh) ) THEN
           xvfrz = rho0(mgs)*qiacr(mgs)/(ciacr(mgs)*900.) ! mean volume of frozen drops; 900. for frozen drop density
           frach = 0.5 *(1. +  Tanh(0.2e12 *( xvfrz - 1.15*xvbiggsnow)))

             qiacrs(mgs) = (1.-frach)*qiacr(mgs)
             ciacrs(mgs) = (1.-frach)*ciacrf(mgs) ! *rzxh(mgs)
           
           ENDIF
           ENDIF

      qiacrf(mgs) = frach*qiacr(mgs)
      ciacrf(mgs) = frach*ciacrf(mgs)

      IF ( lvol(lh) > 1 ) THEN
         viacrf(mgs) = rho0(mgs)*qiacrf(mgs)/rhofrz
      ENDIF
      
      end do
!
!
!
!

! snow aggregation here
      if ( ipconc .ge. 4 ) then !
      do mgs = 1,ngscnt
      csacs(mgs) = 0.0
      IF ( qx(mgs,ls) > qxmin(ls) .and. ess(mgs) .gt. 0.0 ) THEN ! .and. xv(mgs,ls) < 0.25*xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls)))  ) THEN

        IF ( iessec0flag == 0 ) THEN
          ec0(mgs) = 1.0
        ELSE
          tmp = xv(mgs,ls)/(xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls)))) ! fraction of max snow mass
          IF ( tmp .lt. essfrac1 ) THEN
            ec0(mgs) = 1.0
          ELSEIF ( tmp .ge. essfrac2 ) THEN
            ec0(mgs) = 0.0
          ELSE
            ec0(mgs) = (essfrac2 - tmp)/(essfrac2 - essfrac1)
          ENDIF
        ENDIF

      csacs(mgs) = ec0(mgs)*rvt*aa2*ess(mgs)*cx(mgs,ls)**2*Min( xv(mgs,ls), 4.*pii/3.*essrmax**3 ) ! *Min(1.,xdn(mgs,ls)/100. ) ! Min func tries to recalibrate for low diagnosed density 
!      csacs(mgs) = rvt*aa2*ess(mgs)*cx(mgs,ls)**2*Min( xv(mgs,ls), 4.*pii/3.*0.02**3 ) ! *Min(1.,xdn(mgs,ls)/100. ) ! Min func tries to recalibrate for low diagnosed density 
      csacs(mgs) = Min(csacs(mgs),csmxd(mgs))
      ENDIF
      end do
      end if
!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 11'
      if ( ipconc .ge. 2 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      ciacw(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 .and. xmas(mgs,lc) > 0.0 ) THEN
        ciacw(mgs) = qiacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        ciacw(mgs) = min(ciacw(mgs),ccmxd(mgs))
      ENDIF
      end do

      end if

      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 18'
      if ( ipconc .ge. 2 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
       tmp1 = 0.0
       cracw(mgs) = 0.0
       cracr(mgs) = 0.0
       zracr(mgs) = 0.0
       ec0(mgs) = 1.e9
      IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr)    &
     &      .and. qracw(mgs) .gt. 0.0 ) THEN

       IF ( ipconc .lt. 3 ) THEN
        IF ( erw(mgs) .gt. 0.0 ) THEN
        cracw(mgs) =   &
     &   ((0.25)*pi)*erw(mgs)*(cx(mgs,lc) - ccwresv(mgs))*cx(mgs,lr)   &
     &  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf1*xdia(mgs,lc,2)   &
     &    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lr,1)   &
     &    + gf3*xdia(mgs,lr,2) )
        ENDIF
       ELSE ! IF ( ipconc .ge. 3 .and. )
        IF ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) THEN  !{
        IF ( 0.5*xdia(mgs,lr,3) .gt. rh(mgs) ) THEN ! { .or. cx(mgs,lr) .gt. nh(mgs) 
!        IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
          IF ( 0.5*xdia(mgs,lr,3) .gt. rwradmn ) THEN ! r > 50.e-6 
!          DM0CCC=A2*XNC*XNR*(XVC+XVR)                               ! (A11)
!         NOTE: murain drops out, so same result for imurain = 1 and 3
            cracw(mgs) = aa2*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*(xv(mgs,lc) + xv(mgs,lr))
          ELSE
            IF ( imurain == 3 ) THEN
!          DM0CCC=A1*XNC*XNR*(((CNU+2.)/(CNU+1.))*XVC**2+((RNU+2.)/(RNU+1.))*XVR**2) ! (A13)
            cracw(mgs) = aa1*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*   &
     &          ((alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.) +    &
     &          (alpha(mgs,lr) + 2.)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))
            ELSE ! imurain == 1 USE CP00 for rain DSD in diameter
            cracw(mgs) = aa1*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*   &
     &          ((alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.) +    &
     &          (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)*xv(mgs,lr)**2/  &
     &             ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.)) )
            ENDIF ! imurain
          ENDIF
        ENDIF ! } rh
        ENDIF ! } dmrauto
       ENDIF ! ipconc
      ENDIF ! qc > qcmin & qr > qrmin
        
! Rain self collection (cracr) and break-up (factor of ec0)
!
!       
        ec0(mgs) = 1.0 ! 2.e9
        IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        rwrad = 0.5*xdia(mgs,lr,3)
        
        
        ! check median volume diameter
        IF ( icracrthresh > 1 ) THEN
         IF ( imurain == 1 ) THEN
           tmp =  (3.67+alpha(mgs,lr))*xdia(mgs,lr,1) ! median volume diameter; units of m (Ulbrich 1983, JCAM)
         ELSE ! imurain == 3, 
           tmp =  (1.678+alpha(mgs,lr))**(1./3.)*xdia(mgs,lr,1) ! units of mm (using method of Ulbrich 1983. See ventillation_stuff.nb)
         ENDIF
        ELSE
          tmp = xdia(mgs,lr,3) - 0.1e-3
        ENDIF
        tmpdiam = tmp
         
!    Using collection efficiency factor ec0 to simulate break-up that off-sets self-collection (Zieger 1985; Cohard & Pinty 2000)
!    ec0 is 1 for rain diameter < 600 microns and then drop off toward zero until diameter of 2mm to represent passive breakup
!    ec0 does not go negative here (i.e., does not follow other versions that create extra breakup at large rain diameter)
        IF ( ( tmpdiam .gt. 1.9e-3 .and. irainbreak /= 10 .and. irainbreak /= 20 ) .or. icracr <= 0  ) THEN
          ec0(mgs) = 0.0
          cracr(mgs) = 0.0
          IF ( ibincracr == 3 ) THEN
          tmp1 = aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)/ &
     &                  ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.))
          ENDIF
        ELSE
         IF ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) THEN 
          
          IF ( tmpdiam .lt. 6.1e-4 .or. irainbreak == 10 ) THEN
            ec0(mgs) = 1.0
          ELSE
            ec0(mgs) = Exp( -2500.0*(tmpdiam - 6.0e-4) )
          ENDIF

          IF ( rwrad .ge. 50.e-6 ) THEN
              tmp1 = aa2*cx(mgs,lr)**2*xv(mgs,lr)
              cracr(mgs) = ec0(mgs)*tmp1
              IF ( irainbreak == 20 ) THEN
                cracr(mgs) = tmp1
              ENDIF
          ELSE
            IF ( imurain == 3 ) THEN
             cracr(mgs) = ec0(mgs)*aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 2.)/(alpha(mgs,lr) + 1.)
            ELSE ! imurain == 1
              tmp1 = aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)/ &
     &                  ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.))
              cracr(mgs) = ec0(mgs)*tmp1
              IF ( irainbreak == 20 ) THEN
                 cracr(mgs) = tmp1
              ENDIF
            ENDIF
          ENDIF ! rwrad > 50
!          cracr(mgs) = Min(cracr(mgs),crmxd(mgs))
         ENDIF ! dmrauto <= 0
        ENDIF ! tmp > 1.9e-3
          
          IF ( irainbreak == 100 ) THEN ! Morrison breakup
            ec0(mgs) = 1.0
            IF ( xdia(mgs,lr,1) > 300.e-6 ) THEN
              ec0(mgs) = 2. - Exp(2300.*(xdia(mgs,lr,1)-300.e-6))
            ENDIF
            cracr(mgs) = 5.78*ec0(mgs)*cx(mgs,lr)*qx(mgs,lr)
          ENDIF
        
        ENDIF ! ( qx(mgs,lr) .gt. qxmin(lr) )

       ! active breakup option
       crbreak = 0.0
       IF ( irainbreak == 1 .or. irainbreak == 10 ) THEN
                crbreak = Max( 0.0,  rainbreakfac* (rho0(mgs)*qx(mgs,lr))**2 ) ! hand fit to lower range of wkqss output
                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 2 .or. irainbreak == 20 .or. irainbreak == 12 ) THEN
          ! irainbreak == 20 does not work as intended
              IF (  irainbreak == 12 ) THEN
                IF ( xdia(mgs,lr,1) > 300.e-6 ) THEN
                  crbreak = Max( 0.0,  rainbreakfac*(rho0(mgs)*qx(mgs,lr))**2 ) ! hand fit to lower range of wkqss output
                ELSE
                  crbreak = 0.0
                ENDIF
              ELSE
                crbreak = Max( 0.0,  rainbreakfac*(1. - ec0(mgs))*(rho0(mgs)*qx(mgs,lr))**2 ) ! hand fit to lower range of wkqss output
              ENDIF
!                crbreak = Max(0.0, -0.18 + 1.139e6 * (rho0(mgs)*qx(mgs,lr) + 0.00038106)**2)
                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 3 .and. qx(mgs,lr) > qxmin(lr) .and. ipconc >= 5  ) THEN
         
          ! IF ( irainbreak == 1 ) THEN
           
           ratio = Min( maxratiolu, 10.e-3/xdia(mgs,lr,1) )
           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lr),4,1)
            qxd1 = qx(mgs,lr)*(tmp2)

            IF ( ( qxd1 > qxmin(lr)) ) THEN
            
           ! number
              tmp = gaminterp(ratio,alpha(mgs,lr),1,1)
              cxd1 = cx(mgs,lr)*( tmp)
            
              IF ( cxd1 > 2.6e-6 ) THEN
            
                crbreak = 10.0**(2.1967 + 0.8177*log10(cxd1)) ! fit of log-log bin data
            
!               IF ( kgs(mgs) == 1 .and. qx(mgs,lr) > 0.1e-3 ) THEN
!                 write(0,*) 'crbreak: ',crbreak,crbreaksmall,dtpinv*cxd1,cx(mgs,lr),cracr(mgs) - crbreak
!               ENDIF
                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
              ENDIF
            ENDIF
       ELSEIF ( irainbreak == 4 .and. qx(mgs,lr) > qxmin(lr) .and. ipconc >= 5  ) THEN
           crbreak = Exp(14.7863 + 2.2232*log(rho0(mgs)*qx(mgs,lr))) ! "best" fit
           cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 5 .and. qx(mgs,lr) > qxmin(lr) .and. ipconc >= 5  ) THEN
           crbreak = Exp(11.42 + 1.543*log(rho0(mgs)*qx(mgs,lr))) ! hand fit to get closer at large qr
           cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 6 .and. qx(mgs,lr) > qxmin(lr) .and. ipconc >= 5  ) THEN
           crbreak = Exp(14.0 + 1.85*log(rho0(mgs)*qx(mgs,lr))) ! hand fit from wkqsszh3m800cn2km5irtbreak1dsedr1ibinbreak.crbin.txt
           cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 7 ) THEN
                crbreak = Max( 0.0, 2.0e6 * (rho0(mgs)*qx(mgs,lr))**2 ) ! hand fit to lower range of wkqss output
                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 8 ) THEN
                crbreak = Max( 0.0, 2.542e6 * (rho0(mgs)*qx(mgs,lr))**2 ) ! best fit to lower range of wkqss (collision only) output
                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
       ELSEIF ( irainbreak == 11 .and. rho0(mgs)*qx(mgs,lr) > qrbrthresh1 .and. ipconc >= 5  ) THEN
         
          ! Ad hoc method to break up drops in the DSD tail (D > draintail)
           
           ratio = Min( maxratiolu, draintail/xdia(mgs,lr,1) )
           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lr),4,1)
            qxd1 = qx(mgs,lr)*(tmp2)
            qrbreak = dtpinv*qxd1

            crbreaksmall = rho0(mgs)*qrbreak/(xdn(mgs,lr)*pi/6.*drsmall**3)
            IF ( ( qxd1 > qxmin(lr)) ) THEN
            
           ! number
            tmp = gaminterp(ratio,alpha(mgs,lr),1,1)
             IF ( ipconc == 5 ) THEN
          !     tmp = Min( 0.2, tmp )
             ENDIF
            cxd1 = cx(mgs,lr)*( tmp)
            IF ( rho0(mgs)*qx(mgs,lr) > qrbrthresh2 ) THEN
              flim = 1.0
            ELSE
              flim = (rho0(mgs)*qx(mgs,lr) - qrbrthresh1)/(qrbrthresh2 - qrbrthresh1)
            ENDIF
            crbreak = flim*(crbreaksmall - dtpinv*cxd1)
            
!             IF ( kgs(mgs) == 1 .and. qx(mgs,lr) > 0.1e-3 ) THEN
!               write(0,*) 'crbreak: ',crbreak,crbreaksmall,dtpinv*cxd1,cx(mgs,lr),cracr(mgs) - crbreak
!             ENDIF
            cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup

           ! reflectivity -- not used yet: goes into zracr
!            IF ( ipconc >= 6 .and. lzr > 1  ) THEN
!             tmp3 = gaminterp(ratio,alpha(mgs,lr),11,1)
!             zxd1 = zx(mgs,lr)*(tmp3)
!             zrbreak = dtpinv*zxd1
!            ELSE
!             zxd1 = 0
!            ENDIF
!            zrbreak = Max(0.0, zrbreak - crbreaksmall*drsmall**6)
!       ELSEIF ( irainbreak == 12 ) THEN
!                crbreak = Max( 0.0, 3.8098 * (rho0(mgs)*qx(mgs,lr))**1.9416 ) ! best fit to lower range of wkqss (collision only) output
!                cracr(mgs) = cracr(mgs) - crbreak ! cracr is subtracted, so negative value for breakup
             ENDIF
       ENDIF

       IF ( lzr > 0 .and. cracr(mgs) /= 0.0 .and. cx(mgs,lr) > 0.0  ) THEN
!          tmp = qx(mgs,lr)/cx(mgs,lr)
!          zracr(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*1000.))**2*( tmp**2 * cracr(mgs) )
        ! rewrite because original can overestimate zracr if -cracr*dtp is on the order of cx (i.e.,
        !  large increase in the number of drops, which violates differential assumption
        ! Pass -cracr because its meaning is backwards (neg. value ADDS number, positive value SUBTRACTS)
          zracr(mgs) = zraten(dtpinv,dtp,g1x(mgs,lr),rho0(mgs),rho_qr,qx(mgs,lr),cx(mgs,lr),-cracr(mgs))

!          zracr(mgs) = dtpinv*g1x(mgs,lr)*(6.*rho0(mgs)*qx(mgs,lr)/(pi*1000.))**2 &
!                     * ( cracr(mgs) )/((cx(mgs,lr) - dtp*cracr(mgs))*(cx(mgs,lr)))
       ENDIF

!      cracw(mgs) = min(cracw(mgs),cxmxd(mgs,lc)) 

      end do
      end if

!
!  frozen drops
!
      IF ( lf > 1 ) THEN
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22ii'
      cfacw(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt

      IF ( ipconc .ge. 5 ) THEN
       IF ( qfacw(mgs) .gt. 0.0 .and. xmas(mgs,lc) .gt. 0.0 ) THEN

!  This is the explict version of chacw, which turns out to be very close to the
!  approximation that the droplet size does not change, to within a few percent.
!  This may _not_ be the case for cnu other than zero!
!          cfacw(mgs) = (efw(mgs)*cx(mgs,lc)*cx(mgs,lf)*(pi/4.)*
!     :    abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))*
!     :    (2.0*xdia(mgs,lf,1)*(xdia(mgs,lf,1) +
!     :         xdia(mgs,lc,1)*gf43rds) +
!     :      xdia(mgs,lc,2)*gf53rds))

!          cfacw(mgs) = Min( cfacw(mgs), 0.6*cx(mgs,lc)*dtpinv )

!        cfacw(mgs) = qfacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        cfacw(mgs) = qfacw(mgs)*rho0(mgs)/xmascw(mgs)
!        cfacw(mgs) = min(cfacw(mgs),cxmxd(mgs,lc))
        cfacw(mgs) = Min( cfacw(mgs), 0.5*(cx(mgs,lc) - ccwresv(mgs))*dtpinv )
       ELSE
        qfacw(mgs) = 0.0
       ENDIF
      ENDIF
      end do
      end if

     if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      cfaci(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( efi(mgs) .gt. 0.0 .or. ( eficlsn(mgs) > 0.0 .and. ipelec > 0 )) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,li,1) )

          cfaci0(mgs) = 0.25*pi*eficlsn(mgs)*cx(mgs,lf)*cx(mgs,li)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab0lh(mgs,lf,li)*xdia(mgs,lf,3)*xdia(mgs,li,3) +    &
     &            da0(li)*xdia(mgs,li,3)**2 )

        ENDIF

        cfaci(mgs) = min(efi(mgs)*cfaci0(mgs),cimxd(mgs))

       ENDIF
      end do
      end if


!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22nn'
      cfacs(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( efs(mgs) .gt. 0 ) THEN
       IF ( ipconc .ge. 5 .or. ( efsclsn(mgs) > 0.0 .and. ipelec > 0 ) ) THEN

       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,ls,1) )

          cfacs0(mgs) = 0.25*pi*efsclsn(mgs)*cx(mgs,lf)*cx(mgs,ls)*vt*   &
     &         (  da0lx(mgs,lf)*xdia(mgs,lf,3)**2 +     &
     &            dab0lh(mgs,lf,ls)*xdia(mgs,lf,3)*xdia(mgs,ls,3) +    &
     &            da0(ls)*xdia(mgs,ls,3)**2 )

      ENDIF
      cfacs(mgs) = min(efs(mgs)*cfacs0(mgs),csmxd(mgs))
      ENDIF
      end do
      end if


      ELSE
        cfacw(:) = 0.0
        cfaci(:) = 0.0
        cfaci0(:) = 0.0
        cfacis(:) = 0.0
        cfacis0(:) = 0.0
        cfacs(:) = 0.0
        cfacs0(:) = 0.0
      ENDIF ! lf > 1
!
!
!
!  Graupel
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22ii'
      chacw(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt

      IF ( ipconc .ge. 5 ) THEN
       IF ( qhacw(mgs) .gt. 0.0 .and. xmas(mgs,lc) .gt. 0.0 ) THEN

!  This is the explict version of chacw, which turns out to be very close to the
!  approximation that the droplet size does not change, to within a few percent.
!  This may _not_ be the case for cnu other than zero!
!          chacw(mgs) = (ehw(mgs)*cx(mgs,lc)*cx(mgs,lh)*(pi/4.)*
!     :    abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*
!     :    (2.0*xdia(mgs,lh,1)*(xdia(mgs,lh,1) +
!     :         xdia(mgs,lc,1)*gf43rds) +
!     :      xdia(mgs,lc,2)*gf53rds))

!          chacw(mgs) = Min( chacw(mgs), 0.6*cx(mgs,lc)*dtpinv )

!        chacw(mgs) = qhacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        chacw(mgs) = qhacw(mgs)*rho0(mgs)/xmascw(mgs)
!        chacw(mgs) = min(chacw(mgs),cxmxd(mgs,lc))
        chacw(mgs) = Min( chacw(mgs), 0.5*(cx(mgs,lc) - ccwresv(mgs))*dtpinv )
       ELSE
        qhacw(mgs) = 0.0
       ENDIF
      ELSE
      ! single-moment
      chacw(mgs) =   &
     &   ((0.25)*pi)*ehw(mgs)*cx(mgs,lc)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf1*xdia(mgs,lc,2)   &
     &    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lh,1)   &
     &    + gf3*xdia(mgs,lh,2) )
      chacw(mgs) = min(chacw(mgs),0.5*cx(mgs,lc)*dtpinv)
!      chacw(mgs) = min(chacw(mgs),cxmxd(mgs,lc))
!      chacw(mgs) = min(chacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      chaci(:) = 0.0
      chaci0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( ehi(mgs) .gt. 0.0 .or. ( ehiclsn(mgs) > 0.0 .and. ipelec > 0 )) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,li,1) )

          chaci0(mgs) = 0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*cx(mgs,li)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab0lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da0(li)*xdia(mgs,li,3)**2 )

       ELSE
        chaci0(mgs) =   &
     &   ((0.25)*pi)*ehiclsn(mgs)*cx(mgs,li)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))   &
     &  *(  gf1*xdia(mgs,li,2)   &
     &    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lh,1)   &
     &    + gf3*xdia(mgs,lh,2) )
        ENDIF

        chaci(mgs) = min(ehi(mgs)*chaci0(mgs),cimxd(mgs))
       ENDIF
      end do
      end if


!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22nn'
      chacs(:) = 0.0
      chacs0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( ehs(mgs) .gt. 0 ) THEN
       IF ( ipconc .ge. 5 .or. ( ehsclsn(mgs) > 0.0 .and. ipelec > 0 ) ) THEN

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,ls,1) )

          chacs0(mgs) = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*cx(mgs,ls)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab0lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            da0(ls)*xdia(mgs,ls,3)**2 )

       ELSE
      chacs0(mgs) =   &
     &   ((0.25)*pi)*ehsclsn(mgs)*cx(mgs,ls)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))   &
     &  *(  gf3*gf1*xdia(mgs,ls,2)   &
     &    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)   &
     &    + gf1*gf3*xdia(mgs,lh,2) )
      ENDIF
      chacs(mgs) = min(ehs(mgs)*chacs0(mgs),csmxd(mgs))
      ENDIF
      end do
      end if


!
!
!  Hail
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22ii'
      chlacw(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt

      IF ( lhl .gt. 1 .and. ipconc .ge. 5 ) THEN
       IF ( qhlacw(mgs) .gt. 0.0 .and. xmas(mgs,lc) .gt. 0.0 ) THEN

!  This is the explict version of chacw, which turns out to be very close to the
!  approximation that the droplet size does not change, to within a few percent.
!  This may _not_ be the case for cnu other than zero!
!          chlacw(mgs) = (ehlw(mgs)*cx(mgs,lc)*cx(mgs,lhl)*(pi/4.)*
!     :    abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))*
!     :    (2.0*xdia(mgs,lhl,1)*(xdia(mgs,lhl,1) +
!     :         xdia(mgs,lc,1)*gf43rds) +
!     :      xdia(mgs,lc,2)*gf53rds))

!          chlacw(mgs) = Min( chlacw(mgs), 0.6*cx(mgs,lc)*dtpinv )

!        chlacw(mgs) = qhlacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        chlacw(mgs) = qhlacw(mgs)*rho0(mgs)/xmascw(mgs)
!        chlacw(mgs) = min(chlacw(mgs),cxmxd(mgs,lc))
        chlacw(mgs) = Min( chlacw(mgs), 0.5*cx(mgs,lc)*dtpinv )
       ELSE
        qhlacw(mgs) = 0.0
       ENDIF
!      ELSE
!      chlacw(mgs) =
!     >   ((0.25)*pi)*ehlw(mgs)*cx(mgs,lc)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))
!     >  *(  gf1*xdia(mgs,lc,2)
!     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lhl,1)
!     >    + gf3*xdia(mgs,lhl,2) )
!      chlacw(mgs) = min(chlacw(mgs),0.5*cx(mgs,lc)*dtpinv)
!      chlacw(mgs) = min(chlacw(mgs),cxmxd(mgs,lc))
!      chlacw(mgs) = min(chlacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      chlaci(:) = 0.0
      chlaci0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( lhl .gt. 1 .and. ( ehli(mgs) .gt. 0.0 .or. (ipelec > 0 .and. ehliclsn(mgs) > 0.0) )  ) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,li,1) )

          chlaci0(mgs) = 0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*cx(mgs,li)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0(lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            da0(li)*xdia(mgs,li,3)**2 )

!       ELSE
!        chlaci(mgs) =
!     >   ((0.25)*pi)*ehli(mgs)*cx(mgs,li)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))
!     >  *(  gf1*xdia(mgs,li,2)
!     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lhl,1)
!     >    + gf3*xdia(mgs,lhl,2) )
        ENDIF

        chlaci(mgs) = min(ehli(mgs)*chlaci0(mgs),cimxd(mgs))
       ENDIF
      end do
      end if


!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22jj'
      chlacs(:) = 0.0
      chlacs0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      IF ( lhl .gt. 1 .and. ( ehls(mgs) .gt. 0.0 .or. (ipelec > 0 .and. ehlsclsn(mgs) > 0.0) ) ) THEN
       IF ( ipconc .ge. 5 ) THEN

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,ls,1) )

          chlacs0(mgs) = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*cx(mgs,ls)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0(lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            da0(ls)*xdia(mgs,ls,3)**2 )

!       ELSE
!      chlacs(mgs) =
!     >   ((0.25)*pi)*ehls(mgs)*cx(mgs,ls)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))
!     >  *(  gf3*gf1*xdia(mgs,ls,2)
!     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lhl,1)
!     >    + gf1*gf3*xdia(mgs,lhl,2) )
      ENDIF
      chlacs(mgs) = min(ehls(mgs)*chlacs0(mgs),csmxd(mgs))
      ENDIF
      end do
      end if

!
! Ziegler (1985) autoconversion
!
!
      IF ( ipconc .ge. 2 ) THEN
      if (ndebug .gt. 0 ) write(0,*) 'conc 26a'
      
      DO mgs = 1,ngscnt
        zrcnw(mgs) = 0.0
        qrcnw(mgs) = 0.0
        crcnw(mgs) = 0.0
        cautn(mgs) = 0.0
      ENDDO
      
      IF ( dmrauto >= -1 ) THEN !{
      DO mgs = 1,ngscnt
!      qracw(mgs) = 0.0
!      cracw(mgs) = 0.0
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1000. .and. temg(mgs) .gt. tfrh+4.) THEN
       !( .and. w(igs(mgs),jgs,kgs(mgs)) > 5.0) THEN ! DTD: added w threshold for testing                                                                                                            
         volb = xv(mgs,lc)*(1./(1.+alpha(mgs,lc)))**(1./2.)
         cautn(mgs) = Min(ccmxd(mgs),   &
     &      ((alpha(mgs,lc)+2.)/(alpha(mgs,lc)+1.))*aa1*cx(mgs,lc)**2*xv(mgs,lc)**2)
         cautn(mgs) = Max( 0.0d0, cautn(mgs) )
         IF ( rb(mgs) .le. 7.51d-6 .or. dmrauto == -1) THEN
           t2s = 1.d30
!           cautn(mgs) = 0.0
         ELSE
!         XL2P=2.7E-2*XNC*XVC*((1.E12*RB**3*RC)-0.4)
         
!        T2S=3.72E-3/(((1.E4*RB)-7.5)*XNC*XVC) 
!           t2s = 3.72E-3/(((1.e6*rb)-7.5)*cx(mgs,lc)*xv(mgs,lc))
!           t2s = 3.72/(((1.e6*rb(mgs))-7.5)*rho0(mgs)*qx(mgs,lc))
           t2s = 3.72/(1.e6*(rb(mgs)-7.500d-6)*rho0(mgs)*qx(mgs,lc))

           qrcnw(mgs) = Max( 0.0d0, xl2p(mgs)/(t2s*rho0(mgs)) )
           crcnw(mgs) = Max( 0.0d0, Min(3.5e9*xl2p(mgs)/t2s,0.5*cautn(mgs)) )
           
           IF ( dmrauto == 0 ) THEN
             IF ( qx(mgs,lr)*rho0(mgs) > 1.2*xl2p(mgs) .and. cx(mgs,lr) > cxmin ) THEN ! Cohard and Pinty (2000a) switch over from (18) to (19)
               crcnw(mgs) = cx(mgs,lr)/qx(mgs,lr)*qrcnw(mgs)
             ELSEIF ( ( dmropt == 1 .or. dmropt == 3 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               crcnw(mgs) = Min(tmp,crcnw(mgs) )
             ELSEIF ( ( dmropt == 4 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr using converted qc mass
               crcnw(mgs) = (tmp*qrcnw(mgs)+tmp2*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
             ELSEIF ( ( dmropt == 5 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr using full qc mass
               crcnw(mgs) = (tmp*qx(mgs,lc)+tmp2*qx(mgs,lr))/(qx(mgs,lc)+qx(mgs,lr))
             ELSEIF ( ( dmropt == 6 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass*diameter-weighted average of old and new Dmr (using full qc mass)
               crcnw(mgs) = (tmp*xdia(mgs,lc,3)*qx(mgs,lc)+tmp2*xdia(mgs,lr,3)*qx(mgs,lr))/ &
                             (xdia(mgs,lc,3)*qx(mgs,lc)+xdia(mgs,lr,3)*qx(mgs,lr))
             ELSEIF ( ( dmropt == 7 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try diameter-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*xdia(mgs,lc,3)+tmp2*xdia(mgs,lr,3))/(xdia(mgs,lc,3)+xdia(mgs,lr,3))
             ELSEIF ( ( dmropt == 8 ) .and. qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try sqrt(diameter)-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*sqrt(xdia(mgs,lc,3))+tmp2*sqrt(xdia(mgs,lr,3)))/  &
                             (sqrt(xdia(mgs,lc,3))+sqrt(xdia(mgs,lr,3)))
             ENDIF
           ELSEIF ( dmrauto == 1  .and. cx(mgs,lr) > cxmin) THEN
             IF ( qx(mgs,lr) > qxmin(lr) ) THEN
               tmp = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               crcnw(mgs) = Min(tmp,crcnw(mgs) )
             ENDIF
           ELSEIF ( dmrauto == 2  .and. cx(mgs,lr) > cxmin) THEN
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*qrcnw(mgs)+tmp2*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
           ELSEIF ( dmrauto == 3  .and. cx(mgs,lr) > cxmin) THEN ! adapted from MY/CP code
              tmp = Max( 2.d0*rh(mgs), dble( xdia(mgs,lr,3) ) )
              crcnw(mgs) = rho0(mgs)*qrcnw(mgs)/(pi/6.*1000.*tmp**3)
           ENDIF
           
           IF ( crcnw(mgs) < 1.e-30 ) qrcnw(mgs) = 0.0

           IF ( ipconc >= 6 ) THEN
           IF ( lzr > 1 .and. qrcnw(mgs) > 0.0 ) THEN
!            vr = rho0(mgs)*qrcnw(mgs)/(1000.*crcnw(mgs))
!            zrcnw(mgs) = 36.*(xnu(lr)+2.0)*crcnw(mgs)*vr**2/((xnu(lr)+1.0)*pi**2)
             ! DTD: If rain exists at a grid point already either use the alpha-preserving Z-rate eqn. (dmrauto == 1)
             ! or a mass-weighted average of the alpha-preserving Z-rate eqn. and the init. rate eqn. (dmrauto == 2)
             ! or the original initiation rate equation (dmrauto == 0).  Not sure if this is the correct way to go but seems to work ok.
             IF (qx(mgs,lr) .gt. qxmin(lr) .and. ( dmrauto == 1 .or. dmrauto ==2 ) ) THEN
              tmp3 = qx(mgs,lr)/cx(mgs,lr)
              tmp4 =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &                 ( 2.*tmp3 * qrcnw(mgs) - tmp3**2 * crcnw(mgs)  )
              if (imurain == 3) then
                vr = rho0(mgs)*qrcnw(mgs)/(1000.)
                tmp3 = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              else
                tmp3 = galpharaut*(6.*rho0(mgs)*qrcnw(mgs)/(pi*xdn0(lr)))**2/crcnw(mgs)
              endif
              IF ( dmrauto == 1 ) THEN ! Preserve alpha
                zrcnw(mgs) = tmp4
              ELSEIF ( dmrauto == 2 ) THEN ! Mass-weighted average
                zrcnw(mgs) = (tmp3*qrcnw(mgs)+tmp4*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
              ENDIF
             else ! original formulation
              IF ( imurain == 3 ) THEN
                vr = rho0(mgs)*qrcnw(mgs)/(1000.) ! crcnw(mgs) not divided here but is in next line, cancels one factor in the numerator
                zrcnw(mgs) = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              ELSE ! rain in gamma of diameter
                IF ( dmropt <= 1 .or. dmropt >= 4 .or. ( qx(mgs,lr) < qxmin(lr) .and. cx(mgs,lr) < cxmin ) ) THEN
                  zrcnw(mgs) = galpharaut*(6.*rho0(mgs)*qrcnw(mgs)/(pi*xdn0(lr)))**2/crcnw(mgs)
                ELSE
                  tmp3 = qx(mgs,lr)/cx(mgs,lr)
                  zrcnw(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &                 ( 2.*tmp3 * qrcnw(mgs) - tmp3**2 * crcnw(mgs)  )
                ENDIF
!             vr = rho0(mgs)*qrcnw(mgs)/(1000.) ! crcnw(mgs) not divided here but is in next line, cancels one factor in the numerator
!             zrcnw(mgs) = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              ENDIF
             endif
!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
           ENDIF 
           ENDIF ! ipconc >= 6
!           IF (  crcnw(mgs) .gt. cautn(mgs) .and. crcnw(mgs) .gt. 1.0 )
!     :          THEN
!             write(0,*)  'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s,qx(mgs,lr)
!             write(0,*)  '            ',qx(mgs,lc),cx(mgs,lc),0.5e6*xdia(mgs,lc,1)
!             write(0,*)  '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :         1.e6*(( 3/(4.*pi))*rho0(mgs)*qrcnw(mgs)/
!     :       (crcnw(mgs)*xdn(mgs,lr)))**(1./3.),rh(mgs)*1.e6,rwrad(mgs)
!           ELSEIF ( crcnw(mgs) .gt. 1.0 .and. cautn(mgs) .gt. 0.) THEN
!             write(0,*)  'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s
!             write(0,*)  '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :  1.e6*(( 3*pi/4.)*rho0(mgs)*qrcnw(mgs)/
!     :   (crcnw(mgs)*xdn(mgs,lr)))**(1./3.)
!           ENDIF
!           crcnw(mgs) = Min(cautn(mgs),3.5e9*xl2p(mgs)/t2s)

!           IF ( qrcnw(mgs) .gt. 0.3e-2 ) THEN
!            write(0,*)  'QRCNW'
!            write(0,*)  qrcnw(mgs),crcnw(mgs),cautn(mgs)
!            write(0,*)  xl2p,t2s,rho0(mgs),xv(mgs,lc),cx(mgs,lc),qx(mgs,lc)
!            write(0,*)  rb,0.5*xdia(mgs,lc,1),mgs,igs(mgs),kgs(mgs)
!           ENDIF
!           qrcnw(mgs) = Min(qrcnw(mgs),qcmxd(mgs))
         ENDIF


       ENDIF
      ENDDO
      
      ENDIF !} dmrauto >= 0



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
!
!
!  Bigg Freezing of Rain
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 27a'
      qrfrz(:) = 0.0
      qrfrzs(:) = 0.0
      qrfrzf(:) = 0.0
      vrfrzf(:) = 0.0
      crfrz(:) = 0.0
      crfrzs(:) = 0.0
      crfrzf(:) = 0.0
      zrfrz(:)  = 0.0
      zrfrzs(:)  = 0.0
      zrfrzf(:)  = 0.0
      qwcnr(:) = 0.0
      
      IF ( .not. ( ipconc == 0 .and. lwsm6 ) ) THEN
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lr) .gt. qxmin(lr) .and. temcg(mgs) .lt. -5. .and. ibiggopt > 0 ) then
!      brz = 100.0
!      arz = 0.66
       IF ( ipconc .lt. 3 ) THEN
       qrfrz(mgs) =    &
     &  min(   &
     &  (20.0)*(pi**2)*brz*(xdn(mgs,lr)/rho0(mgs))   &
     &   *cx(mgs,lr)*(xdia(mgs,lr,1)**6)   &
     &   *(exp(max(-arz*temcg(mgs), 0.0))-1.0)   &
     &  , qrmxd(mgs))
        qrfrzf(mgs) = qrfrz(mgs)

!       ELSEIF ( ipconc .ge. 3 .and. xv(mgs,lr) .gt. 1.1*xvmn(lr) ) THEN
       ELSEIF ( ipconc .ge. 3 ) THEN
!         tmp = brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
!         crfrz(mgs) = xv(mgs,lr)*tmp

         frach = 1.0d0
         
!         IF ( ibiggopt == 2 .and. imurain == 1 .and. lzr < 1 ) THEN ! lzr check because results are weird for 3-moment
         IF ( ibiggopt == 2 .and. imurain == 1 ) THEN !
         ! integrate from Bigg diameter (for given supercooling Ts) to infinity
           
           volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 (Proc. Phys. Soc. London) 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b, where a = 6.9/6.8, or approx a = 1.0, and b = 16.2
                                               ! volt is given in cm**3, so convert to m**3
           dbigg = (6./pi* volt )**(1./3.) 
           
           ! perhaps should also test that W > V_t_dbigg, i.e., that drops the size of dbigg are being lifted and cooled. 
           IF ( dbigg < 8.e-3 ) THEN !{ only bother if freezing diameter is reasonable
           
             ratio = Min(maxratiolu, dbigg/xdia(mgs,lr,1) )
           
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           crfrz(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)*dtpinv
           crfrzf(mgs) = crfrz(mgs)
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qrfrz(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)*dtpinv
           qrfrzf(mgs) = qrfrz(mgs)

           IF ( qrfrz(mgs)*dtp < qxmin(lh) .or. crfrz(mgs)*dtp < cxmin ) THEN
           
             crfrz(mgs) = 0.0
             qrfrz(mgs) = 0.0
             qrfrzf(mgs) = 0.0
            
           ELSE !{


           IF ( (ipconc >= 5 .or. lzr > 1) ) THEN !{

              cxd1 = crfrz(mgs)*dtp
              qxd1 = qrfrz(mgs)*dtp

           ! interpolate along x, i.e., ratio; 
            tmp1 = ziacrratio(i,j) + delx*dqiacrratioinv*(ziacrratio(ip1,j) - ziacrratio(i,j))
            tmp2 = ziacrratio(i,jp1) + delx*dqiacrratioinv*(ziacrratio(ip1,jp1) - ziacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
            IF ( ipconc >= 6 .and. lzr > 1  ) THEN
              zxd1 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*zx(mgs,lr)
              ! Do the correction for alphamax
              zrfrz(mgs) = zxd1*dtpinv
              ! tmp4 is the Z from the converted particles assuming shape of alphamax
              IF ( icorrectfddbz >= 1 .and. zxd1 > 10.*zxmin ) THEN
              tmp3 = g1xmax*(rho0(mgs)*qxd1)**2/((pi*rhofrz/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new graupel/fd number to match zxd1
                ! increase cxd1 to make z,q,c rates consistent
                ! cxd1 = g1xmax*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lh)/6.0)**2)
                cxd1 = tmp3/zxd1
                crfrzf(mgs) = dtpinv*cxd1
              ENDIF
              ENDIF
            ELSE
              IF ( icorrectfddbz >= 1 ) THEN
            ! tmp5 is rain reflectivity moment
              tmp5 = g1x(mgs,lr)*(rho0(mgs)*qx(mgs,lr))**2/((pi*xdn(mgs,lr)/6.)**2*cx(mgs,lr))
              zxd1 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*tmp5
            ! tmp4 is the reflectivity of the newly-converted graupel particles (use g1x(lh) for loss term)
            ! which we want to match zxd1 to prevent spurious increase in total reflectivity
              IF ( zxd1 > 10.*zxmin ) THEN
              tmp3 =  g1x(mgs,lr)*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lr)/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new FD number to match zxd1
                crfrzf(mgs) = tmp3/zxd1*dtpinv
              ENDIF
              ENDIF
              ENDIF
            ENDIF
           ENDIF !}

           
            IF ( ibiggsmallrain > 0 .and. xv(mgs,lr) < 2.*xvmn(lr) .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) THEN
!            IF ( ibiggsmallrain > 0 .and. xv(mgs,lr) < xvbiggsnow .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) THEN
             ! rain drops are so small that they cannot be pushed smaller, so put into snow (or cloud ice, depending on ifrzs)
              crfrzf(mgs) = 0.0
              qrfrzf(mgs) = 0.0
              crfrzs(mgs) = crfrz(mgs)
              qrfrzs(mgs) = qrfrz(mgs)

              IF ( ipconc >= 6 .and. lzr > 1 ) THEN
                zrfrzs(mgs) = zrfrz(mgs)
                zrfrzf(mgs) = 0.
              ENDIF
           ELSEIF ( dbigg < Max( biggsnowdiam, Max(dfrz,dhmn)) .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) THEN ! { convert some to snow or ice crystals
            ! temporarily store qrfrz and crfrz in snow terms and caclulate new crfrzf, qrfrzf, and zrfrzf. Leave crfrz etc. alone!
            
            crfrzs(mgs) = crfrz(mgs)
            qrfrzs(mgs) = qrfrz(mgs)
            
            IF ( ibiggsmallrain > 0 .and. xv(mgs,lr) < 1.2*xvmn(lr) ) THEN
             ! rain drops are so small that they cannot be pushed smaller, so put into snow (or cloud ice, depending on ifrzs)
            crfrzf(mgs) = 0.0
            qrfrzf(mgs) = 0.0

             IF (ipconc >= 6 .and. lzr > 1 ) THEN
               zrfrzs(mgs) = zrfrz(mgs)
               zrfrzf(mgs) = 0.
             ENDIF
            ELSE !{
            
           ! recalculate using dhmn for ratio
           ratio = Min( maxratiolu, Max(dfrz,dhmn)/xdia(mgs,lr,1) )
           
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
!           j = Int(Max(alphamin,Min(alphamax,alpha(mgs,lr)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(alphamax,alpha(mgs,lr)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))


           ! interpolate along alpha; 
           
           crfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)*dtpinv
           
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qrfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)*dtpinv

           ! now subtract off the difference
            crfrzs(mgs) = crfrzs(mgs) - crfrzf(mgs)
            qrfrzs(mgs) = qrfrzs(mgs) - qrfrzf(mgs)

           IF ( ipconc >= 6 .and. lzr > 1 ) THEN
            zrfrzs(mgs) = zrfrz(mgs)
           ! interpolate along x, i.e., ratio; 
            tmp1 = ziacrratio(i,j) + delx*dqiacrratioinv*(ziacrratio(ip1,j) - ziacrratio(i,j))
            tmp2 = ziacrratio(i,jp1) + delx*dqiacrratioinv*(ziacrratio(ip1,jp1) - ziacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
            zrfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*zx(mgs,lr)*dtpinv
            zrfrzs(mgs) = zrfrzs(mgs) - zrfrzf(mgs)
            zrfrzf(mgs) = (1000./900.)**2*zrfrzf(mgs)
           ENDIF

           IF ( ( ipconc >= 5 .or. lzr > 1 ) ) THEN !{

              cxd1 = crfrzf(mgs)*dtp
              qxd1 = qrfrzf(mgs)*dtp

           ! interpolate along x, i.e., ratio; 
            tmp1 = ziacrratio(i,j) + delx*dqiacrratioinv*(ziacrratio(ip1,j) - ziacrratio(i,j))
            tmp2 = ziacrratio(i,jp1) + delx*dqiacrratioinv*(ziacrratio(ip1,jp1) - ziacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
            IF ( ipconc >= 6 .and. lzr > 1  ) THEN !{
              zxd1 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*zx(mgs,lr)
              ! Do the correction for alphamax
              zrfrz(mgs) = zxd1*dtpinv
              ! tmp4 is the Z from the converted particles assuming shape of alphamax
              IF ( icorrectfddbz >= 2  .and. zxd1 > 10.*zxmin ) THEN
              tmp3 = g1xmax*(rho0(mgs)*qxd1)**2/((pi*rhofrz/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new graupel/fd number to match zxd1
                ! increase cxd1 to make z,q,c rates consistent
                ! cxd1 = g1xmax*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lh)/6.0)**2)
                cxd1 = tmp3/zxd1
                crfrzf(mgs) = dtpinv*cxd1
              ENDIF
              ENDIF
            ELSE ! }{
              IF ( icorrectfddbz >= 2  ) THEN
            ! tmp5 is rain reflectivity moment
              tmp5 = g1x(mgs,lr)*(rho0(mgs)*qx(mgs,lr))**2/((pi*xdn(mgs,lr)/6.)**2*cx(mgs,lr))
              zxd1 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*tmp5
              IF ( zxd1 > 10.*zxmin ) THEN
            ! tmp4 is the reflectivity of the newly-converted graupel particles (use g1x(lh) for loss term)
            ! which we want to match zxd1 to prevent spurious increase in total reflectivity
              tmp3 =  g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lr)/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new FD number to match zxd1
                crfrzf(mgs) = tmp3/zxd1*dtpinv
              ENDIF
              ENDIF
              ENDIF
            ENDIF !}
           ENDIF !}

            ENDIF ! }
           ELSE
            crfrzs(mgs) = 0.0
            qrfrzs(mgs) = 0.0
            zrfrzs(mgs) = 0.0
           ENDIF ! }
           
           ENDIF !}
           
           IF ( (qrfrz(mgs))*dtp > qx(mgs,lr) ) THEN
             fac = ( qrfrz(mgs) )*dtp/qx(mgs,lr)
             qrfrz(mgs) = fac*qrfrz(mgs)
             qrfrzs(mgs) = fac*qrfrzs(mgs)
             qrfrzf(mgs) = fac*qrfrzf(mgs)
             crfrz(mgs) = fac*crfrz(mgs)
             crfrzs(mgs) = fac*crfrzs(mgs)
             crfrzf(mgs) = fac*crfrzf(mgs)
             IF ( ipconc >= 6 .and. lzr > 1 ) THEN
               zrfrz(mgs) = fac*zrfrz(mgs)
               zrfrzf(mgs) = fac*zrfrzf(mgs)
             ENDIF
           ENDIF
           
            ENDIF !}

!           IF ( (crfrzs(mgs) + crfrz(mgs))*dtp > cx(mgs,lr) ) THEN
!             fac = ( crfrzs(mgs) + crfrz(mgs) )*dtp/cx(mgs,lr)
!             crfrz(mgs) = fac*crfrz(mgs)
!             crfrzs(mgs) = fac*crfrzs(mgs)
!           ENDIF
           
!           qrfrzf(mgs) = qrfrz(mgs)
!           crfrzf(mgs) = crfrz(mgs)
           
   !        qrfrz(mgs) = qrfrzf(mgs) + qrfrzs(mgs)
   !        crfrz(mgs) = crfrzf(mgs) + crfrzs(mgs)

           
         ELSEIF ( ibiggopt == 1 ) THEN
         ! Z85, eq. A34
         tmp = xv(mgs,lr)*brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
         IF ( .false. .and. tmp .gt. cxmxd(mgs,lr) ) THEN ! {
!           write(iunit,*) 'Bigg Freezing problem!',mgs,igs(mgs),kgs(mgs)
!           write(iunit,*)  'tmp, cx(lr), xv = ',tmp, cx(mgs,lr), xv(mgs,lr), (Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
!           write(iunit,*)  'qr,temcg = ',qx(mgs,lr)*1000.,temcg(mgs)
           crfrz(mgs) = cxmxd(mgs,lr) ! cx(mgs,lr)*dtpinv
           qrfrz(mgs) = qxmxd(mgs,lr) ! qx(mgs,lr)*dtpinv
!           STOP
         ELSE ! } {
         crfrz(mgs) = tmp
 !        crfrzfmx = cx(mgs,lr)*Exp(-4./3.*pi*(40.e-6)**3/xv(mgs,lr))
 !        IF ( crfrz(mgs) .gt. crfrzmx ) THEN
 !          crfrz(mgs) = crfrzmx
 !          qrfrz(mgs) = bfnu*xmas(mgs,lr)*rhoinv(mgs)*crfrzmx
 !          qwcnr(mgs) = cx(mgs,lr) - crfrzmx
 !        ELSE
         IF ( lzr < 1 ) THEN
           IF ( imurain == 3 ) THEN
             bfnu = bfnu0
           ELSE !imurain == 1
             bfnu = bfnu1
           ENDIF
         ELSE
 !         bfnu = 1.0 ! (alpha(mgs,lr)+2.0)/(alpha(mgs,lr)+1.)
           IF ( imurain == 3 ) THEN
             bfnu = (alpha(mgs,lr)+2.0)/(alpha(mgs,lr)+1.)
           ELSE !imurain == 1
!             bfnu = bfnu1
            bfnu = (4. + alpha(mgs,lr))*(5. + alpha(mgs,lr))*(6. + alpha(mgs,lr))/  &
     &            ((1. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(3. + alpha(mgs,lr)))
!            bfnu = 1.
           ENDIF
         ENDIF 
         qrfrz(mgs) = bfnu*xmas(mgs,lr)*rhoinv(mgs)*crfrz(mgs)

         qrfrz(mgs) = Min( qrfrz(mgs), 1.*qx(mgs,lr)*dtpinv ) ! qxmxd(mgs,lr) 
         crfrz(mgs) = Min( crfrz(mgs), 1.*cx(mgs,lr)*dtpinv ) !cxmxd(mgs,lr) 
         qrfrz(mgs) = Min( qrfrz(mgs), qx(mgs,lr) )
         qrfrzf(mgs) = qrfrz(mgs)
         ENDIF !}

         
         
         
         IF ( crfrz(mgs) .gt. qxmin(lh) ) THEN !{ Yes, it compares cx and qxmin, but this is just to be sure that 
                                                  ! crfrz is greater than zero in the division
!          IF ( xdia(mgs,lr,1) .lt. 200.e-6 ) THEN
!           IF ( xv(mgs,lr) .lt. xvmn(lh) ) THEN
           
           IF ( (ibiggsnow == 1 .or. ibiggsnow == 3 ) .and. ibiggopt /= 2 ) THEN
           xvfrz = rho0(mgs)*qrfrz(mgs)/(crfrz(mgs)*900.) ! mean volume of frozen drops; 900. for frozen drop density
           frach = 0.5 *(1. +  Tanh(0.2e12 *( xvfrz - 1.15*xvmn(lh))))

             qrfrzs(mgs) = (1.-frach)*qrfrz(mgs)
             crfrzs(mgs) = (1.-frach)*crfrz(mgs) ! *rzxh(mgs)
!             qrfrzf(mgs) = frach*qrfrz(mgs)
           
           ENDIF
           
           IF ( ipconc .ge. 14 .and. 1.e-3*rho0(mgs)*qrfrz(mgs)/crfrz(mgs) .lt. xvmn(lh) ) THEN
             qrfrzs(mgs) = qrfrz(mgs)
             crfrzs(mgs) = crfrz(mgs) ! *rzxh(mgs)
           ELSE
!           crfrz(mgs) = Min( crfrz(mgs), 0.1*cx(mgs,lr)*dtpinv ) ! cxmxd(mgs,lr) 
!           qrfrz(mgs) = Min( qrfrz(mgs), 0.1*qx(mgs,lr)*dtpinv ) ! qxmxd(mgs,lr) 
             qrfrzf(mgs) = frach*qrfrz(mgs)
!             crfrzf(mgs) = Min( qrfrz(mgs)*rho0(mgs)/(xdn(mgs,lh)*vgra), crfrz(mgs) )
            IF ( ibfr .le. 1 ) THEN
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*1.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            ELSEIF ( ibfr .eq. 5 ) THEN
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*vfrz*1000.0)*rho0(mgs) )*rzxh(mgs)  !*crfrz(mgs)
            ELSEIF ( ibfr .eq. 2 ) THEN
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*vfrz*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            ELSEIF ( ibfr .eq. 6 ) THEN
             crfrzf(mgs) = frach*Max(crfrz(mgs), qrfrz(mgs)/(bfnu*9.*xv(mgs,lr)*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            ELSE
             crfrzf(mgs) = frach*crfrz(mgs)
            ENDIF 
!             crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*xvmn(lh)*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
!            IF ( lz(lr) > 1 .and. lz(lh) > 1 ) THEN
!              crfrzf(mgs) = crfrz(mgs)
!            ENDIF
            
           ENDIF
!         crfrz(mgs) = Min( cxmxd(mgs,lr), rho0(mgs)*qrfrz(mgs)/xmas(mgs,lr) )
         ELSE
          crfrz(mgs) = 0.0
          qrfrz(mgs) = 0.0
         ENDIF !}

         ENDIF ! ibiggopt

          IF ( lvol(lh) .gt. 1 ) THEN
           vrfrzf(mgs) = rho0(mgs)*qrfrzf(mgs)/rhofrz
          ENDIF

        
        IF ( nsplinter .ne. 0 ) THEN
          IF ( nsplinter .ge. 1000 ) THEN
           ! Lawson et al. 2015 JAS
           ! ave. diam of freezing drops in microns
            tmp = 0
            IF ( qrfrz(mgs)*dtp > qxmin(lh) .and. crfrz(mgs) > 1.e-3 ) THEN
              tmpdiam = 1.e6*( 6.*qrfrz(mgs)/(1000.*pi*crfrz(mgs) ))**(1./3.)  ! avg. diameter of newly frozen drops in microns
              fac = 1.0
              IF ( nsplinter .eq. 1001 ) THEN
             !   fac = 0.2/sqrt(2.0*pi*10.**2)*Exp(-0.5*((258.-temg(mgs))/10.)**2 ) ! temperature dependence from Sullivan et al. 2018 ACP
             ! ELSE
                fac = 0.2*Exp(-0.5*((258.-temg(mgs))/10.)**2 ) ! temperature dependence from Sullivan et al. 2018 ACP
              ENDIF
              tmp = fac*lawson_splinter_fac*tmpdiam**4*crfrz(mgs)
            ENDIF
          ELSEIF ( nsplinter .gt. 0 ) THEN
            tmp = nsplinter*crfrz(mgs)
          ELSE
            tmp = -nsplinter*crfrzf(mgs)
          ENDIF
          csplinter2(mgs) = tmp
          qsplinter2(mgs) = Min(0.1*qrfrz(mgs), tmp*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel

!          csplinter(mgs) = csplinter(mgs) + tmp
!          qsplinter(mgs) = qsplinter(mgs) + Min(0.1*qrfrz(mgs), tmp*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel
        ENDIF
!         IF ( temcg(mgs) .lt. -31.0 ) THEN
!           qrfrz(mgs) = qx(mgs,lr)*dtpinv + qrcnw(mgs)
!           qrfrzf(mgs) = qrfrz(mgs)
!           crfrz(mgs) = cx(mgs,lr)*dtpinv + crcnw(mgs)
!           crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*1.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
!         ENDIF
!         qrfrz(mgs) = 6.0*xdn(mgs,lr)*xv(mgs,lr)**2*tmp*rhoinv(mgs)
!         qrfrz(mgs) = Min( qrfrz(mgs), ffrz*qrmxd(mgs) )
!         crfrz(mgs) = Min( crmxd(mgs), ffrz*crfrz(mgs))
!         crfrz(mgs) = Min(crmxd(mgs),qrfrz(mgs)*rho0(mgs)/xmas(mgs,lr))
       ENDIF
!      if ( temg(mgs) .gt. 268.15 ) then
      else
!      end if
      end if
      end do
      
      ENDIF
!
!  Homogeneous freezing of cloud drops to ice crystals
!  following Bigg (1953) and Ferrier (1994).
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 25b'
      do mgs = 1,ngscnt
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      qwfrzc(mgs) = 0.0
      cwfrzc(mgs) = 0.0
      qwfrzp(mgs) = 0.0
      cwfrzp(mgs) = 0.0
      IF ( ibfc .ge. 1 .and. ibfc /= 3 .and. temg(mgs) < 268.15 ) THEN
!      if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1.  .and.   &
!     &     .not. (ipconc .ge. 2 .and. xdia(mgs,lc,1) .lt. 10.e-6) ) then
      if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. cxmin ) THEN
      IF ( ipconc < 2 ) THEN
      qwfrz(mgs) = ((2.0)*(brz)/(xdn(mgs,lc)*cx(mgs,lc)))   &
     &  *(exp(max(-arz*temcg(mgs), 0.0))-1.0)   &
     &  *rho0(mgs)*(qx(mgs,lc)**2)
      qwfrz(mgs) = max(qwfrz(mgs), 0.0)
      qwfrz(mgs) = min(qwfrz(mgs),qcmxd(mgs))
         cwfrz(mgs) = qwfrz(mgs)*rho0(mgs)/xmas(mgs,li)
       ELSEIF ( ipconc .ge. 2 ) THEN
         IF ( xdia(mgs,lc,3) > 0.e-6 ) THEN
          volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b
                                               ! volt is given in cm**3, so factor of 1.e-6 to convert to m**3
!           dbigg = (6./pi* volt )**(1./3.) 

         IF (  alpha(mgs,lc) == 0.0 ) THEN
         cwfrz(mgs) = cx(mgs,lc)*Exp(-volt/xv(mgs,lc))*dtpinv ! number of droplets with volume greater than volt
!turn off limit so that all can freeze at low temp
!!!       cwfrz(mgs) = Min(cwfrz(mgs),ccmxd(mgs))

         qwfrz(mgs) = cwfrz(mgs)*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
          ELSE
            ratio = (1. + alpha(mgs,lc))*volt/xv(mgs,lc)
            
            IF ( .false. .and. usegamxinfcnu ) THEN
              i = Nint(dgami*(1. + alpha(mgs,lc)))
              gcnup1 = gmoi(i)
              i = Nint(dgami*(2. + alpha(mgs,lc)))
              gcnup2 = gmoi(i)

              cwfrz(mgs) = cx(mgs,lc)*Gamxinf(1.+alpha(mgs,lc), ratio)/(dtp*gcnup1) ! gamxinflu(i,j,1,1)

              qwfrz(mgs) = cx(mgs,lc)*xdn0(lc)*xv(mgs,lc)*rhoinv(mgs)*Gamxinf(2.+alpha(mgs,lc), ratio)/(dtp*gcnup2) !  gamxinflu(i,j,12,1)
            
            ELSE
            
              ratio = Min( maxratiolu, ratio )
!              write(0,*) 'cwfrz: temp,ratio = ',temcg(mgs),ratio
!              write(0,*) 'cwfrz: xv,volt,qx = ',xv(mgs,lc),volt,qx(mgs,lc)
!              write(0,*) 'cwfrz: i,j,k = ',igs(mgs),jgs,kgs(mgs)
              tmp = gaminterp(ratio,alpha(mgs,lc),1,1)
!              write(0,*) 'cwfrz: tmp1 = ',tmp
              cwfrz(mgs) = cx(mgs,lc)*tmp*dtpinv ! Gamxinf(1.+alpha(mgs,lc), ratio)/(dtp*gcnup1) ! gamxinflu(i,j,1,1)

              tmp = gaminterp(ratio,alpha(mgs,lc),12,1)
!              write(0,*) 'cwfrz: tmp2 = ',tmp
              qwfrz(mgs) = cx(mgs,lc)*xdn0(lc)*xv(mgs,lc)*rhoinv(mgs)*dtpinv*tmp ! Gamxinf(2.+alpha(mgs,lc), ratio)/(dtp*gcnup2) !  gamxinflu(i,j,12,1)
            
            ENDIF
          
          ENDIF

         ENDIF
       ENDIF
      if ( temg(mgs) .gt. 268.15 ) then
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      end if
      end if
      ENDIF
!
        if ( xplate(mgs) .eq. 1 ) then
          qwfrzp(mgs) = qwfrz(mgs)
          cwfrzp(mgs) = cwfrz(mgs)
        end if
!  
        if ( xcolmn(mgs) .eq. 1 ) then
          qwfrzc(mgs) = qwfrz(mgs)
          cwfrzc(mgs) = cwfrz(mgs)
        end if
      
!
!     qwfrzp(mgs) = 0.0
!     qwfrzc(mgs) = qwfrz(mgs)
!
      end do
!
!
!  Contact freezing nucleation:  factor is to convert from L-1
!  T < -2C:  via Meyers et al. JAM July, 1992 (31, 708-721)
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 25a'
      do mgs = 1,ngscnt

       ccia(mgs) = 0.0

       cwctfz(mgs) = 0.0
       qwctfz(mgs) = 0.0
       ctfzbd(mgs) = 0.0
       ctfzth(mgs) = 0.0
       ctfzdi(mgs) = 0.0

       cwctfzc(mgs) = 0.0
       qwctfzc(mgs) = 0.0
       cwctfzp(mgs) = 0.0
       qwctfzp(mgs) = 0.0
       IF ( icfn .ge. 1 ) THEN

       IF ( temg(mgs) .lt. 271.15  .and. qx(mgs,lc) .gt. qxmin(lc)) THEN

!       find available # of ice nuclei & limit value to max depletion of cloud water

        IF ( icfn .ge. 2 ) THEN
         ccia(mgs) = exp( 4.11 - (0.262)*temcg(mgs) )  ! in m-3, see Walko et al. 1995; 1000*exp(-2.8 -b*t) = exp(6.91)*exp(-2.8 - b*t) = exp(4.11 -b*t)
         !ccia(mgs) = Min(cwctfz(mgs), ccmxd(mgs) )

!       now find how many of these collect cloud water to form IN
!       Cotton et al 1986

         knud(mgs) = 2.28e-5 * temg(mgs) / ( pres(mgs)*raero ) !Walko et al. 1995
         knuda(mgs) = 1.257 + 0.4*exp(-1.1/knud(mgs))          !Pruppacher & Klett 1997 eqn 11-16
         gtp(mgs) = 1. / ( fai(mgs) + fbi(mgs) )               !Byers 65 / Cotton 72b
         dfar(mgs) = kb*temg(mgs)*(1.+knuda(mgs)*knud(mgs))/(6.*pi*fadvisc(mgs)*raero) !P&K 1997 eqn 11-15
         fn1(mgs) = 2.*pi*xdia(mgs,lc,1)*cx(mgs,lc)*ccia(mgs)
         fn2(mgs) = -gtp(mgs)*(ssw(mgs)-1.)*felv(mgs)/pres(mgs)
         fnft(mgs) = 0.4*(1.+1.45*knud(mgs)+0.4*knud(mgs)*exp(-1./knud(mgs)))*(ftka(mgs)+2.5*knud(mgs)*kaero)      &
     &              / ( (1.+3.*knud(mgs))*(2*ftka(mgs)+5.*knud(mgs)*kaero+kaero) )


!      Brownian diffusion
         ctfzbd(mgs) = fn1(mgs)*dfar(mgs)

!      Thermophoretic contact nucleation
         ctfzth(mgs) = fn1(mgs)*fn2(mgs)*fnft(mgs)/rho0(mgs)

!      Diffusiophoretic contact nucleation
         ctfzdi(mgs) = fn1(mgs)*fn2(mgs)*rw*temg(mgs)/(felv(mgs)*rho0(mgs))

         cwctfz(mgs) = max( ctfzbd(mgs) + ctfzth(mgs) + ctfzdi(mgs) , 0.)

!      Sum of the contact nucleation processes
!         IF ( cx(mgs,lc) .gt. 1.e6) write(0,*) 'ctfzbd,etc = ',cwctfz(mgs),ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs)
!         IF ( wvel(mgs) .lt. -0.05 ) write(6,*) 'ctfzbd,etc = ',ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs),cx(mgs,lc)*1e-6,wvel(mgs)
!         IF ( ssw(mgs) .lt. 1.0 .and. cx(mgs,lc) .gt. 1.e6 .and. cwctfz(mgs) .gt. 1. ) THEN
!          write(6,*) 'ctfzbd,etc = ',ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs),cx(mgs,lc)*1e-6,wvel(mgs),fn1(mgs),fn2(mgs)
!          write(6,*) 'more = ',nstep,ssw(mgs),dfar(mgs),gtp(mgs),felv(mgs),pres(mgs)
!         ENDIF

        ELSEIF ( icfn .eq. 1 ) THEN
         IF ( wvel(mgs) .lt. -0.05 ) THEN ! older kludgy version
           cwctfz(mgs) = cfnfac*exp( (-2.80) - (0.262)*temcg(mgs) )
           cwctfz(mgs) = Min((1.0e3)*cwctfz(mgs), ccmxd(mgs) )  !convert to m-3
         ENDIF
        ENDIF   ! icfn

        IF ( ipconc .ge. 2 ) THEN
         cwctfz(mgs) = Min( cwctfz(mgs)*dtpinv, ccmxd(mgs) )
         qwctfz(mgs) = xmas(mgs,lc)*cwctfz(mgs)/rho0(mgs)
        ELSE
         qwctfz(mgs) = (cimasn)*cwctfz(mgs)/(dtp*rho0(mgs))
         qwctfz(mgs) = max(qwctfz(mgs), 0.0)
         qwctfz(mgs) = min(qwctfz(mgs),qcmxd(mgs))
        ENDIF

!
        if ( xplate(mgs) .eq. 1 ) then
         qwctfzp(mgs) = qwctfz(mgs)
         cwctfzp(mgs) = cwctfz(mgs)
        end if
!
        if ( xcolmn(mgs) .eq. 1 ) then
         qwctfzc(mgs) = qwctfz(mgs)
         cwctfzc(mgs) = cwctfz(mgs)
        end if
        
!        IF ( cwctfz(mgs)*dtp > 0.5 .and. dtp*qwctfz(mgs) > qxmin(li) ) THEN
!          write(91,*) 'cwctfz: ',cwctfz(mgs),qwctfz(mgs) ! ,cwctfzc(mgs),qwctfzc(mgs)
!        ENDIF
       
!
!     qwctfzc(mgs) = qwctfz(mgs)
!     qwctfzp(mgs) = 0.0
!
       end if

       ENDIF ! icfn

      end do
!
!
!
! Hobbs-Rangno ice enhancement (Ferrier, 1994)
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 23a'
      do mgs = 1,ngscnt
      ciihr(mgs) = 0.0
      qiihr(mgs) = 0.0
      cicichr(mgs) = 0.0
      qicichr(mgs) = 0.0
      cipiphr(mgs) = 0.0
      qipiphr(mgs) = 0.0
      ENDDO

      dthr = 300.0
     ! hrifac = (1.e-3)*((0.044)*(0.01**3))
      hrifac = cimas1
      IF ( ihrn .ge. 1 ) THEN
      do mgs = 1,ngscnt
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
      if ( temg(mgs) .lt. 265.15 ) then
!      write(iunit,'(3(1x,i3),3(1x,1pe12.5))')
!     : igs(mgs),jgs,kgs(mgs),cx(mgs,lc),rho0(mgs),qx(mgs,lc)
!      write(iunit,'(1pe15.6)')
!     :  log(cx(mgs,lc)*(1.e-6)/(3.0)),
!     :  ((1.e-3)*rho0(mgs)*qx(mgs,lc)),
!     :  (cx(mgs,lc)*(1.e-6)),
!     : ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)),
!     : (alog(cx(mgs,lc)*(1.e-6)/(3.0)) *
!     >  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)))

      IF ( Log(cx(mgs,lc)*(1.e-6)/(3.0)) .gt. 0.0 ) THEN
      ciihr(mgs) = ((1.69e17))   &
     & *(log(cx(mgs,lc)*(1.e-6)/(3.0)) *   &
     &  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)))**(7./3.)
      ciihr(mgs) = (ciihr(mgs)*(1.0e6) - cx(mgs,li) - cx(mgs,ls))/dthr
      qiihr(mgs) = hrifac*ciihr(mgs)/rho0(mgs)
      qiihr(mgs) = max(qiihr(mgs), 0.0)
      qiihr(mgs) = min(qiihr(mgs),qcmxd(mgs))
      ENDIF
!
      if ( xplate(mgs) .eq. 1 ) then
      qipiphr(mgs) = qiihr(mgs)
      cipiphr(mgs) = ciihr(mgs)
      end if
!
      if ( xcolmn(mgs) .eq. 1 ) then
      qicichr(mgs) = qiihr(mgs)
      cicichr(mgs) = ciihr(mgs)
      end if
!
!     qipiphr(mgs) = 0.0
!     qicichr(mgs) = qiihr(mgs)
!
      end if
      end if
      end do
      ENDIF ! ihrn
!
!
!
!  simple frozen rain to hail conversion.  All of the
!  frozen rain larger than 5.0e-3 m in diameter are converted
!  to hail.  This is done by considering the equation for
!  frozen rain mixing ratio:
!
!
!  qfw = [ cno(lf) * pi * fwdn / (6 rhoair) ]
!
!         /inf
!      *  |     fwdia*3 exp(-dia/fwdia) d(dia)
!         /Do
!
!  The amount to be reclassified as hail is the integral above from
!  Do to inf where Do is 5.0e-3 m.
!
!
!  qfauh = [ cno(lf) * pi * fwdn / (6 rhoair) ]
!
!


      hdia0 = 300.0e-6
      do mgs = 1,ngscnt
      qscnvi(mgs) = 0.0
      cscnvi(mgs) = 0.0
      cscnvis(mgs) = 0.0
!      IF ( .false. ) THEN
!      IF ( temg(mgs) .lt. tfr .and. ssi(mgs) .gt. 1.01 .and. qx(mgs,li) .gt. qxmin(li) ) THEN
      IF ( temg(mgs) .lt. tfr .and. qx(mgs,li) .gt. qxmin(li) ) THEN
        IF ( ipconc .ge. 4 .and. .false. ) THEN
         if ( cx(mgs,li) .gt. 10. .and. xdia(mgs,li,1) .gt. 50.e-6 ) then !{
         cirdiatmp =   &
     &  (qx(mgs,li)*rho0(mgs)   &
     & /(pi*xdn(mgs,li)*cx(mgs,li)))**(1./3.)
          IF ( cirdiatmp .gt. 100.e-6 ) THEN !{
          qscnvi(mgs) =   &
     &  ((pi*xdn(mgs,li)*cx(mgs,li)) / (6.0*rho0(mgs)*dtp))   &
     & *exp(-hdia0/cirdiatmp)   &
     & *( (hdia0**3) + 3.0*(hdia0**2)*cirdiatmp   &
     &  + 6.0*(hdia0)*(cirdiatmp**2) + 6.0*(cirdiatmp**3) )
      qscnvi(mgs) =   &
     &  min(qscnvi(mgs),qimxd(mgs))
          IF ( ipconc .ge. 4 ) THEN
            cscnvi(mgs) = Min( cimxd(mgs), cx(mgs,li)*Exp(-hdia0/cirdiatmp))
          ENDIF
         ENDIF  ! }
        end if ! }

       ELSEIF ( ipconc .lt. 4 ) THEN

        qscnvi(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
        qscnvi(mgs) = min(qscnvi(mgs),qxmxd(mgs,li))
        cscnvi(mgs) = qscnvi(mgs)*rho0(mgs)/xmas(mgs,li)
        cscnvis(mgs) = 0.5*cscnvi(mgs)

       ENDIF
      ENDIF
!      ENDIF
      end do



      IF ( ipelec >= 1 .or. idoniconly ) THEN
!
!  initialize rates to zero
!
      scxacy(:,:,:) = 0.0

      scsacw(:) = 0.0
      scsaci(:) = 0.0
      scsacis(:) = 0.0
      scsacr(:) = 0.0

      schacw(:) = 0.0
      schaci(:) = 0.0
      schacis(:) = 0.0
      schacs(:) = 0.0
      schacr(:) = 0.0

!      IF ( lhl .gt. 1 ) THEN
      schlacw(:) = 0.0
      schlaci(:) = 0.0
      schlacis(:) = 0.0
      schlacs(:) = 0.0
      schlacr(:) = 0.0
!      ENDIF

!      IF ( lf .gt. 1 ) THEN
      scfacw(:) = 0.0
      scfaci(:) = 0.0
      scfacs(:) = 0.0
      scfacr(:) = 0.0
!      ENDIF 
      
      ENDIF


      IF ( elec_ramp_time > 0 .and. time_real > elec_on_time .and. time_real <  ( elec_on_time + elec_ramp_time ) ) THEN
        
        elecfac = (time_real - elec_on_time)/elec_ramp_time
        
      ELSEIF (  time_real >= ( elec_on_time + elec_ramp_time ) ) THEN
        elecfac = 1.0
      ELSEIF ( time_real < elec_on_time ) THEN
        elecfac = 0.0
      ENDIF
      
      altelecfac = 1.0
      ! altelecfac = elecfac



      IF ( ipelec .ge. 2 .or. idoniconly ) THEN

!      write(*,*) 'icezvd_gs: time_real, elec_on_time = ',time_real,elec_on_time
      IF ( nonigrd .ne. -10 .and. time_real > elec_on_time .and.  &
     &    ( (ny .eq. 2 ) .or.    &
     &     ( bcy .eq. 2 ) .or.    &
     &     ( ny .gt. 2 .and. ( jybeg-1+jy .ge. jchgs .and. jybeg-1+jy .lt. nyend-1-jchgn ) ))) THEN
      

      
      
      DO mgs = 1,ngscnt


      exy(mgs,ls,lc) = esw(mgs)
      exy(mgs,ls,li) = esi(mgs)
      exy(mgs,lh,lc) = ehw(mgs)
      exy(mgs,lh,li) = ehi(mgs)
      exy(mgs,lh,ls) = ehs(mgs)
!      IF ( lis > 1 ) exy(mgs,ls,lis) = esi(mgs)
      IF ( lis > 1 ) exy(mgs,lh,lis) = ehis(mgs)

      IF ( lhl .gt. 1 ) THEN
        exy(mgs,lhl,lc) = ehlw(mgs)
        exy(mgs,lhl,li) = ehli(mgs)
        exy(mgs,lhl,ls) = ehls(mgs)
        IF ( lis > 1 ) exy(mgs,lhl,lis) = ehlis(mgs)
      ENDIF

      IF ( lf > 1 ) THEN
        exy(mgs,lf,lc) = efw(mgs)
        exy(mgs,lf,li) = efi(mgs)
        exy(mgs,lf,ls) = efs(mgs)
        IF ( lis > 1 ) exy(mgs,lf,lis) = efis(mgs)
      ENDIF

      cxacy(mgs,ls,lc) = csacw(mgs)
      cxacy(mgs,lh,lc) = chacw(mgs)
      
         ! use ni as the fraction of ice crystals greater than some given diameter
         ni = 1.0
         nis = 1.0
         IF ( cidiamin > 0.0 ) THEN
           ni = Exp(- (cidiamin/xdia(mgs,li,1))**3 )
           IF ( lis > 1 ) nis = Exp(- (cidiamin/xdia(mgs,lis,1))**3 )
         ENDIF

      
      IF ( .true. ) THEN
      cxacy(mgs,ls,li) = ni*csaci0(mgs)
      cxacy(mgs,lh,li) = ni*chaci0(mgs)
      cxacy(mgs,lh,ls) = chacs0(mgs)
        IF ( lis > 1 ) THEN
          cxacy(mgs,lh,lis) = nis*chacis0(mgs)
          cxacy(mgs,ls,lis) = nis*csacis0(mgs)
        ENDIF
      ELSE
      cxacy(mgs,ls,li) = csaci(mgs) ! csaci0(mgs)
      cxacy(mgs,lh,li) = chaci(mgs) ! chaci0(mgs)
      cxacy(mgs,lh,ls) = chacs(mgs) ! chacs0(mgs)
      ENDIF

!         IF ( igs(mgs) == 21 .and. temg(mgs) < 273. .and. temg(mgs) > 233. ) THEN
!           write(0,*) 'k,ni,chaci0 = ',kgs(mgs),ni,chaci0(mgs),cxacy(mgs,lh,li)
!         ENDIF

      IF ( lf > 1 ) THEN
        cxacy(mgs,lf,lc) = cfacw(mgs)
        cxacy(mgs,lf,li) = ni*cfaci0(mgs)
        cxacy(mgs,lf,ls) = cfacs0(mgs)
!        IF ( ny <= 2  ) THEN
!         IF ( igs(mgs) == 21 .and. temg(mgs) < 273. .and. temg(mgs) > 233. ) THEN
!           write(0,*) 'k,ni,cfaci0 = ',kgs(mgs),ni,cfaci0(mgs),cxacy(mgs,lf,li)
!         ENDIF
!        ENDIF

      ENDIF
      IF ( lhl .gt. 1 ) THEN
        cxacy(mgs,lhl,lc) = chlacw(mgs)
        cxacy(mgs,lhl,li) = ni*chlaci0(mgs)
        cxacy(mgs,lhl,ls) = chlacs0(mgs)
      ENDIF


        IF ( iremoveqwfrz == 1 ) THEN
          qcwtmp(mgs) = Max( 0.0, qx(mgs,lc) - dtp*qwfrz(mgs) )
        ELSE
          qcwtmp(mgs) = Max( 0.0, qx(mgs,lc) )
        ENDIF

      ENDDO


!
!  First do loop between large ice and small ice
!
      DO il=ls,lhab ! loop through all ice habits
      DO ic=li,ls ! loop through all mixed phase particles

!      IF ( iexy(il,ic) .eq. 1 .and. il .ne. ic ) THEN ! prevent snow-snow
      IF ( iexy(il,ic) .eq. 1  ) THEN

      do mgs = 1,ngscnt

!        IF (qx(mgs,il) .gt. qxmin(il)  .and. qx(mgs,ic) .gt. qxmin(ic)  .and. exy(mgs,il,ic) .gt. 0.0) THEN  ! make sure there are particle and that they do not stick
        IF (qx(mgs,il) .gt. qxmin(il)  .and. qx(mgs,ic) .gt. qxmin(ic) ) THEN  ! make sure there are particle and that they do not stick
        IF ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) THEN  ! 0.0=tfreezing and -38C= assumed T for homogeneous ice nucleation
!        IF ( .true. ) THEN

!  del cglarge
!
!
      IF (nonigrd .eq. -1 ) THEN ! {

       vt = vtxbar(mgs,il,1) - vtxbar(mgs,ic,1)
       CALL takax(isaund,ntem,nlwc,takalu,temcg(mgs),qcwtmp(mgs),vt,   &
     &     xdia(mgs,ic,1),rho0(mgs),ftelwc,exy(mgs,il,lc),qx(mgs,lr),rarfac)
       scxacy(mgs,il,ic) = ftelwc

        if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .gt. delqxxa(ic) .and. il .ge. lg ) THEN
              ngscmaxy(il,ic) = ngscmaxy(il,ic) + 1
          ENDIF
          scxacy(mgs,il,ic) = min(delqxxa(ic),scxacy(mgs,il,ic))
        end if
        if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .lt. delqnxa(ic) .and. il .ge. lg ) THEN
              ngscminy(il,ic) = ngscminy(il,ic) + 1
          ENDIF
          scxacy(mgs,il,ic) = max(delqnxa(ic),scxacy(mgs,il,ic))
        end if
      ! }
      ELSEIF ( nonigrd .eq. 1 .or. nonigrd .eq. 0 ) THEN ! {

      IF (nonigrd .eq. 0) THEN

      qconkq = 0.0
      qconm = 1.0
      qconn = 1.0
      ftelwc = 0.0
      qsign = 0.0 ! ,ftelwc,qconkq,qconm,qconn

       IF ( qcwtmp(mgs) > qxmin(lc) ) THEN

      IF ( iraropt == 1 .or. il < lh ) THEN
        vt = vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
      ELSEIF ( iraropt == 2 ) THEN ! num-wgt
        vt = vtxbar(mgs,il,2) - vtxbar(mgs,lc,1)
      ELSEIF ( iraropt == 3 ) THEN ! area-wgt
           ! use Milbrandt-MC (2010) eq. 9 for weighted fall speed ratio
           k = 2
           j = 3
           tmp = alpha(mgs,il) + 1.0 + bxx(mgs,il) + k ! 2 for area wgt (k)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + k ! 2 for area wgt (k)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + bxx(mgs,il) + j ! 3 for mass wgt (j)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + j ! 3 for mass wgt (j)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
           fac =  x1*y2/(y1*x2) ! Vk/Vj
!            write(0,*) 'fac,alp,x1,x2,y1,y2 = ',fac,alpha(mgs,il),x1,x2,y1,y2
!            write(0,*) 'v3,v0,v2 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),fac*vtxbar(mgs,il,1)
         vt = fac*vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
      ELSEIF ( iraropt == 4 ) THEN ! diameter-wgt
           ! use Milbrandt-MC (2010) eq. 9 for weighted fall speed ratio
           k = 1
           j = 3
           tmp = alpha(mgs,il) + 1.0 + bxx(mgs,il) + k ! 2 for area wgt (k)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + k ! 2 for area wgt (k)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + bxx(mgs,il) + j ! 3 for mass wgt (j)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,il) + 1.0 + j ! 3 for mass wgt (j)
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
           fac =  x1*y2/(y1*x2) ! Vk/Vj
!            write(0,*) 'fac,alp,x1,x2,y1,y2 = ',fac,alpha(mgs,il),x1,x2,y1,y2
!            write(0,*) 'v3,v0,v1 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),fac*vtxbar(mgs,il,1)
         vt = fac*vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
      ELSE ! IF ( iraropt == 5 ) THEN
        vt = 0.5*(vtxbar(mgs,il,1)+vtxbar(mgs,il,2)) - vtxbar(mgs,lc,1)
      ENDIF
      call saundx(isaund,temcg(mgs),qcwtmp(mgs),exy(mgs,il,lc),vt,xdia(mgs,ic,1),   &
     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac)
     
       ENDIF ! qc > qxmin
!      IF ( ftelwc /= 0.0 ) THEN
!      IF ( igs(mgs) == 42 .and. kgs(mgs) == 43 ) THEN
!      write(0,*) 'called saundx: ',igs(mgs),kgs(mgs),il,lc,isaund,temcg(mgs),qcwtmp(mgs),exy(mgs,il,lc),vt,xdia(mgs,ic,1),   &
!     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac
!      ENDIF
      
!      IF ( lis > 1 .and. ic == lis .and. il == lh ) THEN
!       IF ( cxacy(mgs,il,ic) > 1. .and. qcwtmp(mgs) > 0.05e-3 ) THEN
!       write(0,*) 'called saundx: ',igs(mgs),kgs(mgs),il,lc,isaund,temcg(mgs),qcwtmp(mgs),exy(mgs,il,lc),vt,xdia(mgs,ic,1),   &
!     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac
!       ENDIF
!      ENDIF

      ELSE
        vt = vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
        
      call saundy(isaund,temcg(mgs),qcwtmp(mgs),rarx(mgs,il),vt,xdia(mgs,ic,1),   &
     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq)

!      call saund(temcg(mgs),qcwtmp(mgs),eglw(mgs),xdia(mgs,ic,1),
!     >     rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq)
      ENDIF

      cfce15 = 1.0e-15
      scxacy(mgs,il,ic) = qconkq*(xdia(mgs,ic,1)**qconm)   &
     &  *((abs(vtxbar(mgs,il,1)-vtxbar(mgs,ic,1)))**qconn)   &
     &  *ftelwc*cfce15
      
      tmp = scxacy(mgs,il,ic)
      if ( scxacy(mgs,il,ic) .gt. 0.0 ) then ! {
          IF ( scxacy(mgs,il,ic) .gt. delqxxa(ic) .and. il .ge. lg ) THEN
              ngscmaxy(il,ic) = ngscmaxy(il,ic) + 1
          ENDIF
      scxacy(mgs,il,ic) = min(delqxxa(ic),scxacy(mgs,il,ic))
      end if
      if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .lt. delqnxa(ic) .and. il .ge. lg ) THEN
              ngscminy(il,ic) = ngscminy(il,ic) + 1
          ENDIF
      scxacy(mgs,il,ic) = max(delqnxa(ic),scxacy(mgs,il,ic))
      
!      IF ( ny <= 2 .and. ( il == lh .or. il == lf ) .and. ic == li ) THEN
!       IF ( igs(mgs) == 21 ) THEN
!         IF ( il == lh ) THEN
!         write(0,*) 'lh with li: k = ',kgs(mgs)
!         ELSE
!         write(0,*) 'lf with li: k = ',kgs(mgs)
!         ENDIF
!         write(0,*) 'scxacy,vt,xdia = ',scxacy(mgs,il,ic),vt,xdia(mgs,ic,1)
!       ENDIF
!      ENDIF
      
      end if ! }
      end if ! }
!
      if ( nonigrd .ge. 2 ) then
      IF ( dellwc(mgs) /= 0.0 ) THEN
      scxacy(mgs,il,ic) =    &
     & ( (7.3)*rgard1   &
     & * dellwc(mgs)*sctem(mgs)*(xdia(mgs,ic,1)**4)   &
     & * (abs(vtxbar(mgs,il,1)-vtxbar(mgs,ic,1))**3) )
      ELSE
       scxacy(mgs,il,ic) = 0.0
      ENDIF
      
      if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .gt. delqxxa(ic) .and. il .ge. lg ) THEN
              ngscmaxy(il,ic) = ngscmaxy(il,ic) + 1
          ENDIF
      scxacy(mgs,il,ic) = min(delqxxa(ic),scxacy(mgs,il,ic))
      end if
      if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .lt. delqnxa(ic) .and. il .ge. lg ) THEN
              ngscminy(il,ic) = ngscminy(il,ic) + 1
          ENDIF
      scxacy(mgs,il,ic) = max(delqnxa(ic),scxacy(mgs,il,ic))
      end if
      end if
!
      if ( nonigrd .eq. 3 ) then
      if ( scxacy(mgs,il,ic) .gt. 0 )  scxacy(mgs,il,ic) = delqxxa(ic)
      if ( scxacy(mgs,il,ic) .lt. 0 )  scxacy(mgs,il,ic) = delqnxa(ic)
      end if
!
!      sxxacy(mgs,il,ic) = scxacy(mgs,il,ic)
!
! Formula 7 in M05 to get total space charge per collisions within a grid box
      IF ( .true. ) THEN
      tmp = scxacy(mgs,il,ic)
      scxacy(mgs,il,ic) =   &
     &  cxacy(mgs,il,ic)*ecollmx*(1.0-exy(mgs,il,ic))*scxacy(mgs,il,ic) !   &
!     &  *exy(mgs,il,ic)*(1./(exy(mgs,il,ic)+eps))
!      IF ( igs(mgs) == 43 .and. kgs(mgs) == 47 ) THEN
!      write(0,*) 'scxacy: ',igs(mgs),kgs(mgs),il,ic,isaund,tmp,scxacy(mgs,il,ic),temcg(mgs), &
!     & qcwtmp(mgs),exy(mgs,il,ic),vt,xdia(mgs,ic,1),   &
!     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac
!      ENDIF
      ELSE
      scxacy(mgs,il,ic) =   &
     &  cxacy(mgs,il,ic)*Min(ecollmx,1.0-exy(mgs,il,ic))*scxacy(mgs,il,ic)   &
     &  *(1./(exy(mgs,il,ic)+eps))
      ENDIF

!      IF ( ny <= 2 .and. ( il == lh .or. il == lf ) .and. ic == li ) THEN
!       IF ( igs(mgs) == 21 ) THEN
!         IF ( il == lh ) THEN
!         write(0,*) '2 lh with li: k = ',kgs(mgs)
!         ELSE
!         write(0,*) '2 lf with li: k = ',kgs(mgs)
!         ENDIF
!         write(0,*) 'scxacy,cxacy,1-exy = ',scxacy(mgs,il,ic),cxacy(mgs,il,ic),(1.0-exy(mgs,il,ic))
!       ENDIF
!      ENDIF

     
      IF ( Abs( scxacy(mgs,il,ic) ) .gt. scxacymax ) THEN
           ngscxymax(il,ic) = ngscxymax(il,ic) + 1
           scxacy(mgs,il,ic) = Min(scxacymax, Max( -scxacymax, scxacy(mgs,il,ic) ) )
      ENDIF
!
      end if ! } temperture allows mixed phase


!
! charging with no liquid water (Mitzeva).  Only allow when ssi >= 1, i.e., for the case of deposition growth
!
      IF ( nic_noliq .ne. 0 .and. scxacy(mgs,il,ic) .eq. 0.0 .and. &
     &         temg(mgs) .lt. tfr .and. temg(mgs) > tfrh ) THEN  ! will overwrite pre-existing charge sep rate if no scxacy test

        IF ( ssi(mgs) .ge. 1.0 ) THEN

        IF ( nic_noliq == 1 .and. ssi(mgs) .ge. 1.0 ) THEN

         vt = vtxbar(mgs,il,1) - vtxbar(mgs,ic,1)
        
         call saundmst(temcg(mgs),ssi(mgs),vt,xdia(mgs,ic,1),rho0(mgs),   &
     &     qsign,ftelwc,qcwtmp(mgs),exy(mgs,il,lc),   &
     &     qconkq,qconm,qconn,idelq,rarfac)

!      subroutine saundmst(temcg,ssi(mgs),vt,awdia,rho0,
!                          qsign,ftrar,qcw,exw,
!                          qconkq,qconm,qconn,idelq,rarfac)
        cfce15 = 1.0e-15
        ftelwc = qsign*qc_noliq ! qmst
      
        scxacy(mgs,il,ic) = qconkq*(xdia(mgs,ic,1)**qconm)   &
     &                            *((abs(vt))**qconn)*ftelwc*cfce15     ! MST, eqn 1

        ELSEIF ( nic_noliq >= 2 .and. ssi(mgs) > 1.01 ) THEN
        
         vt = vtxbar(mgs,il,1) - vtxbar(mgs,ic,1)
         IF ( il == ls .and. ic == li ) THEN
           vt = Max( vt, 2.0 )
         ENDIF

        call saundx(14,temcg(mgs),qx(mgs,lc),exy(mgs,il,lc),vt,xdia(mgs,ic,1),   &
     &   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac)

      cfce15 = 1.0e-15
      scxacy(mgs,il,ic) = qconkq*(xdia(mgs,ic,1)**qconm)   &
     &  *abs(vt)**qconn*ftelwc*cfce15
      
        
        ENDIF
        
        if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .gt. delqxxa(ic) .and. il .ge. lg ) THEN
              ngscmaxy(il,ic) = ngscmaxy(il,ic) + 1
          ENDIF
         scxacy(mgs,il,ic) = min(delqxxa(ic),scxacy(mgs,il,ic))
        end if
        if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
          IF ( scxacy(mgs,il,ic) .lt. delqnxa(ic) .and. il .ge. lg ) THEN
              ngscminy(il,ic) = ngscminy(il,ic) + 1
          ENDIF
         scxacy(mgs,il,ic) = max(delqnxa(ic),scxacy(mgs,il,ic))
        end if

!        scxacy(mgs,il,ic) =   &
!     &    cxacy(mgs,il,ic)*Min(ecollmx,1.0-exy(mgs,il,ic))*scxacy(mgs,il,ic)   &
!     &    *(1./(exy(mgs,il,ic)+eps))
        scxacy(mgs,il,ic) =   &
     &    cxacy(mgs,il,ic)*ecollmx*(1.0-exy(mgs,il,ic))*scxacy(mgs,il,ic) 

          
        IF ( ny .eq. 2 .and. scxacy(mgs,il,ic) .ne. 0.0 ) THEN
           cntnic_noliq = cntnic_noliq + 1
           q_noliqmn = Min( q_noliqmn, scxacy(mgs,il,ic) )
           q_noliqmx = Max( q_noliqmx, scxacy(mgs,il,ic) )
           
           IF ( il == ls .and. ic == li ) THEN
             scsacimn = Min( scsacimn, scxacy(mgs,il,ic) )
             scsacimx = Max( scsacimx, scxacy(mgs,il,ic) )
           ENDIF
           
        ENDIF

        tq1(igs(mgs),jy,kgs(mgs),4) =  tq1(igs(mgs),jy,kgs(mgs),4) + scxacy(mgs,il,ic) + scxacy(mgs,il,ic) !temporary to output liquid-free charging rate
        ENDIF !ssi

      ENDIF !nic_noliq
      end if ! } qx(mgs,il) .gt. qxmin(il) .and. qx(mgs,ic) .gt. qxmin(ic)

!
      end do !mgs
      
      ENDIF ! iexy(il,ic) .eq. 1 .and. il .ne. ic 
      
      ENDDO !ic
      ENDDO !il

!
!  End of ice-ice loop
!

!
! Loop for graupel-droplet inductive charging
!  Mansell et al 2005, adapted from Ziegler et al 1991

      IF ( erbnd .gt. 0.0 .and. fdgt .gt. 0.0 .and. ipelec .ge. 3 ) THEN
      
      DO il=ls,lhab
!
      IF ( iexy(il,lc) .eq. 1 ) THEN

      if ( ndebug .ge. 1 ) write(iunit,*) 'elec 19'

      do mgs = 1,ngscnt
!
      if ( qx(mgs,il) .gt. qxmin(il) .and. qcwtmp(mgs) .gt. qxmin(lc)    &
     &     .and. exy(mgs,il,lc) .gt. 0.0 .and. cx(mgs,il) > 1.e-6 ) THEN
      IF ( temg(mgs) .gt. tindmn .and. temg(mgs) .lt. tindmx ) then
!
      scxacy(mgs,il,lc) =  &
     &  0.125*pi**3*exy(mgs,il,lc)*erbnd*fdgt*(cx(mgs,il))   &
     &  *cx(mgs,lc)*vtxbar(mgs,il,1)*6.0/gf4p5*(xdia(mgs,lc,1)**2)   &
     &  *( pi*eperao*costhe*dezcomp(mgs)   &
     &      *gf3p5*(xdia(mgs,il,1)**2)   &
     &     - (1.0/3.0)*gf1p5*scx(mgs,il)/cx(mgs,il) ) ! must use from previous time step
        
        
!        scxacy(mgs,il,lc) = 0.0
        
!    !     makes sure scxacy does nto go over |delqxw|.

      if ( scxacy(mgs,il,lc) .gt. 0.0 ) then
        scxacy(mgs,il,lc) = min(delqxw,scxacy(mgs,il,lc))
      elseif ( scxacy(mgs,il,lc) .lt. 0.0 ) then
!        IF ( scxacy(mgs,il,lc) < delqnw ) THEN
!          write(0,*) 'big inductive: il = ',il,scxacy(mgs,il,lc),temcg(mgs),kgs(mgs)
!          write(0,*) 'cx,vtxbar,scx = ',cx(mgs,il),vtxbar(mgs,il,1),scx(mgs,il),xdia(mgs,lc,1)
!        ENDIF
        scxacy(mgs,il,lc) = max(delqnw,scxacy(mgs,il,lc))
      end if

!       if (abs(scxacy(mgs,il,lc)).gt.0.) write(0,*) 'GT0',cx(mgs,il),cx(mgs,lc),vtxbar(mgs,il,1),xdia(mgs,lc,1),dezcomp(mgs),xdia(mgs,il,1),scx(mgs,il) 

      ELSE

       scxacy(mgs,il,lc) = 0.0

      end if
      ENDIF
!
!      end if
!      end if

      end do  ! mgs

      ENDIF

      ENDDO  ! il

      ENDIF ! erbnd

     ENDIF ! boundary check on y

!      if (MAXVAL(scxacy(:,ls,lc)).gt.0.) write(0,*)'dezcomp',MAXVAL(dezcomp),MAXLOC(dezcomp),MINVAL(dezcomp),MINLOC(dezcomp)
!      if (MAXVAL(scxacy(:,ls,lc)).gt.0.) write(0,*) 'scxacy(mgs,ls,lc)',MAXVAL(scxacy(:,ls,lc)),MAXLOC(scxacy(:,ls,lc)),MINVAL(scxacy(:,ls,lc)),MINLOC(scxacy(:,ls,lc))


!  SET INDUC to ZERO:

!      DO il=ls,lhab
!      DO mgs = 1,ngscnt
!      if (abs(scxacy(mgs,il,lc)).le.delqxw) scxacy(mgs,il,lc)=0.
!      ENDDO
!      ENDDO


      DO mgs = 1,ngscnt

      IF ( .not. ( bcx .ne. 2 .and. ( ixbeg-1+igs(mgs) .le. ichge .or. ixbeg-1+igs(mgs) .ge. nxend-1-ichgw )) ) THEN
      
      scsacw(mgs) = elecfac*scxacy(mgs,ls,lc)
      scsaci(mgs) = elecfac*scxacy(mgs,ls,li)
      scsacr(mgs) = elecfac*scxacy(mgs,ls,lr)

      schacw(mgs) = elecfac*scxacy(mgs,lh,lc)
      schaci(mgs) = elecfac*scxacy(mgs,lh,li)
      schacs(mgs) = elecfac*scxacy(mgs,lh,ls)
      schacr(mgs) = elecfac*scxacy(mgs,lh,lr)
      
      
      IF ( lhl .gt. 1 ) THEN
      schlacw(mgs) = elecfac*scxacy(mgs,lhl,lc)
      schlaci(mgs) = elecfac*scxacy(mgs,lhl,li)
      schlacs(mgs) = elecfac*scxacy(mgs,lhl,ls)
      schlacr(mgs) = elecfac*scxacy(mgs,lhl,lr)
      ENDIF

      IF ( lf .gt. 1 ) THEN
        scfacw(mgs) = elecfac*scxacy(mgs,lf,lc)
        scfaci(mgs) = elecfac*scxacy(mgs,lf,li)
        scfacs(mgs) = elecfac*scxacy(mgs,lf,ls)
        scfacr(mgs) = elecfac*scxacy(mgs,lf,lr)
      ENDIF
      
      ELSE
      
      scsacw(mgs) = 0.0
      scsaci(mgs) = 0.0
      scsacr(mgs) = 0.0

      schacw(mgs) = 0.0
      schaci(mgs) = 0.0
      schacs(mgs) = 0.0
      schacr(mgs) = 0.0
      
      IF ( lhl .gt. 1 ) THEN
      schlacw(mgs) = 0.0
      schlaci(mgs) = 0.0
      schlacs(mgs) = 0.0
      schlacr(mgs) = 0.0
      ENDIF

      IF ( lf .gt. 1 ) THEN
        scfacw(mgs) = 0.0
        scfaci(mgs) = 0.0
        scfacs(mgs) = 0.0
        scfacr(mgs) = 0.0
      ENDIF

      
      ENDIF ! x boundary check

      ENDDO  ! ngscnt


      ENDIF  ! ipelec >= 1


!
!  Ventilation coeficients
!
      do mgs = 1,ngscnt
      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))
      end do
!
!
      if ( ndebug .gt. 0 ) write(0,*) 'civent'
!
      civenta = 1.258e4
      civentb = 2.331
      civentc = 5.662e4
      civentd = 2.373
      civente = 0.8241
      civentf = -0.042
      civentg = 1.70

      do mgs = 1,ngscnt
      IF ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) THEN
      IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
      cireyn =   &
     &  (civenta*xdia(mgs,li,1)**civentb   &
     &  +civentc*xdia(mgs,li,1)**civentd)   &
     &  /   &
     &  (civente*xdia(mgs,li,1)**civentf+civentg)
      xcivent = (fschm(mgs)**(1./3.))*((cireyn/fakvisc(mgs))**0.5)
      if ( xcivent .lt. 1.0 ) then
      civent(mgs) = 1.0 + 0.14*xcivent**2
      end if
      if ( xcivent .ge. 1.0 ) then
      civent(mgs) = 0.86 + 0.28*xcivent
      end if
      ELSE
       civent(mgs) = 0.0
      ENDIF


      ENDIF ! icond .eq. 1
      end do

!
!
      igmrwa = 100.0*2.0
      igmrwb = 100.*((5.0+br)/2.0)
      rwventa = (0.78)*gmoi(igmrwa)  ! 0.78
      rwventb = (0.308)*gmoi(igmrwb) ! 0.562825
      do mgs = 1,ngscnt
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( ipconc .ge. 3 ) THEN
          IF ( imurain == 3 ) THEN
           IF ( izwisventr == 1 ) THEN
            rwvent(mgs) = ventrx(mgs)*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)
           ELSE ! izwisventr = 2
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
          rwvent(mgs) =   &
     &  (0.78*ventrx(mgs) + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
           ENDIF

          ELSE ! imurain == 1
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lr)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        IF ( iferwisventr == 1 ) THEN

  ! Ferrier fall speed in the ventillation term [uses fx(lr) ]
  
        alpr = Min(alpharmax,alpha(mgs,lr) )

        x =  1. + alpha(mgs,lr)

        IF ( ipconc >= 6 .and. lzr > 1 ) THEN ! 3 moment
        tmp = 1. + alpr ! alpha(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lr) + 0.5*bx(lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions
        ELSE
         y = ventrxn(mgs)
        ENDIF

!         vent1 = dble(xdia(mgs,lr,1))**(-2. - alpr) ! Actually OK
!         vent2 = dble(1./xdia(mgs,lr,1) + 0.5*fx(lr))**dble(2.5+alpr+0.5*bx(lr))  ! Actually OK
         vent1 = dble(xdia(mgs,lr,1))**(0.5 + 0.5*bx(lr)) ! 2016.2.26 Changed for consistency with derivation (recast formula -- should be equivalent)
         vent2 = dble(1. + 0.5*fx(lr)*xdia(mgs,lr,1))**dble(2.5+alpr+0.5*bx(lr))
        
        
        rwvent(mgs) =    &
     &    0.78*x +    &
     &    0.308*fvent(mgs)*y*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)
       
        rwventz(mgs) = 0.0

!        rwventz(mgs) =    &
!     &    0.78*x +    &
!     &    0.308*fvent(mgs)*y*   &
!     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)


        ELSEIF ( iferwisventr == 2 ) THEN
          
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
         x =  1. + alpha(mgs,lr)

           rwvent(mgs) =   &
     &  (0.78*x + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )


        IF ( ipconc >= 7 ) THEN
        ! vent coeff. for reflectivity rate from evaporation
        alpr = Min(alpharmax,alpha(mgs,lr) )

           tmp = alpr + 5.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!        rwventz(mgs) =    &
!     &    0.78*(4. + alpha(mgs,lr))*(3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr)) +    &
        rwventz(mgs) =    &
     &    0.78*(4. + alpr)*(3. + alpr)*(2. + alpr)*(1. + alpr) +    &
     &    0.308*fvent(mgs)*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(y/gf1palp(mgs))*(xdia(mgs,lr,1)**((1.0+br)/2.0))

        ENDIF

          
          ENDIF ! iferwisventr
          
          ENDIF ! imurain
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
      igmswa = 100.0*2.0
      igmswb = 100.*((5.0+ds)/2.0)
      swventa = (0.78)*gmoi(igmswa)
      swventb = (0.308)*gmoi(igmswb)
      do mgs = 1,ngscnt
      IF ( qx(mgs,ls) .gt. qxmin(ls) ) THEN
      IF ( ipconc .ge. 4 ) THEN
      swvent(mgs) = 0.65 + 0.44*fvent(mgs)*Sqrt(vtxbar(mgs,ls,1)*xdia(mgs,ls,1))
      ELSE
! 10-ice version:
       swvent(mgs) =   &
     &  (swventa + swventb*fvent(mgs)   &
     &   *Sqrt((cs*rhovt(mgs)))   &
     &   *(xdia(mgs,ls,1)**((1.0+ds)/2.0)) )
      ENDIF
      ELSE
      swvent(mgs) = 0.0
      ENDIF
      end do
!
!

      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
!      hwventc = (4.0*gr/(3.0*cdx(lh)))**(0.25)
      hwvent(:) = 0.0
      hwventy(:) = 0.0

      do mgs = 1,ngscnt
      IF ( qx(mgs,lh) .gt. qxmin(lh) ) THEN
       IF ( icdx /= 6 .and. alpha(mgs,lh) .eq. 0.0 ) THEN
       hwventc = (4.0*gr/(3.0*cdxgs(mgs,lh)))**(0.25)
        hwvent(mgs) =   &
     &  ( hwventa + hwventb*hwventc*fvent(mgs)   &
     &    *((xdn(mgs,lh)/rho0(mgs))**(0.25))   &
     &    *(xdia(mgs,lh,1)**(0.75)))
       ELSE ! Ferrier 1994, eq. B.36
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lh)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
! note that hwvent includes a division by Gamma(1+alpha), so Gamma(2+alpha)/Gamma(1+alpha) = 1 + alpha
! and g1palp = Gamma(1+alpha) divides into y
        x =  1. + alpha(mgs,lh)

        tmp = 1 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        
        hwventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*Sqrt(axx(mgs,lh)*rhovt(mgs)) 
        hwvent(mgs) =    &
     &  ( 0.78*x +  y*hwventy(mgs) ) !   &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       ENDIF
      ELSE
      hwvent(mgs) = 0.0
      hwventy(mgs) = 0.0
      ENDIF
      end do
      
      fwvent(:) = 0.0
      fwventy(:) = 0.0

      IF ( lf .gt. 1 ) THEN
      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
!      hwventc = (4.0*gr/(3.0*cdx(lf)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lf) .gt. qxmin(lf) ) THEN

       IF ( icdx /= 6 .and. alpha(mgs,lf) .eq. 0.0 ) THEN
        hwventc = (4.0*gr/(3.0*cdxgs(mgs,lf)))**(0.25)
        fwvent(mgs) =   &
     &  ( hwventa + hwventb*hwventc*fvent(mgs)   &
     &    *((xdn(mgs,lf)/rho0(mgs))**(0.25))   &
     &    *(xdia(mgs,lf,1)**(0.75)))
       ELSE ! Ferrier 1994, eq. B.36
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lf)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

! note that fwvent includes a division by Gamma(1+alpha), so x = Gamma(2+alpha)/Gamma(1+alpha) = 1 + alpha
! and g1palp = Gamma(1+alpha) divides into y

        x =  1. + alpha(mgs,lf)

        tmp = 1 + alpha(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lf) + 0.5*bxx(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions

        fwventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lf,1)**(0.5 + 0.5*bxx(mgs,lf)))*Sqrt(axx(mgs,lf)*rhovt(mgs)) 
        
        fwvent(mgs) =  0.78*x + y*fwventy(mgs) 

        ENDIF
       ENDIF
      end do
      ENDIF


      hlvent(:) = 0.0
      hlventy(:) = 0.0

      IF ( lhl .gt. 1 ) THEN
      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
!      hwventc = (4.0*gr/(3.0*cdx(lhl)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lhl) .gt. qxmin(lhl) ) THEN

       IF ( icdxhl /= 6 .and. alpha(mgs,lhl) .eq. 0.0 ) THEN
        hwventc = (4.0*gr/(3.0*cdxgs(mgs,lhl)))**(0.25)
        hlvent(mgs) =   &
     &  ( hwventa + hwventb*hwventc*fvent(mgs)   &
     &    *((xdn(mgs,lhl)/rho0(mgs))**(0.25))   &
     &    *(xdia(mgs,lhl,1)**(0.75)))
       ELSE ! Ferrier 1994, eq. B.36
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lhl)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

! note that hlvent includes a division by Gamma(1+alpha), so x = Gamma(2+alpha)/Gamma(1+alpha) = 1 + alpha
! and g1palp = Gamma(1+alpha) divides into y

        x =  1. + alpha(mgs,lhl)

        tmp = 1 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions

        hlventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*Sqrt(axx(mgs,lhl)*rhovt(mgs)) 
        
        hlvent(mgs) =  0.78*x + y*hlventy(mgs)  !   &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*   &
!     &            Sqrt(axx(mgs,lhl)*rhovt(mgs)))
!     :            Sqrt(xdn(mgs,lhl)*ax(lhl)*rhovt(mgs)/rg0))/tmp

        ENDIF
       ENDIF
      end do
      ENDIF

!
!
!
!  Wet growth constants
!
      do mgs = 1,ngscnt
      fwet1(mgs) =   &
     & (2.0*pi)*   &
     & ( felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qss0(mgs)-qx(mgs,lv))   &
     &  -ftka(mgs)*temcg(mgs) )   &
     & / ( rho0(mgs)*(felf(mgs)+fcw(mgs)*temcg(mgs)) )
      fwet2(mgs) =   &
     &  (1.0)-fci(mgs)*temcg(mgs)   &
     & / ( felf(mgs)+fcw(mgs)*temcg(mgs) )
      end do
!
!  Melting constants
!
      do mgs = 1,ngscnt
      fmlt1(mgs) = (2.0*pi)*   &
     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))   &
     &   -ftka(mgs)*temcg(mgs)/rho0(mgs) )    &
     &  / (felf(mgs))
      fmlt2(mgs) = -fcw(mgs)*temcg(mgs)/felf(mgs)
      fmlt1e(mgs) = (2.0*pi)*   &
     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))  ) / (felf(mgs))
      end do
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
!  deposition, sublimation, and melting of snow, graupel and hail
!
      qsmlr(:) = 0.0
      qimlr(:) = 0.0 ! this is not used. qi melts to qc way down in the code.
      qhmlr(:) = 0.0
      qhlmlr(:) = 0.0
      IF ( lhwlg > 1 ) THEN
        qhmlrlg(:) = 0.0
        qhlmlrlg(:) = 0.0
      ENDIF
      qhfzh(:) = 0.0
      qffzf(:) = 0.0
      qhlfzhl(:) = 0.0
      qhfzhlg(:) = 0.0
      qhlfzhllg(:) = 0.0
      vhfzh(:) = 0.0
      vffzf(:) = 0.0
      vhlfzhl(:) = 0.0
      qsfzs(:) = 0.0
!      zsmlr(:) = 0.0
      zhmlr(:) = 0.0
      zhmlrr(:) = 0.0
      zsmlrr(:) = 0.0
      zhshr(:) = 0.0
      zhlmlr(:) = 0.0
      zhlshr(:) = 0.0

      zhshrr(:) = 0.0
      zhlmlrr(:) = 0.0
      zhlshrr(:) = 0.0

      csmlr(:) = 0.0
      csmlrr(:) = 0.0
      chmlr(:) = 0.0
      chmlrr(:) = 0.0
      chlmlr(:) = 0.0
      chlfmlr(:) = 0.0
!      chlmlrsave(:) = 0.0
!      qhlmlrsave(:) = 0.0
!      chlsave(:) = 0.0
!      qhlsave(:) = 0.0
      chlmlrr(:) = 0.0
       qfmlr(:) = 0.0
       cfmlr(:) = 0.0
       cfmul1(:) = 0.0
       cfmlrr(:) = 0.0
       zfmlr(:) = 0.0
       zfmlrr(:) = 0.0
       zfshr(:) = 0.0
       zfshrr(:) = 0.0

      dxshdrate(:,:) = 0.0

      if ( .not. mixedphase ) then !{
      do mgs = 1,ngscnt
!
      IF ( temg(mgs) .gt. tfr ) THEN
      
      IF (  qx(mgs,ls) .gt. qxmin(ls) ) THEN
      qsmlr(mgs) =   &
     &   min(   &
     &  (c1sw*fmlt1(mgs)*cx(mgs,ls)*swvent(mgs)*xdia(mgs,ls,1) ) & ! /rhosm    &
     &   , 0.0 )
      ENDIF

      IF ( lf > 1 ) THEN !{
        IF ( lfw > 1 ) THEN
          tmp = qxw(mgs,lf)
        ELSE
          tmp = 0.0
        ENDIF
      IF (  qx(mgs,lf) - tmp .gt. qxmin(lf) ) THEN !{

!      IF ( ibinhmlr == 0 .or. lzh < 1 ) THEN
       qfmlr(mgs) =   &
     &   meltfac*min(   &
     &  fmlt1(mgs)*cx(mgs,lf)*fwvent(mgs)*xdia(mgs,lf,1)   &
     &  + fmlt2(mgs)*(qfacrmlr(mgs)+qfacwmlr(mgs))    &
     &   , 0.0 )

! Do this later because qfmlr may be adjusted to a smaller value
        cfmlr(mgs)  = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfmlr(mgs)

!       ENDIF

       IF ( ivhmltsoak > 0 .and. qfmlr(mgs) < 0.0 .and. lvol(lf) > 1 .and. xdn(mgs,lf) .lt. xdnmx(lf) ) THEN
         ! act as if 100% of the meltwater were soaked into the graupel
           v1 = (1. - xdn(mgs,lf)/xdnmx(lf))*(vx(mgs,lf) + rho0(mgs)*qfmlr(mgs)/xdn(mgs,lf) )/dtp ! volume available for filling
           v2 = -1.0*rho0(mgs)*qfmlr(mgs)/xdnmx(lf)  ! volume of melted ice if it were refrozen in the matrix
           
           vfsoak(mgs) = Min(v1,v2)
           
       ENDIF
       
       ENDIF !} qf > qxmin
       
       ENDIF !} lf > 1
      
!       IF ( qx(mgs,ls) .gt. 0.1e-4 ) write(0,*) 'qsmlr: ',qsmlr(mgs),qx(mgs,ls),cx(mgs,ls),fmlt1(mgs),
!     :        temcg(mgs),swvent(mgs),xdia(mgs,ls,1),qss0(mgs)-qx(mgs,lv)
!      ELSE
!       qsmlr(mgs) = 0.0
!      ENDIF
! 10ice version:
!     >   min(
!     >  (fmlt1(mgs)*cx(mgs,ls)*swvent(mgs)*xdia(mgs,ls,1) +
!     >   fmlt2(mgs)*(qsacr(mgs)+qsacw(mgs)) )
!     <   , 0.0 )

      IF (  qx(mgs,lh) .gt. qxmin(lh) ) THEN

      IF ( ibinhmlr == 0 .or. lzh < 1 ) THEN
       qhmlr(mgs) =   &
     &   meltfac*min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1)   &
     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacwmlr(mgs))    &
     &   , 0.0 )
       ELSEIF ( ibinhmlr == 1 ) THEN ! use incomplete gamma functions to approximate the bin results


       qhmlr(mgs) =   &
     &   min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1)   &
!     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacw(mgs))    & ! turn off the collection part for now
     &   , 0.0 )

! first part to integrate from mltdiam1 and mltdiam2 to infinity (k loop)

        tmp = 1. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g4palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
 
       DO k = 1,numshedregimes ! do not need mltdiam4 -- treat as infinity; goes to 3 to cover the top 3 diameters (independent of ndiam, the number of smaller diameters)
        
!        ratio = mltdiam1/xdia(mgs,lh,1)
        ratio = Min( maxratiolu, mltdiam(ndiam+k)/xdia(mgs,lh,1) )

        
        IF ( usegamxinf2 ) THEN
           x = gamxinfdp(2. + alpha(mgs,lh), ratio)/g1palp
        ELSE
           x = gaminterp(ratio,alpha(mgs,lh),2,1)
        ENDIF ! usegamxinf

! qxd(ndiam+2), cxd(ndiam+2), qhml(ndiam+2)        
        IF ( usegamxinf2 ) THEN
           cxd1 = cx(mgs,lh)*(1. - gamxinfdp(1. + alpha(mgs,lh), ratio)/g1palp)
        ELSE

          tmp = gaminterp(ratio,alpha(mgs,lh),1,1)

!          IF ( tmp > 0.95 ) THEN
!           tmp = gamxinfdp(1. + alpha(mgs,lh), ratio)/g1palp
!          ENDIF
          cxd1 = cx(mgs,lh)*(1. - tmp)
        ENDIF

         cxd(ndiam+k) = cxd1

        IF ( usegamxinf2 ) THEN
         qxd1 = qx(mgs,lh)*(1. - gamxinfdp(4. + alpha(mgs,lh), ratio)/g4palp)
        ELSE

          tmp = gaminterp(ratio,alpha(mgs,lh),4,1)

!           IF ( tmp > 0.95 ) THEN
!            tmp = gamxinfdp(4. + alpha(mgs,lh), ratio)/g4palp
!           ENDIF
           
           qxd1 = qx(mgs,lh)*(1. - tmp)
        ENDIF

         qxd(ndiam+k) = qxd1

        IF ( usegamxinf2 ) THEN
          y =  gamxinfdp(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp
        ELSE

          y = gaminterp(ratio,alpha(mgs,lh),3,1)

        ENDIF
        
          hwvent1 =  0.78*x + y*hwventy(mgs) 

          qhmlr1 = min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 )

! zero to mltdiam1 is then the difference between qhmlr and qhmlr1
       qhmlr0 = qhmlr(mgs) - qhmlr1
       qhml(ndiam+k) = qhmlr1 ! from diam to infinity
       qhml0(ndiam+k) = qhmlr0 ! from zero to diam

!        IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!          write(iunit,*) 'diam loop k,qhmlr1 = ',k,qhmlr1,qxd1,cxd1
!          write(iunit,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lh),hwvent1,x,y,hwventy(mgs)
!        
!        ENDIF

       
       ENDDO ! k

         IF ( numshedregimes == 2 ) THEN

           dqhml(ndiam+2) = qhml(ndiam+2) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam4
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhmlr1 = qhml(ndiam+1)
           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2

         
         ELSEIF ( numshedregimes == 3 ) THEN
         
           dqhml(ndiam+3) = qhml(ndiam+3) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam4
           dqhml(ndiam+2) = qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhmlr1 = qhml(ndiam+1)
           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           qhmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         ENDIF


         IF ( .not. mixedphase ) qhmlrlg(mgs) = qhmlr1
         

        IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
 !         write(iunit,*) 'k,qhml0(n+1,+2,+3) = ',kgs(mgs),qhml0(ndiam+1),qhml0(ndiam+2),qhml0(ndiam+3)
 !         write(iunit,*) 'k,qhml(n+1,+2,+3) = ',kgs(mgs),qhml(ndiam+1),qhml(ndiam+2),qhml(ndiam+3)
 !         write(iunit,*) 'qhmlr1,qhmlr2,qhmlr12 = ',qhmlr1,qhmlr2,qhmlr12,qhml(ndiam+1)-qhml(ndiam+2)
 !         write(iunit,*) 'qhmlr,qhmlrlg = ',qhmlr(mgs),qhmlrlg(mgs)
        ENDIF


! First calc of number rate from 0 to D1
! numdiam is the number of diameter "bins" for finding the number loss at diameters less than mltdiam1 (9mm)
! If alpha > 1.5, this will be replaced with integrated value

        DO k = 1,numdiam !{
        
!        ratio = mltdiam05/xdia(mgs,lh,1)
        ratio = Min(maxratiolu, mltdiam(k)/xdia(mgs,lh,1) )

        
        IF ( usegamxinf2 ) THEN
          x =  gamxinf(2. + alpha(mgs,lh), ratio)/g1palp
        ELSE

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh))))
!           j = Int(Max(alphamin,Min(alphamax,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,2,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,2,1) - gamxinflu(i,j,2,1))
           tmp2 = gamxinflu(i,jp1,2,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,2,1) - gamxinflu(i,jp1,2,1))
           
           ! interpolate along alpha; 

          x = gaminterp(ratio,alpha(mgs,lh),2,1)

         
           x = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
        ENDIF
        ! number and mass from 0 to mltdiam05

!        IF ( usegamxinf2 .or. k == 1 ) THEN
        IF ( usegamxinf2 ) THEN
          cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)

          IF ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) THEN
           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
           tmp2 = cx(mgs,lh)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lh)*(1. - tmp)
           cxd05 = tmp2
          ENDIF
        ELSE

           ! interpolate along x, i.e., ratio; 
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
!           tmp1 = gamxinflu(i,j,1,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,1,1) - gamxinflu(i,j,1,1))
!           tmp2 = gamxinflu(i,jp1,1,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,1,1) - gamxinflu(i,jp1,1,1))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
          IF ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) THEN
          cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)
           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
           tmp2 = cx(mgs,lh)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 5. ) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
          
          ENDIF

           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)

           IF ( tmp > 0.98 ) THEN ! Interpolated values can be flaky close to 1.0, so revert to expensive gamxinf subroutine
             cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)
           ELSE
             cxd05 = cx(mgs,lh)*(1. - tmp)
           ENDIF

        ENDIF
        
         cxd(k) = cxd05


!        IF ( usegamxinf2 .or. k == 1 ) THEN
        IF ( usegamxinf2 ) THEN
          qxd05 = qx(mgs,lh)*(1. - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp)

!           tmp = gaminterp(ratio,alpha(mgs,lh),4,1)
!           tmp2 = qx(mgs,lh)*(1. - tmp)
!           tmp3 = 2.*100.*(qxd05 - tmp2)/(qxd05 + tmp2)
        
           IF ( Abs(tmp3) > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) THEN
!            write(0,*) 'k,qxd05,lu = ',k,igs(mgs),kgs(mgs),qxd05,tmp2, tmp3
!            write(0,*) 'ratio,alph,gam,interp = ', ratio,alpha(mgs,lh),1.0 - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp, 1.0 - tmp
           ENDIF
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lh)*(1. - tmp)
!          qxd05 = tmp2


        ELSE
           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,4,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,4,1) - gamxinflu(i,j,4,1))
!           tmp2 = gamxinflu(i,jp1,4,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,4,1) - gamxinflu(i,jp1,4,1))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lh),4,1)

           IF ( tmp > 0.98 ) THEN ! Interpolated values can be flaky close to 1.0, so revert to expensive gamxinf subroutine
             qxd05 = qx(mgs,lh)*(1. - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp)
           ELSE
             qxd05 = qx(mgs,lh)*(1. - tmp)
           ENDIF

        ENDIF
        
        qxd(k) = qxd05

        IF ( usegamxinf ) THEN
          y = gamxinf(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp
        ELSE
           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,3,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,3,1) - gamxinflu(i,j,3,1))
           tmp2 = gamxinflu(i,jp1,3,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,3,1) - gamxinflu(i,jp1,3,1))
           
           ! interpolate along alpha; 
           
           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lh),3,1)

            y =  tmp ! 
        ENDIF


        
        ! reuse hwvent1
        hwvent1 =  0.78*x + y*hwventy(mgs)


       qhml(k) = Min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 ) ! from diam to inf.
       qhml0(k) = qhmlr(mgs) - qhml(k)                                    ! from zero to diam
       
!       qhml0(k) = qhmlr05
! zero to 0.5*mltdiam1 is then the difference between qhmlr and qhmlr05
!       qhmlr05 = qhmlr(mgs) - qhmlr05
        IF ( k == 1 ) THEN
         dqhml(k) = qhml0(k) ! qhmlr(mgs) - qhml0(k)
        ELSE 
         dqhml(k) = qhml0(k) - qhml0(k-1)
        ENDIF
       
       ENDDO ! } k loop over diameters less than 9mm

! compute between last diameter and 9mm
       dqhml(numdiam+1) = qhml0(ndiam+1) - qhml0(numdiam)
!       cxd(numdiam+1) = 
       cxd(numdiam+1) = cxd(ndiam+1) ! - cxd(numdiam)
       qxd(numdiam+1) = qxd(ndiam+1) ! - qxd(numdiam)

! compute qx and cx for each 'bin'
! sum the number of melted particle in tmp

     dcxd(1) = cxd(1)
     dqxd(1) = qxd(1)
     tmpcmlt =  dqhml(1)*dcxd(1)/dqxd(1)
     
     IF ( .false. .and. igs(mgs) == 52 ) THEN
       write(0,*) 'HML: kgs, tmpcmlt = ',kgs(mgs),tmpcmlt,dqhml(1),dcxd(1),dqxd(1)
       write(0,*) 'dqhml+1,+2 = ',dqhml(ndiam+1),dqhml(ndiam+2)
       write(0,*) 'qhml0+1,+2 = ',qhml0(ndiam+1),qhml0(ndiam+2)
       write(0,*) 'qhml+1,+2 = ',qhml(ndiam+1),qhml(ndiam+2)
       write(0,*) 'qhmlr(mgs) = ',qhmlr(mgs)
     ENDIF
     
     DO k = 2,numdiam+1
      dcxd(k) = cxd(k) - cxd(k-1)
      dqxd(k) = qxd(k) - qxd(k-1)
!     IF ( .false. .and. igs(mgs) == 52 ) THEN
!       write(0,*) 'HML: kgs, k,delta- tmpcmlt = ',kgs(mgs),k,dqhml(k)*dcxd(k)/dqxd(k), dqhml(k),dcxd(k),dqxd(k),cxd(k),cxd(k-1)
!       write(0,*) cxd(ndiam+1), cxd(ndiam+1)-cxd(k-1)
!       write(0,*) qxd(ndiam+1), qxd(ndiam+1)-qxd(k-1)
!     ENDIF
      IF ( dqxd(k) > qxmin(lh) ) THEN
        tmpcmlt = tmpcmlt + dqhml(k)*dcxd(k)/dqxd(k)
      ENDIF
     ENDDO

     ! now, qhmlr1 and qhmlr2 lose mass but not number
     ! qhmlr0 loses both mass and number

         ! first term is melting for d < d1 divided by mean mass of particles with d < d1
         ! second and third terms are the shedding from larger particles
!         chmlrr(mgs) = qhmlr0/(qxd1/cxd1) + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)



! test integrated N rate for graupel

        IF ( alpha(mgs,lh) > 1.5 ) THEN
        
        ratio = Min( maxratiolu, mltdiam1/xdia(mgs,lh,1) )

        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami



        tmp = -1. + alpha(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf ) THEN
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(-1. + alpha(mgs,lh), ratio))/g1palp
          x = (Gamma_sp(-1. + alpha(mgs,lh)) - gamxinf(-1. + alpha(mgs,lh), ratio))/g1palp
        ELSE

          
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,7,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,7,1) - gamxinflu(i,j,7,1))
           tmp2 = gamxinflu(i,jp1,7,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,7,1) - gamxinflu(i,jp1,7,1))
           
           ! interpolate along alpha; 
           
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lh))
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lh), ratio)/g1palp
!           ENDIF


          !  x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
            
            x = tmp3  ! changed lookup table to be the difference (gam - gamxinf)/g1palp
        
        ENDIF

        tmp = -0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf ) THEN
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
           y = (Gamma_sp(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)) - gamxinf(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
        ELSE
          ! FIX THIS FOR NUMBER!

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,8,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,8,1) - gamxinflu(i,j,8,1))
           tmp2 = gamxinflu(i,jp1,8,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,8,1) - gamxinflu(i,jp1,8,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
        !   y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           
        ENDIF

        hwvent1 =  0.78*x + y*hwventy(mgs)
        
        chmlrtmp =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1/xdia(mgs,lh,1)**(2)   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
         
!          IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!           write(iunit,*) 'k,chmlrtmp',kgs(mgs),chmlrtmp,hwvent1,cx(mgs,lh),fmlt1(mgs),xdia(mgs,lh,1)
!           write(iunit,*) 'x,y,hwventy ',x,y,hwventy(mgs),bxx(mgs,lh),bx(lh)
!         ENDIF
        
        tmp = (chmlrtmp )
        
        tmpcmlt = chmlrtmp
        
!        IF ( igs(mgs) == 40 ) THEN
!          write(iunit,*) 'chmlrtmp,chmlrtmpd1inf = ',chmlrtmp
!          write(iunit,*) 'tmpcmlt old/new, alph = ',tmpcmlt,tmp,alpha(mgs,lh)
!        ENDIF
       
       ENDIF

! test using mltdiam05
!         chmlrr(mgs) = (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05) + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)
!         tmpcmlt = (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05)
         IF ( numshedregimes == 2 ) THEN
           chmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)

         ELSEIF ( numshedregimes == 3 ) THEN
!           qhmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
!           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
!           qhmlr1 = qhml(ndiam+1)
!           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
!           qhmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           chmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhmlr3*mltmass1inv + qhmlr23*mltmass2inv + qhmlr12*mltmass3inv)
         ENDIF

!         IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!           write(iunit,*) 'k,qh,alp,qhmlr,chmlrr(mgs),xdia ',kgs(mgs),1.e3*qx(mgs,lh),alpha(mgs,lh),qhmlr(mgs),chmlrr(mgs),xdia(mgs,lh,3)
!           write(iunit,*) 'zx,tmpcmlt,qhmlr2,qhmlr12,xdn ',zx(mgs,lh),tmpcmlt,qhmlr2,qhmlr12,xdn(mgs,lh),qhmlr12*mltmass2inv,qhmlr2*mltmass1inv
!         ENDIF

         IF ( .true. .or. alpha(mgs,lh) < alphamax - delta_alphamlr .or. xdia(mgs,lh,3) > binmlrmxdia ) THEN ! { allow spectrum to get narrower
     !     chmlr(mgs)  = qhmlr0/(rho0(mgs)*qxd1/cxd1) 
          chmlr(mgs)  =  tmpcmlt ! (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05)
          

! test integrated Z rate for graupel

        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hwventy(mgs)
!        hwvent1 =    &
!     &  ( 0.78*x +    &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       IF ( .false. ) THEN
       ! full shedding version
       zhmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6.*2./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))*2/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 24
         ! factor of rho0 to nullify the 1/rho0 in fmlt1

! alternate zhmlr assumes Z = N D**6, from which dZ/dt = (dN/dt)*D**6 = (1/M) dM/dt D**6
! turns out to be exactly a factor of 2 smaller than the other zhmlr (and matches with bin results when set to remove full particles
! from each bin, i.e., particles do not shrink)
       ELSE
       ! whole particle loss version: factor of 1/2 different from shedding version
       zhmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
       ENDIF
       
       zhmlr0inf = zhmlr(mgs)

! note that fmlt1 has the factor of 2*pi
!        fmlt1(mgs) = (2.0*pi)*   &
!     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))   &
!     &   -ftka(mgs)*temcg(mgs)/rho0(mgs) )    &
!     &  / (felf(mgs))


          IF ( .false. .and. igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
          tmp1 = zhmlr0inf
          
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
          zhmlrtmp =  g1*rho0(mgs)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs)   - tmp**2 * qhmlr(mgs)/tmp  )
          
          write(iunit,*) 'zhmlrtmp,zhmlr,1.5*zhmlr: ',zhmlrtmp,zhmlr(mgs),1.5*zhmlr(mgs),zhmlrtmp/(zhmlr(mgs)),alpha(mgs,lh)
          
         ! try average:
          
    !      zhmlr(mgs) = 0.5*(zhmlr(mgs) + zhmlrtmp)
          
          ENDIF

!        IF ( .false. .and. igs(mgs) == 40 ) THEN !.and. zhmlr(mgs) < -zxmin  THEN
! #include 'bintest.F90'  
!        ENDIF


          !  try using the new N and q with the old Z to compute a new alpha, then use new alpha (assuming it is larger)
          ! to compute delta-Z?

          il = lh
          chw = cx(mgs,il) + chmlr(mgs)
          qr  = qx(mgs,il) + qhmlr(mgs)
          z   = zx(mgs,il)  ! + zhmlr(mgs)

          IF ( .false. .and. z .gt. zxmin .and. qr > qxmin(lh) .and. chw > cxmin ) THEN
           
           il = lh
           
           tmp2 = alpha(mgs,lh) ! save alpha
           
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = Max( alphamin, Min( alphamax, alp ) )
           
         IF ( newton ) THEN
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = alp + ( galpha(alp) - rdi )/dgalpha(alp)
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
           
         ELSE
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
          ENDIF
          
            
            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            
            tmp = qx(mgs,lh)/cx(mgs,lh)
            
            zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
            
            alpha(mgs,lh) = tmp2 ! restore alpha
            
          ENDIF ! everything greater than min values

      ENDIF !}
!
! not sure this is right -- should it use g1x(mgs,lh) or some assumed value for melt rain?
!
! should separately add the reflecitivities from the 2 sizes of shed drops, then add the completely melted particles

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )


! test integrated Z rate for rain for size range up to mltdiam1
! integrate from d1 to infinity, then subtract that from the full value (zhmlr)
        ratio = Min( maxratiolu, mltdiam1/xdia(mgs,lh,1) )


        tmp = 1. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5. + alpha(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf ) THEN
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lh), ratio))/g1palp
          x = (Gamma_sp(5. + alpha(mgs,lh)) - gamxinf(5. + alpha(mgs,lh), ratio))/g1palp
        ELSE

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,5,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,5,1) - gamxinflu(i,j,5,1))
           tmp2 = gamxinflu(i,jp1,5,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,5,1) - gamxinflu(i,jp1,5,1))
           
           ! interpolate along alpha; 
           
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lh))
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
           tmp3 = gaminterp(ratio,alpha(mgs,lh),5,1)
           
!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lh), ratio)/g1palp
!           ENDIF


!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - tmp3)/g1palp
!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)*(1. - tmp3)
       ! x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
         x = tmp3 ! changed lookup table to be the difference (gam - gamxinf)/g1palp
       !  x = Max(0.0,x)
         
        ENDIF

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf ) THEN
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))  - gamxinf(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
        ELSE
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,1) - gamxinflu(i,j,6,1))
           tmp2 = gamxinflu(i,jp1,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,1) - gamxinflu(i,jp1,6,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
           
           tmp3 = gaminterp(ratio,alpha(mgs,lh),6,1)
           
         !  y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3 ! changed lookup table to be the difference (gam - gamxinf)/g1palp
           y = Max(0.0, y)
        ENDIF
        IF ( .false. .and. igs(mgs) == 38 .and. qhmlr(mgs) < 0.0 ) THEN
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,1) - gamxinflu(i,j,6,1))
           tmp2 = gamxinflu(i,jp1,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,1) - gamxinflu(i,jp1,6,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1)) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
          write(iunit,*) 'ibinhmlr = 1: i,k = ',igs(mgs),kgs(mgs)
          write(iunit,*) 'y,ratio,alpha,mltdiam1,xdia(mgs,lh,1) = ',y,ratio,alpha(mgs,lh),mltdiam1,xdia(mgs,lh,1)
          write(iunit,*) 'full gamx = ',gamxinf(5.5 +0.5*bxx(mgs,lh) + alpha(mgs,lh), ratio),bxx(mgs,lh),bx(lh)
           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
          write(iunit,*) 'interpolated y = ', y,tmp1,tmp2,tmp3, &
               gamxinflu(i,j,6,1),gamxinflu(ip1,j,6,1),gamxinflu(i,jp1,6,1),gamxinflu(ip1,jp1,6,1)
           
        ENDIF

        hwvent1 =  0.78*x + y*hwventy(mgs)
       
       zhmlrr(mgs) =   &  !  -- THIS IS THE LOWER INCOMPLETE GAMMA... 0 to ratio
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))/pi comes from removing the mass/volume factors; 

! zhmlrr here is whole particle loss for d < mltdiam1
!
! so now have zhmlrr as the Z loss (whole particle) for graupel for 0 <= d <= mltdiam1, and zhmlr0inf as the integrated whole-particle loss
! We want whole particle loss for 0 < d < mltdiam1 which is (zhmlrr(mgs)), plus shedding loss for d > mltdiam1, which is 2*(zhmlr0inf - zhmlrr(mgs)) 
! Shedding loss is simply double the whole particle loss

       zhmlr(mgs) = zhmlrr(mgs) + binmlrzrrfac*2.*(zhmlr0inf - zhmlrr(mgs)) ! factor of 2 to convert whole particle loss to shedding loss
       ! should that be zhmlr(mgs) = 2.*zhmlrr(mgs) + (zhmlr0inf - zhmlrr(mgs)) = zhmlrr(mgs) + zhmlr0inf
       ! yes, I think it should -- NOPE!
       ! zhmlr(mgs) = 2.*zhmlrr(mgs) + (zhmlr0inf - zhmlrr(mgs))

          IF ( .false. .and. igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
          
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
          zhmlrtmp =  g1*rho0(mgs)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs)   - tmp**2 * qhmlr(mgs)/tmp  )

        IF ( alpha(mgs,lh) >= alphamax - delta_alphamlr ) THEN
        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hwventy(mgs)

       ! whole particle loss version: factor of 1/2 different from shedding version
 !      zhmlr0inf =   &
 !    &   min(   rho0(mgs)* &
 !    &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
 !    &   , 0.0 )/xdn(mgs,lh)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
        
        tmp = zhmlrr(mgs) + 2.*(zhmlr0inf - zhmlrr(mgs)) 
        
       ENDIF
          
!          write(iunit,*) 'zhmlrtmp,zhmlr,zhmlr0inf: ',zhmlrtmp,zhmlr(mgs),zhmlr0inf,zhmlrtmp/(zhmlr(mgs)),alpha(mgs,lh),tmp,zhmlrtmp/(zhmlr0inf)
                    
          ENDIF

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr0 - tmp**2 * chmlr(mgs)  ) )

         IF ( numshedregimes == 2 ) THEN
           zhmlrr(mgs) = (xdn(mgs,lh)/1000.)**2*zhmlrr(mgs) + rho0(mgs)*   &
     &                   (qhmlr2*mltmass1inv*(0.01*takshedsize1)**6 + qhmlr12*mltmass2inv*(0.01*takshedsize2)**6 )
         ELSEIF ( numshedregimes == 3 ) THEN
!           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
          zhmlrr(mgs) = (xdn(mgs,lh)/1000.)**2 *zhmlrr(mgs) + rho0(mgs)*  &
     &              (qhmlr3*mltmass1inv*(0.01*takshedsize1)**6 + qhmlr23*mltmass2inv*(0.01*takshedsize2)**6  + &
     &               qhmlr12*mltmass3inv*(0.01*takshedsize3)**6 )
         ENDIF


         IF ( alpha(mgs,lh) >= alphamax - delta_alphamlr ) THEN ! { cannot allow spectrum to get any narrower, 
         
          chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
         
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
          zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )

!          chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh)))  ! into rain 

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )

!          IF(imltshddmr == 1 ) THEN
!            ! DTD: If Dmg < sheddiam, then assume complete melting into
!            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
!            tmp = -rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh))) ! Min of Maximum raindrop size/mean hail size
!            tmp2 = -rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
!            chmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)
!            chmlrr(mgs) = -Max(tmp,Min(tmp2,chmlrr(mgs)))
!          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
!            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
!            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
!            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
!          ENDIF

         ELSEIF (  xdia(mgs,lh,3) > binmlrmxdia ) THEN ! 
         
          chmlrtmp  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrtmp )
          
          tmp1 = 1.e-3
          IF ( xdia(mgs,lh,3) > binmlrmxdia + tmp1 ) THEN 
             chmlr(mgs) = chmlrtmp
             zhmlr(mgs) = zhmlrtmp
          ELSE ! weighted average
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'k,xdia,chmlr,zhmlr = ',kgs(mgs),xdia(mgs,lh,3)*1.e3,chmlr(mgs),zhmlr(mgs)
!              write(0,*) 'weights 1,2 = ',(xdia(mgs,lh,3) - binmlrmxdia )/tmp1 , (binmlrmxdia + tmp1 - xdia(mgs,lh,3))/tmp1
!            ENDIF
            chmlr(mgs) = chmlrtmp*(xdia(mgs,lh,3) - binmlrmxdia )/tmp1 + chmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lh,3))/tmp1
            zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
            !zhlmlr(mgs) = zhlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + zhlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'tmp,chlmlr,= ',chlmlrtmp,chlmlr(mgs)
!              write(0,*) 'tmp,zhlmlr,= ',zhlmlrtmp,zhlmlr(mgs)
!            ENDIF
          ENDIF
          

         ENDIF ! }

           
!         IF ( igs(mgs) == 38 ) THEN
!           x = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
!           write(iunit,*) 'k,chmlrr: ',kgs(mgs),chmlrr(mgs),qx(mgs,lh),qhmlr(mgs),zhmlrr(mgs),qxd05,qxd1,cxd05,cxd1
!!           write(iunit,*) 'k,zhmlr: ',kgs(mgs),zhmlr(mgs),x,qx(mgs,lh),qhmlr(mgs),zhmlrr(mgs),xdia(mgs,lh,1),zx(mgs,lh)
!!           write(iunit,*) 'zhmlrr parts: ',  2.*tmp * qhmlr(mgs), - tmp**2 * chmlrr(mgs) , 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs) 
!           write(iunit,*) 'chmlrr parts: ',(qhmlr0-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)), qhmlr05/(qxd05/cxd05), rho0(mgs)*(qhmlr2*mltmass1inv), rho0(mgs)*( qhmlr12*mltmass2inv)
!           write(iunit,*) 'qhmlr0: ',qhmlr0,qhmlr05,qhmlr2,qhmlr12
!         ENDIF
         
       ELSEIF ( ibinhmlr == 2 .or. ibinhmlr == 3 ) THEN

       ENDIF
       
       
       IF ( ivhmltsoak > 0 .and. qhmlr(mgs) < 0.0 .and. lvol(lh) > 1 .and. xdn(mgs,lh) .lt. xdnmx(lh) ) THEN
         ! act as if 100% of the meltwater were soaked into the graupel
           v1 = (1. - xdn(mgs,lh)/xdnmx(lh))*(vx(mgs,lh) + rho0(mgs)*qhmlr(mgs)/xdn(mgs,lh) )/(dtp) ! volume available for filling
           v2 = -1.0*rho0(mgs)*qhmlr(mgs)/xdnmx(lh)  ! volume of melted ice if it were refrozen in the matrix
           
           vhsoak(mgs) = Min(v1,v2)
           
       ENDIF

      ENDIF !  qx(mgs,lh) .gt. qxmin(lh)

      
      IF ( lhl .gt. 1  .and. lhlw < 1 ) THEN

       IF ( qx(mgs,lhl) .gt. qxmin(lhl) ) THEN
         IF ( ibinhlmlr == 0  .or. lzhl < 1) THEN
       qhlmlr(mgs) =   &
     &   meltfac*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hlvent(mgs)*xdia(mgs,lhl,1)   &
     &  + fmlt2(mgs)*(qhlacrmlr(mgs)+qhlacwmlr(mgs))    &
     &   , 0.0 )

       ELSEIF ( ibinhlmlr == 1 ) THEN ! use incomplete gamma functions to approximate the bin results


       qhlmlr(mgs) =   &
     &   min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hlvent(mgs)*xdia(mgs,lhl,1)   &
!     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacw(mgs))    & ! turn off the collection part for now
     &   , 0.0 )

! first part to integrate from mltdiam1 and mltdiam2 to infinity (k loop)

        tmp = 1. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g4palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
 
       DO k = 1,numshedregimes  ! do not need k=3, treat as infinity; goes to 3 to cover the top 3 diameters (independent of ndiam, the number of smaller diameters)
        
!        ratio = mltdiam1/xdia(mgs,lhl,1)
        ratio = Min(maxratiolu, mltdiam(ndiam+k)/xdia(mgs,lhl,1) )

        
        IF ( usegamxinf3 ) THEN
        x =  gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
        ELSE

          x = gaminterp(ratio,alpha(mgs,lhl),2,2)
          
          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck01: i,k,xint,x = ',igs(mgs),kgs(mgs),x,tmpgam
            ENDIF
          ENDIF

       ENDIF ! usegamxinf

! qxd(ndiam+2), cxd(ndiam+2), qhml(ndiam+2)        
        IF ( usegamxinf3 ) THEN
         cxd1 = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
        ELSE

          tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)

          
!          IF ( tmp > 0.95 ) THEN
!           tmp = gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp
!          ENDIF

          cxd1 = cx(mgs,lhl)*(1. - tmp)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
              write(6,*) 'gamcheck02: i,k,cxd1,x = ',igs(mgs),kgs(mgs),cxd1,tmpgam
            ENDIF
          ENDIF


        ENDIF

         cxd(ndiam+k) = cxd1

        IF ( usegamxinf3 ) THEN
         qxd1 = qx(mgs,lhl)*(1. - gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp)
        ELSE

          tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)

!           IF ( tmp > 0.95 ) THEN
!            tmp = gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp
!           ENDIF

           qxd1 = qx(mgs,lhl)*(1. - tmp)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = qx(mgs,lhl)*(1. - gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp)
             ! tmp1 =  qx(mgs,lhl)*(1. - gaminterp(ratio,alpha(mgs,lhl),10,2))
              write(6,*) 'gamcheck03: i,k,qxd1,x = ',igs(mgs),kgs(mgs),qxd1,tmpgam !,tmp1
             IF ( .false. ) THEN
              write(6,*) 'alpha,gaminterp = ',alpha(mgs,lhl),tmp
              write(6,*) 'ratio,ratio2,mltdiam,xdia = ',ratio,mltdiam(ndiam+k)/xdia(mgs,lhl,1),mltdiam(ndiam+k),xdia(mgs,lhl,1) 
              tmp1=gamxinfdp(4. + alpha(mgs,lhl), ratio)
              tmp2=gamma_dpr(4. + alpha(mgs,lhl))
              write(6,*) 'gam,gam4,gam4alt,gam/gam4 = ',tmp1,g4palp,tmp2,tmp1/g4palp
           alp = alpha(mgs,lhl)
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           ENDIF
           delx = Min(maxratiolu,ratio) - float(i)*dqiacrratio
           dely = Min(maxalphalu, alp ) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )
           
           luindex = 4
           il = 2

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,luindex,il) + delx*dqiacrratioinv*         &
     &                 (gamxinflu(ip1,j,luindex,il) - gamxinflu(i,j,luindex,il))
           tmp2 = gamxinflu(i,jp1,luindex,il) + delx*dqiacrratioinv*       &
     &                 (gamxinflu(ip1,jp1,luindex,il) - gamxinflu(i,jp1,luindex,il))
           
           ! interpolate along alpha; 
           
           temp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
              write(6,*) 'i,j,temp3,gam1,',i,j,temp3, gamxinflu(i,j,luindex,il)
            ENDIF ! false/true
            
            ENDIF
          ENDIF

        ENDIF

         qxd(ndiam+k) = qxd1

        IF ( usegamxinf3 ) THEN
          y =  gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
        ELSE

          y = gaminterp(ratio,alpha(mgs,lhl),3,2)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck04: i,k,yint,y = ',igs(mgs),kgs(mgs),y,tmpgam
            ENDIF
          ENDIF

        ENDIF
        
          hwvent1 =  0.78*x + y*hlventy(mgs) 

          qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 )

        IF ( igs(mgs) == 40 ) THEN
!          write(0,*) 'kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
!          write(0,*) 'diam loop k,qhmlr1 = ',k,qhmlr1,qxd1,cxd1
!          write(0,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,y,hwventy(mgs)
!          write(0,*) 'ratio0,ratio,xdia(1) = ',mltdiam(ndiam+k)/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1)
        
        ENDIF

! zero to mltdiam1 is then the difference between qhlmlr and qhlmlr1
       qhlmlr0 = qhlmlr(mgs) - qhlmlr1
       qhml(ndiam+k) = qhlmlr1 ! from diam to infinity
       qhml0(ndiam+k) = qhlmlr0 ! from zero to diam
        
!        IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!          write(iunit,*) 'diam loop k,qhlmlr1 = ',k,qhlmlr1,qxd1,cxd1
!          write(iunit,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,tmp3,y,hlventy(mgs)
!          write(iunit,*) 'lookup: ',i,j,gamxinflu(i,j,2,1),gamxinflu(i,j,2,2) !, gaminterp(ratio,alpha(mgs,lhl),2,-2)
!        
!        ENDIF
       
       ENDDO ! k = 1,2

         IF ( numshedregimes == 2 ) THEN
           dqhml(ndiam+2) = qhml(ndiam+2) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhlmlr1 = qhml(ndiam+1)
           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         ELSEIF ( numshedregimes == 3 ) THEN
         
           dqhml(ndiam+3) = qhml(ndiam+3) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+2) = qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhlmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhlmlr1 = qhml(ndiam+1)
           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           qhlmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         ENDIF

         IF ( .not. mixedphase ) qhlmlrlg(mgs) = qhlmlr1
         

!  integrate from mltdiam1/2 to infinity

        DO k = 1,numdiam
        
!        ratio = mltdiam05/xdia(mgs,lhl,1)
        ratio = Min( maxratiolu,  mltdiam(k)/xdia(mgs,lhl,1) )

        
        IF ( usegamxinf3 ) THEN
          x =  gamxinf(2. + alpha(mgs,lhl), ratio)/g1palp
        ELSE

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl))))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,2,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,2,2) - gamxinflu(i,j,2,2))
           tmp2 = gamxinflu(i,jp1,2,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,2,2) - gamxinflu(i,jp1,2,2))
           
           ! interpolate along alpha; 

         !  x = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

          x = gaminterp(ratio,alpha(mgs,lhl),2,2)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinfdp(2.0 + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck05: i,k,xint,x = ',igs(mgs),kgs(mgs),x,tmpgam
            ENDIF
          ENDIF

         
        ENDIF
        ! number and mass from 0 to mltdiam05

!        IF ( usegamxinf2 .or. k == 1 ) THEN
        IF ( usegamxinf3 ) THEN
          cxd05 = cx(mgs,lhl)*(1. - gamxinf(1. + alpha(mgs,lhl), ratio)/g1palp)

          IF ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) THEN
           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)
           tmp2 = cx(mgs,lhl)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
       !    IF ( tmp3 > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lhl)*(1. - tmp)
          cxd05 = tmp2
          ENDIF
        ELSE

           ! interpolate along x, i.e., ratio; 
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           tmp1 = gamxinflu(i,j,1,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,1,2) - gamxinflu(i,j,1,2))
!           tmp2 = gamxinflu(i,jp1,1,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,1,2) - gamxinflu(i,jp1,1,2))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
          IF ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) THEN
          cxd05 = cx(mgs,lhl)*(1. - gamxinf(1. + alpha(mgs,lhl), ratio)/g1palp)
           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)
           tmp2 = cx(mgs,lhl)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 5. ) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
          
          ENDIF

           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)

           IF ( tmp > 0.98 ) THEN
             cxd05 = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
           ELSE
             cxd05 = cx(mgs,lhl)*(1. - tmp)
           ENDIF

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinfdp(1.0 + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck06: i,k,tmp,x = ',igs(mgs),kgs(mgs),tmp,tmpgam
              tmp1 =  cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
              write(6,*) 'cxd05,gam: ',cxd05,tmp1
            ENDIF
          ENDIF
          
          cxd05 = Max(0.0, cxd05)

        ENDIF
        
         cxd(k) = cxd05


!        IF ( usegamxinf2 .or. k == 1 ) THEN
        IF ( usegamxinf3 ) THEN
          qxd05 = qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)

!           tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)
!           tmp2 = qx(mgs,lhl)*(1. - tmp)
!           tmp3 = 2.*100.*(qxd05 - tmp2)/(qxd05 + tmp2)
        
           IF ( Abs(tmp3) > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) THEN
!            write(0,*) 'k,qxd05,lu = ',k,igs(mgs),kgs(mgs),qxd05,tmp2, tmp3
!            write(0,*) 'ratio,alph,gam,interp = ', ratio,alpha(mgs,lhl),1.0 - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp, 1.0 - tmp
           ENDIF
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lhl)*(1. - tmp)
!          qxd05 = tmp2


        ELSE
           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,4,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,4,2) - gamxinflu(i,j,4,2))
!           tmp2 = gamxinflu(i,jp1,4,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,4,2) - gamxinflu(i,jp1,4,2))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)

           IF ( tmp > 0.98 ) THEN
             qxd05 = qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
           ELSE
             qxd05 = qx(mgs,lhl)*(1. - tmp)
           ENDIF

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp
              write(6,*) 'gamcheck07: i,k,tmp,x = ',igs(mgs),kgs(mgs),tmp,tmpgam
              tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
              write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF


        ENDIF
        
        qxd(k) = qxd05

        IF ( usegamxinf3 ) THEN
          y = gamxinf(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
        ELSE
           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,3,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,3,2) - gamxinflu(i,j,3,2))
           tmp2 = gamxinflu(i,jp1,3,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,3,2) - gamxinflu(i,jp1,3,2))
           
           ! interpolate along alpha; 
           
           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lhl),3,2)

            y =  tmp ! 

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = gamxinf(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck08: i,k,tmp,y = ',igs(mgs),kgs(mgs),tmp,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF

        ENDIF


        
        ! reuse hwvent1
        hwvent1 =  0.78*x + y*hlventy(mgs)


       qhml(k) = Min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 ) ! from diam to inf.
       qhml0(k) = qhlmlr(mgs) - qhml(k)                                    ! from zero to diam
 
 
        IF ( igs(mgs) == 40 ) THEN
!          write(0,*) '2 kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
!          write(0,*) '2 diam loop k,qhml,qhml0 = ',k,qhml(k),qhml0(k) ,qxd05,cxd05
!          write(0,*) '2 fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,y,hwventy(mgs)
!          write(0,*) '2 ratio0,ratio,xdia(1) = ',mltdiam(k)/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1),alpha(mgs,lhl)
        
        ENDIF
      
!       qhml0(k) = qhlmlr05
! zero to 0.5*mltdiam1 is then the difference between qhlmlr and qhlmlr05
!       qhlmlr05 = qhlmlr(mgs) - qhlmlr05
        IF ( k == 1 ) THEN
         dqhml(k) = qhml0(k) ! qhlmlr(mgs) - qhml0(k)
        ELSE 
         dqhml(k) = qhml0(k) - qhml0(k-1)
        ENDIF
       
       ENDDO ! k loop over diameters less than 9mm

! compute between last diameter and 9mm
       dqhml(numdiam+1) = qhml0(ndiam+1) - qhml0(numdiam)
!       cxd(numdiam+1) = 
       cxd(numdiam+1) = cxd(ndiam+1) ! - cxd(numdiam)
       qxd(numdiam+1) = qxd(ndiam+1) ! - qxd(numdiam)

! compute qx and cx for each 'bin'
! sum the number of melted particle in tmp

     dcxd(1) = cxd(1)
     dqxd(1) = qxd(1)
     tmpcmlt =  dqhml(1)*dcxd(1)/dqxd(1)
     
     IF ( .false. .and. igs(mgs) == 52 ) THEN
       write(0,*) 'HML: kgs, tmpcmlt = ',kgs(mgs),tmpcmlt,dqhml(1),dcxd(1),dqxd(1)
       write(0,*) 'dqhml+1,+2 = ',dqhml(ndiam+1),dqhml(ndiam+2)
       write(0,*) 'qhml0+1,+2 = ',qhml0(ndiam+1),qhml0(ndiam+2)
       write(0,*) 'qhml+1,+2 = ',qhml(ndiam+1),qhml(ndiam+2)
       write(0,*) 'qhlmlr(mgs) = ',qhlmlr(mgs)
     ENDIF
     
     DO k = 2,numdiam+1
      dcxd(k) = cxd(k) - cxd(k-1)
      dqxd(k) = qxd(k) - qxd(k-1)
     IF ( .false. .and. igs(mgs) == 52 ) THEN
       write(0,*) 'HML: kgs, k,delta- tmpcmlt = ',kgs(mgs),k,dqhml(k)*dcxd(k)/dqxd(k), dqhml(k),dcxd(k),dqxd(k),cxd(k),cxd(k-1)
       write(0,*) cxd(ndiam+1), cxd(ndiam+1)-cxd(k-1)
       write(0,*) qxd(ndiam+1), qxd(ndiam+1)-qxd(k-1)
     ENDIF
      IF ( dqxd(k) > qxmin(lhl) ) THEN
        tmpcmlt = tmpcmlt + dqhml(k)*dcxd(k)/dqxd(k) ! this value of tmpcmlt is used if alpha(lhl) <= 1.5
      ENDIF
     ENDDO

     ! now, qhlmlr1 and qhlmlr2 lose mass but not number
     ! qhlmlr0 loses both mass and number

         ! first term is melting for d < d1 divided by mean mass of particles with d < d1
         ! second and third terms are the shedding from larger particles
!         chlmlrr(mgs) = qhlmlr0/(qxd1/cxd1) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)

! test using mltdiam05
!         chlmlrr(mgs) = (qhml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
!         tmp = (qhml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
!         chlmlrr(mgs) = tmp + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)



! test integrated N rate for Hail

        IF ( alpha(mgs,lhl) > 1.5 ) THEN ! For lower alpha, the estimate of tmpcmlt from above is used
        
        ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )

        tmp = 1.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = -1. + alpha(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf3 ) THEN
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
          x = (Gamma_sp(-1. + alpha(mgs,lhl)) - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
        ELSE

          
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,7,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,7,2) - gamxinflu(i,j,7,2))
           tmp2 = gamxinflu(i,jp1,7,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,7,2) - gamxinflu(i,jp1,7,2))
           
           ! interpolate along alpha; 
           
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
!           IF ( tmp3 > 0.98 ) THEN ! use direct calculation in this nonlinear region
!            tmp3 = gamxinf(5. + alpha(mgs,lhl), ratio)/g1palp
!           ENDIF

          !  x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
            x = tmp3 ! changed lookup table to be the difference already:  Gamma_sp(-1+alp) - Gamma_sp(-1+alp,ratio)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = (Gamma_sp(-1. + alpha(mgs,lhl)) - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
              tmp1 = gaminterp(ratio,alpha(mgs,lhl),7,2)
              write(6,*) 'gamcheck09: i,k,tmp3,interp,x = ',igs(mgs),kgs(mgs),tmp3,tmp1,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF
        
        ENDIF

        tmp = -0.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf3 ) THEN
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(-0.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp( alpha(mgs,lhl) - 0.5 + 0.5*bxx(mgs,lhl) )  - gamxinf( alpha(mgs,lhl) - 0.5 + 0.5*bxx(mgs,lhl), ratio))/g1palp
        ELSE

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,8,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,8,2) - gamxinflu(i,j,8,2))
           tmp2 = gamxinflu(i,jp1,8,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,8,2) - gamxinflu(i,jp1,8,2))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lh))
           
        !   y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           
          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = (Gamma_sp(-0.5 + alpha(mgs,lhl)) - gamxinf(-0.5 + alpha(mgs,lhl), ratio))/g1palp
              tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck10: i,k,tmp3,interp,x = ',igs(mgs),kgs(mgs),tmp3,tmp1,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF
        ENDIF

        hwvent1 =  0.78*x + y*hlventy(mgs)
        
        chlmlrtmp =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1/xdia(mgs,lhl,1)**(2)   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
         
!          IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!           write(iunit,*) 'k,chlmlrtmp',kgs(mgs),chlmlrtmp,hwvent1,cx(mgs,lhl),fmlt1(mgs),xdia(mgs,lhl,1)
!           write(iunit,*) 'x,y,hlventy ',x,y,hlventy(mgs),bxx(mgs,lhl),bx(lhl)
!         ENDIF
       
        tmp = (chlmlrtmp )
        
        tmpcmlt = chlmlrtmp
        
!        IF ( igs(mgs) == 40 ) THEN
!          write(iunit,*) 'chlmlrtmp,chlmlrtmpd1inf = ',chlmlrtmp
!          write(iunit,*) 'tmpcmlt old/new, alph = ',tmpcmlt,tmp,alpha(mgs,lhl)
!        ENDIF
       
       ENDIF

! test using mltdiam05
!         chlmlrr(mgs) = (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
!         tmpcmlt = (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
         IF ( numshedregimes == 2 ) THEN
           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
         ELSEIF ( numshedregimes == 3 ) THEN
!           qhlmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
!           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
!           qhlmlr1 = qhml(ndiam+1)
!           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
!           qhlmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
         ENDIF
         
         ! hack!
         !  chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)*mltmass0inv
         !   IF ( igs(mgs) == nx/2 ) THEN
         !     write(6,*) 'i,k,tmpcmlt,othercmlt = ',igs(mgs),kgs(mgs),tmpcmlt,rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
         !     write(6,*) 'mlr0,1,2,3: ', qhlmlr(mgs)-qhlmlr1,qhlmlr12,qhlmlr23,qhlmlr3  !
         !     tmp = (6.*rho0(mgs)*(qhlmlr(mgs)-qhlmlr1)/(xdn0(lr)*tmpcmlt*pi) )**(1./3.)
         !     tmp1 = (6.*rho0(mgs)*(qhlmlr(mgs))/(xdn0(lr)*chlmlrr(mgs)*pi) )**(1./3.)
         !     write(6,*) 'diam0, diamave = ',tmp*1000.,tmp1*1000.
         !   ENDIF
         
!         IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!           write(iunit,*) 'k,qhl,alp,qhlmlr,chlmlrr(mgs),xdia ',kgs(mgs),1.e3*qx(mgs,lhl),alpha(mgs,lhl),qhlmlr(mgs),chlmlrr(mgs),xdia(mgs,lhl,3)
!           write(iunit,*) 'zx,tmpcmlt,qhlmlr2,qhlmlr12,xdn',zx(mgs,lhl),tmpcmlt,qhlmlr2,qhlmlr12,xdn(mgs,lhl),qhlmlr12*mltmass2inv,qhlmlr2*mltmass1inv
!         ENDIF

         IF ( .true. .or. alpha(mgs,lhl) < alphamax - delta_alphamlr ) THEN ! { allow spectrum to get narrower, BUT still want to calculate chlmlrr
     !     chlmlr(mgs)  = qhlmlr0/(rho0(mgs)*qxd1/cxd1) 
          chlmlr(mgs)  =  tmpcmlt ! (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
          
          ! check for negative value of cx, and, if so revert to typical number loss rate.
          IF ( cx(mgs,lhl) + dtp*chlmlr(mgs) < 0.0 ) THEN
             chlmlr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
          ENDIF
          
! test integrated Z rate for hail

        tmp = 1.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hlventy(mgs)
       
       ! whole particle loss version: factor of 1/2 different from shedding version -- THIS IS THE COMPLETE GAMMA... 0 to infinity
       zhlmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1

       zhlmlr0inf = zhlmlr(mgs)

      ENDIF !}


!
! not sure this is right -- should it use g1x(mgs,lhl) or some assumed value for melt rain?
!
! should separately add the reflecitivities from the 2 sizes of shed drops, then add the completely melted particles

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )


! test integrated Z rate for rain for size range up to mltdiam1
! integrate from d1 to infinity, then subtract that from the full value (zhlmlr)
        ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )


        tmp = 1. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5. + alpha(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        IF ( usegamxinf3 ) THEN
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          x = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
        ELSE

!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           delx = ratio - float(i)*dqiacrratio
!           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
!           ip1 = Min( i+1, nqiacrratio )
!           jp1 = Min( j+1, nqiacralpha )
!
!           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,5,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,5,2) - gamxinflu(i,j,5,2))
!           tmp2 = gamxinflu(i,jp1,5,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,5,2) - gamxinflu(i,jp1,5,2))
!           
!           ! interpolate along alpha; 
!           
!!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lhl))
!           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
!           tmp5 = tmp3

           tmp3 = gaminterp(ratio,alpha(mgs,lhl),5,2)

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
            !  tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck11: i,k,tmp3,x = ',igs(mgs),kgs(mgs),tmp3,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF

!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lhl), ratio)/g1palp
!           ENDIF


!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - tmp3)/g1palp
!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)*(1. - tmp3)
   !     x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
           x = tmp3
           x2 = x

          IF ( x <= 0.0 ) THEN ! check for unphysical negative result from interpolation
       !    x = 0.0
    !       x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
           x = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          ELSE
    !       x2 = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
    !       x2 = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          ENDIF

          IF ( .false. .and. igs(mgs) == 40 ) THEN
            write(0,*) '3 x: tmp3,tmp5 = ',tmp3,tmp5
          ENDIF

        ENDIF

        tmp = 5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        ! problem with interpolation of this one....
        IF ( usegamxinf3 ) THEN
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
        ELSE
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           delx = ratio - float(i)*dqiacrratio
!           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
!           ip1 = Min( i+1, nqiacrratio )
!           jp1 = Min( j+1, nqiacralpha )
!
!           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,2) - gamxinflu(i,j,6,2))
!           tmp2 = gamxinflu(i,jp1,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,2) - gamxinflu(i,jp1,6,2))
!           
!           ! interpolate along alpha; 
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))
!
!           tmp4 = tmp3
           
!           IF ( igs(mgs) == 40 ) THEN
!             tmp3 = gaminterp(ratio,alpha(mgs,lhl),6,-2)
!           ELSE
             tmp3 = gaminterp(ratio,alpha(mgs,lhl),6,2)
!           ENDIF

          IF ( gammacheck ) THEN
            IF ( igs(mgs) == nx/2 ) THEN
              tmpgam = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - &
                        gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
            !  tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck12: i,k,tmp3,y = ',igs(mgs),kgs(mgs),tmp3,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            ENDIF
          ENDIF
           
      !     y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           y2 = y

          IF ( y <= 0.0 ) THEN ! check for unphysical negative result from interpolation
          ! y = 0.0
          ! y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          ELSE
!           y2 = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
!           y2 = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          ENDIF
           
        IF ( .false. .and. igs(mgs) == nx/2 ) THEN
          write(0,*) '3 kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
          write(0,*) '3 y-interp, tmp = ',y,tmp
          tmp = (gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          write(0,*) '3 tmp3,tmp4 gamx = ',tmp3,tmp4,tmp
!          call gaminterpsub(tmp3,ratio,alpha(mgs,lhl),6,2)
          write(0,*) 'tmp3 from subroutine: ',tmp3
          write(0,*) 'i,j,jp1,delx,dely,gam1,gamjp1 = ',i,j,jp1,delx,dely,gamxinflu(i,j,6,2),gamxinflu(i,jp1,6,2) 
          write(0,*) '3 term1,2 = ', (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami),  &
                       gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)
          write(0,*) '3 term1,2/g1palp = ', (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp, &
                       gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
          write(0,*) '3 ratio0,ratio,xdia(1) = ',mltdiam1/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1),alpha(mgs,lhl)
        
        ENDIF
           
        ENDIF
        IF ( .false. .and. igs(mgs) == 38 .and. qhlmlr(mgs) < 0.0 ) THEN
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           IF ( alp0flag ) THEN
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           ELSE
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           ENDIF
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,2) - gamxinflu(i,j,6,2))
           tmp2 = gamxinflu(i,jp1,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,2) - gamxinflu(i,jp1,6,2))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1)) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))
           
          write(iunit,*) 'ibinhmlr = 1: i,k = ',igs(mgs),kgs(mgs)
          write(iunit,*) 'y,ratio,alpha,mltdiam1,xdia(mgs,lhl,1) = ',y,ratio,alpha(mgs,lhl),mltdiam1,xdia(mgs,lhl,1)
          write(iunit,*) 'full gamx = ',gamxinf(5.5 +0.5*bxx(mgs,lhl) + alpha(mgs,lhl), ratio),bxx(mgs,lhl),bx(lhl)
           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
          write(iunit,*) 'interpolated y = ', y,tmp1,tmp2,tmp3, &
               gamxinflu(i,j,6,2),gamxinflu(ip1,j,6,2),gamxinflu(i,jp1,6,2),gamxinflu(ip1,jp1,6,2)
           
        ENDIF

        hwvent1 =  0.78*x + y*hlventy(mgs)
!        hwvent1 =    &
!     &  ( 0.78*x +    &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       zhlmlrr(mgs) =   &
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 

        IF ( .false. .and. igs(mgs) == nx/2 ) THEN
         tmp3 = &
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*(0.78*x2 + y2*hlventy(mgs))*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 

          write(0,*) '3 zhlmlrr(mgs),zhlmlrr2 = ',zhlmlrr(mgs),tmp3
          write(0,*) '3 x,x2,y,y2 = ',x,x2,y,y2,hwvent1,0.78*x + y2*hlventy(mgs)
        
        ENDIF

! so now have zhmlrr as the Z whole-particle loss for hail for 0 to mltdiam1, and zhmlr0inf as the integrated whole-particle loss
! We want whole particle loss for 0 < d < mltdiam1, plus shedding loss for d > mltdiam1. Shedding loss is simple double the whole particle loss

       zhlmlrr(mgs) = Min( 0.0, zhlmlrr(mgs) )
       zhlmlr0inf   = Min( 0.0, zhlmlr0inf)
       zhlmlr(mgs) = zhlmlrr(mgs) + binmlrzrrfac*2.*(zhlmlr0inf - zhlmlrr(mgs)) ! factor of 2 to convert whole particle loss to shedding loss
       ! should that be zhlmlr(mgs) = 2.*zhlmlrr(mgs) + (zhlmlr0inf - zhlmlrr(mgs)) = zhlmlrr(mgs) + zhlmlr0inf
       ! yes, I think it should -- NOPE!
     !  zhlmlr(mgs) = 2.*zhlmlrr(mgs) + Min(0.0, (zhlmlr0inf - zhlmlrr(mgs)) )

        IF ( .false. .and. igs(mgs) == nx/2 ) THEN

          write(0,'(a,i3,10(1pe18.10))') 'kz, zhlmlr0inf,z0-zrr,zhlmlrr,zhlmlr  = ',kgs(mgs),zhlmlr0inf, &
                     zhlmlr0inf - zhlmlrr(mgs),zhlmlrr(mgs),zhlmlr(mgs)
          write(0,'(a,10(1pe18.10))') 'old zhlmlr = ',zhlmlrr(mgs) + 2.*(zhlmlr0inf - zhlmlrr(mgs))
          write(0,'(a,10(1pe18.10))' ) 'test zhlmlr = ',2.*zhlmlrr(mgs) + Min(0.0, (zhlmlr0inf - zhlmlrr(mgs)) )
          write(0,'(a,10(1pe18.10))') 'new - old = ',zhlmlr(mgs) - ( zhlmlrr(mgs) + 2.*(zhlmlr0inf - zhlmlrr(mgs)) )
        !  write(0,*) '3 x,x2,y,y2 = ',x,x2,y,y2,hwvent1,0.78*x + y2*hlventy(mgs)
        
        ENDIF

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr0 - tmp**2 * chlmlr(mgs)  ) )

         tmp1 = (xdn(mgs,lhl)/1000.)**2 * zhlmlrr(mgs) 
         IF ( numshedregimes == 2 ) THEN
          zhlmlrr(mgs) = (xdn(mgs,lhl)/1000.)**2 *zhlmlrr(mgs) + rho0(mgs)*  &
     &              (qhlmlr2*mltmass1inv*(0.01*takshedsize1)**6 + qhlmlr12*mltmass2inv*(0.01*takshedsize2)**6 )
         ELSEIF ( numshedregimes == 3 ) THEN
!           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
          zhlmlrr(mgs) = (xdn(mgs,lhl)/1000.)**2 *zhlmlrr(mgs) + rho0(mgs)*  &
     &              (qhlmlr3*mltmass1inv*(0.01*takshedsize1)**6 + qhlmlr23*mltmass2inv*(0.01*takshedsize2)**6  + &
     &               qhlmlr12*mltmass3inv*(0.01*takshedsize3)**6 )
!     &              (qhlmlr3*mltmass1inv*(0.01*0.5*takshedsize1)**6 + qhlmlr23*mltmass2inv*(0.01*0.5*takshedsize2)**6  + &
!     &               qhlmlr12*mltmass3inv*(0.01*0.5*takshedsize3)**6 )
         ENDIF
         
         IF ( .false. .and. igs(mgs) == nx/2 ) THEN
          write(0,'(a,i3,10(1pe18.10))') 'kz, zhlmlrr1,zhlmlrr,zhlmlr  = ',kgs(mgs),tmp1,zhlmlrr(mgs),zhlmlr(mgs)
          write(0,'(a,10(1pe18.10))') 'z1,z2,z3 = ', rho0(mgs)*qhlmlr3*mltmass1inv*(0.01*takshedsize1)**6 , &
     &               rho0(mgs)*qhlmlr23*mltmass2inv*(0.01*takshedsize2)**6 , &
     &               rho0(mgs)*qhlmlr12*mltmass3inv*(0.01*takshedsize3)**6 
!          write(0,'(a,10(1pe18.10))') 'z1,z2,z3 = ', rho0(mgs)*qhlmlr3*mltmass1inv*(0.01*0.5*takshedsize1)**6 , rho0(mgs)*qhlmlr23*mltmass2inv*(0.01*0.5*takshedsize2)**6 , &
!     &               rho0(mgs)*qhlmlr12*mltmass3inv*(0.01*0.5*takshedsize3)**6 
          write(0,*) 'chlmlrr = ',chlmlrr(mgs)
         ENDIF

         ! chlmlr(mgs)  = Min( chlmlr(mgs), (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)) ! Make sure hail size does not go *down*

!          tmp = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
         IF ( alpha(mgs,lhl) >= alphamax - delta_alphamlr ) THEN ! { cannot allow spectrum to get any narrower, revert to 2-moment method
         
          chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
!          chlmlrsave(mgs) = chlmlr(mgs)
!          qhlmlrsave(mgs) = qhlmlr(mgs)
!          chlsave(mgs) = cx(mgs,lhl)
!          qhlsave(mgs) = qx(mgs,lhl)
 
         

          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )

!        chmlrr(mgs) = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lh)*xv(mgs,lh)) ) ! into rain 
! guess what, this is the same as chmlr: rho0*qhmlr/xmas(lh) --> cx/qx = rho0/xmas
!          IF(imltshddmr == 1 ) THEN
!            ! DTD: If Dmg < sheddiam, then assume complete melting into
!            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
!            tmp = -rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl))) ! Min of Maximum raindrop size/mean hail size
!            tmp2 = -rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
!            chlmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)
!            chlmlrr(mgs) = -Max(tmp,Min(tmp2,chlmlrr(mgs)))
!          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
!            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
!            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
!            chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 
!          ENDIF


!          chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )

         
         ELSEIF (  xdia(mgs,lhl,3) > binmlrmxdia ) THEN ! 
         
          chlmlrtmp  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          zhlmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrtmp )
          
          tmp1 = 1.e-3
          IF ( xdia(mgs,lhl,3) > binmlrmxdia + tmp1 ) THEN 
             chlmlr(mgs) = chlmlrtmp
             zhlmlr(mgs) = zhlmlrtmp
          ELSE ! weighted average
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'k,xdia,chlmlr,zhlmlr = ',kgs(mgs),xdia(mgs,lhl,3)*1.e3,chlmlr(mgs),zhlmlr(mgs)
!              write(0,*) 'weights 1,2 = ',(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 , (binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            ENDIF
            chlmlr(mgs) = chlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + chlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
            zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
            !zhlmlr(mgs) = zhlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + zhlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'tmp,chlmlr,= ',chlmlrtmp,chlmlr(mgs)
!              write(0,*) 'tmp,zhlmlr,= ',zhlmlrtmp,zhlmlr(mgs)
!            ENDIF
          ENDIF

          ! check for negative value of cx, and, if so revert to typical number loss rate.
          IF ( cx(mgs,lhl) + dtp*chlmlr(mgs) < 0.0 ) THEN
             chlmlr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
             zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
          ENDIF
          
         
         ENDIF !} 


       ELSEIF ( ibinhlmlr == -1 ) THEN ! OLD VERSION use incomplete gamma functions to approximate the bin results

        ENDIF ! ibinhlmlr


       IF ( ivhmltsoak > 0 .and.  qhlmlr(mgs) < 0.0 .and. lvol(lhl) > 1 .and. xdn(mgs,lhl) .lt. xdnmx(lhl) ) THEN
         ! act as if 50% of the meltwater were soaked into the graupel
           v1 = (1. - xdn(mgs,lhl)/xdnmx(lhl))*(vx(mgs,lhl) + rho0(mgs)*qhlmlr(mgs)/xdn(mgs,lhl) )/(dtp) ! volume available for filling
           v2 = -1.0*rho0(mgs)*qhlmlr(mgs)/xdnmx(lhl)  ! volume of melted ice if it were refrozen in the matrix
           
           vhlsoak(mgs) = Min(v1,v2)
           
       ENDIF
        
        ENDIF
       ENDIF

      ENDIF
      
!
!      qimlr(mgs)  = max( qimlr(mgs), -qimxd(mgs) ) 
!      qsmlr(mgs)  = max( qsmlr(mgs),  -qsmxd(mgs) ) 
! erm 5/10/2007 changed to next line:
      if ( .not. mixedphase ) qsmlr(mgs)  = max( qsmlr(mgs),  Min( -qsmxd(mgs), -0.7*qx(mgs,ls)*dtpinv ) ) 
      IF ( .not. mixedphase ) THEN
        qhmlr(mgs)  = max( qhmlr(mgs),  Min( -qhmxd(mgs), -0.95*qx(mgs,lh)*dtpinv ) ) 
        chmlr(mgs)  = max( chmlr(mgs),  Min( -chmxd(mgs), -0.95*cx(mgs,lh)*dtpinv ) ) 
        IF ( lf > 0 ) THEN
          qfmlr(mgs)  = max( qfmlr(mgs),  Min( -qfmxd(mgs), -0.95*qx(mgs,lf)*dtpinv ) ) 
          cfmlr(mgs)  = max( cfmlr(mgs),  Min( -cfmxd(mgs), -0.95*cx(mgs,lf)*dtpinv ) ) 
        ENDIF
      ENDIF
!      qhmlr(mgs)  = max( max( qhmlr(mgs),  -qhmxd(mgs) ) , -0.5*qx(mgs,lh)*dtpinv ) !limits to 1/2 qh or max depletion
      qhmlh(mgs)  = 0. ! not used


      ! Rasmussen and Heymsfield say melt water remains on graupel up to 9 mm before shedding


      IF ( lhl .gt. 1 .and. lhlw < 1 ) THEN
        qhlmlr(mgs)  = max( qhlmlr(mgs),  Min( -qxmxd(mgs,lhl), -0.95*qx(mgs,lhl)*dtpinv ) )
        chlmlr(mgs)  = max( chlmlr(mgs),  Min( -cxmxd(mgs,lhl), -0.95*cx(mgs,lhl)*dtpinv ) )
      ENDIF

!
      end do

      endif  ! } not mixedphase
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cimlr(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qimlr(mgs)
      IF ( .not. mixedphase ) THEN !{
        IF ( xdia(mgs,ls,1) .gt. 1.e-6 .and. -qsmlr(mgs) .ge. 0.5*qxmin(ls) .and. ipconc .ge. 4 ) THEN 
!         csmlr(mgs)  = rho0(mgs)*qsmlr(mgs)/(xv(mgs,ls)*rhosm)
         csmlr(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)))*qsmlr(mgs)
        ELSEIF ( qx(mgs,ls) > qxmin(ls) ) THEN
         csmlr(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)))*qsmlr(mgs)
        ENDIF
        
        csmlrr(mgs) = csmlr(mgs)/rzxs(mgs)
         IF ( -csmlrr(mgs)*dtp > cxmin .and. -qsmlr(mgs)*dtp > qxmin(lr) .and. snowmeltdia > 0.0 ) THEN
           rmas = rho0(mgs)*qsmlr(mgs)/csmlrr(mgs)
           IF ( rmas > snowmeltmass ) THEN
             csmlrr(mgs) = rho0(mgs)*qsmlr(mgs)/snowmeltmass
           ENDIF
         ENDIF
           
         IF ( lf > 1 ) THEN !{
!         IF ( lzf < 1 ) THEN
           cfmlr(mgs)  = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfmlr(mgs)
!         ENDIF

      IF ( ihmlt .eq. 1 ) THEN
        cfmlrr(mgs)  = Min( cfmlr(mgs), rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      ELSEIF ( ihmlt .eq. 2 ) THEN
        IF ( xv(mgs,lf) .gt. 0.0 .and. cfmlr(mgs) .lt. 0.0 ) THEN
!        cfmlrr(mgs) = Min( cfmlr(mgs), rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lf)*xv(mgs,lf)) ) ! into rain 
! guess what, this is the same as cfmlr: rho0*qfmlr/xmas(lf) --> cx/qx = rho0/xmas
          IF(imltshddmr == 1) THEN
            ! DTD: If Dmg < sheddiam, then assume complete melting into
            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
            tmp = -rho0(mgs)*qfmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lf)*xv(mgs,lf))) ! Min of Maximum raindrop size/mean hail size
            tmp2 = -rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
            
            cfmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lf,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lf,3)-sheddiam)/(sheddiam0-sheddiam)  ! old version
            cfmlrr(mgs) = -Max(tmp,Min(tmp2,cfmlrr(mgs)))
          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lf))*xdia(mgs,lf,1)
            cfmlrr(mgs) =  rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lf))  ! into rain 
          ELSE ! Old method
            cfmlrr(mgs) =  rho0(mgs)*qfmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lf)*xv(mgs,lf)))  ! into rain 
         ENDIF
        ELSE
        cfmlrr(mgs) = cfmlr(mgs)
        ENDIF
      ELSEIF ( ihmlt .eq. 0 ) THEN
        cfmlrr(mgs) = cfmlr(mgs)
      ENDIF

         ENDIF !} lf > 1


!        IF ( xdia(mgs,lh,1) .gt. 1.e-6 .and. Abs(qhmlr(mgs)) .ge. qxmin(lh) ) THEN
!          chmlr(mgs) = rho0(mgs)*qhmlr(mgs)/(pi*xdn(mgs,lh)*xdia(mgs,lh,1)**3)  ! out of hail
!          chmlr(mgs) = Max( chmlr(mgs), -chmxd(mgs) )
!        ELSE
         IF ( ibinhmlr == 0 .or. lzh < 1 ) THEN
           chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
           IF ( imltshddmr == 3 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
            !  tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            !  
            !  IF ( tmpdiam > sheddiam ) THEN ! let size get smaller until it reaches sheddiam
            !   chmlr(mgs) = 0.0
            !  ENDIF
            
            ! test to remove the part of the melting associated with large ice particles so they get smaller

            tmp = 1. + alpha(mgs,lh)
            i = Int(dgami*(tmp))
            del = tmp - dgam*i
            g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

            ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lh,1) )

            x =  gamxinfdp(2. + alpha(mgs,lh), ratio)/g1palp
            y =  gamxinfdp(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp

            hwvent1 =  0.78*x + y*hwventy(mgs) 

            qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 )

            chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*(qhmlr(mgs) - qhlmlr1)
           
           
           ENDIF
!           IF ( igs(mgs) == 40 ) THEN
!             write(0,*) 'is this running? chmlr = ',kgs(mgs), chmlr(mgs)
!           ENDIF
         ENDIF
!        ENDIF


     IF ( chmlr(mgs) < 0.0 .and. (ibinhmlr < 1 .or. lzh < 1) ) THEN ! { already done if ibinhmlr > 0
      IF (  ipconc >= 6 .and. lzr .gt. 1 .and. lzh < 1  .and. qx(mgs,lh) > qxmin(lh) ) THEN ! Only compute if rain is 3-moment but graupel is not, otherwise is computed later
          tmp = qx(mgs,lh)/cx(mgs,lh)
          alp = alpha(mgs,lh)
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )

      ENDIF
      
      IF ( ibinhmlr == 0 .or. lzh < 1 ) THEN
      IF ( ihmlt .eq. 1 ) THEN
        chmlrr(mgs)  = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      ELSEIF ( ihmlt .eq. 2 ) THEN
        IF ( xv(mgs,lh) .gt. 0.0 .and. chmlr(mgs) .lt. 0.0 ) THEN
!        chmlrr(mgs) = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lh)*xv(mgs,lh)) ) ! into rain 
! guess what, this is the same as chmlr: rho0*qhmlr/xmas(lh) --> cx/qx = rho0/xmas
          IF(imltshddmr == 1) THEN
            ! DTD: If Dmg < sheddiam, then assume complete melting into
            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
            tmp = -rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh))) ! Min of Maximum raindrop size/mean hail size
            tmp2 = -rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
            
            chmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)  ! old version
            chmlrr(mgs) = -Max(tmp,Min(tmp2,chmlrr(mgs)))
          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
          ELSE ! Old method
            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh)))  ! into rain 
         ENDIF
        ELSE
        chmlrr(mgs) = chmlr(mgs)
        ENDIF
      ELSEIF ( ihmlt .eq. 0 ) THEN
        chmlrr(mgs) = chmlr(mgs)
      ENDIF

      ELSE ! ibinhmlr < 0? Already have an outer IF test for ibinhmlr < 1
        chmlrr(mgs)  = Min( chmlrr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*xvmx(lr)) ) ! into rain 
      ENDIF
      
      ENDIF ! } ( chmlr(mgs) < 0.0 .and. ibinhmlr < 1)

      IF ( lhl .gt. 1 .and. lhlw < 1 .and. .not. mixedphase .and. qhlmlr(mgs) < 0.0 ) THEN ! {
      
      IF ( ibinhlmlr == 0 .or. lzhl < 1 ) THEN
!      IF ( xdia(mgs,lhl,1) .gt. 1.e-6 .and. Abs(qhlmlr(mgs)) .ge. qxmin(lhl) ) THEN
!      chlmlr(mgs) = rho0(mgs)*qhlmlr(mgs)/(pi*xdn(mgs,lhl)*xdia(mgs,lhl,1)**3)  ! out of hail
!      chlmlr(mgs) = Max( chlmlr(mgs), -cxmxd(mgs,lhl) )
!      ELSE
      chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
           IF ( imltshddmr == 3 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!           IF ( .false. .and. imltshddmr == 3  ) THEN
!              tmpdiam = (shedalp+alpha(mgs,lhl))*xdia(mgs,lhl,1)
!              
!              IF ( tmpdiam > sheddiam ) THEN ! let size get smaller until it reaches sheddiam
!                chlmlr(mgs) = 0.0
!              ENDIF
 
            ! test to remove the part of the melting associated with large ice particles so they get smaller
!
            tmp = 1. + alpha(mgs,lhl)
            i = Int(dgami*(tmp))
            del = tmp - dgam*i
            g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

            ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )

            x =  gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
            y =  gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp

            hwvent1 =  0.78*x + y*hlventy(mgs) 

            qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 )

            chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*Min(0.0, qhlmlr(mgs) - qhlmlr1)

            ENDIF
!      ENDIF
      ENDIF
      
      IF ( ibinhlmlr == 0 .or. lzhl < 1 ) THEN !{
      IF ( ihmlt .eq. 1 ) THEN
        chlmlrr(mgs)  = Min( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      ELSEIF ( ihmlt .eq. 2 ) THEN
        IF ( xv(mgs,lhl) .gt. 0.0 .and. chlmlr(mgs) .lt. 0.0 ) THEN
!        chlmlrr(mgs) = rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 
!        chlmlrr(mgs) = Min( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lhl)*xv(mgs,lhl)) ) ! into rain 
          IF(imltshddmr == 1 ) THEN
            tmp = -rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl))) ! Min of Maximum raindrop size/mean hail size
            tmp2 = -rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
            chlmlrr(mgs) = tmp*(20.e-3-xdia(mgs,lhl,3))/(20.e-3-sheddiam)+tmp2*(xdia(mgs,lhl,3)-sheddiam)/(20.e-3-sheddiam)
            chlmlrr(mgs) = -Max(tmp,Min(tmp2,chlmlrr(mgs)))
          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 
          ELSE ! old method
            chlmlrr(mgs) = rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 
          ENDIF
        ELSE
        chlmlrr(mgs) = chlmlr(mgs)
        ENDIF
      ELSEIF ( ihmlt .eq. 0 ) THEN
        chlmlrr(mgs) = chlmlr(mgs)
      ENDIF

      ELSE ! } { ibinhlmlr > 0
        chlmlrr(mgs)  = Min( chlmlrr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*xvmx(lr)) ) ! into rain 
      ENDIF !}
      
        
       IF ( ipconc >= 8 .and. lzhl .gt. 1 .and. ibinhlmlr <= 0 ) THEN
        IF ( cx(mgs,lhl) > 0.0 ) THEN

          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          alp = alpha(mgs,lhl)
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( tmp * qhlmlr(mgs) )
       ENDIF
      ENDIF
      ENDIF ! }

      ENDIF ! }.not. mixedphase 

! 10ice versions:
!      chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
!      chmlrr(mgs) = chmlr(mgs)
      end do
      end if

!
!  deposition/sublimation of ice
!
      DO mgs = 1,ngscnt

      rwcap(mgs) = (0.5)*xdia(mgs,lr,1)
      swcap(mgs) = (0.5)*xdia(mgs,ls,1)
      hwcap(mgs) = (0.5)*xdia(mgs,lh,1)
      IF ( lf > 1 ) fwcap(mgs) = (0.5)*xdia(mgs,lf,1)
      IF ( lhl .gt. 1 ) hlcap(mgs) = (0.5)*xdia(mgs,lhl,1)

      if ( qx(mgs,li).gt.qxmin(li) .and. xdia(mgs,li,1) .gt. 0.0 ) then
!
! from Cotton, 1972 (Part II)
!
        cilen(mgs)   = 0.4764*(xdia(mgs,li,1))**(0.958)
        cval = xdia(mgs,li,1)
        aval = cilen(mgs)
        eval = Sqrt(1.0-(aval**2)/(cval**2))
        fval = min(0.99,eval)
        gval = alog( abs( (1.+fval)/(1.-fval) ) )
        cicap(mgs) = cval*fval / gval
      ELSE
       cicap(mgs) = 0.0
      end if
      ENDDO
!
!
      qhdsv(:) = 0.0
      qhldsv(:) = 0.0
      qfdsv(:) = 0.0

      do mgs = 1,ngscnt
      IF ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) THEN
        qidsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,li)*civent(mgs)*cicap(mgs)*depfac
        qsdsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,ls)*swvent(mgs)*swcap(mgs)*depfac

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidsv = ',nstep,kgs(mgs),qidsv(mgs),temg(mgs)-tfrh,100.*(qx(mgs,lv)/qis(mgs) - 1.),1.e6*xdia(mgs,li,1),
!     :            fvds(mgs),civent(mgs),cicap(mgs)
!        ENDIF
      ELSE
        qidsv(mgs) = 0.0
        qsdsv(mgs) = 0.0
      ENDIF
        qhdsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,lh)*hwvent(mgs)*hwcap(mgs)*depfac

        IF ( lf > 1 ) THEN
        qfdsv(mgs) = fvds(mgs)*cx(mgs,lf)*fwvent(mgs)*fwcap(mgs)*depfac
        ENDIF
        IF ( lhl .gt. 1 ) qhldsv(mgs) = fvds(mgs)*cx(mgs,lhl)*hlvent(mgs)*hlcap(mgs)*depfac
!
!
      end do
!


! #include "nssl.qlimit.F"

!
!  Use a test saturation adjustment to set limits on ice deposition/sublimation
!  and rain evaporation
!
!
      IF ( DoSublimationFix ) THEN
      
      do mgs = 1,ngscnt

        qitmp(mgs) = qx(mgs,li) + qx(mgs,ls) + qx(mgs,lh)
        IF ( lis > 1 ) qitmp(mgs) = qitmp(mgs) + qx(mgs,lis)
        IF ( lhl > 1 ) qitmp(mgs) = qitmp(mgs) + qx(mgs,lhl)
        qrtmp(mgs) = qx(mgs,lr)
        qctmp(mgs) = qx(mgs,lc)
        qsimxdep(mgs) = 0.0
        qsimxsub(mgs) = 0.0
        dqcitmp(mgs) = 0.0
        

!      IF ( ( qitmp(mgs) > qxmin(li) .or. qrtmp(mgs) > qxmin(lr) ) ) THEN
      IF ( qitmp(mgs) > qxmin(li)  ) THEN
      
        qitmp1    = qitmp(mgs)
        qctmp1    = qctmp(mgs)
        felvcptmp = felvcp(mgs)
        felscptmp = felscp(mgs)
        qvtmp(mgs) = qx(mgs,lv)
        qss(mgs) = qvs(mgs)
        qsstmp = qvs(mgs)
        qvstmp = qvs(mgs)
        qisstmp = qis(mgs)
        thetatmp  = theta(mgs)
        thetaptmp = thetap(mgs)
        temgtmp   = temg(mgs)
        temcgtmp  = temcg(mgs)
        qvaptmp   = qx(mgs,lv) ! qwvp(mgs) + qv0(mgs)
        qvptmp    = 0.0 ! qwvp(mgs)  ! qv pertubation

        qsstmp = qisstmp

      
       dqwvtmp(mgs) = ( qvtmp(mgs) - qsstmp )

      do itertd = 1,2
      
!
!  calculate super-saturation
!
      IF ( itertd == 1 ) THEN
      
      ELSE
        dqcitmp(mgs) = dqci(mgs)
   !     dqwvtmp(mgs) = dqwv(mgs)
      ENDIF

      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qvtmp(mgs) - qsstmp )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           ! { subsaturated
        if( qitmp(mgs) .gt. -dqwv(mgs) ) then  ! check if qi can make up all the deficit
          dqci(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                  ! otherwise make all ice available for sublimation
          dqci(mgs) = -qitmp(mgs)
          dqwv(mgs) = dqwv(mgs) + qitmp(mgs)
        end if
!
       qvptmp = qvptmp - ( dqcw(mgs) + dqci(mgs) )  ! add to perturbation vapor

       IF ( itertd == 2 .and. eqtset > 1 ) THEN
       ! if eqtset == 2, then need to update the latent heats for change in hydrometeor content
          tmp = qitmp(mgs) !+ qx(mgs,lh)
!          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          cvm = cv+cvv*qvtmp(mgs)+cpl*(qx(mgs,lc)+qrtmp(mgs))   &
                                  +cpigb*(tmp)

          felvcptmp = (felv(mgs)-rw*temg(mgs))/cvm
          felscptmp = (fels(mgs)-rw*temg(mgs))/cvm
       ENDIF


!      qitmp(mgs) = qx(mgs,li)
      qctmp(mgs) = qctmp(mgs) + dqcw(mgs) ! dqcw is zero
      qitmp(mgs) = qitmp(mgs) + dqci(mgs)
      thetaptmp = thetaptmp +   &
     &  1./pi0(mgs)*   &
     &  (felvcp(mgs)*dqcw(mgs) +felscp(mgs)*dqci(mgs))


      end if  ! } dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      IF ( dqwv(mgs) .ge. 0. ) THEN ! {
      
!      write(iunit,*) 'satadj: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
!
!        qitmp(mgs) = qx(mgs,li)
        fracl(mgs) = 0.0
        fraci(mgs) = 1.0
        if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
!          fracl(mgs) = max(min(1.,(temg(mgs)-233.15)/(20.)),0.0)
!          fraci(mgs) = 1.0-fracl(mgs)
        end if
        if ( temg(mgs) .le. thnuc ) then
           fraci(mgs) = 1.0
           fracl(mgs) = 0.0
         end if
!        fraci(mgs) = 1.0-fracl(mgs)

       gamss = (felvcp(mgs)*fracl(mgs) + felscp(mgs)*fraci(mgs))   &
     &      / (pi0(mgs))

          dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv2(mgs)*qsstmp/   &
     &  ((temg(mgs)-cbi)**2))

      if ( temg(mgs) .ge. tfr ) then
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qsstmp/   &
     &  ((temg(mgs)-cbw)**2))
      end if

      delqci1=qx(mgs,li)


      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs) ! is zero
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)

      thetaptmp = thetaptmp +   &
     &   (felvcp(mgs)*dqcw(mgs) + felscp(mgs)*dqci(mgs))   &
     & / (pi0(mgs))

      qvptmp = qvptmp - ( dqvcnd(mgs) )
      qctmp(mgs) = qctmp(mgs) + dqcw(mgs)
      qitmp(mgs) = qitmp(mgs) + dqci(mgs)

       IF ( itertd == 2 .and. eqtset > 1 ) THEN
       ! if eqtset == 2, then need to update the latent heats for change in hydrometeor content
          tmp = qitmp(mgs) ! + qx(mgs,lh)
!          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          cvm = cv+cvv*qvtmp(mgs)+cpl*(qctmp(mgs) +qrtmp(mgs))   &
                                  +cpigb*(tmp)

          felvcptmp = (felv(mgs)-rw*temg(mgs))/cvm
          felscptmp = (fels(mgs)-rw*temg(mgs))/cvm
       ENDIF

!       IF ( eqtset > 2 ) THEN
!         pipert(mgs) = pipert(mgs) + (0   &
!      &  +felspi(mgs)*dqci(mgs)    &
!      &  +felvpi(mgs)*dqcw(mgs)) ! *dtp
!       ENDIF

!
!
      END IF ! } dqwv(mgs) .ge. 0.


!
      IF ( itertd == 1 ) THEN
      ! update temporary saturation values

      thetatmp = thetaptmp + theta0(mgs)
      temgtmp = thetatmp*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      qvaptmp = Max((qvptmp + qv0(mgs)), 0.0)
      temcgtmp = temgtmp - tfr
      tqvcon = temgtmp-cbw
      ltemq = (temgtmp-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      IF ( iqvsopt == 0 ) THEN
        qvstmp = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvstmp = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temgtmp,0))
        qvstmp = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF

      qisstmp = pqs(mgs)*tabqis(ltemq)
      qctmp(mgs) = max( 0.0, qctmp(mgs) )
      qitmp(mgs) = max( 0.0, qitmp(mgs) )
      qvtmp(mgs) = max( 0.0, qvaptmp )
      
!      qsstmp = qvstmp
      qsstmp = qisstmp
      
      ELSE
       ! set max depletion
        qctmp(mgs) = max( 0.0, qctmp(mgs) )
        qitmp(mgs) = max( 0.0, qitmp(mgs) )
       
        IF ( qitmp(mgs) < qitmp1 ) THEN
          qsimxsub(mgs) = (qitmp1 - qitmp(mgs))*dtpinv
        ELSEIF ( qitmp(mgs) > qitmp1 ) THEN
          qsimxdep(mgs) = (qitmp(mgs) - qitmp1)*dtpinv
        ENDIF
       
      
      ENDIF
!      pceds(mgs) = (thetap(mgs) - thsave(mgs))*dtpinv
!      write(iunit,*) 'satadj2: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qxtmp,qctmp(mgs)
!
!  end the saturation adjustment iteration loop
!
      end do ! itertd
      
      ENDIF
      
      end do ! mgs
      
      ELSE
      
       DO mgs = 1,ngscnt
         qsimxdep(mgs) = qvimxd(mgs)
         qsimxsub(mgs) = 1.e20
       ENDDO
      
      ENDIF

! end of qlimit

      qhcev(:) = 0.0
      chcev(:) = 0.0
      qhlcev(:) = 0.0
      chlcev(:) = 0.0
      qfcev(:) = 0.0
      cfcev(:) = 0.0

      do mgs = 1,ngscnt
      qisbv(mgs) = 0.0
      qssbv(mgs) = 0.0
      qidpv(mgs) = 0.0
      qsdpv(mgs) = 0.0
      qhsbv(mgs) = 0.0
      qscev(mgs) = 0.0
      cscev(mgs) = 0.0
      qfsbv(mgs)  = 0.0
      IF ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) THEN ! last condition (qr<qmin & qc<qmin) for case icond=0
!        qisbv(mgs) = max( min(qidsv(mgs), 0.0), -qimxd(mgs) )
!        qssbv(mgs) = max( min(qsdsv(mgs), 0.0), -qsmxd(mgs) )
! erm 5/10/2007:
        qisbv(mgs) = max( min(qidsv(mgs), 0.0),  Min( -qimxd(mgs), -0.5*qx(mgs,li)*dtpinv ) )
        IF ( temg(mgs) < tfr .or. .not. qsmlr(mgs) < 0.0 ) THEN
        qssbv(mgs) = max( min(qsdsv(mgs), 0.0),  Min( -qsmxd(mgs), -0.5*qx(mgs,ls)*dtpinv ) )
        ENDIF
        qidpv(mgs) = Max(qidsv(mgs), 0.0)
        qsdpv(mgs) = Max(qsdsv(mgs), 0.0)
        
        IF ( qsmlr(mgs) < 0.0 .and. .not. mixedphase ) THEN ! switch snow sublimation to evaporation if there is melting

          qscev(mgs) = evapfac*   &
     &  4.0*pi*(qx(mgs,lv)-qss0(mgs))*cx(mgs,ls)*swcap(mgs)*swvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))
          qscev(mgs) = Max( Min(0.0,qscev(mgs)),  Min( -qsmxd(mgs), -0.5*qx(mgs,ls)*dtpinv ) )
        ELSE

        ENDIF



      ELSE
        qisbv(mgs) = 0.0
        qssbv(mgs) = 0.0
        qidpv(mgs) = 0.0
        qsdpv(mgs) = 0.0
      ENDIF

      qhsbv(mgs) = 0.0
      qhdpv(mgs) = 0.0
      IF ( qx(mgs,lh) > qxmin(lh) ) THEN
      IF ( temg(mgs) < tfr .or. .not. qhmlr(mgs) < 0.0 ) THEN
      ! no liquid from melting, so evaporation is greater. Thus can calculate sublimation rate
      qhsbv(mgs) = max( min(qhdsv(mgs), 0.0), -qhmxd(mgs) )
      qhdpv(mgs) = Max(qhdsv(mgs), 0.0)
      ENDIF
      
      IF ( .true. .and. qhmlr(mgs) < 0.0 .and. .not. mixedphase ) THEN
        ! Liquid is forming, so find the evaporation that was subtracted from melting (if it is not condensing)
!       qhcev(mgs) =   &
!     &   evapfac*min(   &
!     &  fmlt1e(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1), 0.0 )
        
        qhcev(mgs) =  evapfac*2.0*pi*(qx(mgs,lv)-qss0(mgs))*  &
     &   cx(mgs,lh)*xdia(mgs,lh,1)*hwvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))

        qhcev(mgs)  = max(qhcev(mgs), -qhmxd(mgs))
        IF ( temg(mgs) > tfr ) qhcev(mgs) = Min(0.0, qhcev(mgs) )
        
      ENDIF
      ENDIF

      IF ( lf > 1 ) THEN
      qfdpv(mgs) = 0.0
      IF ( qx(mgs,lf) > qxmin(lf) ) THEN
        qfsbv(mgs) = max( min(qfdsv(mgs), 0.0), -qxmxd(mgs,lf) )
        qfdpv(mgs) = Max(qfdsv(mgs), 0.0)

        IF ( qfmlr(mgs) < 0.0 .and. .not. mixedphase ) THEN
        ! Liquid is forming, so find the evaporation that was subtracted from melting (if it is not condensing)
         qfcev(mgs) =  evapfac*2.0*pi*(qx(mgs,lv)-qss0(mgs))*   &
     &      cx(mgs,lf)*xdia(mgs,lf,1)*fwvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))

         qfcev(mgs)  = max(qfcev(mgs), -qfmxd(mgs))
         IF ( temg(mgs) > tfr ) qfcev(mgs) = Min(0.0, qfcev(mgs) )

        ENDIF

      ENDIF
      ENDIF

      qhlsbv(mgs) = 0.0
      qhldpv(mgs) = 0.0
      IF ( lhl .gt. 1 ) THEN
      IF ( qx(mgs,lhl) > qxmin(lhl) ) THEN
        IF ( temg(mgs) < tfr .or. .not. qhlmlr(mgs) < 0.0 ) THEN
        qhlsbv(mgs) = max( min(qhldsv(mgs), 0.0), -qxmxd(mgs,lhl) )
        qhldpv(mgs) = Max(qhldsv(mgs), 0.0)
        ENDIF
        IF ( qhlmlr(mgs) < 0.0 .and. .not. mixedphase ) THEN
        ! Liquid is forming, so find the evaporation that was subtracted from melting (if it is not condensing)
         qhlcev(mgs) =  evapfac*2.0*pi*(qx(mgs,lv)-qss0(mgs))*  &
     &      cx(mgs,lhl)*xdia(mgs,lhl,1)*hlvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))

         qhlcev(mgs)  = max(qhlcev(mgs), -qhlmxd(mgs))
         IF ( temg(mgs) > tfr ) qhlcev(mgs) = Min(0.0, qhlcev(mgs) )
        
      ENDIF
      ENDIF
      ENDIF
      
      temp1 = qidpv(mgs) + qsdpv(mgs) + qhdpv(mgs) + qhldpv(mgs)

!      IF ( temp1 .gt. qvimxd(mgs) ) THEN

!      frac = qvimxd(mgs)/temp1

      IF ( temp1 .gt. qsimxdep(mgs) ) THEN
      frac = qsimxdep(mgs)/temp1

      qidpv(mgs) = frac*qidpv(mgs)
      qsdpv(mgs) = frac*qsdpv(mgs)
      qhdpv(mgs) = frac*qhdpv(mgs)
      qhldpv(mgs) = frac*qhldpv(mgs)

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidpv,frac = ',kgs(mgs),qidpv(mgs),frac
!        ENDIF

      ENDIF

      temp1 = qisbv(mgs) + qssbv(mgs) + qhsbv(mgs) + qhlsbv(mgs)
        temp1 = temp1 + qfsbv(mgs)


      IF ( temp1 <  -qsimxsub(mgs) ) THEN
      frac = -qsimxsub(mgs)/temp1

      qisbv(mgs) = frac*qisbv(mgs)
        qfsbv(mgs) = frac*qfsbv(mgs)
      qssbv(mgs) = frac*qssbv(mgs)
      qhsbv(mgs) = frac*qhsbv(mgs)
      qhlsbv(mgs) = frac*qhlsbv(mgs)

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidpv,frac = ',kgs(mgs),qidpv(mgs),frac
!        ENDIF

      ENDIF


      end do
!
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cssbv(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qssbv(mgs)
      cisbv(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qisbv(mgs)
       IF ( lf > 1 )  cfsbv(mgs) = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfsbv(mgs)
      chsbv(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhsbv(mgs)
      IF ( lhl .gt. 1 ) chlsbv(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlsbv(mgs)
      csdpv(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qsdpv(mgs)
      cidpv(mgs) =  0.0 ! (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qidpv(mgs)
      cisdpv(mgs) = 0.0
      chdpv(mgs)  = 0.0 ! (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhdpv(mgs)
      chldpv(mgs) = 0.0
      end do
      end if

!
!  Aggregation or size conversion of small crystals to snow
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 29a'
      do mgs = 1,ngscnt
      qscni(mgs) =  0.0
      cscni(mgs) = 0.0
      cscnis(mgs) = 0.0
      if ( ipconc .ge. 4 .and. iscni .ge. 1 .and. qx(mgs,li) .gt. qxmin(li) ) then
        IF ( iscni .eq. 1 ) THEN
         qscni(mgs) =    &
     &      pi*rho0(mgs)*((0.25)/(6.0))   &
     &      *eii(mgs)*(qx(mgs,li)**2)*(xdia(mgs,li,2))   &
     &      *vtxbar(mgs,li,1)/xmas(mgs,li)
         cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li))
         cscnis(mgs) = 0.5*cscni(mgs)
        ELSEIF ( iscni .eq. 2 .or. iscni .eq. 4 .or. iscni .eq. 5 ) THEN  ! Zeigler 1985/Zrnic 1993, sort of
          IF ( iscni .ne. 5 .and. qidpv(mgs) .gt. 0.0 .and.  xdia(mgs,li,3) .ge. 100.e-6 ) THEN
          ! convert larger crystals to snow
!            IF ( xdia(mgs,ls,3) .gt. xdia(mgs,li,3) ) THEN
!              qscni(mgs) = Max(0.1,xdia(mgs,li,3)/xdia(mgs,ls,3))*qidpv(mgs)
! erm 9/5/08 changed max to min
              qscni(mgs) = Min(0.5, xdia(mgs,li,3)/200.e-6)*qidpv(mgs)
!            ELSE
!              qscni(mgs) = 0.1*qidpv(mgs)
!            ENDIF
            cscni(mgs) = fscni*qscni(mgs)*rho0(mgs)/Max(rho_qs*xvmn(ls),xmas(mgs,li))
!            cscni(mgs) = fscni*Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/Max(xdn(mgs,ls)*xvmn(ls),xmas(mgs,li)))
!            cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li) )
!            IF ( xdia(mgs,ls,3) .le. 200.e-6 ) THEN
              cscnis(mgs) = cscni(mgs)
!            ELSE
!              cscnis(mgs) = 0.0
!            ENDIF
            !  write(91,*) 'qi,qscni = ',igs(mgs),kgs(mgs),qx(mgs,li),qscni(mgs),cscnis(mgs),qidpv(mgs)
          ENDIF
           IF ( iscni .ne. 4 ) THEN
           ! crystal aggregation to become snow
! erm 9/5/08 commented second line and added xv to 1st line (zrnic et al 1993)
             tmp = ess(mgs)*rvt*aa2*cx(mgs,li)*cx(mgs,li)*xv(mgs,li)
!     :         ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,li))

!           csacs(mgs) = rvt*aa2*ess(mgs)*cx(mgs,ls)**2*xv(mgs,ls)

             qscni(mgs) = qscni(mgs) + Min( qxmxd(mgs,li), 2.0*tmp*xmas(mgs,li)*rhoinv(mgs) )
             cscni(mgs) = cscni(mgs) + Min( cxmxd(mgs,li), 2.0*tmp )
             cscnis(mgs) = cscnis(mgs) + Min( cxmxd(mgs,li), tmp )
           ENDIF
        ELSEIF ( iscni .eq. 3 ) THEN ! LFO
           qscni(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
           qscni(mgs) = min(qscni(mgs),qxmxd(mgs,li))
           cscni(mgs) = qscni(mgs)*rho0(mgs)/xmas(mgs,li)
           cscnis(mgs) = 0.5*cscni(mgs)
!           write(iunit,*) 'qscni, qi = ',qscni(mgs),qx(mgs,li),igs(mgs),kgs(mgs)
        ENDIF

      ELSEIF ( ipconc < 4 ) THEN ! LFO
           IF ( lwsm6 ) THEN
             qimax = rhoinv(mgs)*roqimax
             qscni(mgs) = Min(0.90*qx(mgs,li), Max( 0.0, (qx(mgs,li) - qimax)*dtpinv ) )
           ELSE
             qscni(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
             qscni(mgs) = min(qscni(mgs),qxmxd(mgs,li))
           ENDIF
      else ! 10-ice version
      if ( iscni > 0 .and. qx(mgs,li) .gt. qxmin(li) ) then
          qscni(mgs) =    &
     &    pi*rho0(mgs)*((0.25)/(6.0))   &
     &    *eii(mgs)*(qx(mgs,li)**2)*(xdia(mgs,li,2))   &
     &    *vtxbar(mgs,li,1)/xmas(mgs,li)
         cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li))
        end if

      end if
      end do

      IF ( incwet < 1 ) THEN
          dhwet(:) = d1t
          dhlwet(:) = d1t
          dfwet(:) = d1t
      ENDIF

         IF ( incwet >= 1 ) THEN
         ! 'incwet' = incomplete gamma for wet growth
         ! Find diameter where wet growth starts, then compute dry and wet growth 
         ! over [dwet,infinity]. Subtract dry growth from qxacw etc. to get total
         ! dry growth part
         dhwet(:) = dg0thresh + 0.0001
         dhlwet(:) = dg0thresh + 0.0001
         dfwet(:) = dg0thresh + 0.0001
         
         DO mgs = 1,ngscnt

             sqrtrhovt = Sqrt( rhovt(mgs) )
             fventh = sqrtrhovt*(fpndl(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5) 
             fventm = sqrtrhovt*(fschm(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5)
             ltemq = (tfr-163.15)/fqsat+1.5
             qvs0 = pqs(mgs)*tabqvs(ltemq)
             denomdp = felf(mgs) + fcw(mgs)*temcg(mgs)
             denominvdp = 1.d0/(felf(mgs) + fcw(mgs)*temcg(mgs))
         
         IF (((qhacw(mgs) + qhacr(mgs))*dtp > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin .and. &
              temg(mgs) .le. tfr + wetgrthtoffset .and.  temg(mgs) .ge. 243.15 )  ) THEN
!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*ehw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
            x =   1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            IF ( x > 1.e-20 ) THEN
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            ELSE
              dwr = 1.e30
            ENDIF
          d = dwr

           IF ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) THEN

                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = ehi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, ehw(mgs))*qx(mgs,lc) 
                      h4 = ehr(mgs)* qx(mgs,lr)
                      ! iterate to find minimum diameter for wet growth. Start with value of dwr
                      DO n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lh)*d**bxx(mgs,lh) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       IF ( x2 > 1.4 ) THEN
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       ELSE
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       ENDIF


                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)

                        IF ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) EXIT
                        
                      ENDDO
              ENDIF
              
              dhwet(mgs) = Min(dg0thresh + 0.0001, Max( d, dwetmin ))
          ELSE
            dhwet(mgs) = dg0thresh + 0.0001
          ENDIF

         IF (((qhlacw(mgs) + qhlacr(mgs))*dtp > qxmin(lhl) .and. qx(mgs,lhl) > 0.01e-3 &
               .and. temg(mgs) .le. tfr + wetgrthtoffset  .and.  temg(mgs) .ge. 243.15 )  ) THEN
!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*ehlw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(ehlw(mgs)*qx(mgs,lc)+ehlr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
            x =   1.1e4 * rho0(mgs)*(ehlw(mgs)*qx(mgs,lc)+ehlr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            IF ( x > 1.e-20 ) THEN
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            ELSE
              dwr = 1.e30
            ENDIF
          d = dwr
           IF ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) THEN

!                      write(91,*) 'dw,dwr,temcg = ',100.*dw,100.*dwr,temcg(mgs)
                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = ehi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, ehlw(mgs))*qx(mgs,lc) 
                      h4 = ehlr(mgs)* qx(mgs,lr)
                      ! iterate to find minimum diameter for wet growth. Start with value of dwr
                      DO n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lhl)*d**bxx(mgs,lhl) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       IF ( x2 > 1.4 ) THEN
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       ELSE
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       ENDIF


                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)

                        IF ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) EXIT
                        
                      ENDDO
              ENDIF
              
              dhlwet(mgs) = Min(dg0thresh + 0.0001, Max( d, dwetmin ) )
          ELSE
            dhlwet(mgs) = dg0thresh + 0.0001
          ENDIF

         IF ( lf > 0 ) THEN
         IF (((qfacw(mgs) + qfacr(mgs))*dtp > qxmin(lf) .and. qx(mgs,lf) > 0.01e-3 .and. &
               temg(mgs) .le. tfr + wetgrthtoffset .and.  temg(mgs) .ge. 243.15 )  ) THEN
!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*efw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(efw(mgs)*qx(mgs,lc)+efr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
            x =   1.1e4 * rho0(mgs)*(efw(mgs)*qx(mgs,lc)+efr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            IF ( x > 1.e-20 ) THEN
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            ELSE
              dwr = 1.e30
            ENDIF
          d = dwr

           IF ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) THEN

!                      write(91,*) 'dw,dwr,temcg = ',100.*dw,100.*dwr,temcg(mgs)
                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = ehi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, efw(mgs))*qx(mgs,lc) 
                      h4 = efr(mgs)* qx(mgs,lr)
                      ! iterate to find minimum diameter for wet growth. Start with value of dwr
                      DO n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lf)*d**bxx(mgs,lf) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       IF ( x2 > 1.4 ) THEN
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       ELSE
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       ENDIF


                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)

                        IF ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) EXIT
                        
                      ENDDO
              ENDIF
              
              dfwet(mgs) = Min(dg0thresh + 0.0001, Max( d, dwetmin ))
          ELSE
            dfwet(mgs) = dg0thresh + 0.0001
          ENDIF
          
          ENDIF 

            
          ENDDO
          
          ENDIF ! incwet



!
!
!  compute dry growth rate of snow, graupel, and hail
!
      do mgs = 1,ngscnt
!
      qsdry(mgs)  = qsacr(mgs)    + qsacw(mgs)   &
     &            + qsaci(mgs)
!
      qhdry(mgs)  = qhaci(mgs)    + qhacs(mgs)   &
     &            + qhacr(mgs)   &
     &            + qhacw(mgs)
!

      IF ( lf > 1 ) THEN
       qfdry(mgs)  = qfaci(mgs) + qfacs(mgs)   &
     &               + qfacr(mgs)   &
     &               + qfacw(mgs)
      ! qfwet(mgs) = qfdry(mgs)
      ENDIF
      qhldry(mgs) = 0.0
      IF ( lhl .gt. 1 ) THEN
      qhldry(mgs)  = qhlaci(mgs)    + qhlacs(mgs)   &
     &               + qhlacr(mgs)   &
     &               + qhlacw(mgs)
      ENDIF
      end do
!
!  set wet growth and shedding
!
      do mgs = 1,ngscnt
      
      IF ( tfrdry < temg(mgs) .and. temg(mgs) < tfr ) THEN ! {
!
!      qswet(mgs) =
!     >  ( xdia(mgs,ls,1)*swvent(mgs)*cx(mgs,ls)*fwet1(mgs)
!     >  + fwet2(mgs)*(qsaci(mgs)+qsacir(mgs)
!     >               +qsacip(mgs)) )
!      qswet(mgs) = max( 0.0, qswet(mgs))
!
!      IF ( dnu(lh) .ne. 0. ) THEN
!        qhwet(mgs) = qhdry(mgs)
!      ELSE
      ! IF ( incwet == 0 ) THEN
        qhwet(mgs) =   &
     &    ( xdia(mgs,lh,1)*hwvent(mgs)*cx(mgs,lh)*fwet1(mgs)   &
     &   + fwet2(mgs)*(qhaci(mgs) + qhacs(mgs)) )
        qhwet(mgs) = max( 0.0, qhwet(mgs))

       IF ( incwet == 1 .and. qhwet(mgs) < qhdry(mgs) .and. dhwet(mgs) < dg0thresh ) THEN
       !  ELSE
        !   IF ( dhwet(mgs) < dg0thresh ) THEN
             ! find portion of qc and qr collection that are dry/wet growth for d > dwet

               ratio = Min( maxratiolu, dhwet(mgs)/xdia(mgs,lh,1) )
               
               tmp1 = gaminterp(ratio,alpha(mgs,lh),13,1) ! alpha + 3
               tmp2 = gaminterp(ratio,alpha(mgs,lh),12,1) ! alpha + 2
               tmp3 = gaminterp(ratio,alpha(mgs,lh), 9,1) ! alpha + 1

            IF ( qhacw(mgs)*dtp > qxmin(lh) ) THEN
              vt = abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))

          ! dry growth of qc for D > Dwet to substract from qhacw
          qxacwtmp = 0.25*pi*ehw(mgs)*cx(mgs,lh)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab1lh(mgs,lh,lc)*xdia(mgs,lh,3)*xdia(mgs,lc,3) +    &
     &            tmp3*da1lc(mgs)*xdia(mgs,lc,3)**2 )
             ELSE
              qxacwtmp = 0.0
             ENDIF

            IF ( qhacr(mgs)*dtp > qxmin(lh) ) THEN

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lr,1) )

          ! dry growth of qr for D > Dwet to substract from qhacr
         qxacrtmp = 0.25*pi*ehr(mgs)*cx(mgs,lh)*qx(mgs,lr)*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab1lh(mgs,lh,lr)*xdia(mgs,lh,3)*xdia(mgs,lr,3) +    &
     &            tmp3*da1lr(mgs)*xdia(mgs,lr,3)**2 )
             ELSE
               qxacrtmp = 0.0
             ENDIF

        ! hwvent is where the size dependency is, so hxventtmp gives the portion for d > dwet
        x = gaminterp(ratio,alpha(mgs,lh),9,1) ! alpha + 1
        y = gaminterp(ratio,alpha(mgs,lh),3,1) ! alpha + b/2 + 5/2
        
        hxventtmp =  0.78*x + y*hwventy(mgs)  !   &

        ! find the ice and snow collection for d > dwet
        qxacitmp = 0.0
        IF ( qhaci(mgs)*dtp > qxmin(lh) ) THEN
              vt = abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))

          ! note that ehi=1 implicitly here
          qxacitmp = 0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*qx(mgs,li)*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab1lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            tmp3*da1(li)*xdia(mgs,li,3)**2 )
        
          cxacitmp =    &
     &        0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*cx(mgs,li)*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab0lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            tmp3*da0(li)*xdia(mgs,li,3)**2 )
        ENDIF

        qxacstmp = 0.0
        IF ( qhacs(mgs)*dtp > qxmin(lh) ) THEN
              vt = abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))

          ! note that ehs=1 implicitly here
          qxacstmp = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*qx(mgs,ls)*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab1lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da1(ls)*xdia(mgs,ls,3)**2 )

          cxacstmp = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*cx(mgs,ls)*vt*   &
     &         (  tmp1*da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            tmp2*dab0lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da0(ls)*xdia(mgs,ls,3)**2 )
        ENDIF

             qxwettmp =   &
     &       xdia(mgs,lh,1)*hxventtmp*cx(mgs,lh)*fwet1(mgs)   &
     &     + fwet2(mgs)*(qxacitmp + qxacstmp)

          tmp = qhwet(mgs)
          ! as dry growth but subtract part for D > Dw and add wet growth for D > Dw
          qhwet(mgs) = qhacw(mgs) + qhacr(mgs) + qhaci(mgs) + qhacs(mgs) &
                        - ehi(mgs)*qxacitmp - ehs(mgs)*qxacstmp          &
                        -  qxacwtmp - qxacrtmp + qxwettmp

          qhaci(mgs) =  qhaci(mgs) + (1.0 - ehi(mgs))*qxacitmp
          qhacs(mgs) =  qhacs(mgs) + (1.0 - ehs(mgs))*qxacstmp
          chaci(mgs) =  chaci(mgs) + (1.0 - ehi(mgs))*cxacitmp
          chacs(mgs) =  chacs(mgs) + (1.0 - ehs(mgs))*cxacstmp

         ! qhacw(mgs) = Min( qhacw(mgs), 0.5*qx(mgs,lc)*dtpinv )
           
       !    ELSE ! for dwet > 15cm, just assume dry growth
       !      qhwet(mgs) = qhdry(mgs)
       !    ENDIF
         ENDIF

!      ENDIF

      IF ( lf > 1 ) THEN
       !  IF ( incwet == 0 ) THEN
        qfwet(mgs) =   &
     &    ( xdia(mgs,lf,1)*fwvent(mgs)*cx(mgs,lf)*fwet1(mgs)   &
     &   + fwet2(mgs)*(qfaci(mgs) + qfacs(mgs)) )
        qfwet(mgs) = max( 0.0, qfwet(mgs))

       IF ( incwet == 1 .and. qfwet(mgs) < qfdry(mgs) .and. dfwet(mgs) < dg0thresh ) THEN
        ! ELSE
!! || defined (WRFEXTRAS)
        !   IF ( dfwet(mgs) < dg0thresh ) THEN
             ! find portion of qc and qr collection that are dry/wet growth for d > dwet

               ratio = Min( maxratiolu, dfwet(mgs)/xdia(mgs,lf,1) )
               
               tmp1 = gaminterp(ratio,alpha(mgs,lf),13,1) ! alpha + 3
               tmp2 = gaminterp(ratio,alpha(mgs,lf),12,1) ! alpha + 2
               tmp3 = gaminterp(ratio,alpha(mgs,lf), 9,1) ! alpha + 1

            IF ( qfacw(mgs)*dtp > qxmin(lf) ) THEN
              vt = abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))

          qxacwtmp = 0.25*pi*efw(mgs)*cx(mgs,lf)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab1lh(mgs,lf,lc)*xdia(mgs,lf,3)*xdia(mgs,lc,3) +    &
     &            tmp3*da1lc(mgs)*xdia(mgs,lc,3)**2 )
             ELSE
              qxacwtmp = 0.0
             ENDIF

            IF ( qfacr(mgs)*dtp > qxmin(lf) ) THEN

       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,lr,1) )

         qxacrtmp = 0.25*pi*efr(mgs)*cx(mgs,lf)*qx(mgs,lr)*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab1lh(mgs,lf,lr)*xdia(mgs,lf,3)*xdia(mgs,lr,3) +    &
     &            tmp3*da1lr(mgs)*xdia(mgs,lr,3)**2 )
             ELSE
               qxacrtmp = 0.0
             ENDIF

        x = gaminterp(ratio,alpha(mgs,lf),9,1) ! alpha + 1
        y = gaminterp(ratio,alpha(mgs,lf),3,1) ! alpha + b/2 + 5/2
        
        hxventtmp =  0.78*x + y*fwventy(mgs)  !   &

        qxacitmp = 0.0
        IF ( qfaci(mgs)*dtp > qxmin(lf) ) THEN
              vt = abs(vtxbar(mgs,lf,1)-vtxbar(mgs,li,1))

          qxacitmp = 0.25*pi*eficlsn(mgs)*cx(mgs,lf)*qx(mgs,li)*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab1lh(mgs,lf,li)*xdia(mgs,lf,3)*xdia(mgs,li,3) +    &
     &            tmp3*da1(li)*xdia(mgs,li,3)**2 )

          cxacitmp =    &
     &        0.25*pi*eficlsn(mgs)*cx(mgs,lf)*cx(mgs,li)*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab0lh(mgs,lf,li)*xdia(mgs,lf,3)*xdia(mgs,li,3) +    &
     &            tmp3*da0(li)*xdia(mgs,li,3)**2 )
        ENDIF

        qxacstmp = 0.0
        IF ( qfacs(mgs)*dtp > qxmin(lf) ) THEN
              vt = abs(vtxbar(mgs,lf,1)-vtxbar(mgs,ls,1))

          qxacstmp = 0.25*pi*efsclsn(mgs)*cx(mgs,lf)*qx(mgs,ls)*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab1lh(mgs,lf,ls)*xdia(mgs,lf,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da1(ls)*xdia(mgs,ls,3)**2 )

          cxacstmp = 0.25*pi*efsclsn(mgs)*cx(mgs,lf)*cx(mgs,ls)*vt*   &
     &         (  tmp1*da0lf(mgs)*xdia(mgs,lf,3)**2 +     &
     &            tmp2*dab0lh(mgs,lf,ls)*xdia(mgs,lf,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da0(ls)*xdia(mgs,ls,3)**2 )
        ENDIF

             qxwettmp =   &
     &       xdia(mgs,lf,1)*hxventtmp*cx(mgs,lf)*fwet1(mgs)   &
     &     + fwet2(mgs)*(qxacitmp + qxacstmp)

          
!          qfwet(mgs) = qfacw(mgs) + qfacr(mgs) - qxacwtmp - qxacrtmp + qxwettmp
          ! as dry growth but subtract part for D > Dw and add wet growth for D > Dw
          qfwet(mgs) = qfacw(mgs) + qfacr(mgs) + qfaci(mgs) + qfacs(mgs) &
                        - efi(mgs)*qxacitmp - efs(mgs)*qxacstmp          &
                        -  qxacwtmp - qxacrtmp + qxwettmp

          qfaci(mgs) =  qfaci(mgs) + (1.0 - efi(mgs))*qxacitmp
          qfacs(mgs) =  qfacs(mgs) + (1.0 - efs(mgs))*qxacstmp
          cfaci(mgs) =  cfaci(mgs) + (1.0 - efi(mgs))*cxacitmp
          cfacs(mgs) =  cfacs(mgs) + (1.0 - efs(mgs))*cxacstmp

         ! qfacw(mgs) = Min( qfacw(mgs), 0.5*qx(mgs,lc)*dtpinv )
           
       !    ELSE
       !      qfwet(mgs) = qfdry(mgs)
       !    ENDIF
         ENDIF ! incwet
      ENDIF ! lf

       qhlwet(mgs) = 0.0
       IF ( lhl .gt. 1 ) THEN
         !IF ( incwet == 0 ) THEN
         qhlwet(mgs) =   &
     &     ( xdia(mgs,lhl,1)*hlvent(mgs)*cx(mgs,lhl)*fwet1(mgs)   &
     &     + fwet2(mgs)*(qhlaci(mgs) + qhlacs(mgs)) )
         qhlwet(mgs) = max( 0.0, qhlwet(mgs))
         
         IF ( incwet == 1 .and. qhlwet(mgs) < qhldry(mgs) .and. dhlwet(mgs) < dg0thresh ) THEN
         !ELSE
!! || defined (WRFEXTRAS)
         !  IF ( dhlwet(mgs) < dg0thresh ) THEN
             ! find portion of qc and qr collection that are dry/wet growth for d > dwet

               ratio = Min( maxratiolu, dhlwet(mgs)/xdia(mgs,lhl,1) )
               
               tmp1 = gaminterp(ratio,alpha(mgs,lhl),13,2) ! alpha + 3
               tmp2 = gaminterp(ratio,alpha(mgs,lhl),12,2) ! alpha + 2
               tmp3 = gaminterp(ratio,alpha(mgs,lhl), 9,2) ! alpha + 1

            IF ( qhlacw(mgs)*dtp > qxmin(lhl) ) THEN
              vt = abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))

          qxacwtmp = 0.25*pi*ehlw(mgs)*cx(mgs,lhl)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab1lh(mgs,lhl,lc)*xdia(mgs,lhl,3)*xdia(mgs,lc,3) +    &
     &            tmp3*da1lc(mgs)*xdia(mgs,lc,3)**2 )
             ELSE
              qxacwtmp = 0.0
             ENDIF

            IF ( qhlacr(mgs)*dtp > qxmin(lhl) ) THEN

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,lr,1) )

         qxacrtmp = 0.25*pi*ehlr(mgs)*cx(mgs,lhl)*qx(mgs,lr)*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab1lh(mgs,lhl,lr)*xdia(mgs,lhl,3)*xdia(mgs,lr,3) +    &
     &            tmp3*da1lr(mgs)*xdia(mgs,lr,3)**2 )
             ELSE
               qxacrtmp = 0.0
             ENDIF

        x = gaminterp(ratio,alpha(mgs,lhl),9,2) ! alpha + 1
        y = gaminterp(ratio,alpha(mgs,lhl),3,2) ! alpha + b/2 + 5/2
        
        hxventtmp =  0.78*x + y*hlventy(mgs)  !   &

        qxacitmp = 0.0
        IF ( qhlaci(mgs)*dtp > qxmin(lhl) ) THEN
              vt = abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))

          qxacitmp = 0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*qx(mgs,li)*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab1lh(mgs,lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            tmp3*da1(li)*xdia(mgs,li,3)**2 )

          cxacitmp =    &
     &        0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*cx(mgs,li)*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab0lh(mgs,lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            tmp3*da0(li)*xdia(mgs,li,3)**2 )

        ENDIF

        qxacstmp = 0.0
        IF ( qhlacs(mgs)*dtp > qxmin(lhl) ) THEN
              vt = abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))

          qxacstmp = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*qx(mgs,ls)*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab1lh(mgs,lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da1(ls)*xdia(mgs,ls,3)**2 )

          cxacstmp = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*cx(mgs,ls)*vt*   &
     &         (  tmp1*da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            tmp2*dab0lh(mgs,lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            tmp3*da0(ls)*xdia(mgs,ls,3)**2 )
        ENDIF

             qxwettmp =   &
     &       xdia(mgs,lhl,1)*hxventtmp*cx(mgs,lhl)*fwet1(mgs)   &
     &     + fwet2(mgs)*(qxacitmp + qxacstmp)

          ! qhlacw(mgs) + qhlacr(mgs) - qxacwtmp - qxacrtmp is the 'dry' growth
          ! at smaller diameters
!          qhlwet(mgs) = qhlacw(mgs) + qhlacr(mgs) - qxacwtmp - qxacrtmp + qxwettmp
          ! as dry growth but subtract part for D > Dw and add wet growth for D > Dw
          qhlwet(mgs) = qhlacw(mgs) + qhlacr(mgs) + qhlaci(mgs) + qhlacs(mgs) &
                        - ehli(mgs)*qxacitmp - ehls(mgs)*qxacstmp          &
                        -  qxacwtmp - qxacrtmp + qxwettmp

          qhlaci(mgs) =  qhlaci(mgs) + (1.0 - ehli(mgs))*qxacitmp
          qhlacs(mgs) =  qhlacs(mgs) + (1.0 - ehls(mgs))*qxacstmp
          chlaci(mgs) =  chlaci(mgs) + (1.0 - ehli(mgs))*cxacitmp
          chlacs(mgs) =  chlacs(mgs) + (1.0 - ehls(mgs))*cxacstmp

        !   ELSE
        !     qhlwet(mgs) = qhldry(mgs)
        !   ENDIF
         ENDIF ! incwet
       ENDIF
       
       ELSE ! ( tfrdry < temg(mgs) .and. temg(mgs) < tfr )
       
        qhwet(mgs) = qhdry(mgs)
        qhlwet(mgs) = qhldry(mgs)
        qfwet(mgs) = qfdry(mgs)
       ENDIF ! } ( tfrdry < temg(mgs) .and. temg(mgs) < tfr )
!
!      qhlwet(mgs) = qhldry(mgs)

      end do

!
! shedding rate
!
      qsshr(:)  =  0.0
      qhshr(:)  =  0.0
      qhlshr(:) =  0.0
      qhshh(:)  =  0.0
      csshr(:)  =  0.0
      csshrr(:) = 0.0
      chshr(:)  =  0.0
      chlshr(:)  =  0.0
      chshrr(:)  =  0.0
      chlshrr(:)  =  0.0
      vhshdr(:)  = 0.0
      vhlshdr(:)  = 0.0
      wetsfc(:)  = .false.
      wetgrowth(:)  = .false.
      wetsfchl(:)  = .false.
      wetgrowthhl(:)  = .false.
      wetsfcf(:)  = .false.
      qfshr(:)  =  0.0
      vfshdr(:)  = 0.0
      wetgrowthf(:) = .false.
      cfshr(:)  =  0.0
      cfshrr(:)  =  0.0

      do mgs = 1,ngscnt
!
!
!
      qhshr(mgs)  = Min( 0.0, qhwet(mgs) - qhdry(mgs) )  ! water that freezes should never be more than what sheds
      


      qhlshr(mgs)  =  Min( 0.0, qhlwet(mgs) - qhldry(mgs) )

      IF ( lf > 1 ) THEN
       qfshr(mgs)  = Min( 0.0, qfwet(mgs) - qfdry(mgs) )  ! water that freezes should never be more than what sheds
      ENDIF
!
! limit wet growth to only higher density particles
!
      qsshr(mgs)  =  0.0
!
!
!  no shedding for temperatures < 243.15 
!
      if ( temg(mgs) .lt. 243.15 ) then
       qsshr(mgs)  =  0.0
       qhshr(mgs)  =  0.0
       qhlshr(mgs) =  0.0
       vhshdr(mgs)  = 0.0
       vhlshdr(mgs)  = 0.0
       wetsfc(mgs) = .false.
       wetgrowth(mgs) = .false.
       wetsfchl(mgs) = .false.
       wetgrowthhl(mgs) = .false.
       wetsfcf(mgs)  = .false.
       qfshr(mgs)  =  0.0
       vfshdr(mgs)  = 0.0
       wetgrowthf(mgs) = .false.
      end if
!
!  shed all at temperatures > 273.15
!
      if ( temg(mgs) .gt. tfr ) then

       IF ( .false. ) THEN ! old and incorrect -- Thanks to Shaofeng Hua for noticing this error (9/17/2017)
       qsshr(mgs)  = -qsdry(mgs)
       qhshr(mgs)  = -qhdry(mgs)
       qhlshr(mgs) = -qhldry(mgs)
       IF ( lf > 0 ) qfshr(mgs)  = -qfdry(mgs)
       ELSE ! new and correct
       ! note that the qxacr terms should be zero here, so shedding at T > 0 is all from the droplets
       qsshr(mgs)   = - qsacr(mgs) - qsacw(mgs) ! -qsdry(mgs)
       qhlshr(mgs)  = - qhlacw(mgs) - qhlacr(mgs) ! -qhldry(mgs)
       qhshr(mgs)  = - qhacw(mgs) - qhacr(mgs) ! -qhdry(mgs)
       IF ( lf > 0 ) qfshr(mgs)  = - qfacw(mgs) - qfacr(mgs)

       ENDIF

       vhshdr(mgs)  = -vhacw(mgs) - vhacr(mgs)
       vhlshdr(mgs)  = -vhlacw(mgs) - vhlacr(mgs)
       qhwet(mgs)  = 0.0
       qhlwet(mgs) = 0.0
       IF ( lf > 0 ) THEN 
         vfshdr(mgs)  = -vfacw(mgs) - vfacr(mgs)
         qfwet(mgs)  = 0.0
       ENDIF
      end if
!
!      if (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr  ) THEN
        wetsfc(mgs) =  (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qhmlr(mgs) < -qxmin(lh) .and.  temg(mgs) > tfr )
        wetgrowth(mgs) = (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
!      ENDIF
      if (qfshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) THEN
        wetsfcf(mgs) = (qfshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qfmlr(mgs) < -qxmin(lf) .and.  temg(mgs) > tfr )
        wetgrowthf(mgs) = (qfshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
!        wetsfcf(mgs) = .false. ! (qfshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qfmlr(mgs) < -qxmin(lf) .and.  temg(mgs) > tfr )
!        wetgrowthf(mgs) = .false. ! (qfshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
      ENDIF
      if (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) THEN
        wetsfchl(mgs) = (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qhlmlr(mgs) < -qxmin(lhl) .and.  temg(mgs) > tfr )
        wetgrowthhl(mgs) = (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
      ENDIF

      end do
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      csshr(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*Min(0.0,qsshr(mgs))
       
       chshr(mgs) = 0.0 ! no change to graupel number concentration for wet-growth shedding
       
      !   tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
        ! Base the drop size on the shedding regime
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chshrr(mgs) =  rho0(mgs)*qhshr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
      
      
      
      chlshr(mgs) = 0.0
      chlshrr(mgs) = 0.0
      IF ( lhl .gt. 1 ) THEN 
!         chlshr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlshr(mgs)


       chlshr(mgs) = 0.0 ! no change to hail number concentration for wet-growth shedding
       
      !   tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
        ! Base the drop size on the shedding regime
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chlshrr(mgs) =  rho0(mgs)*qhlshr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 

      ENDIF ! ( lhl > 1 )

      IF ( lf > 1 ) THEN 

        cfshr(mgs) = 0.0 ! no change to number concentration for wet-growth shedding
!        cfshrr(mgs) =  rho0(mgs)*qfshr(mgs)/(xdn(mgs,lr)*Min(vshd,vshdgs(mgs,lf)))  ! into rain 
        cfshrr(mgs) =  rho0(mgs)*qfshr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lf))  ! into rain 
       
       ENDIF ! lf > 1
      
      end do
      end if



! count up number of times that wet growth occurs:
      DO mgs = 1,ngscnt

      IF ( ( wetgrowth(mgs) .or. wetgrowthhl(mgs) & 
     &   .or. wetgrowthf(mgs)  &
     &     ) .and. temg(mgs) .lt. 273.15 ) THEN

       iwetg = iwetg + 1
       IF ( qx(mgs,lh) .gt. 0.5e-3 ) iwetg1 = iwetg1 + 1

      
      END IF


      END DO
!
!  final decisions
!
      do mgs = 1,ngscnt
!
!  Snow
!
      if ( qsshr(mgs) .lt. 0.0 ) then
      qsdpv(mgs) = 0.0
      qssbv(mgs) = 0.0
      else
      qsshr(mgs) = 0.0
      end if
!
!     if ( qsdry(mgs) .lt. qswet(mgs) ) then
!     qswet(mgs) = 0.0
!     else
!     qsdry(mgs) = 0.0
!     end if
!

!  graupel
!
!
      if ( wetgrowth(mgs) .or. (mixedphase .and. fhw(mgs) .gt. 0.05 .and. temg(mgs) .gt. 243.15) ) then
      

! soaking (when not advected liquid water film with graupel)

        IF ( lvol(lh) .gt. 1 .and. .not. mixedphase) THEN
        ! rescale volumes to maximum density
         IF ( iwetsoak ) THEN 

         rimdn(mgs,lh) = xdnmx(lh)
         raindn(mgs,lh) = xdnmx(lh)
         vhacw(mgs) = qhacw(mgs)*rho0(mgs)/rimdn(mgs,lh)
         vhacr(mgs) = qhacr(mgs)*rho0(mgs)/raindn(mgs,lh)
!        IF ( lvol(lh) .gt. 1 .and. wetgrowth(mgs) ) THEN
         IF ( xdn(mgs,lh) .lt. xdnmx(lh) ) THEN
         ! soak some liquid into the graupel
!           v1 = xdnmx(lh)*vx(mgs,lh)/(xdn(mgs,lh)*dtp) ! volume available for filling
           v1 = (1. - xdn(mgs,lh)/xdnmx(lh))*vx(mgs,lh)/(dtp) ! volume available for filling
!            tmp = (vx(mgs,lh)/rho0(mgs))*(xdnmx(lh) - xdn(mgs,lh)) ! max mixing ratio of liquid water that can be added
           v2 = rho0(mgs)*qhwet(mgs)/xdnmx(lh)  ! volume of frozen accretion
           
           vhsoak(mgs) = Min(v1,v2)

           
         ENDIF
         
         ENDIF

         vhshdr(mgs) = Min(0.0, rho0(mgs)*qhwet(mgs)/xdnmx(lh) - vhacw(mgs) - vhacr(mgs) )
         
        ELSEIF ( lvol(lh) .gt. 1  .and. mixedphase ) THEN
!         vhacw(mgs) = rho0(mgs)*qhacw(mgs)/xdn0(lr)
!         vhacr(mgs) = rho0(mgs)*qhacr(mgs)/xdn0(lr)
        ENDIF
        

      qhdpv(mgs) = 0.0
!      qhsbv(mgs) = 0.0
      chdpv(mgs) = 0.0
!      chsbv(mgs) = 0.0

      IF ( ipelec > 0 ) THEN
      schacs(mgs) = 0.0
      schaci(mgs) = 0.0
      IF ( lis > 1 ) schacis(mgs) = 0.0
      schacw(mgs) = 0.0
      ENDIF
! collection efficiency modification

      IF ( ehi(mgs) .gt. 0.0 ) THEN
        IF ( incwet == 0 ) THEN
        qhaci(mgs) = Min(qimxd(mgs),qhaci0(mgs))  ! effectively sets collection eff to 1
        chaci(mgs) = Min(cimxd(mgs),chaci0(mgs))  ! effectively sets collection eff to 1
        ENDIF
      ENDIF
      IF ( ehs(mgs) .gt. 0.0 ) THEN
!        qhacs(mgs) = Min(qsmxd(mgs),qhacs(mgs)/ehs(mgs))  ! effectively sets collection eff to 1
        IF ( incwet == 0 ) THEN
        qhacs(mgs) = Min(qsmxd(mgs),qhacs0(mgs)) !/ehs(mgs)                   ! divide out the collection efficiency
        chacs(mgs) = Min(csmxd(mgs),chacs0(mgs)) !/ehs(mgs)                   ! divide out the collection efficiency
        qhacs(mgs) = Min(qsmxd(mgs),qhacs(mgs))   ! plug it back in
        ENDIF
        ehs(mgs) = ehsmax ! 1.0 ! min(ehsfrac*ehs(mgs),ehsmax)            ! modify it
      ENDIF

! be sure to catch particles with wet surfaces but not in wet growth to turn off Hallett-Mossop
      wetsfc(mgs) = .true.

      else
!        qhshr(mgs) = 0.0
      end if
!
!
!  hail
!
!      if ( lhl .gt. 1 .and. qhlshr(mgs) .lt. 0.0 ) then
      if ( lhl > 1 .and. ( wetgrowthhl(mgs) .or. (mixedphase .and. fhlw(mgs) .gt. 0.05 .and. temg(mgs) .gt. 243.15) ) ) then
!      if ( wetgrowthhl(mgs) ) then
       

      qhldpv(mgs) = 0.0
!      qhlsbv(mgs) = 0.0
      chldpv(mgs) = 0.0
!      chlsbv(mgs) = 0.0


      IF ( ipelec > 0 ) THEN
      schlacs(mgs) = 0.0
      schlaci(mgs) = 0.0
      schlacw(mgs) = 0.0
      ENDIF


        IF ( lvol(lhl) .gt. 1  .and. .not. mixedphase ) THEN
!        IF ( lvol(lhl) .gt. 1 .and. wetgrowthhl(mgs) ) THEN

         IF ( iwetsoak ) THEN 

         rimdn(mgs,lhl) = xdnmx(lhl) 
         raindn(mgs,lhl) = xdnmx(lhl) 
         vhlacw(mgs) = qhlacw(mgs)*rho0(mgs)/rimdn(mgs,lhl)
         vhlacr(mgs) = qhlacr(mgs)*rho0(mgs)/raindn(mgs,lhl)

         IF ( xdn(mgs,lhl) .lt. xdnmx(lhl) ) THEN
         ! soak some liquid into the hail
!           v1 = xdnmx(lhl)*vx(mgs,lhl)/(xdn(mgs,lhl)*dtp) ! volume available for filling
           v1 = (1. - xdn(mgs,lhl)/xdnmx(lhl))*vx(mgs,lhl)/(dtp) ! volume available for filling
!            tmp = (vx(mgs,lhl)/rho0(mgs))*(xdnmx(lhl) - xdn(mgs,lhl)) ! max mixing ratio of liquid water that can be added
           v2 = rho0(mgs)*qhlwet(mgs)/xdnmx(lhl)  ! volume of frozen accretion
           IF ( v1 > v2 ) THEN ! all the frozen stuff fits in
             vhlsoak(mgs) = v2
           ELSE  ! fill up the available space
             vhlsoak(mgs) = v1
           ENDIF
!           vhlacw(mgs) = 0.0
!           vhlacr(mgs) = Max( 0.0, v2 - v1 )
         ELSE
           vhlsoak(mgs) = 0.0
!           vhlacw(mgs) = 0.0
!           vhlacr(mgs) = rho0(mgs)*qhlwet(mgs)/raindn(mgs,lhl)
         
         ENDIF
         
         ENDIF

         vhlshdr(mgs) = Min(0.0, rho0(mgs)*qhlwet(mgs)/xdnmx(lhl) - vhlacw(mgs) - vhlacr(mgs) )


        ELSEIF ( lvol(lhl) .gt. 1  .and. mixedphase ) THEN
!         vhlacw(mgs) = rho0(mgs)*qhlacw(mgs)/xdn0(lr)
!         vhlacr(mgs) = rho0(mgs)*qhlacr(mgs)/xdn0(lr)
        ENDIF

      IF ( ehli(mgs) .gt. 0.0 .and. incwet == 0 ) THEN
        qhlaci(mgs) = Min(qimxd(mgs),qhlaci0(mgs))  ! effectively sets collection eff to 1
        chlaci(mgs) = Min(cimxd(mgs),chlaci0(mgs))  ! effectively sets collection eff to 1
      ENDIF

!      IF ( ehls(mgs) .gt. 0.0 ) THEN
!        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs(mgs)/ehls(mgs))
!      ENDIF
      IF ( ehls(mgs) .gt. 0.0 .and. incwet == 0 ) THEN
        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs0(mgs)) !/ehls(mgs)                   ! divide out the collection efficiency
        chlacs(mgs) = Min(csmxd(mgs),chlacs0(mgs)) !/ehls(mgs)                   ! divide out the collection efficiency
        ehls(mgs) = ehsmax ! 1.0 ! min(ehsfrac*ehs(mgs),ehsmax)            ! modify it
!        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs(mgs))   ! plug it back in
      ENDIF

      
!      qhlwet(mgs) = 1.0

! be sure to catch particles with wet surfaces but not in wet growth to turn off Hallett-Mossop
      wetsfchl(mgs) = .true.


      else
!      qhlshr(mgs) = 0.0
!      qhlwet(mgs) = 0.0
      end if

!
!
!  frozen drops
!
      if ( lf > 1 .and. ( wetgrowthf(mgs) .or. (mixedphase .and. ffw(mgs) .gt. 0.05 .and. temg(mgs) .gt. 243.15) ) ) then

      qfdpv(mgs) = 0.0
      cfdpv(mgs) = 0.0


      IF ( ipelec > 0 ) THEN
      scfacs(mgs) = 0.0
      scfaci(mgs) = 0.0
      scfacw(mgs) = 0.0
      ENDIF


        IF ( lvol(lf) .gt. 1  .and. .not. mixedphase ) THEN

         rimdn(mgs,lf) = xdnmx(lf) 
         raindn(mgs,lf) = xdnmx(lf) 
         vfacw(mgs) = qfacw(mgs)*rho0(mgs)/rimdn(mgs,lf)
         vfacr(mgs) = qfacr(mgs)*rho0(mgs)/raindn(mgs,lf)

         IF ( xdn(mgs,lf) .lt. xdnmx(lf) ) THEN
         ! soak some liquid into the ice
           v1 = (1. - xdn(mgs,lf)/xdnmx(lf))*vx(mgs,lf)/(dtp) ! volume available for filling
           v2 = rho0(mgs)*qfwet(mgs)/xdnmx(lf)  ! volume of frozen accretion
           IF ( v1 > v2 ) THEN ! all the frozen stuff fits in
             vfsoak(mgs) = v2
           ELSE  ! fill up the available space
             vfsoak(mgs) = v1
           ENDIF
         ELSE
           vfsoak(mgs) = 0.0
         
         ENDIF

         vfshdr(mgs) = Min(0.0, rho0(mgs)*qfwet(mgs)/xdnmx(lf) - vfacw(mgs) - vfacr(mgs) )


        ELSEIF ( lvol(lf) .gt. 1  .and. mixedphase ) THEN

        ENDIF

      IF ( efi(mgs) .gt. 0.0  .and. incwet == 0 ) THEN
        qfaci(mgs) = Min(qimxd(mgs),qfaci0(mgs))  ! effectively sets collection eff to 1
        cfaci(mgs) = Min(cimxd(mgs),cfaci0(mgs))  ! effectively sets collection eff to 1
      ENDIF

      IF ( efs(mgs) .gt. 0.0 .and. incwet == 0 ) THEN
        qfacs(mgs) = Min(qsmxd(mgs),qfacs0(mgs)) !/efs(mgs)                   ! divide out the collection efficiency
        cfacs(mgs) = Min(csmxd(mgs),cfacs0(mgs)) !/efs(mgs)                   ! divide out the collection efficiency
        efs(mgs) = ehsmax ! 1.0 ! min(ehsfrac*ehs(mgs),ehsmax)            ! modify it
      ENDIF

      
! be sure to catch particles with wet surfaces but not in wet growth to turn off Hallett-Mossop
      wetsfcf(mgs) = .true.

      end if
      end do
!
! Ice -> graupel conversion
!
      DO mgs = 1,ngscnt
      
      qhcni(mgs) = 0.0
      chcni(mgs) = 0.0
      chcnih(mgs) = 0.0
      vhcni(mgs) = 0.0
      
      IF ( iglcnvi .ge. 1 ) THEN
      IF ( temg(mgs) .lt. 273.0 .and. qiacw(mgs) - qidpv(mgs) .gt. 0.0 ) THEN
      
        
        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,li,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
        tmp = Min( Max( rimc3, tmp ), 900.0 )
        
        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)
        
!        write(0,*)  'rime dens = ',tmp
        
        IF ( tmp .ge. 200.0 .or. iglcnvi >= 2 ) THEN
          r = Max( 0.5*(xdn(mgs,li) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcni(mgs) = (qiacw(mgs) - qidpv(mgs)) ! *float(iglcnvi)
          chcni(mgs) = cx(mgs,li)*qhcni(mgs)/qx(mgs,li)
!          chcnih(mgs) = rho0(mgs)*qhcni(mgs)/(1.6e-10)
          chcnih(mgs) = Min(chcni(mgs), rho0(mgs)*qhcni(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcni(mgs) = rho0(mgs)*qhcni(mgs)/r
        ENDIF
      
      ELSEIF ( iglcnvi == 3 ) THEN

       IF ( temg(mgs) .lt. 273.0 .and. qiacw(mgs)*dtp > 2.*qxmin(lh) .and. gamice73fac*xmas(mgs,li) > xdnmn(lh)*xvmn(lh) ) THEN
      
        
        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,li,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
        tmp = Min( Max( rimc3, tmp ), 900.0 )
        
        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)
        
!        write(0,*)  'rime dens = ',tmp
        ! convert to particles with the mass of the mass-weighted diameter
      !  massofmwr = gamice73fac*xmas(mgs,li)
        
        IF ( tmp .ge. xdnmn(lh)  ) THEN
          r = Max( 0.5*(xdn(mgs,li) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcni(mgs) = 0.5*qiacw(mgs)
          chcni(mgs) = qhcni(mgs)/(gamice73fac*xmas(mgs,li))
          chcnih(mgs) = Min(chcni(mgs), rho0(mgs)*qhcni(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcni(mgs) = rho0(mgs)*qhcni(mgs)/r
        ENDIF
      
      ENDIF

      
      ENDIF
      ENDIF
      
      
      ENDDO
      
      
      qhlcnh(:) = 0.0
      chlcnh(:) = 0.0
      chlcnhhl(:) = 0.0
      vhlcnh(:) = 0.0
      vhlcnhl(:) = 0.0
      zhlcnh(:) = 0.0

      qhcnhl(:) = 0.0
      chcnhl(:) = 0.0
      vhcnhl(:) = 0.0
      zhcnhl(:) = 0.0

      qhlcnf(:) = 0.0
      chlcnf(:) = 0.0
      chlcnfhl(:) = 0.0
      fddenfac(:) = 1.0
      vhlcnf(:) = 0.0
      vhlcnfhl(:) = 0.0
      zhlcnf(:) = 0.0

      IF ( lhl .gt. 1  ) THEN
      
      IF ( ihlcnh == 1 .or. ihlcnh == 3 ) THEN

!
!  Graupel (h) conversion to hail (hl) based on Milbrandt and Yau 2005b
!
      DO mgs = 1,ngscnt

!        IF ( lhl .gt. 1 .and. ipconc .ge. 5 .and. qx(mgs,lh) .gt. 1.0e-3 .and.
!     :        xdn(mgs,lh) .gt. 750. .and. qhshr(mgs) .lt. 0.0 .and.
!     :        xdia(mgs,lh,3) .gt. 1.e-3 ) THEN
        IF ( hlcnhdia > 0 ) THEN
          ltest = xdia(mgs,lh,3) .gt. hlcnhdia  ! test on mean volume diameter
        ELSE 
!          ltest =  xdia(mgs,lh,1)*(3. + alpha(mgs,lh)) > Abs( hlcnhdia ) ! test on maximum mass diameter
          ltest =  xdia(mgs,lh,1)*(4. + alpha(mgs,lh)) > Abs( hlcnhdia ) ! test on mass-weighted diameter
        ENDIF

    ! if incwet > 0, then should use dhwet here to avoid calculating again
         IF ( iusedw == 0 .and. ihlcnh == 1 ) THEN
           dg0(mgs) = -1.
         ELSE
         IF ( temg(mgs) .le. tfr+hailcnvtoffset .and. &
              (( (qhacw(mgs) + qhacr(mgs))*dtp > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin &
               .and.  temg(mgs) .gt. dwtempmin ) .or. ( wetgrowth(mgs) .and. qx(mgs,lh) > hlcnhqmin )) ) THEN
!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*ehw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
          IF ( incwet > 0 ) THEN
            d = dhwet(mgs)
          ELSE
            ! First guess for dwet (not that good, but it is something)
            x =   1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            IF ( x > 1.e-20 ) THEN
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            ELSE
              dwr = 1.e30
            ENDIF
            d = Min(dwr, dg0thresh + 0.0001)
           IF ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) THEN
                      sqrtrhovt = Sqrt( rhovt(mgs) )
                      fventh = sqrtrhovt*(fpndl(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5) 
                      fventm = sqrtrhovt*(fschm(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5)
                      ltemq = (tfr-163.15)/fqsat+1.5
                      qvs0 = pqs(mgs)*tabqvs(ltemq)
                      denomdp = felf(mgs) + fcw(mgs)*temcg(mgs)
                      denominvdp = 1.d0/(felf(mgs) + fcw(mgs)*temcg(mgs))

!                      write(91,*) 'dw,dwr,temcg = ',100.*dw,100.*dwr,temcg(mgs)
                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = ehi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, ehw(mgs))*qx(mgs,lc) 
                      h4 = ehr(mgs)* qx(mgs,lr)
                      ! iterate to find minimum diameter for wet growth. Start with value of dwr
                      DO n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lh)*d**bxx(mgs,lh) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       IF ( x2 > 1.4 ) THEN
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       ELSE
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       ENDIF

                       IF ( .false. ) THEN ! this option includes 'am' separate from ah, which makes only small differences. Otherwise equivalent to second option
                        x1 = fventm*sqrtrhovt*Sqrt(d*vth)
                        IF ( x1 > 1.4 ) THEN
                          am = 0.78 + 0.308*x1 ! mass ventillation (Beard and Pruppacher 1971, eq. 8)
                        ELSE
                          am = 1.0 + 0.108*x1**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                        ENDIF
                        
                        d = 8.*denominvdp*( am*felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qvs0 - qx(mgs,lv)) - ah*ftka(mgs)*temcg(mgs)  )/ &
                           (dtp* ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs) +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2*denominvdp))
                       
                        ELSE

                        ! Based on Farley and Orville (1986), eq. 5-9 but neglecting the Ci*(T0-Ts) term in (8) since we want Ts=T0
                        ! Simplified mass rates as dm_w/dt = pi/4*d**2*(Vh - Vc)*rhoair*qc*ehw, etc.
                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)
                            
                        ENDIF
                        IF ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) EXIT
                        
                      ENDDO
                      
                      d = Min( d, dg0thresh + 0.0001 )
              ENDIF ! dwr < 0.2 .and. dwr > 0.0
              ENDIF ! incwet
              
             ! dg0(mgs) = Min( dwmax, Max( d, dwmin ) )
              dg0(mgs) = Max( d, dwmin )
          ELSE
         !   IF ( qx(mgs,lh) > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin .and. temg(mgs) .le. tfr+hailcnvtoffset  ) THEN
         !     dg0(mgs) = dwmax
         !   ELSE
              dg0(mgs) = dg0thresh + 0.0001
          !  ENDIF
          ENDIF
          
            IF ( ihlcnh == 3 .and. (qhacw(mgs) + qhacr(mgs))*dtp > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin &
                   .and. temg(mgs) .le. tfr+hailcnvtoffset .and. temg(mgs) > 238.0 ) THEN
           ! set a secondary condition on to capture large graupel that is riming but not in wet growth
!                dg0(mgs) = Min( dg0(mgs), dg0thresh - 0.0001 )
                dg0(mgs) = Min( dg0(mgs), dwmax )
            ENDIF
            
          ENDIF

        wtest = (dg0(mgs) > 0.0 .and. dg0(mgs) < dg0thresh )
        
        IF ( ihlcnh == 1 ) THEN ! .or. iusedw == 0  THEN
        
        IF ( ( wetgrowth(mgs) .and. (xdn(mgs,lh) .gt. hldnmn .or. lvh < 1 ) .and.  & ! correct this when hail gets turned on
     &        rimdn(mgs,lh) .gt. 800. .and.   &
     &        ltest .and. qx(mgs,lh) .gt. hlcnhqmin ) .or. wtest ) THEN ! {
!     :        xdia(mgs,lh,3) .gt. 2.e-3 .and. qx(mgs,lh) .gt. 1.0e-3  THEN ! 0823.2008 erm test
!        IF ( xdia(mgs,lh,3) .gt. 1.e-3 ) THEN
        IF ( qhacw(mgs) .gt. 0.0 .and. qhacw(mgs) .gt. qhaci(mgs) .and. temg(mgs) .le. tfr+hailcnvtoffset ) THEN ! {
        ! dh0 is the diameter dividing wet growth from dry growth (Ziegler 1985), modified by MY05
!          dh0 = 0.01*(exp(temcg(mgs)/(1.1e4*(qx(mgs,lc)+qx(mgs,lr)) - 
!     :           1.3e3*qx(mgs,li) + 1.0e-3 ) ) - 1.0)
          IF ( wtest ) THEN
            dh0 = dg0(mgs)
          ELSE
            x = (1.1e4*(rho0(mgs)*qx(mgs,lc)) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0e-3 )
            IF ( x > 1.e-20 ) THEN
            arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
            dh0 = 0.01*(exp(arg) - 1.0)
            ELSE
             dh0 = 1.e30
            ENDIF
            dg0(mgs) = Min(dh0, dg0thresh + 0.0001)
          ENDIF ! wtest
!          dh0 = Max( dh0, 5.e-3 )
          
!         IF ( dh0 .gt. 0.0 ) write(0,*) 'dh0 = ',dh0
!         IF ( dh0 .gt. 1.0e-4 ) THEN
         IF ( xdia(mgs,lh,3)/dh0 .gt. 0.1 ) THEN !{
!         IF ( xdia(mgs,lh,3) .lt. 20.*dh0 .and. dh0 .lt. 2.0*xdia(mgs,lh,3) ) THEN 
           tmp = qhacw(mgs) + qhacr(mgs) + qhaci(mgs) + qhacs(mgs)
!           qtmp = Min( 1.0, xdia(mgs,lh,3)/(2.0*dh0) )*(tmp)
           qtmp = Min( 100.0, xdia(mgs,lh,3)/(2.0*dh0) )*(tmp)
           qhlcnh(mgs) = Min(  qxmxd(mgs,lh), qtmp )
           
           IF ( ipconc .ge. 5 ) THEN !{
!           dh0 = Max( xdia(mgs,lh,3), Min( dh0, 5.e-3 ) ) ! do not create hail greater than 5mm diam. unless the graupel is larger
           IF ( .not. wtest ) dh0 = Min( dh0, 10.e-3 ) ! do not create hail greater than 10mm diam., which is the max graupel size
           IF ( qx(mgs,lhl) > 0.1e-3 ) dh0 = Max( dh0, xdia(mgs,lhl,3) ) ! when enough hail is established, do not dilute the size
           chlcnhhl(mgs) = Min( cxmxd(mgs,lh), rho0(mgs)*qhlcnh(mgs)/(pi*xdn(mgs,lh)*dh0**3/6.0) )

           r = rho0(mgs)*qhlcnh(mgs)/(xdn(mgs,lh)*xv(mgs,lh))  ! number of graupel particles at mean volume diameter
           chlcnh(mgs) = Max( chlcnhhl(mgs), r )
           ENDIF !}
           
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))

          ENDIF !}

        ENDIF ! }
        ENDIF ! }
        
        ELSEIF ( ihlcnh == 3 ) THEN !{
         

          IF ( wtest  .and. &
               ( qhacw(mgs)*dtp > qxmin(lh) .and. temg(mgs) .lt. tfr+hailcnvtoffset .and. qx(mgs,lh) > hlcnhqmin ) ) THEN
        ! convert number, mass, and reflectivity for d > dw
           IF ( ipconc == 5 ) THEN
            ! dg0(mgs) = Min( dg0(mgs), hldia1 )
             !dg0(mgs) = hldia1
           ENDIF
           
           ratio = Min( maxratiolu, dg0(mgs)/xdia(mgs,lh,1) )


           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lh),4,1)
           IF ( ipconc == 5 ) THEN
       !      tmp2 = Min( 0.25, tmp2 )
           ENDIF
            qxd1 = qx(mgs,lh)*(tmp2)
            qhlcnh(mgs) = dtpinv*qxd1
            flim = 1.0
            tmp3 = qxmxd(mgs,lh)
            IF (qxd1 > tmp3 ) THEN
!              flim = tmp3/(qxd1)
!              qhlcnh(mgs) = flim*qhlcnh(mgs)
            ENDIF

            
            
            IF ( ( qxd1 > qxmin(lhl) .and. ipconc > 5 ) .or. ( qxd1 > 10.*qxmin(lhl) .and. ipconc == 5) ) THEN
            
           ! number
            tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
             IF ( ipconc == 5 ) THEN
          !     tmp = Min( 0.2, tmp )
             ENDIF
            cxd1 = flim*cx(mgs,lh)*( tmp)
            chlcnh(mgs) = dtpinv*cxd1
            chlcnhhl(mgs) = chlcnh(mgs)

           IF ( qx(mgs,lhl) > qxmin(lhl) .and. dmhlopt > 0 ) THEN
             tmp = rho0(mgs)*qhlcnh(mgs)/chlcnhhl(mgs)
             IF ( tmp < xmas(mgs,lhl) ) THEN
               ! dh0 = ( qxd1*dh0 + qx(mgs,lhl)*xmas(mgs,lhl))/( qxd1 + qx(mgs,lhl))  ! weighted average
               dh0 = (( qxd1*tmp**(1./3.) + qx(mgs,lhl)*xmas(mgs,lhl)**(1./3.))/( qxd1 + qx(mgs,lhl)))**3  ! weighted average
               chlcnhhl(mgs) = Min( chlcnhhl(mgs), rho0(mgs)*qhlcnh(mgs)/dh0 )
             ELSE
!               dh0 = Max( dh0, xmas(mgs,lhl) ) ! when enough hail is established, do not dilute the size
             ENDIF
           ENDIF


           ! reflectivity
           IF ( ipconc >= 6 .and. lzh > 1 .and. lzhl > 1 ) THEN
            tmp3 = gaminterp(ratio,alpha(mgs,lh),11,1)
            zxd1 = flim*zx(mgs,lh)*(tmp3)
            zhlcnh(mgs) = dtpinv*zxd1

              ! tmp4 is the Z from the converted particles assuming shape of alphamax
              IF ( icorrecthaildbz >= 1 .and. zxd1 > 10.*zxmin ) THEN
              tmp3 = g1xmax*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lh)/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new hail number to match zxd1
                ! increase cxd1 to make z,q,c rates consistent
                ! cxd1 = g1xmax*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lh)/6.0)**2)
                cxd1 = tmp3/zxd1
                chlcnhhl(mgs) = dtpinv*cxd1
              ENDIF
              ENDIF
           ELSE
            zxd1 = 0
           ENDIF
           IF ( ipconc == 5 .and. icorrecthaildbz >= 1 ) THEN ! Adjust cxd1 by reflectivity removed from graupel
            tmp3 = gaminterp(ratio,alpha(mgs,lh),11,1)
            ! tmp5 is graupel reflectivity moment
            tmp5 = g1x(mgs,lh)*(rho0(mgs)*qx(mgs,lh))**2/((pi*xdn(mgs,lh)/6.)**2*cx(mgs,lh))
            zxd1 = flim*(tmp3)*tmp5
            IF ( zxd1 > 10.*zxmin ) THEN
            ! tmp4 is the reflectivity of the newly-converted graupel particles (use g1x(lh) for loss term)
            ! which we want to match zxd1 to prevent spurious increase in total reflectivity
              tmp3 =  g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lh)/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new hail number to match zxd1
                ! cxd1 = g1x(mgs,lhl)*(rho0(mgs)*qxd1)**2/(zxd1*pi*xdn(mgs,lh)/6.0) ! trial form results in tiny hail
            ! want the adjust size of the new hail so that Z is conserved, so increase number of
            ! particles to make  qxd1,zxd1, and C consistent.
            ! want zxd1 = g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/(c*(pi*xdn(mgs,lh)/6.0)**2)
                 ! Use g1x(mgs,lh) here instead of g1x(mgs,lhl) because rzxhlh will then multiply
                 ! by g1x(mgs,lhl)/g1x(mgs,lh)
               ! cxd1 = g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lh)/6.0)**2)
                cxd1 = tmp3/zxd1
                chlcnhhl(mgs) = dtpinv*cxd1 ! multiplied later by rzxhlh(mgs)
              ENDIF
              ENDIF
           ENDIF

            ELSE
               qhlcnh(mgs) = 0.0
            ENDIF

           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
           
           ENDIF


        ENDIF !}
      
      ENDDO
      
      ELSEIF ( ihlcnh == 2 ) THEN ! 10-ice type conversion 

!
! Staka and Mansell (2005) type conversion
!
!      hldia1 is set in micro_module and namelist
!      IF ( .true. ) THEN
      
        ! convert number, mass, and reflectivity for d > hldia1,
        ! regardless of wet growth status, but as long as riming > 0
        DO mgs = 1,ngscnt
        IF ( qhacw(mgs)*dtp > qxmin(lh) .and. temg(mgs) .lt. tfr+hailcnvtoffset .and. qx(mgs,lh) > qxmin(lh) ) THEN
           ratio = Min( maxratiolu, hldia1/xdia(mgs,lh,1) )

           ! number
            tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
            cxd1 = cx(mgs,lh)*( tmp)
            chlcnh(mgs) = dtpinv*cxd1
            chlcnhhl(mgs) = chlcnh(mgs)

           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lh),4,1)
            qxd1 = qx(mgs,lh)*(tmp2)
            qhlcnh(mgs) = dtpinv*qxd1

           ! reflectivity
           IF ( lzh > 1 .and. lzhl > 1 ) THEN
            tmp3 = gaminterp(ratio,alpha(mgs,lh),11,1)
            zxd1 = zx(mgs,lh)*(tmp3)
            zhlcnh(mgs) = dtpinv*zxd1
           ELSE
            zxd1 = 0
           ENDIF
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
           
         ENDIF
         
         ENDDO
!        ENDIF
      ELSEIF ( ihlcnh == 0 ) THEN

      do mgs = 1,ngscnt
!      qhlcnh(mgs) = 0.0
!      chlcnh(mgs) = 0.0
      if ( wetgrowth(mgs) .and. temg(mgs) .lt. tfr-5. .and. qx(mgs,lh) > qxmin(lh) ) then
      if ( qhacw(mgs).gt.1.e-6 .and. ( xdn(mgs,lh) > 700. .or. lvh == 0 ) ) then
      qhlcnh(mgs) =                                                   &
        ((pi*xdn(mgs,lh)*cx(mgs,lh)) / (6.0*rho0(mgs)*dtp))           &
       *exp(-hldia1/xdia(mgs,lh,1))                                    &
       *( (hldia1**3) + 3.0*(hldia1**2)*xdia(mgs,lh,1)                  &
        + 6.0*(hldia1)*(xdia(mgs,lh,1)**2) + 6.0*(xdia(mgs,lh,1)**3) )
      qhlcnh(mgs) =   min(qhlcnh(mgs),qhmxd(mgs))
      IF ( ipconc .ge. 5 ) THEN
        chlcnh(mgs) = Min( cxmxd(mgs,lh), cx(mgs,lh)*Exp(-hldia1/xdia(mgs,lh,1)))
        chlcnhhl(mgs) = chlcnh(mgs)
!        chlcnh(mgs) = Min( cxmxd(mgs,lh), rho0(mgs)*qhlcnh(mgs)/(2.0*xmas(mgs,lh) ))
      ENDIF
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
      end if
      end if
      end do
      
!      ENDIF ! true
      
      ENDIF ! ihlcnh options
      
     ! convert low-density hail to graupel
      IF ( icvhl2h >= 1 ) THEN
      DO mgs = 1,ngscnt
        IF (  qx(mgs,lhl) > qxmin(lhl) .and. xdn(mgs,lhl) < 0.5*(xdnmn(lhl) + xdnmx(lhl)) ) THEN
          tmp = Min(0.95, 1. - 0.5*(1. + tanh(0.125*(xdn(mgs,lhl) - 1.01*xdnmn(lhl) )) ))
          qhcnhl(mgs) = tmp*qx(mgs,lhl)*dtpinv
          chcnhl(mgs) = cx(mgs,lhl)*qhcnhl(mgs)/qx(mgs,lhl)
          vhcnhl(mgs) = vx(mgs,lhl)*qhcnhl(mgs)/qx(mgs,lhl)
          
        ENDIF
      ENDDO
      
      ENDIF
      
      ENDIF ! lhl > 1

      IF ( lf > 1 .and. lhl > 1 ) THEN !{
      IF ( ihlcnh == 1 .or. ihlcnh == 3 ) THEN

!
!  frozen drop (f) conversion to hail (hl) based on wet growt
!
      DO mgs = 1,ngscnt

        IF ( hlcnhdia > 0 ) THEN
          ltest = xdia(mgs,lf,3) .gt. hlcnhdia  ! test on mean volume diameter
        ELSE 
          ltest =  xdia(mgs,lf,1)*(4. + alpha(mgs,lf)) > Abs( hlcnhdia ) ! test on mass-weighted diameter
        ENDIF

         IF ( iusedw == 0 .and. ihlcnh == 1 ) THEN
           df0(mgs) = -1.
         ELSE
         IF ( temg(mgs) .le. tfr+hailcnvtoffset .and. &
              (( (qfacw(mgs) + qfacr(mgs))*dtp > qxmin(lf) .and. qx(mgs,lf) > hlcnhqmin &
               .and.  temg(mgs) .gt. dwtempmin ) .or. ( wetgrowth(mgs) .and. qx(mgs,lf) > hlcnhqmin )) ) THEN
!         IF (((qfacw(mgs) + qfacr(mgs))*dtp > qxmin(lf) .and. qx(mgs,lf) > hlcnhqmin .and. temg(mgs) .le. tfr+hailcnvtoffset  &
!               .and.  temg(mgs) .gt. dwtempmin ) .or. ( wetgrowthf(mgs) .and. qx(mgs,lf) > hlcnhqmin ) ) THEN

!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*efw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(efw(mgs)*qx(mgs,lc)+efr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
          IF ( incwet > 0 ) THEN
            d = dfwet(mgs)
          ELSE
            x =   1.1e4 * rho0(mgs)*(efw(mgs)*qx(mgs,lc)+efr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            IF ( x > 1.e-20 ) THEN
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            ELSE
              dwr = 1.e30
            ENDIF
          d = dwr
           IF ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) THEN
                      sqrtrhovt = Sqrt( rhovt(mgs) )
                      fventh = sqrtrhovt*(fpndl(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5) 
                      fventm = sqrtrhovt*(fschm(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5)
                      ltemq = (tfr-163.15)/fqsat+1.5
                      qvs0 = pqs(mgs)*tabqvs(ltemq)
                      denomdp = felf(mgs) + fcw(mgs)*temcg(mgs)
                      denominvdp = 1.d0/(felf(mgs) + fcw(mgs)*temcg(mgs))

!                      write(91,*) 'dw,dwr,temcg = ',100.*dw,100.*dwr,temcg(mgs)
                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = efi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, efw(mgs))*qx(mgs,lc) 
                      h4 = efr(mgs)* qx(mgs,lr)
                      DO n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lf)*d**bxx(mgs,lf) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       IF ( x2 > 1.4 ) THEN
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       ELSE
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       ENDIF

                       IF ( .false. ) THEN ! this option includes 'am' separate from ah, which makes only small differences. Otherwise equivalent to second option
                        x1 = fventm*sqrtrhovt*Sqrt(d*vth)
                        IF ( x1 > 1.4 ) THEN
                          am = 0.78 + 0.308*x1 ! mass ventillation (Beard and Pruppacher 1971, eq. 8)
                        ELSE
                          am = 1.0 + 0.108*x1**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                        ENDIF
                        
                        d = 8.*denominvdp*( am*felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qvs0 - qx(mgs,lv)) - ah*ftka(mgs)*temcg(mgs)  )/ &
                           (dtp* ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs) +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2*denominvdp))
                       
                        ELSE

                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)
                        ENDIF
                        
                        IF ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) EXIT
                        
                      ENDDO
              ENDIF
              ENDIF ! incwet
              
              df0(mgs) = Min( dwmax, Max( d, dwmin ) )

          ELSE
            IF ( qx(mgs,lf) > qxmin(lf) .and. qx(mgs,lf) > hlcnhqmin .and. temg(mgs) .le. tfr+hailcnvtoffset  ) THEN
              df0(mgs) = dwmax
            ELSE
              df0(mgs) = dg0thresh + 0.0001
            ENDIF
!            df0(mgs) = dg0thresh + 0.0001
          ENDIF

           IF ( ihlcnh == 3 .and. (qfacw(mgs) + qfacr(mgs))*dtp > qxmin(lf) .and. qx(mgs,lf) > hlcnhqmin &
                   .and. temg(mgs) .le. tfr+hailcnvtoffset ) THEN
           ! set a secondary condition on to capture large graupel that is riming but not in wet growth
                df0(mgs) = Min( df0(mgs), dwmax ) ! dg0thresh - 0.0001 )
            ENDIF
          
          ENDIF

        wtest = (df0(mgs) > 0.0 .and. df0(mgs) < dg0thresh )
        
        IF ( ihlcnh == 1 ) THEN ! .or. iusedw == 0  THEN
        
        IF ( ( wetgrowthf(mgs) .and. (xdn(mgs,lf) .gt. hldnmn .or. lvf < 1 ) .and.  & ! correct this when hail gets turned on
     &        rimdn(mgs,lf) .gt. 800. .and.   &
     &        ltest .and. qx(mgs,lf) .gt. hlcnhqmin ) .or. wtest ) THEN ! {
!     :        xdia(mgs,lf,3) .gt. 2.e-3 .and. qx(mgs,lf) .gt. 1.0e-3  THEN ! 0823.2008 erm test
!        IF ( xdia(mgs,lf,3) .gt. 1.e-3 ) THEN
        IF ( qfacw(mgs) .gt. 0.0 .and. qfacw(mgs) .gt. qfaci(mgs) .and. temg(mgs) .le. tfr+hailcnvtoffset ) THEN ! {
        ! dh0 is the diameter dividing wet growth from dry growth (Ziegler 1985), modified by MY05
!          dh0 = 0.01*(exp(temcg(mgs)/(1.1e4*(qx(mgs,lc)+qx(mgs,lr)) - 
!     :           1.3e3*qx(mgs,li) + 1.0e-3 ) ) - 1.0)
          IF ( wtest ) THEN
            dh0 = df0(mgs)
          ELSE
            x = (1.1e4*(rho0(mgs)*qx(mgs,lc)) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0e-3 )
            IF ( x > 1.e-20 ) THEN
            arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
            dh0 = 0.01*(exp(arg) - 1.0)
            ELSE
             dh0 = 1.e30
            ENDIF
          ENDIF ! wtest

         IF ( xdia(mgs,lf,3)/dh0 .gt. 0.1 ) THEN !{

           tmp = qfacw(mgs) + qfacr(mgs) + qfaci(mgs) + qfacs(mgs)

           qtmp = Min( 100.0, xdia(mgs,lf,3)/(2.0*dh0) )*(tmp)

           qhlcnf(mgs) = Min(  qxmxd(mgs,lf), qtmp )
           
           IF ( ipconc .ge. 5 ) THEN !{
!           dh0 = Max( xdia(mgs,lf,3), Min( dh0, 5.e-3 ) ) ! do not create hail greater than 5mm diam. unless the graupel is larger
           IF ( .not. wtest ) dh0 = Min( dh0, 10.e-3 ) ! do not create hail greater than 10mm diam., which is the max graupel size
           IF ( qx(mgs,lhl) > 0.1e-3 ) dh0 = Max( dh0, xdia(mgs,lhl,3) ) ! when enough hail is established, do not dilute the size
           chlcnfhl(mgs) = Min( cxmxd(mgs,lf), rho0(mgs)*qhlcnf(mgs)/(pi*xdn(mgs,lf)*dh0**3/6.0) )

           r = rho0(mgs)*qhlcnf(mgs)/(xdn(mgs,lf)*xv(mgs,lf))  ! number of graupel particles at mean volume diameter
           chlcnf(mgs) = Max( chlcnfhl(mgs), r )

           ENDIF !}
           
           vhlcnf(mgs) = rho0(mgs)*qhlcnf(mgs)/xdn(mgs,lf)
           vhlcnfhl(mgs) = rho0(mgs)*qhlcnf(mgs)/Max(xdnmn(lhl), xdn(mgs,lf))

          ENDIF !}

        ENDIF ! }
        ENDIF ! }
        
        ELSEIF ( ihlcnh == 3 ) THEN !{
         

          IF ( wtest  .and. &
               ( qfacw(mgs)*dtp > qxmin(lf) .and. temg(mgs) .lt. tfr+hailcnvtoffset .and. qx(mgs,lf) > hlcnhqmin ) ) THEN
        ! convert number, mass, and reflectivity for d > dw
           IF ( ipconc == 5 ) THEN
            ! df0(mgs) = Min( df0(mgs), hldia1 )
             !df0(mgs) = hldia1
           ENDIF
           
           ratio = Min( maxratiolu, df0(mgs)/xdia(mgs,lf,1) )


           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lf),4,1)
           IF ( ipconc == 5 ) THEN
       !      tmp2 = Min( 0.25, tmp2 )
           ENDIF
            qxd1 = qx(mgs,lf)*(tmp2)
            qhlcnf(mgs) = dtpinv*qxd1
            flim = 1.0
            tmp3 = qxmxd(mgs,lf)
            IF (qxd1 > tmp3 ) THEN
!              flim = tmp3/(qxd1)
!              qhlcnf(mgs) = flim*qhlcnf(mgs)
            ENDIF


            
            
            IF ( ( qxd1 > qxmin(lhl) .and. ipconc > 5 ) .or. ( qxd1 > 10.*qxmin(lhl) .and. ipconc == 5) ) THEN
            
           ! number
            tmp = gaminterp(ratio,alpha(mgs,lf),1,1)
             IF ( ipconc == 5 ) THEN
          !     tmp = Min( 0.2, tmp )
             ENDIF
            cxd1 = flim*cx(mgs,lf)*( tmp)
            chlcnf(mgs) = dtpinv*cxd1
            chlcnfhl(mgs) = chlcnf(mgs)

           IF ( qx(mgs,lhl) > qxmin(lhl) .and. dmhlopt > 0 ) THEN
             tmp = rho0(mgs)*qhlcnf(mgs)/chlcnfhl(mgs) ! mean mass of converted particles
             IF ( tmp < xmas(mgs,lhl) ) THEN
               ! dh0 = ( qxd1*dh0 + qx(mgs,lhl)*xmas(mgs,lhl))/( qxd1 + qx(mgs,lhl))  ! weighted average
               dh0 = (( qxd1*tmp**(1./3.) + qx(mgs,lhl)*xmas(mgs,lhl)**(1./3.))/( qxd1 + qx(mgs,lhl)))**3  ! weighted average
               chlcnfhl(mgs) = Min( chlcnfhl(mgs), rho0(mgs)*qhlcnf(mgs)/dh0 )
             ELSE
!               dh0 = Max( dh0, xmas(mgs,lhl) ) ! when enough hail is established, do not dilute the size
             ENDIF
           ENDIF


           ! reflectivity
           IF ( lzf > 1 .and. lzhl > 1 ) THEN
            tmp = chlcnfhl(mgs) 
            tmp3 = gaminterp(ratio,alpha(mgs,lf),11,1)
            zxd1 = zx(mgs,lf)*(tmp3)
            zhlcnf(mgs) = flim*dtpinv*zxd1
              IF ( icorrecthaildbz >= 1 .and. zxd1 > 10.*zxmin) THEN
              ! tmp4 is the Z from the converted particles assuming shape of alphamax
              tmp3 = g1xmax*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lf)/6.0)**2)
              tmp4 = tmp3/cxd1 ! g1xmax*(rho0(mgs)*qxd1)**2/(cxd1*(pi*xdn(mgs,lf)/6.0)**2)
              IF ( tmp4 > zxd1 ) THEN ! calculate new hail number to match zxd1
                ! increase cxd1 to make z,q,c rates consistent
                ! cxd1 = g1xmax*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lf)/6.0)**2)
                cxd1 = tmp3/zxd1
                chlcnfhl(mgs) = dtpinv*cxd1
              ENDIF
              ENDIF
           ELSE
            zxd1 = 0
           ENDIF

           IF ( ipconc == 5 .and. icorrecthaildbz >= 1 ) THEN ! Adjust cxd1 by reflectivity removed from graupel
            tmp3 = gaminterp(ratio,alpha(mgs,lf),11,1)
            ! tmp5 is FD reflectivity moment (note that alphah is used for FD)
            tmp5 = g1x(mgs,lh)*(rho0(mgs)*qx(mgs,lf))**2/((pi*xdn(mgs,lf)/6.)**2*cx(mgs,lf))
            zxd1 = flim*(tmp3)*tmp5
            IF ( zxd1 > zxmin ) THEN
            ! tmp4 is the reflectivity of the newly-converted graupel particles (use g1x(lh) for loss term)
            ! which we want to match zxd1 to prevent spurious increase in total reflectivity
              tmp3 = g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/((pi*xdn(mgs,lf)/6.0)**2)
              tmp4 = tmp3/cxd1
              IF ( tmp4 > zxd1 ) THEN ! calculate new hail number to match zxd1
                ! cxd1 = g1x(mgs,lhl)*(rho0(mgs)*qxd1)**2/(zxd1*pi*xdn(mgs,lh)/6.0) ! trial form results in tiny hail
            ! want the adjust size of the new hail so that Z is conserved, so increase number of
            ! particles to make  qxd1,zxd1, and C consistent.
            ! want zxd1 = g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/(c*(pi*xdn(mgs,lh)/6.0)**2)
                 ! Use g1x(mgs,lh) here instead of g1x(mgs,lhl) because rzxhlh will then multiply
                 ! by g1x(mgs,lhl)/g1x(mgs,lh)
               ! cxd1 = g1x(mgs,lh)*(rho0(mgs)*qxd1)**2/(zxd1*(pi*xdn(mgs,lf)/6.0)**2)
                cxd1 = tmp3/zxd1
                chlcnfhl(mgs) = dtpinv*cxd1 ! multiplied later by rzxhlh(mgs)
              ENDIF
              ENDIF
           ENDIF

            ELSE
               qhlcnf(mgs) = 0.0
            ENDIF

           vhlcnf(mgs) = rho0(mgs)*qhlcnf(mgs)/xdn(mgs,lf)
           vhlcnfhl(mgs) = rho0(mgs)*qhlcnf(mgs)/Max(xdnmn(lhl), xdn(mgs,lf))
           
           ENDIF


        ENDIF !}
      
      ENDDO
      
      ELSEIF ( ihlcnh == 2 ) THEN ! 10-ice type conversion 

!
! Staka and Mansell (2005) type conversion
!
!      hldia1 is set in micro_module and namelist
!      IF ( .true. ) THEN
      
        ! convert number, mass, and reflectivity for d > hldia1,
        ! regardless of wet growth status, but as long as riming > 0
        DO mgs = 1,ngscnt
        IF ( qfacw(mgs)*dtp > qxmin(lf) .and. temg(mgs) .lt. tfr+hailcnvtoffset .and. qx(mgs,lf) > qxmin(lf) ) THEN
           ratio = Min( maxratiolu, hldia1/xdia(mgs,lf,1) )

           ! number
            tmp = gaminterp(ratio,alpha(mgs,lf),1,1)
            cxd1 = cx(mgs,lf)*tmp
            chlcnf(mgs) = dtpinv*cxd1
            chlcnfhl(mgs) = chlcnf(mgs)

           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lf),4,1)
            qxd1 = qx(mgs,lf)*tmp2
            qhlcnf(mgs) = dtpinv*qxd1

           ! reflectivity
           IF ( lzf > 1 .and. lzhl > 1 ) THEN
            tmp3 = gaminterp(ratio,alpha(mgs,lf),11,1)
            zxd1 = zx(mgs,lf)*(tmp3)
            zhlcnf(mgs) = dtpinv*zxd1
           ELSE
            zxd1 = 0
           ENDIF
           vhlcnf(mgs) = rho0(mgs)*qhlcnf(mgs)/xdn(mgs,lf)
           vhlcnfhl(mgs) = rho0(mgs)*qhlcnf(mgs)/Max(xdnmn(lhl), xdn(mgs,lf))
           
         ENDIF
         
         ENDDO
!        ENDIF
      ELSEIF ( ihlcnh == 0 ) THEN

      do mgs = 1,ngscnt
      if ( wetgrowthf(mgs) .and. temg(mgs) .lt. tfr-5. .and. qx(mgs,lf) > qxmin(lf) ) then
      if ( qfacw(mgs).gt.1.e-6 .and. xdn(mgs,lf) > 700. ) then
      qhlcnf(mgs) =                                                   &
        ((pi*xdn(mgs,lf)*cx(mgs,lf)) / (6.0*rho0(mgs)*dtp))           &
       *exp(-hldia1/xdia(mgs,lf,1))                                    &
       *( (hldia1**3) + 3.0*(hldia1**2)*xdia(mgs,lf,1)                  &
        + 6.0*(hldia1)*(xdia(mgs,lf,1)**2) + 6.0*(xdia(mgs,lf,1)**3) )
      qhlcnf(mgs) =   min(qhlcnf(mgs),qxmxd(mgs,lf))
      IF ( ipconc .ge. 5 ) THEN
        chlcnf(mgs) = Min( cxmxd(mgs,lf), cx(mgs,lf)*Exp(-hldia1/xdia(mgs,lf,1)))
        chlcnfhl(mgs) = chlcnf(mgs)
      ENDIF
           vhlcnf(mgs) = rho0(mgs)*qhlcnf(mgs)/xdn(mgs,lf)
           vhlcnfhl(mgs) = rho0(mgs)*qhlcnf(mgs)/Max(xdnmn(lhl), xdn(mgs,lf))
      end if
      end if
      end do
      
      ENDIF ! ihlcnh options
      
      ENDIF !} lf > 1
  
  

!
! Ziegler snow conversion to graupel
!
      DO mgs = 1,ngscnt

      qhcns(mgs) = 0.0
      chcns(mgs) = 0.0
      chcnsh(mgs) = 0.0
      vhcns(mgs) = 0.0

      qscnh(mgs) = 0.0
      cscnh(mgs) = 0.0
      vscnh(mgs) = 0.0

      IF ( ipconc .ge. 5 ) THEN

        ! test attempt at converting graupel to snow when not riming but growing by deposition
        IF ( temg(mgs) < tfr .and. qx(mgs,lh) .gt. qxmin(lh) .and. qhdpv(mgs) > qxmin(lh)*dtpinv  &
     &       .and. qhacw(mgs) < qxmin(lh)*dtpinv ) THEN
          IF ( xdn(mgs,lh) < 290. ) THEN
!          qscnh(mgs) = 2.*qhdpv(mgs)
!          cscnh(mgs) = cx(mgs,lh)*qscnh(mgs)/qx(mgs,lh)
!          vscnh(mgs) = rho0(mgs)*qscnh(mgs)/xdn(mgs,lh)
          ENDIF
        ENDIF


        IF ( qx(mgs,ls) .gt. qxmin(ls) .and. qsacw(mgs) .gt. 0.0 ) THEN

!      DATA VGRA/1.413E-2/  ! this is the volume (cm**3) of a 3mm diam. sphere
!    vgra = 1.4137e-8 m**3

!      DNNET=DNCNV-DNAGG
!      DQNET=QXCON+QSACC+SDEP
!
!      DNSCNV=EXP(-(ROS*XNS*VGRA/(RO*QI)))*((1.-(XNS*VGRA*ROS/
!     / (RO*QI)))*DNNET + (XNS**2*VGRA*ROS/(RO*QI**2))*DQNET)
!      IF(DNSCNV.LT.0.) DNSCNV=0.
!
!      QIHC=(ROS*VGRA/RO)*DNSCNV
!
!      QH=QH+DT*QIHC
!      QI=QI-DT*QIHC
!      XNH=XNH+DT*DNSCNV
!      XNS=XNS-DT*DNSCNV

        IF ( iglcnvs .eq. 1 ) THEN  ! Zrnic, Ziegler et al (1993)

        dnnet = cscnvis(mgs) + cscnis(mgs) - csacs(mgs)
        dqnet = qscnvi(mgs) + qscni(mgs) + qsacw(mgs) + qsdpv(mgs) + qssbv(mgs)

        a3 = 1./(rho0(mgs)*qx(mgs,ls))
        a1 = Exp( - xdn(mgs,ls)*cx(mgs,ls)*vgra*a3 ) !! EXP(-(ROS*XNS*VGRA/(RO*QI)))
! (1.-(XNS*VGRA*ROS/(RO*QI)))*DNNET
        a2 =  (1.-(cx(mgs,ls)*vgra*xdn(mgs,ls)*a3))*dnnet
! (XNS**2*VGRA*ROS/(RO*QI**2))*DQNET
        a4 = cx(mgs,ls)**2*vgra*xdn(mgs,ls)*a3/qx(mgs,ls)*dqnet

        chcns(mgs) = Max( 0.0, a1*(a2 + a4) )
        chcns(mgs) = Min( chcns(mgs), cxmxd(mgs,ls) )
        chcnsh(mgs) = chcns(mgs)

        qhcns(mgs) = Min( xdn(mgs,ls)*vgra*rhoinv(mgs)*chcns(mgs), qxmxd(mgs,ls) )
        vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),xdnmn(lh))
!        vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),400.)

        ELSEIF ( iglcnvs .ge. 2  ) THEN  ! treat like ice crystals, i.e., check for rime density (ERM)

          IF ( temg(mgs) .lt. 273.0 .and. ( qsacw(mgs) - qsdpv(mgs) .gt. 0.0 .or. &
              ( iglcnvs >= 3 .and. qsacw(mgs)*dtp > 2.*qxmin(lh) .and. gamsnow73fac*xmas(mgs,ls) > xdnmn(lh)*xvmn(lh)  ) ) ) THEN !{


        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,ls,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
!        tmp = Min( Max( rimc3, tmp ), 900.0 )
        tmp = Min( tmp , 900.0 )

        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)

!        write(0,*)  'rime dens = ',tmp

        IF ( iglcnvs == 2 ) THEN !{
        IF ( tmp .ge. 200.0  ) THEN
          r = Max( 0.5*(xdn(mgs,ls) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcns(mgs) = (qsacw(mgs) - qsdpv(mgs))
          chcns(mgs) = cx(mgs,ls)*qhcns(mgs)/qx(mgs,ls)
!          chcnih(mgs) = rho0(mgs)*qhcni(mgs)/(1.6e-10)
          chcnsh(mgs) = Min(chcns(mgs), rho0(mgs)*qhcns(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcns(mgs) = rho0(mgs)*qhcns(mgs)/r
        ENDIF
        
        ELSEIF ( iglcnvs == 3 ) THEN
 
         ! convert to particles with the mass of the mass-weighted diameter
      !  massofmwr = gamice73fac*xmas(mgs,li)
        
        IF ( tmp > xdnmn(lh) ) THEN
          r = Max( 0.5*(xdn(mgs,ls) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcns(mgs) = 0.5*qsacw(mgs)
          chcns(mgs) = qhcns(mgs)/(gamsnow73fac*xmas(mgs,ls))
          chcns(mgs) = Min( chcns(mgs), cx(mgs,ls)*qhcns(mgs)/qx(mgs,ls))
          chcnsh(mgs) = Min(chcns(mgs), rho0(mgs)*qhcns(mgs)/(r*xvmn(lh)) )
          vhcns(mgs) = rho0(mgs)*qhcns(mgs)/r
        ENDIF

        ENDIF !}

      ENDIF !}

        ENDIF


        ENDIF

       ELSE ! single moment lfo

        qhcns(mgs) = 0.001*ehscnv(mgs)*max((qx(mgs,ls)-6.e-4),0.0)
        qhcns(mgs) = min(qhcns(mgs),qxmxd(mgs,ls))
        IF ( lvol(lh) .ge. 1 ) vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),400.)

       ENDIF
      ENDDO
!
!
!  heat budget for rain---not all rain that collects ice can freeze
!
!
!
      if ( irwfrz .gt. 0 .and. .not. mixedphase) then
!
      do mgs = 1,ngscnt
!
!  compute total rain that freeze when it interacts with cloud ice
!
      qrztot(mgs) = qrfrz(mgs) + qiacr(mgs) + qsacr(mgs)
!
!  compute the maximum amount of rain that can freeze
!  Used to limit freezing to 4*qrmxd, but now allow all rain to freeze if possible
!
      qrzmax(mgs) =   &
     &  ( xdia(mgs,lr,1)*rwvent(mgs)*cx(mgs,lr)*fwet1(mgs) )
      qrzmax(mgs) = max(qrzmax(mgs), 0.0)
      qrzmax(mgs) = min(qrztot(mgs), qrzmax(mgs))
      qrzmax(mgs) = min(qx(mgs,lr)*dtpinv, qrzmax(mgs))

      IF ( temcg(mgs) < -30. ) THEN ! allow all to freeze if T < -30 because fwet becomes invalid (negative)
        qrzmax(mgs) = qx(mgs,lr)*dtpinv
      ENDIF
!      qrzmax(mgs) = min(4.*qrmxd(mgs), qrzmax(mgs))
!
!  compute the correction factor
!
!      IF ( qrztot(mgs) .gt. qxmin(lr) ) THEN
      IF ( qrztot(mgs) .gt. qrzmax(mgs) .and. qrztot(mgs) .gt. qxmin(lr) ) THEN
        qrzfac(mgs) = qrzmax(mgs)/(qrztot(mgs))
      ELSE
        qrzfac(mgs) = 1.0
      ENDIF
      qrzfac(mgs) = min(1.0, qrzfac(mgs))
!
      end do
!
!
! now correct the above sources
!
!
      do mgs = 1,ngscnt
      if ( temg(mgs) .le. 273.15 .and. qrzfac(mgs) .lt. 1.0 ) then
      qrfrz(mgs)   = qrzfac(mgs)*qrfrz(mgs)
      qrfrzs(mgs)  = qrzfac(mgs)*qrfrzs(mgs)
      qrfrzf(mgs)  = qrzfac(mgs)*qrfrzf(mgs)
      qiacr(mgs)   = qrzfac(mgs)*qiacr(mgs)
      qsacr(mgs)   = qrzfac(mgs)*qsacr(mgs)
      qiacrf(mgs)  = qrzfac(mgs)*qiacrf(mgs)
      qiacrs(mgs)  = qrzfac(mgs)*qiacrs(mgs)
      crfrz(mgs)   = qrzfac(mgs)*crfrz(mgs)
      crfrzf(mgs)  = qrzfac(mgs)*crfrzf(mgs)
      crfrzs(mgs)  = qrzfac(mgs)*crfrzs(mgs)
      ciacr(mgs)   = qrzfac(mgs)*ciacr(mgs)
      ciacrf(mgs)  = qrzfac(mgs)*ciacrf(mgs)
      ciacrs(mgs)  = qrzfac(mgs)*ciacrs(mgs)

!      IF ( lzh .gt. 1 ) THEN
!        zrfrzf(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.)) * &
!        ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
!      ENDIF
      
       vrfrzf(mgs)  = qrzfac(mgs)*vrfrzf(mgs)
       viacrf(mgs)  = qrzfac(mgs)*viacrf(mgs)
      end if
      end do
!
!
!
      end if
!
!
!
!  evaporation of rain
!
!
!
      qrcev(:) = 0.0
      crcev(:) = 0.0


      do mgs = 1,ngscnt
!
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

      qrcev(mgs) =   &
     &  fvce(mgs)*cx(mgs,lr)*rwvent(mgs)*rwcap(mgs)*evapfac
! this line to allow condensation on rain:
      IF ( rcond .eq. 1 ) THEN
        qrcev(mgs) = min(qrcev(mgs), qxmxd(mgs,lv))
! this line to have evaporation only:
      ELSE
        qrcev(mgs) = min(qrcev(mgs), 0.0)
      ENDIF

      qrcev(mgs) = max(qrcev(mgs), -qrmxd(mgs))
!      if ( temg(mgs) .lt. 273.15 ) qrcev(mgs) = 0.0
      IF ( qrcev(mgs) .lt. 0. .and. lnr > 1 ) THEN
!        qrcev(mgs) =   -qrmxd(mgs)
!        crcev(mgs) = (rho0(mgs)/(xmas(mgs,lr)+1.e-20))*qrcev(mgs)
        IF ( icrcev == 1 ) THEN
          crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)))*qrcev(mgs)
        ELSEIF ( icrcev == 2 ) THEN
          crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)))*qrcev(mgs)*vtxbar(mgs,lr,2)/vtxbar(mgs,lr,1)
        ELSE
          crcev(mgs) = 0.0
        ENDIF
      ELSE
         crcev(mgs) = 0.0
      ENDIF
!      if ( temg(mgs) .lt. 273.15 ) crcev(mgs) = 0.0
!
      ENDIF

      end do
!
! evaporation/condensation of wet graupel and snow
!
      IF ( lhwlg > 1 ) THEN
      qhcevlg(:) = 0.0
      chcevlg(:) = 0.0
      ENDIF
      IF ( lhlwlg > 1 ) THEN
      qhlcevlg(:) = 0.0
      chlcevlg(:) = 0.0
      ENDIF


!
!
!
!  ICE MULTIPLICATION: Two modes (rimpa, and rimpb)
!  (following Cotton et al. 1986)
!
 
      chmul1(:) =  0.0
      chlmul1(:) =  0.0
      csmul1(:) = 0.0
!
      qhmul1(:) =  0.0
      qhlmul1(:) =  0.0
      qsmul1(:) =  0.0
      IF ( lf > 1 ) qfmul1(:) =  0.0
      do mgs = 1,ngscnt
 
       ltest =  qx(mgs,lh) .gt. qxmin(lh)
       IF ( lhl > 1 )  ltest =  ltest .or. qx(mgs,lhl) .gt. qxmin(lhl)
       
      IF ( (itype1 .ge. 1 .or. itype2 .ge. 1 )   &
     &              .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
      if ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 271.15 ) then
       IF ( ipconc .ge. 2 ) THEN
        IF ( xv(mgs,lc) .gt. 0.0     &
     &     .and.  ltest &
!     .and. itype2 .ge. 2    &
     &       ) THEN
!
!  Ziegler et al. 1986 Hallett-Mossop process.  VSTAR = 7.23e-15 (vol of 12micron radius)
!
         IF ( alpha(mgs,lc) == 0.0 ) THEN
           ex1 = (1./250.)*Exp(-7.23e-15/xv(mgs,lc))
         ELSE
           
           ratio = (1. + alpha(mgs,lc))*(7.23e-15)/xv(mgs,lc)

           IF ( usegamxinfcnu ) THEN
            i = Nint(dgami*(1. + alpha(mgs,lc)))
            gcnup1 = gmoi(i)
            ex1 = (1./250.)*Gamxinf(1.+alpha(mgs,lc), ratio)/(gcnup1)
           ELSE
             ratio = Min( maxratiolu, ratio )
             tmp = gaminterp(ratio,alpha(mgs,lc),1,1)
             ex1 = (1./250.)*tmp
           ENDIF
         ENDIF
       IF ( itype2 .le. 2 ) THEN
         ft = Max(0.0,Min(1.0,-0.11*temcg(mgs)**2 - 1.1*temcg(mgs)-1.7))
       ELSE
        IF ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 267.15 ) THEN
          ft = 0.5
        ELSEIF (temg(mgs) .ge. 267.15 .and. temg(mgs) .le. 269.15 ) THEN
          ft = 1.0
        ELSEIF (temg(mgs) .ge. 269.15 .and. temg(mgs) .le. 271.15 ) THEN
          ft = 0.5
        ELSE 
          ft = 0.0
        ENDIF
       ENDIF
!        rhoinv = 1./rho0(mgs)
!        DNSTAR = ex1*cglacw(mgs)
        
       IF ( ft > 0.0 ) THEN
        
        IF ( itype2 > 0 ) THEN
         IF ( qx(mgs,lh) .gt. qxmin(lh) .and. (.not. wetsfc(mgs))  ) THEN
          chmul1(mgs) = ft*ex1*chacw(mgs)
!          chmul1(mgs) = Min( ft*ex1*chacw(mgs), ft*(30.*1.e+06)*rho0(mgs)*qhacw(mgs) ) ! 1.e+6 converts kg to mg; Saunders & Hosseini (2001) average of about 30 crystals per mg
          qhmul1(mgs) = cimas0*chmul1(mgs)*rhoinv(mgs)
         ENDIF
         IF ( lf > 1 ) THEN
         IF ( qx(mgs,lf) .gt. qxmin(lf) .and. (.not. wetsfcf(mgs))  ) THEN
!         IF ( qx(mgs,lf) .gt. qxmin(lf)  ) THEN
          cfmul1(mgs) = ft*ex1*cfacw(mgs)
!          cfmul1(mgs) = Min( ft*ex1*cfacw(mgs), ft*(30.*1.e+06)*rho0(mgs)*qfacw(mgs) ) ! 1.e+6 converts kg to mg; Saunders & Hosseini (2001) average of about 30 crystals per mg
          qfmul1(mgs) = cimas0*cfmul1(mgs)*rhoinv(mgs)
         ENDIF
         ENDIF
         IF ( lhl .gt. 1 ) THEN
           IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs))  ) THEN
            chlmul1(mgs) = (ft*ex1*chlacw(mgs))
            qhlmul1(mgs) = cimas0*chlmul1(mgs)*rhoinv(mgs)
           ENDIF
         ENDIF
        ENDIF ! itype2

        IF ( itype1 > 0 ) THEN
         IF ( qx(mgs,lh) .gt. qxmin(lh) .and. (.not. wetsfc(mgs))  ) THEN
          tmp = ft*(3.5e+08)*rho0(mgs)*qhacw(mgs)
          chmul1(mgs) = chmul1(mgs) + tmp
          qhmul1(mgs) = qhmul1(mgs) + cimas0*tmp*rhoinv(mgs)
         ENDIF
         IF ( lf .gt. 1 ) THEN
           IF ( qx(mgs,lf) .gt. qxmin(lf) .and. (.not. wetsfcf(mgs)) ) THEN
!           IF ( qx(mgs,lf) .gt. qxmin(lf) ) THEN
            tmp = ft*(3.5e+08)*rho0(mgs)*qfacw(mgs)
            cfmul1(mgs) = cfmul1(mgs) + tmp
            qfmul1(mgs) = qfmul1(mgs) + cimas0*tmp*rhoinv(mgs)
           ENDIF
         ENDIF
         IF ( lhl .gt. 1 ) THEN
           IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs)) ) THEN
            tmp = ft*(3.5e+08)*rho0(mgs)*qhlacw(mgs)
            chlmul1(mgs) = chlmul1(mgs) + tmp
            qhlmul1(mgs) = qhlmul1(mgs) + cimas0*tmp*rhoinv(mgs)
           ENDIF
         ENDIF
        ENDIF ! itype1

    ! Hallett and Saunders (1979)
        IF ( ipelec >= 2 .and. ihallettsaund > 0 ) THEN
          schaci(mgs) = schaci(mgs) + elecfac*xhallettsaund*chmul1(mgs)
          IF ( lhl > 1 ) THEN
            schlaci(mgs) = schlaci(mgs) + elecfac*xhallettsaund*chlmul1(mgs)
          ENDIF
          IF ( lf > 1 ) THEN
            scfaci(mgs) = scfaci(mgs) + elecfac*xhallettsaund*cfmul1(mgs)
          ENDIF
        ENDIF
        
        ENDIF ! ft

        ENDIF ! xv(mgs,lc) .gt. 0.0 .and.

       ELSE ! ipconc .lt. 2
!
!  define the temperature function
!
      fimt1(mgs) = 0.0
!
! Cotton et al. (1986) version
!
      if ( temg(mgs) .ge. 268.15 .and. temg(mgs) .le. 270.15 ) then
        fimt1(mgs) = 1.0 -(temg(mgs)-268.15)/2.0
      elseif (temg(mgs) .le. 268.15 .and. temg(mgs) .ge. 265.15 ) then
        fimt1(mgs) = 1.0 +(temg(mgs)-268.15)/3.0
      ELSE 
        fimt1(mgs) = 0.0
      end if
!
! Ferrier (1994) version
!
      if ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 267.15 ) then
        fimt1(mgs) = 0.5
      elseif (temg(mgs) .ge. 267.15 .and. temg(mgs) .le. 269.15 ) then
        fimt1(mgs) = 1.0
      elseif (temg(mgs) .ge. 269.15 .and. temg(mgs) .le. 271.15 ) then
        fimt1(mgs) = 0.5
      ELSE 
        fimt1(mgs) = 0.0
      end if
!
!
!   type I:  350 splinters are formed for every 1e-3 grams of cloud
!            water accreted by graupel/hail (note converted to MKS units)
!            3.5e+8 has units of 1/kg
!
      IF ( itype1 .ge. 1 ) THEN
       fimta(mgs) = (3.5e+08)*rho0(mgs)
      ELSE
       fimta(mgs) = 0.0
      ENDIF

!
!
!   type II:  1 splinter formed for every 250 cloud droplets larger than
!             24 micons in diameter (12 microns in radius) accreted by
!             graupel/hail
!
!
      fimt2(mgs) = 0.0
      xcwmas = xmas(mgs,lc) * 1000.
!
      IF ( itype2 .ge. 1 ) THEN
      if ( xcwmas.lt.1.26e-9 ) then
        fimt2(mgs) = 0.0
      end if
      if ( xcwmas .le. 3.55e-9 .and. xcwmas .ge. 1.26e-9 ) then
        fimt2(mgs) = (2.27)*alog(xcwmas) + 13.39
      end if
      if ( xcwmas .gt. 3.55e-9 ) then
        fimt2(mgs) = 1.0
      end if

      fimt2(mgs) = min(fimt2(mgs),1.0)
      fimt2(mgs) = max(fimt2(mgs),0.0)
      
      ENDIF
!
!     qhmul2 = 0.0
!     qsmul2 = 0.0
!
!     qhmul2 =
!    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qhacw(mgs)
!     qsmul2 =
!    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qsacw(mgs)
!
!      cimas0 = (1.0e-12)
!      cimas0 = 2.5e-10
      IF ( .not. wetsfc(mgs) ) THEN
      chmul1(mgs) =  fimt1(mgs)*(fimta(mgs) +   &
     &                           (4.0e-03)*fimt2(mgs))*qhacw(mgs)
      ENDIF
!
      qhmul1(mgs) =  chmul1(mgs)*(cimas0/rho0(mgs))

         IF ( lf .gt. 1 ) THEN
           IF ( qx(mgs,lf) .gt. qxmin(lf) .and. (.not. wetsfcf(mgs)) ) THEN
            tmp = fimt1(mgs)*(fimta(mgs) +   &
     &                           (4.0e-03)*fimt2(mgs))*qfacw(mgs)
            cfmul1(mgs) =  tmp
            qfmul1(mgs) = cimas0*tmp*rhoinv(mgs)
           ENDIF
         ENDIF
         IF ( lhl .gt. 1 ) THEN
           IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs)) ) THEN
            tmp = fimt1(mgs)*(fimta(mgs) +   &
     &                           (4.0e-03)*fimt2(mgs))*qhlacw(mgs)
            chlmul1(mgs) =  tmp
            qhlmul1(mgs) = cimas0*tmp*rhoinv(mgs)
           ENDIF
         ENDIF

!      qsmul1(mgs) =  csmul1(mgs)*(cimas0/rho0(mgs))
!
      ENDIF ! ( ipconc .ge. 2 )
      
      end if ! (in temperature range)
      
      ENDIF ! ( itype1 .eq. 1 .or. itype2 .eq. 1)
!
      end do
!
!
!
!     end if
!
!     end do
!
!
! ICE MULTIPLICATION FROM SNOW
!   Lo and Passarelli 82 / Willis and Heymsfield 89 / Schuur and Rutledge 00b
!   using kfrag as fragmentation rate (s-1) / 500 microns as char mean diam for max snow mix ratio
!
      csmul(:) = 0.0
      qsmul(:) = 0.0
      
      IF ( isnwfrac /= 0 ) THEN
      do mgs = 1,ngscnt
       IF (temg(mgs) .gt. 265.0) THEN !{
        if (xdia(mgs,ls,1) .gt. 100.e-6 .and. xdia(mgs,ls,1) .lt. 2.0e-3) then  ! equiv diameter 100microns to 2mm

        tmp = rhoinv(mgs)*pi*xdn(mgs,ls)*cx(mgs,ls)*(500.e-6)**3
        qsmul(mgs) = Max( kfrag*( qx(mgs,ls) - tmp ) , 0.0 )

        qsmul(mgs) = Min( qxmxd(mgs,li), qsmul(mgs) )
        csmul(mgs) = Min( cxmxd(mgs,li), rho0(mgs)*qsmul(mgs)/mfrag )

        endif
       ENDIF !}
      enddo
      ENDIF

!
!  frozen rain-rain interaction....
!
!
!
!
!  rain-ice interaction
!
!
      do mgs = 1,ngscnt
      qracif(mgs) = qraci(mgs)
!      cracif(mgs) = craci(mgs)
!      ciacrf(mgs) = ciacr(mgs)
      end do
!
! 
!  vapor to pristine ice crystals   UP
!
!
!
!  compute the nucleation rate
!
!     do mgs = 1,ngscnt
!     idqis = 0
!     if ( ssi(mgs) .gt. 1.0 ) idqis = 1
!     fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
!     dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/
!    >  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
!     qidsvp(mgs) = dqisdt(mgs)
!     cnnt = min(cnit*exp(-temcg(mgs)*bta1),1.0e+09)
!     qiint(mgs) = 
!    >  il5(mgs)*idqis*(1.0*dtpinv)
!    <  *min((6.88e-13)*cnnt/rho0(mgs), 0.25*dqisdt(mgs)) 
!     end do
!
!  Meyers et al. (1992; JAS) and Ferrier (1994) primary ice nucleation
!
      cmassin = cimas1  ! 6.88e-13
      do mgs = 1,ngscnt
      qiint(mgs) = 0.0
      ciint(mgs) = 0.0
      qicicnt(mgs) = 0.0
      cicint(mgs) = 0.0
      qipipnt(mgs) = 0.0
      cipint(mgs) = 0.0
      ccitmp = 0.0
      IF ( icenucopt == 1 .or. icenucopt == -10 .or. icenucopt == -11 ) THEN
      if ( ( temg(mgs) .lt. 268.15 .or.  &
!     : ( imeyers5 .and. temg(mgs) .lt.  273.0) ) .and.    &
     & ( imeyers5 .and. temg(mgs) .lt.  272.0 .and. temgkm2(mgs) .lt. tfr) ) .and.    &
     &    ciintmx .gt. (cx(mgs,li)+ccitmp)  &
!     :    .and. cninm(mgs) .gt. 0.   &
     &     ) then
      fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
      dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/   &
     &  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
!      qidsvp(mgs) = dqisdt(mgs)
      idqis = 0
      if ( ssi(mgs) .gt. 1.0 ) THEN
      idqis = 1 
      dzfacp = max( float(kgsp(mgs)-kgs(mgs)), 0.0 )
      dzfacm = max( float(kgs(mgs)-kgsm(mgs)), 0.0 )
      qiint(mgs) =   &
     &  idqis*il5(mgs)   &
     &  *(cmassin/rho0(mgs))   &
     &  *max(0.0,wvel(mgs))   &
     &  *max((cninp(mgs)-cninm(mgs)),0.0)*gz(kgs(mgs))   &
     &  /((dzfacp+dzfacm))

      qiint(mgs) = min(qiint(mgs), max(0.25*dqisdt(mgs),0.0)) 
      ciint(mgs) = qiint(mgs)*rho0(mgs)/cmassin
      
!
! limit new crystals so it does not increase the current concentration
!  above ciintmx 20,000 per liter (2.e7 per m**3)
!
!      ciintmx = 1.e9
!      ciintmx = 1.e9
      IF ( icenucopt /= -10 ) THEN
      
        IF ( lcin > 1 ) THEN
          ciint(mgs) = Min(ciint(mgs), ccin(mgs)*dtpinv) ! because ciint is a *rate*
          ccin(mgs) = ccin(mgs) - ciint(mgs)*dtp
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
        ELSEIF ( lcina > 1 ) THEN
          ciint(mgs) = Max(0.0, Min( ciint(mgs), Min( cnina(mgs), ciintmx ) - cina(mgs) ))
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
      
        ELSEIF ( icenucopt == 1 .and. ciint(mgs) .gt. Max(0.0, ciintmx - cx(mgs,li) - ccitmp )*dtpinv  ) THEN
          ciint(mgs) = Max(0.0, ciintmx - (cx(mgs,li)) )*dtpinv 
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

        ELSEIF ( icenucopt == -11 .and. dtp*ciint(mgs) .gt. ( cnina(mgs) - (cx(mgs,li) - ccitmp))) THEN
          ciint(mgs) = Max(0.0,  cnina(mgs) - (cx(mgs,li)+ccitmp)*dtpinv )
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

        ENDIF
      ENDIF
      
      end if
      endif

      ELSEIF ( icenucopt == 2 .or. icenucopt == -1 .or. icenucopt == -2 ) THEN
      
        IF ( ( temg(mgs) .lt. 268.15 .and. ssw(mgs) > 1.0 ) .or. ssi(mgs) > 1.25 ) THEN
          IF ( lcin > 1 ) THEN
           ciint(mgs) = Min(cnina(mgs), ccin(mgs))
           ciint(mgs) = Min( ciint(mgs), Max(0.0, ciintmx - (cx(mgs,li) - ccitmp) ) ) ! do not initiate ice beyond concentration of ciintmx
           ccin(mgs) = ccin(mgs) - ciint(mgs)
           ciint(mgs) = ciint(mgs)*dtpinv ! convert total initiation to a rate
          ELSE
           ciint(mgs) = Max( 0.0, cnina(mgs) - cina(mgs) )*dtpinv
          ENDIF
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

          fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
          dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/(1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
          qiint(mgs) = min(qiint(mgs), max(0.25*dqisdt(mgs),0.0))
          ciint(mgs) = qiint(mgs)*rho0(mgs)/cmassin
        ENDIF
      
      
      
      ELSEIF ( icenucopt == 3 .or. icenucopt == 4 .or. icenucopt == 10 ) THEN
        IF (  temg(mgs) .lt. 268.15 ) THEN
          IF ( lcin > 1 ) THEN
           ciint(mgs) = Min(cnina(mgs), ccin(mgs))
           ciint(mgs) = Min( ciint(mgs), Max(0.0, ciintmx - (cx(mgs,li) + ccitmp) ) ) ! do not initiate ice beyond concentration of ciintmx
           ccin(mgs) = ccin(mgs) - ciint(mgs)
           ciint(mgs) = ciint(mgs)*dtpinv ! convert total initiation to a rate
          ELSE
           ciint(mgs) = Max( 0.0, cnina(mgs) - cina(mgs) )*dtpinv
          ENDIF
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
        ENDIF

      ENDIF
!
      if ( xplate(mgs) .eq. 1 ) then
      qipipnt(mgs) = qiint(mgs)
      cipint(mgs) = ciint(mgs)
      end if
!
      if ( xcolmn(mgs) .eq. 1 ) then
      qicicnt(mgs) = qiint(mgs)
      cicint(mgs) = ciint(mgs)
      end if
!
!     qipipnt(mgs) = 0.0
!     qicicnt(mgs) = qiint(mgs)
!
      end do
!
! 

!
!  vapor to cloud droplets   UP
!
      if (ndebug .gt. 0 ) write(0,*) 'dbg = 8'
!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: set 3-component'
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
!  zero some arrays
!
!
      do mgs = 1,ngscnt
      qrshr(mgs) = 0.0
      qwshw(mgs) = 0.0
      cwshw(mgs) = 0.0
      qsshrp(mgs) = 0.0
      qhshrp(mgs) = 0.0
      end do
!
!
!  first sum all of the shed rain
!
!
      do mgs = 1,ngscnt
      qrshr(mgs) = qsshr(mgs) + qhshr(mgs) + qhlshr(mgs)
      crshr(mgs) = chshrr(mgs)/rzxh(mgs) + chlshrr(mgs)/rzxhl(mgs)
      
      IF ( lf > 1 ) THEN
       qrshr(mgs) = qrshr(mgs) + qfshr(mgs)
       crshr(mgs) = crshr(mgs) + cfshrr(mgs)/rzxf(mgs)
      ENDIF
      IF ( ished2cld == 1 ) THEN ! send shed liquid to droplets following Farley and Orville (1986)
        qwshw(mgs) = qrshr(mgs)
        cwshw(mgs) = rho0(mgs)*qwshw(mgs)/cwmas20
        
        crshr(mgs) = 0.0
        qrshr(mgs) = 0.0
        
      ELSEIF ( ished2cld == 2 .or. ( (ished2cld == 3 .or. ished2cld == 5) .and. temcg(mgs) < 0.0 )     &
                              .or. ( ished2cld == 4 .and. temcg(mgs) > 0.0 ) )  THEN
      ! instead of redirecting shed water, just reduce the collection of rain/cloud
      ! 2 : redirect at all temperatures
      ! 3 : redirect only for T < 0
      ! 4 : redirect only for T > 0
        IF ( qhshr(mgs) < 0.0 ) THEN
          IF ( -qhshr(mgs) <= qhacr(mgs) + qhacw(mgs) ) THEN
            IF ( ished2cld == 5 ) THEN
            ENDIF
            frac = -qhshr(mgs)/( qhacr(mgs) + qhacw(mgs) )
            qhacr(mgs) = frac*qhacr(mgs)
            qhacw(mgs) = frac*qhacw(mgs)
            chacr(mgs) = frac*chacr(mgs)
            chacw(mgs) = frac*chacw(mgs)
            chshrr(mgs) = 0.0 ! frac*chshrr(mgs)
            qhshr(mgs) = 0.0
          ELSE
           ! can this happen?
            qhshr(mgs) = 0.0 ! qhshr(mgs) + qhacr(mgs) + qhacw(mgs)
            qhacr(mgs) = 0.0
            qhacw(mgs) = 0.0
            chacr(mgs) = 0.0
            chacw(mgs) = 0.0
            chshrr(mgs) = 0.0
          ENDIF
          
        ENDIF
      
        IF ( lhl > 1 .and. qhlshr(mgs) < 0.0 ) THEN
          IF ( -qhlshr(mgs) <= qhlacr(mgs) + qhlacw(mgs) ) THEN
            frac = -qhlshr(mgs)/( qhlacr(mgs) + qhlacw(mgs) )
            qhlacr(mgs) = frac*qhlacr(mgs)
            qhlacw(mgs) = frac*qhlacw(mgs)
            chlacr(mgs) = frac*chlacr(mgs)
            chlacw(mgs) = frac*chlacw(mgs)
            qhlshr(mgs) = 0.0
            chlshrr(mgs) = 0.0 ! frac*chlshrr(mgs)
          ELSE
            qhlshr(mgs) = 0.0 ! qhlshr(mgs) + qhlacr(mgs) + qhlacw(mgs)
            qhlacr(mgs) = 0.0
            qhlacw(mgs) = 0.0
            chlacr(mgs) = 0.0
            chlacw(mgs) = 0.0
            chlshrr(mgs) = 0.0
          ENDIF
          
        ENDIF

        IF ( lf > 1 .and. qfshr(mgs) < 0.0 ) THEN
          IF ( -qfshr(mgs) <= qfacr(mgs) + qfacw(mgs) ) THEN
            frac = -qfshr(mgs)/( qfacr(mgs) + qfacw(mgs) )
            qfacr(mgs) = frac*qfacr(mgs)
            qfacw(mgs) = frac*qfacw(mgs)
            cfacr(mgs) = frac*cfacr(mgs)
            cfacw(mgs) = frac*cfacw(mgs)
            qfshr(mgs) = 0.0
            cfshrr(mgs) = 0.0 ! frac*cfshrr(mgs)
          ELSE
            qfshr(mgs) = 0.0 ! qfshr(mgs) + qfacr(mgs) + qfacw(mgs)
            qfacr(mgs) = 0.0
            qfacw(mgs) = 0.0
            cfacr(mgs) = 0.0
            cfacw(mgs) = 0.0
            cfshrr(mgs) = 0.0
          ENDIF
          
        ENDIF
       

        ! resum, but should all be zeros
        qrshr(mgs) = qsshr(mgs) + qhshr(mgs) + qhlshr(mgs)
        crshr(mgs) = chshrr(mgs)/rzxh(mgs) + chlshrr(mgs)/rzxhl(mgs)
      
        IF ( lf > 1 ) THEN
         qrshr(mgs) = qrshr(mgs) + qfshr(mgs)
         crshr(mgs) = crshr(mgs) + cfshrr(mgs)/rzxf(mgs)
        ENDIF
      
      ENDIF ! ished2cld

       ! find fraction of collected auto and melt that end up back as rain.
       ! Then pass these through so that they are not reprocessed as shed, then shed
       ! is only from the excess droplets:
       ! 
       ! qhdry
!  T < 0      qhshr(mgs)  = Min( 0.0, qhwet(mgs) - qhdry(mgs) )  ! water that freezes should never be more than what sheds
!  T >= 0     qhshr(mgs)  = - qhacw(mgs) - qhacr(mgs) ! -qhdry(mgs)
      IF (  iraintypes >= 1 .and. ished2cld >= 0 ) THEN
        IF ( qhshr(mgs) < -0.1*qxmin(lr) .and. qhacr(mgs) > 0.0 ) THEN
          frac = -qhshr(mgs)*qhacr(mgs)/( qhacr(mgs) + qhacw(mgs) ) ! fraction of collected rain that is shed
          qxshedfrac(mgs,lh,0) = frac ! qxrainold(mgs,1)/qx(mgs,lr)
        ENDIF

        IF ( lhl > 0 ) THEN
        IF ( qhlshr(mgs) < -0.1*qxmin(lr) .and. qhlacr(mgs) > 0.0 ) THEN
          frac = -qhlshr(mgs)*qhlacr(mgs)/( qhlacr(mgs) + qhlacw(mgs) ) ! fraction of collected rain that is shed
          qxshedfrac(mgs,lhl,0) = frac ! qxrainold(mgs,1)/qx(mgs,lr)
        ENDIF
        ENDIF

        IF ( lf > 0 ) THEN
        IF ( qfshr(mgs) < -0.1*qxmin(lr) .and. qfacr(mgs) > 0.0 ) THEN
          frac = -qfshr(mgs)*qfacr(mgs)/( qfacr(mgs) + qfacw(mgs) ) ! fraction of collected rain that is shed
          qxshedfrac(mgs,lf,0) = frac ! qxrainold(mgs,1)/qx(mgs,lr)
        ENDIF
        ENDIF
        
        IF ( qrshr(mgs) < -0.1*qxmin(lr) .and. qx(mgs,lr) > qxmin(lr) ) THEN
        frac = Sum( qxshedfrac(mgs,lh:lhab,0) ) ! total shed from collected rain
        qxshedfrac(mgs,lr,0) = Abs( frac/qrshr(mgs) )
        qxshedfrac(mgs,lr,1) = qxshedfrac(mgs,lr,0)*qxrain(mgs,1)/qx(mgs,lr) ! collected rain-auto
        qxshedfrac(mgs,lr,3) = qxshedfrac(mgs,lr,0)*qxrain(mgs,3)/qx(mgs,lr) ! collected rain-melt
        qxshedfrac(mgs,lr,2) = 1.0 - ( qxshedfrac(mgs,lr,1) + qxshedfrac(mgs,lr,3) ) ! collected droplets + rain-shed
!         IF ( ny <= 2 .and. frac > 0.0 ) THEN
!           write(0,*) 'qxshedfrac = ',igs(mgs),kgs(mgs),qxshedfrac(mgs,lr,0),qxshedfrac(mgs,lr,1),qxshedfrac(mgs,lr,2),qxshedfrac(mgs,lr,3)
!           write(0,*) 'qrshr,qxacr,qxacw = ',qrshr(mgs),(qhacr(mgs)+qfacr(mgs)+qhlacr(mgs)),(qhacw(mgs)+qfacw(mgs)+qhlacw(mgs))
!         ENDIF
        ELSE
          qxshedfrac(mgs,lr,0:nraintypes) = 0.0
          qxshedfrac(mgs,lr,2) = 1.0
        ENDIF
      ENDIF
        

      
      IF ( ipconc .ge. 3 ) THEN
!       crshr(mgs) = Max(crshr(mgs), rho0(mgs)*qrshr(mgs)/(xdn(mgs,lr)*vr1mm) )
      ENDIF
      end do 
!
!
!
      if (ndebug .gt. 0 ) write(0,*) 'dbg = 8a'

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
!       DO mgs = 1,ngscnt
       pccwi(:) = 0.0
       pccwd(:) = 0.0
       pccwdacc(:) = 0.0
       pccii(:) = 0.0
       pccin(:) = 0.0
       pccid(:) = 0.0
       pcisi(:) = 0.0
       pcisd(:) = 0.0
       pcrwi(:) = 0.0
       pcrwd(:) = 0.0
       pcswi(:) = 0.0
       pcswd(:) = 0.0
       pchwi(:) = 0.0
       pchwd(:) = 0.0
       pcfwi(:) = 0.0
       pcfwd(:) = 0.0
       pchli(:) = 0.0
       pchld(:) = 0.0
!       ENDDO
!
!  Cloud ice (columns)
!
!      IF ( ipconc .ge. 1 ) THEN
      if (ndebug .gt. 0 ) write(0,*) 'cloud ice sum'

      IF ( warmonly < 0.5 ) THEN
      IF ( ffrzs < 1.0 ) THEN
      do mgs = 1,ngscnt
      pccii(mgs) =   &
     &   il5(mgs)*cicint(mgs)  &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))   &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  +cfmul1(mgs)   &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs)
     
       pccii(mgs) = pccii(mgs)*(1.0 - ffrzs)
       
!     >  + nsplinter*(crfrzf(mgs) + crfrz(mgs))
      pccid(mgs) =   &
     &   il5(mgs)*(-cscni(mgs) - cscnvi(mgs) & ! - cwaci(mgs)   &
     &  -craci(mgs)    &
     &  -csaci(mgs)   &
     & - cfaci(mgs) &
     &  -chaci(mgs) - chlaci(mgs)   &
     &  -chcni(mgs))   &
     &  +il5(mgs)*cisbv(mgs)   &
     &  -(1.-il5(mgs))*cimlr(mgs)

      pccin(mgs) = ciint(mgs)
      

      end do
      ENDIF ! ffrzs
      ELSEIF ( warmonly < 0.8 ) THEN
      do mgs = 1,ngscnt
      
!      qiint(mgs) = 0.0
!      cicint(mgs) = 0.0
!      qicicnt(mgs) = 0.0
      
      pccii(mgs) =   &
     &   il5(mgs)*cicint(mgs)   &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))   &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  +cfmul1(mgs)   &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs)

       pccii(mgs) = pccii(mgs)*(1. - ffrzs)
      pccid(mgs) =   &
!     &   il5(mgs)*(-cscni(mgs) - cscnvi(mgs) & ! - cwaci(mgs)   &
!     &  -craci(mgs)    &
!     &  -csaci(mgs)   &
!     &  -chaci(mgs) - chlaci(mgs)   &
!     &  -chcni(mgs))   &
     &  +il5(mgs)*cisbv(mgs)   &
     &  -(1.-il5(mgs))*cimlr(mgs)

      pccin(mgs) = ciint(mgs)

      end do
      ENDIF ! warmonly

      
!      ENDIF ! ( ipconc .ge. 1 )
!
!  Cloud water
!
      IF ( ipconc .ge. 2 ) THEN
      
      do mgs = 1,ngscnt
      pccwi(mgs) =  (0.0) - cwshw(mgs) ! + (1-il5(mgs))*(-cirmlw(mgs))
      
      IF ( warmonly < 0.5 ) THEN
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(-ciacw(mgs)-cwfrz(mgs)-cwctfzp(mgs)   &
     &  -cwctfzc(mgs)   &
     &   )   &
     &  -cfacw(mgs) &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)


      ELSEIF ( warmonly < 0.8 ) THEN
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(  &
     & -ciacw(mgs)-cwfrz(mgs)-cwctfzp(mgs)   &
     &  -cwctfzc(mgs)   &
     &   )   &
     &  -cracw(mgs) -chacw(mgs) -chlacw(mgs) 
      ELSE
      
!       tmp3d(igs(mgs),jy,kgs(mgs)) = crcnw(mgs)

!       cracw(mgs) = 0.0 ! turn off accretion
!       qracw(mgs) = 0.0
!       crcev(mgs) = 0.0 ! turn off evap
!       qrcev(mgs) = 0.0 ! turn off evap
!       cracr(mgs) = 0.0 ! turn off self collection
       
       
!       cautn(mgs) = 0.0 
!       crcnw(mgs) = 0.0
!       qrcnw(mgs) = 0.0

      pccwd(mgs) =    &
     &  - cautn(mgs) -cracw(mgs)
      ENDIF


      IF ( .false. .and. exwmindiam > 0.0 .and. ccwresv(mgs) > 0.0 ) THEN
      pccwdacc(mgs) =    &
     &  il5(mgs)*(-ciacw(mgs)  &
     &   )   &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)

      IF ( -pccwdacc(mgs)*dtp .gt. cx(mgs,lc) - ccwresv(mgs) ) THEN

       frac = -(cx(mgs,lc) - ccwresv(mgs) )/(pccwdacc(mgs)*dtp)
       pccwdacc(mgs) = -(cx(mgs,lc) - ccwresv(mgs) )*dtpinv

        ciacw(mgs)   = frac*ciacw(mgs)
        cracw(mgs)   = frac*cracw(mgs)
        csacw(mgs)   = frac*csacw(mgs)
        chacw(mgs)   = frac*chacw(mgs)
        cautn(mgs)   = frac*cautn(mgs)
        cfacw(mgs)   = frac*cfacw(mgs)
       
        IF ( lhl .gt. 1 ) chlacw(mgs)   = frac*chlacw(mgs)

! resum
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(-ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)   &
     &  -cwfrzc(mgs)-cwctfzc(mgs)   &
     &  -il5(mgs)*(ciihr(mgs))   &
     &   )   &
     &  -cfacw(mgs) &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)

      ENDIF

      ENDIF


      IF ( -pccwd(mgs)*dtp .gt. cx(mgs,lc) ) THEN
!       write(0,*) 'OUCH! pccwd(mgs)*dtp .gt. ccw(mgs) ',pccwd(mgs),cx(mgs,lc)
!       write(0,*) 'qc = ',qx(mgs,lc)
!       write(0,*) -ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)-cwfrzc(mgs)-cwctfzc(mgs)
!       write(0,*)  -cracw(mgs) -csacw(mgs)  -chacw(mgs)
!       write(0,*) - cautn(mgs)

       frac = -cx(mgs,lc)/(pccwd(mgs)*dtp)
       pccwd(mgs) = -cx(mgs,lc)*dtpinv

        ciacw(mgs)   = frac*ciacw(mgs)
        cwfrz(mgs)  = frac*cwfrz(mgs)
        cwfrzp(mgs)  = frac*cwfrzp(mgs)
        cwctfzp(mgs) = frac*cwctfzp(mgs)
        cwfrzc(mgs)  = frac*cwfrzc(mgs)
        cwctfzc(mgs) = frac*cwctfzc(mgs)
        cwctfz(mgs) = frac*cwctfz(mgs)
        cracw(mgs)   = frac*cracw(mgs)
        csacw(mgs)   = frac*csacw(mgs)
        chacw(mgs)   = frac*chacw(mgs)
        cautn(mgs)   = frac*cautn(mgs)
        cfacw(mgs)   = frac*cfacw(mgs)
       
        pccii(mgs) = pccii(mgs) - (1.-frac)*il5(mgs)*(cwfrzc(mgs)+cwctfzc(mgs))*(1. - ffrzs)
        IF ( lhl .gt. 1 ) chlacw(mgs)   = frac*chlacw(mgs)


!       STOP
      ENDIF

      end do

      ENDIF ! ipconc

!
!  Rain
!
      IF ( ipconc .ge. 3 ) THEN

      do mgs = 1,ngscnt

      IF ( warmonly < 0.5 ) THEN
      pcrwi(mgs) = &
!     >   cracw(mgs) +    &
     &   crcnw(mgs)   &
     &  +(1-il5(mgs))*(   &
     &    -cfmlrr(mgs)/rzxf(mgs)   &
     &    -chmlrr(mgs)/rzxh(mgs)   &
     &    -chlmlrr(mgs)/rzxhl(mgs)   &
!     &    -csmlr(mgs)/rzxs(mgs)     &
     &   -csmlrr(mgs)     &
     &   - cimlr(mgs) )   &
     &  - Min(0.0,cracr(mgs)) &  ! cracr is negative if there is enough breakup
     &  -crshr(mgs)             !null at this point when wet snow/graupel included
      pcrwd(mgs) =   &
     &   il5(mgs)*(-ciacr(mgs) - crfrz(mgs) ) & ! - cipacr(mgs))
!     >  -csacr(mgs)   &
     & - cfacr(mgs) &
     &  - chacr(mgs) - chlacr(mgs)   &
     &  +crcev(mgs)   &
     &  - Max(0.0,cracr(mgs))
!     >  -il5(mgs)*ciracr(mgs)

      ELSEIF ( warmonly < 0.8 ) THEN
       pcrwi(mgs) = &
     &   crcnw(mgs)   &
     &  +(1-il5(mgs))*(   &
     &    -chmlrr(mgs)/rzxh(mgs)    &
     &    -chlmlrr(mgs)/rzxhl(mgs)   &
!     &    -csmlr(mgs)     &
     &   -csmlrr(mgs)     &
     &   - cimlr(mgs) )   &
     &  -crshr(mgs)             !null at this point when wet snow/graupel included
      pcrwd(mgs) =   &
     &   il5(mgs)*( - crfrz(mgs) ) & ! - cipacr(mgs))
     &  - chacr(mgs)    &
     &  - chlacr(mgs)    &
     &  +crcev(mgs)   &
     &  - cracr(mgs)
      ELSE
      pcrwi(mgs) =   &
     &   crcnw(mgs)
      pcrwd(mgs) =   &
     &  +crcev(mgs)   &
     &  - cracr(mgs)

!        tmp3d(igs(mgs),jy,kgs(mgs)) = vtxbar(mgs,lr,1) ! crcnw(mgs) ! (pcrwi(mgs) + pcrwd(mgs))
!        pcrwi(mgs) = 0.0
!        pcrwd(mgs) = 0.0
!        qrcnw(mgs) = 0.0

      ENDIF


      frac = 0.0
      IF ( -pcrwd(mgs)*dtp .gt. cx(mgs,lr) ) THEN
!       write(0,*) 'OUCH! pcrwd(mgs)*dtp .gt. crw(mgs) ',pcrwd(mgs)*dtp,cx(mgs,lr),mgs,igs(mgs),kgs(mgs)
!       write(0,*) -ciacr(mgs)
!       write(0,*) -crfrz(mgs)
!       write(0,*) -chacr(mgs)
!       write(0,*)  crcev(mgs)
!       write(0,*)  -cracr(mgs)

       frac =  -cx(mgs,lr)/(pcrwd(mgs)*dtp)
       pcrwd(mgs) = -cx(mgs,lr)*dtpinv

        ciacr(mgs) = frac*ciacr(mgs)
        ciacrf(mgs) = frac*ciacrf(mgs)
        ciacrs(mgs) = frac*ciacrs(mgs)
        crfrz(mgs) = frac*crfrz(mgs)
        crfrzf(mgs) = frac*crfrzf(mgs)
        crfrzs(mgs) = frac*crfrzs(mgs)
        chacr(mgs) = frac*chacr(mgs)
        chlacr(mgs) = frac*chlacr(mgs)
        crcev(mgs) = frac*crcev(mgs)
        cracr(mgs) = frac*cracr(mgs)
        cfacr(mgs) = frac*cfacr(mgs)

!       STOP
      ENDIF

      end do

      ENDIF


      IF ( warmonly < 0.5 ) THEN

!
!  Snow
!
      IF ( ipconc .ge. 4 ) THEN !

      do mgs = 1,ngscnt
      pcswi(mgs) =   &
     &   il5(mgs)*(cscnis(mgs) + cscnvis(mgs) )    &
     &  + cwfrz2snowfrac*cwfrz(mgs)/cwfrz2snowratio  &
     &  + cscnh(mgs)
      
      IF (  ffrzs > 0.0 ) THEN
       pcswi(mgs) =  pcswi(mgs) + ffrzs* (  &
     &   il5(mgs)*cicint(mgs)   &
     &  +il5(mgs)*(cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))  &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  +cfmul1(mgs)   &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs) )
      ENDIF

      
      IF ( ess0 < 0.0 ) THEN
         csacs(mgs) = Max(0.0, csacs(mgs) - (ifrzs)*(crfrzs(mgs) + ciacrs(mgs)))
      ENDIF
      
      pcswd(mgs) = &
!     :  cracs(mgs)     &
     & - cfacs(mgs) &
     &  -chacs(mgs) - chlacs(mgs)   &
     &  -chcns(mgs)   &
     &  +(1-il5(mgs))*csmlr(mgs) + csshr(mgs) & ! + csshrp(mgs)
!     >  +il5(mgs)*(cssbv(mgs))   &
     &   + cssbv(mgs)   &
     &  - csacs(mgs)

      frac = 0.0
      IF ( imixedphase == 0 ) THEN
        IF ( cx(mgs,ls) + dtp*(pcswi(mgs) + pcswd(mgs)) < 0.0 ) THEN
         frac = (-cx(mgs,ls) + pcswi(mgs)*dtp)/(pcswd(mgs)*dtp)
         
           pcswd(mgs) = frac*pcswd(mgs)
           
           chacs(mgs)  = frac*chacs(mgs) 
           chlacs(mgs) = frac*chlacs(mgs)
           chcns(mgs)  = frac*chcns(mgs) 
           csmlr(mgs)  = frac*csmlr(mgs) 
           csshr(mgs)  = frac*csshr(mgs) 
           cssbv(mgs)  = frac*cssbv(mgs) 
           csacs(mgs)  = frac*csacs(mgs)
           cfacs(mgs)  = frac*cfacs(mgs) 
      
        ENDIF
      ENDIF


      
      pccii(mgs) =  pccii(mgs) &
     &  + (1. - ifrzs)*crfrzs(mgs) &
     &  + (1. - ifrzs)*ciacrs(mgs)

      pcswi(mgs) =  pcswi(mgs) &
     &  + (ifrzs)*crfrzs(mgs) &
     &  + (ifrzs)*ciacrs(mgs)

      end do

      ENDIF

!
!  Graupel
!
      IF ( ipconc .ge. 5 ) THEN !
      do mgs = 1,ngscnt
      pchwi(mgs) =   &
     &  +(ffrzh*ifrzg*crfrzf(mgs)   &
     & +il5(mgs)*ffrzh*ifiacrg*(ciacrf(mgs) ))    &
     & + f2h*chcnsh(mgs) + f2h*chcnih(mgs) + chcnhl(mgs)

      pchwd(mgs) =   &
     &  (1-il5(mgs))*chmlr(mgs) &
!     >  + il5(mgs)*chsbv(mgs)   &
     &  + chsbv(mgs)   &
     &  - il5(mgs)*chlcnh(mgs) &
     &  - cscnh(mgs)

        IF ( lnhf > 1 ) THEN
           IF ( cx(mgs,lh) > 0.0 ) THEN
             frac = Min(1.0, chxf(mgs,lh)/cx(mgs,lh))
             frac = Max(0.0, frac)
           ELSE
             frac = 0.0
           ENDIF
 
           tmp = ffrzh*ifrzg*crfrzf(mgs) +il5(mgs)*ffrzh*ifiacrg*(ciacrf(mgs) )

           chxf(mgs,lh) = chxf(mgs,lh) + dtp*( tmp + frac*((1-il5(mgs))*chmlr(mgs) + chsbv(mgs) - il5(mgs)*chlcnh(mgs)) )
           
!           IF ( pchwi(mgs) /= 0.0 .or. pchwd(mgs) /= 0.0 ) THEN
!             write(0,*) 'i,k,tmp,delta-chxf = ',igs(mgs),kgs(mgs),tmp,dtp*( tmp + frac*((1-il5(mgs))*chmlr(mgs) + chsbv(mgs) - il5(mgs)*chlcnh(mgs)) )
!           ENDIF

        ENDIF
      end do


!
!  Frozen drops
!
      IF ( ipconc .ge. 5 .and. lf > 1 .and. lnf > 1 ) THEN !
      do mgs = 1,ngscnt
      pcfwi(mgs) =   &
     &  +((1.0-ffrzh)*crfrzf(mgs)   &
     & +il5(mgs)*(1.0-ffrzh)*(ciacrf(mgs) )) + (1.0-f2h)*chcnsh(mgs) + (1.0-f2h)*chcnih(mgs)


      pcfwd(mgs) =   &
     &  (1-il5(mgs))*cfmlr(mgs) &
     &  + cfsbv(mgs)   &
     &  - il5(mgs)*chlcnf(mgs)
      
      IF ( .false. .and. kgs(mgs) <= 20 .and. ( cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ) > 200. .or. cx(mgs,lf) > 100. )) THEN
        write(0,*) 'ix,jy, kz, cf = ',igs(mgs)+ixbeg,jy+jybeg,kgs(mgs),cx(mgs,lf)
        write(0,*) 'cf_new,pcfwi,pcfwd = ',cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ),pcfwi(mgs) + pcfwd(mgs)
      
      ENDIF
      
!      IF ( cx(mgs,lf) > cxmin .and. cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ) < 0.0 .and.  &
!           qx(mgs,lf) + dtp*( pqfwi(mgs) + pqfwd(mgs) ) > qxmin(lf) ) THEN
!      ! overdepletion
!        ! rescale depletion
!
!         frac = (-cx(mgs,lf) + pcfwi(mgs)*dtp)/(pcfwd(mgs)*dtp)
!         
!         cfmlr(mgs) = frac*cfmlr(mgs)
!         cfsbv(mgs) = frac*cfsbv(mgs)
!         chcnf(mgs) = frac*chcnf(mgs)
!           
!         pcfd(mgs) = frac*pcfd(mgs)
!      
!      
!      ENDIF
      
      end do
      ENDIF


!

!
!  Hail
!
      IF ( lhl .gt. 1 .and. lnhl > 1 ) THEN !
      do mgs = 1,ngscnt
      pchli(mgs) = (ffrzh*(1.0-ifrzg)*crfrzf(mgs) +il5(mgs)*ffrzh*(1.0-ifiacrg)*(ciacrf(mgs) ))  &
     &  + chlcnfhl(mgs)*rzxhlf(mgs) &
     & + chlcnhhl(mgs) *rzxhlh(mgs)

      pchld(mgs) =   &
     &  (1-il5(mgs))*chlmlr(mgs)   &
!     >  + il5(mgs)*chlsbv(mgs)   &
     &  + chlsbv(mgs) - chcnhl(mgs)
      
      IF ( imixedphase == 0 ) THEN
      frac = 0.0
      IF ( cx(mgs,lhl) + dtp*(pchli(mgs) + pchld(mgs)) < 0.0 ) THEN
        ! rescale depletion

         frac = (-cx(mgs,lhl) + pchli(mgs)*dtp)/(pchld(mgs)*dtp)
         
         chlmlr(mgs) = frac*chlmlr(mgs)
         chlsbv(mgs) = frac*chlsbv(mgs)
         chcnhl(mgs) = frac*chcnhl(mgs)
           
         pchld(mgs) = frac*pchld(mgs)
           
      ENDIF
      ENDIF

        IF ( lnhlf > 1 ) THEN
           IF ( cx(mgs,lhl) > 0.0 ) THEN
             frac = Min(1.0, chxf(mgs,lhl)/cx(mgs,lhl))
             frac = Max(0.0, frac)
           ELSE
             frac = 0.0
           ENDIF
 
           IF ( lf > 1 ) THEN ! number of hail from frozen drops
             IF ( ifddenfac >= 1 ) THEN
                IF ( xdn(mgs,lf) <= fddenthresh ) THEN 
                  fddenfac(mgs) = 0.0  ! could make this a "smooth" transition with tanh?
                ELSE
                  fddenfac(mgs) = 1.0
                ENDIF
             ELSE
                fddenfac(mgs) = 1.0
             ENDIF
             chxf(mgs,lhl) = chxf(mgs,lhl) + dtp*( chlcnfhl(mgs)*fddenfac(mgs)*rzxhlf(mgs) + frac*(1-il5(mgs))*chlmlr(mgs) + frac*chlsbv(mgs) )
             chlfmlr(mgs) = frac*(1-il5(mgs))*chlmlr(mgs)
           ELSEIF ( lnhf > 1 ) THEN ! transfer fraction of lnhf
             IF ( cx(mgs,lh) > 0.0 ) THEN
               frach = Min(1.0, chxf(mgs,lh)/cx(mgs,lh))
               frach = Max(0.0, frach)
             ELSE
               frach = 0.0
             ENDIF
             
             chxf(mgs,lhl) = chxf(mgs,lhl) + dtp*( frach*chlcnhhl(mgs)*rzxhlh(mgs) + frac*(1-il5(mgs))*chlmlr(mgs) + frac*chlsbv(mgs) )
        
          ENDIF
        ENDIF
      end do
      
      ENDIF
!

      ENDIF ! (ipconc .ge. 5 )

      ELSEIF ( warmonly < 0.8 ) THEN

!
!  Graupel
!
      IF ( ipconc .ge. 5 ) THEN !
      do mgs = 1,ngscnt
      pchwi(mgs) =   &
     &  +ifrzg*(crfrzf(mgs) ) ! +il5(mgs)*(ciacrf(mgs) ))

      pchwd(mgs) =   &
     &  (1-il5(mgs))*chmlr(mgs) &
     &  - il5(mgs)*chlcnh(mgs)
      end do
!
!  Hail
!
      IF ( lhl .gt. 1 ) THEN !
      do mgs = 1,ngscnt
      pchli(mgs) = (1.0-ifrzg)*(crfrzf(mgs)) & ! +il5(mgs)*(ciacrf(mgs) ))  &
     &  + chlcnfhl(mgs) &
     & + chlcnhhl(mgs) *rzxhl(mgs)/rzxh(mgs)

      pchld(mgs) =   &
     &  (1-il5(mgs))*chlmlr(mgs) !  &
!     >  + il5(mgs)*chlsbv(mgs)   &
!     &  + chlsbv(mgs)

!      IF ( pchli(mgs) .ne. 0. .or. pchld(mgs) .ne. 0 ) THEN
!       write(0,*) 'dr: pchli,pchld = ', pchli(mgs),pchld(mgs), igs(mgs),kgs(mgs)
!      ENDIF
      end do

      ENDIF

      ENDIF ! ipconc >= 5

      ENDIF ! warmonly

!

!
!  Balance and checks for continuity.....within machine precision...
!
      do mgs = 1,ngscnt
      pctot(mgs)   = pccwi(mgs) +pccwd(mgs) +   &
     &               pccii(mgs) +pccid(mgs) +   &
     &               pcrwi(mgs) +pcrwd(mgs) +   &
     &               pcswi(mgs) +pcswd(mgs) +   &
     &               pchwi(mgs) +pchwd(mgs) +   &
     &               pcfwi(mgs) +pcfwd(mgs) +   &
     &               pchli(mgs) +pchld(mgs)
      end do
!
!
      ENDIF ! ( ipconc .ge. 1 )
!
!
!
!
!
!  GOGO
!  production terms for mass
!
!
       pqwvi(:) = 0.0
       pqwvd(:) = 0.0
       pqcwi(:) = 0.0
       pqcwd(:) = 0.0
       pqcwdacc(:) = 0.0
       pqcii(:) = 0.0
       pqcid(:) = 0.0
       pqrwi(:) = 0.0
       pqrwd(:) = 0.0
       pqswi(:) = 0.0
       pqswd(:) = 0.0
       pqhwi(:) = 0.0
       pqhwd(:) = 0.0
       pqhli(:) = 0.0
       pqhld(:) = 0.0
       pqlwsi(:) = 0.0
       pqlwsd(:) = 0.0
       pqlwhi(:) = 0.0
       pqlwhd(:) = 0.0
       pqlwlghi(:) = 0.0
       pqlwlghd(:) = 0.0
       pqlwlghli(:) = 0.0
       pqlwlghld(:) = 0.0
       pqlwhli(:) = 0.0
       pqlwhld(:) = 0.0
!       IF ( lf > 1 ) THEN
       pqfwi(:) = 0.0
       pqfwd(:) = 0.0
       pqlwfi(:) = 0.0
       pqlwfd(:) = 0.0
!       ENDIF
       pqrauto(:) = 0.0
       pqrshed(:) = 0.0
       pqrmelt(:) = 0.0
       pqrother(:) = 0.0
       IF ( ipconc > 5 ) THEN
       pzhwi(:) = 0.0
       pzhwd(:) = 0.0
       pzrwi(:) = 0.0
       pzrwd(:) = 0.0
       pzhli(:) = 0.0
       pzhld(:) = 0.0
       ENDIF


!
!  Vapor
!
      IF ( warmonly < 0.5 ) THEN
      do mgs = 1,ngscnt
      
! NOTE: ANY CHANGES HERE ALSO NEED TO GO INTO THE RESUM FARTHER DOWN!
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))   &
     &  -Min(0.0, qhcev(mgs))   &
     &  -Min(0.0, qfcev(mgs))   &
     &  -Min(0.0, qhlcev(mgs))   &
     &  -Min(0.0, qscev(mgs))   &
!     >  +il5(mgs)*(-qhsbv(mgs) - qhlsbv(mgs) )   &
     &  -qhsbv(mgs) - qhlsbv(mgs)   &
     &  -qssbv(mgs)    &
     &  -qfsbv(mgs)    &
     &  -il5(mgs)*qisbv(mgs)
      
      pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))   &
     &  -Max(0.0, qhcev(mgs))   &
     &  -Max(0.0, qfcev(mgs))   &
     &  -Max(0.0, qhlcev(mgs))   &
     &  -Max(0.0, qscev(mgs))   &
     &  +il5(mgs)*(-qiint(mgs)   &
     &  -qfdpv(mgs)    &
     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -il5(mgs)*qidpv(mgs)  
      
      end do

      ELSEIF ( warmonly < 0.8 ) THEN
      do mgs = 1,ngscnt
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs)) &
     &  -il5(mgs)*qisbv(mgs)
      pqwvd(mgs) =     &
     &  +il5(mgs)*(-qiint(mgs)   &
!     &  -qhdpv(mgs) ) & !- qhldpv(mgs))   &
     &  -qfsbv(mgs)    &
     &  -qhdpv(mgs) - qhldpv(mgs))   &
!     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -Max(0.0, qrcev(mgs))     &
     &  -il5(mgs)*qidpv(mgs)  
      end do

      ELSE
      do mgs = 1,ngscnt
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))
      pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))
      end do

      ENDIF ! warmonly
!
!  Cloud water
!
      do mgs = 1,ngscnt

      pqcwi(mgs) =  (0.0) + qwcnr(mgs) - qwshw(mgs)

      IF ( warmonly < 0.5 ) THEN
      pqcwd(mgs) =    &
     &  il5(mgs)*(-qiacw(mgs)-qwfrz(mgs)-qwctfz(mgs))   &
     &  -il5(mgs)*(qiihr(mgs))   &
     &  - qfacw(mgs) &
     &  -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs) - qhlacw(mgs)  !&
!     &  -il5(mgs)*(qwfrzp(mgs))
      ELSEIF ( warmonly < 0.8 ) THEN
      pqcwd(mgs) =    &
     &  il5(mgs)*(-qiacw(mgs)-qwfrz(mgs)-qwctfz(mgs))   &
     &  -il5(mgs)*(qiihr(mgs))   &
     &  -qracw(mgs) -qrcnw(mgs) -qhacw(mgs) -qhlacw(mgs)
      ELSE
      pqcwd(mgs) =    &
     &  -qracw(mgs) - qrcnw(mgs)
      ENDIF


      IF ( pqcwd(mgs) .lt. 0.0 .and. -pqcwd(mgs)*dtp .gt. qx(mgs,lc) ) THEN

       frac = -Max(0.0,qx(mgs,lc))/(pqcwd(mgs)*dtp)
       pqcwd(mgs) = -qx(mgs,lc)*dtpinv

        qiacw(mgs)   = frac*qiacw(mgs)
!        qwfrzp(mgs)  = frac*qwfrzp(mgs)
!        qwctfzp(mgs) = frac*qwctfzp(mgs)
        qwfrzc(mgs)  = frac*qwfrzc(mgs)
        qwfrz(mgs)  = frac*qwfrz(mgs)
        qwctfzc(mgs) = frac*qwctfzc(mgs)
        qwctfz(mgs) = frac*qwctfz(mgs)
        qracw(mgs)   = frac*qracw(mgs)
        qsacw(mgs)   = frac*qsacw(mgs)
        qhacw(mgs)   = frac*qhacw(mgs)
        qfacw(mgs)   = frac*qfacw(mgs)
        vhacw(mgs)   = frac*vhacw(mgs)
        qrcnw(mgs)   = frac*qrcnw(mgs)
        qwfrzp(mgs)  = frac*qwfrzp(mgs)
        IF ( lhl .gt. 1 ) THEN
          qhlacw(mgs)   = frac*qhlacw(mgs)
          vhlacw(mgs)   = frac*vhlacw(mgs)
        ENDIF
!        IF ( lzh .gt. 1 ) zhacw(mgs) = frac*zhacw(mgs)

!       STOP
      ENDIF
      

      end do
!
!  Cloud ice
!
      IF ( warmonly < 0.5 ) THEN

      do mgs = 1,ngscnt
      IF ( ffrzs < 1.0 ) THEN
      pqcii(mgs) =     &
     &   il5(mgs)*qicicnt(mgs)    &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*qwfrzc(mgs)+qwctfzc(mgs))   &
     &  +il5(mgs)*(qicichr(mgs))  &
     &  +qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     &  + qfmul1(mgs) &
     & + qsplinter(mgs) + qsplinter2(mgs)
!     > + cimas0*nsplinter*(crfrzf(mgs) + crfrz(mgs))/rho0(mgs)
      ENDIF
       
       pqcii(mgs) = pqcii(mgs)*(1.0 - ffrzs) &
     &  +il5(mgs)*qidpv(mgs)    &
     &  +il5(mgs)*qiacw(mgs)
       
      pqcid(mgs) =     &
     &   il5(mgs)*(-qscni(mgs) - qscnvi(mgs)    & ! -qwaci(mgs)    &
     &  -qraci(mgs)    &
     &  -qsaci(mgs) )   &
     & - qfaci(mgs) &
     &  -qhaci(mgs)   &
     &  -qhlaci(mgs)    &
     &  +il5(mgs)*qisbv(mgs)    &
     &  +(1.-il5(mgs))*qimlr(mgs)   &
     &  - qhcni(mgs)
      end do

      
      ELSEIF ( warmonly < 0.8 ) THEN

      do mgs = 1,ngscnt
      pqcii(mgs) =     &
     &   il5(mgs)*qicicnt(mgs)*(1. - ffrzs)    &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*qwfrzc(mgs)+qwctfzc(mgs))*(1. - ffrzs)   &
     &  +il5(mgs)*(qicichr(mgs))*(1. - ffrzs)   &
!     &  +il5(mgs)*(qicichr(mgs))   &
!     &  +qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     & + qsplinter(mgs) + qsplinter2(mgs) &
     &  +il5(mgs)*qidpv(mgs)    &
     &  +il5(mgs)*qiacw(mgs)  ! & ! (qiacwi(mgs)+qwacii(mgs))   &
!     &  +il5(mgs)*(qwfrzc(mgs)+qwctfzc(mgs))   &
!     &  +il5(mgs)*(qicichr(mgs))   &
!     &  +qsmul(mgs)               &
!     &  +qhmul1(mgs) + qhlmul1(mgs)   &
!     & + qsplinter(mgs) + qsplinter2(mgs)

      pqcid(mgs) =     &
!     &   il5(mgs)*(-qscni(mgs) - qscnvi(mgs)    & ! -qwaci(mgs)    &
!     &  -qraci(mgs)    &
!     &  -qsaci(mgs) )   &
!     &  -qhaci(mgs)   &
!     &  -qhlaci(mgs)    &
     &  +il5(mgs)*qisbv(mgs)    &
     &  +(1.-il5(mgs))*qimlr(mgs)  ! &
!     &  - qhcni(mgs)
      end do

      ENDIF
!
!  Rain
!

      do mgs = 1,ngscnt
      IF ( warmonly < 0.5 ) THEN
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))   &
     &  +(1-il5(mgs))*(   &
     &    -qfmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qhmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qsmlr(mgs)  - qhlmlr(mgs)     &
     &    -qimlr(mgs))   &
!     &    -qsshr(mgs)       &                      !null at this point when wet snow/graupel included
!     &    -qfshr(mgs)                 &           !null at this point when wet snow/graupel included
!     &    -qhshr(mgs)       &                      !null at this point when wet snow/graupel included
!     &    -qhlshr(mgs)      &
     & - qrshr(mgs)

      pqrwd(mgs) =     &
     &  il5(mgs)*(-qiacr(mgs)-qrfrz(mgs))    &
     &  - qfacr(mgs) &
     &  - qsacr(mgs) - qhacr(mgs) - qhlacr(mgs) - qwcnr(mgs)   &
     &  + Min(0.0,qrcev(mgs))
      ELSEIF ( warmonly < 0.8 ) THEN
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))   &
     &  +(1-il5(mgs))*(   &
     &    -qfmlr(mgs)                 &           !null at this point when wet snow/graupel included
     &    -qhlmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qhmlr(mgs)  )               &            !null at this point when wet snow/graupel included
     &    -qhshr(mgs)                 &           !null at this point when wet snow/graupel included
     &    -qfshr(mgs)                 &           !null at this point when wet snow/graupel included
     &    -qhlshr(mgs)                            !null at this point when wet snow/graupel included
      pqrwd(mgs) =     &
     &  il5(mgs)*(-qrfrz(mgs))    &
     &   - qhacr(mgs)    &
     &   - qhlacr(mgs)    &
     &  + Min(0.0,qrcev(mgs))
      ELSE
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))
      pqrwd(mgs) =  Min(0.0,qrcev(mgs))
      ENDIF ! warmonly


 !      IF ( pqrwd(mgs) .lt. 0.0 .and. -(pqrwd(mgs) + pqrwi(mgs))*dtp .gt. qx(mgs,lr)  ) THEN
      IF ( pqrwd(mgs) .lt. 0.0 .and. -(pqrwd(mgs) + pqrwi(mgs))*dtp .gt. qx(mgs,lr)  ) THEN

       frac = (-qx(mgs,lr) + pqrwi(mgs)*dtp)/(pqrwd(mgs)*dtp)
!       pqrwd(mgs) = -qx(mgs,lr)*dtpinv  + pqrwi(mgs)

       pqwvi(mgs) = pqwvi(mgs)    &
     &  + Min(0.0, qrcev(mgs))   &
     &  - frac*Min(0.0, qrcev(mgs))
       pqwvd(mgs) =  pqwvd(mgs)   &
     &  + Max(0.0, qrcev(mgs))   &
     &  - frac*Max(0.0, qrcev(mgs))

       qiacr(mgs)  = frac*qiacr(mgs)
       qiacrf(mgs) = frac*qiacrf(mgs)
       qiacrs(mgs) = frac*qiacrs(mgs)
       viacrf(mgs) = frac*viacrf(mgs)
       qrfrz(mgs)  = frac*qrfrz(mgs) 
       qrfrzs(mgs) = frac*qrfrzs(mgs) 
       qrfrzf(mgs) = frac*qrfrzf(mgs)
       vrfrzf(mgs) = frac*vrfrzf(mgs)
       qsacr(mgs)  = frac*qsacr(mgs)
       qhacr(mgs)  = frac*qhacr(mgs)
       vhacr(mgs)  = frac*vhacr(mgs)
       qrcev(mgs)  = frac*qrcev(mgs)
       qhlacr(mgs) = frac*qhlacr(mgs)
       vhlacr(mgs) = frac*vhlacr(mgs)
       qfacr(mgs)  = frac*qfacr(mgs)
       vfacr(mgs)  = frac*vfacr(mgs)
       qfcev(mgs)  = frac*qfcev(mgs)
       qhcev(mgs)  = frac*qhcev(mgs)
       qhlcev(mgs)  = frac*qhlcev(mgs)


      IF ( warmonly < 0.5 ) THEN
       pqrwd(mgs) =     &
     &  il5(mgs)*(-qiacr(mgs)-qrfrz(mgs) - qsacr(mgs))    &
     &  - qhacr(mgs) - qhlacr(mgs) - qwcnr(mgs)   &
     &   - qfacr(mgs)    &
     &  + Min(0.0,qrcev(mgs))
      ELSEIF ( warmonly < 0.8 ) THEN
      pqrwd(mgs) =     &
     &  il5(mgs)*(-qrfrz(mgs))    &
     &   - qfacr(mgs)    &
     &   - qhacr(mgs)    &
     &   - qhlacr(mgs)    &
     &  + Min(0.0,qrcev(mgs))
      ELSE
       pqrwd(mgs) =  Min(0.0,qrcev(mgs))
      ENDIF ! warmonly

!
! Resum for vapor since qrcev has changed
!
      IF ( qrcev(mgs) .ne. 0.0 ) THEN
       pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))   &
     &  -Min(0.0, qhcev(mgs))   &
     &  -Min(0.0, qhlcev(mgs))   &
     &  -Min(0.0, qscev(mgs))   &
!     >  +il5(mgs)*(-qhsbv(mgs)  - qhlsbv(mgs) )   &
     &  -Min(0.0, qfcev(mgs))   &
     &  -qfsbv(mgs)    &
     &  -qhsbv(mgs)  - qhlsbv(mgs)   &
     &  -qssbv(mgs)    &
     &  -il5(mgs)*qisbv(mgs) 
     
       pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))   &
     &  -Max(0.0, qhcev(mgs))   &
     &  -Max(0.0, qhlcev(mgs))   &
     &  -Max(0.0, qscev(mgs))   &
     &  +il5(mgs)*(-qiint(mgs)   &
     &  -Max(0.0, qfcev(mgs))   &
     &  -qfdpv(mgs)    &
     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -il5(mgs)*qidpv(mgs)  

       ENDIF


!       STOP
      ENDIF

       ! rain rates that apply to all source types (collecting qc, evaporation, etc.)
       pqrother(mgs) =     &
     &   qracw(mgs)  + qrcev(mgs)   &
     &  + il5(mgs)*(-qiacr(mgs)-qrfrz(mgs))    &
     &  - qfacr(mgs) &
     &  - qsacr(mgs) - qhacr(mgs) - qhlacr(mgs) - qwcnr(mgs) 


       pqrauto(mgs) = qrcnw(mgs)

       IF ( .not. mixedphase ) THEN
        pqrmelt(mgs) =  &
     &  +(1-il5(mgs))*(   &
     &    -qhmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qfmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qsmlr(mgs)  - qhlmlr(mgs)     &
     &    -qimlr(mgs))

        pqrshed(mgs) = - qrshr(mgs)
!        pqrshed(mgs) = -( qsshr(mgs) + qhshr(mgs) + qhlshr(mgs) ) !  - qrshr(mgs)
!        IF ( lf > 0 ) THEN
!         pqrshed(mgs) = pqrshed(mgs) - qfshr(mgs)*qxshedfrac(mgs,lf,2)
!        ENDIF
      
       ENDIF


      end do

      IF ( warmonly < 0.5 ) THEN

!
!  Snow
!
      do mgs = 1,ngscnt
      pqswi(mgs) =     &
     &   il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs)   &
     &   + qscnvi(mgs)                        &
     &   + ifrzs*(qiacrs(mgs) + qrfrzs(mgs))  &
     &   + il5(mgs)*(( qwfrzc(mgs) + qwctfzc(mgs) + qicichr(mgs) )*ffrzs   &
     &   +  (1.0 - ffrzs)*cwfrz2snowfrac*qwfrz(mgs) ) &
     &   + il2(mgs)*qsacr(mgs))   &
     &   + il5(mgs)*qicicnt(mgs)*ffrzs        &
     &   + il3(mgs)*(qiacrf(mgs)+qracif(mgs)) & ! only applies for ipconc <= 3
     &   + Max(0.0, qscev(mgs))   &
     &   + qsacw(mgs) + qscnh(mgs) &
     &  + ffrzs*(qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     & + qsplinter(mgs) + qsplinter2(mgs))
      pqswd(mgs) =    &
!     >  -qfacs(mgs) ! -qwacs(mgs)   &
     &  -qracs(mgs)*(1-il2(mgs)) -qhacs(mgs) - qhlacs(mgs)   &
     & - qfacs(mgs) &
     &  -qhcns(mgs)   &
     &  +(1-il5(mgs))*qsmlr(mgs) + qsshr(mgs)    &    !null at this point when wet snow included
!     >  +il5(mgs)*(qssbv(mgs))   &
     &  + qssbv(mgs)   &
     &  + Min(0.0, qscev(mgs))  &
     &  -qsmul(mgs)
      
      
      IF ( imixedphase == 0 .and. pqswd(mgs) .lt. 0.0  ) THEN
        IF ( qx(mgs,ls) + dtp*(pqswi(mgs) + pqswd(mgs)) < 0.0 ) THEN
         frac = (-qx(mgs,ls) + pqswi(mgs)*dtp)/(pqswd(mgs)*dtp)
         
           pqswd(mgs) = frac*pqswd(mgs)
           
           qracs(mgs)  = frac*qracs(mgs) ! only used for single moment at this time
           qhacs(mgs)  = frac*qhacs(mgs) 
           qhlacs(mgs) = frac*qhlacs(mgs)
           qhcns(mgs)  = frac*qhcns(mgs) 
           qsmlr(mgs)  = frac*qsmlr(mgs) 
           qsshr(mgs)  = frac*qsshr(mgs) 
           qssbv(mgs)  = frac*qssbv(mgs) 
           qsmul(mgs)  = frac*qsmul(mgs) 
           qfacs(mgs)  = frac*qfacs(mgs) 
           IF ( qscev(mgs) < 0.0 ) qscev(mgs) = frac*qscev(mgs)

        ENDIF
      ENDIF
      
      pqcii(mgs) =  pqcii(mgs) &
     &  + (1. - ifrzs)*qrfrzs(mgs) &
     &  + (1. - ifrzs)*qiacrs(mgs)
      
      end do 
      
!
!  Graupel
!
      do mgs = 1,ngscnt
      pqhwi(mgs) =    &
     &  +il5(mgs)*(ffrzh*ifrzg*qrfrzf(mgs)  + (1-il3(mgs))*ffrzh*ifiacrg*(qiacrf(mgs)+qracif(mgs)))   &
     &  + (1-il2(mgs))*(qracs(mgs) + qsacr(mgs))  & ! only used for ipconc < 3
     &  +il5(mgs)*(qhdpv(mgs))   &
     &  +Max(0.0, qhcev(mgs))   &
     &  +qhacr(mgs)+qhacw(mgs)   &
     &  +qhacs(mgs)+qhaci(mgs)   &
     &  + f2h*qhcns(mgs) + f2h*qhcni(mgs) + qhcnhl(mgs)
      pqhwd(mgs) =     &
     &   qhshr(mgs)                &    !null at this point when wet graupel included
     &  +(1-il5(mgs))*qhmlr(mgs)   &    !null at this point when wet graupel included
!     >  +il5(mgs)*qhsbv(mgs)   &
     &  + qhsbv(mgs)   &
     &  + Min(0.0, qhcev(mgs))   &
     &  -qhmul1(mgs) - qhlcnh(mgs) - qscnh(mgs)  &
     &  - ffrzh*(qsplinter(mgs) + qsplinter2(mgs))
!     > - cimas0*nsplinter*(crfrzf(mgs) + crfrz(mgs))/rho0(mgs)

      end do



!
!  Frozen drops
!
      IF ( lf > 1 ) THEN
      do mgs = 1,ngscnt
      pqfwi(mgs) =    &
     &  +il5(mgs)*((1.0-ffrzh)*qrfrzf(mgs)  + (1-il3(mgs))*(1.0-ffrzh)*ifiacrg*(qiacrf(mgs)+qracif(mgs)))   &
     &  +il5(mgs)*qfdpv(mgs)   &
     &  +Max(0.0, qfcev(mgs))   &
     &  +qfacr(mgs)+qfacw(mgs)   &
     &  +qfacs(mgs)+qfaci(mgs) + (1.0-f2h)*qhcns(mgs) + (1.0-f2h)*qhcni(mgs)

      pqfwd(mgs) =     &
     &   qfshr(mgs)                &    !null at this point when wet graupel included
     &  +(1-il5(mgs))*qfmlr(mgs)   &    !null at this point when wet graupel included
!     >  +il5(mgs)*qhsbv(mgs)   &
     &  + qfsbv(mgs)   &
     &  + Min(0.0, qfcev(mgs))   &
     &  -qfmul1(mgs) - qhlcnf(mgs)   &
     &  - (1.0-ffrzh)*(qsplinter(mgs) + qsplinter2(mgs))
!     > - cimas0*nsplinter*(crfrzf(mgs) + crfrz(mgs))/rho0(mgs)

            IF ( .false. .and.  ny <= 2 .and. ( qfmlr(mgs) ) /= 0.0 ) THEN
             write(91,*)  'i,k,temcg = ',igs(mgs),kgs(mgs),temcg(mgs)
!             write(0,*) 'pzfwi,d = ',pzfwi(mgs),pzfwd(mgs),dtp*( pzfwi(mgs) + pzfwd(mgs) ),zx(mgs,lf)
!             write(0,*) 'pqhwi,d = ',pqhwi(mgs),pqhwd(mgs),dtp*( pqhwi(mgs) + pqhwd(mgs) ),qx(mgs,lf)
!             write(0,*) 'pcfwi,d = ',pcfwi(mgs),pcfwd(mgs),dtp*( pcfwi(mgs) + pcfwd(mgs) ),cx(mgs,lf)
!             write(0,*) 'pzfwi:'
             write(91,*) 'qfmlr,cfmlr = ',qfmlr(mgs) , cfmlr(mgs), fwvent(mgs)
             write(91,*) 'qf,cf,zf = ',qx(mgs,lf),cx(mgs,lf),zx(mgs,lf)
             write(91,*) 'xdia,ax,bx = ',xdia(mgs,lf,1),axx(mgs,lf),bxx(mgs,lf)
!        fwventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lf,1)**(0.5 + 0.5*bxx(mgs,lf)))*Sqrt(axx(mgs,lf)*rhovt(mgs)) 
!        fwvent(mgs) =  0.78*x + y*fwventy(mgs) 
!             write(0,*) zfacw(mgs),zfacr(mgs),zfacs(mgs),zfaci(mgs) 
!             write(0,*) Max( 0.0, zfdsv(mgs) )
!             write(0,*) 'pzfwd:'
!             write(0,*) zfmlr(mgs), zfshr(mgs),zhlcnf(mgs)
              
            ENDIF

      end do
      ENDIF

!
!  Hail
!
      IF ( lhl .gt. 1 ) THEN

      do mgs = 1,ngscnt
      pqhli(mgs) =    &
     &  +il5(mgs)*(qhldpv(mgs) + ((1.0-ifrzg)*qrfrzf(mgs) + (1.0-ifiacrg)*(qiacrf(mgs)+ qracif(mgs))))   &
     &  +Max(0.0, qhlcev(mgs))   &
     &  +qhlacr(mgs)+qhlacw(mgs)   &
     &  +qhlacs(mgs)+qhlaci(mgs)   &
     &  + qhlcnf(mgs) &
     &  + qhlcnh(mgs)
      pqhld(mgs) =     &
     &   qhlshr(mgs)    &
     &  +(1-il5(mgs))*qhlmlr(mgs)    &
!     >  +il5(mgs)*qhlsbv(mgs)   &
     &  + qhlsbv(mgs)   &
     &  + Min(0.0, qhlcev(mgs))   &
     &  -qhlmul1(mgs) - qhcnhl(mgs)

      IF ( imixedphase == 0 ) THEN
      frac = 0.0
      IF ( qx(mgs,lhl) + dtp*(pqhli(mgs) + pqhld(mgs)) < 0.0 ) THEN
        ! rescale depletion

         frac = (-qx(mgs,lhl) + pqhli(mgs)*dtp)/(pqhld(mgs)*dtp)
         
         qhlmlr(mgs) = frac*qhlmlr(mgs)
         qhlsbv(mgs) = frac*qhlsbv(mgs)
         qhcnhl(mgs) = frac*qhcnhl(mgs)
         qhlmul1(mgs) = frac*qhlmul1(mgs)
         IF ( qhlcev(mgs) < 0.0 ) qhlcev(mgs) = frac*qhlcev(mgs)
           
         pqhld(mgs) = frac*pqhld(mgs)
           
      ENDIF
      ENDIF


      end do
      
      ENDIF ! lhl

      ELSEIF ( warmonly < 0.8 ) THEN
!
!  Graupel
!
      do mgs = 1,ngscnt
      pqhwi(mgs) =    &
     &  +il5(mgs)*ifrzg*(qrfrzf(mgs) )   &
     &  +il5(mgs)*(qhdpv(mgs))   &
     &  +qhacr(mgs)+qhacw(mgs)   
      pqhwd(mgs) =     &
     &   qhshr(mgs)                &    !null at this point when wet graupel included
     &  - qhlcnh(mgs)   &
     &  - qhmul1(mgs)   &
     &  - qsplinter(mgs) - qsplinter2(mgs) &
     &  +(1-il5(mgs))*qhmlr(mgs)        !null at this point when wet graupel included
       end do

!
!  Hail
!
      IF ( lhl .gt. 1 ) THEN

      do mgs = 1,ngscnt
      pqhli(mgs) =    &
     &  +il5(mgs)*(qhldpv(mgs) ) & ! + (1.0-ifrzg)*(qiacrf(mgs)+qrfrzf(mgs)  + qracif(mgs)))   &
     &  +il5(mgs)*(1.0-ifrzg)*(qrfrzf(mgs) )  &
     &  +qhlacr(mgs)+qhlacw(mgs)   &
!     &  +qhlacs(mgs)+qhlaci(mgs)   &
     &  + qhlcnf(mgs) &
     &  + qhlcnh(mgs)
      pqhld(mgs) =     &
     &   qhlshr(mgs)    &
     &  +(1-il5(mgs))*qhlmlr(mgs)    &
!     >  +il5(mgs)*qhlsbv(mgs)   &
     &  + qhlsbv(mgs)   &
     &  -qhlmul1(mgs) - qhcnhl(mgs)

      end do

      ENDIF ! lhl

      ENDIF ! warmonly

!
!  Liquid water on snow and graupel 
!

      vhmlr(:) = 0.0
      vhlmlr(:) = 0.0
      vhfzh(:) = 0.0
      vhlfzhl(:) = 0.0
      vfmlr(:) = 0.0
      vffzf(:) = 0.0

      IF ( mixedphase ) THEN
      ELSE ! set arrays for non-mixedphase graupel
      
!        vhshdr(:) = 0.0
        vhmlr(:) = qhmlr(:) ! not actually volume, but treated as q in rate equation
!        vhsoak(:) = 0.0

!        vhlshdr(:) = 0.0
        vhlmlr(:) = qhlmlr(:) ! not actually volume, but treated as q in rate equation
!        vhlmlr(:) = rho0(:)*qhlmlr(:)/xdn(:,lhl) 
!        vhlsoak(:) = 0.0

        IF ( lf > 1 ) vfmlr(:) = qfmlr(:) ! not actually volume, but treated as q in rate equation
      ENDIF  ! mixedphase



!
!  frozen drop reflectivity
!
      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'frozen drop reflectivity'

      do mgs = 1,ngscnt
      
!      zhmlr(mgs) = 0.0
!      zhshr(mgs) = 0.0
!      zhmlrr(mgs) = 0.0
!      zhshrr(mgs) = 0.0
      zfdsv(mgs) = 0.0
      ziacr(mgs) = 0.0
      ziacrf(mgs) = 0.0
      zfacs(mgs) = 0.0
      zfaci(mgs) = 0.0
      pzfwi(mgs) = 0.0
      pzfwd(mgs) = 0.0
      
      ENDDO

      IF ( lzf .gt. 1 ) THEN ! 
      do mgs = 1,ngscnt
      
      
      IF ( qx(mgs,lf) .gt. qxmin(lf) .and. cx(mgs,lf) .gt. 0.0 ) THEN
          tmp = qx(mgs,lf)/cx(mgs,lf)
          alp = Max( alphamin, alpha(mgs,lf) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lf) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

           zfaci(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( tmp ) * qfaci(mgs) )
           zfacs(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( tmp ) * qfacs(mgs) )
        
        IF ( .not. mixedphase  .and. ibinhmlr < 1 ) THEN
         zfmlr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),xdn(mgs,lf),qx(mgs,lf), &
                           cx(mgs,lf),cfmlr(mgs),qfmlr(mgs))
!        zfmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*tmp * qfmlr(mgs) - tmp**2 * cfmlr(mgs)  )
        ENDIF
        
         ! leave zero and combine with zfacr and zfacw because sum is what matters
         zfshr(mgs) = 0.0 ! zrateqn(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),xdn(mgs,lf),qx(mgs,lf), &
                          ! cx(mgs,lf),cfshr(mgs),qfshr(mgs))
!        zfshr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*tmp * qfshr(mgs) - tmp**2 * cfshr(mgs)  )

!        IF ( lzr > 0 .and. qfshr(mgs) /= 0.0 .and. cfshrr(mgs) /= 0.0 .and. ibinhmlr < 1 ) THEN
        IF ( lzr > 0 .and. qfshr(mgs) /= 0.0 .and. cfshrr(mgs) /= 0.0 ) THEN



         IF ( temg(mgs) >= tfr ) THEN
           IF ( (shedalp + alpha(mgs,lf))*xdia(mgs,lf,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of hail
             z1 = g1*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qfshr(mgs)**2/ cfshrr(mgs)  ) 
           ELSE
             z1 = g1shr*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qfshr(mgs)**2/ cfshrr(mgs)  ) ! should this be g1shr?
           ENDIF
           zfshrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qfshr(mgs)**2/ cfshrr(mgs)  ) ! should this be g1shr?
!           zfshrr(mgs) = Max( z1, zfshrr(mgs))
         ELSE
          zfshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qfshr(mgs)**2/ cfshrr(mgs)  )
         ENDIF
         
         zfshrr(mgs) = Min( 0.0, zfshrr(mgs) )
        ENDIF

        IF ( zfshr(mgs) > 0.0 ) THEN
          write(0,*) 'Problem with zfshr! zfshr,qfshr,cfshr = ',zfshr(mgs),qfshr(mgs),cfshr(mgs)
          write(0,*) 'g1,tmp, qx,cx,zx = ',g1,tmp,qx(mgs,lf),cx(mgs,lf),zx(mgs,lf)
          write(0,*) ( 2.*tmp * qfshr(mgs) - tmp**2 * cfshr(mgs)  ),  2.*tmp * qfshr(mgs), - tmp**2 * cfshr(mgs)
          write(0,*) 'temcg = ',temcg(mgs),'cfshr recalc = ',(cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfshr(mgs)
          
          STOP
        ENDIF


!        zfshr(mgs) =  (xdn0(lr)/(xdn(mgs,lf)))**2*( zx(mgs,lf) * qfshr(mgs) )
        
        qtmp = qfdpv(mgs) + qfcev(mgs)
        ctmp = cfdpv(mgs) + cfcev(mgs)

        zfdsv(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

          alp = Max( alphahacx, alpha(mgs,lf) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lf) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

          IF ( .true. ) THEN  ! {
          IF ( qfacr(mgs) + qfacw(mgs) .gt. 0.0 ) THEN

            qtmp = qfacr(mgs) + qfacw(mgs) + qfshr(mgs) - qfmul1(mgs)
            ctmp = cfshr(mgs)

            zfacr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),xdn(mgs,lf),qx(mgs,lf), &
                           cx(mgs,lf),ctmp,qtmp)
!          zfacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( qx(mgs,lf)/cx(mgs,lf)) * qfacr(mgs) )



!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lf)+dtp*qfacr(mgs))**2)/(cx(mgs,lf))

          IF ( z > zx(mgs,lf) ) THEN
!            zfacr(mgs) = (z - zx(mgs,lf))*dtpinv
          ELSE
!            zfacr(mgs) = 0.0
          ENDIF
          ENDIF

!        zfacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qfacr(mgs) )
!        zfacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qfacr(mgs) - tmp**2 * cfacr(mgs) )

!          alp = Max( 1.0, alpha(mgs,lf)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          IF ( qfacw(mgs) .gt. 0.0 ) THEN
!          zfacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lf)/cx(mgs,lf)) * qfacw(mgs) )
           ! zfacw is combined with zfacr (and zfshr)
            zfacw(mgs) = 0.0 ! zrateq(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),xdn(mgs,lf),qx(mgs,lf), &
                          ! cx(mgs,lf),qfacw(mgs))
!          zfacw(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( qx(mgs,lf)/cx(mgs,lf)) * qfacw(mgs) )

!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lf)+dtp*(qfacw(mgs)-qfmul1(mgs)))**2)/(cx(mgs,lf))
          IF ( z > zx(mgs,lf) ) THEN
!            zfacw(mgs) = (z - zx(mgs,lf))*dtpinv
          ENDIF
          ENDIF

          ELSE ! } { ! this is not used because of the 'true' above

          IF ( qfacw(mgs) .gt. 0.0 .or. qfacr(mgs) .gt. 0.0 ) THEN
          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lf)+dtp*(qfacr(mgs) + qfacw(mgs)-qfmul1(mgs)))**2)/(cx(mgs,lf))
!          zfacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lf)/cx(mgs,lf)) * qfacw(mgs) )
          IF ( z > zx(mgs,lf) ) THEN
            zfacw(mgs) = (z - zx(mgs,lf))*dtpinv
          ENDIF
          ENDIF

          ENDIF ! }

          IF ( qhlcnf(mgs) .gt. 0.0 .and. ihlcnh < 2  ) THEN
           zhlcnf(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),xdn(mgs,lf),qx(mgs,lf), &
                           cx(mgs,lf),chlcnf(mgs),qhlcnf(mgs))
!           zhlcnf(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lf)))**2*( 2.*( tmp ) * qhlcnf(mgs) - tmp**2 * chlcnf(mgs) )
           !IF ( zhlcnf(mgs) < 0.0 ) THEN
           !  write(0,*) 'zhlcnf < 0? z,q,c = ', zhlcnf(mgs), qhlcnf(mgs),chlcnf(mgs)
           !  write(0,*) 'tmp,term1,term2: ',tmp,2.*( tmp ) * qhlcnf(mgs),tmp**2 * chlcnf(mgs)
           !ENDIF
           zhlcnf(mgs) = Max( 0.0, zhlcnf(mgs) )
          ENDIF
      ENDIF
! qsplinter(mgs)
      IF ( (1.0-ffrzh)*qiacrf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            IF ( imurain == 3 ) THEN
            ! note that 3.6476 = (6/pi)**2
            ziacr(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.))*  &
     &           ( 2.*tmp * qiacr(mgs) - tmp**2 * ciacr(mgs)  )
            ELSE ! imurain == 1 
            ziacr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lr),rho0(mgs),xdn0(lr),qx(mgs,lr), &
                           cx(mgs,lr),ciacr(mgs),qiacr(mgs))
!            ziacr(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2)*  &
!     &           ( 2.*tmp * qiacr(mgs) - tmp**2 * ciacr(mgs)  )
            ENDIF
            ziacr(mgs) = Min( ziacr(mgs), zxmxd(mgs,lr) )
!            ziacrf(mgs) = (xdn(mgs,lr)/xdn(mgs,lf))**2 * ziacr(mgs)
            ziacrf(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lf),rho0(mgs),rhofrz,qx(mgs,lf), &
                           cx(mgs,lf),ciacrf(mgs),qiacrf(mgs))
!            ziacrf(mgs) = (xdn(mgs,lr)/xdnmx(lf))**2 * ziacr(mgs)
!            z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * (qiacrf(mgs) - qsplinter(mgs)) - tmp**2 * ciacrf(mgs)  )
!            ziacrf(mgs) = Min(  ziacrf(mgs), z )
      ENDIF
      
      IF ( (1.0-ffrzh)*qrfrzf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 ) THEN
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            IF ( imurain == 3 ) THEN
            zrfrz(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.)) * &
     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            zrfrzf(mgs) = (xdn(mgs,lr)/xdn(mgs,lf))**2 * zrfrz(mgs)
            ELSEIF ( imurain == 1 .and. ibiggopt /= 2 ) THEN
!            zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
!     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrz(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lr),rho0(mgs),xdn0(lr),qx(mgs,lr), &
                           cx(mgs,lr),crfrz(mgs),qrfrz(mgs))

!             zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
!      &         ( 2.*tmp * qrfrz(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrzf(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lr),rho0(mgs),rhofrz,qx(mgs,lr), &
                           cx(mgs,lr),crfrzf(mgs),qrfrzf(mgs))
!            zrfrzf(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(rhofrz**2) * &
!     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            ENDIF
            zrfrz(mgs) = Min( zrfrz(mgs), Max(0.4,qrfrz(mgs)/qx(mgs,lr))*zx(mgs,lr)*dtpinv )
      ! change this to be alpha=0?
      
      ENDIF


      pzfwi(mgs) =   &
     &  +(1.0-ffrzh)*(ifrzg*zrfrzf(mgs)   &
     & +il5(mgs)*ifiacrg*(ziacrf(mgs) ) )   &
     & + zfacw(mgs)   &
     & + zfacr(mgs)   &
     & + zfacs(mgs)   &
     & + zfaci(mgs)   &
     &  + (1.0 - f2h)*zhcni(mgs) + (1.0 - f2h)*zhcns(mgs) &
     & + Max( 0.0, zfdsv(mgs) )

      pzfwd(mgs) = 0.0   &
     & + (1-il5(mgs))*zfmlr(mgs)   &
     & + zfshr(mgs)   &
!     >  + il5(mgs)*cfsbv(mgs)   &
     &  + Min( 0.0, zfdsv(mgs) )   &
     &  - il5(mgs)*zhlcnf(mgs)


!           IF ( igs(mgs) == 44 .and. kgs(mgs) == 23 .or. dtp*( pqhwi(mgs) + pqhwd(mgs) ) > qxmin(lf) ) THEN
            IF ( .false. .and. ny <= 2 .and. temcg(mgs) > -5.0 .and. ( pqhwi(mgs) + pqhwd(mgs) ) /= 0.0 ) THEN
!            IF (  ny <= 2 .and. ( qfmlr(mgs) ) /= 0.0 ) THEN
             write(0,*)  'i,k,temcg = ',igs(mgs),kgs(mgs),temcg(mgs)
             write(0,*) 'pzfwi,d = ',pzfwi(mgs),pzfwd(mgs),dtp*( pzfwi(mgs) + pzfwd(mgs) ),zx(mgs,lf)
             write(0,*) 'pqhwi,d = ',pqhwi(mgs),pqhwd(mgs),dtp*( pqhwi(mgs) + pqhwd(mgs) ),qx(mgs,lf)
             write(0,*) 'pcfwi,d = ',pcfwi(mgs),pcfwd(mgs),dtp*( pcfwi(mgs) + pcfwd(mgs) ),cx(mgs,lf)
             write(0,*) 'pzfwi:'
             write(0,*) zrfrzf(mgs) , ziacrf(mgs)
             write(0,*) zfacw(mgs),zfacr(mgs),zfacs(mgs),zfaci(mgs) 
             write(0,*) Max( 0.0, zfdsv(mgs) )
             write(0,*) 'pzfwd:'
             write(0,*) zfmlr(mgs), zfshr(mgs),zhlcnf(mgs)
             write(0,*) Min( 0.0, zfdsv(mgs) ) 

           ENDIF

      end do

      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'end frozen drop reflectivity'
      
      ENDIF


!
!  Graupel reflectivity
!
      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'graupel reflectivity'

      do mgs = 1,ngscnt
      
!      zhmlr(mgs) = 0.0
!      zhshr(mgs) = 0.0
!      zhmlrr(mgs) = 0.0
!      zhshrr(mgs) = 0.0
      zhdsv(mgs) = 0.0
!      IF ( lf < 1 ) THEN
      IF ( ffrzh > 0.0 ) THEN
      ! only initialize if frozen drops are turned off, otherwise is already set above
      ! If ffrzh = 0, then ziacrf is zeroed out for graupel and can leave value set for diagnostics
      ziacr(mgs) = 0.0
      ziacrf(mgs) = 0.0
      ENDIF
!      ENDIF
      zhcns(mgs) = 0.0
      zhcni(mgs) = 0.0
      zhacs(mgs) = 0.0
      zhaci(mgs) = 0.0
      
      ENDDO

      IF ( lzh .gt. 1 ) THEN ! 
      do mgs = 1,ngscnt
      
      
      IF ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) .gt. 0.0 ) THEN
          tmp = qx(mgs,lh)/cx(mgs,lh)
          alp = Max( alphamin, alpha(mgs,lh) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!          g1r = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)

           zhaci(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhaci(mgs) )
           zhacs(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhacs(mgs) )
        
        IF ( .not. mixedphase  .and. ibinhmlr < 1 ) THEN
         zhmlr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),xdn(mgs,lh),qx(mgs,lh), &
                           cx(mgs,lh),chmlr(mgs),qhmlr(mgs))
        ! zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
        ENDIF
        
         ! combined with zhacr
         zhshr(mgs) = 0.0 !zrateqn(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),xdn(mgs,lh),qx(mgs,lh), &
                          ! cx(mgs,lh),chshr(mgs),qhshr(mgs))
!        zhshr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs)  )

!        IF ( lzr > 0 .and. qhshr(mgs) /= 0.0 .and. chshrr(mgs) /= 0.0 .and. ibinhmlr < 1 ) THEN
        IF ( lzr > 0 .and. qhshr(mgs) /= 0.0 .and. chshrr(mgs) /= 0.0 ) THEN
!         IF ( temg(mgs) > tfr + 2.0 ) THEN
!           zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshrr(mgs)  )
!           IF ( zhshrr(mgs) > 0. ) THEN
!             zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs) )
!           ENDIF
!           z1 = g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
!           zhshrr(mgs) = Max( z1, zhshrr(mgs))
!         ELSE
!          zhshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  )


         IF ( temg(mgs) >= tfr ) THEN
 !           zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshrr(mgs)  )
 !           IF ( zhshrr(mgs) > 0.0 ) THEN
 !             zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs)  )
 !           ENDIF
           IF ( (shedalp + alpha(mgs,lh))*xdia(mgs,lh,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of hail
             z1 = g1*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) 
           ELSE
             z1 = g1shr*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
           ENDIF
           zhshrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
!           zhshrr(mgs) = Max( z1, zhshrr(mgs))
         ELSE
          zhshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  )
         ENDIF
         
         zhshrr(mgs) = Min( 0.0, zhshrr(mgs) )
        ENDIF


!        zhshr(mgs) =  (xdn0(lr)/(xdn(mgs,lh)))**2*( zx(mgs,lh) * qhshr(mgs) )
        
        qtmp = qhdpv(mgs) + qhcev(mgs) + qhsbv(mgs)
        ctmp = chdpv(mgs) + chcev(mgs) + chsbv(mgs)

        zhdsv(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

          alp = Max( alphahacx, alpha(mgs,lh) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

          IF ( .true. ) THEN  ! {
          IF ( qhacr(mgs) .gt. 0.0 ) THEN
!          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )

!          g1r = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
!          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )
            qtmp = qhacr(mgs) + qhacw(mgs) + qhshr(mgs) - qhmul1(mgs)
            ctmp = chshr(mgs)
            zhacr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),xdn(mgs,lh),qx(mgs,lh), &
                           cx(mgs,lh),ctmp,qtmp)
!          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )
!          zhacrf(mgs) = g1*zhacr

          ENDIF

!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )

!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          IF ( qhacw(mgs) .gt. 0.0 ) THEN
!          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
            ! combined with zracr
            zhacw(mgs) = 0.0 !zrateq(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),xdn(mgs,lh),qx(mgs,lh), &
                             ! cx(mgs,lh),qhacw(mgs))
!           zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
          ENDIF

          ELSE ! } { ! this is not used because of the 'true' above

!           IF ( qhacw(mgs) .gt. 0.0 .or. qhacr(mgs) .gt. 0.0 ) THEN
!           z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lh)+dtp*(qhacr(mgs) + qhacw(mgs)-qhmul1(mgs)))**2)/(cx(mgs,lh))
! !          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
!           IF ( z > zx(mgs,lh) ) THEN
!             zhacw(mgs) = (z - zx(mgs,lh))*dtpinv
!           ENDIF
!           ENDIF

          ENDIF ! }

          IF ( qhlcnh(mgs) .gt. 0.0 .and. ihlcnh < 2  ) THEN
           zhlcnh(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),xdn(mgs,lh),qx(mgs,lh), &
                           cx(mgs,lh),chlcnh(mgs),qhlcnh(mgs))
          ! zhlcnh(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhlcnh(mgs) - tmp**2 * chlcnh(mgs) )
          ENDIF
      ENDIF
! qsplinter(mgs)
      IF ( ffrzh*qiacrf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            IF ( imurain == 3 ) THEN
            ! note that 3.6476 = (6/pi)**2
            ziacr(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.))*  &
     &           ( 2.*tmp * qiacrf(mgs) - tmp**2 * ciacrf(mgs)  )
            ELSE ! imurain == 1 
             ziacr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lr),rho0(mgs),xdn0(lr),qx(mgs,lr), &
                           cx(mgs,lr),ciacr(mgs),qiacr(mgs))
!             ziacr(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2)*  &
!      &           ( 2.*tmp * qiacrf(mgs) - tmp**2 * ciacrf(mgs)  )
            ENDIF
            ziacr(mgs) = Min( ziacr(mgs), zxmxd(mgs,lr) )
!            ziacrf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * ziacr(mgs)
!            ziacrf(mgs) = (xdn(mgs,lr)/xdnmx(lh))**2 * ziacr(mgs)
            ziacrf(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lh),rho0(mgs),rhofrz,qx(mgs,lh), &
                           cx(mgs,lh),ciacrf(mgs),qiacrf(mgs))
!            z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * (qiacrf(mgs) - qsplinter(mgs)) - tmp**2 * ciacrf(mgs)  )
!            ziacrf(mgs) = Min(  ziacrf(mgs), z )
      ENDIF
      
      
      
      IF ( ffrzh*qrfrzf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 ) THEN
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            IF ( imurain == 3 ) THEN
            zrfrz(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.)) * &
     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            zrfrzf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * zrfrz(mgs)
            ELSEIF ( imurain == 1 .and. ibiggopt /= 2 ) THEN
!            zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
!     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
     &         ( 2.*tmp * qrfrz(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrzf(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(rhofrz**2) * &
     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            ENDIF
            zrfrz(mgs) = Min( zrfrz(mgs), Max(0.4,qrfrz(mgs)/qx(mgs,lr))*zx(mgs,lr)*dtpinv )
!            zrfrzf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * zrfrz(mgs)
!            zrfrzf(mgs) = (xdn(mgs,lr)/xdnmx(lh))**2 * zrfrz(mgs)
!            z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * (qrfrzf(mgs)-qsplinter2(mgs)) - tmp**2 * crfrzf(mgs)  )
!             zrfrzf(mgs) = Min(  zrfrzf(mgs), z )
      ! change this to be alpha=0?
      ENDIF
      
      IF ( lhl > 1 .and. qhcnhl(mgs) .gt. 0.0 ) THEN
        tmp = qx(mgs,lhl)/cx(mgs,lhl)
        zhcnhl(mgs) = g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qhcnhl(mgs) - tmp**2 * chcnhl(mgs) )
        
      ENDIF
      
      IF ( qhcns(mgs) > 0.0 .and. chcns(mgs) > 0.0 .and. cx(mgs,ls) > cxmin .and. vhcns(mgs) > 0 ) THEN
        tmp = qx(mgs,ls)/cx(mgs,ls)
        r = rho0(mgs)*qhcns(mgs)/vhcns(mgs) ! density of new graupel particles
        IF ( imusnow == 3 ) THEN
        zhcns(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,ls)+2.)/(r**2*(alpha(mgs,ls)+1.)) * &
     &         ( 2.*tmp * qhcns(mgs) - tmp**2 * chcnsh(mgs)  )
        ELSE
         write(0,*) 'Value of imusnow not valid. Must be 3 (fix me for =1). imusnow = ',imusnow
        ! STOP
        ENDIF
      ENDIF

      IF ( qhcni(mgs) > 0.0 .and. chcnih(mgs) > 0.0 .and. cx(mgs,li) > cxmin .and. vhcni(mgs) > 0 ) THEN
        tmp = qx(mgs,li)/cx(mgs,li)
        r = rho0(mgs)*qhcni(mgs)/vhcni(mgs) ! density of new graupel particles
        zhcni(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,li)+2.)/(r**2*(alpha(mgs,li)+1.)) * &
     &         ( 2.*tmp * qhcni(mgs) - tmp**2 * chcnih(mgs)  )
      ENDIF
 

      pzhwi(mgs) =   &
     &  +ifrzg*ffrzh*(zrfrzf(mgs)   &
     & +il5(mgs)*ifiacrg*(ziacrf(mgs) ) )   & ! ffrzh turns this off if FD are turned on
!     : + zhcnsh(mgs) + zhcnih(mgs)   &
     & + zhacw(mgs)   &
     & + zhacr(mgs)   &
     & + zhcnhl(mgs)  &
     & + zhacs(mgs)   &
     & + zhaci(mgs)   &
     &  + f2h*zhcni(mgs) + f2h*zhcns(mgs) &
     & + Max( 0.0, zhdsv(mgs) )

      pzhwd(mgs) = 0.0   &
     & + (1-il5(mgs))*zhmlr(mgs)   &
     & + zhshr(mgs)   &
     &  + Min( 0.0, zhdsv(mgs) )   &
     &  - il5(mgs)*zhlcnh(mgs)


!        IF ( zhcnhl(mgs) < 0.0 ) THEN
!          write(0,*) 'Problem with zhcnhl! zhcnhl,qhcnhl,chcnhl = ',zhcnhl(mgs),qhcnhl(mgs),chcnhl(mgs)
!          write(0,*) 'g1,tmp = ',g1x(mgs,lhl),tmp
!          write(0,*) ( 2.*( tmp ) * qhcnhl(mgs) - tmp**2 * chcnhl(mgs) )
!          
!!          STOP
!        ENDIF
      end do

      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'end graupel reflectivity'
      
      ENDIF

!
!  Hail reflectivity
!

      do mgs = 1,ngscnt
      
      zhldsv(mgs) = 0.0
      zhlacr(mgs) = 0.0
      zhlacw(mgs) = 0.0
      
      ENDDO

      IF ( lzhl .gt. 1 .or. ( lzr > 1 .and. lnhl > 1 ) ) THEN ! also run for 2-moment hail for 3-moment rain sources

      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'hail reflectivity'

      do mgs = 1,ngscnt
      
      IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) .gt. 0.0 ) THEN
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          alp = Max( alphamin, alpha(mgs,lhl) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        IF ( .not. mixedphase .and. qhlmlr(mgs) /= 0.0 .and. chlmlr(mgs) /= 0.0 .and. ibinhlmlr < 1 ) THEN
         zhlmlr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lhl),rho0(mgs),xdn(mgs,lhl),qx(mgs,lhl), &
                           cx(mgs,lhl),chlmlr(mgs),qhlmlr(mgs))
!         zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
        ENDIF
        ! combine zhlshr into zhlacr below
        zhlshr(mgs) =  0.0 ! g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshr(mgs)  )
        IF ( lzr > 1 .and. qhlshr(mgs) /= 0.0 .and. chlshrr(mgs) /= 0.0 ) THEN
         IF ( temg(mgs) >= tfr ) THEN
 !           zhlshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshrr(mgs)  )
 !           IF ( zhlshrr(mgs) > 0.0 ) THEN
 !             zhlshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshr(mgs)  )
 !           ENDIF
           IF ( (shedalp + alpha(mgs,lhl))*xdia(mgs,lhl,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of hail
             z1 = g1*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) 
           ELSE
             z1 = g1shr*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) ! should this be g1shr?
           ENDIF
           zhlshrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) ! should this be g1shr?
!           zhlshrr(mgs) = Max( z1, zhlshrr(mgs))
         ELSE
          zhlshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  )
         ENDIF

          zhlshrr(mgs) = Min( 0.0, zhlshrr(mgs) )
        ENDIF

!        zhlshr(mgs) = Min( 0.0, zhlshr(mgs) )

!        zhlshr(mgs) =  (xdn0(lr)/(xdn(mgs,lhl)))**2*( zx(mgs,lhl) * qhlshr(mgs) )
        
        qtmp = qhldpv(mgs) + qhlcev(mgs)
        ctmp = chldpv(mgs) + chlcev(mgs)
        
        zhldsv(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

          alp = Max( alphahacx, alpha(mgs,lhl) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

          IF ( .true. ) THEN ! {
          IF ( qhlacr(mgs) .gt. 0.0 ) THEN
!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*qhlacr(mgs))**2)/(cx(mgs,lhl))
            qtmp = qhlacr(mgs) + qhlacw(mgs) + qhlshr(mgs) - qhlmul1(mgs)
            ctmp = chlshr(mgs)
            zhlacr(mgs) = zrateqn(dtpinv,dtp,g1x(mgs,lhl),rho0(mgs),xdn(mgs,lhl),qx(mgs,lhl), &
                           cx(mgs,lhl),ctmp,qtmp)
!          zhlacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qhlacr(mgs) )
!          zhlacr(mgs) = Min( zxmxd(mgs,lr), zhlacr(mgs) )
          
!          IF ( z > zx(mgs,lhl) ) THEN
!            zhlacr(mgs) = (z - zx(mgs,lhl))*dtpinv
!          ELSE
!            zhlacr(mgs) = 0.0
!          ENDIF
          ENDIF

!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )

!           IF ( qhlacw(mgs) .gt. 0.0 ) THEN
!           alp = Max( 3.0, alpha(mgs,lhl)+1. )
!           g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!           
! !          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*(qhlacw(mgs)-qhlmul1(mgs)))**2)/(cx(mgs,lhl))
! !          zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lhl)/cx(mgs,lhl)) * qhlacw(mgs) )
!           zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlacw(mgs) )
! 
! !          IF ( z > zx(mgs,lhl) ) THEN
! !            zhlacw(mgs) = (z - zx(mgs,lhl))*dtpinv
! !          ENDIF
!           g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!           ENDIF
          
          ELSE ! }  .false. {

!           IF ( qhlacw(mgs) .gt. 0.0 .or. qhlacr(mgs) .gt. 0.0 ) THEN
!           z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*(qhlacr(mgs) + qhlacw(mgs)-qhlmul1(mgs)))**2)/(cx(mgs,lhl))
! !          zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lhl)/cx(mgs,lhl)) * qhlacw(mgs) )
!           IF ( z > zx(mgs,lhl) ) THEN
!             zhlacw(mgs) = (z - zx(mgs,lhl))*dtpinv
!           ENDIF
!           ENDIF
          
          ENDIF ! }
        
      ENDIF
! qsplinter(mgs)
      
      IF ( lzhl > 1 ) THEN
      pzhli(mgs) = ffrzh*(((1.0-ifrzg)*zrfrzf(mgs)   &
     & +il5(mgs)*(1.0-ifiacrg)*ziacrf(mgs) )) &
     &  + il5(mgs)*zhlcnh(mgs)   &
     &  + il5(mgs)*zhlcnf(mgs) &
     & + zhlacw(mgs)   &
     & + zhlacr(mgs)   &
!     : + zhlacs(mgs)   &
     & + Max( 0.0, zhldsv(mgs) )

      pzhld(mgs) = 0.0   &
     & + (1-il5(mgs))*zhlmlr(mgs)   &
     & + zhlshr(mgs)   &
     & - zhcnhl(mgs)   &
     &  + Min( 0.0, zhldsv(mgs) )
      

       IF ( .not. ( -1.0 < pzhli(mgs) .and. pzhli(mgs) < 1.e20 ) ) THEN
         write(iunit,*) 'Problem with pzhli!'
         write(iunit,*) 'zhlcnh,zhlacw,zhlacr,zhldsv = ',zhlcnh(mgs),zhlacw(mgs),zhlacr(mgs),zhldsv(mgs)
       ENDIF

       IF ( .not. ( -1.0e20 < pzhld(mgs) .and. pzhld(mgs) < 1. ) ) THEN
         write(iunit,*) 'Problem with pzhld!'
         write(iunit,*) 'zhlmlr,zhlshr,zhldsv = ',zhlmlr(mgs),zhlshr(mgs),zhldsv(mgs)
       ENDIF
       
      ENDIF ! lzhl > 1
      
      end do
      
      ENDIF

!
!  rain reflectivity
!
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 11'

      IF ( lzr .gt. 1 ) THEN ! 
       
        DO mgs = 1,ngscnt
        
        zracw(mgs) = 0.0
        ! zracr(mgs) = 0.0 ! already set to zero
        zrcev(mgs) = 0.0
        zrach(mgs) = 0.0
        zrachl(mgs) = 0.0
        zsshr(mgs) = 0.0
        zsshrr(mgs) = 0.0
!        zsmlr(mgs) = 0.0
        zsmlrr(mgs) = 0.0
        zracf(mgs) = 0.0

        IF ( qx(mgs,ls) .gt. qxmin(ls) .and. ( csmlr(mgs) /= 0.0 .or. csshr(mgs) /= 0.0 .or. &
              csmlrr(mgs) /= 0.0 .or. csshrr(mgs) /= 0.0) ) THEN !{
         tmp = qx(mgs,ls)/cx(mgs,ls)
         g1 = 36.*(xnu(ls)+2.0)/((xnu(ls)+1.0)*pi**2)
        IF ( .not. mixedphase ) THEN
!          zsmlr(mgs) =  (xdn(mgs,ls)/xdn(mgs,lr))**2*g1*(rho0(mgs)/(xdn(mgs,ls)))**2* &
!     &                 ( 2.*tmp * qsmlr(mgs) - tmp**2 * csmlr(mgs)  )

          IF ( csmlrr(mgs) /= 0.0 ) THEN
            z1 = g1smlr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qsmlr(mgs)**2/ csmlrr(mgs)  )
            zsmlrr(mgs) = z1
          ENDIF
        ENDIF
        
!        zsshr(mgs) =  (xdn(mgs,ls)/xdn(mgs,lr))**2*g1*(rho0(mgs)/(xdn(mgs,ls)))**2*  &
!     &                 ( 2.*tmp * qsshr(mgs) - tmp**2 * csshr(mgs)  )

         IF ( csshrr(mgs) /= 0.0 ) THEN
          z1 = g1smlr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qsshr(mgs)**2/ csshrr(mgs)  )
          zsshrr(mgs) = z1
         ENDIF
        
        ENDIF !}
        
        IF ( .not. mixedphase ) THEN !{
          IF ( zhmlr(mgs) < 0.0 .and. chmlrr(mgs) /= 0.0 .and. ibinhmlr == 0 ) THEN !{
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )
            
!            IF ( zhmlrr(mgs) >= 0. ) THEN
!              zhmlrr(mgs) =  (xdn(mgs,lh)/xdn(mgs,lr))**2 * zhmlr(mgs)
!            ENDIF
           IF ( (shedalp + alpha(mgs,lh))*xdia(mgs,lh,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of graupel
             z1 = g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  ) 
           ELSE ! assume drops are shed off, so use either alpha for shedding or graupel alpha, whichever gives the lower g-factor (i.e., larger alpha)
             z1 = Min(g1x(mgs,lh),g1shr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  )
           ENDIF
           zhmlrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  ) 
!           zhmlrr(mgs) = Max( z1, zhmlrr(mgs))
          ENDIF !}


!          zhshrr(mgs) =  (xdn(mgs,lh)/xdn(mgs,lr))**2 * zhshr(mgs)
         IF ( lf > 1 .and. qfmlr(mgs) /= 0 .and. ibinhlmlr == 0) THEN
          tmp = qx(mgs,lf)/cx(mgs,lf)

           IF ( (shedalp + alpha(mgs,lf))*xdia(mgs,lf,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of hail
             z1 = g1x(mgs,lf)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qfmlr(mgs)**2/ cfmlrr(mgs)  ) 
           ELSE ! assume drops are shed off, so use either alpha for shedding or graupel alpha, whichever gives the lower g-factor (i.e., larger alpha)
             z1 = Min(g1x(mgs,lf),g1shr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qfmlr(mgs)**2/ cfmlrr(mgs)  )
           ENDIF
           zfmlrr(mgs) = z1

         ENDIF

         IF ( lhl > 1 .and. qhlmlr(mgs) /= 0 .and. ibinhlmlr == 0) THEN
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )

!          IF ( zhlmlrr(mgs) >= 0. ) THEN ! should be negative, if not, then use alternate calculation
!           zhlmlrr(mgs) =  (xdn(mgs,lhl)/xdn(mgs,lr))**2 * zhlmlr(mgs)
!          ENDIF

           IF ( (shedalp + alpha(mgs,lhl))*xdia(mgs,lhl,1) < sheddiam ) THEN ! if not shedding small drops, then use alpha of hail
             z1 = g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  ) 
           ELSE ! assume drops are shed off, so use either alpha for shedding or graupel alpha, whichever gives the lower g-factor (i.e., larger alpha)
             z1 = Min(g1x(mgs,lhl),g1shr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
!             z1 = g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
           ENDIF
           zhlmlrr(mgs) = z1

!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
!           zhlmlrr(mgs) = Max( z1, zhlmlrr(mgs))
!         zhlmlr(mgs) =
!          zhlshrr(mgs) =  (xdn(mgs,lhl)/xdn(mgs,lr))**2 * zhlshr(mgs)
         ENDIF
         
         ENDIF ! }

        IF ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) .gt. 0.0 ) THEN

          tmp = qx(mgs,lr)/cx(mgs,lr)
          g1 = g1x(mgs,lr) ! 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)


        IF ( qracw(mgs) > 0.0 .and. cx(mgs,lr) > 0.0 ) THEN
         zracw(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * qracw(mgs) )
        ENDIF
        
! zracr is already done in breakup section
!        IF ( ibincracr /= 2 .and. cracr(mgs) /= 0.0 .and. cx(mgs,lr) > 0.0  ) THEN
        !  zracr(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*1000.))**2*( tmp**2 * cracr(mgs) )
        ! rewrite because original can overestimate zracr if -cracr*dtp is on the order of cx (i.e.,
        !  large increase in the number of drops, which violates differential assumption
!          zracr(mgs) = dtpinv*g1x(mgs,lr)*(6.*rho0(mgs)*qx(mgs,lr)/(pi*1000.))**2 &
!                      * ( cracr(mgs) )/((cx(mgs,lr) - dtp*cracr(mgs))*(cx(mgs,lr)))
!        ENDIF

        qtmp = qrcev(mgs)
        ctmp = crcev(mgs)
        
!        IF ( .false. .or. iferwisventr == 2 ) THEN
!        zrcev(mgs) = Min(0.0, (12./(pii*xdn(mgs,lr)))*xdia(mgs,lr,1)**3*fvce(mgs)*rwcap(mgs)*rwventz(mgs) )
!        ELSE
        zrcev(mgs) = g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

        
        IF (  iferwisventr == 2 ) THEN
          vent1 = Min(0.0, (12./(pii*xdn(mgs,lr)))*xdia(mgs,lr,1)**3*fvce(mgs)*rwcap(mgs)*rwventz(mgs))
          zrcev(mgs) = Max( dble(zrcev(mgs)), vent1 )
        ENDIF
!        IF ( ny == 2 .and. igs(mgs) == 20 ) THEN
!          write(0,*) 'k,zrcevold,new,maxdep : ',kgs(mgs),zrcev(mgs),vent1,-zxmxd(mgs,lr),alpha(mgs,lr),cx(mgs,lr)
!        ENDIF


!        ENDIF
        zrcev(mgs) = Max( zrcev(mgs), -zxmxd(mgs,lr) )

        IF ( qhacr(mgs) > 0.0 ) THEN 
          zrach(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &     ( 2.*( qx(mgs,lr)/cx(mgs,lr)) * qhacr(mgs) - tmp**2 * chacr(mgs) )
          zrach(mgs) = Min( zrach(mgs), zxmxd(mgs,lr) )
         
         ENDIF

        IF ( lhl > 1 .and. qhlacr(mgs) > 0.0 ) THEN 
          zrachl(mgs) = g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*   &
     &     ( 2.*( qx(mgs,lr)/cx(mgs,lr)) * qhlacr(mgs) - tmp**2 * chlacr(mgs) )
          zrachl(mgs) = Min( zrachl(mgs), zxmxd(mgs,lr) )
         ENDIF

        IF ( lf > 1 .and. qfacr(mgs) > 0.0 ) THEN 
          zracf(mgs) = g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*   &
     &     ( 2.*( qx(mgs,lr)/cx(mgs,lr)) * qfacr(mgs) - tmp**2 * cfacr(mgs) )
          zracf(mgs) = Min( zracf(mgs), zxmxd(mgs,lr) )
         ENDIF

        
        ENDIF

         pzrwi(mgs) = zrcnw(mgs) + zracw(mgs) + Max(0.0,zracr(mgs)) &
     &    + Max( 0.,zrcev(mgs) )  &
     &  - (1-il5(mgs))*zsmlrr(mgs)   &
     &  - zsshrr(mgs)   &
     &  - (1-il5(mgs))*zhmlrr(mgs)   &
     &  - zhshrr(mgs)   &
     &  - (1-il5(mgs))*zfmlrr(mgs)   &
     &  - zfshrr(mgs)   &
     &  - (1-il5(mgs))*zhlmlrr(mgs)   &
     &  - zhlshrr(mgs)   


         pzrwd(mgs) = Min(0.0,zracr(mgs))   &
     &   +  Min(0.,zrcev(mgs) )  &
     &    - zrach(mgs)  &
     &    - zrachl(mgs)  &
     &    - zracf(mgs)  &
     &    - zrfrz(mgs)  &
     &    - il5(mgs)*(ziacr(mgs) ) 

         IF (  ny <= 2 .and.                      &
              zx(mgs,lr) + dtp*(pzrwi(mgs)+pzrwd(mgs))  <= 0.0  &
              .and. qx(mgs,lr) > 0.01e-3 .and. alpha(mgs,lr) < alphamax ) THEN
           
           g1 = g1x(mgs,lr) ! 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)

           
           write(0,*) 'zrw negative: pzrwi,pzrwd = ',dtp*pzrwi(mgs),dtp*pzrwd(mgs),igs(mgs),kgs(mgs)
           write(0,*) 'qr,cr = ',1.e3*qx(mgs,lr),cx(mgs,lr),an(igs(mgs),jgs,kgs(mgs),lr),an(igs(mgs),jgs,kgs(mgs),lnr)
           write(0,*) 'zx,z,z0 = ',zx(mgs,lr),g1*rho0(mgs)**2*qx(mgs,lr)**2/(xdn(mgs,lr)**2*cx(mgs,lr)), &
     &           an(igs(mgs),jgs,kgs(mgs),lzr)
           write(0,*) 'alpha_r,rdia1,rdia3 = ', alpha(mgs,lr),xdia(mgs,lr,1),xdia(mgs,lr,3)
           write(0,*) 'zrcev, q ',   dtp*Min(0.,zrcev(mgs) ), dtp*Min(0.0,qrcev(mgs)), dtp*Min(0.0,crcev(mgs))
           write(0,*) 'zrach ', dtp*zrach(mgs)  
           write(0,*) 'zrachl,qhlacr,chlacr ', -dtp*zrachl(mgs),qhlacr(mgs),chlacr(mgs)
           write(0,*) 'zxmxd = ',zxmxd(mgs,lr),zxmxd(mgs,lr)*dtp
           write(0,*) 'zrfrz ', dtp*zrfrz(mgs),dtp*zrfrzf(mgs), qrfrz(mgs),crfrz(mgs),crfrzf(mgs)
           write(0,*) 'ziacr,qiacr,ciacr ',  dtp*il5(mgs)*ziacr(mgs), dtp*qiacr(mgs), dtp*ciacr(mgs)
           write(0,*) 'zrcnw,zracr,zracw = ',dtp*zrcnw(mgs),dtp*zracr(mgs),dtp*zracw(mgs)
           write(0,*) 'temp, qc, rho0 = ',temcg(mgs),1.e3*qx(mgs,lc),rho0(mgs)

          IF ( zx(mgs,lr) > 0.0 ) THEN
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,10
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           write(0,*) 'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
          ENDIF
         ENDIF

         ENDIF ! ny <= 2

         IF ( zx(mgs,lr) + dtp*(pzrwi(mgs)+pzrwd(mgs))  <= 0.0  &
              .and. qx(mgs,lr) > qxmin(lr) ) THEN
           pzrwd(mgs) =  -zx(mgs,lr)*dtpinv - pzrwi(mgs)
         ENDIF

        ENDDO

      ENDIF



!
!  Snow volume
!
      IF ( lvol(ls) .gt. 1 ) THEN
      do mgs = 1,ngscnt
!      pvswi(mgs) = rho0(mgs)*( pqswi(mgs) )/xdn0(ls)

      pvswi(mgs) = rho0(mgs)*(    &
!aps     >   il5*qsfzs(mgs)/xdn(mgs,ls)   &
!aps     >  -il5*qsfzs(mgs)/xdn(mgs,lr)   &
     &  +il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs)   &
     &   + qscnvi(mgs) + (1. - ifrzs)*qiacrs(mgs) &
     &   + (1. - ifrzs)*qrfrzs(mgs)  &
     &  )/xdn0(ls)   &
     &    + (qsacr(mgs))/rimdn(mgs,ls) ) + vsacw(mgs)
!     >   + (qsacw(mgs) + qsacr(mgs))/rimdn(mgs,ls) )
      pvswd(mgs) = rho0(mgs)*( pqswd(mgs) )/xdn0(ls)  &
!     >  -qhacs(mgs)
!     >  -qhcns(mgs)
!     >  +(1-il5(mgs))*qsmlr(mgs) + qsshr(mgs)
!     >  +il5(mgs)*(qssbv(mgs))
     &   -rho0(mgs)*qsmul(mgs)/xdn0(ls)
!aps     >   +rho0(mgs)*(1-il5(mgs))*(
!aps     >             qsmlr(mgs)/xdn(mgs,ls)
!aps     >    +(qscev-qsmlr(mgs))/xdn(mgs,lr) )
      end do

!aps      IF (mixedphase) THEN
!aps        pvswd(mgs) = pvswd(mgs)
!aps     >   + rho0(mgs)*qsshr(mgs)/xdn(mgs,lr)
!aps      ENDIF

      ENDIF
!
!  Graupel volume
!
      IF ( lvol(lh) .gt. 1 ) THEN
      DO mgs = 1,ngscnt
!      pvhwi(mgs) = rho0(mgs)*( (pqhwi(mgs) )/xdn0(lh) )

!      pvhwi(mgs) = rho0(mgs)*( (pqhwi(mgs) - il5(mgs)*qrfrzf(mgs) )/xdn0(lh) !
!     :  +  il5(mgs)*qrfrzf(mgs)/rhofrz )

      pvhwi(mgs) = rho0(mgs)*(   &
     &  +il5(mgs)*( ifiacrg*ffrzh*qracif(mgs))/rhofrz   &
!erm     >  + il5(mgs)*qhfzh(mgs)/rhofrz !aps: or use xdnmx(lh)?   &
     &  + (  il5(mgs)*qhdpv(mgs)/qhdpvdn   &
     &     + (qhacs(mgs) + qhaci(mgs))/qhacidn ) )   &
     &  +   rho0(mgs)*Max(0.0, qhcev(mgs))/1000.   & ! only used in mixed phase: evaporation/condensation of liquid water coating
!     >     + qhacs(mgs) + qhaci(mgs) )/xdn0(ls) )   &
     &  + f2h*vhcns(mgs)   &
     &  + vhacr(mgs) + vhacw(mgs)  + vhfzh(mgs)   & ! qhacw(mgs)/rimdn(mgs,lh)
!     >  + vhfrh(mgs)   &
     &  + f2h*vhcni(mgs) + (ifiacrg*viacrf(mgs) + ifrzg*vrfrzf(mgs))*ffrzh
!     >  +qhacr(mgs)/raindn(mgs,lh) + qhacw(mgs)/rimdn(mgs,lh)
      
!      pvhwd(mgs) = rho0(mgs)*(pqhwd(mgs) )/xdn0(lh)

      pvhwd(mgs) = rho0(mgs)*(   &
!     >   qhshr(mgs)/xdn0(lr)   &
!     >  - il5(mgs)*qhfzh(mgs)/xdn(mgs,lr)   &
     &  +( (1-il5(mgs))*vhmlr(mgs)    &
!     >     +il5(mgs)*qhsbv(mgs)   &
     &     + qhsbv(mgs)   &
     &     + Min(0.0, qhcev(mgs))   &
     &     -qhmul1(mgs) )/xdn(mgs,lh) )   &
     &  - vhlcnh(mgs) + vhshdr(mgs) - vhsoak(mgs) - vscnh(mgs)

!      IF (mixedphase) THEN
!       pvhwd(mgs) = pvhwd(mgs) 
!     >  + rho0(mgs)*qhshr(mgs)/xdn(mgs,lh) !xdn(mgs,lr)
!      ENDIF

       IF ( lzh > 1 .and. qx(mgs,lh) > qxmin(lh) .and.  &
            vx(mgs,lh) + dtp*(pvhwi(mgs) + pvhwd(mgs)) >  rho0(mgs)*qxmin(lh)/900. ) THEN
!       Calculate change in reflectivity due to density changes

        xdn_new = rho0(mgs)*(qx(mgs,lh) + dtp*(pqhwi(mgs) + pqhwd(mgs) ))/   &
     &   (vx(mgs,lh) + dtp*(pvhwi(mgs) + pvhwd(mgs))  )

           IF ( mixedphase ) THEN 
             IF ( qxw(mgs,lh) .gt. 0.0 ) THEN
               dnmx = xdnmx(lr)
             ELSE
               dnmx = xdnmx(lh)
             ENDIF
           ELSE
             dnmx = xdnmx(lh)
           ENDIF

        xdn_new = Max( Min( xdn_new, dnmx ), xdnmn(lh) )
        
        drhodt = (xdn_new - xdn(mgs,lh))*dtpinv
        
        zhwdn(mgs) = -2.*g1x(mgs,lh)*(rho0(mgs)*qx(mgs,lh)*6.*pii )**2/(cx(mgs,lh)*xdn(mgs,lh)**3)*drhodt
        
        pzhwi(mgs) = pzhwi(mgs) + Max(0.0, zhwdn(mgs))
        pzhwd(mgs) = pzhwd(mgs) + Min(0.0, zhwdn(mgs))
        
       
       ENDIF
      IF ( .false. .and. ny .eq. 2 .and. kgs(mgs) .eq. 9 .and. igs(mgs) .eq. 19 ) THEN

      write(iunit,*)
      write(iunit,*)   'Graupel at ',igs(mgs),kgs(mgs)
!
      write(iunit,*)   il5(mgs)*qrfrzf(mgs), qrfrzf(mgs) - qrfrz(mgs)
      write(iunit,*)   il5(mgs)*qiacrf(mgs)
      write(iunit,*)   il5(mgs)*qracif(mgs)
      write(iunit,*)   'qhcns',qhcns(mgs)
      write(iunit,*)   'qhcni',qhcni(mgs)
      write(iunit,*)   il5(mgs)*(qhdpv(mgs))
      write(iunit,*)   'qhacr ',qhacr(mgs)
      write(iunit,*)   'qhacw', qhacw(mgs)
      write(iunit,*)   'qhacs', qhacs(mgs)
      write(iunit,*)   'qhaci', qhaci(mgs)
      write(iunit,*)   'pqhwi = ',pqhwi(mgs)
      write(iunit,*)
      write(iunit,*) 'qhcev',qhcev(mgs)
      write(iunit,*)
      write(iunit,*)   'qhshr',qhshr(mgs)
      write(iunit,*)  'qhmlr', (1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)   'qhsbv', qhsbv(mgs)
      write(iunit,*)   'qhlcnh',-qhlcnh(mgs)
      write(iunit,*)   'qhmul1',-qhmul1(mgs)
      write(iunit,*)   'pqhwd = ', pqhwd(mgs)
      write(iunit,*)
      write(iunit,*)  'Volume'
      write(iunit,*)
      write(iunit,*)  'pvhwi',pvhwi(mgs)
      write(iunit,*)   'vhcns', vhcns(mgs)
      write(iunit,*)  'vhacr,vhacw',vhacr(mgs), vhacw(mgs) ! qhacw(mgs)/rimdn(mgs,lh)
      write(iunit,*)  'vhcni',vhcni(mgs)
      write(iunit,*)
      write(iunit,*)  'pvhwd',pvhwd(mgs)
      write(iunit,*)  'vhlcnh,vhshdr,vhsoak ', vhlcnh(mgs),  vhshdr(mgs), vhsoak(mgs)
      write(iunit,*)  'vhmlr', vhmlr(mgs)
      write(iunit,*)
!      write(iunit,*)
!      write(iunit,*)
!      write(iunit,*)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchwi(mgs),pchwd(mgs)
      write(iunit,*)  crfrzf(mgs)
      write(iunit,*)  chcns(mgs)
      write(iunit,*)  ciacrf(mgs)


      ENDIF


      ENDDO

      ENDIF
!
!
!

!
!  Hail volume
!
      IF ( lhl .gt. 1 ) THEN
      IF ( lvol(lhl) .gt. 1 ) THEN
      DO mgs = 1,ngscnt

      pvhli(mgs) = rho0(mgs)*(   &
     &  + (  il5(mgs)*(((1.0-ifiacrg)*ffrzh*qracif(mgs))/rhofrz  + qhldpv(mgs) )   &
!     &  +    Max(0.0, qhlcev(mgs))   &
!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lhl) )   & ! xdn0(ls) )   &
!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lh) )   &  ! yes, this is 'lh' on purpose
     &     + qhlacs(mgs) + qhlaci(mgs) )/500. )   &  ! changed to 500 instead of min graupel density to keep hail density from dropping too much
     &  +   rho0(mgs)*Max(0.0, qhlcev(mgs))/1000.   &
     &  + vhlcnfhl(mgs) &
     &  + vhlcnhl(mgs) + ((1.0-ifiacrg)*ffrzh*viacrf(mgs) + (1.0-ifrzg)*ffrzh*vrfrzf(mgs))  & 
     &  + vhlacr(mgs) + vhlacw(mgs) + vhlfzhl(mgs) ! qhlacw(mgs)/rimdn(mgs,lhl)
      
      pvhld(mgs) = rho0(mgs)*(   &
     &  +(  qhlsbv(mgs)   &
     &     + Min(0.0, qhlcev(mgs))   &
     &     -qhlmul1(mgs) )/xdn(mgs,lhl) ) &
!     &   + vhlmlr(mgs)                    &
     &   + rho0(mgs)*(1-il5(mgs))*vhlmlr(mgs)/xdn(mgs,lhl)  &
     &   + vhlshdr(mgs) - vhlsoak(mgs)

       IF ( lzhl > 1 .and. qx(mgs,lhl) > qxmin(lhl) .and.  &
            vx(mgs,lhl) + dtp*(pvhli(mgs) + pvhld(mgs)) >  rho0(mgs)*qxmin(lhl)/900. ) THEN
!       Calculate change in reflectivity due to density changes

        xdn_new = rho0(mgs)*(qx(mgs,lhl) + dtp*(pqhli(mgs) + pqhld(mgs) ))/   &
     &   (vx(mgs,lhl) + dtp*(pvhli(mgs) + pvhld(mgs))  )
        
           IF ( mixedphase ) THEN 
             IF ( qxw(mgs,lhl) .gt. 0.0 ) THEN
               dnmx = xdnmx(lr)
             ELSE
               dnmx = xdnmx(lhl)
             ENDIF
           ELSE
             dnmx = xdnmx(lhl)
           ENDIF
        xdn_new = Max( Min( xdn_new, dnmx ), xdnmn(lhl) )
        
        drhodt = (xdn_new - xdn(mgs,lhl))*dtpinv
        
        zhldn(mgs) = -2.*g1x(mgs,lhl)*(rho0(mgs)*qx(mgs,lhl)*6.*pii )**2/(cx(mgs,lhl)*xdn(mgs,lhl)**3)*drhodt
        
        pzhli(mgs) = pzhli(mgs) + Max(0.0, zhldn(mgs))
        pzhld(mgs) = pzhld(mgs) + Min(0.0, zhldn(mgs))
        
       
       ENDIF

      ENDDO
      
      ENDIF
      ENDIF

!
!  frozen drop volume
!
      IF ( lf .gt. 1 ) THEN
      IF ( lvol(lf) .gt. 1 ) THEN
      DO mgs = 1,ngscnt

!      pvfwi(mgs) = rho0(mgs)*(   &
!     &  + (  il5(mgs)*(((1.0-ifiacrg)*qracif(mgs))/rhofrz  + qfdpv(mgs) )   &
!!     &  +    Max(0.0, qhlcev(mgs))   &
!!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lhl) )   & ! xdn0(ls) )   &
!!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lh) )   &  ! yes, this is 'lh' on purpose
!     &     + qfacs(mgs) + qfaci(mgs) )/500. )   &  ! changed to 500 instead of min graupel density to keep hail density from dropping too much
!     &  +   rho0(mgs)*Max(0.0, qfcev(mgs))/1000.   &
!     &  + vfcnhl(mgs) + ((1.0-ifiacrg)*viacrf(mgs) + (1.0-ifrzg)*vrfrzf(mgs))  & 
!     &  + vhlacr(mgs) + vhlacw(mgs) + vhlfzhl(mgs) ! qhlacw(mgs)/rimdn(mgs,lhl)
!      
!      pvfwd(mgs) = rho0(mgs)*(   &
!     &  +(  qfsbv(mgs)   &
!     &     + Min(0.0, qfcev(mgs))   &
!     &     -qfmul1(mgs) )/xdn(mgs,lf) ) &
!!     &   + vfmlr(mgs)                    &
!     &   + rho0(mgs)*(1-il5(mgs))*vfmlr(mgs)/xdn(mgs,lf)  &
!     &   + vfshdr(mgs) - vfsoak(mgs)


      pvfwi(mgs) = rho0(mgs)*(   &
     &  +il5(mgs)*( (1.0-ffrzh)*ifiacrg*qracif(mgs))/rhofrz   &
!erm     >  + il5(mgs)*qhfzh(mgs)/rhofrz !aps: or use xdnmx(lh)?   &
     &  + (  il5(mgs)*qfdpv(mgs)/qhdpvdn   &  ! use qhdpvdn for now
     &     + (qfacs(mgs) + qfaci(mgs))/qhacidn ) )   & ! use qhacidn for now
     &  +   rho0(mgs)*Max(0.0, qfcev(mgs))/1000.   & !  evaporation/condensation of liquid water coating 
!     >     + qhacs(mgs) + qhaci(mgs) )/xdn0(ls) )   &
     &  + vfacr(mgs) + vfacw(mgs)  + vffzf(mgs)   & ! qhacw(mgs)/rimdn(mgs,lh)
     &  + (1.0 - f2h)*(vhcni(mgs) + vhcns(mgs))   &
!     >  + vhfrh(mgs)   &
     &  +  (ifiacrg*viacrf(mgs) + ifrzg*vrfrzf(mgs))*(1.0-ffrzh)
!     >  +qhacr(mgs)/raindn(mgs,lh) + qhacw(mgs)/rimdn(mgs,lh)
      
!      pvhwd(mgs) = rho0(mgs)*(pqhwd(mgs) )/xdn0(lh)

      pvfwd(mgs) = rho0(mgs)*(   &
     &  +( (1-il5(mgs))*vfmlr(mgs)    &
!     >     +il5(mgs)*qhsbv(mgs)   &
     &     + qfsbv(mgs)   &
     &     + Min(0.0, qfcev(mgs))   &
     &     -qfmul1(mgs) )/xdn(mgs,lf) )   &
     &  - vhlcnf(mgs) + vfshdr(mgs) - vfsoak(mgs)


       IF ( lzf > 1 .and. qx(mgs,lf) > qxmin(lf) .and.  &
            vx(mgs,lf) + dtp*(pvfwi(mgs) + pvfwd(mgs)) >  rho0(mgs)*qxmin(lh)/900.) THEN
!       Calculate change in reflectivity due to density changes

        xdn_new = rho0(mgs)*(qx(mgs,lf) + dtp*(pqfwi(mgs) + pqfwd(mgs) ))/   &
     &   (vx(mgs,lf) + dtp*(pvfwi(mgs) + pvfwd(mgs))  )
        
           IF ( mixedphase ) THEN 
             IF ( qxw(mgs,lf) .gt. 0.0 ) THEN
               dnmx = xdnmx(lr)
             ELSE
               dnmx = xdnmx(lf)
             ENDIF
           ELSE
             dnmx = xdnmx(lf)
           ENDIF
        xdn_new = Max( Min( xdn_new, dnmx ), xdnmn(lf) )
        
        drhodt = (xdn_new - xdn(mgs,lf))*dtpinv
        
        zfwdn(mgs) = -2.*g1x(mgs,lf)*(rho0(mgs)*qx(mgs,lf)*6.*pii )**2/(cx(mgs,lf)*xdn(mgs,lf)**3)*drhodt
        
        pzfwi(mgs) = pzfwi(mgs) + Max(0.0, zfwdn(mgs))
        pzfwd(mgs) = pzfwd(mgs) + Min(0.0, zfwdn(mgs))

            IF ( .false. .and. ny <= 2 .and. temcg(mgs) > 3.0 .and. ( pqhwi(mgs) + pqhwd(mgs) ) > 0.0 .and. &
                   (kgs(mgs) == 1 .and. igs(mgs) < 26) ) THEN
             write(0,*)  'i,k,temcg = ',igs(mgs),kgs(mgs),temcg(mgs)
             write(0,*) 'pzfwi,d = ',pzfwi(mgs),pzfwd(mgs),dtp*( pzfwi(mgs) + pzfwd(mgs) ),zx(mgs,lf)
             write(0,*) 'pqhwi,d = ',pqhwi(mgs),pqhwd(mgs),dtp*( pqhwi(mgs) + pqhwd(mgs) ),qx(mgs,lf)
             write(0,*) 'pcfwi,d = ',pcfwi(mgs),pcfwd(mgs),dtp*( pcfwi(mgs) + pcfwd(mgs) ),cx(mgs,lf)
             write(0,*) 'pzfwi:'
             write(0,*) zrfrzf(mgs) , ziacrf(mgs)
             write(0,*) zfacw(mgs),zfacr(mgs),zfacs(mgs),zfaci(mgs) 
             write(0,*) Max( 0.0, zfdsv(mgs) )
             write(0,*) Max(0.0, zfwdn(mgs))
             write(0,*) 'pzfwd:'
             write(0,*) zfmlr(mgs), zfshr(mgs),zhlcnf(mgs)
             write(0,*) Min( 0.0, zfdsv(mgs) ) 
             write(0,*) Min(0.0, zfwdn(mgs))
             write(0,*) 'pvfwi(mgs) + pvfwd(mgs) = ',pvfwi(mgs) + pvfwd(mgs)
             write(0,*) 'xdn,xdn_new,xdnmx = ',xdn(mgs,lf),xdn_new,xdnmx(lf)
             write(0,*) 'pvfwd:'
             write(0,*) rho0(mgs)*(1-il5(mgs))*vfmlr(mgs)/xdn(mgs,lf),vfmlr(mgs),qfmlr(mgs)
             write(0,*) rho0(mgs)*qfsbv(mgs)/xdn(mgs,lf)
             write(0,*) rho0(mgs)*Min(0.0, qfcev(mgs))/xdn(mgs,lf)
             write(0,*) -rho0(mgs)*qfmul1(mgs)/xdn(mgs,lf) 
             write(0,*) - vhlcnf(mgs), vfshdr(mgs), - vfsoak(mgs)

           ENDIF
        
       
       ENDIF

      ENDDO
      
      ENDIF
      ENDIF

      if ( ndebug .ge. -10 ) then
      do mgs = 1,ngscnt
!
      ptotal(mgs) = 0.
      ptotal(mgs) = ptotal(mgs)     &
     &  + pqwvi(mgs) + pqwvd(mgs)   &
     &  + pqcwi(mgs) + pqcwd(mgs)   &
     &  + pqcii(mgs) + pqcid(mgs)   &
     &  + pqrwi(mgs) + pqrwd(mgs)   &
     &  + pqswi(mgs) + pqswd(mgs)   &
     &  + pqhwi(mgs) + pqhwd(mgs)   &
     &  + pqfwi(mgs) + pqfwd(mgs)   &
     &  + pqhli(mgs) + pqhld(mgs)
!

      IF ( .false. .and. kgs(1) > 45 .and. lf > 1 ) THEN
        IF ( pqfwi(mgs) + pqfwd(mgs) /= 0.0 ) THEN
        write(0,*) 'mgs,pqfwi(mgs),pqfwd(mgs),qx(mgs,lf)',mgs,pqfwi(mgs),pqfwd(mgs),qx(mgs,lf)
      write(0,*)   qfshr(mgs)
      write(0,*)   (1-il5(mgs))*qfmlr(mgs)
      write(0,*)   il5(mgs),qfsbv(mgs),qfcev(mgs)
      write(0,*)   -qhlcnf(mgs)
      write(0,*)   -qfmul1(mgs)
      write(0,*)   qsplinter(mgs) + qsplinter2(mgs)
      write(0,*)   'pqfwd = ', pqfwd(mgs)
        ENDIF
      ENDIF
      
      
      ENDDO
      
      do mgs = 1,ngscnt

      if ( ( (ndebug .ge. 0  ) .and. abs(ptotal(mgs)) .gt. eqtot )   &
!      if ( (  abs(ptotal(mgs)) .gt. eqtot )
!     :    .or. pqswi(mgs)*dtp .gt. 1.e-3
!     :    .or. pqhwi(mgs)*dtp .gt. 1.e-3
!     :     .or. dtp*(pqrwi(mgs)+pqrwd(mgs)) .gt. 10.0e-3
!     :     .or. dtp*(pccii(mgs)+pccid(mgs)) .gt. 1.e7
!     :     .or. dtp*(pcipi(mgs)+pcipd(mgs)) .gt. 1.e7    &
     &  .or.  .not. (ptotal(mgs) .lt. 1.0 .and.  ptotal(mgs) .gt. -1.0)   & ! this line is basically checking for NaNs
     &              ) then
      write(iunit,*) 'YIKES! ','ptotal1',mgs,igs(mgs),jgs,   &
     &       kgs(mgs),ptotal(mgs)

      write(iunit,*) 't7: ', t7(igs(mgs),jgs,kgs(mgs))
      write(iunit,*)  'cci,ccw,crw,rdia: ',cx(mgs,li),cx(mgs,lc),cx(mgs,lr),0.5*xdia(mgs,lr,1)
      write(iunit,*)  'qc,qi,qr : ',qx(mgs,lc),qx(mgs,li),qx(mgs,lr)
      write(iunit,*)  'rmas, qrcalc : ',xmas(mgs,lr),xmas(mgs,lr)*cx(mgs,lr)/rho0(mgs)
      write(iunit,*)  'vti,vtc,eiw,vtr: ',vtxbar(mgs,li,1),vtxbar(mgs,lc,1),eiw(mgs),vtxbar(mgs,lr,1)
      write(iunit,*)  'cidia,cwdia,qcmxd: ', xdia(mgs,li,1),xdia(mgs,lc,1),qcmxd(mgs)
      write(iunit,*)  'snow: ',qx(mgs,ls),cx(mgs,ls),swvent(mgs),vtxbar(mgs,ls,1),xdia(mgs,ls,1)
      write(iunit,*)  'graupel: ',qx(mgs,lh),cx(mgs,lh),hwvent(mgs),vtxbar(mgs,lh,1),xdia(mgs,lh,1)
      IF ( lhl .gt. 1 ) write(iunit,*)  'hail: ',qx(mgs,lhl),cx(mgs,lhl),hlvent(mgs),vtxbar(mgs,lhl,1),xdia(mgs,lhl,1)


      write(iunit,*)  'li: ',xdia(mgs,li,1),xdia(mgs,li,2),xmas(mgs,li),qx(mgs,li),   &
     &         vtxbar(mgs,li,1)


      write(iunit,*)  'rain cx,xv : ',cx(mgs,lr),xv(mgs,lr)
      write(iunit,*)  'temcg, w = ', temcg(mgs),wvel(mgs)

      write(iunit,*) 'v ', pqwvi(mgs) ,pqwvd(mgs)
      write(iunit,*) 'c ', pqcwi(mgs) ,pqcwd(mgs)
      write(iunit,*) 'ci', pqcii(mgs) ,pqcid(mgs)
      write(iunit,*) 'r ', pqrwi(mgs) ,pqrwd(mgs)
      write(iunit,*) 's ', pqswi(mgs) ,pqswd(mgs)
      write(iunit,*) 'h ', pqhwi(mgs) ,pqhwd(mgs)
      write(iunit,*) 'hl', pqhli(mgs) ,pqhld(mgs)
      write(iunit,*) 'f ', pqfwi(mgs) ,pqfwd(mgs)
       tmp =  pqwvi(mgs) + pqwvd(mgs)   &
     &  + pqcwi(mgs) + pqcwd(mgs)   &
     &  + pqcii(mgs) + pqcid(mgs)   &
     &  + pqrwi(mgs) + pqrwd(mgs)   &
     &  + pqswi(mgs) + pqswd(mgs)   &
     &  + pqhwi(mgs) + pqhwd(mgs)   &
     &  + pqfwi(mgs) + pqfwd(mgs)   &
     &  + pqhli(mgs) + pqhld(mgs)

      write(iunit,*) 'total = ',tmp
      write(iunit,*) 'END OF OUTPUT OF SOURCE AND SINK'

!
!  print production terms
!
      write(iunit,*)
      write(iunit,*)   'Vapor'
!
      write(iunit,*)   -Min(0.0,qrcev(mgs))
      write(iunit,*)   -il5(mgs)*qhsbv(mgs)
      write(iunit,*)   -il5(mgs)*qhlsbv(mgs)
      write(iunit,*)   -il5(mgs)*qssbv(mgs)
      write(iunit,*)   -il5(mgs)*qisbv(mgs)
      write(iunit,*)    'pqwvi= ', pqwvi(mgs)
      write(iunit,*)   -Max(0.0,qrcev(mgs))
      write(iunit,*)   -Max(0.0,qhcev(mgs))
      write(iunit,*)   -Max(0.0,qhlcev(mgs))
      write(iunit,*)   -Max(0.0,qscev(mgs))
      write(iunit,*)   -il5(mgs)*qiint(mgs)
      write(iunit,*)   -il5(mgs)*qhdpv(mgs)
      write(iunit,*)   -il5(mgs)*qhldpv(mgs)
      write(iunit,*)   -il5(mgs)*qsdpv(mgs)
      write(iunit,*)   -il5(mgs)*qidpv(mgs)
      write(iunit,*)    'pqwvd = ', pqwvd(mgs)
!
      write(iunit,*)
      write(iunit,*)   'Cloud ice'
!
      write(iunit,*)   il5(mgs)*qicicnt(mgs)
      write(iunit,*)   il5(mgs)*qidpv(mgs)
      write(iunit,*)   il5(mgs)*qiacw(mgs)
      write(iunit,*)   il5(mgs)*qwfrzc(mgs)
      write(iunit,*)   il5(mgs)*qwctfzc(mgs)
      write(iunit,*)   il5(mgs)*qicichr(mgs)
      write(iunit,*)   qhmul1(mgs)
      write(iunit,*)   qhlmul1(mgs)
      write(iunit,*)   'pqcii = ', pqcii(mgs)
      write(iunit,*)   -il5(mgs)*qscni(mgs)
      write(iunit,*)   -il5(mgs)*qscnvi(mgs)
      write(iunit,*)   -il5(mgs)*qraci(mgs)
      write(iunit,*)   -il5(mgs)*qsaci(mgs)
      write(iunit,*)   -il5(mgs)*qhaci(mgs)
      write(iunit,*)   -il5(mgs)*qhlaci(mgs)
      write(iunit,*)   il5(mgs)*qisbv(mgs)
      write(iunit,*)   (1.-il5(mgs))*qimlr(mgs)
      write(iunit,*)   -il5(mgs)*qhcni(mgs)
      write(iunit,*)   'pqcid = ', pqcid(mgs)
      write(iunit,*)   ' Conc:'
      write(iunit,*)   pccii(mgs),pccid(mgs)
      write(iunit,*)   il5(mgs),cicint(mgs)
      write(iunit,*)   cwfrzc(mgs),cwctfzc(mgs)
      write(iunit,*)   cicichr(mgs)
      write(iunit,*)   chmul1(mgs)
      write(iunit,*)   cfmul1(mgs)
      write(iunit,*)   chlmul1(mgs)
      write(iunit,*)   csmul(mgs)
!
!
!
!
      write(iunit,*)
      write(iunit,*)   'Cloud water'
!
      write(iunit,*)   'pqcwi =', pqcwi(mgs)
      write(iunit,*)   -il5(mgs)*qiacw(mgs)
      write(iunit,*)   -il5(mgs)*qwfrzc(mgs)
      write(iunit,*)   -il5(mgs)*qwctfzc(mgs)
!      write(iunit,*)   -il5(mgs)*qwfrzp(mgs)
!      write(iunit,*)   -il5(mgs)*qwctfzp(mgs)
      write(iunit,*)   -il5(mgs)*qiihr(mgs)
      write(iunit,*)   -il5(mgs)*qicichr(mgs)
      write(iunit,*)   -il5(mgs)*qipiphr(mgs)
      write(iunit,*)   -qracw(mgs)
      write(iunit,*)   -qsacw(mgs)
      write(iunit,*)   -qrcnw(mgs)
      write(iunit,*)   -qhacw(mgs)
      write(iunit,*)   -qhlacw(mgs)
      write(iunit,*)   'qfacw: ',-qfacw(mgs)
      write(iunit,*)   'pqcwd = ', pqcwd(mgs)


      write(iunit,*)
      write(iunit,*)  'Concentration:'
      write(iunit,*)   -cautn(mgs)
      write(iunit,*)   -cracw(mgs)
      write(iunit,*)   -csacw(mgs)
      write(iunit,*)   -chacw(mgs)
      write(iunit,*)   'cfacw: ',-cfacw(mgs)
      write(iunit,*)  -ciacw(mgs)
      write(iunit,*)  -cwfrzp(mgs)
      write(iunit,*)  -cwctfzp(mgs)
      write(iunit,*)  -cwfrzc(mgs)
      write(iunit,*)  -cwctfzc(mgs)
      write(iunit,*)   pccwd(mgs)
!
      write(iunit,*)
      write(iunit,*)      'Rain '
!
      write(iunit,*)      qracw(mgs)
      write(iunit,*)      qrcnw(mgs)
      write(iunit,*)      Max(0.0, qrcev(mgs))
      write(iunit,*)       -(1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qhlmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qsmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qimlr(mgs)
      write(iunit,*)       -qrshr(mgs)
      write(iunit,*)       'pqrwi = ', pqrwi(mgs)    
      write(iunit,*)        -qsshr(mgs)     
      write(iunit,*)        -qhshr(mgs)     
      write(iunit,*)        -qhlshr(mgs)
      write(iunit,*)        -il5(mgs)*qiacr(mgs),qiacr(mgs), qiacrf(mgs)
      write(iunit,*)        -il5(mgs)*qrfrz(mgs)
      write(iunit,*)        -qsacr(mgs)
      write(iunit,*)        -qhacr(mgs)
      write(iunit,*)        -qhlacr(mgs)
      write(iunit,*)        qrcev(mgs)
      write(iunit,*)       'pqrwd = ', pqrwd(mgs) 
      write(iunit,*)        'qrzfac = ', qrzfac(mgs)
!
      
      write(iunit,*)
      write(iunit,*)  'Rain concentration'
      write(iunit,*)  pcrwi(mgs) 
      write(iunit,*)    crcnw(mgs)
      write(iunit,*)    1-il5(mgs)
      write(iunit,*)   -chmlr(mgs),-csmlr(mgs)
      write(iunit,*)     -crshr(mgs)
      write(iunit,*)  pcrwd(mgs) 
      write(iunit,*)    il5(mgs)
      write(iunit,*)   -ciacr(mgs),-crfrz(mgs) 
      write(iunit,*)   -csacr(mgs),-chacr(mgs)
      write(iunit,*)   +crcev(mgs)
      write(iunit,*)   cracr(mgs)
!      write(iunit,*)   -il5(mgs)*ciracr(mgs)


      write(iunit,*)
      write(iunit,*)   'Snow'
!
      write(iunit,*)        il5(mgs)*qscni(mgs), qscnvi(mgs)
      write(iunit,*)        il5(mgs)*qsaci(mgs)
      write(iunit,*)        il5(mgs)*qrfrzs(mgs), qiacrs(mgs)
      write(iunit,*)        il5(mgs)*qiacrs(mgs),il3(mgs)*(qiacrf(mgs)+qracif(mgs)),il3(mgs),qiacrf(mgs),qracif(mgs)
      write(iunit,*)        il5(mgs)*qsdpv(mgs), qscev(mgs)
      write(iunit,*)        qsacw(mgs),qwfrzc(mgs), qwctfzc(mgs), qicichr(mgs)
      write(iunit,*)        qsacr(mgs), qscnh(mgs)
      write(iunit,*) il2(mgs)*qsacr(mgs) 
      write(iunit,*) il5(mgs)*qicicnt(mgs)*ffrzs        
      write(iunit,*) il3(mgs)*(qiacrf(mgs)+qracif(mgs))  ! only applies for ipconc <= 3
      write(iunit,*) Max(0.0, qscev(mgs))   
      write(iunit,*) qsacw(mgs) + qscnh(mgs) 
      write(iunit,*)        'pqswi = ',pqswi(mgs)
      write(iunit,*)        -qhcns(mgs)
      write(iunit,*)        -qracs(mgs)
      write(iunit,*)        -qhacs(mgs)
      write(iunit,*)        -qhlacs(mgs)
      write(iunit,*)       (1-il5(mgs))*qsmlr(mgs)
      write(iunit,*)       qsshr(mgs)
!      write(iunit,*)       qsshrp(mgs)
      write(iunit,*)       il5(mgs)*(qssbv(mgs))
      write(iunit,*)       'pqswd = ', pqswd(mgs)
      write(iunit,*)   -qracs(mgs)*(1-il2(mgs)) , qhacs(mgs) , qhlacs(mgs)   
      write(iunit,*)   -qhcns(mgs)   
      write(iunit,*)   +(1-il5(mgs))*qsmlr(mgs) , qsshr(mgs)     
      write(iunit,*)    qssbv(mgs)
      write(iunit,*)   Min(0.0, qscev(mgs))  
      write(iunit,*)   -qsmul(mgs)
!
!
      write(iunit,*)
      write(iunit,*)   'Graupel'
!
      write(iunit,*)   il5(mgs)*qrfrzf(mgs), qrfrzf(mgs) - qrfrz(mgs)
      write(iunit,*)   il5(mgs)*qiacrf(mgs)
      write(iunit,*)   il5(mgs)*qracif(mgs)
      write(iunit,*)   qhcns(mgs)
      write(iunit,*)   qhcni(mgs)
      write(iunit,*)   il5(mgs)*(qhdpv(mgs))
      write(iunit,*)   qhacr(mgs)
      write(iunit,*)   qhacw(mgs)
      write(iunit,*)   qhacs(mgs)
      write(iunit,*)   qhaci(mgs)
      write(iunit,*)   'pqhwi = ',pqhwi(mgs)
      write(iunit,*)
      write(iunit,*)   qhshr(mgs)
      write(iunit,*)   (1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)   il5(mgs),qhsbv(mgs)
      write(iunit,*)   -qhlcnh(mgs)
      write(iunit,*)   -qhmul1(mgs)
      write(iunit,*)   'pqhwd = ', pqhwd(mgs)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchwi(mgs),pchwd(mgs)
      write(iunit,*)  crfrzf(mgs)
      write(iunit,*)  chcns(mgs)
      write(iunit,*)  ciacrf(mgs)

      write(iunit,*)
      write(iunit,*)   'Frozen drops'
!
      write(iunit,*)   il5(mgs)*qrfrzf(mgs), qrfrzf(mgs) - qrfrz(mgs)
      write(iunit,*)   il5(mgs)*qiacrf(mgs)
      write(iunit,*)   il5(mgs)*qracif(mgs)
      write(iunit,*)   il5(mgs)*(qfdpv(mgs))
      write(iunit,*)   qfacr(mgs)
      write(iunit,*)   qfacw(mgs)
      write(iunit,*)   qfacs(mgs)
      write(iunit,*)   qfaci(mgs)
      write(iunit,*)   'pfhwi = ',pqfwi(mgs)
      write(iunit,*)
      write(iunit,*)   qfshr(mgs)
      write(iunit,*)   (1-il5(mgs))*qfmlr(mgs)
      write(iunit,*)   il5(mgs),qfsbv(mgs),qfcev(mgs)
      write(iunit,*)   -qhlcnf(mgs)
      write(iunit,*)   -qfmul1(mgs)
      write(iunit,*)   qsplinter(mgs) + qsplinter2(mgs)
      write(iunit,*)   'pqfwd = ', pqfwd(mgs)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pcfwi(mgs),pcfwd(mgs)
      write(iunit,*)  crfrzf(mgs)
      write(iunit,*)  ciacrf(mgs)
!
      write(iunit,*)
      write(iunit,*)   'Hail'
!
      write(iunit,*)   qhlcnh(mgs)
      write(iunit,*)   il5(mgs)*(qhldpv(mgs))
      write(iunit,*)   qhlacr(mgs)
      write(iunit,*)   qhlacw(mgs)
      write(iunit,*)   qhlacs(mgs)
      write(iunit,*)   qhlaci(mgs)
      write(iunit,*)   pqhli(mgs)
      write(iunit,*)
      write(iunit,*)   qhlshr(mgs)
      write(iunit,*)   (1-il5(mgs))*qhlmlr(mgs)
      write(iunit,*)   il5(mgs)*qhlsbv(mgs)
      write(iunit,*)   pqhld(mgs)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchli(mgs),pchld(mgs)
      write(iunit,*)  chlcnh(mgs)
!
!  Balance and checks for continuity.....within machine precision...
!
!
      write(iunit,*) 'END OF OUTPUT OF SOURCE AND SINK'
      write(iunit,*) 'PTOTAL',ptotal(mgs)
!
      end if ! ptotal out of bounds or NaN
!
      end do
!

      DO mgs = 1,ngscnt
      IF ( .not. (ptotal(mgs) .lt. eqtot .and.   &
     &            ptotal(mgs) .gt. -eqtot) ) THEN
! cmic$ guard
! C$PAR CRITICAL SECTION    
!! c$omp  critical
        iptotal = iptotal + 1
      END IF
!        write(iunit,*) igs(mgs),jy,kgs(mgs),ptotal(mgs)
        ptotalmx = Max( ptotalmx, ptotal(mgs) )
        ptotalmn = Min( ptotalmn, ptotal(mgs) )
!! c$omp  end critical
! C$PAR END CRITICAL SECTION   
! cmic$ endguard
      END DO
      end if ! ( nstep/12*12 .eq. nstep )

!
!  latent heating from phase changes (except qcw, qci cond, and evap)
!
      do mgs = 1,ngscnt
      IF ( warmonly < 0.5 ) THEN
      pfrz(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+    &
     &   qfmlr(mgs) +  &
     &   qsmlr(mgs)+qhlmlr(mgs))   & !+qhmlh(mgs))   &
     &  +il5(mgs)*(1-imixedphase)*(   &
     &   qsacw(mgs)+qhacw(mgs) + qhlacw(mgs)   &
     &  +qsacr(mgs)+qhacr(mgs) + qhlacr(mgs)   &
     &  +qsshr(mgs)   &
     &  +qhshr(mgs)   &
     &  +qfshr(mgs) + qfacr(mgs) + qfacw(mgs)   & ! normally offsets qfacw for wet growth, so is flagged by (1-imixedphase) 
     &  +qhlshr(mgs)  &
     &  +qrfrz(mgs)+qiacr(mgs)  &
     &  )  &
     &  +il5(mgs)*(qwfrz(mgs)    &
     &  +qwctfz(mgs)+qiihr(mgs)   &
     &  +qiacw(mgs))
      pmlt(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+qsmlr(mgs)+  &
     &   qfmlr(mgs) +  & 
     &   qhlmlr(mgs))    !+qhmlh(mgs))   
      ! NOTE: psub is sum of sublimation and deposition
      psub(mgs) =    &
     &   il5(mgs)*(   &
     &  + qsdpv(mgs) + qhdpv(mgs)   &
     &  +qfdpv(mgs)    &
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs) + qisbv(mgs) )   &
     &   + qssbv(mgs)  + qhsbv(mgs) &
     &  + qfsbv(mgs) &
     &  + qhlsbv(mgs)   &
     &  +il5(mgs)*(qiint(mgs))
      pvap(mgs) =    &
     &   qrcev(mgs) + qhcev(mgs) + qscev(mgs) + qhlcev(mgs) + qfcev(mgs)
      pevap(mgs) =    &
     &   Min(0.0,qrcev(mgs)) + Min(0.0,qhcev(mgs)) + Min(0.0,qscev(mgs)) + Min(0.0,qhlcev(mgs)) &
         +  Min(0.0,qfcev(mgs))
      ! NOTE: pdep is the deposition part only
      pdep(mgs) =    &
     &   il5(mgs)*(   &
     &  + qsdpv(mgs) + qhdpv(mgs)   &
     &  + qfdpv(mgs)  & 
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs)  )   & 
     &  +il5(mgs)*(qiint(mgs))
      ELSEIF ( warmonly < 0.8 ) THEN
      pfrz(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+qhlmlr(mgs))   & !+qhmlh(mgs))   &
     &  +il5(mgs)*(qhfzh(mgs)+qhlfzhl(mgs))   &
     &  +il5(mgs)*(   &
     &  +qhshr(mgs)   &
     &  +qhlshr(mgs)   &
     &  +qrfrz(mgs)+qwfrz(mgs)   &
     &  +qwctfz(mgs)+qiihr(mgs)   &
     &  +qiacw(mgs)                &
     & +qhacw(mgs) + qhlacw(mgs)   &
     & +qhacr(mgs) + qhlacr(mgs)  ) 
      psub(mgs) =  0.0 +  &
     &   il5(mgs)*(   &
     &  + qsdpv(mgs)   &
     &  + qhdpv(mgs)   &
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs) + qisbv(mgs) )   &
     &  +il5(mgs)*(qiint(mgs))
      pvap(mgs) =    &
     &   qrcev(mgs) + qhcev(mgs) + qhlcev(mgs) + qfcev(mgs) 
      ELSE
      pfrz(mgs) = 0.0
      psub(mgs) = 0.0
      pvap(mgs) = qrcev(mgs)
      ENDIF ! warmonly
      ptem(mgs) =    &
     &  (1./pi0(mgs))*   &
     &  (felfcp(mgs)*pfrz(mgs)   &
     &  +felscp(mgs)*psub(mgs)    &
     &  +felvcp(mgs)*pvap(mgs))
      thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)
      ptem2(mgs) = ptem(mgs)
      IF ( eqtset > 2 ) THEN
        pipert(mgs) = pipert(mgs) + (felfpi(mgs)*pfrz(mgs)   &
     &  +felspi(mgs)*psub(mgs)    &
     &  +felvpi(mgs)*pvap(mgs))*dtp
      ENDIF
      end do



!  production of space charge ! use values from previous time step
!
!
!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscrw(mgs) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      fscrw(mgs) = scx(mgs,lr)/qx(mgs,lr)
!      if (abs(fscrw(mgs)).gt.1) write(0,*) 'fscrw(mgs),scx(mgs,lr),qx(mgs,lr)=',fscrw(mgs),scx(mgs,lr),qx(mgs,lr)
      end if
      end do
      end if

!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fsccw(mgs) = 0.0
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
      fsccw(mgs) = scx(mgs,lc)/qx(mgs,lc)
      end if
      end do
      end if

!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscci(mgs) = 0.0
      if ( qx(mgs,li) .gt. qxmin(li) ) then
      fscci(mgs) = scx(mgs,li)/qx(mgs,li)
      end if
      end do
      end if


!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscsw(mgs) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      fscsw(mgs) = scx(mgs,ls)/qx(mgs,ls)
      end if
      end do
      end if

!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fschw(mgs) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      fschw(mgs) = scx(mgs,lh)/qx(mgs,lh)
      end if
      end do
      end if

!
!  space charge production
!  do transfer during conversion
!  accretion, etc
!
      fschl(:) = 0.
      if ( lhl .gt. 1 .and. ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      fschl(mgs) = scx(mgs,lhl)/qx(mgs,lhl)
      end if
      end do
      end if

      fscfw(:) = 0.
      if ( lf .gt. 1 .and. ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
      fscfw(mgs) = scx(mgs,lf)/qx(mgs,lf)
      end if
      end do
      end if

      schmlr(:) = 0.0
      scsmlr(:) = 0.0
      schshr(:) = 0.0
      scfshr(:) = 0.0
      scsshr(:) = 0.0
      scsdep(:) = 0.0

! "Alternative" electrification mechanisms

     IF (ipelec .ge. 1) THEN

! Drake electrification mechanism 
!   for melting hail between 1mm and 5mm diameter
!   put after melting (instead of with other charging code) in order to have a consistent measure 
!   of melt water to convert Drake estimation of charge remaining on the mixed-phase particle
!
!   as a first guess, use Drake estimate of 8 esu cm-3, or 2.67e-3 C m-3 (charge per vol of melt water).
!   if it affects anything significantly, might want to improve the scheme by using temperature/velocity info
!
!   aps (30nov08): Will be overestimate if water is shed every time step; no retention.

      IF ( ieopt .eq. 1 .or. ieopt .ge. 11 ) THEN
      ! consideration: rasmussen and heymsfield say hail up to 9 mm before shedding.
      
      
      do mgs = 1,ngscnt
      
       if ( xdia(mgs,lh,1) .gt. 1.0e-03 .and. (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1) .lt. sheddiam ) then 
       ! noting here that qhmlr <= 0, so schmlr >= 0.
        schmlr(mgs) = -altelecfac*qhmlr(mgs)*rho0(mgs)/xdn(mgs,lr)*2.67e-3 
!        schmlr(mgs) = (qx(mgs,lh)*dtpinv)*rho0(mgs)/xdn(mgs,lr)*1.5e-3
       endif ! xdia

       if ( xdia(mgs,ls,1) .gt. 1.0e-03 .and. xdia(mgs,ls,3) .lt. sheddiam ) then 
        scsmlr(mgs) = -altelecfac*qsmlr(mgs)*rho0(mgs)/xdn(mgs,lr)*2.67e-3
       end if

      end do ! mgs

!debug

      if (ieopt .eq. 1) then
      tmp=0.
      do mgs=1,ngscnt
       if ( scsmlr(mgs).ne.0. ) then
       if (tmp==0.) write(0,"('-------------- ICEZVD_GS: after Drake --------------')")
       if (tmp==0.) write(0,"(1x,'i',2x,'j',2x,'k',5x,'qsmlr',10x,'rho0',10x,'xdn',10x,'scsmlr')")
          write(0,"(3(i2,1x),4(es12.5,2x))") igs(mgs),jy,kgs(mgs),qsmlr(mgs),rho0(mgs),xdn(mgs,lr),scsmlr(mgs)
        tmp=1.
       endif
     enddo
     endif

      ENDIF ! ieopt=1/drake


      
! Dong and Hallett (1992)
!   ice undergoing depositional growth
!   like Schuur and Rutledge eqns B7-B11 for inverse exponential distribution

      IF ( ieopt .eq. 2 ) THEN
       do mgs = 1,ngscnt
        if ( qscev(mgs) .gt. 0.0 ) then
         fqt = 0.0008*(temcg(mgs)**3.)-0.0025*(temcg(mgs)**2.)-0.4477*(temcg(mgs))-1.286 ! C cm-2 s-1
         scsdep(mgs) = 1.e-12*altelecfac*xsfca(mgs,ls)*fqt ! C s-1

         if (scsdep(mgs).gt.0) scsdep(mgs) = min(scsdep(mgs), ec*cionpmxd(mgs))
         if (scsdep(mgs).lt.0) scsdep(mgs) = max(scsdep(mgs),-ec*cionnmxd(mgs))
        endif
       enddo

      ENDIF ! ieopt=2/donghallett

! Canosa and List (1993)
!   Inductive charge separation for drop collisions that result in disjection
!   First pass: use multiplier with inductive charge already calculated,
!               since we now include wet snow/graupel, disjection may be more frequent

      IF ( ieopt .eq. 3 .or. ieopt .ge. 13 ) THEN
      do il = ls,lhab    ! lr also?
      do mgs = 1,ngscnt

       if ( qx(mgs,il) .gt. qxmin(il) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( temg(mgs) .gt. tindmn .and. temg(mgs) .lt. tindmx ) then
          if ( (xdia(mgs,lr,1) .gt. 5.e-4 .and. xdia(mgs,lr,1) .lt. 2.5e-3) .and.     &
     &         (xdia(mgs,ls,1) .gt. 5.e-4 .and. xdia(mgs,ls,1) .lt. 2.5e-3) ) then


!          tmp = 2.5*pii*rhovt(mgs)*xdia(mgs,lr,1)*xdia(mgs,il,1)**2                   &
!     &             *rvt*0.10*0.36*fdgt*cx(mgs,il)                                     &  ! 10% collision eff / 36% disject eff
!     &             *( dezcomp(mgs)*costhe*(xdia(mgs,lr,1)**2)*cx(mgs,lr)*pii**2/3.0  ) ! &
!     &                 - 4.*scx(mgs,lr) )

           tmp = (pii/4.)*abs(vtxbar(mgs,lr,1)-vtxbar(mgs,il,1))*0.38*0.36*fdgt   &       ! letting erbnd=0.36 and exy(mgs,ls,lr)=0.38
      &                  *cx(mgs,lr)*cx(mgs,il)*xdia(mgs,lr,1)**2*11.*costhe*1.e-12*dezcomp(mgs)/50.e3 

!    The following did NOT work despite using the same assumptions (fall speed, distros, etc) as rest of model
!           tmp = (pii/4.)*(1./-5.73855)*(vtxbar(mgs,lr,1)-vtxbar(mgs,ls,1))*0.10*0.36*rvt*fdgt   &
!      &                  *cx(mgs,lr)*cx(mgs,il)*xdia(mgs,il,1)**2                                & ! assumes rnu=-0.8
!      &              *( (pi**2)*(1./24.)*0.908639*costhe*dezcomp(mgs)*(xdia(mgs,lr,1)**2)       &
!      &                  + 4.32685*scx(mgs,lr)/cx(mgs,lr) )

!      scxacy(mgs,il,lc) =   &
!     &  0.125*pi**3*exy(mgs,il,lc)*erbnd*fdgt*(cx(mgs,il))   &
!     &  *cx(mgs,lc)*vtxbar(mgs,il,1)*6.0/gf4p5*(xdia(mgs,lc,1)**2)   &
!     &  *( pi*eperao*costhe*dezcomp(mgs)   &
!     &      *gf3p5*(xdia(mgs,il,1)**2)   &
!     &     - (1.0/3.0)*gf1p5*scx(mgs,il)/cx(mgs,il) )

         if ( fsw(mgs) .gt. 0.00 .and. il .eq. ls ) then
          scsshr(mgs) = altelecfac*tmp
!         else if ( fhw(mgs) .gt. 0.00 .and. il .eq. lh ) then
!          schshr(mgs) = tmp
         endif
         endif !xdia
        ENDIF
       endif

      end do ! il
      end do ! mgs

      ENDIF ! ieopt=3/canosalist

!debug

      if (ieopt .eq. 3) then
      tmp=0.
      do mgs=1,ngscnt
       if ( scsshr(mgs).ne.0. ) then
       if (tmp==0.) write(0,"('-------------- ICEZVD_GS: after CL93s --------------')")
       if (tmp==0.) write(0, &
     & "(1x,'i',2x,'j',2x,'k',4x,'snow conc',4x,'rain conc',6x,'delta vt',5x,'rain diam',5x,'snow diam',6x,'dezcomp',7x,'scsshr')")
          write(0,"(3(i2,1x),7(es12.5,2x))") igs(mgs),jy,kgs(mgs),cx(mgs,ls),cx(mgs,lr),vtxbar(mgs,lr,1)-vtxbar(mgs,ls,1), &
          xdia(mgs,lr,1),xdia(mgs,ls,1),dezcomp(mgs),scsshr(mgs)
        tmp=1.
       endif
     enddo
     endif

     ENDIF ! ipelec
!
!
!
!SSSS
!
!  Charge production from mass transfer
!
!  Cloud ice
!
      IF ( ipelec > 0 ) THEN
        psccimi(1:ngscnt) = 0.0
        psccimd(1:ngscnt) = 0.0
        psccwmi(1:ngscnt) = 0.0
        psccwmd(1:ngscnt) = 0.0
        pscswmi(1:ngscnt) = 0.0
        pscswmd(1:ngscnt) = 0.0
        pscrwmi(1:ngscnt) = 0.0
        pscrwmd(1:ngscnt) = 0.0
        pschwmi(1:ngscnt) = 0.0
        pschwmd(1:ngscnt) = 0.0
        pschlmi(1:ngscnt) = 0.0
        pschlmd(1:ngscnt) = 0.0
        pscfwmi(1:ngscnt) = 0.0
        pscfwmd(1:ngscnt) = 0.0
      ENDIF
      
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccimi(mgs) =    &
     &   il5(mgs)*(fsccw(mgs)*qiacw(mgs)+fsccw(mgs)*qwfrzc(mgs)   &
     &  +fsccw(mgs)*qwctfzc(mgs)+fsccw(mgs)*qicichr(mgs))   &
!     >  +il5(mgs)*fscci(mgs)*qwacii(mgs))   &
     &  +il5(mgs)*fscrw(mgs)*(1. - ifrzs)*qrfrzs(mgs)   &
     &  +fschw(mgs)*qhmul1(mgs)   &
     &  +fschl(mgs)*qhlmul1(mgs)  &
     &  +fscfw(mgs)*qfmul1(mgs)   &
     &  +fscsw(mgs)*qsmul(mgs)
      psccimd(mgs) =    &
     &   il5(mgs)*(-fscci(mgs)*qscni(mgs)   & ! -fscci(mgs)*qwaci(mgs)   &
     &  -fscci(mgs)*qsaci(mgs))  -fscci(mgs)*qraci(mgs)   &
     &  -fscci(mgs)*qfaci(mgs)   &
     &  -fscci(mgs)*qhaci(mgs) -fscci(mgs)*qhlaci(mgs)
      end do
      end if


!
!
!  Cloud water
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccwmi(mgs) = 0.0
      psccwmd(mgs) =    &
     &  -fsccw(mgs)*qracw(mgs)-fsccw(mgs)*qrcnw(mgs)   &
     &  +il5(mgs)*(-fsccw(mgs)*qiacw(mgs)   &
     &  -fsccw(mgs)*qwctfz(mgs)    & ! -fsccw(mgs)*qwctfzp(mgs)
     &  -fsccw(mgs)*qwfrz(mgs) )   & ! -fsccw(mgs)*qwfrzp(mgs))
     &  -fsccw(mgs)*qiihr(mgs)   &
     &  - fsccw(mgs)*qfacw(mgs)   &
     &  -fsccw(mgs)*qsacw(mgs) -fsccw(mgs)*qhacw(mgs)  -fsccw(mgs)*qhlacw(mgs)
      end do
      end if

      mixedphasefac = 1.0
      IF ( mixedphase ) mixedphasefac = 0.0
!
!  Snow
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscswmi(mgs) =    &
     &  +il5(mgs)*(fscci(mgs)*qscni(mgs)+fscci(mgs)*qsaci(mgs)   &
     &  +fscrw(mgs)*ifrzs*qrfrzs(mgs))   &
     &  +fsccw(mgs)*qsacw(mgs)+fscrw(mgs)*qsacr(mgs)
      pscswmd(mgs) =    &
!     >  -fscsw(mgs)*qracs(mgs) ! -fscsw(mgs)*qwacs(mgs)   &
     &  - fscsw(mgs)*qfacs(mgs)   &
     &  - fscsw(mgs)*qhacs(mgs) - fscsw(mgs)*qhlacs(mgs)   &
     &  +(1-il5(mgs))*fscsw(mgs)*qsmlr(mgs)*mixedphasefac &
!     &   + fscsw(mgs)*qsshr(mgs) &
     &  -fscsw(mgs)*qsmul(mgs)  &
!       + fscsw(mgs)*max(0.0,qscev(mgs))  &
        + fscsw(mgs)*min(0.0,qscev(mgs))

      fscsw2(mgs) = 0.0
      tmp =  qx(mgs,ls) + dtp*(pqswi(mgs) + pqswd(mgs) - qsshr(mgs)) ! should this add back the shed rain mass? probably!
      IF ( tmp > qxmin(ls) .and. qsshr(mgs) < 0.0 ) THEN
        chgtmp = scx(mgs,ls) + dtp*(pscswmi(mgs) + pscswmd(mgs))
        fscsw2(mgs) = chgtmp/tmp
        pscswmd(mgs) =  pscswmd(mgs) + fscsw2(mgs)*qsshr(mgs)
      ENDIF

        il = ls
        IF ( Abs( fscsw2(mgs)*qsshr(mgs)) > 1.e-8  ) THEN
          write(0,*) 'Problem with scs = ',il,scx(mgs,il),pqswi(mgs),pqswd(mgs),qsshr(mgs),mgs,igs(mgs),kgs(mgs)
      write(0,*) 'qs,pqswi(mgs),pqswd(mgs),qsnew ',qx(mgs,ls),pqswi(mgs),pqswd(mgs),qx(mgs,ls) + dtp*(pqswi(mgs)+pqswd(mgs))
      write(0,*) ' -qracs(mgs)*(1-il2(mgs)) , qhacs(mgs) , qhlacs(mgs)',  -qracs(mgs)*(1-il2(mgs)) , qhacs(mgs) , qhlacs(mgs)   
      write(0,*) ' -qhcns(mgs)',  -qhcns(mgs)   
      write(0,*)  '+(1-il5(mgs))*qsmlr(mgs) , qsshr(mgs) ', (1-il5(mgs))*qsmlr(mgs) , qsshr(mgs)        !null at this point when wet snow included
      write(0,*)  'qssbv', (qssbv(mgs))   
      write(0,*)   'qscev', Min(0.0, qscev(mgs))  
      write(0,*)   '-qsmul',-qsmul(mgs)
 
      write(0,*)   'pqswi',pqswi(mgs)
!      write(0,*)   il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs)   &
!      write(0,*)    + qscnvi(mgs)                        &
!     write(0,*)   + ifrzs*(qiacrs(mgs) + qrfrzs(mgs))  &
!      write(0,*)   + il2(mgs)*qsacr(mgs))   
!      write(0,*)   + il3(mgs)*(qiacrf(mgs)+qracif(mgs)) 
      write(0,*)    'qscev',Max(0.0, qscev(mgs))   
      write(0,*)     'qsacw(mgs) , qscnh(mgs)',qsacw(mgs) , qscnh(mgs)

        write(0,*) 'rank,temperature = ',my_rank,temcg(mgs)
          write(0,*) 'tmp, chgtmp, qs, scs,fscsw = ',tmp,chgtmp,qx(mgs,ls),scx(mgs,ls),fscsw(mgs)
        write(0,*) 'tmp,chgtmp, fscsw2, qsshr,fscsw2(mgs)*qsshr(mgs) ',tmp,chgtmp,fscsw2(mgs), qsshr(mgs), fscsw2(mgs)*qsshr(mgs)
        write(0,*) 'pscswmd new,old = ',pscswmd(mgs), pscswmd(mgs) - fscsw2(mgs)*qsshr(mgs)
        write(0,*)  psccwi(mgs), psccwd(mgs)
        write(0,*)  psccii(mgs), psccid(mgs)
        write(0,*)  pscrwi(mgs), pscrwd(mgs)
        write(0,*)  pscswi(mgs), pscswd(mgs)
        write(0,*)  pschwi(mgs), pschwd(mgs)
        write(0,*)  pschli(mgs), pschld(mgs)
        write(0,*) 'FD: ',pscfwi(mgs), pscfwd(mgs)
        write(0,*)  pscpii(mgs) ,pscpid(mgs)
        write(0,*)  pscnii(mgs) ,pscnid(mgs)
        write(0,*)
        write(0,*)  psccwmi(mgs), psccwmd(mgs)
        write(0,*)  psccimi(mgs), psccimd(mgs)
        write(0,*)  pscrwmi(mgs), pscrwmd(mgs)
        write(0,*)  pscswmi(mgs), pscswmd(mgs), pscswmd(mgs) - fscsw2(mgs)*qsshr(mgs)
        write(0,*)  pschwmi(mgs), pschwmd(mgs)
        write(0,*) 'FD: ',pscfwmi(mgs), pscfwmd(mgs)
        write(0,*)  pschlmi(mgs), pschlmd(mgs)

     write(0,*)  'pscswmd(mgs)'
     write(0,*) - fscsw(mgs)*qhacs(mgs), - fscsw(mgs)*qhlacs(mgs)
     write(0,*)  +(1-il5(mgs))*fscsw(mgs)*qsmlr(mgs), fscrw(mgs)*qsshr(mgs)
     write(0,*)  -fscsw(mgs)*qsmul(mgs)
     write(0,*) qhacs(mgs), qhlacs(mgs),qsmlr(mgs),qsshr(mgs) ,qsmul(mgs)
     write(0,*)
     write(0,*)  fschw(mgs),fscsw(mgs),fschl(mgs)
     write(0,*)  fsccw(mgs),fscci(mgs),fscrw(mgs)

!      tmp =  qx(mgs,ls) + dtp*(pqswi(mgs) + pqswd(mgs) - qsshr(mgs)) ! should this add back the shed rain mass? probably!
!      chgtmp = 0.0
!      IF ( tmp > qxmin(ls) .and. qsshr(mgs) < 0.0 ) THEN
!        chgtmp = scx(mgs,ls) + dtp*(pscswmi(mgs) - pscswmd(mgs))
!       ! fscsw2(mgs) = chgtmp/tmp
!       ! pscswmd(mgs) =  pscswmd(mgs) + fscsw2(mgs)*qsshr(mgs)
!      ENDIF
!      write(0,*) 'tmp,chgtmp, fscsw2, qsshr,fscsw2(mgs)*qsshr(mgs) ',tmp,chgtmp,fscsw2(mgs), qsshr(mgs), fscsw2(mgs)*qsshr(mgs)

       ENDIF
      end do
      end if

!
!
!  Graupel
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschwmi(mgs) =    &
     &  +il5(mgs)*(ifiacrg*ffrzh*(fscci(mgs)*qracif(mgs)+fscrw(mgs)*qiacrf(mgs))   &
     &  +fscrw(mgs)*ifrzg*ffrzh*qrfrzf(mgs))   &
     &  +fscrw(mgs)*qhacr(mgs)+fsccw(mgs)*qhacw(mgs)   &
     &  +fscsw(mgs)*qhacs(mgs)+fscci(mgs)*qhaci(mgs) 
      pschwmd(mgs) =    &
     &  - fschw(mgs)*qhlcnh(mgs)   &
!     &  + fschw(mgs)*qhshr(mgs)    &
     &  +(1-il5(mgs))*fschw(mgs)*qhmlr(mgs)*mixedphasefac   &
     &  -fschw(mgs)*qhmul1(mgs)

      fschw2(mgs) = 0.0
      tmp =  qx(mgs,lh) + dtp*(pqhwi(mgs) + pqhwd(mgs) - qhshr(mgs)) ! should this add back the shed rain mass? probably!
      IF ( tmp > qxmin(lh) .and. qhshr(mgs) < 0.0 ) THEN
        chgtmp = scx(mgs,lh) + dtp*(pschwmi(mgs) + pschwmd(mgs)) ! sign error bug fix 3.4.2021 Thanks to Yichen (Jade) Cai.
        IF ( sign(1.0,chgtmp) == sign(1.0,scx(mgs,lh) ) ) THEN
          fschw2(mgs) = chgtmp/tmp
          pschwmd(mgs) =  pschwmd(mgs) + fschw2(mgs)*qhshr(mgs)
        ENDIF
      ENDIF

      end do
      end if

!
!
!  Frozen drops
!
!      &  +il5(mgs)*((1.0-ffrzh)*qrfrzf(mgs)  + (1-il3(mgs))*(1.0-ffrzh)*(qiacrf(mgs)+qracif(mgs)))   &

      pscfwmi(:) = 0.0
      pscfwmd(:) = 0.0
      fscfw2(:) = 0.0

      if ( ipelec .ge. 1 .and. itest .eq. 1 .and. lf > 1 ) then
      do mgs = 1,ngscnt
      pscfwmi(mgs) =    &
     &  +il5(mgs)*(ifiacrg*(1.0-ffrzh)*(fscci(mgs)*qracif(mgs)+fscrw(mgs)*qiacrf(mgs))   &
     &  +fscrw(mgs)*(1.0-ffrzh)*qrfrzf(mgs))   &
     &  +fscrw(mgs)*qfacr(mgs)+fsccw(mgs)*qfacw(mgs)   &
     &  +fscsw(mgs)*qfacs(mgs)+fscci(mgs)*qfaci(mgs)
      pscfwmd(mgs) =    &
     &  - fscfw(mgs)*qhlcnf(mgs)   &
!     &  + fschw(mgs)*qhshr(mgs)    &
     &  +(1-il5(mgs))*fscfw(mgs)*qfmlr(mgs)*mixedphasefac   &
     &  -fscfw(mgs)*qfmul1(mgs)

      tmp =  qx(mgs,lf) + dtp*(pqfwi(mgs) + pqfwd(mgs) - qfshr(mgs)) ! should this add back the shed rain mass? probably!
      IF ( tmp > qxmin(lf) .and. qfshr(mgs) < 0.0 ) THEN
        chgtmp = scx(mgs,lf) + dtp*(pscfwmi(mgs) + pscfwmd(mgs)) ! sign error bug fix 3.4.2021 Thanks to Yichen (Jade) Cai.
        IF ( sign(1.0,chgtmp) == sign(1.0,scx(mgs,lf) ) ) THEN
          fscfw2(mgs) = chgtmp/tmp
          pscfwmd(mgs) =  pscfwmd(mgs) + fscfw2(mgs)*qfshr(mgs)
        ENDIF
      ENDIF

      end do
      end if
!
!
!  Hail
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then

      pschlmi(:) = 0.0
      pschlmd(:) = 0.0

      IF ( lhl .gt. 1 ) THEN
      do mgs = 1,ngscnt
      pschlmi(mgs) =    &
     &  +il5(mgs)*((1.0-ifiacrg)*ffrzh*(fscci(mgs)*qracif(mgs)+fscrw(mgs)*qiacrf(mgs))   &
     &  +fscrw(mgs)*(1.0-ifrzg)*qrfrzf(mgs))   &
     &  +fscrw(mgs)*qhlacr(mgs)+fsccw(mgs)*qhlacw(mgs)   &
     &  +fscsw(mgs)*qhlacs(mgs)+fscci(mgs)*qhlaci(mgs)   &
     &  + fschw(mgs)*qhlcnh(mgs)
      pschlmd(mgs) =    &
!     &  + fschl(mgs)*qhlshr(mgs)    &
     &  +(1-il5(mgs))*fschl(mgs)*qhlmlr(mgs)*mixedphasefac    &
     &  -fschl(mgs)*qhlmul1(mgs)

       fschl2(mgs) = 0.0
      tmp =  qx(mgs,lhl) + dtp*(pqhli(mgs) + pqhld(mgs) - qhlshr(mgs)) ! should this add back the shed rain mass? probably!
      IF ( tmp > qxmin(lhl) .and. qhlshr(mgs) < 0.0 ) THEN
        chgtmp = scx(mgs,lhl) + dtp*(pschlmi(mgs) + pschlmd(mgs))
        IF ( sign(1.0,chgtmp) == sign(1.0,scx(mgs,lhl) ) ) THEN
          fschl2(mgs) = chgtmp/tmp
          pschlmd(mgs) =  pschlmd(mgs) + fschl2(mgs)*qhlshr(mgs)
        ENDIF
      ENDIF

      end do
      end if
      ENDIF


!
!  Rain
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscrwmi(mgs) =    &
     &    +fsccw(mgs)*qracw(mgs) +fsccw(mgs)*qrcnw(mgs)    &
     &  +(1-il5(mgs))*   &
     &   (-fschl(mgs)*qhlmlr(mgs) - fschw(mgs)*qhmlr(mgs)-fscsw(mgs)*qsmlr(mgs))*mixedphasefac   &
!     &    -fschl2(mgs)*qhlshr(mgs) - fschw2(mgs)*qhshr(mgs) - fscsw(mgs)*qsshr(mgs) ! -fscsw2(mgs)*qsshr(mgs)
     &  -(1-il5(mgs))*fscfw(mgs)*qfmlr(mgs)*mixedphasefac - fscfw2(mgs)*qfshr(mgs)  &
     &    -fschl2(mgs)*qhlshr(mgs) - fschw2(mgs)*qhshr(mgs) -fscsw2(mgs)*qsshr(mgs)
      pscrwmd(mgs) =     &
     &  + il5(mgs)*(-fscrw(mgs)*qiacr(mgs)-fscrw(mgs)*qrfrz(mgs))    &
     &  -fscrw(mgs)*qsacr(mgs)-fscrw(mgs)*qhacr(mgs) - fscrw(mgs)*qhlacr(mgs)  &
     &  -fscrw(mgs)*qfacr(mgs) &
!       + fscrw(mgs)*max(0.0,qrcev(mgs))  &
        + fscrw(mgs)*min(0.0,qrcev(mgs))
      end do
      end if

      pscpii(:) = 0.0 ! positive ion charge production here
     pscpid(:) = 0.0 
     pscnii(:) = 0.0 ! negative ion charge production here
     pscnid(:) = 0.0 
!
!
! Collisional Charging transfers (inductive/noninductive)
!
!
!  rain water
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscrwi(mgs) = 0.0
      pscrwd(mgs) = 0.0

      pscrwi(mgs) =  pscrwi(mgs)  &
     & + max(-schacr(mgs),0.0)    &
     & + max(-scfacr(mgs),0.0)    &
     & + max(-scsacr(mgs),0.0)    &
     & + max(scsshr(mgs),0.0)     & !Canosa & List 93
     & + max(schshr(mgs),0.0)       !Canosa & List 93

      pscrwd(mgs) =  pscrwd(mgs)  &
     & + min(-schacr(mgs),0.0)    &
     & + min(-scfacr(mgs),0.0)    &
     & + min(-scsacr(mgs),0.0)    &
     & + min(scsshr(mgs),0.0)     & !Canosa & List 93
     & + min(schshr(mgs),0.0)       !Canosa & List 93

      end do
      end if

!
!  cloud water
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccwi(mgs) = 0.0
      psccwd(mgs) = 0.0

      psccwi(mgs) =  psccwi(mgs)   &
     & + max(-schacw(mgs),0.0)   &
     & + max(-scfacw(mgs),0.0)  &
     & + max(-schlacw(mgs),0.0)   &
     & + max(-scsacw(mgs),0.0)  

      psccwd(mgs) =  psccwd(mgs)   &
     & + min(-schacw(mgs),0.0)   &
     & + min(-scfacw(mgs),0.0)  &
     & + min(-schlacw(mgs),0.0)   &
     & + min(-scsacw(mgs),0.0)   

!      end do
!      end if

       qcwtmp(mgs) = Max(0.0,qx(mgs,lc) +   &
     &   dtp*(pqcwi(mgs)+pqcwd(mgs)) )

       IF ( qcwtmp(mgs) .gt. qxmin(lc) ) then
       
      IF ( ipconc .lt. 2 ) THEN
       ccwtmp = cwccn
       ccwtmp = Max(1.0,qcwtmp(mgs)*rho0(mgs)/xmas(mgs,lc))
      ELSE
 
       ccwtmp = Max(0.001, cx(mgs,lc) + dtp*(pccwi(mgs)+pccwd(mgs)) )
      
      ENDIF ! ( ipconc .lt. 2 )

      
       cpqc = (scx(mgs,lc)   &
     &   +dtp*(psccwi(mgs)+psccwd(mgs)   &
     &   +psccwmi(mgs)+psccwmd(mgs)))/ccwtmp
      
         IF ( Abs(cpqc) .gt. scwppmx .and.   &
     &       ( psccwi(mgs)+psccwd(mgs) ) .ne. 0.0 .and.   &
     &    Sign(1.0,scx(mgs,lc)) .eq.   &
     &        Sign(1.0,( psccwi(mgs)+psccwd(mgs) )) ) THEN
          
           cpqc0 = Sign(scwppmx,cpqc)
          
           scfac = (cpqc0*ccwtmp - scx(mgs,lc) -   &
     &       dtp*(psccwmi(mgs)+psccwmd(mgs)))/   &
     &              ( dtp*(psccwi(mgs)+psccwd(mgs)) )
!           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: scwfac = ',scfac
             write(0,*) 'OUCH: scwfac = ',scfac
           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)
        
!        IF ( scfac .lt. 1.0 .and. scfac .gt. 0.0 ) THEN
!        write(iunit,*) 'LOOK: sccwfac = ',scfac,cpqc,scx(mgs,lc),ccwtmp,
!     : kgs(mgs),temg(mgs)
!        ENDIF
         
        ELSE
          scfac = 1.0
        ENDIF
      
      ELSE
       scfac = 0.0
      ENDIF ! ( qcwtmp .gt. qxmin(lc) ) 
      
      IF ( scfac .lt. 1.0 .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
       schacw(mgs) = scfac*schacw(mgs)
       schlacw(mgs) = scfac*schlacw(mgs)
       scsacw(mgs) = scfac*scsacw(mgs)
       scfacw(mgs) = scfac*scfacw(mgs)

      psccwi(mgs) = 0.0
      psccwd(mgs) = 0.0
      psccwi(mgs) =  psccwi(mgs)   &
     & + max(-schacw(mgs),0.0)   &
     & + max(-scfacw(mgs),0.0)  &
     & + max(-schlacw(mgs),0.0)   &
     & + max(-scsacw(mgs),0.0)
      psccwd(mgs) =  psccwd(mgs)   &
     & + min(-schacw(mgs),0.0)   &
     & + min(-scfacw(mgs),0.0)  &
     & + min(-schlacw(mgs),0.0)   &
     & + min(-scsacw(mgs),0.0)
      ENDIF
      end do
      
      
      end if

!
!  cloud ice
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccii(mgs) = 0.0
      psccid(mgs) = 0.0
      psccii(mgs) =  psccii(mgs)   &
     & + max(-schaci(mgs),0.0)   &
     & + max(-scfaci(mgs),0.0)  &
     & + max(-schlaci(mgs),0.0)   &
     & + max(-scsaci(mgs),0.0)
      psccid(mgs) =  psccid(mgs)   &
     & + min(-schaci(mgs),0.0)   &
     & + min(-scfaci(mgs),0.0)  &
     & + min(-schlaci(mgs),0.0)   &
     & + min(-scsaci(mgs),0.0)

      qcitmp = Max(0.0,qx(mgs,li) +   &
     &   dtp*(pqcii(mgs)+pqcid(mgs)) )

      IF ( qcitmp .gt. qxmin(li) ) then

      IF ( ipconc .le. 0 ) THEN
      
      ccitmp = cnina(mgs)
      IF ( cimn .gt. 1.0 ) THEN
        cx(mgs,li) = Max(cimn,cx(mgs,li))
      ENDIF
      IF ( cimx .gt. 1.0 ) THEN
        cx(mgs,li) = Min(cimx,cx(mgs,li))
      ENDIF
      IF ( itype1 .ge. 1 .or. itype2 .ge. 1 ) THEN 
!        ccitmp = Max(ccitmp,qcitmp*rho0(mgs)/cimasx)
        ccitmp = Min(ccitmp,qcitmp*rho0(mgs)/cimasn)
      ENDIF
      ccitmp = max(1.0,ccitmp)

      ELSE
 
       ccitmp = Max(0.001, cx(mgs,li) + dtp*(pccii(mgs)+pccid(mgs)) )
      
      ENDIF ! ( ipconc .le. 0 )
      
      
      cpci = (scx(mgs,li)   &
     &   +dtp*( psccii(mgs)+psccid(mgs)   &
     &         +psccimi(mgs)+psccimd(mgs)))/ccitmp
      
          IF ( Abs(cpci) .gt. scippmx .and.   &
     &       ( psccii(mgs)+psccid(mgs) ) .ne. 0.0 .and.   &
     &    Sign(1.0,scx(mgs,li)) .eq.   &
     &        Sign(1.0,( psccii(mgs)+psccid(mgs) )) ) THEN
          
           cpci0 = Sign(scippmx,cpci)
          
           scfac = (cpci0*ccitmp - scx(mgs,li) -    &
     &       dtp*(psccimi(mgs)+psccimd(mgs)))/   &
     &              ( dtp*(psccii(mgs)+psccid(mgs)) )
!           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: sccifac = ',scfac,   &
     & dtp*(psccimi(mgs)+psccimd(mgs)),cpci,scx(mgs,li),ccitmp
             write(0,*) 'OUCH: sccifac = ',scfac,   &
     & dtp*(psccimi(mgs)+psccimd(mgs)),cpci,scx(mgs,li),ccitmp,   &
     & kgs(mgs),temg(mgs)
           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)

!        IF ( scfac .lt. 1.0 ) THEN
!        write(iunit,*) 'LOOK: sccifac = ',scfac,cpci,scx(mgs,li),ccitmp,
!     : kgs(mgs),temg(mgs)
!        ENDIF
         
        ELSE
          scfac = 1.0
        ENDIF
      ELSE
       scfac = 0.0
      ENDIF ! ( qcitmp .gt. qxmin(li) ) 

      IF ( scfac .lt. 1.0 ) THEN


       schaci(mgs) = scfac*schaci(mgs)
       schlaci(mgs) = scfac*schlaci(mgs)
       scsaci(mgs) = scfac*scsaci(mgs)
       scfaci(mgs) = scfac*scfaci(mgs)

      psccii(mgs) = 0.0
      psccid(mgs) = 0.0
      psccii(mgs) =  psccii(mgs)   &
     & + max(-schaci(mgs),0.0)   &
     & + max(-scfaci(mgs),0.0)  &
     & + max(-schlaci(mgs),0.0)   &
     & + max(-scsaci(mgs),0.0)
      psccid(mgs) =  psccid(mgs)   &
     & + min(-schaci(mgs),0.0)   &
     & + min(-scfaci(mgs),0.0)  &
     & + min(-schlaci(mgs),0.0)   &
     & + min(-scsaci(mgs),0.0)
      
      ENDIF

      end do
      end if


!
!  snow
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscswi(mgs) = 0.0
      pscswd(mgs) = 0.0

      pscswi(mgs) =  pscswi(mgs)   &
     & + max(scsaci(mgs),0.0)   &
     & + max(scsacw(mgs),0.0)   &
     & + max(scsacr(mgs),0.0)   &
     & + max(-schacs(mgs),0.0)   &
     & + max(-schlacs(mgs),0.0)  &
     & + max(-scfacs(mgs),0.0)  &
     & + max(scsmlr(mgs),0.0)    &  !Drake
     & + max(scsdep(mgs),0.0)    &  !Dong&Hallett
     & - min(scsshr(mgs),0.0)       !Canosa & List 93
      pscswd(mgs) =  pscswd(mgs)   &
     & + min(scsaci(mgs),0.0)   &
     & + min(scsacw(mgs),0.0)   &
     & + min(scsacr(mgs),0.0)   &
     & + min(-schacs(mgs),0.0)   &
     & + min(-schlacs(mgs),0.0)  &
     & + min(-scfacs(mgs),0.0)   &
     & + min(scsmlr(mgs),0.0)    &   !Drake, but always 0.? (since scsmlr >= 0)
     & + min(scsdep(mgs),0.0)    &   !Dong&Hallett
     & - max(scsshr(mgs),0.0)        !Canosa & List 93
      end do
      end if

!
!  Graupel
!
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschwi(mgs) = 0.0
      pschwd(mgs) = 0.0

      pschwi(mgs) =  pschwi(mgs)   &
     & + max(schaci(mgs),0.0)   &
!     > + max(schacip(mgs),0.0)   &
!     > + max(schacir(mgs),0.0)   &
     & + max(schacs(mgs),0.0)   &
     & + max(schacw(mgs),0.0)   &
     & + max(schacr(mgs),0.0)   &
     & + max(schmlr(mgs),0.0)   & !Drake
     & - min(schshr(mgs),0.0)     !Canosa & List 93
      pschwd(mgs) =  pschwd(mgs)   &
     & + min(schaci(mgs),0.0)   &
!     > + min(schacip(mgs),0.0)   &
!     > + min(schacir(mgs),0.0)   &
     & + min(schacs(mgs),0.0)   &
     & + min(schacw(mgs),0.0)   &
     & + min(schacr(mgs),0.0)   &
     & + min(schmlr(mgs),0.0)   & !Drake, but always 0.? (since schmlr >= 0)
     & - max(schshr(mgs),0.0)     !Canosa & List 93
      
      
      end do
      end if

!
!  Frozen drops
!
      pscfwi(:) = 0.0
      pscfwd(:) = 0.0

      IF ( lf .gt. 1 ) THEN
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt

      pscfwi(mgs) =  pscfwi(mgs)   &
     & + max(scfaci(mgs),0.0)   &
     & + max(scfacs(mgs),0.0)   &
     & + max(scfacw(mgs),0.0)   &
     & + max(scfacr(mgs),0.0)

      pscfwd(mgs) =  pscfwd(mgs)   &
     & + min(scfaci(mgs),0.0)   &
     & + min(scfacs(mgs),0.0)   &
     & + min(scfacw(mgs),0.0)   &
     & + min(scfacr(mgs),0.0)

      IF ( lis > 1 ) THEN
        pscfwi(mgs) =  pscfwi(mgs) + max(scfacis(mgs),0.0)   
        pscfwd(mgs) =  pscfwd(mgs) + min(scfacis(mgs),0.0)
      ENDIF
      
!      IF ( ny == 2 .and. igs(mgs) == 21 .and. temg(mgs) < 273. .and. temg(mgs) > 233. ) THEN
!        write(0,*) 'k, pscfwi, pscfwd = ',kgs(mgs),pscfwi(mgs),pscfwd(mgs)
!        write(0,*) ' pschwi, pschwd = ',kgs(mgs),pschwi(mgs),pschwd(mgs)
!        write(0,*) 'chaci0,cfaci0 = ',chaci0(mgs),cfaci0(mgs)
!      ENDIF

      end do
      end if
      ENDIF

!
!  Hail
!
      pschli(:) = 0.0
      pschld(:) = 0.0

      IF ( lhl .gt. 1 ) THEN
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt

      pschli(mgs) =  pschli(mgs)   &
     & + max(schlaci(mgs),0.0)   &
     & + max(schlacs(mgs),0.0)   &
     & + max(schlacw(mgs),0.0)   &
     & + max(schlacr(mgs),0.0)

      pschld(mgs) =  pschld(mgs)   &
     & + min(schlaci(mgs),0.0)   &
     & + min(schlacs(mgs),0.0)   &
     & + min(schlacw(mgs),0.0)   &
     & + min(schlacr(mgs),0.0)


      end do
      end if
      ENDIF

!
! Ions
!
      pscpii(:) = 0.0 ! positive ion charge production here
      pscpid(:) = 0.0 
      pscnii(:) = 0.0 ! negative ion charge production here
      pscnid(:) = 0.0
      pscplii(:) = 0.0
      pscplid(:) = 0.0
      pscnlii(:) = 0.0
      pscnlid(:) = 0.0      
      
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
       IF (largeion) THEN
        ! Drop break-off from melting snow/graupel goes to large ions (Drake 1968)
        ! Evaporating drops create large ions (Takahashi 1979)
        ! Depositional growth through small ions (Dong & Hallett 1992)
        ! Condensing drops remove small ions, but handled in ion attachment
        pscnlii(mgs) =  pscnlii(mgs)  &
     &               + min(fscrw(mgs),0.0)*max(-qrcev(mgs),0.0)  &
     &               + min(fscsw(mgs),0.0)*max(-qscev(mgs),0.0)  &
     &               + min(-schmlr(mgs),0.0)    &
     &               + min(-scsmlr(mgs),0.0)
        pscplii(mgs) =  pscplii(mgs)  &
     &               + max(fscrw(mgs),0.0)*max(-qrcev(mgs),0.0)  &
     &               + max(fscsw(mgs),0.0)*max(-qscev(mgs),0.0)  &
     &               + max(-schmlr(mgs),0.0)    &
     &               + max(-scsmlr(mgs),0.0)
       ELSE
        pscnii(mgs) =  pscnii(mgs)  &
     &               - min(fscrw(mgs),0.0)*min(qrcev(mgs),0.0)  &
     &               - min(fscsw(mgs),0.0)*min(qscev(mgs),0.0)  &
     &               + min(-schmlr(mgs),0.0)    &
     &               + min(-scsmlr(mgs),0.0)
        pscpii(mgs) =  pscpii(mgs)  &
     &               - max(fscrw(mgs),0.0)*min(qrcev(mgs),0.0)  &
     &               - max(fscsw(mgs),0.0)*min(qscev(mgs),0.0)  &
     &               + max(-schmlr(mgs),0.0)    &
     &               + max(-scsmlr(mgs),0.0)
      ENDIF

        pscnid(mgs) =  pscnid(mgs)  &
     &               - max(-scsdep(mgs),0.0)
        pscpid(mgs) =  pscpid(mgs)  &
     &               - min(-scsdep(mgs),0.0)


      end do

      end if

!
!  continuity of space charge
!
      if ( ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      psctot(mgs) = psccwi(mgs) +psccwd(mgs) +   &
     &              psccii(mgs) +psccid(mgs) +   &
     &              pscrwi(mgs) +pscrwd(mgs) +   &
     &              pscswi(mgs) +pscswd(mgs) +   &
     &              pschwi(mgs) +pschwd(mgs) +       &
     &              pschli(mgs) +pschld(mgs) +       &
     &              pscfwi(mgs) +pscfwd(mgs) +       &
     &              pscpii(mgs) +pscpid(mgs) +       &
     &              pscnii(mgs) +pscnid(mgs) +       &
     &              pscplii(mgs) +pscplid(mgs) +       &
     &              pscnlii(mgs) +pscnlid(mgs) +       &
!
     &              psccwmi(mgs) +psccwmd(mgs) +   &
     &              psccimi(mgs) +psccimd(mgs) +   &
     &              pscrwmi(mgs) +pscrwmd(mgs) +   &
     &              pscswmi(mgs) +pscswmd(mgs) +   &
     &              pschwmi(mgs) +pschwmd(mgs) +     &
     &              pscfwmi(mgs) +pscfwmd(mgs) +     &
     &              pschlmi(mgs) +pschlmd(mgs)
      IF ( .not. ( psctot(mgs) .gt. -1. .and. psctot(mgs) .lt. 1.0 ) .or.    &
     &      Abs( psctot(mgs) ) .gt. esctot .or. Abs(pschwmi(mgs)) > 0.1 )  THEN

!        write(iunit,*) 'yikes! psctot = ',psctot(mgs),'and esctot =',esctot,'at',mgs

      IF ( .not. ( psctot(mgs) .gt. -1.0 .and. psctot(mgs) .lt. 1.0 ) .or.  &
            Abs( psctot(mgs) ) .gt. 1.0 ) THEN
        write(0,*)  psccwi(mgs), psccwd(mgs)
        write(0,*)  psccii(mgs), psccid(mgs)
        write(0,*)  pscrwi(mgs), pscrwd(mgs)
        write(0,*)  pscswi(mgs), pscswd(mgs)
        write(0,*)  pschwi(mgs), pschwd(mgs)
        write(0,*)  pscfwi(mgs), pscfwd(mgs)
        write(0,*)  pschli(mgs), pschld(mgs)
        write(0,*)  pscpii(mgs) ,pscpid(mgs)
        write(0,*)  pscnii(mgs) ,pscnid(mgs)
        write(0,*)
        write(0,*)  psccwmi(mgs), psccwmd(mgs)
        write(0,*)  psccimi(mgs), psccimd(mgs)
        write(0,*)  pscrwmi(mgs), pscrwmd(mgs)
        write(0,*)  pscswmi(mgs), pscswmd(mgs)
        write(0,*)  pschwmi(mgs), pschwmd(mgs)
        write(0,*)  pscfwmi(mgs), pscfwmd(mgs)
        write(0,*)  pschlmi(mgs), pschlmd(mgs)

      write(0,*) 'Problem with psctot,pschwmi: ',mgs,psctot(mgs),pschwmi(mgs)
      DO il = lc,lhab
       write(0,*) 'il,scx,qx = ',il,scx(mgs,il),qx(mgs,il)
      ENDDO

      write(0,*)     fsccw(mgs)*qracw(mgs) 
      write(0,*) fsccw(mgs)*qrcnw(mgs) 
      write(0,*) -fschl(mgs)*qhlmlr(mgs)
      write(0,*) -fschw(mgs)*qhmlr(mgs)
      write(0,*) -fscsw(mgs)*qsmlr(mgs)
      write(0,*)   -(1-il5(mgs))*fscfw(mgs)*qfmlr(mgs)*mixedphasefac 
      write(0,*)  - fscfw2(mgs)*qfshr(mgs)
      write(0,*)     -fschl2(mgs)*qhlshr(mgs) 
      write(0,*) - fschw2(mgs)*qhshr(mgs) 
      write(0,*) -fscsw2(mgs)*qsshr(mgs)

       call commasmpi_abort()
       ENDIF
      ENDIF

      IF ( Abs( psctot(mgs) ) .gt. esctot .and. ndebug > 0  )  THEN
      write(0,*) 'psctot: ',   psccwi(mgs) +psccwd(mgs) +   &
     &              psccii(mgs) +psccid(mgs) +   &
     &              pscrwi(mgs) +pscrwd(mgs) +   &
     &              pscswi(mgs) +pscswd(mgs) +   &
     &              pschwi(mgs) +pschwd(mgs) +       &
     &              pschli(mgs) +pschld(mgs) +       &
     &              pscfwi(mgs) +pscfwd(mgs) +       &
     &              pscpii(mgs) +pscpid(mgs) +       &
     &              pscnii(mgs) +pscnid(mgs) +       &
     &              pscplii(mgs) +pscplid(mgs) +       &
     &              pscnlii(mgs) +pscnlid(mgs) +       &
!
     &              psccwmi(mgs) +psccwmd(mgs) +   &
     &              psccimi(mgs) +psccimd(mgs) +   &
     &              pscrwmi(mgs) +pscrwmd(mgs) +   &
     &              pscswmi(mgs) +pscswmd(mgs) +   &
     &              pschwmi(mgs) +pschwmd(mgs) +     &
     &              pscfwmi(mgs) +pscfwmd(mgs) +     &
     &              pschlmi(mgs) +pschlmd(mgs)


      write(0,*) 'psctot1: ',   psccwi(mgs) +psccwd(mgs) +   &
     &              psccii(mgs) +psccid(mgs) +   &
     &              pscrwi(mgs) +pscrwd(mgs) +   &
     &              pscswi(mgs) +pscswd(mgs) +   &
     &              pschwi(mgs) +pschwd(mgs) +       &
     &              pschli(mgs) +pschld(mgs) +       &
     &              pscfwi(mgs) +pscfwd(mgs) +       &
     &              pscpii(mgs) +pscpid(mgs) +       &
     &              pscnii(mgs) +pscnid(mgs) +       &
     &              pscplii(mgs) +pscplid(mgs) +       &
     &              pscnlii(mgs) +pscnlid(mgs) 

       write(0,*) 'psctotx: ', psccwmi(mgs) +psccwmd(mgs) +   &
     &              psccimi(mgs) +psccimd(mgs) +   &
     &              pscrwmi(mgs) +pscrwmd(mgs) +   &
     &              pscswmi(mgs) +pscswmd(mgs) +   &
     &              pschwmi(mgs) +pschwmd(mgs) +     &
     &              pscfwmi(mgs) +pscfwmd(mgs) +     &
     &              pschlmi(mgs) +pschlmd(mgs)
        
        write(0,*)  psccwi(mgs), psccwd(mgs)
        write(0,*)  psccii(mgs), psccid(mgs)
        write(0,*)  pscrwi(mgs), pscrwd(mgs)
        write(0,*)  pscswi(mgs), pscswd(mgs)
        write(0,*)  pschwi(mgs), pschwd(mgs)
        write(0,*)  pscfwi(mgs), pscfwd(mgs)
        write(0,*)  pschli(mgs), pschld(mgs)
        write(0,*)  pscpii(mgs) ,pscpid(mgs)
        write(0,*)  pscnii(mgs) ,pscnid(mgs)
        write(0,*)
        write(0,*)  psccwmi(mgs), psccwmd(mgs)
        write(0,*)  psccimi(mgs), psccimd(mgs)
        write(0,*)  pscrwmi(mgs), pscrwmd(mgs)
        write(0,*)  pscswmi(mgs), pscswmd(mgs)
        write(0,*)  pschwmi(mgs), pschwmd(mgs)
        write(0,*)  pscfwmi(mgs), pscfwmd(mgs)
        write(0,*)  pschlmi(mgs), pschlmd(mgs)


      
         write(0,*) 'scfaci', scfaci(mgs),scfacs(mgs),scfacw(mgs),scfacr(mgs)

         write(0,*) 'pscfwmi'
         write(0,*)  +il5(mgs)*(ifiacrg*(1.0-ffrzh)*(fscci(mgs)*qracif(mgs)) )
         write(0,*)  +il5(mgs)*(ifiacrg*(1.0-ffrzh)*(fscrw(mgs)*qiacrf(mgs)) )
         write(0,*)   il5(mgs)*fscrw(mgs)*(1.0-ffrzh)*qrfrzf(mgs)
         write(0,*)   +fscrw(mgs)*qfacr(mgs),fsccw(mgs)*qfacw(mgs)   
         write(0,*)   +fscsw(mgs)*qfacs(mgs),fscci(mgs)*qfaci(mgs)   
         write(0,*) 'fscrw(mgs)*qrfrz(mgs) ',fscrw(mgs)*qrfrz(mgs),qrfrz(mgs),fscrw(mgs)*qrfrzs(mgs)
         write(0,*) 'hail: ', fscrw(mgs)*(1.0-ifrzg)*qrfrzf(mgs) 
          write(0,*) 'pscfwmd'
         write(0,*)   - fscfw(mgs)*qhlcnf(mgs)   
         write(0,*)   +(1-il5(mgs))*fscfw(mgs)*qfmlr(mgs)*mixedphasefac   
         write(0,*)   -fscfw(mgs)*qfmul1(mgs)
      ENDIF

      end do
      
      ENDIF

      do mgs = 1,ngscnt
        psctotmx = Max( psctotmx, psctot(mgs) )
        psctotmn = Min( psctotmn, psctot(mgs) )
      end do
! compute space charge arrays:

    if (ipelec.ge.1) then

     do mgs = 1,ngscnt

      DO il = lc,lhab
!      IF ( .not. (scx(mgs,il) > -1.e-6 .and. scx(mgs,il) < 1.e-6 ) .or. scx(mgs,il) > 1.e-9 ) THEN ! DEBUGTED
      IF ( .not. (scx(mgs,il) > -1.e-6 .and. scx(mgs,il) < 1.e-6 ) ) THEN 
        write(0,*) 'Problem1a with scx il = ',il,scx(mgs,il),mgs,igs(mgs),kgs(mgs),jyslab
        write(0,*) 'psccw: ',psccwi(mgs),psccwd(mgs),psccwmi(mgs),psccwmd(mgs)
        call commasmpi_abort()
      ENDIF
      ENDDO

     scx(mgs,lc) = scx(mgs,lc) +   &
    &   dtp*(psccwi(mgs)+psccwd(mgs)+psccwmi(mgs)+psccwmd(mgs))
     scx(mgs,li) = scx(mgs,li) +   &
    &   dtp*(psccii(mgs)+psccid(mgs)+psccimi(mgs)+psccimd(mgs))
     scx(mgs,lr) = scx(mgs,lr) +   &
    &   dtp*(pscrwi(mgs)+pscrwd(mgs)+pscrwmi(mgs)+pscrwmd(mgs))
     scx(mgs,ls) = scx(mgs,ls) +   &
    &   dtp*(pscswi(mgs)+pscswd(mgs)+pscswmi(mgs)+pscswmd(mgs))
     scx(mgs,lh) = scx(mgs,lh) +   &
    &   dtp*(pschwi(mgs)+pschwd(mgs)+pschwmi(mgs)+pschwmd(mgs))

     IF ( lf .gt. 1 ) THEN
     scx(mgs,lf) = scx(mgs,lf) +   &
    &   dtp*(pscfwi(mgs)+pscfwd(mgs)+pscfwmi(mgs)+pscfwmd(mgs))
     ENDIF

     IF ( lhl .gt. 1 ) THEN
     scx(mgs,lhl) = scx(mgs,lhl) +   &
    &   dtp*(pschli(mgs)+pschld(mgs)+pschlmi(mgs)+pschlmd(mgs))
     ENDIF

      cionp(mgs) = cionp(mgs) +  &
     &   dtp*eci*(pscpii(mgs)+pscpid(mgs)) ! #/m3 * C / C/m3. ! each ion has electron charge

      cionn(mgs) = cionn(mgs) -  &
     &   dtp*eci*(pscnii(mgs)+pscnid(mgs))

      if ( largeion ) then
      clionp(mgs) = clionp(mgs) +  &
     &   dtp*eci*(pscplii(mgs)+pscplid(mgs))
      clionn(mgs) = clionn(mgs) -  &
     &   dtp*eci*(pscnlii(mgs)+pscnlid(mgs))
      endif 

!
!  sum total charge
!
        sctot = 0.0
        DO il = lc,lhab
        sctot = sctot + scx(mgs,il)
        ENDDO
        
! C$PAR CRITICAL SECTION    
!! c$omp  critical
       sctot1p = sctot1p + Max(0.0d0,sctot)*dvmgs(mgs)
       sctot1n = sctot1n + Min(0.0d0,sctot)*dvmgs(mgs)
!! c$omp  end critical
! C$PAR END CRITICAL SECTION   
! cmic$ endguard

     end do
     end if ! (ipelec .ge. 1 )
    



!
!  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
!
!
      do mgs = 1,ngscnt

! do *before* updating qr
      IF ( iraintypes >= 1 ) THEN
        IF ( qx(mgs,lr) > qxmin(lr) ) THEN

!          write(0,*) 'update qxrain: pqrauto,shd,mlt = ',qx(mgs,lr),pqrauto(mgs),pqrshed(mgs),pqrmelt(mgs)
          tmpqxrain = 0.0
          DO il = 1,nraintypes

            IF ( lrain(il) == lrauto ) THEN
            
              qxrain(mgs,il) = qxrain(mgs,il)*(1.0 + dtp*pqrother(mgs)/qx(mgs,lr)) + dtp*(pqrauto(mgs) + pqrshed(mgs)*qxshedfrac(mgs,lr,1))
              tmpqxrain = tmpqxrain + qxrain(mgs,il) 
            ENDIF
            IF ( lrain(il) == lrshed ) THEN
              qxrain(mgs,il) = qxrain(mgs,il)*(1.0 + dtp*pqrother(mgs)/qx(mgs,lr)) + dtp*pqrshed(mgs)*qxshedfrac(mgs,lr,2)
              tmpqxrain = tmpqxrain + qxrain(mgs,il) 
            ENDIF
            IF ( lrain(il) == lrmelt ) THEN
              qxrain(mgs,il) = qxrain(mgs,il)*(1.0 + dtp*pqrother(mgs)/qx(mgs,lr)) + dtp*(pqrmelt(mgs) + pqrshed(mgs)*qxshedfrac(mgs,lr,3))
              tmpqxrain = tmpqxrain + qxrain(mgs,il) 
            ENDIF

          ENDDO
          
          ! Ensure that parts add up to sum
          tmpqr = qx(mgs,lr) + dtp*(pqrwi(mgs)+pqrwd(mgs))
          IF ( tmpqxrain > 0.d0 ) THEN
            frac = tmpqr/tmpqxrain
!            IF ( ndebug > 1 ) write(0,*) 'qxrain check: mgs, frac = ', mgs, frac, tmpqr, tmpqxrain
            DO il = 1,nraintypes
              qxrain(mgs,il) = frac*qxrain(mgs,il)
            ENDDO
          ENDIF
          
        ENDIF
      ENDIF

      qwvp(mgs) = qwvp(mgs) +        &
     &   dtp*(pqwvi(mgs)+pqwvd(mgs))
   !   qcwresv(mgs) = qx(mgs,lc) ! temporary save of old qc value
      qx(mgs,lc) = qx(mgs,lc) +   &
     &   dtp*(pqcwi(mgs)+pqcwd(mgs))
      qx(mgs,lr) = qx(mgs,lr) +   &
     &   dtp*(pqrwi(mgs)+pqrwd(mgs))
      qx(mgs,li) = qx(mgs,li) +   &
     &   dtp*(pqcii(mgs)+pqcid(mgs))
      qx(mgs,ls) = qx(mgs,ls) +   &
     &   dtp*(pqswi(mgs)+pqswd(mgs))
      qx(mgs,lh) = qx(mgs,lh) +    &
     &   dtp*(pqhwi(mgs)+pqhwd(mgs))
      IF ( lf > 1 ) THEN
      qx(mgs,lf) = qx(mgs,lf) +    &
     &   dtp*(pqfwi(mgs)+pqfwd(mgs))
!       IF ( qx(mgs,lf) < -1.e8 ) THEN
!         write(iunit,*) 'negative fd: ',qx(mgs,lf),pqfwi(mgs),pqfwd(mgs)
!         write(iunit,*) 'i,k = ',igs(mgs),kgs(mgs)
!       ENDIF
      ENDIF

      IF ( lhl .gt. 1 ) THEN
      qx(mgs,lhl) = qx(mgs,lhl) +    &
     &   dtp*(pqhli(mgs)+pqhld(mgs))
      ENDIF


      end do

! sum sources for particle volume

      IF ( ldovol ) THEN

      do mgs = 1,ngscnt

      IF ( lvol(ls) .gt. 1 ) THEN
      vx(mgs,ls) = vx(mgs,ls) +    &
     &   dtp*(pvswi(mgs)+pvswd(mgs))
      ENDIF

      IF ( lvol(lh) .gt. 1 ) THEN
      vx(mgs,lh) = vx(mgs,lh) +    &
     &   dtp*(pvhwi(mgs)+pvhwd(mgs))
!     >   rho0(mgs)*dtp*(pqhwi(mgs)+pqhwd(mgs))/xdn0(lh)
      ENDIF
      IF ( lf > 1 ) THEN
      IF ( lvol(lf) > 1 ) THEN
      vx(mgs,lf) = vx(mgs,lf) +    &
     &   dtp*(pvfwi(mgs)+pvfwd(mgs))
      ENDIF
      ENDIF

      IF ( lhl .gt. 1 ) THEN
      IF ( lvol(lhl) .gt. 1 ) THEN
      vx(mgs,lhl) = vx(mgs,lhl) +    &
     &   dtp*(pvhli(mgs)+pvhld(mgs))
!     >   rho0(mgs)*dtp*(pqhwi(mgs)+pqhwd(mgs))/xdn0(lh)
      ENDIF
      ENDIF

      ENDDO

      ENDIF  ! ldovol

!
!
!
! concentrations
!
      if ( ipconc .ge. 1  ) then
      do mgs = 1,ngscnt
      cx(mgs,li) = cx(mgs,li) +   &
     &   dtp*(pccii(mgs)+pccid(mgs)) 
      cina(mgs) = cina(mgs) + pccin(mgs)*dtp
      IF ( ipconc .ge. 2 ) THEN
      cx(mgs,lc) = cx(mgs,lc) +   &
     &   dtp*(pccwi(mgs)+pccwd(mgs))
      ENDIF
      IF ( ipconc .ge. 3 ) THEN
      cx(mgs,lr) = cx(mgs,lr) +   &
     &   dtp*(pcrwi(mgs)+pcrwd(mgs))
      ENDIF
      IF ( ipconc .ge. 4 ) THEN
      cx(mgs,ls) = cx(mgs,ls) +   &
     &   dtp*(pcswi(mgs)+pcswd(mgs))
      ENDIF
      IF ( ipconc .ge. 5 ) THEN
      cx(mgs,lh) = cx(mgs,lh) +    &
     &   dtp*(pchwi(mgs)+pchwd(mgs))
      IF ( lf > 1 ) THEN
      cx(mgs,lf) = cx(mgs,lf) +    &
     &   dtp*(pcfwi(mgs)+pcfwd(mgs))
      ENDIF
       IF ( lhl .gt. 1 ) THEN
        cx(mgs,lhl) = cx(mgs,lhl) +    &
     &     dtp*(pchli(mgs)+pchld(mgs))


        
!        IF ( (cx(mgs,lhl) > 5000. .or. cx(mgs,lhl) < 0. .or. dtp*(pchli(mgs)+pchld(mgs)) > 5000.) &
!              .and. qx(mgs,lhl) > qxmin(lhl) ) THEN
!          write(0,*) 'GS: cx,pci,pcd,time: ',cx(mgs,lhl), pchli(mgs),pchld(mgs),time_real
!          write(0,*) 'qx,pqi,pqd: ',qx(mgs,lhl), pqhli(mgs),pqhld(mgs)
!          write(0,*) 'alphl, xdia = ',alpha(mgs,lhl),xdia(mgs,lhl,3)
!          h1 = qx(mgs,lhl) - dtp*(pqhli(mgs)+pqhld(mgs))
!          h2 = cx(mgs,lhl) - dtp*(pchli(mgs)+pchld(mgs)) 
!          write(0,*) 'orig qx,cx = ',h1,h2
!          write(0,*) 'orig qx,cx = ',an(igs(mgs),jgs,kgs(mgs),lhl),an(igs(mgs),jgs,kgs(mgs),ln(lhl))
!          write(0,*) 'qx,cxsave = ',qhlsave(mgs),chlsave(mgs)
!          write(0,*) 'chlmlr, qmlr*cx/qx, save = ',chlmlr(mgs), qhlmlr(mgs)*h2/h1,chlmlrsave(mgs)
!          write(0,*) 'qhlmlrsave = ',qhlmlrsave(mgs)
!          write(0,*) 'qhlshr = ' ,  qhlshr(mgs)  
!          write(0,*) 'qhlmlr = ' , (1-il5(mgs))*qhlmlr(mgs)    
!          write(0,*)  qhlsbv(mgs) 
!          write(0,*) -qhlmul1(mgs), -qhcnhl(mgs)
!          write(0,*) 'chlmlr : ', (1-il5(mgs))*chlmlr(mgs) 
!          write(0,*) 'chlsbv,chcnhl : ', chlsbv(mgs), - chcnhl(mgs)
!
!        ENDIF
        
       ENDIF
      ENDIF
      IF ( ipconc .ge. 6 ) THEN
       IF ( lzr .gt. 1 ) THEN
       zx(mgs,lr) = zx(mgs,lr) +    &
     &   dtp*(pzrwi(mgs)+pzrwd(mgs))
       ENDIF
       IF ( lzs .gt. 1 ) THEN
       zx(mgs,ls) = zx(mgs,ls) +    &
     &   dtp*(pzswi(mgs)+pzswd(mgs))
       ENDIF
       IF ( lzh .gt. 1 ) THEN
       zx(mgs,lh) = zx(mgs,lh) +    &
     &   dtp*(pzhwi(mgs)+pzhwd(mgs))
       ENDIF
       IF ( lf > 1 ) THEN
        zx(mgs,lf) = zx(mgs,lf) +    &
     &   dtp*(pzfwi(mgs)+pzfwd(mgs))
       ENDIF
       IF ( lzhl .gt. 1 ) THEN
        zx(mgs,lhl) = zx(mgs,lhl) +    &
     &     dtp*(pzhli(mgs)+pzhld(mgs))
!      IF ( pchli(mgs) .ne. 0. .or. pchld(mgs) .ne. 0 ) THEN
!       write(0,*) 'dr: cx,pchli,pchld = ', cx(mgs,lhl),pchli(mgs),pchld(mgs), igs(mgs),kgs(mgs)
!      ENDIF
       ENDIF
      ENDIF
      end do
      end if

!
!
!
! start saturation adjustment
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 30a'
!      include 'sam.jms.satadj.sgi'
!
!
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!
!
!  set up temperature and vapor arrays
!
      do mgs = 1,ngscnt
      pqs(mgs) = (380.0)/(pres(mgs))
      theta(mgs) = thetap(mgs) + theta0(mgs)
      qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      end do
!
!  melting of cloud ice
!
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptimlw(mgs) = 0.0
      end do
!
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li)
      if( temg(mgs) .gt. tfr .and.   &
     &    qitmp(mgs) .gt. 0.0 ) then
      qx(mgs,lc) = qx(mgs,lc) + qitmp(mgs)
!      pfrz(mgs) = pfrz(mgs) - qitmp(mgs)*dtpinv
      ptem(mgs) =  ptem(mgs) +   &
     &  (1./pi0(mgs))*   &
     &  felfcp(mgs)*(- qitmp(mgs)*dtpinv)  
      IF ( eqtset > 2 ) THEN
        pipert(mgs) = pipert(mgs) - (felfpi(mgs)*qitmp(mgs))
      ENDIF
      pmlt(mgs) = pmlt(mgs) - qitmp(mgs)*dtpinv
      scx(mgs,lc) = scx(mgs,lc) + scx(mgs,li)
      thetap(mgs) = thetap(mgs) -   &
     &  fcc3(mgs)*qitmp(mgs)
      ptimlw(mgs) = -fcc3(mgs)*qitmp(mgs)*dtpinv
      cx(mgs,lc) = cx(mgs,lc) + cx(mgs,li)
      qx(mgs,li) = 0.0
      cx(mgs,li) = 0.0
      scx(mgs,li) = 0.0
      vx(mgs,li) = 0.0
      qitmp(mgs) = 0.0
      end if
      end do

!
!


!      do mgs = 1,ngscnt
!      qimlw(mgs) = (qcwtmp(mgs)-qx(mgs,lc))*dtpinv
!      end do
!
!  homogeneous freezing of cloud water
!
      IF ( warmonly < 0.8 ) THEN

      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptwfzi(mgs) = 0.0
      end do
!
      do mgs = 1,ngscnt

!      if( temg(mgs) .lt. tfrh ) THEN
!       write(0,*) 'GS: mgs,temp,qc,qi = ',mgs,temg(mgs),temcg(mgs),qx(mgs,lc),qx(mgs,li)
!      ENDIF

      ctmp = 0.0
      frac = 0.0
      qtmp = 0.0
      
!      if( ( temg(mgs) .lt. thnuc + 2. .or. (ibfc == 2 .and. temg(mgs) < thnuc + 10. ) ) .and.    &
!     &  qx(mgs,lc) .gt. qxmin(lc) .and. (ipconc < 2 .or. ibfc == 0 .or. ibfc == 2 )) then
! commented for test (12/01/2015):
!      if( temg(mgs) .lt. thnuc + 0. .and.    &
!     &  qx(mgs,lc) .gt. 0.0 .and. (ipconc < 2 .or. ibfc == 0 )) then
      if( ( ( temg(mgs) .lt. thnuc + 0.) .or. (temg(mgs) .lt. thnuc + 2. .and. ibfc >= 3) ) .and.    &
     &  qx(mgs,lc) .gt. 0.0 .and. (ipconc < 2 .or. ibfc == 0 .or. ibfc == 2)) then

      IF ( ibfc >= 3 ) THEN
        frac = Max( 0.25, Min( 1., ((thnuc + 2.) - temg(mgs) )/4.0 ) )
      ELSEIF ( ibfc /= 2 .or. ipconc < 2 ) THEN
        frac = Max( 0.25, Min( 1., ((thnuc + 1.) - temg(mgs) )/4.0 ) )
      ELSE
          volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b
                                               ! volt is given in cm**3, so factor of 1.e-6 to convert to m**3
         
         cwfrz(mgs) = cx(mgs,lc)*Exp(-volt/xv(mgs,lc)) ! number of droplets with volume greater than volt

         qtmp = cwfrz(mgs)*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
         frac = qtmp/qx(mgs,lc) ! reset number frozen to same fraction as mass. This makes 
                                                       ! sure that cwfrz and qwfrz are consistent and prevents 
                                                       ! spurious creation of ice crystals.
      
      ENDIF
      qtmp = frac*qx(mgs,lc)

      IF ( ibfc == 4 .and. lis >= 1 ) THEN
        qx(mgs,lis) = qx(mgs,lis) + qtmp
      ELSE
        qx(mgs,li) = qx(mgs,li) + qtmp ! qx(mgs,lc)
      ENDIF
      pfrz(mgs) = pfrz(mgs) + qtmp*dtpinv
      ptem(mgs) =  ptem(mgs) +   &
     &  (1./pi0(mgs))*   &
     &  felfcp(mgs)*(qtmp*dtpinv)  

      IF ( eqtset > 2 ) THEN
        pipert(mgs) = pipert(mgs) + felfpi(mgs)*qtmp
      ENDIF

!      IF ( lvol(li) .gt. 1 ) vx(mgs,li) = vx(mgs,li) + rho0(mgs)*qx(mgs,lc)/xdn0(li)
      IF ( lvol(li) .gt. 1 ) vx(mgs,li) = vx(mgs,li) + rho0(mgs)*qtmp/xdn0(li)

      IF ( ipconc .ge. 2 ) THEN
        ctmp = frac*cx(mgs,lc)
!        cx(mgs,li) = cx(mgs,li) + cx(mgs,lc)
        IF ( ibfc == 4 .and. lis >= 1 ) THEN
          cx(mgs,lis) = cx(mgs,lis) + ctmp
        ELSE
          cx(mgs,li) = cx(mgs,li) + ctmp
        ENDIF
      ELSE ! (ipconc .lt. 2 )
        ctmp = 0.0
        IF ( t9(igs(mgs),jgs,kgs(mgs)-1) .gt. qx(mgs,lc) ) THEN
           qtmp = frac*t9(igs(mgs),jgs,kgs(mgs)-1)  

!           cx(mgs,lc) = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
           ctmp = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
        ELSE
           cx(mgs,lc) = Max(0.0,wvel(mgs))*dtp*cwccn   &
     &      *gz(kgs(mgs))
          cx(mgs,lc) = cwccn
        ENDIF

       IF ( ipconc .ge. 1 ) cx(mgs,li) = Min(ccimx, cx(mgs,li) + cx(mgs,lc))
      ENDIF

      sctmp = frac*scx(mgs,lc)
!      scx(mgs,li) = scx(mgs,li) + scx(mgs,lc)
      scx(mgs,li) = scx(mgs,li) + sctmp
!      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qx(mgs,lc)
!      ptwfzi(mgs) = fcc3(mgs)*qx(mgs,lc)*dtpinv
!      qx(mgs,lc) = 0.0
!      cx(mgs,lc) = 0.0
!      scx(mgs,lc) = 0.0
      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qtmp
      ptwfzi(mgs) = fcc3(mgs)*qtmp*dtpinv
      qx(mgs,lc) = qx(mgs,lc) - qtmp
      cx(mgs,lc) = cx(mgs,lc) - ctmp
      scx(mgs,lc) = scx(mgs,lc) - sctmp
      end if
      end do

      ENDIF ! warmonly
!
!      do mgs = 1,ngscnt
!      qwfzi(mgs) = (qcwtmp(mgs)-qx(mgs,lc))*dtpinv   ! Not used?? (ERM)
!      end do
!
!  reset temporaries for cloud particles and vapor
!
      qcond(:) = 0.0
      
      IF ( ipconc .le. 1 .and.  lwsm6 ) THEN ! Explicit cloud condensation/evaporation (Rutledge and Hobbs 1983)
       DO mgs = 1,ngscnt

        qcwtmp(mgs) = qx(mgs,lc)
        theta(mgs) = thetap(mgs) + theta0(mgs)
        temgtmp = temg(mgs)
!        temg(mgs) = theta(mgs)*(p2(igs(mgs),jgs,kgs(mgs)) ) ! *pk(mgs) ! ( pres(mgs) / poo ) ** cap
!        temsav = temg(mgs)
!        thsave(mgs) = thetap(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
        temcg(mgs) = temg(mgs) - tfr
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )

        IF ( iqvsopt == 0 ) THEN
          qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        ELSEIF ( iqvsopt == 1 ) THEN
          qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
        ELSE
          tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
          qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
        ENDIF

        IF ( ( qvap(mgs) > qvs(mgs) .or. qx(mgs,lc) > qxmin(lc) ) .and. temg(mgs) > tfrh ) THEN
          tmp = (qvap(mgs) - qvs(mgs))/(1. + qvs(mgs)*felv(mgs)**2/(cp*rw*temg(mgs)**2) )
          qcond(mgs) = Min( Max( 0.0, tmp ), (qvap(mgs)-qvs(mgs)) )
          IF ( qx(mgs,lc) > qxmin(lc) .and. tmp < 0.0 ) THEN ! evaporation
            qcond(mgs) = Max( tmp, -qx(mgs,lc) )
          ENDIF
          qwvp(mgs) = qwvp(mgs) - qcond(mgs)
          qvap(mgs) = qvap(mgs) - qcond(mgs)
          qx(mgs,lc) = Max( 0.0, qx(mgs,lc) + qcond(mgs) )
          thetap(mgs) = thetap(mgs) + felvcp(mgs)*qcond(mgs)/(pi0(mgs))
          
        ENDIF
        
        ENDDO
      
      ENDIF
      
      
      IF ( ipconc .le. 1 .and. .not. lwsm6 ) THEN
!      IF ( ipconc .le. 1  ) THEN
      
      do mgs = 1,ngscnt
      qx(mgs,lv) = max( 0.0, qvap(mgs) )
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qx(mgs,li) = max( 0.0, qx(mgs,li) )
      qitmp(mgs) = qx(mgs,li)
      end do
!
!
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      qitmp(mgs) = qx(mgs,li)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temgtmp = temg(mgs)
      temg(mgs) = theta(mgs)*(pinit(kgs(mgs)) + p2(igs(mgs),jgs,kgs(mgs)) ) ! *pk(mgs) ! ( pres(mgs) / poo ) ** cap
      temsav = temg(mgs)
      thsave(mgs) = thetap(mgs)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      IF ( iqvsopt == 0 ) THEN
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
        qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )   &
     &  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /   &
     &   (qx(mgs,lc) + qitmp(mgs))
      end if
      end do
!
!  iterate  adjustment
!
      do itertd = 1,2
!
      do mgs = 1,ngscnt
!
!  calculate super-saturation
!
      qitmp(mgs) = qx(mgs,li)
      fcci(mgs) = 0.0
      fcip(mgs) = 0.0
      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qx(mgs,lv) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           !  subsaturated
        if( qx(mgs,lc) .gt. -dqwv(mgs) ) then  ! check if qc can make up all of the deficit
          dqcw(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                 !  otherwise make all qc available for evap
          dqcw(mgs) = -qx(mgs,lc)
          dqwv(mgs) = dqwv(mgs) + qx(mgs,lc)
        end if
!
        if( qitmp(mgs) .gt. -dqwv(mgs) ) then  ! check if qi can make up all the deficit
          dqci(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                  ! otherwise make all ice available for sublimation
          dqci(mgs) = -qitmp(mgs)
          dqwv(mgs) = dqwv(mgs) + qitmp(mgs)
        end if
!
       qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs) + dqci(mgs) )  ! add to perturbation vapor
!
! This next line removed 3/19/2003 thanks to Adam Houston,
!  who found the bug in the 3-ICE code
!      qwvp(mgs) = max(qwvp(mgs), 0.0)
      qitmp(mgs) = qx(mgs,li)
      IF ( qitmp(mgs) .ge. qxmin(li) ) THEN
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
      ELSE
        fcci(mgs) = 1.0
      ENDIF
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
      qx(mgs,li) = qx(mgs,li) + dqci(mgs) * fcci(mgs)
      thetap(mgs) = thetap(mgs) +   &
     &  1./pi0(mgs)*   &
     &  (felvcp(mgs)*dqcw(mgs) +felscp(mgs)*dqci(mgs))

      IF ( eqtset > 2 ) THEN
        pipert(mgs) = pipert(mgs)   &
     &  +(felspi(mgs)*dqci(mgs)    &
     &  +felvpi(mgs)*dqcw(mgs)) ! *dtp (remove dtp since dqxx are not rates)
      ENDIF

      end if  ! dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      IF ( dqwv(mgs) .ge. 0. ) THEN
      
!      write(iunit,*) 'satadj: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
!
        qitmp(mgs) = qx(mgs,li)
        fracl(mgs) = 1.0
        fraci(mgs) = 0.0
        if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
          fracl(mgs) = max(min(1.,(temg(mgs)-233.15)/(20.)),0.0)
          fraci(mgs) = 1.0-fracl(mgs)
        end if
        if ( temg(mgs) .le. thnuc ) then
           fraci(mgs) = 1.0
           fracl(mgs) = 0.0
         end if
        fraci(mgs) = 1.0-fracl(mgs)
!
       gamss = (felvcp(mgs)*fracl(mgs) + felscp(mgs)*fraci(mgs))   &
     &      / (pi0(mgs))
!
      IF ( temg(mgs) .lt. tfr ) then
        IF (qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) ) then
         dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbw)**2))
        END IF
        IF ( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
          dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv2(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbi)**2))
        END IF
        IF ( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
         cdw = caw*pi0(mgs)*tfrcbw/((temg(mgs)-cbw)**2)
         cdi = cai*pi0(mgs)*tfrcbi/((temg(mgs)-cbi)**2)
         denom1 = qx(mgs,lc) + qitmp(mgs)
         denom2 = 1.0 + gamss*   &
     &    (qx(mgs,lc)*qvs(mgs)*cdw + qitmp(mgs)*qis(mgs)*cdi) / denom1
         dqvcnd(mgs) =  dqwv(mgs) / denom2
        END IF 

      ENDIF  !  temg(mgs) .lt. tfr
!
      if ( temg(mgs) .ge. tfr ) then
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbw)**2))
      end if
!
      delqci1=qx(mgs,li)
!
      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
      ELSE
        fcci(mgs) = 1.0
      ENDIF
!
      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs)
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)
!
      thetap(mgs) = thetap(mgs) +   &
     &   (felvcp(mgs)*dqcw(mgs) + felscp(mgs)*dqci(mgs))   &
     & / (pi0(mgs))

      IF ( eqtset > 2 ) THEN
        pipert(mgs) = pipert(mgs) + (0   &
     &  +felspi(mgs)*dqci(mgs)    &
     &  +felvpi(mgs)*dqcw(mgs)) ! *dtp (remove dtp since dqxx are not rates)
      ENDIF

      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
!      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        qx(mgs,li) = qx(mgs,li) + dqci(mgs)*fcci(mgs)
        qitmp(mgs) = qx(mgs,li)
!      ENDIF
!
!      delqci(mgs) =  dqci(mgs)*fcci(mgs)
!
      END IF !  dqwv(mgs) .ge. 0.
      end do
!
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      qvap(mgs) = Max((qwvp(mgs) + qv0(mgs)), 0.0)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      IF ( iqvsopt == 0 ) THEN
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
        qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qitmp(mgs) = max( 0.0, qitmp(mgs) )
      qx(mgs,lv) = max( 0.0, qvap(mgs))
!      if ( temg(mgs) .lt. tfr ) then
!      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
!     >  qss(mgs) = qvs(mgs)
!c      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
!      if( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
!     >  qss(mgs) = qis(mgs)
!c      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
!      if( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
!     >  qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
!     > (qx(mgs,lc) + qitmp(mgs))
!      else
!      qss(mgs) = qvs(mgs)
!      end if
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )   &
     &  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /   &
     &   (qx(mgs,lc) + qitmp(mgs))
      end if
!      pceds(mgs) = (thetap(mgs) - thsave(mgs))*dtpinv
!      write(iunit,*) 'satadj2: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
      end do
!
!  end the saturation adjustment iteration loop
!
      end do

     ENDIF ! ( ipconc .le. 1 )

!
!  spread the growth owing to vapor diffusion onto the
!  ice crystal categories using the
!
!  END OF SATURATION ADJUSTMENT
!

      if (ndebug .gt. 0 ) write(0,*) 'conc 30b'
!
!
!  end of saturation adjustment

       IF ( allocated( microrates ) ) THEN
        DO mgs = 1,ngscnt
         microrates(igs(mgs),jy,kgs(mgs),1) = pfrz(mgs)*felf(mgs)*3600./(cp*pi0(mgs))
         microrates(igs(mgs),jy,kgs(mgs),2) = psub(mgs)*fels(mgs)*3600./(cp*pi0(mgs))
         microrates(igs(mgs),jy,kgs(mgs),3) = pvap(mgs)*felv(mgs)*3600./(cp*pi0(mgs))
         microrates(igs(mgs),jy,kgs(mgs),4) = qrcev(mgs)*felv(mgs)/(cp*pi0(mgs))
         microrates(igs(mgs),jy,kgs(mgs),5) = (qhmlr(mgs)+qhlmlr(mgs))*felf(mgs)/(cp*pi0(mgs))
!         microrates(igs(mgs),jy,kgs(mgs),6) = pcond2(mgs)
        ENDDO
       ENDIF

     IF ( ipelec >= 1 ) THEN
      DO il = lc,lhab
       DO mgs = 1,ngscnt
        psctot2 = psctot2 + scx(mgs,il)*dvmgs(mgs)
       ENDDO
      ENDDO



       IF ( allocated( elecrates ) ) THEN
        DO mgs = 1,ngscnt
         elecrates(igs(mgs),jy,kgs(mgs), 9) = ehw(mgs)*qx(mgs,lc)*rho0(mgs)*1.e3*(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*rarfac
         elecrates(igs(mgs),jy,kgs(mgs),10) = schaci(mgs)
         elecrates(igs(mgs),jy,kgs(mgs),11) = ehw(mgs)
         elecrates(igs(mgs),jy,kgs(mgs),11) = ehw(mgs)
         IF ( lhl > 1 ) THEN
         elecrates(igs(mgs),jy,kgs(mgs),12) = ehlw(mgs)*qx(mgs,lc)*rho0(mgs)*1.e3*(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))*rarfac
         ELSE
         elecrates(igs(mgs),jy,kgs(mgs),12) = ehw(mgs)*qx(mgs,lc)*rho0(mgs)*1.e3*(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*rarfac
         ENDIF
         elecrates(igs(mgs),jy,kgs(mgs),13) = schlaci(mgs)
         elecrates(igs(mgs),jy,kgs(mgs),14) = ehlw(mgs)
!         elecrates(igs(mgs),jy,kgs(mgs),4) = qrcev(mgs)*felv(mgs)/(cp*pi0(kgs(mgs)))
!         elecrates(igs(mgs),jy,kgs(mgs),5) = (qhmlr(mgs)+qhlmlr(mgs))*felf(mgs)/(cp*pi0(mgs))
        ENDDO
       ENDIF


     ENDIF
!
!
! !DIR$ IVDEP
      do mgs = 1,ngscnt
      t0(igs(mgs),jy,kgs(mgs)) =  temg(mgs)
      end do
!
! Load the save arrays
!
      IF ( numproc > 1 ) THEN
      DO mgs = 1,ngscnt
       dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
!       tq1(igs(mgs),jy,kgs(mgs),1) = qsmlr(mgs)
       tq1(igs(mgs),jy,kgs(mgs),1) = crfrzf(mgs) + il5(mgs)*ciacrf(mgs)
       tq1(igs(mgs),jy,kgs(mgs),2) = chcnsh(mgs) + chcnih(mgs)
       tq1(igs(mgs),jy,kgs(mgs),3) = schmlr(mgs) + scsmlr(mgs) + scsdep(mgs) + scsshr(mgs) + schshr(mgs) !"alternative" charging mechanisms

       IF ( iraintypes > 0 ) THEN
         IF ( qxrainold(mgs,0) > qxmin(lr) ) THEN
          DO il = 1,nraintypes
           rate2d(igs(mgs),jy,il) = rate2d(igs(mgs),jy,il) + &
                (qrfrzf(mgs) + qiacrf(mgs))*qxrainold(mgs,il)/qxrainold(mgs,0)*rho0(mgs)*dtp*dv
          ENDDO
         ENDIF
       ENDIF

       IF ( ipconc > 2 ) THEN
       thproc(kzbeg-1+kgs(mgs),1) = thproc(kzbeg-1+kgs(mgs),1) + crfrzf(mgs)*dtp*dv
       ELSE
       thproc(kzbeg-1+kgs(mgs),1) = thproc(kzbeg-1+kgs(mgs),1) + qrfrzf(mgs)*rho0(mgs)*dtp*dv
       ENDIF
       thproc(kzbeg-1+kgs(mgs),2) = thproc(kzbeg-1+kgs(mgs),2) + il5(mgs)*ciacrf(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),3) = thproc(kzbeg-1+kgs(mgs),3) + chcnsh(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),4) = thproc(kzbeg-1+kgs(mgs),4) + chcnih(mgs)*dtp*dv
       IF (  qhacw(mgs)+qhacr(mgs) > 0.0 .and. temg(mgs) < tfr ) THEN
       thproc(kzbeg-1+kgs(mgs),5) = thproc(kzbeg-1+kgs(mgs),5) + (qhacw(mgs)+qhacr(mgs)+qhshr(mgs))*rho0(mgs)*dtp*dv
       ENDIF
       thproc(kzbeg-1+kgs(mgs),6) = thproc(kzbeg-1+kgs(mgs),6) + qracw(mgs)*rho0(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),7) = thproc(kzbeg-1+kgs(mgs),7) + qrcnw(mgs)*rho0(mgs)*dtp*dv
       IF ( qhacw(mgs) > 0.0 .and. temg(mgs) < tfr ) THEN
       thproc(kzbeg-1+kgs(mgs),8) = thproc(kzbeg-1+kgs(mgs),8) + (vhacw(mgs)+vhacr(mgs)+vhshdr(mgs))*dtp*dv
!       thproc(kzbeg-1+kgs(mgs),8) = thproc(kzbeg-1+kgs(mgs),8) + qhacw(mgs)*rho0(mgs)/rimdn(mgs,lh)*dtp*dv
       ENDIF
       thproc(kzbeg-1+kgs(mgs),9) = thproc(kzbeg-1+kgs(mgs),9) + pi0(mgs)*ptem(mgs)*dtp*dv  ! latent heating -- mult by pi0 to get T rather than theta
       thproc(kzbeg-1+kgs(mgs),10) = thproc(kzbeg-1+kgs(mgs),10) +  &
     &                             ( chmul1(mgs) + chlmul1(mgs)  )*dtp*dv
       IF ( lf > 1 ) THEN
       thproc(kzbeg-1+kgs(mgs),11) = thproc(kzbeg-1+kgs(mgs),11) +  &
     &                             ( cfmul1(mgs) )*dtp*dv
       ELSE
       thproc(kzbeg-1+kgs(mgs),11) = thproc(kzbeg-1+kgs(mgs),11) +  &
     &                             ( csplinter(mgs) + csplinter2(mgs) )*dtp*dv
       ENDIF
       thproc(kzbeg-1+kgs(mgs),12) = thproc(kzbeg-1+kgs(mgs),12) + qrfrzf(mgs)*rho0(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),13) = thproc(kzbeg-1+kgs(mgs),13) + il5(mgs)*qiacrf(mgs)*rho0(mgs)*dtp*dv ! mass of rain freezing by ice crystal capture
       thproc(kzbeg-1+kgs(mgs),14) = thproc(kzbeg-1+kgs(mgs),14) + crcnw(mgs)*dtp*dv    ! rain drop prod. by autoconv.
       thproc(kzbeg-1+kgs(mgs),15) = thproc(kzbeg-1+kgs(mgs),15) + (pcrwi(mgs)-crcnw(mgs))*dtp*dv ! rain drop prod by melting/shedding (i.e., everything but autoconv.)
!       thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate
       thproc(kzbeg-1+kgs(mgs),19) = thproc(kzbeg-1+kgs(mgs),19) + pmlt(mgs)*rho0(mgs)*dv  ! melting rate
       thproc(kzbeg-1+kgs(mgs),20) = thproc(kzbeg-1+kgs(mgs),20) + pdep(mgs)*rho0(mgs)*dv  ! deposition rate
       thproc(kzbeg-1+kgs(mgs),21) = thproc(kzbeg-1+kgs(mgs),21) + (psub(mgs)-pdep(mgs))*rho0(mgs)*dv  ! sublimation rate
       thproc(kzbeg-1+kgs(mgs),22) = thproc(kzbeg-1+kgs(mgs),22) + (pfrz(mgs)-pmlt(mgs))*rho0(mgs)*dv  ! freezing rate

!       thproc(kzbeg-1+kgs(mgs),20) = thproc(kzbeg-1+kgs(mgs),20) + (1./pi0(mgs))*felfcp(mgs)*pvap(mgs)*rho0(mgs)*dv  ! deposition rate
!       thproc(kzbeg-1+kgs(mgs),21) = thproc(kzbeg-1+kgs(mgs),21) + (1./pi0(mgs))*felscp(mgs)*psub(mgs)*rho0(mgs)*dv  ! sublimation rate
!       thproc(kzbeg-1+kgs(mgs),22) = thproc(kzbeg-1+kgs(mgs),22) + (1./pi0(mgs))*felfcp(mgs)*pfrz(mgs)*rho0(mgs)*dv ! (pfrz(mgs)-pmlt(mgs))*rho0(mgs)*dv  ! freezing rate

       thproc(kzbeg-1+kgs(mgs),23) = thproc(kzbeg-1+kgs(mgs),23) + crfrzs(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),24) = thproc(kzbeg-1+kgs(mgs),24) + il5(mgs)*ciacrs(mgs)*dtp*dv

       thproc(kzbeg-1+kgs(mgs),25) = thproc(kzbeg-1+kgs(mgs),25) + qhmlr(mgs)*rho0(mgs)*dv  ! melting rate
       thproc(kzbeg-1+kgs(mgs),26) = thproc(kzbeg-1+kgs(mgs),26) + qhlmlr(mgs)*rho0(mgs)*dv  ! melting rate

       IF (  qhlacw(mgs)+qhlacr(mgs) > 0.0 .and. temg(mgs) < tfr ) THEN
        thproc(kzbeg-1+kgs(mgs),27) = thproc(kzbeg-1+kgs(mgs),27) + (qhlacw(mgs)+qhlacr(mgs)+qhlshr(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),28) = thproc(kzbeg-1+kgs(mgs),28) + (qhlacw(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),29) = thproc(kzbeg-1+kgs(mgs),29) + (qhlacr(mgs))*rho0(mgs)*dtp*dv
       ENDIF

       IF ( temg(mgs) < tfr ) THEN
        thproc(kzbeg-1+kgs(mgs),30) = thproc(kzbeg-1+kgs(mgs),30) + (qhacw(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),31) = thproc(kzbeg-1+kgs(mgs),31) + (qhacr(mgs))*rho0(mgs)*dtp*dv
       ENDIF

       thproc(kzbeg-1+kgs(mgs),32) = thproc(kzbeg-1+kgs(mgs),32) + qhlcnh(mgs)*rho0(mgs)*dtp*dv ! graupel mass conversion to hail

       IF ( ihrn > 0 ) THEN
       thproc(kzbeg-1+kgs(mgs),33) = thproc(kzbeg-1+kgs(mgs),33) + ciihr(mgs)*dtp*dv ! contact freezing of droplets
       ELSE
       IF ( qwctfz(mgs)*dtp >= qxmin(li) ) THEN
       thproc(kzbeg-1+kgs(mgs),33) = thproc(kzbeg-1+kgs(mgs),33) + cwctfz(mgs)*dtp*dv ! contact freezing of droplets
       ENDIF
       ENDIF
       thproc(kzbeg-1+kgs(mgs),34) = thproc(kzbeg-1+kgs(mgs),34) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate
       IF ( qiint(mgs)*dtp >= qxmin(li) ) THEN
       thproc(kzbeg-1+kgs(mgs),35) = thproc(kzbeg-1+kgs(mgs),35) + ciint(mgs)*dtp*dv ! primary ice initiation
       ENDIF
       IF ( lf > 1 ) THEN
          thproc(kzbeg-1+kgs(mgs),38) = thproc(kzbeg-1+kgs(mgs),38) + (vfacw(mgs)+vfacr(mgs)+vfshdr(mgs))*dtp*dv
       ELSE
          thproc(kzbeg-1+kgs(mgs),38) = thproc(kzbeg-1+kgs(mgs),38) + (vhacw(mgs)+vhacr(mgs)+vhshdr(mgs))*dtp*dv
        ENDIF
      IF ( lhl > 1 ) THEN
        thproc(kzbeg-1+kgs(mgs),36) = thproc(kzbeg-1+kgs(mgs),36) + chlcnhhl(mgs)*dtp*dv
        IF ( lf > 1 ) THEN 
          thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + chlcnfhl(mgs)*dtp*dv
        ELSEIF ( lnhf > 1 ) THEN 
             IF ( cx(mgs,lh) > 0.0 ) THEN
               frach = Min(1.0, chxf(mgs,lh)/cx(mgs,lh))
               frach = Max(0.0, frach)
               chlcnfhl(mgs) = frach*chlcnhhl(mgs)
             ELSE
               chlcnfhl(mgs) = 0.0
             ENDIF
        thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + chlcnfhl(mgs)*dtp*dv
        ENDIF
!        thproc(kzbeg-1+kgs(mgs),38) = thproc(kzbeg-1+kgs(mgs),38) + qhlcnh(mgs)*rho0(mgs)*dtp*dv
        IF ( lf > 1 ) thproc(kzbeg-1+kgs(mgs),39) = thproc(kzbeg-1+kgs(mgs),39) + qhlcnf(mgs)*rho0(mgs)*dtp*dv

        IF ( numproc >= 41 ) THEN
          thproc(kzbeg-1+kgs(mgs),40) = thproc(kzbeg-1+kgs(mgs),40) + pchld(mgs)*dtp*dv
          thproc(kzbeg-1+kgs(mgs),41) = thproc(kzbeg-1+kgs(mgs),41) + chlfmlr(mgs)*dtp*dv
        ENDIF

       ELSE
         IF ( lf > 1 ) THEN 
           thproc(kzbeg-1+kgs(mgs),36) = thproc(kzbeg-1+kgs(mgs),36) + (pcfwi(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),39) = thproc(kzbeg-1+kgs(mgs),39) + (pcfwd(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + (cfmlr(mgs))*dtp*dv
         ELSE
           thproc(kzbeg-1+kgs(mgs),36) = thproc(kzbeg-1+kgs(mgs),36) + (pchwi(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),39) = thproc(kzbeg-1+kgs(mgs),39) + (pchwd(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + (chmlr(mgs))*dtp*dv
         ENDIF
       ENDIF

!        IF ( numproc >= 44 ) THEN
          thproc(kzbeg-1+kgs(mgs),42) = thproc(kzbeg-1+kgs(mgs),42) - qhshr(mgs)*dtp*dv
          thproc(kzbeg-1+kgs(mgs),43) = thproc(kzbeg-1+kgs(mgs),43) - qfshr(mgs)*dtp*dv
          thproc(kzbeg-1+kgs(mgs),44) = thproc(kzbeg-1+kgs(mgs),44) - qhlshr(mgs)*dtp*dv
!        ENDIF
       
       IF ( iraintypes >= 1 ) THEN
         IF ( qxrainold(mgs,0) > qxmin(lr) ) THEN
          DO il = 1,nraintypes
          thproc(kzbeg-1+kgs(mgs),numproc-nraintypes+il) = thproc(kzbeg-1+kgs(mgs),numproc-nraintypes+il) + &
                (crfrzf(mgs) + il5(mgs)*ciacrf(mgs))*qxrainold(mgs,il)/qxrainold(mgs,0)*dtp*dv
          ENDDO
         ENDIF


       ENDIF
       
!       thproc(kzbeg-1+kgs(mgs),35) = thproc(kzbeg-1+kgs(mgs),35) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate


!      ptem(mgs) =    &
!     &  (1./pi0(mgs))*   &
!     &  (felfcp(mgs)*pfrz(mgs)   &
!     &  +felscp(mgs)*psub(mgs)    &
!     &  +felvcp(mgs)*pvap(mgs))

      ENDDO
      ENDIF


 1234 Continue

      IF ( io_flag ) THEN
        DO mgs = 1,ngscnt

          IF ( lcwnu >= 1 )  axtra(igs(mgs),jy,kgs(mgs),lcwnu ) = alpha(mgs,lc)
          IF ( lmvr  >= 1 .and. qx(mgs,lr) > 100.*qxmin(lr) ) axtra(igs(mgs),jy,kgs(mgs),lmvr ) = vtxbar(mgs,lr,1)
          IF ( lmvh  >= 1 .and. qx(mgs,lh) > 100.*qxmin(lh) ) axtra(igs(mgs),jy,kgs(mgs),lmvh ) = vtxbar(mgs,lh,1)
          IF ( lf > 1 ) THEN
            IF ( lmvf  >= 1 .and. qx(mgs,lf) > 100.*qxmin(lh) ) axtra(igs(mgs),jy,kgs(mgs),lmvf ) = vtxbar(mgs,lf,1)
          ENDIF
          IF ( lmvhl >= 1 .and. lhl > 1 .and. lzvhl == 0 ) THEN
            IF ( qx(mgs,lhl) > 100.*qxmin(lhl) ) axtra(igs(mgs),jy,kgs(mgs),lmvhl) = vtxbar(mgs,lhl,1)
          ENDIF
          IF ( lmvs  >= 1 .and. ( qx(mgs,ls) > 10.*qxmin(ls) ) ) axtra(igs(mgs),jy,kgs(mgs),lmvs ) = vtxbar(mgs,ls,1)
          IF ( lmvi  >= 1 .and. ( qx(mgs,li) > 10.*qxmin(li) ) ) axtra(igs(mgs),jy,kgs(mgs),lmvi ) = vtxbar(mgs,li,1)
          IF ( idfw > 1 ) THEN 
            axtra(igs(mgs),jy,kgs(mgs),idfw  ) =  1000.*xdia(mgs,lf,3)
          ENDIF
          IF ( idhw > 1 ) axtra(igs(mgs),jy,kgs(mgs),idhw  ) =  1000.*xdia(mgs,lh,3)
           !hack   IF ( lmnvhl > 1 )  axtra(igs(mgs),jy,kgs(mgs),lmnvhl) = qcwresv(mgs)
           !hack   IF ( lzvhl > 1 )   axtra(igs(mgs),jy,kgs(mgs),lzvhl) = qhacw(mgs)
          IF ( lhl > 1 ) THEN
            IF ( idhl > 1 ) axtra(igs(mgs),jy,kgs(mgs),idhl  ) =  1000.*xdia(mgs,lhl,3)
            IF ( idnhl > 1 ) axtra(igs(mgs),jy,kgs(mgs),idnhl  ) =  1000.*xdia(mgs,lhl,1)
            IF ( idmhl > 1 ) axtra(igs(mgs),jy,kgs(mgs),idmhl  ) =  1000.*(4.0 + alpha(mgs,lhl))*xdia(mgs,lhl,1)
            IF ( ialphahl > 1 ) axtra(igs(mgs),jy,kgs(mgs),ialphahl  ) =  alpha(mgs,lhl)
            IF ( mixedphase .and. qx(mgs,lhl) > qxmin(lhl) )THEN
        !      IF ( lmnvhl > 1 )  axtra(igs(mgs),jy,kgs(mgs),lmnvhl) = fhlw(mgs)
        !      IF ( lzvhl > 1 )   axtra(igs(mgs),jy,kgs(mgs),lzvhl) = qhlshr(mgs)
            ENDIF
          ELSE
            IF ( idhl > 1 .and. idhw == 0 ) axtra(igs(mgs),jy,kgs(mgs),idhl  ) =  1000.*xdia(mgs,lh,3)
            IF ( idnhl > 1 ) axtra(igs(mgs),jy,kgs(mgs),idnhl  ) =  1000.*xdia(mgs,lh,1)
            IF ( idmhl > 1 ) axtra(igs(mgs),jy,kgs(mgs),idmhl  ) =  1000.*(4.0 + alpha(mgs,lh))*xdia(mgs,lh,1)
            IF ( ialphahl > 1 ) axtra(igs(mgs),jy,kgs(mgs),ialphahl  ) =  alpha(mgs,lh)
          ENDIF

          IF ( ialphah > 1 .and. lh > 1 ) axtra(igs(mgs),jy,kgs(mgs),ialphah  ) =  alpha(mgs,lh)
          IF ( ialphaf > 1 .and. lf > 1 ) axtra(igs(mgs),jy,kgs(mgs),ialphaf  ) =  alpha(mgs,lf)
          IF ( ialphar > 1 .and. lr > 1 ) axtra(igs(mgs),jy,kgs(mgs),ialphar  ) =  alpha(mgs,lr)

!          IF ( lmvs  >= 1 .and. ( qx(mgs,ls) > qxmin(ls) .or. qx(mgs,li) > qxmin(li) ) )             &
!              axtra(igs(mgs),jy,kgs(mgs),lmvs ) = (qx(mgs,ls)*vtxbar(mgs,ls,1) +  &
!     &                    qx(mgs,li)*vtxbar(mgs,li,1))/(qx(mgs,ls) + qx(mgs,li))
          IF ( ld0 >= 1 .and. idiagnosecnu > 0 ) THEN
               axtra(igs(mgs),jy,kgs(mgs),ld0  ) = alpha(mgs,lc)
!          ELSEIF ( ld0 >= 1 .and.  dmrauto <= 1 ) THEN
!               axtra(igs(mgs),jy,kgs(mgs),ld0  ) = 1.2*xl2p(mgs)*1.e3
          ELSEIF ( ld0   >= 1 .and. qx(mgs,lr) > 100.*qxmin(lr) ) THEN
            IF ( qx(mgs,lr) > 1.e-8 ) THEN
              IF ( imurain == 1 ) THEN
              axtra(igs(mgs),jy,kgs(mgs),ld0  ) =  1000.*(3.67+alpha(mgs,lr))*xdia(mgs,lr,1) ! units of mm (Ulbrich 1983, JCAM)
              ELSE ! imurain == 3, 
              axtra(igs(mgs),jy,kgs(mgs),ld0  ) =  1000.*(1.678+alpha(mgs,lr))**(1./3.)*xdia(mgs,lr,1) ! units of mm (using method of Ulbrich 1983. See ventillation_stuff.nb)
! Below is the maximum mass diameter, which turns out to be very very close to the median volume diameter!
!              axtra(igs(mgs),jy,kgs(mgs),ld0  ) =  1000.*(5./3. + alpha(mgs,lr))**(1./3.)*xdia(mgs,lr,1) ! units of mm. This is maximum mass diam (see ventillation_stuff.nb)
              ENDIF
            ELSE
              axtra(igs(mgs),jy,kgs(mgs),ld0  ) = 0.0
            ENDIF
          ENDIF
!          IF ( ld0   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),ld0  ) =  xdia(mgs,lr,3)
          IF ( lqaut   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqaut  ) = qrcnw(mgs)
          IF ( lqacc   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqacc  ) = qracw(mgs)
          IF ( lqrevap >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqrevap) = qrcev(mgs) + qscev(mgs) + qhcev(mgs) + qhlcev(mgs)
          IF ( lqdep   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqdep  ) = pdep(mgs)
!          IF ( ny <= 2 ) THEN
!           IF ( lqmelt  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqmelt ) = qsmlr(mgs)
!          ELSE
           IF ( lqmelt  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqmelt ) = pmlt(mgs)
!          ENDIF
!          IF ( lqaut   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqaut  ) = pqlwfd(mgs)  !qfshr(mgs)
!          IF ( lqacc   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqacc  ) = qfmlr(mgs)
          IF ( lqrevap >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqrevap) = axtra(igs(mgs),jy,kgs(mgs),lqrevap) + qfcev(mgs)
!          IF ( lqdep   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqdep  ) = qx(mgs,lf)-qxw(mgs,lf)
!          IF ( lqmelt  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqmelt ) = pqlwfi(mgs) ! qfcev(mgs)
!          IF ( lqcevap >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqcevap) = 
!          IF ( lqcond  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqcond ) = 
!          IF ( ny <= 2 ) THEN
!            IF ( lqsub   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqsub  ) = qssbv(mgs)
!          ELSE
            IF ( lqsub   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqsub  ) = psub(mgs)-pdep(mgs)
!          ENDIF

!          IF ( lqrevap >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqrevap) = ptem(mgs)
!          IF ( lqdep   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqdep  ) = psub(mgs)
!          IF ( lqmelt  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqmelt ) = pfrz(mgs)
!          IF ( lqsub   >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lqsub  ) = pvap(mgs)
          IF ( lchlcnh >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lchlcnh ) = chlcnh(mgs)
!          IF ( lchlcnh >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lchlcnh ) = &
!                       axtra(igs(mgs),jy,kgs(mgs),lchlcnh )  + chlcnh(mgs)*dtp/thistory ! chlcnhhl(mgs)*dtp/thistory
          IF ( ldhlcnh >= 1 .and. dg0(mgs) > 0. ) axtra(igs(mgs),jy,kgs(mgs),ldhlcnh ) = dg0(mgs)*1000.

          IF ( lswagg  >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lswagg  ) = csacs(mgs)
          IF ( lswdn  >= 1 )  axtra(igs(mgs),jy,kgs(mgs),lswdn  ) = xdn(mgs,ls)
          IF ( lhwdn  >= 1  .and. lh > 1 )  THEN
             IF ( qx(mgs,lh) > qxmin(lh) ) axtra(igs(mgs),jy,kgs(mgs),lhwdn  ) = xdn(mgs,lh)
!             IF ( qx(mgs,lf) > qxmin(lf) ) axtra(igs(mgs),jy,kgs(mgs),lhwdn  ) = ffw(mgs)
!             IF ( qx(mgs,lh) > qxmin(lh) ) axtra(igs(mgs),jy,kgs(mgs),lhwdn  ) = fhw(mgs)
          ENDIF
          IF ( lhldn  >= 1  .and. lhl > 1 )  THEN
             IF ( qx(mgs,lhl) > qxmin(lhl) ) axtra(igs(mgs),jy,kgs(mgs),lhldn  ) = xdntmp(mgs,lhl)
          ENDIF
          IF ( lfwdn  >= 1  .and. lf > 1 )  THEN
             IF ( qx(mgs,lf) > qxmin(lf) ) axtra(igs(mgs),jy,kgs(mgs),lfwdn  ) = xdn(mgs,lf)
             IF ( lchlcnf >= 1 ) axtra(igs(mgs),jy,kgs(mgs),lchlcnf ) = &
                     axtra(igs(mgs),jy,kgs(mgs),lchlcnf ) + chlcnfhl(mgs)*dtp/thistory
            IF ( ldhlcnf >= 1 .and. df0(mgs) > 0. ) axtra(igs(mgs),jy,kgs(mgs),ldhlcnf ) = df0(mgs)*1000.
          ENDIF
          
          IF ( lswdia >= 1 )  THEN
!          write(0,*) 'write to lswdia,nxtra ',lswdia,nxtra,xdia(mgs,ls,1)
            axtra(igs(mgs),jy,kgs(mgs),lswdia ) = xdia(mgs,ls,1)*1.e3
!          write(0,*) 'wrote to lswdia,nxtra ',mgs
          ENDIF
          IF ( lswmass >= 1 )  THEN
!          write(0,*) 'write to lswdia,nxtra ',lswdia,nxtra,xdia(mgs,ls,1)
            axtra(igs(mgs),jy,kgs(mgs),lswmass ) = xmas(mgs,ls)*1.e6 ! units of mg
!          write(0,*) 'wrote to lswdia,nxtra ',mgs
          ENDIF
          
        ENDDO
      ENDIF



      IF ( ipelec >= 1 ) THEN

      DO mgs = 1,ngscnt
        t1(igs(mgs),jy,kgs(mgs)) = scsaci(mgs)
        t2(igs(mgs),jy,kgs(mgs)) = schacs(mgs)
        t3(igs(mgs),jy,kgs(mgs)) = schaci(mgs) 
        IF ( lf > 1 ) THEN
          t3(igs(mgs),jy,kgs(mgs)) =  t3(igs(mgs),jy,kgs(mgs)) + scfaci(mgs) 
          t2(igs(mgs),jy,kgs(mgs)) = t2(igs(mgs),jy,kgs(mgs)) + scfacs(mgs)
        ENDIF
        t4(igs(mgs),jy,kgs(mgs)) = scsacw(mgs)
        t5(igs(mgs),jy,kgs(mgs)) = schacw(mgs) + schlacw(mgs)
        t6(igs(mgs),jy,kgs(mgs)) = schlacs(mgs)
        t8(igs(mgs),jy,kgs(mgs)) = schlaci(mgs)

      ENDDO

      ENDIF


      if (ndebug .gt. 0 ) write(0,*) 'gs 11'

      do mgs = 1,ngscnt
!
      an(igs(mgs),jy,kgs(mgs),lt) =    &
     &  theta0(mgs) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) = qwvp(mgs) + qv0(mgs) !

      IF ( eqtset > 2 ) THEN
        p2(igs(mgs),jy,kgs(mgs)) = pipert(mgs)
      ENDIF
!
      
      DO il = lc,lhab
        IF ( ido(il) .eq. 1 ) THEN
        IF ( lf > 1 .and. il == lf ) THEN 
           lfsave(mgs,1) = an(igs(mgs),jy,kgs(mgs),il)
           lfsave(mgs,2) = qx(mgs,il)
        ENDIF
         an(igs(mgs),jy,kgs(mgs),il) = qx(mgs,il) +   &
     &     min( an(igs(mgs),jy,kgs(mgs),il), 0.0 )
         qx(mgs,il) = an(igs(mgs),jy,kgs(mgs),il)
        ENDIF
      ENDDO

      IF ( lcina > 1 ) THEN
        an(igs(mgs),jy,kgs(mgs),lcina) = cina(mgs)
      ENDIF

! do *before* updating qr
      IF ( iraintypes >= 1 ) THEN

          DO il = 1,nraintypes
            an(igs(mgs),jy,kgs(mgs),lrain(il)) = Max( 0.0, qxrain(mgs,il) )
          ENDDO
          
      ENDIF




!
!  6th moments
!

      IF ( ipconc .ge. 6 ) THEN
       DO il = lr,lhab
        IF ( lz(il) .gt. 1 ) THEN
        IF ( lf > 1 .and. il == lf ) THEN 
           lfsave(mgs,3) = an(igs(mgs),jy,kgs(mgs),lz(il))
           lfsave(mgs,4) = zx(mgs,il)
        ENDIF

         an(igs(mgs),jy,kgs(mgs),lz(il)) = zx(mgs,il) +   &
     &     min( an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0 )
         zx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lz(il))
         
        ENDIF
       ENDDO
       
      ENDIF
!
      end do
!

      if ( ipconc .ge. 1 ) then
      DO il = lc,lhab !{

!        write(0,*) 'limiter loop: il,ipc,lz: ',il,ipc(il),lz(il),ipconc

       IF ( ipconc .ge. ipc(il) .and. ido(il) > 0 ) THEN ! {

         IF (  ipconc .ge. 4 .and. ipc(il) .ge. 1 ) THEN ! {

!            write(0,*) 'MY limiter: il,ipc,lz: ',il,ipc(il),lz(il),lr,lzr
!            STOP

          IF ( lz(il) <= 1 .or. ioldlimiter == 1 ) THEN ! { { is a two-moment category so dont worry about reflectivity
          

           DO mgs = 1,ngscnt
            IF ( qx(mgs,il) .le. 0.0 ) THEN
              cx(mgs,il) = 0.0
            ELSE !{
              IF ( cx(mgs,il) .gt. cxmin .and. qx(mgs,il) > qxmin(il) ) THEN !{ only do this if mass is sufficient
!              xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
!              xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(cxmin,cx(mgs,il)))
                xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*cx(mgs,il))
              
!              IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!               write(0,*) 'dr: xv,cx,qx,xdn,ln = ',xv(mgs,il),cx(mgs,il),qx(mgs,il),xdn(mgs,il),ln(il)
!              ENDIF

               ! 8/26/2015 erm: apply imaxdiaopt for 2-moment also
               IF ( imaxdiaopt == 1 .or. il == lc .or. il == li .or. (il == lr .and. imurain == 3) .or. &
     &              (il == ls .and. imusnow == 3 ) .or. ( il >= lh .and. lh > 0 ) ) THEN
!              IF ( imaxdiaopt == 1 .or. (il == lr .and. imurain == 3) .or. .not. (il == lr .and. imurain == 1) ) THEN
                 xvbarmax = xvmx(il)
               ELSEIF ( imaxdiaopt == 2 ) THEN ! test against maximum mass diameter
                 xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
               ELSEIF ( imaxdiaopt == 3 ) THEN ! test against mass-weighted diameter
                 xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
               ELSE
                 xvbarmax = xvmx(il)
               ENDIF

               tmp = 1.0
               IF ( il == ls ) THEN
                 xvbarmax = xvbarmax*Max(1.,100./Min(100.,xdn(mgs,ls)))
               ENDIF
               
               IF ( xv(mgs,il) .lt. xvmn(il) .or. xv(mgs,il) .gt. xvbarmax ) THEN
                xv(mgs,il) = Min( xvbarmax, xv(mgs,il) )
                xv(mgs,il) = Max( xvmn(il), xv(mgs,il) )
                cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xv(mgs,il)*xdn(mgs,il))
               ENDIF
              
             ENDIF !}

!              IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!               write(0,*) 'dr: xv,cx,= ',xv(mgs,il),cx(mgs,il)
!              ENDIF

            ENDIF !}
           ENDDO ! mgs
          
          ELSE ! } { is three-moment, so have to adjust Z if size is too large
           IF ( il == lr .and. imurain == 3 ) THEN ! { { RAIN

!          rdmx = 
!          rdmn = 

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
           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
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
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

! check for artificial breakup (rain larger than allowed max size)
        IF (  xv(mgs,il) .gt. xvmx(il) .or. (ioldlimiter == 2 .and. xv(mgs,il) .gt. xvmx(il)/8.) ) THEN
          tmp = cx(mgs,il)
!            write(0,*) 'MY limiter: xv: ',xv(mgs,il), xv(mgs,il)/(xvmx(il)/8.)
!            STOP
          IF ( ioldlimiter == 2 ) THEN ! MY-style active breakup
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          ELSE ! simple cutoff 
            xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          ENDIF
            !xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            !cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          
          
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
           

           
           ENDIF
          ENDIF
          
          ENDIF
          
          ENDDO
!        CALL cld_cpu('Z-MOMENT-1r')  
           
           
           ELSEIF ( il == lh .or. il == lhl .or. il == lf .or. (il == lr .and. imurain == 1 )) THEN ! } { Rain, GRAUPEL OR HAIL

        
        
        DO mgs = 1,ngscnt

        IF ( lf > 1 .and. il == lf ) THEN 
           lfsave(mgs,5) = an(igs(mgs),jy,kgs(mgs),ln(il))
           lfsave(mgs,6) = cx(mgs,il)
        ENDIF
        
        IF ( il == lhl .and. lnhlf > 1 ) THEN
          IF ( cx(mgs,lhl) > cxmin ) THEN
            frac = chxf(mgs,lhl)/cx(mgs,lhl)
          ELSE
            frac = 0.0
          ENDIF
        ENDIF

        IF ( il == lh .and. lnhf > 1 ) THEN
          IF ( cx(mgs,lh) > cxmin ) THEN
            frach = chxf(mgs,lh)/cx(mgs,lh)
          ELSE
            frach = 0.0
          ENDIF
        ENDIF

!             IF ( lhl > 1 .and. il == lhl ) THEN
!              h1 = cx(mgs,lhl)
!              IF ( cx(mgs,lhl) > 5000. ) THEN
!                write(0,*) 'GS: cx check-pre, chl,q,z time = ',cx(mgs,lhl),qx(mgs,lhl),zx(mgs,lhl),time_real
!              ENDIF
!             ENDIF


         IF ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) THEN ! { .or. qx(mgs,il) <= qxmin(il) 
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3 
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ELSEIF ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3  
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
         ELSE
            IF ( zx(mgs,il) < 0.0 ) THEN !  .and. qx(mgs,il) > 0.05e-3 
               zx(mgs,il) = 0.0
             ENDIF
         ENDIF !}

!             IF ( lhl > 1 .and. il == lhl ) THEN
!              IF ( cx(mgs,lhl) > 5000. ) THEN
!                write(0,*) 'GS: cx check1, chl,ha,q,z time = ',cx(mgs,lhl),h1,qx(mgs,lhl),zx(mgs,lhl),time_real
!              ENDIF
!             ENDIF

         IF (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
        
        IF ( qx(mgs,il) .gt. qxmin(il) ) THEN !{

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        IF ( xv(mgs,il) .lt. xvmn(il) ) THEN
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
!             IF ( lhl > 1 .and. il == lhl ) THEN
!              IF ( cx(mgs,lhl) > 5000. ) THEN
!                write(0,*) 'GS: cx check2, chl,ha,q,z time = ',cx(mgs,lhl),h1,qx(mgs,lhl),zx(mgs,lhl),time_real
!              ENDIF
!             ENDIF
        ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN !{
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
!            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)

!             IF ( lhl > 1 .and. il == lhl ) THEN
!              IF ( cx(mgs,lhl) > 5000. ) THEN
!                write(0,*) 'GS: cx check3, chl,h1,q,z time = ',cx(mgs,lhl),h1,qx(mgs,lhl),zx(mgs,lhl),time_real
!                write(0,*) 'pchli,pchld = ',pchli(mgs),pchld(mgs),dtp*(pchli(mgs)+pchld(mgs))
!                write(0,*) 'an,alpha: = ', an(igs(mgs),jgs,kgs(mgs),lnhl),alpha(mgs,lhl)
!              ENDIF
!             ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
!            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            g1 = (6.0 + alphamax)*(5.0 + alphamax)*(4.0 + alphamax)/ &
     &            ((3.0 + alphamax)*(2.0 + alphamax)*(1.0 + alphamax))
            zx(mgs,il) = Max(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2

!               write(0,*) 'GS: moment problem! il,c,z,q = ',il,cx(mgs,il),zx(mgs,il),qx(mgs,il)
               
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
!               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
               
!               write(0,*) 'GS: moment problem! reset il,c,z,q = ',il,cx(mgs,il),zx(mgs,il),qx(mgs,il)
               
           ELSE
          ! have all valid moments, so find shape parameter
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          IF ( zx(mgs,il) .gt. zxmin .and. qr > qxmin(il) .and. chw > cxmin ) THEN !{

!            rdi = z*(pi/6.*1000.)**2*chw/((rho0(mgs)*qr)**2)
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           DO i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,il))/(Abs(alpha(mgs,il))+1.e-5) .lt. 1. ) EXIT
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO


! check for artificial breakup (graupel/hail larger than allowed max size)
        IF (  xv(mgs,il) .gt. xvmx(il) ) THEN !{
          tmp = cx(mgs,il)


          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          IF ( tmp < cx(mgs,il) ) THEN ! breakup
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF !}

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
 
           IF ( ( lrescalelow(il) .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) THEN !{

            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( lrescalelow(il) .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) .and. &
                     .not. ( il == lr .and. .not. rescale_low_alphar ) ) THEN ! alpha = alphamin, so reset Z to prevent growth in C
             
             wtest = .false.
             IF ( irescalerainopt == 0 ) THEN
               wtest = .false.
             ELSEIF ( irescalerainopt == 1 ) THEN
               wtest = qx(mgs,lc) > qxmin(lc) 
             ELSEIF ( irescalerainopt == 2 ) THEN
               wtest = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             ELSEIF ( irescalerainopt == 3 ) THEN
               wtest = temcg(mgs) > rescale_tempthresh .and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             ENDIF
             
             IF ( il == lr .and. ( wtest .or. .not. rescale_low_alphar ) ) THEN
             ! certain situations where rain number is adjusted instead of Z. Helps avoid rain being 'zapped' by autoconverted 
             ! drops (i.e., favor preserving Z when alpha tries to go negative)
             chw = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2 ! g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1
             cx(mgs,il) = chw
             an(igs(mgs),jy,kgs(mgs),ln(il)) = chw
             ELSE
             ! Usual resetting of reflectivity moment to force consisntency between Q, N, Z, and alpha when alpha = alphamin
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             ENDIF

!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!             z  = z1*(6./(pi*xdn(mgs,il)))**2
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
            ENDIF

           ENDIF !}
          
          
           ENDIF !}
          
           
           ENDIF ! !}
 
          
          
          ENDIF !}

        IF ( lzr > 1 ) THEN
          alpha2d(igs(mgs),kgs(mgs),1) = Max(alphamin, Min(alphamax, alpha(mgs,lr) ))
        ENDIF
        IF ( lzh > 1 ) THEN
          alpha2d(igs(mgs),kgs(mgs),2) = Max(alphamin, Min(alphamax, alpha(mgs,lh) ))
        ENDIF
        IF ( lzhl > 1 ) THEN
          alpha2d(igs(mgs),kgs(mgs),3) = Max(alphamin, Min(alphamax, alpha(mgs,lhl) ))
        ENDIF
        IF ( lzf > 1 ) THEN
          alpha2d(igs(mgs),kgs(mgs),4) = Max(alphamin, Min(alphamax, alpha(mgs,lf) ))
        ENDIF

        IF ( il == lhl .and. lnhlf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lhl) = frac*cx(mgs,lhl)
        ENDIF
        IF ( il == lh .and. lnhf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lh) = frach*cx(mgs,lh)
        ENDIF

!             IF ( lhl > 1 .and. il == lhl ) THEN
!              IF ( cx(mgs,lhl) > 5000. ) THEN
!                write(0,*) 'GS: cx check-post, chl,ha,q,z time = ',cx(mgs,lhl),h1,qx(mgs,lhl),zx(mgs,lhl),time_real
!              ENDIF
!             ENDIF

!      IF ( lf > 0 .and. il == lf .and. kgs(mgs) <= 20 .and. ( cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ) > 200. .or. cx(mgs,lf) > 400. )) THEN
!        write(0,*) 'ix,jy, kz, cf = ',igs(mgs)+ixbeg,jy+jybeg,kgs(mgs), an(igs(mgs),jy,kgs(mgs),ln(lf)),lfsave(mgs,5),lfsave(mgs,6)
!        write(0,*) 'qold,qxold,zold,zxold = ',lfsave(mgs,1),lfsave(mgs,2),lfsave(mgs,3),lfsave(mgs,4)
!        write(0,*) 'cf_new,pcfwi,pcfwd = ',cx(mgs,lf),cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ),pcfwi(mgs) + pcfwd(mgs)
!      
!      ENDIF
        
        ENDDO ! mgs

!         CALL cld_cpu('Z-DELABK')  
        

!         CALL cld_cpu('Z-DELABK')  
        
        
 
           
           ENDIF ! } }

          ENDIF ! }}
          ENDIF ! }

          DO mgs = 1,ngscnt

            IF ( il == lh ) THEN
            IF ( lnhf > 1 ) THEN ! number of graupel from frozen drops
              an(igs(mgs),jy,kgs(mgs),lnhf) = Max( chxf(mgs,lh), 0.0)
            ENDIF
            ENDIF

            IF ( il == lhl ) THEN
            
            IF ( lnhlf > 1 ) THEN ! number of hail from frozen drops
!              an(igs(mgs),jy,kgs(mgs),lnhlf) = Min( cx(mgs,lhl), Max( chxf(mgs,lhl), 0.0) )
              an(igs(mgs),jy,kgs(mgs),lnhlf) = Max( chxf(mgs,lhl), 0.0)
            ENDIF
!              IF ( cx(mgs,il) > 5000. ) THEN
!                write(0,*) 'GS: cx->an, chl, time = ',cx(mgs,il),an(igs(mgs),jy,kgs(mgs),ln(il)) ,pchli(mgs),pchld(mgs),time_real
!              ENDIF
            ENDIF
            an(igs(mgs),jy,kgs(mgs),ln(il)) = Max(cx(mgs,il), 0.0)
          ENDDO
        ENDIF ! }
      ENDDO ! il }

      IF ( lcin > 1 ) THEN
      do mgs = 1,ngscnt
        an(igs(mgs),jy,kgs(mgs),lcin) = Max(0.0, ccin(mgs))
      end do
      ENDIF

      IF ( ipconc .ge. 2 ) THEN
      do mgs = 1,ngscnt
        IF ( lss > 1 ) THEN
          an(igs(mgs),jy,kgs(mgs),lss) = Max(0.0, ssmax(mgs) )
        ENDIF

        IF ( lccn > 1 ) THEN
          an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0, ccnc(mgs) )
        ENDIF
      end do
      ENDIF
      
      ELSEIF ( ipconc .eq. 0 .and. lni .gt. 1 ) THEN
      
          DO mgs = 1,ngscnt
            an(igs(mgs),jy,kgs(mgs),lni) = Max(cx(mgs,li), 0.0)
          ENDDO


      end if

      IF ( ldovol ) THEN

       DO il = li,lhab

        IF ( lvol(il) .ge. 1 ) THEN

          DO mgs = 1,ngscnt

           an(igs(mgs),jy,kgs(mgs),lvol(il)) = Max( 0.0, vx(mgs,il) )
          ENDDO
          
        ENDIF
      
       ENDDO
      
      ENDIF
!
!
!
!
!
      if (ndebug .gt. 0 ) write(0,*) 'gs 12'



! put scx arrays back into an array
      IF ( lscw .gt. 1 .and. ipelec >= 1 ) THEN
      do mgs = 1,ngscnt
!     total space charge ! dont forget to repass array back to AN() as this is used in the driver
         an(igs(mgs),jy,kgs(mgs),lscpi) = cionp(mgs)
         an(igs(mgs),jy,kgs(mgs),lscni) = cionn(mgs)

      if ( largeion ) then

       if ( clionp(mgs) .gt. clionpmx ) then
        an(igs(mgs),jy,kgs(mgs),lscpli) = clionpmx
        an(igs(mgs),jy,kgs(mgs),lscpi)  = an(igs(mgs),jy,kgs(mgs),lscpi) + ( clionp(mgs) - clionpmx )
       else
        an(igs(mgs),jy,kgs(mgs),lscpli) = clionp(mgs)
       endif
       if ( clionn(mgs) .gt. clionnmx ) then
        an(igs(mgs),jy,kgs(mgs),lscnli) = clionnmx
        an(igs(mgs),jy,kgs(mgs),lscni)  = an(igs(mgs),jy,kgs(mgs),lscni) + ( clionn(mgs) - clionnmx )
       else
        an(igs(mgs),jy,kgs(mgs),lscnli) = clionn(mgs)
       endif

      end if

      ENDDO
      DO il = lc,lhab
      do mgs = 1,ngscnt
        tmp = an(igs(mgs),jy,kgs(mgs),lsc(il))
        an(igs(mgs),jy,kgs(mgs),lsc(il)) = scx(mgs,il)
!        IF ( .not. (scx(mgs,il) > -1.e-6 .and. scx(mgs,il) < 1.e-6 ) ) THEN
!        IF ( Abs (scx(mgs,il) ) > 1000.e-9 .or. Abs(scx(mgs,il)) > 1.e-9 ) THEN ! DEBUGTED
        IF ( Abs (scx(mgs,il) ) > 1000.e-9  ) THEN
          write(0,*) 'Problem2 with scx il = ',il,scx(mgs,il),tmp,mgs,tmp-scx(mgs,il),(tmp-scx(mgs,il))*dtpinv
          write(0,*) 'temper = ',temcg(mgs)
          write(0,*) 'qx,cx = ',qx(mgs,il),cx(mgs,il)
        write(0,*)  psccwi(mgs), psccwd(mgs)
        write(0,*)  psccii(mgs), psccid(mgs)
        write(0,*)  pscrwi(mgs), pscrwd(mgs)
        write(0,*)  pscswi(mgs), pscswd(mgs)
        write(0,*)  pschwi(mgs), pschwd(mgs)
        write(0,*)  pschli(mgs), pschld(mgs)
        write(0,*)  pscpii(mgs) ,pscpid(mgs)
        write(0,*)  pscnii(mgs) ,pscnid(mgs)
        write(0,*)
        write(0,*)  psccwmi(mgs), psccwmd(mgs)
        write(0,*)  psccimi(mgs), psccimd(mgs)
        write(0,*)  pscrwmi(mgs), pscrwmd(mgs)
        write(0,*)  pscswmi(mgs), pscswmd(mgs)
        write(0,*)  pschwmi(mgs), pschwmd(mgs)
        write(0,*)  pschlmi(mgs), pschlmd(mgs)

        write(0,*) 'ci: tot, parts:', dtp*(psccii(mgs)+psccid(mgs)+psccimi(mgs)+psccimd(mgs)), &
     &             dtp*psccii(mgs),dtp*psccid(mgs),dtp*psccimi(mgs),dtp*psccimd(mgs)

        write(0,*) 'rain: tot, parts:', dtp*(pscrwi(mgs)+pscrwd(mgs)+pscrwmi(mgs)+pscrwmd(mgs)), &
     &             dtp*pscrwi(mgs), dtp*pscrwd(mgs), dtp*pscrwmi(mgs), dtp*pscrwmd(mgs)
        
        write(0,*) +fsccw(mgs)*qracw(mgs) ,fsccw(mgs)*qrcnw(mgs) 
        write(0,*)  fschl(mgs)*qhlmlr(mgs), fschw(mgs)*qhmlr(mgs), fscsw(mgs)*qsmlr(mgs)
        write(0,*)  fschl2(mgs)*qhlshr(mgs), fschw2(mgs)*qhshr(mgs) , fscsw2(mgs)*qsshr(mgs)
        write(0,*) fschw2(mgs), qhshr(mgs)
        write(0,*) 'qh,sch = ',an(igs(mgs),jgs,kgs(mgs),lh),an(igs(mgs),jgs,kgs(mgs),lsch)
        write(0,*) 'scx = ',scx(mgs,lh)

        write(0,*) 'graupel chg: new, old = ',scx(mgs,lh),scx(mgs,lh) -   &
    &   dtp*(pschwi(mgs)+pschwd(mgs)+pschwmi(mgs)+pschwmd(mgs))
        write(0,*) 'graupel rates: ',pschwi(mgs),pschwd(mgs),pschwmi(mgs),pschwmd(mgs)
        write(0,*) 'pschwmi parts'
        write(0,*)  +il5(mgs)*(ifiacrg*(fscci(mgs)*qracif(mgs)+fscrw(mgs)*qiacrf(mgs))   &
     &  +fscrw(mgs)*ifrzg*qrfrz(mgs))   
        write(0,*)   +fscrw(mgs)*qhacr(mgs)+fsccw(mgs)*qhacw(mgs)   
        write(0,*)   +fscsw(mgs)*qhacs(mgs)+fscci(mgs)*qhaci(mgs)   
        write(0,*) 'pschwmd parts'
        write(0,*)   - fschw(mgs)*qhlcnh(mgs),fschw(mgs),qhlcnh(mgs)  
!     &  + fschw(mgs)*qhshr(mgs)    &
        write(0,*)   +(1-il5(mgs))*fschw(mgs)*qhmlr(mgs)*mixedphasefac   
        write(0,*)   -fschw(mgs)*qhmul1(mgs), fschw(mgs), qhmul1(mgs)


        write(0,*) 'pschlmi parts'
        write(0,*)  +fscrw(mgs)*qhlacr(mgs),fsccw(mgs)*qhlacw(mgs)   
        write(0,*)  +fscsw(mgs)*qhlacs(mgs),fscci(mgs)*qhlaci(mgs)   
        write(0,*)  + fschw(mgs)*qhlcnh(mgs) ,fschw(mgs),qhlcnh(mgs) 

     write(0,*) 'qs,qsold: ',qx(mgs,ls), qx(mgs,ls) - dtp*(pqswi(mgs)+pqswd(mgs))
     write(0,*) 'cs,csold: ',cx(mgs,ls), cx(mgs,ls) - dtp*(pcswi(mgs)+pcswd(mgs))

     
     write(0,*)  'pscswmd(mgs)',fscsw(mgs),fscsw2(mgs)
     write(0,*) - fscsw(mgs)*qhacs(mgs), - fscsw(mgs)*qhlacs(mgs)
     write(0,*)  +(1-il5(mgs))*fscsw(mgs)*qsmlr(mgs)*mixedphasefac , fscrw(mgs)*qsshr(mgs)
     write(0,*)  -fscsw(mgs)*qsmul(mgs)
     write(0,*) qhacs(mgs), qhlacs(mgs),qsmlr(mgs),qsshr(mgs) ,qsmul(mgs), mixedphasefac
     write(0,*)
     write(0,*)  fschw(mgs),fscsw(mgs),fschl(mgs)
     write(0,*)  fsccw(mgs),fscci(mgs),fscrw(mgs)

     write(0,*)  'fscsw2(mgs)*qsshr(mgs), ', fscsw2(mgs)*qsshr(mgs)

      tmp =  qx(mgs,ls) + dtp*(pqswi(mgs) + pqswd(mgs) - qsshr(mgs)) ! should this add back the shed rain mass? probably!
      chgtmp = 0.0
      IF ( tmp > qxmin(ls) .and. qsshr(mgs) < 0.0 ) THEN
        chgtmp = scx(mgs,ls) + dtp*(pscswmi(mgs) + pscswmd(mgs))
       ! fscsw2(mgs) = chgtmp/tmp
       ! pscswmd(mgs) =  pscswmd(mgs) + fscsw2(mgs)*qsshr(mgs)
      ENDIF
      write(0,*) 'tmp,chgtmp, fscsw2, qsshr,fscsw2(mgs)*qsshr(mgs) ',tmp,chgtmp,fscsw2(mgs), qsshr(mgs), fscsw2(mgs)*qsshr(mgs)


        STOP
        ENDIF
      end do
      ENDDO

      ENDIF

!      write(0,*) 'SCW',MAXVAL(scx(:,lc)),MAXLOC(scx(:,lc)),MINVAL(scx(:,lc)),MINLOC(scx(:,lc))



      if (ndebug .gt. 0 ) write(0,*) 'gs 13'

 9998 continue

      if ( kz .gt. nz-1 .and. ix .ge. itile) then
        if ( ix .ge. itile ) then
         go to 1200 ! exit gather scatter
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz
      end if

      if ( ix .ge. itile ) then
        nxmpb = 1
        nzmpb = kz+1
      else
       nxmpb = ix+1
      end if

 1000 continue
 1200 continue
!
!  end of gather scatter (for this jy slice)
!
!

      return
      end subroutine nssl_2mom_gs
!
!--------------------------------------------------------------------------
!

      real function  galpha(a_in) 
      implicit none
      real :: a_in
        galpha = ((4. + a_in)*(5. + a_in)*(6. + a_in))/((1. + a_in)*(2. + a_in)*(3. + a_in))
      end function galpha
!
!--------------------------------------------------------------------------
!

      real function dgalpha(a_in) 
      real :: a_in
        dgalpha = (876. + 1260.*a_in + 621.*a_in**2 + 126.*a_in**3 + 9.*a_in**4)/            &
     &  (36. + 132.*a_in + 193.*a_in**2 + 144.*a_in**3 + 58.*a_in**4 + 12.*a_in**5 + a_in**6)
      end function dgalpha
!
!--------------------------------------------------------------------------
!
!  Calculate reflectivity change when only number changes
!  Differential version can have large error when crate is big and time step is big
      real function zraten(dtpinv,dtp,g1x,rho0,xdn,qx,cx,crate)
      implicit none
      real, intent(in) :: dtpinv,dtp,g1x,rho0,qx,cx,crate,xdn
      real, parameter :: pi = 3.141592653589793

        IF ( cx > 1.e-8 ) THEN
          zraten = (6./pi)**2*dtpinv*g1x*(rho0*qx/xdn)**2 &
                     *  crate /((cx + dtp*crate)*cx)
        ELSE
          zraten = 0.0
        ENDIF
      end function zraten
!
!--------------------------------------------------------------------------
!
!  Calculate reflectivity change when only mass changes
      real function zrateq(dtpinv,dtp,g1x,rho0,xdn,qx,cx,qrate)
      implicit none
      real, intent(in) :: dtpinv,dtp,g1x,rho0,qx,cx,qrate,xdn
      real :: tmp1,tmp2
      real, parameter :: pi = 3.141592653589793

        IF ( cx > 1.e-8 ) THEN
          tmp1 = qx**2
          tmp2 = (qx+dtp*qrate)**2
          zrateq = (6./pi)**2*dtpinv*g1x*(rho0/xdn)**2*(tmp2 - tmp1)/cx
        ELSE
          zrateq = 0.0
        ENDIF

      end function zrateq
!
!--------------------------------------------------------------------------
!
!  Calculate reflectivity change when both mass and  number change
      real function zrateqn(dtpinv,dtp,g1x,rho0,xdn,qx,cx,crate,qrate)
      implicit none
      real, intent(in) :: dtpinv,dtp,g1x,rho0,qx,cx,crate,xdn,qrate
      integer :: ioldnew = 1
      real :: tmp1,tmp2
      real, parameter :: pi = 3.141592653589793

        IF ( cx > 1.e-8 ) THEN
          IF ( ioldnew == 1 .and. cx + dtp*crate > 1.e-8 .and. qx+dtp*qrate > 0. ) THEN
          ! final-initial
          tmp1 = qx**2/cx
          tmp2 = (qx+dtp*qrate)**2/(cx + dtp*crate)
          zrateqn = (6./pi)**2*dtpinv*g1x*(rho0/xdn)**2*(tmp2 - tmp1)
          ELSE
          ! differential form
           tmp1 = qx/cx
           zrateqn = (6./pi)**2*g1x*(rho0/xdn)**2*( 2.*tmp1*qrate - tmp1**2 * crate )
          ENDIF
        ELSE
          zrateqn = 0.0
        ENDIF

      end function zrateqn
!
!--------------------------------------------------------------------------
!
