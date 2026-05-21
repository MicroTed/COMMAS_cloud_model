#if defined(COMMAS)
#define CHGELEC
#define Z3MOM
#define USEFROZENDROPS 1
#define USEICESPHERES 1
#define USERAINTYPES 1
#define NUWRFMODS 1
#endif
! #####################################################################
! #####################################################################
!
! Subroutine for explicit cloud condensation and droplet nucleation
!
! 11/30/2022: Fixed droplet evaporation heating term for CM1 eqtset=2 (was only doing eqtset=1)
!
   SUBROUTINE NUCOND    &
     &  (nx,ny,nz,na,jyslab & 
     &  ,nor,norz,dtp,nxi & 
#ifdef COMMAS
     &  ,dz1d & 
#else
     &  ,dz3d & 
#endif
     &  ,t0,t9 & 
     &  ,an,dn,p2 & 
     &  ,pn,w & 
#ifdef COMMAS
     &  ,ventr,ventrn,ventc,qxmin, xdn0, xvmn, xvmx, cno  &
     &  ,cckm,ccne,ccnefac,cnexp,ccne0  &
     &  ,ido,lsc,ln,ipc,lvol,lz,lliq,lrain &
     &  ,xdnmn,xdnmx  &
     &  ,pb, pinit, thproc,numproc, dxx,dyy,dzz    &
#else
     &  ,ngs   &
#endif
#ifdef WRFCODE
     &  ,thproc,numproc, dx1,dy1,gz    &
#endif
     &  ,axtra,io_flag &
     &  ,ssfilt,t00,t77,flag_qndrop  &
     & )

#ifdef COMMAS
       USE GRID_MODULE
       USE INDEX_MODULE !, cwmasn_index => cwmasn
       USE MICRO_MODULE
       USE COMMASMPI_MODULE, only : nzend_commas => nzend, kzbeg_commas => kzbeg,ixbeg,jybeg,my_rank
       USE TRAJ_MODULE
!       USE PARAM_MODULE, only : gr => g
       use param_module, only : gr => g, rd, cp, cap => rcp, poo => p00, rdorv => epsilon
#endif

   implicit none

!      real :: cwmasn = 1000.*0.523599*(2.*2.e-6)**3 
      integer :: nx,ny,nz,na,nxi
      integer :: nor,norz, jyslab ! ,nht,ngt,igsr
      real    :: dtp  ! time step
      logical :: flag_qndrop

      integer, parameter :: ng1 = 1


!
! external temporary arrays
!
      real t00(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t77(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t3(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t4(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t5(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t6(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t8(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t9(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      

      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)  ! perturbation Pi
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real qv(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real ssfilt(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      

      real pb(-norz+ng1:nz+norz)
      real pinit(-norz+ng1:nz+norz)

#ifdef COMMAS
      real dz1d(-norz+1:nz+norz)
#else
      real dz3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
#endif

      
    ! local

#ifdef COMMAS
      ! passed in
      real ventr,ventrn,ventc
      real qxmin(lc:lqmx)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real cno(lc:lhab)
      real cckm,ccne,ccnefac,cnexp,ccne0
      integer ido(lc:lqmx)
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      integer lvol(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
#if USERAINTYPES > 0
      integer lrain(nraintypes)
#endif
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      
      ! local
      
      real tfr,tfrh
      parameter ( tfr = 273.15, tfrh = 233.15)
      
!       real cp, 
!       parameter ( cp = 1004.0, rd = 287.04 )
      
      real cpi
      parameter ( cpi = 1./cp )
      
!       real poo,cap
!       parameter ( cap = rd/cp, poo = 1.0e+05 )

      real, parameter :: ec = 1.602e-19 ! fundamental unit of charge
      real, parameter :: eci = 1.0/ec

      real, parameter :: rw = 461.5              ! gas const. for water vapor
      real, parameter :: advisc0 = 1.832e-05     ! reference dynamic viscosity (SMT; see Beard & Pruppacher 71)
      real, parameter :: advisc1 = 1.718e-05     ! dynamic viscosity constant used in thermal conductivity calc
      real, parameter :: tka0 = 2.43e-02         ! reference thermal conductivity
      
      integer, parameter :: eqtset = 0

      real, parameter ::      cv = 717.0             ! specific heat at constant volume - air
      REAL, parameter ::      cvv = 1408.5
      REAL, parameter ::      cpl = 4190.0
      REAL, parameter ::      cpigb = 2106.0

      REAL, parameter :: rho00 = 1.225          ! reference/MSL air density

      real, parameter :: pi = 3.141592653589793
      real, parameter :: piinv = 1./pi
      real, parameter :: c1f3 = 1.0/3.0
      real, parameter :: cwc1 = 6.0/(pi*1000.)

      real, parameter :: ar = 841.99666         ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: br = 0.8               ! rain terminal velocity power law coefficient (LFO)

      integer numproc
      real thproc(nzend_commas,numproc)

      real dzz(nz),dxx(nx),dyy(ny)          ! dz(k),dx(i),dy(j)
      
      logical, parameter :: is_aerosol_aware = .false.

      real, parameter :: rovcp = rd/cp
      
      real, external :: beta, POLYSVP1

#endif

#ifdef WRFCODE
      integer, intent(in) :: numproc
      real, intent(inout) :: thproc(nz,numproc)
      real, intent(in) :: dx1,dy1, gz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
#endif
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)
      logical :: io_flag
      
      real :: dv
      real :: ccnefactwo, sstmp, cn1, cnuctmp

! 
!  declarations microphysics and for gather/scatter
!
      real, parameter :: cwmas30 = 1000.*0.523599*(2.*30.e-6)**3 ! mass of 30-micron radius droplet, for sat. adj.
      real, parameter :: cwmas20 = 1000.*0.523599*(2.*20.e-6)**3 ! mass of 20-micron radius droplet, for sat. adj.
      integer nxmpb,nzmpb,nxz
      integer mgs,ngs,numgs,inumgs
#if defined( COMMAS )
      parameter (ngs=64)
#endif
      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs)
      integer nsvcnt
      
      integer ix,kz,i,n, kp1, km1
      integer :: jy, jgs
      integer ixb,ixe,jyb,jye,kzb,kze
    
      integer itile,jtile,ktile
      integer ixend,jyend,kzend,kzbeg
      integer nxend,nyend,nzend,nzbeg

!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs), ccna(ngs), cnuc(ngs), cwnccn(ngs), ccnaco(ngs), ccnanu(ngs)
      real :: ccnc_nu(ngs), ccnc_ac(ngs), ccnc_co(ngs)
      real ccncuf(ngs)
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
      real :: ssmax(ngs)      ! maximum SS experienced by a parcel
      real ssmx
      real dnnet,dqnet
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real ventrx(ngs)
      real ventrxn(ngs)
      real volb, t2s
      real, parameter :: aa1 = 9.44e15, aa2 = 5.78e3  ! a1 in Ziegler

      real rhoinv(ngs)
      
      real chw, g1, rd1

      real ac1,bc, taus, c1,d1,e1,f1,p380,tmp,tmp2 ! , sstdy, super
      real tmpmx, fw, qctmp
      real x,y,del,r,alpr
      double precision :: vent1,vent2
      real g1palp
      real bs
      real v1, v2
      real d1r, d1i, d1s, e1i
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
      real dcrit
      real cn(ngs), cnuf(ngs)
      real :: ccwmax
      

      integer ltemq
      
      integer il

      real  es(ngs) ! ss(ngs),
!      real  eis(ngs)
      real ssf(ngs),ssfkp1(ngs),ssfkm1(ngs),ssat0(ngs)
      real, parameter :: ssfcut = 4.0
      real ssfjp1(ngs),ssfjm1(ngs)
      real ssfip1(ngs),ssfim1(ngs)

      real supcb, supmx
      parameter (supcb=0.5,supmx=238.0)
      real r2dxm, r2dym, r2dzm
      real dssdz, dssdy, dssdx
!      real tqvcon
      real epsi,d
      parameter (epsi = 0.622, d = 0.266)
      real r1,qevap ! ,slv
      
      real vr,nrx,qr,z1,z2,rdi,alp,xnutmp,xnuc
      real ctmp, ccwtmp
      real f5, qvs0  ! Kessler condensation factor
      real    :: t0p1, t0p3
      real qvex
      
!      real, dimension(ngs) :: temp, tempc, elv, elf, els, pqs, theta, temg, temcg
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
      real temp(ngs),tempc(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) ! ,tembzg(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real felv(ngs),felf(ngs),fels(ngs)
      real felvcp(ngs),felvpi(ngs)
      real gamw(ngs),gams(ngs)   !   qciavl(ngs),
      real tsqr(ngs),ssi(ngs),ssw(ngs)
      real cc3(ngs),cqv1(ngs),cqv2(ngs)
      real qcwtmp(ngs),qtmp

      real fvent(ngs) !,fraci(ngs),fracl(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs)
      real fschm(ngs),fpndl(ngs)

      real pres(ngs),pipert(ngs)
      real pk(ngs)
      real rho0(ngs),pi0(ngs)
      real rhovt(ngs)
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
      real qss0(ngs)
      real fcqv1(ngs)
      real wvel(ngs),wvelkm1(ngs)

      real wvdf(ngs),tka(ngs)
      real advisc(ngs)

      real rwvent(ngs)
      

      real :: qx(ngs,lv:lhab)
      real :: cx(ngs,lc:lhab)
      real :: xv(ngs,lc:lhab)
      real :: xmas(ngs,lc:lhab)
      real :: xdn(ngs,lc:lhab)
      real :: xdia(ngs,lc:lhab,3)
      real :: alpha(ngs,lc:lhab)
      real :: zx(ngs,lr:lhab)


      logical zerocx(lc:lqmx)
      
      logical :: lprint

      integer, parameter :: iunit = 0
      
      real :: frac, hwdn, tmpg, xdia1, xdia3, cwch,xvol
      
      real :: cvm,cpm,rmm

      real, parameter ::      cpv = 1885.0       ! specific heat of water vapor at constant pressure
      real, parameter ::      Mair = 0.0284        ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)


      integer :: kstag
      
      integer :: count
      
!     Addtion T.Iguchi Y2021 Update
      real, parameter :: mwwater = 0.01801528  ! Molecular weight of water (kg/mol)
      real, parameter :: rhowater = 997.0  ! Density of liquid water (kg/m3)
      real, parameter :: gasconst = 8.3144598   ! Gas constant (m2 kg s-2 K-1 mol-1)
      real :: sswater  ! unit change supersaturation from percentage to n/a
      real :: sigvl, aact

      real :: alpha_ar, gamma_ar, G_ar, evs, zeta, smax
      real :: f_ac, g_ac, eta_ac
      real :: f_nu, g_nu, eta_nu
      real :: f_co, g_co, eta_co

      real :: sm_nu, sm_ac, sm_co, ss_ac, ss_nu, ss_co
      real :: uu_nu, uu_ac, uu_co
      
      real :: cn_ac, cn_co, cn_nu
#ifdef COMMAS
      real, external :: Derf
#endif

! -------------------------------------------------------------------------------
      itile = nxi
      jtile = ny
      ktile = nz
      ixend = nxi
      jyend = ny
      kzend = nz
      nxend = nxi + 1
      nyend = ny + 1
      nzend = nz
#ifdef COMMAS
      kzbeg = kzbeg_commas
#else
      kzbeg = 1
#endif
      nzbeg = 1

      IF ( ac_opt > 0 )  ccnefactwo =  (1.63e-3/(cck * beta(3./2., cck/2.)))**(1.0/(cck + 2.0))
      f5 = 237.3 * 17.27 * 2.5e6 / cp ! combined constants for rain condensation (Soong and Ogura 73)

#ifdef COMMAS
       jy = jyslab
       kstag = 1
#else
       jy = 1
       kstag = 0
       pb(:) = 0.0
       pinit(:) = 0.0
#endif
      
      IF ( ipconc <= 1 .or. isedonly == 2 ) GOTO 2200

!
!  Ziegler nucleation 
!

!      ssfilt(:,:,:) = 0.0
      ssmx = 0
      count = 0

      do kz = 1,nz-kstag
        do ix = 1,nxi

         temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz)
          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

!          c1 = t00(ix,jy,kz)*tabqvs(ltemq)
          IF ( iqvsopt == 0 ) THEN
            c1 = t00(ix,jy,kz)*tabqvs(ltemq)
          ELSEIF ( iqvsopt == 1 ) THEN
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pn(ix,jy,kz) + pb(kz) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
          ELSE
            tmp = min(0.99*(pn(ix,jy,kz)+pb(kz)),POLYSVP1(temp1,0))
            c1 = rdorv*tmp/(pn(ix,jy,kz) + pb(kz) - tmp)
#endif
          ENDIF

          IF ( c1 > 0. ) THEN
            ssfilt(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
          ELSE
            ssfilt(ix,jy,kz) = -100.
          ENDIF

        ENDDO
      ENDDO


!
!     jy = 1 ! working on a 2d slab
!!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

!
!..Gather microphysics
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Gather stage'

      nxmpb = 1
      nzmpb = 1
      nxz = nxi*nz
      numgs = nxz/ngs + 1


      do 2000 inumgs = 1,numgs

      ngscnt = 0


      kzb = nzmpb
      kze = nz-kstag
 !     if (kzbeg .le. nzmpb .and. kzend .gt. nzmpb) kzb = nzmpb

      ixb = nxmpb
      ixe = itile

      do kz = kzb,kze
      do ix = nxmpb,nxi

      pres(1) = pn(ix,jy,kz) + pb(kz)
      pqs(1) = 380.0/(pn(ix,jy,kz) + pb(kz))
      theta(1) = an(ix,jy,kz,lt)
      temg(1) = t0(ix,jy,kz)

      temcg(1) = temg(1) - tfr
      ltemq = (temg(1)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
     ! qvs(1) = pqs(1)*tabqvs(ltemq)
      IF ( iqvsopt == 0 ) THEN
        qvs(1) = pqs(1)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(1) = rdorv*esbolton*tabqvs(ltemq)/(pres(1) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
      ELSE
        tmp = min(0.99*pres(1),POLYSVP1(temg(1),0))
        qvs(1) = rdorv*tmp/(pres(1) - tmp)
#endif
      ENDIF
      qis(1) = pqs(1)*tabqis(ltemq)

      qss(1) = qvs(1)


      if ( temg(1) .lt. tfr ) then
      end if
!
      if ( (temg(1) .gt. tfrh .or. an(ix,jy,kz,lv)/qvs(1) > maxlowtempss ) .and.  &
     &   ( an(ix,jy,kz,lv)  .gt. qss(1) .or. &
     &     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.  &
     &     ( an(ix,jy,kz,lr)  .gt. qxmin(lr) .and. rcond == 2 )  &
     &     )) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 2100
      end if

      end do  !ix

      nxmpb = 1
      end do  !kz
!      if ( jy .eq. (ny-jstag) ) iend = 1
 2100 continue

      if ( ngscnt .eq. 0 ) go to 29998

      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: dbg = 8'
      
!      write(0,*) 'NUCOND: dbg = 8, ngscnt,ssmx = ',ngscnt,ssmx

      
      qx(:,:) = 0.0
      cx(:,:) = 0.0
#ifdef Z3MOM
      zx(:,:) = 0.0
#endif

      xv(:,:) = 0.0
      xmas(:,:) = 0.0

      IF ( imurain == 1 ) THEN
        alpha(:,lr) = alphar
      ELSEIF ( imurain == 3 ) THEN
        alpha(:,lr) = xnu(lr)
      ENDIF

!
!  define temporaries for state variables to be used in calculations
!
      DO mgs = 1,ngscnt
      qx(mgs,lv) = an(igs(mgs),jy,kgs(mgs),lv)
       DO il = lc,lhab
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0)
       ENDDO

       qcwtmp(mgs) = qx(mgs,lc)


      theta0(mgs) = an(igs(mgs),jy,kgs(mgs),lt) !
      thetap(mgs) = 0.0
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      qv0(mgs) =  qx(mgs,lv)
      qwvp(mgs) = qx(mgs,lv) - qv0(mgs)

       pres(mgs) = pn(igs(mgs),jy,kgs(mgs)) + pb(kgs(mgs))
       pipert(mgs) = p2(igs(mgs),jy,kgs(mgs))
       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhoinv(mgs) = 1.0/rho0(mgs)
       rhovt(mgs) = Sqrt(rho00/rho0(mgs))
       pi0(mgs) = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
!       pk(mgs) = t77(igs(mgs),jy,kgs(mgs)) ! ( pres(mgs) / poo ) ** cap
       pk(mgs)   = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs)) ! t77(igs(mgs),jy,kgs(mgs))
       temcg(mgs) = temg(mgs) - tfr
       qss0(mgs) = (380.0)/(pres(mgs))
       pqs(mgs) = (380.0)/(pres(mgs))
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
!       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       IF ( iqvsopt == 0 ) THEN
         qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       ELSEIF ( iqvsopt == 1 ) THEN
         qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
       ELSE
         tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
         qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
#endif
       ENDIF
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
!
        qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
        IF ( iqvsopt == 0 ) THEN
          es(mgs) = 6.1078e2*tabqvs(ltemq)
        ELSEIF ( iqvsopt == 1 ) THEN
          es(mgs) = esbolton*tabqvs(ltemq)
#if defined(COMMAS) || defined(COMMASTMP)
        ELSE
          es(mgs) = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
          ! qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
#endif
        ENDIF
!        es(mgs) = 6.1078e2*tabqvs(ltemq)
        qss(mgs) = qvs(mgs)


        temgx(mgs) = min(temg(mgs),313.15)
        temgx(mgs) = max(temgx(mgs),233.15)
        felv(mgs) = 2500837.367 * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
!
        IF ( eqtset <= 1 ) THEN
          felvcp(mgs) = felv(mgs)*cpi
        ELSE ! equation set 2 in cm1
          tmp = qx(mgs,li)+qx(mgs,ls)+qx(mgs,lh)
          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          IF ( lf > 1 ) tmp = tmp + qx(mgs,lf)
          cvm = cv+cvv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          cpm = cp+cpv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          rmm=rd+rw*qx(mgs,lv)
          
          IF ( eqtset == 2 ) THEN

           felvcp(mgs) = (felv(mgs)-rw*temg(mgs))/cvm

          ELSE
            felvcp(mgs) = (felv(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
            felvpi(mgs) = pi0(mgs)*rovcp*(felv(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm
          ENDIF

        ENDIF

        temcgx(mgs) = min(temg(mgs),273.15)
        temcgx(mgs) = max(temcgx(mgs),223.15)
        temcgx(mgs) = temcgx(mgs)-273.15
        felf(mgs) = 333690.6098 + (2030.61425)*temcgx(mgs) - (10.46708312)*temcgx(mgs)**2
!
        fels(mgs) = felv(mgs) + felf(mgs)
        fcqv1(mgs) = 4098.0258*felv(mgs)*cpi

      wvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)* &
     &  (101325.0/(pb(kgs(mgs)) + pn(igs(mgs),jgs,kgs(mgs))))                            ! diffusivity of water vapor, Hall and Pruppacher (76)
      advisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &  (temg(mgs)/296.0)**(1.5)                         ! dynamic viscosity (SMT; see Beard & Pruppacher 71)
      tka(mgs) = tka0*advisc(mgs)/advisc1                 ! thermal conductivity


      ENDDO



!
! load concentrations
!
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cwnccn(mgs) = cwccn*rho0(mgs)/rho00 ! background ccn count
        cn(mgs) = 0.0
        IF ( lss > 1 ) THEN 
          ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        ELSE
          ssmax(mgs) = 0.0
        ENDIF
        IF ( lccn .gt. 1 .and. ac_opt == 0 ) THEN
          IF ( lccnuf .gt. 1 .and. i_uf_or_ccn > 0 ) THEN
             ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn) + an(igs(mgs),jy,kgs(mgs),lccnuf)
          ELSE
             ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
             IF ( lccna > 1 ) THEN
               cnuc(mgs) = ccnc(mgs)
             ENDIF
          ENDIF
#ifdef NUWRFMODS
        ELSEIF ( lcn_ac > 1 .and. ( ac_opt == 1 .or. ac_opt == 11 ) ) THEN
             ccnc_ac(mgs) = an(igs(mgs),jy,kgs(mgs),lcn_ac)
             cnuc(mgs) = ccnc_ac(mgs)
           !  write(0,*) 'ccnc_ac,mgs = ', ccnc_ac(mgs),mgs,igs(mgs),jy,kgs(mgs)
        ELSEIF ( lcn_ac > 1 .and. ( ac_opt == 2 .or. ac_opt == 22 ) ) THEN
          ccnc_nu(mgs) = an(igs(mgs),jy,kgs(mgs),lcn_nu)
          ccnc_ac(mgs) = an(igs(mgs),jy,kgs(mgs),lcn_ac)
          ccnc_co(mgs) = an(igs(mgs),jy,kgs(mgs),lcn_co)
#endif
        ELSE
          ccnc(mgs) = cwnccn(mgs)
        ENDIF
        IF ( lccnuf .gt. 1 .and. i_uf_or_ccn == 0 ) THEN
          ccncuf(mgs) = an(igs(mgs),jy,kgs(mgs),lccnuf)
        ELSE
          ccncuf(mgs) = 0.0
        ENDIF
        cnuf(mgs) = 0.0
        IF ( lccna > 1 ) THEN
          ccna(mgs) = an(igs(mgs),jy,kgs(mgs),lccna) ! predicted count of activated ccn
          IF ( ac_opt == 22 ) THEN
            IF ( lccnaco > 1 ) THEN 
              ccnaco(mgs) = an(igs(mgs),jy,kgs(mgs),lccnaco)
            ELSE
              ccnaco(mgs) = 0.0
            ENDIF
            IF ( lccnanu > 1 ) THEN
              ccnanu(mgs) = an(igs(mgs),jy,kgs(mgs),lccnanu)
            ELSE
              ccnanu(mgs) = 0.0
            ENDIF
          ENDIF
        ELSE
          IF ( lccn > 1 ) THEN
#ifdef COMMAS
            ccna(mgs) = cwnccn(mgs) - ccnc(mgs) ! diagnose activated ccn as background value - remaining unactivated ccn
#else
            ccna(mgs) = 0.0 ! WRF driver interface already has ccw subtracted from ccnc
#endif
          ELSE
            ccna(mgs) = cx(mgs,lc) ! approximation of number of activated ccn
          ENDIF
        ENDIF
       end do
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
       end do
      end if

!        cnuc(1:ngscnt) = cwccn*rho0(mgs)/rho00*(1. - renucfrac) + ccnc(1:ngscnt)*renucfrac
       DO mgs = 1,ngscnt
       ! default value of renucfrac is 0.0
        IF ( irenuc /= 6 ) THEN
          IF ( irenuc == 2 ) THEN
            cnuc(mgs) = Max(ccnc(mgs),cwnccn(mgs))*(1. - renucfrac) + ccnc(mgs)*renucfrac
          ELSE
            cnuc(mgs) = ccnc(mgs)*(1. - renucfrac) + ccnc(mgs)*renucfrac
          ENDIF
        ELSE
        cnuc(mgs) = ccnc(mgs)*(1. - renucfrac) + Max(0.0,ccnc(mgs) - ccna(mgs))*renucfrac
        ENDIF
        IF ( renucfrac >= 0.999 ) THEN
          IF ( temg(mgs) < 265. ) THEN
            IF ( qx(mgs,lc) > 10.*qxmin(lc) .and. w(igs(mgs),jgs,kgs(mgs)) > 2.0 ) THEN
             cnuc(mgs) = 0.0 !  Min(cnuc(mgs), 0.5*cx(mgs,lc) ) ! Hack to reduce nucleation at low temp in updraft when ccn are not predicted
            ELSE
             cnuc(mgs) = 0.1*cnuc(mgs)
            ENDIF
          ENDIF
        ENDIF
       ENDDO

!  Set density
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Set density'

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
      end do

      ventrx(:) = ventr
      ventrxn(:) = ventrn
      
#ifdef Z3MOM

!  Find shape parameter rain

      IF ( lzr > 1 .and. rcond == 2 ) THEN ! { RAIN SHAPE PARAM
      DO mgs = 1,ngscnt
         zx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lzr), 0.0)
      ENDDO

!      CALL cld_cpu('Z-MOMENT-1r2')
          il = lr
          DO mgs = 1,ngscnt

         IF ( iresetmoments == 1 .or. iresetmoments == il  .or. iresetmoments == -1 ) THEN
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3 ) THEN
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
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3  ) THEN
!!            write(91,*) 'cx=0; qx,zx = ',1000.*qx(mgs,il),1.e18*zx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
         ENDIF

         IF (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
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
            xv(mgs,lr) = xvmx(lr)
            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
          ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
            xv(mgs,lr) = xvmn(lr)
            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            IF ( imurain == 3 ) THEN
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z1   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
            ELSE
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
            z1   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
            
            ENDIF
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            IF ( imurain == 3 ) THEN
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(chw*1000.*1000)
            ELSE
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(chw*1000.*1000)
            
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
              IF ( imurain == 3 ) THEN
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z1   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              ELSEIF ( imurain == 1 ) THEN
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z1   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z1*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              
              ENDIF
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

          IF ( imurain == 3 ) THEN
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           write(0,*) 'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z1,xv
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

         ELSE ! imurain == 1
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)

            rd1 = z1*(pi/6.*xdn(mgs,il))**2*nrx/(rho0(mgs)*qr)**2

           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd1) - 1.0

           DO i = 1,10
            IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )

             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd1) - 1.0

             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

         
         ENDIF
!         ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
          IF ( imurain == 3 ) THEN
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
           
          ELSEIF ( imurain == 1 ) THEN
          
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

           IF ( (rescale_low_alpha .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) THEN



            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*rho0(mgs)**2*(qr)*qr/zx(mgs,lr)*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alpha .and. alp <= alphamin ) THEN ! alpha = alphamin, so reset Z to prevent growth in C
             z1 = g1*rho0(mgs)**2*(qr)*qr/nrx
             z2  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z2
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z2
            ENDIF
          ENDIF ! imurain

          ENDIF ! z > 0

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

           IF ( imurain == 3 .and. izwisventr == 2 ) THEN

           tmp = alpha(mgs,lr) + 1.5 + br/6.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 1.5 + br/6.)/Gamma(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/(y*(alpha(mgs,lr) + 1.)**((1.+br)/6. + 1./3.))
           
           ELSEIF ( imurain == 1 .and.  iferwisventr == 2 ) THEN

           tmp = alpha(mgs,lr) + 2.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 1.5 + br/6.)/Gamma(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/y
           
           
           ENDIF

           
           ENDIF
          ENDIF
          
          ENDIF
          
          ENDDO
!        CALL cld_cpu('Z-MOMENT-1r2')  
        ENDIF ! }
#endif


!       write(0,*) 'NUCOND: Set ssf variables, ssmxinit =',ssmxinit
      ssmx = 0.0
      DO mgs = 1,ngscnt
      
      kp1 = Min(nz, kgs(mgs)+1 )
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1) & 
     &                  +w(igs(mgs),jgs,kgs(mgs)))
      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) & 
     &                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))

      ssat0(mgs)  = ssfilt(igs(mgs),jgs,kgs(mgs))
      ssf(mgs)    = ssfilt(igs(mgs),jgs,kgs(mgs))
!      ssmx = Max( ssmx, ssf(mgs) )

      
      ssfkp1(mgs) = ssfilt(igs(mgs),jgs,Min(nz-1,kgs(mgs)+1))
      ssfkm1(mgs) = ssfilt(igs(mgs),jgs,Max(1,kgs(mgs)-1))

!        IF ( wvel(mgs) /= 0.0 ) write(0,*) 'mgs,wvel1,ssf = ',mgs,wvel(mgs),ssf(mgs)


      ENDDO



!
!  cloud water variables
!

      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Set cloud water variables'

      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      IF ( ipconc .ge. 2 .and. cx(mgs,lc) .gt. 1.0e6 ) THEN
        xmas(mgs,lc) = &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
      ELSE
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. cxmin ) THEN
        xmas(mgs,lc) = &
     &     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)), &
     &      xdn(mgs,lc)*xvmx(lc) )

        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)

       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. cxmin ) THEN
!        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
!        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)

       ELSE
        xmas(mgs,lc) = cwmasn
       ENDIF
      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3


      end do
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
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.) ! xdia(mgs,lr,1)
        IF ( imurain == 3 ) THEN
!          xdia(mgs,lr,1) = (6.*pii*xv(mgs,lr)/(alpha(mgs,lr)+1.))**(1./3.)
          xdia(mgs,lr,1) = xdia(mgs,lr,3) ! formulae for Ziegler (1985) use mean volume diameter, not lambda**(-1)
        ELSE ! imurain == 1, Characteristic diameter (1/lambda)
          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
        ENDIF
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
!  Ventilation coefficients

      do mgs = 1,ngscnt


      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* & 
     &  (temg(mgs)/296.0)**(1.5)

      fakvisc(mgs) = fadvisc(mgs)*rhoinv(mgs)

      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)* & 
     &  (101325.0/(pres(mgs)))
      
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))

      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))

      end do
!
!
!  Ziegler nucleation 
!
!
! cloud evaporation, condensation, and nucleation
!  sqsat -> qss(mgs)

      DO mgs=1,ngscnt
        dcloud = 0.0
        ! Skip points at low temperature if SS stays less than 1.08, 
        ! otherwise allow nucleation at low temp (will freeze at next time step)
        IF ( temg(mgs) .le. tfrh .and. qx(mgs,lv)/qvs(mgs) < maxlowtempss ) THEN 
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


      IF ( qx(mgs,lc) <= QEVAP ) THEN ! GO TO 63
        qwvp(mgs) = qwvp(mgs) + qx(mgs,lc)
        thetap(mgs) = thetap(mgs) - felvcp(mgs)*qx(mgs,lc)/(pi0(mgs))
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felvcp(mgs)*qx(mgs,lc)/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - qx(mgs,lc)*rho0(mgs)*dv/dtp  ! evaporation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felvcp(mgs)*qx(mgs,lc)/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - qx(mgs,lc)*rho0(mgs)*dv/dtp  ! evaporation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) - qx(mgs,lc)/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcevap >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcevap) = -qx(mgs,lc)/dtp
        ENDIF
#else
        IF ( io_flag .and. nxtra > 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),1) = -qx(mgs,lc)/dtp
        ENDIF
#endif
        qx(mgs,lc) = 0.
        IF ( restoreccn ) THEN
           IF ( lccna > 1 ) THEN
              ccna(mgs) = ccna(mgs) - restoreccnfrac*cx(mgs,lc)
           ELSEIF ( irenuc <= 2 ) THEN
              IF ( .not. invertccn ) THEN
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*cx(mgs,lc) ) )
              ELSE
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*cx(mgs,lc)
              ENDIF
          ENDIF
        ENDIF
        cx(mgs,lc) = 0.
      ELSE
        qctmp = qx(mgs,lc)
        qwvp(mgs) = qwvp(mgs) + QEVAP
        qx(mgs,lc) = qx(mgs,lc) - QEVAP
        IF ( qx(mgs,lc) .le. 0. ) THEN
          IF ( restoreccn ) THEN
            IF ( lccna > 1 ) THEN
              ccna(mgs) = ccna(mgs) - restoreccnfrac*cx(mgs,lc)
            ELSEIF ( irenuc <= 2 ) THEN
!              ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + cx(mgs,lc) ) )
!              ccnc(mgs) = ccnc(mgs) + cx(mgs,lc)
              IF ( .not. invertccn ) THEN
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*cx(mgs,lc) ) )
              ELSE
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*cx(mgs,lc)
              ENDIF
            ENDIF
          ENDIF
          cx(mgs,lc) = 0.
        ELSE
          tmp = 0.9*QEVAP*cx(mgs,lc)/qctmp ! let droplets get smaller but also remove some. A factor of 1.0 would maintain same size
          IF ( restoreccn ) THEN
            IF ( lccna > 1 ) THEN
              ccna(mgs) = ccna(mgs) - restoreccnfrac*tmp
            ELSEIF ( irenuc <= 2 ) THEN
 !             ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + tmp ) )
!              ccnc(mgs) = ccnc(mgs) + tmp
              IF ( .not. invertccn ) THEN
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*tmp ) )
              ELSE
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*tmp
              ENDIF
            ENDIF
          ENDIF
          cx(mgs,lc) = cx(mgs,lc) - tmp
        ENDIF
        thetap(mgs) = thetap(mgs) - felvcp(mgs)*QEVAP/(pi0(mgs))
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felvcp(mgs)*QEVAP/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - QEVAP*rho0(mgs)*dv/dtp  ! evaporation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) - felvcp(mgs)*QEVAP/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) - QEVAP*rho0(mgs)*dv/dtp  ! evaporation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) - QEVAP/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcevap >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcevap) = -QEVAP/dtp
        ENDIF
#else
        IF ( io_flag .and. nxtra > 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),1) = -QEVAP/dtp
        ENDIF
#endif

      ENDIF

      GO TO 631


  620 CONTINUE

!.... CLOUD CONDENSATION

        IF ( qx(mgs,lc) .GT. qxmin(lc) .and. cx(mgs,lc) .ge. 1. ) THEN



!       ac1 =  xdn(mgs,lc)*elv(kgs(mgs))**2*epsi/
!     :        (tka(kgs(mgs))*rw*temg(mgs)**2)
! took out xdn factor because it cancels later...
       ac1 =  felv(mgs)**2/(tka(mgs)*rw*temg(mgs)**2)


!       bc = xdn(mgs,lc)*rw*temg(mgs)/
!     :       (epsi*wvdf(kgs(mgs))*es(mgs))
! took out xdn factor because it cancels later...
       bc =   rw*temg(mgs)/(wvdf(mgs)*es(mgs))

!       bs = rho0(mgs)*((rd*temg(mgs)/(epsi*es(mgs)))+
!     :             (epsi*elv(kgs(mgs))**2/(pres(mgs)*temg(mgs)*cp)))

!       taus = Min(dtp, xdn(mgs,lc)*rho0(mgs)*(ac1+bc)/
!     :        (4*pi*0.89298*BS*0.5*xdia(mgs,lc,1)*cx(mgs,lc)*xdn(mgs,lc)))

!
      IF ( ssf(mgs) .gt. 0.0 .or. ssat0(mgs) .gt. 0.0 ) THEN
       IF ( ny .le. 2 ) THEN
!        write(0,*)  'undershoot: ',ssf(mgs),
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

       IF ( rcond .eq. 2 .and. qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > 1.e-9 ) THEN
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

           IF ( iferwisventr == 1 ) THEN
             alpr = Min(alpharmax,alpha(mgs,lr) )
!             alpr = alpha(mgs,lr)
             x =  1. + alpr

              tmp = 1 + alpr
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = 2.5 + alpr + 0.5*bx(lr)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions

!         vent1 = dble(xdia(mgs,lr,1))**(-2. - alpr) ! Actually OK
!         vent2 = dble(1./xdia(mgs,lr,1) + 0.5*fx(lr))**dble(2.5+alpr+0.5*bx(lr))  ! Actually OK
         vent1 = dble(xdia(mgs,lr,1))**(0.5 + 0.5*bx(lr)) ! 2016.2.26 Changed for consistency with derivation (recast formula)
         vent2 = dble(1. + 0.5*fx(lr)*xdia(mgs,lr,1))**dble(2.5+alpr+0.5*bx(lr))
        
        
        rwvent(mgs) =    &
     &    0.78*x +    &
     &    0.308*fvent(mgs)*y*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)

           ELSEIF ( iferwisventr == 2 ) THEN
          
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
            x =  1. + alpha(mgs,lr)

            rwvent(mgs) =   &
     &        (0.78*x + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &         *Sqrt((ar*rhovt(mgs)))   &
     &         *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )

          
          ENDIF ! iferwisventr
          
       ENDIF ! imurain

       d1r = (1./(ac1 + bc))*4.0*pi*rwvent(mgs) & 
     &        *0.5*xdia(mgs,lr,1)*cx(mgs,lr)*rhoinv(mgs)
       ELSE
       d1r = 0.0
       ENDIF
       
       
       e1  = felvcp(mgs)/(pi0(mgs))
       f1 = pk(mgs) ! (pres(mgs)/poo)**cap

!
!  fifth trial to see what happens:
!
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       ltemq1 = ltemq
       temp1 = temg(mgs)
       IF ( iqvsopt == 0 ) THEN
         p380 = 380.0/pres(mgs)
       ELSE
         p380 = esbolton*rdorv/(pres(mgs) - es(mgs))
       ENDIF

!       taus = Max( 0.05*dtp, Min(taus, 0.25*dtp ) )
!       nc = NInt(dtp/Min(1.0,0.5*taus))
!       dtcon = dtp/float(nc)
       ss1 = qx(mgs,lv)/qvs(mgs)
       ss2 = ss1
       temp2 = temp1
       qv1 = qx(mgs,lv)
       qvs1 = qvs(mgs)
       qis1 = qis(mgs)
       dt1 = 0.0


!          dtcon = Max(dtcon,0.2)
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

       n = 1
       dt1 = 0.0
       nc = 0
       dqc = 0.0
       dqr = 0.0
       dqi = 0.0
       dqs = 0.0
       dqvii = 0.0
       dqvis = 0.0

#ifdef COMMAS
       lprint = .false.
!       IF ( ixbeg + igs(mgs) == 29 .and. jybeg + jyslab == 39 .and. ss1 > 1.04 ) THEN
!         lprint = .true.
!         write(0,*) 'rk2c start: ', kgs(mgs), qv1, ss1, 100.*(qv1/qvs1-1.0), nc, dtcon1, dtcon2
!         write(0,*) 'cx,cdia = ', 1.e-6*cx(mgs,lc),1.e6*xdia(mgs,lc,1),temp1,thetap(mgs)+theta0(mgs)
!       ENDIF
#endif
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
!          write(0,*) 'RK2c dqv1 = ',dqv
! calculate midpoint values:
     !      ltemq1m = ltemq1 + Nint(dtemp*fqsat + 0.5)

         ! 7.6.2016: Test full calc of ltemq
           ltemq1m = (temp1+dtemp-163.15)*fqsati+1.5
           ltemq1m = Min( nqsat, Max(1,ltemq1m) )

           IF ( ltemq1m .lt. 1 .or. ltemq1m .gt. nqsat ) THEN
             write(0,*) 'STOP in nucond line 1192 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,f1,dqv,dqvr = ', dtemp,e1,f1,dqv,dqvr
             write(0,*) ' d1,d1r,dtcon,ss1 = ',d1,d1r,dtcon,ss1
             write(0,*) ' dqc, dqr = ',dqc,dqr
             write(0,*) ' qv,qc,qr = ',qx(mgs,lv)*1000.,qx(mgs,lc)*1000.,qx(mgs,lr)*1000.
             write(0,*) ' i, j, k = ',igs(mgs),jy,kgs(mgs)
             write(0,*) ' dtcon1,dtcon2,delta = ',dtcon1,dtcon2,delta
             write(0,*) ' nc,dtp = ',nc,dtp
             write(0,*) ' rwvent,xdia,crw,ccw = ', rwvent(mgs),xdia(mgs,lr,1),cx(mgs,lr),cx(mgs,lc)
             write(0,*) ' fvent,alphar = ',fvent(mgs),alpha(mgs,lr)
             write(0,*) ' xvr,xmasr,xdnr,cwc1 = ',xv(mgs,lr),xmas(mgs,lr),xdn(mgs,lr),cwc1
           ENDIF
            dqvs = dtemp*p380*dtabqvs(ltemq1m)
            qv1m = qv1 + dqv + dqvr
!          qv1mr = qv1r + dqvr

            qvs1m = qvs1 + dqvs
            ss1m = qv1m/qvs1m

    ! check for undersaturation when no ice is present, if so, then reduce time step
          IF ( ss1m .lt. 1.  .and. (dqvii + dqvis) .eq. 0.0 ) THEN
            dtcon = (0.5*dtcon)
            IF ( dtcon .ge. dtcon1 ) THEN
             GOTO 609
            ELSE
             EXIT
            ENDIF
          ENDIF
! calculate full step:
          dqv  = -(ss1m - 1.)*d1*dtcon
          dqvr = -(ss1m - 1.)*d1r*dtcon


!          write(0,*) 'RK2a dqv1m = ',dqv
          dtemp = -e1*f1*(dqv + dqvr)
          
         ! ltemq1 = ltemq1 + Nint(dtemp*fqsat + 0.5)

         ! 7.6.2016: Test full calc of ltemq
           ltemq1 = (temp1+dtemp-163.15)*fqsati+1.5
           ltemq1 = Min( nqsat, Max(1,ltemq1) )

           IF ( ltemq1 .lt. 1 .or. ltemq1 .gt. nqsat ) THEN
             write(0,*) 'STOP in nucond line 1230 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,dqv,dqvr = ', dtemp,e1,dqv,dqvr
           ENDIF
          dqvs = dtemp*p380*dtabqvs(ltemq1)

          qv1 = qv1 + dqv + dqvr

          dqc = dqc - dqv
          dqr = dqr - dqvr

          qvs1 = qvs1 + dqvs
          ss1 = qv1/qvs1
          temp1 = temp1 + dtemp
          IF ( temp2 .eq. temp1 .or. ss2 .eq. ss1 .or.  &
     &           ss1 .eq. 1.00 .or.  &
     &      ( n .gt. 10 .and. ss1 .lt. 1.0005 ) ) THEN
!           write(0,*) 'RK2c break'
           EXIT
          ELSE
           ss2 = ss1
           temp2 = temp1
           dt1 = dt1 + dtcon
           n = n + 1
          ENDIF
       ENDDO RK2c


        dcloud = dqc ! qx(mgs,lv) - qv1
        thetap(mgs) = thetap(mgs) + e1*(DCLOUD + dqr)

#ifdef COMMAS
!       IF ( lprint ) THEN

!         temp1 = (theta0(mgs)+thetap(mgs))*t77(ix,jy,kz)
!          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
!         ltemq = Min( nqsat, Max(1,ltemq) )
!
!          c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
!
!          ssf(mgs) = 0.0
!          IF ( c1 > 0. ) THEN
!            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
!          ENDIF

!         write(0,*) 'rk2c end: ', qv1, ss1, dqc, e1*(DCLOUD + dqr), temp1,thetap(mgs)+theta0(mgs),  &
!                (thetap(mgs)+theta0(mgs))*t77(igs(mgs),jy,kgs(mgs))
!       ENDIF
#endif

        IF ( eqtset > 2 ) THEN
           pipert(mgs) = pipert(mgs) + felvpi(mgs)*(DCLOUD + dqr)
        ENDIF
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + e1*(DCLOUD + dqr)*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + (DCLOUD + dqr)*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + e1*(DCLOUD + dqr)*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + (DCLOUD + dqr)*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + &
                  (DCLOUD + dqr)/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = DCLOUD/dtp
        ENDIF
#else
        IF ( io_flag .and. nxtra > 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),1) = DCLOUD/dtp
           axtra(igs(mgs),jy,kgs(mgs),2) = axtra(igs(mgs),jy,kgs(mgs),2) + dqr/dtp
        ENDIF
#endif
        qwvp(mgs) = qwvp(mgs) - (DCLOUD + dqr)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        qx(mgs,lr) = qx(mgs,lr) + dqr
!        t9(igs(mgs),jy,kgs(mgs)) = t9(igs(mgs),jy,kgs(mgs)) + (DCLOUD + dqr)/dtp*felv(mgs)/(cp*pi0(mgs)) !* &
!!     &                 dx*dy*dz3d(igs(mgs),jy,kgs(mgs))

#ifdef Z3MOM

        IF ( lzr > 1 .and. rcond == 2 .and. qx(mgs,lr) .gt. qxmin(lr)   &
     &       .and. cx(mgs,lr) .gt. 1.e-9 ) THEN
          tmp = qx(mgs,lr)/cx(mgs,lr)
          IF ( imurain == 3 ) THEN
          g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
          ELSE
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
          
          ENDIF
          zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(xdn(mgs,lr)))**2*( 2.*( tmp ) * dqr )
        ENDIF
#endif

        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*f1
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        IF ( iqvsopt == 0 ) THEN
          qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        ELSEIF ( iqvsopt == 1 ) THEN
          qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
        ELSE
          tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
          qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
#endif
        ENDIF
!        es(mgs) = 6.1078e2*tabqvs(ltemq)

!

      ENDIF  ! dcloud .gt. 0.


      ELSE  ! qc .le. qxmin(lc)

!        IF ( ssf(mgs) .gt. 0.0 .and. .not. flag_qndrop ) THEN ! flag_qndrop turns off primary nucleation when using wrf-chem with progn=1
        IF ( ssf(mgs) .gt. 0.0 ) THEN ! .and.  ssmax(mgs) .lt. sscb ) THEN ! except that wrf-chem does not seem to initialize qc for activated aerosols, so keep this, after all

          IF ( iqcinit == 1 ) THEN

         qvs0   = 380.*exp(17.27*(temg(mgs)-273.)/(temg(mgs)- 36.))/pk(mgs)

         dcloud = Max(0.0, (qx(mgs,lv)-qvs0) / (1.+qvs0*f5/(temg(mgs)-36.)**2) )

          ELSEIF ( iqcinit == 3 ) THEN
              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felvcp(mgs)/ & 
     &             ((temg(mgs) - cbw)**2))
            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
          
          ELSEIF ( iqcinit == 2 ) THEN
!              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/
!     :             (cp*(temg(mgs) - cbw)**2))
!            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
         ssmx = ssmxinit

!          IF ( ssf(mgs) > ssmx  .and. ssmax(mgs) < 3.0 ) THEN
!          IF ( ssf(mgs) > ssmx  .and. ccnc(mgs) > 1.0 ) THEN
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 5.0 .and. ccnc(mgs) > 0.1*cwnccn(mgs) ) THEN ! this one works
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 ) THEN ! test -- fails
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 .and. ccnc(mgs) > 0.1*cwnccn(mgs)) THEN ! test -- is OK
          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 .and.  &
             ( ccnc(mgs) > 0.05*cwnccn(mgs) .or. ( ac_opt > 0 .and. ccnc_ac(mgs) - cx(mgs,lc) > 0.0 ) ) ) THEN ! test
!          IF ( ssf(mgs) > ssmx ) THEN ! original condition
           CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,dcloud, & 
     &      pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felvcp,ssmx,pk,ngscnt)
          ELSE
            dcloud = 0.0
          ENDIF
         ENDIF
        ELSE
            dcloud = 0.0
        ENDIF
        
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = DCLOUD/dtp
        ENDIF
#else
        IF ( io_flag .and. nxtra > 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),1) = DCLOUD/dtp
        ENDIF
#endif
        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) !( pres(mgs) / poo ) ** cap
!        temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
       ! qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        IF ( iqvsopt == 0 ) THEN
          qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        ELSEIF ( iqvsopt == 1 ) THEN
          qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
        ELSE
          tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
          qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
#endif
        ENDIF
!        es(mgs) = 6.1078e2*tabqvs(ltemq)

!.... S. TWOMEY (1959)
! Note: get here if there is no previous cloud water and w > 0.
      cn(mgs) = 0.0
      
      IF ( ncdebug .ge. 1 ) THEN
        write(iunit,*) 'at 613: ',qx(mgs,lc),cx(mgs,lc),wvel(mgs),ssmax(mgs),kgs(mgs)
      ENDIF
      
      IF (  .not. flag_qndrop ) THEN ! { do not calculate number of droplets if using wrf-chem

      IF ( ac_opt == 0 ) THEN
        cnuctmp = cnuc(mgs)
      ELSE
        cnuctmp = ccnc_ac(mgs)
      ENDIF
      
!      IF ( ssmax(mgs) .lt. sscb .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
      IF ( dcloud .gt. qxmin(lc) .and. wvel(mgs) > 0.0) THEN
!       CN(mgs) =   CCNE*wvel(mgs)**cnexp ! *Min(1.0,1./dtp) ! 0.3465
       CN(mgs) =   CCNE0*cnuctmp**(2./(2.+cck))*wvel(mgs)**cnexp ! *Min(1.0,1./dtp) ! 0.3465
        IF ( ny .le. 2 .and. cn(mgs) .gt. 0.0    &
     &                    .and. ncdebug .ge. 1 ) THEN 
          write(iunit,*) 'CN: ',cn(mgs)*1.e-6, cx(mgs,lc)*1.e-6, qx(mgs,lc)*1.e3,   &
     &       wvel(mgs), dcloud*1.e3
          IF ( cn(mgs) .gt. 1.0 ) write(iunit,*) 'cwrad = ',   &
     &       1.e6*(rho0(mgs)*qx(mgs,lc)/cn(mgs)*cwc1)**c1f3,   &
     &   igs(mgs),kgs(mgs),temcg(mgs),    &
     &   1.e3*an(igs(mgs),jgs,kgs(mgs)-1,lc)
        ENDIF
        IF ( iccwflg .eq. 1 ) THEN
          cn(mgs) = Min(cwccn*rho0(mgs)/rho00, Max(cn(mgs),   &
     &       rho0(mgs)*qx(mgs,lc)/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3)))
        ENDIF
      ELSE
       cn(mgs) = 0.0
       dcloud = 0.0
!          cn(mgs) = Min(cwccn,    &
!     &       rho0(mgs)*dcloud/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3) )
      ENDIF

      IF ( cn(mgs) .gt. 0.0 ) THEN
       IF ( ac_opt == 0 ) THEN
         IF ( cn(mgs) .gt. ccnc(mgs) ) THEN
           cn(mgs) = ccnc(mgs)
!          ccnc(mgs) = 0.0
         ENDIF
       ELSE 
         cn(mgs) = Min( cn(mgs), ccnc_ac(mgs) )
       ENDIF
!      cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       IF ( irenuc <= 2 .and. lccna < 1  ) ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
       ccna(mgs) = ccna(mgs) + cn(mgs)
      ENDIF

!       write(91,*) 'nuc1: cn, ix, kz = ',cn(mgs),igs(mgs),kgs(mgs),wvel(mgs),cnexp,ccnc(mgs)

      IF( CN(mgs) .GT. cx(mgs,lc) ) cx(mgs,lc) = CN(mgs)
      IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .le. qxmin(lc) ) THEN
        cx(mgs,lc) = 0.
      ELSE
        cx(mgs,lc) = Min(cx(mgs,lc),rho0(mgs)*Max(0.0,qx(mgs,lc))/cwmasn)
      ENDIF
      
      ENDIF ! }.not. flag_qndrop

        GOTO 613
        
        END IF ! qc .gt. 0.

!        ES=EES(PIB(K)*PT)
!        SQSAT=EPSI*ES/(PB(K)*1000.-ES)

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

      if (ndebug .gt. 0) write(0,*) "ICEZVD_DR: Entered Ziegler Cloud Nucleation" !mpidebug

      DSSDZ=0.
#ifdef COMMAS
      r2dzm=0.50/dz1d(kgs(mgs))
#else
      r2dzm=0.50/dz3d(igs(mgs),jy,kgs(mgs))
#endif

      IF ( irenuc >= 0 .and. ac_opt == 0 .and. .not. flag_qndrop ) THEN ! turn off nucleation when flag_qndrop (using WRF-CHEM for activation)

      IF ( irenuc < 2 ) THEN !{

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
!        IF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1 & !)
     &  .and. ssf(mgs) .gt. 0.0  .and. wvel(mgs) .gt. 0.0 &
     &  .and. ssfkp1(mgs) .gt. 0.0   &
     &  .AND. ssfkm1(mgs) .le. 0.0 .and. wvelkm1(mgs) .gt. 0.0 &
     &  .AND. ssf(mgs) .gt. ssfkm1(mgs)  &
     &  .and. t0p1 .gt. 233.2) THEN
         DSSDZ = 2.*(ssf(mgs) - ssfkm1(mgs))*R2DZM  ! 1-sided difference
        ENDIF

       ENDIF
!
!CLZ  IF(wijk.LE.0.) CN=CCN*ssfilt(ix,jy,kz)**CCK
! note: CCN -> cwccn, DELT -> dtp
      c1 = Max(0.0, rho0(mgs)*(qx(mgs,lv) - qss(mgs))/ &
     &        (xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3))
      IF ( lccn .lt. 1 ) THEN
       CN(mgs) = cwccn*rho0(mgs)/rho00*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         (wvel(mgs)*DSSDZ) )      ! probably the vertical gradient dominates
      ELSE
      CN(mgs) =  &
     &    Min(ccnc(mgs), cnuc(mgs)*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         ( wvel(mgs)*DSSDZ) )  )
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

      ELSEIF ( irenuc == 2 ) THEN !} { 
      ! simple Twomey scheme
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       CN(mgs) =   CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp ! *Min(1.0,1./dtp) ! 0.3465
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       CN(mgs) = Min( CN(mgs), Max(0.0, (cnuc(mgs) - ccna(mgs) )) )
       
        IF ( .false. .and. ny <= 2 ) THEN
          write(0,*) 'i,k, cwmasn = ',igs(mgs),kgs(mgs),cwmasn
          write(0,*) 'wvel, cnuc, cn = ',wvel(mgs),cnuc(mgs),cn(mgs)
          write(0,*) 'ccne0,cnexp,cck = ',ccne0,cnexp,cck
          write(0,*) 'part1, part2 = ',CCNE0*cnuc(mgs)**(2./(2.+cck)), Max(0.0,wvel(mgs))**cnexp
          write(0,*) 'ccnc, dqc, dqc/cwmasn = ',ccnc(mgs), dqc, 0.5*dqc/cwmasn
        ENDIF
       
       IF ( icnuclimit > 0 ) THEN 
#ifdef XCOMMENTX
! max droplet conc. based on Chandrakar et al. (2016) and Konwar et al. (2012)
! see chandrakar2016.nb
#endif
         tmp = ccnc(mgs) + cx(mgs,lc)
         IF ( tmp < 330.34e6 ) THEN
           ccwmax = 1.1173e6 * (1.e-6*tmp)**0.9504
         ELSE
           ccwmax = 21.57e6 * (1.e-6*tmp)**0.44
         ENDIF
         
!         IF ( cn(mgs) > 0. ) THEN
!          write(0,*) 'cn,tmp,ccwmax,cx,c-cx = ',cn(mgs),tmp,ccwmax,cx(mgs,lc),ccwmax - cx(mgs,lc) 
!         ENDIF
         
        cn(mgs) = Max( 0.0, Min( cn(mgs), ccwmax - cx(mgs,lc) ) )
       
       ENDIF
       
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       IF ( lccna < 1 ) ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))

#ifdef Z3MOM
      ELSEIF ( irenuc == 3 ) THEN !} { 
      ! Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!       CN(mgs) =   cwccn*Min(ssf(mgs),ssfcut)**cck 

! Need to calculate new ssf since condensation has happened:
         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1= pqs(mgs)*tabqvs(ltemq)

          ssf(mgs) = 0.0
          IF ( c1 > 0. ) THEN
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          ENDIF
       CN(mgs) =   cnuc(mgs)*Min(1.0, (ssf(mgs))**cck ) ! 

       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
        ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs)) 
       
      ELSEIF ( irenuc == 4 ) THEN !} { 
      ! modification of Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!       CN(mgs) =   CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp
!       cn(mgs) = Min( cn(mgs), cnuc(mgs) )
! Need to calculate new ssf since condensation has happened:
         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1= pqs(mgs)*tabqvs(ltemq)
          IF ( c1 > 0. ) THEN
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          ELSE
            ssf(mgs) = 0.0
          ENDIF
       CN(mgs) =   cnuc(mgs)*Min(ssf2kmax, ssf(mgs)**cck) ! this allows cn(mgs) > cnuc(mgs)

       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       IF ( cn(mgs) > 0.0 ) THEN
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs)) 
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
        ENDIF
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
!        ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
       


      ELSEIF ( irenuc == 6 ) THEN !} { 

      ! simple Twomey scheme but limit activation to try to do most activation near cloud base, but keep some CCN available for renuclation
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       cn(mgs) = 0.0
!       IF ( ccna(mgs) < 0.7*cnuc(mgs) .and. ccnc(mgs) > 0.69*cnuc(mgs) - ccna(mgs)) THEN ! here, assume we are near cloud base and use Twomey formulation
       IF ( ccna(mgs) < 0.7*cnuc(mgs) ) THEN ! here, assume we are near cloud base and use Twomey formulation
         CN(mgs) =  Min( 0.9*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
!         IF ( cn(mgs) + ccna(mgs) > 0.71*cnuc ) THEN
         ! prevent this branch from activating more than 70% of CCN
#ifdef COMMAS
           IF ( ny == 2 .and. igs(mgs) == 29 .and. kgs(mgs) < 14 ) THEN
!             write(0,*) 'iren6: cn,cnuc,ccnc,ccna,max =',igs(mgs),kgs(mgs),cn(mgs),cnuc(mgs),ccnc(mgs),ccna(mgs), &
!                 Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ),0.71*cnuc(mgs) - ccna(mgs),ccna(mgs)/cnuc(mgs)
           ENDIF
#endif
           CN(mgs) = Min( CN(mgs), Max(0.0, (0.7*cnuc(mgs) - ccna(mgs) )) )
!           CN(mgs) = Min( CN(mgs), Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ) )
           
       ELSE
        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

!          c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          c1= pqs(mgs)*tabqvs(ltemq)
          IF ( c1 > 0. ) THEN
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          ELSE
            ssf(mgs) = 0.0
          ENDIF

!        CN(mgs) = cnuc(mgs)*Min(0.99, Min(ssf(mgs),ssfcut)**cck ) ! 
         CN(mgs) =   cnuc(mgs)*Min(2.0, Max(0.0,ssf(mgs))**cck ) ! 
!         CN(mgs) =   cnuc(mgs)*Min(ssf(mgs),ssfcut)**cck ! 


#ifdef COMMAS
!           IF ( ny == 2 .and. igs(mgs) == 29 .and. kgs(mgs) < 14 ) THEN
           IF ( (ny == 2 .and. ssf(mgs) > 0.5 ) .or. ssf(mgs) > 0.5 ) THEN
!         write(0,*) 'iren6b: ',igs(mgs),kgs(mgs),cn(mgs),ssf(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ,cnuc(mgs),ccnc(mgs),ccna(mgs)
           ENDIF
#endif
        CN(mgs) = Min(0.01*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
!        cn(mgs) = 0.0
       ENDIF
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       IF ( cn(mgs) > 0.0 ) THEN
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
        ENDIF
#endif /* Z3MOM */
      ELSEIF ( irenuc == 5 ) THEN !} { 

      ! modification of Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!      CN(mgs) =  Min( 0.91*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
       CN(mgs) =  Min( cnuc(mgs),  CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )

         
        IF ( ccna(mgs) >= cnuc(mgs) ) THEN ! apply limit after all "base" CCN have been depleted
        temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          IF ( iqvsopt == 0 ) THEN
            c1 = pqs(mgs)*tabqvs(ltemq)
          ELSEIF ( iqvsopt == 1 ) THEN
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
          ELSE
            tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
            c1 = rdorv*tmp/(pres(mgs) - tmp)
#endif
          ENDIF
          IF ( c1 > 0. ) THEN
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          ELSE
            ssf(mgs) = 0.0
          ENDIF
          

       CN(mgs) =   Max( cn(mgs), cnuc(mgs)*Min(ssf2kmax, ssf(mgs)**cck) ) ! this allows cn(mgs) > cnuc(mgs)

   !    cn(mgs) = Min( cn(mgs), cnuc(mgs) )

!       IF ( ccna(mgs) >= cnuc(mgs) ) THEN ! apply limit after all "base" CCN have been depleted
       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
       
       ELSE
         CN(mgs) =  Min( cn(mgs), cnuc(mgs) - ccna(mgs) ) ! no more than remaining "base" CCN
       ENDIF
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       dcrit = 2.0*2.0e-6
       dcloud = 1000.*dcrit**3*Pi/6.
 !      cn(mgs) = Min(cn(mgs), 0.5*dqc/dcloud) ! limit the nucleation mass to half of the condensation mass
       ! check new droplet size:
         ! tmp is number of droplets at diameter dcrit
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         cn(mgs) = Min(tmp, cn(mgs) )

      
       IF ( cn(mgs) > 0.0 ) THEN
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       dcrit = 2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
        ENDIF
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      ELSEIF ( irenuc == 7 .or. irenuc == 17 ) THEN !} { 

      ! simple Twomey scheme but limit activation to try to do most activation near cloud base, but keep some CCN available for renuclation
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       cn(mgs) = 0.0
       IF ( irenuc == 7 ) THEN
         frac = 0.9
       ELSE
         frac = 0.98
       ENDIF
!       IF ( ccna(mgs) < 0.7*cnuc(mgs) .and. ccnc(mgs) > 0.69*cnuc(mgs) - ccna(mgs)) THEN ! here, assume we are near cloud base and use Twomey formulation
       IF ( ccna(mgs) < frac*cnuc(mgs) ) THEN ! { here, assume we are near cloud base and use Twomey formulation
         CN(mgs) =  Min( (frac+0.01)*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
!         IF ( cn(mgs) + ccna(mgs) > 0.71*cnuc ) THEN
         ! prevent this branch from activating more than 70% of CCN
#ifdef COMMAS
           IF ( ny == 2 .and. igs(mgs) == 29 .and. kgs(mgs) < 14 ) THEN
!             write(0,*) 'iren6: cn,cnuc,ccnc,ccna,max =',igs(mgs),kgs(mgs),cn(mgs),cnuc(mgs),ccnc(mgs),ccna(mgs), &
!                 Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ),0.71*cnuc(mgs) - ccna(mgs),ccna(mgs)/cnuc(mgs)
           ENDIF
#endif
           CN(mgs) = Min( CN(mgs), Max(0.0, (frac*cnuc(mgs) - ccna(mgs) )) )
!           CN(mgs) = Min( CN(mgs), Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ) )
         !  write(0,*) '1: k,cn = ',kgs(mgs),cn(mgs),ssf(mgs)
!!           IF ( ccncuf(mgs) > 0.0 .and. cn(mgs) < 1.e-3 .and. ssmax(mgs) > 1.0 ) THEN
!           IF ( ccncuf(mgs) > 0.0 .and. ssf(mgs) > ssmxuf .and. ssmax(mgs) > ssmxuf ) THEN
!            CNuf(mgs) =  Min( ccncuf(mgs), CCNE0*ccncuf(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
          !  IF ( cnuf(mgs) >= 0.0 ) write(0,*) '1: cnuf, k = ',cnuf(mgs),ccncuf(mgs),kgs(mgs)
!           ENDIF

           
       ELSE ! }{
        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

        !  c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          IF ( iqvsopt == 0 ) THEN
            c1 = pqs(mgs)*tabqvs(ltemq)
          ELSEIF ( iqvsopt == 1 ) THEN
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
          ELSE
            tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
            c1 = rdorv*tmp/(pres(mgs) - tmp)
#endif
          ENDIF

          ssf(mgs) = 0.0
          IF ( c1 > 0. ) THEN
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          ENDIF

!          IF ( ssf(mgs) <= 1.0 .or. cnuc(mgs) > ccna(mgs) ) THEN
          IF ( ssf(mgs) <= 1.0 ) THEN
          CN(mgs) =   cnuc(mgs)*Min(1.0, Max(0.0,ssf(mgs))**cck ) ! 
          ELSE
          CN(mgs) =   cnuc(mgs)*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) !           
!          write(0,*) 'iren7: ssf,ssmx = ',ssf(mgs),ssmax(mgs),cn(mgs),ccna(mgs),cnuc(mgs)
!          write(0,*) 'c1,qv = ',c1,qx(mgs,lv),temp1,ltemq
          ENDIF

         !  write(0,*) 'k,cn = ',kgs(mgs),cn(mgs),ssf(mgs)
         !  write(0,*) 'ccn-ccna = ',cnuc(mgs) - ccna(mgs),ccnc(mgs) - ccna(mgs)
!           IF ( ccncuf(mgs) > 0.0 .and. cn(mgs) < 1.e-3 .and. ssmax(mgs) > 1.0 ) THEN
           IF ( ccncuf(mgs) > 0.0 .and. ssf(mgs) > ssmxuf .and. ( ssmax(mgs) > ssmxuf .or. lss < 1 ) ) THEN
            CNuf(mgs) =  Min( ccncuf(mgs), CCNE0*ccncuf(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
          !  IF ( cnuf(mgs) >= 0.0 ) write(0,*) 'cnuf, k = ',cnuf(mgs),ccncuf(mgs),kgs(mgs)
           ENDIF
          

#ifdef COMMAS
!           IF ( ny == 2 .and. igs(mgs) == 29 .and. kgs(mgs) < 14 ) THEN
           IF ( (ny == 2 .and. ssf(mgs) > 0.5 ) .or. ssf(mgs) > 0.5 ) THEN
!         write(0,*) 'iren6b: ',igs(mgs),kgs(mgs),cn(mgs),ssf(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ,cnuc(mgs),ccnc(mgs),ccna(mgs)
           ENDIF
#endif
!        CN(mgs) = Min( Min(0.1,ssf(mgs)-1.)*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
!        CN(mgs) = Min( Min(0.5*cx(mgs,lc), Min(0.1,ssf(mgs)/100.)*cnuc(mgs)), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
        
        CN(mgs) = Min(0.01*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from

       ENDIF ! }
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       

        IF ( icnuclimit > 0 ) THEN
! max droplet conc. based on Chandrakar et al. (2016) and Konwar et al. (2012)
#ifdef XCOMMENTX
! see chandrakar2016.nb
#endif
           tmp = ccnc(mgs) - ccna(mgs) + cx(mgs,lc)
           IF ( tmp < 330.34e6 ) THEN
             ccwmax = 1.1173e6 * (1.e-6*tmp)**0.9504
           ELSE
             ccwmax = 21.57e6 * (1.e-6*tmp)**0.44
           ENDIF
          
           cn(mgs) = Max( 0.0, Min( cn(mgs), ccwmax - cx(mgs,lc) ) )
           
        ENDIF

       IF ( cn(mgs) + cnuf(mgs) > 0.0 ) THEN

       dcrit = 2.0*2.0e-6
       dcloud = 1000.*dcrit**3*Pi/6.
 !      cn(mgs) = Min(cn(mgs), 0.5*dqc/dcloud) ! limit the nucleation mass to half of the condensation mass
       ! check new droplet size:
         ! tmp is number of droplets at diameter dcrit
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         cn(mgs) = Min(tmp, cn(mgs) )

        cx(mgs,lc) = cx(mgs,lc) + cn(mgs) + cnuf(mgs)
 
 
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       
       dcrit = 2.0*2.5e-7
       dcloud = 1000.*dcrit**3*Pi/6.*(cn(mgs) + cnuf(mgs) )
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
         ccncuf(mgs) = Max(0.0, ccncuf(mgs) - cnuf(mgs))
        ENDIF

      ELSEIF ( irenuc == 8 ) THEN !} { 
      ! simple Twomey scheme
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       
       cn(mgs) = 0.0
       
       IF ( ccnc(mgs) > 0. ) THEN
       CN(mgs) =   CCNE0*ccnc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp ! *Min(1.0,1./dtp) ! 0.3465
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       
       ELSEIF ( cx(mgs,lc) < 0.01e9 ) THEN

        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

        !  c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          IF ( iqvsopt == 0 ) THEN
            c1 = pqs(mgs)*tabqvs(ltemq)
          ELSEIF ( iqvsopt == 1 ) THEN
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
#if defined(COMMAS) || defined(COMMASTMP)
          ELSE
            tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
            c1 = rdorv*tmp/(pres(mgs) - tmp)
#endif
          ENDIF

          ssf(mgs) = 0.0
          IF ( c1 > 0. ) THEN
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          ENDIF

!          IF ( ssf(mgs) <= 1.0 .or. cnuc(mgs) > ccna(mgs) ) THEN
          IF ( ssf(mgs) <= 1.0 ) THEN
          CN(mgs) = 0.0
          ELSE
!           CN(mgs) = 0.01e9*rho0(mgs)/rho00*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) - cx(mgs,lc) !           
           CN(mgs) = 0.01e9*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) - cx(mgs,lc) !           
          ENDIF
       
       ENDIF

       IF ( cn(mgs) > 0.0 ) THEN
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
       
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
        ENDIF
       

      ELSEIF ( irenuc == 9 .or. irenuc == 10 ) THEN ! } {

#ifdef NUWRFMODS
           ! following Abdul-Razzak et al. (1998, 'part 1') and Abdul-Razzak and Ghan (2000, 'part 2') for
           ! estimated maximum SS and just one aerosol mode (accumulation)
          temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
          ltemq = Min( nqsat, Max(1,ltemq) )

          IF ( iqvsopt == 0 ) THEN
            c1 = pqs(mgs)*tabqvs(ltemq)
          ELSE !  ( iqvsopt == 1 ) THEN
            c1 = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
          ENDIF

          ssf(mgs) = 0.0
          IF ( c1 > 0. ) THEN
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          ENDIF

          IF ( ssf(mgs) > 0.0 .and. wvel(mgs) > 0.0 ) THEN

           cpm = cp + cpv*qx(mgs,lv)
           evs = pres(mgs)*qvs(mgs)/0.622
           alpha_ar = gr*mwwater*felv(mgs)/(cpm*gasconst*temg(mgs)**2)-gr*Mair/(gasconst*temg(mgs))  ! AR98 eq. 11 (alpha)
           gamma_ar = gasconst*temg(mgs)/(evs*mwwater) + mwwater*felv(mgs)**2/(cpm*pres(mgs)*Mair*temg(mgs)) ! AR98 eq. 12 (gamma)

           G_ar = 1./(rhowater*gasconst*temg(mgs)/(evs*wvdf(mgs)*mwwater)+ felv(mgs)*rhowater/(tka(mgs)*temg(mgs))*(felv(mgs)*mwwater/ &
              (temg(mgs)*gasconst)-1.)) ! AR98 eq. 16

!            write(91,*) 'sstmp,ssf,cn1,wvel = ',sstmp,ssf(mgs),cn1,wvel(mgs)
!            sswater = Min(sstmp, ssf(mgs) * 0.01)  ! unit change from percentage to n/a
            sswater = ssf(mgs) * 0.01  ! unit change from percentage to n/a
            sigvl = 0.0761 - 1.55E-4*temcg(mgs)
            aact  = 2.0*mwwater/(rhowater*gasconst*temg(mgs))*sigvl  ! AR98 eq. 5 (A), where sigvl = tau = surface tension (wrongly defined as time in AR98)
            ! could add a subgrid TKE-based W-pert to wvel at some point
            zeta = 2./3.*(alpha_ar*wvel(mgs)/G_ar)**0.5*aact ! wvel is vertical velocity ! ARG2000 eq. 10

              f_ac = 0.5*Exp(2.5*(ac_pgw)**2)     ! ARG2000 eq. 7 (f_i), but size distribution in Cheng et al. uses sigma instead of log(sigma)
              g_ac = 1.0 + 0.25*ac_pgw           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
              eta_ac = (alpha_ar*wvel(mgs)/G_ar)**1.5/(2.*pi*rhowater*gamma_ar*ccnc(mgs)) ! na zero ! ARG2000 eq. 11

              sm_ac = 2.0 / ac_kappa**0.5 * (aact / (3.0*ac_pmr) )**1.5    ! AR98 eq. 9 (S_mi)
              ss_ac = 1./sm_ac**2*(f_ac*(zeta/eta_ac)**1.5 + g_ac*(sm_ac**2/(eta_ac + 3.*zeta))**0.75)  ! ARG2000 eq. 6 (1 term in sum); comes from AR98 eq. 31
            
            
            IF  ( sswater > 0.0 .or. ss_ac > 0.0 )  THEN  ! 
!              F11 = 0.5*EXP(2.5*(LOG(SIG1))**2) ! ARG2000 eq. 7 (f_i)
!              F21 = 1.+0.25*LOG(SIG1)           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
!              ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1) ! na zero ! ARG2000 eq. 11
!     sm_ac ->    SM1  = 2./BACT1**0.5*(AACT/(3.*RM1))**1.5  ! AR98 eq. 8 (or 7?) and ARG2000 eq. 9 for particle critical supersaturation
!              DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)  ! ARG2000 eq. 6 (1 term in sum); comes from AR98 eq. 31
!             SMAX = 1./(DUM1 + ... )**0.5
!      uu_ac ->      UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
!    NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
!             DUM1 = NANEW1/2.*(1.-DERF1(UU1))
!  BACT -> kappa
!  RM1 -> pmr

!-              sm_ac = 2.0 / ac_kappa**0.5 * (aact / (3.0*ac_pmr) )**1.5
              
              IF ( ss_ac > 0.0 ) THEN
               IF ( .true. ) THEN
                 smax = 1./sqrt(ss_ac)
                 uu_ac = 2.0 * log(sm_ac / smax) / (3.0 * 2**0.5 * ac_pgw )  ! AR98 eq. 15 (u_i)
               ELSE
                 uu_ac = 2.0 * log(sm_ac / sswater) / (3.0 * 2**0.5 * ac_pgw)
               ENDIF
!              cn(mgs) = Max( cn1, ccnc_ac(mgs) * 0.5 * (1.0 - Derf(uu_ac))  ) ! number concentration
              cn(mgs) = Max( 0.0, ccnc(mgs) * 0.5 * (1.0 - Derf(uu_ac))  ) ! number concentration of newly-activated droplets
              
              ENDIF
!              IF ( sswater > 0.0 ) &
!              write(91,*) 'cn,c_ac,cn-cx,uuac,sswater = ',cn(mgs),ccnc_ac(mgs),cn(mgs)- cx(mgs,lc),uu_ac,sswater
              IF ( irenuc == 10 ) THEN
                cn(mgs) = Min( cn(mgs) - ccna(mgs), ccnc(mgs) - ccna(mgs) )
              ENDIF
             ENDIF
            ELSE
              cn(mgs) = 0.0
            ENDIF

       IF ( cn(mgs)  > 0.0 ) THEN

           CN(mgs) = Min( CN(mgs), Max(0.0, (cnuc(mgs) - ccna(mgs) )) )

         dcrit = 2.0*2.0e-6
         dcloud = 1000.*dcrit**3*Pi/6.
       ! check new droplet size:
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         ! tmp is number of droplets at diameter dcrit. If cn is too large and makes avg droplet diam too small, then replace with tmp
          cn(mgs) = Min(tmp, cn(mgs) )

          cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
 
 
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       
         dcrit = 2.0*2.5e-7
         dcloud = 1000.*dcrit**3*Pi/6.*(cn(mgs) )
         qx(mgs,lc) = qx(mgs,lc) + DCLOUD
         thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
         qwvp(mgs) = qwvp(mgs) - DCLOUD
       ENDIF
#else
         write(0,*) 'irenuc=9 requires nuwrfmods=1'
#endif
      ENDIF ! }
      

      ccna(mgs) = ccna(mgs) + cn(mgs)



      ENDIF ! irenuc >= 0 .and. .not. flag_qndrop

#ifdef NUWRFMODS
      IF ( ac_opt > 0 ) THEN
!          IF ( dcloud > qxmin(lc) .and. wvel(mgs) > 0.0 ) THEN
!          IF ( ssf(mgs) > 0.0 .and. wvel(mgs) > 0.0 ) THEN

         !  write(0,*) 'nucond: ac_opt = ',ac_opt
          temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
          ltemq = Min( nqsat, Max(1,ltemq) )
          c1= pqs(mgs)*tabqvs(ltemq)

          ssf(mgs) = 0.0
          IF ( c1 > 0. ) THEN
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          ENDIF

          IF ( ssf(mgs) > 0.0 .and. wvel(mgs) > 0.0 ) THEN

!      ccnefac =  (1.63/(cck * beta(3./2., cck/2.)))**(cck/(cck + 2.0))
!      cnexp   = (3./2.)*cck/(cck+2.0)
! ccne is all the factors with w in eq. A7 in Mansell et al. 2010 (JAS).  The constant changes
! if k (cck) is changed!
!      ccne0 = ccnefac*1.e6*(1.e-6)**(2./(2.+cck))
!         CN(mgs) =  Min( (frac+0.01)*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465

            sstmp = ccnefactwo*( (100.0*wvel(mgs))**(1.5)/(1.e-6*ccnc_ac(mgs)) )**(1.0/(cck+2.0))
            cn1 =  Min( ccnc_ac(mgs), CCNE0*ccnc_ac(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )

!           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
!           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

!           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
!              (T3D(K)*RR)-1.))
!     REAL, DIMENSION(KTS:KTE) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
!     REAL, DIMENSION(KTS:KTE) ::   EVS   ! SATURATION VAPOR PRESSURE
!     REAL, DIMENSION(KTS:KTE) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR

           ! following Abdul-Razzak et al. (1998, 'part 1') and Abdul-Razzak and Ghan (2000, 'part 2') for
           ! estimated maximum SS and one or more aerosol modes
           cpm = cp + cpv*qx(mgs,lv)
           evs = pres(mgs)*qvs(mgs)/0.622
           alpha_ar = gr*mwwater*felv(mgs)/(cpm*gasconst*temg(mgs)**2)-gr*Mair/(gasconst*temg(mgs))  ! AR98 eq. 11 (alpha)
           gamma_ar = gasconst*temg(mgs)/(evs*mwwater) + mwwater*felv(mgs)**2/(cpm*pres(mgs)*Mair*temg(mgs)) ! AR98 eq. 12 (gamma)

           G_ar = 1./(rhowater*gasconst*temg(mgs)/(evs*wvdf(mgs)*mwwater)+ felv(mgs)*rhowater/(tka(mgs)*temg(mgs))*(felv(mgs)*mwwater/ &
              (temg(mgs)*gasconst)-1.)) ! AR98 eq. 16

!            write(91,*) 'sstmp,ssf,cn1,wvel = ',sstmp,ssf(mgs),cn1,wvel(mgs)
!            sswater = Min(sstmp, ssf(mgs) * 0.01)  ! unit change from percentage to n/a
            sswater = ssf(mgs) * 0.01  ! unit change from percentage to n/a
            sigvl = 0.0761 - 1.55E-4*temcg(mgs)
            aact  = 2.0*mwwater/(rhowater*gasconst*temg(mgs))*sigvl  ! AR98 eq. 5 (A), where sigvl = tau = surface tension (wrongly defined as time in AR98)
            ! could add a subgrid TKE-based W-pert to wvel at some point
            zeta = 2./3.*(alpha_ar*wvel(mgs)/G_ar)**0.5*aact ! wvel is vertical velocity ! ARG2000 eq. 10

              f_ac = 0.5*Exp(2.5*(ac_pgw)**2)     ! ARG2000 eq. 7 (f_i), but size distribution in Cheng et al. uses sigma instead of log(sigma)
              g_ac = 1.0 + 0.25*ac_pgw           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
              eta_ac = (alpha_ar*wvel(mgs)/G_ar)**1.5/(2.*pi*rhowater*gamma_ar*ccnc_ac(mgs)) ! na zero ! ARG2000 eq. 11

              sm_ac = 2.0 / ac_kappa**0.5 * (aact / (3.0*ac_pmr) )**1.5    ! AR98 eq. 8 (S_mi)
!              IF ( wvel(mgs) < ac_wthresh ) THEN
                ss_ac = 1./sm_ac**2*(f_ac*(zeta/eta_ac)**1.5 + g_ac*(sm_ac**2/(eta_ac + 3.*zeta))**0.75)  ! ARG2000 eq. 6 (1 term in sum); comes from AR98 eq. 31
!              ELSE
!                ss_ac = 1.e20
!              ENDIF
            
            IF ( (ac_opt == 1 .or. ac_opt == 11) .and. ( sswater > 0.0 .or. ss_ac > 0.0 ) ) THEN  ! 
!              F11 = 0.5*EXP(2.5*(LOG(SIG1))**2) ! ARG2000 eq. 7 (f_i)
!              F21 = 1.+0.25*LOG(SIG1)           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
!              ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1) ! na zero ! ARG2000 eq. 11
!     sm_ac ->    SM1  = 2./BACT1**0.5*(AACT/(3.*RM1))**1.5  ! AR98 eq. 8 (or 7?) and ARG2000 eq. 9 for particle critical supersaturation
!              DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)  ! ARG2000 eq. 6 (1 term in sum); comes from AR98 eq. 31
!             SMAX = 1./(DUM1 + ... )**0.5
!      uu_ac ->      UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
!    NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
!             DUM1 = NANEW1/2.*(1.-DERF1(UU1))
!  BACT -> kappa
!  RM1 -> pmr

!-              sm_ac = 2.0 / ac_kappa**0.5 * (aact / (3.0*ac_pmr) )**1.5
              
              IF ( ss_ac > 0.0 ) THEN
               IF ( .true. ) THEN
                ! IF ( wvel(mgs) < ac_wthresh ) THEN
                   smax = Max(sswater, 1./sqrt(ss_ac) )
                ! ELSE
                !   smax = Min(sswater, 1./sqrt(ss_ac) )
                ! ENDIF
                 uu_ac = 2.0 * log(sm_ac / smax) / (3.0 * 2**0.5 * ac_pgw )  ! AR98 eq. 15 (u_i)
!                 tmp = 2.0 * log(sm_ac / sswater) / (3.0 * 2**0.5 * ac_pgw)
               ELSE
                 uu_ac = 2.0 * log(sm_ac / sswater) / (3.0 * 2**0.5 * ac_pgw)
               ENDIF
                cn(mgs) = Max( cn1, ccnc_ac(mgs) * 0.5 * (1.0 - Derf(uu_ac))  ) ! number concentration
                
               IF ( ac_opt == 11 ) THEN
                   cn_ac = Max( 0.0, cn(mgs) - ccna(mgs) )
                   cn(mgs) = cn_ac
                   !ccna(mgs) = ccna(mgs) + cn_ac
               ELSE
                  cn(mgs) = Max( 0.0, cn(mgs) - cx(mgs,lc)  ) ! number concentration of newly-activated droplets
               ENDIF
            !    write(0,*) 'i,k,smax,sswater = ',igs(mgs),kgs(mgs),smax,sswater,uu_ac, tmp, cn(mgs), Max( 0.0, ccnc_ac(mgs) * 0.5 * (1.0 - Derf(tmp))  ) 
              ENDIF
!              IF ( sswater > 0.0 ) &
!              write(91,*) 'cn,c_ac,cn-cx,uuac,sswater = ',cn(mgs),ccnc_ac(mgs),cn(mgs)- cx(mgs,lc),uu_ac,sswater
!              cn(mgs) = Min( cn(mgs)- cx(mgs,lc), ccnc_ac(mgs) - cx(mgs,lc) )
!              cn(mgs) = Min( cn(mgs), ccnc_ac(mgs) - cx(mgs,lc) )
            ELSEIF ( (ac_opt == 2 .or. ac_opt == 22 ) .and. ( sswater > 0.0 .or. ss_ac + ss_nu + ss_co > 0. ) ) THEN  ! multi-mode, use ccnc_xx(mgs) = an(igs(mgs),jy,kgs(mgs),lcn_xx) as number concentration

              f_nu = 0.5*Exp(2.5*(nu_pgw)**2)     ! ARG2000 eq. 7 (f_i), but size distribution in Cheng et al. uses sigma instead of log(sigma)
              g_nu = 1.0 + 0.25*nu_pgw           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
              eta_nu = (alpha_ar*wvel(mgs)/G_ar)**1.5/(2.*pi*rhowater*gamma_ar*ccnc_nu(mgs)) ! na zero ! ARG2000 eq. 11

              sm_nu = 2.0 / nu_kappa**0.5 * (aact / (3.0*nu_pmr) )**1.5    ! AR98 eq. 8 (S_mi)
!              IF ( wvel(mgs) < ac_wthresh ) THEN
                ss_nu = 1./sm_nu**2*(f_nu*(zeta/eta_nu)**1.5 + g_nu*(sm_nu**2/(eta_nu + 3.*zeta))**0.75)  ! ARG2000 eq. 6; comes from AR98 eq. 31
!              ELSE
!                ss_nu = 1.e20
!              ENDIF

              f_co = 0.5*Exp(2.5*(co_pgw)**2)     ! ARG2000 eq. 7 (f_i), but size distribution in Cheng et al. uses sigma instead of log(sigma)
              g_co = 1.0 + 0.25*co_pgw           ! ARG2000 eq. 8 (g_i); sig is sigma of the aerosol distribution
              eta_co = (alpha_ar*wvel(mgs)/G_ar)**1.5/(2.*pi*rhowater*gamma_ar*ccnc_co(mgs)) ! na zero ! ARG2000 eq. 11

              sm_co = 2.0 / co_kappa**0.5 * (aact / (3.0*co_pmr) )**1.5    ! AR98 eq. 8 (S_mi)
!              IF ( wvel(mgs) < ac_wthresh ) THEN
                ss_co = 1./sm_co**2*(f_co*(zeta/eta_co)**1.5 + g_co*(sm_co**2/(eta_co + 3.*zeta))**0.75)  ! ARG2000 eq. 6 (1 term in sum); comes from AR98 eq. 31
!              ELSE
!                ss_co = 1.e20
!              ENDIF

          !     IF ( wvel(mgs) < ac_wthresh ) THEN
                 smax = Max(sswater, 1./sqrt(ss_ac + ss_nu + ss_co) )
          !     ELSE
          !       smax = Min(sswater, 1./sqrt(ss_ac + ss_nu + ss_co) )
          !     ENDIF


             ! sm_nu = 2.0 / nu_kappa**0.5 * (aact / (3.0*nu_pmr) )**1.5
             ! sm_ac = 2.0 / ac_kappa**0.5 * (aact / (3.0*ac_pmr) )**1.5
             ! sm_co = 2.0 / co_kappa**0.5 * (aact / (3.0*co_pmr) )**1.5
              uu_nu = 2.0 * log(sm_nu / smax) / (3.0 * 2**0.5 * nu_pgw)
              uu_ac = 2.0 * log(sm_ac / smax) / (3.0 * 2**0.5 * ac_pgw)
              uu_co = 2.0 * log(sm_co / smax) / (3.0 * 2**0.5 * co_pgw)
!               uu_nu = 2.0 * log(sm_nu / sswater) / (3.0 * 2**0.5 * nu_pgw)
!               uu_ac = 2.0 * log(sm_ac / sswater) / (3.0 * 2**0.5 * ac_pgw)
!               uu_co = 2.0 * log(sm_co / sswater) / (3.0 * 2**0.5 * co_pgw)
              cn_co = Min( ccnc_co(mgs), ccnc_co(mgs) * 0.5 * (1.0 - Derf(uu_co)) ) 
              cn_ac = Min( ccnc_ac(mgs), ccnc_ac(mgs) * 0.5 * (1.0 - Derf(uu_ac)) ) 
              cn_nu = Min( ccnc_nu(mgs), ccnc_nu(mgs) * 0.5 * (1.0 - Derf(uu_nu)) ) 
              ! should try to initiate small rain drops from coarse mode?
              ! cn(mgs) = Min( cn(mgs)- cx(mgs,lc), ccnc_ac(mgs) + ccnc_nu(mgs) + ccnc_co(mgs) - cx(mgs,lc) )
              IF ( ac_opt == 2 ) THEN
              cn(mgs) = Min( ccnc_nu(mgs), ccnc_nu(mgs) * 0.5 * (1.0 - Derf(uu_nu)) ) + &
                      & Min( ccnc_ac(mgs), ccnc_ac(mgs) * 0.5 * (1.0 - Derf(uu_ac)) ) + &
                      & Min( ccnc_co(mgs), ccnc_co(mgs) * 0.5 * (1.0 - Derf(uu_co)) )     ! number concentration
            !  cn(mgs) = Min( cn(mgs), Max(0.0, ccnc_ac(mgs) + ccnc_nu(mgs) + ccnc_co(mgs) - cx(mgs,lc)) )
              cn(mgs) = Min( cn(mgs) - cx(mgs,lc), Max(0.0, ccnc_ac(mgs) + ccnc_nu(mgs) + ccnc_co(mgs) - cx(mgs,lc)) )
              ELSEIF ( ac_opt == 22 ) THEN
              !  write(0,*) 'cn_ac1,co,nu = ',igs(mgs),kgs(mgs),cn_ac,cn_co,cn_nu,ccnanu(mgs),ccnc_nu(mgs)
                cn_co = Max( 0.0, cn_co - ccnaco(mgs) ) ! should initiate rain drops from this??
                cn_ac = Max( 0.0, cn_ac - ccna(mgs) )
                cn_nu = Max( 0.0, cn_nu - ccnanu(mgs) )
                ! write(0,*) 'cn_ac2,co,nu = ',cn_ac,cn_co,cn_nu
                ! cn(mgs) = cx(mgs,lc) + cn_co + cn_ac + cn_nu ! add cx because of check method used below
                cn(mgs) =  cn_co + cn_ac + cn_nu ! add cx because of check method used below
              
              ENDIF
            ENDIF


       IF ( cn(mgs) > 0.0 ) THEN

       dcrit = 2.0*2.0e-6
       dcloud = 1000.*dcrit**3*Pi/6.
       ! check new droplet size:
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         ! tmp is number of droplets at diameter dcrit. If cn is too large and makes avg droplet diam too small, then replace with tmp
         IF ( tmp < cn(mgs) ) THEN
           IF ( lccna > 1 ) cn_ac = cn_ac*tmp/cn(mgs)
           IF ( lccnaco > 1 ) cn_co = cn_co*tmp/cn(mgs)
           IF ( lccnanu > 1 ) cn_nu = cn_nu*tmp/cn(mgs)
       !    write(0,*) 'cn_ac3,co,nu = ',igs(mgs),kgs(mgs),cn_ac,cn_co,cn_nu
          
           cn(mgs) = tmp ! Min(tmp, cn(mgs) )
         
         ENDIF

!          cx(mgs,lc) = Min(cn(mgs), cx(mgs,lc)+0.5*dqc/cwmasn)  ! limit the nucleation mass to half of the condensation mass
        IF ( ac_opt == 11 ) THEN
!          write(0,*) 'i,k,cn_ac2,co,nu = ',igs(mgs),kgs(mgs),cn_ac,cn_co,cn_nu,ccna(mgs),ccnc_ac(mgs),cn(mgs),cx(mgs,lc)
           ccna(mgs) = ccna(mgs) + cn_ac
        ELSEIF ( ac_opt == 22 ) THEN
           ccnanu(mgs) = ccnanu(mgs) + cn_nu
           ccnaco(mgs) = ccnaco(mgs) + cn_co
           ccna(mgs) = ccna(mgs) + cn_ac
        ENDIF

!!          cx(mgs,lc) = Max(cx(mgs,lc), cn(mgs) ) 
         cx(mgs,lc) = cx(mgs,lc) + cn(mgs)

 
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       
       dcrit = 2.0*2.5e-7
       dcloud = 1000.*dcrit**3*Pi/6.*Max(0.0,cn(mgs) + cnuf(mgs) )
        IF ( qx(mgs,lc) + dcloud <= qxmin(lc) ) dcloud = 1.001*( qxmin(lc) - qx(mgs,lc) )
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD

       ! IF ( ac_opt >= 11 ) THEN
       !   write(0,*) 'i,k,cn_ac3,cx,qx = ',igs(mgs),kgs(mgs),cx(mgs,lc),qx(mgs,lc)
       ! ENDIF

#ifdef WRFCODE
        IF ( numproc > 1 ) THEN
         dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
         thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
         thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        ENDIF
#endif
#ifdef COMMAS
        dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
        thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*DCLOUD/(pi0(mgs))*dv  ! latent heating
        thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + DCLOUD*rho0(mgs)*dv/dtp  ! condensation rate
        IF ( allocated( microrates ) ) THEN
          microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + DCLOUD/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
        ENDIF
        IF ( io_flag .and. lqcond >= 1 ) THEN
           axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + DCLOUD/dtp
        ENDIF
#endif
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
         ccncuf(mgs) = Max(0.0, ccncuf(mgs) - cnuf(mgs))
        ENDIF

            ! cn(mgs) = Max( cn(mgs), cnnumber_minimum )
!            IF ( iccwflg == 1 ) THEN  ! sets max size of first droplets in parcel to 4 micron radius (in two-moment liquid)
!              cn(mgs) = Max(cn(mgs), rho0(mgs)*qx(mgs,lc)/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3))
!            ENDIF
!          ELSE
!            cn(mgs) = 0.0
!            dcloud = 0.0
          ENDIF
      
      ENDIF
#endif
      ! IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .LE. qxmin(lc)) cx(mgs,lc)=0.
      GO TO 631
!.... NUCLEATION ON CLOUD INFLOW BOUNDARY POINT

  613 CONTINUE

  631  CONTINUE

!
! Check for supersaturation greater than ssmx and adjust down
!
       ssmx = maxsupersat
       qv1 = qv0(mgs) + qwvp(mgs)
       qvs1 = qvs(mgs)
       
!       IF ( flag_qndrop .and. do_satadj_for_wrfchem ) ssmx = 1.04 ! set lower threshold for progn=1 when using WRF-CHEM

       IF ( qv1 .gt. (ssmx*qvs1) ) THEN
! use line below to disable saturation adjustment when flag_qndrop is true
!       IF ( qv1 .gt. (ssmx*qvs1) .and. .not. flag_qndrop ) THEN
        
         ss1 = qv1/qvs1

        ssmx = 100.*(ssmx - 1.0)
        
        qvex = 0.0

        CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,qvex,   &
     &    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felvcp,ssmx,pk,ngscnt)



        IF ( qvex .gt. 0.0 ) THEN
          thetap(mgs) = thetap(mgs) + felvcp(mgs)*qvex/(pi0(mgs))
#ifdef WRFCODE
          IF ( numproc > 1 ) THEN
            dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
            thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*qvex/(pi0(mgs))*dv  ! latent heating
            thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + qvex*rho0(mgs)*dv/dtp  ! condensation rate
          ENDIF
#endif
#ifdef COMMAS
          dv = dxx(igs(mgs))*dyy(jgs)*dzz(kgs(mgs))
          thproc(kzbeg-1+kgs(mgs),16) = thproc(kzbeg-1+kgs(mgs),16) + felvcp(mgs)*qvex/(pi0(mgs))*dv  ! latent heating
          thproc(kzbeg-1+kgs(mgs),17) = thproc(kzbeg-1+kgs(mgs),17) + qvex*rho0(mgs)*dv/dtp  ! condensation rate
          IF ( allocated( microrates ) ) THEN
            microrates(igs(mgs),jy,kgs(mgs),3) = microrates(igs(mgs),jy,kgs(mgs),3) + qvex/dtp*felvcp(mgs)/(pi0(mgs))*3600.0
          ENDIF
          IF ( io_flag .and. lqcond >= 1 ) THEN
             axtra(igs(mgs),jy,kgs(mgs),lqcond) = axtra(igs(mgs),jy,kgs(mgs),lqcond) + qvex/dtp
          ENDIF
#else
          IF ( io_flag .and. nxtra > 1 ) THEN
             axtra(igs(mgs),jy,kgs(mgs),1) = axtra(igs(mgs),jy,kgs(mgs),1) + qvex/dtp
          ENDIF
#endif
          qwvp(mgs) = qwvp(mgs) - qvex
          qx(mgs,lc) = qx(mgs,lc) + qvex
          IF ( .not. flag_qndrop) THEN
            IF ( imaxsupopt == 1 ) THEN
              cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, xmas(mgs,lc) )  )
            ELSEIF ( imaxsupopt == 2 ) THEN
              cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmas30,xmas(mgs,lc)) )  )
            ELSEIF ( imaxsupopt == 3 ) THEN
              cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmasx,xmas(mgs,lc)) )  )
!            cn(mgs) = 1.5*cxmin
            ELSEIF ( imaxsupopt == 4 ) THEN
              cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmas20,xmas(mgs,lc)) )  )
            ENDIF

            IF ( lccna > 1 ) THEN
             !IF ( ac_opt == 0 ) THEN
               ccna(mgs) = ccna(mgs) + cn(mgs)
             !ENDIF
            ELSE
              ccnc(mgs) = Max( 0.0, ccnc(mgs) - cn(mgs) )
            ENDIF

            cx(mgs,lc) = cx(mgs,lc) + cn(mgs)

          ENDIF ! flag_qndrop

        ENDIF ! ( qvex .gt. 0.0 )

       ENDIF ! ( qv1 .gt. (ssmx*qvs1) )

!
! Calculate droplet volume and check if it is within bounds.
!  Adjust if necessary
!  
!      if (ndebug .gt. 0) write(0,*) "ICEZVD_DR: check droplet volume" 


!      cx(mgs,lc) = Min( cwnccn(mgs), cx(mgs,lc) )
      IF( cx(mgs,lc) > cxmin .AND. qx(mgs,lc) .GT. qxmin(lc)) THEN
!        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc)*xdn(mgs,lc))
        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc))
        
       IF (  xmas(mgs,lc) < cwmasn .or.  xmas(mgs,lc) > cwmasx ) THEN
        tmp = cx(mgs,lc)
        xmas(mgs,lc) = Min( xmas(mgs,lc), cwmasx )
        xmas(mgs,lc) = Max( xmas(mgs,lc), cwmasn )
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
!        IF ( cx(mgs,lc) > tmp*1.1 ) THEN
!        ENDIF
       ENDIF
      ENDIF


!      IF( cx(mgs,lc) .GT. 10.e6 .AND. qx(mgs,lc) .GT. qxmin(lc) ) GO TO 681
!        ccwtmp = cx(mgs,lc)
!        cwmastmp = xmas(mgs,lc)
!       xmas(mgs,lc) = Max(xmas(mgs,lc), cwmasn)
!       IF (qx(mgs,lc) .GT. qxmin(lc) .AND. cx(mgs,lc) .le. 0.) THEN
!          cx(mgs,lc) = Min(0.5*cwccn,rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc))
!          xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!       ENDIF
!      IF (cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. qxmin(lc))    &
!     &        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!      IF (qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .LT. cwmasn)    &
!     &          xmas(mgs,lc) = cwmasn
!      IF (qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .GT. cwmasx)    &
!     &    xmas(mgs,lc) = cwmasx
!      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN
!        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/Max(cwmasn,xmas(mgs,lc))
!      ENDIF
!        
!
! 681  CONTINUE

        
      IF ( ipconc .ge. 3 .and. rcond == 2 ) THEN

        
        IF (cx(mgs,lr) .GT. 0. .AND. qx(mgs,lr) .GT. qxmin(lr))    &
     &       xv(mgs,lr)=rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
        IF (xv(mgs,lr) .GT. xvmx(lr)) xv(mgs,lr) = xvmx(lr)
        IF (xv(mgs,lr) .LT. xvmn(lr)) xv(mgs,lr) = xvmn(lr)

      ENDIF



      ENDDO ! mgs


! ################################################################
      DO mgs=1,ngscnt
      IF ( lss > 1 .and. ssf(mgs) .gt. ssmax(mgs)    &
     &  .and. ( idecss .eq. 0 .or. qx(mgs,lc) .gt. qxmin(lc)) ) THEN
        ssmax(mgs) = ssf(mgs)
      ENDIF
      ENDDO
!

      do mgs = 1,ngscnt
      an(igs(mgs),jy,kgs(mgs),lt) = theta0(mgs) + thetap(mgs)
      an(igs(mgs),jy,kgs(mgs),lv) =  qv0(mgs) + qwvp(mgs)
!      tmp3d(igs(mgs),jy,kgs(mgs)) = tmp3d(igs(mgs),jy,kgs(mgs)) + t9(igs(mgs),jy,kgs(mgs)) !  pi0(mgs) ! wvdf(mgs) ! ssf(mgs) ! cn(mgs)
!
      IF ( eqtset > 2 ) THEN
        p2(igs(mgs),jy,kgs(mgs)) = pipert(mgs)
      ENDIF

       if ( ido(lc) .eq. 1 )  then
        an(igs(mgs),jy,kgs(mgs),lc) = qx(mgs,lc) +    &
     &    min( an(igs(mgs),jy,kgs(mgs),lc), 0.0 )
!        qx(mgs,lc) = an(igs(mgs),jy,kgs(mgs),lc)
       end if
!

       if ( ido(lr) .eq. 1 .and. rcond == 2 )  then
        an(igs(mgs),jy,kgs(mgs),lr) = qx(mgs,lr) +    &
     &    min( an(igs(mgs),jy,kgs(mgs),lr), 0.0 )
!        qx(mgs,lr) = an(igs(mgs),jy,kgs(mgs),lr)
       end if

#ifdef Z3MOM
        IF ( lzr > 1 .and. rcond == 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lzr) = zx(mgs,lr) +  &
     &    min( an(igs(mgs),jy,kgs(mgs),lzr), 0.0 )
        ENDIF
#endif


       IF (  ipconc .ge. 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnc) = Max(cx(mgs,lc) , 0.0)
       ! IF ( ac_opt > 10 .and. (cx(mgs,lc) > 0. .or. ccna(mgs) > 0. ) ) THEN
       !   write(0,*) 'i,k final cx/cna = ',igs(mgs),kgs(mgs),cx(mgs,lc),ccna(mgs)
       ! ENDIF
        
        IF ( lss > 1 ) an(igs(mgs),jy,kgs(mgs),lss) = Max( 0.0, ssmax(mgs) )
        IF ( ac_opt == 0 ) THEN
          IF ( lccn .gt. 1 .and. lccna .lt. 1  ) THEN
            an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0,  ccnc(mgs) )
          ENDIF
#ifdef NUWRFMODS
        ELSEIF ( ac_opt == 1 .and. lcn_ac > 1) THEN
            an(igs(mgs),jy,kgs(mgs),lcn_ac) = Max( 0.0, ccnc_ac(mgs) )
        ELSEIF ( ac_opt == 11 .and. lccna > 1) THEN
            an(igs(mgs),jy,kgs(mgs),lccna) = Max( 0.0, ccna(mgs) )
        ELSEIF ( ac_opt == 2 .and. lcn_ac > 1) THEN
            an(igs(mgs),jy,kgs(mgs),lcn_ac) = Max( 0.0, ccnc_ac(mgs) )
            an(igs(mgs),jy,kgs(mgs),lcn_nu) = Max( 0.0, ccnc_nu(mgs) )
            an(igs(mgs),jy,kgs(mgs),lcn_co) = Max( 0.0, ccnc_co(mgs) )
        ELSEIF ( ac_opt == 22 .and. lccna > 1) THEN
            an(igs(mgs),jy,kgs(mgs),lccna) = Max( 0.0, ccna(mgs) )
            an(igs(mgs),jy,kgs(mgs),lccnanu) = Max( 0.0, ccnanu(mgs) )
            an(igs(mgs),jy,kgs(mgs),lccnaco) = Max( 0.0, ccnaco(mgs) )
#endif
        ENDIF
        IF ( lccnuf .gt. 1 .and. .not. ( lccna .gt. 1 .and. i_uf_or_ccn > 0 ) ) THEN
          an(igs(mgs),jy,kgs(mgs),lccnuf) = Max(0.0,  ccncuf(mgs) )
        ENDIF
        IF ( lccna .gt. 1 ) THEN
          an(igs(mgs),jy,kgs(mgs),lccna) = Max(0.0, ccna(mgs) )
        ENDIF
       ENDIF
       IF (  ipconc .ge. 3 .and. rcond == 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnr) = Max(cx(mgs,lr) , 0.0)
       ENDIF
      end do


29998 continue


      if ( kz .gt. nz-1 .and. ix .ge. nxi) then
        if ( ix .ge. nxi ) then
         go to 2200 ! exit gather scatter
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz
      end if

      if ( ix .ge. nxi ) then
        nxmpb = 1
        nzmpb = kz+1
      else
       nxmpb = ix+1
      end if

 2000 continue ! inumgs
 2200 continue
!
!  end of gather scatter (for this jy slice)


! Redistribute inappreciable cloud particles and charge
!
! Redistribution everywhere in the domain...
!
! moved to separate subroutine (below)
!

   
   
   9999 RETURN
   
   END SUBROUTINE NUCOND


! #####################################################################
! #####################################################################
! Clean up tiny values of mixing ratio
! Redistribute inappreciable cloud particles and charge
!
! Redistribution everywhere in the domain...
!
     subroutine smallvalues   &
     &  (nx,ny,nz,na,jyslab & 
     &  ,nor,norz,dtp,nxi & 
     &  ,t0 & 
     &  ,an,dn, w & 
#ifdef COMMAS
     &  ,qxmin, xdn0, xvmn, xvmx, cno  &
     &  ,ido,lsc,ln,ipc,lvol,lz,lliq,lrain &
     &  ,xdnmn,xdnmx  &
#endif
     &  ,t77,flag_qndrop  &
     & )

#ifdef COMMAS
       USE INDEX_MODULE !, cwmasn_index => cwmasn
       USE MICRO_MODULE
#endif

      implicit none

      integer :: nx,ny,nz,na,nxi
      integer :: nor,norz, jyslab ! ,nht,ngt,igsr
      real    :: dtp  ! time step
      logical,intent(in) :: flag_qndrop

!
! external temporary arrays
!
      real t77(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

    ! local

#ifdef COMMAS
      ! passed in
      real qxmin(lc:lqmx)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real cno(lc:lhab)
      integer ido(lc:lqmx)
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      integer lvol(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
#if USERAINTYPES > 0
      integer lrain(nraintypes)
#endif
      real xdnmx(lc:lhab), xdnmn(lc:lhab)

      real tfr,tfrh
      parameter ( tfr = 273.15, tfrh = 233.15)

      real, parameter :: ec = 1.602e-19 ! fundamental unit of charge
      real, parameter :: eci = 1.0/ec

      real, parameter :: pi = 3.141592653589793
      real, parameter :: piinv = 1./pi

#endif
 
      logical zerocx(lc:lqmx)

      real :: frac, hwdn, tmpg, xdia1, xdia3, cwch,xvol

      integer ix,kz,i,n, km1
      integer :: il
      integer :: jy, jgs
      real :: chw, g1, z1, tmp, fw, tmpmx, qr
    

! Redistribute inappreciable cloud particles and charge
!
! Redistribution everywhere in the domain...
!
#ifdef COMMAS
       jy = jyslab
#else
       jy = 1
#endif

      frac = 1.0 ! 0.25 ! 1.0 ! 0.2

      cwch = ((3. + alphah)*(2. + alphah)*(1.0 + alphah))**(-1./3.)
!
!  alternate test version for ipconc .ge. 3
!  just vaporize stuff to prevent noise in the number concentrations


      do kz = 1,nz
!      do jy = 1,1
      do ix = 1,nxi
      
      t0(ix,jy,kz) = an(ix,jy,kz,lt)*t77(ix,jy,kz)
      
      zerocx(:) = .false.
      DO il = lc,lhab
       IF ( iresetmoments == 1 .or. iresetmoments == il ) THEN
        IF ( ln(il) > 1 ) zerocx(il) = ( an(ix,jy,kz,ln(il)) < cxmin )
        IF ( lz(il) > 1 ) zerocx(il) = ( zerocx(il) .or. (an(ix,jy,kz,lz(il)) < zxmin) )
       ELSE
        IF ( il == lc ) THEN
          IF ( ln(il) > 1 ) THEN
           zerocx(il) = ( an(ix,jy,kz,ln(il)) < cxmin ) .and. .not. flag_qndrop ! do not reset if progn=1 (WRF-CHEM)
          ENDIF
        ELSE
         IF ( ln(il) > 1 ) zerocx(il) = ( an(ix,jy,kz,ln(il)) < cxmin )
        ENDIF
       ENDIF
      ENDDO

      IF ( lhl .gt. 1 ) THEN
      
#ifdef Z3MOM
      IF ( lzhl .gt. 1 ) THEN

        an(ix,jy,kz,lzhl) = Max(0.0, an(ix,jy,kz,lzhl) )
        
        IF ( an(ix,jy,kz,lhl) .ge. frac*qxmin(lhl) .and. rescale_low_alpha ) THEN ! check 6th moment
          
          IF ( an(ix,jy,kz,lnhl) .gt. 0.0 ) THEN

           IF ( lvhl .gt. 1 ) THEN
             IF ( an(ix,jy,kz,lvhl) .gt. 0.0 ) THEN
               hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
             ELSE
               hwdn = xdn0(lhl)
             ENDIF
             hwdn = Max( xdnmn(lhl), hwdn )
           ELSE
             hwdn = xdn0(lhl)
           ENDIF

             chw = an(ix,jy,kz,lnhl)
             g1 = (6.0+alphamin)*(5.0+alphamin)*(4.0+alphamin)/   &
     &            ((3.0+alphamin)*(2.0+alphamin)*(1.0+alphamin))
             z1 = g1*dn(ix,jy,kz)**2*( an(ix,jy,kz,lhl) )*an(ix,jy,kz,lhl)/chw
             z1 = z1*(6./(pi*hwdn))**2
          ELSE
             z1 = 0.0
          ENDIF
          
          an(ix,jy,kz,lzhl) = Min( z1, an(ix,jy,kz,lzhl) )
          
          IF (  an(ix,jy,kz,lnhl) .lt. 1.e-5 ) THEN
!            an(ix,jy,kz,lzhl) = 0.9*an(ix,jy,kz,lzhl)
          ENDIF
        ENDIF
        
      ENDIF !lzhl
#endif
      
      if ( (an(ix,jy,kz,lhl) .lt. frac*qxmin(lhl)) .or. zerocx(lhl) ) then

!        IF ( an(ix,jy,kz,lhl) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lhl)
          an(ix,jy,kz,lhl) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnhl) = 0.0
        ENDIF

        IF ( lvhl .gt. 1 ) THEN
           an(ix,jy,kz,lvhl) = 0.0
        ENDIF

        IF ( lhlw .gt. 1 ) THEN
           an(ix,jy,kz,lhlw) = 0.0
        ENDIF

        IF ( lnhlf .gt. 1 ) THEN
           an(ix,jy,kz,lnhlf) = 0.0
        ENDIF
      
        IF ( lzhl .gt. 1 ) THEN
           an(ix,jy,kz,lzhl) = 0.0
        ENDIF
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lschl .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lschl) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lschl)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lschl)*eci
          ENDIF
          
           an(ix,jy,kz,lschl) = 0.0
        ENDIF
#endif

      ELSE
       IF ( lvol(lhl) .gt. 1 ) THEN  ! check density
        IF ( an(ix,jy,kz,lvhl) .gt. 0.0 ) THEN
         tmp = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
        ELSE ! in case volume is zero but mass is above threshold (should not happen, of course)
          tmp = rho_qhl
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        ENDIF

        IF (  tmp .lt. xdnmn(lhl) ) THEN
          tmp = Max( xdnmn(lhl), tmp )
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        ENDIF

        IF ( tmp .gt. xdnmx(lhl) .and. lhlw .le. 0 ) THEN ! no liquid allowed on hail
          tmp = Min( xdnmx(lhl), tmp )
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        ELSEIF ( tmp .gt. xdnmx(lhl) .and. lhlw .gt. 1 ) THEN  ! allow for liquid on hail
          fw = an(ix,jy,kz,lhlw)/an(ix,jy,kz,lhl)
!          tmpmx = xdnmx(lhl) + fw*(xdnmx(lr) - xdnmx(lhl)) ! maximum possible average density
                                                           ! it is not exactly linear, but approx. is close enough for this
!          tmpmx = 1./( (1. - fw)/900. + fw/1000. ) is exact max, where 900 is xdnmx

          tmpmx = xdnmx(lhl)/( 1. - fw*(1. - xdnmx(lhl)/xdnmx(lr) )) 

          IF ( tmp .gt. tmpmx  ) THEN
            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmpmx
          ENDIF

!          IF ( tmp .gt. xdnmx(lhl) .and. an(ix,jy,kz,lhlw) .lt. qxmin(lhl) ) THEN
!            tmp = Min( xdnmx(lhl), tmp )
!            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
!          ELSEIF ( tmp .gt. xdnmx(lr) ) THEN
!            tmp =  xdnmx(lr)
!            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
!          ENDIF
        ENDIF

        IF ( lhlw .gt. 1 ) THEN ! check if basically pure water
          IF ( an(ix,jy,kz,lhlw) .gt. 0.98*an(ix,jy,kz,lhl) ) THEN
           tmp = xdnmx(lr)
           an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
          ENDIF
        ENDIF
        
       ENDIF
       
         IF ( lvhl .gt. 1 ) THEN
           IF ( an(ix,jy,kz,lvhl) .gt. 0.0 ) THEN
             hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
           ELSE
             hwdn = xdn0(lhl)
           ENDIF
           hwdn = Max( xdnmn(lhl), hwdn )
         ELSE
           hwdn = xdn0(lhl)
         ENDIF

         IF ( ipconc >= 5 .and.  an(ix,jy,kz,lhl) .gt. qxmin(lhl) ) THEN
            qr = an(ix,jy,kz,lhl)
            xvol = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/(hwdn*an(ix,jy,kz,lnhl))
            chw = an(ix,jy,kz,lnhl)

             IF ( xvol .lt. xvmn(lhl) .or. xvol .gt. xvmx(lhl) ) THEN
              xvol = Min( xvmx(lhl), Max( xvmn(lhl),xvol ) )
              chw = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/(xvol*hwdn)
              an(ix,jy,kz,lnhl) = chw
             ENDIF
          ENDIF
       
!  CHECK INTERCEPT
       IF ( ipconc == 5 .and.  an(ix,jy,kz,lhl) .gt. qxmin(lhl) .and.  alphahl .le. 0.1 .and. lnhl .gt. 1 .and. lzhl == 0 ) THEN
       
         IF ( lvhl .gt. 1 ) THEN
           hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
         ELSE
           hwdn = xdn0(lhl)
         ENDIF
           tmp = (hwdn*an(ix,jy,kz,lnhl))/(dn(ix,jy,kz)*an(ix,jy,kz,lhl))
           tmpg = an(ix,jy,kz,lnhl)*(tmp*pi)**(1./3.)
           IF ( tmpg .lt. cnohlmn ) THEN
             tmp = ( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lhl))*pi)**(1./3.)
              an(ix,jy,kz,lnhl) = (cnohlmn/tmp)**(3./4.)
           ENDIF
       
       ENDIF
!      ELSE  ! check mean size here?

      end if

      ENDIF !lhl

#if USEFROZENDROPS > 0
        IF ( lf .gt. 1 ) THEN
      
#ifdef Z3MOM
      IF ( lzf .gt. 1 ) THEN

        an(ix,jy,kz,lzf) = Max(0.0, an(ix,jy,kz,lzf) )
        
        IF ( .false. .and. an(ix,jy,kz,lf) .ge. frac*qxmin(lf) .and. rescale_low_alpha  ) THEN
          
          IF ( an(ix,jy,kz,lnf) .gt. 0.0 ) THEN

           IF ( lvf .gt. 1 ) THEN
             IF ( an(ix,jy,kz,lvf) .gt. 0.0 ) THEN
               hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
             ELSE
               hwdn = xdn0(lf)
             ENDIF
             hwdn = Max( xdnmn(lf), hwdn )
           ELSE
             hwdn = xdn0(lf)
           ENDIF

             chw = an(ix,jy,kz,lnf)
             g1 = (6.0+alphamin)*(5.0+alphamin)*(4.0+alphamin)/   &
     &            ((3.0+alphamin)*(2.0+alphamin)*(1.0+alphamin))
             z1 = g1*dn(ix,jy,kz)**2*( an(ix,jy,kz,lf) )*an(ix,jy,kz,lf)/chw
             z1  = z1*(6./(pi*hwdn))**2
          ELSE
             z1 = 0.0
          ENDIF
          
          an(ix,jy,kz,lzf) = Min( z1, an(ix,jy,kz,lzf) )
          
          IF (  an(ix,jy,kz,lnf) .lt. 1.e-5 ) THEN
!            an(ix,jy,kz,lzf) = 0.9*an(ix,jy,kz,lzf)
          ENDIF
        ENDIF
        
      ENDIF !lzf

#endif
      
      if ( (an(ix,jy,kz,lf) .lt. frac*qxmin(lf)) .or. zerocx(lf) ) then

!        IF ( an(ix,jy,kz,lf) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lf)
          an(ix,jy,kz,lf) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnf) = 0.0
        ENDIF

        IF ( lvf .gt. 1 ) THEN
           an(ix,jy,kz,lvf) = 0.0
        ENDIF

        IF ( lfw .gt. 1 ) THEN
           an(ix,jy,kz,lfw) = 0.0
        ENDIF
      
        IF ( lzf .gt. 1 ) THEN
           an(ix,jy,kz,lzf) = 0.0
        ENDIF
      
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscf .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lscf) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscf)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscf)*eci
          ENDIF
           an(ix,jy,kz,lscf) = 0.0
        ENDIF
#endif

      ELSE
       IF ( lvol(lf) .gt. 1 ) THEN  ! check density
        IF ( an(ix,jy,kz,lvf) .gt. 0.0 ) THEN
         tmp = dn(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
        ELSE 
         tmp = 0.5*( xdnmn(lf) + xdnmx(lf) )
          an(ix,jy,kz,lvf) = dn(ix,jy,kz)*an(ix,jy,kz,lf)/tmp
        ENDIF
        
        IF ( tmp .lt. xdnmn(lf) ) THEN
          tmp =  Max( xdnmn(lf) , tmp )
          an(ix,jy,kz,lvf) = dn(ix,jy,kz)*an(ix,jy,kz,lf)/tmp
        ENDIF
        
        IF ( tmp .gt. xdnmx(lf) .and. lfw .le. 0 ) THEN ! no liquid allowed on hail
          tmp = Min( xdnmx(lf), tmp )
          an(ix,jy,kz,lvf) = dn(ix,jy,kz)*an(ix,jy,kz,lf)/tmp
        ELSEIF ( tmp .gt. xdnmx(lf) .and. lfw .gt. 1 ) THEN  ! allow for liquid on hail
          fw = an(ix,jy,kz,lfw)/an(ix,jy,kz,lf)
          tmpmx = xdnmx(lf) + fw*(xdnmx(lr) - xdnmx(lf)) ! maximum possible average density
                                                           ! it is not exactly linear, but approx. is close enough for this
!          tmpmx = 1./( (1. - fw)/900. + fw/1000. ) is exact max, where 900 is xdnmx
         
          IF ( tmp .gt. tmpmx  ) THEN
!            tmp = Min( xdnmx(lf), tmp )
            an(ix,jy,kz,lvf) = dn(ix,jy,kz)*an(ix,jy,kz,lf)/tmpmx
          ENDIF
        ENDIF
        
        IF ( lfw .gt. 1 ) THEN ! check if basically pure water
          IF ( an(ix,jy,kz,lfw) .gt. 0.98*an(ix,jy,kz,lf) ) THEN
           tmp = xdnmx(lr)
           an(ix,jy,kz,lvf) = dn(ix,jy,kz)*an(ix,jy,kz,lf)/tmp
          ENDIF
        ENDIF

        
       ENDIF

         IF ( lvf .gt. 1 ) THEN
           IF ( an(ix,jy,kz,lvf) .gt. 0.0 ) THEN
             hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
           ELSE
             hwdn = xdn0(lf)
           ENDIF
           hwdn = Max( xdnmn(lf), hwdn )
         ELSE
           hwdn = xdn0(lf)
         ENDIF

         IF ( ipconc >= 5 .and.  an(ix,jy,kz,lf) .gt. qxmin(lh) ) THEN
            qr = an(ix,jy,kz,lf)
            xvol = dn(ix,jy,kz)*an(ix,jy,kz,lf)/(hwdn*an(ix,jy,kz,lnf))
            chw = an(ix,jy,kz,lnf)

             IF ( xvol .lt. xvmn(lf) .or. xvol .gt. xvmx(lf) ) THEN
              xvol = Min( xvmx(lf), Max( xvmn(lf),xvol ) )
              chw = dn(ix,jy,kz)*an(ix,jy,kz,lf)/(xvol*hwdn)
              an(ix,jy,kz,lnf) = chw
             ENDIF
          ENDIF

       
       
!  CHECK INTERCEPT
#if 0
       IF ( ipconc == 5 .and. an(ix,jy,kz,lf) .gt. qxmin(lf) .and.  alphah .le. 0.5 .and. lnf .gt. 1 ) THEN
       
         IF ( lvf .gt. 1 ) THEN
           hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
         ELSE
           hwdn = xdn0(lf)
         ENDIF
           tmp = (hwdn*an(ix,jy,kz,lnf))/(dn(ix,jy,kz)*an(ix,jy,kz,lf))
           tmpg = an(ix,jy,kz,lnf)*(tmp*pi)**(1./3.)
           IF ( tmpg .lt. cnohmn ) THEN
!           tmpg = an(ix,jy,kz,lnh)*( (hwdn*an(ix,jy,kz,lnh))/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
!           tmpg = an(ix,jy,kz,lnh)**(4./3.)*( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
             tmp = ( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lf))*pi)**(1./3.)
              an(ix,jy,kz,lnf) = (cnohmn/tmp)**(3./4.)
           ENDIF
       
       ENDIF
#endif
!      ELSE  ! check mean size here?
        
      end if
      
      ENDIF !lf

#endif /* frozen drops */


#ifdef Z3MOM
      IF ( lzh .gt. 1 ) THEN

        an(ix,jy,kz,lzh) = Max(0.0, an(ix,jy,kz,lzh) )
        
        IF ( .false. .and. an(ix,jy,kz,lh) .ge. frac*qxmin(lh) .and. rescale_low_alpha ) THEN
          
          IF ( an(ix,jy,kz,lnh) .gt. 0.0 ) THEN

           IF ( lvh .gt. 1 ) THEN
             IF ( an(ix,jy,kz,lvh) .gt. 0.0 ) THEN
               hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
             ELSE
               hwdn = xdn0(lh)
             ENDIF
             hwdn = Max( xdnmn(lh), hwdn )
           ELSE
             hwdn = xdn0(lh)
           ENDIF

             chw = an(ix,jy,kz,lnh)
             g1 = (6.0+alphamin)*(5.0+alphamin)*(4.0+alphamin)/   &
     &            ((3.0+alphamin)*(2.0+alphamin)*(1.0+alphamin))
             z1 = g1*dn(ix,jy,kz)**2*( an(ix,jy,kz,lh) )*an(ix,jy,kz,lh)/chw
             z1  = z1*(6./(pi*hwdn))**2
          ELSE
             z1 = 0.0
          ENDIF
          
          an(ix,jy,kz,lzh) = Min( z1, an(ix,jy,kz,lzh) )
          
          IF (  an(ix,jy,kz,lnh) .lt. 1.e-5 ) THEN
!            an(ix,jy,kz,lzh) = 0.9*an(ix,jy,kz,lzh)
          ENDIF
        ENDIF
        
      ENDIF
#endif

      if ( (an(ix,jy,kz,lh) .lt. frac*qxmin(lh)) .or. zerocx(lh) ) then

!        IF ( an(ix,jy,kz,lh) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lh)
          an(ix,jy,kz,lh) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnh) = 0.0
        ENDIF

        IF ( lvh .gt. 1 ) THEN
           an(ix,jy,kz,lvh) = 0.0
        ENDIF
      
        IF ( lhw .gt. 1 ) THEN
           an(ix,jy,kz,lhw) = 0.0
        ENDIF

        IF ( lnhf .gt. 1 ) THEN
           an(ix,jy,kz,lnhf) = 0.0
        ENDIF
      
        IF ( lzh .gt. 1 ) THEN
           an(ix,jy,kz,lzh) = 0.0
        ENDIF
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lsch .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lsch) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lsch)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lsch)*eci
          ENDIF
           an(ix,jy,kz,lsch) = 0.0
        ENDIF
#endif

      ELSE
       IF ( lvol(lh) .gt. 1 ) THEN  ! check density
        IF ( an(ix,jy,kz,lvh) .gt. 0.0 ) THEN
         tmp = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
        ELSE
         tmp = rho_qh
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        ENDIF

        IF (  tmp .lt. xdnmn(lh) ) THEN
          tmp = Max( xdnmn(lh), tmp )
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        ENDIF

        IF ( tmp .gt. xdnmx(lh) .and. lhw .le. 0 ) THEN ! no liquid allowed on graupel
          tmp = Min( xdnmx(lh), tmp )
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        ELSEIF ( tmp .gt. xdnmx(lh) .and. lhw .gt. 1 ) THEN  ! allow for liquid on graupel
          fw = an(ix,jy,kz,lhw)/an(ix,jy,kz,lh)
!          tmpmx = xdnmx(lh) + fw*(xdnmx(lr) - xdnmx(lh)) ! maximum possible average density
                                                           ! it is not exactly linear, but approx. is close enough for this
!          tmpmx = 1./( (1. - fw)/900. + fw/1000. ) is exact max, where 900 is xdnmx
          tmpmx = xdnmx(lh)/( 1. - fw*(1. - xdnmx(lh)/xdnmx(lr) )) 

          IF ( tmp .gt. tmpmx  ) THEN
            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmpmx
          ENDIF

!          IF ( tmp .gt. xdnmx(lh) .and. an(ix,jy,kz,lhw) .lt. qxmin(lh) ) THEN
!            tmp = Min( xdnmx(lh), tmp )
!            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
!          ELSEIF ( tmp .gt. xdnmx(lr) ) THEN
!            tmp =  xdnmx(lr)
!            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
!          ENDIF

        ENDIF

        IF ( lhw .gt. 1 ) THEN ! check if basically pure water
          IF ( an(ix,jy,kz,lhw) .gt. 0.98*an(ix,jy,kz,lh) ) THEN
           tmp = xdnmx(lr)
           an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
          ENDIF
        ENDIF
        
       ENDIF

         IF ( lvh .gt. 1 ) THEN
           IF ( an(ix,jy,kz,lvh) .gt. 0.0 ) THEN
             hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
           ELSE
             hwdn = xdn0(lh)
           ENDIF
           hwdn = Max( xdnmn(lh), hwdn )
         ELSE
           hwdn = xdn0(lh)
         ENDIF

         IF ( ipconc >= 5 .and.  an(ix,jy,kz,lh) .gt. qxmin(lh) ) THEN
            qr = an(ix,jy,kz,lh)
            xvol = dn(ix,jy,kz)*an(ix,jy,kz,lh)/(hwdn*an(ix,jy,kz,lnh))
            chw = an(ix,jy,kz,lnh)

             IF ( xvol .lt. xvmn(lh) .or. xvol .gt. xvmx(lh) ) THEN
              xvol = Min( xvmx(lh), Max( xvmn(lh),xvol ) )
              chw = dn(ix,jy,kz)*an(ix,jy,kz,lh)/(xvol*hwdn)
              an(ix,jy,kz,lnh) = chw
             ENDIF
          ENDIF

!  CHECK INTERCEPT
       IF ( ipconc == 5 .and.  an(ix,jy,kz,lh) .gt. qxmin(lh) .and.  alphah .le. 0.1 .and. lnh .gt. 1 .and. lzh == 0 ) THEN
       
           tmp = (hwdn*an(ix,jy,kz,lnh))/(dn(ix,jy,kz)*an(ix,jy,kz,lh))
           tmpg = an(ix,jy,kz,lnh)*(tmp*pi)**(1./3.)
           IF ( tmpg .lt. cnohmn ) THEN
!           tmpg = an(ix,jy,kz,lnh)*( (hwdn*an(ix,jy,kz,lnh))/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
!           tmpg = an(ix,jy,kz,lnh)**(4./3.)*( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
             tmp = ( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*pi)**(1./3.)
              an(ix,jy,kz,lnh) = (cnohmn/tmp)**(3./4.)
           ENDIF
       
       ENDIF

         IF (  ipconc == 5 .and. imorrgdnglimit == 1 ) THEN
           ! limit on characteristic diameter (i.e., 1/slope)
            xdia3 = (xvol*6.*piinv)**(1./3.)
            xdia1 = cwch*xdia3
            IF ( xdia1 > morrdnglimit ) THEN
               xdia1 = morrdnglimit
               xvol = pi/6.0*(xdia1/cwch)**3
               chw = dn(ix,jy,kz)*qr/(xvol*hwdn)
               an(ix,jy,kz,lnh) = chw
               xdia3 = (xvol*6.*piinv)**(1./3.)
            ENDIF
         
         ENDIF
        
      end if


      if ( (an(ix,jy,kz,ls) .lt.  frac*qxmin(ls))  .or. zerocx(ls) ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF
      
        IF ( ipconc .ge. 4 ) THEN ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
!          an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF
        
        IF ( lvs .gt. 1 ) THEN
           an(ix,jy,kz,lvs) = 0.0
        ENDIF

        IF ( lsw .gt. 1 ) THEN
           an(ix,jy,kz,lsw) = 0.0
        ENDIF

      ELSE
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF

        IF ( lvs .gt. 1 ) THEN
           an(ix,jy,kz,lvs) = 0.0
        ENDIF

        IF ( lsw .gt. 1 ) THEN
           an(ix,jy,kz,lsw) = 0.0
        ENDIF

        IF ( ipconc .ge. 4 ) THEN ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF

      ENDIF
      
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscs > 1 ) THEN
          IF ( an(ix,jy,kz,lscs) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscs)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscs)*eci
          ENDIF
          
          an(ix,jy,kz,lscs) = 0.0
          
        ENDIF
#endif

      ELSEIF ( lvol(ls) .gt. 1 ) THEN  ! check density
        IF ( an(ix,jy,kz,lvs) .gt. 0.0 ) THEN
          tmp = dn(ix,jy,kz)*an(ix,jy,kz,ls)/an(ix,jy,kz,lvs)
          IF ( tmp .gt. xdnmx(ls) .or. tmp .lt. xdnmn(ls) ) THEN
            tmp = Min( xdnmx(ls), Max( xdnmn(ls), tmp ) )
            an(ix,jy,kz,lvs) = dn(ix,jy,kz)*an(ix,jy,kz,ls)/tmp
          ENDIF
        ELSE
          tmp = rho_qs
          an(ix,jy,kz,lvs) = dn(ix,jy,kz)*an(ix,jy,kz,ls)/tmp
        ENDIF


      end if

#ifdef Z3MOM
        IF ( lzr > 1 ) THEN
          an(ix,jy,kz,lzr) = Max(0.0, an(ix,jy,kz,lzr) )
        ENDIF
#endif

      if ( (an(ix,jy,kz,lr) .lt. frac*qxmin(lr))  .or. zerocx(lr) ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        IF ( ipconc .ge. 3 ) THEN
!          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
        ENDIF
        
#ifdef Z3MOM
        IF ( lzr > 1 ) THEN
          an(ix,jy,kz,lzr) = 0.0
        ENDIF
#endif
#if USERAINTYPES > 0
        IF ( iraintypes >= 1 ) THEN
          DO il = 1,nraintypes
            an(ix,jy,kz,lrain(il)) = 0.0
          ENDDO
        ENDIF
#endif

#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscr .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lscr) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscr)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscr)*eci
          ENDIF

           an(ix,jy,kz,lscr) = 0.0

        ENDIF
#endif
      end if

!
!  for qi
!
      IF ( (an(ix,jy,kz,li) .le. frac*qxmin(li)) .or. zerocx(li) ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
      an(ix,jy,kz,li)= 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lni) = 0.0
       ENDIF
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lsci .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lsci) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lsci)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lsci)*eci
          ENDIF
           an(ix,jy,kz,lsci) = 0.0
        ENDIF
#endif
      ENDIF

#ifdef USEICESPHERES
!
!  for qis
!
      IF ( lis > 1 ) THEN ! {
      IF ( (an(ix,jy,kz,lis) .le. frac*qxmin(lis)) .or. zerocx(lis)   & ! .or.  an(ix,jy,kz,lni) .lt. 0.1
     &    ) THEN ! { {
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lis)
      an(ix,jy,kz,lis)= 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lnis) = 0.0
       ENDIF
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscis .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lscis) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscis)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscis)*eci
          ENDIF
           an(ix,jy,kz,lscis) = 0.0
        ENDIF
#endif
      
      ELSEIF ( icespheres >= 2 ) THEN ! } {
       km1 = Max(1, kz-1)
       IF ( 0.5*( w(ix,jy,kz) + w(ix,jy,kz+1)) < -1.0 .or.    &
     &      (icespheres == 3 .and. ( t0(ix,jy,kz) < 232.15 .or. an(ix,jy,kz,lc) < qxmin(lc) ) ) .or. &
     &      (icespheres == 5 .and. ( t0(ix,jy,kz) < 232.15 .or. &
     &         ( an(ix,jy,kz,lc) < qxmin(lc) .and. an(ix,jy,km1,lc) < qxmin(lc)  )) ) .or. &
     &      (icespheres == 4 .and. ( t0(ix,jy,kz) < 235.15 )) ) THEN ! transfer to regular ice crystals in downdraft or at low temp
         an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,lis)
         an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lnis)
         an(ix,jy,kz,lis)= 0.0
         an(ix,jy,kz,lnis)= 0.0
#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscis .gt. 1 ) THEN
         an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscis)
         an(ix,jy,kz,lscis)= 0.0
        ENDIF
#endif
         
       ENDIF
       
      ENDIF ! } }
      ENDIF ! }
#endif /* USEICESPHERES */
!
!  for qcw
!

      IF ( (an(ix,jy,kz,lc) .le. frac*qxmin(lc)) .or. zerocx(lc) ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0
       IF ( ipconc .ge. 2 ) THEN
        IF ( lccn .gt. 1 .or. ac_opt == 1 ) THEN
          IF ( irenuc < 5 .and. lccna <= 1 ) THEN
            IF ( ac_opt == 0 ) THEN
               an(ix,jy,kz,lccn) = an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc))
#ifdef NUWRFMODS
            ELSEIF ( lcn_ac > 1 ) THEN
               an(ix,jy,kz,lcn_ac) = an(ix,jy,kz,lcn_ac) + Max(0.0,an(ix,jy,kz,lnc))
#endif
            ENDIF
          ELSEIF ( lccna > 1 ) THEN
            an(ix,jy,kz,lccna) = Max( 0.0, an(ix,jy,kz,lccna) - Max(0.0,an(ix,jy,kz,lnc)) )
          ENDIF
        ENDIF
         an(ix,jy,kz,lnc) = 0.0
         IF ( lccn > 1 ) an(ix,jy,kz,lccn) = Max( 0.0, an(ix,jy,kz,lccn) )
         
         IF ( lccna > 0 .and. ac_opt == 0  ) THEN ! apply exponential decay to activated CCN to restore to environmental value
           IF ( restoreccn ) THEN
           tmp = an(ix,jy,kz,li) + an(ix,jy,kz,ls)  
           
           IF ( an(ix,jy,kz,lccna) > 1. .and. tmp < qxmin(li) ) an(ix,jy,kz,lccna) = an(ix,jy,kz,lccna)*Exp(-dtp/ccntimeconst)
           ENDIF
         ELSEIF ( lccn > 1 .and. restoreccn .and. ac_opt == 0  ) THEN
           ! in this case, we are treating the ccn field as ccna
           tmp = an(ix,jy,kz,li) + an(ix,jy,kz,ls)  
!           IF ( ny == 2 .and. ix == nx/2 ) THEN
!             write(0,*) 'restore: k, qccn,exp = ',kz,qccn,dn(ix,jy,kz)*qccn,Exp(-dtp/ccntimeconst)
!             write(0,*) 'ccn1,ccn2 = ',an(ix,jy,kz,lccn),dn(ix,jy,kz)*qccn - Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*Exp(-dtp/ccntimeconst)
!           ENDIF
           IF ( an(ix,jy,kz,lccn) > 1. .and. tmp < qxmin(li) .and. ( an(ix,jy,kz,lccn) < dn(ix,jy,kz)*qccn .or. .not. invertccn ) ) THEN 
        !      an(ix,jy,kz,lccn) =  &
        !            an(ix,jy,kz,lccn) +  Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*(1.0 - Exp(-dtp/ccntimeconst))
        ! Equivalent form after expanding last term:
               an(ix,jy,kz,lccn) =  &
                    dn(ix,jy,kz)*qccn - Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*Exp(-dtp/ccntimeconst)
           ENDIF
         
         ENDIF

#ifdef CHGELEC
        IF ( ipelec > 0 .and. lscw .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lscw) .gt. 0.0 ) THEN
            an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscw)*eci
          ELSE
            an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscw)*eci
          ENDIF

           an(ix,jy,kz,lscw) = 0.0

        ENDIF
#endif
       ENDIF

      ENDIF

      end do
!      end do
      end do
      
      
      end subroutine smallvalues
