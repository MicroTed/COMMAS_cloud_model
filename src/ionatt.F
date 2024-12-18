c
c  subroutine to calculate fall speeds of hydrometeors and ion attachment
c
      subroutine ionattdd(nx,ny,nz,nor,na,nba,dtp,dz,jgs,
     :  elecez,attach,uz,diff,cond,conc,
     :  an,dn,db,t0,t7,
     :  qmin)
c
c TO DO LIST:
c
c    Note that ngs=1 and only one point is done at a time.  Eventually
c    may reimplement the full gather-scatter as in the microphysics
c
c
      USE MICRO_MODULE
      USE INDEX_MODULE
      implicit none
      
c
c  include file for mix ratio and charge indices
c
!      include 'sam.index.ion.h'
      
      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag,istag
      parameter (kstag=1, istag=1)
c
      integer ndebug2
      parameter (ndebug2 = 0)

      integer nx,ny,nz,nor,ngt,jgs,na,nba
c      real gt(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,ngt)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)

c
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)


      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
c      real vt(nx,nz)
      real dtp,dz,dtz1
c      real qtmp1(nx,nz),qtmp2(nx,nz)
c      real cmax
c      real xfall(nx,ny)
      

      real       elecez(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      integer idx

      integer ix,jy,kz,i,j,k,item,ia,n, il, ic
      
      integer ngs1
c      real cinccn(nz),
      real cwc1,cwc0

c      real xvcmn, xvcmx  ! min, max droplet volumes
c      real xvrmn, xvrmx  ! min, max rain volumes
c      real xvsmn, xvsmx  ! min, max snow volumes
c      real xvfmn, xvfmx  ! min, max frozen drop volumes
c      real xvgmn, xvgmx  ! min, max graupel volumes
c      real xvhmn, xvhmx  ! min, max hail volumes
c      real xvhlmn, xvhlmx  ! min, max lg hail volumes
      real rwmasn,rwmasx
c      PARAMETER(XVCMN=4.188E-12,XVCMX=6.54E-8)   ! CGS
c      PARAMETER(XVRMN=2.8866E-7,XVRMX=4.1887E-3) ! CGS
c      parameter( xvcmn=4.188e-18,  xvcmx=6.54e-14 )    ! mks
c      parameter( xvrmn=2.8866e-13, xvrmx=4.1887e-9 )  ! mks
c      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(2.e-3)**3 ) !( was 4.1887e-9 )  ! mks
c      parameter( xvsmn=2.8866e-13, xvsmx=0.523599*(3.e-3)**3 ) !( was 4.1887e-9 )  ! mks
c      parameter( xvfmn=2.8866e-13, xvfmx=0.523599*(3.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvgmn=2.8866e-13, xvgmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvhmn=0.523599*(1.e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvhlmn=0.523599*(10.e-3)**3, 
c     :           xvhlmx=0.523599*(100.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      
c      real cwvt(nx,nz)
c      real civt(nx,nz)
c      real cirvt(nx,nz)
c      real cipvt(nx,nz)
c      real rwvt(nx,nz)
c      real fwvt(nx,nz)
c      real swvt(nx,nz)
c      real glvt(nx,nz)
c      real gmvt(nx,nz)
c      real ghvt(nx,nz)
c      real hlvt(nx,nz)

      real attach(nx,nz,lc:lhab,2)

c      real vt(nx,nz,lc:lhab)
      
c      real rad(nx,nz,lc:lhab)
      
      real conc(nx,nz,lc:lhab) ! concentration

      real diff(nx,nz,lc:lhab,2)
      real cond(nx,nz,lc:lhab,2)
 
       real qmin(lqb:lhab)
!       integer lqb1
!       real qmin(lqb1:lqb1+50)

      integer ngs,ngscnt,mgs
      parameter ( ngs=500 )
      integer igs(ngs),kgs(ngs)
      
      real
     :  cwmas(ngs),vtwbar(ngs),
     :  rho0(ngs),temcg(ngs),cip(ngs),cci(ngs),ccw(ngs),
     :  qcw(ngs),qrw(ngs),qci(ngs),qip(ngs),qir(ngs),
     :  qsw(ngs),qgl(ngs),qgm(ngs),
     :  qgh(ngs),qhw(ngs),qhl(ngs),qfw(ngs),qwv(ngs),
     :  qvs(ngs)

      real temg(ngs),ez(ngs)
      
      real rhovt(ngs), cirdn(ngs),cwdn(ngs),rwdn(ngs),
     : swdn(ngs), gldn(ngs),gmdn(ngs),ghdn(ngs),
     : fwdn(ngs), hwdn(ngs), hldn(ngs), cipdn(ngs),cidn(ngs)

      real cirdia(ngs), rwdia(ngs), swdia(ngs),gldia(ngs),
     : gmdia(ngs), ghdia(ngs), hwdia(ngs),hldia(ngs),fwdia(ngs)
      real cipdia(ngs),cidia(ngs)

      real rwrad(ngs)

      real vtirbar(ngs),vtrbar(ngs),vtsbar(ngs),vtglbar(ngs),
     :  vtgmbar(ngs), vtghbar(ngs),vtfbar(ngs),vthbar(ngs),
     :  vthlbar(ngs) ,vtipbar(ngs),vtibar(ngs)
      
!      real vtxbar(ngs,lc:lhab)
      real, allocatable :: vtxbar(:,:,:)
      real, allocatable :: xmas(:,:)
      real, allocatable :: xdn(:,:)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real, allocatable :: xdia(:,:,:)

      real xrad(ngs,lc:lhab)  ! mean diameter
!      real scx(ngs,lc:lhab)   ! charge density of category
      real scppx(ngs,lc:lhab) ! charge per particle
!      real cx(ngs,lc:lhab)    ! number density
!      real qx(ngs,lc:lhab)    ! mixing ratio

      real,allocatable :: qx(:,:) ! qx(ngs,lv:lhab)
      real,allocatable :: cx(:,:) ! cx(ngs,lc:lhab)
      real,allocatable :: scx(:,:) ! scx(ngs,lc:lhab)
      real,allocatable :: xv(:,:)

      real cion(ngs,2)  ! 1=pos ions; 2=neg. ions
      real scxaciond(ngs,lc:lhab,2)
      real scxacionc(ngs,lc:lhab,2)
      
      real diffco(ngs,2)  ! diffusion coefficient for pos (1) and neg (2) ions

      real uz(nz,4)  ! ion mobilities
      double precision s
c      parameter (s(1) = 1.0, s(2) = -1.0)
      
      double precision x,p1,p2,p3
      double precision qd,qm
      
      real kb   ! Boltzman constant
      parameter (kb = 1.3807e-23)
      
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)

      real eperao
      parameter (eperao  = 8.8592e-12 )

      real cir(ngs), crw(ngs) ,csw(ngs), cgl(ngs), cgm(ngs),
     :  cgh(ngs), chw(ngs),chl(ngs), cfw(ngs)
     
      real ciplen(ngs),cilen(ngs)
      
      real cwnc(ngs) ! ,cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
      real cimasn,cimasx,cnina(ngs),cimas(ngs)
      real xvr(ngs),rwmas(ngs)
      real xvs(ngs),xvgl(ngs),xvgm(ngs),xvgh(ngs),xvf(ngs)
      real xvh(ngs),xvhl(ngs)
      real swmas(ngs),glmas(ngs),gmmas(ngs),ghmas(ngs)
      real fwmas(ngs),hwmas(ngs),hlmas(ngs)
      real mwfac

c
c  Ice-Ice: Collection (cxc) efficiencies
c
      real eii(ngs),eipi(ngs),eipip(ngs),eipir(ngs),eirir(ngs),
     :     eiri(ngs),eirw(ngs),ehxw(ngs),ehxr(ngs),ehxi(ngs),
     :     ehxs(ngs)
c
c  Ice-cloud water: Collection (cxc) efficiencies
c
      real eiw(ngs),eipw(ngs)
c
c  Rain: Collection (cxc) efficiencies
c
      real erw(ngs),err(ngs),ers(ngs)
      real erir(ngs),eri(ngs),erip(ngs)
c
c  Snow aggregates: Collection (cxc) efficiencies
c
      real cwrad
      real esw(ngs),esr(ngs),ess(ngs),esi(ngs)
      real esir(ngs),esip(ngs)
c
c  Graupel-low density: Collection (cxc) efficiencies
c
      real eglw(ngs),eglr(ngs),egls(ngs),egli(ngs)
      real eglip(ngs),eglir(ngs)
c
c  Graupel-medium density: Collection (cxc) efficiencies
c
      real egmw(ngs),egmr(ngs), egms(ngs)
      real egmi(ngs),egmip(ngs),egmir(ngs)
c
c  Graupel-high density: Collection (cxc) efficiencies
c
      real eghw(ngs),eghr(ngs), eghs(ngs)
      real eghi(ngs),eghir(ngs),eghip(ngs)
c
c  Frozen drops: Collection (cxc) efficiencies
c
      real efw(ngs),efr(ngs),efs(ngs),efi(ngs),efip(ngs),efir(ngs)
c
c  Hail: Collection (cxc) efficiencies
c
      real ehw(ngs),ehr(ngs),ehs(ngs)
      real ehi(ngs),ehip(ngs),ehir(ngs)
c
c  Hail-large: Collection (cxc) efficiencies
c
      real ehlw(ngs),ehlr(ngs),ehls(ngs),ehli(ngs)
      real ehlip(ngs),ehlir(ngs)

      real poo, cp608, cp, cv,ar
      real aradcw, bradcw, cradcw, dradcw
      real br, bta1, cnit, dragh, dnz00, rho00, cs, ds
      real pi, pii, pid4, qccrit, qscrit
      real gf1, gf1p5, gf2,  gf3, gf3p5, gf4, gf5
      real gf6, gf7, gf4br, gf4ds, gf4p5, gf3ds, gf1ds, gr
      real gamma
c
c intercepts
c
      real cnow, cnoi, cnoip, cnoir 
      real cnogl, cnogm, cnogh
c
c  density maximums and minimums
c
      real rwdnmx, cwdnmx, cidnmx, xidnmx
      real swdnmx, gldnmx, gmdnmx, ghdnmx, fwdnmx, hwdnmx, hldnmx
c
      real rwdnmn, cwdnmn, xidnmn, cidnmn
      real swdnmn, gldnmn, gmdnmn, ghdnmn, fwdnmn, hwdnmn, hldnmn0
c
c  constants
c
      real c1f3
c
c  general constants for microphysics
c
      real brz, arz, rd
      real qcmin, qimin, qrmin, qsmin, qhmin
      real qfmin, qghmin, qgmmin, qglmin, qirmin, qipmin, tdtol, tfr
      real thnuc, rw, advisc0, advisc1, tka0, cpi, cap, tfrcbw, tfrcbi
      real qxmin(lc:lhab)
c
c   drag coefficients
c
      real cdrw, cdhw, cdhl, cdfw, cdgh, cdgm,  cdgl, cdir

c 
c Miscellaneous
c
c      character*80 infile, outfile
c      integer ifile
      integer ihabdo
      parameter(ihabdo = 0)
      integer ierr

      integer lrt(lc:lhab)
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)

      real cdx(lc:lhab)
      real cno(lc:lhab), cnox
      real xdn0(lc:lhab)

      real xvmn(lc:lhab), xvmx(lc:lhab)
      
      real :: tmp, tmp2

      real cwnccn(ngs)

c #####################################################################
c BEGIN EXECUTABLE
c #####################################################################
c


      DO kz = 1,nz-kstag 
      cwnccn(kz) = cwccn*db(kz)/db(1)  
      ENDDO


      lsc(lc) = lscw
      lsc(lr) = lscr
      lsc(li) = lsci
      lsc(lip) = lscip
      lsc(lir) = lscir
      lsc(ls) = lscs
      lsc(lgl) = lscgl
      lsc(lgm) = lscgm
      lsc(lgh) = lscgh
      lsc(lf) = lscf
      lsc(lh) = lsch
      lsc(lhl) = lschl
      

      ln(lc) = lnc
      ln(lr) = lnr
      ln(li) = lni
      ln(lip) = lnip
      ln(lir) = lnir
      ln(ls) = lns
      ln(lgl) = lngl
      ln(lgm) = lngm
      ln(lgh) = lngh
      ln(lf) = lnf
      ln(lh) = lnh
      ln(lhl) = lnhl

      ipc(lc) = 2
      ipc(lr) = 3
      ipc(li) = 1
      ipc(lip) = 1
      ipc(lir) = 1
      ipc(ls) = 4
      ipc(lgl) = 5
      ipc(lgm) = 5
      ipc(lgh) = 5
      ipc(lf) = 5
      ipc(lh) = 5
      ipc(lhl) = 5

      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn
      xvmn(lgl) = xvgmn
      xvmn(lgm) = xvgmn
      xvmn(lgh) = xvgmn
      xvmn(lf) = xvfmn
      xvmn(lh) = xvhmn
      xvmn(lhl) = xvhlmn

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx
      xvmx(ls) = xvsmx
      xvmx(lgl) = xvgmx
      xvmx(lgm) = xvgmx
      xvmx(lgh) = xvgmx
      xvmx(lf) = xvfmx
      xvmx(lh) = xvhmx
      xvmx(lhl) = xvhlmx

      IF ( ipconc .ge. 2 ) THEN
!        cwmasn = xvmn(lc)*1000.
!        cwradn = 1.0e-6
!        cwmasx = xvmx(lc)*1000.
      ENDIF

      call setqxmin(qxmin)
      CALL setcno(cno)

c
c  density maximums and minimums
c
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(lip) =  917.0
      xdnmx(lir) =  400.0
      xdnmx(ls) =  200.0
      xdnmx(lgl) =  400.0
      xdnmx(lgm) =  600.0
      xdnmx(lgh) =  800.0
      xdnmx(lf) =  900.0
      xdnmx(lh) =  900.0
      xdnmx(lhl) =  900.0
c
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(lip) =  100.0
      xdnmn(lir) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lgl) =  200.0
      xdnmn(lgm) =  400.0
      xdnmn(lgh) =  600.0
      xdnmn(lf) =  700.0
      xdnmn(lh) =  700.0
      xdnmn(lhl) =  900.0

      xdn0(lir) = 300.0
      xdn0(lc) = 1000.0
      xdn0(lr) = 1000.0
      xdn0(ls) = 100.0
      
      call setdn10(xdn0)
      
!      xdn0(lgl) = (0.5)*(xdnmn(lgl)+xdnmx(lgl))
!      xdn0(lgm) = (0.5)*(xdnmn(lgm)+xdnmx(lgm))
!      xdn0(lgh) = (0.5)*(xdnmn(lgh)+xdnmx(lgh))
!      xdn0(lf) = (0.5)*(xdnmn(lf)+xdnmx(lf))
!      xdn0(lh) = (0.5)*(xdnmn(lh)+xdnmx(lh))
!      xdn0(lhl) = (0.5)*(xdnmn(lhl)+xdnmx(lhl))

c
c  Set terminal velocities...
c    also set drag coefficients
c
      cdx(lr) = 0.60
      cdx(lh) = 0.45
      cdx(lhl) = 0.45
      cdx(lf) = 0.45
      cdx(lgh) = 0.60
      cdx(lgm) = 0.80
      cdx(lgl) = 0.80
      cdx(lir) = 2.00
      cdx(ls) = 2.00

        rwmasn = xvrmn*1000.
        rwmasx = xvrmx*1000.

c  constants
c
      mwfac = 6.0**(1./3.) 
      poo = 1.0e+05
      cp608 = 0.608
      cp = 1004.0
      cv = 717.0
      ar = 841.99666  
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      br = 0.8
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
      rho00 = 1.225
      cs = 4.83607122
      ds = 0.25
c  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pi = 4.0*atan(1.0)
      pii = 1./pi
      pid4 = pi/4.0 
      qccrit = 2.0e-03
      qscrit = 6.0e-04
      gf1 = 1.0 ! gamma(1.0)
      gf1p5 = 0.8862269255 ! gamma(1.5)
      gf2 = 1.0 ! gamma(2.0)
      gf3 = 2.0 ! gamma(3.0)
      gf3p5 = 3.32335097 ! gamma(3.5)
      gf4 = 6.0 ! gamma(4.0)
      gf5 = 24.0 ! gamma(5.0)
      gf6 = 120.0 ! gamma(6.0)
      gf7 = 720.0 ! gamma(7.0)
      gf4br = gamma(4.0+br)
      gf4ds = gamma(4.0+ds)
      gf4p5 = gamma(4.0+0.5)
      gf3ds = gamma(3.0+ds)
      gf1ds = gamma(1.0+ds)
      gr = 9.8
      cwc0 = pii ! 6.0*pii
      cwc1 = 6.0/(pi*1000.)
c
c  slope intercepts
c
      cnow  = 1.0e+08
      cnoi  = 1.0e+08
      cnoip = 1.0e+08 
      cnoir = 1.0e+08 
      cnor  = 8.0e+06 
      cnos  = 8.0e+06 
      cnogl = 4.0e+05 
      cnogm = 4.0e+05 
      cnogh = 4.0e+05 
!      cnof  = 4.0e+05
      cnoh  = 4.0e+04
!      cnohl = 1.0e+03
c
c  density maximums and minimums
c
      rwdnmx = 1000.0
      cwdnmx = 1000.0
      cidnmx =  917.0
      xidnmx =  917.0
      swdnmx =  200.0
      gldnmx =  400.0
      gmdnmx =  600.0
      ghdnmx =  800.0
      fwdnmx =  900.0
      hwdnmx =  900.0
      hldnmx =  900.0
c
      rwdnmn = 1000.0
      cwdnmn = 1000.0
      xidnmn =  001.0
      cidnmn =  001.0
      swdnmn =  001.0
      gldnmn =  200.0
      gmdnmn =  400.0
      ghdnmn =  600.0
      fwdnmn =  700.0
      hwdnmn =  700.0
      hldnmn0 =  900.0
c
c  constants
c
      c1f3 = 1.0/3.0
c
c  general constants for microphysics
c


      brz = 100.0
      arz = 0.66
      rd = 287.04
      cp = 1004.0
      
c      qxmin = 1.0e-09
c      qcmin = qmin(lc) ! 1.0e-09
c      qimin = qmin(li) ! 1.0e-09
c      qrmin =  qmin(lr) !1.0e-07
c      qsmin =  qmin(ls) !1.0e-09
c      qhmin = qmin(lh) ! 1.0e-07
c      qfmin = qmin(lf) ! 1.0e-07
c      qghmin = qmin(lgh) ! 1.0e-07
c      qgmmin = qmin(lgm) ! 1.0e-07
c      qglmin = qmin(lgl) ! 1.0e-07
c      qirmin = qmin(lir) ! 1.0e-08
c      qipmin = qmin(lip) ! 1.0e-09

c
c  set minimum mass mixing ratios
c
c      qmin(lc) = 1.0e-09
c      qmin(li) = 1.0e-09
c      qmin(lr) = 1.0e-07
c      qmin(ls) = 1.0e-07
c      qmin(lf) = 1.0e-07
c      qmin(lgl) = 1.0e-07
c      qmin(lgm) = 1.0e-07
c      qmin(lgh) = 1.0e-07
c      qmin(lh) = 1.0e-07
c      qmin(lhl) = 1.0e-07
c      qmin(lir) = 1.0e-09
c      qmin(lip) = 1.0e-09


      tdtol = 1.0e-05
      tfr = 273.15
      thnuc = 233.15
      rw = 461.5
      advisc0 = 1.832e-05
      advisc1 = 1.718e-05
      tka0 = 2.43e-02
      cpi = 1.0/cp
      cap = rd/cp 
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
c
c  ci constants in mks units
c
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
c
c  Set terminal velocities...
c    also set drag coefficients
c
      cdrw = 0.60
      cdhw = 0.45
      cdhl = 0.45
      cdfw = 0.45
      cdgh = 0.60
      cdgm = 0.80
      cdgl = 0.80
      cdir = 2.00
      
      jy = jgs

      DO ia=lc,lhab
       DO kz = 1,nz-2
        DO ix = 1,nx-1
c          vt(ix,kz,ia) = 0.0
c          rad(ix,kz,ia) = 0.0
          conc(ix,kz,ia) = 0.0
          attach(ix,kz,ia,1) = 0.0
          attach(ix,kz,ia,2) = 0.0
        ENDDO ! ia
       ENDDO ! kz
      ENDDO ! ix

c      DO kz = 1,nz-2
c        DO ix = 1,nx-1

      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do 1000 inumgs = 1,numgs 
      ngscnt = 0
      do 1001 kz = nzmpb,nz-kstag-1 
      do 1002 ix = nxmpb,nx-istag

      if ( 
     >     an(ix,jy,kz,lc)  .gt. qmin(lc)   .or. 
     >     an(ix,jy,kz,li)  .gt. qmin(li)   .or.
     >     an(ix,jy,kz,lip) .gt. qmin(lip)   .or.
     >     an(ix,jy,kz,lir) .gt. qmin(lir)   .or.
     >     an(ix,jy,kz,lr)  .gt. qmin(lr)   .or.
     >     an(ix,jy,kz,ls)  .gt. qmin(ls)   .or.
     >     an(ix,jy,kz,lgl) .gt. qmin(lgl)   .or.
     >     an(ix,jy,kz,lgm) .gt. qmin(lgm)   .or.
     >     an(ix,jy,kz,lgh) .gt. qmin(lgh)   .or.
     >     an(ix,jy,kz,lf)  .gt. qmin(lf)   .or.
     >     an(ix,jy,kz,lhl) .gt. qmin(lh)   .or.
     <     an(ix,jy,kz,lh)  .gt. qmin(lhl) ) then
       ! load temp quantities
       
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 1100
      end if
 1002 continue
      nxmpb = 1
 1001 continue
c      if ( jy .eq. (ny-jstag) ) iend = 1
 1100 continue
      if ( ngscnt .eq. 0 ) go to 9998

c      ngscnt = ngscnt + 1 
c      mgs = 1
c      igs(mgs) = ix
c      kgs(mgs) = kz


      allocate (  qx(ngscnt,lv:lhab) )
      allocate (  cx(ngscnt,lc:lhab) )
      allocate (  scx(ngscnt,lc:lhab) )
      allocate (  xv(ngscnt,lc:lhab) )
      allocate (  vtxbar(ngscnt,lc:lhab,2) )
      allocate (  xmas(ngscnt,lc:lhab) )
      allocate (  xdn(ngscnt,lc:lhab) )
      allocate (  xdia(ngscnt,lc:lhab,3) )
      
      xv(:,:)       = 0.0
      xmas(:,:)     = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:)   = 0.0
      scx(:,:)      = 0.0

c
c  set temporaries for microphysics variables
c
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

      DO il = lc,lhab
      do mgs = 1,ngscnt
      scx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lsc(il))
      ENDDO
      ENDDO

      
      do mgs = 1,ngscnt
     
C      qip(mgs) = max(an(igs(mgs),jy,kgs(mgs),lip), 0.0) 
C      qir(mgs) = max(an(igs(mgs),jy,kgs(mgs),lir), 0.0) 
C      qcw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lc), 0.0) 
C      qrw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lr), 0.0) 
C      qsw(mgs) = max(an(igs(mgs),jy,kgs(mgs),ls), 0.0) 
C      qgl(mgs) = max(an(igs(mgs),jy,kgs(mgs),lgl), 0.0) 
C      qgm(mgs) = max(an(igs(mgs),jy,kgs(mgs),lgm), 0.0) 
C      qgh(mgs) = max(an(igs(mgs),jy,kgs(mgs),lgh), 0.0) 
C      qfw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lf), 0.0) 
C      qhw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lh), 0.0) 
C      qhl(mgs) = max(an(igs(mgs),jy,kgs(mgs),lhl), 0.0) 
C      qci(mgs) = max(an(igs(mgs),jy,kgs(mgs),li), 0.0) 
C      
C      scx(mgs,li)  = an(igs(mgs),jy,kgs(mgs),lsci)
C      scx(mgs,lip) = an(igs(mgs),jy,kgs(mgs),lscip)
C      scx(mgs,lir) = an(igs(mgs),jy,kgs(mgs),lscir)
C      scx(mgs,lc)  = an(igs(mgs),jy,kgs(mgs),lscw)
C      scx(mgs,lr)  = an(igs(mgs),jy,kgs(mgs),lscr)
C      scx(mgs,lf)  = an(igs(mgs),jy,kgs(mgs),lscf)
C      scx(mgs,ls)  = an(igs(mgs),jy,kgs(mgs),lscs)
C      scx(mgs,lgl) = an(igs(mgs),jy,kgs(mgs),lscgl)
C      scx(mgs,lgm) = an(igs(mgs),jy,kgs(mgs),lscgm)
C      scx(mgs,lgh) = an(igs(mgs),jy,kgs(mgs),lscgh)
C      scx(mgs,lh)  = an(igs(mgs),jy,kgs(mgs),lsch)
C      scx(mgs,lhl) = an(igs(mgs),jy,kgs(mgs),lschl)
      
      cion(mgs,1) = an(igs(mgs),jy,kgs(mgs),lscpi)
      cion(mgs,2) = an(igs(mgs),jy,kgs(mgs),lscni)
      
      ez(mgs) = elecez(igs(mgs),jy,kgs(mgs))


      end do

      cx(:,:) = 0.0
      
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
        cx(mgs,lir) = Max(an(igs(mgs),jy,kgs(mgs),lnir), 0.0)
        cx(mgs,lip) = Max(an(igs(mgs),jy,kgs(mgs),lnip), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
       end do
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
          cx(mgs,lr) = 0.0
        ELSEIF ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) THEN
          qx(mgs,lr) = 0.0
        ELSE
          cx(mgs,lr) = Max( 1.e-3, cx(mgs,lr) )
        ENDIF
       end do
      end if
      if ( ipconc .ge. 4 ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
          cx(mgs,ls) = 0.0
        ELSEIF ( cx(mgs,ls) .eq. 0.0 .and. qx(mgs,ls) .lt. 3.0*qxmin(ls) ) THEN
          qx(mgs,ls) = 0.0
        ELSE
          cx(mgs,ls) = Max( 1.e-3, cx(mgs,ls) )

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
        cx(mgs,lgl) = Max(an(igs(mgs),jy,kgs(mgs),lngl), 0.0)
        IF ( qx(mgs,lgl) .le. qxmin(lh) ) THEN
          cx(mgs,lgl) = 0.0
        ELSEIF ( cx(mgs,lgl) .eq. 0.0 .and. qx(mgs,lgl) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgl) = 0.0
        ELSE
          cx(mgs,lgl) = Max( 1.e-3, cx(mgs,lgl) )
c test check for intercept set max of 4.e6
            
         IF ( ilimit .ge. ipc(lgl) ) THEN
            tmp = (xdn0(lgl)*cx(mgs,lgl))/(rho0(mgs)*qx(mgs,lgl))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgl)*(tmp2)
c           IF ( qx(mgs,lgl) .gt. 0.5e-3 ) THEN
c           write(6,*) 'nogl1: ',igs(mgs),kgs(mgs),cx(mgs,lgl),rho0(mgs),
c     :      Log10(cnox),cnox*1.e-6,tmp2,qx(mgs,lgl)*1000.
c           ENDIF
         IF ( cnox .gt. 3.0*cno(lgl) ) THEN
           cx(mgs,lgl) = 3.0*cno(lgl)/tmp2
c           IF ( qx(mgs,lgl) .gt. 0.5e-3 ) THEN
c           write(6,*) 'nogl2: ',cx(mgs,lgl)
c           ENDIF
         ENDIF
         ENDIF
        ENDIF

        cx(mgs,lgm) = Max(an(igs(mgs),jy,kgs(mgs),lngm), 0.0)
        IF ( qx(mgs,lgm) .le. qxmin(lh) ) THEN
          cx(mgs,lgm) = 0.0
        ELSEIF ( cx(mgs,lgm) .eq. 0.0 .and. qx(mgs,lgm) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgm) = 0.0
        ELSE
          cx(mgs,lgm) = Max( 1.e-3, cx(mgs,lgm) )

         IF ( ilimit .ge. ipc(lgm) ) THEN
            tmp = (xdn0(lgm)*cx(mgs,lgm))/(rho0(mgs)*qx(mgs,lgm))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgm)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lgm) ) THEN
           cx(mgs,lgm) = 3.0*cno(lgm)/tmp2
         ENDIF
         ENDIF
        ENDIF

        cx(mgs,lgh) = Max(an(igs(mgs),jy,kgs(mgs),lngh), 0.0)
        IF ( qx(mgs,lgh) .le. qxmin(lh) ) THEN
          cx(mgs,lgh) = 0.0
        ELSEIF ( cx(mgs,lgh) .eq. 0.0 .and. qx(mgs,lgh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgh) = 0.0
        ELSE
          cx(mgs,lgh) = Max( 1.e-3, cx(mgs,lgh) )
         IF ( ilimit .ge. ipc(lf) ) THEN
            tmp = (xdn0(lgh)*cx(mgs,lgh))/(rho0(mgs)*qx(mgs,lgh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgh)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lgh) ) THEN
           cx(mgs,lgh) = 3.0*cno(lgh)/tmp2
         ENDIF
         ENDIF
        ENDIF

        cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)
        IF ( qx(mgs,lf) .le. qxmin(lh) ) THEN
          cx(mgs,lf) = 0.0
        ELSEIF ( cx(mgs,lf) .eq. 0.0 .and. qx(mgs,lf) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lf) = 0.0
        ELSE
          cx(mgs,lf) = Max( 1.e-3, cx(mgs,lf) )
        ENDIF

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
          cx(mgs,lh) = 0.0
        ELSEIF ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lh) = 0.0
        ELSE
          cx(mgs,lh) = Max( 1.e-3, cx(mgs,lh) )
         IF ( ilimit .ge. ipc(lh) ) THEN
            tmp = (xdn0(lh)*cx(mgs,lh))/(rho0(mgs)*qx(mgs,lh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lh)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lh) ) THEN
           cx(mgs,lh) = 3.0*cno(lh)/tmp2
         ENDIF
         ENDIF
        ENDIF

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        IF ( qx(mgs,lhl) .le. qxmin(lh) ) THEN
          cx(mgs,lhl) = 0.0
        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lhl) = 0.0
        ELSE
          cx(mgs,lhl) = Max( 1.e-3, cx(mgs,lhl) )
        ENDIF
       end do
      end if
        
       
c
c
c  Reconstruct various quantities 
c

      
      do mgs = 1,ngscnt
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhovt(mgs) = Sqrt(rho00/rho0(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
c      temg(mgs) = temcg(mgs) + tfr
      end do
c

c
c  Set density
c
      if (ndebug .gt. 0 ) print*,'Set density'
c
      DO mgs=1,ngscnt
      cirdn(mgs) = 300.0
      cwdn(mgs) = 1000.0
      rwdn(mgs) = 1000.0
      swdn(mgs) = 100.0
      gldn(mgs) = (0.5)*(gldnmn+gldnmx)
      gmdn(mgs) = (0.5)*(gmdnmn+gmdnmx)
      ghdn(mgs) = (0.5)*(ghdnmn+ghdnmx)
      fwdn(mgs) = (0.5)*(fwdnmn+fwdnmx)
      hwdn(mgs) = (0.5)*(hwdnmn+hwdnmx)
      hldn(mgs) = (0.5)*(hldnmn0+hldnmx)
      cipdn(mgs)  = 900.0
      END DO
c
c  Set density
c
      if (ndebug .gt. 0 ) print*,'Set density'
c

      do mgs = 1,ngscnt
        xdn(mgs,lir) = xdn0(lir)
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
        xdn(mgs,ls) = xdn0(ls)
        xdn(mgs,lgl) = xdn0(lgl)
        xdn(mgs,lgm) = xdn0(lgm)
        xdn(mgs,lgh) = xdn0(lgh)
        xdn(mgs,lf) = xdn0(lf)
        xdn(mgs,lh) = xdn0(lh)
        xdn(mgs,lhl) = xdn0(lhl)
      end do

c
c  Set a couple of cloud variables...
c
      do mgs = 1,ngscnt
      cwnc(mgs) = cwccn ! cwnccn(kgs(mgs))
c      cinc(mgs) = cinccn(kgs(mgs))
      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*
     >  (temg(mgs)/296.0)**(1.5)
c
      end do
c

      call setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
     :                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,
     :                 ipconc,ndebug,ngs,nz,kgs,cwnccn,fadvisc,
     :                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,
     :                 itype1,itype2,temcg,0)


c
c  Set number concentrations (need xdia from setvt)
c
      if ( ndebug .gt. 0 ) print*,'Set concentration'
      if ( ipconc .lt. 5 ) then
      do mgs = 1,ngscnt

      IF ( ipconc .eq. 0 ) THEN
      cx(mgs,lir) = 0.0
      if ( qx(mgs,lir) .gt. qxmin(li) )  then
      cx(mgs,lir) = cno(lir)*xdia(mgs,lir,1)
      end if
      ENDIF

      IF ( ipconc .lt. 3 ) THEN
      cx(mgs,lr) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) )  then
      cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
      end if
      ENDIF

      IF ( ipconc .lt. 4 ) THEN
      tmp = cx(mgs,ls)
      cx(mgs,ls) = 0.0
      xv(mgs,ls) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) )  then
      cx(mgs,ls) = cno(ls)*xdia(mgs,ls,1)
      xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
      end if
      ENDIF ! ( ipconc .lt. 4 )

      IF ( ipconc .lt. 5 ) THEN 

      cx(mgs,lgl) = 0.0
      xv(mgs,lgl) = 0.0
      if ( qx(mgs,lgl) .gt. qxmin(lh) )  then
      cx(mgs,lgl) = cno(lgl)*xdia(mgs,lgl,1)
      xv(mgs,lgl) = rho0(mgs)*qx(mgs,lgl)/(xdn(mgs,lgl)*cx(mgs,lgl))
      end if

      cx(mgs,lgm) = 0.0
      xv(mgs,lgm) = 0.0
      if ( qx(mgs,lgm) .gt. qxmin(lh) )  then
      cx(mgs,lgm) = cno(lgm)*xdia(mgs,lgm,1)
      xv(mgs,lgm) = rho0(mgs)*qx(mgs,lgm)/(xdn(mgs,lgm)*cx(mgs,lgm))
      end if

      cx(mgs,lgh) = 0.0
      xv(mgs,lgh) = 0.0
      if ( qx(mgs,lgh) .gt. qxmin(lh) )  then
      cx(mgs,lgh) = cno(lgh)*xdia(mgs,lgh,1)
      xv(mgs,lgh) = rho0(mgs)*qx(mgs,lgh)/(xdn(mgs,lgh)*cx(mgs,lgh))
      end if

      cx(mgs,lf) = 0.0
      xv(mgs,lf) = 0.0
      if ( qx(mgs,lf) .gt. qxmin(lh) )  then
      cx(mgs,lf) = cno(lf)*xdia(mgs,lf,1)
      xv(mgs,lf) = rho0(mgs)*qx(mgs,lf)/(xdn(mgs,lf)*cx(mgs,lf))
      end if

      cx(mgs,lh) = 0.0
      xv(mgs,lh) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) )  then
      cx(mgs,lh) = cno(lh)*xdia(mgs,lh,1)
      xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))
      end if

      cx(mgs,lhl) = 0.0
      xv(mgs,lhl) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lh) )  then
      cx(mgs,lhl) = cno(lhl)*xdia(mgs,lhl,1)
      xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))
      end if

      ENDIF ! ( ipconc .lt. 5 )

      end do
      end if
      
      
c  Set a couple of cloud variables...
c
c      do mgs = 1,ngscnt
c      cwnc(mgs) = cwnccn(kgs(mgs))
c      cinc(mgs) = cinccn(kgs(mgs))
c      end do
c
c  cloud water variables
c
c cloud ice variables
c plates
c
c      if ( ndebug .gt. 0 )  print*,'Set ci'
c      if ( ihabdo .lt. 1 ) then
c      do mgs = 1,ngscnt
cc      ciplen(mgs) = 1.25e-3*(cipdia(mgs)**(0.4740))
cc      cipdia2(mgs) = cipdia(mgs)**2
c      vtipbar(mgs) = (131.6)*(cipdia(mgs)**(0.824))
c      vtipbar(mgs) = vtipbar(mgs)*rhovt(mgs)
c      end do
c      endif
cc
cc cloud ice variables
cc columns
cc
c      if ( ndebug .gt. 0 ) print*,'Set cip'
c      if ( ihabdo .lt. 1 ) then
c      do mgs = 1,ngscnt
c      cilen(mgs)   = 0.4764*(cidia(mgs))**(0.958)
c      vtibar(mgs) = (4.942e4)*(cidia(mgs)**(1.4150))
c      cidn(mgs)   = min(max(769.8*cidia(mgs)**(-0.0140),300.0),900.0)
cc      cidia2(mgs) = cidia(mgs)**2
c      vtibar(mgs) = vtibar(mgs)*rhovt(mgs)
cc      cicap(mgs) = 0.0
c      end do
c      endif

      do mgs = 1,ngscnt
            
C      vtxbar(mgs,li)  = vtibar(mgs)
C      vtxbar(mgs,lip) = vtipbar(mgs)
C      vtxbar(mgs,lir) = vtirbar(mgs)
C      vtxbar(mgs,lc)  = vtwbar(mgs)
C      vtxbar(mgs,lr)  = vtrbar(mgs)
C      vtxbar(mgs,lf)  = vtfbar(mgs)
C      vtxbar(mgs,ls)  = vtsbar(mgs)
C      vtxbar(mgs,lgl) = vtglbar(mgs)
C      vtxbar(mgs,lgm) = vtgmbar(mgs)
C      vtxbar(mgs,lgh) = vtghbar(mgs)
C      vtxbar(mgs,lh)  = vthbar(mgs)
C      vtxbar(mgs,lhl) = vthlbar(mgs)
C
C
C      cx(mgs,li) =  cci(mgs)
C      cx(mgs,lip) = cip(mgs)
C      cx(mgs,lir) = cir(mgs)
C      cx(mgs,lc) =  ccw(mgs)
C      cx(mgs,lr) =  crw(mgs)
C      cx(mgs,lf) =  cfw(mgs)
C      cx(mgs,ls) =  csw(mgs)
C      cx(mgs,lgl) = cgl(mgs)
C      cx(mgs,lgm) = cgm(mgs)
C      cx(mgs,lgh) = cgh(mgs)
C      cx(mgs,lh) =  chw(mgs)
C      cx(mgs,lhl) = chl(mgs)
C
C      qx(mgs,li)  = qci(mgs)
C      qx(mgs,lip) = qip(mgs)
C      qx(mgs,lir) = qir(mgs)
C      qx(mgs,lc)  = qcw(mgs)
C      qx(mgs,lr)  = qrw(mgs)
C      qx(mgs,lf)  = qfw(mgs)
C      qx(mgs,ls)  = qsw(mgs)
C      qx(mgs,lgl) = qgl(mgs)
C      qx(mgs,lgm) = qgm(mgs)
C      qx(mgs,lgh) = qgh(mgs)
C      qx(mgs,lh)  = qhw(mgs)
C      qx(mgs,lhl) = qhl(mgs)

      DO il=lc,lhab
        xrad(mgs,il) = 0.5*xdia(mgs,il,1)
      ENDDO

C      xrad(mgs,li) =  cidia(mgs)*0.5
C      xrad(mgs,lip) = cipdia(mgs)*0.5
C      xrad(mgs,lir) = cirdia(mgs)*0.5
C      xrad(mgs,lc) =  cwdia(mgs)*0.5
C      xrad(mgs,lr) =  rwdia(mgs)*0.5
C      xrad(mgs,lf) =  fwdia(mgs)*0.5
C      xrad(mgs,ls) =  swdia(mgs)*0.5
C      xrad(mgs,lgl) = gldia(mgs)*0.5
C      xrad(mgs,lgm) = gmdia(mgs)*0.5
C      xrad(mgs,lgh) = ghdia(mgs)*0.5
C      xrad(mgs,lh) =  hwdia(mgs)*0.5
C      xrad(mgs,lhl) = hldia(mgs)*0.5
      
CVD$L SKIP
      DO ia=lc,lhab
        scppx(mgs,ia) = 0.0
      ENDDO
      
      end do
      
      
c
c ion diffusion/conduction
c      
      do mgs = 1,ngscnt
c        diffco(mgs,1) = 2.0*pi*uz(kgs(mgs),1)*kb*temg(mgs)/ec
c        diffco(mgs,2) = 2.0*pi*uz(kgs(mgs),2)*kb*temg(mgs)/ec
        diffco(mgs,1) = 5.415227187e-4*uz(kgs(mgs),1)*temg(mgs)
        diffco(mgs,2) = 5.415227187e-4*uz(kgs(mgs),2)*temg(mgs)
        
          DO ia=lc,lhab
            IF ( qx(mgs,ia) .gt. qmin(ia) ) THEN
            DO n=1,2
             scxaciond(mgs,ia,n) = 0.0
             scxacionc(mgs,ia,n) = 0.0
             s = 3 - 2*n   ! 1 for n=1, -1 for n=2

c
c diffusion
c            
c            qd = 4.0*pi*eperao*xrad(mgs,ia)*kb*temg(mgs)/ec
            qd = 9.59491613906812e-15*xrad(mgs,ia)*temg(mgs)
            
            IF ( qd .gt. 1.0e-20 .and. cx(mgs,ia) .gt. 1.0e-6 ) THEN
            scppx(mgs,ia) = scx(mgs,ia)/cx(mgs,ia)
            
            x = scppx(mgs,ia)/qd
            
            p1 = 2.0*xrad(mgs,ia)*diffco(mgs,n)*cion(mgs,n)
            
            IF ( s*x .gt. 1.0e-5 .and. (s*x) .lt. 40.0 ) THEN
c              p2 = Min(1.0, s*x/(Exp(s*x) - 1))
              p2 = Min(1.0D0, s*x/(Exp(s*x) - 1.0D0))
            ELSEIF ( s*x .ge. 40.0 ) THEN
              p2 = 0.0
            ELSE
              p2 = 1
            ENDIF
            
            p3 = 1.0 + Sqrt(xrad(mgs,ia)*vtxbar(mgs,ia,1)/diffco(mgs,n))
            
            scxaciond(mgs,ia,n) = p1*p2*p3
            ELSE
             write(6,*) 'Huh? ia,qd,cx,q = ',ia,qd,cx(mgs,ia),qx(mgs,ia)
             write(6,*) xrad(mgs,ia),temg(mgs),vtxbar(mgs,ia,1),qmin(ia)
             scxaciond(mgs,ia,n) = 0.0
            ENDIF
            
            ENDDO ! n
c
c  conduction
c
            qm = 12.0*pi*eperao*Abs(ez(mgs))*xrad(mgs,ia)**2
            
              IF ( scppx(mgs,ia) .ge. qm ) THEN
                scxacionc(mgs,ia,1) = 0.0
                scxacionc(mgs,ia,2) = cion(mgs,2)
     :            *uz(kgs(mgs),2)*scppx(mgs,ia)/eperao
              
              ELSEIF ( scppx(mgs,ia) .le. -qm ) THEN
                scxacionc(mgs,ia,1) = -cion(mgs,1)
     :            *uz(kgs(mgs),1)*scppx(mgs,ia)/eperao
                scxacionc(mgs,ia,2) = 0.0
              
              ELSEIF ( ez(mgs) .lt. 0.0 ) THEN ! E parallel to vt
                
                scxacionc(mgs,ia,2) = 3.0*pi*xrad(mgs,ia)**2
     :            *cion(mgs,2)*Abs(ez(mgs))
     :            *uz(kgs(mgs),2)*(1.0 + scppx(mgs,ia)/qm)**2
                
                IF ( vtxbar(mgs,ia,1) .gt. uz(kgs(mgs),1)*Abs(ez(mgs)) )
     :                   THEN
                  IF ( scppx(mgs,ia) .gt. 0.0 ) THEN
                    scxacionc(mgs,ia,1) = 0.0
                  ELSE
                    scxacionc(mgs,ia,1) = -uz(kgs(mgs),1)*cion(mgs,1)*
     :                scppx(mgs,ia)/eperao
                  ENDIF
                ELSE  
                scxacionc(mgs,ia,1) = 3.0*pi*xrad(mgs,ia)**2
     :            *cion(mgs,1)*Abs(ez(mgs))
     :            *uz(kgs(mgs),1)*(1.0 - scppx(mgs,ia)/qm)**2
                  
                ENDIF
              ELSEIF ( ez(mgs) .gt. 0.0 ) THEN ! E anti-parallel to vt

                scxacionc(mgs,ia,1) = 3.0*pi*xrad(mgs,ia)**2
     :            *cion(mgs,1)*Abs(ez(mgs))
     :            *uz(kgs(mgs),1)*(1.0 - scppx(mgs,ia)/qm)**2

                IF ( vtxbar(mgs,ia,1) .gt. uz(kgs(mgs),2)*Abs(ez(mgs)) )
     :               THEN
                  IF ( scppx(mgs,ia) .lt. 0.0 ) THEN
                    scxacionc(mgs,ia,2) = 0.0
                  ELSE
                    scxacionc(mgs,ia,2) = uz(kgs(mgs),2)*cion(mgs,2)*
     :                scppx(mgs,ia)/eperao
                  ENDIF
                ELSE  
                scxacionc(mgs,ia,2) = 3.0*pi*xrad(mgs,ia)**2
     :            *cion(mgs,2)*Abs(ez(mgs))
     :            *uz(kgs(mgs),2)*(1.0 + scppx(mgs,ia)/qm)**2
                  
                ENDIF
              
              ELSE 
                scxacionc(mgs,ia,1) = 0.0
                scxacionc(mgs,ia,2) = 0.0
              ENDIF
            ELSE
            DO n=1,2
             scxaciond(mgs,ia,n) = 0.0
             scxacionc(mgs,ia,n) = 0.0
            ENDDO ! n
           ENDIF ! q .gt. qmin
          ENDDO ! ia -- diffusion
          
       
      end do
      
      
c
c load some stuff (may not be needed...)
c
      
      do mgs = 1,ngscnt
c      vt(igs(mgs),kgs(mgs),li) = vtibar(mgs)
c      vt(igs(mgs),kgs(mgs),lip) = vtipbar(mgs)
c      vt(igs(mgs),kgs(mgs),lir) = vtirbar(mgs)
c      vt(igs(mgs),kgs(mgs),lc) = vtwbar(mgs)
c      vt(igs(mgs),kgs(mgs),lr) = vtrbar(mgs)
c      vt(igs(mgs),kgs(mgs),lf) = vtfbar(mgs)
c      vt(igs(mgs),kgs(mgs),ls) = vtsbar(mgs)
c      vt(igs(mgs),kgs(mgs),lgl) = vtglbar(mgs)
c      vt(igs(mgs),kgs(mgs),lgm) = vtgmbar(mgs)
c      vt(igs(mgs),kgs(mgs),lgh) = vtghbar(mgs)
c      vt(igs(mgs),kgs(mgs),lh) = vthbar(mgs)
c      vt(igs(mgs),kgs(mgs),lhl) = vthlbar(mgs)

!      conc(igs(mgs),kgs(mgs),li)  = cci(mgs)
!      conc(igs(mgs),kgs(mgs),lip) = cip(mgs)
!      conc(igs(mgs),kgs(mgs),lir) = cir(mgs)
!      conc(igs(mgs),kgs(mgs),lc)  = ccw(mgs)
!      conc(igs(mgs),kgs(mgs),lr)  = crw(mgs)
!      conc(igs(mgs),kgs(mgs),lf)  = cfw(mgs)
!      conc(igs(mgs),kgs(mgs),ls)  = csw(mgs)
!      conc(igs(mgs),kgs(mgs),lgl) = cgl(mgs)
!      conc(igs(mgs),kgs(mgs),lgm) = cgm(mgs)
!      conc(igs(mgs),kgs(mgs),lgh) = cgh(mgs)
!      conc(igs(mgs),kgs(mgs),lh)  = chw(mgs)
!      conc(igs(mgs),kgs(mgs),lhl) = chl(mgs)

      
c      rad(igs(mgs),kgs(mgs),li) =  cidia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lip) = cipdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lir) = cirdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lc) =  cwdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lr) =  rwdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lf) =  fwdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),ls) =  swdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lgl) = gldia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lgm) = gmdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lgh) = ghdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lh) =  hwdia(mgs)*0.5
c      rad(igs(mgs),kgs(mgs),lhl) = hldia(mgs)*0.5
      
      DO n=1,2
      DO ia=lc,lhab
        attach(igs(mgs),kgs(mgs),ia,n) = cx(mgs,ia)*
     :     (scxaciond(mgs,ia,n) + scxacionc(mgs,ia,n)) 
      IF ( attach(igs(mgs),kgs(mgs),ia,n) .ne. 0.0 .and.
     :      ndebug2 .ge. 1  ) THEN
       write(6,*) 'Attach: ',ia,igs(mgs),kgs(mgs),cx(mgs,ia),
     :   qx(mgs,ia), attach(igs(mgs),kgs(mgs),ia,n),n,
     :  scxaciond(mgs,ia,n),scxacionc(mgs,ia,n)
      ENDIF

      IF ( .not.  ( attach(igs(mgs),kgs(mgs),ia,n) .lt. 1.0e30 .and.
     :     attach(igs(mgs),kgs(mgs),ia,n) .gt. -1.0e30 )) THEN
       write(6,*) 'Attach: ',ia,igs(mgs),kgs(mgs),cx(mgs,ia),
     :   qx(mgs,ia), attach(igs(mgs),kgs(mgs),ia,n),n,
     :  scxaciond(mgs,ia,n),scxacionc(mgs,ia,n)
      write(6,*) 'xdia,vt = ',xdia(mgs,ia,1),vtxbar(mgs,ia,1),xv(mgs,ia)
      write(6,*) 'cno,xdn = ',cno(ia),xdn(mgs,ia)
        STOP
      ENDIF
        
          diff(igs(mgs),kgs(mgs),ia,n) = cx(mgs,ia)*scxaciond(mgs,ia,n)
          cond(igs(mgs),kgs(mgs),ia,n) = cx(mgs,ia)*scxacionc(mgs,ia,n)

      ENDDO
      ENDDO
      
      end do

c
c
c


 123  CONTINUE

      deallocate ( qx )
      deallocate ( cx )
      deallocate ( scx )
      deallocate ( xv )
      deallocate ( vtxbar )
      deallocate ( xmas )
      deallocate ( xdn )
      deallocate ( xdia )
      

 9998 continue
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
c
      go to 1200
      else
      nzmpb = kz 
      end if
      if ( ix+1 .gt. nx-1 ) then
      nxmpb = 1
      else
      nxmpb = ix+1
      end if
 1000 continue
 1200 continue

c       ENDDO ! ix
c      ENDDO ! kz
      
      
      
      RETURN
      END

