#define USEMIXEDPHASE
!
! ##############################################################################
!
      SUBROUTINE setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno, &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,      &
     &                 ipconc1,ndebug1,ngs,nz,igs,kgs,cwnccn,fadvisc,   &
     &                 cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
     &                 itype1a,itype2a,temcg,infdo,alpha,alphan,axx,bxx,ildo)
!     &                 itype1a,itype2a,temcg,infdo,alpha,axh,bxh,axhl,bxhl)
      
      USE COMMASMPI_MODULE, only: commasmpi_abort
      USE INDEX_MODULE, only: lt,lc,lr,li,lis,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lfw,lhw,lhlw, &
                               lnc,lnr,lni,lns,lnh,lnf,lnhl,cinu,dmuh,dnu,xnu,xmu,dmuhl,snu,alphamax,alphahl
      USE MICRO_MODULE
      
      implicit none
      
!      include 'sam.index.ion.h'
!      include 'swm.index.zieg.h'
      
      integer ngscnt,ngs,nz
!      integer infall    ! whether to calculate number-weighted fall speeds
      
      real xv(ngs,lc:lhab)
      real qx(ngs,lv:lhab)
      real qxw(ngs,ls:lhab)
      real cx(ngs,lc:lhab)
      real vtxbar(ngs,lc:lhab,3)
      real xmas(ngs,lc:lhab)
      real xdn(ngs,lc:lhab)
      real xdia(ngs,lc:lhab,3)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real qxmin(lc:lhab)
      real cdx(lc:lhab)
      real cdxgs(ngs,lc:lhab)
      real alpha(ngs,lc:lhab)
      real alphan(ngs,lc:lhab)
      
      real rho0(ngs),rhovt(ngs),temcg(ngs)
      real cno(lc:lhab)
      
      real cwc1, cwnccn(nz), cimna, cimxa
      real cnina(ngs)
      integer igs(ngs),kgs(ngs)
      real fadvisc(ngs)
      real fsw
      
      integer ipconc1
      integer ndebug1
      
      integer itype1a,itype2a,infdo
      real :: axx(ngs,lh:lhab), bxx(ngs,lh:lhab)
      integer, intent(in) :: ildo
!      real :: axh(ngs),bxh(ngs)
!      real :: axhl(ngs),bxhl(ngs)
      
! Local vars

      
      real swmasmx, dtmp
      real cd
      real cwc0 ! ,cwc1
      real :: cwch(ngs), cwchl(ngs), cwcf(ngs)
      real :: cwchtmp,cwcftmp,cwchltmp,xnutmp
      real cimasx,cimasn
      real cwmasn,cwmasx,cwradn
      real cwrad
      real vr,rnux
      real alp
      real myfac(ngs,lc:lhab)
      
      real ccimx

      integer mgs
      
      real arx,frx,vtrain,fw,vtrainnum,vtrainz
      real fwlo,fwhi,rfwdiff
      real ar,br,cs,ds
      real gf4p5, gf4ds, gf4br, ifirst, gf1ds
      real gfcinu1, gfcinu1p47, gfcinu2p47
      real gfcinu1p22,gfcinu2p22
      real gfcinu1p18,gfcinu2p18
      real gsnow1, gsnow53, gsnow73
      real gr
      real rwrad,rwdia
      real frac
      real gamma, gamma_sp
      real mwfac
      integer il
      real alphahltmp

!      save gf4p5, gf4ds, gf4br, ifirst, gf1ds
!      save gfcinu1, gfcinu1p47, gfcinu2p47
!      data ifirst /0/
      
      real bta1,cnit
      parameter ( bta1 = 0.6, cnit = 1.0e-02 )
      real x,y,tmp,tmp2,del,x2,y2
      real aax,bbx,delrho
      integer :: indxr
      real mwt, nwt, zwt
      real, parameter :: rho00 = 1.225
      integer i
      real cnotmp
      real xvbarmax,xvbar,volfac,xvmxtmp

      real, parameter ::  pi = 3.14159265 ! 4.0*atan(1.0)
      real, parameter ::  piinv = 1.0/pi



!
! set values
!
!      cwmasn = 5.23e-13  ! radius of 5.0e-6
!      cwradn = 5.0e-6
!      cwmasx = 5.25e-10  ! radius of 50.0e-6

      fwlo = 0.2                ! water fraction to start weighting toward rain fall speed
      fwhi = 0.4                ! water fraction at which rain fall speed only is used
      rfwdiff = 1./(fwhi - fwlo)
      
      arx = 10.
      frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.

      ar = 841.99666  
      br = 0.8
      gr = 9.8
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42

!      IF ( ifirst .eq. 0 ) THEN
!        ifirst = 1
        gf4br = gamma(4.0+br)
        gf4ds = gamma(4.0+ds)
!        gf1ds = gamma(1.0+ds)
        gf4p5 = gamma(4.0+0.5)
        gfcinu1 = gamma(cinu + 1.0)
        gfcinu1p47 = gamma(cinu + 1.47167)
        gfcinu2p47 = gamma(cinu + 2.47167)
        gfcinu1p22 = gamma(cinu + 1.22117)
        gfcinu2p22 = gamma(cinu + 2.22117)
        gfcinu1p18 = gamma(cinu + 1.18333)
        gfcinu2p18 = gamma(cinu + 2.18333)
        gsnow1 = gamma(snu + 1.0)
        gsnow53 = gamma(snu + 5./3.)
        gsnow73 = gamma(snu + 7./3.)
        
        IF ( lh  .gt. 1 ) THEN
          IF ( dmuh == 1.0 ) THEN
            cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
          ELSE
            cwchtmp = 6.0*piinv*gamma( (xnu(lh) + 1.)/xmu(lh) )/gamma( (xnu(lh) + 2.)/xmu(lh) )
          ENDIF
        ENDIF
        IF ( lf .gt. 1 ) THEN
          IF ( dmuh == 1.0 ) THEN
            cwcftmp = ((3. + dnu(lf))*(2. + dnu(lf))*(1.0 + dnu(lf)))**(-1./3.)
          ELSE
            cwcftmp = 6.0*piinv*gamma( (xnu(lf) + 1.)/xmu(lf) )/gamma( (xnu(lf) + 2.)/xmu(lf) )
          ENDIF
        ENDIF
        IF ( lhl .gt. 1 ) THEN
          IF ( dmuhl == 1.0 ) THEN
            cwchltmp = ((3. + dnu(lhl))*(2. + dnu(lhl))*(1.0 + dnu(lhl)))**(-1./3.)
          ELSE
            cwchltmp = 6.0*piinv*gamma( (xnu(lhl) + 1)/xmu(lhl) )/gamma( (xnu(lhl) + 2)/xmu(lhl) )
          ENDIF
        ENDIF

        IF ( ipconc .le. 5 .and. imydiagalpha <= 0 ) THEN
          IF ( lh  .gt. 1 ) cwch(:) =  cwchtmp 
          IF ( lf  .gt. 1 ) cwcf(:) =  cwcftmp 
          IF ( lhl .gt. 1 ) cwchl(:) = cwchltmp
        ELSE
          DO mgs = 1,ngscnt
          
          IF ( lh  .gt. 1 .and. ( ildo == 0 .or. ildo == lh )) THEN
           IF ( qx(mgs,lh) .gt. qxmin(lh) ) THEN
            IF ( dmuh == 1.0 ) THEN
              cwch(mgs) = ((3. + alpha(mgs,lh))*(2. + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))**(-1./3.)
             ELSE
             xnutmp = (alpha(mgs,lh) - 2.0)/3.0
             cwch(mgs) =  6.0*piinv*gamma( (xnutmp + 1.)/xmu(lh) )/gamma( (xnutmp + 2.)/xmu(lh) )
            ENDIF
           ELSE
             cwch(mgs) = cwchtmp
           ENDIF
          ENDIF

          IF ( lf  .gt. 1 .and. ( ildo == 0 .or. ildo == lf ) ) THEN
           IF ( qx(mgs,lf) .gt. qxmin(lf) ) THEN
            IF ( dmuh == 1.0 ) THEN
              cwcf(mgs) = ((3. + alpha(mgs,lf))*(2. + alpha(mgs,lf))*(1.0 + alpha(mgs,lf)))**(-1./3.)
             ELSE
             xnutmp = (alpha(mgs,lf) - 2.0)/3.0
             cwcf(mgs) =  6.0*piinv*gamma( (xnutmp + 1.)/xmu(lf) )/gamma( (xnutmp + 2.)/xmu(lf) )
            ENDIF
           ELSE
             cwcf(mgs) = cwcftmp
           ENDIF
          ENDIF

          IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN
           IF ( qx(mgs,lhl) .gt. qxmin(lhl) ) THEN
            IF ( dmuhl == 1.0 ) THEN
              cwchl(mgs) = ((3. + alpha(mgs,lhl))*(2. + alpha(mgs,lhl))*(1.0 + alpha(mgs,lhl)))**(-1./3.)
             ELSE
             xnutmp = (alpha(mgs,lhl) - 2.0)/3.0
             cwchl(mgs) = 6.0*piinv*gamma( (xnutmp + 1)/xmu(lhl) )/gamma( (xnutmp + 2)/xmu(lhl) )
            ENDIF
           ELSE
             cwchl(mgs) = cwchltmp
           ENDIF
          ENDIF
          
          ENDDO
        
        ENDIF
       

      cimasn = Min( cimas0, 6.88e-13)
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter

      cwc1 = 6.0/(pi*1000.)
      cwc0 = piinv ! 6.0*piinv
      mwfac = 6.0**(1./3.)

      
      if (ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set scale diameter, ndebug = ',ndebug1,ngscnt
!


!
!  cloud water variables
! ################################################################
!
!  DROPLETS
!
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cloud water variables'

      IF ( ildo == 0 .or. ildo == lc ) THEN

      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      
      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN !{
      tmp = cx(mgs,lc)
      IF ( ipconc .ge. 2 ) THEN
       IF  ( cx(mgs,lc) .gt. cxmin ) THEN !{
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
       
       ELSE
        ! 02/02/2021: Get here when droplet number has been depleted but still have mass above minimum.
        ! Used to have code here that reset the droplet radius to 5 microns, which can spike
        ! spurious large maximum condensation rates, which
        ! ERM traced back to this size resetting. Instead, now reset ccw < cxmin according
        ! cxmin and cwmasx (maximum droplet mass). Thanks to Xiaofei Li for showing maximum
        ! condensation rates that displayed unexpected trends.
        
        cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
       
       ENDIF
      ELSE
       IF ( ipconc .lt. 2 ) THEN
         cx(mgs,lc) = rho0(mgs)*ccn/rho00 ! scales to local density, relative to standard air density
       ENDIF
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 0.01 ) THEN !{
        xmas(mgs,lc) =  &
     &     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)), &
     &      xdn(mgs,lc)*xvmx(lc) )
        
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)
        
       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 1.0e-9 ) THEN
        cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)

       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 0.01 ) THEN
        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        
       ELSE
        xmas(mgs,lc) = cwmasn
        xv(mgs,lc) = xmas(mgs,lc)/1000.
! do not define ccw here! it can feed back to ccn!!!    cx(mgs,lc) = 0.0 ! cwnc(mgs)
       ENDIF !}
      ENDIF !}
      
      
!      IF ( ipconc .lt. 2 ) THEN
!        xmas(mgs,lc) = &
!     &    min( max(qx(mgs,lc)*rho0(mgs)/cwnc(mgs),cwmasn),cwmasx )
!        cx(mgs,lc) = Max(1.0,qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc))
!      ELSE
!        cwnc(mgs) = an(igs(mgs),jgs,kgs(mgs),lnc)
!        cx(mgs,lc) = cwnc(mgs)
!      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**(1./3.)
      xdia(mgs,lc,2) = xdia(mgs,lc,1)**2
      xdia(mgs,lc,3) = xdia(mgs,lc,1)
      cwrad = 0.5*xdia(mgs,lc,1)
      IF ( fadvisc(mgs) > 0.0 ) THEN
      vtxbar(mgs,lc,1) =  &
     &   (2.0*gr*xdn(mgs,lc) *(cwrad**2)) &
     &  /(9.0*fadvisc(mgs))
      ELSE
       vtxbar(mgs,lc,1) = 0.0
      ENDIF

      
      ELSE
       xmas(mgs,lc) = cwmasn
       xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
       IF ( qx(mgs,lc) <= 0.0 ) cx(mgs,lc) = 0.0
       IF ( ipconc .le. 1 ) cx(mgs,lc) = 0.01
       xdia(mgs,lc,1) = 2.*cwradn
       xdia(mgs,lc,2) = 4.*cwradn**2
       xdia(mgs,lc,3) = xdia(mgs,lc,1)
       vtxbar(mgs,lc,1) = 0.0
       
      ENDIF !} qcw .gt. qxmin(lc)
      
      end do

      ENDIF


!
! cloud ice variables
! columns
!
! ################################################################
!
!  CLOUD ICE
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cip'
      
      IF ( li .gt. 1 .and. ( ildo == 0 .or. ildo == li ) ) THEN
      do mgs = 1,ngscnt
       xdn(mgs,li)  = 900.0
      IF ( ipconc .eq. 0 ) THEN
!       cx(mgs,li) = min(cnit*exp(-temcg(mgs)*bta1),1.e+09)
        cx(mgs,li) = cnina(mgs)
       IF ( cimna .gt. 1.0 ) THEN
         cx(mgs,li) = Max(cimna,cx(mgs,li))
       ENDIF
       IF ( cimxa .gt. 1.0 ) THEN
         cx(mgs,li) = Min(cimxa,cx(mgs,li))
       ENDIF
! erm 3/28/2002
       IF ( itype1a .ge. 1 .or. itype2a .ge. 1 ) THEN
        cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
        cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
       ENDIF
!
       cx(mgs,li) = max(1.0e-20,cx(mgs,li))
!       cx(mgs,li) = Min(ccimx, cx(mgs,li))

      
      ELSEIF ( ipconc .ge. 1 ) THEN
        IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
         cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
         cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
!         cx(mgs,li) = Max(1.0,cx(mgs,li))
        ENDIF
      ENDIF
      
      IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
      xmas(mgs,li) = &
     &     max( qx(mgs,li)*rho0(mgs)/cx(mgs,li), cimasn )
!     &  min( max(qx(mgs,li)*rho0(mgs)/cx(mgs,li),cimasn),cimasx )
      
!      if ( temcg(mgs) .gt. 0.0 ) then
!      xdia(mgs,li,1) = 0.0
!      else
      if ( xmas(mgs,li) .gt. 0.0 ) THEN ! cimasn ) then
!c      xdia(mgs,li,1) = 0.4892*(xmas(mgs,li)**(0.4554))
!       xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))

!       xdia(mgs,li,1) = (132.694*5.40662/xmas(mgs,li))**(-1./2.9163)  ! for inverse exponential distribution
       IF ( ixtaltype == 1 ) THEN ! column
       xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))
       xdia(mgs,li,3) = 0.1871*(xmas(mgs,li)**(0.3429))
       ELSEIF  ( ixtaltype == 2 ) THEN ! disk
        xdia(mgs,li,1) = 0.277823*xmas(mgs,li)**0.359971
        xdia(mgs,li,3) = 0.277823*xmas(mgs,li)**0.359971
       ENDIF
      end if
!      end if
!      xdia(mgs,li,1) = max(xdia(mgs,li,1), 5.e-6)
!      xdia(mgs,li,1) = min(xdia(mgs,li,1), 1000.e-6)

       IF ( ipconc .ge. 0 ) THEN
!      vtxbar(mgs,li,1) = rhovt(mgs)*49420.*40.0005/5.40662*xdia(mgs,li,1)**(1.415) ! mass-weighted
!      vtxbar(mgs,li,1) = (4.942e4)*(xdia(mgs,li,1)**(1.4150))
        xv(mgs,li) = xmas(mgs,li)/xdn(mgs,li)
        IF ( icefallopt == 1 ) THEN ! default ice fall
          IF ( ixtaltype == 1 ) THEN ! column
          tmp = (67056.6300748612*rhovt(mgs))/  &
     &    (((1.0 + cinu)/xv(mgs,li))**0.4716666666666667*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p47
          vtxbar(mgs,li,1) = tmp*gfcinu2p47/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
        ELSEIF  ( ixtaltype == 2 ) THEN ! disk -- but just use Ferrier (1994) snow fall speeds for now
            vtxbar(mgs,li,1) = 11.9495*rhovt(mgs)*(xv(mgs,li))**(0.14)
            vtxbar(mgs,li,2) = 7.02909*rhovt(mgs)*(xv(mgs,li))**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
           vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
        
        ENDIF
          
       ELSEIF ( icefallopt == 2 ) THEN !   ! Ferrier ice fall speed
          tmp = (82.3166*rhovt(mgs))/  &
     &     (((1.0 + cinu)/xv(mgs,li))**0.22117*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p22
          vtxbar(mgs,li,1) = tmp*gfcinu2p22/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 

       ELSEIF ( icefallopt == 3 ) THEN !   ! Adjusted Ferrier (smaller exponent)
       
          tmp = (47.6273*rhovt(mgs))/  &
     &     (((1.0 + cinu)/xv(mgs,li))**0.18333*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p18
          vtxbar(mgs,li,1) = tmp*gfcinu2p18/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
       
       ENDIF
!      vtxbar(mgs,li,1) = vtxbar(mgs,li,2)*(1.+cinu)/(1. + cinu)
!      xdn(mgs,li)   = min(max(769.8*xdia(mgs,li,1)**(-0.0140),300.0),900.0)
!      xdn(mgs,li) = 900.0
        xdia(mgs,li,2) = xdia(mgs,li,1)**2
!      vtxbar(mgs,li,1) = vtxbar(mgs,li,1)*rhovt(mgs)
       ELSE
         xdia(mgs,li,1) = max(xdia(mgs,li,1), 10.e-6)
         xdia(mgs,li,1) = min(xdia(mgs,li,1), 1000.e-6)
         vtxbar(mgs,li,1) = (4.942e4)*(xdia(mgs,li,1)**(1.4150))
!      xdn(mgs,li)   = min(max(769.8*xdia(mgs,li,1)**(-0.0140),300.0),900.0)
         xdn(mgs,li) = 900.0
         xdia(mgs,li,2) = xdia(mgs,li,1)**2
         vtxbar(mgs,li,1) = vtxbar(mgs,li,1)*rhovt(mgs)
         xv(mgs,li) = xmas(mgs,li)/xdn(mgs,li)
       ENDIF ! ipconc gt 3
      ELSE
       xmas(mgs,li) = 1.e-13
       xdn(mgs,li)  = 900.0
       xdia(mgs,li,1) = 1.e-7
       xdia(mgs,li,2) = (1.e-14)
       xdia(mgs,li,3) = 1.e-7
       vtxbar(mgs,li,1) = 0.0
!       cicap(mgs) = 0.0
!       ciat(mgs) = 0.0
      ENDIF
      
      IF ( icefallfac /= 1.0 ) THEN
        vtxbar(mgs,li,1) = icefallfac*vtxbar(mgs,li,1)
        vtxbar(mgs,li,2) = icefallfac*vtxbar(mgs,li,2)
        vtxbar(mgs,li,3) = icefallfac*vtxbar(mgs,li,3)
      ENDIF

      
      
      end do
      
      ENDIF ! li .gt. 1

! ################################################################
! ICE SPHERES

!
      IF ( lis > 1 .and. ( ildo == 0 .or. ildo == lis ) ) THEN
       IF ( ( ndebug1 .gt. 0 ) ) write(0,*) 'SETVTZ: Set ice spheres variables'

      IF ( .false. ) THEN ! treat like cloud droplets
      
      do mgs = 1,ngscnt
      xv(mgs,lis) = 0.0
      
      IF ( qx(mgs,lis) .gt. qxmin(lis) ) THEN !{
      
      IF ( ipconc .ge. 2 .and. (cx(mgs,lis) .gt. 1.0e-9) ) THEN !{
        xmas(mgs,lis) =  &
     &    min( max(qx(mgs,lis)*rho0(mgs)/cx(mgs,lis),cwmasn),cimasx )
        xv(mgs,lis) = xmas(mgs,lis)/xdn(mgs,lis)
      ELSE
       IF ( ipconc .lt. 2 ) THEN
         cx(mgs,lis) = rho0(mgs)*ccn/rho00 ! scales to local density, relative to standard air density
       ENDIF
       IF ( qx(mgs,lis) .gt. qxmin(lis) .and. cx(mgs,lis) .gt. 0.01 ) THEN !{
        xmas(mgs,lis) =  &
     &     min( max(qx(mgs,lis)*rho0(mgs)/cx(mgs,lis),xdn(mgs,lis)*xvmn(lis)), &
     &      xdn(mgs,lis)*xvmx(lis) )
        
        xv(mgs,lis) = xmas(mgs,lis)/xdn(mgs,lis)
        cx(mgs,lis) = qx(mgs,lis)*rho0(mgs)/xmas(mgs,lis)
        
       ELSEIF ( qx(mgs,lis) .gt. qxmin(lis) .and. cx(mgs,lis) .le. 0.01 ) THEN
        xmas(mgs,lis) = xdn(mgs,lis)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lis) = rho0(mgs)*qx(mgs,lis)/xmas(mgs,lis)
        xv(mgs,lis) = xmas(mgs,lis)/xdn(mgs,lis)
        
       ELSE
        xmas(mgs,lis) = cwmasn
        xv(mgs,lis) = xmas(mgs,lis)/900.
! do not define ccw here! it can feed back to ccn!!!    cx(mgs,lis) = 0.0 ! cwnc(mgs)
       ENDIF !}
      ENDIF !}
      xdia(mgs,lis,1) = (xmas(mgs,lis)*cwc1)**(1./3.)
      xdia(mgs,lis,2) = xdia(mgs,lis,1)**2
      xdia(mgs,lis,3) = xdia(mgs,lis,1)
      cwrad = 0.5*xdia(mgs,lis,1)
      IF ( fadvisc(mgs) > 0.0 ) THEN
      vtxbar(mgs,lis,1) =  &
     &   (2.0*gr*xdn(mgs,lis) *(cwrad**2)) &
     &  /(9.0*fadvisc(mgs))
      ELSE
       vtxbar(mgs,lis,1) = 0.0
      ENDIF

      
      ELSE
       xmas(mgs,lis) = cwmasn
       IF ( ipconc .le. 1 ) cx(mgs,lis) = 0.01
       xdia(mgs,lis,1) = 2.*cwradn
       xdia(mgs,lis,2) = 4.*cwradn**2
       xdia(mgs,lis,3) = xdia(mgs,lis,1)
       vtxbar(mgs,lis,1) = 0.0
       
      ENDIF !} qis .gt. qxmin(lis)
      
      end do
      
      
      ELSE ! treat as columns
      
      do mgs = 1,ngscnt
       xdn(mgs,lis)  = 900.0

        IF ( qx(mgs,lis) .gt. qxmin(lis) ) THEN
         cx(mgs,lis) = Max(cx(mgs,lis),qx(mgs,lis)*rho0(mgs)/cimasx)
!         cx(mgs,lis) = Min(cx(mgs,lis),qx(mgs,lis)*rho0(mgs)/cwmasn) ! this is cwmasn to be frozen droplets
         cx(mgs,lis) = Min(cx(mgs,lis),qx(mgs,lis)*rho0(mgs)/cimasn) ! this is cwmasn to be frozen droplets
        ENDIF
      
      IF ( qx(mgs,lis) .gt. qxmin(lis) ) THEN
      xmas(mgs,lis) = &
     &     max( qx(mgs,lis)*rho0(mgs)/cx(mgs,lis), cimasn )
!     &     max( qx(mgs,lis)*rho0(mgs)/cx(mgs,lis), cwmasn )
      if ( xmas(mgs,lis) .gt. 0.0 ) THEN
       IF ( ixtaltype == 1 ) THEN ! column
       xdia(mgs,lis,1) = 0.1871*(xmas(mgs,lis)**(0.3429))
       xdia(mgs,lis,3) = 0.1871*(xmas(mgs,lis)**(0.3429))
       ELSEIF  ( ixtaltype == 2 ) THEN ! disk
        xdia(mgs,lis,1) = 0.277823*xmas(mgs,lis)**0.359971
        xdia(mgs,lis,3) = 0.277823*xmas(mgs,lis)**0.359971
       ENDIF
      end if


       IF ( ipconc .ge. 0 ) THEN
        xv(mgs,lis) = xmas(mgs,lis)/xdn(mgs,lis)
        IF ( ixtaltype == 1 ) THEN ! column
        tmp = (67056.6300748612*rhovt(mgs))/  &
     &   (((1.0 + cinu)/xv(mgs,lis))**0.4716666666666667*gfcinu1)
        vtxbar(mgs,lis,2) = tmp*gfcinu1p47
        vtxbar(mgs,lis,1) = tmp*gfcinu2p47/(1. + cinu)
        vtxbar(mgs,lis,3) = vtxbar(mgs,lis,1) 
        ELSEIF  ( ixtaltype == 2 ) THEN ! disk -- but just use column fall speed for now
        tmp = (67056.6300748612*rhovt(mgs))/  &
     &   (((1.0 + cinu)/xv(mgs,lis))**0.4716666666666667*gfcinu1)
        vtxbar(mgs,lis,2) = tmp*gfcinu1p47
        vtxbar(mgs,lis,1) = tmp*gfcinu2p47/(1. + cinu)


        vtxbar(mgs,lis,3) = vtxbar(mgs,lis,1) 
        
        ENDIF
        xdia(mgs,lis,2) = xdia(mgs,lis,1)**2
       ELSE
!         xdia(mgs,lis,1) = max(xdia(mgs,lis,1), 5.e-6)
         xdia(mgs,lis,1) = min(xdia(mgs,lis,1), 1000.e-6)
         vtxbar(mgs,lis,1) = (4.942e4)*(xdia(mgs,lis,1)**(1.4150))
         xdn(mgs,lis) = 900.0
         xdia(mgs,lis,2) = xdia(mgs,lis,1)**2
         vtxbar(mgs,lis,1) = vtxbar(mgs,lis,1)*rhovt(mgs)
         xv(mgs,lis) = xmas(mgs,lis)/xdn(mgs,lis)
       ENDIF ! ipconc gt 3
      ELSE
       xmas(mgs,lis) = 1.e-13
       xdn(mgs,lis)  = 900.0
       xdia(mgs,lis,1) = 1.e-7
       xdia(mgs,lis,2) = (1.e-14)
       xdia(mgs,lis,3) = 1.e-7
       vtxbar(mgs,lis,1) = 0.0
      ENDIF
      end do
      
      ENDIF
      
      ENDIF ! lis > 1

! ################################################################
!
!  RAIN
!
      
!
      IF ( ildo == 0 .or. ildo == lr ) THEN

      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      
!      IF ( qx(mgs,lr) .gt. 10.0e-3 ) &
!     &  print*, 'RAIN1: ',igs(mgs),kgs(mgs),qx(mgs,lr)
      
      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        xvbarmax = xvmx(lr)
        IF ( imaxdiaopt == 1 ) THEN
          xvbarmax = xvmx(lr)
        ELSEIF ( imaxdiaopt == 2 ) THEN ! test against maximum mass diameter
         IF ( imurain == 1 ) THEN
         ! maximum mass diameter is xdia(mgs,lr,3)*(3.0 + alpha(mgs,lr))*(gamma(1.+alpha(mgs,lr))/Gamma(4.0 + alpha(mgs,lr)))**(1./3.)
         ! expand Gamma(4+alp) = (3+alp)*Gamma(3+alp) = (3+alp)*(2+alp)*Gamma(2+alp) = (3+alp)*(2+alp)*(1+alp)*Gamma(1+alp)
         ! to get xdiamaxmass =  xdia(mgs,lr,3)*(3.0 + alpha(mgs,lr))*((3.+alpha(mgs,lr))*(2.+alpha(mgs,lr))*(1. + alpha(mgs,lr)) )**(-1./3.)
         ! cube that to do volume ratio
           xvbarmax = xvmx(lr)/((3. + alpha(mgs,lr))**3/((3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr))))
         ELSEIF ( imurain == 3 ) THEN
            ! no code for this yet
         ENDIF
        ELSEIF ( imaxdiaopt == 3 ) THEN ! test against mass-weighted diameter
         IF ( imurain == 1 ) THEN
           xvbarmax = xvmx(lr)/((4. + alpha(mgs,lr))**3/((3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr))))
         ELSEIF ( imurain == 3 ) THEN
            ! no code for this yet
         ENDIF
        ENDIF
        
        ! write(0,*) 'setvt: k,xvbarmax,alpr = ',kgs(mgs),xvbarmax,xv(mgs,lr)
        ! xv(mgs,lr) = 1.01*xvbarmax
       
        IF ( xv(mgs,lr) .gt. xvbarmax ) THEN
          xv(mgs,lr) = xvbarmax
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvbarmax*xdn(mgs,lr))
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF


        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.) ! xdia(mgs,lr,1)
        IF ( imurain == 3 ) THEN
!          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/(alpha(mgs,lr)+1.))**(1./3.)
          xdia(mgs,lr,1) = xdia(mgs,lr,3) ! formulae for Ziegler (1985) use mean volume diameter, not lambda**(-1)
        ELSE ! imurain == 1, Characteristic diameter (1/lambda)
          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
        ENDIF
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

! Inverse exponential version:
!        xdia(mgs,lr,1) =
!     &  (qx(mgs,lr)*rho0(mgs)
!     & /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      ELSE
        xdia(mgs,lr,1) = &
     &  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.)
        cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
        xdia(mgs,lr,3) = 1.e-9
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if
      xdia(mgs,lr,2) = xdia(mgs,lr,1)**2
!      xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
      end do

      ENDIF

! ################################################################
!
!  SNOW
!

      IF ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) THEN
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      if ( ipconc .ge. 4 ) then ! two-moment, gamma of volume

        xmas(mgs,ls) =  rho0(mgs)*qx(mgs,ls)/(Max(1.0e-9,cx(mgs,ls)))
        swmasmx = 13.7e-6
!       IF ( xmas(mgs,ls) > swmasmx ) THEN
!          xmas(mgs,ls) = swmasmx
!          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
!        ENDIF
          
        
        IF ( isnowdens == 2 ) THEN ! Set values according to Cox relationship
        
         ! xdn(mgs,ls) = 0.0346159*Sqrt(cx(mgs,ls)/(qx(mgs,ls)*rho0(mgs)) )
          xdn(mgs,ls) = 0.0346159/Sqrt(Min( swmasmx, xmas(mgs,ls) ) )
        !  xdn(mgs,ls) = Max( 100.0, xdn(mgs,ls) )  ! limit snow to 100. to keep other equations in line
          
          IF ( xdn(mgs,ls) <= 900. ) THEN
             dtmp = Sqrt( xmas(mgs,ls)/0.069 ) ! diameter (meters) of mean mass particle using Cox 1998 relation (m = p d^2)
             xv(mgs,ls) = xmas(mgs,ls)/xdn(mgs,ls) ! 28.8887*xmas(mgs,ls)**(3./2.)
          ELSE ! at small sizes, assume ice spheres
             xdn(mgs,ls) = 900.
             xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
             dtmp = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
          ENDIF
          
        ELSE ! leave xdn(ls) at default value
             xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
             dtmp = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
        ENDIF

!        xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
!      parameter( xvmn(lr)=2.8866e-13, xvmx(lr)=4.1887e-9 )  ! mks
!        xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)

        xdia(mgs,ls,1) = dtmp ! (xv(mgs,ls)*cwc0*6.0)**(1./3.)

        IF ( xv(mgs,ls) .lt. xvmn(ls) .and. isnowdens == 1) THEN
          xv(mgs,ls) = Max( xvmn(ls),xv(mgs,ls) )
          xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
        ENDIF

        IF ( xv(mgs,ls) .gt. xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls))) ) THEN
          xv(mgs,ls) = Min( xvmx(ls), Max( xvmn(ls),xv(mgs,ls) ) )
          xmas(mgs,ls) = 0.106214*xv(mgs,ls)**(2./3.)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdn(mgs,ls) = 0.0346159*Sqrt(cx(mgs,ls)/(qx(mgs,ls)*rho0(mgs)) )
          xdia(mgs,ls,1) = Sqrt( xmas(mgs,ls)/0.069 ) 
        ENDIF

#ifdef USEMIXEDPHASE
        IF ( mixedphase ) THEN
            
            fsw = Min(1.0, qxw(mgs,ls)/qx(mgs,ls) )

            IF ( temcg(mgs) .gt. 0.0 .and. fsw .ge. 0.10) then
             frac = (fsw - 0.1)/0.9 ! scale range from 0.1-1 to 0-1.
             rwdia = (xmas(mgs,ls)*cwc1)**(1./3.)                       !melted diameter
             xdia(mgs,ls,1) = (1.0 - frac)*xdia(mgs,ls,1) + frac*rwdia
            ENDIF
            
         ENDIF

#endif
        xdia(mgs,ls,3) = xdia(mgs,ls,1)

      ELSE ! single moment, inverse exponential of diameter
        xdia(mgs,ls,1) =  &
     &    (qx(mgs,ls)*rho0(mgs)/(pi*xdn(mgs,ls)*cno(ls)))**(0.25) 
        cx(mgs,ls) = cno(ls)*xdia(mgs,ls,1)
        xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
        xdia(mgs,ls,3) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
      end if
      else
      xdia(mgs,ls,1) = 1.e-9
      xdia(mgs,ls,3) = 1.e-9
      cx(mgs,ls) = 0.0
      
       IF ( isnowdens == 2 ) THEN ! Set values according to Cox relationship
         xdn(mgs,ls) = 90.
       ENDIF
      
      end if
      xdia(mgs,ls,2) = xdia(mgs,ls,1)**2
!      swdia3(mgs) = xdia(mgs,ls,2)*xdia(mgs,ls,1)
!      xmas(mgs,ls) = xdn(mgs,ls)*(pi/6.)*swdia3(mgs)
      end do
      
      ENDIF ! ls .gt 1
!
!
! ################################################################
!
!  GRAUPEL
!

      IF ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) THEN
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*Max(1.0e-9,cx(mgs,lh)))
        xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)

        tmp = 0.
        IF ( xv(mgs,lh) .lt. xvmn(lh) .or. xv(mgs,lh) .gt. xvmx(lh) ) THEN
          tmp = xv(mgs,lh)
          xv(mgs,lh) = Min( xvmx(lh), Max( xvmn(lh),xv(mgs,lh) ) )
          xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
          cx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xmas(mgs,lh))
        ENDIF

         IF ( ipconc == 5 .and. imorrgdnglimit == 1 ) THEN
           ! limit on characteristic diameter (i.e., 1/slope)
           xvmxtmp = pi/6.0*(morrdnglimit/cwch(mgs))**3
            IF ( xv(mgs,lh) .gt. xvmxtmp ) THEN
               xdia(mgs,lh,1) = morrdnglimit
               tmp = xv(mgs,lh)
               xdia(mgs,lh,3) = xdia(mgs,lh,1)/cwch(mgs)
               xv(mgs,lh) = Min( xvmxtmp, Max( xvmn(lh),xv(mgs,lh) ) )
               xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
               cx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xmas(mgs,lh))
            ENDIF
         
         ENDIF

         xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
         IF ( dmuh == 1.0 ) THEN
           xdia(mgs,lh,1) = cwch(mgs)*xdia(mgs,lh,3)
         ELSE
           xdia(mgs,lh,1) = (xv(mgs,lh)*cwch(mgs))**(1./3.)
         ENDIF

      ELSE
      IF ( ithompsoncnoh == 0 ) THEN
       cnotmp = cno(lh)
      ELSE
       cnotmp = Max( 1.e4, Min( 200./qx(mgs,lh), 5.e6 ) )
      ENDIF
      xdia(mgs,lh,1) =  &
     &  (qx(mgs,lh)*rho0(mgs)/(pi*xdn(mgs,lh)*cnotmp))**(0.25) 
      cx(mgs,lh) = cnotmp*xdia(mgs,lh,1)
      xv(mgs,lh) = Max(xvmn(lh), rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh)) )
      xdia(mgs,lh,3) = (xv(mgs,lh)*6./pi)**(1./3.) 
      end if
      else
      xdia(mgs,lh,1) = 1.e-9
      xdia(mgs,lh,3) = 1.e-9
      end if
      xdia(mgs,lh,2) = xdia(mgs,lh,1)**2
!      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
!      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
      
      ENDIF


!
!
! ################################################################
!
!  Frozen Drops
!

      IF ( lf .gt. 1 .and. ( ildo == 0 .or. ildo == lf ) ) THEN
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lf) = rho0(mgs)*qx(mgs,lf)/(xdn(mgs,lf)*Max(1.0e-9,cx(mgs,lf)))
        xmas(mgs,lf) = xv(mgs,lf)*xdn(mgs,lf)

        IF ( xv(mgs,lf) .lt. xvmn(lf) .or. xv(mgs,lf) .gt. xvmx(lf) ) THEN
          xv(mgs,lf) = Min( xvmx(lf), Max( xvmn(lf),xv(mgs,lf) ) )
          xmas(mgs,lf) = xv(mgs,lf)*xdn(mgs,lf)
          cx(mgs,lf) = rho0(mgs)*qx(mgs,lf)/(xmas(mgs,lf))
        ENDIF

         xdia(mgs,lf,3) = (xv(mgs,lf)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lf,1) ! (xv(mgs,lf)*cwc0*6.0)**(1./3.)
         IF ( dmuh == 1.0 ) THEN
           xdia(mgs,lf,1) = cwcf(mgs)*xdia(mgs,lf,3)
         ELSE
           xdia(mgs,lf,1) = (xv(mgs,lf)*cwcf(mgs))**(1./3.)
         ENDIF

      ELSE
      IF ( ithompsoncnoh == 0 ) THEN
       cnotmp = cno(lf)
      ELSE
       cnotmp = Max( 1.e4, Min( 200./qx(mgs,lf), 5.e6 ) )
      ENDIF
      xdia(mgs,lf,1) =  &
     &  (qx(mgs,lf)*rho0(mgs)/(pi*xdn(mgs,lf)*cnotmp))**(0.25) 
      cx(mgs,lf) = cnotmp*xdia(mgs,lf,1)
      xv(mgs,lf) = Max(xvmn(lf), rho0(mgs)*qx(mgs,lf)/(xdn(mgs,lf)*cx(mgs,lf)) )
      xdia(mgs,lf,3) = (xv(mgs,lf)*6./pi)**(1./3.) 
      end if
      else
      xdia(mgs,lf,1) = 1.e-9
      xdia(mgs,lf,3) = 1.e-9
      end if
      xdia(mgs,lf,2) = xdia(mgs,lf,1)**2
!      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
!      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
      
      ENDIF


!
! ################################################################
!
!  HAIL
!

      IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*Max(1.0e-9,cx(mgs,lhl)))
        xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xmas(mgs,lhl),qx(mgs,lhl)

        IF ( xv(mgs,lhl) .lt. xvmn(lhl) .or. xv(mgs,lhl) .gt. xvmx(lhl) ) THEN
          xv(mgs,lhl) = Min( xvmx(lhl), Max( xvmn(lhl),xv(mgs,lhl) ) )
          xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
          cx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xmas(mgs,lhl))
        ENDIF

        xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
         IF ( dmuhl == 1.0 ) THEN
           xdia(mgs,lhl,1) = cwchl(mgs)*xdia(mgs,lhl,3)
         ELSE
           xdia(mgs,lhl,1) = (xv(mgs,lhl)*cwchl(mgs))**(1./3.)
         ENDIF
        
!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xdia(mgs,lhl,3)
      ELSE
      xdia(mgs,lhl,1) = &
     &  (qx(mgs,lhl)*rho0(mgs)/(pi*xdn(mgs,lhl)*cno(lhl)))**(0.25) 
      cx(mgs,lhl) = cno(lhl)*xdia(mgs,lhl,1)
      xv(mgs,lhl) = Max(xvmn(lhl), rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl)) )
      xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.) 
      end if
      else
      xdia(mgs,lhl,1) = 1.e-9
      xdia(mgs,lhl,3) = 1.e-9
      end if
      xdia(mgs,lhl,2) = xdia(mgs,lhl,1)**2
!      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
!      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
      
      ENDIF
!      
!
!
!  Set terminal velocities...
!    also set drag coefficients (moved to start of subroutine)
!
!      cdx(lr) = 0.60
!      cdx(lh) = 0.45
!      cdx(lhl) = 0.45
!      cdx(lf) = 0.45
!      cdx(lgh) = 0.60
!      cdx(lgm) = 0.80
!      cdx(lgl) = 0.80
!      cdx(lir) = 2.00
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set terminal velocities'
!
!
! ################################################################
!
!  RAIN FALL SPEED
!
      IF ( ildo == 0 .or. ildo == lr ) THEN

      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      IF ( ipconc .lt. 3 ) THEN
        vtxbar(mgs,lr,1) = rainfallfac*(ar*gf4br/6.0)*(xdia(mgs,lr,1)**br)*rhovt(mgs)
!        write(91,*) 'vtxbar: ',vtxbar(mgs,lr,1),mgs,gf4br,xdia(mgs,lr,1),rhovt(mgs)
      ELSE
        
        IF ( imurain == 1 ) THEN ! DSD of Diameter
        
        ! using functional form of  arx*(1 - Exp(-frx*diameter) ), with arx = 10.
        !  and frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.
        ! Similar form as in Atlas et al. (1973), who had 9.65 - 10.3*Exp[-600 * d]

        
          alp = alpha(mgs,lr)
          
          vtxbar(mgs,lr,1) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 4.0) ) ! mass weighted
          
          IF ( infdo .ge. 1 .and. rssflg == 1 ) THEN
            vtxbar(mgs,lr,2) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alpha(mgs,lr) - 1.0) ) ! number weighted
          ELSEIF ( infdo .ge. 1 .and. rssflg == 2 ) THEN ! test code to try to reduce size sorting at the melting layer
!            IF ( temcg(mgs) > 0. .and. temcg(mgs) < 10. .and. qx(mgs,ls) > 1.e-4 .and. qx(mgs,ls) >= qx(mgs,lr) ) THEN
            IF ( temcg(mgs) > 0. .and. temcg(mgs) < 10. .and. qx(mgs,ls) > qxmin(ls) ) THEN
              vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
            ELSE
              vtxbar(mgs,lr,2) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 1.0) ) ! number weighted
            ENDIF
          ELSE
            vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
          ENDIF
          
          IF ( infdo .ge. 2 .and. rssflg == 1 ) THEN
            vtxbar(mgs,lr,3) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 7.0) ) ! z-weighted
          ELSEIF ( infdo .ge. 2 .and. rssflg == 2 ) THEN
            IF ( temcg(mgs) > 0. .and. temcg(mgs) < 10. .and. qx(mgs,ls) > 1.e-4 .and. qx(mgs,ls) >= qx(mgs,lr) ) THEN
              vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
            ELSE
              vtxbar(mgs,lr,3) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 7.0) ) ! z-weighted
            ENDIF
          ELSE
            vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
          ENDIF
          
!          write(91,*) 'setvt: alp,vn,vm,vz = ',alp,vtxbar(mgs,lr,2), vtxbar(mgs,lr,1), vtxbar(mgs,lr,3),alphan(mgs,lr)

      IF ( .false. .and. infdo >= 1 .and. imydiagalpha == 3 ) THEN
       ! diagnostic shape parameter for Vn (Milbrandt and M-C 2010)
       ! Mk = n0*gamma(k+alpha+1)*xdia(mgs,lh,1)**(k+alpha+1)
       ! M3/M0 = gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3)
       !  tmp = (gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3))**(1./3.)

        tmp = 4. + alpha(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alpha(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = (x/y)**(1./3.)*xdia(mgs,lr,1)

        alphan(mgs,lr) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)

        tmp = 4. + alphan(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alphan(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alphan(mgs,lr) + bxx(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alphan(mgs,lr) + bxx(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
        ! vn = vm*gamma(0+1+alpha+b)*gamma(3+1+alpha)/(gamma(3+1+alpha+b)*gamma(0+1+alpha))
        
        vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)*y2*x/(x2*y)
        
        ! write(0,*) 'alphan: ',alphan(mgs,lh),qx(mgs,lh)*1.e3,1.e3*xdia(mgs,lh,1),x,y,(x/y)**(1./3.)

      ELSEIF ( infdo >= 1 .and. imydiagalpha == 4 ) THEN
        ! Thompson hack
   !     alphan(mgs,lr) = alpha(mgs,lr) + 1.5
   !     vtxbar(mgs,lr,2) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alpha(mgs,lr) - 1.0) ) ! number weighted

      ENDIF

        ELSEIF ( imurain == 3 ) THEN ! DSD of Volume
        
        IF ( lzr < 1 ) THEN ! not 3-moment rain
        rwdia = Min( xdia(mgs,lr,1), 8.0e-3 )
        
         vtxbar(mgs,lr,1) = rhovt(mgs)*6.0*piinv*( 0.04771 + 3788.0*rwdia -  &
     &        1.105e6*rwdia**2 + 1.412e8*rwdia**3 - 6.527e9*rwdia**4)
        
        IF ( infdo .ge. 1 ) THEN
          IF (  rssflg >= 1 ) THEN
         vtxbar(mgs,lr,2) = (0.09112 + 2714.0*rwdia - 4.872e5*rwdia**2 +  &
     &            4.495e7*rwdia**3 - 1.626e9*rwdia**4)*rhovt(mgs)
          ELSE
            vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
          ENDIF
        ENDIF
        
        IF ( infdo .ge. 2 ) THEN ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  = rhovt(mgs)*(  &
     &       0.0911229 +                  &
     &  9246.494*(rwdia) -               &
     &  3.2839926e6*(rwdia**2) +          &
     &  4.944093e8*(rwdia**3) -          &
     &  2.631718e10*(rwdia**4) )
        ENDIF
        
        ELSE ! 3-moment rain, gamma-volume

        vr = xv(mgs,lr)
        rnux = alpha(mgs,lr)
        
        IF ( infdo .ge. 1 .and. rssflg == 1) THEN ! number-weighted; DTD: added size-sorting flag
        vtxbar(mgs,lr,2) = rhovt(mgs)*                             &
     &     (((1. + rnux)/vr)**(-1.333333)*                         &
     &    (0.0911229*((1. + rnux)/vr)**1.333333*Gamma(1. + rnux) + &
     &      (5430.3131*(1. + rnux)*Gamma(4./3. + rnux))/           &
     &       vr - 1.0732802e6*((1. + rnux)/vr)**0.6666667*         &
     &       Gamma(1.666667 + rnux) +                              &
     &      8.584110982429507e7*((1. + rnux)/vr)**(1./3.)*         &
     &       Gamma(2. + rnux) -                                    &
     &      2.3303765697228556e9*Gamma(7./3. + rnux)))/            &
     &  Gamma(1. + rnux)
        ENDIF

!  mass-weighted
       vtxbar(mgs,lr,1)  = rhovt(mgs)*                                                 &
     &   (0.0911229*(1 + rnux)**1.3333333333333333*Gamma(2. + rnux) +                  &
     &    5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                         &
     &     Gamma(2.333333333333333 + rnux) -                                           &
     &    1.0732802065650471e6*(1 + rnux)**0.6666666666666666*vr**0.6666666666666666*  &
     &     Gamma(2.6666666666666667 + rnux) +                                          &
     &    8.584110982429507e7*(1 + rnux)**0.3333333333333333*vr*Gamma(3 + rnux) -      &
     &    2.3303765697228556e9*vr**1.3333333333333333*                                 &
     &     Gamma(3.333333333333333 + rnux))/                                           &
     &  ((1 + rnux)**2.333333333333333*Gamma(1 + rnux)) 
     
        IF(infdo .ge. 1 .and. rssflg == 0) THEN ! No size-sorting, set N-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
        ENDIF     
      
        IF ( infdo .ge. 2 .and. rssflg == 1) THEN ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  =   rhovt(mgs)*                                          &
     &  ((1. + rnux)*(0.0911229*(1 + rnux)**1.3333333333333333*Gamma(3. + rnux) +  &
     &      5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                   &
     &       Gamma(3.3333333333333335 + rnux) -                                    &
     &      1.0732802065650471e6*(1 + rnux)**0.6666666666666666*                   &
     &       vr**0.6666666666666666*Gamma(3.6666666666666665 + rnux) +             &
     &      8.5841109824295e7*(1 + rnux)**0.3333333333333333*vr*Gamma(4. + rnux) - &
     &      2.3303765697228556e9*vr**1.3333333333333333*                           &
     &       Gamma(4.333333333333333 + rnux)))/                                    &
     &  ((1 + rnux)**3.3333333333333335*(2 + rnux)*Gamma(1 + rnux))
        
!         write(0,*) 'setvt: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
        
        ELSEIF (infdo .ge. 2) THEN ! No size-sorting, set Z-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
        ENDIF
        
        
        ENDIF
       ENDIF ! imurain

!        IF ( rwrad*mwfac .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,1) = 20.1*Sqrt(100.*rwrad*mwfac)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,1) = 80.0e2*rwrad*rhovt(mgs)*mwfac
!        ENDIF
!        IF ( rwrad .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,2) = 20.1*Sqrt(100.*rwrad)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,2) = 80.0e2*rwrad*rhovt(mgs)
!        ENDIF
      ENDIF ! ipconc
      else  ! qr < qrmin
      vtxbar(mgs,lr,1) = 0.0
      vtxbar(mgs,lr,2) = 0.0
      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set rain vt'

      ENDIF
!
! ################################################################
!
!  SNOW  FALL SPEED ! Zrnic et al. (1993) or Ferrier (1994, JAS)
!
      IF ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) THEN

      do mgs = 1,ngscnt


      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
        IF ( ipconc .ge. 4 ) THEN
         if ( mixedphase .and. qsvtmod ) then
#ifdef USEMIXEDPHASE
          IF ( qx(mgs,ls) .gt. qxmin(ls) ) fsw = Min(1.0, qxw(mgs,ls)/qx(mgs,ls) )
          IF ( isnowfall == 1 ) THEN
           ! original (Zrnic et al. 1993)
           vtxbar(mgs,ls,1) = 5.72462*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
          ELSEIF ( isnowfall == 2 ) THEN
          ! Ferrier:
            vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls))**(0.14)
          ELSEIF ( isnowfall == 3 ) THEN
          ! Cox, mass distrib:
            vtxbar(mgs,ls,1) = 50.092*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
          ENDIF
          IF(Abs(sssflg) >= 1) THEN
            IF ( isnowfall == 1 ) THEN
              vtxbar(mgs,ls,2) = 4.04091*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
            ELSEIF ( isnowfall == 2 ) THEN
            ! Ferrier:
              vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls))**(0.14) ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
            ELSEIF ( isnowfall == 3 ) THEN
            ! Cox, mass distrib:
              vtxbar(mgs,ls,2) = 21.6147*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            ENDIF
          ELSE
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
          ENDIF
          if ( temcg(mgs) .gt. 0.0 .and. fsw .gt. 0.10) then
!           vtxbar(mgs,ls,1) = (1.-fsw)*vtxbar(mgs,ls,1) + fsw*vtxbar(mgs,lr,1)
           rwdia = (xmas(mgs,ls)*cwc1)**(1./3.)                       !melted diameter

           IF ( infdo  >= 2 ) THEN
            IF ( isnowfall == 1 ) THEN
             vtxbar(mgs,ls,3) = 6.12217*rhovt(mgs)*(xv(mgs,ls))**(1./12.) ! Zrnic et al 93
            ELSEIF ( isnowfall == 2 ) THEN
             vtxbar(mgs,ls,3) = 13.3436*rhovt(mgs)*(xv(mgs,ls))**(0.14)   ! Ferrier 94
            ELSEIF ( isnowfall == 3 ) THEN
            ! Cox, mass distrib:
              vtxbar(mgs,ls,3) = 61.0914*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            ENDIF
           ENDIF
          
           IF ( imusnow == 3 ) THEN ! gamma volume
           
           ! mass weighted fall speed
           mwt = rhovt(mgs)*6.0*piinv*( 0.04771 + 3788.0*rwdia -   &
     &        1.105e6*rwdia**2 + 1.412e8*rwdia**3 - 6.527e9*rwdia**4) !fall speed for melted snow
           tmp = (1.-fsw)*vtxbar(mgs,ls,1) + fsw*mwt        !linearly wtd by liq water frac
           vtxbar(mgs,ls,1) = max(tmp,vtxbar(mgs,ls,1))               !liq mass should incr fall speed
          
          ! number weighted fall speed (added 08/20/2015)
            nwt = (0.09112 + 2714.0*rwdia - 4.872e5*rwdia**2 +  &
     &            4.495e7*rwdia**3 - 1.626e9*rwdia**4)*rhovt(mgs)
           tmp = (1.-fsw)*vtxbar(mgs,ls,2) + fsw*nwt        !linearly wtd by liq water frac
           vtxbar(mgs,ls,2) = max(tmp,vtxbar(mgs,ls,2))               !liq mass should incr fall speed

           IF ( infdo  >= 2 ) THEN
             zwt  = rhovt(mgs)*(  &
     &              0.0911229 +                  &
     &              9246.494*(rwdia) -               &
     &              3.2839926e6*(rwdia**2) +          &
     &              4.944093e8*(rwdia**3) -          &
     &              2.631718e10*(rwdia**4) )
           tmp = (1.-fsw)*vtxbar(mgs,ls,3) + fsw*zwt        !linearly wtd by liq water frac
           vtxbar(mgs,ls,3) = max(tmp,vtxbar(mgs,ls,3))               !liq mass should incr fall speed
            
           ENDIF

           ELSEIF ( imusnow == 1 ) THEN ! gamma diameter
           
             write(0,*) 'SETVT: no code yet for imusnow == 1 for mixed phase snow with qsvtmod=true.'
             call commasmpi_abort()
           
           ENDIF

         IF ( sssflg < 0 .and. temcg(mgs) > Abs(sssflg) ) THEN ! above a given temperature, effectively turn off size sorting
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
            vtxbar(mgs,ls,3) = vtxbar(mgs,ls,1)
         ENDIF

          
          endif
#endif
         else
          IF ( isnowfall == 1 ) THEN
           ! original (Zrnic et al. 1993)
           vtxbar(mgs,ls,1) = 5.72462*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
          ELSEIF ( isnowfall == 2 ) THEN
          ! Ferrier:
            IF ( isnowdens == 1 ) THEN
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls))**(0.14)
            ELSE
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14) 
            ENDIF
          ELSEIF ( isnowfall == 3 ) THEN
          ! Cox, mass distrib:
            vtxbar(mgs,ls,1) = 50.092*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
          ENDIF
          
          IF(Abs(sssflg) >= 1) THEN
            IF ( isnowfall == 1 ) THEN
              vtxbar(mgs,ls,2) = 4.04091*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
            ELSEIF ( isnowfall == 2 ) THEN
            ! Ferrier:
              IF ( isnowdens == 1 ) THEN
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls))**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              ELSE
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              ENDIF
            ELSEIF ( isnowfall == 3 ) THEN
            ! Cox, mass distrib:
              vtxbar(mgs,ls,2) = 21.6147*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            ENDIF
          ELSE
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
          ENDIF
           IF ( infdo  >= 2 ) THEN
            IF ( isnowfall == 1 ) THEN
             vtxbar(mgs,ls,3) = 6.12217*rhovt(mgs)*(xv(mgs,ls))**(1./12.) ! Zrnic et al 93
            ELSEIF ( isnowfall == 2 ) THEN
             vtxbar(mgs,ls,3) = 13.3436*rhovt(mgs)*(xv(mgs,ls))**(0.14)   ! Ferrier 94
            ELSEIF ( isnowfall == 3 ) THEN
            ! Cox, mass distrib:
              vtxbar(mgs,ls,3) = 61.0914*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            ENDIF
           ENDIF
         
         IF ( sssflg < 0 .and. temcg(mgs) > Abs(sssflg) ) THEN ! above a given temperature, effectively turn off size sorting
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
            vtxbar(mgs,ls,3) = vtxbar(mgs,ls,1)
         ENDIF
         
         endif
        ELSE ! single-moment:
         vtxbar(mgs,ls,1) = (cs*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        ENDIF
      else
      vtxbar(mgs,ls,1) = 0.0
      end if

      IF ( snowfallfac /= 1.0 ) THEN
        vtxbar(mgs,ls,1) = snowfallfac*vtxbar(mgs,ls,1)
        vtxbar(mgs,ls,2) = snowfallfac*vtxbar(mgs,ls,2)
        vtxbar(mgs,ls,3) = snowfallfac*vtxbar(mgs,ls,3)
      ENDIF

        IF ( .false. .and. .not. mixedphase ) THEN
        ! hack to slow down rain to snow speed in melting layer
        IF ( qx(mgs,ls) .gt. qxmin(ls) .and. qx(mgs,lr) .gt. qxmin(lr) &
            .and. temcg(mgs) .gt. 0.0 ) THEN
          tmp = Min( qx(mgs,ls), qx(mgs,lr) )/qx(mgs,lr)
          vtxbar(mgs,lr,1) = (tmp*vtxbar(mgs,ls,1) + (1.0 - tmp)*vtxbar(mgs,lr,1))
          vtxbar(mgs,lr,2) = (tmp*vtxbar(mgs,ls,2) + (1.0 - tmp)*vtxbar(mgs,lr,2))
          vtxbar(mgs,lr,3) = (tmp*vtxbar(mgs,ls,3) + (1.0 - tmp)*vtxbar(mgs,lr,3))
          
        ENDIF
        ENDIF


      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set snow vt'
      
      ENDIF ! ls .gt. 1

!
!
! ################################################################
!
!  Frozen Drops  FALL SPEED ! Wisner et al. (1972)
!
      IF ( lf .gt. 1 .and. ( ildo == 0 .or. ildo == lf ) ) THEN
      
      do mgs = 1,ngscnt
      vtxbar(mgs,lf,1) = 0.0
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
       IF ( icdx .eq. 1 ) THEN
         cd = cdx(lf)
       ELSEIF ( icdx .eq. 2 ) THEN
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lf) - xdn(mgs,lf))/(xdnmx(lf)-xdnmn(lf)) ) )
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lf))/(900. - 300.) ) )
         cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 500.) ) )
!         cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 500.) ) )
       ELSEIF ( icdx .eq. 3 ) THEN
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 170.0) ) )
       ELSEIF ( icdx .eq. 4 ) THEN
         cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
     &        (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lf) ) ) )/(cdhdnmax - cdhdnmin) ) )
       ELSEIF ( icdx .eq. 5 ) THEN
         cd = cdx(lf)*(xdn(mgs,lf)/rho_qf)**(2./3.)
       ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lf)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )
         
         
         delrho = Max( 0.0, 0.01*(xdn(mgs,lf) - mmgraupvt(indxr,1)) )
         IF ( indxr < ngdnmm ) THEN
          ! interpolate from MM13 table
          axx(mgs,lf) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxx(mgs,lf) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )

          
         ELSE
          axx(mgs,lf) = mmgraupvt(indxr,2)
          bxx(mgs,lf) = mmgraupvt(indxr,3)
         ENDIF
         
         aax = axx(mgs,lf)
         bbx = bxx(mgs,lf)

         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 170.0) ) )
         
       ELSEIF ( icdx <= 0 ) THEN ! 
         aax = ax(lf)
         bbx = bx(lf)
          cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 170.0) ) )
      ELSE
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 170.0) ) )
      ENDIF
       
       cdxgs(mgs,lf) = cd

      IF ( alpha(mgs,lf) .eq. 0.0 .and. icdx > 0 .and. icdx /= 6 ) THEN
!      axx(mgs,lf) =  (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lf)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      axx(mgs,lf) = Sqrt(4.0*xdn(mgs,lf)*gr/(3.0*cd*rho00))
      bxx(mgs,lf) = 0.5
      vtxbar(mgs,lf,1) = (gf4p5/6.0)* rhovt(mgs)*axx(mgs,lf) * Sqrt(xdia(mgs,lf,1)) 

      ELSE
        IF ( icdx /= 6 ) bbx = bx(lf)
        tmp = 4. + alpha(mgs,lf) + bbx
!        IF ( .not. ( -0.5 < alpha(mgs,lf) .and. alpha(mgs,lf) < 16. ) ) THEN
!          write(0,*) 'setvt: problem with alphaf: ',alpha(mgs,lf),qx(mgs,lf),cx(mgs,lf)
!          STOP
!        ENDIF
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
        
        IF ( icdx > 0 .and. icdx /= 6 .and. bbx == 0.5 ) THEN
          aax = Sqrt(4.0*xdn(mgs,lf)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lf,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lf,1)) * x/y
!          vtxbar(mgs,lf,1) =  rhovt(mgs)*aax* xdia(mgs,lf,1)**bbx * x/y
          axx(mgs,lf) = aax
          bxx(mgs,lf) = bbx
        ELSEIF (icdx == 6 ) THEN
          vtxbar(mgs,lf,1) =  rhovt(mgs)*aax* xdia(mgs,lf,1)**bbx * x/y
        ELSE
          axx(mgs,lf) = ax(lf)
          bxx(mgs,lf) = bx(lf)
          vtxbar(mgs,lf,1) =  rhovt(mgs)*ax(lf)*(xdia(mgs,lf,1)**bx(lf)*x)/y
        ENDIF

      ENDIF

      IF ( lwsm6 .and. ipconc == 0 ) THEN
         vtxbar(mgs,lf,1) = (330.*gf4br/6.0)*(xdia(mgs,lf,1)**br)*rhovt(mgs)
      ENDIF

      IF ( infdo >= 1 .and. imydiagalpha == 3 ) THEN
       ! diagnostic shape parameter for Vn (Milbrandt and M-C 2010)
       ! Mk = n0*gamma(k+alpha+1)*xdia(mgs,lh,1)**(k+alpha+1)
       ! M3/M0 = gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3)
       !  tmp = (gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3))**(1./3.)

        tmp = 4. + alpha(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alpha(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = (x/y)**(1./3.)*xdia(mgs,lf,1)

!        alphan(mgs,lf) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
        alp = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)

        tmp = 4. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alp + bxx(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp + bxx(mgs,lf)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        ! vn = vm*gamma(0+1+alpha+b)*gamma(3+1+alpha)/(gamma(3+1+alpha+b)*gamma(0+1+alpha))

!        myfac(mgs,lf) = y2*x/(x2*y)
        ! vn = vm*gamma(0+1+alpha+b)*gamma(3+1+alpha)/(gamma(3+1+alpha+b)*gamma(0+1+alpha))

        vtxbar(mgs,lf,2) = vtxbar(mgs,lf,1)*y2*x/(x2*y)

        ! write(0,*) 'alphan: ',alphan(mgs,lf),qx(mgs,lf)*1.e3,1.e3*xdia(mgs,lf,1),x,y,(x/y)**(1./3.)

      ELSEIF ( infdo >= 1 .and. imydiagalpha == 4 ) THEN
        ! Thompson hack
        alphan(mgs,lf) = alpha(mgs,lf) + 1.5

      ENDIF
      
#ifdef USEMIXEDPHASE
! adjust for liquid fraction, if predicted
!      IF ( ifwmfall > 0 .and. lhw > 1  .and. (xdia(mgs,lf,3) < sheddiam .or. ifwmfall > 1 ) ) THEN ! have liquid fraction
      IF ( ifwmfall > 0 .and. lhw > 1 ) THEN ! have liquid fraction
       IF ( qxw(mgs,lf) > qxmin(lf) ) THEN
         fw = qxw(mgs,lf)/qx(mgs,lf)  ! liquid mass fraction
         IF ( fw > fwlo ) THEN ! adjust toward rain fall speed
           vtrain = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lf,1))**(-alpha(mgs,lf) - 4.0) )
           
!           vtrainnum = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lf,1))**(-alpha(mgs,lf) - 1.0) ) ! number weighted
           
!           vtrainz = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lf,1))**(-alpha(mgs,lf) - 7.0) ) ! z-weighted
           
!          IF ( vtrain <  vtxbar(mgs,lf,1) ) THEN
           IF ( xdia(mgs,lf,1)*(3. + alpha(mgs,lf)) < 9.0e-3 ) THEN ! test on maximum mass diameter against largest rain drop
           
!           IF ( my_rank == 3 ) THEN
!             write(0,*) 'setvtz: qh,qhw,vtr,vtrn,vtrz=',qx(mgs,lf),qx(mgs,lhw),vtrain,vtrainnum,vtrainz
!           ENDIF
           
            IF ( fw > fwhi ) THEN
              vtxbar(mgs,lf,1) = vtrain
!              vtxbar(mgs,lf,2) = vtrainnum
!              vtxbar(mgs,lf,3) = Max(vtrain, vtrainz)
            ELSE
               vtxbar(mgs,lf,1) = vtrain + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lf,1) - vtrain)
!               vtxbar(mgs,lh,2) = vtrainnum + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lf,2) - vtrainnum)
!               vtxbar(mgs,lf,3) = vtrainz + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lf,3) - vtrainz)
!               vtxbar(mgs,lf,3) = Max( vtxbar(mgs,lf,3), vtxbar(mgs,lf,1) )
            ENDIF
          ENDIF
          
         ENDIF
       ENDIF
      ENDIF
#endif
      end if



      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'
      
      ENDIF ! lf .gt. 1

!
!
! ################################################################
!
!  GRAUPEL  FALL SPEED ! Wisner et al. (1972)
!
      IF ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) THEN

      do mgs = 1,ngscnt
      vtxbar(mgs,lh,1) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
       IF ( icdx .eq. 1 ) THEN
         cd = cdx(lh)
       ELSEIF ( icdx .eq. 2 ) THEN
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
         cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!         cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
       ELSEIF ( icdx .eq. 3 ) THEN
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
       ELSEIF ( icdx .eq. 4 ) THEN
         cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
     &        (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
       ELSEIF ( icdx .eq. 5 ) THEN
         cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
       ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lh)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )
         
         
         delrho = Max( 0.0, 0.01*(xdn(mgs,lh) - mmgraupvt(indxr,1)) )
         IF ( indxr < ngdnmm ) THEN
          
          axx(mgs,lh) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxx(mgs,lh) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )

          
         ELSE
          axx(mgs,lh) = mmgraupvt(indxr,2)
          bxx(mgs,lh) = mmgraupvt(indxr,3)
         ENDIF
         
         aax = axx(mgs,lh)
         bbx = bxx(mgs,lh)

         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
         
       ELSEIF ( icdx <= 0 ) THEN ! 
         aax = ax(lh)
         bbx = bx(lh)
          cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
      ELSE
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
      ENDIF
       
       cdxgs(mgs,lh) = cd

      IF ( alpha(mgs,lh) .eq. 0.0 .and. icdx > 0 .and. icdx /= 6 ) THEN
!      axx(mgs,lh) =  (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lh)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      axx(mgs,lh) = Sqrt(4.0*xdn(mgs,lh)*gr/(3.0*cd*rho00))
      bxx(mgs,lh) = 0.5
      vtxbar(mgs,lh,1) = (gf4p5/6.0)* rhovt(mgs)*axx(mgs,lh) * Sqrt(xdia(mgs,lh,1)) 
!      vtxbar(mgs,lh,1) = (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lh)*xdia(mgs,lh,1)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      ELSE
        IF ( icdx /= 6 ) bbx = bx(lh)
        tmp = 4. + alpha(mgs,lh) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
!        aax = Max( 1.0, Min(2.0, (xdn(mgs,lh)/400.) ) )
!        vtxbar(mgs,lh,1) =  rhovt(mgs)*aax*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y
        
        IF ( icdx > 0 .and. icdx /= 6 .and. bbx == 0.5 ) THEN
          aax = Sqrt(4.0*xdn(mgs,lh)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lh,1)) * x/y
          axx(mgs,lh) = aax
          bxx(mgs,lh) = bbx
        ELSEIF (icdx == 6 ) THEN
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* xdia(mgs,lh,1)**bbx * x/y
        ELSE
          axx(mgs,lh) = ax(lh)
          bxx(mgs,lh) = bx(lh)
          vtxbar(mgs,lh,1) =  rhovt(mgs)*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y
        ENDIF

!     &    Gamma(4.0 + dnu(lh) + 0.6))/Gamma(4. + dnu(lh))
      ENDIF

      IF ( lwsm6 .and. ipconc == 0 ) THEN
!         vtxbar(mgs,lh,1) = (330.*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         vtxbar(mgs,lh,1) = (330.*gf4br/6.0)*(xdia(mgs,lh,1)**br)*rhovt(mgs)
      ENDIF
      
      IF ( infdo >= 1 .and. imydiagalpha == 3 ) THEN
       ! diagnostic shape parameter for Vn (Milbrandt and M-C 2010)
       ! Mk = n0*gamma(k+alpha+1)*xdia(mgs,lh,1)**(k+alpha+1)
       ! M3 = n0*gamma(alpha+4)*xdia(mgs,lh,1)**(alpha+4)
       ! M0 = n0*gamma(alpha+1)*xdia(mgs,lh,1)**(alpha+1)
       ! M3/M0 = gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3)
       !  tmp = (gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3))**(1./3.)
       !  tmp = (gamma(alpha+4)/gamma(alpha+1))**(1./3.)*xdia(mgs,lh,1)

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp2 = (x/y)**(1./3.)*xdia(mgs,lh,1)

!        alphan(mgs,lh) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
        alp = Min(15., 11.8*(1000.*tmp2 - 0.7)**2 + 2.)

!        IF ( igs(mgs) == 10 ) THEN
!          write(0,*) 'k,xdia,mu,tmp = ',kgs(mgs),xdia(mgs,lh,1),alp,tmp
!        ENDIF

        tmp = 4. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alp + bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp + bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        ! vn = vm*gamma(0+1+alpha+b)*gamma(3+1+alpha)/(gamma(3+1+alpha+b)*gamma(0+1+alpha))

        vtxbar(mgs,lh,2) = vtxbar(mgs,lh,1)*y2*x/(x2*y)
!        vtxbar(mgs,lh,3) = vtxbar(mgs,lh,1)
        myfac(mgs,lh) = y2*x/(x2*y)
!         IF ( igs(mgs) == 10 ) THEN
!           write(0,*) 'k,xdia,mu,tmp2,myfac = ',kgs(mgs),xdia(mgs,lh,1),alp,tmp2,myfac(mgs,lh)
!         ENDIF
!         write(0,*) 'alphan,vm,vn,qh,tem: ',alphan(mgs,lh),vtxbar(mgs,lh,1),vtxbar(mgs,lh,2), &
!                       qx(mgs,lh)*1.e3,temcg(mgs),1.e3*xdia(mgs,lh,1),x,y,(x/y)**(1./3.)

      ELSEIF ( infdo >= 1 .and. imydiagalpha == 4 ) THEN
        ! Thompson hack
        alphan(mgs,lh) = alpha(mgs,lh) + 1.5

      ENDIF
      
#ifdef USEMIXEDPHASE
! adjust for liquid fraction, if predicted
!      IF ( ifwmfall > 0 .and. lhw > 1  .and. (xdia(mgs,lh,3) < sheddiam .or. ifwmfall > 1 ) ) THEN ! have liquid fraction
      IF ( ifwmfall > 0 .and. lhw > 1 ) THEN ! have liquid fraction
       IF ( qxw(mgs,lh) > qxmin(lh) ) THEN
         fw = qxw(mgs,lh)/qx(mgs,lh)  ! liquid mass fraction
         IF ( fw > fwlo ) THEN ! adjust toward rain fall speed
           vtrain = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lh,1))**(-alpha(mgs,lh) - 4.0) )
           
!           vtrainnum = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lh,1))**(-alpha(mgs,lh) - 1.0) ) ! number weighted
           
!           vtrainz = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lh,1))**(-alpha(mgs,lh) - 7.0) ) ! z-weighted
           
!          IF ( vtrain <  vtxbar(mgs,lh,1) ) THEN
           IF ( xdia(mgs,lh,1)*(3. + alpha(mgs,lh)) < 9.0e-3 ) THEN ! test on maximum mass diameter against largest rain drop
           
!           IF ( my_rank == 3 ) THEN
!             write(0,*) 'setvtz: qh,qhw,vtr,vtrn,vtrz=',qx(mgs,lh),qx(mgs,lhw),vtrain,vtrainnum,vtrainz
!           ENDIF
           
            IF ( fw > fwhi ) THEN
              vtxbar(mgs,lh,1) = vtrain
!              vtxbar(mgs,lh,2) = vtrainnum
!              vtxbar(mgs,lh,3) = Max(vtrain, vtrainz)
            ELSE
               vtxbar(mgs,lh,1) = vtrain + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lh,1) - vtrain)
!               vtxbar(mgs,lh,2) = vtrainnum + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lh,2) - vtrainnum)
!               vtxbar(mgs,lh,3) = vtrainz + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lh,3) - vtrainz)
!               vtxbar(mgs,lh,3) = Max( vtxbar(mgs,lh,3), vtxbar(mgs,lh,1) )
            ENDIF
          ENDIF
          
         ENDIF
       ENDIF
      ENDIF
#endif
      end if



      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'
      
      ENDIF ! lh .gt. 1
!
!
! ################################################################
!
!  HAIL
!
      IF ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) THEN
      
      do mgs = 1,ngscnt
      vtxbar(mgs,lhl,1) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then

       IF ( icdxhl .eq. 1 ) THEN
         cd = cdx(lhl)
       ELSEIF ( icdxhl .eq. 3 ) THEN
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
       ELSEIF ( icdxhl .eq. 4 ) THEN
         cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
     &       (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
       ELSEIF ( icdxhl .eq. 5 ) THEN
         cd = cdx(lh)*(xdn(mgs,lhl)/rho_qh)**(2./3.)
       ELSEIF ( icdxhl .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lhl)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )
         
         
         delrho = Max( 0.0, 0.01*(xdn(mgs,lhl) - mmgraupvt(indxr,1)) )
         IF ( indxr < ngdnmm ) THEN
          
          axx(mgs,lhl) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxx(mgs,lhl) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )

          
         ELSE
          axx(mgs,lhl) = mmgraupvt(indxr,2)
          bxx(mgs,lhl) = mmgraupvt(indxr,3)
         ENDIF
         
         aax = axx(mgs,lhl)
         bbx = bxx(mgs,lhl)

         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
         
       ELSE
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lhl))/(900. - 300.) ) )
!        cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lhl) - xdn(mgs,lhl))/(xdnmx(lhl)-xdnmn(lhl)) ) )
!         cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
       ENDIF

       cdxgs(mgs,lhl) = cd

      IF ( alpha(mgs,lhl) .eq. 0.0 .and. icdxhl > 0 .and. icdxhl /= 6) THEN
!      axx(mgs,lhl) =  (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lhl)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      axx(mgs,lhl) = Sqrt(4.0*xdn(mgs,lhl)*gr/(3.0*cd*rho00))
      bxx(mgs,lhl) = 0.5
      vtxbar(mgs,lhl,1) = (gf4p5/6.0)* rhovt(mgs)*axx(mgs,lhl) * Sqrt(xdia(mgs,lhl,1)) 
      ELSE
        IF ( icdxhl /= 6 ) bbx = bx(lhl)
        tmp = 4. + alpha(mgs,lhl) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
!        IF ( i+1 > ngm0 ) THEN
!         write(0,*) 'SETVTZ: i out of bounds for gmoi',i,del,tmp,alpha(mgs,lhl),bbx
!        ENDIF
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        IF ( icdxhl > 0 .and. icdxhl /= 6 .and. bbx == 0.5 ) THEN
          aax = Sqrt(4.0*xdn(mgs,lhl)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lhl,1)) * x/y
          axx(mgs,lhl) = aax
          bxx(mgs,lhl) = bbx
        ELSEIF ( icdxhl == 6 ) THEN
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* (xdia(mgs,lhl,1))**bbx * x/y
        ELSE
          axx(mgs,lhl) = ax(lhl)
          bxx(mgs,lhl) = bx(lhl)
         vtxbar(mgs,lhl,1) =  rhovt(mgs)*(ax(lhl)*xdia(mgs,lhl,1)**bx(lhl)*x)/y
        ENDIF
        
!     &    Gamma(4.0 + dnu(lh) + 0.6))/Gamma(4. + dnu(lh))
      ENDIF

      IF ( infdo >= 1 .and. imydiagalpha == 3 ) THEN
       ! diagnostic shape parameter for Vn (Milbrandt and M-C 2010)
       ! Mk = n0*gamma(k+alpha+1)*xdia(mgs,lh,1)**(k+alpha+1)
       ! M3/M0 = gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3)
       !  tmp = (gamma(alpha+4)/gamma(alpha+1)*xdia(mgs,lh,1)**(3))**(1./3.)

        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = (x/y)**(1./3.)*xdia(mgs,lhl,1)

!        alphan(mgs,lhl) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
        alp = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)

        tmp = 4. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alp + bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 1. + alp + bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        ! vn = vm*gamma(0+1+alpha+b)*gamma(3+1+alpha)/(gamma(3+1+alpha+b)*gamma(0+1+alpha))

        vtxbar(mgs,lhl,2) = vtxbar(mgs,lhl,1)*y2*x/(x2*y)
!        vtxbar(mgs,lhl,3) = vtxbar(mgs,lhl,1)
        myfac(mgs,lhl) = y2*x/(x2*y)

      ELSEIF ( infdo >= 1 .and. imydiagalpha == 4 ) THEN
        ! Thompson hack
        alphan(mgs,lhl) = alpha(mgs,lhl) + 1.5

      ENDIF

#ifdef USEMIXEDPHASE
! adjust for liquid fraction, if predicted
!      IF ( ifwmfall > 0 .and. lhlw > 1 .and. (xdia(mgs,lhl,3) < sheddiam .or. ifwmfall > 1 ) ) THEN ! have liquid fraction
      IF ( ifwmfall > 1 .and. lhlw > 1 ) THEN ! have liquid fraction
!       IF ( ifwmfall == 2 ) THEN
!       
!       ENDIF
       IF ( qxw(mgs,lhl) > qxmin(lhl) ) THEN
         fw = qxw(mgs,lhl)/qx(mgs,lhl)  ! liquid mass fraction
         IF ( fw > fwlo ) THEN ! adjust toward rain fall speed
           vtrain = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lhl,1))**(-alpha(mgs,lhl) - 4.0) )

!           vtrainnum = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lhl,1))**(-alpha(mgs,lhl) - 1.0) ) ! number weighted
           
!           vtrainz = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lhl,1))**(-alpha(mgs,lhl) - 7.0) ) ! z-weighted
!          IF ( vtrain <  vtxbar(mgs,lhl,1) ) THEN
           IF ( xdia(mgs,lhl,1)*(3. + alpha(mgs,lhl)) < 9.0e-3 ) THEN ! test on maximum mass diameter against largest rain drop
            IF ( fw > fwhi ) THEN
              vtxbar(mgs,lhl,1) = vtrain
!              vtxbar(mgs,lhl,2) = vtrainnum
!              vtxbar(mgs,lhl,3) = Max(vtrain, vtrainz)
            ELSE
               vtxbar(mgs,lhl,1) = vtrain + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lhl,1) - vtrain)
!               vtxbar(mgs,lhl,2) = vtrainnum + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lhl,2) - vtrainnum)
!               vtxbar(mgs,lhl,3) = vtrainz + rfwdiff*(fwhi - fw)*(vtxbar(mgs,lhl,3) - vtrainz)
            ENDIF
          ENDIF
          
         ENDIF
       ENDIF
      ENDIF
#endif

      end if


      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'
      
      ENDIF ! lhl .gt. 1


      IF ( infdo .ge. 1 ) THEN

!      DO il = lc,lhab
!      IF ( il .ne. lr ) THEN
        DO mgs = 1,ngscnt
          IF ( ildo == 0 .or. ildo == lc ) THEN
            vtxbar(mgs,lc,2) = vtxbar(mgs,lc,1)
          ENDIF
        IF ( li .gt. 1 ) THEN
!          vtxbar(mgs,li,2) = rhovt(mgs)*49420.*1.25447*xdia(mgs,li,1)**(1.415) ! n-wgt (Ferrier 94)
!          vtxbar(mgs,li,2) = vtxbar(mgs,li,1)

! test print stuff...
!          IF ( xdia(mgs,li,1) .gt. 200.e-6 ) THEN
!            tmp = (xv(mgs,li)*cwc0)**(1./3.)
!            x = rhovt(mgs)*49420.*40.0005/5.40662*tmp**(1.415)
!            y = rhovt(mgs)*49420.*1.25447*tmp**(1.415)
!            write(6,*) 'Ice fall: ',vtxbar(mgs,li,1),x,y,tmp,xdia(mgs,li,1)
!          ENDIF
        ENDIF
!          vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        ENDDO

!        IF ( lg .gt. lr ) THEN 
        IF ( lg .gt. lr .and. .not. ( imydiagalpha == 3 .and. ipconc == 5 )) THEN

        DO il = lg,lhab
         IF ( ildo == 0 .or. ildo == il ) THEN

                alphahltmp = alphahl

            DO mgs = 1,ngscnt
             IF ( qx(mgs,il) .gt. qxmin(il) ) THEN
             
               IF ( lhl > 1 ) alphahltmp = alphan(mgs,lhl)
              
              IF ( il==lf .or. (il .eq. lh .and. ( hssflg == 1 .or. (hssflg == -1 .and. temcg(mgs) <= -0.5)  ) ) .or.   &
                ( lhl .gt. 1 .and. il .eq. lhl .and. ( hlssflg == 1 .or. (hlssflg == -1 .and. temcg(mgs) <= -0.5) .or. &
                  (hlssflg == 2 .and. alphahltmp < alphamax ) ) ) ) THEN ! DTD: added flag for size-sorting
              
              ! DTD: allow for setting of number-weighted and z-weighted fall speeds to the mass-weighted value,
              ! effectively turning off size-sorting

              IF ( il .eq. lh ) THEN ! {
             
               IF ( icdx .eq. 1 ) THEN
                 cd = cdx(lh)
               ELSEIF ( icdx .eq. 2 ) THEN
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
                 cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!                 cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
               ELSEIF ( icdx .eq. 3 ) THEN
!                 cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 170.0, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
                 cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
               ELSEIF ( icdx .eq. 4 ) THEN
                 cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
     &            (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
               ELSEIF ( icdx .eq. 5 ) THEN
                 cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
               ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
                  aax = axx(mgs,lh)
                  bbx = bxx(mgs,lh)
               ELSEIF ( icdx <= 0 ) THEN ! 
                  aax = ax(lh)
                  bbx = bx(lh)
               ENDIF
               
              ELSEIF ( lf .gt. 1 .and. il .eq. lf ) THEN
             
               IF ( icdx .eq. 1 ) THEN
                 cd = cdx(lf)
               ELSEIF ( icdx .eq. 3 ) THEN
!               cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 300.) ) )
                cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 170.0) ) )
               ELSEIF ( icdx .eq. 4 ) THEN
                cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)*  &
     &               (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lf) ) ) )/(cdhdnmax - cdhdnmin) ) )
               ELSEIF ( icdx == 5 ) THEN
!                cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lf))/(900. - 300.) ) )
!                cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lf) - xdn(mgs,lf))/(xdnmx(lf)-xdnmn(lf)) ) )
                 cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lf) ) ) )/(800. - 500.) ) )
               ELSEIF ( icdx .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
                  aax = axx(mgs,lf)
                  bbx = bxx(mgs,lf)
               ELSEIF ( icdx <= 0 ) THEN ! 
                  aax = ax(lf)
                  bbx = bx(lf)
               ENDIF

              ELSEIF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
             
               IF ( icdxhl .eq. 1 ) THEN
                 cd = cdx(lhl)
               ELSEIF ( icdxhl .eq. 3 ) THEN
!               cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 300.) ) )
                cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
               ELSEIF ( icdxhl .eq. 4 ) THEN
                cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
     &               (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
               ELSEIF ( icdxhl == 5 ) THEN
!                cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lhl))/(900. - 300.) ) )
!                cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lhl) - xdn(mgs,lhl))/(xdnmx(lhl)-xdnmn(lhl)) ) )
                 cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
               ELSEIF ( icdxhl .eq. 6 ) THEN ! Milbrandt and Morrison (2013)
                  aax = axx(mgs,lhl)
                  bbx = bxx(mgs,lhl)
               ELSEIF ( icdxhl <= 0 ) THEN ! 
                  aax = ax(lhl)
                  bbx = bx(lhl)
               ENDIF
               
              ENDIF ! }

               IF ( alphan(mgs,il) .eq. 0. .and. infdo .lt. 2 .and.   &
               ( ( (il==lh .or. il==lf) .and. icdx > 0 .and. icdx /= 6) .or. ( il==lhl .and. icdxhl > 0 .and. icdxhl /= 6 ) ) ) THEN ! {
                 vtxbar(mgs,il,2) =   &
     &              Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) / &
     &                (3.0*cd*rho0(mgs)) )

               ELSE
               IF ( ( il == lh .or. il == lf )  .and. icdx   /= 6 ) bbx = bx(il)
               IF ( il == lhl .and. icdxhl /= 6 ) bbx = bx(il)
               tmp = 1. + alphan(mgs,il) + bbx
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
  
               tmp = 1. + alphan(mgs,il)
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

                 IF ( il .eq. lh .or. il == lf .or. il .eq. lhl) THEN ! {
                   IF ( ( ( il==lh .or. il==lf ) .and. icdx > 0 ) ) THEN
                     IF ( icdx /= 6 ) THEN
                      aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                      vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bx(il) * x/y
                     ELSE !  (icdx == 6 ) THEN
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     ENDIF

                   ELSEIF ( ( il==lhl .and. icdxhl > 0 ) ) THEN
                     IF ( icdxhl /= 6 ) THEN
                       aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bx(il) * x/y
                     ELSE ! ( icdxhl == 6 )
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     ENDIF
                   
                   ELSE  ! get here if il==lh and icdx <= 0 -- or -- il==lhl and icdxhl <= 0
                     aax = ax(il)
                     vtxbar(mgs,il,2) =  rhovt(mgs)*ax(il)*(xdia(mgs,il,1)**bx(il)*x)/y
                   ENDIF
                   
!                   IF ( il == lh .and. imydiagalpha == 3 .and. ipconc == 5 ) THEN
!                    vtxbar(mgs,il,2) = myfac(mgs,lh)*vtxbar(mgs,il,1)
!                   ENDIF

!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(75.715*xdia(mgs,il,1)**0.6* &
!     &               x)/y
!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
!     &               x)/y
                  IF ( infdo .ge. 2 ) THEN ! Z-weighted

               tmp = 7. + alpha(mgs,il) + bbx
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
  
               tmp = 7. + alpha(mgs,il)
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

                   vtxbar(mgs,il,3) = rhovt(mgs)*                 &
     &                (aax*(xdia(mgs,il,1) )**bbx *  &
     &                 x)/y
!     &                 Gamma(7.0 + alpha(mgs,il) + bbx))/Gamma(7. + alpha(mgs,il))

!            IF ( igs(mgs) == 10 .and. il == lhl ) THEN
!              write(91,*) 'k,xdia,mu,vt3 = ',kgs(mgs),xdia(mgs,lhl,1),alpha(mgs,il), &
!                     vtxbar(mgs,il,3)
!            ENDIF

          IF ( .not. (vtxbar(mgs,il,1) > -1. .and. vtxbar(mgs,il,1) < 250. ) .or. &
               .not. (vtxbar(mgs,il,3) > -1. .and. vtxbar(mgs,il,3) < 250. ) ) THEN
           write(0,*) 'Setvtz: problem with vtxbar1/3: ',il,igs(mgs),kgs(mgs),vtxbar(mgs,il,1),vtxbar(mgs,il,3),aax,bbx,x,y
           write(0,*) 'qx,cx,alpha,xdia1,3 = ',qx(mgs,il),cx(mgs,il),alpha(mgs,il),xdia(mgs,il,1),xdia(mgs,il,3)
           write(0,*) 'cd = ',cd,gamma_sp(1. + alpha(mgs,il) + bbx),gamma_sp(1. + alpha(mgs,il)),x/y
          ! call commasmpi_abort()
          ENDIF
!     &                (aax*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &                 Gamma(7.0 + alpha(mgs,il) + bx(il)))/Gamma(7. + alpha(mgs,il))
                  ENDIF

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt3'

                 ELSE ! hail -- I think the code should never get here?
                  vtxbar(mgs,il,2) =  &
     &               rhovt(mgs)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
     &               x)/y

                 IF ( infdo .ge. 2 ) THEN ! Z-weighted
                  vtxbar(mgs,il,3) = rhovt(mgs)*                 &
     &              (aax*(1.0/xdia(mgs,il,1) )**(- bbx)*  &
     &               Gamma(7.0 + alpha(mgs,il) + bbx))/Gamma(7. + alpha(mgs,il))
!     &              (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &               Gamma(7.0 + alpha(mgs,il) + bx(il)))/Gamma(7. + alpha(mgs,il))
                  ENDIF

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt4'

                 ENDIF ! }
!     &             Gamma(1.0 + dnu(il) + 0.6)/Gamma(1. + dnu(il))
               ENDIF ! }

#ifdef USEMIXEDPHASE
! adjust for liquid fraction, if predicted
      IF ( ( il == lh .and. lhw > 1 .and. ifwmfall > 0) .or.  &
     &     ( il == lhl .and. lhlw > 1 .and. ifwmfall > 1 ) ) THEN ! have liquid fraction
!       IF ( ifwmfall > 0 .and. qxw(mgs,il) > qxmin(il) .and. (xdia(mgs,il,3) < sheddiam .or. ifwmfall > 1 ) ) THEN
       IF ( ifwmfall > 0 .and. qxw(mgs,il) > qxmin(il) ) THEN
         fw = qxw(mgs,il)/qx(mgs,il)  ! liquid mass fraction
         IF ( fw > fwlo ) THEN ! adjust toward rain fall speed
          vtrain = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,il,1))**(-alphan(mgs,il) - 1.0) ) ! number weighted
!          IF ( vtrain <  vtxbar(mgs,il,2) ) THEN
          IF ( xdia(mgs,il,1)*(3. + alpha(mgs,il)) < 9.0e-3 ) THEN ! test on maximum mass diameter against largest rain drop
            IF ( fw > fwhi ) THEN
              vtxbar(mgs,il,2) = vtrain
            ELSE
               vtxbar(mgs,il,2) = vtrain + rfwdiff*(fwhi - fw)*(vtxbar(mgs,il,2) - vtrain)
            ENDIF
          ENDIF
          IF ( infdo .ge. 2 ) THEN ! Z-weighted
!          IF ( vtrain <  vtxbar(mgs,il,3) ) THEN
           IF ( xdia(mgs,il,1)*(3. + alpha(mgs,il)) < 9.0e-3 ) THEN !{ test on maximum mass diameter against largest rain drop
            
            vtrainz = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,il,1))**(-alpha(mgs,il) - 7.0) ) ! z-weighted
            
            IF ( fw > fwhi ) THEN
              vtxbar(mgs,il,3) = vtrainz
            ELSE
               vtxbar(mgs,il,3) = vtrainz + rfwdiff*(fwhi - fw)*(vtxbar(mgs,il,3) - vtrainz)
            ENDIF
            vtxbar(mgs,il,3) = Max( vtxbar(mgs,il,3), vtxbar(mgs,il,1) )
          
           ENDIF !}
          ELSE
            vtxbar(mgs,il,3) = vtxbar(mgs,il,1) 
          
          ENDIF
          
         ENDIF
       ENDIF
      ENDIF
#endif
!              IF ( infdo .ge. 2 ) THEN ! Z-weighted
!               vtxbar(mgs,il,3) = rhovt(mgs)*                 &
!     &            (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &             Gamma(7.0 + alpha(mgs,il) + bx(il)))/Gamma(7. + alpha(mgs,il))
!              ENDIF

!               IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!                write(0,*) 'setvt: ',qx(mgs,il),xdia(mgs,il,1),xdia(mgs,il,3),dnu(il),ax(il),bx(il)
!               ENDIF
             ELSEIF ( ( (il .eq. lh .or. il==lf) .and. ( hssflg == 0 .or. (hssflg == -1 .and. temcg(mgs) > -0.5) ) ) .or.   &
             ( lhl .gt. 1 .and. il .eq. lhl .and. (hlssflg == 0 .or. (hlssflg == -1 .and. temcg(mgs) > -0.5) .or. &
               ( hlssflg == 2 .and. alphahltmp >= alphamax )) ) ) THEN ! no size-sorting for graupel or hail
              vtxbar(mgs,il,2) = vtxbar(mgs,il,1)
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)
             ELSE ! not lh or lhl
              vtxbar(mgs,il,2) = &
     &            Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) /  &
     &              (3.0*cdx(il)*rho0(mgs)) )
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt5'


              ENDIF
             ELSE ! qx < qxmin
              vtxbar(mgs,il,2) = 0.0

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt6'

             ENDIF
         
           ENDDO ! mgs

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt7'

        ENDIF ! ildo
        ENDDO ! il

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt8'

        ENDIF ! lg .gt. 1 

        

        
!      ENDIF
!      ENDDO

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt9'

!       DO mgs = 1,ngscnt
!        IF ( qx(mgs,lr) > qxmin(lr) ) THEN
!         write(0,*) 'setvt2: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
!        ENDIF
!       ENDDO

      ENDIF ! infdo .ge. 1 


       IF ( ipconc == 5 .and. infall == 0 .and. infalln > 0 .and. infdo > 0 ) THEN
           DO mgs = 1,ngscnt
              vtxbar(mgs,lh,2) = 1.7*vtxbar(mgs,lh,2)
           ENDDO
           IF ( irfall == 0 ) THEN
           DO mgs = 1,ngscnt
              vtxbar(mgs,lr,2) = 1.7*vtxbar(mgs,lr,2)
           ENDDO
           ENDIF

           IF ( lhl > 0 ) THEN
           DO mgs = 1,ngscnt
              vtxbar(mgs,lhl,2) = 1.7*vtxbar(mgs,lhl,2)
           ENDDO
           ENDIF
       ENDIF

       IF ( ( ifallsedonly == 0 ) .or. (  ifallsedonly == 1 .and. infdo > 0 ) .or. &
             ( ifallsedonly == 2 .and. infdo == 0 ) ) THEN
         IF (  lh > 0 .and. graupelfallfac /= 1.0 ) THEN
           DO mgs = 1,ngscnt
              vtxbar(mgs,lh,1) = graupelfallfac*vtxbar(mgs,lh,1)
              vtxbar(mgs,lh,2) = graupelfallfac*vtxbar(mgs,lh,2)
              vtxbar(mgs,lh,3) = graupelfallfac*vtxbar(mgs,lh,3)
              axx(mgs,lh) = graupelfallfac*axx(mgs,lh)
            ENDDO
          ENDIF

         IF (  lf > 0 .and. fdfallfac /= 1.0 ) THEN
           DO mgs = 1,ngscnt
              vtxbar(mgs,lf,1) = fdfallfac*vtxbar(mgs,lf,1)
              vtxbar(mgs,lf,2) = fdfallfac*vtxbar(mgs,lf,2)
              vtxbar(mgs,lf,3) = fdfallfac*vtxbar(mgs,lf,3)
              axx(mgs,lf) = fdfallfac*axx(mgs,lf)
            ENDDO
          ENDIF

          IF ( lhl > 0 .and. hailfallfac /= 1.0 ) THEN
            DO mgs = 1,ngscnt
              vtxbar(mgs,lhl,1) = hailfallfac*vtxbar(mgs,lhl,1)
              vtxbar(mgs,lhl,2) = hailfallfac*vtxbar(mgs,lhl,2)
              vtxbar(mgs,lhl,3) = hailfallfac*vtxbar(mgs,lhl,3)
              axx(mgs,lhl) = hailfallfac*axx(mgs,lhl)
            ENDDO
          ENDIF

        ENDIF
     
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: END OF ROUTINE'

!############ SETVTZ ############################

      RETURN
      END
!--------------------------------------------------------------------------

!
! ##############################################################################
!
      SUBROUTINE setqxminz(qxmin)
      
      USE INDEX_MODULE

      implicit none
      
      real qxmin(lc:lhab)

      
      qxmin(:) = 1.0e-12
      
      qxmin(lc) = 1.e-9
      IF ( lis > 1 ) qxmin(lis) = 1.e-13
      qxmin(lr) = 1.e-7
      IF ( li > 1 ) qxmin(li) = 1.e-12
      IF ( ls > 1 ) qxmin(ls) = 1.e-7
      IF ( lh > 1 ) qxmin(lh) = 1.e-7
      IF ( lf > 1 ) qxmin(lf) = 1.e-7
      IF ( lhl .gt. 1 ) qxmin(lhl) = 1.e-7
      
      IF ( lc .gt. 1 .and. lnc .gt. 1 ) qxmin(lc) = 1.0e-13
      IF ( lr .gt. 1 .and. lnr .gt. 1 ) qxmin(lr) = 1.0e-12

      IF ( li .gt. 1 .and. lni .gt. 1 ) qxmin(li ) = 1.0e-13
      IF ( ls .gt. 1 .and. lns .gt. 1 ) qxmin(ls ) = 1.0e-13
      IF ( lh .gt. 1 .and. lnh .gt. 1 ) qxmin(lh ) = 1.0e-12
      IF ( lf .gt. 1 .and. lnf .gt. 1 ) qxmin(lf ) = 1.0e-12
      IF ( lhl.gt. 1 .and. lnhl.gt. 1 ) qxmin(lhl) = 1.0e-12
!      IF ( lzhl > 1 ) qxmin(lhl) = 1.0e-9
      
      RETURN
      END

!
! ##############################################################################
!
      SUBROUTINE setqxmin(qxmin)

      USE INDEX_MODULE

      implicit none

      real qxmin(lc:lhab)


      qxmin(:) = 1.0e-12

!      qxmin(lc) = 1.0e-09
!      IF ( lr  .ge. lc ) qxmin(lr) = 1.0e-10
!      IF ( li  .ge. lc ) qxmin(li) = 1.0e-12
!      IF ( lir .ge. lc ) qxmin(lir) = 1.0e-12
!      IF ( ls  .ge. lc ) qxmin(ls) = 1.0e-11
!      IF ( lgm .ge. lc ) qxmin(lgm) = 1.0e-09
!      IF ( lh  .ge. lc ) qxmin(lh) = 1.0e-09
!      IF ( lgl .ge. lc ) qxmin(lgl) = 1.0e-09
!      IF ( lgh .ge. lc ) qxmin(lgh) = 1.0e-09
!      IF ( lf  .ge. lc ) qxmin(lf) = 1.0e-09
!      IF ( lhl .ge. lc ) qxmin(lhl) = 1.0e-09
!      IF ( lip .ge. lc ) qxmin(lip) = 1.0e-09

      RETURN
      END


!
! ##############################################################################
!


!
! ##############################################################################
!
      SUBROUTINE sett09s(nx,ny,nz,nor,norz,na,lt,   &
     &             t77,pinit,p2,an,pn,pb,           &
     &          t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00 )

! #ifdef MPI
      USE COMMASMPI_MODULE
! #endif
      implicit none

      integer nx,ny,nz,nor,norz,na,lt

      real t00(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t77(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t3(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t4(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t5(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t6(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t8(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t9(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real pb(-norz+1:nz+norz)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)

      real pinit(-norz+1:nz+norz)
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real cp, rd
      parameter ( cp = 1004.0, rd = 287.04 )

      real poo,cap
      parameter ( cap = rd/cp, poo = 1.0e+05 )

      integer ix,jy,kz

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

! #ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.
! #endif

!      print*,'sett09s: lt = ',lt

      t1(:,:,:) = 0.0
      t2(:,:,:) = 0.0
      t3(:,:,:) = 0.0
      t4(:,:,:) = 0.0
      t5(:,:,:) = 0.0
      t6(:,:,:) = 0.0
      t7(:,:,:) = 0.0
      t8(:,:,:) = 0.0
      t9(:,:,:) = 0.0
      t77(:,:,:) = 0.0

      kzb = -1
      kze = ktile+1
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO kz = kzb,kze
      DO jy = jyb,jye
      DO ix = ixb,ixe
        t77(ix,jy,kz) = pinit(kz) + p2(ix,jy,kz)
        t0(ix,jy,kz) = an(ix,jy,kz,lt)*t77(ix,jy,kz)
        t00(ix,jy,kz) = 380.0/(pn(ix,jy,kz)+pb(kz))
      END DO
      END DO
      END DO
      

      RETURN
      END
      

!
! #####################################################################
!
!
! #####################################################################
      Function delbk(bb,nu,mu,k)
!   
!  Purpose: Caluculates collection coefficients following Siefert (2006)
!
!  delbk is equation (90) (b collecting b -- self-collection)
!  mass-diameter relationship: D = a*x**(b), where x = particle mass
!  general distribution: n(x) = A*x**(nu)*Exp(-lam*x**(mu))
!  where
!      A = mu*N/(Gamma((nu+1)/mu)) *lam**((nu+1)/mu)
!
!      lam = ( Gamma((nu+1)/mu)/Gamma((nu+2)/mu) * xbar )**(-mu)
!
!     where  xbar = L/N  (mass content)/(number concentration) = q*rhoa/N
!

      USE MICRO_MODULE

      implicit none
      real delbk, gamma
      real nu, mu, bb
      integer k
      
      real tmp, del
      real x1, x2, x3, x4
      integer i

        tmp = ((1.0 + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((2.0 + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1.0 + 2.0*bb + k + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x3 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
      
!      delbk =  &
!     &  ((Gamma((1.0 + nu)/mu)/Gamma((2.0 + nu)/mu))**(2.0*bb + k)* &
!     &    Gamma((1.0 + 2.0*bb + k + nu)/mu))/Gamma((1.0 + nu)/mu)

      delbk =  &
     &  ((x1/x2)**(2.0*bb + k)* &
     &    x3)/x1
      
      RETURN
      END
      
! #####################################################################
!
!
! #####################################################################
! Equation (91) in Seifert and Beheng (2006) ("a" collecting "b")
      Function delabk(ba,bb,nua,nub,mua,mub,k)

      USE COMMASMPI_MODULE, only: commasmpi_abort
      USE MICRO_MODULE
      
      implicit none
      real delabk, gamma
      real nua, mua, ba
      integer k
      real nub, mub, bb
      
      integer i
      real tmp,del
      
      real g1pnua, g2pnua, g1pbapnua, g1pbbpk, g1pnub, g2pnub
      
        tmp = (1. + nua)/mua
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
!        IF ( i+1 > ngm0 ) THEN
!          write(0,*) 'delabk: i+1 > ngm0!!!!',i,ngm0,nua,mua,tmp
!          call commasmpi_abort()
!        ENDIF
        g1pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!        write(91,*) 'delabk: g1pnua,gamma = ',g1pnua,Gamma((1. + nua)/mua)

        tmp = ((2. + nua)/mua)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g2pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + ba + nua)/mua)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pbapnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((2 + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g2pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + bb + k + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pbbpk = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

      delabk =  &
     &  (2.*(g1pnua/g2pnua)**ba*     &
     &    g1pbapnua*                                               &
     &    (g1pnub/g2pnub)**(bb + k)*                                &
     &    g1pbbpk)/                                                &
     &  (g1pnua*g1pnub)              
      
!      delabk =  &
!     &  (2.*(Gamma((1. + nua)/mua)/Gamma((2. + nua)/mua))**ba*     &
!     &    Gamma((1. + ba + nua)/mua)*                              &
!     &    (Gamma((1. + nub)/mub)/Gamma((2 + nub)/mub))**(bb + k)*  &
!     &    Gamma((1. + bb + k + nub)/mub))/                         &
!     &  (Gamma((1. + nua)/mua)*Gamma((1. + nub)/mub))              
! old with reversed a and b:
!!     :   (2*(Gamma((1 + nua)/mua)/Gamma((2 + nua)/mua))**(ba + k)* &
!!     :    Gamma((1 + ba + k + nua)/mua)*                           &
!!     :    (Gamma((1 + nub)/mub)/Gamma((2 + nub)/mub))**bb*         &
!!     :    Gamma((1 + bb + nub)/mub))/                              &
!!     :  (Gamma((1 + nua)/mua)*Gamma((1 + nub)/mub))    
     
      RETURN
      END
      
! #####################################################################
!
!
! #####################################################################
      Function delvbk(beta,b,nu,mu,k)
!   
!  Purpose: Caluculates (velocity) collection coefficients following Siefert (2006)
!  Equation (92) (self-collection: b collecting b)
!  mass-diameter relationship: D = a*x**(b), where x = particle mass
!  velocity relationship: v = alpha*x**beta
!  general distribution: n(x) = A*x**(nu)*Exp(-lam*x**(mu))
!  where
!      A = mu*N/(Gamma((nu+1)/mu)) *lam**((nu+1)/mu)
!
!      lam = ( Gamma((nu+1)/mu)/Gamma((nu+2)/mu) * xbar )**(-mu)
!
!     where  xbar = L/N  (mass content)/(number concentration) = q*rhoa/N
!
      implicit none
      real delvbk, gamma
      real nu, mu, b, beta
      integer k
      
      delvbk =  &
     &   ((Gamma((1 + nu)/mu)/Gamma((2 + nu)/mu))**(2*beta + k)*  &
     &    Gamma((1 + 2*b + 2*beta + k + nu)/mu))/                 &
     &    Gamma((1 + 2*b + k + nu)/mu)
      RETURN
      END
      
! #####################################################################
!
!
! #####################################################################
! Equation (93) in Seifert and Beheng (2006) ("a" collecting "b")
      Function delvabk(betaa,betab,ba,bb,nua,nub,mua,mub,k)

      implicit none
      real delvabk, gamma
      real nua, mua, ba, k
      real nub, mub, bb
      real betaa, betab
      
      
      delvabk =  &
     &  (2*(Gamma((1 + nua)/mua)/Gamma((2 + nua)/mua))**betaa*   &
     &    Gamma((1 + ba + betaa + nua)/mua)*                     &
     &    (Gamma((1 + nub)/mub)/Gamma((2 + nub)/mub))**betab*    &
     &    Gamma((1 + bb + betab + k + nub)/mub))/                &
     &  (Gamma((1 + ba + nua)/mua)*Gamma((1 + bb + k + nub)/mub))
      
      RETURN
      END

! #####################################################################
!
!
! #####################################################################

     REAL FUNCTION cnudiag(ccw)
     
     implicit none
     
     real :: ccw
     
     ! local
     real :: vard, tmp
     
     real :: vard_nu(9,2)  ! table of variance vs. shape parameter
     real :: vard_ccw(4,3) ! table of Charndrakar et al. 2016

      DATA vard_nu(:,1) / 0.406873, 0.381485, 0.360077, 0.34174, 0.325827, 0.311862, 0.299489, 0.288433, 0.278483 /
      DATA vard_nu(:,2) / -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5 /
     
      DATA vard_ccw(:,1) / 0.40, 0.38, 0.35, 0.28 /
      DATA vard_ccw(:,2) / 21.3e6, 76.9e6, 201.2e6, 564.6e6 /
      DATA vard_ccw(:,3) / -0.274661, -0.193569, -0.0468808, 0.484082 /
      
      IF ( ccw <= vard_ccw(1,2) ) THEN
        vard = 0.4
        cnudiag = -0.274661
      ELSEIF ( ccw >= vard_ccw(4,2) ) THEN
        vard = 0.28
        cnudiag = 0.484082
      ELSE
        ! interpolate
        IF ( ccw <= vard_ccw(2,2) ) THEN
         ! range of 21.3e6 to 76.9e6
          
          cnudiag = vard_ccw(1,3) + (vard_ccw(2,3)-vard_ccw(1,3))*(ccw - vard_ccw(1,2))/( vard_ccw(2,2) - vard_ccw(1,2) )
          
        ELSEIF ( ccw <= vard_ccw(3,2) ) THEN
         ! range of 76.9e6 to 201.2e6

          cnudiag = vard_ccw(2,3) + (vard_ccw(3,3)-vard_ccw(2,3))*(ccw - vard_ccw(2,2))/( vard_ccw(3,2) - vard_ccw(2,2) )
        
        ELSE
         ! range of 201.2e6 to 564.6e6

          cnudiag = vard_ccw(3,3) + (vard_ccw(4,3)-vard_ccw(3,3))*(ccw - vard_ccw(3,2))/( vard_ccw(4,2) - vard_ccw(3,2) )
        
        ENDIF
      
      
      ENDIF
     RETURN
     END FUNCTION cnudiag



