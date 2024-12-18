c #include "sam.def.h"
c
c ##############################################################################
c
      SUBROUTINE setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
     :                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,
     :                 ipconc,ndebug,ngs,nz,kgs,cwnccn,fadvisc,
     :                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,
     :                 itype1,itype2,temcg,infdo)

       USE INDEX_MODULE, only: lt,lc,lr,li,ls,lh,lhl,lip,lir,lgl,lgm,lgh,lf,lv,lhab,lg,
     :                          lnc,lnr,lni,lns,lnh,lngl,lngm,lngh,lnhl,lnf,lnir,lnip

      implicit none

c      include 'sam.index.ion.h'
c#ifdef ICE4
c      include 'sam.index.ice4.h'
c#else
c      include 'sam.index.ion.h'
c#endif
      
      integer ngscnt,ngs,nz
      integer infall    ! whether to calculate number-weighted fall speeds
      
      real xv(ngscnt,lc:lhab)
      real qx(ngscnt,lv:lhab)
      real cx(ngscnt,lc:lhab)
      real vtxbar(ngscnt,lc:lhab,2)
      real xmas(ngscnt,lc:lhab)
      real xdn(ngscnt,lc:lhab)
      real xdia(ngscnt,lc:lhab,3)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real qxmin(lc:lhab)
      real cdx(lc:lhab)
      
      real rho0(ngs),rhovt(ngs),temcg(ngs)
      real cno(lc:lhab)
      
      real cwc1, cwnccn(nz), cimn, cimx
      real cnina(ngs)
      integer kgs(ngs)
      real fadvisc(ngs)
      
      integer ipconc
      integer ndebug
      
      integer itype1,itype2,infdo
      
c Local vars

      real cwc0 ! ,cwc1
      real pi,pii
      real cimasx,cimasn
      real cwmasn,cwmasx,cwradn
      real cwrad
      
      real ccimx

      integer mgs
      
      real ar,br,cs,ds
      real gf4p5, gf4ds, gf4br, ifirst, gf1ds
      real gr
      real rwrad,rwdia
      real gamma
      real mwfac
      integer il

      save gf4p5, gf4ds, gf4br, ifirst, gf1ds
      data ifirst /0/

c
c set values
c
c      cwmasn = 5.23e-13  ! radius of 5.0e-6
c      cwradn = 5.0e-6
c      cwmasx = 5.25e-10  ! radius of 50.0e-6

      ar = 841.99666  
      br = 0.8
      gr = 9.8
c  new values for  cs and ds
      cs = 12.42
      ds = 0.42

      IF ( ifirst .eq. 0 ) THEN
        ifirst = 1
        gf4br = gamma(4.0+br)
        gf4ds = gamma(4.0+ds)
        gf1ds = gamma(1.0+ds)
        gf4p5 = gamma(4.0+0.5)
      ENDIF
      

      pi = 4.0*atan(1.0)
      pii = 1.0/pi
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter

      cwc1 = 6.0/(pi*1000.)
      cwc0 = pii ! 6.0*pii
      mwfac = 6.0**(1./3.)

      
      if (ndebug .gt. 0 ) print*,'Set scale diameter'
c


c
c  cloud water variables
c
      if ( ndebug .gt. 0 )print*,'Set cloud water variables'
      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      
      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN !{
      
      IF ( ipconc .ge. 2 .and. cx(mgs,lc) .gt. 1.0e6 ) THEN !{
        xmas(mgs,lc) =
     >    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
      ELSE
       IF ( ipconc .lt. 2 ) THEN
         cx(mgs,lc) = cwnccn(kgs(mgs))
       ENDIF
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 0.01 ) THEN !{
        xmas(mgs,lc) =
     >     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)),
     :      xdn(mgs,lc)*xvmx(lc) )
        
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)
        
       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 0.01 ) THEN
        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        
       ELSE
        xmas(mgs,lc) = cwmasn
c do not define ccw here! it can feed back to ccn!!!    cx(mgs,lc) = 0.0 ! cwnc(mgs)
       ENDIF !}
      ENDIF !}
c      IF ( ipconc .lt. 2 ) THEN
c        xmas(mgs,lc) =
c     >    min( max(qx(mgs,lc)*rho0(mgs)/cwnc(mgs),cwmasn),cwmasx )
c        cx(mgs,lc) = Max(1.0,qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc))
c      ELSE
c        cwnc(mgs) = an(igs(mgs),jgs,kgs(mgs),lnc)
c        cx(mgs,lc) = cwnc(mgs)
c      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**(1./3.)
      xdia(mgs,lc,2) = xdia(mgs,lc,1)**2
      cwrad = 0.5*xdia(mgs,lc,1)
      vtxbar(mgs,lc,1) =
     >   (2.0*gr*xdn(mgs,lc) *(cwrad**2))
     >  /(9.0*fadvisc(mgs))

      
      ELSE
       xmas(mgs,lc) = cwmasn
       IF ( ipconc .le. 1 ) cx(mgs,lc) = 0.01
       xdia(mgs,lc,1) = 2.*cwradn
       xdia(mgs,lc,2) = 4.*cwradn**2
       vtxbar(mgs,lc,1) = 0.0
       
      ENDIF !} qcw .gt. qxmin(lc)
      
      end do



c
c cloud ice variables
c plates
c
c #ifndef ICE4
      if ( ndebug .gt. 0 )  print*,'Set cip'
      do mgs = 1,ngscnt
      IF ( ipconc .eq. 0 ) THEN
c       cx(mgs,lip) = min(cnit*exp(-temcg(mgs)*bta1),1.e+09)
       cx(mgs,lip) = cnina(mgs)
       IF ( cimn .gt. 1.0 ) THEN
         cx(mgs,lip) = Max(cimn,cx(mgs,lip))
       ENDIF
       IF ( cimx .gt. 1.0 ) THEN
         cx(mgs,lip) = Min(cimx,cx(mgs,lip))
       ENDIF
c erm 3/28/2002
       IF ( itype1 .ge. 1 .or. itype2 .ge. 1 ) THEN 
        cx(mgs,lip) = Max(cx(mgs,lip),qx(mgs,lip)*rho0(mgs)/cimasx)
        cx(mgs,lip) = Min(cx(mgs,lip),qx(mgs,lip)*rho0(mgs)/cimasn)
       ENDIF
c
       cx(mgs,lip) = max(1.0,cx(mgs,lip))
       cx(mgs,lip) = Min(ccimx, cx(mgs,lip))

      ELSEIF ( ipconc .ge. 1 ) THEN
        IF ( qx(mgs,lip) .gt. qxmin(lip) ) THEN
         cx(mgs,lip) = Max(cx(mgs,lip),qx(mgs,lip)*rho0(mgs)/cimasx)
         cx(mgs,lip) = Min(cx(mgs,lip),qx(mgs,lip)*rho0(mgs)/cimasn)
c         cx(mgs,lip) = Max(1.0,cx(mgs,lip))
        ENDIF
      ENDIF

      IF ( qx(mgs,lip) .gt. qxmin(lip) ) THEN
      xmas(mgs,lip) =
     >  min( max(qx(mgs,lip)*rho0(mgs)/cx(mgs,lip),cimasn),cimasx )
      if ( temcg(mgs) .gt. 0 ) then
      xdia(mgs,lip,1) = 0.0
      else
      if ( xmas(mgs,lip) .gt. 0.0 ) THEN ! cimasn ) then
      xdia(mgs,lip,1) = 1.188*(xmas(mgs,lip)**(0.4040))
      end if
      end if
      xdia(mgs,lip,1) = max(xdia(mgs,lip,1), 10.e-6)
      xdia(mgs,lip,1) = min(xdia(mgs,lip,1), 3000.e-6)
c ciplen not used!
c      ciplen(mgs) = 1.25e-3*(xdia(mgs,lip,1)**(0.4740))
      xdn(mgs,lip)  = 900.0
      xdia(mgs,lip,2) = xdia(mgs,lip,1)**2
      vtxbar(mgs,lip,1) = (131.6)*(xdia(mgs,lip,1)**(0.824))
      vtxbar(mgs,lip,1) = vtxbar(mgs,lip,1)*rhovt(mgs)
c      cipat(mgs) = (0.25)*pi*xdia(mgs,lip,2)*cx(mgs,lip)
      ELSE
       xmas(mgs,lip) = 1.e-13
       xdn(mgs,lip)  = 900.0
       xdia(mgs,lip,1) = 1.e-7
       xdia(mgs,lip,2) = (1.e-14)
       vtxbar(mgs,lip,1) = 0.0
c       cipat(mgs) = 0.0
      ENDIF
      
      end do
c #endif
! ICE4 cut
c
c cloud ice variables
c columns
c
      if ( ndebug .gt. 0 ) print*,'Set ci'
      do mgs = 1,ngscnt
      IF ( ipconc .eq. 0 ) THEN
c       cx(mgs,li) = min(cnit*exp(-temcg(mgs)*bta1),1.e+09)
        cx(mgs,li) = cnina(mgs)
       IF ( cimn .gt. 1.0 ) THEN
         cx(mgs,li) = Max(cimn,cx(mgs,li))
       ENDIF
       IF ( cimx .gt. 1.0 ) THEN
         cx(mgs,li) = Min(cimx,cx(mgs,li))
       ENDIF
c erm 3/28/2002
       IF ( itype1 .ge. 1 .or. itype2 .ge. 1 ) THEN
        cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
        cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
       ENDIF
c
       cx(mgs,li) = max(1.0e-20,cx(mgs,li))
c       cx(mgs,li) = Min(ccimx, cx(mgs,li))

      
      ELSEIF ( ipconc .ge. 1 ) THEN
        IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
         cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
         cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
c         cx(mgs,li) = Max(1.0,cx(mgs,li))
        ENDIF
      ENDIF
      
      IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
      xmas(mgs,li) =
     >  min( max(qx(mgs,li)*rho0(mgs)/cx(mgs,li),cimasn),cimasx )
      
      if ( temcg(mgs) .gt. 0.0 ) then
      xdia(mgs,li,1) = 0.0
      else
      if ( xmas(mgs,li) .gt. 0.0 ) THEN ! cimasn ) then
c      xdia(mgs,li,1) = 0.4892*(xmas(mgs,li)**(0.4554))
      xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))
      end if
      end if
      xdia(mgs,li,1) = max(xdia(mgs,li,1), 10.e-6)
      xdia(mgs,li,1) = min(xdia(mgs,li,1), 1000.e-6)
      vtxbar(mgs,li,1) = (4.942e4)*(xdia(mgs,li,1)**(1.4150))
c      xdn(mgs,li)   = min(max(769.8*xdia(mgs,li,1)**(-0.0140),300.0),900.0)
      xdn(mgs,li) = 900.0
      xdia(mgs,li,2) = xdia(mgs,li,1)**2
      vtxbar(mgs,li,1) = vtxbar(mgs,li,1)*rhovt(mgs)
c      ciat(mgs) = (0.25)*pi*xdia(mgs,li,2)*cx(mgs,li)
      ELSE
       xmas(mgs,li) = 1.e-13
       xdn(mgs,li)  = 900.0
       xdia(mgs,li,1) = 1.e-7
       xdia(mgs,li,2) = (1.e-14)
       vtxbar(mgs,li,1) = 0.0
c       cicap(mgs) = 0.0
c       ciat(mgs) = 0.0
      ENDIF
      end do


      do mgs = 1,ngscnt

      IF  ( ipconc .ge. 1 ) THEN
        IF ( qx(mgs,lir) .gt. qxmin(lir) ) THEN
         cx(mgs,lir) = Max(cx(mgs,lir),qx(mgs,lir)*rho0(mgs)/cimasx)
         cx(mgs,lir) = Min(cx(mgs,lir),qx(mgs,lir)*rho0(mgs)/cimasn)
        ENDIF
      ENDIF


      if ( qx(mgs,lir) .gt. qxmin(lir) ) then
      if ( ipconc .eq. 0 ) then
      xdia(mgs,lir,1) =
     >  (qx(mgs,lir)*rho0(mgs)/(pi*xdn(mgs,lir)*cno(lir)))**(0.25) 
      elseif ( ipconc .ge. 1 ) then
c      IF ( cx(mgs,lir) .lt. 10. ) THEN 
c         cx(mgs,lir) = Max(1.0, qx(mgs,lir)*rho0(mgs)/(cimasx+cimasn))
c      ENDIF
      xdia(mgs,lir,1) =
     >  (qx(mgs,lir)*rho0(mgs)
     > /(pi*xdn(mgs,lir)*cx(mgs,lir)))**(0.333333)
      end if
      else
c      xdia(mgs,lir,1) = 1.e-20
      xdia(mgs,lir,1) = 1.e-9
      end if
      xdia(mgs,lir,2) = xdia(mgs,lir,1)**2
      xmas(mgs,lir) = Max(cimasn,xdn(mgs,lir)*(pi/6.)*xdia(mgs,lir,1)**3)
      xmas(mgs,lir) = Min(cimasx, xmas(mgs,lir))
      end do
      
c #ifndef ICE4
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

        xdia(mgs,lf,1) = (xv(mgs,lf)*cwc0)**(1./3.)

      ELSE
      xdia(mgs,lf,1) =
     >  (qx(mgs,lf)*rho0(mgs)/(pi*xdn(mgs,lf)*cno(lf)))**(0.25) 
      end if
      else
      xdia(mgs,lf,1) = 1.e-20
      end if
      xdia(mgs,lf,2) = xdia(mgs,lf,1)**2
c      fwdia3(mgs) = xdia(mgs,lf,2)*xdia(mgs,lf,1)
c      xmas(mgs,lf) = xdn(mgs,lf)*(pi/6.)*fwdia3(mgs)
      end do
c #endif
c
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      
c      IF ( qx(mgs,lr) .gt. 10.0e-3 )
c     :  print*, 'RAIN1: ',igs(mgs),kgs(mgs),qx(mgs,lr)
      
      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-9,cx(mgs,lr)))
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
          xv(mgs,lr) = xvmx(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF

        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,1) = (xmas(mgs,lr)*cwc1)**(1./3.)
        xdia(mgs,lr,3) = xdia(mgs,lr,1)
c        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

c Inverse exponential version:
c        xdia(mgs,lr,1) =
c     >  (qx(mgs,lr)*rho0(mgs)
c     > /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      ELSE
        xdia(mgs,lr,1) =
     >  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.)
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
c        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if
      xdia(mgs,lr,2) = xdia(mgs,lr,1)**2
c      xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
      end do
c
      do mgs = 1,ngscnt 
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      if ( ipconc .ge. 4 ) then ! 

        xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
c      parameter( xvmn(lr)=2.8866e-13, xvmx(lr)=4.1887e-9 )  ! mks
        xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)

        IF ( xv(mgs,ls) .lt. xvmn(ls) .or. xv(mgs,ls) .gt. xvmx(ls) ) THEN
          xv(mgs,ls) = Min( xvmx(ls), Max( xvmn(ls),xv(mgs,ls) ) )
          xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
        ENDIF

        xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0)**(1./3.)

      ELSE
      xdia(mgs,ls,1) = 
     >  (qx(mgs,ls)*rho0(mgs)/(pi*xdn(mgs,ls)*cno(ls)))**(0.25) 
      end if
      else
      xdia(mgs,ls,1) = 1.e-9
      cx(mgs,ls) = 0.0
      end if
      xdia(mgs,ls,2) = xdia(mgs,ls,1)**2
c      swdia3(mgs) = xdia(mgs,ls,2)*xdia(mgs,ls,1)
c      xmas(mgs,ls) = xdn(mgs,ls)*(pi/6.)*swdia3(mgs)
      end do
c
c #ifndef ICE4
      do mgs = 1,ngscnt 
      if ( qx(mgs,lgl) .gt. qxmin(lgl) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lgl) = rho0(mgs)*qx(mgs,lgl)/(xdn(mgs,lgl)*Max(1.0e-9,cx(mgs,lgl)))
        xmas(mgs,lgl) = xv(mgs,lgl)*xdn(mgs,lgl)

        IF ( xv(mgs,lgl) .lt. xvmn(lgl) .or. xv(mgs,lgl) .gt. xvmx(lgl) ) THEN
          xv(mgs,lgl) = Min( xvmx(lgl), Max( xvmn(lgl),xv(mgs,lgl) ) )
          xmas(mgs,lgl) = xv(mgs,lgl)*xdn(mgs,lgl)
          cx(mgs,lgl) = rho0(mgs)*qx(mgs,lgl)/(xmas(mgs,lgl))
        ENDIF

        xdia(mgs,lgl,1) = (xv(mgs,lgl)*cwc0)**(1./3.)


      ELSE
      xdia(mgs,lgl,1) = 
     >  (qx(mgs,lgl)*rho0(mgs)/(pi*xdn(mgs,lgl)*cno(lgl)))**(0.25) 
      end if
      else
      xdia(mgs,lgl,1) = 1.e-9
      end if
      xdia(mgs,lgl,2) = xdia(mgs,lgl,1)**2
c      gldia3(mgs) = xdia(mgs,lgl,2)*xdia(mgs,lgl,1)
c      xmas(mgs,lgl) = xdn(mgs,lgl)*(pi/6.)*gldia3(mgs)
      end do
c #endif
c
      do mgs = 1,ngscnt 
      if ( qx(mgs,lgm) .gt. qxmin(lgm) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lgm) = rho0(mgs)*qx(mgs,lgm)/(xdn(mgs,lgm)*Max(1.0e-9,cx(mgs,lgm)))
        xmas(mgs,lgm) = xv(mgs,lgm)*xdn(mgs,lgm)

        IF ( xv(mgs,lgm) .lt. xvmn(lgm) .or. xv(mgs,lgm) .gt. xvmx(lgm) ) THEN
          xv(mgs,lgm) = Min( xvmx(lgm), Max( xvmn(lgm),xv(mgs,lgm) ) )
          xmas(mgs,lgm) = xv(mgs,lgm)*xdn(mgs,lgm)
          cx(mgs,lgm) = rho0(mgs)*qx(mgs,lgm)/(xmas(mgs,lgm))
        ENDIF

        xdia(mgs,lgm,1) = (xv(mgs,lgm)*cwc0)**(1./3.)

      ELSE
       xdia(mgs,lgm,1) = 
     >  (qx(mgs,lgm)*rho0(mgs)/(pi*xdn(mgs,lgm)*cno(lgm)))**(0.25) 
      end if
      else
      xdia(mgs,lgm,1) = 1.e-9
      end if
      xdia(mgs,lgm,2) = xdia(mgs,lgm,1)**2
c      gmdia3(mgs) = xdia(mgs,lgm,2)*xdia(mgs,lgm,1)
c      xmas(mgs,lgm) = xdn(mgs,lgm)*(pi/6.)*gmdia3(mgs)
      end do
c
c #ifndef ICE4
      do mgs = 1,ngscnt  
      if ( qx(mgs,lgh) .gt. qxmin(lgh) ) then
      if ( ipconc .ge. 5 ) then

        cx(mgs,lgh) = Max(1.0e-9,cx(mgs,lgh))
        xv(mgs,lgh) = rho0(mgs)*qx(mgs,lgh)/(xdn(mgs,lgh)*cx(mgs,lgh))
        xmas(mgs,lgh) = xv(mgs,lgh)*xdn(mgs,lgh)

        IF ( xv(mgs,lgh) .lt. xvmn(lgh) .or. xv(mgs,lgh) .gt. xvmx(lgh) ) THEN
          xv(mgs,lgh) = Min( xvmx(lgh), Max( xvmn(lgh),xv(mgs,lgh) ) )
          xmas(mgs,lgh) = xv(mgs,lgh)*xdn(mgs,lgh)
          cx(mgs,lgh) = rho0(mgs)*qx(mgs,lgh)/(xmas(mgs,lgh))
        ENDIF

        xdia(mgs,lgh,1) = (xv(mgs,lgh)*cwc0)**(1./3.)

      ELSE
      xdia(mgs,lgh,1) = 
     >  (qx(mgs,lgh)*rho0(mgs)/(pi*xdn(mgs,lgh)*cno(lgh)))**(0.25) 
      end if
      else
      xdia(mgs,lgh,1) = 1.e-9
      end if
      xdia(mgs,lgh,2) = xdia(mgs,lgh,1)**2
c      ghdia3(mgs) = xdia(mgs,lgh,2)*xdia(mgs,lgh,1)
c      xmas(mgs,lgh) = xdn(mgs,lgh)*(pi/6.)*ghdia3(mgs)
      end do
c #endif
c
      do mgs = 1,ngscnt 
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*Max(1.0e-9,cx(mgs,lh)))
        xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)

        IF ( xv(mgs,lh) .lt. xvmn(lh) .or. xv(mgs,lh) .gt. xvmx(lh) ) THEN
          xv(mgs,lh) = Min( xvmx(lh), Max( xvmn(lh),xv(mgs,lh) ) )
          xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
          cx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xmas(mgs,lh))
        ENDIF

        xdia(mgs,lh,1) = (xv(mgs,lh)*cwc0)**(1./3.)

      ELSE
      xdia(mgs,lh,1) = 
     >  (qx(mgs,lh)*rho0(mgs)/(pi*xdn(mgs,lh)*cno(lh)))**(0.25) 
      end if
      else
      xdia(mgs,lh,1) = 1.e-9
      end if
      xdia(mgs,lh,2) = xdia(mgs,lh,1)**2
c      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
c      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
c      
c
c #ifndef ICE4
      do mgs = 1,ngscnt 
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*Max(1.0e-9,cx(mgs,lhl)))
        xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)

        IF ( xv(mgs,lhl) .lt. xvmn(lhl) .or. xv(mgs,lhl) .gt. xvmx(lhl) ) THEN
          xv(mgs,lhl) = Min( xvmx(lhl), Max( xvmn(lhl),xv(mgs,lhl) ) )
          xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
          cx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xmas(mgs,lhl))
        ENDIF

        xdia(mgs,lhl,1) = (xv(mgs,lhl)*cwc0)**(1./3.)

      ELSE
      xdia(mgs,lhl,1) = 
     >  (qx(mgs,lhl)*rho0(mgs)/(pi*xdn(mgs,lhl)*cno(lhl)))**(0.25) 
      end if
      else
      xdia(mgs,lhl,1) = 1.e-9
      end if
      xdia(mgs,lhl,2) = xdia(mgs,lhl,1)**2
c      hldia3(mgs) = xdia(mgs,lhl,2)*xdia(mgs,lhl,1)
c      xmas(mgs,lhl) = xdn(mgs,lhl)*(pi/6.)*hldia3(mgs)
      end do
c #endif
c
c  Set terminal velocities...
c    also set drag coefficients (moved to start of subroutine)
c
c      cdx(lr) = 0.60
c      cdx(lh) = 0.45
c      cdx(lhl) = 0.45
c      cdx(lf) = 0.45
c      cdx(lgh) = 0.60
c      cdx(lgm) = 0.80
c      cdx(lgl) = 0.80
c      cdx(lir) = 2.00
c
      if ( ndebug .gt. 0 ) print*,'Set terminal velocities'
c
      do mgs = 1,ngscnt
      if ( qx(mgs,lir) .gt. qxmin(lir) ) then
      vtxbar(mgs,lir,1) = (cs*gf4ds/6.0)*(xdia(mgs,lir,1)**ds)*rhovt(mgs)
      else
      vtxbar(mgs,lir,1) = 0.0
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set rimed ice vt'
c
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      IF ( ipconc .lt. 3 ) THEN
        vtxbar(mgs,lr,1) = (ar*gf4br/6.0)*(xdia(mgs,lr,1)**br)*rhovt(mgs)
      ELSE
        rwdia = Min( xdia(mgs,lr,1), 8.0e-3 )
        
         vtxbar(mgs,lr,1) = rhovt(mgs)*6.0*pii*( 0.04771 + 3788.0*rwdia - 
     :        1.105e6*rwdia**2 + 1.412e8*rwdia**3 - 6.527e9*rwdia**4)
        
        IF ( infdo .ge. 1 ) THEN
         vtxbar(mgs,lr,2) = (0.09112 + 2714.0*rwdia - 4.872e5*rwdia**2 + 
     :            4.495e7*rwdia**3 - 1.626e9*rwdia**4)*rhovt(mgs)
        ENDIF

c        IF ( rwrad*mwfac .gt. 6.0e-4  ) THEN
c          vtxbar(mgs,lr,1) = 20.1*Sqrt(100.*rwrad*mwfac)*rhovt(mgs)
c        ELSE
c          vtxbar(mgs,lr,1) = 80.0e2*rwrad*rhovt(mgs)*mwfac
c        ENDIF
c        IF ( rwrad .gt. 6.0e-4  ) THEN
c          vtxbar(mgs,lr,2) = 20.1*Sqrt(100.*rwrad)*rhovt(mgs)
c        ELSE
c          vtxbar(mgs,lr,2) = 80.0e2*rwrad*rhovt(mgs)
c        ENDIF
      ENDIF
      else
      vtxbar(mgs,lr,1) = 0.0
      vtxbar(mgs,lr,2) = 0.0
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set rain vt'
c
      do mgs = 1,ngscnt
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      vtxbar(mgs,ls,1) = (cs*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
      else
      vtxbar(mgs,ls,1) = 0.0
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set snow vt'
c
c #ifndef ICE4
      do mgs = 1,ngscnt
      if ( qx(mgs,lgl) .gt. qxmin(lgl) ) then
      vtxbar(mgs,lgl,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lgl)*xdia(mgs,lgl,1)*4.0*gr) /
     <    (3.0*cdx(lgl)*rho0(mgs)) )
      else
      vtxbar(mgs,lgl,1) = 0.0
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set low graupel vt'
c #endif
c
      do mgs = 1,ngscnt
      vtxbar(mgs,lgm,1) = 0.0
      if ( qx(mgs,lgm) .gt. qxmin(lgm) ) then
      vtxbar(mgs,lgm,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lgm)*xdia(mgs,lgm,1)*4.0*gr) /
     <    (3.0*cdx(lgm)*rho0(mgs)) )
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set medium graupel vt'
c
c #ifndef ICE4
      do mgs = 1,ngscnt
      vtxbar(mgs,lgh,1) = 0.0
      if ( qx(mgs,lgh) .gt. qxmin(lgh) ) then
      vtxbar(mgs,lgh,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lgh)*xdia(mgs,lgh,1)*4.0*gr) /
     <    (3.0*cdx(lgh)*rho0(mgs)) )
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set heavy graupel vt'
c
      do mgs = 1,ngscnt
      vtxbar(mgs,lf,1) = 0.0
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
      vtxbar(mgs,lf,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lf)*xdia(mgs,lf,1)*4.0*gr) /
     <    (3.0*cdx(lf)*rho0(mgs)) )
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set frozen drop vt'
c #endif
c
      do mgs = 1,ngscnt
      vtxbar(mgs,lh,1) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      vtxbar(mgs,lh,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lh)*xdia(mgs,lh,1)*4.0*gr) /
     <    (3.0*cdx(lh)*rho0(mgs)) )
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set hail vt'
c
c #ifndef ICE4
      do mgs = 1,ngscnt
      vtxbar(mgs,lhl,1) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      vtxbar(mgs,lhl,1) = (gf4p5/6.0)*
     >  Sqrt( (xdn(mgs,lhl)*xdia(mgs,lhl,1)*4.0*gr) /
     <    (3.0*cdx(lhl)*rho0(mgs)) )
      end if
      end do
      if ( ndebug .gt. 0 ) print*,'Set hail vt'
c #endif
      
      
      
      IF ( infdo .ge. 1 ) THEN

      
c      DO il = lc,lhab
c      IF ( il .ne. lr ) THEN
        DO mgs = 1,ngscnt
          vtxbar(mgs,lc,2) = vtxbar(mgs,lc,1)
          vtxbar(mgs,li,2) = vtxbar(mgs,li,1)
c #ifndef ICE4
          vtxbar(mgs,lip,2) = vtxbar(mgs,lip,1)
c #endif
          vtxbar(mgs,lir,2) = vtxbar(mgs,lir,1)
          vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        ENDDO

        DO mgs = 1,ngscnt
         IF ( qx(mgs,ls) .gt. qxmin(ls) ) THEN
          vtxbar(mgs,ls,2) = (cs*gf1ds)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         ENDIF
        ENDDO

        
!#ifndef ICE4
        DO il = lg,lhab
!#else
!        DO il = lgm,lhab
!#endif
          IF ( il .ne. lip ) THEN
            DO mgs = 1,ngscnt
             IF ( qx(mgs,il) .gt. qxmin(il) ) THEN
              vtxbar(mgs,il,2) = 
     >            Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) /
     <              (3.0*cdx(il)*rho0(mgs)) )
        
             ENDIF
           ENDDO
          ENDIF
        ENDDO
        
c      ENDIF
c      ENDDO
      
      ENDIF
      
      

c############ SETVT ############################

      RETURN
      END
c--------------------------------------------------------------------------

