
      subroutine raintable(mqcw,mqxw,mrho,qfrwax,qfxarw) 
      
      implicit none
      
      integer mqcw,mqxw,mrho
      real, dimension(mqcw,mqxw,mrho) :: qfrwax,qfxarw
      
      real sqcw,sqfw, qcw, qfw
      real fqcw, fqfw
      real temg, p0s, rho0
      
      integer lt
      integer lqcw,lqfw

! ####################################################################

      sqcw = 1.e-13
      sqfw = 1.e-13
      
      fqcw = 0.0005
      fqfw = 0.0005
      
      temg = 273.16
      p0s = 67383.0
      rho0 = p0s/(287.04*temg)
      print*, 'Rain lookup table: ',temg,p0s,rho0
      
      DO lqcw = 1,mqcw
      DO lqfw = 1,mqxw

      qcw = sqcw + (lqcw-1)*fqcw
      qfw = sqfw + (lqfw-1)*fqfw


      
      ENDDO
      ENDDO
      
      RETURN
      END
      

!
!2345678901234567890123456789012345678901234567890123456789012345678912
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  STRAKA'S ATMOSPHERIC MODEL  (SAM)
!    Designed by Jerry M. Straka
!
!  SUBROUTINE RAINMCONV
!
!  A subroutine version of sam.qfxacr.f adapted by ERM to calculate
!  conversions of graupel and frozen drops on the fly.  Main reason 
!  is to eliminate the need for the lookup tables.
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine rainconv(ipconc,qrw,qfw,crw,cfw,                  &
                 cdfw,cnof,fwdia,rwmas,rwdia,                      &
                 cnor,xvr,rhovt,                                   &
                 rfrh,rfrr,intyp,temg,rho0,cd,cs,ds,crfrh,crfrr)

!
      implicit none
      
      integer ndebug
      parameter ( ndebug = 0 )
      
      integer ipconc  !  which categories are two-moment

      real qrw,qfw,crw,cfw
      real crwt, cfwt
      real rwdia,rwmas
      real fwdia
      real rtfw    !  riming time
      real cnof,cnor
      real cdfw
      integer intyp
      real rnu,xvr
      parameter ( rnu = -0.8 )

!      character*40 filouth, filoutg
!      parameter(ngs=100,nh=100,ng=100)
!      parameter(mqcw=21,mqfw=21,mt=1,mr=5)
      integer mr,nh,nnr
      parameter(nh=11,nnr=nh)
      parameter(mr=5)
      
      
      real rfrh(mr)     !  fraction of mass to each category
      real rfrr(mr)
      real crfrh(mr), crfrr(mr)
      real rfrh2(mr), rfrh3(mr)  ! test fractions
      
      real rhop(nh,nnr), rhog(nh)
      real rhopr(nh,nnr)
      real rn(nnr) !,rd(nnr),rm(nnr)
      real rq(nnr),vtr(nnr) !,rdrd(nnr)
      real hn(nh),hd(nh),hm(nh) ! ,dmdt(nh)
!      real dmhdt(nh,nnr),dmrdt(nh,nnr)
      real dmhdt0(nh,nnr,mr),dmrdt0(nh,nnr,mr),hmrmainv(nh,nnr,mr)
      save dmhdt0, dmrdt0, hmrmainv
      real hq(nh),hqp(nh),hdhd(nh),vth(nh) ! ,dqdt(nh)
      real dqhdt(nh,nnr),dqrdt(nh,nnr)
      real dchdt(nh,nnr),dcrdt(nh,nnr)
      real rhonew(nh)
      real rwdn
      real fwmas
      real fwdn,fslp,wslp
      real rho0,temg,rho0inv,rhovt
      real rimf,qcmxd,qfacw,efw,rhobar,qtotp
      real qhtot(mr),qrtot(mr),qhtot2(mr),qhtot3(mr)
      real chtot(mr),crtot(mr)
      real chtot2(mr)
      real qhtotp,qrtotp
      real qhtotp2,qhtotp3
      real chtotp,crtotp
      real chtotp2

      real rhos(mr),rhoe(mr),qtot(mr)
      real rhoss(mr),rhoes(mr)
      data rhoss /000,150,350,600,770/
      data rhoes /150,350,600,770,1001/
      real rhosg(5),rhoeg(5)
      data rhosg /000,350,600,770,9000/
      data rhoeg /350,600,770,1001,9999/

       integer nr, n
       real hrho
       real gr, ds, cs, pi, ar, br
       real cd(5)
       real fac,x
       
       integer ifirst
       data ifirst/1/
       save ifirst
       
       integer nn,m
       parameter (nn=nh)
       real hma(nn,mr),hda(nn,mr),hdhda(nn,mr),hrhoa(mr)
       real hvai(nn,mr)  ! inverse of particle volume
       save hma,hda,hdhda,hrhoa,hvai
       real vtha(nn,mr)
       data hrhoa/100.,300.,500.,700.,800./
       
       real vtra(nn),rda(nn),rma(nn),rdrda(nn),rva(nn),rvna(nn)
       save vtra,rda,rma,rdrda,rva,rvna,vtha
       
       
       real nt
       parameter ( nt = 1 )
       
       parameter ( rwdn = 1000. )
       
       integer l,lr,i
       
       integer idebug
       parameter ( idebug = 0 )
       
       real, parameter :: deltat = 30.
       

!      real rfrh(mqcw,mqfw,mt,mr)
!      real rfrg(mqcw,mqfw,mt,mr)
!
!
!     read
!
!     constants
!
       integer ivt
       
       real hmmin,hjo, hrmin, rjo
       parameter ( hjo = 7.5*nh/(41.) )
       parameter ( rjo = 7.5*nnr/(41.) )
       
       real sumh, sumg
       real vrnu
       
       real gfrnu,xvri
       parameter ( gfrnu = 4.5908437 )  ! Gamma(1 + rnu)
       
       real rho00
       parameter ( rho00 = 1.225 )

       parameter ( ar = 841.99666 )
       parameter ( br = 0.8 )
       parameter ( gr = 9.8 )

! ####################################################################
       n = nh
      hmmin = 4.7e-10
      hrmin = 2.192e-11  ! give top diameter of 5mm
      pi = 4.*atan(1.)
      

      rho0inv = 1./rho0
      
      xvri = 1./xvr
      vrnu = (xvri)**rnu/(xvr * gfrnu )
!
       IF ( ifirst .eq. 1 ) THEN
        ifirst = 0
       
        DO m=1,mr
          hrho = hrhoa(m)
        DO l = 1,nh
         hma(l,m) = hmmin*exp(3.0*(l-1)/hjo)           ! mass of particle type m in bin l
         hvai(l,m) = hrho/hma(l,m)                     ! particle volume
         hda(l,m) = (6.*hma(l,m)/(pi*hrho))**(1./3.)   ! diameter of particle in bin l
         hdhda(l,m) = hda(l,m)/hjo                     ! width of bin (length) i.e., d(diameter)
        ENDDO
        ENDDO
                
        DO l = 1,nh
         vtha(l,1) = cs*(hda(l,1)**ds) ! vth(l) = vtha(l,1)*Sqrt(1.225/rho0)
        ENDDO
        
        IF ( ipconc .ge. 3 ) THEN
        
        DO l = 1,nnr
         
         rma(l)  = hrmin*exp(3.0*(l-1)/rjo)
         rva(l)  = 1.e-3*rma(l)
         rvna(l) = (rnu + 1)**(rnu + 1) * rva(l)**rnu
         rda(l)  = (6.*rma(l)/(pi*1000.))**(1./3.)
         x = rda(l)
         vtra(l) = 0.09112 + 4376.81*x - 697234.5*x**2 +    &
                     4.4946e7*x**3 - 9.8346e8*x**4
!          IF ( 0.5*rda(l) .gt. 6.0e-4  ) THEN
!            vtra(l) = 20.1*Sqrt(100.*0.5*rda(l))
!          ELSE
!            vtra(l) = 80.0e2*0.5*rda(l)
!          ENDIF
!         vtra(l) = ar*rda(l)**br
         rdrda(l) = 3.0*rva(l)/rjo
        ENDDO
        
        ELSE ! single-moment rain for lookup tables
        
        DO l = 1,nnr
         
         rma(l)  = hrmin*exp(3.0*(l-1)/rjo)
         
         rva(l)  = 1.e-3*rma(l)
         
         rvna(l) = 0.0 ! (rnu + 1)**(rnu + 1) * rva(l)**rnu
         rda(l)  = (6.*rma(l)/(pi*1000.))**(1./3.) 
         x = rda(l)
         vtra(l) = 0.09112 + 4376.81*x - 697234.5*x**2 +    &
                     4.4946e7*x**3 - 9.8346e8*x**4
         rdrda(l) = rda(l)/rjo
        ENDDO
        
        ENDIF ! ipconc

        DO m = 2,mr
        DO l = 1,nh
          vtha(l,m) = Sqrt(hrhoa(m)*(hda(l,m)*4.0*gr) /       &
                  (3.0*cd(m)*rho00) ) ! *Sqrt(hrho/rho0)
        ENDDO
        ENDDO

      DO m = 1,mr
      do l = 1,nh
      do lr = 1,nnr
       hd(l) = hda(l,m)  
       dmhdt0(l,lr,m) = .25*pi*((hd(l)+rda(lr))**2)*      & ! collision kernel for graupel 
                             max((vtha(l,m)-vtra(lr)),0.)   ! is zero if rain particle is faster
       dmrdt0(l,lr,m) = .25*pi*((hd(l)+rda(lr))**2)*      & ! collision kernel for rain 
                             max((vtra(lr)-vtha(l,m)),0.)   ! is zero if graupel particle is faster
      end do
      end do

      do l = 1,nh
      do lr = 1,nnr
      hmrmainv(l,lr,m) = 1./(hma(l,m)+rma(lr))    ! inverse of combined mass of graupel & rain particle
      end do
      end do

      ENDDO
       
      IF ( ndebug .ge. 3 ) THEN
      open(unit=85,file='rainconv.txt',status='unknown',form='formatted')

      DO m = 1,mr
      write(85,*) 'm = ',m
      do l = 1,nh
      write(85,*) 'l = ',l
      write(85,'(11(1x,e9.4))') (dmhdt0(l,lr,m), lr=1,nnr)
      write(85,'(11(1x,e9.4))') (dmrdt0(l,lr,m), lr=1,nnr)
!      do lr = 1,n
       
!      end do
      end do
      ENDDO
      
      close(85)
      stop

      ENDIF

      ENDIF ! ( ifirst )

       
       hrho = hrhoa(intyp)
      
       IF ( intyp .ne. 1 ) THEN
         nr = 4
         ivt = 1
         fac = rhovt ! Sqrt(hrho*rho0inv)
         DO i = 1,nr
           rhos(i) = rhosg(i)
           rhoe(i) = rhoeg(i)
         ENDDO
       ELSE !  ( intyp .eq. 1 ) THEN  ! SA
!         hrho = 100.
         nr = 5
         ivt = 2
         fac = rhovt ! Sqrt(1.225*rho0inv)
         DO i = 1,nr
           rhos(i) = rhoss(i)
           rhoe(i) = rhoes(i)
         ENDDO
       ENDIF

!      do lqcw = 1,nqcw
!      do lqfw = 1,nqfw
!
      if ( idebug .gt. 0 ) then
      write(*,*)
      write(*,*)'l,hm(g),hd(cm)'
      end if
      if ( idebug .gt. 0 ) then
      write(*,*) '=1 to write distribution d, m'
      end if


      fwdn = hrho
      qcmxd = 1000000.

!
!  compute slope
!
!      if ( qfw .gt. qfmin ) then
!      fslp = ( (rho0*qfw)/(pi*cnof*fwdn) )**(-0.25)
      fslp = 1./fwdia
!
!  write slope if debug
!
      if ( idebug .gt. 0 ) then
!      write(*,*) '1/fslp',1./fslp
      end if

      do l = 1,nh
!      rq(l) = 0.0
!      hq(l) = 0.0
!      rn(l) = 0.0
!      hn(l) = 0.0
! hm is 'graupel' particle mass
      hm(l) =  hma(l,intyp) ! hmmin*exp(3.0*(l-1)/hjo)
! hg is rain particle mass
!      rm(l) = hmmin*exp(3.0*(l-1)/hjo)
      hd(l) = hda(l,intyp) ! (6.*hm(l)/(pi*hrho))**.33333333333
!      rd(l) = (6.*rma(l)/(pi*1000))**.33333333333
      hdhd(l) = hdhda(l,intyp) ! hd(l)/hjo
!      rdrd(l) = rdrda(l) ! rda(l)/hjo
! vth is 'graupel' fall speed
       vth(l) = vtha(l,intyp)*fac

!      hn(l) = cnof*exp(-hd(l)*fslp)*hdhd(l)
      hn(l) = (cfw*fslp)*exp(-hd(l)*fslp)*hdhd(l)
      hq(l) = (hn(l)*hm(l)*rho0inv)
!      hqtot =  hqtot + hq(l)
      end do 


      IF ( ipconc .ge. 3 ) THEN
      do l = 1,nnr
! vtr is rain fall speed
      vtr(l) = vtra(l)*rhovt ! ar*(rda(l)**br)*(1.225/rho0)**.5
!      hn(l) = cnof*exp(-hd(l)*fslp)*hdhd(l)

!      rn(l) = cnor*exp(-rda(l)*wslp)*rdrda(l)
      rn(l) = crw*rvna(l)*vrnu* Exp( -(rnu+1)*rva(l)*xvri )*rdrda(l)
      rq(l) = (rn(l)*rma(l)*rho0inv)
!      rqtot =  rqtot + rq(l)
      end do 
      
      ELSE ! single-moment rain
     
      IF ( qrw .gt. 0. ) THEN
      wslp = ( (rho0*qrw)/(pi*cnor*1000.) )**(-0.25)
      do l = 1,nnr
! vtr is rain fall speed
      vtr(l) = vtra(l)*rhovt ! ar*(rda(l)**br)*(1.225/rho0)**.5
      rn(l) = cnor*exp(-rda(l)*wslp)*rdrda(l)
      rq(l) = (rn(l)*rma(l)*rho0inv)
!      rqtot =  rqtot + rq(l)
      end do 
      ELSE
       rn(:) = 0.0
       rq(:) = 0.0
      ENDIF
      
      ENDIF ! ipconc
!
!  compute mass growth rate of discretized solution
!
!      hqtot = 0.0
!      rqtot = 0.0
      efw = 1.
!
!
!      IF ( ndebug .ge. 1 ) THEN
!       write(6,*) 'rhovt = ', rho
!      ENDIF
      dqhdt(:,:) = 0.0
      dchdt(:,:) = 0.0
      dqrdt(:,:) = 0.0
      dcrdt(:,:) = 0.0
      do l = 1,nh
      rhog(l) = hm(l)
      do lr = 1,nnr
      
      IF ( rda(lr) .gt. 200.e-6 ) THEN
      
!      dmhdt(l,lr) = dmhdt0(l,lr,intyp)*efw*rq(lr)
      
!      dmrdt(l,lr) = dmrdt0(l,lr,intyp)*efw*hq(l)

! collection by 1 graupel particle
      dqhdt(l,lr) = 1.0*dmhdt0(l,lr,intyp)*efw*rq(lr)*rhovt

      dchdt(l,lr) = 1.0*dmhdt0(l,lr,intyp)*efw*rn(lr)*rhovt

! changed rn(l) to rn(lr)
      dqrdt(l,lr) = 1.0*dmrdt0(l,lr,intyp)*efw*hq(l)*rhovt

      dcrdt(l,lr) = 1.0*dmrdt0(l,lr,intyp)*efw*hn(l)*rhovt
      
      ENDIF

      IF ( dqhdt(l,lr) .gt. 0.0 ) THEN
        rhop(l,lr) = (hm(l)+dchdt(l,lr)*rma(lr))*hvai(l,intyp) !hmrmainv(l,lr,intyp)
!        rhop(l,lr) = (hm(l)*fwdn+dchdt(l,lr)*rma(lr)*rwdn)*hmrmainv(l,lr,intyp)

      dqhdt(l,lr) = hn(l)*dqhdt(l,lr)

      dchdt(l,lr) = hn(l)*dchdt(l,lr)

      ELSEIF ( dqrdt(l,lr) .gt. 0.0 ) THEN
        rhopr(l,lr) = (hm(l)*fwdn*dcrdt(l,lr)+rma(lr)*rwdn)*hmrmainv(l,lr,intyp)
      dqrdt(l,lr) = rn(lr)*dqrdt(l,lr)

      dcrdt(l,lr) = rn(lr)*dcrdt(l,lr)
        
      ELSE
        rhop(l,lr) = fwdn
      ENDIF
!      rhop(l,lr) = (hm(l)*fwdn+rma(lr)*rwdn)*hmrmainv(l,lr,intyp)
      
! assume water drop soaks into graupel, add mass of collected droplets
      rhog(l) = rhog(l) + dchdt(l,lr)*rma(lr)*deltat

      IF ( ndebug .ge. 2  ) THEN
       write(6,*) 'rainconv,l,lr ',l,lr,dqhdt(l,lr), dqrdt(l,lr),       &
         dchdt(l,lr), dcrdt(l,lr), rhop(l,lr)
!       write(6,*) hmrmainv(l,lr,intyp)
      ENDIF
      end do
      rhog(l) = rhog(l)*hvai(l,intyp)
      end do 

      IF ( ndebug .ge. 1 .and. qrw .gt. 5.e-4 .and. qfw .gt. 5.e-4 ) THEN

      crwt = 0.0
      DO l=1,nnr
        crwt = crwt + rn(l)
      ENDDO
      cfwt = 0.0
      DO l=1,nh
        cfwt = cfwt + hn(l)
      ENDDO
      
      write(6,*) 'intyp,qrw,qfw,crw,cfw = ', intyp, qrw*1.e3, qfw*1.e3, crw,cfw,crwt,cfwt
      write(6,'(11(1x,e9.4))') (hn(l), l=1,nh)
      write(6,'(11(1x,f9.4))') (hd(l)*1000., l=1,nh)
      write(6,*) 
      write(6,'(11(1x,e9.4))') (rn(l), l=1,nnr)
      write(6,'(11(1x,f9.4))') (rda(l)*1000., l=1,nnr)
      write(6,'(11(1x,e9.4))') (rn(l)*rma(l)*1000., l=1,nnr)
      write(6,*) 'density:'
      write(6,'(11(1x,e9.4))') (rhog(l),l=1,nh)
      write(6,*) 'mass in bin:'
      write(6,'(11(1x,e9.4))') (hn(l)*hm(l),l=1,nh)
!      do l = 1,nh
!      write(6,*) 'l = ',l
!      write(6,'(11(1x,e9.4))') (rhop(l,lr), lr=1,nnr)
!      write(6,'(11(1x,e9.4))') (rhopr(l,lr), lr=1,nnr)
!      write(6,'(11(1x,e9.4))') (dqhdt(l,lr), lr=1,nnr)
!      write(6,'(11(1x,e9.4))') (dchdt(l,lr), lr=1,nnr)
!      write(6,'(11(1x,e9.4))') (dqrdt(l,lr), lr=1,nnr)
!      write(6,'(11(1x,e9.4))') (dcrdt(l,lr), lr=1,nnr)
!!      do lr = 1,n
!       
!!      end do
!      end do
      
      ENDIF
!
!  compute rhop
!
      do l = 1,n
      do lr = 1,n
      end do
      end do

!
!
!  find fractions of q in each density category
!
      do i = 1,nr
      qhtot(i) = 0.
      qhtot2(i) = 0.0
      qhtot3(i) = 0.0
      qrtot(i) = 0.
      chtot(i) = 0.
      chtot2(i) = 0.
      crtot(i) = 0.
      do l = 1,nh
! test sum
      IF ( rhog(l).gt.rhos(i).and.rhog(l).lt.rhoe(i)) THEN
        DO lr = 1,nnr
          qhtot2(i) = qhtot2(i) + dqhdt(l,lr)
        ENDDO
        qhtot3(i) = qhtot3(i) + hm(l)*hn(l)
        chtot2(i) = chtot2(i) + hn(l)
      ENDIF

      do lr = 1,nnr


      if ( rhop(l,lr).gt.rhos(i).and.rhop(l,lr).lt.rhoe(i)) then
      qhtot(i) = qhtot(i) + dqhdt(l,lr)
      chtot(i) = chtot(i) + dchdt(l,lr)
      end if
      if ( rhopr(l,lr).gt.rhos(i).and.rhopr(l,lr).lt.rhoe(i)) then
      qrtot(i) = qrtot(i) + dqrdt(l,lr)
      crtot(i) = crtot(i) + dcrdt(l,lr)
      end if
      end do
      end do
      end do
!
      qhtotp  = 0.0
      qrtotp  = 0.0
      qhtotp2 = 0.0
      qhtotp3 = 0.0
      chtotp2 = 0.0
      chtotp = 0.0
      crtotp = 0.0
      do i = 1,nr
      qhtotp = qhtotp + qhtot(i) 
      qhtotp2 = qhtotp2 + qhtot2(i) 
      qhtotp3 = qhtotp3 + qhtot3(i) 
      qrtotp = qrtotp + qrtot(i) 
      chtotp = chtotp + chtot(i) 
      chtotp2 = chtotp2 + chtot2(i) 
      crtotp = crtotp + crtot(i) 
      end do
!
      sumh = 0
      sumg = 0
      do i = 1,nr
      if ( qhtotp .gt. 0 ) then
      rfrh(i) = qhtot(i)/qhtotp
      else
      rfrh(i) = 0.0
      crfrh(i) = 0.0
      end if

      if ( qhtotp2 .gt. 0 ) then
      rfrh2(i) = qhtot2(i)/qhtotp2
      else
      rfrh2(i) = 0.0
      end if

      if ( qhtotp3 .gt. 0 ) then
      rfrh3(i) = qhtot3(i)/qhtotp3
      crfrh(i) = rfrh3(i) ! chtot2(i)/chtotp2
      else
      rfrh3(i) = 0.0
      end if

      if ( qrtotp .gt. 0 ) then
      rfrr(i) = qrtot(i)/qrtotp
      crfrr(i) = crtot(i)/crtotp
      else
      rfrr(i) = 0.0
      crfrr(i) = 0.0
      end if
      sumh = sumh + rfrh(i)
      sumg = sumg + rfrr(i)
      end do
      
      IF ( ndebug .ge. 1 ) THEN
        write(6,*) 'intyp,rfrh   = ',intyp,(rfrh(i),i=1,nr)
        write(6,*) 'intyp,rfrh2  = ',intyp,(rfrh2(i),i=1,nr)
        write(6,*) 'intyp,rfrh3  = ',intyp,(rfrh3(i),i=1,nr)
        write(6,*) 'intyp,crfrh  = ',intyp,(crfrh(i),i=1,nr)
        write(6,*) 'intyp,rfrr   = ',intyp,(rfrr(i),i=1,nr)
        write(6,*) 'qf,qr,qhtotp,qrtotp = ',                  &
                    1.e3*qfw,1.e3*qrw,qhtotp,qrtotp
      ENDIF
      
      DO i = 1,nr
       rfrh(i) = rfrh3(i)
       rfrr(i) = rfrh3(i)
      ENDDO
!
      if ( sumg .lt. .9999 .and. sumg .ne. 0. ) then
      write(*,*) 'OUCH-G',sumg
        write(6,*) 'intyp,rfrh = ',intyp,(rfrh(i),i=1,nr)
        write(6,*) 'intyp,rfrr = ',intyp,(rfrr(i),i=1,nr)
        write(6,*) 'qf,qr,qhtotp,qrtotp = ',                  &
                    1.e3*qfw,1.e3*qrw,qhtotp,qrtotp
      end if
      if ( sumh .lt. .9999 .and. sumh .ne. 0. ) then
      write(*,*) 'OUCH-H',sumh
        write(6,*) 'intyp,rfrh = ',intyp,(rfrh(i),i=1,nr)
        write(6,*) 'intyp,rfrr = ',intyp,(rfrr(i),i=1,nr)
        write(6,*) 'qf,qr,qhtotp,qrtotp = ',                  &
                    1.e3*qfw,1.e3*qrw,qhtotp,qrtotp
      end if
!
      if ( idebug .ge. -1 ) then
      end if
!
!
!      end do
!      end do
!
!
      RETURN
      END
!
!
!


