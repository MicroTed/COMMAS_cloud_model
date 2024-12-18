c
c2345678901234567890123456789012345678901234567890123456789012345678912
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  STRAKA'S ATMOSPHERIC MODEL  (SAM)
c    Designed by Jerry M. Straka
c
c  SUBROUTINE RIMCONV
c
c  A subroutine version of sam.qfxacwt.f adapted by ERM to calculate
c  conversions of graupel and frozen drops on the fly.  The look-up 
c  table method requires too many variables when cloud droplet 
c  concentration is predicted or varies with air density.  
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine rimconv(qcw,qfw,ccw,cfw,rtfw,
     :           cdfw,cnof,fwdia,cwmas,cwdia,
     :           cwdia2,vtwbar,
     :           rfr,intyp,temg,rho0,cd,cs,ds,crfr)

c
c Things needed to pass in:
c    number concentration or intercept (type of distribution?)
c    density ranges for categories?  Could hard-wire it here for now
c    Assume that qcw > qcmin and qfw > qfmin!!!!
c
c  Other considerations:
c    assume gamma distibution for droplets
c    Must decide whether to put mgs loop in here or put call inside a loop (do the latter for now)
c    efw should be variable with respect to collector size and droplet size
c
c  Should calculate cnof internally from cfw for the time when cfw is predicted.
c
c     declare
c
      implicit none
      
c      integer ngscnt,ngs
c      real qcw(ngs),qfw(ngs),ccw(ngs),cfw(ngs)
c      real cwdia(ngs),cwdia2(ngs),cwdia3(ngs),cwmas(ngs)
c      real fwdia(ngs)
c      real rtfw(ngs)    !  riming time

      real qcw,qfw,ccw,cfw
      real cwdia,cwdia2,cwmas
      real fwdia
      real rtfw    !  riming time
      real cnof
      real cdfw
      integer intyp
      
c      character*40 filout
      integer mr,nh
      parameter(nh=21) ! 41)
      parameter(mr=5)
      real hn(nh),hd(nh),hm(nh),dmdt(nh),hrim(nh),hmp(nh)
      real hq(nh),hqp(nh),hdhd(nh),dqdt(nh),vth(nh),hdp(nh)
      real rhonew(nh)
      real vtwbar,cwdn
      real vtfbar,fwdn,fslp
      real rho0,temg
      real rimf,qcmxd,qfacw,efw(nh),rhobar,qtotp,ctotp
      double precision qtotp2
      real rhos(mr),rhoe(mr),qtot(mr), ctot(mr)
      real qtot2(mr)
      real rhoss(mr),rhoes(mr)
      data rhoss /000,150,350,600,770/
      data rhoes /150,350,600,770,1001/
      real rhosg(5),rhoeg(5)
      data rhosg /000,350,600,770,9000/
      data rhoeg /350,600,770,1001,9999/

c      real rfr(mqcw,mqfw,mt,mr,mtim)
      real rfr(mr)     !  fraction of mass to each category
      double precision rfr2(mr)    ! test fraction
      real crfr(mr)     !  fraction of mass to each category
      real crfr2(mr)
c      real qffr(mqcw,mqfw,mt,mtim)
c
c
c     read
c
c      open(7,file='inputtim',form='formatted')
c      read(7,*) filout
c      read(7,*) iplay
c      read(7,*) n, ivt
c      read(7,*) nqcw, nqfw, nt, nr,ntim
c      read(7,*) fqcw, fqfw, ft
c      read(7,*) sqcw, sqfw, st
c      read(7,*) rho0s, p0s, dtdz
c      read(7,*) rimc1, rimc2
c      read(7,*) hjo
c      read(7,*) cwccn
c      read(7,*) cnof,cdfw,hrho,dtrim
c      read(7,*) (rhos(i),i=1,nr)
c      read(7,*) (rhoe(i),i=1,nr)
c      read(7,*) idebug
       
       integer idebug
       parameter ( idebug = 0 )
       real rimc1, rimc2
       parameter ( rimc1 = 300., rimc2 = 0.44 )
       real hmmin,hjo
       parameter ( hjo = 7.5*nh/41. )
       integer l,i,iread
       real dt
c       real gf1,gf2,gf3,gf4,gf5,gf6
c       real fadvisc0
       integer ivt
       real hqtot
       real sum, rim
c       real gf4ds, gf4p5, gf4br
       integer nr, n
       real hrho
       real gr, eps, ds, cs, pi
       real cd(5)
       real fac
       
       integer ifirst
       data ifirst/1/
       save ifirst
       
       integer nn,m
       parameter (nn=nh)
       real hma(nn,mr),hda(nn,mr),hdhda(nn,mr),hrhoa(mr)
       save hma,hda,hdhda,hrhoa
       real vtha(nn,mr)
       save vtha
       data hrhoa/100.,300.,500.,700.,800./

c ####################################################################
       n = nh

      pi = 4.*atan(1.)
      hmmin = 4.7e-10
c
      eps = 1.0e-20
      gr = 9.8
c      ar = 841.99666
c      br = 0.8
c      cs = 12.42
c      ds = 0.42 
       dt = Min(120., rtfw) ! ltim*dtrim
c       dt = 120.
       
       IF ( ifirst .eq. 1 ) THEN
        ifirst = 0
       
        DO m=1,mr
          hrho = hrhoa(m)
        DO l = 1,n
         hma(l,m) = hmmin*exp(3.0*(l-1)/hjo)
         hda(l,m) = (6.*hma(l,m)/(pi*hrho))**.33333333333
         hdhda(l,m) = hda(l,m)/hjo
        ENDDO
        ENDDO
                
        DO l = 1,n
         vtha(l,1) = cs*(hda(l,1)**ds) ! vth(l) = vtha(l,1)*Sqrt(1.225/rho0)
        ENDDO

        DO m = 2,mr
        DO l = 1,n
          vtha(l,m) = Sqrt( (hda(l,m)*4.0*gr) /(3.0*cd(m)) ) ! *Sqrt(hrho/rho0)
        ENDDO
        ENDDO
       
       ENDIF ! ( ifirst )
       
       hrho = hrhoa(intyp)
      
       IF ( intyp .ne. 1 ) THEN
         nr = 4
         ivt = 1
         fac = Sqrt(hrho/rho0)
         DO i = 1,nr
           rhos(i) = rhosg(i)
           rhoe(i) = rhoeg(i)
         ENDDO
       ELSE !  ( intyp .eq. 1 ) THEN  ! SA
c         hrho = 100.
         nr = 5
         ivt = 2
         fac = Sqrt(1.225/rho0)
         DO i = 1,nr
           rhos(i) = rhoss(i)
           rhoe(i) = rhoes(i)
         ENDDO
       ENDIF
       
c
c     constants
c
      if ( idebug .gt. 2 ) then
      write(*,*)
      write(*,*)'l,hm(g),hd(cm)'
      end if
c      if ( idebug .gt. 0 ) then
c      write(*,*) '=1 to write distribution d, m'
c      read(*,*) iread
c      end if
      do l = 1,n
      hm(l) = hma(l,intyp) ! hmmin*exp(3.0*(l-1)/hjo)
      hd(l) = hda(l,intyp) ! (6.*hm(l)/(pi*hrho))**.33333333333
      hdhd(l) = hdhda(l,intyp) ! hd(l)/hjo
c      if ( ivt .eq. 1 ) then
       vth(l) = vtha(l,intyp)*fac
c      vth(l) = 
c     >  Sqrt( (hrho*hd(l)*4.0*gr) /
c     <    (3.0*cdfw*rho0) )
c      end if
c      if ( ivt .eq. 2 ) then
c       vth(l) = vtha(l,intyp)*fac
c      end if
c      if(iread.eq.1)write(*,100)l,hm(l)*1000,hd(l)*100,vth(l)
      end do 
c
      do l = 1,n
      dqdt(l) = 0.0
      dmdt(l) = 0.0
      end do
c
c      qcw = sqcw + (lqcw-1)*fqcw
c      qfw = sqfw + (lqfw-1)*fqfw
      if ( idebug .gt. 2 ) then
      print*,'intype,qcw,qfw,temg ',intyp,qcw*1000.,qfw*1000.,temg
      print*,qcw,qfw,temg
      end if
c
c      fadvisc = fadvisc0*(416.16/(temg+120.0))*
c     >  (temg/296.0)**(1.5)
c
      if ( qcw .eq. 0.0 ) go to 999
      if ( qfw .eq. 0.0 ) go to 999
c
      fwdn = hrho
      cwdn = 1000.
      efw = 1.0
100   format(1x,i4,6f12.7)
c
c  compute slope
c
c      if ( qfw .gt. qfmin ) then
c      fslp = ( (rho0*qfw)/(pi*cnof*fwdn) )**(-0.25)
c      fslp = Sqrt(Sqrt(( (pi*cnof*fwdn)/(rho0*qfw))))
      fslp = 1./fwdia
c      else
c      fslp = 0.0
c      end if        
c
c  write slope if debug
c
      if ( idebug .gt. 2 ) then
      write(*,*) '1/fslp',1./fslp
      end if
c    
c  compute number conc
c  NEED TO GENERALIZE FOR GAMMA FUNCTION????
c
      hqtot = 0.0
      do l = 1,n
c      hn(l) = cnof*exp(-hd(l)*fslp)*hdhd(l)
      hn(l) = (cfw*fslp)*exp(-hd(l)*fslp)*hdhd(l)
      hq(l) = (hn(l)*hm(l)/rho0)
      end do 
c
c  compute mass growth rate of discretized solution
c
!      if ( idebug .gt. 0 ) then
c      write(*,*) 'l,hd(l)*100,dmdt(l),dqdt(l),hqtot'
c      iread = 0
c      write(*,*) '=1 to write growth stuff'
c      read(*,*) iread
!      end if
c
c      hqtot = 0.0
      do l = 1,n
      efw(l) = 1.0
      dmdt(l) = .25*pi*((cwdia+hd(l))**2)
     >  *efw(l)*rho0*qcw*abs(vtwbar-vth(l))
      dqdt(l) = (1./rho0)*hn(l)*dmdt(l)
      end do 
c
!      if ( idebug .ge. 0 ) then
!      iread = 0
c      write(*,*) '=1 to write density stuff'
c      read(*,*) iread
!      end if
c
c  compute m = m+dtrim*dmdt and hdp from m hmp
c
      do l = 1,n
      hmp(l) = hm(l) + dt*dmdt(l)      
c      hdp(l) = (6.*hmp(l)/(pi*hrho))**.33333333333
      rim = 1000.0
      if ( temg .lt. 273.15 ) then
      rim = rimc1*(-((0.5)*(1.e+06)*cwdia)
     >                *((0.60)*vth(l))
     >                /(temg-273.15))**(rimc2)
      rim = min( max(100.0, rim), 900.0)
      end if
      rhonew(l) = (hm(l)*hrho + dt*dmdt(l)*rim)
     >  / (hmp(l))
      end do 
c
      go to 998
999   continue
998   continue
c
c  find fractions of q in each density category
c

      qtot(:) = 0.
      qtot2(:) = 0.
      ctot(:) = 0.
      do i = 1,nr
      do l = 1,n
      if ( rhonew(l) .gt. rhos(i) .and. rhonew(l) .lt. rhoe(i) ) then
      qtot(i) = qtot(i) + dqdt(l)
      qtot2(i) = qtot2(i) + hm(l)*hn(l)
      ctot(i) = ctot(i) + hn(l)
      end if
      end do
      end do
c
      qtotp = 0.0
      qtotp2 = 0.0
      ctotp = 0.0
      do i = 1,nr
      qtotp = qtotp + qtot(i)
      qtotp2 = qtotp2 + qtot2(i)
      ctotp = ctotp + ctot(i)
      end do
c
      do i = 1,nr
      if ( qtotp .gt. 0.0 ) then
      rfr(i) = qtot(i)/qtotp
      crfr(i) = ctot(i)/ctotp
      else
      rfr(i) = 0.0
      crfr(i) = 0.0
      end if

      IF ( qtotp2 .gt. 0.0 ) THEN
       rfr2(i) = qtot2(i)/qtotp2
      ELSE
       rfr2(i) = 0.0
      ENDIF
      

c      if ( idebug .ge. -1 ) then
c      write(*,100) i,dt,rfr(lqcw,lqfw,lt,i,ltim),qcw,qfw,hqtot,qfacw
c      end if
      end do
c
! rescale number conc. fractions to prevent large changes
c      IF ( crfr(intyp) .lt. 0.9*rfr2(intyp) ) THEN
c       IF ( crfr(intyp) .gt. 0.0 ) THEN
c         x = 0.9*rfr2(intyp)/crfr(intyp)
c         y = (1.0 - x*crfr(intyp))/(1.0 - crfr(intyp))
c         crfr(intyp) = x*crfr(intyp)
c         DO i = 1,nr
c          IF ( i .ne. intyp ) crfr(i) = y * crfr(i)
c         ENDDO
c       ENDIF
c      ELSEIF ( crfr(intyp) .gt. 1.1*rfr2(intyp) ) THEN
c        
c      ENDIF

      IF ( idebug .ge. 1 .and. qfw .gt. 1.0e-3 .and. qcw .gt. 0.1e-3) THEN
        write(6,*) 'intyp,rfr   = ',intyp,(rfr(i),i=1,nr)
        write(6,*) 'intyp,rfr2  = ',intyp,(rfr2(i),i=1,nr)
        write(6,*) 'intyp,crfr  = ',intyp,(crfr(i),i=1,nr)
        write(6,*) 'qf,qc,qhtotp,qrtotp = ',
     :              1.e3*qfw,1.e3*qcw,cfw,ccw
      write(6,'(11(1x,e9.4))') (hn(l), l=1,nh)
      write(6,'(11(1x,f9.4))') (hd(l)*1000., l=1,nh)
      write(6,'(11(1x,f9.3))') (rhonew(l), l=1,nh)
      ENDIF
      
      sum = 0.0
      DO i = 1,nr
        rfr(i) = rfr2(i)
!        crfr(i) = rfr(i) 
        sum = sum + rfr(i)
      ENDDO

      if ( sum .lt. .9999 .and. sum .ne. 0. ) then
      write(*,*) 'OUCH, sum = ',sum,nr,rfr(1),rfr(2),rfr(3),rfr(4),rfr(5)
      end if
c
      RETURN
      END
c
c
c


