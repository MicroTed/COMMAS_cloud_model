#define TMINTEMP -43.

!
! 08.13.04: version 02b: added missing low-temp roll-off to "saund" (S91)
!           scheme and changed all roll-offs to cut off at -40 
!           instead of -43.
!
! 11.11.03: version 02a: Fixes to the S91 (saund) subroutine
!
! ####################################################################
!                SUBROUTINE SAUNDX
! ####################################################################
      subroutine saundx(isaund,temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,rarfac)
      
! ####################################################################
!
! 05.04.2004  Fixed another error in the S91 code (in region S8 of 
!             Helsdon et al.)
!             Added isaund = -2 option to replace anom. pos. zone with normal 
!              negative charging
!             Added isaund = -3 option to use same equations as Wojcik 1994
!
!  11/9/2001  Fixed a discrepancy between the crystal diameter used
!             to calculate charge separation and the diameter used
!             to find the size range.  Previously used awdia for finding
!             the size range, but then used fac*awdia for calculating the
!             charge separation.  (OK if fac = 1, but not for fac=3.67)
!
!
!
!  Purpose:
!   Calls the appropriate version of the Saunders et al. scheme
!
!  -5 : Saunders etal 1991 (following Helsdon et al. 2001, but use normal charging instead of 'anomalous' zones)
!  -4 : Saunders etal 1991 (following Helsdon et al. 2001, but remove positive 'anomalous' zone)
!  -3 : Saunders etal 1991 (following Helsdon et al. 2001)
!  -2 : Saunders etal 1991 (remove positive 'anomalous' zone)
!  -1 : Saunders etal 1991 (unmodified 'anomalous' regions) (DO NOT USE)
!   0 : Saunders 1991 (modified as in Wojcik 1994) 
!   1 : RR scheme (with extra factor 3.67**qconm)
!   2 : RR scheme ( no extra factor )
!   3 : Saunders and Peck Scheme (with extra factor 3.67**qconm) (DO NOT USE)
!   4 : Saunders and Peck Scheme ( no extra factor )
!   5 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.5*rar) (DO NOT USE)
!   6 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.75*rar) (DO NOT USE)
!   7 : Saunders and Peck Scheme ( no extra factor, 0.5*rar)
!   8 : RR scheme ( no extra factor, , 0.5*rar )
!   9 : Saunders and Peck Scheme ( no extra factor, cutoff at -32.47 as orig eq. from sp98 )
!  10 : Brooks et al. RARcrit for T > -15 using saund2 (otherwise same as isaund=2) (set rarfac to negative in saund2)
!  11 : Brooks et al. RARcrit for T > -15 using saund6 (otherwise same as isaund=4)
!  12 : Brooks et al. RARcrit for T > -15 using saund6 (otherwise same as isaund=9)
!  13 : Brooks et al. RARcrit for T > -15 using saund8 
!  14 : Brooks et al. RARcrit for T > -15 using saund8 and minimum RAR (on no liquid assumption)
! ####################################################################
      
      implicit none
      integer isaund  ! scheme choice
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real awdia1
      real rarfac

! ####################################################################
      
      idelq = 0
      awdia1 = awdia
      
      IF ( isaund .eq. 1 ) THEN
        awdia1 = 3.67*awdia
        call saund2(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,3.67,1.0)
      ELSEIF ( isaund .eq. 2 ) THEN
        call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
      ELSEIF ( isaund .eq. 8 ) THEN
        call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,0.5)
      ELSEIF ( isaund .eq. 3 ) THEN
        awdia1 = 3.67*awdia
        call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,3.67,1.0)
      ELSEIF ( isaund .eq. 4 ) THEN
        call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
      ELSEIF ( isaund .eq. 5 ) THEN
        awdia1 = 3.67*awdia
        call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,3.67,0.5)
      ELSEIF ( isaund .eq. 6 ) THEN
        awdia1 = 3.67*awdia
        call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,3.67,0.75)
      ELSEIF ( isaund .eq. 7 ) THEN
        call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,0.5)
      ELSEIF ( isaund .eq. 9 ) THEN
        idelq = 1
        call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
      ELSEIF ( isaund .eq. 10 ) THEN
        call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
      ELSEIF ( isaund .eq. 11 ) THEN
        call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
      ELSEIF ( isaund .eq. 12 ) THEN
       idelq = 1
        call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
      ELSEIF ( isaund .eq. 13 ) THEN
        call saund8(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
      ELSEIF ( isaund .eq. 14 ) THEN
        call saund8(temcg,-100.0,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
      ELSEIF ( isaund .eq. 0 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0)
      ELSEIF ( isaund .eq. -1 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,1)
      ELSEIF ( isaund .eq. -2 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,2)
      ELSEIF ( isaund .eq. -3 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,3)
      ELSEIF ( isaund .eq. -4 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,4)
      ELSEIF ( isaund .eq. -5 ) THEN
        call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,5)
      ELSE
        write(6,*) 'STOP! unsupported value of isaund = ',isaund
        STOP
      ENDIF
      RETURN
      END subroutine saundx

! ####################################################################
!                SUBROUTINE SAUNDY
! ####################################################################
      subroutine saundy(isaund,temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq)
      
! ####################################################################
!
!  4.1.2002   altered version of saundx to take rar as input instead
!             of exw.
!
!  11/9/2001  Fixed a discrepancy between the crystal diameter used
!             to calculate charge separation and the diameter used
!             to find the size range.  Previously used awdia for finding
!             the size range, but then used fac*awdia for calculating the
!             charge separation.  (OK if fac = 1, but not for fac=3.67)
!
!
!
!  Purpose:
!   Calls the appropriate version of the Saunders et al. scheme
!
!   1 : RAR scheme (with extra factor 3.67**qconm)
!   2 : RAR scheme ( no extra factor )
!   3 : Saunders and Peck Scheme (with extra factor 3.67**qconm)
!   4 : Saunders and Peck Scheme ( no extra factor )
!   5 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.5*rar)
!   6 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.75*rar)
!
!
! ####################################################################
      
      implicit none
      integer isaund  ! scheme choice
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real rar
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real awdia1

! ####################################################################
      
      awdia1 = awdia
      
      IF ( isaund .eq. 4 ) THEN
        call saund7(temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,0.0,1.0)
      ENDIF
      RETURN
      END subroutine saundy
! ####################################################################
! ####################################################################
!                SUBROUTINE TAKAX
! ####################################################################
      subroutine takax(isaund0,nt,nlc,lookup,temcg,qcw,vt,awdia,     &
     &                 rho0,ftlwc,exw,qrw,rarfac)
      
! ####################################################################
!
!  Purpose: For calling different versions of Takahashi charging 
!           routine
!  
!
! ####################################################################
      
      implicit none
      
      integer isaund0,isaund,nt,nlc
      real lookup(0:nt,0:nlc) 
      real qcw, qrw, lwc
      real rho0
      real ftlwc
      real temcg
      real vt
      real awdia,awdia1
      real exw
      real rarfac
      
      lwc = qcw
      isaund = Abs(isaund0)
      IF ( isaund0 .eq. -1 ) lwc = (qcw + qrw)
      IF ( isaund0 .eq. -2 ) lwc = (qcw + qrw)*rarfac
       
!  BEGIN EXECUTABLE CODE
      IF ( isaund .eq. 1 ) THEN ! taka size/vel. depend.
        call taka(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
      ELSEIF ( isaund .eq. 2 ) THEN ! saunders size dependence
        awdia1 = awdia !  3.67*awdia
        call taka2(nt,nlc,lookup,temcg,lwc,vt,awdia1,rho0,ftlwc)
      ELSEIF ( isaund .eq. 3 ) THEN ! no size or velocity dependence
        call taka3(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
      ELSEIF ( isaund .eq. 4 ) THEN ! taka size/vel. depend. w/ D = 3.67/lambda
        call taka4(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
      ELSEIF ( isaund .eq. 5 ) THEN ! taka size/vel. depend. and using RAR
        call takarar(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc,exw)
      ELSE
        call taka(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
      ENDIF

      END subroutine takax

! ####################################################################
!                SUBROUTINE SAUND7
! ####################################################################
      subroutine saund7(temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,qfac,rarfac)
      
! ####################################################################
!
! 4.1.2002  New version that uses the model-calculated RAR (derived
!           from the collection rate qxacw)
!
!cPurpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the 
!           drop-off starting at -23.7
!
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40)
!
! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
!         approximate using the mass-weighted mean crystal diameter
!
! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
!          for all temperatures (took out alternate eq. for T<-34)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
!      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real qfac
      real rarfac  ! factor to reduce rar for purposes of charge calculation
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin
      parameter (tema = -7.0 , tmin = TMINTEMP)
      
      real fac
      
! ####################################################################
! Begin Executable code
! ####################################################################

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      ftrar = 0.0
      
!      rar = exw*qcw*1.0e3*rho0*vt*rarfac
      t = temcg

      IF (t .gt. tc ) THEN
      rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &  t*(7.4754e-3 + t*(5.4686e-4 +     &
     &  t*(1.6737e-5 + t*1.7613e-7)))))
!      ELSE ! IF ( t .le. -23.0) THEN
!       rarc = 3.27
!      ELSE ! IF ( t .le. -30.0) THEN
!       rarc = 1.795
      ELSE ! IF ( t .le. tc) THEN
       rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc + 40) )**3) )
      END IF


      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunders and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
!         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
!        ftrar = Max(0.0, ftrar)

!        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
!         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
!        ELSE
!         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
!        END IF
        
!        IF ( temcg .gt. -34.0 ) THEN
          ftrar = 6.74*( rar - rarc )
!        ELSE
!          ftrar = Max ( 0.0, 6.74*( rar - 
!     :     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
!        END IF
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF

      ELSEIF (qsign .lt. -0.5) THEN

!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
!        ELSEIF (rar .lt. 0.4) THEN
!         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
!         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if

        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF
        
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund7


! ####################################################################
!                SUBROUTINE TAKA2
! ####################################################################
      subroutine TAKA2(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)
      
! ####################################################################
!
!  Purpose:
!           Noninductive charging per Takahashi 1978 lab results.  
!           Using lookup table from Wojcik (1994).  For temperatures
!           lower than -30 C, the values at -30 C are used.
!
!  7/31/99 Now limit low-temp charging rates by a parabolic function
!          with a value of zero at -30 C and goes to zero at -43 C.
!         
!  8/10/99 Fixed bug in the interpolation: was accidentally 
!          _extrapolating_.   Oops.
!
!  10/26/99  Test using a Saunders et al. size dependence on the crystal size
!
! ####################################################################
      
      implicit none
      
      integer nt,nlc
      real lookup(0:nt,0:nlc) 
      real qlwc
      real qcw
      real rho0
      real ftlwc
      real temt,temcg
      integer item,ilwc
      real fact,facl
      real a1,a2,a3
      real alf
      real vt
      real awdia
      real d0
      parameter(d0 = 100.0e-6)
      real v0 
      parameter(v0 = 8.0)
      real fac
      real qconkq, qconm, qconn
      
!  BEGIN EXECUTABLE CODE

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      
      temt = Max (temcg, -30.0)
      item = Int(-temt) + 1
      qlwc = qcw*(1.e3)*rho0
      
      qlwc = Min ( qlwc, 30.0)
      ilwc = 0
      ftlwc = 0.0
      
      IF (qlwc .ge. 10.0) THEN
        ilwc = Int( qlwc/10.0 ) + 27
        facl = (qlwc - lookup(0,ilwc) )*0.1
      ELSEIF (qlwc .ge. 1.0) THEN
        ilwc = Int( qlwc ) + 18
        facl = (qlwc - lookup(0,ilwc) )*1.0
      ELSEIF (qlwc .ge. 0.1) THEN
        ilwc = Int( qlwc*10.0 ) + 9
        facl = (qlwc - lookup(0,ilwc) )*10.0
      ELSEIF (qlwc .ge. 0.01) THEN
        ilwc = Int ( qlwc*100.0 )
        facl = (qlwc - lookup(0,ilwc) )*100.0
      ELSE 
         ilwc = 0
         ftlwc = 0.0
      END IF

      IF (ilwc .gt. 0) THEN
        fact = -temt - Float(item) + 1.0
!        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
!        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
        a1 = (1.0-fact)*lookup(item,ilwc) +      &
     &          fact*lookup(Min(nt,item+1),ilwc)
        a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
     &          fact*lookup(Min(nt,item+1),ilwc+1)
        ftlwc = fac*(a1 + facl*(a2-a1))
        IF ( ftlwc .gt. 0.0 ) THEN
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
        ELSEIF ( ftlwc .lt. 0.0 ) THEN
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if
        
        ENDIF
        
!        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
!        alf = qconkq*(3.67*awdia)**qconm*(0.333*Abs(vt))**qconn ! 11/9/01 erm
!        alf = qconkq*(awdia)**qconm*(Abs(vt))**qconn        ! 11/9/01 erm
! Need to figure out what to do with qconkq, if anything...
        alf = (awdia/d0)**qconm*(Abs(vt)/v0)**qconn        ! 12/9/04 erm
        ftlwc = 1.0e-15*alf*ftlwc
      END IF
      
      RETURN
      END subroutine TAKA2

! ####################################################################
!                SUBROUTINE TAKARAR
! ####################################################################
      subroutine TAKARAR(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc,     &
     &                   exw)
      
! ####################################################################
!
!  Purpose:
!           Noninductive charging per Takahashi 1978 lab results.  
!           Using lookup table from Wojcik (1994).  For temperatures
!           lower than -30 C, the values at -30 C are used.
!
!  7/31/99 Now limit low-temp charging rates by a parabolic function
!          with a value of zero at -30 C and goes to zero at -43 C.
!         
!  8/10/99 Fixed bug in the interpolation: was accidentally 
!          _extrapolating_.   Oops.
!
! ####################################################################
      
      implicit none
      
      integer nt,nlc
      real lookup(0:nt,0:nlc) 
      real qlwc
      real qcw
      real rho0
      real ftlwc
      real exw
      real temt,temcg
      integer item,ilwc
      real fact,facl
      real a1,a2,a3
      real alf
      real vt
      real awdia
      real d0
      parameter(d0 = 100.0e-6)
      real v0 
      parameter(v0 = 8.0)
      real fac
      real rar
      
!  BEGIN EXECUTABLE CODE

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
!
! here, qlwc is scaled by vt/9.0, where 9.0m/s is the riming rod speed in Takahashi (1978)
!
      qlwc = exw*qcw*1.0e3*rho0*vt/9.0
      
      temt = Max (temcg, -30.0)
      item = Int(-temt) + 1
!      qlwc = qcw*(1.e3)*rho0
      
      qlwc = Min ( qlwc, 30.0)
      ilwc = 0
      ftlwc = 0.0
      
      IF (qlwc .ge. 10.0) THEN
        ilwc = Int( qlwc/10.0 ) + 27
        facl = (qlwc - lookup(0,ilwc) )*0.1
      ELSEIF (qlwc .ge. 1.0) THEN
        ilwc = Int( qlwc ) + 18
        facl = (qlwc - lookup(0,ilwc) )*1.0
      ELSEIF (qlwc .ge. 0.1) THEN
        ilwc = Int( qlwc*10.0 ) + 9
        facl = (qlwc - lookup(0,ilwc) )*10.0
      ELSEIF (qlwc .ge. 0.01) THEN
        ilwc = Int ( qlwc*100.0 )
        facl = (qlwc - lookup(0,ilwc) )*100.0
      ELSE 
         ilwc = 0
         ftlwc = 0.0
      END IF
      
      IF ( ilwc .eq. nlc ) THEN
       write(0,*) 'Warning: ilwc = nlc! setting to nlc-1'
       ilwc = nlc - 1
      ENDIF 
      IF ( item .gt. nt ) THEN
       write(0,*) 'Warning: item > nt! setting to nt'
       item = nt
      ENDIF 
      IF ( ilwc .lt. 0 ) THEN
       write(0,*) 'Warning: ilwc < 0! setting to 0'
       ilwc = 0
      ENDIF 
      IF ( item .lt. 0 ) THEN
       write(0,*) 'Warning: item < 0! setting to 0'
       item = 0
      ENDIF 

      IF (ilwc .gt. 0) THEN
        fact = -temt - Float(item) + 1.0
        a1 = (1.0-fact)*lookup(item,ilwc) +      &
     &          fact*lookup(Min(nt,item+1),ilwc)
        a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
     &          fact*lookup(Min(nt,item+1),ilwc+1)
        ftlwc = fac*(a1 + facl*(a2-a1))
        
        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
        ftlwc = 1.0e-15*alf*ftlwc
      END IF
      
      RETURN
      END subroutine takarar
! ####################################################################
!                SUBROUTINE TAKA
! ####################################################################
      subroutine TAKA(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)
      
! ####################################################################
!
!  Purpose:
!           Noninductive charging per Takahashi 1978 lab results.  
!           Using lookup table from Wojcik (1994).  For temperatures
!           lower than -30 C, the values at -30 C are used.
!
!  7/31/99 Now limit low-temp charging rates by a parabolic function
!          with a value of zero at -30 C and goes to zero at -43 C.
!         
!  8/10/99 Fixed bug in the interpolation: was accidentally 
!          _extrapolating_.   Oops.
!
! ####################################################################
      
      implicit none
      
      integer nt,nlc
      real lookup(0:nt,0:nlc) 
      real qlwc
      real qcw
      real rho0
      real ftlwc
      real temt,temcg
      integer item,ilwc
      real fact,facl
      real a1,a2,a3
      real alf
      real vt
      real awdia
      real d0
      parameter(d0 = 100.0e-6)
      real v0 
      parameter(v0 = 8.0)
      real fac
      
!  BEGIN EXECUTABLE CODE

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      
      temt = Max (temcg, -30.0)
      item = Int(-temt) + 1
      qlwc = qcw*(1.e3)*rho0
      
      qlwc = Min ( qlwc, 30.0)
      ilwc = 0
      ftlwc = 0.0
      
      IF (qlwc .ge. 10.0) THEN
        ilwc = Int( qlwc/10.0 ) + 27
        facl = (qlwc - lookup(0,ilwc) )*0.1
      ELSEIF (qlwc .ge. 1.0) THEN
        ilwc = Int( qlwc ) + 18
        facl = (qlwc - lookup(0,ilwc) )*1.0
      ELSEIF (qlwc .ge. 0.1) THEN
        ilwc = Int( qlwc*10.0 ) + 9
        facl = (qlwc - lookup(0,ilwc) )*10.0
      ELSEIF (qlwc .ge. 0.01) THEN
        ilwc = Int ( qlwc*100.0 )
        facl = (qlwc - lookup(0,ilwc) )*100.0
      ELSE 
         ilwc = 0
         ftlwc = 0.0
      END IF
      
      IF ( ilwc .eq. nlc ) THEN
       write(0,*) 'Warning: ilwc = nlc! setting to nlc-1'
       ilwc = nlc - 1
      ENDIF 
      IF ( item .gt. nt ) THEN
       write(0,*) 'Warning: item > nt! setting to nt'
       item = nt
      ENDIF 
      IF ( ilwc .lt. 0 ) THEN
       write(0,*) 'Warning: ilwc < 0! setting to 0'
       ilwc = 0
      ENDIF 
      IF ( item .lt. 0 ) THEN
       write(0,*) 'Warning: item < 0! setting to 0'
       item = 0
      ENDIF 

      IF (ilwc .gt. 0) THEN
        fact = -temt - Float(item) + 1.0
        a1 = (1.0-fact)*lookup(item,ilwc) +      &
     &          fact*lookup(Min(nt,item+1),ilwc)
        a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
     &          fact*lookup(Min(nt,item+1),ilwc+1)
        ftlwc = fac*(a1 + facl*(a2-a1))
        
        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
        ftlwc = 1.0e-15*alf*ftlwc
      END IF
      
      RETURN
      END subroutine TAKA

! ####################################################################
!                SUBROUTINE TAKA3
! ####################################################################
      subroutine TAKA3(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)
      
! ####################################################################
!
!  Purpose:
!           Noninductive charging per Takahashi 1978 lab results.  
!           Using lookup table from Wojcik (1994).  For temperatures
!           lower than -30 C, the values at -30 C are used.
!
!  7/31/99 Now limit low-temp charging rates by a parabolic function
!          with a value of zero at -30 C and goes to zero at -43 C.
!         
!  8/10/99 Fixed bug in the interpolation: was accidentally 
!          _extrapolating_.   Oops.
!
! ####################################################################
      
      implicit none
      
      integer nt,nlc
      real lookup(0:nt,0:nlc) 
      real qlwc
      real qcw
      real rho0
      real ftlwc
      real temt,temcg
      integer item,ilwc
      real fact,facl
      real a1,a2,a3
      real alf
      real vt
      real awdia
      real d0
      parameter(d0 = 100.0e-6)
      real v0 
      parameter(v0 = 8.0)
      real fac
      
!  BEGIN EXECUTABLE CODE

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      
      temt = Max (temcg, -30.0)
      item = Int(-temt) + 1
      qlwc = qcw*(1.e3)*rho0
      
      qlwc = Min ( qlwc, 30.0)
      ilwc = 0
      ftlwc = 0.0
      
      IF (qlwc .ge. 10.0) THEN
        ilwc = Int( qlwc/10.0 ) + 27
        facl = (qlwc - lookup(0,ilwc) )*0.1
      ELSEIF (qlwc .ge. 1.0) THEN
        ilwc = Int( qlwc ) + 18
        facl = (qlwc - lookup(0,ilwc) )*1.0
      ELSEIF (qlwc .ge. 0.1) THEN
        ilwc = Int( qlwc*10.0 ) + 9
        facl = (qlwc - lookup(0,ilwc) )*10.0
      ELSEIF (qlwc .ge. 0.01) THEN
        ilwc = Int ( qlwc*100.0 )
        facl = (qlwc - lookup(0,ilwc) )*100.0
      ELSE 
         ilwc = 0
         ftlwc = 0.0
      END IF

      IF (ilwc .gt. 0) THEN
        fact = -temt - Float(item) + 1.0
!        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
!        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
        a1 = (1.0-fact)*lookup(item,ilwc) +      &
     &          fact*lookup(Min(nt,item+1),ilwc)
        a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
     &          fact*lookup(Min(nt,item+1),ilwc+1)
        ftlwc = fac*(a1 + facl*(a2-a1))
        
!        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
        ftlwc = 1.0e-15*ftlwc
      END IF
      
      RETURN
      END subroutine TAKA3

! ####################################################################
!                SUBROUTINE TAKA4
! ####################################################################
      subroutine TAKA4(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)
      
! ####################################################################
!
!  Purpose:
!           Noninductive charging per Takahashi 1978 lab results.  
!           Using lookup table from Wojcik (1994).  For temperatures
!           lower than -30 C, the values at -30 C are used.
!
!  7/31/99 Now limit low-temp charging rates by a parabolic function
!          with a value of zero at -30 C and goes to zero at -43 C.
!         
!  8/10/99 Fixed bug in the interpolation: was accidentally 
!          _extrapolating_.   Oops.
!
! ####################################################################
      
      implicit none
      
      integer nt,nlc
      real lookup(0:nt,0:nlc) 
      real qlwc
      real qcw
      real rho0
      real ftlwc
      real temt,temcg
      integer item,ilwc
      real fact,facl
      real a1,a2,a3
      real alf
      real vt
      real awdia
      real d0
      parameter(d0 = 100.0e-6)
      real v0 
      parameter(v0 = 8.0)
      real fac
      
!  BEGIN EXECUTABLE CODE

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      
      temt = Max (temcg, -30.0)
      item = Int(-temt) + 1
      qlwc = qcw*(1.e3)*rho0
      
      qlwc = Min ( qlwc, 30.0)
      ilwc = 0
      ftlwc = 0.0
      
      IF (qlwc .ge. 10.0) THEN
        ilwc = Int( qlwc/10.0 ) + 27
        facl = (qlwc - lookup(0,ilwc) )*0.1
      ELSEIF (qlwc .ge. 1.0) THEN
        ilwc = Int( qlwc ) + 18
        facl = (qlwc - lookup(0,ilwc) )*1.0
      ELSEIF (qlwc .ge. 0.1) THEN
        ilwc = Int( qlwc*10.0 ) + 9
        facl = (qlwc - lookup(0,ilwc) )*10.0
      ELSEIF (qlwc .ge. 0.01) THEN
        ilwc = Int ( qlwc*100.0 )
        facl = (qlwc - lookup(0,ilwc) )*100.0
      ELSE 
         ilwc = 0
         ftlwc = 0.0
      END IF

      IF (ilwc .gt. 0) THEN
        fact = -temt - Float(item) + 1.0
!        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
!        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
        a1 = (1.0-fact)*lookup(item,ilwc) +      &
     &          fact*lookup(Min(nt,item+1),ilwc)
        a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
     &          fact*lookup(Min(nt,item+1),ilwc+1)
        ftlwc = fac*(a1 + facl*(a2-a1))
        
        alf = Min(5.0*(3.67*awdia/d0)**2*Abs(vt)/v0,10.0)
        ftlwc = 1.0e-15*alf*ftlwc
      END IF
      
      RETURN
      END subroutine TAKA4

! ####################################################################
!                SUBROUTINE SAUND6
! ####################################################################
      subroutine saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,qfac,rarfac)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
!              values of rarc for -43 < T < -40
!
!
! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the 
!           drop-off starting at -23.7
!
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40) 
!         (correction, 10/12/2003: charging was set zero below -43, not -40)
!
! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
!         approximate using the mass-weighted mean crystal diameter
!
! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
!          for all temperatures (took out alternate eq. for T<-34)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real qfac
      real rarfac  ! factor to reduce rar for purposes of charge calculation
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin
      parameter (tema = -7.0 , tmin = -37.0)
      integer ibs
      
      real fac
      
! ####################################################################
! Begin Executable code
! ####################################################################

      IF ( temcg .gt. -30. ) THEN
        fac = 1.0
!      ELSEIF ( temcg .gt. -40.0 ) THEN
!        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSEIF ( temcg .gt. tmin ) THEN
        fac = 1.0 - ((temcg+30.0)/(-tmin - 30.0))**2
      ELSE
        fac = 0.0
      END IF

      IF ( rarfac .gt. 0.0 ) THEN
        ibs = 0
      ELSE
        ibs = 1
      ENDIF
      
      ftrar = 0.0
      
      rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
      t = temcg

      IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
       rarc = Max( 0.0, Min( 3.29, -1.47 - 0.2*t ) )
      ELSE
        IF (t .gt. tc ) THEN
        rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7)))))
!      ELSE ! IF ( t .le. -23.0) THEN
!       rarc = 3.27
!      ELSE ! IF ( t .le. -30.0) THEN
!       rarc = 1.795
        ELSE ! IF ( t .le. tc) THEN
         rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc - tmin) )**3) )
        ENDIF
      ENDIF
      
       rarc = Max(rarc, rar0)  ! ERM 10/12/2003; changed 0 to rar0 5/8/2004

!
! New option (2/9/2005) to use original SP98 equation and cut off charging where
!  it hits zero
!
      IF ( idelq .eq. 1 ) THEN  
      
        IF ( t .le. tc) THEN
        rarc = Max(0.0, 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7))))) )
        END IF
        
        IF ( t .lt. -33. ) rarc = rar0
      
        IF ( temcg .gt. -25.0 ) THEN
          fac = 1.0
        ELSEIF ( temcg .gt. -32.47 ) THEN
          fac = 1.0 + ((temcg+25.0)/(32.47 - 25.00))
        ELSE
          fac = 0.0
          qsign = 0.0
          ftrar = 0.0
          qconkq = 0.0
          qconn = 1.0
          qconm = 1.0
          GOTO 999
        END IF
      
      ENDIF

      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunders and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
          ftrar = 6.74*( rar - rarc )
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF

      ELSEIF (qsign .lt. -0.5) THEN


        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if

        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF
        
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund6


! ####################################################################
!                SUBROUTINE SAUND8
! ####################################################################
      subroutine saund8(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,qfac,rarfac)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
!              values of rarc for -43 < T < -40
!
!
! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the 
!           drop-off starting at -23.7
!
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40) 
!         (correction, 10/12/2003: charging was set zero below -43, not -40)
!
! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
!         approximate using the mass-weighted mean crystal diameter
!
! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
!          for all temperatures (took out alternate eq. for T<-34)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real, intent(in) :: qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real qfac
      real rarfac  ! factor to reduce rar for purposes of charge calculation
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
!      parameter (tc = -23.7)
      parameter (tc = -24.5)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin,tcc
      parameter (tema = -7.0 , tmin = -37.0)
      real :: tmincutoff
      integer ibs
      real, parameter :: pi = 3.141592654
      
      real fac
      
! ####################################################################
! Begin Executable code
! ####################################################################

      tcc = -30.
      tmincutoff = -37.
      IF ( temcg .gt. tcc ) THEN
        fac = 1.0
!      ELSEIF ( temcg .gt. -40.0 ) THEN
!        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSEIF ( temcg .gt. tmincutoff ) THEN
!        fac = 1.0 - ((temcg+30.0)/(-tmin - 30.0))**2
! Cosine funtion roll-off
         fac = 0.5 *(1. + Cos(pi *(((temcg - tcc)/(-tmincutoff + tcc)))))
      ELSE
        fac = 0.0
      END IF

      IF ( rarfac .gt. 0.0 ) THEN
        ibs = 0
      ELSE
        ibs = 1
      ENDIF
      
      ftrar = 0.0
      
      IF ( qcw > 0.0 ) THEN
        rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
      ELSE
        rar = 0.0 ! 1.5*rar0
      ENDIF
      t = temcg

      IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
       rarc = Max( 0.0, Min( 3.29, -1.47 - 0.2*t ) )
      ELSE
        IF (t .gt. tc ) THEN
        rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7)))))
!      ELSE ! IF ( t .le. -23.0) THEN
!       rarc = 3.27
!      ELSE ! IF ( t .le. -30.0) THEN
!       rarc = 1.795
        ELSEIF ( t > tmin) THEN
!         rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc - tmin) )**3) )
         rarc = 3.42533* 0.5 *(1. + Cos(pi *(((t - tc)/(-tmin + tc))))) ! 3.42533 is the value on the curve at T = -24.5284
        ELSE
          rarc = 0
        ENDIF
      ENDIF
      
       rarc = Max(rarc, rar0)  ! ERM 10/12/2003; changed 0 to rar0 5/8/2004

!
! New option (2/9/2005) to use original SP98 equation and cut off charging where
!  it hits zero
!
      IF ( idelq .eq. 1 ) THEN  
      
        IF ( t .le. tc) THEN
        rarc = Max(0.0, 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7))))) )
        END IF
        
        IF ( t .lt. -33. ) rarc = rar0
      
        IF ( temcg .gt. -25.0 ) THEN
          fac = 1.0
        ELSEIF ( temcg .gt. -32.47 ) THEN
          fac = 1.0 + ((temcg+25.0)/(32.47 - 25.00))
        ELSE
          fac = 0.0
          qsign = 0.0
          ftrar = 0.0
          qconkq = 0.0
          qconn = 1.0
          qconm = 1.0
          GOTO 999
        END IF
      
      ENDIF

      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunders and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
          ftrar = 6.74*( rar - rarc )
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF

      ELSEIF (qsign .lt. -0.5) THEN


        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if

        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF
        
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund8

! ####################################################################
!                SUBROUTINE SAUND2
! ####################################################################
      subroutine saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq,qfac,rarfac)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 8.27.2002 added rarfac input
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40)
!  version for small droplet curve (Saunders, ICAE 1999)
! 8/15/99 use higher RARcrit at low temp and keep small droplet
!
! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
!         approximate using the mass-weighted mean crystal diameter
!
! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
!          for all temperatures (took out alternate eq. for T<-34)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real qfac
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin
      parameter (tema = -7.0 , tmin = TMINTEMP)
      integer ibs
      
      real fac,rarfac
      
! ####################################################################
! Begin Executable code
! ####################################################################

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      IF ( rarfac .gt. 0.0 ) THEN
        ibs = 0
      ELSE
        ibs = 1
      ENDIF
      
      ftrar = 0.0
      
      rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
      t = temcg
      IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
       rarc = Max( rar0, Min( 3.29, -1.47 - 0.2*t ) )
      ELSE
        IF (t .gt. -7.0 ) THEN
         rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7)))))
        ELSEIF ( t .gt. -16.0 ) THEN
         rarc =  8.0*Abs((tema-t)/(10.0+tema))*     &
     &      Exp(- Abs((tema-t)/(10.0+tema))) +     &
     &     1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &    t*(7.4754e-3 + t*(5.4686e-4 +     &
     &    t*(1.6737e-5 + t*1.7613e-7)))))
        ELSEIF ( t .gt. -21.7 ) THEN
          rarc = 4*( 1.0 - Abs( ( (t + 25.0 )/(-25.0 - tmin) )**2) )
        ELSE 
          rarc = 4*( 1.0 - Abs( ( (-21.7 + 25.0 )/(-25.0 - tmin) )**2) )
        ENDIF
      ENDIF ! ibs .eq. 1
      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunder and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
!         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
!        ftrar = Max(0.0, ftrar)

!        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
!         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
!        ELSE
!         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
!        END IF
        
!        IF ( temcg .gt. -34.0 ) THEN
          ftrar = 6.74*( rar - rarc )
!        ELSE
!          ftrar = Max ( 0.0, 6.74*( rar - 
!     :     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
!        END IF
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF

      ELSEIF (qsign .lt. -0.5) THEN

!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
!        ELSEIF (rar .lt. 0.4) THEN
!         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
!         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if

        IF ( qfac .gt. 1.0 ) THEN
         ftrar = fac*ftrar*qfac**qconm
        ELSE
         ftrar = fac*ftrar
        ENDIF
        
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund2


! ####################################################################
!                SUBROUTINE SAUND5
! ####################################################################
      subroutine saund5(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40)
!  version for small droplet curve (Saunders, ICAE 1999)
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin
      parameter (tema = -7.0 , tmin = TMINTEMP)
      
! ####################################################################
! Begin Executable code
! ####################################################################
      
      ftrar = 0.0
      
      rar = exw*qcw*1.0e3*rho0*vt
      t = temcg
      IF (t .gt. -7.0 ) THEN
      rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &  t*(7.4754e-3 + t*(5.4686e-4 +     &
     &  t*(1.6737e-5 + t*1.7613e-7)))))
      ELSEIF ( t .gt. -16.0 ) THEN
       rarc =  8.0*Abs((tema-t)/(10.0+tema))*     &
     &      Exp(- Abs((tema-t)/(10.0+tema))) +     &
     &     1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &  t*(7.4754e-3 + t*(5.4686e-4 +     &
     &  t*(1.6737e-5 + t*1.7613e-7)))))
      ELSE
       rarc = 4*( 1.0 - Abs( ( (t + 25.0 )/(-25.0 - tmin) )**2) )
      END IF
      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunder and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
!         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
!        ftrar = Max(0.0, ftrar)

!        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
!         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
!        ELSE
!         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
!        END IF
        
        IF ( temcg .gt. -34.0 ) THEN
          ftrar = 6.74*( rar - rarc )
        ELSE
          ftrar = Max ( 0.0, 6.74*( rar -      &
     &     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
        END IF
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if

      ELSEIF (qsign .lt. -0.5) THEN

!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
!        ELSEIF (rar .lt. 0.4) THEN
!         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
!         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund5

! ####################################################################
!                SUBROUTINE SAUND4
! ####################################################################
      subroutine saund4(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
!              rarc for -43 < T < -40
!
!
! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
!          from -23.7 to -40 (and zero charging below -40)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      
! ####################################################################
! Begin Executable code
! ####################################################################
      
      ftrar = 0.0
      
      rar = exw*qcw*1.0e3*rho0*vt
      t = temcg
!      IF (t .gt. -23.0 ) THEN
!      IF (t .gt. -30.0 ) THEN
      IF (t .gt. tc ) THEN
      rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &  t*(7.4754e-3 + t*(5.4686e-4 +     &
     &  t*(1.6737e-5 + t*1.7613e-7)))))
!      ELSE ! IF ( t .le. -23.0) THEN
!       rarc = 3.27
!      ELSE ! IF ( t .le. -30.0) THEN
!       rarc = 1.795
      ELSE ! IF ( t .le. tc) THEN
       rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc + 40) )**3) )
      END IF

      
      rarc = Max(rarc, 0.0)

      
! check for RAR below threshold:
      IF ( rar .le. rar0 .or. temcg .le. -40.0) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunder and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
!         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
!        ftrar = Max(0.0, ftrar)

        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
         ftrar = Max( 0.0, delri*Abs(rar - rarc)*     &
     &           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
        ELSE
         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
        END IF
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if

      ELSEIF (qsign .lt. -0.5) THEN

!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
!        ELSEIF (rar .lt. 0.4) THEN
!         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
!         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

        ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(7.0)*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        ftrar = Min( 0.0, ftrar )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund4
      
! ####################################################################
!                SUBROUTINE SAUND3
! ####################################################################
      subroutine saund3(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
     &                 qconkq,qconm,qconn,idelq)
      
! ####################################################################
!
!  Purpose:
!   Saunders charging scheme based on rime accretion rate (RAR)
!   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
!  8/2/99  Put in quadratic drop-off in transferred charge 
!          starting at -30 (fac=1.0) down to -43 (fac=0.0)
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real fac
      
! ####################################################################
! Begin Executable code
! ####################################################################
      
      ftrar = 0.0

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF
      
      rar = exw*qcw*1.0e3*rho0*vt
      t = temcg
      IF (t .gt. -23.0 ) THEN
!      IF (t .gt. -30.0 ) THEN
      rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
     &  t*(7.4754e-3 + t*(5.4686e-4 +     &
     &  t*(1.6737e-5 + t*1.7613e-7)))))
      ELSE ! IF ( t .le. -23.0) THEN
       rarc = 3.27
!      ELSE ! IF ( t .le. -30.0) THEN
!       rarc = 1.795
      END IF
      
      
! check for RAR below threshold:
      IF ( rar .le. rar0 ) THEN
        qsign = 0.0
        ftrar = 0.0
        qconkq = 0.0
        qconn = 1.0
        qconm = 1.0
        GOTO 999
      END IF

! for now, use charging values at -8.0 for warmer temps
! and values at -23.0 for colder temps (Saunder and Peck, 1998)      
      IF (temcg .gt. -8.0 ) THEN
         t = -8.0
      ELSEIF (temcg .lt. -23.0) THEN
         t = -23.0
      END IF
      
      qsign = -1.0
      IF (rar .gt. rarc) qsign = 1.0
      
      IF (qsign .gt. 0.5) THEN
        
!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
!         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
!        ftrar = Max(0.0, ftrar)

         ftrar = fac*Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if

      ELSEIF (qsign .lt. -0.5) THEN

!        IF ( Abs(rar - rarc) .lt. delr) THEN
!         ftrar = delri*Abs(rar - rarc)*
!     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
!        ELSEIF (rar .lt. 0.4) THEN
!         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
!         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

        ftrar = fac*q0*(-1.0 +      &
     &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saund3
      
!
! 6/21/99 
!         Also put in a linear drop-off for -7.4 < T < 0
!       try setting elwc = 0.5*lwc
!
!  11.11.03  ERM fixed some bugs that allowed charging at EW below
!            threshold for T < -20.  Also corrected charging at low EW
!            at -16 > T > -20.
!
! ####################################################################
!                SUBROUTINE SAUND
! ####################################################################
      subroutine saund(temcg,qcw,exw,awdia,rho0,qsign,ftelwc,     &
     &                 qconkq,qconm,qconn,idelq,ianom)
     
      implicit none
      
      integer ianom  ! =0 to reduce 'anomalous' regions
                     ! =1 to keep 'anonmalous' regions unchanged
                     ! =2 to reduce neg. anom. zone and
                     !     remove the pos. anom. zone (use normal neg. instead)
                     ! =3 to run same as Helsdon et al. 2001
                     ! =4 to run Helsdon et al. 2001 without pos. anom. zone
                     ! =5 to run Helsdon et al. 2001 without pos. or neg. anom. zones (use normal instead)
                     ! =(anything else) to remove the anomalous zones altogether
      real anom
      integer ieq, iseq

      real temcg ! temperature
      real qcw ! cloud water mixing ratio 
      real exw ! cloud water collection efficiency
!      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
      real ftelwc, fac
      
      real elwc, trevsau, qlwc, pi, cew
!
! ####################################################
!

      IF ( temcg .gt. -30 ) THEN
        fac = 1.0
      ELSEIF ( temcg .gt. -40.0 ) THEN
        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
      ELSE
        fac = 0.0
      END IF

      ieq = 0
      pi    = 4.*atan(1.0)
!
      qlwc = qcw*(1.e3)*rho0
      elwc = qlwc*Min(1.0,exw)
!      elwc = qlwc*0.5
!
!
!   Saunders et al. 1991 non-inductive parameterizations
!
!
      idelq = 0
      ftelwc = 0.0
      qsign = 0.0
      qconkq = 0.0
      qconm = 1.0
      qconn = 1.0

      IF ( ( temcg .lt. -20.0 .and. elwc .le. 0.061 ) .or.     &
     &     (  elwc .le. 0.026 ) ) THEN
       RETURN
      ENDIF
!
!   case I  (positive anomalous zone)
!
      IF ( ( temcg.lt.(-20.0) .or.      &
     &    ( temcg.lt.(-18.0) .and. ianom .ne. 3 ) )     &
     &      .and. elwc.lt.(0.1596) .and. ianom .ne. 2      &
     &       .and. ianom .ne. 4 .and. ianom .ne. 5 ) then
!
      ieq = 1
      idelq = 1
!
!  sign
!
      qsign = 1.0
!
!  f(t,lwc)
!
      ftelwc = 0.0
      
      if (elwc.lt.(0.12).and.elwc.gt.(0.061)) then
      ftelwc = ( (2041.76)*elwc-(128.70) )
      end if
      if (elwc.lt.(0.1596).and.elwc.ge.(0.12)) then
      ftelwc = ( (-2900.22)*elwc+(462.91) )
      end if
      
      IF ( temcg .gt. -20.0 ) THEN
       ftelwc = ftelwc*( -temcg - 18.0 )*0.5
      ENDIF
!
! f(awdia)
!
      
      IF ( ianom .eq. 0 .or. ianom .eq. 3 ) THEN
         anom = 0.1
      ELSEIF ( ianom .eq. 1 ) THEN
         anom = 1.0
      ELSE
         anom = 0.0
      ENDIF
      
      if ( awdia*1.e6.lt.155. ) then
      qconkq = anom*4.9e13
      qconm  = 3.76
      qconn  = 2.5
      end if
      if ( awdia*1.e6.ge.155. .and.      &
     &     awdia*1.e6.le.452. ) then
      qconkq = anom*4.9e6
      qconm  = 1.9
      qconn  = 2.5
      end if
      if ( awdia*1.e6.gt.452. ) then
      qconkq = anom*52.8
      qconm  = 0.44
      qconn  = 2.5
      end if
!
!  end case I
!
      end if
!
!
!   case II  (zero charging at low EW)
!
      if ( ( (temcg.ge.(-20.0) .or. (ianom.eq.2) ) .and.      &
     &      temcg.lt.(-16.0)) .and.      &
     &     elwc.lt.(0.06) ) then
!
      ieq = 2
      idelq = 1
!
!  sign
!
      qsign = 0.0
!
!  f(t,lwc)
!
      ftelwc = 0.0
!
! f(awdia)
!
      qconkq = 0.0
      qconm  = 1.0
      qconn  = 1.0
!
!  end case II
!
      end if
!
!   case III  (negative anomalous zone)
!
      if ( ( temcg.ge.(-16.0) .or.      &
     &      ( temcg.ge.(-18.0) .and. ianom .ne. 3 .and. ianom .ne. 4) )     &
     &        .and. elwc.lt.(0.22) .and. ianom .ne. 5 ) then
!
      ieq = 3
      idelq = 1
!
!  sign
!
      qsign = -1.
!
!  f(t,lwc)
!
!  If elwc.lt.0.03:
      ftelwc = 0.0                                  ! ERM bug fix 11.11.2003

      if (elwc.lt.(0.14).and.elwc.gt.(0.026)) then
       ftelwc = ( (-314.40)*elwc+(7.92) ) !iseq=2
      end if
      if (elwc.lt.(0.22).and.elwc.ge.(0.14)) then
       ftelwc = ( (419.4)*elwc-(92.64) ) !iseq=3
      end if

      ftelwc = Min ( 0.0, ftelwc ) 
      
      IF ( temcg .lt. -16.0 ) THEN
        ftelwc = ftelwc*( temcg + 18.0 )*0.5 
      ENDIF
      
      IF ( temcg .gt. -7.38 ) ftelwc = ftelwc*Abs(temcg/7.38)
!
! f(awdia)
!
      
      IF ( ianom .eq. 0 .or. (ianom .ge. 2 .and. ianom .le. 4) ) THEN
         anom = 0.2
      ELSEIF ( ianom .eq. 1 ) THEN
         anom = 1.0
      ELSE
         anom = 0.0
      ENDIF

      if ( awdia*1.e6.lt.253. ) then
      qconkq = anom*5.24e8
      qconm  = 2.54
      qconn  = 2.8
      end if
      if ( awdia*1.e6.gt.253. ) then
      qconkq = anom*24.0
      qconm  = 0.50
      qconn  = 2.8
      end if
!
!  end case III
!
      end if
!
!  ASSIGN NORMAL CHARGING ZONES
!
      trevsau = -15.06*elwc-7.38
      
      iseq = 7 ! set default for [something]
      IF ( temcg .gt. -10.69 .and. elwc .gt. 0.22 ) THEN
        cew = 0.22
        iseq = 1
      ELSEIF (temcg .lt. -24. .and. elwc .gt. 1.1 ) THEN
        cew = 1.1
        iseq = 8
      ELSEIF (temcg .lt. -24.0 .and. elwc .le. 1.1 ) THEN
        cew = 1.1
        iseq = 7
      ELSEIF ( temcg .le. -10.69 .and. temcg .ge. -24.0 ) THEN
         cew = -0.49 - (6.64e-2)*temcg
         IF ( elwc .gt. cew ) THEN
           iseq = 4
         ELSE
           iseq = 7
         ENDIF
      ENDIF
 
      IF ( temcg .le. -7.38 .and. temcg .ge. -24.0 .and. ianom .eq. 5 ) THEN
         cew = -0.49 - (6.64e-2)*temcg
         IF ( elwc .gt. cew ) THEN
           iseq = 4
         ELSE
           iseq = 7
         ENDIF
      ENDIF


!
!   case IV (normal positive zone, above the critical line)
!

      IF ( iseq .eq. 1 .or. iseq .eq. 8 .or. iseq .eq. 4 ) THEN
!
      ieq = 4
      idelq = 1
!
!  sign
!
      qsign = 1.0
!
!  f(t,lwc)
!
      IF ( iseq .eq. 4 ) THEN
        ftelwc = (20.22)*elwc+(1.36)*temcg+(10.05)
      ELSEIF ( iseq .eq. 1 .and. ianom /= 5 ) THEN
        ftelwc = 20.22*elwc - 4.4484 ! value of T set to -10.69
!cccc      ENDIF
      ELSEIF ( iseq .eq. 1 .and. ianom == 5 ) THEN
        ftelwc = 20.22*(elwc - cew)
      ELSEIF ( iseq .eq. 8 ) THEN
        IF ( ianom .eq. 3 .or. ianom .eq. 4 ) THEN 
          ftelwc = Max(0.0, (20.22)*elwc - 22.24 )
        ELSE  ! create region of zero charging
          ftelwc = Max(0.0, (20.22)*elwc+(1.36)*temcg+(10.05))
        ENDIF
      ENDIF
      IF ( temcg .gt. -7.38 )     &
     &   ftelwc = Abs(temcg/7.38)*     &
     &            ((20.22)*elwc + (1.36)*(-7.38) + 10.05)

!      IF ( ftelwc .lt. 0 ) ftelwc = 0.0
      ftelwc = Max ( 0.0 , ftelwc )
!
! f(awdia)
!
      if ( awdia*1.e6.lt.155. ) then
      qconkq = 4.9e13
      qconm  = 3.76
      qconn  = 2.5
      end if
      if ( awdia*1.e6.ge.155. .and.      &
     &     awdia*1.e6.le.452. ) then
      qconkq = 4.9e6
      qconm  = 1.9
      qconn  = 2.5
      end if
      if ( awdia*1.e6.gt.452. ) then
      qconkq = 52.8
      qconm  = 0.44
      qconn  = 2.5
      end if
!
!  end case IV
!
      END IF
!
!
!   case V (normal negative zone, below the critical line)
!
!      if ( temcg.le.(trevsau) .and. 
!     >     elwc.lt.(1.1) .and. 
!     >     idelq .eq. 0 ) then
        IF ( idelq .eq. 0 .and. iseq .eq. 7 ) THEN
!
      ieq = 7
      idelq = 1
!
!  sign
!
      qsign = -1.
!
!  f(t,lwc)
!
      ftelwc = Min( 0.0, (3.02)-(31.76)*elwc+(26.53)*elwc**2 )
!
! f(awdia)
!
      if ( awdia*1.e6.lt.253. ) then
      qconkq = 5.24e8
      qconm  = 2.54
      qconn  = 2.8
      end if
      if ( awdia*1.e6.gt.253. ) then
      qconkq = 24.0
      qconm  = 0.50
      qconn  = 2.8
      end if
!
!  end case V
!
      end if
      
      ftelwc = fac*ftelwc
!
!     print*,'qsign',qsign
!     print*,'ftelwc',ftelwc
!     print*,'qconkq',qconkq
!     print*,'qconm',qconm
!     print*,'qconn',qconn
!     print*,'idelq',idelq
!     print*,'awdia',awdia
!     print*,'temcg',temcg
!     print*,'qcw',qcw
!
      return
      end subroutine saund
!
!
! ####################################################################
!                SUBROUTINE SAUNDMST
! ####################################################################
      subroutine saundmst(temcg,ssi,vt,awdia,rho0,     &
     &                 qsign,ftrar,qcw,exw,     &
     &                 qconkq,qconm,qconn,idelq,rarfac)
      
! ####################################################################
!
!  Purpose:
!    Charge separation based on having just ice supersaturation but little
!    or no riming. (Mitzeva, Saunders, and Tsenova, 2006, Atmos. Res.)
!
!  NOTE: 'q' here does nothing -- ftrar is really set in the calling subroutine
!     using qsign for polarity.
!
!
! ####################################################################
      
      implicit none
      real temcg ! temperature
      real ssi ! ice supersaturation
      real qcw !  cloud liquid water mixing ratio
      real exw ! cloud water collection efficiency
      real vt ! terminal speed difference between x and cw 
      real awdia ! crystal diameter
      real rho0  ! air density
      real qsign ! sign of charge acquired by rimer (not used by calling prog)
      real ftrar ! charge factor based on temp and RAR 
      real qconkq ! factor kq ( or 'B') 
      real qconm  ! exponent on crystal diameter ('a')
      real qconn  ! exponent on speed ('b')
      integer idelq ! charge sign, not used by calling program
!      real qfac
      real rarfac  ! factor to reduce rar for purposes of charge calculation
      
      real rar ! rime accretion rate
      real rarc ! critical RAR
      real t,tc
      parameter (tc = -23.7)
      real delr,delri
      parameter (delr=0.5, delri=1.0/delr)
      real rar0  ! lower limit of RAR where charging goes to zero
      parameter (rar0 = 0.1)
      real q0
      parameter (q0 = 6.48)
      real tema,tmin
      parameter (tema = -7.0 , tmin = -37.0)
      integer ibs
      
      real fac
      real q
      
! ####################################################################
! Begin Executable code
! ####################################################################


      fac = 1.0
      q = 1.
      
      rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
      t = temcg

      IF ( rar .gt. rar0 ) THEN ! should not be here if there is riming
          fac = 0.0
          qsign = 0.0
          ftrar = 0.0
          qconkq = 0.0
          qconn = 1.0
          qconm = 1.0
          GOTO 999
        
      ENDIF
      
      
      IF ( ssi .ge. 1.0 ) THEN
        qsign = -1.0
      ELSE
        qsign =  1.0
      ENDIF
      
      IF (qsign .gt. 0.5) THEN
        
        
        if ( awdia*1.e6.lt.155. ) then
          qconkq = 4.9e13
          qconm  = 3.76
          qconn  = 2.5
        end if
        if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
          qconkq = 4.0e6
          qconm  = 1.9
          qconn  = 2.5
        end if
        if ( awdia*1.e6.gt.452. ) then
          qconkq = 52.8
          qconm  = 0.44
          qconn  = 2.5
        end if
        
         ftrar = q

      ELSEIF (qsign .lt. -0.5) THEN


        
        if ( awdia*1.e6.lt.253. ) then
          qconkq = 5.24e8
          qconm  = 2.54
          qconn  = 2.8
        end if
        if ( awdia*1.e6.gt.253. ) then
          qconkq = 24.0
          qconm  = 0.50
          qconn  = 2.8
        end if

         ftrar = -q
        
      END IF
      
      
 999  CONTINUE      
      RETURN
      END subroutine saundmst
      
