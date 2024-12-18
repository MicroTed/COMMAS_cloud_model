!--------------------------------------------------------------------------
      subroutine gaminterpsub(gaminterp, ratio, alp, luindex, ilh)
      
      USE MICRO_MODULE, only: dqiacrratioinv,dqiacralphainv, &
                              gamxinflu,nqiacrratio,nqiacralpha, &
                              dqiacralpha,dqiacrratio,maxratiolu,maxalphalu

      implicit none

      real, intent(out) :: gaminterp
     
      real, intent(in) :: ratio, alp
      integer, intent(in) :: ilh  ! 1 = graupel, 2 = hail
      integer, intent(in) :: luindex ! which argument: 
                         ! gamxinflu(i,j,1,1) = x/y
                          ! gamxinflu(i,j,2,1) = gamxinf( 2.0+alp, ratio )/y
                          ! gamxinflu(i,j,3,1) = gamxinf( 2.5+alp+0.5*bxh, ratio )/y
                          ! gamxinflu(i,j,5,1) = gamxinf( 5.0+alp, ratio )/y
                          ! gamxinflu(i,j,6,1) = gamxinf( 5.5+alp+0.5*bxh, ratio )/y

      
      real :: delx, dely, tmp1, tmp2, temp3
      integer :: i,j,ip1,jp1, il
      
      il = Abs(ilh)


           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           delx = Min(maxratiolu,ratio) - float(i)*dqiacrratio
           dely = alp - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,luindex,il) + delx*dqiacrratioinv*         &
     &                 (gamxinflu(ip1,j,luindex,il) - gamxinflu(i,j,luindex,il))
           tmp2 = gamxinflu(i,jp1,luindex,il) + delx*dqiacrratioinv*       &
     &                 (gamxinflu(ip1,jp1,luindex,il) - gamxinflu(i,jp1,luindex,il))
           
           ! interpolate along alpha; 
           
           temp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           gaminterp = temp3
           
           ! debug
!           IF ( ilh < 0 ) THEN
!             write(0,*) 'gaminterp: ',i,j,il,ratio,delx,dely,gamxinflu(i,j,luindex,il),tmp1,tmp2,temp3
!           ENDIF
           
        END SUBROUTINE gaminterpsub


!--------------------------------------------------------------------------
      real function gaminterp(ratio, alp, luindex, ilh)
      
      USE MICRO_MODULE, only: dqiacrratioinv,dqiacralphainv, &
                              gamxinflu,nqiacrratio,nqiacralpha, &
                              dqiacralpha,dqiacrratio,maxratiolu,maxalphalu,minalphalu

      implicit none

      real, intent(in) :: ratio, alp
      integer, intent(in) :: ilh  ! 1 = graupel, 2 = hail
      integer, intent(in) :: luindex ! which argument: 
                         ! gamxinflu(i,j,1,1) = x/y
                          ! gamxinflu(i,j,2,1) = gamxinf( 2.0+alp, ratio )/y
                          ! gamxinflu(i,j,3,1) = gamxinf( 2.5+alp+0.5*bxh, ratio )/y
                          ! gamxinflu(i,j,5,1) = gamxinf( 5.0+alp, ratio )/y
                          ! gamxinflu(i,j,6,1) = gamxinf( 5.5+alp+0.5*bxh, ratio )/y

      
      real :: delx, dely, tmp1, tmp2, temp3, tmpratio
      integer :: i,j,ip1,jp1, il
      
      il = Abs(ilh)


           tmpratio = Min(maxratiolu,ratio)
           IF ( ratio > 1.e6 ) THEN
             write(0,*) 'gaminterp: ratio = ',ratio
             tmp1 = 0.
             tmp2 = ratio/tmp1
           ENDIF
           i = Min(nqiacrratio,Int(tmpratio*dqiacrratioinv))
           j = Int(Max(minalphalu,Min(maxalphalu,alp))*dqiacralphainv)
           delx = Min(maxratiolu,ratio) - float(i)*dqiacrratio
           dely = Min(maxalphalu, alp ) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,luindex,il) + delx*dqiacrratioinv*         &
     &                 (gamxinflu(ip1,j,luindex,il) - gamxinflu(i,j,luindex,il))
           tmp2 = gamxinflu(i,jp1,luindex,il) + delx*dqiacrratioinv*       &
     &                 (gamxinflu(ip1,jp1,luindex,il) - gamxinflu(i,jp1,luindex,il))
           
           ! interpolate along alpha; 
           
           temp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           gaminterp = temp3
           
           ! debug
!           IF ( ilh < 0 ) THEN
!             write(0,*) 'gaminterp: ',i,j,il,ratio,delx,dely,gamxinflu(i,j,luindex,il),tmp1,tmp2,temp3
!           ENDIF
           
        END FUNCTION gaminterp


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Routine from Numerical Recipes to replace other gamma function
!  using 32-bit reals, this is accurate to 6th decimal place.
! This copy is just renamed to avoid conflict with intrinic in newer Fortran standards
      REAL FUNCTION GAMMA_SP(xx)

      implicit none
      real xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      IF ( xx <= 0.0 ) THEN
        write(0,*) 'Argument to gamma_sp must be > 0!! xx = ',xx
        STOP
      ENDIF
      
      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      DO j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      END DO
      gamma_sp = Exp(tmp + log(stp*ser/x))

      RETURN
      END FUNCTION GAMMA_SP

!--------------------------------------------------------------------------
! Routine from Numerical Recipes to replace other gamma function
!  using 32-bit reals, this is accurate to 6th decimal place.
#ifndef CM1
      REAL FUNCTION GAMMA(xx)

      implicit none
      real xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      IF ( xx <= 0.0 ) THEN
        write(0,*) 'Argument to gamma must be > 0!! xx = ',xx
        STOP
      ENDIF
      
      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      DO j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      END DO
      gamma = Exp(tmp + log(stp*ser/x))

      RETURN
      END
#endif
      DOUBLE PRECISION FUNCTION GAMMA_DPR(x)
      ! dp gamma with real input
        implicit none
        real :: x
        double precision :: xx, gamma_dp
        
        xx = x
        
        gamma_dpr = gamma_dp(xx)
        
        return
        end
        

      DOUBLE PRECISION FUNCTION GAMMA_DP(xx)

      implicit none
      double precision xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      DO j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      END DO
      gamma_dp = Exp(tmp + log(stp*ser/x))

      RETURN
      END

! #####################################################################
!
!  FUNCTION GAMMA
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
      function gammastirling(x)
      implicit none
      real x, xx, xxx, xxxx
      real gammastirling
      real aa1, aa2, aa3, aa4
      real bb1, bb2, bb3, cc1
      xx = x*x
      xxx = x*xx
      xxxx = xx*xx
      aa1 = 1./12.
      aa2 = 1./288.
      aa3 = 139./5184.
      aa4 = 571./2488320.
      bb1 = x**x
      bb2 = exp(-x)
      bb3 = sqrt(2.*3.1415936526/x)
      cc1 = bb1*bb2*bb3
      gammastirling = cc1*(1. + aa1/x + aa2/xx- aa3/xxx - aa4/xxxx)
      return
      end
!--------------------------------------------------------------------------
!  Gamma function interpolated from lookup table
!
      REAL FUNCTION GAMMAL(x)
      
      implicit none
      real x
      real gamma
      real dx
      integer ngm0,igam, i
      parameter (ngm0=3001)
      real, save :: gmoi(0:ngm0)
      real f, arg
      integer, save :: imake = 0
      
      IF ( imake .eq. 0 ) THEN
       imake = 1
       do igam = 1,ngm0
        arg = 0.01*igam
        gmoi(igam) = gamma(arg)
       end do
      ENDIF
      
      i = Int(100.*x)
      IF ( i .lt. ngm0 ) THEN
       dx = x - 0.01*float(i)
       arg = gmoi(i) + (gmoi(i+1) - gmoi(i))*dx*100.
      ELSE
       arg = gamma(x)
      ENDIF
       gammal = arg
        
      RETURN
      END
!--------------------------------------------------------------------------


!**************************** GAML ************************** 
!  From Conrad Ziegler
! **********************************************************
      real FUNCTION GAML(ig,x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      real gamxg(10,4), xg(10)
      DATA xg/.01,.05,.1,.5,1.,2.,4.,6.,8.,10./ 
      DATA gamxg/9.95E-3,4.87E-2,9.51E-2,.393,.632,.864,  &
     & .981,.997,.999,.999,                                &
     & 4.96E-5,1.2E-3,4.67E-3,9.02E-2,.264,.593,.908,      &
     & .982,.996,.999,                                    &
     & 2.83E-8,7.68E-6,8.36E-5,1.71E-2,.133,.731,2.21,      &
     & 2.98,3.23,3.3,                                      &
     & 1.29E-8,4.54E-6,5.52E-5,1.46E-2,.126,.767,2.53,      &
     & 3.52,3.86,3.95/
      IF ( ig .eq. 1 ) gaml = 1.00    ! Gamma[1,0,x]
      IF ( ig .eq. 2 ) gaml = 1.00    ! Gamma[2,0,x]
      IF ( ig .eq. 3 ) gaml = 3.323351 ! Gamma[3.5,0,x]
      IF ( ig .eq. 4 ) gaml = 4.00 ! Gamma[11/3] = 4.012201 ! Gamma[3.667,0,x]
      IF ( x .gt. xg(10) ) RETURN
      DO ii = 1,9 
        i = 10 - ii
        n = i
        np = n + 1
        IF ( x .ge. xg(i) ) THEN
!         GOTO 2 
          gaml = gamxg(N,IG)+((X-XG(N))/(XG(NP)-XG(N)))*  &
     &            ( gamxg(NP,IG) - gamxg(N,IG) ) 
          RETURN
        ENDIF
      ENDDO
      gaml = (x/xg(1))*gamxg(1,ig)  ! linear approx. from zero to 0.01
      RETURN
      END 

!**************************** GAML02 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 40 micron drops)
! **********************************************************
      real FUNCTION GAML02(x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=12)
      real gamxg(ng), xg(ng)
      DATA xg/0.01,0.02,0.025,0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      DATA gamxg/  &
     &  7.391019203578037e-8,0.02212726874591478,0.06959352407989682, &
     &  0.2355654024970809,0.46135930387500346,0.545435791452399,     &
     &  0.7371571313308203,                                           &
     &  0.8265676632204345,0.8640182781845841,0.8855756211304151,     &
     &  0.9245079225301251,                                           &
     &  0.9712578342732681/
      IF ( x .ge. xg(ng) ) THEN
        gaml02 = xg(ng)
        RETURN
      ENDIF
      IF ( x .lt. xg(1) ) THEN
        gaml02 = 0.0
        RETURN
      ENDIF
      DO ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        IF ( x .ge. xg(i) ) THEN
!         GOTO 2 
          gaml02 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))* &
     &            ( gamxg(NP) - gamxg(N) ) 
          RETURN
        ENDIF
      ENDDO
      RETURN
      END 

!**************************** GAML02d300 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 300 micron drops) (see zieglerstuff.nb)
! **********************************************************
      real FUNCTION GAML02d300(x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=9)
      real gamxg(ng), xg(ng)
      DATA xg/0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      DATA gamxg/                           &
     &  0.0,                                  &
     &  7.391019203578011e-8,0.0002260640810600053,  &
     &  0.16567071824457152,                         &
     &  0.4231369044918005,0.5454357914523988,       &
     &  0.6170290936864555,                           &
     &  0.7471346054110058,0.9037156157718299 /
      IF ( x .ge. xg(ng) ) THEN
        GAML02d300 = xg(ng)
        RETURN
      ENDIF
      IF ( x .lt. xg(1) ) THEN
        GAML02d300 = 0.0
        RETURN
      ENDIF
      DO ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        IF ( x .ge. xg(i) ) THEN
!         GOTO 2 
          GAML02d300 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))*  &
     &            ( gamxg(NP) - gamxg(N) ) 
          RETURN
        ENDIF
      ENDDO
      RETURN
      END 
!c

! #####################################################################
! #####################################################################

!**************************** GAML02 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 500 micron drops) (see zieglerstuff.nb)
! **********************************************************
      real FUNCTION GAML02d500(x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=9)
      real gamxg(ng), xg(ng)
      DATA xg/0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      DATA gamxg/  &
     &  0.0,0.0,   &
     &  2.2346039e-13, 0.0221272687459,  &
     &  0.23556540,  0.38710348,         &
     &  0.48136183,0.6565833,            &
     &  0.86918315 /
      IF ( x .ge. xg(ng) ) THEN
        GAML02d500 = xg(ng)
        RETURN
      ENDIF
      IF ( x .lt. xg(1) ) THEN
        GAML02d500 = 0.0
        RETURN
      ENDIF
      DO ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        IF ( x .ge. xg(i) ) THEN
!         GOTO 2 
          GAML02d500 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))*  &
     &            ( gamxg(NP) - gamxg(N) ) 
          RETURN
        ENDIF
      ENDDO
      RETURN
      END 
!c

! #####################################################################
! #####################################################################

        subroutine GAMINC(A1,X1,ingam)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                r(a,x), ‚(a,x) and P(a,x)
!       Input :  a   --- Parameter ( a Û 170 )
!                x   --- Argument 
!       Output:  GIN --- g(a,x) t=0,x
!                GIM --- G(a,x) t=x,Infinity
!       Routine called: GAMMA for computing G(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1,gamma,ingam
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        double precision :: GAMMA_DP
        
        a = a1
        x = x1
!        print*,'a,x = ',a,x
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           GIM = GAMMA_DP(A)
        ELSE IF (X.LE.1.0+A) THEN
!        print*,'branch 1'
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           ga = GAMMA(A1)
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
!        print*,'branch 2, xam =',xam,dexp(xam)
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
!              print*,'k,t0 = ',k,t0
20         CONTINUE
!           print*,'x+t0 = ',x+t0
           GIM=DEXP(XAM)/(X+T0)
           GA = GAMMA(A1)
           GIN=GA-GIM
!           print*,'gim,ga,gin,a1 = ',gim,ga,gin,a1
        ENDIF
        
        ingam = GIN
        return
        END

! #####################################################################

        real function GAMZEROX(A1,X1)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                from 0 to x
!       Input :  a   --- Parameter ( a Û 170 )
!                x   --- Argument 
!       Output:  GIN --- r(a,x) t=0,x 
!       Routine called: GAMMA for computing ‚(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1,gamma
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        
        a = a1
        x = x1
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           GIM = GAMMA(A1)
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
!           ga = GAMMA(A1)
!           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           GA = GAMMA(A1)
           GIN=GA-GIM
        ENDIF
        
        gamzerox = GIN
        return
        END

! #####################################################################

        real function GAMXINF(A1,X1)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                from x to infinity
!       Input :  a   --- Parameter ( a Û 170 )
!                x   --- Argument 
!       Output:  GIM --- ‚(a,x) t=x,Infinity
!       Routine called: GAMMA for computing ‚(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1,gamma
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        
        a = a1
        x = x1
        IF ( x1 <= 0.0 ) THEN
           gamxinf = GAMMA(A1)
           return
        ENDIF
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           GIM = GAMMA(A1)
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           ga = GAMMA(A1)
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
!           GA = GAMMA(A1)
!           GIN=GA-GIM
        ENDIF
        
        gamxinf = GIM
        return
        END


! #####################################################################

        double precision function GAMXINFDP(A1,X1)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                from x to infinity
!       Input :  a   --- Parameter ( a Û 170 )
!                x   --- Argument 
!       Output:  GIM --- ‚(a,x) t=x,Infinity
!       Routine called: GAMMA for computing ‚(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1
        double precision :: gamma_dp
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        
        a = a1
        x = x1
        IF ( x1 <= 0.0 ) THEN
           gamxinfdp = GAMMA_DP(A)
           return
        ENDIF
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           GIM = GAMMA_dp(A)
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           ga = GAMMA_DP(A)
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
!           GA = GAMMA_dp(A)
!           GIN=GA-GIM
        ENDIF
        
        gamxinfdp = GIM
        return
        END


! #####################################################################

!       ==================================================

        real function GAMMAZJ(X1)

!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!       Code from "mincog.f" from the book Computation of Special Functions by
!        S. Zhang and J. M. Jin (John Wiley & Sons) 1996 (copyrighted)
        
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        real x1
        DIMENSION G(26)
        x = x1
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO K=2,M1
                 GA=GA*K
              ENDDO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO K=1,M
                R=R*(Z-K)
              ENDDO
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,                 &
     &          -0.6558780715202538D0, -0.420026350340952D-1, &
     &          0.1665386113822915D0,-.421977345555443D-1,    &
     &          -.96219715278770D-2, .72189432466630D-2,      &
     &          -.11651675918591D-2, -.2152416741149D-3,      &
     &          .1280502823882D-3, -.201348547807D-4,         &
     &          -.12504934821D-5, .11330272320D-5,            &
     &          -.2056338417D-6, .61160950D-8,                &
     &          .50020075D-8, -.11812746D-8,                  &
     &          .1043427D-9, .77823D-11,                      &
     &          -.36968D-11, .51D-12,                         &
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO K=25,1,-1
             GR=GR*Z+G(K)
           ENDDO
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        gammazj = ga
        RETURN
        END

! #####################################################################


        real function BETA(P,Q)
!
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA for computing ‚(x)
!       ==========================================
!
        IMPLICIT real (A-H,O-Z)
        double precision p1,gp,q1,gq, ppq,gpq
        
        p1 = p
        q1 = q
        CALL GAMMADP(P1,GP)
        CALL GAMMADP(Q1,GQ)
        PPQ=P1+Q1
        CALL GAMMADP(PPQ,GPQ)
        beta=GP*GQ/GPQ
        RETURN
        END

! #####################################################################

        SUBROUTINE GAMMADP(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO K=2,M1
                GA=GA*K
              ENDDO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO K=1,M
                R=R*(Z-K)
              ENDDO
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,                  &
     &          -0.6558780715202538D0, -0.420026350340952D-1,  &
     &          0.1665386113822915D0,-.421977345555443D-1,     &
     &          -.96219715278770D-2, .72189432466630D-2,       &
     &          -.11651675918591D-2, -.2152416741149D-3,       &
     &          .1280502823882D-3, -.201348547807D-4,          &
     &          -.12504934821D-5, .11330272320D-5,             &
     &          -.2056338417D-6, .61160950D-8,                 &
     &          .50020075D-8, -.11812746D-8,                   &
     &          .1043427D-9, .77823D-11,                       &
     &          -.36968D-11, .51D-12,                          &
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO K=25,1,-1
             GR=GR*Z+G(K)
           ENDDO
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
! #####################################################################
!
! #####################################################################
!--------------------------------------------------------------------------
      real function DERF(X)

      implicit none
     
      real :: X
      double precision, dimension(0 : 64) :: A, B
      double precision :: W,T,Y
      integer :: K,I
           data A/                                                 &
              0.00000000005958930743E0, -0.00000000113739022964E0, &
              0.00000001466005199839E0, -0.00000016350354461960E0, &
              0.00000164610044809620E0, -0.00001492559551950604E0, &
              0.00012055331122299265E0, -0.00085483269811296660E0, &
              0.00522397762482322257E0, -0.02686617064507733420E0, &
              0.11283791670954881569E0, -0.37612638903183748117E0, &
              1.12837916709551257377E0,                            &
              0.00000000002372510631E0, -0.00000000045493253732E0, &
              0.00000000590362766598E0, -0.00000006642090827576E0, &
              0.00000067595634268133E0, -0.00000621188515924000E0, &
              0.00005103883009709690E0, -0.00037015410692956173E0, &
              0.00233307631218880978E0, -0.01254988477182192210E0, &
              0.05657061146827041994E0, -0.21379664776456006580E0, &
              0.84270079294971486929E0,                            &
              0.00000000000949905026E0, -0.00000000018310229805E0, &
              0.00000000239463074000E0, -0.00000002721444369609E0, &
              0.00000028045522331686E0, -0.00000261830022482897E0, &
              0.00002195455056768781E0, -0.00016358986921372656E0, &
              0.00107052153564110318E0, -0.00608284718113590151E0, &
              0.02986978465246258244E0, -0.13055593046562267625E0, &
              0.67493323603965504676E0,                            &
              0.00000000000382722073E0, -0.00000000007421598602E0, &
              0.00000000097930574080E0, -0.00000001126008898854E0, &
              0.00000011775134830784E0, -0.00000111992758382650E0, &
              0.00000962023443095201E0, -0.00007404402135070773E0, &
              0.00050689993654144881E0, -0.00307553051439272889E0, &
              0.01668977892553165586E0, -0.08548534594781312114E0, &
              0.56909076642393639985E0,                            &
              0.00000000000155296588E0, -0.00000000003032205868E0, &
              0.00000000040424830707E0, -0.00000000471135111493E0, &
              0.00000005011915876293E0, -0.00000048722516178974E0, &
              0.00000430683284629395E0, -0.00003445026145385764E0, &
              0.00024879276133931664E0, -0.00162940941748079288E0, &
              0.00988786373932350462E0, -0.05962426839442303805E0, &
              0.49766113250947636708E0 /
           data (B(I), I = 0, 12) /                                 &
              -0.00000000029734388465E0,  0.00000000269776334046E0, &
              -0.00000000640788827665E0, -0.00000001667820132100E0, &
              -0.00000021854388148686E0,  0.00000266246030457984E0, &
               0.00001612722157047886E0, -0.00025616361025506629E0, &
               0.00015380842432375365E0,  0.00815533022524927908E0, &
              -0.01402283663896319337E0, -0.19746892495383021487E0, &
               0.71511720328842845913E0 /
           data (B(I), I = 13, 25) /                                &
              -0.00000000001951073787E0, -0.00000000032302692214E0, &
               0.00000000522461866919E0,  0.00000000342940918551E0, &
              -0.00000035772874310272E0,  0.00000019999935792654E0, &
               0.00002687044575042908E0, -0.00011843240273775776E0, &
              -0.00080991728956032271E0,  0.00661062970502241174E0, &
               0.00909530922354827295E0, -0.20160072778491013140E0, &
               0.51169696718727644908E0 /
           data (B(I), I = 26, 38) /                                &
              0.00000000003147682272E0, -0.00000000048465972408E0,  &
              0.00000000063675740242E0,  0.00000003377623323271E0,  &
             -0.00000015451139637086E0, -0.00000203340624738438E0,  &
              0.00001947204525295057E0,  0.00002854147231653228E0,  &
             -0.00101565063152200272E0,  0.00271187003520095655E0,  &
              0.02328095035422810727E0, -0.16725021123116877197E0,  &
              0.32490054966649436974E0 /
           data (B(I), I = 39, 51) /                                &
              0.00000000002319363370E0, -0.00000000006303206648E0,  &
             -0.00000000264888267434E0,  0.00000002050708040581E0,  &
              0.00000011371857327578E0, -0.00000211211337219663E0,  &
              0.00000368797328322935E0,  0.00009823686253424796E0,  &
             -0.00065860243990455368E0, -0.00075285814895230877E0,  &
              0.02585434424202960464E0, -0.11637092784486193258E0,  &
              0.18267336775296612024E0 /
           data (B(I), I = 52, 64) /                                &
             -0.00000000000367789363E0,  0.00000000020876046746E0,  &
             -0.00000000193319027226E0, -0.00000000435953392472E0,  &
              0.00000018006992266137E0, -0.00000078441223763969E0,  &
             -0.00000675407647949153E0,  0.00008428418334440096E0,  &
             -0.00017604388937031815E0, -0.00239729611435071610E0,  &
              0.02064129023876022970E0, -0.06905562880005864105E0,  &
              0.09084526782065478489E0 /
           W = ABS(X)
           if (W .LT. 2.2D0) then
               T = W * W
               K = INT(T)
               T = T - K
               K = K * 13
               Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
                   A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
                   A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
                   A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T +    &
                   A(K + 11)) * T + A(K + 12)) * W
           elseif (W .LT. 6.9D0) then
               K = INT(W)
               T = W - K
               K = 13 * (K - 2)
               Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
                   B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T +     &
                   B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T +     &
                   B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T +    &
                   B(K + 11)) * T + B(K + 12)
               Y = Y * Y
               Y = Y * Y
               Y = Y * Y
               Y = 1 - Y * Y
           else
               Y = 1
           endif
           if (X .LT. 0) Y = -Y
           DERF = Y
     
      end function DERF

! #####################################################################

!-------------------------------------------------------------------------
!
!  subroutine mlint1 linearly interpolates, in three dimensions,
!  for any variable 'ac'
!
      subroutine mlint2                             &
       (idebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,  &
        nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
!
      integer  idebug
      integer  ic
      integer  kc
      integer  jc
      integer  lc
!
      real     acx1
      real     acx2
      real     acx3
      real     acx4
      real     acy1
      real     acy2
      real     aint
      real     facx
      real     facy
      real     facz
      real     ac(n1b:n1e,n2b:n2e,n3b:n3e,n4b:n4e)
      
      
      acx1 = 0.0
      acx2 = 0.0
      acx3 = 0.0
      acx4 = 0.0
      acy1 = 0.0
      acy2 = 0.0
      
!
!  x-direction interpolation
!
      if( nx .gt. 1 ) then
!
      acx1 = facx*(ac(ic+1,jc,kc,lc)-ac(ic,jc,kc,lc)) +  &
             ac(ic,jc,kc,lc)
      acx2 = facx*(ac(ic+1,jc,kc+1,lc)-ac(ic,jc,kc+1,lc)) + &
             ac(ic,jc,kc+1,lc)
!
      else if( nx .le. 1 .and. ny .gt. 1 ) then
!
      acy1 = facy*(ac(ic,jc+1,kc,lc)-ac(ic,jc,kc,lc)) +       &
             ac(ic,jc,kc,lc)
      acy2 = facy*(ac(ic,jc+1,kc+1,lc)-ac(ic,jc+1,kc+1,lc)) + &
            ac(ic,jc+1,kc+1,lc)
!
      
      else  ! 1-d vertical interpolation

      acy1 = ac(ic,jc,kc,lc)
      acy2 = ac(ic,jc,kc+1,lc)
      
      end if
!
!  y-direction interpolation of data interpolated in x-direction
!  for bottem of cube if in 3-d; eg x:y:z trajectory
!
      if( nx .gt. 1 .and. ny .gt.1 ) then
!
      acx3 = facx*(ac(ic+1,jc+1,kc,lc)-ac(ic,jc+1,kc,lc)) +     &
             ac(ic,jc+1,kc,lc)
      acx4 = facx*(ac(ic+1,jc+1,kc+1,lc)-ac(ic,jc+1,kc+1,lc)) + &
             ac(ic,jc+1,kc+1,lc)
!
! y-direction interpolation of x interpolated data
!
      acy1 = facy*(acx3-acx1) + acx1
      acy2 = facy*(acx4-acx2) + acx2
!
!  z-direction interpolation of data interpolated in x:y-directions
!
      aint = facz*(acy2-acy1) + acy1
!
!  z-direction interpolation of data interpolated in y-direction for
!  2-d y:z trajectory
!
      else if ( nx .le. 1 ) then
      aint = facz*(acy2-acy1) + acy1
!
!  z-direction interpolation of data interpolated in x-direction for
!  2-d x:z trajectory
!
      else if ( ny .le. 1 ) then
      aint = facz*(acx2-acx1) + acx1
      end if
!
      if ( idebug .ge. 2 ) then
      IF ( nx .gt. 1 ) print*,'acx1,acx2: ',acx1, acx2, facx
      print*,'acx3,acx4: ',acx3, acx4, facx
      print*,'acy1,acy2: ',acy1, acy2, facy
      print*,'aint: ', aint, facz
      end if
!
      return
      end subroutine mlint2



!
!-------------------------------------------------------------------------
!  Borrowed from ARPS ADAS code

!      SUBROUTINE adaslcl(nx,ny,nz,zs,pr,pt,qv,hgtlcl)
      SUBROUTINE adaslcl(nx,ny,nz,gzt,pipert,piinit1d,theta,qv,hgtlcl)

      USE PARAM_MODULE, only: rcp,ng
      
      IMPLICIT NONE
      INTEGER :: nx,ny,nz
!      REAL :: zs(nx,ny,nz)
      REAL :: gzt(-ng+1:nz+ng,4)
!      REAL :: pr(nx,ny,nz)
      REAL :: theta(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      REAL :: pipert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: piinit1d(-ng+1:nz+ng)
      REAL :: qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      REAL :: hgtlcl(nx,ny)

      REAL :: pres(nz)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
      INTEGER :: i,j,k
      INTEGER :: imid,jmid
    
      REAL :: pmb,tc,tk,td,wmr,thepcl,plcl,tlcl
      REAL :: plnhi,plnlo,plnlcl,whi,wlo
    
      REAL :: oe,wmr2td

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Find pressure of lcl
!  Using pressure of lcl, find height of lcl
!
!-----------------------------------------------------------------------
!
        imid=nx/2
        jmid=ny/2
      
        DO j=1,ny-1
          DO i=1,nx-1
            pres(1) = (1.e5*(piinit1d(1) + pipert(i,j,1))**3.509)
            pres(2) = (1.e5*(piinit1d(2) + pipert(i,j,2))**3.509)
            pmb=0.01*pres(1) ! pr(i,j,2)
!            tk=pt(i,j,2)*((1000./pmb)**rddcp)
            tk=theta(i,j,1)*((1000./pmb)**rcp)
            tc=tk-273.15
!            wmr=1000.*(qv(i,j,1)/(1.-qv(i,j,1)))  ! why was this divided by (1 - qv)? wmr should just be qv, unless it is assuming qv is specific humidity (vapor mass per kg of moist air) to mixing ratio (per kg of dry air)?
            wmr=1000.*qv(i,j,1)
            td=wmr2td(pmb,wmr)
!            thepcl=oe(tc,td,pmb)
            CALL ptlcl(pmb,tc,td,plcl,tlcl)
            plcl=plcl*100.
            DO k=2,nz-2
              pres(k) = (1.e5*(piinit1d(k) + pipert(i,j,k))**3.509)
!              IF(pr(i,j,k) < plcl) EXIT
              IF(pres(k) < plcl) EXIT
            END DO
      !      81   CONTINUE
!            plnhi=LOG(pr(i,j,k))
!            plnlo=LOG(pr(i,j,k-1))
            plnhi=LOG(pres(k))
            plnlo=LOG(pres(k-1))
            plnlcl=LOG(plcl)
            whi=(plnlo-plnlcl)/(plnlo-plnhi)
            wlo=1.-whi
            hgtlcl(i,j)=whi*gzt(k,1)+wlo*gzt(k-1,1)
            IF(i == imid .AND. j == jmid) THEN
      !        PRINT *, 'arpslcl: ',pmb,tc,td,wmr,plcl,hgtlcl(i,j)
            END IF
          END DO
        END DO
      !
        RETURN
      END SUBROUTINE adaslcl

!-----------------------------------------------------------------------
!

      SUBROUTINE ptlcl(p,t,td,pc,tc)

      implicit none
      real, intent(IN) :: p,t,td
      real, intent(OUT) :: pc,tc
      real, parameter :: akap = 0.28541
      real, parameter :: cta  = 273.16
      
      real c1,c2
   
!
!   this subroutine estimates the pressure pc (mb) and the temperature
!   tc (celsius) at the lifted condensation level (lcl), given the
!   initial pressure p (mb), temperature t (celsius) and dew point
!   (celsius) of the parcel.  the approximation is that lines of
!   constant potential temperature and constant mixing ratio are
!   straight on the skew t/log p chart.
!
!    baker,schlatter   17-may-1982   original version
!
!   teten's formula for saturation vapor pressure as a function of
!   pressure was used in the derivation of the formula below.  for
!   additional details, see math notes by t. schlatter dated 8 sep 81.
!   t. schlatter, noaa/erl/profs program office, boulder, colorado,
!   wrote this subroutine.
!
!   akap = (gas constant for dry air) / (specific heat at constant
!       pressure for dry air)
!   cta = difference between kelvin and celsius temperatures
!
!  DATA akap,cta/0.28541,273.16/
         c1 = 4098.026/(td+237.3)**2
         c2 = 1./(akap*(t+cta))
         pc = p*EXP(c1*c2*(t-td)/(c2-c1))
         tc = t+c1*(t-td)/(c2-c1)
         RETURN
       END SUBROUTINE ptlcl

!-----------------------------------------------------------------------
!

  real FUNCTION oe(t,td,p)
  implicit none
  real t,td,p
  real tw,os,atw
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns equivalent potential temperature oe (celsius)
!   of a parcel of air given its temperature t (celsius), dew point
!   td (celsius) and pressure p (millibars).
!   find the wet bulb temperature of the parcel.

  atw = tw(t,td,p)

!   find the equivalent potential temperature.

  oe = os(atw,p)
  RETURN
  END FUNCTION oe

!-----------------------------------------------------------------------
!

  real FUNCTION os(t,p)
  implicit none 
  real t,p
  real tk,b,osk
  real w_adas
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the equivalent potential temperature os
!   (celsius) for a parcel of air saturated at temperature t (celsius)
!   and pressure p (millibars).
  DATA b/2.6518986/
!   b is an empirical constant approximately equal to the latent heat
!   of vaporization for water divided by the specific heat at constant
!   pressure for dry air.

  tk = t+273.15
  osk= tk*((1000./p)**.286)*(EXP(b*w_adas(t,p)/tk))
  os= osk-273.15
  RETURN
  END FUNCTION os

!-----------------------------------------------------------------------
!

  real FUNCTION tw(t,td,p)

!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the wet-bulb temperature tw (celsius)
!   given the temperature t (celsius), dew point td (celsius)
!   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
!   above, for a description of the technique.
!
!
!   determine the mixing ratio line thru td and p.

  aw = w_adas(td,p)
!
!   determine the dry adiabat thru t and p.

  ao = o(t,p)
  pi = p

!   iterate to locate pressure pi at the intersection of the two
!   curves .  pi has been set to p for the initial guess.

  DO i= 1,10
    x= .02*(tmr(aw,pi)-tda(ao,pi))
    IF (ABS(x) < 0.01) EXIT
    pi= pi*(2.**(x))
  END DO

!   find the temperature on the dry adiabat ao at pressure pi.

  ti= tda(ao,pi)

!   the intersection has been located...now, find a saturation
!   adiabat thru this point. function os returns the equivalent
!   potential temperature (c) of a parcel saturated at temperature
!   ti and pressure pi.

  aos= os(ti,pi)

!   function tsa returns the wet-bulb temperature (c) of a parcel at
!   pressure p whose equivalent potential temperature is aos.

  tw = tsa(aos,p)
  RETURN
  END FUNCTION tw

!-----------------------------------------------------------------------
!

  real FUNCTION w_adas(t,p)
!
!    g.s. stipanuk     1973              original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!  this function returns the mixing ratio (grams of water vapor per
!  kilogram of dry air) given the dew point (celsius) and pressure
!  (millibars). if the temperture  is input instead of the
!  dew point, then saturation mixing ratio (same units) is returned.
!  the formula is found in most meteorological texts.

  x= esat_fast(t)
  w_adas= 622.*x/(p-x)
  RETURN
  END FUNCTION w_adas

!-----------------------------------------------------------------------
!

  real FUNCTION esat_fast(t)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!   this function returns the saturation vapor pressure over
!   water (mb) given the temperature (celsius).
!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
!   tions of selected meteorlolgical parameters for cloud physics prob-
!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
!   electronics command, white sands missile range, new mexico 88002.

  tk = t+273.15
  p1 = 11.344-0.0303998*tk
  p2 = 3.49149-1302.8844/tk
  c1 = 23.832241-5.02808*ALOG10(tk)
  esat_fast = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tk)
  RETURN
  END FUNCTION esat_fast

!-----------------------------------------------------------------------
!

  real FUNCTION tda(o,p)
!
!    g.s. stipanuk     1973           original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns the temperature tda (celsius) on a dry adiabat
!   at pressure p (millibars). the dry adiabat is given by
!   potential temperature o (celsius). the computation is based on
!   poisson's equation.

  ok= o+273.15
  tdak= ok*((p*.001)**.286)
  tda= tdak-273.15
  RETURN
  END FUNCTION tda

!-----------------------------------------------------------------------
!
  real FUNCTION o(t,p)
!
!    g.s. stipanuk     1973          original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982

!   this function returns potential temperature (celsius) given
!   temperature t (celsius) and pressure p (mb) by solving the poisson
!   equation.

  tk= t+273.15
  ok= tk*((1000./p)**.286)
  o= ok-273.15
  RETURN
  END FUNCTION o

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION WMR2TD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

  real FUNCTION wmr2td(pres,wmr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate Dewpoint Temperature (TD) from vapor mixing ratio (WMR)
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  April, 1995   Based on GEMPAK routine of same name.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!  REAL :: wmr2td
  REAL :: pres
  REAL :: wmr
  REAL :: wkgkg,e,evap
!
  wkgkg = 0.001 * wmr
  wkgkg = AMAX1(wkgkg,0.00005)
  e= (pres*wkgkg) / (0.62197 + wkgkg)  ! from w = eps*e/(p - e)
  evap = e /(1.001 + (( pres - 100.) /900.) * 0.0034) ! Looks like some kind of interpolation to get e_sat from e and p.
  wmr2td = ALOG(evap/6.112) * 243.5 /( 17.67 - ALOG (evap/6.112)) ! Bolton (1980) formula, solve for T

  RETURN
  END FUNCTION wmr2td

!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                  FUNCTION RAN0                       ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

      FUNCTION ran0(idum)

!  ##################################################################
!  PURPOSE: to generate a random number.  To continue a sequence,
!           do not change idum between calls to ran0
! Minimal random number generator of Park and Miller. Returns a
! uniform random deviate between 0.0 and 1.0. Set or reset idum to
! any integer value (except the unlikely value MASK) to initialize
! the sequence; idum must not be altered between calls for
! successive deviates in a sequence.
!
!  AUTHOR:  lifted straight out of Numerical Recipes in Fortran 77
!
!  ##################################################################
      
      
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL ran0,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,  &
                 IQ=127773,IR=2836,MASK=123459876)
     
      INTEGER k
      
      idum=ieor(idum,MASK) 
!      XORing with MASK allows use of zero and other simple
!      bit patterns for idum. 
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k 
! Compute idum=mod(IA*idum,IM) without overflows by
!      Schrage's method. 
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum                !   Convert idum to afloating result.
      idum=ieor(idum,MASK)        !    Unmask before return.
      return
      END

!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                  FUNCTION ran0mpi                       ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

      FUNCTION ran0mpi(idum)

!  ##################################################################
!  PURPOSE: to generate a random number.  To continue a sequence,
!           do not change idum between calls to ran0mpi
! Minimal random number generator of Park and Miller. Returns a
! uniform random deviate between 0.0 and 1.0. Set or reset idum to
! any integer value (except the unlikely value MASK) to initialize
! the sequence; idum must not be altered between calls for
! successive deviates in a sequence.
!
!  AUTHOR:  lifted straight out of Numerical Recipes in Fortran 77
!
!  ##################################################################
      
      
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL ran0mpi,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,      &
     &           IQ=127773,IR=2836,MASK=123459876)
     
      INTEGER k
      
      idum=ieor(idum,MASK) 
!      XORing with MASK allows use of zero and other simple
!      bit patterns for idum. 
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k 
! Compute idum=mod(IA*idum,IM) without overflows by
!      Schrage's method. 
      if (idum.lt.0) idum=idum+IM
      ran0mpi=AM*idum                !   Convert idum to afloating result.
      idum=ieor(idum,MASK)        !    Unmask before return.
      return
      END



