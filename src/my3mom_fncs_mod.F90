module my3mom_fncs_mod

!==============================================================================!
!  The following functions are used the 'tmom' version of the multimoment      !
!  package.  Note: 'tmom' also uses functions in 'my_fncs_mod.cdk90'           !
!                                                                              !
!  Package version:  2.20.0      (internal bookkeeping)                        !
!  Last modified  :  2011-03-03                                                !
!==============================================================================!

   implicit none

   private
   public  :: NccnFNC_v33,SxFNC_v33,diagAlpha_v33,solveAlpha_v33,gammaDP

   contains

!==============================================================================!

 REAL FUNCTION NccnFNC_v33(Win,Tin,Pin,AIRTYPE)

!---------------------------------------------------------------------------!
! This function returns number concentration (activated aerosols) as a
! function of w,T,p, based on polynomial approximations of detailed
! approach using a hypergeometric function, following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

  IMPLICIT NONE

! PASSING PARAMETERS:
  real,    INTENT(IN) :: Win, Tin, Pin
  integer, INTENT(IN) :: AIRTYPE

! LOCAL PARAMETERS:
  real :: T,p,x,y,a,b,c,d,e,f,g,h,T2,T3,T4,x2,x3,x4,p2

  x= log10(Win*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin - 273.15;      T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;          p2= p*p

  if (AIRTYPE==1) then  !Maritime

     a= 1.47e-9*T4 -6.944e-8*T3 -9.933e-7*T2 +2.7278e-4*T -6.6853e-4
     b=-1.41e-8*T4 +6.662e-7*T3 +4.483e-6*T2 -2.0479e-3*T +4.0823e-2
     c= 5.12e-8*T4 -2.375e-6*T3 +4.268e-6*T2 +3.9681e-3*T -3.2356e-1
     d=-8.25e-8*T4 +3.629e-6*T3 -4.044e-5*T2 +2.1846e-3*T +9.1227e-1
     e= 5.02e-8*T4 -1.973e-6*T3 +3.944e-5*T2 -9.0734e-3*T +1.1256e0
     f= -1.424e-6*p2 +3.631e-3*p -1.986
     g= -0.0212*x4 +0.1765*x3 -0.3770*x2 -0.2200*x +1.0081
     h= 2.47e-6*T3 -3.654e-5*T2 +2.3327e-3*T +0.1938
     y= a*x4 + b*x3 + c*x2 + d*x + e + f*g*h
     NccnFNC_v33= 10.**min(2.,max(0.,y)) *1.e6                ![m-3]

  else if (AIRTYPE==2) then  !Continental

     a= 0.
     b= 0.
     c=-2.112e-9*T4 +3.9836e-8*T3 +2.3703e-6*T2 -1.4542e-4*T -0.0698
     d=-4.210e-8*T4 +5.5745e-7*T3 +1.8460e-5*T2 +9.6078e-4*T +0.7120
     e= 1.434e-7*T4 -1.6455e-6*T3 -4.3334e-5*T2 -7.6720e-3*T +1.0056
     f= 1.340e-6*p2 -3.5114e-3*p  +1.9453
     g= 4.226e-3*x4 -5.6012e-3*x3 -8.7846e-2*x2 +2.7435e-2*x +0.9932
     h= 5.811e-9*T4 +1.5589e-7*T3 -3.8623e-5*T2 +1.4471e-3*T +0.1496
     y= a*x4 +b*x3 +c*x2 + d*x + e + (f*g*h)
     NccnFNC_v33= 10.**max(0.,y) *1.e6

  else

    print*, '*** STOPPED in MODULE ### NccnFNC  *** '
    print*, '    Parameter AIRTYPE incorrectly specified'
    stop

  endif

 END FUNCTION NccnFNC_v33
!======================================================================!

   real*8 FUNCTION SxFNC_v33(Win,Tin,Pin,Qsw,Qsi,AIRTYPE,WRT)

!---------------------------------------------------------------------------!
! This function returns the peak supersaturation achieved during ascent with
! activation of CCN aerosols as a function of w,T,p, based on polynomial
! approximations of detailed approach using a hypergeometric function,
! following Cohard and Pinty (2000a).
!---------------------------------------------------------------------------!

 IMPLICIT NONE

! PASSING PARAMETERS:
  integer, INTENT(IN) :: WRT
  integer, INTENT(IN) :: AIRTYPE
  real,    INTENT(IN) :: Win, Tin, Pin, Qsw, Qsi

! LOCAL PARAMETERS:
  real   ::  FOQSA,FOQST,Si,Sw,Qv,T,p,x,a,b,c,d,f,g,h,Pcorr,T2corr,   &
             T2,T3,T4,x2,x3,x4,p2
  real, parameter :: TRPL= 273.15

  x= log10(max(Win,1.e-20)*100.);   x2= x*x;  x3= x2*x;  x4= x2*x2
  T= Tin;                           T2= T*T;  T3= T2*T;  T4= T2*T2
  p= Pin*0.01;                      p2= p*p

  if (AIRTYPE==1) then  !Maritime

     a= -5.109e-7*T4 -3.996e-5*T3 -1.066e-3*T2 -1.273e-2*T +0.0659
     b=  2.014e-6*T4 +1.583e-4*T3 +4.356e-3*T2 +4.943e-2*T -0.1538
     c= -2.037e-6*T4 -1.625e-4*T3 -4.541e-3*T2 -5.118e-2*T +0.1428
     d=  3.812e-7*T4 +3.065e-5*T3 +8.795e-4*T2 +9.440e-3*T +6.14e-3
     f= -2.012e-6*p2 + 4.1913e-3*p    - 1.785e0
     g=  2.832e-1*x3 -5.6990e-1*x2 +5.1105e-1*x -4.1747e-4
     h=  1.173e-6*T3 +3.2174e-5*T2 -6.8832e-4*T +6.7888e-2
     Pcorr= f*g*h
     T2corr= 0.9581-4.449e-3*T-2.016e-4*T2-3.307e-6*T3-1.725e-8*T4

  else if (AIRTYPE==2) then  !Continental [computed for -35<T<-5C]

     a=  3.80e-5*T2 +1.65e-4*T +9.88e-2
     b= -7.38e-5*T2 -2.53e-3*T -3.23e-1
     c=  8.39e-5*T2 +3.96e-3*T +3.50e-1
     d= -1.88e-6*T2 -1.33e-3*T -3.73e-2
     f= -1.9761e-6*p2 + 4.1473e-3*p - 1.771e0
     g=  0.1539*x4 -0.5575*x3 +0.9262*x2 -0.3498*x -0.1293
     h=-8.035e-9*T4+3.162e-7*T3+1.029e-5*T2-5.931e-4*T+5.62e-2
     Pcorr= f*g*h
     T2corr= 0.98888-5.0525e-4*T-1.7598e-5*T2-8.3308e-8*T3

  else

    print*, '*** STOPPED in MODULE ### SxFNC  *** '
    print*, '    Parameter AIRTYPE incorrectly specified'
    stop

  endif

  Sw= (a*x3 + b*x2 +c*x + d) + Pcorr
  Sw= 1. + 0.01*Sw
  Qv= Qsw*Sw
  Si= Qv/Qsi
  Si= Si*T2corr
  if (WRT.eq.1) then
     SxFNC_v33= Sw
  else
     SxFNC_v33= Si
  endif
  if (Win.le.0.) SxFNC_v33= 1.

 END function SxFNC_v33
!======================================================================!

 FUNCTION diagAlpha_v33(Dm,x)

  IMPLICIT NONE

  integer :: x
  real*8  :: diagAlpha_v33,Dm
  real*8, dimension(5) :: c1,c2,c3,c4
  real*8, parameter    :: pi = 3.14159265d0
  real*8, parameter  :: alphaMAX= 80.d0
  data c1 /19.0d0, 12.0d0, 4.5d0, 5.5d0, 3.7d0/
  data c2 / 0.6d0,  0.7d0, 0.5d0, 0.7d0, 0.3d0/
  data c3 / 1.8d0,  1.7d0, 5.0d0, 4.5d0, 9.0d0/
  data c4 /17.0d0, 11.0d0, 5.5d0, 8.5d0, 6.5d0/
  diagAlpha_v33= c1(x)*tanh(c2(x)*(1.d3*Dm-c3(x)))+c4(x)
  if (x==5.and.Dm>0.008d0) diagAlpha_v33= 1.d3*Dm-2.6d0
  diagAlpha_v33= min(diagAlpha_v33, alphaMAX)

 END function diagAlpha_v33
!======================================================================!

 FUNCTION solveAlpha_v33(Q,N,Z,Cx,rho)

 IMPLICIT NONE

! PASSING PARAMETERS:
  real, INTENT(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real*8 :: solveAlpha_v33
  real   :: a,g,a1,g1,g2,tmp1
  integer :: i
  real, parameter :: alphaMax= 40.
  real, parameter :: epsQ    = 1.e-14
  real, parameter :: epsN    = 1.e-3
  real, parameter :: epsZ    = 1.e-32

!  Q         mass mixing ratio
!  N         total concentration
!  Z         reflectivity
!  Cx        (pi/6)*RHOx
!  rho       air density
!  a         alpha (returned as solveAlpha)
!  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
!              where g = (Cx/(rho*Q))**2.*(Z*N)


  if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
  ! For testing/debugging only; this module should never be called
  ! if the above condition is true.
    print*,'*** STOPPED in MODULE ### solveAlpha *** '
    print*,'*** : ',Q,N,Z,Cx*1.9099,rho
    stop
  endif

  IF (Q>epsQ .and. N>epsN .and. Z>epsZ ) THEN

     tmp1= Cx/(rho*Q)
     g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2

 !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force; for testing only)
!      a= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            a = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of g(a) to solve for a:  [2004-11-29]
     if (g>=20.) then
       a= 0.
     else
       g2= g*g
       if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
       if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
       if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
       if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
       if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
       if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
       if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
       if (g<1.230) a= alphaMax
     endif

     solveAlpha_v33= max(0.,min(a,alphaMax))

  ELSE

     solveAlpha_v33= 0.

  ENDIF

 END FUNCTION solveAlpha_v33
!======================================================================!

 FUNCTION gammaDP(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gammaDP
  INTEGER  :: j
  DOUBLE PRECISION  :: ser,stp,tmp,x,y,cof(6)


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammaDP=tmp+log(stp*ser/x)
  gammaDP= exp(gammaDP)

 END FUNCTION gammaDP
!======================================================================!

end module my3mom_fncs_mod
