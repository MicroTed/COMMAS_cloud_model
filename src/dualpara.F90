#ifdef COMMAS
#define IJSTAG 1
#else
#define IJSTAG 0
#endif
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module dualpara                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE DUALPARA

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declare some constants used for calculaion of dual polarization
! parameters such as Zhh, Zdr, and Kdp. (It can be expanded...)
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
! MODIFIED: Dan Dawson, 05/26/2015
!           General code cleanup. Removed unneeded/unused functions and subroutines.
!
!-----------------------------------------------------------------------
! Declare parameters.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

  REAL, PARAMETER :: pi = 3.141592   ! pi

  REAL :: lambda           ! wavelength of radar (mm)

  REAL,PARAMETER :: Kw2 = 0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: alphaa = 4.28e-4   ! backscattering amplitude constant
                                       ! along major axis for rain
  REAL,PARAMETER :: beta_ra = 3.04
  REAL,PARAMETER :: alphab = 4.28e-4   ! backscattering amplitude constant
                                       ! along minor axis for rain
  REAL,PARAMETER :: beta_rb = 2.77
  REAL,PARAMETER :: alphak = 3.88e-4   ! differential forward scattering
                                       ! amplitude for rain
  REAL,PARAMETER :: alphask = 8.53e-7   ! differential forward scattering
                                        ! amplitude for snow
  REAL,PARAMETER :: alphaa_ds = 1.94e-5 ! for dry snow at horz plane
  REAL,PARAMETER :: alphab_ds = 1.91e-5 ! for dry snow at vert plane
  REAL,PARAMETER :: alphaa_dh = 1.91e-4 ! for dry hail at horz plane
  REAL,PARAMETER :: alphab_dh = 1.65e-4 ! for dry hail at vert plane
  REAL,PARAMETER :: alphaa_dg = 0.81e-4 ! for dry graupel at horz plane
  REAL,PARAMETER :: alphab_dg = 0.76e-4 ! for dry graupel at vert plane

  REAL,PARAMETER :: alphak_ds = 0.03e-5 ! alphaa_ds - alphab_ds
  REAL,PARAMETER :: alphak_dh = 0.26e-4 ! alphaa_dh - alphab_dh
  REAL,PARAMETER :: alphak_dg = 0.05e-4 ! alphaa_dh - alphab_dh

  REAL,PARAMETER :: rho_0r = 1.0      ! rho_0 for rain
  REAL,PARAMETER :: rho_0s = 1.0      ! rho_0 for snow
  REAL,PARAMETER :: rho_0h = 0.97     ! rho_0 for hail
  REAL,PARAMETER :: rho_0g = 0.40     ! rho_0 for graupel
  REAL,PARAMETER :: rho_0rsi = 0.82   ! lower limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rsf = 0.95   ! upper limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rhi = 0.85   ! lower limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rhf = 0.95   ! upper limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rgi = 0.82   ! lower limit of rho_0rg (rain-graupel mixture)
  REAL,PARAMETER :: rho_0rgf = 0.95   ! upper limit of rho_0rg (rain-graupel mixture)

  REAL,PARAMETER :: Coef = 5.6212976E+26
                    ! 4.*(lambda)**4.*gamma(7)/(pi**4.*Kw2)*(6*10**x/(pi*gamma(4)))^1.75
  REAL,PARAMETER :: kdpCoef = 1.1709e10    ! 180*lambda*1.e6*6/(pi**2)

  REAL,PARAMETER :: degKtoC=273.15 ! Conversion factor from degrees K to
                                   !   degrees C

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)

  REAL,PARAMETER :: mm3todBZ=1.0E+9 ! Conversion factor from mm**3 to
                                    !   mm**6 m**-3.

! The following parameters are for the Thompson microphysics scheme bimodal snow distribution
  REAL,PARAMETER :: tN0s1 = 490.6 ! 1st snow intercept parameter (k0)
  REAL,PARAMETER :: tN0s2 = 17.46 ! 2nd snow intercept parameter (k1)
  REAL,PARAMETER :: tlamdas1 = 20.78 ! 1st snow slope parameter
  REAL,PARAMETER :: tlamdas2 = 3.29 ! 2nd snow slope parameter
  REAL,PARAMETER :: talphas2 = 0.6357 ! 2nd snow shape parameter (1st is always 0)

! The following parameters are for the KMA implementation of the UM microphysics of snow
! which assumes a bimodal distribution with two sets of intercept, slope, and shape parameters
! The values are from the midlatitude option in Table 2 of Field et al. (2007)
  REAL,PARAMETER :: umN0s1 = 141.     ! 1st snow intercept parameter
  REAL,PARAMETER :: umN0s2 = 102.     ! 2nd snow intercept parameter
  REAL,PARAMETER :: umalphas2 = 2.07   ! 2nd snow shape parameter (1st is always 0)
  REAL,PARAMETER :: umlamdas1 = 16.8  ! 1st snow slope parameter
  REAL,PARAMETER :: umlamdas2 = 4.82  ! 2nd snow slope parameter

! The following parameters are for the KMA implementation of the UM microphysics for 
! rain and graupel diagnostic N0

  REAL,PARAMETER :: nar = 0.22
  REAL,PARAMETER :: nbr = 2.2
  REAL,PARAMETER :: nag = 5e25
  REAL,PARAMETER :: nbg = -4.0

  REAL :: tmpN0s1,tmpN0s2,tmplamdas1,tmplamdas2,tmpalphas1,tmpalphas2
  REAL :: N0s1,N0s2,alphas1,alphas2,lamdas1,lamdas2

! Missing value
  REAL,PARAMETER :: missing = -9999.0
  REAL :: grpl_miss
  REAL :: hl_miss 
  
  LOGICAL :: firstcall = .true.
  INTEGER :: grpl_ON, snow_ON
  INTEGER :: hl_ON 
  INTEGER :: qgh_opt

!-----------------------------------------------------------------------
! Precalculated complete gamma function values
!-----------------------------------------------------------------------
  REAL,PARAMETER :: gamma7_08 = 836.7818
  REAL,PARAMETER :: gamma6_81 = 505.8403
  REAL,PARAMETER :: gamma6_54 = 309.3308
  REAL,PARAMETER :: gamma5_63 = 64.6460
  REAL,PARAMETER :: gamma4_16 = 7.3619
  REAL,PARAMETER :: gamma3_97 = 5.7788

!-----------------------------------------------------------------------
! Variables to can be changed by parameter retrieval
!-----------------------------------------------------------------------
  REAL :: N0r        ! Intercept parameter in 1/(m^4) for rain
  REAL :: N0s        ! Intercept parameter in 1/(m^4) for snow
  REAL :: N0h        ! Intercept parameter in 1/(m^4) for hail
  REAL :: N0g        ! Intercept parameter in 1/(m^4) for graupel

  REAL :: rhor=1000. ! Density of rain (kg m**-3)
  REAL :: rhoh       ! Density of hail (kg m**-3)
  REAL :: rhos       ! Density of snow (kg m**-3)
  REAL :: rhog       ! Density of graupel (kg m**-3)

!-----------------------------------------------------------------------
! Variables to can be changed for meling ice
!-----------------------------------------------------------------------
  REAL :: fos        ! Maximum fraction of rain-snow mixture
  REAL :: foh        ! Maximum fraction of rain-hail mixture
  REAL :: fog        ! Maximum fraction of rain-graupel mixture

!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------

  REAL :: ZerhfGamma, ZervfGamma, ZerhvfGamma
  REAL :: constKdpr

!-----------------------------------------------------------------------
! Scattering matrix coefficient for snow
!
! phi=0.       (Mean orientation)
! sigmas=pi/9
! As=1/8*(3+4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Bs=1/8*(3-4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Cs=1/8*(1-cos(4*phi)*exp(-8*sigmas**2))
! Ds=1/8*(3+cos(4*phi)*exp(-8*sigmas**2))
! Cks=cos(2*phi)*exp(-2*sigmas**2)
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmas = 0.3491
  REAL,PARAMETER :: As = 0.8140
  REAL,PARAMETER :: Bs = 0.0303
  REAL,PARAMETER :: Cs = 0.0778
  REAL,PARAMETER :: Ds = 0.4221
  REAL,PARAMETER :: Cks = 0.7837

!-----------------------------------------------------------------------
! Scattering matrix coefficient for hail
!
! phi=0.     (Mean orientation)
! sigmah=pi/3*(1-sf*fw)
! Ah=1/8*(3+4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Bh=1/8*(3-4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Ch=1/8*(1-cos(4*phi)*exp(-8*sigmah**2))
! Dh=1/8*(3+cos(4*phi)*exp(-8*sigmah**2))
! Ckh=cos(2*phi)*exp(-2*sigmah**2)
!
! corresponding coefficient for dry hail: Ahd, Bhd, Chd, Dhd, Ckhd
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmahd = 1.0472
  REAL,PARAMETER :: Ahd = 0.4308
  REAL,PARAMETER :: Bhd = 0.3192
  REAL,PARAMETER :: Chd = 0.1250
  REAL,PARAMETER :: Dhd = 0.3750
  REAL,PARAMETER :: Ckhd = 0.1116

  REAL,PARAMETER :: q_threshold = 2.e-4
  REAL :: sf
  REAL :: sigmah, Ah, Bh, Ch, Dh, Ckh

!-----------------------------------------------------------------------
! Scattering matrix coefficient for graupel
! 
! phi=0.     (Mean orientation)
! sigmag=pi/3*(1-sf*fw)
! Ag=1/8*(3+4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Bg=1/8*(3-4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Cg=1/8*(1-cos(4*phi)*exp(-8*sigmag**2))
! Dg=1/8*(3+cos(4*phi)*exp(-8*sigmag**2))
! Ckg=cos(2*phi)*exp(-2*sigmag**2)
! 
! corresponding coefficient for dry graupel: Agd, Bgd, Cgd, Dgd, Ckgd
!-----------------------------------------------------------------------
  
  REAL,PARAMETER :: sigmagd = 1.0472
  REAL,PARAMETER :: Agd = 0.4308
  REAL,PARAMETER :: Bgd = 0.3192
  REAL,PARAMETER :: Cgd = 0.1250
  REAL,PARAMETER :: Dgd = 0.3750
  REAL,PARAMETER :: Ckgd = 0.1116
  
  REAL :: sigmag, Ag, Bg, Cg, Dg, Ckg

#ifndef COMMAS
!-----------------------------------------------------------------------
!  Declare new observation type
!-----------------------------------------------------------------------

! TYPE T_obs_dual
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw
#endif
!-----------------------------------------------------------------------
!  Parameters related to water fraction diagnostics
!-----------------------------------------------------------------------

  REAL :: wgfac = 1.0 ! Reduction factor for water shell for wet growth conditions (valid for MFflg = 1)
  LOGICAL :: tosoak   = .true. ! Soak water into graupel and hail to raise density to 910 before adding to shell?
  LOGICAL :: freezesoak=.true. ! Assume soaked water is frozen?  Recommended true for now since fractional water lookup table
                    ! assumes axis ratios are as if the water is all on the shell
  LOGICAL :: conserveice=.true. ! Conserve ice mass during shell iteration? (If false, liquid water mass is still conserved)
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SUBROUTINES AND FUNCTIONS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  CONTAINS

  SUBROUTINE setgrplhl(graupel_ON, hail_ON, wetsnow_ON)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set graupel and hail parameters
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 02/11/2011
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: graupel_ON, hail_ON, wetsnow_ON

  grpl_ON = graupel_ON
  hl_ON = hail_ON
  snow_ON = wetsnow_ON

  RETURN

  END SUBROUTINE setgrplhl

  SUBROUTINE init_fox()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can vary depend on whether graupel/hail exists. 
!-----------------------------------------------------------------------
  fos = 0.35             ! Maximum fraction of rain-snow mixture
  foh = 0.2              ! Maximum fraction of rain-hail mixture
  fog = 0.25             ! Maximum fraction of rain-graupel mixture

  END SUBROUTINE init_fox

  SUBROUTINE init_fox_no_grpl()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice 
! when graupel is suppressed.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  fos = 0.5              ! Maximum fraction of rain-snow mixture
  foh = 0.3              ! Maximum fraction of rain-hail mixture
  fog = 0.0              ! Maximum fraction of rain-graupel mixture
      
  END SUBROUTINE init_fox_no_grpl

  SUBROUTINE init_fox_no_hail() 

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
!  Setup default maximum fraction of water in the melting ice 
!  when hail is suppressed. 
!
!-----------------------------------------------------------------------
!
! AUTHOR: Bryan Putnam, 12/14/10
!
!-----------------------------------------------------------------------
! Force explicit declarations. 
!-----------------------------------------------------------------------

  IMPLICIT NONE 

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval 
!-----------------------------------------------------------------------

  fos = 0.5             ! Maximum fraction of rain-snow mixture
  foh = 0.0             ! Maximum fraction of rain-hail mixture
  fog = 0.4             ! Maximum fraction of rain-graupel mixture

  END SUBROUTINE init_fox_no_hail

  SUBROUTINE init_dsd()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default dsd values or reinialize default dsd values
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  N0r=8.0E+06 ! Intercept parameter in 1/(m^4) for rain.
  N0h=4.0E+04 ! Intercept parameter in 1/(m^4) for hail.
  N0s=3.0E+06 ! Intercept parameter in 1/(m^4) for snow.
  N0g=4.0E+05 ! Intercept parameter in 1/(m^4) for graupel.

  rhoh=913.  ! Density of hail (kg m**-3)
  rhos=100.  ! Density of snow (kg m**-3)
  rhog=400.  ! Density of graupel (kg m**-3)

  END SUBROUTINE init_dsd

  SUBROUTINE model_dsd(MPflg,n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set dsd values to those used in the arps forecasts
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  INTEGER :: MPflg
  REAL :: n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r=n0rain
  N0s=n0snow
  N0h=n0hail
  N0g=n0grpl

  rhos=rhosnow
  rhoh=rhohail
  rhog=rhogrpl

  ! Set some additional parameters based on choice of particular microphysics scheme
  SELECT CASE (MPflg)
  CASE(3) ! UM scheme
    tmpN0s1 = umN0s1
    tmpN0s2 = umN0s2
    tmplamdas1 = umlamdas1
    tmplamdas2 = umlamdas2
    tmpalphas2 = umalphas2
  CASE(4) ! Thompson scheme
    tmpN0s1 = tN0s1
    tmpN0s2 = tN0s2
    tmplamdas1 = tlamdas1
    tmplamdas2 = tlamdas2
    tmpalphas2 = talphas2
  END SELECT

  END SUBROUTINE model_dsd

  SUBROUTINE coeff_hail(fw,qml)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for hail
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: fw, qml

! DTD: test 02/16/2012, do not adjust coefficient downward for small hail
!  IF(qml < q_threshold) THEN
!     sf = 4*qml*1.e3
!  ELSE
     !sf = 0.6
     !sf = 0.8
     sf = 1.0
     !sf = 1.2
!  ENDIF
!
!  !sigmah=2.*pi/9.*(1-sf*fw) ! Was pi/3.
!	IF(fw > 0.5) THEN
!		sigmah = 0.
!	ELSEIF(fw > 0.1) THEN
!		sigmah = pi/3.*(1.-(1/.4)*sf*(fw-0.1))
!	ELSE
!		sigmah = pi/3.
!	ENDIF
	
! The following curve is used for the ZDR paper:
	
  IF(fw > 0.5) THEN
    sigmah = 0.
  ELSE
    sigmah = pi/3.*(1.-2.*sf*fw)
  ENDIF
	
	
  !sigmah=pi/3.*max(0.4,1.-sf*fw)
  !sigmah = 0.0
  Ah=.125*(3+4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Bh=.125*(3-4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Ch=.125*(1-exp(-8*sigmah**2))
  Dh=.125*(3+exp(-8*sigmah**2))
  Ckh=exp(-2*sigmah**2)

  END SUBROUTINE coeff_hail
    
  SUBROUTINE coeff_grpl(fw,qml)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for graupel
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
! MODIFIED: Dan Dawson, 02/16/2012
!           Made separate version for graupel.
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: fw, qml

! DTD: test 02/16/2012, do not adjust coefficient downward for small hail
!  IF(qml < q_threshold) THEN
!     sf = 4*qml*1.e3
!  ELSE
     !sf = 0.6
     !sf = 0.8
     sf = 1.0
     !sf = 1.2
!  ENDIF
    
!  sigmag=pi/3.*(1-sf*fw)
!	IF(fw > 0.5) THEN
!		sigmag = 0. ! pi/24. ! 0.
!	ELSEIF(fw > 0.1) THEN
!		sigmag = pi/3.*(1.-(1/.4)*sf*(fw-0.1))! pi/3. + (0.-pi/3.)*(fw-0.1)/(0.5-0.1) ! pi/3.*(1.-(1/.4)*sf*(fw-0.1))
!	ELSE
!		sigmag = pi/3.
!	ENDIF
	
! The following curve is used for the ZDR paper:
	
  IF(fw > 0.5) THEN
    sigmag = 0.
  ELSE
    sigmag = pi/3.*(1.-2.*sf*fw)
  ENDIF

	
  !sigmag=pi/3.*max(0.4,1.-sf*fw)
  Ag=.125*(3+4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
  Bg=.125*(3-4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
  Cg=.125*(1-exp(-8*sigmag**2))
  Dg=.125*(3+exp(-8*sigmag**2))
  Ckg=exp(-2*sigmag**2)

  END SUBROUTINE coeff_grpl

SUBROUTINE cal_N0(rhoa,q,Ntx,rhox,alpha,N0,mu)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  (03/26/2008)
!  Recast N0 as a double precision variable, and used double precision for
!  all intermediate calculations.  The calling subroutine should
!  also define it as double precision.  For situations with large alpha,
!  N0 can become very large, and loss of precision can result.
!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
!  With Jason Milbrandt's calculation of N0 just before evaporation in
!  the multi-moment code.
!
!  (02/14/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q,Ntx
  REAL*8 :: N0,alpha,mu
!JYS  REAL*8 :: N0_eff
  REAL :: rhox
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma_r8

  DOUBLE PRECISION :: lamda

  gamma1 = gamma_r8((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma_r8((4.d0+alpha)/(3.d0*mu))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
        dble(q)))**mu
  ELSE
    lamda = 0.d0
  END IF

  N0 = 3*mu*dble(Ntx)*lamda**(0.5d0*((1.d0+alpha)/(3*mu)))*                         &
              (1.d0/gamma1)*lamda**(0.5d0*((1.d0+alpha))/(3*mu))

!JYS  IF(lamda /= 0.d0) THEN
!JYS    N0_eff = N0*(((4.d0+alpha)/lamda)**alpha)*gamma4*(128.d0/3.d0)/   &
!JYS                ((4.d0+alpha)**(4.d0+alpha))
!JYS  ELSE
!JYS    N0_eff = 0.d0
!JYS  END IF

END SUBROUTINE cal_N0

SUBROUTINE cal_lamda(rhoa,q,Ntx,rhox,alpha,lamda,mu)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter lamda
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Converted intermediate calculations and arrays alpha and lamda to
!  double precision.
!
!  (02/14/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q
  REAL*8 :: lamda,alpha,mu
  REAL :: rhox
  REAL :: Ntx
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma_r8

  gamma1 = gamma_r8((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma_r8((4.d0+alpha)/(3.d0*mu))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
          dble(q)))**mu
  ELSE
    lamda = 0.d0
  END IF

END SUBROUTINE cal_lamda

SUBROUTINE cal_Nt(rhoa,q,N0,cx,alpha,Ntx,mu)

!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration at scalar points
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/06/2008) 
!   
!  MODIFICATION HISTORY:
!  
!  03/31/08 - converted intermediate calculations to double precision
!             as well as a few of the input arguments.
!
!  (02/15/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
! 
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  REAL :: rhoa,q
  REAL*8 :: N0,alpha,mu
  REAL :: cx
  REAL :: Ntx
  REAL*8 :: gamma1,gamma4

  REAL*8 :: gamma_r8

  gamma1 = gamma_r8((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma_r8((4.d0+alpha)/(3.d0*mu))

  Ntx = sngl((N0*gamma1/(3.d0*mu))**(3.d0/(4.d0+alpha))*   &
                   ((gamma1/gamma4)*dble(rhoa)* &
                   dble(q)/dble(cx))**((1.d0+alpha)/(4.d0+alpha)))

END SUBROUTINE cal_Nt

SUBROUTINE power_mom(MPflg,power,cx,t_c,q,moment)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates moments of the PSD based on the Field et al. 2005 power law
! relations. Used for Thompson scheme.
!
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
! MODIFIED: Dan Dawson, 08/19/2014 for use in python version.  Switched
!           to Field et al. 2007 version for use with "global version"
!           of UM microphysics for aggregates in the KMA model.
!
!           Dan Dawson, 09/18/2014.  Added MPflg to switch between Field et al. 
!           (2005) for Thompson scheme and Field et al. (2007) for UM scheme.
!
!-----------------------------------------------------------------------

!  USE DUALPARA

  IMPLICIT NONE

!  INCLUDE 'globcst.inc'
!  INCLUDE 'phycst.inc'
!  INCLUDE 'arpsenkf.inc'

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------


  INTEGER, INTENT(IN) :: MPflg,power
  REAL, INTENT(IN) :: t_c,q,cx
  REAL, INTENT(OUT) :: moment

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: a,b,d,e,f
  REAL*8 :: log_a
  REAL :: second_moment,test

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  second_moment = (q/cx)

  IF(power == 2) THEN
    moment = second_moment
  ELSE
    SELECT CASE (MPflg)
    CASE(3) ! UM scheme
      d = exp(13.6 - 7.76*power + 0.479*power**2.)
      e = -0.0361 + 0.0151*power + 0.00149*power**2.
      f = 0.807 + 0.00581*power + 0.0457*power**2.
     
      moment = d*exp(e*t_c)*second_moment**f
    CASE(4) ! Thompson scheme
      log_a = dble(5.065339-.062659*t_c - 3.032362*power +                 &
                   0.029469*t_c*power -  &
        0.000285*(t_c**2) + 0.312550*(power**2) + 0.0000204*(t_c**2)*power + &
        0.003199*t_c*(power**2) + 0.000000*(t_c**3) - 0.015952*(power**3))

      a = sngl(10.d0**log_a)

      b = 0.476221 - 0.015896*t_c + 0.165977*power + 0.007468*t_c*power -   &
        0.000141*(t_c**2) + 0.060366*(power**2) + 0.000079*(t_c**2)*power + &
        0.000594*t_c*(power**2) + 0.000000*(t_c**3) - 0.003577*(power**3)
        
      moment = a*(second_moment)**b
    END SELECT
  END IF

END SUBROUTINE

SUBROUTINE cal_D0(rhoa,q,Ntx,rhox,alpha,D0)
  !   
  !-----------------------------------------------------------------------
  !  PURPOSE:  Calculates median volume diameter (m) for a gamma distribution
  !            Note, assumes spherical particles.
  !-----------------------------------------------------------------------
  !
  !  AUTHOR: Dan Dawson 
  !  (02/06/2008) 
  !   
  !  MODIFICATION HISTORY:
  !   
  !-----------------------------------------------------------------------
  !  Variable Declarations:
  !-----------------------------------------------------------------------
  !    
  
  REAL :: rhoa,q,Ntx,rhox,D0
  REAL*8 :: alpha,lamda,mu
  
  mu = 1./3. ! Assumes standard gamma distribution!
  
  ! Compute lamda
  CALL cal_lamda(rhoa,q,Ntx,rhox,alpha,lamda,mu)
  IF(lamda > 0) THEN
    D0 = (3.672 + alpha)/lamda
  ELSE
    D0 = 0.0
  ENDIF

END SUBROUTINE cal_D0

SUBROUTINE calN0g_Thompson(tair,rhoa,qg,cg,alphag,rhog,qr,mvd_r,N0g)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter for graupel in the 
!            Thompson scheme. Based on code from module_mp_thompson.F 
!            in WRFV3.
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/04/2015) 
!   
!  MODIFICATION HISTORY:
!  
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     
  IMPLICIT NONE

  REAL :: tair,rhoa,qg,cg,rhog,qr,mvd_r
  REAL :: gonv_min, gonv_max
  REAL*8 :: alphag, N0_min, N0_exp, N0g, lam_exp, lamg, ilamg
  REAL :: D0r,xslw1,rg,ygra1,zans1,am_g,bm_g,oge1,cgg_1,cgg_2,cgg_3,ogg1,ogg2,obmg
  REAL :: cge_2

  gonv_min = 1.e4
  gonv_max = 3.e6
  N0_min = gonv_max
  D0r = 50.e-6
  !CALL cal_D0(rhoa,qr,Ntr,cr,alphar)
  mvd_r = MIN(mvd_r,2.5e-3)
  mvd_r = MAX(mvd_r,0.75*D0r)
  if (tair < 270.65 .and. qr > 1.e-12 .and. mvd_r > 100.E-6) then
    xslw1 = 4.01 + alog10(mvd_r)
  else
    xslw1 = 0.01
  endif
  rg = qg*rhoa
  ygra1 = 4.31 + alog10(max(5.E-5, rg))
  zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
  N0_exp = 10.**(zans1)
  N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
  am_g = cg
  bm_g = 3.
  oge1 = 1./(bm_g + 1.)
  cge_2 = SNGL(alphag) + 1.
  cgg_1 = WGAMMA(bm_g + 1.)
  cgg_2 = WGAMMA(cge_2)
  cgg_3 = WGAMMA(bm_g + SNGL(alphag) + 1.)
  ogg1 = 1./cgg_1
  ogg2 = 1./cgg_2
  obmg = 1./bm_g
  lam_exp = (N0_exp*am_g*cgg_1/rg)**oge1
  lamg = lam_exp*(cgg_3*ogg2*ogg1)**obmg
  ilamg = 1./lamg
    
  N0g = N0_exp/(cgg_2*lam_exp)*lamg**cge_2

END SUBROUTINE calN0g_Thompson

SUBROUTINE callamda_UM(rhoa,q,cx,na,nb,alpha,lamda)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter for graupel in the 
!            UM scheme. 
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/04/2015) 
!   
!  MODIFICATION HISTORY:
!  
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     
  IMPLICIT NONE

  REAL :: rhoa,q,cx,na,nb
  REAL*8 :: alpha, lamda,gamma4

  REAL*8 :: gamma_r8
  
  gamma4 = gamma_r8(4.+alpha)
  lamda = (gamma4*cx*na/(rhoa*q))**(1./(4.-nb+alpha))

END SUBROUTINE callamda_UM

SUBROUTINE calN0_UM(na,nb,lamda,N0)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter for graupel or rain in the 
!            UM scheme. 
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/04/2015) 
!   
!  MODIFICATION HISTORY:
!  
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     
  IMPLICIT NONE

  REAL :: na,nb
  REAL*8 :: lamda,N0
  
  REAL*8 :: gamma

  N0 = na*lamda**nb

END SUBROUTINE calN0_UM

!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
REAL FUNCTION GAMMLN(XX)
  !     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
  IMPLICIT NONE
  REAL, INTENT(IN) :: XX
  REAL*8, PARAMETER :: STP = 2.5066282746310005D0
  REAL*8, DIMENSION(6), PARAMETER:: &
           COF = (/76.18009172947146D0, -86.50532032941677D0, &
                   24.01409824083091D0, -1.231739572450155D0, &
                  .1208650973866179D-2, -.5395239384953D-5/)
  REAL*8 :: SER,TMP,X,Y
  INTEGER:: J

  X=XX
  Y=X
  TMP=X+5.5D0
  TMP=(X+0.5D0)*LOG(TMP)-TMP
  SER=1.000000000190015D0
  DO 11 J=1,6
    Y=Y+1.D0
    SER=SER+COF(J)/Y
  11    CONTINUE
  GAMMLN=TMP+LOG(STP*SER/X)
END FUNCTION GAMMLN

REAL FUNCTION WGAMMA(y)

  IMPLICIT NONE
  REAL, INTENT(IN):: y

  WGAMMA = EXP(GAMMLN(y))

END FUNCTION WGAMMA

SUBROUTINE fractionWater(qr,qi,fo,density_ice,fracqr,fracqi,fm,fw,rhom)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, fo, density_ice
  REAL,INTENT(OUT) :: fracqr, fracqi, fm, fw, rhom

  REAL :: fr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fr = 0.
  fw = 0.
  fracqr = 0.
  fracqi = 0.
  fm = 0.
  rhom = 0.

!-----------------------------------------------------------------------
! Calculate the fraction of mleting ice (fr) based on the ratio between
! qr and qi. fo is the maximum allowable fraction of melting snow.
!-----------------------------------------------------------------------
  IF (qr > 0. .AND. qi > 0.) THEN
    fr = fo*(MIN(qi/qr,qr/qi))**.3
  ENDIF

! DTD: test, set fr to 1.0
! IF (qr > 0. .AND. qi > 0.) THEN
!   fr = 1.0
! ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

!JYS   fracqr = MIN(fr*qr,qi)
!JYS   fracqi = MIN(fr*qi,qr)
   fracqr = fr * qr
   fracqi = fr * qi

  fm = fracqr + fracqi

  IF (fm .EQ. 0. .AND. qr > 0.) THEN
    fw = 1.
  ELSE IF (fm > 0.) THEN
    fw = fracqr/fm
  ENDIF

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater

SUBROUTINE fractionWater_temperature(qr,qi,density_ice,fm,fw,rhom,tair_C)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture based on the air temperature. 
! It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 7/25/2014
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, density_ice, tair_C
  REAL,INTENT(OUT) :: fm, fw, rhom
  REAL :: layer_tmax, layer_tmin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fw = 0.
  fm = 0.
  rhom = 0.

!-----------------------------------------------------------------------
! Calculate the faction of water.
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  fm = qi

! Compute the degree of wet in percentage based on air temperature
  layer_tmax = 2.5
  layer_tmin = -2.5
  if(tair_C >= layer_tmin .and. tair_C < layer_tmax) then
    fw = (tair_C - layer_tmin)/(layer_tmax-layer_tmin)
  else if(tair_C >= layer_tmax) then
    fw = 1.
  else
    fm = 0.
    fw = 0.
  endif

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater_temperature

SUBROUTINE fractionWaterD14(MPflg,ii,jj,kk,fos,rhoa,tairC,qr,qs,qg,qh,ntr,nts,ntg,nth,alfr,alfg,alfh,mug,muh,rhos,rhog,rhoh, &
                          fracqrs,fracqs,fracqrg,fracqrh,qsoakg,qsoakh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
                          rhomh,Ndrg,Ndrh,fwgbin,fwhbin,rhomgbin,rhomhbin)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 6/27/2012
!
! Based on Youngsun Jung's fractionWater subroutine, this subroutine
! uses an alternate method to compute the liquid water fraction on melting
! ice.  It also handles the case where there are multiple ice species at a given
! grid point.
!
! MODIFIED: Dan Dawson 08/03/2012
!
! Rewrote the logic of water fraction diagnostics for graupel and hail.
! The critical water mass on melting graupel and hail is now integrated over the distribution
! based on the formula of Rasmussen and Heymsfield (1987, eqn. 6)
! for the assumed discrete wet graupel/hail distribution and compared
! with the available water from the rainfield.  The minimum of the available
! water and the critical water fraction is then used for graupel and hail.  This requires
! an iterative technique using a first guess for the wet graupel/hail mass.
! Any excess, which should be slight, is added back to the rain field later in the 
! code.  For snow, we are just punting for now and allowing up to a maximum of 75%
! water fraction distributed evenly amongst sizes.  NOTE, this is the version used
! in Dawson et al. (2014)
!
! MODIFIED: Dan Dawson 04/11/2014
!
! Changed method for how much rain is transferred to graupel and hail.  Instead of
! adding almost all the rain to the graupel/hail and then iterating until convergence,
! the graupel/hail distribution is modified "in place" to determine how much water is
! needed to "saturate".  This is then compared with the available water from the rain 
! field.  If there's enough, the appropriate amount is subtracted from rain.  If not, the
! amount of water in each graupel/hail bin is reduced by the appropriate factor.
!
! In the future we may want to do this in mass space and convert back to diameter afterwards.
!
! MODIFIED: Dan Dawson Fall 2014
! 
! Further updates to fix some bugs and some conceptual errors. An iteration like the original is performed but is based on
! a first guess that starts with the amount of rain needed to saturate the initial graupel distribution.  Also reverted to original
! Jung et al. (2008) method for snow only.
!
! MODIFIED: Dan Dawson 02/04/15
!
! More minor bug fixes.  Now makes sure wet graupel and hail distributions are consistent with individual
! microphysics schemes during water fraction iteration.
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  IMPLICIT NONE
 
!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  INTEGER :: ii,jj,kk,MPflg
  REAL :: rhoa, tairC, rhos, rhog, rhoh, rhogs, rhohs, rhogd, rhohd,volg,volh,fos
  REAL*8 :: alfr,alfg,mug,alfh,muh
  REAL,INTENT(INOUT) :: qr,qs,qg,qh,ntr,nts,ntg,nth
  REAL,INTENT(OUT) :: fracqrs,fracqs,fracqrg,fracqrh,qsoakg,qsoakh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,rhomh
  REAL,INTENT(INOUT) :: fwgbin(:),fwhbin(:),rhomgbin(:),rhomhbin(:)
  REAL*8, INTENT(INOUT) :: Ndrg(:)
  REAL*8, INTENT(INOUT) :: Ndrh(:)
  
  REAL :: qgwfac,qhwfac
  REAL :: frs,frg,frh
  REAL :: qti ! Total ice mixing ratio (sum of snow, graupel, and hail)
  REAL :: qswcrit,qgwcrit,qhwcrit,qtwcrit,tmp1,tmp2,tmp3
  REAL :: qra
  REAL*8 :: db_N0,lamrg,lamrh,Ntw
  REAL :: intv,cx,cr,mtot,mice,mices,miced,micesw,mwcrit,tmpfwg,tmpfwh,diff,tmpntg,tmpnth,rhomgtmp,rhomhtmp
  REAL :: alfgtmp,alfhtmp,ratiog,ratioh,dmr,dmg,dmh,D0r,micetot,qgsave
  INTEGER :: i,iter
  LOGICAL :: adjustNt
  LOGICAL :: stopflag
!  LOGICAL :: tosoak, freezesoak,conserveice
    
!  REAL, PARAMETER :: qmin = 1.e-6 ! test higher value (erm, 3/25/2016) Threshold for minimum mixing ratio for diagnostic water routine
  REAL, PARAMETER :: qmin = 1.e-12 ! Threshold for minimum mixing ratio for diagnostic water routine
  REAL :: lwscoeff ! Liquid water shell coefficient: set to 1 above freezing, to wgfac below freezing.
  integer :: iend

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	adjustNt = .true. ! Adjust number concentration of hail and graupel during water transfers
										! to maintain mean volume diameter?
    
  IF(MPflg >= 3) adjustNt = .false. ! Don't maintain mean volume diameter for single-moment schemes
  
  frs = 0.
  frg = 0.
  frh = 0.
  fws = 0.
  fwg = 0.
  fwh = 0.
  fracqrs = 0.
  fracqrg = 0.
  fracqrh = 0.
  qra = 0.
  fms = 0.
  fmg = 0.
  fmh = 0.
  rhoms = rhos
  rhomg = rhog
  rhomh = rhoh
  qgwfac = 0.
  qhwfac = 0.
  tmpntg = ntg
  tmpnth = nth
  fwgbin = 0.
  fwhbin = 0.
  alfhtmp = alfh
  alfgtmp = alfg
  rhomgbin = 0.0
  rhomhbin = 0.0
  tmpfwg = 0.0
  tmpfwh = 0.0
  ratiog = 0.0
  ratioh = 0.0
  micesw = 0.0
  rhogs = 910.
  rhohs = 910.
  qsoakg = 0.0
  qsoakh = 0.0
  dmr = 0.0
  dmg = 0.0
  dmh = 0.0
  cr = 1000.*(pi/6.)

  ! Reduction factor only operative for T < 0 C for now
  IF(tairC > 0.) THEN
    lwscoeff = 1.0
  ELSE
    lwscoeff = wgfac
  ENDIF

  qti = qs
  IF(grpl_ON == 1) qti = qti+qg
  IF(hl_ON == 1) qti = qti+qh 
  !qti = qs+qg+qh
  
!-----------------------------------------------------------------------
! Calculate the first guess fraction of rain water to be transfered to the ice 
! surface (fr), based on the ratio of qs/qg/qh and qti up to a maximum of 100%
! of the rain.
!-----------------------------------------------------------------------

  IF (qr > qmin .AND. qti > qmin) THEN   
    frs = qs/qti
    IF(grpl_ON == 1) frg = qg/qti
    IF(hl_ON == 1) frh = qh/qti
  ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  ! For snow, use old Jung et al. (2008) method
  IF (qs > qmin .and. snow_ON == 1) CALL fractionWater(frs*qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
  IF (qg > qmin .and. frg > 0. .and. grpl_ON == 1) THEN ! Iterate to find the water fraction on graupel
    dmr = (6.*rhoa*qr/(pi*1000.*ntr))**(1./3.)
    dmg = (6.*rhoa*qg/(pi*rhog*ntg))**(1./3.)
    intv = (dsg(2) - dsg(1))*1.e-3
    diff = 1000. ! fracqrg-qgwcrit
!    qra = Min(0.99,lwscoeff)*frg*qr
    qra = 0.99*frg*qr
    qgwcrit = 0.0
    qgsave = qg
    ! First, add enough water to raise the bulk density to 910 kg m^-3
    IF (tosoak) THEN
      rhogd = rhogs
      tmp1 = qg*rhoa/rhog ! Total volume of graupel/unit volume air
      qsoakg = tmp1*rhogs/rhoa ! Mass with same total volume as above but with density = 910 kg/m^3
      qsoakg = qsoakg-qg ! Amount of water needed to raise density to 910 kg/m^3
!        IF ( ii==15 ) THEN ! .and. jj==90 .and. kk==17 ) THEN
!          write(6,*) 'qsoakg,qra,qg',qsoakg,qra,qg
!        ENDIF
      IF (qsoakg >= qra) THEN ! Not enough water to spill onto shell, skip shell iteration
        qsoakg = qra
        qg = qg+qsoakg
        fmg = qg
        tmpntg = ntg
        rhomg = rhog*qg/(qg-qsoakg)
        rhogd = rhomg
        rhomgbin = rhomg
        IF(freezesoak) THEN
          fwg = 0.0
        ELSE
          fwg = qsoakg/fmg
        ENDIF
        fwgbin = fwg

!        IF ( ii==15 ) THEN ! .and. jj==90 .and. kk>=6 .and. kk<=18 ) THEN
!          write(6,*) 'kk,qsoakg,qra,rhomg,rhog,ntg,fmg',kk,qsoakg,qra,rhomg,rhog,ntg,fmg,fmg-qsoakg
!        ENDIF

        ! Set up new graupel distribution
        if(tmpntg > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment graupel
          CALL cal_N0(rhoa,fmg,tmpntg,rhomg,alfg,db_N0,mug)
          CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
        else ! Single-moment
          cx = rhomg*(pi/6.)
          IF(MPflg == 3) THEN ! Diagnostic N0 for UM scheme
            CALL callamda_UM(rhoa,fmg,cx,nag,nbg,alfg,lamrg)
            CALL calN0_UM(nag,nbg,lamrg,db_N0)
          ELSEIF(MPflg == 4) THEN ! Diagnostic N0 for Thompson scheme
            CALL cal_D0(rhoa,qr,ntr,1000.,alfr,D0r) ! Compute median volume diameter for rain (using original rain distribution)
            CALL calN0g_Thompson(tairC+273.15,rhoa,fmg,cx,alfg,rhomg,qr,D0r,db_N0)
          ELSE
            db_N0 = dble(N0g)
          ENDIF
          CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,tmpntg,mug)
          IF(MPflg /= 3) CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
        endif
        Ntw = tmpntg
        DO i=1,nd_g
          Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
          Ntw = Ntw - Ndrg(i)
          if(Ntw <= 0.d0) Ndrg(i) = 0.d0
        ENDDO
        stopflag = .true.
      ELSEIF(qsoakg > 0.) THEN ! First add the water for soaking
        qg = qg+qsoakg
        fmg = qg
        tmpntg = ntg
        rhomg = 910.
        rhogd = rhomg
        rhomgbin = rhomg
        qra = qra-qsoakg ! available rainwater minus that used for soaking
        fracqrg = qra
        stopflag = .false.
!        IF ( ii==15 ) THEN ! 121 .and. jj==90 .and. kk>=6 .and. kk<=18 ) THEN
!          write(6,*) 'kk,qsoakg,qr,qgsave,rhomg,rhog,',kk,qsoakg,qr,qgsave,rhomg,rhog
!        ENDIF
      ELSE ! Don't soak
        qsoakg = 0.0
        rhomg = rhog
        rhogd = rhog
        rhomgbin = rhomg
        fmg = qg
        tmpntg = ntg
        fracqrg = qra ! Available rain water
        stopflag = .false.
      ENDIF
    ELSE ! Don't soak
      qsoakg = 0.0
      rhomg = rhog
      rhogd = rhog
      rhomgbin = rhog
      fmg = qg
      tmpntg = ntg
      fracqrg = qra ! Available rain water
      stopflag = .false.
    ENDIF
!    IF ( ii==15 ) THEN
!     write(6,*) 'i,j,k,qsoakg,qra,qg,ntg,rhogd,rhog = ',ii,jj,kk,qsoakg,qra,qg,ntg,rhogd,rhog
!     write(6,*) 'fractionrain = ',qra/(qra+qg)
!    ENDIF 
    ! Iterate to find how much water graupel distribution can "hold" on the shell.
    iter = 0
    tmp1 = qsoakg
    DO WHILE(diff > 1.e-3 .and. .not. stopflag) ! Iterate until convergence
      qgwcrit = 0.0
      ! Compute N0 and lambda for graupel distribution
      if(tmpntg > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment graupel
        CALL cal_N0(rhoa,fmg,tmpntg,rhomg,alfg,db_N0,mug)
        CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
      else ! Single-moment
        cx = rhomg*(pi/6.)
        IF(MPflg == 3) THEN ! Diagnostic N0 for UM scheme
          CALL callamda_UM(rhoa,fmg,cx,nag,nbg,alfg,lamrg)
          CALL calN0_UM(nag,nbg,lamrg,db_N0)
        ELSEIF(MPflg == 4) THEN ! Diagnostic N0 for Thompson scheme
          CALL cal_D0(rhoa,qr,ntr,1000.,alfr,D0r) ! Compute median volume diameter for rain (using original rain distribution)
          CALL calN0g_Thompson(tairC+273.15,rhoa,fmg,cx,alfg,rhomg,qr,D0r,db_N0)
        ELSE
          db_N0 = dble(N0g)
        ENDIF
        CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,tmpntg,mug)
        IF(MPflg /= 3) CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
      endif
      Ntw = tmpntg
      rhomgtmp = 0.0
      tmp3 = 0.0
      iend = 1
      micetot = 0.0
      DO i = 1,nd_g
        volg = (1./6.)*pi*(dsg(i)*1.0e-3)**3.
        mice = volg*rhog   ! Mass of original dry ice sphere
        miced = volg*rhogd ! Mass of "dry" ice sphere after soaking (if any, otherwise will be equal to mice)
        mices = volg*rhogs ! Mass of "dry" ice sphere as if it were soaked 
        micesw = miced-mice ! Mass of soaked water in ice sphere (if any, otherwise will be 0) 
      
        IF(dsg(i) <= 8.0) THEN ! Try to completely melt graupel less than or equal to 8 mm (EDIT: Now do so up to factor set by lwscoeff)
          mwcrit = lwscoeff*1000.*volg
          miced = rhogd*(volg-mwcrit/1000.) ! is 0 if lwscoeff=1
          mtot = mwcrit+miced
        ELSE ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula (multiplied by lwscoeff factor)
          ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
          mtot = volg*rhomgbin(i)
          mwcrit = lwscoeff*(0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
          mwcrit = MIN(mwcrit,mtot)
          miced = mtot-mwcrit ! New ice mass
        ENDIF
        mice = miced*(rhog/rhogd) ! New ice mass assuming original density
        micesw = MAX(0.0,miced-mice) ! New soaked water mass
        tmp2 = miced/rhogd ! Volume of ice assuming 910 density
        tmp2 = MAX(0.0,volg-tmp2) ! New volume of water
        mwcrit = tmp2*1000. ! Adjusted water mass
        
        rhomgbin(i) = MAX(0.0,(miced+mwcrit)/volg) ! New bulk density
        mtot = mwcrit+miced ! New total mass 
        IF(tosoak .and. .not. freezesoak) THEN
          tmpfwg = (mwcrit+micesw)/mtot
        ELSE
          tmpfwg = mwcrit/mtot
        ENDIF
        
        Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
        Ntw = Ntw - Ndrg(i)
        if(Ntw <= 0.d0) Ndrg(i) = 0.d0 
        fwgbin(i) = tmpfwg
        rhomgtmp = rhomgtmp + Ndrg(i)*rhomgbin(i)
        tmp3 = tmp3 + Ndrg(i)
        micetot = micetot + mice*Ndrg(i)/rhoa 
        qgwcrit = qgwcrit + mwcrit*Ndrg(i)/rhoa ! tmpfwg*Ndrg(i)*volg*rhomgbin(i)/rhoa ! 2nd term is mass mixing ratio of liquid water 
                                                           ! in current diameter bin    
!        IF ( ii==15 ) THEN ! 121 .and. jj==90 .and. kk==17 ) THEN
!          write(6,*) 'i,dsg,N,tmp3,qgwcrit,qra',i,dsg(i),Ndrg(i),tmp3,qgwcrit,qra,micetot,miced,rhog,rhogd
!        ENDIF

        iend = i
        IF(qgwcrit > qra) THEN
          stopflag = .true.
          qgwcrit = qra
          EXIT
        ENDIF
        
        IF(Ntw <= 0.0) EXIT
            
      ENDDO
      ! Compute new rhomg and update qg/ntg
      rhomg = rhomgtmp/tmp3
!        IF ( ii==15 ) THEN ! 121 .and. jj==90 .and. kk==17 ) THEN
!          write(6,*) 'iend,rhomg',iend,rhomg,rhomgtmp,tmp3,micetot
!        ENDIF
      IF(qgwcrit > qra) THEN
        stopflag = .true.
        qgwcrit = qra
      ENDIF
      IF ( qgwcrit > 0.0 ) THEN
       diff = ABS(fracqrg-qgwcrit)/qgwcrit
      qgwfac = MIN(fracqrg/qgwcrit,1.0)
      fracqrg = qgwcrit
      ELSE
        write(0,*) 'Problem with qgwcrit=0, stopping'
        stop
        diff = 0.0
        qgwfac = 1.0
        ! what to do here?
      ENDIF

      IF(conserveice) fmg = qg + fracqrg ! Note, any soaked water has already been added to qg!
      IF(adjustNt .and. conserveice) THEN
        tmpntg = ntg+ntg*fracqrg/qg
      ELSEIF(MPflg < 3 .or. .not. conserveice) THEN
        tmpntg = ntg
      ENDIF
      IF(fmg > 0.) THEN
        IF(tosoak .and. .not. freezesoak) THEN
          fwg = (fracqrg+tmp1)/fmg
        ELSE
          fwg = fracqrg/fmg
        ENDIF
      ELSE
        fwg = 0.0
      END IF
      iter = iter+1

!        IF ( ii==15 ) THEN ! 121 .and. jj==90 .and. kk==17 ) THEN
!          write(6,*) 'iter,ntw,fwg',iter,ntw,fwg
!        ENDIF
    
      IF(iter > 1000) stopflag = .true. ! Catchall condition to avoid infinite loop from lack of convergence
    ENDDO
    qsoakg = tmp1

    ! Recompute graupel distribution
    if(tmpntg > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment graupel
      CALL cal_N0(rhoa,fmg,tmpntg,rhomg,alfg,db_N0,mug)
      CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
    else ! Single-moment
      cx = rhomg*(pi/6.)
      IF(MPflg == 3) THEN ! Diagnostic N0 for UM scheme
        CALL callamda_UM(rhoa,fmg,cx,nag,nbg,alfg,lamrg)
        CALL calN0_UM(nag,nbg,lamrg,db_N0)
      ELSEIF(MPflg == 4) THEN ! Diagnostic N0 for Thompson scheme
        CALL cal_D0(rhoa,qr,ntr,1000.,alfr,D0r) ! Compute median volume diameter for rain (using original rain distribution)
        CALL calN0g_Thompson(tairC+273.15,rhoa,fmg,cx,alfg,rhomg,qr,D0r,db_N0)
      ELSE
        db_N0 = dble(N0g)
      ENDIF
      CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,tmpntg,mug)
      IF(MPflg /= 3) CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
    endif

    ! Loop through the bins one more time and set fwgbin to 0 for those beyond which we have water.
    Ntw = tmpntg
    tmp2 = 0.0
    rhomgtmp = 0.0
    IF(tosoak .and. .not. freezesoak) THEN
      tmp1 = fracqrg+qsoakg
    ELSE
      tmp1 = fracqrg
    ENDIF
    DO i = 1,nd_g
      Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
      Ntw = Ntw - Ndrg(i)
      if(Ntw <= 0.d0) Ndrg(i) = 0.d0 
      tmp2 = tmp2 + Ndrg(i)
      if(tmp1 <= 0.) THEN
        fwgbin(i) = 0.
        rhomgbin(i) = rhogd ! Original (possibly after soaking) density
      ENDIF
      tmp1 = tmp1 - fwgbin(i)*Ndrg(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhomgbin(i)/rhoa
      rhomgtmp = rhomgtmp + Ndrg(i)*rhomgbin(i)
    ENDDO
    rhomg = rhomgtmp/tmp2 ! New average bulk density for distribution
  ENDIF

!    IF ( ii==15 ) write(6,*) 'i,j,k,rhog,rhomg,rhomgtmp,tmp2 = ',ii,jj,kk,rhog,rhomg,rhomgtmp,tmp2 
  
  IF(qh > qmin .and. frh > 0. .and. hl_ON == 1) THEN ! Iterate to find the water fraction on hail
    dmr = (6.*rhoa*qr/(pi*1000.*ntr))**(1./3.)
    dmh = (6.*rhoa*qh/(pi*rhoh*nth))**(1./3.)
    intv = (dsh(2) - dsh(1))*1.e-3
    diff = 1000. ! fracqrh-qhwcrit
!    qra = Min(0.99,lwscoeff)*frh*qr
    qra = 0.99*frh*qr
    qhwcrit = 0.0
    ! First, add enough water to raise the bulk density to 910 kg m^-3
    IF(tosoak) THEN
      rhohd = rhohs
      tmp1 = qh*rhoa/rhoh ! Total volume of hail/unit volume air
      qsoakh = tmp1*rhohs/rhoa ! Mass with same total volume as above but with density = 910 kg/m^3
      qsoakh = qsoakh-qh ! Amount of water needed to raise density to 910 kg/m^3
      IF(qsoakh >= qra) THEN ! Not enough water to spill onto shell, skip shell iteration
        qsoakh = qra
        qh = qh+qsoakh
        fmh = qh
        tmpnth = nth
        rhomh = rhoh*qh/(qh-qsoakh)
        rhohd = rhomh
        rhomhbin = rhomh
        IF(freezesoak) THEN
          fwh = 0.0
        ELSE
          fwh = qsoakh/fmh
        ENDIF
        fwhbin = fwh
        ! Set up new hail distribution
        if(tmpnth > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment hail
          CALL cal_N0(rhoa,fmh,tmpnth,rhomh,alfh,db_N0,muh)
          CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
        else ! Single-moment
          cx = rhomh*(pi/6.)
          db_N0 = dble(N0h)
          CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,tmpnth,muh)
          CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
        endif
        Ntw = tmpnth
        DO i=1,nd_h
          Ndrh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamrh*(dsh(i)*1.e-3)**(3.0*muh))*intv
          Ntw = Ntw - Ndrh(i)
          if(Ntw <= 0.d0) Ndrh(i) = 0.d0
        ENDDO
        stopflag = .true.
      ELSEIF(qsoakh > 0.) THEN ! First add the water for soaking
        qh = qh+qsoakh
        fmh = qh
        tmpnth = nth
        rhomh = 910.
        rhohd = rhomh
        rhomhbin = rhomh
        qra = qra-qsoakh ! available rainwater minus that used for soaking
        fracqrh = qra
        stopflag = .false.
      ELSE ! Don't soak
        qsoakh = 0.0
        rhomh = rhoh
        rhohd = rhoh
        rhomhbin = rhomh
        fmh = qh
        tmpnth = nth
        fracqrh = qra ! Available rain water
        stopflag = .false.
      ENDIF
    ELSE ! Don't soak
      qsoakh = 0.0
      rhomh = rhoh
      rhohd = rhoh
      rhomhbin = rhoh
      fmh = qh
      tmpnth = nth
      fracqrh = qra ! Available rain water
      stopflag = .false.
    ENDIF
    ! Iterate to find how much water hail distribution can "hold" on the shell.
    iter = 0
    tmp1 = qsoakh
    DO WHILE(diff > 1.e-3 .and. .not. stopflag) ! Iterate until convergence
      qhwcrit = 0.0
      ! Compute N0 and lambda for hail distribution
      if(tmpnth > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment hail
        CALL cal_N0(rhoa,fmh,tmpnth,rhomh,alfh,db_N0,muh)
        CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
      else ! Single-moment
        cx = rhomh*(pi/6.)
        db_N0 = dble(N0h)
        CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,tmpnth,muh)
        CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
      endif
      Ntw = tmpnth
      rhomhtmp = 0.0
      tmp3 = 0.0
      DO i = 1,nd_h
        volh = (1./6.)*pi*(dsh(i)*1.0e-3)**3.
        mice = volh*rhoh   ! Mass of original dry ice sphere
        miced = volh*rhohd ! Mass of "dry" ice sphere after soaking (if any, otherwise will be equal to mice)
        mices = volh*rhohs ! Mass of "dry" ice sphere as if it were soaked 
        micesw = miced-mice ! Mass of soaked water in ice sphere (if any, otherwise will be 0) 
      
        IF(dsh(i) <= 8.0) THEN ! Try to completely melt hail less than or equal to 8 mm (EDIT: Now do so up to factor set by lwscoeff)
          mwcrit = lwscoeff*1000.*volh
          miced = rhohd*(volh-mwcrit/1000.)
          mtot = mwcrit+miced
        ELSE ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula (multiplied by lwscoeff factor)
          ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
          mtot = volh*rhomhbin(i)
          mwcrit = lwscoeff*(0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
          mwcrit = MIN(mwcrit,mtot)
          miced = mtot-mwcrit ! New ice mass
        ENDIF
        mice = miced*(rhoh/rhohd) ! New ice mass assuming original density
        micesw = MAX(0.0,miced-mice) ! New soaked water mass
        tmp2 = miced/rhohd ! Volume of ice assuming 910 density
        tmp2 = MAX(0.0,volh-tmp2) ! New volume of water
        mwcrit = tmp2*1000. ! Adjusted water mass
        
        rhomhbin(i) = MAX(0.0,(miced+mwcrit)/volh) ! New bulk density
        mtot = mwcrit+miced ! New total mass 
        IF(tosoak .and. .not. freezesoak) THEN
          tmpfwh = (mwcrit+micesw)/mtot
        ELSE
          tmpfwh = mwcrit/mtot
        ENDIF
        
        Ndrh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamrh*(dsh(i)*1.e-3)**(3.0*muh))*intv
        Ntw = Ntw - Ndrh(i)
        if(Ntw <= 0.d0) Ndrh(i) = 0.d0 
        fwhbin(i) = tmpfwh
        rhomhtmp = rhomhtmp + Ndrh(i)*rhomhbin(i)
        tmp3 = tmp3 + Ndrh(i)
        qhwcrit = qhwcrit + mwcrit*Ndrh(i)/rhoa ! tmpfwh*Ndrh(i)*volh*rhomhbin(i)/rhoa ! 2nd term is mass mixing ratio of liquid water 
                                                           ! in current diameter bin    
        IF(qhwcrit > qra) THEN
          stopflag = .true.
          qhwcrit = qra
          EXIT
        ENDIF
        IF(Ntw <= 0.0) EXIT
            
      ENDDO
      ! Compute new rhomh and update qh/nth
      rhomh = rhomhtmp/tmp3
      IF(qhwcrit > qra) THEN
        stopflag = .true.
        qhwcrit = qra
      ENDIF
      diff = ABS(fracqrh-qhwcrit)/qhwcrit
      qhwfac = MIN(fracqrh/qhwcrit,1.0)
      fracqrh = qhwcrit
      IF(conserveice) fmh = qh + fracqrh ! Note, any soaked water has already been added to qh!
      IF(adjustNt .and. conserveice) THEN
        tmpnth = nth+nth*fracqrh/qh
      ELSEIF(MPflg < 3 .or. .not. conserveice) THEN
        tmpnth = nth
      ENDIF
      IF(fmh > 0.) THEN
        IF(tosoak .and. .not. freezesoak) THEN
          fwh = (fracqrh+tmp1)/fmh
        ELSE
          fwh = fracqrh/fmh
        ENDIF
      ELSE
        fwh = 0.0
      END IF
      iter = iter+1
      IF(iter > 1000) stopflag = .true. ! Catchall condition to avoid infinite loop from lack of convergence
    ENDDO
    qsoakh = tmp1
    
    ! Recompute hail distribution
    if(tmpnth > 0.0 .and. MPflg < 3) then ! Multi-moment case, MPflg >= 3 have single-moment hail
      CALL cal_N0(rhoa,fmh,tmpnth,rhomh,alfh,db_N0,muh)
      CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
    else ! Single-moment
      cx = rhomh*(pi/6.)
      db_N0 = dble(N0h)
      CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,tmpnth,muh)
      IF(MPflg /= 3) CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
    endif

    ! Loop through the bins one more time and set fwhbin to 0 for those beyond which we have water.
    Ntw = tmpnth
    tmp2 = 0.0
    rhomhtmp = 0.0
    IF(tosoak .and. .not. freezesoak) THEN
      tmp1 = fracqrh+qsoakh
    ELSE
      tmp1 = fracqrh
    ENDIF
    DO i = 1,nd_h
      Ndrh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamrh*(dsh(i)*1.e-3)**(3.0*muh))*intv
      Ntw = Ntw - Ndrh(i)
      if(Ntw <= 0.d0) Ndrh(i) = 0.d0 
      tmp2 = tmp2 + Ndrh(i)
      if(tmp1 <= 0.) THEN
        fwhbin(i) = 0.
        rhomhbin(i) = rhohd ! Original (possibly after soaking) density
      ENDIF
      tmp1 = tmp1 - fwhbin(i)*Ndrh(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhomhbin(i)/rhoa
      rhomhtmp = rhomhtmp + Ndrh(i)*rhomhbin(i)
    ENDDO
    rhomh = rhomhtmp/tmp2 ! New average bulk density for distribution
  ENDIF
  
  IF(grpl_ON == 1) THEN
  	IF(conserveice) THEN 
  	  fmg = qg+fracqrg
  	ELSE
  	  fmg = qg
  	ENDIF
  	ntg = tmpntg
  	IF(fmg > 0.) THEN
  	  IF(tosoak .and. .not. freezesoak) THEN
        fwg = (fracqrg+qsoakg)/fmg
      ELSE
        fwg = fracqrg/fmg
      ENDIF
    ELSE
      fwg = 0.0
    END IF
  ENDIF
  
  IF(hl_ON == 1) THEN
  	IF(conserveice) THEN
  	  fmh = qh+fracqrh
  	ELSE
  	  fmh = qh
  	ENDIF
  	nth = tmpnth
  	IF(fmh > 0.) THEN
  	  IF(tosoak .and. .not. freezesoak) THEN
        fwh = (fracqrh+qsoakh)/fmh
      ELSE
        fwh = fracqrh/fmh
      ENDIF
    ELSE
      fwh = 0.0
    END IF
  ENDIF

END SUBROUTINE fractionWaterD14

SUBROUTINE adjustWater(ii,jj,kk,rhoa,qr,qs,qg,qh,ntr,nts,ntg,nth,alfg,alfh,mug,muh,rhos,rhog,rhoh, &
                          fracqrs,fracqrg,fracqrh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
                          rhomh,fwgbin,fwhbin)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine is similar to fractionWater3 but instead adjusts the predicted
! fraction of water on graupel and hail (from the ZVDM scheme) to be consistent with the experimental
! findings of Rasmussen and Heymsfield.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 02/07/2013
!
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  IMPLICIT NONE
 
!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  INTEGER :: ii,jj,kk
  INTEGER, PARAMETER :: nd_gtmp = 625
  INTEGER, PARAMETER :: nd_htmp = 875
  REAL :: rhoa, qr, qs, qg, qh, rhos, rhog, rhoh
  REAL*8 :: alfg,mug,alfh,muh
  REAL,INTENT(IN) :: rhoms,rhomg,rhomh
  REAL,INTENT(INOUT) :: ntr,nts,ntg,nth,fracqrs,fracqrg,fracqrh
  REAL,INTENT(INOUT) :: fms,fmg,fmh,fws,fwg,fwh
  real :: fmg0,fmh0,fwg0,fwh0
  REAL:: fwgbin(:),fwhbin(:)
  
  REAL :: qgwfac,qhwfac
  REAL :: frs,frg,frh
  REAL :: qti ! Total ice mixing ratio (sum of snow, graupel, and hail)
  REAL :: qswcrit,qgwcrit,qhwcrit,qtwcrit,tmp1
  real :: qhtot,qgtot
  REAL*8 :: db_N0,lamrg,lamrh,Ntw
  REAL :: intv,cx,mtot,mice,mwcrit,tmpfwg,diff,tmpntg,tmpnth
  INTEGER :: i,iter
!  INTEGER, PARAMETER :: nd_r = 112 ! 100
!  INTEGER, PARAMETER :: nd_s = 112 
!  INTEGER, PARAMETER :: nd_g = 112 ! 625
!  INTEGER, PARAMETER :: nd_h = 112 ! 875
  !REAL*8, DIMENSION (nd) :: Ndr, Nds, Ndh, Ndg, Ndrs, Ndrh, Ndrg
  REAL*8, DIMENSION (nd_r) :: Ndr
  REAL*8, DIMENSION (nd_s) :: Nds, Ndrs
  REAL*8, DIMENSION (nd_g) :: Ndg, Ndrg
  REAL*8, DIMENSION (nd_h) :: Ndh, Ndrh
  
  real, parameter :: d1t = 1000.*(6.0 * 0.268e-3/(917.* pi))**(1./3.) ! d1t is the diameter of the ice sphere with the mass (0.268e-3 kg) of an 8mm spherical drop
  real :: f
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	fwgbin = 0.0
	fwhbin = 0.0
  frs = 0.
  frg = 0.
  frh = 0.
!  fws = 0.
!  fwg = 0.
!  fwh = 0.
!  fracqrs = 0.
!  fracqrg = 0.
!  fracqrh = 0.
!  fms = 0.
!  fmg = 0.
!  fmh = 0.
!  rhoms = 0.
!  rhomg = 0.
!  rhomh = 0.
  qgwfac = 0.
  qhwfac = 0.
  tmpntg = ntg
  tmpnth = nth
  
  
  fmg0 = fmg
  fmh0 = fmh
  fwg0 = fwg
  fwh0 = fwh
  
!  write(0,*) 'start of adjustwater: i,j,k = ',ii,jj,kk

  IF(fmg > 0.) THEN
    fwg = Min(1.0,fracqrg/fmg)
  ELSE
  	fwg = 0.0
  END IF
  
  ! turn off graupel temporarily and work on hail first
  
    fwgbin(:) = fwg
  
  IF ( .true. ) THEN
  
! Special treatment for 100% water graupel
  IF(qg > 0.0 .and. fwg == 1.0) THEN

    fwgbin(:) = 1.0

  ELSEIF(qg > 0.0 .and. grpl_ON == 1) THEN ! Iterate to find the water fraction on graupel
  
    ! First recompute qg as the mass of dry graupel.  This will be shortly added back during the iteration, which precedes identically to fractionWater3
!    qg = qg - fracqrg ! fracqrg is the *original predicted* water fraction here
!    ntg = ntg - ntg*fracqrg/fmg ! Adjust ntg to maintain mean diameter
    qgwcrit = 0.0
    diff = fracqrg-qgwcrit
    iter=0
!    DO WHILE(diff > 1.0e-9)
      qgwcrit = 0.0
      qgtot = 0.0
      db_N0 = 0.d0
      intv = (dsg(2) - dsg(1))*1.e-3
      Ntw = 0.d0
  !    fmg = qg+fracqrg
  !    tmpntg = ntg+ntg*fracqrg/qg
      !tmpntg = ntg
      IF(fmg > 0.) THEN
        fwg = fracqrg/fmg
      ELSE
        fwg = 0.0
      END IF
      if(ntg > 0.0) then
        CALL cal_N0(rhoa,fmg,tmpntg,rhomg,alfg,db_N0,mug)
        CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
        Ntw = tmpntg
      else
        db_N0 = dble(N0g)
        cx = rhomg*(pi/6.)
        CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,tmpntg,mug)
        CALL cal_lamda(rhoa,fmg,tmpntg,rhomg,alfg,lamrg,mug)
        Ntw = tmpntg
      endif
      Do i = 1,nd_g
        if(dsg(i) <= d1t) THEN ! Try to completely melt graupel less than or equal to 8 mm
          !tmpfwg = fwh
          tmpfwg = 1.0	
        else ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula
          ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
          mtot = (1./6.)*pi*(dsg(i)*1.0e-3)**3.*MAX(917.,rhomg) ! Mass of wet hailstone of diameter dsg(i): treat it as high density even if it isn't
          !mice = (1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhog ! rhog
          !mwcrit = 0.268e-3 + 0.1389*mice*rhog/910.
          !mtot = mwcrit+mice
          mwcrit = (0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
          tmpfwg = mwcrit/mtot  ! Max water fraction for wet hailstone of diameter dsg(i)
        endif
        Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
        Ntw = Ntw - Ndrg(i)
        if(Ntw <= 0.d0) Ndrg(i) = 0.d0 
        qgwcrit = qgwcrit + tmpfwg*Ndrg(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhomg/rhoa ! rhomg ! 2nd term is mass mixing ratio of liquid water in 
                                                                  ! current diameter bin
        qgtot = qgtot + 1.0*Ndrg(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhomg/rhoa ! rhomg ! 2nd term is mass mixing ratio of liquid water in 
        fwgbin(i) = tmpfwg

      ENDDO 

      diff = fracqrg-qgwcrit

        
        ! IF fracqrh > qhwcrit, then the model is allowing more liquid than calculated here
        ! and we will just leave that be for now. (And eventually figure out why the model has higher fraction)
        ! For the opposite case, just scale fwhbin:
        IF ( fracqrg < qgwcrit ) THEN
           f = fracqrg/qgwcrit
           DO i = 1,nd_g
             fwgbin(i) = f*fwgbin(i)
           ENDDO
        ENDIF

!     IF(ii == 28 .and. jj == 54 .and. kk == 1) THEN
!       write(*,*) 'iter,fracqrg,qgwcrit,fwg',iter,fracqrg,qgwcrit,fwg
!     ENDIF
!      IF(diff <= 0.0 .and. iter == 0) THEN ! water fraction to be added less than critical, so exit loop
!          EXIT
!      ELSEIF ( iter > 20 ) THEN
!!        write(0,*) 'qgw: too many iterations. iter = ',iter
!        EXIT
!      ELSE
!        fracqrg = qgwcrit
!        iter=iter+1
!      ENDIF
!    ENDDO
    ! Loop through the bins one more time and set fwgbin to 0 for those beyond which we have water.
!    tmp1 = fracqrg
!    DO i = 1,nd_g
!      if(tmp1 <= 0.) fwgbin(i) = 0.
!      tmp1 = tmp1 - fwgbin(i)*Ndrg(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhomg/rhoa
!    ENDDO
!    
!    qgwfac = MIN(fracqrg/qgwcrit,1.0)
!    
!    fmg = qg+fracqrg
!    ntg = tmpntg
!    qg = fmg ! Now qg is just fmg again (full mass of dry+wet graupel)
    
  ENDIF
  
  ENDIF ! .false./.true.
    
  IF(fmh > 0.) THEN
    fwh = fracqrh/fmh
  ELSE
  	fwh = 0.0
  END IF
  
  ! Special treatment for 100% water hail
  
  IF(qh > 0.0 .and. fwh == 1.0) THEN
  
    fwhbin(:) = 1.0
  
  ELSEIF(qh > 0. .and. hl_ON == 1) THEN ! Iterate to find the water fraction on hail
    ! First recompute qh as the mass of dry hail.  This will be shortly added back during the iteration, which precedes identically to fractionWater3
   ! qh = qh - fracqrh ! fracqrh is the *original predicted* water fraction here
   ! nth = nth - nth*fracqrh/fmh ! Adjust nth to maintain mean diameter
    qhwcrit = 0.0
    diff = fracqrh-qhwcrit
    iter=0

!    DO WHILE(diff > 1.0e-9)

      qhwcrit = 0.0
      qhtot = 0.0
      db_N0 = 0.d0
      intv = (dsh(2) - dsh(1))*1.e-3
      Ntw = 0.d0
!      fmh = qh+fracqrh
!      tmpnth = nth+nth*fracqrh/qh
!     fmh is qh (total: ice + liquid)
!     fwh is hail water fraction
!     fracqrh is qhw
      tmpnth = nth
      IF(fmh > 0.) THEN
        fwh = fracqrh/fmh
      ELSE
        fwh = 0.0
      END IF
      if(tmpnth > 0.0) then
        CALL cal_N0(rhoa,fmh,tmpnth,rhomh,alfh,db_N0,muh)
        CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
        Ntw = tmpnth
!        IF(kk == 1 .and. ii == 45 .and. jj == 62) THEN
!          write(*,*) 'db_N0,lamrh,Ntw',db_N0,lamrh,Ntw
!        ENDIF
      else
        db_N0 = dble(N0h)
        cx = rhomh*(pi/6.)
        CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,tmpnth,muh)
        CALL cal_lamda(rhoa,fmh,tmpnth,rhomh,alfh,lamrh,muh)
        Ntw = tmpnth
      endif
      Do i = 1,nd_h
        if(dsh(i) <= d1t) THEN ! Try to completely melt hail less than or equal to 9 mm
          tmpfwg = 1.0
        else ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula
          ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
          mtot = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*MAX(917.,rhomh) ! Mass of wet hailstone of diameter dsh(i): treat it as high density even if it isn't
          !mice = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*910. ! rhoh
!          mwcrit = 0.268 + 0.1389*mice
!          mtot = mwcrit+mice
          !mice = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhoh
          mwcrit = (0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
          tmpfwg = mwcrit/mtot  ! Max water fraction for wet hailstone of diameter dsg(i)
          IF ( tmpfwg > 1.0 ) THEN
!            write(0,*) 'bin,mtot,mwcrit,tmpfwg = ',i,mtot,mwcrit,tmpfwg 
          ENDIF
        endif
        Ndrh(i) = db_N0*(dble(dsh(i)*1.e-3))**alfh*exp(-lamrh*(dble(dsh(i)*1.e-3))**(3.d0*muh))*dble(intv)
        Ntw = Ntw - Ndrh(i)
        if(Ntw <= 0.d0) Ndrh(i) = 0.d0 
        qhwcrit = qhwcrit + tmpfwg*sngl(Ndrh(i))*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhomh/rhoa ! rhomh ! 2nd term is mass mixing ratio of liquid water in 
        qhtot = qhtot + 1.0*sngl(Ndrh(i))*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhomh/rhoa ! rhomh ! 2nd term is mass mixing ratio of liquid water in 
                                                                  ! current diameter bin
!        IF(kk == 1 .and. ii == 45 .and. jj == 62) THEN                 
!          	write(*,*) 'i,dsh(i),Ndrh(i),tmpfwg,qhwcrit',i,dsh(i),Ndrh(i),tmpfwg,qhwcrit
!        ENDIF
         fwhbin(i) = tmpfwg
      ENDDO 
      diff = fracqrh-qhwcrit

!        write(0,*) 'i,j,fracqrh,qhwcrit',ii,jj,fracqrh,qhwcrit
!        write(0,*) 'qhtot,fmh = ',qhtot,fmh
!        write(0,*) 'fwh0,fwh,fmax =', fwh0,fwh, qhwcrit/qhtot
!        write(0,*) 'fmh,tmpnth,rhomh',fmh,tmpnth,rhomh
        
        
        
        ! IF fracqrh > qhwcrit, then the model is allowing more liquid than calculated here
        ! and we will just leave that be for now. (And eventually figure out why the model has higher fraction)
        ! For the opposite case, just scale fwhbin:
        IF ( fracqrh < qhwcrit ) THEN
           f = fracqrh/qhwcrit
           DO i = 1,nd_h
             fwhbin(i) = f*fwhbin(i)
           ENDDO
        ENDIF

!      ENDIF

!      IF(diff <= 0.0 .and. iter == 0) THEN ! water fraction to be added less than critical, so exit loop
!        EXIT
!      ELSEIF ( iter > 20 ) THEN
!!        write(0,*) 'qhlw: too many iterations. iter = ',iter
!        EXIT
!      ELSE
!        fracqrh = qhwcrit
!        iter=iter+1
!      ENDIF
!    ENDDO ! do while


! -------
    IF ( .false. ) THEN
    ! Loop through the bins one more time and set fwgbin to 0 for those beyond which we have water.
    tmp1 = fracqrh
    DO i = 1,nd_h
      if(tmp1 <= 0.) fwhbin(i) = 0.
      tmp1 = tmp1 - fwhbin(i)*Ndrh(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhomh/rhoa
    ENDDO
    
    qhwfac = MIN(fracqrh/qhwcrit,1.0)
    
    fmh = qh+fracqrh
    nth = tmpnth
    qh = fmh ! Now qh is just fmh again (full mass of dry+wet hail)
    
    ENDIF
! ----------

  ENDIF

  IF(grpl_ON == 1) THEN
!  	fmg = qg+fracqrg
!  	ntg = tmpntg
!	qg = fmg ! Now qg is just fmg again (full mass of dry+wet graupel)
  ENDIF
  IF(hl_ON == 1) THEN
!  	fmh = qh+fracqrh
!  	nth = tmpnth
!  	qh = fmh ! Now qh is just fmh again (full mass of dry+wet hail)
  ENDIF

  ! Check water fraction of graupel, hail, and snow. If it is nearly 100% add everything
  ! to the rain field, preserving the mean size.
  
!  IF(fws >= 0.90) THEN
!    fws = 0.0
!    fms = 0.0
!    nts = 0.0
!    fracqrs = 0.0
!    qr = qr+qs
!    ntr = ntr+qs*ntr/qr
!  ENDIF
!  
!  IF(fwg >= 0.90) THEN
!    fwg = 0.0
!    fmg = 0.0
!    ntg = 0.0
!    fracqrg = 0.0
!    qr = qr+qg
!    ntr = ntr+qg*ntr/qr
!  ENDIF
!  
!  IF(fwh >= 0.999) THEN
!    fwh = 0.0
!    fmh = 0.0
!    nth = 0.0
!    fracqrh = 0.0
!    qr = qr+qh
!    ntr = ntr+qh*ntr/qr
!  ENDIF 
   
!   write(0,*) 'leaving adjustwater'

END SUBROUTINE adjustWater

SUBROUTINE fractionWaterTAK(ii,jj,kk,rhoa,qr,qs,qg,qh,Ndr_in,Nds_in,Ndg_in,Ndh_in,rhos,rhog,rhoh, &
                         fracqrs,fracqrg,fracqrh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
                         rhomh,fwgbin,fwhbin)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 6/27/2012
!
! Based on Youngsun Jung's fractionWater subroutine, this subroutine
! uses an alternate method to compute the liquid water fraction on melting
! ice.  It also handles the case where there are multiple ice species at a given
! grid point.
!
! MODIFIED: Dan Dawson 08/03/2012
!
! Rewrote the logic of water fraction diagnostics for graupel and hail.
! The critical water mass on melting graupel and hail is now integrated over the distribution
! based on the formula of Rasmussen and Heymsfield (1987, eqn. 6)
! for the assumed discrete wet graupel/hail distribution and compared
! with the available water from the rainfield.  The minimum of the available
! water and the critical water fraction is then used for graupel and hail.  This requires
! an iterative technique using a first guess for the wet graupel/hail mass.
! Any excess, which should be slight, is added back to the rain field later in the 
! code.  For snow, we are just punting for now and allowing up to a maximum of 75%
! water fraction distributed evenly amongst sizes.
!
! MODIFIED: Dan Dawson 12/18/2012
!
! Split off the subroutine to handle the Takahashi bin microphysics scheme.
! All handling of partitioning of water between rain and snow/graupel/hail is now handled
! exclusively in this subroutine so that no further adjustments should have to be made
! later in the code.  The same approach as before in which the mean volume diameter is
! preserved for rain during transfer from rain to hail is used.
! Note, mixing ratios are assumed to have been pre-calculated and passed in to the subroutine
!
! NOTE: Dan Dawson 06/10/2014
!
! This subroutine needs revisiting. It is not recommended to use it at the present time.
! A different approach to that of the bulk scheme is needed.
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
 USE rsa_table
 IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

 INTEGER :: ii,jj,kk
 REAL :: rhoa, qr, qs, qg, qh, rhos, rhog, rhoh
!  INTEGER, PARAMETER :: nd_r = 34         ! Number of bins in Takahashi rain category
!  INTEGER, PARAMETER :: nd_s = 21         ! Number of bins in Takahashi ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nk_s =  5         ! Number of thickness categories for ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nd_g = 45         ! Number of bins in Takahashi graupel category
!  INTEGER, PARAMETER :: nd_h = 45         ! Number of bins in Takahashi hail category
 REAL*8, DIMENSION (nd_r) :: Ndr_in 
 REAL*8, DIMENSION (nk_s,nd_s) :: Nds_in
 REAL*8, DIMENSION (nd_g) :: Ndg_in
 REAL*8, DIMENSION (nd_h) :: Ndh_in
 
 REAL*8, DIMENSION (nd_r) :: Ndrtmp
 REAL*8, DIMENSION (nk_s,nd_s) :: Ndstmp
 REAL*8, DIMENSION (nd_g) :: Ndgtmp
 REAL*8, DIMENSION (nd_h) :: Ndhtmp
 
 REAL,INTENT(OUT) :: fracqrs,fracqrg,fracqrh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,rhomh
 REAL :: qgwfac,qhwfac
 REAL,INTENT(OUT) :: fwgbin(nd_g),fwhbin(nd_h)
 
 REAL :: frs,frg,frh
 REAL :: qti ! Total ice mixing ratio (sum of snow, graupel, and hail)
 REAL :: qswcrit,qgwcrit,qhwcrit,qtwcrit,tmp1
 REAL :: intv,cx,mtot,mwcrit,diff
 INTEGER :: i,iter

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 frs = 0.
 frg = 0.
 frh = 0.
 fws = 0.
 fwg = 0.
 fwh = 0.
 fracqrs = 0.
 fracqrg = 0.
 fracqrh = 0.
 fms = 0.
 fmg = 0.
 fmh = 0.
 rhoms = 0.
 rhomg = 0.
 rhomh = 0.
 qgwfac = 0.
 qhwfac = 0.
 fwgbin = 0.
 fwhbin = 0.

 Ndrtmp = 0.0
 Ndstmp = 0.0
 Ndgtmp = 0.0

 qti = qs+qg+qh

!-----------------------------------------------------------------------
! Calculate the first guess fraction of rain water to be transfered to the ice 
! surface (fr), based on the ratio of qs/qg/qh and qti up to a maximum of 100%
! of the rain.
!-----------------------------------------------------------------------

 IF (qr > 0. .AND. qti > 0.0) THEN   
   frs = qs/qti
   frg = qg/qti
   frh = qh/qti
 ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------
 IF(qs > 0.) fracqrs = MIN(frs*qr,3.*frs*qs) ! This gives a maximum of 75% water fraction for snow
 
 IF(qg > 0.0 .and. frg > 0.) THEN ! Iterate to find the water fraction on graupel
   fracqrg = frg*qr
   qgwcrit = 0.0
   diff = fracqrg-qgwcrit
   iter=0
   DO WHILE(diff > 0.0)
     qgwcrit = 0.0
     fmg = qg+fracqrg
     IF(fmg > 0.) THEN
       fwg = fracqrg/fmg
     ELSE
       fwg = 0.0
     END IF
     rhomg = 1000.*fwg**2. + (1.-fwg**2.)*rhog
     DO i = 1,nd_g
       intv = (dsg(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
       ! For each graupel diameter bin, first add water equally, preserving mean-volume diameter, and then adjust for maximum allowed
       Ndgtmp(i) = Ndg_in(i) + Ndg_in(i)*fracqrg/qg
       !Ndgtmp(i) = Ndg_in(i) ! Need to decided how to adjust number concentration in each bin...
       if(dsg(i) <= 9.0) THEN ! Try to completely melt graupel less than or equal to 8 mm
         !tmpfwg = fwh
         fwgbin(i) = 1.0
       elseif(dsg(i) <= 20.0) THEN ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula
         ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
         !mtot = (1./6.)*pi*(dsg(i)*1.0e-3)**3.*MAX(910.,rhomg) ! Mass of wet graupel of diameter dsg(i): treat it as high density even if it isn't
         mtot = (1./6.)*pi*(dsg(i)*1.0e-3)**3.*910.	! Rasmussen and Heymsfield (1987) assume 910 kg/m^3 density for this formula
         !mice = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhoh
         mwcrit = (0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
         fwgbin(i) = mwcrit/mtot ! Max water fraction for wet hailstone of diameter dsg(i)
       else
       	fwgbin(i) = 0.1
       endif
       qgwcrit = qgwcrit + intv*fwgbin(i)*Ndgtmp(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*1000./rhoa ! 2nd term is mass mixing ratio of liquid water in 
                                                                                      ! current diameter bin
     ENDDO 
     diff = fracqrg-qgwcrit
     !IF(ii == 60 .and. jj == 7) THEN
     !  write(*,*) 'iter,fracqrg,qgwcrit,fwg',iter,fracqrg,qgwcrit,fwg
     !ENDIF
     IF(iter == 0 .and. diff < 0) THEN ! water fraction to be added less than critical
       ! Redistribute available water evenly across diameter bins
       qgwfac = fracqrg/qgwcrit
       DO i = 1,nd_g
         fwgbin(i) = qgwfac*fwgbin(i) 
       ENDDO
       EXIT
     ELSE
       fracqrg = qgwcrit
       iter=iter+1
     ENDIF
   ENDDO
   qgwfac = MIN(fracqrg/qgwcrit,1.0)
   DO i = 1,nd_g
     Ndg_in(i) = Ndgtmp(i)
   ENDDO
 ENDIF
 
 IF(qh > 0. .and. frh > 0.) THEN ! Iterate to find the water fraction on hail
   fracqrh = frh*qr
   qhwcrit = 0.0
   diff = fracqrh-qhwcrit
   iter=0
   DO WHILE(diff > 0.0)
     qhwcrit = 0.0
     fmh = qh+fracqrh
     IF(fmh > 0.) THEN
       fwh = fracqrh/fmh
     ELSE
       fwh = 0.0
     END IF
     rhomh = 1000.*fwh**2. + (1.-fwh**2.)*rhoh
     Do i = 1,nd_h
       intv = (dsh(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
       ! For each hail diameter bin, first add water equally, preserving mean-volume diameter, and then adjust for maximum allowed
       Ndhtmp(i) = Ndh_in(i) + Ndh_in(i)*fracqrh/qh
       !Ndhtmp(i) = Ndh_in(i) ! Need to decided how to adjust number concentration in each bin...
       if(dsh(i) <= 9.0) THEN ! Try to completely melt hail less than or equal to 8 mm
         fwhbin(i) = 1.0
       else ! Use the Rasmussen and Heymsfield (1987) critical water mass on ice formula
         ! Mwcrit = 0.268 + 0.1389*Mice in grams, where 0.268 g is the mass of an 8mm spherical water drop
         mtot = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*MAX(910.,rhomh) ! Mass of wet hailstone of diameter dsh(i): treat it as high density even if it isn't
         mtot = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*910
         !mice = (1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhoh
         mwcrit = (0.268e-3 + 0.1389*mtot)/1.1389  ! Critical water mass before shedding
         fwhbin(i) = mwcrit/mtot ! Max water fraction for wet hailstone of diameter dsg(i)
       endif
       qhwcrit = qhwcrit + intv*fwhbin(i)*Ndhtmp(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*1000./rhoa ! rhomh ! 2nd term is mass mixing ratio of liquid water in 
                                                                                         ! current diameter bin
     ENDDO 
     diff = fracqrh-qhwcrit
!     IF(diff > 0.0) THEN
!       write(*,*) 'iter,fracqrg,qgwcrit,fwg',iter,fracqrg,qgwcrit,fwg
!     ENDIF
     IF(iter == 0 .and. diff < 0) THEN ! water fraction to be added less than critical
       ! Redistribute available water evenly across diameter bins
       qhwfac = fracqrh/qhwcrit
       DO i = 1,nd_h
         fwhbin(i) = qhwfac*fwhbin(i)
       ENDDO
       EXIT
     ELSE
       fracqrh = qhwcrit
       iter=iter+1
     ENDIF
   ENDDO 
   qhwfac = MIN(fracqrh/qhwcrit,1.0)
   DO i = 1,nd_h
     Ndh_in(i) = Ndhtmp(i)
   ENDDO
 ENDIF
 fms = qs+fracqrs
 fmg = qg+fracqrg
 fmh = qh+fracqrh
 
 IF(fms > 0.) THEN
   fws = fracqrs/fms
 ELSE
   fws = 0.
 ENDIF
 
 IF(fmg > 0.) THEN
   fwg = fracqrg/fmg
 ELSE
   fwg = 0.0
 END IF
 
 IF(fmh > 0.) THEN
   fwh = fracqrh/fmh
 ELSE
   fwh = 0.0
 ENDIF

 ! Remove water that was added above to snow,graupel, and hail from rain
 
 DO i = 1,nd_r
   Ndr_in(i) = Ndr_in(i) - Ndr_in(i)*(fracqrs+fracqrg+fracqrh)/qr
 ENDDO
 
 qr = qr-fracqrs-fracqrg-fracqrh

 ! Check water fraction of graupel, hail, and snow in each diameter bin shared with rain.  If it is (nearly) 100%
 ! add that bin to the corresponding rain bin
 
! DO i = 1,nd_r
!   IF(fwgbin(i) > 0.99) THEN
!     intv = (dsg(i)*1.0e-3)/4.329 
!     Ndr_in(i) = Ndr_in(i) + Ndg_in(i)
!     fmg = fmg - intv*Ndg_in(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*1000./rhoa
!     qr = qr + intv*Ndg_in(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*1000./rhoa
!     Ndg_in(i) = 0.0
!   ENDIF
!   IF(fwhbin(i) > 0.99) THEN
!     intv = (dsh(i)*1.0e-3)/4.329 
!     Ndr_in(i) = Ndr_in(i) + Ndh_in(i)
!     Ndh_in(i) = 0.0
!     fmh = fmh - intv*Ndh_in(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*1000.
!     qr = qr + intv*Ndh_in(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*1000.
!   ENDIF
! END DO
 
 rhoms = 1000.*fws**2. + (1.-fws**2.)*rhos
 rhomg = 1000.*fwg**2. + (1.-fwg**2.)*rhog
 rhomh = 1000.*fwh**2. + (1.-fwh**2.)*rhoh

END SUBROUTINE fractionWaterTAK

!########################################################################
!########################################################################
!#########                                                      #########
!#########              SUBROUTINE refl_rsa                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

#ifdef COMMAS
SUBROUTINE refl_rsa (addtorain,ii,jj,kk,MPflg,MFflg,rhoa,fws,fwg,fwh,qs,qg,qh,qrf,   &
                     nts,ntg,nth,ntr,alfr,alfs,alfg,alfh,rhoms,rhomg,rhomh,tair_C,   &
                     T_log_ref, T_sum_ref_h, T_sum_ref_v, &
                     T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
                     T_Ahh,     T_Avv,  &
                     T_qr, T_qs, T_qh, T_qg, &
                     T_Ntr, T_Nts, T_Nth, T_Ntg, &
                     T_alfr,T_alfs,T_alfh,T_alfg, &
                     T_mur, T_mus, T_mug, T_muh, &
                     T_qsw, T_qgw, T_qhw, rhog, rhoh, &
                     nrbin,nsbin,rhosbin,ngbin,rhogbin,nhbin,rhohbin,fwsbinout,      &
                     fwgbinout,fwhbinout)

#else
SUBROUTINE refl_rsa (addtorain,ii,jj,kk,MPflg,MFflg,rhoa,fws,fwg,fwh,qs,qg,qh,qrf,   &
                     nts,ntg,nth,ntr,alfr,alfs,alfg,alfh,rhoms,rhomg,rhomh,tair_C,   &
                     nrbin,nsbin,rhosbin,ngbin,rhogbin,nhbin,rhohbin,fwsbinout,      &
                     fwgbinout,fwhbinout)
#endif
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates specific differential phase.
! Compute radar observations by integrating radar scattering amplitudes
! over the drop size range. Radar scattering amplitues were calculated
! using T-matrix method and stored in the table.
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
! MODIFIED: Dan Dawson, 2/08/2011
!           Changed from a function to a subroutine.
!           Put derived type variables explicitly in argument list.
!           Removed references to derived type and instead refer 
!           explicitly to the variables within (this is done because f2py doesn't
!           work with fortran derived types).
!
!           Dan Dawson, 2/14/2011
!           Modified DSD integration to allow for gamma distribution of
!           the form n(D) = N0*D^alpha*exp(-lambda*D^(3*mu)),
!           that is, an additional shape parameter, mu, as in the ZVD
!           scheme. Also included the case (through the MFflg variable)
!           where the melted water fraction on ice is explicitly predicted
!           (the ZVDM scheme variants).  
!
!           Dan Dawson, 01/17/2012
!           Added MFhail flag to allow for special treatment of small melting
!           hail, by converting it to rain.
!
!           Dan Dawson, 12/04/2012
!           Removed MFhail flag and associated parameters and code, modified MFflg 
!           to control the method by which water fraction is diagnosed:
!             MFflg = 0: Water fraction on snow, graupel, and hail is diagnosed
!                        using the original method of Jung et al. (2008)
!             MFflg = 1: Water fraction on snow, graupel, and hail is diagnosed
!                        using the method of Dawson et al. (2012)
!             MFflg = 2: Water fraction on snow, graupel, and hail are explicitly
!                        predicted by the microphysics or at least passed into the
!                        subroutine from the external calling routine.
!
!           Youngsun Jung, 8/12/2014
!           Added temperature-based melting layer option (MFflg = 3 and 4)
!             MFflg = 3 Water fraction on snow, graupel, and hail are diagnosed
!                       based on the air temperature. 
!                       See subroutine fractionWater_temperature for details.
!             MFflg = 4 Similar to MFflg = 3 but dry and wet scattering amplitudes
!                       are averaged as a quadratic function of fractionWater.
!
!           Dan Dawson, 01/26/2015-02/05/2015
!           More modifications to fractionWater4.  Also renamed to fractionWater14, and
!           removed fractionWater2.  Renamed fractionWater3 to fractionWaterD14_old.  
!           Added several optional output arrays holding binned PSD information for diagnostics. WIP.  
!
!           Dan Dawson, 05/13/2015
!           General code cleanup
!    
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER :: ii,jj,kk     ! Grid indices for debugging/diagnostics
  INTEGER :: MPflg,MFflg
  REAL, PARAMETER :: pi4 = 97.409           ! pi^4
  REAL*8, DIMENSION (nd_r) :: Ndr, Ndrtmp 
  REAL*8, DIMENSION (nd_s) :: Nds, Ndrs 
  REAL*8, DIMENSION (nd_h) :: Ndh, Ndrh  
  REAL*8, DIMENSION (nd_g) :: Ndg, Ndrg
	REAL :: fwgbin(nd_g),rhomgbin(nd_g)
  REAL :: fwhbin(nd_h),rhomhbin(nd_h)

  REAL*8, EXTERNAL :: gamma_r8
  INTEGER, EXTERNAL :: get_qgh_opt

  REAL :: rhoa,cx,tair_C
  REAL :: qr,qs,qh,qg
  REAL :: ntr,nts,nth,ntg,ntrold
  REAL*8 :: alfr,alfs,alfh,alfg
  REAL*8 :: mur,mus,mug,muh
  REAL*8 :: qsw,qgw,qhw
  REAL*8 :: db_N0,db_N02,Ntw,Ntd
  REAL :: fracqrs,fracqrh,fracqrg,qsoakg,qsoakh,qgwfac,qhwfac
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg
  REAL :: fws,fwh,fwg,tmpqgw,tmpfwg,mtot,mwcrit,mice
  REAL :: rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf,qrfold,qrtmp
  REAL*8 :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL :: tsar_h,tsas_h,tsah_h,tsag_h,tsars_h,tsarh_h,tsarg_h
  REAL :: tsar_v,tsas_v,tsah_v,tsag_v,tsars_v,tsarh_v,tsarg_v
  COMPLEX :: tsar_hv,tsas_hv,tsah_hv,tsag_hv,tsars_hv,tsarh_hv,tsarg_hv
  REAL :: tfsar,tfsas,tfsah,tfsag,tfsars,tfsarh,tfsarg
  REAL :: D,intv,far,far2,lambda4
  REAL :: fa2,fb2,temph,tempv,temphv,tempk,temp,temp1,temp2,fa2h,fb2h
  real :: temph_rain, tempv_rain
  COMPLEX :: fab,fba,fconj, fabh,fbah,fconjh
  INTEGER :: i, idx
  REAL :: tempAhh,tempAvv
  REAL :: Ar_h,As_h,Ag_h,Ah_h,Ars_h,Arg_h,Arh_h
  REAL :: Ar_v,As_v,Ag_v,Ah_v,Ars_v,Arg_v,Arh_v
  
  REAL :: tem1,tem2,D0r
  
  real :: frac,finterp

  LOGICAL :: ranout, samesizebinsg, samesizebinsh, addtorain, adjustNt
  
  ! New optional output arrays containing binned PSD information
  REAL, OPTIONAL, INTENT(INOUT) :: nrbin(:,:),nsbin(:,:),ngbin(:,:),nhbin(:,:)
  REAL, OPTIONAL, INTENT(INOUT) :: rhosbin(:,:),rhogbin(:,:),rhohbin(:,:)
  REAL, OPTIONAL, INTENT(INOUT) :: fwsbinout(:),fwgbinout(:),fwhbinout(:)
#ifdef COMMAS
! TYPE T_obs_dual
  real :: rhog, rhoh
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw
#endif
  REAL, PARAMETER :: qmin = 1.e-12 ! Threshold for minimum mixing ratio for diagnostic water routine
  logical :: lflag
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

  adjustNt = .true. ! Adjust number concentration of hail and graupel during water transfers
		    ! to maintain mean volume diameter?

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  lambda4 = lambda**4.

  qr = T_qr
  qs = T_qs
  qh = T_qh
  qg = T_qg
  ntr = T_Ntr
  nts = T_Nts
  nth = T_Nth
  ntg = T_Ntg
  alfr = T_alfr
  alfs = T_alfs
  alfh = T_alfh
  alfg = T_alfg

  ntrold = 0.0
  qrfold = 0.0

  mur = T_mur  !
  mus = T_mus  ! New shape parameters (to accomodate gamma-in-volume distribution for ZVD)
  mug = T_mug  !
  muh = T_muh  !

  qsw = T_qsw  !
  qgw = T_qgw  ! New variables holding liquid water fraction
  qhw = T_qhw  ! mixing ratios

  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  Ndr = 0.; Nds = 0.; Ndh = 0.; Ndg = 0.
  Ndrs = 0.; Ndrh = 0.; Ndrg = 0.
  Ndrtmp = 0.; qrtmp = 0.
  temph = 0.; tempv = 0.; temphv = 0.; temp = 0.; tempk = 0.
  temph_rain = 0. ; tempv_rain = 0.0
  tempAhh = 0.; tempAvv = 0.
  fracqrs = 0.; fracqs = 0.; fms = 0.; fws = 0.; rhoms = rhos  ! 100.
  fracqrh = 0.; fracqh = 0.; fmh = 0.; fwh = 0.; rhomh = rhoh  ! 913.
  fracqrg = 0.; fracqg = 0.; fmg = 0.; fwg = 0.; rhomg = rhog  ! 400.
  T_sum_ref_h = 0.
  T_sum_ref_v = 0.
  T_log_zdr = missing
  T_log_zdr_rain = missing
  T_log_ref = 0.
  T_sum_ref_hv = 0.
  T_kdp = missing
  T_Ahh = 0.
  T_Avv = 0.

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

  ! Check to see if we are using the same size bins for graupel and hail
  ! The diagnosed water fraction is handled differently for this case
  temp1 = (dsg(2)-dsg(1))*1.0e-3
  temp2 = (dsr(2)-dsr(1))*1.0e-3
  samesizebinsg = (temp1 == temp2)
	
  temp1 = (dsh(2)-dsh(1))*1.0e-3
  temp2 = (dsr(2)-dsr(1))*1.0e-3
  samesizebinsh = (temp1 == temp2)
  
  fwgbin = 0.
  fwhbin = 0.
  rhomgbin = 0.
  rhomhbin = 0.

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
! For the variable definition, see "FUNCTION rainIceRefl".
! DTD: In the case of the ZVDM scheme variants, the liquid water mixing
! ratio on ice is explicitly predicted
!-----------------------------------------------------------------------

  lflag = .true.
!  IF ( (qg+qh) > qmin .and. qr > qmin ) lflag = (qr)/(qr+qg+qh) > 0.7 ! mf1b70
! IF ( qg > qmin .and. qr > qmin ) lflag = (qr)/(qr+qg) > 0.90 ! mf1b90g
!  IF ( (qg+qh) > qmin .and. qr > qmin ) lflag = (qr)/(qr+qg+qh) > 0.5 ! mf1b50
!  lflag = ( qg > qmin .and. rhog > 870.0 ) .or. ( qh > qmin .and. rhoh > 870.0 )
  
  IF (MFflg == 0) THEN       ! Diagnose water fraction on snow, graupel, and hail
                             ! based on method of Jung et al. (2008)

    addtorain = .false.

    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hl_ON == 1)  &
    CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
    IF(grpl_ON == 1) &
    CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)
    qrf = qr - fracqrs - fracqrh - fracqrg
    ! DTD: adjust rain number concentration to preserve mean diameter
    IF(ntr > 0.0) ntr = ntr - ntr*(fracqrs+fracqrh+fracqrg)/qr
    if(qrf < 0.0) qrf = 0.0

    qsf = qs - fracqs
    if(qsf < 0.0) qsf = 0.0
    qhf = qh - fracqh
    if(qhf < 0.0) qhf = 0.0
    qgf = qg - fracqg
    if(qgf < 0.0) qgf = 0.0
    
    qr = qrf
    qs = qsf+fms
    qg = qgf+fmg
    qh = qhf+fmh
    
    fwgbin(:) = fwg
    fwhbin(:) = fwh
    rhomgbin(:) = rhomg
    rhomhbin(:) = rhomh

  ELSEIF (MFflg == 1 .and. lflag) THEN       ! New diagnostic method by DTD (06/27/2012) as used in
                                 ! Dawson et al. (2014)
    
    ! First store "untouched" graupel and hail PSD bin numbers and densities (rain is handled later)
    IF(present(ngbin)) THEN ! Graupel
      ngbin = 0.0
      if(qg > 0.) then
        intv = (dsg(2) - dsg(1))*1.e-3
        Ntd = 0.
        if(ntg > 0.0) then
          CALL cal_N0(rhoa,qg,ntg,rhog,alfg,db_N0,mug)
          CALL cal_lamda(rhoa,qg,ntg,rhog,alfg,lamg,mug)
          Ntd = ntg
        else
          db_N0 = dble(N0g)
          cx = rhog*(pi/6.)
          CALL cal_Nt(rhoa,qg,db_N0,cx,alfg,ntg,mug)
          CALL cal_lamda(rhoa,qg,ntg,rhog,alfg,lamg,mug)
          Ntd = ntg
        endif
        DO i = 1,nd_g
          ngbin(1,i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamg*(dsg(i)*1.e-3)**(3.0*mug))*intv
          Ntd = Ntd - ngbin(1,i)
          if(Ntd <= 0) ngbin(1,i) = 0.
          rhogbin(1,i) = rhog
        ENDDO
      endif
    ENDIF
    
    IF(present(nhbin)) THEN ! Hail
      nhbin = 0.0
      if(qh > 0.) then
        intv = (dsh(2) - dsh(1))*1.e-3
        Ntd = 0.
        if(nth > 0.0) then
          CALL cal_N0(rhoa,qh,nth,rhoh,alfh,db_N0,muh)
          CALL cal_lamda(rhoa,qh,nth,rhoh,alfh,lamh,muh)
          Ntd = nth
        else
          db_N0 = dble(N0h)
          cx = rhoh*(pi/6.)
          CALL cal_Nt(rhoa,qh,db_N0,cx,alfh,nth,muh)
          CALL cal_lamda(rhoa,qh,nth,rhoh,alfh,lamh,muh)
          Ntd = nth
        endif
        DO i = 1,nd_h
          nhbin(1,i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamh*(dsh(i)*1.e-3)**(3.0*muh))*intv
          Ntd = Ntd - nhbin(1,i)
          if(Ntd <= 0) nhbin(1,i) = 0.
          rhohbin(1,i) = rhoh
        ENDDO
      endif
    ENDIF
        
!    IF ( qg > qmin .and. qr > qmin .and. (qr+0.1*qmin)/(qr+qg+0.1*qmin) > 0.90 ) THEN
    
    CALL fractionWaterD14(MPflg,ii,jj,kk,fos,rhoa,tair_C,qr,qs,qg,qh,ntr,nts,ntg,nth,alfr,alfg,alfh,mug,muh,rhos,rhog,rhoh, &
                          fracqrs,fracqs,fracqrg,fracqrh,qsoakg,qsoakh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
                          rhomh,Ndrg,Ndrh,fwgbin,fwhbin,rhomgbin,rhomhbin)

!    ENDIF
    
    qrfold = qr ! Save old mixing ratio and number conc. of rain in case it gets wiped out and needs
                ! to be added back later, else we lose the mean diameter information
    if(ntr > 0.0) ntrold = ntr 
    qrf = qr - fracqrs - fracqrg - fracqrh - qsoakg - qsoakh
    if(qrf <= 0.0) then
      qrf = 0.0
      ntr = 0.0
    endif
    ! DTD: adjust rain number concentration to preserve mean diameter
    IF(ntr > 0.0) then
      ntr = ntr - ntr*(fracqrs+fracqrh+fracqrg+qsoakg+qsoakh)/qr
    endif
 
    ! Now reverting to Jung et al. (2008) method for snow only
    qsf = qs - fracqs
    if(qsf < 0.0) qsf = 0.0
    qs = qsf+fms
 
    if(fwg == 0.0) then
      qgf = fmg
      fmg = 0.0
      qg = qgf
    else
      qgf = 0.0
      qg = fmg
    endif
    if(fwh == 0.0) then
      qhf = fmh
      fmh = 0.0
      qh = qhf
    else
      qhf = 0.0
      qh = fmh
    endif
    qr = qrf


  ELSEIF (MFflg == 3) THEN      ! Temperature-based melting.
    
    CALL fractionWater_temperature(qr,qs,rhos,fms,fws,rhoms,tair_C)
    IF(hl_ON == 1)  &
    CALL fractionWater_temperature(qr,qh,rhoh,fmh,fwh,rhomh,tair_C)
    IF(grpl_ON == 1)  &
    CALL fractionWater_temperature(qr,qg,rhog,fmg,fwg,rhomg,tair_C)

    qrf = qr
    qsf = qs-fms
    qhf = qh-fmh
    qgf = qg-fmg
    
    fwgbin(:) = fwg
    fwhbin(:) = fwh
    rhomgbin(:) = rhomg
    rhomhbin(:) = rhomh

  ELSEIF (MFflg == 4) THEN      ! Temperature-based melting option 2.

    CALL fractionWater_temperature(qr,qs,rhos,fms,fws,rhoms,tair_C)
    IF(hl_ON == 1)  &
    CALL fractionWater_temperature(qr,qh,rhoh,fmh,fwh,rhomh,tair_C)
    IF(grpl_ON == 1)  &
    CALL fractionWater_temperature(qr,qg,rhog,fmg,fwg,rhomg,tair_C)

    qrf = qr
    qsf = qs
    qhf = qh
    qgf = qg
    
    fwgbin(:) = fwg
    fwhbin(:) = fwh
    rhomgbin(:) = rhomg
    rhomhbin(:) = rhomh

  ELSE  ! MFflg == 2 : Melted fraction explicitly predicted (passed into subroutine)
   
    ! Note, the handling of this case is somewhat different than the above
    ! Above, there is an explicit mixing ratio of a rain/snow mixture = fms
    ! while qsf represents the mixing ratio of the nonmelting dry snow.  In the ZVD
    ! case, it is assumed that if melting is occuring, all of the snow
    ! at a given point is melting: thus, fms = qs, and qsf = 0.
    ! Similarly for the other ice hydrometeors.
    ! For rain, since the water fraction on ice is carried separately,
    ! qrf = qr

!    IF(kk == 1 .and. ii == 130 .and. jj == 142) THEN
!        write(*,*) 'before: qr,ntr,qh,fmh,nth,rhomh',qr,ntr,qh,fmh,nth,rhomh
!    ENDIF

    qrf = qr
    qgwfac = 1.0
    qhwfac = 1.0

    qrfold = qr ! Save old mixing ratio and number conc. of rain in case it gets wiped out and needs
                ! to be added back later, else we lose the mean diameter information
    if(ntr > 0.0) ntrold = ntr 

    IF(qsw > 0.0) THEN
      qsf = 0.0
      fms = qs
      fracqrs = qsw
      fracqs = qs-qsw
      IF(qs > 0.0) THEN
        fws = qsw/qs
      ELSE 
        fws = 0.0
      END IF
      ! For now just assume rhoms is the same as original dry snow (i.e. constant)
      rhoms = rhos        
    ELSE
      qsf = qs
      fms = 0.0
      fracqrs = 0.0
      fracqs = qs
      fws = 0.0
      rhoms = rhos
    END IF
    IF(qgw > 0.0) THEN
      qgf = 0.0
      fmg = qg
      fracqrg = qgw
      fracqg = qg-qgw
      IF(qg > 0.0) THEN
        fwg = Min(1.0, qgw/qg)
      ELSE
        fwg = 0.0
      END IF
      ! The (variable) density of graupel should already take into account melting
      rhomg = rhog
    ELSE
      qgf = qg
      fmg = 0.0
      fracqrg = 0.0
      fracqg = qg
      fwg = 0.0
      rhomg = rhog
    END IF
    IF(qhw > 0.0) THEN
      IF(qh > 0.0) THEN
        fwh = qhw/qh
      ELSE
        fwh = 0.0
      END IF
      qhf = 0.0
      fmh = qh
      fracqrh = qhw
      fracqh = qh-qhw
      ! The (variable) density of hail should already take into account melting
      rhomh = rhoh
    ELSE
      qhf = qh
      fmh = 0.0
      fracqrh = 0.0
      fracqh = qh
      fwh = 0.0
      rhomh = rhoh
    END IF
    
    qrtmp = 0.0 ! Will hold excess water mixing ratio
    temp1 = fracqrg
    temp2 = fracqrh
    

    addtorain = ( fwg < 1.0 .and. fwg > 0.0 ) .or. ( fwh < 1.0 .and. fwh > 0.0 ) ! 

!    addtorain = .false.
    IF(addtorain) THEN    
      CALL adjustWater(ii,jj,kk,rhoa,qr,qs,qg,qh,ntr,nts,ntg,nth,alfg,alfh,mug,muh,rhos,rhog,rhoh, &
                          fracqrs,fracqrg,fracqrh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
                          rhomh,fwgbin(:),fwhbin(:))

!      fwgbin(:) = fwg
!      fwhbin(:) = fwh

    ELSE
      fwgbin(:) = fwg
      fwhbin(:) = fwh
    ENDIF
    
    rhomgbin(:) = rhomg
    rhomhbin(:) = rhomh
    
    ! Add up excess water from wet graupel and hail.  We'll deal with it later.
    ! The mixing ratio/number concentration of graupel/hail should already have been adjusted in adjustWater
    qrtmp = temp1-fracqrg
    qrtmp = qrtmp+(temp2-fracqrh)
    
  END IF

!-----------------------------------------------------------------------
! Calculate N(d) for each diameter bin for each species based on input
! DSD information.
!-----------------------------------------------------------------------

  db_N0 = 0.d0

  ! Dry snow
  if(qsf > 0.) then
    intv = (dss(2) - dss(1))*1.e-3
    Ntd = 0.
    if(MPflg >= 3 .and. MPflg <= 4) then ! UM (3) or Thompson (4) microphysics case -- bimodal snow distribution
      ! Compute shape and slope parameters for snow
      cx = 0.069
      CALL power_mom(MPflg,2,cx,tair_C,qsf,tem1)
      CALL power_mom(MPflg,3,cx,tair_C,qsf,tem2)
      lamdas1 = (tem1/tem2)*tmplamdas1
      lamdas2 = (tem1/tem2)*tmplamdas2
      N0s1 = sngl(((dble(tem1)**4)/(dble(tem2)**3))*dble(tmpN0s1))
      N0s2 = sngl(((dble(tem1)**4)/(dble(tem2)**3))*dble(tmpN0s2)*      &
                  ((dble(tem1)/dble(tem2))**dble(tmpalphas2)))
      db_N0 = dble(N0s1)
      db_N02 = dble(N0s2)
      Ntd = sngl((db_N0/lamdas1) + (db_N02/(lamdas2**tmpalphas2))*gamma_r8(1.d0+tmpalphas2))
      nts = Ntd 
      DO i = 1,nd
        Nds(i) = db_N0*exp(-lamdas1*(dss(i)*1.e-3))*  &
                    intv + db_N02*(dss(i)*1.e-3)**tmpalphas2*exp(-lamdas2*       &
                    (dss(i)*1.e-3))*intv
        Ntd = Ntd - Nds(i)
        if(Ntd <= 0) Nds(i) = 0. 
     enddo
    else 
      if(nts > 0.0) then
        CALL cal_N0(rhoa,qsf,nts,rhos,alfs,db_N0,mus)
        CALL cal_lamda(rhoa,qsf,nts,rhos,alfs,lams,mus)
        Ntd = nts
      else
        db_N0 = dble(N0s)
        cx = rhos*(pi/6.)
        CALL cal_Nt(rhoa,qsf,db_N0,cx,alfs,nts,mus)
        CALL cal_lamda(rhoa,qsf,nts,rhos,alfs,lams,mus)
        Ntd = nts
      endif
      Do i = 1,nd_s
        Nds(i) = db_N0*(dss(i)*1.e-3)**alfs*exp(-lams*(dss(i)*1.e-3)**(3.0*mus))*intv
        Ntd = Ntd - Nds(i)
        if(Ntd <= 0) Nds(i) = 0.
      ENDDO
      ! write(0,*) 'qs,nts,ntd = ',qs,nts,ntd
      nts = ntd
    endif
  endif

  db_N0 = 0.d0

  ! Wet snow
  if(fms > 0.) then
    intv = (dss(2) - dss(1))*1.e-3
    Ntw = 0.
    if(MPflg >= 3 .and. MPflg <= 4) then ! UM (3) or Thompson (4) microphysics case -- bimodal snow distribution
      ! Compute shape and slope parameters for snow
      cx = 0.069
      CALL power_mom(MPflg,2,cx,tair_C,fms,tem1)
      CALL power_mom(MPflg,3,cx,tair_C,fms,tem2)
      lamdas1 = (tem1/tem2)*tmplamdas1
      lamdas2 = (tem1/tem2)*tmplamdas2
      N0s1 = sngl(((dble(tem1)**4)/(dble(tem2)**3))*dble(tmpN0s1))
      N0s2 = sngl(((dble(tem1)**4)/(dble(tem2)**3))*dble(tmpN0s2)*      &
                  ((dble(tem1)/dble(tem2))**dble(tmpalphas2)))
      db_N0 = dble(N0s1)
      db_N02 = dble(N0s2)
      Ntw = sngl((db_N0/lamdas1) + (db_N02/(lamdas2**tmpalphas2))*gamma_r8(1.d0+tmpalphas2))
      nts = Ntw 
      DO i = 1,nd
        Ndrs(i) = db_N0*exp(-lamdas1*(dss(i)*1.e-3))*  &
                    intv + db_N02*(dss(i)*1.e-3)**tmpalphas2*exp(-lamdas2*       &
                    (dss(i)*1.e-3))*intv
        Ntw = Ntw - Nds(i)
        if(Ntw <= 0) Nds(i) = 0. 
     enddo
    else 
      if(nts > 0.0) then
        CALL cal_N0(rhoa,fms,nts,rhoms,alfs,db_N0,mus)
        CALL cal_lamda(rhoa,fms,nts,rhoms,alfs,lamrs,mus)
        Ntw = nts
      else
        db_N0 = dble(N0s)
        cx = rhoms*(pi/6.)
        CALL cal_Nt(rhoa,fms,db_N0,cx,alfs,nts,mus)
        CALL cal_lamda(rhoa,fms,nts,rhoms,alfs,lamrs,mus)
        Ntw = nts
      endif
      Do i = 1,nd_s
        Ndrs(i) = db_N0*(dss(i)*1.e-3)**alfs*exp(-lamrs*(dss(i)*1.e-3)**(3.0*mus))*intv
        Ntw = Ntw - Ndrs(i)
        if(Ntw <= 0) Ndrs(i) = 0.
      ENDDO
      nts = ntw
    endif
  endif

  IF(hl_ON == 1) THEN  
    db_N0 = 0.d0

    ! Dry hail
    if(qhf > 0.) then
      intv = (dsh(2) - dsh(1))*1.e-3
      Ntd = 0.
      if(nth > 0.0) then ! Multi-moment case (or Nt already computed for single-moment case)
        CALL cal_N0(rhoa,qhf,nth,rhoh,alfh,db_N0,muh)
        CALL cal_lamda(rhoa,qhf,nth,rhoh,alfh,lamh,muh)
        Ntd = nth
      else ! Single-moment case
        db_N0 = dble(N0h)
        cx = rhoh*(pi/6.)
        CALL cal_Nt(rhoa,qhf,db_N0,cx,alfh,nth,muh)
        CALL cal_lamda(rhoa,qhf,nth,rhoh,alfh,lamh,muh)
        Ntd = nth
      endif
      Do i = 1,nd_h
        Ndh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamh*(dsh(i)*1.e-3)**(3.0*muh))*intv
        Ntd = Ntd - Ndh(i)
        if(Ntd <= 0) Ndh(i) = 0.
      ENDDO
    endif

    !write(*,*) 'Computing hail DSD'
    db_N0 = 0.d0
    ! Wet hail
    if(fmh > 0.) then
      intv = (dsh(2) - dsh(1))*1.e-3
      Ntw = 0.
      IF(MFflg /= 1) THEN ! Don't need to do this for MFflg == 1, as it's already computed
        if(nth > 0.0) then ! Multi-moment case (or Nt already computed for single-moment case)
          CALL cal_N0(rhoa,fmh,nth,rhomh,alfh,db_N0,muh)
          CALL cal_lamda(rhoa,fmh,nth,rhomh,alfh,lamrh,muh)
        else ! Single-moment case
          db_N0 = dble(N0h)
          cx = rhomh*(pi/6.)
          CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,nth,muh)
          CALL cal_lamda(rhoa,fmh,nth,rhomh,alfh,lamrh,muh)
        endif
      ENDIF
      Ntw = nth
      Do i = 1,nd_h
        IF(MFflg /= 1) THEN
          Ndrh(i) = db_N0*(dble(dsh(i)*1.e-3))**alfh*exp(-lamrh*(dble(dsh(i)*1.e-3))**(3.d0*muh))*dble(intv)
          Ntw = Ntw - Ndrh(i)
          if(Ntw <= 0.d0) Ndrh(i) = 0.d0
        ENDIF
        
        temp = sngl(Ndrh(i))*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhomhbin(i)/rhoa ! Mass mixing ratio of bin
        
        ! Now figure out what to do with 100% water bins
        IF(addtorain .and. samesizebinsh .and. dsh(i) <= 8.0 .and. fwhbin(i) == 1.0) THEN ! Add directly to corresponding rain bin
        	Ndrtmp(i) = Ndrtmp(i) + Ndrh(i)
        ELSEIF(addtorain .and. fwhbin(i) == 1.0) THEN ! Add to excess water mixing ratio
        	qrtmp = qrtmp + temp
        ENDIF
        
        IF(fwhbin(i) == 1.0 .and. addtorain) THEN ! Zero out this bin since we just added it to rain
        	Ndrh(i) = 0.0
        	! Adjust number concentration and mixing ratio to account for this
        	nth = nth-Ndrh(i)
        	fmh = fmh-temp
        	qh = fmh
        ENDIF
        
        ! Now store the adjusted wet hail distribution
    
        IF(present(nhbin) .and. MFflg == 1) THEN
          nhbin(2,i) = Ndrh(i)
          IF(present(fwhbinout)) fwhbinout(i) = fwhbin(i)
          rhohbin(2,i) = rhomhbin(i)
        ENDIF
        
      ENDDO
      ! If we have Ntw left over, add it to the last bin to preserve total number concentration
!    	if(Ntw > 0.) then
!    		Ndrh(nd_h) = Ndrh(nd_h) + dble(Ntw)
!    	endif
    endif
  ENDIF 

  IF(grpl_ON == 1) THEN
    db_N0 = 0.d0

    ! Dry graupel
    if(qgf > 0.) then
      intv = (dsg(2) - dsg(1))*1.e-3
      Ntd = 0.
      if(ntg > 0.0) then ! Multi-moment case (or Nt already computed for single-moment case)
        CALL cal_N0(rhoa,qgf,ntg,rhog,alfg,db_N0,mug)
        CALL cal_lamda(rhoa,qgf,ntg,rhog,alfg,lamg,mug)
        Ntd = ntg
      else ! Single-moment case
        cx = rhog*(pi/6.)
        IF(MPflg == 3) THEN ! Diagnostic N0 for UM scheme
          CALL callamda_UM(rhoa,qgf,cx,nag,nbg,alfg,lamg)
          CALL calN0_UM(nag,nbg,lamg,db_N0)
        ELSEIF(MPflg == 4) THEN ! Diagnostic N0 for Thompson scheme
          CALL cal_D0(rhoa,qrf,ntr,1000.,alfr,D0r) ! Compute median volume diameter for rain (using original rain distribution)
          CALL calN0g_Thompson(tair_C+273.15,rhoa,qgf,cx,alfg,rhog,qrf,D0r,db_N0)
        ELSE ! Assume fixed N0
          db_N0 = dble(N0g)
        ENDIF
        CALL cal_Nt(rhoa,qgf,db_N0,cx,alfg,ntg,mug)
        IF(MPflg /= 3) CALL cal_lamda(rhoa,qgf,ntg,rhog,alfg,lamg,mug)
        Ntd = ntg
      endif
      Do i = 1,nd_g
        Ndg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamg*(dsg(i)*1.e-3)**(3.0*mug))*intv
        Ntd = Ntd - Ndg(i)
        if(Ntd <= 0) Ndg(i) = 0.
      ENDDO
    endif

    db_N0 = 0.d0
    ! Wet graupel
    if(fmg > 0.) then
      intv = (dsg(2) - dsg(1))*1.e-3
      Ntw = 0.
      IF(MFflg /= 1) THEN ! Don't need to do this for MFflg == 1, as it's already computed
        if(ntg > 0.0) then ! Multi-moment case (or Nt already computed for single-moment case)
          CALL cal_N0(rhoa,fmg,ntg,rhomg,alfg,db_N0,mug)
          CALL cal_lamda(rhoa,fmg,ntg,rhomg,alfg,lamrg,mug)
        else ! Single-moment case
          cx = rhomg*(pi/6.)
          IF(MPflg == 3) THEN ! Diagnostic N0 for UM scheme
            CALL callamda_UM(rhoa,fmg,cx,nag,nbg,alfg,lamrg)
            CALL calN0_UM(nag,nbg,lamrg,db_N0)
          ELSEIF(MPflg == 4) THEN ! Diagnostic N0 for Thompson scheme
            CALL cal_D0(rhoa,qrf,ntr,1000.,alfr,D0r) ! Compute median volume diameter for rain (using original rain distribution)
            CALL calN0g_Thompson(tair_C+273.15,rhoa,fmg,cx,alfg,rhomg,qrf,D0r,db_N0)
          ELSE ! Assume fixed N0
            db_N0 = dble(N0g)
          ENDIF
          CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,ntg,mug)
          IF(MPflg /= 3) CALL cal_lamda(rhoa,fmg,ntg,rhomg,alfg,lamrg,mug)
        endif
      ENDIF
      Ntw = ntg
      DO i = 1,nd_g 
        IF(MFflg /= 1) THEN 
          Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
          Ntw = Ntw - Ndrg(i)
          if(Ntw <= 0.d0) Ndrg(i) = 0.d0 
        ENDIF
 	  
        temp = sngl(Ndrg(i))*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhomgbin(i)/rhoa ! Mass mixing ratio of bin
        
        ! Now figure out what to do with 100% water bins
        IF(addtorain .and. samesizebinsg .and. dsg(i) <= 8.0 .and. fwgbin(i) == 1.0) THEN ! Add directly to corresponding rain bin
        	Ndrtmp(i) = Ndrtmp(i) + Ndrg(i)
        ELSEIF(addtorain .and. fwgbin(i) == 1.0) THEN ! Add to excess water mixing ratio
        	qrtmp = qrtmp + temp
        ENDIF
        
        IF(fwgbin(i) == 1.0 .and. addtorain) THEN ! Zero out this bin since we just added it to rain
        	Ndrg(i) = 0.0
        	! Adjust number concentration and mixing ratio to account for this
        	ntg = ntg-Ndrg(i)
        	fmg = fmg-temp
        	qg = fmg
        	!write(*,*) 'ii,jj,kk,qg',ii,jj,kk,qg
        ENDIF
        
        ! Now store the adjusted wet graupel distribution
    
        IF(present(ngbin) .and. MFflg == 1) THEN
          ngbin(2,i) = Ndrg(i)
          IF(present(fwgbinout)) fwgbinout(i) = fwgbin(i)
          rhogbin(2,i) = rhomgbin(i)
        ENDIF
        
      ENDDO
      ! If we have Ntw left over, add it to the last bin to preserve total number concentration
!   	if(Ntw > 0.) then
!   		Ndrg(nd_g) = Ndrg(nd_g) + Ntw
!   	endif
		
    endif
  ENDIF
  
  lamr = 0.d0; lams = 0.d0; lamh = 0.d0; lamrs = 0.d0; lamrh = 0.d0;
  db_N0 = 0.d0

  ! Compute rain N(d)
  
  ! First deal with the case with no existing rain at a grid point but potentially
  ! excess water to add to rain.
  ! We'll build the rain distribution using the 100%-water graupel and hail bins
  ! if they exist at this point.
    
  if(qrf == 0. .and. ntr == 0) then
    temp1 = 0.0
    DO i=1,nd_r
      IF(Ndrtmp(i) > 0.) THEN
        temp1 = temp1 + 1.
        Ndr(i) = Ndrtmp(i)
        ntr = ntr + Ndr(i)
        qrf = qrf + sngl(Ndr(i))*(1./6.)*pi*(dsr(i)*1.0e-3)**3.*1000./rhoa ! Mass mixing ratio of bin
      ENDIF
    ENDDO
    temp2 = 0.0
    ! Now add any additional excess water proportionally to each bin
    IF(qrtmp > 0. .and. ntr > 0. .and. qrf > 0.) THEN
      DO i = 1,nd_r
        Ndr(i) = Ndr(i)*(qrtmp+qrf)/qrf
        temp2 = temp2 + sngl(Ndr(i))*(1./6.)*pi*(dsr(i)*1.0e-3)**3.*1000./rhoa
      ENDDO
      ntr = ntr + qrtmp*ntr/qrf
      qrf = qrf + qrtmp
    ENDIF

  elseif(qrf > 0.) then ! if we already have rain
    intv = (dsr(2) - dsr(1))*1.e-3
    Ntw = 0.
    if(ntr > 0.0 .and. (MPflg < 3 .or. MPflg == 4 .or. MPflg == 5)) then ! Multi-moment case
      if(qrtmp > 0.) THEN ! should be zero if samesizebinsg/h is true
        ntr = ntr + qrtmp*ntr/qrf
        qrf = qrf+qrtmp
      endif
      CALL cal_N0(rhoa,qrf,ntr,rhor,alfr,db_N0,mur)
      CALL cal_lamda(rhoa,qrf,ntr,rhor,alfr,lamr,mur)
      Ntw = ntr
    else ! Single-moment case
      cx = rhor*(pi/6.)
      if(qrtmp > 0.) qrf = qrf+qrtmp
      IF(MPflg == 3) THEN ! UM case, diagnostic N0
        CALL callamda_UM(rhoa,qrf,cx,nar,nbr,alfr,lamr)
        CALL calN0_UM(nar,nbr,lamr,db_N0)
      ELSE
        db_N0 = dble(N0r)
      ENDIF
      CALL cal_Nt(rhoa,qrf,db_N0,cx,alfr,ntr,mur)  
      CALL cal_lamda(rhoa,qrf,ntr,rhor,alfr,lamr,mur)
      Ntw = ntr
    endif
    DO i = 1,nd_r
      Ndr(i) = db_N0*(dsr(i)*1.e-3)**alfr*exp(-lamr*(dsr(i)*1.e-3)**(3.0*mur))*intv
      Ntw = Ntw - Ndr(i)
      if(Ntw <= 0) Ndr(i) = 0.
      ! Now store the initial rain distribution (note this isn't true if samesizebinsg/h is false, the rain
      ! might have had qrtmp added to it...)
      IF(present(nrbin) .and. MFflg == 1) nrbin(1,i) = Ndr(i)
    ENDDO
    DO i = 1,nd_r
    ! Add extra water bins (if samesizebinsg/h is true)
      IF(Ndrtmp(i) > 0.) THEN
        Ndr(i) = Ndr(i) + Ndrtmp(i)
        ntr = ntr + Ndrtmp(i)
        qrf = qrf + sngl(Ndrtmp(i))*(1./6.)*pi*(dsr(i)*1.0e-3)**3.*1000./rhoa ! Mass mixing ratio of bin
      ENDIF
    ENDDO
    ! If we have Ntw left over, add it to the last bin to preserve total number concentration
    ! EDIT 02/04/2013: This is problematic for two moment ZVD with alpha = -0.4.  Need to revisit.
    ! Perhaps add extra Ntr evenly to all bins?
!    if(Ntw > 0.) then
!    	Ndr(nd) = Ndr(nd) + Ntw
!    endif

  endif

  ! Now store the adjusted rain distribution
  
  IF(present(nrbin) .and. MFflg == 1) THEN
    DO i=1,nd_r
      nrbin(2,i) = Ndr(i)
    ENDDO
  ENDIF

  qr = qrf ! Just to be safe!

!-----------------------------------------------------------------------
! Calculate radar observations.
!-----------------------------------------------------------------------
  tsar_h=0.; tsas_h=0.; tsah_h=0.; tsag_h=0.
  tsars_h=0.; tsarh_h=0.; tsarg_h=0.
  tsar_v=0.; tsas_v=0.; tsah_v=0.; tsag_v=0.
  tsars_v=0.; tsarh_v=0.; tsarg_v=0.
  tsar_hv=0.; tsas_hv=0.; tsah_hv=0.; tsag_hv = 0.
  tsars_hv=0.; tsarh_hv=0.; tsarg_hv=0.
  tfsar=0.; tfsas=0.; tfsah=0.; tfsag=0.
  tfsars=0.; tfsarh=0.; tfsarg=0.
  fa2=0.; fb2=0.; fab=0.; far=0.
  Ar_h=0.; As_h=0.; Ag_h=0.; Ah_h=0.
  Ars_h=0.; Arg_h=0.; Arh_h=0.
  Ar_v=0.; As_v=0.; Ag_v=0.; Ah_v=0.
  Ars_v=0.; Arg_v=0.; Arh_v=0.

  fa2=0.; fb2=0.; fab=0.; far=0.

  ! Dry snow
  if(qsf > 0.) then
    do i=1,nd_s
      fa2 = ABS(fas_b(i,1))**2
      fb2 = ABS(fbs_b(i,1))**2
      fab = fas_b(i,1)*CONJG(fbs_b(i,1))
      fba = fbs_b(i,1)*CONJG(fas_b(i,1))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsas_h = tsas_h + far*Nds(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsas_v = tsas_v + far*Nds(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsas_hv = tsas_hv + fconj*Nds(i)
      far=Cks*REAL(fas_f(i,1) - fbs_f(i,1))
      tfsas = tfsas + far*Nds(i)
    enddo
  endif

  fa2=0.; fb2=0.; fab=0.; far=0.

  ! Wet snow
  if(fms > 0.) then
    idx = INT(fws * 20 + 0.5) + 1
    do i=1,nd_s
      fa2 = ABS(fas_b(i,idx))**2
      fb2 = ABS(fbs_b(i,idx))**2
      fab = fas_b(i,idx)*CONJG(fbs_b(i,idx))
      fba = fbs_b(i,idx)*CONJG(fas_b(i,idx))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsars_h = tsars_h + far*Ndrs(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsars_v = tsars_v + far*Ndrs(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsars_hv = tsars_hv + fconj*Ndrs(i)
      far=Cks*REAL(fas_f(i,idx)-fbs_f(i,idx))
      tfsars = tfsars + far*Ndrs(i)
    enddo
  endif

  IF(hl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    ! Dry hail
    if(qhf > 0.) then

       frac = (rhoh - 500.)/(913.0 - 500.) ! yes, it can be negative

      IF ( .true. ) THEN
      do i=1,nd_g
        fa2 = ABS(fag_b(i,1))**2
        fb2 = ABS(fbg_b(i,1))**2
        fab = fag_b(i,1)*CONJG(fbg_b(i,1))
        fba = fbg_b(i,1)*CONJG(fag_b(i,1))

        fa2h = ABS(fah_b(i,1))**2
        fb2h = ABS(fbh_b(i,1))**2
        fabh = fah_b(i,1)*CONJG(fbh_b(i,1))
        fbah = fbh_b(i,1)*CONJG(fah_b(i,1))

        fa2 = fa2 + frac*(fa2h-fa2)
        fb2 = fb2 + frac*(fb2h-fb2)
        temp2 = REAL(fab) + frac*( REAL(fabh) - REAL(fab))
        
        far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*temp2)
        tsah_h = tsah_h + far*Ndh(i)
        far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*temp2)
        tsah_v = tsah_v + far*Ndh(i)
        
        
        fab = fab + frac*(fabh-fab)
        fba = fba + frac*(fbah-fba)
        
        fconj = Cgd*(fa2+fb2)+ Agd*fab+Bgd*fba 
        
        tsah_hv = tsah_hv + fconj*Ndh(i)
        far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
        far2=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsah = tfsah + (far + frac*(far2-far))*Ndh(i)
      enddo
       
     ELSE
      do i=1,nd_h
        fa2 = ABS(fah_b(i,1))**2
        fb2 = ABS(fbh_b(i,1))**2
        fab = fah_b(i,1)*CONJG(fbh_b(i,1))
        fba = fbh_b(i,1)*CONJG(fah_b(i,1))
        far=(Ahd*fa2 + Bhd*fb2 + 2*Chd*REAL(fab))
        tsah_h = tsah_h + far*Ndh(i)
        far=(Bhd*fa2 + Ahd*fb2 + 2*Chd*REAL(fab))
        tsah_v = tsah_v + far*Ndh(i)
        fconj=(Chd*(fa2+fb2)+Ahd*fab+Bhd*fba)
        tsah_hv = tsah_hv + fconj*Ndh(i)
        far=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsah = tfsah + far*Ndh(i)
      enddo
    endif
    
    ENDIF

    fa2=0.; fb2=0.; fab=0.; far=0.
    ! Wet hail
    if(fmh > 0.) then
      do i=1,nd_h  
        idx = INT(fwhbin(i) * 20 + 0.5) + 1
        CALL coeff_hail(fwhbin(i),fmh)      
              
        fa2 = ABS(fah_b(i,idx))**2
        fb2 = ABS(fbh_b(i,idx))**2
        fab = fah_b(i,idx)*CONJG(fbh_b(i,idx))
        fba = fbh_b(i,idx)*CONJG(fah_b(i,idx))
        far=(Ah*fa2 + Bh*fb2 + 2*Ch*REAL(fab))
        tsarh_h = tsarh_h + far*sngl(Ndrh(i))
        far=(Bh*fa2 + Ah*fb2 + 2*Ch*REAL(fab))
        tsarh_v = tsarh_v + far*sngl(Ndrh(i))
        fconj=(Ch*(fa2+fb2)+Ah*fab+Bh*fba)
        tsarh_hv = tsarh_hv + fconj*sngl(Ndrh(i))
        far=Ckh*REAL(fah_f(i,idx)-fbh_f(i,idx))
        tfsarh = tfsarh + far*sngl(Ndrh(i))
      enddo
    endif
  ENDIF 

  IF(grpl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    ! Dry graupel
    if(qgf > 0.) then
      
       frac = (rhog - 500.)/(913.0 - 500.) ! yes, it can be negative

      IF ( .true. ) THEN
      do i=1,nd_g
        fa2 = ABS(fag_b(i,1))**2
        fb2 = ABS(fbg_b(i,1))**2
        fab = fag_b(i,1)*CONJG(fbg_b(i,1))
        fba = fbg_b(i,1)*CONJG(fag_b(i,1))

        fa2h = ABS(fah_b(i,1))**2
        fb2h = ABS(fbh_b(i,1))**2
        fabh = fah_b(i,1)*CONJG(fbh_b(i,1))
        fbah = fbh_b(i,1)*CONJG(fah_b(i,1))

        fa2 = fa2 + frac*(fa2h-fa2)
        fb2 = fb2 + frac*(fb2h-fb2)
        temp2 = REAL(fab) + frac*( REAL(fabh) - REAL(fab))
        
        far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*temp2)
        tsag_h = tsag_h + far*Ndg(i)
        far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*temp2)
        tsag_v = tsag_v + far*Ndg(i)
        
        
        fab = fab + frac*(fabh-fab)
        fba = fba + frac*(fbah-fba)
        
        fconj = Cgd*(fa2+fb2)+ Agd*fab+Bgd*fba 
        
        tsag_hv = tsag_hv + fconj*Ndg(i)
        far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
        far2=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsag = tfsag + (far + frac*(far2-far))*Ndg(i)
      enddo
       
     ELSE
      IF ( rhog < 700. ) THEN
      do i=1,nd_g
        fa2 = ABS(fag_b(i,1))**2
        fb2 = ABS(fbg_b(i,1))**2
        fab = fag_b(i,1)*CONJG(fbg_b(i,1))
        fba = fbg_b(i,1)*CONJG(fag_b(i,1))
        far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*REAL(fab))
        tsag_h = tsag_h + far*Ndg(i)
        far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*REAL(fab))
        tsag_v = tsag_v + far*Ndg(i)
        fconj=(Cgd*(fa2+fb2)+Agd*fab+Bgd*fba)
        tsag_hv = tsag_hv + fconj*Ndg(i)
        far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
        tfsag = tfsag + far*Ndg(i)
      enddo
      ELSE ! use hail for high density
      do i=1,nd_g
        fa2 = ABS(fah_b(i,1))**2
        fb2 = ABS(fbh_b(i,1))**2
        fab = fah_b(i,1)*CONJG(fbh_b(i,1))
        fba = fbh_b(i,1)*CONJG(fah_b(i,1))
        far=(Ahd*fa2 + Bhd*fb2 + 2*Chd*REAL(fab))
        tsag_h = tsag_h + far*Ndg(i)
        far=(Bhd*fa2 + Ahd*fb2 + 2*Chd*REAL(fab))
        tsag_v = tsag_v + far*Ndg(i)
        fconj=(Chd*(fa2+fb2)+Ahd*fab+Bhd*fba)
        tsag_hv = tsag_hv + fconj*Ndg(i)
        far=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsag = tfsag + far*Ndg(i)
      enddo
      ENDIF
      
      ENDIF
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.

    !Wet graupel
    if(fmg > 0.) then
      do i=1,nd_g
        idx = INT(fwgbin(i) * 20 + 0.5) + 1
        CALL coeff_grpl(fwgbin(i),fmg)
                
        fa2 = ABS(fag_b(i,idx))**2
        fb2 = ABS(fbg_b(i,idx))**2
        fab = fag_b(i,idx)*CONJG(fbg_b(i,idx))
        fba = fbg_b(i,idx)*CONJG(fag_b(i,idx))
        far=(Ag*fa2 + Bg*fb2 + 2*Cg*REAL(fab))
        tsarg_h = tsarg_h + far*Ndrg(i)
        far=(Bg*fa2 + Ag*fb2 + 2*Cg*REAL(fab))
        tsarg_v = tsarg_v + far*Ndrg(i)
        fconj=(Cg*(fa2+fb2)+Ag*fab+Bg*fba)
        tsarg_hv = tsarg_hv + fconj*Ndrg(i)
        far=Ckg*REAL(fag_f(i,idx)-fbg_f(i,idx))
        tfsarg = tfsarg + far*Ndrg(i)
      enddo
    endif
  ENDIF

  ! Compute rain radar variables
  if(qrf > 0.) then
    do i=1,nd_r
      fa2 = ABS(far_b(i))**2
      fb2 = ABS(fbr_b(i))**2
      fab = far_b(i)*CONJG(fbr_b(i))
      fba = fbr_b(i)*CONJG(far_b(i))
      tsar_h = tsar_h + fa2*Ndr(i)
      tsar_v = tsar_v + fb2*Ndr(i)
      tsar_hv = tsar_hv + fab*Ndr(i)
      far=REAL(far_f(i) - fbr_f(i))
      tfsar = tfsar + far*Ndr(i)
    enddo
  endif

  IF (MFflg == 4) THEN      ! Temperature-based melting.
    tsas_h = (1-fws**2)*tsas_h + tsars_h*fws**2
    tsars_h = 0.0
    tsag_h = (1-fwg**2)*tsag_h + tsarg_h*fwg**2
    tsarg_h = 0.0
    tsah_h = (1-fwh**2)*tsah_h + tsarh_h*fwh**2
    tsarh_h = 0.0
    tsas_v = (1-fws**2)*tsas_v + tsars_v*fws**2
    tsars_v = 0.0
    tsag_v = (1-fwg**2)*tsag_v + tsarg_v*fwg**2
    tsarg_v = 0.0
    tsah_v = (1-fwh**2)*tsah_v + tsarh_v*fwh**2
    tsarh_v = 0.0
    tsas_hv = (1-fws**2)*tsas_hv + tsars_hv*fws**2
    tsars_h = 0.0
    tsag_hv = (1-fwg**2)*tsag_hv + tsarg_hv*fwg**2
    tsarg_hv = 0.0
    tsah_hv = (1-fwh**2)*tsah_hv + tsarh_hv*fwh**2
    tsarh_hv = 0.0
    tfsas = (1-fws**2)*tfsas + tfsars*fws**2
    tfsag = (1-fwg**2)*tfsag + tfsarg*fwg**2
    tfsah = (1-fwh**2)*tfsah + tfsarh*fwh**2
    tfsars = 0.0
    tfsarg = 0.0
    tfsarh = 0.0
  END IF

  temph = 4*lambda4/(pi4*Kw2)*(tsar_h+tsas_h+tsah_h+tsag_h+tsars_h+tsarh_h+tsarg_h)
  tempv = 4*lambda4/(pi4*Kw2)*(tsar_v+tsas_v+tsah_v+tsag_v+tsars_v+tsarh_v+tsarg_v)
 
  temph_rain = 4*lambda4/(pi4*Kw2)*(tsar_h)
  tempv_rain = 4*lambda4/(pi4*Kw2)*(tsar_v)
 
  T_sum_ref_h = temph
  T_sum_ref_v = tempv

  temphv = 4*lambda4/(pi4*Kw2)*ABS(tsar_hv+tsas_hv+tsah_hv+tsag_hv+tsars_hv+tsarh_hv+tsarg_hv) 
  T_sum_ref_hv = temphv

  if(temph > 0.) T_log_ref = 10*log10(temph)


  if(tempv > 0.) then
    T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
  endif

  if(tempv_rain > 0.) then
    T_log_zdr_rain = 10.*LOG10(MAX(1.0,temph_rain/tempv_rain))
  endif

!    IF ( MFflg == 0 ) T_log_zdr = T_log_zdr_rain
  
!JYS  if(tempk < 0.) tempk = 0.0
  tempk = 180.*lambda/pi*(tfsar+tfsas+tfsah+tfsag+tfsars+tfsarh+tfsarg)*1.e-3
  T_kdp = tempk

!   IF ( tair_C < 4.0 .and. T_qs > 1.e-4 ) THEN
!     write(0,*) 't,qs,cx,tref = ',tair_C,t_qs,qs,nts,ntd,t_log_ref,tfsas,rhoms,rhos,rhoh,rhomh
!   ENDIF

  ! For Jung et al. 2008 melting, adjust water fractions so that they represent the 
  ! water fraction with respect to the total ice, and not just the "melting ice", for 
  ! consistency of interpretation with the Dawson et al. 2014 formulation (MFflg = 1)
  IF ( qs > 1.e-13 ) THEN
    fws = fws*fms/qs
  ELSE
    fws = 0
  ENDIF
  IF(MFflg == 0) THEN 
    fwg = fwg*fmg/qg
    fwh = fwh*fmh/qh
  ENDIF
  
END SUBROUTINE refl_rsa

#ifdef COMMAS
SUBROUTINE refl_rsa_tak (ii,jj,kk,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,MFflg,rhoa,fws,fwg,fwh,qs,qg,qh,qrf,   &
                     Nds_in,Ndg_in,Ndh_in,Ndr_in,rhoms,rhomg,rhomh, &
                     T_log_ref, T_sum_ref_h, T_sum_ref_v, &
                     T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
                     T_Ahh,     T_Avv,  &
                     T_qr, T_qs, T_qh, T_qg, &
                     T_Ntr, T_Nts, T_Nth, T_Ntg, &
                     T_alfr,T_alfs,T_alfh,T_alfg, &
                     T_mur, T_mus, T_mug, T_muh, &
                     T_qsw, T_qgw, T_qhw, rhog, rhoh)
#else
SUBROUTINE refl_rsa_tak (ii,jj,kk,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,MFflg,rhoa,fws,fwg,fwh,qs,qg,qh,qrf,   &
                     Nds_in,Ndg_in,Ndh_in,Ndr_in,rhoms,rhomg,rhomh)
#endif
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates specific differential phase.
! Compute radar observations by integrating radar scattering amplitudes
! over the drop size range. Radar scattering amplitues were calculated
! using T-matrix method and stored in the table.
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
! MODIFIED: Dan Dawson, 2/08/2011
!           Changed from a function to a subroutine.
!           Put derived type variables explicitly in argument list.
!           Removed references to derived type and instead refer 
!           explicitly to the variables within (this is done because f2py doesn't
!           work with fortran derived types).
!
!           Dan Dawson, 2/14/2011
!           Modified DSD integration to allow for gamma distribution of
!           the form n(D) = N0*D^alpha*exp(-lambda*D^(3*mu)),
!           that is, an additional shape parameter, mu, as in the ZVD
!           scheme. Also included the case (through the MFflg variable)
!           where the melted water fraction on ice is explicitly predicted
!           (the ZVDM scheme variants).  
!
!           Dan Dawson, 01/17/2012
!           Added MFhail flag to allow for special treatment of small melting
!           hail, by converting it to rain.
!
!           Dan Dawson, 12/04/2012
!           Removed MFhail flag and associated parameters and code, modified MFflg 
!           to control the method by which water fraction is diagnosed:
!             MFflg = 0: Water fraction on snow, graupel, and hail is diagnosed
!                        using the original method of Jung et al. (2008)
!             MFflg = 1: Water fraction on snow, graupel, and hail is diagnosed
!                        using the method of Dawson et al. (2012)
!             MFflg = 2: Water fraction on snow, graupel, and hail are explicitly
!                        predicted by the microphysics or at least passed into the
!                        subroutine from the external calling routine.
!
!           Dan Dawson, 12/12/2012
!           Split off from subroutine refl_rsa.  This version handles specifically
!           the Takahashi bin scheme
!    
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER :: ii,jj,kk     ! Grid indices for debugging/diagnostics
  INTEGER :: nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,MFflg
!  INTEGER, PARAMETER :: nd_r = 34         ! Number of bins in Takahashi rain category
!  INTEGER, PARAMETER :: nd_s = 21         ! Number of bins in Takahashi ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nk_s =  5         ! Number of thickness categories for ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nd_g = 45         ! Number of bins in Takahashi graupel category
!  INTEGER, PARAMETER :: nd_h = 45         ! Number of bins in Takahashi hail category
  REAL, PARAMETER :: pi4 = 97.409           ! pi^4
  REAL*8, DIMENSION (nd_r_in) :: Ndr_in, Ndr 
  REAL*8, DIMENSION (nk_s_in,nd_s_in) :: Nds_in ! , Nds, Ndrs Not using these at the moment, but later need to update the code below
  REAL*8, DIMENSION (nd_s_in) :: Nds,Ndrs        ! For now keep them as 1D 
  REAL*8, DIMENSION (nd_g_in) :: Ndg_in, Ndg, Ndrg 
  REAL*8, DIMENSION (nd_h_in) :: Ndh_in, Ndh, Ndrh 
  REAL, DIMENSION (nd_g_in) :: fwgbin
  REAL, DIMENSION (nd_h_in) :: fwhbin
!  REAL*8, EXTERNAL :: gamma
  INTEGER, EXTERNAL :: get_qgh_opt

#ifdef COMMAS
! TYPE T_obs_dual
  real :: rhog, rhoh
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw
#endif

  REAL :: rhoa,cx
  REAL :: qr,qs,qh,qg
  REAL :: ntr,nts,nth,ntg,ntrold
  REAL*8 :: alfr,alfs,alfh,alfg
  !DTD add mu shape parameter
  REAL*8 :: mur,mus,mug,muh
  !DTD add liquid water fraction variables
  REAL*8 :: qsw,qgw,qhw
  REAL*8 :: db_N0,Ntw,Ntd
  REAL :: fracqrs,fracqrh,fracqrg,qgwfac,qhwfac
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg
  REAL :: fws,fwh,fwg,tmpqgw,tmpfwg,mtot,mwcrit
  REAL :: rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf,qrfold
  REAL*8 :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL :: tsar_h,tsas_h,tsah_h,tsag_h,tsars_h,tsarh_h,tsarg_h
  REAL :: tsar_v,tsas_v,tsah_v,tsag_v,tsars_v,tsarh_v,tsarg_v
  COMPLEX :: tsar_hv,tsas_hv,tsah_hv,tsag_hv,tsars_hv,tsarh_hv,tsarg_hv
  REAL :: tfsar,tfsas,tfsah,tfsag,tfsars,tfsarh,tfsarg
  REAL :: D,intv,far,lambda4
  REAL :: fa2,fb2,temph,tempv,temphv,tempk,temp,temp2,temphr,tempvr
  real :: temph_rain, tempv_rain
  COMPLEX :: fab,fba,fconj
  INTEGER :: i, idx
  REAL :: tempAhh,tempAvv
  REAL :: Ar_h,As_h,Ag_h,Ah_h,Ars_h,Arg_h,Arh_h
  REAL :: Ar_v,As_v,Ag_v,Ah_v,Ars_v,Arg_v,Arh_v

  LOGICAL :: ranout

!  TYPE(T_para_dsd) :: var_dsd

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    write(*,*) 'grpl_ON,hl_ON',grpl_ON,hl_ON
    write(*,*) 'qgh_opt',qgh_opt

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  lambda4 = lambda**4.

! qr = var_dsd%T_qr
! qs = var_dsd%T_qs
! qh = var_dsd%T_qh
! qg = var_dsd%T_qg
! ntr = var_dsd%T_Ntr
! nts = var_dsd%T_Nts
! nth = var_dsd%T_Nth
! ntg = var_dsd%T_Ntg
! alfr = var_dsd%T_alfr
! alfs = var_dsd%T_alfs
! alfh = var_dsd%T_alfh
! alfg = var_dsd%T_alfg

  qr = T_qr
  qs = T_qs
  qh = T_qh
  qg = T_qg
!  ntr = T_Ntr
!  nts = T_Nts
!  nth = T_Nth
!  ntg = T_Ntg
!  alfr = T_alfr
!  alfs = T_alfs
!  alfh = T_alfh
!  alfg = T_alfg

  ntr = 0.0
  nts = 0.0
  nth = 0.0
  ntg = 0.0
  alfr = 0.0
  alfs = 0.0
  alfh = 0.0
  alfg = 0.0

  ntrold = 0.0
  qrfold = 0.0

!  mur = T_mur  !
!  mus = T_mus  ! New shape parameters (to accomodate gamma-in-volume distribution for ZVD)
!  mug = T_mug  !
!  muh = T_muh  !

  mur = 0.0
  mus = 0.0
  mug = 0.0
  muh = 0.0

!  qsw = T_qsw  !
!  qgw = T_qgw  ! New variables holding liquid water fraction
!  qhw = T_qhw  ! mixing ratios

  qsw = 0.0
  qgw = 0.0
  qhw = 0.0

  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  Ndr = 0.; Nds = 0.; Ndh = 0.; Ndg = 0.
  Ndrs = 0.; Ndrh = 0.; Ndrg = 0.
  temph = 0.; tempv = 0.; temphv = 0.; temp = 0.; tempk = 0.
  temphr = 0.; tempvr = 0.
  temph_rain = 0. ; tempv_rain = 0.0
  tempAhh = 0.; tempAvv = 0.
  fracqrs = 0.; fracqs = 0.; fms = 0.; fws = 0.; rhoms = rhos  ! 100.
  fracqrh = 0.; fracqh = 0.; fmh = 0.; fwh = 0.; rhomh = rhoh  ! 913.
  fracqrg = 0.; fracqg = 0.; fmg = 0.; fwg = 0.; rhomg = rhog  ! 400.
  !refl_rsa = init_Refl()
  T_sum_ref_h = 0.
  T_sum_ref_v = 0.
  T_log_zdr = missing
  T_log_zdr_rain = missing
  T_log_ref = 0.
  T_sum_ref_hv = 0.
  T_kdp = missing
  T_Ahh = 0.
  T_Avv = 0.

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
! For the variable definition, see "FUNCTION rainIceRefl".
! NOT YET IMPLEMENTED for the Takahashi scheme.  All ice assumed dry for
! the moment.
!-----------------------------------------------------------------------


IF(MFflg == 2 .or. .true.) THEN

  ! Rain
  qrf = qr
  ! Snow
  fms = 0.0
  fws = 0.0
  ! Set snow mixing ratio to zero here for now
  qsf = 0.0
  ! Graupel
  fmg = 0.0
  fwg = 0.0
  qgf = qg
  ! Hail
  fmh = 0.0
  fwh = 0.0
  qhf = qh
  
!  IF ( MFflg == 2 ) THEN
!    grpl_ON = 0
!    hl_ON = 0
!  ENDIF

ELSE

  qrf = qr

  IF(ii == 60 .and. jj == 42) THEN
    write(*,*) 'ii,kk',ii,kk
    write(*,*) 'Before, qr,qg',qrf,qg
    write(*,*) 'Before, Ndr',(Ndr_in(i),i=1,nd_r)
    write(*,*) 'Before, Ndg',(Ndg_in(i),i=1,nd_g)
    write(*,*) 'Before, Ndh',(Ndh_in(i),i=1,nd_h)
    ntg = 0.0
    nth = 0.0
    DO i=1,nd_g
      intv = (dsg(i)*1.0e-3)/4.329
      ntg = ntg+Ndg_in(i)*intv
      nth = nth+Ndh_in(i)*intv
    ENDDO
    write(*,*) 'Before, dmg',(qg/ntg)**(1./3.)
    write(*,*) 'Before, dmh',(qh/nth)**(1./3.)
  ENDIF

!  CALL fractionWaterTAK(ii,jj,kk,rhoa,qrf,qs,qg,qh,Ndr_in,Nds_in,Ndg_in,Ndh_in,rhos,rhog,rhoh, &
!                          fracqrs,fracqrg,fracqrh,fms,fmg,fmh,fws,fwg,fwh,rhoms,rhomg,   &
!                          rhomh,fwgbin,fwhbin)
    
  IF(ii == 60 .and. jj == 42) THEN
    write(*,*) 'After, qr,qg',qrf,fmg
    temp = 0.0
    temp2 = 0.0
    ntg = 0.0
    nth = 0.0
    DO i=1,nd_g
      intv = (dsg(i)*1.0e-3)/4.329
      temp = temp + intv*Ndg_in(i)*(1./6.)*pi*(dsg(i)*1.0e-3)**3.*rhog
			temp2 = temp2 + intv*Ndh_in(i)*(1./6.)*pi*(dsh(i)*1.0e-3)**3.*rhoh
			ntg = ntg + Ndg_in(i)*intv
      nth = nth + Ndh_in(i)*intv
    END DO
    write(*,*) 'After, qg2,qh2',temp,temp2
    write(*,*) 'After, Ndr',(Ndr_in(i),i=1,nd_r)
    write(*,*) 'After, Ndg',(Ndg_in(i),i=1,nd_g)
    write(*,*) 'After, Ndh',(Ndh_in(i),i=1,nd_h)
    write(*,*) 'After, dmg',(qg/ntg)**(1./3.)
    write(*,*) 'After, dmh',(qh/nth)**(1./3.)
    write(*,*) 'fwg = ',fwg
    write(*,*) 'fwh = ',fwh
    write(*,*) 'fwgbin',(fwgbin(i),i=1,nd_g)
    write(*,*) 'fwhbin',(fwhbin(i),i=1,nd_h)
    write(*,*) 'fmg',fmg
    write(*,*) 'fmh',fmh
    write(*,*) 'rhomg,rhomh',rhomg,rhomh
  ENDIF
    
  if(fws == 0.0) then
    qsf = fms
    fms = 0.0
    qs = qsf
  else
    qsf = 0.0
    qs = fms
  endif
  if(fwg == 0.0) then
    qgf = fmg
    fmg = 0.0
    qg = qgf
  else
    qgf = 0.0
    qg = fmg
  endif
  if(fwh == 0.0) then
    qhf = fmh
    fmh = 0.0
    qh = qhf
  else
    qhf = 0.0
    qh = fmh
  endif
  qr = qrf
ENDIF


!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
! DTD: This is now done for each diameter bin below due to allowing
! water fraction to vary for each bin
!-----------------------------------------------------------------------
!  IF(hl_ON == 1) CALL coeff_hail(fwh,fmh)
!  IF(grpl_ON == 1) CALL coeff_grpl(fwg,fmg)

! Calculate diameter bins for graupel, hail, and rain based on input spectra
! Ignore snow for now since it's a bit more complicated (5 thickness bins in addition to
! 21 diameter bins)

  IF(hl_ON == 1) THEN
    ! Dry hail
    if(qhf > 0.) then
      DO i = 1,nd_h
        intv = (dsh(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
        !write(*,*) 'i,dsh(i),Ndh_in(i)',i,dsh(i),Ndh_in(i)
        Ndh(i) = Ndh_in(i)*intv
      END DO
    endif
    ! Wet Hail
    if(fmh > 0.) then
      DO i = 1,nd_g
        intv = (dsh(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
        Ndrh(i) = Ndh_in(i)*intv
      ENDDO
    endif
  ENDIF

  IF(grpl_ON == 1) THEN
    ! Dry graupel
    if(qgf > 0.) then
      DO i = 1,nd_g
        intv = (dsg(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
        !write(*,*) 'i,dsg(i),Ndg_in(i)',i,dsg(i),Ndg_in(i)
        Ndg(i) = Ndg_in(i)*intv
      END DO
    endif
    ! Wet Graupel
    if(fmg > 0.) then
      DO i = 1,nd_g
        intv = (dsg(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
        Ndrg(i) = Ndg_in(i)*intv
      ENDDO
    endif
  ENDIF

  ! Rain
  !write(*,*) 'qrf = ',qrf
  if(qrf > 0.) then
    DO i = 1,nd_r
      intv = (dsr(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
      !write(*,*) 'i,dsr(i),Ndr_in(i)',i,dsr(i),Ndr_in(i)
      Ndr(i) = Ndr_in(i)*intv
    END DO
  endif

!-----------------------------------------------------------------------
! Calculate radar observations.
!-----------------------------------------------------------------------
  tsar_h=0.; tsas_h=0.; tsah_h=0.; tsag_h=0.
  tsars_h=0.; tsarh_h=0.; tsarg_h=0.
  tsar_v=0.; tsas_v=0.; tsah_v=0.; tsag_v=0.
  tsars_v=0.; tsarh_v=0.; tsarg_v=0.
  tsar_hv=0.; tsas_hv=0.; tsah_hv=0.; tsag_hv = 0.
  tsars_hv=0.; tsarh_hv=0.; tsarg_hv=0.
  tfsar=0.; tfsas=0.; tfsah=0.; tfsag=0.
  tfsars=0.; tfsarh=0.; tfsarg=0.
  fa2=0.; fb2=0.; fab=0.; far=0.
  Ar_h=0.; As_h=0.; Ag_h=0.; Ah_h=0.
  Ars_h=0.; Arg_h=0.; Arh_h=0.
  Ar_v=0.; As_v=0.; Ag_v=0.; Ah_v=0.
  Ars_v=0.; Arg_v=0.; Arh_v=0.

  fa2=0.; fb2=0.; fab=0.; far=0.

  ! rain
  if(qrf > 0.) then
    do i=1,nd_r
!    	IF(ii == 60 .and. jj == 42) THEN
!        write(*,*) 'i,dsr(i),Ndr(i)',i,dsg(i)*1.0e-3,Ndr(i) 
!      ENDIF
      fa2 = ABS(far_b(i))**2
      fb2 = ABS(fbr_b(i))**2
      fab = far_b(i)*CONJG(fbr_b(i))
      fba = fbr_b(i)*CONJG(far_b(i))
      tsar_h = tsar_h + fa2*Ndr(i)
      tsar_v = tsar_v + fb2*Ndr(i)
      tsar_hv = tsar_hv + fab*Ndr(i)
      far=REAL(far_f(i) - fbr_f(i))
      tfsar = tfsar + far*Ndr(i)
    enddo
  endif

  ! Dry snow
  if(qsf > 0.) then
    do i=1,nd_s
      fa2 = ABS(fas_b(i,1))**2
      fb2 = ABS(fbs_b(i,1))**2
      fab = fas_b(i,1)*CONJG(fbs_b(i,1))
      fba = fbs_b(i,1)*CONJG(fas_b(i,1))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsas_h = tsas_h + far*Nds(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsas_v = tsas_v + far*Nds(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsas_hv = tsas_hv + fconj*Nds(i)
      far=Cks*REAL(fas_f(i,1) - fbs_f(i,1))
      tfsas = tfsas + far*Nds(i)
    enddo
  endif

  fa2=0.; fb2=0.; fab=0.; far=0.

  ! Wet snow
  if(fms > 0.) then
    idx = INT(fws * 20 + 0.5) + 1
    do i=1,nd_s
      fa2 = ABS(fas_b(i,idx))**2
      fb2 = ABS(fbs_b(i,idx))**2
      fab = fas_b(i,idx)*CONJG(fbs_b(i,idx))
      fba = fbs_b(i,idx)*CONJG(fas_b(i,idx))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsars_h = tsars_h + far*Ndrs(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsars_v = tsars_v + far*Ndrs(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsars_hv = tsars_hv + fconj*Ndrs(i)
      far=Cks*REAL(fas_f(i,idx)-fbs_f(i,idx))
      tfsars = tfsars + far*Ndrs(i)
    enddo
  endif

  IF(hl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    ! Dry hail
    if(qhf > 0.) then
      do i=1,nd_h
        fa2 = ABS(fah_b(i,1))**2
        fb2 = ABS(fbh_b(i,1))**2
        fab = fah_b(i,1)*CONJG(fbh_b(i,1))
        fba = fbh_b(i,1)*CONJG(fah_b(i,1))
        far=(Ahd*fa2 + Bhd*fb2 + 2*Chd*REAL(fab))
        tsah_h = tsah_h + far*Ndh(i)
        far=(Bhd*fa2 + Ahd*fb2 + 2*Chd*REAL(fab))
        tsah_v = tsah_v + far*Ndh(i)
        fconj=(Chd*(fa2+fb2)+Ahd*fab+Bhd*fba)
        tsah_hv = tsah_hv + fconj*Ndh(i)
        far=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsah = tfsah + far*Ndh(i)
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.
    ! Wet hail
    if(fmh > 0.) then
      !idx = INT(fwh * 20 + 0.5) + 1
      ! Greatly simplified this section.  Now uses fwhbin, which is the diagnosed water fraction
      ! computed earlier for each diameter bin in fractionWaterTAK.  No further adjustment
      ! is done.  This approach will be applied soon to the bulk scheme version (refl_rsa and fractionwater3)
      do i=1,nd_h
        idx = INT(fwhbin(i) * 20 + 0.5) + 1
        
        CALL coeff_hail(fwhbin(i),fmh)      
!        IF(ii == 60 .and. jj == 42) THEN
!          write(*,*) 'i,dsh(i),fwhbin(i),Ndrh(i)',i,dsh(i)*1.0e-3,fwhbin(i),Ndrh(i) 
!        ENDIF
        
        fa2 = ABS(fah_b(i,idx))**2
        fb2 = ABS(fbh_b(i,idx))**2
        fab = fah_b(i,idx)*CONJG(fbh_b(i,idx))
        fba = fbh_b(i,idx)*CONJG(fah_b(i,idx))
        far=(Ah*fa2 + Bh*fb2 + 2*Ch*REAL(fab))
        tsarh_h = tsarh_h + far*Ndrh(i)
        far=(Bh*fa2 + Ah*fb2 + 2*Ch*REAL(fab))
        tsarh_v = tsarh_v + far*Ndrh(i)
        fconj=(Ch*(fa2+fb2)+Ah*fab+Bh*fba)
        tsarh_hv = tsarh_hv + fconj*Ndrh(i)
        far=Ckh*REAL(fah_f(i,idx)-fbh_f(i,idx))
        tfsarh = tfsarh + far*Ndrh(i)
      enddo
    endif
  ENDIF 

  IF(grpl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    ! Dry graupel
    if(qgf > 0.) then
      do i=1,nd_g
        fa2 = ABS(fag_b(i,1))**2
        fb2 = ABS(fbg_b(i,1))**2
        fab = fag_b(i,1)*CONJG(fbg_b(i,1))
        fba = fbg_b(i,1)*CONJG(fag_b(i,1))
        far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*REAL(fab))
        tsag_h = tsag_h + far*Ndg(i)
        far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*REAL(fab))
        tsag_v = tsag_v + far*Ndg(i)
        fconj=(Cgd*(fa2+fb2)+Agd*fab+Bgd*fba)
        tsag_hv = tsag_hv + fconj*Ndg(i)
        far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
        tfsag = tfsag + far*Ndg(i)
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.

    !Wet graupel
    if(fmg > 0.) then
      !idx = INT(fwg * 20 + 0.5) + 1
      ! Greatly simplified this section.  Now uses fwgbin, which is the diagnosed water fraction
      ! computed earlier for each diameter bin in fractionWaterTAK.  No further adjustment
      ! is done.  This approach will be applied soon to the bulk scheme version (refl_rsa and fractionwater3)
      do i=1,nd_g
        idx = INT(fwgbin(i) * 20 + 0.5) + 1
        CALL coeff_grpl(fwgbin(i),fmg)
!        IF(ii == 60 .and. jj == 42) THEN
!          write(*,*) 'i,dsg(i),fwgbin(i),Ndrg(i)',i,dsg(i)*1.0e-3,fwgbin(i),Ndrg(i) 
!        ENDIF
        fa2 = ABS(fag_b(i,idx))**2
        fb2 = ABS(fbg_b(i,idx))**2
        fab = fag_b(i,idx)*CONJG(fbg_b(i,idx))
        fba = fbg_b(i,idx)*CONJG(fag_b(i,idx))
        far=(Ag*fa2 + Bg*fb2 + 2*Cg*REAL(fab))
        tsarg_h = tsarg_h + far*Ndrg(i)
        far=(Bg*fa2 + Ag*fb2 + 2*Cg*REAL(fab))
        tsarg_v = tsarg_v + far*Ndrg(i)
        fconj=(Cg*(fa2+fb2)+Ag*fab+Bg*fba)
        tsarg_hv = tsarg_hv + fconj*Ndrg(i)
        far=Ckg*REAL(fag_f(i,idx)-fbg_f(i,idx))
        tfsarg = tfsarg + far*Ndrg(i)
      enddo
    endif
  ENDIF

  temph = 4*lambda4/(pi4*Kw2)*(tsar_h+tsas_h+tsah_h+tsag_h+tsars_h+tsarh_h+tsarg_h)
  tempv = 4*lambda4/(pi4*Kw2)*(tsar_v+tsas_v+tsah_v+tsag_v+tsars_v+tsarh_v+tsarg_v)

  temph_rain = 4*lambda4/(pi4*Kw2)*(tsar_h)
  tempv_rain = 4*lambda4/(pi4*Kw2)*(tsar_v)
  !refl_rsa%T_sum_ref_h = temph
  !refl_rsa%T_sum_ref_v = tempv
  T_sum_ref_h = temph
  T_sum_ref_v = tempv

  !write(*,*) 'temph',temph
  !write(*,*) 'tempv',tempv

  temphv = 4*lambda4/(pi4*Kw2)*ABS(tsar_hv+tsas_hv+tsah_hv+tsag_hv+tsars_hv+tsarh_hv+tsarg_hv) 
  !refl_rsa%T_sum_ref_hv = temphv
  T_sum_ref_hv = temphv

  if(temph > 0.) T_log_ref = 10*log10(temph)
  
  if(tempv > 0.) then
!    refl_rsa%T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
    T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
!    if(kk == 1 .and. ii == 38 .and. jj == 31) then
!     write(*,*) 'T_log_zdr',T_log_zdr
!    endif
  endif

  if(tempv_rain > 0.) then
!    refl_rsa%T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
    T_log_zdr_rain = 10.*LOG10(MAX(1.0,temph_rain/tempv_rain))
!    IF ( MFflg == 2 ) T_log_zdr = T_log_zdr_rain
!    if(kk == 1 .and. ii == 38 .and. jj == 31) then
!     write(*,*) 'T_log_zdr',T_log_zdr
!    endif
  endif
  
!JYS  if(tempk < 0.) tempk = 0.0
  tempk = 180.*lambda/pi*(tfsar+tfsas+tfsah+tfsag+tfsars+tfsarh+tfsarg)*1.e-3
  !refl_rsa%T_kdp = tempk
  T_kdp = tempk
    
END SUBROUTINE refl_rsa_tak

END MODULE DUALPARA

#ifdef COMMAS
SUBROUTINE refl_rsa_array_takahashi(MFflg,nx,ny,nz,ng,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,ibgn,iend,jbgn,jend,kbgn,kend,rsafndir,wavelen,rhoa_arr,qr,qs,qg,qh,  &
           Ndr_in_arr,Nds_in_arr,Ndg_in_arr,Ndh_in_arr,qsw,qgw,qhw,rhogrpl,rhohail, &
           logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh,dpmask)
#else
SUBROUTINE refl_rsa_array_takahashi(MFflg,nx,ny,nz,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,ibgn,iend,jbgn,jend,kbgn,kend,rsafndir,wavelen,rhoa_arr,qr,qs,qg,qh,  &
           Ndr_in_arr,Nds_in_arr,Ndg_in_arr,Ndh_in_arr,qsw,qgw,qhw,rhogrpl,rhohail, &
           logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh)
#endif
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This is a wrapper subroutine for the subroutine refl_rsa_tak 
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 12/12/2012
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
#ifdef COMMAS
  USE dualpara, only: refl_rsa_tak,lambda
#else
  USE dualpara
#endif
  IMPLICIT NONE

  INTEGER :: MFflg
  INTEGER :: nx,ny,nz,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,ibgn,iend,jbgn,jend,kbgn(nx,ny),kend(nx,ny)
  REAL :: fws,fwg,fwh
  CHARACTER (LEN=256) :: rsafndir
  REAL, INTENT(IN) :: wavelen
  REAL, INTENT(IN) :: rhoa_arr(nx,ny,nz)

#ifdef COMMAS
  integer, intent(in) :: ng
#else
  integer, parameter :: ng = 0
#endif

!  REAL, INTENT(INOUT) :: qr(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng),qs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng), &
!                         qg(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng),qh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
  REAL, INTENT(INOUT) :: qr(nx,ny,nz),qs(nx,ny,nz), &
                         qg(nx,ny,nz),qh(nx,ny,nz)
!  INTEGER, PARAMETER :: nd_r = 34         ! Number of bins in Takahashi rain category
!  INTEGER, PARAMETER :: nd_s = 21         ! Number of bins in Takahashi ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nk_s =  5         ! Number of thickness categories for ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nd_g = 45         ! Number of bins in Takahashi graupel category
!  INTEGER, PARAMETER :: nd_h = 45         ! Number of bins in Takahashi hail category
#ifdef COMMAS
  REAL :: Ndr_in_arr(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nd_r_in) 
  REAL :: Nds_in_arr(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nd_s_in,nk_s_in)
  REAL :: Ndg_in_arr(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nd_g_in)
  REAL :: Ndh_in_arr(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nd_h_in) 
#else
  REAL :: Ndr_in_arr(nd_r_in,nx,ny,nz) 
  REAL :: Nds_in_arr(nk_s_in,nd_s_in,nx,ny,nz)
  REAL :: Ndg_in_arr(nd_g_in,nx,ny,nz)
  REAL :: Ndh_in_arr(nd_h_in,nx,ny,nz) 
#endif

  REAL*8, ALLOCATABLE :: Ndr_in(:),Nds_in(:,:),Ndg_in(:),Ndh_in(:)
  REAL, INTENT(INOUT) :: qsw(nx,ny,nz),qgw(nx,ny,nz),qhw(nx,ny,nz)
  REAL, INTENT(IN) :: rhogrpl(nx,ny,nz),rhohail(nx,ny,nz)
!f2py intent(in,out) qr,qs,qg,qh,qsw,qgw,qhw
  REAL, INTENT(OUT) :: logZ(nx,ny,nz),sumZh(nx,ny,nz),sumZv(nx,ny,nz),logZdr(nx,ny,nz),logZdrrain(nx,ny,nz)
  REAL, INTENT(OUT) :: sumZhv(nx,ny,nz),Kdp(nx,ny,nz),Ahh(nx,ny,nz),Avv(nx,ny,nz),rhv(nx,ny,nz)
  REAL, INTENT(OUT) :: rhoms(nx,ny,nz),rhomg(nx,ny,nz),rhomh(nx,ny,nz)
  REAL :: qsout,qgout,qhout,qrout,nsout,ngout,nhout,nrout,rhomsout,rhomgout,rhomhout
  REAL :: rhoa

#ifdef COMMAS
  integer, intent(inout) :: dpmask(nx,ny,nz)
  real :: rhog, rhoh
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw
#endif

  INTEGER :: i,j,k,a,b
  !INTEGER :: bin_opt
  
! Initialize radar wavelength

  lambda = wavelen

! Initialize DSD parameters

  !CALL init_dsd()

! DTD: Read RSA table here

  CALL set_nbins(1)

  ! Sanity check to make sure bin sizes match up with input arrays.
  
  IF(nd_r /= nd_r_in) THEN
    write(*,*) "Something is wrong!  Bin sizes don't match up! Exiting now!"
    write(*,*) "nd_r,nd_r_in",nd_r,nd_r_in
    STOP
  ENDIF
  
  ALLOCATE(Ndr_in(nd_r_in))
  ALLOCATE(Nds_in(nk_s,nd_s_in))
  ALLOCATE(Ndg_in(nd_g_in))
  ALLOCATE(Ndh_in(nd_h_in))

  CALL read_table(rsafndir)

  ! Set variables to zero that aren't needed
  
  T_mur = 0.0
  T_mus = 0.0
  T_mug = 0.0
  T_muh = 0.0
  T_Ntr = 0.0
  T_Nts = 0.0
  T_Ntg = 0.0
  T_Nth = 0.0
  T_alfr = 0.0
  T_alfs = 0.0
  T_alfg = 0.0
  T_alfh = 0.0
  T_qsw = 0.0
  T_qgw = 0.0
  T_qhw = 0.0

  fws = 0.0
  fwg = 0.0
  fwh = 0.0

#if 1
!$OMP PARALLEL DO DEFAULT(SHARED) schedule(dynamic, 1) &
!$OMP PRIVATE(T_log_ref, T_sum_ref_h, T_sum_ref_v, T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
!$OMP   T_Ahh, T_Avv, T_qr, T_qs, T_qh, T_qg, T_Ntr, T_Nts, T_Nth, T_Ntg, &
!$OMP   T_qsw, T_qgw, T_qhw, T_alfr,T_alfs,T_alfh,T_alfg, i, j, k, rhog, rhoh, rhoa, &
!$OMP   qsout,qgout, qhout,qrout,nsout,ngout,nhout,nrout,rhomsout, &
!$OMP   rhomgout,rhomhout,fws,fwg,fwh ,Ndr_in, Nds_in, Ndg_in, Ndh_in, a, b)
#endif
  DO j=jbgn,jend
    DO i=ibgn,iend
      DO k=kbgn(i,j),kend(i,j)
#ifdef COMMAS
        IF ( dpmask(i,j,k) == 0 ) THEN
           cycle
        ELSE
            dpmask(i,j,k) = 0 ! done with this point
        ENDIF
#endif

        ! Assign input variables
        T_qr = qr(i,j,k)
        T_qs = qs(i,j,k) 
        T_qh = qh(i,j,k)
        T_qg = qg(i,j,k)
        
        rhog = rhogrpl(i,j,k)
        rhoh = rhohail(i,j,k)
        rhoa = rhoa_arr(i,j,k)
        
      
#ifdef COMMAS
        DO b=1,nd_r
          Ndr_in(b) = dble(Ndr_in_arr(i,j,k,b))
          !write(*,*) 'Ndr_in(b)',b,Ndr_in(b)
        ENDDO
        
        DO a=1,nk_s
          DO b=1,nd_s
            Nds_in(a,b) = dble(Nds_in_arr(i,j,k,b,a))
            !write(*,*) 'Nds_in(a,b)',a,b,Nds_in(a,b)
          ENDDO
        ENDDO
        
        DO b=1,nd_g
          Ndg_in(b) = dble(Ndg_in_arr(i,j,k,b))
          !write(*,*) 'Ndg_in(b)',b,Ndg_in(b)
        ENDDO
        
        DO b=1,nd_h
          Ndh_in(b) = dble(Ndh_in_arr(i,j,k,b))
          !write(*,*) 'Ndh_in(b)',b,Ndh_in(b)
        ENDDO
        
        !write(*,*) 'i,j,k',i,j,k
        
        ! Compute the polarimetric variables
        CALL refl_rsa_tak (i,j,k,nd_r,nd_s,nk_s,nd_g,nd_h,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,qhout,qrout,   &
                     Nds_in,Ndg_in,Ndh_in,Ndr_in,rhomsout,rhomgout,rhomhout, &
                     T_log_ref, T_sum_ref_h, T_sum_ref_v, &
                     T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
                     T_Ahh,     T_Avv,  &
                     T_qr, T_qs, T_qh, T_qg, &
                     T_Ntr, T_Nts, T_Nth, T_Ntg, &
                     T_alfr,T_alfs,T_alfh,T_alfg, &
                     T_mur, T_mus, T_mug, T_muh, &
                     T_qsw, T_qgw, T_qhw, rhog, rhoh )
#else
        DO b=1,nd_r
          Ndr_in(b) = dble(Ndr_in_arr(b,i,j,k))
          !write(*,*) 'Ndr_in(b)',b,Ndr_in(b)
        ENDDO
        
        DO a=1,nk_s
          DO b=1,nd_s
            Nds_in(a,b) = dble(Nds_in_arr(a,b,i,j,k))
            !write(*,*) 'Nds_in(a,b)',a,b,Nds_in(a,b)
          ENDDO
        ENDDO
        
        DO b=1,nd_g
          Ndg_in(b) = dble(Ndg_in_arr(b,i,j,k))
          !write(*,*) 'Ndg_in(b)',b,Ndg_in(b)
        ENDDO
        
        DO b=1,nd_h
          Ndh_in(b) = dble(Ndh_in_arr(b,i,j,k))
          !write(*,*) 'Ndh_in(b)',b,Ndh_in(b)
        ENDDO
        
        !write(*,*) 'i,j,k',i,j,k
        
        ! Compute the polarimetric variables
        CALL refl_rsa_tak (i,j,k,nd_r,nd_s,nk_s,nd_g,nd_h,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,qhout,qrout,   &
                     Nds_in,Ndg_in,Ndh_in,Ndr_in,rhomsout,rhomgout,rhomhout)
#endif
        
        
        IF(MFflg < 2) THEN ! Compute water fraction mixing ratios from those calculated with fractionWater subroutine
          !if(i == 60 .and. j == 7) write(*,*) 'fwg,qg',fwg,qgout
          !write(*,*) 'i,j,k,fwh',i,j,k,fwh
          qsw(i,j,k) = fws*qsout
          qgw(i,j,k) = fwg*qgout
          qhw(i,j,k) = fwh*qhout
          
          ! Also adjust qr, qs, qg, qh appropriately EDIT: This is not needed because they are recalculated later in the python code
          qr(i,j,k) = qrout
          qs(i,j,k) = qsout
          qg(i,j,k) = qgout
          qh(i,j,k) = qhout
          
          rhoms(i,j,k) = rhomsout
          rhomg(i,j,k) = rhomgout
          rhomh(i,j,k) = rhomhout
        ELSE
          rhoms(i,j,k) = 100. ! Hardcoded for now.
          rhomg(i,j,k) = rhog
          rhomh(i,j,k) = rhoh
        ENDIF
        ! Store result in arrays for output
        logZ(i,j,k) = T_log_ref
        sumZh(i,j,k) = T_sum_ref_h
        sumZv(i,j,k) = T_sum_ref_v
        logZdr(i,j,k) = T_log_zdr
        logZdrrain(i,j,k) = T_log_zdr_rain
        sumZhv(i,j,k) = T_sum_ref_hv
        Kdp(i,j,k) = T_kdp
        Ahh(i,j,k) = T_Ahh
        Avv(i,j,k) = T_Avv
        ! Compute rho_hv
        IF(T_sum_ref_h*T_sum_ref_v > 0.) THEN
          rhv(i,j,k) = T_sum_ref_hv/                      &
            SQRT(T_sum_ref_h*T_sum_ref_v)
        ELSE
          rhv(i,j,k) = 0.0
        ENDIF
        
      END DO
    END DO
  END DO
      
  RETURN

END SUBROUTINE refl_rsa_array_takahashi

SUBROUTINE d0r_rsa_array_takahashi(MFflg,nx,ny,nz,ndr,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,ibgn,iend,jbgn,jend,slice, &
           sr,sx,sxdx, &
           rhoa_arr, &
           Ndr_in_arr,Nds_in_arr,Ndg_in_arr,Ndh_in_arr, &
           d0r)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculate the median volume diameter of rain
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Ted Mansell (7/17/2019), borrowing from Dan Dawson
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  USE dualpara
  IMPLICIT NONE

  INTEGER :: MFflg
  INTEGER :: nx,ny,nz,ndr,nd_r_in,nd_s_in,nk_s_in,nd_g_in,nd_h_in,ibgn,iend,jbgn,jend,slice ! ,kbgn(nx,ny),kend(nx,ny)
  REAL :: fws,fwg,fwh
!  CHARACTER (LEN=256) :: rsafndir
  REAL, INTENT(IN) :: rhoa_arr(nx,ny,nz)
!  REAL, INTENT(INOUT) :: qr(nx,ny,nz),qs(nx,ny,nz),qg(nx,ny,nz),qh(nx,ny,nz)
!  INTEGER, PARAMETER :: nd_r = 34         ! Number of bins in Takahashi rain category
!  INTEGER, PARAMETER :: nd_s = 21         ! Number of bins in Takahashi ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nk_s =  5         ! Number of thickness categories for ice/snow category (not used currently)
!  INTEGER, PARAMETER :: nd_g = 45         ! Number of bins in Takahashi graupel category
!  INTEGER, PARAMETER :: nd_h = 45         ! Number of bins in Takahashi hail category
  real, intent( in ) :: sr(ndr), sx(nd_r_in), sxdx(nd_r_in)
  REAL :: Ndr_in_arr(nd_r_in,nx,ny) 
  REAL :: Nds_in_arr(nk_s_in,nd_s_in,nx,ny)
  REAL :: Ndg_in_arr(nd_g_in,nx,ny)
  REAL :: Ndh_in_arr(nd_h_in,nx,ny) 
  REAL*8, ALLOCATABLE :: Ndr_in(:),Nds_in(:,:),Ndg_in(:),Ndh_in(:)
!f2py intent(out) d0r
  REAL, INTENT(OUT) :: d0r(nx,ny)
  REAL :: qsout,qgout,qhout,qrout,nsout,ngout,nhout,nrout,rhomsout,rhomgout,rhomhout
  REAL :: rhoa
  INTEGER :: i,j,k,a,b,n,il
  !INTEGER :: bin_opt
  real :: tmp,tmpn,d0rmax
  double precision :: masscum(45),qr,qrmax,chw
  

  ! Sanity check to make sure bin sizes match up with input arrays.
  
!  IF(nd_r /= nd_r_in) THEN
!    write(*,*) "Something's wrong!  Bin sizes don't match up! Exiting now!"
!    write(0,*) 'ndr, nd_r,nd_r_in',ndr,nd_r,nd_r_in
!    write(0,*) 'slice = ',slice
!    STOP
!  ENDIF
  
!  ALLOCATE(Ndr_in(nd_r_in))
!  ALLOCATE(Nds_in(nk_s,nd_s_in))
!  ALLOCATE(Ndg_in(nd_g_in))
!  ALLOCATE(Ndh_in(nd_h_in))


  ! Set variables to zero that aren't needed
  
    DO i = 1,nd_r_in
!      intv = (dsr(i)*1.0e-3)/4.329 ! Exponentially-varying bin width
!      write(0,*) 'i,sr,sx,sxdx',i,sr(i) ,sx(i),sxdx(i)
!      Ndr(i) = Ndr_in(i)*intv
    END DO

      qrmax = 0.0
      d0rmax = 0.0
      k=slice+1
  DO j=jbgn,jend
    DO i=ibgn,iend

                qr = 0.0d0
                chw = 0.0d0
                masscum(:) = 0.0

               DO n = 16,nd_r_in
!               DO n = 1,lfmax
                qr  = qr  +  Ndr_in_arr(n,i,j)*sxdx(n)*sx(n) ! total mass
                chw = chw +  Ndr_in_arr(n,i,j)*sxdx(n)
                masscum(n) = qr ! cumulative mass
               ENDDO
               
               qrmax = Max(qrmax,qr)

             IF ( qr > 0.0 ) THEN
                n = 0
          
                DO il = 16,nd_r_in
                  IF ( qr*0.5d0 < masscum(il) ) THEN
                    n = il
                    EXIT
                  ENDIF
                ENDDO

                IF ( n <= 16 ) THEN
                   d0r(i,j)  = 0.0 ! sr(n)
                ELSE
                ! n is index of masscum just above d0
                   tmp =  sr(n-1) + (0.5*qr - masscum(n-1))/(masscum(n) - masscum(n-1))*(sr(n) - sr(n-1)) ! radius
                  ! pltdp(i,j,k) = 2.*tmp*10. ! *2 for diam, *10 for units of mm
                   d0r(i,j) = tmp
                ENDIF
                
                d0rmax = Max( d0rmax, d0r(i,j) )

                
              ELSE
                 d0r(i,j) = 0.0
              ENDIF
              
!              IF ( (i/10)*10 == i .and. (j/10)*10 == j .and. qr /= 0.0d0 ) THEN
!                write(0,*) 'i,j,qr,Nr,d0r = ',i,j,qr,chw,d0r(i,j)
!                IF ( qr < 0.0 .and. chw < 0.0 ) THEN
!                  DO n = 1,nd_r_in
!                    write(0,*) 'n,nr = ',n,Ndr_in_arr(n,i,j),Ndr_in_arr(n,i,j)*sxdx(n)
!                  ENDDO
!                ENDIF
!              ENDIF

        ! Assign input variables
        
        rhoa = rhoa_arr(i,j,k)
        
        !write(*,*) 'i,j,k',i,j,k
        
        ! Compute the polarimetric variables
      
        
!      END DO
    END DO
  END DO
  
  write(0,*) 'qrmax,d0rmax = ',qrmax,d0rmax
      
  RETURN

END SUBROUTINE d0r_rsa_array_takahashi


#ifdef COMMAS
SUBROUTINE refl_rsa_array(addtorain,MPflg,MFflg,nx,ny,nz,rsafndir,wavelen,ibgn,iend,jbgn,jend,kbgn,kend,      &
           rhoa_arr,qr,qs,qg,qh,Ntr,Nts,Ntg,Nth,alphar,alphas,alphag,alphah,qsw,qgw,qhw,rhogrpl,rhohail,      &
           logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh,tair_arr,dpmask,dbz)
#else
SUBROUTINE refl_rsa_array(addtorain,MPflg,MFflg,nx,ny,nz,rsafndir,wavelen,ibgn,iend,jbgn,jend,kbgn,kend,      &
           rhoa_arr,qr,qs,qg,qh,Ntr,Nts,Ntg,Nth,alphar,alphas,alphag,alphah,qsw,qgw,qhw,rhogrpl,rhohail,      &
           logZ,sumZh,sumZv,logZdr,logZdrrain,sumZhv,Kdp,Ahh,Avv,rhv,rhoms,rhomg,rhomh,tair_arr)
#endif
!           logZ,sumZh,sumZv,logZdr,sumZhv,Kdp,Ahh,Avv,rhv,logZdrrain,rhoms,rhomg,rhomh,tair_arr)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This is a wrapper subroutine for the subroutine refl_rsa, to allow iteration
! over a 3D array, and to accomodate direct passing of microphysics variables
! through the argument list
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 2/10/2011
!
! MODIFIED: Dan Dawson and others, a bunch of times between 2011-2015
!
! MODIFIED: Dan Dawson, 05/13/15
!           General code cleanup.
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
#ifdef COMMAS
  USE dualpara, only: refl_rsa,lambda
#else
  USE dualpara
#endif
  IMPLICIT NONE

  INTEGER :: MPflg,MFflg
  LOGICAL :: addtorain
  INTEGER :: nx,ny,nz,ibgn,iend,jbgn,jend,kbgn(nx,ny),kend(nx,ny)
  REAL :: fws,fwg,fwh
  CHARACTER (LEN=256) :: rsafndir
  REAL, INTENT(IN) :: wavelen
  REAL, INTENT(IN) :: rhoa_arr(nx,ny,nz)
  REAL, INTENT(IN) :: tair_arr(nx,ny,nz)
  REAL, INTENT(INOUT) :: rhogrpl(nx,ny,nz),rhohail(nx,ny,nz)
  REAL :: qr(nx,ny,nz),qs(nx,ny,nz),qg(nx,ny,nz),qh(nx,ny,nz)
  REAL :: Ntr(nx,ny,nz),Nts(nx,ny,nz),Ntg(nx,ny,nz),Nth(nx,ny,nz)
  REAL :: alphar(nx,ny,nz),alphas(nx,ny,nz),alphag(nx,ny,nz),alphah(nx,ny,nz)
  REAL :: qsw(nx,ny,nz),qgw(nx,ny,nz),qhw(nx,ny,nz)
!f2py intent(in,out) qr,qs,qg,qh,Ntr,Nts,Ntg,Nth,alphar,alphas,alphag,alphah,qsw,qgw,qhw  
  REAL, INTENT(OUT) :: logZ(nx,ny,nz),sumZh(nx,ny,nz),sumZv(nx,ny,nz),logZdr(nx,ny,nz)
  real, intent(out) :: logZdrrain(nx,ny,nz)
  REAL, INTENT(OUT) :: sumZhv(nx,ny,nz),Kdp(nx,ny,nz),Ahh(nx,ny,nz),Avv(nx,ny,nz),rhv(nx,ny,nz)
  REAL, INTENT(INOUT) :: rhoms(nx,ny,nz),rhomg(nx,ny,nz),rhomh(nx,ny,nz)

  REAL :: qsout,qgout,qhout,qrout,nsout,ngout,nhout,nrout,rhomsout,rhomgout,rhomhout
  REAL*8 :: alfrout,alfsout,alfgout,alfhout

#ifdef COMMAS
  integer, intent(inout) :: dpmask(nx,ny,nz)
  real, intent(in) :: dbz(nx,ny,nz) ! model dBZ
  real :: rhog, rhoh
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw
#endif

  REAL :: rhoa,tair_C
  INTEGER :: i,j,k
  logical :: mask_logz
  !INTEGER :: bin_opt
  
!  write(0,*) 'start refl_rsa_array'
! Initialize radar wavelength

  lambda = wavelen

! DTD: Read RSA table here
  
  CALL set_nbins(0)

  CALL read_table(rsafndir)

  ! Assume that mu is constant for each category, and assign it here based on the value of MPflg

  IF(MPflg == 2) THEN ! We are using ZVD scheme with gamma-diameter for rain

    T_mur = 1.0d0/3.0d0
    T_mus = 1.0d0       ! Only poor snow is still gamma volume!
    T_mug = 1.0d0/3.0d0
    T_muh = 1.0d0/3.0d0
    
  ELSEIF(MPflg == 1) THEN ! We are using the ZVD scheme, which has an additional shape parameter mu for
                       ! the gamma distribution (setting it to 1/3 simplifies it to the standard gamma)
    T_mur = 1.0d0
    T_mus = 1.0d0
    T_mug = 1.0d0/3.0d0
    T_muh = 1.0d0/3.0d0

  ELSE                 ! Milbrandt and Yau scheme, LFO, etc. (based on standard gamma or exponential, no mu)

    T_mur = 1.0d0/3.0d0
    T_mus = 1.0d0/3.0d0
    T_mug = 1.0d0/3.0d0
    T_muh = 1.0d0/3.0d0

  END IF

  

  IF ( jbgn /= jend ) THEN ! horizontal x-section
#ifdef COMMAS
!$OMP PARALLEL DO DEFAULT(SHARED) schedule(dynamic, 1) &
!$OMP PRIVATE(T_log_ref, T_sum_ref_h, T_sum_ref_v, T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
!$OMP   T_Ahh, T_Avv, T_qr, T_qs, T_qh, T_qg, T_Ntr, T_Nts, T_Nth, T_Ntg, &
!$OMP   T_qsw, T_qgw, T_qhw, T_alfr,T_alfs,T_alfh,T_alfg, i, j, k, rhog, rhoh, rhoa, tair_C, &
!$OMP   qsout,qgout, qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout, &
!$OMP   rhomgout,rhomhout,fws,fwg,fwh )
#endif
  DO j=jbgn,jend
    DO i=ibgn,iend
      DO k=kbgn(i,j),kend(i,j)

#ifdef COMMAS
        IF ( dpmask(i,j,k) == 0 ) THEN
           cycle
        ELSE
            dpmask(i,j,k) = 0 ! done with this point
        ENDIF
#endif
        ! Assign input variables
        T_qr = qr(i,j,k)
        T_qs = qs(i,j,k) 
        T_qh = qh(i,j,k)
        T_qg = qg(i,j,k)
        T_Ntr = Ntr(i,j,k)
        T_Nts = Nts(i,j,k)
        T_Ntg = Ntg(i,j,k)
        T_Nth = Nth(i,j,k)
        T_alfr = alphar(i,j,k)
        T_alfs = alphas(i,j,k)
        T_alfg = alphag(i,j,k)
        T_alfh = alphah(i,j,k)

        IF(MFflg == 2) THEN ! Explicit melted water fraction
          T_qsw = qsw(i,j,k)
          T_qgw = qgw(i,j,k)
          T_qhw = qhw(i,j,k)
        ELSE ! No explicit melted water fraction and/or use fractionWater subroutine
          T_qsw = 0.0
          T_qgw = 0.0
          T_qhw = 0.0
        END IF

        rhog = rhogrpl(i,j,k)
        rhoh = rhohail(i,j,k)

        rhoa = rhoa_arr(i,j,k)
        tair_C = tair_arr(i,j,k)-273.15
        ! Compute the polarimetric variables
#ifdef COMMAS
        CALL refl_rsa(addtorain,i,j,k,MPflg,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,     &
                      qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout,rhomgout,rhomhout,tair_C, &
                     T_log_ref, T_sum_ref_h, T_sum_ref_v, &
                     T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
                     T_Ahh,     T_Avv,  &
                     T_qr, T_qs, T_qh, T_qg, &
                     T_Ntr, T_Nts, T_Nth, T_Ntg, &
                     T_alfr,T_alfs,T_alfh,T_alfg, &
                     T_mur, T_mus, T_mug, T_muh, &
                     T_qsw, T_qgw, T_qhw, rhog, rhoh )
#else
        CALL refl_rsa(addtorain,i,j,k,MPflg,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,     &
                      qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout,rhomgout,rhomhout,tair_C)
#endif
        IF(MFflg <= 4) THEN ! Compute water fraction mixing ratios from those calculated with fractionWater subroutine
          qsw(i,j,k) = fws*qsout
          qgw(i,j,k) = fwg*qgout
          qhw(i,j,k) = fwh*qhout
          
          IF ( tair_C < 4.0 .and. qs(i,j,k) > 1.e-4 ) THEN
!           write(0,*) 't,qs,cx,tref = ',tair_C,qs(i,j,k),qsout,nts(i,j,k),nsout,t_log_ref,rhomsout
          ENDIF
          
          ! Also adjust qr, qs, qg, qh, nr, ns, ng, and nh appropriately
          qr(i,j,k) = qrout
          qs(i,j,k) = qsout
          qg(i,j,k) = qgout
          qh(i,j,k) = qhout
          Ntr(i,j,k) = nrout
          Nts(i,j,k) = nsout
          Ntg(i,j,k) = ngout
          Nth(i,j,k) = nhout
          alphar(i,j,k) = alfrout
          alphas(i,j,k) = alfsout
          alphag(i,j,k) = alfgout
          alphah(i,j,k) = alfhout
          rhoms(i,j,k) = rhomsout
          rhomg(i,j,k) = rhomgout
          rhomh(i,j,k) = rhomhout
          
        END IF

! #ifdef COMMAS
!        IF ( i==15 ) THEN ! 121 .and. j==90 .and. k>=6 .and. k<=18 ) THEN
!          write(6,*) 'i,j,k,dbzmodel,dbzdp = ',i,j,k, dbz(i,j,k), T_log_ref
!        ENDIF
! #endif
        ! Store result in arrays for output
        logZ(i,j,k) = T_log_ref
        sumZh(i,j,k) = T_sum_ref_h
        sumZv(i,j,k) = T_sum_ref_v
        logZdr(i,j,k) = T_log_zdr
!        logZdr(i,j,k) = T_log_zdr_rain
        logZdrrain(i,j,k) = T_log_zdr_rain
        sumZhv(i,j,k) = T_sum_ref_hv
        Kdp(i,j,k) = T_kdp
        Ahh(i,j,k) = T_Ahh
        Avv(i,j,k) = T_Avv
        ! Compute rho_hv
        IF(T_sum_ref_h*T_sum_ref_v > 0.) THEN
          rhv(i,j,k) = T_sum_ref_hv/                      &
            SQRT(T_sum_ref_h*T_sum_ref_v)
        ELSE
          rhv(i,j,k) = 0.0
        ENDIF
        
      END DO
    END DO
  END DO

  ELSE ! vertical x-section

#ifdef COMMAS
!$OMP PARALLEL DO DEFAULT(SHARED) schedule(dynamic, 1) &
!$OMP PRIVATE(T_log_ref, T_sum_ref_h, T_sum_ref_v, T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
!$OMP   T_Ahh, T_Avv, T_qr, T_qs, T_qh, T_qg, T_Ntr, T_Nts, T_Nth, T_Ntg, &
!$OMP   T_qsw, T_qgw, T_qhw, T_alfr,T_alfs,T_alfh,T_alfg, i, j, k, rhog, rhoh, rhoa, tair_C, &
!$OMP   qsout,qgout, qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout, &
!$OMP   rhomgout,rhomhout,fws,fwg,fwh )
#endif
    DO i=ibgn,iend
     DO j=jbgn,jend
      DO k=kbgn(i,j),kend(i,j)

#ifdef COMMAS
        IF ( dpmask(i,j,k) == 0 ) THEN
           cycle
        ELSE
            dpmask(i,j,k) = 0 ! done with this point
        ENDIF
#endif
        ! Assign input variables
        T_qr = qr(i,j,k)
        T_qs = qs(i,j,k) 
        T_qh = qh(i,j,k)
        T_qg = qg(i,j,k)
        T_Ntr = Ntr(i,j,k)
        T_Nts = Nts(i,j,k)
        T_Ntg = Ntg(i,j,k)
        T_Nth = Nth(i,j,k)
        T_alfr = alphar(i,j,k)
        T_alfs = alphas(i,j,k)
        T_alfg = alphag(i,j,k)
        T_alfh = alphah(i,j,k)

        IF(MFflg == 2) THEN ! Explicit melted water fraction
          T_qsw = qsw(i,j,k)
          T_qgw = qgw(i,j,k)
          T_qhw = qhw(i,j,k)
        ELSE ! No explicit melted water fraction and/or use fractionWater subroutine
          T_qsw = 0.0
          T_qgw = 0.0
          T_qhw = 0.0
        END IF

        rhog = rhogrpl(i,j,k)
        rhoh = rhohail(i,j,k)

        rhoa = rhoa_arr(i,j,k)
        tair_C = tair_arr(i,j,k)-273.15
        ! Compute the polarimetric variables
#ifdef COMMAS
        CALL refl_rsa(addtorain,i,j,k,MPflg,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,     &
                      qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout,rhomgout,rhomhout,tair_C, &
                     T_log_ref, T_sum_ref_h, T_sum_ref_v, &
                     T_log_zdr, T_sum_ref_hv, T_kdp, T_log_zdr_rain, & 
                     T_Ahh,     T_Avv,  &
                     T_qr, T_qs, T_qh, T_qg, &
                     T_Ntr, T_Nts, T_Nth, T_Ntg, &
                     T_alfr,T_alfs,T_alfh,T_alfg, &
                     T_mur, T_mus, T_mug, T_muh, &
                     T_qsw, T_qgw, T_qhw, rhog, rhoh )
#else
        CALL refl_rsa(addtorain,i,j,k,MPflg,MFflg,rhoa,fws,fwg,fwh,qsout,qgout,     &
                      qhout,qrout,nsout,ngout,nhout,nrout,alfrout,alfsout,alfgout,alfhout,rhomsout,rhomgout,rhomhout,tair_C)
#endif
        IF(MFflg <= 4) THEN ! Compute water fraction mixing ratios from those calculated with fractionWater subroutine
          qsw(i,j,k) = fws*qsout
          qgw(i,j,k) = fwg*qgout
          qhw(i,j,k) = fwh*qhout
          
          ! Also adjust qr, qs, qg, qh, nr, ns, ng, and nh appropriately
          qr(i,j,k) = qrout
          qs(i,j,k) = qsout
          qg(i,j,k) = qgout
          qh(i,j,k) = qhout
          Ntr(i,j,k) = nrout
          Nts(i,j,k) = nsout
          Ntg(i,j,k) = ngout
          Nth(i,j,k) = nhout
          alphar(i,j,k) = alfrout
          alphas(i,j,k) = alfsout
          alphag(i,j,k) = alfgout
          alphah(i,j,k) = alfhout
          rhoms(i,j,k) = rhomsout
          rhomg(i,j,k) = rhomgout
          rhomh(i,j,k) = rhomhout
          
        END IF

#ifdef COMMAS
!        IF ( i==15 ) THEN ! 121 .and. j==90 .and. k>=6 .and. k<=18 ) THEN
!          write(6,*) 'i,j,k,dbzmodel,dbzdp = ',i,j,k, dbz(i,j,k), T_log_ref
!        ENDIF
#endif

        ! Store result in arrays for output
        logZ(i,j,k) = T_log_ref
        sumZh(i,j,k) = T_sum_ref_h
        sumZv(i,j,k) = T_sum_ref_v
        logZdr(i,j,k) = T_log_zdr
!        logZdr(i,j,k) = T_log_zdr_rain
        logZdrrain(i,j,k) = T_log_zdr_rain
        sumZhv(i,j,k) = T_sum_ref_hv
        Kdp(i,j,k) = T_kdp
        Ahh(i,j,k) = T_Ahh
        Avv(i,j,k) = T_Avv
        ! Compute rho_hv
        IF(T_sum_ref_h*T_sum_ref_v > 0.) THEN
          rhv(i,j,k) = T_sum_ref_hv/                      &
            SQRT(T_sum_ref_h*T_sum_ref_v)
        ELSE
          rhv(i,j,k) = 0.0
        ENDIF
        
      END DO
    END DO
  END DO

  ENDIF

  RETURN

END SUBROUTINE refl_rsa_array

SUBROUTINE solve_alpha(nx,ny,nz,rhoa,cx,q,Ntx,Z,alpha)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Changed alpha array to double precision
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  INTEGER :: nx,ny,nz
  REAL    :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL*8,INTENT(OUT)  :: alpha(nx,ny,nz)

  REAL*8  :: solveAlpha

  REAL :: cx

  INTEGER i,j,k

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN

          alpha(i,j,k) = solveAlpha(q(i,j,k),Ntx(i,j,k),Z(i,j,k),cx,rhoa(i,j,k))

        ELSE
          alpha(i,j,k) = 0.d0
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE solve_alpha

SUBROUTINE solve_alpha_iter(nx,ny,nz,rhoa,mu,q,Ntx,Z,rhox,alpha)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (06/09/2011)
!
!  Calculates the shape parameter (alpha) for a gamma distribution
!  of the form N(D) = N0*D^alpha*exp(-lambda*D^3mu), using an iterative
!  approach.  When mu=1/3 the standard gamma distribution is recovered.
!  When mu = 1.0, you get a gamma distribution in terms of volume.  Right
!  now the subroutine only allows a fixed mu=1/3 or 1.
!
!  MODIFICATION HISTORY:
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: pi = 3.141592   ! pi
  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL, INTENT(IN) :: mu
  REAL, INTENT(IN)    :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL, INTENT(IN)    :: rhox(nx,ny,nz)
  REAL, INTENT(OUT)  :: alpha(nx,ny,nz)

  REAL :: temp,t1,t2,t3,t4,t5
  REAL :: nu

  INTEGER :: i,j,k,iter

  REAL, PARAMETER :: rnumin = -0.8
  REAL, PARAMETER :: rnumax = 8.0
  REAL, PARAMETER :: alpmin = 0.0
  REAL, PARAMETER :: alpmax = 15.0

  alpha = 0.0

  IF( Int(mu+0.001) == 1) THEN
    alpha = rnumin
#ifdef COMMAS
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,temp,nu)
#endif
    DO k=1,nz-IJSTAG
      DO j=1,ny-IJSTAG
        DO i=1,nx-IJSTAG
!          IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN
          IF(q(i,j,k) > 1.e-12 .and. Ntx(i,j,k) > 1.e-19 .and. Z(i,j,k) > 1.e-23) THEN
            ! Compute mean volume
            temp = rhoa(i,j,k)*q(i,j,k)/(1000.*max(1.0e-9,Ntx(i,j,k)))
            ! Compute initial value of nu before iteration
            nu = 36.*(alpha(i,j,k)+2.0)*Ntx(i,j,k)*temp**2/(Z(i,j,k)*pi**2)-1.
            ! Perform the iteration on nu
            DO iter=1,10 ! Ten iterations should be enough, according to Ted
              IF(abs(nu - alpha(i,j,k)) .lt. 0.01) EXIT
              alpha(i,j,k) = max(rnumin, min(rnumax,nu))
              nu = 36.*(alpha(i,j,k)+2.0)*Ntx(i,j,k)*temp**2/(Z(i,j,k)*pi**2)-1.
              nu = max(rnumin,min(rnumax,nu))
            END DO
            ! Convert alpha to "standard gamma" alpha (i.e. alpha = 3*nu+2)
            alpha(i,j,k) = max(alpmin,min(alpmax,3.*alpha(i,j,k)+2.))
          ELSE
            alpha(i,j,k) = 0.0
          END IF   
        END DO
      END DO
    END DO
  ELSE IF (NInt(1./mu) == 3 ) THEN
    write(0,*) 'in solve_alpha_iter'
    alpha = alpmin
#ifdef COMMAS
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,t1,t2,t3,t4,t5,temp,nu)
#endif
    DO k=1,nz-IJSTAG
      DO j=1,ny-IJSTAG
        DO i=1,nx-IJSTAG
!          IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN
          IF(q(i,j,k) > 1.e-12 .and. Ntx(i,j,k) > 1.e-12 .and. Z(i,j,k) > 1.e-23 .and. rhoa(i,j,k) > 0.) THEN
            t1 = q(i,j,k)
            t2 = Ntx(i,j,k)
            t3 = rhox(i,j,k)
            t4 = rhoa(i,j,k)
            t5 = Z(i,j,k)
           ! temp = Z(i,j,k)*(pi/6.*rhox(i,j,k))**2*Ntx(i,j,k)/((rhoa(i,j,k)*q(i,j,k))**2)
            temp = t5*(pi/6.*t3)**2*t2/((t4*t1)**2)
            ! Compute initial value of nu before iteration (note, nu is actually alpha)
            nu = (6.+alpha(i,j,k))*(5.+alpha(i,j,k))*(4.+alpha(i,j,k))/  &
                 ((3.+alpha(i,j,k))*(2.+alpha(i,j,k))*temp) - 1.0
            DO iter=1,10
              IF(abs(nu - alpha(i,j,k)) .lt. 0.01) EXIT
              alpha(i,j,k) = max(alpmin,min(alpmax,nu))
              nu = (6.+alpha(i,j,k))*(5.+alpha(i,j,k))*(4.+alpha(i,j,k))/  &
                 ((3.+alpha(i,j,k))*(2.+alpha(i,j,k))*temp) - 1.0
              nu = max(alpmin,min(alpmax,nu))
            END DO
          ELSE
            alpha(i,j,k) = 0.0
          END IF
        END DO
      END DO
    END DO
  ELSE
    write(*,*) 'Sorry, only mu = 1/3 or mu = 1 accepted for now. Try again later!'
    RETURN
  END IF

  RETURN

END SUBROUTINE solve_alpha_iter

INTEGER FUNCTION get_qgh_opt (graupel_ON, hail_ON)

  INTEGER :: graupel_ON,hail_ON

  get_qgh_opt = 0

  IF(graupel_ON == 0 .and. hail_ON == 0) THEN
    get_qgh_opt = 1
  ELSE IF(graupel_ON == 0 .and. hail_ON == 1) THEN
    get_qgh_opt = 2
  ELSE IF(graupel_ON == 1 .and. hail_ON == 0) THEN
    get_qgh_opt = 3
  ELSE IF(graupel_ON == 1 .and. hail_ON == 1) THEN 
    get_qgh_opt = 4
  ENDIF 

END FUNCTION get_qgh_opt

FUNCTION gamma_r8(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gamma_r8
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
  gamma_r8=tmp+log(stp*ser/x)
  gamma_r8= exp(gamma_r8)

END FUNCTION gamma_r8

FUNCTION solveAlpha(Q,N,Z,Cx,rho)

 IMPLICIT NONE

! PASSING PARAMETERS:
  real, INTENT(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real*8 :: solveAlpha
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
    write(*,*) '*** STOPPED in MODULE ### solveAlpha *** '
    write(*,*) '*** : ',Q,N,Z,Cx*1.9099,rho
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

!Piecewise-polynomial approximation of g(a) to solve for a:
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

     solveAlpha= max(0.,min(a,alphaMax))

  ELSE

     solveAlpha= 0.

  ENDIF

END FUNCTION solveAlpha
