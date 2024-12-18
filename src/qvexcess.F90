      SUBROUTINE QVEXCESS(ngs,mgs,qwvp0,qv0,qcw1,pres,thetap0,theta0, &
     &    qvex,pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felv,ss1,pk,ngscnt)
      
!#####################################################################
!  Purpose: find the amount of vapor that can be condensed to liquid
!#####################################################################

      use micro_module, only : iqvsopt, esbolton
      use param_module, only : rd, cp, cap => rcp, poo => p00, rdorv => epsilon
      
      implicit none

      integer ngs,mgs,ngscnt
      
      real theta2temp
      
      real qvex
      
      integer nqsat
      real fqsat, cbw
      
      real ss1  ! 'target' supersaturation
!
!  input arrays
!
      real qv0(ngs), qcw1(ngs), pres(ngs), qwvp0(mgs)
      real thetap0(ngs), theta0(ngs)
      real fcqv1(ngs), felv(ngs), pi0(ngs)
      real pk(ngs)
      
      real tabqvs(nqsat)
!
! Local stuff
!
      
      integer itertd
      integer ltemq
      real gamss
      real theta(ngs), qvap(ngs), pqs(ngs), qcw(ngs), qwv(ngs)
      real qcwtmp(ngs), qss(ngs), qvs(ngs), qwvp(ngs)
      real dqcw(ngs), dqwv(ngs), dqvcnd(ngs)
      real temg(ngs), temcg(ngs), thetap(ngs)
      
      real tfr
      parameter ( tfr = 273.15 )
            
      real, external :: POLYSVP1
      
      real :: tmp

!
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!
!
!  set up temperature and vapor arrays
!
      pqs(mgs) = (380.0)/(pres(mgs))
      thetap(mgs) = thetap0(mgs)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      qwvp(mgs) = qwvp0(mgs)
      qvap(mgs) = max( (qwvp0(mgs) + qv0(mgs)), 0.0 )
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
!      temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
!
!
!
!  reset temporaries for cloud particles and vapor
!
      
      qwv(mgs) = max( 0.0, qvap(mgs) )
      qcw(mgs) = max( 0.0, qcw1(mgs) )
!
!
      qcwtmp(mgs) = qcw(mgs)
!      temgtmp = temg(mgs)
!      temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap
!      thsave(mgs) = thetap(mgs)
      temcg(mgs) = temg(mgs) - tfr
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      IF ( iqvsopt == 0 ) THEN
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
        qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF
      qss(mgs) = (0.01*ss1 + 1.0)*qvs(mgs)
!
!  iterate  adjustment
!
      do itertd = 1,2
!
!
!  calculate super-saturation
!
      dqcw(mgs) = 0.0
      dqwv(mgs) = ( qwv(mgs) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           !  subsaturated
        if( qcw(mgs) .gt. -dqwv(mgs) ) then  ! check if qc can make up all of the deficit
          dqcw(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                 !  otherwise make all qc available for evap
          dqcw(mgs) = -qcw(mgs)
          dqwv(mgs) = dqwv(mgs) + qcw(mgs)
        end if
!
        qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs)  )  ! add to perturbation vapor
!
        qcw(mgs) = qcw(mgs) + dqcw(mgs)

        thetap(mgs) = thetap(mgs) +  &
     &                felv(mgs)*dqcw(mgs)/(cp*pi0(mgs))

      end if  ! dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      IF ( dqwv(mgs) .ge. 0. ) THEN
!
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/  &
     &  ((temg(mgs)-cbw)**2))
!
!
      dqcw(mgs) = dqvcnd(mgs)
!
      thetap(mgs) = thetap(mgs) +  &
     &   (felv(mgs)*dqcw(mgs) )    &
     & / (pi0(mgs)*cp)
      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )
      qcw(mgs) = qcw(mgs) + dqcw(mgs)
!
      END IF !  dqwv(mgs) .ge. 0.

      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
!      temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
      qvap(mgs) = Max((qwvp(mgs) + qv0(mgs)), 0.0)
      temcg(mgs) = temg(mgs) - tfr
!      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      IF ( iqvsopt == 0 ) THEN
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      ELSEIF ( iqvsopt == 1 ) THEN
        qvs(mgs) = rdorv*esbolton*tabqvs(ltemq)/(pres(mgs) - esbolton*tabqvs(ltemq))
      ELSE
        tmp = min(0.99*pres(mgs),POLYSVP1(temg(mgs),0))
        qvs(mgs) = rdorv*tmp/(pres(mgs) - tmp)
      ENDIF
      qcw(mgs) = max( 0.0, qcw(mgs) )
      qwv(mgs) = max( 0.0, qvap(mgs))
      qss(mgs) = (0.01*ss1 + 1.0)*qvs(mgs)
      end do
!
!  end the saturation adjustment iteration loop
!
!
      qvex = Max(0.0, qcw(mgs) - qcw1(mgs) )

      RETURN
      END

!==========================================================================================!

 real function polysvp1(T,i_type)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      implicit none

      real    :: T
      integer :: i_type

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

!-------------------------------------------

      if (i_type.EQ.1 .and. T.lt.273.15) then
! ICE

! use Goff-Gratch for T < 195.8 K and Flatau et al. equal or above 195.8 K
         if (t.ge.195.8) then
            dt=t-273.15
            polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
            polysvp1 = polysvp1*100.
         else
            polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
                alog10(273.16/t)+0.876793*(1.-t/273.16)+ &
                alog10(6.1071))*100.
         end if

      elseif (i_type.EQ.0 .or. T.ge.273.15) then
! LIQUID

! use Goff-Gratch for T < 202.0 K and Flatau et al. equal or above 202.0 K
         if (t.ge.202.0) then
            dt = t-273.15
            polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
            polysvp1 = polysvp1*100.
         else
! note: uncertain below -70 C, but produces physical values (non-negative) unlike flatau
            polysvp1 = 10.**(-7.90298*(373.16/t-1.)+ &
                5.02808*alog10(373.16/t)- &
                1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+ &
                8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
                alog10(1013.246))*100.
         end if

         endif


 end function polysvp1

!------------------------------------------------------------------------------------------!
