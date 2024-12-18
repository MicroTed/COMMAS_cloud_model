!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      subroutine takinitbin


        use takcommon,          only: pi,lmax,lfmax,kfmax,iimax,kimax          &
     &                              ,dji,sh0,sr0,sri0,dj,lr41,lr100,ls250,rhoi,shedsmall &
     &                              ,lmax1,lfmax1,dmu,accel,rhoa0,ceta,ck
         use comm1,           only: rhof
        use comm2,           only: sr,sx,sf,srf,fw,vw,vf,bew,cdw,szr,szf,szi      &
     &                            ,sxfmlr,isxfmlr,sxfinv,sxinv,sxfmlrinv             &
     &                            ,sri,sfi,di,dichg,srfchg,sxf,sxi,sh,sh00,ff
        use comm3,           only: workwp,workwm                               &
     &                            ,termkp,termkm,termfp,termfm,ref,rew,sqrtrefnummks
        use micro_module,    only: takgrmassopt,rho_qh_tak,rho_qhl_tak,takvf,icdx,icdxhl, &
     &                             cdhmin,cdhmax,cdhdnmin,cdhdnmax,cdhlmin,cdhlmax, &
     &                             cdhldnmin,cdhldnmax,takshedsmall,takshedsize1,takshedsize2,takshedsize3, &
     &                             takicethickopt,takicethickness,numshedregimes,mltdiam1, mltdiam2, mltdiam3 
        use takmpi,             only: rank,root
        USE COMMASMPI_MODULE,only: my_rank
        
        implicit none

       real, parameter :: rho00 = 1.225          ! reference/MSL air density
       real, parameter :: arx = 10.
       real, parameter :: frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.

          integer :: i,j,k,l,n,ii,kf,ki,kmm,jp,imt,jmt
          
          double precision :: ribar,sxibar,sim,total,dhj,bef,cdf
          real :: cd
          logical, parameter :: newice = .true.

!---


         RETURN
         END SUBROUTINE TAKINITBIN

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
