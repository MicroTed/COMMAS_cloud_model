      subroutine vegpbl(an,u,v,t,qv,p,piinit,  &
           tcanp, wcanp, qav, eflx, fflx, precip, veg,  &
           g12,glat,glon,  &
           zsfc,stype,tsfc,wsfc,tsoil,wsoil,vlai,  &
           uflx,vflx,tflx,qflx,radsw,radlw,albedo,rough,  &
           sfctke,sgz,wgz,mfc,mfe,time,nx,ny,nz,ns,dt)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      USE PARAM_MODULE
      USE SPHYS_MODULE
      USE RUN_ATT_NML
      USE COMMASMPI_MODULE

      implicit none
      
      integer nx,ny,nz,ns
      
      real  an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      real  u (-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  v (-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  t (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  p (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  piinit(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

      real precip(-ng+1:nx+ng,-ng+1:ny+ng)
      real sfctke(-ng+1:nx+ng,-ng+1:ny+ng)
      
      real glat(-ng+1:nx+ng,-ng+1:ny+ng)
      real glon(-ng+1:nx+ng,-ng+1:ny+ng)
      real zsfc(-ng+1:nx+ng,-ng+1:ny+ng)
      real g12 (-ng+1:nx+ng,-ng+1:ny+ng)

      real stype(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real wsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real tsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real tcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real wcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real qav  (-ng+1:nx+ng,-ng+1:ny+ng)
      real veg  (-ng+1:nx+ng,-ng+1:ny+ng)
      real vlai (-ng+1:nx+ng,-ng+1:ny+ng)

      real uflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real vflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real tflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real qflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real eflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real fflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real radsw (-ng+1:nx+ng,-ng+1:ny+ng)
      real radlw (-ng+1:nx+ng,-ng+1:ny+ng)
      real albedo(-ng+1:nx+ng,-ng+1:ny+ng)
      real rough (-ng+1:nx+ng,-ng+1:ny+ng)


      real sgz(nz), wgz(nz), mfc(nz), mfe(nz)

      real, dimension(:,:,:), allocatable :: thv

      real, dimension(:,:), allocatable :: cpm, rho
      real, dimension(:,:), allocatable :: qvssfc, thsfc, thvsfc
      real, dimension(:,:), allocatable :: qvscanp, qcanp, tav
      real, dimension(:,:), allocatable :: rpp, rars, wspd, vmag
      real, dimension(:,:), allocatable :: gz1z0, bri, mavail
      real, dimension(:,:), allocatable :: radlwdn
      real, dimension(:,:), allocatable :: cun, ctn, cu, ct, cf
      real, dimension(:,:), allocatable :: ustar, wstar, zpbl, molen

      real, dimension(13) :: wwlt, wsat

      real dt, mbar, nsum, npsum
      integer i,j,k
      integer time

      logical debug, pbldump
      real :: third = 1./3.
      real :: czo  = 0.032
      real :: ozo  = 1.E-4
      real :: kzo  = 1.0
      real :: c1   = 0.2721655
      real :: c2   = -0.33333
!.....msb 5/8/08  try critical soil moisture value of 0.40
!      real :: wcrit =  0.30= 
      real :: wcrit =  0.435
      real :: Ric   =  3.05      
      integer :: ivtype =   3              ! vegetation type (3 = grass)
      real :: ctv    =  0.005           ! vegetation heat transisivity
      real :: roughv =  0.05            ! roughness of vegetation (meters)

      real :: pi1, pres, tv, pib, pip
      real :: qsfc, dth, dum, vc, delz, roughn, azol
      real :: psi, fac, ra, solvar, rs, wcanpmax
      real :: delc, rsra, ta, bot, radswv, faclw, radlwg
      real :: radnet
      real :: fac1, fac2, fac3
      real :: dqvsdt, top, efpot, phis, cd, vel, tau
      real :: vtflx, vqflx, wtmp
      
      integer :: ixe,jye
!-----------------------------------------------------------------------

!
! --- Soil moisture saturation values by soil types
!
! 1 - sand
! 2 - loamy sand
! 3 - sandy loam
! 4 - silt loam
! 5 - silt
! 6 - loam
! 7 - clay
! 8 - clay loam
! 9 - 
!10 - 
!11 - 
!12 - 
!13 - water
!
      data wsat / .395, .410, .435, .485, .451, .420,  &
                 .477, .476, .426, .482, .482, 1.0E-20, 1.00/
!
! --- Vegetation wilting points by vegetation types
!
      data wwlt / .068, .075, .114, .179, .155, .175,  &
                 .218, .250, .219, .283, .286, 1.0E-20, 1.00/

!-----------------------------------------------------------------------
      debug = .false.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' Enter vegpbl ',nx,ny,nz

      allocate( thv(nz,nx,ny) )

      allocate(    cpm(nx,ny) ); allocate(   rho(nx,ny) )
      allocate( qvssfc(nx,ny) ); allocate( thsfc(nx,ny) )
      allocate( thvsfc(nx,ny) ); allocate( qvscanp(nx,ny) )
      allocate( qcanp(nx,ny) ); allocate(    tav(nx,ny) )
      allocate(   rpp(nx,ny) ); allocate(   rars(nx,ny) )
      allocate(  wspd(nx,ny) ); allocate(  gz1z0(nx,ny) )
      allocate(   bri(nx,ny) )
      allocate(mavail(nx,ny) )
      allocate(radlwdn(nx,ny) )
      allocate(  vmag(nx,ny) )

      allocate(  cun(nx,ny) ); allocate(  ctn(nx,ny) )

      allocate(   cu(nx,ny) ); allocate(   ct(nx,ny) )
      allocate(   cf(nx,ny) )
      allocate(ustar(nx,ny) ); allocate(wstar(nx,ny) )
      allocate( zpbl(nx,ny) ); allocate(molen(nx,ny) )


      ixe = itile
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend - jybeg

!     write(0,"('--------------VEGPBL: start --------------')")
!     do i=1,ixe
!     do j=1,ny
!     do k=1,2
!     IF( .not. ( u(i,j,k) .lt. 500. .and.  u(i,j,k) .gt. -500. ) ) then
!     print*, u(i,j,k), t(i,j,k), qv(i,j,k),i,j,k,time
!     end if
!     end do
!     end do
!     end do
     
     
      do i=1,ixe
      do j=1,jye
      uflx(i,j) = 0.
      vflx(i,j) = 0.
      eflx(i,j) = 0.
      fflx(i,j) = 0.
      tflx(i,j) = 0.
      qflx(i,j) = 0.
      end do
      end do


!     if(debug) write(*,*) ' vegpbl 1'
!$omp parallel do default(shared) &
!$omp private(i,j,k,pi1,pres,pib,pip,qsfc,dth,dum,vc, &
!$omp         delz,roughn,azol,psi,fac,tv)
      do k=1,nz-1
       do j=1,jye
        do i=1,ixe
          pi1 = (piinit(i,j,k) + p(i,j,k))
          pres = (psl*100.) * pi1**rcpinv
          thv(k,i,j) = t(i,j,k) * (1. + 0.61 * qv(i,j,k) )
        enddo
       enddo
      enddo

      do j=1,jye
      do i=1,ixe
        pi1 = (piinit(i,j,1) + p(i,j,1))
        tv = thv(1,i,j) * pi1
        pres = (psl*100.) * pi1**rcpinv
        rho(i,j) = pres / (rd * tv )
!        if(j.eq.1) print *,rho, pres, tv, rd
        cpm(i,j) = cp * (1. + 0.8 * qv(i,j,1) )
!
!--- Eliminate negative numbers from regridding
!
         wsfc(i,j) = MAX( wsfc(i,j),0.)
        wsoil(i,j) = MAX(wsoil(i,j),0.)
        wcanp(i,j) = MAX(wcanp(i,j),0.)
!     enddo
!
!----convert ground temperature to potential temperature:
!    set the ground mixing ratio from a moisture availability
!    determined from Deardorff 1978 (see Sun and Wu 1992)
!
!     if(debug) write(*,*) ' vegpbl 3'
!     do j=1,ny-1
        pib = .5 * (3.*piinit(i,j,1) - piinit(i,j,2))
        pip = .5 * (3.*p(i,j,1) - p(i,j,2))
        pi1 = (pib + pip)
        pres  = (psl*100.)* pi1**rcpinv
        thsfc(i,j)= tsfc(i,j) / pi1


!.......msb 5/7/10 - recalculate qvssfc to reduce latent heat fluxes
!       due to extreme ground temps on 5/22/02
!       try using just saturation at air temp
!       multiply by pi1 to get Temp from Pot Temp


!        qvssfc(i,j)= 380.*exp(17.27*(tsfc(i,j)-273.)/   &
!                       (tsfc(i,j) - 36.))/pres


        qvssfc(i,j)= 380.*exp(17.27*(t(i,j,1)*pi1 - 273.)/  &
                       (t(i,j,1)*pi1 - 36.))/pres
     

!        print*,qvssfc(i,j),j,i,t(i,j,1)*pi1,tsfc(i,j)
     
     
        qvscanp(i,j)= 380.*exp(17.27*(tcanp(i,j)-273.)/  &
                       (tcanp(i,j) - 36.))/pres

        mavail(i,j) = MIN(1., wsfc(i,j) / wcrit)

        qsfc = mavail(i,j) * qvssfc(i,j) +   &
                (1.-mavail(i,j)) * qv(i,j,1)
        qsfc = MIN( qvssfc(i,j), qsfc )
        thvsfc(i,j)= thsfc(i,j) * (1. + 0.61 * qsfc )
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get Bulk Ri '
!-----calculate bulk richardson no. of surface layer
!
!        vmag(i,j)=SQRT( u(i,j,1) * u(i,j,1) + v(i,j,1) * v(i,j,1) )
        vmag(i,j)=Max(5., SQRT( u(i,j,1) * u(i,j,1) + v(i,j,1) * v(i,j,1) ) )
        wspd(i,j) = ((1. - veg(i,j)) +   &
                    0.83 * veg(i,j) * SQRT(ctv)) * vmag(i,j)

!        print*,wspd(i,j),SQRT(ctv)

        dth = thv(1,i,j) - thvsfc(i,j)
        dum = - dth * ( 0.5 - sign(0.5, dth) )
        vc = 2. * SQRT(dum)

!        print*,u(i,j,1),v(i,j,1),wspd(i,j)

        wspd(i,j) = MAX( SQRT(wspd(i,j) * wspd(i,j) + vc * vc), 1.)

!        print*,wspd(i,j),vc

        delz = sgz(1) / g12(i,j)
        bri(i,j) = g / t(i,j,1) * (dth / delz) /   &
                  ( (wspd(i,j) * wspd(i,j)) / (delz * delz) )
        bri(i,j) = MIN(bri(i,j), 0.9 * Ric)

        roughn = (1.-veg(i,j))*rough(i,j) + veg(i,j) * roughv
        
        gz1z0(i,j) = ALOG( delz / roughn )

        azol = ALOG(0.025 * delz / roughn )
        cun(i,j) = 1. / ( (azol/        vonk ) + 8.4 )
        ctn(i,j) = 1. / ( (azol/(0.74 * vonk)) + 7.3 )
!-----------------------------------------------------------------------
!     Diagnose basic parameters for the appropriated stability class:
!     (Deardorff, 1972)
!-----------------------------------------------------------------------
        if (bri(i,j) .GE. 0. ) then                ! stable case
          cu(i,j) = cun(i,j) * (1. - bri(i,j) / Ric)
          ct(i,j) = ctn(i,j) * (1. - bri(i,j) / Ric)
        else                                     ! unstable case
          psi = alog10(-bri(i,j)) - 3.5

          fac = (1./cun(i,j)) - 25. * exp( psi * (0.26 - 0.03 * psi) )
          cu(i,j) = 1. / fac
          
          fac = (1./ctn(i,j)) + (1./cu(i,j)) - (1./cun(i,j))
          ct(i,j) = 1. / fac
        endif

      enddo
      enddo


!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get radiation  '
!-----compute radiation and soil thermal capacity:
!
        radsw(:,:) = 0.0
        radlw(:,:) = 0.0
        radlwdn(:,:) = 0.0
      IF ( isfcphys > 1 ) THEN
      
      call radflx(an,radsw,radlw,albedo,p,t,qv,piinit,stype,veg, &
                 tsfc,wsfc,tcanp,glat,glon,zsfc,wsat,radlwdn,  &
                 time,nx,ny,nz,ns)
                 
      ELSE
      ENDIF

!      print*,'after radflx'
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get vegetation stuff  '
! --- vegetation stuff
!


      do j=1,jye
      do i=1,ixe


!......msb 5/9/08  recalculate wind speed in foliage layer (uaf)

!        wspd(i,j)=SQRT( u(i,j,1) * u(i,j,1) + v(i,j,1) * v(i,j,1) )
        wspd(i,j) = ((1. - veg(i,j)) +  &
                   0.83 * veg(i,j) * SQRT(ctv)) * vmag(i,j)
                     


        ct(i,j) = (1. - veg(i,j)) * ct(i,j) + veg(i,j) * ctv/cu(i,j)
        cf(i,j) = 0.01 * ( 1. + 0.3 / (wspd(i,j)+1.E-10) )

        ra = 1. / (cf(i,j) * wspd(i,j) + 1.E-10 )
        fac = wwlt(ivtype) / ( 0.9 * wsoil(i,j) + 0.1 * wsfc(i,j) )
        solvar = (1. - albedo(i,j) ) * solarc /   &
                (radsw(i,j) + 0.03 * (1. - albedo(i,j)) * solarc)

!......msb 4/18/08 .....add in to make radiation at top of canpoy

!        solvar = (1. - albedov) * solarc /  &
!                (radsw(i,j)/((1-veg(i,j))*(1-albedo(i,j))) & 
!                + 0.03 * (1. - albedov)*solarc)

        rs = 200. * ( solvar +  0. + fac**2 )      ! s m**-1
        
        rars(i,j) = ra / (ra + rs)

        wcanpmax = MAX(.1 * veg(i,j), 1.E-15)
        delc = 0.5 + sign(0.5, qvscanp(i,j) - qav(i,j) )
        rsra = rs / (ra + rs)
        rpp(i,j) = 1. - rsra * delc *  &
                  (1. - ( wcanp(i,j) / wcanpmax)**(2.*third) )
        rpp(i,j) = MAX(1.E-10, rpp(i,j) )
        rpp(i,j) = MIN(1., rpp(i,j) )

!        print*,rpp(i,j),radsw(i,j),solvar,solarc,rs

        pi1 = (piinit(i,j,1) + p(i,j,1))
        ta = t(i,j,1) * pi1
        tav(i,j) = (1. - veg(i,j)) * ta + veg(i,j) * (  &
                   0.3 * ta + 0.6 * tcanp(i,j) + 0.1 * tsfc(i,j) )
     
        qcanp(i,j) = rpp(i,j) * qvscanp(i,j) + (1.-rpp(i,j)) * qav(i,j)
	    qcanp(i,j) = MIN(qcanp(i,j), qvscanp(i,j))

        bot = 1. / ( 1. - 0.1 * veg(i,j) * (1. - mavail(i,j)) )  
!        qav(i,j) = bot * ( (1. - 0.7 * veg(i,j)) * qv(i,j,1) +  & 
!                    veg(i,j) * (0.6 * qcanp(i,j) + &
!                    0.1 * mavail(i,j) * qvssfc(i,j)) )


!......msb 5/14/10...recalculate qav from eq 21b (Deardorff 1978)

        qav(i,j) = (1. - veg(i,j)) * qv(i,j,1) + veg(i,j) * (  &
                   0.3 * qv(i,j,1) + 0.6 * qcanp(i,j) + 0.1 * qvssfc(i,j) )

!         print*,qcanp(i,j),qav(i,j),j,i

      enddo
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get Tcanp  '
! --- diagnose the canapy temperature (eqn 32) 
!     Note that veg is not needed as it cancels from all terms!
!
!      do j=1,jye
      do i=1,ixe
!erm: 10/19/2010: restored original without veg., which keeps tcanp stable.
         radswv = radsw(i,j) * (1. - albedov) / (1. - albedo(i,j))
!         radswv = radsw(i,j)*(1. - albedov)/((1. - albedo(i,j))*(1.-veg(i,j)))         
 
!          print*,radsw(i,j),radswv,albedov,albedo(i,j),veg(i,j),j,i
         
         faclw = emissv * emissg / (emissv + emissg - emissv * emissg)
         radlwg = faclw * sbconst * tsfc(i,j)**4

         radnet = radswv + emissv * radlwdn(i,j) + radlwg

         fac3 = 7.7 * rho(i,j) * cpm(i,j) * cf(i,j) * wspd(i,j)
     
         fac1 = (emissv + 2.*emissg - emissv * emissg) /   &
               (emissv +    emissg - emissv * emissg)
         fac1 = fac1 * (emissv * sbconst)
         fac2 = 7. * rho(i,j) * cf(i,j) * wspd(i,j) *   &
                    Lav * rpp(i,j)

         dqvsdt = Lav * qvscanp(i,j) / (rv * tcanp(i,j)**2)

         top = ( radnet + fac3 * tav(i,j) +  &
                3. * fac1 * tcanp(i,j)**4 -  &
               fac2 * (qvscanp(i,j) - qav(i,j) - tcanp(i,j) * dqvsdt))

         bot = 4. * fac1 * tcanp(i,j)**3 + fac3 + fac2 * dqvsdt

         tcanp(i,j) = top / bot

!        print*,radsw(i,j),radswv,albedo(i,j),tcanp(i,j),t(i,j,1)*pi1,j,i

      enddo





!.......msb 4/8/08  diagnose canapy temperature (eq 32) (Tf)

!      do j=1,jye
    
!.....radsw is (1-albedog)(1-veg)*Sh_down in (eq 36)

!      radswv = radsw(i,j) * (1. - albedov)/((1. - albedo(i,j))*(1-veg(i,j)))  
        
!      faclw = emissv * emissg / (emissv + emissg - emissv * emissg)  
!      radlwg = faclw * sbconst * tsfc(i,j)**4        

!      radnet = radswv + emissv * radlwdn(i,j) + radlwg        
       
!      fac1 = (emissv + 2.*emissg - emissv * emissg) /   &
!            (emissv +    emissg - emissv * emissg)
      
!      fac1 = fac1 * (emissv * sbconst)

      
!      enddo








!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get fluxes  '
!-----compute radiation and soil thermal capacity:
!-----compute the surface sensible and latent heat fluxes
!     NOTE: over water, alter roughness length (zl, or rough)
!
!      do j=1,jye
      do i=1,ixe
!
!-----compute surface moist flux:
!

!......msb 04/21/08 add in (1-veg) (Chen + Dudhia 2001; WRF model)
!......msb 5/9/08 add in dirunal dependence of sfc evaporation.
!......msb 5/7/10 take out (1-veg), change calulation of qvssfc

!        qflx(i,j) = rho(i,j) * ct(i,j)*cu(i,j)*mavail(i,j)*wspd(i,j)

        qflx(i,j) = rho(i,j) * ct(i,j)*cu(i,j)*mavail(i,j)*  &
       wspd(i,j)*(1.0-veg(i,j))  &
       *sin((2.0*pii/24.0)*((hour+minute/60.0+time/3600.0)-13.5))

!        qflx(i,j) = rho(i,j) * ct(i,j)*cu(i,j)*mavail(i,j)*wspd(i,j) &
!       *sin((2*3.1415/24.0)*((hour+minute/60.0+time/3600.0)-13.5))
     
        qflx(i,j) = MAX( qflx(i,j) * (qvssfc(i,j) - qav(i,j)), 0.)

!-----compute evaporation moist flux (eqn 25c) :

        fflx(i,j) = (7.*veg(i,j)) * rho(i,j) * cf(i,j) * wspd(i,j) * &
                   (qvscanp(i,j) - qav(i,j)) * rpp(i,j)
!
!-----compute vegetation moist flux (eqn 26) :
!
        efpot = fflx(i,j) / rpp(i,j)
        delc = 0.5 + sign(0.5, qvscanp(i,j) - qav(i,j) )

        eflx(i,j) = delc * efpot * rars(i,j) *  &
                   (1. - ( wcanp(i,j) / wcanpmax)**(2.*third) )
                   
!        print*,sin((2*3.1415/24.0)*((hour+minute/60.0)-13.5))
!        print*,qflx(i,j)
             
!         print*,qflx(i,j),qvssfc(i,j),qav(i,j),ct(i,j)*cu(i,j),  &
!        mavail(i,j),wspd(i,j),rho(i,j),(1.0-veg(i,j))

!
!-----compute surface sensible heat flux with a lower cooling limit
!
        pi1 = (piinit(i,j,1) + p(i,j,1))
!        tflx(i,j) = cpm(i,j) * rho(i,j) * ct(i,j) * cu(i,j) *  &  ! hack
!                   4.0*wspd(i,j) * (tsfc(i,j) - tav(i,j))             ! hack
        tflx(i,j) = cpm(i,j) * rho(i,j) * ct(i,j) * cu(i,j) *  &
                   wspd(i,j) * (tsfc(i,j) - tav(i,j))
!       tflx(i,j) = MAX( tflx(i,j), -100.0 )

!        print*,tflx(i,j),ct(i,j),cu(i,j),tsfc(i,j),tav(i,j),wspd(i,j),rho(i,j),cpm(i,j),j,i
!
!-----compute surface momentum fluxes
!

!        print*,fflx(i,j),qflx(i,j),qvssfc(i,j),qvscanp(i,j),qav(i,j),wspd(i,j), &
!       cu(i,j),ct(i,j),mavail(i,j),rpp(i,j),tcanp(i,j),tsfc(i,j),eflx(i,j) 

        phis = g * zsfc(i,j)
        cd = cu(i,j)**2 + 3.0E-03 * ( phis / (phis + 9800.) )
        vel = vmag(i,j) + 1.e-9 ! SQRT( u(i,j,1) * u(i,j,1) + v(i,j,1) * v(i,j,1) )+1.E-9
        tau = rho(i,j) * cd * vel**2

        uflx(i,j) = -tau * u(i,j,1) / vel
        vflx(i,j) = -tau * v(i,j,1) / vel
        ustar(i,j) = SQRT( MAX(tau / rho(i,j), 0.) )


!        print*,phis,zsfc(i,j),cd,uflx(i,j),vflx(i,j),j,i,vel**2,rho(i,j)

      enddo ! i
      enddo ! j
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get soil ebm  '
!
! --- given the surface fluxes compute the ground temperature
!
      call vegebm(tsfc,eflx,fflx,tflx,qflx,tsoil,wsfc,wsoil,tcanp, &
                   wcanp,precip,radsw,radlw,veg,nx,ny,dt)

!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get tendencies  '
!
! --- total surface layer sensible heat and moisture fluxes
!
      do j=1,jye
      do i=1,ixe
!         delz1 = sgz(1) / g12(i,j)
!         delz = g12(i,j) / (sgz(2)-sgz(1))

          pi1 = (piinit(i,j,1) + p(i,j,1))
          vtflx = 1.1 * vlai(i,j)*rho(i,j)*cpm(i,j)*cf(i,j)*wspd(i,j)* &
                (tcanp(i,j) - tav(i,j))
          tflx(i,j) = (tflx(i,j) + vtflx) / pi1     ! THETA not T

!          print*,vlai(i,j),cf(i,j),tcanp(i,j),tav(i,j),j,i

          vqflx = vlai(i,j)*rho(i,j)*cf(i,j)*wspd(i,j)*rpp(i,j)* &
                (qvscanp(i,j) - qav(i,j))
          qflx(i,j) = qflx(i,j) + vqflx


!          print*,vtflx,tflx(i,j),vqflx,qflx(i,j),qvscanp(i,j), &
!         qav(i,j),mavail(i,j),tav(i,j),tcanp(i,j),tsfc(i,j)

      enddo
!-----------------------------------------------------------------------
!     Diagnose the sfc TKE for mixing (from Troen and Marht?)
!-----------------------------------------------------------------------
!     if(debug) write(*,*) ' get sfctke ' 
!      do j=1,jye
      do i=1,ixe
! --- PBL depth from theta-v
        zpbl(i,j) = wgz(2)
        do k=2,nz-3
          if( thv(k,i,j) .LT. thv(1,i,j)) zpbl(i,j) = wgz(k+1)
!         if( i .EQ. 400 .AND. j .EQ. 2 ) then
!      write(*,*)k,zpbl(i,j),thv(1,i,j),thv(k,i,j),wgz(k+1),thvsfc(i,j)
!         endif
        enddo

         wtmp = MAX( g * zpbl(i,j) * tflx(i,j) /   &
                    (rho(i,j) * cpm(i,j) * thvsfc(i,j)), 0.)
         wstar(i,j) = wtmp**third

       molen(i,j)=-1.*(ustar(i,j)**3)*rho(i,j)*cpm(i,j)*thvsfc(i,j)/ &
                     (vonk * g * tflx(i,j))
 
! --- convective velocity
        cd = cu(i,j) * cu(i,j)
        if (bri(i,j) .GE. 0. ) then                ! stable case
          sfctke(i,j) = 3.75 * cd * wspd(i,j)**2
        else                                     ! unstable case
          sfctke(i,j) = (3.75 + (-bri(i,j))**.66667) * cd*wspd(i,j)**2 &
                      + wstar(i,j)**2
        endif
        sfctke(i,j) = MIN(MAX(sfctke(i,j),1.E-8),9.)

      enddo ! i
       IF ( myprocj == 1 ) sfctke(i,1) = sfctke(i,2)
       IF ( myprocj == nprocj )  sfctke(i,ny-1) = sfctke(i,ny-2)
      enddo ! j

!     if(debug) write(*,*) ' VEGPBL: deallocate  '
      pbldump = .false.
      if(pbldump) then
        write(55,'(2i6)') nx,ny
        do j=1,jye
        do i=1,ixe
          write(55,'(2i6,4f20.9)') i,j,rho(i,j),wstar(i,j), &
                    zpbl(i,j),molen(i,j)
        enddo
        enddo
        stop
      endif
      
!-----------------------------------------------------------------------
      deallocate( thv )

      deallocate(    cpm ); deallocate(   rho )
      deallocate( qvssfc ); deallocate( thsfc )
      deallocate( thvsfc ); deallocate( qvscanp )
      deallocate(  qcanp ); deallocate(    tav )
      deallocate(    rpp ); deallocate(   rars )
      deallocate(   wspd ); deallocate(  gz1z0 )
      deallocate(    bri )
      deallocate( mavail )
      deallocate( radlwdn)

      deallocate( cun ); deallocate( ctn ); deallocate( cu )
      deallocate( ct ); deallocate( cf )
      deallocate( ustar ); deallocate( wstar )
      deallocate( zpbl ) ; deallocate( molen )
!     if(debug) write(*,*) ' VEGPBL: exit  '
!-----------------------------------------------------------------------
      return
      end

