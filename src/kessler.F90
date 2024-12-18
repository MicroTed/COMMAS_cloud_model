!------------------------------------------------------------------------------
!
!   /////////////////////          BEGIN         \\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE KESSLER   //////////////////////
!
! KESSLER does a 3 water catagory microphysical parameterization
!
! Created by:  Louis Wicker, July 18, 1988
! Latest Revision: 06-10-00
!
! Update Notes:
!------------------------------------------------------------------------------
      subroutine kessler(t,                 &    ! potential temperature
                         qv,                &    ! water vapor mixing ratio
                         qc,                &    ! cloud water mixing ratio at t 
                         qr,                &    ! rain  water mixing ratio at t 
                         precip,            &    ! 2D precip array 
                         tb, pb,            &    ! base state theta and pi
                         dzc,               &    ! 1 / vertical grid spacing
                         dt,                &    ! time step and grid size
                         qrprod, prod, vt, fq,  &
                         nx, ny, nz)            ! grid dimensions


      USE PARAM_MODULE
      USE COMMASMPI_MODULE
      USE TRAJ_MODULE

      implicit none

! Passed variables

      integer nx, ny, nz
      real     t    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qv    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qc    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qr    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qrprod(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    prod  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    vt    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    fq    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    precip(-ng+1:nx+ng,-ng+1:ny+ng,4)
      real    dzc   (-ng+1:nz+ng)
      real    tb    (-ng+1:nz+ng)
      real    pb    (-ng+1:nz+ng)
      real    dz
      real    dt

! Local variables

      integer i, j, k, nfall, n
      double precision ::  product, factor, f5, qvs, qv0, cr
      double precision ::  prod0
      double precision ::  ern
      double precision ::  rcgs(-ng+1:nz+ng)
      double precision ::  pinit(-ng+1:nz+ng)
      double precision ::  gam(-ng+1:nz+ng)
      double precision ::  temp, dtfall, mxfall
      double precision ::  rcgsi(-ng+1:nz+ng)
      double precision ::  vtden(-ng+1:nz+ng)

      real c1, c2, c3, c4
      parameter( c1 = .001, c2 = .001, c3 = 2.2, c4 = .875 )
      parameter( f5 = 237.3 * 17.27 * 2.5e6 / cp )
      parameter( mxfall = 10.0 )

!     real pow, a, b
!     pow(a,b) = exp( log(a) )

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer ixb, jyb, kzb
      integer ixe, jye, kze

!------------------------------------------------------------------------------
! Parameters for the time split terminal advection

      nfall  = 1 + int( mxfall * dt * dzc(1) )
      dtfall = dt / float(nfall)
      
!      IF ( nprock > 1 ) RETURN

!------------------------------------------------------------------------------
! Create constants

#ifdef MPI
      kzb = 0
      kze = ktile+1
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg
      
      do k = kzb,kze
#else
      DO k = 1,nz-1
#endif
       rcgs(k)   = 1.0e2*pb(k)**2.509/(287.04*tb(k))
       rcgsi(k)  = 1.0 / rcgs(k)
       pinit(k)  = 1.0e5*pb(k)**3.509
       gam(k)    = Lav/(Cp*pb(k))
       vtden(k)  = sqrt(0.0011225/rcgs(k))

      ENDDO

      if(kzend .eq. nzend) rcgs(nz) = 0.0
!------------------------------------------------------------------------------
! Outer JI LOOP

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,factor,temp,qvs,prod0,ern,product)
      do j = jyb,jye
       do k = kzb,kze
        do i = ixb,ixe
#else
! THIS OMP LINE HAS NOT BEEN TESTED!
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,factor,temp,qvs,prod0,ern,product)
      DO j = 1,ny-1

! Autoconversion + Collision/Collection
! Terminal velocity calculation and advection by 3rd order upwind Crowley

       DO k = 1,nz-1
        DO i = 1,nx-1
#endif

!        factor        = pow(1.+c3*dt*qr(i,j,k),-c4)
         factor        = 1.0 / (1.+c3*dt*qr(i,j,k)**c4)
         qrprod(i,j,k) = qc(i,j,k) * (1.0 - factor) + factor*c1*dt*max(qc(i,j,k)-c2,0.)
         prod(i,j,k)   = qr(i,j,k)
         vt(i,j,k)     = 36.34*(Max( 0.0,qr(i,j,k) )*rcgs(k))**0.1364 * vtden(k)
!         vt(i,j,k) = 0.0
!         vt(i,j,k)     = 36.34*pow(qr(i,j,k)*rcgs(k),0.1364) * vtden(k)

!         cr            = vt(i,j,k) * dt * dzc(k)

!         IF( k .GT. 1 .AND. k .LT. nz-1) THEN
!          fq(i,j,k) = -     (qr(i,j,k+1) - 5.*qr(i,j,k) - 2.*qr(i,j,k-1)) / 6.0    &
!                      -    cr*(            -    qr(i,j,k) +    qr(i,j,k-1)) / 2.0  &
!                      + cr*cr*(qr(i,j,k+1) - 2.*qr(i,j,k) +    qr(i,j,k-1)) / 6.0
!         ELSE
          fq(i,j,k)   = Max( 0.0, qr(i,j,k) )
!         ENDIF

!         fq(i,j,k) = max(min(fq(i,j,k), qr(i,j,k)), 0.0)

        ENDDO
       ENDDO

#ifdef MPI
   ixb = 1
   ixe = itile

   IF ( myprock == nprock ) THEN
   do i = ixb,ixe
     fq(i,j,nz) = 0.0
     vt(i,j,nz) = 0.0
   enddo
   ENDIF
   
#else
        fq(1:nx-1,j,nz) = 0.0
        vt(1:nx-1,j,nz) = 0.0
#endif
!------------------------------------------------------------------------------
! Fallout done with flux upstream
#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do k = kzb,kze
         factor = dt*dzc(k)*rcgsi(k)
         do i = ixb,ixe
#else
        DO k = 1,nz-1
         factor = dt*dzc(k)*rcgsi(k)
         DO i = 1,nx-1
#endif
          prod(i,j,k) = prod(i,j,k) - factor                 &
                      * (rcgs(k  )*fq(i,j,k  )*vt(i,j,k  )   &
                        -rcgs(k+1)*fq(i,j,k+1)*vt(i,j,k+1))

         ENDDO
        ENDDO

!------------------------------------------------------------------------------
! Compute instantaneous precip rate on grid
!    den of air    mix rat.   fall velo x  time step  x   1/den of water
!     g /cm^3   x    g/g    x    cm/s   x                   cm^3 / g    = precip rate (mm/s)
!     g /cm^3   x    g/g    x    cm/s   x     DT      x     cm^3 / g    = depth of precip (mm)

       IF ( myprock == 1 ) THEN
#ifdef MPI
        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do i = ixb,ixe
#else
        DO i = 1,nx-1
#endif
         precip(i,j,1) = rcgs(1)*prod(i,j,1)*1000.*vt(i,j,1)
         precip(i,j,3) = precip(i,j,3) + precip(i,j,1)*dt
        ENDDO
       ENDIF

!------------------------------------------------------------------------------
! Production of rain and deletion of qc
! Production of qc from supersaturation
! Evaporation of qr
#ifdef MPI
        kzb = 1
        kze = ktile+1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        do i = ixb,ixe
#else
       DO k = 1,nz-1
        DO i = 1,nx-1
#endif
         qc(i,j,k) = max(qc  (i,j,k) - qrprod(i,j,k),0.)
         qr(i,j,k) = max(prod(i,j,k) + qrprod(i,j,k),0.)

         temp  = pb(k) * t(i,j,k)
         qvs   = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pinit(k)
         qv0  = qv(i,j,k)

         prod0 = (qv0-qvs) / (1.+qvs*f5/(temp-36.)**2)

!        ern   = amin1(dt*(((1.6+124.9*pow(rcgs(k)*qr(i,j,k),.2046))
!    $         * pow(rcgs(k)*qr(i,j,k),.525))/(2.55e8/(pinit(k)*qvs)
!    $          +5.4e5))*(dim(qvs,qv(i,j,k))/(rcgs(k)*qvs)),
!    $           Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

!
! note: dim(x,y) = The difference X-Y if it is positive; otherwise zero.
         ern   = Min(dt*(((1.6+124.9*(rcgs(k)*qr(i,j,k))**.2046)  &
               *(rcgs(k)*qr(i,j,k))**.525)/(2.55e8/(pinit(k)*qvs)   &
                +5.4e5))*(ddim(qvs,qv0)/(rcgs(k)*qvs)),        &
                 Max(-prod0-qc(i,j,k),0.d0),qr(i,j,k))

! 
! Formula from Sun and Crook evap
!        ern = amin1(0.0486*dt*(qvs-qv(i,j,k))
!    $              *89.1*(rcgs(k)*qr(i,j,k))**0.65,
!    $         Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

! Next line shuts off conversion at low temperatures - keeps qc up at anvil level
! This was used in storm video

         IF(prod0 .lt. 0.0 .and. temp .lt. 253.) prod0 = 0.0

! Update all variables

         product = max(prod0,-qc(i,j,k))
         t (i,j,k) = t(i,j,k) + gam(k)*(product - ern)
         qv(i,j,k) = max(qv(i,j,k) - product + ern,0.)
         qc(i,j,k) = qc(i,j,k) + product
         qr(i,j,k) = qr(i,j,k) - ern

       IF ( allocated( microrates ) ) THEN
!         microrates(i,j,k,1) = 0.0
!         microrates(i,j,k,2) = 0.0
         microrates(i,j,k,3) = gam(k)*(product - ern)*3600. ! gam(k) = Lav/(Cp*pb(k))
!         microrates(i,j,k,4) = -qvevr(kz)*elv/(cp*piz(kz))*3600.
!         microrates(i,j,k,5) = -qrmlh(kz)*elf/(cp*piz(kz))*3600.
!         microrates(i,j,k,6) = pcond2(kz)
       ENDIF
         

        ENDDO
       ENDDO

      ENDDO ! j

      RETURN  
      END SUBROUTINE KESSLER




!------------------------------------------------------------------------------
!
!   /////////////////////          BEGIN         \\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE ZHANG   //////////////////////
!
! ZHANG does a 3 water catagory microphysical parameterization
!
! Created by:  Louis Wicker, Feb 1, 2008
!
! Update Notes:  APS: Not updated for MPI yet
!------------------------------------------------------------------------------
      subroutine zhang(t,                 &    ! potential temperature
                         qv,                &    ! water vapor mixing ratio
                         qc,                &    ! cloud water mixing ratio at t 
                         qr,                &    ! rain  water mixing ratio at t 
                         precip,            &    ! 2D precip array 
                         tb, pb,            &    ! base state theta and pi
                         dzc,               &    ! 1 / vertical grid spacing
                         dt,                &    ! time step and grid size
                         qrprod, prod, vt, fq,  &
                         nx, ny, nz)            ! grid dimensions

      
      USE PARAM_MODULE
      
      implicit none

!      include 'param.h'

! Passed variables

      integer nx, ny, nz
      real     t    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qv    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qc    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qr    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qrprod(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    prod  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    vt    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    fq (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    precip(-ng+1:nx+ng,-ng+1:ny+ng,4)
      real    dzc   (-ng+1:nz+ng)
      real    tb    (-ng+1:nz+ng)
      real    pb    (-ng+1:nz+ng)
      real    dz
      real    dt

! Local variables

      integer i, j, k, nfall, n
      real product, factor, f5, qvs, cr
      real prod0
      real ern
      real lw
      real rcgs(nz)
      real pinit(nz)
      real gam(nz)
      real temp, dtfall, mxfall
      real rcgsi(nz)
      real vtden(nz)
      real rinit(nz)
      real riniti(nz)
      
      real, allocatable :: Rc(:,:,:)
      real, allocatable :: Re(:,:,:)


      real c1, c2, c3, c4
      parameter( c1 = .001, c2 = .001, c3 = 2.2, c4 = .875 )
      parameter( f5 = 237.3 * 17.27 * 2.5e6 / cp )
      parameter( mxfall = 10.0 )

!     real pow, a, b
!     pow(a,b) = exp( log(a) )

!------------------------------------------------------------------------------
! Parameters for the time split terminal advection

      nfall  = 1 + ifix( mxfall * dt * dzc(1) )
      dtfall = dt / float(nfall)

!------------------------------------------------------------------------------
! Create constants

      DO k = 1,nz-1

       rcgs(k)   = 1.0e2*pb(k)**2.509/(287.04*tb(k))
       rcgsi(k)  = 1.0 / rcgs(k)
       rinit(k)  = 1.0e5*pb(k)**2.509/(287.04*tb(k))
       riniti(k) = 1.0 / rinit(k)
       pinit(k)  = 1.0e5*pb(k)**3.509
       gam(k)    = Lav/(Cp*pb(k))
       vtden(k)  = sqrt(0.0011225/rcgs(k))

      ENDDO

      rcgs(nz) = 0.0

      allocate( Rc(nx,ny,nz) )
      allocate( Re(nx,ny,nz) )

!------------------------------------------------------------------------------
! Outer JI LOOP

! THIS OMP LINE HAS NOT BEEN TESTED!
!!! !$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,factor,temp,qvs,prod0,ern,product,lw)
      DO j = 1,ny-1

! Autoconversion + Collision/Collection
! Terminal velocity calculation and advection by 3rd order upwind Crowley

       DO k = 1,nz-1
        DO i = 1,nx-1

! Kessler:
!         factor        = 1.0 / (1.+c3*dt*qr(i,j,k))**c4
!         qrprod(i,j,k) = qc(i,j,k) * (1.0 - factor) + factor*c1*dt*max(qc(i,j,k)-c2,0.)
!         prod(i,j,k)   = qr(i,j,k)
!         vt(i,j,k)     = 36.34*(qr(i,j,k)*rcgs(k))**0.1364 * vtden(k)
!         cr            = vt(i,j,k) * dt * dzc(k)
!          fq(i,j,k)   = qr(i,j,k)

! Zhang:
         IF ( qr(i,j,k) .gt. 1.e-12) THEN
!          lw = rinit(k)*qr(i,j,k)
          lw = 1.e3*rinit(k)*Max( 0.0, qr(i,j,k) )
         Re(i,j,k)= 4.63e-4 * lw**0.875
         Rc(i,j,k)= 4.92e-3 * lw**0.955
         Vt(i,j,k)= 5.499   * lw**0.0447 * vtden(k)

         ELSE
!          lw = 0.1
         Re(i,j,k)= 0.0
         Rc(i,j,k)= 0.0
         Vt(i,j,k)= 0.0
         ENDIF

         fq(i,j,k)   = qr(i,j,k)

!         factor        = 1.0 / (1.+dt*Rc(i,j,k))
!         qrprod(i,j,k) = qc(i,j,k) * (1.0 - factor) + factor*c1*dt*max(qc(i,j,k)-c2,0.)

         qrprod(i,j,k) = amin1( dt*(qc(i,j,k)*Rc(i,j,k) + c1*max(qc(i,j,k)-c2,0.0)), qc(i,j,k))
         prod(i,j,k)   = qr(i,j,k)


        ENDDO


       ENDDO

        prod(1:nx-1,j,nz) = 0.0
        fq(1:nx-1,j,nz) = 0.0
        vt(1:nx-1,j,nz) = 0.0

!------------------------------------------------------------------------------
! Fallout done with flux upstream

        DO k = 1,nz-1
         factor = dt*dzc(k)*rcgsi(k)
         DO i = 1,nx-1
          prod(i,j,k) = prod(i,j,k) - factor                 &
                      * (rcgs(k  )*fq(i,j,k  )*vt(i,j,k  )   &
                        -rcgs(k+1)*fq(i,j,k+1)*vt(i,j,k+1))
         ENDDO
        ENDDO

!------------------------------------------------------------------------------
! Compute instantaneous precip rate on grid
!    den of air    mix rat.   fall velo x  time step  x   1/den of water
!     g /cm^3   x    g/g    x    cm/s   x                   cm^3 / g    = precip rate (mm/s)
!     g /cm^3   x    g/g    x    cm/s   x     DT      x     cm^3 / g    = depth of precip (mm)

        DO i = 1,nx-1
         precip(i,j,1) = rcgs(1)*prod(i,j,1)*1000.*vt(i,j,1)
         precip(i,j,3) = precip(i,j,3) + precip(i,j,1)*dt
        ENDDO

!------------------------------------------------------------------------------
! Production of rain and deletion of qc
! Production of qc from supersaturation
! Evaporation of qr

       DO k = 1,nz-1
        DO i = 1,nx-1

         qc(i,j,k) = max(qc  (i,j,k) - qrprod(i,j,k),0.)
         qr(i,j,k) = max(prod(i,j,k) + qrprod(i,j,k),0.)

         temp  = pb(k) * t(i,j,k)
         qvs   = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pinit(k)

         prod0 = (qv(i,j,k)-qvs) / (1.+qvs*f5/(temp-36.)**2)

!        ern   = amin1(dt*(((1.6+124.9*pow(rcgs(k)*qr(i,j,k),.2046))
!    $         * pow(rcgs(k)*qr(i,j,k),.525))/(2.55e8/(pinit(k)*qvs)
!    $          +5.4e5))*(dim(qvs,qv(i,j,k))/(rcgs(k)*qvs)),
!    $           Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

! note: dim(x,y) = The difference X-Y if it is positive; otherwise zero.

!         ern   = amin1(dt*(((1.6+124.9*(rcgs(k)*qr(i,j,k))**.2046)  &
!               *(rcgs(k)*qr(i,j,k))**.525)/(2.55e8/(pinit(k)*qvs)   &
!                +5.4e5))*(dim(qvs,qv(i,j,k))/(rcgs(k)*qvs)),        &
!                 Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

! 
! Formula from Sun and Crook evap
!        ern = amin1(0.0486*dt*(qvs-qv(i,j,k))
!    $              *89.1*(rcgs(k)*qr(i,j,k))**0.65,
!    $         Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

! Formulas from Zhang et al.

!         ern = amin1( dt*Re(i,j,k), Max(-prod0-qc(i,j,k),0.), qr(i,j,k) )

!         ern   = amin1(dt*Re(i,j,k)*(dim(qvs,qv(i,j,k))/qvs),        &
!                 Max(-prod0-qc(i,j,k),0.),qr(i,j,k))

        ern = amin1(dt*(qvs-qv(i,j,k))*Re(i,j,k), &
               Max(-prod0-qc(i,j,k),0.),qr(i,j,k))


! Next line shuts off conversion at low temperatures - keeps qc up at anvil level
! This was used in storm video

         IF(prod0 .lt. 0.0 .and. temp .lt. 253.) prod0 = 0.0

! Update all variables

         product = max(prod0,-qc(i,j,k))
         t (i,j,k) = t(i,j,k) + gam(k)*(product - ern)
         qv(i,j,k) = max(qv(i,j,k) - product + ern,0.)
         qc(i,j,k) = qc(i,j,k) + product
         qr(i,j,k) = qr(i,j,k) - ern
         
        ENDDO
       ENDDO

      ENDDO
      
      deallocate( Rc )
      deallocate( Re )

      RETURN  
      END SUBROUTINE ZHANG
      
      
!------------------------------------------------------------------------------
!
!   /////////////////////          BEGIN         \\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE KESSLERQC  //////////////////////
!
! KESSLERQC does a 2 water catagory microphysical parameterization (no rain)
!
! Created by:  Louis Wicker, July 18, 1988
! Latest Revision: 06-10-00
!
!
! Update Notes:
!------------------------------------------------------------------------------
      subroutine kesslerqc(t,                 &    ! potential temperature
                         qv,                &    ! water vapor mixing ratio
                         qc,                &    ! cloud water mixing ratio at t 
                         precip,            &    ! 2D precip array 
                         tb, pb,            &    ! base state theta and pi
                         dzc,               &    ! 1 / vertical grid spacing
                         dt,                &    ! time step and grid size
                         qrprod, prod, vt, fq,  &
                         nx, ny, nz)            ! grid dimensions


      USE PARAM_MODULE
      USE COMMASMPI_MODULE

      implicit none

! Passed variables

      integer nx, ny, nz
      real     t    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qv    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qc    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qrprod(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    prod  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    vt    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    fq    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    precip(-ng+1:nx+ng,-ng+1:ny+ng,4)
      real    dzc   (-ng+1:nz+ng)
      real    tb    (-ng+1:nz+ng)
      real    pb    (-ng+1:nz+ng)
      real    dz
      real    dt

! Local variables

      integer i, j, k, nfall, n
      real product, factor, f5, qvs, cr
      real prod0
      real ern
      real rcgs(-ng+1:nz+ng)
      real pinit(-ng+1:nz+ng)
      real gam(-ng+1:nz+ng)
      real temp, dtfall, mxfall
      real rcgsi(-ng+1:nz+ng)
      real vtden(-ng+1:nz+ng)

      real c1, c2, c3, c4
      parameter( c1 = .001, c2 = .001, c3 = 2.2, c4 = .875 )
      parameter( f5 = 237.3 * 17.27 * 2.5e6 / cp )
      parameter( mxfall = 10.0 )

!     real pow, a, b
!     pow(a,b) = exp( log(a) )

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer ixb, jyb, kzb
      integer ixe, jye, kze
      
      logical :: debug_mpi = .true.

!------------------------------------------------------------------------------
! Parameters for the time split terminal advection

!------------------------------------------------------------------------------
! Create constants

!------------------------------------------------------------------------------
! Create constants

#ifdef MPI
      kzb = 0
      kze = ktile+1
      if(kzbeg .eq. nzbeg) kzb = 1
      if(kzend .eq. nzend) kze = kzend-kzbeg
      
      do k = kzb,kze
#else
      DO k = 1,nz-1
#endif
       rcgs(k)   = 1.0e2*pb(k)**2.509/(287.04*tb(k))
       rcgsi(k)  = 1.0 / rcgs(k)
       pinit(k)  = 1.0e5*pb(k)**3.509
       gam(k)    = Lav/(Cp*pb(k))
       vtden(k)  = sqrt(0.0011225/rcgs(k))

      ENDDO

      if(kzend .eq. nzend) rcgs(nz) = 0.0

!------------------------------------------------------------------------------
! Outer JI LOOP

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,factor,temp,qvs,prod0,ern,product)
      do j = jyb,jye
#else
! THIS OMP LINE HAS NOT BEEN TESTED!
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,factor,temp,qvs,prod0,ern,product)
      DO j = 1,ny-1

! Autoconversion + Collision/Collection
! Terminal velocity calculation and advection by 3rd order upwind Crowley

#endif


!------------------------------------------------------------------------------
! Production of rain and deletion of qc
! Production of qc from supersaturation
! Evaporation of qr
#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        do i = ixb,ixe
#else
       DO k = 1,nz-1
        DO i = 1,nx-1
#endif
         qc(i,j,k) = max(qc  (i,j,k) ,0.)
         temp  = pb(k) * t(i,j,k)
         qvs   = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pinit(k)

         prod0 = (qv(i,j,k)-qvs) / (1.+qvs*f5/(temp-36.)**2)


! Next line shuts off conversion at low temperatures - keeps qc up at anvil level
! This was used in storm video

         IF(prod0 .lt. 0.0 .and. temp .lt. 253.) prod0 = 0.0

! Update all variables

         product = max(prod0,-qc(i,j,k))
         t (i,j,k) = t(i,j,k) + gam(k)*(product)
         qv(i,j,k) = max(qv(i,j,k) - product,0.)
         qc(i,j,k) = qc(i,j,k) + product

        ENDDO
       ENDDO

      ENDDO

      RETURN  
      END SUBROUTINE KESSLERQC



