!-----------------------------------------------------------------------
!
!   /////////////////////            BEGIN        \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE RADFLX     ////////////////////
!
! Compute the surface short wave radiation
!-----------------------------------------------------------------------
      SUBROUTINE RADFLX(an,radsw,radlw,albedo,p,t,qv,piinit,stype,veg,  &
                tsfc,wsfc,tcanp,glat,glon,zsfc,wsat,radlwdn,    &
                time, nx,ny,nz,ns)


      USE PARAM_MODULE
      USE SPHYS_MODULE
      USE RUN_ATT_NML
      USE COMMASMPI_MODULE
      USE INDEX_MODULE
      
      implicit none
      
      integer, parameter :: npw=13, npl=46

      integer nx,ny,nz,ns

      real an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      real radsw (-ng+1:nx+ng,-ng+1:ny+ng)
      real radlw (-ng+1:nx+ng,-ng+1:ny+ng)
      real albedo(-ng+1:nx+ng,-ng+1:ny+ng)
      
      real  t (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  p (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)      
      
      real  piinit(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 


      real stype(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real wsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real tcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real veg  (-ng+1:nx+ng,-ng+1:ny+ng)


      real glat(-ng+1:nx+ng,-ng+1:ny+ng)
      real glon(-ng+1:nx+ng,-ng+1:ny+ng)
      real zsfc(-ng+1:nx+ng,-ng+1:ny+ng)

      
      real radlwdn(nx,ny)

      real wsat(13)
      real tac(3), tsc(3), clw(3)
      real tran(npw, npl)
      
      integer daymonth(12)
      integer time

      real, dimension(:,:,:), allocatable :: cldfr
      real :: rh

      real, dimension(:,:), allocatable :: radlwup
      real, dimension(:,:), allocatable :: zenith, thour, cjdayloc
      real, dimension(:,:), allocatable :: prw, psfc
      real, dimension(:,:), allocatable :: ftabs , ftsca , fbsc
      real, dimension(:,:), allocatable :: ftabsd, ftscad, fbscd
      real, dimension(:,:), allocatable :: path, tacc, tscc

      logical debug
      logical pbl_run
      
      integer i,j,k
      
      real    :: pi1, pi2, pi3, pres, temp, qvs
      real    :: cld1, cld2, cld3
      real    :: p1, p2, p3
      real    :: fac, t1
      real    :: prwpt, epsair
      real    :: depcl
      real    :: emissa
      real    :: sunup
      real    :: radmax
      
      integer :: m, im1, ip1, jm1, jp1
      
      real    :: dzdx, dzdy
      real    :: slpmag, slpdir
      real    :: term1, anncyc, adjust, etau, anglehr, soldec
      real    :: coszen, sinpsi, cospsi, azimuth, dipsi, cosi
      real    :: temp1, albz, temp2
      
      integer :: iradloc, jradloc
      
      real    :: bdbar, xser, hi

      integer :: ixe,jye

 DATA tran(:, 1) /  0.926, 0.868, 0.855, 0.846, 0.838, 0.832, 0.827, 0.822, 0.818, 0.814, 0.811, 0.776, 0.314/
 DATA tran(:, 2) /  0.915, 0.854, 0.840, 0.831, 0.823, 0.817, 0.811, 0.806, 0.802, 0.798, 0.794, 0.740, 0.313/
 DATA tran(:, 3) /  0.903, 0.841, 0.826, 0.816, 0.808, 0.802, 0.796, 0.791, 0.787, 0.782, 0.779, 0.707, 0.312/
 DATA tran(:, 4) /  0.892, 0.828, 0.813, 0.803, 0.795, 0.788, 0.782, 0.777, 0.772, 0.768, 0.764, 0.676, 0.311/
 DATA tran(:, 5) /  0.881, 0.815, 0.800, 0.790, 0.782, 0.775, 0.769, 0.763, 0.758, 0.754, 0.750, 0.648, 0.310/
 DATA tran(:, 6) /  0.870, 0.803, 0.788, 0.777, 0.769, 0.762, 0.756, 0.750, 0.745, 0.741, 0.737, 0.620, 0.310/
 DATA tran(:, 7) /  0.860, 0.792, 0.776, 0.765, 0.757, 0.750, 0.743, 0.738, 0.733, 0.728, 0.724, 0.595, 0.309/
 DATA tran(:, 8) /  0.850, 0.781, 0.765, 0.754, 0.745, 0.738, 0.731, 0.726, 0.721, 0.716, 0.712, 0.571, 0.308/
 DATA tran(:, 9) /  0.839, 0.770, 0.753, 0.742, 0.733, 0.726, 0.720, 0.714, 0.709, 0.704, 0.700, 0.549, 0.308/
 DATA tran(:,10) /  0.830, 0.759, 0.743, 0.731, 0.722, 0.715, 0.709, 0.703, 0.698, 0.693, 0.689, 0.527, 0.307/
 DATA tran(:,11) /  0.820, 0.748, 0.732, 0.721, 0.712, 0.704, 0.698, 0.692, 0.687, 0.682, 0.678, 0.507, 0.307/
 DATA tran(:,12) /  0.810, 0.738, 0.722, 0.710, 0.701, 0.694, 0.687, 0.681, 0.676, 0.671, 0.667, 0.488, 0.307/
 DATA tran(:,13) /  0.801, 0.728, 0.712, 0.700, 0.691, 0.683, 0.677, 0.671, 0.666, 0.661, 0.656, 0.470, 0.306/
 DATA tran(:,14) /  0.791, 0.719, 0.702, 0.690, 0.681, 0.674, 0.667, 0.661, 0.656, 0.651, 0.646, 0.453, 0.306/
 DATA tran(:,15) /  0.782, 0.709, 0.692, 0.681, 0.671, 0.664, 0.657, 0.651, 0.646, 0.641, 0.636, 0.437, 0.306/
 DATA tran(:,16) /  0.773, 0.700, 0.683, 0.671, 0.662, 0.654, 0.648, 0.642, 0.636, 0.631, 0.627, 0.422, 0.305/
 DATA tran(:,17) /  0.764, 0.691, 0.674, 0.662, 0.653, 0.645, 0.638, 0.632, 0.627, 0.622, 0.618, 0.407, 0.305/
 DATA tran(:,18) /  0.756, 0.682, 0.665, 0.653, 0.644, 0.636, 0.629, 0.623, 0.618, 0.613, 0.608, 0.393, 0.305/
 DATA tran(:,19) /  0.747, 0.673, 0.656, 0.644, 0.635, 0.627, 0.621, 0.615, 0.609, 0.604, 0.600, 0.380, 0.305/
 DATA tran(:,20) /  0.738, 0.665, 0.647, 0.636, 0.626, 0.619, 0.612, 0.606, 0.600, 0.596, 0.591, 0.367, 0.305/
 DATA tran(:,21) /  0.730, 0.656, 0.639, 0.627, 0.618, 0.610, 0.603, 0.597, 0.592, 0.587, 0.582, 0.355, 0.305/
 DATA tran(:,22) /  0.722, 0.648, 0.631, 0.619, 0.610, 0.602, 0.595, 0.589, 0.584, 0.579, 0.574, 0.344, 0.304/
 DATA tran(:,23) /  0.714, 0.640, 0.623, 0.611, 0.602, 0.594, 0.587, 0.581, 0.576, 0.571, 0.566, 0.333, 0.304/
 DATA tran(:,24) /  0.706, 0.632, 0.615, 0.603, 0.594, 0.586, 0.579, 0.573, 0.568, 0.563, 0.558, 0.322, 0.304/
 DATA tran(:,25) /  0.698, 0.624, 0.607, 0.595, 0.586, 0.578, 0.571, 0.565, 0.560, 0.555, 0.550, 0.312, 0.304/
 DATA tran(:,26) /  0.690, 0.616, 0.599, 0.588, 0.578, 0.571, 0.564, 0.558, 0.552, 0.547, 0.543, 0.302, 0.304/
 DATA tran(:,27) /  0.683, 0.609, 0.592, 0.580, 0.571, 0.563, 0.556, 0.550, 0.545, 0.540, 0.535, 0.293, 0.304/
 DATA tran(:,28) /  0.675, 0.602, 0.585, 0.573, 0.564, 0.556, 0.549, 0.543, 0.538, 0.533, 0.528, 0.284, 0.304/
 DATA tran(:,29) /  0.668, 0.594, 0.577, 0.566, 0.556, 0.549, 0.542, 0.536, 0.531, 0.526, 0.521, 0.276, 0.304/
 DATA tran(:,30) /  0.661, 0.587, 0.570, 0.559, 0.549, 0.542, 0.535, 0.529, 0.524, 0.519, 0.514, 0.268, 0.304/
 DATA tran(:,31) /  0.653, 0.580, 0.563, 0.552, 0.542, 0.535, 0.528, 0.522, 0.517, 0.512, 0.507, 0.260, 0.304/
 DATA tran(:,32) /  0.646, 0.573, 0.556, 0.545, 0.536, 0.528, 0.521, 0.515, 0.510, 0.505, 0.501, 0.253, 0.304/
 DATA tran(:,33) /  0.639, 0.567, 0.550, 0.538, 0.529, 0.521, 0.515, 0.509, 0.503, 0.499, 0.494, 0.245, 0.304/
 DATA tran(:,34) /  0.633, 0.560, 0.543, 0.532, 0.522, 0.515, 0.508, 0.502, 0.497, 0.492, 0.488, 0.238, 0.304/
 DATA tran(:,35) /  0.626, 0.553, 0.537, 0.525, 0.516, 0.508, 0.502, 0.496, 0.491, 0.486, 0.481, 0.232, 0.304/
 DATA tran(:,36) /  0.619, 0.547, 0.530, 0.519, 0.510, 0.502, 0.496, 0.490, 0.484, 0.480, 0.475, 0.225, 0.304/
 DATA tran(:,37) /  0.613, 0.541, 0.524, 0.512, 0.503, 0.496, 0.489, 0.484, 0.478, 0.473, 0.469, 0.219, 0.304/
 DATA tran(:,38) /  0.606, 0.534, 0.518, 0.506, 0.497, 0.490, 0.483, 0.477, 0.472, 0.467, 0.463, 0.213, 0.304/
 DATA tran(:,39) /  0.600, 0.528, 0.512, 0.500, 0.491, 0.484, 0.477, 0.472, 0.466, 0.462, 0.457, 0.208, 0.304/
 DATA tran(:,40) /  0.594, 0.522, 0.506, 0.494, 0.486, 0.478, 0.472, 0.466, 0.461, 0.456, 0.451, 0.202, 0.304/
 DATA tran(:,41) /  0.587, 0.516, 0.500, 0.489, 0.480, 0.472, 0.466, 0.460, 0.455, 0.450, 0.446, 0.197, 0.304/
 DATA tran(:,42) /  0.581, 0.511, 0.494, 0.483, 0.474, 0.467, 0.460, 0.454, 0.449, 0.444, 0.440, 0.192, 0.304/
 DATA tran(:,43) /  0.575, 0.505, 0.488, 0.477, 0.468, 0.461, 0.455, 0.449, 0.444, 0.439, 0.435, 0.187, 0.304/
 DATA tran(:,44) /  0.569, 0.499, 0.483, 0.472, 0.463, 0.456, 0.449, 0.444, 0.438, 0.434, 0.429, 0.182, 0.304/
 DATA tran(:,45) /  0.563, 0.494, 0.477, 0.466, 0.458, 0.450, 0.444, 0.438, 0.433, 0.428, 0.424, 0.178, 0.304/
 DATA tran(:,46) /  0.558, 0.488, 0.472, 0.461, 0.452, 0.445, 0.439, 0.433, 0.428, 0.423, 0.419, 0.173, 0.304/
      
!-----------------------------------------------------------------------

      ixe = itile
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend - jybeg

      allocate(cldfr(3 ,nx,ny) )

      allocate(  radlwup(nx,ny) )
      allocate(   zenith(nx,ny) ); allocate( thour(nx,ny) )
      allocate( cjdayloc(nx,ny) )
      allocate(      prw(nx,ny) ); allocate(  psfc(nx,ny) )
      allocate(    ftabs(nx,ny) ); allocate( ftsca(nx,ny) )
      allocate(     fbsc(nx,ny) )
      allocate(   ftabsd(nx,ny) ); allocate( ftscad(nx,ny) )
      allocate(    fbscd(nx,ny) )
      allocate(     path(nx,ny) ); allocate( tacc(nx,ny) )
      allocate(     tscc(nx,ny) )
!-----------------------------------------------------------------------
      debug = .false.
      pbl_run = .false.
!-----------------------------------------------------------------------
!
! --- Soil moisture saturation values by soil types
!
!     data wsat / .395, .410, .435, .485, .451, .420,
!    :            .477, .476, .426, .482, .482, 1.0E-20, 1.00/
!
! --- Cloud transmisivities (high, mid, low)
!
      data tac / 0.98, 0.85, 0.80 /          ! absorption
      data tsc / 0.80, 0.60, 0.48 /          ! scattering
      data clw / 0.06, 0.22, 0.26 /          ! enhancement coeff.

!     data csw /4.183E+06/

! --- Pressures differentiating cloud levels (mb)
!      DATA pcllo /800.0/                   ! low clouds
!      DATA pclhi /450.0/                   ! high clouds
!-----------------------------------------------------------------------
     if(debug) write(0,*)'Entering RADFLX ',jday, year, month,day,hour
     if(debug) write(0,*)'RADFLX time ',time, 'rcpinv = ',rcpinv
!-----------------------------------------------------------------------
! Initialize atmospheric transmissivity
!-----------------------------------------------------------------------
!      open(15,file='tran.data',form='unformatted',status='old',convert='big_endian')
!      read(15) tran
!      close(15)
!     write(*,*)
!     write(*,*) ' atmospheric transmissivity'
!     do k=1,npl
!       write(*,2000) (tran(i,k), i=1,npw)
!     enddo
!     write(*,*)
!2000 format(15f6.3)
!     stop
!-----------------------------------------------------------------------

      do j=1,jye
      IF ( debug ) write(0,*) ' j = ',j
      do i=1,ixe
      IF ( debug ) write(0,*) ' i = ',i
        cldfr(1,i,j) = 0.
        cldfr(2,i,j) = 0.
        cldfr(3,i,j) = 0.
        prw(i,j)     = 0.
        radlw(i,j)   = 0.
        radlwdn(i,j)   = 0.
        radlwup(i,j)   = 0.
        radsw(i,j)   = 0.
        tacc(i,j)    = 1.
        tscc(i,j)    = 1.
        psfc(i,j)    = (psl*100.)*(piinit(i,j,1) + p(i,j,1))**rcpinv
!
! --- Determine the cloud fraction
!
      do k=1,nz-1
        pi1   = (piinit(i,j,k) + p(i,j,k))**rcpinv
        pres  = (psl*100.)* pi1
        temp  = (piinit(i,j,k)+p(i,j,k))*t(i,j,k)
        qvs   = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pres
        IF ( .not. ns .ge. lc ) THEN
        rh = qv(i,j,k) / qvs
        ELSE
           rh = 0.0
           IF ( an(i,j,k,lc) > 1.0e-5 ) rh = 1.0
           IF ( li > 1 ) THEN
              IF ( an(i,j,k,li) > 1.0e-5 ) rh = 1.0
           ENDIF
           IF ( ls > 1 ) THEN
              IF ( an(i,j,k,ls) > 1.0e-4 ) rh = 1.0
           ENDIF
           IF ( lr > 1 ) THEN
              IF ( an(i,j,k,lr) > 0.3e-3 ) rh = 1.0
           ENDIF
        ENDIF
 
        cld1 = 0.5 + sign(0.5,pclhi - pres*.01)        ! hi cloud
        cld3 = 0.5 - sign(0.5,pcllo - pres*.01)        ! lo cloud
        cld2 = 1. - cld1 - cld3                  ! md cloud

        cldfr(3,i,j) = MAX(cld3 * rh, cldfr(3,i,j))
        cldfr(2,i,j) = MAX(cld2 * rh, cldfr(2,i,j))
        cldfr(1,i,j) = MAX(cld1 * rh, cldfr(1,i,j))
      enddo

        cldfr(3,i,j) = 4.0 * cldfr(3,i,j) - 3.0
        cldfr(2,i,j) = 4.0 * cldfr(2,i,j) - 3.0
        cldfr(1,i,j) = 2.5 * cldfr(1,i,j) - 1.5
!
! --- Limit the cloud fraction
!
      do m=1,3
        cldfr(m,i,j) = MIN( cldfr(m,i,j), 1. )
        cldfr(m,i,j) = MAX( cldfr(m,i,j), 0. )
      enddo
!-----------------------------------------------------------------------
! --- Compute the precipital water (cm)
!
        pi1 = (piinit(i,j,nz-1) + p(i,j,nz-1))**rcpinv
        pi2 = (piinit(i,j,nz-2) + p(i,j,nz-2))**rcpinv
        p1 = (psl*100.)* pi1
        p2 = (psl*100.)* pi2

        prw(i,j) = prw(i,j) -                                &
           0.5 * (qv(i,j,nz-1)+qv(i,j,nz-2)) * (p1-p2) * (0.1 / g)

      do k=nz-2,2,-1
        pi1 = (piinit(i,j,k+1) + p(i,j,k+1))**rcpinv
        pi2 = (piinit(i,j,k  ) + p(i,j,k  ))**rcpinv
        pi3 = (piinit(i,j,k-1) + p(i,j,k-1))**rcpinv
        p1 = (psl*100.)* pi1
        p2 = (psl*100.)* pi2
        p3 = (psl*100.)* pi3

        prw(i,j) = prw(i,j) - 0.5 * (    &
           0.5 * (qv(i,j,k)+qv(i,j,k+1)) * (p1-p2) * (0.1 / g) +  &
           0.5 * (qv(i,j,k)+qv(i,j,k-1)) * (p2-p3) * (0.1 / g) )
      enddo

        pi1 = (piinit(i,j,2) + p(i,j,2))**rcpinv
        pi2 = (piinit(i,j,1) + p(i,j,1))**rcpinv
        p1 = (psl*100.)* pi1
        p2 = (psl*100.)* pi2

        prw(i,j) = prw(i,j) -   &
           0.5 * (qv(i,j,1)+qv(i,j,2)) * (p1-p2) * (0.1 / g)
      enddo
      enddo
!-----------------------------------------------------------------------
!
! --- Calculate the longwave irradiance incident on the ground
!
!     if no clouds depcl = 1.   => emissa = 1.
!

!       year=coards(1)
!       month=coards(2)
!       day=coards(3)
!       hour=coards(4)
!       minute=coards(5)
!       second=coards(6)
       
!.....calculate julian day
      daymonth(1)=31
      daymonth(2)=28
      daymonth(3)=31
      daymonth(4)=30
      daymonth(5)=31
      daymonth(6)=30
      daymonth(7)=31
      daymonth(8)=31
      daymonth(9)=30
      daymonth(10)=31
      daymonth(11)=30
      daymonth(12)=31
      
      jday=0
      do i=month-1,1,-1
      jday=jday + daymonth(i)
      end do
      jday=jday+day

     if(debug) write(0,*) ' RADFLX: LW calc 2'

      fac = emissv + emissg - emissv * emissg ! emissv + emissg* (1.0 - emissv)
      do j=1,jye
      do i=1,ixe
        prwpt = max(prw(i,j), 0.1)
        epsair = 0.725 + 0.17*alog10(prwpt)

        depcl = 1. + ( clw(1)*cldfr(1,i,j) + clw(2)*cldfr(2,i,j)   &
                  +   clw(3)*cldfr(3,i,j) )
        emissa = MIN(epsair * depcl, 1.)

! - include vegetation in calculation
        pi1 = (piinit(i,j,1) + p(i,j,1))
        t1 = t(i,j,1) * pi1
        radlwdn(i,j) = (1. - veg(i,j)) * emissa * sbconst * t1**4
        radlwdn(i,j) = radlwdn(i,j) + veg(i,j) *  (   &
                    emissv * sbconst * tcanp(i,j)**4 +   &
                    (1. - emissv) * emissg * sbconst * tsfc(i,j)**4 ) &
                    / fac
!
! --- Calculate the longwave irradiance from the ground
!
        radlwup(i,j) = (1. - veg(i,j)) *  &
                          (emissg * sbconst * tsfc(i,j)**4 + &
                    (1. - emissg) * emissa * sbconst * t1**4 )
        radlwup(i,j) = radlwup(i,j) +  &
                 veg(i,j) * ( emissg * sbconst * tsfc(i,j)**4 + &
                 (1. - emissg) * emissv * sbconst * tcanp(i,j)**4) /  &
                  fac

        radlw(i,j) = radlwup(i,j) - radlwdn(i,j)
!      print*,radlw(i,j),radlwup(i,j),radlwdn(i,j),i,j
      enddo
      enddo

     if(debug) write(0,*) ' RADFLX: past LW calc'
!-----------------------------------------------------------------------
! --- Calculate the shortwave irradiance incident on the ground
!
      sunup = 0.
      radmax = 0.
!#ifndef MPI
!$omp  parallel do default(shared)  &
!$omp  private(i,j,dzdx,dzdy,slpdir,slpmag,anncyc,etau,adjust, &
!$omp        anglehr,soldec,coszen,sunup,sinpsi,cospsi,term1, &
!$omp        azimuth,dipsi,cosi,albz,temp1,temp2) &
!$omp reduction(+:radmax)
!#endif
      do j=1,jye
        jm1 = MAX(j-1,1)
        jp1 = MIN(j+1,jye)

      do i=1,ixe
        im1 = MAX(i-1,1)
        ip1 = MIN(i+1,ixe)

!     terrain slope magnitude and angle
!
! --- zsfc(*,2) defined at u-velocity locations
!
        dzdx = (zsfc(ip1,j  ) - zsfc(im1,j  ))*rdx*(1./float(ip1-im1))
        dzdy = (zsfc(i  ,jp1) - zsfc(i  ,jm1))*rdx*(1./float(jp1-jm1))
        slpmag = acos(1./SQRT(1. + dzdx * dzdx + dzdy * dzdy))
        slpdir = 0.5 * pii * sign(1.,dzdx)
!       slpdir = atan2( dzdx / dzdy )

!     hour in the day, hour => initial hour !!


!!!! msb 02/13/08 check calculation
!!!  added minutes and seconds into calculation for solar hour

!        thour(i,j) = float(hour) + time/3600. +  glon(i,j)/15.
        thour(i,j) = float(hour) + float(minute)/60. +  &
       float(second)/3600. + time/3600. +glon(i,j)/15.
     
        thour(i,j) = MOD( thour(i,j) , hrday )

!         print*,float(hour),float(minute),float(second),glon(i,j),time

        cjdayloc(i,j) = jday + int ( thour(i,j)/hrday )

        term1 = 0.5 - sign(0.5, thour(i,j))
        thour(i,j) = thour(i,j) + term1 * hrday
        cjdayloc(i,j) = MOD( cjdayloc(i,j)-term1,daysyr)

!.......msb 6/2/08 try keeping solar rad const.
!        thour(i,j)=16.03556

!        print*,hour,minute,second,term1,thour(i,j),time

!     Compute the zenith angle

         anncyc = 2.*pii*( cjdayloc(i,j) -1.)/daysyr

         adjust = 1.000110 + 0.034221 * cos(anncyc) &
                          + 0.001280 * sin(anncyc)  &
                          + 0.000719 * cos(2.*anncyc) &
                          + 0.000077 * sin(2.*anncyc)

         etau = 0.158 * sin(pii*(cjdayloc(i,j)+10.)/91.25) +  &
               0.125 * sin(pii*cjdayloc(i,j) / 182.5 )

         anglehr = 15. * (degtorad) * (thour(i,j) + etau - 12.)

         soldec = 23.5 * (degtorad)*  &
                 cos(2.*pii*(cjdayloc(i,j)-173.)/daysyr)

!        coszen = cos(glat(i,j)*degtorad)*cos(soldec)*cos(anglehr) &
!              + sin(glat(i,j)*degtorad)*sin(soldec)

         coszen = cos(glat(i,ny/2)*degtorad)*cos(soldec)*cos(anglehr) &
               + sin(glat(i,ny/2)*degtorad)*sin(soldec)

         sunup = 0.5 + sign(0.5, coszen)
         zenith(i,j) = acos(coszen)

!         print*,zenith(i,j),coszen,glat(i,ny/2),thour(i,j),anglehr,i,j

!
! --- Path length calculation
!
         path(i,j) = 1. / MAX( coszen, 0.00005)
!-----------------------------------------------------------------------
! --- Effects of terrain slope and angle on radiation

         sinpsi = cos(soldec)*sin(anglehr)/sin(zenith(i,j))

!        cospsi = (coszen * sin(glat(i,j)*degtorad) - sin(soldec)) &
!              / (sin(zenith(i,j))*cos(glat(i,j)*degtorad))

         cospsi = (coszen * sin(glat(i,ny/2)*degtorad) - sin(soldec)) &
               / (sin(zenith(i,j))*cos(glat(i,ny/2)*degtorad))

         azimuth = atan2(sinpsi,cospsi)

         dipsi = azimuth - slpdir

         cosi = cos(slpmag) * coszen  &
             + sin(slpmag) * sin(zenith(i,j)) * cos(dipsi)
!        if( i .EQ. 1 ) then
!        write(0,*) coszen,anglehr,thour(i,j),cosi/coszen
!        endif

!        radsw(i,j) = solarc * adjust * (cosi / coszen) * sunup
         radsw(i,j) = solarc * adjust * coszen * sunup
         radmax = radmax + radsw(i,j)

!-----------------------------------------------------------------------
!
! --- Albedo calculations  - zenith effects 
!    
         temp1 = SQRT( (zenith(i,j)*180./pii)**3 )
         albz = 0.01 * ( exp(0.003286 * temp1) -1. )

!     Albedo calculations  - soil moisture effects
         if( Nint(stype(i,j)) .eq. 12 ) then
           temp1 = 0.
         else
           temp1 = wsfc(i,j) / wsat(int(stype(i,j)))
         endif

         temp2 = 0.5 + sign(.5,temp1-.5)
         albedo(i,j) = albz + temp2 * 0.14  &
               + (1.-temp2)*(.31-.34*temp1)

!          print*,albedo(i,j),albz,temp2,temp1,zenith(i,j)

!        albedo(i,j) = (1.-veg(i,j)) * albedo(i,j) + veg(i,j) * albedov

!.......msb 4/9/08 albedo such that (1-albedo)=(1-albedog)(1-veg)
!       and radsw=(1-albedog)(1-veg)*sh_down (eq. 36)
!        albedo(i,j) = (1.-veg(i,j)) * albedo(i,j) + veg(i,j)


!-----------------------------------------------------------------------
      enddo
      enddo
     if(debug) write(0,*) ' RADFLX: past SW calc'
!-----------------------------------------------------------------------
! --- Uniform surface heating
      iradloc=nx/2
      jradloc=ny/2
      if( pbl_run ) then
        do j=1,jye
        do i=1,ixe
          radsw(i,j) = radsw(iradloc,jradloc)
          path(i,j) = path(iradloc,jradloc)
        enddo
        enddo
      else
     if(debug) write(0,*) ' RADFLX: radsw changes in x '
        do j=1,jye
        do i=1,ixe
          radsw(i,j) = radsw(i,jradloc)
          path(i,j) = path(i,jradloc)
        enddo
        enddo
      endif
!-----------------------------------------------------------------------
! --- If the sun is not up skip solar radiation calc
!
     if(debug) write(0,*) ' RADFLX: radmax: ',radmax
      if( radmax .GT. 0. ) then
!
! --- SW transmissivity determined from look-up table from
!     shortwave model of carlson and boland (JAM, 1978)
!
!     direct radiation : use actual solar path length
!
      CALL TRANSM(path, prw, ftabs , ftsca , fbsc , psfc, tran, &
                 npw, npl, nx, ny, ixe, jye)
!
!     diffuse radiation : use path length (diffusivity factor) = 1.67
!       from rodgers and walshaw, 1967
!
      do j=1,jye
      do i=1,ixe
        path(i,j) = 1.67
      enddo
      enddo
      CALL TRANSM(path, prw, ftabsd, ftscad, fbscd, psfc, tran, &
                 npw, npl, nx, ny, ixe, jye)

      do j=1,jye
      do i=1,ixe
      do m=1,3
        tacc(i,j) = tacc(i,j) * ( 1. - ((1.-tac(m)) * cldfr(m,i,j)))
        tscc(i,j) = tscc(i,j) * ( 1. - ((1.-tsc(m)) * cldfr(m,i,j)))
      enddo
!
! --- minimum absorption and scattering transmissivities through cloud
!     from drummond and hickey 1971
!
        tacc(i,j) = MAX( tacc(i,j), 0.70 )
        tscc(i,j) = MAX( tscc(i,j), 0.44 )
!-----------------------------------------------------------------------
     if(debug) write(0,*) ' Set net sfc rad '
! --- Compute the net radiation at the surface
!
        bdbar = ( fbscd(i,j) * (1. - ftscad(i,j)) + (1. - tscc(i,j)) ) / &
                           ((1. - ftscad(i,j)) + (1. - tscc(i,j)) )
        xser  = bdbar * albedo(i,j) * (1. - ftscad(i,j) * tscc(i,j) ) * &
                                        ftabsd(i,j) * tacc(i,j)
        hi    = radsw(i,j) * ftabs(i,j) * tacc(i,j) * tscc(i,j) * &
               ( ftsca(i,j) + (1.-ftsca(i,j)) * (1.-fbsc(i,j)) )

        radsw(i,j) = hi * (1. - albedo(i,j)) / (1.-xser)


!.......msb 4/9/08  add in effects of vegetation
        radsw(i,j)=radsw(i,j)*(1.0-veg(i,j))
!........msb 5/14/10 delete line        



!      print*,radsw(i,j),i,j
 
      enddo
      enddo
!
! --- fix circulation @ western boundary
!
!     if( .not. pbl_run ) then
!     do j=1,ny-1
!         radsw(j,1) = radsw(j,3)
!         radsw(j,2) = radsw(j,3)
!     enddo
!     endif

      endif   ! endif for radmax
!---------------------------------------------------------------------- 
!     if(debug) then
!     call amaxmin( radlw,v1,v2,i1,j1,k1,i2,j2,k2,0,0,0,nx,ny,1)
!     write(0,*) ' RADLW: ',v1,v2
!     call amaxmin( radlwup,v1,v2,i1,j1,k1,i2,j2,k2,0,0,0,nx,ny,1)
!     write(0,*) ' RADLWUP: ',v1,v2
!     call amaxmin( radsw,v1,v2,i1,j1,k1,i2,j2,k2,0,0,0,nx,ny,1)
!     write(0,*) ' RADSW: ',v1,v2,i1,j1,i2,j2
!     call amaxmin( path,v1,v2,i1,j1,k1,i2,j2,k2,0,0,0,nx,ny,1)
!     write(0,*) ' PATH: ',v1,v2,i1,j1,i2,j2
!     endif
!-----------------------------------------------------------------------
      deallocate( cldfr )

      deallocate(  radlwup )
      deallocate(   zenith ); deallocate( thour )
      deallocate( cjdayloc )
      deallocate(      prw ); deallocate(  psfc )
      deallocate(    ftabs ); deallocate( ftsca )
      deallocate(     fbsc )
      deallocate(   ftabsd ); deallocate( ftscad )
      deallocate(    fbscd )
      deallocate(     path ); deallocate( tacc )
      deallocate(     tscc )
!-----------------------------------------------------------------------
      if(debug) write(0,*) ' Exiting SWRAD '
!     if(debug) stop
      return
      END

      SUBROUTINE TRANSM(path,prw,ftabs,ftscat,fbscat,psfc,tran, &
                       npw,npl,nx,ny,ixe,jye)
!-----------------------------------------------------------------------
!     This subroutine obtains sw transmissivity as a function of      
!       1) path length (1-10) (path length = 1 when zenith angle = 0),
!       2) precipitable water (0-5 cm)                                
!     using bilinear interpolation from look-up tables determined     
!     a shortwave model of carlson and boland (1978, JAM)             
!
!     path   : path length.                                - input    
!     prw    : precipitable water in cm.                   - input    
!
!     ftabs  : absorption transmissivity (0 - 1).          - output   
!     ftscat : scattering transmissivity (0 - 1).          - output   
!     fbscat : backscattering coefficient (0 - 1).         - output   
!              (percent of scattered radiation which                  
!               is backscattered. the rest is forward                 
!               scattered.)
!-----------------------------------------------------------------------
      
      implicit none
      
      integer nx,ny,npw,npl,ixe,jye
      
      real path(nx,ny), ftabs(nx,ny)
      real ftscat(nx,ny), fbscat(nx,ny)
      real prw(nx,ny), psfc(nx,ny)
      real tran(npw,npl)
      
      integer :: i,j
      real    :: fract,fract1, fract2
      integer :: ipath, jpath, iom, jom
      real    :: gract, gract1, gract2
      real    :: fta, fts, sp
      
!---------------------------------------------------------------------- 
      do j = 1,jye
      do i = 1,ixe
        path(i,j) = MIN(path(i,j),9.99)
        path(i,j) = MAX(path(i,j),1.01)
        prw(i,j)  = MIN(prw(i,j),4.99)
        prw(i,j)  = MAX(prw(i,j),0.01)
      enddo
      enddo
!                                                                       
!-----46 path lengths, 1-10, interval=0.2                              
!                                                                      
      do j = 1,jye
      do i = 1,ixe
       fract = 5.*(path(i,j)-1.) + 1.
       ipath = IFIX(fract+.001)
       jpath = ipath + 1
       fract1 = fract-ipath
       fract2 = 1.-fract1
!                                                                    
!-----11 precip water amts., 0-5 cm, interval=0.5                   
!                                                                  
      gract = 2.*prw(i,j) + 1.
      iom = IFIX(gract+.001)
      jom = iom + 1
      gract1 = gract-iom
      gract2 = 1.-gract1
!                                                                
!-----interpolate to get abs, scat, bscat values                
!        - bilinear interpolation for abs (f(path,omega)       
!        - linear interpolation for scat, bscat (f(path only))
!                                                            
      fta    = fract2*gract2*tran(iom,ipath)  &
            + fract1*gract1*tran(jom,jpath)   &
            + fract1*gract2*tran(iom,jpath)   &
            + fract2*gract1*tran(jom,ipath)
      fts         = fract2*TRAN(12,ipath) + fract1*tran(12,jpath)
      fbscat(i,j) = fract2*TRAN(13,ipath) + fract1*tran(13,jpath)
!                                                                
!-----transmissions from table are valid for surface pressure of
!     1013.25 mb.  we correct below for the actual surface pressure
!     assuming that attenuative effects of water vapor are        
!     secondary to those of other active constituents          
!     (air molecules, aerosols, ozone)                        
!                                                            
      sp =  psfc(i,j) / 101325.
      ftabs(i,j) = sp*(fta-1.) + 1.
      ftscat(i,j)= sp*(fts-1.) + 1.
      enddo
      enddo
!---------------------------------------------------------------------- 
      return
      END
