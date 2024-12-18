!-----------------------------------------------------------------------------------------------------------------------------------
!  
!  SPHYS_MODULE:  parameters for surface and radiation physics 
!                 replaces stuff from const.h
!

MODULE SPHYS_MODULE
   
  implicit none
      
      real, parameter :: psl = 1000.  ! sea-level pressure
      real, parameter :: b61 = 0.608
      real, parameter :: rv = 461.5
      real, parameter :: tzero = 273.15
      real, parameter :: rhow = 1000.              ! water density per m^3
      real, parameter :: eradius = 6.371E+06       ! Radius of Earth in meters
      real, parameter :: rcpinv = 1004./287.04
!      real, parameter :: cpi = 3.141592654
!      real, parameter :: deg2r = cpi / 180.
! note rcvinv replaced by cvr


      real, parameter :: cloudthreshold2 = 2.0e-8

      real, parameter :: emissg = 0.920  ! Emissivity of the ground (assume a sandy soil)

      real, parameter :: emissv = 0.940, albedov = 0.25   ! Emissivity / Albedo of the vegetation (grass)

      real, parameter :: sbconst = 5.67E-08     ! Stephan-Boltzman const

      real, parameter :: omega = 7.27220521664E-05

      real, parameter :: daysec = 86400.      ! seconds per day

      real, parameter :: solarc  = 1365.       ! solar constant
!      parameter(solarc = 1200.)                    ! Model too hot
!      parameter(solarc = 1250.)                    ! Model too hot
!      parameter(solarc = 1300.)                    ! Model too hot
!      parameter(solarc = 1365.)                    ! 100 m vert. resolution
!      parameter(solarc = 1395.)                    !  10 m vert. resolution

      real, parameter :: Pr0 = 0.74           ! Neutral stability Prandtl number

      real, parameter :: vonk = 0.40          ! Von Karman constant

      real, parameter :: qrthres = 1.E-06
      real, parameter :: vlimit = 1.0

      real, parameter :: daysyr = 365.

      real, parameter :: hrday = 24.
      
      integer      nsm, isfcphy_opt, icor_opt, imoist_opt,  &
                   jday !,  &
!                    year, month, day, hour, minute, second

      real  dts, rdx, rdy, rdz, udom, vdom,  &
                   udomain, vdomain, cs, cssq, alowk,  &
                   f2dmpx, f2dmpz, f4dmpx, f4dmpz, &
                   tstart, tstop, &
                   tplot, tprint, tsave
  
      real :: pcllo = 800.0
      real :: pclhi = 200.0 ! ERM changed from 450 to 200 to be closer to typical anvil level 2/22/2012
      
      real :: albedo_init = 0.20
      real :: veg_init = 0.2 ! 0.84
                    
END MODULE SPHYS_MODULE


!==================================================================================================================
!==================================================================================================================


      subroutine sfcphy(an,u,fu,v,fv,t,ft,qv,fqv,p,km,kh,piinit,  & 
           tcanp, wcanp, qav, eflx, fflx, precip, veg,         &
           g12,sgz,wgz,glat,glon,                              &
           zsfc,stype,tsfc,wsfc,tsoil,wsoil,vlai,              &
           uflx,vflx,tflx,qflx,radsw,radlw,albedo,rough,       &
           sfctke,time,nx,ny,nz,ns,dt)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      USE GRID_MODULE
      USE PARAM_MODULE
      USE SPHYS_MODULE
      USE COMMASMPI_MODULE
               
      implicit none
      
      integer nx,ny,nz,ns
      integer time
      
      real dt

      real  an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      real  u(-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  fu(-ng+1:nx+ng,-ng+1:ny+ng) 
      real  v(-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  fv(-ng+1:nx+ng,-ng+1:ny+ng) 
      real  t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  ft(-ng+1:nx+ng,-ng+1:ny+ng)
      real  qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  fqv(-ng+1:nx+ng,-ng+1:ny+ng)
      real  p(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  km(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  kh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real  piinit(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

      real precip(-ng+1:nx+ng,-ng+1:ny+ng)
      real sfctke(-ng+1:nx+ng,-ng+1:ny+ng)
      
      real g12(-ng+1:nx+ng,-ng+1:ny+ng)

      real sgz(nz), wgz(nz), mfc(nz), mfe(nz)

      real glat(-ng+1:nx+ng,-ng+1:ny+ng)
      real glon(-ng+1:nx+ng,-ng+1:ny+ng)
      real zsfc(-ng+1:nx+ng,-ng+1:ny+ng)

      real stype(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsfc(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsfc(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real tcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real wcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real qav(-ng+1:nx+ng,-ng+1:ny+ng)
      real veg(-ng+1:nx+ng,-ng+1:ny+ng)
      real vlai(-ng+1:nx+ng,-ng+1:ny+ng)

      real uflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real vflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real tflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real qflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real eflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real fflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real radsw(-ng+1:nx+ng,-ng+1:ny+ng)
      real radlw(-ng+1:nx+ng,-ng+1:ny+ng)
      real albedo(-ng+1:nx+ng,-ng+1:ny+ng)
      real rough(-ng+1:nx+ng,-ng+1:ny+ng)

      real mbar, nsum, npsum
!      real lmax
      real  kmsfc, kmtop

      real tflx1, qflx1, uflx1, vflx1
      
      integer i,j
      
      real pi1, pres, thv, tv, rho, cpm
      real delz1, delz


      integer :: ixe,jye

!-----------------------------------------------------------------------
      logical debugsfc
!-----------------------------------------------------------------------
!     real, dimension(:,:), allocatable :: dum
!     allocate( dum(ny,nx) )
!-----------------------------------------------------------------------
      debugsfc = .false.
!-----------------------------------------------------------------------
      if(debugsfc) write(*,*) ' Enter sfcphy '


      ixe = itile
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend - jybeg

      eflx(:,:) = 0.
      fflx(:,:) = 0.


!       year=coards(1)
!       month=coards(2)
!       day=coards(3)
!       hour=coards(4)
!       minute=coards(5)
!       second=coards(6)
                                          
! --- tst is the number of time steps in 1 minute. 
! --- Changing tst will change the frequenct of the convective
!       initiation check.

!     tst =  dt / dt
!     ntst = INT(tst)
!     ntst = MAX(ntst,1)
!     iternum = INT( time / dt )
!     icheck=MOD(iternum, ntst)
!     if(icheck .EQ. 0) then

!       do j=1,ny
!       do i=1,nx
!       print*,tsfc(i,j),t(1,i,j),qv(1,i,j),i,j,time
!       end do
!       end do


!     write(0,"('--------------SFCPHY: before vegpbl --------------')")
!     do i=1,nx
!     do j=1,ny
!     do k=1,nz
!     IF( .not. ( u(k,i,j) .lt. 500. .and.  u(k,i,j) .gt. -500. ) ) then
!     print*, u(k,i,j), t(k,i,j), qv(k,i,j),k,i,j,time
!     end if
!     end do
!     end do
!     end do


        call vegpbl(an,u,v,t,qv,p,piinit,    &
           tcanp, wcanp, qav, eflx, fflx, precip, veg,  &
           g12,glat,glon,                               &
           zsfc,stype,tsfc,wsfc,tsoil,wsoil,vlai,       &
           uflx,vflx,tflx,qflx,radsw,radlw,albedo,rough, &
           sfctke,sgz,wgz,mfc,mfe,time,nx,ny,nz,ns,dt)

!    endif


!      print*,'after vegpbl'
!-----------------------------------------------------------------------
      do j=1,jye
      do i=1,ixe
        pi1 = (piinit(i,j,1) + p(i,j,1))
        pres = (psl*100.) * pi1**cvr
        thv = t(i,j,1) * (1. + 0.61 * qv(i,j,1) )
        tv = thv * pi1
        
        rho = pres / (rd * tv )
        cpm = cp * (1. + 0.8 * qv(i,j,1) )

!        print*,rho,pres,rd,tv,pi1,thv,t(i,j,1),qv(i,j,1),i,j
!
! --- surface layer tendencies 
!
          delz1 = sgz(1) / g12(i,j)
          delz = g12(i,j) / (sgz(2)-sgz(1))
!         delz =       1. / (sgz(2)-sgz(1))    ! mix on computational grid

!           delz1=100.0/g12(i,j)
!           delz=g12(i,j)/100.0


!         /msb/ tflx(i,j) is potential temp flux
!
!         msb 3/21/08 multiply by cpm*rho to get units of w m-2 s-1
!                     tflx(i,j) is already in units w m-2 s-1

!          tflx1 = kh(i,j,2) * delz * (t(i,j,1)-t(i,j,2))
          tflx1 = kh(i,j,2) * delz * cpm * rho * (t(i,j,1)-t(i,j,2))
!          tflx1 = Max(50.0, kh(i,j,2)) * delz * cpm * rho * (t(i,j,1)-t(i,j,2))
          ft(i,j)  = ft(i,j) +                   &
                  (tflx(i,j) - tflx1)/(cpm*rho*delz1)


!          print*,'fluxes',tflx1,tflx(i,j),t(i,j,1),t(i,j,2),ft(i,j)

!..........msb 4/12/08 multiply qflx1, uflx1, vflx1 by rho to correct units

!          qflx1 = kh(i,j,2) * delz * (qv(i,j,1)-qv(i,j,2))
          qflx1 = kh(i,j,2) * delz * rho * (qv(i,j,1)-qv(i,j,2))
          
!         dum(i,j) = qflx1
          fqv(i,j) = fqv(i,j) +   &
                  ( qflx(i,j) - qflx1 ) / (rho*delz1)

!          uflx1 = km(i,j,2) * delz * (u(i,j,1) - u(i,j,2))
          uflx1 = km(i,j,2) * delz * rho * (u(i,j,1) - u(i,j,2))
          
          fu(i,j)  =  fu(i,j) +    &
                  (uflx(i,j) - uflx1 ) / (rho*delz1)

!          vflx1 = km(i,j,2) * delz * (v(i,j,1) - v(i,j,2))
          vflx1 = km(i,j,2) * delz * rho * (v(i,j,1) - v(i,j,2))
          
          fv(i,j)  =  fv(i,j) +   &
                  (vflx(i,j) - vflx1 ) / (rho*delz1)


!.......test turning off sfc fluxes

!          ft(i,j)=0.0
!          fqv(i,j)=0.0



!        print*,qflx1,qflx(i,j),fqv(i,j,1),qv(i,j,1),qv(i,j,2)
!        print*,kh(i,j,2),delz,rho,cpm,i,j
        
 1015 format('sfcflx: ',7(f7.2,1x),i5)
 
!         if((mod(time,4).eq.0).and.i.eq.30.and.j.eq.25) then
!         IF ( my_rank == 0 .and. i == ixend-1 .and. j == jyend-1 ) THEN
!         write(91,1015) tflx(i,j),qflx(i,j)*Lav,     &
!        tsfc(i,j),tsoil(i,j),t(i,j,1),wsfc(i,j),wsoil(i,j),time
!         end if
 
!........msb 5/9/08  turn off momentum fluxes

         fu(i,j)=0.0
         fv(i,j)=0.0

      enddo ! j
!        fv(1,ny,i)  =  0.5 * (3.*fv(1,ny-1,i) - fv(1,ny-2,i))
!         fv(1,ny,i)=0.0
      enddo
!      do j=1,jye
!        fu(1,j,nx)  =  0.5 * (3.*fu(1,j,nx-1) - fu(1,j,nx-2))
!         fu(1,j,nx)=0.0
!      enddo

!........msb 5/9/08  turn off momentum fluxes

         fu(:,:)=0.0
         fv(:,:)=0.0



!     i0 = 140
!     i1 = 150
!     write(*,*) ' FQV '
!     write(*,999) (fqv(1,10,i),i=i0,i1)   
!     write(*,*) ' QFLX (k=2) '
!     write(*,999) (dum(10,i),i=i0,i1)   
!     write(*,*) ' SFC QFLX '
!     write(*,999) (qflx(10,i),i=i0,i1)   
!     write(*,*) ' KH (k=2) '
!     write(*,999) (kh(2,10,i),i=i0,i1)   
!     write(*,*) ' KH (k=1) '
!     write(*,999) (kh(1,10,i),i=i0,i1)   
!     write(*,*) ' QV (k=2) '
!     write(*,999) (qv(2,10,i),i=i0,i1)   
!     write(*,*) ' QV (k=1) '
!     write(*,999) (qv(1,10,i),i=i0,i1)   
!999  format(12e13.5)
!     stop
!-----------------------------------------------------------------------
!     deallocate( dum )
!-----------------------------------------------------------------------
      if(debugsfc) write(*,*) ' past sfcphy  '
      return
      END


!==================================================================================================================
!==================================================================================================================

      subroutine set_var_sfcphy(tcanp,   &
                 wcanp,qav,veg,stype,    &
                 tsrfc,wsfc,tsoil,wsoil, &
                 vlai, albedo,rough,     &
                 nx,ny,nz,ns,s,          &
!                 nxend,nyend,            &
                 eflx, fflx, uflx,       &
                 vflx, tflx, qflx,       &
                 radsw, radlw,piinit,pi,tsfc,qsfc) 
                

      USE COMMASMPI_MODULE, only : itile,jtile, nxend, nyend, myproci, myprocj, &
                                   nproci, nprocj, ixend, ixbeg, jyend, jybeg
      USE GRID_MODULE
      USE SPHYS_MODULE, only: albedo_init, veg_init

      implicit none

      integer i,j,k,nx,ny,nz,ns

!      real tsrfc(nx,ny), wsfc(nx,ny)
!      real tsoil(nx,ny), wsoil(nx,ny)
!      real tcanp(nx,ny), wcanp(nx,ny), qav(nx,ny)
!      real veg(nx,ny), vlai(nx,ny)
!      real albedo(nx,ny), rough(nx,ny)
!      real stype(nx,ny)
!      real t(nx,ny,nz),qv(nx,ny,nz)
 
      TYPE(VARIABLE)     :: tcanp, wcanp
      TYPE(VARIABLE)     :: qav, veg
      TYPE(VARIABLE)     :: stype, tsrfc
      TYPE(VARIABLE)     :: wsfc, tsoil
      TYPE(VARIABLE)     :: wsoil, vlai
      TYPE(VARIABLE)     :: albedo, rough 
      TYPE(VARIABLE)     :: s(ns)
      TYPE(VARIABLE)     :: eflx, fflx
      TYPE(VARIABLE)     :: uflx, vflx
      TYPE(VARIABLE)     :: tflx, qflx
      TYPE(VARIABLE)     :: radsw, radlw
      TYPE(VARIABLE)     :: piinit, pi
      
      real               :: tsfc,qsfc
      integer            :: ixe,jye,i1,j1,ix,jy
      integer  :: iranseed
      real     :: rndnum

      real, allocatable :: ranarray(:,:),ranarraysmth(:,:)
                                                                
      real pib,pip,pi1
      


      allocate ( ranarray(1:nxend,1:nyend) )
      allocate ( ranarraysmth(nxend,1:nyend) )

!.....fill 2D arrays

      iranseed = 1
  ranarray(:,:) = 0.0
  ranarraysmth(:,:) = 0.0
      DO j = 1,nyend-1
       DO i = 1,nxend-1
         ranarray(i,j) = 2.0*(rndnum(iranseed) - 0.5) ! values of -1 to 1
       ENDDO
      ENDDO

    DO j = 1,nyend-2,4
     DO i = 1,nxend-2,4
       DO  j1 = 0,3
       DO  i1 = 0,3
         ix = Min(nxend,i+i1)
         jy = Min(nyend,j+j1)
         ranarraysmth(ix,jy) = ranarray(i,j)
       ENDDO
       ENDDO
     ENDDO
    ENDDO
      
      ixe = itile
      IF ( myproci == nproci ) ixe = ixend-ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend-jybeg
      
      do j=1,jye
      do i=1,ixe

      pib = .5 * (3.*piinit%flt1d(1) - piinit%flt1d(2))
      pip = .5 * (3.*pi%flt3d(i,j,1) - pi%flt3d(i,j,1))
      pi1 = (pib + pip)   

      stype%flt2d(i,j)=3.
      tsrfc%flt2d(i,j)= tsfc !314.75     !T, not THETA
!      tsrfc%flt2d(i,j)=s(1)%flt3d(i,j,1)*pi1
!      wsfc%flt2d(i,j)=s(2)%flt3d(i,j,1)
      wsfc%flt2d(i,j)=0.11
      tsoil%flt2d(i,j)= tsfc ! 305.25 
!      wsoil%flt2d(i,j)=0.25                  
      wsoil%flt2d(i,j)=0.11
      tcanp%flt2d(i,j)=tsfc !s(1)%flt3d(i,j,1)*pi1   ! T, not THETA
!      print*,tcanp%flt2d(i,j),s(1)%flt3d(i,j,1),pi1,i,j
!      wcanp%flt2d(i,j)=s(2)%flt3d(i,j,1)
      wcanp%flt2d(i,j) = 0.0
      qav%flt2d(i,j)=s(2)%flt3d(i,j,1)
!      veg%flt2d(i,j) = veg_init + (1.0 - veg_init)*ranarraysmth(ixbeg-1+i,jybeg-1+j)
      veg%flt2d(i,j) = veg_init + 0.05*ranarraysmth(ixbeg-1+i,jybeg-1+j)
!.....msb 5/6/08.....try LAI 2*veg, not 7*veg           
!      vlai%flt2d(i,j)=7.0*veg%flt2d(i,j)
      vlai%flt2d(i,j)=0.7*veg%flt2d(i,j)
      albedo%flt2d(i,j) = albedo_init
      rough%flt2d(i,j)=.05

      eflx%flt2d(i,j)=0.0
      fflx%flt2d(i,j)=0.0
      uflx%flt2d(i,j)=0.0
      vflx%flt2d(i,j)=0.0
      tflx%flt2d(i,j)=0.0
      qflx%flt2d(i,j)=0.0      
      radsw%flt2d(i,j)=0.0
      radlw%flt2d(i,j)=0.0      
                              

!      print*,tsrfc%flt2d(i,j),wsfc%flt2d(i,j),i,j

      end do
      end do

      Deallocate( ranarray, ranarraysmth )
      return
      END           
                
!==================================================================================================================
!==================================================================================================================
                
      subroutine gf(var,nx,ny,nz,ng)

      implicit none
      
      integer nx,ny,nz,ng
      real var(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      integer i,j,k
      
      do i=1,nx
      do j=1,ny
      do k=1,nz
      print*,var(i,j,k),k,j,i,'b'
      end do
      end do
      end do

      
      return
      END      
                      

!==================================================================================================================
! msb 2/8/08 Added surfcae layer parameterization
!
!PART-VI: SURFACE LAYER
  
      subroutine sfc_interface(gd,     &
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    u, v, w,                &
                    s, ns,                  &
                    an,                     &
                    precip,                 &
                    lat,lon, gtx,gty,       &
                    time,nx,ny,nz,dt,ugrid,vgrid,nst,dx,dy)

   USE GRID_MODULE
   USE PARAM_MODULE
   USE CPUTIME_MODULE
   USE COMMASMPI_MODULE
   
   implicit none
   
   integer            :: nx,ny,nz
   integer            :: ns
   
   TYPE(GRID)         :: gd 
   TYPE(VARIABLE)     :: u
   TYPE(VARIABLE)     :: v
   TYPE(VARIABLE)     :: w
   TYPE(VARIABLE)     :: pi, piinit
   TYPE(VARIABLE)     :: km, kminit
   TYPE(VARIABLE)     :: s(ns)     ! We are going to try and create arrays of the scalar variables needed here
   real               :: precip(-ng+1:nx+ng,-ng+1:ny+ng,5)

   real an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
   
   
   real    :: dt
   integer :: time
   real    :: lat, lon
   real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng)
   real    :: ugrid,vgrid
   integer :: nst
   real    :: dx,dy
   
!   real    :: glat(-ng+1:nx+ng,-ng+1:ny+ng)
!   real    :: glon(-ng+1:nx+ng,-ng+1:ny+ng)
   real    :: sgz(nz),wgz(nz)

   real, allocatable,dimension(:,:)     ::  zsfc,g12,sfctke
   real, allocatable,dimension(:,:,:)   ::  kh,piinit3d,tmp3d

   real, allocatable   ::  ji6(:,:)
   real, allocatable   ::  ji10(:,:),ji11(:,:)

   real, allocatable :: mfu(:,:)
   real, allocatable :: mfv(:,:) 
   real, allocatable :: mft(:,:) 
   real, allocatable :: mfq(:,:) 

   real, allocatable :: utmp(:,:,:)
   real, allocatable :: vtmp(:,:,:)
   
   real, pointer :: zc(:)

  integer :: tcanp
  integer :: wcanp
  integer :: qav
  integer :: veg
  integer :: stype
  integer :: tsrfc
  integer :: wsfc
  integer :: tsoil
  integer :: wsoil
  integer :: vlai
  integer :: albedo
  integer :: rough
  integer :: qv
  integer :: eflx
  integer :: fflx
  integer :: uflx
  integer :: vflx
  integer :: tflx
  integer :: qflx
  integer :: radsw
  integer :: radlw
   
   
   integer i,j,k,n,num2d,nb
   integer, parameter :: ngtmp=3
   
   integer m(50)
   
   integer, save :: istart = 0
   real :: pmax,zfac

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag

      real, parameter :: zlen = 100.0

! -------------------------------------------------------

      CALL GET_VARIABLE(gd, 'ZC',   zc)

 tcanp   = GET_VARIABLE_INDEX(gd, 'TCANP')
 wcanp   = GET_VARIABLE_INDEX(gd, 'WCANP')
 qav     = GET_VARIABLE_INDEX(gd, 'QAV')
 veg     = GET_VARIABLE_INDEX(gd, 'VEG')
 stype   = GET_VARIABLE_INDEX(gd, 'STYPE')
 tsrfc   = GET_VARIABLE_INDEX(gd, 'TSRFC')
 wsfc    = GET_VARIABLE_INDEX(gd, 'WSFC')
 tsoil   = GET_VARIABLE_INDEX(gd, 'TSOIL')
 wsoil   = GET_VARIABLE_INDEX(gd, 'WSOIL')
 vlai    = GET_VARIABLE_INDEX(gd, 'VLAI')
 albedo  = GET_VARIABLE_INDEX(gd, 'ALBEDO')
 rough   = GET_VARIABLE_INDEX(gd, 'ROUGH')

 eflx    = GET_VARIABLE_INDEX(gd, 'EFLX')
 fflx    = GET_VARIABLE_INDEX(gd, 'FFLX')
 uflx    = GET_VARIABLE_INDEX(gd, 'UFLX')
 vflx    = GET_VARIABLE_INDEX(gd, 'VFLX')
 tflx    = GET_VARIABLE_INDEX(gd, 'TFLX')
 qflx    = GET_VARIABLE_INDEX(gd, 'QFLX')  
 radsw   = GET_VARIABLE_INDEX(gd, 'RADSW')  
 radlw   = GET_VARIABLE_INDEX(gd, 'RADLW') 
    
    allocate( mfu(-ng+1:nx+ng,-ng+1:ny+ng), &
              mfv(-ng+1:nx+ng,-ng+1:ny+ng), &
              mft(-ng+1:nx+ng,-ng+1:ny+ng), &
              mfq(-ng+1:nx+ng,-ng+1:ny+ng) )
              
      mfu(:,:) = 0.
      mfv(:,:) = 0.
      mft(:,:) = 0.
      mfq(:,:) = 0.   ! set fw to be fqv                                

      allocate( utmp(-ng+1:nx+ng,-ng+1:ny+ng,2) )
      allocate( vtmp(-ng+1:nx+ng,-ng+1:ny+ng,2) )
      
      DO k = 1,2
        DO j = -ng+1,ny+ng
          DO i = -ng+1,nx+ng
           utmp(i,j,k) = u%flt3d(i,j,k) + ugrid
           vtmp(i,j,k) = v%flt3d(i,j,k) + vgrid
          ENDDO
        ENDDO
      ENDDO

!....set piinit,glat,glon,g12,sgz,wgz,zsfc,kh


 !     CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
 !     CALL GET_VARIABLE(gd, 'ZEDX',   zedx)

   allocate( g12(-ng+1:nx+ng,-ng+1:ny+ng) )
   allocate( zsfc(-ng+1:nx+ng,-ng+1:ny+ng) )
   allocate( kh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
   allocate( piinit3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
   allocate( sfctke(-ng+1:nx+ng,-ng+1:ny+ng) )


   allocate ( ji10(-ng+1:nx+ng,-ng+1:ny+ng) )
   allocate ( ji11(-ng+1:nx+ng,-ng+1:ny+ng) )
   allocate ( ji6(-ng+1:nx+ng,-ng+1:ny+ng) )

      do i=1,nx
      do j=1,ny
    
!      glat(i,j)=35.00
!      glon(i,j)=-100.00
      zsfc(i,j)=zc(1)

!.....msb 5/15/08 change for stretched grid
!      g12(i,j)=1.017     !g12(j,i)=ztop/(ztop-zsfc(j,i))
!      g12(i,j)=zcdx(1,nz-1)/(zcdx(1,nz-1)-zcdx(1,1))
      g12(i,j)=zc(nz-1)/(zc(nz-1)-zc(1))
       
      do k=1,nz

      piinit3d(i,j,k)=piinit%flt1d(k)
      kh(i,j,k)=km%flt3d(i,j,k)*3.0
      sgz(k)=zc(k)
      wgz(k)=zc(k)


!       print*,sgz(k),zc(k)

  
      end do

!      print*,tcanp%flt2d(i,j),wcanp%flt2d(i,j),i,j

      end do
      end do

! 2D advection, if needed
       IF ( ( ugrid /= 0.0 .or. vgrid /= 0.0 ) .and. isfcphys <= 3) THEN

         i = 0
         DO n = tcanp,radlw
           IF ( gd%var(n)%dyntype == 1 ) THEN
            i = i + 1
            m(i) = n
           ENDIF
         ENDDO
         
         num2d = i+2
         
! MPI: Send 2D data for advection
#ifdef MPI
    IF ( number_of_processes .gt. 1 ) THEN
    
    CALL cld_cpu('MPI-COMM-SFC')
         
         allocate( tmp3d(-ng+1:nx+ng,-ng+1:ny+ng,num2d) )
         
         DO n = 1,num2d-2
            tmp3d(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp,n) = gd%var(m(n))%flt2d(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp)
         ENDDO

         tmp3d(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp,num2d-1) = precip(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp,3) ! RAIN_ACC
         tmp3d(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp,num2d  ) = precip(-ngtmp+1:nx+ngtmp,-ngtmp+1:ny+ngtmp,4) ! HAIL_ACC
        
        nb = ngtmp ! number of ghost zones needed to comm. 3 zones for 6th-order 2D Crowley advection
        
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,num2d,ng,ng,0,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,tmp3d)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,num2d,ng,ng,0,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,tmp3d)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,num2d,ng,ng,0,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,tmp3d)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,num2d,ng,ng,0,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,tmp3d)

        ENDIF

         DO n = 1,num2d-2
            gd%var(m(n))%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = tmp3d(-ng+1:nx+ng,-ng+1:ny+ng,n) 
         ENDDO
         precip(-ng+1:nx+ng,-ng+1:ny+ng,3) = tmp3d(-ng+1:nx+ng,-ng+1:ny+ng,num2d-1) ! RAIN_ACC
         precip(-ng+1:nx+ng,-ng+1:ny+ng,4) = tmp3d(-ng+1:nx+ng,-ng+1:ny+ng,num2d  ) ! HAIL_ACC

         deallocate( tmp3d )

    CALL cld_cpu('MPI-COMM-SFC')
    
    ENDIF

#endif

! 2D advection for grid motion
         DO n = 1,num2d-2
           i = m(n)
           call ADVECTCRW2D         &
                (nx,ny,nst,             &
                 dt,dx,dy, gtx,gty,   &
                 gd%var(i)%flt2d, ugrid,vgrid,           &
                 ji10) 
         ENDDO
 
         
         DO n = 3,4
           call ADVECTCRW2D          &
                (nx,ny,nst,             &
                 dt,dx,dy, gtx,gty,   &
                 precip(-ng+1,-ng+1,n), ugrid,vgrid,           &
                 ji10) 
         ENDDO
       
       ENDIF ! ( ugrid /= 0.0 .or. vgrid /= 0.0 )


       ji10(:,:) = lat
       ji11(:,:) = lon



      ji6(-ng+1:nx+ng,-ng+1:ny+ng) = precip(-ng+1:nx+ng,-ng+1:ny+ng,1)

!      subroutine sfcphy(an,u,fu,v,fv,t,ft,qv,fqv,p,km,kh,piinit,  & 
!           tcanp, wcanp, qav, eflx, fflx, precip, veg,         &
!           g12,sgz,wgz,glat,glon,                              &
!           zsfc,stype,tsfc,wsfc,tsoil,wsoil,vlai,              &
!           uflx,vflx,tflx,qflx,radsw,radlw,albedo,rough,       &
!           sfctke,time,nx,ny,nz,ns,dt)
                  
!      write(91,*) 'call sfcphy'
      IF ( isfcphys <= 3 ) THEN
      call sfcphy(an,utmp,mfu,vtmp,mfv,s(1)%flt3d,mft,s(2)%flt3d, &
           mfq,pi%flt3d,km%flt3d,kh,piinit3d,gd%var(tcanp)%flt2d,    &
              gd%var(wcanp)%flt2d,gd%var(qav)%flt2d,                 &
           gd%var(eflx)%flt2d,gd%var(fflx)%flt2d,ji6,                &
                gd%var(veg)%flt2d,g12,sgz,wgz,ji10,ji11,             &
           zsfc,gd%var(stype)%flt2d,gd%var(tsrfc)%flt2d,             &
             gd%var(wsfc)%flt2d,gd%var(tsoil)%flt2d,                 &
             gd%var(wsoil)%flt2d,gd%var(vlai)%flt2d,gd%var(uflx)%flt2d,   &
           gd%var(vflx)%flt2d,gd%var(tflx)%flt2d,gd%var(qflx)%flt2d,gd%var(radsw)%flt2d,  &
              gd%var(radlw)%flt2d,gd%var(albedo)%flt2d,gd%var(rough)%flt2d,   &
           sfctke,time,nx,ny,nz,ns,dt)  


      u%flt3d(:,:,1)    = u%flt3d(:,:,1)    + dt * mfu(:,:)
      v%flt3d(:,:,1)    = v%flt3d(:,:,1)    + dt * mfv(:,:)
      s(1)%flt3d(:,:,1) = s(1)%flt3d(:,:,1) + dt * mft(:,:)
      s(2)%flt3d(:,:,1) = s(2)%flt3d(:,:,1) + dt * mfq(:,:)
      
      ELSEIF ( isfcphys == 4 .or. isfcphys == 5 ) THEN ! Grabowski (2006, QJ) heat flux (hard-wired values)
      
      call sfcheat(an,utmp,mfu,vtmp,mfv,s(1)%flt3d,mft,s(2)%flt3d, &
           mfq,pi%flt3d,km%flt3d,kh,piinit3d,gd%var(tcanp)%flt2d,    &
              gd%var(wcanp)%flt2d,gd%var(qav)%flt2d,                 &
           gd%var(eflx)%flt2d,gd%var(fflx)%flt2d,ji6,                &
                gd%var(veg)%flt2d,g12,sgz,wgz,ji10,ji11,             &
           zsfc,gd%var(stype)%flt2d,gd%var(tsrfc)%flt2d,             &
             gd%var(wsfc)%flt2d,gd%var(tsoil)%flt2d,                 &
             gd%var(wsoil)%flt2d,gd%var(vlai)%flt2d,gd%var(uflx)%flt2d,   &
           gd%var(vflx)%flt2d,gd%var(tflx)%flt2d,gd%var(qflx)%flt2d,gd%var(radsw)%flt2d,  &
              gd%var(radlw)%flt2d,gd%var(albedo)%flt2d,gd%var(rough)%flt2d,   &
           sfctke,time,nx,ny,nz,ns,dt)  

       DO k = 1,10
        ! zfac = exp(-(zc(k)-zc(1))/zlen) ! test 1
         zfac = exp(-zc(k)/zlen) ! test 2
         s(1)%flt3d(:,:,k) = s(1)%flt3d(:,:,k) + dt * mft(:,:)*zfac
         s(2)%flt3d(:,:,k) = s(2)%flt3d(:,:,k) + dt * mfq(:,:)*zfac
       ENDDO
      
      ELSE 
      
        write(0,*) 'invalid value of isfcphys ',isfcphys
      ENDIF
       

   
   deallocate ( ji6, ji10, ji11 )
   
   deallocate ( mfu,mfv,mft,mfq)
   
   deallocate ( g12,zsfc,kh,piinit3d,sfctke )


   
   RETURN
   END

!==================================================================================================================
!==================================================================================================================


      subroutine sfcheat(an,u,fu,v,fv,t,ft,qv,fqv,p,km,kh,piinit,  & 
           tcanp, wcanp, qav, eflx, fflx, precip, veg,         &
           g12,sgz,wgz,glat,glon,                              &
           zsfc,stype,tsfc,wsfc,tsoil,wsoil,vlai,              &
           uflx,vflx,tflx,qflx,radsw,radlw,albedo,rough,       &
           sfctke,time,nx,ny,nz,ns,dt)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      USE GRID_MODULE
      USE PARAM_MODULE
      USE SPHYS_MODULE
      USE COMMASMPI_MODULE
               
      implicit none
      
      integer nx,ny,nz,ns
      integer time
      
      real dt

      real  an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      real  u(-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  fu(-ng+1:nx+ng,-ng+1:ny+ng) 
      real  v(-ng+1:nx+ng,-ng+1:ny+ng,2)
      real  fv(-ng+1:nx+ng,-ng+1:ny+ng) 
      real  t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  ft(-ng+1:nx+ng,-ng+1:ny+ng)
      real  qv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  fqv(-ng+1:nx+ng,-ng+1:ny+ng)
      real  p(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  km(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real  kh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real  piinit(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

      real precip(-ng+1:nx+ng,-ng+1:ny+ng)
      real sfctke(-ng+1:nx+ng,-ng+1:ny+ng)
      
      real g12(-ng+1:nx+ng,-ng+1:ny+ng)

      real sgz(nz), wgz(nz), mfc(nz), mfe(nz)

      real glat(-ng+1:nx+ng,-ng+1:ny+ng)
      real glon(-ng+1:nx+ng,-ng+1:ny+ng)
      real zsfc(-ng+1:nx+ng,-ng+1:ny+ng)

      real stype(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsfc(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsfc(-ng+1:nx+ng,-ng+1:ny+ng)
      real tsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real tcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real wcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real qav(-ng+1:nx+ng,-ng+1:ny+ng)
      real veg(-ng+1:nx+ng,-ng+1:ny+ng)
      real vlai(-ng+1:nx+ng,-ng+1:ny+ng)

      real uflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real vflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real tflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real qflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real eflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real fflx(-ng+1:nx+ng,-ng+1:ny+ng)
      real radsw(-ng+1:nx+ng,-ng+1:ny+ng)
      real radlw(-ng+1:nx+ng,-ng+1:ny+ng)
      real albedo(-ng+1:nx+ng,-ng+1:ny+ng)
      real rough(-ng+1:nx+ng,-ng+1:ny+ng)

      real mbar, nsum, npsum
!      real lmax
      real  kmsfc, kmtop

      real tflx1, qflx1, uflx1, vflx1
      
      integer i,j
      
      real pi1, pres, thv, tv, rho, cpm
      real delz1, delz


      integer :: ixe,jye,ix,jy,i1,j1,iranseed
      real :: rndnum
      real  :: ffunc,fluxs,fluxl,thour
      real, allocatable :: ranarray(:,:),ranarraysmth(:,:)

!-----------------------------------------------------------------------
      logical debugsfc
!-----------------------------------------------------------------------
!     real, dimension(:,:), allocatable :: dum
!     allocate( dum(ny,nx) )
!-----------------------------------------------------------------------
      debugsfc = .false.
!-----------------------------------------------------------------------
      if(debugsfc) write(*,*) ' Enter sfcphy '


      ixe = itile
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend - jybeg

      eflx(:,:) = 0.
      fflx(:,:) = 0.


!       year=coards(1)
!       month=coards(2)
!       day=coards(3)
!       hour=coards(4)
!       minute=coards(5)
!       second=coards(6)
                                          
! --- tst is the number of time steps in 1 minute. 
! --- Changing tst will change the frequenct of the convective
!       initiation check.

!     tst =  dt / dt
!     ntst = INT(tst)
!     ntst = MAX(ntst,1)
!     iternum = INT( time / dt )
!     icheck=MOD(iternum, ntst)
!     if(icheck .EQ. 0) then

!       do j=1,ny
!       do i=1,nx
!       print*,tsfc(i,j),t(1,i,j),qv(1,i,j),i,j,time
!       end do
!       end do


!     write(0,"('--------------SFCPHY: before vegpbl --------------')")
!     do i=1,nx
!     do j=1,ny
!     do k=1,nz
!     IF( .not. ( u(k,i,j) .lt. 500. .and.  u(k,i,j) .gt. -500. ) ) then
!     print*, u(k,i,j), t(k,i,j), qv(k,i,j),k,i,j,time
!     end if
!     end do
!     end do
!     end do


!    endif

      allocate ( ranarray(1:nxend,1:nyend) )
      allocate ( ranarraysmth(nxend,1:nyend) )

!.....fill 2D arrays

      iranseed = Nint( time/dt )
  ranarray(:,:) = 0.0
  ranarraysmth(:,:) = 0.0
      DO j = 1,nyend-1
       DO i = 1,nxend-1
         ranarray(i,j) = 1.0 + 0.1*2.0*(rndnum(iranseed) - 0.5) ! values of  0.9 to 1.1
       ENDDO
      ENDDO

    DO j = 1,nyend-2,4
     DO i = 1,nxend-2,4
       DO  j1 = 0,3
       DO  i1 = 0,3
         ix = Min(nxend,i+i1)
         jy = Min(nyend,j+j1)
         ranarraysmth(ix,jy) = ranarray(i,j)
       ENDDO
       ENDDO
     ENDDO
    ENDDO

!      print*,'after vegpbl'

      ! prescribed fluxes from Grabowski 2006 (QJRMS) Appendix A
      IF ( isfcphys == 4 ) THEN
        thour = time/3600.
        ffunc = Max(0.0, Cos(0.5*pii*(5.25 - thour)/5.25) )
      ELSEIF ( isfcphys == 5 ) THEN
        IF ( time <= 3600. ) THEN
          thour = 1.0
        ELSE
          thour = 0.25
        ENDIF
          ffunc = Max(0.0, Cos(0.5*pii*(5.25 - thour)/5.25) )
      ENDIF
      fluxs = 270.0*ffunc**1.5
      fluxl = 554.0*ffunc**1.3
!-----------------------------------------------------------------------
      do j=1,jye
      do i=1,ixe
        pi1 = (piinit(i,j,1) + p(i,j,1))
        pres = (psl*100.) * pi1**cvr
        thv = t(i,j,1) * (1. + 0.61 * qv(i,j,1) )
        tv = thv * pi1
        
        rho = pres / (rd * tv )
        cpm = cp * (1. + 0.8 * qv(i,j,1) )

!        print*,rho,pres,rd,tv,pi1,thv,t(i,j,1),qv(i,j,1),i,j
!
! --- surface layer tendencies 
!
          delz1 = sgz(1) / g12(i,j)
          delz = g12(i,j) / (sgz(2)-sgz(1))
!         delz =       1. / (sgz(2)-sgz(1))    ! mix on computational grid

!           delz1=100.0/g12(i,j)
!           delz=g12(i,j)/100.0


!         /msb/ tflx(i,j) is potential temp flux
!
!         msb 3/21/08 multiply by cpm*rho to get units of w m-2 s-1
!                     tflx(i,j) is already in units w m-2 s-1

!          tflx1 = kh(i,j,2) * delz * (t(i,j,1)-t(i,j,2))
          tflx1 = kh(i,j,2) * delz * cpm * rho * (t(i,j,1)-t(i,j,2))
          tflx1 = 0.0
          tflx(i,j) = fluxs*ranarray(ixbeg-1+i,jybeg-1+j)
!          tflx1 = Max(50.0, kh(i,j,2)) * delz * cpm * rho * (t(i,j,1)-t(i,j,2))
          ft(i,j)  = ft(i,j) +                   &
                  (fluxs - tflx1)/(cpm*rho*delz1)


!          print*,'fluxes',tflx1,tflx(i,j),t(i,j,1),t(i,j,2),ft(i,j)

!..........msb 4/12/08 multiply qflx1, uflx1, vflx1 by rho to correct units

!          qflx1 = kh(i,j,2) * delz * (qv(i,j,1)-qv(i,j,2))
          qflx1 = kh(i,j,2) * delz * rho * (qv(i,j,1)-qv(i,j,2))
          qflx1 = 0.0
          qflx(i,j) = fluxl/2.5e6*ranarray(ixbeg-1+i,jybeg-1+j)
!         dum(i,j) = qflx1
          fqv(i,j) = fqv(i,j) +   &
                  ( fluxl/2.5e6 - qflx1 ) / (rho*delz1)


!.......test turning off sfc fluxes

!          ft(i,j)=0.0
!          fqv(i,j)=0.0



!        print*,qflx1,qflx(i,j),fqv(i,j,1),qv(i,j,1),qv(i,j,2)
!        print*,kh(i,j,2),delz,rho,cpm,i,j
        
 1015 format('sfcflx: ',7(f7.2,1x),i5)
 
!         if((mod(time,4).eq.0).and.i.eq.30.and.j.eq.25) then
!         IF ( my_rank == 0 .and. i == ixend-1 .and. j == jyend-1 ) THEN
!         write(91,1015) tflx(i,j),qflx(i,j)*Lav,     &
!        tsfc(i,j),tsoil(i,j),t(i,j,1),wsfc(i,j),wsoil(i,j),time
!         end if
 
!........msb 5/9/08  turn off momentum fluxes

         fu(i,j)=0.0
         fv(i,j)=0.0

      enddo ! j
!        fv(1,ny,i)  =  0.5 * (3.*fv(1,ny-1,i) - fv(1,ny-2,i))
!         fv(1,ny,i)=0.0
      enddo
!      do j=1,jye
!        fu(1,j,nx)  =  0.5 * (3.*fu(1,j,nx-1) - fu(1,j,nx-2))
!         fu(1,j,nx)=0.0
!      enddo

!........msb 5/9/08  turn off momentum fluxes

         fu(:,:)=0.0
         fv(:,:)=0.0

      deallocate( ranarray, ranarraysmth )
      if(debugsfc) write(*,*) ' past sfcheat  '
      return
      END



!2345678901234567890123456789012345678901234567890123456789012345678912
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!  SUBROUTINE  ADVECTCRW
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!  Issues:
!       need time step or time for xy switching: should we switch 
!         directions between EnKF steps?
!
! 12.2005: Conversion from SAM to SWM model (erm)
!
! Original code by Jerry Straka
!
      subroutine ADVECTCRW2D         &
       (nx,ny,nst,             &
        dt,dx,dy, gtx,gty,   &
        an, ugrid,vgrid,           &
        t0) 

      USE GRID_MODULE
      USE PARAM_MODULE
      USE COMMASMPI_MODULE

      implicit none


!  declare all variables
!
!  INTEGERS
!
!      integer    ng,nor
      integer    ip
      integer    ipass
      integer    mdnstp
      integer    n1
      integer    n2
      integer    np
      integer    npass
      integer    nx
      integer    ny
      parameter (np=6)
      integer    i,j,k,n
      
      integer  nst
!
!  REALS
!
      real  :: ugrid, vgrid
      
      real       qu, qd, qc, qr, qdel, qxs, qxp, qsr, qsl
      real       aqdel, qcurv, aqcurv
      
      real   ba1, ba2
!
      real       dt
      real       dx
      real       dy

      real       dtdx
      real       dtdy

      real       dtdx4
      real       dtdy4
      real       dtfac
      real       vnorme
      real       vnormw
      real       vnormn
      real       vnorms

      real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng)
!
!      real       cx(nxl,nzl,np),  
!      real  px(nx,np,np),py(ny,np,np),pz(nz,np,np)
!      real  px(0:nx-1,np,np),py(0:ny-1,np,np),pz(0:nz-1,np,np)
      real  px(-ng+1:nx+ng,np,np),py(-ng+1:ny+ng,np,np)
!      real       cy(nyl,nzl,np),  py(nyl,np,np)
!      real       cz(nxl,nzl,np),  pz(nzl,np,np)
!
!
!
      real       ab
!
!      real       dc(0:nz)

      real       an(-ng+1:nx+ng,-ng+1:ny+ng) ! st
!
!      real u(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real v(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real w(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

!      real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real fu, fv
      real u,v
      real fus,fvs

      real t0(-ng+1:nx+ng,-ng+1:ny+ng)
!      real t7(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real t8(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real t9(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
      real t0s,t1s,t2s,t3s,t4s,t5s,t6s
      
      integer imn,imx
      integer jmn,jmx
      integer kmn,kmx
      
      
! params from sam.param.h
!      integer icrwmp
      integer istag,jstag,kstag
      parameter ( istag = 1, jstag = 1, kstag = 1 )
!
!  Grid parameters:  dimensions (on=1)
!
      integer    id1, id2, id3, id4, id5
      parameter (id1=1)
      parameter (id2=id1*2)
      parameter (id3=id1*3)
      parameter (id4=id1*4)
      parameter (id5=id1*5)
      integer    jd1, jd2, jd3, jd4, jd5
      parameter (jd1=1)
      parameter (jd2=jd1*2)
      parameter (jd3=jd1*3)
      parameter (jd4=jd1*4)
      parameter (jd5=jd1*5)
      integer    kd1, kd2, kd3, kd4, kd5
      parameter (kd1=1)
      parameter (kd2=kd1*2)
      parameter (kd3=kd1*3)
      parameter (kd4=kd1*4)
      parameter (kd5=kd1*5)

!      logical debug
!      parameter ( debug = .true. ) 

      real binflo
      parameter ( binflo = 0.0 )
      
      real ainflotmp
      logical  relaxscalar

!      integer icrwmn
!      parameter ( icrwmn = 1 )

      integer   ibc,jbc

      logical DBG
      parameter ( DBG = .false. )

      real, parameter :: eps = 1.0e-25

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
#endif

      ab = 0.0
!      icrwmp = 2
      
!      IF ( DBG ) print*, 'nx,ny,nz,dt,ia = ',nx,ny,nz,dt,ia
!
!
!  CHECKS
!
!
!
!  check to see that work arrays wont overflow...
!
      if ( ny .gt. nx ) then
!      write(*,*) 'nyl is greater than nxl;  ',ny,nx
!      write(*,*) 'stop in sam.a3.cx.f -- this may be outdated'
!      stop
      end if
!
!
!
!  SET RUNTIME PARAMETER TO DETERMINE IMPLEMENTATION OF CROWLEY SCHEME
!
!    icrwmp = 1;   at step 1, do x-pass, y-pass
!                  at step 2, do y-pass, x-pass
!    icrwmp = 2;   at step 1, do x-pass, y-pass, y-pass, x-pass
!                  at step 2, do y-pass, x-pass, x-pass, y-pass
!                    each part uses only 1/2dt for the time step giving a
!                    cumulative of dt; e.g., two half passes for each direction
! 
!  set type of crowley implementation ( sam.param.h )
!
!      IF ( bcx .eq. 2 .or. bcy .eq. 2 ) THEN
!        write(0,*) 'STOP! Crowley scheme not set up yet for periodic BCs!!'
!        STOP
!      ENDIF

!      if ( icrwmp .eq. 1 ) then
       dtfac = 1.0
       npass = 2
       mdnstp = mod(nst,2)
!      end if
!      if ( icrwmp .eq. 2 ) then
!       dtfac = 0.5
!       npass = 4
!       mdnstp = mod(nst,2)
!      end if
      
      
      ainflotmp = 0.0
      
      IF (  ainflo .gt. 0.0 ) THEN

         relaxscalar = .true.
  
   IF ( relaxscalar ) THEN
     ainflotmp = ainflo
   ENDIF
  
  ENDIF
!
!
!
!  CONSTANTS
!
!
!
!      dxi(ix)   = (1.0)/dx
!      dyi   = (1.0)/dy
!      dzi   = (1.0)/dz
!      dxt   = -dt/dx
!      dyt   = -dt/dy
!      dzt   = -dt/dz
!      dxt4  = -dt/(4.0*dx)
!      dyt4  = -dt/(4.0*dy)
!      dzt4  = -dt/(4.0*dz)

      dtdx  = -dtfac*dt/dx
      dtdy  = -dtfac*dt/dy
      dtdx4 =  (0.25)*dtfac*dt/dx
      dtdy4 =  (0.25)*dtfac*dt/dy
      
      u = -ugrid
      v = -vgrid
      
#ifdef MPI
        imn = 1
        imx = nxend-istag
        jmn = 1
        jmx = nyend-jstag
!        kmn = 1
!        kmx = nzend-kstag
#else
        imn = 1
        imx = nx-istag
        jmn = 1
        jmx = ny-jstag
!        kmn = 1
!        kmx = nz-kstag
#endif
       
       
       IF ( .false. .and. nst .le. 2 ) THEN
        do i=1,nx
         print*,'ix,gtx = ',i,gtx(i)
        ENDDO
        do j=1,ny
         print*,'jy,gty = ',j,gty(j)
        ENDDO
        print*, 'dtdx, etc: ',dtdx,dx,dy
        STOP
       ENDIF


!  IF( bcy .eq. 2 ) THEN    ! periodic
!
!   DO n = 1,ng
!    an(1:nx-1,ny-1+n,1:nz-1,ia) = an(1:nx-1,n   ,1:nz-1,ia)
!    an(1:nx-1,1-n   ,1:nz-1,ia) = an(1:nx-1,ny-n,1:nz-1,ia)
!    ad(1:nx-1,ny-1+n,1:nz-1) = ad(1:nx-1,n   ,1:nz-1)
!    ad(1:nx-1,1-n   ,1:nz-1) = ad(1:nx-1,ny-n,1:nz-1)
!   ENDDO
!   
!  ENDIF


!      IF ( DBG .and. ia .eq. 2 ) print*, 'imn,imx,jmn,jmx= ',imn,imx,jmn,jmx
      
!      print*, 'nst, ia, imn,imx,kmn,kmx= ',nst,ia,imn,imx,kmn,kmx
       
!       IF ( debug ) print*,'ia,mdnstp',ia,mdnstp
!
!
!  CREATE COEFFICIENTS FOR CROWELY FLUX SCHEME
!
!
!
!  create x-direction polynomial coefficients
!
!      if ( nx .gt. 2 ) then
!      call crwcff(np,nx,nx,ng,istag,id1,itile,ixbeg,ixend,nxbeg,nxend,px)
!      end if
!
!  create y-direction polynomial coefficients
!
!      if ( ny .gt. 2 ) then
!      call crwcff(np,ny,ny,ng,jstag,jd1,jtile,jybeg,jyend,nybeg,nyend,py)
!      end if

!  create x-direction polynomial coefficients
!
      ibc = 0
      jbc = 0
      if ( nx .gt. 2 ) then
       IF ( bcx .ne. 2 ) THEN
!        write(0,*) 'call crwcff for x'
        call crwcff(np,nx,nx,ng,istag,id1,itile,ixbeg,ixend,nxbeg,nxend,px)
       ELSE  !  periodic BC:
        call crwcff1(np,nx,nx,ng,istag,id1,px)
        ibc = id1 ! + istag
#ifdef MPI
        IF ( ixbeg .eq. nxbeg .and. imn .eq. 1    )    imx = nxend - 1
        IF ( ixend .eq. nxend .and. imx .eq. nxend-1 ) imn = 1
#else
        IF ( imn .eq. 1    ) imx = nx - 1
        IF ( imx .eq. nx-1 ) imn = 1
#endif
       ENDIF
      end if
!
!  create y-direction polynomial coefficients
!
      if ( ny .gt. 2 ) then
       IF ( bcy .ne. 2 ) THEN
!        write(0,*) 'call crwcff for y'
        call crwcff(np,ny,ny,ng,jstag,jd1,jtile,jybeg,jyend,nybeg,nyend,py)
       ELSE
        call crwcff1(np,ny,ny,ng,jstag,jd1,py)
        jbc = jd1 ! + jstag
#ifdef MPI
        IF ( jybeg .eq. nybeg .and. jmn .eq. 1    )    jmx = nyend - 1
        IF ( jyend .eq. nyend .and. jmx .eq. nyend-1 ) jmn = 1
#else
        IF ( jmn .eq. 1    ) jmx = ny - 1
        IF ( jmx .eq. ny-1 ) jmn = 1
#endif
       ENDIF
      end if
!
!
!  create z-direction polynomial coefficients
!
!
!
!  SIXTH ORDER CROWLEY FLUX-ANTI FLUX SCHEME
!
!
!
      do ipass = 1,npass
!
!  x-direction crowley corrected flux computations
!
      if ( ipass .eq. (1+mdnstp) .or. ipass .eq. (4-mdnstp) ) then 
#ifdef MPI
      if ( nxend .gt. 5 ) then
#else
      if ( nx .gt. 5 ) then
#endif
!
!
! C$DOACROSS LOCAL(k,j,i)
!!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
!      do 2011 k = 1,     nz-kstag
!      do 2012 j = 1,     ny-jstag
!      do 2013 i = 1,     nx-istag-id1
! 2013 continue
! 2012 continue
! 2011 continue

!      if ( icrwmn .eq. 1 ) then
!      do  k = 1,     nz-kstag
!      do  j = 1,     ny-jstag
!      do  i = 0,     ng-1 
!      an(-i,j,k,ia)   = an(i+1,j,k,ia)
!      an(nx+i,j,k,ia) = an(nx-istag-i,j,k,ia)
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDIF

      IF ( bcx .eq. 2 ) THEN  ! periodic
#ifdef MPI
    IF ( nproci == 1 ) THEN
       ixb=1
       ixe=itile
       if(ixend.eq.nxend) ixe=ixend-ixbeg
       
       jyb=1
       jye=jtile
       if(jyend.eq.nyend) jye=jyend-jybeg

       
       do n=1,ng
        do j=jyb,jye ; do i=ixb,ixe
        an( nxend-1+n,j) = an(n   ,j)
        an( 1-n   ,j) = an(nxend-n,j)
        enddo ; enddo
       enddo
     ENDIF
#else
      DO n = 1,ng
        an( nx-1+n,1:ny) = an(n   ,1:ny)
        an( 1-n   ,1:ny) = an(nx-n,1:ny)
       ENDDO
#endif
      ENDIF

      
      t0(:,:) = 0.0
!
! C$DOACROSS LOCAL(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,k,j,i,ia), 
! C$&       share(t1,t2,t3,t4,t5,t6,an,ad,fu)
! !$omp  PARALLEL DO DEFAULT(SHARED),  &
! !$omp  PRIVATE(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,k,j,i,ia) &
! !$omp shared(t1,t2,t3,t4,t5,t6,an,ad,fu)
!  do 2020 ia = 1,     iadva
#ifdef MPI
  ixb = 0
  ixe = itile
  if (ixbeg .le. imn) ixb = imn-ibc
  if (ixend .ge. imx) ixe = imx-ixbeg+1-id1+ibc

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
  if (jybeg .le. jmn) jyb = jmn
  if (jyend .ge. jmx) jye = jmx-jybeg+1


   do j = jyb,jye ; do i = ixb,ixe
  
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i,  &
!$OMP   fus,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)  
      do j = jmn,jmx ! 1,     ny-jstag
!      do 2023 i = Max(1,imn-1),Min(nx-istag-id1,imx)  ! 1,     nx-istag-id1
      do i = imn-ibc, imx-id1+ibc ! 1,     nx-istag-id1
#endif
      fus =                             &
        dtdx4                                &
       *(2.0*u)      &
       /(gtx(Max(1,i)) + gtx(i+id1))

      IF ( Abs(fus) .lt. eps ) THEN
       t0(i,j) = 0.0
       CYCLE
      ENDIF

      t1s =    &
        fus*(px(i,1,1) +   &
        fus*(px(i,1,2) +   &
        fus*(px(i,1,3) +   &
        fus*(px(i,1,4) +   &
        fus*(px(i,1,5) +   &
        fus*(px(i,1,6)))))))
      t2s =    &
        fus*(px(i,2,1) +   &
        fus*(px(i,2,2) +   &
        fus*(px(i,2,3) +   &
        fus*(px(i,2,4) +   &
        fus*(px(i,2,5) +   &
        fus*(px(i,2,6)))))))
      t3s =    &
        fus*(px(i,3,1) +   &
        fus*(px(i,3,2) +   &
        fus*(px(i,3,3) +   &
        fus*(px(i,3,4) +   &
        fus*(px(i,3,5) +   &
        fus*(px(i,3,6)))))))
      t4s =    &
        fus*(px(i,4,1) +   &
        fus*(px(i,4,2) +   &
        fus*(px(i,4,3) +   &
        fus*(px(i,4,4) +   &
        fus*(px(i,4,5) +   &
        fus*(px(i,4,6)))))))
      t5s =    &
        fus*(px(i,5,1) +   &
        fus*(px(i,5,2) +   &
        fus*(px(i,5,3) +   &
        fus*(px(i,5,4) +   &
        fus*(px(i,5,5) +   &
        fus*(px(i,5,6)))))))
      t6s =    &
        fus*(px(i,6,1) +   &
        fus*(px(i,6,2) +   &
        fus*(px(i,6,3) +   &
        fus*(px(i,6,4) +   &
        fus*(px(i,6,5) +   &
        fus*(px(i,6,6)))))))
      t0s =    &
        t1s*an(i-id2,j)+   &
        t2s*an(i-id1,j)+   &
        t3s*an(i    ,j)+   &
        t4s*an(i+id1,j)+   &
        t5s*an(i+id2,j)+   &
        t6s*an(i+id3,j)
      IF ( icrwmn .eq. 1  ) THEN
      t0s = -t0s / (fus)
      qsr  = max(sign((1.0),fus), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(i-id1*1,j)   &
           + qsl*an(i+id1*2,j)
      qd   = qsr*an(i+id1*1,j)   &
           + qsl*an(i,      j)
      qc   = qsr*an(i,      j)   &
           + qsl*an(i+id1*1,j)
      qr   = qu + (qc - qu) / (abs(fus))
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),   &
                 abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
!      IF ( DBG ) print*, 'i,j,k,imn,imx,ia = ',i,j,k,imn,imx,ia
      t0(i,j) =   &
        -fus *    &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
       t0(i,j) = t0s
      ENDIF
      
      enddo
      enddo

!       

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = max(2,imn)
      if (ixend .eq. nxend) ixe = min(ixend-ixbeg+1-istag-id1,imx-ixbeg+1)

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
      if (jybeg .eq. nybeg) jyb = jmn
      if (jyend .eq. nyend) jye = jmx-jybeg+1

!      kzb = 1
!      kze = ktile
!      if (kzbeg .lt. kmn) kzb = kmn
!      if (kzend .gt. kmx) kze = kmx-kzbeg+1

      do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do  j = jmn,jmx ! 1,     ny-jstag
      do i = imn+id1-ibc, imx-id1+ibc ! Max(2,imn), Min(nx-istag-id1,imx) ! 1+id1, nx-id1-istag
#endif
      an(i,j) = an(i,j)    &
        + (t0(i,j)-t0(i-id1,j))   &
        * (gtx(i)**2)
      ENDDO
      ENDDO
      
!
! 2020 continue
!
!  x-direction boundaries
!

      IF ( .FALSE. ) THEN
#ifdef MPI

      IF ( bcx .ne. 2 .and. ( ixbeg .eq. nxbeg .or. ixend .eq. nxend ) ) THEN

       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag


       do j = jyb,jye
#else

      IF ( bcx .ne. 2 .and. ( imn .eq. 1 .or. imx .eq. nx-istag ) ) THEN

! !$omp  PARALLEL DO DEFAULT(SHARED), &
! !$omp  PRIVATE(vnormw,vnorme,k,j,ba1,ba2,ia),SHARED(vc,an,ad) 
      do j = 1,ny-jstag

#endif
        ba1 = ab
        ba2 = ab
#ifdef MPI
       IF ( ixbeg .eq. nxbeg ) THEN
#else
       IF ( imn .eq. 1 ) THEN
#endif
       vnormw = u
       an(1,j) = an(1,j)    &
        + (ainflotmp*min(vnormw*dtdx, 0.0) - dt*dtfac*binflo)   &
          *(an(1,j) - ba1)
       ENDIF
#ifdef MPI
       IF ( ixend .eq. nxend ) THEN
#else
       IF ( imx .eq. nx-istag ) THEN
#endif
       vnorme = u
       an(nx-istag,j) = an(nx-istag,j)    &
        - (ainflotmp*max(vnorme*dtdx, 0.0) + dt*dtfac*binflo)   &
          *(an(nx-istag,j) - ba2)
       ENDIF

      enddo
      ENDIF !( bcy .ne. 2 .and. ( imn .eq. 1 .or. imx .eq. nx-istag ) )
      ENDIF
      end if ! nx.gt.5
!
      end if ! ipass
!
!  y-direction crowley corrected flux computation
!
      if ( ipass .eq. (2-mdnstp) .or. ipass .eq. (3+mdnstp) ) then
!
      if ( ny .gt. 5 ) then
!
! C$DOACROSS LOCAL(k,j,i)
!
! C$DOACROSS LOCAL(k,j,i)
! !$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i) 

!
  IF( bcy .eq. 2 ) THEN    ! periodic
#ifdef MPI
    IF ( nprocj == 1 ) THEN
    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg


     do n = 1,ng
     do j = jyb,jye ; do i = ixb,ixe
      an(i,nyend-1+n) = an(i,n      )
      an(i,1-n      ) = an(i,nyend-n)
     enddo ; enddo
     enddo
    ENDIF
#else
   DO n = 1,ng
    an(1:nx-1,ny-1+n) = an(1:nx-1,n   )
    an(1:nx-1,1-n   ) = an(1:nx-1,ny-n)
   ENDDO
#endif
  ENDIF
!      if ( icrwmn .eq. 1 ) then 
!      do kz = 1,     nz-kstag
!      do jy = 0,     ng-1
!      do ix = 1,     nx-istag
!      an(ix,-jy,kz,ia)   = an(ix,jy+1,kz,ia)
!      an(ix,ny+jy,kz,ia) = an(ix,ny-jstag-jy,kz,ia)
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDIF
!

      t0(:,:) = 0.0

! C$DOACROSS LOCAL(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,k,j,i,ia), 
! C$&       share(t1,t2,t3,t4,t5,t6,an,ad,fv)
! !$omp PARALLEL DO DEFAULT(SHARED), &
! !$omp PRIVATE(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,k,j,i,ia), &
! !$omp SHARED(t1,t2,t3,t4,t5,t6,an,ad,fv)  
!      do 3020 ia = 1,     iadva

#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = imn
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 0
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn-jbc
      if (jyend .ge. jmx) jye = jmx-jybeg+1-jd1+jbc


       do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(k,j,i,     &
!$OMP   fvs,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)
      do j = jmn-jbc, jmx-jd1+jbc ! Max(1,jmn-1),Min(ny-jstag-jd1,jmx)  ! 1, ny-jstag-jd1
      do i = imn,imx ! 1, nx-istag
#endif
      fvs =                             &
        dtdy4                                &
       *(2.0*v)      &
       /(gty(Max(1,j)) + gty(j+jd1))

      IF ( Abs(fvs) .lt. eps ) THEN
       t0(i,j) = 0.0
       CYCLE
      ENDIF

      t1s =     &
        fvs*(py(j,1,1) +   &
        fvs*(py(j,1,2) +   &
        fvs*(py(j,1,3) +   &
        fvs*(py(j,1,4) +   &
        fvs*(py(j,1,5) +   &
        fvs*(py(j,1,6)))))))
      t2s =    &
        fvs*(py(j,2,1) +   &
        fvs*(py(j,2,2) +   &
        fvs*(py(j,2,3) +   &
        fvs*(py(j,2,4) +   &
        fvs*(py(j,2,5) +   &
        fvs*(py(j,2,6)))))))
      t3s =    &
        fvs*(py(j,3,1) +   &
        fvs*(py(j,3,2) +   &
        fvs*(py(j,3,3) +   &
        fvs*(py(j,3,4) +   &
        fvs*(py(j,3,5) +   &
        fvs*(py(j,3,6)))))))
      t4s =    &
        fvs*(py(j,4,1) +   &
        fvs*(py(j,4,2) +   &
        fvs*(py(j,4,3) +   &
        fvs*(py(j,4,4) +   &
        fvs*(py(j,4,5) +   &
        fvs*(py(j,4,6)))))))
      t5s =    &
        fvs*(py(j,5,1) +   &
        fvs*(py(j,5,2) +   &
        fvs*(py(j,5,3) +   &
        fvs*(py(j,5,4) +   &
        fvs*(py(j,5,5) +   &
        fvs*(py(j,5,6)))))))
      t6s =    &
        fvs*(py(j,6,1) +   &
        fvs*(py(j,6,2) +   &
        fvs*(py(j,6,3) +   &
        fvs*(py(j,6,4) +   &
        fvs*(py(j,6,5) +   &
        fvs*(py(j,6,6)))))))
      t0s =    &
        t1s*an(i,j-jd2) +   &
        t2s*an(i,j-jd1) +   &
        t3s*an(i,j    ) +   &
        t4s*an(i,j+jd1) +   &
        t5s*an(i,j+jd2) +   &
        t6s*an(i,j+jd3)
      IF ( icrwmn .eq. 1 ) THEN
      t0s = -t0s / (fvs)
      qsr  = max(sign((1.0),fvs), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(i,j-jd1*1)   &
           + qsl*an(i,j+jd1*2)
      qd   = qsr*an(i,j+jd1*1)   &
           + qsl*an(i,j      )
      qc   = qsr*an(i,j      )   &
           + qsl*an(i,j+jd1*1)
      qr   = qu + (qc - qu) / (abs(fvs))
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),   &
                 abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
      t0(i,j) =   &
       -fvs *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
        t0(i,j) = t0s
      ENDIF

     ENDDO
     ENDDO

#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = imn
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = max(2,jmn)
      if (jyend .eq. nyend) jye = min(nyend-jybeg-jstag-jd1,jmx-jybeg+1)


      do j = jyb,jye ; do i = ixb,ixe
#else
      do  j = Max(2,jmn), Min(ny-jstag-jd1,jmx) ! 1+jd1, ny-jstag-jd1
      do  i = imn,imx                           ! 1, nx-istag
#endif
       an(i,j) = an(i,j)      &
         + (t0(i,j)-t0(i,j-jd1))    &
         * (gty(j)**2)
      ENDDO
      ENDDO
      

!
! 3020 continue
!
!  y-direction boundaries
!  scalars on south and north boundaries
!
     IF ( .FALSE. ) THEN
#ifdef MPI
      IF ( bcy .ne. 2 .and. ( jybeg .eq. nybeg .or. jyend .eq. nyend ) ) THEN
#else
      IF ( bcy .ne. 2 .and. ( jmn .eq. 1 .or. jmx .eq. ny-jstag) ) THEN
#endif
! C$DOACROSS LOCAL(vnorms,vnormn,kz,ix,ia), 
! C$&       share(vc,an,ad)
! !$omp PARALLEL DO DEFAULT(SHARED), &
! !$omp PRIVATE(vnorms,vnormn,kz,ix,ia,ba1,ba2), SHARED(vc,an,ad)  
!      do 3051 ia = 1,     iadva
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      do i = ixb,ixe
#else
      do i = 1, nx-istag
#endif
        ba1 = ab
        ba2 = ab

#ifdef MPI
      IF ( jybeg .eq. nybeg ) THEN
#else
      IF ( jmn .eq. 1 ) THEN
#endif
      vnorms = v
      an(i,1) = an(i,1)    &
       + (ainflotmp*min(vnorms*dtdy, 0.0) - dt*dtfac*binflo)   &
         *(an(i,1) - ba1)
      ENDIF

#ifdef MPI
      IF ( jyend .eq. nyend ) THEN
#else
      IF ( jmx .eq. ny-jstag ) THEN
#endif
      vnormn = v
      an(i,ny-jstag) = an(i,ny-jstag)    &
       - (ainflotmp*max(vnormn*dtdy, 0.0) + dt*dtfac*binflo)   &
         *(an(i,ny-jstag) - ba2)
      ENDIF
     enddo


      ENDIF  !( bcy .ne. 2 .and. ( jmn .eq. 1 .or. jmx .eq. nyend-jstag) )
      ENDIF ! false
!     enddo
      end if
!
      end if
      
!
!
!  end of crowley flux scheme for scalars
!
!
   end do
!
!
!
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!
!
      return
      end
!

