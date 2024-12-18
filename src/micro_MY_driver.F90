#ifndef RKIND
#define RKIND 4
#endif

!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE MY_drive                        ######
!######                                                      ######
!##################################################################
!##################################################################

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  The driver for the multimoment microphysics scheme developed by
!  Jason Milbrandt (McGill University / Meteorological Services of Canada).
!
!  REFERENCES:
!
!  Milbrandt, J. A. and M. K. Yau, 2005: A multi-moment bulk microphysics
!  parameterization. Part I: Analysis of the role of the spectral shape
!  parameter. J. Atmos. Sci., 62, 3051-3064.
!
!  Milbrandt, J. A. and M. K. Yau, 2005: A multi-moment bulk microphysics
!  parameterization. Part II: A proposed three-moment closure and scheme
!  description. J. Atmos. Sci., 62, 3065-3081.
!
!-----------------------------------------------------------------------
!
!   AUTHOR: Dan Dawson
!   08/02/2011
!
!   MODIFICATION HISTORY:
!     Dan Dawson 08/02/2011
!       First written.  Based on the similar micro_MY.f90 driver in ARPS 5.3
!       and adapted for the COMMAS model.
!
!-----------------------------------------------------------------------

SUBROUTINE MY_DRIVE(gd,mscheme,nscalar,            &
                      pibase,                               & ! base state pi
                      pim3d,                                & ! pi_prime at t
                      sm4d,                                 & ! scalars at t
                      pi3d,                                 & ! pi_prime at t+1
                      s4d,                                  & ! scalars at t+1
                      xfalltot,                             & ! precip
                      dt,z1d,                               & ! time step, height
                      nx,ny,nz,nor,norz,na,iunit,           &
                      w,                                    & ! vertical velocity
                      t0,t1,                                & ! temporary arrays
                      dbz,io_flag)                            ! reflectivity

  !USE FILE_MODULE
  !USE VIS5D_MODULE
  USE GRID_MODULE
  USE MICRO_MODULE, only:  ntc_my,         &
                           cnor_my,        &
                           cnos_my,        &
                           cnog_my,        &
                           cnoh_my,        &
                           rho_qr_my,      &
                           rho_qi_my,      &
                           rho_qs_my,      &
                           rho_qg_my,      &
                           rho_qh_my,      &
                           alphar_my,      &
                           alphai_my,      &
                           alphas_my,      &
                           alphag_my,      &
                           alphah_my,      &
                           ndebug

                            

  USE INDEX_MODULE
  USE CPUTIME_MODULE
  USE COMMASMPI_MODULE, only: itile,jtile,ktile, &
                              ixbeg,ixend,nxend, &
                              jybeg,jyend,nyend, &
                              kzbeg,kzend,nzend  
  USE MY_TMOM_MOD
  
  implicit none

  integer nx,ny,nz,na,nor,norz

  TYPE(GRID)         :: gd 

  INTEGER            :: nscalar

  integer, parameter    :: ng1 = 1

  real               :: w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

!
! external temporary arrays
!
  real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
  real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

  real pibase(-norz+ng1:nz+norz)  ! base state Pi
  real pim3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi (time t)
  real pi3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi (time t+1)
  real sm4d(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
  real  s4d(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
  real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

  real xfalltot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,5)

  real dt

  real z1d(-norz+ng1:nz+norz,4)

  integer mscheme

  integer i,j,k,n
  integer iunit
  
  logical io_flag

  integer, parameter :: istag = 1, jstag = 1, kstag = 1
  integer :: ixe, jye, kze
  
  !
  ! local arrays, Memory order in i-k-j
  !
  REAL     :: tk (nx,nz,ny)   ! temperature at time tfuture
  REAL     :: tkm(nx,nz,ny)   ! temperature at time tpast

  REAL     :: p (nx,nz,ny)    ! total pressure at time tfuture
  REAL     :: pm(nx,nz,ny)    ! total pressure at time tpast

  REAL     :: q (nx,nz,ny)    ! total mixing ratio at time tfuture
  REAL     :: qm(nx,nz,ny)    ! total mixing ratio at time tpast

  REAL     :: ww  (nx,nz)
  REAL     :: zp2d(nx,nz)

  REAL     :: lpr(nx,ny)       ! liquid precipitation rate
  REAL     :: spr(nx,ny)       ! soild precipitation rate

  REAL     :: qnz (nx,nz,nscalar) ! hydrometeor arrays at time tfuture
  REAL     :: qnzm(nx,nz,nscalar) ! at time tpresent, qc, qr, qi, qs, qg, qh
                                  !                   nc, nr, ni, ns, ng, nh
                                  !                       zr, zi, zs, zg, zh
  REAL     :: tem(nx,nz,12)       ! used only when scheme == 1
  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc Begin Execute
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  
  CALL cld_cpu('MICROPHYSICS')
  !WRITE(0,*) 'MY_DRIVE: Inside MY_DRIVE'
  !write(iunit,*) 'rho_qi,rho_qs,rho_qg,rho_qh = ',rho_qi_my,rho_qs_my,rho_qg_my,rho_qh_my
  !write(iunit,*) 'lc,lr,li,ls,lh,lhl,lnc,lnr,lni,lns,lnh,lnhl,lzr,lzi,lzs,lzh,lzhl'
  !write(iunit,*) lc,lr,li,ls,lh,lhl,lnc,lnr,lni,lns,lnh,lnhl,lzr,lzi,lzs,lzh,lzhl
  
  !-----------------------------------------------------------------------
  !
  ! Compute total temperature and total pressure
  !
  !-----------------------------------------------------------------------
  
  !write(0,*) 'nx,ny,nz',nx,ny,nz

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag
      
      kze = nz-1 ! assume no vertical tiling for now. Will fail otherwise.
  
  DO k=1,kze
    DO j=1,jye
      DO i=1,ixe

        p  (i,k,j) =  1.0e5*(pibase(k)+pi3d(i,j,k))**3.509
        tk (i,k,j) =  s4d(i,j,k,lt)*(pibase(k)+pi3d(i,j,k))
        q  (i,k,j) =  s4d(i,j,k,lv)

        pm (i,k,j) =  1.0e5*(pibase(k)+pim3d(i,j,k))**3.509
        tkm(i,k,j) =  sm4d(i,j,k,lt)*(pibase(k)+pim3d(i,j,k))
        qm (i,k,j) =  sm4d(i,j,k,lv)
        
        !if(i == 50 .and. j == 50) THEN
        !  WRITE(0,*) 'p,tk,q',p(i,k,j),tk(i,k,j),q(i,k,j)
        !END IF

      END DO
    END DO
  END DO

  DO j = 1,jye
    DO k = 1,kze
      DO i = 1,ixe

        ww  (i,k) = w(i,j,k)
        zp2d(i,k) = z1d(k,2)
                
        IF(mscheme >= 1) THEN 
          qnz(i,k,1) = s4d(i,j,k,lc) ! QC
          qnz(i,k,2) = s4d(i,j,k,lr) ! QR
          qnz(i,k,3) = s4d(i,j,k,li) ! QI
          qnz(i,k,4) = s4d(i,j,k,ls) ! QS
          qnz(i,k,5) = s4d(i,j,k,lh) ! QG
          qnz(i,k,6) = s4d(i,j,k,lhl) ! QH
          
          qnzm(i,k,1) = sm4d(i,j,k,lc) ! QC
          qnzm(i,k,2) = sm4d(i,j,k,lr) ! QR
          qnzm(i,k,3) = sm4d(i,j,k,li) ! QI
          qnzm(i,k,4) = sm4d(i,j,k,ls) ! QS
          qnzm(i,k,5) = sm4d(i,j,k,lh) ! QG
          qnzm(i,k,6) = sm4d(i,j,k,lhl) ! QH
        END IF
        
        IF(mscheme >= 2) THEN
          qnz(i,k,7) = s4d(i,j,k,lnc) ! NC
          qnz(i,k,8) = s4d(i,j,k,lnr) ! NR
          qnz(i,k,9) = s4d(i,j,k,lni) ! NI
          qnz(i,k,10) = s4d(i,j,k,lns) ! NS
          qnz(i,k,11) = s4d(i,j,k,lnh) ! NG
          qnz(i,k,12) = s4d(i,j,k,lnhl) ! NH
          
          qnzm(i,k,7) = sm4d(i,j,k,lnc) ! NC
          qnzm(i,k,8) = sm4d(i,j,k,lnr) ! NR
          qnzm(i,k,9) = sm4d(i,j,k,lni) ! NI
          qnzm(i,k,10) = sm4d(i,j,k,lns) ! NS
          qnzm(i,k,11) = sm4d(i,j,k,lnh) ! NG
          qnzm(i,k,12) = sm4d(i,j,k,lnhl) ! NH
        END IF
        
        IF(mscheme == 4) THEN
          qnz(i,k,13) = s4d(i,j,k,lzr) ! ZR
          qnz(i,k,14) = s4d(i,j,k,lzi) ! ZI
          qnz(i,k,15) = s4d(i,j,k,lzs) ! ZS
          qnz(i,k,16) = s4d(i,j,k,lzh) ! ZG
          qnz(i,k,17) = s4d(i,j,k,lzhl) ! ZH
          
          qnzm(i,k,13) = sm4d(i,j,k,lzr) ! ZR
          qnzm(i,k,14) = sm4d(i,j,k,lzi) ! ZI
          qnzm(i,k,15) = sm4d(i,j,k,lzs) ! ZS
          qnzm(i,k,16) = sm4d(i,j,k,lzh) ! ZG
          qnzm(i,k,17) = sm4d(i,j,k,lzhl) ! ZH
        END IF
          
      END DO
    END DO
    
    ! Microphysics calls
  
    if (ndebug .gt. 0) write(iunit,*) "DRIVER: ABOUT TO MAKE MICROPHYSICS CALL "
    
   
    CALL MYTMOM_MAIN(ww,tk(:,:,j),q(:,:,j),qnz,p(:,:,j),tkm(:,:,j),qm(:,:,j),qnzm,pm(:,:,j),lpr(:,j), &
                     spr(:,j),tem,zp2d,dt,nx,nz,j,mscheme,nscalar,ntc_my,cnor_my,cnos_my,  &
                     cnog_my,cnoh_my,rho_qi_my,rho_qs_my,rho_qg_my,rho_qh_my,alphar_my,         &
                     alphai_my,alphas_my,alphag_my,alphah_my,ixe)

  
    ! Copy back microphysics scalars
    
    DO k = 1,kze
      DO i = 1,ixe
      
        IF(mscheme >= 1) THEN 
          s4d(i,j,k,lc) = qnz(i,k,1)  ! QC
          s4d(i,j,k,lr) = qnz(i,k,2)  ! QR
          s4d(i,j,k,li) = qnz(i,k,3)  ! QI
          s4d(i,j,k,ls) = qnz(i,k,4)  ! QS
          s4d(i,j,k,lh) = qnz(i,k,5)  ! QG
          s4d(i,j,k,lhl) = qnz(i,k,6) ! QH
          
          sm4d(i,j,k,lc) = qnzm(i,k,1) ! QC
          sm4d(i,j,k,lr) = qnzm(i,k,2) ! QR
          sm4d(i,j,k,li) = qnzm(i,k,3) ! QI
          sm4d(i,j,k,ls) = qnzm(i,k,4) ! QS
          sm4d(i,j,k,lh) = qnzm(i,k,5) ! QG
          sm4d(i,j,k,lhl) = qnzm(i,k,6) ! QH
        END IF
        
        IF(mscheme >= 2) THEN
          s4d(i,j,k,lnc) = qnz(i,k,7) ! NC
          s4d(i,j,k,lnr) = qnz(i,k,8) ! NR
          s4d(i,j,k,lni) = qnz(i,k,9) ! NI
          s4d(i,j,k,lns) = qnz(i,k,10) ! NS
          s4d(i,j,k,lnh) = qnz(i,k,11) ! NG
          s4d(i,j,k,lnhl) = qnz(i,k,12) ! NH
          
          sm4d(i,j,k,lnc) = qnzm(i,k,7) ! NC
          sm4d(i,j,k,lnr) = qnzm(i,k,8) ! NR
          sm4d(i,j,k,lni) = qnzm(i,k,9) ! NI
          sm4d(i,j,k,lns) = qnzm(i,k,10) ! NS
          sm4d(i,j,k,lnh) = qnzm(i,k,11) ! NG
          sm4d(i,j,k,lnhl) = qnzm(i,k,12) ! NH
        END IF
        
        IF(mscheme == 4) THEN
          s4d(i,j,k,lzr) = qnz(i,k,13) ! ZR
          s4d(i,j,k,lzi) = qnz(i,k,14) ! ZI
          s4d(i,j,k,lzs) = qnz(i,k,15) ! ZS
          s4d(i,j,k,lzh) = qnz(i,k,16) ! ZG
          s4d(i,j,k,lzhl) = qnz(i,k,17) ! ZH
          
          sm4d(i,j,k,lzr) = qnzm(i,k,13) ! ZR
          sm4d(i,j,k,lzi) = qnzm(i,k,14) ! ZI
          sm4d(i,j,k,lzs) = qnzm(i,k,15) ! ZS
          sm4d(i,j,k,lzh) = qnzm(i,k,16) ! ZG
          sm4d(i,j,k,lzhl) = qnzm(i,k,17) ! ZH
        END IF
        
        ! Convert temperature back to potential temperature (first store temperature back in
        ! i,j,k order array for use in computing reflectivity below, although it is not used
        ! at this time, it may be in the future for determining wet vs. dry snow/graupel/hail).
        
        t1(i,j,k) = tk(i,k,j)
        s4d(i,j,k,lt) =  tk(i,k,j)/(pibase(k)+pi3d(i,j,k))
        sm4d(i,j,k,lt) = tkm(i,k,j)/(pibase(k)+pi3d(i,j,k))
        
        ! Copy back modified qv
        
        s4d(i,j,k,lv) = q(i,k,j)
        sm4d(i,j,k,lv) = qm(i,k,j)
        
        ! Compute density (needed for reflectivity calculation below) -- store in t0
        
        t0(i,j,k) = 1.0e5*(pibase(k)+pi3d(i,j,k))**2.509/(287.04*s4d(i,j,k,lt) )
        
        ! TODO: Need to figure out how to fill in precip rate and total precip arrays
      
        
      END DO
    END DO
  END DO ! Outer j-loop over slices

  ! Compute reflectivity
   
  IF (io_flag) THEN
    CALL reflec_MY(mscheme,nx,ny,nz,na,nor,norz,ntc_my,cnor_my,cnos_my,  &
                  cnog_my,cnoh_my,rho_qi_my,rho_qs_my,rho_qg_my,rho_qh_my,alphar_my,         &
                  alphai_my,alphas_my,alphag_my,alphah_my,                                   &
                  t0,                                                        &                                                                         
                  s4d,                                                     &
                  t1,                                                        &
                  dbz)
  END IF
  
  CALL cld_cpu('MICROPHYSICS')

  RETURN
END SUBROUTINE MY_DRIVE

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REFLEC_MY                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE reflec_MY(mscheme,nx,ny,nz,na,nor,norz,ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,      &
                      rhogrpl,rhohail,alpharain,alphaice,alphasnow,alphagrpl,           &
                      alphahail,rho,qscalar,t,rff)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates logarithmic radar reflectivity using
! equations based on the formulation described in Milbrandt and Yau (2005 Part I)
! It is designed for use with any of the 1, 2, or 3-moment versions of
! of the MY multi-moment microphysics scheme within ARPS.
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, Fall 2006.
!
! MODIFICATION HISTORY:
!   Dan Dawson 08/04/2011
!     Modified for use in the NSSL COMMAS model system.  Changed name to
!     reflec_MY.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  USE INDEX_MODULE
  USE COMMASMPI_MODULE, only: itile,jtile,ktile, &
                              ixbeg,ixend,nxend, &
                              jybeg,jyend,nyend, &
                              kzbeg,kzend,nzend  

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  integer, parameter    :: ng1 = 1

  INTEGER, INTENT(IN) :: mscheme
  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid
  INTEGER, INTENT(IN) :: na ! number of scalars
  INTEGER, INTENT(IN) :: nor,norz

  REAL, INTENT(IN) :: rho(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: qscalar(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz,na) ! Scalar microphysics array
  REAL, INTENT(IN) :: t(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz) ! Temperature (K)

  REAL, INTENT(OUT) :: rff(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz) ! Reflectivity (dBZ)

!-----------------------------------------------------------------------
! Declare local parameters.
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nqscalar=6       ! Number of microphysics species (6 for MY)
  REAL,PARAMETER :: pi = 3.14159265
  REAL,PARAMETER :: deratio = 0.224     ! Ratio of dielectric constants for ice and water
  REAL,PARAMETER :: cr = (pi/6.)*1000.  ! constant in mass power law relation for water particle

  REAL :: ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,alpharain,alphaice
  REAL :: alphasnow,alphagrpl,alphahail

!  REAL, PARAMETER :: rhor = 1000        ! Density of liquid water
!  REAL, PARAMETER :: rhoi = 500         ! Density of ice particles
!  REAL, PARAMETER :: rhos = 100         ! Density of snow
!  REAL, PARAMETER :: rhog = 400         ! Density of graupel
!  REAL, PARAMETER :: rhoh = 900         ! Density of hail

  REAL, PARAMETER :: rhor = 1000        ! Density of liquid water
  REAL :: rhoi          ! Density of ice particles
  REAL :: rhos         ! Density of snow
  REAL :: rhog         ! Density of graupel
  REAL :: rhoh         ! Density of hail

!  REAL, PARAMETER :: N0rfix = 1.0e6
!  REAL, PARAMETER :: N0sfix = 1.0e7
!  REAL, PARAMETER :: N0gfix = 4.0e5
!  REAL, PARAMETER :: N0hfix = 1.0e5

  REAL :: N0rfix
  REAL :: N0sfix
  REAL :: N0gfix
  REAL :: N0hfix

  ! Constants in diagnostic alpha relations

  REAL, PARAMETER :: c1r = 19.0
  REAL, PARAMETER :: c2r = 0.6
  REAL, PARAMETER :: c3r = 1.8
  REAL, PARAMETER :: c4r = 17.0
  REAL, PARAMETER :: c1i = 12.0
  REAL, PARAMETER :: c2i = 0.7
  REAL, PARAMETER :: c3i = 1.7
  REAL, PARAMETER :: c4i = 11.0
  REAL, PARAMETER :: c1s = 4.5
  REAL, PARAMETER :: c2s = 0.5
  REAL, PARAMETER :: c3s = 5.0
  REAL, PARAMETER :: c4s = 5.5
  REAL, PARAMETER :: c1g = 5.5
  REAL, PARAMETER :: c2g = 0.7
  REAL, PARAMETER :: c3g = 4.5
  REAL, PARAMETER :: c4g = 8.5
  REAL, PARAMETER :: c1h = 3.7
  REAL, PARAMETER :: c2h = 0.3
  REAL, PARAMETER :: c3h = 9.0
  REAL, PARAMETER :: c4h = 6.5
  REAL, PARAMETER :: c5h = 1.0
  REAL, PARAMETER :: c6h = 6.5


  REAL :: Ntr
  REAL :: Nti
  REAL :: Nts
  REAL :: Ntg
  REAL :: Nth
  REAL :: Zr                 ! rain reflectivity
  REAL :: Zi                 ! ice reflectivity
  REAL :: Zs                 ! snow reflectivity
  REAL :: Zsd                ! dry snow reflectivity
  REAL :: Zsw                ! wet snow reflectivity
  REAL :: Zg                 ! graupel reflectivity
  REAL :: Zgd                ! dry graupel reflectivity
  REAL :: Zgw                ! wet graupel reflectivity
  REAL :: Zh                 ! hail reflectivity
  REAL :: Zhd                ! dry hail reflectivity
  REAL :: Zhw                ! wet hail reflectivity
  REAL :: Zt                 ! total reflectivity
  REAL :: alpha              ! shape parameter in microphysical gamma distribution
  REAL :: dmx                ! mean-mass diameter of distribution
  REAL :: Gx                 ! parameter in radar equation

  REAL,PARAMETER :: epsQ = 1.0e-14
  REAL,PARAMETER :: epsN = 1.0e-3
  REAL,PARAMETER :: epsZ = 1.0e-32

  ! Note: for now, these must be set to the same as the values in the original
  ! run of the 2-moment version, or the reflectivity values will be incorrect.
  ! EDIT 08/04/08: values now passed in through common block in phycst.inc

!  REAL, PARAMETER :: alpharfix = 0.0
!  REAL, PARAMETER :: alphaifix = 0.0
!  REAL, PARAMETER :: alphasfix = 0.0
!  REAL, PARAMETER :: alphagfix = 0.0
!  REAL, PARAMETER :: alphahfix = 0.0

  REAL :: alpharfix
  REAL :: alphaifix
  REAL :: alphasfix
  REAL :: alphagfix
  REAL :: alphahfix

  INTEGER :: i,j,k,nq

  integer :: ixe, jye, kze
  integer, parameter :: istag = 1, jstag = 1, kstag = 1

! Beginning of executable code

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag
      
      kze = nz-1 ! assume no vertical tiling for now. Will fail otherwise.

  IF (ntcloud <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: ntcloud is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         ntcloud was reset to the default value of 1.0e8'

    ntcloud = 1.0e8
  END IF
  IF (n0rain <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0rain is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0rain was reset to the default value of 8.0e6'

    n0rain = 8.0e6
  END IF
  IF (n0snow <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0snow is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0snow was reset to the default value of 3.0e6'

    n0snow = 3.0e6
  END IF
  IF (n0grpl <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0grpl is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0grpl was reset to the default value of 4.0e5'

    n0snow = 4.0e5
  END IF
  IF (n0hail <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0hail is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0hail was reset to the default value of 4.0e4'

    n0hail = 4.0e4
  END IF
  IF (rhoice <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhoice is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhoice was reset to the default value of 500.'

    rhoice = 500.
  END IF
  IF (rhosnow <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhosnow is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhosnow was reset to the default value of 100.'

    rhosnow = 100.
  END IF
  IF (rhogrpl <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhogrpl is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhogrpl was reset to the default value of 400.'

    rhogrpl = 400.
  END IF
  IF (rhohail <= 0.0) THEN
    WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhohail is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhohail was reset to the default value of 913.'

    rhohail = 913.
  END IF

  rhoi = rhoice
  rhos = rhosnow
  rhog = rhogrpl
  rhoh = rhohail
  N0rfix = n0rain
  N0sfix = n0snow
  N0gfix = n0grpl
  N0hfix = n0hail
  alpharfix = alpharain
  alphaifix = alphaice
  alphasfix = alphasnow
  alphagfix = alphagrpl
  alphahfix = alphahail

  Zr = 0.0
  Zi = 0.0
  Zs = 0.0
  Zg = 0.0
  Zh = 0.0

  rff = 0.0

  DO k = 1,kze
    DO j = 1,jye
      DO i = 1,ixe

      Zr = 0.0
      Zi = 0.0
      Zs = 0.0
      Zg = 0.0
      Zh = 0.0

        IF(mscheme == 4) THEN   ! 3-moment case: reflectivity predicted
          IF(qscalar(i,j,k,lzr) >= epsZ) THEN
            Zr = (((pi/6.)*rhor/cr)**2.)*qscalar(i,j,k,lzr)
            !print*,'rain block,qr,nr,Zr_raw,Zr',qscalar(i,j,k,nq),qscalar(i,j,k,nq+nqscalar),qscalar(i,j,k,nq+2*nqscalar-1),Zr
          END IF
          IF(qscalar(i,j,k,lzi) >= epsZ) THEN
            Zi = deratio*((440./cr)**2.)*qscalar(i,j,k,lzi)
          END IF
          IF(qscalar(i,j,k,lzs) >= epsZ) THEN
  !              IF(t(i,j,k) <= 273.15) THEN
              Zs = deratio*(((pi/6.)*rhos/cr)**2.)*qscalar(i,j,k,lzs)
  !              ELSE
  !                Zs = (((pi/6.)*rhos/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
  !              END IF
          END IF
          IF(qscalar(i,j,k,lzh) >= epsZ) THEN
            Zgd = deratio*(((pi/6.)*rhog/cr)**2.)*qscalar(i,j,k,lzh)
  !              Zgw = (((pi/6.)*rhog/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
  !              IF(t(i,j,k) <= 275.65 .and. t(i,j,k) >= 270.65) THEN
  !                Zg = Zgw*(t(i,j,k)-270.65)/5.0 + Zgd*(275.65-t(i,j,k))/5.0
  !              ELSE IF (t(i,j,k) < 270.65) THEN
            Zg = Zgd
  !              ELSE
  !                Zg = Zgw
  !              END IF
          END IF
          IF(qscalar(i,j,k,lzhl) >= epsZ) THEN
            Zhd = deratio*(((pi/6.)*rhoh/cr)**2.)*qscalar(i,j,k,lzhl)
  !              Zhw = (((pi/6.)*rhoh/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
  !              IF(t(i,j,k) <= 275.65 .and. t(i,j,k) >= 270.65) THEN
  !                Zh = Zhw*(t(i,j,k)-270.65)/5.0 + Zhd*(275.65-t(i,j,k))/5.0
  !              ELSE IF (t(i,j,k) < 270.65) THEN
            Zh = Zhd
  !              ELSE
  !                Zh = Zhw
  !              END IF
          END IF
        ELSE IF(mscheme == 3) THEN     ! 2-moment with diagnostic alpha
          IF(qscalar(i,j,k,lr) >= epsQ .and. qscalar(i,j,k,lnr) >= epsN) THEN
            dmx = (rho(i,j,k)*qscalar(i,j,k,lr)/((pi/6.)*rhor*qscalar(i,j,k,lnr)))**(1./3.)
            alpha = c1r*TANH(c2r*(dmx-c3r))+c4r
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zr = ((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lr))**2.)/(qscalar(i,j,k,lnr))
          END IF
          IF(qscalar(i,j,k,li) >= epsQ .and. qscalar(i,j,k,lni) >= epsN) THEN
            dmx = (rho(i,j,k)*qscalar(i,j,k,li)/((pi/6.)*rhoi*qscalar(i,j,k,lni)))**(1./3.)
            alpha = c1i*TANH(c2i*(dmx-c3i))+c4i
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zi = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,li))**2.)/(qscalar(i,j,k,lni))
          END IF
          IF(qscalar(i,j,k,ls) >= epsQ .and. qscalar(i,j,k,lns) >= epsN) THEN
            dmx = (rho(i,j,k)*qscalar(i,j,k,ls)/((pi/6.)*rhos*qscalar(i,j,k,lns)))**(1./3.)
            alpha = c1s*TANH(c2s*(dmx-c3s))+c4s
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zs = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,ls))**2.)/(qscalar(i,j,k,lns))
          END IF
          IF(qscalar(i,j,k,lh) >= epsQ .and. qscalar(i,j,k,lnh) >= epsN) THEN
            dmx = (rho(i,j,k)*qscalar(i,j,k,lh)/((pi/6.)*rhog*qscalar(i,j,k,lnh)))**(1./3.)
            alpha = c1g*TANH(c2g*(dmx-c3g))+c4g
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zg = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lh))**2.)/(qscalar(i,j,k,lnh))
          END IF
          IF(qscalar(i,j,k,lhl) >= epsQ .and. qscalar(i,j,k,lnhl) >= epsN) THEN
            dmx = (rho(i,j,k)*qscalar(i,j,k,lhl)/((pi/6.)*rhoh*qscalar(i,j,k,lnhl)))**(1./3.)
            IF(dmx < 8e-3) THEN
              alpha = c1h*TANH(c2h*(dmx-c3h))+c4h
            ELSE
              alpha = c5h*dmx-c6h
            END IF
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zh = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lhl))**2.)/(qscalar(i,j,k,lnhl))
          END IF
        ELSE IF(mscheme == 2) THEN     ! 2-moment scheme with fixed alpha
          IF(qscalar(i,j,k,lr) >= epsQ .and. qscalar(i,j,k,lnr) >= epsN) THEN
            alpha = alpharfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zr = ((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lr))**2.)/(qscalar(i,j,k,lnr))
            !print*,'rain block,i,j,k,alpha,Gx,rho,qr,nr',i,j,k,alpha,Gx,rho(i,j,k),qscalar(i,j,k,nq),qscalar(i,j,k,nq+nqscalar),10.0*LOG10(1.0e18*Zr)
          END IF
          IF(qscalar(i,j,k,li) >= epsQ .and. qscalar(i,j,k,lni) >= epsN) THEN
            alpha = alphaifix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zi = deratio*((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,li))**2.)/(qscalar(i,j,k,lni))
            !  print*,'ice block'
          END IF
          IF(qscalar(i,j,k,ls) >= epsQ .and. qscalar(i,j,k,lns) >= epsN) THEN
            alpha = alphasfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zs = deratio*((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,ls))**2.)/(qscalar(i,j,k,lns))
            !  print*,'snow block'
          END IF
          IF(qscalar(i,j,k,lh) >= epsQ .and. qscalar(i,j,k,lnh) >= epsN) THEN
            alpha = alphagfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zg = deratio*((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lh))**2.)/(qscalar(i,j,k,lnh))
            !  print*,'graupel block'
          END IF
          IF(qscalar(i,j,k,lhl) >= epsQ .and. qscalar(i,j,k,lnhl) >= epsN) THEN
            alpha = alphahfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            Zh = deratio*((1./cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lhl))**2.)/(qscalar(i,j,k,lnhl))
            !  print*,'hail block'
          END IF
        ELSE IF(mscheme == 1) THEN     ! 1-moment scheme with fixed N0x (and alpha)
          IF(qscalar(i,j,k,lr) >= epsQ) THEN
            Ntr = N0rfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,lr)/(pi*rhor))**(1./4.)
            alpha = alpharfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            IF(Ntr >= epsN) THEN
              Zr = ((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lr))**2.)/Ntr
            END IF
          END IF
          IF(qscalar(i,j,k,li) >= epsQ) THEN
            Nti = 5.*exp(0.304*(273.15-max(233.,t(i,j,k))))
            alpha = alphaifix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            IF(Nti >= epsN) THEN
              Zi = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,li))**2.)/Nti
            END IF
          END IF
          IF(qscalar(i,j,k,ls) >= epsQ) THEN
            Nts = N0sfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,ls)/(pi*rhos))**(1./4.)
            alpha = alphasfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            IF(Nts >= epsN) THEN
              Zs = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,ls))**2.)/Nts
            END IF
          END IF
          IF(qscalar(i,j,k,lh) >= epsQ) THEN
            Ntg = N0gfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,lh)/(pi*rhog))**(1./4.)
            alpha = alphagfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            IF(Ntg >= epsN) THEN
              Zg = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lh))**2.)/Ntg
            END IF
          END IF
          IF(qscalar(i,j,k,lhl) >= epsQ) THEN
            Nth = N0hfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,lhl)/(pi*rhoh))**(1./4.)
            alpha = alphahfix
            Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
            IF(Nth >= epsN) THEN
            Zh = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,lhl))**2.)/Nth
            END IF
          END IF
        END IF ! which microphysics scheme option

        ! Sum up contributions from all hydrometeor species and convert
        ! to logarithmic reflectivity

        Zt = 1.0e18*(Zr+Zi+Zs+Zg+Zh)   ! Now in units of mm^6*m^-3

        !IF(Zt >= epsZ) THEN
        !  print*,'Zt = ',Zt
        !END IF
        rff(i,j,k) = 10.0*LOG10(MAX(Zt,1.0))   ! Now in units of dBZ


        !IF(rff(i,j,k) >= 0.0) THEN
        !  print*,rff(i,j,k)
        !END IF

      END DO  ! k loop
    END DO  ! j loop
  END DO  ! i loop
  RETURN
END SUBROUTINE reflec_MY
