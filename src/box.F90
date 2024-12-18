!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN          \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE FOLLOW    ////////////////////
!
!
! This subroutine either adds or subtracts storm motion from the horizontal 
! velocity field.
!-----------------------------------------------------------------------------
  MODULE FOLLOW_MODULE
  
  CONTAINS
  
  SUBROUTINE FOLLOW(gd,sign,nxend1)

   USE GRID_MODULE

   USE COMMASMPI_MODULE

   implicit none

   TYPE(GRID) :: gd
   real       :: sign
   integer, optional :: nxend1
   
   real, pointer :: u(:,:,:), v(:,:,:)
   integer       :: nx, ny, nz, ng
   real          :: ugrid, vgrid
   
   integer :: i,j,k

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   IF ( .not. present(nxend1) ) THEN
     CALL GET_VARIABLE(gd, 'NX', nx)
     CALL GET_VARIABLE(gd, 'NY', ny)
     CALL GET_VARIABLE(gd, 'NZ', nz)
   ELSE
     CALL GET_VARIABLE(gd, 'NXEND', nx)
     CALL GET_VARIABLE(gd, 'NYEND', ny)
     CALL GET_VARIABLE(gd, 'NZEND', nz)
     
     IF ( ncxe .ge. 1 ) THEN
     nx = Min( nx, ncxe )
     ny = Min( ny, ncye )
     nz = Min( nz, ncze )
     ENDIF
     
   ENDIF
   CALL GET_VARIABLE(gd, 'NG', ng)
   CALL GET_VARIABLE(gd, 'UGRID', ugrid)
   CALL GET_VARIABLE(gd, 'VGRID', vgrid)
   
   CALL GET_VARIABLE(gd, 'U', u)
   CALL GET_VARIABLE(gd, 'V', v)

   IF( sign .gt. 0.0 ) THEN

#ifdef MPI
       ixb = -ng+1
       ixe = itile+ng

       jyb = -ng+1
       jye = jtile+ng

       kzb = -ng+1
       kze = ktile+ng

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =1,nz-1 ; DO j = 1,ny ; DO i = 1,nx
#endif
         u(i,j,k) = u(i,j,k) + ugrid
         v(i,j,k) = v(i,j,k) + vgrid
       ENDDO ; ENDDO ; ENDDO

   ELSE

#ifdef MPI
       ixb = -ng+1
       ixe = itile+ng

       jyb = -ng+1
       jye = jtile+ng

       kzb = -ng+1
       kze = ktile+ng

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =1,nz-1 ; DO j = 1,ny ; DO i = 1,nx 
#endif
         u(i,j,k) = u(i,j,k) - ugrid
         v(i,j,k) = v(i,j,k) - vgrid

       ENDDO ; ENDDO ; ENDDO

   ENDIF


  RETURN
  END SUBROUTINE FOLLOW
  
  END MODULE FOLLOW_MODULE

!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN        \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE BOX     ////////////////////
!
!
! Implements the storm tracking algorithm suggested by Proctor for TASS
!  Proctor, F. H., The Terminal Area Simulation System / Volume 1: 
!       Theoretical Formulation, NASA Technical Report 4046. 1987. (See page 86)
!-----------------------------------------------------------------------------
!
! Created by LJW: 07-19-02
! APS JUN07: Not converted to MPI since currently not being called
! ERM March 2019: Updated for COMMAS data structure and converted to MPI
!-----------------------------------------------------------------------------
      SUBROUTINE BOX(gd, u,v,w,gx,gy,gz,                                   &
                     xold,yold,uold,vold,ugrid,vgrid, piinit,sinit,      &
                     xmid,ymid, track_tau, track_alpha, track_zeta0,     &
                     track_xmid_delta, track_ymid_delta,luno)

      USE GRID_MODULE
      USE PARAM_MODULE, only : rd
      USE COMMASMPI_MODULE
#ifdef MPI
    USE mpi
#endif


      implicit none 

! Passed variable declarations

   TYPE(GRID) :: gd
!   real       :: sign
!   integer, optional :: nxend1

   TYPE(VARIABLE)     :: u
   TYPE(VARIABLE)     :: v
   TYPE(VARIABLE)     :: w
   TYPE(VARIABLE)     :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE)     :: piinit
   TYPE(VARIABLE)     :: sinit(2)     ! We are going to try and create arrays of the scalar variables needed here
   
   integer :: luno ! stdout file unit
   
   integer       :: nx, ny, nz, ng

!      integer is, js, ks,nx1d,ny1d,nz1d
!      real u(nx,ny,nz)
!      real v(nx,ny,nz)
!      real w(nx,ny,nz)

!      real x1d(nx,nx1d) 
!      real y1d(ny,ny1d) 
!      real z1d(nz,nz1d) 

      real xold,yold,uold,vold,ugrid,vgrid,time
      real :: xmid, ymid, track_tau, track_alpha, track_zeta0, track_xmid_delta, track_ymid_delta

! Local variable declarations
      
      integer i, j, k, n, im1, jm1, ip1, jp1

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      real, allocatable :: den(:)

      real :: zeta, zeta0, alpha, tau, unew, vnew
      real :: dv
      double precision :: xbar, qxbar, ybar, qybar, qbar, q, eps

      parameter( eps   = 1.0e-15 )
!      parameter( zeta0 = 0.001 )
!      parameter( xmid  = 50000., ymid = 55000. )
!      parameter( alpha = 0.75 )
!      parameter( tau   = 200. )


      logical, parameter :: debugme = .false.
      
      integer, parameter :: ntot = 5
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      real :: wzz
      
      
! Vorticity function statements

   wzz(i,j,k,im1,jm1)=(v%flt3d(i,j,k)-v%flt3d(im1,j,k))*gx(4)%flt1d(i)-(u%flt3d(i,j,k)-u%flt3d(i,jm1,k))*gy(4)%flt1d(j)
      

     CALL GET_VARIABLE(gd, 'NX', nx)
     CALL GET_VARIABLE(gd, 'NY', ny)
     CALL GET_VARIABLE(gd, 'NZ', nz)
     
     IF (debugme) write(0,*) 'allocate den, nz = ',nz
     allocate( den(nz) )
     IF (debugme) write(0,*) 'allocated den'

!   CALL GET_VARIABLE(gd, 'NG', ng)
!   CALL GET_VARIABLE(gd, 'UGRID', ugrid)
!   CALL GET_VARIABLE(gd, 'VGRID', vgrid)
   
!   CALL GET_VARIABLE(gd, 'U', u)
!   CALL GET_VARIABLE(gd, 'V', v)
!   CALL GET_VARIABLE(gd, 'W', w)


        tau   = track_tau
        alpha = track_alpha
        zeta0 = track_zeta0
! Init sum variables

      qbar  = 0.0
      qxbar = 0.0
      qybar = 0.0

! Compute where storm is above the boundary layer through mid levels

   IF (debugme) write(0,*) 'set den'
   DO k = 1,nz-1
    IF (debugme) write(0,*) 'k = ',k
    IF (debugme) write(0,*) 'pinit = ',piinit%flt1d(k)
    IF (debugme) write(0,*) 'sinit(1) = ',sinit(1)%flt1d(k)
    IF (debugme) write(0,*) 'sinit(2) = ',sinit(2)%flt1d(k)
    
    den(k) = 1.0e5*piinit%flt1d(k)**2.509/(rd*sinit(1)%flt1d(k) * (1.0+0.61*sinit(2)%flt1d(k)))
   ENDDO
! Compute vorticity

!      DO k = 1,nz/2
!       DO j = 3,ny-3
!        DO i = 3,nx-3

      kzb = 1
      kze = ktile/2
      if (kzend .eq. nzend) kze = (kzend-kzbeg)/2
      
      IF ( gz(3)%flt1d(kze) < 10000. .and. gz(3)%flt1d(ktile-1) > 13000. ) THEN
        DO k = 1,ktile-1
          IF ( gz(3)%flt1d(k) <= 10100. ) THEN 
            kze = k
          ELSE
            exit
          ENDIF
        ENDDO
      ENDIF

      jyb = 1
!      IF ( jybeg == nybeg ) jyb = 3
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
!      IF ( ixbeg == nxbeg ) ixb = 3
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      IF (debugme) write(0,*) 'zeta loop'
      
      do k = kzb,kze

       do j = jyb,jye
        jm1 = j-1
        jp1 = j+1
        if (jybeg .eq. nybeg) jm1 = max(j-1,1)
        if (jyend .eq. nyend) jp1 = min(j+1,ny-1)

        do i = ixb,ixe
          im1 = i-1
          ip1 = i+1
          if (ixbeg .eq. nxbeg) im1 = max(i-1,1)
          if (ixend .eq. nxend) ip1 = min(i+1,nx-1)
          
          IF ( track_xmid_delta > 0. ) THEN
             IF ( .not. ( (gx(1)%flt1d(i) > xmid - track_xmid_delta) .and. &
                ( gx(1)%flt1d(i) < xmid + track_xmid_delta) ) ) THEN
               CYCLE
             ENDIF
          ENDIF

          IF ( track_ymid_delta > 0. ) THEN
             IF ( .not. ( (gy(1)%flt1d(j) > ymid - track_ymid_delta) .and. &
                ( gy(1)%flt1d(j) < ymid + track_ymid_delta) ) ) THEN
               CYCLE
             ENDIF
          ENDIF
          
          dv = 1.0/(gx(3)%flt1d(i)*gy(3)%flt1d(j)*gz(3)%flt1d(k) )

!         zeta = (v(j,  k,i  )-v(j,  k,i-1))*x1d(i,  4) - (u(j,  k,i  )-u(j-1,k,i  ))*y1d(j,  4)
!     $        + (v(j+1,k,i  )-v(j+1,k,i-1))*x1d(i,  4) - (u(j+1,k,i  )-u(j,  k,i  ))*y1d(j+1,4)
!     $        + (v(j,  k,i+1)-v(j,  k,i  ))*x1d(i+1,4) - (u(j,  k,i+1)-u(j-1,k,i+1))*y1d(j,  4)
!     $        + (v(j+1,k,i+1)-v(j+1,k,i  ))*x1d(i+1,4) - (u(j+1,k,i+1)-u(j,  k,i+1))*y1d(j+1,4)

!         zeta = (v(i,  j,  k)-v(im1,j,  k))*x1d(i,  4) - (u(i,  j,  k)-u(i,  j-1,k))*y1d(j,  4)  &
!              + (v(i,  j+1,k)-v(im1,j+1,k))*x1d(i,  4) - (u(i,  j+1,k)-u(i,  j,  k))*y1d(j+1,4)  &
!              + (v(i+1,j,  k)-v(i,  j,  k))*x1d(i+1,4) - (u(i+1,j,  k)-u(i+1,j-1,k))*y1d(j,  4)  &
!              + (v(i+1,j+1,k)-v(i,  j+1,k))*x1d(i+1,4) - (u(i+1,j+1,k)-u(i+1,j,  k))*y1d(j+1,4)

!         zeta = zeta/4.0

          zeta  = 0.25 * (wzz(i,j,k,im1,jm1) + wzz(i,jp1,k,im1,j)+wzz(ip1,j,k,i,jm1) + wzz(ip1,jp1,k,i,j))

         IF( zeta .gt. zeta0 ) THEN

          q = 1.0 + alog10( zeta0 ) + alog10( zeta ) * den(k) * w%flt3d(i,j,k)

         ELSE

          q = den(k) * w%flt3d(i,j,k)

         ENDIF

         qxbar = qxbar + q*q*gx(1)%flt1d(i)*dv ! x1d(i,1)
         qybar = qybar + q*q*gy(1)%flt1d(j)*dv  ! y1d(j,1)
         qbar  = qbar  + q*q*dv

        ENDDO
       ENDDO

      ENDDO
      
      
#ifdef MPI

    IF (debugme) write(0,*) 'allreduce'
    
    mpitotindp(1)  = qxbar
    mpitotindp(2)  = qybar
    mpitotindp(3)  = qbar

!   write(luno,*) 'tile tmp1,tmp2 = ',tmp1,tmp2,my_rank
  CALL MPI_Allreduce(mpitotindp, mpitotoutdp, 3, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)
  
  qxbar = mpitotoutdp(1)
  qybar = mpitotoutdp(2)
  qbar  = mpitotoutdp(3)

#endif
      xbar = qxbar / (qbar+eps)
      ybar = qybar / (qbar+eps)

      IF ( xold == 0.0 ) xold = xbar
      IF ( yold == 0.0 ) yold = ybar
      unew = uold + (alpha*(xbar-xold) + (1.0-alpha)*(xbar-xmid)) / tau
      vnew = vold + (alpha*(ybar-yold) + (1.0-alpha)*(ybar-ymid)) / tau

!     unew = uold + ((xbar-xold) + 0.25*(xbar-xmid)) / tau
!     vnew = vold + ((ybar-yold) + 0.25*(ybar-ymid)) / tau

    IF ( my_rank == 0 ) THEN
      write(luno,*) 
      write(luno,*)  'SUBROUTINE BOX'
      write(luno,*) 
      write(luno,*)  'XOLD  = ',xold
      write(luno,*)  'YOLD  = ',yold
      write(luno,*)  'UOLD  = ',uold
      write(luno,*)  'VOLD  = ',vold
      write(luno,*) 
      write(luno,*)  'XBAR  = ',xbar
      write(luno,*)  'YBAR  = ',ybar
      write(luno,*)  'UNEW  = ',unew
      write(luno,*)  'VNEW  = ',vnew
    ENDIF

      ugrid = unew
      vgrid = vnew

      xold = xbar
      yold = ybar

      IF (debugme) write(0,*) 'deallocate'
      deallocate( den )
      IF (debugme) write(0,*) 'return'

      RETURN
      END SUBROUTINE BOX


!-----------------------------------------------------------------------------------------------------------------------------------
!  
!  Base state substitution (BSS) NAMELIST
!

MODULE BSS_NML

!
! Base State Substition (BSS) (Letkewicz et al. 2013, MWR)
!
!  This sets up BSS with up to maxnumbss different soundings. Here, BSS is done as 
!  the model runs, without restarts. BSS is done fractionally on each time
!  step over the time period(s) specified by 'startbss' and 'endbss'.  Note that 
!  overlapping time periods are not checked. For example, if the first BSS goes from 
!  1000 s to 1500s, and the second goes from 1400 s to 1600s, the first one will end
!  at 1400 s (and be incomplete).
!
!  Default values are shown. For multiple BSS, simply add to the appropriate namelist variables for
!  non-default values. Mean U and V can be adjusted with 'umeanbss' and 'vmeanbss' if they 
!  are set to values > -100 (m/s).

  implicit none
  
  integer, parameter :: maxnumbss = 10
  
  integer :: numbss = 0                   ! number of BSS
  integer :: startbss(maxnumbss) = 0      ! starting time for each sounding
  integer :: endbss(maxnumbss)   = -1     ! ending time for each sounding 
  integer :: sndtypebss(maxnumbss) = 1    ! as for 'sndtype' in homog_init: 1 = WK, 2 = read in sounding
  integer :: wtypebss(maxnumbss) = 1      ! WK analytical sounding type (as for 'wtype' in homog_init)
  integer :: idoqbss(maxnumbss) = 0       ! substitute QV (0 = off, 1 = on)
  integer :: idopbss(maxnumbss) = 0       ! substitute pressure (pi) (0 = off, 1 = on)
  integer :: idotbss(maxnumbss) = 0       ! substitute theta (0 = off, 1 = on)
  integer :: idoubss(maxnumbss) = 1       ! substitute U (0 = off, 1 = on)
  integer :: idovbss(maxnumbss) = 1       ! substitute V (0 = off, 1 = on)
  real    :: usbss(maxnumbss) = 0.0       ! as for 'Us' in homog_init (sndtype=1)
  real    :: uslbss(maxnumbss) = 0.0      ! as for 'Usl' in homog_init (sndtype=1)
  real    :: usmvbss(maxnumbss) = 0.0      ! as for 'Usl' in homog_init (sndtype=1)
  real    :: uzbss(maxnumbss) = 2500.     ! as for 'Uz' in homog_init (sndtype=1)
  real    :: ubasebss(maxnumbss) = 0.0    ! as for 'Ubase' in homog_init (sndtype=1)
  real    :: psfcbss(maxnumbss+1) = 100000. ! as for 'psfc' in homog_init (sndtype=1)
  real    :: tsfcbss(maxnumbss+1)   = 300.  ! as for 'tsfc' in homog_init (sndtype=1)
  real    :: tshiftbss(maxnumbss+1)   = 0.  ! temperature shift in entire column
  real    :: qsfcbss(maxnumbss+1)   = 14.   ! as for 'qsfc' in homog_init (sndtype=1)
  real    :: rhmaxbss(maxnumbss)  = 0.90  ! as for 'rhmax' in homog_init (sndtype=1)
  real    :: rhmax2bss(maxnumbss) = 0.25  ! as for 'rhmax2' in homog_init (sndtype=1)
  real    :: zrhmax2bss(maxnumbss) = -100. ! as for 'zrhmax2' in homog_init (sndtype=1)
  integer :: shapebss(maxnumbss) = 0       ! as for 'shape' in homog_init (sndtype=1)
  real    :: dudz0bss(maxnumbss) = 0.0     ! as for 'dudz0' in homog_init (sndtype=1)
  real    :: dudz1bss(maxnumbss) = 0.0     ! as for 'dudz1' in homog_init (sndtype=1)
  real    :: dudz2bss(maxnumbss) = 0.0     ! as for 'dudz2' in homog_init (sndtype=1)
  real    :: z0bss(maxnumbss)    = 99000.  ! as for 'z0' in homog_init (sndtype=1)
  real    :: z1bss(maxnumbss)    = 98000.  ! as for 'z1' in homog_init (sndtype=1)
  real    :: z2bss(maxnumbss)    = 97000.  ! as for 'z2' in homog_init (sndtype=1)
  real    :: umeanbss(maxnumbss) = -1000.0 ! set mean U in depth of domain
  real    :: vmeanbss(maxnumbss) = -1000.0 ! ! adjust mean U in depth of domain mean V in depth of domain
  real    :: rhnewbss(maxnumbss) = -1.     ! specified RH (e.g., 0.5 for 50%) for a layer defined by rhz1bss and rhz2bss
  real    :: rhz1bss(maxnumbss)  = -1.     ! lower bound of RH layer
  real    :: rhz2bss(maxnumbss)  = -1.     !  upper bound of RH layer
  real    :: ireplaceuvbss(maxnumbss) = 0     ! whether to replace U and V in an input sounding with U and V from the WK settings
  real    :: ireplaceqbss(maxnumbss)  = 0     ! whether to replace qv in an input sounding with qv from the WK settings
  character(LEN = 100) :: sndfilebss(maxnumbss)   = ' ' ! as for 'sndfile' in homog_init (sndtype=2)


  
  
  logical :: bssinit = .false.
  integer :: i_current_bss = 1
  
  character(len=255) :: sndnamebss(maxnumbss)
  
  real, allocatable :: u1dbss(:,:),v1dbss(:,:),t1dbss(:,:),p1dbss(:,:),q1dbss(:,:)
  

  NAMELIST /bss_init/              &
                 numbss,      &
                 startbss,    &
                 endbss,      &
                 sndtypebss,  &
                 wtypebss,  &
                 usbss,     &
                 uslbss,    &
                 uzbss,  &
                 ubasebss,  &
                 psfcbss,  &
                 tsfcbss,  &
                 tshiftbss,  &
                 qsfcbss,  &
                 rhmaxbss,  &
                 rhmax2bss,  &
                 zrhmax2bss,  &
                 shapebss,  &
                 dudz0bss,  &
                 dudz1bss,  &
                 dudz2bss,  &
                 z0bss,  &
                 z1bss,  &
                 z2bss,  &
                 umeanbss, &
                 vmeanbss, &
                 sndfilebss, &
                 idoqbss, &
                 idopbss, &
                 idotbss, &
                 idoubss, &
                 idovbss, &
                 rhnewbss, &
                 rhz1bss,  &
                 rhz2bss,  &
                 ireplaceuvbss, &
                 ireplaceqbss


                 
 CONTAINS
 
   SUBROUTINE BSS_MODULE_INIT()
   END SUBROUTINE BSS_MODULE_INIT


END MODULE BSS_NML




!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN        \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE BSS     ////////////////////
!
!-----------------------------------------------------------------------------

   SUBROUTINE BSS(gd, time_real)

   USE GRID_MODULE
   USE BSS_NML
   USE COMMASMPI_MODULE, only: my_rank,nzbeg,nzend,kzbeg,kzend,commasmpi_abort, &
                               itile,jtile,ktile
   USE INIT_MODULE, only: z1dinit,zrhdel,wkumax1,wkumax2,wkvmax,wkh1,wkh2
   USE PARAM_MODULE, only: cp,rd

   implicit none

   TYPE(GRID) :: gd
   real       :: time_real
   
   real, pointer :: u(:,:,:), v(:,:,:), qv3(:,:,:), th3(:,:,:), qc3(:,:,:),pi3(:,:,:)
   integer       :: nx, ny, nz, ng
   real          :: dt
   real          :: ugrid, vgrid
   
   TYPE(ATTRIBUTE), pointer   :: prefix
   real, pointer :: zc(:)
   real, pointer :: ze(:)
   real, pointer :: dzc(:)
   real, pointer :: dze(:)
   real, pointer :: u1d1(:)
   real, pointer :: v1d1(:)
   real, pointer :: t1d1(:)
   real, pointer :: p1d1(:)
   real, pointer :: q1d1(:)

   real, pointer :: u1d1b(:)
   real, pointer :: v1d1b(:)
   real, pointer :: tz1b(:)
   real, pointer :: qz1b(:)
   real, pointer :: pz1b(:)

   integer :: i,j,k,n
   logical :: work_to_do
   character(LEN=255)         :: outsoundfile, outsoundfiletmp
   character(len=2)           :: numstr
   real    :: rotdeg = 0
   integer :: ntr = 1
   real    :: sum, height,dudz
   real    :: deltatstep,deltat, deltas
   real    :: psfc,pres,qvs,rh
   integer :: indx
   
!   real, allocatable :: upert(:,:,:),vpert(:,:,:),tpert(:,:,:),qpert(:,:,:),pipert(:,:,:)
   
   real, allocatable, save :: udel(:),vdel(:),tdel(:),qdel(:),pidel(:)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   
   IF ( numbss == 0 ) Return

   CALL GET_VARIABLE(gd, 'NX', nx)
   CALL GET_VARIABLE(gd, 'NY', ny)
   CALL GET_VARIABLE(gd, 'NZ', nz)
   CALL GET_VARIABLE(gd, 'DT', dt)
     
   CALL GET_VARIABLE(gd, 'NG', ng)
   CALL GET_VARIABLE(gd, 'UGRID', ugrid)
   CALL GET_VARIABLE(gd, 'VGRID', vgrid)
   
   IF ( .not. bssinit ) THEN !{
! initialize BSS arrays 
    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', prefix)

   CALL GET_VARIABLE(gd,'ZC',   zc)
   CALL GET_VARIABLE(gd,'ZE',   ze)
   CALL GET_VARIABLE(gd,'DZC',  dzc)
   CALL GET_VARIABLE(gd,'DZE',  dze)


      bssinit = .true.
      
      allocate( u1dbss(-ng+1:nzend+ng,0:numbss+1),v1dbss(-ng+1:nzend+ng,0:numbss+1), &
                t1dbss(-ng+1:nzend+ng,0:numbss+1),p1dbss(-ng+1:nzend+ng,0:numbss+1), &
                q1dbss(-ng+1:nzend+ng,0:numbss+1) )

      allocate( udel(-ng+1:nzend+ng), vdel(-ng+1:nzend+ng), tdel(-ng+1:nzend+ng), &
                qdel(-ng+1:nzend+ng), pidel(-ng+1:nzend+ng) )

! get original base state and put into slot 0

      CALL GET_VARIABLE(gd,'UINIT0', u1d1b)
      CALL GET_VARIABLE(gd,'VINIT0', v1d1b)
      CALL GET_VARIABLE(gd,'THINIT0',tz1b)
      CALL GET_VARIABLE(gd,'QVINIT0',qz1b)
      CALL GET_VARIABLE(gd,'PIINIT0',pz1b)

      DO k = 1,nzend
       u1dbss(k,0) = u1d1b(k)
       v1dbss(k,0) = v1d1b(k)
       t1dbss(k,0) = tz1b(k) 
       q1dbss(k,0) = qz1b(k) 
       p1dbss(k,0) = pz1b(k) 
      ENDDO

      DO n = 1,numbss
       
       write(numstr,'(i2.2)') n
       outsoundfile = trim(prefix%str)//'.bss.'//numstr//'.sound'
        IF ( sndtypebss(n) == 1 ) THEN
          call snd1x(gd,t1dbss(-ng+1,n),p1dbss(-ng+1,n),q1dbss(-ng+1,n),dzc,dze,zc, &
                  nzend,ng,psfcbss(n),tsfcbss(n),qsfcbss(n),rhmaxbss(n),rhmax2bss(n),  &
                   zrhmax2bss(n),zrhdel,wtypebss(n),Usbss(n),Uzbss(n),Uslbss(n),Usmvbss(n),Ubasebss(n),  &
                   ntr,u1dbss(-ng+1,n),v1dbss(-ng+1,n),z0bss(n),z1bss(n),dudz0bss(n),rotdeg, &
                   outsoundfile,wkumax1,wkumax2,wkvmax,wkh1,wkh2) 
        ELSEIF ( sndtypebss(n) == 2 ) THEN
           IF ( ireplaceuvbss(n) >= 1 .or. ireplaceqbss(n) >= 1 ) THEN
            outsoundfiletmp = '/dev/null'
!            outsoundfiletmp = 'dum.sound'
            call snd1x(gd,t1dbss(-ng+1,numbss+1),p1dbss(-ng+1,numbss+1),q1dbss(-ng+1,numbss+1),dzc,dze,zc, &
                  nzend,ng,psfcbss(n),tsfcbss(n),qsfcbss(n),rhmaxbss(n),rhmax2bss(n),  &
                   zrhmax2bss(n),zrhdel,wtypebss(n),Usbss(n),Uzbss(n),Uslbss(n),Usmvbss(n),Ubasebss(n),  &
                   ntr,u1dbss(-ng+1,numbss+1),v1dbss(-ng+1,numbss+1),z0bss(n),z1bss(n),dudz0bss(n),rotdeg, &
                   outsoundfiletmp,wkumax1,wkumax2,wkvmax,wkh1,wkh2) 
!            call snd1x(gd,t1dbss(-ng+1,numbss+1),p1dbss(-ng+1,numbss+1),q1dbss(-ng+1,numbss+1),dzc,dze,zc, &
!                  nzend,ng,psfcbss(n),tsfcbss(n),qsfcbss(n),rhmaxbss(n),rhmax2bss(n),  &
!                   zrhmax2bss(n),zrhdel,wtypebss(n),Usbss(n),Uzbss(n),Uslbss(n),Ubasebss(n),  &
!                   ntr,u1dbss(-ng+1,n),v1dbss(-ng+1,n),z0bss(n),z1bss(n),dudz0bss(n),rotdeg, &
!                   outsoundfiletmp) 
!     IF ( my_rank == 0 ) THEN
!     open(21,file=outsoundfiletmp,form='formatted', status='unknown', position='append')
!     DO k = 1,nzend-1
!      write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), t1dbss(k,numbss+1), q1dbss(k,numbss+1)*1000., &
!        u1dbss(k,numbss+1), v1dbss(k,numbss+1), 1000.*p1dbss(k,numbss+1)**3.508
!     ENDDO
!     close(21)
!     ENDIF
!            DO k = 1,nzend
!               IF ( my_rank == 0 ) write(0,*) 'k,u,v = ',k,u1dbss(k,numbss+1),v1dbss(k,numbss+1)
!            ENDDO
          ENDIF
          
         call snd2x(gd,t1dbss(-ng+1,n),p1dbss(-ng+1,n),q1dbss(-ng+1,n),dzc,dze,zc, &
                  nzend,ng,psfcbss(n),tsfcbss(n),qsfcbss(n),rhmaxbss(n),rhmaxbss(n),rhmax2bss(n),  &
                   zrhmax2bss(n),zrhdel,wtypebss(n),Usbss(n),Uzbss(n),Uslbss(n),Usmvbss(n),Ubasebss(n),  &
                   ntr,u1dbss(-ng+1,n),v1dbss(-ng+1,n),z0bss(n),z1bss(n),dudz0bss(n),rotdeg, &
                   outsoundfile,sndfilebss(n)) 
          
            IF ( ireplaceuvbss(n) >= 1 ) THEN
            DO k = 1,nzend
!               IF ( my_rank == 0 ) write(0,*) 'k,u,v = ',k,u1dbss(k,numbss+1),v1dbss(k,numbss+1)
               u1dbss(k,n) = u1dbss(k,numbss+1) 
               v1dbss(k,n) = v1dbss(k,numbss+1) 
            ENDDO
            ENDIF
            
            IF ( ireplaceqbss(n) >= 1 ) THEN
            DO k = 1,nzend
               q1dbss(k,n) = q1dbss(k,numbss+1) 
            ENDDO
            ENDIF
          
          
        ELSE
          IF ( my_rank == 0 ) THEN
            write(0,*) 'Sorry, sndtypebss = ',sndtypebss(n), ' is not yet supported'
          ENDIF
          call commasmpi_abort()
        ENDIF

       IF ( rhnewbss(n) > 0.0 .and. rhz1bss(n) > 0.0 .and. rhz2bss(n) > rhz1bss(n) ) THEN
       ! specify RH in layer from rhz1bss to rhz2bss
          psfc = psfcbss(n)        
        DO k = 1,nzend
          IF ( zc(k) >= rhz1bss(n) .and. zc(k) <= rhz2bss(n) ) THEN
          pres  = psfc*p1dbss(k,n)**(cp/rd)
          qvs   = 380.*exp(17.27*(p1dbss(k,n)*t1dbss(k,n)-273.16) / (p1dbss(k,n)*t1dbss(k,n)- 35.86)) / pres
          rh    = q1dbss(k,n)/qvs
          q1dbss(k,n) = qvs*Min(rh, rhnewbss(n))
         ! qz(k) = min( rh(k) * qvs, qsfc )
          ENDIF
        ENDDO
       ENDIF
        
        
   ! adjust mean winds
      IF ( umeanbss(n) > -100. ) THEN
        
        sum = 0.0
        height = 0.0 !  ze(nzend)
        DO k = 1,nzend-1
          IF ( zc(k) >= ubasebss(n) .and. zc(k) <= uzbss(n) ) THEN
            sum = sum + u1dbss(k,n)/dzc(k)
            height = height + 1.0/dzc(k)
          ENDIF
        ENDDO
        
        dudz = sum/height - umeanbss(n)
!       write(0,*) 'sum,height,dudz = ',sum,height,dudz

        DO k = 1,nzend-1
!          write(0,*) 'k,u1dbss,  = ',k,u1dbss(k,n),dudz/dzc(k),1.0/dzc(k)
          u1dbss(k,n) = u1dbss(k,n) - dudz
        ENDDO
        
      
      ENDIF

      IF ( vmeanbss(n) > -100. ) THEN

        sum = 0.0
        height = 0.0 !  ze(nzend)
        DO k = 1,nzend-1
          sum = sum + v1dbss(k,n)/dzc(k)
          height = height + 1.0/dzc(k)
        ENDDO
        
        dudz = sum/height - vmeanbss(n)
!       write(0,*) 'sum,height,dudz = ',sum,height,dudz

        DO k = 1,nzend-1
!          write(0,*) 'k,v1dbss,  = ',k,v1dbss(k,n),dudz/dzc(k),1.0/dzc(k)
          v1dbss(k,n) = v1dbss(k,n) - dudz
        ENDDO
      
      
      ENDIF

     IF ( my_rank == 0 ) THEN
     open(21,file=outsoundfile,form='formatted', status='unknown', position='append')
     DO k = 1,nzend-1
      write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), t1dbss(k,n), q1dbss(k,n)*1000., &
        u1dbss(k,n), v1dbss(k,n), 1000.*p1dbss(k,n)**3.508
     ENDDO
     close(21)
     ENDIF
      
      
      ENDDO
   
   ENDIF !} ( .not. bssinit )
   
   
   
  ! check if time is within a BSS time period
   work_to_do = .false.
   
   DO n = 1,numbss

    IF ( time_real > startbss(n) .and. time_real .le. endbss(n) ) THEN
      work_to_do = .true.
      i_current_bss = n
    ENDIF

   ENDDO
   
   IF ( .not. work_to_do ) RETURN
   
   deltat = endbss(i_current_bss) - startbss(i_current_bss)
   deltatstep = time_real - startbss(i_current_bss)
   
   IF ( my_rank == 0 ) THEN
!     write(91,*) 'BSS: time,delt,deltstep : ',time_real,deltat,deltatstep
   
   ENDIF
   

   CALL GET_VARIABLE(gd,'UINIT', u1d1)
   CALL GET_VARIABLE(gd,'VINIT', v1d1)
   CALL GET_VARIABLE(gd,'QVINIT', q1d1)
   CALL GET_VARIABLE(gd,'THINIT', t1d1)
   CALL GET_VARIABLE(gd,'PIINIT', p1d1)
   
   
!       allocate( upert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( vpert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( tpert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( pipert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( qpert(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
       
       ixb = -ng+1
       ixe = itile+ng

       jyb = -ng+1
       jye = jtile+ng

       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       IF ( idoubss(i_current_bss) /= 0 ) THEN
         CALL GET_VARIABLE(gd, 'U', u)
         indx = GET_VARIABLE_INDEX(gd,  'U')
       DO k = kzb,kze
         udel(k) = dt*(u1dbss(k,i_current_bss) - u1dbss(k,i_current_bss-1))/deltat
         u1d1(k) = u1dbss(k,i_current_bss-1) + deltatstep*(u1dbss(k,i_current_bss) - u1dbss(k,i_current_bss-1))/deltat
         gd%var(indx)%base1d(k) = u1d1(k)
!         IF ( my_rank == 0 ) write(91,*) 'udel,uinit = ',k,udel(k),u1d1(k)
       ENDDO
       
       DO k = kzb,kze
       DO j = jyb,jye
       DO i = ixb,ixe
          u(i,j,k) = u(i,j,k) + udel(k)
       ENDDO
       ENDDO
       ENDDO
       
       ENDIF


       IF ( idovbss(i_current_bss) /= 0 ) THEN
         CALL GET_VARIABLE(gd, 'V', v)
         indx = GET_VARIABLE_INDEX(gd,  'V')
       DO k = kzb,kze
         vdel(k) = dt*(v1dbss(k,i_current_bss) - v1dbss(k,i_current_bss-1))/deltat
         v1d1(k) = v1dbss(k,i_current_bss-1) + deltatstep*(v1dbss(k,i_current_bss) - v1dbss(k,i_current_bss-1))/deltat
         gd%var(indx)%base1d(k) = v1d1(k)
!         IF ( my_rank == 0 ) write(91,*) 'vdel,vinit = ',k,vdel(k),v1d1(k)
       ENDDO
       
       DO k = kzb,kze
       DO j = jyb,jye
       DO i = ixb,ixe
          v(i,j,k) = v(i,j,k) + vdel(k)
       ENDDO
       ENDDO
       ENDDO
       
       ENDIF
       
       IF ( idotbss(i_current_bss) /= 0 ) THEN
         CALL GET_VARIABLE(gd, 'TH', th3)
         indx = GET_VARIABLE_INDEX(gd,  'TH')
       DO k = kzb,kze
         tdel(k) = dt*(t1dbss(k,i_current_bss) - t1dbss(k,i_current_bss-1) + tshiftbss(i_current_bss))/deltat
         t1d1(k) = t1dbss(k,i_current_bss-1) +  &
              deltatstep*(t1dbss(k,i_current_bss) - t1dbss(k,i_current_bss-1) + tshiftbss(i_current_bss))/deltat
         gd%var(indx)%base1d(k) = t1d1(k)
        ! IF ( my_rank == 0 ) write(91,*) 'tdel,thinit = ',k,tdel(k),t1d1(k)
       ENDDO
       
       
       DO k = kzb,kze
       DO j = jyb,jye
       DO i = ixb,ixe
          th3(i,j,k) = th3(i,j,k) + tdel(k)
          IF ( my_rank == 0 .and. i == (ixe-ixb)/2 .and. j == (jye-jyb)/2 ) THEN
!            write(91,*) 'tdel,thinit,th3 = ',k,tdel(k),t1d1(k),th3(i,j,k)
          ENDIF
       ENDDO
       ENDDO
       ENDDO
       
       ENDIF


       IF ( idopbss(i_current_bss) /= 0 ) THEN
         CALL GET_VARIABLE(gd, 'PI', pi3)
         indx = GET_VARIABLE_INDEX(gd,  'PI')
       DO k = kzb,kze
         pidel(k) = dt*(p1dbss(k,i_current_bss) - p1dbss(k,i_current_bss-1))/deltat
         p1d1(k) = p1dbss(k,i_current_bss-1) + deltatstep*(p1dbss(k,i_current_bss) - p1dbss(k,i_current_bss-1))/deltat
         gd%var(indx)%base1d(k) = p1d1(k)
       ENDDO
       
       
       DO k = kzb,kze
       DO j = jyb,jye
       DO i = ixb,ixe
!          th3(i,j,k) = th3(i,j,k) + tdel(k)
       ENDDO
       ENDDO
       ENDDO
       
       ENDIF


       IF ( idoqbss(i_current_bss) /= 0 ) THEN
         CALL GET_VARIABLE(gd, 'QV', qv3)
         CALL GET_VARIABLE(gd, 'QC', qc3)
         indx = GET_VARIABLE_INDEX(gd,  'QV')

       DO k = kzb,kze
         qdel(k) = dt*(q1dbss(k,i_current_bss) - q1dbss(k,i_current_bss-1))/deltat
         q1d1(k) = q1dbss(k,i_current_bss-1) + deltatstep*(q1dbss(k,i_current_bss) - q1dbss(k,i_current_bss-1))/deltat
         gd%var(indx)%base1d(k) = q1d1(k)
       ENDDO
       
       DO k = kzb,kze
       DO j = jyb,jye
       DO i = ixb,ixe
          ! not applying within liquid clouds
          IF ( qc3(i,j,k) < 1.e-8 ) qv3(i,j,k) = Max(0.0, qv3(i,j,k) + qdel(k) )
       ENDDO
       ENDDO
       ENDDO
       
       ENDIF
       
!       deallocate ( upert, vpert, tpert, pipert, qpert )
   
   RETURN
   END SUBROUTINE BSS

!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN        \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE BOX     ////////////////////
!
!
! Implements the storm tracking algorithm suggested by Proctor for TASS
!-----------------------------------------------------------------------------
!
! Created by LJW: 07-19-02
!
!-----------------------------------------------------------------------------
!  SUBROUTINE BOX(u,v,w,x1d,y1d,z1d,den,xold,yold,uold,vold,ugrid,vgrid,nx,ny,nz,ng,nx1d,ny1d,nz1d)
!
!		implicit none
!
!! Passed variable declarations
!
!		integer nx, ny, nz, ng, is, js, ks,nx1d,ny1d,nz1d
!		real u(ny,nz,nx)
!		real v(ny,nz,nx)
!		real w(ny,nz,nx)
!
!		real x1d(nx,nx1d)
!		real y1d(ny,ny1d)
!		real z1d(nz,nz1d)
!
!		real den(nz)
!
!
!		real xold,yold,uold,vold,ugrid,vgrid,time
!
!! Local variable declarations
!
!		integer i, j, k, n
!
!		real xmid, ymid
!		real zeta, zeta0, alpha, tau, unew, vnew
!		real*8 xbar, qxbar, ybar, qybar, qbar, q, eps
!
!		parameter( eps   = 1.0e-15 )
!		parameter( zeta0 = 0.001 )
!		parameter( xmid  = 50000., ymid = 55000. )
!		parameter( alpha = 0.75 )
!		parameter( tau   = 200. )
!
!! Init sum variables
!
!		qbar  = 0.0
!		qxbar = 0.0
!		qybar = 0.0
!
!! Compute where storm is above the boundary layer through mid levels
!
!! Compute vorticity
!
!		DO i = 3,nx-3
!		 DO k = 1,nz/2
!			DO j = 3,ny-3
!
!			 zeta = (v(j,  k,i  )-v(j,  k,i-1))*x1d(i,  4) - (u(j,  k,i  )-u(j-1,k,i  ))*y1d(j,  4)
!	 $        + (v(j+1,k,i  )-v(j+1,k,i-1))*x1d(i,  4) - (u(j+1,k,i  )-u(j,  k,i  ))*y1d(j+1,4)
!	 $        + (v(j,  k,i+1)-v(j,  k,i  ))*x1d(i+1,4) - (u(j,  k,i+1)-u(j-1,k,i+1))*y1d(j,  4)
!	 $        + (v(j+1,k,i+1)-v(j+1,k,i  ))*x1d(i+1,4) - (u(j+1,k,i+1)-u(j,  k,i+1))*y1d(j+1,4)
!
!			 zeta = zeta/4.0
!
!			 IF( zeta .gt. zeta0 ) THEN
!
!				q = 1.0 + alog10( zeta0 ) + alog10( zeta ) * den(k) * w(j,k,i)
!
!			 ELSE
!
!				q = den(k) * w(j,k,i)
!
!			 ENDIF
!
!			 qxbar = qxbar + q*q*x1d(i,1)
!			 qybar = qybar + q*q*y1d(j,1)
!			 qbar  = qbar  + q*q
!
!			ENDDO
!		 ENDDO
!
!		ENDDO
!
!		xbar = qxbar / (qbar+eps)
!		ybar = qybar / (qbar+eps)
!
!		unew = uold + (alpha*(xbar-xold) + (1.0-alpha)*(xbar-xmid)) / tau
!		vnew = vold + (alpha*(ybar-yold) + (1.0-alpha)*(ybar-ymid)) / tau
!
!!     unew = uold + ((xbar-xold) + 0.25*(xbar-xmid)) / tau
!!     vnew = vold + ((ybar-yold) + 0.25*(ybar-ymid)) / tau
!
!		write(6,*)
!		print *, 'SUBROUTINE BOX'
!		write(6,*)
!		print *, 'XOLD  = ',xold
!		print *, 'YOLD  = ',yold
!		print *, 'UOLD  = ',uold
!		print *, 'VOLD  = ',vold
!		write(6,*)
!		print *, 'XBAR  = ',xbar
!		print *, 'YBAR  = ',ybar
!		print *, 'UNEW  = ',unew
!		print *, 'VNEW  = ',vnew
!
!		ugrid = unew
!		vgrid = vnew
!
!		xold = xbar
!		yold = ybar
!
!	RETURN
!	END SUBROUTINE BOX
