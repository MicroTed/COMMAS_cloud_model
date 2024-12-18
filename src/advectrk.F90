#define EPSVAL 1.0d-40
#define EPSVAL9 1.0d-40
#define EPSWVAL 1.0d+30

!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN          \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE ADVECT    ////////////////////
!
! ADVECT computes the fluxes for advection for a three-dimensional
!      variable using standard Langrange interpolation
!   
!      Eight different algorithms are available:
!
!        IF( ATYPE .eq. -1 ) CALL ADVECT5N()
!          ==> 5th-order upwind-biased in a fused loop with no collaping of the stencil
!
!        IF( ATYPE .eq. 0 ) CALL ADVECT5()
!          ==> 5th-order upwind-biased
!
!       IF( ATYPE .eq. 1 ) CALL ADVECT_WENO() or WENOC()
!          ==> 5th-order WENO differencing
!
!       IF( ATYPE .eq. 2 ) CALL ADVECT6()
!          ==> 6th order horizontal advection
!          ==> 6th order vertical advection
!
!       IF( ATYPE .eq. 3 ) CALL ADVECT_MONO()
!          ==> 6th order horizontal advection with Leonard 1D monotonic limiter (no spatial filter -- faster than option 5)
!          ==> 5th order vertical advection with Leonard 1D monotonic limiter
!
!       IF( ATYPE .eq. 4 ) CALL ADVECT_MONOF(damph,0.0,.false.)
!          ==> 6th order horizontal advection with Monotonic 6th order filter with damph coefficient
!          ==> 5th order vertical advection
!
!       IF( ATYPE .eq. 5 ) CALL ADVECT_MONOF1(damph)
!          ==> 6th order horizontal advection with Monotonic 6th order filter + Leonard limiter
!          ==> 5th order vertical advection with Leonard 1D monotonic limiter
!
!       IF( ATYPE .eq. 6 ) CALL ADVECT_MONOF(0.0,1.0,.true.)
!          ==> 5th order horizontal advection with with Leonard limiter
!          ==> 5th order vertical advection with Leonard 1D monotonic limiter
!
!       IF( ATYPE .eq. 7 ) CALL ADVECT_MONOF(damph,1.0,.false.)
!          ==> 6th order horizontal advection with Monotonic 6th order filter, 
!              filter coefficent is max of 5th order and damph
!          ==> 5th order vertical advection
!
!       IF( ATYPE .eq. 8 ) CALL ADVECT_MONOF(damph,1.0,.true.)
!          ==> 6th order horizontal advection with Monotonic 6th order filter, 
!              filter coefficent is max of 5th order and damph + Leonard limiter
!          ==> 5th order vertical advection + Leonard limiter
!
!       IF( ATYPE .eq. 9 ) CALL ADVECT_WENOZ() ! not working!!!
!          ==> enhanced 5th-order WENO differencing Borges et al. (2008) J. Comp. Phys.

!       IF( ATYPE .eq. 10 ) CALL WENO7C() 
!          ==> 7th-order WENO differencing Balsara and Shu (2000) J. Comp. Phys.
!
!        IF( ATYPE .eq. 13 ) CALL ADVECT5(.true.) 
!          ==> 5th-order upwind-biased with Leonard monotonic limiter
!
!        IF( ATYPE .eq. 14 ) CALL ADVECT3()
!          ==> 3rd-order upwind-biased
!
!       IF( ATYPE .eq. 15 ) CALL WENOCZ()
!          ==> 5th-order WENO-Z differencing
!
!       IF( ATYPE .eq. 16 ) CALL WENO9C() 
!          ==> 9th-order WENO differencing Balsara and Shu (2000) J. Comp. Phys.

!       IF( ATYPE .eq. 17 ) CALL ADVECT9()
!          ==> 9th order upwind horizontal advection
!          ==> 5th order upwind vertical advection
!

! Standard choice
! 
! For RK iterations 1 & 2:
!                          use ATYPE = 0, e.g, ADVECT5
! For RK iteration  3:
!                          u,v,w,theta use ATYPE = 1, e.g., ADVECT_WENO
!                          All other positive definite variables, use ATYPE = 3, ADVECT_MONO(0.0,0.0,.true)
!
!-----------------------------------------------------------------------------
!
! Created by LJW: 11-01-99
! Latest update:  06-30-06
!
!-----------------------------------------------------------------------------
! DT   = time step - needed for flux limiting procedures
!
! NX   = number of grid points in x
! NY   = number of grid points in y
! NZ   = number of grid points in z
! NG   = number of boundary zones
!
!-----------------------------------------------------------------------------
! NOTES!!
!
! Flux-anti flux formulation used for divergence correction
!
!-----------------------------------------------------------------------------
  SUBROUTINE ADVECT(gd,s,s0,fs,u,v,w,           &
                    gxt,gyt,gzt,rrp1,rrm1,dt, dtlarge, &
                    nx,ny,nz,                &
                    svar,atype,izero0,time,time_real,  &
                    ugrid,vgrid,x_sw_loc,y_sw_loc)

   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE PARAM_MODULE, only: ng, bcx, bcy, damph, ainflo, ainfloqv, ainflom, nocollapse, &
                           vert_adv_scheme, izero_opt,  &
                           zinflo_u_s, zinflo_u_n, zinflo_v_w, zinflo_v_e, zinflo_w_s, &
                           zinflo_w_n, zinflo_w_e, zinflo_w_w
   USE MICRO_MODULE, only: denscale
   USE INIT_MODULE, only: inhom, khomog


   USE COMMASMPI_MODULE

   implicit none 

#ifdef MPI
      INCLUDE "mpif.h"
#endif

  TYPE(GRID) :: gd

   integer, INTENT(IN)    :: nx, ny, nz, atype, izero0
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: rrp1(-ng+1:nz+ng), rrm1(-ng+1:nz+ng)
   real,    INTENT(IN)    :: dt, dtlarge

   TYPE(VARIABLE) :: svar

   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)

   integer :: time
   real    :: time_real
   real    :: ugrid, vgrid
   real    :: x_sw_loc,y_sw_loc

! Local vars

   integer is, js, ks
   integer n, i, j, k, lz
   integer izero
   real    uw, ue, vn, vs, wt, wb, unorm, vnorm
   real    qx1, ux1, qz0, qz1, wz1, qx0, vy0, vy1, qy0, qy1
   real,  allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
   real,  allocatable :: ndgvarns(:,:,:), ndgvarew(:,:,:)

   real    div, divplus, divminus
   real    sijk
   real    sim3, sim2, sim1, sip1, sip2
   real    sjm3, sjm2, sjm1, sjp1, sjp2
   real    skm3, skm2, skm1, skp1, skp2
   real    xmask(-ng+1:nx+ng), ymask(-ng+1:ny+ng), zmask(-ng+1:nz+ng)

   real  :: rrp(-ng+1:nz+ng), rrm(-ng+1:nz+ng)

   integer imn,imx,jmn,jmx,kmn,kmx,imx1,jmx1,kmx1,kmax
   integer imxb, jmxb
   integer i1,i2,i3,i4,i5
   integer j1,j2,j3,j4,j5
   integer k1,k2,k3,k4,k5
         
   real damp1
   real, parameter :: damp0 = 0.0
   real dt1
   character(LEN=10)  :: varname
   
   double precision :: divmult ! = DIVM  ! should be set to 1.0d0 except for special tests
   
   logical  relaxscalar
   logical :: firstcall = .true.
   
   integer bcx1,bcy1
   
   integer :: mpiintin(3),mpiintout(3)

#ifdef MPI
!   logical, parameter :: nocollapse = .true. ! true to turn off stencil collapse at boundaries
#else
!   logical, parameter :: nocollapse = .false. ! true to turn off stencil collapse at boundaries
#endif
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   
   logical :: debug_mpi = .false.

! Parameter flags

! Relaxation at the boundaries parameters
   
!   real, parameter    :: ainflo = 0.5                 !  Leave as 0.0, test value is 0.5
   real               :: ainflow(-ng+1:nz+ng), ainfloe(-ng+1:nz+ng)
   real               :: ainflos(-ng+1:nz+ng), ainflon(-ng+1:nz+ng)

!-----------------------------------------------------------------------------

!   IF ( firstcall .and. my_rank == 0 ) THEN
!     firstcall = .false.
!     write(0,*) 'ADVECTRK: eps for weno is ',EPSVAL
!   ENDIF
   
! Diffusion (only applied in the horizontal) DAMPH IS IN PARAM_MODULE!!

!-----------------------------------------------------------------------------
! Determine staggering
      
   is      = svar%istag
   js      = svar%jstag
   ks      = svar%kstag

   varname = svar%name
   
   bcy1 = 1
   
   rrp(:) = rrp1(:)
   rrm(:) = rrm1(:)
   
   IF ( ( denscale == 2 .and. (varname(1:1) == 'C' .or. varname(1:1) == 'Q' .or. varname(1:1) == 'V') ) .or.  &
        ( denscale == -2 .and. ( svar%dyntype .eq. 0 ) )  ) THEN
!        ( denscale == -2 .and. (varname(1:1) == 'C' .or. varname(1:1) == 'V' .or. varname(1:1) == 'Z' .or. varname(1:2) == 'SC') )  ) THEN
      divmult = 0.0d0
   ELSE
      divmult = 1.0d0
   ENDIF

!   IF ( ( (denscale == -1 .or. denscale == -2) .and. ( varname(1:1) == 'C' .or. varname(1:1) == 'V' .or. varname(1:2) == 'SC') ) .or. &
   IF ( ( (denscale == -1 .or. denscale == -2) .and. ( svar%dyntype .eq. 0 ) ) .or. &
       ( ( denscale  == 3 .and. (varname(1:1) == 'C' ) ) ) .or. denscale == 4 ) THEN
     rrp(:) = 1.0
     rrm(:) = 1.0
   ENDIF

!-----------------------------------------------------------------------------

   if (debug_mpi)  write(0,"('ADVECTRK: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank

  IF ( ny .le. 2 .and. svar%name .eq. 'V' ) RETURN
  IF ( nx .le. 2 .and. svar%name .eq. 'U' ) RETURN

#ifdef MPI
!  if ( bcx .eq. 2 .or. bcy .eq. 2) then
!    write(luno,*) 'ADVECTRK: MPI not tested for periodic bc yet'
!    write(luno,*) 'ADVECTRK: stopping run ...'
!    STOP
!  endif
#endif

!-----------------------------------------------------------------------------
! Set bounds for advection if izero=1 to avoid advecting zeroes.

  ! debug:
  IF ( izero_opt == 0 ) THEN
    izero = 0
  ELSE
    izero = izero0
  ENDIF

  IF ( bcx == 2 .or. bcy == 2 ) izero = 0

#ifdef MPI
!   izero = 0
#endif

  IF ( izero .ge. 1 ) THEN
    CALL cld_cpu('ADVECT-IZERO1')
#ifdef MPI
        imn = nxend-is
        imx = 1
        jmn = nyend-js
        jmx = 1
        kmn = nzend-ks
        kmx = 1
       ixb=-ng+1
       if(ixbeg == nxbeg ) ixb = 1
       ixe=itile+ng-1
       if(ixend.eq.nxend) ixe=ixend-ixbeg-is
       
       jyb=-ng+1
       if ( jybeg == nybeg ) jyb = 1
       jye=jtile+ng-1
       if(jyend.eq.nyend) jye=jyend-jybeg-js

       kzb=1
       kze=ktile
       if(kzend.eq.nzend) kze=kzend-kzbeg-ks
       
!       write(0,*) 'check for nonzero in ',varname,my_rank
       do k=kzb,kze ; do j=jyb,jye ; do i=ixb,ixe
          IF ( s(i,j,k) .ne. 0.0 .or. s0(i,j,k) .ne. 0.0 ) THEN
            imn = Min(imn,Max(1,ixbeg+i-1))
            imx = Max(imx,Min(ixbeg+i-1,nxend-1))
            jmn = Min(jmn,Max(1,jybeg+j-1))
            jmx = Max(jmx,Min(jybeg-1+j,nyend-1))

            kmn = Min(kmn,Max(1,kzbeg-1+k))
            kmx = Max(kmx,Min(kzbeg-1+k,nzend-1))
          ENDIF
       enddo ; enddo ; enddo
!       write(0,*) 'done check',varname,my_rank,imn,imx,jmn,jmx,kmn,kmx

      mpiintin(1) = imn
      mpiintin(2) = jmn
      mpiintin(3) = kmn

      CALL MPI_AllReduce(mpiintin, mpiintout, 3, MPI_INTEGER, MPI_MIN, my_comm, mpi_error_code)
      
      imn = mpiintout(1)
      jmn = mpiintout(2)
      kmn = mpiintout(3)
      
      mpiintin(1) = imx
      mpiintin(2) = jmx
      mpiintin(3) = kmx

      CALL MPI_AllReduce(mpiintin, mpiintout, 3, MPI_INTEGER, MPI_MAX, my_comm, mpi_error_code)

      imx = mpiintout(1)
      jmx = mpiintout(2)
      kmx = mpiintout(3)

! Bounds for the flux calculations

    imn = MAX(imn - ng, 1+is)
    imx = MIN(imx + ng + 1,nxend-1+is)
    jmn = MAX(jmn - ng, 1+js)
    jmx = MIN(jmx + ng + 1,nyend-1+js)
    kmn = MAX(kmn - ng, 1+ks)
    kmx = MIN(kmx + ng + 1,nzend-1)

!     imn = 1+is
!     imx = nxend-1+is
!     jmn = 1+js
!     jmx = nyend-1+js
!     kmn = 1+ks
!     kmx = nzend-1 
!     imx1 = nxend-1
!     jmx1 = nyend-1
!     kmx1 = nzend-1
    
    imx1 = Min(imx,nxend-1)
    jmx1 = Min(jmx,nyend-1)
    kmx1 = Min(kmx,nzend-1)

    CALL cld_cpu('ADVECT-IZERO1')

    IF ( imn .gt. imx ) RETURN
    
    IF ( imn > ixend .or. imx < ixbeg .or. jmn > jyend .or. jmx < jybeg ) RETURN

#else
    imn = nx-1
    imx = 1
    jmn = ny-1
    jmx = 1
    kmn = nz-1
    kmx = 1

    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1

       IF ( s(i,j,k) .ne. 0.0 .or. s0(i,j,k) .ne. 0.0 ) THEN
       
         imn = Min(imn,Max(1,i))
         imx = Max(imx,Min(i,nx-1))
         jmn = Min(jmn,Max(1,j))
         jmx = Max(jmx,Min(j,ny-1))
         kmn = Min(kmn,Max(1,k))
         kmx = Max(kmx,Min(k,nz-1))
            
       ENDIF

      ENDDO
     ENDDO
    ENDDO

! Bounds for the flux calculations

    imn = MAX(imn - ng, 1+is)
    imx = MIN(imx + ng + 1,nx-1+is)
    jmn = MAX(jmn - ng, 1+js)
    jmx = MIN(jmx + ng + 1,ny-1+js)
    kmn = MAX(kmn - ng, 1)
    kmx = MIN(kmx + ng + 1,nz-1)

    CALL cld_cpu('ADVECT-IZERO1')

    IF ( imn .gt. imx ) RETURN

    imx1 = Min(imx,nx-1)
    jmx1 = Min(jmx,ny-1)
    kmx1 = Min(kmx,nz-1)
#endif

! Check if upper bound is below nz-1; if so, need to reduce kmx1 so that
! the flux at kmx1+1 is defined (i.e., at kmx)

    IF ( kmx .lt. nz-1 ) kmx1 = kmx-1

  ELSE ! izero is set to force the whole domain (e.g., for qv, u, v, etc.)

#ifdef MPI
    imn = 1+is
    imx = nxend-1+is
    jmn = 1+js
    jmx = nyend-1+js
    kmn = 1+ks
!    kmn = 1
    kmx = nzend-1 
    imx1 = nxend-1
    jmx1 = nyend-1
    kmx1 = nzend-1
    IF ( bcy == 2 ) THEN
      jmn = 1
!      jmx1 = nyend
      jmx = nyend+js
    ENDIF
    IF ( bcx == 2 ) THEN
      imn = 1
!      imx1 = nxend
      imx = nxend+is
    ENDIF
#else
    imn = 1+is
    imx = nx-1+is
    jmn = 1+js
    jmx = ny-1+js
!    kmn = 1+ks
    kmn = 1
    kmx = nz-1
    imx1 = nx-1
    jmx1 = ny-1
    kmx1 = nz-1
#endif

  ENDIF

  IF ( ny .eq. 2 ) THEN
   jmn = 1
   jmx = 1
  ENDIF

#ifdef MPI
    i5 = nxend-5+is
    i4 = nxend-4+is
    i3 = nxend-3+is
    i2 = nxend-2+is
    i1 = nxend-1+is

    j5 = nyend-5+js
    j4 = nyend-4+js
    j3 = nyend-3+js
    j2 = nyend-2+js
    j1 = nyend-1+js

    k5 = nzend-5+ks
    k4 = nzend-4+ks
    k3 = nzend-3+ks
    k2 = nzend-2+ks
    k1 = nzend-1+ks
    
    imxb = imx
    jmxb = jmx
#else
    i5 = nx-5+is
    i4 = nx-4+is
    i3 = nx-3+is
    i2 = nx-2+is
    i1 = nx-1+is

    j5 = ny-5+js
    j4 = ny-4+js
    j3 = ny-3+js
    j2 = ny-2+js
    j1 = ny-1+js

    k5 = nz-5+ks
    k4 = nz-4+ks
    k3 = nz-3+ks
    k2 = nz-2+ks
    k1 = nz-1

    imxb = imx
    jmxb = jmx
#endif

  CALL cld_cpu('ADVECT-SETBC')

!----------------------------------------------------------------------
! Boundary conditions
!
!  IF( bcx .ge. 0 ) THEN    ! zero gradient
!
!   DO n = 1,ng
!    s( 1-n      ,1:ny,1:nz) = s0(1      ,1:ny,1:nz)
!    s( nx-1+n+is,1:ny,1:nz) = s0(nx-1+is,1:ny,1:nz)
!   ENDDO
!
!  ENDIF

!! x-direction

  IF ( bcx .eq. 1 .or. bcx .ge. 3 ) THEN   ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug: no need to fill ghost zones
   DO n = 1,ng
    s( 1-n      ,1:ny,1:nz) = s(1      ,1:ny,1:nz)
    s( nx-1+n+is,1:ny,1:nz) = s(nx-1+is,1:ny,1:nz)
   ENDDO
#else
   IF ( ixbeg .eq. nxbeg ) THEN
   DO n = 1,ng
    s( 1-n      ,1:ny,1:nz) = s(1      ,1:ny,1:nz)
   ENDDO
   ENDIF
   
   IF ( ixend .eq. nxend ) THEN
   DO n = 1,ng
    s( nx-1+n+is,1:ny,1:nz) = s(nx-1+is,1:ny,1:nz)
   ENDDO
   ENDIF
#endif

  ELSEIF( bcx .eq. 2 ) THEN    ! periodic

   DO n = 1,ng
#ifdef MPI
!    kzb = 1
!    kze = ktile
!    if (kzend .eq. nzend) kze = kzend-kzbeg+1
!
!    jyb = 1
!    jye = jtile
!    if (jyend .eq. nyend) jye = jyend-jybeg+1
!
!    do k = kzb,kze ; do j = jyb, jye
!     s( nxend-1+n+is,j,k) = s(n+is   ,j,k)
!     s( 1-n         ,j,k) = s(nxend-n,j,k)
!    enddo ; enddo
#else
     s( nx-1+n+is,1:ny,1:nz) = s(n+is    ,1:ny,1:nz)
     s( 1-n      ,1:ny,1:nz) = s(nx-n,1:ny,1:nz)
#endif
   ENDDO
#ifdef MPI
    IF ( imn .eq. 1 + is ) THEN
      imn = 1
      imx = nxend-1+is
    ENDIF
    IF ( imx .eq. nxend-1+is ) THEN
      imxb = nxend + is
      imn  = 1
    ENDIF
#else
    IF ( imn .eq. 1 + is ) THEN
      imn = 1
      imx = nx-1+is
    ENDIF
    IF ( imx .eq. nx-1+is ) THEN
      imxb = nx + is
      imn  = 1
    ENDIF
#endif

  ELSEIF ( bcx .eq. 0 ) THEN  ! mirror plane (solid wall) on western boundary, open (zero gradient) on east

   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg+1

    do k = kzb,kze ; do j = jyb, jye
      s(    1-n   ,j,k) = s(   n   ,j,k)
!      s( nx-1+n+is,j,k) = s(nx-n+is,j,k)
      s( nx-1+n+is,j,k) = s(nx-1+is,j,k)
    enddo ; enddo
#else
    s( 1-n      ,1:ny,1:nz) = s(n      ,1:ny,1:nz)
!    s( nx-1+n+is,1:ny,1:nz) = s(nx-n+is,1:ny,1:nz)
    s( nx-1+n+is,1:ny,1:nz) = s(nx-1+is,1:ny,1:nz)
#endif
   ENDDO

  ENDIF !bcx

!! y-direction

  IF ( bcy .eq. 1 .or. bcy .ge. 3 ) THEN    ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug:  no need to fill ghost zones
   DO n = 1,ng
    s(1:nx, 1-n     ,1:nz) = s(1:nx,1      ,1:nz)
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,ny-1+js,1:nz) 
   ENDDO
#else
   IF ( jybeg .eq. nybeg ) THEN
   DO n = 1,ng
    s(1:nx, 1-n     ,1:nz) = s(1:nx,1      ,1:nz)
   ENDDO
   ENDIF

   IF ( jyend .eq. nyend ) THEN
   DO n = 1,ng
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,ny-1+js,1:nz) 
   ENDDO
   ENDIF

#endif

  ELSEIF( bcy .eq. 2 ) THEN    ! periodic

!  IF ( jyend == nyend .and. svar%name .eq. 'CCCN' ) THEN
!    k = 3
!    write(91,*) 'CCCN at N boundary'
!    DO j = ny-2,ny+2
!     write(91,*) 'j = ',j,'s= ', ( s(i,j,k), i=nx/2,nx-1)
!    ENDDO
!  
!  ENDIF
!
!  IF ( jybeg == nybeg .and. svar%name .eq. 'CCCN' ) THEN
!    k = 3
!    write(91,*) 'CCCN at S boundary'
!    DO j = -2,2
!     write(91,*) 'j = ',j,'s= ', ( s(i,j,k), i=nx/2,nx-1)
!    ENDDO
!  
!  ENDIF

   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    ixb = 1
    ixe = jtile
    if (ixend .eq. nxend) ixe = ixend-ixbeg+1
!mpi need processor communication here?  No, already have done it??
!    do k = kzb,kze ; do i = ixb,ixe
!     s(i,nyend-1+n+js,k) = s(i,n+js,k) 
!     s(i,1-n      ,k) = s(i,nyend-n,k)
!    enddo ; enddo
#else
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,n+js,1:nz) 
    s(1:nx,1-n      ,1:nz) = s(1:nx,ny-n,1:nz)
#endif
   ENDDO

#ifdef MPI
      jmn = 1
      jmx = nyend+js
    IF ( jyend == nyend ) THEN
      jmxb = nyend + js
!      jmx = nxend - 1 + js
    ENDIF
#else
    IF ( jmn .eq. 1 + js ) THEN
      jmn = 1
      jmx = ny-1+js
    ENDIF
    IF ( jmx .eq. ny-1+js ) THEN 
      jmxb = ny + js
      jmn = 1
    ENDIF
#endif
    
  ELSEIF ( bcy .eq. 0 ) THEN    ! mirror (solid wall) on north boundary, open on south:
   DO n = 1,ng
!    s(1:nx, 1-n     ,1:nz) = s(1:nx,n      ,1:nz)
    s(1:nx, 1-n     ,1:nz) = s(1:nx,1      ,1:nz)
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,ny-n+js,1:nz) 
   ENDDO

  ENDIF ! bcy

!
!  IF( bcy .ge. 0 ) THEN    ! zero gradient
!
!   DO n = 1,ng
!    s(1:nx,1-n,      1:nz)      = s0(1:nx,1,      1:nz)
!    s(1:nx,ny-1+n+js,1:nz)      = s0(1:nx,ny-1+js,1:nz)
!   ENDDO
!
!  ENDIF
!

! Upper and lower bounary conditions

   IF( ks .eq. 0 ) THEN    ! reflective
 
    DO n = 1,ng
     s (1:nx,1:ny,   1-n) = s0(1:nx,1:ny,   n)
     s (1:nx,1:ny,nz-1+n) = s0(1:nx,1:ny,nz-n)

! Need for SL WENO

     s0(1:nx,1:ny,   1-n) = s0(1:nx,1:ny,   n)
     s0(1:nx,1:ny,nz-1+n) = s0(1:nx,1:ny,nz-n)
    ENDDO
 
   ENDIF
 
   IF( ks .eq. 1 ) THEN    ! anti - reflective
 
    DO n = 1,ng
     s (1:nx,1:ny, 1-n) = -s0(1:nx,1:ny, n+1)
     s (1:nx,1:ny,nz+n) = -s0(1:nx,1:ny,nz-n)

! Need for SL WENO

     s0(1:nx,1:ny, 1-n) = -s0(1:nx,1:ny, n+1)
     s0(1:nx,1:ny,nz+n) = -s0(1:nx,1:ny,nz-n)
    ENDDO
 
   ENDIF

  CALL cld_cpu('ADVECT-SETBC')

!--------------------------------------------------------------------------
! DO THE VARIOUS TYPES OF ADVECTION


    CALL cld_cpu('ADVECT-IZERO2')

    allocate ( ndgvarns(nx,nz,2) )
    allocate ( ndgvarew(ny,nz,2) )

  IF ( .not. ( ATYPE .eq. 0 .or. ATYPE .eq. 1 .or. ATYPE .eq. 51 .or. & 
               ATYPE .eq. 10 .or. ATYPE .eq. 16 ) ) THEN

#ifdef MPI
  allocate (  fx(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  allocate (  fy(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  allocate (  fz(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  
  IF ( izero >= 1 .and. atype >= 2 .and. ( atype <= 9 )  ) THEN ! need to have the imx+1, etc., planes set to zero if not filling the whole array
  i = imx1-ixbeg+1
  fx(i:i+1,: ,: ) = 0.0 
  j = jmx1-jybeg+1
  fy(: ,j:j+1,: ) = 0.0
  k = kmx1-kzbeg+1
  fz(: ,: ,k:k+1) = 0.0 
  ENDIF
  
  IF ( izero >= 1 .and. atype == 12 ) THEN
   fx(:,:,:) = 0.0
   fy(:,:,:) = 0.0
   fz(:,:,:) = 0.0
  ENDIF
  
!  write(0,*) 'advectrk: rank,ny,ny+ng = ',my_rank,ny,ny+ng,jmx,jmxb
  
#else
  allocate (  fx(0:nx+1,ny,nz) )
  allocate (  fy(nx,0:ny+1,nz) )
  allocate (  fz(nx,ny,nz) )

! #ifdef __ia64__

  fx(0:1,:  ,:) = 0.0 
  fy(:  ,0:1,:) = 0.0 
  fz(:  ,:  ,1) = 0.0 

  fx(nx-1:nx,: ,: ) = 0.0 
  fy(: ,ny-1:ny,: ) = 0.0 
  fz(: ,: ,nz-1:nz) = 0.0 

  fz(:,:,kmx1+1) = 0.0 ! for non-mpi, set z-flux to zero at one point above kmx1

! #endif

#endif

  IF ( nx .eq. 2 ) fx(:,:,:) = 0.0
  IF ( ny .eq. 2 ) fy(:,:,:) = 0.0

!   fz(:,:,:) = -1.e32
!   fx(:,:,:) = -1.e32
!   fy(:,:,:) = -1.e32

   ENDIF
   
!   s(:,:,nz+ks:nz+ng) = -1.e32
!   s0(:,:,nz+ks:nz+ng) = -1.e32

!#ifdef __ia64__
!   fz(:,:,:) = 0.0 
!   fx(:,:,:) = 0.0 
!   fy(:,:,:) = 0.0 
!#endif



! Changed xmask,ymask to kill fluxes at sides (imx,jmx) and added zmax for top (kmx1+1) of subdomain box
  xmask(:) = 1.0 
  IF ( bcx .ne. 2 ) THEN
#ifdef MPI
    if (ixbeg .eq. nxbeg) xmask(-ng+1:1) = 0.0
    if (ixend .eq. nxend) then
      ixb = imx-ixbeg+1
      ixe = ixend-ixbeg+1+ng

      xmask(ixb:ixe) = 0.0
    endif
#else
    xmask(-ng+1:1) = 0.0 ; xmask(imx:nx+ng) = 0.0
#endif
  ELSE
#ifndef MPI
    imx1 = imx1 + is
!    fx(nx+is,1:ny-1+js,1:nz-1+ks) = fx(1,1:ny-1+js,1:nz-1+ks)
#endif
  ENDIF

  ymask(:) = 1.0 
  IF ( bcy .ne. 2 ) THEN
#ifdef MPI
!    if (jybeg .eq. nybeg)  ymask(-ng+1:1) = 0.0
    if (jyend .eq. nyend) then
      jyb = jmx-jybeg+1
      jye = jyend-jybeg+1+ng

      ymask(jyb:jye) = 0.0
    endif
#else
    ymask(-ng+1:1) = 0.0 ; ymask(jmx:ny+ng) = 0.0
#endif
  ENDIF

#ifdef MPI
  IF ( bcy .ne. 2 )  THEN
   if (jybeg .eq. nybeg)  ymask(-ng+1:1) = 0.0
  ENDIF
#endif

  zmask(:) = 1.0
#ifdef MPI
  if (kzbeg .eq. nzbeg)  zmask(-ng+1:1) = 0.0
  if (kzend .eq. nzend) then
   kzb = kmx1-kzbeg+2
   kze = kzend-kzbeg+1+ng

   zmask(kzb:kze) = 0.0
  endif
#else
  zmask(kmx1+1:nz+ng) = 0.0 ; zmask(-ng+1:1) = 0.0
#endif

  IF ( ny .le. 2 ) ymask(:) = 0.0
  IF ( nx .le. 2 ) xmask(:) = 0.0


    CALL cld_cpu('ADVECT-IZERO2')
  damp1 = damph

  IF( ATYPE .eq.-1 ) CALL ADVECT5N()
  IF( ATYPE .eq. 0 ) CALL ADVECT5()

  IF ( ATYPE .eq. 51 ) CALL WENOC() ! WENO5C() !  WENOC()
  IF ( ATYPE .eq. 1 ) CALL WENO5C() !  WENOC()
! There is something weird about ADVECT_WENOCF -- gives same answer as WENOC on intel Mac with -O0,
! but takes 3 times longer.  At higher opt it gives different results than WENOC, which gives about
! same for both -O0 and -O3.
!      CALL ADVECT_WENOCF(s,s0,fs,u,v,w,rrp,rrm,dt,nx,ny,nz,      &
!                         imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1, &
!                         imxb, jmxb, &
!                         i1,i2,i3,j1,j2,j3,k1,k2,k3,is,js,ks,  &
!                         gxt,gyt,gzt,xmask,ymask,varname)

  IF ( ATYPE .eq. 02 ) CALL ADVECT6()
!  IF ( ATYPE .eq. 03 ) CALL ADVECT_MONO()
  IF ( ATYPE .eq. 03 ) CALL ADVECT_MONOF1(damp0)
  IF ( ATYPE .eq. 04 ) CALL ADVECT_MONOF(damp1,0.0,.false.)
  IF ( ATYPE .eq. 05 ) CALL ADVECT_MONOF1(damp1)
  IF ( ATYPE .eq. 06 ) CALL ADVECT_MONOF(0.0,1.0,.true.)
  IF ( ATYPE .eq. 07 ) CALL ADVECT_MONOF(damp1,1.0,.false.)
  IF ( ATYPE .eq. 08 ) CALL ADVECT_MONOF(damp1,1.0,.true.)
  IF ( ATYPE .eq. 09 ) CALL ADVECT_WENOZ()  
  IF ( ATYPE .eq. 10 ) CALL WENO7C()
  IF ( ATYPE .eq. 11 ) CALL ADVECT_WENO()  ! old, non-cache-friendly version
  IF ( ATYPE .eq. 12 ) CALL ADVECT_SLWENO()
  IF ( ATYPE .eq. 13 ) CALL ADVECT5(.true.) ! turn on 1D Leonard filter
  IF ( ATYPE .eq. 14 ) CALL ADVECT3()       ! 3rd order upwind advection
  IF ( ATYPE .eq. 15 ) CALL WENOCZ()
  IF ( ATYPE .eq. 16 ) CALL WENO9C()
  IF ( ATYPE .eq. 17 ) CALL ADVECT9()

!-----------------------------------------------------------------------------
! Set boundary masks

! xmask(:) = 1.0 ; xmask(1) = 0.0 ; xmask(nx-1+is:nx) = 0.0
! ymask(:) = 1.0 ; ymask(1) = 0.0 ; ymask(ny-1+js:ny) = 0.0

  IF ( .not. ( ATYPE .eq. 0 .or. ATYPE .eq. 1 .or. ATYPE .eq. 10 .or. ATYPE == 13 .or. &
               atype == 14 .or. atype == 15 .or. atype == 16 .or. atype == 17  .or. atype == 51 ) ) THEN
  
   CALL cld_cpu('ADVECT-FLUX')

  IF (myprock == 1 ) fz(:,:, 1) = 0.0 
  IF (myprock == nprock ) fz(:,:,nz) = 0.0 

#ifdef MPI
  ixb = 1
  ixe = itile
!  if (ixbeg .le. imn)  ixb = imn -ixbeg+1
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1
  
!  write(0,*) 'advectrk: jye,my_rank = ',jye,my_rank

  kzb = 1
  kze = ktile
  if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
  if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,ue,vs,vn,wb,wt,div, divplus, divminus)
  do k = kzb,kze
   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,ue,vs,vn,wb,wt,div, divplus, divminus)
  DO k = kmn,kmx1 
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
        uw        = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
        ue        = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
        vs        = 0.5*(v(i  ,j,  k  ) + v(i-is  ,j-js,  k-ks  ))
        vn        = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))
        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)
        wt        = 0.5*(w(i  ,j,  k+1) + w(i-is  ,j-js,  k-ks+1))*rrp(k)
        div       = (gxt(i,3+is)*(ue-uw) + gyt(j,3+js)*(vn-vs) + gzt(k,3+ks)*(wt-wb))

    IF ( .not. ( ATYPE .eq. 11 .or. ATYPE .eq. 12 ) ) THEN
! Old version:
      fs(i,j,k) =  -            (fz(i  ,j  ,k+1)*rrp(k) - fz(i,j,k)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy(i  ,j+1,k  )        - fy(i,j,k)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx(i+1,j  ,k  )        - fx(i,j,k)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div * divmult

    ELSE

! Eq 4.3 & 4.4 from Shchepetkin & Williams, 2009
! Needed because of the SLWENO, plus I think its more accurate in time.

        divplus  = 1.0 - 0.5*dt*div
        divminus = 1.0 + 0.5*dt*div
 
        fs(i,j,k) =  (s0(i,j,k)*divminus - dt*(          (wt*fz(i  ,j  ,k+1) - wb*fz(i,j,k))*gzt(k,3+ks) &
                                              + ymask(j)*(vn*fy(i  ,j+1,k  ) - vs*fy(i,j,k))*gyt(j,3+js) &
                                              + xmask(i)*(ue*fx(i+1,j  ,k  ) - uw*fx(i,j,k))*gxt(i,3+is) ) ) / divplus
 
        fs(i,j,k) = (fs(i,j,k) - s0(i,j,k)) / dt
      
     ENDIF

     ENDDO
    ENDDO
   ENDDO
   
    CALL cld_cpu('ADVECT-FLUX')
   
   ENDIF ! ( ATYPE .ne. 0 ) 

!-----------------------------------------------------------------------------
! Add in outflow boundary conditions using the WRF split formulation at the boundaries
!     for the flux and also have code to set up boundary relaxation terms


    CALL cld_cpu('ADVECT-INFLOW')

  ainflow(:) = 0.0
  ainfloe(:) = 0.0
  ainflos(:) = 0.0
  ainflon(:) = 0.0

  ndgvarew(:,:,:) = 0.0
  ndgvarns(:,:,:) = 0.0
  
!
! MSB 7/5/11 add nudging to boundary u accelerations for density current simulations
!

  IF ((bcx.ne.1).and.(bcy.ne.1) .and. (bcx.ne.3).and.(bcy.ne.3) .and. (bcx.ne.4).and.(bcy.ne.4)) then

  IF( ixbeg.eq.nxbeg.or.ixend.eq.nxend.or.jybeg.eq.nybeg.or.jyend.eq.nyend ) THEN  
!  CALL SET_BASE(gd,varname,ndgvarew,ndgvarns,nx,ny,nz,time)
  ENDIF
  
  ENDIF


  DO k = 1,nz
    IF ( nx .gt. 2 ) THEN
    ndgvarew(1:ny,k,1:2) = 0.0 ! svar%base1d(k)
    ENDIF
    IF ( ny .gt. 2 ) THEN
    ndgvarns(1:nx,k,1:2) = 0.0 ! svar%base1d(k)
    ENDIF
  ENDDO
  
! here, turn on inflow nudging on scalars that have base states of all zero.
  IF (  ainflo .gt. 0.0 .and. is == 0 .and. js == 0 .and. ks == 0 ) THEN

   relaxscalar = .true.

   DO k = 1,nz-1
      relaxscalar = (relaxscalar .and. svar%base1d(k) == 0.0)
   ENDDO
  
   IF ( relaxscalar ) THEN
     IF ( bcx /= 2 ) THEN ! do not do this for periodic EW 
       ainflow(:) = ainflo
       ainfloe(:) = ainflo
     ENDIF
     IF ( bcy /= 2 ) THEN ! do not do this for periodic NS
       ainflos(:) = ainflo
       ainflon(:) = ainflo
     ENDIF

     ndgvarew(:,:,:) = 0.0
     IF ( ny .gt. 2 ) THEN
     ndgvarns(:,:,:) = 0.0
     ENDIF

   ENDIF
  
  ENDIF


  IF (  ainfloqv > 0.0 .and. varname == 'QV' ) THEN
    
     IF ( bcx /= 2 ) THEN ! do not do this for periodic EW 
       ainflow(:) = ainfloqv
       ainfloe(:) = ainfloqv
     ENDIF
     IF ( bcy /= 2 ) THEN ! do not do this for periodic NS
       ainflos(:) = ainfloqv
       ainflon(:) = ainfloqv
     ENDIF

     DO k = 1,nz-1
     ! write(0,*) 'k, base1d(qv) = ',k,svar%base1d(k),my_rank
     ndgvarew(:,k,:) = svar%base1d(k)
     IF ( ny .gt. 2 ) THEN
     ndgvarns(:,k,:) = svar%base1d(k)
     ENDIF
     ENDDO
    
  ENDIF
  
  IF ( ainflom > 0.0 .and. ( is == 1 .or. js == 1 ) ) THEN
  
    IF (  zinflo_u_s > 0.0 .and. is == 1 ) THEN
      DO k = 1,nz-1
        IF ( zinflo_u_s > gzt(k,1) ) THEN
          ainflos(k) = ainflom
          ndgvarns(:,k,1) = svar%base1d(k) - ugrid
        ENDIF
      ENDDO
    ENDIF

  
  ENDIF

  !  zinflo_w_s to nudge W on south inflow boundary
  IF ( ainflom > 0.0 .and.  ks == 1 .and. zinflo_w_s > 0.0 ) THEN
  
      DO k = 1,nz-1
        IF ( zinflo_w_s > gzt(k,1) ) THEN
          ainflos(k) = ainflom
          ndgvarns(:,k,1) = svar%base1d(k) 
        ENDIF
      ENDDO
 
  ENDIF

  !  zinflo_w_n to nudge W on north inflow boundary
  IF ( ainflom > 0.0 .and.  ks == 1 .and. zinflo_w_n > 0.0 ) THEN
  
      DO k = 1,nz-1
        IF ( zinflo_w_n > gzt(k,1) ) THEN
          ainflon(k) = ainflom
          ndgvarns(:,k,2) = svar%base1d(k) 
        ENDIF
      ENDDO
 
  ENDIF

  !  zinflo_w_w to nudge W on west inflow boundary
  IF ( ainflom > 0.0 .and.  ks == 1 .and. zinflo_w_w > 0.0 ) THEN
  
      DO k = 1,nz-1
        IF ( zinflo_w_s > gzt(k,1) ) THEN
          ainflow(k) = ainflom
          ndgvarew(:,k,1) = svar%base1d(k) 
        ENDIF
      ENDDO
 
  ENDIF

  !  zinflo_w_e to nudge W on east inflow boundary
  IF ( ainflom > 0.0 .and.  ks == 1 .and. zinflo_w_e > 0.0 ) THEN
  
      DO k = 1,nz-1
        IF ( zinflo_w_s > gzt(k,1) ) THEN
          ainfloe(k) = ainflom
          ndgvarew(:,k,2) = svar%base1d(k) 
        ENDIF
      ENDDO
 
  ENDIF





  IF ( varname .eq. 'TH' .or. varname .eq. 'QV' .or. varname.eq.'U' .or. varname.eq.'V') THEN

  IF (ixbeg.eq.nxbeg .or. ixend.eq.nxend .or. jybeg.eq.nybeg .or. jyend.eq.nyend ) THEN  !{

  IF ( bcx == 3 .or. bcy == 3 ) THEN !{

  ndgvarew(:,:,:)=0.0
  ndgvarns(:,:,:)=0.0

  kmax = 0
   CALL cld_cpu('TDLBC')
   CALL INT_TDLBC(gd,varname,ndgvarew,ndgvarns,time,nx,ny,nz,kmax)
   CALL cld_cpu('TDLBC')
 
!
!.... CLZ (3-26-13): experiment with ainflo as function of height using new index lz
!

  do lz = 1, nz
!
!.... CLZ (4-30-13): experiment with ainflo values for scalars and parallel wind components
!
  if (lz .lt. khomog) then
  ainflow(lz) = 0.5
  ainfloe(lz) = 0.5
  ainflos(lz) = 0.0
!  ainflos(lz) = 0.5
  ainflon(lz) = 0.0

  elseif (lz .ge. khomog) then
  ainflow(lz) = 0.0
  ainfloe(lz) = 0.0
  ainflos(lz) = 0.0
  ainflon(lz) = 0.0

  endif
  enddo

!
!.... CLZ (3-6-13): moving grid in fixed mesoscale domain
!

  ELSEIF ( (bcx .eq. 4) .or. (bcy .eq. 4) ) THEN ! 

!
! MSB 8/22/11 add BCs for moving grid within larger mesoscale domain
!  bc=4 analytic mesoscale grid

  ndgvarew(:,:,:) = 0.0
  ndgvarns(:,:,:) = 0.0
  
   CALL cld_cpu('TDLBC')
  CALL INT_TDLBC_MOV(gd,varname,ndgvarew,ndgvarns,   &
                     time,nx,ny,nz,kmax,ugrid,vgrid,x_sw_loc,y_sw_loc)
   CALL cld_cpu('TDLBC')

!
!.... CLZ (3-26-13): experiment with ainflo as function of height using new index lz
!

  do lz = 1, nz

!
!.... CLZ (5-1-13): experiment with CUSTOM ainflo values for scalars and wind components
!


  IF (lz .lt. khomog) THEN

!.... west face
  IF (varname.eq.'TH'.or.varname.eq.'QV') THEN
  ainflow(lz) = 0.5
  ELSEIF (varname.eq.'U') THEN
  ainflow(lz) = 0.0
  ELSEIF (varname.eq.'V') THEN
  ainflow(lz) = 0.0
  ENDIF

!.... east face
  IF (varname.eq.'TH'.or.varname.eq.'QV') THEN
  ainfloe(lz) = 0.5
  ELSEIF (varname.eq.'U') THEN
  ainfloe(lz) = 0.0
  ELSEIF (varname.eq.'V') THEN
  ainfloe(lz) = 0.0
  ENDIF

!.... south face
  IF (varname.eq.'TH'.or.varname.eq.'QV') THEN
  ainflos(lz) = 0.5
  ELSEIF (varname.eq.'U') THEN
  ainflos(lz) = 0.5
  ELSEIF (varname.eq.'V') THEN
  ainflos(lz) = 0.0
  ENDIF

!.... north face
  IF (varname.eq.'TH'.or.varname.eq.'QV') THEN
  ainflon(lz) = 0.0
  ELSEIF (varname.eq.'U') THEN
  ainflon(lz) = 0.0
  ELSEIF (varname.eq.'V') THEN
  ainflon(lz) = 0.0
  ENDIF

  ELSEIF (lz .ge. khomog) THEN

  ainflow(lz) = 0.0
  ainfloe(lz) = 0.0
  ainflos(lz) = 0.0
  ainflon(lz) = 0.0

  ENDIF

!
!....enddo lz = 1, nz
!

  enddo


  
  ENDIF !}
  
  ENDIF !}
  
  ENDIF ! ( varname.eq.'TH'.or.varname.eq.'QV' )


  IF ( nx .gt. 2 .and. bcx .ne. 2 ) THEN
  
  IF ( varname == 'U' ) THEN
    IF ( ixbeg .eq. nxbeg ) THEN
      fs( 1,-ng+1:ny+ng,0:nz) = 0.0
    ENDIF
    IF ( ixend .eq. nxend ) THEN
      fs(nx,-ng+1:ny+ng,0:nz) = 0.0
    ENDIF
  ENDIF
  
   IF ( ixbeg .eq. nxbeg ) THEN ! have western boundary

    jyb = 1
    jye = jtile
    if (jybeg .le. jmn) jyb = jmn
    if (jyend .ge. jmx1) jye = jmx1-jybeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn) kzb = kmn
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1
   
   IF ( imn == 1 .and. is == 0 ) THEN

    do k = kzb,kze 
      do j = jyb,jye

       uw = 0.5*(u(1,j,k) + u(1,j-js,k-ks))
       ue = 0.5*(u(2,j,k) + u(2,j-js,k-ks))
       unorm = uw ! 0.5*(ue + uw)
       fs(1,j,k)   = fs(1,j,k)  - (amin1(unorm,0.0)*(s0(2,j,k)-s0(1,j,k))     &
                                   +  s0(1,j,k)*(ue-uw)       &
                                   +  (ainflow(k)*Max(unorm, 0.0))*(s0(1,j,k) - ndgvarew(j,k,1))) * gxt(1,3)
      ENDDO
     ENDDO
  
   ELSEIF ( varname == 'U' .and. (bcx == 3 .or. bcx == 4) ) THEN

    do k = kzb,kze ; do j = jyb,jye

       uw = u(1,j,k)
       unorm = uw !
       fs(1,j,k)   =   - (amin1(unorm,0.0)*(s0(2,j,k)-s0(1,j,k))     &
                       +  (ainflow(k)*Max(unorm, 0.0))*(s0(1,j,k) - ndgvarew(j,k,1))) * gxt(2,4)
      
      ENDDO
     ENDDO

   
   ENDIF
   ENDIF

   IF ( ixend .eq. nxend ) THEN ! have eastern boundary

    jyb = 1
    jye = jtile
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmx1) jye = jmx1-jybeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

   IF ( imx .eq. nxend-1 .and. is == 0 ) THEN

     do k = kzb,kze
      do j = jyb,jye

       uw = 0.5*(u(nx-1,j,k) + u(nx-1,j-js,k-ks))
       ue = 0.5*(u(nx,  j,k) + u(nx,  j-js,k-ks))
       unorm = ue ! 0.5*(ue + uw)
       fs(nx-1,j,k) = fs(nx-1,j,k) - (amax1(unorm,0.0)*(s0(nx-1,j,k)-s0(nx-2,j,k)) &
                                      +  s0(nx-1,j,k)*(ue-uw)       &
                                      -  (ainfloe(k)*Min(unorm, 0.0)) *(s0(nx-1,j,k) - ndgvarew(j,k,2))) * gxt(nx-1,3)
      ENDDO
     ENDDO

   ELSEIF ( varname == 'U' .and. ( bcx == 3 .or. bcx == 4 ) ) THEN
   
     do k = kzb,kze
      do j = jyb,jye

       ue = u(nx,j,k)
       unorm = ue 

       fs(nx,j,k) =    - (amax1(unorm,0.0)*(s0(nx,j,k)-s0(nx-1,j,k)) &
                       -  (ainfloe(k)*Min(unorm, 0.0)) *(s0(nx,j,k) - ndgvarew(j,k,2))) * gxt(nx-1,4)

      ENDDO
     ENDDO
   
   ENDIF
   
   ENDIF ! ( ixend .eq. nxend )
   
  ENDIF ! ( nx .gt. 2 .and. bcx .ne. 2 )

  IF ( ny .gt. 2 .and. bcy .ne. 2 ) THEN

  IF ( varname == 'V' ) THEN
    IF ( jybeg .eq. nybeg ) THEN
      fs(-ng+1:nx+ng, 1,0:nz) = 0.0
    ENDIF
    IF ( jyend .eq. nyend ) THEN
      fs(-ng+1:nx+ng,ny,0:nz) = 0.0
    ENDIF
  ENDIF

   IF( jybeg .eq. nybeg ) THEN ! have southern boundary

    ixb = 1
    ixe = itile
    if (ixbeg .le. imn)  ixb = imn -ixbeg+1
    if (ixend .ge. imx1) ixe = imx1-ixbeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

   IF( jmn .eq. 1 .and. js == 0 ) THEN

    do k = kzb,kze
      do i = ixb,ixe

       vs = 0.5*(v(i,1,k) + v(i-is,1,k-ks))
       vn = 0.5*(v(i,2,k) + v(i-is,2,k-ks))
       unorm = vs ! 0.5*(vs + vn)
       fs(i,1,k) = fs(i,1,k) - (amin1(unorm,0.0)*(s0(i,2,k)-s0(i,1,k)) &
                                 + s0(i,1,k)*(vn-vs)       &
                                 +  (ainflos(k)*Max(unorm, 0.0))*(s0(i,1,k) - ndgvarns(i,k,1) )) * gyt(1,3)
      ENDDO
     ENDDO

   ELSEIF  ( varname == 'V' .and. (bcy == 3 .or. bcy == 4) ) THEN

    do k = kzb,kze
      do i = ixb,ixe

       vs = v(i,1,k)
       unorm = vs ! 0.5*(vs + vn)

       fs(i,1,k) = fs(i,1,k) - (amin1(unorm,0.0)*(s0(i,2,k)-s0(i,1,k)) &
                                 +  (ainflos(k)*Max(unorm, 0.0))*(s0(i,1,k) - ndgvarns(i,k,1) )) * gyt(1,3)
      ENDDO
     ENDDO
   
   ENDIF

   ENDIF ! jybeg .eq. nybeg

   IF( jyend .eq. nyend ) THEN ! have northern boundary

    ixb = 1
    ixe = itile
    if (ixbeg .le. imn)  ixb = imn -ixbeg+1
    if (ixend .ge. imx1) ixe = imx1-ixbeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1


   IF( jmx .eq. nyend-1 .and. js == 0 ) THEN

    do k = kzb,kze 
      do i = ixb,ixe

       vs = 0.5*(v(i,ny-1,k) + v(i-is,ny-1,k-ks))
       vn = 0.5*(v(i,ny,  k) + v(i-is,ny,  k-ks))
       unorm = vn ! 0.5*(vs + vn)
       fs(i,ny-1,k) = fs(i,ny-1,k) - (amax1(unorm,0.0)*(s0(i,ny-1,k)-s0(i,ny-2,k)) &
                                   + s0(i,ny-1,k)*(vn-vs)       &
                                   -  (ainflon(k)*Min(unorm, 0.0)) *(s0(i,ny-1,k) - ndgvarns(i,k,2) )) * gyt(ny-1,3)
      ENDDO
     ENDDO
   
   ELSEIF ( varname == 'V' .and. ( bcy == 3 .or. bcy == 4 ) ) THEN
   
     do k = kzb,kze
      do i = ixb,ixe

       vn = v(i,ny,k)
       unorm = vn ! 0.5*(vs + vn)
       fs(i,ny,k) = fs(i,ny,k) - (amax1(unorm,0.0)*(s0(i,ny,k)-s0(i,ny-1,k)) &
                                   -  (ainflon(k)*Min(unorm, 0.0)) *(s0(i,ny,k) - ndgvarns(i,k,2) )) * gyt(ny-1,3)


      ENDDO
     ENDDO
   
   ENDIF

   ENDIF
  ENDIF

   CALL cld_cpu('ADVECT-INFLOW')

 IF ( allocated( fx ) ) THEN
  deallocate ( fx )
  deallocate ( fy )
  deallocate ( fz )
 ENDIF
 
 deallocate ( ndgvarns )
 deallocate ( ndgvarew )

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECTRK: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

  RETURN

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM DEFINITIONS (Fortran90 stuff....)

  CONTAINS

! this routine somehow prevents optimization on itanium2 with ifort, so
! a separate subroutine is used instead
#if !defined(__ia64__)
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE ADVECT5(mono)

    implicit none

    logical, optional :: mono
    
    integer         :: i, kz
    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1, km1m,km2m,kp1m
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real, parameter :: eps = 1.0e-25
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

   double precision :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: div2(-ng+1:nx+ng,-ng+1:ny+ng)
   logical :: domono

   integer kt,kb

!-----------------------------------------------------------------------------
!  LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: ENTERING SUBROUTINE, var= ',A,' , my_rank= ',i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-5TH')

    domono = .false.
    IF ( present( mono ) ) domono = mono
      
!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      fx2(:,:) = 0.0d0
      fy2(:,:) = 0.0d0
      fz2(:,:,:) = 0.0d0
      kb = 2
      kt = 1

      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = Max(kzb,kmn-kzbeg+1)
      if (kzend .ge. kmx) kze = Min(kze,kmx1-kzbeg+1)

   if(debug_mpi)  write(0,*) 'Advect5: kzb,kze = ',kzb,kze

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 .and. k > 0 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,ue,uw,dir,km1,i0,im1,im2,ip1,qD,qC,qU,qinmin,qinmax,qi,cr,qoutmin,qoutmax)
       do j = jyb,jye

!DIR$ IVDEP
        DO i = ixb,ixe
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))

         div2(i,j) = gxt(i,3+is)*(ue-uw)
!         div2x(i,j) = gxt(i,3+is)*(ue-uw)

        ENDDO
!DIR$ IVDEP
        do i = ixb,ixe
        im1 = i-1
         IF ( domono ) THEN
         i0   = i
         im2  = i-2
         ip1  = i+1
         ENDIF
         km1 = max(k-1,1)
         vv  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
!         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
!         uw = vv

 !        div2(i,j) = gxt(i,3+is)*(ue-uw)
!         div2x(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

!         IF( nocollapse .or. (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN

!          IF ( dir .ge. 0.0 ) THEN
!          fx2(i,j) = vv * ( 27./60. * s(i,  j,k) + 47./60. * s(i-1,j,k)  &
!                           - 3./60. * s(i+1,j,k) - 13./60. * s(i-2,j,k)  &
!                            + 2./60. * s(i-3,j,k)   )
!                           
!          ELSE
!          fx2(i,j) = vv * ( 47./60. * s(i,  j,k) +  27./60. *s(i-1,j,k)  &
!                           - 13./60. * s(i+1,j,k) - 3./60. * s(i-2,j,k)    &
!                           + 2./60. * s(i+2,j,k)  )
!          ENDIF

          qi =  ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )

         IF ( domono ) THEN
         IF( vv .ge. 0.0 ) THEN
          qD = s0(i0, j,k  )
          qC = s0(im1,j,k)
          qU = s0(im2,j,k)
         ELSE
          qD = s0(im1,j,k)
          qC = s0(i0, j,k)
          qU = s0(ip1,j,k)
         ENDIF

           
         qinmin  = min(qD, qC)
         qinmax  = max(qD, qC)
         qi      = max(qi, qinmin)
         qi      = min(qi, qinmax)
  
         cr      = min(dt*gxt(i,4-is)*abs(vv), 1.0)
         qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
         qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

         qi      = min(max(qi, qoutmin), qoutmax)
         
         ENDIF

          fx2(i,j) = vv * qi

!         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
!          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
!                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
!                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
!                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )
!
!        
!        ELSE
!        
!          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
!         
!          fx2(i,j) = 0.0
         
!         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1) 
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vs,vn,vv,dir,km1,j0,jm1,jm2,jp1,qD,qC,qU,qinmin,qinmax,qi,cr,qoutmin,qoutmax)
      do j = jyb,jye
       jm1 = j-1
          IF ( domono ) THEN
            j0     = j 
            jm1    = j - 1
            jm2    = j - 2
            jp1    = j + 1
          ENDIF

!DIR$ IVDEP
      DO i = ixb,ixe
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)
      ENDDO
!DIR$ IVDEP
      do i = ixb,ixe
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

!         vs  = vv
!         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

!         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)
!!         div2y(i,j) = gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

!         IF( nocollapse .or. (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN

!          IF ( dir .ge. 0.0 ) THEN
!          fy2(i,j) = vv * ( 27./60. * s(i,  j,k) + 47./60. * s(i,j-1,k)  &
!                           - 3./60. * s(i,j+1,k)  - 13./60. * s(i,j-2,k)  &
!                            + 2./60. * s(i,j-3,k)  )
!                          
!          ELSE
!
!          fy2(i,j) = vv * ( 47. * s(i,  j,k) + 27.*s(i,j-1,k)   &
!                           - 13. * s(i,j+1,k) - 3. * s(i,j-2,k)  &
!                           + 2. * s(i,j+2,k)   )/60.
!
!          ENDIF
          
          qi =     ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )

           IF ( domono ) THEN
           IF( vv .ge. 0.0 ) THEN
              qD = s0(i,j0, k)
              qC = s0(i,jm1,k)
              qU = s0(i,jm2,k)
           ELSE
              qD = s0(i,jm1,k)
              qC = s0(i,j0, k)
              qU = s0(i,jp1,k)
           ENDIF
           
             qinmin  = min(qD, qC)
             qinmax  = max(qD, qC)
             qi      = max(qi, qinmin)
             qi      = min(qi, qinmax)
  
             cr      = min(dt*gyt(j,4)*abs(vv), 1.0)
             qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
             qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

             qi      = min(max(qi, qoutmin), qoutmax)
            ENDIF
            
            fy2(i,j)  = vv * qi 

!         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 .or. bcy == 2 ) THEN
!          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
!                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
!                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
!                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )
!
!
!        ELSE
!        
!          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))
!
!
!         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
          km1m    = max(kz-1,1)
          km2m    = max(kz-2,1)
          kp1m    = min(kz+1,k1+ks)
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,wt,wb,dir,qD,qC,qU,qinmin,qinmax,qi,cr,qoutmin,qoutmax)
      do j = jyb,jye
!DIR$ IVDEP
      DO i = ixb,ixe
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)
      ENDDO
!DIR$ IVDEP
      do i = ixb,ixe
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

!        wt  = vv*rrp(kz-1)

!        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

!         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)
!!         div2z(i,j) = gzt(k,3+ks)*(wt-wb)

         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
!          IF ( dir .ge. 0.0 ) THEN
!          fz2(i,j,kt) = vv * ( 27./60. * s(i,j,kz  ) + 47./60.*s(i,j,kz-1)  &
!                           - 3./60. * s(i,j,kz+1) - 13./60. * s(i,j,kz-2)  &
!                           + 2./60. *  s(i,j,kz-3) )
!
!          ELSE
!          
!          fz2(i,j,kt) = vv * ( 47./60. * s(i,j,kz  ) + 27./60.*s(i,j,kz-1)  &
!                           - 13./60. * s(i,j,kz+1) - 3./60. * s(i,j,kz-2)  &
!                           + 2./60. * s(i,j,kz+2)    )
!
!          ENDIF
          
          IF ( Abs(vert_adv_scheme) >= 5 ) THEN
          
          qi =  ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
                           

          IF ( domono ) THEN
!          IF ( .false. ) THEN
           IF( vv .ge. 0.0 ) THEN
              qD = s0(i,j,kz )
              qC = s0(i,j,km1m)
              qU = s0(i,j,km2m)
           ELSE
              qD = s0(i,j,km1m)
              qC = s0(i,j,kz )
              qU = s0(i,j,kp1m)
           ENDIF
           
             qinmin  = min(qD, qC)
             qinmax  = max(qD, qC)
             qi      = max(qi, qinmin)
             qi      = min(qi, qinmax)
  
             cr      = min(dt*gzt(k,4)*abs(vv), 1.0)
             qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
             qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

             qi      = min(max(qi, qoutmin), qoutmax)
            ENDIF
            
            fz2(i,j,kt)  = vv * qi


          ELSE
          
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
          
          ENDIF


         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   do j = jyb,jye
    do i = ixb,ixe
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult
       
       ENDIF
      ENDDO
      ENDDO
      
!      IF ( svar%name == 'W' ) THEN
!        write(luno,*) 'ADVECT5: k,w,fw = ',kzbeg-1+k,w(nx/2,ny/2,k),fs(nx/2,ny/2,k)
!        write(luno,*) 'ADVECT5: fzt,fzb,rrp,rrm,gzt = ', fz2(nx/2,ny/2,kt),fz2(nx/2,ny/2,kb),rrp(k),rrm(k),gzt(k,3+ks)
!      ENDIF
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var= ',A,' , my_rank= ',i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-5TH')

    RETURN
    END SUBROUTINE ADVECT5

#endif

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  3rd-order upwind-biased advection

    SUBROUTINE ADVECT3

    implicit none
    
    integer         :: i, kz
    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1, km1m,km2m,kp1m
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real, parameter :: eps = 1.0e-25
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

   double precision :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng)
   double precision :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: ENTERING SUBROUTINE, var= ',A,' , my_rank= ',i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-5TH')

      
!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      fx2(:,:) = 0.0d0
      fy2(:,:) = 0.0d0
      fz2(:,:,:) = 0.0d0
      kb = 2
      kt = 1

      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = Max(kzb,kmn-kzbeg+1)
      if (kzend .ge. kmx) kze = Min(kze,kmx1-kzbeg+1)

   if(debug_mpi)  write(0,*) 'Advect5: kzb,kze = ',kzb,kze

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 .and. k > 0 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,ue,uw,im1,dir)
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
     km1    = max(k-1,1)
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)
!         div2x(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

         IF( nocollapse .or. (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN


          qi  =    ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

          fx2(i,j) = vv * qi

         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

        
        ELSE
        
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1) 
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vs,vn,vv,jm1,dir)
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)
!         div2y(i,j) = gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

         IF( nocollapse .or. (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN

          qi =       ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )
            
            fy2(i,j)  = vv * qi 

         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 .or. bcy == 2 ) THEN
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )


        ELSE
        
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))


         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
          km1m    = max(kz-1,1)
          km2m    = max(kz-2,1)
          kp1m    = min(kz+1,k1+ks)
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,wt,wb,dir)
      do j = jyb,jye
      do i = ixb,ixe
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)
!         div2z(i,j) = gzt(k,3+ks)*(wt-wb)

         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
           
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )


         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   do j = jyb,jye
    do i = ixb,ixe
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult
       
       ENDIF
      ENDDO
      ENDDO
      
!      IF ( svar%name == 'W' ) THEN
!        write(luno,*) 'ADVECT5: k,w,fw = ',kzbeg-1+k,w(nx/2,ny/2,k),fs(nx/2,ny/2,k)
!        write(luno,*) 'ADVECT5: fzt,fzb,rrp,rrm,gzt = ', fz2(nx/2,ny/2,kt),fz2(nx/2,ny/2,kb),rrp(k),rrm(k),gzt(k,3+ks)
!      ENDIF
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var= ',A,' , my_rank= ',i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-5TH')

    RETURN
    END SUBROUTINE ADVECT3



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECTN:  5th-order upwind-biased advection with No stencil collapse

    SUBROUTINE ADVECT5N()
    
    integer         :: i, im1, jm1, km1
    real            :: dirx,diry,dirz, uu, vv, ww
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    CALL cld_cpu('ADVECT-5N')


!#ifdef MPI
!   if(debug_mpi)  write(0,"('ADVECT5N: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
!   write(0,*) "ADVECT5N: not updated for MPI yet"
!#endif


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES


    IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,jm1,dirx,diry,dirz)
      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
!DIR$ IVDEP
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,jm1,dirx,diry,dirz)
      DO k = kmn,kmx1
       km1 = max(k-1,1)
       DO j = jmn,jmxb 
        jm1 = max(j-1,1)
!DIR$ IVDEP
        DO i = imn,imxb  
         im1 = max(i-1,1)
#endif
        
!     x-direction:
        
          uu  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dirx = sign(1.0,uu)
         
          fx(i,j,k) = uu * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dirx )
!     y-direction:
         
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         diry = sign(1.0,vv)

          fy(i,j,k) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*diry )

!     z-direction:

         ww  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dirz = sign(1.0,ww)

          fz(i,j,k) = ww * ( f50 * (s(i,j,k  ) + s(i,j,k-1))  &
                           - f51 * (s(i,j,k+1) + s(i,j,k-2))  &
                           + f52 * (s(i,j,k+2) + s(i,j,k-3))  &
                           - f52 * (s(i,j,k+2) - s(i,j,k-3)   &
                           - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
                           + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dirz )

     
        ENDDO
       ENDDO
      ENDDO

    ELSEIF ( nx .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,dirx,dirz)
      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,dirx,dirz)
      DO k = kmn,kmx1
       km1 = max(k-1,1)
       DO j = jmn,jmxb 
        DO i = imn,imxb  
         im1 = max(i-1,1)
#endif
        
!     x-direction:
        
          uu  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dirx = sign(1.0,uu)
         
          fx(i,j,k) = uu * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dirx )


         ww  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dirz = sign(1.0,ww)

          fz(i,j,k) = ww * ( f50 * (s(i,j,k  ) + s(i,j,k-1))  &
                           - f51 * (s(i,j,k+1) + s(i,j,k-2))  &
                           + f52 * (s(i,j,k+2) + s(i,j,k-3))  &
                           - f52 * (s(i,j,k+2) - s(i,j,k-3)   &
                           - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
                           + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dirz )

     
        ENDDO
       ENDDO
      ENDDO

    ELSEIF ( ny .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = imn-ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,jm1,diry,dirz)
      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,jm1,diry,dirz)
      DO k = kmn,kmx1
       km1 = max(k-1,1)
       DO j = jmn,jmxb 
        jm1 = max(j-1,1)
        DO i = imn,imxb  
#endif
        
!     y-direction:
         
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         diry = sign(1.0,vv)

          fy(i,j,k) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*diry )

!     z-direction:

         ww  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dirz = sign(1.0,ww)

          fz(i,j,k) = ww * ( f50 * (s(i,j,k  ) + s(i,j,k-1))  &
                           - f51 * (s(i,j,k+1) + s(i,j,k-2))  &
                           + f52 * (s(i,j,k+2) + s(i,j,k-3))  &
                           - f52 * (s(i,j,k+2) - s(i,j,k-3)   &
                           - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
                           + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dirz )

     
        ENDDO
       ENDDO
      ENDDO

    ENDIF
    

    CALL cld_cpu('ADVECT-5N')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5N: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    
    RETURN
    END SUBROUTINE ADVECT5N



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  6th-order advection

    SUBROUTINE ADVECT6()
    
    implicit none
    
    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
    real            :: uw, vs, wb

    CALL cld_cpu('ADVECT-6TH')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT6: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

    IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,vs,wb)
     DO k = kmn,kmx1
      DO j = jmn,jmxb
       DO i = imn,imxb
#endif
          uw = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
          vs = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
          wb = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))

         fx(i,j,k)   = uw * (f60 * (s(i,  j,k) + s(i-1,j,k))  &
                            -f61 * (s(i+1,j,k) + s(i-2,j,k))  &
                            +f62 * (s(i+2,j,k) + s(i-3,j,k)))

         fy(i,j,k)   = vs * (f60 * (s(i,j,  k) + s(i,j-1,k))  &
                            -f61 * (s(i,j+1,k) + s(i,j-2,k))  &
                            +f62 * (s(i,j+2,k) + s(i,j-3,k)))

         fz(i,j,k)   = wb * (f60 * (s(i,j,k  ) + s(i,j,k-1))  &
                            -f61 * (s(i,j,k+1) + s(i,j,k-2))  &
                            +f62 * (s(i,j,k+2) + s(i,j,k-3)))

         ENDDO
        ENDDO
       ENDDO

    ELSEIF ( nx .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,wb)
     DO k = kmn,kmx1
      DO j = jmn,jmxb
       DO i = imn,imxb
#endif
          uw = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
          wb = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))

         fx(i,j,k)   = uw * (f60 * (s(i,  j,k) + s(i-1,j,k))  &
                            -f61 * (s(i+1,j,k) + s(i-2,j,k))  &
                            +f62 * (s(i+2,j,k) + s(i-3,j,k)))

         fy(i,j,k)   = 0

         fz(i,j,k)   = wb * (f60 * (s(i,j,k  ) + s(i,j,k-1))  &
                            -f61 * (s(i,j,k+1) + s(i,j,k-2))  &
                            +f62 * (s(i,j,k+2) + s(i,j,k-3)))

         ENDDO
        ENDDO
       ENDDO

    ELSEIF ( ny .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,vs,wb)
     DO k = kmn,kmx1
      DO j = jmn,jmxb
       DO i = imn,imxb
#endif
          vs = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
          wb = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))

         fx(i,j,k)   = 0

         fy(i,j,k)   = vs * (f60 * (s(i,j,  k) + s(i,j-1,k))  &
                            -f61 * (s(i,j+1,k) + s(i,j-2,k))  &
                            +f62 * (s(i,j+2,k) + s(i,j-3,k)))

         fz(i,j,k)   = wb * (f60 * (s(i,j,k  ) + s(i,j,k-1))  &
                            -f61 * (s(i,j,k+1) + s(i,j,k-2))  &
                            +f62 * (s(i,j,k+2) + s(i,j,k-3)))

         ENDDO
        ENDDO
       ENDDO

    ENDIF

    CALL cld_cpu('ADVECT-6TH')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT6: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    RETURN
    END SUBROUTINE ADVECT6

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_MONOF1:  Rather complicated general advection algorithm.  See above

    SUBROUTINE ADVECT_MONOF1(fdamp1)
    
    implicit none
    
    integer :: mono = 2
    real    :: fdamp1

! Local vars

    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real            :: dtbydx, dtbydy, dtbydz, del
    real, parameter :: f40 =  7./12., f41 = 1./12
    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
    real, parameter :: eps = 1.0e-28
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1
    real :: fd1,fd2
!    real :: gxs(nx), gxw(nx)
    real :: dt2
!    real :: s0x(-ng+1:nx+ng)
!    real :: s0y(-ng+1:ny+ng)
!    real,allocatable :: s0t(:,:,:)
    
    real :: qr, qdel, qxs, qxp, t0s, fvs

!-----------------------------------------------------------------------------

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONOF1: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT_MONOF1')

!    allocate( s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )

!    s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) = s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

    dt2 = dt ! 5.0

!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

    IF ( nx .gt. 2 ) THEN
        IF (  bcx .ne. 2 ) THEN
#ifdef MPI
          imn1 = 1
#else
          imn1 = Max(2,imn)
!          imx1 = Min(nx-1,imx)
#endif

#ifndef MPI
!no need to fill ghost zones when using MPI
         DO n = 1,ng
          s0(1-n   ,1:ny-1,1:nz-1) = s0(1   ,1:ny-1,1:nz-1)
          s0(nx+n-1,1:ny-1,1:nz-1) = s0(nx-1,1:ny-1,1:nz-1)
         ENDDO
#endif
        ELSE ! periodic
          imn1 = imn
#ifndef MPI
         DO n = 1,ng
          s0(1-n   ,1:ny-1,1:nz-1) = s0(nx-n,1:ny-1,1:nz-1)
          s0(nx+n-1,1:ny-1,1:nz-1) = s0(n   ,1:ny-1,1:nz-1)
         ENDDO
#endif
        ENDIF

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn)  ixb = Max( ixb, imn-ixbeg+1 )
    if (ixend .ge. imxb) ixe = imxb-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2, &
!$OMP         qdel, fvs, qr, t0s)
	DO k = kmn,kmx1 
	 DO j = jmn,jmx 
	 
!       s0x(1:nx-1) = s0(1:nx-1,j,k)
	  
	  DO i = imn,imxb
#endif
            dtbydx = dt2*gxt(i,3+is)

              i0     = i
              im1    = i-1
              im2    = i-2
              ip1    = i+1

!            IF ( bcx .ne. 2 ) THEN
              del    = s(i,j,k) - s(im1,j,k)
!            ELSE
!              del    = s(i,j,k) - s(i-1,j,k)
!            ENDIF

		qi     = 0.0
		qid    = 0.0
		damp   = 0.0
		vv     = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
  
#ifdef MPI
        IF( (nocollapse .and. ixbeg-1+i .le. i1) .or. (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
        IF( (nocollapse .and. i .le. i1) .or. (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2 ) THEN
#endif
                 qi   =    f60 * (s(i,  j,k) + s(i-1,j,k))  &
                         - f61 * (s(i+1,j,k) + s(i-2,j,k))  &
                         + f62 * (s(i+2,j,k) + s(i-3,j,k))  
		 
		 qid  = (-s(i-3,j,k)+5.*s(i-2,j,k)-10.*s(i-1,j,k)+10.*s(i,j,k)-5.*s(i+1,j,k)+s(i+2,j,k))
		 qid  = qid*Max(0.0, sign(1.0,qid*del) )
		 damp = fdamp1 / (64.*dtbydx) 

#ifdef MPI
        ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
		ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
                 qi   =   f40 * (s(i,  j,k) + s(i-1,j,k))  &
                        - f41 * (s(i+1,j,k) + s(i-2,j,k))
                  
		 qid  = -(s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)) )
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
             damp =  fdamp1 / (16.*dtbydx)

#ifdef MPI
        ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
		ELSEIF ( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
		 qi   = 0.5 * (s(i,j,k) + s(i-1,j,k))
		  
		 qid  = (s(i,j,k)-s(i-1,j,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp =  fdamp1 / (4.*dtbydx) 


          ELSE
          
           qi = 0.0
           fx(i,j,k) = 0.0
           CYCLE
          
		ENDIF

        
	     IF( vv .ge. 0.0 ) THEN
		  qD = s0(i0, j,k  )
		  qC = s0(im1,j,k)
		  qU = s0(im2,j,k)
	     ELSE
		  qD = s0(im1,j,k)
		  qC = s0(i0, j,k)
		  qU = s0(ip1,j,k)
	     ENDIF

           IF ( mono == 1 ) THEN 
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gxt(i,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)

           ELSEIF ( mono == 2 ) THEN
           
             IF ( Abs(vv) > 1.e-6 ) THEN

               qdel = qd - qu
             
               IF ( abs(qdel) <= abs(qu-(2.0)*qc+qd) ) THEN
                 qi = qc
               ELSE
                 fvs = vv*dt*gxt(i,4)
                 qr   = qu + (qc - qu) / (abs(fvs))  ! phi_ref
                 t0s = qi
                 IF ( qdel > 0 ) THEN
                   qi = Min( Min(qr,qd) ,Max(qc,qi) )
                 ELSE
                   qi = Max( Max(qr,qd) ,Min(qc,qi) )
                 ENDIF
!                 qxs  = max(sign((1.0), qdel), (0.0))
!                 qxp  = max(sign((1.0), abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))

!                qi =    -dt*gxt(i,4) *   &
!                (   &
!                ((1.0)-qxp) * qc   &
!              +        qxp   &
!              * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
!              + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
           
             ENDIF
             ENDIF
           ENDIF

		fx(i,j,k)  = vv * qi - damp * qid

	   ENDDO
	  ENDDO
	 ENDDO

     ENDIF 

!-----------------------------------------------------------------------------
! COMPUTE Y-INTERFACE FLUXES

    IF ( ny .gt. 2 ) THEN

     IF (  bcy .ne. 2 ) THEN
#ifdef MPI
       jmn1 = 1
#else
       jmn1 = Max(2,jmn)
#endif
#ifndef MPI
       DO n = 1,ng
        s0(1:nx-1,1-n   ,1:nz-1) = s0(1:nx-1,1   ,1:nz-1)
        s0(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,ny-1,1:nz-1)
       ENDDO
#endif
     ELSE ! periodic
#ifndef MPI
       jmn1 = jmn
       DO n = 1,ng
        s0(1:nx-1,1-n   ,1:nz-1) = s0(1:nx-1,ny-n,1:nz-1)
        s0(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,n   ,1:nz-1)
       ENDDO
#endif
     ENDIF

!   write(0,*) 'mono, y-dir , var = ', varname

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmxb) jye = jmxb-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye

            j0     = j 
            jm1    = j-1 
            jm2    = j-2 
            jp1    = j+1 

      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,k,j0,jm1,jm2,jp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax, &
!$OMP         qdel, fvs, qr, t0s)
	DO k = kmn,kmx1 
	 DO j = jmn,jmxb

            j0     = j 
            jm1    = j - 1 
            jm2    = j - 2 
            jp1    = j + 1 

	  DO i = imn,imx 
#endif
		dtbydx = dt2*gyt(j,3+js)
!            IF ( bcy .ne. 2 ) THEN
		del    = s(i,j,k) - s(i,jm1,k)   
!            ELSE 
!                del    = s(i,j,k) - s(i,j-1,k)
!            ENDIF
		qid    = 0.0
		qi     = 0.0
		damp   = 0.0
		vv     = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
  
#ifdef MPI
        IF( (nocollapse .and. jybeg-1+j .le. j1 ) .or. (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
        IF( (nocollapse .and. j .le. j1 ) .or. ( j .ge. 4 .and. j .le. j3 ) .or. bcy .eq. 2 ) THEN
#endif		
         qi   = f60 * (s(i,j,  k) + s(i,j-1,k))  &
              - f61 * (s(i,j+1,k) + s(i,j-2,k))  &
              + f62 * (s(i,j+2,k) + s(i,j-3,k))  
		 
		 qid  = (-s(i,j-3,k)+5.*s(i,j-2,k)-10.*s(i,j-1,k)+10.*s(i,j,k)-5.*s(i,j+1,k)+s(i,j+2,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp = fdamp1 / (64.*dtbydx) 

#ifdef MPI
        ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
		ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
         qi   = f40 * (s(i,j,  k) + s(i,j-1,k))  &
              - f41 * (s(i,j+1,k) + s(i,j-2,k))
                      
		 qid  = -(s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)) )
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp = fdamp1 / (16.*dtbydx) 

#ifdef MPI
        ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
		ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
		 qi   = 0.5 * (s(i,j,k) + s(i,j-1,k))
		  
		 qid  = (s(i,j,k)-s(i,j-1,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
 		 damp =  fdamp1 / (4.*dtbydx) 

		ENDIF
        
           IF( vv .ge. 0.0 ) THEN
              qD = s0(i,j0, k)
              qC = s0(i,jm1,k)
              qU = s0(i,jm2,k)
           ELSE
              qD = s0(i,jm1,k)
              qC = s0(i,j0, k)
              qU = s0(i,jp1,k)
           ENDIF
     
           IF ( mono == 1 ) THEN
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gyt(j,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)

           ELSEIF ( mono == 2 ) THEN
           
             IF ( Abs(vv) > 1.e-6 ) THEN

               qdel = qd - qu
             
               IF ( abs(qdel) <= abs(qu-(2.0)*qc+qd) ) THEN
                 qi = qc
               ELSE
                 fvs = vv*dt*gyt(j,4)
                 qr   = qu + (qc - qu) / (abs(fvs))  ! phi_ref
                 t0s = qi
                 IF ( qdel > 0 ) THEN
                   qi = Min( Min(qr,qd) ,Max(qc,qi) )
                 ELSE
                   qi = Max( Max(qr,qd) ,Min(qc,qi) )
                 ENDIF
!                 qxs  = max(sign((1.0), qdel), (0.0))
!                 qxp  = max(sign((1.0), abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))

!                qi =    -dt*gxt(i,4) *   &
!                (   &
!                ((1.0)-qxp) * qc   &
!              +        qxp   &
!              * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
!              + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
           
             ENDIF
             ENDIF

           ENDIF

		fy(i,j,k)  = vv * qi - damp * qid

	   ENDDO
	  ENDDO
	 ENDDO
     
      ENDIF

!-----------------------------------------------------------------------------
! COMPUTE Z-INTERFACE FLUXES

    IF ( nz .gt. 2 ) THEN

!     kmn1 = Max(2,kmn)

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze

     km1    = k-1
     km2    = k-2
     kp1    = k+1
     
     IF ( myprock == 1 ) THEN
       km1    = max(k-1,1)
       km2    = max(k-2,1)
     ENDIF
     
     IF ( myprock == nprock ) THEN
        kp1    = min(k+1,k1+ks)
     ENDIF

     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax, &
!$OMP         qdel, fvs, qr, t0s)
     DO k = kmn,kmx1

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,k1+ks)

	 DO j = jmn,jmx 
	  DO i = imn,imx 
#endif
         qi  = 0.0
         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)


         IF( k+kzbeg-1 .ge. 4 .and. k+kzbeg-1 .le. k3 ) THEN

          IF ( Abs(vert_adv_scheme) >= 5 ) THEN

           IF ( mono /= 12 ) THEN
! 5th order:
           qi = ( f60 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f61 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f62 * (s(i,j,k+2) + s(i,j,k-3))  &
               - f62 * (s(i,j,k+2) - s(i,j,k-3)   &
               - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
               + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dir )
            ELSE
! 6th order:
            qi = (f60 * (s(i,j,k  ) + s(i,j,k-1))  &
                 -f61 * (s(i,j,k+1) + s(i,j,k-2))  &
                 +f62 * (s(i,j,k+2) + s(i,j,k-3)))
           ENDIF
 
          ELSE

! 3rd order
          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f41 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f41 * (s(i,j,k+1) - s(i,j,k-2)   &
               - 3.0 * (s(i,j,k  ) - s(i,j,k-1)))*dir )
! 4th order
!          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
!               - f41 * (s(i,j,k+1) + s(i,j,k-2)))
              

          ENDIF
         ELSEIF( k+kzbeg-1 .eq. 3 .or. k+kzbeg-1 .eq. k2 ) THEN

! 3rd order
          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f41 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f41 * (s(i,j,k+1) - s(i,j,k-2)   &
               - 3.0 * (s(i,j,k  ) - s(i,j,k-1)))*dir )

! 4th order
!          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
!               - f41 * (s(i,j,k+1) + s(i,j,k-2)))

         ELSEIF( k+kzbeg-1 .le. 2 .or. k+kzbeg-1 .ge. k1 ) THEN

          qi = 0.5 * (s(i,j,k) + s(i,j,km1))
    
!         ELSE
         
!           fz(i,j,k) = 0.0
!           CYCLE
         
         ENDIF
        
	     IF( vv .ge. 0.0 ) THEN
		  qD = s0(i,j,k  )
		  qC = s0(i,j,km1)
		  qU = s0(i,j,km2)
	     ELSE
		  qD = s0(i,j,km1)
		  qC = s0(i,j,k  )
		  qU = s0(i,j,kp1)
	     ENDIF


          IF ( mono == 1 ) THEN
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gzt(k,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)

           ELSEIF ( mono == 2 ) THEN
           
             IF ( Abs(vv) > 1.e-6 ) THEN

               qdel = qd - qu
             
               IF ( abs(qdel) <= abs(qu-(2.0)*qc+qd) ) THEN
                 qi = qc
               ELSE
                 fvs = vv*dt*gzt(k,4)
                 qr   = qu + (qc - qu) / (abs(fvs))  ! phi_ref
                 t0s = qi
                 IF ( qdel > 0 ) THEN
                   qi = Min( Min(qr,qd) ,Max(qc,qi) )
                 ELSE
                   qi = Max( Max(qr,qd) ,Min(qc,qi) )
                 ENDIF
!                 qxs  = max(sign((1.0), qdel), (0.0))
!                 qxp  = max(sign((1.0), abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))

!                qi =    -dt*gxt(i,4) *   &
!                (   &
!                ((1.0)-qxp) * qc   &
!              +        qxp   &
!              * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
!              + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
           
               ENDIF
             ENDIF
           ENDIF
		fz(i,j,k)  = vv * qi

	   ENDDO
	  ENDDO
	 ENDDO

     ENDIF
     
!    deallocate( s0t )

    CALL cld_cpu('ADVECT_MONOF1')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONOF1: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    RETURN
    END SUBROUTINE ADVECT_MONOF1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_MONOF:  Rather complicated general advection algorithm.  See above

    SUBROUTINE ADVECT_MONOF(fdamp1,fdamp2,mono)
    
    implicit none
    
    logical :: mono
    real    :: fdamp1, fdamp2

! Local vars

    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real            :: dtbydx, dtbydy, dtbydz, del
    real, parameter :: f40 =  7./12., f41 = 1./12
    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
    real, parameter :: eps = 1.0e-25
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1
    real :: fd1,fd2
!    real :: gxs(nx), gxw(nx)
    real :: dt2
!    real :: s0x(-ng+1:nx+ng)
!    real :: s0y(-ng+1:ny+ng)
    real,allocatable :: s0t(:,:,:)

!-----------------------------------------------------------------------------


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONOF: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT_MONOF')

    allocate( s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )

    s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) = s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

    dt2 = dt ! 5.0

!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

    IF ( nx .gt. 2 ) THEN
        IF (  bcx .ne. 2 ) THEN
!          imn1 = Max(2,imn)
!          imx1 = Min(nx-1,imx)
#ifndef MPI
         DO n = 1,ng
          s0t(1-n   ,1:ny-1,1:nz-1) = s0(1   ,1:ny-1,1:nz-1)
          s0t(nx+n-1,1:ny-1,1:nz-1) = s0(nx-1,1:ny-1,1:nz-1)
         ENDDO
#endif
        ELSE ! periodic
          imn1 = imn
         DO n = 1,ng
          s0t(1-n   ,1:ny-1,1:nz-1) = s0(nx-n,1:ny-1,1:nz-1)
          s0t(nx+n-1,1:ny-1,1:nz-1) = s0(n   ,1:ny-1,1:nz-1)
         ENDDO
        ENDIF

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
    if (ixend .ge. imxb) ixe = imxb-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2)
	DO k = kmn,kmx1 
	 DO j = jmn,jmx 
	 
!       s0x(1:nx-1) = s0(1:nx-1,j,k)
	  
	  DO i = imn,imxb
#endif
            dtbydx = dt2*gxt(i,3+is)
              i0     = i
              im1    = i-1
              im2    = i-2
              ip1    = i+1

            IF ( bcx .ne. 2 ) THEN
!              im1    = max(i-1,1)
!              im2    = max(i-2,1)
!              ip1    = min(i+1,i1)
              del    = s(i,j,k) - s(im1,j,k)
            ELSE
              del    = s(i,j,k) - s(i-1,j,k)
            ENDIF

		qi     = 0.0
		qid    = 0.0
		damp   = 0.0
		vv     = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))

#ifdef MPI
        IF( (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else  
		IF( (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2 ) THEN
#endif
                 qi   =    f60 * (s(i,  j,k) + s(i-1,j,k))  &
                         - f61 * (s(i+1,j,k) + s(i-2,j,k))  &
                         + f62 * (s(i+2,j,k) + s(i-3,j,k))  
		 
		 qid  = (-s(i-3,j,k)+5.*s(i-2,j,k)-10.*s(i-1,j,k)+10.*s(i,j,k)-5.*s(i+1,j,k)+s(i+2,j,k))
		 qid  = qid*Max(0.0, sign(1.0,qid*del) )
		 damp = Max( fdamp2*abs(vv) / 60., fdamp1 / (64.*dtbydx) )

#ifdef MPI
        ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
		ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
                 qi   =   f40 * (s(i,  j,k) + s(i-1,j,k))  &
                        - f41 * (s(i+1,j,k) + s(i-2,j,k))
                  
		 qid  = -(s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)) )
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
                 damp = Max( fdamp2*abs(vv) / 12., fdamp1 / (16.*dtbydx) )
#ifdef MPI
        ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
		ELSEIF ( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
		 qi   = 0.5 * (s(i,j,k) + s(i-1,j,k))
		  
		 qid  = (s(i,j,k)-s(i-1,j,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp = Max( fdamp2*abs(vv) / 2., fdamp1 / (4.*dtbydx) )

		ENDIF

           IF( mono ) THEN
        
	     IF( vv .ge. 0.0 ) THEN
		  qD = s0t(i0, j,k  )
		  qC = s0t(im1,j,k)
		  qU = s0t(im2,j,k)
	     ELSE
		  qD = s0t(im1,j,k)
		  qC = s0t(i0, j,k)
		  qU = s0t(ip1,j,k)
	     ENDIF

	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gxt(i,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		ENDIF
		
		fx(i,j,k)  = vv * qi - damp * qid

	   ENDDO
	  ENDDO
	 ENDDO
	 
    ENDIF

!-----------------------------------------------------------------------------
! COMPUTE Y-INTERFACE FLUXES

    IF ( ny .gt. 2 ) THEN
        

     IF (  bcy .ne. 2 ) THEN
!       jmn1 = Max(2,jmn)
#ifndef MPI
       DO n = 1,ng
        s0t(1:nx-1,1-n,   1:nz-1) = s0(1:nx-1,1   ,1:nz-1)
        s0t(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,ny-1,1:nz-1)
       ENDDO
#endif
     ELSE ! periodic
       jmn1 = jmn
#ifndef MPI
       DO n = 1,ng
        s0t(1:nx-1,1-n  , 1:nz-1) = s0(1:nx-1,ny-n,1:nz-1)
        s0t(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,n   ,1:nz-1)
       ENDDO
#endif
     ENDIF

!   write(0,*) 'mono, y-dir , var = ', varname

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmxb) jye = jmxb-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye

            j0     = j 
            jm1    = j-1 
            jm2    = j-2 
            jp1    = j+1 

      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,k,j0,jm1,jm2,jp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
	DO k = kmn,kmx1 
	 DO j = jmn,jmxb

!       s0y(1:ny-1) = s0(i,1:ny-1,k)

            j0     = j 
            jm1    = j - 1 
            jm2    = j - 2 
            jp1    = j + 1 
!            IF ( bcy .ne. 2 ) THEN
!              jm1    = max(j-1,1)
!              jm2    = max(j-2,1)
!              jp1    = min(j+1,j1)
!            ENDIF

	  DO i = imn,imx 
#endif
		dtbydx = dt2*gyt(j,3+js)
            IF ( bcy .ne. 2 ) THEN
		del    = s(i,j,k) - s(i,jm1,k)   
            ELSE 
                del    = s(i,j,k) - s(i,j-1,k)
            ENDIF
		qid    = 0.0
		qi     = 0.0
		damp   = 0.0
		vv     = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

#ifdef MPI
        IF( (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else  
		IF( ( j .ge. 4 .and. j .le. j3 ) .or. bcy .eq. 2 ) THEN
#endif
         qi   = f60 * (s(i,j,  k) + s(i,j-1,k))  &
              - f61 * (s(i,j+1,k) + s(i,j-2,k))  &
              + f62 * (s(i,j+2,k) + s(i,j-3,k))  
		 
		 qid  = (-s(i,j-3,k)+5.*s(i,j-2,k)-10.*s(i,j-1,k)+10.*s(i,j,k)-5.*s(i,j+1,k)+s(i,j+2,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp = Max( fdamp2*abs(vv) / 60., fdamp1 / (64.*dtbydx) )

#ifdef MPI
        ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
		ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
         qi   = f40 * (s(i,j,  k) + s(i,j-1,k))  &
              - f41 * (s(i,j+1,k) + s(i,j-2,k))
                      
		 qid  = -(s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)) )
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
		 damp = Max( fdamp2*abs(vv) / 12., fdamp1 / (16.*dtbydx) )

#ifdef MPI
        ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
		ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
		 qi   = 0.5 * (s(i,j,k) + s(i,j-1,k))
		  
		 qid  = (s(i,j,k)-s(i,j-1,k))
		 qid  = qid*amax1(0.0, sign(1.0,qid*del) )
 		 damp = Max( fdamp2*abs(vv) / 2., fdamp1 / (4.*dtbydx) )

		ENDIF

        IF( mono ) THEN
        
           IF( vv .ge. 0.0 ) THEN
              qD = s0t(i,j0, k)
              qC = s0t(i,jm1,k)
              qU = s0t(i,jm2,k)
           ELSE
              qD = s0t(i,jm1,k)
              qC = s0t(i,j0, k)
              qU = s0t(i,jp1,k)
           ENDIF
	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gyt(j,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		ENDIF
		
		fy(i,j,k)  = vv * qi - damp * qid

	   ENDDO
	  ENDDO
	 ENDDO
	 
!   write(0,*) 'mono2, y-dir done'
    ENDIF

!-----------------------------------------------------------------------------
! COMPUTE Z-INTERFACE FLUXES

    IF ( nz .gt. 2 ) THEN

!     kmn1 = Max(2,kmn)
#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx1-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze

     km1 = max(k-1,1)
     km2 = max(k-2,1)
     kp1 = min(k+1,k1+ks)

     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
     DO k = kmn,kmx1

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,nz-1+ks)

	 DO j = jmn,jmx 
	  DO i = imn,imx 
#endif
		 qi  = 0.0
         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF( kzbeg-1+k .ge. 4 .and. kzbeg-1+k .le. k3 ) THEN
#else
         IF( k .ge. 4 .and. k .le. k3 ) THEN
#endif
          qi = ( f60 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f61 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f62 * (s(i,j,k+2) + s(i,j,k-3))  &
               - f62 * (s(i,j,k+2) - s(i,j,k-3)   &
               - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
               + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 3 .or. kzbeg-1+k .eq. k2 ) THEN
#else
         ELSEIF( k .eq. 3 .or. k .eq. k2 ) THEN
#endif
          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f41 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f41 * (s(i,j,k+1) - s(i,j,k-2)   &
               - 3.0 * (s(i,j,k  ) - s(i,j,k-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 2 .or. kzbeg-1+k .eq. k1 ) THEN
#else
         ELSEIF( k .eq. 2 .or. k .eq. k1 ) THEN
#endif
          qi = 0.5 * (s(i,j,k) + s(i,j,km1))
    
!         ELSE
!         
!           fz(i,j,k) = 0.0
!           CYCLE
         
         ENDIF

        IF( mono ) THEN
        
	     IF( vv .ge. 0.0 ) THEN
		  qD = s0(i,j,k  )
		  qC = s0(i,j,km1)
		  qU = s0(i,j,km2)
	     ELSE
		  qD = s0(i,j,km1)
		  qC = s0(i,j,k  )
		  qU = s0(i,j,kp1)
	     ENDIF
	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gzt(k,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		ENDIF
		
		fz(i,j,k)  = vv * qi

	   ENDDO
	  ENDDO
	 ENDDO
	 
    ENDIF

    deallocate( s0t )

    CALL cld_cpu('ADVECT_MONOF')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONOF: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    RETURN
    END SUBROUTINE ADVECT_MONOF

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_MONO:  Rather complicated general advection algorithm.  See above

    SUBROUTINE ADVECT_MONO()
    
    implicit none

! Local vars

    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real            :: dtbydx, dtbydy, dtbydz, del
    real, parameter :: f40 =  7./12., f41 = 1./12
    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
    real, parameter :: eps = 1.0e-25
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1
    real :: fd1,fd2
    real :: dt2
    real,allocatable :: s0t(:,:,:)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONO: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

!-----------------------------------------------------------------------------

    CALL cld_cpu('ADVECT_MONO')

    allocate( s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )

    s0t(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) = s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

    dt2 = dt ! 5.0

!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

    IF ( nx .gt. 2 ) THEN
        IF (  bcx .ne. 2 ) THEN
!#ifdef MPI
!          imn1 = Max(2,imn)
!          imx1 = Min(nxend-1,imx)
!#else
!          imn1 = Max(2,imn)
!          imx1 = Min(nx-1,imx)
!#endif

#ifndef MPI
         DO n = 1,ng  !!mpidebug: no need to fill ghost zones if filled at start of ADVECT_MONO
          s0t(1-n   ,1:ny-1,1:nz-1) = s0t(1   ,1:ny-1,1:nz-1)
          s0t(nx+n-1,1:ny-1,1:nz-1) = s0t(nx-1,1:ny-1,1:nz-1)
         ENDDO
#endif
        ELSE ! periodic

          imn1 = imn
         DO n = 1,ng
#ifdef MPI
          jyb = 1
          jye = jtile+1
          if (jyend .eq. nyend) jye = jyend-jybeg

          kzb = 1
          kze = ktile+1
          if (kzbeg .eq. nzbeg) kzb = 1
          if (kzend .eq. nzend) kze = kzend-kzbeg

          do k = kzb,kze ; do j = jyb,jye
           s0t(1-n   ,j,k) = s0t(nxend-n,j,k)
           s0t(nxend+n-1,j,k) = s0t(n   ,j,k)
          enddo;enddo
#else
          s0t(1-n   ,1:ny-1,1:nz-1) = s0t(nx-n,1:ny-1,1:nz-1)
          s0t(nx+n-1,1:ny-1,1:nz-1) = s0t(n   ,1:ny-1,1:nz-1)
#endif
         ENDDO
        ENDIF !! (  bcx .ne. 2 )

#ifdef MPI
       ixb = 1
       ixe = itile+1
       if (ixbeg .le. imn ) ixb = Max(ixb,imn -ixbeg+1)
       if (ixend .ge. imxb) ixe = imxb-ixbeg+1

       jyb = 1
       jye = jtile+1
       if (jybeg .le. jmn) jyb = jmn-jybeg+1
       if (jyend .ge. jmx) jye = jmx-jybeg+1

       kzb = 1
       kze = ktile+1
       if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
       if (kzend .ge. kmx) kze = kmx1-kzbeg+1

       do k = kzb,kze
        do j = jyb,jye
         do i = ixb,ixe
#else

!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2)
       DO k = kmn,kmx1 
        DO j = jmn,jmx 
         DO i = imn,imxb
#endif

         i0   = i
         im1  = i-1
         im2  = i-2
         ip1  = i+1

	   vv   = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
  
         qi   =    f60 * (s(i,  j,k) + s(i-1,j,k))  &
                 - f61 * (s(i+1,j,k) + s(i-2,j,k))  &
                 + f62 * (s(i+2,j,k) + s(i-3,j,k))  
		 
         IF( vv .ge. 0.0 ) THEN
	    qD = s0t(i0, j,k  )
	    qC = s0t(im1,j,k)
	    qU = s0t(im2,j,k)
	   ELSE
	    qD = s0t(im1,j,k)
	    qC = s0t(i0, j,k)
	    qU = s0t(ip1,j,k)
	   ENDIF

	     
	   qinmin  = min(qD, qC)
	   qinmax  = max(qD, qC)
	   qi      = max(qi, qinmin)
	   qi      = min(qi, qinmax)
  
	   cr      = min(dt2*gxt(i,4-is)*abs(vv), 1.0)
	   qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
	   qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

	   qi      = min(max(qi, qoutmin), qoutmax)
		 
	   fx(i,j,k)  = vv * qi 

         ENDDO
        ENDDO
       ENDDO
     
    ENDIF

!-----------------------------------------------------------------------------
! COMPUTE Y-INTERFACE FLUXES

    IF ( ny .gt. 2 ) THEN
        

     IF (  bcy .ne. 2 ) THEN
!       jmn1 = Max(2,jmn)
#ifndef MPI
       DO n = 1,ng  !!mpidebug: no need to fill ghost zones if filled above
        s0t(1:nx-1,1-n, 1:nz) = s0t(1:nx-1,1,1:nz)
        s0t(1:nx-1,ny+n-1,1:nz) = s0t(1:nx-1,ny-1,1:nz)
       ENDDO
#endif
     ELSE ! periodic
       jmn1 = jmn
       DO n = 1,ng
#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixend .eq. nxend) ixe = ixend-ixbeg
      
       kzb = 1
       kze = ktile
       if (kzend .eq. nzend) kze = kzend-kzbeg+1
      
!       do k = kzb,kze ; do i = ixb,ixe
!        s0t(i,1-n   ,k) = s0t(i,ny-n,k)
!        s0t(i,ny+n-1,k) = s0t(i,n   ,k)
!       enddo;enddo
#else
        s0t(1:nx-1,1-n, 1:nz) = s0t(1:nx-1,ny-n,1:nz)
        s0t(1:nx-1,ny+n-1,1:nz) = s0t(1:nx-1,n,1:nz)
#endif
      ENDDO
     ENDIF

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn ) jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
             j0     = j 
             jm1    = j - 1 
             jm2    = j - 2 
             jp1    = j + 1 
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,k,j0,jm1,jm2,jp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
    DO k = kmn,kmx1 
     DO j = jmn,jmxb

            j0     = j 
            jm1    = j - 1 
            jm2    = j - 2 
            jp1    = j + 1 

      DO i = imn,imx 
#endif

		dtbydx = dt2*gyt(j,3+js)
		vv     = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
  
         qi   = f60 * (s(i,j,  k) + s(i,j-1,k))  &
              - f61 * (s(i,j+1,k) + s(i,j-2,k))  &
              + f62 * (s(i,j+2,k) + s(i,j-3,k))  
		 
           IF( vv .ge. 0.0 ) THEN
              qD = s0t(i,j0, k)
              qC = s0t(i,jm1,k)
              qU = s0t(i,jm2,k)
           ELSE
              qD = s0t(i,jm1,k)
              qC = s0t(i,j0, k)
              qU = s0t(i,jp1,k)
           ENDIF
	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gyt(j,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		fy(i,j,k)  = vv * qi 

       ENDDO
      ENDDO
     ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! COMPUTE Z-INTERFACE FLUXES

    IF ( nz .gt. 2 ) THEN

!     kmn1 = Max(2,kmn)

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
         km1    = max(k-1,1)
         km2    = max(k-2,1)
         kp1    = min(k+1,k1+ks)
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
     DO k = kmn,kmx1

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,k1+ks)

	 DO j = jmn,jmx 
	  DO i = imn,imx 
#endif
         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

          qi = ( f60 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f61 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f62 * (s(i,j,k+2) + s(i,j,k-3))  &
               - f62 * (s(i,j,k+2) - s(i,j,k-3)   &
               - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
               + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dir )

	     IF( vv .ge. 0.0 ) THEN
		  qD = s0(i,j,k  )
		  qC = s0(i,j,km1)
		  qU = s0(i,j,km2)
	     ELSE
		  qD = s0(i,j,km1)
		  qC = s0(i,j,k  )
		  qU = s0(i,j,kp1)
	     ENDIF
	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gzt(k,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		fz(i,j,k)  = vv * qi

       ENDDO
      ENDDO
     ENDDO
     
    ENDIF

    deallocate( s0t )

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONO: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT_MONO')

    RETURN
    END SUBROUTINE ADVECT_MONO

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_WENO:  5th-order WENO interpolation

    SUBROUTINE ADVECT_WENO()
    
    implicit none

    integer            :: i, im1, jm1, km1
    double precision, parameter    :: f30 =  7.d0/12.d0, f31 = 1.d0/12.d0
    real               :: qim2, qim1, qi, qip1, qip2, vv, dir
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 3

!-----------------------------------------------------------------------------
!  LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

!--------------------------------------------------------------------------
! X-INTERFACE

    CALL cld_cpu('ADVECT-WENO')

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    IF ( nx .gt. 2 ) THEN
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
         vv  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF ( nocollapse .or. ( ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3 ) .or. bcx .eq. 2 ) THEN

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k  )
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fx(i,j,k) = (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF ( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN

         fx(i,j,k) = ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                          - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                          + f31 * (s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)))*dir )

         ELSE

         fx(i,j,k) = 0.5 * (s(i,j,k) + s(im1,j,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Y-INTERFACE

    IF ( ny .gt. 2 ) THEN

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1
    
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      do k = kzb,kze
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe

         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF( nocollapse .or. ( jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3 ) .or. bcy .eq. 2 ) THEN

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fy(i,j,k) = (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF ( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN

         fy(i,j,k) = ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                          - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                          + f31 * (s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)))*dir )

         ELSE

         fy(i,j,k) = 0.5 * (s(i,j,k) + s(i,jm1,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Z-INTERFACE

    IF ( nz .gt. 2 ) THEN
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      do k = kzb,kze
       km1 = k-1
       if (kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe

         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF( kzbeg-1+k .ge. 4 .and. kzbeg-1+k .le. k3 ) THEN

         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,k+1)
            qip1 = s(i,j,k  )
            qi   = s(i,j,k-1)
            qim1 = s(i,j,k-2)
            qim2 = s(i,j,k-3)
          ELSE
            qip2 = s(i,j,k-2)
            qip1 = s(i,j,k-1)
            qi   = s(i,j,k  )
            qim1 = s(i,j,k+1)
            qim2 = s(i,j,k+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fz(i,j,k) = (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF( kzbeg-1+k .eq. 3 .or. kzbeg-1+k .eq. k2 ) THEN

         fz(i,j,k) = ( f30 * (s(i,j,k  ) + s(i,j,k-1))  &
                          - f31 * (s(i,j,k+1) + s(i,j,k-2))  &
                          + f31 * (s(i,j,k+1) - s(i,j,k-2) - 3.*(s(i,j,k)-s(i,j,k-1)))*dir )

         ELSEIF( kzbeg-1+k .eq. 2 .or. kzbeg-1+k .eq. k1 ) THEN

         fz(i,j,k) = 0.5 * (s(i,j,k) + s(i,j,km1))

         ELSE
         
          fz(i,j,k) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENO')

    RETURN
    END SUBROUTINE ADVECT_WENO

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_SLWENO:  5th-order SEMI-LAGRANGIAN WENO interpolation

    SUBROUTINE ADVECT_SLWENO()
    
    implicit none

    integer(kind=4)            :: i, im1, jm1, km1, kshift, j, k, kz
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.
    real(kind=8),    parameter :: f30 =  7.d0/12.d0, f31 = 1.d0/12.d0
    real(kind=4)               :: qim3, qim2, qim1, qi, qip1, qip2, vv, dir
    real(kind=4)               :: h1, h2, h3, cr, crp1, crm1, crl, qtmp, qsum, qres
    real(kind=8)               :: beta0, beta1, beta2, beta3, f0, f1, f2, f3, w1, w2, w3, sumwk
    double precision           :: wi0, wi1, wi2
    real(kind=8),    parameter :: g1 = 1.d0/10.d0, g2 = 6.d0/10.d0, g3 = 3.d0/10.d0, eps=EPSVAL
    double precision, parameter  :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0
    real(kind=8)               :: g11, g22, g33
    integer(kind=4), parameter :: pw5 = 2, pw5b = 3

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

!--------------------------------------------------------------------------
! X-INTERFACE

    CALL cld_cpu('ADVECT-SLWENO')

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    IF ( nx .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1,qip2,qip1,qi,qim1,qim2,f1,f2,f3,beta1,beta2,beta3,w1,w2,sumwk)
      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f1,f2,f3,beta1,beta2,beta3,w1,w2,sumwk)
      DO k = kmn,kmx1
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         vv  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF ( nocollapse .or. ( ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3 ) .or. bcx .eq. 2 ) THEN
#else
         IF ( nocollapse .or. ( i .ge. 4 .and. i .le. i3 ) .or. bcx .eq. 2 ) THEN
#endif
         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k  )
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f1 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f2 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f3 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta1 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta2 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta3 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         w1 = g1 / (eps + beta1)**pw5
         w2 = g2 / (eps + beta2)**pw5
         w3 = g3 / (eps + beta3)**pw5
    
         sumwk = w1 + w2 + w3
    
         fx(i,j,k) = (w1*f1 + w2*f2 + w3*f3) / sumwk

#ifdef MPI
         ELSEIF ( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF ( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
         fx(i,j,k) = ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                     - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                     + f31 * (s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)))*dir )

         ELSE

         fx(i,j,k) = 0.5 * (s(i,j,k) + s(im1,j,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Y-INTERFACE

    IF ( ny .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
      
!      IF ( jye > jtile+1 ) write(0,*) 'advectsl: jye,rank = ',jye,my_rank
    
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1
    
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1,qip2,qip1,qi,qim1,qim2,f1,f2,f3,beta1,beta2,beta3,w1,w2,sumwk)
      do k = kzb,kze
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f1,f2,f3,beta1,beta2,beta3,w1,w2,sumwk)
      DO k = kmn,kmx1
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
!        IF ( bcy .eq. 2 ) jm1 = j-1
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. ( jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3 ) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. ( j .ge. 4 .and. j .le. j3 ) .or. bcy .eq. 2 ) THEN
#endif
         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f1 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f2 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f3 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta1 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta2 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta3 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         w1 = g1 / (eps + beta1)**pw5
         w2 = g2 / (eps + beta2)**pw5
         w3 = g3 / (eps + beta3)**pw5
    
         sumwk = w1 + w2 + w3
    
         fy(i,j,k) = (w1*f1 + w2*f2 + w3*f3) / sumwk

#ifdef MPI
         ELSEIF ( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF ( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
         fy(i,j,k) = ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                     - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                     + f31 * (s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)))*dir )

         ELSE

         fy(i,j,k) = 0.5 * (s(i,j,k) + s(i,jm1,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Z-INTERFACE

    IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile
      if (kzbeg .le. kmn) kzb = Max(kzb, kmn-kzbeg+1 )
      if (kzend .ge. kmx) kze = Min(kze,kmx1-kzbeg+1 )


     ! tests
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. 1+is) ixb = Max( ixb, 1+is-ixbeg+1 )
      if (ixend .ge. nxend-1+is) ixe = nxend-1+is-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. 1+js) jyb = 1+js-jybeg+1
      if (jyend .ge. nyend-1+js) jye = nyend-1+js-jybeg+1

      kzb = 1
      kze = ktile
      if (kzbeg .le. 1+ks) kzb = Max(kzb, 1+ks-kzbeg+1 )
      if (kzend .ge. nzend-1 ) kze = Min(kze,nzend-1 -kzbeg+1 )


!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir, kshift, cr, crm1, crp1, crl, &
!$OMP im1,qip2,qip1,qi,qim1,qim2,qim3,f1,f2,f3,beta1,beta2,beta3,w1,w2,w3, &
!$OMP g11, g22, g33, qtmp, qres, sumwk)
      do k = kzb,kze
       km1 = k-1
       if (kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir, kshift, cr, crm1, crp1, crl, &
!$OMP im1,qip2,qip1,qi,qim1,qim2,qim3,f1,f2,f3,beta1,beta2,beta3,w1,w2,w3, &
!$OMP g11, g22, g33, qtmp, qres, sumwk)
      DO k = kmn,kmx1
       km1 = max(k-1,1)
       DO j = jmn,jmx 
        DO i = imn,imx  
#endif
         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF( kzbeg-1+k .ge. 4 .and. kzbeg-1+k .le. k3 ) THEN
#else
         IF( k .ge. 4 .and. k .le. k3 ) THEN
#endif
           cr     = dtlarge*vv*gzt(k,3+ks)



!           cr     = dt*vv*gzt(k,3+ks)

      ! test:
       
!          IF ( cr >= 2.0 ) THEN
!           crp1     = dtlarge*0.5*(w(i,j,k+1) + w(i-is,j-js,k-ks+1))*gzt(k+1,3+ks)
!           IF ( crp1 < 2.0 ) cr = Min(cr,  1.95)
!          ELSEIF ( cr <= -2.0 ) THEN
!           crm1     = dtlarge*0.5*(w(i,j,k-1) + w(i-is,j-js,k-ks-1))*gzt(k-1,3+ks)
!           IF ( crm1 > -2.0 ) cr = Max(cr, -1.95)
!          ENDIF
!           cr = Min(cr,  1.95)
!           cr = Max(cr, -1.95)

           kshift = floor(abs(cr)) * dir
           crl    = abs(cr - float(kshift))

           IF( abs(kshift) .gt. 2 ) THEN
            write(6,*)  "ADVECT_WENO-SLW5 PROBLEM, BIG COURANT#: (i,j,k,ng)        ", i, j, k, ng
            write(6,*)  "ADVECT_WENO-SLW5 PROBLEM, BIG COURANT#: (is,js,ks)        ", is, js, ks 
            write(6,*)  "ADVECT_WENO-SLW5 PROBLEM, BIG COURANT#: (kmin,k-,k+,kmax) ", -ng+1, k-3-kshift, k+2-kshift, ng
            write(6,*)  "ADVECT_WENO-SLW5 PROBLEM, BIG COURANT#: (cr, crl, kshift) ", cr, crl, kshift
           ENDIF
           
!           IF (  ( Abs( cr ) <= 1.0 .or. Abs( dt*vv*gzt(k,3+ks) ) <= 1 ) ) THEN
           IF (   Abs( cr ) <= 1.0  ) THEN
           
             kz = k
             km1 = max(kz-1,1)

            IF ( dt > dtlarge*0.99 ) THEN ! { assume loop == rkscheme

             IF( dir .ge. 0.0 ) THEN
                 qip2 = s(i,j,kz+1)
                 qip1 = s(i,j,kz  )
                 qi   = s(i,j,kz-1)
                 qim1 = s(i,j,kz-2)
                 qim2 = s(i,j,kz-3)
               ELSE
                 qip2 = s(i,j,kz-2)
                 qip1 = s(i,j,kz-1)
                 qi   = s(i,j,kz  )
                 qim1 = s(i,j,kz+1)
                 qim2 = s(i,j,kz+2)
              ENDIF
         
              f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
              f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
              f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
         
              beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
              beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
              beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
         
              wi0 = gi0 / (eps + beta0)**pw5
              wi1 = gi1 / (eps + beta1)**pw5
              wi2 = gi2 / (eps + beta2)**pw5
         
              sumwk = wi0 + wi1 + wi2
         
               fz(i,j,k) = (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

            
            ELSE ! plain 5th for loop < rkscheme
            
            
            fz(i,j,k)  =  ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
              
            
            ENDIF !}
           
           ELSE

           qim3  = s0(i,j,k-3-kshift)
           qim2  = s0(i,j,k-2-kshift)
           qim1  = s0(i,j,k-1-kshift)
           qi    = s0(i,j,k  -kshift)
           qip1  = s0(i,j,k+1-kshift)
           qip2  = s0(i,j,k+2-kshift)

           IF( cr > 0. ) THEN

            f1 = qim3*(     (crl)/03. -     (crl**2)/02. + (crl**3)/06.) &
               + qim2*(-07.*(crl)/06. + 03.*(crl**2)/02. - (crl**3)/03.) &
               + qim1*( 11.*(crl)/06. -     (crl**2)     + (crl**3)/06.)

            f2 = qim2*(    -(crl)/06.                    + (crl**3)/06.) &
               + qim1*( 05.*(crl)/06. +     (crl**2)/02. - (crl**3)/03.) &
               + qi  *(     (crl)/03. -     (crl**2)/02. + (crl**3)/06.)

            f3 = qim1*(     (crl)/03. +     (crl**2)/02. + (crl**3)/06.) &
               + qi  *( 05.*(crl)/06. -     (crl**2)/02. - (crl**3)/03.) &
               + qip1*(    -(crl)/06.                    + (crl**3)/06.)

            g11 = (2.0 + 03.*crl + crl**2)/20.0
            g22 = (6.0 +     crl - crl**2)/10.0
            g33 = (6.0 - 05.*crl + crl**2)/20.0

            beta1 = 13./12.*(qim3 - 2.*qim2 + qim1)**2 + 1./4.*(   qim3 - 4.*qim2 + 3.*qim1)**2
            beta2 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(   qim2           -    qi  )**2
            beta3 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(3.*qim1 - 4.*qi   +    qip1)**2

           ELSE
  
            f1 = - qim2*(    -(crl)/06.                    + (crl**3)/06.) &
                 - qim1*( 05.*(crl)/06. -     (crl**2)/02. - (crl**3)/03.) &
                 - qi  *(     (crl)/03. +     (crl**2)/02. + (crl**3)/06.)

            f2 = - qim1*(     (crl)/03. -     (crl**2)/02. + (crl**3)/06.) &
                 - qi  *( 05.*(crl)/06. +     (crl**2)/02. - (crl**3)/03.) &
                 - qip1*(    -(crl)/06.                    + (crl**3)/06.)

            f3 = - qi  *( 11.*(crl)/06. -     (crl**2)     + (crl**3)/06.) &
                 - qip1*(-07.*(crl)/06. + 03.*(crl**2)/02. - (crl**3)/03.) &
                 - qip2*(     (crl)/03. -     (crl**2)/02. + (crl**3)/06.)

            g11 = (6.0 - 05.*crl + crl**2)/20.0
            g22 = (6.0 +     crl - crl**2)/10.0
            g33 = (2.0 + 03.*crl + crl**2)/20.0

            beta1 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(   qim2 - 4.*qim1 + 3.*qi  )**2
            beta2 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(   qim1           -    qip1)**2
            beta3 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(3.*qi   - 4.*qip1 +    qip2)**2

          ENDIF

          w1   = g11 / (eps + beta1)**2
          w2   = g22 / (eps + beta2)**2
          w3   = g33 / (eps + beta3)**2

          qres = (w1*f1 + w2*f2 + w3*f3 ) / (w1 + w2 + w3)

          IF( cr .gt. 1.0 ) THEN

            qtmp = SUM(s0(i,j,k-kshift:k-1))
            fz(i,j,k) = (qtmp + qres) / max(eps,cr)

          ELSEIF( cr .lt. -1.0 ) THEN

            qtmp = -SUM(s0(i,j,k:k-kshift-1))
            fz(i,j,k) = (qtmp + qres) / min(eps,cr)

          ELSE

            fz(i,j,k) = qres / (eps+cr)

          ENDIF
          
          ENDIF
          
#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 3 .or. kzbeg-1+k .eq. k2 ) THEN
#else
         ELSEIF( k .eq. 3 .or. k .eq. k2 ) THEN
#endif
           fz(i,j,k) = ( f30 * (s(i,j,k  ) + s(i,j,k-1))  &
                       - f31 * (s(i,j,k+1) + s(i,j,k-2))  &
                       + f31 * (s(i,j,k+1) - s(i,j,k-2) - 3.*(s(i,j,k)-s(i,j,k-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 2 .or. kzbeg-1+k .eq. k1 ) THEN
#else
         ELSEIF( k .eq. 2 .or. k .eq. k1 ) THEN
#endif
           fz(i,j,k) = 0.5 * (s(i,j,k) + s(i,j,km1))

         ELSE
         
           fz(i,j,k) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_SLWENO: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-SLWENO')

    RETURN
    END SUBROUTINE ADVECT_SLWENO

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_WENOZ:  5th-order WENO interpolation; 
! version for enhanced accuracy following Borges et al. (2008) J. Comp. Phys.

    SUBROUTINE ADVECT_WENOZ()
    
    implicit none

    integer            :: i, im1, jm1, km1
    double precision, parameter    :: f30 =  7.d0/12.d0, f31 = 1.d0/12.d0
    real               :: qim2, qim1, qi, qip1, qip2, vv, dir
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk, tau5
    double precision               :: a0, a1, a2, suma
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 2, pw5b = 3

!-----------------------------------------------------------------------------
! LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

!--------------------------------------------------------------------------
! X-INTERFACE

    CALL cld_cpu('ADVECT-WENO')

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    IF ( nx .gt. 2 ) THEN

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)

         vv  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF ( nocollapse .or. ( ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3 ) .or. bcx .eq. 2 ) THEN

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k  )
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
         
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5b)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5b)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5b)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fx(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF ( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN

         fx(i,j,k) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                          - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                          + f31 * (s(i+1,j,k) - s(i-2,j,k) - 3.*(s(i,j,k)-s(i-1,j,k)))*dir )

         ELSE

         fx(i,j,k) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Y-INTERFACE

    IF ( ny .gt. 2 ) THEN

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1
    
      do k = kzb,kze
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe

         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF( nocollapse .or. ( jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3 ) .or. bcy .eq. 2 ) THEN

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5b)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5b)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5b)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fy(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF ( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN

         fy(i,j,k) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                          - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                          + f31 * (s(i,j+1,k) - s(i,j-2,k) - 3.*(s(i,j,k)-s(i,j-1,k)))*dir )

         ELSE

         fy(i,j,k) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
     ENDIF

!--------------------------------------------------------------------------
! Z-INTERFACE

    IF ( nz .gt. 2 ) THEN

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze
       km1 = k-1
       if (kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe

         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

         IF( kzbeg-1+k .ge. 4 .and. kzbeg-1+k .le. k3 ) THEN

         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,k+1)
            qip1 = s(i,j,k  )
            qi   = s(i,j,k-1)
            qim1 = s(i,j,k-2)
            qim2 = s(i,j,k-3)
          ELSE
            qip2 = s(i,j,k-2)
            qip1 = s(i,j,k-1)
            qi   = s(i,j,k  )
            qim1 = s(i,j,k+1)
            qim2 = s(i,j,k+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5b)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5b)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5b)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fz(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ELSEIF( kzbeg-1+k .eq. 3 .or. kzbeg-1+k .eq. k2 ) THEN

         fz(i,j,k) = vv * ( f30 * (s(i,j,k  ) + s(i,j,k-1))  &
                          - f31 * (s(i,j,k+1) + s(i,j,k-2))  &
                          + f31 * (s(i,j,k+1) - s(i,j,k-2) - 3.*(s(i,j,k)-s(i,j,k-1)))*dir )

         ELSEIF( kzbeg-1+k .eq. 2 .or. kzbeg-1+k .eq. k1 ) THEN

         fz(i,j,k) = vv * 0.5 * (s(i,j,k) + s(i,j,km1))

         ELSE
         
          fz(i,j,k) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO
      ENDDO
      
    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENO')

    RETURN
    END SUBROUTINE ADVECT_WENOZ

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE WENOC()

    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

!    real               :: qim2, qim1, qi, qip1, qip2
    double precision   :: qim2, qim1, qi, qip1, qip2
!    real               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
!    real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-8
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 2, pw5b = 3
    integer, parameter :: iweight = 1 ! 0 for old smooth, 1 for Borges et al. 2008

!   double precision :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   double precision :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENOC')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       do j = jyb,jye
        do i = ixb,ixe
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         div2(i,j) = gxt(i,3+is)*(ue-uw)
!          enddo
!          
!          do i = ixb,ixe

        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         vv = uw
         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
         IF( nocollapse .or. (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2) THEN
#endif

!          fx2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
!                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
!                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
!                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
!                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
!                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2

         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2

         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
         ELSE
         
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. (j .ge. 4 .and. j .le. j3) .or. bcy .eq. 2 ) THEN
#endif

!          fy2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
!                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
!                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
!                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
!                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
!                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2

         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
         ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
#else
         IF( kz .ge. 4 .and. kz .le. k3 ) THEN
#endif
!          fz2(i,j,kt) = vv * ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
!                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
!                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
!                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
!                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
!                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

          IF ( Abs(vert_adv_scheme) >= 5 ) THEN

         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

          
          ELSE
          
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
          
          ENDIF

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENOC')

    RETURN
    END SUBROUTINE WENOC

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  7th-order upwind-biased advection

    SUBROUTINE WENO7C()

! 7th-order WENO from Balsara and Shu (2000) J. Comp. Phys.
    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    real               :: qim2, qim1, qi, qip1, qip2, qim3, qip3
    double precision               :: beta0, beta1, beta2, beta3, f0, f1, f2, f3, wi0, wi1, wi2, wi3, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    double precision, parameter    :: c04 = 1.d0/35.d0, c14 = 12.d0/35.d0, c24 = 18.d0/35.d0, c34 = 4.d0/35.d0
    integer, parameter :: pw5 = 2, pw7 = 2, pw5b = 3, pw7b = 4
    integer, parameter :: iweight = 1 ! 0 for old smooth, 1 for Borges et al. 2008


     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
!  LOCAL VARIABLES

    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
   if(debug_mpi)  write(0,"('WENO7C: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENO7C')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      
      IF ( ng < 4 .and. number_of_processes > 1 ) THEN
       write(0,*) 'Increase the number of ghost zones for 7th order WENO! Must have ng=4, but ng = ',ng
       call commasmpi_abort()
      ENDIF
      
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
!      fx2(imx1+1,:) = 0.0
!      fy2(:,jmx1+1) = 0.0
!      fz2(:,:,kt) = 0.0
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,f3,beta0,beta1,beta2,beta3,wi0,wi1,wi2,wi3,sumwk, &
!$OMP         qim3,qim2,qim1,qi,qip1,qip2,qip3)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (ixbeg-1+i .ge. 5 .and. ixbeg-1+i .le. i4) .or. bcx .eq. 2) THEN
#else
         IF( nocollapse .or. (i .ge. 5 .and. i .le. i4) .or. bcx .eq. 2) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip3 = s(i+2,j,k)
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
            qim3 = s(i-4,j,k)
          ELSE
            qip3 = s(i-3,j,k)
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
            qim3 = s(i+3,j,k)
         ENDIF
    
         IF ( .true. ) THEN
         f0 = -1./4.*qim3  + 13./12.*qim2 - 23./12.*qim1 + 25./12.*qi
         f1 =  1./12.*qim2 -  5./12.*qim1 + 13./12.*qi   + 1./4.*qip1
         f2 = -1./12.*qim1 +  7./12.*qi   +  7./12.*qip1 - 1./12.*qip2
         f3 =  1./4.*qi    + 13./12.*qip1 -  5./12.*qip2 + 1./12.*qip3
    
         beta0 = qim3*(547.*qim3 - 3882.*qim2 + 4642.*qim1 - 1854.*qi) +  &
                 qim2*(7043.*qim2 - 17246.*qim1 + 7042.*qi) + qim1*(11003.*qim1 - 9402.*qi) + 2107.*qi**2
         beta1 = qim2*(267.*qim2 - 1642.*qim1 +1602.*qi - 494.*qip1) +   &
                  qim1*(2843.*qim1 - 5966.*qi + 1922.*qip1) + qi*(3443.*qi - 2522.*qip1)+547.*qip1**2
         beta2 = qim1*(547.*qim1 - 2522.*qi + 1922.*qip1 - 494.*qip2) +   &
                 qi*(3443.*qi - 5966.*qip1 + 1602.*qip2) + qip1*(2843.*qip1 - 1642.*qip2) + 267.*qip2**2
         beta3 = qi*(2107.*qi - 9402.*qip1 + 7042.*qip2 - 1854.*qip3) +  &
                 qip1*(11003.*qip1 - 17246.*qip2 + 4642.*qip3) + qip2*(7043.*qip2 - 3882.*qip3) + 547.*qip3**2

         IF ( iweight == 0 ) THEN
           wi0 = c04 / (eps + beta0)**pw7
           wi1 = c14 / (eps + beta1)**pw7
           wi2 = c24 / (eps + beta2)**pw7
           wi3 = c34 / (eps + beta3)**pw7
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c04*(1.0 + (Abs(beta0 - beta3)/(beta0 + eps))**pw7b)
           wi1 = c14*(1.0 + (Abs(beta0 - beta3)/(beta1 + eps))**pw7b)
           wi2 = c24*(1.0 + (Abs(beta0 - beta3)/(beta2 + eps))**pw7b)
           wi3 = c34*(1.0 + (Abs(beta0 - beta3)/(beta3 + eps))**pw7b)
         ENDIF

         sumwk = wi0 + wi1 + wi2 + wi3
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3) / sumwk
         
         ELSE

         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2

         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
         
         ENDIF


#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 4 .or. ixbeg-1+i .eq. i3 ) THEN
#else
         ELSEIF( i .eq. 4 .or. i .eq. i3 ) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
         ELSE
         
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,f3,beta0,beta1,beta2,beta3,wi0,wi1,wi2,wi3,sumwk, &
!$OMP         qim3,qim2,qim1,qi,qip1,qip2,qip3)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (jybeg-1+j .ge. 5 .and. jybeg-1+j .le. j4) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. (j .ge. 5 .and. j .le. j4) .or. bcy .eq. 2 ) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip3 = s(i,j+2,k)
            qip2 = s(i,j+1,k)
            qip1 = s(i,j  ,k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
            qim3 = s(i,j-4,k)
          ELSE
            qip3 = s(i,j-3,k)
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j  ,k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
            qim3 = s(i,j+3,k)
         ENDIF
    
         IF ( .true. ) THEN
         f0 = -1./4.*qim3  + 13./12.*qim2 - 23./12.*qim1 + 25./12.*qi
         f1 =  1./12.*qim2 -  5./12.*qim1 + 13./12.*qi   + 1./4.*qip1
         f2 = -1./12.*qim1 +  7./12.*qi   +  7./12.*qip1 - 1./12.*qip2
         f3 =  1./4.*qi    + 13./12.*qip1 -  5./12.*qip2 + 1./12.*qip3
    
         beta0 = qim3*(547.*qim3 - 3882.*qim2 + 4642.*qim1 - 1854.*qi) +  &
                 qim2*(7043.*qim2 - 17246.*qim1 + 7042.*qi) + qim1*(11003.*qim1 - 9402.*qi) + 2107.*qi**2
         beta1 = qim2*(267.*qim2 - 1642.*qim1 +1602.*qi - 494.*qip1) +   &
                  qim1*(2843.*qim1 - 5966.*qi + 1922.*qip1) + qi*(3443.*qi - 2522.*qip1)+547.*qip1**2
         beta2 = qim1*(547.*qim1 - 2522.*qi + 1922.*qip1 - 494.*qip2) +   &
                 qi*(3443.*qi - 5966.*qip1 + 1602.*qip2) + qip1*(2843.*qip1 - 1642.*qip2) + 267.*qip2**2
         beta3 = qi*(2107.*qi - 9402.*qip1 + 7042.*qip2 - 1854.*qip3) +  &
                 qip1*(11003.*qip1 - 17246.*qip2 + 4642.*qip3) + qip2*(7043.*qip2 - 3882.*qip3) + 547.*qip3**2
    
          IF ( iweight == 0 ) THEN
           wi0 = c04 / (eps + beta0)**pw7
           wi1 = c14 / (eps + beta1)**pw7
           wi2 = c24 / (eps + beta2)**pw7
           wi3 = c34 / (eps + beta3)**pw7
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c04*(1.0 + (Abs(beta0 - beta3)/(beta0 + eps))**pw7b)
           wi1 = c14*(1.0 + (Abs(beta0 - beta3)/(beta1 + eps))**pw7b)
           wi2 = c24*(1.0 + (Abs(beta0 - beta3)/(beta2 + eps))**pw7b)
           wi3 = c34*(1.0 + (Abs(beta0 - beta3)/(beta3 + eps))**pw7b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2 + wi3
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3) / sumwk
         
         ELSE
         
         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
         
         ENDIF

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 4 .or. jybeg-1+j .eq. j3 ) THEN
#else
         ELSEIF( j .eq. 4 .or. j .eq. j3 ) THEN
#endif


         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
         ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,f3,beta0,beta1,beta2,beta3,wi0,wi1,wi2,wi3,sumwk, &
!$OMP         qim3,qim2,qim1,qi,qip1,qip2,qip3)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .gt. 4 .and. kzbeg-1+kz .lt. k3 ) THEN
#else
         IF( kz .gt. 4 .and. kz .lt. k3 ) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip3 = s(i,j,kz+2)
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
            qim3 = s(i,j,kz-4)
          ELSE
            qip3 = s(i,j,kz-3)
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
            qim3 = s(i,j,kz+3)
         ENDIF
    
         f0 = -1./4.*qim3  + 13./12.*qim2 - 23./12.*qim1 + 25./12.*qi
         f1 =  1./12.*qim2 -  5./12.*qim1 + 13./12.*qi   + 1./4.*qip1
         f2 = -1./12.*qim1 +  7./12.*qi   +  7./12.*qip1 - 1./12.*qip2
         f3 =  1./4.*qi    + 13./12.*qip1 -  5./12.*qip2 + 1./12.*qip3
    
         beta0 = qim3*(547.*qim3 - 3882.*qim2 + 4642.*qim1 - 1854.*qi) +  &
                 qim2*(7043.*qim2 - 17246.*qim1 + 7042.*qi) + qim1*(11003.*qim1 - 9402.*qi) + 2107.*qi**2
         beta1 = qim2*(267.*qim2 - 1642.*qim1 +1602.*qi - 494.*qip1) +   &
                  qim1*(2843.*qim1 - 5966.*qi + 1922.*qip1) + qi*(3443.*qi - 2522.*qip1)+547.*qip1**2
         beta2 = qim1*(547.*qim1 - 2522.*qi + 1922.*qip1 - 494.*qip2) +   &
                 qi*(3443.*qi - 5966.*qip1 + 1602.*qip2) + qip1*(2843.*qip1 - 1642.*qip2) + 267.*qip2**2
         beta3 = qi*(2107.*qi - 9402.*qip1 + 7042.*qip2 - 1854.*qip3) +  &
                 qip1*(11003.*qip1 - 17246.*qip2 + 4642.*qip3) + qip2*(7043.*qip2 - 3882.*qip3) + 547.*qip3**2

         IF ( iweight == 0 ) THEN
           wi0 = c04 / (eps + beta0)**pw7
           wi1 = c14 / (eps + beta1)**pw7
           wi2 = c24 / (eps + beta2)**pw7
           wi3 = c34 / (eps + beta3)**pw7
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c04*(1.0 + (Abs(beta0 - beta3)/(beta0 + eps))**pw7b)
           wi1 = c14*(1.0 + (Abs(beta0 - beta3)/(beta1 + eps))**pw7b)
           wi2 = c24*(1.0 + (Abs(beta0 - beta3)/(beta2 + eps))**pw7b)
           wi3 = c34*(1.0 + (Abs(beta0 - beta3)/(beta3 + eps))**pw7b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2 + wi3

         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3) / sumwk

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 4 .or. kzbeg-1+kz .eq. k3 ) THEN
#else
         ELSEIF( kz .eq. 4 .or. kz .eq. k3 ) THEN
#endif

!          IF ( Abs(vert_adv_scheme) == 5 ) THEN
!
         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
!
!          
!          ELSE
!          
!          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
!                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
!                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
!                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
!          
!          ENDIF

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('WENO7C: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENO7C')

    RETURN
    END SUBROUTINE WENO7C
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  9th-order upwind-biased advection

    SUBROUTINE WENO9C()

! Horizontal 9th-order WENO from Balsara and Shu (2000) J. Comp. Phys.
! Vertical is 9th order
    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    double precision   :: qim4, qim3, qim2, qim1, qi, qip1, qip2, qip3, qip4
    double precision   :: beta0, beta1, beta2, beta3, beta4
    double precision   ::  f0, f1, f2, f3, f4, wi0, wi1, wi2, wi3, wi4, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0
    double precision, parameter    :: eps=EPSVAL
    double precision   :: eps9
    double precision, parameter    :: c04 = 1.d0/35.d0, c14 = 12.d0/35.d0, c24 = 18.d0/35.d0, c34 = 4.d0/35.d0
    double precision, parameter    :: c05 = 1.d0/126.d0, c15 = 20.d0/126.d0, c25 = 60.d0/126.d0
    double precision, parameter    :: c35 = 40.d0/126.d0, c45 = 5.d0/126.d0
    integer, parameter :: pw5 = 2, pw7 = 2, pw9 = 2, pw5b = 3, pw7b = 4, pw9b = 2
    integer, parameter :: iweight = 1 ! 0 for old smooth, 1 for Borges et al. 2008
    
    real :: weno9

     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
!  LOCAL VARIABLES

    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
   if(debug_mpi)  write(0,"('WENO9C: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENO9C')
    
    eps9=EPSVAL9
! Attempt to get pw9b=5 to work, but still blows up, so stick with 2 for now
!     IF ( pw9b > 2 ) THEN
!     IF ( trim(varname) == 'U' .or. trim(varname) == 'V' .or. trim(varname) == 'W' &
!          .or. trim(varname) == 'TH' ) THEN 
!        eps9 = 1.d-6
!     ELSEIF ( trim(varname) == 'KM' .or. trim(varname) == 'TKE'  .or. trim(varname) == 'SSMX' ) THEN
!        eps9 = 1.d-6
!     ELSEIF ( varname(1:1) == 'Q' ) THEN
!         eps9 = 1.d-16
!     ELSEIF ( varname(1:1) == 'C' ) THEN
!         eps9 = 1.d-10
!     ELSEIF ( varname(1:1) == 'Z' ) THEN
!         eps9 = 1.d-30
!     ELSEIF ( varname(1:1) == 'V' ) THEN
!         eps9 = 1.d-16
!     ENDIF
!     ENDIF


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      
      IF ( ng < 5 .and. number_of_processes > 1 ) THEN
       write(0,*) 'Increase the number of ghost zones for 9th order WENO! Must have ng=5, but ng = ',ng
       call commasmpi_abort()
      ENDIF
      
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
!      fx2(imx1+1,:) = 0.0
!      fy2(:,jmx1+1) = 0.0
!      fz2(:,:,kt) = 0.0
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,f3,f4,beta0,beta1,beta2,beta3,beta4,wi0,wi1,wi2,wi3,wi4,sumwk, &
!$OMP         qim4,qim3,qim2,qim1,qi,qip1,qip2,qip3,qip4)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (ixbeg-1+i .ge. 6 .and. ixbeg-1+i .le. i5) .or. bcx .eq. 2) THEN
#else
         IF( nocollapse .or. (i .ge. 6 .and. i .le. i5) .or. bcx .eq. 2) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip4 = s(i+3,j,k)
            qip3 = s(i+2,j,k)
            qip2 = s(i+1,j,k)
            qip1 = s(i  ,j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
            qim3 = s(i-4,j,k)
            qim4 = s(i-5,j,k)
          ELSE
            qip4 = s(i-4,j,k)
            qip3 = s(i-3,j,k)
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i  ,j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
            qim3 = s(i+3,j,k)
            qim4 = s(i+4,j,k)
         ENDIF
    
         IF ( .true. ) THEN

!         IF ( .true. ) THEN
         f0 =   1./5.*qim4  - 21./20.*qim3 + 137./60.*qim2 - 163./60.*qim1  + 137./60.*qi
         f1 = -1./20.*qim3  + 17./60.*qim2 -  43./60.*qim1 +  77./60.*qi    +    1./5.*qip1
         f2 =  1./30.*qim2  - 13./60.*qim1 +  47./60.*qi   +   9./20.*qip1  -   1./20.*qip2
         f3 = -1./20.*qim1  +  9./20.*qi   +  47./60.*qip1 -  13./60.*qip2  +   1./30.*qip3
         f4 =   1./5.*qi    + 77./60.*qip1 -  43./60.*qip2 +  17./60.*qip3  -   1./20.*qip4


      beta0 = qim4*(22658.0*qim4  - 208501.0*qim3 + 364863.0*qim2  - 288007.0*qim1 + 86329.0*qi ) &
            + qim3*(482963.0*qim3  - 1704396.0*qim2 + 1358458.0*qim1 - 411487.0*qi ) &
            + qim2*(1521393.0*qim2 - 2462076.0*qim1 + 758823.0*qi )  &
            + qim1*(1020563.0*qim1  - 649501.0*qi ) + 107918.0*qi**2

      beta1 = qim3*(6908.0*qim3 - 60871.0*qim2 + 99213.0*qim1 - 70237.0*qi + 18079.0*qip1) &
            + qim2*(138563.0*qim2 - 464976.0*qim1 + 337018.0*qi - 88297.0*qip1) &
            + qim1*(406293.0*qim1 - 611976.0*qi + 165153.0*qip1) &
            + qi  *(242723.0*qi - 140251.0*qip1 ) + 22658.0*qip1**2

     beta2 = qim2* (6908.0*qim2 - 51001.0*qim1 + 67923.0*qi - 38947.0*qip1 + 8209.0*qip2) &
           + qim1*(104963.0*qim1 - 299076.0*qi + 179098.0*qip1 - 38947.0*qip2) &
           + qi  *(231153.0*qi - 299076.0*qip1 + 67923.0*qip2) &
           + qip1*(104963.0*qip1 - 51001.0*qip2 ) + 6908.0*qip2**2

     beta3 = qim1* (22658.0*qim1 - 140251.0*qi + 165153.0*qip1 - 88297.0*qip2 + 18079.0*qip3 )  &
           + qi * (242723.0*qi - 611976.0*qip1 + 337018.0*qip2 - 70237.0*qip3) &
           + qip1*(406293.0*qip1 - 464976.0*qip2 + 99213.0*qip3) &
           + qip2*(138563.0*qip2 - 60871.0*qip3 ) + 6908.0*qip3**2

     beta4 = qi*  (107918.0*qi - 649501.0*qip1 + 758823.0*qip2 - 411487.0*qip2 + 86329.0*qip4) &
           + qip1*(1020563.0*qip1 - 2462076.0*qip2 + 1358458.0*qip3 - 288007.0*qip4) &
           + qip2*(1521393.0*qip2 - 1704396.0*qip3 + 364863.0*qip4) &
           + qip3*(482963.0*qip3 - 208501.0*qip4 ) + 22658.0*qip4**2

         IF ( iweight == 0 ) THEN
           wi0 = c05 * Min( EPSWVAL, (eps + beta0)**(-pw9))
           wi1 = c15 * Min( EPSWVAL, (eps + beta1)**(-pw9))
           wi2 = c25 * Min( EPSWVAL, (eps + beta2)**(-pw9))
           wi3 = c35 * Min( EPSWVAL, (eps + beta3)**(-pw9))
           wi4 = c45 * Min( EPSWVAL, (eps + beta4)**(-pw9))
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c05*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta0 + eps9))**pw9b)
           wi1 = c15*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta1 + eps9))**pw9b)
           wi2 = c25*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta2 + eps9))**pw9b)
           wi3 = c35*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta3 + eps9))**pw9b)
           wi4 = c45*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta4 + eps9))**pw9b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2 + wi3 + wi4
    
         fx2(i,j) = vv * ((wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3 + wi4*f4) / sumwk)
         
!           ELSE
!           
!              fx2(i,j) = vv *  weno9(qim4,qim3,qim2,qim1,qi,qip1,qip2,qip3,qip4,eps)
!           ENDIF
         
         ELSE
! 5th order
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
         
         ENDIF


         ELSEIF( ixbeg-1+i .eq. 4 .or. ixbeg-1+i .eq. i3 .or. ixbeg-1+i .eq. 5 .or. ixbeg-1+i .eq. i4 ) THEN

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
         ELSE
         
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (jybeg-1+j .ge. 6 .and. jybeg-1+j .le. j3-2) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. (j .ge. 6 .and. j .le. j3-2) .or. bcy .eq. 2 ) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip4 = s(i,j+3,k)
            qip3 = s(i,j+2,k)
            qip2 = s(i,j+1,k)
            qip1 = s(i,j  ,k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
            qim3 = s(i,j-4,k)
            qim4 = s(i,j-5,k)
          ELSE
            qip4 = s(i,j-4,k)
            qip3 = s(i,j-3,k)
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j  ,k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
            qim3 = s(i,j+3,k)
            qim4 = s(i,j+4,k)
         ENDIF
    
         IF ( .true. ) THEN
         f0 =   1./5.*qim4  - 21./20.*qim3 + 137./60.*qim2 - 163./60.*qim1  + 137./60.*qi
         f1 = -1./20.*qim3  + 17./60.*qim2 -  43./60.*qim1 +  77./60.*qi    +    1./5.*qip1
         f2 =  1./30.*qim2  - 13./60.*qim1 +  47./60.*qi   +   9./20.*qip1  -   1./20.*qip2
         f3 = -1./20.*qim1  +  9./20.*qi   +  47./60.*qip1 -  13./60.*qip2  +   1./30.*qip3
         f4 =   1./5.*qi    + 77./60.*qip1 -  43./60.*qip2 +  17./60.*qip3  -   1./20.*qip4

      beta0 = qim4*(22658.0*qim4  - 208501.0*qim3 + 364863.0*qim2  - 288007.0*qim1 + 86329.0*qi ) &
            + qim3*(482963.0*qim3  - 1704396.0*qim2 + 1358458.0*qim1 - 411487.0*qi ) &
            + qim2*(1521393.0*qim2 - 2462076.0*qim1 + 758823.0*qi )  &
            + qim1*(1020563.0*qim1  - 649501.0*qi ) + 107918.0*qi**2

      beta1 = qim3*(6908.0*qim3 - 60871.0*qim2 + 99213.0*qim1 - 70237.0*qi + 18079.0*qip1) &
            + qim2*(138563.0*qim2 - 464976.0*qim1 + 337018.0*qi - 88297.0*qip1) &
            + qim1*(406293.0*qim1 - 611976.0*qi + 165153.0*qip1) &
            + qi  *(242723.0*qi - 140251.0*qip1 ) + 22658.0*qip1**2

     beta2 = qim2* (6908.0*qim2 - 51001.0*qim1 + 67923.0*qi - 38947.0*qip1 + 8209.0*qip2) &
           + qim1*(104963.0*qim1 - 299076.0*qi + 179098.0*qip1 - 38947.0*qip2) &
           + qi  *(231153.0*qi - 299076.0*qip1 + 67923.0*qip2) &
           + qip1*(104963.0*qip1 - 51001.0*qip2 ) + 6908.0*qip2**2

     beta3 = qim1* (22658.0*qim1 - 140251.0*qi + 165153.0*qip1 - 88297.0*qip2 + 18079.0*qip3 )  &
           + qi * (242723.0*qi - 611976.0*qip1 + 337018.0*qip2 - 70237.0*qip3) &
           + qip1*(406293.0*qip1 - 464976.0*qip2 + 99213.0*qip3) &
           + qip2*(138563.0*qip2 - 60871.0*qip3 ) + 6908.0*qip3**2

     beta4 = qi*  (107918.0*qi - 649501.0*qip1 + 758823.0*qip2 - 411487.0*qip2 + 86329.0*qip4) &
           + qip1*(1020563.0*qip1 - 2462076.0*qip2 + 1358458.0*qip3 - 288007.0*qip4) &
           + qip2*(1521393.0*qip2 - 1704396.0*qip3 + 364863.0*qip4) &
           + qip3*(482963.0*qip3 - 208501.0*qip4 ) + 22658.0*qip4**2

    
         IF ( iweight == 0 ) THEN
           wi0 = c05 * Min( EPSWVAL, (eps + beta0)**(-pw9))
           wi1 = c15 * Min( EPSWVAL, (eps + beta1)**(-pw9))
           wi2 = c25 * Min( EPSWVAL, (eps + beta2)**(-pw9))
           wi3 = c35 * Min( EPSWVAL, (eps + beta3)**(-pw9))
           wi4 = c45 * Min( EPSWVAL, (eps + beta4)**(-pw9))
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c05*(1.0 + (Abs(beta0 - beta4)/(beta0 + eps9))**pw9b)
           wi1 = c15*(1.0 + (Abs(beta0 - beta4)/(beta1 + eps9))**pw9b)
           wi2 = c25*(1.0 + (Abs(beta0 - beta4)/(beta2 + eps9))**pw9b)
           wi3 = c35*(1.0 + (Abs(beta0 - beta4)/(beta3 + eps9))**pw9b)
           wi4 = c45*(1.0 + (Abs(beta0 - beta4)/(beta4 + eps9))**pw9b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2 + wi3 + wi4
    
         fy2(i,j) = vv * ( (wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3 + wi4*f4) / sumwk )
         
         ELSE
         
         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
         
         ENDIF

         ELSEIF( jybeg-1+j .eq. 4 .or. jybeg-1+j .eq. j3 .or. jybeg-1+j .eq. 5 .or. jybeg-1+j .eq. j4 ) THEN


         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
         ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

!#ifdef MPI
!         IF( nocollapse .or. (jybeg-1+j .ge. 6 .and. jybeg-1+j .le. j3-2) .or. bcy .eq. 2 ) THEN
!#else
!         IF( nocollapse .or. (j .ge. 6 .and. j .le. j3-2) .or. bcy .eq. 2 ) THEN
!#endif

         IF( kzbeg-1+kz .ge. 6 .and. kzbeg-1+kz .le. k4 - 1 ) THEN
         
         IF ( .true. .and. vert_adv_scheme >= 5 ) THEN
         
         IF ( dir .ge. 0.0 ) THEN
            qip4 = s(i,j,kz+3)
            qip3 = s(i,j,kz+2)
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
            qim3 = s(i,j,kz-4)
            qim4 = s(i,j,kz-5)
          ELSE
            qip4 = s(i,j,kz-4)
            qip3 = s(i,j,kz-3)
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
            qim3 = s(i,j,kz+3)
            qim4 = s(i,j,kz+4)
         ENDIF
    
!         IF ( .true. ) THEN
         f0 =   1./5.*qim4  - 21./20.*qim3 + 137./60.*qim2 - 163./60.*qim1  + 137./60.*qi
         f1 = -1./20.*qim3  + 17./60.*qim2 -  43./60.*qim1 +  77./60.*qi    +    1./5.*qip1
         f2 =  1./30.*qim2  - 13./60.*qim1 +  47./60.*qi   +   9./20.*qip1  -   1./20.*qip2
         f3 = -1./20.*qim1  +  9./20.*qi   +  47./60.*qip1 -  13./60.*qip2  +   1./30.*qip3
         f4 =   1./5.*qi    + 77./60.*qip1 -  43./60.*qip2 +  17./60.*qip3  -   1./20.*qip4


      beta0 = qim4*(22658.0*qim4  - 208501.0*qim3 + 364863.0*qim2  - 288007.0*qim1 + 86329.0*qi ) &
            + qim3*(482963.0*qim3  - 1704396.0*qim2 + 1358458.0*qim1 - 411487.0*qi ) &
            + qim2*(1521393.0*qim2 - 2462076.0*qim1 + 758823.0*qi )  &
            + qim1*(1020563.0*qim1  - 649501.0*qi ) + 107918.0*qi**2

      beta1 = qim3*(6908.0*qim3 - 60871.0*qim2 + 99213.0*qim1 - 70237.0*qi + 18079.0*qip1) &
            + qim2*(138563.0*qim2 - 464976.0*qim1 + 337018.0*qi - 88297.0*qip1) &
            + qim1*(406293.0*qim1 - 611976.0*qi + 165153.0*qip1) &
            + qi  *(242723.0*qi - 140251.0*qip1 ) + 22658.0*qip1**2

     beta2 = qim2* (6908.0*qim2 - 51001.0*qim1 + 67923.0*qi - 38947.0*qip1 + 8209.0*qip2) &
           + qim1*(104963.0*qim1 - 299076.0*qi + 179098.0*qip1 - 38947.0*qip2) &
           + qi  *(231153.0*qi - 299076.0*qip1 + 67923.0*qip2) &
           + qip1*(104963.0*qip1 - 51001.0*qip2 ) + 6908.0*qip2**2

     beta3 = qim1* (22658.0*qim1 - 140251.0*qi + 165153.0*qip1 - 88297.0*qip2 + 18079.0*qip3 )  &
           + qi * (242723.0*qi - 611976.0*qip1 + 337018.0*qip2 - 70237.0*qip3) &
           + qip1*(406293.0*qip1 - 464976.0*qip2 + 99213.0*qip3) &
           + qip2*(138563.0*qip2 - 60871.0*qip3 ) + 6908.0*qip3**2

     beta4 = qi*  (107918.0*qi - 649501.0*qip1 + 758823.0*qip2 - 411487.0*qip2 + 86329.0*qip4) &
           + qip1*(1020563.0*qip1 - 2462076.0*qip2 + 1358458.0*qip3 - 288007.0*qip4) &
           + qip2*(1521393.0*qip2 - 1704396.0*qip3 + 364863.0*qip4) &
           + qip3*(482963.0*qip3 - 208501.0*qip4 ) + 22658.0*qip4**2

         IF ( iweight == 0 ) THEN
           wi0 = c05 / (eps + beta0)**pw9
           wi1 = c15 / (eps + beta1)**pw9
           wi2 = c25 / (eps + beta2)**pw9
           wi3 = c35 / (eps + beta3)**pw9
           wi4 = c45 / (eps + beta4)**pw9
         ELSEIF ( iweight == 1 ) THEN
           wi0 = c05*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta0 + eps9))**pw9b)
           wi1 = c15*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta1 + eps9))**pw9b)
           wi2 = c25*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta2 + eps9))**pw9b)
           wi3 = c35*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta3 + eps9))**pw9b)
           wi4 = c45*(1.0 + Min( EPSWVAL, Abs(beta0 - beta4)/(beta4 + eps9))**pw9b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2 + wi3 + wi4
    
         fz2(i,j,kt) = vv * ( (wi0*f0 + wi1*f1 + wi2*f2 + wi3*f3 + wi4*f4) / sumwk )
         
!           ELSE
!           
!              fz2(i,j,kt) = vv *  weno9(qim4,qim3,qim2,qim1,qi,qip1,qip2,qip3,qip4,eps)
!           ENDIF
         
         ELSE ! vert_adv_scheme
         
         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
         
         ENDIF ! vert_adv_scheme


         ELSEIF( kzbeg-1+kz .eq. 4 .or. kzbeg-1+kz .eq. k3 .or. kzbeg-1+kz .eq. 5 .or. kzbeg-1+kz .eq. k4) THEN


         IF ( .false. ) THEN ! plain 5th-order -- does not seem to make much difference using this or the WENO 5th
         
          fz2(i,j,kt)  =  vv * ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
         
         ELSE ! 5th-order WENO
         
         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

         ENDIF
         
         

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('WENO9C: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENO9C')

    RETURN
    END SUBROUTINE WENO9C
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  7th-order upwind-biased advection

    SUBROUTINE ADVECT9()

! Horizontal 9th-order upwind from Li (1997) J. Comp. Phys. 133, 235-255
! Vertical is 5th order
    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv, qq
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    real               :: qim4, qim3, qim2, qim1, qi, qip1, qip2, qip3, qip4
    double precision   :: beta0, beta1, beta2, beta3, beta4
    double precision   ::  f0, f1, f2, f3, f4, wi0, wi1, wi2, wi3, wi4, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    double precision, parameter    :: c04 = 1.d0/35.d0, c14 = 12.d0/35.d0, c24 = 18.d0/35.d0, c34 = 4.d0/35.d0
    double precision, parameter    :: c05 = 1.d0/126.d0, c15 = 10.d0/63.d0, c25 = 10.d0/21.d0
    double precision, parameter    :: c35 = 20.d0/63.d0, c45 = 5.d0/126.d0
    integer, parameter :: pw = 2


     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
!  LOCAL VARIABLES

    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT9: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-ADVECT9')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      
      IF ( ng < 5 .and. number_of_processes > 1 ) THEN
       write(0,*) 'Increase the number of ghost zones for 9th order! Must have ng=5, but ng = ',ng
       call commasmpi_abort()
      ENDIF
      
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
!      fx2(imx1+1,:) = 0.0
!      fy2(:,jmx1+1) = 0.0
!      fz2(:,:,kt) = 0.0
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,f3,f4,beta0,beta1,beta2,beta3,beta4,wi0,wi1,wi2,wi3,wi4,sumwk, &
!$OMP         qim4,qim3,qim2,qim1,qi,qip1,qip2,qip3,qip4)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (ixbeg-1+i .ge. 6 .and. ixbeg-1+i .le. i5) .or. bcx .eq. 2) THEN
#else
         IF( nocollapse .or. (i .ge. 6 .and. i .le. i5) .or. bcx .eq. 2) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip4 = s(i+3,j,k)
            qip3 = s(i+2,j,k)
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
            qim3 = s(i-4,j,k)
            qim4 = s(i-5,j,k)
          ELSE
            qip4 = s(i-4,j,k)
            qip3 = s(i-3,j,k)
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
            qim3 = s(i+3,j,k)
            qim4 = s(i+4,j,k)
         ENDIF
    

          qq =  1./362880.*(576*qim4 - 5904*qim4  + 28656*qim2  &
                - 92304*qim1 + 270576*qi + 198000*qip1 - 43920*qip2 + 7920*qip3 - 720*qip4)
         
         fx2(i,j) = vv * qq


         ELSEIF( ixbeg-1+i .eq. 4 .or. ixbeg-1+i .eq. i3 .or. ixbeg-1+i .eq. 5 .or. ixbeg-1+i .eq. i4 ) THEN

          qq =  ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )
    
         fx2(i,j) = vv * qq

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
         ELSE
         
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (jybeg-1+j .ge. 6 .and. jybeg-1+j .le. j3-2) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. (j .ge. 6 .and. j .le. j3-2) .or. bcy .eq. 2 ) THEN
#endif

         IF ( dir .ge. 0.0 ) THEN
            qip4 = s(i,j+3,k)
            qip3 = s(i,j+2,k)
            qip2 = s(i,j+1,k)
            qip1 = s(i,j  ,k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
            qim3 = s(i,j-4,k)
            qim4 = s(i,j-5,k)
          ELSE
            qip4 = s(i,j-4,k)
            qip3 = s(i,j-3,k)
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j  ,k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
            qim3 = s(i,j+3,k)
            qim4 = s(i,j+4,k)
         ENDIF
    
          qq =  1./362880.*(576*qim4 - 5904*qim4  + 28656*qim2  &
                - 92304*qim1 + 270576*qi + 198000*qip1 - 43920*qip2 + 7920*qip3 - 720*qip4)
         
         fy2(i,j) = vv * qq

         ELSEIF( jybeg-1+j .eq. 4 .or. jybeg-1+j .eq. j3 .or. jybeg-1+j .eq. 5 .or. jybeg-1+j .eq. j4 ) THEN


          qq =     ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )
    
         fy2(i,j) = vv * qq

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
         ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
#else
         IF( kz .ge. 4 .and. kz .le. k3 ) THEN
#endif

          qq =  ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
    
         fz2(i,j,kt) = vv * qq

          

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT9: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-ADVECT9')

    RETURN
    END SUBROUTINE ADVECT9
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE WENOCZ()

    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    real               :: qim2, qim1, qi, qip1, qip2
!    real               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
!    real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-8
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    double precision               :: a0, a1, a2, suma, tau5
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 2

!   double precision :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   double precision :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

#ifdef MPI
   if(debug_mpi)  write(0,"('WENOCZ: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENOCZ')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
!      fx2(imx1+1,:) = 0.0
!      fy2(:,jmx1+1) = 0.0
!      fz2(:,:,kt) = 0.0
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
         IF( nocollapse .or. (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2) THEN
#endif

!          fx2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
!                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
!                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
!                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
!                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
!                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
         ELSE
         
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( nocollapse .or. (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
         IF( nocollapse .or. (j .ge. 4 .and. j .le. j3) .or. bcy .eq. 2 ) THEN
#endif

!          fy2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
!                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
!                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
!                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
!                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
!                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
         ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
#else
         IF( kz .ge. 4 .and. kz .le. k3 ) THEN
#endif
!          fz2(i,j,kt) = vv * ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
!                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
!                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
!                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
!                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
!                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

          IF ( Abs(vert_adv_scheme) == 5 ) THEN

         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         tau5 = abs(beta0-beta2)
         
         a0 = gi0*(1.+(tau5/(beta0+eps))**pw5)
         a1 = gi1*(1.+(tau5/(beta1+eps))**pw5)
         a2 = gi2*(1.+(tau5/(beta2+eps))**pw5)

         suma = a0 + a1 + a2

!  set w0 through w2 and sum in sumw

         wi0 = a0 / suma
         wi1 = a1 / suma
         wi2 = a2 / suma
    
!         wi0 = gi0 / (eps + beta0)**pw5
!         wi1 = gi1 / (eps + beta1)**pw5
!         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

          
          ELSE
          
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
          
          ENDIF

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENOCZ')

    RETURN
    END SUBROUTINE WENOCZ

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM WENO5C:  5th-order upwind-biased WENO advection

    SUBROUTINE WENO5C()

    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

!    real               :: qim2, qim1, qi, qip1, qip2
    double precision   :: qim2, qim1, qi, qip1, qip2
!    real               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
!    real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-8
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 2, pw5b = 3
    integer, parameter :: iweight = 1 ! 0 for old smooth, 1 for Borges et al. 2008

!   double precision :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
!   double precision :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   double precision :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze

#ifdef MPI
  if(debug_mpi) write(0,"('WENO5C: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENO5C')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = Max(kzb,kmn-kzbeg+1)
      if (kzend .ge. kmx) kze = Min(kze,kmx1-kzbeg+1)

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

      
    IF ( nx .gt. 2  .and. k > 0 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,ue,uw,im1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
       do j = jyb,jye
!DIR$ IVDEP
        DO i = ixb,ixe
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))

         div2(i,j) = gxt(i,3+is)*(ue-uw)

        ENDDO
!DIR$ IVDEP
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
         km1 = max(k-1,1)
         vv  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))

         dir = sign(1.0,vv)

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k)
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2

         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2

         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1) 
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vs,vn,vv,jm1,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)

!DIR$ IVDEP
      DO i = ixb,ixe
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)
      ENDDO
!DIR$ IVDEP
      do i = ixb,ixe
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         dir = sign(1.0,vv)

         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2

         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF

         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
!$OMP PARALLEL DO DEFAULT(SHARED),  &
!$OMP PRIVATE(i,j,vv,wt,wb,dir,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk, &
!$OMP         qim2,qim1,qi,qip1,qip2)
      do j = jyb,jye
!DIR$ IVDEP
      DO i = ixb,ixe
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))

        wt  = vv*rrp(kz-1)
        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)
      ENDDO
!DIR$ IVDEP
      do i = ixb,ixe
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN

          IF ( Abs(vert_adv_scheme) >= 5 ) THEN

         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         IF ( iweight == 0 ) THEN
           wi0 = gi0 / (eps + beta0)**pw5
           wi1 = gi1 / (eps + beta1)**pw5
           wi2 = gi2 / (eps + beta2)**pw5
         ELSEIF ( iweight == 1 ) THEN
           wi0 = gi0*(1.0 + (Abs(beta0 - beta2)/(beta0 + eps))**pw5b)
           wi1 = gi1*(1.0 + (Abs(beta0 - beta2)/(beta1 + eps))**pw5b)
           wi2 = gi2*(1.0 + (Abs(beta0 - beta2)/(beta2 + eps))**pw5b)
         ENDIF
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

          
          ELSE
          
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )
          
          ENDIF

         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN

          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN

          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   do j = jyb,jye
    do i = ixb,ixe
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult

       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('WENO5C: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENO5C')

    RETURN
    END SUBROUTINE WENO5C

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT_TEST:  Test program which uses the 1D agorithms
!          This program can be used to test a new algorithm which initially
!          only needs a 1D representation.  Example routines are at the bottom
!          of the file.
!
!    SUBROUTINE ADVECT_TEST()
!    
!      real, dimension(:), allocatable :: u1, q1, f1, dtbydx
!
!      CALL cld_cpu('ADVECT-TEST')
!
!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES
!
!    IF ( nx .gt. 2 ) THEN
!    
!      allocate(u1(nx+is))
!      allocate(f1(nx+is))
!      allocate(q1(nx-1+is))
!      allocate(dtbydx(nx+is))
!
!      u1(:)           = 0.0
!      f1(:)           = 0.0
!      dtbydx(:)       = 0.0
!      dtbydx(imn:imx) = dt * gx(3+is)%flt1d(imn:imx)
!    
!      DO k = kmn,kmx
!       DO j = jmn,jmx 
!       
!        DO i = 1,nx-1+is  
!
!         u1(i)     = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
!         q1(i)     = s(i,j,k)
!         
!        ENDDO
!        
!       CALL ADVECT5_1D(q1,u1,f1,imn,imx,0,nx-1+is)
!       CALL ADVECT6_FILTER6_1D(q1,u1,f1,dtbydx,0.05,imn,imx,0,nx-1+is)
!       CALL ADVECT6_MONO_1D(q1,u1,f1,dtbydx,imn,imx,0,nx-1+is)
!       CALL ADVECT_WENO_1D(q1,u1,f1,imn,imx,0,nx-1+is)
!
!        fx(imn:imx,j,k) = f1(imn:imx)
!       
!       ENDDO
!      ENDDO
!
!      deallocate(u1)
!      deallocate(f1)
!      deallocate(q1)
!      deallocate(dtbydx)
!
!    ENDIF
!
!-----------------------------------------------------------------------------
! COMPUTE Y-INTERFACE FLUXES
!
!    IF( ny .gt. 2 ) THEN
!    
!      allocate(u1(ny+js))
!      allocate(f1(ny+js))
!      allocate(q1(ny-1+js))
!      allocate(dtbydx(ny+js))
!      
!      u1(:)           = 0.0
!      f1(:)           = 0.0
!      dtbydx(:)       = 0.0
!      dtbydx(jmn:jmx) = dt * gy(3+js)%flt1d(jmn:jmx)
!
!      DO k = kmn,kmx
!       DO i = imn,imx      
!
!        DO j = 1,ny-1+js  
!
!         u1(j)     = 0.5*(v(i,j,k) +v(i-is,j-js,k-ks))
!         q1(j)     = s(i,j,k)
!         
!        ENDDO
!       
!       CALL ADVECT5_1D(q1,u1,f1,jmn,jmx,0,ny-1+js)
!       CALL ADVECT6_FILTER6_1D(q1,u1,f1,dtbydx,0.05,jmn,jmx,0,ny-1+js)
!       CALL ADVECT6_MONO_1D(q1,u1,f1,dtbydx,jmn,jmx,0,ny-1+js)
!        CALL ADVECT_WENO_1D(q1,u1,f1,jmn,jmx,0,ny-1+js)
!        
!        fy(i,jmn:jmx,k) = f1(jmn:jmx)
!       
!       ENDDO
!      ENDDO
!
!      deallocate(u1)
!      deallocate(f1)
!      deallocate(q1)
!      deallocate(dtbydx)
!         
!    ENDIF
!
!-----------------------------------------------------------------------------
! COMPUTE Z-INTERFACE FLUXES
!
!    allocate(u1(nz+ks))
!    allocate(f1(nz+ks))
!    allocate(q1(nz-1+ks))
!    allocate(dtbydx(nz+ks))
!    
!    u1(:)           = 0.0
!    f1(:)           = 0.0
!    dtbydx(:)       = 0.0
!    dtbydx(kmn:kmx) = dt * gz(3+ks)%flt1d(kmn:kmx)
!
!    DO i = imn,imx      
!     DO j = jmn,jmx      
!        
!      DO k = 1,nz-1+ks
!        
!       u1(k)     = 0.5*(w(i,j,k) +w(i-is,j-js,k-ks))
!       q1(k)     = s(i,j,k)
!    
!      ENDDO
!       
!     CALL ADVECT5_1D(q1,u1,f1,kmn,kmx,0,nz-1+ks)
!     CALL ADVECT6_FILTER6_1D(q1,u1,f1,dtbydx,0.05,kmn,kmx,0,nz-1+ks)
!     CALL ADVECT6_MONO_1D(q1,u1,f1,dtbydx,kmn,kmx,0,nz-1+ks)
!      CALL ADVECT_WENO_1D(q1,u1,f1,kmn,kmx,0,nz-1+ks)
!      
!      fz(i,j,kmn:kmx) = f1(kmn:kmx)
!       
!     ENDDO
!    ENDDO
!    
!    deallocate(u1)
!    deallocate(f1)
!    deallocate(q1)
!    deallocate(dtbydx)
!    
!    CALL cld_cpu('ADVECT-TEST')
!
!    RETURN
!    END SUBROUTINE ADVECT_TEST
!
 END SUBROUTINE ADVECT


#if defined(__ia64__)
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE ADVECT5TH(s,s0,fs,u,v,w,rrp,rrm,dt,nx,ny,nz,      &
                         imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1, &
                         imxb, jmxb, &
                         i1,i2,i3,j1,j2,j3,k1,k2,k3,is,js,ks,  &
                         gxt,gyt,gzt,xmask,ymask,varname,divmult)

!   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE PARAM_MODULE

#ifdef MPI
   USE COMMASMPI_MODULE
#endif

    implicit none

    integer imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1,imxb,jmxb,is,js,ks
    integer i1,i2,i3,j1,j2,j3,k1,k2,k3

   real    xmask(-ng+1:nx+ng), ymask(-ng+1:ny+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)

   integer, INTENT(IN)    :: nx, ny, nz
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: rrp(-ng+1:nz+ng), rrm(-ng+1:nz+ng)
   real,    INTENT(IN)    :: dt
   character(LEN=10)  :: varname
   double precision divmult

    integer         :: i, j, k, im1, jm1, km1, kz
    real            :: dir, vv
    real            :: ue,uw,vs,vn,wb,wt
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

   real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
   real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
   real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
   real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   integer kt,kb
   
   logical, parameter :: debug_mpi = .false.

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-5TH')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      fx2(:,:) = 0.0
      fy2(:,:) = 0.0
      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,ue,uw,im1,dir)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
         IF( (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2) THEN
#endif

          fx2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

         ELSE
         
          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vs,vn,vv,jm1,dir)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
         IF( (j .ge. 4 .and. j .le. j3) .or. bcy .eq. 2 ) THEN
#endif

          fy2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

         ELSE

          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,wt,wb,dir)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
#else
         IF( kz .ge. 4 .and. kz .le. k3 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult
       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-5TH')

    RETURN
    END ! SUBROUTINE ADVECT5

#endif
! \\ __ia64__ 


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


    SUBROUTINE ADVECT_WENOCF(s,s0,fs,u,v,w,rrp,rrm,dt,nx,ny,nz,      &
                         imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1, &
                         imxb, jmxb, &
                         i1,i2,i3,j1,j2,j3,k1,k2,k3,is,js,ks,  &
                         gxt,gyt,gzt,xmask,ymask,varname)

!   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE PARAM_MODULE
   USE MICRO_MODULE

#ifdef MPI
   USE COMMASMPI_MODULE
#endif

    implicit none

    integer imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1,imxb,jmxb,is,js,ks
    integer i1,i2,i3,j1,j2,j3,k1,k2,k3

   real    xmask(-ng+1:nx+ng), ymask(-ng+1:ny+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)

   integer, INTENT(IN)    :: nx, ny, nz
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: rrp(-ng+1:nz+ng), rrm(-ng+1:nz+ng)
   real,    INTENT(IN)    :: dt


   character(LEN=10)  :: varname

    integer         :: i, j, k, im1, jm1, km1, kz
    real            :: dir, vv
    real            :: ue,uw,vs,vn,wb,wt
    real, parameter :: f30 =  7./12., f31 = 1./12.
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.


   integer kt,kb
   
    
    real               :: qim2, qim1, qi, qip1, qip2
    double precision               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    double precision, parameter    :: gi0 = 1.d0/10.d0, gi1 = 6.d0/10.d0, gi2 = 3.d0/10.d0, eps=EPSVAL
    integer, parameter :: pw5 = 2

     real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
     real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
     real :: div2(-ng+1:nx+ng,-ng+1:ny+ng)

   
   logical, parameter :: debug_mpi = .false.
   double precision :: divmult ! = DIVM  ! should be set to 1.0d0 except for special tests

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    CALL cld_cpu('ADVECT-WENOCF')

   IF ( Abs( denscale ) > 1 .and. (varname(1:1) == 'C' .or. varname(1:1) == 'Q' &
                            .or. varname(1:1) == 'V'  .or. varname(1:6) == 'TRACER') ) THEN
      divmult = 0.0d0
   ELSE
      divmult = 1.0d0
   ENDIF

!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      fz2(:,:,:) = 0.0
      IF ( nx .le. 2 ) fx2(:,:) = 0.0
      IF ( ny .le. 2 ) fy2(:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      kzb = 0  ! need to start at zero to fill the kb level to use at k=1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1-kzbeg+1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

!      fx2(:,:) = 0.0
!      fy2(:,:) = 0.0
!      fz2(:,:,kt) = 0.0
      
    IF ( nx .gt. 2 ) THEN
       do j = jyb,jye
        do i = ixb,ixe
        im1 = i-1
        if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
      DO k = kmn,kmx1

      kb = 3 - kb
      kt = 3 - kt
      div2(:,:) = 0.0

    IF ( nx .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,ue,uw,im1,dir)
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         uw  = 0.5*(u(i  ,j,  k  ) + u(i-is  ,j-js,  k-ks  ))
         ue  = 0.5*(u(i+1,j,  k  ) + u(i-is+1,j-js,  k-ks  ))
         vv = uw

         div2(i,j) = gxt(i,3+is)*(ue-uw)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
         IF( (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2) THEN
#endif

!          fx2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i-1,j,k))  &
!                           - f51 * (s(i+1,j,k) + s(i-2,j,k))  &
!                           + f52 * (s(i+2,j,k) + s(i-3,j,k))  &
!                           - f52 * (s(i+2,j,k) - s(i-3,j,k)   &
!                           - 5.0 * (s(i+1,j,k) - s(i-2,j,k))  &
!                           + 10. * (s(i,  j,k) - s(i-1,j,k)))*dir )


         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i+1,j,k)
            qip1 = s(i,  j,k  )
            qi   = s(i-1,j,k)
            qim1 = s(i-2,j,k)
            qim2 = s(i-3,j,k)
          ELSE
            qip2 = s(i-2,j,k)
            qip1 = s(i-1,j,k)
            qi   = s(i,  j,k)
            qim1 = s(i+1,j,k)
            qim2 = s(i+2,j,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fx2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
          fx2(i,j) = vv * ( f30 * (s(i,j,k)   + s(i-1,j,k))  &
                           - f31 * (s(i+1,j,k) + s(i-2,j,k))  &
                           + f31 * (s(i+1,j,k) - s(i-2,j,k)   &
                           - 3.0 * (s(i,j,k)   - s(i-1,j,k)))*dir )

#ifdef MPI
!         ELSEIF( ixbeg-1+i .eq. 2 .or. ixbeg-1+i .eq. i1 ) THEN
#else
!         ELSEIF( i .eq. 2 .or. i .eq. i1 ) THEN
#endif
        ELSE
          fx2(i,j) = vv * 0.5 * (s(i,j,k) + s(im1,j,k))

!         ELSE
!         
!          fx2(i,j) = 0.0
         
         ENDIF
     
        ENDDO
       ENDDO

    ENDIF

!-----------------------------------------------------------------------------
! Y-INTERFACE

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
    IF ( ny .gt. 2 ) THEN
      do j = jyb,jye
       jm1 = j-1
       if(jybeg.eq.nybeg) jm1 = max(j-1,1)
      do i = ixb,ixe
#else

    IF ( ny .gt. 2 .and. k .ge. kmn ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vs,vn,vv,jm1,dir)
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))

         vs  = vv
         vn  = 0.5*(v(i  ,j+1,k  ) + v(i-is  ,j-js+1,k-ks  ))

         div2(i,j) = div2(i,j) + gyt(j,3+js)*(vn-vs)

         dir = sign(1.0,vv)

#ifdef MPI
         IF( (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
         IF( (j .ge. 4 .and. j .le. j3) .or. bcy .eq. 2 ) THEN
#endif

!          fy2(i,j) = vv * ( f50 * (s(i,  j,k) + s(i,j-1,k))  &
!                           - f51 * (s(i,j+1,k) + s(i,j-2,k))  &
!                           + f52 * (s(i,j+2,k) + s(i,j-3,k))  &
!                           - f52 * (s(i,j+2,k) - s(i,j-3,k)   &
!                           - 5.0 * (s(i,j+1,k) - s(i,j-2,k))  &
!                           + 10. * (s(i,j,  k) - s(i,j-1,k)))*dir )


         IF ( dir .ge. 0.0 ) THEN
            qip2 = s(i,j+1,k)
            qip1 = s(i,j,  k)
            qi   = s(i,j-1,k)
            qim1 = s(i,j-2,k)
            qim2 = s(i,j-3,k)
          ELSE
            qip2 = s(i,j-2,k)
            qip1 = s(i,j-1,k)
            qi   = s(i,j,  k)
            qim1 = s(i,j+1,k)
            qim2 = s(i,j+2,k)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fy2(i,j) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk


#ifdef MPI
         ELSEIF( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
          fy2(i,j) = vv * ( f30 * (s(i,j,  k) + s(i,j-1,k))  &
                           - f31 * (s(i,j+1,k) + s(i,j-2,k))  &
                           + f31 * (s(i,j+1,k) - s(i,j-2,k)   &
                           - 3.0 * (s(i,j,  k) - s(i,j-1,k)))*dir )

#ifdef MPI
!         ELSEIF ( jybeg-1+j .eq. 2 .or. jybeg-1+j .eq. j1 ) THEN
#else
!         ELSEIF ( j .eq. 2 .or. j .eq. j1 ) THEN
#endif
        ELSE
          fy2(i,j) = vv * 0.5 * (s(i,j,k) + s(i,jm1,k))

!         ELSE

!          fy2(i,j) = 0.0

         ENDIF

        ENDDO
       ENDDO

     ENDIF

! Z-INTERFACE

     IF ( nz .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = Max( ixb, imn-ixbeg+1 )
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

       kz =  k + 1
       km1 = kz-1
       if(kzbeg.eq.nzbeg) km1 = max(kz-1,1)
      do j = jyb,jye
      do i = ixb,ixe
#else
       kz = k + 1
       km1 = max(kz-1,1)

       
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,vv,wt,wb,dir)
       DO j = jmn,jmx
        DO i = imn,imx
#endif
         vv  = 0.5*(w(i,j,kz) + w(i-is,j-js,kz-ks))
         dir = sign(1.0,vv)

        wt  = vv*rrp(kz-1)

        wb        = 0.5*(w(i  ,j,  k  ) + w(i-is  ,j-js,  k-ks  ))*rrm(k)

         div2(i,j) = div2(i,j) + gzt(k,3+ks)*(wt-wb)

#ifdef MPI
         IF( kzbeg-1+kz .ge. 4 .and. kzbeg-1+kz .le. k3 ) THEN
#else
         IF( kz .ge. 4 .and. kz .le. k3 ) THEN
#endif
!          fz2(i,j,kt) = vv * ( f50 * (s(i,j,kz  ) + s(i,j,kz-1))  &
!                           - f51 * (s(i,j,kz+1) + s(i,j,kz-2))  &
!                           + f52 * (s(i,j,kz+2) + s(i,j,kz-3))  &
!                           - f52 * (s(i,j,kz+2) - s(i,j,kz-3)   &
!                           - 5.0 * (s(i,j,kz+1) - s(i,j,kz-2))  &
!                           + 10. * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )


         IF( dir .ge. 0.0 ) THEN
            qip2 = s(i,j,kz+1)
            qip1 = s(i,j,kz  )
            qi   = s(i,j,kz-1)
            qim1 = s(i,j,kz-2)
            qim2 = s(i,j,kz-3)
          ELSE
            qip2 = s(i,j,kz-2)
            qip1 = s(i,j,kz-1)
            qi   = s(i,j,kz  )
            qim1 = s(i,j,kz+1)
            qim2 = s(i,j,kz+2)
         ENDIF
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw5
         wi1 = gi1 / (eps + beta1)**pw5
         wi2 = gi2 / (eps + beta2)**pw5
    
         sumwk = wi0 + wi1 + wi2
    
         fz2(i,j,kt) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk


#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 3 .or. kzbeg-1+kz .eq. k2 ) THEN
#else
         ELSEIF( kz .eq. 3 .or. kz .eq. k2 ) THEN
#endif
          fz2(i,j,kt) = vv * ( f30 * (s(i,j,kz  ) + s(i,j,kz-1))  &
                           - f31 * (s(i,j,kz+1) + s(i,j,kz-2))  &
                           + f31 * (s(i,j,kz+1) - s(i,j,kz-2)   &
                           - 3.0 * (s(i,j,kz  ) - s(i,j,kz-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+kz .eq. 2 .or. kzbeg-1+kz .eq. k1 ) THEN
#else
         ELSEIF( kz .eq. 2 .or. kz .eq. k1 ) THEN
#endif
          fz2(i,j,kt) = vv * 0.5 * (s(i,j,kz) + s(i,j,km1))

         ELSE
          
          fz2(i,j,kt) = 0.0
          
         ENDIF
     
        ENDDO
       ENDDO

       ENDIF ! ( nz .gt. 2 )




#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = Max(ixb,imn -ixbeg+1)
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1



   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j)
   DO j = jmn,jmx1 
    DO i = imn,imx1 
#endif
      IF ( k .ge. kmn-ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j) * divmult
       ENDIF
      ENDDO
      ENDDO
      
      
      ENDDO
      
!    ENDIF

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    CALL cld_cpu('ADVECT-WENOCF')

    RETURN
    END ! SUBROUTINE ADVECT-WENOCF

!
!--------------------------------------------------------------------------
! SUBROUTINE ADVECT5_1D:  5th-order upwind-biased advection for 1D
! 
!  q(:) is dimensioned the ACTUAL size of the computational grid
!  u(:) is dimensioned size(q) + 1
! fq(:) is dimensioned size(q) + 1, the flux is returned
!
!  SUBROUTINE ADVECT5_1D(q,u,fq,imin,imax,bc,ni)
!  
!    implicit none
!    
!    real,    intent(in)  ::  q(ni) 
!    real,    intent(in)  ::  u(ni+1)
!    real,    intent(out) :: fq(ni+1)
!    integer, intent(in)  :: imin, imax, bc, ni
!    
!    integer         :: i
!    real            :: dir
!    real, parameter :: f30 =  7./12., f31 = 1./12
!    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.
!
!-----------------------------------------------------------------------------
! COMPUTE INTERFACE FLUXES
!
!   DO i = imin,imax 
!
!     dir = sign(1.0,u(i))
!
!     IF( i .ge. 4 .and. i .le. ni-2 ) THEN
!
!       fq(i) = u(i) * ( f50 * (q(i)   + q(i-1))  &
!                      - f51 * (q(i+1) + q(i-2))  &
!                      + f52 * (q(i+2) + q(i-3))  &
!                      - f52 * (q(i+2) - q(i-3) - 5.*(q(i+1)-q(i-2)) + 10.*(q(i)-q(i-1)))*dir )
!
!     ELSEIF( i .eq. 3 .or. i .eq. ni-1 ) THEN
!
!       fq(i) = u(i) * ( f30 * (q(i)   + q(i-1))  &
!                      - f31 * (q(i+1) + q(i-2))  &
!                      + f31 * (q(i+1) - q(i-2) - 3.*(q(i)-q(i-1)))*dir )
!
!     ELSE
!
!       fq(i) = u(i) * 0.5 * (q(i) + q(i-1))
!
!     ENDIF
!        
!   ENDDO
!
!  END SUBROUTINE ADVECT5_1D
!  
!---------------------------------------------------------------------------------------------
! SUBROUTINE ADVECT6_FILTER6_1D:  6th order advection plus (Del)^6 filter with Xue monotonic limiter
! 
!  q(:) is dimensioned the ACTUAL size of the computational grid
!  u(:) is dimensioned size(q) + 1
! fq(:) is dimensioned size(q) + 1, the flux is returned
!
! fdamp is % of 2 dx filter per time step
! fdamp == 0 is monotonic filter with same coefficient as 5th order scheme [|u|*dt/(60*dx)]
! fdamp != 0 is max { [|u|*dt/(60*dx)], % of 2 dx filter coeff }
!
!  SUBROUTINE ADVECT6_FILTER6_1D(q,u,fq,dtbydx,fdamp,imin,imax,bc,ni)
!  
!    implicit none
!    
!    real,    intent(in)  ::  q(ni) 
!    real,    intent(in)  ::  u(ni+1)
!    real,    intent(in)  ::  dtbydx(ni+1)
!    real,    intent(out) :: fq(ni+1)
!    integer, intent(in)  :: imin, imax, bc, ni
!    real,    intent(in)  :: fdamp
!    
!    integer         :: i
!    real, parameter :: f40 =  7./12., f41 = 1./12
!    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
!    real            :: qi, qid, damp
!
! COMPUTE INTERFACE FLUXES
!
!   DO i = imin,imax 
!
!     IF( i .ge. 4 .and. i .le. ni-2 ) THEN
!
!       qi  = f60 * (q(i)   + q(i-1))  &
!           - f61 * (q(i+1) + q(i-2))  &
!           + f62 * (q(i+2) + q(i-3))  
!                      
!       qid = (-q(i-3)+5.*q(i-2)-10.*q(i-1)+10.*q(i)-5.*q(i+1)+q(i+2))
!       qid = qid*amax1(0.0, sign(1.0,qid*(q(i)-q(i-1)) ) )
!       
!       damp   = amax1( abs(u(i)) / 60., fdamp / (64.*dtbydx(i)) )
!       fq(i)  = u(i) * qi - damp * qid
!
!     ELSEIF( i .eq. 3 .or. i .eq. ni-1 ) THEN
!
!       qi  = f40 * (q(i)   + q(i-1))  &
!           - f41 * (q(i+1) + q(i-2))
!                      
!       qid = -(q(i+1) - q(i-2) - 3.*(q(i)-q(i-1)) )
!       qid = qid*amax1(0.0, sign(1.0,qid*(q(i)-q(i-1)) ) )
!       
!       damp   = amax1( abs(u(i)) / 12., fdamp / (16.*dtbydx(i)) )
!       fq(i)  = u(i) * qi - damp * qid
!
!     ELSE
!
!       qi  = 0.5 * (q(i) + q(i-1))
!       
!       qid = (q(i)-q(i-1))/4.0
!       qid = qid*amax1(0.0, sign(1.0,qid*(q(i)-q(i-1)) ) )
!       
!       damp   = amax1( abs(u(i)) / 2., fdamp / (4.*dtbydx(i)) )
!       fq(i)  = u(i) * qi - damp * qid
!       
!     ENDIF
!
!   ENDDO
!
!  END SUBROUTINE ADVECT6_FILTER6_1D
!--------------------------------------------------------------------------
! SUBROUTINE ADVECT_WENO_1D:  5th-order WENO advection for 1D
! 
!  q(:) is dimensioned the ACTUAL size of the computational grid
!  u(:) is dimensioned size(q) + 1
! fq(:) is dimensioned size(q) + 1, the flux is returned
!
!  SUBROUTINE ADVECT_WENO_1D(q,u,fq,imin,imax,bc,ni)
!  
!    implicit none
!    
!    real,    intent(in)  ::  q(ni) 
!    real,    intent(in)  ::  u(ni+1)
!    real,    intent(out) :: fq(ni+1)
!    integer, intent(in)  :: imin, imax, bc, ni
!    
!    integer            :: i, im1
!    real               :: dir
!    real, parameter    :: f30 =  7./12., f31 = 1./12
!    real               :: qim2, qim1, qi, qip1, qip2
!    real               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
!    real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-15
!    integer, parameter :: pw = 1
!
!-----------------------------------------------------------------------------
! COMPUTE INTERFACE FLUXES
!
!   DO i = imin,imax 
!
!     im1 = max(i-1,1)
!     dir = sign(1.0,u(i))
!
!     IF( i .ge. 4 .and. i .le. ni-2 ) THEN
!
!       IF( dir .ge. 0.0 ) THEN
!          qip2 = q(i+1)
!          qip1 = q(i  )
!          qi   = q(i-1)
!          qim1 = q(i-2)
!          qim2 = q(i-3)
!        ELSE
!          qip2 = q(i-2)
!          qip1 = q(i-1)
!          qi   = q(i  )
!          qim1 = q(i+1)
!          qim2 = q(i+2)
!       ENDIF
!    
!       f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
!       f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
!       f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
!    
!       beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
!       beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
!       beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
!    
!       wi0 = gi0 / (eps + beta0)**pw
!       wi1 = gi1 / (eps + beta1)**pw
!       wi2 = gi2 / (eps + beta2)**pw
!    
!       sumwk = wi0 + wi1 + wi2
!    
!       fq(i) = u(i) * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk
!
!     ELSEIF( i .eq. 3 .or. i .eq. ni-1 ) THEN
!
!       fq(i) = u(i) * ( f30 * (q(i)   + q(i-1))  &
!                      - f31 * (q(i+1) + q(i-2))  &
!                      + f31 * (q(i+1) - q(i-2) - 3.*(q(i)-q(i-1)))*dir )
!
!     ELSE
!
!       fq(i) = u(i) * 0.5 * (q(i) + q(im1))
!
!     ENDIF
!        
!   ENDDO
!
!  END SUBROUTINE ADVECT_WENO_1D
!
!--------------------------------------------------------------------------
! SUBROUTINE ADVECT6_MONO_1D:  6th-order advection with Leonard monotonicity
! 
!  q(:) is dimensioned the ACTUAL size of the computational grid
!  u(:) is dimensioned size(q) + 1
! fq(:) is dimensioned size(q) + 1, the flux is returned
!
!  SUBROUTINE ADVECT6_MONO_1D(q,u,fq,dtbydx,imin,imax,bc,ni)
!  
!    implicit none
!    
!    real,    intent(in)  ::  q(ni) 
!    real,    intent(in)  ::  u(ni+1)
!    real,    intent(in)  ::  dtbydx(ni+1)
!    real,    intent(out) :: fq(ni+1)
!    integer, intent(in)  :: imin, imax, bc, ni
!    
!    integer         :: nx, i, idir, im1, im2, ip1
!    real, parameter :: f40 =  7./12., f41 = 1./12
!    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
!    real            :: qi, qoutmin, qoutmax, qinmax, qinmin, qD, qC, qU, cr
!    real, parameter :: eps = 1.0e-15 
!
!-----------------------------------------------------------------------------
! COMPUTE INTERFACE FLUXES  
!  
!    DO i = imin,imax 
!
!     im1 = max(i-1,1)
!     im2 = max(i-2,1)
!     ip1 = min(i+1,nx)
!
!     idir = 1
!     IF( u(i) .lt. 0.0 ) idir = -1
!
!     qD =  0.5*((1 + idir)*q(i)   - (idir - 1)*q(im1))
!     qC =  0.5*((1 + idir)*q(im1) - (idir - 1)*q(i  ))
!     qU =  0.5*((1 + idir)*q(im2) - (idir - 1)*q(ip1))
!
!     IF( i .ge. 4 .and. i .le. ni-2 ) THEN
!
!       qi  = f60 * (q(i)   + q(i-1))  &
!           - f61 * (q(i+1) + q(i-2))  &
!           + f62 * (q(i+2) + q(i-3))  
!                      
!     ELSEIF( i .eq. 3 .or. i .eq. ni-1 ) THEN
!
!       qi  = f40 * (q(i)   + q(i-1))  &
!           - f41 * (q(i+1) + q(i-2))
!                      
!     ELSE
!
!       qi  = 0.5 * (q(i) + q(i-1))
!       
!     ENDIF
!
!     qinmin  = min(qD, qC)
!     qinmax  = max(qD, qC)
!     qi      = max(qi, qinmin)
!     qi      = min(qi, qinmax)
!
!     cr      = min(dtbydx(i) * abs(u(i)), 1.0)
!     qoutmin = (qC - (1. - cr) * amax1(qC, qU)) / (cr + eps)
!     qoutmax = (qC - (1. - cr) * amin1(qC, qU)) / (cr + eps)
!
!     qi      = amax1(qi, qoutmin)
!     fq(i)   = u(i) * amin1(qi, qoutmax)
!
!    ENDDO
!
!  END SUBROUTINE ADVECT6_MONO_1D
