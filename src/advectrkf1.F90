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
!       IF( ATYPE .eq. 1 ) CALL ADVECT_WENO()
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
  SUBROUTINE ADVECT(s,s0,fs,u,v,w,           &
                    gxt,gyt,gzt,rrp,rrm,dt,  &
                    nx,ny,nz,                &
                    svar,atype,izero0)

   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE PARAM_MODULE

#ifdef MPI
   USE COMMASMPI_MODULE
#endif

   implicit none 

   integer, INTENT(IN)    :: nx, ny, nz, atype, izero0
   real,    INTENT(INOUT) :: s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: s0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(INOUT) :: fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real,    INTENT(IN)    :: rrp(-ng+1:nz+ng), rrm(-ng+1:nz+ng)
   real,    INTENT(IN)    :: dt

!   TYPE(VARIABLE) :: gx(4) , gy(4), gz(4)
   TYPE(VARIABLE) :: svar

   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)

! Local vars

   integer is, js, ks
   integer n, i, j, k
   integer izero
   real    uw, ue, vn, vs, wt, wb, unorm, vnorm
   real    qx1, ux1, qz0, qz1, wz1, qx0, vy0, vy1, qy0, qy1
   real,  allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
   real    div
   real    sijk
   real    sim3, sim2, sim1, sip1, sip2
   real    sjm3, sjm2, sjm1, sjp1, sjp2
   real    skm3, skm2, skm1, skp1, skp2
   real    xmask(-ng+1:nx+ng), ymask(-ng+1:ny+ng), zmask(-ng+1:nz+ng)

   integer imn,imx,jmn,jmx,kmn,kmx,imx1,jmx1,kmx1
   integer imxb, jmxb
   integer i1,i2,i3
   integer j1,j2,j3
   integer k1,k2,k3
         
   real damp1
   real dt1
   character(LEN=10)  :: varname

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   
   logical :: debug_mpi = .false.
#endif

! Parameter flags

! Relaxation at the boundaries parameters
   
   real, parameter    :: ainflo = 0.0                 !  Leave as 0.0, test value is 0.5
   real               :: ainflow(nz), ainfloe(nz), ainflos(nz), ainflon(nz)
   integer            :: kbl = 7
   
! Diffusion (only applied in the horizontal) DAMPH IS IN PARAM_MODULE!!

!-----------------------------------------------------------------------------
! Determine staggering
      
   is      = svar%istag
   js      = svar%jstag
   ks      = svar%kstag

   varname = svar%name

!-----------------------------------------------------------------------------

#ifdef MPI
   if (debug_mpi)  write(0,"('ADVECTRK: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

  IF ( ny .le. 2 .and. svar%name .eq. 'V' ) RETURN
  IF ( nx .le. 2 .and. svar%name .eq. 'U' ) RETURN

#ifdef MPI
  if ( bcx .eq. 2 .or. bcy .eq. 2) then
    write(luno,*) 'ADVECTRK: MPI not tested for periodic bc yet'
    write(luno,*) 'ADVECTRK: stopping run ...'
    STOP
  endif
#endif

!-----------------------------------------------------------------------------
! Set bounds for advection if izero=1 to avoid advecting zeroes.

  izero = izero0

#ifdef MPI
  izero = 0
#endif

  IF ( izero .ge. 1 ) THEN
    CALL cld_cpu('ADVECT-IZERO1')
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

    imn = MAX(imn - 3, 1)
    imx = MIN(imx + 4,nx)
    jmn = MAX(jmn - 3, 1)
    jmx = MIN(jmx + 4,ny)
    kmn = MAX(kmn - 3, 1)
    kmx = MIN(kmx + 4,nz)

    CALL cld_cpu('ADVECT-IZERO1')

    IF ( imn .gt. imx ) RETURN

    imx1 = Min(imx,nx-1)
    jmx1 = Min(jmx,ny-1)
    kmx1 = Min(kmx,nz-1)

! Check if upper bound is below nz-1; if so, need to reduce kmx1 so that
! the flux at kmx1+1 is defined (i.e., at kmx)

    IF ( kmx .lt. nz-1 ) kmx1 = kmx-1

  ELSE ! izero

#ifdef MPI
    imn = 1+is
    imx = nxend-1+is
    jmn = 1+js
    jmx = nyend-1+js
!    kmn = 1+ks
    kmn = 1
    kmx = nzend-1+ks
    imx1 = nxend-1
    jmx1 = nyend-1
    kmx1 = nzend-1
#else
    imn = 1+is
    imx = nx-1+is
    jmn = 1+js
    jmx = ny-1+js
!    kmn = 1+ks
    kmn = 1
    kmx = nz-1+ks
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
    i3 = nxend-3+is
    i2 = nxend-2+is
    i1 = nxend-1+is

    j3 = nyend-3+js
    j2 = nyend-2+js
    j1 = nyend-1+js

    k3 = nzend-3+ks
    k2 = nzend-2+ks
    k1 = nzend-1+ks
    
    imxb = imx
    jmxb = jmx
#else
    i3 = nx-3+is
    i2 = nx-2+is
    i1 = nx-1+is

    j3 = ny-3+js
    j2 = ny-2+js
    j1 = ny-1+js

    k3 = nz-3+ks
    k2 = nz-2+ks
    k1 = nz-1+ks

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

  IF ( bcx .eq. 1 ) THEN   ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug: no need to fill ghost zones
   DO n = 1,ng
    s( 1-n      ,1:ny,1:nz) = s(1      ,1:ny,1:nz)
    s( nx-1+n+is,1:ny,1:nz) = s(nx-1+is,1:ny,1:nz)
   ENDDO
#endif

  ELSEIF( bcx .eq. 2 ) THEN    ! periodic

   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg+1

    do k = kzb,kze ; do j = jyb, jye
     s( nxend-1+n+is,j,k) = s(n+is   ,j,k)
     s( 1-n         ,j,k) = s(nxend-n,j,k)
    enddo ; enddo
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

  IF ( bcy .eq. 1 ) THEN    ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug:  no need to fill ghost zones
   DO n = 1,ng
    s(1:nx, 1-n     ,1:nz) = s(1:nx,1      ,1:nz)
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,ny-1+js,1:nz) 
   ENDDO
#endif

  ELSEIF( bcy .eq. 2 ) THEN    ! periodic

   DO n = 1,ng
#ifdef MPI
    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg+1

    ixb = 1
    ixe = jtile
    if (ixend .eq. nxend) ixe = ixend-ixbeg+1
!mpi need processor communication here
    do k = kzb,kze ; do i = ixb,ixe
     s(i,nyend-1+n+js,k) = s(i,n+js,k) 
     s(i,1-n      ,k) = s(i,nyend-n,k)
    enddo ; enddo
#else
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,n+js,1:nz) 
    s(1:nx,1-n      ,1:nz) = s(1:nx,ny-n,1:nz)
#endif
   ENDDO

    IF ( jmn .eq. 1 + js ) THEN
      jmn = 1
      jmx = ny-1+js
    ENDIF
    IF ( jmx .eq. ny-1+js ) THEN 
      jmxb = ny + js
      jmn = 1
    ENDIF
    
  ELSEIF ( bcy .eq. 0 ) THEN    ! mirror (solid wall) on north boundary, open on south:
   DO n = 1,ng
!    s(1:nx, 1-n     ,1:nz) = s(1:nx,n      ,1:nz)
    s(1:nx, 1-n     ,1:nz) = s(1:nx,1      ,1:nz)
    s(1:nx,ny-1+n+js,1:nz) = s(1:nx,ny-n+js,1:nz) 
   ENDDO

  ENDIF !bcy

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
!  IF( bcz .ge. 0 .and. ks .eq. 0 ) THEN    ! reflective
!
!   DO n = 1,ng
!   s(1:nx,1:ny,   1-n   )      = s0(1:nx,1:ny,n+ks)
!   s(1:nx,1:ny,nz-1+n+ks)      = s0(1:nx,1:ny,nz-n)
!    s(1:nx,1:ny,   1-n   )      = s(1:nx,1:ny,n+ks)
!    s(1:nx,1:ny,nz-1+n+ks)      = s(1:nx,1:ny,nz-n)
!   ENDDO
!
!  ENDIF
!
!  IF( bcz .ge. 0 .and. ks .eq. 1 ) THEN    ! anti - reflective
!
!   DO n = 1,ng
!    s (1:nx,1:ny,   1-n   )      = -s0(1:nx,1:ny,n+ks)
!    s (1:nx,1:ny,nz-1+n+ks)      = -s0(1:nx,1:ny,nz-n)
!   ENDDO
!
!  ENDIF

  CALL cld_cpu('ADVECT-SETBC')

!--------------------------------------------------------------------------
! DO THE VARIOUS TYPES OF ADVECTION


    CALL cld_cpu('ADVECT-IZERO2')

    IF ( ATYPE .ne. 0 ) THEN

#ifdef MPI
  allocate (  fx(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  allocate (  fy(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
  allocate (  fz(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )


!   fx(:,:,:) = 0.0
!   fy(:,:,:) = 0.0
!   fz(:,:,:) = 0.0

#ifdef __ia64__
  fx(-ng+1:1,:      ,:      ) = 0.0 
  fy(:      ,-ng+1:1,:      ) = 0.0 
  fz(:      ,:      ,-ng+1:1) = 0.0 

  fx(nx:nx+ng,:       ,:       ) = 0.0 
  fy(:       ,ny:ny+ng,:       ) = 0.0 
  fz(:       ,:       ,nz:nz+ng) = 0.0 
#endif

#else
  allocate (  fx(0:nx+1,ny,nz) )
  allocate (  fy(nx,0:ny+1,nz) )
  allocate (  fz(nx,ny,nz) )

#ifdef __ia64__
  fx(0:1,:  ,:) = 0.0 
  fy(:  ,0:1,:) = 0.0 
  fz(:  ,:  ,1) = 0.0 

  fx(nx,: ,: ) = 0.0 
  fy(: ,ny,: ) = 0.0 
  fz(: ,: ,nz) = 0.0 
#endif

#endif

   ENDIF

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
    imx1 = imx1 + is
!    fx(nx+is,1:ny-1+js,1:nz-1+ks) = fx(1,1:ny-1+js,1:nz-1+ks)
  ENDIF

  ymask(:) = 1.0 
  IF ( bcy .ne. 2 ) THEN
#ifdef MPI
    if (jybeg .eq. nybeg)  ymask(-ng+1:1) = 0.0
    if (jyend .eq. nyend) then
      jyb = jmx-jybeg+1
      jye = jyend-jybeg+1+ng

      ymask(jyb:jye) = 0.0
    endif
#else
    ymask(-ng+1:1) = 0.0 ; ymask(jmx:ny+ng) = 0.0
#endif
  ENDIF
  

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

!  write(0,*) 'Advect ',varname
  IF( ATYPE .eq.-1 ) CALL ADVECT5N()
!  IF( ATYPE .eq. 0 ) CALL ADVECT5()
  IF( ATYPE .eq. 0 ) CALL ADVECT5TH(s,s0,fs,u,v,w,rrp,rrm,dt,nx,ny,nz,      &
                         imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1, &
                         imxb, jmxb, &
                         i1,i2,i3,j1,j2,j3,k1,k2,k3,is,js,ks,  &
                         gxt,gyt,gzt,xmask,ymask,varname)

  IF( ATYPE .eq. 1 ) CALL ADVECT_WENO()
  IF( ATYPE .eq. 2 ) CALL ADVECT6()
  IF( ATYPE .eq. 3 ) CALL ADVECT_MONO()
  IF( ATYPE .eq. 4 ) CALL ADVECT_MONOF(damp1,0.0,.false.)
  IF( ATYPE .eq. 5 ) CALL ADVECT_MONOF1(damp1)
  IF( ATYPE .eq. 6 ) CALL ADVECT_MONOF(0.0,1.0,.true.)
  IF( ATYPE .eq. 7 ) CALL ADVECT_MONOF(damp1,1.0,.false.)
  IF( ATYPE .eq. 8 ) CALL ADVECT_MONOF(damp1,1.0,.true.)

!-----------------------------------------------------------------------------
! Set boundary masks

! xmask(:) = 1.0 ; xmask(1) = 0.0 ; xmask(nx-1+is:nx) = 0.0
! ymask(:) = 1.0 ; ymask(1) = 0.0 ; ymask(ny-1+js:ny) = 0.0


  IF ( ATYPE .ne. 0 ) THEN
  
   CALL cld_cpu('ADVECT-FLUX')

#ifdef MPI
  ixb = 1
  ixe = itile
  if (ixbeg .le. imn)  ixb = imn -ixbeg+1
  if (ixend .ge. imx1) ixe = imx1-ixbeg+1

  jyb = 1
  jye = jtile
  if (jybeg .le. jmn)  jyb = jmn -jybeg+1
  if (jyend .ge. jmx1) jye = jmx1-jybeg+1

  kzb = 1
  kze = ktile
  if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
  if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

  do k = kzb,kze
   do j = jyb,jye
    do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,ue,vs,vn,wb,wt,div)
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

!      fs(i,j,k) =  - (zmask(k+1)*fz(i  ,j  ,k+1)*rrp(k) - zmask(k)*fz(i,j,k)*rrm(k))*gzt(k,3+ks) &
      fs(i,j,k) =  -            (fz(i  ,j  ,k+1)*rrp(k) - fz(i,j,k)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy(i  ,j+1,k  )        - fy(i,j,k)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx(i+1,j  ,k  )        - fx(i,j,k)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div

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
  
  IF ( nx .gt. 2 .and. bcx .ne. 2 ) THEN
#ifdef MPI
   IF( imn .eq. 1 .and. ixbeg .eq. nxbeg) THEN

    jyb = 1
    jye = jtile
    if (jybeg .le. jmn) jyb = jmn
    if (jyend .ge. jmx1) jye = jmx1-jybeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn) kzb = kmn
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

    do k = kzb,kze ; do j = jyb,jye

#else
   IF( imn .eq. 1 ) THEN

     DO k = kmn,kmx1 
      DO j = jmn,jmx1 
#endif
       uw = 0.5*(u(1,j,k) + u(1,j-js,k-ks))
       ue = 0.5*(u(2,j,k) + u(2,j-js,k-ks))
       unorm = uw ! 0.5*(ue + uw)
       fs(1,j,k)   = fs(1,j,k)  - (amin1(unorm,0.0)*(s0(2,j,k)-s0(1,j,k))     &
                                   +  s0(1,j,k)*(ue-uw)       &
                                   +  (ainflow(k)*Max(unorm, 0.0))*(s0(1,j,k) - svar%base1d(k))) * gxt(1,3)
      ENDDO
     ENDDO
   ENDIF

#ifdef MPI
   IF( imx .eq. nxend-1 .and. ixend .eq. nxend                                &
                        .and. bcx .ne. 2 ) THEN

    jyb = 1
    jye = jtile
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmx1) jye = jmx1-jybeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

     do k = kzb,kze
      do j = jyb,jye
#else
   IF( imx .eq. nx-1 .and. bcx .ne. 2 ) THEN

     DO k = kmn,kmx1 
      DO j = jmn,jmx1 
#endif
       uw = 0.5*(u(nx-1,j,k) + u(nx-1,j-js,k-ks))
       ue = 0.5*(u(nx,  j,k) + u(nx,  j-js,k-ks))
       unorm = ue ! 0.5*(ue + uw)
       fs(nx-1,j,k) = fs(nx-1,j,k) - (amax1(unorm,0.0)*(s0(nx-1,j,k)-s0(nx-2,j,k)) &
                                      +  s0(nx-1,j,k)*(ue-uw)       &
                                      -  (ainfloe(k)*Min(unorm, 0.0)) *(s0(nx-1,j,k) - svar%base1d(k))) * gxt(nx-1,3)
      ENDDO
     ENDDO

   ENDIF
   
  ENDIF

  IF ( ny .gt. 2 .and. bcy .ne. 2 ) THEN
#ifdef MPI
   IF( jmn .eq. 1 .and. jybeg .eq. nybeg) THEN

    ixb = 1
    ixe = itile
    if (ixbeg .le. imn)  ixb = imn -ixbeg+1
    if (ixend .ge. imx1) ixe = imx1-ixbeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

    do k = kzb,kze ; do i = ixb,ixe

#else
   IF( jmn .eq. 1 ) THEN

     DO k = kmn,kmx1 
      DO i = imn,imx1 
#endif
       vs = 0.5*(v(i,1,k) + v(i-is,1,k-ks))
       vn = 0.5*(v(i,2,k) + v(i-is,2,k-ks))
       unorm = vs ! 0.5*(vs + vn)
       fs(i,1,k) = fs(i,1,k) - (amin1(unorm,0.0)*(s0(i,2,k)-s0(i,1,k)) &
                                 + s0(i,1,k)*(vn-vs)       &
                                 +  (ainflos(k)*Max(unorm, 0.0))*(s0(i,1,k) - svar%base1d(k))) * gyt(1,4)
      ENDDO
     ENDDO
   ENDIF

#ifdef MPI
   IF( jmx .eq. nyend-1 .and. jyend .eq. nyend                                &
                        .and. bcy .ne. 2) THEN

    ixb = 1
    ixe = itile
    if (ixbeg .le. imn)  ixb = imn -ixbeg+1
    if (ixend .ge. imx1) ixe = imx1-ixbeg+1

    kzb = 1
    kze = ktile
    if (kzbeg .le. kmn)  kzb = kmn -kzbeg+1
    if (kzend .ge. kmx1) kze = kmx1-kzbeg+1

    do k = kzb,kze ; do i = ixb,ixe

#else
   IF( jmx .eq. ny-1 .and. bcy .ne. 2) THEN

     DO k = kmn,kmx1 
      DO i = imn,imx1 
#endif
       vs = 0.5*(v(i,ny-1,k) + v(i-is,ny-1,k-ks))
       vn = 0.5*(v(i,ny,  k) + v(i-is,ny,  k-ks))
       unorm = vn ! 0.5*(vs + vn)
       fs(i,ny-1,k) = fs(i,ny-1,k) - (amax1(unorm,0.0)*(s0(i,ny-1,k)-s0(i,ny-2,k)) &
                                   + s0(i,ny-1,k)*(vn-vs)       &
                                   -  (ainflon(k)*Min(unorm, 0.0)) *(s0(i,ny-1,k) - svar%base1d(k))) * gyt(ny-1,4)
      ENDDO
     ENDDO
   ENDIF
  ENDIF

   CALL cld_cpu('ADVECT-INFLOW')



 IF ( allocated( fx ) ) THEN
  deallocate ( fx )
  deallocate ( fy )
  deallocate ( fz )
 ENDIF


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECTRK: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

  RETURN

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM DEFINITIONS (Fortran90 stuff....)

  CONTAINS

!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE ADVECT5()

    implicit none

    integer         :: i, im1, jm1, km1, kz
    real            :: dir, vv
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

   real :: fx2(-ng+1:nx+ng,-ng+1:ny+ng)
   real :: fy2(-ng+1:nx+ng,-ng+1:ny+ng)
   real :: fz2(-ng+1:nx+ng,-ng+1:ny+ng,2)
!   real :: div2x(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2y(-ng+1:nx+ng,-ng+1:ny+ng)
!   real :: div2z(-ng+1:nx+ng,-ng+1:ny+ng,2)
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

    CALL cld_cpu('ADVECT-5TH')


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES

      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = imn -ixbeg+1
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

    IF ( nx .gt. 2 .and. k .ge. kmn+ks ) THEN
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
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
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

    IF ( ny .gt. 2 .and. k .ge. kmn+ks ) THEN
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
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
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
  if (ixbeg .le. imn)  ixb = imn -ixbeg+1
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
      IF ( k .ge. kmn+ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j)
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
    END SUBROUTINE ADVECT5



!--------------------------------------------------------------------------
! INTERNAL SUBPROGRAM ADVECTN:  5th-order upwind-biased advection with No stencil collapse

    SUBROUTINE ADVECT5N()
    
    integer         :: i, im1, jm1, km1
    real            :: dirx,diry,dirz, uu, vv, ww
    real, parameter :: f30 =  7./12., f31 = 1./12
    real, parameter :: f50 = 37./60., f51 = 2./15., f52 = 1./60.

    CALL cld_cpu('ADVECT-5N')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT5N: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
   write(0,*) "ADVECT5N: not updated for MPI yet"
#endif


!-----------------------------------------------------------------------------
! COMPUTE X-INTERFACE FLUXES


    IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN
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
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,jm1,dirx,diry,dirz)
      DO k = kmn,kmx
       km1 = max(k-1,1)
       DO j = jmn,jmxb 
        jm1 = max(j-1,1)
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
      if (ixbeg .le. imn)  ixb = imn-ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,im1,dirx,dirz)
      DO k = kmn,kmx
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
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       km1 = k-1
       if(kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uu,vv,ww,jm1,diry,dirz)
      DO k = kmn,kmx
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
      if (ixbeg .le. imn)  ixb = imn-ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,vs,wb)
     DO k = kmn,kmx
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
      if (ixbeg .le. imn)  ixb = imn-ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,uw,wb)
     DO k = kmn,kmx
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
      if (ixbeg .le. imn)  ixb = imn-ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn-jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,vs,wb)
     DO k = kmn,kmx
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
! INTERNAL SUBPROGRAM ADVECT_MONOF1:  Rather complicated general advection algorithm.  See above

    SUBROUTINE ADVECT_MONOF1(fdamp1)
    
    implicit none
    
    logical :: mono
    real    :: fdamp1

! Local vars

    integer         :: i0, im1, im2, ip1, j0, jm1, jm2, jp1, km1, km2, kp1
    real            :: qi, qid, damp, dir
    real            :: vv, cr, qinmax, qoutmax, qoutmin, qinmin, qC, qD, qU
    real            :: dtbydx, dtbydy, dtbydz, del
    real, parameter :: f40 =  7./12., f41 = 1./12
    real, parameter :: f60 = 37./60., f61 = 2./15., f62 = 1./60.
    real, parameter :: eps = 1.0e-15
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1,kmx1
    real :: fd1,fd2
!    real :: gxs(nx), gxw(nx)
    real :: dt2
!    real :: s0x(-ng+1:nx+ng)
!    real :: s0y(-ng+1:ny+ng)
!    real,allocatable :: s0t(:,:,:)

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
         DO n = 1,ng
          s0(1-n   ,1:ny-1,1:nz-1) = s0(nx-n,1:ny-1,1:nz-1)
          s0(nx+n-1,1:ny-1,1:nz-1) = s0(n   ,1:ny-1,1:nz-1)
         ENDDO
        ENDIF

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn)  ixb = imn-ixbeg+1
    if (ixend .ge. imxb) ixe = imxb-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2)
	DO k = kmn,kmx 
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
        IF( (ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3) .or. bcx .eq. 2) THEN
#else
		IF( (i .ge. 4 .and. i .le. i3) .or. bcx .eq. 2 ) THEN
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

	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gxt(i,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		
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
       jmn1 = jmn
       DO n = 1,ng
        s0(1:nx-1,1-n   ,1:nz-1) = s0(1:nx-1,ny-n,1:nz-1)
        s0(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,n   ,1:nz-1)
       ENDDO
     ENDIF

!   write(0,*) 'mono, y-dir , var = ', varname

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmxb) jye = jmxb-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = imn-ixbeg+1
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
	DO k = kmn,kmx 
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
        IF( (jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3) .or. bcy .eq. 2 ) THEN
#else
		IF( ( j .ge. 4 .and. j .le. j3 ) .or. bcy .eq. 2 ) THEN
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
	     
		 qinmin  = min(qD, qC)
		 qinmax  = max(qD, qC)
		 qi      = max(qi, qinmin)
		 qi      = min(qi, qinmax)
  
		 cr      = min(dt2*gyt(j,4)*abs(vv), 1.0)
		 qoutmin = (qC - (1.-cr)*max(qC,qU)) / (cr+eps)
		 qoutmax = (qC - (1.-cr)*min(qC,qU)) / (cr+eps)

		 qi      = min(max(qi, qoutmin), qoutmax)
		 
		
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
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = imn-ixbeg+1
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze

     km1 = max(k-1,1)
     km2 = max(k-2,1)
     kp1 = min(k+1,k1)

     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
     DO k = kmn,kmx

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,k1)

	 DO j = jmn,jmx 
	  DO i = imn,imx 
#endif
         qi  = 0.0
         vv  = 0.5*(w(i,j,k) + w(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

!         IF( k .ge. 4 .and. k .le. k3 ) THEN

          qi = ( f60 * (s(i,j,k  ) + s(i,j,k-1))  &
               - f61 * (s(i,j,k+1) + s(i,j,k-2))  &
               + f62 * (s(i,j,k+2) + s(i,j,k-3))  &
               - f62 * (s(i,j,k+2) - s(i,j,k-3)   &
               - 5.0 * (s(i,j,k+1) - s(i,j,k-2))  &
               + 10. * (s(i,j,k  ) - s(i,j,k-1)))*dir )


!         ELSEIF( k .eq. 3 .or. k .eq. k2 ) THEN
!
!          qi = ( f40 * (s(i,j,k  ) + s(i,j,k-1))  &
!               - f41 * (s(i,j,k+1) + s(i,j,k-2))  &
!               + f41 * (s(i,j,k+1) - s(i,j,k-2)   &
!               - 3.0 * (s(i,j,k  ) - s(i,j,k-1)))*dir )
!
!         ELSEIF( k .eq. 2 .or. k .eq. k1 ) THEN
!
!          qi = 0.5 * (s(i,j,k) + s(i,j,km1))
!    
!!         ELSE
!!         
!!           fz(i,j,k) = 0.0
!!           CYCLE
!         
!         ENDIF
        
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
     
!    deallocate( s0t )

    CALL cld_cpu('ADVECT_MONOF1')


#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_MONOF1: EXITING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif


    RETURN
    END SUBROUTINE ADVECT_MONOF1

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
    real, parameter :: eps = 1.0e-15
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1,kmx1
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
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn)  ixb = imn -ixbeg+1
    if (ixend .ge. imxb) ixe = imxb-ixbeg+1

    DO k = kzb,kze
     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2)
	DO k = kmn,kmx 
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
       DO n = 1,ng
        s0t(1:nx-1,1-n  , 1:nz-1) = s0(1:nx-1,ny-n,1:nz-1)
        s0t(1:nx-1,ny+n-1,1:nz-1) = s0(1:nx-1,n   ,1:nz-1)
       ENDDO
     ENDIF

!   write(0,*) 'mono, y-dir , var = ', varname

#ifdef MPI
    kzb = 1
    kze = ktile+1
    if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn)  jyb = jmn -jybeg+1
    if (jyend .ge. jmxb) jye = jmxb-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = imn-ixbeg+1
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
	DO k = kmn,kmx 
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
    if (kzend .ge. kmx) kze = kmx-kzbeg+1

    jyb = 1
    jye = jtile+1
    if (jybeg .le. jmn) jyb = jmn-jybeg+1
    if (jyend .ge. jmx) jye = jmx-jybeg+1

    ixb = 1
    ixe = itile+1
    if (ixbeg .le. imn) ixb = imn-ixbeg+1
    if (ixend .ge. imx) ixe = imx-ixbeg+1

    DO k = kzb,kze

     km1 = max(k-1,1)
     km2 = max(k-2,1)
     kp1 = min(k+1,k1)

     DO j = jyb,jye
      DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
     DO k = kmn,kmx

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,k1)

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
    real, parameter :: eps = 1.0e-15
    integer         :: imn1,imx1
    integer         :: jmn1,jmx1
    integer         :: kmn1,kmx1
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
       if (ixbeg .le. imn ) ixb = imn -ixbeg+1
       if (ixend .ge. imxb) ixe = imxb-ixbeg+1

       jyb = 1
       jye = jtile+1
       if (jybeg .le. jmn) jyb = jmn-jybeg+1
       if (jyend .ge. jmx) jye = jmx-jybeg+1

       kzb = 1
       kze = ktile+1
       if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
       if (kzend .ge. kmx) kze = kmx-kzbeg+1

       do k = kzb,kze
        do j = jyb,jye
         do i = ixb,ixe
#else

!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,i0,im1,im2,ip1,dtbydx,del,qi,qid,damp,  &
!$OMP         vv,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax,fd1,fd2)
       DO k = kmn,kmx 
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
      
       do k = kzb,kze ; do i = ixb,ixe
        s0t(i,1-n   ,k) = s0t(i,ny-n,k)
        s0t(i,ny+n-1,k) = s0t(i,n   ,k)
       enddo;enddo
#else
        s0t(1:nx-1,1-n, 1:nz) = s0t(1:nx-1,ny-n,1:nz)
        s0t(1:nx-1,ny+n-1,1:nz) = s0t(1:nx-1,n,1:nz)
#endif
      ENDDO
     ENDIF

#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn ) jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

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
    DO k = kmn,kmx 
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
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
         km1    = max(k-1,1)
         km2    = max(k-2,1)
         kp1    = min(k+1,k1)
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,km1,km2,kp1,dtbydx,del,qi,qid,damp, &
!$OMP         vv,dir,qD,qC,qU,qinmin,qinmax,cr,qoutmin,qoutmax)
     DO k = kmn,kmx

     km1    = max(k-1,1)
     km2    = max(k-2,1)
     kp1    = min(k+1,k1)

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
! INTERNAL SUBPROGRAM ADVECT_WENO:  5th-order WENO interpolation

    SUBROUTINE ADVECT_WENO()
    
    implicit none

    integer            :: i, im1, jm1, km1
    real, parameter    :: f30 =  7./12., f31 = 1./12
    real               :: qim2, qim1, qi, qip1, qip2, vv, dir
    real               :: beta0, beta1, beta2, f0, f1, f2, wi0, wi1, wi2, sumwk
    real, parameter    :: gi0 = 1./10., gi1 = 6./10., gi2 = 3./10., eps=1.0e-8
    integer, parameter :: pw = 2

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
    integer :: ixb, jyb, kzb
    integer :: ixe, jye, kze
#endif

!--------------------------------------------------------------------------
! X-INTERFACE

    CALL cld_cpu('ADVECT-WENO')

#ifdef MPI
   if(debug_mpi)  write(0,"('ADVECT_WENO: ENTERING SUBROUTINE, var=',1x,A,1x,', my_rank=',1x,i2)") varname,my_rank
#endif

    IF ( nx .gt. 2 ) THEN
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = imn -ixbeg+1
      if (ixend .ge. imxb) ixe = imxb-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
         im1 = i-1
         if(ixbeg.eq.nxbeg) im1 = max(i-1,1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      DO k = kmn,kmx
       DO j = jmn,jmx
        DO i = imn,imxb
         im1 = max(i-1,1)
#endif
         vv  = 0.5*(u(i,j,k) + u(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF ( ( ixbeg-1+i .ge. 4 .and. ixbeg-1+i .le. i3 ) .or. bcx .eq. 2 ) THEN
#else
         IF ( ( i .ge. 4 .and. i .le. i3 ) .or. bcx .eq. 2 ) THEN
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
    
         f0 =  1./3.*qim2 - 7./6.*qim1 + 11./6.*qi
         f1 = -1./6.*qim1 + 5./6.*qi   + 1./3. *qip1
         f2 =  1./3.*qi   + 5./6.*qip1 - 1./6. *qip2
    
         beta0 = 13./12.*(qim2 - 2.*qim1 + qi  )**2 + 1./4.*(qim2 - 4.*qim1 + 3.*qi)**2
         beta1 = 13./12.*(qim1 - 2.*qi   + qip1)**2 + 1./4.*(qim1 - qip1)**2
         beta2 = 13./12.*(qi   - 2.*qip1 + qip2)**2 + 1./4.*(qip2 - 4.*qip1 + 3.*qi)**2
    
         wi0 = gi0 / (eps + beta0)**pw
         wi1 = gi1 / (eps + beta1)**pw
         wi2 = gi2 / (eps + beta2)**pw
    
         sumwk = wi0 + wi1 + wi2
    
         fx(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF ( ixbeg-1+i .eq. 3 .or. ixbeg-1+i .eq. i2 ) THEN
#else
         ELSEIF ( i .eq. 3 .or. i .eq. i2 ) THEN
#endif
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
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
      if (ixend .ge. imx) ixe = imx-ixbeg+1
    
      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn)  jyb = jmn -jybeg+1
      if (jyend .ge. jmxb) jye = jmxb-jybeg+1
    
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1
    
      do k = kzb,kze
       do j = jyb,jye
        jm1 = j-1
        if(jybeg.eq.nybeg) jm1 = max(j-1,1)
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      DO k = kmn,kmx
       DO j = jmn,jmxb
        jm1 = max(j-1,1)
!        IF ( bcy .eq. 2 ) jm1 = j-1
        DO i = imn,imx  
#endif
         vv  = 0.5*(v(i,j,k) + v(i-is,j-js,k-ks))
         dir = sign(1.0,vv)

#ifdef MPI
         IF( ( jybeg-1+j .ge. 4 .and. jybeg-1+j .le. j3 ) .or. bcy .eq. 2 ) THEN
#else
         IF( ( j .ge. 4 .and. j .le. j3 ) .or. bcy .eq. 2 ) THEN
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
    
         wi0 = gi0 / (eps + beta0)**pw
         wi1 = gi1 / (eps + beta1)**pw
         wi2 = gi2 / (eps + beta2)**pw
    
         sumwk = wi0 + wi1 + wi2
    
         fy(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF ( jybeg-1+j .eq. 3 .or. jybeg-1+j .eq. j2 ) THEN
#else
         ELSEIF ( j .eq. 3 .or. j .eq. j2 ) THEN
#endif
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
#ifdef MPI
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile+1
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze
       km1 = k-1
       if (kzbeg.eq.nzbeg) km1 = max(k-1,1)
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,vv,dir,  &
!$OMP im1, qip2,qip1,qi,qim1,qim2,f0,f1,f2,beta0,beta1,beta2,wi0,wi1,wi2,sumwk)
      DO k = kmn,kmx
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
    
         wi0 = gi0 / (eps + beta0)**pw
         wi1 = gi1 / (eps + beta1)**pw
         wi2 = gi2 / (eps + beta2)**pw
    
         sumwk = wi0 + wi1 + wi2
    
         fz(i,j,k) = vv * (wi0*f0 + wi1*f1 + wi2*f2) / sumwk

#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 3 .or. kzbeg-1+k .eq. k2 ) THEN
#else
         ELSEIF( k .eq. 3 .or. k .eq. k2 ) THEN
#endif
         fz(i,j,k) = vv * ( f30 * (s(i,j,k  ) + s(i,j,k-1))  &
                          - f31 * (s(i,j,k+1) + s(i,j,k-2))  &
                          + f31 * (s(i,j,k+1) - s(i,j,k-2) - 3.*(s(i,j,k)-s(i,j,k-1)))*dir )

#ifdef MPI
         ELSEIF( kzbeg-1+k .eq. 2 .or. kzbeg-1+k .eq. k1 ) THEN
#else
         ELSEIF( k .eq. 2 .or. k .eq. k1 ) THEN
#endif
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
    END SUBROUTINE ADVECT_WENO

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


!--------------------------------------------------------------------------
! SUBPROGRAM ADVECT:  5th-order upwind-biased advection

    SUBROUTINE ADVECT5TH(s,s0,fs,u,v,w,rrp,rrm,dt,nx,ny,nz,      &
                         imn,imx,imx1,jmn,jmx,jmx1,kmn,kmx,kmx1, &
                         imxb, jmxb, &
                         i1,i2,i3,j1,j2,j3,k1,k2,k3,is,js,ks,  &
                         gxt,gyt,gzt,xmask,ymask,varname)

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

      fz2(:,:,:) = 0.0
      kb = 2
      kt = 1

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx1

      do k = kzb,kze 

!
!  Note that x and y loop bounds need to be set _inside_ the k loop
!  because they change for x-, y-, and z-direction fluxes
!
      ixb = 1
      ixe = itile+1
      if (ixbeg .le. imn)  ixb = imn -ixbeg+1
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

    IF ( nx .gt. 2 .and. k .ge. kmn+ks ) THEN
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
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
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

    IF ( ny .gt. 2 .and. k .ge. kmn+ks ) THEN
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
      if (ixbeg .le. imn) ixb = imn-ixbeg+1
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
  if (ixbeg .le. imn)  ixb = imn -ixbeg+1
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
      IF ( k .ge. kmn+ks ) THEN
       fs(i,j,k) = -            (fz2(i  ,j  ,kt)*rrp(k) - fz2(i,j,kb)*rrm(k))*gzt(k,3+ks) &
                   -   ymask(j)*(fy2(i  ,j+1)        - fy2(i,j)       )*gyt(j,3+js) &
                   -   xmask(i)*(fx2(i+1,j  )        - fx2(i,j)       )*gxt(i,3+is) &
                   +  s0(i,j,k) * div2(i,j)
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
