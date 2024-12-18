!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN          \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE TRMMDATDRIVE ////////////////////
!     
!-----------------------------------------------------------------------------
  SUBROUTINE TRMMDATDRIVE(gd,                 &
                    u,  uinit ,         &            ! U,  UINIT
                    v,  vinit ,         &            ! V,  VINIT
                    w,  winit ,         &            ! W,  WINIT
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    s,  sinit,          &            ! S,  SINIT
                    precip,             &            ! PRECIP
                    gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                    gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                    gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                 &            ! DT
                    ugrid, vgrid,       &            ! GRID MOTION
                    nsmall,             &            ! NSMALL
                    microphys,          &            ! microphysical scheme
                    nx, ny, nz, ns,     &            ! NX,NY,NZ,NS
                    io_flag,            &            ! IO_FLAG
                    st, ut, vt, wt, pt, &
                    fu, fv, fw, fp,     &
                    fs, kt, t0,         &
                    dbz, wz,            &
                    elec, cion,         &
                    muz,                &
                    time, tstat, lstt)
 
   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE MICRO_MODULE
   USE PARAM_MODULE
   USE INDEX_MODULE
   USE FORCE_MODULE
   USE TRMM_MODULE
   
!-----------------------------------------------------------------------------
   
   implicit none
 
!-----------------------------------------------------------------------------
! GRID DEFINITIONS

   TYPE(GRID)         :: gd 

   character(LEN = *) :: microphys
   integer            :: nx, ny, nz, ns
   integer            :: dt, nsmall    
   real               :: ugrid, vgrid     
   logical            :: io_flag

   TYPE(VARIABLE)     :: u, uinit
   TYPE(VARIABLE)     :: v, vinit
   TYPE(VARIABLE)     :: w, winit
   TYPE(VARIABLE)     :: pi, piinit
   TYPE(VARIABLE)     :: km, kminit
   TYPE(VARIABLE)     :: s(ns), sinit(2)     ! We are going to try and create arrays of the scalar variables needed here
   TYPE(VARIABLE)     :: precip(5)
   TYPE(VARIABLE)     :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE)     :: dbz, wz
   TYPE(VARIABLE)     :: elec(neelec)
   TYPE(VARIABLE)     :: cion(2)
   TYPE(VARIABLE)     :: muz(4)


   real :: preciptmp(nx,ny,5)
   
   real :: preciptot(5)
      
   real :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)

   real :: ut (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: vt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: wt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: pt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: kt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

   real :: fu (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fv (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fw (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fp (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: fs (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   
   integer :: tstat
   logical :: lstt
      
! Local variables
      
   logical, parameter :: debugsolver = .false.
   integer :: loop, ge, i, j, k, n, nstep, atype
   real    :: dx, dy, dz
   real    :: dts, sdt
   real    :: den(0:nz,2), rrp(nz,2), rrm(nz,2)
   integer :: nsrk,nsfwd
      
   real    :: gtx(nx), gty(ny), gtz(nz)
   real    :: gxt(nx,4), gyt(ny,4), gzt(nz,4)
   integer :: nst = 1                                ! dummy timestep; used to flip xy|yx ordering in Crowley advection
   save nst
   integer :: iadiv
   real    :: dt1
   integer :: nsub

   integer :: ix,jy,kz,ia
   integer :: i1,i2, j1,j2, k1,k2
   integer :: is, js, ks
   integer :: im1, jm1, ip1, jp1
   real    :: wmax, wmin, vmax, a
   integer :: time

   real    :: q(2*lqmx), pii_total, t_total
   real    :: qtodbz
      
   real    :: ice_param(12)
      
   real, allocatable :: sbase(:,:)
      
   real wzz
   
   logical, parameter :: filter = .false.
   integer, parameter :: ihole  = 2

! microphysics flags for number of sub-steps and minimum dbz

      
   real, parameter    :: mindbz = 0.0
!   integer, parameter :: nphys = 1
   integer, parameter :: iuvwadv = 1  ! 1 = normal; 2 = box scheme (Now WORKING)
   
   logical, parameter :: ltkeall = .true. ! .false. 
   
! Set advection options for rksteps:  
!  atypes1 (scalars) and atypem1 (momentum) are for the first steps
!  atypes2 (scalars) and atypem2 (momentum) are for the last step
!
!  IF( ATYPE .eq. 0 ) CALL ADVECT5()
!  IF( ATYPE .eq. 1 ) CALL ADVECT_WENO()
!  IF( ATYPE .eq. 2 ) CALL ADVECT_MONO(0.0,0.0,.false.)
!  IF( ATYPE .eq. 3 ) CALL ADVECT_MONO(0.0,0.0,.true.)
!  IF( ATYPE .eq. 4 ) CALL ADVECT_MONO(damph,0.0,.false.)
!  IF( ATYPE .eq. 5 ) CALL ADVECT_MONO(damph,0.0,.true.)
!  IF( ATYPE .eq. 6 ) CALL ADVECT_MONO(0.0,1.0,.true.)
!  IF( ATYPE .eq. 7 ) CALL ADVECT_MONO(damph,1.0,.false.)
!  IF( ATYPE .eq. 8 ) CALL ADVECT_MONO(damph,1.0,.true.)

!   integer, parameter :: atypes1 = 0
!   integer, parameter :: atypes2 = 5

!   integer, parameter :: atypem1 = 0
!   integer, parameter :: atypem2 = 1

! Vorticity function statements

   wzz(i,j,k,im1,jm1)=(v%flt3d(i,j,k)-v%flt3d(im1,j,k))*gx(4)%flt1d(i)-(u%flt3d(i,j,k)-u%flt3d(i,jm1,k))*gy(4)%flt1d(j)
      
!==================================================================================================================================

  IF ( debugsolver ) print *,'SOLVER3D:  BEGIN '
  IF ( debugsolver ) print *,''
  IF ( debugsolver ) print *,'SOLVER3D:  NX/NY/NZ/NS/NG ', nx, ny, nz, ns, ng
  IF ( debugsolver ) print *,'RKSCHEME = ',RKSCHEME
      
!==================================================================================================================================
! COMPUTE locally need variables

  CALL cld_cpu('TIMESTEP')  

   CALL GET_VARIABLE (gd, 'DX', dx)
   CALL GET_VARIABLE (gd, 'DY', dy)
   CALL GET_VARIABLE (gd, 'DZ', dz)

   IF ( debugsolver ) print *,'SOLVER3D:  DX/DY/DZ ', dx, dy, dz

   allocate ( sbase(nz,ns) )
   sbase(:,:) = 0.0

   DO n = 1,2
     DO k = 1,nz
       sbase(k,n) = sinit(n)%flt1d(k)
     ENDDO
   ENDDO

   DO n = 3,ns
     DO k = 1,nz
      sbase(k,n) = s(n)%base1d(k)
     ENDDO
   ENDDO

   DO k = 1,5
   DO j = 1,ny
     DO i = 1,nx
       preciptmp(i,j,k) = precip(k)%flt2d(i,j)
     ENDDO
   ENDDO 
   ENDDO

! FOR RK/Forward:

   IF ( lfwds ) THEN
     nsrk  = 2
     nsfwd = ns
   ELSE
     nsrk = ns
     nsfwd = 0
   ENDIF

! Constants needed
      
   gtx(:) = gx(3)%flt1d(:)*dx
   gty(:) = gy(3)%flt1d(:)*dy
   gtz(:) = gz(3)%flt1d(:)*dz
       
   DO i = 1,4
     gxt(:,i) = gx(i)%flt1d(:)
     gyt(:,i) = gy(i)%flt1d(:)
     gzt(:,i) = gz(i)%flt1d(:)
   ENDDO
       
   ge  = -ng + 1
   dts = dt / float(nsmall)

! Density variables (den(:,1) = scalar-pt, den(:,2) = w-pt)

   den(1:nz-1,1) = 1.0e5*piinit%flt1d(1:nz-1)**2.509/(rd*sbase(1:nz-1,1) * (1.0+0.61*sbase(1:nz-1,2)))
   den(0,1)      = den(1,1)
   den(nz,1)     = den(nz-1,1)
   den(0:1,2)    = den(0:1,1)
   den(nz,2)     = den(nz,1)
   den(2:nz-1,2) = 0.5*(den(2:nz-1,1) + den(1:nz-2,1))

! Reciprocal density variables for flux calculations

   rrp(:,:)      = 1.0  ! init to 1.0 for all undefined reciprocals
   rrm(:,:)      = 1.0  ! init to 1.0 for all undefined reciprocals
                                                                                                                                
   rrp(1:nz-1,1) = den(2:nz,  2) / den(1:nz-1,1)   ! scalar correction top
   rrm(1:nz-1,1) = den(1:nz-1,2) / den(1:nz-1,1)   ! scalar correction bottom
   rrp(2:nz-1,2) = den(2:nz-1,1) / den(2:nz-1,2)   ! w      correction top
   rrm(2:nz-1,2) = den(1:nz-2,1) / den(2:nz-1,2)   ! w      correction bottom

!==================================================================================================================================
! COPY data into tmp arrays for RK time stepping

   DO k = 1,nz
    ut(1:nx,1:ny,k) = u%flt3d(1:nx,1:ny,k) 
    vt(1:nx,1:ny,k) = v%flt3d(1:nx,1:ny,k) 
    wt(1:nx,1:ny,k) = w%flt3d(1:nx,1:ny,k) 
    pt(1:nx,1:ny,k) = pi%flt3d(1:nx,1:ny,k) 
    kt(1:nx,1:ny,k) = km%flt3d(1:nx,1:ny,k) 
   ENDDO

   DO n = 1,ns
!    IF ( s(n)%name(1:1) .eq. 'V' .or. s(n)%name(1:1) .eq. 'C' .or. s(n)%name(1:2) .eq. 'SC' ) THEN
!
!      s(n)%dyntype = 1
!      
!      DO k = 1,nz-1
!      s(n)%flt3d(1:nx-1,1:ny-1,k) = s(n)%flt3d(1:nx-1,1:ny-1,k)/den(k,1)
!      ENDDO
!    ENDIF
    
    st(1:nx-1,1:ny-1,1:nz-1,n) = s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1)

   ENDDO

  CALL cld_cpu('TIMESTEP')  




!=================================================================================================================================
! PART-VI:  MICROPHYSICS

  
      IF ( microphys .eq. 'KESSLER'  ) THEN
  
       CALL cld_cpu('MICROPHYSICS')
  
        CALL KESSLER(s(1)%flt3d,                         & ! potential temperature
                     s(2)%flt3d,                         & ! water vapor mixing ratio
                     s(3)%flt3d,                         & ! cloud water mixing ratio at t
                     s(4)%flt3d,                         & ! rain    water mixing ratio at t
                     preciptmp,                          & ! precip array
                     sbase(1,1), piinit%flt1d,           & ! base state theta and pi
                     gz(3)%flt1d,                        & ! vertical grid spacing
                     dt,                                 & ! time step
                     fu, fv, fw, pt,                     & ! scratch arrays
                     nx, ny, nz)                           ! grid dimensions
  
        IF ( io_flag ) THEN
  
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
  
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(1)%flt3d(i,j,k)
 
              q(lr)             = s(4)%flt3d(i,j,k) ! st(i,j,k,4)
              dbz%flt3d(i,j,k) = qtodbz(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                                 mindbz, microphys, 0)
  
            ENDDO
           ENDDO
          ENDDO
      
        ENDIF
  
        IF ( debugsolver ) print *,'SOLVER3D:   PAST KESSLER'

        CALL cld_cpu('MICROPHYSICS')
  
      ELSEIF ( microphys .eq. 'WARMLFO' ) THEN
  
       CALL cld_cpu('MICROPHYSICS')

         CALL WARMLFO(s(1)%flt3d,                         & ! potential temperature
                      s(2)%flt3d,                         & ! water vapor mixing ratio
                      s(3)%flt3d,                         & ! cloud water mixing ratio at t
                      s(4)%flt3d,                         & ! rain    water mixing ratio at t
                      preciptmp,                          & ! precip array
                      sbase(1,1), piinit%flt1d,           & ! base state theta and pi
                      gz(3)%flt1d,                        & ! vertical grid spacing
                      dt,                                 & ! time step
                      nx, ny, nz)                           ! grid dimensions

        IF ( io_flag ) THEN
  
          DO k = 1,nz-1
           DO j = 1,ny-1
            DO i = 1,nx-1
  
              pii_total        = piinit%flt1d(k) + pi%flt3d(i,j,k)
              t_total          = s(1)%flt3d(i,j,k) ! st(i,j,k,1)
  
              q(lr)             = s(4)%flt3d(i,j,k) ! st(i,j,k,4)
              dbz%flt3d(i,j,k) = qtodbz(1, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                                 mindbz, microphys, 0)
  
            ENDDO
           ENDDO
          ENDDO
      
        ENDIF
  
        IF ( debugsolver ) print *,'SOLVER3D:   PAST WARMLFO'
  
       CALL cld_cpu('MICROPHYSICS')
 
 
      ELSEIF ( microphys .eq. 'LFO' ) THEN
   
     CALL cld_cpu('MICROPHYSICS')
      ice_param(01) = hole_fill
      ice_param(02) = autoconversion
      ice_param(03) = rho_qr
      ice_param(04) = cnor
      ice_param(05) = rho_qs
      ice_param(06) = cnos
      ice_param(07) = rho_qh
      ice_param(08) = cnoh
      ice_param(09) = qcmincwrn
      ice_param(10) = cwdiap
      ice_param(11) = cwdisp
      ice_param(12) = ccn
  
      DO n = 1,nphys
  
        CALL LFO_ICE_DRIVE(ice_param,                                            &
                           piinit%flt1d, sbase(1,1), sbase(1,2),                 &   ! base state
                           pi%flt3d,                                             &   ! pi'
                           st,                                                   &   ! scalars
                           preciptmp,                                            &   ! precip
                           dt/float(nphys), gz(3)%flt1d,                         &   ! time step, vertical grid spacing
                           nx, ny, nz, ns, ng )
  
      ENDDO
  
      DO n = 1,ns
        s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1) = st(1:nx-1,1:ny-1,1:nz-1,n)
      ENDDO 
  
      IF ( io_flag ) THEN
  
        DO k = 1,nz-1
          DO j = 1,ny-1
          DO i = 1,nx-1
  
            pii_total = piinit%flt1d(k) + pi%flt3d(i,j,k)
            t_total     = st(i,j,k,1)
  
            q(lr) = st(i,j,k,4)
            q(li) = st(i,j,k,5)
            q(ls) = st(i,j,k,6)
            q(lh) = st(i,j,k,7)
            dbz%flt3d(i,j,k) = qtodbz(4, q, pii_total, t_total, cnor, rho_qr, cnos, rho_qs, cnoh, rho_qh, &
                               mindbz, microphys, 0)
  
          ENDDO
          ENDDO
        ENDDO
      
      ENDIF
  
      IF ( debugsolver ) print *,'SOLVER3D:   PAST LFO_ICE_DRIVE'
      CALL cld_cpu('MICROPHYSICS')
  
   ELSEIF( microphys(1:5) .eq. 'ICE10' .or.        &
           microphys(1:4) .eq. 'ZIEG'  .or.        &
           microphys(1:3) .eq. 'ZVD'   .or.        &
           microphys(1:8) .eq. 'WARMZIEG' ) THEN
  
      IF ( debugsolver ) print *,'SOLVER3D:   CALL ICE10_DRIVE'
      
      dt1 = dt/float(nphys)

   
      IF ( microphys(1:5) .eq. 'ICE10' ) THEN
      CALL ICE10_DRIVE(gd,                                                  &
                       piinit%flt1d, sbase,                                 &     ! base state
                       pi%flt3d,                                            &     ! pi'
                       st,                                                  &     ! scalars
                       preciptmp,                                           &     ! precip
                       dt1,                                                 &     ! time step
                       gx(3)%flt1d,                                         &     ! horizontal grid spacing
                       gy(3)%flt1d,                                         &     ! horizontal grid spacing
                       gz(3)%flt1d,                                         &     ! vertical grid spacing
                       gzt,                                                 &     ! vertical grid height
                       nx, ny, nz, ng, ns,                                  &
                       u%flt3d(1,1,1),v%flt3d(1,1,1),w%flt3d(1,1,1),        &     ! vn (u=1,v=2,w=3)
                       t0,ut,vt,wt,pt,fu,fv,fw,fp,fs,                       &     ! temporary arrays
                       dx,dy,dz,                                            &     ! average dx,dy,dz
                       nphys,luno,1,1,1,                                    &     ! istag,jstag,kstag
                       dbz%flt3d,km%flt3d, io_flag, time, dt,               &
                       elec, cion, muz, gxt, gyt, gzt,                      &
                       tstat, lstt, bcx, bcy)

      ELSEIF ( microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) .eq. 'ZVD' ) THEN
      
      CALL ZIEG_DRIVE(gd,                                                   &
                      piinit%flt1d, sbase,                                  &     ! base state
                      pi%flt3d,                                             &     ! pi'
                      st,                                                   &     ! scalars
                      preciptmp,                                            &     ! precip
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,                                                  &     ! vertical grid height
                      nx, ny, nz, ng, ns,                                   &
                      u%flt3d(1,1,1),v%flt3d(1,1,1),w%flt3d(1,1,1),         &     ! vn (u=1,v=2,w=3)
                      t0,ut,vt,wt,pt,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d,km%flt3d, io_flag, time, dt,                &
                      elec, cion, muz, gxt, gyt, gzt,                       &
                      tstat, lstt, bcx, bcy)
!                      tstat, lstt, 1, 1)

      ELSEIF  ( microphys(1:8) .eq. 'WARMZIEG' ) THEN
      
      CALL WARMZIEG_DRIVE(gd,                                               &
                      piinit%flt1d, sbase,                                  &     ! base state
                      pi%flt3d,                                             &     ! pi'
                      st,                                                   &     ! scalars
                      preciptmp,                                            &     ! precip
                      dt1,                                                  &     ! time step
                      gx(3)%flt1d,                                          &     ! horizontal grid spacing
                      gy(3)%flt1d,                                          &     ! horizontal grid spacing
                      gz(3)%flt1d,                                          &     ! vertical grid spacing
                      gzt,                                                  &     ! vertical grid height
                      nx, ny, nz, ng, ns,                                   &
                      u%flt3d(1,1,1),v%flt3d(1,1,1),w%flt3d(1,1,1),         &     ! vn (u=1,v=2,w=3)
                      t0,ut,vt,wt,pt,fu,fv,fw,fp,fs,                        &     ! temporary arrays
                      dx,dy,dz,                                             &     ! average dx,dy,dz
                      nphys,luno,u%istag,v%jstag,w%kstag,                   &     ! istag,jstag,kstag
                      dbz%flt3d,km%flt3d, io_flag, time, dt,                &
                      elec, cion, muz, gxt, gyt, gzt,                       &
                      tstat, lstt, bcx, bcy)

      ENDIF
      
      DO n = 1,ns
        s(n)%flt3d(1:nx-1,1:ny-1,1:nz-1) = st(1:nx-1,1:ny-1,1:nz-1,n)
      ENDDO 
      
  
      IF ( debugsolver ) print *,'SOLVER3D:   PAST ICE10_DRIVE'
   ELSE
   
     write(0,*) 'Unsupported microphysics option: ',microphys
     write(0,*) 'STOP'
     STOP
     
   ENDIF
  
  
   IF ( debugsolver ) print *,'SOLVER3D:    EXITING...'
  
!=================================================================================================================================
! Clean up some stuff...

   deallocate ( sbase )

   preciptot(:) = 0.0
   DO j = ibsd,iesd ! 1,ny-1
     DO i = jbsd,jesd ! 1,nx-1
       a = 1./(gx(3)%flt1d(i) * gy(3)%flt1d(j))
       DO k = 1,5
       precip(k)%flt2d(i,j) = preciptmp(i,j,k)
       preciptot(k) = preciptot(k) + preciptmp(i,j,k)*a
       ENDDO
     ENDDO 
   ENDDO
   
   IF ( io_flag ) THEN
     write(luno,*) 'Total rainfall = ',preciptot(3)
     write(luno,*) 'Total hailfall = ',preciptot(4)
   ENDIF
 
 
  CALL emagmax  (nx,ny,nz,elec,gxt,gyt,gzt,luno)  


  RETURN
    
!-----------------------------------------------------------------------------
!
!   /////////////////////           END           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE TRMMDATDRIVE ////////////////////
!     
!-----------------------------------------------------------------------------

   END SUBROUTINE TRMMDATDRIVE
   
   
!
!
!  START OF EFIELDSTAG: Calculates efield on staggered C-grid
!
!
      subroutine emagmax   &
         (nx,ny,nz   &
         ,elec   &
         ,gxt,gyt,gzt   &
         ,iunit)  
     
      USE INDEX_MODULE
      USE GRID_MODULE
!
!
! 7.23.2002: replaced 1st order centered difference for Ez at k=1
!            with a 2nd order forward difference.
!
! 5.13.2002: update for vertical stretching
!
      implicit none
!
      integer nx,ny,nz

      TYPE(VARIABLE)     :: elec(neelec)

      real gxt(nx,4)
      real gyt(ny,4)
      real gzt(nz,4)



      integer ix,jy,kz
      integer, parameter :: istag = 1, jstag = 1, kstag = 1
      
      integer iunit
      
      integer ibg  ! =1 for phi=0 at k=1, -1 for phi=0 at k=0
      
      integer imx,jmx,kmx, imx2, jmx2, kmx2
!
      real emag
      real emax, emax2
      real emin
!      real dx,dy,dz
!      real dxi,dyi,dzi
!      real dxi2,dyi2,dzi2
      
      integer im1,ip1
      integer jm1,jp1
      integer km1,kp1


!
      emax  = -1.0
      emax2 = -1.0
      emin  = 0.0
      do kz = 1,nz-kstag
      do jy = 1,ny-jstag
      do ix = 1,nx-istag
      IF ( emax .lt. elec(iemag)%flt3d(ix,jy,kz) ) THEN
        imx = ix
        jmx = jy
        kmx = kz
      ENDIF
      emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
      emin = min(emin,elec(iemag)%flt3d(ix,jy,kz))

       emag =    &
        Sqrt( (( elec(iex)%flt3d(ix,jy,kz) )**2   &
            +  ( elec(iey)%flt3d(ix,jy,kz) )**2   &
            +  ( elec(iez)%flt3d(ix,jy,kz) )**2 )) 

      IF ( emax2 .lt. emag ) THEN
        emax2 = emag
        imx2 = ix
        jmx2 = jy
        kmx2 = kz
      ENDIF

!        elec(iemag)%flt3d(ix,jy,kz) = Max( elec(iemag)%flt3d(ix,jy,kz), emag )

      end do
      end do
      end do
!
      write(iunit,*) 'EFIELD--E-MAX=',emax
      write(iunit,*) 'EFIELD--E-MIN=',emin
      write(iunit,*) 'ijk,Ez at emax: ',imx,jmx,kmx,elec(iez)%flt3d(imx,jmx,kmx),   &
                      0.5*( elec(iez)%flt3d(imx,jmx,kmx) + elec(iez)%flt3d(imx,jmx,kmx+1) )

      write(iunit,*) 'EFIELD2--E-MAX=',emax2
      write(iunit,*) 'ijk,Ez at emax2: ',imx2,jmx2,kmx2,elec(iez)%flt3d(imx2,jmx2,kmx2)
!
!
      return
      end
!
!
!  END OF EFIELD
!
!
