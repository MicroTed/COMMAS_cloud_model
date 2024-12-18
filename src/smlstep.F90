!-----------------------------------------------------------------------------
! 
! SUBROUTINE SMLSTEP
!
! Latest update:  01-18-00
!
!
!-----------------------------------------------------------------------------
      SUBROUTINE SMLSTEP(ut, u, fu,           &
                         vt, v, fv,           &
                         wt, w, fw,           &
                         pt, p, fp,           &
                         pinit,tinit,qinit,uinit,vinit,     &
                         gx,gy,gz,dts,nx,ny,nz,nsmall)

      USE PARAM_MODULE
      USE CPUTIME_MODULE
      USE COMMASMPI_MODULE
      implicit none

! Passed variables

      integer nx, ny, nz, nsmall

      real ut(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real vt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real wt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real pt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real p (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fp(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real gx(-ng+1:nx+ng,4), gy(-ng+1:ny+ng,4), gz(-ng+1:nz+ng,4)
      real pinit(-ng+1:nz+ng), tinit(-ng+1:nz+ng), qinit(-ng+1:nz+ng)
      real uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng)
      real dts

! Local declarations

      integer i, j, k, n, ip1, im1, jm1, km1
      integer nx1, ny1
      real dz0, cs2, eps0, eps1,tvtmp,rtmp,tvinit,tvhalf
      real atri, btri, ctri(-ng+1:nz+ng), etri(-ng+1:nz+ng), gtri(-ng+1:nz+ng)
      real uw, ue, vn, vs
      real kdivx(-ng+1:nx+ng), kdivy(-ng+1:ny+ng), kdivz(-ng+1:nz+ng)
!      real rhalfi(-ng+1:nz+ng), rinit(-ng+1:nz+ng), rhalf(-ng+1:nz+ng)
      real,allocatable,save,dimension(:) :: rinit, rhalf, rhalfi
      real, allocatable :: fpp(:,:,:), div(:,:,:)
!      real coefu(-ng+1:nz+ng), coefv(-ng+1:nz+ng), coefw(-ng+1:nz+ng)
!      real coefpx(-ng+1:nz+ng)
!      real coefpy(-ng+1:nz+ng)
!      real coefpz(-ng+1:nz+ng)
      real,allocatable,save,dimension(:) :: coefu, coefv, coefw
      real,allocatable,save,dimension(:) :: coefpx
      real,allocatable,save,dimension(:) :: coefpy
      real,allocatable,save,dimension(:) :: coefpz
      real,allocatable :: temp2d(:,:)
      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng)
      
      logical :: ctest = .false.
      logical :: fucomm = .false.
      logical :: explicit = .true.
!      logical :: explicit = .false.
      
      integer nb ! number of ghost zones to communicate
      integer ia

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      integer :: downward_tag, upward_tag

      logical :: debug_mpi = .false.
      
   if (debug_mpi) write(0,"('SMLSTEP: ENTERING SUBROUTINE, my_rank=',1x,i2)") my_rank

!-----------------------------------------------------------------------------
! Bnd flags
       
       explicit = ( .not. vert_implicit )
       IF ( nprock > 1 ) explicit = .true.
       

#ifdef MPI
      allocate ( fpp(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
      allocate ( div(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
#else
      allocate ( fpp(0:nx,0:ny,nz) )
      allocate ( div(0:nx,0:ny,0:nz) )
#endif

      nx1 = nx - 1
      ny1 = ny - 1
      
  ! For periodic BC
      IF ( bcx .eq. 2 ) nx1 = nx
      IF ( bcy .eq. 2 ) ny1 = ny
      
#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      DO i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO 

      if(ixbeg.eq.nxbeg) xbnd(-ng+1:1)  = 0.0 
      if(ixend.eq.nxend) xbnd(nx:nx+ng) = 0.0

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg

      DO j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO 

      if(jybeg.eq.nybeg) ybnd(-ng+1:1)  = 0.0
      if(jyend.eq.nyend) ybnd(ny:ny+ng) = 0.0
#else
      DO i = 2,nx-1 ; xbnd(i) = 1.0 ; ENDDO 
      xbnd(1) = 0.0 ; xbnd(nx) = 0.0

      DO j = 2,ny-1 ; ybnd(j) = 1.0 ; ENDDO 
      ybnd(1) = 0.0 ; ybnd(ny) = 0.0
#endif

      cs2  = cspd**2
      eps0 = 0.5 * (1.0 + alpha)
      eps1 = 0.5 * (1.0 - alpha)

!------------------------------------------------------------------------------
! Initialize need vertical coefficients

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: set coefficients'

      IF ( .not. allocated( rinit ) ) THEN
      
      allocate( rinit(-ng+1:nz+ng), rhalf(-ng+1:nz+ng), rhalfi(-ng+1:nz+ng))
      allocate( coefu(-ng+1:nz+ng), coefv(-ng+1:nz+ng), coefw(-ng+1:nz+ng))
      allocate( coefpx(-ng+1:nz+ng))
      allocate( coefpy(-ng+1:nz+ng))
      allocate( coefpz(-ng+1:nz+ng))

      ENDIF

          rinit(:)  = 0
          rhalf(:)  = 0
          rhalfi(:) = 0
          coefu(:)  = 0
          coefv(:)  = 0
          coefw(:)  = 0
          coefpx(:) = 0
          coefpy(:) = 0
          coefpz(:) = 0
      
      rhalf(:) = 1.0
#ifdef MPI
      kzb = -ng+1
      kze = ktile+ng-1
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
      do k=kzb,kze
       km1 = Max( k-1, -ng+1 )
       if(kzbeg .eq. nzbeg) km1 = max(k-1,1)
#else
      DO k = 1,nz-1
       km1 = max(k-1,1)
#endif

       tvinit   = tinit(k)   * (1.0 + 0.61*qinit(k)  )
       tvtmp    = tinit(km1) * (1.0 + 0.61*qinit(km1))
       tvhalf   = 0.5*(tvtmp + tvinit)

       rinit(k) = pinit(k)  **cvr * p00 / (rd*tvinit)
       rtmp     = pinit(km1)**cvr * p00 / (rd*tvtmp)
       rhalf(k) = 0.5*(rtmp + rinit(k))
       rhalfi(k) = 1.0/rhalf(k)

! Quasi compressible coefficients

!      coefu(k) = dts
!      coefv(k) = dts
!      coefw(k) = dts*gz(k,4)
!      coefpx(k)= dts*cs2
!      coefpy(k)= dts*cs2
!      coefpz(k)= dts*cs2*gz(k,3)

! Pi coefficients

       coefu(k) = cp * tvinit * dts 
       coefv(k) = cp * tvinit * dts 
       coefw(k) = cp * tvhalf * rhalf(k) * dts * gz(k,4)
       coefpx(k)= dts*cs2/(cp * tvinit)
       coefpy(k)= dts*cs2/(cp * tvinit)
       coefpz(k)= dts*cs2/(cp * tvinit * rinit(k))*gz(k,3)

!      write(91,*) 'k,coefpz,rinit,pinit = ',k,coefpz(k),rinit(k),pinit(k)
      ENDDO
      
#ifdef MPI
        IF ( nprock > 1 ) THEN
          allocate( temp2d(-ng+1:nz+ng,9) )
          
          temp2d(:,1) = rinit(:)
          temp2d(:,2) = rhalf(:)
          temp2d(:,3) = rhalfi(:)
          temp2d(:,4) = coefu(:)
          temp2d(:,5) = coefv(:)
          temp2d(:,6) = coefw(:)
          temp2d(:,7) = coefpx(:)
          temp2d(:,8) = coefpy(:)
          temp2d(:,9) = coefpz(:)
          
       downward_tag = 10001
       CALL sendrecv_downward(1,1,nz,0,0,ng,ng,9,  &
         d_proc(my_rank),u_proc(my_rank),downward_tag, temp2d(-ng+1,1) )

       upward_tag = 10002
       CALL sendrecv_upward(1,1,nz,0,0,ng,ng,9,  &
         d_proc(my_rank),u_proc(my_rank),upward_tag, temp2d(-ng+1,1) )
        
          rinit(:)  = temp2d(:,1)
          rhalf(:)  = temp2d(:,2)
          rhalfi(:) = temp2d(:,3)
          coefu(:)  = temp2d(:,4)
          coefv(:)  = temp2d(:,5)
          coefw(:)  = temp2d(:,6)
          coefpx(:) = temp2d(:,7)
          coefpy(:) = temp2d(:,8)
          coefpz(:) = temp2d(:,9)
        
          deallocate( temp2d )
          
!          DO k = -ng+1,nz+ng
!            write(luno,*) 'k,rhalf = ',k,rhalf(k)
!          ENDDO
          
        ENDIF
#endif
      
      
!      ENDIF ! allocated

      etri(1) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
      do k=kzb,kze
#else
      DO k = 2,nz-1
#endif
       atri    =-coefw(k)*coefpz(k)    * eps0**2
       ctri(k) =-coefw(k)*coefpz(k-1)  * eps0**2
       btri    = 1.0 - (atri + ctri(k))
       gtri(k) = 1.0 / (btri + ctri(k)*etri(k-1))
       etri(k) = -atri * gtri(k)
      ENDDO 

!------------------------------------------------------------------------------
!   kdiv  = non-dimensional coefficient for divergent damping (x/y/z)
!           kdiv(k) = kdiv0 * dx / dts so the div-damping = kdiv0 per
!           small time step.

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (myproci == 1 .and. bcx .ne. 2) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2) ixe = ixend-ixbeg

      jyb = -ng+1
      jye = jtile+ng
      if (myprocj == 1 .and. bcy .ne. 2) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2) jye = jyend-jybeg

      kzb = -ng+1
      kze = ktile+ng
      if (myprock == 1) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

!      kdivx(:) = kdiv0/(gx(1,3)*dts)
!      kdivy(:) = kdiv0/(gy(1,3)*dts)
!      kdivz(-ng+1:nz+ng) = kdiv0/(gz(-ng+1:nz+ng,3)*dts)

      DO i = ixb,ixe ; kdivx(i)  = kdiv0 / (gx(i,3)*dts) ; ENDDO
      DO j = jyb,jye ; kdivy(j)  = kdiv0 / (gy(j,3)*dts) ; ENDDO
      DO k = kzb,kze ; kdivz(k)  = kdiv0 / (gz(k,3)*dts) ; ENDDO
#else
      DO i = 1,nx-1
       kdivx(i)  = kdiv0 / (gx(i,3)*dts)
      ENDDO
                                                                                                                          
      DO j = 1,ny-1
       kdivy(j)  = kdiv0 / (gy(j,3)*dts)
      ENDDO
                                                                                                                          
      DO k = 1,nz-1
       kdivz(k)  = kdiv0 / (gz(k,3)*dts)
      ENDDO
#endif


   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: done coefficients'

! Initialize edges of DIV = 0

      fpp(:,:,:) = 0.0
#ifdef MPI
#else
!      div(:, 0,:) = 0.0
!      div(:,ny,:) = 0.0
!      div(:,:, 0) = 0.0
      div(:,:,nz) = 0.0
!      div( 0,:,:) = 0.0
!      div(nx,:,:) = 0.0
#endif

#ifdef MPI
    IF ( number_of_processes .gt. 1 .and. fucomm ) THEN
    
    CALL cld_cpu('MPI-COMM-SMLSTEP')


! fu,fv,fw velocity tendencies -- note that they are contiguous in memory

        
        nb = 1
        ia = 3

        IF ( nproci > 1 ) THEN
        
        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,fu)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,fu)
        
        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,fu)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,fu)

        ENDIF

        IF ( nprock > 1 ) THEN

        downward_tag = 2005
        CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),downward_tag,fu)

        upward_tag = 2006
        CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),upward_tag,fu)
   
        ENDIF

    CALL cld_cpu('MPI-COMM-SMLSTEP')
    
    ENDIF

#endif



!-----------------------------------------------------------------
! Set up tendency arrays

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: set fu'

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg+1

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg+1

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
     DO k = 1,nz-1
       DO j = 1,ny
        DO i = 1,nx
#endif
         fu(i,j,k) = fu(i,j,k) * dts
         fv(i,j,k) = fv(i,j,k) * dts
         wt(i,j,k) = wt(i,j,k) * rhalf(k)
         fw(i,j,k) = fw(i,j,k) * rhalf(k) * dts

        ENDDO
       ENDDO
      ENDDO



   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: set u'

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng

      jyb = -ng+1
      jye = jtile+ng

      kzb = -ng+1
      kze = ktile+ng

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      DO k = 1,nz
       DO j = 1,ny
        DO i = 1,nx
#endif

          u(i,j,k) = 0.0
          v(i,j,k) = 0.0
          w(i,j,k) = 0.0
          p(i,j,k) = 0.0
        div(i,j,k) = 0.0

        ENDDO
       ENDDO
      ENDDO

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: set fpp'

!-----------------------------------------------------------------------
! DO SMLSTEP USING PERTUBATIONS FROM V(t) 
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb, jye ; do i = ixb, ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          fpp(i,j,k)  = coefpx(k)*(ut(i+1,j,  k  ) - ut(i,j,k))*gx(i,3)  &
                      + coefpy(k)*(vt(i,  j+1,k  ) - vt(i,j,k))*gy(j,3)  &
                      + coefpz(k)*(wt(i  ,j,  k+1) - wt(i,j,k))

         ENDDO
        ENDDO
       ENDDO

       IF ( bcx .eq. 2 .and. nproci == 1 ) THEN
#ifdef MPI
        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg

        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        do k = kzb,kze ; do j = jyb, jye
         fpp(0    ,j,k) = fpp(nxend-1,j,k)
         fpp(nxend,j,k) = fpp(1      ,j,k)
        enddo ; enddo
#else
         fpp(0 ,1:ny-1,1:nz-1) = fpp(nx-1,1:ny-1,1:nz-1)
         fpp(nx,1:ny-1,1:nz-1) = fpp(1   ,1:ny-1,1:nz-1)
#endif
       ENDIF

       IF ( bcy .eq. 2 .and. nprocj == 1 ) THEN
#ifdef MPI
        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        do k = kzb,kze ; do i = ixb, ixe
         fpp(i,0    ,k) = fpp(i,nyend-1,k)
         fpp(i,nyend,k) = fpp(i,1      ,k)
        enddo ; enddo
#else
         fpp(1:nx-1,0 ,1:nz-1) = fpp(1:nx-1,ny-1,1:nz-1)
         fpp(1:nx-1,ny,1:nz-1) = fpp(1:nx-1,1   ,1:nz-1)
#endif
       ENDIF


   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: update fu'

#ifdef MPI
      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = 0
      kze = ktile+1
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb, jye ; do i = ixb, ixe
          IF( ixbeg-1+i .ne. 1 .or. bcx .eq. 2 ) fu(i,j,k) = fu(i,j,k)                           &
                                   - coefu(k) * (pt(i,j,k)-pt(i-1,j,k))*gx(i,4)  &
                                   + kdivx(i) * (fpp(i,j,k)-fpp(i-1,j,k))

          IF( jybeg-1+j .ne. 1 .or. bcy .eq. 2  ) fv(i,j,k) = fv(i,j,k)                           &
                                   - coefv(k) * (pt(i,j,k)-pt(i,j-1,k))*gy(j,4)  &
                                   + kdivy(j) * (fpp(i,j,k)-fpp(i,j-1,k))

          IF( kzbeg-1+k .ne. 1 ) fw(i,j,k) = fw(i,j,k)                           &
                                   - coefw(k) * (pt(i,j,k) - pt(i,j,k-1))        &
                                   + kdivz(k) * (fpp(i,j,k) - fpp(i,j,k-1))*rhalf(k)
      enddo ; enddo ; enddo
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       DO k = 1,nz-1 
        DO j = 1,ny-1
         DO i = 1,nx-1

          IF( i .ne. 1 ) fu(i,j,k) = fu(i,j,k)                                   &
                                   - coefu(k) * (pt(i,j,k)-pt(i-1,j,k))*gx(i,4)  &
                                   + kdivx(i) * (fpp(i,j,k)-fpp(i-1,j,k))

          IF( j .ne. 1 ) fv(i,j,k) = fv(i,j,k)                                   &
                                   - coefv(k) * (pt(i,j,k)-pt(i,j-1,k))*gy(j,4)  &
                                   + kdivy(j) * (fpp(i,j,k)-fpp(i,j-1,k))

          IF( k .ne. 1 ) fw(i,j,k) = fw(i,j,k)                                   &
                                   - coefw(k) * (pt(i,j,k) - pt(i,j,k-1))        &
                                   + kdivz(k) * (fpp(i,j,k) - fpp(i,j,k-1))*rhalf(k)
         ENDDO
        ENDDO
       ENDDO
#endif



       IF ( bcx .eq. 2 .and. nproci == 1 ) THEN
#ifdef MPI
      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb, jye
#else
       DO k = 1,nz-1 
        DO j = 1,ny-1
#endif

#ifdef MPI
         if (ixbeg .eq. nxbeg) then
#endif
         i = 1

           fu(i,j,k) = fu(i,j,k)      &
                                   - coefu(k) * (pt(i,j,k)-pt(nx-1,j,k))*gx(i,4)  &
                                   + kdivx(i) * (fpp(i,j,k)-fpp(nx-1,j,k))
#ifdef MPI
         endif
         
         if (ixend .eq. nxend) then
         i = nxend
#else
         i = nx
#endif
            fu(i,j,k) = fu(1,j,k)
!           fu(i,j,k) = fu(i,j,k)      &
!                                   - coefu(k) * (pt(1,j,k)-pt(i-1,j,k))*gx(1,4)  &
!                                   + kdivx(1) * (fpp(i,j,k)-fpp(i-1,j,k))
#ifdef MPI
         endif
#endif
        ENDDO
       ENDDO

       ENDIF !! ( bcx .eq. 2 )

       IF ( bcy .eq. 2 .and. nprocj == 1 ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do i = ixb, ixe
#else
       DO k = 1,nz-1 
        DO i = 1,nx-1
#endif

#ifdef MPI
          if (jybeg .eq. nybeg) then
#endif
          j = 1

          fv(i,j,k) = fv(i,j,k)      &
                                   - coefv(k) * (pt(i,j,k)-pt(i,ny-1,k))*gy(j,4)  &
                                   + kdivy(j) * (fpp(i,j,k)-fpp(i,ny-1,k))

#ifdef MPI
         endif

         if (jyend .eq. nyend) then
          j = nyend
#else
          j = ny
#endif
          fv(i,j,k) = fv(i,1,k)
!          fv(i,j,k) = fv(i,j,k)      &
!                                   - coefv(k) * (pt(i,1,k)-pt(i,j-1,k))*gy(1,4)  &
!                                   + kdivy(1) * (fpp(i,j,k)-fpp(i,j-1,k))
#ifdef MPI
         endif
#endif
        ENDDO
       ENDDO

       ENDIF !! ( bcy .eq. 2 )


! Compute outflow boundaries tendencies

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze

       IF ( nx .gt. 2 .and. bcx .ne. 2 ) THEN

        jyb = -1
        jye = jtile+2
        IF (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        do j = jyb, jye
          if(ixbeg.eq.nxbeg) uw = dts * amin1( ( ut(1 ,j,k) - dxt ), 0.0 ) *  gx( 2  ,4)
          if(ixend.eq.nxend) ue = dts * amax1( ( ut(nx,j,k) + dxt ), 0.0 ) *  gx(nx-1,4)
!          if(ixbeg.eq.nxbeg) fu(1 ,j,k) = - uw * (ut(2 ,j,k) - ut( 1  ,j,k))
!          if(ixend.eq.nxend) fu(nx,j,k) = - ue * (ut(nx,j,k) - ut(nx-1,j,k))
          if(ixbeg.eq.nxbeg) fu(1 ,j,k) = fu(1 ,j,k) - uw * (ut(2 ,j,k) - ut( 1  ,j,k))
          if(ixend.eq.nxend) fu(nx,j,k) = fu(nx,j,k) - ue * (ut(nx,j,k) - ut(nx-1,j,k))
        enddo

       ENDIF

       IF ( ny .gt. 2 .and. bcy .ne. 2 ) THEN
        ixb = -1
        ixe = itile+2
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do i = ixb, ixe
          if(jybeg.eq.nybeg) vs = dts * amin1( ( vt(i,1, k) - dxt ), 0.0 ) * gy(2,   4)
          if(jyend.eq.nyend) vn = dts * amax1( ( vt(i,ny,k) + dxt ), 0.0 ) * gy(ny-1,4)
!          if(jybeg.eq.nybeg) fv(i, 1,k) = - vs * (vt(i,2, k) - vt(i,  1 ,k))
!          if(jyend.eq.nyend) fv(i,ny,k) = - vn * (vt(i,ny,k) - vt(i,ny-1,k))
          if(jybeg.eq.nybeg) fv(i, 1,k) = fv(i, 1,k) - vs * (vt(i,2, k) - vt(i,  1 ,k))
          if(jyend.eq.nyend) fv(i,ny,k) = fv(i,ny,k) - vn * (vt(i,ny,k) - vt(i,ny-1,k))
        enddo

       ENDIF

      ENDDO

!-----------------------------------------------------------------------
! Small step loop

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: start nsmall loop'

      DO n = 1,nsmall


      IF ( explicit ) THEN

      
! Compute divergence; need div in the ghost zone, too, for MPI, to get gradient.
#ifdef MPI
       ixb = 0
       ixe = itile+1
       if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
       if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

       jyb = 0
       jye = jtile+1
       if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
       if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

       kzb = 0
       kze = ktile+1
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          div(i,j,k)  = coefpx(k)*(u(i+1,j,  k  ) - u(i,j,k))*gx(i,3)  &
                      + coefpy(k)*(v(i  ,j+1,k  ) - v(i,j,k))*gy(j,3)  &
                      + coefpz(k)*(w(i  ,j,  k+1) - w(i,j,k))

          p(i,j,k)    = p(i,j,k) - fpp(i,j,k) - div(i,j,k) 

         ENDDO
        ENDDO
       ENDDO

      
      ELSE
! Compute divergence; need div in the ghost zone, too, for MPI, to get gradient.
#ifdef MPI
       ixb = 0
       ixe = itile+1
       if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
       if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

       jyb = 0
       jye = jtile+1
       if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
       if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

       kzb = 0
       kze = ktile+1
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          div(i,j,k)  = coefpx(k)*(u(i+1,j,  k  ) - u(i,j,k))*gx(i,3)  &
                      + coefpy(k)*(v(i  ,j+1,k  ) - v(i,j,k))*gy(j,3)  &
                      + coefpz(k)*(w(i  ,j,  k+1) - w(i,j,k))

          fp(i,j,k)   = p(i,j,k) - fpp(i,j,k) - div(i,j,k) + coefpz(k)*(w(i,j,k+1) - w(i,j,k))*eps0

         ENDDO
        ENDDO
       ENDDO
       
       ENDIF

!-----------------------------------------------------------------------
! Update pressure using w and u data


      
      IF ( .not. explicit ) THEN

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO k = kzb,kze
       DO j = jyb,jye
        DO i = ixb,ixe
#else
!$OMP PARALLEL DO IF ( ny .gt. 2 ), DEFAULT(SHARED), PRIVATE(i,j,k)
      DO j = 1,ny-1
       DO k = 2,nz-1
        DO i = 1,nx-1
#endif
! Compute RHS for tridiagonal system

          w(i,j,k) = w(i,j,k) + fw(i,j,k)   &
                    - coefw(k) * (eps1*( p(i,j,k) -  p(i,j,k-1))   &
                                 +eps0*(fp(i,j,k) - fp(i,j,k-1)))  &
                    + kdivz(k)*rhalf(k)*(div(i,j,k) - div(i,j,k-1))

! Upward sweep for forward substitution

          w(i,j,k) = (w(i,j,k)-ctri(k)*w(i,j,k-1))*gtri(k)

         ENDDO
        ENDDO
       ENDDO


! Downward sweep for backward substitution + update p and t

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = ktile
      kze = 1
      if (kzbeg .eq. nzbeg) kze = 1
      if (kzend .eq. nzend) kzb = kzend-kzbeg

      do k = kzb,kze,-1 ; do j = jyb, jye ; do i = ixb, ixe
#else
!$OMP PARALLEL DO IF ( ny .gt. 2 ), DEFAULT(SHARED), PRIVATE(i,j,k)
       DO j = 1,ny-1
        DO k = nz-1,1,-1
         DO i = 1,nx-1
#endif
          w(i,j,k) = w(i,j,k+1) * etri(k) + w(i,j,k)
          p(i,j,k) = fp(i,j,k) - eps0*coefpz(k)*(w(i,j,k+1)-w(i,j,k))

         ENDDO
        ENDDO

       ENDDO
       
       ENDIF ! .not. explicit

! MPI: Send pressure update to other processes

#ifdef MPI

    IF ( number_of_processes .gt. 1 ) THEN

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: comms on p'
    
    CALL cld_cpu('MPI-COMM-SMLSTEP')


        nb = 1
        IF ( n == nsmall ) nb = 1
        
        IF ( nproci > 1 ) THEN

        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,p)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,1,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,p)

        ENDIF

        IF ( nprocj > 1 ) THEN
        
        southward_tag = 203
        CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,p)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,1,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,p)

        ENDIF

        IF ( nprock > 1 ) THEN

        downward_tag = 2005
        CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,1,  &
           d_proc(my_rank),u_proc(my_rank),downward_tag,p)

        upward_tag = 2006
        CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,1,  &
           d_proc(my_rank),u_proc(my_rank),upward_tag,p)
   
        ENDIF

    CALL cld_cpu('MPI-COMM-SMLSTEP')
    
    ENDIF

#endif

      IF ( explicit ) THEN ! {


      
      
#ifdef MPI
         if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: done check w'
        IF ( debug_mpi ) CALL MPI_BARRIER(my_comm, mpi_error_code)
         if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: update w'

      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = 0
      kze = ktile+1
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO k = kzb,kze
      if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: update w: k,kdivz,rhalf,coefw = ',k,kdivz(k),rhalf(k),coefw(k)
       DO j = jyb,jye
        DO i = ixb,ixe
#else
!$OMP PARALLEL DO IF ( ny .gt. 2 ), DEFAULT(SHARED), PRIVATE(i,j,k)
      DO j = 1,ny-1
       DO k = 2,nz-1
        DO i = 1,nx-1
#endif
! Compute RHS for w (explicit)


          w(i,j,k) = w(i,j,k) + fw(i,j,k)   &
                    - coefw(k) * (p(i,j,k) -  p(i,j,k-1))   &
                    + kdivz(k)*rhalf(k)*(div(i,j,k) - div(i,j,k-1))

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF ! } explicit
!-----------------------------------------------------------------------
! Compute u-equation
! Compute v-equation

#ifdef MPI
        IF ( debug_mpi ) CALL MPI_BARRIER(my_comm, mpi_error_code)
        if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: update u and v'
#endif


#ifdef MPI
      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = 1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg
      if (ixend .eq. nxend .and. bcx .eq. 2 ) ixe = ixend-ixbeg+1

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 1
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg
      if (jyend .eq. nyend .and. bcy .eq. 2 ) jye = jyend-jybeg+1

      kzb = 0
      kze = ktile+1
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb, jye ; do i = ixb, ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       DO k = 1,nz-1 
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
#ifdef MPI
           IF(ixbeg-1+i .gt. 1 .or. bcx == 2 ) THEN
#else
           IF( i .ne. 1 ) THEN
#endif
             u(i,j,k) = u(i,j,k) + fu(i,j,k)   &
                      - coefu(k)*(p(i,j,k)-p(i-1,j,k))*gx(i,4) + kdivx(i)*(div(i,j,k)-div(i-1,j,k))
           ENDIF
#ifdef MPI
           IF(jybeg-1+j .gt. 1 .or. bcy == 2 ) THEN
#else
           IF( j .ne. 1 ) THEN
#endif
             v(i,j,k) = v(i,j,k) + fv(i,j,k)   &
                      - coefv(k)*(p(i,j,k)-p(i,j-1,k))*gy(j,4) + kdivy(j)*(div(i,j,k)-div(i,j-1,k))
           ENDIF

         ENDDO
        ENDDO  
       ENDDO

       IF ( bcx .eq. 2 .and. nproci == 1 ) THEN
#ifdef MPI
      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do j = jyb, jye
#else
       DO k = 1,nz-1 
        DO j = 1,ny-1
#endif

#ifdef MPI
         if (ixbeg .eq. nxbeg) then
#endif
           i = 1

           u(i,j,k) = u(i,j,k) + fu(i,j,k)   &
                      - coefu(k)*(p(i,j,k)-p(nx-1,j,k))*gx(i,4) + kdivx(i)*(div(i,j,k)-div(nx-1,j,k))
#ifdef MPI
         endif
         
         if(ixend .eq. nxend) then
           i = nxend
#else
           i = nx
#endif
           u(i,j,k) = u(1,j,k) 
!           u(i,j,k) = u(i,j,k) + fu(i,j,k)   &
!                      - coefu(k)*(p(1,j,k)-p(i-1,j,k))*gx(1,4) + kdivx(1)*(div(1,j,k)-div(i-1,j,k))
#ifdef MPI
         endif
#endif

        ENDDO
       ENDDO
       ENDIF

       IF ( bcy .eq. 2 .and. nprocj == 1 ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze ; do i = ixb, ixe
#else
       DO k = 1,nz-1 
        DO i = 1,nx-1
#endif

#ifdef MPI
         if(jybeg .eq. nybeg) then
#endif
            j = 1

            v(i,j,k) = v(i,j,k) + fv(i,j,k)   &
                      - coefv(k)*(p(i,j,k)-p(i,ny-1,k))*gy(j,4) + kdivy(j)*(div(i,j,k)-div(i,ny-1,k))
#ifdef MPI
         endif
         
         if(jyend .eq. nyend) then
          j = nyend
#else
          j = ny
#endif
             v(i,j,k) = v(i,1,k)
!             v(i,j,k) = v(i,j,k) + fv(i,j,k)   &
!                      - coefv(k)*(p(i,1,k)-p(i,j-1,k))*gy(1,4) + kdivy(1)*(div(i,1,k)-div(i,j-1,k))
#ifdef MPI
         endif
#endif
        ENDDO
       ENDDO
       ENDIF

! Update boundaries

#ifdef MPI
      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        IF ( nx .gt. 2 .and. bcx .ne. 2 ) THEN
#ifdef MPI
        jyb = -1
        jye = jtile+2
        if (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        do j = jyb, jye
          if(ixbeg.eq.nxbeg) u(1 ,j,k) = u( 1,j,k) + fu( 1,j,k)
          if(ixend.eq.nxend) u(nx,j,k) = u(nx,j,k) + fu(nx,j,k)
        enddo
#else
        DO j = 1,ny-1
          u(1, j,k) = u( 1,j,k) + fu( 1,j,k)
          u(nx,j,k) = u(nx,j,k) + fu(nx,j,k)
        ENDDO
#endif
        ENDIF

        IF ( ny .gt. 2 .and. bcy .ne. 2  ) THEN
#ifdef MPI
        ixb = -1
        ixe = itile+2
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do i = ixb, ixe
          if(jybeg.eq.nybeg) v(i,1 ,k) = v(i,1 ,k) + fv(i,1 ,k)
          if(jyend.eq.nyend) v(i,ny,k) = v(i,ny,k) + fv(i,ny,k)
        enddo
#else
        DO i = 1,nx-1
          v(i, 1,k) = v(i, 1,k) + fv(i, 1,k)
          v(i,ny,k) = v(i,ny,k) + fv(i,ny,k)
        ENDDO
#endif
        ENDIF

       ENDDO


! MPI: Send velocity data at end of last small step

#ifdef MPI
    IF ( number_of_processes .gt. 1 ) THEN
    
    CALL cld_cpu('MPI-COMM-SMLSTEP')

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: comm u,v,w'

! u,v,w velocity -- note that ut,vt,wt are contiguous in memory

        nb = 1
        ia = 3
        IF ( n == nsmall ) nb = ng
        IF ( .not. ctest .or.  n == nsmall ) THEN
        
        IF ( nproci > 1 ) THEN

        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,u)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,u)

        ENDIF

        IF ( nprocj > 1 ) THEN

        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,u)

        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,u)

        ENDIF

        IF ( nprock > 1 ) THEN

        downward_tag = 2005
        CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),downward_tag,u)

        upward_tag = 2006
        CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),upward_tag,u)
   
        ENDIF

        ELSEIF ( .not. fucomm ) THEN  ! if n < nsmall, only need to communicate U to E-W and V to N-S.

        ia = 1

        IF ( nproci > 1 ) THEN
        westward_tag = 2001
        CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),westward_tag,u)

        eastward_tag = 2002
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             w_proc(my_rank),e_proc(my_rank),eastward_tag,u)
        ENDIF

        IF ( nprocj > 1 ) THEN
        southward_tag = 2003
        CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),southward_tag,v)
             
        northward_tag = 2004
        CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,ia,  &
             n_proc(my_rank),s_proc(my_rank),northward_tag,v)
        ENDIF

        IF ( nprock > 1 ) THEN
        downward_tag = 2005
        CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),downward_tag,w)

        upward_tag = 2006
        CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,ia,  &
           d_proc(my_rank),u_proc(my_rank),upward_tag,w)
        ENDIF

        ENDIF

    CALL cld_cpu('MPI-COMM-SMLSTEP')
    
    ENDIF

#endif

! End small step iteration loop

   if (debug_mpi) write(0,*) my_rank, 'SMLSTEP: end of nsmall loop'

      ENDDO

! Rewrite w-array and add back in V(t) values

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng

      jyb = -ng+1
      jye = jtile+ng

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
!       DO k = kzb,kze
!        write(91,*) my_rank, 'u,w,p,fpp ,div= ',k,u(nx/2+2,ny/2,k),w(nx/2,ny/2,k),p(nx/2,ny/2,k),fpp(nx/2+2,ny/2,k), &
!          div(nx/2+2,ny/2,k),fw(nx/2+2,ny/2,k),fu(nx/2+2,ny/2,k)
!       ENDDO

      do k = kzb,kze ; do j = jyb, jye ; do i = ixb, ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      DO k = 1,nz-1
       DO j = 1,ny
        DO i = 1,nx
#endif
          u(i,j,k) = u(i,j,k) + ut(i,j,k)
          v(i,j,k) = v(i,j,k) + vt(i,j,k)
          w(i,j,k) = w(i,j,k) + wt(i,j,k)
          p(i,j,k) = p(i,j,k) + pt(i,j,k)

          w(i,j,k) =   w(i,j,k) / rhalf(k)
         wt(i,j,k) =  wt(i,j,k) / rhalf(k)

        ENDDO
       ENDDO
      ENDDO


      deallocate ( fpp )
      deallocate ( div )

#ifdef MPI
   if (debug_mpi)  write(0,"('SMLSTEP: EXITING SUBROUTINE, my_rank=',1x,i2)") my_rank
#endif

      RETURN
      END SUBROUTINE SMLSTEP

