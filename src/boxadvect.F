c2345678901234567890123456789012345678901234567890123456789012345678912
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  STRAKA'S ATMOSPHERIC MODEL  (SAM)
c    Designed by Jerry M. Straka
c
c  SUBROUTINE VELADV
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c2345678901234567890123456789012345678901234567890123456789012345678912
c
c
c
      subroutine boxadvect
     : (nx,ny,nz,nor,
     :  dt,dx,dy,dz,gtz,
     :  vc,vd,
     :  fu,fv,fw )
c
c      include 'sam.impnone.h'
      USE PARAM_MODULE

      USE COMMASMPI_MODULE

      implicit none
!      include 'sam.param.h'
c
c  INTEGERS
c
      integer :: ix, jy, kz
      integer :: nx, ny, nz
      
      integer  nv, mu, mv, mw
      parameter ( nv = 3, mu = 1, mv = 2, mw = 3 )
      integer istag,jstag,kstag
      parameter ( istag = 1, jstag = 1, kstag = 1 )
      integer nor
      

      integer    id1
      parameter (id1=1)
      integer    jd1
      parameter (jd1=1)
      integer    kd1
      parameter (kd1=1)
c
c  REALS
c
      real         dt,   grdt    
      real         dx,     dy,     dz
      real        dxi,    dyi,    dzi
      real       dxis,   dyis,   dzis
      real       dxt4,   dyt4,   dzt4
      real       dxi4,   dyi4,   dzi4
      real       dxi2,   dyi2,   dzi2
      real               dyt,    dzt
      real      dxtca,  dytca,  dztca
      real      dxist,  dyist,  dzist   
      real        aml,  amlsq,  ddspt,   dsht 
      real     vnorme, vnormw, vnorms, vnormn
c
c
      real gtz(nz)

c      real, allocatable :: vc(:,:,:,:) ! vc(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv)
c      real, allocatable :: vd(:,:,:,:) ! vd(nx,ny,nz,nv)
c
c

!      real fo(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)


      real fu(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fw(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      real  vc(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv)
      real  vd(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv)

      real, allocatable :: t0(:,:,:)
      real, allocatable :: t1(:,:,:)
      real, allocatable :: t2(:,:,:)
      real, allocatable :: t3(:,:,:)
      real, allocatable :: t4(:,:,:)
      
      
      real base, ainflobx,binflo
      parameter ( base = 1.0, ainflobx = 0.0, binflo = 0.0 )
      
      integer iadvx, iadvy, iadvz
      
      integer ibc,jbc

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!
      ibc = 0
      jbc = 0
      
      IF ( bcx == 2 ) ibc = 1
      IF ( bcy == 2 ) jbc = 1
!      IF ( bcx .eq. 2 .or. bcy .eq. 2 ) THEN
!        write(0,*) 'STOP! boxadvect scheme not set up yet for periodic BCs!!'
!        STOP
!      ENDIF

     
!      allocate ( vc(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv) )
!      allocate ( vd(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv) )
      allocate ( t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )
      allocate ( t1(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )
      allocate ( t2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )
      allocate ( t3(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )
      allocate ( t4(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )


      iadvx = 0
      iadvy = 0
      iadvz = 0
      
      IF ( nx .gt. 5 ) iadvx = 1
      IF ( ny .gt. 5 ) iadvy = 1
      IF ( nz .gt. 5 ) iadvz = 1
c
c
c  constants
c
      dxt4  = -1./(4.0*dx)
      dyt4  = -1./(4.0*dy)
      dzt4  = -1./(4.0*dz)
!      dxi   =  (1.0)/dx
!      dxi2  =  (0.5)/dx
!      dxi4  =  (0.25)/dx
!      dyi   =  (1.0)/dy
!      dyi2  =  (0.5)/dy
!      dyi4  =  (0.25)/dy
!      dzi   =  (1.0)/dz
!      dzi2  =  (0.5)/dz
!      dzi4  =  (0.25)/dz
      dxi   = -1./dx
      dyt   = -1./dy
      dzt   = -1./dz
c
c
c
c  ADVECTION
c
c
c
c  x-direction advection -box method- for u velocity.
c
      if ( iadvx .ge. 1 ) then
      if ( nx .gt. 5 ) then
#ifdef MPI
      ixb = 0
      IF (ixbeg .eq. nxbeg) ixb = 1-ibc
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1,     ny-jstag
      do ix = 1-ibc,     nx-id1+ibc
#endif
      t1(ix,jy,kz) =
     >  dxt4
     >  *(vc(ix+id1,jy,kz,mu)+vc(ix,jy,kz,mu))
     <  *(vc(ix+id1,jy,kz,mu)-vc(ix,jy,kz,mu))
      enddo
      enddo
      enddo
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1,     ny-jstag
      do ix = 1+id1-ibc, nx-id1+ibc
#endif
      fu(ix,jy,kz) =  fu(ix,jy,kz) 
     > + t1(ix,jy,kz) + t1(ix-id1,jy,kz) 
      enddo
      enddo
      enddo
      end if
c
c  x-direction advection -box method- for v velocity.
c
      if ( nx .gt. 1 ) then
#ifdef MPI
      ixb = 0
      IF (ixbeg .eq. nxbeg) ixb = 1-ibc
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      jyb = 1
      jye = jtile
      if (jybeg .eq. nxbeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1+jbc
      do ix = 1-ibc,     nx-id1+ibc
#endif
      t2(ix,jy,kz) = 
     > dxt4
     > *(vc(ix+id1,jy,kz,mu)+vc(ix+id1,jy-jd1,kz,mu))
     > *(vc(ix+id1,jy,kz,mv)-vc(ix,jy,kz,mv))
      enddo
      enddo
      enddo
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag-id1+ibc

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1+jbc
      do ix = 1+id1-ibc, nx-istag-id1+ibc
#endif
      fv(ix,jy,kz) =  fv(ix,jy,kz) 
     > + t2(ix,jy,kz) + t2(ix-id1,jy,kz) 
      enddo
      enddo
      enddo
      end if
c
c  east - west boundary condition for v in x-direction
c
      if ( nx .gt. 5 .and. ( myproci == 1 .or. myproci == nproci) .and. bcx /= 2 ) then
#ifdef MPI
      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-kd1

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1
#endif
#ifdef MPI
      if ( myproci == 1) then
#endif
      vnormw = 
     >  0.25*(vd(1,jy,kz,mu)+vd(1+istag,jy,kz,mu)
     >      + vd(1,jy-jd1,kz,mu)+vd(1+istag,jy-jd1,kz,mu))
      fv(1,jy,kz) = fv(1,jy,kz) 
     > + max(vnormw*dxi, 0.0)
     >   *((vd(1+id1,jy,kz,mv) - vd(1,jy,kz,mv))
!     >   -base*(vb(1+id1,jy,kz,mv) - vb(1,jy,kz,mv))
     >   )
!     > + (ainflo*min(vnormw*dxi, 0.0) - dt*binflo)
!     >   *(vd(1,jy,kz,mv)-base*vb(1,jy,kz,mv))
#ifdef MPI
      endif

      if ( myproci == nproci ) then
#endif
      vnorme = 
     >  0.25*(vd(nx,jy,kz,mu)+vd(nx-istag,jy,kz,mu)
     >      + vd(nx,jy-jd1,kz,mu)+vd(nx-istag,jy-jd1,kz,mu))
      fv(nx-istag,jy,kz) = fv(nx-istag,jy,kz) 
     > + min(vnorme*dxi, 0.0)
     >   *((vd(nx-istag,jy,kz,mv) -vd(nx-istag-id1,jy,kz,mv))
!     >   -base*(vb(nx-istag,jy,kz,mv) -vb(nx-istag-id1,jy,kz,mv))
     >  )
!     > - (ainflo*max(vnorme*dxi, 0.0) + dt*binflo)
!     >   *(vd(nx-istag,jy,kz,mv)-base*vb(nx-istag,jy,kz,mv))
#ifdef MPI
      endif
#endif
      enddo
      enddo
      end if
      end if
c  
c  u advection -box method- of w-momentum
c
      if ( iadvx .ge. 1 ) then
      if ( nx .gt. 1 ) then
#ifdef MPI
      ixb = 0
      IF (ixbeg .eq. nxbeg) ixb = 1-ibc
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do jy = 1,     ny-jstag
      do ix = 1-ibc,     nx-istag-id1+ibc
#endif
      t3(ix,jy,kz) = 
     > dxt4
     > *(vc(ix+id1,jy,kz,mu)+vc(ix+id1,jy,kz-kd1,mu))
     > *(vc(ix+id1,jy,kz,mw)-vc(ix,jy,kz,mw))
      enddo
      enddo
      enddo
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do jy = 1,     ny-jstag
      do ix = 1+id1-ibc, nx-istag-id1+ibc
#endif
      fw(ix,jy,kz) =  fw(ix,jy,kz) 
     > + t3(ix,jy,kz) + t3(ix-id1,jy,kz) 
      enddo
      enddo
      enddo
      end if
c
c  east - west boundary condition for w in x-direction
c
      if ( nx .gt. 5 .and. ( myproci == 1 .or. myproci == nproci) .and. bcx /= 2 ) then
#ifdef MPI
      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
#else
      do jy = 1,     ny-jstag
      do kz = 1+kd1, nz-kd1
#endif
      if ( myproci == 1) then
      vnormw = 
     >   0.25*(vd(1+istag,jy,kz,mu)+vd(1+istag,jy,kz-kd1,mu)
     >  +      vd(1,jy,kz,mu)+vd(1,jy,kz-kd1,mu))
      fw(1,jy,kz) = fw(1,jy,kz) 
     > + max(vnormw*dxi, 0.0)
     >   *((vd(1+istag,jy,kz,mw) - vd(1,jy,kz,mw))
!     >   -base*(vb(1+istag,jy,kz,mw) - vb(1,jy,kz,mw))
     > )
!     > + (ainflo*min(vnormw*dxi, 0.0) - dt*binflo)
!     >   *(vd(1,jy,kz,mw)-base*vb(1,jy,kz,mw))
      ENDIF
      IF ( myproci == nproci) then
      vnorme = 
     >   0.25*(vd(nx-istag,jy,kz,mu)+vd(nx-istag,jy,kz-kd1,mu)
     >  +      vd(nx,jy,kz,mu)+vd(nx,jy,kz-kd1,mu))
      fw(nx-istag,jy,kz) = fw(nx-istag,jy,kz) 
     > + min(vnorme*dxi, 0.0)
     >   *((vd(nx-istag,jy,kz,mw) - vd(nx-istag-id1,jy,kz,mw))
!     >   -base*(vb(nx-istag,jy,kz,mw) - vb(nx-istag-id1,jy,kz,mw))
     > )
!     > - (ainflo*max(vnorme*dxi, 0.0) + dt*binflo)
!     >   *(vd(nx-istag,jy,kz,mw)-base*vb(nx-istag,jy,kz,mw))
      ENDIF
      enddo
      enddo
      end if
      end if
c
c
c  y-direction
c
c 
c
c  v advection -box method- of u momentum
c
      if ( iadvy .ge. 1 ) then
      if ( ny .gt. 5 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      jyb = 0
      jye = jtile
      IF (jybeg .eq. nybeg) jyb = 1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1-jbc,     ny-jstag
      do ix = 1+id1-ibc, nx-id1
#endif
      t1(ix,jy,kz) = 
     >  dyt4
     >  *(vc(ix,jy+jd1,kz,mv)+vc(ix-id1,jy+jd1,kz,mv))
     >  *(vc(ix,jy+jd1,kz,mu)-vc(ix,jy,kz,mu))
      enddo
      enddo
      enddo
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag-id1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jstag-jd1+jbc
      do ix = 1+id1-ibc, nx-id1
#endif
      fu(ix,jy,kz) =  fu(ix,jy,kz) 
     > + t1(ix,jy,kz) + t1(ix,jy-jd1,kz) 
      enddo
      enddo
      enddo
      end if
c
c  north - south boundary condition for u in y-direction
c
      if ( ny .gt. 5 .and. ( myprocj == 1 .or. myprocj == nprocj) .and. bcy /= 2 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do ix = 1+id1-ibc, nx-id1+ibc
#endif
      IF ( myprocj == 1) THEN
      vnorms = 
     >  0.25*(vd(ix-istag,1,kz,mv) + vd(ix-istag,1+jstag,kz,mv)
     >      + vd(ix,1,kz,mv) + vd(ix,1+jstag,kz,mv))
      fu(ix,1,kz) = fu(ix,1,kz) 
     > + max(vnorms*dyt, 0.0)
     >   *((vd(ix,1+jd1,kz,mu) - vd(ix,1,kz,mu))
!     >   -base*(vb(ix,1+jd1,kz,mu) - vb(ix,1,kz,mu))
     > )
!     > + (ainflo*min(vnorms*dyt, 0.0) - dt*binflo)
!     >   *(vd(ix,1,kz,mu)-base*vb(ix,1,kz,mu))
      ENDIF
      IF ( myprocj == nprocj) THEN
      vnormn = 
     >  0.25*(vd(ix-istag,ny-jstag,kz,mv) + vd(ix-istag,ny,kz,mv)
     >       +vd(ix,ny-jstag,kz,mv) + vd(ix,ny,kz,mv))
      fu(ix,ny-jstag,kz) = fu(ix,ny-jstag,kz) 
     >  + min(vnormn*dyt, 0.0)
     >   *((vd(ix,ny-jstag,kz,mu) - vd(ix,ny-jstag-jd1,kz,mu))
!     >   -base*(vb(ix,ny-jstag,kz,mu) - vb(ix,ny-jstag-jd1,kz,mu))
     >  )
!     > - (ainflo*max(vnormn*dyt, 0.0) + dt*binflo)
!     >   *(vd(ix,ny-jstag,kz,mu)-base*vb(ix,ny-jstag,kz,mu))
      ENDIF
      
      enddo
      enddo
      end if
      end if
c
c  v advection -box method- of v-component
c    
      if ( iadvy .ge. 1 ) then
      if ( ny .gt. 5 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 0
      jye = jtile
      IF (jybeg .eq. nybeg) jyb = 1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1-jbc,     ny-jstag+jbc
      do ix = 1,     nx-istag
#endif
      t2(ix,jy,kz) = 
     >  dyt4
     >  *(vc(ix,jy+jd1,kz,mv)+vc(ix,jy,kz,mv))
     >  *(vc(ix,jy+jd1,kz,mv)-vc(ix,jy,kz,mv))      
      enddo
      enddo
      enddo
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1, ny-jd1
      do ix = 1,     nx-istag
#endif
      fv(ix,jy,kz) =  fv(ix,jy,kz)
     > + t2(ix,jy,kz) + t2(ix,jy-jd1,kz)
      enddo
      enddo
      enddo
      end if
      end if
c
c  v advection -box-method- of w wind
c
      if ( iadvy .ge. 1 ) then
      if ( ny .gt. 5 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 0
      jye = jtile
      IF (jybeg .eq. nybeg) jyb = 1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do jy = 1-jbc, ny-jstag
      do ix = 1,     nx-istag
#endif
      t3(ix,jy,kz) = 
     >   dyt4
     >   *(vc(ix,jy+jd1,kz,mv)+vc(ix,jy+jd1,kz-kd1,mv))
     >   *(vc(ix,jy+jd1,kz,mw)-vc(ix,jy,kz,mw))
      enddo
      enddo
      enddo
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do jy = 1+jd1-jbc, ny-jstag-jd1+jbc
      do ix = 1,     nx-istag
#endif
      fw(ix,jy,kz) =  fw(ix,jy,kz) 
     > + t3(ix,jy,kz) + t3(ix,jy-jd1,kz) 
      enddo
      enddo
      enddo
      end if
c
c  north - south boundary condition for w in y-direction
c
      if ( ny .gt. 5 .and. ( myprocj == 1 .or. myprocj == nprocj) .and. bcy /= 2 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do ix = 1,     nx-istag
#endif
      IF ( myprocj == 1 ) THEN
      vnorms = 
     <   0.25*(vd(ix,1,kz,mv)+vd(ix,1+jstag,kz,mv)
     <       + vd(ix,1,kz-kstag,mv)+vd(ix,1+jstag,kz-kstag,mv))
      fw(ix,1,kz) = fw(ix,1,kz) 
     < + max(vnorms*dyt, 0.0)
     <   *((vd(ix,1+jd1,kz,mw) - vd(ix,1,kz,mw))
!     <   -base*(vb(ix,1+jd1,kz,mw) - vb(ix,1,kz,mw))
     < )
!     < + (ainflo*min(vnorms*dyt, 0.0) - dt*binflo)
!     >   *(vd(ix,1,kz,mw)-base*vb(ix,1,kz,mw))
      ENDIF
      IF ( myprocj == nprocj ) THEN
      vnormn =
     >  0.25*(vd(ix,ny,kz,mv)+vd(ix,ny-jstag,kz,mv)
     <      + vd(ix,ny,kz-kstag,mv)+vd(ix,ny-jstag,kz-kstag,mv))
      fw(ix,ny-jstag,kz) = fw(ix,ny-jstag,kz) 
     > + min(vnormn*dyt, 0.0)
     >   *((vd(ix,ny-jstag,kz,mw) - vd(ix,ny-jstag-jd1,kz,mw))
!     >   -base*(vb(ix,ny-jstag,kz,mw) - vb(ix,ny-jstag-jd1,kz,mw))
     > )
!     > - (ainflo*max(vnormn*dyt, 0.0) + dt*binflo)
!     >   *(vd(ix,ny-jstag,kz,mw)-base*vb(ix,ny-jstag,kz,mw))
      ENDIF
      
      enddo
      enddo
      endif
      endif
c
c
c
c
c  z-direction
c
c  set fo if itopo ge 1 
c
c      if ( itopo .ge. 1 ) then
c      call setfo(nx,ny,nz,dx,dy,dz,ht,gt,dc,vc)
c      end if
c
c
c  w advection -box method- of u momentum
c
c 
c  no topography (itopo .eq. 0 ) just these terms need be computed
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1,     ny-jstag
      do ix = 1+id1-ibc, nx-id1+ibc
#endif
      t3(ix,jy,kz) =
     >  (0.50)*(vc(ix,jy,kz+kd1,mw)+vc(ix-id1,jy,kz+kd1,mw))
      t4(ix,jy,kz) =  
     >  (0.50)
     > *(gtz(kz+kd1)+gtz(kz))
     > *(vc(ix,jy,kz+kd1,mu)-vc(ix,jy,kz,mu))
      enddo
      enddo
      enddo
c
c  extra terms for topography (itopo .eq. 1 ) 
c  redefine t3 = average of fo 
c
c      if ( itopo .ge. 1 ) then
cc
c      do 71111 kz = 1,     nz-kstag
c      do 71112 jy = 1,     ny-jstag
c      do 71113 ix = 1+id1, nx-id1
c      t3(ix,jy,kz) = 
c     <  (0.50)*(fo(ix,jy,kz+kd1)+fo(ix-id1,jy,kz+kd1))
c71113 continue
c71112 continue
c71111 continue
cc
c      end if
c
c  compute wdu/dz (t0)
c
      if ( iadvz .ge. 1 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1,     ny-jstag
      do ix = 1+id1, nx-id1
#endif
      t0(ix,jy,kz) = t3(ix,jy,kz) * t4(ix,jy,kz)
      enddo
      enddo
      enddo
c
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      IF ( myprock == 1 ) THEN
      do jy = jyb,jye
      do ix = ixb,ixe
      t0(ix,jy,00) = 0.0  
      enddo
      enddo
      ENDIF

      IF ( myprock == nprock ) THEN
      do jy = jyb,jye
      do ix = ixb,ixe
      t0(ix,jy,nz-kstag) = 0.0
      enddo
      enddo
      ENDIF
      
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1,     ny-jstag
      do ix = 1+id1-ibc, nx-id1+ibc
#endif
      fu(ix,jy,kz) =  fu(ix,jy,kz)
     > + (0.50)*dzt*(t0(ix,jy,kz) + t0(ix,jy,kz-kd1))
      enddo
      enddo
      enddo
      end if
c
c  w advection -box method- of v momentum
c
c
c  no topography (itopo .eq. 0 ) just these terms need be computed
c
      if ( iadvz .ge. 1 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1+jbc
      do ix = 1,     nx-istag
#endif
      t3(ix,jy,kz) =
     >  (0.50)*(vc(ix,jy,kz+kd1,mw)+vc(ix,jy-jd1,kz+kd1,mw))
      t4(ix,jy,kz) = 
     >  (0.50)
     > *(gtz(kz+kd1)+gtz(kz))
     > *(vc(ix,jy,kz+kd1,mv)-vc(ix,jy,kz,mv))
      enddo
      enddo
      enddo
c
c  extra terms for topography (itopo .eq. 1 )
c  redefine t3 = average of fo
c
c      if ( itopo .eq. 1 ) then
cc
c      do 72111 kz = 1,     nz-kstag
c      do 72112 jy = 1+jd1, ny-jd1
c      do 72113 ix = 1,     nx-istag
c      t3(ix,jy,kz) = 
c     <  (0.50)*(fo(ix,jy,kz+kd1)+fo(ix,jy-jd1,kz+kd1))
c72113 continue
c72112 continue
c72111 continue
cc
c      end if
      end if
c
c  compute wdv/dz (t0)
c
      if ( iadvz .ge. 1 ) then
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1+jbc
      do ix = 1,     nx-istag
#endif
      t0(ix,jy,kz) = t3(ix,jy,kz) * t4(ix,jy,kz)
      enddo
      enddo
      enddo
c
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      IF ( myprock == 1 ) THEN
      do jy = jyb,jye
      do ix = ixb,ixe
      t0(ix,jy,00) = 0.0  
      enddo
      enddo
      ENDIF

      IF ( myprock == nprock ) THEN
      do jy = jyb,jye
      do ix = ixb,ixe
      t0(ix,jy,nz-kstag) = 0.0
      enddo
      enddo
      ENDIF

c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1+jd1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kstag
      do jy = 1+jd1-jbc, ny-jd1+jbc
      do ix = 1,     nx-istag
#endif
      fv(ix,jy,kz) =  fv(ix,jy,kz)
     > + (0.50)*dzt*(t0(ix,jy,kz) + t0(ix,jy,kz-kd1))
      enddo
      enddo
      enddo
      end if
c
c  w advection -box method- of w momentum
c
      if ( iadvz .ge. 1 ) then
c
c  no topography (itopo .eq. 0 )
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kd1
      do jy = 1,     ny-jstag
      do ix = 1,     nx-istag
#endif
      t3(ix,jy,kz) =
     >  (0.50)*(vc(ix,jy,kz+kd1,mw)+vc(ix,jy,kz,mw))
      t4(ix,jy,kz) = 
     >  (gtz(kz))
     > *(vc(ix,jy,kz+kd1,mw)-vc(ix,jy,kz,mw))
      enddo
      enddo
      enddo
c
c  extra terms for topography (itopo .eq. 1 )
c  redefine t3 = average of fo  
c
c      if ( itopo .eq. 1 ) then
cc
c      do 73111 kz = 1,     nz-kd1
c      do 73112 jy = 1,     ny-jstag
c      do 73113 ix = 1,     nx-istag
c      t3(ix,jy,kz) = 
c     <  (0.50)*(fo(ix,jy,kz+kd1)+fo(ix,jy,kz))
c73113 continue
c73112 continue
c73111 continue
cc
c      end if
c
c  compute wdw/dz (t0)
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 0
      kze = ktile
      if (kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1,     nz-kd1
      do jy = 1,     ny-jstag
      do ix = 1,     nx-istag
#endif
      t0(ix,jy,kz) = t3(ix,jy,kz) * t4(ix,jy,kz)
      enddo
      enddo
      enddo
c
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1+kd1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
#else
      do kz = 1+kd1, nz-kd1
      do jy = 1,     ny-jstag
      do ix = 1,     nx-istag
#endif
      fw(ix,jy,kz) =  fw(ix,jy,kz)
     > + (0.50)*dzt*(t0(ix,jy,kz) + t0(ix,jy,kz-kd1))
      enddo
      enddo
      enddo
c
      end if
c
c
c  end of routine
c
c
!      deallocate ( vc )
!      deallocate ( vd )
      deallocate ( t0 )
      deallocate ( t1 )
      deallocate ( t2 )
      deallocate ( t3 )
      deallocate ( t4 )
      
      return
      end
c
c
c
