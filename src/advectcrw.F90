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
      subroutine ADVECTCRW           &
       (nx,ny,nz,ns,nst,             &
        dt,dx,dy,dz, gtx,gty,gtz,   &
        ab,a,an,dc, u,v,w,           &
        t0, fu,fv,fw,iadiv,ia,izero,ad, svar,jpass) 

      USE GRID_MODULE
      USE PARAM_MODULE
      USE COMMASMPI_MODULE

      implicit none

   TYPE(VARIABLE) :: svar

!  declare all variables
!
!  INTEGERS
!
      integer    ia,ns,nba
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
      integer    nz
      parameter (np=6)
      integer    i,j,k,n
      
      integer  iadiv
      integer  nst
      integer  izero
      integer  jpass
      integer  ipass0
!
!  REALS
!
      real       qu, qd, qc, qr, qdel, qxs, qxp, qsr, qsl
      real       aqdel, qcurv, aqcurv
      
      real   ba1, ba2
!
      real       dt
      real       dx
      real       dy
      real       dz

      real       dtdx
      real       dtdy
      real       dtdz

      real       dtdx4
      real       dtdy4
      real       dtdz4
      real       dtfac
      real       vnorme
      real       vnormw
      real       vnormn
      real       vnorms

      real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng), gtz(-ng+1:nz+ng)
!
!      real       cx(nxl,nzl,np),  
!      real  px(nx,np,np),py(ny,np,np),pz(nz,np,np)
!      real  px(0:nx-1,np,np),py(0:ny-1,np,np),pz(0:nz-1,np,np)
      real  px(-ng+1:nx+ng,np,np),py(-ng+1:ny+ng,np,np),pz(-ng+1:nz+ng,np,np)
!      real       cy(nyl,nzl,np),  py(nyl,np,np)
!      real       cz(nxl,nzl,np),  pz(nzl,np,np)
!
!
!
      real       ab(-ng+1:nz+ng)
!
      real       dc(-ng+1:nz+ng)

      real       a(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)  ! s
      real       ad(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real       an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns) ! st
!
      real u(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real v(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real w(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real fus,fvs,fws

      real t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
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
!    icrwmp = 1;   at step 1, do x-pass, y-pass, z-pass
!                  at step 2, do y-pass, x-pass, z-pass
!    icrwmp = 2;   at step 1, do x-pass, y-pass, z-pass, z-pass, y-pass, x-pass
!                  at step 2, do y-pass, x-pass, z-pass, z-pass, x-pass, y-pass
!                    each part uses only 1/2dt for the time step giving a
!                    cumulative of dt; e.g., two half passes for each direction
! 
!  set type of crowley implementation ( sam.param.h )
!
!      IF ( bcx .eq. 2 .or. bcy .eq. 2 ) THEN
!        write(0,*) 'STOP! Crowley scheme not set up yet for periodic BCs!!'
!        STOP
!      ENDIF

      if ( icrwmp .eq. 1 ) then
       ipass0 = 1
       dtfac = 1.0
       npass = 3
       mdnstp = mod(nst,2)
      end if
      if ( icrwmp .eq. 2 ) then
       dtfac = 0.5
       npass = 3*jpass
       ipass0 = 1 + 3*(jpass-1) ! 1 for jpass=1 or 4 for jpass=2
       mdnstp = mod(nst,2)
      end if
      
      
      ainflotmp = 0.0
      
      IF (  ainflo .gt. 0.0 ) THEN

         relaxscalar = .true.

   DO k = 1,nz-1
      relaxscalar = (relaxscalar .and. svar%base1d(k) == 0.0)
   ENDDO
  
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
      dtdz  = -dtfac*dt/dz
      dtdx4 =  (0.25)*dtfac*dt/dx
      dtdy4 =  (0.25)*dtfac*dt/dy
      dtdz4 =  (0.25)*dtfac*dt/dz
      
!      IF ( ia .ge. 3 .and. izero .ne. 0 ) THEN
      IF ( izero .ne. 0 ) THEN
#ifdef MPI
        imn = nxend-istag
        imx = 1
        jmn = nyend-jstag
        jmx = 1
        kmn = nzend-kstag
        kmx = 1
#else
        imn = nx-istag
        imx = 1
        jmn = ny-jstag
        jmx = 1
        kmn = nz-kstag
        kmx = 1
#endif
!      ENDDO

#ifdef MPI
       ixb=-ng+1
       if(ixbeg == 1 ) ixb = 1
       ixe=itile+ng-1
       if(ixend.eq.nxend) ixe=ixend-ixbeg-istag
       
       jyb=-ng+1
       if ( jybeg == 1 ) jyb = 1
       jye=jtile+ng-1
       if(jyend.eq.nyend) jye=jyend-jybeg-jstag

       kzb=1
       kze=ktile
       if(kzend.eq.nzend) kze=kzend-kzbeg-kstag
       
!       write(0,*) 'check zeros for ia = ',ia,my_rank,ixbeg,ixb,ixe
       
       do k=kzb,kze ; do j=jyb,jye ; do i=ixb,ixe
          IF ( an(i,j,k,ia) .ne. 0.0 ) THEN
            imn = Min(imn,Max(1,ixbeg-1+i-3))
            imx = Max(imx,Min(ixbeg-1+i+3,nxend-1))
            jmn = Min(jmn,Max(1,jybeg-1+j-3))
            jmx = Max(jmx,Min(jybeg-1+j+3,nyend-1))
            kmn = Min(kmn,Max(1,kzbeg-1+k-3))
            kmx = Max(kmx,Min(kzbeg-1+k+3,nzend-1))
          ENDIF
       enddo ; enddo ; enddo

        IF ( imn <= ixbeg ) imn = 1
        IF ( imx >= ixend ) imx = nxend - 1
!        imn = 1
!        imx = nxend-istag
!        jmn = 1
!        jmx = nyend-jstag
        kmn = 1
        kmx = nzend-kstag
!       write(0,*) 'imn, imx, jmn, jmx for ia = ',ia,my_rank,imn,imx,jmn,jmx
#else
       DO k=1,nz-kstag
        DO j=1,ny-jstag
         DO i=1,nx-istag
          IF ( an(i,j,k,ia) .ne. 0.0 ) THEN
            imn = Min(imn,Max(1,i-3))
            imx = Max(imx,Min(i+3,nx-1))
            jmn = Min(jmn,Max(1,j-3))
            jmx = Max(jmx,Min(j+3,ny-1))
            kmn = Min(kmn,Max(1,k-3))
            kmx = Max(kmx,Min(k+3,nz-1))
          ENDIF
         ENDDO
        ENDDO
       ENDDO
#endif
      ELSE
#ifdef MPI
        imn = 1
        imx = nxend-istag
        jmn = 1
        jmx = nyend-jstag
        kmn = 1
        kmx = nzend-kstag
#else
        imn = 1
        imx = nx-istag
        jmn = 1
        jmx = ny-jstag
        kmn = 1
        kmx = nz-kstag
#endif
      ENDIF
       
       IF ( kmn .gt. kmx ) THEN 
         RETURN
       ENDIF
       
       IF ( .false. .and. nst .le. 2 ) THEN
        do i=1,nx
         print*,'ix,gtx = ',i,gtx(i)
        ENDDO
        do j=1,ny
         print*,'jy,gty = ',j,gty(j)
        ENDDO
        do k=1,nz
         print*,'kz,gtz = ',k,gtz(k)
        ENDDO
        print*, 'dtdx, etc: ',dtdx,dtdz,dx,dy,dz
        STOP
       ENDIF

#ifdef MPI
       ixb=-ng+1
       ixe=itile+ng
!       if(ixbeg.eq.nxbdg) ixe=1
!       if(ixend.eq.nxend) ixe=ixend-ixbeg
       
       jyb=-ng+1
       jye=jtile+ng
!       if(jyend.eq.nyend) jye=jyend-jybeg

       kzb=-ng+1
       kze=ktile+ng
!       if(kzend.eq.nzend) kze=kzend-kzbeg
       
       do k=kzb,kze ; do j=jyb,jye ; do i=ixb,ixe
         ad(i,j,k) = a(i,j,k)
       enddo ; enddo ; enddo
#else
       DO k = 1,nz-1
         ad(1:nx-1,1:ny-1,k) = a(1:nx-1,1:ny-1,k)
       ENDDO
#endif

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
!        write(0,*) 'call crwcff for z'
      call crwcff(np,nz,nz,ng,kstag,kd1,ktile,kzbeg,kzend,nzbeg,nzend,pz)
!
!
!
!  SIXTH ORDER CROWLEY FLUX-ANTI FLUX SCHEME
!
!
!
      do ipass = ipass0,npass
!
!  x-direction crowley corrected flux computations
!
      if ( ipass .eq. (1+mdnstp) .or. ipass .eq. (6-mdnstp) ) then 
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

       kzb=1
       kze=ktile
       if(kzend.eq.nzend) kze=kzend-kzbeg
       
       do n=1,ng
       do k=kzb,kze ; do j=jyb,jye ; do i=ixb,ixe
        an( nxend-1+n,j,k,ia) = an(n   ,j,k,ia)
        an( 1-n   ,j,k,ia) = an(nxend-n,j,k,ia)
        ad( nxend-1+n,j,k) = ad(n   ,j,k)
        ad( 1-n   ,j,k) = ad(nxend-n,j,k)
       enddo ; enddo ; enddo
       enddo
     ENDIF
#else
      DO n = 1,ng
        an( nx-1+n,1:ny,1:nz,ia) = an(n   ,1:ny,1:nz,ia)
        an( 1-n   ,1:ny,1:nz,ia) = an(nx-n,1:ny,1:nz,ia)
        ad( nx-1+n,1:ny,1:nz) = ad(n   ,1:ny,1:nz)
        ad( 1-n   ,1:ny,1:nz) = ad(nx-n,1:ny,1:nz)
       ENDDO
#endif
      ENDIF

      
      t0(:,:,:) = 0.0
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
  if (ixbeg .le. imn) ixb = Max(0,imn-ixbeg+1-ibc)
  if (ixend .ge. imx) ixe = Min(itile,imx-ixbeg+1-id1+ibc)

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
  if (jybeg .le. jmn) jyb = jmn-jybeg+1
  if (jyend .ge. jmx) jye = jmx-jybeg+1

  kzb = -ng+1
  kze = ktile+ng
  if (kzbeg .le. kmn) kzb = kmn
  if (kzend .ge. kmx) kze = kmx-kzbeg+1

  do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
  
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i,  &
!$OMP   fus,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)  
      do k = kmn,kmx ! 1,     nz-kstag
      do j = jmn,jmx ! 1,     ny-jstag
!      do 2023 i = Max(1,imn-1),Min(nx-istag-id1,imx)  ! 1,     nx-istag-id1
      do i = imn-ibc, imx-id1+ibc ! 1,     nx-istag-id1
#endif
      fus = fu(i,j,k)

      IF ( Abs(fus) .lt. eps ) THEN
       t0(i,j,k) = 0.0
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
        t1s*an(i-id2,j,k,ia)+   &
        t2s*an(i-id1,j,k,ia)+   &
        t3s*an(i    ,j,k,ia)+   &
        t4s*an(i+id1,j,k,ia)+   &
        t5s*an(i+id2,j,k,ia)+   &
        t6s*an(i+id3,j,k,ia)
      IF ( icrwmn .eq. 1  ) THEN
      t0s = -t0s / (fus)
      qsr  = max(sign((1.0),fus), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(i-id1*1,j,k,ia)   &
           + qsl*an(i+id1*2,j,k,ia)
      qd   = qsr*an(i+id1*1,j,k,ia)   &
           + qsl*an(i,      j,k,ia)
      qc   = qsr*an(i,      j,k,ia)   &
           + qsl*an(i+id1*1,j,k,ia)
      qr   = qu + (qc - qu) / (abs(fus))
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),   &
                 abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
!      IF ( DBG ) print*, 'i,j,k,imn,imx,ia = ',i,j,k,imn,imx,ia
      t0(i,j,k) =   &
        -fus *    &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
       t0(i,j,k) = t0s
      ENDIF
      
      enddo
      enddo
      enddo

!  monotonic face values
!
!      if ( icrwmn .eq. 1 ) then
      if ( .false. ) then
#ifdef MPI
       ixb = 0
       ixe = ng-1

       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg-jstag

       kzb = 1
       kze = ktile
       if (kzend .eq. nxend) kze = kzend-kzbeg-kstag

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
        an(-i,j,k,ia)   = an(i+1,j,k,ia)
        an(nxend+i,j,k,ia) = an(nxend-istag-i,j,k,ia)
       enddo ; enddo ; enddo
#else
       do k = 1,     nz-kstag
       do j = 1,     ny-jstag
       do i = 0,     ng-1 
        an(-i,j,k,ia)   = an(i+1,j,k,ia)
        an(nx+i,j,k,ia) = an(nx-istag-i,j,k,ia)
       enddo
       enddo
       enddo
#endif

#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixbeg .lt. imn) ixb = Max(1,imn-ixbeg+1)
       if (ixend .gt. imx) ixe = min(imx,ixend-ixbeg+1-istag-id1)

       jyb = 1
       jye = jtile
       if (jybeg .lt. jmn) jyb = jmn-jybeg+1
       if (jyend .gt. jmx) jye = jmx-jybeg+1

       kzb = 1
       kze = ktile
       if (kzbeg .lt. kmn) kzb = kmn-kzbeg+1
       if (kzend .gt. kmx) kze = kmx-kzbeg+1

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       do k = kmn,kmx !  1, nz-kstag
       do j = jmn,jmx !  1, ny-jstag
       do i = imn,Min(nx-istag-id1,imx)  ! 1, nx-istag-id1
#endif
        t0(i,j,k) = -t0(i,j,k) / (fu(i,j,k) + 1.0e-20)
        qsr  = max(sign((1.0),fu(i,j,k)), (0.0))
        qsl  = (1.0) - qsr
        qu   = qsr*an(i-id1*1,j,k,ia)   &
             + qsl*an(i+id1*2,j,k,ia)
        qd   = qsr*an(i+id1*1,j,k,ia)   &
             + qsl*an(i,      j,k,ia)
        qc   = qsr*an(i,      j,k,ia)   &
             + qsl*an(i+id1*1,j,k,ia)
        qr   = qu + (qc - qu) / (abs(fu(i,j,k)) + 1.0e-20)
        qdel = qd - qu
        qxs  = max(sign((1.0), qdel), (0.0))
        qxp  = max(sign((1.0),   &
                   abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
        t0(i,j,k) =   &
          -fu(i,j,k) *    &
         (   &
           ((1.0)-qxp) * qc   &
         +        qxp   &
         * (      qxs  * min(max(t0(i,j,k), qc), min(qr, qd))   &
         + ((1.0)-qxs) * max(min(t0(i,j,k), qc), max(qr, qd)) ) )
       enddo
       enddo
       enddo
      endif
!       
      IF ( iadiv .eq. 1 .or. iadiv .eq. 2  ) THEN
#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1+id1-ibc)
       if (ixend .ge. imx) ixe = Min(itile,imx-ixbeg+1-id1+ibc)

       IF ( ipass .eq. 1 ) THEN
       jyb = -ng+1
       jye = jtile+ng
       ELSE
       jyb = 1
       jye = jtile
       ENDIF
       if (jybeg .le. jmn) jyb = jmn-jybeg+1
       if (jyend .ge. jmx) jye = jmx-jybeg+1

       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .le. kmn) kzb = kmn
       if (kzend .ge. kmx) kze = kmx-kzbeg+1

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
       do k = kmn , kmx ! 1, nz-kstag
       do j = jmn , jmx ! 1, ny-jstag
       do i = imn+id1-ibc, imx-id1+ibc ! Max(2,imn), Min(nx-istag-id1,imx) ! 1+id1, nx-id1-istag
#endif
        an(i,j,k,ia) = an(i,j,k,ia)     &
          + (t0(i,j,k)-t0(i-id1,j,k))   &
          * (gtx(i)**2)                      &
          + (ad(i,j,k))                    &
          * (fu(i,j,k)-fu(i-id1,j,k))   &
          * (gtx(i)**2)
       enddo
       enddo
       enddo

      ELSE

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .lt. imn) ixb = max(2,imn-ixbeg+1)
      if (ixend .gt. imx) ixe = min(ixend-ixbeg+1-istag-id1,imx-ixbeg+1)

      IF ( ipass .eq. 1 ) THEN
       jyb = -ng+1
       jye = jtile+ng
      ELSE
       jyb = 1
       jye = jtile
      ENDIF
      if (jybeg .lt. jmn) jyb = jmn-jybeg+1
      if (jyend .gt. jmx) jye = jmx-jybeg+1

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .lt. kmn) kzb = kmn
      if (kzend .gt. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do  k = kmn,kmx ! 1,     nz-kstag
      do  j = jmn,jmx ! 1,     ny-jstag
      do i = imn+id1-ibc, imx-id1+ibc ! Max(2,imn), Min(nx-istag-id1,imx) ! 1+id1, nx-id1-istag
#endif
      an(i,j,k,ia) = an(i,j,k,ia)    &
        + (t0(i,j,k)-t0(i-id1,j,k))   &
        * (gtx(i)**2)
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
!
! 2020 continue
!
!  x-direction boundaries
!
#ifdef MPI

      IF ( bcx .ne. 2 .and. ( ixbeg .eq. nxbeg .or. ixend .eq. nxend ) ) THEN

       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

       kzb = 1
       kze = ktile
       if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

       do k = kzb,kze ; do j = jyb,jye
#else

      IF ( bcx .ne. 2 .and. ( imn .eq. 1 .or. imx .eq. nx-istag ) ) THEN

! !$omp  PARALLEL DO DEFAULT(SHARED), &
! !$omp  PRIVATE(vnormw,vnorme,k,j,ba1,ba2,ia),SHARED(vc,an,ad) 
      do k = 1,nz-kstag
      do j = 1,ny-jstag

#endif
        ba1 = ab(k)
        ba2 = ab(k)
#ifdef MPI
       IF ( ixbeg .eq. nxbeg ) THEN
#else
       IF ( imn .eq. 1 ) THEN
#endif
       vnormw =    &
         0.5*(u(1,j,k) +u(1+istag,j,k))
       an(1,j,k,ia) = an(1,j,k,ia)    &
        + max(vnormw*dtdx, 0.0)*gtx(1)   &
          *(ad(1+id1,j,k) - ad(1,j,k))   &
        + (ainflotmp*min(vnormw*dtdx, 0.0) - dt*dtfac*binflo)   &
          *(ad(1,j,k) - ba1)
       ENDIF
#ifdef MPI
       IF ( ixend .eq. nxend ) THEN
#else
       IF ( imx .eq. nx-istag ) THEN
#endif
       vnorme =    &
         0.5*(u(nx,j,k) + u(nx-istag,j,k))
       an(nx-istag,j,k,ia) = an(nx-istag,j,k,ia)    &
        + min(vnorme*dtdx, 0.0)*gtx(nx-istag)   &
          *(ad(nx-istag,j,k) - ad(nx-istag-id1,j,k))   &
        - (ainflotmp*max(vnorme*dtdx, 0.0) + dt*dtfac*binflo)   &
          *(ad(nx-istag,j,k) - ba2)
       ENDIF

      enddo
      enddo
      ENDIF !( bcy .ne. 2 .and. ( imn .eq. 1 .or. imx .eq. nx-istag ) )
      end if ! nx.gt.5
!
      end if ! ipass
!
!  y-direction crowley corrected flux computation
!
      if ( ipass .eq. (2-mdnstp) .or. ipass .eq. (5+mdnstp) ) then
!
      if ( ny .gt. 5 ) then
!
! C$DOACROSS LOCAL(k,j,i)
!
! C$DOACROSS LOCAL(k,j,i)
! !$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i) 
      do k = 1, nz-kstag
      do j = 1, ny-jstag-jd1
      do i = 1, nx-istag
      enddo
      enddo
      enddo
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

    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

     do n = 1,ng
     do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
      an(i,nyend-1+n,k,ia) = an(i,n   ,k,ia)
      an(i,1-n   ,k,ia)    = an(i,nyend-n,k,ia)
      ad(i,nyend-1+n,k)    = ad(i,n   ,k)
      ad(i,1-n   ,k)       = ad(i,nyend-n,k)
     enddo ; enddo ; enddo
     enddo
    ENDIF
#else
   DO n = 1,ng
    an(1:nx-1,ny-1+n,1:nz-1,ia) = an(1:nx-1,n   ,1:nz-1,ia)
    an(1:nx-1,1-n   ,1:nz-1,ia) = an(1:nx-1,ny-n,1:nz-1,ia)
    ad(1:nx-1,ny-1+n,1:nz-1) = ad(1:nx-1,n   ,1:nz-1)
    ad(1:nx-1,1-n   ,1:nz-1) = ad(1:nx-1,ny-n,1:nz-1)
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

      t0(:,:,:) = 0.0

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
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .ge. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 0
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn-jybeg+1-jbc
      if (jyend .ge. jmx) jye = jmx-jybeg+1-jd1+jbc

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .le. kmn) kzb = kmn
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(k,j,i,     &
!$OMP   fvs,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)
      do k = kmn,kmx ! 1, nz-kstag
      do j = jmn-jbc, jmx-jd1+jbc ! Max(1,jmn-1),Min(ny-jstag-jd1,jmx)  ! 1, ny-jstag-jd1
      do i = imn,imx ! 1, nx-istag
#endif
      fvs = fv(i,j,k)

      IF ( Abs(fvs) .lt. eps ) THEN
       t0(i,j,k) = 0.0
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
        t1s*an(i,j-jd2,k,ia) +   &
        t2s*an(i,j-jd1,k,ia) +   &
        t3s*an(i,j    ,k,ia) +   &
        t4s*an(i,j+jd1,k,ia) +   &
        t5s*an(i,j+jd2,k,ia) +   &
        t6s*an(i,j+jd3,k,ia)
      IF ( icrwmn .eq. 1 ) THEN
      t0s = -t0s / (fvs)
      qsr  = max(sign((1.0),fvs), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(i,j-jd1*1,k,ia)   &
           + qsl*an(i,j+jd1*2,k,ia)
      qd   = qsr*an(i,j+jd1*1,k,ia)   &
           + qsl*an(i,j      ,k,ia)
      qc   = qsr*an(i,j      ,k,ia)   &
           + qsl*an(i,j+jd1*1,k,ia)
      qr   = qu + (qc - qu) / (abs(fvs))
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),   &
                 abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
      t0(i,j,k) =   &
       -fvs *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
        t0(i,j,k) = t0s
      ENDIF

     ENDDO
     ENDDO
     ENDDO
!
!  monotonic face values
!
!      if ( icrwmn .eq. 1 ) then 
      if ( .false. ) then 

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      jyb = 0
      jye = ng-1

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
       an(i,-j,k,ia)   = an(i,j+1,k,ia)
       an(i,nyend+j,k,ia) = an(i,nyend-jstag-j,k,ia)
      enddo ; enddo ; enddo
#else    
      do k = 1, nz-kstag
      do j = 0, ng-1
      do i = 1, nx-istag
       an(i,-j,k,ia)   = an(i,j+1,k,ia)
       an(i,ny+j,k,ia) = an(i,ny-jstag-j,k,ia)
      enddo
      enddo
      enddo
#endif

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .lt. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .gt. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .lt. jmn) jyb = jmn-jybeg+1
      if (jyend .gt. jmx) jye = min(nyend-jstag-jd1,jmx-jybeg+1)

      kzb = 1
      kze = ktile
      if (kzbeg .lt. kmn) kzb = kmn-kzbeg+1
      if (kzend .gt. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
      do k = kmn,kmx !  1, nz-kstag
      do j = jmn,Min(ny-jstag-jd1,jmx)  ! 1, ny-jstag-jd1
      do i = imn,imx ! 1, nx-istag
#endif
       t0(i,j,k) = -t0(i,j,k) / (fv(i,j,k) + 1.0e-20)
       qsr  = max(sign((1.0),fv(i,j,k)), (0.0))
       qsl  = (1.0) - qsr
       qu   = qsr*an(i,j-jd1*1,k,ia)   &
            + qsl*an(i,j+jd1*2,k,ia)
       qd   = qsr*an(i,j+jd1*1,k,ia)   &
            + qsl*an(i,j      ,k,ia)
       qc   = qsr*an(i,j      ,k,ia)   &
            + qsl*an(i,j+jd1*1,k,ia)
       qr   = qu + (qc - qu) / (abs(fv(i,j,k)) + 1.0e-20)
       qdel = qd - qu
       qxs  = max(sign((1.0), qdel), (0.0))
       qxp  = max(sign((1.0),   &
                  abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
       t0(i,j,k) =   &
        -fv(i,j,k) *   &
        (   &
          ((1.0)-qxp) * qc   &
        +        qxp   &
        * (      qxs  * min(max(t0(i,j,k), qc), min(qr, qd))   &
        + ((1.0)-qxs) * max(min(t0(i,j,k), qc), max(qr, qd)) ) )
      enddo
      enddo
      enddo
      end if
!
      IF ( iadiv .eq. 1 .or. iadiv .eq. 2 ) THEN


#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .ge. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn-jybeg+1+jd1-jbc
      if (jyend .ge. jmx) jye = jmx-jybeg+1-jd1+jbc

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do k = kmn,kmx                  ! 1, nz-kstag
      do j = jmn+jd1-jbc, jmx-jd1+jbc ! Max(2,jmn), Min(ny-jstag-jd1,jmx) ! 1+jd1, ny-jstag-jd1
      do i = imn,imx                  ! 1, nx-istag
#endif
       an(i,j,k,ia) = an(i,j,k,ia)      &
         + (t0(i,j,k)-t0(i,j-jd1,k))    &
         * (gty(j)**2)                  &
         + (ad(i,j,k))                  &
         * (fv(i,j,k)-fv(i,j-jd1,k))    &
         * (gty(j)**2)
      enddo
      enddo
      enddo

!      ix = nx/2
!      kz = 16
!      print*,'nstep,limits = ',nst,jmn+id1-jbc, jmx-id1+jbc
!      DO jy = -1,1
!        print*,'t0 ',jy,t0(ix,jy+1,kz),t0(ix,ny+jy-2,kz)
!      ENDDO
!      DO jy = 0,2
!        print*,'an ',jy,an(i,j,k,ia),an(ix,ny-2+jy,kz,ia)
!        print*,'ab ',jy,ad(i,j,k),ad(ix,ny-2+jy,kz)
!      ENDDO

      ELSE

#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .ge. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = max(2,jmn-jybeg+1)
      if (jyend .ge. jmx) jye = min(nyend-jstag-jd1,jmx-jybeg+1)

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
      do  k = kmn,kmx                           ! 1, nz-kstag
      do  j = Max(2,jmn), Min(ny-jstag-jd1,jmx) ! 1+jd1, ny-jstag-jd1
      do  i = imn,imx                           ! 1, nx-istag
#endif
       an(i,j,k,ia) = an(i,j,k,ia)      &
         + (t0(i,j,k)-t0(i,j-jd1,k))    &
         * (gty(j)**2)
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF

!
! 3020 continue
!
!  y-direction boundaries
!  scalars on south and north boundaries
!

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

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      do k = kzb,kze
      do i = ixb,ixe
#else
      do k = 1, nz-kstag
      do i = 1, nx-istag
#endif
        ba1 = ab(k)
        ba2 = ab(k)

#ifdef MPI
      IF ( jybeg .eq. nybeg ) THEN
#else
      IF ( jmn .eq. 1 ) THEN
#endif
      vnorms =    &
        0.5*(v(i,1+jstag,k)+ v(i,1,k))
      an(i,1,k,ia) = an(i,1,k,ia)    &
       + max(vnorms*dtdy, 0.0)*gty(1)   &
         *(ad(i,1+jd1,k) - ad(i,1,k))   &
       + (ainflotmp*min(vnorms*dtdy, 0.0) - dt*dtfac*binflo)   &
         *(ad(i,1,k) - ba1)
      ENDIF

#ifdef MPI
      IF ( jyend .eq. nyend ) THEN
#else
      IF ( jmx .eq. ny-jstag ) THEN
#endif
      vnormn =    &
        0.5*(v(i,ny-jstag,k) + v(i,ny,k))
      an(i,ny-jstag,k,ia) = an(i,ny-jstag,k,ia)    &
       + min(vnormn*dtdy, 0.0)*gty(ny-jstag)   &
         *(ad(i,ny-jstag,k) - ad(i,ny-jstag-jd1,k))   &
       - (ainflotmp*max(vnormn*dtdy, 0.0) + dt*dtfac*binflo)   &
         *(ad(i,ny-jstag,k) - ba2)
      ENDIF
     enddo
     enddo


      ENDIF  !( bcy .ne. 2 .and. ( jmn .eq. 1 .or. jmx .eq. nyend-jstag) )
!     enddo
      end if
!
      end if
!
!  z-direction crowley corrected flux computation
!
!     if ( ipass .eq. (3+mdnstp) .or. ipass .eq. (4-mdnstp) ) then
      if ( ipass .eq. 3 .or. ipass .eq. 4 ) then
!
!
      if ( nz .gt. 5 ) then
!
!  compute polynomials
!
! C$DOACROSS LOCAL(k,j,i)
!!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)
!      do 4011 kz = 1,     nz-kstag-kd1
!      do 4012 jy = 1,     ny-jstag
!      do 4013 ix = 1,     nx-istag
! 4013 continue
! 4012 continue
! 4011 continue
!
! C$DOACROSS LOCAL(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,i,j,k,ia), 
! C$&       share(t1,t2,t3,t4,t5,t6,an,ad,fw,dc,gt)
!!$OMP  PARALLEL DO DEFAULT(SHARED), &
!!$omp PRIVATE(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,i,j,k,ia) 
!! c$omp+ ,shared(t1,t2,t3,t4,t5,t6,an,ad,fw,dc,gt)
!      do ia = 1,     iadva
!
!      if ( icrwmn .eq. 1 ) then 
!      do  kz = 0,     ng-1
!      do  jy = 1,     ny-jstag
!      do  ix = 1,     nx-istag
!      an(ix,jy,-kz,ia)   = an(i,j,k+1,ia)
!      an(ix,jy,nz+kz,ia) = an(ix,jy,nz-kstag-kz,ia)
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDIF

      t0(:,:,:) = 0.0

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .ge. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn-jybeg+1
      if (jyend .ge. jmx) jye = jmx-jybeg+1

      kzb = 0
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn-kzbeg+1
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k,   &
!$OMP   fws,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp) 
      do k = kmn,kmx ! Max(1,kmn-1),Min(kmx,nz-kstag-kd1) ! 1, nz-kstag-kd1
      do j = jmn,jmx ! 1, ny-jstag
      do i = imn,imx ! 1, nx-istag
#endif
      fws = fw(i,j,k)

      IF ( Abs(fws) .lt. eps ) THEN
       t0(i,j,k) = 0.0
       CYCLE
      ENDIF

      t1s =    &
        fws*(pz(k,1,1) +   &
        fws*(pz(k,1,2) +   &
        fws*(pz(k,1,3) +   &
        fws*(pz(k,1,4) +   &
        fws*(pz(k,1,5) +   &
        fws*(pz(k,1,6)))))))
      t2s =    &
        fws*(pz(k,2,1) +   &
        fws*(pz(k,2,2) +   &
        fws*(pz(k,2,3) +   &
        fws*(pz(k,2,4) +   &
        fws*(pz(k,2,5) +   &
        fws*(pz(k,2,6)))))))
      t3s =    &
        fws*(pz(k,3,1) +   &
        fws*(pz(k,3,2) +   &
        fws*(pz(k,3,3) +   &
        fws*(pz(k,3,4) +   &
        fws*(pz(k,3,5) +   &
        fws*(pz(k,3,6)))))))
      t4s =    &
        fws*(pz(k,4,1) +   &
        fws*(pz(k,4,2) +   &
        fws*(pz(k,4,3) +   &
        fws*(pz(k,4,4) +   &
        fws*(pz(k,4,5) +   &
        fws*(pz(k,4,6)))))))
      t5s =    &
        fws*(pz(k,5,1) +   &
        fws*(pz(k,5,2) +   &
        fws*(pz(k,5,3) +   &
        fws*(pz(k,5,4) +   &
        fws*(pz(k,5,5) +   &
        fws*(pz(k,5,6)))))))
      t6s =    &
        fws*(pz(k,6,1) +   &
        fws*(pz(k,6,2) +   &
        fws*(pz(k,6,3) +   &
        fws*(pz(k,6,4) +   &
        fws*(pz(k,6,5) +   &
        fws*(pz(k,6,6)))))))
      t0s =    &
        t1s*an(i,j,k-kd2,ia)+   &
        t2s*an(i,j,k-kd1,ia)+   &
        t3s*an(i,j,k    ,ia)+   &
        t4s*an(i,j,k+kd1,ia)+   &
        t5s*an(i,j,k+kd2,ia)+   &
        t6s*an(i,j,k+kd3,ia)


      IF ( icrwmn .eq. 1 .and. k .gt. 1 ) THEN
      t0s = -t0s / (fws)
      qsr  = max(sign((1.0),fws), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(i,j,k-kd1*1,ia)   &
           + qsl*an(i,j,k+kd1*2,ia)
      qd   = qsr*an(i,j,k+kd1*1,ia)   &
           + qsl*an(i,j,k      ,ia)
      qc   = qsr*an(i,j,k      ,ia)   &
           + qsl*an(i,j,k+kd1*1,ia)
      qr   = qu + (qc - qu) / (abs(fws))
      qdel = qd - qu
      qxs  = max(sign((1.0), qdel), (0.0))
      qxp  = max(sign((1.0),   &
                 abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
      t0(i,j,k) =   &
       -fws *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      ELSE
       t0(i,j,k) = t0s
      ENDIF

!      IF ( DBG .and. ia .eq. 2 .and. fws .gt. 0.01 )  &
!            print*, 'i,j,k,ia,fws,t0 = ',i,j,k,ia,fws,w(i,j,k),t0(i,j,k),an(i,j,k,ia)

     enddo
     enddo
     enddo
!
!  monotonic face values
!
      if ( .false. ) then 
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg-jstag

      kzb = 0
      kze = ng-1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
       an(i,j,-k,ia)   = an(i,j,k+1,ia)
       an(i,j,nzend+k,ia) = an(i,j,nzend-kstag-k,ia)
      enddo ; enddo ; enddo
#else
      do k = 0,     ng-1
      do j = 1,     ny-jstag
      do i = 1,     nx-istag
       an(i,j,-k,ia)   = an(i,j,k+1,ia)
       an(i,j,nz+k,ia) = an(i,j,nz-kstag-k,ia)
      enddo
      enddo
      enddo
#endif

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .lt. imn) ixb = imn-ixbeg+1
      if (ixend .gt. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile
      if (jybeg .lt. jmn) jyb = jmn-jybeg+1
      if (jyend .gt. jmx) jye = jmx-jybeg+1

      kzb = 1
      kze = ktile
      if (kzbeg .lt. kmn) kzb = kmn-kzbeg+1
      if (kzend .gt. kmx) kze = kmx-kzbeg+1 ! min(kmx-kzbeg+1,nzend-kstag-kd1)

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
      do k = kmn,Min(kmx,nz-kstag-kd1) ! 1, nz-kstag-kd1
      do j = jmn,jmx                   ! 1, ny-jstag
      do i = imn,imx                   ! 1, nx-istag
#endif
       t0(i,j,k) = -t0(i,j,k) / (fw(i,j,k) + 1.0e-20)
       qsr  = max(sign((1.0),fw(i,j,k)), (0.0))
       qsl  = (1.0) - qsr
       qu   = qsr*an(i,j,k-kd1*1,ia)   &
            + qsl*an(i,j,k+kd1*2,ia)
       qd   = qsr*an(i,j,k+kd1*1,ia)   &
            + qsl*an(i,j,k      ,ia)
       qc   = qsr*an(i,j,k      ,ia)   &
            + qsl*an(i,j,k+kd1*1,ia)
       qr   = qu + (qc - qu) / (abs(fw(i,j,k)) + 1.0e-20)
       qdel = qd - qu
       qxs  = max(sign((1.0), qdel), (0.0))
       qxp  = max(sign((1.0),   &
                  abs(qdel)-abs(qu-(2.0)*qc+qd)), (0.0))
       t0(i,j,k) =   &
        -fw(i,j,k) *   &
        (   &
          ((1.0)-qxp) * qc   &
        +        qxp   &
        * (      qxs  * min(max(t0(i,j,k), qc), min(qr, qd))   &
        + ((1.0)-qxs) * max(min(t0(i,j,k), qc), max(qr, qd)) ) )
      enddo
      enddo
      enddo
      end if

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg-jstag

      do j = jyb,jye ; do i = ixb,ixe
       IF ( myprock == 1 ) THEN
       t0(i,j, 0) = 0.0
       fw(i,j, 0) = 0.0
       ENDIF
       IF ( myprock == nprock ) THEN
       t0(i,j,nz-kstag) = 0.0
       fw(i,j,nz-kstag) = 0.0
       ENDIF
      enddo ; enddo
#else
      do j = 1, ny-jstag
      do i = 1, nx-istag
       t0(i,j, 0) = 0.0
       t0(i,j,nz-kstag) = 0.0
       fw(i,j, 0) = 0.0
       fw(i,j,nz-kstag) = 0.0
      enddo
      enddo
#endif
      IF ( iadiv .eq. 1 ) THEN

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .le. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .ge. imx) ixe = Min(ixe,imx-ixbeg+1)
!      if (ixend .eq. nxend) ixe = Min(ixe,ixend-ixbeg-istag)

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = Max(jyb,jmn-jybeg+1)
      if (jyend .ge. jmx) jye = Min(jye,jmx-jybeg+1)

      kzb = 1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do k = kmn,kmx ! 1, nz-kstag
      do j = jmn,jmx ! 1, ny-jstag
      do i = imn,imx ! 1, nx-istag
#endif
!     an(i,j,k,ia) = an(i,j,k,ia) 
!    >  + ( (t0(i,j,k)-t0(i,j,k-kd1))
!    >     +(ad(i,j,k,ia))
!    >     *(fw(i,j,k)-fw(i,j,k-kd1)) ) 
!    >  / (dc(i,j,k))
      an(i,j,k,ia) = an(i,j,k,ia)    &
        + ( (0.5)                    &
          *(t0(i,j,k)                &
          *(dc(k+kd1)+dc(k))         &
          -(dc(k-kd1)+dc(k))         &
           *t0(i,j,k-kd1))           &
           +(ad(i,j,k))*0.5          &
           *(fw(i,j,k)               &
          *(dc(k+kd1)+dc(k))         &
          -(dc(k-kd1)+dc(k))         &
            *fw(i,j,k-kd1)) )        &
        * (gtz(k))                   &
        / (dc(k))
      enddo
      enddo
      enddo

      ELSEIF (  iadiv .eq. 2 ) THEN

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .lt. imn) ixb = Max(ixb,imn-ixbeg+1)
      if (ixend .gt. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .lt. jmn) jyb = Max(jyb,jmn-ixbeg+1)
      if (jyend .gt. jmx) jye = Min(jye,jmx-jybeg+1)

      kzb = 1
      kze = ktile
      if (kzbeg .lt. kmn) kzb = kmn
      if (kzend .gt. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do  k = kmn,kmx ! 1, nz-kstag
      do  j = jmn,jmx ! 1, ny-jstag
      do  i = imn,imx ! 1, nx-istag
#endif
      an(i,j,k,ia) = an(i,j,k,ia)       &
        + (t0(i,j,k) - t0(i,j,k-kd1))   &
        * (gtz(k))                      &
         + (ad(i,j,k))                  &
        * (fw(i,j,k)-fw(i,j,k-kd1))     &
        * (gtz(k))

      
      ENDDO
      ENDDO
      ENDDO

      ELSE

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .lt. imn) ixb = Max(1,imn-ixbeg+1)
      if (ixend .gt. imx) ixe = Min(ixe,imx-ixbeg+1)

      jyb = 1
      jye = jtile
      if (jybeg .lt. jmn) jyb = Max(1,jmn-jybeg+1)
      if (jyend .gt. jmx) jye = Min(jye,jmx-jybeg+1)

      kzb = 1
      kze = ktile
      if (kzbeg .lt. kmn) kzb = kmn
      if (kzend .gt. kmx) kze = kmx-kzbeg+1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(i,j,k) 
      do  k = kmn,kmx ! 1, nz-kstag
      do  j = jmn,jmx ! 1, ny-jstag
      do  i = imn,imx ! 1, nx-istag
#endif
!      an(i,j,k,ia) = an(i,j,k,ia) 
!     >  + (t0(i,j,k) - t0(i,j,k-kd1))
!     >  * (gtz(k))
      an(i,j,k,ia) = an(i,j,k,ia)    &
        + ( (0.5)   &
          *(t0(i,j,k)   &
          *(dc(k+kd1)+dc(k))   &
          -(dc(k-kd1)+dc(k))   &
           *t0(i,j,k-kd1))  )    &
        * (gtz(k))   &
        / (dc(k))
      
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
!  end do ! ia
!
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
!
!
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!
!
!
!
!  subroutine crwcff
!
      subroutine crwcff(norder,nap,nm,ng,nstag,id1,  &
                       itile,ixbeg,ixend,nxbeg,nxend,p)
!


      implicit none

      integer      n1, n2, ia, nstag
      integer      nap, nm, nd1, norder, ng
      integer      itile, ixbeg, ixend, nxbeg, nxend
      
      integer ixb, ixe, ix

      integer    id1
!      parameter (id1=1)

!      real p(0:nap-1,norder,norder)
      real p(-ng+1:nap+ng,norder,norder)
!
!  zero coeficients
!      
      p(:,:,:) = 0.0

!      do n2 = 1,norder
!      do n1 = 1,norder
!      do ia = 0,nm-1
!       p(ia,n1,n2) = 0.0
!      enddo
!      enddo
!      enddo
!
!  second order coeficients (4 values)
!

      ixb = -ng+1
      if ( ixbeg .eq. nxbeg ) ixb = 1
      ixe = itile+ng
      if ( ixend .eq. nxend ) ixe = ixend-ixbeg-id1

      do ix = ixb, ixe

      if ( ((ixbeg.eq.nxbeg) .and. (ix.eq.1)) .or.             &
           ((ixend.eq.nxend) .and. (ix.eq.ixend-ixbeg-id1)) ) then

      p(ix,3,1) = -(1.0/2.0)
      p(ix,4,1) = -(1.0/2.0)

      p(ix,3,2) = -(1.0/2.0)
      p(ix,4,2) =  (1.0/2.0)
 
!      end if
!
!  fourth order coefficients (16 values)
!
      elseif ( ((ixbeg.eq.nxbeg) .and. (ix.eq.id1*2))  .or.        &
           ((ixend.eq.nxend) .and. (ix.eq.ixend-ixbeg-id1*2)) ) then

      p(ix,2,1) =  (1.0/12.0)
      p(ix,3,1) = -(7.0/12.0)
      p(ix,4,1) = -(7.0/12.0)
      p(ix,5,1) =  (1.0/12.0)

      p(ix,2,2) =  (1.0/24.0)
      p(ix,3,2) = -(15.0/24.0)
      p(ix,4,2) =  (15.0/24.0)
      p(ix,5,2) = -(1.0/24.0)

      p(ix,2,3) = -(1.0/12.0)
      p(ix,3,3) =  (1.0/12.0)
      p(ix,4,3) =  (1.0/12.0)
      p(ix,5,3) = -(1.0/12.0)

      p(ix,2,4) = -(1.0/24.0)
      p(ix,3,4) =  (3.0/24.0)
      p(ix,4,4) = -(3.0/24.0)
      p(ix,5,4) =  (1.0/24.0)

!      end if
!
!  sixth order coefficients (36 values)
!
!      if ( ((ixend .ne. nxend) .and. (ix .ge. id1*3))  .or.             &
!           ((ixbeg .ne. nxbeg) .and. (ix .le. ixend-ixbeg-id1*3)) .or.  &
!           ((ixbeg .ne. nxbeg) .and. (ixend .ne. nxend)) ) then

      else

      p(ix,1,1) = -(1.0/60.0)
      p(ix,2,1) =  (8.0/60.0)
      p(ix,3,1) = -(37.0/60.0)
      p(ix,4,1) = -(37.0/60.0)
      p(ix,5,1) =  (8.0/60.0)
      p(ix,6,1) = -(1.0/60.0)

      p(ix,1,2) = -(2.0/360.0)
      p(ix,2,2) =  (25.0/360.0)
      p(ix,3,2) = -(245.0/360.0)
      p(ix,4,2) =  (245.0/360.0)
      p(ix,5,2) = -(25.0/360.0)
      p(ix,6,2) =  (2.0/360.0)

      p(ix,1,3) =  (1.0/48.0)
      p(ix,2,3) = -(7.0/48.0)
      p(ix,3,3) =  (6.0/48.0)
      p(ix,4,3) =  (6.0/48.0)
      p(ix,5,3) = -(7.0/48.0)
      p(ix,6,3) =  (1.0/48.0)

      p(ix,1,4) =  (1.0/144.0)
      p(ix,2,4) = -(11.0/144.0)
      p(ix,3,4) =  (28.0/144.0)
      p(ix,4,4) = -(28.0/144.0)
      p(ix,5,4) =  (11.0/144.0)
      p(ix,6,4) = -(1.0/144.0)

      p(ix,1,5) = -(1.0/240.0)
      p(ix,2,5) =  (3.0/240.0)
      p(ix,3,5) = -(2.0/240.0)
      p(ix,4,5) = -(2.0/240.0)
      p(ix,5,5) =  (3.0/240.0)
      p(ix,6,5) = -(1.0/240.0)

      p(ix,1,6) = -(1.0/720.0)
      p(ix,2,6) =  (5.0/720.0)
      p(ix,3,6) = -(10.0/720.0)
      p(ix,4,6) =  (10.0/720.0)
      p(ix,5,6) = -(5.0/720.0)
      p(ix,6,6) =  (1.0/720.0)

      end if

      end do

!
!  return
!
      return
      end
!
!  end of subroutine
!
!
!

!c
!

! #############################################################
!
!  subroutine crwcff1 (for periodic BC)
!
      subroutine crwcff1(norder,nap,nm,ng,nstag,nd1,p)
!
      implicit none
!
      integer      n1,     n2,     ia,    nstag
      integer      nap,    nm,    nd1, norder
      integer      ng
      real p(-ng+1:nap+ng,norder,norder)
!
!  zero coeficients
!      
      p(:,:,:) = 0.0

!      do 1000 n2 = 1,norder
!      do 1001 n1 = 1,norder
!      do 1002 ia = 0,nm-1
!      p(ia,n1,n2) = 0.0
! 1002 continue
! 1001 continue
! 1000 continue
!
!  second order coeficients (4 values)
!
!      do 2000 ia = 1, nm-nstag-nd1, nm-nstag-nd1-1
!!
!      p(ia,3,1) = -(1.0/2.0)
!      p(ia,4,1) = -(1.0/2.0)
!!
!      p(ia,3,2) = -(1.0/2.0)
!      p(ia,4,2) =  (1.0/2.0)
!
! 2000 continue
!
!  fourth order coefficients (16 values)
!
!      do 4000 ia = 2,nm-nstag-nd1-1,nm-nstag-nd1-3
!!
!      p(ia,2,1) =  (1.0/12.0)
!      p(ia,3,1) = -(7.0/12.0)
!      p(ia,4,1) = -(7.0/12.0)
!      p(ia,5,1) =  (1.0/12.0)
!!
!      p(ia,2,2) =  (1.0/24.0)
!      p(ia,3,2) = -(15.0/24.0)
!      p(ia,4,2) =  (15.0/24.0)
!      p(ia,5,2) = -(1.0/24.0)
!!
!      p(ia,2,3) = -(1.0/12.0)
!      p(ia,3,3) =  (1.0/12.0)
!      p(ia,4,3) =  (1.0/12.0)
!      p(ia,5,3) = -(1.0/12.0)
!!
!      p(ia,2,4) = -(1.0/24.0)
!      p(ia,3,4) =  (3.0/24.0)
!      p(ia,4,4) = -(3.0/24.0)
!      p(ia,5,4) =  (1.0/24.0)
!!
! 4000 continue
!
!  sixth order coefficients (36 values)
!
      do 6000 ia = 0,nap-1 ! 1,nm-nstag-nd1 !  3,nm-nstag-nd1-2
!
      p(ia,1,1) = -(1.0/60.0)
      p(ia,2,1) =  (8.0/60.0)
      p(ia,3,1) = -(37.0/60.0)
      p(ia,4,1) = -(37.0/60.0)
      p(ia,5,1) =  (8.0/60.0)
      p(ia,6,1) = -(1.0/60.0)
!
      p(ia,1,2) = -(2.0/360.0)
      p(ia,2,2) =  (25.0/360.0)
      p(ia,3,2) = -(245.0/360.0)
      p(ia,4,2) =  (245.0/360.0)
      p(ia,5,2) = -(25.0/360.0)
      p(ia,6,2) =  (2.0/360.0)
!
      p(ia,1,3) =  (1.0/48.0)
      p(ia,2,3) = -(7.0/48.0)
      p(ia,3,3) =  (6.0/48.0)
      p(ia,4,3) =  (6.0/48.0)
      p(ia,5,3) = -(7.0/48.0)
      p(ia,6,3) =  (1.0/48.0)
!
      p(ia,1,4) =  (1.0/144.0)
      p(ia,2,4) = -(11.0/144.0)
      p(ia,3,4) =  (28.0/144.0)
      p(ia,4,4) = -(28.0/144.0)
      p(ia,5,4) =  (11.0/144.0)
      p(ia,6,4) = -(1.0/144.0)
!
      p(ia,1,5) = -(1.0/240.0)
      p(ia,2,5) =  (3.0/240.0)
      p(ia,3,5) = -(2.0/240.0)
      p(ia,4,5) = -(2.0/240.0)
      p(ia,5,5) =  (3.0/240.0)
      p(ia,6,5) = -(1.0/240.0)
!
      p(ia,1,6) = -(1.0/720.0)
      p(ia,2,6) =  (5.0/720.0)
      p(ia,3,6) = -(10.0/720.0)
      p(ia,4,6) =  (10.0/720.0)
      p(ia,5,6) = -(5.0/720.0)
      p(ia,6,6) =  (1.0/720.0)
!
 6000 continue
!
!  return
!
      return
      end
!
!  end of subroutine

!
!
!

! #############################################################
       subroutine SETFUVW   &
          (nx,ny,nz,dx,dy,dz,dt,u,ut,v,vt,w,wt, gtx,gty,gtz,fu,fv,fw)

! ################################################################

      USE PARAM_MODULE
#ifdef MPI
      USE COMMASMPI_MODULE
#endif

      implicit none

      integer    i,j,k, ix,jy,kz
      integer    nx,ny,nz
      
      real dx,dy,dz,dt

      real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      
      real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng), gtz(-ng+1:nz+ng)


      real :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real ut (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real vt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real wt (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

      real       cintx(-ng+1:nx+ng,0:5)
      real       cinty(-ng+1:ny+ng,0:5)
      real       cintz(-ng+1:nz+ng,0:5)

      integer istag,jstag,kstag
      parameter ( istag = 1, jstag = 1, kstag = 1 )

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

      real       dtdx4
      real       dtdy4
      real       dtdz4
      real       dtfac
      
!      integer   icrwmp
      
      logical debug
!      parameter ( debug = .true. )
      parameter ( debug = .false. )
      
      real wmax
      integer imx,jmx,kmx
      integer n
      
      integer ibc,jbc

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
      integer ixb, jyb, kzb
      integer ixe, jye, kze
#endif

! ################################################################

!      icrwmp = 2


      if ( icrwmp .eq. 1 ) then
       dtfac = 1.0
       elseif ( icrwmp .eq. 2 ) then
      dtfac = 0.5
      endif

      ibc = 0
      jbc = 0
      IF ( bcx .eq. 2 ) ibc = 1
      IF ( bcy .eq. 2 ) jbc = 1



      dtdx4 =  (0.25)*dtfac*dt/dx
      dtdy4 =  (0.25)*dtfac*dt/dy
      dtdz4 =  (0.25)*dtfac*dt/dz

!      DO ix=1,nx-1
!       print*,'i,gtx : ',ix,gtx(ix)
!      ENDDO
!      print*,'dx,dtdx4,nx = ',dx,dtdx4,nx

!      STOP


!
!
!  SET UP INTERPOLATION (MIDPOINT) COEFFICIENTS (FROM PURSER AND LESLIE 1987)
!
      do ix = 1,nx,nx-1
      cintx(ix,3) = 0.0
      cintx(ix,2) = 0.0
      cintx(ix,1) = 0.5
      cintx(ix,0) = 0.0
      end do
!
      do ix =  1+id1,nx-id1,nx-id2-id1
      cintx(ix,3) = 0.0
      cintx(ix,2) = -1.0/16.0
      cintx(ix,1) = 9.0/16.0
      cintx(ix,0) = 0.0
      end do
!
      do ix = 1+id2,nx-id2
      cintx(ix,3) = 3.0/256.0
      cintx(ix,2) = -25.0/256.0
      cintx(ix,1) = 75.0/128.0
      cintx(ix,0) = 0.0
      end do
!
! y-direction
!
      do jy = 1,ny,ny-1
      cinty(jy,3) = 0.0
      cinty(jy,2) = 0.0
      cinty(jy,1) = 0.5
      cinty(jy,0) = 0.0
      end do
!
      do jy =  1+id1,ny-id1,ny-id2-id1
      cinty(jy,3) = 0.0
      cinty(jy,2) = -1.0/16.0
      cinty(jy,1) = 9.0/16.0
      cinty(jy,0) = 0.0
      end do
!
      do jy = 1+id2,ny-id2
      cinty(jy,3) = 3.0/256.0
      cinty(jy,2) = -25.0/256.0
      cinty(jy,1) = 75.0/128.0
      cinty(jy,0) = 0.0
      end do
!
!  z-direction
!
      do kz = 1,nz,nz-id1
      cintz(kz,3) = 0.0
      cintz(kz,2) = 0.0
      cintz(kz,1) = 0.5
      cintz(kz,0) = 0.0
      end do
!
      do kz =  1+id1,nz-id1,nz-id3
      cintz(kz,3) = 0.0
      cintz(kz,2) = -1.0/16.0
      cintz(kz,1) = 9.0/16.0
      cintz(kz,0) = 0.0
      end do
!
      do kz = 1+id2,nz-id2
      cintz(kz,3) = 3.0/256.0
      cintz(kz,2) = -25.0/256.0
      cintz(kz,1) = 150.0/256.0
      cintz(kz,0) = 0.0
      end do
!
!


      if ( nx .gt. 5 ) then

      IF ( debug ) print*,'calculate fu'

#ifdef MPI
      ixb = 0
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
      DO k = 1, nz-kstag
      DO j = 1, ny-jstag
      DO i = 0, nx-istag !-id1
#endif
      fu(i,j,k) =                            &
        dtdx4                                &
       *(ut(i+id1,j,k)+ut(i+istag,j,k)       &
        + u(i+id1,j,k)+ u(i+istag,j,k))      &
       /((0.50)*(gtx(Max(1,i)) + gtx(i+id1)))
      
      ENDDO
      ENDDO
      ENDDO

      ENDIF
!
      if ( ny .gt. 5 ) then
!
#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 0
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 0
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
      DO k = 1, nz-kstag
      DO j = 0, ny-jstag !-jd1
      DO i = 1, nx-istag
#endif
      fv(i,j,k) =    &
       dtdy4   &
       *(vt(i,j+jd1,k)+vt(i,j+jstag,k)   &
        + v(i,j+jd1,k)+ v(i,j+jstag,k))   &
       /((0.50)*(gty(Max(1,j))+gty(j+jd1)))
      ENDDO
      ENDDO
      ENDDO
      
!      IF ( bcy .eq. 2 ) THEN
!        DO n = 1,2
!         fv(1:nx,ny-1+n ,1:nz)   = fv(1:nx,n-1,1:nz) 
!         fv(1:nx, 0-n ,1:nz)   = fv(1:nx,ny-n,1:nz) 
!        ENDDO
!      ENDIF

      ENDIF
!
!
      if ( nz .gt. 5 ) then

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 0
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
      DO k = 1, nz-kstag-kd1
      DO j = 1, ny-jstag
      DO i = 1, nx-istag
#endif
      fw(i,j,k) = dtdz4   &
       * ( (wt(i,j,k+kd1)+wt(i,j,k+kstag))   &
          +( w(i,j,k+kd1)+ w(i,j,k+kstag)) )
      ENDDO
      ENDDO
      ENDDO

!      if ( itopo .ge. 1 ) then
!c
!      call setfo
!     > (nx,ny,nz,dx,dy,dz,ht,gt,dc,vc,
!     :  fo,
!     :  t0,t1,t2,t3,t4,t5,t6,t7,t8,t9)
!c
!      DO 29121 k = 1,     nz-kstag
!      DO 29122 j = 1,     ny-jstag
!      DO 29123 i = 1,     nx-istag
!      t7(i,j,k) = fo(i,j,k)
!29123 continue
!29122 continue
!29121 continue
!c
!      call setfo
!     > (nx,ny,nz,dx,dy,dz,ht,gt,dc,vn,
!     :  fo,
!     :  t0,t1,t2,t3,t4,t5,t6,t7,t8,t9)
!c
!      DO 29131 k = 1,     nz-kstag
!      DO 29132 j = 1,     ny-jstag
!      DO 29133 i = 1,     nx-istag
!      fw(i,j,k) =
!     > dtdz4*
!     > ((t7(i,j,k+kd1)+t7(i,j,k+kstag))
!     < +(fo(i,j,k+kd1)+fo(i,j,k+kstag)))
!29133 continue
!29132 continue
!29131 continue
!
!      end if

      ENDIF
      
      IF ( debug ) THEN
      wmax = 0.0
      imx = 1
      jmx = 1
      kmx = 1
      DO k = 1, nz-kstag-kd1
      DO j = 1, ny-jstag
      DO i = 1, nx-istag
       IF ( w(i,j,k) .gt. wmax ) THEN
       wmax = w(i,j,k)
       imx = i
       jmx = j
       kmx = k
       ENDIF
      ENDDO
      ENDDO
      ENDDO
      print*,'imx,jmx,kmx =',imx,jmx,kmx
      DO k = 1,nz-1
       IF ( fw(imx,jmx,k) .gt. 0.01 ) THEN
       print*,'k,fu,fv,fw',k,fu(imx,jmx,k)/dtdx4/4., u(imx+1,jmx,k),   &
            fv(imx,jmx,k)/dtdy4/4.,fw(imx,jmx,k)/dtdz4/4.
       print*,'k,fu,fv,fw',k,fu(imx-1,jmx,k)/dtdx4/4., u(imx,jmx,k)
            
       STOP
       ENDIF
      ENDDO
      
      ENDIF
      
      RETURN
      END
!

