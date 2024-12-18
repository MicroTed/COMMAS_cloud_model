#define KUSFCX 0
#define KVSFCX -1
!-----------------------------------------------------------------------------
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE COMPUTETKE    ////////////////////
!     
!-----------------------------------------------------------------------------
!       
! Created by LJW: 11-01-99
! Latest update:  05-06-04
! tke_type == 2 'km' is SqrtE
!-----------------------------------------------------------------------------
!
! km  => TKE:  The RK3 integration variable for Km
!     => SMG:  The array where Km is stored
!
! kt  => TKE:  The 'n' time level of km, used in RK3 TKE integration
!     => SMG:  Not used.
!
! ft  => TKE:  Holds mixing and advection time tendencies for km.
!     => SMG:  Not used.
!
!-----------------------------------------------------------------------------
      SUBROUTINE COMPUTETKE(km,kt,shear,buoy,ft,u,v,w,kmt,kht,tke_diss,khh,khv,kmh,kmv,     &
                            uinit,vinit,pinit,tinit,kmbaserm,               &
                            gx,gy,gz,dt,ugrid,vgrid,nx,ny,nz,ns)

      USE PARAM_MODULE, only : ng, bcx, bcy, bcz, Cm, Ce, Pr, len_type, mix_type, &
     &                         zPBLhgt, drag, ahighk, Lmax, rmbasekm, g, kmgen_thresh,kmbasefac
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile, my_rank
      USE FORCE_MODULE, only : tt23sl, tt13sl, td13sl, td23sl, twt3sl
      USE MICRO_MODULE, only : iturbenhance

      implicit none

! Passed variables

      integer, INTENT(IN)    :: nx, ny, nz, ns
      real,    INTENT(INOUT) :: km    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: kt    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: kmt   (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: kht   (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: ft    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: u     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: v     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: w     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: shear (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: buoy  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: uinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: vinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: pinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: tinit(-ng+1:nz+ng,2)
      real,    INTENT(INOUT) :: kmbaserm(-ng+1:nz+ng)
      real,    INTENT(IN)    :: gx(-ng+1:nx+ng,4)
      real,    INTENT(IN)    :: gy(-ng+1:ny+ng,4)
      real,    INTENT(IN)    :: gz(-ng+1:nz+ng,4)
      real,    INTENT(IN)    :: ugrid, vgrid, dt
      real,    INTENT(INOUT) :: tke_diss(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: khh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: khv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: kmh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: kmv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

! Local variables
       
      integer i, j, k, kPBL, im1, ip1, jm1, jp1, km1, kp1
      real dx, dy, dz
      real dte, tm, rc, moist, temp, qvs, aterm
      real hgt, dz0, term1d, bb, bsh, ush, vsh 
!      double precision :: mlen       
      real :: mlen       
      real ustar, kmPBL, kmSFC, zPBL, zSFC, dzPBL, z
      real vsq0, angle0, ubase_flux, vbase_flux
!      double precision :: div
      real :: div
      real dxl,dyl,dzl

      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng), zbnd(-ng+1:nz+ng)

      real vsq(-ng+1:nx+ng,-ng+1:ny+ng), angle(-ng+1:nx+ng,-ng+1:ny+ng), uavg, vavg

      real kmbase(-ng+1:nz+ng)
      real kmbasesh(-ng+1:nz+ng)
      real kmbasebb(-ng+1:nz+ng)
      real prz(-ng+1:nz+ng)
      real prztmp(-ng+1:nz+ng)

      integer, parameter :: kusfc = KUSFCX
      integer, parameter :: kvsfc = KVSFCX
      
      double precision :: tmpdp
      
      logical, save :: firstcall = .true.
      
      real :: slen, sdelta,cetmp,kmtmin,kmtmax,shearmax,buoymax,shbumax

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!-----------------------------------------------------------------------

      xbnd(:) = 1.0
      ybnd(:) = 1.0
!      zbnd(:) = 0.0

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg-1

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

      DO i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO
      DO j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO
!      DO k = kzb,kze ; zbnd(k) = 1.0 ; ENDDO
#else
      DO i = 2,nx-2 ; xbnd(i) = 1.0 ; ENDDO
      DO j = 2,ny-2 ; ybnd(j) = 1.0 ; ENDDO
!      DO k = 2,nz-2 ; zbnd(k) = 1.0 ; ENDDO
#endif

      IF ( ny .le. 2 ) THEN
       ybnd(:) = 1.0
      ENDIF
      
      IF ( bcy .eq. 2 ) THEN
        ybnd(:) = 1.0
      ENDIF

      IF ( nx .le. 2 ) THEN
        xbnd(:) = 1.0
      ENDIF

      IF ( bcx .eq. 2 ) THEN
        xbnd(:) = 1.0
      ENDIF

!-----------------------------------------------------------------------
! Compute base state mixing

      kmbase(:) = 0.0
      kmbasesh(:) = 0.0
      kmbasebb(:) = 0.0
      prz(:) = Pr

      IF( rmbasekm .eq. 1 ) THEN

        prz(:) = 0.0
        
        kzb = -ng+1
        kze = ktile+ng
        if (kzbeg .eq. nzbeg) kzb = 1 
        if (kzend .eq. nzend) kze = kzend-kzbeg
        
        DO k = kzb,kze
         km1 = k-1
         kp1 = k+1
         if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
         if (kzend .eq. nzend) kp1 = min(k+1,nz-1)

!          dte = 0.5*((tinit(kp1,1)-tinit(k,1))*gz(k+1,4) + (tinit(k,1)-tinit(km1,1))*gz(k,4))
!          bb  = g * dte / tinit(k,1)
          dte = 0.5*((tinit(kp1,1)*(1.0 + 0.61*tinit(kp1,2))-tinit(k,1)*(1.0 + 0.61*tinit(k,2)))*gz(k+1,4) + & 
                    (tinit(k,1)*(1.0 + 0.61*tinit(k,2))-tinit(km1,1)*(1.0 + 0.61*tinit(km1,2)))*gz(k,4))
          bb  = g * dte / ( tinit(k,1)*(1.0 + 0.61*tinit(k,2)) )
          ush = 0.5*((uinit(kp1)-uinit(k))*gz(k+1,4) + (uinit(k)-uinit(km1))*gz(k,4))
          vsh = 0.5*((vinit(kp1)-vinit(k))*gz(k+1,4) + (vinit(k)-vinit(km1))*gz(k,4))
          bsh = ush**2 + vsh**2

          IF( mix_type .eq. 0 ) THEN
           kmbase(k) = kmbasefac*sqrt(max(bsh - Pr*bb,0.))
          ELSE
!           kmbase(k) = max(bsh - Pr*bb,0.0)
           IF ( kmbasefac > 1.0 ) THEN
             kmbase(k) = kmbasefac* max(kmbasefac*bsh - (1.0 - Min(1.0,(kmbasefac - 1.0)))*Pr*bb,0.0) 
           ELSE
             kmbase(k) = kmbasefac*max(bsh - Pr*bb,0.0) 
           ENDIF
!           prztmp(k) = Max(Pr, 1./(0.95*min(bb/(bsh+1.0e-10),50.) ))
           
          ENDIF
          
          kmbaserm(k) = kmbase(k)
          
       ENDDO !! k
       
!       DO k = kzb,kze
!         km1 = k-1
!         kp1 = k+1
!         if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
!         if (kzend .eq. nzend) kp1 = min(k+1,nz-1)
!         
!         prz(k) = Max(prztmp(k), Max( prztmp(km1),prztmp(kp1) ) )
!         
!       ENDDO



           IF ( firstcall ) THEN
             firstcall = .false.
           ENDIF

      ENDIF !! ( rmbasekm .eq. 1 )

! Base state O'brien profile for boundary layer fluxes.. 

      IF( rmbasekm .eq. 1 .and. bcz .eq. 2 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
        
       DO k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        IF( gz(k,1) .le. zPBLhgt ) kPBL = k
       ENDDO

       ustar = sqrt(drag*(uinit(1)**2 + vinit(1)**2))
       zPBL  = gz(kPBL+1,1)
       zSFC  = gz(1,1)
       dzPBL = zPBL-zSFC
       kmSFC = zSFC * ustar

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .ge. kPBL)  kze = kPBL
        
       do k = kzb,kze
#else
       DO k = 1,kPBL
#endif
        z = gz(k,1)
        kmPBL = kmbase(kPBL+1)
        kmbase(k) = kmPBL + ((zPBL-z)/dzPBL)**2 * (kmSFC-kmPBL + (z-zSFC)*(ustar + 2.*(kmSFC-kmPBL)/dzPBL))
       ENDDO !! k

      ENDIF !! ( rmbasekm .eq. 1 .and. bcz .eq. 2 )

!-----------------------------------------------------------------------
! SFC fluxes

      IF ( kzbeg == nzbeg ) THEN
        km(:,:,kusfc) = 0.0
        km(:,:,kvsfc) = 0.0

      IF( bcz .ge. 1 ) THEN

! Base state surface fluxes

       angle0     = atan2(vinit(1),uinit(1))
       vsq0       = uinit(1)**2 + vinit(1)**2

       kmbase(-1) = drag * vsq0 * cos(angle0)  ! U-base state friction flux
       kmbase( 0) = drag * vsq0 * sin(angle0)  ! V-base state friction flux

#ifdef MPI
       ixb = -ng+1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -ng+1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg
       
       DO j = jyb,jye
        DO i = ixb,ixe
#else
       DO j = 1,ny-1
        DO i = 1,nx-1
#endif
        uavg       = 0.5 * (u(i,j,1) + u(i+1,j,  1)) + ugrid
        vavg       = 0.5 * (v(i,j,1) + v(i  ,j+1,1)) + vgrid
        angle(i,j) = atan2(vavg,uavg)
        vsq(i,j)   = uavg**2 + vavg**2

        ENDDO
       ENDDO

#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       DO j = jyb,jye
        DO i = ixb,ixe
         IF( ixbeg-1+i .gt. 1 ) km(i,j,kusfc) = drag*0.5*(vsq(i,j)*cos(angle(i,j))+vsq(i-1,j)*cos(angle(i-1,j))) - kmbase(-1) ! U - Ubar flux
         IF( jybeg-1+j .gt. 1 ) km(i,j,kvsfc) = drag*0.5*(vsq(i,j)*sin(angle(i,j))+vsq(i,j-1)*sin(angle(i,j-1))) - kmbase( 0) ! V - Vbar flux
        ENDDO
       ENDDO
#else
       DO j = 1,ny-1
        DO i = 1,nx-1

         IF( i .ne. 1 ) km(i,j,kusfc) = drag*0.5*(vsq(i,j)*cos(angle(i,j))+vsq(i-1,j)*cos(angle(i-1,j))) - kmbase(-1) ! U - Ubar flux
         IF( j .ne. 1 ) km(i,j,kvsfc) = drag*0.5*(vsq(i,j)*sin(angle(i,j))+vsq(i,j-1)*sin(angle(i,j-1))) - kmbase( 0) ! V - Vbar flux

        ENDDO
       ENDDO
#endif

      ENDIF !! ( bcz .ge. 1 )

        IF ( (bcz == -1 .or. bcz == -2 ) .and. allocated( tt13sl ) .and. allocated( tt23sl ) ) THEN

       ixb = -1
       ixe = itile+1
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+1
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       DO j = jyb,jye
        DO i = ixb,ixe
         IF( ixbeg-1+i .gt. 1 ) km(i,j,kusfc) = tt13sl(i,j) ! /gz(1,3) ! U  flux
         IF( jybeg-1+j .gt. 1 ) km(i,j,kvsfc) = tt23sl(i,j) ! /gz(1,3) ! V  flux
        ENDDO
       ENDDO
        
        ENDIF

      ENDIF !! kzbeg == nzbeg


!-----------------------------------------------------------------------
! Special case, IF Mix_type < 0, set Km to a value of ahighk

      IF( mix_type .eq. -1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       do k = kzb,kze
        km(:,:,k) = ahighk
       enddo
#else
        km(:,:,1:nz-1) = ahighk
#endif

!-----------------------------------------------------------------------
! Mix_type = 0, compute Smagorinsky formula

      ELSEIF( mix_type .eq. 0 ) THEN

#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -1
       kze = ktile+2
       if(kzbeg .eq. nzbeg) kzb = 1
       if(kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,dx,dy,dz,mlen)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(k,3)
        
          IF( len_type .eq. 0 ) mlen = Cm*(dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dz
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(k,1)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
    
          km(i,j,k) = (mlen**2)*(sqrt(max(shear(i,j,k)-Pr*buoy(i,j,k),0.))-kmbase(k))*xbnd(i)*ybnd(j)
          km(i,j,k) = max(min(km(i,j,k), ahighk),0.0) 

         ENDDO
        ENDDO
       ENDDO

! O'brien profile for boundary layer fluxes..

       IF( bcz .eq. 2 ) THEN
#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -1
       kze = ktile+2
       if(kzbeg .eq. nzbeg) kzb = 1
       if(kzend .ge. kPBL) kze=kPBL

       DO k = kzb,kze
        z = gz(k,1)
        DO j = jyb,jye
         DO i = ixb,ixe
#else
        DO k = 1,kPBL
         z = gz(k,1)
         DO j = 1,ny-1
          DO i = 1,nx-1
#endif
           ustar = sqrt(drag*vsq(i,j))
           kmSFC = zSFC * ustar
           kmPBL = kmbase(kPBL+1)
           km(i,j,k) = kmPBL + ((zPBL-z)/dzPBL)**2 * (kmSFC-kmPBL + (z-zSFC)*(ustar + 2.*(kmSFC-kmPBL)/dzPBL))
           km(i,j,k) = (km(i,j,k) - kmbase(k))*xbnd(i)*ybnd(j)

          ENDDO
         ENDDO
        ENDDO

       ENDIF !! ( bcz .eq. 2 )

!-----------------------------------------------------------------------
! Mix_type = 1, compute RHS of TKE EQ AND UPDATE!

      ELSEIF( mix_type .eq. 1 ) THEN
      
      kmtmax = 0.0
      kmtmin = 0.0
      shearmax = 0.0
      buoymax = 0.0
      shbumax = 0.0

#ifdef MPI
       ixb = -ng+1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -ng+1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -ng+1
       kze = ktile+2
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,dx,dy,dz,mlen,dyl,div)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(k,3)
        
          ! Deardorff 1980 variable mixing length
          sdelta = (dx*dy*dz)**(1./3.)
          
          IF ( .true. .and. buoy(i,j,k) > 0.0 .and. km(i,j,k) > 1.e-3 ) THEN
            slen = 0.76*km(i,j,k)/Sqrt(buoy(i,j,k))
            slen = Min( slen, dz )
            !slen = Min( slen, sdelta )
          ELSE
            slen = dz ! Min(dx, dy)
          ENDIF
          
          cetmp = Ce !  0.19 + 0.51*slen/sdelta
          
          kmt(i,j,k) = cm*slen*km(i,j,k) ! momentum KM
          kmv(i,j,k) = cm*slen*km(i,j,k) ! momentum KM - vertical
          kmh(i,j,k) = cm*Sqrt(dx*dy)*km(i,j,k) ! momentum KM - horizontal

          kht(i,j,k) = (1.0 + 2.0*slen/sdelta)*kmt(i,j,k) ! scalar KM
          khh(i,j,k) = (1.0 + 2.0*slen/sdelta)*kmh(i,j,k) ! scalar KM
          khv(i,j,k) = (1.0 + 2.0*slen/sdelta)*kmv(i,j,k) ! scalar KM

          mlen = Cm*slen
          
          div =  gx(i,3)*(u(i+1,j  ,k  ) - u(i,j,k))      &
               + gy(j,3)*(v(i  ,j+1,k  ) - v(i,j,k))      &
               + gz(k,3)*(w(i  ,j  ,k+1) - w(i,j,k))      


          
            ft(i,j,k) = ft(i,j,k) + 0.5*mlen*(shear(i,j,k) - Pr*buoy(i,j,k) - kmbase(k) )*xbnd(i)*ybnd(j)  &
                        - 0.5*cetmp*Cm*km(i,j,k)**2/mlen - (1./3.)*div*km(i,j,k)

           IF ( iturbenhance > 0 ) THEN ! tke dissipation for E instead of sqrtE
             tke_diss(i,j,k) =  cetmp*Cm*km(i,j,k)**3/mlen 
           ENDIF

!           kmtmin = Max(kmtmin, ft(i,j,k) )
!           shearmax = Max(shearmax, shear(i,j,k) )
!           buoymax = Max( buoymax, buoy(i,j,k)*Pr )
!           shbumax = Max( shbumax, shear(i,j,k) - Pr*buoy(i,j,k))

         ENDDO
        ENDDO
       ENDDO
       
!       write(91,*) 'kmtmax, ftmax, sh, bu = ', kmtmax,kmtmin,shearmax,buoymax,shbumax

! Fluxes from surface layer if BCZ == 1 and mixing at first grid pt

       IF( bcz .eq. 1 ) THEN

#ifdef MPI
        ixb = -ng+1
        ixe = itile+ng
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = -ng+1
        jye = jtile+ng
        if (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO j = jyb,jye
         DO i = ixb,ixe
#else
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(1,3)

          IF( len_type .eq. 0 ) mlen = Cm*(dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dz
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(1,1)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
          IF( len_type .eq. 4 ) THEN
            dyl = (dx*dy*dz)**(1./3.)
            IF ( gz(k,1) .gt. 1400. ) THEN
             mlen = Cm*dyl
            ELSE
             mlen = Cm / (1./(0.4*gz(k,1)) + 1./dyl)
            ENDIF
          ENDIF
 
          ustar     = max(drag*(vsq(i,j) - vsq0),0.0)
          km(i,j,1) = km(i,j,1) + dt                                   &
                   * (0.5*(mlen**2)*ustar/((0.35*gz(1,1))**2)          &
                   +  (km(i,j,2)+km(i,j,1))*(km(i,j,2)-km(i,j,1))*gz(2,4)*gz(1,3))
 
         ENDDO
        ENDDO
 
       ENDIF !! ( bcz .eq. 1 )

      ENDIF !! ( mix_type .eq. 1 )

      RETURN

      END SUBROUTINE COMPUTETKE

!-----------------------------------------------------------------------------
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE COMPUTEKM    ////////////////////
!     
!-----------------------------------------------------------------------------
!       
! Created by LJW: 11-01-99
! Latest update:  05-06-04
! Renamed from computetke to computekm
! tke_type == 1 'km' is Km (original default)
!-----------------------------------------------------------------------------
!
! km  => TKE:  The RK3 integration variable for Km
!     => SMG:  The array where Km is stored
!
! kt  => TKE:  The 'n' time level of km, used in RK3 TKE integration
!     => SMG:  Not used.
!
! ft  => TKE:  Holds mixing and advection time tendencies for km.
!     => SMG:  Not used.
!
!-----------------------------------------------------------------------------
      SUBROUTINE COMPUTEKM(km,kt,shear,buoy,ft,u,v,w,tke_diss,              &
                            uinit,vinit,pinit,tinit,kmbaserm,               &
                            gx,gy,gz,dt,ugrid,vgrid,nx,ny,nz,ns)

      USE PARAM_MODULE, only : ng, bcx, bcy, bcz, Cm, Ce, Pr, len_type, mix_type, &
     &                         zPBLhgt, drag, ahighk, Lmax, rmbasekm, g, kmgen_thresh,kmbasefac
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile, my_rank
      USE FORCE_MODULE, only : tt23sl, tt13sl, td13sl, td23sl, twt3sl
      USE MICRO_MODULE, only : iturbenhance

      implicit none

! Passed variables

      integer, INTENT(IN)    :: nx, ny, nz, ns
      real,    INTENT(INOUT) :: km    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: kt    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(INOUT) :: ft    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: u     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: v     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: w     (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: shear (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: buoy  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: uinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: vinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: pinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: tinit(-ng+1:nz+ng,2)
      real,    INTENT(INOUT) :: kmbaserm(-ng+1:nz+ng)
      real,    INTENT(IN)    :: gx(-ng+1:nx+ng,4)
      real,    INTENT(IN)    :: gy(-ng+1:ny+ng,4)
      real,    INTENT(IN)    :: gz(-ng+1:nz+ng,4)
      real,    INTENT(IN)    :: ugrid, vgrid, dt
      real,    INTENT(INOUT) :: tke_diss(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

! Local variables
       
      integer i, j, k, kPBL, im1, ip1, jm1, jp1, km1, kp1
      real dx, dy, dz
      real dte, tm, rc, moist, temp, qvs, aterm
      real hgt, dz0, term1d, bb, bsh, ush, vsh 
!      double precision :: mlen       
      real :: mlen       
      real ustar, kmPBL, kmSFC, zPBL, zSFC, dzPBL, z
      real vsq0, angle0, ubase_flux, vbase_flux
!      double precision :: div
      real :: div
      real dxl,dyl,dzl

      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng), zbnd(-ng+1:nz+ng)

      real vsq(-ng+1:nx+ng,-ng+1:ny+ng), angle(-ng+1:nx+ng,-ng+1:ny+ng), uavg, vavg

      real kmbase(-ng+1:nz+ng)
      real kmbasesh(-ng+1:nz+ng)
      real kmbasebb(-ng+1:nz+ng)
      real prz(-ng+1:nz+ng)
      real prztmp(-ng+1:nz+ng)

      integer, parameter :: kusfc = KUSFCX
      integer, parameter :: kvsfc = KVSFCX
      
      double precision :: tmpdp
      
      logical, save :: firstcall = .true.
      
      real :: slen, sdelta,cetmp,kmtmin,kmtmax,shearmax,buoymax,shbumax

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!-----------------------------------------------------------------------

      xbnd(:) = 1.0
      ybnd(:) = 1.0
!      zbnd(:) = 0.0

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg-1

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

      DO i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO
      DO j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO
!      DO k = kzb,kze ; zbnd(k) = 1.0 ; ENDDO
#else
      DO i = 2,nx-2 ; xbnd(i) = 1.0 ; ENDDO
      DO j = 2,ny-2 ; ybnd(j) = 1.0 ; ENDDO
!      DO k = 2,nz-2 ; zbnd(k) = 1.0 ; ENDDO
#endif

      IF ( ny .le. 2 ) THEN
       ybnd(:) = 1.0
      ENDIF
      
      IF ( bcy .eq. 2 ) THEN
        ybnd(:) = 1.0
      ENDIF

      IF ( nx .le. 2 ) THEN
        xbnd(:) = 1.0
      ENDIF

      IF ( bcx .eq. 2 ) THEN
        xbnd(:) = 1.0
      ENDIF

!-----------------------------------------------------------------------
! Compute base state mixing

      kmbase(:) = 0.0
      kmbasesh(:) = 0.0
      kmbasebb(:) = 0.0
      prz(:) = Pr

      IF( rmbasekm .eq. 1 ) THEN

        prz(:) = 0.0
        
        kzb = -ng+1
        kze = ktile+ng
        if (kzbeg .eq. nzbeg) kzb = 1 
        if (kzend .eq. nzend) kze = kzend-kzbeg
        
        DO k = kzb,kze
         km1 = k-1
         kp1 = k+1
         if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
         if (kzend .eq. nzend) kp1 = min(k+1,nz-1)

!          dte = 0.5*((tinit(kp1,1)-tinit(k,1))*gz(k+1,4) + (tinit(k,1)-tinit(km1,1))*gz(k,4))
!          bb  = g * dte / tinit(k,1)
          dte = 0.5*((tinit(kp1,1)*(1.0 + 0.61*tinit(kp1,2))-tinit(k,1)*(1.0 + 0.61*tinit(k,2)))*gz(k+1,4) + & 
                    (tinit(k,1)*(1.0 + 0.61*tinit(k,2))-tinit(km1,1)*(1.0 + 0.61*tinit(km1,2)))*gz(k,4))
          bb  = g * dte / ( tinit(k,1)*(1.0 + 0.61*tinit(k,2)) )
          ush = 0.5*((uinit(kp1)-uinit(k))*gz(k+1,4) + (uinit(k)-uinit(km1))*gz(k,4))
          vsh = 0.5*((vinit(kp1)-vinit(k))*gz(k+1,4) + (vinit(k)-vinit(km1))*gz(k,4))
          bsh = ush**2 + vsh**2

          IF( mix_type .eq. 0 ) THEN
           kmbase(k) = kmbasefac*sqrt(max(bsh - Pr*bb,0.))
          ELSE
!           kmbase(k) = max(bsh - Pr*bb,0.0)
           IF ( kmbasefac > 1.0 ) THEN
             kmbase(k) = kmbasefac* max(kmbasefac*bsh - (1.0 - Min(1.0,(kmbasefac - 1.0)))*Pr*bb,0.0) 
           ELSE
             kmbase(k) = kmbasefac*max(bsh - Pr*bb,0.0) 
           ENDIF
!           prztmp(k) = Max(Pr, 1./(0.95*min(bb/(bsh+1.0e-10),50.) ))
           
          ENDIF
          
          kmbaserm(k) = kmbase(k)
          
       ENDDO !! k
       
!       DO k = kzb,kze
!         km1 = k-1
!         kp1 = k+1
!         if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
!         if (kzend .eq. nzend) kp1 = min(k+1,nz-1)
!         
!         prz(k) = Max(prztmp(k), Max( prztmp(km1),prztmp(kp1) ) )
!         
!       ENDDO



           IF ( firstcall ) THEN
             firstcall = .false.
           ENDIF

      ENDIF !! ( rmbasekm .eq. 1 )

! Base state O'brien profile for boundary layer fluxes.. 

      IF( rmbasekm .eq. 1 .and. bcz .eq. 2 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
        
       DO k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        IF( gz(k,1) .le. zPBLhgt ) kPBL = k
       ENDDO

       ustar = sqrt(drag*(uinit(1)**2 + vinit(1)**2))
       zPBL  = gz(kPBL+1,1)
       zSFC  = gz(1,1)
       dzPBL = zPBL-zSFC
       kmSFC = zSFC * ustar

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .ge. kPBL)  kze = kPBL
        
       do k = kzb,kze
#else
       DO k = 1,kPBL
#endif
        z = gz(k,1)
        kmPBL = kmbase(kPBL+1)
        kmbase(k) = kmPBL + ((zPBL-z)/dzPBL)**2 * (kmSFC-kmPBL + (z-zSFC)*(ustar + 2.*(kmSFC-kmPBL)/dzPBL))
       ENDDO !! k

      ENDIF !! ( rmbasekm .eq. 1 .and. bcz .eq. 2 )

!-----------------------------------------------------------------------
! SFC fluxes

      IF ( kzbeg == nzbeg ) THEN
        km(:,:,kusfc) = 0.0
        km(:,:,kvsfc) = 0.0

      IF( bcz .ge. 1 ) THEN

! Base state surface fluxes

       angle0     = atan2(vinit(1),uinit(1))
       vsq0       = uinit(1)**2 + vinit(1)**2

       kmbase(-1) = drag * vsq0 * cos(angle0)  ! U-base state friction flux
       kmbase( 0) = drag * vsq0 * sin(angle0)  ! V-base state friction flux

#ifdef MPI
       ixb = -ng+1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -ng+1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg
       
       DO j = jyb,jye
        DO i = ixb,ixe
#else
       DO j = 1,ny-1
        DO i = 1,nx-1
#endif
        uavg       = 0.5 * (u(i,j,1) + u(i+1,j,  1)) + ugrid
        vavg       = 0.5 * (v(i,j,1) + v(i  ,j+1,1)) + vgrid
        angle(i,j) = atan2(vavg,uavg)
        vsq(i,j)   = uavg**2 + vavg**2

        ENDDO
       ENDDO

#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       DO j = jyb,jye
        DO i = ixb,ixe
         IF( ixbeg-1+i .gt. 1 ) km(i,j,kusfc) = drag*0.5*(vsq(i,j)*cos(angle(i,j))+vsq(i-1,j)*cos(angle(i-1,j))) - kmbase(-1) ! U - Ubar flux
         IF( jybeg-1+j .gt. 1 ) km(i,j,kvsfc) = drag*0.5*(vsq(i,j)*sin(angle(i,j))+vsq(i,j-1)*sin(angle(i,j-1))) - kmbase( 0) ! V - Vbar flux
        ENDDO
       ENDDO
#else
       DO j = 1,ny-1
        DO i = 1,nx-1

         IF( i .ne. 1 ) km(i,j,kusfc) = drag*0.5*(vsq(i,j)*cos(angle(i,j))+vsq(i-1,j)*cos(angle(i-1,j))) - kmbase(-1) ! U - Ubar flux
         IF( j .ne. 1 ) km(i,j,kvsfc) = drag*0.5*(vsq(i,j)*sin(angle(i,j))+vsq(i,j-1)*sin(angle(i,j-1))) - kmbase( 0) ! V - Vbar flux

        ENDDO
       ENDDO
#endif

      ENDIF !! ( bcz .ge. 1 )

        IF ( (bcz == -1 .or. bcz == -2 ) .and. allocated( tt13sl ) .and. allocated( tt23sl ) ) THEN

       ixb = -1
       ixe = itile+1
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+1
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       DO j = jyb,jye
        DO i = ixb,ixe
         IF( ixbeg-1+i .gt. 1 ) km(i,j,kusfc) = tt13sl(i,j) ! /gz(1,3) ! U  flux
         IF( jybeg-1+j .gt. 1 ) km(i,j,kvsfc) = tt23sl(i,j) ! /gz(1,3) ! V  flux
        ENDDO
       ENDDO
        
        ENDIF

      ENDIF !! kzbeg == nzbeg


!-----------------------------------------------------------------------
! Special case, IF Mix_type < 0, set Km to a value of ahighk

      IF( mix_type .eq. -1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       do k = kzb,kze
        km(:,:,k) = ahighk
       enddo
#else
        km(:,:,1:nz-1) = ahighk
#endif

!-----------------------------------------------------------------------
! Mix_type = 0, compute Smagorinsky formula

      ELSEIF( mix_type .eq. 0 ) THEN

#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -1
       kze = ktile+2
       if(kzbeg .eq. nzbeg) kzb = 1
       if(kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,dx,dy,dz,mlen)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(k,3)
        
          IF( len_type .eq. 0 ) mlen = Cm*(dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dz
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(k,1)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
    
          km(i,j,k) = (mlen**2)*(sqrt(max(shear(i,j,k)-Pr*buoy(i,j,k),0.))-kmbase(k))*xbnd(i)*ybnd(j)
          km(i,j,k) = max(min(km(i,j,k), ahighk),0.0) 

         ENDDO
        ENDDO
       ENDDO

! O'brien profile for boundary layer fluxes..

       IF( bcz .eq. 2 ) THEN
#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -1
       kze = ktile+2
       if(kzbeg .eq. nzbeg) kzb = 1
       if(kzend .ge. kPBL) kze=kPBL

       DO k = kzb,kze
        z = gz(k,1)
        DO j = jyb,jye
         DO i = ixb,ixe
#else
        DO k = 1,kPBL
         z = gz(k,1)
         DO j = 1,ny-1
          DO i = 1,nx-1
#endif
           ustar = sqrt(drag*vsq(i,j))
           kmSFC = zSFC * ustar
           kmPBL = kmbase(kPBL+1)
           km(i,j,k) = kmPBL + ((zPBL-z)/dzPBL)**2 * (kmSFC-kmPBL + (z-zSFC)*(ustar + 2.*(kmSFC-kmPBL)/dzPBL))
           km(i,j,k) = (km(i,j,k) - kmbase(k))*xbnd(i)*ybnd(j)

          ENDDO
         ENDDO
        ENDDO

       ENDIF !! ( bcz .eq. 2 )

!-----------------------------------------------------------------------
! Mix_type = 1, compute RHS of TKE EQ AND UPDATE!

      ELSEIF( mix_type .eq. 1 ) THEN

      kmtmax = 0.0
      kmtmin = 0.0
      shearmax = 0.0
      buoymax = 0.0
      shbumax = 0.0

#ifdef MPI
       ixb = -ng+1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -ng+1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -ng+1
       kze = ktile+2
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       DO k = kzb,kze
        DO j = jyb,jye
         DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,dx,dy,dz,mlen,dyl,div)
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(k,3)
        
          IF( len_type .eq. 0 ) mlen = Cm * (dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dz
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(k,1)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
          IF( len_type .eq. 4 ) THEN
            dyl = (dx*dy*dz)**(1./3.)
            IF ( gz(k,1) .gt. 1400. ) THEN
             mlen = Cm*dyl
            ELSE
             mlen = Cm / (1./(0.4*gz(k,1)) + 1./dyl)
            ENDIF
          ENDIF
          div =  gx(i,3)*(u(i+1,j  ,k  ) - u(i,j,k))      &
               + gy(j,3)*(v(i  ,j+1,k  ) - v(i,j,k))      &
               + gz(k,3)*(w(i  ,j  ,k+1) - w(i,j,k))      

!          km(i,j,k) = kt(i,j,k) + dt*(0.5*(mlen**2)*(shear(i,j,k) - Pr*buoy(i,j,k) - kmbase(k))*xbnd(i)*ybnd(j)  &
!                  - 0.5*Ce*Cm*(km(i,j,k)/mlen)**2 - (1./3.)*div*km(i,j,k) )

          IF ( .true. ) THEN
          
            ft(i,j,k) = ft(i,j,k) + 0.5*(mlen**2)*(shear(i,j,k) - Pr*buoy(i,j,k) - kmbase(k))*xbnd(i)*ybnd(j)  &
                        - 0.5*Ce*Cm*(km(i,j,k)/mlen)**2 - (1./3.)*div*km(i,j,k) 

!           kmtmin = Max(kmtmin, ft(i,j,k) )
!           shearmax = Max(shearmax, shear(i,j,k) )
!           buoymax = Max( buoymax, buoy(i,j,k)*Pr )
!           shbumax = Max( shbumax, shear(i,j,k) - Pr*buoy(i,j,k) )

          ELSE
          ! test code
            tmpdp = 0.0
          
            IF ( xbnd(i)*ybnd(j) > 0.0 ) THEN
              tmpdp = 0.5*(mlen**2)*(shear(i,j,k) - Pr*buoy(i,j,k) - kmbase(k)) - 0.5*Ce*Cm*(km(i,j,k)/mlen)**2 - &
                      (1./3.)*div*km(i,j,k) 
!              tmpdp = 0.5*(mlen**2)*(shear(i,j,k) - Prz(k)*buoy(i,j,k) ) - 0.5*Ce*Cm*(km(i,j,k)/mlen)**2 - (1./3.)*div*km(i,j,k) 
            ELSE
              tmpdp =  Min( 0.0,  - 0.5*Ce*Cm*(km(i,j,k)/mlen)**2 - (1./3.)*div*km(i,j,k) )
            ENDIF

            IF ( tmpdp > 0.0 .and. tmpdp < kmgen_thresh ) tmpdp = 0.0
                  
            ft(i,j,k) = ft(i,j,k) + tmpdp
          
          ENDIF
          
           IF ( iturbenhance > 0 ) THEN ! Need dissipation for E, not sqrtE
             tke_diss(i,j,k) =  Ce*Cm*km(i,j,k)**3/mlen**4 
           ENDIF
          
!          km(i,j,k) = max(km(i,j,k), 0.0)
!          km(i,j,k) = min(km(i,j,k), ahighk)
          
! Limit on km is now done on kt in solver, so no need to do it here
!          IF ( iturbenhance >= 1 ) THEN
!            km(i,j,k) = max(min(km(i,j,k), mlen*ahighk),0.0) ! (ERM) essentially treats ahighk as max value of Sqrt(tke), 
!                                                             ! effectively unlimiting KM, but km is then limited by ahighk
!                                                             ! in the actual fluxes. This is done to preserve gradients in
!                                                             ! TKE (i.e., prevent TKE from being limited) for the rain 
!                                                             ! autoconversion turbulence-enhancement option
!
!          ELSE
!!            km(i,j,k) = max(min(km(i,j,k), ahighk),0.0)
!          ENDIF
         ENDDO
        ENDDO
       ENDDO

!       write(91,*) 'kmtmax, ftmax, sh, bu = ', kmtmax,kmtmin,shearmax,buoymax,shbumax

! Fluxes from surface layer if BCZ == 1 and mixing at first grid pt

       IF( bcz .eq. 1 ) THEN

#ifdef MPI
        ixb = -ng+1
        ixe = itile+ng
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = -ng+1
        jye = jtile+ng
        if (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        DO j = jyb,jye
         DO i = ixb,ixe
#else
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          dz = 1.0/gz(1,3)

          IF( len_type .eq. 0 ) mlen = Cm*(dx*dy*dz)**(1./3.)
          IF( len_type .eq. 1 ) mlen = Cm * dz
          IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(1,1)) + 1./Lmax)
          IF( len_type .eq. 3 ) mlen = Cm * Lmax
          IF( len_type .eq. 4 ) THEN
            dyl = (dx*dy*dz)**(1./3.)
            IF ( gz(k,1) .gt. 1400. ) THEN
             mlen = Cm*dyl
            ELSE
             mlen = Cm / (1./(0.4*gz(k,1)) + 1./dyl)
            ENDIF
          ENDIF
 
          ustar     = max(drag*(vsq(i,j) - vsq0),0.0)
          km(i,j,1) = km(i,j,1) + dt                                   &
                   * (0.5*(mlen**2)*ustar/((0.35*gz(1,1))**2)          &
                   +  (km(i,j,2)+km(i,j,1))*(km(i,j,2)-km(i,j,1))*gz(2,4)*gz(1,3))
 
         ENDDO
        ENDDO
 
       ENDIF !! ( bcz .eq. 1 )

      ENDIF !! ( mix_type .eq. 1 )

      RETURN

      END SUBROUTINE COMPUTEKM

!-----------------------------------------------------------------------------
! 
! SUBROUTINE TKE_SHEAR
!
!-----------------------------------------------------------------------------
      SUBROUTINE TKE_SHEAR(shear,u,v,w,gx,gy,gz,nx,ny,nz,ng1)

      USE GRID_MODULE
      USE PARAM_MODULE, only : ng , bcx, bcy, bcz, kmshear_thresh
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile
      USE FORCE_MODULE, only : td13sl, td23sl

      implicit none

! Passed variables

      integer, INTENT(IN)    :: nx, ny, nz, ng1
      real,    INTENT(INOUT) :: shear(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real,    INTENT(IN)    :: u    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: v    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: w    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      TYPE(VARIABLE) :: gx(4) , gy(4), gz(4)

! Local variables

      integer i, j, k, im1, ip1, jm1, jp1, km1, kp1
      
      real :: dx,dy, fac

      double precision :: def12_0, def12_1, def12_2, def12_3
      double precision :: def13_0, def13_1, def13_2, def13_3
      double precision :: def23_0, def23_1, def23_2, def23_3
      double precision :: ux, vy, wz

      integer, parameter :: kusfc = KUSFCX
      integer, parameter :: kvsfc = KVSFCX

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
#endif

! Compute MAG(DEF) 

#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       kzb = -1
       kze = ktile+2
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       do k = kzb,kze
        km1 = k-1
        kp1 = k+1
        if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
        if (kzend .eq. nzend) kp1 = min(k+1,nz-1)

       do j = jyb,jye
        jm1 = j-1
        jp1 = j+1
        if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jm1 = max(j-1,1)
        if (jyend .eq. nyend .and. bcy .ne. 2 ) jp1 = min(j+1,ny-1)

       do i = ixb,ixe
        im1 = i-1
        ip1 = i+1
        if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) im1 = max(i-1,1)
        if (ixend .eq. nxend .and. bcx .ne. 2 ) ip1 = min(i+1,nx-1)
#else

!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,im1,ip1,jm1,jp1,km1,kp1,ux,vy,wz, &
!$OMP  def12_0, def12_1, def12_2, def12_3, &
!$OMP  def13_0, def13_1, def13_2, def13_3, &
!$OMP  def23_0, def23_1, def23_2, def23_3  )
      DO k = 1,nz-1

       km1 = max(k-1,1)
       kp1 = min(k+1,nz-1)

       DO j = 1,ny-1

       IF ( bcy .ne. 2 ) THEN
        jm1 = max(j-1,1)
        jp1 = min(j+1,ny-1)
        IF ( bcy .eq. 0 ) jp1 = j+1
       ELSE
        jm1 = j-1
        jp1 = j+1
       ENDIF

        DO i = 1,nx-1

         IF ( bcx .ne. 2 ) THEN
          im1 = max(i-1,1)
          ip1 = min(i+1,nx-1)
          IF ( bcx .eq. 0 ) im1 = i-1
         ELSE
          im1 = i-1
          ip1 = i+1
         ENDIF
#endif


         IF ( k == 1 ) THEN
           dx = 1.0/gx(3)%flt1d(i)
           dy = 1.0/gy(3)%flt1d(j) 
         ENDIF

         ux = (u(i+1,j,  k  )-u(i,j,k))*gx(3)%flt1d(i)

         vy = (v(i,  j+1,k  )-v(i,j,k))*gy(3)%flt1d(j)

         wz = (w(i,  j,  k+1)-w(i,j,k))*gz(3)%flt1d(k)

         def12_0 = (gx(4)%flt1d(i)   * ( v(i,j,k) - v(im1,j,  k) )   &
                + gy(4)%flt1d(j)   * ( u(i,j,k) - u(i,jm1,k) ) )

         def12_1 = gx(4)%flt1d(i+1) * ( v(ip1,j,k) - v(i,j,  k  ) )  &
                + gy(4)%flt1d(j) * ( u(i+1,j,k) - u(i+1,jm1,k) )

         def12_2 = gx(4)%flt1d(i) * ( v(i,j+1,k) - v(im1,j+1,k) )  &
                 + gy(4)%flt1d(j+1) * ( u(i,jp1,k) - u(i,j,  k  ) )

         def12_3 = gx(4)%flt1d(i+1) * ( v(ip1,j+1,k) - v(i,j+1,k  ) ) &
                 + gy(4)%flt1d(j+1) * ( u(i+1,jp1,k) - u(i+1,j,  k) )

         def13_0 = gx(4)%flt1d(i)   * ( w(i,j,k) - w(im1, j,k  ) )  &
                 + gz(4)%flt1d(k)   * ( u(i,j,k) - u(i,j,km1  ) )
         

         def13_1 = gx(4)%flt1d(i+1) * ( w(ip1,j,k) - w(i,j,k  ) )  &
                 + gz(4)%flt1d(k)   * ( u(i+1,j,k) - u(i+1,j,km1) )

         def13_2 = gx(4)%flt1d(i) * ( w(i,j,k+1) - w(im1,j,k+1) )  &
                 + gz(4)%flt1d(k+1) * ( u(i,j,kp1) - u(i,j,k   ) )

         def13_3 = gx(4)%flt1d(i+1) * ( w(ip1,j,k+1) - w(i,  j,k+1) )  &
                 + gz(4)%flt1d(k+1) * ( u(i+1,j,kp1) - u(i+1,j,k  ) )
         IF ( k == 1 .and. (bcz == -1 .or. bcz == -3) .and. allocated( td13sl ) ) THEN
           ! def13_0 will be zero at k=1, anyway, because w(k=1) is zero and u(i,j,1) - u(i,j,1) = 0
           def13_0 = td13sl(i,j)
           def13_1 = td13sl(ip1,j)
           
         ENDIF


         def23_0 = gy(4)%flt1d(j)   * ( w(i,j,k) - w(i,jm1,k  ) )  &
                 + gz(4)%flt1d(k)   * ( v(i,j,k) - v(i,j,  km1) )


         def23_1 = gy(4)%flt1d(j+1) * ( w(i,jp1,k) - w(i,j,  k  ) )  &
                 + gz(4)%flt1d(k)   * ( v(i,j+1,k) - v(i,j+1,km1) )

         def23_2 = gy(4)%flt1d(j) * ( w(i,j,k+1) - w(i,jm1,k+1) )  &
                 + gz(4)%flt1d(k+1) * ( v(i,j,kp1) - v(i,j,  k  ) )

         def23_3 = gy(4)%flt1d(j+1) * ( w(i,jp1,k+1) - w(i,j,  k+1) )  &
                 + gz(4)%flt1d(k+1) * ( v(i,j+1,kp1) - v(i,j+1,k  ) )

         IF ( k == 1 .and. (bcz == -1 .or. bcz == -3) .and. allocated( td23sl ) ) THEN
           ! def23_0 will be zero at k=1, anyway, because w(k=1) is zero and v(i,j,1) - v(i,j,1) = 0
           def23_0 = td23sl(i,j)
           def23_1 = td23sl(i,jp1)

         ENDIF

         IF ( k == 1 .and. ( dx > 5000. .or. dy > 5000. ) .and. (bcz == -1 .or. bcz == -3) ) THEN 
         
           IF ( Max(dx,dy) > 10000. ) THEN
             fac = 0.0
           ELSE
             fac = (10000. - Max(dx,dy) )/5000.
           ENDIF
           def13_0 = fac*def13_0
           def13_1 = fac*def13_1
           def23_0 = fac*def23_0
           def23_1 = fac*def23_1
         ENDIF

         shear(i,j,k) = 2.0d0*ux**2  &
                      + 2.0d0*vy**2  &
                      + 2.0d0*wz**2  &
                      + (.25d0*(def12_0 + def12_1 + def12_2 + def12_3))**2  &
                      + (.25d0*(def13_0 + def13_1 + def13_2 + def13_3))**2  &
                      + (.25d0*(def23_0 + def23_1 + def23_2 + def23_3))**2

         IF ( shear(i,j,k) < kmshear_thresh ) shear(i,j,k) = 0.0

!         IF ( k == 1 .and. ( dx > 7500. .or. dy > 7500. ) ) THEN 
!           shear(i,j,k) = 0.0
!         ENDIF

        ENDDO
       ENDDO
      ENDDO

      RETURN

      END SUBROUTINE TKE_SHEAR

!-----------------------------------------------------------------------------
! 
! SUBROUTINE TKE_BUOY
!
!-----------------------------------------------------------------------------
      SUBROUTINE TKE_BUOY(buoy,s,st,pinit,tinit,gz,nx,ny,nz,ns)

       USE GRID_MODULE
      USE PARAM_MODULE, only : ng , bcx, bcy, bcz, g, cthres, Lav, rd, epsilon, Cp, iusebvfreq
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile
      USE FORCE_MODULE, only : twt3sl

      implicit none

! Passed variables

      integer, INTENT(IN)    :: nx, ny, nz, ns
      real,    INTENT(INOUT) :: buoy  (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real,    INTENT(IN)    :: st    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
      TYPE(VARIABLE)         :: s(ns)
      real,    INTENT(IN)    :: pinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: tinit(-ng+1:nz+ng)
      real,    INTENT(IN)    :: gz(-ng+1:nz+ng,4)

! Local variables

      integer i, j, k, im1, ip1, jm1, jp1, km1, kp1
      integer n
      double precision :: pk, dte, tm, rc, moist, temp, qvs, aterm
      double precision :: dq3
      integer lv, lc, li, lr

      integer, parameter :: kusfc = KUSFCX
      integer, parameter :: kvsfc = KVSFCX
      real, allocatable :: thetal(:,:,:)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
! #ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
! #endif

       lv = 0
       lc = 0
       li = 0
       lr = 0
       
       DO n = 1,ns
       
        IF ( s(n)%buotype .eq. 2 .or. s(n)%name .eq. 'QV ' ) lv = n
        IF ( s(n)%name .eq. 'QC ' ) lc = n
        IF ( s(n)%name .eq. 'QR ' ) lr = n
        IF ( s(n)%name .eq. 'QI ' ) li = n
       
       ENDDO
       
!       print*, 'lv,lc,li = ',lv,lc,li
!       IF ( li .eq. 0 ) li = lc
       
       IF ( lv .eq. 0 ) THEN 
! Buoyancy term for dry model

#ifdef MPI
        ixb = -1
        ixe = itile+2
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = -1
        jye = jtile+2
        if (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        kzb = -1
        kze = ktile+2
        if (kzbeg .eq. nzbeg) kzb = 1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        DO k = kzb,kze
          km1 = k-1
          kp1 = k+1
          if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
          if (kzend .eq. nzend) kp1 = min(k+1,nz-1)
         DO j = jyb,jye
          DO i = ixb,ixe
#else
        DO k = 1,nz-1

          km1 = max(k-1,1)
          kp1 = min(k+1,nz-1)

!$OMP PARALLEL DO IF ( ny .gt. 2 ), DEFAULT(SHARED), PRIVATE(i,j,dte)
         DO j = 1,ny-1

          DO i = 1,nx-1
#endif
           IF ( iusebvfreq == 0 ) THEN
             dte  = 0.5*((st(i,j,kp1,1)-st(i,j,k,1))*gz(k+1,4) + (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4))
             buoy(i,j,k) = g * dte / st(i,j,k,1) 
           ELSEIF ( iusebvfreq == -1 ) THEN
             dte  = (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)
             buoy(i,j,k) = g * dte / st(i,j,k,1) 
           ELSE
             dte = alog( st(i,j,k,1)/st(i,j,km1,1) )*gz(k,4)
             buoy(i,j,k) = g * dte
           ENDIF

          ENDDO
         ENDDO
        ENDDO

       ELSE

         IF ( iusebvfreq >= 1 ) THEN
           allocate( thetal (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
        
         ixb = -1
         ixe = itile+2
         if (ixbeg .eq. nxbeg) ixb = 1
         if (ixend .eq. nxend) ixe = ixend-ixbeg
 
         jyb = -1
         jye = jtile+2
         if (jybeg .eq. nybeg) jyb = 1
         if (jyend .eq. nyend) jye = jyend-jybeg

         kzb = -1
         kze = ktile+2
         if (kzbeg .eq. nzbeg) kzb = 1
         if (kzend .eq. nzend) kze = kzend-kzbeg

         DO k = kzb,kze
          DO j = jyb,jye
           DO i = ixb,ixe
             rc = 0.0
             IF ( lc > 0 ) rc = rc + st(i,j,k,lc)
             IF ( lr > 0 ) rc = rc + st(i,j,k,lr)

             temp  = pinit(k)*st(i,j,k,1)

             thetal(i,j,k) = st(i,j,k,1) - rc*Lav/(Cp*pinit(k))
           ENDDO
          ENDDO
         ENDDO
        
         ENDIF

#ifdef MPI
        ixb = -1
        ixe = itile+2
        if (ixbeg .eq. nxbeg) ixb = 1
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        jyb = -1
        jye = jtile+2
        if (jybeg .eq. nybeg) jyb = 1
        if (jyend .eq. nyend) jye = jyend-jybeg

        kzb = -1
        kze = ktile+2
        if (kzbeg .eq. nzbeg) kzb = 1
        if (kzend .eq. nzend) kze = kzend-kzbeg

        DO k = kzb,kze
         km1 = k-1
         kp1 = k+1
         if (kzbeg .eq. nzbeg) km1 = max(k-1,1)
         if (kzend .eq. nzend) kp1 = min(k+1,nz-1)

         pk   = 1.0e5*pinit(k)**3.509

         DO j = jyb,jye
          DO i = ixb,ixe
#else
        DO k = 1,nz-1

         km1 = max(k-1,1)
         kp1 = min(k+1,nz-1)

         pk   = 1.0d5*pinit(k)**3.509d0

!$OMP PARALLEL  DO IF ( ny .gt. 2 ), DEFAULT(SHARED), PRIVATE(i,j,moist,rc,temp,qvs,dq3,aterm,dte,n)
         DO j = 1,ny-1
          DO i = 1,nx-1
#endif
           moist = 0.0
           IF ( lc .ne. 0 ) THEN
             IF ( li == 0 ) THEN
               rc = st(i,j,k,lc) 
               IF ( rc .ge. cthres ) moist = 1.0
             ELSE
               rc = st(i,j,k,lc) + st(i,j,k,li)
               IF ( rc .ge. cthres ) moist = 1.0
             ENDIF
           ENDIF

           temp  = pinit(k)*st(i,j,k,1)
           qvs   = 380.*exp(17.27*(temp-273.)/(temp-36.)) / pk

           dq3 = 0.0
           
           IF ( moist .eq. 1.0 ) THEN
           aterm = (1.0 + moist * Lav*qvs/(rd*temp)) / (1.0 + moist * epsilon*Lav**2*qvs/(Cp*rd*temp**2))

           IF ( iusebvfreq == 0 ) THEN
           dte  = 0.5*( (st(i,j,kp1,1)-st(i,j,k,1))*gz(k+1,4) + (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)  &
                + moist*Lav/Cp*((st(i,j,kp1,lv)/pinit(kp1)-st(i,j,k,  lv)/pinit(k  ))*gz(kp1,4)      &
                               +(st(i,j,k,  lv)/pinit(k  )-st(i,j,km1,lv)/pinit(km1))*gz(k  ,4)) ) 
           ELSEIF ( iusebvfreq == -1 ) THEN
           dte  =  (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)  &
                + moist*Lav/Cp*((st(i,j,k,  lv)/pinit(k  )-st(i,j,km1,lv)/pinit(km1))*gz(k  ,4)) 
           ELSE
           ! need to update, and then also below
           dte  =  (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)  &
                + moist*Lav/Cp*((st(i,j,k,  lv)/pinit(k  )-st(i,j,km1,lv)/pinit(km1))*gz(k  ,4)) 
!           dte  = 0.5*( (st(i,j,kp1,1)-st(i,j,k,1))*gz(k+1,4) + (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)  &
!                + moist*Lav/Cp*((st(i,j,kp1,lv)/pinit(kp1)-st(i,j,k,  lv)/pinit(k  ))*gz(k+1,4)      &
!                               +(st(i,j,k,  lv)/pinit(k  )-st(i,j,km1,lv)/pinit(km1))*gz(k  ,4)) )            
           ENDIF

             DO n=3,ns
              IF ( s(n)%buotype .eq. 3 ) THEN
               dq3 = dq3 + 0.5*((st(i,j,kp1,n)-st(i,j,k,n))*gz(k+1,4) + (st(i,j,k,n)-st(i,j,km1,n))*gz(k,4))
              ENDIF
             ENDDO

!           buoy(i,j,k) = g * (aterm*dte/st(i,j,k,1) - (dq3))

           ELSE
             dq3 = 0.0
             aterm = (1.0)

             IF ( iusebvfreq == 0 ) THEN
                dte  = 0.5*((st(i,j,kp1,1)-st(i,j,k,1))*gz(k+1,4) + (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4))
             ELSEIF ( iusebvfreq == -1 ) THEN
                dte  = (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)
             ELSE
                dte  = (st(i,j,k,1)-st(i,j,km1,1))*gz(k,4)
!                dte = alog( st(i,j,k,1)/st(i,j,km1,1) )*gz(k,4)
             ENDIF

           ENDIF
           
           IF ( (bcz == -1 .or. bcz == -4 .or. bcz == -2 ) .and. k == 1 .and. allocated( twt3sl ) ) THEN
             aterm = 1.0
             dq3 = 0.0
             dte = twt3sl(i,j)
           ENDIF

           IF ( iusebvfreq <= 0 ) THEN
             buoy(i,j,k) = g *  (aterm*dte/st(i,j,k,1) - dq3 )
           ELSE
             buoy(i,j,k) = g *  (aterm*dte/st(i,j,k,1) - dq3 )
!             buoy(i,j,k) = g * (aterm*dte - dq3)
           ENDIF

          ENDDO
         ENDDO
        ENDDO

         IF ( iusebvfreq >= 1 ) THEN
           deallocate( thetal )
         ENDIF

       ENDIF
       

      RETURN
     
      END SUBROUTINE TKE_BUOY

!-----------------------------------------------------------------------------
! 
! SUBROUTINE MIX_VELO - does a deformational-based mixing of momentum
!
!-----------------------------------------------------------------------------
      SUBROUTINE MIX_VELO(u,v,w,fu,fv,fw,kht,kmh,kmv,km,e,ubase,vbase,gx,gy,gz, &
                          ugrid,vgrid,nx,ny,nz,                     &
                          t0,den)

      USE PARAM_MODULE, only: ng, ahighk, len_type, tke_type, Lmax, Cm, bcx, bcy, bcz
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile, nproci, my_rank

      implicit none

      integer nx, ny, nz
      real dx, dy, dz, ugrid, vgrid
      real kht(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real kmh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real kmv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real km(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real e (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real ubase(-ng+1:nz+ng), vbase(-ng+1:nz+ng)
      real gx(-ng+1:nx+ng,4), gy(-ng+1:ny+ng,4), gz(-ng+1:nz+ng,4)
      real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: den(-ng+1:nz+ng,2)

! Local variables

      integer i, j, k, kbnd, i1, j1
      integer km1, kp1
      double precision :: rkm0, rkm1, rkm2, rkm0b
      real uavg, vavg, mlen
      real deninv, deninvw
      integer jm1, jp1, im1, ip1
      integer ibc, jbc
      real dyl
      real, allocatable :: t1(:,:,:)
!      real, allocatable :: d12(:,:,:), d13(:,:,:), d23(:,:,:)
!      real, allocatable :: d11(:,:,:), d22(:,:,:), d33(:,:,:)

      double precision :: def13ik, def13ip, def13kp
      double precision :: def23jk, def23jp, def23kp
      double precision :: def12ij, def12ip, def12jp
      
      double precision :: tmpdp,tmpdp1,tmpdp2,tmpdp3,tmpdp4

      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng), zbnd(-ng+1:nz+ng)
      
      logical, parameter :: ifudo = .true.
      logical, parameter :: ifvdo = .true.
      logical, parameter :: ifwdo = .true.
      
      integer :: k1

      integer, parameter :: kusfc = KUSFCX
      integer, parameter :: kvsfc = KVSFCX

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
!
!  #########################################################################
!

      allocate( t1(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d11(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d22(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d33(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d12(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d13(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
!       allocate( d23(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )

      t1(:,:,:) = 0.0

      xbnd(:) = 0.0
      ybnd(:) = 0.0
!      zbnd(:) = 0.0

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg-1

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

!      do k = kzb,kze ; zbnd(k) = 1.0 ; ENDDO
      do j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO
      do i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO
#else
      DO i = 2,nx-2 ; xbnd(i) = 1.0 ; ENDDO
      DO j = 2,ny-2 ; ybnd(j) = 1.0 ; ENDDO
!      DO k = 2,nz-2 ; zbnd(k) = 1.0 ; ENDDO
#endif
      IF ( ny .le. 2 ) THEN
        j1 = 1
      ELSE
        j1 = 1
      ENDIF
      
      jbc = 0
      IF ( bcy .eq. 2 ) THEN
        j1 = 1
        jbc = 1
        ybnd(:) = 1.0
      ENDIF

      IF ( nx .le. 2 ) THEN
        i1 = 1
      ELSE
        i1 = 1
      ENDIF

      ibc = 0
      IF ( bcx .eq. 2 ) THEN
        i1 = 1
        ibc = 1
        xbnd(:) = 1.0
      ENDIF
      
! Compute E from Km for diagonal terms

! #ifdef MPI
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

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,dx,dy,dz,mlen,dyl,im1,jm1)
      DO k = kzb,kze
        km1 = Max(1,k-1)
       DO j = jyb,jye
         jm1 = j-1
         if (jybeg .eq. nybeg) jm1 = max(j-1,1)
         IF ( bcy .eq. 2 ) THEN
          jm1 = j-1
         ENDIF
        DO i = ixb,ixe
        im1 = i-1
        if (ixbeg .eq. nxbeg .and. bcx /= 2 ) im1 = max(i-1,1)        

! #else
!       DO k = 1,nz-1
!        DO j = 1,ny-1
!         DO i = 1,nx-1
! #endif
         dx = 1.0/gx(i,3)
         dy = 1.0/gy(j,3)
         dz = 1.0/gz(k,3)

         IF ( tke_type == 1 ) THEN
         IF( len_type .eq. 0 ) mlen = Cm*(dx*dy*dz)**(1./3.)
         IF( len_type .eq. 1 ) mlen = Cm * dz
         IF( len_type .eq. 2 ) mlen = Cm / (1./(0.4*gz(k,1)) + 1./Lmax)
         IF( len_type .eq. 3 ) mlen = Cm * Lmax
         IF( len_type .eq. 4 ) THEN
           dyl = (dx*dy*dz)**(1./3.)
           IF ( gz(k,1) .gt. 1400. ) THEN
             mlen = Cm*dyl
           ELSE
             mlen = Cm / (1./(0.4*gz(k,1)) + 1./dyl)
           ENDIF
         ENDIF
 
           e(i,j,k) = den(k,1)* (2./3.) * (Min(ahighk,km(i,j,k)) / mlen)**2
          
          ELSE
          ! km is sqrt(tke)
            e(i,j,k) = den(k,1)* (2./3.) * km(i,j,k)**2
          
          ENDIF
         
         
!            d13(i,j,k) = gx(i,4)   * ( w(i,j,k) - w(im1,j,k  ) )   &
!                   + gz(k,4)   * ( u(i,j,k) - u(i,  j,km1) )
! 
!            d23(i,j,k) = gy(j,4)   * ( w(i,j,k) - w(i,jm1,k  ) )    &
!                   + gz(k,4)   * ( v(i,j,k) - v(i,j,  km1) )
! 
         ENDDO
        ENDDO
       ENDDO

      k1 = -ng+1
      if (kzbeg .eq. nzbeg) k1 = 1

#ifdef MPI
      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = -1
      if (kzend .eq. nzend) kze = kzend-kzbeg
      

      DO k = kzb,kze
#else
      DO k = -1,nz-1
#endif
         IF ( tke_type == 1 ) THEN
           t0(:,:,k) = den(Max(k1,k),1)*Min(ahighk,km(:,:,k))
           t1(:,:,k) = den(Max(k1,k),1)*Min(ahighk,km(:,:,k))
         ELSE
           t0(:,:,k) = den(Max(k1,k),1)*Min(ahighk,kmh(:,:,k))
           t1(:,:,k) = den(Max(k1,k),1)*Min(ahighk,kmv(:,:,k))
         ENDIF
      ENDDO

! U-momentum --> d/dx (Km D11)

      IF ( nx .gt. 2 .and. ifudo ) THEN !{
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2) ixb = 2
      if (ixend .eq. nxend .and. bcx .ne. 2) ixe = ixend-ixbeg

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = j1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-j1

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO k = kzb,kze
       deninv = 1.0/den(k,1)
       DO j = jyb,jye
        DO i = ixb,ixe 
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv)
      DO k = 1,nz-1
       deninv = 1.0/den(k,1)
       DO j = j1,ny-j1
        DO i = 2,nx-1
#endif
         ip1 = Min(i+1,ixend-ixbeg+1*ibc)

        tmpdp1 = (e(i,j,k) - e(i-1,j,k))
        tmpdp2 = t0(i,  j,k) * (u(i+1,j,k)-u(i,j,k)) * gx(i+1,3)
        tmpdp3 = t0(i-1,j,k) * (u(i-1,j,k)-u(i,j,k)) * gx(i,  3)

         tmpdp = gx(i,4)* (tmpdp2 + tmpdp3  - tmpdp1 )* deninv

        fu(i,j,k) = fu(i,j,k) + tmpdp

!          IF (  k == 27 .and. i == 18 ) THEN
!            write(0,*) 'j,tmp = ',tmpdp,tmpdp2,tmpdp3, tmpdp1
!          ENDIF

!          fu(i,j,k) = fu(i,j,k) + gx(i,4)  &
!                    * (t0(i,  j,k) * (u(i+1,j,k)-u(i,j,k)) * gx(ip1,3)  &
!                      +t0(i-1,j,k) * (u(i-1,j,k)-u(i,j,k)) * gx(i,  3)  &
!                    - (e(i,j,k) - e(i-1,j,k)) ) * deninv

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

      do k = kzb,kze
       deninv = 1.0/den(k,1)
      do j = jyb,jye
#else
      DO k = 1,nz-1
       deninv = 1.0/den(k,1)
       DO j = 1,ny-1
#endif
         i = 1
         fu(i,j,k) = fu(i,j,k) + gx(i,4)  &
                   * (t0(i,  j,k) * (u(i+1,j,k)-u(i,j,k)) * gx(i+1,3)  &
                     +t0(nx-1,j,k) * (u(i-1,j,k)-u(i,j,k)) * gx(i,  3)  &
                   - (e(i,j,k) - e(nx-1,j,k)) ) * deninv

#ifdef MPI
         i = nxend
#else
         i = nx
#endif
         fu(i,j,k) = fu(1,j,k)
!         fu(i,j,k) = fu(i,j,k) + gx(i,4)  &
!                   * (t0(1,  j,k) * (u(i+1,j,k)-u(i,j,k)) * gx(1,3)  &
!                     +t0(i-1,j,k) * (u(i-1,j,k)-u(i,j,k)) * gx(i,  3)  &
!                   - (e(1,j,k) - e(i-1,j,k)) ) 
       ENDDO
      ENDDO
      ENDIF
      
      ENDIF   !}


! V-momentum --> d/dy (Km D22)

     IF ( ny .gt. 2 .and. ifvdo ) THEN !{
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg .and. bcx .ne. 2 ) ixb = i1
      if (ixend .eq. nxend .and. bcx .ne. 2 ) ixe = ixend-ixbeg+1-i1

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg .and. bcy .ne. 2 ) jyb = 2
      if (jyend .eq. nyend .and. bcy .ne. 2 ) jye = jyend-jybeg

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO k = kzb,kze
       deninv = 1.0/den(k,1)
       DO j = jyb,jye
       jp1 = Min(j+1, jyend-jybeg+1*jbc)
        DO i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv,jp1)
      DO k = 1,nz-1
       deninv = 1.0/den(k,1)
       DO j = 2,ny-1
       jp1 = Min(j+1, jyend-jybeg+1*jbc)
        DO i = i1,nx-i1
#endif
        tmpdp1 = e(i,j,k) - e(i,j-1,k)
        tmpdp2 = t0(i,j,  k) * (v(i,j+1,k)-v(i,j,k)) * gy(j+1,3)
        tmpdp3 = t0(i,j-1,k) * (v(i,j-1,k)-v(i,j,k)) * gy(j,  3)

         tmpdp = gy(j,4)* (tmpdp2 + tmpdp3  - tmpdp1 )* deninv

        fv(i,j,k) = fv(i,j,k) + tmpdp

!         fv(i,j,k) = fv(i,j,k) + gy(j,4)  &
!                   * (t0(i,j,  k) * (v(i,j+1,k)-v(i,j,k)) * gy(jp1,3)  &
!                     +t0(i,j-1,k) * (v(i,j-1,k)-v(i,j,k)) * gy(j,  3)  &
!                   - (e(i,j,k) - e(i,j-1,k)) ) * deninv

        ENDDO
       ENDDO
      ENDDO

#ifndef MPI
      IF ( bcy .eq. 2 ) THEN

      DO k = 1,nz-1
       deninv = 1.0/den(k,1)
        DO i = 1,nx-1
         j = 1
         fv(i,j,k) = fv(i,j,k) + gy(j,4)  &
                  * (t0(i,j,   k) * (v(i,j+1,k)-v(i,j,k)) * gy(j+1,3)  &
                    +t0(i,ny-1,k) * (v(i,j-1,k)-v(i,j,k)) * gy(j,  3)  &
                  - (e(i,j,k) - e(i,ny-1,k)) ) * deninv
         j = ny

         fv(i,j,k) = fv(i,1,k)

       ENDDO
      ENDDO

      ENDIF
#endif
      
      ELSE ! } ( ny .gt. 2 .and. ifvdo )
       ! fv(:,:,:) = 0.0
      ENDIF
      

! W-momentum --> d/dz (Km D33)

      IF ( ifwdo ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = i1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-i1

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = j1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-j1

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
       deninvw = 1.0/den(k,2)
      do j = jyb,jye
      do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninvw)
      DO k = 2,nz-1
       deninvw = 1.0/den(k,2)
       DO j = j1,ny-j1
        DO i = i1,nx-i1
#endif
        fw(i,j,k) = fw(i,j,k) + gz(k,4)   &
                  * (t1(i,j,k  ) * (w(i,j,k+1)-w(i,j,k)) * gz(k+1,3)  &
                    +t1(i,j,k-1) * (w(i,j,k-1)-w(i,j,k)) * gz(k,  3) &
                  - (e(i,j,k) - e(i,j,k-1)) )*deninvw

        ENDDO
       ENDDO
      ENDDO
      ENDIF

! U-equation needs  d/dy(du/dy + dv/dx),  d/dz(du/dz + dw/dx)    
!                                  ^                     ^

! U-momentum --> d/dz (Km D13)
! W-momentum --> d/dx (Km D13)

      IF ( nx .gt. 2 ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1*ibc

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = j1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-j1

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
       deninv = 1.0/den(k,1)
       deninvw = 1.0/den(k,2)
      do j = jyb,jye
      do i = ixb,ixe
        im1 = i-1
        ip1 = i+1
        if (ixbeg .eq. nxbeg) im1 = max(i-1,1)
        if (ixend .eq. nxend) then
          ip1 = min(i+1,nx-1)
        endif
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,im1,ip1,deninv,deninvw,def13ik,def13ip,def13kp,rkm0,rkm1,rkm2)
      DO k = 2,nz-1
       deninv = 1.0/den(k,1)
       deninvw = 1.0/den(k,2)
       DO j = j1,ny-j1
        DO i = 2-ibc,nx-1+1*ibc
         im1 = Max( 1, i - 1 )
         ip1 = Min( nx-1, i + 1 )
#endif
         IF ( bcx .eq. 2 ) THEN
           im1 = i - 1
           ip1 = i + 1
         ENDIF

         def13ik = gx(i,4)   * ( w(i,j,k) - w(im1,j,k  ) )   &
                 + gz(k,4)   * ( u(i,j,k) - u(i,  j,k-1) )

         def13ip = gx(ip1,4) * ( w(ip1,j,k) - w(i,  j,k  ) )  &
                 +         gz(k,4)   * ( u(ip1,j,k) - u(ip1,j,k-1) )

         def13kp = gx(i,  4) * ( w(i,j,k+1) - w(im1,j,k+1) )  &
                 + gz(k+1,4) * ( u(i,j,k+1) - u(i,  j,k  ) )

        tmpdp1 = t1(i,j,k)+t1(im1,j,k)
        tmpdp2 = t1(i,j,k)+t1(i,j,k-1)
        rkm0 = 0.25 * (tmpdp1 + t1(i,j,k-1) + t1(im1,j,k-1))
        rkm0b = 0.25 * (tmpdp2+t1(im1,j,k) + t1(im1,j,k-1))
        rkm1 = 0.25 * (tmpdp1 + t1(i,j,k+1) + t1(im1,j,k+1))
        rkm2 = 0.25 * (tmpdp2+t1(ip1,j,k) + t1(ip1,j,k-1))

! d/dz(dw/dx)
!        IF( k .ne. nz-1 ) fu(i,j,k) = fu(i,j,k) + gz(k,3) * gx(i,4)   &
!                                    * (km1*(w(i,j,k+1) - w(i-1,j,k+1)) - km0*(w(i,j,k) - w(i-1,j,k)))

!        IF( k .ne. nz-1 .and. i.ne. im1 ) fu(i,j,k) = fu(i,j,k) + gz(k,3) * (km1*def13kp - km0*def13ik)
#ifdef MPI
        IF( kzbeg-1+k .ne. nzend-1 .and. ifudo ) fu(i,j,k) = fu(i,j,k) + gz(k,3) * (rkm1*def13kp - rkm0*def13ik)*deninv
#else
        IF( k .ne. nz-1 .and. ifudo ) fu(i,j,k) = fu(i,j,k) + gz(k,3) * (rkm1*def13kp - rkm0*def13ik)*deninv
#endif
! d/dx(du/dz)
!           fw(i,j,k) = fw(i,j,k) + gx(i,3) * gz(k,4)  &
!                            * (km2*(u(i+1,j,k) - u(i+1,j,k-1)) - km0*(u(i,j,k) - u(i,j,k-1)))

!        IF( i .ne. ip1 ) fw(i,j,k) = fw(i,j,k) + gx(i,3) * (km2*def13ip - km0*def13ik)

         IF ( ifwdo )  fw(i,j,k) = fw(i,j,k) + gx(i,3) * (rkm2*def13ip - rkm0b*def13ik)*xbnd(i)*deninvw
!         IF ( ifwdo )  fw(i,j,k) = fw(i,j,k) + gx(i,3) * (rkm2*d13(i+1,j,k) - rkm0*d13(i,j,k))*xbnd(i)*deninvw

!        IF( i .eq. ip1 ) fw(i,j,k) = fw(i,j,k) + gx(i,3) * (km2*def13ip - km0*(gz(k,4)*( u(i,j,k) - u(i,j,k-1) )))

        ENDDO
       ENDDO
      ENDDO

! Lower boundary condition if bcz = 1

      IF( bcz /= 0 .and. ifudo .and. kzbeg == nzbeg ) THEN
#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = 2-ibc
       if (ixend .eq. nxend) ixe = ixend-ixbeg+ibc

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = j1
       if (jyend .eq. nyend) jye = jyend-jybeg+1-j1

       deninv = 1.0/den(1,1)
       DO j = jyb,jye
        DO i = ixb,ixe
#else
       deninv = 1.0/den(1,1)
       DO j = j1,ny-j1
        DO i = 2-ibc,nx-1+ibc
#endif
         def13kp = gx(i,4) * ( w(i,j,2) - w(i-1,j,2) )  &
                 + gz(2,4) * ( u(i,j,2) - u(i,  j,1) )

        rkm1 = 0.25 * (t1(i,j,1)+t1(i-1,j,1) + t1(i,j,2) + t1(i-1,j,2))
!        fu(i,j,1) = fu(i,j,1) + gz(1,3)*(km1*(w(i,j,2)-w(i-1,j,2))*gx(i,4) - km(i,j,-1))   ! km(j,-1,i) contain (U-Ubar) surface flux
        fu(i,j,1) = fu(i,j,1) + gz(1,3)*(rkm1*def13kp - t0(i,j,kusfc)) *deninv      ! rkm(i,j,-1) contain (U-Ubar) surface flux

        ENDDO
       ENDDO

      ENDIF
      
      ENDIF ! nx .gt. 2
      
! V-momentum --> d/dz (Km D23)
! W-momentum --> d/dy (Km D23)

      IF ( ny .gt. 2 ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = i1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-i1

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1*jbc

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
       deninv = 1.0/den(k,1)
       deninvw = 1.0/den(k,2)

       kp1 = k+1
       if (kzend .eq. nzend) kp1 = min(k+1,nz-1)

      do j = jyb,jye
       jm1 = j-1
       jp1 = j+1
       if (jybeg .eq. nybeg) jm1 = max(j-1,1)
       if (jyend .eq. nyend) then
         jp1 = min(j+1,ny-1)
       endif
        
       IF ( bcy .eq. 2 ) THEN
        jm1 = j-1
        jp1 = j+1
       ENDIF
      do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,kp1,jm1,jp1,deninv,deninvw,def23jk,def23jp,def23kp,rkm0,rkm1,rkm2)
      DO k = 2,nz-1
       deninv = 1.0/den(k,1)
       deninvw = 1.0/den(k,2)
       kp1 = Min( k + 1, nz - 1 )
       DO j = 2-jbc,ny-1 + 1*jbc
         jm1 = Max( 1, j - 1 )
         jp1 = Min( ny-1, j + 1 )
         IF ( bcy .eq. 2 ) THEN
           jm1 = j - 1
           jp1 = j + 1
         ENDIF
        DO i = i1,nx-i1
#endif
         def23jk = gy(j,4)   * ( w(i,j,k) - w(i,jm1,k  ) )    &
                 + gz(k,4)   * ( v(i,j,k) - v(i,j,  k-1) )

         def23jp = gy(jp1,4) * ( w(i,jp1,k) - w(i,j,  k  ) )  &
                 +         gz(k,4)   * ( v(i,jp1,k) - v(i,jp1,k-1) )

         def23kp = gy(j,  4) * ( w(i,j,k+1) - w(i,j-1,k+1) )  &
                 + gz(k+1,4) * ( v(i,j,k+1) - v(i,j,  k  ) )

        tmpdp1 = t1(i,j,k)+t1(i,jm1,k)
        tmpdp2 = t1(i,j,k-1) + t1(i,jm1,k-1)
        tmpdp3 = t1(i,j,kp1) + t1(i,jm1,kp1)
        rkm0 = 0.25 * (tmpdp1 + tmpdp2)
        rkm1 = 0.25 * (tmpdp1 + tmpdp3)
        rkm2 = 0.25 * (t1(i,j,k)+t1(i,jp1,k) + t1(i,j,k-1) + t1(i,jp1,k-1))
        rkm0b = 0.25 * (t1(i,j,k)+t1(i,jm1,k) + t1(i,j,k-1) + t1(i,jm1,k-1))
!         rkm0 = 0.25 * (t0(i,j,k)+t0(i,jm1,k) + t0(i,j,k-1) + t0(i,jm1,k-1))
!         rkm1 = 0.25 * (t0(i,j,k)+t0(i,jm1,k) + t0(i,j,kp1) + t0(i,jm1,kp1))
!         rkm2 = 0.25 * (t0(i,j,k)+t0(i,jp1,k) + t0(i,j,k-1) + t0(i,jp1,k-1))

! d/dz(dw/dy)
!         IF( k .ne. nz-1 ) fv(i,j,k) = fv(i,j,k) + gz(k,3) * gy(j,4)  &
!                                     * (rkm1*(w(i,j,k+1)-w(i,j-1,k+1)) - rkm0*(w(i,j,k)-w(i,j-1,k)))
!        IF( k .ne. nz-1 .and. j .ne. jm1 ) fv(i,j,k) = fv(i,j,k) + gz(k,3) * (rkm1*def23kp - rkm0*def23jk)

#ifdef MPI
        IF( kzbeg-1+k .ne. nzend-1 .and. ifvdo ) THEN
!            tmpdp = gz(k,3) * (rkm1*def23kp - rkm0*def23jk)*deninv
           fv(i,j,k) = fv(i,j,k) + gz(k,3) * (rkm1*def23kp - rkm0*def23jk)*deninv
        ENDIF
#else
        IF( k .ne. nz-1 .and. ifvdo ) fv(i,j,k) = fv(i,j,k) + gz(k,3) * (rkm1*def23kp - rkm0*def23jk)*deninv
#endif
!          IF ( ifvdo .and. k == 27 .and. i == 18 ) THEN
!            write(0,*) 'j,tmp = ',tmpdp,rkm1*def23kp,rkm0*def23jk,rkm1,rkm0,def23kp,def23jk
!          ENDIF
! d/dy(dv/dz)
!          fw(i,j,k) = fw(i,j,k) + gy(j,3) * gz(k,4)  &
!                                      * (rkm2*(v(i,j+1,k)-v(i,j+1,k-1))- rkm0*(v(i,j,k)-v(i,j,k-1)))

!        IF ( j .ne. jp1 ) fw(i,j,k) = fw(i,j,k) + gy(j,3) * (rkm2*def23jp - rkm0*def23jk)

         IF ( ifwdo )  fw(i,j,k) = fw(i,j,k) + gy(j,3) * (rkm2*def23jp - rkm0b*def23jk)*ybnd(j)*deninvw
!         IF ( ifwdo )  fw(i,j,k) = fw(i,j,k) + gy(j,3) * (rkm2*d23(i,j+1,k) - rkm0*d23(i,j,k))*ybnd(j)*deninvw

!        IF ( j .eq. jp1 ) fw(i,j,k) = fw(i,j,k) + gy(j,3) * (rkm2*def23jp - rkm0*gz(k,4)*( v(i,j,k) - v(i,j, k-1)))

        ENDDO
       ENDDO
      ENDDO

! Lower boundary condition if bcz = 1

      IF( bcz /= 0 .and. ifvdo .and. kzbeg == nzbeg ) THEN
#ifdef MPI
       ixb = -1
       ixe = itile+2
       if (ixbeg .eq. nxbeg) ixb = i1
       if (ixend .eq. nxend) ixe = ixend-ixbeg+1-i1

       jyb = -1
       jye = jtile+2
       if (jybeg .eq. nybeg) jyb = 2-jbc
       if (jyend .eq. nyend) jye = jyend-jybeg+jbc

       deninv = 1.0/den(1,1)
       DO j = jyb,jye
        DO i = ixb,ixe
#else
       deninv = 1.0/den(1,1)
       DO j = 2-jbc,ny-1+jbc
        DO i = i1,nx-i1
#endif
         def23kp = gy(j,4) * ( w(i,j,2) - w(i,j-1,2) )  &
                 + gz(2,4) * ( v(i,j,2) - v(i,j,  1) )

         rkm1 = 0.25 * (t1(i,j,1)+t1(i,j-1,1) + t1(i,j,2) + t1(i,j-1,2))
!         fv(i,j,1) = fv(i,j,1) + gz(1,3)*(rkm1*(w(i,j,2)-w(i,j-1,2))*gy(j,4) - rkm(i,j,0))   ! km(j,0,i) contain (V-Vbar) surface flux
         fv(i,j,1) = fv(i,j,1) + gz(1,3)*(rkm1*def23kp - t0(i,j,kvsfc)*den(1,1)) * deninv          ! km(i,j,0) contain (V-Vbar) surface flux

        ENDDO
       ENDDO

      ENDIF

      ENDIF ! ny .gt. 2

! U-momentum --> d/dy (Km D12)
! V-momentum --> d/dx (Km D12)

      IF ( ny .gt. 2 .and. nx .gt. 2 ) THEN
#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+ibc

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+jbc

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
       deninv = 1.0/den(k,1)
      do j = jyb,jye
        jm1 = j-1
        jp1 = j+1
        if (jybeg .eq. nybeg) jm1 = max(j-1,1)
        if (jyend .eq. nyend) jp1 = min(j+1,ny-1)
         IF ( bcy .eq. 2 ) THEN
           jm1 = j - 1
           jp1 = j + 1
         ENDIF
      do i = ixb,ixe
        im1 = i-1
        ip1 = i+1
        if (ixbeg .eq. nxbeg) im1 = max(i-1,1)
        if (ixend .eq. nxend) ip1 = min(i+1,nx-1)
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,im1,ip1,jm1,jp1,deninv,def12ij,def12ip,def12jp,rkm0,rkm1,rkm2)
      DO k = 1,nz-1
       deninv = 1.0/den(k,1)
       DO j = 2-jbc,ny-1+jbc
         jm1 = Max( 1, j - 1 )
         jp1 = Min( ny-1, j + 1 )
         IF ( bcy .eq. 2 ) THEN
           jm1 = j - 1
           jp1 = j + 1
         ENDIF
        DO i = 2-ibc,nx-1+ibc
         im1 = Max( 1, i - 1 )
         ip1 = Min( nx-1, i + 1 )
         IF ( bcx .eq. 2 ) THEN
           im1 = i - 1
           ip1 = i + 1
         ENDIF
#endif
         def12ij = gx(i,4)   * ( v(i,j,k) - v(i-1,j,  k) )  &
                 + gy(j,4)   * ( u(i,j,k) - u(i,j-1,k) )

         def12ip = gx(ip1,4) * ( v(ip1,j,k) - v(i,j,  k  ) )  &
                 +         gy(j,  4) * ( u(i+1,j,k) - u(i+1,j-1,k) )

         def12jp =         gx(i,  4) * ( v(i,j+1,k) - v(i-1,j+1,k) )  &
                 + gy(jp1,4) * ( u(i,jp1,k) - u(i,j,  k  ) )

          tmpdp1 = t0(i,j,k)+t0(i-1,j,k)
          tmpdp2 = t0(i,j,k)+t0(i,j-1,k)
         
         rkm0 = 0.25*(tmpdp1+t0(i,j-1,k)+t0(i-1,j-1,k))
         rkm1 = 0.25*(tmpdp2+t0(ip1,j,k)+t0(ip1,j-1,k))
         rkm2 = 0.25*(tmpdp1+t0(i,jp1,k)+t0(i-1,jp1,k))
         rkm0b = 0.25*(tmpdp2 +t0(i-1,j,k) +t0(i-1,j-1,k))
 
 ! d/dx(du/dy)
!         IF( i .ne. nx-1 ) fv(i,j,k) = fv(i,j,k) + gx(i,3) * gy(j,4)  &
!          fv(i,j,k) = fv(i,j,k) + gx(i,3) * gy(j,4)  &
!                                     * (rkm1*(u(i+1,j,k)-u(i+1,j-1,k)) - rkm0*(u(i,j,k)-u(i,j-1,k)))
 
!        IF ( j .ne. jm1 .and. i .ne. im1) THEN
!         IF( i .ne. ip1 ) fv(i,j,k) = fv(i,j,k) + gx(i,3) * (rkm1*def12ip - rkm0*def12ij)
        IF ( ifvdo ) fv(i,j,k) = fv(i,j,k) + gx(i,3) * (rkm1*def12ip - rkm0b*def12ij)*xbnd(i) * deninv
!         IF( i .eq. ip1 ) fv(i,j,k) = fv(i,j,k) + gx(i,3) * (rkm1*def12ip - rkm0*gy(j,4) * ( u(i,j,k) - u(i,j-1,k) ))
!        ENDIF

 ! d/dy(dv/dx)
!         IF( j .ne. ny-1 ) fu(i,j,k) = fu(i,j,k) + gy(j,3) * gx(i,4)   &
!          fu(i,j,k) = fu(i,j,k) + gy(j,3) * gx(i,4)   &
!                                     * (rkm2*(v(i,j+1,k)-v(i-1,j+1,k)) - rkm0*(v(i,j,k)-v(i-1,j,k)))

!         IF ( i .ne. im1 .and. j .ne. jm1) THEN
!           IF( j .ne. jp1 ) fu(i,j,k) = fu(i,j,k) + gy(j,3) * (rkm2*def12jp - rkm0*def12ij)
            IF ( ifudo ) fu(i,j,k) = fu(i,j,k) + gy(j,3) * (rkm2*def12jp - rkm0*def12ij)*ybnd(j) * deninv
!           IF( j .eq. jp1 ) fu(i,j,k) = fu(i,j,k) + gy(j,3) * (rkm2*def12jp - rkm0*gx(i,4)*( v(i,j,k) - v(i-1,j,k) ) )
!         ENDIF

        ENDDO
       ENDDO
      ENDDO

      ENDIF

      deallocate( t1 )
!       deallocate (d11, d22, d33, d12, d13, d23 )

      RETURN
      END SUBROUTINE MIX_VELO

!-----------------------------------------------------------------------------
! 
! SUBROUTINE RAYDAMP
!
!-----------------------------------------------------------------------------
      SUBROUTINE RAYDAMP(s,fs,sinit,gz,sconst,nx,ny,nz,svar)

      USE GRID_MODULE
      USE PARAM_MODULE, only : ng , bcx, bcy, rayd_hgt, rayd_mag, hrayd_mag, pii, hsponge
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile

      implicit none

! Passed variables

      integer nx, ny, nz, is, js, ks
      real dt, sconst
      real s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real sinit(-ng+1:nz+ng), gz(-ng+1:nz+ng,4)

   TYPE(VARIABLE) :: svar

! Local variables

      integer i, j, k
      real z, zdif, zfac, top
      real rayd, raydZ(-ng+1:nz+ng), raydX(-ng+1:nx+ng), raydY(-ng+1:ny+ng)
      real xdis, ydis 

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

! Compute vertical gravity wave damping coefficient for this variable

      is = svar%istag
      js = svar%jstag
      ks = svar%kstag

      top  = gz(nz,2)
      zdif = top - rayd_hgt

      IF( zdif .lt. 1.0 ) RETURN

#ifdef MPI
      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
#else
      DO k = 1,nz-1
#endif
       raydZ(k) = 0.0
       z        = gz(k,1+ks)
       zfac     = ( z - rayd_hgt ) / zdif

       IF( z .gt. rayd_hgt ) THEN

        IF( zfac .le. 0.5 ) THEN

        raydZ(k) = 0.5 * (1.0 - cos(zfac*pii))

        ENDIF

        IF( zfac .gt. 0.5 .and. zfac .le. 1.0 ) THEN

        raydZ(k) = 0.5 * (1.0 + (zfac - 0.5)*pii)

        ENDIF

       ENDIF

      ENDDO

! Compute horizontal sponge layer near grid edge

      raydX(:) = 0.0
      raydY(:) = 0.0
 
      IF( hsponge .gt. 0 .and. ( is /= 0 .or. js /= 0 .or. ks /= 0 ) ) THEN
#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = nxend-ixbeg+is ! max(hsponge-ixbeg+1,1)

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = nyend-jybeg+js ! max(hsponge-jybeg+1,1)

      IF ( bcx .ne. 2 ) THEN

!      do i = ixb,ixe
!        xdis = 0.5*pii * float(ixbeg-1+i-1) / float(hsponge)
!        if(ixbeg-1+i .le. hsponge)             raydX(i) = cos(xdis)
!        if(ixbeg-1+i .ge. nxend-1-hsponge+is)  raydX(nx-i+is) = cos(xdis)
!      enddo

       DO i = 1,max(hsponge,1)
        xdis = 0.5*pii * float(i-1) / float(hsponge)
        IF (ixbeg-1+i .le. hsponge )  raydX(i) = cos(xdis)
        IF (ixend-1+is-hsponge .ge. nxend-1-hsponge+is) raydX(nx-i+is) = cos(xdis)
!        write(0,*) 'i,nx-i+is,xdis,raydx = ',i,nx-i+is,ixbeg-1+i,xdis,raydX(i),raydX(nx-i+is),cos(xdis)
       ENDDO

      ENDIF

      IF ( bcy .ne. 2 ) THEN

!      do j = jyb,jye
!        ydis = 0.5*pii * float(jybeg-1+j-1) / float(hsponge)
!        if(jybeg-1+j .le. hsponge)            raydY(j) = cos(ydis)
!        if(jybeg-1+j .gt. nyend-1-hsponge+js) raydY(ny-j+is) = cos(ydis)
!      enddo

       DO j = 1,max(hsponge,1)
        ydis = 0.5*pii * float(j-1) / float(hsponge)
        IF (jybeg-1+j .le. hsponge)     raydY(j) = cos(ydis)
        IF (jyend-hsponge-1+js+j .gt. nyend-1-hsponge+js)   raydY(ny-j+js) = cos(ydis)
!        write(0,*) 'j,ny-j+is,ydis,raydy = ',j,ny-j+is,jybeg-1+j,ydis,raydy(j),raydY(ny-j+js),cos(ydis)
       ENDDO

      ENDIF

#else

      IF ( bcx .ne. 2 ) THEN

       DO i = 1,max(hsponge,1)
        xdis = 0.5*pii * float(i-1) / float(hsponge)
        raydX(i)       = cos(xdis)
        raydX(nx-i+is) = cos(xdis)
       ENDDO

      ENDIF

      IF ( bcy .ne. 2 ) THEN

       DO j = 1,max(hsponge,1)
        ydis = 0.5*pii * float(j-1) / float(hsponge)
        raydY(j)       = cos(ydis)
        raydY(ny-j+js) = cos(ydis)
       ENDDO
       
      ENDIF

#endif
      ENDIF

! Do damping

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1+is
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1+js
      if (jyend .eq. nyend) jye = jyend-jybeg

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb = 1+ks
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze
       do j = jyb,jye
        do i = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$OMP PRIVATE(i,j,k,rayd)
      DO k = 1+ks,nz-1
       DO j = 1+js,ny-1
        DO i = 1+is,nx-1
#endif
         rayd = rayd_mag*raydZ(k) + hrayd_mag*(raydX(i) + raydY(j))

         fs(i,j,k) = fs(i,j,k) - rayd * (s(i,j,k) - sinit(k) + sconst)

        ENDDO
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RAYDAMP
!-----------------------------------------------------------------------
!        
! SUBROUTINE MIX_S
!    
!-----------------------------------------------------------------------------
      SUBROUTINE MIX_SCAL(s,fs,kh,sbase,gx,gy,gz,Prtl,nx,ny,nz,t0,den,tend3d,sname,tendflag) !)

      USE PARAM_MODULE, only: ahighk, bcx, bcy, ng
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile, my_rank
      USE FORCE_MODULE, only : thf3sl, tmf3sl

      implicit none

      integer nx, ny, nz
      real kh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real tend3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real sbase(-ng+1:nz+ng) 
      real gx(-ng+1:nx+ng,4), gy(-ng+1:ny+ng,4), gz(-ng+1:nz+ng,4)
      real Prtl
      integer, intent(in) :: tendflag

      real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: den(-ng+1:nz+ng,2)

      character(LEN=*), intent(IN) :: sname

! Local

      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng), zbnd(-ng+1:nz+ng)
      real deninv
      integer k, j, i, n
      integer ibc, jbc
      double precision tend,tmp
      real :: dx,dy,fac

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!-----------------------------------------------------------------------------
! Bnd flags

      ibc = 0
      jbc = 0
      xbnd(:) = 0.0
      ybnd(:) = 0.0
      zbnd(:) = 0.0

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb=2
      if (ixend .eq. nxend) ixe=ixend-ixbeg

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb=2
      if (jyend .eq. nyend) jye=jyend-jybeg

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb=2
      if (kzend .eq. nzend) kze=kzend-kzbeg-1

      DO i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO
      DO j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO
      DO k = kzb,kze ; zbnd(k) = 1.0 ; ENDDO
#else
      DO i = 2,nx-1 ; xbnd(i) = 1.0 ; ENDDO
      DO j = 2,ny-1 ; ybnd(j) = 1.0 ; ENDDO
      DO k = 2,nz-2 ; zbnd(k) = 1.0 ; ENDDO
#endif

! Check for periodic boundary:
      IF ( bcx .eq. 2 ) THEN
        xbnd(:) = 1.0
        ibc = 1
      ENDIF
      IF ( bcy .eq. 2 ) THEN
        ybnd(:) = 1.0
        jbc = 1
      ENDIF

!-----------------------------------------------------------------------------

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
      
       DO k = kzb,kze
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k)
       DO k = 1,nz-1
#endif
         t0(:,:,k) = den(k,1)*Min(ahighk,kh(:,:,k))
       ENDDO

!#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1+ibc

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg-1+jbc

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv,tend)
      do k=kzb,kze
       deninv = 1.0/den(k,1)
      do j=jyb,jye
      do i=ixb,ixe!
!#else
! !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv)
!       DO k = 1,nz-1
!       deninv = 1.0/den(k,1)
!        DO j = 2-jbc,ny-2+jbc
!         DO i = 2-ibc,nx-2+ibc
!#endif

        tend = 0.0d0
        
        IF ( nx .gt. 2 ) THEN

!        fs(i,j,k) = fs(i,j,k)+.25*gx(i,3)*Prtl   &
!          *((t0(i+1,j,k)+t0(i,j,k)+t0(i+1,j,k+1)+t0(i,j,k+1))  &
!          *(s(i+1,j,k)-s(i,j,k))*gx(i+1,4)*xbnd(i+1)           &
!           +(t0(i-1,j,k)+t0(i,j,k)+t0(i-1,j,k+1)+t0(i,j,k+1))  &
!           *(s(i-1,j,k)-s(i,j,k))*gx(i,  4)*xbnd(i  ))*deninv

!        fs(i,j,k) = fs(i,j,k)+0.5*gx(i,3)*Prtl                  &
        tend = tend + 0.5*gx(i,3)*Prtl                  &
          *((t0(i+1,j,k)+t0(i,j,k))                             &
          *(s(i+1,j,k)-s(i,j,k))*gx(i+1,4)*xbnd(i+1)            &
           +(t0(i-1,j,k)+t0(i,j,k))                             &
           *(s(i-1,j,k)-s(i,j,k))*gx(i,  4)*xbnd(i  ))*deninv
        
        ENDIF

        IF ( ny .gt. 2 ) THEN
         
!         fs(i,j,k) = fs(i,j,k)+.25*gy(j,3)*Prtl  &
!          *((t0(i,j+1,k)+t0(i,j,k)+t0(i,j+1,k+1)+t0(i,j,k+1))  &
!          *(s(i,j+1,k)-s(i,j,k))*gy(j+1,4)*ybnd(j+1)           &
!           +(t0(i,j-1,k)+t0(i,j,k)+t0(i,j-1,k+1)+t0(i,j,k+1))  &
!           *(s(i,j-1,k)-s(i,j,k))*gy(j,  4)*ybnd(j  ))*deninv

!         fs(i,j,k) = fs(i,j,k)+0.5*gy(j,3)*Prtl  &
         tend = tend + 0.5*gy(j,3)*Prtl  &
          *((t0(i,j+1,k)+t0(i,j,k))  &
          *(s(i,j+1,k)-s(i,j,k))*gy(j+1,4)*ybnd(j+1)           &
           +(t0(i,j-1,k)+t0(i,j,k))  &
           *(s(i,j-1,k)-s(i,j,k))*gy(j,  4)*ybnd(j  ))*deninv
         
        ENDIF

!         fs(i,j,k) = fs(i,j,k)+gz(k,3)*zbnd(k)*Prtl  &
         
         IF ( zbnd(k) > 0 ) tend = tend + gz(k,3)*zbnd(k)*Prtl  &
                * (0.5*(t0(i,j,k+1) + t0(i,j,k))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                  +0.5*(t0(i,j,k-1) + t0(i,j,k))*(s(i,j,k-1)-s(i,j,k))*gz(k  ,4))*deninv

         fs(i,j,k) = fs(i,j,k) + tend
         
         IF ( tendflag > 0 ) THEN
           tend3d(i,j,k) = tend
         ENDIF
         
         ENDDO
        ENDDO
       ENDDO
       
!      ixb = 1
!      ixe = itile
!      if (ixend .eq. nxend) ixe = ixend-ixbeg

!      jyb = 1
!      jye = jtile
!      if (jyend .eq. nyend) jye = jyend-jybeg


!      if (ixbeg .eq. nxbeg .and. bcx == 1) ixb = 2-ibc
!      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-2+ibc

!      jyb = -1
!      jye = jtile+2
!      if (jybeg .eq. nybeg) jyb = 2-jbc
!      if (jyend .eq. nyend) jye = jyend-jybeg+1-2+jbc

      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-2+ibc

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-2+jbc

       k = 1
       deninv = 1.0/den(1,1)

       IF ( .true. .and. sname .eq. 'TH' .and. allocated (thf3sl) ) THEN
        DO j = jyb,jye
         DO i = ixb,ixe
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)

          IF ( ( dx > 5000. .or. dy > 5000. )  ) THEN 
         
           IF ( Max(dx,dy) > 10000. ) THEN
             fac = 0.0
           ELSE
             fac = (10000. - Max(dx,dy) )/5000.
           ENDIF
         
         ELSE
           fac = 1.0
         ENDIF

           fac = 1.0

          tmp = fac*thf3sl(i,j)

          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl  &
                * ((kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv
         
         ENDDO
        ENDDO
       ENDIF

       IF ( sname .eq. 'QV' .and. allocated (tmf3sl) ) THEN
        DO j = jyb,jye
         DO i = ixb,ixe

          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          
          IF ( ( dx > 5000. .or. dy > 5000. )  ) THEN 
         
           IF ( Max(dx,dy) > 10000. ) THEN
             fac = 0.0
           ELSE
             fac = (10000. - Max(dx,dy) )/5000.
           ENDIF
         
         ELSE
           fac = 1.0
         ENDIF

           fac = 1.0

          tmp = fac*tmf3sl(i,j)


          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
                * ((kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv
                   

!          tmp = (kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))
!          IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmp =',i,j,fs(i,j,k)
         ENDDO
        ENDDO

       ELSEIF ( .false. .and. sname .eq. 'QV'  ) THEN
        ! hack code to test sfc flux of qv

      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1 ! 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1 ! -2+ibc

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1 ! 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1 ! -2+jbc

        DO j = jyb,jye
         DO i = ixb,ixe
!          tmp = tmf3sl(i,j)
          tmp = 0
!          IF ( gx(i,1) .gt. 30.e3 .and. gx(i,1) < 90.e3 .and. gy(j,1) > 30.e3 .and. gy(j,1) < 90.e3 ) THEN
!          write(0,*) 'mix: my_rank,i,j =',my_rank,i,j
          tmp = -6.5870481e-8
!          ENDIF
     !     IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmf3sl =',i,j,tmp
          
!          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
!                * ((kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
!                  - ((tmp/Prtl)*gz(k  ,4)) )
          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
                * ((kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv

!          tmp = (kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))
!          IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmp =',i,j,fs(i,j,k)
         ENDDO
        ENDDO
       
       ENDIF

      RETURN
      END SUBROUTINE MIX_SCAL

!-----------------------------------------------------------------------
!        
! SUBROUTINE MIX_SCAL_HV
!    Separate horizontal and vertical: kh -> (khh,khv)
!-----------------------------------------------------------------------------
      SUBROUTINE MIX_SCAL_HV(s,fs,khh,khv,sbase,gx,gy,gz,Prtl,nx,ny,nz,t0,den,tend3d,sname,tendflag) !)

      USE PARAM_MODULE, only: ahighk, bcx, bcy, ng
      USE COMMASMPI_MODULE, only : nxbeg, nxend, ixbeg, ixend, &
     &                             nybeg, nyend, jybeg, jyend, &
     &                             nzbeg, nzend, kzbeg, kzend, &
     &                             itile, jtile, ktile, my_rank
      USE FORCE_MODULE, only : thf3sl, tmf3sl

      implicit none

      integer nx, ny, nz
      real khh(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real khv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real fs(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real tend3d(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real sbase(-ng+1:nz+ng) 
      real gx(-ng+1:nx+ng,4), gy(-ng+1:ny+ng,4), gz(-ng+1:nz+ng,4)
      real Prtl
      integer, intent(in) :: tendflag

      real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: den(-ng+1:nz+ng,2)

      character(LEN=*), intent(IN) :: sname

! Local

      real xbnd(-ng+1:nx+ng), ybnd(-ng+1:ny+ng), zbnd(-ng+1:nz+ng)
      real deninv
      integer k, j, i, n
      integer ibc, jbc
      double precision tend,tmp
      real :: dx,dy,fac
      real, allocatable :: t1(:,:,:)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!-----------------------------------------------------------------------------
! Bnd flags

      ibc = 0
      jbc = 0
      xbnd(:) = 0.0
      ybnd(:) = 0.0
      zbnd(:) = 0.0

#ifdef MPI
      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb=2
      if (ixend .eq. nxend) ixe=ixend-ixbeg

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb=2
      if (jyend .eq. nyend) jye=jyend-jybeg

      kzb = -ng+1
      kze = ktile+ng
      if (kzbeg .eq. nzbeg) kzb=2
      if (kzend .eq. nzend) kze=kzend-kzbeg-1

      DO i = ixb,ixe ; xbnd(i) = 1.0 ; ENDDO
      DO j = jyb,jye ; ybnd(j) = 1.0 ; ENDDO
      DO k = kzb,kze ; zbnd(k) = 1.0 ; ENDDO
#else
      DO i = 2,nx-1 ; xbnd(i) = 1.0 ; ENDDO
      DO j = 2,ny-1 ; ybnd(j) = 1.0 ; ENDDO
      DO k = 2,nz-2 ; zbnd(k) = 1.0 ; ENDDO
#endif

! Check for periodic boundary:
      IF ( bcx .eq. 2 ) THEN
        xbnd(:) = 1.0
        ibc = 1
      ENDIF
      IF ( bcy .eq. 2 ) THEN
        ybnd(:) = 1.0
        jbc = 1
      ENDIF

!-----------------------------------------------------------------------------

      allocate( t1(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
      
#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
      
       DO k = kzb,kze
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k)
       DO k = 1,nz-1
#endif
         t0(:,:,k) = den(k,1)*Min(ahighk*Prtl,khh(:,:,k))
         t1(:,:,k) = den(k,1)*Min(ahighk*Prtl,khv(:,:,k))
       ENDDO

!#ifdef MPI
      ixb = -1
      ixe = itile+2
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1+ibc

      jyb = -1
      jye = jtile+2
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg-1+jbc

      kzb = -1
      kze = ktile+2
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv,tend)
      do k=kzb,kze
       deninv = 1.0/den(k,1)
      do j=jyb,jye
      do i=ixb,ixe!
!#else
! !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,deninv)
!       DO k = 1,nz-1
!       deninv = 1.0/den(k,1)
!        DO j = 2-jbc,ny-2+jbc
!         DO i = 2-ibc,nx-2+ibc
!#endif

        tend = 0.0d0
        
        IF ( nx .gt. 2 ) THEN

!        fs(i,j,k) = fs(i,j,k)+.25*gx(i,3)*Prtl   &
!          *((t0(i+1,j,k)+t0(i,j,k)+t0(i+1,j,k+1)+t0(i,j,k+1))  &
!          *(s(i+1,j,k)-s(i,j,k))*gx(i+1,4)*xbnd(i+1)           &
!           +(t0(i-1,j,k)+t0(i,j,k)+t0(i-1,j,k+1)+t0(i,j,k+1))  &
!           *(s(i-1,j,k)-s(i,j,k))*gx(i,  4)*xbnd(i  ))*deninv

!        fs(i,j,k) = fs(i,j,k)+0.5*gx(i,3)*Prtl                  &
        tend = tend + 0.5*gx(i,3)                  &
          *((t0(i+1,j,k)+t0(i,j,k))                             &
          *(s(i+1,j,k)-s(i,j,k))*gx(i+1,4)*xbnd(i+1)            &
           +(t0(i-1,j,k)+t0(i,j,k))                             &
           *(s(i-1,j,k)-s(i,j,k))*gx(i,  4)*xbnd(i  ))*deninv
        
        ENDIF

        IF ( ny .gt. 2 ) THEN
         
!         fs(i,j,k) = fs(i,j,k)+.25*gy(j,3)*Prtl  &
!          *((t0(i,j+1,k)+t0(i,j,k)+t0(i,j+1,k+1)+t0(i,j,k+1))  &
!          *(s(i,j+1,k)-s(i,j,k))*gy(j+1,4)*ybnd(j+1)           &
!           +(t0(i,j-1,k)+t0(i,j,k)+t0(i,j-1,k+1)+t0(i,j,k+1))  &
!           *(s(i,j-1,k)-s(i,j,k))*gy(j,  4)*ybnd(j  ))*deninv

!         fs(i,j,k) = fs(i,j,k)+0.5*gy(j,3)*Prtl  &
         tend = tend + 0.5*gy(j,3)  &
          *((t0(i,j+1,k)+t0(i,j,k))  &
          *(s(i,j+1,k)-s(i,j,k))*gy(j+1,4)*ybnd(j+1)           &
           +(t0(i,j-1,k)+t0(i,j,k))  &
           *(s(i,j-1,k)-s(i,j,k))*gy(j,  4)*ybnd(j  ))*deninv
         
        ENDIF

!         fs(i,j,k) = fs(i,j,k)+gz(k,3)*zbnd(k)*Prtl  &
         
         IF ( zbnd(k) > 0 ) tend = tend + gz(k,3)*zbnd(k)  &
                * (0.5*(t1(i,j,k+1) + t1(i,j,k))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                  +0.5*(t1(i,j,k-1) + t1(i,j,k))*(s(i,j,k-1)-s(i,j,k))*gz(k  ,4))*deninv

         fs(i,j,k) = fs(i,j,k) + tend
         
         IF ( tendflag > 0 ) THEN
           tend3d(i,j,k) = tend
         ENDIF
         
         ENDDO
        ENDDO
       ENDDO
       
!      ixb = 1
!      ixe = itile
!      if (ixend .eq. nxend) ixe = ixend-ixbeg

!      jyb = 1
!      jye = jtile
!      if (jyend .eq. nyend) jye = jyend-jybeg


!      if (ixbeg .eq. nxbeg .and. bcx == 1) ixb = 2-ibc
!      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-2+ibc

!      jyb = -1
!      jye = jtile+2
!      if (jybeg .eq. nybeg) jyb = 2-jbc
!      if (jyend .eq. nyend) jye = jyend-jybeg+1-2+jbc

      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-2+ibc

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-2+jbc

       k = 1
       deninv = 1.0/den(1,1)

       IF ( .true. .and. sname .eq. 'TH' .and. allocated (thf3sl) ) THEN
        DO j = jyb,jye
         DO i = ixb,ixe
          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)

          IF ( ( dx > 5000. .or. dy > 5000. )  ) THEN 
         
           IF ( Max(dx,dy) > 10000. ) THEN
             fac = 0.0
           ELSE
             fac = (10000. - Max(dx,dy) )/5000.
           ENDIF
         
         ELSE
           fac = 1.0
         ENDIF

           fac = 1.0

          tmp = fac*thf3sl(i,j)

          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl  &
                * ((khv(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv
         
         ENDDO
        ENDDO
       ENDIF

       IF ( sname .eq. 'QV' .and. allocated (tmf3sl) ) THEN
        DO j = jyb,jye
         DO i = ixb,ixe

          dx = 1.0/gx(i,3)
          dy = 1.0/gy(j,3)
          
          IF ( ( dx > 5000. .or. dy > 5000. )  ) THEN 
         
           IF ( Max(dx,dy) > 10000. ) THEN
             fac = 0.0
           ELSE
             fac = (10000. - Max(dx,dy) )/5000.
           ENDIF
         
         ELSE
           fac = 1.0
         ENDIF

           fac = 1.0

          tmp = fac*tmf3sl(i,j)


          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
                * ((khv(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv
                   

!          tmp = (kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))
!          IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmp =',i,j,fs(i,j,k)
         ENDDO
        ENDDO

       ELSEIF ( .false. .and. sname .eq. 'QV'  ) THEN
        ! hack code to test sfc flux of qv

      ixb = -ng+1
      ixe = itile+ng
      if (ixbeg .eq. nxbeg) ixb = 1 ! 2-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1 ! -2+ibc

      jyb = -ng+1
      jye = jtile+ng
      if (jybeg .eq. nybeg) jyb = 1 ! 2-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1 ! -2+jbc

        DO j = jyb,jye
         DO i = ixb,ixe
!          tmp = tmf3sl(i,j)
          tmp = 0
!          IF ( gx(i,1) .gt. 30.e3 .and. gx(i,1) < 90.e3 .and. gy(j,1) > 30.e3 .and. gy(j,1) < 90.e3 ) THEN
!          write(0,*) 'mix: my_rank,i,j =',my_rank,i,j
          tmp = -6.5870481e-8
!          ENDIF
     !     IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmf3sl =',i,j,tmp
          
!          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
!                * ((kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
!                  - ((tmp/Prtl)*gz(k  ,4)) )
          fs(i,j,k) = fs(i,j,k)+gz(k,3)*Prtl                      &
                * ((khv(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))*gz(k+1,4)  &
                   ) - gz(k,3)*tmp*deninv

!          tmp = (kh(i,j,k+1))*(s(i,j,k+1)-s(i,j,k))
!          IF ( i == nx/2 .and. my_rank == 2 ) write(0,*) 'mix: i,j,tmp =',i,j,fs(i,j,k)
         ENDDO
        ENDDO
       
       ENDIF
       
       deallocate( t1 )

      RETURN
      END SUBROUTINE MIX_SCAL_HV
