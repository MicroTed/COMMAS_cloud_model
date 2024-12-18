!
! Originally written by Jerry Straka
! Adapted for COMMAS by ERM 4/2007
!
!
      subroutine stats(          &
         nx,ny,nz,ns             &
        ,time,dx,dy,dz           &
        ,gx,gy,gz                &
        ,t0,sinit,s              &
        ,u,v,w,uinit,vinit,winit &
        ,pi,piinit,km,kminit      &
        ,z1d4                     &
        ,ibsd,iesd,jbsd,jesd)
!

      USE GRID_MODULE
      USE PARAM_MODULE
      USE COMMASMPI_MODULE
      
      implicit none
!
#ifdef MPI
      INCLUDE "mpif.h"
#endif
      integer ns
      
      TYPE(VARIABLE)     :: u, uinit
      TYPE(VARIABLE)     :: v, vinit
      TYPE(VARIABLE)     :: w, winit
      TYPE(VARIABLE)     :: pi, piinit
      TYPE(VARIABLE)     :: km, kminit
      TYPE(VARIABLE)     :: gx(4), gy(4), gz(4)
      TYPE(VARIABLE)     :: s(ns), sinit(2)
      
      real               :: z1d4(nzend,4)
      
      integer ibsd,iesd,jbsd,jesd
      integer, parameter :: igrd=1
      integer, parameter :: istag=1,jstag=1,kstag=1
      integer ic,ip,icnt,ncntmx
      integer jypu,jymu,ixpv,ixmv
      integer jypw,jymw,kzpv,kzmv
      integer kzpu,kzmu,ixpw,ixmw
      integer kzm,kzp
      integer ia,iv,ix,jy,kz,imnmx
      integer nst,nsv,nsp,nss,nsc,nmxmn
!      parameter (nst=50,nsa=75,nsv=3)
!      parameter (nsp=1,nss=14,nsc=10,nlvl=42,nmxmn=4)
      parameter (nst=1) ! (183)
      parameter  (nsv=3)
      parameter (nsp=1,nss=15,nsc=10,nmxmn=4)
      integer ixsmax(nst,nss),jysmax(nst,nss),kzsmax(nst,nss)
      integer ixsmin(nst,nss),jysmin(nst,nss),kzsmin(nst,nss)
      integer ixpmax(nst,nsp),jypmax(nst,nsp),kzpmax(nst,nsp)
      integer ixpmin(nst,nsp),jypmin(nst,nsp),kzpmin(nst,nsp)
      integer ixamax(nst,ns),jyamax(nst,ns),kzamax(nst,ns)
      integer ixamin(nst,ns),jyamin(nst,ns),kzamin(nst,ns)
      integer ixvmax(nst,nsv),jyvmax(nst,nsv),kzvmax(nst,nsv)
      integer ixcmin(nst,nsc),jycmin(nst,nsc),kzcmin(nst,nsc)
      integer ixcmax(nst,nsc),jycmax(nst,nsc),kzcmax(nst,nsc)
      integer ixvmin(nst,nsv),jyvmin(nst,nsv),kzvmin(nst,nsv)
      integer nx,ny,nz
!
      real bvnum,thmn,xval,yval,zval
      real awrk,pwrk,dwrk,akewrk,uwrk,vwrk,wwrk
      real ptim(nst),fac,vor,div,acon,temp
      real vori,vorj,vork,ens
      real divi,divj,divk,dsq
      real facxv,facyv,faczv
      real dx,dy,dz,ake
      integer time
      real cmax(nst,nsc),cmin(nst,nsc),cfac(nsc),vkeavl(nz)
      real smax(nst,nss),smin(nst,nss),sfac(ns),davl(nz)
      real pmax(nst,nsp),pmin(nst,nsp),pfac(nsp),pavl(nz)
      real amax(nst,ns),amin(nst,ns)
      real, allocatable, save :: aavl(:,:) ! (nlvl,nsa)
      real vmax(nst,nsv),vmin(nst,nsv),vfac(nsv),vavl(nz,nsv)
      real pmaxz(nst,nsp,nz),pminz(nst,nsp,nz)
      real, allocatable, save :: amaxz(:,:,:) !(nst,nsa,nlvl)
      real, allocatable, save :: aminz(:,:,:) !(nst,nsa,nlvl)
      real vmaxz(nz,nsv,nst),vminz(nz,nsv,nst)
      real atot(nst,ns),vketot(nst),vave(nst,nsv)
      real, allocatable, save ::  atavl(:,:) ! (nlvl,nsa)
      real trumfx(nz), commfx(nz)
!
      real upup(nz),vpvp(nz),wpwp(nz),tptp(nz),wptp(nz) 
!
!
      real t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!
!      real gt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ngt)
!
!      real ab(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nba)
!      real vb(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nv)
!      real dc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
!      real, allocatable :: pc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real, allocatable :: dc(:,:,:), pc(:,:,:)
      real :: db(nz), pb(nz)
!      real ac(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)
!      real vc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nv)
      
!      real    :: den(0:nz-1)
      real    :: uvw
      integer, parameter :: nv = 3
      integer :: ip1, jp1, kp1
      integer :: n

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
!  integer :: mpi_error_code

   logical :: debug_mpi = .FALSE.

   integer, parameter :: ntot = 50
   real  mpitotin(ntot), mpitotout(ntot)
   double precision, allocatable :: mpitotinth(:,:),mpitotoutth(:,:)
   real, allocatable :: mpitotinthr(:,:),mpitotoutthr(:,:)
      
   double precision  mpitotindp(ntot), mpitotoutdp(ntot)


! ##################################################################

      IF ( .not. allocated( amaxz ) ) THEN
        allocate( amaxz(nst,ns,nzend) )
        allocate( aminz(nst,ns,nzend) )
        allocate( aavl(nzend,ns) )
        allocate( atavl(nzend,ns) ) 
      ENDIF

      allocate ( dc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
      allocate ( pc(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
      
#ifdef MPI
      div = 1./float(max(nxend-1,1)*max(nyend-1,1))
#else
!      div = 1./float(max(nx-1,1)*max(ny-1,1))
      div = 1./float(max(iesd-ibsd+1,1)*max(jesd-jbsd+1,1))
#endif
!      div = 1./float(max(nx-1,1)*max(ny-1,1))

!      den(1:nz-1) = 1.0e5*piinit%flt1d(1:nz-1)**2.509/(rd*sinit(1)%flt1d(1:nz-1) * (1.0+0.61*sinit(2)%flt1d(1:nz-1)))

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
      DO kz = kzb,kze

       db(kz)= 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sinit(1)%flt1d(kz))
       pb(kz)= 1.0e5*piinit%flt1d(kz)**3.509
        
!        write(6,*) kz,db(kz),pb(kz),ab(kz,lt)

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg
        
        jyb = 1
        jye = jtile
        if (jyend .eq. nyend) jye = jyend-jybeg
        
        DO ix=ixb,ixe
          DO jy=jyb,jye
!          IF ( Abs(p2(ix,jy,kz)) .gt. 0.2 ) THEN
!           print*,'ice driver: piinit,p2,an(lt) = ',piinit(kz),p2(ix,jy,kz),an(ix,jy,kz,lt)
!     :        ,u(ix,jy,kz),w(ix,jy,kz)
!          ENDIF
            dc(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/ &
              (287.04*s(1)%flt3d(ix,jy,kz) )
            pc(ix,jy,kz)= 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509 - pb(kz)
          ENDDO
        ENDDO
!        write(6,*) 'pc(nx/2,ny/2,kz),pb = ',p2(nx/2,ny/2,kz),pb(kz),
!     :     db(kz),dc(nx/2,ny/2,kz),an(nx/2,ny/2,kz,lt)
      ENDDO

!
!
!
!  set consts
!
!      icnt = icnt + 1
      icnt = 1
      ptim(icnt) = Float(time)
      pfac(1) = 1.0e-02
      vfac(1) = 1.0
      vfac(2) = 1.0
      vfac(3) = 1.0
!

      aavl(:,:) = 0.0
      atavl(:,:) = 0.0

!#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1
      
      DO kz = kzb,kze
!#else
!      do 01200 kz = 1,nz
!#endif
      davl(kz) = 0.0
      pavl(kz) = 0.0
        DO ia = 1,ns
          aavl(kzbeg-1+kz,ia) = 0.0
          atavl(kzbeg-1+kz,ia) = 0.0
        ENDDO
      
        DO  iv = 1,nv
         vavl(kz,iv) = 0.0
        ENDDO
      vkeavl(kz) = 0.0
      
      ENDDO


!
!     compute max and mins of 'p'
!     also compute max and mins of vorticity and divergence
!
      pmax(icnt,1) = 0.0
      pmin(icnt,1) = 0.0
      ixpmax(icnt,1) = 0
      jypmax(icnt,1) = 0
      kzpmax(icnt,1) = 0
      ixpmin(icnt,1) = 0
      jypmin(icnt,1) = 0
      kzpmin(icnt,1) = 0
      do 71001 ic = 1,nsc 
      cmax(icnt,1) = 0.0
      cmin(icnt,1) = 0.0
      ixcmax(icnt,1) = 0
      jycmax(icnt,1) = 0
      kzcmax(icnt,1) = 0
      ixcmin(icnt,1) = 0
      jycmin(icnt,1) = 0
      kzcmin(icnt,1) = 0
71001 continue

!      IF ( number_of_processes == 1 ) THEN
!      write(0,*) 'u(nz-1),u(nz),uinit(nz) = ',u%flt3d(1,1,nz-1),u%flt3d(1,1,nz),  &
!             uinit%flt1d(nz)
!      ENDIF
             
!
!#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
      do 10001 kz = kzb,kze
!#else
!      do 10001 kz = 1,nz
!#endif
      pmaxz(icnt,01,kz) = -1.0e10
      pminz(icnt,01,kz) =  1.0e10

#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .gt. jesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 10002 jy = jyb,jye
      do 10003 ix = ixb,ixe
#else
!      do 10002 jy = 1,ny-1
!      do 10003 ix = 1,nx-1
      do 10002 jy = jbsd,jesd ! 1,ny-1
      do 10003 ix = ibsd,iesd ! 1,nx-1
#endif
      if ( pc(ix,jy,kz) .gt. pmaxz(icnt,1,kz) ) then
      pmaxz(icnt,1,kz) = pc(ix,jy,kz)
      end if
      if ( pc(ix,jy,kz) .lt. pminz(icnt,1,kz) ) then
      pminz(icnt,1,kz) = pc(ix,jy,kz)
      end if
      if ( pc(ix,jy,kz) .gt. pmax(icnt,1) ) then
      pmax(icnt,1) = pc(ix,jy,kz)
      ixpmax(icnt,1) = ix
      jypmax(icnt,1) = jy
      kzpmax(icnt,1) = kz
      end if
      if ( pc(ix,jy,kz) .lt. pmin(icnt,1) ) then
      pmin(icnt,1) = pc(ix,jy,kz)
      ixpmin(icnt,1) = ix
      jypmin(icnt,1) = jy
      kzpmin(icnt,1) = kz
      end if
      pavl(kz) = pavl(kz) + pc(ix,jy,kz)
      davl(kz) = davl(kz) + dc(ix,jy,kz)
!
!  vorticity, enstrophy and divergence
!
!
!  total tke
!
! #ifdef MPI
      ip1 = ix+istag
      if (ixend .eq. nxend) ip1 = Min(nx,ix+istag)

      jp1 = jy+istag
      if (jyend .eq. nyend) jp1 = Min(ny,jy+jstag)

      kp1 = kz+istag
      if (kzend .eq. nzend) kp1 = Min(nz,kz+kstag)
!#else
!      ip1 = Min(nx,ix+istag)
!      jp1 = Min(ny,jy+istag)
!      kp1 = Min(nz,kz+istag)
!#endif
!      write(0,*) 'stats: ix,jy,kz = ',ix,jy,kz,nx,ny,nz
      uwrk = ((u%flt3d(ix,jy,kz)+u%flt3d(ip1,jy,kz))   &
             -(uinit%flt1d(kz)))
      vwrk = ((v%flt3d(ix,jy,kz)+v%flt3d(ix,jp1,kz))   &
             -(vinit%flt1d(kz)))
      wwrk = ((w%flt3d(ix,jy,kz)+w%flt3d(ix,jy,kp1))   &
             -(winit%flt1d(kz)))
      ake = 0.125*(uwrk**2 + vwrk**2 + wwrk**2)
      vkeavl(kz) = vkeavl(kz) + ake
      vketot(icnt) = vketot(icnt) + ake
!
10003 continue
10002 continue
10001 continue
!
!     compute max and mins of 'a'
!
      amaxz(icnt,:,:) = 0.0
      aminz(icnt,:,:) = 0.0

      do 20000 ia = 1,ns
      amax(icnt,ia) = 0.0
      amin(icnt,ia) = 0.0
      ixamax(icnt,ia) = 0
      jyamax(icnt,ia) = 0
      kzamax(icnt,ia) = 0
      ixamin(icnt,ia) = 0
      jyamin(icnt,ia) = 0
      kzamin(icnt,ia) = 0
!#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag
      
      do 20001 kz = kzb,kze
!#else
!      do 20001 kz = 1,nz-kstag
!#endif
      amaxz(icnt,ia,kzbeg-1+kz) = -1.0e10
      aminz(icnt,ia,kzbeg-1+kz) =  1.0e10

#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .gt. jesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 20002 jy = jyb,jye
      do 20003 ix = ixb,ixe
#else
!      do 20002 jy = 1,ny-jstag
!      do 20003 ix = 1,nx-istag
      do 20002 jy = jbsd,jesd ! 1,ny-jstag
      do 20003 ix = ibsd,iesd ! 1,nx-istag
#endif
       IF ( s(ia)%name == 'TH' ) THEN
        awrk = s(ia)%flt3d(ix,jy,kz)-sinit(ia)%flt1d(kz)
       ELSE 
        awrk = s(ia)%flt3d(ix,jy,kz)
       ENDIF
      if ( s(ia)%flt3d(ix,jy,kz) .gt. amaxz(icnt,ia,kzbeg-1+kz) ) then
      amaxz(icnt,ia,kzbeg-1+kz) = s(ia)%flt3d(ix,jy,kz)
      end if
      if ( s(ia)%flt3d(ix,jy,kz) .lt. aminz(icnt,ia,kzbeg-1+kz) ) then
      aminz(icnt,ia,kzbeg-1+kz) = s(ia)%flt3d(ix,jy,kz)
      end if 
      if ( awrk .gt. amax(icnt,ia) ) then
      amax(icnt,ia) = awrk
      ixamax(icnt,ia) = ix
      jyamax(icnt,ia) = jy
      kzamax(icnt,ia) = kz
      end if
      if ( awrk .lt. amin(icnt,ia) ) then
      amin(icnt,ia) = awrk
      ixamin(icnt,ia) = ix
      jyamin(icnt,ia) = jy
      kzamin(icnt,ia) = kz
      end if
      if ( s(ia)%name == 'TH' .or. s(ia)%name == 'QV'  ) then
        aavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia) + awrk
        atavl(kzbeg-1+kz,ia) = atavl(kzbeg-1+kz,ia) + s(ia)%flt3d(ix,jy,kz)
        atot(icnt,ia) = atot(icnt,ia) + awrk
!      end if
      elseif ( s(ia)%name == 'CPION' .or. s(ia)%name == 'CNION') THEN
        aavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)    &
          +    s(ia)%flt3d(ix,jy,kz) *div
        atavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)
        atot(icnt,ia) = atot(icnt,ia)    &
           + s(ia)%flt3d(ix,jy,kz)
      elseif (  s(ia)%name(1:1) == 'Q' ) then
        aavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)    &
          + max( dc(ix,jy,kz)*s(ia)%flt3d(ix,jy,kz), 0.0 ) 
!          + max( s(ia)%flt3d(ix,jy,kz), 0.0 ) 
        atavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)
        atot(icnt,ia) = atot(icnt,ia)    &
          + max( dc(ix,jy,kz)*s(ia)%flt3d(ix,jy,kz), 0.0 ) 
      else ! if ( ia .gt. 15 ) THEN
        aavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)    &
         +    s(ia)%flt3d(ix,jy,kz)  
        atavl(kzbeg-1+kz,ia) = aavl(kzbeg-1+kz,ia)
        atot(icnt,ia) = atot(icnt,ia)    &
         + s(ia)%flt3d(ix,jy,kz)
      end if
20003 continue
20002 continue
20001 continue
20000 continue
!
!     compute max and mins of 'v'
!
      do 30000 iv = 1,nv
      vmax(icnt,iv) = 0.0
      vmin(icnt,iv) = 0.0
      ixvmax(icnt,iv) = 0
      jyvmax(icnt,iv) = 0
      kzvmax(icnt,iv) = 0
      ixvmin(icnt,iv) = 0
      jyvmin(icnt,iv) = 0
      kzvmin(icnt,iv) = 0

!#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg
      
      do 30001 kz = kzb,kze
!#else
!      do 30001 kz = 1,nz
!#endif
      vmaxz(kz,iv,icnt) = -1.0e10
      vminz(kz,iv,icnt) =  1.0e10

#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .gt. jesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 30002 jy = jyb,jye
      do 30003 ix = ixb,ixe
#else
!      do 30002 jy = 1,ny
!      do 30003 ix = 1,nx
      do 30002 jy = jbsd,jesd ! 1,ny
      do 30003 ix = ibsd,iesd ! 1,nx
#endif
      IF ( iv .eq. 1 ) THEN
         uvw = u%flt3d(ix,jy,kz)
         vwrk = u%flt3d(ix,jy,kz) - uinit%flt1d(kz)
      ELSEIF ( iv .eq. 2 ) THEN
         uvw = v%flt3d(ix,jy,kz)
         vwrk = v%flt3d(ix,jy,kz) - vinit%flt1d(kz)
      ELSEIF ( iv .eq. 3 ) THEN
         uvw = w%flt3d(ix,jy,kz)
         vwrk = w%flt3d(ix,jy,kz) - winit%flt1d(kz)
      ENDIF
!      vwrk = vc(ix,jy,kz,iv) - vb(00,00,kz,iv)
      if ( vwrk .gt. vmaxz(kz,iv,icnt) ) then
      vmaxz(kz,iv,icnt) = uvw ! vc(ix,jy,kz,iv)
      end if
      if ( vwrk .lt. vminz(kz,iv,icnt) ) then
      vminz(kz,iv,icnt) = uvw ! vc(ix,jy,kz,iv)
      end if 
      if ( vwrk .gt. vmax(icnt,iv) ) then
      vmax(icnt,iv) = vwrk
      ixvmax(icnt,iv) = ix
      jyvmax(icnt,iv) = jy
      kzvmax(icnt,iv) = kz
      end if
      if ( vwrk .lt. vmin(icnt,iv) ) then
      vmin(icnt,iv) = vwrk
      ixvmin(icnt,iv) = ix
      jyvmin(icnt,iv) = jy
      kzvmin(icnt,iv) = kz
      end if
      vave(icnt,iv) = vave(icnt,iv) + vwrk
      vavl(kz,iv) = vavl(kz,iv) + vwrk
30003 continue
30002 continue
30001 continue

#ifdef MPI

! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
       ip = iv
        mpitotin(:) = 0
      
        mpitotin(1) = vmax(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        vmax(icnt,ip) = mpitotout(1)
      ENDIF

        mpitotin(1) = vmin(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        vmin(icnt,ip) = mpitotout(1)
      ENDIF

      ENDIF

#endif


30000 continue
!
!
#ifdef MPI
! communicate wind component and pressure info
      IF ( number_of_processes .gt. 1 ) THEN
        allocate(  mpitotinth(nzend,nsv) )
        allocate( mpitotoutth(nzend,nsv) )
       
! first the winds
        mpitotinth(:,:) = 0
        mpitotinth(kzbeg:kzend,1:nsv) = vavl(kzbeg:kzend,1:nsv) ! atavl(kzbeg:kzend,1:ns)
      
       n = nzend*nsv

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        vavl(1:nzend,1:nsv) = mpitotoutth(1:nzend,1:nsv)
      ENDIF

      mpitotinth(:,:) = 0
      DO ia = 1,nsv
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = vmaxz(kz,ia,icnt)
        ENDDO
      ENDDO

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,nsv
        DO kz = 1,nzend
          vmaxz(kz,ia,icnt) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

      mpitotinth(:,:) = 1.e20
      DO ia = 1,nsv
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = vminz(kz,ia,icnt)
        ENDDO
      ENDDO

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,nsv
        DO kz = 1,nzend
          vminz(kz,ia,icnt) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

! now do the pressure
      mpitotinth(:,:) = 0
      DO ia = 1,1
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = pmaxz(icnt,1,kz)
        ENDDO
      ENDDO
       n = nzend*1

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,1
        DO kz = 1,nzend
          pmaxz(icnt,1,kz) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

      mpitotinth(:,:) = 1.e20
      DO ia = 1,1
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = pminz(icnt,1,kz)
        ENDDO
      ENDDO

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,1
        DO kz = 1,nzend
          pminz(icnt,1,kz) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

       
        mpitotin(:) = 0
      
        mpitotin(1) = pmax(icnt,1)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        pmax(icnt,1) = mpitotout(1)
      ENDIF

        mpitotin(1) = pmin(icnt,1)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        pmin(icnt,1) = mpitotout(1)
      ENDIF


        deallocate( mpitotinth )
        deallocate( mpitotoutth )

      ENDIF

#endif


!

      IF ( my_rank == 0 ) THEN
      write(3,*)    &
       '======================================================'
      write(3,*) 'TIME HISTORY STATISTICS: GROUP I'
      write(3,*)   &
       '======================================================'
      ip = 1
      write(3,*)
      write(3,*)   &
        'Max and Mins for pert pressure (mb) on grid: ', igrd
      write(3,*)   &
        '______________________________________________________'
      write(3,90101)    &
        ptim(icnt),pfac(1)*pmax(icnt,1),ip ! ,ixpmax(icnt,1)   &
!        ,jypmax(icnt,1),kzpmax(icnt,1)
      write(3,90102)    &
        ptim(icnt),pfac(1)*pmin(icnt,1),ip,ixpmin(icnt,1)   &
        ,jypmin(icnt,1),kzpmin(icnt,1)
      write(3,*)
      write(3,*)   &
        'Layer average, min & max for pert pressure (mb) on grid: ', igrd
      ENDIF
#ifdef MPI
      div = 1./float(max(nxend-1,1)*max(nyend-1,1))
#else
      div = 1./float(max(nx-1,1)*max(ny-1,1))
#endif
!#ifdef MPI
      kzb = nz-1 ! should be 'nzend' but need to fix pavl to be nzend ! ktile
      kze = 1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg+1
      
      do 91001 kz = kzb,kze,-1
!#else
!      do 91001 kz = nz,1,-1
!#endif
!      acon = gz(1)%flt1d(kz) ! dz*(kz-1)+(0.5)*dz
      acon = z1d4(kz,1)
      IF ( my_rank == 0 ) THEN
      write(3,90103)    &
        ptim(icnt),kz,acon,pfac(1)*pavl(kz)*div,pfac(1)*pminz(icnt,01,kz)   &
       ,pfac(1)*pmaxz(icnt,01,kz), pb(kz)
      ENDIF
91001 continue
90101 format(1x,'pmax  ',1f10.2,1f12.6,i6)
! 90101 format(1x,'pmax  ',1f10.2,1f12.6,4i6)
90102 format(1x,'pmin  ',1f10.2,1f12.6,4i6)
90103 format(1x,'p(anx)',1f10.2,1x,i6,1x,1f10.2,4(1x,f12.3))
!
      IF ( my_rank == 0 ) THEN
      write(3,*)   &
        '______________________________________________________'
      write(3,*)
      write(3,*)   &
        'Max and Mins for grid-relative velocity (m/s) components on grid: '
      write(3,*)   &
        '______________________________________________________'
      do iv = 1,nv
      write(3,90201)    &
        ptim(icnt),vfac(iv)*vmax(icnt,iv),iv,ixvmax(icnt,iv)   &
        ,jyvmax(icnt,iv),kzvmax(icnt,iv)
      write(3,90202)    &
        ptim(icnt),vfac(iv)*vmin(icnt,iv),iv,ixvmin(icnt,iv)   &
        ,jyvmin(icnt,iv),kzvmin(icnt,iv) 
      write(3,*)   &
        '______________________________________________________'

      ENDDO
      ENDIF
      
      do 81001 iv = 1,nv

      IF ( my_rank == 0 ) THEN

      write(3,*)
      write(3,*)   &
        'Layer average, min & max for velocity (m/s) on grid: ',igrd 
      ENDIF
#ifdef MPI
      if ( iv .eq. 1 ) div = 1./float(max(nxend,1)*max(nyend-1,1))
      if ( iv .eq. 2 ) div = 1./float(max(nxend-1,1)*max(nyend,1))
      if ( iv .eq. 3 ) div = 1./float(max(nxend-1,1)*max(nyend-1,1))
!      div = 1./float(max(iesd-ibsd+1,1)*max(jesd-jbsd+1,1))
#else
!      if ( iv .eq. 1 ) div = 1./float(max(nx,1)*max(ny-1,1))
!      if ( iv .eq. 2 ) div = 1./float(max(nx-1,1)*max(ny,1))
!      if ( iv .eq. 3 ) div = 1./float(max(nx-1,1)*max(ny-1,1))
      div = 1./float(max(iesd-ibsd+1,1)*max(jesd-jbsd+1,1))
#endif
!#ifdef MPI
      kzb = nz-1 ! +Max(0,nv-2) ! ! should be 'nzend' but need to fix vavl to be nzend ! ktile
      kze = 1
!      if (kzend .eq. nzend) kzb = kzend-kzbeg+1
      
      do 81101 kz = kzb,kze,-1
!#else
!      do 81101 kz = nz,1,-1
!#endif
!      if ( iv .eq. 1 ) acon = dz*(kz-1)+.5*dz
!      if ( iv .eq. 2 ) acon = dz*(kz-1)+.5*dz
!      if ( iv .eq. 3 ) acon = dz*(kz-1)
      if ( iv .eq. 1 ) acon = z1d4(kz,1) ! gz(1)%flt1d(kz)
      if ( iv .eq. 2 ) acon = z1d4(kz,1) ! gz(1)%flt1d(kz)
      if ( iv .eq. 3 ) acon = z1d4(kz,2) ! gz(2)%flt1d(kz) 

      IF ( my_rank == 0 ) THEN
      write(3,90203)    &
        ptim(icnt),iv,kz,acon,vfac(iv)*vavl(kz,iv)*div,vminz(kz,iv,icnt)   &
       ,vmaxz(kz,iv,icnt)
      ENDIF
81101 continue
      IF ( my_rank == 0 ) THEN
      write(3,*)   &
        '______________________________________________________'
      ENDIF
81001 continue
      IF ( my_rank == 0 ) THEN
      write(3,*)
      write(3,*)   &
        'Averages for resolved TKE  ',igrd
      ENDIF
#ifdef MPI
      kzb = ktile
      kze = 1
      if (kzend .eq. nzend) kzb = kzend-kzbeg
      
      do 82001 kz = kzb,kze,-1
#else
      do 82001 kz = nz,1,-1
#endif
!      acon = dz*(kz-1)+.5*dz
      acon = z1d4(kz,1) ! gz(1)%flt1d(kz)
      IF ( my_rank == 0 ) THEN
      write(3,90204)    &
        ptim(icnt),kz,acon,vkeavl(kz)
      ENDIF
82001 continue
      IF ( my_rank == 0 ) THEN
      write(3,*)   &
        '______________________________________________________'
      ENDIF
90201 format(1x,'vmax  ',1f10.2,1f12.6,4i6)
90202 format(1x,'vmin  ',1f10.2,1f12.6,4i6)
90203 format(1x,'v(anx)',1f10.2,2i6,1f10.2,3f12.6)
90204 format(1x,'kavl  ',1f10.2,1i4,1f10.2,1f16.6)
!
!
!
!
!
      IF ( my_rank == 0 ) THEN

      IF ( ns .ge. 100 ) THEN
      write(3,*)
      write(3,*) 'UNITS FOR SCALARS TAKAHASHI'
      write(3,*) 'Scalar 1, potential temperature: K'

      ELSEIF ( ns .ge. 40 ) THEN
      write(3,*)
      write(3,*) 'UNITS FOR SCALARS'
      write(3,*) 'Scalar 1, potential temperature: K'
      write(3,*) 'Scalar 2, TKE: m2 s-2'
      write(3,*) 'Scalar =>3, mixing ratios: kg/kg'
      write(3,*) '  qv=3, qc=4, qr=5, qi=6, qir=7 qs=8, qgl=9 qgm=10'
      write(3,*) '  qgh=11, qf=12, qh=13, qip=14, qhl=15 '
      write(3,*) '  lscw=17 lschl=28'
      ELSEIF ( ns .gt. 15 ) THEN
      write(3,*)
      write(3,*) 'UNITS FOR SCALARS'
      write(3,*) 'Scalar 1, potential temperature: K'
      write(3,*) 'Scalar 2, TKE: m2 s-2'
      write(3,*) 'Scalar =>3, mixing ratios: kg/kg'
      write(3,*) '  qv=3, qc=4, qr=5, qi=6, qs=7 qh=8'
      write(3,*) ' lscw=9, lscr=10, lsci=11, lscs=12 lsch=13'
      write(3,*) ' lscpi=14, lscni=15'
      write(3,*) ' lncw=14, lni=15'      
      ELSEIF ( ns .eq. 15 ) THEN
      write(3,*)
      write(3,*) 'UNITS FOR SCALARS'
      write(3,*) 'Scalar 1, potential temperature: K'
      write(3,*) 'Scalar 2, TKE: m2 s-2'
      write(3,*) 'Scalar =>3, mixing ratios: kg/kg'
      write(3,*) '  qv=3, qc=4, qr=5, qi=6, qs=7 qh=8'
      write(3,*) ' lscw=9, lscr=10, lsci=11, lscs=12 lsch=13'
      write(3,*) ' lscpi=14, lscni=15'
      ELSE
      write(3,*)
      write(3,*) 'UNITS FOR SCALARS'
      write(3,*) 'Scalar 1, potential temperature: K'
      write(3,*) 'Scalar 2, TKE: m2 s-2'
      write(3,*) 'Scalar =>3, mixing ratios: kg/kg'
      write(3,*) '  qv=3, qc=4, qr=5, qi=6, qs=7 qh=8'
      ENDIF
       write(3,*)    &
       'Max and Mins for scalars on grid: ', igrd
       write(3,*)    &
       '______________________________________________________'

! don't care about max/mins for MPI because they are already output in print_mpi
      do  ia = 1,ns
      write(3,90301)    &
        ptim(icnt),amax(icnt,ia),ia,ixamax(icnt,ia)   &
        ,jyamax(icnt,ia),kzamax(icnt,ia), s(ia)%name
      write(3,90302)    &
        ptim(icnt),amin(icnt,ia),ia,ixamin(icnt,ia)   &
        ,jyamin(icnt,ia),kzamin(icnt,ia)
      write(3,*)   &
        '______________________________________________________'
      ENDDO

      write(3,*)
      ENDIF

#ifdef MPI
! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
        allocate(  mpitotinth(nzend,ns) )
        allocate( mpitotoutth(nzend,ns) )
       
        mpitotinth(:,:) = 0
        mpitotinth(kzbeg:kzend,1:ns) = atavl(kzbeg:kzend,1:ns)
      
       n = nzend*ns

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        atavl(1:nzend,1:ns) = mpitotoutth(1:nzend,1:ns)
      ENDIF

      mpitotinth(:,:) = 0
      DO ia = 1,ns
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = amaxz(icnt,ia,kz)
        ENDDO
      ENDDO

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,ns
        DO kz = 1,nzend
          amaxz(icnt,ia,kz) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

      mpitotinth(:,:) = 1.e20
      DO ia = 1,ns
        DO kz = kzbeg,kzend
          mpitotinth(kz,ia) = aminz(icnt,ia,kz)
        ENDDO
      ENDDO

      CALL MPI_Reduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_MIN, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
      DO ia = 1,ns
        DO kz = 1,nzend
          aminz(icnt,ia,kz) = mpitotoutth(kz,ia)
        ENDDO
      ENDDO
      ENDIF

        deallocate( mpitotinth )
        deallocate( mpitotoutth )

      ENDIF

#endif
      do 83001 ia = 1,ns
      IF ( my_rank == 0 ) THEN
      write(3,*)
      write(3,*)   &
        'Layer average, max & min for scalars ',s(ia)%name,' on grid: ',igrd
      ENDIF
#ifdef MPI
      div = 1./float(max(nxend-1,1)*max(nyend-1,1))
#else
!      div = 1./float(max(nx-1,1)*max(ny-1,1))
      div = 1./float(max(iesd-ibsd+1,1)*max(jesd-jbsd+1,1))
#endif
#ifdef MPI
      kze = nzend
      kzb = 1
!      if (kzend .eq. nzend) kze = kzend-kzbeg+1
      
      do 83101 kz = kze,kzb,-1
#else
      do 83101 kz = nz,1,-1
#endif
!      acon = dz*(kz-1)+.5*dz
!      acon = gz(1)%flt1d(kz)
      acon = z1d4(kz,1)
      IF ( my_rank == 0 ) THEN
      IF ( .not. ( s(ia)%name == 'CPION' .or. s(ia)%name == 'CNION') ) THEN
      write(3,90303)    &
        ptim(icnt),ia,kz,acon,   &
        atavl(kz,ia)*div   &
       ,amaxz(icnt,ia,kz)   &
       ,aminz(icnt,ia,kz)
      ELSE
      write(3,90303)    &
        ptim(icnt),ia,kz,acon,   &
        atavl(kz,ia)   &
       ,amaxz(icnt,ia,kz)   &
       ,aminz(icnt,ia,kz)
      ENDIF
      ENDIF
83101 continue
      IF ( my_rank == 0 ) THEN
      write(3,*)   &
        '______________________________________________________'
      ENDIF
83001 continue
90301 format(1x,'amax  ',1f10.2,1x,1pe14.7,4i6,2x,a)
90302 format(1x,'amin  ',1f10.2,1x,1pe14.7,4i6)
90303 format(1x,'a(anx)',1f10.2,2i5,1f10.2,3(1x,1pe12.5))
!
!  special values
!
!
!     compute max and mins of 'theta' at first level
!
      smax(icnt,1) = -1.0e-10
      smin(icnt,1) =  1.0e-10
      ixsmax(icnt,1) = 0
      jysmax(icnt,1) = 0
      kzsmax(icnt,1) = 0
      ixsmin(icnt,1) = 0
      jysmin(icnt,1) = 0
      kzsmin(icnt,1) = 0
      do 40000 ia = 1,1
#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      kzb = 1
      kze = 1
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 40001 jy = jyb,jye
      do 40002 kz = kzb,kze
      do 40003 ix = ixb,ixe
#else
      do 40001 jy = jbsd,jesd ! 1,ny-1
      do 40002 kz = 1,1
      do 40003 ix = ibsd,iesd ! 1,nx-1
#endif
      if ( s(ia)%flt3d(ix,jy,kz)-sinit(ia)%flt1d(kz) .gt. smax(icnt,1) ) then
      smax(icnt,1) = s(ia)%flt3d(ix,jy,kz)-sinit(ia)%flt1d(kz)
      ixsmax(icnt,1) = ix
      jysmax(icnt,1) = jy
      kzsmax(icnt,1) = kz
      end if
      if ( s(ia)%flt3d(ix,jy,kz)-sinit(ia)%flt1d(kz)  &
       .lt. smin(icnt,1) ) then
      smin(icnt,1) = s(ia)%flt3d(ix,jy,kz)-sinit(ia)%flt1d(kz)
      ixsmin(icnt,1) = ix
      jysmin(icnt,1) = jy
      kzsmin(icnt,1) = kz
      end if
40003 continue
40002 continue
40001 continue              
40000 continue              
      ip = 1
#ifdef MPI
! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
       
        mpitotin(:) = 0
      
        mpitotin(1) = smax(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        smax(icnt,ip) = mpitotout(1)
      ENDIF

        mpitotin(1) = smin(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        smin(icnt,ip) = mpitotout(1)
      ENDIF

      ENDIF

#endif
      IF ( my_rank == 0 ) THEN
      write(3,*)
      write(3,*)   &
        'Max and Mins for theta prime at level 1: '
      write(3,*)   &
        '______________________________________________________'
      write(3,94001)    &
        ptim(icnt),smax(icnt,1),ip,ixsmax(icnt,1)   &
        ,jysmax(icnt,1),kzsmax(icnt,1)
      write(3,94002)    &
        ptim(icnt),smin(icnt,1),ip,ixsmin(icnt,1)   &
        ,jysmin(icnt,1),kzsmin(icnt,1)
      write(3,*)   &
        '______________________________________________________'
      ENDIF
94001 format(1x,'tmax  ',f11.3,1x,f12.6,4i6)
94002 format(1x,'tmin  ',f11.3,1x,f12.6,4i6)
!
!     compute max and mins of 'u' at level 1
!
      smax(icnt,2) = -1.0e-10
      smin(icnt,2) =  1.0e-10
      ixsmax(icnt,2) = 0
      jysmax(icnt,2) = 0
      kzsmax(icnt,2) = 0
      ixsmin(icnt,2) = 0
      jysmin(icnt,2) = 0
      kzsmin(icnt,2) = 0
      do 40010 iv = 01,01
#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      kzb = 1
      kze = 1
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 40011 jy = jyb,jye
      do 40012 kz = kzb,kze
      do 40013 ix = ixb,ixe
#else
      do 40011 jy = jbsd,jesd ! 1,ny-1
      do 40012 kz = 1,1
      do 40013 ix = ibsd,iesd ! 1,nx
#endif
      vwrk = u%flt3d(ix,jy,kz)-uinit%flt1d(kz)
      if ( vwrk .gt. smax(icnt,2) ) then
      smax(icnt,2) = vwrk
      ixsmax(icnt,2) = ix
      jysmax(icnt,2) = jy
      kzsmax(icnt,2) = kz
      end if
      if ( vwrk .lt. smin(icnt,2) ) then
      smin(icnt,2) = vwrk
      ixsmin(icnt,2) = ix
      jysmin(icnt,2) = jy
      kzsmin(icnt,2) = kz
      end if
40013 continue
40012 continue
40011 continue
40010 continue
      ip = 2
#ifdef MPI
! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
       
        mpitotin(:) = 0
      
        mpitotin(1) = smax(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        smax(icnt,ip) = mpitotout(1)
      ENDIF

        mpitotin(1) = smin(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        smin(icnt,ip) = mpitotout(1)
      ENDIF

      ENDIF

#endif
      IF ( my_rank == 0 ) THEN
      write(3,*)
      write(3,*)   &
        'Max and Mins for u at level 1:'
      write(3,*)   &
        '______________________________________________________'
      write(3,94011)   &
        ptim(icnt),smax(icnt,2),ip,ixsmax(icnt,2)   &
        ,jysmax(icnt,2),kzsmax(icnt,2)
      write(3,94012)   &
        ptim(icnt),smin(icnt,2),ip,ixsmin(icnt,2)   &
        ,jysmin(icnt,2),kzsmin(icnt,2)
      write(3,*)   &
        '______________________________________________________'
      ENDIF
94011 format(1x,'u-wind max  ',f11.3,1x,f12.6,4i6)
94012 format(1x,'u-wind min  ',f11.3,1x,f12.6,4i6)     
!
!     compute max and mins of 'v' at level 1
!
      smax(icnt,3) = -1.0e-10
      smin(icnt,3) =  1.0e-10
      ixsmax(icnt,3) = 0
      jysmax(icnt,3) = 0
      kzsmax(icnt,3) = 0
      ixsmin(icnt,3) = 0
      jysmin(icnt,3) = 0
      kzsmin(icnt,3) = 0
      do 40020 iv = 02,02
#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      kzb = 1
      kze = 1
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 40021 jy = jyb,jye
      do 40022 kz = kzb,kze
      do 40023 ix = ixb,ixe
#else
      do 40021 jy = jbsd,jesd ! 1,ny
      do 40022 kz = 1,1
      do 40023 ix = ibsd,iesd ! 1,nx-1
#endif
      vwrk = v%flt3d(ix,jy,kz)-vinit%flt1d(kz)
      if ( vwrk .gt. smax(icnt,3) ) then
      smax(icnt,3) = vwrk
      ixsmax(icnt,3) = ix
      jysmax(icnt,3) = jy
      kzsmax(icnt,3) = kz
      end if
      if ( vwrk .lt. smin(icnt,3) ) then
      smin(icnt,3) = vwrk
      ixsmin(icnt,3) = ix
      jysmin(icnt,3) = jy
      kzsmin(icnt,3) = kz
      end if
40023 continue
40022 continue
40021 continue
40020 continue
      ip = 3
#ifdef MPI
! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
       
        mpitotin(:) = 0
      
        mpitotin(1) = smax(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        smax(icnt,ip) = mpitotout(1)
      ENDIF

        mpitotin(1) = smin(icnt,ip)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        smin(icnt,ip) = mpitotout(1)
      ENDIF

      ENDIF

#endif
      IF ( my_rank == 0 ) THEN
      write(3,*)
      write(3,*)   &
        'Max and Mins for v at level 1:'
      write(3,*)   &
        '______________________________________________________'
      write(3,94021)   &
        ptim(icnt),smax(icnt,3),ip,ixsmax(icnt,3)   &
        ,jysmax(icnt,3),kzsmax(icnt,3)
      write(3,94022)   &
        ptim(icnt),smin(icnt,3),ip,ixsmin(icnt,3)   &
        ,jysmin(icnt,3),kzsmin(icnt,3)
      write(3,*)   &
        '______________________________________________________'
      ENDIF
94021 format(1x,'v-wind max  ',f11.3,1x,f12.6,4i6)
94022 format(1x,'v-wind min  ',f11.3,1x,f12.6,4i6)     
!
!     compute max and mins of 'w' at level 1
!
      smax(icnt,4) = -1.0e-10
      smin(icnt,4) =  1.0e-10
      ixsmax(icnt,4) = 0
      jysmax(icnt,4) = 0
      kzsmax(icnt,4) = 0
      ixsmin(icnt,4) = 0
      jysmin(icnt,4) = 0
      kzsmin(icnt,4) = 0
      do 40030 iv = 03,03
#ifdef MPI
      jyb = 1
      jye = jtile
!      if (jybeg .le. jbsd) jyb = jbsd
!      if (jyend .ge. jesd) jye = jesd-jybeg+1
      IF ( jyend .eq. nyend ) jye = jyend - jybeg
      
      kzb = 1
      kze = 2
      
      ixb = 1
      ixe = itile
!      if (ixbeg .le. ibsd) ixb = ibsd
!      if (ixend .ge. iesd) ixe = iesd-ixbeg+1
      IF ( ixend .eq. nxend ) ixe = ixend-ixbeg
      
      do 40031 jy = jyb,jye
      do 40032 kz = kzb,kze
      do 40033 ix = ixb,ixe
#else
      do 40031 jy = jbsd,jesd ! 1,ny-1
      do 40032 kz = 1,2
      do 40033 ix = ibsd,iesd ! 1,nx-1
#endif
      if ( w%flt3d(ix,jy,kz) .gt. smax(icnt,4) ) then
      smax(icnt,4) = w%flt3d(ix,jy,kz)
      ixsmax(icnt,4) = ix
      jysmax(icnt,4) = jy
      kzsmax(icnt,4) = kz
      end if
      if ( w%flt3d(ix,jy,kz) .lt. smin(icnt,4) ) then
      smin(icnt,4) = w%flt3d(ix,jy,kz)
      ixsmin(icnt,4) = ix
      jysmin(icnt,4) = jy
      kzsmin(icnt,4) = kz
      end if
40033 continue
40032 continue
40031 continue
40030 continue
      ip = 4
#ifdef MPI
! communicate scalar info
!        atavl(kz,ia)*div   &
!       ,amaxz(icnt,ia,kz)   &
!       ,aminz(icnt,ia,kz)
      IF ( number_of_processes .gt. 1 ) THEN
       
        mpitotin(:) = 0
      
        mpitotin(1) = smax(icnt,4)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
      
      IF ( my_rank == 0 ) THEN
        smax(icnt,4) = mpitotout(1)
      ENDIF

        mpitotin(1) = smin(icnt,4)
        
        n = 1

      CALL MPI_AllReduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
        smin(icnt,4) = mpitotout(1)
      ENDIF

      ENDIF

#endif
      IF ( my_rank == 0 ) THEN ! need to fix for MPI
      write(3,*)
      write(3,*)   &
        'Max and Mins for w at level 1:'
      write(3,*)   &
        '______________________________________________________'
      write(3,94031)   &
        ptim(icnt),smax(icnt,4),ip ! ,ixsmax(icnt,4)   &
!        ,jysmax(icnt,4),kzsmax(icnt,4)
      write(3,94032)   &
        ptim(icnt),smin(icnt,4),ip ! ,ixsmin(icnt,4)   &
!        ,jysmin(icnt,4),kzsmin(icnt,4)
      write(3,*)   &
        '______________________________________________________'
      ENDIF
94031 format(1x,'wmax  ',f11.3,1x,f12.6,i6)
94032 format(1x,'wmin  ',f11.3,1x,f12.6,i6)     
! 94031 format(1x,'wmax  ',f11.3,1x,f12.6,4i6)
! 94032 format(1x,'wmin  ',f11.3,1x,f12.6,4i6)     
!
!
      deallocate ( dc )
      deallocate ( pc )
      
      RETURN
      END
