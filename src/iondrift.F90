!2345678901234567890123456789012345678901234567890123456789012345678912
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  STRAKA'S ATMOSPHERIC MODEL  (SAM)
!    Designed by Jerry M. Straka
!
!  SUBROUTINE IONDRIFT
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!
!  Crowley advection scheme adapted for ions
!
!  call separately for pos,neg ions with respective mobilities (uz)
!  
!
! think about boundaries (vnormx)... and how efield drives ions.
!  maybe make uzn negative to make uz*elec positive (like velocity)?
!
      subroutine iondrift &
       (nx,ny,nz,nor,na,ia,s,emaxc,icorona,nstep, &
        dt,dx,dy,dz,cion,scionfx, &
        ab,ad,an,                 & ! ad is t7 or t8 &
        elec,uz,ezfair,           &
!        ex,ey,ez,uz,ezfair,           &
        t0,ezw,                   & ! ezw is t9 &
        fu,fv,fw,                 & ! t1, t2, t3 &
        dxx,dyy,dzz,gx,gy,gz,iestag,bcx,bcy) 
     
       USE GRID_MODULE
       USE INDEX_MODULE
       USE COMMASMPI_MODULE
!
!  Revisions:
!   2.17.2003  Now using ezw(ix,jy,1) to calculate corona current
!              because it is more responsive to feedback of the 
!              new corona ions.
!
!
!  include file :  contains implicit none
!
      implicit none
!      include 'sam.param.h'
!      include 'sam.index.f'
!
!  declare all variables

      real dxx(nx),dyy(ny),dzz(nz)         ! dx(i),dy(j),dz(k)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)

      integer ng1, nor, nstep, ng
      parameter(ng1 = 1)

!      real ab(nz,na)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)

!
!      real pb(nz)
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
!      real db(nz)
      
!      real ex(nx,ny,nz)
!      real ey(nx,ny,nz)
!      real ez(nx,ny,nz)

!

      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ad(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ezw(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real fu(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fw(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      integer bcx, bcy

      real fus,fvs,fws

      double precision t0s,t1s,t2s,t3s,t4s,t5s,t6s

      integer imn,imx
      integer jmn,jmx
      integer kmn,kmx
!
!  INTEGERS
!
      integer    ia,na
      integer    ip
      integer    ipass
      integer    i,j,k
      integer    ix
      integer    jy
      integer    kz
      integer    mdnstp
      integer    n1
      integer    n2
      integer    np
      integer    npass
      integer    nx
      integer    ny
      integer    nz
      parameter (np=6)
      
      integer icorona
      
      integer ndebug
      parameter (ndebug=0)
      
      integer iestag
!
!  REALS
!
      double precision scionfx(6)
      
      real dv
      
      real delc,ez1
      
      real jc
      parameter ( jc = -2.0e-12 )
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)

      real eperao
      parameter (eperao  = 8.8592e-12 )
      
      real alpha ! rate for forcing ions at ground back to fair weather value
      parameter (alpha = 0.0577/5.0)

      real emaxc ! point discharge (corona) threshold

      double precision qu, qd, qc, qr, qdel, qxs, qxp, qsr, qsl
      
      real  s   ! sign of ion charge
      real ezfair(nz)
      real cion(nz)  !  base state ion density (pos or neg depending)
!
      real       dt
      real       dx
      real       dy
      real       dz
      real       dxi
      real       dyi
      real       dzi
      real       dxt
      real       dyt
      real       dzt
      real       dtdx
      real       dtdy
      real       dtdz
      real       dxt4
      real       dyt4
      real       dzt4
      real       dtdx4
      real       dtdy4
      real       dtdz4
      real       dtfac
      real       vnorme
      real       vnormw
      real       vnormn
      real       vnorms
      real vnormb,vnormt
      real      tmp
!
!      real       cx(nxl,nzl,np),  
      real  px(-nor+1:nx+nor,np,np),py(-nor+1:ny+nor,np,np),pz(-nor+1:nz+nor,np,np)
!      real       cy(nyl,nzl,np),  py(nyl,np,np)
!      real       cz(nxl,nzl,np),  pz(nzl,np,np)
!
!

      TYPE(VARIABLE), INTENT(IN)     :: elec(neelec)


      real   uz(nz,3)  ! ion mobility

      integer icrwmp
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
      
      integer iadvx, iadvy, iadvz
      parameter ( iadvx = 1, iadvy = 1, iadvz = 1 ) 

      integer icrwmn
      parameter ( icrwmn = 1 )
      
      integer   ibc,jbc
      
      integer n

      real, parameter :: eps = 1.0e-25

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!
!
!
!  CHECKS
!
!
!
!  check to see that work arrays wont overflow...
!
!      if ( ny .gt. nx ) then
!      write(*,*) 'nyl is greater than nxl;  ',ny,nx
!      write(*,*) 'stop in sam.a3.cx.f -- this may be outdated'
!      stop
!      end if
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
      icrwmp = 1
      IF ( ndebug .gt. 0 ) print*,'IONDRIFT',nstep
      if ( icrwmp .eq. 1 ) then
      dtfac = 1.0
      npass = 3
      mdnstp = mod(nstep,2)
      end if
      if ( icrwmp .eq. 2 ) then
      dtfac = 0.5
      npass = 6
      mdnstp = mod(nstep,2)
      end if
!
!
!
!  CONSTANTS
!
!
!

      ng = nor
      
      dxi   = (1.0)/dx
      dyi   = (1.0)/dy
      dzi   = (1.0)/dz
      dxt   = -dt/dx
      dyt   = -dt/dy
      dzt   = -dt/dz
      dxt4  = -dt/(4.0*dx)
      dyt4  = -dt/(4.0*dy)
      dzt4  = -dt/(4.0*dz)
      dtdx  = -dtfac*dt/dx
      dtdy  = -dtfac*dt/dy
      dtdz  = -dtfac*dt/dz
      dtdx4 =  (0.25)*dtfac*dt/dx
      dtdy4 =  (0.25)*dtfac*dt/dy
      dtdz4 =  (0.25)*dtfac*dt/dz
!
#ifdef MPI
    imn = 1
    imx = nxend-1
    jmn = 1
    jmx = nyend-1
    kmn = 1
    kmx = nzend-1 
#else
    imn = 1
    imx = nx-1
    jmn = 1
    jmx = ny-1
    kmn = 1
    kmx = nz-1
#endif

!
!
!  CREATE COEFFICIENTS FOR CROWELY FLUX SCHEME
!
!
!
!  create x-direction polynomial coefficients
!
      ibc = 0
      jbc = 0
      if ( nx .gt. 2 ) then
       IF ( bcx .ne. 2 ) THEN
!        call crwcff(np,nx,nx,istag,id1,px)
        call crwcff(np,nx,nx,ng,istag,id1,itile,ixbeg,ixend,nxbeg,nxend,px)
       ELSE
!      subroutine crwcff1(norder,nap,nm,ng,nstag,nd1,p)
        call crwcff1(np,nx,nx,nor,0,0,px)
        ibc = id1 ! + istag
       ENDIF
      end if
!
!  create y-direction polynomial coefficients
!
!      subroutine crwcff(norder,nap,nm,ng,nstag,id1,  &
!                       itile,ixbeg,ixend,nxbeg,nxend,p)

      if ( ny .gt. 2 ) then
       IF ( bcy .ne. 2 ) THEN
!        call crwcff(np,ny,ny,jstag,jd1,py)
        call crwcff(np,ny,ny,ng,jstag,jd1,jtile,jybeg,jyend,nybeg,nyend,py)
       ELSE
        call crwcff1(np,ny,ny,nor,0,0,py)
        jbc = jd1 ! + jstag
       ENDIF
      end if
!
!  create z-direction polynomial coefficients
!
!      call crwcff(np,nz,nz,kstag,kd1,pz)
      call crwcff(np,nz,nz,ng,kstag,kd1,ktile,kzbeg,kzend,nzbeg,nzend,pz)
!
!
!
!  SIXTH ORDER CROWLEY FLUX
!
!
!
      scionfx(1) = 0.0d0
      scionfx(2) = 0.0d0
      scionfx(3) = 0.0d0
      scionfx(4) = 0.0d0
      scionfx(5) = 0.0d0
      scionfx(6) = 0.0d0


  IF ( bcx .eq. 1 .or. bcx .ge. 3 ) THEN   ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug: no need to fill ghost zones
   DO n = 1,ng
    an( 1-n   ,1:ny,1:nz,ia) = an(1      ,1:ny,1:nz,ia)
    an( nx-1+n,1:ny,1:nz,ia) = an(nx-1,1:ny,1:nz,ia)
   ENDDO
#else
   IF ( ixbeg .eq. nxbeg ) THEN
   DO n = 1,ng
    an( 1-n      ,-ng+1:ny+ng,1:nz,ia) = an(1      ,-ng+1:ny+ng,1:nz,ia)
   ENDDO
   ENDIF
   
   IF ( ixend .eq. nxend ) THEN
   DO n = 1,ng
    an( nx-1+n,-ng+1:ny+ng,1:nz,ia) = an(nx-1,-ng+1:ny+ng,1:nz,ia)
   ENDDO
   ENDIF
#endif

   ENDIF

  IF ( bcy .eq. 1 .or. bcy .ge. 3 ) THEN    ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug:  no need to fill ghost zones
   DO n = 1,ng
    an(1:nx, 1-n  ,1:nz,ia) = an(1:nx,1   ,1:nz,ia)
    an(1:nx,ny-1+n,1:nz,ia) = an(1:nx,ny-1,1:nz,ia)
   ENDDO
#else
   IF ( jybeg .eq. nybeg ) THEN
   DO n = 1,ng
    an(-ng+1:nx+ng, 1-n     ,1:nz,ia) = an(-ng+1:nx+ng,1      ,1:nz,ia)
   ENDDO
   ENDIF

   IF ( jyend .eq. nyend ) THEN
   DO n = 1,ng
    an(-ng+1:nx+ng,ny-1+n,1:nz,ia) = an(-ng+1:nx+ng,ny-1,1:nz,ia) 
   ENDDO
   ENDIF

#endif
   ENDIF
      
      do 9000 ipass = 1,npass
!
!  x-direction crowley corrected flux computations
!
      if ( ipass .eq. (1+mdnstp) .or. ipass .eq. (6-mdnstp) ) then 
!
      if ( iadvx .ge. 1 ) then 
      if ( nx .gt. 5 ) then

!      fu(:,:,:) = 0.0

#ifdef MPI
      ixb = 0
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
      DO kz = 1, nz-kstag
      DO jy = 1, ny-jstag
      DO ix = 1-ibc,     nx-istag-id1+ibc
#endif
!
!c C$DOACROSS LOCAL(kz,jy,ix)
! !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)  
!      do 2001 kz = 1,     nz-kstag
!      do 2002 jy = 1,     ny-jstag
!      do 2003 ix = 1-ibc,     nx-istag-id1+ibc
      IF ( iestag .eq. 1 ) THEN
       
       fu(ix,jy,kz) =  &
        dtdx4 &
       *4.0*s*uz(kz,1)*(elec(iex)%flt3d(ix+id1,jy,kz)) &
       /((0.50)*dx*(gx(ix) + gx(ix+id1)))
       
!       fu(ix,jy,kz) =  &
!        -dt*s*uz(kz,1)*elec(iex)%flt3d(ix+1,jy,kz) &
!       /((0.50)*(gx(ix) + gx(ix+id1)))
      
      ELSE
       
       fu(ix,jy,kz) =  &
        dtdx4 &
       *2.0*s*uz(kz,1)*(elec(iex)%flt3d(ix+id1,jy,kz)+ &
                elec(iex)%flt3d(ix,jy,kz)) &
       /((0.50)*dx*(gx(ix) + gx(ix+id1)))
       ENDIF
       
       ENDDO
       ENDDO
       ENDDO
!
       t0(:,:,:) = 0.0

#ifdef MPI
      ixb = 0
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,  &
!$OMP   fus,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)  
      do  kz = 1,     nz-kstag
      do  jy = 1,     ny-jstag
      do  ix = 1-ibc, nx-istag-id1+ibc
#endif

      fus = fu(ix,jy,kz)

      IF ( Abs(fus) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF

      t1s =    &
        fus*(px(ix,1,1) +   &
        fus*(px(ix,1,2) +   &
        fus*(px(ix,1,3) +   &
        fus*(px(ix,1,4) +   &
        fus*(px(ix,1,5) +   &
        fus*(px(ix,1,6)))))))
      t2s =    &
        fus*(px(ix,2,1) +   &
        fus*(px(ix,2,2) +   &
        fus*(px(ix,2,3) +   &
        fus*(px(ix,2,4) +   &
        fus*(px(ix,2,5) +   &
        fus*(px(ix,2,6)))))))
      t3s =    &
        fus*(px(ix,3,1) +   &
        fus*(px(ix,3,2) +   &
        fus*(px(ix,3,3) +   &
        fus*(px(ix,3,4) +   &
        fus*(px(ix,3,5) +   &
        fus*(px(ix,3,6)))))))
      t4s =    &
        fus*(px(ix,4,1) +   &
        fus*(px(ix,4,2) +   &
        fus*(px(ix,4,3) +   &
        fus*(px(ix,4,4) +   &
        fus*(px(ix,4,5) +   &
        fus*(px(ix,4,6)))))))
      t5s =    &
        fus*(px(ix,5,1) +   &
        fus*(px(ix,5,2) +   &
        fus*(px(ix,5,3) +   &
        fus*(px(ix,5,4) +   &
        fus*(px(ix,5,5) +   &
        fus*(px(ix,5,6)))))))
      t6s =    &
        fus*(px(ix,6,1) +   &
        fus*(px(ix,6,2) +   &
        fus*(px(ix,6,3) +   &
        fus*(px(ix,6,4) +   &
        fus*(px(ix,6,5) +   &
        fus*(px(ix,6,6)))))))
      t0s =    &
        t1s*an(ix-id2,jy,kz,ia)+   &
        t2s*an(ix-id1,jy,kz,ia)+   &
        t3s*an(ix    ,jy,kz,ia)+   &
        t4s*an(ix+id1,jy,kz,ia)+   &
        t5s*an(ix+id2,jy,kz,ia)+   &
        t6s*an(ix+id3,jy,kz,ia)
      IF ( icrwmn .eq. 1 ) THEN ! .and. px(ix,2,1) .ne. 0.0 .and. px(ix,5,1) .ne. 0.0 ) THEN
      t0s = -t0s / (fus)
      qsr  = max(sign((1.0),fus), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix-id1*1,jy,kz,ia)   &
           + qsl*an(ix+id1*2,jy,kz,ia)
      qd   = qsr*an(ix+id1*1,jy,kz,ia)   &
           + qsl*an(ix,      jy,kz,ia)
      qc   = qsr*an(ix,      jy,kz,ia)   &
           + qsl*an(ix+id1*1,jy,kz,ia)
      qr   = qu + (qc - qu) / (abs(fus) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
!      IF ( DBG ) print*, 'ix,jy,kz,imn,imx,ia = ',ix,jy,kz,imn,imx,ia
      t0(ix,jy,kz) =   &
        -fus *    &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
       t0(ix,jy,kz) = t0s
      ENDIF
     ENDDO
     ENDDO
     ENDDO

!
! C$DOACROSS LOCAL(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,kz,jy,ix,ia), 
! C$&       share(t1,t2,t3,t4,t5,t6,an,ad,fu)
! C$OMP PARALLEL DO DEFAULT(SHARED),  
! c$omp+  PRIVATE(t0,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp,kz,jy,ix,ia)
! c$omp+ shared(t1,t2,t3,t4,t5,t6,an,ad,fu)
!       
!      IF ( bcx .eq. 2 ) THEN
!        t0(0,1:ny-1,1:nz-1) = t0(nx-1,1:ny-1,1:nz-1)
!      ENDIF
      
#ifdef MPI

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1-istag

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF

      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

       do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      do kz = 1,         nz-kstag
      do jy = 1,         ny-jstag
      do ix = 1+id1-ibc, nx-id1-istag+ibc
#endif
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)  &
        + (t0(ix,jy,kz)-t0(ix-id1,jy,kz)) &
        * (gx(ix)*dx)**2
!     >  + (ad(ix,jy,kz,ia))
!     >  * (fu(ix,jy,kz)-fu(ix-id1,jy,kz))
!     >  * (gt(ix,jy,kz,imapx)**2)
       enddo
       enddo
       enddo
!
!
!  x-direction boundaries
!
#ifdef MPI

      IF ( bcy .ne. 2 .and. ( ixbeg .eq. nxbeg .or. ixend .eq. nxend ) ) THEN

       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

       kzb = 1
       kze = ktile
       if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

       do kz = kzb,kze ; do jy = jyb,jye
#else
      IF ( bcx .ne. 2 ) THEN
! C$DOACROSS LOCAL(vnormw,vnorme,kz,jy), 
! C$&       share(elec,uz,an,ad)
! !$OMP  PARALLEL DO DEFAULT(SHARED), &
! !$omp PRIVATE(vnormw,vnorme,kz,jy,dv),SHARED(elec,s,uz,an,ad) 
!      do 2050 ia = na1,na1
      do kz = 1,nz-kstag
      do jy = 1,ny-jstag

#endif
      IF ( myproci == 1 ) THEN
        dv = dxx(1)*dyy(jy)*dzz(kz)
        scionfx(1) = scionfx(1) - t0(1,jy,kz)*dv ! /gt(1,jy,kz,imapz)
      ENDIF
      
      IF ( myproci == nproci ) THEN
        dv = dxx(nx-2)*dyy(jy)*dzz(kz)
        scionfx(2) = scionfx(2) + t0(nx-2,jy,kz)*dv !/gt(nx-2,jy,kz,imapz)
      ENDIF
      
      vnormw = 0.0
!c     >  0.5*s*uz(kz)*(elec(1,jy,kz,iex) +elec(1,jy,kz,iex))
!      an(1,jy,kz,ia) = an(1,jy,kz,ia) 
!     > + max(vnormw*dtdx, 0.0)*gt(1,jy,kz,imapx)
!     >   *(ad(1+id1,jy,kz,ia) - ad(1,jy,kz,ia))
!     > + (ainflo*min(vnormw*dtdx, 0.0) - dt*dtfac*binflo)
!     >   *(ad(1,jy,kz,ia) - ab(1,jy,kz,ia))
      vnorme = 0.0
!c     >  0.5*s*uz(kz)*(elec(nx,jy,kz,iex) +elec(nx,jy,kz,iex))
!      an(nx-istag,jy,kz,ia) = an(nx-istag,jy,kz,ia) 
!     > + min(vnorme*dtdx, 0.0)*gt(nx-istag,jy,kz,imapx)
!     >   *(ad(nx-istag,jy,kz,ia) - ad(nx-istag-id1,jy,kz,ia))
!     > - (ainflo*max(vnorme*dtdx, 0.0) + dt*dtfac*binflo)
!     >   *(ad(nx-istag,jy,kz,ia) - ab(nx-istag,jy,kz,ia))
      enddo
      enddo
 
      ENDIF

      end if
      end if
!
      end if
!
!  y-direction crowley corrected flux computation
!
      if ( ipass .eq. (2-mdnstp) .or. ipass .eq. (5+mdnstp) ) then
!
      if ( iadvy .ge. 1 ) then
      if ( ny .gt. 5 ) then

!      fv(:,:,:) = 0.0
!
#ifdef MPI
!      ixb = -ng+1
!      ixe = itile+ng
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 0
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1,     nz-kstag
      DO jy = 1-jbc, ny-jstag-jd1+jbc
      DO ix = 1,     nx-istag
#endif
      IF ( iestag .eq. 1 ) THEN
      
!      fv(ix,jy,kz) =  &
!       -dt*s*uz(kz,1)*(elec(iey)%flt3d(ix,jy+1,kz)) &
!       /((0.50)*(gy(jy)+gy(jy+jd1)))

      fv(ix,jy,kz) =  &
       dtdy4 &
       *4.0*s*uz(kz,1)* &
        (elec(iey)%flt3d(ix,jy+jd1,kz)) &
       /((0.50)*dy*(gy(jy)+gy(jy+jd1)))
      
      ELSE
      fv(ix,jy,kz) =  &
       dtdy4 &
       *2.0*s*uz(kz,1)* &
        (elec(iey)%flt3d(ix,jy+jd1,kz)+elec(iey)%flt3d(ix,jy,kz)) &
       /((0.50)*dy*(gy(jy)+gy(jy+jd1)))
      ENDIF
     ENDDO
     ENDDO
     ENDDO
!
  
      t0(:,:,:) = 0.0

#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      jyb = 0
      jye = jtile
      if (jybeg .le. nybeg) jyb = 1-jbc
      if (jyend .ge. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzbeg .le. nzbeg) kzb = 1
      if (kzend .ge. nzend) kze = kzend-kzbeg+1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix,     &
!$OMP   fvs,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)
      DO kz = 1,     nz-kstag
      DO jy = 1-jbc,     ny-jstag-jd1+jbc
      DO ix = 1,     nx-istag
#endif
      fvs = fv(ix,jy,kz)

      IF ( Abs(fvs) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF
      
      t1s =     &
        fvs*(py(jy,1,1) +   &
        fvs*(py(jy,1,2) +   &
        fvs*(py(jy,1,3) +   &
        fvs*(py(jy,1,4) +   &
        fvs*(py(jy,1,5) +   &
        fvs*(py(jy,1,6)))))))
      t2s =    &
        fvs*(py(jy,2,1) +   &
        fvs*(py(jy,2,2) +   &
        fvs*(py(jy,2,3) +   &
        fvs*(py(jy,2,4) +   &
        fvs*(py(jy,2,5) +   &
        fvs*(py(jy,2,6)))))))
      t3s =    &
        fvs*(py(jy,3,1) +   &
        fvs*(py(jy,3,2) +   &
        fvs*(py(jy,3,3) +   &
        fvs*(py(jy,3,4) +   &
        fvs*(py(jy,3,5) +   &
        fvs*(py(jy,3,6)))))))
      t4s =    &
        fvs*(py(jy,4,1) +   &
        fvs*(py(jy,4,2) +   &
        fvs*(py(jy,4,3) +   &
        fvs*(py(jy,4,4) +   &
        fvs*(py(jy,4,5) +   &
        fvs*(py(jy,4,6)))))))
      t5s =    &
        fvs*(py(jy,5,1) +   &
        fvs*(py(jy,5,2) +   &
        fvs*(py(jy,5,3) +   &
        fvs*(py(jy,5,4) +   &
        fvs*(py(jy,5,5) +   &
        fvs*(py(jy,5,6)))))))
      t6s =    &
        fvs*(py(jy,6,1) +   &
        fvs*(py(jy,6,2) +   &
        fvs*(py(jy,6,3) +   &
        fvs*(py(jy,6,4) +   &
        fvs*(py(jy,6,5) +   &
        fvs*(py(jy,6,6)))))))
      t0s =    &
        t1s*an(ix,jy-jd2,kz,ia) +   &
        t2s*an(ix,jy-jd1,kz,ia) +   &
        t3s*an(ix,jy    ,kz,ia) +   &
        t4s*an(ix,jy+jd1,kz,ia) +   &
        t5s*an(ix,jy+jd2,kz,ia) +   &
        t6s*an(ix,jy+jd3,kz,ia)
      IF ( icrwmn .eq. 1 ) THEN ! .and. py(jy,2,1) .ne. 0.0 .and. py(jy,5,1) .ne. 0.0 ) THEN
      t0s = -t0s / (fvs)
      qsr  = max(sign((1.0),fvs), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix,jy-jd1*1,kz,ia)   &
           + qsl*an(ix,jy+jd1*2,kz,ia)
      qd   = qsr*an(ix,jy+jd1*1,kz,ia)   &
           + qsl*an(ix,jy      ,kz,ia)
      qc   = qsr*an(ix,jy      ,kz,ia)   &
           + qsl*an(ix,jy+jd1*1,kz,ia)
      qr   = qu + (qc - qu) / (abs(fvs) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
      t0(ix,jy,kz) =   &
       -fvs *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
        t0(ix,jy,kz) = t0s
      ENDIF

     ENDDO
     ENDDO
     ENDDO


!      IF ( bcy .eq. 2 ) THEN
!        t0(1:nx-1, 0, 1:nz-1) = t0(1:nx-1, ny-1, 1:nz-1)
!      ENDIF
!
!
!
#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = imn
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn+jd1-jbc
      if (jyend .ge. jmx) jye = jmx-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1,         nz-kstag
      DO jy = 1+jd1-jbc, ny-jstag-jd1+jbc
      DO ix = 1,         nx-istag
#endif
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)    &
        + (t0(ix,jy,kz)-t0(ix,jy-jd1,kz))    &
        * (gy(jy)*dy)**2
!     >  + (ad(ix,jy,kz,ia))
!     >  * (fv(ix,jy,kz)-fv(ix,jy-jd1,kz))
!     >  * (gt(ix,jy,kz,imapy)**2)
     ENDDO
     ENDDO
     ENDDO
!
! 3020 continue
!
!  y-direction boundaries
!  scalars on south and north boundaries
!
#ifdef MPI
      IF ( bcy .ne. 2 .and. ( jybeg .eq. nybeg .or. jyend .eq. nyend ) ) THEN
#else
      IF ( bcy .ne. 2 ) THEN
#endif


#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(vnorms,vnormn,kz,ix,ia,dv), SHARED(elec,s,uz,an,ad)  
!      do 3051 ia = na1,na1
      DO kz = 1,     nz-kstag
      DO ix = 1,     nx-istag
#endif
#ifdef MPI
      IF ( jybeg .eq. nybeg ) THEN
#else
      IF ( jmn .eq. 1 ) THEN
#endif
       dv = dxx(ix)*dyy(1)*dzz(kz)
      scionfx(3) = scionfx(3) - t0(ix,1,kz)*dv !/gt(ix,1,kz,imapz)
      ENDIF

#ifdef MPI
      IF ( jyend .eq. nyend ) THEN
#else
      IF ( jmx .eq. ny-jstag ) THEN
#endif
       dv = dxx(ix)*dyy(ny-2)*dzz(kz)
      scionfx(4) = scionfx(4) + t0(ix,ny-2,kz)*dv !/gt(ix,ny-2,kz,imapz)
      ENDIF
      vnorms = 0.0
!c     >  0.5*s*uz(kz)*(elec(ix,1+jstag,kz,iey)+elec(ix,1,kz,iey))
!      an(ix,1,kz,ia) = an(ix,1,kz,ia) 
!     > + max(vnorms*dtdy, 0.0)*gt(ix,1,kz,imapy)
!     >   *(ad(ix,1+jd1,kz,ia) - ad(ix,1,kz,ia))
!     > + (ainflo*min(vnorms*dtdy, 0.0) - dt*dtfac*binflo)
!     >   *(ad(ix,1,kz,ia)- ab(ix,1,kz,ia))
      vnormn = 0.0
!c     >  0.5*s*uz(kz)*(elec(ix,ny-jstag,kz,iey)+elec(ix,ny,kz,iey))
!      an(ix,ny-jstag,kz,ia) = an(ix,ny-jstag,kz,ia) 
!     > + min(vnormn*dtdy, 0.0)*gt(ix,ny-jstag,kz,imapy)
!     >   *(ad(ix,ny-jstag,kz,ia) - ad(ix,ny-jstag-jd1,kz,ia))
!     > - (ainflo*max(vnormn*dtdy, 0.0) + dt*dtfac*binflo)
!     >   *(ad(ix,ny-jstag,kz,ia) - ab(ix,ny-jstag,kz,ia))
     ENDDO
     ENDDO

      ENDIF
      
      
      end if
      end if
!
      end if
!
!  z-direction crowley corrected flux computation
!
!     if ( ipass .eq. (3+mdnstp) .or. ipass .eq. (4-mdnstp) ) then
      
      
      IF ( ndebug .gt. 0 ) print*,'ipass,iadvz,nz = ',ipass,iadvz,nz
      if ( ipass .eq. 3 .or. ipass .eq. 4 ) then
!
      if ( iadvz .ge. 1 ) then
!
      if ( nz .gt. 5 ) then

!      fw(:,:,:) = 0.0

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      DO kz = 1,     nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jstag
      DO ix = 1,     nx-istag
#endif
!      IF ( iestag .eq. 1 ) THEN
!       fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz+kd1,3)*elec(iez)%flt3d(ix,jy,kz+1)
!      ELSE
       fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz+kd1,3)*ezw(ix,jy,kz+kd1) 
!      ENDIF
!      fws          = dtdz4*4.0*s*uz(kz+kd1,3)*ezw(ix,jy,kz+kd1)
!      IF ( ndebug .ge. 2 .and. ix .eq. nx/2 .and. jy .eq. ny/2 ) THEN
!       write(6,*) 'fw: ',kz,fws, fw(ix,jy,kz), dtdz4,s,uz(kz+kd1,3),ezw(ix,jy,kz+kd1)
!      ENDIF

      ENDDO
      ENDDO
      ENDDO
!
!      if ( itopo .ge. 1 ) then
!c
!      write(6,*) 'No ions with itopo .ge. 1 !!!'
!      write(0,*) 'No ions with itopo .ge. 1 !!!'
!      STOP
!      
!      call setfo
!     > (nx,ny,nz,dx,dy,dz,ht,gt,dc,elec,
!     :  fo,
!     :  t0,t1,t2,t3,t4,t5,t6,t7,t8,t9)
!c
!C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
!      do 29121 kz = 1,     nz-kstag
!      do 29122 jy = 1,     ny-jstag
!      do 29123 ix = 1,     nx-istag
!      t7(ix,jy,kz) = fo(ix,jy,kz)
!29123 continue
!29122 continue
!29121 continue
!c
!      call setfo
!     > (nx,ny,nz,dx,dy,dz,ht,gt,dc,elec,
!     :  fo,
!     :  t0,t1,t2,t3,t4,t5,t6,t7,t8,t9)
!c
!C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
!      do 29131 kz = 1,     nz-kstag
!      do 29132 jy = 1,     ny-jstag
!      do 29133 ix = 1,     nx-istag
!      fw(ix,jy,kz) =
!     > dtdz4*
!     > ((t7(ix,jy,kz+kd1)+t7(ix,jy,kz+kstag))
!     < +(fo(ix,jy,kz+kd1)+fo(ix,jy,kz+kstag)))
!29133 continue
!29132 continue
!29131 continue
!c
!      end if
!
!  compute polynomials
!
! C$DOACROSS LOCAL(kz,jy,ix)
! !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
!      do 4011 kz = 1,     nz-kstag-kd1
!      do 4012 jy = 1,     ny-jstag
!      do 4013 ix = 1,     nx-istag
! 4013 continue
! 4012 continue
! 4011 continue
!
!      do 4020 ia = na1,na1
!

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(ix,jy,kz,   &
!$OMP   fws,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp) 
      DO kz = 1,     nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jstag
      DO ix = 1,     nx-istag
#endif
      fws = fw(ix,jy,kz)

      IF ( Abs(fws) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF

      t1s =    &
        fws*(pz(kz,1,1) +   &
        fws*(pz(kz,1,2) +   &
        fws*(pz(kz,1,3) +   &
        fws*(pz(kz,1,4) +   &
        fws*(pz(kz,1,5) +   &
        fws*(pz(kz,1,6)))))))
      t2s =    &
        fws*(pz(kz,2,1) +   &
        fws*(pz(kz,2,2) +   &
        fws*(pz(kz,2,3) +   &
        fws*(pz(kz,2,4) +   &
        fws*(pz(kz,2,5) +   &
        fws*(pz(kz,2,6)))))))
      t3s =    &
        fws*(pz(kz,3,1) +   &
        fws*(pz(kz,3,2) +   &
        fws*(pz(kz,3,3) +   &
        fws*(pz(kz,3,4) +   &
        fws*(pz(kz,3,5) +   &
        fws*(pz(kz,3,6)))))))
      t4s =    &
        fws*(pz(kz,4,1) +   &
        fws*(pz(kz,4,2) +   &
        fws*(pz(kz,4,3) +   &
        fws*(pz(kz,4,4) +   &
        fws*(pz(kz,4,5) +   &
        fws*(pz(kz,4,6)))))))
      t5s =    &
        fws*(pz(kz,5,1) +   &
        fws*(pz(kz,5,2) +   &
        fws*(pz(kz,5,3) +   &
        fws*(pz(kz,5,4) +   &
        fws*(pz(kz,5,5) +   &
        fws*(pz(kz,5,6)))))))
      t6s =    &
        fws*(pz(kz,6,1) +   &
        fws*(pz(kz,6,2) +   &
        fws*(pz(kz,6,3) +   &
        fws*(pz(kz,6,4) +   &
        fws*(pz(kz,6,5) +   &
        fws*(pz(kz,6,6)))))))
      t0s =    &
        t1s*an(ix,jy,kz-kd2,ia)+   &
        t2s*an(ix,jy,kz-kd1,ia)+   &
        t3s*an(ix,jy,kz    ,ia)+   &
        t4s*an(ix,jy,kz+kd1,ia)+   &
        t5s*an(ix,jy,kz+kd2,ia)+   &
        t6s*an(ix,jy,kz+kd3,ia)


      IF ( icrwmn .eq. 1 .and. kz .gt. 1 ) THEN
      t0s = -t0s / (fws)
      qsr  = max(sign((1.0),fws), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix,jy,kz-kd1*1,ia)   &
           + qsl*an(ix,jy,kz+kd1*2,ia)
      qd   = qsr*an(ix,jy,kz+kd1*1,ia)   &
           + qsl*an(ix,jy,kz      ,ia)
      qc   = qsr*an(ix,jy,kz      ,ia)   &
           + qsl*an(ix,jy,kz+kd1*1,ia)
      qr   = qu + (qc - qu) / (abs(fws) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
      t0(ix,jy,kz) =   &
       -fws *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      ELSE
       t0(ix,jy,kz) = t0s
      ENDIF

!      IF ( DBG .and. ia .eq. 2 .and. fws .gt. 0.01 )  &
!            print*, 'ix,jy,kz,dt,ia,fws,t0s = ',ix,jy,kz,ia,fws,t0s

      ENDDO
      ENDDO
      ENDDO
!
!  monotonic face values
!
!
      IF ( ndebug .ge. 2 ) THEN
      write(6,*) 'nx,ny,nz,id1,jd1,kd1 = ',nx,ny,nz,id1,jd1,kd1
      write(6,*) 'ia,kz,t0,t0,ex,ey,ez,ezw,fw'
      ix = nx/2
      jy = ny/2
      DO kz =  nz-kstag-kd1,1,-1
       write(6,'(2i3,8(1x,1pe12.5))') ia,kz, &
         t0(nx/2,ny/2,kz),t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1), &
         elec(iex)%flt3d(nx/2,ny/2,kz),elec(iey)%flt3d(nx/2,ny/2,kz),elec(iez)%flt3d(nx/2,ny/2,kz),ezw(nx/2,ny/2,kz), &
         fw(nx/2,ny/2,kz)
!        fw(nx/2,ny/2,kz),fw(nx/2,ny/2,kz) - fw(nx/2,ny/2,kz-kd1), &
!        ad(nx/2,ny/2,kz)*(fw(nx/2,ny/2,kz) - fw(nx/2,ny/2,kz-kd1)), &
!        ezfair(kz),elec(iez)%flt3d(ix,jy,kz), &
!       ezfair(kz)*uz(kz,1)*cion(1)*dt/dz, &
!       elec(iez)%flt3d(ix,jy,kz)*uz(kz,1)*cion(1)*dt/dz
!     :  ( (t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1))
!     >     +(ad(ix,jy,kz,ia))
!     >     *(fw(ix,jy,kz) - fw(ix,jy,kz-kd1)) ) 
      ENDDO
      ENDIF
      
      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'Ion drift, s= ',s
      ENDIF
      
      
      
      
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      do jy = jyb,jye ; do ix = ixb,ixe
#else
      do jy = 1,     ny-id1
      do ix = 1,     nx-jd1
#endif
      IF ( s .eq. 1.0 ) THEN
!        IF ( elec(ix,jy,1,iez) .gt. emaxc .and.
!     :       elec(ix,jy,2,iez) .ge. emaxc .and.
!     :       elec(ix,jy,3,iez) .ge. emaxc  ) THEN
! 2/17/03 tm subst. ezw for elec
        IF ( ezw(ix,jy,1) .gt. emaxc .and. &
!             ezw(ix,jy,2) .gt. emaxc .and. &
!             ezw(ix,jy,3) .gt. emaxc .and. &
             icorona .eq. 1  ) THEN
!
! positive corona current
!
!         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Min(15.e3,elec(ix,jy,1,iez))*
!     :      (Min(15.e3,elec(ix,jy,1,iez)) - emaxc)**2 *dt/dz
! 2/17/03 tm subst. ezw for elec
         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Min(15.e3,ezw(ix,jy,1))* &
            (Min(15.e3,ezw(ix,jy,1)) - emaxc)**2 *dt/dz

         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'pos corona current at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
!
! check limit (use 1-D Gauss's Law to determine ion density that will
!   bring ez down to emaxc
!
!         ez1 = elec(ix,jy,1,iez)
         ez1 = ezw(ix,jy,1)
         delc = - (eperao/ec)*(emaxc*Sign(1.0,ez1) - ez1)/dz
         IF ( -t0(ix,jy,0) .gt. delc ) THEN
           t0(ix,jy,0) = -delc
         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'pos corona current REDUCED at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
         ENDIF
        
!        IF ( elec(ix,jy,1,iez) .lt. 0.0 ) THEN
!
!  First check if ion density needs to be relaxed toward fair-weather value
!
        ELSEIF ( (elec(iez)%flt3d(ix,jy,2) .lt. elec(iez)%flt3d(ix,jy,1) .or.  &
         Sign(1.0,elec(iez)%flt3d(ix,jy,4)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,1))) &
          .and. Abs(elec(iez)%flt3d(ix,jy,1)-ezfair(1)) .gt. 100.0  ) THEN
        
        t0(ix,jy,0) = Max( 0.0, t0(ix,jy,1) ) +   &
            alpha*dt*(an(ix,jy,1,ia) - cion(1))/(dz*gz(1)) ! /gt(ix,jy,1,imapz)
        
!        ENDIF
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .lt. 0.0 ) THEN
!            t0(ix,jy, 0) =  -s*uz(1,1)*ezw(ix,jy,1)*
            t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)* &
            an(ix,jy,1,ia)*dt/dz
          
         
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1) .and. &
             elec(iez)%flt3d(ix,jy,2) .le. ezfair(2)) THEN
!
! pass through fair weather (or downward positive) current.
!
!         IF ( t0(ix,jy,1) .lt. 0.0 ) THEN
           t0(ix,jy, 0)  =  Max( 0.0, t0(ix,jy,1) )
!         ENDIF
!
! allow at most the fair weather current out of ground for Ez .le. ezfair
!
!         t0(ix,jy, 0) =  -s*uz(1,1)*ezfair(1)*
!     :   Min(an(ix,jy,1,na1), cion(1))*dt/dz
!         t0(ix,jy, 0) =  s*uz(1,1)*
!     :   Max(ezfair(1), elec(ix,jy,1,iez))*
!     :          Sign(1.0,elec(ix,jy,1,iez))*ab(ix,jy,1,na1)*dt/dz
!        ELSEIF ( elec(ix,jy,1,iez) .gt. emaxc ) THEN
        ELSE
         t0(ix,jy, 0) = 0.0
        ENDIF
      ELSEIF ( s .eq. -1.0 ) THEN
!        IF ( elec(ix,jy,1,iez) .lt. -emaxc .and.
!     :       elec(ix,jy,2,iez) .le. -emaxc .and.
!     :       elec(ix,jy,3,iez) .le. -emaxc  ) THEN
! tm 2/17/03
        IF ( ezw(ix,jy,1) .lt. -emaxc  .and. &
!             ezw(ix,jy,2) .lt. -emaxc .and. &
!             ezw(ix,jy,3) .lt. -emaxc .and. &
             icorona .eq. 1  ) THEN
!
! negative corona
!
!         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Max(-15.e3,elec(ix,jy,1,iez))*
!     :      (Max(-15.e3,elec(ix,jy,1,iez)) + emaxc)**2 *dt/dz
! tm 2/17/03
         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Max(-15.e3, ezw(ix,jy,1))* &
            (Max(-15.e3,ezw(ix,jy,1)) + emaxc)**2 *dt/dz
          IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'neg corona current at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
!
! check limit (use 1-D Gauss's Law to determine ion density that will
!   bring ez down to emaxc
!
!         ez1 = elec(ix,jy,1,iez)
! tm 2/17/03
         ez1 = ezw(ix,jy,1)
         delc =  (eperao/ec)*(emaxc*Sign(1.0,ez1) - ez1)/dz
         IF ( -t0(ix,jy,0) .gt. delc ) THEN
           t0(ix,jy,0) = -delc
         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'neg corona current REDUCED at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
         ENDIF
!
!  First check if ion density needs to be relaxed toward fair-weather value
!
        ELSEIF ( (elec(iez)%flt3d(ix,jy,2) .gt. elec(iez)%flt3d(ix,jy,1) .or.  &
         Sign(1.0,elec(iez)%flt3d(ix,jy,4)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,1))) &
          .and. Abs(elec(iez)%flt3d(ix,jy,1)-ezfair(1)) .gt. 100.0  ) THEN
        
        t0(ix,jy,0) = Max( 0.0, t0(ix,jy,1) )+ &
            alpha*dt*(an(ix,jy,1,ia) - cion(1))/(dz*gz(1)) ! /gt(ix,jy,1,imapz)
        
!        ENDIF
!        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1)) THEN
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1) .and. &
                 elec(iez)%flt3d(ix,jy,2) .le. ezfair(2) ) THEN
!
! allow at most the fair weather current out of ground for Ez .le. ezfair
!
         t0(ix,jy, 0) =  -s*uz(1,1)*ezfair(1)* &
         cion(1)*dt/dz
!     :   Min(an(ix,jy,1,na1), cion(1))*dt/dz
         
!         t0(ix,jy, 0) =  s*uz(1,1)*ezfair(1)*
!     :          Sign(1.0,elec(iez)%flt3d(ix,jy,1))*an(ix,jy,1,na1)*dt/dz
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .lt. 0.0 .and. &
                 elec(iez)%flt3d(ix,jy,2) .lt. 0.0 ) THEN
!
!  allow a smaller current (into ground) for small negative Ez
!
         t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)*cion(1)*dt/dz
         
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .gt. 0.0 .and. &
                 elec(iez)%flt3d(ix,jy,2) .gt. 0.0 ) THEN
!
!  pass through downward negative ion current
!
!         t0(ix,jy, 0) =  t0(ix,jy,1)

           t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)* &
            an(ix,jy,1,ia)*dt/dz
                   

!         t0(ix,jy, 0) =  s*uz(1,1)*
!     :    Min(-ezfair(1), elec(iez)%flt3d(ix,jy,1))*
!     :          Sign(1.0,elec(iez)%flt3d(ix,jy,1))*ab(ix,jy,1,na1)*dt/dz
        ELSE
         t0(ix,jy, 0) = 0.0
        ENDIF
      ENDIF
      fw(ix,jy, 0) = 0.0
!      t0(ix,jy,nz-kstag) = 0.0
!      fw(ix,jy,nz-kstag) = 0.0

! not multiplying these by gt because it will get multiplied by
! dz later, canceling out the 1/dz factor in t0.
      IF ( myprock == 1 ) THEN
        dv = dxx(ix)*dyy(jy)*dzz(1)
        scionfx(5) = scionfx(5) - t0(ix,jy,0)*dv ! *gt(ix,jy,1,imapz)
      ENDIF
      
      IF ( myprock == nprock ) THEN
        dv = dxx(ix)*dyy(jy)*dzz(nz-2)
        scionfx(6) = scionfx(6) + t0(ix,jy,nz-2)*dv ! *gt(ix,jy,nz-2,imapz)
      ENDIF

      ENDDO
      ENDDO

#ifdef MPI

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1


      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1, nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jd1
      DO ix = 1,     nx-id1
#endif

! without zero divergence
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)                  &
      + (t0(ix,jy,kz) - t0(ix,jy,kz-kd1))*(gz(kz)*dz)

      ENDDO
      ENDDO
      ENDDO

      IF ( ndebug .ge. 3 ) THEN
      write(6,*) 'ia,kz,t0,t0,fw,fw'
      ix = nx/2
      jy = ny/2
      DO kz =  nz-kstag-kd1,1,-1
       write(6,'(2i3,8(1x,1pe12.5))') ia,kz, &
         t0(nx/2,ny/2,kz),t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1), &
         elec(iex)%flt3d(nx/2,ny/2,kz),elec(iey)%flt3d(nx/2,ny/2,kz),elec(iez)%flt3d(nx/2,ny/2,kz),ezw(nx/2,ny/2,kz)
!        fw(nx/2,ny/2,kz),fw(nx/2,ny/2,kz) - fw(nx/2,ny/2,kz-kd1), &
!        ad(nx/2,ny/2,kz)*(fw(nx/2,ny/2,kz) - fw(nx/2,ny/2,kz-kd1)), &
!        ezfair(kz),elec(iez)%flt3d(ix,jy,kz), &
!       ezfair(kz)*uz(kz,1)*cion(1)*dt/dz, &
!       elec(iez)%flt3d(ix,jy,kz)*uz(kz,1)*cion(1)*dt/dz
!     :  ( (t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1))
!     >     +(ad(ix,jy,kz,ia))
!     >     *(fw(ix,jy,kz) - fw(ix,jy,kz-kd1)) ) 
      ENDDO
      ENDIF

!
! 4020 continue
!
      IF (.false.) THEN
!      IF (.true.) THEN

      IF ( s .gt. 0.0 ) THEN
        vnormb = -5.385811127460374e6*dt/dz
!        vnormb = ezfair(2)*uz(2,3)*ab(1,1,1,na1)*dt/dz
      ELSE
        vnormb = 7.098583379406049e6*dt/dz
      ENDIF
!      do  ia = na1,na1
      do  jy = 1,     ny-kstag
      do  ix = 1,     nx-istag
!      vnormb = s*uz(1)*elec(iez)%flt3d(ix,jy,1)
!      vnormb = s*jc/ec
      
      an(ix,jy,1,ia) = an(ix,jy,1,ia) + vnormb
    
!c     > + max(vnormb*dtdz, 0.0)*gt(ix,jy,1,imapz)
!     > + vnormb*dtdz*gt(ix,jy,1,imapz)
!     >   *(ad(ix,jy,1+kd1,ia) - ad(ix,jy,1,ia))
      IF ( jy .eq. ny/2 .and. ix.eq.nx/2 .and. ndebug.ge.2) THEN
        write(6,*) 'vnormb at nx/2,ny/2 = ',vnormb,  &
!     : max(vnormb*dtdz, 0.0)*gt(ix,jy,1,imapz)
       + vnormb*dtdz*dz*gz(1)                &  ! *gt(ix,jy,1,imapz)               &
         *(ad(ix,jy,1+kd1) - ad(ix,jy,1))   
      ENDIF
!      vnormn = 0.0
!     >  0.5*s*uz(kz)*(elec(ix,ny-jstag,kz,iey)+elec(ix,ny,kz,iey))
!      an(ix,ny-jstag,kz,ia) = an(ix,ny-jstag,kz,ia) 
!     > + min(vnormn*dtdy, 0.0)*gt(ix,ny-jstag,kz,imapy)
!     >   *(ad(ix,ny-jstag,kz,ia) - ad(ix,ny-jstag-jd1,kz,ia))
!     > - (ainflo*max(vnormn*dtdy, 0.0) + dt*dtfac*binflo)
!     >   *(ad(ix,ny-jstag,kz,ia) - ab(ix,ny-jstag,kz,ia))
      ENDDO
      ENDDO
!      ENDDO
      
      ENDIF

      end if
      end if
!
      end if
!
!
!  end of crowley flux scheme for scalars
!
!
 9000 continue
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
! 1st order upwind version of drift
!
      subroutine iondrift1st &
       (nx,ny,nz,nor,na,ia,s,emaxc,icorona,nstep, &
        dt,dx,dy,dz,cion,scionfx, &
        ab,ad,an,                 & ! ad is t7 or t8 &
        elec,uz,ezfair,           &
!        ex,ey,ez,uz,ezfair,           &
        t0,ezw,                   & ! ezw is t9 &
        fu,fv,fw,                 & ! t1, t2, t3 &
        dxx,dyy,dzz,gx,gy,gz,iestag,bcx,bcy) 
     
       USE GRID_MODULE
       USE INDEX_MODULE
       USE COMMASMPI_MODULE
       USE PARAM_MODULE, only: icrwmn1st
       
!
!  Revisions:
!   2.17.2003  Now using ezw(ix,jy,1) to calculate corona current
!              because it is more responsive to feedback of the 
!              new corona ions.
!
!
!  include file :  contains implicit none
!
      implicit none
!      include 'sam.param.h'
!      include 'sam.index.f'
!
!  declare all variables

      real dxx(nx),dyy(ny),dzz(nz)         ! dx(i),dy(j),dz(k)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)

      integer ng1, nor, nstep, ng
      parameter(ng1 = 1)

!      real ab(nz,na)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)

!
!      real pb(nz)
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
!      real db(nz)
      
!      real ex(nx,ny,nz)
!      real ey(nx,ny,nz)
!      real ez(nx,ny,nz)

!

      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ad(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ezw(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

      real fu(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real fw(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      integer bcx, bcy

      real fus,fvs,fws

      double precision t0s,t1s,t2s,t3s,t4s,t5s,t6s

      integer imn,imx
      integer jmn,jmx
      integer kmn,kmx
!
!  INTEGERS
!
      integer    ia,na
      integer    ip
      integer    ipass
      integer    i,j,k
      integer    ix
      integer    jy
      integer    kz
      integer    mdnstp
      integer    n1
      integer    n2
      integer    np
      integer    npass
      integer    nx
      integer    ny
      integer    nz
      parameter (np=6)
      
      integer icorona
      
      integer ndebug
      parameter (ndebug=0)
      
      integer iestag
!
!  REALS
!
      double precision scionfx(6)
      
      real dv
      
      real delc,ez1
      
      real jc
      parameter ( jc = -2.0e-12 )
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)

      real eperao
      parameter (eperao  = 8.8592e-12 )
      
      real alpha ! rate for forcing ions at ground back to fair weather value
      parameter (alpha = 0.0577/5.0)

      real emaxc ! point discharge (corona) threshold

      double precision qu, qd, qc, qr, qdel, qxs, qxp, qsr, qsl
      
      real  s   ! sign of ion charge
      real ezfair(nz)
      real cion(nz)  !  base state ion density (pos or neg depending)
!
      real       dt
      real       dx
      real       dy
      real       dz
      real       dxi
      real       dyi
      real       dzi
      real       dxt
      real       dyt
      real       dzt
      real       dtdx
      real       dtdy
      real       dtdz
      real       dxt4
      real       dyt4
      real       dzt4
      real       dtdx4
      real       dtdy4
      real       dtdz4
      real       dtfac
      real       vnorme
      real       vnormw
      real       vnormn
      real       vnorms
      real vnormb,vnormt
      real      tmp
!
!      real       cx(nxl,nzl,np),  
!      real  px(-nor+1:nx+nor,np,np),py(-nor+1:ny+nor,np,np),pz(-nor+1:nz+nor,np,np)
!      real       cy(nyl,nzl,np),  py(nyl,np,np)
!      real       cz(nxl,nzl,np),  pz(nzl,np,np)
!
!

      TYPE(VARIABLE), INTENT(IN)     :: elec(neelec)


      real   uz(nz,3)  ! ion mobility

      integer icrwmp
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
      
      integer iadvx, iadvy, iadvz
      parameter ( iadvx = 1, iadvy = 1, iadvz = 1 ) 

      integer icrwmn
      parameter ( icrwmn = 1 )
      
      integer   ibc,jbc
      
      integer n

      real, parameter :: eps = 1.0e-25

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

!
!
!
!  CHECKS
!
!
!
!  check to see that work arrays wont overflow...
!
!      if ( ny .gt. nx ) then
!      write(*,*) 'nyl is greater than nxl;  ',ny,nx
!      write(*,*) 'stop in sam.a3.cx.f -- this may be outdated'
!      stop
!      end if
!
!
!
!  SET RUNTIME PARAMETER TO DETERMINE IMPLEMENTATION OF CROWLEY SCHEME
!
!    icrwmp = 1;   at step 1, do x-pass, y-pass, z-pass
!                  at step 2, do y-pass, x-pass, z-pass
! 
!  set type of crowley implementation ( sam.param.h )
!
      icrwmp = 1
      IF ( ndebug .gt. 0 ) print*,'IONDRIFT',nstep
      
      dtfac = 1.0
      npass = 3
      mdnstp = mod(nstep,2)
!
!
!
!  CONSTANTS
!
!
!

      ng = nor
      
      dxi   = (1.0)/dx
      dyi   = (1.0)/dy
      dzi   = (1.0)/dz
      dxt   = -dt/dx
      dyt   = -dt/dy
      dzt   = -dt/dz
      dxt4  = -dt/(4.0*dx)
      dyt4  = -dt/(4.0*dy)
      dzt4  = -dt/(4.0*dz)
      dtdx  = -dtfac*dt/dx
      dtdy  = -dtfac*dt/dy
      dtdz  = -dtfac*dt/dz
      dtdx4 =  (0.25)*dtfac*dt/dx
      dtdy4 =  (0.25)*dtfac*dt/dy
      dtdz4 =  (0.25)*dtfac*dt/dz
!
#ifdef MPI
    imn = 1
    imx = nxend-1
    jmn = 1
    jmx = nyend-1
    kmn = 1
    kmx = nzend-1 
#else
    imn = 1
    imx = nx-1
    jmn = 1
    jmx = ny-1
    kmn = 1
    kmx = nz-1
#endif

!
!
!  CREATE COEFFICIENTS FOR CROWELY FLUX SCHEME
!
!
!
!  create x-direction polynomial coefficients
!
      ibc = 0
      jbc = 0
      if ( nx .gt. 2 ) then
       IF ( bcx .ne. 2 ) THEN
       ELSE
        ibc = id1 ! + istag
       ENDIF
      end if


      if ( ny .gt. 2 ) then
       IF ( bcy .ne. 2 ) THEN
       ELSE
        jbc = jd1 ! + jstag
       ENDIF
      end if
!
!  First ORDER CROWLEY FLUX
!
!
!
      scionfx(1) = 0.0d0
      scionfx(2) = 0.0d0
      scionfx(3) = 0.0d0
      scionfx(4) = 0.0d0
      scionfx(5) = 0.0d0
      scionfx(6) = 0.0d0


  IF ( bcx .eq. 1 .or. bcx .ge. 3 ) THEN   ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug: no need to fill ghost zones
   DO n = 1,ng
    an( 1-n   ,1:ny,1:nz,ia) = an(1      ,1:ny,1:nz,ia)
    an( nx-1+n,1:ny,1:nz,ia) = an(nx-1,1:ny,1:nz,ia)
   ENDDO
#else
   IF ( ixbeg .eq. nxbeg ) THEN
   DO n = 1,ng
    an( 1-n      ,-ng+1:ny+ng,1:nz,ia) = an(1      ,-ng+1:ny+ng,1:nz,ia)
   ENDDO
   ENDIF
   
   IF ( ixend .eq. nxend ) THEN
   DO n = 1,ng
    an( nx-1+n,-ng+1:ny+ng,1:nz,ia) = an(nx-1,-ng+1:ny+ng,1:nz,ia)
   ENDDO
   ENDIF
#endif

   ENDIF

  IF ( bcy .eq. 1 .or. bcy .ge. 3 ) THEN    ! zero gradient, set ghost zones to boundary value

#ifndef MPI
!mpidebug:  no need to fill ghost zones
   DO n = 1,ng
    an(1:nx, 1-n  ,1:nz,ia) = an(1:nx,1   ,1:nz,ia)
    an(1:nx,ny-1+n,1:nz,ia) = an(1:nx,ny-1,1:nz,ia)
   ENDDO
#else
   IF ( jybeg .eq. nybeg ) THEN
   DO n = 1,ng
    an(-ng+1:nx+ng, 1-n     ,1:nz,ia) = an(-ng+1:nx+ng,1      ,1:nz,ia)
   ENDDO
   ENDIF

   IF ( jyend .eq. nyend ) THEN
   DO n = 1,ng
    an(-ng+1:nx+ng,ny-1+n,1:nz,ia) = an(-ng+1:nx+ng,ny-1,1:nz,ia) 
   ENDDO
   ENDIF

#endif
   ENDIF
      
      do 9000 ipass = 1,npass
!
!  x-direction crowley corrected flux computations
!
      if ( ipass .eq. (1+mdnstp) .or. ipass .eq. (6-mdnstp) ) then 
!
      if ( iadvx .ge. 1 ) then 
      if ( nx .gt. 5 ) then


#ifdef MPI
      ixb = 0
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1+ibc

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,j,i)  
      DO kz = 1, nz-kstag
      DO jy = 1, ny-jstag
      DO ix = 1-ibc,     nx-istag-id1+ibc
#endif

      IF ( iestag .eq. 1 ) THEN
       
       fu(ix,jy,kz) =  &
        dtdx4 &
       *4.0*s*uz(kz,1)*(elec(iex)%flt3d(ix+id1,jy,kz)) &
       /((0.50)*dx*(gx(ix) + gx(ix+id1)))
       
      
      ELSE
       
       fu(ix,jy,kz) =  &
        dtdx4 &
       *2.0*s*uz(kz,1)*(elec(iex)%flt3d(ix+id1,jy,kz)+ &
                elec(iex)%flt3d(ix,jy,kz)) &
       /((0.50)*dx*(gx(ix) + gx(ix+id1)))
       ENDIF
       
       ENDDO
       ENDDO
       ENDDO
!
       t0(:,:,:) = 0.0

#ifdef MPI
      ixb = 0
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,  &
!$OMP   fus,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)  
      do  kz = 1,     nz-kstag
      do  jy = 1,     ny-jstag
      do  ix = 1-ibc, nx-istag-id1+ibc
#endif

      fus = fu(ix,jy,kz)

      IF ( Abs(fus) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF
      
       t0s = -fus*an(ix,jy,kz,ia)

      IF ( icrwmn1st >= 2 ) THEN ! .and. px(ix,2,1) .ne. 0.0 .and. px(ix,5,1) .ne. 0.0 ) THEN
      t0s = -t0s / (fus)
      qsr  = max(sign((1.0),fus), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix-id1*1,jy,kz,ia)   &
           + qsl*an(ix+id1*2,jy,kz,ia)
      qd   = qsr*an(ix+id1*1,jy,kz,ia)   &
           + qsl*an(ix,      jy,kz,ia)
      qc   = qsr*an(ix,      jy,kz,ia)   &
           + qsl*an(ix+id1*1,jy,kz,ia)
      qr   = qu + (qc - qu) / (abs(fus) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
!      IF ( DBG ) print*, 'ix,jy,kz,imn,imx,ia = ',ix,jy,kz,imn,imx,ia
      t0(ix,jy,kz) =   &
        -fus *    &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
       t0(ix,jy,kz) = t0s
      ENDIF

     ENDDO
     ENDDO
     ENDDO

!
      
#ifdef MPI

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1+id1-ibc
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1-istag

  IF ( ipass .eq. 1 ) THEN
  jyb = -ng+1
  jye = jtile+ng
  ELSE
  jyb = 1
  jye = jtile
  ENDIF

      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

       do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      do kz = 1,         nz-kstag
      do jy = 1,         ny-jstag
      do ix = 1+id1-ibc, nx-id1-istag+ibc
#endif
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)  &
        + (t0(ix,jy,kz)-t0(ix-id1,jy,kz)) &
        * (gx(ix)*dx)**2

       enddo
       enddo
       enddo
!
!
!  x-direction boundaries
!
#ifdef MPI

      IF ( bcy .ne. 2 .and. ( ixbeg .eq. nxbeg .or. ixend .eq. nxend ) ) THEN

       jyb = 1
       jye = jtile
       if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

       kzb = 1
       kze = ktile
       if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

       do kz = kzb,kze ; do jy = jyb,jye
#else
      IF ( bcx .ne. 2 ) THEN
! C$DOACROSS LOCAL(vnormw,vnorme,kz,jy), 
! C$&       share(elec,uz,an,ad)
! !$OMP  PARALLEL DO DEFAULT(SHARED), &
! !$omp PRIVATE(vnormw,vnorme,kz,jy,dv),SHARED(elec,s,uz,an,ad) 
!      do 2050 ia = na1,na1
      do kz = 1,nz-kstag
      do jy = 1,ny-jstag

#endif
      IF ( myproci == 1 ) THEN
        dv = dxx(1)*dyy(jy)*dzz(kz)
        scionfx(1) = scionfx(1) - t0(1,jy,kz)*dv ! /gt(1,jy,kz,imapz)
      ENDIF
      
      IF ( myproci == nproci ) THEN
        dv = dxx(nx-2)*dyy(jy)*dzz(kz)
        scionfx(2) = scionfx(2) + t0(nx-2,jy,kz)*dv !/gt(nx-2,jy,kz,imapz)
      ENDIF
      
      vnormw = 0.0
      vnorme = 0.0

      enddo
      enddo
 
      ENDIF

      end if
      end if
!
      end if
!
!  y-direction crowley corrected flux computation
!
      if ( ipass .eq. (2-mdnstp) .or. ipass .eq. (5+mdnstp) ) then
!
      if ( iadvy .ge. 1 ) then
      if ( ny .gt. 5 ) then

!      fv(:,:,:) = 0.0
!
#ifdef MPI

      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 0
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1-jbc
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1,     nz-kstag
      DO jy = 1-jbc, ny-jstag-jd1+jbc
      DO ix = 1,     nx-istag
#endif
      IF ( iestag .eq. 1 ) THEN
      
!     fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz+kd1,3)*ezw(ix,jy,kz+kd1) 

      fv(ix,jy,kz) =  &
       dtdy4 &
       *4.0*s*uz(kz,1)* &
        (elec(iey)%flt3d(ix,jy+jd1,kz)) &
       /((0.50)*dy*(gy(jy)+gy(jy+jd1)))
      
      ELSE
      fv(ix,jy,kz) =  &
       dtdy4 &
       *2.0*s*uz(kz,1)* &
        (elec(iey)%flt3d(ix,jy+jd1,kz)+elec(iey)%flt3d(ix,jy,kz)) &
       /((0.50)*dy*(gy(jy)+gy(jy+jd1)))
      ENDIF
     ENDDO
     ENDDO
     ENDDO
!
  
      t0(:,:,:) = 0.0

#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1

      jyb = 0
      jye = jtile
      if (jybeg .le. nybeg) jyb = 1-jbc
      if (jyend .ge. nyend) jye = jyend-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzbeg .le. nzbeg) kzb = 1
      if (kzend .ge. nzend) kze = kzend-kzbeg+1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix,     &
!$OMP   fvs,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp)
      DO kz = 1,     nz-kstag
      DO jy = 1-jbc,     ny-jstag-jd1+jbc
      DO ix = 1,     nx-istag
#endif
      fvs = fv(ix,jy,kz)

      IF ( Abs(fvs) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF
      
        t0s = -fvs*an(ix,jy,kz,ia)
        
      IF ( icrwmn1st >= 2 ) THEN ! .and. py(jy,2,1) .ne. 0.0 .and. py(jy,5,1) .ne. 0.0 ) THEN
      t0s = -t0s / (fvs)
      qsr  = max(sign((1.0),fvs), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix,jy-jd1*1,kz,ia)   &
           + qsl*an(ix,jy+jd1*2,kz,ia)
      qd   = qsr*an(ix,jy+jd1*1,kz,ia)   &
           + qsl*an(ix,jy      ,kz,ia)
      qc   = qsr*an(ix,jy      ,kz,ia)   &
           + qsl*an(ix,jy+jd1*1,kz,ia)
      qr   = qu + (qc - qu) / (abs(fvs) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
      t0(ix,jy,kz) =   &
       -fvs *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      
      ELSE
        t0(ix,jy,kz) = t0s
      ENDIF

     ENDDO
     ENDDO
     ENDDO



!
#ifdef MPI
      IF ( ipass .eq. 1 ) THEN
      ixb = -ng+1
      ixe = itile+ng
      ELSE
      ixb = 1
      ixe = itile
      ENDIF
      if (ixbeg .le. imn) ixb = imn
      if (ixend .ge. imx) ixe = imx-ixbeg+1

      jyb = 1
      jye = jtile
      if (jybeg .le. jmn) jyb = jmn+jd1-jbc
      if (jyend .ge. jmx) jye = jmx-jybeg+1-jd1+jbc

      kzb = 1
      kze = ktile
      if (kzbeg .le. kmn) kzb = kmn
      if (kzend .ge. kmx) kze = kmx-kzbeg+1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1,         nz-kstag
      DO jy = 1+jd1-jbc, ny-jstag-jd1+jbc
      DO ix = 1,         nx-istag
#endif
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)    &
        + (t0(ix,jy,kz)-t0(ix,jy-jd1,kz))    &
        * (gy(jy)*dy)**2

     ENDDO
     ENDDO
     ENDDO
!
! 3020 continue
!
!  y-direction boundaries
!  scalars on south and north boundaries
!
#ifdef MPI
      IF ( bcy .ne. 2 .and. ( jybeg .eq. nybeg .or. jyend .eq. nyend ) ) THEN
#else
      IF ( bcy .ne. 2 ) THEN
#endif


#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg-istag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      do kz = kzb,kze
      do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(vnorms,vnormn,kz,ix,ia,dv), SHARED(elec,s,uz,an,ad)  
!      do 3051 ia = na1,na1
      DO kz = 1,     nz-kstag
      DO ix = 1,     nx-istag
#endif
#ifdef MPI
      IF ( jybeg .eq. nybeg ) THEN
#else
      IF ( jmn .eq. 1 ) THEN
#endif
       dv = dxx(ix)*dyy(1)*dzz(kz)
      scionfx(3) = scionfx(3) - t0(ix,1,kz)*dv !/gt(ix,1,kz,imapz)
      ENDIF

#ifdef MPI
      IF ( jyend .eq. nyend ) THEN
#else
      IF ( jmx .eq. ny-jstag ) THEN
#endif
       dv = dxx(ix)*dyy(ny-2)*dzz(kz)
      scionfx(4) = scionfx(4) + t0(ix,ny-2,kz)*dv !/gt(ix,ny-2,kz,imapz)
      ENDIF
      vnorms = 0.0

      vnormn = 0.0

     ENDDO
     ENDDO

      ENDIF
      
      
      end if
      end if
!
      end if
!
!  z-direction crowley corrected flux computation
!
!     if ( ipass .eq. (3+mdnstp) .or. ipass .eq. (4-mdnstp) ) then
      
      
      IF ( ndebug .gt. 0 ) print*,'ipass,iadvz,nz = ',ipass,iadvz,nz
      if ( ipass .eq. 3 .or. ipass .eq. 4 ) then
!
      if ( iadvz .ge. 1 ) then
!
      if ( nz .gt. 5 ) then

!      fw(:,:,:) = 0.0

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kd1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
! C$DOACROSS LOCAL(kz,jy,ix)
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      DO kz = 1,     nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jstag
      DO ix = 1,     nx-istag
#endif

       fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz+kd1,3)*ezw(ix,jy,kz+kd1) 
!       fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz,1)*ezw(ix,jy,kz+kd1)  ! test 1
!       fw(ix,jy,kz) = dtdz4*4.0*s*uz(kz+1,1)*ezw(ix,jy,kz+kd1)  ! test 2


      ENDDO
      ENDDO
      ENDDO


#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP  PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(ix,jy,kz,   &
!$OMP   fws,t0s,t1s,t2s,t3s,t4s,t5s,t6s,qsr,qsl,qu,qd,qc,qr,qdel,qxs,qxp) 
      DO kz = 1,     nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jstag
      DO ix = 1,     nx-istag
#endif
      fws = fw(ix,jy,kz)

      IF ( Abs(fws) .lt. eps ) THEN
       t0(ix,jy,kz) = 0.0
       CYCLE
      ENDIF



       t0s = -fws*an(ix,jy,kz,ia)
      
      IF ( icrwmn1st >= 1 .and. kz .gt. 1 ) THEN
      t0s = -t0s / (fws)
      qsr  = max(sign((1.0),fws), (0.0))
      qsl  = (1.0) - qsr
      qu   = qsr*an(ix,jy,kz-kd1*1,ia)   &
           + qsl*an(ix,jy,kz+kd1*2,ia)
      qd   = qsr*an(ix,jy,kz+kd1*1,ia)   &
           + qsl*an(ix,jy,kz      ,ia)
      qc   = qsr*an(ix,jy,kz      ,ia)   &
           + qsl*an(ix,jy,kz+kd1*1,ia)
      qr   = qu + (qc - qu) / (abs(fws) + 1.0e-20)
      qdel = qd - qu
      qxs  = max(sign((1.0d0), qdel), (0.0d0))
      qxp  = max(sign((1.0d0),   &
                 abs(qdel)-abs(qu-(2.0d0)*qc+qd)), (0.0d0))
      t0(ix,jy,kz) =   &
       -fws *   &
       (   &
         ((1.0)-qxp) * qc   &
       +        qxp   &
       * (      qxs  * min(max(t0s, qc), min(qr, qd))   &
       + ((1.0)-qxs) * max(min(t0s, qc), max(qr, qd)) ) )
      ELSE
       t0(ix,jy,kz) = t0s
      ENDIF



      ENDDO
      ENDDO
      ENDDO
!
!  monotonic face values
!
!
      IF ( ndebug .ge. 2 ) THEN
      write(6,*) 'nx,ny,nz,id1,jd1,kd1 = ',nx,ny,nz,id1,jd1,kd1
      write(6,*) 'ia,kz,t0,t0,ex,ey,ez,ezw,fw'
      ix = nx/2
      jy = ny/2
      DO kz =  nz-kstag-kd1,1,-1
       write(6,'(2i3,8(1x,1pe12.5))') ia,kz, &
         t0(nx/2,ny/2,kz),t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1), &
         elec(iex)%flt3d(nx/2,ny/2,kz),elec(iey)%flt3d(nx/2,ny/2,kz),elec(iez)%flt3d(nx/2,ny/2,kz),ezw(nx/2,ny/2,kz), &
         fw(nx/2,ny/2,kz)

      ENDDO
      ENDIF
      
      IF ( ndebug .ge. 1 ) THEN
       write(6,*) 'Ion drift, s= ',s
      ENDIF
      
      
      
      
#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      do jy = jyb,jye ; do ix = ixb,ixe
#else
      do jy = 1,     ny-id1
      do ix = 1,     nx-jd1
#endif
      IF ( s .eq. 1.0 ) THEN
!        IF ( elec(ix,jy,1,iez) .gt. emaxc .and.
!     :       elec(ix,jy,2,iez) .ge. emaxc .and.
!     :       elec(ix,jy,3,iez) .ge. emaxc  ) THEN
! 2/17/03 tm subst. ezw for elec
        IF ( ezw(ix,jy,1) .gt. emaxc .and. &
!             ezw(ix,jy,2) .gt. emaxc .and. &
!             ezw(ix,jy,3) .gt. emaxc .and. &
             icorona .eq. 1  ) THEN
!
! positive corona current
!
!         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Min(15.e3,elec(ix,jy,1,iez))*
!     :      (Min(15.e3,elec(ix,jy,1,iez)) - emaxc)**2 *dt/dz
! 2/17/03 tm subst. ezw for elec
         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Min(15.e3,ezw(ix,jy,1))* &
            (Min(15.e3,ezw(ix,jy,1)) - emaxc)**2 *dt/dz

         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'pos corona current at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
!
! check limit (use 1-D Gauss's Law to determine ion density that will
!   bring ez down to emaxc
!
!         ez1 = elec(ix,jy,1,iez)
         ez1 = ezw(ix,jy,1)
         delc = - (eperao/ec)*(emaxc*Sign(1.0,ez1) - ez1)/dz
         IF ( -t0(ix,jy,0) .gt. delc ) THEN
           t0(ix,jy,0) = -delc
         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'pos corona current REDUCED at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
         ENDIF
        
!        IF ( elec(ix,jy,1,iez) .lt. 0.0 ) THEN
!
!  First check if ion density needs to be relaxed toward fair-weather value
!
        ELSEIF ( (elec(iez)%flt3d(ix,jy,2) .lt. elec(iez)%flt3d(ix,jy,1) .or.  &
         Sign(1.0,elec(iez)%flt3d(ix,jy,4)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,1))) &
          .and. Abs(elec(iez)%flt3d(ix,jy,1)-ezfair(1)) .gt. 100.0  ) THEN
        
        t0(ix,jy,0) = Max( 0.0, t0(ix,jy,1) ) +   &
            alpha*dt*(an(ix,jy,1,ia) - cion(1))/(dz*gz(1)) ! /gt(ix,jy,1,imapz)
        
!        ENDIF
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .lt. 0.0 ) THEN
!            t0(ix,jy, 0) =  -s*uz(1,1)*ezw(ix,jy,1)*
            t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)* &
            an(ix,jy,1,ia)*dt/dz
          
         
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1) .and. &
             elec(iez)%flt3d(ix,jy,2) .le. ezfair(2)) THEN
!
! pass through fair weather (or downward positive) current.
!
!         IF ( t0(ix,jy,1) .lt. 0.0 ) THEN
           t0(ix,jy, 0)  =  Max( 0.0, t0(ix,jy,1) )
!         ENDIF
!
! allow at most the fair weather current out of ground for Ez .le. ezfair
!
!         t0(ix,jy, 0) =  -s*uz(1,1)*ezfair(1)*
!     :   Min(an(ix,jy,1,na1), cion(1))*dt/dz
!         t0(ix,jy, 0) =  s*uz(1,1)*
!     :   Max(ezfair(1), elec(ix,jy,1,iez))*
!     :          Sign(1.0,elec(ix,jy,1,iez))*ab(ix,jy,1,na1)*dt/dz
!        ELSEIF ( elec(ix,jy,1,iez) .gt. emaxc ) THEN
        ELSE
         t0(ix,jy, 0) = 0.0
        ENDIF
      ELSEIF ( s .eq. -1.0 ) THEN
!        IF ( elec(ix,jy,1,iez) .lt. -emaxc .and.
!     :       elec(ix,jy,2,iez) .le. -emaxc .and.
!     :       elec(ix,jy,3,iez) .le. -emaxc  ) THEN
! tm 2/17/03
        IF ( ezw(ix,jy,1) .lt. -emaxc  .and. &
!             ezw(ix,jy,2) .lt. -emaxc .and. &
!             ezw(ix,jy,3) .lt. -emaxc .and. &
             icorona .eq. 1  ) THEN
!
! negative corona
!
!         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Max(-15.e3,elec(ix,jy,1,iez))*
!     :      (Max(-15.e3,elec(ix,jy,1,iez)) + emaxc)**2 *dt/dz
! tm 2/17/03
         t0(ix,jy,0) = (-s)*(2.0e-20/ec)*Max(-15.e3, ezw(ix,jy,1))* &
            (Max(-15.e3,ezw(ix,jy,1)) + emaxc)**2 *dt/dz
          IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'neg corona current at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
!
! check limit (use 1-D Gauss's Law to determine ion density that will
!   bring ez down to emaxc
!
!         ez1 = elec(ix,jy,1,iez)
! tm 2/17/03
         ez1 = ezw(ix,jy,1)
         delc =  (eperao/ec)*(emaxc*Sign(1.0,ez1) - ez1)/dz
         IF ( -t0(ix,jy,0) .gt. delc ) THEN
           t0(ix,jy,0) = -delc
         IF ( ndebug .ge. 1 ) THEN
           write(6,*) 'neg corona current REDUCED at ix,jy',ix,jy, &
           elec(iez)%flt3d(ix,jy,1),elec(iez)%flt3d(ix,jy,2),elec(iez)%flt3d(ix,jy,3), &
            t0(ix,jy,0)
          ENDIF
         ENDIF
!
!  First check if ion density needs to be relaxed toward fair-weather value
!
        ELSEIF ( (elec(iez)%flt3d(ix,jy,2) .gt. elec(iez)%flt3d(ix,jy,1) .or.  &
         Sign(1.0,elec(iez)%flt3d(ix,jy,4)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,1))) &
          .and. Abs(elec(iez)%flt3d(ix,jy,1)-ezfair(1)) .gt. 100.0  ) THEN
        
        t0(ix,jy,0) = Max( 0.0, t0(ix,jy,1) )+ &
            alpha*dt*(an(ix,jy,1,ia) - cion(1))/(dz*gz(1)) ! /gt(ix,jy,1,imapz)
        
!        ENDIF
!        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1)) THEN
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .le. ezfair(1) .and. &
                 elec(iez)%flt3d(ix,jy,2) .le. ezfair(2) ) THEN
!
! allow at most the fair weather current out of ground for Ez .le. ezfair
!
         t0(ix,jy, 0) =  -s*uz(1,1)*ezfair(1)* &
         cion(1)*dt/dz
!     :   Min(an(ix,jy,1,na1), cion(1))*dt/dz
         
!         t0(ix,jy, 0) =  s*uz(1,1)*ezfair(1)*
!     :          Sign(1.0,elec(iez)%flt3d(ix,jy,1))*an(ix,jy,1,na1)*dt/dz
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .lt. 0.0 .and. &
                 elec(iez)%flt3d(ix,jy,2) .lt. 0.0 ) THEN
!
!  allow a smaller current (into ground) for small negative Ez
!
         t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)*cion(1)*dt/dz
         
        ELSEIF ( elec(iez)%flt3d(ix,jy,1) .gt. 0.0 .and. &
                 elec(iez)%flt3d(ix,jy,2) .gt. 0.0 ) THEN
!
!  pass through downward negative ion current
!
!         t0(ix,jy, 0) =  t0(ix,jy,1)

           t0(ix,jy, 0) =  -s*uz(1,1)*elec(iez)%flt3d(ix,jy,1)* &
            an(ix,jy,1,ia)*dt/dz
                   

!         t0(ix,jy, 0) =  s*uz(1,1)*
!     :    Min(-ezfair(1), elec(iez)%flt3d(ix,jy,1))*
!     :          Sign(1.0,elec(iez)%flt3d(ix,jy,1))*ab(ix,jy,1,na1)*dt/dz
        ELSE
         t0(ix,jy, 0) = 0.0
        ENDIF
      ENDIF
      fw(ix,jy, 0) = 0.0
!      t0(ix,jy,nz-kstag) = 0.0
!      fw(ix,jy,nz-kstag) = 0.0

! not multiplying these by gt because it will get multiplied by
! dz later, canceling out the 1/dz factor in t0.
      IF ( myprock == 1 ) THEN
        dv = dxx(ix)*dyy(jy)*dzz(1)
        scionfx(5) = scionfx(5) - t0(ix,jy,0)*dv ! *gt(ix,jy,1,imapz)
      ENDIF
      
      IF ( myprock == nprock ) THEN
        dv = dxx(ix)*dyy(jy)*dzz(nz-2)
        scionfx(6) = scionfx(6) + t0(ix,jy,nz-2)*dv ! *gt(ix,jy,nz-2,imapz)
      ENDIF

      ENDDO
      ENDDO

#ifdef MPI

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-id1

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jd1

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag-kd1


      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!
!$OMP PARALLEL DO DEFAULT(SHARED), &
!$omp PRIVATE(kz,jy,ix)
      DO kz = 1, nz-kstag-kd1 ! +iestag
      DO jy = 1,     ny-jd1
      DO ix = 1,     nx-id1
#endif

! without zero divergence
      an(ix,jy,kz,ia) = an(ix,jy,kz,ia)                  &
      + (t0(ix,jy,kz) - t0(ix,jy,kz-kd1))*(gz(kz)*dz)

      ENDDO
      ENDDO
      ENDDO

      IF ( ndebug .ge. 3 ) THEN
      write(6,*) 'ia,kz,t0,t0,fw,fw'
      ix = nx/2
      jy = ny/2
      DO kz =  nz-kstag-kd1,1,-1
       write(6,'(2i3,8(1x,1pe12.5))') ia,kz, &
         t0(nx/2,ny/2,kz),t0(nx/2,ny/2,kz) - t0(nx/2,ny/2,kz-kd1), &
         elec(iex)%flt3d(nx/2,ny/2,kz),elec(iey)%flt3d(nx/2,ny/2,kz),elec(iez)%flt3d(nx/2,ny/2,kz),ezw(nx/2,ny/2,kz)
      ENDDO
      ENDIF

!
! 4020 continue
!

      end if
      end if
!
      end if
!
!
!  end of crowley flux scheme for scalars
!
!
 9000 continue
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
