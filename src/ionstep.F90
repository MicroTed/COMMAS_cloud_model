#ifndef RKIND
#define RKIND 4
#endif

!VD$F SKIP
!c
!
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
!  SUBROUTINE IONSTEP  (E. Mansell)
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!
! Purpose: subroutine to handle ion processes of ion drift,
!  attachment (diffusion and conduction) cosmic ray generation, 
!  recombination
!
!  7.14.2003: version optimized for OMP and memory
!
!  2.17.2003: Removed the fair-weather field 'correction' and changed
!             the way ezw is calculated (used to be an average of Ez from
!             the neighboring scalar points, but now calculate directly
!             from the potential at the neighboring scalar points).
!             The result was the elimination of ion oscillations near
!             the surface and above cloud top (and probably elsewhere).
!
!             Also added some extra net charge (and net ion charge) counters
!             in double precision.
!
!  12.27.2002 version 01b: eliminate use of sam.index.ion.h by passing
!   in the range of indices for hydrometeors and space charges
!
!  10.2.2002: version 01a: completely remove charge from hydrometeors
!    below qxmin if charge density magnitude is less than 1.0e-12 (still
!    goes into ions)
!
!    Also set ezw = ez at kz=1 (instead of fair-weather field)
!
!
      subroutine ionstep       &
     & (nx,ny,nz,na,nba,nor,nstep,       &
     &  dt1,dx,dy,dz,gxt,gyt,gzt,dxx,dyy,dzz,       &
     &  an,dn,pn,       &
     &  db,pb,       &
     &  elec,uz,cion,       &
     &  cwccn,ipconc,cimn,cimx,  &
     &  t0,t1,t2,t3,t4,t5,t6,t7,t8,ezw,       &
     &  tt0,tt7,       &
     &   nnxs,nnys,nnzs,nztop,       &
     &   iixps,jjyqs,kkzrs,iiexs,jjeys,kkezs,       &
     &   nbw,nbe,nbs,nbn,       &
     &   llworks,istretch,       &
     &   ibg,iunit0,       &
     &   id1,jd1,kd1,istag,jstag,kstag,       &
     &   microp,iestag,bcx,bcy,       &
     &   epot)
     
       USE ELEC_MODULE
       USE INDEX_MODULE
       USE GRID_MODULE
       USE CPUTIME_MODULE
       USE COMMASMPI_MODULE
       USE MICRO_MODULE, only: takcxmin, chaffconc


!
! TO DO:
!
!  deal with swapping ions in an,ac,ad after each subtime step
!  if more than one step.  Dont do after last step.  Maybe use ab as 
!  intermediary
!
!

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"

#ifdef BOXMG
#include        "BMG_constants.h"
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'
#endif


#endif

      TYPE(VARIABLE)     :: elec(neelec)
      
      integer            :: iestag  ! whether e-field components are staggered or not
      
!      real               :: gxt(nx,4)
!      real               :: gyt(ny,4)
!      real               :: gzt(nz,4)

      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)
      
      real dxx(nx),dyy(ny),dzz(nz)         ! dx(i),dy(j),dz(k)
!      real gx(nx)
!      real gy(ny)
!      real gz(nz)
      
      character(len=*) microp
      
      integer   ipconc
      
      integer   :: bcx, bcy
!      include 'sam.index.ion.h'

!
!  INTEGERS
!
      integer    ia,na,na1,nba
      integer    ix
      integer    jy
      integer    kz
      integer    nx
      integer    ny
      integer    nz
      integer    nor, ng
      integer    jgs
      integer    i,j,k
      
      
      integer id1,jd1,kd1,istag,jstag,kstag,nstep
      
      integer ndebug
      parameter ( ndebug = 0 )
      
      
      integer nionstep,n
      
      logical istretch
      

      real cwccn
      real cimn,cimx
!      real cinccn(nz)
      
      integer ithomson
      parameter(ithomson=0)
      
      real alpha,alpha0
      parameter(alpha0 = 1.95e-12, alpha = 1.6e-12)
      real eta
      parameter(eta = 1.0e-11) ! coeff for large ion and small ion recombination, Takahashi 79
      real ezfairo
      parameter(ezfairo = -80.00)
      real efa1, efa2, efa3, efb1, efb2, efb3
      parameter(efa1 = 4.5e-3, efa2 = 3.8e-4, efa3 = 1.0e-4,       &
     &          efb1 = 0.50,   efb2 = 0.65,   efb3 = 0.10    )
     
      real ezfairs        ! fair field at surface
      parameter(  ezfairs = -100.0 )

      real ec,ecinv  ! fundamental unit of charge
      parameter (ec = 1.602e-19, ecinv=1./ec)
      
!      real emaxc ! point discharge (corona) threshold

      real ez,dezdz,zgrd,dezdz2

      real eperao
      parameter (eperao  = 8.8592e-12 )

      real poo, rd, cp, cap
      parameter( poo = 1.0e+05,       &
     &            rd = 287.04,       &
     &            cp = 1004.0,       &
     &           cap = rd/cp    )
      
      double precision temc,x,ffx,x1
      double precision xp,xn,tmp
      real cmax,climit,c1,c2,c3,ds,vel,s,c1t,c2t,c3t
      real cmaxall
      double precision pionpd,pionnd,attot
      double precision plionpd,plionnd
      


      integer ipass, mdnstp
      
      integer iixps(2),jjyqs(2),kkzrs(2)
      integer iiexs(2),jjeys(2),kkezs(2)
      integer*8 llworks(2)
      integer nnxs(2),nnys(2),nnzs(2)
      integer nbw,nbe,nbs,nbn,nztop
      integer ibg
!      real potfair(nztop),rhofair(nztop)
      

      integer ng1
      integer iunit,iunit0

      parameter(ng1 = 1)
      
      integer istop
      
      logical, parameter :: lcalce = .true.
!
!  REALS
!
!
      real       dt,dt1
      real       dx
      real       dy
      real       dz
      double precision       dv

      real       dtfac
      double precision       fac,fac1
      
      double precision delt1,delt2
!
!
!      real ab(nz,na)
      real ab(-nor+1:nz+nor,na)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real epot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,neelec)

!
!      real pb(nz)
!      real db(nz)
      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

!

      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ezw(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)


     
      real uz(nz,4)  ! ion mobilities
      real cion(nz,2) ! fairweather (base) ion concentrations

!

!      real fo(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real fu(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real fv(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)
!      real fw(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor)

      real tt0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) ! air temp (kelvin)
      real tt7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) ! cloud ice number conc.
      
      
!
! temporary arrays
!
      
!      real cwvt(nx,nz)
!      real civt(nx,nz)
!      real cirvt(nx,nz)
!      real cipvt(nx,nz)
!      real rwvt(nx,nz)
!      real fwvt(nx,nz)
!      real swvt(nx,nz)
!      real glvt(nx,nz)
!      real gmvt(nx,nz)
!      real ghvt(nx,nz)
!      real hlvt(nx,nz)
!      real hwvt(nx,nz)
      
      
      
      integer nzmax
      parameter (nzmax=410)

!      real vt(nx,nz,lqb:lqe)   ! terminal fall speed
      real alph(nzmax)
      double precision crg(nzmax)
      double precision recomb(nx,nz)
      double precision recomb12(nx,nz), recomb21(nx,nz)
      real chaffion(nx,nz)
      double precision chaffiontot,chaffioncount
      real, parameter :: achaff   = 1.4e-11 ! chaff corona constant (F m-1)
      real, parameter :: lchaff   = 0.1     ! chaff length (m)
      real, parameter :: ethchaff = 30.e3   ! chaff corona threshold (V m-1)
      real            :: chaffuzfac(nzmax)
      real            :: chaffionrate

      real,allocatable :: attach(:,:,:,:) ! (nx,nz,lqb:lqe,2)
!      real rad(nx,nz,lqb:lqe)  ! mean radius
      real,allocatable :: conc(:,:,:) ! (nx,nz,lqb:lqe) ! concentration
      real :: qmin(lqb:lqe)
!      real qmin(lqb:lqb+50)
      
      real,allocatable :: diff(:,:,:,:) ! (nx,nz,lqb:lqe,2)
      real,allocatable :: cond(:,:,:,:) ! (nx,nz,lqb:lqe,2)
      
      
      real test(nz,2)
      
      double precision scionfx(6,2), scionfxp,scionfxn
      double precision scionn1,scionn2,scionp1,scionp2,sciont1,sciont2
      double precision :: chglossp, chglossn
      double precision sciont0,scionp0,scionn0
      real sciont0s, sciont1s, sciont2s
      
      double precision chgneg,chgpos,chgnet ! total net negative/positive charge in storm (not including corona)
      
      integer ifirstcall
      data ifirstcall/0/
      
      save alph,ifirstcall,crg,chaffuzfac
      
      integer km1,kp1
      integer ioffset,count,countp,countn

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      
      real :: smin,smax
      double precision :: dpmin,dpmax
      
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      real*8  :: dpt1, dpt2

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag

      logical :: debug_mpi = .false.
      integer       :: nampi

! ################################################################
      
      ng = nor
      
#ifdef TAKON
!       return
#endif

      allocate( attach(nx,nz,lqb:lqe,2) )
      allocate( conc(nx,nz,lqb:lqe) )
      allocate( diff(nx,nz,lqb:lqe,2) )
      allocate( cond(nx,nz,lqb:lqe,2) )

      CALL cld_cpu('ION-MISC')

! erm: 2/13/2008 t7 and t8 are not used in iondrift anymore, so no need to copy
!      t7(1:nx-1,1:ny-1,1:nz-1) = an(1:nx-1,1:ny-1,1:nz-1,lscpi)
!      t8(1:nx-1,1:ny-1,1:nz-1) = an(1:nx-1,1:ny-1,1:nz-1,lscni)
      
      iunit = iunit0
!      iunit = 0
      IF ( ifirstcall .eq. 0 ) THEN
      
       ifirstcall = 1
       
       IF ( nz .gt. nzmax ) THEN
         write(6,*) 'STOP: increase nzmax in ionstep.F90'
         STOP
       ENDIF
!
! set up recombination coefficient (alpha) and cosmic ray generation
!
      
!      write(6,*) 'cosmic ray generation rates (kz,crg,x1)'
      IF ( my_rank == 0 ) THEn
       write(iunit,*) 'kz,crg(kz),ez,dezdz'
      ENDIF

      DO kz=1,nz-1
       IF ( ithomson .eq. 1 ) THEN
          temc = ab(kz,lt)       &
     &     *(pb(kz)/poo)**cap
        x = 2.43*(273.0/temc)**2        &
     &       * (pb(kz))/101300.0
        ffx = 1 - 4.0/(x**4)*(1 - Exp(-x)*(x+1))**2
        alph(kz) = alpha0*Sqrt((273.0/temc)**3)*ffx
       ELSE
        alph(kz) = alpha
       ENDIF
       
       chaffuzfac(kz) = achaff*lchaff*uz(kz,1)*uz(kz,2)/(uz(kz,1) + uz(kz,2))/ec

       zgrd = gzt(kz,1) 
       ez = ezfairo *       &
     &           ( efb1*exp(-efa1*zgrd)       &
     &            +efb2*exp(-efa2*zgrd)       &
     &            +efb3*exp(-efa3*zgrd) )
       dezdz = ezfairo *       &
     &           ( -efa1*efb1*exp(-efa1*zgrd)       &
     &             -efa2*efb2*exp(-efa2*zgrd)       &
     &             -efa3*efb3*exp(-efa3*zgrd) )
       dezdz2 = ezfairo *       &
     &           (   efa1**2*efb1*exp(-efa1*zgrd)       &
     &             + efa2**2*efb2*exp(-efa2*zgrd)       &
     &             + efa3**2*efb3*exp(-efa3*zgrd) )

        x1 = (1.4e-4)*ez*dezdz + dezdz**2 + ez*dezdz2
        
!        crg(kz) = uz(kz,1)*uz(kz,2)*eperao/ec
!     :   /( uz(kz,1) + uz(kz,2) )*x1
!     :    + alph(kz)*cion(kz,1)*cion(kz,2)

      IF ( fairweather ) THEN
        crg(kz) = alph(kz)*cion(kz,1)*cion(kz,2) ! * crgfac ! use crgfac instead in ion init.
      ELSE ! set to zero if 'fairweather' is false
        crg(kz) = 0.0
      ENDIF
      
      IF ( my_rank == 0 ) THEN
       write(iunit,*) kz,crg(kz),ez,dezdz
      ENDIF
        

!         write(6,*) kz,crg(kz),x1,ez,dezdz,dezdz2,
!     :   uz(kz,1)*uz(kz,2)*eperao/ec/( uz(kz,1) + uz(kz,2) )*x1 ,
!     :    alph(kz),cion(kz,1),cion(kz,2)
       ENDDO

      
      ENDIF ! ifirstcall

!
!  set minimum mass mixing ratios
!
      IF ( microp(1:5) .eq. 'ICE10' ) THEN
      call setqxmin(qmin)
!      qmin(lc) = 1.0e-09
!      qmin(li) = 1.0e-09
!      qmin(lr) = 1.0e-07
!      qmin(ls) = 1.0e-07
!      qmin(lf) = 1.0e-07
!      qmin(lgl) = 1.0e-07
!      qmin(lgm) = 1.0e-07
!      qmin(lgh) = 1.0e-07
!      qmin(lh) = 1.0e-07
!      qmin(lhl) = 1.0e-07
!      qmin(lir) = 1.0e-08
!      qmin(lip) = 1.0e-09
      ELSEIF ( microp(1:1) .eq. 'Z' ) THEN
       call setqxminz(qmin)
!       qmin(:) = 1.0e-7
!      qmin(4) = 1.0e-09
!      qmin(5) = 1.0e-07
!      qmin(6) = 1.0e-12
!      qmin(7) = 1.0e-07
!      qmin(8) = 1.0e-07
      ELSEIF ( microp(1:3) .eq. 'TAK' ) THEN
       qmin(:) = 5.*takcxmin
      ENDIF

      sciont0 = 0.0d0
      sciont0s = 0.0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      dpmax = 0.0d0
      dpmin = 0.0d0
      

#ifdef MPI
      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg

      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do kz = kzb,kze ; do jy = jyb,jye ; do ix = ixb,ixe
#else
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : sciont0,sciont0s,scionp0,scionn0)
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
#endif
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
       IF ( largeion ) THEN
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
       ENDIF
       
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv

       IF ( largeion ) THEN
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
       ENDIF
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv
       IF ( largeion ) THEN
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv
       ENDIF
       
       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
       IF ( largeion ) THEN
       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
       ENDIF
       
       dpmax = Max( dpmax, ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) )
       dpmin = Min( dpmin, ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) )
       
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont0
       mpitotindp(2)  = sciont0s
       mpitotindp(3)  = scionp0
       mpitotindp(4)  = scionn0

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont0 = mpitotoutdp(1)
       sciont0s = mpitotoutdp(2)
       scionp0 = mpitotoutdp(3)
       scionn0 = mpitotoutdp(4)
      ENDIF

       mpitotindp(1)  = dpmax
       mpitotindp(2)  = -dpmin

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, my_comm, mpi_error_code)

      IF ( my_rank == 0 ) THEN
       dpmax = mpitotoutdp(1)
       dpmin = -mpitotoutdp(2)
      ENDIF


#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'preion ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0


      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'preion single precision ion net charge (C):',       &
     &      sciont0s
     
      write(iunit,'(a,2(2x,1pe15.8))')        &
     &     'preion ion max/min charge (nC/m^3):',       &
     &      1.d9*dpmax,1.d9*dpmin
      ENDIF

!  Determine if multiple time steps are needed
!
      
      cmax = 0.0
#ifdef MPI
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg
#else
      ixe = nx-1
      jye = ny-1
      kze = nz-1
#endif
      IF ( iestag .eq. 1 ) THEN
       DO kz=1,kze
        DO jy=1,jye
         DO ix=1,ixe
          ds = dzz(kz)
          vel = Abs(Max(uz(kz,1),uz(kz,2))*elec(iez)%flt3d(ix,jy,kz))
          c1 = vel ! *dt1/ds
          c1t = vel*dt1/ds
          cmax = Max(cmax,c1*dt1/ds)

          ds = dyy(jy)
          vel = Abs(Max(uz(kz,1),uz(kz,2))*elec(iey)%flt3d(ix,jy,kz))
          c2 = vel ! *dt1/ds
          c2t = vel*dt1/ds
          cmax = Max(cmax,c2*dt1/ds)

          ds = dxx(ix)
          vel = Abs(Max(uz(kz,1),uz(kz,2))*elec(iex)%flt3d(ix,jy,kz))
          c3 = vel !*dt1/ds
          c3t = vel*dt1/ds
          cmax = Max(cmax,c3*dt1/ds)

          ds = Min(dzz(kz),Min(dyy(jy),dxx(ix)))
          cmax = Max( cmax, Sqrt(c1**2 + c2**2 + c3**2)*dt1/ds )
          ! cmax = Max( cmax, c1t+c2t+c3t )

         ENDDO
        ENDDO
       ENDDO
      ELSE
       DO kz=1,kze
        DO jy=1,jye
         DO ix=1,ixe
          ds = Min(dzz(kz),Min(dyy(jy),dxx(ix)))
          vel = Abs(uz(kz,2)*elec(iemag)%flt3d(ix,jy,kz))
          c1 = vel*dt1/ds
          cmax = Max(cmax,c1)
         ENDDO
        ENDDO
       ENDDO
      ENDIF

#ifdef MPI
!
! Get the global cmax
!
      CALL MPI_Allreduce(cmax, cmaxall, 1, MPI_REAL, MPI_MAX, my_comm, mpi_error_code)
  
      cmax = cmaxall


#endif 
      climit = 0.41
      IF ( cmax .lt. climit ) THEN
        nionstep = 1
        dt = dt1
      ELSE
        nionstep = Max(2, Int(cmax/climit) + 1)
        dt = dt1/nionstep
      ENDIF

      IF ( my_rank == 0 ) write(iunit,*) 'ION subtimesteps = ',nionstep,cmax


!      IF ( microp(1:3) /= 'TAK' ) THEN
!
!  remove charge from below-minimum mixing ratios (bulk schemes!) and put into ions
!  aps: since downscaling is handled in icezvd_dr, this is zero
!  aps: need to update using large ions...
!
      IF ( .false. ) THEN ! ions already created from small q
!      IF ( .true. ) THEN
      
      xp = 0.0d0
      xn = 0.0d0
      
      fac = Exp(-0.105*dt/5.0)
      
      
      count = 0
      countp = 0
      countn = 0
      chgneg = 0.0
      chgpos = 0.0
#ifndef MPI
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,ia,fac1,fac), REDUCTION(+ : xp,xn)
#endif
      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
        DO ia = lqb,lqe
          IF ( an(ix,jy,kz,ia) .lt. qmin(ia) ) THEN
          IF ( an(ix,jy,kz,lscb+ia-lqb) /= 0.0 ) THEN
!            count = count + 1
!            chgneg = chgneg + Min(0.0,  an(ix,jy,kz,lscb+ia-lqb))
!            chgpos = chgpos + Max(0.0,  an(ix,jy,kz,lscb+ia-lqb))
          ENDIF
            IF ( an(ix,jy,kz,lscb+ia-lqb) .gt. 0.0 ) THEN
!            countp = countp + 1
!             chgpos = chgpos + Max(0.0,  an(ix,jy,kz,lscb+ia-lqb))
              IF ( an(ix,jy,kz,lscb+ia-lqb) .gt. 1.0e-12 ) THEN
                fac1 = fac
              ELSE
                fac1 = 0.0
              ENDIF
              an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) +       &
     &             (1.0-fac1)*an(ix,jy,kz,lscb+ia-lqb)/ec
              xp = xp + (1.0-fac1)*an(ix,jy,kz,lscb+ia-lqb)/ec
              an(ix,jy,kz,lscb+ia-lqb) = fac1*an(ix,jy,kz,lscb+ia-lqb)
            ELSEIF ( an(ix,jy,kz,lscb+ia-lqb) .lt. 0.0 ) THEN
!            countn = countn + 1
!             chgneg = chgneg + Min(0.0,  an(ix,jy,kz,lscb+ia-lqb))
              IF ( an(ix,jy,kz,lscb+ia-lqb) .lt. -1.0e-12 ) THEN
                fac1 = fac
              ELSE
                fac1 = 0.0
              ENDIF
              an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) -       &
     &             (1.0-fac1)*an(ix,jy,kz,lscb+ia-lqb)/ec
              xn = xn + (1.0-fac1)*an(ix,jy,kz,lscb+ia-lqb)/ec
              an(ix,jy,kz,lscb+ia-lqb) = fac1*an(ix,jy,kz,lscb+ia-lqb)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = xp
       mpitotindp(2)  = xn

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       xp = mpitotoutdp(1)  
       xn = mpitotoutdp(2)  
      ENDIF
#endif

      
      IF ( my_rank == 0 ) write(iunit,'(a,4(2x,1pe13.5))') 'ion redistribute: xp,xn = ',xp,xn,xp*ec,xn*ec !,count,countp,countn,chgpos,chgneg
      
      ENDIF
      
!      ENDIF ! microp(1:3) /= 'TAK'

      
      chgneg = 0.0d0
      chgpos = 0.0d0
      t0(:,:,:) = 0.0
!      call zerond(nx,ny,nz,1,nor,t0)
      
      DO ia = lscb,lsceq
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO
 
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : chgpos, chgneg)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  
       
       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ENDIF
       
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv
      END IF

     ! lightning charge tendency, assuming iscnet has value from start of time step
      IF ( ichgtndlgt > 1 ) THEN
        elec(ichgtndlgt)%flt3d(ix,jy,kz) = elec(ichgtndlgt)%flt3d(ix,jy,kz) + (t0(ix,jy,kz) - elec(iscnet)%flt3d(ix,jy,kz) )
        elec(iscnet)%flt3d(ix,jy,kz) = t0(ix,jy,kz)
      ENDIF

      

      end do
      end do
      end do
#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &   'ionredistrib pos/neg/net charge (C):',       &
     &    chgpos,chgneg,(chgpos+chgneg)
      ENDIF
      
!      ENDIF

      CALL cld_cpu('ION-MISC')
      
      DO n=1,nionstep

      an(:,:,nz,lscpi) = cion(nz,1)
      an(:,:,nz,lscni) = cion(nz,2)

! check for periodic boundary conditions
      IF ( bcx .eq. 2 ) THEN

        elec(iex)%flt3d(nx,1:ny-1,1:nz-1)  = elec(iex)%flt3d(1,1:ny-1,1:nz-1) 

      DO ia=1,nor+1
      DO kz=1,kze
      DO jy=1,jye
        an(nx-1+ia,jy,kz,lscpi) = an(ia   ,jy,kz,lscpi)

        an(nx-1+ia,jy,kz,lscni) = an(ia   ,jy,kz,lscni)
      ENDDO
      ENDDO
      ENDDO

      DO ia=1,nor
      DO kz=1,kze
      DO jy=1,jye
        an(1-ia   ,jy,kz,lscpi) = an(nx-ia,jy,kz,lscpi)

        an(1-ia   ,jy,kz,lscni) = an(nx-ia,jy,kz,lscni)
      ENDDO
      ENDDO
      ENDDO

      ENDIF

#ifndef MPI
      IF ( bcy .eq. 2 ) THEN
      
        elec(iey)%flt3d(1:nx-1,ny,1:nz-1) = elec(iey)%flt3d(1:nx-1,1,1:nz-1) 
      
      DO ia=1,nor
      DO kz=1,kze
      DO ix=1,ixe
        an(ix,ny-1+ia,kz,lscpi) = an(ix,ia   ,kz,lscpi)

        an(ix,ny-1+ia,kz,lscni) = an(ix,ia   ,kz,lscni)
      ENDDO
      ENDDO
      ENDDO

      DO ia=1,nor
      DO kz=1,kze
      DO ix=1,ixe
        an(ix,1-ia   ,kz,lscpi) = an(ix,ny-ia,kz,lscpi)

        an(ix,1-ia   ,kz,lscni) = an(ix,ny-ia,kz,lscni)
      ENDDO
      ENDDO
      ENDDO

      ENDIF
#endif
      
!      IF ( .false. ) THEN
!      IF ( .true. ) THEN
      
      CALL cld_cpu('ION-DRIFT')
      
      IF ( ndebug .ge. 1 ) THEN
      DO kz=1,nz
       test(kz,1) = an(nx/2,ny/2,kz,lscpi) 
       test(kz,2) = an(nx/2,ny/2,kz,lscni) 
       write(iunit,'(a,i3,5(2x,1pe12.5))')        &
     &  'kz,test: ',kz,gzt(kz,1),       &
     &   an(nx/2,ny/2,kz,lscpi)-cion(kz,1), cion(kz,1),       &
     &   an(nx/2,ny/2,kz,lscni)-cion(kz,2), cion(kz,2)
      ENDDO
      ENDIF


      IF ( ndebug .ge. 1 ) THEN
      DO kz=nz-1,1,-1
        write(iunit,'(i3,10(2x,1pe12.5))') kz,       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*an(nx/2,ny/2,kz,lscpi)*ec,       &
     &  uz(kz,2)*elec(iez)%flt3d(nx/2,ny/2,kz)*an(nx/2,ny/2,kz,lscni)*ec,       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &   (an(nx/2,ny/2,kz,lscpi))*ec,       &
     &  uz(kz,2)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &   (an(nx/2,ny/2,kz,lscni) )*ec,       &
     & 0.5*(uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz) +        &
     &     uz(kz+1,1)*elec(iez)%flt3d(nx/2,ny/2,kz+1))*       &
     &     (an(nx/2,ny/2,kz,lscpi)),       &
     & 0.25*(uz(kz,1)+uz(kz+1,1))*(elec(iez)%flt3d(nx/2,ny/2,kz) +        &
     &     elec(iez)%flt3d(nx/2,ny/2,kz+1))*       &
     &    (an(nx/2,ny/2,kz,lscpi)),       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &  (an(nx/2,ny/2,kz,lscpi)),       &
     &  ezfair(kz)*cion(kz,1)*uz(kz,1),        &
     &  ezfair(kz)*       &
     &  (an(nx/2,ny/2,kz,lscpi))*uz(kz,1),       &
     &  ezfairw(kz)*       &
     &  (an(nx/2,ny/2,kz,lscpi))*uz(kz,3)
      ENDDO
      ENDIF


      DO jy=1,jye
      DO ix=1,ixe

          kz = 0
          an(ix,jy,kz,lscpi) = cion(1,1) 
          an(ix,jy,kz,lscni) = cion(1,2)

!          kz = nz-1
!          an(ix,jy,kz,lscpi) = cion(nz,1)
!          an(ix,jy,kz,lscni) = cion(nz,2)
!          ad(ix,jy,kz,lscpi) = cion(nz,1)
!          ad(ix,jy,kz,lscni) = cion(nz,2)

      ENDDO
      ENDDO

      mdnstp = mod(nstep,2)

      scionn1 = 0.0d0
      scionp1 = 0.0d0
      sciont1 = 0.0d0
      sciont1s = 0.0


!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : sciont1,sciont1s,scionp1,scionn1)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       
       IF ( largeion ) THEN
       sciont1 = sciont1 + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv        &
     &                   + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv 
      
       sciont1s = sciont1s + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &                     + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv 

       scionp1 = scionp1 + ec*(an(ix,jy,kz,lscpi))*dv       &
     &                   + ec*(an(ix,jy,kz,lscpli))*dv

       scionn1 = scionn1 + ec*( - an(ix,jy,kz,lscni))*dv       &
     &                   + ec*(an(ix,jy,kz,lscnli))*dv
      ELSE
       sciont1 = sciont1 + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv 
      
       sciont1s = sciont1s + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv

       scionp1 = scionp1 + ec*(an(ix,jy,kz,lscpi))*dv

       scionn1 = scionn1 + ec*( - an(ix,jy,kz,lscni))*dv
      ENDIF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont1
       mpitotindp(2)  = sciont1s
       mpitotindp(3)  = scionp1
       mpitotindp(4)  = scionn1

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont1 = mpitotoutdp(1)
       sciont1s = mpitotoutdp(2)
       scionp1 = mpitotoutdp(3)
       scionn1 = mpitotoutdp(4)
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'predrift ion pos/neg/net charge (C):',       &
     &      scionp1,scionn1,sciont1

      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'predrift single pres. ion net charge (C):',       &
     &      sciont1s
      
      ENDIF


#if defined ( BOXMG ) && defined ( MPI )
! run a few fine-grid relaxations to try to reduce noise at the boundaries.
      IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep call putf'
         
       IF ( .false. ) THEN
       
!         IF ( my_rank == 0 ) THEN
!           DO k = 1,NLzdg
!             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!           ENDDO
!         ENDIF
      CALL cld_cpu('BOXMG')
         
!         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,     
!     :               nbw,nbe,nbs,nbn,
!     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!     :               iGsdg, jGsdg, kGsdg, dx, dy, dz, 0
!     :              )
!         IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep: done putf'


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dpt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

       
       
         IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep: call BMG3_SymStd_SOLVE_boxmg'

        i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)
        BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 1
        BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 3

!         BMG_iPARMSdg(id_BMG3_NRELAX_FG) = 2
         CALL commas_SymStd_SOLVE_boxmgdg(                                          &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg,       &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                          &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,        &
     &             SOdg, NSOdg,                                                     &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                          &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                            &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,                      &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                                 &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                                 &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                                 &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                                 &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                         &
     &             pMSGdg, pMSGSOdg,                                                &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,                &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg,       &
     &             my_comm                                               &
     &             )

!         BMG_iPARMSdg(id_BMG3_NRELAX_FG) = 0

        BMG_iPARMSdg(id_BMG3_NRELAX_UP )   = 1
        BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i

         IF ( ndebug .ge. 1 ) write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dpt2 = MPI_Wtime()

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for boxmg = ', dpT2 - dpT1
            WRITE(0,*) ' time for boxmg (ionstep) = ', dpT2 - dpT1
         ENDIF

!      IF ( my_rank .eq. 0 ) THEN
!          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!     :              )
!      ENDIF

!
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,            &
     &               gxt,gyt,gzt,       &
     &               nbw,nbe,nbs,nbn,iunit,       &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz       &
     &              )

      CALL cld_cpu('BOXMG')

      ENDIF ! true or false
#endif

#ifdef MPI

      IF ( number_of_processes .gt. 1 ) THEN

       CALL cld_cpu('MPICOM-IONSTEP')

       nampi = 2


       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,an(-ng+1,-ng+1,-ng+1,lscpi) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,an(-ng+1,-ng+1,-ng+1,lscpi))

       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,an(-ng+1,-ng+1,-ng+1,lscpi))

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,an(-ng+1,-ng+1,-ng+1,lscpi))


       nampi = 3

       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,epot(-ng+1,-ng+1,-ng+1,1) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,epot(-ng+1,-ng+1,-ng+1,1))

       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,epot(-ng+1,-ng+1,-ng+1,1))

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,         &
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,epot(-ng+1,-ng+1,-ng+1,1))


        CALL cld_cpu('MPICOM-IONSTEP')

       ENDIF

#endif

      IF ( iestag .eq. 0 ) THEN
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,km1)
      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
       km1 = Max( 1, kz-1 )
!       IF ( .true. ) THEN
       IF ( kz .gt. 1 ) THEN
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-elec(ipot)%flt3d(ix,jy,km1))*gzt(kz,3)
       ELSE
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-0.0)*gzt(kz,3)
       
       ENDIF

!       ELSE
!       IF ( kz .eq. 1 ) THEN
!c       ezw(ix,jy,kz) = ezfairw(kz)
!        ezw(ix,jy,kz) = elec(ix,jy,kz,iez)
!       ELSE
!       ezw(ix,jy,kz) = 0.5*(elec(ix,jy,kz,iez)+elec(ix,jy,kz-1,iez)) -
!     :     0.5*(elec(ix,jy,kz,iezb)+elec(ix,jy,kz-1,iezb)) +
!     :     ezfairw(kz)
!       ENDIF
!          elec(ix,jy,kz,iez) = elec(ix,jy,kz,iez) - elec(ix,jy,kz,iezb)
!     :       + ezfair(kz)
!
!        ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE
       DO kz = -ng+1,nz+ng
!         km1 = Max( 1   ,kz-1 )
!         kp1 = Min( nz-1,kz   )
         DO jy = -ng+1,ny+ng
         DO ix = -ng+1,nx+ng
! erm 7.31.2007: restored staggered value of ezw
         ezw(ix,jy,kz) = elec(iez)%flt3d(ix,jy,kz)
!         ezw(ix,jy,kz) = 0.5*(elec(iez)%flt3d(ix,jy,km1) + elec(iez)%flt3d(ix,jy,kp1))
!         IF ( Sign(1.0,ezw(ix,jy,kz)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,kz)) ) THEN
!           ezw(ix,jy,kz) = 0.5*elec(iez)%flt3d(ix,jy,kz)
!         ENDIF
         ENDDO
         ENDDO
       ENDDO
      ENDIF

      
!      DO ipass = 1,0
      DO ipass = 1,2

!      IF ( .true. .and. nstep .gt. 5 ) THEN



      IF ( ipass .eq. (1+mdnstp) ) THEN
!
!  ion drift motion, positive ions
! 
      s = 1.0
      na1 = lscpi
      
      IF ( iondriftorder == 6 ) THEN
      call iondrift       &
     & (nx,ny,nz,nor,na,na1,s,emaxc,icorona,nstep,       &
     &  dt,dx,dy,dz,cion(1,1),scionfx(1,1),       &
     &  ab,t7,an,       &
     &  elec,uz(1,1),ezfair,  &
     &  t0,ezw,       &
     &  t1,t2,t3,       &
     &  dxx,dyy,dzz,gxt(-nor+1,3),gyt(-nor+1,3),gzt(-nor+1,3),iestag,bcx,bcy) 

      ELSEIF ( iondriftorder < 0 ) THEN
!         IF ( my_rank == 0 ) THEN
!           write(iunit,*) 'ion drift is off: iondriftorder = ',iondriftorder
!         ENDIF
      ELSE
      
      call iondrift1st       &
     & (nx,ny,nz,nor,na,na1,s,emaxc,icorona,nstep,       &
     &  dt,dx,dy,dz,cion(1,1),scionfx(1,1),       &
     &  ab,t7,an,       &
     &  elec,uz(1,1),ezfair,  &
     &  t0,ezw,       &
     &  t1,t2,t3,       &
     &  dxx,dyy,dzz,gxt(-nor+1,3),gyt(-nor+1,3),gzt(-nor+1,3),iestag,bcx,bcy) 

      ENDIF
      
      ELSE
!
!  ion drift motion, negative ions
! 
      s = -1.0
      na1 = lscni
      
      IF ( iondriftorder == 6 ) THEN
      call iondrift       &
     & (nx,ny,nz,nor,na,na1,s,emaxc,icorona,nstep,       &
     &  dt,dx,dy,dz,cion(1,2),scionfx(1,2),       &
     &  ab,t8,an,       &
     &  elec,uz(1,2),ezfair,  &
     &  t0,ezw,       &
     &  t1,t2,t3,       &
     &  dxx,dyy,dzz,gxt(-nor+1,3),gyt(-nor+1,3),gzt(-nor+1,3),iestag,bcx,bcy) 

      ELSEIF ( iondriftorder < 0 ) THEN
!         IF ( my_rank == 0 ) THEN
!           write(iunit,*) 'ion drift is off: iondriftorder = ',iondriftorder
!         ENDIF

      ELSE
      call iondrift1st       &
     & (nx,ny,nz,nor,na,na1,s,emaxc,icorona,nstep,       &
     &  dt,dx,dy,dz,cion(1,2),scionfx(1,2),       &
     &  ab,t8,an,       &
     &  elec,uz(1,2),ezfair,  &
     &  t0,ezw,       &
     &  t1,t2,t3,       &
     &  dxx,dyy,dzz,gxt(-nor+1,3),gyt(-nor+1,3),gzt(-nor+1,3),iestag,bcx,bcy) 

      ENDIF
      
      ENDIF
      
!      ENDIF
! solve potential/electric field 
!
!

!      IF ( .true. ) THEN
!      IF ( ipass .eq. 1 .or.
!     :     (nionstep .gt. 1 .and. n .lt. nionstep) ) THEN 
!      IF ( ipass .eq. 2 .or. ipass .eq. 0 .or. dz .gt. 250. ) THEN
      IF ( ipass .eq. 2 .or. ipass .eq. 0 ) THEN

      CALL cld_cpu('ION-DRIFT')

      CALL cld_cpu('ION-EFIELD')

!
!  preserve ion densities at nz-1
!
!      DO jy=1,ny-1
!      DO ix=1,nx-1
!          kz = nz-1
!          an(ix,jy,kz,lscpi) = cion(kz,1) 
!          an(ix,jy,kz,lscni) = cion(kz,2) 
!      ENDDO
!      ENDDO

      t0(:,:,:) = 0.0
      
      DO ia = lscb,lsceq
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO
 
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,ia)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
!        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
!        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))

       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))         &
     &  + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  
       ENDIF
      end do
      end do
      end do


      IF ( ndebug .ge. 1 ) write(iunit,*) 'call mud'

#if defined ( BOXMG ) && defined ( MPI )

      IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep call putf'
         
!         IF ( my_rank == 0 ) THEN
!           DO k = 1,NLzdg
!             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!           ENDDO
!         ENDIF
      CALL cld_cpu('BOXMG')
         
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,            &
     &               nbw,nbe,nbs,nbn,       &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 0       &
     &              )
         IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep: done putf'


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dpt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

         IF ( ndebug .ge. 1 ) write(iunit,*) 'ionstep: call BMG3_SymStd_SOLVE_boxmg'

         CALL commas_SymStd_SOLVE_boxmgdg(                                          &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg,       &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                          &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,        &
     &             SOdg, NSOdg,                                                     &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                          &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                            &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,                      &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                                 &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                                 &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                                 &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                                 &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                         &
     &             pMSGdg, pMSGSOdg,                                                &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,                &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg,       &
     &             my_comm                                               &
     &             )

         IF ( ndebug .ge. 1 ) write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dpt2 = MPI_Wtime()

         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) ' time for boxmg = ', dpT2 - dpT1
         ENDIF

!      IF ( my_rank .eq. 0 ) THEN
!          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!     :              )
!      ENDIF

!
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,            &
     &               gxt,gyt,gzt,       &
     &               nbw,nbe,nbs,nbn,iunit,       &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz       &
     &              )

      CALL cld_cpu('BOXMG')


#else

      IF ( number_of_processes .eq. 1 ) THEN

      CALL mud(nx,ny,nz,nor,       &
     &   dx,dy,dz,gxt,gyt,gzt,       &
     &   nnxs(1),nnys(1),nnzs(1),       &
     &   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1),       &
     &   nbw,nbe,nbs,nbn,       &
     &   llworks(1),istretch,       &
     &   t0,ibg,       &
     &   elec, iestag, iunit, lcalce, bcx, bcy)
!
!
      IF ( ndebug .ge. 1 ) write(iunit,*) 'call efield'

      
      IF ( iestag .eq. 0 ) THEN
      call efield       &
     &  (nx,ny,nz       &
     &  ,id1,jd1,kd1,istag,jstag,kstag       &
     &  ,dx,dy,dz,elec       &
     &  ,gxt(1,3),gyt(1,3),gzt(1,3)       &
     &  ,iunit)    
       ENDIF

       
         IF ( ndebug .ge. 1 ) write(iunit,*) 'done efield'

      ENDIF ! number of processes

#endif

      IF ( iex2 > 1 ) THEN
      call efield2       &
     &  (nx,ny,nz       &
     &  ,id1,jd1,kd1,istag,jstag,kstag       &
     &  ,dx,dy,dz,elec       &
     &  ,gxt(1,3),gyt(1,3),gzt(1,3)       &
     &  ,iunit)   
      ENDIF

      IF ( iestag .eq. 0 ) THEN 
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
      
!       IF ( .true. ) THEN
       IF ( kz .gt. 1 ) THEN
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-elec(ipot)%flt3d(ix,jy,kz-1))*gzt(kz,3)
       ELSE
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-0.0)*gzt(kz,3)
       
       ENDIF
       
      ENDDO
      ENDDO
      ENDDO
      
      ELSE

!       ezw(1:nx-1,1:ny-1,1:nz) = elec(iez)%flt3d(1:nx-1,1:ny-1,1:nz)

       DO kz = -ng+1,nz+ng
!         km1 = Max( 1   ,kz-1 )
!         kp1 = Min( nz-1,kz   )
         DO jy = -ng+1,ny+ng
         DO ix = -ng+1,nx+ng
! erm 7.31.2007: restored staggered value of ezw
         ezw(ix,jy,kz) = elec(iez)%flt3d(ix,jy,kz)
!         ezw(ix,jy,kz) = 0.5*(elec(iez)%flt3d(ix,jy,km1) + elec(iez)%flt3d(ix,jy,kp1))
!         IF ( Sign(1.0,ezw(ix,jy,kz)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,kz)) ) THEN
!           ezw(ix,jy,kz) = 0.5*elec(iez)%flt3d(ix,jy,kz)
!         ENDIF
         ENDDO
         ENDDO
       ENDDO

      ENDIF

      IF ( ndebug .ge. 1 ) THEN
      write(iunit,*) 'print efield'
      DO kz=nz-1,1,-1
        write(iunit,'(i3,3(2x,1pe12.5))') kz,       &
     &   elec(iez)%flt3d(nx/2,ny/2,kz),ezfair(kz),        &
     &   elec(iez)%flt3d(nx/2,ny/2,kz)-ezfair(kz)
      ENDDO
      ENDIF

      CALL cld_cpu('ION-EFIELD')

      CALL cld_cpu('ION-DRIFT')

      ENDIF ! ipass .eq. 1 .or. ...


      ENDDO ! ipass

      scionn2 = 0.0d0
      scionp2 = 0.0d0
      sciont2 = 0.0d0
      sciont2s = 0.0

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : sciont2,sciont2s,scionp2,scionn2)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)

       IF ( largeion ) THEN
       sciont2 = sciont2        &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
       sciont2s = sciont2s        &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
       scionp2 = scionp2        &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv
       scionn2 = scionn2       &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      ELSE
       sciont2 = sciont2        &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
       sciont2s = sciont2s        &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
       scionp2 = scionp2        &
     &  + ec*(an(ix,jy,kz,lscpi))*dv
       scionn2 = scionn2       &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF

      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont2
       mpitotindp(2)  = sciont2s
       mpitotindp(3)  = scionp2
       mpitotindp(4)  = scionn2
       k = 4
       DO j=1,2
        DO i=1,6
          k = k + 1
          mpitotindp(k) = scionfx(i,j)
        ENDDO
       ENDDO

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, k, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont2 = mpitotoutdp(1)
       sciont2s = mpitotoutdp(2)
       scionp2 = mpitotoutdp(3)
       scionn2 = mpitotoutdp(4)

       k = 4
       DO j=1,2
        DO i=1,6
          k = k + 1
          scionfx(i,j) = mpitotoutdp(k) 
        ENDDO
       ENDDO

      ENDIF
#endif

      CALL cld_cpu('ION-DRIFT')

      CALL cld_cpu('ION-MISC')

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postdrift ion pos/neg/net charge (C):',       &
     &      scionp2,scionn2,sciont2
      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'postdrift single pres. ion net charge (C):',       &
     &      sciont2s

      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'difference ion pos/neg/net charge (C):',       &
     &      (scionp2-scionp1),       &
     &      (scionn2-scionn1),       &
     &      (sciont2-sciont1)

      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'difference single pres. ion net charge (C):',       &
     &      (sciont2s-sciont1s)
     


      
      scionfxp = ec*       &
     &         (scionfx(1,1) +  scionfx(2,1) + scionfx(3,1) +       &
     &          scionfx(4,1) +  scionfx(5,1) + scionfx(6,1) )

      scionfxn = ec*       &
     &         (scionfx(1,2) +  scionfx(2,2) + scionfx(3,2) +       &
     &          scionfx(4,2) +  scionfx(5,2) + scionfx(6,2) )
      
      write(iunit,'(a,5(1x,1pe15.8))')        &
     & 'Net Ion charge fluxes (domain,top,bottom): ',       &
     &  scionfxp, scionfxn, scionfxp - scionfxn,       &
     &  ec*(scionfx(6,1)-scionfx(6,2)),       &
     &  ec*(scionfx(5,1)-scionfx(5,2))

      write(iunit,'(a,5(1x,1pe15.8))')        &
     & 'Net Lateral Ion charge fluxes (w,e,s,n): ',       &
     &  ec*(scionfx(1,1)-scionfx(1,2)),       &
     &  ec*(scionfx(2,1)-scionfx(2,2)),       &
     &  ec*(scionfx(3,1)-scionfx(3,2)),       &
     &  ec*(scionfx(4,1)-scionfx(4,2))      
      
      write(iunit,'(a,5(1x,1pe15.8))')        &
     & 'Neg Lateral Ion charge fluxes (w,e,s,n): ',       &
     &  ec*(scionfx(1,2)),       &
     &  ec*(scionfx(2,2)),       &
     &  ec*(scionfx(3,2)),       &
     &  ec*(scionfx(4,2))      

      write(iunit,'(a,5(1x,1pe15.8))')        &
     & 'Pos Lateral Ion charge fluxes (w,e,s,n): ',       &
     &  ec*(scionfx(1,1)),       &
     &  ec*(scionfx(2,1)),       &
     &  ec*(scionfx(3,1)),       &
     &  ec*(scionfx(4,1))      

      ENDIF

!      write(iunit,*) 'Ex(west),Ex(east),Ey(south),Ey(North)'
!      DO kz = nz-kd1,1,-1
!      write(iunit,'(i4,4(1x,1pe12.4))') kz,elec(1,ny/2,kz,iex),
!     :  elec(nx-1,ny/2,kz,iex), elec(nx/2,1,kz,iey),
!     :  elec(nx/2,ny-1,kz,iey)
!      
!      ENDDO
      
      DO jy=1,jye
      DO ix=1,ixe
         IF ( kzbeg .eq. nzbeg ) THEN
          kz = 0
          an(ix,jy,kz,lscpi) = cion(1,1) ! 0.0 
          an(ix,jy,kz,lscni) = cion(1,2) ! 0.0 
         ENDIF

         IF ( kzend .eq. nzend ) THEN
          kz = kze
          an(ix,jy,kz,lscpi) = cion(nz,1) 
          an(ix,jy,kz,lscni) = cion(nz,2) 
         ENDIF
      ENDDO
      ENDDO


      sciont0 = 0.0d0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      sciont0s = 0.0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : sciont0,sciont0s,scionp0,scionn0)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       
       IF ( largeion ) THEN
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv

       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      
      ELSE
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv

       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF
      
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont0
       mpitotindp(2)  = sciont0s
       mpitotindp(3)  = scionp0
       mpitotindp(4)  = scionn0

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont0 = mpitotoutdp(1)
       sciont0s = mpitotoutdp(2)
       scionp0 = mpitotoutdp(3)
       scionn0 = mpitotoutdp(4)
      ENDIF
#endif
      
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postnzreset ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0
       
      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'postnzreset single pres. ion net charge (C):',       &
     &      sciont0s
      ENDIF
       
       
      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
      ENDDO
      ENDDO
      ENDDO


      sciont0 = 0.0d0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      sciont0s = 0.0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv), &
!$OMP REDUCTION (+ : sciont0,sciont0s,scionp0,scionn0)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       
       IF ( largeion ) THEN
       
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv

       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      
      ELSE
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv

       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF
      
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont0
       mpitotindp(2)  = sciont0s
       mpitotindp(3)  = scionp0
       mpitotindp(4)  = scionn0

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont0 = mpitotoutdp(1)
       sciont0s = mpitotoutdp(2)
       scionp0 = mpitotoutdp(3)
       scionn0 = mpitotoutdp(4)
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postmaxzero ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0

      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postmaxzero single pres. ion net charge (C):',       &
     &      sciont0s

      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'difference2 ion pos/neg/net charge (C):',       &
     &      (scionp0-scionp1),       &
     &      (scionn0-scionn1),       &
     &      (sciont0-sciont1)

      write(iunit,'(a,(2x,1pe15.8))')        &
     &     'difference2 single pres. ion net charge (C):',       &
     &      (sciont0s-sciont1s)
      ENDIF

      IF ( ndebug .ge. 1 ) THEN
      
      DO kz=nz-1,1,-1
       write(iunit,'(a,i3,7(2x,1pe12.5))') '       &
     &  kz,test: ',kz,gzt(kz,1),       &
     &   test(kz,1),test(kz,2),       &
     &  -test(kz,1)+(an(nx/2,ny/2,kz,lscpi)),       &
     &  -test(kz,2)+(an(nx/2,ny/2,kz,lscni)),       &
     &  ec*(an(nx/2,ny/2,kz,lscpi)-an(nx/2,ny/2,kz,lscni)),       &
     &  elec(iez)%flt3d(nx/2,ny/2,kz)
      ENDDO

      DO kz=nz-1,1,-1
        write(iunit,'(i3,9(2x,1pe12.5))') kz,       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*an(nx/2,ny/2,kz,lscpi)*ec,       &
     &  uz(kz,2)*elec(iez)%flt3d(nx/2,ny/2,kz)*an(nx/2,ny/2,kz,lscni)*ec,       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &   (an(nx/2,ny/2,kz,lscpi))*ec,       &
     &  uz(kz,2)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &   (an(nx/2,ny/2,kz,lscni) )*ec,       &
     & 0.5*(uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz) +        &
     &     uz(kz+1,1)*elec(iez)%flt3d(nx/2,ny/2,kz+1))*       &
     &     (an(nx/2,ny/2,kz,lscpi)),       &
     & 0.25*(uz(kz,1)+uz(kz+1,1))*(elec(iez)%flt3d(nx/2,ny/2,kz) +        &
     &     elec(iez)%flt3d(nx/2,ny/2,kz+1))*       &
     &    (an(nx/2,ny/2,kz,lscpi)),       &
     &  uz(kz,1)*elec(iez)%flt3d(nx/2,ny/2,kz)*       &
     &  (an(nx/2,ny/2,kz,lscpi)),       &
     &  ezfair(kz)*cion(kz,1)*uz(kz,1),        &
     &  ezfair(kz)*       &
     &  (an(nx/2,ny/2,kz,lscpi))*uz(kz,1)
      ENDDO
      
      ENDIF ! ndebug

!      DO kz=1,nz-1
!      DO jy=1,ny-1
!      DO ix=1,nx-1
!          elec(ix,jy,kz,iez) = elec(ix,jy,kz,iez) + elec(ix,jy,kz,iezb)
!     :       - ezfair(kz)
!      ENDDO
!      ENDDO
!      ENDDO

      IF ( ndebug .ge. 1 ) THEN
      
      DO kz=nz-1,1,-1
        write(iunit,'(i3,3(2x,1pe12.5))') kz,       &
     &   elec(iez)%flt3d(nx/2,ny/2,kz),ezfair(kz),       &
     &   elec(iez)%flt3d(nx/2,ny/2,kz)-ezfair(kz)
      ENDDO
      
      ENDIF ! ndebug

       
!       ENDIF ! .false.

       IF ( ndebug .ge. 1 ) THEN
       
       write(iunit,*) 'post drift values'
      DO kz=nz,1,-1
       write(iunit,'(a,i3,5(2x,1pe12.5))')        &
     &  'kz,test: ',kz,gzt(kz,1),       &
     &   an(nx/2,ny/2,kz,lscpi)-cion(kz,1), cion(kz,1),       &
     &   an(nx/2,ny/2,kz,lscni)-cion(kz,2), cion(kz,2)
      ENDDO
      
      ENDIF ! ndebug
!
! Check net charge
!
      chgneg = 0.0d0
      chgpos = 0.0d0
      t0(:,:,:) = 0.0
      sciont0 = 0.0d0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      sciont0s = 0.0
      chglossp = 0.0d0
      chglossn = 0.0d0
      
      DO ia = lscb,lsceq
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      end do
      ENDDO
 
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv),   &
!$OMP  REDUCTION (+ : chgpos, chgneg, sciont0, sciont0s, scionp0, scionn0)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)
       
        chglossp = chglossp + Min( 0.0, an(ix,jy,kz,lscpi))*dv
        chglossn = chglossn + Min( 0.0, an(ix,jy,kz,lscni))*dv
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
       
       IF ( largeion ) THEN
       
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))         &
     &  + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  

       ENDIF
       
      IF ( t0(ix,jy,kz) .gt. 0.0 ) THEN
        chgpos = chgpos + t0(ix,jy,kz)*dv
      ELSE
        chgneg = chgneg + t0(ix,jy,kz)*dv
      END IF

       IF ( largeion ) THEN
       
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) *dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) *dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      
      ELSE
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF

      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg
       mpitotindp(3)  = sciont0
       mpitotindp(4)  = sciont0s
       mpitotindp(5)  = scionp0
       mpitotindp(6)  = scionn0
       mpitotindp(7)  = chglossp
       mpitotindp(8)  = chglossn

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 8, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  
       sciont0 = mpitotoutdp(3)
       sciont0s = mpitotoutdp(4)
       scionp0 = mpitotoutdp(5)
       scionn0 = mpitotoutdp(6)
       chglossp = mpitotindp(7)
       chglossn = mpitotindp(8)
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &   'postdrift pos/neg/net charge (C):',       &
     &    chgpos,chgneg,(chgpos+chgneg)
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postdrift ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postdrift ion charge loss (C):',       &
     &      chglossp*ec,chglossn*ec, (chglossp+chglossn)*ec
      ENDIF


      CALL cld_cpu('ION-MISC')

!
! do other work on x-z slices (for smaller temporary arrays)
!
      chaffiontot = 0.0d0
      
      CALL cld_cpu('ION-ATTACH')

       istop = 0
! C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(ix,jy,kz,jgs,ia,attach,diff,cond,conc,
! C$OMP+      recomb, x, pionpd, pionnd, attot )
      DO jy=1,jye

      jgs = jy

      DO ia=lqb,lqe
       DO kz=1,kze
        DO ix=1,ixe
          attach(ix,kz,ia,1) = 0.0
          attach(ix,kz,ia,2) = 0.0
!          diff(ix,kz,ia,1) = 0.0
!          diff(ix,kz,ia,2) = 0.0
!          cond(ix,kz,ia,1) = 0.0
!          cond(ix,kz,ia,2) = 0.0
          conc(ix,kz,ia) = 0.0
        ENDDO
       ENDDO
      ENDDO
!
!  ion attachment by diffusion and conduction
!
       IF ( microp(1:5) .eq. 'ICE10' ) THEN

       call ionattdd(nx,ny,nz,nor,na,nba,dt,dz,jgs,       &
     &  elec(iez)%flt3d,attach,uz,diff,cond,conc,       &
     &  an,dn,db,tt0,tt7,       &
     &  qmin)

       ELSEIF ( microp(1:1) .eq. 'Z' ) THEN

       call ionattz(nx,ny,nz,nor,na,nba,dt,dz,jgs,       &
     &  elec(iez)%flt3d,attach,uz,diff,cond,conc,       &
     &  an,dn,db,tt0,tt7,       &
     &  qmin)
#ifdef TAKON
       ELSEIF ( microp(1:3) .eq. 'TAK' ) THEN

       call ionatttak(nx,ny,nz,nor,na,nba,dt,dz,jgs,       &
     &  elec(iez)%flt3d,attach,uz,diff,cond,conc,       &
     &  an,dn,db,tt0,tt7,       &
     &  qmin)
#endif
       ENDIF


!
!  recombination
!

      DO kz=1,nz-2
      DO ix=1,ixe
        recomb(ix,kz) = alph(kz)*       &
     &         Max(0.0,an(ix,jy,kz,lscpi)*an(ix,jy,kz,lscni))
       IF ( largeion ) THEN
        recomb21(ix,kz) = eta*       &
     &         Max(0.0,an(ix,jy,kz,lscpli)*an(ix,jy,kz,lscni))
        recomb12(ix,kz) = eta*       &
     &         Max(0.0,an(ix,jy,kz,lscnli)*an(ix,jy,kz,lscpi))
       ELSE
        recomb12(ix,kz) = 0.0 ; recomb21(ix,kz) = 0.0
       ENDIF
      ENDDO
      ENDDO
      
      chaffioncount = 0.0d0
      IF ( lnchaff .gt. 1 .or. chaffconc > 0.0) THEN
        DO kz=1,nz-2
          DO ix=1,ixe
           IF ( elec(iemag)%flt3d(ix,jy,kz) .gt. ethchaff .and. ( an(ix,jy,kz,lnchaff) .gt. 0. .or. chaffconc > 0.) ) THEN
             dv = dxx(ix)*dyy(jy)*dzz(kz)
             chaffionrate = chaffuzfac(kz)*(elec(iemag)%flt3d(ix,jy,kz) - ethchaff)**2
             IF ( lnchaff > 1 ) THEN
               chaffion(ix,kz) = chaffionrate*an(ix,jy,kz,lnchaff) ! already divided by ec, and lnchaff is concentration, not mass density
             ELSE
               chaffionrate = chaffuzfac(kz)*( Min(5.e3, Abs(elec(iemag)%flt3d(ix,jy,kz) - ethchaff) )  )**2
               chaffion(ix,kz) = chaffionrate*chaffconc
             ENDIF
!             write(iunit,*) 'ix,kz,chaffion = ',ix,kz,chaffion(ix,kz),0.001*elec(iemag)%flt3d(ix,jy,kz)
             an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + chaffion(ix,kz)*dt
             an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) + chaffion(ix,kz)*dt
             chaffioncount = chaffioncount + chaffion(ix,kz)*dt*dv
           ELSE
             chaffion(ix,kz) = 0.0
           ENDIF
          ENDDO
        ENDDO
      
        chaffiontot = chaffiontot + chaffioncount
      
      ELSE
        chaffion(:,:) = 0.0
      ENDIF
      
      
      DO kz=1,nz-2
       IF ( ndebug .ge. 1 ) THEN
       IF ( jy .eq. ny/2) THEN
!$OMP CRITICAL
       write(iunit,'(a,i3,5(2x,1pe12.5))')        &
     &  'nz-kz,recomb,crg,diff: ',nz-kz,recomb(nx/2,nz-kz),crg(nz-kz),       &
     &  recomb(nx/2,nz-kz) - crg(nz-kz),       &
     &  attach(nx/2,nz-kz,lqb,1),attach(nx/2,nz-kz,lqb,2)
!$OMP END CRITICAL
       ENDIF
       ENDIF
      DO ix=1,ixe
        pionpd = recomb(ix,kz) + recomb12(ix,kz) - crg(kz) !small pos ions
        pionnd = recomb(ix,kz) + recomb21(ix,kz) - crg(kz) !small neg ions
        
        plionpd = recomb21(ix,kz) !large pos ions
        plionnd = recomb12(ix,kz) !large neg ions
        
!        IF ( jy .eq. ny/2 .and. ix.eq.nx/2 ) THEN
!         write(iunit,*) 'kz,pionpd,recomb,crg',kz,pionpd,recomb(ix,kz),
!     :     crg(kz)
!        ENDIF
        
        DO ia=lqb,lqe
          IF( ndebug .ge. 1 ) THEN
          IF (  attach(ix,kz,ia,1) .gt. 0.0 ) THEN
!$OMP CRITICAL
           write(iunit,*) 'ix,jy,kz,ia,1,attach= ',ix,jy,kz,ia,       &
     &        attach(ix,kz,ia,1)
!$OMP END CRITICAL
          ENDIF
          ENDIF
          pionpd = pionpd + attach(ix,kz,ia,1)
          IF ( ndebug .ge. 1 ) THEN
          IF (  attach(ix,kz,ia,2) .gt. 0.0 ) THEN
!$OMP CRITICAL
           write(iunit,*) 'ix,jy,kz,ia,2,attach= ',ix,jy,kz,ia,       &
     &        attach(ix,kz,ia,2)
!$OMP END CRITICAL
          ENDIF
          ENDIF
          pionnd = pionnd + attach(ix,kz,ia,2)
        ENDDO ! ia 
        
!
! make sure small ion numbers dont go below zero
!
        IF ( pionpd*dt .gt.        &
     &       Max(0.0,an(ix,jy,kz,lscpi)  )) THEN
          IF ( ndebug .ge. 1 ) THEN
!$OMP CRITICAL
          write(iunit,*) 'pionpd too big! must reduce from ',pionpd,       &
     &     ' to ', (an(ix,jy,kz,lscpi))/dt,       &
     &    ' at ix,jy,kz ',ix,jy,kz,       &
     &   an(ix,jy,kz,lscpi),       &
     &   an(ix,jy,kz,lscni)
          
          write(iunit,*) 'crg,recomb = ',crg(kz),recomb(ix,kz),       &
     &     an(ix,jy,kz,lscpi),an(ix,jy,kz,lscni),pionpd*dt

           DO ia=lqb,lqe
             write(iunit,*) 'ia,attach(1,2)',ia,attach(ix,kz,ia,1),       &
     &        attach(ix,kz,ia,2),an(ix,jy,kz,ia),       &
     &        an(ix,jy,kz,ia+lscb-lqb),       &
!     &         diff(ix,kz,ia,1),diff(ix,kz,ia,2),       &
!     &         cond(ix,kz,ia,1),cond(ix,kz,ia,2),elec(iez)%flt3d(ix,jy,kz),       &
     &       conc(ix,kz,ia)
           ENDDO
!$OMP END CRITICAL
           ENDIF

           attot = 0.0d0
           DO ia=lqb,lqe
             attot = attot + attach(ix,kz,ia,1) 
           ENDDO

          x = (an(ix,jy,kz,lscpi)/dt + crg(kz))/       &
     &          (attot + recomb(ix,kz) + recomb12(ix,kz))
           
           IF ( x .gt. 1.0 ) THEN
           IF ( x - 1.0 .lt. 1.0e-4 ) THEN  ! eliminate rare round-off/truncation problem
             x = 1.0
           ELSE
!$OMP CRITICAL
             write(iunit,*) 'WARNING! (pos ion) x>1.0 x = ',x
             write(iunit,*) 'crg,recomb = ',crg(kz),recomb(ix,kz),       &
     &             an(ix,jy,kz,lscpi),an(ix,jy,kz,lscni),pionnd*dt

           DO ia=lqb,lqe
             write(iunit,*) 'ia,attach(1,2)',ia,attach(ix,kz,ia,1),       &
     &        attach(ix,kz,ia,2),an(ix,jy,kz,ia),       &
     &        an(ix,jy,kz,ia+lscb-lqb),       &
!     &         diff(ix,kz,ia,1),diff(ix,kz,ia,2),       &
!     &         cond(ix,kz,ia,1),cond(ix,kz,ia,2),elec(iez)%flt3d(ix,jy,kz),       &
     &       conc(ix,kz,ia)
            ENDDO
             istop = 1 ! STOP
!$OMP END CRITICAL
           ENDIF
           ENDIF
           

           DO ia=lqb,lqe
             attach(ix,kz,ia,1) = x*attach(ix,kz,ia,1)
           ENDDO
           
           an(ix,jy,kz,lscpi) = 0.0 
           
         ELSE
           an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) - pionpd*dt
         ENDIF !pionpd

        IF ( pionnd*dt .gt.        &
     &       Max(0.0,an(ix,jy,kz,lscni) )) THEN
          IF ( ndebug .ge. 1 ) THEN
!$OMP CRITICAL
          write(iunit,*) 'pionnd too big! must reduce from ',pionnd,       &
     &     ' to ', (an(ix,jy,kz,lscni))/dt,       &
     &    ' at ix,jy,kz ',ix,jy,kz
          write(iunit,*) 'crg,recomb = ',crg(kz),recomb(ix,kz),       &
     &     an(ix,jy,kz,lscpi),        &
     &      an(ix,jy,kz,lscni),pionnd*dt
!           DO ia=lscwi,lscwi+nhab
           DO ia=lqb,lqe
             write(iunit,*) 'ia,attach(1,2)',ia,attach(ix,kz,ia,1),       &
     &        attach(ix,kz,ia,2),an(ix,jy,kz,ia),       &
     &        an(ix,jy,kz,ia+lscb-lqb),       &
!     &         diff(ix,kz,ia,1),diff(ix,kz,ia,2),       &
!     &         cond(ix,kz,ia,1),cond(ix,kz,ia,2),elec(iez)%flt3d(ix,jy,kz),       &
     &       conc(ix,kz,ia)
           ENDDO
!$OMP END CRITICAL
           
           ENDIF
          
          
           attot = 0.0d0
           DO ia=lqb,lqe
             attot = attot + attach(ix,kz,ia,2) 
           ENDDO

          x = (an(ix,jy,kz,lscni)/dt + crg(kz))/       &
     &          (attot + recomb(ix,kz) + recomb21(ix,kz))
           
           IF ( x .gt. 1.0 ) THEN
           IF ( x - 1.0 .lt. 1.0e-4 ) THEN  ! eliminate rare round-off/truncation problem
             x = 1.0
           ELSE
!$OMP CRITICAL
             write(iunit,*) 'WARNING! (neg ion) x>1.0 x = ',x
          write(iunit,*) 'crg,recomb = ',crg(kz),recomb(ix,kz),       &
     &     an(ix,jy,kz,lscpi),an(ix,jy,kz,lscni),pionnd*dt

           DO ia=lqb,lqe
             write(iunit,*) 'ia,attach(1,2)',ia,attach(ix,kz,ia,1),       &
     &        attach(ix,kz,ia,2),an(ix,jy,kz,ia),       &
     &        an(ix,jy,kz,ia+lscb-lqb),       &
!     &         diff(ix,kz,ia,1),diff(ix,kz,ia,2),       &
!     &         cond(ix,kz,ia,1),cond(ix,kz,ia,2),elec(iez)%flt3d(ix,jy,kz),       &
     &       conc(ix,kz,ia)
           ENDDO
            istop = 1 ! STOP
!$OMP END CRITICAL
           ENDIF
           ENDIF


           DO ia=lqb,lqe
             attach(ix,kz,ia,2) = x*attach(ix,kz,ia,2)
           ENDDO
           
           an(ix,jy,kz,lscni) = 0.0 
           
         ELSE
           an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - pionnd*dt
         ENDIF !pionnd

!
! make sure large ion numbers dont go below zero
!
         IF ( largeion ) THEN
         IF (  plionpd*dt .gt.        &
     &       Max(0.0,an(ix,jy,kz,lscpli)) ) THEN
          IF ( ndebug .ge. 1 ) THEN
!$OMP CRITICAL
          write(iunit,*) 'plionpd too big! must reduce from ',plionpd,       &
     &     ' to ', (an(ix,jy,kz,lscpli))/dt,       &
     &    ' at ix,jy,kz ',ix,jy,kz,       &
     &   an(ix,jy,kz,lscpli),       &
     &   an(ix,jy,kz,lscnli)
          
          write(iunit,*) 'crg,recomb21 = ',crg(kz),recomb21(ix,kz),       &
     &     an(ix,jy,kz,lscpli),an(ix,jy,kz,lscnli),plionpd*dt
!$OMP END CRITICAL
           ENDIF

           x = (an(ix,jy,kz,lscpli)/dt) / (recomb21(ix,kz))
           
           IF ( x .gt. 1.0 ) THEN
           IF ( x - 1.0 .lt. 1.0e-4 ) THEN  ! eliminate rare round-off/truncation problem
             x = 1.0
           ELSE
!$OMP CRITICAL
             write(iunit,*) 'WARNING! (pos lg ion) x>1.0 x = ',x
             write(iunit,*) 'crg,recomb21 = ',crg(kz),recomb21(ix,kz),       &
     &             an(ix,jy,kz,lscpli),an(ix,jy,kz,lscnli),plionpd*dt

             istop = 1 ! STOP
!$OMP END CRITICAL
         ENDIF
           ENDIF
           
           an(ix,jy,kz,lscpli) = 0.0 
           
         ELSE
           an(ix,jy,kz,lscpli) = an(ix,jy,kz,lscpli) - plionpd*dt
         ENDIF !plionpd
         ENDIF !largeion

        IF ( largeion ) THEN
        IF ( plionnd*dt .gt.        &
     &       Max(0.0,an(ix,jy,kz,lscnli) )) THEN
          IF ( ndebug .ge. 1 ) THEN
!$OMP CRITICAL
          write(iunit,*) 'plionnd too big! must reduce from ',plionnd,       &
     &     ' to ', (an(ix,jy,kz,lscnli))/dt,       &
     &    ' at ix,jy,kz ',ix,jy,kz
          write(iunit,*) 'crg,recomb12 = ',crg(kz),recomb12(ix,kz),       &
     &     an(ix,jy,kz,lscpi),        &
     &      an(ix,jy,kz,lscni),pionnd*dt
!$OMP END CRITICAL
           
          ENDIF
          
          x = (an(ix,jy,kz,lscnli)/dt) / (recomb12(ix,kz))
           
           IF ( x .gt. 1.0 ) THEN
           IF ( x - 1.0 .lt. 1.0e-4 ) THEN  ! eliminate rare round-off/truncation problem
             x = 1.0
           ELSE
!$OMP CRITICAL
             write(iunit,*) 'WARNING! (neg lg ion) x>1.0 x = ',x
          write(iunit,*) 'crg,recomb12 = ',crg(kz),recomb12(ix,kz),       &
     &     an(ix,jy,kz,lscpli),an(ix,jy,kz,lscnli),pionnd*dt

            istop = 1 ! STOP
!$OMP END CRITICAL
           ENDIF
           ENDIF

           an(ix,jy,kz,lscnli) = 0.0 
           
         ELSE
           an(ix,jy,kz,lscnli) = an(ix,jy,kz,lscnli) - plionnd*dt
         ENDIF !plionnd
         ENDIF !largeion

         

           DO ia=lqb,lqe
             an(ix,jy,kz,ia-lqb+lscb) = an(ix,jy,kz,ia-lqb+lscb) +        &
     &            ec*dt*(attach(ix,kz,ia,1)-attach(ix,kz,ia,2))
            IF ( (attach(ix,kz,ia,1).gt.0.0 .or.        &
     &           attach(ix,kz,ia,2).gt.0.0 ) .and.       &
     &            ndebug .ge. 1     ) THEN
!$OMP CRITICAL
             write(iunit,*) 'attach2: ix,jy,kz,ia ',ix,jy,kz,ia,       &
     &        attach(ix,kz,ia,1),       &
     &    attach(ix,kz,ia,2),attach(ix,kz,ia,1)-attach(ix,kz,ia,2),       &
     &     an(ix,jy,kz,ia-lqb+lscb),       &
     &     ec*dt*(attach(ix,kz,ia,1)-attach(ix,kz,ia,2))
!$OMP END CRITICAL
            ENDIF
          ENDDO
         
      ENDDO
      ENDDO
       
      
      ENDDO  ! jy loop

      CALL cld_cpu('ION-ATTACH')
      
      IF ( istop .eq. 1 ) STOP

!
! if multiple steps, put an(ions) into ad(ions)
!
!      IF ( nionstep .gt. 1 .and. n .lt. nionstep ) THEN
!         call swap(nx,ny,nz,2,nor,
!     :     an(-nor+1,-nor,-nor,lscpi),an(-nor,-nor,-nor,lscpi))
!
! recompute electric field?  Taken care of after drift, if necessary.
!
!      ENDIF


!      IF (  n .lt. nionstep ) THEN
      IF (  .false. ) THEN
      
      CALL cld_cpu('ION-EFIELD')

      t0(:,:,:) = 0.0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,ia)
      do kz = 1, kze
      DO ia = lscb,lsceq
      do jy = 1, jye
      do ix = 1, ixe
        t0(ix,jy,kz) = t0(ix,jy,kz) + an(ix,jy,kz,ia)
      end do
      end do
      ENDDO
      end do
 
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       IF ( largeion ) THEN
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))         &
     &  + ec*(an(ix,jy,kz,lscpli)- an(ix,jy,kz,lscnli))  
       ELSE
       t0(ix,jy,kz) = t0(ix,jy,kz)        &
     &  + ec*(an(ix,jy,kz,lscpi)- an(ix,jy,kz,lscni))  
       ENDIF
      
      end do
      end do
      end do

      IF ( ndebug .ge. 1 ) write(iunit,*) 'call mud'

#if defined (MPI) && defined (BOXMG) 

         write(iunit,*) 'ionstep2: call putf'
         
!         IF ( my_rank == 0 ) THEN
!           DO k = 1,NLzdg
!             write(iunit,*) 'k,rhofair = ',k,rhofair(k)
!           ENDDO
!         ENDIF
      CALL cld_cpu('BOXMG')
         
         CALL PUTF( SOdg, QFdg, QFdg, Qdg, Qdg, t0, elec(ipot)%flt3d, nx,ny,nz,nor,            &
     &               nbw,nbe,nbs,nbn,       &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz, 0       &
     &              )
!         write(iunit,*) 'done putf'


         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
        
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dpt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg'

         CALL commas_SymStd_SOLVE_boxmgdg(                                          &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg,       &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                          &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,        &
     &             SOdg, NSOdg,                                                     &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                          &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                            &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,                      &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                                 &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                                 &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                                 &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                                 &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                         &
     &             pMSGdg, pMSGSOdg,                                                &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,                &
     &             number_of_processes, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg,       &
     &             my_comm                                               &
     &             )

!         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dpt2 = MPI_Wtime()
!
         IF (my_rank .eq. 0) THEN
            WRITE(iunit,*) 'Total Time = ', dpt2 - dpt1
         ENDIF

!      IF ( my_rank .eq. 0 ) THEN
!          CALL PRINTQ( SOdg, QFdg, Qdg,                      
!     :               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   
!     :               iGsdg, jGsdg, kGsdg, dx, dy, dlz, 0    
!     :              )
!      ENDIF

!
         
         CALL PUTPHI( Qdg, elec, nx,ny,nz,nor,            &
     &               gxt,gyt,gzt,       &
     &               nbw,nbe,nbs,nbn,iunit,       &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,          &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz       &
     &              )

      CALL cld_cpu('BOXMG')

#else

#ifndef BOXMG
      IF ( number_of_processes .eq. 1 ) THEN
#endif
      CALL mud(nx,ny,nz,nor,       &
     &   dx,dy,dz,gxt,gyt,gzt,       &
     &   nnxs(1),nnys(1),nnzs(1),       &
     &   iixps(1),jjyqs(1),kkzrs(1),iiexs(1),jjeys(1),kkezs(1),       &
     &   nbw,nbe,nbs,nbn,       &
     &   llworks(1),istretch,       &
     &   t0,ibg,       &
     &   elec,iestag,iunit,.true., bcx, bcy)
!
!
#endif

      IF ( ndebug .ge. 1 ) write(iunit,*) 'call efield'



      IF ( iestag .eq. 0 ) THEN
      call efield       &
     &  (nx,ny,nz       &
     &  ,id1,jd1,kd1,istag,jstag,kstag       &
     &  ,dx,dy,dz,elec       &
     &  ,gxt(1,3),gyt(1,3),gzt(1,3)       &
     &  ,iunit)    
       
         IF ( ndebug .ge. 1 ) write(iunit,*) 'done efield'
         


!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,ia)
      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
      
!       IF ( .true. ) THEN
       IF ( kz .gt. 1 ) THEN
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-elec(ipot)%flt3d(ix,jy,kz-1))*gzt(kz,3)
       ELSE
       ezw(ix,jy,kz) =        &
     & -(elec(ipot)%flt3d(ix,jy,kz)-0.0)*gzt(kz,3)
       
       ENDIF
       
!       ELSE
!       IF ( kz .eq. 1 ) THEN
!       ezw(ix,jy,kz) = ezfairw(kz)
!       ELSE
!       ezw(ix,jy,kz) = 0.5*(elec(ix,jy,kz,iez)+elec(ix,jy,kz-1,iez)) -
!     :     0.5*(elec(ix,jy,kz,iezb)+elec(ix,jy,kz-1,iezb)) +
!     :     ezfairw(kz)
!       ENDIF
!          elec(ix,jy,kz,iez) = elec(ix,jy,kz,iez) - elec(ix,jy,kz,iezb)
!     :       + ezfair(kz)
!       ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE

       DO kz = 1,kze
         km1 = Max( 1   ,kz-1 )
         kp1 = Min( nz-1,kz   )
         DO jy = 1,jye
         DO ix = 1,ixe
         ezw(ix,jy,kz) = 0.5*(elec(iez)%flt3d(ix,jy,km1) + elec(iez)%flt3d(ix,jy,kp1))
         IF ( Sign(1.0,ezw(ix,jy,kz)) .ne. Sign(1.0,elec(iez)%flt3d(ix,jy,kz) ) ) THEN
           ezw(ix,jy,kz) = 0.5*elec(iez)%flt3d(ix,jy,kz)
         ENDIF
         ENDDO
         ENDDO
       ENDDO

      
      ENDIF ! iestag

#ifndef BOXMG
      ENDIF
#endif
      CALL cld_cpu('ION-EFIELD')

      ENDIF ! false hack

      ENDDO ! nionstep

      CALL cld_cpu('ION-MISC')

      sciont0 = 0.0d0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      sciont0s = 0.0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv),   &
!$OMP  REDUCTION (+ : sciont0,sciont0s,scionp0,scionn0)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)

       IF ( largeion ) THEN
       
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) *dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli)) *dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      
      ELSE
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv
      
       sciont0s = sciont0s       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) *dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont0
       mpitotindp(2)  = sciont0s
       mpitotindp(3)  = scionp0
       mpitotindp(4)  = scionn0
       mpitotindp(5)  = chaffiontot
       
      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 5, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont0 = mpitotoutdp(1)
       sciont0s = mpitotoutdp(2)
       scionp0 = mpitotoutdp(3)
       scionn0 = mpitotoutdp(4)
       chaffiontot = mpitotoutdp(5)
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postrecomb ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0

      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'postrecomb single pres. ion net charge (C):',       &
     &      sciont0s
      
       IF ( lnchaff .gt. 1 .or. chaffconc > 0. ) THEN
      write(iunit,'(a,1(2x,1pe15.8))')        &
     &     'chaff net ion charge production (C):',       &
     &      chaffiontot*ec
       ENDIF
      
      ENDIF

      DO kz=1,kze
      DO jy=1,jye
      DO ix=1,ixe
        an(ix,jy,kz,lscpi) = Max(0.0,an(ix,jy,kz,lscpi))
        an(ix,jy,kz,lscni) = Max(0.0,an(ix,jy,kz,lscni))
      ENDDO
      ENDDO
      ENDDO
      

      sciont0 = 0.0d0
      scionp0 = 0.0d0
      scionn0 = 0.0d0
      
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix,dv),   &
!$OMP  REDUCTION (+ : sciont0,scionp0,scionn0)
      do kz = 1, kze
      do jy = 1, jye
      do ix = 1, ixe
       dv = dxx(ix)*dyy(jy)*dzz(kz)

       IF ( largeion ) THEN
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli) - an(ix,jy,kz,lscnli))*dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv       &
     &  + ec*(an(ix,jy,kz,lscpli))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv       &
     &  + ec*( - an(ix,jy,kz,lscnli))*dv
      ELSE
       sciont0 = sciont0       &
     &  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*dv
      
       scionp0 = scionp0       &
     &  + ec*(an(ix,jy,kz,lscpi))*dv

       scionn0 = scionn0        &
     &  + ec*( - an(ix,jy,kz,lscni))*dv
      
      ENDIF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = sciont0
       mpitotindp(2)  = sciont0s
       mpitotindp(3)  = scionp0
       mpitotindp(4)  = scionn0

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       sciont0 = mpitotoutdp(1)
       sciont0s = mpitotoutdp(2)
       scionp0 = mpitotoutdp(3)
       scionn0 = mpitotoutdp(4)
      ENDIF
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe15.8))')        &
     &     'postzero ion pos/neg/net charge (C):',       &
     &      scionp0,scionn0,sciont0
      ENDIF

      CALL cld_cpu('ION-MISC')

      deallocate( attach )
      deallocate( conc )
      deallocate( diff )
      deallocate( cond )

      RETURN
      END
      
      
!
! ####################################################################
! ####################################################################
