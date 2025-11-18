! #include "sam.def.h"
!
!  subroutine to calculate fall speeds of hydrometeors
!
      subroutine ziegfall(nx,ny,nz,nor,norz,na,dtp,dz,jgs,ngs1, &
     &  xvt,                                                    &
     &  an,dn,ipconc0,t0,t7,cwccn0,cwmasn,cwmasx,cimn0,cimx0,       &
     &  rwmasn,rwmasx,cwradn,                                   &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx0,xvmn,xvmx,cwnccn,  &
     &  itype1x,itype2x,infdo)

       USE INDEX_MODULE, only: lt,lc,lr,li,lis,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lhw,lfw,lhlw, &
                               lnc,lnr,lni,lnis,lns,lnh,lnf,lnhl,cinu,dmuh,dnu,dmu,xnu,xmu,dmuhl,rnu,cnu,snu, &
                               lss,lsat,lsati,xvcmx,xvcmn,xvrmn,xvrmx,lccn,rnumin,rnumax, &
                               lqmx,nxtra,lvi,lvs,lvh,lvf,lvhl,lzi,lzs,lzr,lzh,lzf,lzhl,lsw, &
                               alphar,alphamin,alphamax, lnhf,lnhlf
      USE COMMASMPI_MODULE
      USE CPUTIME_MODULE
      USE MICRO_MODULE

! 12.16.2005: .F version use in transitional SWM model
!
! 10.10.2003: Added cimn and cimx to setting for cci and cip.
!
! TO DO LIST:
!
! need to set up values for:
!     :  cipdia,cidia,cwdia,cwmas,vtwbar,
!     :  rho0,temcg,cip,cci
!
! and need to put fallspeed values in cwvt etc.
!
      
      implicit none

      integer ng1
      parameter(ng1 = 1)
      
      integer nx,ny,nz,nor,norz,ngt,jgs,na
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
!      real vt(nx,nz)
      real dtp,dz,dtz1
      
      real ccwmx0
      real cwnccn(nz)
      real, parameter :: pi1 = 3.14159265 ! 4.0*atan(1.0)
      real, parameter :: piinv = 1.0/pi1
      real, parameter :: cwc1 = 6.0/(pi1*1000.)
!      real qtmp1(nx,nz),qtmp2(nx,nz)
!      real cmax
!      real xfall(nx,ny)
      
      integer ndebugzf
      parameter (ndebugzf = 0)
      integer idx

      integer ix,jy,kz,i,j,k,item,il
      integer itype1x,itype2x,infdo
!
!  include file for mix ratio and charge indices
!
!      include 'swm.index.zieg.h'
!
      real xvt(nx,nz+1,3,lc:lhab) ! 1=mass-weighted, 2=number-weighted

      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      
      real xdnmx(lc:lhab), xdnmn(lc:lhab)

!
!   drag coefficients
!
      real cdx(lc:lhab)
!
! Fixed intercept values for single moment scheme
!
      real cno(lc:lhab)
      
      real cwccn0,cwmasn,cwmasx,cimn0,cimx0,cwradn
!      real xvcmn, xvcmx  ! min, max droplet volumes
!      real xvrmn, xvrmx  ! min, max rain volumes
!      real xvsmn, xvsmx  ! min, max snow volumes
!      real xvfmn, xvfmx  ! min, max frozen drop volumes
!      real xvgmn, xvgmx  ! min, max graupel volumes
!      real xvhmn, xvhmx  ! min, max hail volumes
!      real xvhlmn, xvhlmx  ! min, max lg hail volumes
      real rwmasn,rwmasx
!      PARAMETER(XVCMN=4.188E-12,XVCMX=6.54E-8)   ! CGS
!      PARAMETER(XVRMN=2.8866E-7,XVRMX=4.1887E-3) ! CGS
!      parameter( xvcmn=4.188e-18,  xvcmx=6.54e-14 )    ! mks
!c      parameter( xvrmn=2.8866e-13, xvrmx=4.1887e-9 )  ! mks
!      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(2.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvsmn=2.8866e-13, xvsmx=0.523599*(3.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvfmn=2.8866e-13, xvfmx=0.523599*(3.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvgmn=2.8866e-13, xvgmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhmn=0.523599*(1.e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhlmn=0.523599*(10.e-3)**3, 
!           xvhlmx=0.523599*(100.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      integer ngs1
!      real cinccn(ngs1)
      real cwc0

      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag,istag
      parameter (kstag=1, istag=1)

      integer ngs,ngscnt,mgs,ipconc0
      parameter ( ngs=500 )
      integer igs(ngs),kgs(ngs)
      
      real :: qx(ngs,lv:lhab)
      real :: qxw(ngs,ls:lhab)
      real :: cx(ngs,lc:lhab)
      real :: chxf(ngs,lh:lhab)
      real :: xv(ngs,lc:lhab)
      real :: vtxbar(ngs,lc:lhab,3)
      real :: xmas(ngs,lc:lhab)
      real :: xdn(ngs,lc:lhab)
      real :: cdxgs(ngs,lc:lhab)
      real :: xdia(ngs,lc:lhab,3)
      real :: vx(ngs,li:lhab)
      real :: alpha(ngs,lc:lhab)
      real :: alphan(ngs,lc:lhab)
      real :: zx(ngs,lr:lhab)

!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)
      real :: axx(ngs,lh:lhab),bxx(ngs,lh:lhab)

      real rho0(ngs),temcg(ngs)
!     :  qcw(ngs),qrw(ngs),qci(ngs),qip(ngs),qir(ngs),
!     :  qsw(ngs),qgl(ngs),qgm(ngs),
!     :  qgh(ngs),qhw(ngs),qhl(ngs),qfw(ngs),qwv(ngs),
!     :  qvs(ngs)

!      real ccw(ngs),cci(ngs),crw(ngs),csw(ngs),chw(ngs) 
!      real cgl(ngs),cgh(ngs),cfw(ngs)
!      real cgm(ngs),chl(ngs),cir(ngs)
!      real cip(ngs)

      real temg(ngs)
      
      real rhovt(ngs)
!     : , cirdn(ngs),cwdn(ngs),rwdn(ngs),
!     : swdn(ngs), gldn(ngs),gmdn(ngs),ghdn(ngs),
!     : fwdn(ngs), hwdn(ngs), hldn(ngs), cipdn(ngs),cidn(ngs)

!      real cirdia(ngs), rwdia(ngs), swdia(ngs),gldia(ngs),
!     : gmdia(ngs), ghdia(ngs), hwdia(ngs),hldia(ngs),fwdia(ngs)
!      real cipdia(ngs),cidia(ngs)
      
!      real rwrad(ngs)

!      real vtirbar(ngs),vtrbar(ngs),vtsbar(ngs),vtglbar(ngs),
!     :  vtgmbar(ngs), vtghbar(ngs),vtfbar(ngs),vthbar(ngs),
!     :  vthlbar(ngs) ,vtipbar(ngs),vtibar(ngs)
     
!      real vtrnbar(ngs)

      
!      real ciplen(ngs),cilen(ngs)
      
      real cwnc(ngs),cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
      real cimasn,cimasx,cnina(ngs),cimas(ngs)
!      real xvr(ngs),rwmas(ngs)
!      real xvs(ngs),xvgl(ngs),xvgm(ngs),xvgh(ngs),xvf(ngs)
!      real xvh(ngs),xvhl(ngs)
!      real swmas(ngs),glmas(ngs),gmmas(ngs),ghmas(ngs)
!      real fwmas(ngs),hwmas(ngs),hlmas(ngs)
!      real mwfac
!      parameter ( mwfac = 6.0**(1./3.) ) ! factor for mass-weighted rain volume diameter

!

      real poo, cp608, cp, cv
      real dnz00, rho00, cs, ds
      real pi, pii, pid4, qccrit, qscrit
!
! intercepts
!
!
!  density maximums and minimums
!
      real rwdnmx, cwdnmx, cidnmx, xidnmx
      real swdnmx, gldnmx, gmdnmx, ghdnmx, fwdnmx, hwdnmx, hldnmx
!
      real rwdnmn, cwdnmn, xidnmn, cidnmn
      real swdnmn, gldnmn, gmdnmn, ghdnmn, fwdnmn, hwdnmn
!
!  constants
!
      real c1f3
!
!  general constants for microphysics
!
       real tfr, advisc0

! 
! Miscellaneous
!
!      character*80 infile, outfile
!      integer ifile
      integer ihabdo
      parameter(ihabdo = 0)
      integer ierr
      
      logical flag
      logical ldovol, ldoliq
      integer lvol(lc:lhab)
      integer ln(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
      logical lrescalelow(lc:lhab)
      
      real chw, qr, z, rd, alp, z1, g1, vr, nrx, tmp, frac, frach
      real x,y,del,cwchtmp,cwcrtmp
      
      real vtmax
      real xvbarmax

      real, parameter ::  c1r=19.0, c2r=0.6, c3r=1.8, c4r=17.0   ! rain
      real, parameter ::  c1h=5.5, c2h=0.7, c3h=4.5, c4h=8.5   ! Graupel
      real, parameter ::  c1hl=3.7, c2hl=0.3, c3hl=9.0, c4hl=6.5, c5hl=1.0, c6hl=6.5 ! Hail

! inline functions for Newton method
       real :: galpha, dgalpha
       real :: a_in
       logical, parameter :: newton = .false.

      galpha(a_in) = ((4. + a_in)*(5. + a_in)*(6. + a_in))/((1. + a_in)*(2. + a_in)*(3. + a_in))
      dgalpha(a_in) = (876. + 1260.*a_in + 621.*a_in**2 + 126.*a_in**3 + 9.*a_in**4)/            &
     &  (36. + 132.*a_in + 193.*a_in**2 + 144.*a_in**3 + 58.*a_in**4 + 12.*a_in**5 + a_in**6)


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.


#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: ENTERED SUBROUTINE, my_rank=",my_rank
#else
      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL: ENTERED SUBROUTINE"
#endif

! #####################################################################
! BEGIN EXECUTABLE
! #####################################################################
!

!  constants
!
      ldovol = .false.
      lvol(:) = 0
      IF ( lvi .gt. 1 ) lvol(li) = lvi
      IF ( lvs .gt. 1 ) lvol(ls) = lvs
      IF ( lvh .gt. 1 ) lvol(lh) = lvh
      IF ( lhl .gt. 1 .and. lvhl .gt. 1 ) lvol(lhl) = lvhl
      IF ( lf .gt. 1 .and. lvf .gt. 1 ) lvol(lf) = lvf
      
      
! #ifdef Z3MOM
      lrescalelow(:) = rescale_low_alpha
      lrescalelow(lr) = rescale_low_alphar .and. rescale_low_alpha
      lrescalelow(lh) = rescale_low_alphah .and. rescale_low_alpha
      IF ( lf > 1 ) lrescalelow(lf) = rescale_low_alphah .and. rescale_low_alpha
      IF ( lhl > 1 ) lrescalelow(lhl) = rescale_low_alphahl .and. rescale_low_alpha
!      write(91,*) 'my_rank, lrescalelow: ',my_rank,lrescalelow(lr),lrescalelow(lh)
! #endif
      
      IF ( li .gt. 1 ) THEN
      DO il = li,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO
      ENDIF

      ln(:) = 0
      ln(lc) = lnc
      ln(lr) = lnr
      IF ( li > 0 ) THEN
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      ENDIF
      IF ( lhl .gt. 1 ) ln(lhl) = lnhl
      IF ( lf .gt. 1 ) ln(lf) = lnf
      IF ( lis > 1 ) ln(lis) = lnis

      lz(:) = 0
      lz(lr) = lzr
      IF ( li > 0 ) THEN
      lz(li) = lzi
      lz(ls) = lzs
      lz(lh) = lzh
      ENDIF
      IF ( lhl .gt. 1 .and. lzhl > 1 ) lz(lhl) = lzhl
      IF ( lf .gt. 1 .and. lzf > 1 ) lz(lf) = lzf

      lliq(:) = 0
      IF ( lsw .gt. 1 ) lliq(ls) = lsw
      IF ( lhw .gt. 1 ) lliq(lh) = lhw
      IF ( lhl .gt. 1 .and. lhlw .gt. 1 ) lliq(lhl) = lhlw
      IF ( lf .gt. 1 ) lliq(lf) = lfw

      ldoliq = .false.
      IF ( ls .gt. 1 ) THEN
      DO il = ls,lhab
        ldoliq = ldoliq .or. ( lliq(il) .gt. 1 )
      ENDDO
      ENDIF
      
      poo = 1.0e+05
      cp608 = 0.608
      cp = 1004.0
      cv = 717.0
      dnz00 = 1.225
      rho00 = 1.225
      cs = 4.83607122
      ds = 0.25
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pi = 4.0*atan(1.0)
      pii = 1./pi
      pid4 = pi/4.0 
      qccrit = 2.0e-03
      qscrit = 6.0e-04
      cwc0 = pii
      advisc0 = 1.832e-05
!
!  constants
!
      c1f3 = 1.0/3.0
!
!  general constants for microphysics
!
      tfr = 273.15
!
!  ci constants in mks units
!
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
!
!  Set terminal velocities...
!    also set drag coefficients
!

#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: start loop"
#endif
      jy = jgs
!      DO kz = 1,nz-1
!        DO ix = 1,nx-1
      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      if (ixend .eq. nxend) nxz = (ixend-ixbeg+1-istag)*nz
      numgs = nxz/ngs + 1

!      ixb = 1
!      ixe = itile
!      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

!#ifdef MPI
!      numgs = ixe
!      do ix = 1,ixe
!#else
      do inumgs = 1,numgs
!#endif
       ngscnt = 0

#ifdef MPI

      kzb = nzmpb
      kze = ktile+1
      if (kzend .eq. nzend) kze = kzend-kzbeg-kstag

      ixb = nxmpb
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

       do kz = kzb,kze
        do ix = nxmpb,ixe
#else
       do kz = nzmpb,nz-kstag-1
        do ix = nxmpb,nx-istag
#endif
        flag = .false.

        DO il = lc,lhab
          flag =  flag .or. ( an(ix,jy,kz,il)  .gt. qxmin(il) ) 
        ENDDO

        if ( flag ) then
! load temp quantities
!        ngscnt = 1
!        mgs = 1
!        igs(mgs) = ix
!        kgs(mgs) = kz

        ngscnt = ngscnt + 1
        igs(ngscnt) = ix
        kgs(ngscnt) = kz
        if ( ngscnt .eq. ngs ) goto 1100
        end if
!#ifndef MPI
        end do !!ix
!#endif
        nxmpb = 1
       end do !! kz

!      if ( jy .eq. (ny-jstag) ) iend = 1

 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998
!
!  set temporaries for microphysics variables
!

#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: after allocate"
#endif

!
!  Reconstruct various quantities 
!
      do mgs = 1,ngscnt

       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhovt(mgs) = Sqrt(rho00/rho0(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
       temcg(mgs) = temg(mgs) - tfr

       fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &   (temg(mgs)/296.0)**(1.5)
!
      end do
!
      IF ( ipconc .eq. 0 ) THEN
      do mgs = 1,ngscnt
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      end do
      ENDIF


      vtxbar(:,:,:) = 0.0
      
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

!
!  set concentrations
!
      cx(:,:) = 0.0
      IF ( lnhf > 1 .or. lnhlf > 1  ) chxf(:,:) = 0.0
      
      if ( ipconc .ge. 1 .and. li .gt. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
        IF ( lis > 1 ) cx(mgs,lis) = Max(an(igs(mgs),jy,kgs(mgs),lnis), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 .and. lc .gt. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
!        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
!        IF ( na .ge. lccn ) THEN
!         ccn(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
!        ELSE
!         ccn(mgs) = 0.0
!        ENDIF
       end do
      end if
      if ( ipconc .ge. 3 .and. lr .gt. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
!          cx(mgs,lr) = 0.0
!        ELSEIF ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) THEN
!          qx(mgs,lr) = 0.0
        ELSE
          cx(mgs,lr) = Max( 0.0, cx(mgs,lr) )
        ENDIF
       end do
      end if
      if ( ipconc .ge. 4  .and. ls .gt. 1) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
!          cx(mgs,ls) = 0.0
!        ELSEIF ( cx(mgs,ls) .eq. 0.0 .and. qx(mgs,ls) .lt. 3.0*qxmin(ls) ) THEN
!          qx(mgs,ls) = 0.0
        ELSE
          cx(mgs,ls) = Max( 0.0, cx(mgs,ls) )
        ENDIF
       end do
      end if

      if ( ipconc .ge. 5  .and. lh .gt. 1) then
       do mgs = 1,ngscnt

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
!          cx(mgs,lh) = 0.0
!        ELSEIF ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) THEN
!          qx(mgs,lh) = 0.0
        ELSE
          cx(mgs,lh) = Max( 0.0, cx(mgs,lh) )
        ENDIF

        IF ( lnhf > 1 ) THEN
           chxf(mgs,lh) = Min(cx(mgs,lh), Max(an(igs(mgs),jy,kgs(mgs),lnhf), 0.0))
        ENDIF

       end do
      ENDIF

      if ( ipconc .ge. 5  .and. lf .gt. 1) then
       do mgs = 1,ngscnt

        cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)
        IF ( qx(mgs,lf) .le. qxmin(lf) ) THEN
!          cx(mgs,lf) = 0.0
!        ELSEIF ( cx(mgs,lf) .eq. 0.0 .and. qx(mgs,lf) .lt. 3.0*qxmin(lf) ) THEN
!          qx(mgs,lf) = 0.0
        ELSE
          cx(mgs,lf) = Max( 0.0, cx(mgs,lf) )
        ENDIF


       end do
      ENDIF

      if ( ipconc .ge. 5  .and. lhl .gt. 1) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        IF ( qx(mgs,lhl) .le. qxmin(lhl) ) THEN
!          cx(mgs,lhl) = 0.0
!        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) THEN
!          qx(mgs,lhl) = 0.0
        ELSE
          cx(mgs,lhl) = Max( 0.0, cx(mgs,lhl) )
        ENDIF

        IF ( lnhlf > 1 ) THEN ! number of hail from frozen drops
           chxf(mgs,lhl) = Min(cx(mgs,lhl), Max(an(igs(mgs),jy,kgs(mgs),lnhlf), 0.0))
        ENDIF

       end do
      end if
       
      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
!        IF ( ls .gt. 1 .and. lvs .eq. 0 ) xdn(mgs,ls) = xdn0(ls)
!        IF ( lh .gt. 1 .and. lvh .eq. 0 ) xdn(mgs,lh) = xdn0(lh)
        IF ( li .gt. 1 )  xdn(mgs,li) = xdn0(li)
        IF ( lis > 1 )  xdn(mgs,lis) = xdn0(lis)
        IF ( ls .gt. 1 )  xdn(mgs,ls) = xdn0(ls)
        IF ( lh .gt. 1 )  xdn(mgs,lh) = xdn0(lh)
        IF ( lf .gt. 1 )  xdn(mgs,lf) = xdn0(lf)
        IF ( lhl .gt. 1 ) xdn(mgs,lhl) = xdn0(lhl)
      end do

!
! Set mean particle volume
!
      IF ( ldovol ) THEN
      
      vx(:,:) = 0.0
      
       DO il = li,lhab
        
        IF ( lvol(il) .ge. 1 ) THEN
        
          DO mgs = 1,ngscnt
            vx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lvol(il)), 0.0)
            IF ( vx(mgs,il) .gt. rho0(mgs)*qxmin(il)*1.e-3 .and. qx(mgs,il) .gt. qxmin(il) ) THEN
              xdn(mgs,il) = Min( xdnmx(il), Max( xdnmn(il), rho0(mgs)*qx(mgs,il)/vx(mgs,il) ) )
            ENDIF
          ENDDO
          
        ENDIF
      
       ENDDO
      
      ENDIF

      DO il = lg,lhab
      DO mgs = 1,ngscnt
        alpha(mgs,il) = dnu(il)
      ENDDO
      ENDDO
      
      IF ( imurain == 1 ) THEN
        alpha(:,lr) = alphar
      ELSEIF ( imurain == 3 ) THEN
        alpha(:,lr) = xnu(lr)
      ENDIF
      
      alphan(:,:) = alpha(:,:)

      
      IF ( ipconc == 5 .and. imydiagalpha > 0 ) THEN
        IF ( imydiagalpha == 1 .or. imydiagalpha == 2 ) THEN

            cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
            cwcrtmp = ((3. + dnu(lr))*(2. + dnu(lr))*(1.0 + dnu(lr)))**(-1./3.)

        DO mgs = 1,ngscnt
          IF ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > cxmin ) THEN
          ! MY 2005:
             xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))            ! 
             xdia(mgs,lr,3) = (xv(mgs,lr)*6.0*cwc1)**(1./3.) 
          !   alpha(mgs,lr) = Min(alphamax, c1r*tanh(c2r*(xdia(mgs,lr,3)*1000. - c3r)) + c4r)

            ! M&M-C 2010:
             tmp = 4. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*cwcrtmp*xdia(mgs,lr,3)

             alpha(mgs,lr) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
             alphan(mgs,lr) = alpha(mgs,lr)

          ENDIF
          IF ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) > cxmin ) THEN
            ! MY 2005:
             xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))            ! 
             xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
            ! alpha(mgs,lh) = Min(alphamax, c1h*tanh(c2h*(xdia(mgs,lh,3)*1000. - c3h)) + c4h)
             
            ! M&M-C 2010:
             tmp = 4. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*cwchtmp*xdia(mgs,lh,3)

             alpha(mgs,lh) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
             alphan(mgs,lh) = alpha(mgs,lh)
             !alp = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
        
!          write(0,*) 'zfall: alp = ',igs(mgs),kgs(mgs),alphan(mgs,lh),qx(mgs,lh)*1.e3,1.e3*xdia(mgs,lh,3)*cwchtmp

          ENDIF
!        alpha(:,lr) = 0. ! 10.
!        alpha(:,lh) = 0. ! 10.
          IF ( lhl > 0 ) THEN
          IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) > cxmin ) THEN
             xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))            ! 
             xdia(mgs,lhl,3) = (xv(mgs,lhl)*6.*piinv)**(1./3.)
             IF ( xdia(mgs,lhl,3) < 0.008 ) THEN
               alpha(mgs,lhl) = Min(alphamax, c1hl*tanh(c2hl*(xdia(mgs,lhl,3)*1000. - c3hl)) + c4hl)
             ELSE
               alpha(mgs,lhl) = Min(alphamax, c5hl*xdia(mgs,lhl,3)*1000. + c6hl)
             ENDIF
          ENDIF
          ENDIF
        ENDDO
        ELSEIF ( imydiagalpha == 3 ) THEN ! Milbrandt 2010
        
        ENDIF
      ENDIF

!
! Set 6th moments
!
      IF ( ipconc .ge. 6 .or. lzr > 1) THEN
      
      zx(:,:) = 0.0
      
       DO il = lr,lhab
        
        IF ( lz(il) .ge. 1 ) THEN
          DO mgs = 1,ngscnt
            zx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0)
          ENDDO
          
        
        ENDIF
      
       ENDDO
      
      ENDIF
       

!
! Set liquid water fractions
!
      IF ( ldoliq ) THEN
      
      DO il = ls,lhab
      IF ( lliq(il) .gt. 1 ) THEN
        DO mgs = 1,ngscnt
          qxw(mgs,il) = max(min(qx(mgs,il),an(igs(mgs),jy,kgs(mgs),lliq(il))), 0.0) 
        ENDDO
      ENDIF
      ENDDO
      
      ENDIF


       
!      CALL cld_cpu('Z-MOMENT-ZFAll')  

!  Find shape parameter rain


     IF ( lz(lr) > 1 .and. imurain == 3 ) THEN ! { RAIN SHAPE PARAM
          il = lr
          DO mgs = 1,ngscnt
         
         IF ( iresetmoments == 1 .or. iresetmoments == il  ) THEN
!         IF (  .false. .and. zx(mgs,lr) <= zxmin ) THEN
         IF ( zx(mgs,lr) <= zxmin ) THEN
           qx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
!         ELSEIF ( zx(mgs,lr) <= 0.0 .and. cx(mgs,lr) > 0.0 .and. qx(mgs,il) .gt. qxmin(il)) THEN
!           write(91,*) 'ZF: overdepletion of Zr: z,c,q = ',zx(mgs,il),cx(mgs,il),qx(mgs,il)
         ELSEIF ( cx(mgs,lr) <= cxmin ) THEN
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         ENDIF
         ENDIF
         
          
         
         IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
!          tmp = cx(mgs,lr)
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
!          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!             g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
!!             zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(1000.))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!!             an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
!          ENDIF
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
        ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)

            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)

!            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
!            vr = xv(mgs,lr)

!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z

            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              write(91,*) 'alpha = ',alpha(mgs,il)
             IF ( qx(mgs,il) < 1.e-8 ) THEN
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             ELSE
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
             ENDIF
          ENDIF
          
          IF ( zx(mgs,lr) > 0.0 ) THEN
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
            vr = xv(mgs,lr)
!            z = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
!            IF ( 100.*Abs(alp - alpha(mgs,lr))/Abs(alpha(mgs,lr)) .lt. 1. ) EXIT
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           print*,'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(mgs,lr),qr*1000,z*1.e18,vr,nrx


! check for artificial breakup (rain larger than allowed max size)
!        massfactor = ((3+alpha(mgs,lr))*xdia(mgs,lr,1)/xdia(mgs,lr,3))**3
        IF (  xv(mgs,il) .gt. xvmx(il) ) THEN
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          IF ( tmp < cx(mgs,il) ) THEN ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
!           IF ( alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax ) THEN
           IF ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) THEN

            IF ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alphar .and. alp <= rnumin ) THEN

             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             
             ENDIF
           ENDIF

           ENDIF
          ENDIF
           
          ELSE
          
           zx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
          
          ENDIF
          
          ENDDO
        ENDIF ! }
        

      IF ( ipconc .ge. 6 ) THEN

!  Find shape parameters for graupel,hail


        DO il = lr,lhab

        
        IF (  lz(il) .gt. 1 .and. ( .not. ( il == lr .and. imurain == 3 )) ) THEN
        
        DO mgs = 1,ngscnt

        IF ( il == lhl .and. lnhlf > 1 ) THEN
          IF ( cx(mgs,lhl) > cxmin ) THEN
            frac = chxf(mgs,lhl)/cx(mgs,lhl)
          ELSE
            frac = 0.0
          ENDIF
        ENDIF
        IF ( il == lh .and. lnhf > 1 ) THEN
          IF ( cx(mgs,lh) > cxmin ) THEN
            frach = chxf(mgs,lh)/cx(mgs,lh)
          ELSE
            frach = 0.0
          ENDIF
        ENDIF

         IF ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) THEN ! .or. qx(mgs,il) <= qxmin(il) ) THEN
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3 ) THEN
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ELSEIF ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3  ) THEN
!!            write(91,*) 'cx=0; qx,zx = ',1000.*qx(mgs,il),1.e18*zx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
         ENDIF

         IF (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF

         IF ( qx(mgs,il) .gt. qxmin(il) ) THEN

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        IF ( xv(mgs,il) .lt. xvmn(il)  ) THEN
!          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
!             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
!            
!          ENDIF
        ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            g1 = (6.0 + alphamax)*(5.0 + alphamax)*(4.0 + alphamax)/ &
     &            ((3.0 + alphamax)*(2.0 + alphamax)*(1.0 + alphamax))
            zx(mgs,il) = Max(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
!              write(91,*) 'ziegfall: something screwy with moments: il = ',il
!              write(91,*) 'q,n,z = ', 1.e3*qx(mgs,il),cx(mgs,il),zx(mgs,il)
!              write(91,*) 'alpha = ',alpha(mgs,il)

             IF ( qx(mgs,il) < 1.e-8 ) THEN
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             ELSE
         ! set values according to dBZ of -10
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
!               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
!              write(0,*) 'ziegfall: values of reset moments: il = ',il
!              write(0,*) 'q,n,z = ', 1.e3*qx(mgs,il),cx(mgs,il),zx(mgs,il)
            ENDIF
          ENDIF
         ENDIF

        IF ( qx(mgs,il) .gt. qxmin(il) .and. cx(mgs,il) .gt. 0.0 ) THEN
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          IF ( zx(mgs,il) .gt. 0. ) THEN
           
!            rd = z*(pi/6.*1000.)**2*chw/(0.224*(dn(igs(mgs),jy,kgs(mgs))*qr)**2)
            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((dn(igs(mgs),jy,kgs(mgs))*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rd,z,xv
           alp = Max( alphamin, Min( alphamax, alp ) )

         IF ( newton ) THEN
           DO i = 1,10
             IF ( i > 1 .and. Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = alp + ( galpha(alp) - rd )/dgalpha(alp)
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
           IF ( .not. ( alp < alphamax+0.001 .and. alp > alphamin - 0.001 ) ) THEN
             write(0,*) 'Zfall: problem with alp! myrank,alp,chw,qr,z,dn,xdn = ',my_rank,alp,chw,qr,z, &
     &                   dn(igs(mgs),jy,kgs(mgs)), xdn(mgs,il),il,rd

           alpha(mgs,il) = 0
           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           alp = Max( alphamin, Min( alphamax, alp ) )
           DO i = 1,10
            IF ( i > 1 .and. Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
!     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
             write(0,*) 'Zfall: old iteration, alp = ',alp
             call commasmpi_abort()
            ENDDO
           ENDIF
           
         ELSE
           DO i = 1,10
            IF ( i > 1 .and. Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
!     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
          ENDIF



! check for artificial breakup (graupel/hail larger than allowed max size)
        
        IF ( imaxdiaopt == 1 .or. il /= lr ) THEN
          xvbarmax = xvmx(il) 
        ELSEIF ( imaxdiaopt == 2 ) THEN ! test against maximum mass diameter
          xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ELSEIF ( imaxdiaopt == 3 ) THEN ! test against mass-weighted diameter
          xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ENDIF
        
        IF (  xv(mgs,il) .gt. xvbarmax ) THEN
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvbarmax, Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          IF ( tmp < cx(mgs,il) ) THEN ! breakup
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF
           
!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
           IF ( ( lrescalelow(il) .or. rescale_high_alpha ) .and.  &
     &        ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) THEN

             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( lrescalelow(il) .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) ) THEN

!!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*( 0.224*qr)*qr/chw
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!!!             z  = 1.e18*z1*(6./(pi*1000.))**2
!!             z  = z1*(6./(pi*1000.))**2
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
            ENDIF
           ENDIF
          ELSE
!             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*( 0.224*qr)*qr/chw
!!             z  = 1.e18*z1*(6./(pi*1000.))**2
!             z  = z1*(6./(pi*1000.))**2
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
          ENDIF
        ENDIF

        IF ( il == lhl .and. lnhlf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lhl) = frac*cx(mgs,lhl)
          an(igs(mgs),jy,kgs(mgs),lnhlf) = chxf(mgs,lhl)
        ENDIF
        IF ( il == lh .and. lnhf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lh) = frach*cx(mgs,lh)
          an(igs(mgs),jy,kgs(mgs),lnhf) = chxf(mgs,lh)
        ENDIF

        ENDDO ! mgs
        
        ENDIF ! lz(il) .gt. 1
        
        ENDDO ! il

        alphan(:,lr) = alpha(:,lr)
        alphan(:,lh) = alpha(:,lh)
        IF ( lf > 1 ) alphan(:,lf) = alpha(:,lf)
        IF ( lhl > 1 ) alphan(:,lhl) = alpha(:,lhl)

!      CALL cld_cpu('Z-MOMENT-ZFAll')  
          
      ENDIF

!      IF ( .not. ( ipconc == 5 .and. imydiagalpha > 0 ) ) THEN
!        alphan(:,:) = alpha(:,:)
!      ENDIF


!
!  Set density
!
#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: done setvtz"
#else
      if (ndebugzf .gt. 0 ) write(0,*)  my_rank,'ZIEGFALL: call setvtz'
#endif
!
      
      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,       &
     &                 ipconc,ndebugzf,ngs,nz,igs,kgs,cwnccn,fadvisc, &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,    &
     &                 itype1,itype2,temcg,infdo,alpha,alphan,axx,bxx,0)
!     &                 itype1,itype2,temcg,infdo,alpha,axh,bxh,axhl,bxhl)

#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: done setvtz"
#endif

!
! put fall speeds into the x-z arrays
!
      DO il = lc,lhab
      do mgs = 1,ngscnt
       
       vtmax = 150.0

!       IF ( lf > 1 .and. il == lf .and. ny <= 2 .and. igs(mgs) == 3 ) THEN
!         IF ( qx(mgs,il) > 1.e-6 ) THEN
!          write(0,*) 'FD vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          write(0,*) 'z,alp,xdn = ',zx(mgs,il),alpha(mgs,il),xdn(mgs,il)
!          write(0,*) 'axx,bxx,xmas = ',axx(mgs,il),bxx(mgs,il),xmas(mgs,il)
!          write(0,*) 'ax,bx,icdx = ',ax(il),bx(il),icdx
!         ENDIF
!       ENDIF
!       IF ( lh > 1 .and. il == lh .and. ny <= 2 .and. igs(mgs) == 3 ) THEN
!         IF ( qx(mgs,il) > 1.e-6 ) THEN
!          write(0,*) 'Graup vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          write(0,*) 'z,alp,xdn = ',zx(mgs,il),alpha(mgs,il),xdn(mgs,il)
!          write(0,*) 'axx,bxx,xmas = ',axx(mgs,il),bxx(mgs,il),xmas(mgs,il)
!          write(0,*) 'ax,bx,icdx = ',ax(il),bx(il),icdx
!         ENDIF
!       ENDIF
       
       
       IF ( vtxbar(mgs,il,2) .gt. vtxbar(mgs,il,1)  .or. &
     &      ( vtxbar(mgs,il,1) .gt. vtxbar(mgs,il,3) .and. vtxbar(mgs,il,3) > 0.0) ) THEN
          
          
          IF ( qx(mgs,il) > 1.e-4 .and.  &
     &        .not. ( il == lr .and. 1.e3*xdia(mgs,il,3) > 5.0 ) ) THEN
          write(0,*) 'my_rank, infdo,mgs = ',my_rank,infdo,lzr,mgs
          write(0,*) 'Moment problem with vtxbar for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
          IF ( il .ge. lg .or. il == lr ) THEN
            write(0,*) 'alpha = ',alpha(mgs,il)
          ENDIF
          ENDIF
          
          vtxbar(mgs,il,1) = Max( vtxbar(mgs,il,1), vtxbar(mgs,il,2) )
          vtxbar(mgs,il,3) = Max( vtxbar(mgs,il,3), vtxbar(mgs,il,1) )
          
       ENDIF

       
       IF ( vtxbar(mgs,il,1) .gt. vtmax .or. vtxbar(mgs,il,2) .gt. vtmax .or. &
     &      vtxbar(mgs,il,3) .gt. vtmax ) THEN
       
        IF ( ndebugzf >= 0 .and.  1.e3*qx(mgs,il) > 0.1 ) THEN
          write(0,*) 'my_rank,infdo = ',my_rank,infdo
          write(0,*) 'vtxbar exceeds vtmax for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
          write(0,*) 'vtmax = ',vtmax
!          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
          IF ( il .ge. lg ) THEN
            write(0,*) 'alpha = ',alpha(mgs,il)
          ENDIF
        ENDIF
        vtxbar(mgs,il,1) = Min(vtmax,vtxbar(mgs,il,1) )
        vtxbar(mgs,il,2) = Min(vtmax,vtxbar(mgs,il,2) )
        vtxbar(mgs,il,3) = Min(vtmax,vtxbar(mgs,il,3) )
        
!        call commasmpi_abort()
       ENDIF


       xvt(igs(mgs),kgs(mgs),1,il) = vtxbar(mgs,il,1)
       xvt(igs(mgs),kgs(mgs),2,il) = vtxbar(mgs,il,2)
       IF ( infdo .ge. 2 ) THEN
       xvt(igs(mgs),kgs(mgs),3,il) = vtxbar(mgs,il,3)
       ELSE
       xvt(igs(mgs),kgs(mgs),3,il) = 0.0
       ENDIF

!       xvt(igs(mgs),kgs(mgs),2,il) = xvt(igs(mgs),kgs(mgs),1,il)

      enddo
      ENDDO

      if (ndebugzf .gt. 0 ) write(0,*)  my_rank,'ZIEGFALL: COPIED FALL SPEEDS'


 9998 continue

      if (ndebugzf .gt. 0 ) write(0,*)  my_rank,'ZIEGFALL: DONE WITH LOOP'

#ifdef MPI
!      IF ( .false. ) THEN
      
      if ( kzbeg-1+kz .gt. nzend-kstag-1 .and. ixbeg-1+ix .gt. nxend-istag ) then
!      if ( kz .gt. nz-kstag-1 .and. ix .ge. nx-istag) then
        if ( ixbeg-1+ix .eq. nxend-1 ) then
         go to 1200
        elseif ( ix .ge. nx ) then
         go to 1200
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz 
      end if
      
!      ENDIF
!      if ( kz .gt. kze-1 .and. ix .gt. ixe ) then
!        go to 1200
!      else
!        nzmpb = kz 
!      end if
#else
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
        go to 1200
      else
        nzmpb = kz 
      end if
#endif

      if (ndebugzf .gt. 0 ) print*,'ZIEGFALL: SET NZMPB'

#ifdef MPI
!      IF ( .false. ) THEN
      
      if ( ix .ge. nx-1 ) then
       if ( ixbeg-1+ix .eq. nxend-1 ) then
        nxmpb = 1
       elseif ( ix .ge. nx ) then
        nxmpb = 1
       else
        nxmpb = ix+1
       endif
      else
       nxmpb = ix+1
      end if
      
!      ENDIF
!      if ( ix+1 .gt. ixe ) then
!       nxmpb = 1
!      else
!       nxmpb = ix+1
!      end if
#else
      if ( ix+1 .gt. nx-1 ) then
       nxmpb = 1
      else
       nxmpb = ix+1
      end if
#endif

      end do !! inumgs

      if (ndebugzf .gt. 0 ) print*,'ZIEGFALL: SET NXMPB'

 1200 continue


!       ENDDO ! ix
!      ENDDO ! kz


#ifdef MPI
      if (debug_mpi) write(0,*) "ZIEGFALL: EXITING SUBROUTINE, my_rank=",my_rank,' jgs = ',jgs
#else
      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL: EXITING SUBROUTINE"
#endif


      RETURN
      END subroutine ziegfall

! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------

      subroutine ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ixcol, &
     &  xvt, rhovtzx,                                           &
     &  an,dn,ipconc0,t0,t7,cwmasn,cwmasx,       &
     &  cwradn,                                   &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx,  &
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,ln,lz,lvol,lliq,                     &
     &  infdo,ildo,cwnccn)

! 12.16.2005: .F version use in transitional SWM model
!
! 10.10.2003: Added cimn and cimx to setting for cci and cip.
!
! TO DO LIST:
!
! need to set up values for:
!     :  cipdia,cidia,cwdia,cwmas,vtwbar,
!     :  rho0,temcg,cip,cci
!
! and need to put fallspeed values in cwvt etc.
!
       USE INDEX_MODULE, only: lt,lc,lr,li,lis,ls,lh,lhl,lf,lv,lg,lhab,lzr,ax,bx,lhw,lfw,lhlw, &
                               lnc,lnr,lni,lnis,lns,lnh,lnf,lnhl,cinu,dmuh,dnu,dmu,xnu,xmu,dmuhl,rnu,cnu,snu, &
                               lss,lsat,lsati,xvcmx,xvcmn,xvrmn,xvrmx,lccn,rnumin,rnumax, &
                               lqmx,nxtra,lvi,lvs,lvh,lvf,lvhl,lzi,lzs,lzr,lzh,lzf,lzhl,lsw, &
                               alphar,alphamin,alphamax, lnhf,lnhlf
      USE COMMASMPI_MODULE, only: my_rank
      USE CPUTIME_MODULE
      USE MICRO_MODULE

      implicit none
      integer ng1
      parameter(ng1 = 1)
      
      integer, intent(in) :: ixcol ! which column to return
      integer, intent(in) :: ildo
      real,    intent(in) :: cwnccn(nz)
      
      integer nx,ny,nz,nor,norz,ngt,jgs,na
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real dtp,dtz1
      
      real :: rhovtzx(nz,nx)
      
      integer ndebugzf
      parameter (ndebugzf = 0)

      integer ix,jy,kz,i,j,k,il
      integer infdo
!
!
      real xvt(nz+1,nx,3,lc:lhab) ! 1=mass-weighted, 2=number-weighted

      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)

      integer :: ngs
      integer :: ngscnt,mgs,ipconc0
!      parameter ( ngs=200 )
      
      real ::  qx(ngs,lv:lhab) 
      real ::  qxw(ngs,ls:lhab) 
      real ::  cx(ngs,lc:lhab) 
      real ::  xv(ngs,lc:lhab) 
      real ::  vtxbar(ngs,lc:lhab,3) 
      real ::  xmas(ngs,lc:lhab) 
      real ::  xdn(ngs,lc:lhab) 
      real ::  cdxgs(ngs,lc:lhab) 
      real ::  xdia(ngs,lc:lhab,3) 
      real ::  vx(ngs,li:lhab) 
      real ::  alpha(ngs,lc:lhab) 
      real ::  zx(ngs,lr:lhab) 
      real :: alphan(ngs,lc:lhab)

      real ::  chxf(ngs,lh:lhab)

      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real :: axx(ngs,lh:lhab), bxx(ngs,lh:lhab)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)

!
!   drag coefficients
!
      real cdx(lc:lhab)
!
! Fixed intercept values for single moment scheme
!
      real cno(lc:lhab)
      
      real cwccn0,cwmasn,cwmasx,cwradn
!      real cwc0

      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag
      parameter (kstag=1)

      integer igs(ngs),kgs(ngs)
      
      real rho0(ngs),temcg(ngs)

      real temg(ngs)
      
      real rhovt(ngs)
      
      real cwnc(ngs),cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
!      real cimasn,cimasx,
      real :: cnina(ngs),cimas(ngs)
      
      real :: cnostmp(ngs)

      logical ldovol, ldoliq
      integer, intent(in) :: lvol(lc:lhab)
      integer, intent(in) :: ln(lc:lhab)
      integer, intent(in) :: lz(lc:lhab)
      integer, intent(in) :: lliq(li:lhab)
      
      real, parameter :: advisc0 = 1.832e-05
      real, parameter :: tfr = 273.15
      real, parameter :: pi1 = 3.14159265 ! 4.0*atan(1.0)
      real, parameter :: pi = 3.14159265 ! 4.0*atan(1.0)
      real, parameter :: piinv = 1.0/pi1
      real, parameter :: cwc1 = 6.0/(pi1*1000.)

!      real pii
!
!
!  general constants for microphysics
!

! 
! Miscellaneous
!
      
      logical flag
      
    
      real chw, qr, z, rd, alp, z1, g1, vr, nrx, tmp, tmp2, tmpz,tmpmas,tmpc, frac, frach
      
      real vtmax
      real xvbarmax

      real, parameter ::  c1r=19.0, c2r=0.6, c3r=1.8, c4r=17.0   ! rain
      real, parameter ::  c1h=5.5, c2h=0.7, c3h=4.5, c4h=8.5   ! Graupel
      real, parameter ::  c1hl=3.7, c2hl=0.3, c3hl=9.0, c4hl=6.5, c5hl=1.0, c6hl=6.5 ! Hail

      integer l1, l2
      
      double precision :: dpt1, dpt2


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.


      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL1D: ENTERED SUBROUTINE"

! #####################################################################
! BEGIN EXECUTABLE
! #####################################################################
!

!  constants
!

      ldovol = .false.
            
      IF ( li .gt. 1 ) THEN
      DO il = li,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO
      ENDIF

      ldoliq = .false.
      IF ( ls .gt. 1 ) THEN
      DO il = ls,lhab
        ldoliq = ldoliq .or. ( lliq(il) .gt. 1 )
      ENDDO
      ENDIF
!
      jy = jgs
      nxmpb = ixcol
      nzmpb = 1
      nxz = 1*nz
!      ngs = nz
      numgs = 1

      IF ( ildo == 0 ) THEN
        l1 = lc
        l2 = lhab
      ELSE
        l1 = ildo
        l2 = ildo
      ENDIF


      do inumgs = 1,numgs
       ngscnt = 0


       do kz = nzmpb,nz-1
        do ix = ixcol,ixcol
        flag = .false.

        
        DO il = l1,l2
          flag =  flag .or. ( an(ix,jy,kz,il)  .gt. qxmin(il) ) 
        ENDDO

        if ( flag ) then
! load temp quantities

        ngscnt = ngscnt + 1
        igs(ngscnt) = ix
        kgs(ngscnt) = kz
        if ( ngscnt .eq. ngs ) goto 1100
        end if
        end do !!ix
        nxmpb = 1
       end do !! kz

!      if ( jy .eq. (ny-jstag) ) iend = 1

 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998
!
!  set temporaries for microphysics variables
!


!
!  Reconstruct various quantities 
!
      do mgs = 1,ngscnt

       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhovt(mgs) = rhovtzx(kgs(mgs),ixcol) !  Sqrt(rho00/rho0(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
       temcg(mgs) = temg(mgs) - tfr

        
!
      end do
!
! only need fadvisc for 
      IF ( lc .gt. 1 .and. (ildo == 0 .or. ildo == lc ) ) then
        do mgs = 1,ngscnt
         fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &   (temg(mgs)/296.0)**(1.5)
        end do
      ENDIF

      IF ( ipconc .eq. 0 ) THEN
      do mgs = 1,ngscnt
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      end do
      ENDIF


      IF ( ildo > 0 ) THEN
        vtxbar(:,ildo,:) = 0.0
      ELSE
        vtxbar(:,:,:) = 0.0
      ENDIF
      
!      do mgs = 1,ngscnt
!        qx(mgs,lv) = max(an(igs(mgs),jy,kgs(mgs),lv), 0.0) 
!      ENDDO
      DO il = l1,l2
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do
      
      cnostmp(:) = cno(ls)
      IF ( ipconc < 1 .and. lwsm6 .and. (ildo == 0 .or. ildo == ls )) THEN
        DO mgs = 1,ngscnt
          tmp = Min( 0.0, temcg(mgs) )
          cnostmp(mgs) = Min( 2.e8, 2.e6*exp(0.12*tmp) )
        ENDDO
      ENDIF


!
!  set concentrations
!
      cx(:,:) = 0.0
       IF ( lnhf > 1 .or. lnhlf > 1  ) chxf(:,:) = 0.0

      if ( ipconc .ge. 1 .and. li .gt. 1 .and. (ildo == 0 .or. ildo == li ) ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
        IF ( lis > 1 ) cx(mgs,lis) = Max(an(igs(mgs),jy,kgs(mgs),lnis), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 .and. lc .gt. 1 .and. (ildo == 0 .or. ildo == lc ) ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
!        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
       end do
      end if
      if ( ipconc .ge. 3 .and. lr .gt. 1 .and. (ildo == 0 .or. ildo == lr ) ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
!        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
!        ELSE
!          cx(mgs,lr) = Max( 0.0, cx(mgs,lr) )
!        ENDIF
       end do
      end if
      if ( ipconc .ge. 4  .and. ls .gt. 1 .and. (ildo == 0 .or. ildo == ls ) ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
!        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
!        ELSE
!          cx(mgs,ls) = Max( 0.0, cx(mgs,ls) )
!        ENDIF
       end do
      end if

      if ( ipconc .ge. 5  .and. lh .gt. 1 .and. (ildo == 0 .or. ildo == lh ) ) then
       do mgs = 1,ngscnt

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
!        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
!        ELSE
!          cx(mgs,lh) = Max( 0.0, cx(mgs,lh) )
!        ENDIF


        IF ( lnhf > 1 ) THEN
           chxf(mgs,lh) = Min(cx(mgs,lh), Max(an(igs(mgs),jy,kgs(mgs),lnhf), 0.0))
        ENDIF

       end do
      ENDIF

      if ( ipconc .ge. 5  .and. lf .gt. 1 .and. (ildo == 0 .or. ildo == lf ) ) then
       do mgs = 1,ngscnt

        cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)

       end do
      ENDIF

      if ( ipconc .ge. 5  .and. lhl .gt. 1 .and. (ildo == 0 .or. ildo == lhl ) ) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
!        IF ( qx(mgs,lhl) .le. qxmin(lhl) ) THEN
!          cx(mgs,lhl) = 0.0
!        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) THEN
!          qx(mgs,lhl) = 0.0
!        ELSE
!          cx(mgs,lhl) = Max( 0.0, cx(mgs,lhl) )
!        ENDIF


        IF ( lnhlf > 1 ) THEN ! number of hail from frozen drops
           chxf(mgs,lhl) = Min(cx(mgs,lhl), Max(an(igs(mgs),jy,kgs(mgs),lnhlf), 0.0))
        ENDIF

       end do
      end if

       
      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
!        IF ( ls .gt. 1 .and. lvs .eq. 0 ) xdn(mgs,ls) = xdn0(ls)
!        IF ( lh .gt. 1 .and. lvh .eq. 0 ) xdn(mgs,lh) = xdn0(lh)
        IF ( li .gt. 1 )  xdn(mgs,li) = xdn0(li)
        IF ( lis .gt. 1 )  xdn(mgs,lis) = xdn0(lis)
        IF ( ls .gt. 1 )  xdn(mgs,ls) = xdn0(ls)
        IF ( lh .gt. 1 )  xdn(mgs,lh) = xdn0(lh)
        IF ( lf .gt. 1 )  xdn(mgs,lf) = xdn0(lf)
        IF ( lhl .gt. 1 ) xdn(mgs,lhl) = xdn0(lhl)
      end do

!
! Set mean particle volume
!
      IF ( ldovol .and. (ildo == 0 .or. ildo >= li ) ) THEN
      
      vx(:,:) = 0.0
      
       DO il = l1,l2
        
        IF ( lvol(il) .ge. 1 ) THEN
        
          DO mgs = 1,ngscnt
            vx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lvol(il)), 0.0)
            IF ( vx(mgs,il) .gt. rho0(mgs)*qxmin(il)*1.e-3 .and. qx(mgs,il) .gt. qxmin(il) ) THEN
              xdn(mgs,il) = Min( xdnmx(il), Max( xdnmn(il), rho0(mgs)*qx(mgs,il)/vx(mgs,il) ) )
            ENDIF
          ENDDO
          
        ENDIF
      
       ENDDO
      
      ENDIF

      DO il = lg,lhab
      DO mgs = 1,ngscnt
        alpha(mgs,il) = dnu(il)
      ENDDO
      ENDDO
      
      IF ( imurain == 1 ) THEN
        alpha(:,lr) = alphar
      ELSEIF ( imurain == 3 ) THEN
        alpha(:,lr) = xnu(lr)
      ENDIF


      IF ( ipconc == 5 .and. imydiagalpha > 0 ) THEN
        DO mgs = 1,ngscnt
          IF ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > cxmin ) THEN
             xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))            ! 
             xdia(mgs,lr,3) = (xv(mgs,lr)*6.0*cwc1)**(1./3.) 
             alpha(mgs,lr) = Min(alphamax, c1r*tanh(c2r*(xdia(mgs,lr,3)*1000. - c3r)) + c4r)
          ENDIF
          IF ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) > cxmin ) THEN
             xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))            ! 
             xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
             alpha(mgs,lh) = Min(alphamax, c1h*tanh(c2h*(xdia(mgs,lh,3)*1000. - c3h)) + c4h)
          ENDIF
!        alpha(:,lr) = 0. ! 10.
!        alpha(:,lh) = 0. ! 10.
          IF ( lhl > 0 ) THEN
          IF ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) > cxmin ) THEN
             xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))            ! 
             xdia(mgs,lhl,3) = (xv(mgs,lhl)*6.*piinv)**(1./3.)
             IF ( xdia(mgs,lhl,3) < 0.008 ) THEN
               alpha(mgs,lhl) = Min(alphamax, c1hl*tanh(c2hl*(xdia(mgs,lhl,3)*1000. - c3hl)) + c4hl)
             ELSE
               alpha(mgs,lhl) = Min(alphamax, c5hl*xdia(mgs,lhl,3)*1000. + c6hl)
             ENDIF
          ENDIF
          ENDIF
        ENDDO
      ENDIF


!
! Set 6th moments
!
      IF ( ipconc .ge. 6 .or. lzr > 1) THEN
      
      zx(:,:) = 0.0
      
!      DO il = lr,lhab
       DO il = l1,l2
        
        IF ( lz(il) .ge. 1 ) THEN
        
          DO mgs = 1,ngscnt
            zx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0)
          ENDDO
          
        
        ENDIF
      
       ENDDO
      
      ENDIF
       

!
! Set liquid water fractions
!
      IF ( ldoliq ) THEN
      
!     DO il = ls,lhab
      DO il = l1,l2
      IF ( il >= ls ) THEN
      IF ( lliq(il) .gt. 1 ) THEN
        DO mgs = 1,ngscnt
          qxw(mgs,il) = max(min(qx(mgs,il),an(igs(mgs),jy,kgs(mgs),lliq(il))), 0.0) 
        ENDDO
      ENDIF
      ENDIF
      ENDDO
      
      ENDIF
       
!  Find shape parameter rain


     IF ( lz(lr) > 1 .and. (ildo == 0 .or. ildo == lr ) .and. imurain == 3  ) THEN ! { RAIN SHAPE PARAM
          il = lr
          DO mgs = 1,ngscnt
         
         IF ( iresetmoments == 1 .or. iresetmoments == il  ) THEN
!         IF (  .false. .and. zx(mgs,lr) <= zxmin ) THEN
         IF ( zx(mgs,lr) <= zxmin ) THEN
           qx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
!         ELSEIF ( zx(mgs,lr) <= 0.0 .and. cx(mgs,lr) > 0.0 .and. qx(mgs,il) .gt. qxmin(il)) THEN
!           write(91,*) 'ZF: overdepletion of Zr: z,c,q = ',zx(mgs,il),cx(mgs,il),qx(mgs,il)
         ELSEIF ( cx(mgs,lr) <= cxmin ) THEN
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         ENDIF
         ENDIF
         
          
         
         IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
!          tmp = cx(mgs,lr)
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
!          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!             g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
!!             zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(1000.))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!!             an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
!          ENDIF
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
        ENDIF

          IF ( zx(mgs,il) > zxmin .and. cx(mgs,il) <= cxmin ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)

            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

            IF ( cx(mgs,lr) <= cxmin ) THEN
            ! if resulting concentration is still too small, then zero out
              zx(mgs,lr) = 0.0
              qx(mgs,lr) = 0.0
              an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
              an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
              an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)

!            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
!            vr = xv(mgs,lr)

!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z

            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              write(91,*) 'alpha = ',alpha(mgs,il)
             IF ( qx(mgs,il) < 1.e-8 ) THEN
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             ELSE
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
             ENDIF
          ENDIF
          
          IF ( zx(mgs,lr) > 0.0 ) THEN
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
            vr = xv(mgs,lr)
!            z = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           IF ( z .gt. 0.0 ) THEN
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
!            IF ( 100.*Abs(alp - alpha(mgs,lr))/Abs(alpha(mgs,lr)) .lt. 1. ) EXIT
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
!           write(0,*) 'kz, alp, alpha(kz) = ',kz,alp,alpha(mgs,lr),qr*1000,z*1.e18,vr,nrx


! check for artificial breakup (rain larger than allowed max size)
        IF (  xv(mgs,il) .gt. xvmx(il) ) THEN
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          IF ( tmp < cx(mgs,il) ) THEN ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           DO i = 1,20
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
!           IF ( alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax ) THEN
           IF ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) THEN

            IF ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alphar .and. alp <= rnumin ) THEN

             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             
             ENDIF
           ENDIF

           ENDIF
          ENDIF
           
          ELSE
          
           zx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
          
          ENDIF
          
          ENDDO
        ENDIF ! }
        

      IF ( ipconc .ge. 6 ) THEN

!  Find shape parameters for graupel,hail

        DO il = lr,lhab

        IF ( lz(il) .gt. 1 .and. (ildo == 0 .or. ildo == il ) .and.  &
                ( .not. ( il == lr .and. imurain == 3 )) ) THEN
        
        DO mgs = 1,ngscnt

        IF ( il == lhl .and. lnhlf > 1 ) THEN
          IF ( cx(mgs,lhl) > cxmin ) THEN
            frac = chxf(mgs,lhl)/cx(mgs,lhl)
          ELSE
            frac = 0.0
          ENDIF
        ENDIF
        IF ( il == lh .and. lnhf > 1 ) THEN
          IF ( cx(mgs,lh) > cxmin ) THEN
            frach = chxf(mgs,lh)/cx(mgs,lh)
          ELSE
            frach = 0.0
          ENDIF
        ENDIF

!          IF ( igs(mgs) == 10 .and. il == lhl ) THEN
!            write(91,*) 'zf1d: k,q,c,z = ',kgs(mgs),qx(mgs,il),cx(mgs,il),zx(mgs,il)
!          ENDIF
         IF ( iresetmoments == 1 .or. iresetmoments == il  .or. iresetmoments == -1 ) THEN
         IF ( zx(mgs,il) <= zxmin ) THEN !  .and. qx(mgs,il) > 0.05e-3 ) THEN
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         ELSEIF ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         ELSEIF ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) THEN !  .and. qx(mgs,il) > 0.05e-3  ) THEN
!!            write(91,*) 'cx=0; qx,zx = ',1000.*qx(mgs,il),1.e18*zx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF
         ENDIF

         IF (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) THEN
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         ENDIF

         IF ( qx(mgs,il) .gt. qxmin(il) ) THEN

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        IF ( xv(mgs,il) .lt. xvmn(il)  ) THEN
!          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
!             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
!            
!          ENDIF
        ENDIF

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= cxmin ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

            IF ( cx(mgs,il) <= cxmin ) THEN
            ! if resulting number is still too small, then zero out
              cx(mgs,il) = 0.0
              zx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

              qx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

            IF ( zx(mgs,il) <= zxmin ) THEN
            ! if resulting reflectivity is still too small, then zero out
              cx(mgs,il) = 0.0
              zx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

              qx(mgs,il) = 0.0
              an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
            ENDIF

           ELSEIF ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
!              write(91,*) 'ziegfall: something screwy with moments: il = ',il
!              write(91,*) 'q,n,z = ', 1.e3*qx(mgs,il),cx(mgs,il),zx(mgs,il)
!              write(91,*) 'alpha = ',alpha(mgs,il)

             IF ( qx(mgs,il) < 1.e-8 ) THEN
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             ELSE
!              write(0,*) 'alpha = ',alpha(mgs,il)
         ! set values according to dBZ of -10
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
            ENDIF
          ENDIF
         ENDIF

        IF ( qx(mgs,il) .gt. qxmin(il) .and. cx(mgs,il) .gt. cxmin ) THEN
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          IF ( zx(mgs,il) .gt. 0. ) THEN
           
!            rd = z*(pi/6.*1000.)**2*chw/(0.224*(dn(igs(mgs),jy,kgs(mgs))*qr)**2)
            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((dn(igs(mgs),jy,kgs(mgs))*qr)**2)

           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           DO i = 1,10
            IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO



! check for artificial breakup (graupel/hail larger than allowed max size)
        
        IF ( imaxdiaopt == 1 .or. il /= lr ) THEN
          xvbarmax = xvmx(il) 
        ELSEIF ( imaxdiaopt == 2 ) THEN ! test against maximum mass diameter
          xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ELSEIF ( imaxdiaopt == 3 ) THEN ! test against mass-weighted diameter
          xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        ENDIF
        
        IF (  xv(mgs,il) .gt. xvbarmax ) THEN
          tmp = cx(mgs,il)
          tmp2 = xv(mgs,il)
          xv(mgs,il) = Min( xvbarmax, Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          tmpc = cx(mgs,il) ! rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))

          IF ( tmp < tmpc ) THEN ! breakup
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
            ! zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            ! check if incoming zx is consistent
            ! Z from incoming cx, qx, and alpha
            tmpz = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/tmp
            IF ( tmpz > zx(mgs,il) ) THEN
              tmpc = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/zx(mgs,il)
              cx(mgs,il) = Max(cx(mgs,il), tmpc)
              ! find cx that gives zx
            ENDIF
            zx(mgs,il) = g1/(pi/6.*xdn(mgs,il))**2 * ((rho0(mgs)*qx(mgs,il))**2)/cx(mgs,il)
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

            qr  = qx(mgs,il)
            chw = cx(mgs,il)
            z   = zx(mgs,il)

            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           DO i = 1,10
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

            
          ENDIF
        ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
           IF ( (rescale_low_alpha .or. rescale_high_alpha ) .and.  &
     &        ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) THEN

             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            ELSEIF ( rescale_low_alpha .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) ) THEN

!!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*( 0.224*qr)*qr/chw
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
            ENDIF
           ENDIF
          ELSE
          ENDIF
        ENDIF


        IF ( il == lhl .and. lnhlf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lhl) = frac*cx(mgs,lhl)
          an(igs(mgs),jy,kgs(mgs),lnhlf) = chxf(mgs,lhl)
        ENDIF
        IF ( il == lh .and. lnhf > 1 ) THEN
        ! update chxf in case cx has changed
          chxf(mgs,lh) = frach*cx(mgs,lh)
          an(igs(mgs),jy,kgs(mgs),lnhf) = chxf(mgs,lh)
        ENDIF

        ENDDO ! mgs
        
        ENDIF ! lz(il) .gt. 1
        
        ENDDO ! il

!      CALL cld_cpu('Z-MOMENT-ZFAll')  
          
      ENDIF

      IF ( lzhl > 1 ) THEN
        IF ( lhl .gt. 1 ) THEN
        
        ENDIF
      ENDIF


!        alphan(:,:) = alpha(:,:)
        alphan(:,lr) = alpha(:,lr)
        alphan(:,lh) = alpha(:,lh)
        IF ( lf > 1 ) alphan(:,lf) = alpha(:,lf)
        IF ( lhl > 1 ) alphan(:,lhl) = alpha(:,lhl)

!
!  Set density
!
      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: call setvtz'
!
!       call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,   &
!      &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,       &
!      &                 ipconc,ndebugzf,ngs,nz,igs,kgs,cwnccn,fadvisc, &
!      &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,    &
!      &                 itype1,itype2,temcg,infdo,alpha,alphan,axx,bxx)
      
      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,        &
     &                 ipconc,ndebugzf,ngs,nz,igs,kgs,cwnccn,fadvisc, &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,    &
     &                 itype1,itype2,temcg,infdo,alpha,alphan,axx,bxx,ildo)
!     &                 itype1,itype2,temcg,infdo,alpha,ildo,axh,bxh,axhl,bxhl)



!
! put fall speeds into the x-z arrays
!
      DO il = l1,l2
      do mgs = 1,ngscnt
       
       vtmax = 150.0

       
       IF ( vtxbar(mgs,il,2) .gt. vtxbar(mgs,il,1)  .or. &
     &      ( vtxbar(mgs,il,1) .gt. vtxbar(mgs,il,3) .and. vtxbar(mgs,il,3) > 0.0) ) THEN
          
          
#ifdef Z3MOM
!          IF ( qx(mgs,il) > 1.e-4 .and.  &
!     &        .not. ( il == lr .and. 1.e3*xdia(mgs,il,3) > 5.0 ) ) THEN
!          write(0,*) 'infdo,mgs = ',infdo,lzr,mgs
!          write(0,*) 'Moment problem with vtxbar for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
!          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
!          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
!          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
!          IF ( il .ge. lg .or. il == lr ) THEN
!            write(0,*) 'alpha = ',alpha(mgs,il)
!          ENDIF
!          ENDIF
#endif
          
          vtxbar(mgs,il,1) = Max( vtxbar(mgs,il,1), vtxbar(mgs,il,2) )
          vtxbar(mgs,il,3) = Max( vtxbar(mgs,il,3), vtxbar(mgs,il,1) )
          
       ENDIF

       
       IF ( vtxbar(mgs,il,1) .gt. vtmax .or. vtxbar(mgs,il,2) .gt. vtmax .or. &
     &      vtxbar(mgs,il,3) .gt. vtmax ) THEN
       
#ifdef Z3MOM
!        IF ( ndebugzf >= 0 .and.  1.e3*qx(mgs,il) > 0.1 ) THEN
!          write(0,*) 'infdo = ',infdo
!          write(0,*) 'Problem with vtxbar for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
!          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
!          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
!          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
!          IF ( il .ge. lg ) THEN
!            write(0,*) 'alpha = ',alpha(mgs,il)
!          ENDIF
!        ENDIF
#endif
        vtxbar(mgs,il,1) = Min(vtmax,vtxbar(mgs,il,1) )
        vtxbar(mgs,il,2) = Min(vtmax,vtxbar(mgs,il,2) )
        vtxbar(mgs,il,3) = Min(vtmax,vtxbar(mgs,il,3) )
        
!        call commasmpi_abort()
       ENDIF


       xvt(kgs(mgs),igs(mgs),1,il) = vtxbar(mgs,il,1)
       xvt(kgs(mgs),igs(mgs),2,il) = vtxbar(mgs,il,2)
       IF ( infdo .ge. 2 ) THEN
       xvt(kgs(mgs),igs(mgs),3,il) = vtxbar(mgs,il,3)
       ELSE
       xvt(kgs(mgs),igs(mgs),3,il) = 0.0
       ENDIF

!       xvt(kgs(mgs),igs(mgs),2,il) = xvt(kgs(mgs),igs(mgs),1,il)

      enddo
      ENDDO


      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: COPIED FALL SPEEDS'



 9998 continue

      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: DONE WITH LOOP'

      if ( kz .gt. nz-1 ) then
        go to 1200
      else
        nzmpb = kz 
      end if

      if (ndebugzf .gt. 0 ) write(0,*) 'ZIEGFALL: SET NZMPB'

      end do !! inumgs

      if (ndebugzf .gt. 0 ) write(0,*) 'ZIEGFALL: SET NXMPB'

 1200 continue


!       ENDDO ! ix
!      ENDDO ! kz


      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL: EXITING SUBROUTINE"


      RETURN
      END subroutine ziegfall1d

! #####################################################################
