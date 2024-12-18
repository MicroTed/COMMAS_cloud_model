#include "sam.def.h"
c
c  subroutine to calculate fall speeds of hydrometeors
c
      subroutine ddfallnx(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
     :  xvt,
     :  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx,
     :  rwmasn,rwmasx,cwradn,
     :  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn,
     :  itype1,itype2,infdo)
     
       USE INDEX_MODULE, only: lt,lc,lr,li,ls,lh,lhl,lip,lir,lgl,lgm,lgh,lf,lv,lhab,
     :                          lnc,lnr,lni,lns,lnh,lngl,lngm,lngh,lnhl,lnf,lnir,lnip
       USE COMMASMPI_MODULE
c
c 12.16.2005: .F version use in transitional SWM model
c
c 10.10.2003: Added cimn and cimx to setting for cci and cip.
c
c TO DO LIST:
c
c need to set up values for:
c     :  cipdia,cidia,cwdia,cwmas,vtwbar,
c     :  rho0,temcg,cip,cci
c
c and need to put fallspeed values in cwvt etc.
c
      
      implicit none

      integer ng1
      parameter(ng1 = 1)
      
      integer nx,ny,nz,nor,ngt,jgs,na

      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
c      real vt(nx,nz)
      real dtp,dz,dtz1
      
      real ccwmx
      real cwc1(nz), cwnccn(nz)
c      real qtmp1(nx,nz),qtmp2(nx,nz)
c      real cmax
c      real xfall(nx,ny)
      
      integer ndebug
      parameter (ndebug = -1)
      integer idx

      integer ix,jy,kz,i,j,k,item,il
      integer itype1,itype2,infdo
c
c  include file for mix ratio and charge indices
c
c
      real xvt(nx,nz,2,lc:lhab) ! 1=mass-weighted, 2=number-weighted

      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      
      real, allocatable :: qx(:,:) ! qx(ngs,lv:lhab)
      real, allocatable :: cx(:,:) ! qx(ngs,lv:lhab)
      real, allocatable :: xv(:,:)
      real, allocatable :: vtxbar(:,:,:)
      real, allocatable :: xmas(:,:)
      real, allocatable :: xdn(:,:)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real, allocatable :: xdia(:,:,:)

c
c   drag coefficients
c
      real cdx(lc:lhab)
c
c Fixed intercept values for single moment scheme
c
      real cno(lc:lhab)
      
      real cwccn,cwmasn,cwmasx,cimn,cimx,cwradn
      real rwmasn,rwmasx
      integer ngs1
      real cwc0

      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag,istag
      parameter (kstag=1, istag=1)
      
      

      integer ngs,ngscnt,mgs,ipconc
      parameter ( ngs=500 )
      integer igs(ngs),kgs(ngs)
      
      real rho0(ngs),temcg(ngs)

      real temg(ngs)
      
      real rhovt(ngs)
      
      real cwnc(ngs),cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
      real cimasn,cimasx,cnina(ngs),cimas(ngs)

      real poo, cp608, cp, cv,ar
      real aradcw, bradcw, cradcw, dradcw
      real br, bta1, cnit, dragh, dnz00, rho00, cs, ds
      real pi, pii, pid4, qccrit, qscrit
      real gf1, gf1p5, gf2,  gf3, gf3p5, gf4, gf5
      real gf6, gf7, gf4br, gf4ds, gf4p5, gf3ds, gf1ds, gr
      real gamma
c
c intercepts
c
c
c
c  constants
c
      real c1f3
c
c  general constants for microphysics
c
       real tfr, advisc0

c 
c Miscellaneous
c
c      character*80 infile, outfile
c      integer ifile
      integer ihabdo
      parameter(ihabdo = 0)
      integer ierr
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.

c #####################################################################
c BEGIN EXECUTABLE
c #####################################################################
c

c  constants
c
      poo = 1.0e+05
      cp608 = 0.608
      cp = 1004.0
      cv = 717.0
      ar = 841.99666  
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      br = 0.8
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
      rho00 = 1.225
      cs = 4.83607122
      ds = 0.25
c  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pi = 4.0*atan(1.0)
      pii = 1./pi
      pid4 = pi/4.0 
      qccrit = 2.0e-03
      qscrit = 6.0e-04
      gf1 = 1.0 ! gamma(1.0)
      gf1p5 = 0.8862269255 ! gamma(1.5)
      gf2 = 1.0 ! gamma(2.0)
      gf3 = 2.0 ! gamma(3.0)
      gf3p5 = 3.32335097 ! gamma(3.5)
      gf4 = 6.0 ! gamma(4.0)
      gf5 = 24.0 ! gamma(5.0)
      gf6 = 120.0 ! gamma(6.0)
      gf7 = 720.0 ! gamma(7.0)
      gf4br = gamma(4.0+br)
      gf4ds = gamma(4.0+ds)
      gf4p5 = gamma(4.0+0.5)
      gf3ds = gamma(3.0+ds)
      gf1ds = gamma(1.0+ds)
      gr = 9.8
      cwc0 = pii
      advisc0 = 1.832e-05
c
c  constants
c
      c1f3 = 1.0/3.0
c
c  general constants for microphysics
c
      tfr = 273.15
c
c  ci constants in mks units
c
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
c
c  Set terminal velocities...
c    also set drag coefficients
c
      
      jy = jgs
c      DO kz = 1,nz-1
c        DO ix = 1,nx-1

      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do 1000 inumgs = 1,numgs 
      ngscnt = 0
      do 1001 kz = nzmpb,nz-kstag-1 
      do 1002 ix = nxmpb,nx-istag

      if ( 
     >     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or. 
     >     an(ix,jy,kz,li)  .gt. qxmin(li)   .or.
     >     an(ix,jy,kz,lip) .gt. qxmin(li)   .or.
     >     an(ix,jy,kz,lir) .gt. qxmin(li)   .or.
     >     an(ix,jy,kz,lr)  .gt. qxmin(lr)   .or.
     >     an(ix,jy,kz,ls)  .gt. qxmin(ls)   .or.
     >     an(ix,jy,kz,lgl) .gt. qxmin(lh)   .or.
     >     an(ix,jy,kz,lgm) .gt. qxmin(lh)   .or.
     >     an(ix,jy,kz,lgh) .gt. qxmin(lh)   .or.
     >     an(ix,jy,kz,lf)  .gt. qxmin(lh)   .or.
     >     an(ix,jy,kz,lhl) .gt. qxmin(lh)   .or.
     <     an(ix,jy,kz,lh)  .gt. qxmin(lh) ) then
       ! load temp quantities
       
c      ngscnt = 1
c      mgs = 1
c      igs(mgs) = ix
c      kgs(mgs) = kz
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 1100
      end if
 1002 continue
      nxmpb = 1
 1001 continue
c      if ( jy .eq. (ny-jstag) ) iend = 1
 1100 continue
      if ( ngscnt .eq. 0 ) go to 9998
c
c  set temporaries for microphysics variables
c
      allocate (  qx(ngscnt,lv:lhab) )
      allocate (  cx(ngscnt,lc:lhab) )
      allocate (  xv(ngscnt,lc:lhab) )
      allocate (  vtxbar(ngscnt,lc:lhab,2) )
      allocate (  xmas(ngscnt,lc:lhab) )
      allocate (  xdn(ngscnt,lc:lhab) )
      allocate (  xdia(ngscnt,lc:lhab,3) )

      vtxbar(:,:,:) = 0.0
      
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

c
c  set concentrations
c
      cx(:,:) = 0.0
      
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
        cx(mgs,lir) = Max(an(igs(mgs),jy,kgs(mgs),lnir), 0.0)
        cx(mgs,lip) = Max(an(igs(mgs),jy,kgs(mgs),lnip), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
c        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
c        IF ( na .ge. lccn ) THEN
c         ccn(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
c        ELSE
c         ccn(mgs) = 0.0
c        ENDIF
       end do
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(ls) ) THEN
          cx(mgs,lr) = 0.0
        ELSEIF ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(ls) ) THEN
          qx(mgs,lr) = 0.0
        ELSE
          cx(mgs,lr) = Max( 1.e-3, cx(mgs,lr) )
        ENDIF
       end do
      end if
      if ( ipconc .ge. 4 ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
          cx(mgs,ls) = 0.0
        ELSEIF ( cx(mgs,ls) .eq. 0.0 .and. qx(mgs,ls) .lt. 3.0*qxmin(ls) ) THEN
          qx(mgs,ls) = 0.0
        ELSE
          cx(mgs,ls) = Max( 1.e-3, cx(mgs,ls) )
        ENDIF
       end do
      end if
      if ( ipconc .ge. 5 ) then
       do mgs = 1,ngscnt
        cx(mgs,lgl) = Max(an(igs(mgs),jy,kgs(mgs),lngl), 0.0)
        IF ( qx(mgs,lgl) .le. qxmin(lh) ) THEN
          cx(mgs,lgl) = 0.0
        ELSEIF ( cx(mgs,lgl) .eq. 0.0 .and. qx(mgs,lgl) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgl) = 0.0
        ELSE
          cx(mgs,lgl) = Max( 1.e-3, cx(mgs,lgl) )
        ENDIF

        cx(mgs,lgm) = Max(an(igs(mgs),jy,kgs(mgs),lngm), 0.0)
        IF ( qx(mgs,lgm) .le. qxmin(lh) ) THEN
          cx(mgs,lgm) = 0.0
        ELSEIF ( cx(mgs,lgm) .eq. 0.0 .and. qx(mgs,lgm) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgm) = 0.0
        ELSE
          cx(mgs,lgm) = Max( 1.e-3, cx(mgs,lgm) )
        ENDIF

        cx(mgs,lgh) = Max(an(igs(mgs),jy,kgs(mgs),lngh), 0.0)
        IF ( qx(mgs,lgh) .le. qxmin(lh) ) THEN
          cx(mgs,lgh) = 0.0
        ELSEIF ( cx(mgs,lgh) .eq. 0.0 .and. qx(mgs,lgh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgh) = 0.0
        ELSE
          cx(mgs,lgh) = Max( 1.e-3, cx(mgs,lgh) )
        ENDIF

        cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)
        IF ( qx(mgs,lf) .le. qxmin(lh) ) THEN
          cx(mgs,lf) = 0.0
        ELSEIF ( cx(mgs,lf) .eq. 0.0 .and. qx(mgs,lf) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lf) = 0.0
        ELSE
          cx(mgs,lf) = Max( 1.e-3, cx(mgs,lf) )
        ENDIF

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
          cx(mgs,lh) = 0.0
        ELSEIF ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lh) = 0.0
        ELSE
          cx(mgs,lh) = Max( 1.e-3, cx(mgs,lh) )
        ENDIF

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        IF ( qx(mgs,lhl) .le. qxmin(lh) ) THEN
          cx(mgs,lhl) = 0.0
        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lhl) = 0.0
        ELSE
          cx(mgs,lhl) = Max( 1.e-3, cx(mgs,lhl) )
        ENDIF
       end do
      end if
       
       
       
      do mgs = 1,ngscnt
        xdn(mgs,lir) = xdn0(lir)
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
        xdn(mgs,ls) = xdn0(ls)
        xdn(mgs,lgl) = xdn0(lgl)
        xdn(mgs,lgm) = xdn0(lgm)
        xdn(mgs,lgh) = xdn0(lgh)
        xdn(mgs,lf) = xdn0(lf)
        xdn(mgs,lh) = xdn0(lh)
        xdn(mgs,lhl) = xdn0(lhl)
      end do


c
c
c  Reconstruct various quantities 
c

      
      do mgs = 1,ngscnt
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhovt(mgs) = Sqrt(rho00/rho0(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
c      temg(mgs) = temcg(mgs) + tfr

      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*
     >  (temg(mgs)/296.0)**(1.5)
c
      end do
c
      IF ( ipconc .eq. 0 ) THEN
      do mgs = 1,ngscnt
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      end do
      ENDIF

c
c  Set density
c
      if (ndebug .gt. 0 ) print*,'Set density'
c
      call setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
     :                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,
     :                 ipconc,ndebug,ngs,nz,kgs,cwnccn,fadvisc,
     :                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,
     :                 itype1,itype2,temcg,infdo)



c
c put fall speeds into the x-z arrays
c
      DO il = lc,lhab
      do mgs = 1,ngscnt
      
       xvt(igs(mgs),kgs(mgs),1,il) = vtxbar(mgs,il,1)
       xvt(igs(mgs),kgs(mgs),2,il) = vtxbar(mgs,il,2)

c       xvt(igs(mgs),kgs(mgs),2,il) = xvt(igs(mgs),kgs(mgs),1,il)
      
      ENDDO
      end do

      deallocate ( qx )
      deallocate ( cx )
      deallocate ( xv )
      deallocate ( vtxbar )
      deallocate ( xmas )
      deallocate ( xdn )
      deallocate ( xdia )

 9998 continue
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
c
      go to 1200
      else
      nzmpb = kz 
      end if
      if ( ix+1 .gt. nx-1 ) then
      nxmpb = 1
      else
      nxmpb = ix+1
      end if
 1000 continue
 1200 continue


c       ENDDO ! ix
c      ENDDO ! kz
      
      
      


      RETURN
      END
