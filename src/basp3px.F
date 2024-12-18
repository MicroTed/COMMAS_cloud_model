c
c
c      subroutine getslb
c
c
c 11/3/99 modify to calculate values for an oblique slice
c
      subroutine getslb
     > (gd,idx,ib,jb,kb,ivint,nxp,nyp,nzp
     >  ,nxslb,nyslb,nzslb,istag,jstag,kstag
     >  ,id1,jd1,kd1,imx,jmx,kmx,imn,jmn,kmn
     <  ,xmx,xmn,cstag,plt2dp,plt2dn,uvec,vvec,wvec,uava,vava,
     :   nx,ny,nz,na,nor,
     :   st,u,v,w,km,dbz,xtra,precip,precip_old,
     :   uinit,vinit,sinit,elec,rho,
     :   pb,db,pn,dn,rh,the,tt0,ng,
     :   dx,dy,dz,gx,gy,gz)


c     :   nx,ny,nz,na,nor,
c     :   st,u%flt3d,v%flt3d,w%flt3d,km%flt3d,dbz%flt3d
c     :   uinit%flt1d,vinit%flt1d,elec)
c
      USE PLOT_MODULE
      USE MICRO_MODULE
      USE INDEX_MODULE
      USE GRID_MODULE
      use takcommon,       only: lmax,lfmax,kfmax,iimax,kimax,lr100,ls250
      use comm1,           only: rhof
      use comm2,           only: sr,sx,srf,fw,vw,bew,cdw,szr,szf,szi      
     &                            ,sri,sfi,di,sxf,sxi,sh,ff,vf,fi,vi       
      
      implicit none

      integer na,nor,nx,ny,nz,ng

      TYPE(GRID) :: gd
      TYPE(VARIABLE) :: xtra(nxtra)
      TYPE(VARIABLE) :: elec(neelec)
      TYPE(VARIABLE) :: precip(nprecip+neelec2d)
      TYPE(VARIABLE) :: dbz
      TYPE(VARIABLE) :: u, uinit
      TYPE(VARIABLE) :: v, vinit
      TYPE(VARIABLE) :: w
      TYPE(VARIABLE) :: gx(4), gy(4), gz(4)
      real           :: sinit(nz,na)

      real :: tt0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real :: precip_old(-ng+1:nx+ng,-ng+1:ny+ng,nprecip)
      
      integer lamcir,lamcs,lamcgl
      integer lamcgm,lamcgh,lamcf,lamch,lamcip,lamchl
      integer lamci,lamcr
      real ftol
      integer kz,jy,ix,i,j,k,k1,n,il,ip1,ip2
      real  maxdbz, x
      integer idx
      real cstag,uava,vava
      integer kstag,istag,jstag,ivint
      real density,densi
      integer ivb,ive,jvb,jve,kvb,kve,ii,jj,kk
      real xi,yj,delx,dely,f2,f3,f4
      integer imn,jmn,kmn,imx,jmx,kmx,jypu,jymu
      real facyv,facxv
      integer id1,jd1,kd1
      integer ixpv,ixmv,irdpre,iread,nxx,nyx,ityp
      real delh_w
      real dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2
      real xmx,xmn
      real xib,xie,yjb,yje,zkb,zke,phi,dd,d1pobl
      integer ipltyp,nl
      integer ia, cnt
      real pi,d2r

      character*40 filrad
      character*40 flnmpr
      integer mu,mv,mw
c      include 'sam.plot.dims'
c      parameter (nx=26,ny=26,nz=26,na=75,nv=3,nor=3,nh=1)
c      parameter (nxm=42,nym=42,nzm=42)
c      parameter (nx2=2,ny2=2,nz2=2,na2=5,nv2=3,nor2=3)
c      parameter (lt=1,le=2)
      real, allocatable :: pltdp(:,:,:)

c      include 'sam.index.ion.h'
c
c      integer lt,le
c      parameter (lt=1,le=2)
c
c   mixing ratio
c
c      integer lv,lc,lr,li,lir,ls
c      integer lgl,lgm,lgh,lf,lh
c      integer lip,lhl,lhab
c      
c      integer lqb,lqe

c
c  number concentration (was 16 to 27)
c
c      integer lnc,lnr,lni,lnir,lns
c      integer lngl,lngm,lngh,lnf
c      integer lnh,lnip,lnhl
c      integer lnshi
c      parameter (lnshi = 1) ! 64)
c      parameter (lnc=0+lnshi,lnr=1+lnshi,lni=2+lnshi,lnir=3+lnshi)
c      parameter (lns=4+lnshi,lngl=5+lnshi,lngm=6+lnshi,lngh=7+lnshi)
c      parameter (lnf=8+lnshi,lnh=9+lnshi,lnip=10+lnshi,lnhl=11+lnshi)

c  space charge (28 to 39)

c      integer lscw,lscr,lsci,lscir,lscs
c      integer lscgl,lscgm,lscgh,lscf,lsch
c      integer lscip,lschl

c      integer lscb,lsce

c      parameter (lscw=28,lscr=29,lsci=30,lscir=31,lscs=32)
c      parameter (lscgl=33,lscgm=34,lscgh=35,lscf=36,lsch=37)
c      parameter (lscip=38,lschl=39)
c
c  temperature (was 40 to 51)
c
      integer ltpc,ltpr,ltpi,ltpir,ltps
      integer ltpgl,ltpgm,ltpgh,ltpf,ltph
      integer ltpip,ltphl
      integer ltshi
      parameter (ltshi = 1) ! 52)
      parameter (ltpc=0+ltshi,ltpr=1+ltshi,ltpi=2+ltshi)
      parameter (ltpir=3+ltshi,ltps=4+ltshi,ltpgl=5+ltshi)
      parameter (ltpgm=6+ltshi,ltpgh=7+ltshi,ltpf=8+ltshi)
      parameter (ltph=9+ltshi,ltpip=10+ltshi,ltphl=11+ltshi)
c
c  mean cloud content exposure (was 52 to 63)
c
      integer lamcw ! ,lamcr,lamci,lamcir,lamcs
c      integer lamcgl,lamcgm,lamcgh,lamcf,lamch
c      integer lamcip,lamchl
      integer lashi
      parameter (lashi = 16)
      parameter (lamcw=0+lashi) ! ,lamcr=1+lashi,lamci=2+lashi)
      parameter (lamcir=3+lashi,lamcs=4+lashi,lamcgl=5+lashi)
      parameter (lamcgm=6+lashi,lamcgh=7+lashi,lamcf=8+lashi)
      parameter (lamch=9+lashi,lamcip=10+lashi,lamchl=11+lashi)
c
c
c  mean cloud content exposure (was 52 to 63)
c
c      integer lamcw,lamcr,lamci,lamcir,lamcs
c      integer lamcgl,lamcgm,lamcgh,lamcf,lamch
c      integer lamcip,lamchl
c      integer lashi
c      parameter (lashi = 16)
c      parameter (lamcw=0+lashi,lamcr=1+lashi,lamci=2+lashi)
c      parameter (lamcir=3+lashi,lamcs=4+lashi,lamcgl=5+lashi)
c      parameter (lamcgm=6+lashi,lamcgh=7+lashi,lamcf=8+lashi)
c      parameter (lamch=9+lashi,lamcip=10+lashi,lamchl=11+lashi)
c
c  tau  (was 64 to 75)
c
c      integer lrtc,lrtr,lrti,lrtir,lrts
c      integer lrtgl,lrtgm,lrtgh,lrtf,lrth
c      integer lrtip,lrthl
c      integer lrshi
c      parameter (lrshi = 31)
c      parameter (lrtc=0+lrshi) !,lrtr=1+lrshi,
c      parameter (lrti=1+lrshi,lrtip=2+lrshi)
c      parameter (lrtir=3+lrshi,lrts=4+lrshi,lrtgl=5+lrshi)
c      parameter (lrtgm=6+lrshi,lrtgh=7+lrshi,lrtf=8+lrshi)
      
      parameter (mu=1,mv=2,mw=3)
      integer nxp,nyp,nzp,nxslb,nyslb,nzslb,ib,jb,kb
      real cntz(nz)
      real plt2dp(nxslb,nyslb,nzslb)
      real plt2dn(nxslb,nyslb,nzslb)
      real uvec(nxslb,nyslb,nzslb)
      real vvec(nxslb,nyslb,nzslb)
      real wvec(nxslb,nyslb,nzslb)
      real st(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real km(nx,ny,nz)
c      real u(nx,ny,nz)
c      real v(nx,ny,nz)
c      real w(nx,ny,nz)
c      real uinit(nz)
c      real vinit(nz)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
c      real a2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2,na2)
c      real v2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2,nv2)
c      real p2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2)
c      real d2(-nor2:nx2+nor2,-nor2:ny2+nor2,-nor2:nz2+nor2)
c      real ab(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na+1)
c      real vb(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,nv)
      real pb(nz)
      real db(nz)
c      real temc(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real, allocatable :: tema(:,:,:) ! tema(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
c      real elec(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,6)
c      integer  iex,iey,iez,iemag,ipot,iezb
c      parameter (iex=1,iey=2,iez=3,iemag=4,ipot=5,iezb=6)
c      real algt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

c      real flashm(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) !  flash initiation points
c      real flashp(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) ! number of 'positive' segments (negative leaders)
c      real flashn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) ! number of 'negative' segments (positive leaders)
c      real drhom(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real rho(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      real rh(nx,ny,nz),the(nx,ny,nz)
c      real dbz(nx,ny,nz)      
      real, allocatable :: dum(:,:,:) ! (nx,ny,nz)
c      real divxy(nx,ny,nz),vorz(nx,ny,nz)
      real prt2d(nx,ny)
      real prti2d(nx,ny)
      real prtl2d(nx,ny)
      real cwrt2d(nx,ny)
      real cirt2d(nx,ny)
      real cirrt2d(nx,ny)
      real ciprt2d(nx,ny)
      real rwrt2d(nx,ny)
      real fwrt2d(nx,ny)
      real swrt2d(nx,ny)
      real glrt2d(nx,ny)
      real gmrt2d(nx,ny)
      real ghrt2d(nx,ny)
      real hlrt2d(nx,ny)
      real hxrt2d(nx,ny)
      real hwrt2d(nx,ny)
      real ptot2d(nx,ny)
      real ptoti2d(nx,ny)
      real ptotl2d(nx,ny)
      real cwtot2d(nx,ny)
      real citot2d(nx,ny)
      real cirtot2d(nx,ny)
      real ciptot2d(nx,ny)
      real rwtot2d(nx,ny)
      real fwtot2d(nx,ny)
      real swtot2d(nx,ny)
      real gltot2d(nx,ny)
      real gmtot2d(nx,ny)
      real ghtot2d(nx,ny)
      real hltot2d(nx,ny)
      real hxtot2d(nx,ny)
      real hwtot2d(nx,ny)
      logical fenele,feht1,femic,femica,febase
c      character*60 inname,outname,inname2
      character*60 inname3,inname4,basename
      real scmax,scmin

       integer :: l
!      integer nqsat
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati
!      parameter (fqsat=0.002,fqsati=1./fqsat)
!      parameter (fqsat=0.01,fqsati=1./fqsat)
!      real, save :: tabqvs(nqsat),tabqis(nqsat)
      integer, save :: itabflag = 0

      real cap,cp,cpi,rd,poo
     
      real pqs,qwv,temg,qvs,qis,temq,tmp,tmpn,qr,chw,area,areaw,z
      real :: masscum(45)
      real :: alpha1,alp,rdi,pi1,diam
      real qxmin(lc:lhab)
!      real    :: cxmin = 1.e-4  ! threshold cutoff for number concentration
!      real    :: zxmin = 1.e-28 ! threshold cutoff for reflectivity moment
      integer,save :: lf75 = 0 ! 75 micron radius; use as dividing line between ice crystal and graupel/frozen drop
      integer imx1,jmx1,kmx1
      real  :: buoymax
      integer ltemq
      logical test
      real cld
      real hwdn
      
      real tmp856, tmp857, tmp858

      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      
      integer :: istart, iskip  


c
c  common 'plot'
c
c     common /ppppp/ ab,db,vb,a,d,p,v,a2,d2,p2,v2,elec,sch,algt,hdnx,
c    > divxy,vorz,rh,the,dbz,dum,cntz,umove,vmove,ftol,
c    > dx,dy,dz,dxi2,dyi2,dzi2
c      common /ppppp/ ab,db,pb,vb,a,d,p,v,elec,
c     > umove,vmove,dx,dy,dz,dxi2,dyi2,dzi2,temc,flashp,flashn,
c     > drhom,flashm,rho,isum,dbz,the,rh

      common /oblique/xib,xie,yjb,yje,zkb,zke,ipltyp,phi,dd,nl,d1pobl
       
c       integer ipconc
      common /bastuff/femic,inname3


c      common /indices/lt,le,lv,lc,lr,li,lir,ls,lgl,lgm,lgh,
c     :    lf,lh,lip,lhl,lhab,lqb,lqe,
c     :    lscw,lscr,lsci,lscir,lscs,
c     :    lscgl,lscgm,lscgh,lscf,lsch,
c     :    lscip,lschl,lscb,lsce

c
      pi = 4.0*ATan(1.0)
      pi1 = pi
      d2r = pi/180.0

        IF ( lf75 < 1 .and. microp(1:3) == 'TAK') THEN
          IF ( taknicfd <= 1 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0075 ) lf75 = l
          enddo
          ELSEIF ( taknicfd == 2 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0025 ) lf75 = l
          enddo
          ELSEIF ( taknicfd == 3 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0125 ) lf75 = l
          enddo
          ENDIF
        ENDIF

      call setqxminz(qxmin)
      
!      qxmin(lh) = 1.e-8
!      IF ( lhl > 1 ) qxmin(lhl) = 1.e-8
      
      allocate( pltdp(nx+1,ny+1,nz+1) )
      
      DO kz=1,nz+1
      DO jy=1,ny+1
      DO ix=1,nx+1
        pltdp(ix,jy,kz) = 0.0
      ENDDO
      ENDDO
      ENDDO

      poo = 1.0e+05
      rd = 287.04
      cp = 1004.0
      cap = rd/cp 

      IF ( itabflag == 0 ) THEN
        itabflag = 1

        do l = 1,nqsat
          temq = 163.15 + (l-1)*fqsat
          tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
          tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
        end do
      
      ENDIF
      
      write(*,*) 'sp3p: dx,dy,dz',dx,dy,dz
      if ( dx .gt. 0. ) then
      dxi = 1.0/dx
      dyi = 1.0/dy
      dzi = 1.0/dz
      dxi2 = 0.5*dxi
      dyi2 = 0.5*dyi
      dzi2 = 0.5*dzi
      end if



c
c      do kz = 1,nz-1
c      do jy = 1,ny-1
c      do ix = 1,nx-1
c      xmax = max(xmax,algt(ix,jy,kz))
c      xmin = min(xmin,algt(ix,jy,kz))
c      end do
c      end do
c      end do
c      print*,'XMAX XMIN LIGHT ',xmax, xmin
c
      xmx = -1.e+5
      xmn = 1.e+5
      IF ( idx .eq. 856 .or. idx .eq. 857 .or. idx .eq. 858 ) THEN
       IF ( jyave .eq. 0 .and. ixave .eq. 0 ) THEN

       do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
         kz = (k-1) + kb
         do i = 1,nxp
          ix = (i-1) + ib
           IF ( idx .eq. 856)  pltdp(i,j,k) = 1.0e12*elec(icghis)%flt3d(ix,jy,kz)
           IF ( idx .eq. 857)  pltdp(i,j,k) = 1.0e12*elec(icghw)%flt3d(ix,jy,kz)
           IF ( idx .eq. 858)  pltdp(i,j,k) = 1.0e12*(elec(icghis)%flt3d(ix,jy,kz) + elec(icghw)%flt3d(ix,jy,kz))
         end do
        end do
       end do
      
       ELSE

       do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
         kz = (k-1) + kb
         do i = 1,nxp
          ix = (i-1) + ib

           pltdp(i,j,k) = 0.0
           tmp856 = 0.0
           tmp857 = 0.0
           tmp858 = 0.0
           IF (nyp .eq. 1) THEN
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
              tmp856 = tmp856 + elec(icghis)%flt3d(ix,ia,kz)
              tmp857 = tmp857 + elec(icghw)%flt3d(ix,ia,kz)
              tmp858 = tmp858 + (elec(icghis)%flt3d(ix,ia,kz) + elec(icghw)%flt3d(ix,ia,kz))
               cnt = cnt + 1
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
              tmp856 = tmp856 + elec(icghis)%flt3d(ia,jy,kz)
              tmp857 = tmp857 + elec(icghw)%flt3d(ia,jy,kz)
              tmp858 = tmp858 + (elec(icghis)%flt3d(ia,jy,kz) + elec(icghw)%flt3d(ia,jy,kz))
               cnt = cnt + 1
             ENDDO
           ELSE

              cnt = 1
              tmp856 = tmp856 + elec(icghis)%flt3d(ix,ia,kz)
              tmp857 = tmp857 + elec(icghw)%flt3d(ix,ia,kz)
              tmp858 = tmp858 + (elec(icghis)%flt3d(ix,ia,kz) + elec(icghw)%flt3d(ix,ia,kz))
           ENDIF

           IF ( cnt .gt. 0 ) THEN
            IF ( idx .eq. 856)  pltdp(i,j,k) = 1.0e12*tmp856/Float(cnt)
            IF ( idx .eq. 857)  pltdp(i,j,k) = 1.0e12*tmp857/Float(cnt)
            IF ( idx .eq. 858)  pltdp(i,j,k) = 1.0e12*tmp858/Float(cnt)
           ENDIF
           
         end do
        end do
       end do
       
       ENDIF
      ELSEIF ( idx .ge. 800 .and. idx .le. 858 .and. femic ) THEN
      allocate ( tema(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor) )
      write (*,*) idx
c      CALL flush(*)
c      CALL baread(nx,ny,nz,nor,na,st,tema,idx,inname3,ipconc)
      do j = 1,nyp
       jy = (j-1) + jb
       do k = 1,nzp
        kz = (k-1) + kb
        do i = 1,nxp
         ix = (i-1) + ib
         IF ( idx .lt. 840 .or. idx .gt. 846 ) THEN
           pltdp(i,j,k) = 1.0e12*tema(ix,jy,kz)
         ELSE
           pltdp(i,j,k) = tema(ix,jy,kz)
         ENDIF
        end do
       end do
      end do
      deallocate( tema )
      END IF

c
c      title(1) = 'POT. TEMP. PERURBATION (K)              '
      if ( idx .eq. 1 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + (st(ix,ia,kz,lt)- sinit(kz,lt) )
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + (st(ia,jy,kz,lt)- sinit(kz,lt) )
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSE
            pltdp(i,j,k) = (st(ix,jy,kz,lt)- sinit(kz,lt) )
           ENDIF
!      pltdp(i,j,k) = (st(ix,jy,kz,lt)- sinit(kz,lt) ) ! ab(ix,jy,kz,lt))
      enddo
      enddo
      enddo
      end if

c      title(33) = 'POT. TEMP. (K)              '
      if ( idx .eq. 33 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = (st(ix,jy,kz,lt) )
      ENDDO
      ENDDO
      ENDDO
      
      end if
c
c      title(2) = 'TEMPERATURE (C)              '
      if ( idx .eq. 2 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = tt0(ix,jy,kz) - 273.15
!      pltdp(i,j,k) = st(ix,jy,kz,lt)*
!     :     ((pn(ix,jy,kz)+pb(kz))/100000.0)**(287./1004.) - 273.15
c     pltdp(i,j,k) = temc(ix,jy,kz) - 273.15
c     >  (st(ix,jy,kz,lt))
c     > *((p(ix,jy,kz)+pb(ix,jy,kz))/1.0e5)**(287.04/1004.)
c     > -273.15
      enddo
      enddo
      enddo
      end if

c      title(34) = 'TEMPERATURE (K)              '
      if ( idx .eq. 34 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = tt0(ix,jy,kz)
!      pltdp(i,j,k) = st(ix,jy,kz,lt)*
!     :     ((pn(ix,jy,kz)+pb(kz))/100000.0)**(287./1004.)
      ENDDO
      ENDDO
      ENDDO
      end if
c
c      title(3) = 'PERTURBATION U -WIND (M/S)'
      if ( idx .eq. 3 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + cstag*(u%flt3d(ix,ia,kz)+u%flt3d(ix+istag,ia,kz))
     > -(uinit%flt1d(kz))
     > + umove - uava
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + cstag*(u%flt3d(ia,jy,kz)+u%flt3d(ia+istag,jy,kz))
     > -(uinit%flt1d(kz))
     > + umove - uava
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSE
      pltdp(i,j,k) =
     >  cstag*(u%flt3d(ix,jy,kz)+u%flt3d(ix+istag,jy,kz))
     > -(uinit%flt1d(kz))
     > + umove - uava
           ENDIF
!      pltdp(i,j,k) =
!     >  cstag*(u%flt3d(ix,jy,kz)+u%flt3d(ix+istag,jy,kz))
!     > -(uinit%flt1d(kz))
!     > + umove - uava
      enddo
      enddo
      enddo
      end if
c
c      title(4) = 'PERTURBATION V -WIND (M/S)'
      if ( idx .eq. 4 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) =
     >  cstag*(v%flt3d(ix,jy,kz)+v%flt3d(ix,jy+jstag,kz))
     > -(vinit%flt1d(kz))
     > + vmove - vava
      enddo
      enddo
      enddo
      end if
c
c      title(5) = 'PERTURBATION W -WIND (M/S)'
      if ( idx .eq. 5 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) =
     >  cstag*(w%flt3d(ix,jy,kz)+w%flt3d(ix,jy,kz+kstag))
      enddo
      enddo
      enddo
      end if
c
c      title(6) = 'TOTAL U -WIND (M/S)  '
      if ( idx .eq. 6 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + cstag*(u%flt3d(ix,ia,kz)+u%flt3d(ix+istag,ia,kz))
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + cstag*(u%flt3d(ia,jy,kz)+u%flt3d(ia+istag,jy,kz))
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSE
      pltdp(i,j,k) =
     >  cstag*(u%flt3d(ix,jy,kz)+u%flt3d(ix+istag,jy,kz))
           ENDIF
c     > + umove
      enddo
      enddo
      enddo
      end if
c
c      title(7) = 'TOTAL V -WIND (M/S)  '
      if ( idx .eq. 7 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) =
     >  cstag*(v%flt3d(ix,jy,kz)+v%flt3d(ix,jy+jstag,kz))
c     > + vmove
      enddo
      enddo
      enddo
      end if
c
c      title(8) = 'PERTURBATION PRESSURE (MB)              '
      if ( idx .eq. 8 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = pn(ix,jy,kz)*.01
      enddo
      enddo
      enddo
      end if

c      title(56) = 'PRESSURE (MB)'
      if ( idx .eq. 56 .or. idx == 57 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = (pb(kz)+pn(ix,jy,kz))*.01
      ENDDO
      ENDDO
      ENDDO
      end if
c
cylylylylylylylylylylylylylylylylylylyllylylylylylylylylylylylylylylyl
c
c      title(251) = 'HORIZ. W GRAD. (*1.0E+3 M/S**2)'
	if ( idx .eq. 251 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 1000.0*((w%flt3d(ix+1,jy,kz)-w%flt3d(ix-1,jy,kz))*dxi2
     ^		    +(w%flt3d(ix,jy+1,kz)-w%flt3d(ix,jy-1,kz))*dyi2)
      enddo
      enddo
      enddo
      do  j = 1,nyp
      do  k = 1,nzp
      do  i = 1,nxp
	if (i.eq.1) pltdp(i,j,k)=pltdp(i+1,j,k)
	if (i.eq.nxp) pltdp(i,j,k)=pltdp(i-1,j,k)
	if (j.eq.1) pltdp(i,j,k)=pltdp(i,j+1,k)
	if (j.eq.nyp) pltdp(i,j,k)=pltdp(i,j-1,k)
      enddo
      enddo
      enddo
	end if
c
c
c      title(252) = 'VERT. PRESSURE GRAD. (*1.0E+2 N/KG)'
	if ( idx .eq. 252 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
	density=(dn(ix,jy,kz)+dn(ix,jy,kz-1))/2.
	densi=1./(density*dz)
c
      pltdp(i,j,k) = -100.0*((pn(ix,jy,kz)-pn(ix,jy,kz-1))*densi)
      enddo
      enddo
      enddo

      do  j = 1,nyp
      do  k = 1,nzp
      do  i = 1,nxp
        if (i.eq.1) pltdp(i,j,k)=pltdp(i+1,j,k)
        if (i.eq.nxp) pltdp(i,j,k)=pltdp(i-1,j,k)
        if (j.eq.1) pltdp(i,j,k)=pltdp(i,j+1,k)
        if (j.eq.nyp) pltdp(i,j,k)=pltdp(i,j-1,k)
      enddo
      enddo
      enddo

      end if
c	
cylylylylylylylylylylylylylylylylylylylyllylylyllylylylylyllylylylyllyly
c      title(9) = 'MIXING COEF.(M2/S2) FROM TKE EQ.'
      if ( idx .eq. 9 ) then
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = km(ix,jy,kz)
      enddo
      enddo
      enddo
      end if
c
c      title(10) = 'WATER VAPOR MIXING RATIO (G/KG)         '
      if ( idx .eq. 10) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = (st(ix,jy,kz,lv)-ab(ix,jy,kz,lv))*1000.
      pltdp(i,j,k) = (st(ix,jy,kz,lv)-sinit(kz,lv))*1000.
      enddo
      enddo
      enddo
      end if
c
c      title(11) = 'CLOUD WATER MIXING RATIO (G/KG)         '
      if ( idx .eq. 11 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.,0.0)
      enddo
      enddo
      enddo
      end if
c
c      title(12) = 'RAIN WATER MIXING RATIO (G/KG)          '
      if ( idx .eq. 12 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lr)*1000.,0.0)
      enddo
      enddo
      enddo
      end if
c
c      title(13) = 'CLOUD ICE MIXING RATIO (G/KG)           '
      if ( idx .eq. 13 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
       enddo
      enddo
      enddo

      end if
c
c      title(14) = 'SNOW AGGREGATE MIXING RATIO (G/KG)      '
      if ( idx .eq. 14 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,ls)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(15) = 'TOTAL HAIL MIXING RATIO (G/KG)  '
      if ( idx .eq. 15 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(16) = 'CLOUD WATER-ICE (G/KG)                  '
      if ( idx .eq. 16 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max((st(ix,jy,kz,lc)+st(ix,jy,kz,li))*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(17) = 'VISIBLE CLOUD (G/KG)'
      if ( idx .eq. 17 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( microp(1:3) == 'TAK' ) THEN
      pltdp(i,j,k) = max((st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                   st(ix,jy,kz,ls))*1000.,0.0)      
      ELSEIF ( na .ge. 40 ) THEN
      pltdp(i,j,k) = max((st(ix,jy,kz,lc)+st(ix,jy,kz,li)
     :                   +st(ix,jy,kz,lip)
     >                   +st(ix,jy,kz,ls))*1000.,0.0)
      ELSE
      pltdp(i,j,k) = max((st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                   st(ix,jy,kz,ls))*1000.,0.0)      
      ENDIF
      enddo
      enddo
      enddo

      end if
c
c      title(18) = 'TOTAL ICE (G/KG)'
      if ( idx .eq. 18 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      tmp = st(ix,jy,kz,li)+st(ix,jy,kz,ls)+st(ix,jy,kz,lh)
       IF ( lhl > 1 ) tmp = tmp + st(ix,jy,kz,lhl)
      pltdp(i,j,k) = max(tmp*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(19) = 'TOTAL LIQUID WATER (G/KG)'
      if ( idx .eq. 19 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max((st(ix,jy,kz,lc)+
     >                   st(ix,jy,kz,lr))*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(20) = 'TOTAL PRECIPITATION WATER (G/KG)'
      if ( idx .eq. 20 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max((st(ix,jy,kz,ls)+st(ix,jy,kz,lr)+
     >                   st(ix,jy,kz,lh))*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
c      title(21) = 'TOTAL CONDENSATE (G/KG)'
      if ( idx .eq. 21 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      tmp = st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                    st(ix,jy,kz,lr)+st(ix,jy,kz,ls)+
     >                    st(ix,jy,kz,lh)
      IF ( lhl > 1 ) tmp = tmp + st(ix,jy,kz,lhl)
      pltdp(i,j,k) = max(tmp*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 22 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + Max(0.0, dbz%flt3d(ix,ia,kz) )
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + Max(0.0, dbz%flt3d(ia,jy,kz) )
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)
           ELSE
             pltdp(i,j,k) = dbz%flt3d(ix,jy,kz)
           ENDIF
      enddo
      enddo
      enddo

      end if


      if ( idx .eq. 35 ) then ! composite reflectivity
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
       maxdbz = -100.
        DO k1 = 1,nz-1
         maxdbz = Max( maxdbz, dbz%flt3d(ix,jy,k1) )
        ENDDO
      
      pltdp(i,j,k) = maxdbz
      
      ENDDO
      ENDDO
      ENDDO
      
      end if

      if ( idx .eq. 38 ) then ! iwp ice water path (graupel)
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
        x = 0.0
        DO k1 = 1,nz-1
          IF ( tt0(ix,jy,k1) <= 263.15 ) THEN
            x = x + st(ix,jy,k1,lh)*dn(ix,jy,k1)/gz(3)%flt1d(k1)
          ENDIF
        ENDDO
      
      pltdp(i,j,k) = x
      
      ENDDO
      ENDDO
      ENDDO
      
      end if


c
c      title(39) = 'buoyancy'
      if ( idx .eq. 39 ) then
      imx1 = 0
      jmx1 = 0
      kmx1 = 0
      buoymax = -1.0e30
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      tmp = st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                    st(ix,jy,kz,lr)+st(ix,jy,kz,ls)+
     >                    st(ix,jy,kz,lh)
      IF ( lhl > 1 ) tmp = tmp + st(ix,jy,kz,lhl)
      
      pltdp(i,j,k) = 100.*9.8*( (st(ix,jy,kz,lt) - sinit(kz,lt))/sinit(kz,lt) +  
     >  0.61*(st(ix,jy,kz,lv) - sinit(kz,lv)) - max(tmp,0.0) )

        IF ( pltdp(i,j,k) > buoymax .and. tt0(ix,jy,kz) < 253.15 ) THEN
          imx1 = ix
          jmx1 = jy
          kmx1 = kz
          buoymax = pltdp(i,j,k)
        ENDIF

      ENDDO
      ENDDO
      ENDDO

        ix = ihack
        IF ( jhack > 0 ) THEN
          jy = jhack
        ELSE
          jy = jb
        ENDIF
        kz = khack

          pqs = 380.0/(pn(ix,jy,kz)+pb(kz))
          qwv = max(st(ix,jy,kz,lv),0.0)
c          temg = temc(ix,jy,kz)
          temg = tt0(ix,jy,kz)
!          temg = st(ix,jy,kz,lt)*
!     :     ((pn(ix,jy,kz)+pb(kz))/100000.0)**(287./1004.) ! temc(ix,jy,kz)

          ltemq = (temg-163.15)/fqsat+1.5
          ltemq = Min( nqsat, Max(1,ltemq) )
          qvs = pqs*exp(caw*(temg-273.15)/(temg-cbw)) ! tabqvs(ltemq)
          qis = pqs*tabqis(ltemq)
          hwdn = dn(ix,jy,kz)

        tmp = st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                    st(ix,jy,kz,lr)+st(ix,jy,kz,ls)+
     >                    st(ix,jy,kz,lh)
       IF ( lhl > 1 ) tmp = tmp + st(ix,jy,kz,lhl)
        write(0,*) 'point output hack at ix,jy,kz = ',ix,jy,kz
        write(0,*) 'B,B1,Load = ',
     >   100.*9.8*( (st(ix,jy,kz,lt) - sinit(kz,lt))/sinit(kz,lt) +  
     >  0.61*(st(ix,jy,kz,lv) - sinit(kz,lv)) - max(tmp,0.0) ),
     >   100.*9.8*( (st(ix,jy,kz,lt) - sinit(kz,lt))/sinit(kz,lt) +  
     >  0.61*(st(ix,jy,kz,lv) - sinit(kz,lv))),
     >   100.*9.8*(  - max(tmp,0.0) )
        write(0,*) 'qc,qr,qi+qs,qh,tot = ',1000.*st(ix,jy,kz,lc),1000.*st(ix,jy,kz,lr),
     >   1000.*(st(ix,jy,kz,li)+st(ix,jy,kz,ls)), 1000.*st(ix,jy,kz,lh),1000.*tmp
        write(0,*) 'content: qc,qr,qi+qs,qh,tot = ',1000.*st(ix,jy,kz,lc)*hwdn,1000.*st(ix,jy,kz,lr)*hwdn,
     >   1000.*(st(ix,jy,kz,li)+st(ix,jy,kz,ls))*hwdn, 1000.*st(ix,jy,kz,lh)*hwdn,1000.*tmp*hwdn
        write(0,*) 'pert theta, theta-v =',(st(ix,jy,kz,lt) - sinit(kz,lt)),
     >     st(ix,jy,kz,lt)*(1.0 + 0.61*st(ix,jy,kz,lv)) - 
     >                     sinit(kz,lt)*(1.0 + 0.61*sinit(kz,lv))
        write(0,*) 'w at w-point = ', w%flt3d(ix,jy,kz)
        write(0,*) 'w at th-point = ',0.5*(w%flt3d(ix,jy,kz)+w%flt3d(ix,jy,kz+kstag)) ! w%flt3d(ix,jy,kz)
        write(0,*) 'SSw, SSi = ', 100.*qwv/qvs - 100., 100.*qwv/qis - 100.
        write(0,*) 'ccw = ', 1.e-6*st(ix,jy,kz,lnc)
      end if


c
c      title(40) = 'therm buoyancy'
      if ( idx .eq. 40 ) then
      imx1 = 0
      jmx1 = 0
      kmx1 = 0
      buoymax = -1.0e30
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      pltdp(i,j,k) = 100.*9.8*( (st(ix,jy,kz,lt) - sinit(kz,lt))/sinit(kz,lt) +  
     >  0.61*(st(ix,jy,kz,lv) - sinit(kz,lv)) )


      ENDDO
      ENDDO
      ENDDO

      end if

c
c      title(41) = 'load buoyancy'
      if ( idx .eq. 41 ) then
      imx1 = 0
      jmx1 = 0
      kmx1 = 0
      buoymax = -1.0e30
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      tmp = st(ix,jy,kz,lc)+st(ix,jy,kz,li)+
     >                    st(ix,jy,kz,lr)+st(ix,jy,kz,ls)+
     >                    st(ix,jy,kz,lh)
      IF ( lhl > 1 ) tmp = tmp + st(ix,jy,kz,lhl)
      
      pltdp(i,j,k) = 100.*9.8*( - max(tmp,0.0) )


      ENDDO
      ENDDO
      ENDDO

      end if

c
c      title(28) = 'hack point'
      if ( idx .eq. 28) then
      imx1 = 0
      jmx1 = 0
      kmx1 = 0
      buoymax = -1.0e30
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      pltdp(i,j,k) = 0.0
      
      IF ( ix == ihack .and. kz == khack ) THEN
         pltdp(i,j,k) = 10.
      ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF

      if ( idx .eq. 332 .or. idx == 333) then ! flash extent density (8 km)
      IF ( iflshr < 1 ) THEN
        write(0,*) 'flshr not available! iflshr = ',iflshr
        STOP
      ELSE
      
      IF ( iflshr8km > 1 .and. idx == 332 ) THEN
        maxdbz = 0
      DO jy = 4,ny-1,8
       DO ix = 4,nx-1,8
         DO j = Max(1,jy-4),Min(ny-1,jy+3)
         DO i = Max(1,ix-4),Min(nx-1,ix+3)
          pltdp(i,j,nz)  = precip(iflshr8km)%flt2d(ix,jy)
          maxdbz = max(maxdbz,pltdp(i,j,nz))
!          pltdp(i,j,nz)  = precip(iflshr)%flt2d(i,j)
         ENDDO
         ENDDO
       ENDDO
       ENDDO
       write(0,*) 'max value for iflshr8km = ',maxdbz
     
      ELSE
        istart = Nint(4000./dx)
        iskip  = Nint(8000./dx)
      ! assuming 1km for this hack!
!      DO jy = 3,ny-1,8
!       DO ix = 3,nx-1,8
      DO jy = istart,ny-1,iskip
       DO ix = istart,nx-1,iskip
        maxdbz = 0
!         DO j = Max(1,jy-3),Min(ny-1,jy+4)
!         DO i = Max(1,ix-3),Min(nx-1,ix+4)
         DO j = Max(1,jy-istart),Min(ny-1,jy+istart-1)
         DO i = Max(1,ix-istart),Min(nx-1,ix+istart-1)
          maxdbz = Max( maxdbz,  precip(iflshr)%flt2d(i,j) )
         ENDDO
         ENDDO
         
!         DO j = jy-3,jy+4
!         DO i = ix-3,ix+4
!         DO j = Max(1,jy-3),Min(ny-1,jy+4)
!         DO i = Max(1,ix-3),Min(nx-1,ix+4)
         DO j = Max(1,jy-istart),Min(ny-1,jy+istart-1)
         DO i = Max(1,ix-istart),Min(nx-1,ix+istart-1)
          pltdp(i,j,nz)  = maxdbz
         ENDDO
         ENDDO
         
       ENDDO
      ENDDO
      
      ENDIF
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      
      pltdp(i,j,k) = pltdp(ix,jy,nz) ! precip(iflshr)%flt2d(ix,jy)
      
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      end if


      if ( idx .eq. 331 ) then ! flash extent density
      IF ( iflshr < 1 ) THEN
        write(0,*) 'flshr not available! iflshr = ',iflshr
        STOP
      ELSE
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      
      pltdp(i,j,k) = precip(iflshr)%flt2d(ix,jy)
!      IF ( pltdp(i,j,k) <= 0.0 ) THEN
!        IF ( j == 1 ) THEN
!         pltdp(i,j,k) = 5.*(Mod(i,5))
!        ELSEIF ( i == 1  ) THEN
!         pltdp(i,j,k) = 5.*(Mod(j,5))
!        ENDIF
!      ENDIF
      
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      end if

c      title(36) = 'THETA-V (K)'
      if ( idx .eq. 36 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lt)*(1.0 + 0.61*st(ix,jy,kz,lv)) ! ab(ix,jy,kz,lt))
      
      enddo
      enddo
      enddo
      end if

c      title(37) = 'THETA-V Perturbation (K)'
      if ( idx .eq. 37 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,Min(nz-1,nzp)
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lt)*(1.0 + 0.61*st(ix,jy,kz,lv)) - 
     :                     sinit(kz,lt)*(1.0 + 0.61*sinit(kz,lv))
      
      enddo
      enddo
      enddo
      end if

c
      if(idx.eq.23 .or. idx.eq.32) then
      uvec(:,:,:) = 0.0
      vvec(:,:,:) = 0.0
      wvec(:,:,:) = 0.0
      IF ( ipltyp .ne. 5 ) THEN
        IF ( idx .eq. 23 ) THEN
          write(*,*) 'UAVA,VAVA,IVINT:,',uava,vava,ivint
          write(*,*) 'umove, vmove = ',umove,vmove
        ENDIF
      if ( nxp .eq. 1 ) then
      ivb = 1
      ive = 1
      else
      ivb = 2
      ive = nxp-1
      end if
      if ( nyp .eq. 1 ) then
      jvb = 1
      jve = 1
      else
      jvb = 2
      jve = nyp-1
      end if
      if ( nzp .eq. 1 ) then
      kvb = 1
      kve = 1
      else
      kvb = 1
      kve = nzp-1
      end if
      do j = jvb,jve
      jy = (j-1) + jb
      jj = jy
      if ( nyp .eq. 1 ) jj = 1
      do k = kvb,kve
      kz = (k-1) + kb
      kk = kz
      if ( nzp .eq. 1 ) kk = 1
      do i = ivb,ive
      ix = (i-1) + ib
      ii = ix
      if ( nxp .eq. 1 ) ii = 1
c      uvec(i,j,k) = 999.0
c      vvec(i,j,k) = 999.0
c      wvec(i,j,k) = 999.0
      uvec(i,j,k) = 0.0
      vvec(i,j,k) = 0.0
      wvec(i,j,k) = 0.0
c      if ( 
c     >     ((ii-1)/ivint*ivint .eq. (ii-1) .and.
cc    >      (ii-1) .ne. 1  .and.
c     >      (ii-1) .ne. nxp) .and.
c     >     ((jj-1)/ivint*ivint .eq. (jj-1) .and.
cc    >      (jj-1) .ne. 1  .and.
c     >      (jj-1) .ne. nyp) .and.
c     >     ((kk-1)/ivint*ivint .eq. (kk-1) .and.
cc    >      (kk-1) .ne. 1  .and.
c     >      (kk-1) .ne. 1) 
c     <   )  then
cc
      if ( idx.eq.23 ) then
      uvec(i,j,k) = 0.5*(u%flt3d(ix,jy,kz)+u%flt3d(ix+istag,jy,kz))
c     >            - 0.5*(vb(ix,jy,kz,1)+vb(ix+istag,jy,kz,1))
     >            + umove - uava - umvnew
      vvec(i,j,k) = 0.5*(v%flt3d(ix,jy,kz)+v%flt3d(ix,jy+jstag,kz))
c     >            - 0.5*(vb(ix,jy,kz,2)+vb(ix,jy+jstag,kz,2))
     >            + vmove - vava - vmvnew
      wvec(i,j,k) = 0.5*(w%flt3d(ix,jy,kz)+w%flt3d(ix,jy,kz+kstag))
      
      IF ( radarmask .eq. 1 ) THEN
        IF ( dbz%flt3d(ix,jy,kz) .lt. lowdbz ) THEN
c     :       dbz(i+1,j,k) .lt. lowdbz .and.
c     :       dbz(i-1,j,k) .lt. lowdbz .and.
c     :       dbz(i,j+1,k) .lt. lowdbz .and.
          uvec(i,j,k) = 0.0
          vvec(i,j,k) = 0.0
          wvec(i,j,k) = 0.0
        ENDIF
      ELSEIF ( radarmask .eq. 2 ) THEN
        cld = st(ix,jy,kz,lc)
        IF ( microp(1:3) .eq. 'LFO' .or. microp(1:1) .eq. 'Z' .or. microp(1:1) .eq. 'T') THEN
          cld = cld + st(ix,jy,kz,li)
        ELSEIF ( microp(1:5) .eq. 'ICE10' ) THEN
          cld = cld + st(ix,jy,kz,li) + st(ix,jy,kz,lip) + st(ix,jy,kz,lir)
        ENDIF
        test = dbz%flt3d(ix,jy,kz) .ge. lowdbz .or. cld .gt. 0.01e-3

        IF ( .not. test ) THEN
          uvec(i,j,k) = 0.0
          vvec(i,j,k) = 0.0
          wvec(i,j,k) = 0.0
        ENDIF
      
      ENDIF
      
c      pltdp(i,j,k) = Sqrt(uvec(i,j,k)**2 + vvec(i,j,k)**2
c     :  + wvec(i,j,k)**2  )
c      pltdp(i,j,k) = wvec(i,j,k)
      end if
c
      if ( idx.eq.32 ) then
      uvec(i,j,k) = 0.5*(elec(iex)%flt3d(ix,jy,kz) + elec(iex)%flt3d(ix+istag,jy,kz))*1.e-3
      vvec(i,j,k) = 0.5*(elec(iey)%flt3d(ix,jy,kz) + elec(iey)%flt3d(ix,jy+jstag,kz))*1.e-3
      wvec(i,j,k) = 0.5*(elec(iez)%flt3d(ix,jy,kz) + elec(iez)%flt3d(ix,jy,kz+kstag))*1.e-3

c      vvec(i,j,k) = elec(iey)%flt3d(ix,jy,kz)*1.e-3
c      wvec(i,j,k) = elec(iez)%flt3d(ix,jy,kz)*1.e-3
c      pltdp(i,j,k) = Sqrt(elec(iex)%flt3d(ix,jy,kz)**2 
c     :  + elec(iey)%flt3d(ix,jy,kz)**2
c     :  + elec(iez)%flt3d(ix,jy,kz)**2  )*1.e-3
c      pltdp(i,j,k) = elec(iez)%flt3d(ix,jy,kz)*1.e-3
      end if
      
c      ENDIF

      enddo
      enddo
      enddo
c

      ELSEIF ( ipltyp .eq. 5 )  THEN
c interpolate z-component into wvec
       IF ( idx .eq. 23 ) THEN
         DO kz=1,nz
         DO jy=1,ny
         DO ix=1,nx
          pltdp(ix,jy,kz) = 0.5*(w%flt3d(ix,jy,kz)+w%flt3d(ix,jy,kz+kstag))
         ENDDO
         ENDDO
         ENDDO
       ELSEIF ( idx.eq.32 ) THEN
         DO kz=1,nz
         DO jy=1,ny
         DO ix=1,nx
          pltdp(ix,jy,kz) = elec(iez)%flt3d(ix,jy,kz)*1.e-3
         ENDDO
         ENDDO
         ENDDO
       ENDIF
       DO j=1,1
       DO k=1,nzp
         kz = (k-1) + kb
       DO i=1,nl
         xi = xib + (i-1 + d1pobl/dx)*Cos(phi*d2r)
         yj = yjb + (i-1 + d1pobl/dx)*Sin(phi*d2r)*(dx/dy)
         ix = Int(xi)
         jy = Int(yj)
         delx = xi - Float(ix)
         dely = yj - Float(jy)
         f2 = (1.0 - delx)*pltdp(ix,jy,kz) + delx*pltdp(ix+1,jy,kz) 
         f3 = (1.0 - delx)*pltdp(ix,jy+1,kz) + delx*pltdp(ix+1,jy+1,kz)
         f4 = (1.0 - dely)*f2 + dely*f3
         wvec(i,j,k) = f4 
         if (wvec(i,j,k).lt.xmn) then
          xmn = wvec(i,j,k)
          imn = ix
          jmn = jy
          kmn = kz
         end if
         if (wvec(i,j,k).gt.xmx) then
           xmx = wvec(i,j,k)
           imx = ix
           jmx = jy
           kmx = kz
         end if
       END DO
       END DO
       END DO


c rotate axes and interpolate horizontal component into uvec
       IF ( idx .eq. 23 ) THEN
         DO kz=1,nz
         DO jy=1,ny
         DO ix=1,nx
          pltdp(ix,jy,kz) = 
     :     0.5*(u%flt3d(ix,jy,kz)+u%flt3d(ix,jy,kz+kstag))*Cos(phi*d2r)
     :  +  0.5*(v%flt3d(ix,jy,kz)+v%flt3d(ix,jy,kz+kstag))*Sin(phi*d2r)
         ENDDO
         ENDDO
         ENDDO
       ELSEIF ( idx.eq.32 ) THEN
         DO kz=1,nz
         DO jy=1,ny
         DO ix=1,nx
         pltdp(ix,jy,kz) = elec(iex)%flt3d(ix,jy,kz)*1.e-3*Cos(phi*d2r)
     :                   + elec(iey)%flt3d(ix,jy,kz)*1.e-3*Sin(phi*d2r)
         ENDDO
         ENDDO
         ENDDO
       ENDIF
       DO j=1,1
       DO k=1,nzp
         kz = (k-1) + kb
       DO i=1,nl
         xi = xib + (i-1 + d1pobl/dx)*Cos(phi*d2r)
         yj = yjb + (i-1 + d1pobl/dx)*Sin(phi*d2r)*(dx/dy)
         ix = Int(xi)
         jy = Int(yj)
         delx = xi - Float(ix)
         dely = yj - Float(jy)
         f2 = (1.0 - delx)*pltdp(ix,jy,kz) + delx*pltdp(ix+1,jy,kz) 
         f3 = (1.0 - delx)*pltdp(ix,jy+1,kz) + delx*pltdp(ix+1,jy+1,kz)
         f4 = (1.0 - dely)*f2 + dely*f3
         uvec(i,j,k) = f4 

      
      IF ( radarmask .ge. 1 ) THEN
         f2 = (1.0 - delx)*dbz%flt3d(ix,jy,kz) + delx*dbz%flt3d(ix+1,jy,kz) 
         f3 = (1.0 - delx)*dbz%flt3d(ix,jy+1,kz) + delx*dbz%flt3d(ix+1,jy+1,kz)
         f4 = (1.0 - dely)*f2 + dely*f3
       IF ( radarmask .eq. 1 ) THEN
        IF ( f4 .lt. lowdbz ) THEN
c     :       dbz(i+1,j,k) .lt. lowdbz .and.
c     :       dbz(i-1,j,k) .lt. lowdbz .and.
c     :       dbz(i,j+1,k) .lt. lowdbz .and.
          uvec(i,j,k) = 0.0
          wvec(i,j,k) = 0.0
        ENDIF
      ELSEIF ( radarmask .eq. 2 ) THEN
        cld = st(ix,jy,kz,lc)
        IF ( microp(1:3) .eq. 'LFO' .or. microp(1:1) .eq. 'Z' ) THEN
          cld = cld + st(ix,jy,kz,li)
        ELSEIF ( microp(1:5) .eq. 'ICE10' ) THEN
          cld = cld + st(ix,jy,kz,li) + st(ix,jy,kz,lip) + st(ix,jy,kz,lir)
        ENDIF
        test = f4 .ge. lowdbz .or. cld .gt. 0.01e-3

        IF ( .not. test ) THEN
          uvec(i,j,k) = 0.0
          wvec(i,j,k) = 0.0
        ENDIF
      
       ENDIF
       ENDIF

       END DO
       END DO
       END DO

      END IF ! ( ipltyp .ne. 5 )

      end if  ! (idx.eq.23 .or. idx.eq.32)


c
      if ( idx .eq. 24 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = rh(ix,jy,kz)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 25 ) then ! horizontal divergence/convergence
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 
     >  (u%flt3d(ix+1,jy,kz)-u%flt3d(ix,jy,kz))*dxi
     > +(v%flt3d(ix,jy+1,kz)-v%flt3d(ix,jy,kz))*dyi
      pltdp(i,j,k) = pltdp(i,j,k)*1000
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 27) then ! windspeed
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = sqrt(
     >  ((0.5*(u%flt3d(ix+istag,jy,kz)+u%flt3d(ix,jy,kz)))**2) +
     >  ((0.5*(v%flt3d(ix,jy,kz)+v%flt3d(ix,jy+jstag,kz)))**2) )
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 26) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      jypu = min(jy+1,ny-jstag)
      jymu = max(jy-1,1)
      facyv = 1.0
      if ( jy .eq. 1 .or. jy .eq. ny-jd1 ) facyv = 2.0
      ixpv = min(ix+1,nx-istag)
      ixmv = max(ix-1,1)
      facxv = 1.0
      if ( ix .eq. 1 .or. ix .eq. nx-id1 ) facxv = 2.0
      pltdp(i,j,k) =
     >  0.5*( (v%flt3d(ixpv,jy,kz)+v%flt3d(ixpv,jy+1,kz))
     >       -(v%flt3d(ixmv,jy,kz)+v%flt3d(ixmv,jy+1,kz)) )*dxi2*facxv
     > -0.5*( (u%flt3d(ix,jypu,kz)+u%flt3d(ix+1,jypu,kz))
     >       -(u%flt3d(ix,jymu,kz)+u%flt3d(ix+1,jymu,kz)) )*dyi2*facyv
      pltdp(i,j,k) = pltdp(i,j,k)*1000.
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 29 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = the(ix,jy,kz)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 30 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      delh_w = 
     >  (w%flt3d(ix+1,jy,kz)-w%flt3d(ix-1,jy,kz))*dxi2
     > +(w%flt3d(ix,jy+1,kz)-w%flt3d(ix,jy-1,kz))*dyi2
      pltdp(i,j,k) = delh_w*1.e3
      enddo
      enddo
      enddo

      end if
c
c    
c
      irdpre=0
      if ( (idx .ge. 261 .and. idx .le. 275) .or. (idx .ge. 581 .and. idx .le. 595) ) then ! {
      read(22,*) iread,flnmpr
      if ( iread .eq. 1 ) then ! {
      open(unit=210,file=flnmpr,form='unformatted')
      read(210) nxx,nyx,rwtot2d,rwrt2d
      close (unit=210)
      end if  ! }
      if ( iread .eq. 2 ) then ! {
      open(unit=210,file=flnmpr,form='unformatted')
      read(210) nxx,nyx,citot2d,cirt2d,swtot2d,swrt2d
     >                  ,cwtot2d,cwrt2d,rwtot2d,rwrt2d
     >                  ,swtot2d,swrt2d,hwtot2d,hwrt2d
     >                  ,prt2d,ptot2d,prti2d,ptoti2d
     >                  ,prtl2d,ptotl2d
      close (unit=210)
      end if  ! }
      if ( iread .eq. 3 ) then  ! {
      open(unit=210,file=flnmpr,form='unformatted')
      read(210) nxx,nyx,citot2d,cirt2d,ciptot2d,ciprt2d
     >                  ,cirtot2d,cirrt2d,swtot2d,swrt2d
     >                  ,cwtot2d,cwrt2d,rwtot2d,rwrt2d
     >                  ,swtot2d,swrt2d,gltot2d,glrt2d
     >                  ,gmtot2d,gmrt2d,ghtot2d,ghrt2d
     >                  ,fwtot2d,fwrt2d,hwtot2d,hwrt2d
     >                  ,prt2d,ptot2d,prti2d,ptoti2d
     >                  ,prtl2d,ptotl2d
      close (unit=210)
      end if  ! }
      end if  ! }
c
c
!     title(281)  = 'PRECIP PRECIP RATE (MM/HR)'
!     title(282)  = 'PRECIP-ICE PRECIP RATE (MM/HR)'
!     title(283)  = 'Rain Precip Rate (MM/HR)'
!     title(284)  = 'Graupel Precip Rate (MM/HR)'
!     title(285)  = 'Hail Precip Rate (MM/HR)'
!     title(286)  = 'Frozen drop Precip Rate (MM/HR)'
!     title(287)  = 'Hail % from Frozen drops Precip Rate fraction (mass)'
!
!     title(291)  = 'PRECIP PRECIP Accum. (mm)'
!     title(292)  = 'PRECIP-ICE PRECIP RATE (mm)'
!     title(293)  = 'Rain Accum. (mm)'
!     title(294)  = 'Graupel Precip Accum. (mm)'
!     title(295)  = 'Hail Precip Accum. (mm)'
!     title(296)  = 'Hail from Frozen drop Precip Accum. (mm)'
!     title(297)  = 'Hail % from Frozen drops Precip Accum. fraction (#)'

!    iprecip = Max( 0, GET_VARIABLE_INDEX(gd, 'RAIN_RAT', 1)  )
!    iprecip + 1 = hail_rate (hail + graupel + fd)
!    iprecip + 2 = rain_acc
!    iprecip + 3 = hail_acc (hail + graupel + fd)
!    iphailacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACCH', 1) - iprecip + 1  )
!    iphailfacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_ACCHF', 1) - iprecip + 1  )
!    iphailnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACC', 1) - iprecip + 1  )
!     iphailfnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_NACCF', 1) - iprecip + 1  )
!    iphaildiam = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_DIAM', 1) - iprecip + 1  )
!
!    ipgracc = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_ACC', 1) - iprecip + 1  )
!    ipgrnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_NACC', 1) - iprecip + 1  )
!    ipgrdiam = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_DIAM', 1) - iprecip + 1  )
!
!    ipfdacc = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_ACC', 1) - iprecip + 1  )
!    ipfdnumacc = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_NACC', 1) - iprecip + 1  )
!    ipfddiam = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_DIAM', 1) - iprecip + 1  )

      if ( idx .eq. 283 ) then ! rain rate
      IF ( iprecip < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       pltdp(i,j,k) = precip(1)%flt2d(ix,jy)*3600.0 ! convert 1/s to 1/hr
      
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if

      if ( idx .eq. 285 ) then ! hail rate
      IF ( iprecip < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       pltdp(i,j,k) = precip(2)%flt2d(ix,jy)*3600.0
      
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if

      if ( idx .eq. 291 ) then ! rain accumulation
      IF ( iprecip < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
       IF ( precip(3)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(3)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif

      if ( idx == 293 ) then ! rain accumulation
      IF ( iprainacc2 < 1 ) THEN
        write(0,*) 'rainacc2 not available! iprainacc2 = ',iprainacc2
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
       IF ( precip(iprainacc2)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(iprainacc2)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif

      if ( idx .eq. 288 ) then ! rain-auto accumulation
      IF ( iprecip < 1 .or. iprainaccauto < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       IF ( precip(iprainaccauto)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(iprainaccauto)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
      
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif

      if ( idx .eq. 289 ) then ! rain-melt accumulation
      IF ( iprecip < 1 .or. iprainaccmelt < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       IF ( precip(iprainaccmelt)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(iprainaccmelt)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
      
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif

      if ( idx .eq. 290 ) then ! rain-shed accumulation
      IF ( iprecip < 1 .or. iprainaccshed < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       IF ( precip(iprainaccshed)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(iprainaccshed)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
      
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif


      if ( idx .eq. 295 .or. idx == 301 ) then ! hail accum
      IF ( iprecip < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      write(0,*) 'iphailacc,iprecip = ',iphailacc,iprecip
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
         IF ( idx == 295 ) THEN
           pltdp(i,j,k) = precip(iphailacc)%flt2d(ix,jy)
         ELSE
           pltdp(i,j,k) = precip(iphailacc)%flt2d(ix,jy) - precip_old(ix,jy,iphailacc)
         ENDIF
      
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if

      if ( idx == 308 ) then ! graupel accumulation (grid motion)
      IF ( ipgracc2 < 1 ) THEN
        write(0,*) 'ipgracc2 not available! ipgracc2 = ',ipgracc2
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
       IF ( precip(ipgracc2)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(ipgracc2)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif

      if ( idx == 309 ) then ! hail accumulation (grid motion)
      IF ( iphailacc2 < 1 ) THEN
        write(0,*) 'iphailacc2 not available! iphailacc2 = ',iphailacc2
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
       IF ( precip(iphailacc2)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = precip(iphailacc2)%flt2d(ix,jy) ! *3600.0 ! convert 1/s to 1/hr
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF
      endif


      if ( idx .eq. 296 ) then ! hail from fd accum
      IF ( iprecip < 1 ) THEN
        write(0,*) 'precip not available! iprecip = ',iprecip
        STOP
      ELSE
      write(0,*) 'iphailfacc,iprecip = ',iphailfacc,iprecip
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       pltdp(i,j,k) = precip(iphailfacc)%flt2d(ix,jy)
      
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if

      if ( idx .eq. 287 ) then ! Mass fraction of hail from FD accum
       write(0,*) 'iphailacc,iphailfacc = ',iphailacc,iphailfacc
      IF ( iphailacc < 1 .or. iphailfacc < 1 ) THEN
        write(0,*) 'Cannot do hail fraction from FD!  '
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       IF ( precip(iphailacc)%flt2d(ix,jy) > 0.01 ) THEN
         pltdp(i,j,k) = 100.*precip(iphailfacc)%flt2d(ix,jy)/precip(iphailacc)%flt2d(ix,jy)
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if

      if ( idx .eq. 297 .or. idx == 300 ) then ! Number fraction of hail from FD accum
       IF ( precipgroundrelative ) THEN
         write(0,*) 'iphailnumacc2,iphailfnumacc2 = ',iphailnumacc2,iphailfnumacc2
         ip1 = iphailnumacc2
         ip2 = iphailfnumacc2
       ELSE
         write(0,*) 'iphailnumacc,iphailfnumacc = ',iphailnumacc,iphailfnumacc
         ip1 = iphailnumacc
         ip2 = iphailfnumacc
       ENDIF
      IF ( ip1 < 1 .or. ip2 < 1 ) THEN
        write(0,*) 'Cannot do hail fraction from FD!  '
        STOP
      ELSE
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib
      
      
       IF ( idx == 297 ) THEN
        IF ( precip(ip1)%flt2d(ix,jy) > 1.0 .and. precip(ip2)%flt2d(ix,jy) >= 0.0 ) THEN
          pltdp(i,j,k) = 100.*precip(ip2)%flt2d(ix,jy)/precip(ip1)%flt2d(ix,jy)
        ELSE
          pltdp(i,j,k) = -0.01
!         pltdp(i,j,k) = 0.0
        ENDIF
       ELSE
         tmp  = precip(ip1)%flt2d(ix,jy) - precip_old(ix,jy,ip1)
         tmpn = precip(ip2)%flt2d(ix,jy) - precip_old(ix,jy,ip2)
!         tmpn = precip_old(ix,jy,ip1)
!         tmp  = precip_old(ix,jy,ip2)
        IF ( tmp > 1.0 .and. tmp > 0.0 ) THEN
          pltdp(i,j,k) = 100.*tmpn/tmp !  precip(ip2)%flt2d(ix,jy)/precip(ip1)%flt2d(ix,jy)
        ELSE
          pltdp(i,j,k) = -0.01
        ENDIF
       ENDIF
       ENDDO; ENDDO; ENDDO
      ENDIF

      end if
      
      IF  ( idx .eq. 298 ) THEN ! 'Hail mean diameter at sfc. (mm) (thistory)
        IF ( iphaildiam < 1 .and. ipgrdiam < 1 ) THEN
          write(0,*) 'Cannot plot hail diameter: idx,iphaildiam,ipgrdiam = ',idx, iphaildiam,ipgrdiam
        ELSE
        
        do j = 1,nyp
          jy = (j-1) + jb
          do k = 1,nzp
            kz = (k-1) + kb
            do i = 1,nxp
             ix = (i-1) + ib
      
      
            IF ( iphaildiam > 0 ) THEN
              IF ( precip(iphaildiam)%flt2d(ix,jy) > 0.0 ) THEN
                pltdp(i,j,k) = precip(iphaildiam)%flt2d(ix,jy) ! already in units of mm
              ELSE
                pltdp(i,j,k) = -0.01
!               pltdp(i,j,k) = 0.0
              ENDIF
            ELSEIF (ipgrdiam > 0 ) THEn
              IF ( precip(ipgrdiam)%flt2d(ix,jy) > 0.0 ) THEN
                pltdp(i,j,k) = precip(ipgrdiam)%flt2d(ix,jy) ! already in units of mm
              ELSE
                pltdp(i,j,k) = -0.01
!               pltdp(i,j,k) = 0.0
              ENDIF
            ENDIF
            ENDDO
           ENDDO
          ENDDO
        
        ENDIF
      ENDIF

      IF  ( idx .eq. 299 ) THEN ! 'Hail mean diameter at sfc. (mm) (total)'
      IF ( iphailnumacc < 1 .or. iphailacc < 1 ) THEN
        write(0,*) 'Cannot do hail diameter (299)!  '
        ELSE
        
        do j = 1,nyp
          jy = (j-1) + jb
          do k = 1,nzp
            kz = (k-1) + kb
            do i = 1,nxp
             ix = (i-1) + ib
      
      
       IF ( precip(iphailnumacc)%flt2d(ix,jy)*50.**2 > 10. .and. precip(iphailacc)%flt2d(ix,jy) > 0.01 ) THEN
         tmp = precip(iphailacc)%flt2d(ix,jy)/precip(iphailnumacc)%flt2d(ix,jy)
         pltdp(i,j,k) = 1000.*(tmp*6.0/(3.14159*900.0) )**(1.0/3.0)
       ELSE
         pltdp(i,j,k) = -0.01
!         pltdp(i,j,k) = 0.0
       ENDIF
            ENDDO
           ENDDO
          ENDDO
        
        ENDIF
      ENDIF

      IF  ( idx .eq. 310 ) THEN ! 'graupel mean diameter at sfc. (mm) (thistory)
        IF ( ipgrdiam < 1 ) THEN
          write(0,*) 'Cannot plot mean hail diameter: idx,ipgrdiam = ',idx,ipgrdiam
        ELSE
        
        do j = 1,nyp
          jy = (j-1) + jb
          do k = 1,nzp
            kz = (k-1) + kb
            do i = 1,nxp
             ix = (i-1) + ib
      
      
              IF ( precip(ipgrdiam)%flt2d(ix,jy) > 0.0 ) THEN
                pltdp(i,j,k) = precip(ipgrdiam)%flt2d(ix,jy) ! already in units of mm
              ELSE
                pltdp(i,j,k) = -0.01
              ENDIF

            ENDDO
           ENDDO
          ENDDO
        
        ENDIF
      ENDIF

      IF  ( idx .eq. 311 ) THEN ! graupel mean diameter at sfc. (mm) (total)'
        IF ( ipgrnumacc < 1 ) THEN
          write(0,*) 'Cannot do graupel diameter (311)!  '
        ELSE
        
        do j = 1,nyp
          jy = (j-1) + jb
          do k = 1,nzp
            kz = (k-1) + kb
            do i = 1,nxp
             ix = (i-1) + ib
      
            IF ( precip(ipgrnumacc)%flt2d(ix,jy)*50.**2 > 10. .and. precip(ipgracc)%flt2d(ix,jy) > 0.01 ) THEN
              tmp = precip(ipgracc)%flt2d(ix,jy)/precip(ipgrnumacc)%flt2d(ix,jy)
              pltdp(i,j,k) = 1000.*(tmp*6.0/(3.14159*900.0) )**(1.0/3.0)
            ELSE
              pltdp(i,j,k) = -0.01
            ENDIF
            
            ENDDO
           ENDDO
          ENDDO
        
        ENDIF
      ENDIF

      IF  ( idx .eq. 314 ) THEN ! 'Hail mean diameter at sfc. (mm) (total)'
      IF ( iphailnumacc2 < 1 .or. iphailacc2 < 1 ) THEN
        write(0,*) 'Cannot do hail diameter2 (314)!  '
        ELSE
        
        do j = 1,nyp
          jy = (j-1) + jb
          do k = 1,nzp
            kz = (k-1) + kb
            do i = 1,nxp
             ix = (i-1) + ib
      
      
       IF ( precip(iphailnumacc2)%flt2d(ix,jy)*50.**2 > 10. .and. precip(iphailacc2)%flt2d(ix,jy) > 0.01 ) THEN
         tmp = precip(iphailacc2)%flt2d(ix,jy)/precip(iphailnumacc2)%flt2d(ix,jy)
         pltdp(i,j,k) = 1000.*(tmp*6.0/(3.14159*900.0) )**(1.0/3.0)
       ELSE
         pltdp(i,j,k) = -0.01
       ENDIF
            ENDDO
           ENDDO
          ENDDO
        
        ENDIF
      ENDIF

      
      IF ( idx .eq. 45 ) THEN ! SSI
      
      write(0,*) 'DO SSI in basp3p'

      
c      xmx = 0.0

      DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,Min(nz-1,nzp)
         kz = (k-1) + kb
          DO i = 1,nxp
          ix = (i-1) + ib

          pqs = 380.0/(pn(ix,jy,kz)+pb(kz))
          qwv = max(st(ix,jy,kz,lv),0.0)
c          xmx = Max(xmx,qwv)
          temg = tt0(ix,jy,kz)
!          temg = st(ix,jy,kz,lt)*
!     :     ((pn(ix,jy,kz)+pb(kz))/100000.0)**(287./1004.) ! temc(ix,jy,kz)

          ltemq = (temg-163.15)/fqsat+1.5
          ltemq = Min( nqsat, Max(1,ltemq) )
          qis = pqs*tabqis(ltemq)

          pltdp(i,j,k) = 100.*qwv/qis - 100.
          
         ENDDO
        ENDDO
      ENDDO

c      write(0,*) 'Max qv = ',xmx
      write(0,*) 'DONE SSI in basp3p'

      ENDIF ! idx 45
c
      IF ( idx .eq. 46 ) THEN  ! SSW


      DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,Min(nz-1,nzp)
         kz = (k-1) + kb
          DO i = 1,nxp
          ix = (i-1) + ib

          pqs = 380.0/(pn(ix,jy,kz)+pb(kz))
          qwv = max(st(ix,jy,kz,lv),0.0)
c          temg = temc(ix,jy,kz)
          temg = tt0(ix,jy,kz)
!          temg = st(ix,jy,kz,lt)*
!     :     ((pn(ix,jy,kz)+pb(kz))/100000.0)**(287./1004.) ! temc(ix,jy,kz)

          ltemq = (temg-163.15)/fqsat+1.5
          ltemq = Min( nqsat, Max(1,ltemq) )
          qvs = pqs*tabqvs(ltemq)

          pltdp(i,j,k) = 100.*qwv/qvs - 100.
         ENDDO
        ENDDO
      ENDDO

      ENDIF ! idx 46

c
      IF ( idx .eq. 49 .and. lr .gt. 1 .and. lnr .gt. 1 ) THEN ! rwdia

       DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
          ix = (i-1) + ib
              IF ( st(ix,jy,kz,lr) .gt. 1.e-6 .and. st(ix,jy,kz,lnr) .gt. 1.e-3 ) THEN
                pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lr)/(1000.*st(ix,jy,kz,lnr))
                pltdp(i,j,k) = 1.e6*(6.*pltdp(i,j,k)/(3.14159))**(1./3.)
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF
         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

      IF ( idx .eq. 50 .and. ls .gt. 1 .and. lns .gt. 1 ) THEN ! swdia

       DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
          ix = (i-1) + ib
              IF ( st(ix,jy,kz,ls) .gt. 1.e-6 .and. st(ix,jy,kz,lns) .gt. 1.e-3 ) THEN
                pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,ls)/(100.*st(ix,jy,kz,lns))
                pltdp(i,j,k) = 1.e6*(6.*pltdp(i,j,k)/(3.14159))**(1./3.)
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF
         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

      IF ( idx .eq. 51 .and. lh .gt. 1 .and. lnh .gt. 1 ) THEN ! hwdia

       DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
          ix = (i-1) + ib
              IF ( st(ix,jy,kz,lh) .gt. 1.e-6 .and. st(ix,jy,kz,lnh) .gt. 1.e-3 ) THEN
                IF ( lvh .gt. 0 ) THEN
                  hwdn = dn(ix,jy,kz)* st(ix,jy,kz,lh)/st(ix,jy,kz,lvh)
                  hwdn = Max( 170., hwdn )
                ELSE
                  hwdn = rho_qh
                ENDIF
                
                pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lh)/(hwdn*st(ix,jy,kz,lnh))
                pltdp(i,j,k) = 1.e6*(6.*pltdp(i,j,k)/(3.14159))**(1./3.)
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF
         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

      IF ( idx .eq. 52 .and. lhl .gt. 1 .and. lnhl .gt. 1 ) THEN ! hldia

       DO j = 1,nyp
        jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
          ix = (i-1) + ib
              IF ( st(ix,jy,kz,lhl) .gt. 1.e-6 .and. st(ix,jy,kz,lnhl) .gt. 1.e-3 ) THEN
                IF ( lvhl .gt. 0 ) THEN
                  hwdn = dn(ix,jy,kz)* st(ix,jy,kz,lhl)/st(ix,jy,kz,lvhl)
                  hwdn = Max( 170., hwdn )
                ELSE
                  hwdn = rho_qh
                ENDIF
                pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lhl)/(hwdn*st(ix,jy,kz,lnhl))
                pltdp(i,j,k) = 1.e6*(pltdp(i,j,k)/(3.14159))**(1./3.)
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF
         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

c
      if ( idx .eq. 261 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = prt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 262 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = prti2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 263 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = prtl2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 264 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) =  cwrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 265 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = rwrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 266 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = cirt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 267 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = swrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 268 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = glrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 269 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = gmrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 270 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ghrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 271 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = fwrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 272 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = hwrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 273 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = hlrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 274 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = cirrt2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 275 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ciprt2d(ix,jy)
      enddo
      enddo

      end if
c

      IF ( idx .eq. 570 ) THEN ! ZRW 

      DO j = 1,nyp
       jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
          DO i = 1,nxp
           ix = (i-1) + ib

           IF ( microp(1:3) == 'TAK' ) THEN
              tmp = 0.0
               DO l = 24,ntakrd
                 tmp = tmp + 1.d6*szr(l)*st(ix,jy,kz,lnr-1+l)
               ENDDO

              pltdp(i,j,k) =  Max(-10., 10*Log10(Max(0.001, tmp ) ) )
           ELSEIF ( lzr .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! predicted moment
                
              pltdp(i,j,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*( st(ix,jy,kz,lzr) ) ) ) )
            
           ELSEIF ( lnr .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! two-moment
        
           
           ENDIF
            
          ENDDO
        ENDDO
       ENDDO
        
        
       ENDIF


      IF ( idx .eq. 576 ) THEN ! ZHW

      DO j = 1,nyp
       jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
          DO i = 1,nxp
           ix = (i-1) + ib

           IF ( microp(1:3) == 'TAK' ) THEN
              tmp = 0.0
               DO l = ls250+1,ntakpd
                 tmp = tmp + 0.224*1.d6*szf(l,1)*st(ix,jy,kz,lnh-1+l)
               ENDDO

              pltdp(i,j,k) = Max(-10., 10*Log10(Max(0.001, tmp ) ) )
           ELSEIF ( lzh .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! predicted moment

              IF ( lvhl .gt. 1 ) THEN
               IF ( st(ix,jy,kz,lh) > qxmin(lh) .and. 
     &              st(ix,jy,kz,lvh) > 1.0e-28 ) THEN
               hwdn = dn(ix,jy,kz)*       
     &                st(ix,jy,kz,lh)/   
     &                st(ix,jy,kz,lvh)
               ELSE
                hwdn = rho_qh
               ENDIF
              ELSE
               hwdn = rho_qh
              ENDIF

                
              pltdp(i,j,k) =  Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(  
     &                            st(ix,jy,kz,lzh)*(hwdn/1000.)**2  ) ) ) )
            
           ELSEIF ( lnhl .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! two-moment
        
           
           ENDIF
            
          ENDDO
        ENDDO
       ENDDO
        
        
       ENDIF


      IF ( idx .eq. 577 ) THEN ! ZHL

      DO j = 1,nyp
       jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
          DO i = 1,nxp
           ix = (i-1) + ib

           IF ( microp(1:3) == 'TAK' ) THEN
              tmp = 0.0
               DO l = ls250+1,ntakpd
                 tmp = tmp + 0.224*1.d6*szf(l,2)*st(ix,jy,kz,lnhl-1+l)
               ENDDO

              pltdp(i,j,k) = Max(-10., 10*Log10(Max(0.001, tmp ) ) )
           ELSEIF ( lzhl .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! predicted moment

              IF ( lvhl .gt. 1 ) THEN
               IF ( st(ix,jy,kz,lhl) > qxmin(lhl) .and. 
     &              st(ix,jy,kz,lvhl) > 1.0e-28 ) THEN
               hwdn = dn(ix,jy,kz)*       
     &                st(ix,jy,kz,lhl)/   
     &                st(ix,jy,kz,lvhl)
               ELSE
                hwdn = rho_qhl
               ENDIF
              ELSE
               hwdn = rho_qhl
              ENDIF

                
              pltdp(i,j,k) =  Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(  
     &                            st(ix,jy,kz,lzhl)*(hwdn/1000.)**2  ) ) ) )
            
           ELSEIF ( lnhl .gt. 1 .and. microp(1:1) == 'Z' ) THEN ! two-moment
        
           
           ENDIF
            
          ENDDO
        ENDDO
       ENDDO
        
        
       ENDIF







      if ( idx .eq. 581 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ptot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 582 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ptoti2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 583 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ptotl2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 584 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = cwtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 585 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = rwtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 586 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = citot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 587 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + i
      pltdp(i,j,k) = swtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 588 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = gltot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 589 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = gmtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 590 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ghtot2d(ix,jy) 
      enddo
      enddo

      end if
c
      if ( idx .eq. 591 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = fwtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 592 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = hwtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 593 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = hltot2d(ix,jy)
      enddo
      enddo

      end if
c
      k=1
      if ( idx .eq. 594 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = cirtot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 595 ) then
      k=1
      do j = 1,nyp
      jy = (j-1) + jb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ciptot2d(ix,jy)
      enddo
      enddo

      end if
c
      if ( idx .eq. 503 ) then
      DO j = 1,nyp
       jy = (j-1) + jb
        DO k = 1,nzp
         kz = (k-1) + kb
          DO i = 1,nxp
           ix = (i-1) + ib
           pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.*dn(ix,jy,kz),0.0)
          ENDDO
        ENDDO
      ENDDO
      end if

      if ( idx .eq. 504 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,lc)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,lc)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.,0.0)
           ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 505 ) then ! qr
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,lr)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,lr)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,lr)*1000.,0.0)
           ENDIF
      enddo
      enddo
      enddo

      end if


c
      if ( idx .eq. 506 ) then ! cloud ice
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,li)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,li)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
           ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 507 ) then ! snow
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,ls)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,ls)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,ls)*1000.,0.0)
           ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 508 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lgl)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 509 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lgm)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 510 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lgh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 511 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lf)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 512 ) then ! graupel
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,lh)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,lh)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
           ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 513 ) then ! hail
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        IF ( lhl > 1 ) THEN
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ix,ia,kz,lhl)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + st(ia,jy,kz,lhl)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e3
           ELSE
            pltdp(i,j,k) = max(st(ix,jy,kz,lhl)*1000.,0.0)
           ENDIF
         ELSE
           pltdp(i,j,k) = 0
         ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 514 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( lis >= 1 ) THEN
      pltdp(i,j,k) = max(st(ix,jy,kz,lis)*1000.,0.0)
      ELSE
      pltdp(i,j,k) = max(st(ix,jy,kz,lir)*1000.,0.0)
      ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 515 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lip)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 516 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lip)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 517 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 518 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 519 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 520 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 521 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 522 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 523 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 524 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 526 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 528 ) then
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
     >             + max(st(ix,jy,kz,lip)*1000.,0.0)
      ENDDO
      ENDDO
      ENDDO
      end if

      if ( idx .eq. 529 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
     >             + max(st(ix,jy,kz,lir)*1000.,0.0)
     >             + max(st(ix,jy,kz,lip)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 530 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lgl)*1000.,0.0)
     >             + max(st(ix,jy,kz,lgm)*1000.,0.0)
     >             + max(st(ix,jy,kz,lgh)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 531 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.,0.0)
     >             + max(st(ix,jy,kz,li)*1000.,0.0)
     >             + max(st(ix,jy,kz,lir)*1000.,0.0)
     >             + max(st(ix,jy,kz,lip)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 532 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.,0.0)
     >             + max(st(ix,jy,kz,li)*1000.,0.0)
     >             + max(st(ix,jy,kz,lir)*1000.,0.0)
     >             + max(st(ix,jy,kz,lip)*1000.,0.0)
     >             + max(st(ix,jy,kz,ls)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 533 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lc)*1000.,0.0)
     >             + max(st(ix,jy,kz,lr)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 534 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,li)*1000.,0.0)
     >             + max(st(ix,jy,kz,lir)*1000.,0.0)
     >             + max(st(ix,jy,kz,lip)*1000.,0.0)
     >             + max(st(ix,jy,kz,ls)*1000.,0.0)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 535 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0

      IF ( ls > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,ls)*1000.,0.0)

      IF ( lgl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgl)*1000.,0.0)
      IF ( lgm > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgm)*1000.,0.0)
      IF ( lgh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgh)*1000.,0.0)

      IF ( lf > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lf)*1000.,0.0)
      IF ( lh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lh)*1000.,0.0)
      IF ( lhl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lhl)*1000.,0.0)

      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 536 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

      pltdp(i,j,k) = 0.0

      IF ( ls > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,ls)*1000.,0.0)

      IF ( li > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,li)*1000.,0.0)
      IF ( lir > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lir)*1000.,0.0)
      IF ( lip > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lip)*1000.,0.0)

      IF ( lgl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgl)*1000.,0.0)
      IF ( lgm > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgm)*1000.,0.0)
      IF ( lgh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgh)*1000.,0.0)

      IF ( lf > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lf)*1000.,0.0)
      IF ( lh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lh)*1000.,0.0)
      IF ( lhl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lhl)*1000.,0.0)

      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 537 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

      pltdp(i,j,k) = 0.0

      IF ( lc > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lc)*1000.,0.0)
      IF ( lr > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lr)*1000.,0.0)

      IF ( ls > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,ls)*1000.,0.0)

      IF ( li > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,li)*1000.,0.0)
      IF ( lir > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lir)*1000.,0.0)
      IF ( lip > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lip)*1000.,0.0)

      IF ( lgl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgl)*1000.,0.0)
      IF ( lgm > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgm)*1000.,0.0)
      IF ( lgh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgh)*1000.,0.0)

      IF ( lf > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lf)*1000.,0.0)
      IF ( lh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lh)*1000.,0.0)
      IF ( lhl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lhl)*1000.,0.0)

      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 538 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

      pltdp(i,j,k) = 0.0

      IF ( lr > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lr)*1000.,0.0)

      IF ( ls > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,ls)*1000.,0.0)


      IF ( lgl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgl)*1000.,0.0)
      IF ( lgm > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgm)*1000.,0.0)
      IF ( lgh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgh)*1000.,0.0)

      IF ( lf > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lf)*1000.,0.0)
      IF ( lh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lh)*1000.,0.0)
      IF ( lhl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lhl)*1000.,0.0)

      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 539 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

      pltdp(i,j,k) = 0.0

      IF ( lgl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgl)*1000.,0.0)
      IF ( lgm > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgm)*1000.,0.0)
      IF ( lgh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lgh)*1000.,0.0)

      IF ( lf > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lf)*1000.,0.0)
      IF ( lh > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lh)*1000.,0.0)
      IF ( lhl > 1 ) pltdp(i,j,k) = pltdp(i,j,k) + max(st(ix,jy,kz,lhl)*1000.,0.0)

      enddo
      enddo
      enddo

      end if
c
c
      if ( idx .eq. 540 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1000.,0.0)
     >             + max(st(ix,jy,kz,lhl)*1000.,0.0)
      enddo
      enddo
      enddo

      end if

c graupel particle density
      if ( idx .eq. 541 .or. idx == 547 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( idx == 547 .and. do_dualpol ) THEN
        pltdp(i,j,k) = 0.0
        IF ( st(ix,jy,kz,lh) .gt. 1.e-7 ) THEN
           pltdp(i,j,k) = Min( 1000., Max( 50., rhomg(ix,jy,kz) ) )
        ENDIF
      ELSE
       pltdp(i,j,k) = 0.0
        IF ( st(ix,jy,kz,lh) .gt. 1.e-7 ) THEN
         IF ( lvh > 0 ) THEN
          IF ( st(ix,jy,kz,lvh) .gt. 1.e-20  ) THEN
            pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lh)/st(ix,jy,kz,lvh)
            pltdp(i,j,k) = Min( 1000., Max( 50., pltdp(i,j,k) ) )
          ENDIF
         ELSE
           pltdp(i,j,k) = rho_qh
         ENDIF
        ENDIF
      
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

c frozen drop particle density
      if ( idx .eq. 545 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( st(ix,jy,kz,lvf) .gt. 1.e-20 .and. st(ix,jy,kz,lf) .gt. 1.e-5 ) THEN
        pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lf)/st(ix,jy,kz,lvf)
        pltdp(i,j,k) = Min( 1000., Max( 50., pltdp(i,j,k) ) )
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

c hail particle density
      if ( idx .eq. 546 .or. idx == 548 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        pltdp(i,j,k) = 0.0
      IF ( idx == 548 .and. do_dualpol ) THEN
        IF ( st(ix,jy,kz,lvhl) .gt. 1.e-20 .and. st(ix,jy,kz,lhl) .gt. 1.e-5 ) THEN
          pltdp(i,j,k) =  Max(50., Min( 1000., rhomh(ix,jy,kz) ) )
        ENDIF
      ELSE
      IF ( st(ix,jy,kz,lvhl) .gt. 1.e-20 .and. st(ix,jy,kz,lhl) .gt. 1.e-5 ) THEN
        pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lhl)/st(ix,jy,kz,lvhl)
        pltdp(i,j,k) = Min( 1000., Max( 50., pltdp(i,j,k) ) )
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

c Liquid water on snow (QSW)
      if ( idx .eq. 542 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( lsw .gt. 1 ) THEN
          pltdp(i,j,k) = 1000.*Max(0.0, st(ix,jy,kz,lsw) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

c Liquid water on graupel (QHW)
      if ( Abs(idx) .eq. 543 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( do_dualpol .and. idx == 543 ) THEN
          pltdp(i,j,k) = 1000.*Max(0.0, qgw(ix,jy,kz) )
      ELSEIF ( lhw .gt. 1 ) THEN
          pltdp(i,j,k) = 1000.*Max(0.0, st(ix,jy,kz,lhw) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

c Liquid water on Hail (QHLW)
      if ( idx .eq. 544 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( lhlw .gt. 1 ) THEN
          pltdp(i,j,k) = 1000.*Max(0.0, st(ix,jy,kz,lhlw) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

! dbz (simulator)
      if ( idx .eq. 601 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( do_dualpol ) THEN
          pltdp(i,j,k) = Max(-1.0, logz(ix,jy,kz) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

! Zdr (simulator)
      if ( idx .eq. 603 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( do_dualpol ) THEN
          pltdp(i,j,k) = Max(-1.,logzdr(ix,jy,kz) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

! rhohv (simulator)
      if ( idx .eq. 607 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( do_dualpol ) THEN
          pltdp(i,j,k) = Max(-1.,rhv(ix,jy,kz) )
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if

! Zdr-rain (simulator)
      if ( idx .eq. 611 ) then  
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 0.0
      IF ( do_dualpol ) THEN
          pltdp(i,j,k) = logZdrrain(ix,jy,kz) 
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      end if


! #########################################################################################
! #########################################################################################
c hail shape parameter
      if ( idx == 711  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
       IF (lhl > 1) THEN
      hwdn = rho_qhl
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      IF ( lzhl < 1 ) THEN
        pltdp(i,j,k) = alphahl
        CYCLE
      ENDIF

      IF ( lnhl < 1 ) THEN
        pltdp(i,j,k) = alphahl
        CYCLE
      ENDIF
      
      IF ( st(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. st(ix,jy,kz,lnhl) .gt. 1.e-19 .and.  
     &         10.*log10(Max(1.e-30,1.e18*st(ix,jy,kz,lzhl))) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnhl)
                qr = st(ix,jy,kz,lhl)
                z = st(ix,jy,kz,lzhl)
                alpha1 = 0

                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO

               pltdp(i,j,k) = alp
                
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
       ENDIF ! lhl > 1
      
      ELSE ! takhashi bin
      
      hwdn = rho_qhl_tak
      
      write(0,*) 'BASP3P: lf75,szf(40,2),rhof',lf75,szf(40,2),rhof(2),rho_qhl_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

               IF (  st(ix,jy,kz,lhl) .gt. qxmin(lhl) ) THEN

          ! get Z
              z = 0.0
               DO l = lf75+1,lfmax
                 z = z + szf(l,2)*st(ix,jy,kz,lnhl-1+l)
               ENDDO
              z = z*1.0e-12/rhof(2)**2   ! divide by rhof because it is contained in szf
          ! get concentration
              chw = 0.0
               DO l = lf75+1,lfmax
                 chw = chw + st(ix,jy,kz,lnhl-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(ix,jy,kz,lhl)
               hwdn = rho_qhl_tak
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-19 .and. 10.*log10(Max(1.e-30,1.e18*z)) > 0.0 ) THEN

            rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
                 
                pltdp(i,j,k) = alp

               ELSE
                pltdp(i,j,k) = -1.
               ENDIF
               
               ELSE
                pltdp(i,j,k) = -1.
               ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if ! shape parameter

! #########################################################################################
! #########################################################################################
c hail content
      if ( idx == 712  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qhl
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( st(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. st(ix,jy,kz,lnhl) .gt. 1.e-19 ) THEN
        pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lhl)
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qhl_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

                qr = 0.0

               DO n = ls250+1,lfmax
                qr = qr + st(ix,jy,kz,lnhl+n-1)*sxf(n,2)  ! total mass 
               ENDDO
              IF ( 1.e6*chw > 1.e-3 ) THEN
              
                  qr = 1.e6*qr
                  pltdp(i,j,k) = qr
                
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if


! #########################################################################################
! #########################################################################################
c hail number
      if ( idx == 713  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qhl
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( st(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. st(ix,jy,kz,lnhl) .gt. 1.e-19 ) THEN
        pltdp(i,j,k) = st(ix,jy,kz,lnhl)
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qhl_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

                chw = 0.0

               DO n = ls250+1,lfmax
                chw = chw + st(ix,jy,kz,lnhl+n-1)  ! total number 
               ENDDO
              IF ( 1.e6*chw > 1.e-3 ) THEN
              
                  chw = 1.e6*chw
                  pltdp(i,j,k) = chw
                
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if

! #########################################################################################
! #########################################################################################

c hail mean volume diameter, max mass diameter, mass-weighted diam
      if ( idx .eq. 714 .or. idx == 716 .or. idx == 717 ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qhl
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        
        IF ( st(ix,jy,kz,lhl) <= qxmin(lhl) ) THEN
         pltdp(i,j,k) = 0.0
         CYCLE
        ENDIF
        
        IF ( lvhl .gt. 1 ) THEN
          hwdn = dn(ix,jy,kz)*st(ix,jy,kz,lhl)/st(ix,jy,kz,lvhl)
        ENDIF
          
        IF ( idx == 714 ) THEN
          IF ( st(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. st(ix,jy,kz,lnhl) .gt. cxmin ) THEN
            pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lhl)/(hwdn*st(ix,jy,kz,lnhl))
            pltdp(i,j,k) = 1.e6*(6.0*pltdp(i,j,k)/(3.14159))**(1./3.)
          ELSE
            pltdp(i,j,k) = 0.0
          ENDIF
        
        ELSEIF ( idx == 716 .or. idx == 717 ) THEN ! need shape parameter for other diameters
 
         IF ( st(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. st(ix,jy,kz,lnhl) .gt. cxmin ) THEN
         IF ( lzhl < 1 .or. lnhl < 1 ) THEN
           alp = alphahl
         ENDIF
      
         IF ( lzhl > 1 ) THEN ! find shape parameter
            IF ( st(ix,jy,kz,lzhl) > zxmin ) THEN

                chw = st(ix,jy,kz,lnhl)
                qr = st(ix,jy,kz,lhl)
                z = st(ix,jy,kz,lzhl)
                alpha1 = 0
                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
             ENDIF
           ENDIF
           
           ! now have shape parameter
                tmp = dn(ix,jy,kz)* st(ix,jy,kz,lhl)/(hwdn*st(ix,jy,kz,lnhl))
                diam = (6.0*tmp/(3.14159))**(1./3.)
            
           !  alp = alphahl ! test!
            
           IF ( idx == 716 ) THEN ! max mass diam
                  pltdp(i,j,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           
           ELSEIF ( idx == 717 ) THEN ! mass-weighted
                  pltdp(i,j,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           ENDIF
           
           ELSE 
             pltdp(i,j,k) = 0.0
           ENDIF
          
          ENDIF ! idx
        
        
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qhl_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
              
              tmp = 0.0
              tmpn = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0

               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                qr  = qr  + st(ix,jy,kz,lnhl+n-1)*sxf(n,2) ! total mass
                chw = chw + st(ix,jy,kz,lnhl+n-1)  ! total number 
                tmp = tmp + st(ix,jy,kz,lnhl+n-1)*sxf(n,2)*srf(n,2)*2. ! D*N(D)*m(D)
                area = area + st(ix,jy,kz,lnhl+n-1)*srf(n,2)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(ix,jy,kz,lnhl+n-1)*2.*srf(n,2)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO
!               qr = 1.e6*1.e-3*tmp ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
              
!              chw = 1.e6*chw
              IF (  1.e6*1.e-3*qr > qxmin(lhl) .and. 1.e6*chw > cxmin ) THEN
              
                IF ( idx == 714 ) THEN ! mean vol. diam
                   chw = 1.e6*chw
                   qr = 1.e6*1.e-3*qr ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
                  pltdp(i,j,k) = 1.e6*(6.*qr/(1000.*3.14159*chw))**(1./3.)
                ELSEIF ( idx == 716 ) THEN ! max mass diam
                   pltdp(i,j,k) =  1.e6*0.01*areaw/area  ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ELSEIF ( idx == 717 ) THEN ! mass-weighted diam
                   pltdp(i,j,k) =   0.01*1.e6*tmp/qr 
                ENDIF
                
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if


! #########################################################################################
! #########################################################################################

c frozen drop shape parameter
      if ( idx == 721  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qh
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      IF ( lzh < 1 ) THEN
        pltdp(i,j,k) = alphah
        CYCLE
      ENDIF

      IF ( lnh < 1 ) THEN
        pltdp(i,j,k) = alphah
        CYCLE
      ENDIF
      
      IF ( st(ix,jy,kz,lf) .gt. qxmin(lf) .and. st(ix,jy,kz,lnf) .gt. 1.e-19 .and. st(ix,jy,kz,lzf) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnf)
                qr = st(ix,jy,kz,lf)
                z = st(ix,jy,kz,lzf)
                alpha1 = 0

                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO

               pltdp(i,j,k) = alp
                
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      
      ENDIF
      
      end if ! shape parameter

! #########################################################################################
! #########################################################################################
c
c frozen drop content
      if ( idx .eq. 722  ) then

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        pltdp(i,j,k) = 1.e3*dn(ix,jy,kz)*st(ix,jy,kz,lf)
      ENDDO
      ENDDO
      ENDDO

       end if
! #########################################################################################
! #########################################################################################
c frozen drop number
      if ( idx == 723  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qh
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( st(ix,jy,kz,lf) .gt. qxmin(lf) .and. st(ix,jy,kz,lnf) .gt. 1.e-19 ) THEN
        pltdp(i,j,k) = st(ix,jy,kz,lnf)
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      
      ENDIF
      
      end if


! #########################################################################################
! #########################################################################################
c
c frozen drop mean volume diameter
      if ( idx .eq. 724 .or. idx == 726 .or. idx == 727 ) then

      IF ( microp(1:3) /= 'TAK' ) THEN
      
        hwdn = rho_qh
      
        do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
        kz = (k-1) + kb
        do i = 1,nxp
        ix = (i-1) + ib

        IF ( st(ix,jy,kz,lf) <= qxmin(lf) ) THEN
         pltdp(i,j,k) = 0.0
         CYCLE
        ENDIF
        
        IF ( lvf .gt. 1 ) THEN
          hwdn = dn(ix,jy,kz)*st(ix,jy,kz,lf)/st(ix,jy,kz,lvf)
        ENDIF
          
        IF ( idx == 724 ) THEN
          IF ( st(ix,jy,kz,lf) .gt. qxmin(lf) .and. st(ix,jy,kz,lnf) .gt. 0.0 ) THEN
            pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lf)/(hwdn*st(ix,jy,kz,lnf))
            pltdp(i,j,k) = 1.e6*(6.0*pltdp(i,j,k)/(3.14159))**(1./3.)
          ELSE
            pltdp(i,j,k) = 0.0
          ENDIF
        
        ELSEIF ( idx == 726 .or. idx == 727 ) THEN ! need shape parameter for other diameters
 
         IF ( lzf < 1 .or. lnf < 1 ) THEN
           alp = alphah
         ENDIF
      
         IF ( lzf > 1 ) THEN ! find shape parameter
           IF ( st(ix,jy,kz,lf) .gt. qxmin(lf) .and. st(ix,jy,kz,lnf) .gt. 0.0 .and. st(ix,jy,kz,lzf) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnf)
                qr = st(ix,jy,kz,lf)
                z = st(ix,jy,kz,lzf)
                alpha1 = 0
                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
        
             ENDIF
           ENDIF
           
           ! now have shape parameter
                tmp = dn(ix,jy,kz)* st(ix,jy,kz,lf)/(hwdn*st(ix,jy,kz,lnf))
                diam = (6.0*tmp/(3.14159))**(1./3.)
            
           !  alp = alphah ! test!
            
           IF ( idx == 726 ) THEN ! max mass diam
                  pltdp(i,j,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           
           ELSEIF ( idx == 727 ) THEN ! mass-weighted
                  pltdp(i,j,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           ENDIF
           
           
          ENDIF ! idx


        ENDDO
        ENDDO
        ENDDO


      ELSE ! takhashi bin
      

      
      ENDIF

      end if

! #########################################################################################


c graupel shape parameter
      if ( idx == 731  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qh
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      IF ( lzh < 1 ) THEN
        pltdp(i,j,k) = alphah
        CYCLE
      ENDIF

      IF ( lnh < 1 ) THEN
        pltdp(i,j,k) = alphah
        CYCLE
      ENDIF
      
      IF ( st(ix,jy,kz,lh) .gt. qxmin(lh) .and. st(ix,jy,kz,lnh) .gt. 1.e-19 .and. st(ix,jy,kz,lzh) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnh)
                qr = st(ix,jy,kz,lh)
                z = st(ix,jy,kz,lzh)
                alpha1 = 0

                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO

               pltdp(i,j,k) = alp
                
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qh_tak
      
    !  write(0,*) 'BASP3P: lf75,szf(40,1),rhof',lf75,szf(40,1),rhof(1),rho_qh_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

               IF (  st(ix,jy,kz,lh) .gt. qxmin(lh) ) THEN

          ! get Z
              z = 0.0
               DO l = lf75+1,lfmax
                 z = z + szf(l,1)*st(ix,jy,kz,lnh-1+l)
               ENDDO
              z = z*1.0e-12/rhof(1)**2   ! divide by rhof because it is contained in szf
          ! get concentration
              chw = 0.0
               DO l = lf75+1,lfmax
                 chw = chw + st(ix,jy,kz,lnh-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(ix,jy,kz,lh)
               hwdn = rho_qh_tak
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-9 .and.  z .gt. 0.0 ) THEN

            rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
                 
                pltdp(i,j,k) = alp

               ELSE
                pltdp(i,j,k) = -1.
               ENDIF
               
               ELSE
                pltdp(i,j,k) = -1.
               ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if ! shape parameter

! #########################################################################################
! #########################################################################################
c
c graupel content
      if ( idx .eq. 732  ) then

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        pltdp(i,j,k) = 1.e3*dn(ix,jy,kz)*st(ix,jy,kz,lh)
      ENDDO
      ENDDO
      ENDDO

       end if
! #########################################################################################
! #########################################################################################
c graupel number
      if ( idx == 733  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qh
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( st(ix,jy,kz,lh) .gt. qxmin(lh) .and. st(ix,jy,kz,lnh) .gt. 1.e-19 ) THEN
        pltdp(i,j,k) = st(ix,jy,kz,lnh)
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qh_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

                chw = 0.0

               DO n = ls250+1,lfmax
                chw = chw + st(ix,jy,kz,lnh+n-1)  ! total number 
               ENDDO
              IF ( 1.e6*chw > 1.e-3 ) THEN
              
                  chw = 1.e6*chw
                  pltdp(i,j,k) = chw
                
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if


! #########################################################################################
! #########################################################################################
c
c graupel mean volume diameter
      if ( idx .eq. 734 .or. idx == 736 .or. idx == 737 ) then

      IF ( microp(1:3) /= 'TAK' ) THEN
      
        hwdn = rho_qh
      
        do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
        kz = (k-1) + kb
        do i = 1,nxp
        ix = (i-1) + ib

        IF ( st(ix,jy,kz,lh) <= qxmin(lh) ) THEN
         pltdp(i,j,k) = 0.0
         CYCLE
        ENDIF
        
        IF ( lvh .gt. 1 ) THEN
          hwdn = dn(ix,jy,kz)*st(ix,jy,kz,lh)/st(ix,jy,kz,lvh)
        ENDIF
          
        IF ( idx == 734 ) THEN
          IF ( st(ix,jy,kz,lh) .gt. qxmin(lh) .and. st(ix,jy,kz,lnh) .gt. 0.0 ) THEN
            pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lh)/(hwdn*st(ix,jy,kz,lnh))
            pltdp(i,j,k) = 1.e6*(6.0*pltdp(i,j,k)/(3.14159))**(1./3.)
          ELSE
            pltdp(i,j,k) = 0.0
          ENDIF
        
        ELSEIF ( idx == 736 .or. idx == 737 ) THEN ! need shape parameter for other diameters
 
         IF ( lzh < 1 .or. lnh < 1 ) THEN
           alp = alphah
         ENDIF
      
         IF ( lzh > 1 ) THEN ! find shape parameter
           IF ( st(ix,jy,kz,lh) .gt. qxmin(lh) .and. st(ix,jy,kz,lnh) .gt. 0.0 .and. st(ix,jy,kz,lzh) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnh)
                qr = st(ix,jy,kz,lh)
                z = st(ix,jy,kz,lzh)
                alpha1 = 0
                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
        
             ENDIF
           ENDIF
           
           ! now have shape parameter
                tmp = dn(ix,jy,kz)* st(ix,jy,kz,lh)/(hwdn*st(ix,jy,kz,lnh))
                diam = (6.0*tmp/(3.14159))**(1./3.)
            
           !  alp = alphah ! test!
            
           IF ( idx == 736 ) THEN ! max mass diam
                  pltdp(i,j,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           
           ELSEIF ( idx == 737 ) THEN ! mass-weighted
                  pltdp(i,j,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           ENDIF
           
           
          ENDIF ! idx


        ENDDO
        ENDDO
        ENDDO


      ELSE ! takhashi bin
      
      hwdn = rho_qh_tak

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

              tmp = 0.0
              tmpn = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0

               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                qr  = qr  + st(ix,jy,kz,lnh+n-1)*sxf(n,1) ! total mass
                chw = chw + st(ix,jy,kz,lnh+n-1)  ! total number 
                tmp = tmp + st(ix,jy,kz,lnh+n-1)*sxf(n,1)*srf(n,1)*2. ! D*N(D)*m(D)
                area = area + st(ix,jy,kz,lnh+n-1)*srf(n,1)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(ix,jy,kz,lnh+n-1)*2.*srf(n,1)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO
              
              
              IF (  1.e6*1.e-3*qr > qxmin(lh) .and. 1.e6*chw > 1.e-19 ) THEN
              
                IF ( idx == 734 ) THEN ! mean vol. diam
                   chw = 1.e6*chw
                   qr = 1.e6*1.e-3*qr ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
                  pltdp(i,j,k) = 1.e6*(6.*qr/(1000.*3.14159*chw))**(1./3.)
                ELSEIF ( idx == 736 ) THEN ! max mass diam
                   pltdp(i,j,k) =  1.e6*0.01*areaw/area  ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ELSEIF ( idx == 737 ) THEN ! mass-weighted diam
                   pltdp(i,j,k) =   0.01*1.e6*tmp/qr 
                ENDIF
                
              ELSE
                pltdp(i,j,k) = 0.0
              ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF

      end if

! #########################################################################################
! #########################################################################################

c rain shape parameter
      if ( idx == 751  ) then
      
      IF ( microp(1:3) /= 'TAK' ) THEN
      hwdn = rho_qr
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      IF ( lzr < 1 ) THEN
        pltdp(i,j,k) = alphar
        CYCLE
      ENDIF

      IF ( lnr < 1 ) THEN
        pltdp(i,j,k) = alphar
        CYCLE
      ENDIF
      
      IF ( st(ix,jy,kz,lr) .gt. qxmin(lr) .and. st(ix,jy,kz,lnr) .gt. 1.e-19 .and. st(ix,jy,kz,lzr) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnr)
                qr = st(ix,jy,kz,lr)
                z = st(ix,jy,kz,lzr)
                alpha1 = 0

                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO

               pltdp(i,j,k) = alp
                
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE ! takhashi bin
      
      hwdn = rho_qr
      

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

           IF (  st(ix,jy,kz,lr) .gt. qxmin(lr) ) THEN

          ! get Z
              z = 0.0
               DO l = lr100,lmax
                 z = z + szr(l)*st(ix,jy,kz,lnr-1+l)
               ENDDO
              z = z*1.0e-12
          ! get concentration
              chw = 0.0
               DO l = lr100,lmax
                 chw = chw + st(ix,jy,kz,lnr-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(ix,jy,kz,lr)
               hwdn = rho_qr
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-9 .and.  z .gt. 0.0 ) THEN

            rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
                 
                pltdp(i,j,k) = alp

               ELSE
                pltdp(i,j,k) = -1.
               ENDIF
               
               ELSE
                pltdp(i,j,k) = -1.
               ENDIF

      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      
      end if ! shape parameter

! #########################################################################################
! #########################################################################################
c
c rain content
      if ( idx .eq. 752 ) then
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
        pltdp(i,j,k) = 1.e3*dn(ix,jy,kz)*st(ix,jy,kz,lr)
      ENDDO
      ENDDO
      ENDDO
      end if
! #########################################################################################
! #########################################################################################
c
c rain mean volume diameter
      if ( idx .eq. 754 .or. idx == 756 .or. idx == 757 .or. idx == 758 ) then
      hwdn = rho_qr
      
      write(0,*) 'BASP3P: idx = ',idx,', Micro = ',microp
      
      IF ( microp(1:3) /= 'TAK' ) THEN

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
          hwdn = 1000.

        IF ( st(ix,jy,kz,lr) <= qxmin(lr) ) THEN
         pltdp(i,j,k) = 0.0
         CYCLE
        ENDIF
        
          
        IF ( idx == 754 ) THEN
          IF ( st(ix,jy,kz,lr) .gt. qxmin(lr) .and. st(ix,jy,kz,lnr) .gt. 0 ) THEN
            pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lr)/(hwdn*st(ix,jy,kz,lnr))
            pltdp(i,j,k) = 1.e6*(6.0*pltdp(i,j,k)/(3.14159))**(1./3.)
          ELSE
            pltdp(i,j,k) = 0.0
          ENDIF
        
        ELSEIF ( idx == 756 .or. idx == 757 .or. idx == 758 ) THEN ! need shape parameter for other diameters
 
         IF ( lzr < 1 .or. lnr < 1 ) THEN
           alp = alphar
         ENDIF
      
         IF ( lzr > 1 ) THEN ! find shape parameter
           IF ( st(ix,jy,kz,lr) .gt. qxmin(lr) .and. st(ix,jy,kz,lnr) .gt. 0 .and. st(ix,jy,kz,lzr) > 0.0 ) THEN

                chw = st(ix,jy,kz,lnr)
                qr = st(ix,jy,kz,lr)
                z = st(ix,jy,kz,lzr)
                alpha1 = 0
                
                rdi = z*(pi/6.*hwdn)**2*chw/((dn(ix,jy,kz)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
        
             ENDIF
           ENDIF
           
           ! now have shape parameter
                tmp = dn(ix,jy,kz)* st(ix,jy,kz,lr)/(hwdn*st(ix,jy,kz,lnr)) ! mean volume, xv0
                diam = (6.0*tmp/(3.14159))**(1./3.) ! mean volume diameter
            
            
           IF ( idx == 756 ) THEN ! max mass diam
                  pltdp(i,j,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           
           ELSEIF ( idx == 757 ) THEN ! mass-weighted
                  pltdp(i,j,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           ELSEIF ( idx == 758 ) THEN ! median volume diameter
                  pltdp(i,j,k) =  1000.*(3.67+alp)*diam*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)! units of mm (Ulbrich 1983, JCAM) ! 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
           ENDIF
           
           
          ENDIF ! idx

      
      ENDDO
      ENDDO
      ENDDO


      ELSE ! Takahashi
 
       hwdn = rho_qr
 
       write(0,*) 'BASP3P: lr100 = ',lr100,', lnr = ',lnr, ', sx10 = ',sx(10),',lmax = ',lmax

      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib

              tmp = 0.0
              tmpn = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0

               DO n = lr100,lmax
!               DO n = 1,lfmax
                qr  = qr  + st(ix,jy,kz,lnr+n-1)*sx(n) ! total mass
                masscum(n) = qr ! cumulative mass
                chw = chw + st(ix,jy,kz,lnr+n-1)  ! total number 
                tmp = tmp + st(ix,jy,kz,lnr+n-1)*sx(n)*sr(n)*2. ! D*N(D)*m(D)
                area = area + st(ix,jy,kz,lnr+n-1)*sr(n)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(ix,jy,kz,lnr+n-1)*2.*sr(n)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO

!               qr = 1.e6*1.e-3*qr ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
              
!              write(0,*) 'i,j,qr,chw = ',ix,jy,qr,chw
              
!              chw = 1.e6*chw
              
               IF (  1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
              
                IF ( idx == 754 ) THEN ! mean vol. diam
                   chw = 1.e6*chw
                   qr = 1.e6*1.e-3*qr ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
                  pltdp(i,j,k) = 1.e6*(6.*qr/(1000.*3.14159*chw))**(1./3.)
                ELSEIF ( idx == 756 ) THEN ! max mass diam
                   pltdp(i,j,k) =  1.e6*0.01*areaw/area  ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ELSEIF ( idx == 757 ) THEN ! mass-weighted diam
                   pltdp(i,j,k) =   0.01*1.e6*tmp/qr 
                ELSEIF ( idx == 758 ) THEN ! median volume diameter
                
                n = 0
          
                DO il = 1,ntakrd
                  IF ( qr*0.5d0 < masscum(il) ) THEN
                    n = il
                    EXIT
                  ENDIF
                ENDDO

                IF ( n == 1 ) THEN
                   pltdp(i,j,k)  = sr(n)
                ELSE
                ! n is index of masscum just above d0
                   tmp =  sr(n-1) + (0.5*qr - masscum(n-1))/(masscum(n) - masscum(n-1))*(sr(n) - sr(n-1)) ! radius
                   pltdp(i,j,k) = 2.*tmp*10. ! *2 for diam, *10 for units of mm
                ENDIF

                ENDIF
                
              ELSE
                pltdp(i,j,k) = 0.
              ENDIF

      ENDDO
      ENDDO
      ENDDO

      
      ENDIF ! ( microp(1:3) /= 'TAK' ) 
      
      end if

! #########################################################################################
! #########################################################################################


      IF ( radarmask .eq. 1 .and. ( idx >= 711 .and. idx <= 759  ) ) THEN
        do j = 1,nyp
        jy = (j-1) + jb
          do k = 1,nzp
          kz = (k-1) + kb
            do i = 1,nxp
            ix = (i-1) + ib
             
              IF ( dbz%flt3d(ix,jy,kz) .lt. lowdbz ) THEN
                pltdp(i,j,k) = 0.0
              ENDIF
         
            ENDDO
          ENDDO
        ENDDO
      
      ENDIF

! #########################################################################################
! #########################################################################################
c
c cloud droplet mean volume diameter
      if ( idx .eq. 47 ) then
      hwdn = rho_qh
      
      IF ( lnc .gt. 1 ) THEN
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
          hwdn = 1000.
      IF ( st(ix,jy,kz,lc) .gt. 1.e-6 .and. st(ix,jy,kz,lnc) .gt. 1.e-3 ) THEN
        pltdp(i,j,k) = dn(ix,jy,kz)*st(ix,jy,kz,lc)/(hwdn*st(ix,jy,kz,lnc))
        pltdp(i,j,k) = 1.e6*(6.0*pltdp(i,j,k)/(3.14159))**(1./3.)
      ELSE
        pltdp(i,j,k) = 0.0
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      ELSE
       write(*,*) 'Cannot calc droplet diameter'
       pltdp(:,:,:) = 0.0
      ENDIF
      
      end if

      if ( na .gt. 8 ) then
      if ( idx .eq. 140 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,9)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 141 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,10)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 142 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,11)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 143 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,12)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 144 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,13)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 145 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,14)*1000., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 239 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,15)*1., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 240 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,16)*1., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 241 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lgm)*1., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 242 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lh)*1., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 243 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = max(st(ix,jy,kz,lhl)*1., 0.0)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 146 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,1)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 147 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,2)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 148 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,3)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 149 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,4)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 150 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,5)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 151 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,6)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( na .gt. 8 ) then
      if ( idx .eq. 152 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = hdnx(ix,jy,kz,7)
      enddo
      enddo
      enddo

      end if
      end if
c
      if ( idx .eq. 180 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      if(idx.eq.180) pltdp(i,j,k) = st(ix,jy,kz,2)
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 181 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,3)*1000.
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 182 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,4)*100.
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 183 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,5)
      enddo
      enddo
      enddo

      end if
c
c      if ( idx .eq. 332 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,1)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 333 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = v2(ix,jy,kz,1)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 334 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = v2(ix,jy,kz,2)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 335 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = v2(ix,jy,kz,3)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 336 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = p2(ix,jy,kz)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 337 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,02)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 338 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,03)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 339 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,04)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 340 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,05)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 341 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,06)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 342 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,07)
c      enddo
c      enddo
c      enddo

c      end if 
cc
c      if ( idx .eq. 343 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = a2(ix,jy,kz,08)
c      enddo
c      enddo
c      enddo

c      end if 
c
c  RADAR STUFF
c
c      if ( idx .ge. 600 .and. idx .le. 800 ) then
cc
c      write(0,*) 'BASP3P: radar stuff;  idx = ',idx
c      read(22,*) ityp,filrad
c      if(ityp.eq.1) open(44,file=filrad,form='formatted')
cc
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      if(ityp.eq.1) read(44,*) pltdp(i,j,k)
c      enddo
c      enddo
c      enddo

 
cc     close(44)
cc
c      end if
c
c      if ( idx .eq. 320 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = algt(ix,jy,kz)
c      enddo
c      enddo
c      enddo

c      end if
c
c      if ( idx .eq. 321 ) then
c      do 13211 j = 1,nyp
c      do 13211 k = 1,nzp
c      do 13211 i = 1,nxp
c      pltdp(i,j,k) = 0.0
c13211 continue
c      j = 1
c      do jy = 1,ny
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      if (  algt(ix,jy,kz) .lt.  -0.5 .and.
c     >      algt(ix,jy,kz) .gt. -1.5 )  then
c      pltdp(i,j,k) = pltdp(i,j,k) + algt(ix,jy,kz)
c      end if
c      enddo
c      enddo
c      enddo

c      end if
c
c      if ( idx .eq. 322 ) then
c      do 13221 j = 1,nyp
c      do 13221 k = 1,nzp
c      do 13221 i = 1,nxp
c      pltdp(i,j,k) = 0.0
c13221 continue
c      j = 1
c      do jy = 1,ny
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      if (  algt(ix,jy,kz) .gt.  0.5 .and.
c     >      algt(ix,jy,kz) .lt. 1.5 )  then
c      pltdp(i,j,k) = pltdp(i,j,k) + algt(ix,jy,kz)
c      end if
c      enddo
c      enddo
c      enddo

c      end if
cc
c      if ( idx .eq. 323 ) then
c      do 13231 j = 1,nyp
c      do 13231 k = 1,nzp
c      do 13231 i = 1,nxp
c      pltdp(i,j,k) = 0.0
c13231 continue
c      j = 1
c      do jy = 1,ny
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      if (  algt(ix,jy,kz) .gt.  1.5 )  then
c      pltdp(i,j,k) = pltdp(i,j,k) + algt(ix,jy,kz)
c      end if
c      enddo
c      enddo
c      enddo

c      end if


      IF ( idx .eq. 324 ) THEN ! POSITIVE LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           pltdp(i,j,k) = Abs( elec(ieflshp)%flt3d(ix,jy,kz) )
         END DO
       END DO
      END DO
      end if
c
c
      IF ( idx .eq. 325 ) THEN ! NEGATIVE LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           pltdp(i,j,k) = Abs(elec(ieflshn)%flt3d(ix,jy,kz))
         END DO
       END DO
      END DO
      end if
c
c
      IF ( idx .eq. 326 ) THEN ! INITIATION POINTS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           pltdp(i,j,k) = elec(ieinit)%flt3d(ix,jy,kz)
         END DO
       END DO
      END DO
      end if
c
c
c
c      isum = 3
      IF ( idx .eq. 328 ) THEN ! NegaTIVE LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshn)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO
           ELSEIF ( nzp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(kz-isum,1),Min(kz+isum,nz)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshn)%flt3d(ix,jy,ia) )
             ENDDO
           ELSEIF ( ipltyp .eq. 5 ) THEN
             pltdp(i,j,k) = 0.0
             IF ( phi .gt. -22.5 .and. phi .lt. 22.5 ) THEN
             ! treat as x-z slice
               DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
                pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshn)%flt3d(ix,ia,kz) )
               ENDDO
             ELSEIF ( phi .le. -22.5 .and. phi .ge. -67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshn)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
                            
             ELSEIF ( phi .ge. 22.5 .and. phi .lt. 67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshn)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
             
             ELSEIF ( phi .gt. 67.5 .or. phi .lt. -67.5 ) THEN
             ! treat as y-z slice
             
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO 
             ENDIF            
           ELSE
           pltdp(i,j,k) = Abs( elec(ieflshn)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if
c
c
      IF ( idx .eq. -327 ) THEN ! positive LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ia,jy,kz) )
             ENDDO
           ELSE
           pltdp(i,j,k) = Abs( elec(ieflshp)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if

      IF ( idx .eq. 327 ) THEN ! positive LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ia,jy,kz) )
             ENDDO
           ELSEIF ( nzp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(kz-isum,1),Min(kz+isum,nz)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,jy,ia) )
             ENDDO
           ELSEIF ( ipltyp .eq. 5 ) THEN
             pltdp(i,j,k) = 0.0
             IF ( phi .gt. -22.5 .and. phi .lt. 22.5 ) THEN
             ! treat as x-z slice
               DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
                pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,ia,kz) )
               ENDDO
             ELSEIF ( phi .le. -22.5 .and. phi .ge. -67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshp)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
                            
             ELSEIF ( phi .ge. 22.5 .and. phi .lt. 67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshp)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
             
             ELSEIF ( phi .gt. 67.5 .or. phi .lt. -67.5 ) THEN
             ! treat as y-z slice
             
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ia,jy,kz) )
             ENDDO 
             ENDIF            
           ELSE
           pltdp(i,j,k) = Abs( elec(ieflshp)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if
c
c
c
c      isum = 10
      IF ( idx .eq. 329 ) THEN ! INITIATION POINTS
      write(6,*) '329 -- nyp,jb,nxp,ib,nzp,kb=',nyp,jb,nxp,ib,nzp,kb
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieinit)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieinit)%flt3d(ia,jy,kz) )
             ENDDO
           ELSEIF ( nzp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(kz-isum,1),Min(kz+isum,nz)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieinit)%flt3d(ix,jy,ia) )
             ENDDO
           ELSEIF ( ipltyp .eq. 5 ) THEN
             pltdp(i,j,k) = 0.0
             IF ( phi .gt. -22.5 .and. phi .lt. 22.5 ) THEN
             ! treat as x-z slice
               DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
                pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieinit)%flt3d(ix,ia,kz) )
               ENDDO
             ELSEIF ( phi .le. -22.5 .and. phi .ge. -67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieinit)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
                            
             ELSEIF ( phi .ge. 22.5 .and. phi .lt. 67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieinit)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
             
             ELSEIF ( phi .gt. 67.5 .or. phi .lt. -67.5 ) THEN
             ! treat as y-z slice
             
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieinit)%flt3d(ia,jy,kz) )
             ENDDO 
             ENDIF            
           ELSE
           pltdp(i,j,k) = Abs( elec(ieinit)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if
c
      IF ( idx .eq. 330 ) THEN ! total LEADERS
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,ia,kz) ) + Abs( elec(ieflshn)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ia,jy,kz) ) + Abs( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO
           ELSEIF ( nzp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(kz-isum,1),Min(kz+isum,nz)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,jy,ia) ) + Abs( elec(ieflshn)%flt3d(ix,jy,ia) )
             ENDDO
           ELSEIF ( ipltyp .eq. 5 ) THEN
             pltdp(i,j,k) = 0.0
             IF ( phi .gt. -22.5 .and. phi .lt. 22.5 ) THEN
             ! treat as x-z slice
               DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
                pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ix,ia,kz) ) + Abs( elec(ieflshn)%flt3d(ix,ia,kz) )
               ENDDO
             ELSEIF ( phi .le. -22.5 .and. phi .ge. -67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshp)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) ) + 
     :         Abs( elec(ieflshn)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )               
               ENDDO
                            
             ELSEIF ( phi .ge. 22.5 .and. phi .lt. 67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         Abs( elec(ieflshp)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) ) + 
     :         Abs( elec(ieflshn)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
             
             ELSEIF ( phi .gt. 67.5 .or. phi .lt. -67.5 ) THEN
             ! treat as y-z slice
             
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + Abs( elec(ieflshp)%flt3d(ia,jy,kz) )  + Abs( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO 
             ENDIF            
           ELSE
           pltdp(i,j,k) = Abs( elec(ieflshp)%flt3d(ix,jy,kz) ) + Abs( elec(ieflshn)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if


      IF ( idx .eq. -330 ) THEN ! total LEADERS with negative values for NEG Leaders
      DO j = 1,nyp
        jy = (j-1) + jb
       DO k = 1,nzp
         kz = (k-1) + kb
         DO i = 1,nxp
           ix = (i-1) + ib
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
             pltdp(i,j,k) = pltdp(i,j,k) + ( elec(ieflshp)%flt3d(ix,ia,kz) ) + ( elec(ieflshn)%flt3d(ix,ia,kz) )
             ENDDO
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + ( elec(ieflshp)%flt3d(ia,jy,kz) ) + ( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO
           ELSEIF ( nzp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             DO  ia=Max(kz-isum,1),Min(kz+isum,nz)
             pltdp(i,j,k) = pltdp(i,j,k) + ( elec(ieflshp)%flt3d(ix,jy,ia) ) + ( elec(ieflshn)%flt3d(ix,jy,ia) )
             ENDDO
           ELSEIF ( ipltyp .eq. 5 ) THEN
             pltdp(i,j,k) = 0.0
             IF ( phi .gt. -22.5 .and. phi .lt. 22.5 ) THEN
             ! treat as x-z slice
               DO  ia=Max(jy-isum,1),Min(jy+isum,ny)
                pltdp(i,j,k) = pltdp(i,j,k) + ( elec(ieflshp)%flt3d(ix,ia,kz) ) + ( elec(ieflshn)%flt3d(ix,ia,kz) )
               ENDDO
             ELSEIF ( phi .le. -22.5 .and. phi .ge. -67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         ( elec(ieflshp)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) ) + 
     :         ( elec(ieflshn)%flt3d( Max(1,Min((ix+ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )               
               ENDDO
                            
             ELSEIF ( phi .ge. 22.5 .and. phi .lt. 67.5 ) THEN
             pltdp(i,j,k) = 0.0
               DO  ia=-isum,isum
                pltdp(i,j,k) = pltdp(i,j,k) + 
     :         ( elec(ieflshp)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) ) + 
     :         ( elec(ieflshn)%flt3d( Max(1,Min((ix-ia),nx)) ,
     :                      Max(1,Min((jy+ia),ny)) , kz ) )
               ENDDO
             
             ELSEIF ( phi .gt. 67.5 .or. phi .lt. -67.5 ) THEN
             ! treat as y-z slice
             
             DO  ia=Max(ix-isum,1),Min(ix+isum,nx)
             pltdp(i,j,k) = pltdp(i,j,k) + ( elec(ieflshp)%flt3d(ia,jy,kz) )  + ( elec(ieflshn)%flt3d(ia,jy,kz) )
             ENDDO 
             ENDIF            
           ELSE
           pltdp(i,j,k) = ( elec(ieflshp)%flt3d(ix,jy,kz) ) + ( elec(ieflshn)%flt3d(ix,jy,kz) )
           ENDIF
         END DO
       END DO
      END DO
      end if

c
      if ( idx .eq. 393 ) then ! cloud + rain charge
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = (st(ix,jy,kz,lscw)+st(ix,jy,kz,lscr))*1.0e9
      ENDDO
      ENDDO
      ENDDO
      end if

      if ( idx .eq. 394 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) =
     >              ( st(ix,jy,kz,lsci)
     >              + st(ix,jy,kz,lscir)
     >              + st(ix,jy,kz,lscip))*1.e9
      enddo
      enddo
      enddo

      endif
c
c
      if ( idx .eq. 395 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( microp(1:5) == 'ICE10') THEN
      pltdp(i,j,k) =
     >              ( st(ix,jy,kz,lscgl)
     >              + st(ix,jy,kz,lscgm)
     >              + st(ix,jy,kz,lscgh)
     >              + st(ix,jy,kz,lscf)
     >              + st(ix,jy,kz,lschl)
     >              + st(ix,jy,kz,lsch))*1.0e9
      ELSE
      pltdp(i,j,k) =
     >               (st(ix,jy,kz,lschl)
     >              + st(ix,jy,kz,lsch))*1.0e9
      ENDIF
      enddo
      enddo
      enddo

      endif
c
      if ( idx .eq. 396 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lschl)*1.e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 399 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscw)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 400 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lsci)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 401 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscip)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 402 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscir)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 403 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscr)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 404 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscs)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 405 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscgl)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 406 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscgm)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 407 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscgh)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 408 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lscf)*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 409 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lsch)*1.0e9
      enddo
      enddo
      enddo

      end if

      if ( idx .eq. 984 ) then ! positive ion charge
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ec*(st(ix,jy,kz,lscpi) )*1.0e9
      ENDDO
      ENDDO
      ENDDO
      end if

      if ( idx .eq. 985 ) then ! negative ion charge
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ec*( - st(ix,jy,kz,lscni) )*1.0e9
      ENDDO
      ENDDO
      ENDDO
      end if

      if ( idx .eq. 986 ) then ! net ion charge
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = ec*(st(ix,jy,kz,lscpi) - st(ix,jy,kz,lscni) )*1.0e9
      ENDDO
      ENDDO
      ENDDO
      end if
c
      if ( idx .eq. 410 ) then
      scmax = 0.0
      scmin = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
!      pltdp(i,j,k) = rho(ix,jy,kz)*1.0e9
      scmax = Max(scmax,rho(ix,jy,kz))
      scmin = Min(scmin,rho(ix,jy,kz))
           IF (nyp .eq. 1) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(jy-jyave,1),Min(jy+jyave,ny)
               pltdp(i,j,k) = pltdp(i,j,k) + rho(ix,ia,kz)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e9
           ELSEIF ( nxp .eq. 1 ) THEN
             pltdp(i,j,k) = 0.0
             cnt = 0
             DO  ia=Max(ix-ixave,1),Min(ix+ixave,nx)
               pltdp(i,j,k) = pltdp(i,j,k) + rho(ia,jy,kz)
               cnt = cnt + 1
             ENDDO
             IF ( cnt .gt. 0 ) pltdp(i,j,k) = pltdp(i,j,k)/Float(cnt)*1.0e9
           ELSE
            pltdp(i,j,k) = rho(ix,jy,kz)*1.0e9
           ENDIF

c      (st(ix,jy,kz,lscw)
c     >              + st(ix,jy,kz,lsci)
c     >              + st(ix,jy,kz,lscir)
c     >              + st(ix,jy,kz,lscip)
c     >              + st(ix,jy,kz,lscr)
c     >              + st(ix,jy,kz,lscs)
c     >              + st(ix,jy,kz,lscgl)
c     >              + st(ix,jy,kz,lscgm)
c     >              + st(ix,jy,kz,lscgh)
c     >              + st(ix,jy,kz,lscf)
c     >              + st(ix,jy,kz,lschl)
c     >              + st(ix,jy,kz,lsch))*1.0e9
      enddo
      enddo
      enddo

      write(6,*) 'scmin,scmax,ix,jy,kz',scmin,scmax,ix,jy,kz
      endif
c
      if ( idx .eq. 416 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = 
     >              ( st(ix,jy,kz,lscgl)
     >              + st(ix,jy,kz,lscgm)
     >              + st(ix,jy,kz,lscgh)
     >              + st(ix,jy,kz,lscf))*1.0e9
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 417 ) then
      IF ( lscir .gt. 1 .and. lscip .gt. 1 ) THEN
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = (st(ix,jy,kz,lsci) + st(ix,jy,kz,lscs)
     >              + st(ix,jy,kz,lscir)
     >              + st(ix,jy,kz,lscip))*1.e9
      enddo
      enddo
      enddo

      ELSE
      DO j = 1,nyp
      jy = (j-1) + jb
      DO k = 1,nzp
      kz = (k-1) + kb
      DO i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = (st(ix,jy,kz,lsci) + st(ix,jy,kz,lscs))*1.0e9
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF
      end if
c
c      if ( idx .eq. 418 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = sch(ix,jy,kz)*1.0e9
c      enddo
c      enddo
c      enddo

c      end if
c
      if ( idx .eq. 411 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = elec(iex)%flt3d(ix,jy,kz)*1.e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 412 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = elec(iey)%flt3d(ix,jy,kz)*1.e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 413 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = elec(iez)%flt3d(ix,jy,kz)*1.e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 419 ) then
       do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
         kz = (k-1) + kb
         do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = elec(iez)%flt3d(ix,jy,kz)*1.e-3
         ENDDO
        ENDDO
       ENDDO
      end if
c
      if ( idx .eq. 414 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = elec(iemag)%flt3d(ix,jy,kz)*0.001
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 415 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = elec(ipot)%flt3d(ix,jy,kz)*1.0e-6
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 420 ) then
      IF ( microp(1:3) == 'TAK' ) THEN
      do  j = 1,nyp
      jy = (j-1) + jb
      do  k = 1,nzp
      kz = (k-1) + kb
      do  i = 1,nxp
      ix = (i-1) + ib

       tmp = 0.0
         DO l = 1,lr100-1
           tmp = tmp + st(ix,jy,kz,lnr-1+l)
         ENDDO

      pltdp(i,j,k) = tmp
      
      enddo
      enddo
      enddo
      
      ELSE
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnc)*1.0e-6
      enddo
      enddo
      enddo

      ENDIF
      end if
c
      if ( idx .eq. 421 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lni)*1.0e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 422 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnip)*1.0e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 423 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      IF ( lnis >= 1 ) THEN
      pltdp(i,j,k) = st(ix,jy,kz,lnis)*1.0e-3
      ELSE
      pltdp(i,j,k) = st(ix,jy,kz,lnir)*1.0e-3
      ENDIF
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 424 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnr)*1.0e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 425 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lns)*1.0e-3
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 426 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lngl)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 427 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lngm)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 428 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lngh)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 429 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnf)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 430 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnh)*1.0
      enddo
      enddo
      enddo

      end if

      if ( idx == 434 .and. lsat >= 1 ) then
      xmx = 0.0
      xmn = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
       do k = 1,nzp
       kz = (k-1) + kb
        do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = xtra(lsat)%flt3d(ix,jy,kz)
        enddo
       enddo
      enddo
      end if

      if ( idx == 435 .and. lsati >= 1 ) then
      xmx = 0.0
      xmn = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
       do k = 1,nzp
       kz = (k-1) + kb
        do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = xtra(lsati)%flt3d(ix,jy,kz)
        enddo
       enddo
      enddo
      end if
c
      if ( idx .eq. 440 .and. lnchaff .gt. 1 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lnchaff)*1000.
      enddo
      enddo
      enddo

      end if
c
      if ( idx == 441 .and. lnox > 1 ) then
      write(6,*) 'LNOX: factor is 1.e9, lnox = ',lnox,na
      xmx = 0.0
      xmn = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
       do k = 1,nzp
       kz = (k-1) + kb
        do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = st(ix,jy,kz,lnox)*1.e9
          xmx = Max(xmx, pltdp(i,j,k))
          xmn = Min(xmn, pltdp(i,j,k))
        enddo
       enddo
      enddo
      write(6,*) 'LNOX: xmn,xmx = ',xmn,xmx
      end if

      if ( idx == 442 .and. lanox > 1 ) then
      write(6,*) 'LANOx: factor is 1.e9, lanox = ',lanox,na
      xmx = 0.0
      xmn = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
       do k = 1,nzp
       kz = (k-1) + kb
        do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = st(ix,jy,kz,lanox)*1.e9
          xmx = Max(xmx, pltdp(i,j,k))
          xmn = Min(xmn, pltdp(i,j,k))
        enddo
       enddo
      enddo
      write(6,*) 'LANOX: xmn,xmx = ',xmn,xmx
      end if

      if ( idx == 443 .and. lco > 1 ) then
      write(6,*) 'CO: factor is 1.e9, lco = ',lco,na
      xmx = 0.0
      xmn = 0.0
      do j = 1,nyp
      jy = (j-1) + jb
       do k = 1,nzp
       kz = (k-1) + kb
        do i = 1,nxp
          ix = (i-1) + ib
          pltdp(i,j,k) = st(ix,jy,kz,lco)*1.e9
          xmx = Max(xmx, pltdp(i,j,k))
          xmn = Min(xmn, pltdp(i,j,k))
        enddo
       enddo
      enddo
      write(6,*) 'CO: xmn,xmx = ',xmn,xmx
      end if


!      if ( idx .eq. 441 ) then
!      do j = 1,nyp
!      jy = (j-1) + jb
!      do k = 1,nzp
!      kz = (k-1) + kb
!      do i = 1,nxp
!      ix = (i-1) + ib
!      pltdp(i,j,k) = st(ix,jy,kz,ltpi)*1.0
!      enddo
!      enddo
!      enddo

!      end if
c
      if ( idx .eq. 442 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpip)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 443 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpir)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 444 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpr)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 445 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltps)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 446 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpgl)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 447 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpgm)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 448 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpgh)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 449 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,ltpf)*1.0
      enddo
      enddo
      enddo

      end if
c
!      if ( idx .eq. 450 ) then
!      do j = 1,nyp
!      jy = (j-1) + jb
!      do k = 1,nzp
!      kz = (k-1) + kb
!      do i = 1,nxp
!      ix = (i-1) + ib
!      pltdp(i,j,k) = st(ix,jy,kz,ltph)*1.0
!      enddo
!      enddo
!      enddo

!      end if
c
c
      IF ( idx == 450 ) THEN ! cloud top temp
      
      DO jy = 1,ny-1
       DO ix = 1,nx-1
        pltdp(ix,jy,nz) = tt0(ix,jy,1) - 273.15
         DO kz = nz-1,1,-1
          IF ( st(ix,jy,kz,li) + st(ix,jy,kz,lc) > 0.04e-3 ) THEN
            pltdp(ix,jy,nz) = tt0(ix,jy,Max(1,kz-1)) - 273.15
            EXIT
          ENDIF
         ENDDO
         
       ENDDO
      ENDDO
      
      
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      
      
      pltdp(i,j,k) = pltdp(ix,jy,nz) ! precip(iflshr)%flt2d(ix,jy)
      
      ENDDO
      ENDDO
      ENDDO
      ENDIF


      if ( idx .eq. 451 ) then
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib

            pltdp(i,j,k) = st(ix,jy,kz,lccn)*1.0e-6
      
          enddo
        enddo
      enddo
      end if

      if ( idx .eq. 452 .and. lccna > 1 ) then
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib

            pltdp(i,j,k) = st(ix,jy,kz,lccna)*1.0e-6
      
          enddo
        enddo
      enddo
      end if

      if ( idx .eq. 453 .and. lcina > 1 ) then
      do j = 1,nyp
        jy = (j-1) + jb
        do k = 1,nzp
          kz = (k-1) + kb
          do i = 1,nxp
          ix = (i-1) + ib

            pltdp(i,j,k) = st(ix,jy,kz,lcina)*1.0e-3
      
          enddo
        enddo
      enddo
      end if



      if ( idx .eq. 460 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lamcw)*1.0
      enddo
      enddo
      enddo

      end if
c
c      if ( idx .eq. 461 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamci)*1.0
c      enddo
c      enddo
c      enddo

c      end if
c
c      if ( idx .eq. 462 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcip)*1.0
c      enddo
c      enddo
c      enddo

c      end if
c
      if ( idx .eq. 463 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcir)*1.0
      enddo
      enddo
      enddo

      end if
c
c      if ( idx .eq. 464 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcr)*1.0
c      enddo
c      enddo
c      enddo

c      end if
c
c      if ( idx .eq. 465 ) then
c      do j = 1,nyp
c      jy = (j-1) + jb
c      do k = 1,nzp
c      kz = (k-1) + kb
c      do i = 1,nxp
c      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcs)*1.0
c      enddo
c      enddo
c      enddo

c      end if
c
      if ( idx .eq. 466 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcgl)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 467 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcgm)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 468 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcgh)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 469 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamcf)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 470 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
c      pltdp(i,j,k) = st(ix,jy,kz,lamch)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 480 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtc)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 481 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrti)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 482 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtip)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 483 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtir)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 484 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtr)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 485 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrts)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 486 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtgl)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 487 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtgm)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 488 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtgh)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 489 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrtf)*1.0
      enddo
      enddo
      enddo

      end if
c
      if ( idx .eq. 490 ) then
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = st(ix,jy,kz,lrth)*1.0
      enddo
      enddo
      enddo

      end if
c
c  micro-physics stuff
c 
      if ( .false. .and. idx .ge. 47 .and. idx .le. 139 ) then
      allocate( dum(nx,ny,nz) )
      write(0,*) 'BASP3P: Special microphysics readmp'
      call readmp( nx, ny, nz, dum )
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      pltdp(i,j,k) = dum(ix,jy,kz)
      enddo
      enddo
      enddo

      deallocate( dum )
      end if


c
c  transfer values to the 2d slab
c
      xmx = -1.0e29
      xmn =  1.0e29
      IF ( ipltyp .ne. 5 ) THEN
      do j = 1,nyp
      jy = (j-1) + jb
      do k = 1,nzp
      kz = (k-1) + kb
      do i = 1,nxp
      ix = (i-1) + ib
      plt2dp(i,j,k) = Max( 0.0, pltdp(i,j,k) )
      plt2dn(i,j,k) = Min( 0.0, pltdp(i,j,k) )
c
c  find min and max
c
      if (pltdp(i,j,k).lt.xmn) then
      xmn = pltdp(i,j,k)
      imn = ix
      jmn = jy
      kmn = kz
      end if
      if (plt2dp(i,j,k).gt.xmx) then
      xmx = pltdp(i,j,k)
      imx = ix
      jmx = jy
      kmx = kz
      end if
c separate positive and negative values      
!      plt2dn(i,j,k) = 0.0
!      if ( plt2dp(i,j,k) .lt. 0.0 ) then
!      plt2dn(i,j,k) = plt2dp(i,j,k)
!      plt2dp(i,j,k) = 0.0 
!      end if
      
      END DO
      END DO
      END DO
      
      ELSEIF ( ipltyp .eq. 5 )  THEN
c interpolate slice into pltdp2           
       DO j=1,1
       DO k=1,nzp
         kz = k  ! dont add (kb-1) because pltdp already starts at kb
       DO i=1,nl
         xi = xib + (i-1 + d1pobl/dx)*Cos(phi*d2r)
         yj = yjb + (i-1 + d1pobl/dx)*Sin(phi*d2r)*(dx/dy)
         ix = Int(xi)
         jy = Int(yj)
         delx = xi - Float(ix)
         dely = yj - Float(jy)
         f2 = (1.0 - delx)*pltdp(ix,jy,kz) + delx*pltdp(ix+1,jy,kz) 
         f3 = (1.0 - delx)*pltdp(ix,jy+1,kz) + delx*pltdp(ix+1,jy+1,kz)
         f4 = (1.0 - dely)*f2 + dely*f3
         plt2dp(i,j,k) = f4 
         if (plt2dp(i,j,k).lt.xmn) then
          xmn = plt2dp(i,j,k)
          imn = ix
          jmn = jy
          kmn = kz
         end if
         if (plt2dp(i,j,k).gt.xmx) then
           xmx = plt2dp(i,j,k)
           imx = ix
           jmx = jy
           kmx = kz
         end if
c separate positive and negative values      
        plt2dn(i,j,k) = Min( 0.0, f4 )
        plt2dp(i,j,k) = Max( 0.0, plt2dp(i,j,k) )
!        if ( plt2dp(i,j,k) .lt. 0.0 ) then
!         plt2dn(i,j,k) = plt2dp(i,j,k)
!         plt2dp(i,j,k) = 0.0 
!        end if
       END DO
       END DO
       END DO
      END IF ! ( ipltyp .ne. 5 )
      
      write(6,*) 'ib,jb,kb,nxp, nyp, nzp', ib,jb,kb,nxp,nyp,nzp
      write(6,*) 'xib,yjb,zkb,zke',xib,yjb,zkb,zke
      write(6,'(a,1(1pe12.5),a,4i6)') 'xmn: ',xmn,' at ',imn,jmn,kmn
      write(6,'(a,1(1pe12.5),a,4i6)') 'xmx: ',xmx,' at ',imx,jmx,kmx
c
c  end
c
      deallocate( pltdp )

      return
      end 
c
c
c

