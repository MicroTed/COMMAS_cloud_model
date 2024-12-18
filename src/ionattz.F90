!
!  subroutine to calculate fall speeds of hydrometeors and ion attachment
!
      subroutine ionattz(nx,ny,nz,nor,na,nba,dtp,dz,jgs,     &
     &  elecez,attach,uz,diff,cond,conc,     &
     &  an,dn,db,t0,t7,     &
     &  qxmin)
!
! 8/17/2024: Removed unused variables; switched GS arrays to allocatable 
!            and convert to column-wise calculation (removed complicated GS logic)
!
!
      USE ELEC_MODULE, only: do_ionatt_conduction, do_ionatt_diffusion, &
                            ionatt_conduction_factor, ionatt_diffusion_factor
      USE INDEX_MODULE
      USE MICRO_MODULE
      USE COMMASMPI_MODULE
      implicit none

      integer nxmpb,ixtmp,nzmpb,nxz,numgs,inumgs
      integer kstag,istag
      parameter (kstag=1, istag=1)
!
      integer ndebug2
      parameter (ndebug2 = 0)

      integer nx,ny,nz,nor,ngt,jgs,na,nba

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)

      real pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)


      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dtp,dz,dtz1

      real elecez(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      integer idx

      integer ix,jy,kz,i,j,k,item,ia,n
      
      integer ngs1
      real cwc1,cwc0
      real chw, qr, z, rd, alp, z1, g1, vr, nrx, frac, frach

      real attach(nx,nz,lc:lhab,2)

      real conc(nx,nz,lc:lhab) ! concentration

      real diff(nx,nz,lc:lhab,2)
      real cond(nx,nz,lc:lhab,2)
 
      integer ngs,ngscnt,mgs,count
      double precision :: atotpos,atotneg

      real, allocatable :: axx(:,:),bxx(:,:)

      real, allocatable ::  qx(:,:)
      real, allocatable ::  qxw(:,:)
      real, allocatable ::  cx(:,:)
      real, allocatable ::  xv(:,:)
      real, allocatable ::  vtxbar(:,:,:)
      real, allocatable ::  xmas(:,:)
      real, allocatable ::  xdn(:,:)
      real, allocatable ::  xdia(:,:,:)
      real, allocatable ::  vx(:,:)
      real, allocatable ::  alpha(:,:)
      real, allocatable ::  zx(:,:), scx(:,:)

      integer, allocatable :: igs(:),kgs(:)
      
      real, allocatable :: rho0(:),temcg(:)

      real, allocatable :: temg(:),ez(:)
      
      real, allocatable :: rhovt(:)
      
      real, allocatable :: cwnc(:),cinc(:),cdxgs(:,:)
      real, allocatable :: fadvisc(:),cwdia(:),cipmas(:)
      
      real, allocatable :: cnina(:)
      real, allocatable :: cnostmp(:)

      real xdnmx(lc:lhab), xdnmn(lc:lhab)

      real qxmin(lc:lhab)
      real cdx(lc:lhab)

      real cno(lc:lhab)
      real xdn0(lc:lhab)

      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)

      real, allocatable :: xrad(:,:), scppx(:,:)
      real, allocatable :: cion(:,:)  ! 1=pos ions; 2=neg. ions
      real, allocatable :: scxaciond(:,:,:), scxacionc(:,:,:)
      
      double precision, allocatable :: diffco(:,:)  ! diffusion coefficient for pos (1) and neg (2) ions
      real uz(nz,4)  ! ion mobilities
      double precision s(2)
      
      double precision x,y,p1,p2,p3
      double precision qd,qm
      
      real kb   ! Boltzman constant
      parameter (kb = 1.3807e-23)
      
      real ec  ! fundamental unit of charge
      parameter (ec = 1.602e-19)

      real eperao
      parameter (eperao  = 8.8592e-12 )
      
      real, parameter :: pi = 3.141592653589793, pii = 1.0d0/pi
      real, parameter :: difffac = 2.0d0*pi*kb/ec

      real rho00
      real xvmn(lc:lhab), xvmx(lc:lhab)
!
!  general constants for microphysics
!
      real tfr
      real advisc0, advisc1

! 
! Miscellaneous
!
      
      real cnox, tmp, tmp2
      integer il
      real cwnccn(nz)

      logical flag
      logical ldovol, ldoliq
      integer lvol(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze


! #####################################################################
! BEGIN EXECUTABLE
! #####################################################################
!
      s(1) = 1.0d0 ; s(2) = -1.0d0
      
      ngs = nz


      allocate( xrad(ngs,lc:lhab), scppx(ngs,lc:lhab) )
      allocate( diffco(ngs,2), cion(ngs,2), scxaciond(ngs,lc:lhab,2), scxacionc(ngs,lc:lhab,2) )
      allocate( qx(ngs,lv:lhab),  &
                qxw(ngs,ls:lhab),  &
                cx(ngs,lc:lhab),  &
                xv(ngs,lc:lhab),  &
                vtxbar(ngs,lc:lhab,3),  &
                xmas(ngs,lc:lhab),  &
                xdn(ngs,lc:lhab),  &
                xdia(ngs,lc:lhab,3),  &
                vx(ngs,li:lhab),  &
                alpha(ngs,lc:lhab),  &
                zx(ngs,lr:lhab),     &
                scx(ngs,lc:lhab),    &
                igs(ngs),kgs(ngs), &
                ez(ngs),axx(ngs,lh:lhab),bxx(ngs,lh:lhab), &
                rho0(ngs),temcg(ngs),temg(ngs), rhovt(ngs), &
                fadvisc(ngs), cdxgs(ngs,lc:lhab), &
                cnina(ngs), &
                cnostmp(ngs) )

      ldovol = .false.
      lvol(:) = 0
      IF ( lvi .gt. 1 ) lvol(li) = lvi
      IF ( lvs .gt. 1 ) lvol(ls) = lvs
      IF ( lvh .gt. 1 ) lvol(lh) = lvh
      
      
      IF ( li .gt. 1 ) THEN
      DO il = li,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO
      ENDIF

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
      IF ( lfw .gt. 1 .and. lf > 1 ) lliq(lf) = lfw
      IF ( lhl .gt. 1 .and. lhlw .gt. 1 ) lliq(lhl) = lhlw

      ldoliq = .false.
      IF ( ls .gt. 1 ) THEN
      DO il = ls,lhab
        ldoliq = ldoliq .or. ( lliq(il) .gt. 1 )
      ENDDO
      ENDIF
 
      cwnccn(:) = cwccn
      
      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn
      xvmn(lh) = xvhmn
      IF ( lf .gt. 1 ) xvmn(lf) = xvhmn ! (yes, this is right -- not large hail here)
      IF ( lhl .gt. 1 ) xvmn(lhl) = xvhmn ! (yes, this is right -- not large hail here)

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx
      xvmx(ls) = xvsmx
      xvmx(lh) = xvhmx
      IF ( lf .gt. 1 ) xvmx(lf) = xvhmx
      IF ( lhl .gt. 1 ) xvmx(lhl) = xvhlmx

      xvmn(li) = xvimn
      xvmx(li) = xvimx

      lsc(lc) = lscw
      lsc(lr) = lscr
      lsc(li) = lsci
      lsc(ls) = lscs
      lsc(lh) = lsch
      IF ( lf .gt. 1 ) lsc(lf) = lscf
      IF ( lhl .gt. 1 ) lsc(lhl) = lschl
      IF ( lis .gt. 1 ) lsc(lis) = lscis

      ln(lc) = lnc
      ln(lr) = lnr
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      IF ( lf .gt. 1 ) ln(lf) = lnf
      IF ( lhl .gt. 1 ) ln(lhl) = lnhl
      IF ( lis .gt. 1 ) ln(lis) = lnis

      ipc(lc) = 2
      ipc(lr) = 3
      ipc(li) = 1
      ipc(ls) = 4
      ipc(lh) = 5
      IF ( lf .gt. 1 ) ipc(lf) = 5 
      IF ( lhl .gt. 1 ) ipc(lhl) = 5 
      IF ( lis .gt. 1 ) ipc(lis) = 1
 
!  constants
!
      rho00 = 1.225
!
!  slope intercepts
!
      CALL setcnoz(cno)
!
!  density maximums and minimums
!
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(ls) =  200.0
      xdnmx(lh) =  900.0
!
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lh) =  170.0

      xdn0(lc) = 1000.0
      xdn0(lr) = 1000.0
      xdn0(li) = 900.0
      xdn0(ls) = rho_qs ! 100.0
      xdn0(lh) = rho_qh ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      IF ( lf .gt. 1 ) THEN
        xdnmx(lf) =  rho_qf_max
        xdnmn(lf) =  fdnmn
        xdn0(lf) = rho_qf
      ENDIF
      IF ( lhl .gt. 1 ) THEN
        xdnmx(lhl) =  900.0
        xdnmn(lhl) =  hldnmn
        xdn0(lhl) = rho_qhl
      ENDIF

      IF ( lis > 1 ) THEN
        xvmn(lis) = xvimn
        xvmx(lis) = xvimx
        xdn0(lis) = 917.0
        xdnmx(lis) =  917.0
        xdnmn(lis) =  100.0
      ENDIF

!
!  Set terminal velocities...
!    also set drag coefficients
!
      cdx(lr) = 0.60
      cdx(lh) = 0.8 ! 1.0 ! 0.45
      IF ( lhl .gt. 1 ) cdx(lhl) = 1.0 ! 1.0 ! 0.45
      cdx(ls) = 2.00
!
!  general constants for microphysics
!
      tfr = 273.15
      advisc0 = 1.832e-05
      advisc1 = 1.718e-05


      jy = jgs

      DO ia=lc,lhab
       DO kz = 1,nz
        DO ix = 1,nx
!          vt(ix,kz,ia) = 0.0
!          rad(ix,kz,ia) = 0.0
!          conc(ix,kz,ia) = 0.0
!          attach(ix,kz,ia,1) = 0.0 ! already zero from ionstep.F90
!          attach(ix,kz,ia,2) = 0.0
        ENDDO ! ia
       ENDDO ! kz
      ENDDO ! ix


      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

!       atotpos = 0
!       atotneg = 0
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do ix = 1,ixe
       ngscnt = 0
       do kz = 1,kze-1 

      flag = .false.
      
      DO il = lc,lhab
        flag =  flag .or. ( an(ix,jy,kz,il) .gt. qxmin(il) ) 
      ENDDO

      if ( flag ) then
       ! load temp quantities
        ngscnt = ngscnt + 1
        igs(ngscnt) = ix
        kgs(ngscnt) = kz
      end if
      
      enddo ! kz

      if ( ngscnt .eq. 0 ) CYCLE ! skip to next ix column

      xv(:,:) = 0.0
      qx(:,:) = 0.0
      scx(:,:) = 0.0
      zx(:,:) = 0.0
      alpha(:,:) = 0.0
      xmas(:,:) = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:) = 0.0
      xdn(:,:) = 0.0
      cion(:,:) = 0.0
      diffco(:,:) = 0.0
      scxaciond(:,:,:) = 0.0
      scxacionc(:,:,:) = 0.0

!
!  set temporaries for microphysics variables
!
      DO il = lv,lhab
        DO mgs = 1,ngscnt
          qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
        ENDDO
      ENDDO

      
      DO il = lc,lhab
        DO mgs = 1,ngscnt
         scx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lsc(il))
        ENDDO
      ENDDO
      
      DO mgs = 1,ngscnt
        cion(mgs,1) = an(igs(mgs),jy,kgs(mgs),lscpi)
        cion(mgs,2) = an(igs(mgs),jy,kgs(mgs),lscni)
        ez(mgs) = elecez(igs(mgs),jy,kgs(mgs))
      ENDDO

       
      DO il = lg,lhab
      DO mgs = 1,ngscnt
        alpha(mgs,il) = dnu(il)
      ENDDO
      ENDDO
       
      alpha(:,lr) = xnu(lr)
       
!
!  set concentrations
!
      cx(:,:) = 0.0
      
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
        IF ( lis > 1 ) cx(mgs,lis) = Max(an(igs(mgs),jy,kgs(mgs),lnis), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
       end do

      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
          cx(mgs,lr) = 0.0
        ELSEIF ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) THEN
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

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
          cx(mgs,lh) = 0.0
        ELSEIF ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lh) = 0.0
        ELSE
          cx(mgs,lh) = Max( 1.e-3, cx(mgs,lh) )
        ENDIF
       end do
      end if

      if ( lf .gt. 1 .and. ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lf) = Max(an(igs(mgs),jy,kgs(mgs),lnf), 0.0)
        IF ( qx(mgs,lf) .le. qxmin(lf) ) THEN
          cx(mgs,lf) = 0.0
        ELSEIF ( cx(mgs,lf) .eq. 0.0 .and. qx(mgs,lf) .lt. 3.0*qxmin(lf) ) THEN
          qx(mgs,lf) = 0.0
        ELSE
          cx(mgs,lf) = Max( 1.e-3, cx(mgs,lf) )
        ENDIF
       end do
      end if

      if ( lhl .gt. 1 .and. ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        IF ( qx(mgs,lhl) .le. qxmin(lhl) ) THEN
          cx(mgs,lhl) = 0.0
        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) THEN
          qx(mgs,lhl) = 0.0
        ELSE
          cx(mgs,lhl) = Max( 1.e-3, cx(mgs,lhl) )
        ENDIF
       end do
      end if

       
      do mgs = 1,ngscnt
        xdn(mgs,li) = xdn0(li)
        IF ( lis > 1 ) xdn(mgs,lis) = xdn0(lis)
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
        xdn(mgs,ls) = xdn0(ls)
        xdn(mgs,lh) = xdn0(lh)
        IF ( lf .gt. 1 ) xdn(mgs,lf) = xdn0(lf)
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

!
! Set 6th moments
!
      IF ( ipconc .ge. 6 .or. lzr > 1) THEN
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

!
!  Reconstruct various quantities 
!
      
      do mgs = 1,ngscnt
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhovt(mgs) = Sqrt(rho00/rho0(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
      end do

      ! viscosity for droplet fall speed
      do mgs = 1,ngscnt
       fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*(temg(mgs)/296.0)**(1.5)
      end do
!


      IF ( ipconc .ge. 6 ) THEN

!  Find shape parameters by iteration for rain, graupel,hail, fd 
!  Assuming that imurain = 1

        DO il = lr,lhab

        
        IF (  lz(il) .gt. 1 .and. ( .not. ( il == lr .and. imurain == 3 )) ) THEN
        
        DO mgs = 1,ngscnt

        IF ( qx(mgs,il) .gt. qxmin(il) .and. cx(mgs,il) .gt. 0.0 .and. zx(mgs,il) > 0.0 ) THEN
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          IF ( zx(mgs,il) .gt. 0. ) THEN

           rd = z*(pi/6.*xdn(mgs,il))**2*chw/((dn(igs(mgs),jy,kgs(mgs))*qr)**2)

           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           alp = Max( alphamin, Min( alphamax, alp ) )

           DO i = 1,10
            IF ( i > 1 .and. Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

          ENDIF
          ENDIF


        ENDDO ! mgs
        
        ENDIF ! lz(il) .gt. 1
        
        ENDDO ! il

      ENDIF



      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,     &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,     &
     &                 ipconc,ndebug,ngs,nz,igs,kgs,cwnccn,fadvisc,     &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,     &
     &                 itype1,itype2,temcg,0,alpha,alpha,axx,bxx,0)


            


      DO il = lc,lhab
        xrad(1:ngscnt,il) = 0.5*xdia(1:ngscnt,il,3)
      ENDDO
      
      DO ia=lc,lhab
      DO mgs = 1,ngscnt
        scppx(mgs,ia) = 0.0
      ENDDO
      ENDDO

!
! ion diffusion/conduction
!      
      do mgs = 1,ngscnt
!        diffco(mgs,1) = 2.0*pi*uz(kgs(mgs),1)*kb*temg(mgs)/ec
!        diffco(mgs,2) = 2.0*pi*uz(kgs(mgs),2)*kb*temg(mgs)/ec
        diffco(mgs,1) = difffac*uz(kgs(mgs),1)*temg(mgs)
        diffco(mgs,2) = difffac*uz(kgs(mgs),2)*temg(mgs)
        
          DO ia=lc,lhab
             scxaciond(mgs,ia,:) = 0.0
             scxacionc(mgs,ia,:) = 0.0

            IF ( qx(mgs,ia) .gt. qxmin(ia) ) THEN

!
! diffusion
!            
            IF ( do_ionatt_diffusion ) THEN

            DO n=1,2
            qd = 4.0*pi*eperao*xrad(mgs,ia)*kb*temg(mgs)/ec
!            qd = 9.59491613906812d-15*xrad(mgs,ia)*temg(mgs)
            
            IF ( qd .gt. 1.0d-20 .and. cx(mgs,ia) .gt. 1.0e-6 ) THEN
            scppx(mgs,ia) = scx(mgs,ia)/cx(mgs,ia)
            
            x = scppx(mgs,ia)/qd ! is the X in eq. A9 of Mansell et al. 2005, where qd is the "balance" charge
            
            p1 = 2.0d0*xrad(mgs,ia)*diffco(mgs,n)*cion(mgs,n) ! leading term of A8, except for f(X) and cx,
                                                              ! Note that diffco includes factor of 2*pi
            
            y = s(n)*x
            
!            IF ( s(n)*x .gt. 1.0d-5 .and. (s(n)*x) .lt. 40.0d0 ) THEN
!c              p2 = Min(1.0, s*x/(Exp(s*x) - 1))
!              p2 = Min(1.0D0, s(n)*x/(Exp(s(n)*x) - 1.0D0))
!            ELSEIF ( s(n)*x .ge. 40.0d0 ) THEN
            ! p2 is the f(X) factor (A9) used in A8 (Mansell et al. 2005). Paper says only using y > 0, but
            ! numerically need to have it a little bit larger
            IF ( y .gt. 1.0d-5 .and. y .lt. 10.0d0 ) THEN
!              p2 = Min(1.0, s*x/(Exp(s*x) - 1))
              p2 = Min(1.0D0, y/(Exp(y) - 1.0D0))
            ELSEIF ( y .ge. 10.0d0 ) THEN
              p2 = 0.0d0
            ELSE
              p2 = 1.0d0
            ENDIF
            
            
            p3 = 1.0d0 + Sqrt(xrad(mgs,ia)*vtxbar(mgs,ia,1)/diffco(mgs,n))

            scxaciond(mgs,ia,n) = p1*p2*p3*ionatt_diffusion_factor
            ELSE
!             write(6,*) 'Huh? ia,qd,cx,q = ',ia,qd,cx(mgs,ia),qx(mgs,ia)
!             write(6,*) xrad(mgs,ia),temg(mgs),vtxbar(mgs,ia,1),qxmin(ia)
             scxaciond(mgs,ia,n) = 0.0
            ENDIF
            
            ENDDO ! n
            
            ENDIF ! do_ionatt_diffusion
!
!  conduction per Chiu (1978)
!
!            scxacionc(mgs,ia,1) = 0.0
!            scxacionc(mgs,ia,2) = 0.0

            IF ( do_ionatt_conduction ) THEN

            qm = 12.0d0*pi*eperao*Abs(ez(mgs))*xrad(mgs,ia)**2 ! eq. 55 from Chiu (1978)
            
              IF ( scppx(mgs,ia) .ge. qm ) THEN ! Chiu eq. 51
                scxacionc(mgs,ia,1) = 0.0
                scxacionc(mgs,ia,2) = cion(mgs,2)     &
     &            *uz(kgs(mgs),2)*scppx(mgs,ia)/eperao
              
              ELSEIF ( scppx(mgs,ia) .le. -qm ) THEN  ! Chiu eq. 52
                scxacionc(mgs,ia,1) = -cion(mgs,1)     &
     &            *uz(kgs(mgs),1)*scppx(mgs,ia)/eperao
                scxacionc(mgs,ia,2) = 0.0
              
              ELSEIF ( ez(mgs) .lt. 0.0 ) THEN ! E parallel to vt; Chiu eq. 53a-d
              
                scxacionc(mgs,ia,2) = 3.0*pi*xrad(mgs,ia)**2     &
     &            *cion(mgs,2)*Abs(ez(mgs))     &
     &            *uz(kgs(mgs),2)*(1.0 + scppx(mgs,ia)/qm)**2
                
                IF ( vtxbar(mgs,ia,1) .gt. uz(kgs(mgs),1)*Abs(ez(mgs)) )     &
     &                   THEN
                  IF ( scppx(mgs,ia) .gt. 0.0 ) THEN
                    scxacionc(mgs,ia,1) = 0.0
                  ELSE
                    scxacionc(mgs,ia,1) = -uz(kgs(mgs),1)*cion(mgs,1)*     &
     &                scppx(mgs,ia)/eperao
                  ENDIF
                ELSE  
                scxacionc(mgs,ia,1) = 3.0*pi*xrad(mgs,ia)**2     &
     &            *cion(mgs,1)*Abs(ez(mgs))     &
     &            *uz(kgs(mgs),1)*(1.0 - scppx(mgs,ia)/qm)**2
                  
                ENDIF
              ELSEIF ( ez(mgs) .gt. 0.0 ) THEN ! E anti-parallel to vt; Chiu eq. 54a-d

                scxacionc(mgs,ia,1) = 3.0*pi*xrad(mgs,ia)**2     &
     &            *cion(mgs,1)*Abs(ez(mgs))     &
     &            *uz(kgs(mgs),1)*(1.0 - scppx(mgs,ia)/qm)**2

                IF ( vtxbar(mgs,ia,1) .gt. uz(kgs(mgs),2)*Abs(ez(mgs)) )     &
     &               THEN
                  IF ( scppx(mgs,ia) .lt. 0.0 ) THEN
                    scxacionc(mgs,ia,2) = 0.0
                  ELSE
                    scxacionc(mgs,ia,2) = uz(kgs(mgs),2)*cion(mgs,2)*     &
     &                scppx(mgs,ia)/eperao
                  ENDIF
                ELSE  
                scxacionc(mgs,ia,2) = 3.0*pi*xrad(mgs,ia)**2     &
     &            *cion(mgs,2)*Abs(ez(mgs))     &
     &            *uz(kgs(mgs),2)*(1.0 + scppx(mgs,ia)/qm)**2
                  
                ENDIF
              
              ENDIF
              
              IF ( ionatt_conduction_factor /= 1.0 ) THEN
                 scxacionc(mgs,ia,1) = scxacionc(mgs,ia,1)*ionatt_conduction_factor
                 scxacionc(mgs,ia,2) = scxacionc(mgs,ia,2)*ionatt_conduction_factor
              ENDIF
              
             ENDIF ! do_ionatt_conduct
              
            ELSE ! q .le. qxmin
            DO n=1,2
             scxaciond(mgs,ia,n) = 0.0
             scxacionc(mgs,ia,n) = 0.0
            ENDDO ! n
           ENDIF ! q .gt. qxmin
          ENDDO ! ia -- diffusion
          
       
      end do
      
      
!
! load some stuff (may not be needed...)
!
      
      do mgs = 1,ngscnt

!      DO il = lc,lhab
!         conc(igs(mgs),kgs(mgs),il)  = cx(mgs,il)
!      ENDDO
      
      DO n=1,2
      DO il=lc,lhab
        attach(igs(mgs),kgs(mgs),il,n) = cx(mgs,il)*     &
     &     (scxaciond(mgs,il,n) + scxacionc(mgs,il,n)) 
!       IF ( attach(igs(mgs),kgs(mgs),il,n) .ne. 0.0 .and.     &
!      &      ndebug2 .ge. 1  ) THEN
!        write(6,*) 'Attach: ',il,igs(mgs),kgs(mgs),cx(mgs,il),     &
!      &   qx(mgs,il), attach(igs(mgs),kgs(mgs),il,n),n,     &
!      &  scxaciond(mgs,il,n),scxacionc(mgs,il,n)
!       ENDIF
        
!          diff(igs(mgs),kgs(mgs),il,n) = cx(mgs,il)*scxaciond(mgs,il,n)
!          cond(igs(mgs),kgs(mgs),il,n) = cx(mgs,il)*scxacionc(mgs,il,n)

!         IF ( n == 1 ) atotpos = atotpos + attach(igs(mgs),kgs(mgs),il,n)
!         IF ( n == 2 ) atotneg = atotneg + attach(igs(mgs),kgs(mgs),il,n)

      ENDDO
      ENDDO
      
      end do ! mgs

      enddo ! ix

      RETURN
      END

