!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        BALLOON_MODULE        ////////////////////
!
!
!
!===========================================================================
 MODULE BALLOON_MODULE
 
 implicit none
!
! Microphysics Parameters
!
  
      integer, parameter ::  maxsamp = 5000
      integer, parameter ::  numvar  = 16
      integer, parameter ::  maxsnd  = 10
      
!      real, allocatable  :: sound(:,:,:)
      real :: sound(numvar,maxsamp,maxsnd)
      
      integer :: numsnd = 0
      real    :: brise(maxsnd) = 5.0
       
      integer :: istpos(maxsnd)   = 0
      integer :: ijkst(3,maxsnd)  = 0
      real    :: xyzst0(3,maxsnd) = 0.0
      real    :: xyzst(3,maxsnd)  = 0.0
      real    :: tstsnd(maxsnd)   = 0.0
      integer :: icont(maxsnd)    = 0
      character(len=80) :: sndname(maxsnd) = ' '
  
 END MODULE BALLOON_MODULE

!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE BALLOONF                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

! CHECK interpolations for correct array limits!!!!!

      SUBROUTINE BALLOONF        &
       ( ntmul,time,tstop,nstep,nx,ny,nz,na,iunit,  &
        ng,istag,jstag,kstag,                     &
        elec,                                      &
        dtp1,dx,dy,dz,                             &
        dxx,dyy,dzz,                               &
        gxt,gyt,gzt,                               &
        t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,             &
        an,pb,pn,                                  &
        u,v,w)                                     

!     >  (ntmul,nstart,nstop,nstep,nx,ny,nz,na,nba,nv
!     >  ,nor,nht,ngt,imapz,istag,jstag,kstag,itopo,mzdist
!     >  ,neelec,nxelec,nyelec,nzelec
!     >  ,elec
!     >  ,dtp1,dx,dy,dz,ht,gt
!     >  ,nxl,nyl,nzl
!     >  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
!     >  ,ab,ac,ad,an,db,dc,dd,dn
!     <  ,pb,pc,pd,pn,vb,vc,vd,vn)
!
!#######################################################################
!
!     PURPOSE:
!
!     Make a series of virtual balloon EFM flights during a storm
!     simulation.  Behavior is controlled by input file 'inballoon'
!     Sounding data includes charge density, 3D E-field, electric 
!     potential, winds, temperature, pressure, altitude, water substance
!     mixing ratios.
!
!     Assume that the balloon rises/falls at a rate of w + 5 m/s
!
!     Assumes staggered C-grid.
!
!
!   5.13.2002 Updated to use vertically stretched grid
!
!#######################################################################
!

       USE GRID_MODULE
       USE INDEX_MODULE
       USE ELEC_MODULE
       USE BALLOON_MODULE
       USE FILE_MODULE
      
      implicit none
      
!      include 'sam.files.h'

      integer, parameter    :: ng1 = 1
      
      integer ntmul
      real dtp1
      integer iunit ! text IO unit number
      
      integer nstart, tstop, nstep
      integer nx,ny,nz,na,nba,nv
      integer ng,istag,jstag,kstag,itopo
      real dtp,dx,dy,dz
      integer nxl,nyl,nzl

      integer id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      
      integer i,j,k,ix,jy,kz
      integer idebug
      parameter(idebug = 0)
      
      integer ierr
!
! include file for mixing ratio and charge density indices
!
!      include 'sam.index.ion.h'

!
!  Velocity parameters
!
      TYPE(VARIABLE)     :: elec(neelec)

      real               :: u(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real               :: v(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real               :: w(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)


!
! external temporary arrays
!

      real t0(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t1(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t2(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t3(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t4(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t5(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t6(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t7(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t8(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)
      real t9(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng)

      real pb(-ng+1:nz+ng)
      real db(-ng+1:nz+ng)

      real dn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real pn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

      real an(-ng+ng1:nx+ng,-ng+ng1:ny+ng,-ng+ng1:nz+ng,na)

      real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
      real dzz(nz),dxx(nx),dyy(ny)         ! dz(k),dx(i),dy(j)

!
!
!
!  electricity
!
!      integer  neelec,nxelec,nyelec,nzelec
!      integer  iex,iey,iez,iemag,ipot,ie
       integer ie
!      parameter (iex=1,iey=2,iez=3,iemag=4,ipot=5)
!      real elec(nxelec,nyelec,nzelec,neelec)
!      real elec(-ng:nx+ng,-ng:ny+ng,-ng:nz+ng,neelec)
!
!  Interpolation vars
!
      integer ic,jc,kc
      integer iu,ju,ku
      integer iv,jv,kv
      integer iw,jw,kw
      real    eps,facx,facy,facz
      parameter (eps = 1.0e-7)
      real facxu,facyu,faczu
      real facxv,facyv,faczv
      real facxw,facyw,faczw
      integer iwmax,jwmax
      real dumint,telec
      real pbint,pnint   ! pressure 'new' and 'base'
      real uint,vint,wint
      real x,y,z
      
      integer iwmx,jwmx,kwmx
      real wmax   ! max w at k=2
      real emax
      
      real tempsc(2,2,2)
      real scramt
      real qvap
      
!
! sounding
!
      integer llx,lly,llz,llu,llv,llw
      integer llex,lley,llez,llemag,llphi
      integer llsc  ! charge density
      integer llp   ! pressure
      integer llt   ! temp (celsius)
      integer llpt  ! potential temperature (K)
      parameter (llx=1,lly=2,llz=3,llu=4,llv=5,llw=6)
      parameter (llex=7,lley=8,llez=9,llemag=10,llphi=11)
      parameter (llsc=12,llp=13,llt=14,llpt=15)
      integer llqv
      parameter (llqv=15)
      integer lltim
      parameter (lltim=16)
      
      integer ilev(maxsnd)
      save ilev
      data ilev /maxsnd*0/
      
      integer iflight(maxsnd)  ! flight number for a current sounding
      save iflight
      data iflight /maxsnd*0/
      
      integer nflight  ! flight number for next sounding
      save nflight
      
      real times  ! current model time
      integer time
!      real brise  ! balloon rise rate (m/s)
!      parameter(brise = 5.0)
      
! vars to set position of balloon at next timestep
      real nextx(maxsnd),nexty(maxsnd),nextz(maxsnd)
      save nextx,nexty,nextz

!
!  temperature
!
!      t0(ix,jy,kz) = an(ix,jy,kz,lt)
!     >  *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
!

      real poo,rd,cp,cap


!
!  File vars
!
      integer ifile1,ifile2,ifile3
      parameter (ifile1=54,ifile2=55,ifile3=56)
      character*80 fname
      integer iwrite    ! flag to write sounding to diskfile
!      character*ihttima  rstime1, rstime2
      character*6  rstime1, rstime2
      character*2  nmflight
      
!      integer numsnd
!      save numsnd
      
!      integer istpos(maxsnd)
!      integer ijkst(3,maxsnd)
!      real xyzst0(3,maxsnd)
!      real xyzst(3,maxsnd)
!      real tstsnd(maxsnd)
!      integer icont(maxsnd)
!      character*80 sndname

!      save istpos
!      save ijkst
!      save xyzst
!      save tstsnd
!      save icont
      
      integer nevergo
      data nevergo /0/
      
      integer ifirst
      data ifirst /0/

! misc vars
!      integer i,j,k,ix,jy,kz
      integer isnd
           

! #####################################################################
!   Begin Executable Code
! #####################################################################

      IF (numsnd .le. 0) RETURN
!      IF (nevergo .eq. 1) RETURN

      dtp = ntmul*dtp1
      
      poo = 1.0e+05
      rd = 287.04
      cp = 1004.0
      cap = rd/cp 
      iwrite = 0

      times = Float(time) ! dtp1*(nstep-1)

! initialization
      IF ( ifirst .eq. 0) THEN
      
         ifirst = 1
         
         nflight = 1
         
!         DO k=1,maxsnd
!           iflight(k) = k
!         END DO



!        open(unit=ifile3,file='inballoon',status='old',form='formatted')
!        read(ifile3,*) numsnd
        write(iunit,*) 'numsnd = ',numsnd
         IF (numsnd .le. 0) THEN  ! no soundings
            nevergo = 1
            RETURN
         END IF
         IF ( numsnd .gt. maxsnd ) THEN
          write(6,*) 'WARNING! You want ',numsnd,' soundings, ', &
           ' but I can only do ',maxsnd,' at a time.'
           numsnd = maxsnd
         END IF     
        DO k=1,numsnd
           ilev(k) = 0
!           read(ifile3,*) istpos(k)
           write(iunit,*) 'istpos(',k,')= ',istpos(k)
!           read(ifile3,*) ijkst(1,k),ijkst(2,k),ijkst(3,k)
!           read(ifile3,*) xyzst0(1,k),xyzst0(2,k),xyzst0(3,k)
           DO j=1,3
             xyzst0(j,k) = 1000*xyzst0(j,k)
           END DO
           IF (istpos(k) .eq. 1) THEN
           write(iunit,*) ijkst(1,k),ijkst(2,k),ijkst(3,k)
             xyzst0(1,k) = 0.5*dx + (ijkst(1,k)-1)*dx
             xyzst0(2,k) = 0.5*dy + (ijkst(2,k)-1)*dy
!             xyzst0(3,k) = 0.5*dz + (ijkst(3,k)-1)*dz
             xyzst0(3,k) = gzt(ijkst(3,k),1)
           END IF
           
           IF ( istpos(k) .le. 2 ) THEN
             write(iunit,*) 'Ball.sound #',k,' starts at ',xyzst0(1,k),xyzst0(2,k),xyzst0(3,k)
           ENDIF
            xyzst(1,k) = xyzst0(1,k)
            xyzst(2,k) = xyzst0(2,k)
            xyzst(3,k) = xyzst0(3,k)
            nextx(k) = xyzst(1,k)
            nexty(k) = xyzst(2,k) 
            nextz(k) = xyzst(3,k) 
!           read(ifile3,*) tstsnd(k)
!           read(ifile3,*) icont(k)
!           read(ifile3,*) sndname
           IF (icont(k) .eq. 1) THEN
            tstsnd(k) = times - 0.01   ! set start time to current time
!        CALL asnctl ('NEWLOCAL', 1, ierr)
!        CALL asnfile(sndname, '-F f77 -N ieee', ierr)
            open(unit=ifile2,file=sndname(k),form='unformatted')
            read(ifile2) ilev(k)
            call readsound(numvar,ilev(k),ifile2,sound(1,1,k))
            close(ifile2)

!     extrapolate to new balloon position if last time
!     in the sounding file is the current time
            IF (times+0.001 .ge. sound(lltim,ilev(k),k) ) THEN
            xyzst(1,k) = sound(llx,ilev(k),k)
            xyzst(2,k) = sound(lly,ilev(k),k)
            xyzst(3,k) = sound(llz,ilev(k),k)
            nextx(k) = 2.0*xyzst(1,k) - sound(llx,ilev(k)-1,k)
            nexty(k) = 2.0*xyzst(2,k) - sound(lly,ilev(k)-1,k)
            nextz(k) = 2.0*xyzst(3,k) - sound(llz,ilev(k)-1,k)
            ELSE
            DO j=ilev(k)-1,1,-1
              IF ( times+0.001 .ge. sound(lltim,j,k) ) THEN
               xyzst(1,k) = sound(llx,j-1,k)
               xyzst(2,k) = sound(lly,j-1,k)
               xyzst(3,k) = sound(llz,j-1,k)
               nextx(k) = sound(llx,j,k)
               nexty(k) = sound(lly,j,k)
               nextz(k) = sound(llz,j,k)
               write(iunit,*) 'restarting sounding at ', nextx(k),nexty(k),nextz(k)
               ilev(k) = j - 1
               GOTO 100
              END IF
            END DO ! j
            write(6,*) 'hey, this sounding file isnt right!'
            write(6,*) 'Ill just start a new one....'   
            ilev(k) = 0         
  100      CONTINUE 
        
            END IF ! (times .ge. sound(lltim,ilev(k),k)
                 
           END IF ! (icont(k) .eq. 1) 
        
        END DO ! k
        
!        close(ifile3)
      
      END IF ! (nstep .eq. nstart)
      

!
! Loop through all soundings
!      
      DO isnd=1,numsnd
      
      IF (times .lt. tstsnd(isnd) ) GOTO 900
      
        
      IF ( iflight(isnd) .eq. 0 ) THEN
         iflight(isnd) = nflight
         nflight = nflight + 1
      END IF  

      IF ( ilev(isnd) .eq. 0) THEN
        ilev(isnd) = 1
        
        IF ( istpos(isnd) .eq. 3 ) THEN
         IF ( idebug .ge. 1 ) write(6,*) 'here at 1, istpos = 3'
! if start under updraft: find max w in lower half of domain
          wmax = 0.0001
          iwmx = nx/2
          jwmx = ny/2
          kwmx = 2
          DO kz=2,nz/2
          DO jy=2,ny-1
          DO ix=2,nx-1
           IF ( w(ix,jy,kz) .gt. wmax ) THEN
             wmax = w(ix,jy,kz)
             iwmx = ix
             jwmx = jy
             kwmx = 1
           END IF
         END DO
         END DO
         ENDDO
         x = 0.5*dx + (iwmx-1)*dx
         y = 0.5*dy + (jwmx-1)*dy
!         z = 0.5*dz + (kz-1)*dz
         z = gzt(kwmx,1)
        ELSEIF ( istpos(isnd) .eq. 4 ) THEN
         IF ( associated( elec(ipot)%flt3d ) ) THEN
          emax = 0.02
          iwmx = nx/2
          jwmx = ny/2
          kwmx = 1
          DO kz=2,nz-2
          DO jy=2,ny-1
          DO ix=2,nx-1
           IF ( elec(iemag)%flt3d(ix,jy,kz) .gt. emax ) THEN
             emax = elec(iemag)%flt3d(ix,jy,kz) 
             iwmx = ix
             jwmx = jy
             kwmx = kz
           END IF
         ENDDO
         ENDDO
         ENDDO
         x = 0.5*dx + (iwmx-1)*dx
         y = 0.5*dy + (jwmx-1)*dy
!         z = 0.5*dz
         z = gzt(1,1)
         ELSE
           write(0,*) 'Electrification not turned on! Cannot use istpos=4!  Stop!'
         ENDIF
        ELSEIF ( istpos(isnd) .eq. 5 .or. istpos(isnd) .eq. 6 ) THEN
!        write(6,*) 'here at 1'
! if start under updraft: find max w at 2km
          wmax = 0.0001
          iwmx = nx/2
          jwmx = ny/2
          kwmx = 2
          DO kz = 1,nz-1
           IF ( gzt(kz,1) .le. 2000.0 ) THEN
            kwmx = kz
           ENDIF
          ENDDO
          
          DO kz=kwmx,kwmx
          DO jy=2,ny-1
          DO ix=2,nx-1
           IF ( w(ix,jy,kz) .gt. wmax ) THEN
             wmax = w(ix,jy,kz)
             iwmx = ix
             jwmx = jy
           END IF
         END DO
         END DO
         ENDDO
          kwmx = 1
         x = 0.5*dx + (iwmx-1)*dx
         y = 0.5*dy + (jwmx-1)*dy
         IF ( istpos(isnd) .eq. 6 ) THEN
           x = Min(0.5*dx + (nx-2)*dx,Max(1.5*dx, x + xyzst0(1,isnd)))
           y = Min(0.5*dy + (ny-2)*dy,Max(1.5*dy, y + xyzst0(2,isnd)))
         ENDIF
!         z = 0.5*dz + (kz-1)*dz
         z = gzt(kwmx,1)
        
         ELSE 
!         write(6,*) 'here at 2'
         x = xyzst(1,isnd)
         y = xyzst(2,isnd)
         z = xyzst(3,isnd)
         
         END IF
         
        
      ELSE
       ilev(isnd) = ilev(isnd) + 1
       x = nextx(isnd)
       y = nexty(isnd)
       z = nextz(isnd)
      END IF
      
      sound(llx,ilev(isnd),isnd) = x
      sound(lly,ilev(isnd),isnd) = y
      sound(llz,ilev(isnd),isnd) = z

       write(iunit,'(a,i2,a,3(1x,e12.5))') 'Balloon #',iflight(isnd),' is at x,y,z:',x,y,z
      
      sound(lltim,ilev(isnd),isnd) = times

      IF (time .ge. tstop .or. ilev(isnd) .eq. maxsamp) THEN
         iwrite = 1
      ELSE
         iwrite = 0
      END IF

!
!  find lower left corner of scalar grid cell that trajectory
!  end point is currently in for interpolation of scalars
!
      ic = ifix((x+eps)/dx + 0.5)
      jc = ifix((y+eps)/dy + 0.5)
!      kc = ifix((z+eps)/dz + 0.5)
      kc = 1
      DO kz=1,nz-1
        IF ( gzt(kz,1) .le. z ) THEN
          kc = kz
        ENDIF
      ENDDO
      facx = (x/dx - Real(ic) + 0.5)
      facy = (y/dy - Real(jc) + 0.5)
!      facz = (z/dz - Real(kc) + 0.5)
      facz = (z - gzt(kc,1))*gzt(kc,3)
      
      IF ( idebug .ge. 1 ) write(6,*) 'ic,jc,kc = ',ic,jc,kc,facx,facy,facz
      IF ( idebug .ge. 1 ) write(6,*) 'facz: ',z,gzt(kc,1),gzt(kc,3)
      
      IF (kc .eq. nz-2 .or. z .gt. 15.0e3 .or.  &
          ic .le. 0 .or. ic .ge. nx-1 .or.       &
          jc .le. 0 .or. jc .ge. ny-1) THEN
        iwrite = 1
      END IF

!
! interpolate electric field components and potential
!
      IF ( associated( elec(ipot)%flt3d ) ) THEN
      
      do ie = 1,5
      IF (ie .ne. iemag) THEN
      call mlint2   &
      (idebug, elec(ie)%flt3d, -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
        1, 1,nx, ny, nz, telec, facx, facy, facz,ic,jc,kc, 1)
      ELSE
       telec = Sqrt(sound(llex,ilev(isnd),isnd)**2 + &
                    sound(lley,ilev(isnd),isnd)**2 + &
                    sound(llez,ilev(isnd),isnd)**2 )
      END IF
      sound((llex-1+ie),ilev(isnd),isnd) = telec              
      end do
      
      ENDIF

! temperature in celsius:
      call mlint2 &
      (idebug, pb,   1, 1, 1, 1, -ng+1, nz+ng, &
       1, 1, 1, 1, nz, pbint, facx, facy, facz,1,1,kc, 1)

      call mlint2   &
      (idebug, pn,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, pnint, facx, facy, facz,ic,jc,kc, 1)

      call mlint2  &
      (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lt)
      
      sound(llt,ilev(isnd),isnd) = dumint*((pbint + pnint)/poo)**cap
      sound(llp,ilev(isnd),isnd) = pbint + pnint
! substituting qv for potential temp
!      sound(llpt,ilev(isnd),isnd) = dumint
      
      call mlint2  &
      (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, na, nx, ny, nz, qvap, facx, facy, facz,ic,jc,kc, lv)
      
       sound(llqv,ilev(isnd),isnd) = qvap

!
! Note:  Assuming here that t0 has net charge.  Otherwise use tempsc to sum
!        up the charge variables at the 8 points around the balloon -- no need
!        to sum over all points.
!

      DO kz=kc,kc+1
      DO jy=jc,jc+1
      DO ix=ic,ic+1

      tempsc(ix-ic+1,jy-jc+1,kz-kc+1) = t0(ix,jy,kz)
      
      END DO
      END DO
      END DO

      call mlint2 &
      (idebug, tempsc,   1, 2, 1, 2, 1, 2,  &
       1, 1, 2, 2, 2, dumint, facx, facy, facz,1,1,1, 1)
      
      sound(llsc,ilev(isnd),isnd) = dumint


!
!  find lower left corner of u component grid cell 
!
      IF (facx .lt. 0.5) THEN
        iu = ic
        facxu = facx + 0.5
      ELSE
        iu = ic + 1
        facxu = facx - 0.5
      END IF
      
      ju = jc
      ku = kc   
      facyu = facy
      faczu = facz
      
      uint = 0.0

      IF ( idebug .ge. 1 ) write(6,*) 'iu,ju,ku = ',iu,ju,ku,facxu,facyu,faczu,uint

      call mlint2  &
      (idebug, u, -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, uint, facxu, facyu, faczu, iu,ju,ku, 1)

      
      sound(llu,ilev(isnd),isnd) = uint
      
      IF ( idebug .ge. 1 ) write(6,*) 'iu,ju,ku = ',iu,ju,ku,facxu,facyu,faczu,uint
      
!
!  find lower left corner of v component grid cell 
!
      IF (facy .lt. 0.5) THEN
        jv = jc
        facyv = facy + 0.5
      ELSE
        jv = jc + 1
        facyv = facy - 0.5
      END IF
      
      iv = ic
      kv = kc   
      facxv = facx
      faczv = facz      
      
      vint = 0.0

      IF ( idebug .ge. 1 ) write(6,*) 'iv,jv,kv = ',iv,jv,kv,facxv,facyv,faczv,vint

      call mlint2  &
      (idebug, v,  -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, vint, facxv, facyv, faczv, iv,jv,kv, 1)
      
      sound(llv,ilev(isnd),isnd) = vint

      IF ( idebug .ge. 1 ) write(6,*) 'iv,jv,kv = ',iv,jv,kv,facxv,facyv,faczv,vint
      
!
!  find lower left corner of w component grid cell 
!
      IF (facz .lt. 0.5) THEN
        kw = kc
        faczw = facz + 0.5
      ELSE
        kw = kc + 1
        faczw = facz - 0.5
      END IF
      
      iw = ic   
      jw = jc
      facxw = facx
      facyw = facy
      
      call mlint2  &
      (idebug, w, -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, wint, facxw, facyw, faczw, iw,jw,kw, 1)
      
      sound(llw,ilev(isnd),isnd) = wint

      IF ( idebug .ge. 1 ) write(6,*) 'iw,jw,kw = ',iw,jw,kw,facxw,facyw,faczw
      
      nextx(isnd) = x + dtp*uint
      nexty(isnd) = y + dtp*vint
      nextz(isnd) = z + dtp*(wint + brise(isnd))
      
      IF ( idebug .ge. 1 ) write(6,*) 'u,v,w = ',uint,vint,wint,isnd,nextx(isnd),nexty(isnd),nextz(isnd)
      
      IF (iwrite .eq. 1) THEN
      
      iwrite = 0
      
      write(nmflight,915) iflight(isnd)
 915  format(i2.2)
      write(rstime1,timfmt) int(sound(lltim,1,isnd))
      write(rstime2,timfmt) int(sound(lltim,ilev(isnd),isnd))
! 912  format(i5.5)
      
      IF ( lhdrcfs .ge. 1 ) THEN
      fname = hdrcfs(1:lhdrcfs)//filehead(1:lfilehead)//'.bal.'// &
        nmflight//'.'//rstime1(1:ihttima)//'.'//rstime2(1:ihttima)
      ELSE
      fname = filehead(1:lfilehead)//'.bal.'// &
        nmflight//'.'//rstime1(1:ihttima)//'.'//rstime2(1:ihttima)
      ENDIF
      
!        CALL asnctl ('NEWLOCAL', 1, ierr)
!        CALL asnfile(fname, '-F f77 -N ieee', ierr)
      open(unit=ifile1,file=fname,status='unknown',form='unformatted')
      write(ifile1) ilev(isnd)
      call savesound(numvar,ilev(isnd),ifile1,sound(1,1,isnd))
      
      close(ifile1)
      ilev(isnd) = 0  ! setting level to zero signals a new flight
      icont(isnd) = 0  ! sounding finished, so don't continue anymore
      iflight(isnd) = 0  ! a new flight number will be assigned next time
!      iflight = iflight + 1
! set new starting position from inballoon:
            xyzst(1,isnd) = xyzst0(1,isnd)
            xyzst(2,isnd) = xyzst0(2,isnd)
            xyzst(3,isnd) = xyzst0(3,isnd)
            nextx(isnd) = xyzst(1,isnd)
            nexty(isnd) = xyzst(2,isnd) 
            nextz(isnd) = xyzst(3,isnd) 
      END IF ! (iwrite .eq. 1) 
      
 900  CONTINUE
      END DO ! isnd=1,numsnd
      
      RETURN
      END
      
! #####################################################################      
      SUBROUTINE SAVESOUND(numvar,nlev,ifile,sound)
! #####################################################################      

      implicit none

      integer nlev,ifile,numvar
      real sound(numvar,nlev)
      
      write(ifile) sound
      
      RETURN
      END

! #####################################################################      
      SUBROUTINE READSOUND(numvar,nlev,ifile,sound)
! #####################################################################      

      implicit none

      integer nlev,ifile,numvar
      real sound(numvar,nlev)
      
      read(ifile) sound
      
      RETURN
      END


      
