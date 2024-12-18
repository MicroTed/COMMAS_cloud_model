! ##############################################################################
!
!  Subroutine to calculate graupel or rain reflectivity in 2-moment ZXX scheme.
!  Calculation is in a slab (constant jgs)
!
!      call calczgr(nx,ny,nz,nor,na,an,
!     :    z,db1,jgs,ipconc, microp, alpha, lh, ln(lh) )

      subroutine calczgr(nx,ny,nz,nor,na,a,ixe,kze,              &
     &    z,db,jgs,ipconc, alpha, l,ln, qmin, xvmn,xvmx, lvol, rho_qx, ixcol )

      USE INDEX_MODULE
      USE COMMASMPI_MODULE
      USE MICRO_MODULE, only: imurain,imusnow,imorrgdnglimit,morrdnglimit
      
      implicit none

      integer nx,ny,nz,nor,na,ngt,jgs
      integer, parameter :: norz = 3
#ifdef CM1
      real a(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real z(-nor+1:nx+nor,-norz+1:nz+norz,lr:lhab)   ! reflectivity
#else
      real a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real z(-nor+1:nx+nor,-nor+1:nz+nor,lr:lhab)   ! reflectivity
#endif
      real db(nx,nz+1)  ! air density
!      real gt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,ngt)
      
      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer l   ! index for q
      integer ln  ! index for N
      integer lvol ! index for volume
      real    rho_qx
      integer, intent(in) :: ixcol

      integer ix,jy,kz,i1,i2
      real vr,qr,nrx,rd,xv,g1,zx,chw,xdn,xdia1,xdia3
!      real pi
      real :: ynu, cwch, xvmxtmp
      real, parameter ::  pi = 3.14159265 ! 4.0*atan(1.0)
      real, parameter ::  piinv = 1.0/pi
      

!      pi = 4*ATan(1.0)
      
      jy = jgs

      IF ( ixcol > 0 ) THEN
        i1 = ixcol
        i2 = ixcol
      ELSE
        i1 = 1
        i2 = ixe
      ENDIF
      
      IF ( l .eq. lh .or. l .eq. lhl .or. ( l .eq. lr .and. imurain == 1 ) .or. ( l .eq. ls .and. imusnow == 1 ) ) THEN
      
        IF ( l == lh ) cwch = ((3. + alpha)*(2. + alpha)*(1.0 + alpha))**(-1./3.)
        
      DO kz = 1,kze
        DO ix = i1,i2
          
          
          
          IF (  a(ix,jy,kz,l) .gt. qmin .and. a(ix,jy,kz,ln) .gt. 1.e-15 ) THEN
            
            IF ( lvol .gt. 1 ) THEN
                IF ( a(ix,jy,kz,lvol) .gt. 0.0 ) THEN
                  xdn = db(ix,kz)*a(ix,jy,kz,l)/a(ix,jy,kz,lvol)
                  xdn = Min( 900., Max( 170., xdn ) )
                ELSE 
                  xdn = rho_qx
                ENDIF
            ELSE
                xdn = rho_qx
            ENDIF

            IF ( l == lr ) xdn = 1000.

            qr = a(ix,jy,kz,l)
            xv = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
            chw = a(ix,jy,kz,ln)

             IF ( xv .lt. xvmn .or. xv .gt. xvmx ) THEN
              xv = Min( xvmx, Max( xvmn,xv ) )
              chw = db(ix,kz)*a(ix,jy,kz,l)/(xv*xdn)
             ENDIF

         IF (  l == lh .and. ipconc == 5 .and. imorrgdnglimit == 1 ) THEN
           ! limit on characteristic diameter (i.e., 1/slope)
            xdia3 = (xv*6.*piinv)**(1./3.)
            xdia1 = cwch*xdia3
            IF ( xdia1 > morrdnglimit ) THEN
               xdia1 = morrdnglimit
              ! xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
               xv = pi/6.0*(xdia1/cwch)**3
               chw = db(ix,kz)*qr/(xv*xdn)
               a(ix,jy,kz,ln) = chw
               xdia3 = (xv*6.*piinv)**(1./3.)
            ENDIF
         ENDIF


             g1 = (6.0 + alpha)*(5.0 + alpha)*(4.0 + alpha)/  &
     &            ((3.0 + alpha)*(2.0 + alpha)*(1.0 + alpha))
             zx = g1*db(ix,kz)**2*(a(ix,jy,kz,l))*a(ix,jy,kz,l)/chw
!             z(ix,kz,l)  = 1.e18*zx*(6./(pi*1000.))**2
             z(ix,kz,l)  = zx*(6./(pi*1000.))**2
            
          
!          IF ( ny.eq.2 .and. kz .ge. 25 .and. kz .le. 29 .and. z(ix,kz,l) .gt. 0. ) THEN
!             write(*,*) 'calczgr: z,dbz,xdn = ',ix,kz,z(ix,kz,l),10*log10(z(ix,kz,l)),xdn
!          ENDIF
          
          ELSE
           
            z(ix,kz,l) = 0.0
           
          ENDIF
          
        ENDDO
      ENDDO
      
      ELSEIF ( (l == ls .and. imusnow == 3) .or. (l .eq. lr .and. imurain == 3) .or. &
               (l == li ) ) THEN

      xdn = rho_qx ! 1000.
      IF ( l == ls ) ynu = snu
      IF ( l == lr ) ynu = rnu
      IF ( l == li ) ynu = cinu
      
      DO kz = 1,kze
        DO ix = i1,i2
          IF (  a(ix,jy,kz,l) .gt. qmin .and. a(ix,jy,kz,ln) .gt. 1.e-15 ) THEN

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
!            z(ix,kz,l) = 3.6e18*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
            z(ix,kz,l) = 3.6*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
!            qr = a(ix,jy,kz,lr)
!            nrx = a(ix,jy,kz,lnr)
          
          ELSE
           
            z(ix,kz,l) = 0.0
           
          ENDIF
      
          
        ENDDO
      ENDDO
      
      ENDIF
      
      RETURN
      
      END

! ##############################################################################
! ##############################################################################
!
!  Subroutine to correct number concentration to prevent reflectivity growth by 
!  sedimentation in 2-moment ZXX scheme.
!  Calculation is in a slab (constant jgs)
!

      subroutine calcnfromz(nx,ny,nz,nor,na,a,t0,ixe,kze,    &
     &    z0,db,jgs,ipconc, alpha, l,ln, qmin, xvmn,xvmx,t1, &
     &    lvol, rho_qx, infall, ixcol )

      USE INDEX_MODULE
      USE COMMASMPI_MODULE
      USE MICRO_MODULE, only: imurain,imusnow
      
      implicit none

      integer nx,ny,nz,nor,na,ngt,jgs

#ifdef CM1
      integer, parameter :: norz = 3
      real a(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)  ! sedimented N and q
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)    ! sedimented reflectivity
      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)    ! sedimented N (by Vm)
!      real gt(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,ngt)
      real z0(-nor+1:nx+nor,-norz+1:nz+norz,lr:lhab)   ! initial reflectivity
#else
      real a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)  ! sedimented N and q
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)    ! sedimented reflectivity
      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)    ! sedimented N (by Vm)
!      real gt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,ngt)
      real z0(-nor+1:nx+nor,-nor+1:nz+nor,lr:lhab)   ! initial reflectivity
#endif

      real db(nx,nz+1)  ! air density
      
      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer l   ! index for q
      integer ln  ! index for N
      integer lvol ! index for volume
      real    rho_qx
      integer infall
      integer, intent(in) :: ixcol
      
      
      integer ix,jy,kz,i1,i2
      double precision vr,qr,nrx,rd,g1,zx,chw,z,znew,zt,zxt
      real xv,xdn
      real pi
      integer :: ndbz, nmwgt, nnwgt, nwlessthanz
      real :: ynu
      
      ndbz = 0
      nmwgt = 0
      nnwgt = 0
      nwlessthanz = 0
      

      pi = 4*ATan(1.0)
      
      jy = jgs
      
      IF ( ixcol > 0 ) THEN
        i1 = ixcol
        i2 = ixcol
      ELSE
        i1 = 1
        i2 = ixe
      ENDIF
      
      IF ( l .eq. lh .or. l .eq. lhl .or. ( l == lr .and. imurain == 1 )) THEN
      
      
      DO kz = 1,kze
        DO ix = i1,i2

         
          IF (   t0(ix,jy,kz) .gt. 0. ) THEN ! {
!          IF (   a(ix,jy,kz,ln) .gt. 0. ) THEN
            
            IF ( lvol .gt. 1 ) THEN
               IF ( a(ix,jy,kz,lvol) .gt. 0.0 ) THEN
                 xdn = db(ix,kz)*a(ix,jy,kz,l)/a(ix,jy,kz,lvol)
                 xdn = Min( 900., Max( 170., xdn ) )
               ELSE 
                 xdn = rho_qx
               ENDIF
            ELSE
               xdn = rho_qx
            ENDIF
          
            qr = a(ix,jy,kz,l)
            xv = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
            chw = a(ix,jy,kz,ln)

             IF ( xv .lt. xvmn .or. xv .gt. xvmx ) THEN
              xv = Min( xvmx, Max( xvmn,xv ) )
              chw = db(ix,kz)*a(ix,jy,kz,l)/(xv*xdn)
             ENDIF

             g1 = (6.0 + alpha)*(5.0 + alpha)*(4.0 + alpha)/  &
     &            ((3.0 + alpha)*(2.0 + alpha)*(1.0 + alpha))
             zx = g1*db(ix,kz)**2*( a(ix,jy,kz,l))*a(ix,jy,kz,l)/chw
!             z  = 1.e18*zx*(6./(pi*1000.))**2
             z  = zx*(6./(pi*1000.))**2

!          IF ( jy.eq.1 .and. ix .eq. 39  ) THEN
!          IF ( ny.eq.2 .and. kz .ge. 25 .and. kz .le. 29 .and. z .gt. 100. ) THEN
!            write(*,*) 'calcnfromz 1: dbzold/new = ',ix,kz,10*log10(Max(1.e-15,z0(ix,kz,l))),
!     &             10.*log10(t0(ix,jy,kz)),10.*log10(Max(1.e-15,z))
!          ENDIF
            
!           IF ( .true. .or. (z .gt. t0(ix,jy,kz) .and. z .gt. 0.0 .and.
           IF ( (z .gt. t0(ix,jy,kz) .and. z .gt. 0.0 .and.  &
!     &          t0(ix,jy,kz) .gt. 0.0 .and. ! THEN
     &           t0(ix,jy,kz) .gt. z0(ix,kz,l) )) THEN !{

!           IF ( t0(ix,jy,kz) .gt. z0(ix,kz,l) ) THEN
           
!            zx = t0(ix,jy,kz)/(1.e18*(6./(pi*1000.))**2)
            zx = t0(ix,jy,kz)/((6./(pi*1000.))**2)
            
            nrx =  g1*db(ix,kz)**2*( a(ix,jy,kz,l))*a(ix,jy,kz,l)/zx
            IF ( infall .eq. 3 ) THEN
              ! Method I
              IF ( nrx .gt. a(ix,jy,kz,ln) ) THEN
                ndbz = ndbz + 1
                IF ( t1(ix,jy,kz) .lt. ndbz ) nwlessthanz = nwlessthanz + 1
              ELSE
                nnwgt = nnwgt + 1
              ENDIF
              a(ix,jy,kz,ln) = Max( real(nrx), a(ix,jy,kz,ln) )
            ELSE
             ! method I+II
             IF (  nrx .gt. a(ix,jy,kz,ln) .and. t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) THEN
              IF ( nrx .lt. t1(ix,jy,kz)  ) THEN
                ndbz = ndbz + 1
              ELSE
                nmwgt = nmwgt + 1
                IF ( t1(ix,jy,kz) .lt. ndbz ) nwlessthanz = nwlessthanz + 1
              ENDIF
             ELSE
              nnwgt = nnwgt + 1
             ENDIF
              
              a(ix,jy,kz,ln) = Max(Min( real(nrx), t1(ix,jy,kz) ), a(ix,jy,kz,ln) )
            ENDIF


!          IF ( ny.eq.2 .and. kz .ge. 25 .and. kz .le. 29 .and. z .gt. 100. ) THEN
!             zx = g1*db(ix,kz)**2*( 0.224*a(ix,jy,kz,l))*a(ix,jy,kz,l)/nrx
!             znew  = 1.e18*zx*(6./(pi*1000.))**2
!
!             zxt = g1*db(ix,kz)**2*( 0.224*a(ix,jy,kz,l))*a(ix,jy,kz,l)/t1(ix,jy,kz)
!             zt  = 1.e18*zxt*(6./(pi*1000.))**2
!             write(*,*) 'calcnfromz: meth I: dbzold/new = ',ix,kz,10*log10(Max(1.e-15,z)),
!     &             10.*log10(znew),10.*log10(zt),xdn
!             write(*,*) 'calcnfromz: meth I: N old/new = ',ix,kz,chw,nrx,t1(ix,jy,kz),alpha
!           ENDIF

           ELSE ! } {
           
            IF ( t1(ix,jy,kz) .gt. 0 .and. a(ix,jy,kz,ln) .gt. 0 ) THEN
              IF ( t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) THEN
                nmwgt = nmwgt + 1
              ELSE
                nnwgt = nnwgt + 1
              ENDIF
            ENDIF
            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )
            nrx = a(ix,jy,kz,ln)
            


!            IF ( ny.eq.2 .and. kz .ge. 25 .and. kz .le. 29 .and. z .gt. 100. ) THEN
!               write(*,*) 'calcnfromz: method II = ',ix,kz,10*log10(Max(1.e-15,z)),
!     &               10.*log10(t0(ix,jy,kz)),xdn
!               write(*,*) 'calcnfromz: meth II: N old/new = ',ix,kz,chw,nrx,t1(ix,jy,kz),alpha
!             ENDIF
           
           ENDIF ! }
          
           ! }
          ELSE ! {
            IF ( t1(ix,jy,kz) .gt. 0 .and. a(ix,jy,kz,ln) .gt. 0 ) THEN
              IF ( t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) THEN
                nmwgt = nmwgt + 1
              ELSE
                nnwgt = nnwgt + 1
              ENDIF
            ENDIF
!            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )
          ENDIF! }
          
        ENDDO
      ENDDO
      
!      IF ( ny .eq. 2 .and. l .eq. lh ) THEN
!        write(*,'(a,4(1x,i6))') 'ndbz,nmwgt,nnwgt,nwlessthanz = ',ndbz,nmwgt,nnwgt,nwlessthanz
!      ENDIF
      
      ELSEIF ( (l .eq. lr .and. imurain == 3) .or. (l .eq. ls .and. imusnow == 3) ) THEN

      xdn = rho_qx ! 1000.
      IF ( l == ls ) ynu = snu
      IF ( l == lr ) ynu = rnu
      
      DO kz = 1,kze
        DO ix = i1,i2
          IF (  t0(ix,jy,kz) .gt. 0. .and. a(ix,jy,kz,ln) > 1.e-5 ) THEN

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
!            z = 3.6e18*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
            z = 3.6*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
          
             IF ( z .gt. t0(ix,jy,kz) .and. z .gt. 0.0 .and.  &
     &          t0(ix,jy,kz) .gt. 0.0                         &
     &          .and. t0(ix,jy,kz) .gt. z0(ix,kz,l) ) THEN

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn)
            
            chw =  a(ix,jy,kz,ln)
!            nrx =   3.6e18*(ynu+2.0)*vr**2/((ynu+1.0)*t0(ix,jy,kz))
            nrx =   3.6*(ynu+2.0)*vr**2/((ynu+1.0)*t0(ix,jy,kz))
            
            IF ( infall .eq. 3 ) THEN
              a(ix,jy,kz,ln) = Max( real(nrx), a(ix,jy,kz,ln) )
            ELSEIF ( infall .eq. 4 ) THEN
              a(ix,jy,kz,ln) = Max( Min( real(nrx), t1(ix,jy,kz)), a(ix,jy,kz,ln) )
            ENDIF
!            a(ix,jy,kz,ldbzr) = z
!             write(*,*) 'calcnfromz: dbzold/new = ',ix,kz,10*log10(Max(1.e-15,z)),
!     &             10.*log10(znew),10.*log10(zt)
!             write(*,*) 'calcnfromz: Nr old/new = ',ix,kz,chw,nrx,t1(ix,jy,kz)

           ELSE
           
            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )
            
           ENDIF
            
          ELSE
           
            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )
            
          ENDIF
      
          
        ENDDO
      ENDDO
      
      ENDIF
      
      RETURN
      
      END

! ##############################################################################

      subroutine radardd02(nx,ny,nz,nor,na,an,temk,         &
     &    dbz,db,nzdbz,cnoh0t,hwdn1t,ipconc, iunit, microp, printyn, vzflag0, vzf)

      USE INDEX_MODULE
      USE COMMASMPI_MODULE
      USE MICRO_MODULE, only: ithompsoncnoh,imurain, cxmin, iuseferrier,    &
                              idbzci,iusewetgraupel,iusewethail,iusewetsnow

!
! 11.13.2005: Changed values of indices for reordering of lip
!
! 07.13.2005: Fixed an error where cnoh was being used for graupel and frozen drops
!
! 01.24.2005: add ice crystal reflectivity using parameterization of
!             Heymsfield (JAS, 1977).  Could also try Ferrier for this, too.
!
!  09.28.2002 Test alterations for dry ice following Ferrier (1994)
!      for equivalent melted diameter reflectivity.
!      Converted to Fortran by ERM.
!      
!Date: Tue, 21 Nov 2000 10:13:36 -0600 (CST)
!From: Matthew Gilmore <gilmore@hesston.met.tamu.edu>
!
!PRO RF_SPEC ; Computes Radar Reflectivity
!COMMON MAINB, data, x1d, y1d, z1d, iconst, rconst, labels, nx, ny, nz, dshft
!
!;MODIFICATION HISTORY
!; 5/99  -Svelta Veleva introduces variable dielf (const_ki_x) as a (weak)
!;   function of density.  This leads to slight modification of dielf such
!;   that the snow reflectivity is slightly increased - not a big effect.
!;   This is believed to be more accurate than assuming the dielectric
!;   constant for snow is the same as for hail in previous versions.
!
!;On 6/13/99 I added the VIL computation (k=0 in vil array)
!;On 6/15/99 I removed the number concentration dependencies as a function
!;           of temperature (only use for ferrier!)
!;On 6/15/99 I added the Composite reflectivity (k=1 in VIL array)
!;On 6/15/99 I added the Severe Hail Index computation (k=2 in vil array)
!;
!; 6/99 - Veleva and Seo argue that since graupel is more similar to
!;   snow (in number conc and size density) than it is to hail, we
!;   should not weight wetted graupel with the .95 exponent correction
!;   factor as in the case of hail.  An if-statement checks the size
!;   density for wet hail/graupel and treats them appropriately.
!;
!; 6/22/99 - Added function to compute height of max rf and 40 dbz echo top
!;           Also added vilqr which is the model vertical integrated liquid only
!;           using qr.  Will need to check...doesn't seem consistent with vilZ
!;


      implicit none
      
      character(LEN=*) :: microp
      integer nx,ny,nz,nor,na,ngt
      integer nzdbz    !  how many levels actually to process
      
      integer ng1,n10
      integer iunit
      integer :: printyn

      parameter( ng1 = 1, n10 = 23)
      real cno(3:50)
!      real db(nz)
!      real a(nx,ny,nz,na)
      
      real cnoh0t,hwdn1t
      integer ipconc
      real vr
!      parameter ( rnu = -0.8, snu = -0.8 )

!
!   mixing ratio
!
!      integer lv,lc,lr,li,lir,ls
!      integer lgl,lgm,lgh,lf,lh
!      integer lip,lhl,lhab
!      parameter (lv=3,lc=4,lr=5,li=6,lir=7,ls=8)
!      parameter (lgl=9,lgm=10,lgh=11,lf=12,lh=13)
!      parameter (lip=14,lhl=15,lhab=15)
      
!      integer lnr,lns,lnh
!      parameter ( lnr = 44 ) 
      
!      integer lngl,lngm,lngh,lnf
!      integer lnip,lnhl
!      integer lnshi

      integer imapz,mzdist
      
      integer vzflag
      integer, intent(in) :: vzflag0
      
      integer, parameter :: norz = 3
#ifdef CM1
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
      real db(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)  ! air density
!      real gt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,ngt)
      real temk(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)  ! air temperature (kelvin)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)   ! reflectivity
      real vzf(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)   ! power-weighted fall velocity
      real gz(-norz+1:nz+norz) ! ,z1d(-nor+1:nz+nor,4)
#else
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real db(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)  ! air density
!      real gt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,ngt)
      real temk(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)  ! air temperature (kelvin)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)   ! reflectivity
      real, intent(out) :: vzf(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)   ! power-weighted fall velocity
      real gz(-nor+1:nz+nor) ! ,z1d(-nor+1:nz+nor,4)
#endif
      
      real g,cv,cp,rgas,rcp,eta,inveta,rcpinv,cpr,cvr
      real cr1, cr2 ,  hwdnsq,swdnsq
      real rwdnsq, dhmin, qrmin, qsmin, qhmin, qhlmin, tfr, tfrh, zrc
      real reflectmin,  kw_sq
      real const_ki_sn, const_ki_h, ki_sq_sn
      real ki_sq_h, dielf_sn, dielf_h
      real pi
      logical ltest

!  Other data arrays
       real gtmp     (nx,nz)
       real dtmp     (nx,nz)

       real*8 dtmps, dtmpr, dtmph, dtmpf, dtmphl, g1, zx, ze, tmp, x

       integer i,j,k,l,ix,jy,kz,ihcnt

        real*8 xcnoh, xcnos, dadh, dads, zhdryc, zsdryc, zhwetc,zswetc
        real*8 dadr
        real dbzmax,dbzmin
        parameter ( dbzmin = -10 )

      real cnow,cnoi,cnoip,cnoir,cnor,cnos
      real cnogl,cnogm,cnogh,cnof,cnoh,cnohl

      real swdn, rwdn ,hwdn,gldn,gmdn,ghdn,fwdn,hldn
      real swdn0

      real rwdnmx,cwdnmx,cidnmx,xidnmx,swdnmx,gldnmx,gmdnmx
      real ghdnmx,fwdnmx,hwdnmx,hldnmx
      real rwdnmn,cwdnmn,cidnmn,xidnmn,swdnmn,gldnmn,gmdnmn
      real ghdnmn,fwdnmn,hwdnmn,hldnmn
 
      real gldnsq,gmdnsq,ghdnsq,fwdnsq,hldnsq

      real dadgl,dadgm,dadgh,dadhl,dadf
      real zgldryc,zglwetc,zgmdryc, zgmwetc,zghdryc,zghwetc
      real zhldryc,zhlwetc,zfdryc,zfwetc

      real dielf_gl,dielf_gm,dielf_gh,dielf_hl,dielf_fw
      
      integer imx,jmx,kmx
      
      real swdia,gldia,gmdia,ghdia,fwdia,hwdia,hldia
      
      real csw,cgl,cgm,cgh,cfw,chw,chl
      real xvs,xvgl,xvgm,xvgh,xvf,xvh,xvhl
      
!      real xvrmn, xvrmx  ! min, max rain volumes
!      real xvfmn, xvfmx  ! min, max frozen drop volumes
!      real xvsmn, xvsmx  ! min, max snow volumes
!      real xvfmn, xvfmx  ! min, max frozen drop volumes
!      real xvgmn, xvgmx  ! min, max graupel volumes
!      real xvhmn, xvhmx  ! min, max hail volumes
!      real xvhlmn, xvhlmx  ! min, max lg hail volumes
      
!      parameter( xvrmn=2.8866e-13, xvrmx=4.1887e-9 )  ! mks
!      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(1.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvfmn=2.8866e-13, xvfmx=0.523599*(3.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(2.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvsmn=2.8866e-13, xvsmx=0.523599*(3.e-3)**3 ) !( was 4.1887e-9 )  ! mks
!      parameter( xvfmn=2.8866e-13, xvfmx=0.523599*(3.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvgmn=2.8866e-13, xvgmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhmn=0.523599*(1.e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhmn=0.523599*(0.3e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
!      parameter( xvhlmn=0.523599*(10.e-3)**3, 
!     &           xvhlmx=0.523599*(100.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      
      real cwc0
      integer izieg
      integer ice10
      real rhos
      parameter ( rhos = 0.1 )
      
      real qxw, qxw1    ! temp value for liquid water on ice mixing ratio
      real :: dnsnow
      real :: qh
      real, allocatable ::  xvt(:,:,:,:) ! xvt(nx,nz,3,lc:lhab) ! 1=mass-weighted, 2=number-weighted, 3=Z-weighted

!      real, parameter :: cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
!      real, parameter :: cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
!      real, parameter :: cwradn = 5.0e-6     ! minimum radius

      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real cdx(lc:lhab)
      real cwnccn(nz)
      
      real :: vzsnow, vzrain, vzgraupel, vzhail, vzfd
      
      real :: ksq
      
      real :: rdamelt(nch)

      real gsnow1, gsnow53, gsnow73, gamma

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.
#endif

! #########################################################################      

      vzflag = 0
      
#ifndef CM1
      IF ( vzflag0 >= 1 ) THEN
        vzflag = vzflag0
!        write(0,*) 'radardd02: vzflag = ',vzflag
        allocate( xvt(nx,nz+1,3,lc:lhab) )
        xvt(:,:,:,:) = 0.0
        vzf(:,:,:) = 0.0
      ENDIF

#endif

      izieg = 0
      ice10 = 0
      g=9.806                 ! g: gravity constant
      cv=717.0                ! cv: specific heat at constant volume
      cp=1004.0               ! cp: specific heat at constant pressure
      rgas=287.04             ! rgas: gas constant for dry air
      rcp=rgas/cp             ! rcp: gamma constant
      eta=0.622
      inveta = 1./eta
      rcpinv = 1./rcp
      cpr=cp/rgas
      cvr=cv/rgas
      pi = 4.0*ATan(1.)
      cwc0 = 1./pi ! 6.0/pi
      
!      cnoh = 4.0e+04  ! Hail (supercells)
!      hwdn = 900.   ! Hail
      
      cnoh = cnoh0t
      hwdn = hwdn1t

      rwdn = 1000.0
      swdn = 100.0

      qrmin = 1.0e-05
      qsmin = 1.0e-06
      qhmin = 1.0e-05

!
!  default slope intercepts
!
      cnow  = 1.0e+08
      cnoi  = 1.0e+08
      cnoip = 1.0e+08 
      cnoir = 1.0e+08 
      cnor  = 8.0e+06 
      cnos  = 8.0e+06 
      cnogl = 4.0e+05 
      cnogm = 4.0e+05 
      cnogh = 4.0e+05 
      cnof  = 4.0e+05
!      cnoh  = 4.0e+04
      cnohl = 1.0e+03

        gsnow1 = gamma(snu + 1.0)
        gsnow53 = gamma(snu + 5./3.)
        gsnow73 = gamma(snu + 7./3.)


!      cnoh = 4.0e+06 ! GRAUPEL  (Squall lines)
!      hwdn = 400.0   ! GRAUPEL
      

      imx = 1
      jmx = 1
      kmx = 1
!      lv=3
!      lc=4
!      lr=5
      i = 1

      IF ( microp(1:5) .eq. 'ICE10' ) THEN
      
       ice10 = 1

!        print*, 'Set reflectivity for 10ICE'
!         lv=3 - i
!         lc=4 - i
!         lr=5 - i
!         li=6 - i
!         lir=7
!         ls=8
!         lgl=9
!         lgm=10
!         lgh=11
!         lf=12
!         lh=13
!         lip=14
!         lip=7 - i
!         lir=8 - i
!         ls=9 - i
!         lgl=10 - i
!         lgm=11 - i
!         lgh=12 - i
!         lf=13 - i
!         lh=14 - i
!         lhl=15 - i
!         lhab=15 - i
!!         lnr = 44
         
!         cnoh  = 4.0e+04
!         lnshi = 25
!         lnip=3+lnshi
!         lnr=1+lnshi
!         lns=5+lnshi
!         lngl=6+lnshi
!         lngm=7+lnshi
!         lngh=8+lnshi
!         lnf=9+lnshi
!         lnh=10+lnshi
!         lnhl=11+lnshi

         CALL setcno(cno)
         call setqxmin(qxmin)
         CALL setdn10(xdn0)

         cnor  = cno(lr)
         cnos  = cno(ls)
         cnogl = cno(lgl)
         cnogm = cno(lgm)
         cnogh = cno(lgh)
         cnof  = cno(lf)
         cnoh  = cno(lh)
         cnohl = cno(lhl)
         
         
       ELSEIF ( microp(1:8) .eq. 'WARMZIEG' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 

         izieg = 1
         CALL setcnoz(cno)
         call setqxminz(qxmin)

         cnor  = cno(lr)
         qrmin = qxmin(lr)

       ELSEIF ( microp(1:4) .eq. 'ZIEG' .or. microp(1:3) == 'HCM' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 

!        print*, 'Set reflectivity for ZIEG'
         IF ( microp(1:1) .eq. 'Z' ) izieg = 1
         
         IF ( microp(1:3) == 'HCM' ) THEN
           DO l = 1,nch
             rdamelt(l)  = (6.*hm(l)/(pi*1000.))**(1./3.)
           ENDDO
         ENDIF

!         lv = 2
!         lc = 3
!         lr = 4
!         li = 5
!         ls = 6
!         lh = 7
!         lnr= 10
!         lns= 12
!         lnh= 13
         hwdn = hwdn1t ! 500.

         CALL setcnoz(cno)
         call setqxminz(qxmin)

         cnor  = cno(lr)
         cnos  = cno(ls)
         cnoh  = cno(lh)
         qrmin = qxmin(lr)
         qsmin = qxmin(ls)
         qhmin = qxmin(lh)
         IF ( lhl .gt. 1 ) THEN
            cnohl  = cno(lhl)
            qhlmin = qxmin(lhl)
         ENDIF

       ELSEIF ( microp(1:3) .eq. 'ZVD' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 

         izieg = 1
         CALL setcnoz(cno)
         call setqxminz(qxmin)
         
         swdn0 = swdn

         cnor  = cno(lr)
         cnos  = cno(ls)
         cnoh  = cno(lh)
         
         qrmin = qxmin(lr)
         qsmin = qxmin(ls)
         qhmin = qxmin(lh)
         IF ( lhl .gt. 1 ) THEN
            cnohl  = cno(lhl)
            qhlmin = qxmin(lhl)
         ENDIF
!         write(*,*) 'radardbz: ',db(1,1,1),temk(1,1,1),an(1,1,1,lr),an(1,1,1,ls),an(1,1,1,lh)

       ELSE

!         lv=3 - i
!         lc=4 - i
!         lr=5 - i
!         li=6 - i
!         ls=7 - i
!         lh=8 - i
!         lnr=19 - i
!         lns=20 - i
!         lnh=21 - i

        ENDIF


      cdx(lr) = 0.60
      
      IF ( lh > 1 ) THEN
      cdx(lh) = 0.8 ! 1.0 ! 0.45
      cdx(ls) = 2.00
      ENDIF

      IF ( lhl .gt. 1 ) cdx(lhl) = 0.45
      IF ( lf .gt. 1 ) cdx(lf) = 0.45

      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx

      IF ( xvhmn == 0. ) THEN
      
          xvhmn = 0.523599*(dhmn)**3
         IF ( dhmx <= 0.0 ) THEN
           xvhmx = xvhmx0
         ELSE
           xvhmx = 0.523599*(dhmx)**3
         ENDIF
      
      ENDIF
      
      IF ( lh > 1 ) THEN
      xvmn(ls) = xvsmn
      xvmn(lh) = xvhmn
      xvmx(ls) = xvsmx
      xvmx(lh) = xvhmx
      ENDIF

      IF ( lhl .gt. 1 ) THEN
      xvmn(lhl) = xvhlmn
      xvmx(lhl) = xvhlmx
      ENDIF

      IF ( lf .gt. 1 ) THEN
        xvmn(lf) = xvfmn
        xvmx(lf) = xvfmx
      ENDIF

      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      IF ( lh > 1 ) THEN
      xdnmx(li) =  917.0
      xdnmx(ls) =  300.0
      xdnmx(lh) =  900.0
      ENDIF
      IF ( lhl .gt. 1 ) xdnmx(lhl) = 900.0
      IF ( lf .gt. 1 ) xdnmx(lf) = 900.0
!
      xdnmn(:) = 900.0
      
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      IF ( lh > 1 ) THEN
      xdnmn(li) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lh) =  170.0
      ENDIF
      IF ( lhl .gt. 1 ) xdnmn(lhl) = 500.0
      IF ( lf .gt. 1 ) xdnmn(lf) = 150.0

      xdn0(:) = 900.0
      
      xdn0(lc) = 1000.0
      xdn0(lr) = 1000.0
      IF ( lh > 1 ) THEN
      xdn0(li) = 900.0
      xdn0(ls) = 100.0 ! 100.0
      xdn0(lh) = hwdn1t ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      ENDIF
      IF ( lhl .gt. 1 ) xdn0(lhl) = 800.0
      IF ( lf .gt. 1 ) xdn0(lf) = hwdn1t ! 800.0

!
!  slope intercepts
!
!      cnow  = 1.0e+08
!      cnoi  = 1.0e+08
!      cnoip = 1.0e+08 
!      cnoir = 1.0e+08 
!      cnor  = 8.0e+06 
!      cnos  = 8.0e+06 
!      cnogl = 4.0e+05 
!      cnogm = 4.0e+05 
!      cnogh = 4.0e+05 
!      cnof  = 4.0e+05
!c      cnoh  = 4.0e+04
!      cnohl = 1.0e+03
!
!
!  density maximums and minimums
!
      rwdnmx = 1000.0
      cwdnmx = 1000.0
      cidnmx =  917.0
      xidnmx =  917.0
      swdnmx =  200.0
      gldnmx =  400.0
      gmdnmx =  600.0
      ghdnmx =  800.0
      fwdnmx =  900.0
      hwdnmx =  900.0
      hldnmx =  900.0
!
      rwdnmn = 1000.0
      cwdnmn = 1000.0
      xidnmn =  001.0
      cidnmn =  001.0
      swdnmn =  001.0
      gldnmn =  200.0
      gmdnmn =  400.0
      ghdnmn =  600.0
      fwdnmn =  700.0
      hwdnmn =  700.0
!      hldnmn =  900.0

      
      gldn = (0.5)*(gldnmn+gldnmx)  ! 300.
      gmdn = (0.5)*(gmdnmn+gmdnmx)  ! 500.
      ghdn = (0.5)*(ghdnmn+ghdnmx)  ! 700.
      fwdn = (0.5)*(fwdnmn+fwdnmx)  ! 800.
      IF ( ice10 .ge. 1 ) THEN
        hwdn = xdn0(lh) !(0.5)*(hwdnmn+hwdnmx)  ! 800.
        gldn = xdn0(lgl)
        gmdn = xdn0(lgm)
        ghdn = xdn0(lgh)
        fwdn = xdn0(lf)
        hldn = xdn0(lhl)
      ENDIF
      hldn = 900. ! (0.5)*(hldnmn+hldnmx)  ! 900.


      cr1  = 7.2e+20
      cr2  = 7.295e+19
      hwdnsq = hwdn**2
      swdnsq = swdn**2
      rwdnsq = rwdn**2

      gldnsq = gldn**2
      gmdnsq = gmdn**2
      ghdnsq = ghdn**2
      fwdnsq = fwdn**2
      hldnsq = hldn**2
      
      dhmin = 0.005
      tfr   = 273.16
      tfrh  = tfr - 8.0
      zrc   = cr1*cnor
      reflectmin = 0.0
      kw_sq = 0.93
      dbzmax = dbzmin
      
      ihcnt=0

            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Dielectric Factor  - Formulas implemented by Svetla Veleva
!                       following Battan, "Radar Meteorology" - p. 40
!  The result of these calculations is that the dielf numerator (ki_sq) without
!  the density ratio is  .2116 for hail if using 917 density and .25 for
!  snow if using 220 density.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      const_ki_sn = 0.5 - (0.5-0.46)/(917.-220.)*(swdn-220.)
      const_ki_h  = 0.5 - (0.5-0.46)/(917.-220.)*(hwdn-220.)
      ki_sq_sn = (swdnsq/rwdnsq) * const_ki_sn**2
      ki_sq_h  = (hwdnsq/rwdnsq) * const_ki_h**2
      dielf_sn = ki_sq_sn / kw_sq
      dielf_h  = ki_sq_h  / kw_sq
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Use the next line if you want to hardwire dielf for dry hail for both dry
!  snow and dry hail.
!  This would be equivalent to what Straka had originally. (i.e, .21/.93)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dielf_sn = (swdnsq/rwdnsq)*.21/ kw_sq
      dielf_h  = (hwdnsq/rwdnsq)*.21/ kw_sq

      dielf_gl  = (gldnsq/rwdnsq)*.21/ kw_sq
      dielf_gm  = (gmdnsq/rwdnsq)*.21/ kw_sq
      dielf_gh  = (ghdnsq/rwdnsq)*.21/ kw_sq
      dielf_hl  = (hldnsq/rwdnsq)*.21/ kw_sq
      dielf_fw  = (fwdnsq/rwdnsq)*.21/ kw_sq

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Notes on dielectric factors  - from Eun-Kyoung Seo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! constants for both snow and hail would be (x=s,h).....
!       xwdnsq/rwdnsq *0.21/kw_sq   ! Straka/Smith - the original
!       xwdnsq/rwdnsq *0.224        ! Ferrier - for particle sizes in equiv. drop diam
!       xwdnsq/rwdnsq *0.176/kw_sq  ! =0.189 in Smith - for particle sizes in equiv 
!                       ice spheres
!       xwdnsq/rwdnsq *0.208/kw_sq  ! Smith '84 - for particle sizes in equiv melted drop diameter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! VIL algorithm constants
!      Ztop = 10.**(56./10)           !56 dbz is the max rf used by WATADS in cell vil


! Hail detection algorithm constants
!      ZL = 40.
!      ZU = 50.
!      Ho = 3400.  !WATADS Defaults
!      Hm20 = 6200.      !WATADS Defaults

!      DO kz = 1,Min(nzdbz,nz-1)

#ifdef MPI
      kzb = 1
      kze = ktile
      IF ( kzend == nzend ) kze = Max(1,Min(nzdbz,nz-1))
!      if (kzend .le. Max(1,Min(nzdbz,nzend-1))) kze = Max(1,Min(nzdbz,nzend-1)) - kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend - jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend - ixbeg

      DO jy = jyb,jye
#ifndef CM1
        IF ( vzflag == 1 ) THEN
        cwnccn(:) = 1.e9

! a number of arguments are dummy or reasonable values for single-moment droplets/cloud ice, 
! which we don't care about since they're not included for reflectivity.

      call ziegfall(nx,ny,nz,nor,norz,na,1.0,500.,jy,1, & 
     &  xvt, & 
     &  an,db,ipconc,temk,db,1.e9,cwmasn,cwmasx,0.0,1.e8, &  ! substituted db for t7, which doesn't matter, anyway
     &  1.0,1.0,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,2.0e9,xvmn,xvmx,cwnccn, & 
     &  1,1,2)

        ENDIF
#endif
        DO kz = kzb,kze
!         dz = 1./gz(kz)
          DO ix = ixb,ixe
#else
      DO jy=1,Max(1,ny-1)

#ifndef CM1
        IF ( vzflag == 1 ) THEN
        cwnccn(:) = 1.e9

! a number of arguments are dummy or reasonable values for single-moment droplets/cloud ice, 
! which we don't care about since they're not included for reflectivity.

      call ziegfall(nx,ny,nz,nor,norz,na,1.0,500.,jy,1, & 
     &  xvt, & 
     &  an,db,ipconc,temk,db,1.e9,cwmasn,cwmasx,0.0,1.e8, &  ! substituted db for t7, which doesn't matter, anyway
     &  1.0,1.0,cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,2.0e9,xvmn,xvmx,cwnccn, & 
     &  1,1,2)

        ENDIF
#endif
        DO kz = 1,Max(1,Min(nzdbz,nz-1))
!         dz = 1./gz(kz)
!         write(iunit,*) 'kz,db: ',kz,db(1,1,kz)
!         write(iunit,*) 'kz,temk: ',kz,temk(3,3,kz)
         
          DO ix=1,Max(1,nx-1)
#endif
            dbz(ix,jy,kz) = 0.0
                      
          vzsnow = 0.0
          vzrain = 0.0
          vzgraupel = 0.0
          vzhail = 0.0
          vzfd = 0.0
          
          dtmph = 0.0
          dtmpf = 0.0
          dtmps = 0.0
          dtmphl = 0.0
          dtmpr = 0.0
           dadr = (db(ix,jy,kz)/(pi*rwdn*cnor))**(0.25)
!-----------------------------------------------------------------------
! Compute Rain Radar Reflectivity
!-----------------------------------------------------------------------
           
           dtmp(ix,kz) = 0.0
           gtmp(ix,kz) = 0.0
           IF ( an(ix,jy,kz,lr) .ge. qrmin ) THEN
             IF ( ipconc .le. 2 ) THEN
               gtmp(ix,kz) = dadr*an(ix,jy,kz,lr)**(0.25)
               dtmp(ix,kz) = zrc*gtmp(ix,kz)**7
             ELSEIF ( lzr .gt. 1 ) THEN
               dtmp(ix,kz) = 1e18*an(ix,jy,kz,lzr)
             ELSEIF ( an(ix,jy,kz,lnr) .gt. cxmin ) THEN
               IF ( imurain == 3 ) THEN
                 vr = db(ix,jy,kz)*an(ix,jy,kz,lr)/(1000.*an(ix,jy,kz,lnr))
                 dtmp(ix,kz) = 3.647e18*(rnu+2.)*an(ix,jy,kz,lnr)*vr**2/(rnu+1.)
               ELSE ! imurain == 1
                g1 = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
                zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lr))**2/an(ix,jy,kz,lnr)
                ze =1.e18*zx*(6./(pi*1000.))**2 ! note: using 1000. here for water density
                dtmp(ix,kz) = ze
               ENDIF
             ENDIF
             dtmpr = dtmp(ix,kz)
             IF ( vzflag >= 1 ) THEN
               vzrain = xvt(ix,kz,3,lr)*dtmp(ix,kz)
             ENDIF
           ENDIF
           
!-----------------------------------------------------------------------
! Compute snow and graupel reflectivity
!
! Lou modified to look at parcel temperature rather than base state
!-----------------------------------------------------------------------

          IF( lhab .gt. lr ) THEN

!    qs2d   = reform(data[*,*,k,10],[nx*ny])
!    qh2d   = reform(data[*,*,k,11],[nx*ny])

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Only use the following lines if running Straka's GEMS microphysics
!  (Sam 1-d version modified by L Wicker does not use this)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ;xcnoh    = cnoh*exp(-0.025*(temp-tfr))
!    ;xcnos    = cnos*exp(-0.038*(temp-tfr))
!    ;good = where(temp GT tfr, n_elements)
!    ;IF n_elements NE 0 THEN xcnoh(good) = cnoh*exp(-0.075*(temp(good)-tfr))
!    ;IF n_elements NE 0 THEN xcnos(good) = cnos*exp(-0.088*(temp(good)-tfr))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Only use the following lines if running Ferrier micro with No=No(T)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ;  NOSE = -.15
!    ;  NOGE =  .0
!    ;  xcnoh = cnoh*(1.>exp(NOGE*(temp-tfr)) )
!    ;  xcnos = cnos*(1.>exp(NOSE*(temp-tfr)) )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Use the following lines if Nos and Noh are constant
!  (As in Svetla's version of Ferrier, GCE Tao, and SAM 1-d)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xcnoh    = cnoh
        IF ( ipconc < 5 .and. lh > 1 .and. ithompsoncnoh > 0 .and. microp(1:1) == 'Z') THEN
          IF ( an(ix,jy,kz,lh) > qxmin(lh) ) THEN
           xcnoh = Max( 1.e4, Min( 200./an(ix,jy,kz,lh), 5.e6 ) )
          ENDIF
        ENDIF
        xcnos    = cnos

!
! Temporary fix for predicted number concentration -- need a 
! more appropriate reflectivity equation!
!
!        IF ( an(ix,jy,kz,lns) .lt. 0.1 ) THEN
!         swdia = (xvrmn*cwc0)**(1./3.)
!         xcnos = an(ix,jy,kz,ls)*db(ix,jy,kz)/(xvrmn*swdn*swdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         swdia =
!     >  (an(ix,jy,kz,ls)*db(ix,jy,kz)
!     > /(pi*swdn*an(ix,jy,kz,lns)))**(1./3.)
!
!        xcnos = an(ix,jy,kz,lns)/swdia
!        ENDIF

        IF ( ls .gt. 1 ) THEN ! {
        
        IF ( lvs .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lvs) .gt. 0.0 ) THEN
            swdn = db(ix,jy,kz)*an(ix,jy,kz,ls)/an(ix,jy,kz,lvs)
            swdn = Min( 300., Max( 100., swdn ) )
          ELSE 
            swdn = swdn0
          ENDIF
        
        ENDIF 
        
        IF ( ipconc .ge. 5 ) THEN ! {

        xvs = db(ix,jy,kz)*an(ix,jy,kz,ls)/  &
     &      (swdn*Max(1.0e-3,an(ix,jy,kz,lns)))
        IF ( xvs .lt. xvsmn .or. xvs .gt. xvsmx ) THEN
          xvs = Min( xvsmx, Max( xvsmn,xvs ) )
          csw = db(ix,jy,kz)*an(ix,jy,kz,ls)/(xvs*swdn)
        ENDIF

         swdia = (xvs*cwc0)**(1./3.)
         xcnos = an(ix,jy,kz,ls)*db(ix,jy,kz)/(xvs*swdn*swdia)
         
         ENDIF ! }
         ENDIF  ! }

!        IF ( an(ix,jy,kz,lnh) .lt. 0.1 ) THEN
!         hwdia = (xvrmn*cwc0)**(1./3.)
!         xcnoh = an(ix,jy,kz,lh)*db(ix,jy,kz)/(xvrmn*hwdn*hwdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         hwdia =
!     >  (an(ix,jy,kz,lh)*db(ix,jy,kz)
!     > /(pi*hwdn*an(ix,jy,kz,lnh)))**(1./3.)
!        
!         xcnoh = an(ix,jy,kz,lnh)/hwdia
!        ENDIF

        IF ( lh .gt. 1 ) THEN ! {

        IF ( lvh .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lvh) .gt. 0.0 ) THEN
            hwdn = db(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
            hwdn = Min( 900., Max( 170., hwdn ) )
          ELSE 
            hwdn = 500. ! hwdn1t
          ENDIF
        ELSE
          hwdn = hwdn1t
        ENDIF 
        
        IF ( ipconc .ge. 5 ) THEN ! {

        xvh = db(ix,jy,kz)*an(ix,jy,kz,lh)/       &
     &      (hwdn*Max(1.0e-9,an(ix,jy,kz,lnh)))
        IF ( xvh .lt. xvhmn .or. xvh .gt. xvhmx ) THEN
          xvh = Min( xvhmx, Max( xvhmn,xvh ) )
          chw = db(ix,jy,kz)*an(ix,jy,kz,lh)/(xvh*hwdn)
        ENDIF

         hwdia = (xvh*cwc0)**(1./3.)
         xcnoh = an(ix,jy,kz,lh)*db(ix,jy,kz)/(xvh*hwdn*hwdia)
         
        ENDIF ! } ipconc .ge. 5
 
        ENDIF ! }

        dadh = 0.0
        dadhl = 0.0
        dads = 0.0
        IF ( xcnoh .gt. 0.0 ) THEN 
          dadh = ( db(ix,jy,kz) /(pi*hwdn*xcnoh) )**(.25)
          zhdryc = 0.224*cr2*(db(ix,jy,kz)/rwdn)**2/xcnoh ! dielf_h*cr1*xcnoh          ! SV - equiv formula as before but
                                        ! ratio of densities included in
                                        ! dielf_h rather than here following
                                        ! Battan.
        ELSE
          dadh = 0.0
          zhdryc = 0.0
        ENDIF
        
        IF ( xcnos .gt. 0.0 ) THEN
          dads = ( db(ix,jy,kz) /(pi*swdn*xcnos) )**(.25)
          zsdryc = 0.224*cr2*(db(ix,jy,kz)/rwdn)**2/xcnos ! dielf_sn*cr1*xcnos         ! SV - similar change as above
        ELSE
          dads = 0.0
          zsdryc = 0.0
        ENDIF
        zhwetc = zhdryc ! cr1*xcnoh      !Hail/graupel version with .95 power bug removed
        zswetc = zsdryc ! cr1*xcnos
!           
! snow contribution
!
          IF ( ls .gt. 1 ) THEN
          
          gtmp(ix,kz) = 0.0 
          qxw = 0.0 
          qxw1 = 0.0
          dtmps = 0.0
           IF ( an(ix,jy,kz,ls) .ge. qsmin ) THEN !{
            IF ( ipconc .ge. 4 ) THEN  ! (Ferrier 94) !{

             if (lsw .gt. 1) THEN 
               qxw = an(ix,jy,kz,lsw)
               qxw1 = 0.0
             ELSEIF ( Abs(iusewetsnow) == 1 .and. temk(ix,jy,kz) .gt. tfr+1. .and. an(ix,jy,kz,ls) > an(ix,jy,kz,lr) &
     &              .and. an(ix,jy,kz,lr) > qsmin) THEN
               qxw = Min(0.5*an(ix,jy,kz,ls), an(ix,jy,kz,lr))
               qxw1 = qxw
             ENDIF

             vr = xvs ! db(ix,jy,kz)*an(ix,jy,kz,lr)/(1000.*an(ix,jy,kz,lnr))
!             gtmp(ix,kz) = 3.6e18*(0.243*rhos**2/0.93)*(snu+2.)*an(ix,jy,kz,lns)*vr**2/(snu+1.)
             
             ksq = 0.189 ! Smith (1984, JAMC) for equiv. ice sphere
             IF ( an(ix,jy,kz,lns) .gt. 1.e-7 ) THEN
             ksq = 0.224 ! *(xdn0(ls)/917.)**2
!               IF ( .true. ) THEN
               IF ( qxw > qsmin .or. iusewetsnow <= -1) THEN ! old version
!                gtmp(ix,kz) = 3.6e18*(snu+2.)*( 0.224*an(ix,jy,kz,ls) + 0.776*qxw)*an(ix,jy,kz,ls)/ &
!     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2
                gtmp(ix,kz) = 3.6e18*(snu+2.)*( 0.224*(an(ix,jy,kz,ls)+qxw1) + 0.776*qxw)*(an(ix,jy,kz,ls)+qxw1)/ &
     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2

               ELSE ! new form using a mass relationship m = p d^2 (instead of d^3 -- Cox 1988 QJRMS) so that density depends on size
! This one was incorrect but gave a reasonable result thanks to compensating errors
!                 gtmp(ix,kz) = 1.e18* 1.06214**2*(ksq*an(ix,jy,kz,ls) + (1.-ksq)*qxw)*an(ix,jy,kz,ls)*db(ix,jy,kz)**2*gsnow1* gsnow73/    &
!     &                   (an(ix,jy,kz,lns)*(917.)**2* gsnow53**2)
                 dnsnow = 0.0346159*Sqrt(an(ix,jy,kz,lns)/(an(ix,jy,kz,ls)*db(ix,jy,kz)) )
                 IF ( .true. .or. dnsnow < 900. ) THEN 
#if 0
   ! see notebook zrnic1993.nb for details
#endif
                 gtmp(ix,kz) = 1.e18*323.3226* 0.106214**2*(ksq*an(ix,jy,kz,ls) + &
     &                 (1.-ksq)*qxw)*an(ix,jy,kz,ls)*db(ix,jy,kz)**2*gsnow73/    &
     &                   (an(ix,jy,kz,lns)*(917.)**2* gsnow1*(1.0+snu)**(4./3.))
                 ELSE ! otherwise small enough to assume ice spheres?
                gtmp(ix,kz) = (36./pi**2) * 1.e18*(snu+2.)*( 0.224*(an(ix,jy,kz,ls)+qxw1) + 0.776*qxw)*(an(ix,jy,kz,ls)+qxw1)/ &
     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2
                 ENDIF
               
               ENDIF
             ENDIF
             
         !    tmp = Min(1.0,1.e3*(an(ix,jy,kz,ls))*db(ix,jy,kz))
         !    gtmp(ix,kz) = Max( 1.0d0*gtmp(ix,kz), 750.0*(tmp)**1.98)
             dtmps = gtmp(ix,kz)
             dtmp(ix,kz) = dtmp(ix,kz) + gtmp(ix,kz)
            ELSE
             gtmp(ix,kz) = dads*an(ix,jy,kz,ls)**(0.25)
             
             IF ( gtmp(ix,kz) .gt. 0.0 ) THEN !{
             dtmps = zsdryc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                   zsdryc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             ELSE
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zswetc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             ENDIF
             ENDIF !}
            ENDIF !}
           
            IF ( vzflag >= 1 ) THEN
               vzsnow = xvt(ix,kz,3,ls)*dtmps
             ENDIF

           ENDIF !}
           
           ENDIF


         IF ( li .gt. 1 .and. idbzci .ne. 0 ) THEN
          
          
          IF ( idbzci == 1 .and. lni > 0 ) THEN
          ! assume spherical ice with density of 900 for dbz calc
            IF ( an(ix,jy,kz,li) > qxmin(li) .and. an(ix,jy,kz,lni) > 1.0 ) THEN
                 vr = db(ix,jy,kz)*an(ix,jy,kz,li)/(900.*an(ix,jy,kz,lni))
                 dtmp(ix,kz) = dtmp(ix,kz) +  &
     &                 0.224*3.6e18*(cinu+2.)*an(ix,jy,kz,lni)*vr**2/(cinu+1.)*(900./1000.)**2
            ENDIF

          ELSEIF ( idbzci == 2 ) THEN
!
! ice crystal contribution (Heymsfield, 1977, JAS)
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,li) .ge. 0.1e-3 ) THEN
            IF ( ice10 .ge. 1 ) THEN
             gtmp(ix,kz) = Min(1.0,                   &
     &        1.e3*(an(ix,jy,kz,li)+an(ix,jy,kz,lip))*db(ix,jy,kz) )
            ELSE
             gtmp(ix,kz) = Min(1.0,1.e3*(an(ix,jy,kz,li))*db(ix,jy,kz))
            ENDIF
             dtmp(ix,kz) = dtmp(ix,kz) + 750.0*(gtmp(ix,kz))**1.98
           ENDIF

          ELSEIF ( idbzci == 3 .and. lni > 1 ) THEN
             ksq = 0.224 ! *(xdn0(ls)/917.)**2
             qxw1 = 0.0
             qxw = 0.0
            IF ( an(ix,jy,kz,li) > qxmin(li) .and. an(ix,jy,kz,lni) > 1.0 ) THEN
             !  IF ( qxw > qsmin .or. iusewetsnow <= -1) THEN ! old version
!                gtmp(ix,kz) = 3.6e18*(snu+2.)*( 0.224*an(ix,jy,kz,ls) + 0.776*qxw)*an(ix,jy,kz,ls)/ &
!     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2
                gtmp(ix,kz) = 3.6e18*(cinu+2.)*( 0.224*(an(ix,jy,kz,li)+qxw1) + 0.776*qxw)*(an(ix,jy,kz,li)+qxw1)/ &
     &              (an(ix,jy,kz,lni)*(cinu+1.)*rwdn**2)*db(ix,jy,kz)**2

               
             
         !    tmp = Min(1.0,1.e3*(an(ix,jy,kz,ls))*db(ix,jy,kz))
         !    gtmp(ix,kz) = Max( 1.0d0*gtmp(ix,kz), 750.0*(tmp)**1.98)
         !    dtmps = gtmp(ix,kz)
             dtmp(ix,kz) = dtmp(ix,kz) + gtmp(ix,kz)
             ENDIF
           
          
          ENDIF
          
         ENDIF

! frozen drops
           dtmpf = 0.0
         IF ( lf .gt. 1 ) THEN ! {
           gtmp(ix,kz) = 0.0 
           dtmpf = 0.0
           qxw = 0.0

          IF ( izieg .ge. 1 .and. ipconc .ge. 5 ) THEN

           ltest = .false.
           IF ( lzf > 1 ) THEN
             IF ( an(ix,jy,kz,lzf) > 0.0 .and. an(ix,jy,kz,lf) > qhmin .and. &
                  an(ix,jy,kz,lnf) > 0.001*cxmin ) ltest = .true.
           ENDIF
           
           IF ( ltest .or. (an(ix,jy,kz,lf) .ge. qhmin .and. an(ix,jy,kz,lnf) .gt. 0.001*cxmin )) THEN
            
             IF ( lvf .gt. 1 ) THEN
             
              IF ( an(ix,jy,kz,lvf) .gt. 0.0 ) THEN
               hwdn = db(ix,jy,kz)*an(ix,jy,kz,lf)/an(ix,jy,kz,lvf)
               hwdn = Min( 900., Max( 100., hwdn ) )
              ELSE 
               hwdn = 500. ! hwdn1t
              ENDIF

             ENDIF

             chw = an(ix,jy,kz,lnf)
            IF ( chw .gt. 0.0 ) THEN                                         ! (Ferrier 94)
             xvh = db(ix,jy,kz)*an(ix,jy,kz,lf)/(hwdn*Max(1.0e-9,chw))
             IF ( xvh .lt. xvhmn .or. xvh .gt. xvhmx ) THEN
              xvh = Min( xvhmx, Max( xvhmn,xvh ) )
              chw = db(ix,jy,kz)*an(ix,jy,kz,lf)/(xvh*hwdn)
             ENDIF

             qh = an(ix,jy,kz,lf)
             
             IF ( lfw .gt. 1 ) THEN
               IF ( iusewetgraupel .eq. 1 ) THEN
                  qxw = an(ix,jy,kz,lfw)
               ELSEIF ( iusewetgraupel .eq. 2 ) THEN
                  IF ( hwdn .lt. 300. ) THEN
                    qxw = an(ix,jy,kz,lfw)
                  ENDIF
               ENDIF
             ELSEIF ( iusewetgraupel .eq. 3 ) THEN
                  IF ( hwdn .lt. 300. .and. temk(ix,jy,kz) > tfr .and. an(ix,jy,kz,lr) > qhmin ) THEN
                    qxw = Min( an(ix,jy,kz,lf), an(ix,jy,kz,lr))
                    qh = qh + qxw
                  ENDIF
             ELSEIF ( iusewetgraupel == 4 .and. temk(ix,jy,kz) .gt. tfr+0.25 .and. an(ix,jy,kz,lf) > an(ix,jy,kz,lr) &
     &              .and. an(ix,jy,kz,lr) > qhmin) THEN
               qxw = Min(0.5*an(ix,jy,kz,lf), an(ix,jy,kz,lr))
               qh = qh + qxw
             ENDIF
             
              ksq = 0.224 ! *(hwdn/917.)**2
             IF ( lzf .gt. 1 ) THEN
              x = (ksq*qh +  (1.0-ksq)*qxw)/an(ix,jy,kz,lf)  ! weighted average of dielectric const
              dtmpf = 1.e18*x*an(ix,jy,kz,lzf)*(hwdn/rwdn)**2
              dtmp(ix,kz) = dtmp(ix,kz) + dtmpf
             ELSE
             g1 = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))
!             zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lf))**2/chw
!             ze = 0.224*1.e18*zx*(6./(pi*1000.))**2
             zx = g1*db(ix,jy,kz)**2*( ksq*an(ix,jy,kz,lf) + (1.0-ksq)*qxw)*an(ix,jy,kz,lf)/chw
             ze =1.e18*zx*(6./(pi*1000.))**2 ! note: using 1000. here instead of hwdn due to multiplication by (hwdn/rwdn)**2
!             ze =1.e18*zx*(6./(pi))**2 ! note: removed 1000. to go with using ksq adjustment to ice index of refraction
             dtmp(ix,kz) = dtmp(ix,kz) + ze
             dtmpf = ze
             ENDIF
             
            ENDIF
             
        !     IF ( an(ix,jy,kz,lf) .gt. 1.0e-3 ) print*, 'Graupel Z : ',dtmpf,ze
           ENDIF
                    
          ELSE
          
          dtmpf = 0.0
          
           IF ( an(ix,jy,kz,lf) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadh*an(ix,jy,kz,lf)**(0.25)
             IF ( gtmp(ix,kz) .gt. 0.0 ) THEN
             dtmpf =  zhdryc*an(ix,jy,kz,lf)**2/gtmp(ix,kz)
             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lf)**2/gtmp(ix,kz)
             ELSE
!               IF ( hwdn .gt. 700.0 ) THEN
                 dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lf)**2/gtmp(ix,kz)
! 
!     &                               (zhwetc*gtmp(ix,kz)**7)**0.95
!               ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zhwetc*gtmp(ix,kz)**7
!               ENDIF
             ENDIF
             ENDIF
           ENDIF
          
         
          
          ENDIF
 
             IF ( vzflag >= 1 ) THEN
               vzfd = xvt(ix,kz,3,lf)*dtmpf
             ENDIF

          ENDIF ! }         
!           
! graupel/hail contribution
!
         IF ( lh .gt. 1 ) THEN ! {
           gtmp(ix,kz) = 0.0 
           dtmph = 0.0
           qxw = 0.0

          IF ( izieg .ge. 1 .and. ipconc .ge. 5 ) THEN

           ltest = .false.
           IF ( lzh > 1 ) THEN
             IF ( an(ix,jy,kz,lzh) > 0.0 .and. an(ix,jy,kz,lh) > qhmin .and. &
                  an(ix,jy,kz,lnh) > 0.001*cxmin ) ltest = .true.
           ENDIF
           
           IF ( ltest .or. (an(ix,jy,kz,lh) .ge. qhmin .and. an(ix,jy,kz,lnh) .gt. 0.001*cxmin )) THEN
            
             IF ( lvh .gt. 1 ) THEN
             
              IF ( an(ix,jy,kz,lvh) .gt. 0.0 ) THEN
               hwdn = db(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
               hwdn = Min( 900., Max( 100., hwdn ) )
              ELSE 
               hwdn = 500. ! hwdn1t
              ENDIF

             ENDIF

             chw = an(ix,jy,kz,lnh)
            IF ( chw .gt. 0.0 ) THEN                                         ! (Ferrier 94)
             xvh = db(ix,jy,kz)*an(ix,jy,kz,lh)/(hwdn*Max(1.0e-9,chw))
             IF ( xvh .lt. xvhmn .or. xvh .gt. xvhmx ) THEN
              xvh = Min( xvhmx, Max( xvhmn,xvh ) )
              chw = db(ix,jy,kz)*an(ix,jy,kz,lh)/(xvh*hwdn)
             ENDIF

             qh = an(ix,jy,kz,lh)
             
             IF ( lhw .gt. 1 ) THEN
               IF ( iusewetgraupel .eq. 1 ) THEN
                  qxw = an(ix,jy,kz,lhw)
               ELSEIF ( iusewetgraupel .eq. 2 ) THEN
                  IF ( hwdn .lt. 300. ) THEN
                    qxw = an(ix,jy,kz,lhw)
                  ENDIF
               ENDIF
             ELSEIF ( iusewetgraupel .eq. 3 ) THEN
                  IF ( hwdn .lt. 300. .and. temk(ix,jy,kz) > tfr .and. an(ix,jy,kz,lr) > qhmin ) THEN
                    qxw = Min( an(ix,jy,kz,lh), an(ix,jy,kz,lr))
                    qh = qh + qxw
                  ENDIF
             ELSEIF ( iusewetgraupel == 4 .and. temk(ix,jy,kz) .gt. tfr+0.25 .and. an(ix,jy,kz,lh) > an(ix,jy,kz,lr) &
     &              .and. an(ix,jy,kz,lr) > qhmin) THEN
               qxw = Min(0.5*an(ix,jy,kz,lh), an(ix,jy,kz,lr))
               qh = qh + qxw
             ENDIF
             
              ksq = 0.224 ! *(hwdn/917.)**2
             IF ( lzh .gt. 1 ) THEN
              x = (ksq*qh +  (1.0-ksq)*qxw)/an(ix,jy,kz,lh)  ! weighted average of dielectric const
              dtmph = 1.e18*x*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
              dtmp(ix,kz) = dtmp(ix,kz) + dtmph
             ELSE
             g1 = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))
!             zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lh))**2/chw
!             ze = 0.224*1.e18*zx*(6./(pi*1000.))**2
             zx = g1*db(ix,jy,kz)**2*( ksq*an(ix,jy,kz,lh) + (1.0-ksq)*qxw)*an(ix,jy,kz,lh)/chw
             ze =1.e18*zx*(6./(pi*1000.))**2 ! note: using 1000. here instead of hwdn due to multiplication by (hwdn/rwdn)**2
!             ze =1.e18*zx*(6./(pi))**2 ! note: removed 1000. to go with using ksq adjustment to ice index of refraction
             dtmp(ix,kz) = dtmp(ix,kz) + ze
             dtmph = ze
             ENDIF
             
            ENDIF
             
        !     IF ( an(ix,jy,kz,lh) .gt. 1.0e-3 ) print*, 'Graupel Z : ',dtmph,ze
           ENDIF
                    
          ELSE
          
          dtmph = 0.0
          
           IF ( an(ix,jy,kz,lh) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadh*an(ix,jy,kz,lh)**(0.25)
             IF ( gtmp(ix,kz) .gt. 0.0 ) THEN
             dtmph =  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
             ELSE
!               IF ( hwdn .gt. 700.0 ) THEN
                 dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
! 
!     &                               (zhwetc*gtmp(ix,kz)**7)**0.95
!               ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zhwetc*gtmp(ix,kz)**7
!               ENDIF
             ENDIF
             ENDIF
           ENDIF
          
         
          
          ENDIF
 
             IF ( vzflag >= 1 ) THEN
               vzgraupel = xvt(ix,kz,3,lh)*dtmph
             ENDIF

          ENDIF ! }
          
          ENDIF ! na .gt. 5


        
        IF ( izieg .ge. 1 .and. lhl .gt. 1 ) THEN

        hldn = 900.0
        gtmp(ix,kz) = 0.0
        dtmphl = 0.0
        qxw = 0.0
        

        IF ( lvhl .gt. 1 ) THEN
          IF ( an(ix,jy,kz,lvhl) .gt. 0.0 ) THEN
            hldn = db(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
            IF ( lhlw > 1 ) THEN
              hldn = Min( 1000., Max( 300., hldn ) )
            ELSE
              hldn = Min( 900., Max( 300., hldn ) )
            ENDIF
          ELSE 
            hldn = 900. 
          ENDIF
        ELSE
          hldn = 900.
        ENDIF 


        IF ( ipconc .ge. 5 ) THEN

           ltest = .false.
           IF ( lzhl > 1 ) THEN
             IF ( an(ix,jy,kz,lzhl) > 0.0 .and. an(ix,jy,kz,lhl) > qhlmin .and. &
                  an(ix,jy,kz,lnhl) > 0.001*cxmin ) ltest = .true.
           ENDIF

          IF ( ltest .or. ( an(ix,jy,kz,lhl) .ge. qhlmin .and. an(ix,jy,kz,lnhl) .gt. 0.001*cxmin) ) THEN !{
            chl = an(ix,jy,kz,lnhl)
            IF ( chl .gt. 0.0 ) THEN !{
             xvhl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/         &
     &        (hldn*Max(1.0e-9,an(ix,jy,kz,lnhl)))
            IF ( xvhl .lt. xvhlmn .or. xvhl .gt. xvhlmx ) THEN ! {
              xvhl = Min( xvhlmx, Max( xvhlmn,xvhl ) )
              chl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/(xvhl*hldn)
              ! do not update state in dbz calc. ! an(ix,jy,kz,lnhl) = chl
            ENDIF ! }

             qxw = 0
             IF ( lhlw .gt. 1 ) THEN
               IF ( iusewethail .eq. 1 ) THEN
                  qxw = an(ix,jy,kz,lhlw)
               ELSEIF ( iusewethail .eq. 2 ) THEN
                  IF ( hldn .lt. 300. ) THEN
                    qxw = an(ix,jy,kz,lhlw)
                  ENDIF
               ENDIF
             ENDIF
            
              ksq = 0.224 ! *(hldn/917.)**2
             IF ( lzhl .gt. 1 ) THEN !{
              x = (ksq*an(ix,jy,kz,lhl) +  (1.0-ksq)*qxw)/an(ix,jy,kz,lhl)  ! weighted average of dielectric const
              dtmphl = 1.e18*x*an(ix,jy,kz,lzhl)*(hldn/rwdn)**2
              dtmp(ix,kz) = dtmp(ix,kz) + dtmphl
             ELSE !} {

             g1 = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))
!             zx = g1*db(ix,jy,kz)**2*( 0.224*an(ix,jy,kz,lhl) + 0.776*qxw)*an(ix,jy,kz,lhl)/chl
             zx = g1*db(ix,jy,kz)**2*( ksq*an(ix,jy,kz,lhl) + (1.0-ksq)*qxw)*an(ix,jy,kz,lhl)/chl
!             zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lhl))**2/chl
!             ze = 0.224*1.e18*zx*(6./(pi*1000.))**2
             ze = 1.e18*zx*(6./(pi*1000.))**2 ! Argh, had extra factor of 0.224 when already had factor of ksq inside zx
             dtmp(ix,kz) = dtmp(ix,kz) + ze
             dtmphl = ze
             
             ENDIF !}
            ENDIF!}
        !     IF ( an(ix,jy,kz,lh) .gt. 1.0e-3 ) print*, 'Graupel Z : ',dtmph,ze
           ENDIF

          
          ELSE
          
          
           IF ( an(ix,jy,kz,lhl) .ge. qhlmin ) THEN ! {
            dadhl = ( db(ix,jy,kz) /(pi*hldn*cnohl) )**(.25)
             gtmp(ix,kz) = dadhl*an(ix,jy,kz,lhl)**(0.25)
             IF ( gtmp(ix,kz) .gt. 0.0 ) THEN ! {

              zhldryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnohl 

             dtmphl =  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)

             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)
             ELSE
!               IF ( hwdn .gt. 700.0 ) THEN
                 dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)
! 
!     :                               (zhwetc*gtmp(ix,kz)**7)**0.95
!               ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zhwetc*gtmp(ix,kz)**7
!               ENDIF
             ENDIF
             ENDIF ! }
           
           ENDIF ! }
          
         ENDIF ! ipconc .ge. 5

            IF ( vzflag >= 1 ) THEN
               vzhail = xvt(ix,kz,3,lhl)*dtmphl
             ENDIF

        ENDIF ! izieg .ge. 1 .and. lhl .gt. 1 

          IF ( microp(1:3) == 'HCM' ) THEN
          
          dtmph = 0.0
          DO l = 1,nch
            dtmph = dtmph + 0.224*(1000.*rdamelt(l))**6*db(ix,jy,kz)*an(ix,jy,kz,lh-1+l)/hm(l)
          ENDDO
          dtmp(ix,kz) = dtmp(ix,kz) + dtmph
          ENDIF

          
          IF( ice10 .ge. 1 .and. .true. ) THEN

!
! Temporary fix for predicted number concentration -- need a 
! more appropriate reflectivity equation!
!
        IF ( ipconc .ge. 5 ) THEN


!        IF ( an(ix,jy,kz,lngl) .lt. 0.1 ) THEN
!         gldia = (xvrmn*cwc0)**(1./3.)
!         cnogl = an(ix,jy,kz,lgl)*db(ix,jy,kz)/(xvrmn*gldn*gldia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         gldia =
!     >  (an(ix,jy,kz,lgl)*db(ix,jy,kz)
!     > /(pi*gldn*an(ix,jy,kz,lngl)))**(1./3.)
!        
!         cnogl = an(ix,jy,kz,lngl)/gldia
!        ENDIF

        xvgl = db(ix,jy,kz)*an(ix,jy,kz,lgl)/         &
     &      (gldn*Max(1.0e-9,an(ix,jy,kz,lngl)))
        IF ( xvgl .lt. xvgmn .or. xvgl .gt. xvgmx ) THEN
          xvgl = Min( xvgmx, Max( xvgmn,xvgl ) )
          cgl = db(ix,jy,kz)*an(ix,jy,kz,lgl)/(xvgl*gldn)
        ENDIF

         gldia = (xvgl*cwc0)**(1./3.)
         cnogl = an(ix,jy,kz,lgl)*db(ix,jy,kz)/(xvgl*gldn*gldia)


!        IF ( an(ix,jy,kz,lngm) .lt. 0.1 ) THEN
!         gmdia = (xvrmn*cwc0)**(1./3.)
!         cnogm = an(ix,jy,kz,lgm)*db(ix,jy,kz)/(xvrmn*gmdn*gmdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         gmdia =
!     >  (an(ix,jy,kz,lgm)*db(ix,jy,kz)
!     > /(pi*gmdn*an(ix,jy,kz,lngm)))**(1./3.)
!        
!         cnogm = an(ix,jy,kz,lngm)/gmdia
!        ENDIF

        xvgm = db(ix,jy,kz)*an(ix,jy,kz,lgm)/         &
     &      (gmdn*Max(1.0e-9,an(ix,jy,kz,lngm)))
        IF ( xvgm .lt. xvgmn .or. xvgm .gt. xvgmx ) THEN
          xvgm = Min( xvgmx, Max( xvgmn,xvgm ) )
          cgm = db(ix,jy,kz)*an(ix,jy,kz,lgm)/(xvgm*gmdn)
        ENDIF

         gmdia = (xvgm*cwc0)**(1./3.)
         cnogm = an(ix,jy,kz,lgm)*db(ix,jy,kz)/(xvgm*gmdn*gmdia)

!        IF ( an(ix,jy,kz,lngh) .lt. 0.1 ) THEN
!         ghdia = (xvrmn*cwc0)**(1./3.)
!         cnogh = an(ix,jy,kz,lgh)*db(ix,jy,kz)/(xvrmn*ghdn*ghdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         ghdia =
!     >  (an(ix,jy,kz,lgh)*db(ix,jy,kz)
!     > /(pi*ghdn*an(ix,jy,kz,lngh)))**(1./3.)
!        
!         cnogh = an(ix,jy,kz,lngh)/ghdia
!        ENDIF

        xvgh = db(ix,jy,kz)*an(ix,jy,kz,lgh)/         &
     &      (ghdn*Max(1.0e-9,an(ix,jy,kz,lngh)))
        IF ( xvgh .lt. xvgmn .or. xvgh .gt. xvgmx ) THEN
          xvgh = Min( xvgmx, Max( xvgmn,xvgh ) )
          cgh = db(ix,jy,kz)*an(ix,jy,kz,lgh)/(xvgh*ghdn)
        ENDIF

         ghdia = (xvgh*cwc0)**(1./3.)
         cnogh = an(ix,jy,kz,lgh)*db(ix,jy,kz)/(xvgh*ghdn*ghdia)

!        IF ( an(ix,jy,kz,lnf) .lt. 0.1 ) THEN
        xvf = db(ix,jy,kz)*an(ix,jy,kz,lf)/         &
     &      (fwdn*Max(1.0e-9,an(ix,jy,kz,lnf)))
        IF ( xvf .lt. xvfmn .or. xvf .gt. xvfmx ) THEN
          xvf = Min( xvfmx, Max( xvfmn,xvf ) )
          cfw = db(ix,jy,kz)*an(ix,jy,kz,lf)/(xvf*fwdn)
        ENDIF

         fwdia = (xvf*cwc0)**(1./3.)
         cnof = an(ix,jy,kz,lf)*db(ix,jy,kz)/(xvf*fwdn*fwdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         fwdia =
!     >  (an(ix,jy,kz,lf)*db(ix,jy,kz)
!     > /(pi*fwdn*an(ix,jy,kz,lnf)))**(1./3.)
!        
!         cnof = an(ix,jy,kz,lnf)/fwdia
!        ENDIF

!        IF ( an(ix,jy,kz,lnhl) .lt. 0.1 ) THEN
!         hldia = (xvrmn*cwc0)**(1./3.)
!         cnohl = an(ix,jy,kz,lhl)*db(ix,jy,kz)/(xvrmn*hldn*hldia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         hldia =
!     >  (an(ix,jy,kz,lhl)*db(ix,jy,kz)
!     > /(pi*hldn*an(ix,jy,kz,lnhl)))**(1./3.)
!        
!         cnohl = an(ix,jy,kz,lnhl)/hldia
!        ENDIF

        xvhl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/         &
     &      (hldn*Max(1.0e-9,an(ix,jy,kz,lnhl)))
        IF ( xvhl .lt. xvhlmn .or. xvhl .gt. xvhlmx ) THEN
          xvhl = Min( xvhlmx, Max( xvhlmn,xvhl ) )
          chl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/(xvhl*hldn)
        ENDIF

         hldia = (xvhl*cwc0)**(1./3.)
         cnohl = an(ix,jy,kz,lhl)*db(ix,jy,kz)/(xvhl*hldn*hldia)

        ENDIF
!        dadgl = ( db(ix,jy,kz) /(pi*gldn*xcnoh) )**(.25)
!        dadgm = ( db(ix,jy,kz) /(pi*gmdn*xcnoh) )**(.25)
!        dadgh = ( db(ix,jy,kz) /(pi*ghdn*xcnoh) )**(.25)
!        dadhl = ( db(ix,jy,kz) /(pi*hldn*xcnoh) )**(.25)
!        dadf = ( db(ix,jy,kz) /(pi*fwdn*xcnoh) )**(.25)

        dadgl = ( db(ix,jy,kz) /(pi*gldn*cnogl) )**(.25)
        dadgm = ( db(ix,jy,kz) /(pi*gmdn*cnogm) )**(.25)
        dadgh = ( db(ix,jy,kz) /(pi*ghdn*cnogh) )**(.25)
        dadhl = ( db(ix,jy,kz) /(pi*hldn*cnohl) )**(.25)
        dadf  = ( db(ix,jy,kz) /(pi*fwdn*cnof) )**(.25)

        zgldryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnogl 
        zgmdryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnogm 
        zghdryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnogh 
        zhldryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnohl 
        zfdryc  = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnof 
        
!        zgldryc = dielf_gl*cr1*cnogl 
        zglwetc = zgldryc ! cr1*cnogl                

!        zgmdryc = dielf_gm*cr1*cnogm 
        zgmwetc = zgmdryc ! cr1*cnogm               

!        zghdryc = dielf_gh*cr1*cnogh 
        zghwetc = zghdryc ! cr1*cnogh                

!        zhldryc = dielf_hl*cr1*cnohl
        zhlwetc = zhldryc ! cr1*cnohl                

!        zfdryc = dielf_fw*cr1*cnof 
        zfwetc = zfdryc  !  cr1*cnof                

!           
! low dens. graupel contribution
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,lgl) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadgl*an(ix,jy,kz,lgl)**(0.25)
             
!             IF ( temk(ix,jy,kz) .lt. tfr .or.
!     &             temk(ix,jy,kz+1) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zgldryc*an(ix,jy,kz,lgl)**2/gtmp(ix,kz)

!               + zgldryc*gtmp(ix,kz)**7
!             ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zglwetc*gtmp(ix,kz)**7
!             ENDIF
           ENDIF
!           
! med. dens. graupel contribution
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,lgm) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadgm*an(ix,jy,kz,lgm)**(0.25)
             
!             IF ( temk(ix,jy,kz) .lt. tfr .or.
!     &             temk(ix,jy,kz+1) .lt. tfr) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zgmdryc*an(ix,jy,kz,lgm)**2/gtmp(ix,kz)

!            zgmdryc*gtmp(ix,kz)**7
!             ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + 
!     &                  zgmwetc*(gtmp(ix,kz)*(gmdn/900.)**0.25)**7
!             ENDIF
           ENDIF
!           
! high dens. graupel contribution
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,lgh) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadgh*an(ix,jy,kz,lgh)**(0.25)
             
!             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zghdryc*an(ix,jy,kz,lgh)**2/gtmp(ix,kz)

               
!c                zghdryc*gtmp(ix,kz)**7
!             ELSE
!c               IF ( ghdn .gt. 700.0 ) THEN
!                 dtmp(ix,kz) = dtmp(ix,kz) + 
!     &                               (zghwetc*gtmp(ix,kz)**7)**0.95
!c               ELSE
!c                 dtmp(ix,kz) = dtmp(ix,kz) + zghwetc*gtmp(ix,kz)**7
!c               ENDIF
!             ENDIF
           ENDIF
!           
! frozen drop contribution
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,lf) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadf*an(ix,jy,kz,lf)**(0.25)
             
!             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +         &
     &                  zfdryc*an(ix,jy,kz,lf)**2/gtmp(ix,kz)

               
!c                zfdryc*gtmp(ix,kz)**7
!             ELSE
!c               IF ( fwdn .gt. 700.0 ) THEN
!                 dtmp(ix,kz) = dtmp(ix,kz) + 
!     &                               (zfwetc*gtmp(ix,kz)**7)**0.95
!c               ELSE
!c                 dtmp(ix,kz) = dtmp(ix,kz) + zfwetc*gtmp(ix,kz)**7
!c               ENDIF
!             ENDIF
           ENDIF
!           
! hail large contribution
!
          gtmp(ix,kz) = 0.0 
           IF ( an(ix,jy,kz,lhl) .ge. qhmin ) THEN
             gtmp(ix,kz) = dadhl*an(ix,jy,kz,lhl)**(0.25)
             
!             IF ( temk(ix,jy,kz) .lt. tfr ) THEN
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)

               
!c                zhldryc*gtmp(ix,kz)**7
!             ELSE
!c               IF ( hldn .gt. 700.0 ) THEN
!                 dtmp(ix,kz) = dtmp(ix,kz) + 
!     &                               (zhlwetc*gtmp(ix,kz)**7)**0.95
!c               ELSE
!c                 dtmp(ix,kz) = dtmp(ix,kz) + zhlwetc*gtmp(ix,kz)**7
!c               ENDIF
!             ENDIF
           ENDIF
          
          ENDIF ! na .ge. 40
           
          IF ( dtmp(ix,kz) .gt. 0.0 ) THEN
            dbz(ix,jy,kz) = Max(dbzmin, 10.0*Log10(dtmp(ix,kz)) )
            IF ( vzflag >= 1 .and. dbz(ix,jy,kz) > dbzmin ) THEN
              IF ( lf > 0 ) THEN
               vzf(ix,jy,kz) = (vzsnow + vzrain + vzgraupel + vzhail + vzfd)/dtmp(ix,kz)
              ELSE
               vzf(ix,jy,kz) = (vzsnow + vzrain + vzgraupel + vzhail)/dtmp(ix,kz)
              ENDIF

!            IF ( vzf(ix,jy,kz) > 0.2 .and. dbz(ix,jy,kz) > 45.0 ) THEN
!            write(*,*) 'qtovzf: vzf, s,r,h ',vzf(ix,jy,kz),vzsnow, vzrain, vzgraupel, vzhail
!            write(*,*) 'z: s,r,h,0: ', dtmps, dtmpr, dtmph, dtmphl,dtmp(ix,kz)
!            write(*,*) 'q: s,r,h: ',1.e3*an(ix,jy,kz,ls),1.e3*an(ix,jy,kz,lr),1.e3*an(ix,jy,kz,lh)
!            write(*,*) 'dbz: ',dbz(ix,jy,kz)
!            write(*,*) 'xvt2: r,s,h: ',xvt(ix,kz,2,lr),xvt(ix,kz,2,ls),xvt(ix,kz,2,lh)
!            write(*,*) 'xvt1: r,s,h: ',xvt(ix,kz,1,lr),xvt(ix,kz,1,ls),xvt(ix,kz,1,lh)
!            write(*,*) 'xvt3: r,s,h: ',xvt(ix,kz,3,lr),xvt(ix,kz,3,ls),xvt(ix,kz,3,lh)
!            IF ( lhl > 1 ) THEN
!            write(*,*) 'lhl: q,xv2,1,3 = ',1.e3*an(ix,jy,kz,lhl),xvt(ix,kz,2,lhl),xvt(ix,kz,1,lhl),xvt(ix,kz,3,lhl)
!            ENDIF
!            write(*,*)
!            ENDIF
            
            ENDIF
            
            IF ( dbz(ix,jy,kz) .gt. dbzmax ) THEN
              dbzmax = Max(dbzmax,dbz(ix,jy,kz))
              imx = ix
              jmx = jy
              kmx = kz
            ENDIF
          ELSE 
             dbz(ix,jy,kz) = dbzmin
             IF ( lh > 1 .and. lhl > 1) THEN
               IF ( an(ix,jy,kz,lh) > 1.0e-3 .or. an(ix,jy,kz,lr) > 1.0e-3 ) THEN
                 write(0,*) 'radardbz: qr,qh,qhl,ipconc = ',an(ix,jy,kz,lr), an(ix,jy,kz,lh),an(ix,jy,kz,lhl),ipconc
                 write(0,*) 'radardbz: dtmps,dtmph,dadh,dadhl,dtmphl = ',dtmps,dtmph,dadh,dadhl,dtmphl
                 IF ( lnr > 1 .and. lnh > 1 .and. lnhl > 1 ) THEN
                 write(0,*) 'cr,ch,chl = ',an(ix,jy,kz,lnr), an(ix,jy,kz,lnh),an(ix,jy,kz,lnhl)
                 write(0,*) 'lvh,lvhl,hwdn = ',lvh,lvhl,hwdn
                 ENDIF
                 
                 IF ( lzh>1 .and. lzhl>1 ) write(0,*) 'radardbz: zh, zhl = ',an(ix,jy,kz,lzh),an(ix,jy,kz,lzhl)
               ENDIF
             ENDIF
          ENDIF

!         IF ( an(ix,jy,kz,lh) .gt. 1.e-4 .and. 
!     &        dbz(ix,jy,kz) .le. 0.0 ) THEN
!          print*,'dbz = ',dbz(ix,jy,kz)
!          print*,'Hail intercept: ',xcnoh,ix,kz
!          print*,'Hail,snow q: ',an(ix,jy,kz,lh),an(ix,jy,kz,ls)
!          print*,'Hail,snow c: ',an(ix,jy,kz,lnh),an(ix,jy,kz,lns)
!          print*,'dtmps,dtmph = ',dtmps,dtmph
!         ENDIF
        IF ( .not. dtmp(ix,kz) .lt. 1.e30 .or. dbz(ix,jy,kz) > 90.0 ) THEN
!        IF ( ix == 31 .and. kz == 20 .and. jy == 23 ) THEN
          write(0,*) 'my_rank = ',my_rank
          write(0,*) 'ix,jy,kz = ',ix,jy,kz
          write(0,*) 'dbz = ',dbz(ix,jy,kz)
          write(0,*) 'Hail intercept: ',xcnoh,ix,kz
          write(0,*) 'Hail,snow q: ',an(ix,jy,kz,lh),an(ix,jy,kz,ls)
          write(0,*) 'rain q: ',an(ix,jy,kz,lr)
          write(0,*) 'ice q: ',an(ix,jy,kz,li)
          IF ( lhl .gt. 1 ) write(0,*) 'Hail (lhl): ',an(ix,jy,kz,lhl)
          IF (ipconc .ge. 3 ) write(0,*) 'rain c: ',an(ix,jy,kz,lnr)
          IF ( lzr > 1 ) write(0,*) 'rain Z: ',an(ix,jy,kz,lzr)
          IF ( ipconc .ge. 5 ) THEN
          write(0,*) 'Hail,snow c: ',an(ix,jy,kz,lnh),an(ix,jy,kz,lns)
          IF ( lhl .gt. 1 ) write(0,*) 'Hail (lnhl): ',an(ix,jy,kz,lnhl)
          IF ( lzhl .gt. 1 ) THEN 
            write(0,*) 'Hail (lzhl): ',an(ix,jy,kz,lzhl)
            write(0,*) 'chl,xvhl,dhl = ',chl,xvhl,(xvhl*6./3.14159)**(1./3.)
            write(0,*) 'xvhlmn,xvhlmx = ',xvhlmn,xvhlmx
          ENDIF
          ENDIF
          write(0,*) 'chw,xvh = ', chw,xvh
          write(0,*) 'dtmps,dtmph,dadh,dadhl,dtmphl = ',dtmps,dtmph,dadh,dadhl,dtmphl
          write(0,*) 'dtmpr = ',dtmpr
          write(0,*) 'gtmp = ',gtmp(ix,kz),dtmp(ix,kz)
          IF ( .not. (dbz(ix,jy,kz) .gt. -100 .and. dbz(ix,jy,kz) .lt. 200 ) ) THEN
            write(0,*) 'dbz out of bounds! STOP!'
!            STOP
          ENDIF
         ENDIF

           
          ENDDO ! ix
         ENDDO ! jy
      ENDDO ! kz
            
      
      
      
!      print*, 'na,lr = ',na,lr
      IF ( printyn .eq. 1 .and. my_rank == 0 .and. number_of_processes <= 2 ) THEN
!      IF ( dbzmax .gt. dbzmin ) THEN
        write(iunit,*) 'maxdbz,ijk = ',dbzmax,imx,jmx,kmx
        write(iunit,*) 'qrw = ',an(imx,jmx,kmx,lr)
        
        IF ( lh .gt. 1 ) THEN
          write(iunit,*) 'qi  = ',an(imx,jmx,kmx,li)
          write(iunit,*) 'qsw = ',an(imx,jmx,kmx,ls)
          write(iunit,*) 'qhw = ',an(imx,jmx,kmx,lh)
          IF ( lhl .gt. 1 ) write(iunit,*) 'qhl = ',an(imx,jmx,kmx,lhl)
        ENDIF

        IF ( ice10 .ge. 1 ) THEN
          write(iunit,*) 'qgl = ',an(imx,jmx,kmx,lgl)
          write(iunit,*) 'qgm = ',an(imx,jmx,kmx,lgm)
          write(iunit,*) 'qgh = ',an(imx,jmx,kmx,lgh)
          write(iunit,*) 'qhl = ',an(imx,jmx,kmx,lhl)
          write(iunit,*) 'qfw = ',an(imx,jmx,kmx,lf)
        ENDIF
      
      ENDIF
      
      
      IF ( allocated( xvt ) ) deallocate( xvt )
      
      RETURN
      END
      


!-----------------------------------------------------------------------
!
!     ##################################################################
!     ######                                                      ######
!     ######                REAL FUNCTION QTODBZ                  ######
!     ######                                                      ######
!     ##################################################################
!
! Computes effective radar-reflectivity factor corresponding to model's
! hydrometeor variables.
!
! Algorithm is derived from Jerry Straka's microphysics code, which is
! based on the work of Smith (1975).
!
! NOTE:  An important trick to get these calculations to work is to
! use double precision variables for all the intermediate calculations.
! If you dont, you get overflow and underflow error when computing
! the logarithms.  Ugh. 
!
! Units are MKS, and for most accurate results, make sure that the
! number concentrations and densities of rain, snow and hail are the
! same as the model used in producing the fields.
!
! To be fully consistent with the microphysics, hail/graupel particle
! conditions (wet or dry) would need to be passed to this routine.  This
! is because hail/graupel can be wet at temperatures colder than freezing.
! Instead, here we assume that all particles are dry at T<0.  (M. Gilmore)
!
!--------------------------------------------------------------------------
!
! Code obtained from Lou Wicker, 30 August 2004
! Modified by David Dowell, 7 September 2004, after input from Matt Gilmore
!
!
! 2005.07.18:  (erm) Added option for Ferrier (1994) version of dBZ 
!             calculation, which uses equivalent melted diameter.
!             Here it is assumed that all ice particles are dry, which 
!             may not be realistic in that regard.
!
!             Also added dBZ calculation for 10-ice.
!
!--------------------------------------------------------------------------
! INPUTS:
!
!           nc:     number of hydrometeor categories
!
!           q(1):   rainwater mixing ratio (kg/kg) from model
!           q(2):   ice crystal mixing ratio (kg/kg) from model
!           q(3):   snow mixing ratio (kg/kg) from model
!           q(4):   graupel/hail mixing ratio (kg/kg) from model
!
!           pb:     Exner function (total or base-state pressure)
!           tb:     potential temperature (K)
!
!           cnor:   slope intercept of rainwater
!           rho_qr: density of rainwater
!
!           cnos:   slope intercept of snow
!           rho_qs: density of snow
!
!           cnoh:   slope intercept of hail
!           rho_qh: density of hail
!
! OUTPUTS:
!
!           qtodbz: reflectivity (dBZ)
!
!-----------------------------------------------------------------------
      REAL FUNCTION qtodbz(nc, q,                  &
     &                     pb, tb,                 &
     &                     cnor, rho_qr,           &
     &                     cnos, rho_qs,           &
     &                     cnoh, rho_qh,           &
     &                     mindbz, microp,ipconc)

      USE INDEX_MODULE
      USE MICRO_MODULE, only : idbzci, iuseferrier

      implicit none

!---- Passed Variables

      integer nc
      real q(2*lqmx)
      real pb, tb
      real cno(3:50)
      real cnoh, cnos, cnor, rho_qh, rho_qs, rho_qr
      real mindbz
      character(len=*)  :: microp
      integer           :: ipconc
!---  Local Variables

!      integer lv,lc,lr,li,lir,ls
!      integer lgl,lgm,lgh,lf,lh
!      integer lip,lhl,lhab
      

      real*8 qr8, qs8, qh8, qi8, den
      real*8 qgl8, qgm8, qgh8, qf8, qhl8
      real*8 xcnoh, xcnos, zsdryc, zswetc, zhdryc, zhwetc, dadr, dads, dadh
      real*8 cnow,cnoi,cnoip,cnoir
      real*8 cnogl,cnogm,cnogh,cnof,cnohl
      real*8 cr1, cr2, rho_qhsq, rho_qssq, rho_qrsq
      real*8 hwdnsq,swdnsq,rwdnsq,gldnsq,gmdnsq,ghdnsq,fwdnsq,hldnsq
      real*8 swdn,rwdn,hwdn,gldn,gmdn,ghdn,fwdn,hldn
      real*8 dhmin, dielf, pie, qrmin, qsmin, qhmin, tfr, zrc
      real*8 gmh, gms, gmr, zrain, zswet, zhwet, zsdry, zhdry
      real*8 temp, tmp
      real*8 dadgl,dadgm,dadgh,dadhl,dadf
      real*8 zgldryc,zglwetc,zgmdryc, zgmwetc,zghdryc,zghwetc
      real*8 zhldryc,zhlwetc,zfdryc,zfwetc
      real*8 dielf_gl,dielf_gm,dielf_gh,dielf_hl,dielf_fw
      real z0
      
      logical ice10
      
!      integer iuseferrier  ! set = 1 to use alternate dBZ based on Ferrier 1994
                           ! NOTE that snow is treated as always dry.  Might want to
                           ! change that....
!      integer idbzci       ! set = 1 to include dBZ contribution of cloud ice 
                           ! in 3ice (Heymsfield JAS, 1977).  Used for 10-ice by default
!      parameter ( iuseferrier = 1, idbzci = 0 )

      integer nx,ny,nz,nor,na, k
      parameter ( nx = 1, ny = 1, nz = 1, nor = 0, na = 2*lqmx)
      real z1d(1,4), gz(1)
      real an(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real vzf(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real temk(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dbz(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

!-----------------------------------------------------------------------

      cr1  = 7.2d+20
      cr2  = 7.295d+19
      rho_qhsq = rho_qh**2
      rho_qssq = rho_qs**2
      rho_qrsq = rho_qr**2
      dhmin = 0.005
      dielf = 0.21/0.93
      pie   = 4.0*atan(1.0)
      qrmin = 1.0e-10
      qsmin = 1.0e-10
      qhmin = 1.0e-7
      tfr   = 273.16
      zrc   = cr1*cnor
      
      swdn = rho_qs
      rwdn = rho_qr
      hwdn = rho_qh
      hwdnsq = rho_qhsq
      swdnsq = rho_qssq
      rwdnsq = rho_qrsq
      
!      IF ( nc .ge. 12 ) THEN
!        ice10 = .true.
!      ELSE
!        ice10 = .false.
!      ENDIF
      
! Define double precision mixing ratios

      den = 1.0e5*pb**2.509/(287.04*tb)
      temp = tb*pb


      IF ( microp(1:5) .eq. 'ICE10' .or. microp(1:1) .eq. 'Z' .or. microp .eq. 'WARMZIEG' ) THEN
      
      CALL setmicro(cnoh,rho_qh,cnor,rho_qr,cnos,rho_qs)
      
      z1d(:,:) = 1.0
      gz(1) = 1.0
      an(:,:,:,:) = 0.0
      dn = den
      temk = temp
      
      DO k = 1,2*lqmx
        an(1,1,1,k) = q(k)
      ENDDO
      
!      write(*,*) 'qtodbz: den,temp = ',den,temp
      
! assume ipconc = 0 for now....
! assume print unit=6
           call radardd02(nx,ny,nz,nor,na,an,temk,dbz,dn, 1,  &
     &                    cnoh,rho_qh,ipconc, 6, microp, 0, 0, vzf)
         
!         IF ( dbz .gt. 1.0 ) write(*,*) 'qtodbz: dbz = ', dbz
         qtodbz = Max( dbz(1,1,1), mindbz )
         
       RETURN
      
      ENDIF


      xcnoh = cnoh
      xcnos = cnos

      dadh = ( den / (pie*rho_qh*xcnoh) )**.25
      dads = ( den / (pie*rho_qs*xcnos) )**.25
      dadr = ( den / (pie*rho_qr*cnor)  )**.25

      IF ( iuseferrier .eq. 1 ) THEN
        zhdryc = 0.224*cr2*(den/rwdn)**2/xcnoh 
        zsdryc = 0.224*cr2*(den/rwdn)**2/xcnos 
        zhwetc = zhdryc
        zswetc = zsdryc
      ELSE
        zhdryc = dielf*cr1*(rho_qhsq/rho_qrsq)*xcnoh
        zhwetc = cr1*xcnoh
        zsdryc = dielf*cr1*rho_qssq/rho_qrsq*xcnos
        zswetc = cr1*xcnos
      ENDIF

      zrain = 0.
      zswet = 0.
      zsdry = 0.
      zhwet = 0.
      zhdry = 0.
      gmr   = 0.0
      gms   = 0.0
      gmh   = 0.0

!-----------------------------------------------------------------------
! If nc=1, then water only microphysics
!-----------------------------------------------------------------------
      IF( nc.eq.1 ) THEN

        qr8 = q(lr)

!  Slope for rain

        IF ( qr8 .gt. qrmin ) THEN
         IF ( microp .ne. 'SCG'  ) THEN
           gmr = dadr*(qr8)**.25
           zrain = zrc*gmr**7
         ELSE
           zrain = 2.46e4*(1000.*den*qr8)**(1.27)
         ENDIF
        ENDIF

!-----------------------------------------------------------------------
! If nc=4, then rain plus 3 ice categories
!-----------------------------------------------------------------------
      ELSEIF (nc.eq.4 ) THEN

        qr8 = q(lr)
        qs8 = q(ls)
        qh8 = q(lh)
        
        IF ( idbzci .ne. 0  ) THEN
            qi8 = q(li)
          IF ( qi8 .gt. qsmin ) THEN
            tmp = Min(1.5d0,1.d3*qi8*den)
            zsdry = 750.0*tmp**1.98
          ENDIF
        ENDIF

!  Slope for rain

        IF ( qr8 .gt. qrmin ) THEN
         gmr = dadr*(qr8)**.25 
         zrain = zrc*gmr**7
        ENDIF

!  Slope for snow

        IF ( qs8 .gt. qsmin ) THEN
         gms = dads*(qs8)**.25

!  Computation for dry and wet snow reflectivity

        IF ( temp .lt. tfr ) THEN
          IF ( iuseferrier .eq. 1  ) THEN
            zsdry = zsdry + zsdryc*qs8**2/gms
          ELSE
            zsdry = zsdry + zsdryc*gms**7
          ENDIF
        ELSE
          IF ( iuseferrier .eq. 1 ) THEN
            zswet = zswetc*qs8**2/gms
          ELSE
            zswet = zswetc*gms**7
          ENDIF
        END IF
        
        ENDIF ! ( qs8 .gt. qsmin )

!  Slope for hail/graupel
 
        IF ( qh8 .gt. qhmin ) THEN
        
          gmh = dadh*(qh8)**.25

!  Computation for dry and wet hail reflectivity including mie scattering
!  via parameterization (Smith, 1975) 

        IF ( temp .lt. tfr ) THEN
          IF ( iuseferrier .eq. 1 ) THEN
            zhdry = zhdryc*qh8**2/gmh
          ELSE
            zhdry = zhdryc*gmh**7
          ENDIF
        ELSE
          IF ( iuseferrier .eq. 1 ) THEN
            zhwet = zhwetc*qh8**2/gmh
          ELSE
            zhwet = (zhwetc*(gmh**7))**.95
          ENDIF
        ENDIF
        
        ENDIF ! ( qh8 .gt. qhmin )
      
      
      ELSE ! ( nc value not recognized )

        write(*,*) 'qtodbz: *** ERROR *** nc =', nc
        stop
      
      ENDIF

 

!-----------------------------------------------------------------------
! Compute final reflectivity value
!-----------------------------------------------------------------------

      z0 = zrain + zswet + zsdry + zhwet + zhdry

      IF( z0 .gt. 0. ) THEN
        qtodbz = max(mindbz, 10.0 * log10(z0))
      ELSE
        qtodbz = mindbz
      ENDIF

      RETURN
      END FUNCTION QTODBZ

!
! ##############################################################################
!

!-----------------------------------------------------------------------
!
!     ##################################################################
!     ######                                                      ######
!     ######                REAL FUNCTION QTODBZ                  ######
!     ######                                                      ######
!     ##################################################################
!
! Computes effective radar-reflectivity factor corresponding to model's
! hydrometeor variables.
!
! Algorithm is derived from Jerry Straka's microphysics code, which is
! based on the work of Smith (1975).
!
! NOTE:  An important trick to get these calculations to work is to
! use double precision variables for all the intermediate calculations.
! If you dont, you get overflow and underflow error when computing
! the logarithms.  Ugh. 
!
! Units are MKS, and for most accurate results, make sure that the
! number concentrations and densities of rain, snow and hail are the
! same as the model used in producing the fields.
!
! To be fully consistent with the microphysics, hail/graupel particle
! conditions (wet or dry) would need to be passed to this routine.  This
! is because hail/graupel can be wet at temperatures colder than freezing.
! Instead, here we assume that all particles are dry at T<0.  (M. Gilmore)
!
!--------------------------------------------------------------------------
!
! Code obtained from Lou Wicker, 30 August 2004
! Modified by David Dowell, 7 September 2004, after input from Matt Gilmore
!
!
! 2005.07.18:  (erm) Added option for Ferrier (1994) version of dBZ 
!             calculation, which uses equivalent melted diameter.
!             Here it is assumed that all ice particles are dry, which 
!             may not be realistic in that regard.
!
!             Also added dBZ calculation for 10-ice.
!
!--------------------------------------------------------------------------
! INPUTS:
!
!           nc:     number of hydrometeor categories
!
!           q(1):   rainwater mixing ratio (kg/kg) from model
!           q(2):   ice crystal mixing ratio (kg/kg) from model
!           q(3):   snow mixing ratio (kg/kg) from model
!           q(4):   graupel/hail mixing ratio (kg/kg) from model
!
!           pb:     Exner function (total or base-state pressure)
!           tb:     potential temperature (K)
!
!           cnor:   slope intercept of rainwater
!           rho_qr: density of rainwater
!
!           cnos:   slope intercept of snow
!           rho_qs: density of snow
!
!           cnoh:   slope intercept of hail
!           rho_qh: density of hail
!
! OUTPUTS:
!
!           qtodbz: reflectivity (dBZ)
!
!-----------------------------------------------------------------------
      SUBROUTINE    QTOVZF(nc, q,                  &
     &                     pb, tb,                 &
     &                     cnor, rho_qr,           &
     &                     cnos, rho_qs,           &
     &                     cnoh, rho_qh,           &
     &                     mindbz, microp,ipconc,qtodbz,vzf)

      USE INDEX_MODULE
      USE MICRO_MODULE, only : idbzci, iuseferrier

      implicit none

!---- Passed Variables

      integer nc
      real q(2*lqmx)
      real pb, tb
      real cno(3:50)
      real cnoh, cnos, cnor, rho_qh, rho_qs, rho_qr
      real mindbz
      character(len=*)  :: microp
      integer           :: ipconc
      real              :: qtodbz,vzf
!---  Local Variables

!      integer lv,lc,lr,li,lir,ls
!      integer lgl,lgm,lgh,lf,lh
!      integer lip,lhl,lhab
      

      real*8 qr8, qs8, qh8, qi8, den
      real*8 qgl8, qgm8, qgh8, qf8, qhl8
      real*8 xcnoh, xcnos, zsdryc, zswetc, zhdryc, zhwetc, dadr, dads, dadh
      real*8 zhdrycsmith,zsdrycsmith
      real*8 cnow,cnoi,cnoip,cnoir
      real*8 cnogl,cnogm,cnogh,cnof,cnohl
      real*8 cr1, cr2, rho_qhsq, rho_qssq, rho_qrsq
      real*8 hwdnsq,swdnsq,rwdnsq,gldnsq,gmdnsq,ghdnsq,fwdnsq,hldnsq
      real*8 swdn,rwdn,hwdn,gldn,gmdn,ghdn,fwdn,hldn
      real*8 dhmin, dielf, pie, qrmin, qsmin, qhmin, tfr, zrc
      real*8 gmh, gms, gmr, zrain, zswet, zhwet, zsdry, zhdry
      real*8 temp, tmp
      real*8 dadgl,dadgm,dadgh,dadhl,dadf
      real*8 zgldryc,zglwetc,zgmdryc, zgmwetc,zghdryc,zghwetc
      real*8 zhldryc,zhlwetc,zfdryc,zfwetc
      real*8 dielf_gl,dielf_gm,dielf_gh,dielf_hl,dielf_fw
      real*8 vzsnow, vzrain, vzhail
      real z0
      
      logical ice10
      
!      integer iuseferrier  ! set = 1 to use alternate dBZ based on Ferrier 1994
                           ! NOTE that snow is treated as always dry.  Might want to
                           ! change that....
!      integer idbzci       ! set = 1 to include dBZ contribution of cloud ice 
                           ! in 3ice (Heymsfield JAS, 1977).  Used for 10-ice by default
!      parameter ( iuseferrier = 1, idbzci = 0 )

      integer nx,ny,nz,nor,na, k
      parameter ( nx = 1, ny = 1, nz = 1, nor = 0, na = 2*lqmx)
      real z1d(1,4), gz(1)
      real an(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real zv(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real temk(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real dbz(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ar,br,cs,ds
      real :: gf7br = -1.0
      real :: gf7ds = -1.0
      real :: gf7p5 = -1.0
      real :: gf4ds = -1.0
      save gf7br, gf7ds, gf7p5, gf4ds
      real :: gamma
      real, parameter :: rho00 = 1.204 ! 1.225
      real, parameter :: gr = 9.8
      real  ::  cd

!-----------------------------------------------------------------------

      cr1  = 7.2d+20
      cr2  = 7.295d+19
      rho_qhsq = rho_qh**2
      rho_qssq = rho_qs**2
      rho_qrsq = rho_qr**2
      dhmin = 0.005
      dielf = 0.21/0.93
      pie   = 4.0*atan(1.0)
      qrmin = 1.0e-10
      qsmin = 1.0e-10
      qhmin = 1.0e-7
      tfr   = 273.16
      zrc   = cr1*cnor
      cd = 0.6
      
      swdn = rho_qs
      rwdn = rho_qr
      hwdn = rho_qh
      hwdnsq = rho_qhsq
      swdnsq = rho_qssq
      rwdnsq = rho_qrsq
      
      vzsnow = 0.0d0
      vzrain = 0.0d0
      vzhail = 0.0d0
      vzf    = 0.0

      
!      IF ( nc .ge. 12 ) THEN
!        ice10 = .true.
!      ELSE
!        ice10 = .false.
!      ENDIF
      
! Define double precision mixing ratios

      den = 1.0e5*pb**2.509/(287.04*tb)
      temp = tb*pb


      IF ( microp(1:5) .eq. 'ICE10' .or. microp(1:1) .eq. 'Z' .or. microp .eq. 'WARMZIEG' ) THEN
      
      CALL setmicro(cnoh,rho_qh,cnor,rho_qr,cnos,rho_qs)
      
      z1d(:,:) = 1.0
      gz(1) = 1.0
      an(:,:,:,:) = 0.0
      dn = den
      temk = temp
      
      DO k = 1,2*lqmx
        an(1,1,1,k) = q(k)
      ENDDO
      
!      write(*,*) 'qtodbz: den,temp = ',den,temp
      
! assume ipconc = 0 for now....
! assume print unit=6
           call radardd02(nx,ny,nz,nor,na,an,temk,dbz,dn, 1,  &
     &                    cnoh,rho_qh,ipconc, 6, microp, 0, 1, zv)
         
!         IF ( dbz .gt. 1.0 ) write(*,*) 'qtodbz: dbz = ', dbz
         qtodbz = Max( dbz(1,1,1), mindbz )
         vzf = zv(1,1,1)
       RETURN
      
      ENDIF


      xcnoh = cnoh
      xcnos = cnos

      dadh = ( den / (pie*rho_qh*xcnoh) )**.25
      dads = ( den / (pie*rho_qs*xcnos) )**.25
      dadr = ( den / (pie*rho_qr*cnor)  )**.25

      IF ( iuseferrier .eq. 1 ) THEN
        zhdryc = 0.224*cr2*(den/rwdn)**2/xcnoh 
        zsdryc = 0.224*cr2*(den/rwdn)**2/xcnos 
        zhwetc = zhdryc
        zswetc = zsdryc
        zhdrycsmith = dielf*cr1*(rho_qhsq/rho_qrsq)*xcnoh
        zsdrycsmith = dielf*cr1*rho_qssq/rho_qrsq*xcnos
      ELSE
        zhdryc = dielf*cr1*(rho_qhsq/rho_qrsq)*xcnoh
        zhwetc = cr1*xcnoh
        zsdryc = dielf*cr1*rho_qssq/rho_qrsq*xcnos
        zswetc = cr1*xcnos
      ENDIF

      zrain = 0.
      zswet = 0.
      zsdry = 0.
      zhwet = 0.
      zhdry = 0.
      gmr   = 0.0
      gms   = 0.0
      gmh   = 0.0

      ar    = 841.99666        ! In mks units....(identical to LFO83 of 2115 in cgs)
      br    = 0.8
      cs    = 4.83607122
      ds    = 0.25  
!      IF ( gf7br <= 0.0 ) gf7br = 3376.92 !gamma(7.0+br)
!      IF ( gf7ds <= 0.0 ) gf7ds = 1155.38 !  gamma(7.0+ds)
!      IF ( gf7p5 <= 0.0 ) gf7p5 = 1871.25 ! gamma(7.5)
!      IF ( gf4ds <= 0.0 ) gf4ds = 8.28509 ! gamma(4+ds)

      IF ( gf7br <= 0.0 ) gf7br = gamma(7.0+br)
      IF ( gf7ds <= 0.0 ) gf7ds = gamma(7.0+ds)
      IF ( gf7p5 <= 0.0 ) gf7p5 = gamma(7.5)
      IF ( gf4ds <= 0.0 ) gf4ds = gamma(4+ds)

!-----------------------------------------------------------------------
! If nc=1, then water only microphysics
!-----------------------------------------------------------------------
      IF( nc.eq.1 ) THEN

        qr8 = q(lr)

!  Slope for rain

        IF ( qr8 .gt. qrmin ) THEN
         IF ( microp .ne. 'SCG'  ) THEN
           gmr = dadr*(qr8)**.25
           zrain = zrc*gmr**7
           vzrain = 1.e18*cnor*ar*Sqrt(rho00/den)*gf7br*((den*qr8)/(pie*cnor*rwdn))**(0.25*(br+7.0))
         ELSE
           zrain = 2.46e4*(1000.*den*qr8)**(1.27)
         ENDIF
        ENDIF

!-----------------------------------------------------------------------
! If nc=4, then rain plus 3 ice categories
!-----------------------------------------------------------------------
      ELSEIF (nc.eq.4 ) THEN

        qr8 = q(lr)
        qs8 = q(ls)
        qh8 = q(lh)
        
        IF ( idbzci .ne. 0  ) THEN
            qi8 = q(li)
          IF ( qi8 .gt. qsmin ) THEN
            tmp = Min(1.5d0,1.d3*qi8*den)
            zsdry = 750.0*tmp**1.98
          ENDIF
        ENDIF

!  Slope for rain

        IF ( qr8 .gt. qrmin ) THEN
         gmr = dadr*(qr8)**.25 
         zrain = zrc*gmr**7
        ENDIF

!  Slope for snow

        IF ( qs8 .gt. qsmin ) THEN
         gms = dads*(qs8)**.25

!  Computation for dry and wet snow reflectivity

        IF ( temp .lt. tfr ) THEN
          IF ( iuseferrier .eq. 1  ) THEN
            zsdry = zsdry + zsdryc*qs8**2/gms
            vzsnow = (zsdryc*qs8**2/gms)*(cs*gf4ds/6.0)*(gms**ds)*Sqrt(rho00/den)
!           vzsnow = 1.e18*cnos*cs*Sqrt(rho00/den)*gf7ds*((den*qs8)/(pie*cnos*swdn))*(0.25*(ds+7.0))
          ELSE
            zsdry = zsdry + zsdryc*gms**7
           vzsnow = 1.e18*dielf*rho_qssq/rho_qrsq*cnos*cs*Sqrt(rho00/den)*gf7ds*((den*qs8)/(pie*cnos*swdn))**(0.25*(ds+7.0))
          ENDIF
        ELSE
          IF ( iuseferrier .eq. 1 ) THEN
            zswet = zswetc*qs8**2/gms
            vzsnow = (zsdryc*qs8**2/gms)*(cs*gf4ds/6.0)*(gms**ds)*Sqrt(rho00/den)
          ELSE
            zswet = zswetc*gms**7
            vzsnow = 1.e18*cnos*cs*Sqrt(rho00/den)*gf7ds*((den*qs8)/(pie*cnos*swdn))**(0.25*(ds+7.0))
          ENDIF
        END IF
        
        ENDIF ! ( qs8 .gt. qsmin )

!  Slope for hail/graupel
 
        IF ( qh8 .gt. qhmin ) THEN
        
          gmh = dadh*(qh8)**.25

!  Computation for dry and wet hail reflectivity including mie scattering
!  via parameterization (Smith, 1975) 

        IF ( temp .lt. tfr .or. iuseferrier .eq. 1 ) THEN
          IF ( iuseferrier .eq. 1 ) THEN
            zhdry = zhdryc*qh8**2/gmh
            tmp = (zhdrycsmith*gmh**7)
!            vzhail = (zhdry/(zhdrycsmith*gmh**7))* 1.e18*dielf*rho_qhsq/rho_qrsq*cnoh*Sqrt(rho00/den)*gf7p5*((den*qs8)/ &
!                   (pie*cnoh*hwdn))**(0.25*(7.5))
            vzhail = (zhdry/(zhdrycsmith*gmh**7))*1.e18*dielf*rho_qhsq/rho_qrsq*cnoh*Sqrt(4.0*rho_qh*gr/(3*cd*den)) &
                      *gf7p5*((den*qh8)/(pie*cnoh*hwdn))**(0.25*(7.5))
          ELSE
            zhdry = zhdryc*gmh**7
            vzhail = 1.e18*dielf*rho_qhsq/rho_qrsq*cnoh*Sqrt(rho00/den)*gf7p5*((den*qh8)/(pie*cnoh*hwdn))**(0.25*(7.5))
          ENDIF
        ELSE
          IF ( iuseferrier .eq. 1 ) THEN
            zhwet = zhwetc*qh8**2/gmh
            vzhail = (zhwet/(zhdrycsmith*gmh**7)) *1.e18*dielf*rho_qhsq/rho_qrsq*cnoh*Sqrt(rho00/den)*gf7p5*((den*qh8)/ &
                      (pie*cnoh*hwdn))**(0.25*(7.5))
          ELSE
            zhwet = (zhwetc*(gmh**7))**.95
            vzhail = 1.e18*dielf*rho_qhsq/rho_qrsq*cnoh*Sqrt(rho00/den)*gf7p5*((den*qh8)/(pie*cnoh*hwdn))**(0.25*(7.5))
          ENDIF
        ENDIF
        
        ENDIF ! ( qh8 .gt. qhmin )
      
      
      ELSE ! ( nc value not recognized )

        write(*,*) 'qtodbz: *** ERROR *** nc =', nc
        stop
      
      ENDIF

 

!-----------------------------------------------------------------------
! Compute final reflectivity value
!-----------------------------------------------------------------------

      z0 = zrain + zswet + zsdry + zhwet + zhdry

      IF( z0 .gt. 0. ) THEN
        qtodbz = max(mindbz, 10.0 * log10(z0))
        IF ( qtodbz > mindbz ) THEN
          vzf = (vzsnow + vzrain + vzhail)/z0
          IF ( vzf > 90.0  ) THEN
!            write(*,*) 'qtovzf: vzf, s,r,h ',vzf,vzsnow, vzrain, vzhail
!            write(*,*) 'z: s,r,h,0: ', zswet + zsdry, zrain, zhwet, zhdry,z0
!            write(*,*) 'q: s,r,h: ',1.e3*qs8,1.e3*qr8,1.e3*qh8
!            write(*,*) 'dbz: ',qtodbz, tmp, gf7p5, gamma(7.5)
          ENDIF
          vzf = Min(100.0,vzf)
        ELSE
         vzf = 0.0
        ENDIF
      ELSE
        qtodbz = mindbz
        vzf = 0.0
      ENDIF

      RETURN
      END SUBROUTINE QTOVZF

!
! ##############################################################################
!


      SUBROUTINE setmicro(xcnoh,xrho_qh,xcnor,xrho_qr,xcnos,xrho_qs)      


      USE MICRO_MODULE
      USE INDEX_MODULE

      implicit none
      
      real xcnoh,xrho_qh,xcnor,xrho_qr,xcnos,xrho_qs
      
      cnoh = xcnoh
      cnor = xcnor
      cnos = xcnos
      
      rho_qh = xrho_qh
      rho_qr = xrho_qr
      rho_qs = xrho_qs
      
      
      RETURN
      END
      
!
! ##############################################################################
!
      SUBROUTINE setcnoz(cno)

      USE MICRO_MODULE
      USE INDEX_MODULE

      implicit none
      
      real cno(lc:lhab)


      cno(lc)  = 1.0e+08
      IF ( li .gt. 1 ) cno(li)  = 1.0e+08
!      cno(lr)  = 4.0e+05 
      cno(lr)  = cnor ! 8.0e+06 
      IF ( ls .gt. 1 ) cno(ls)  = cnos ! 8.0e+06 
!      cno(lgm) = 4.0e+04 
!      cno(lgm) = 2.0e+05 
      IF ( lh .gt. 1 ) cno(lh)  = cnoh ! 4.0e+05
      IF ( lhl .gt. 1 ) cno(lhl)  = cnohl ! 4.0e+05
      
      RETURN
      END

!
! ##############################################################################
!
      SUBROUTINE setcno(cno)
      
      USE MICRO_MODULE
      USE INDEX_MODULE

      implicit none
      logical, parameter :: iclz = .false. ! use CLZ intercepts for Binger
      
      real cno(lc:lhab)


      IF ( lc .ge. lc )  cno(lc)  = 1.0e+08
      IF ( li .ge. lc )  cno(li)  = 1.0e+08
      IF ( lir .ge. lc ) cno(lir) = 1.0e+08 
!      IF ( lr .ge. lc ) cno(lr)  = 4.0e+05 
      IF ( lr .ge. lc )  cno(lr)  = cnor_i10
      IF ( ls .ge. lc )  cno(ls)  = cnos_i10
!      IF ( lgm .ge. lc ) cno(lgm) = 4.0e+04 
      IF ( lgm .ge. lc ) cno(lgm) = cnogm_i10 
      IF ( lh .ge. lc )  cno(lh)  = cnoh_i10
!      IF ( lgl .ge. lc ) cno(lgl) = 4.0e+04 
!      IF ( lgh .ge. lc ) cno(lgh) = 4.0e+04 
      IF ( lgl .ge. lc ) cno(lgl) = cnogl_i10
      IF ( lgh .ge. lc ) cno(lgh) = cnogh_i10
      IF ( lf .ge. lc )  cno(lf)  = cnof_i10
      IF ( lhl .ge. lc ) cno(lhl) = cnohl_i10
      IF ( lip .ge. lc ) cno(lip) = 1.0e+08 
      
      IF ( iclz ) THEN
        
      IF ( lr .ge. lc )  cno(lr)  = 4.0e+05

      IF ( lgl .ge. lc ) cno(lgl) = 4.0e+05 
      IF ( lgm .ge. lc ) cno(lgm) = 9.0e+04
      IF ( lgh .ge. lc ) cno(lgh) = 4.0e+04 
      IF ( lf .ge. lc )  cno(lf)  = 4.0e+05

      IF ( lh .ge. lc )  cno(lh)  = 4.0e+04
      IF ( lhl .ge. lc ) cno(lhl) = 1.0e+03
        
      ENDIF
      
      RETURN
      END
!
! ##############################################################################
!
      SUBROUTINE setdn10(xdn0)
      
      USE MICRO_MODULE
      USE INDEX_MODULE

      implicit none
      
      real xdn0(lc:lhab)
      
      xdn0(lr)  = rho_qr_i10 
      xdn0(ls)  = rho_qs_i10 
      xdn0(lf)  = rho_qf_i10 
      xdn0(lgl) = rho_qgl_i10
      xdn0(lgm) = rho_qgm_i10
      xdn0(lgh) = rho_qgh_i10
      xdn0(lh)  = rho_qh_i10 
      xdn0(lhl) = rho_qhl_i10

      
      RETURN
      END
      


!
! ##############################################################################

!
!  subroutine to check for consistency of 3 moments after EnKF adjustments
!
      subroutine chkalpha(nx,ny,nz,nor,na, &
     &  an,dn,ipconc0,member )

      USE INDEX_MODULE
      USE COMMASMPI_MODULE
      USE MICRO_MODULE

!
      
      implicit none

      integer ng1
      parameter(ng1 = 1)
      
      integer nx,ny,nz,nor,ngt,jgs,na
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real dtp,dz,dtz1
      integer :: member
      
      integer ndebug0
      parameter (ndebug0 = -1)
      integer idx

      integer ix,jy,kz,i,j,k,item,il
!
!  include file for mix ratio and charge indices
!
!      include 'swm.index.zieg.h'
!
      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      
      real, allocatable :: qx(:,:) ! qx(ngs,lv:lhab)
      real, allocatable :: zx(:,:) ! zx(ngs,lr:lhab)
      real, allocatable :: cx(:,:) ! cx(ngs,lv:lhab)
      real, allocatable :: xv(:,:)
      real, allocatable :: xmas(:,:)
      real, allocatable :: xdn(:,:)
      real, allocatable :: xdia(:,:,:)
      real, allocatable :: vx(:,:)
      real, allocatable :: alpha(:,:)  ! alpha(ngs,lc:lhab)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)

!
! Fixed intercept values for single moment scheme
!
      
!      real cwccn,cwmasn,cwmasx,cimn,cimx,cwradn
      real rwmasn,rwmasx
      integer ngs1
      real cwc0

      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag,istag
      parameter (kstag=1, istag=1)
      

      integer ngs,ngscnt,mgs,ipconc0
      parameter ( ngs=500 )
      integer igs(ngs),kgs(ngs)
      
      real rho0(ngs),temcg(ngs)

      real temg(ngs)
      
      real rhovt(ngs)
      
      real cwnc(ngs),cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
      real cimasn,cimasx,cnina(ngs),cimas(ngs)

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
      integer ierr
      
      logical flag
      logical ldovol, ldoliq
      integer lvol(lc:lhab)
      integer ln(lc:lhab)
      integer lz(lc:lhab)
      integer lliq(li:lhab)
      
      real chw, qr, z, rd, alp, z1, g1, vr, nrx
      
      real zmin, ctmp


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.


#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: ENTERED SUBROUTINE, my_rank=",my_rank
#else
      if (ndebug0 .gt. 0 ) write(0,*) "ZIEGFALL: ENTERED SUBROUTINE"
#endif

! #####################################################################
! BEGIN EXECUTABLE
! #####################################################################
!

      call setqxmin(qxmin)

!  constants
!
      ldovol = .false.
      lvol(:) = 0
      IF ( lvi .gt. 1 ) lvol(li) = lvi
      IF ( lvs .gt. 1 ) lvol(ls) = lvs
      IF ( lvh .gt. 1 ) lvol(lh) = lvh
      IF ( lhl .gt. 1 .and. lvhl .gt. 1 ) lvol(lhl) = lvhl
      
      
      IF ( li .gt. 1 ) THEN
      DO il = li,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      ENDDO
      ENDIF

      ln(lc) = lnc
      ln(lr) = lnr
      IF ( li > 0 ) THEN
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      ENDIF
      IF ( lhl .gt. 1 ) ln(lhl) = lnhl

      lz(:) = 0
      lz(lr) = lzr
      IF ( li > 0 ) THEN
      lz(li) = lzi
      lz(ls) = lzs
      lz(lh) = lzh
      ENDIF
      IF ( lhl .gt. 1 .and. lzhl > 1 ) lz(lhl) = lzhl

      lliq(:) = 0
      IF ( lsw .gt. 1 ) lliq(ls) = lsw
      IF ( lhw .gt. 1 ) lliq(lh) = lhw
      IF ( lhl .gt. 1 .and. lhlw .gt. 1 ) lliq(lhl) = lhlw

      ldoliq = .false.
      IF ( ls .gt. 1 ) THEN
      DO il = ls,lhab
        ldoliq = ldoliq .or. ( lliq(il) .gt. 1 )
      ENDDO
      ENDIF

!
!  density maximums and minimums
!
      xdnmx(:) = 900.0
      
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(ls) =  300.0
      xdnmx(lh) =  900.0
      IF ( lhl .gt. 1 ) xdnmx(lhl) = 900.0
!
      xdnmn(:) = 900.0
      
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lh) =  170.0
      IF ( lhl .gt. 1 ) xdnmn(lhl) = 500.0

      xdn0(:) = 900.0
      
      xdn0(lc) = 1000.0
      xdn0(li) = 900.0
      xdn0(lr) = 1000.0
      xdn0(ls) = 100. ! rho_qs ! 100.0
      xdn0(lh) = 500. !rho_qh ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      IF ( lhl .gt. 1 ) xdn0(lhl) =  800. ! rho_qhl ! 800.0
      
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
      jye = ny
      IF ( jyend == nyend ) jye = ny-1
      DO jy = 1,jye
      jgs = jy
      
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
      allocate (  qx(ngscnt,lv:lhab) )
      allocate (  cx(ngscnt,lc:lhab) )
      allocate (  xdn(ngscnt,lc:lhab) )
      allocate (  vx(ngscnt,li:lhab) )
      allocate (  alpha(ngscnt,lr:lhab) )
      IF ( ipconc .ge. 6 ) THEN
      allocate (  zx(ngscnt,lr:lhab) )
      ENDIF

#ifdef MPI
      if (debug_mpi) write(0,*) my_rank, "ZIEGFALL: after allocate"
#endif

!
!  Reconstruct various quantities 
!
!

      do mgs = 1,ngscnt
        rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      ENDDO

      
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

!
!  set concentrations
!
      cx(:,:) = 0.0

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

       end do
      end if
       
      do mgs = 1,ngscnt
        xdn(mgs,lr) = xdn0(lr)
        IF ( lh .gt. 1 )  xdn(mgs,lh) = xdn0(lh)
        IF ( lhl .gt. 1 ) xdn(mgs,lhl) = xdn0(lhl)
      end do

!
! Set mean particle volume
!
      IF ( ldovol ) THEN
      
      vx(:,:) = 0.0
      
       DO il = lg,lhab
        
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

      dnu(lh) = alphah
      IF ( lhl .gt. 1 ) dnu(lhl) = alphahl

      DO il = lg,lhab
       DO mgs = 1,ngscnt
        alpha(mgs,il) = dnu(il)
      ENDDO
      ENDDO
      
      alpha(:,lr) = rnumin
!
! Set 6th moments
!
      IF ( ipconc .ge. 6 ) THEN
      
      zx(:,:) = 0.0
      
       DO il = lr,lhab
        
        IF ( lz(il) .ge. 1 ) THEN
        
          DO mgs = 1,ngscnt
            zx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0)
          ENDDO
          
        
        ENDIF
      
       ENDDO
      
      ENDIF
       


!  Find shape parameter rain

     IF ( lz(lr) > 1 ) THEN ! { RAIN SHAPE PARAM
          il = lr
          DO mgs = 1,ngscnt
         
         IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN

          IF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)

            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
      IF ( my_rank == 0 .and. .false.) THEN
      write(0,*) 'new cr is',cx(mgs,lr)*1000.
      ENDIF

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(1000*1000*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = cx(mgs,il)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
          ENDIF
          
          IF ( zx(mgs,lr) > 0.0 ) THEN
            vr = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-11,cx(mgs,lr)))
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
           DO i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,lr))/Abs(alpha(mgs,lr)) .lt. 1. ) EXIT
            IF ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) EXIT
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           print*,'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(mgs,lr),qr*1000,z*1.e18,vr,nrx

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
!           IF ( alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax ) THEN
           IF ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) THEN

            IF ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z

!            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)

!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z

              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)

            
            ELSEIF ( rescale_low_alpha .and. alp <= rnumin ) THEN

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

       
!  Find shape parameters for graupel,hail
      IF ( ipconc .ge. 6 ) THEN
        DO il = lg,lhab
        
        IF ( lz(il) .gt. 1 ) THEN
        
        DO mgs = 1,ngscnt

              IF ( cx(mgs,il) > 1.e10 ) THEN
                write(0,*) 'chkalpha0: large CHW: ',member,il,cx(mgs,il),qx(mgs,il),zx(mgs,il),xdn(mgs,il),igs(mgs),jgs,kgs(mgs)
              ENDIF

         IF ( qx(mgs,il) .gt. qxmin(il) ) THEN
          zmin = 1.e-19/0.224*(xdn0(lr)/xdn(mgs,il))**2
          IF ( zx(mgs,il) <=  0.0 .and. cx(mgs,il) <= 0.0 ) THEN
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = 0.0
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
          ELSEIF ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
!            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) > 0.0 ) THEN
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2)
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
           ELSEIF ( zx(mgs,il) <= 0.0 .and. cx(mgs,il) <= 0.0 ) THEN
!   How did this happen?
!              write(0,*) 'ziegfall: something screwy with moments: il = ',il
!              write(0,*) 'q,n,z = ', 1.e3*qx(mgs,il),cx(mgs,il),zx(mgs,il)
!              write(0,*) 'alpha = ',alpha(mgs,il)
         ! set values according to dBZ of -10
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn(mgs,il))**2
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
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = 0.0
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
          ENDIF
         
         ELSE
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = 0.0
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
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
           DO i = 1,10
            IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
!     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
           
!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
           IF ( .true. .or. alpha(mgs,il) <= alphamin+0.1 .or. alp == alphamin .or. alp >= alphamax-0.1 ) THEN

             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

            IF ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) THEN  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
              
              IF ( cx(mgs,il) > 1.e10 ) THEN
                write(0,*) 'chkalpha: large CHW: ',member,il,cx(mgs,il),chw,alpha(mgs,il),g1,qr,z,xdn(mgs,il),igs(mgs),jgs,kgs(mgs)
              ENDIF
            
            ELSEIF ( rescale_low_alpha .and. alp <= alphamin ) THEN

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

              IF ( cx(mgs,il) > 1.e10 ) THEN
                write(0,*) 'chkalpha1: large CHW: ',member,il,cx(mgs,il),qx(mgs,il),zx(mgs,il),xdn(mgs,il),igs(mgs),jgs,kgs(mgs)
              ENDIF
        ENDDO ! mgs
        
        ENDIF ! lz(il) .gt. 1
        
        ENDDO ! il
          
      ENDIF

      IF ( ipconc .ge. 7 ) THEN
        IF ( lhl .gt. 1 ) THEN
        
        ENDIF
      ENDIF

!
! put fall speeds into the x-z arrays
!


      deallocate ( qx )
      deallocate ( cx )
      deallocate ( xdn )
      deallocate ( vx )
      deallocate ( alpha )
      IF ( allocated(zx) ) deallocate ( zx )

 9998 continue

      if (ndebug0 .gt. 0 ) write(0,*)  my_rank,'ZIEGFALL: DONE WITH LOOP'

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

      if (ndebug0 .gt. 0 ) print*,'ZIEGFALL: SET NZMPB'

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

      if (ndebug0 .gt. 0 ) print*,'ZIEGFALL: SET NXMPB'

 1200 continue


!       ENDDO ! ix
!      ENDDO ! kz
      ENDDO ! jy

#ifdef MPI
      if (debug_mpi) write(0,*) "ZIEGFALL: EXITING SUBROUTINE, my_rank=",my_rank,' jgs = ',jgs
#else
      if (ndebug0 .gt. 0 ) write(0,*) "ZIEGFALL: EXITING SUBROUTINE"
#endif


      RETURN
      END

