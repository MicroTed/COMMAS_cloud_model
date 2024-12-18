c------------------------------------------------------------------------------
c
c   /////////////////////          BEGIN         \\\\\\\\\\\\\\\\\\\\\\
c   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE WARMLFO   //////////////////////
c
c WARMLFO is the warm rain portion of the LFO microphysics scheme from
c         Gilmore et al. 2003.
c
c Created by:  Louis Wicker, January, 2006
c Latest Revision: 01-20-06
c
c Converted to ijk indexing by ERM, March 2006
c
c------------------------------------------------------------------------------
      SUBROUTINE WARMLFO(t,             ! potential temperature
     $                   qv,            ! water vapor mixing ratio
     $                   qc,            ! cloud water mixing ratio at t 
     $                   qr,            ! rain  water mixing ratio at t 
     $                   precip,        ! 2D precip array 
     $                   tb, pb,        ! base state theta and pi
     $                   dzc,           ! vertical map factor
     $                   dt,            ! time step and grid size
     $                   nx, ny, nz)    ! grid dimensions
      
      USE MICRO_MODULE
      USE PARAM_MODULE
#ifdef MPI
      USE COMMASMPI_MODULE
#endif
      USE TRAJ_MODULE

      implicit none

c------------------------------------------------------------------------------
c Passed variables

      integer nx, ny, nz
      real     t    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qv    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qc    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    qr    (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
      real    precip(-ng+1:nx+ng,-ng+1:ny+ng,4)
      real    dzc  (-ng+1:nz+ng)
      real    tb   (-ng+1:nz+ng)
      real    pb   (-ng+1:nz+ng)
      real    dt
      
c------------------------------------------------------------------------------
c LFO PARAMETERS

      real dnz00      ; parameter ( dnz00   = 1.225      )
      real br         ; parameter ( br      = 0.8        )
      real ar         ; parameter ( ar      = 841.9966   )
      real cs         ; parameter ( cs      = 4.83607122 )
      real ds         ; parameter ( ds      = 0.25       )
      real dragh      ; parameter ( dragh   = 0.6        )
      real elv        ; parameter ( elv     = 2500300.   )
      real tfr        ; parameter ( tfr     = 273.16     )
      real advisc0    ; parameter ( advisc0 = 1.832e-05  )
      real advisc1    ; parameter ( advisc1 = 1.718e-05  )
      real tka0       ; parameter ( tka0    = 2.43e-02   )
      real cwdn       ; parameter ( cwdn    = 1000.0     )
      real cwmasn     ; parameter ( cwmasn  = 4.25e-15   )
      real cwmasx     ; parameter ( cwmasx  = 5.25e-10   )
      real pid4       ; parameter ( pid4    = pii / 4.0  )
      real c1f3       ; parameter ( c1f3    = 1.0/3.0    )
      real cwc1       ; parameter ( cwc1    = 6.0/(pii*cwdn) )
      real erw        ; parameter ( erw     = 1.0        )
      real qrmin      ; parameter ( qrmin   = 1.0e-10    )
      real qcmin      ; parameter ( qcmin   = 1.0e-10    )

      real temp(nz)
      real xvtr(nz) 
      real xvts(nz) 
      real xvth2(nz)
      real xslop(nz)
      real wvdf(nz) 
      real advisc(nz)
      real akvisc(nz)
      real tka(nz)  
      real schm(nz)  
      real rwd    
      real ssw   
      real xbv   
      real xav    
      real xxcev  
      real xce    
      real xvevr2 
      real xrv    
      real crw
      real xvevr1
      real gf4br 
      real gf4ds 
      real qr5dt
      real qc5dt
      real gf1
      real gf2
      real gf3
      real cwd
      real rwd2
      real cwd2
      real vtrbar
      real vtwbar
      real cwmas
      real qccrit
      real qdiff
      real argrcnw1
      real argrcnw
      real dqr
      real dqc
      real dqv
      

      real gamma
      external gamma

c------------------------------------------------------------------------------
c Local variables

      integer i, j, k, nfall, n
      real product, factor, f5, qvs
      real prod0
      real ern
      real cr
      real rhoa(nz), rcgs(nz), rcgsi(nz)
      real pinit(nz)
      real gam(nz)
      real vtden(nz)
      real tempg, dtfall
      real esat
      real satrat
      real e1,e2,e3,e4
      real gammae,psi

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      logical :: debug_mpi = .TRUE.
#endif

c These used to be passed arrays for diagnostics, here, simply allocate the memory locally

      real fq    (0:nx+1,0:nz+1)
      real vt    (0:nx+1,0:nz+1)
      real qracw (nx,nz)
      real qrcnw (nx,nz)
      real qvevr (nx,nz)

      parameter( f5 = 237.3 * 17.27 * 2.5e6 / cp )

c------------------------------------------------------------------------------
c Compute some constants

      xvevr1 = 0.78*gamma(2.0)
      gf1    = gamma(1.0)
      gf2    = gamma(2.0)
      gf3    = gamma(3.0)
      gf4ds  = gamma(4.0+ds)
      gf4br  = gamma(4.0+br)

c------------------------------------------------------------------------------
c Create constants

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      do k = kzb,kze 
#else
      DO k = 1,nz-1
#endif
       rhoa(k)   = 1.0e5*pb(k)**2.509/(287.04*tb(k))
       pinit(k)  = 1.0e5*pb(k)**3.509
       gam(k)    = Lav/(Cp*pb(k))
       temp(k)   = tb(k)*pb(k)

       rcgs(k)   = 1.0e2*pb(k)**2.509/(287.04*tb(k))
       rcgsi(k)  = 1.0 / rcgs(k)
       vtden(k)  = sqrt(0.0011225/rcgs(k))

       xvtr(k)   = (ar*gf4br*(dnz00/rhoa(k))**(0.5))/(6.0)
       xslop(k)  = (rhoa(k)/(pii))**(0.25)
       wvdf(k)   = (2.11e-05)*((temp(k)/tfr)**1.94)* (101325.0/(pinit(k)))
       advisc(k) = advisc0*(416.16/(temp(k)+120.0))*(temp(k)/296.0)**(1.5)
       akvisc(k) = advisc(k)/rhoa(k)
       tka(k)    = tka0*advisc(k)/advisc1
       schm(k)   = akvisc(k) / wvdf(k)

      ENDDO

c------------------------------------------------------------------------------
c Outer JI LOOP
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do j = jyb,jye
       fq    (:,:) = 0.0
       vt    (:,:) = 0.0

       do k = kzb,kze 
        do i = ixb,ixe
#else
      DO j = 1,ny-1

c Set some stuff to zero
       fq    (:,:) = 0.0
       vt    (:,:) = 0.0

c Compute autoconversion, coalesense using LFO equations

       DO k = 1,nz-1
        DO i = 1,nx-1
#endif
         qracw(i,k) = 0.0
         qrcnw(i,k) = 0.0

c Coalescence

         vtrbar = 0.0
         rwd    = 0.0

         IF ( qr(i,j,k) .gt. qrmin )  THEN
          rwd    = xslop(k)* (qr(i,j,k)/(rho_qr*cnor))**(0.25)
          vtrbar = xvtr(k)*(rwd**br)
         ENDIF

         vtwbar = 0.0
         cwd    = 0.0

         IF ( qc(i,j,k) .gt. qcmin )  THEN
          cwmas  = min( max(qc(i,j,k)*rhoa(k)/ccn,cwmasn),cwmasx )
          cwd    = (cwmas*cwc1)**c1f3                           ! min diam gives 2e-6, max diam gives 100e-6
          vtwbar = (ar*(cwd**br))*(dnz00/rhoa(k))**0.5          ! MSG Extrapolation of Liu and Orville to qcw droplet sizes (Gilmore et al. 2004)
c         vtwbar = g*cwdn*(cwd**2)/( advisc(k) * 18. )          ! MSG Alternative equation - Stokes Law
         ENDIF     

         crw     = cnor*rwd
         rwd2    = rwd**2
         cwd2    = cwd**2

         qc5dt = 0.20 * qc(i,j,k) / dt

         qracw(i,k) = min(erw*qc(i,j,k)*pid4*crw*abs(vtrbar-vtwbar)*(gf3*rwd2+2.*gf2*rwd*cwd+gf1*cwd2),qc5dt)

c AutoConversion qcw -> qrw
c Autoconversion == 0 --> do standard autoconversion theshold
c Autoconversion =/ 1 --> Berry 1968 auto conversion for rain (Ferrier 1994) with time delay


!         IF ( qc(i,j,k) .gt. 1.e-3 ) THEN
!         print*, 'autoconversion = ',autoconversion,qcmincwrn,qc(i,j,k)
!         ENDIF
         IF( autoconversion .eq. 0 ) THEN
          qdiff  = max( qc(i,j,k) - qcmincwrn, 0.)
         ELSE
          qccrit = (pii/6.)*((ccn*cwdiap**3)*cwdn)/rhoa(k)
          qdiff  = max( qc(i,j,k) - qccrit,0.)
         ENDIF

         IF ( qdiff .gt. 0.0 ) THEN
!          print*, 'Make rain'
          argrcnw1   = (1.596e-12)*(1e-6) / (cwdisp*1e-3)
          argrcnw    = 1.2e-4 + (argrcnw1*ccn) / (rhoa(k)*qdiff)
          qrcnw(i,k) = min(max(rhoa(k)*1e-3*(qdiff**2)/argrcnw,0.0),qc5dt)
         ENDIF

c BEGIN LFO EVAP / COND

         tempg      = pb(k) * t(i,j,k)

         qr5dt      = 0.20*qr(i,j,k)/dt
         ssw        = qv(i,j,k)/qvs
         xbv        = (1.0 / (rhoa(k)*wvdf(k)))
         xav        = (elv**2) / (tka(k)*rw)
         xxcev      = 2.0 * pii / rhoa(k)
         xce        = xxcev * (ssw - 1.0) * ( 1.0 / (xav / tempg**2 + xbv / qvs) )
         xvevr2     = 0.308 * (gamma((5.0+br)/2.0)) * (schm(k)**(1./3.)) * (akvisc(k)**(-0.5))*(ar**(0.5))*((dnz00/rhoa(k))**(0.25))
         xrv        = (xvevr1*rwd + xvevr2*(rwd**((3.0+br)/2.0)))
         qvevr(i,k) = -1.0 * max(min(xce * crw * xrv, 0.0), -qr5dt)

c END   LFO EVAP / COND

c Adjust mixing ratios and account for evap cooling...

         dqr        = dt*(qrcnw(i,k) + qracw(i,k) - qvevr(i,k))
         dqc        = dt*( - qrcnw(i,k) - qracw(i,k) )
         dqv        = dt * qvevr(i,k)
 

         qv(i,j,k)  = max(qv(i,j,k) + dqv,0.0)
         qc(i,j,k)  = max(qc(i,j,k) + dqc,0.0)
         qr(i,j,k)  = max(qr(i,j,k) + dqr,0.0)
         t (i,j,k)  = t(i,j,k) - gam(k)*dqv

c Latent heat of vaporization/condensation

         tempg      = pb(k) * t(i,j,k)
         qvs        = 380.*exp(17.27*(tempg - 273.)/(tempg - 36.))/pinit(k)
         product    = (qv(i,j,k) - qvs) / (1. + qvs*f5/(tempg - 36.)**2)

         product    = max(product,-qc(i,j,k))
         t (i,j,k)  = t(i,j,k) + gam(k)*product
         qv(i,j,k)  = max(qv(i,j,k) - product,0.)
         qc(i,j,k)  = max(qc(i,j,k) + product,0.)

         IF ( allocated( microrates ) ) THEN
!         microrates(i,j,k,1) = 0.0
!         microrates(i,j,k,2) = 0.0
          microrates(i,j,k,3) = gam(k)*(product)*3600. ! gam(k) = Lav/(Cp*pb(k))
          microrates(i,j,k,4) = -qvevr(i,k)*gam(k)*3600.
!         microrates(i,j,k,5) = -qrmlh(kz)*elf/(cp*piz(kz))*3600.
!         microrates(i,j,k,6) = pcond2(kz)
         ENDIF

         IF (qc(i,j,k).lt.qcmin) THEN
           qv(i,j,k) = qv(i,j,k) + qc(i,j,k)
           qc(i,j,k)=0.0
         ENDIF
         IF (qr(i,j,k).lt.qrmin) THEN
           qv(i,j,k) = qv(i,j,k) + qc(i,j,k)
           qr(i,j,k)=0.0
         ENDIF

        ENDDO !! i
       ENDDO !! k
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze 
        do i = ixb,ixe
#else
       DO k = 1,nz-1
        DO i = 1,nx-1
#endif
         IF ( qr(i,j,k) .gt. qrmin ) THEN
          rwd       = xslop(k)* (qr(i,j,k)/(rho_qr*cnor))**(0.25)
          vt(i,k) = xvtr(k)*(rwd**br)
          fq(i,k) = qr(i,j,k)
!         IF ( qr(i,j,k) .gt. 3.e-3 ) THEN
!          print*,'vt: qr,vt,rwd: ',
!     $         qr(i,j,k),vt(i,k),rwd,xvtr,rho_qr,cnor
!         ENDIF
         
         ENDIF

        ENDDO !! i
       ENDDO !! k

c------------------------------------------------------------------------------
c Compute instantaneous precip rate on grid
c    den of air    mix rat.   fall velo x  time step  x   1/den of water
c     g /cm^3   x    g/g    x    cm/s   x                   cm^3 / g    = precip rate (mm/s)
c     g /cm^3   x    g/g    x    cm/s   x     DT      x     cm^3 / g    = depth of precip (mm)
#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do i = ixb,ixe
#else
       DO i = 1,nx-1
#endif
        precip(i,j,1) = rcgs(1)*qr(i,j,1)*1000.*vt(i,1)
        precip(i,j,3) = precip(i,j,3) + precip(i,j,1)*dt
       ENDDO
c------------------------------------------------------------------------------
c Fallout
#ifdef MPI
        kzb = 1
        kze = ktile
        if (kzend .eq. nzend) kze = kzend-kzbeg

        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do k = kzb,kze
         factor = dt*dzc(k)*rcgsi(k)
         do i = ixb,ixe
#else
        DO k = 1,nz-1
         factor = dt*dzc(k)*rcgsi(k)
         DO i = 1,nx-1
#endif
!         IF ( qr(i,j,k) .gt. 3.e-3 ) THEN
!          print*,'qfall: qr,vt,factor: ',qr(i,j,k),vt(i,k),factor
!         ENDIF
          qr(i,j,k) = qr(i,j,k) - factor
     $                * (rcgs(k  )*fq(i,k  )*vt(i,k  )
     $                  -rcgs(k+1)*fq(i,k+1)*vt(i,k+1))

         ENDDO  !! i
        ENDDO  !! k

      ENDDO  !! j-loop

      RETURN  
      END SUBROUTINE WARMLFO
