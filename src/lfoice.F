c-------------------------------------------------------------------------------
c 
c   /////////////////////         BEGIN            \\\\\\\\\\\\\\\\\\\\
c   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE LFO_DRIVER  ////////////////////
c
c LFO_ICE_DRIVE interfaces COMMAS model with Straka 3-class LFO ice scheme
c-------------------------------------------------------------------------------
      SUBROUTINE LFO_ICE_DRIVE(pb, tb, qb,                   ! pibar, tbar, qvbar
     $                         pi,                           ! pi'
     $                         s,                            ! scalars
     $                         precip,                       ! precip
     $                         dt, dzc,                      ! constants
     $                         nx, ny, nz, ns, ng )
     
     
      USE INDEX_MODULE, only: lt,lc,lr,li,ls,lh,lhl,lf,lv
      USE MICRO_MODULE
      USE COMMASMPI_MODULE

c-------------------------------------------------------------------------------
c     PURPOSE
c
c     calculate and apply the microphysical contributions to the water,
c     ice and temperature fields, using an ice microphysics parameterization
c     scheme.
c-------------------------------------------------------------------------------

c INPUT

c       nx       number of grid points in the x-direction (east/west)
c       ny       number of grid points in the y-direction (north/south)
c       nz       number of grid points in the vertical

c       dt       the large time step size for this call.

c       pb       base state Exner function
c       tb       base state potential temperature (K)
c       qb       base state mixing ratio (kg/kg)
c       p        perturbation Exner function
c       t        potential temperature  (K)
c       qv       water vapor specific humidity at all time levels (kg/kg)
c       qc       cloud water mixing ratio at all time levels (kg/kg)
c       qr       rainwater mixing ratio at all time levels (kg/kg)
c       qi       cloud ice mixing ratio at all time levels (kg/kg)
c       qs       snow mixing ratio at all time levels (kg/kg)
c       qg       hail mixing ratio at all time levels (kg/kg)
c       precip   accumulated grid-scale rainfall (mm)

c-------------------------------------------------------------------------------

c     variable declarations.

c-------------------------------------------------------------------------------

      implicit none

      integer nsint, nsreal
      integer nx,ny,nz,ng,ns       ! number of grid points in 3 directions
      real dt                ! the large time step size for this call.
      real dzc(-ng+1:nz+ng)
      
      real precip(-ng+1:nx+ng,-ng+1:ny+ng,4)   ! precipitation

      real tb (-ng+1:nz+ng)           ! base state potential temperature (k)
      real qb (-ng+1:nz+ng)           ! base state water vapor specific humidity (kg/kg)
      real pb (-ng+1:nz+ng)           ! base state pressure

      real  pi(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)     ! perturbation Exner
      real  s (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns)  !  (t,qv,qc,qr,qi,qs,qh)

c-------------------------------------------------------------------------------

c     temporary arrays

c-------------------------------------------------------------------------------

      real sb(-ng+1:nz+ng,ns)                       ! base state array
      real rinit(-ng+1:nz+ng)                       ! base state density
      real pinit(-ng+1:nz+ng)                       ! base state pi
  
c-------------------------------------------------------------------------------

c     misc. local variables

c-------------------------------------------------------------------------------

      integer i,j,k
      

      integer constset
      save constset
      data constset /0/

      real*8 cqtotn, cqvtotn, cqctotn, cqrtotn, cqitotn, cqstotn, cqhtotn
      real*8 cqtotp, cqvtotp, cqctotp, cqrtotp, cqitotp, cqstotp, cqhtotp
      real*8 cqfac,  cqvfac,  cqcfac,  cqrfac,  cqifac,  cqsfac,  cqhfac

c Species threshold,  mass,    mass constraints

      real qcmin,   cwmas, cwmasn, cwmasx
      real qimin,   cimas, cimasn, cimasx
      real qrmin
      real qsmin
      real qhmin
      parameter (qcmin   = 1.0e-09, cwmasn = 4.25e-15, cwmasx = 5.25e-10)   
      parameter (qimin   = 1.0e-12, cimasn = 4.25e-15, cimasx = 3.23e-8)      ! cimasx corresponds to Di=3000 microns.
      parameter (qrmin   = 1.0e-07)
      parameter (qsmin   = 1.0e-07)
      parameter (qhmin   = 1.0e-07)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      logical :: debug_mpi = .FALSE.

#ifdef MPI
      if(debug_mpi) write(0,*) "LFO_ICE_DRIVE: ENTERED SUBROUTINE, my_rank=",my_rank
#else
!      write(0,*) "LFO_ICE_DRIVE: ENTERED SUBROUTINE"
#endif

c-------------------------------------------------------------------------------

c     IF (constset .eq. 0) THEN
c
c       write(6,*)
c       write(6,*) '2004 Gilmore et al. microphysics'
c       write(6,*)
c       IF( hole_fill .eq. 0 ) write(6,*) 'HOLE_FILLING OFF '
c       IF( hole_fill .eq. 1 ) write(6,*) 'GLOBAL HOLE FILLING TURN ON'
c       IF( hole_fill .eq. 2 ) write(6,*) 'SPECIES HOLE FILLING TURN ON'
c       write(6,*)
c
c       constset = 1
c
c     ENDIF

      sb(:,:)  = 0.0
      rinit(:) = 0.0
      pinit(:) = 0.0


c Set up 1D vertical arrays

#ifdef MPI
      kzb = 0
      kze = ktile+1
      if( kzbeg .eq. nzbeg ) kzb = 1
      if( kzend .eq. nzend ) kze = kzend - kzbeg
      
      do k = kzb,kze
#else
      DO k = 1,nz-1
#endif
       sb(k,1) = tb(k)
       sb(k,2) = qb(k)
       rinit(k)= 1.0e5*pb(k)**2.509/(287.04*tb(k))
       pinit(k)= 1.0e5*pb(k)**3.509
      ENDDO


c DO MICROPHYSICS

       CALL LFO_ICE (nx,ny,nz,ns,ng,dt,dzc,sb,s,precip,rinit,pinit,pi,pb)

c-------------------------------------------------------------------------------

c HOLE FILLING SCHEME:  WATER MASS IS CONSERVED

c-------------------------------------------------------------------------------
c HOLE_FILL == 1 --> GLOBAL METHOD ACROSS ALL SPECIES
#ifdef MPI
       hole_fill = 0 ! only do hole-filling in solver under MPI
#endif

      IF( hole_fill .eq. 1 ) THEN

       cqtotn = 0.0
       cqtotp = 0.0

#ifdef MPI
       kzb = 1
       kze = ktile
       if(kzend .eq. nzend) kze = kzend - kzbeg
       
       jyb = 1
       jye = jtile
       if(jyend .eq. nyend) jye = jyend - jybeg
       
       ixb = 1
       ixe = itile
       if(ixend .eq. nxend) ixe = ixend - ixbeg
       
       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =  1,nz-1
        DO j =  1,ny-1
         DO i =  1,nx-1
#endif
          cqtotn = cqtotn + rinit(k)*(s(i,j,k,2)+s(i,j,k,3)+s(i,j,k,4)+s(i,j,k,5)+s(i,j,k,6)+s(i,j,k,7))

          s(i,j,k,2) = max(s(i,j,k,2),0.0)
          s(i,j,k,3) = max(s(i,j,k,3),0.0)
          s(i,j,k,4) = max(s(i,j,k,4),0.0)
          s(i,j,k,5) = max(s(i,j,k,5),0.0)
          s(i,j,k,6) = max(s(i,j,k,6),0.0)
          s(i,j,k,7) = max(s(i,j,k,7),0.0)

          cqtotp = cqtotp + rinit(k)*(s(i,j,k,2)+s(i,j,k,3)+s(i,j,k,4)+s(i,j,k,5)+s(i,j,k,6)+s(i,j,k,7))

         ENDDO
        ENDDO
       ENDDO

c CREATE ratio

       cqfac = (cqtotn+1.0e-20)/(cqtotp+1.0e-20)

c ADJUST fields

#ifdef MPI
       kzb = 1
       kze = ktile
       if(kzend .eq. nzend) kze = kzend - kzbeg
       
       jyb = 1
       jye = jtile
       if(jyend .eq. nyend) jye = jyend - jybeg
       
       ixb = 1
       ixe = itile
       if(ixend .eq. nxend) ixe = ixend - ixbeg
       
       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =  1,nz-1
        DO j =  1,ny-1
         DO i =  1,nx-1
#endif
          s(i,j,k,2) = s(i,j,k,2) * cqfac
          s(i,j,k,3) = s(i,j,k,3) * cqfac
          s(i,j,k,4) = s(i,j,k,4) * cqfac
          s(i,j,k,5) = s(i,j,k,5) * cqfac
          s(i,j,k,6) = s(i,j,k,6) * cqfac
          s(i,j,k,7) = s(i,j,k,7) * cqfac

         ENDDO
        ENDDO
       ENDDO

      ENDIF

c-------------------------------------------------------------------------------
c HOLE_FILL == 2 --> LOCAL FILL METHOD FROM EACH SPECIES

      IF( hole_fill .eq. 2 ) THEN

       cqvtotn = 0.0
       cqvtotp = 0.0
       cqctotn = 0.0
       cqctotp = 0.0
       cqrtotn = 0.0
       cqrtotp = 0.0
       cqitotn = 0.0
       cqitotp = 0.0
       cqstotn = 0.0
       cqstotp = 0.0
       cqhtotn = 0.0
       cqhtotp = 0.0

#ifdef MPI
       kzb = 1
       kze = ktile
       if(kzend .eq. nzend) kze = kzend - kzbeg
       
       jyb = 1
       jye = jtile
       if(jyend .eq. nyend) jye = jyend - jybeg
       
       ixb = 1
       ixe = itile
       if(ixend .eq. nxend) ixe = ixend - ixbeg
       
       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =  1,nz-1
        DO j =  1,ny-1
         DO i =  1,nx-1
#endif
          cqvtotn = cqvtotn + rinit(k)*s(i,j,k,2)
          cqctotn = cqctotn + rinit(k)*s(i,j,k,3)
          cqrtotn = cqrtotn + rinit(k)*s(i,j,k,4)
          cqitotn = cqitotn + rinit(k)*s(i,j,k,5)
          cqstotn = cqstotn + rinit(k)*s(i,j,k,6)
          cqhtotn = cqhtotn + rinit(k)*s(i,j,k,7)

          s(i,j,k,2) = max(s(i,j,k,2),0.0)
          s(i,j,k,3) = max(s(i,j,k,3),0.0)
          s(i,j,k,4) = max(s(i,j,k,4),0.0)
          s(i,j,k,5) = max(s(i,j,k,5),0.0)
          s(i,j,k,6) = max(s(i,j,k,6),0.0)
          s(i,j,k,7) = max(s(i,j,k,7),0.0)

          cqvtotp = cqvtotp + rinit(k)*s(i,j,k,2)
          cqctotp = cqctotp + rinit(k)*s(i,j,k,3)
          cqrtotp = cqrtotp + rinit(k)*s(i,j,k,4)
          cqitotp = cqitotp + rinit(k)*s(i,j,k,5)
          cqstotp = cqstotp + rinit(k)*s(i,j,k,6)
          cqhtotp = cqhtotp + rinit(k)*s(i,j,k,7)

         ENDDO
        ENDDO
       ENDDO

c CREATE ratio

       cqvfac = (cqvtotn+1.0e-20)/(cqvtotp+1.0e-20)
       cqcfac = (cqctotn+1.0e-20)/(cqctotp+1.0e-20)
       cqrfac = (cqrtotn+1.0e-20)/(cqrtotp+1.0e-20)
       cqifac = (cqitotn+1.0e-20)/(cqitotp+1.0e-20)
       cqsfac = (cqstotn+1.0e-20)/(cqstotp+1.0e-20)
       cqhfac = (cqhtotn+1.0e-20)/(cqhtotp+1.0e-20)

c ADJUST fields

#ifdef MPI
       kzb = 1
       kze = ktile
       if(kzend .eq. nzend) kze = kzend - kzbeg
       
       jyb = 1
       jye = jtile
       if(jyend .eq. nyend) jye = jyend - jybeg
       
       ixb = 1
       ixe = itile
       if(ixend .eq. nxend) ixe = ixend - ixbeg
       
       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k =  1,nz-1
        DO j =  1,ny-1
         DO i =  1,nx-1
#endif
          s(i,j,k,2) = s(i,j,k,2) * cqvfac
          s(i,j,k,3) = s(i,j,k,3) * cqcfac
          s(i,j,k,4) = s(i,j,k,4) * cqrfac
          s(i,j,k,5) = s(i,j,k,5) * cqifac
          s(i,j,k,6) = s(i,j,k,6) * cqsfac
          s(i,j,k,7) = s(i,j,k,7) * cqhfac

         ENDDO
        ENDDO
       ENDDO

      ENDIF

c      DO k = 1,nz-1
c       DO j = 1,ny-1
c        DO kz = 1,nz-1
c 
c          IF ( an(ix,jy,kz,lr) .le. qrmin ) THEN
c           an(ix,jy,kz,lc) = an(ix,jy,kz,lc) + an(ix,jy,kz,lr)
c           an(ix,jy,kz,lr) = 0.0
c          ENDIF
c
c          IF ( an(ix,jy,kz,lc) .le. qcmin ) THEN
c           an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
c           an(ix,jy,kz,lc) = 0.0
c          ENDIF
c
c          IF ( an(ix,jy,kz,lh) .le. qhmin ) THEN
c           an(ix,jy,kz,ls) = an(ix,jy,kz,ls)+ an(ix,jy,kz,lh)
c           an(ix,jy,kz,lh) = 0.0
c          ENDIF
c 
c          IF ( an(ix,jy,kz,ls) .le. qsmin ) THEN
c           an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,ls)
c           an(ix,jy,kz,ls) = 0.0
c          ENDIF
c
c          IF ( an(ix,jy,kz,li) .le. qimin ) THEN
c           an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
c           an(ix,jy,kz,li) = 0.0
c          ENDIF
c
c        ENDDO
c       ENDDO
c      ENDDO

#ifdef MPI
      if(debug_mpi) write(0,*) "LFO_ICE_DRIVE: EXITING SUBROUTINE, my_rank=",my_rank
#else
!      write(0,*) "LFO_ICE_DRIVE: EXITING SUBROUTINE"
#endif

      RETURN
      END
c--------------------------------------------------------------------------
c
c
c  LFO 3-ICE MODEL 
c
c  VERSION:  RELEASE 1.1.1 (internal version #20 from Gilmore)
c
c  
c  LIN-FARELY-ORVILLE-like "Simple Ice and Liquid Microphysics Scheme"
c   based on Gilmore et al. (2004) Mon. Wea. Rev.
c
c   Copyright (C) <2004>  <Jerry Straka and Matthew Gilmore>
c
c   This library is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Lesser General Public
c   License as published by the Free Software Foundation; either
c   version 2.1 of the License, or (at your option) any later version.
c
c   This library is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c   Lesser General Public License for more details.
c
c   If you find this code useful and publish results using it, please reference:
c
c     Gilmore M. S., J. M. Straka, and E. N. Rasmussen,
c           Monthly Weather Review: Vol. 132, No. 8, pp. 1897-1916.
c--------------------------------------------------------------------------

       SUBROUTINE LFO_ICE(nx, ny, nz, na, nor, dtp, dzc, ab, an, precip,
     $                    db, pb, pin, pib)

c--------------------------------------------------------------------------
c
c RELEASE NOTES
c
c  MSG - 01/03/03  Deleted all unnecessary variables.
c        Also grouped varables by process.  Made all variables explicit.
c
c  MSG - 11/26/03  sped up code by omitting micro and fallout for cells with no activity.
c        This saved between 25 and 30% of total model wallclock time on several test runs.
c        Saves approx 50% of the micro CPU that was being used in previous versions.
c        Very small differences noted after 3 hour simulation.
c
c  ERM - 11/08/03  Reordered arrays to be j,k,i to avoid cost of copying arrays
c
c  MSG - 03/21/04  Corrected a few bugs which were making -qinitauto flag necessary in version 1.0. 
c        That flag is now no longer needed during compilation for version 1.1.
c
c  LJW - 03/22/04  Found a few more uninitialized variables and corrected them.  Also
c        commented out all IF tests relating to debugging and other diagnostics.
c        Also moved downscaling to end of routine and now checks ALL grid zones.
c
c  LJW - 06/18/04  Added 3rd order PD Crowley fluxes for fallout, added precip calculation.
c        !!!NOT USED!!!!
c
c  LJW - 08/01/04  Added Berry autoconversion formula option based on conversation with Ted Mansell
c
c  LJW - 11/30/04  Introduce "ice_param" array to pass in parameters set in initialization about 
c                  number concentrations and other microphysics flags that get changed...
c
c  5/30/2006 - ERM - Switch to sat. adj. as in Kessler for T > 0C.
c------------------------------------------------------------------------------------------------------

c  Declarations

c     USE IEEE_ARITHMETIC   ! USED WITH THE IBM COMPILER & see code below

      USE PARAM_MODULE
      USE MICRO_MODULE
      USE COMMASMPI_MODULE
      USE TRAJ_MODULE

      implicit none

!	include 'const.h'

c-------------------------------------------------------------------------------------------------------
c PARAMETERS FOR MICROPHYSICS CONSTANTS ARE NOW IN MICRO_MODULE
c-------------------------------------------------------------------------------------------------------
 
c Passed Variables

      integer nx, ny, nz, na, nor
      real dtp                                                ! time step in seconds
      real dzc(-nor+1:nz+nor)                                            ! 1/dz (inverse of physical vertical grid spacing)
      real ab(-nor+1:nz+nor,na), pb(-nor+1:nz+nor)
      real db(-nor+1:nz+nor), pib(-nor+1:nz+nor), pnz(-nor+1:nz+nor)        ! base state fields ab is base state of scalar fields
      real pin(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)                                      ! perturbation Exner function

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)                                    ! scalar array containing theta and hydrometeor mixing ratios

      real precip(-nor+1:nx+nor,-nor+1:ny+nor,4)            ! precip array
	                                                        ! precip(:,:,1) == instantaneous rain rate
	                                                        ! precip(:,:,2) == instantaneous hail rate
	                                                        ! precip(:,:,3) == accumulated (on moving grid) rainfall
	                                                        ! precip(:,:,4) == accumulated (on moving grid) hailfall

c Specification of location of scalar variables (so they can be stored in any order, or the scalar array can contain other
c                                                fields such as TKE, etc, in which case na > 7)

	integer lt ; parameter( lt = 1 )				  ! memory location of theta field: an(:,:,:,lt)
	integer lv ; parameter( lv = 2 )				  ! memory location of vapor field: an(:,:,:,lv)
	integer lc ; parameter( lc = 3 )				  ! memory location of cloud field: an(:,:,:,lc)
	integer lr ; parameter( lr = 4 )				  ! memory location of rain  field: an(:,:,:,lr)
	integer li ; parameter( li = 5 )				  ! memory location of ice   field: an(:,:,:,li)
	integer ls ; parameter( ls = 6 )				  ! memory location of snow  field: an(:,:,:,ls)
	integer lh ; parameter( lh = 7 )				  ! memory location of hail  field: an(:,:,:,lh)

c Local variables

      integer ix, jy, kz
      integer kzfstart, kzfend    
      integer nzm; parameter (nzm = 201)                      ! USED TO CREATE SCRATCH SPACE, NZ < NZM
	
c Misc constants

      integer gcnt(nzm), gcntcolumn
      logical gsat(nzm),gqcw(nzm),gqci(nzm),gqrw(nzm),gqsw(nzm),gqhw(nzm)
      real advisc0,advisc1
      real ar,br,bta1,cs,ds
      real cpi,cnit,cwc1,c1f3
      real dragh
      real elv,elf,els   
      real gf6,gf5,gf4,gf3,gf2,gf1,gf4br,gf4ds,gf4p5
      real pi,pid4
      real tfr,tfrcbw,tfrcbi,thnuc,tka0
      parameter (advisc0 = 1.832e-05, advisc1 = 1.718e-05)
      parameter (cpi = 1.0/cp)
!      parameter (cbi  = 7.66   , cbw  =   35.86)
      parameter (tka0 = 2.43e-02)
      parameter (tfr  = 273.16,  thnuc= 233.16, tfrcbw= tfr - cbw, tfrcbi= tfr - cbi)

c Misc Variables

      real wvdf(nzm),akvisc(nzm),ci,cw,tka(nzm),advisc(nzm),schm(nzm)
      real cc3

c Species mixing ratio, # conc., fall velocity, density, intercept params
      real qwv(nzm), qv0
      real qcw(nzm),  ccw,  vtwbar,  cwdn
      real qci(nzm),  cci,  vtibar
      real qrw(nzm),  crw,  vtrbar,  rwdn(nzm)
      real qsw(nzm),  csw,  vtsbar,  swdn(nzm)
      real qhw(nzm),  chw,  vthbar,  hwdn(nzm)

c Species threshold,  mass,    mass constraints

      real qcmin,   cwmas, cwmasn, cwmasx
      real qimin,   cimas, cimasn, cimasx
      real qrmin
      real qsmin
      real qhmin
      parameter (qcmin   = 1.0e-09, cwmasn = 4.25e-15, cwmasx = 5.25e-10)   
      parameter (qimin   = 1.0e-12, cimasn = 4.25e-15, cimasx = 3.23e-8)      ! cimasx corresponds to Di=3000 microns.
      parameter (qrmin   = 1.0e-07)
      parameter (qsmin   = 1.0e-07)
      parameter (qhmin   = 1.0e-07)
	real qccrit

c Species inverse slope, inverse slope^2, part of slope 

      real cwd,  cwd2 
      real cid,  cid2
      real rwd,  rwd2,  xslop(nzm)
      real swd,  swd2
      real hwd,  hwd2

c Fallout, fall velocity parameters/vars

      real dtz1               
      real vtrdzmax           ! input courant                         ! MSG
      real dtsplit            ! small timestep for rain sedimentation ! MSG
      integer nrnstp,inrnstp  ! number of small timesteps, counter    ! MSG
      real cflux(nzm), pflux(nzm), rflux(nzm), sflux(nzm), hflux(nzm)
      real xvtr(nzm),xvts(nzm),xvth1,xvth2(nzm)
	real eps ; parameter( eps = 1.0e-15 )
      real crh, crs, cri, crr, crc
      real qi, qoutmin, qoutmax, qinmax, qinmin, qD, qC, qU

C Air and particle temperatures (equivalent to each other in LFO version)

      real temp(nzm),theta(nzm),thetap(nzm),theta0(nzm)                        
      real temg(nzm),temgsqr,temcg(nzm)

C Air pressure, density, & Exner function

      real pbz(nzm), pres(nzm),presp(nzm),   poo
      real rho0(nzm), dnz(nzm), dnz00, rho00
      real piz(nzm)
      real gam(nzm)
      parameter ( dnz00 = 1.225, rho00 = 1.225 )
      parameter ( poo   = 1.e+05 )

c Saturation Adjustment

      real qvap(nzm)
      real dqvcnd(nzm),dqwv(nzm),dqcw(nzm),dqci(nzm)
      real cqv1,cqv2
      real gamss,denom1,denom2
      real qwvp(nzm)
      real fraci(nzm),fracl(nzm)
      real cdw, cdi, product, f5
      parameter( f5 = 237.3 * 17.27 * 2.5e6 / cp )
      integer itertd
      
c Saturation lookup table, vapor pressures, ratio's     

      integer ltemq, l
!      parameter(nqsat=20001)
       real temq
!      parameter (fqsat=0.01)
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati,temq
!      parameter (fqsat=0.002,fqsati=1./fqsat)

!      real tabqvs(nqsat),tabqis(nqsat)
!      real  cai, caw
!      parameter (cai  = 21.87, caw = 17.27)   ! rounded to hundredth's place
      real qvs(nzm),qis(nzm),qss(nzm),pqs(nzm)
      real ssi,ssw
      integer :: imaktab = 0
      save imaktab

c Accretions

      real qiacr(nzm)
      real qracw(nzm),qraci(nzm),qracs(nzm)           
      real qsacw(nzm),qsaci(nzm),           qsacr(nzm)
      real qhacw(nzm),qhaci(nzm),qhacs(nzm),qhacr(nzm)
      real xacx
      
      real   erw,  eri,  ers       ! Efficiencies 
      real   esw,  esi
      real   ehw,  ehi,  ehs,  ehr
 
c Biggs Freezing

      real qhfzr(nzm), xrfrz, arz, brz
      parameter (arz = 0.66, brz = 100.0)

c Bergeron Process

      real qsfw(nzm),qsfi(nzm), eic, bsfw
      real cs10(32),cs11(32),cs9
      real bfa1(32),bfa2(32),cbtim(32)
      real cmn, cmi40, cmi50, ri50, vti50, a, cm50a, cm40b, cm50b
      integer ib
      parameter (cmn = 1.05e-15, cmi40 = 2.4546e-07, cmi50 = 4.8e-07)
      parameter (ri50 = 5.0e-03, vti50 = 100.0)

c    Note that the following table of a1 values from Koenig (1971) has an error at T=-14 C.
c    A reviewer noted that 0.1725e-6 should really be 0.1725e-4 (fifteenth index of bfa1).  
c    This was corrected below.

      data bfa1/0., 0.7939E-7, 0.7841E-6, 0.3369E-5, 0.4336E-5,
     $   0.5285E-5, 0.3728E-5, 0.1852E-5, 0.2991E-6, 0.4248E-6,
     $   0.7434E-6, 0.1812E-5, 0.4394E-5, 0.9145E-5, 0.1725E-4,
     $   0.3348E-4, 0.1725E-4, 0.9175E-5, 0.4412E-5, 0.2252E-5,
     $   0.9115E-6, 0.4876E-6, 0.3473E-6, 0.4758E-6, 0.6306E-6,
     $   0.8573E-6, 0.7868E-6, 0.7192E-6, 0.6513E-6, 0.5956E-6,
     $   0.5333E-6, 0.4834E-6/

      data bfa2/0., 0.4006, 0.4831, 0.5320, 0.5307, 0.5319,
     $      0.5249, 0.4888, 0.3894, 0.4047, 0.4318, 0.4771,
     $      0.5183, 0.5463, 0.5651, 0.5813, 0.5655, 0.5478,
     $      0.5203, 0.4906, 0.4447, 0.4126, 0.3960, 0.4149,
     $      0.4320, 0.4506, 0.4483, 0.4460, 0.4433, 0.4413,
     $      0.4382, 0.4361/
      
c Conversions

      real qhcns(nzm),qrcnw(nzm),qscni(nzm), qdiff, argrcnw, argrcnw1
      real ehscnv, esicnv     !Efficiencies

c Deposition/Sublimation

      real xas,xbs, xds, xxdsv
      real qhdsv(nzm),qhdpv(nzm),qvsbh(nzm)
      real qsdsv(nzm),qsdpv(nzm),qvsbs(nzm)
      
c Evaporation

      real qvevr(nzm), xav,xbv
      real xvevr1,xvevr2,xce,xrv,xxcev
 
c Initiation of cloud ice

      real dqisdt,qiint(nzm),xiint,cnnt

c Melting/Wet growth of hail, Melting snow

      real qrmls(nzm),qrmlh(nzm)
      real xvth3      
      real xmlt1,xmlt2,xmlt3
      real xhmlt1,xhmlt2
      real xsmlt1,xsmlt2
      real xsv, xhv, xhsw

c Wet/dry growth, shedding

      real qhdry(nzm),qhwet(nzm),qrshh(nzm)
      real xhwet1,xhwet2,xcwt,xwt1
      real qhacip,qhacsp           !Accretions if wet growth

c Water Budgets

c     real ptotal(nzm)
      real pqcwi(nzm),pqcii(nzm),pqrwi(nzm)
      real pqswi(nzm),pqhwi(nzm),pqwvi(nzm)
      real pqcwd(nzm),pqcid(nzm),pqrwd(nzm)
      real pqswd(nzm),pqhwd(nzm),pqwvd(nzm)
      integer il2(nzm),  il3(nzm),  il5(nzm)    !boolean operators

c Latent Heating Computation

      real psub(nzm),pvap(nzm),pfrz(nzm),ptem(nzm),pcond2(nzm)

c Maximum Depletion Tendencies

      real qc5dt,qi5dt,qr5dt,qs5dt,qh5dt 

c Functions

      real gamma
      external gamma

c Other declarations

c      integer autoconversion 
c      real    rho_qr, cnor 
c      real    rho_qs, cnos
c      real    rho_qh, cnoh 
c      real    qcmincwrn        
c      real    cwdiap          
c      real    cwdisp         
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      
      logical :: debug_mpi = .FALSE.

c Flags   

      integer imake !, ndebug     !prints intercept/density to output file
      data imake / 0 /
      
c      ndebug = 0        !prints debugging flags
#ifdef MPI
      if ( debug_mpi )  write(0,*) "LFO_ICE: ENTERED SUBROUTINE, my_rank=",my_rank
#else
!      write(0,*) "LFO_ICE: ENTERED SUBROUTINE"
#endif

c      IF ( ndebug .gt. 0 ) print*,'dbg = 0a'


c  Other Constants

      ar    = 841.99666        ! In mks units....(identical to LFO83 of 2115 in cgs)
      br    = 0.8
      bta1  = 0.6
      cs    = 4.83607122
      ds    = 0.25  
      cnit  = 1.0e-02
      c1f3  = 1.0/3.0
      pi    = 4.0*atan(1.0)
      pid4  = pi/4.0
      gf1   = gamma(1.0)
      gf2   = gamma(2.0)
      gf3   = gamma(3.0)
      gf4   = gamma(4.0)
      gf5   = gamma(5.0)
      gf6   = gamma(6.0)
      gf4br = gamma(4.0+br)
      gf4ds = gamma(4.0+ds)
      gf4p5 = gamma(4.0+0.5)

      cwdn  = 1000.0
      cwc1  = 6.0/(pi*cwdn)
      elv   = 2500300.             ! From LFO83 (uses a constant to avoid the problem of a singularity)
      elf   = 335717.              ! From LFO83 (ditto)
      els   = elv + elf
      xacx  = pid4/gf4
      xvth1 = gf4p5/6.
      xmlt3 = -4218./elf
      xhmlt1= 0.78*gamma(2.0)
      xsmlt1= 0.78*gamma(2.0)
      xvevr1= 0.78*gamma(2.0)
      cqv1  = 4098.0258*elv*cpi
      cqv2  = 5807.6953*els*cpi

c Set up passed in constants

      ccw            = Abs(ccn)
      
      IF (rho_qh .lt. 600.) THEN        ! MSG (5/21/02) Pruppacher and Klett (1978) 
        dragh = 1.0                     ! Hail experiences more drag than graupel.
      ELSE                              ! This was used in Gilmore et al. (2004).
        dragh = 0.6
      ENDIF

      IF ( nz .gt. nzm-1 )  THEN
       write(6,*) 'STRAKA/GILMORE LFO MICROPHYSICS:  PROBLEM!!!'
       write(6,*) 'STOP, NZM= ',nzm,' must be => NZ= ',nz
       stop
      END IF

c     IF ( imake .eq. 0 ) THEN
c
c      write(6,*) '-----------------------------------------------------------------------'
c      write(6,*) 'LFO MICROPHYSICAL CONSTANTS'
c
c101    format(1x,a,6(g10.5,2x))
c
c      write(6,101) 'CNOR/DENR:        ',cnor, rho_qr
c
c      write(6,101) 'CNOS/DENS:        ',cnos, rho_qs
c
c      write(6,101) 'CNOH/DENH/DRAGH:  ',cnoh, rho_qh,dragh
c
c      write(6,*)   'TIME STEP:        ', dtp
c
c      IF( autoconversion .eq. 0 ) THEN
c       write(6,*) 'Simple qc -> qr autoconversion, QC threshold = ',qcmincwrn
c      ELSE
c       write(6,*) 'Berry (1968) AC / CCW = ',ccw
c      ENDIF
c 
c      write(6,*) '-----------------------------------------------------------------------'
c
c      IF ( ndebug .gt. 0 ) print*,'dbg = 0b'
c      imake = 1
c
c     ENDIF

c Generate lookup table for vapor pressures
      IF ( imaktab .eq. 0 ) THEN
      imaktab = 1
      DO l = 1,nqsat
       temq = 163.15 + (l-1)*fqsat
       tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
       tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
c      IF (IEEE_IS_NAN(tabqvs(l))) stop 'L'
c      IF (IEEE_IS_NAN(tabqis(l))) stop 'M'
      ENDDO
      ENDIF
c  Define lookup table for Bergeron process, note in cgs units

       bsfw     = (ri50**2)*vti50*pi
       cbtim(1) = 1.0e+30
       cs10(1)  = 0.0                   ! LJW bug fix 03/22/04
       cs11(1)  = 0.0                   ! LJW bug fix 03/22/04
       DO ib = 2,32
        cm50a     = cmi50**bfa2(ib)
        a         = 1.0-bfa2(ib)
        cm40b     = cmi40**a
        cm50b     = cmi50**a
        cbtim(ib) = (cm50b-cm40b)/(bfa1(ib)*a)
        cs10(ib)  = bfa1(ib)*cm50a
        cs11(ib)  = dtp/(cbtim(ib)*cmi50)
       ENDDO

c  Define altitude-dependent variables

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend - kzbeg
!      IF ( kzbeg == nzbeg ) kzb = 1
      do kz = kzb,kze
#else
      DO kz = 1,nz-1
#endif

      IF (ndebug .gt. 0 ) print*,'dbg = 0c'

        piz(kz)    = pib(kz) ! (pb(kz)/poo)**rcp ! = pib(kz)
        dnz(kz)    = db(kz)
        pbz(kz)    = pb(kz)
        temp(kz)   = piz(kz)*ab(kz,lt)
        gam(kz)    = (elv) / (pib(kz)*cp)
c       if (IEEE_IS_NAN(piz(kz))) stop 'A'
c       if (IEEE_IS_NAN(dnz(kz))) stop 'B'
c       if (IEEE_IS_NAN(pbz(kz))) stop 'C'
c       if (IEEE_IS_NAN(temp(kz))) stop 'D'

        xvtr(kz)   = (ar*gf4br*(dnz00/dnz(kz))**(0.5))/(6.0)
        xvts(kz)   = (cs*gf4ds*(dnz00/dnz(kz))**(0.5))/(6.0)
        xvth2(kz)  = (4.0*g/(3.0*dragh*dnz(kz)))**(0.5)
        xslop(kz)  = (dnz(kz)/(pi))**(0.25)
        wvdf(kz)   = (2.11e-05)*((temp(kz)/tfr)**1.94)* (101325.0/(pbz(kz)))
        advisc(kz) = advisc0*(416.16/(temp(kz)+120.0))*(temp(kz)/296.0)**(1.5)
        akvisc(kz) = advisc(kz)/dnz(kz)
        tka(kz)    = tka0*advisc(kz)/advisc1
        schm(kz)   = (akvisc(kz)/wvdf(kz))
      ENDDO

c      IF (ndebug .gt. 0 ) print*,'dbg = 1'

c Start ix loop

#ifdef MPI
      jyb = 1
      jye = jtile
      if( jyend .eq. nyend ) jye = jyend - jybeg

      ixb = 1
      ixe = itile
      if( ixend .eq. nxend ) ixe = ixend - ixbeg

      do jy = jyb,jye ; do ix = ixb,ixe
#else
      DO jy = 1,ny-1                            !  doing XZ slabs here
      DO ix = 1,nx-1
#endif

c      IF (ndebug .gt. 0 ) print*,'start of xz slab'

c      Reset all tendency arrays for the column

       il2(:) =0
       il3(:) =0
       il5(:) =0
       gsat(:)=.false.
       gqcw(:)=.false.
       gqci(:)=.false.
       gqrw(:)=.false.
       gqsw(:)=.false.
       gqhw(:)=.false.
       gcnt(:)  = 0   
       gcntcolumn = 0
#ifdef MPI
       kzfstart  = -1
       kzfend    = nzend-1
#else
       kzfstart  = -1
       kzfend    = nz-1
#endif
       vtrdzmax = 0.0
       qracw(:) = 0.0
       qsacw(:) = 0.0
       qhacw(:) = 0.0
       qraci(:) = 0.0
       qsaci(:) = 0.0
       qhaci(:) = 0.0
       qiacr(:) = 0.0
       qhacs(:) = 0.0
       qhacr(:) = 0.0
       qracs(:) = 0.0
       qsacr(:) = 0.0
       qsfw(:)  = 0.0                                                                                
       qsfi(:)  = 0.0
       qscni(:) = 0.0
       qhcns(:) = 0.0
       qhfzr(:) = 0.0
       qrcnw(:) = 0.0
       qrmls(:) = 0.0
       qrmlh(:) = 0.0
       qsdsv(:) = 0.0
       qhdsv(:) = 0.0
       qvsbs(:) = 0.0
       qvsbh(:) = 0.0
       qsdpv(:) = 0.0
       qhdpv(:) = 0.0
       qvevr(:) = 0.0
       qiint(:) = 0.0
       qrshh(:) = 0.0
       pfrz(:) = 0.0
       psub(:) = 0.0
       pvap(:) = 0.0

c      IF (ndebug .gt. 0 ) print*,'dbg = 4'

c Define temporaries to be used in calculations...done for all k even if no microphysics

#ifdef MPI
      kzb = 1
      kze = ktile+1
      if(kzend .eq. nzend) kze = kzend - kzbeg
      do kz = kzb,kze
#else
      DO kz = 1,nz-1                                       ! Main k loop for microphysics
#endif

c         IF (ndebug .gt. 0 ) print*,'start of mp k-loop'

         pnz(kz)   = 1.0e5*((pin(ix,jy,kz)+pib(kz))**3.509) - pb(kz)
         pqs(kz)   = 380.0/(pnz(kz)+pb(kz))
         qwv(kz)   = max(an(ix,jy,kz,lv), 0.0)
         qcw(kz)   = max(an(ix,jy,kz,lc), 0.0)
         qci(kz)   = max(an(ix,jy,kz,li), 0.0)
         qrw(kz)   = max(an(ix,jy,kz,lr), 0.0)
         qsw(kz)   = max(an(ix,jy,kz,ls), 0.0)
         qhw(kz)   = max(an(ix,jy,kz,lh), 0.0)
         theta0(kz)= an(ix,jy,kz,lt) ! ab(kz,lt)
         thetap(kz)= 0.0 ! an(ix,jy,kz,lt) - ab(kz,lt)
         theta(kz) = an(ix,jy,kz,lt)
         presp(kz) = pnz(kz)
         pres(kz)  = presp(kz) + pb(kz)
!         temg(kz)  = theta(kz)* (pres(kz) / poo ) ** rcp
         
         temg(kz) = theta(kz)*(pib(kz) + pin(ix,jy,kz))
         
         temcg(kz) = temg(kz) - tfr
         rho0(kz)  = db(kz)
         qwvp(kz)  = 0.0 ! an(ix,jy,kz,lv) - ab(kz,lv)
         rwdn(kz)  = rho_qr     ! was 1000.0....now defined above
         swdn(kz)  = rho_qs     ! was 100.0
         hwdn(kz)  = rho_qh     ! was 900.0

c        IF (IEEE_IS_NAN(temg(kz))) stop '8a9'

         ltemq     = nint((temg(kz)-163.15)/fqsat+1.5)
         ltemq     = min(max(ltemq,1),nqsat)
         qvs(kz)   = pqs(kz)*tabqvs(ltemq)
         qis(kz)   = pqs(kz)*tabqis(ltemq)

!      IF ( jy .eq. 1 .and. ix .eq. 24 ) THEN
!       write(91,*) 'kz,qv,th: ',kz,an(ix,jy,kz,lv),an(ix,jy,kz,lt),pqs(kz),tabqvs(ltemq),qvs(kz)
!      ENDIF

         IF( temg(kz) .lt. tfr ) THEN
          IF( qcw(kz) .ge. 0.0 .and. qci(kz) .eq. 0.0) qss(kz) = qvs(kz)
          IF( qcw(kz) .eq. 0.0 .and. qci(kz) .gt. 0.0) qss(kz) = qis(kz)
          IF( qcw(kz) .gt. 0.0 .and. qci(kz) .gt. 0.0) qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz))/(qcw(kz) + qci(kz))
         ELSE
          qss(kz) = qvs(kz)
!           IF ( an(ix,jy,kz,lv)  .gt. qss(kz) ) write(91,*) 'qss exceeded at ',ix,jy,kz,qss(kz),an(ix,jy,kz,lv),temg(kz)
         ENDIF

         IF (ndebug .gt. 0 ) print*,'dbg = 5a'

c  Check whether to do microphysics for the gridcell or not

         IF ( qwv(kz) .gt. qss(kz)) gsat(kz) = .true.
         IF ( qcw(kz) .gt. qcmin  ) gqcw(kz) = .true.
         IF ( qci(kz) .gt. qimin  ) gqci(kz) = .true.
         IF ( qrw(kz) .gt. qrmin  ) gqrw(kz) = .true.
         IF ( qsw(kz) .gt. qsmin  ) gqsw(kz) = .true.
         IF ( qhw(kz) .gt. qhmin  ) gqhw(kz) = .true.
         IF (gsat(kz) .or. gqcw(kz) .or. gqci(kz) .or. gqrw(kz) .or. gqsw(kz) .or. gqhw(kz) ) gcnt(kz) = 1

         IF (ndebug .gt. 0 ) print*,'dbg = 5b'

         IF (gcnt(kz) .gt. 0) THEN      ! Do microphysics for this grid cell!

         IF ( ndebug .gt. 0)  THEN
c          print *, ix,jy,kz,'Micro. gcnt',gcnt(kz),'gsat:',gsat(kz),'gqcw:',gqcw(kz),'gqci:',gqci(kz), 
c   $              'gqrw:',gqrw(kz),'gqsw:',gqsw(kz),'gqhw:',gqhw(kz)
          print *, 'ICE ON ', ix, jy, kz, gcnt(kz)
         ENDIF

c Set Boolean operators

          IF ( temg(kz) .lt. tfr ) THEN
           il5(kz) = 1
           if ( qrw(kz) .lt. 1.0e-04 .and. qsw(kz) .lt. 1.0e-04 ) il2(kz) = 1
           if ( qrw(kz) .lt. 1.0e-04 ) il3(kz) = 1
          ENDIF

         IF (ndebug .gt. 0 ) print*,'dbg = 6'

c Set misc variables

          cwd    = 0.0
          cid    = 0.0
          rwd    = 0.0
          swd    = 0.0
          hwd    = 0.0
          vtibar = 0.0
          vtwbar = 0.0
          vtrbar = 0.0
          vtsbar = 0.0
          vthbar = 0.0
           cci   = 0.0

c COMPUTE AVERAGE DIAMETERS, AVERAGE TERMINAL VELOCITIES, AND NUMBER CONCENTRATIONS
c Cloud water variables

          IF ( gqcw(kz) ) THEN
            cwmas= min( max(qcw(kz)*rho0(kz)/ccw,cwmasn),cwmasx )
            cwd  = (cwmas*cwc1)**c1f3                           ! min diam gives 2e-6, max diam gives 100e-6
            vtwbar = (ar*(cwd**br))*(rho00/rho0(kz))**0.5       ! MSG Extrapolation of Liu and Orville to qcw droplet sizes (Gilmore et al. 2004)
c           vtwbar = g*cwdn*(cwd**2)/( advisc(kz) * 18. )      ! MSG Alternative equation - Stokes Law
          ENDIF 

c Cloud ice variables

          IF ( gqci(kz) ) THEN
            cci = max(min(cnit*exp(-temcg(kz)*bta1),1.e+09),1.0)
            cimas = min( max(qci(kz)*rho0(kz)/cci,cimasn),cimasx )
            IF ( temcg(kz) .le. 0 ) cid = min(max(16.7*(cimas**(0.5)),1.e-5),1.e-2) ! MSG thresholding not really needed since applied cimasn already
            vtibar = (cs*(cid**ds))*(rho00/rho0(kz))**0.5      ! MSG Snow fall velocity equation applied for cloudice.  (Gilmore et al. 2004)
c           vtibar = (58.51*(cid**.824))*(rho00/rho0(kz))**0.5 ! MSG Alternative "slower" formula from Pruppacher and Klett (1997) in mks units.
          ENDIF

c Marshall-Palmer distribution information for rain, snow aggs, and graupel/hail
c Definitions for marshall palmer distribution variables (rain, snow, hail) when mixing ratio only is predicted
c MSG - JMS used to use a threshold of 1.0e-10 for q mins.  This was inconsistent with the q mins set elsewhere.
c LJW - model blows up if dont put thresholding mixing ratios going into the diameters here. 1/2/03

        IF (ndebug .gt. 0 ) print*,'dbg = 7a'

         erw = 1.0
         IF ( gqrw(kz) ) THEN 
           rwd = xslop(kz)* (qrw(kz)/(rwdn(kz)*cnor))**(0.25)
           vtrbar = xvtr(kz)*(rwd**br)
           if ( rwd .le. 200.0e-6 ) erw = 0.0
           IF ( ierw .eq. 0 ) erw = erw0
         ENDIF
         IF ( gqsw(kz) ) THEN
           swd = xslop(kz)* (qsw(kz)/(swdn(kz)*cnos))**(0.25)
           vtsbar = xvts(kz)*(swd**ds)
         ENDIF
         IF ( gqhw(kz) ) THEN 
           hwd = xslop(kz)* (qhw(kz)/(hwdn(kz)*cnoh))**(0.25)
           vthbar = xvth1*xvth2(kz)*((hwdn(kz)*hwd)**(0.5))
         ENDIF
         
         rwd2 = rwd**2
         swd2 = swd**2
         hwd2 = hwd**2
         cwd2 = cwd**2
         cid2 = cid**2

         crw  = cnor*rwd
         csw  = cnos*swd
         chw  = cnoh*hwd

c Maximum depletion tendency by any one source

         qc5dt = 0.20*qcw(kz)/dtp
         qi5dt = 0.20*qci(kz)/dtp
         qr5dt = 0.20*qrw(kz)/dtp
         qs5dt = 0.20*qsw(kz)/dtp
         qh5dt = 0.20*qhw(kz)/dtp

c Collection/Conversion efficiencies 

         eic = 1.0
         eri = 1.0
         esi = esilfo0*exp(0.025*min(temcg(kz),0.0))        
         esicnv = esi                                               ! ERM includes esicnv instead of esi for snow AC
         esw = 1.0
         ers = 1.0
         ehw = 1.0
         ehr = 1.0
         ehs = 1.0
         IF ( temcg(kz) .le. 0.0 ) ehs = ehslfo0*exp(0.09*min(temcg(kz),0.0))
         ehscnv = ehs                                               ! LJW & ERM includes ehscnv instead of ehs for hail AC
         ehi = 0.1

c Constants for hydrometeor-vapor interactions

        ssi = qwv(kz)/qis(kz)
        temgsqr = temg(kz)**2

c Variables for Melting & deposition/sublimation of snow and hail

        IF ( gqsw(kz) .or. gqhw(kz) ) THEN      
          xbs       = (1.0/(dnz(kz)*wvdf(kz)))
          xas       = (els**2)/(tka(kz)*rw)
          xxdsv     = 2.0*pi/dnz(kz)
          xmlt2     = wvdf(kz)*elv*dnz(kz)
          xds       = xxdsv*(ssi-1.0) * (1.0/(xas/temgsqr+xbs/qis(kz)))
          xhsw      = (tka(kz)*temcg(kz) + xmlt2*(qwv(kz)-pqs(kz)))
          xmlt1     = -2.0*pi/(elf*dnz(kz))
        ENDIF

        IF (  gqhw(kz)) THEN     !HAIL ACTIVE

c Accretions

          IF (gqcw(kz)) qhacw(kz)=min(ehw*qcw(kz)*pid4*chw*abs(vthbar-vtwbar)*(gf3*hwd2+2.*gf2*hwd*cwd+gf1*cwd2),qc5dt)
          IF (gqci(kz)) qhaci(kz)=min(ehi*qci(kz)*pid4*chw*abs(vthbar-vtibar)*(gf3*hwd2+2.*gf2*hwd*cid+gf1*cid2),qi5dt)
          IF (gqsw(kz) .and. temcg(kz) .le. 0.0 ) 
     &        qhacs(kz)=min(ehs*qsw(kz)*xacx*chw*abs(vthbar-vtsbar)*(gf6*gf1*swd2+2.*gf5*gf2*swd*hwd+gf4*gf3*hwd2),qs5dt)
          IF (gqrw(kz)) qhacr(kz)=min(ehr*qrw(kz)*xacx*chw*abs(vthbar-vtrbar)*(gf6*gf1*rwd2+2.*gf5*gf2*rwd*hwd+gf4*gf3*hwd2),qr5dt)

c Misc variables needed

          xvth3   = (4.0*g/(3.0*dragh*dnz(kz)))**(0.25)
          xhmlt2  = 0.308*gamma(2.75)*(schm(kz)**(1./3.))*(akvisc(kz)**(-0.5))
          xhv     = (xhmlt1*(hwd) + (hwdn(kz)**(0.25))*xvth3* xhmlt2*(hwd**(1.75)))

c Melting, sublimation, deposition of hail

          qrmlh(kz) = min( (xmlt1*chw*xhv*xhsw + temcg(kz)*xmlt3*(qhacr(kz)+qhacw(kz)) ), 0.0 )
          qrmlh(kz) = -1.0*max( qrmlh(kz), -qh5dt )
          qhdsv(kz) = xds*chw*xhv  ! LJW & ERM fix - explicit thresholding for dep/sub of hail (MSG AGREES)
          qvsbh(kz) = abs(max( min(qhdsv(kz), 0.0), -qh5dt ))
          qhdpv(kz) =     max(qhdsv(kz), 0.0)

c=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=
c Wet growth/Shedding evaluation

          qhacsp = min((qhacs(kz)/ehs),qs5dt)  !MSG - ehs can be small....threshold by qs5dt for consistency  3/21/04
          qhacip = min((qhaci(kz)/ehi),qi5dt)  !MSG - ehi is never zero...still threshold by qi5dt for consistency 3/21/04

          IF ( temg(kz) .gt. tfr ) THEN            !Melting - Special shedding
          
            qhdry(kz) = 0.0
            qhwet(kz) = 0.0
            qrshh(kz) = qhacr(kz) +qhacw(kz) +qhacsp +qhacip      ! MSG - changed to qhacsp and qhacip for consistency 3/21/04

          ELSEIF ( temg(kz) .lt. 243.15 ) THEN     !Dry - No shedding

            qhdry(kz) = qhacr(kz) + qhacw(kz) + qhaci(kz) + qhacs(kz)
            qhwet(kz) = 0.0
            qrshh(kz) = 0.0

          ELSE                                     !Possible wet growth shedding
          
            qhdry(kz) = qhacr(kz) +qhacw(kz) + qhaci(kz) +qhacs(kz) 
            xcwt      = 1.0 / ( elf + 4218.*temcg(kz) )
            xhwet1    = 2.0*pi/dnz(kz)
            xhwet2    = dnz(kz)*elv*wvdf(kz)
            xwt1      = xhwet2*(pqs(kz)-qwv(kz)) -tka(kz)*temcg(kz)
            ci        = (2.118636 + 0.007371*(temp(kz)-tfr) )*(1.0e+03)
            qhwet(kz) = max(0.0,  (xhv*chw*xwt1*xhwet1*xcwt) + (1.0-ci*temcg(kz)*xcwt) * (qhacip + qhacsp)  )

            IF ( (qhwet(kz) .lt. qhdry(kz)) .and. (qhwet(kz) .gt. 0.0) ) THEN
              qhdry(kz) = 0.0                                 ! Wet growth case
              qrshh(kz) = -qhwet(kz) +(qhacw(kz) +qhacr(kz))  ! MSG -  qhacip,qhacsp are implicitly included in qhwet (see above)
            ELSE                                              ! Dry growth (defaults here if qhwet<0 due to a singularity at low T)
              qhwet(kz) = 0.0
              qrshh(kz) = 0.0
            END IF

           END IF    !end of shedding evaluation

           IF ( qrshh(kz) .gt. 0.0 ) THEN   ! Reset some variables if wet qh surface
            qhaci(kz) = qhacip
            qhacs(kz) = qhacsp
            qhdpv(kz) = 0.0
            qvsbh(kz) = 0.0
           END IF

        ENDIF                    !END OF ACTIVE HAIL GRIDCELLS

c=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=-=-=-=-0-=

        IF (  gqsw(kz)) THEN     !SNOW ACTIVE

c Accretions

         IF (gqcw(kz)) qsacw(kz)=min(esw*qcw(kz)*pid4*csw*abs(vtsbar-vtwbar)*(gf3*swd2+2.*gf2*swd*cwd+gf1*cwd2),qc5dt)
         IF (gqci(kz)) qsaci(kz)=min(esi*qci(kz)*pid4*csw*abs(vtsbar-vtibar)*(gf3*swd2+2.*gf2*swd*cid+gf1*cid2),qi5dt)
         IF (gqrw(kz)) qsacr(kz)=min(ers*qrw(kz)*xacx*csw*abs(vtsbar-vtrbar)*(gf6*gf1*rwd2+2.*gf5*gf2*rwd*swd+gf4*gf3*swd2),qr5dt)


c Melting, sublimation, deposition of snow

         xsmlt2    = 0.308*(gamma((5.0+ds)/2.0))*(schm(kz)**(1./3.))*(akvisc(kz)**(-0.5))*(cs**(0.5))*((dnz00/dnz(kz))**(0.25))
         xsv       = (xsmlt1*(swd) + xsmlt2*(swd**((3.0+ds)/2.0)))
         qrmls(kz) = min( (xmlt1*csw*xsv*xhsw + temcg(kz)*xmlt3*(qsacr(kz)+qsacw(kz)) ), 0.0 )
!        qrmls(kz) = -1.0*max( qrmls(kz), -qs5dt )
         qrmls(kz) = -1.0*max( qrmls(kz), -qsw(kz)/dtp ) ! ERM: allow melting to remove all snow
         qsdsv(kz) = xds*csw*xsv  ! LJW & ERM fix - explicit thresholding for dep/sub of snow (MSG AGREES)
         qvsbs(kz) = abs(max( min(qsdsv(kz), 0.0), -qs5dt ))
         qsdpv(kz) =     max(qsdsv(kz), 0.0)

c Conversion qsw->qhw 

         IF ( iglcnvs .lt. 3 ) THEN
           qhcns(kz) = min(0.001*ehscnv*max((qsw(kz)-6.e-4),0.0),qs5dt)   ! LJW & ERM includes ehscnv instead of ehs for hail AC
         ELSE
           IF ( qsacw(kz) .gt. 1.0e-12 .and. qsacw(kz) .gt. qsdpv(kz) .and. temcg(kz) < 0.0 ) THEN
            qhcns(kz) = min(0.001*ehscnv*max((qsw(kz)-6.e-4),0.0),qs5dt)   ! LJW & ERM includes ehscnv instead of ehs for hail AC
           ENDIF
         ENDIF
        ENDIF                     !END OF ACTIVE SNOW GRIDCELLS

        IF (  gqrw(kz)) THEN      !RAIN ACTIVE

c Accretions

         IF (gqcw(kz)) qracw(kz)=min(erw*qcw(kz)*pid4*crw*abs(vtrbar-vtwbar)*(gf3*rwd2+2.*gf2*rwd*cwd+gf1*cwd2),qc5dt)
c         IF ( qracw(kz) .gt. 0.0 ) THEN
c           write(91,*) 'qracw,cx =',qracw(kz),ccn,kz,crw,1.e6*rwd,vtrbar,vtwbar
c           write(91,*) 'vtr: ',vtrbar, (ar*gf4br/6.0)*(rwd**br)* (dnz00/dnz(kz))**(0.5), dnz(kz), ar, br,qrw(kz)     !  xvtr(kz)*(xslop(kz)* (qrw(kz)/(rwdn(kz)*cnor))**(0.25))**br
c         ENDIF
         IF (gqci(kz)) qraci(kz)=min(eri*qci(kz)*pid4*crw*abs(vtrbar-vtibar)*(gf3*rwd2+2.*gf2*rwd*cid+gf1*cid2),qi5dt)
         IF (gqci(kz)) qiacr(kz)=min(eri*qrw(kz)*xacx*cci*abs(vtrbar-vtibar)*(gf6*rwd2+2.*gf5*rwd*cid+gf4*cid2),qr5dt)
         IF (gqsw(kz)) qracs(kz)=min(ers*qsw(kz)*xacx*crw*abs(vtrbar-vtsbar)*(gf6*gf1*swd2+2.*gf5*gf2*swd*rwd+gf4*gf3*rwd2),qs5dt)

c Biggs Freezing

         xrfrz     = (20.0)*(pi**2)*brz*(cwdn/dnz(kz))
         qhfzr(kz) = min(xrfrz*crw*(rwd**6)*(exp(max(-arz*temcg(kz), 0.0))-1.0),qr5dt)

c Evaporation

          ssw = qwv(kz)/qvs(kz)
          xbv = (1.0/(dnz(kz)*wvdf(kz)))
          xav = (elv**2)/(tka(kz)*rw)
          xxcev= 2.0*pi/dnz(kz)
          xce = xxcev*(ssw-1.0)*(1.0/(xav/temgsqr+xbv/qvs(kz)))
          xvevr2 = 0.308*(gamma((5.0+br)/2.0))*(schm(kz)**(1./3.))
     $             * (akvisc(kz)**(-0.5))*(ar**(0.5))*((dnz00/dnz(kz))**(0.25))
          xrv = (xvevr1*(rwd) + xvevr2*(rwd**((3.0+br)/2.0)))
          qvevr(kz) = -1.0*max(min(xce*crw*xrv, 0.0), -qr5dt)

        ENDIF                     !END ACTIVE RAIN GRIDCELLS

c Bergeron process for snow, Conversion into snow of ice
        
        IF ( gqci(kz) ) THEN                       !(ERM thresholding on 12.31.02)
          ib       = min(max(1,int(-temcg(kz))),32)
          qsfi(kz) = qci(kz)/cbtim(ib)
          qsfi(kz) = min(qsfi(kz),qi5dt)

          IF ( gqcw(kz)) THEN                      !(ERM Thresholding 12/31/02)
            cs9      = bsfw*dnz(kz)*(0.001)
            qsfw(kz) = qci(kz)*cs11(ib)*(cs10(ib) + eic*cs9*qcw(kz))
            qsfw(kz) = min(qsfw(kz),qc5dt)
          ENDIF 

          qscni(kz) = min(0.001*esicnv*max((qci(kz)-1.e-3),0.0),qi5dt)   ! ERM includes esicnv instead of esi for snow AC

        ENDIF        !End of Bergeron
        
c Cloud ice initiation (Vapor to pristine ice crystals)

        IF ( ssi .gt. 1.0 ) THEN
          xiint  = (elv**2)/(cp*rw)
          dqisdt = (qwv(kz)-qis(kz))/ (1.0 + xiint*qis(kz)/temgsqr)
          cnnt   = min(cnit*exp(-temcg(kz)*bta1),1.0e+09)
          qiint(kz) = (1.0/dtp)*min((1.0e-12)*cnnt/rho0(kz),0.50*dqisdt)
        ENDIF

c AutoConversion qcw -> qrw
c Autoconversion == 0 --> do standard autoconversion theshold
c Autoconversion == 1 --> Berry 1968 auto conversion for rain with Ferrier (1994) ccn dependence

	  IF( autoconversion .eq. 0 ) THEN
         qdiff  = max((qcw(kz)-qcmincwrn),0.)
	  ELSE
         qccrit = (pi/6.)*((ccw*cwdiap**3)*cwdn)/rho0(kz)
         qdiff  = max((qcw(kz)-qccrit),0.)
	  ENDIF

        IF ( qdiff .gt. 0.0 ) THEN
          argrcnw1  = (1.596e-12)*(1e-6) / (cwdisp*1e-3)
!          argrcnw   = (1.2e-4 + argrcnw1*ccw) / (rho0(kz)*qdiff)
          argrcnw   = 1.2e-4 + ( argrcnw1*ccw) / (rho0(kz)*qdiff)

! version from G Bryan's code  (buggy, missing 1e-3*rho0 in denominator) :
!         argrcnw = (1.2e-4)+(1.596e-12)*ccw*(1e-6)/(cwdisp*qdiff)
         argrcnw = (1.2e-4)+(1.596e-12)*ccw*(1e-6)/(rho0(kz)*1.e-3*cwdisp*qdiff)

          qrcnw(kz) = min(max(rho0(kz)*1e-3*(qdiff**2)/argrcnw,0.0),qc5dt)
c          write(91,*) 'qrcnw,cx =',qrcnw(kz),ccw,kz,1.e3*qcw(kz)
        ENDIF

c      IF (ndebug .gt. 0 ) print*,'dbg = 8'

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c Production    ! MSG i=source terms (all positive), d=sink terms (all negative)

        pqwvi(kz)  =  qvevr(kz) +il5(kz)*( qvsbh(kz) +qvsbs(kz) )
        pqwvd(kz)  =  il5(kz)*( -qhdpv(kz) -qsdpv(kz) -qiint(kz) )

        pqcii(kz)  =  il5(kz)*qiint(kz)
        pqcid(kz)  =  il5(kz)*( -qscni(kz) -qsaci(kz) -qraci(kz) -qsfi(kz)) -qhaci(kz)

        pqcwi(kz)  =  0.0
        pqcwd(kz)  = -(qracw(kz)) -qsacw(kz) -qrcnw(kz) -qhacw(kz) - il5(kz)*qsfw(kz)

        pqrwi(kz)  =  (qracw(kz)) + qrcnw(kz) + (1-il5(kz)) * (qsacw(kz)+qrmlh(kz) + qrmls(kz))+qrshh(kz)
        pqrwd(kz)  =  il5(kz)*(-qiacr(kz) -qhfzr(kz) -qsacr(kz)) -qvevr(kz) -qhacr(kz)

        pqswi(kz)  =  il5(kz)*(qsacw(kz) +qscni(kz) +qsaci(kz) + qsfi(kz) +qsfw(kz)
     $             +  il3(kz)*(qraci(kz) +qiacr(kz))+il2(kz)*qsacr(kz)+qsdpv(kz) )

        pqswd(kz)  = -qhcns(kz) -qhacs(kz) - (1-il5(kz))*qrmls(kz) + il5(kz)*(-qvsbs(kz) -(1-il2(kz))*qracs(kz))

        pqhwi(kz)  =  qhcns(kz) + il5(kz)* ( (1-il3(kz))*(qraci(kz) +qiacr(kz))+(1-il2(kz))*(qsacr(kz) +qracs(kz))
     $             +  qhdpv(kz) + qhfzr(kz)) + qhacr(kz) +qhacw(kz) +qhacs(kz) +qhaci(kz)

        pqhwd(kz)  =  -qrshh(kz) -(1-il5(kz))*qrmlh(kz) - il5(kz)*qvsbh(kz)

c       ptotal(kz) = pqwvi(kz) + pqwvd(kz) + pqcwi(kz) +pqcwd(kz) + pqcii(kz) +pqcid(kz)
c    $             + pqrwi(kz) + pqrwd(kz) + pqswi(kz) +pqswd(kz) + pqhwi(kz) +pqhwd(kz)


c      IF (ndebug .gt. 0 ) print*,'dbg = 9'

c Latent heating from phase changes (except qcw, qci cond, and evap)

        pfrz(kz)   = (1.-il5(kz))*(-qrmlh(kz) - qrmls(kz))+il5(kz)*(qiacr(kz)+qsacr(kz)+qsfw(kz)
     $                                        + qhfzr(kz)+qsacw(kz)+qhacw(kz)+qhacr(kz)-qrshh(kz))
        psub(kz)   = il5(kz)*(qhdpv(kz)-qvsbh(kz) + qiint(kz)+qsdpv(kz)-qvsbs(kz))
        pvap(kz)   = -qvevr(kz)
        ptem(kz)   = (cpi/piz(kz))*(elf*pfrz(kz) + els*psub(kz) + elv*pvap(kz))
        thetap(kz) = thetap(kz) + dtp*ptem(kz)
        theta(kz)  = thetap(kz) + theta0(kz)
!        temg(kz)   = theta(kz)*( pres(kz) / poo ) ** rcp
        temg(kz)   = theta(kz)*(pib(kz) + pin(ix,jy,kz))

c       IF (IEEE_IS_NAN(thetap(kz))) stop 'Z'
c       IF (IEEE_IS_NAN(qci(kz))) stop 'Y'
c       IF (IEEE_IS_NAN(qcw(kz))) stop 'W'
c       IF (IEEE_IS_NAN(qwvp(kz))) stop 'V'

c      IF (ndebug .gt. 0 ) print*,'dbg = 9b'



c Update q for all sources and sinks up until now

        qwvp(kz) = qwvp(kz) + dtp*(pqwvi(kz)+pqwvd(kz))    !updated vapor perturbation
        qcw(kz)  = qcw(kz)  + dtp*(pqcwi(kz)+pqcwd(kz))
        qci(kz)  = qci(kz)  + dtp*(pqcii(kz)+pqcid(kz))
        qrw(kz)  = qrw(kz)  + dtp*(pqrwi(kz)+pqrwd(kz))
        qsw(kz)  = qsw(kz)  + dtp*(pqswi(kz)+pqswd(kz))
        qhw(kz)  = qhw(kz)  + dtp*(pqhwi(kz)+pqhwd(kz))

      IF (ndebug .gt. 0 ) print*,'dbg = 9c'


c Melting of cloud ice   *based upon updated temperature*

       IF( temg(kz) .gt. tfr .and. qci(kz) .gt. 0.0 ) THEN
         qcw(kz)    = qcw(kz) + qci(kz)
         pfrz(kz)   = pfrz(kz) - qci(kz)/dtp
         cc3        = cpi*elf/piz(kz)
         thetap(kz) = thetap(kz) - cc3*qci(kz)
         qci(kz)    = 0.0
       ENDIF

      IF (ndebug .gt. 0 ) print*,'dbg = 9d'



c Homogeneous freezing of cloud water T < -40  *based upon updated temperature*

       IF( temg(kz) .lt. thnuc .and. qcw(kz) .gt. 0.0 ) THEN      
         qci(kz)    = qci(kz) + qcw(kz)
         pfrz(kz)   = pfrz(kz) + qcw(kz)/dtp
         cc3        = cpi*elf/piz(kz)
         thetap(kz) = thetap(kz) + cc3*qcw(kz)
         qcw(kz)    = 0.0
       ENDIF

      IF (ndebug .gt. 0 ) print*,'dbg = 9e'


c---------------------------------------------------------------------
c Modified Straka adjustment (similar to Tao et al. 1989 MWR)
c---------------------------------------------------------------------

c Iterate  adjustment   !Note that Gilmore et al. (2004) found no improvement iterating twice versus once.

        theta(kz)  = thetap(kz) + theta0(kz)
!        temg(kz)   = theta(kz)*( pres(kz) / poo ) ** rcp 
        temg(kz) = theta(kz)*(pib(kz) + pin(ix,jy,kz))
        temcg(kz)  = temg(kz) - tfr                       

      IF (ndebug .gt. 0 ) print*,'dbg = 9f'


      IF ( .false. .and. temg(kz) .ge. tfr ) THEN
      

         qvs(kz)    = 380.*exp(17.27*(temg(kz) - 273.)/(temg(kz) - 36.))/pb(kz)
         product    = (qwv(kz) - qvs(kz)) / (1. + qvs(kz)*f5/(temg(kz) - 36.)**2)

         product    = max(product,-qcw(kz))
         thetap(kz)  = thetap(kz) + gam(kz)*product
         qwvp(kz)  = qwvp(kz) - product
         qcw(kz)  = max(qcw(kz) + product,0.)
      
      
      ELSE
      
       DO itertd = 1,2

c       Update temperature due to melting/freezing of cloud (or adjustment if on 2nd iteration)
       IF ( itertd .gt. 1 ) THEN
        theta(kz) = thetap(kz) + theta0(kz)
!        temg(kz)   = theta(kz)*( pres(kz) / poo ) ** rcp 
        temg(kz)  = theta(kz)*(pib(kz) + pin(ix,jy,kz))
        temcg(kz)  = temg(kz) - tfr
       ENDIF

c Reset temporaries 

        qv0      = qwv(kz) !  ab(kz,lv)
        qvap(kz) = max( (qwvp(kz) + qv0), 0.0 )    !updated total vapor 
        qwv(kz)  = max( 0.0, qvap(kz))
        qwvp(kz) = 0.0
        qcw(kz)  = max( 0.0, qcw(kz) )
        qci(kz)  = max( 0.0, qci(kz) )             !MSG 10/25/03 - Was originally unfairly omitted in 2nd iteration loop

c Calculate Saturation

         ltemq     = nint((temg(kz)-163.15)/fqsat+1.5)
         ltemq     = min(max(ltemq,1),nqsat)
         qvs(kz)   = pqs(kz)*tabqvs(ltemq)
         qis(kz)   = pqs(kz)*tabqis(ltemq)

         IF ( temg(kz) .lt. tfr ) THEN     
           IF( qcw(kz) .ge. 0.0 .and. qci(kz) .eq. 0.0) qss(kz) = qvs(kz)
           IF( qcw(kz) .eq. 0.0 .and. qci(kz) .gt. 0.0) qss(kz) = qis(kz)
           IF( qcw(kz) .gt. 0.0 .and. qci(kz) .gt. 0.0) qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz))/(qcw(kz) + qci(kz))
         ELSE
           qss(kz) = qvs(kz)
         ENDIF

c Calculate super-saturation/under-saturation

         dqcw(kz) = 0.0
         dqci(kz) = 0.0
         dqwv(kz) = ( qwv(kz) - qss(kz) )

c Evaporation and sublimation adjustment

         IF( dqwv(kz) .lt. 0. ) THEN

           IF( qcw(kz) .gt. -dqwv(kz) ) THEN
             dqcw(kz) = dqwv(kz)
             dqwv(kz) = 0.
           ELSE
             dqcw(kz) = -qcw(kz)
             dqwv(kz) = dqwv(kz) + qcw(kz)
           ENDIF
 
           IF( qci(kz) .gt. -dqwv(kz) ) THEN
             dqci(kz) = dqwv(kz)
             dqwv(kz) = 0.
           ELSE
             dqci(kz) = -qci(kz)
             dqwv(kz) = dqwv(kz) + qci(kz)
           ENDIF

           qwvp(kz)   = qwvp(kz) - ( dqcw(kz) + dqci(kz) )                         !cumulative change in water vapor due to evap/sub
           qcw(kz)    = qcw(kz) + dqcw(kz)                                         !cumulative change in qcw due to evap
           qci(kz)    = qci(kz) + dqci(kz)                                         !cumulative change in qci due to sub
           pvap(kz) = pvap(kz) + dqcw(kz)/dtp
           psub(kz) = psub(kz) + dqci(kz)/dtp
           thetap(kz) = thetap(kz) + cpi/piz(kz)*(elv*dqcw(kz)+els*dqci(kz))       !Update TH' due to evap/sub
         ENDIF

c Condensation/deposition

         IF( dqwv(kz) .ge. 0. ) THEN

!         write(91,*) 'satadj: mgs,iter = ',kz,itertd,dqwv(kz),qss(kz),qwv(kz)

           fracl(kz) = 1.0
           fraci(kz) = 0.0

           IF ( temg(kz) .lt. tfr .and. temg(kz) .gt. thnuc ) THEN
             fracl(kz) = max(min(1.,(temg(kz)-233.15)/(20.)),0.0)
             fraci(kz) = 1.0-fracl(kz)
           ENDIF

           IF ( temg(kz) .le. thnuc ) THEN
             fraci(kz) = 1.0
             fracl(kz) = 0.0
           ENDIF

           fraci(kz) = 1.0-fracl(kz)

           gamss = (elv*fracl(kz) + els*fraci(kz)) / (piz(kz)*cp)

           IF ( temg(kz) .lt. tfr ) THEN

             IF (qcw(kz) .ge. 0.0 .and. qci(kz) .le. 0.0 ) dqvcnd(kz) = dqwv(kz)/(1. + cqv1*qss(kz)/ ((temg(kz)-cbw)**2))
             IF( qcw(kz) .eq. 0.0 .and. qci(kz) .gt. 0.0 ) dqvcnd(kz) = dqwv(kz)/(1. + cqv2*qss(kz)/ ((temg(kz)-cbi)**2))

             IF ( qcw(kz) .gt. 0.0 .and. qci(kz) .gt. 0.0 ) THEN           !MSG- Note that we are not thresholding qcwmin, qcimin.
                                                                           !MSG- In other words, negligible qi, qw are dealt with
               cdw = caw*piz(kz)*tfrcbw/((temg(kz)-cbw)**2)
               cdi = cai*piz(kz)*tfrcbi/((temg(kz)-cbi)**2)
               denom1 = qcw(kz) + qci(kz)
               denom2 = 1.0 + gamss*(qcw(kz)*qvs(kz)*cdw + qci(kz)*qis(kz)*cdi) / denom1
               dqvcnd(kz) =  dqwv(kz) / denom2
             ENDIF

           ENDIF

           IF ( temg(kz) .ge. tfr ) dqvcnd(kz) = dqwv(kz)/(1. + cqv1*qss(kz) / ((temg(kz)-cbw)**2))

           dqcw(kz) = dqvcnd(kz)*fracl(kz)
           dqci(kz) = dqvcnd(kz)*fraci(kz)

           qwvp(kz) = qwvp(kz) - ( dqvcnd(kz) )                                   !cumulative change in vapor, accounting now for cond/dep
           qcw(kz)  = qcw(kz) + dqcw(kz)                                          !cumulative change in qcw, accounting now for cond
           qci(kz)  = qci(kz) + dqci(kz)                                          !cumulative change in qci, accounting now for dep
           pvap(kz) = pvap(kz) + dqcw(kz)/dtp
           psub(kz) = psub(kz) + dqci(kz)/dtp
           thetap(kz) = thetap(kz) + (elv*dqcw(kz) + els*dqci(kz)) / (piz(kz)*cp) !cumulative change in TH', accounting now for cond/dep

          ENDIF        ! End of condensation/deposition within sat.adj

         ENDDO    ! Sat adj loop
       ENDIF !! ( .false. .and. temg(kz) .ge. tfr )

      IF (ndebug .gt. 0 ) print*,'dbg = 9g'

c  MSG - Before leaving the microphysics, we determine which boxes will be affected by fallout
c  Saves computation time later!  (10/25/03)
c
c  Update booleans indicating whether enough precip remains to fallout

         gcnt(kz) = 0
         gqcw(kz) = .false.
         gqci(kz) = .false.
         gqrw(kz) = .false.
         gqsw(kz) = .false.
         gqhw(kz) = .false.

         IF ( qcw(kz) .gt. qcmin  ) gqcw(kz) = .true.
         IF ( qci(kz) .gt. qimin  ) gqci(kz) = .true.
         IF ( qrw(kz) .gt. qrmin  ) gqrw(kz) = .true.
         IF ( qsw(kz) .gt. qsmin  ) gqsw(kz) = .true.
         IF ( qhw(kz) .gt. qhmin  ) gqhw(kz) = .true.
         IF ( gqcw(kz) .or. gqci(kz) .or. gqrw(kz) .or. gqsw(kz) .or. gqhw(kz) ) gcnt(kz) = 1

c COMPUTE SPLIT-EXPLICIT TIME STEP USING UPDATED qrw and qhw (fastest falling species)

         IF ( gqrw(kz) ) THEN 
           rwd = xslop(kz)* (qrw(kz)/(rwdn(kz)*cnor))**(0.25)
           vtrbar = xvtr(kz)*(rwd**br)
           IF ( nprock > 1 ) THEN ! limit fall speeds for vertical tiling to avoid time splitting.
             vtrbar = Min( vtrbar, 0.85*dzc(kz)/dtp )
           ENDIF
           vtrdzmax  = max(vtrdzmax,dtp*dzc(kz)*vtrbar)
         ENDIF
         IF ( gqrw(kz) .or. gqhw(kz) ) THEN 
           hwd = xslop(kz)* (qhw(kz)/(hwdn(kz)*cnoh))**(0.25)
           vthbar = xvth1*xvth2(kz)*((hwdn(kz)*hwd)**(0.5))
           IF ( nprock > 1 ) THEN ! limit fall speeds for vertical tiling to avoid time splitting.
             vthbar = Min( vthbar, 0.85*dzc(kz)/dtp )
           ENDIF
           vtrdzmax  = max(vtrdzmax,dtp*dzc(kz)*vthbar)
         ENDIF

c DETERMINE KZ LAYER FOR FALLOUT

        IF ((gcnt(kz) .eq. 1) .and. (kzfstart .eq. -1)) kzfstart = kz

        IF (gcnt(kz) .eq. 1) THEN
          kzfend   = kz
          gcntcolumn = 1              !"Switch on" the column
        ENDIF

      IF (ndebug .gt. 0 ) print*,'dbg = 9h'


c       IF (ndebug .eq. 1) print *, ix,jy,kz,'Fall   gcnt',gcnt(kz),' col:',gcntcolumn,'gqcw:',gqcw(kz),
c    $                              'gqci:',gqci(kz),'gqrw:',gqrw(kz),  'gqsw:',gqsw(kz),'gqhw:',gqhw(kz)


c END OF MICROPHYSICS FOR GRIDBOX --> INSERT NEW VALUES OF THETA AND HYDROMETEORS 
c                                     (MSG - q values important to update here since might have fallen below qmin threshold
c                                            and thus would not be updated properly after the fallout routine)

         an(ix,jy,kz,lt) = theta0(kz)+ thetap(kz)
         an(ix,jy,kz,lv) = qwv(kz) + qwvp(kz)
         an(ix,jy,kz,li) = qci(kz) + Min( 0.0, an(ix,jy,kz,li) )
         an(ix,jy,kz,lc) = qcw(kz) + Min( 0.0, an(ix,jy,kz,lc) )
         an(ix,jy,kz,lr) = qrw(kz) + Min( 0.0, an(ix,jy,kz,lr) )
         an(ix,jy,kz,ls) = qsw(kz) + Min( 0.0, an(ix,jy,kz,ls) )
         an(ix,jy,kz,lh) = qhw(kz) + Min( 0.0, an(ix,jy,kz,lh) )

        ENDIF                                   ! Done with micro in box/gridcell


      IF (ndebug .gt. 0 ) print*,'dbg = 10'

       IF ( allocated( microrates ) ) THEN
         microrates(ix,jy,kz,1) = pfrz(kz)*elf/(cp*piz(kz))*3600.
         microrates(ix,jy,kz,2) = psub(kz)*els/(cp*piz(kz))*3600.
         microrates(ix,jy,kz,3) = pvap(kz)*elv/(cp*piz(kz))*3600.
         microrates(ix,jy,kz,4) = -qvevr(kz)*elv/(cp*piz(kz))*3600.
         microrates(ix,jy,kz,5) = -qrmlh(kz)*elf/(cp*piz(kz))*3600.
!         microrates(ix,jy,kz,6) = pcond2(kz)
       ENDIF

       ENDDO                                    ! End of main k loop for Microphysics


c Vertical fluxes of hydrometers
                                         ! LJW 12/31/02:  MAJOR CHANGE -> TIME SPLIT RECOMPUTES FLUXES EACH TIME STEP 
                                         ! MSG 1/1/03  :  USES NEW FALL VELOCITIES
                                         ! MSG 10/25/03:  Only computed for those gridcells needing it

       IF (gcntcolumn .gt. 0) THEN       ! Precip present in column!

        nrnstp = 1 + aint(2.*vtrdzmax)   ! MSG Define number of small steps based upon
                                         ! fastest falling particles (rain or hail/graupel)
                                         ! nrnstp for 1 + aint(2*0.45) = 1 (use original timestep)
                                         ! whereas 1+aint(2*0.55) = 2 (requiring dtp be split into
                                         ! two smaller timesteps.) Thus if courant number, vtrdzmax>=0.5,
                                         ! then automatically uses 2 small steps (0.25 courant each).
                                         ! If vtrdzmax=1.0, then uses 3 small steps (0.33 courant each).
                                         ! If vtrdzmax=1.5, then uses 4 small steps (0.25 courant, Etc.)
                                         ! Thus the actual fallout courant ranges 0.25 to 0.5

        IF ( nprock > 1 ) nrnstp = 1
        
        dtsplit = dtp / float( nrnstp )  ! MSG Define small timestep (sec)

        precip(ix,jy,1) = 0.0            ! Set instantaneous rates to zero for time averaging the precip rate
        precip(ix,jy,2) = 0.0            ! due to time splitting

#ifdef MPI
        kzfstart = max(kzfstart-1,    1)    ! Added in 3RD order vertical advection with limiting (LJW, 06/18/04)
        kzfend   = min(kzfend+1  ,  kze)    ! Because the domain of influence is now 1 grid pt higher and lower, 
                                            ! we need to expand the calculation of fluxes by one grid point.
#else
        kzfstart = max(kzfstart-1,   1)    ! Added in 3RD order vertical advection with limiting (LJW, 06/18/04)
        kzfend   = min(kzfend+1,  nz-1)    ! Because the domain of influence is now 1 grid pt higher and lower, 
                                           ! we need to expand the calculation of fluxes by one grid point.
#endif

        DO inrnstp = 1,nrnstp            ! MSG BEGIN OF TIMESPLITTING

         pflux(:) = 0.0                  ! LJW 06/24/04:  Fixed another time-split bug, fluxes should be zero at beginning of each iteration.
         cflux(:) = 0.0
         rflux(:) = 0.0
         sflux(:) = 0.0
         hflux(:) = 0.0

         DO kz = kzfstart,kzfend           ! Only compute fall velocity where there is precip (between kzfstart and kzfend)
                                           ! Compute fall velocities & fluxes for the split step
                                           ! Thresholding of mixing ratio is needed to prevent arithmetic underflow (LJW, 01/2/03).
c Cloud water 

          IF ( gqcw(kz) ) THEN
            cwmas     = min( max(qcw(kz)*rho0(kz)/ccw,cwmasn),cwmasx )
            cwd       = (cwmas*cwc1)**c1f3
            vtwbar    = (ar*(cwd**br))*(rho00/rho0(kz))**0.5               ! MSG Extrapolation of Liu and Orville
            crc       = vtwbar * dzc(kz) * dtsplit
            
c           IF( kz .GT. 1 .AND. kz .LT. nz-1) THEN
c            cflux(kz) = -         (an(jy,kz+1,ix,lc) - 5.*an(ix,jy,kz,lc) - 2.*an(jy,kz-1,ix,lc)) / 6.0
c    $                   - crc*    (                  -    an(ix,jy,kz,lc) +    an(jy,kz-1,ix,lc)) / 2.0
c    $                   - crc*crc*(an(jy,kz+1,ix,lc) - 2.*an(ix,jy,kz,lc) +    an(jy,kz-1,ix,lc)) / 6.0
c
c            qD        = an(jy,kz-1,ix,lc)                                 ! Thuburn monotonic limiter
c            qC        = an(jy,kz  ,ix,lc)
c            qU        = an(jy,kz+1,ix,lc)
c            qinmin    = min(qD, qC)
c            qinmax    = max(qD, qC)
c            qi        = max(cflux(kz), qinmin)
c            qi        = min(qi,        qinmax)
c            qoutmin   = (qC - (1.-crc)*max(qC,qU)) / (crc+eps)
c            qoutmax   = (qC - (1.-crc)*min(qC,qU)) / (crc+eps)
c            qi        = max(qi, qoutmin)
c            cflux(kz) = min(qi, qoutmax)
c
c           ELSE  

             cflux(kz) = an(ix,jy,kz,lc)

c           ENDIF  

            cflux(kz) = rho0(kz) * cflux(kz) * vtwbar

          ENDIF

c Cloud ice 

          IF ( gqci(kz) .and. temcg(kz) .le. 0.0 ) THEN
            cci       = max(min(cnit*exp(-temcg(kz)*bta1),1.e+09),1.0)
            cimas     = min( max(qci(kz)*rho0(kz)/cci,cimasn),cimasx )
            cid       = min(max(16.7*(cimas**(0.5)),1.e-5),1.e-2)          ! MSG thresholding not really needed since applied cimasn already
            vtibar    = (cs*(cid**ds))*(rho00/rho0(kz))**0.5               ! MSG Snow fall velocity equation applied for cloudice.  (Gilmore et al. 2004)
            cri       = vtibar * dzc(kz) * dtsplit
            
c           IF( kz .GT. 1 .AND. kz .LT. nz-1) THEN
c            pflux(kz) = -         (an(jy,kz+1,ix,li) - 5.*an(jy,kz,ix,li) - 2.*an(jy,kz-1,ix,li)) / 6.0
c    $                   - cri*    (                  -    an(jy,kz,ix,li) +    an(jy,kz-1,ix,li)) / 2.0
c    $                   - cri*cri*(an(jy,kz+1,ix,li) - 2.*an(jy,kz,ix,li) +    an(jy,kz-1,ix,li)) / 6.0
c
c            qD        = an(jy,kz-1,ix,li)                                 ! Thuburn monotonic limiter
c            qC        = an(jy,kz  ,ix,li)
c            qU        = an(jy,kz+1,ix,li)
c            qinmin    = min(qD, qC)
c            qinmax    = max(qD, qC)
c            qi        = max(pflux(kz), qinmin)
c            qi        = min(qi,        qinmax)
c            qoutmin   = (qC - (1.-cri)*max(qC,qU)) / (cri+eps)
c            qoutmax   = (qC - (1.-cri)*min(qC,qU)) / (cri+eps)
c            qi        = max(qi, qoutmin)
c            pflux(kz) = min(qi, qoutmax)
c
c           ELSE  

             pflux(kz) = an(ix,jy,kz,li)

c           ENDIF  

            pflux(kz) = rho0(kz) * pflux(kz) * vtibar

          ENDIF

c Rain

          IF ( gqrw(kz) ) THEN 
            rwd       = xslop(kz)* (qrw(kz)/(rwdn(kz)*cnor))**(0.25)
            vtrbar    = xvtr(kz)*(rwd**br)
            crr       = vtrbar * dzc(kz) * dtsplit
            
c           IF( kz .GT. 1 .AND. kz .LT. nz-1) THEN
c            rflux(kz) = -         (an(jy,kz+1,ix,lr) - 5.*an(jy,kz,ix,lr) - 2.*an(jy,kz-1,ix,lr)) / 6.0
c    $                   - crr*    (                  -    an(jy,kz,ix,lr) +    an(jy,kz-1,ix,lr)) / 2.0
c    $                   - crr*crr*(an(jy,kz+1,ix,lr) - 2.*an(jy,kz,ix,lr) +    an(jy,kz-1,ix,lr)) / 6.0
c
c            qD        = an(jy,kz-1,ix,lr)                        	! Thuburn monotonic limiter
c            qC        = an(jy,kz  ,ix,lr)
c            qU        = an(jy,kz+1,ix,lr)
c            qinmin    = min(qD, qC)
c            qinmax    = max(qD, qC)
c            qi        = max(rflux(kz), qinmin)
c            qi        = min(qi,        qinmax)
c            qoutmin   = (qC - (1.-crr)*max(qC,qU)) / (crr+eps)
c            qoutmax   = (qC - (1.-crr)*min(qC,qU)) / (crr+eps)
c            qi        = max(qi, qoutmin)
c            rflux(kz) = min(qi, qoutmax)
c
c           ELSE  

             rflux(kz) = an(ix,jy,kz,lr)

c           ENDIF  

            rflux(kz) = rho0(kz) * rflux(kz) * vtrbar

          ENDIF

c Snow
          IF ( gqsw(kz) ) THEN
            swd       = xslop(kz)* (qsw(kz)/(swdn(kz)*cnos))**(0.25)
            vtsbar    = xvts(kz)*(swd**ds)
            crs       = vtsbar * dzc(kz) * dtsplit
            
c           IF( kz .GT. 1 .AND. kz .LT. nz-1) THEN
c            sflux(kz) = -         (an(jy,kz+1,ix,ls) - 5.*an(jy,kz,ix,ls) - 2.*an(jy,kz-1,ix,ls)) / 6.0
c    $                   - crs*    (                  -    an(jy,kz,ix,ls) +    an(jy,kz-1,ix,ls)) / 2.0
c    $                   - crs*crs*(an(jy,kz+1,ix,ls) - 2.*an(jy,kz,ix,ls) +    an(jy,kz-1,ix,ls)) / 6.0
c
c            qD        = an(jy,kz-1,ix,ls)                        	! Thuburn monotonic limiter
c            qC        = an(jy,kz  ,ix,ls)
c            qU        = an(jy,kz+1,ix,ls)
c            qinmin    = min(qD, qC)
c            qinmax    = max(qD, qC)
c            qi        = max(sflux(kz), qinmin)
c            qi        = min(qi,        qinmax)
c            qoutmin   = (qC - (1.-crs)*max(qC,qU)) / (crs+eps)
c            qoutmax   = (qC - (1.-crs)*min(qC,qU)) / (crs+eps)
c            qi        = max(qi, qoutmin)
c            sflux(kz) = min(qi, qoutmax)
c
c           ELSE  

             sflux(kz) = an(ix,jy,kz,ls)

c           ENDIF  

            sflux(kz) = rho0(kz) * sflux(kz) * vtsbar

          ENDIF

c Hail/graupel

          IF ( gqhw(kz) ) THEN 
            hwd       = xslop(kz)* (qhw(kz)/(hwdn(kz)*cnoh))**(0.25)
            vthbar    = xvth1*xvth2(kz)*((hwdn(kz)*hwd)**(0.5))
            crh       = vthbar * dzc(kz) * dtsplit
            
c           IF( kz .GT. 1 .AND. kz .LT. nz-1) THEN
c            hflux(kz) = -         (an(jy,kz+1,ix,lh) - 5.*an(jy,kz,ix,lh) - 2.*an(jy,kz-1,ix,lh)) / 6.0
c    $                   - crh*    (                  -    an(jy,kz,ix,lh) +    an(jy,kz-1,ix,lh)) / 2.0
c    $                   - crh*crh*(an(jy,kz+1,ix,lh) - 2.*an(jy,kz,ix,lh) +    an(jy,kz-1,ix,lh)) / 6.0
c
c            qD        = an(jy,kz-1,ix,lh)                        	! Thuburn monotonic limiter
c            qC        = an(jy,kz  ,ix,lh)
c            qU        = an(jy,kz+1,ix,lh)
c            qinmin    = min(qD, qC)
c            qinmax    = max(qD, qC)
c            qi        = max(hflux(kz), qinmin)
c            qi        = min(qi,        qinmax)
c            qoutmin   = (qC - (1.-crh)*max(qC,qU)) / (crh+eps)
c            qoutmax   = (qC - (1.-crh)*min(qC,qU)) / (crh+eps)
c            qi        = max(qi, qoutmin)
c            hflux(kz) = min(qi, qoutmax)
c
c           ELSE  

             hflux(kz) = an(ix,jy,kz,lh)

c           ENDIF  

            hflux(kz) = rho0(kz) * hflux(kz) * vthbar

          ENDIF

         ENDDO                              ! End computing fallout fluxes for the kzfstart-kzfend layer

c Update hydrometeors and always start 1 zone lower than precip level in order to get fallout from above

         kzfstart = max(kzfstart-1, 1)     

        IF (ndebug .eq. 1) print *, ix, jy, 'kzfstart,kzfend for fallout:',kzfstart,kzfend

         DO kz = kzfstart,kzfend
          dtz1            = dzc(kz)*dtsplit/db(kz)
          an(ix,jy,kz,li) = an(ix,jy,kz,li) + dtz1*(pflux(kz+1)-pflux(kz))
          an(ix,jy,kz,lc) = an(ix,jy,kz,lc) + dtz1*(cflux(kz+1)-cflux(kz))
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + dtz1*(rflux(kz+1)-rflux(kz))
          an(ix,jy,kz,ls) = an(ix,jy,kz,ls) + dtz1*(sflux(kz+1)-sflux(kz))
          an(ix,jy,kz,lh) = an(ix,jy,kz,lh) + dtz1*(hflux(kz+1)-hflux(kz))
         ENDDO

c------------------------------------------------------------------------------
c PRECIPITATION 
c Compute instantaneous precip rate on grid and accumulated precip on the moving grid
c    den of air    mix rat.   fall velo  x  1/den of water  x  mm / meter
c     kg /m^3   x    g/g    x    m/s     x     m^3 / kg     x  1000 / m    = precip rate (mm/s)

c Average the instantaneous precipitation rate over time to account for time split fallout.
        
        IF ( myprock == 1 ) THEN
         precip(ix,jy,1) = precip(ix,jy,1) + (rflux(1) / 1000.) * 1000. / float(nrnstp)
         precip(ix,jy,2) = precip(ix,jy,2) + (hflux(1) / 1000.) * 1000. / float(nrnstp)
        ENDIF

        ENDDO                                                      ! End timesplit loop

	ENDIF                                                          ! END FALLOUT

c Accumulate precipitation based on the the average instantaneous precipitation rate

      IF ( myprock == 1 ) THEN
      precip(ix,jy,3) = precip(ix,jy,3) + precip(ix,jy,1)*dtp
      precip(ix,jy,4) = precip(ix,jy,4) + precip(ix,jy,2)*dtp
      ENDIF

c  At all grid points, check for minimum values and cascade diminishing values of mixing ratio downscale (Proctor 1987).  
c  Added by MSG 09/14/03, moved here by LJW & ERM 03/22/04

#ifdef MPI
       kzb = 1
       kze = ktile
       if( kzend .eq. nzend) kze = kzend - kzbeg
       
       do kz = kzb,kze
#else
       DO kz = 1,nz-1
#endif

          IF ( an(ix,jy,kz,lr) .le. qrmin ) THEN
           an(ix,jy,kz,lc) = an(ix,jy,kz,lc) + an(ix,jy,kz,lr)
           an(ix,jy,kz,lr) = 0.0
          ENDIF

          IF ( an(ix,jy,kz,lc) .le. qcmin ) THEN
           an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
           an(ix,jy,kz,lc) = 0.0
          ENDIF

          IF ( an(ix,jy,kz,lh) .le. qhmin ) THEN
           an(ix,jy,kz,ls) = an(ix,jy,kz,ls)+ an(ix,jy,kz,lh)
           an(ix,jy,kz,lh) = 0.0
          ENDIF
 
          IF ( an(ix,jy,kz,ls) .le. qsmin ) THEN
           an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls) ! ERM: put snow into vapor instead of cloud ice (can get residual cloud ice at T > 0, which is weird)
           an(ix,jy,kz,ls) = 0.0
          ENDIF

          IF ( an(ix,jy,kz,li) .le. qimin ) THEN
           an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
           an(ix,jy,kz,li) = 0.0
          ENDIF

        ENDDO  !! kz

       ENDDO    !for ix loop

       ENDDO    !for jy loop

#ifdef MPI
      if(debug_mpi) write(0,*) "LFO_ICE: EXITING SUBROUTINE, my_rank=",my_rank
#else
      IF ( ndebug .ge. 1 ) write(0,*) "LFO_ICE: EXITING SUBROUTINE"
#endif

       RETURN
       END
c
c  END OF SUBROUTINE
c
