#define DEBUG 0
! #DEFINE HSTR 0 ! to be used for future horizontal stretching option
! CVD$F SKIP
!#######################################################################
!  SUBROUTINE SLPM
!
!  written by ERM      
!c
!
!
! 'phi2' for initial solution, copied into 'phi' which is used for 
!   propagation (values changed by lightning channels)
! 'pot2' was same as 'phi2' but extended one point below ground.
! 'phi2w' and 'pot2w' were set up so that (1,1,1) corresponds (1,1,1) on the
!  storm grid.
!
!  so now:
!    pot2 -> phi
!    pot2w -> phi0
!    f -> rhs
!    fw -> rhs0
!    phi -> phip  (phi-propagate)
!
!   can get rid of phi2 and phi2w
!
! phi, phi0, and phi2 are all the same array expressed differently (pot2)
!       they are all the 'background' (i.e. preflash) potential
! phip (phi-propagate) is a copy of phi that can be updated with 
!       channel potentials
!
!#######################################################################
!##################################################################
!

      subroutine mudslpmmpi      &
     &  (time_real,time,nstep1,iliter,lgtstp,pchnl,nchnl,ldumpsp      &
     &  ,nor,nx,ny,nz,NLx,NLy,NLz,pot      &
     &  ,icldtop,zcldtop      &
     &  ,icgyn,icgpn,ebrkds      &
     &  ,dtp,dx1,dy1,dz1,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
     &  ,dzmax, istretch      &
     &  ,gxt,gyt,gzt,dxx,dyy,dzz      &
     &  ,ixst,jyst,kzst,outname,nkzmn,icheckdims,      &
     &   nnx,nny,nnz,nnzt,      &
     &   iixp,jjyq,kkzr,iiex,jjey,kkez,      &
     &   nbw,nbe,nbs,nbn,nbz2,  &
!     &   llwork,      &
     &   nxs,nys,nzs,nb,nxb,nyb,nzbt,nzb,nxyz,nwslv,      &
     &   phi,phi0,phi2,rhs,rhs0,rho,rhow,mask,phip,      &
     &   numic,numcgn,numcgp,rhssor,rhs0sor,pickmask,      &
     &   iunit0, bcx, bcy)
     

! Note: the end of this first #if is at the very end of the subroutine.
#if defined (MPI) && defined (BOXMG)

!       USE GRID_MODULE
       USE ELEC_MODULE, only :                         intfg,iseed, igrid,         &
                        intfgmpi,       &
                        mdel,           &
                        ibal,           &
                        cgthres,        &
                        dslight, dslightz,        &
                        l2nde,          &
                        nohold,         &
                        maxpicks,       &
                        multipickmod,   &
                        pickmaskwidth,  &
                        cgfr,           &
                        soreps,         &
                        intsormpi,      &
                        lightintx,      &
                        lightintz, overvolt, elgtfdel, eint, ieint, feint, fecrit, fefac, &
                        dlzlgt, dslightz, z1d2lgt, ltgeta, nliter, isa, rhofairlgt,       &
                        BMG_iPARMSlg, nlxlg, nlylg, nlzlg, ngxlg, ngylg, ngzlg,           &
                        igslg, jgslg, kgslg, BMG_rPARMSlg, BMG_IOFLAGlg, Qlg, Qflg,       &
                        BMG_rWORKlg, NFlg, NCbmglg, SOlg, NSOlg, NSORlg, NCIlg,           &
                        BMG_iWORKlg, NOGlg, NOGclg, &
                        id_BMG3_MAX_ITERSdg, BMG_pWORKlg, BMG_iWORK_PLlg, NBMG_iWORK_PLlg, &
                        BMG_rWORK_PLlg, NBMG_rWORK_PLlg, BMG_iWORK_CSlg, NBMG_iWORK_CSlg,  &
                        BMG_rWORK_CSlg, NBMG_rWORK_CSlg, pMSGlg, pMSGSOlg, &
                        BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg, NMSGrlg,  &
                        id_BMG3_SETUPdg, nzmax1lgt, efracinitwire, irand, igslg0, jgslg0, kgslg0, &
                        ip_RESdg, ip_SORdg, ip_CIdg, ip_iGdg, ip_MSGdg, NMSGilg, ip_MSG_BUFdg, &
                        x1d2lgt, y1d2lgt, &
                        hstretch

       USE COMMASMPI_MODULE

      implicit none

#ifdef MPI
      INCLUDE "mpif.h"

#ifdef BOXMG
!#include        "BMG_constants_f77.h"
!      INCLUDE   'BMG_workspace_f77.h'
!      INCLUDE   'BMG_parameters_f77.h'
#endif

#endif

!#######################################################################
            
!##################################################################
!  PURPOSE: To generate a lightning channel using the Dielectric
!       Breakdown Model of L. Niemeyer, L. Pietronero, and 
!       H.J. Wiesmann (Phys. Rev. Let., 1984, v. 52, p. 1033)
!
!
!  TO DO list:  vertical refinement on stretched grid
!
!##################################################################
!
! 2006.9:  version 2.6:  Allow for 2-D grid and use staggered E for emag
!
! 2005.10.3: Added input variable cgfr to control whether ibal=2 resets
!            the channel potential along the channel or not.  (A value 
!            of cgfr < 1.0 turns this off an cgfrac = cgfr for the fraction
!            of charge retained on the lower channels.)
!
! 2005.05.25  CHANGED the grid refinement so that the smaller boxes fit inside
!    the larger boxes.  This will preserve the total charge, but grid 
!    points will not be co-located with any coarse grid points.  This
!    requires that dxl = n*ds (integer multiple).
!      In other words, now do
!
!     |   *   |   *   |   *   |
! to
!     | *   * | *   * | *   * |  (igrid=1)
!
!   instead of 
! 
!     |   *   |   *   |   *   |
! to
!     |   * | * | * | * | *   | (igrid=0)
!
!  2007.03.22: erm; changed back to igrid = 0.
!
!  This also makes putting charge back on the coarse grid much simpler!
!
! 2005.01.14: Try normalizing phistart hi and lo to a length of 1km by adding
!             or subtracting the quantity  1.e3*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
!             to phistart: 
!         phistlow = phistart - 5.e2*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
!         phisthi  = phistart + 5.e2*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
!
!         Also change count of one-sided growth to be normalized to a dslight=500m:
!
!            IF ( ibal .eq. 2 .and. istopcnt(iend) .gt. 50*500./ds ) THEN
!          
!         so that a flash can go about the same total length without 
!         growth on the other end
!
!
! 8.2004: some tests with setting channel potential to zero upon becoming a CG
!         if ibal = 2
!
! 6.30.2004: Fixed a bug in the final calculations of the channel charge. 
!            buggy version had
!                IF ( k1 .gt. 1 .and. k1 .lt. kz-2 ) THEN
!         changed to
!                IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN
!         Didn't see an obvious problem with this until trying to rerun
!         the WK Us=50 case.  It could have been a culprit in causing 
!         a sizable fraction of a channel not to be counted, sometimes resulting
!         in a retry of the flash.  It would also affect the actual flash 
!         charge.
!           The problem appeared after a change was made to the z-limit of
!         the initiation point list, which was change from kz=1,nzs-4 to 
!         kz=1,k2 where k2 is the point just below cloud top.  This left kz
!         defined as k2+1, so the channels above k2-2 were being cut off, which
!         is a big problem for flashes that initiate close to cloud top, as 
!         they apparently do in the Us=50 WK supercell case.
!
!         Also removed 'the top' condition of numchan(iend).gt.numchan(3-iend) 
!         for ibal=2.  Not much point to that, really.
!
! 5.06.2004: For ibal=2, added a constraint that the channel reference potential
!            (phistart) cannot be shifted beyond the values of the first 
!            two channel points.
!
! 4.17.2004: Found strange logic in the ibal=5 scheme.  It wasn't really
!            doing ihold like it should and was not stopping the flash
!            if it couldn't force growth on the less-charged end.  The
!            changes cause more total steps in the process, but those are
!            pretty cheap compared to the poisson solver.  Also causes a higher
!            flash rate for the same fecrit
!
!
! 4.8.2004: increased the number of times that charge balance is checked
!           for ibal=5.
!      
!            was:
!      IF ( (( nstep/5*5 .eq. nstep .or. nstep .le. 3 )  .and.
!
!           now:
!      IF ( (( nstep/4*4 .eq. nstep .or. nstep .le. 30 )  .and.
!
!
! 2.14.2004: Two bug fixes: cofz was not set for k=0 and k=1, causing 
!            possible problems with strzsor.  Other bug was incorrect
!            setting of chandist when a new channel point is added.
!
! 1.27.2004: fixed a bug in the calculation of Ez in the initiation point
!            list generation code.  It thought that dz was always dzmax
!            (wrong subscript in z1d2(kz,3))
!
! There was also a bug in the calculation of qtot for the charge balance
! checks (though not at the very end where charge is put into the 'coarse' grid).
! The charge density was being scaled by the stretch factor, but the wrong 
! index was being used (k2 instead of k1).  Not a huge problem, but good to
! have it fixed.
!
! 01.05.2004 Update/bug fixes for ibal=2 option, which seems to be behaving
!            better since other bugs were fixed.  Apparently need to use a 
!            relatively low internal field (1kV/m or less).  Added full flash SOR
!            after adjusting channel potential.  Had a bug that caused internal
!            field to be applied only to the first points, with zero drop after that.
!            Applied the channel potential criterion for ground flashes, but
!            it seems questionable to do so in this case, since the channel is
!            artificially widened, thus diluting the value of the channel potential.
!            Must think about this some more!
!
!            Also leave charge on CG channels as for ibal=5.
!
!
!
! 07.25.2003 sam.spark.ion.2.5.f: update to allow unstretched grid with 
!            dslight < dzl
!
! 01.31.2003 update to documentation: ibal=5 acts like ibal=4 and also
!            feeds back charge from the ground leaders (zeroing points that
!            have the wrong sign of charge, which can happen on the part of the 
!            channel that is extended to ground when a threshold is used.)
!
!            Both ibal=4 and ibal=5 add a second criterion for deciding ground
!            flashes: the potential of the downward leader tip must still have
!            the same sign it started with (i.e. a negative leader tip must have
!            a negative channel potential.)
!
!            Also removed some old commented code to clean things up a bit.
!
! 10.22.2002 version 2.4c:
!            Added ibal = 4 option that is like ibal=3 but does not stop
!            an unbalanced flash if the the downward channels have more 
!            charge.  
!            planned to add: preserve charge on downward CG channels
!
!              In sll, ibal=4 will limit charge on ground channels to 
!            -0.25*(charge on upper channels) (1/31/2003 Now only ibal=5 does this)
!
! 9.16.2002 increased search range for initiation points to get closer
!           to the edge of the domain:
!
!      DO jy = nbs+3,nny-nbn-3
!      DO ix = nbw+3,nnx-nbe-3
!  becomes
!
!      DO jy = nbs+2,nny-nbn-1
!      DO ix = nbw+2,nnx-nbe-1
!
!  also
!     :         i1+i .gt. nbw+1 .and. i1+i .lt. nxs+nbw-1 .and.
!     :         j1+j .gt. nbs+1 .and. j1+j .lt. nys+nbs-1 .and.
!  to
!     :         i1+i .gt. nbw+1 .and. i1+i .le. nxs+nbw .and.
!     :         j1+j .gt. nbs+1 .and. j1+j .le. nys+nbs .and.
!
! 8.26.2002 Had a problem with the growth on the upward channels of a CG
!           also getting to ground.  Suppressed as follows: 
!            1. the flash is stopped if the upper end gets within 0.5*cgthres
!               of the ground
!            2. the upper branches are allowed to grow a maximum of 50 extra 
!               channel segments
!           Also (mostly) fixed an occasional problem where the propagation
!           thresholds of both ends could go below fecrit.  Now, when one end
!           lowered, the other end is raised back toward fecrit (to a maximum
!           of fecrit).
!
! 8.20.2002 Originally left all charge on channel for CGs when ibal=3 but
!           now changed back to the original scheme of setting charge to 
!           zero on downward leaders.  Otherwise can get weird stuff like
!           net negative charge to the storm by a negative CG!
!
! 8.6.2002  Increased CG threshold to 2500m (was 2000m).  
!           Also changed ibal=3 scheme to adjust propagation if imbalance
!           exceeds 5% (was 10%)
!
!
! 5.22.2002  New version (sparkmud.2.4b.f) for vertically stretched grid.
!            For the moment must have dslight = dzmax (i.e. lightning uses the
!            same vertical levels as on the dynamics grid and refinement allowed
!            in horizontal only)  There is no longer a choice about whether to set
!            phi=0 at kz=1 or kz = 0.  Now hard-wired for phi=0 at kz = 0 (i.e.
!            boundary is below ground).
!
!            Also fixed a bug in the code that puts lightning channel charge
!            back into the coarse (dynamics) grid.
!
! 4.8.2002  Fixed bug where still had mseed used for ichoose instead of iseed
!
! 1/17/2002  Now have more variables that are from inlightdbm,
!             including 
!                   ieint (selects constant or fractional internal field)
!                   eint1 (input value for constant internal field)
!                   ibal  (selects channel neutrality balancing method)
!                   fefac (individual propagation settings for pos/neg ends)
!              
!              Added code for the flag 'ieint' to set internal field values
!              
!              Made some changes to ibal=2 method:
!                If an initial potential shift causes imbalance to 'go the wrong way,'
!                then dont try to adjust it on the current step and instead reset 
!                ihold = 0 and let the channel grow a little more.  This seems to be
!                working well.
!                
!                 Also now check channel balance on every step up to 40 
!                 (instead of 4 previously) and after that only every 3rd step
!                 
!
! 1/7/2002  spark.2.4b.f   Try using a second-order difference to calculate e-field
!           next to channel when possible (otherwise revert to 1st order difference)
!
! 1/4/2002  spark.2.4a.f   Local arrays are now automatic or passed in from lightning
!            driver.  This allows changing dimensions without resetting the lightning
!            grid dimensions.
!
! 12/3/2001  Minor fix to reject small flashes (< 4 segments total)
!            when growth cant be forced.
!
! 11/26/2001 version sam.spark.2.4.f:
!            Combine regular and stand-alone versions into one.  Stand-alone
!            behavior is turned on by isa = 1
!
! 10/3/2001 2.3f.f
!    Fix up to allow arbitrary ds (i.e. anything other than the storm dz).

! 9/26/2001
!   increased eint to 5kV/m (from 500 V/m)
!
! 9/25/2001  spark.2.3e.f
!    Removed test for channel outside cloud and put in a difference in prop. threshold
!    for neg. vs. pos. leaders.  Start out with neg. leaders needing 20% higher field.
!
!  9/12/2001 sam.spark.2.3d.f:
!      Increased domain height to improve solution.
!  
!  9/10/2001  sam.spark.2.3c.f:
!        Correct the transfer of lightning charge to the coarse grid.   
!  9/8/2001 sam.spark.2.3b.f
!
!        Changed array t4 to have its logical meaning (i.e. t4 is positive if
!          t3 is positive)
!  8/22/2001
!       put in stronger smoothing (horizontal only) when the lightning grid dx is
!       close to the dynamics dx (dx1) (i.e. when dx1/dx .lt. 1.5 ).  If dx1/dx .gt. 1.5 
!       then the injection into the coarser grid has enough horizontal smoothing already.
!       (The reason for doing this is that a branch from one flash went up above the
!        cloud, and the charge it deposited up there, as 'ions' then instigated another
!        flash.  The second flash initiated _above_ the storm -- not very realistic!)
!
!  version sam.spark.2.3a.f:
!     Still not quite working, so grafted in old code for initial channel...
!     The problem turned out to be accessing the wrong array for ambient charge
!       when calculating the channel charge (should have used rho instead of rhow,
!       resulting in an offset of (nb,nb) gridpoints horizontally.  Oops.)
!
!  version sam.spark.2.3.f:
!     Returned to 2-point initial channel, which seemed to help.
!     
!  7/31/2001 version sam.spark.2.2.f
!            Now trying the 'proper' way to balance charge: by adjusting
!            the reference potential.
!
!            Start using a 3-point initial channel instead of 2-point.  The
!            positive and negative leaders now 'share' the middle point and
!            half of that charge.
!
!            Also fixed up 'signchan' to actually be the sign of charge
!            carried by that end.  (Was reversed before for arcane reasons)
!
!            Also fixed up accounting of the charge on the channel by
!            subtracting the ambient charge density.  This should actually
!            enhance slightly the charge on the channel.
!
!  7/31/2001 version sam.spark.2.1.f
!            Upon suggestion by D. MacGorman, try a new charge balancing 
!            procedure by which the lopsided end is kept from growing 
!            at all until the other end catches up.  (Works pretty well, too.)
!
!
!  Expanded version to feed potential back into SOR  
!  spark2r uses single precision arrays instead of double
!
! added in a call to the FACR routine in FISHPACK to compare with SOR
!   The FACR solution to the initial charge field is fed to the SOR as
!   a first guess.  Then SOR can work on the short wavelengths.
!
! 4/15/99  added diagonal propagation possibility
!
! 4/15/99  Now finally try it with a dipole distribution!
!
! 4/16/99  Now read in a space charge field from the cloud model
!          interpolates from a 1500x1500x500m  grid to 500x500x500m
!
! 4/18/99  Go down to 250m resolution for lightning propagation
!
! 4/21/99  Altered version of spark6b to simulate the Williams (1985)
!          experiments with space charge embedded in PMMA
! 4/25     change to allow for reversal of charge regions
!
! 5/1/99  version sparkw2 for simulating a rectangular slab of charge
!         with the discharge started at one end (Williams et al. (1985),
!         figure 1.)
!
! 5/2/99  Add calculation of charge density at each point in the
!         channel and just add up all the adjustments from the 
!         initial values.
!
!   5/5/99  rename from sparkw3 to sparkc1 and try again with a 
!           cloud model field
!
!  5/9/99  sparkc2
!          Try doing 'local' relaxation: just around the added bond
!
!  5/10/99 sparkc3: 500m grid
!
!  5/25/99 sparkc6: add break-even field option
!
!  6/15/99 sparkc7: testing speed-up on 500m grid
!
!  6/26/99 convert to a subroutine to be called from SAM lightning
!          subroutine
!
!  7/5/99  put a buffer zone extension to the lateral boundaries
!          so there is a charge-free zone next to the boundaries
!
!  8/8/99 allow propagation thresholds to adjust to promote growth 
!         on an end with substantially less charge and retard growth
!         of the end with more charge.  Idea is to check the channel 
!         charges every so many steps and alter the propagation  
!         threshold fraction fecrit for each end.  This should have
!         a similar effect as adjusting the channel potential.
!
!  8/11/99 increase minimum fecrit from 0.5 to 0.6
!  8/12/99 output flash bounding box, put in logical switch 'lsgi'
!          to easily switch between sgi and cray timing routines
!
!          If a CG occurs, don't record the channel charge of the 
!          ground-reaching end in the t3 array (but do calculate the
!          original charge on that end, anyway, and write to output
!          files)  Also, if fecritn of the upward-propagating end
!          back was greater than fecrit, then reset it to fecrit.  
!
! 8/22/99  When calculating probabilities, now check that the local
!          charge density of the bond point is > 1.0e-16.  If it is
!          not, then set Ecrit = Min ( 150kV/m, 2*Ecrit(z) ) to
!          simulate the increased breakdown field for clear air vs.
!          cloudy/hydrometeor-filled air.
!
! 9/12/99  Now after 3 flashes in the same time step we stop using
!          the offset and only start lightning where E exceeds 
!          inititation threshold 
!
!          Also don't relax after the initial starting point is
!          chosen (hope this will help avoid 1-segment channels)
!          (changed back again 9/17)
!
! 9/17/99  Now pick a new starting point if propagation stops with
!          4 or fewer total channel points (do this up to maxtry times-- it's 
!          a lot cheaper than doing the potential all over for a trivial
!          flash)  
!
!          Also reinstate the initial relaxation around the 
!          starting points.
!
!          Also treat the initial 2-point breakdown as a perfect 
!          conductor.  i.e. no internal field for the first points.
!
! 10/2/99  Change to set zero potential below ground level as in
!          sam.08.driver.sgi.f
!
! 10/4/99  set mask array values back to 1 when trying a new point after
!          a small flash
!
! 10/26/99 Stop IC development after 40 steps if charge imbalance
!          is greater than a facter of 2. (CGs are left alone)
!
! 2/12/2000 now stop flashes if level nz-2 is reached (used to be nz-1)
!           Also, changed
!              ecomp = Min( 150.0e3, 2.0*ecritz(k1) )
!           to
!              ecomp = Min( 150.0e3, Max(fecritn(iend),2.0)*ecritz(k1) )
!           just in case fecritn(iend) happens to be greater than 2.
!
! 3/15/2000 Changed lightning file naming convention.  Just reordered
!           the model time in front of the iteration number so that a normal 
!           file listing will put the files in chronological order.
!
! TO DO:
!
!        1. put in check that grid dimensions are OK
!
!#######################################################################
!
!  Variable declarations
!
!#######################################################################
      

!      include 'sam.files.h'

      integer nx,ny,nz,nb,nor  ! nb is number of boundary buffer zone points

      integer, parameter    :: ng1 = 1
      
      integer               :: iunit0,iunit
      integer               :: bcx, bcy
      integer               :: bcxtmp, bcytmp

      real dzz(nz),dxx(nx),dyy(ny)         ! dz(k),dx(i),dy(j)
!      real    :: gxt(nx,4), gyt(ny,4), gzt(nz,4)
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      real ex,ey,ez,emag,dxinv2,dyinv2,dzinv2
      real ex1,ey1,ez1,emag1, dxinv, dyinv, dzinv
      
      double precision tmp
      real tmp4
      integer ntry
      
!      real cgfr
      
      integer idown ! which end is going downward
      real phi1
      integer numic,numcgn,numcgp
      
      integer icnt1(2),icnt2(2),icnt3(2),icnt4(2)
      double precision chgdis(2),chgdis0(2)
      
!      real cgthres       ! altitude threshold to declare a CG
      real chanlen       !  channel length from initiation point to ground
      integer iextra     ! counter for max allowed extra segments after CG
      integer maxextra, ngroundstep
      integer imapz, mzdist, dzmax
      logical istretch  ! whether grid is stretched or not
      logical ldumpsp   ! flag to make special dump when something goes screwy
      logical lct       ! above-cloud energy test
      integer iixp,jjyq,kkzr,iiex,jjey,kkez
!      integer*8 llwork
      integer nnx,nny,nnz
      integer nxs,nys,nzs   ! indices of 'interior' of slm grid
      integer nbw,nbe,nbs,nbn,nbz2
      
      real eold,enew,eoldct,enewct
!
! note that nnzt = nnz-1
!
      integer nnzt
!      real phi(nnx,nny,nnz),rhs(nnx,nny,nnz)
!      real phi0(-nbw:nxs+nbe,-nbs:nys+nbn,0:nnzt)
!      real rhs0(-nbw:nxs+nbe,-nbs:nys+nbn,0:nnzt)
      integer NLx, NLy, NLz
! Note that NLx = nnx, NLy = nny, NLz = nnzt = nnz-1
      real phi(0:NLx+1,0:NLy+1,0:NLz+1) ! (NLx+2,NLy+2,NLz+2)
      real rhs(0:NLx+1,0:NLy+1,0:NLz+1) ! (NLx+2,NLy+2,NLz+2)
      real rhssor(0:NLx+1,0:NLy+1,0:NLz+1) ! (NLx+2,NLy+2,NLz+2)
      real phi0(-nbw:nxs+1+nbe,-nbs:nys+1+nbn,0:NLz+1)
      real rhs0(-nbw:nxs+1+nbe,-nbs:nys+1+nbn,0:NLz+1)
      real rhs0sor(-nbw:nxs+1+nbe,-nbs:nys+1+nbn,0:NLz+1)

!      real phi2(nnx,nny,nnzt)
      real phi2(NLx+2,NLy+2,NLz+2)
      
!      real potfair(nnz),rhofair(nnz)

!      real rho(nnx,nny,nnzt)
!      real rhow(-nbw:nxs+nbe,-nbs:nys+nbn,nnzt)

      real rho(0:NLx+1,0:NLy+1,0:NLz+1)
      real rhow(-nbw:nxs+1+nbe,-nbs:nys+1+nbn,0:NLz+1)

      integer nzmax
!      parameter (nzmax1 = 400)
!      real z1d2(nzmax1,4)
      real dlzz,dlz2
!      common/strz/z1d2,dlz,nz1d

      real eztop, pottop
      common/bndcz/eztop, pottop
      
      real,allocatable :: hgt(:),zmap(:)
      real,allocatable :: cofx(:,:),cofy(:,:),cofz(:,:)
      real,allocatable :: cz(:),czz(:)
      integer,allocatable :: footprint(:,:,:)

!      real hgt(nnz),zmap(nnz)
!      real cofx(nnx,3),cofy(nny,3),cofz(0:nnzt,3)
!      real cz(0:nnz-1),czz(0:nnz-1)

!
!     put integer and floating point arguments names in contiguous
!     storeage for labelling
!
      integer*8 iprm(22)
      integer mgopt(4)
      real fprm(8)


      integer*8 intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,      &
     &         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
      common/itmsp3/intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,      &
     &         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc0,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc0,yd,ze,zf,tolmax,relmax
      real dlx,dly ! ,dlz2,x,y,z,cxx,cyy,czz,cx,cy,cz,cex,cey,cez,ce,errm

      equivalence(intl,iprm)
      equivalence(xa,fprm)

      integer*8 iprm2d(16)
      real fprm2d(6)

      integer*8 intl2d,nxa2d,nxb12d,nyc2d,nyd2d,ixp2d,jyq2d,iex2d,jey2d,      &
     &              nxm2d,nym2d,      &
     &              iguess2d,maxcy2d,method2d,nwork2d,lwrkqd2d,itero2d
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
       common/itmsp2d/ intl2d,nxa2d,nxb12d,nyc2d,nyd2d,ixp2d,jyq2d,iex2d,jey2d,      &
     &              nxm2d,nym2d,      &
     &              iguess2d,maxcy2d,method2d,nwork2d,lwrkqd2d,itero2d
      real xa2d,xb2d,yc2d,yd2d,tolmax2d,relmax2d
      common/ftmsp2d/xa2d,xb2d,yc2d,yd2d,tolmax2d,relmax2d

      equivalence(intl2d,iprm2d)
      equivalence(xa2d,fprm2d)
!
!     declare coefficient and boundary condition input subroutines external
!
      external cfx,cfy,cfz,cfzlgt,cfzu,bndc,bndc2d


      double precision e,e1     ! temp value of efield next to channel
!      integer l2nde  !  =2 to use 2nd order difference for e
                     !  =1 for 1st order
!      parameter (l2nde = 1)
      integer icheckdims
      integer icldtop
      real    zcldtop
      integer ihold
      real fx,fz
      real scintmx, scintmn 
      real phidiff
      real bfac

!      integer ibal  ! flag for which balancing procedure to use

      double precision qtot1,qtot2
      integer iadd,iaddz,istop1
      double precision chg,chg0
      integer ilast, iflip
      real adjlast
      
      integer iw, ie, js, jn, ku, kd
      real    xw, xe, yn, zu, zd
      real   xc, yc, zc
      real frs, frn, fre, frw, fru, frd


      
! input stuff
      integer pchnl,nchnl  ! number of segments in neg/pos ends
      integer icgyn  ! flag for cg flash
      integer icgpn  ! flag for positive(1)/negative(-1) cg 
      real dx1,dy1,dz1  !  
      integer iliter  !  flash counter for the timestep
      integer lgtstp  !  flag to exit lightning
!      real elgtfdel  !  fraction of E that must be exceeded for possible
                     !  flash initiation point
      real dtp      ! time step size
      integer nstep1
      integer ismall ! flag for small flash and to try another


                 
! nsavephi:  =2 or 3 read facr calculated phi from file,1=write facr phi to file
!   =3 write phi after sor, =4 read sor-relaxed phi (skips first relaxation)
!  nsavephi =1 1 causes call to fishpack
      integer ndebug
      parameter (ndebug=DEBUG)
      integer interphi
      parameter (interphi = 0)  ! 0 = interpolate rho, solve phi
                                ! 1 = solve phi, interpolate phi
      integer ireverse
!      parameter (ireverse = 0)    ! =1 to swap density regions
      integer ntest
      parameter (ntest=0)  ! =2 stops code somewhere....
      integer idphi
      parameter (idphi=0)  ! =1 to dump phi from fishpack and quit

      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)


!
!  vars for stand-alone call
!
      integer ixst,jyst,kzst
      character(len=80) outname
      integer nkzmn(nzs)

      integer nxb,nyb,nzbt
      integer jslab
! use dx = dz = ds 'cause it's easier that way       
      real ds
      double precision ds2,ds3,dls3
      real eps0,rjac,pi,omega,soreps0
!      real  mask(nnx,nny,nnzt)
      real mask(0:NLx+1,0:NLy+1,0:NLz+1)
      real phip(0:NLx+1,0:NLy+1,0:NLz+1)
      integer pickmask(0:NLx+1,0:NLy+1,0:NLz+1)
      logical :: pickmaskcheck
!      integer :: pickmaskwidth = 2
!      integer :: maxpicks = 2 ! 10
!      integer :: multipickmod = 20 ! add an extra pick for each number of added steps
      integer :: numpick, pickloop, addpick, badpick
      real maskmx,maskmn
!      real,allocatable ::  mask(:,:,:)
!      real,allocatable ::  phip(:,:,:)   !   (nnx,nny,0:nnzt)
      integer i,j,k,i1,j1,k1,n,ix,jy,kz,nstop,nstep,time,m
      real    time_real
      integer iter,maxiter
      integer idif, jdif, kdif
      integer i2,j2,k2,i3,j3,k3
      real r2
      integer ibc  !  value for SOR bound. conds: 1=Neum. 0=Diri.
      integer ibcn,ibcs,ibce,ibcw
      integer itype
      real resid
      real temp
      integer, parameter :: nend=2  ! 1 or 2 for uni or bidirectional
      integer :: numextrapicks(nend)
      real    :: chanlentot(nend)
      real    :: chanareatot(nend+1)
      integer nxyz
      integer, parameter :: nxyz2 = 30000/nend
      integer chanijk(nxyz2,3,nend) ! i,j,k of channel points
      integer ibonds(nxyz2,3,2,nend) ! list of connected points
      integer linksrc(nxyz2,nend) ! list of source points, i.e., the bond number from which this one branched.
      integer linkstep(nxyz2,nend) ! list of source points, i.e., the bond number from which this one branched.
      integer linkrank(nxyz2,nend) ! list of tile ranks, i.e., the rank of the tile with this bond
      integer linkoffset(nxyz2,2,nend) ! offset to get global location of point
      integer ioffset, joffset, koffset
      integer ishift, jshift, kshift
      integer NLxchoose, NLychoose, NLzchoose
      real    phichan(nxyz2,nend)   ! potential of channel point
      real signchan(2)           ! -1.0 or +1.0 depending on direction
      integer numchan(2)         ! number of points in channel
      integer numchancheck(2)         ! number of points in channel
      integer :: numchanmin  ! min of the number of channel points on each end at the beginning of each nstep (lightloop) integration
      integer npbond(nxyz2,nend)    ! number of potential bonds from point nxyz
      real bondinfoin(5,nxyz2) ! info at each point for MPI comms
      real bondinfoout(5,nxyz2) ! info at each point for MPI comms

!
! with diagonals, have 9+8+8=25 possible unbonded neighbors (at least
!  one neighbor is already bonded).  So change 5 to 25
!
      integer bondijk(25,nxyz2,3,nend) ! locations of pot. bond points
!      real  bondphi(25,nxyz2,2)     ! phi of potential bond
      real  bondprob(25,nxyz2,nend)    ! probability of bond
      double precision  probtot, probcum
      double precision  probtotall, probcumall
      real  phistart,phistlow,phisthi,phistart0
      real  phistart1, phistart2
      real  chandist(nxyz2,nend)
      real  cgchandist(nxyz2,nend)
      integer istop(2)  ! flag for no propagation on one end or the other
      integer istopcnt(2)
!
!  need arrays for distance between neighboring points
!
      real,allocatable,save :: gd(:,:,:,:,:)
!      real gd(-1:1,-1:1,-1:1,0:nnzt,2)  !  1: distance, 2: inverse of distance
      real :: gdstr(-1:1,-1:1,-1:1,2)
      real,allocatable,save :: gdxy(:,:,:),gdxz(:,:,:), gdyz(:,:,:),gdxyz(:,:,:,:) ! factors for horizontal stretch
      real,allocatable,save :: gdx(:,:), gdy(:,:)
      
      real ran0mpi   ! random number function
      double precision pick
      integer npick,bpick
      integer eta    !  power law for probability
      parameter (eta=1)
!      integer ieint   ! =1 then eint(k) = eint
                      ! =2 then eint(k) = feint*ecritz(k)
!      real  eint     ! value of constant internal field
      real eintz(nzs)     ! internal field of channel
      integer iend,iend0
      integer iendgnd ! says which end hit ground
!      real  feint
      integer mstop   ! input value for number of iterations on bond-adding
                      ! before stopping. If mstop < 0, then use max value. 
      integer intype  ! whether input charge field is (=1) single array or
                      !   (=anything else) a full mdirt file
      integer icheck  ! set to 1 to activate stopping 2nd direction of 
                      ! leader after propagation stops in the other end.
      integer nstop2
      
!      real tem1(nnx,nny,nnzt)
       real emax,emax2,emagmax,ezmax,ehmax
       real dx,dy,dz,dzi
       real x,y,z

      real dzl,dyl,dxl
      real tdx,tdy,tdz
!      integer neelec
      real pot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
!      integer ipot
!      parameter (ipot=5)

      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

!      character*40 fname   ! input filename with charge density field
!      character*40 outfile  ! name for binary output
!      character*40 tfile    ! name for text file with link points
      character*6  rstime
      character*2  nmliter
!      integer idump   ! input flag for binary output
      real eps1
      parameter (eps1 = 1.e-2)
      real facx,facy,facz
      real scint,potint
      integer ifile,istat
!
!  electrical permitivity of air C / (N m**2) -  check the units
!
      real eperao
      parameter (eperao  = 8.8592e-12 )
!
!  arrays used in FISHPAK for solver for electric potential field
!
      integer id1,jd1,kd1
      parameter (id1=0,jd1=0,kd1=0)

      integer nwslv

      integer nzb

! other stuff for poisson solver
      integer ierror
      real xs, xf, ys, yf, zs ! , zf  ! domain limits in meters
      real pertrb,elmbda

      integer ist,jst,kst,ist1  ! endpoint of needle from groundplane
      integer jend
      integer isp,jsp
!
! stuff for charge removed
!

      double precision qtot(2)
      double precision rhotemp, rhox, phix, energy
      double precision rhotot(2),rhotot1(2)
      logical ltest
      real arhopos, arhoneg, rhopos1, rhopos2, rhoneg1, rhoneg2  
      
! vars for localsor
      integer i0,j0,k0
      integer iloc0(nend),jloc0(nend),kloc0(nend)
      integer  mdelf
      real omegal
      integer imn,jmn,kmn,imx,jmx,kmx
      integer imn1,jmn1,kmn1,imx1,jmx1,kmx1
      integer imn2,jmn2,kmn2,imx2,jmx2,kmx2
      
!
! vars for allowing for breakeven field threshold
!
!      integer ibrkd  ! = 1 for constant breakdown threshold
                     ! = 2 for breakeven field
      real ebrkd(nzs)  !  this array is passed in -- don't change it!
      real ebrkds(nzs)
      real ecritz(nzs)
      real zlev(nzs)
      
       integer ninit  ! number of possible initiation points
       integer npick1
       integer,allocatable :: ipick(:,:)
!       integer ipick(nxs*nys*nzs/10,3)
       integer ichoose
       integer idum, mseed
       data mseed /678/
       real etest
!       real fecrit
       real fecritn(nend)  ! 
!       real fefac(nend)
       real fecinc   ! increment for adjusting fecrit
       real fecmin   ! minimum value for fecrit
!
! timing vars
!
      real tim1,tim2,timarr(2)
      real second,commas_dtime
      logical lsgi
      parameter( lsgi = .true. )
      
      real times
      double precision volbi,vol,volratio
!
! Misc. vars
!
      integer istart
      data istart /0/
      save istart
      integer pstop
      real ecomp
      integer numthr
      integer mp_numthreads
      
      integer numtry  ! counts number of tries to get a non-trivial flash
      integer maxtry  ! maximum number of tries to get a non-trivial flash
      parameter (maxtry = 6)
      
      integer ifac
      real    delta
      
      double precision :: dt1, dt2, delt, dtSTRZSORmpi
      integer nSTRZSORmpi, nboxmgiter,numitstot,numits
      integer ixe, jye, kze

      integer,allocatable,save :: init_send_buffer(:)
      integer,allocatable,save :: init_recv_buffer(:,:)
      integer,allocatable,save :: mpi_init0(:,:)

      double precision,allocatable,save :: prob_send_buffer(:)
      double precision,allocatable,save :: prob_recv_buffer(:,:)
      double precision,allocatable,save :: mpi_prob0(:,:)
      
      integer,parameter :: ntot = 50
      integer :: ichoose_send_buffer(ntot)
      integer :: ichoose_recv_buffer(ntot)
      real    :: rchoose_send_buffer(ntot)
      
      integer :: mpitotinint(ntot), mpitotoutint(ntot)
      real    :: mpitotin(ntot), mpitotout(ntot)
      double precision :: mpitotindp(ntot), mpitotoutdp(ntot)
      
      integer :: init_send_size, init_recv_size, nz1, ninittmp, ninittot
      integer :: prob_send_size, prob_recv_size
      
      integer :: ichooserank, ichoosenumber
      integer :: ichoosenextrank
      
!      logical, parameter :: usebox = .true.
      logical, parameter :: usebox = .false.
      
      logical edgecheck 

            
! #################################################      
! Begin Executable Code      
! #################################################      


      IF ( ndebug > 0 ) write(0,*) 'Start of spark. my_rank = ',my_rank
      soreps0 = soreps
      
      IF ( ndebug .gt. 0 ) THEN
         iunit = 0
      ELSE
         iunit = iunit0
      ENDIF
      
      IF ( lightintz .gt. 1 ) igrid = 1
      
#ifdef MPI
        l2nde = 1  ! turn off 2nd order for now. 
                   ! Can allow later with some logic to revert to 1st order at tile boundaries
        igrid = 1
#endif
      
      allocate ( hgt(nnz) )
      allocate ( zmap(nnz) )
      allocate ( cofx(nnx,3) )
      allocate ( cofy(nny,3) )
      allocate ( cofz(0:nnzt,3) )
      allocate ( cz(0:nnz-1) )
      allocate ( czz(0:nnz-1) )
!      allocate ( footprint(-nbw:nxs+1+nbe,-nbs:nys+1+nbn,nend) )
      allocate ( footprint( 0:NLx+1,0:NLy+1, nend) )

       npick1 = number_of_processes*nxs*nys*nzs/Min( 5, 11-lightintz )

      allocate ( ipick(npick1,3) )

       IF ( ibal .eq. 3 .or. ibal .ge. 4) THEN
          fecmin = fecrit/2.0
       ELSE
          fecmin = 0.5
       ENDIF
       
!       allocate ( mask(-1:nnx+2,-1:nny+2,nnzt) )
!       allocate ( phip(-1:nnx+2,-1:nny+2,0:nnzt+1) )
       

!
! value of ds is now passed in!
!
      ds  = dslight
      maxextra = Int(500*(Max(1.,500./ds)))/2
      ngroundstep = maxextra
      ds2 = ds**2
      dx = ds
      dy = ds
      dxl = dx1
      dyl = dy1
      dzl = dz1
      IF ( istretch ) THEN
        dz = dlzlgt
        dls3 = dlzlgt*ds**2
        ds3 = ds2*dzl
      ELSE
        dz = dslightz
        dls3 = dz*ds**2
        ds3 = ds2*dz
      ENDIF
      
      numextrapicks(:) = 0

      dzi = 1.0/dz
      dxinv2 = 0.5/dx
      dyinv2 = 0.5/dy
      dzinv2 = 0.5/dz

      dxinv = 1./dx
      dyinv = 1./dy
      dzinv = 1./dz
      
      volratio = 1.d0/(lightintx**2*lightintz)
      
      DO kz = 1,nnz
        IF ( istretch ) THEN
          hgt(kz) = z1d2lgt(kz,1) 
          zmap(kz) = 1.0/z1d2lgt(kz,3) ! vertical map factor ( dz_ave/dz(kz) )
        ELSE
          hgt(kz) = 0.5*dz*(dz1/dz) + (kz - 1)*dz
          zmap(kz) = 1.0
        ENDIF
      ENDDO
      
       DO i=1,nnx
         cofx(i,1) = 1.0/x1d2lgt(i)**2 ! 1.0/ds2
         cofx(i,2) = 1.0/x1d2lgt(i)**2 ! 1.0/ds2
         cofx(i,3) = -(cofx(i,1)+cofx(i,2))
       ENDDO

       DO j=1,nny
         cofy(j,1) = 1.0/y1d2lgt(j)**2 ! 1.0/ds2
         cofy(j,2) = 1.0/y1d2lgt(j)**2 ! 1.0/ds2
         cofy(j,3) = -(cofy(j,1)+cofy(j,2))
       ENDDO
       
       dlzz = dlzlgt*dlzlgt
       dlz2 = dlzlgt+dlzlgt
       cz(0) = z1d2lgt(1,3)*(z1d2lgt(2,3)-z1d2lgt(1,3))/dlzlgt
       czz(0) = z1d2lgt(1,3)**(2)
       cz(1) = cz(0)
       czz(1) = czz(0)
       
        DO k = 0,1
           cofz(k,1) = czz(k)/dlzz-cz(k)/dlz2
           cofz(k,2) = czz(k)/dlzz+cz(k)/dlz2
           cofz(k,3) = -(cofz(k,1)+cofz(k,2))
        ENDDO

       IF ( istretch ) THEN
       DO k=2,nnz-1
         IF ( k .lt. nnz-1) THEN
          cz(k) = z1d2lgt(k,3)*(z1d2lgt(k+1,3)-z1d2lgt(k-1,3))/dlzlgt*0.5
         ELSE
          cz(k) = z1d2lgt(nnz-1,3)*(z1d2lgt(nnz-1,3)-z1d2lgt(nnz-2,3))/dlzlgt
         ENDIF
         czz(k) = z1d2lgt(k,3)**(2)
         cofz(k,1) = czz(k)/dlzz-cz(k)/dlz2
         cofz(k,2) = czz(k)/dlzz+cz(k)/dlz2
         cofz(k,3) = -(cofz(k,1)+cofz(k,2))
         
       ENDDO
       
       ELSE  ! unstretched grid

       DO k=0,nnzt
         cofz(k,1) = 1.0/ds2
         cofz(k,2) = 1.0/ds2
         cofz(k,3) = -(cofz(k,1)+cofz(k,2))
       ENDDO
       
       
       ENDIF ! istretch

  45  CONTINUE

       icgyn = 0
       icgpn = 0
       iendgnd = 0
       fecinc = 0.05
       mdelf = Max(mdel,15)
       imn = nxs-1
       jmn = Max(1, nys-1)
       kmn = nzs-1
       imx = 2
       jmx = 2
       kmx = 2
       imn1 = nxs-1
       jmn1 = Max(1, nys-1)
       kmn1 = nzs-1
       imx1 = 2
       jmx1 = 2
       kmx1 = 2
       ibc = 0  !  1 for Neumann condition on top and sides
       ibcn = 0
       ibcs = 0
       ibce = 0
       ibcw = 0
       IF ( ixbeg == nxbeg ) ibcw = 1
       IF ( ixend == nxend ) ibce = 1
       IF ( jybeg == nybeg ) ibcs = 1
       IF ( jyend == nyend ) ibcn = 1

       IF ( ny .eq. 2 ) THEN
        jmx = 1
        jmx1 = 1
       ENDIF

      probtot = 0.0
      DO n=1,2
      numchan(n) = 0
      signchan(n) = 0.0
      END DO
      nstop = nxyz2
      phistart = 0.0
      
      pi = 4.00*ATan(1.0d0)
      r2 = (3.0/(4.0*pi))**(2.0/3.0)
      rjac = 0.86
!      rjac = 0.745356
      omega = 2.00/(1.00 + Sqrt(1.00 - rjac**2))
#if DEBUG > 0
      IF (ndebug .ge. 2) write(iunit,*) 'pi,rjac,r2',pi,rjac,r2,omega
#endif
      

       mstop = -1
       icheck = 1

#if DEBUG > 0
      IF (ndebug .ge. 1) THEN
        write(iunit,*) 'eta = ',ltgeta
        write(iunit,*) 'intfgmpi = ',intfgmpi
        write(iunit,*) 'ds = ',ds
      END IF
#endif
      
      IF (mstop .ge. 0) THEN
        write(iunit,'(a)') 'mstop .ge. 0, so setting nstop'
        nstop = Min(mstop,nstop)
      ENDIF   !  
      IF (icheck .ne. 0) icheck = 1

      

       fecritn(1) = fecrit
       fecritn(2) = fecrit

!
!      IF (ibrkd .le. 2) THEN
!       DO kz = 1,nzs
!
!c      eintz(kz) = feint*ebrkds(kz)
!c       eintz(kz) = 5000.0
!c      ecritz(kz) = ebrkds(kz)
!       fecritn(1) = fecrit
!       fecritn(2) = fecrit
!       ENDDO
!      ELSE
!        DO kz=1,nzs
!c        ebrkds(kz) = 150.0e3
!c        fecrit = ecrit/ebrkds(1)
!        fecritn(1) = fecrit
!        fecritn(2) = fecrit
!c        ecritz(kz) = ebrkds(kz)
!c        eintz(kz) = feint*ecrit
!        END DO
!      END IF  !  (ibrkd .le. 2)


      IF ( .not. allocated( gd ) ) THEN
      allocate ( gd(-1:1,-1:1,-1:1,0:nnzt,2) )
      allocate ( gdx(0:nnx,2) )
      allocate ( gdy(0:nny,2) )
      allocate ( gdxy(0:nnx,0:nny,2) )
      allocate ( gdxz(0:nnx,0:nnzt,2) )
      allocate ( gdyz(0:nny,0:nnzt,2) )
      allocate ( gdxyz(0:nnx,0:nny,0:nnzt,2) )
      
      gdx(:,:) = 1.0
      gdy(:,:) = 1.0
      gdxy(:,:,:) = 1.0
      gdxz(:,:,:) = 1.0
      gdyz(:,:,:) = 1.0
      gdxyz(:,:,:,:) = 1.0
      

      IF ( istretch ) THEN
      
      DO kz=0,nnzt
      DO k=-1,1
      DO j=-1,1
      DO i=-1,1
        IF (i.eq.0 .and. j.eq.0 .and. k.eq.0) THEN
         gd(i,j,k,kz,1) = 0.0
         gd(i,j,k,kz,2) = 0.0
         cycle
        ENDIF

        gd(i,j,k,kz,1) =       &
     &     Sqrt(dx**2*Abs(i) + dy**2*Abs(j) +       &
     &           (dz/z1d2lgt(Max(1,kz+Max(0,k)),3))**2*Abs(k))

          gd(i,j,k,kz,2) = 1.0/gd(i,j,k,kz,1)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      ELSE
      
      DO kz=0,nnzt
      DO k=-1,1
      DO j=-1,1
      DO i=-1,1
        IF (i.eq.0 .and. j.eq.0 .and. k.eq.0) THEN
         gd(i,j,k,kz,1) = 0.0
         gd(i,j,k,kz,2) = 0.0
         cycle
        ENDIF

        gd(i,j,k,kz,1) =       &
     &     Sqrt(dx**2*Abs(i) + dy**2*Abs(j) +  dz**2*Abs(k) )

          gd(i,j,k,kz,2) = 1.0/gd(i,j,k,kz,1)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF ! istretch

      IF ( hstretch ) THEN
        ! stretch factors in 1,0,0 direction
        write(iunit,*) 'hstretch,dx = ',hstretch,dx
        DO i = 0,nnx
          gdx(i,1) = Max(dx,x1d2lgt(i))/dx
          gdx(i,2) = 1.0/gdx(i,1)
          write(iunit,*) 'i,gdx = ',i,gdx(i,1),gdx(i,2)
        ENDDO
        ! stretch factors in 0,1,0 direction
        DO j = 0,nny
          gdy(j,1) = Max(dy,y1d2lgt(j) )/dy
          gdy(j,2) = 1.0/gdy(j,1)
!          write(iunit,*) 'j,gdy = ',j,gdy(j,1)
        ENDDO

        ! diagonals in XY plane: 1,1,0
        DO j = 1,nny
          DO i = 0,nnx
            
            gdxy(i,j,1) = Sqrt((dx*gdx(i,1))**2 + (dy*gdy(j,1))**2) / &
              Sqrt(dx**2 + dy**2 )
            IF ( gdxy(i,j,1) /= 0.0 ) THEN
              gdxy(i,j,2) = 1.0/gdxy(i,j,1)
            ELSE
              gdxy(i,j,2) = 1.
            ENDIF
            
          ENDDO
        ENDDO
      
        ! diagonals in XZ plane: 1,0,1
        DO k = 1,nnzt
          DO i = 0,nnx
            
            gdxz(i,k,1) = Sqrt((dx*gdx(i,1))**2 + (dz/z1d2lgt(k,3))**2) / &
              Sqrt(dx**2 + (dz/z1d2lgt(k,3))**2 )
            IF ( gdxz(i,k,1) /= 0.0 ) THEN
              gdxz(i,k,2) = 1.0/gdxz(i,k,1)
            ELSE
              gdxz(i,k,2) = 1.
            ENDIF
            
            IF ( k == nnzt/3 ) THEN
!              write(iunit,*) 'i,gdxz = ',i,gdxz(i,k,1),gdx(i,1),(dz/z1d2lgt(k,3)),Sqrt(gdx(i,1)**2 + (dz/z1d2lgt(k,3))**2), &
!              Sqrt(dx**2 + (dz/z1d2lgt(k,3))**2 )
            ENDIF
            
          ENDDO
        ENDDO

        ! diagonals in YZ plane: 0,1,1
        DO k = 1,nnzt
          DO j = 0,nny
            
            gdyz(j,k,1) = Sqrt((dy*gdy(j,1))**2 + (dz/z1d2lgt(k,3))**2) / &
              Sqrt(dy**2 +  (dz/z1d2lgt(k,3))**2 )
            IF ( gdyz(j,k,1) /= 0.0 ) THEN
              gdyz(j,k,2) = 1.0/gdyz(j,k,1)
            ELSE
              gdyz(j,k,2) = 1.
            ENDIF

            IF ( k == nnzt/3 ) THEN
!              write(iunit,*) 'j,gdyz = ',j,gdyz(j,k,1)
            ENDIF
            
          ENDDO
        ENDDO

        ! long diagonals: 1,1,1
        DO k = 1,nnzt
          DO j = 0,nny
          DO i = 0,nnx
            
            gdxyz(i,j,k,1) = Sqrt((dx*gdx(i,1))**2 + (dy*gdy(j,1))**2 + (dz/z1d2lgt(k,3))**2) / &
              Sqrt(dx**2 + dy**2 +  (dz/z1d2lgt(k,3))**2 )
            IF ( gdxyz(i,j,k,1) /= 0.0 ) THEN
              gdxyz(i,j,k,2) = 1.0/gdxyz(i,j,k,1)
            ELSE
              gdxyz(i,j,k,2) = 1.
            ENDIF

            IF ( k == nnzt/3 .and. i == j ) THEN
!              write(iunit,*) 'i,gdxyz = ',i,gdxyz(i,j,k,1)
            ENDIF
          
          ENDDO
          ENDDO
        ENDDO

      ENDIF
      
      ENDIF ! allocated
      
      DO m=1,nend
      DO n=1,nxyz2
      chandist(n,m) = 0.0
      chanijk(n,1,m) = 0
      chanijk(n,2,m) = 0
      chanijk(n,3,m) = 0
      npbond(n,m) = 0
      phichan(n,m) = 0.0
      ibonds(n,1,1,m) = 0
      ibonds(n,2,1,m) = 0
      ibonds(n,3,1,m) = 0
      ibonds(n,1,2,m) = 0
      ibonds(n,2,2,m) = 0
      ibonds(n,3,2,m) = 0
       DO i=1,25
        bondijk(i,n,1,m) = 0
        bondijk(i,n,2,m) = 0
        bondijk(i,n,3,m) = 0
        bondprob(i,n,m) = 0.0   
       END DO ! i
      END DO ! n
      END DO ! m

      numextrapicks(:) = 0
      
! set up charge distribution

!  45  CONTINUE
! IF THIS IS NOT THE FIRST TIME (stand-alone), THEN 
!   SKIP OVER THE CHARGE/POTENTIAL SETUP, but still reset the mask and tem1

      
      IF ( nliter .gt. 1 .and. isa .ge. 1 ) THEN ! {
        iliter = iliter + 1
        IF ( iliter .gt. nliter ) THEN ! {
       
         IF ( istop(1) .eq. 2 ) THEN
           numcgn = numcgn + 1
         ELSEIF ( istop(2) .eq. 2 ) THEN
           numcgp = numcgp + 1
         ELSE
           numic = numic + 1
         ENDIF
          GOTO 999
        ENDIF ! }

      IF ( my_rank == 0 ) THEN
       write(iunit,*) '*********************************************'
       write(iunit,*) '****** iliter = ',iliter,'*************'
       write(iunit,*) '*********************************************'
      ENDIF
      ENDIF ! }

      IF ( iliter .gt. 2 .and. isa .ge. 1 ) THEN
!           tem1(i,j,k) = 0.0
! (0:NLxlg+1,0:NLylg+1,0:NLzlg+1)
        mask(:,:,:) = 1.0e30 ! set a very large number (larger than any valid potential value) to represent "no channel". 
                             ! Then when there is a channel, the value is set to an initial value of potential
                             ! This allows us to use the BoxMG subroutines for a potential solve
        pickmask(:,:,:) = 0
!       DO k = 1,nnzt
!         DO j = -1,nny+2
!          DO i = -1,nnx+2
!           mask(i,j,k) = 1.0 
!          END DO
!         END DO
!        END DO
       
       IF ( istop(1) .eq. 2 ) THEN
         numcgn = numcgn + 1
       ELSEIF ( istop(2) .eq. 2 ) THEN
         numcgp = numcgp + 1
       ELSE
         numic = numic + 1
       ENDIF
       
       
       GOTO 55

      ENDIF ! ( iliter .gt. 2 .and. isa .ge. 1 )

! set up charge distribution
#if DEBUG > 0
      IF (ndebug .ge. 1) write(iunit,*) 'setup charge: nxs,nys,nzs=',nxs,nys,nzs
#endif
! C$DOACROSS LOCAL(i,j,k)
! C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)

        mask(:,:,:) = 1.0e30
        pickmask(:,:,:) = 0
!        DO k = 1,nnzt
!         DO j = -1,nny+2
!          DO i = -1,nnx+2
!           mask(i,j,k) = 1.0 
!          END DO
!         END DO
!        END DO

!      DO k = 2,nnzt
!      DO j = 1,nny
!      DO i = 1,nnx
!        rho(i,j,k) = rhofair(k) ! 0.0
!        rhs(i,j,k) = rhofair(k-1)/eperao ! 0.0
!        phi(i,j,k) = 0.0 ! potfair(k-1)
!      END DO
!      END DO
!      END DO


! (0:NLx+1,0:NLy+1,0:NLz+1)

      DO k = 1,nnzt-1
        rho(1:NLx,1:NLy,k) = rhofairlgt(k)
        DO j = 1,Nly
        DO i = 1,Nlx
        rhs(1:NLx,1:NLy,k) = rhofairlgt(k)/eperao*x1d2lgt(i)*y1d2lgt(j)*(dlzlgt/z1d2lgt(k,3)) ! volume multiplied for boxmg solver
        ENDDO
        ENDDO
        rhssor(1:NLx,1:NLy,k) = -rhofairlgt(k)/eperao
      ENDDO
        phi(:,:,:) = 0.0

      DO k = 1,1
      DO j = 1,nny
      DO i = 1,nnx
        rho(i,j,1) = rhofairlgt(1) ! 0.0
        rhs(i,j,1) =  0.0
        phi(i,j,1) = 0.0

        phi(i,j,nnz-1:nnz) = pottop ! potfairlgt(nztop)
        rhs(i,j,nnz) = rhofairlgt(nnz-1)/eperao*x1d2lgt(i)*y1d2lgt(j)*(dlzlgt/z1d2lgt(k,3))
        rhssor(i,j,nnz) = -rhofairlgt(nnz-1)/eperao
      END DO
      END DO
      END DO

       
      

!     subroutine mlint2
!    > (ndebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,
!    >  nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
!
      ixe = nxs
      IF ( myproci == nproci ) ixe = nxs - 1
      
      jye = nys
      IF ( myprocj == nprocj ) jye = Max(1, nys - 1)

      scintmx = 0.0
      scintmn = 0.0
      ifac = lightintx ! Nint(dxl/ds)
!      IF ( ifac /= lightintx ) THEN
!        write(0,*) 'ifac,lightinx = ',ifac,lightintx,dxl,ds
!      ENDIF
!      jend = Max(1, nys-1)
      DO k=1,nzs-1
        IF ( istretch .or. dz .eq. dz1 ) THEN
         IF ( lightintz .eq. 1 ) THEN
          k1 = k
          facz = 0.0
         ELSEIF  ( lightintz .eq. 2 ) THEN
          k1 = (k + 1)/lightintz
          facz = 0.0 + 0.5*(k + 1 - k1*lightintz) ! facz is 0 or 0.5
         ENDIF
        ELSE
          k1 = Int( ((k-1)*dz + eps1)/dzl ) + 1
          facz = ( (k-1)*dz - (k1-1)*dzl )/dzl
        ENDIF
      DO j=1,jye
        j1 = Int( ((j-1)*dy + eps1)/dyl ) + 1
        facy = ( (j-1)*dy - (j1-1)*dyl )/dyl
        IF ( ny .eq. 2 ) THEN ; facy = 0.0 ; j1 = 1 ; ENDIF
      DO i=1,ixe
        i1 = Int( ((i-1)*dx + eps1)/dxl ) + 1
        facx = ( (i-1)*dx - (i1-1)*dxl )/dxl
       IF ( igrid .eq. 0 ) THEN
       IF ( i1 .ge. 1   .and. j1 .ge. 1   .and. k1 .ge.  1 .and.      &
     &     i1 .le. nx-1 .and. j1 .le. ny-1 .and. k1 .le. nz-1 ) THEN
        call mlint2      &
     & (0, t0,  -nor+1, nx+nor, -nor+1, ny+nor, -nor+1, nz+nor,       &
     &  1, 1, nx, ny, nz, scint, facx, facy, facz,i1,j1,k1, 1)
        call mlint2      &
     & (0, pot,   -nor+1, nx+nor, -nor+1, ny+nor, -nor+1, nz+nor,      &
     &  1, 1, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,1)
!       scint = t0(i1,j1,k1)
!       potint = pot(i1,j1,k1)
       rhow(i,j,k) = scint
       phi0(i,j,k) = potint
       scintmx = Max(scintmx,scint)
       scintmn = Min(scintmn,scint)
       ENDIF
       ELSEIF ( igrid .eq. 1 ) THEN
       IF ( i1 .ge. 1  .and. j1 .ge. 1   .and. k1 .ge.  1 .and.      &
     &     i1 .le. ixe .and. j1 .le. jye .and. k1 .le. nz-1 ) THEN
!        call mlint2
!     > (0, t0,  -nor, nx+nor, -nor, ny+nor, -nor, nz+nor, 
!     >  1, 1, nx, ny, nz, scint, facx, facy, facz,i1,j1,k1, 1)
!        call mlint2
!     > (0, elec,  -nor, nx+nor, -nor, ny+nor, -nor, nz+nor, 
!     >  1, neelec, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,ipot)
       scint = t0(i1,j1,k1)
       potint = pot(i1,j1,k1)
       rhow(i,j,k) = scint
       phi0(i,j,k) = potint
       scintmx = Max(scintmx,scint)
       scintmn = Min(scintmn,scint)
       ENDIF
       ENDIF
      END DO ! i
      END DO ! j
        
        ecritz(k) = ebrkds(k1) + facz*(ebrkds(k1+1) - ebrkds(k1))
        IF ( ieint .eq. 1 ) THEN
          eintz(k) = eint
        ELSE
          eintz(k) = feint*ecritz(k)
        ENDIF
        
      END DO ! k
      
!      IF ( isa .lt. 1 ) THEN
#ifdef MPI
         CALL MPI_AllReduce(scintmx, tmp4, 1, MPI_real, MPI_Max, my_comm, mpi_error_code)
         scintmx = tmp4
         CALL MPI_AllReduce(scintmn, tmp4, 1, MPI_real, MPI_Min, my_comm, mpi_error_code)
         scintmn = tmp4
#endif
        IF ( my_rank == 0 ) THEN
         write(iunit,*) 'scintmx,scintmn = ', scintmx, scintmn 
        ENDIF
!      ENDIF
            
      ixe = nxs
      IF ( myproci == nproci ) ixe = nxs - 1
      
      jye = nys
      IF ( myprocj == nprocj ) jye = Max(1, nys - 1)


#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'setup forcing'
#endif
! C$DOACROSS LOCAL(i,j,k)
!      jend = Max(1, nys-1)
#ifndef MPI
!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
#endif
      DO k = 1,nzs-1
      DO j = 1,jye
      DO i = 1,ixe 
        rhs0(i,j,k)  = rhow(i,j,k)/eperao*x1d2lgt(i+nbw)*y1d2lgt(j+nbs)*(dlzlgt/z1d2lgt(k,3))
        rhs0sor(i,j,k)  = -rhow(i,j,k)/eperao
      END DO
      END DO
      END DO
      
! turn this off for MPI (2015/02/23)
!      IF ( bcx .eq. 2 ) rhs0(nxs,1:jend,1:nzs-1) = 
!     :                          rhs0(1,1:jend,1:nzs-1)
!      IF ( bcy .eq. 2 ) rhs0(1:nxs-1,nys,1:nzs-1) = 
!     :                          rhs0(1:nxs-1,1,1:nzs-1)

      
!      tim1 = commas_dtime(timarr)  ! for SGI
!      IF ( .not. lsgi ) tim1 = second()  !  for cray

! do fft solve on work grid    
#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'call boxmg, ny = ', ny 
#endif
      IF ( ny .gt. 2 ) THEN ! 3D case

#ifdef MPI
      
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt2 = 0
         dt1 = MPI_Wtime()

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!         write(iunit,*) 'call BMG3_SymStd_SOLVE_boxmg'
!         write(iunit,*) 'sizes of rho: ',size(rho,1),size(rho,2),size(rho,3)
!         write(iunit,*) 'sizes of rhow: ',size(rhow,1),size(rhow,2),size(rhow,3)
!         write(iunit,*) 'sizes of rhs: ',size(rhs,1),size(rhs,2),size(rhs,3)
!         write(iunit,*) 'sizes of rhs0: ',size(rhs0,1),size(rhs0,2),size(rhs0,3)

!         CALL BMG3_SymStd_SOLVE_boxmg(                                   

!         BMG_iPARMSlg(id_BMG3_SETUPdg) = 0

       BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 2

         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :              my_comm                                              &
     &             my_comm,mask                                              &
     &             )

       BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 1

!      IF ( .not. usebox ) THEN
      DO k = 1,nnzt
      DO j = 1,nny
      DO i = 1,nnx
        rhssor(i,j,k) = -rhs(i,j,k)/(x1d2lgt(i)*y1d2lgt(j)*(dlzlgt/z1d2lgt(k,3)))
      END DO
      END DO
      END DO
!      ENDIF

!
! Turn off future setup:
         BMG_iPARMSlg(id_BMG3_SETUPdg) = 2

!         write(iunit,*) 'done BMG3_SymStd_SOLVE_boxmg'

         CALL MPI_Barrier(my_comm, mpi_error_code)
         dt2 = dt2 + ( MPI_Wtime() - dt1 )

          IF ( my_rank == 0 ) THEN
            WRITE(iunit,*) ' time for boxmg lightning initial solution = ', dT2
!            WRITE(0,*) ' time for boxmg lightning initial solution = ', dT2
            dt2 = 0
          ENDIF

!         IF (my_rank .eq. 1) THEN
!          CALL PRINTQ( SOlg, QFlg, Qlg,                      
!     :               NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg,   
!     :               iGslg, jGslg, kGslg, dx, dy, dlz, 0    
!     :              )

!         ENDIF
         
#endif

      ELSE ! 2D code

       write(0,*) '2D lightning not available yet with BOXMG -- might never be here'
       call commasmpi_abort()
       
      ENDIF ! ny .gt. 2
      
  55  CONTINUE


#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'copy phi2 into phip'
#endif
      
      DO k = 1,nnzt
      DO j = 1,nny
      DO i = 1,nnx
!        phip(i,j,k) =  phi2(i,j,k)
        phi2(i,j,k) = phip(i,j,k)
      END DO
      END DO
      END DO

#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'set phip = 0 at ground'
#endif

      DO j = 1,nny
      DO i = 1,nnx
        phip(i,j,0) =  0.0
      END DO
      END DO
      
! turn this off for MPI (2015/02/23)
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF



       pstop = 0

#if DEBUG > 0
      IF ( ndebug .ge. 2 ) THEN
      j = Max(1, nys/2)
      write(iunit,*) 'k,z,phi,rho through nx/2,ny/2, i,j=',nxs/2,j
      DO k=nzs-1,2,-1
         
        write(iunit,*) k,z1d2lgt(k,1),phi0(nxs/2,j,k),      &
     &   -0.5/dlzlgt*(phi0(nxs/2,j,k+1)-phi0(nxs/2,j,k-1))*      &
     &     z1d2lgt(k,3),      &
     &   rhow(nxs/2,j,k)
      ENDDO

      write(iunit,*)
      
      IF ( ny .gt. 2 ) THEN
      ix = nxs/2
      jy = Min(ny,54)
      i1 = ix
      j1 = jy
      write(iunit,*) 'k,z,phi,rho through ix,jy= ',ix,jy
      gdstr(:,:,:,:) = 1.0
      IF ( hstretch ) THEN
      ! do not need gdstr for this bit because only vertical is considered
        gdstr(1,0,0,2) = gdx(i1,2)
        gdstr(0,1,0,2) = gdy(j1,2)
        gdstr(1,1,0,2) = gdxy(i1,j1,2)
        gdstr(1,0,1,2) = gdxz(i1,k1,2)
        gdstr(0,1,1,2) = gdyz(j1,k1,2)
        gdstr(1,1,1,2) = gdxyz(i1,j1,k1,2)
        ! gdstr(Abs(i),Abs(j),Abs(j),2)
      ENDIF
      DO k=nzs-1,2,-1
        k1 = k
        rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
         k1 = k
         i = 0
         j = 0
         kz = 1
            e = gd(i,j,kz,k1+kz,2)/2.0*(       &
     &                - 1.0*phip(i1+2*i,j1+2*j,k1+2*kz)       &
     &                + 4.0*phip(i1+i,j1+j,k1+kz)       &
     &                - 3.0*phip(i1,j1,k1) )
           e1 = gd(i,j,kz,k1+kz,2)*      &
     &            (phip(i1+i,j1+j,k1+kz) -  phip(i1,j1,k1))


        write(iunit,*) k,z1d2lgt(k,1),phip(ix,jy,k),      &
     &   -0.5/dlzlgt*(phip(ix,jy,k+1)-phip(ix,jy,k-1))*      &
     &     z1d2lgt(k,3),e,e1,      &
     &   rho(ix,jy,k),rhotemp
      ENDDO
      ENDIF
      
      ENDIF ! ndebug
#endif

#if DEBUG > 0
       IF ( ndebug >= 1 ) write(iunit,*) 'calling mxfldstrz, dx,dy,dlzlgt=',dx,dy,dlzlgt
#endif
       
!      SUBROUTINE MXFLDSTRZ(nx,ny,nz,
!     :    phi,gd,pstop,nbw,nbs,
!     :    emax,emagmax,i1,j1,k1,i2,j2,k2,dx,dy,dz,istretch,
!     :    z1d2,nz1d)

#if DEBUG > 0
      IF (ndebug .ge. 1) write(iunit,*) 'call mxfldstrz'
#endif

       CALL mxfldstrzmpi(nnx,nny,nnzt,      &
     &    phip,gd,gdx,gdy,pstop,nbw,nbs,      &
     &    emax,emagmax,i1,j1,k1,i2,j2,k2,ds,ds,dz,istretch,      &
     &    z1d2lgt,nzmax1lgt,ezmax,iunit,bcx,bcy)
      
!      write(0,*) 'done call mxfldstrz',my_rank

#ifdef MPI

#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,1pe12.5,a,3(1x,i5))')      &
     &   'emagmax = ',emagmax,' at ijk ',i1,j1,k1
      write(iunit,'(a,1pe12.4)')      &
     &   'ezmax = ',ezmax
#ifndef MPI
      write(iunit,'(a,1pe11.5,a,3(1x,i5))')      &
     &   'emax = ',emax,' at ijk ',i2,j2,k2
#endif
      ENDIF
      
      IF ( emagmax .gt. 1.0e6 .or. pstop .eq. 1 ) THEN
        IF ( pstop .eq. 0 ) THEN
           write(iunit,*) 'Efield too large!  Stopping now...'
        ENDIF
        DO kz=nzs-1,1,-1
          write(iunit,*) kz,phi(i1,j1,kz),rho(i1,j1,kz),      &
     &     rhs(i1,j1,kz)
        ENDDO
        i1=nnx/2
        j1=nny/2
        DO kz=nzs-1,1,-1
          write(iunit,*) kz,phi(i1,j1,kz),rho(i1,j1,kz),      &
     &     rhs(i1,j1,kz)
        ENDDO
       
        call commasmpi_abort()
      ENDIF
      

! make list of possible initiation points
      ninit = 0
      
      IF ( icldtop .ne. 0 ) THEN
        DO kz=1,nzs
          IF ( hgt(kz) .lt. zcldtop ) THEN
            k2 = kz
            icldtop = kz
          ENDIF
        ENDDO
        k2 = Min( k2, nzs-4 )
!        k2 = Min( icldtop, nzs-4 )
      ELSE
        k2 = nzs-4
      ENDIF
      
      IF ( iliter .lt. 40 .or. isa .ge. 1 ) THEN
      
      ist = Max(1,nbw+1)
      jst = Max(1,nbs+1)
      jend = Max(1, nny-nbn-1)

      ixe = NLx
      IF ( myproci == nproci ) ixe = nnx - nbe - 1
      
      jye = NLy
      IF ( myprocj == nprocj ) jye = Max(1, nny - nbn - 1)

      DO kz = 3, k2 ! nzs-4
       emagmax = 0
       ezmax = 0
       ehmax = 0
       x = 0
       y = 0
      DO jy = jst,jye
      DO ix = ist,ixe ! nnx-nbe-1
         IF ( ix .gt. nbw+1 ) THEN
           ex  = (phip(ix-1,jy,kz) - phip(ix+1,jy,kz))*dxinv2
         ELSE
           ex  = (phip(ix,jy,kz) - phip(ix+1,jy,kz))*dxinv
         ENDIF
          ex1 = (phip(ix-1,jy  ,kz) - phip(ix,jy,kz))*dxinv
         
         IF ( ny .gt. 2 ) THEN
          IF ( jy .gt. nbs+1 ) THEN
            ey  = (phip(ix,jy-1,kz) - phip(ix,jy+1,kz))*dyinv2
          ELSE
            ey  = (phip(ix,jy,kz) - phip(ix,jy+1,kz))*dyinv
          ENDIF
          ey1 = (phip(ix  ,jy-1,kz) - phip(ix,jy,kz))*dyinv
         ELSE
           ey = 0.0
         ENDIF
         IF ( istretch ) THEN
           ez  = (phip(ix,jy,kz-1) - phip(ix,jy,kz))*dzinv*z1d2lgt(kz,3)
!           ez1 = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2*z1d2lgt(kz,3)
           ez1 = ez
         ELSE
           ez = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2
           ez1 = (phip(ix,jy,kz-1 ) - phip(ix,jy,kz))*dzinv
         ENDIF
         emag  =  Sqrt(ex**2 + ey**2 + ez**2) 
!         emag1 = Sqrt(ex1**2 + ey1**2 + ez1**2)
!         emag  = Max( emag1, emag )
!         IF ( emag .gt. emagmax ) THEN
!           emagmax = emag
!           ezmax = ez
!           ehmax = Sqrt(ex**2 + ey**2)
!           x = (ds*Real(ix-nbw-1 - 1))/dx1 + 1
!           y =  (ds*Real(jy-nbs-1 - 1))/dy1 + 1
!         ENDIF

       IF (emag .gt. efracinitwire*overvolt*elgtfdel*ecritz(kz) ) THEN
!        write(iunit,*) 'kz, emag, elgtfdel,ectritz(kz) = ',kz,emag, elgtfdel,ecritz(kz) 
!        write(0,*) 'kz, emag, elgtfdel,ectritz(kz) = ',kz,emag, elgtfdel,ecritz(kz) 
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(0,*) 'Ee-gad! Over-ran the ipick array! ABORT!',npick1,ninit,my_rank
        call commasmpi_abort()
      END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
!       IF ( iliter .eq. 1 ) THEN
!        write(iunit,'(2(f4.1,1x),i4,2x,3(1x,1pe12.4))') 
!     :   x,y,kz,emagmax,ezmax,ehmax
!       ENDIF
      ENDDO  ! k loop
      ELSE
      ist = Max(1,nbw+1)
      jst = Max(1,nbs+1)
      jend = Max(1, nny-nbn-1)
      DO kz = 3,k2  !  nzs-4
      DO jy = jst,jend
      DO ix = ist,nnx-nbe-1
      
         IF ( ix .gt. nbw+1 ) THEN
           ex  = (phip(ix-1,jy,kz) - phip(ix+1,jy,kz))*dxinv2
         ELSE
           ex  = (phip(ix,jy,kz) - phip(ix+1,jy,kz))*dxinv
         ENDIF

         IF ( ny .gt. 2 ) THEN
          IF ( jy .gt. nbs+1 ) THEN
            ey  = (phip(ix,jy-1,kz) - phip(ix,jy+1,kz))*dyinv2
          ELSE
            ey  = (phip(ix,jy,kz) - phip(ix,jy+1,kz))*dyinv
          ENDIF
         ELSE
            ey = 0.0
         ENDIF

         IF ( istretch ) THEN
           ez = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2*z1d2lgt(kz,3)
         ELSE
           ez = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2
         ENDIF
         emag = Sqrt(ex**2 + ey**2 + ez**2)

       IF (emag .gt. efracinitwire*overvolt*ecritz(kz) ) THEN
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(0,'(a)') 'Ee-gad! Over-ran the ipick array! ABORT!'
        call commasmpi_abort()
      END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
      ENDDO  ! k loop
      END IF

#ifndef MPI 
      write(iunit,'(1x,i10,a)') ninit, ' possible starting points'
#endif

#ifdef MPI

! send e-field components of emag max
      nz1 = 2

!      integer,allocatable,save :: init_send_buffer(:)
!      integer,allocatable,save :: init_recv_buffer(:,:)
!      integer,allocatable,save :: mpi_init0(:,:)
      
!      integer :: init_send_size, init_recv_size, nz1
      
      IF ( .not. allocated(init_recv_buffer) ) THEN
      allocate(init_send_buffer(2))
      allocate(init_recv_buffer(number_of_processes,2))

      allocate(mpi_init0(2,0:number_of_processes-1))
      
      ENDIF
      
!      DO kz = 1,nz1
        init_send_buffer(1) = my_rank
        init_send_buffer(2) = ninit
!      ENDDO

      init_send_size = nz1
      init_recv_size = nz1

      CALL MPI_Gather(init_send_buffer,init_send_size,MPI_INTEGER,      &
     &               init_recv_buffer,init_recv_size,MPI_INTEGER,       &
     &               0,my_comm,mpi_error_code)

       IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpi_init0 = reshape(init_recv_buffer,(/2,number_of_processes/))

           ninittot = 0
           DO kz = 0,number_of_processes-1
       
!            irank1 = int(mpi_init0(1,kz))

            ninittot = ninittot +  mpi_init0(2, kz)

          ENDDO
          
!          write(iunit,*) 'ninittot = ',ninittot
       write(iunit,'(1x,i10,a,i6)') ninittot, ' possible starting points' ! on rank =',my_rank
          

       ENDIF
!       write(0,*) 'spark: starting points ',my_rank

       mpitotinint(1)  = ninit

      CALL MPI_AllReduce(mpitotinint, mpitotoutint, 1, MPI_INTEGER, MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       ninittot = mpitotoutint(1)  
!       write(iunit,*) 'reduced ninittot = ',ninittot
          

#else

       ninittot = ninit

#endif

          IF ( ninittot == 0 ) THEN
!            write(iunit,*) 'SPARKMPI: Something wrong, ninittot = 0'
!            write(0,*) 'SPARKMPI: Something wrong, ninittot = 0'
!         write(iunit,'(a)') 'No more points to choose!'
         lgtstp = 1
         GOTO 999
            
!            lgtstp = 1
!            return
          ENDIF

       
       numtry = 0
 199  CONTINUE
       numtry = numtry + 1
       IF ( numtry .gt. 1 ) THEN ! {
       icgyn = 0
       icgpn = 0
       iendgnd = 0
       imn = nxs-1
       jmn = Max(1, nys-1)
       kmn = nzs-1
       imx = 2
       jmx = 2
       kmx = 2
       imn1 = nxs-1
       jmn1 = Max(1, nys-1)
       kmn1 = nzs-1
       imx1 = 2
       jmx1 = 2
       kmx1 = 2
       IF ( ny .eq. 2 ) THEN
        jmx = 1
        jmx2 = 1
       ENDIF
! clear out mask array 
           DO iend = 1, nend
             DO n=1,numchan(iend)
              IF ( linkrank(n,iend) == my_rank ) THEN
               i1 = chanijk(n,1,iend)
               j1 = chanijk(n,2,iend)
               k1 = chanijk(n,3,iend)
               mask(i1,j1,k1 ) = 1.0e30
               phip(i1,j1,k1) =  phi2(i1,j1,k1) 
               IF ( numchan(iend) .le. 5 .and. my_rank == 0 ) THEN
                 write(iunit,'(5i5)') iend,n, i1,j1,k1
               ENDIF
              ENDIF
             END DO
           END DO
           pickmask(:,:,:) = 0
      DO m=1,nend
      numchan(m) = 0
      DO n=1,nxyz2
      linksrc(n,m) = 0
      linkstep(n,m) = 0
      chandist(n,m) = 0.0
      chanijk(n,1,m) = 0
      chanijk(n,2,m) = 0
      chanijk(n,3,m) = 0
      npbond(n,m) = 0
      phichan(n,m) = 0.0
      ibonds(n,1,1,m) = 0
      ibonds(n,2,1,m) = 0
      ibonds(n,3,1,m) = 0
      ibonds(n,1,2,m) = 0
      ibonds(n,2,2,m) = 0
      ibonds(n,3,2,m) = 0
       DO i=1,25
        bondijk(i,n,1,m) = 0
        bondijk(i,n,2,m) = 0
        bondijk(i,n,3,m) = 0
        bondprob(i,n,m) = 0.0   
       END DO ! i
      END DO ! n
      END DO ! m

            mask(:,:,:) = 1.0e30
            pickmask(:,:,:) = 0

        DO k = 1,nnzt
        DO j = 1,nny
        DO i = 1,nnx
          phip(i,j,k) =  phi2(i,j,k)
!          mask(i,j,k) = 1.0
        END DO
        END DO
        END DO

        DO j = 1,nny
        DO i = 1,nnx
          phip(i,j,0) =  0.0
        END DO
        END DO
         IF ( my_rank == 0 ) THEN
           write(iunit,'(a,i4)') 'Try number ',numtry
         ENDIF
       END IF ! } ( numtry .gt. 1 )
       
       IF ( ninittot .le. 0 ) THEN
        IF ( my_rank == 0 ) THEN
          write(iunit,'(a)') 'No more points to choose!'
        ENDIF
         lgtstp = 1
         GOTO 999
       END IF
!      
! here either choose random point (irand = 0) or use
! a specified starting point (irand = 1)
!
      IF ( irand .ne. 1 .or. isa .lt. 1 ) THEN
        ichoose = int(real(ninittot)*ran0mpi(iseed)) + 1
        IF ( irand .eq. -1 ) ichoose = Min(iliter,ninittot)
        
!        ichoose = 1
        
      IF ( my_rank == 0 ) THEN
        write (iunit,'(2(a,i10))') 'ichoose = ',ichoose
        write(iunit,*) 'iseed = ',iseed
      ENDIF
!       write(0,*) 'spark: ichoose ',my_rank

#ifdef MPI
! find which processor has the initiation point
      IF ( my_rank == 0 ) THEN
       ninittmp = 0
           ichooserank = -1
           ichoosenumber = 0
           DO kz = 0,number_of_processes-1
       
!            write(iunit,*) 'kz = ',kz
            
            IF ( ichoose .le. ninittmp + mpi_init0(2, kz) ) THEN
               write(iunit,*) 'ichooserank = ',mpi_init0(1,kz)
               ichooserank = mpi_init0(1,kz)
               ichoosenumber = ichoose - ninittmp
               EXIT
!            ELSE
!              write(iunit,*) 'init point not on processor ',mpi_init0(1,kz)
            ENDIF
            
            ninittmp = ninittmp + mpi_init0(2, kz)

          ENDDO
          
          IF ( ichooserank == -1 ) THEN
            write(0,*) 'Problem choosing init point! ichooserank = ',ichooserank
          ENDIF

          ichoose_send_buffer(1) = ichooserank
          ichoose_send_buffer(2) = ichoosenumber
          
      ELSE
          ichoose_send_buffer(1) = -1
          ichoose_send_buffer(2) = -1
      ENDIF ! (my_rank == 0 )

       
       CALL MPI_Bcast(ichoose_send_buffer, 2, MPI_INTEGER, 0, my_comm, mpi_error_code)
       ichooserank = ichoose_send_buffer(1)
       ichoose = ichoose_send_buffer(2)
       
#if DEBUG > 0
       IF ( ndebug .ge. 1 )  write(iunit,*) 'my_rank: ichooserank,ichoose = ',my_rank,ichooserank,ichoose
#endif
#endif

       
      IF ( my_rank == ichooserank ) THEN
      
      i1 = ipick(ichoose,1) ! + Igslg - 1
      j1 = ipick(ichoose,2) ! + Jgslg - 1
      k1 = ipick(ichoose,3)
      ixst = i1
      jyst = j1
      kzst = k1
!      write(iunit,'(a,3(1x,i3))') 'Picked starting point ',i1,j1,k1
! update the points array (move up the points lower on the list)     
      DO idum=ichoose,ninit-1
        ipick(idum,1) = ipick(idum+1,1)
        ipick(idum,2) = ipick(idum+1,2)
        ipick(idum,3) = ipick(idum+1,3)
      END DO

          ichoose_send_buffer(1) = i1
          ichoose_send_buffer(2) = j1
          ichoose_send_buffer(3) = k1
          ichoose_send_buffer(4) = igslg0
          ichoose_send_buffer(5) = jgslg0
          ichoose_send_buffer(6) = kgslg0
          ichoose_send_buffer(7) = NLx
          ichoose_send_buffer(8) = NLy
          ichoose_send_buffer(9) = NLz
      
       ninit = ninit - 1
      ENDIF ! my_rank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
#endif
       i1 = ichoose_send_buffer(1)
       j1 = ichoose_send_buffer(2)
       k1 = ichoose_send_buffer(3)
       ixst = i1 + ichoose_send_buffer(4)
       jyst = j1 + ichoose_send_buffer(5)
       kzst = k1 + ichoose_send_buffer(6)
       ioffset = ichoose_send_buffer(4)
       joffset = ichoose_send_buffer(5)
       koffset = ichoose_send_buffer(6)
       NLxchoose = ichoose_send_buffer(7)
       NLychoose = ichoose_send_buffer(8)
       NLzchoose = ichoose_send_buffer(9)
      
      
      ELSE ! either irand = 1 or isa = 1
      
       IF ( my_rank == 0 ) THEN
        write(iunit,*) 'Starting at predetermined point: ',ixst,jyst,kzst
       ENDIF

          ichoose_send_buffer(2) = my_rank

        IF ( igslg <= ixst .and. ixst < igslg + nlx .and.      &
     &       jgslg <= jyst .and. jyst < jgslg + nly .and.      &
     &       kgslg <= kzst .and. kzst < kgslg + nlz  ) THEN

          ichoose_send_buffer(1) = 1
          write(0,*) 'the point is in my_rank = ',my_rank
        ELSE
          ichoose_send_buffer(1) = -1
        ENDIF

#ifdef MPI
        call MPI_ALLREDUCE( ichoose_send_buffer, ichoose_recv_buffer, 1, MPI_2INTEGER,      &
     &                      MPI_MAXLOC, my_comm, mpi_error_code)
#endif
       
        ichooserank = ichoose_recv_buffer(2)
        
        write(0,*) my_rank, 'ichooserank = ',ichooserank,ichoose_recv_buffer(1)
        
      IF ( my_rank == ichooserank ) THEN
      
      i1 = ixst - igslg0 ! ipick(ichoose,1) ! + Igslg - 1
      j1 = jyst - jgslg0 ! ipick(ichoose,2) ! + Jgslg - 1
      k1 = kzst - kgslg0 ! ipick(ichoose,3)
      
      ixst = i1
      jyst = j1
      kzst = k1
!      write(iunit,'(a,3(1x,i3))') 'Picked starting point ',i1,j1,k1

          ichoose_send_buffer(1) = i1
          ichoose_send_buffer(2) = j1
          ichoose_send_buffer(3) = k1
          ichoose_send_buffer(4) = igslg0
          ichoose_send_buffer(5) = jgslg0
          ichoose_send_buffer(6) = kgslg0
          ichoose_send_buffer(7) = NLx
          ichoose_send_buffer(8) = NLy
          ichoose_send_buffer(9) = NLz
      
      ENDIF ! my_rank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
#endif
       i1 = ichoose_send_buffer(1)
       j1 = ichoose_send_buffer(2)
       k1 = ichoose_send_buffer(3)
       ixst = i1 + ichoose_send_buffer(4)
       jyst = j1 + ichoose_send_buffer(5)
       kzst = k1 + ichoose_send_buffer(6)
       ioffset = ichoose_send_buffer(4)
       joffset = ichoose_send_buffer(5)
       koffset = ichoose_send_buffer(6)
       NLxchoose = ichoose_send_buffer(7)
       NLychoose = ichoose_send_buffer(8)
       NLzchoose = ichoose_send_buffer(9)


        ninittot = Max(ninittot,2)
      ENDIF  ! ( irand .ne. 1 .or. isa .ne. 1 )
      
! one fewer init point now      
      ninittot = ninittot - 1
! find neighboring point with largest E and make it the first bond.

      IF ( my_rank == ichooserank ) THEN
      emax = 0.0
      ist = -1 ; isp = 1
      jst = -1 ; jsp = 1
       IF ( ny .eq. 2 ) THEN
        jst = 0
        jsp = 0
       ENDIF
!      IF ( bcy .eq. 2 .and. ( j1 .eq. 1     ) ) jst = 0
!      IF ( bcy .eq. 2 .and. ( j1 .eq. nny-1 ) ) jsp = 0
!      IF ( bcx .eq. 2 .and. ( i1 .eq. 1     ) ) ist = 0
!      IF ( bcx .eq. 2 .and. ( i1 .eq. nnx-1 ) ) isp = 0
 
      gdstr(:,:,:,:) = 1.0
      IF ( hstretch ) THEN
        gdstr(1,0,0,2) = gdx(i1,2)
        gdstr(0,1,0,2) = gdy(j1,2)
        gdstr(1,1,0,2) = gdxy(i1,j1,2)
        gdstr(1,0,1,2) = gdxz(i1,k1,2)
        gdstr(0,1,1,2) = gdyz(j1,k1,2)
        gdstr(1,1,1,2) = gdxyz(i1,j1,k1,2)
        ! gdstr(Abs(i),Abs(j),Abs(j),2)
      ENDIF

      DO k= -1,1
      DO j=  jst,jsp
      DO i=  ist,isp
         etest = Abs(phip(i1+i,j1+j,k1+k) - phip(i1,j1,k1))*      &
     &              gd(i,j,k,k1,2)*gdstr(Abs(i),Abs(j),Abs(j),2)
        
!         write(0,*) 'i,j,k,gd,gdstr = ',i,j,k,gd(i,j,k,k1,2),gdstr(Abs(i),Abs(j),Abs(j),2)
!         IF ( k .ne. 0 .and. istretch ) THEN
!          etest = etest*z1d2lgt(k1,1)
!         ENDIF
         IF (etest .gt. emax) THEN
           emax = etest
           i2 = i1+i
           j2 = j1+j
           k2 = k1+k
         END IF
      END DO ! i
      END DO ! j
      END DO ! k

!
! check whether second point goes into another tile
!
       ichoosenextrank = ichooserank
       IF ( i2 .eq. 0 .or. j2 .eq. 0 .or.       &
     &      i2 .eq. NLx+1 .or. j2 .eq. NLy+1 ) THEN
         ishift = 0
         jshift = 0
         IF ( i2 .eq. 0 )     ishift = -1
         IF ( j2 .eq. 0 )     jshift = -1
         IF ( i2 .eq. NLx+1 ) ishift = 1
         IF ( j2 .eq. NLy+1 ) jshift = 1
         
         ichoosenextrank = procmap(myproci+ishift, myprocj+jshift,myprock)
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
           write(0,*) 'Second lightning point is in another tile! '
           write(0,*) 'i1,j1,i2,j2,ishift,jshift = ',i1,j1,i2,j2,ishift,jshift
           DO j=1,nprocj
           DO i=1,nproci
             write(0,*) 'i,j, procmap = ',procmap(i,j,myprock)
           ENDDO
           ENDDO
           
           write(0,*) 'ichooserank,ichoosenextrank = ',ichooserank,ichoosenextrank
           
         ENDIF
#endif
         IF ( ichoosenextrank .lt. 0 .or. ichoosenextrank .gt. number_of_processes-1) THEN
           write(0,*) 'problem with ichoosenextrank!'
           CALL commasmpi_abort()
           STOP
         ENDIF
       ENDIF
       
        IF (  .not.      &
     &      (i2 + ioffset .ge. igslg .and. i2 + ioffset .le. igslg0 + nxs .and.      &
     &       j2 + joffset .ge. jgslg .and. j2 + joffset .le. jgslg0 + nys ) ) THEN

        ENDIF

         mask(i2,j2,k2) = phip(i2,j2,k2)
         IF ( my_rank == ichoosenextrank) THEN
          rhotemp = -eperao*(      &
     &      cofx(i2,1)*phip(i2-1,j2,k2)+cofx(i2,2)*phip(i2+1,j2,k2) +      &
     &      cofy(j2,1)*phip(i2,j2-1,k2)+cofy(j2,2)*phip(i2,j2+1,k2) +      &
     &      cofz(k2,1)*phip(i2,j2,k2-1)+cofz(k2,2)*phip(i2,j2,k2+1) +      &
     &      (cofx(i2,3)+cofy(j2,3)+cofz(k2,3))*phip(i2,j2,k2)  )
          
          rhs(i2,j2,k2) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k2,3))*gdx(i2,1)*gdy(j1,1)
          
         ENDIF
          

          phistart1 = phip(i1,j1,k1)
          phistart2 = phip(i2,j2,k2) 
          rchoose_send_buffer(1) = i2
          rchoose_send_buffer(2) = j2
          rchoose_send_buffer(3) = k2
          rchoose_send_buffer(4) = phistart1
          rchoose_send_buffer(5) = phistart2
          rchoose_send_buffer(6) = ichoosenextrank
 
        ! gdstr(:,:,:,:) = 1.0
        IF ( hstretch ) THEN
          DO k = 1,2
          rchoose_send_buffer(7+(k-1)*6) = gdx(i2,k)
          rchoose_send_buffer(8+(k-1)*6) = gdy(j2,k)
          rchoose_send_buffer(9+(k-1)*6) = gdxy(i2,j2,k)
          rchoose_send_buffer(10+(k-1)*6) = gdxz(i2,k2,k)
          rchoose_send_buffer(11+(k-1)*6) = gdyz(j2,k2,k)
          rchoose_send_buffer(12+(k-1)*6) = gdxyz(i2,j2,k2,k)
          ENDDO
          ! gdstr(Abs(i),Abs(j),Abs(j),2)
        ELSE
          rchoose_send_buffer(7:18) = 1.0
        ENDIF

      ENDIF ! my_rank == ichooserank

#ifdef MPI
       CALL MPI_Bcast(rchoose_send_buffer, 18, MPI_REAL, ichooserank, my_comm, mpi_error_code)
#endif
       IF ( my_rank /= ichooserank ) THEN
       i2 = Nint(rchoose_send_buffer(1))
       j2 = Nint(rchoose_send_buffer(2))
       k2 = Nint(rchoose_send_buffer(3))
       phistart1 = rchoose_send_buffer(4)
       phistart2 = rchoose_send_buffer(5)
       ichoosenextrank = Nint(rchoose_send_buffer(6))
       ENDIF
      

       idif = i2 - i1
       jdif = j2 - j1
       
      gdstr(:,:,:,:) = 1.0
      IF ( hstretch ) THEN
        DO k = 1,2
        gdstr(1,0,0,k) = rchoose_send_buffer(7+(k-1)*6) 
        gdstr(0,1,0,k) = rchoose_send_buffer(8+(k-1)*6) 
        gdstr(1,1,0,k) = rchoose_send_buffer(9+(k-1)*6) 
        gdstr(1,0,1,k) = rchoose_send_buffer(10+(k-1)*6)
        gdstr(0,1,1,k) = rchoose_send_buffer(11+(k-1)*6)
        gdstr(1,1,1,k) = rchoose_send_buffer(12+(k-1)*6)
        ENDDO
        ! gdstr(Abs(i),Abs(j),Abs(j),2)
      ENDIF

#ifndef MPI       
! checks for periodic case (not being used)
       IF ( bcy .eq. 2 .and. ( j2 .le. 1 .or. j2 .ge. nny-1) )       &
     &    print*,'0 j1,j2,nny = ',j1,j2,nny
       IF ( bcy .eq. 2 .and. ( j2 .le. 1 .or. j2 .ge. nny-1 ) ) THEN
        IF ( j2 .le. 1 .or. j2 .ge. nny-1 ) print*,'1 j2,nny = ',j2,nny
        IF ( j2 .eq. 0   ) j2 = nny - 1
         IF ( j2 .eq. nny ) j2 = 1
        IF ( j2 .le. 1 .or. j2 .ge. nny-1 ) print*,'2 j2,nny = ',j2,nny
       ENDIF

       IF ( bcx .eq. 2 .and. ( i2 .eq. 0 .or. i2 .eq. nnx ) ) THEN
         IF ( i2 .eq. 0   ) i2 = nnx - 1
         IF ( i2 .eq. nnx ) i2 = 1
       ENDIF
#endif
!      IF ( phip(i2,j2,k2) .gt.  phip(i1,j1,k1) ) THEN
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(1x,i5))') 'Picked starting point ',i1+ioffset,j1+joffset,k1+koffset
      ENDIF
        
#if DEBUG > 0
      IF ( ndebug .ge. 1 ) THEN
        write(iunit,'(a,3(1x,i5),a)') 'First bond point      ',      &
     &       i2,j2,k2, ' on my_rank = ',my_rank
      ENDIF
#endif
      
      IF ( phistart2 .gt.  phistart1 ) THEN
       IF ( my_rank == 0 ) THEN
        write(iunit,'(a,3(1x,i5),a)') 'First bond point      ',      &
     &       i2+ioffset,j2+joffset,k2+koffset,' (neg end)'
       ENDIF
        IF ( k2 .gt. k1 ) THEN
         idown = 2 ! pos end going down
        ELSE 
         idown = 1 ! neg end going down
        ENDIF
      ELSE
       IF ( my_rank == 0 ) THEN
        write(iunit,'(a,3(1x,i5),a)') 'First bond point      ',      &
     &       i2+ioffset,j2+joffset,k2+koffset,' (pos end)'
       ENDIF
        IF ( k2 .lt. k1 ) THEN
         idown = 2 ! pos end going down
        ELSE 
         idown = 1 ! neg end going down
        ENDIF
      ENDIF
      
      IF ( k1 .eq. k2 ) idown = 0
      
      IF ( isa .ge. 1 .and. my_rank == 0 ) THEN
        IF ( k1 .eq. k2 ) THEN
          write(iunit,*) 'Neither channel downward'
        ELSEIF ( idown .eq. 1 ) THEN
          write(iunit,*) 'Neg channel downward'
        ELSEIF ( idown .eq. 2 ) THEN
          write(iunit,*) 'Pos channel downward'
        ENDIF
      ENDIF

        tim1 = commas_dtime(timarr)
!       IF ( lsgi ) tim1 = commas_dtime(timarr)
!       IF ( .not. lsgi ) tim1 = second()
        
! start the leader at the location of max efield
      IF (nend .eq. 2) THEN
!       phistart = ( phip(i1,j1,k1) + phip(i2,j2,k2) )/2.0
!       phidiff = Abs( phip(i1,j1,k1) - phip(i2,j2,k2) )
       phistart = ( phistart1 + phistart2 )/2.0
       phidiff = Abs( phistart1 - phistart2 )
       IF ( my_rank == 0 ) THEN
         write(iunit,'(a,3(1x,1pe13.5))') 'phistart,phi1,phi2 = ',      &
     &       phistart, phistart1, phistart2
       ENDIF
       phistart0 = phistart
!         phistlow = Min( phip(i1,j1,k1),phip(i2,j2,k2) )
         phistlow = phistart - 5.e2*gd(idif,jdif,k2-k1,k2,2)* &
                gdstr(Abs(idif),Abs(jdif),Abs(k2-k1),2)*phidiff
!         phisthi  = Max( phip(i1,j1,k1),phip(i2,j2,k2) )
         phisthi  = phistart + 5.e2*gd(idif,jdif,k2-k1,k2,2)* &
                gdstr(Abs(idif),Abs(jdif),Abs(k2-k1),2)*phidiff
         signchan(1) = -1.0
         signchan(2) =  1.0
!       IF (phip(i1,j1,k1) .gt. phistart ) THEN
       IF ( phistart1 .gt. phistart ) THEN
          m = 1
          n = 2
       ELSE
         m = 2
         n = 1 
       END IF
       numchan(m) = 1
       chandist(1,m) = gd(idif,jdif,k2-k1,k2,1)*0.5*gdstr(Abs(idif),Abs(jdif),Abs(k2-k1),1)
       phichan(1,m) = phistart ! - chandist(1,m)*signchan(m)*eintz(k1)
       
       IF ( my_rank == ichooserank ) THEN 
         phip(i1,j1,k1) = phichan(1,m)

!         IF ( usebox ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(iunit,*) 'calc rhotemp at i1,j1,k1,my_rank = ',i1,j1,k1,my_rank
         ENDIF
#endif
          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          
          rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdstr(1,0,0,1)*gdstr(0,1,0,1)
          
!         ENDIF
       
       ENDIF

       chanijk(numchan(m),1,m) = i1
       chanijk(numchan(m),2,m) = j1
       chanijk(numchan(m),3,m) = k1
       ibonds(numchan(m),1,1,m) = i1
       ibonds(numchan(m),2,1,m) = j1
       ibonds(numchan(m),3,1,m) = k1
       ibonds(numchan(m),1,2,m) = i1
       ibonds(numchan(m),2,2,m) = j1
       ibonds(numchan(m),3,2,m) = k1
       linksrc(1,m) = 0
       linkstep(1,m) = 0
       
       linkrank(1,m) = ichooserank
       linkoffset(1,1,m) = ioffset
       linkoffset(1,2,m) = joffset
       
       numchan(n) = 1
       chandist(1,n) = gd(idif,jdif,k2-k1,k2,1)*0.5*gdstr(Abs(idif),Abs(jdif),Abs(k2-k1),1)
       phichan(1,n) = phistart ! - chandist(1,n)*signchan(n)*eintz(k2)

      IF ( my_rank == ichooserank ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(iunit,*) 'calc rhotemp at i1,j1,k1,my_rank = ',i1,j1,k1,my_rank
         ENDIF
#endif
        mask(i1,j1,k1) = phistart
!         IF ( usebox ) THEN
          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          
          rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdstr(1,0,0,1)*gdstr(0,1,0,1)
          
!         ENDIF
      ENDIF

       IF (       &
     &      (i1 + ioffset .ge. igslg0 .and. i1 + ioffset .le. igslg0 + NLx+1 .and.      &
     &       j1 + joffset .ge. jgslg0 .and. j1 + joffset .le. jgslg0 + NLy+1 ) ) THEN
        
        mask(i1 + ioffset - igslg0, j1 + joffset - jgslg0, k1) = phichan(1,n)
        phip(i1 + ioffset - igslg0, j1 + joffset - jgslg0, k1) = phichan(1,n)
        
      ENDIF

!
! check whether second point goes into another tile or its ghost zone
!
!      IF ( my_rank == ichooserank ) THEN
       IF (       &
     &      (i2 + ioffset .ge. igslg0 .and. i2 + ioffset .le. igslg0 + NLx+1 .and.      &
     &       j2 + joffset .ge. jgslg0 .and. j2 + joffset .le. jgslg0 + NLy+1 ) ) THEN
        
        mask(i2 + ioffset - igslg0, j2 + joffset - jgslg0, k2) = phichan(1,n)
        phip(i2 + ioffset - igslg0, j2 + joffset - jgslg0, k2) = phichan(1,n)
        
      ENDIF

       IF ( my_rank == ichooserank .and. ichoosenextrank == ichooserank ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(iunit,*) 'calc rhotemp at i2,j2,k2,my_rank = ',i2,j2,k2,my_rank
         ENDIF
#endif
        phip(i2,j2,k2) = phichan(1,n)
        mask(i2,j2,k2) = phichan(1,n)
!         IF ( usebox ) THEN
          rhotemp = -eperao*(      &
     &      cofx(i2,1)*phip(i2-1,j2,k2)+cofx(i2,2)*phip(i2+1,j2,k2) +      &
     &      cofy(j2,1)*phip(i2,j2-1,k2)+cofy(j2,2)*phip(i2,j2+1,k2) +      &
     &      cofz(k2,1)*phip(i2,j2,k2-1)+cofz(k2,2)*phip(i2,j2,k2+1) +      &
     &      (cofx(i2,3)+cofy(j2,3)+cofz(k2,3))*phip(i2,j2,k2)  )
          
          rhs(i2,j2,k2) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k2,3))*gdstr(1,0,0,1)*gdstr(0,1,0,1)
          
!         ENDIF
      
       ENDIF
!
! Check whether first point goes across tile boundary into local tile:
!
       IF ( my_rank /= ichoosenextrank .and.      &
     &      i2 + ioffset .ge. igslg0 .and. i2 + ioffset .le. igslg0 + NLx+1 .and.      &
     &      j2 + joffset .ge. jgslg0 .and. j2 + joffset .le. jgslg0 + NLy+1 ) THEN
       
       phip(i2+ioffset-igslg0,j2+joffset-jgslg0,k2) = phichan(1,n)
       mask(i2+ioffset-igslg0,j2+joffset-jgslg0,k2) = phichan(1,n)
       
       ENDIF


      IF ( ichooserank .ne. ichoosenextrank ) THEN
      
       IF ( my_rank .eq. ichoosenextrank ) THEN
          ichoose_send_buffer(1) = i2 + ioffset - igslg0
          ichoose_send_buffer(2) = j2 + joffset - jgslg0
          ichoose_send_buffer(3) = k2
          ichoose_send_buffer(4) = igslg0
          ichoose_send_buffer(5) = jgslg0
          ichoose_send_buffer(6) = kgslg0
          ichoose_send_buffer(7) = NLx
          ichoose_send_buffer(8) = NLy
          ichoose_send_buffer(9) = NLz
      
       ENDIF ! my_rank == ichoosenextrank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichoosenextrank, my_comm, mpi_error_code)
#endif

          i2 = ichoose_send_buffer(1)
          j2 = ichoose_send_buffer(2)
          k2 = ichoose_send_buffer(3)

       IF ( my_rank ==  ichoosenextrank ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(iunit,*) 'calc rhotemp at i2,j2,k2,my_rank = ',i2,j2,k2,my_rank
         ENDIF
#endif
        phip(i2,j2,k2) = phichan(1,n)
        mask(i2,j2,k2) = phichan(1,n)
!         IF ( usebox ) THEN
          rhotemp = -eperao*(      &
     &      cofx(i2,1)*phip(i2-1,j2,k2)+cofx(i2,2)*phip(i2+1,j2,k2) +      &
     &      cofy(j2,1)*phip(i2,j2-1,k2)+cofy(j2,2)*phip(i2,j2+1,k2) +      &
     &      cofz(k2,1)*phip(i2,j2,k2-1)+cofz(k2,2)*phip(i2,j2,k2+1) +      &
     &      (cofx(i2,3)+cofy(j2,3)+cofz(k2,3))*phip(i2,j2,k2)  )
          
          rhs(i2,j2,k2) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k2,3))*gdx(i2,1)*gdy(j2,1)
          
!         ENDIF
      
       ENDIF
       
       ENDIF
       
       chanijk(numchan(n),1,n) = i2
       chanijk(numchan(n),2,n) = j2 
       chanijk(numchan(n),3,n) = k2
       ibonds(numchan(n),1,1,n) = i2
       ibonds(numchan(n),2,1,n) = j2
       ibonds(numchan(n),3,1,n) = k2
       ibonds(numchan(n),1,2,n) = i2
       ibonds(numchan(n),2,2,n) = j2
       ibonds(numchan(n),3,2,n) = k2
       linksrc(1,n) = 0
       linkstep(1,n) = 0

       linkrank(1,n) = ichoosenextrank
       linkoffset(1,1,n) = ioffset
       linkoffset(1,2,n) = joffset

!      IF ( my_rank == ichooserank ) THEN
!        mask(i1,j1,k1) = 0.0
!        mask(i2,j2,k2) = 0.0
!      ENDIF
       
       iadd = NInt(4000./ds)
       iaddz = NInt(3000./dz)
       imn = Min(i1,i2)
       jmn = Min(j1,j2)
       kmn = Min(k1,k2)
       imx = Max(i1,i2)
       jmx = Max(j1,j2)
       kmx = Max(k1,k2)

       imn = Min(imn,Max(1,i1-iadd))
       jmn = Min(jmn,Max(1,j1-iadd))
       kmn = Min(kmn,Max(1,k1-iaddz))
       imx = Max(imx,Min(nxs+nbw,i1+iadd))
       jmx = Max(jmx,Min(nys+nbs,j1+iadd))
       kmx = Max(kmx,Min(nzs,k1+iaddz))
       
       IF ( ny .eq. 2 ) THEN ; jmn = 1 ; jmx = 1 ; ENDIF
       
       IF ( bcx .eq. 2 ) THEN
         IF ( imn .eq. 1   ) imx = nxs+nbw
         IF ( imx .eq. nxs+nbw ) imn = 1
       ENDIF

       IF ( bcy .eq. 2 ) THEN
!         IF ( jmn .le. 1 .or. jmx .ge. nys ) THEN
!           print*, '1 jmn,jmx = ',jmn,jmx,nys,nbs
!         ENDIF
         IF ( jmn .eq. 1   ) jmx = nys+nbs
         IF ( jmx .ge. nys+nbs ) jmn = 1
!           print*, '2 jmn,jmx = ',jmn,jmx,nys,nbs
       ENDIF

       imn1 = Min(i1+ioffset,i2+ioffset)
       jmn1 = Min(j1+joffset,j2+joffset)
       kmn1 = Min(k1,k2)
       imx1 = Max(i1+ioffset,i2+ioffset)
       jmx1 = Max(j1+joffset,j2+joffset)
       kmx1 = Max(k1,k2)
      
      END IF ! (nend .eq. 2)
      
!      ENDIF ! my_rank == ichooserank
#ifdef MPI
        CALL MPI_Barrier(my_comm, mpi_error_code)
        
        nSTRZSORmpi = 0
        nboxmgiter = 0
        numitstot = 0

!         CALL BMG3_SymStd_SOLVE_boxmg(                                   
        IF ( usebox ) THEN
#if DEBUG > 0
        IF ( ndebug .ge. 1 ) THEN
         write(iunit,*) 'calling commas_SymStd_SOLVE_boxmg on my_rank = ',my_rank
         write(iunit,*) 'rank, size of mask = ',my_rank, size(mask,1),size(mask,2),size(mask,3)
         maskmn = 1.e31
         maskmx = 0
      do k=1, NLz-1
         do j=1, NLy-1
            do i=1, NLx-1
               maskmn = Min( maskmn, mask(i,j,k) )
               maskmx = Max( maskmx, mask(i,j,k) )
               
               IF ( mask(i,j,k) .lt. 0.9e30 ) THEN
                 write(iunit,*) 'i,j,k,mask = ',i,j,k,mask(i,j,k)
               ENDIF
            enddo
         enddo
      enddo
         
         write(iunit,*) 'rank,maskmx,maskmn = ',my_rank,maskmx,maskmn
        ENDIF ! ndebug
#endif
         
         dt1 = MPI_Wtime()

         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :             my_comm                                              &
     &             my_comm,mask                                              &
     &             )

         nboxmgiter = nboxmgiter + 1
!         dt2 = MPI_Wtime()
         dt2 = dt2 + ( MPI_Wtime() - dt1 )

!            WRITE(iunit,*) ' time for boxmg lightning = ', dT2 - dT1

         ELSE
         dt1 = MPI_Wtime()

         i0 = Min(i1,i2)
         j0 = Min(j1,j2)
         k0 = Min(k1,k2) 
        itype = 0
         dt1 = MPI_Wtime()
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(0,*) 'call STRZSORmpi: nnx,nny,nnzt,nzs = ',nnx,nny,nnzt,nzs
          write(0,*) 'call STRZSORmpi: nlx,nly,nlz = ',nlx,nly,nlz
         ENDIF
#endif
        CALL STRZSORmpi(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,2*mdel,      &
     &       cofx,cofy,cofz, bcx, bcy,ibcn,ibcs,ibce,ibcw,      &
     &       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps0,ndebug,numits)
         
         nSTRZSORmpi = nSTRZSORmpi + 1
         numitstot = numitstot + numits
         
         dtSTRZSORmpi = ( MPI_Wtime() - dt1 )
         delt =  MPI_Wtime() - dt1
         
         IF ( my_rank == 0 ) THEN
           write(iunit,*) 'time for initial STRZSORmpi = ',delt
         ENDIF

       ENDIF
         
         CALL MPI_Barrier(my_comm, mpi_error_code)
#endif
             BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 1
       
      
       
!      omega=1.45
!
!      i0 = Min(i1,i2)
!      j0 = Min(j1,j2)
!      k0 = Min(k1,k2) 
!      itype = 0
!      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 1'
!      IF ( ny .gt. 2 ) THEN
!      CALL strzsor(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,2*mdel,
!     :       cofx,cofy,cofz, bcx, bcy,
!     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
!      ELSE
!      CALL strzsor2d(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,2*mdel,
!     :       cofx,cofy,cofz, bcx, bcy,
!     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
!      ENDIF
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF

      IF (ntest .eq. 2) call commasmpi_abort()
      

! start loop for adding points to channel, relaxing, doing it again
! nstop indicates max number of times to repeat
#if DEBUG > 0
      IF (ndebug .ge. 1) write(iunit,*) 'start propagation'
#endif

      DO iend = 1,nend
        istop(iend) = 0
        istopcnt(iend) = 0
      END DO
      
      ihold = 0
      iextra = 0
      
!      DO nstep = 1,nstop
      nstep = 1
!      LightLoop: DO WHILE ( nstep .lt. nstop - maxpicks )
      LightLoop: DO WHILE ( Maxval(numchan(1:2)) .lt. nstop - maxpicks )
      
      
      
      IF (icheck .eq. 2) THEN
        IF (nstep .ge. nstop2) THEN
         IF ( my_rank == 0 ) THEN
           write(iunit,*) 'stopping second end now....'
         ENDIF
         GOTO 201
        END IF
      END IF
            
      numchanmin = Minval( numchan )
      
      DO iend = 1,nend
! look for possible bonds....
! At 'positive' end (end that is going toward higher potential, which is
! really the negative leader since it carries negative charge), the 
! potential at a neighboring point should be _higher_ in order to be considered
! and vice versa for the 'negative' end (seeks neg. charge and lower potential).
!  So, look for   -signchan(iend)*(bondphi-chanphi) .gt. 0.0, and actually want
!  -signchan(iend)*(bondphi-chanphi)/ds .gt. ecrit.  
!  This should prevent travel in an unfavored direction.
!
!       IF ( bcy .eq. 2 ) THEN
!         mask(imn1:imx1,0,kmn1:kmx1) = mask(imn1:imx1,nys-1,kmn1:kmx1)
!       ENDIF

       IF (istop(iend) .le. 1 .and. ihold .ne. iend) THEN
       probtot = 0.0
       DO n=1,numchan(iend)
         npbond(n,iend) = 0
!
! check neighboring points if channel point is on the domain, otherwise skip it
!
         IF ( my_rank == linkrank(n,iend) ) THEN
         i1 = chanijk(n,1,iend)
         j1 = chanijk(n,2,iend)
         k1 = chanijk(n,3,iend)
#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,*)      &
     &     'i1,j1,k1: ',i1,j1,k1
#endif
         
         jst = -1 ; jend = 1
         IF ( ny .eq. 2 ) THEN ; jst = 0 ; jend = 0 ; ENDIF

!         IF ( .not. btest( pickmask(i1,j1,k1), iend) ) THEN ! do not allow growth from points near a current pick (for multiple picks before relaxation)
      gdstr(:,:,:,:) = 1.0
      IF ( hstretch ) THEN
        gdstr(1,0,0,2) = gdx(i1,2)
        gdstr(0,1,0,2) = gdy(j1,2)
        gdstr(1,1,0,2) = gdxy(i1,j1,2)
        gdstr(1,0,1,2) = gdxz(i1,k1,2)
        gdstr(0,1,1,2) = gdyz(j1,k1,2)
        gdstr(1,1,1,2) = gdxyz(i1,j1,k1,2)
        ! gdstr(Abs(i),Abs(j),Abs(j),2)
      ENDIF

         DO k= -1,1
         DO j= jst,jend
         DO i= -1,1
!           IF ( j1 +j  .le. 0 .and. ndebug .ge. 1 ) THEN
!             write(iunit,*) 'Wait, j1+j = ',j1+j, j1,j
!             write(iunit,*) 'my_rank, linkrank(n,iend)= ',my_rank, linkrank(n,iend) 
!           ENDIF

!
!  add condition for pickmask here: check for growth toward a current pick-
!

       pickmaskcheck = .true.
       
       ! need to set some threshold size of flash before allowing multiple picks per channel end
       IF ( numchan(iend) > 10 ) THEN
!          pickmaskcheck = .not. btest( pickmask(i1+i,j1+j,k1+k), iend) ! btest is true when bit flag is set, so the .not. makes a zero bit OK for propagation
       ENDIF

#ifdef MPI
           IF ( pickmaskcheck .and. mask(i1+i,j1+j,k1+k) .gt. 0.9e30 .and.      &
     &        ( ( i1+i .ge. nbw+1 .or. myproci .gt. 1 )  .and.       &
     &           (i1+i .lt. nxs+nbw .or. myproci .lt. nproci )) .and.      &
     &        ( ( j1+j .ge. nbs+1 .or. myprocj .gt. 1 ) .and.       &
     &          ( j1+j .lt. nys+nbs .or. myprocj .lt. nprocj )) .and.      &
     &         k1+k .gt. 0 .and. k1+k .lt. nzs-1 )    THEN
#else
           IF (mask(i1+i,j1+j,k1+k) .gt. 0.9e30 .and.      &
     &        (bcx .eq. 2 .or. ( i1+i .ge. nbw+1 .and. i1+i .le. nxs+nbw)) .and.      &
     &        ( (bcy .eq. 2 .or. ny .eq. 2) .or. ( j1+j .ge. nbs+1 .and. j1+j .le. nys+nbs)) .and.      &
     &         k1+k .gt. 0 .and. k1+k .lt. nzs-1 )    THEN
#endif
!            IF ( ( Abs(rho(i1+i,j1+j,k1+k)) .gt. 1.0e-16 ) .or.
!     :           k1 .lt. 10  ) THEN
              ecomp = fefac(iend)*fecritn(iend)*ecritz(k1)
!            ELSE
!              ecomp = Max(fecritn(iend)*ecritz(k1),
!     :                     Min(125.0e3, 2.0*ecritz(k1) ))
!            ENDIF
            
            edgecheck = ( i1+2*i .ge. 0 .and. i1+2*i .le. NLx+1 .and.      &
     &                    j1+2*j .ge. 0 .and. j1+2*j .le. NLy+1 )
            
            IF ( edgecheck ) THEN
              edgecheck = edgecheck .and. mask(i1+2*i,j1+2*j,k1+2*k) .gt. 0.9e30 
            ENDIF
!           IF ( j1 +2*j  .le. 0 ) THEN
!             write(iunit,*) 'Wait, j1+2*j = ',j1+2*j, j1,j,edgecheck
!             write(iunit,*) 'my_rank, linkrank(n,iend)= ',my_rank, linkrank(n,iend) 
!           ENDIF
            
           IF ( l2nde .ge. 2 .and. edgecheck   &
!     :           .and.  mask(i1+2*i,j1+2*j,k1+2*k) .gt. 0.9e30       &
     &            ) THEN
            e = gd(i,j,k,k1+k,2)*gdstr(Abs(i),Abs(j),Abs(j),2)/2.0*(       &
     &                - 1.0*phip(i1+2*i,j1+2*j,k1+2*k)       &
     &                + 4.0*phip(i1+i,j1+j,k1+k)       &
     &                - 3.0*phichan(n,iend)  )
            IF ( l2nde .eq. 3 ) THEN
              e1 = gd(i,j,k,k1+k,2)*gdstr(Abs(i),Abs(j),Abs(j),2)*      &
     &            (phip(i1+i,j1+j,k1+k) -  phichan(n,iend))
!              e = 0.5*( e + e1 )
            ENDIF
           ELSE
            e = gd(i,j,k,k1+k,2)*gdstr(Abs(i),Abs(j),Abs(j),2)*      &
     &            (phip(i1+i,j1+j,k1+k) -  phichan(n,iend))
           ENDIF
            
            temp = Max(0.0d0,       &
     &      (  - signchan(iend)*e  - ecomp )  )
     
           IF ( temp .ne. 0.0 .and. l2nde .eq. 3 ) THEN
             write(iunit,*) 'e,e1: ',e,e1
           ENDIF

           IF ( temp .gt. 0.0 ) THEN
             npbond(n,iend) = npbond(n,iend) + 1
             bondijk(npbond(n,iend),n,1,iend) = i1+i
             bondijk(npbond(n,iend),n,2,iend) = j1+j
             bondijk(npbond(n,iend),n,3,iend) = k1+k
           bondprob(npbond(n,iend),n,iend) = (temp)**ltgeta
           probtot = probtot + bondprob(npbond(n,iend),n,iend)
           ENDIF
#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,*)      &
     &     'n,npbond: ',n,npbond(n,iend),phip(i1+i,j1+j,k1+k)
#endif
           END IF
         END DO ! i
         END DO ! j
         END DO ! k
!         ENDIF ! pickmask check
         
         ENDIF ! my_rank == linkrank
         
         ! check last point for this end of the channel. If all bond probabilities are zero, then
         ! update to at least one bond with zero bondprob.
         ! but why? bondprob should be initiated to zero
         IF ( npbond(n,iend) .eq. 0 .and. n .eq. numchan(iend) ) THEN
          npbond(n,iend) = 1
          bondprob(npbond(n,iend),n,iend) = 0.0
          IF ( .false. ) THEN
          write(iunit,*) 'Lightning debug: numchan,n,iend',      &
     &           numchan(iend),n,iend
          write(iunit,*) 'i1,j1,k1',i1,j1,k1

          DO k=1,numchan(iend)
           write(iunit,'(3(1x,i4))')       &
     &      (ibonds(k,i,1,iend), i=1,3)
           write(iunit,'(3(1x,i4),1(1x,1pe13.5))')       &
     &      (ibonds(k,i,2,iend), i=1,3),phi(i1,j1,k1)
          ENDDO
          ENDIF
          
         ENDIF
     
       END DO ! n


! calculate probabilities
!       DO n=1,numchan(iend)
!         DO m=1,npbond(n,iend)
!           temp = Max(0.0, 
!     :           Abs(bondphi(m,n,iend) - phichan(n,iend))/ds - ecrit)
!           bondprob(m,n,iend) = (temp)**eta
!           probtot = probtot + bondprob(m,n,iend)
!         END DO ! m
!       END DO ! n

#ifdef MPI
!
!  Now gather the probabilities as for choosing the first point.  
!  Rank 0 process chooses which tile the next point will be on.
!

       mpitotindp(1)  = probtot

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       probtotall = mpitotoutdp(1)  
!       write(iunit,*) 'probtotall = ',probtotall
#else
       probtotall = probtot
#endif
       
#if DEBUG > 0
       IF (ndebug .ge. 1) write(iunit,'(a,2(1x,e12.5))') 'probtot,probtotall = ',probtot,probtotall
#endif
       
       IF (probtotall .lt. 1.0e-10) THEN
#if DEBUG > 0
         IF (ndebug .ge. 1) write(iunit,*)  'No propagation at iend,nstep,rank: ',iend,nstep, my_rank
#endif
         istop(iend) = 1
         istopcnt(iend) = istopcnt(iend) + 1
           IF ( ibal .eq. 2 .and. istopcnt(iend) .gt. 50*500./ds ) THEN
            IF ( my_rank == 0 ) THEN
              write(iunit,"(a,i2,a,i3)")       &
     &          "istopcnt(",iend," ) = ",istopcnt(iend)
            ENDIF
            GOTO 201
           ENDIF
         IF (icheck .eq. 0) THEN
           icheck = 2
           nstop2 = Min(nstop, nstep + nstep/3)
           write(iunit,*) 'New nstop is ',nstop2
         END IF
          IF (nend .eq. 1) GOTO 201
       ELSE
         istop(iend) = 0
         istopcnt(iend) = 0
       END IF
       IF (nend .gt. 1) THEN
       IF (istop(1) .ge. 1 .and. istop(2) .ge. 1 ) THEN
!       IF (istop(1) .ge. 1 .and. istop(2) .ge. 1 .and. 
!     :       ihold .eq. 0) THEN
!         IF (ndebug .ge. 1) 
         IF ( my_rank == 0 ) THEN
         write(iunit,*)  'No propagation at either end...stopping loop.',      &
     &       ' nstep=',nstep
         ENDIF
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) write(iunit,*) ' ihold = ',ihold,       &
     &        'istop1,2 = ',istop(1),istop(2)
#endif
         IF ( ( ((numchan(1) + numchan(2) ) .le. 4 .and. ibal.ne.2)      &
     &       .or. (numchan(1) + numchan(2) ) .le. 2 ) .and.       &
     &       numtry .lt. maxtry  .and. ninittot .ge. 1 ) THEN
           GOTO 199
         ELSE
          GOTO 201
         END IF
       END IF
       END IF 
       
      IF (istop(iend) .eq. 0) THEN 
       
       probcum = 0.0
       DO n=1,numchan(iend)
         DO m=1,npbond(n,iend)
!           bondprob(m,n,iend) = bondprob(m,n,iend)/probtot
           probcum = probcum + bondprob(m,n,iend)
           bondprob(m,n,iend) = probcum
#if DEBUG > 0
       IF (ndebug .ge. 3) write(iunit,'(a,2i5,2x,e12.5)')       &
     &     'n,bond,probcum = ',n,m,probcum
#endif
         END DO ! m
       END DO ! n

#ifdef MPI
!
      nz1 = 2
      IF ( .not. allocated(prob_recv_buffer) ) THEN
      allocate(prob_send_buffer(2))
      allocate(prob_recv_buffer(number_of_processes,2))

      allocate(mpi_prob0(2,0:number_of_processes-1))
      
      ENDIF
      

        prob_send_buffer(1) = my_rank
        prob_send_buffer(2) = probcum
!      ENDDO

      prob_send_size = nz1
      prob_recv_size = nz1

      CALL MPI_Gather(prob_send_buffer,prob_send_size,MPI_DOUBLE_PRECISION,      &
     &               prob_recv_buffer,prob_recv_size,MPI_DOUBLE_PRECISION,       &
     &               0,my_comm,mpi_error_code)

       IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)


       numpick = Min( maxpicks, numchanmin/multipickmod + 1 )
       
 !      numpick = 1
       
       numextrapicks(iend) = -1
       addpick = 0
       badpick = 0
       
       pickloop1 : DO pickloop = 1,numpick
       
        numextrapicks(iend) = numextrapicks(iend) + 1

       
! pull a random number to maintain same sequence on all tiles
        tmp = ran0mpi(iseed) 

       IF (my_rank == 0) THEN

        mpi_prob0 = reshape(prob_recv_buffer,(/2,number_of_processes/))

!        tmp = ran0mpi(iseed) 
        pick = probtotall*tmp

           probcumall = 0
           ichooserank = -1
           DO kz = 0,number_of_processes-1
             IF ( pick .le. probcumall + mpi_prob0(2, kz) .and. ichooserank == -1) THEN
               ichooserank = Nint(mpi_prob0(1,kz))
               
#if DEBUG > 0
               IF ( ndebug .ge. 1 ) THEN
                 write(iunit,*) 'ichooserank for propagation = ',ichooserank,numchan(iend)
                 write(iunit,*) 'iend,numchan1,2 = ',iend,numchan(1),numchan(2)
               ENDIF
#endif
               
             ENDIF

             probcumall = probcumall +  mpi_prob0(2, kz)

           ENDDO
          
!          write(iunit,*) 'probcumall = ',probcumall

          ichoose_send_buffer(1) = ichooserank
       
       ENDIF

       CALL MPI_Bcast(ichoose_send_buffer, 1, MPI_INTEGER, 0, my_comm, mpi_error_code)
       
       ichooserank = ichoose_send_buffer(1)


       mpitotindp(1)  = probcum

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       probcumall = mpitotoutdp(1)  
!       write(iunit,*) 'probcumall = ',probcumall
#else
       probcumall = probcum
       ichooserank = my_rank
#endif


#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,'(a,e12.5)') 'probcum = ',probcum
#endif
     
       IF (probcumall .lt. 0.98*probtotall .or.       &
     &     probcumall .gt. 1.1*probtotall) THEN
         write(iunit,*) 'Somethings gone wrong! probcum out of range: ',      &
     &       probcumall,probtotall
         GOTO 201
       END IF

       tmp = ran0mpi(iseed) 

! 'roll the die' and get a bond within the range of the points within tile
       IF ( my_rank == ichooserank ) THEN !{
       n = numchan(iend)
       n = Max(1,n)
       m = npbond(n,iend)
       m = Max(1,m)
       ntry = 0
 99    CONTINUE
       IF ( ntry > 0 ) tmp = ran0mpi(iseed) 
       pick = probtot*tmp
!       write(0,*) 'tmp,pick = ', tmp,pick,bondprob(m,n,iend)
       
       IF ( .false. ) THEN
       IF ( tmp > 1.0 .or. pick > bondprob(m,n,iend) ) THEN
         ntry = ntry + 1
         IF ( ntry .lt. 100 ) THEN
           write(0,*) 'ntry, myrank = ',ntry,my_rank,tmp,pick,bondprob(m,n,iend)
           GOTO 99
         ELSE
           tmp = 0.0  
           pick = 0.0
         ENDIF
       ENDIF 
       ENDIF 
       
       IF ( ntry .gt. 0 ) THEN
        write(iunit,*) 'Lightning debug: ntry,tmp,pick,probtot = ',      &
     &     ntry,tmp,pick,probtot
       ENDIF
       
       
       IF ( pick .gt. bondprob(m,n,iend) .and. probtot .gt. 0.0 .and.      &
     &   bondprob(m,n,iend) .gt. 0.0   ) THEN
        write(iunit,*)       &
     &   'pick,probtot,probcum,bondprob(m,n,iend),m,n,iend = ',      &
     &   pick,probtot,probcum,bondprob(m,n,iend),m,n,iend
!        GOTO 99
       ELSEIF ( bondprob(m,n,iend) .eq. 0.0 ) THEN
!       write(iunit,*) 'm,n,iend,bondprob(m,n,iend),probcum,probtot = ',      &
!     &   m,n,iend,bondprob(m,n,iend),probcum,probtot
!       write(iunit,*) 'numchan(iend) =', numchan(iend)
!!       DO n=1,numchan(iend)
!!        write(iunit,*) 'n,npbond = ', n, npbond(n,iend)
!!         DO m=1,npbond(n,iend)
!!         write(iunit,*) m,n,bondprob(m,n,iend)
!!         END DO ! m
!!       END DO ! n
       ENDIF
       
       npick = 0
       bpick = 0
#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,'(a,e12.5)') 'pick = ',pick
#endif
       IF ( tmp >= 1.0 ) THEN

         npick = numchan(iend)
         bpick = npbond(n,iend)
       
       ELSE
       
       bondloop1 : DO n=1,numchan(iend)
         DO m=1,npbond(n,iend)
           IF (pick .lt. bondprob(m,n,iend) ) THEN
             npick = n
             bpick = m
             EXIT bondloop1
             ! GOTO 101
           END IF
         END DO ! m
       END DO bondloop1 ! n

       ENDIF

! 101  CONTINUE      

         ichoosenextrank = ichooserank
         IF ( npick .gt. 0 .and. bpick .gt. 0 ) THEN ! {
         
         ! i2,j2,k2 is originating point
          i2 = chanijk(npick,1,iend)
          j2 = chanijk(npick,2,iend)
          k2 = chanijk(npick,3,iend)
          
         ! i1,j1,k1 is new point
          i1 = bondijk(bpick,npick,1,iend)
          j1 = bondijk(bpick,npick,2,iend)
          k1 = bondijk(bpick,npick,3,iend)

       IF ( pickloop > 1  ) THEN
          pickmaskcheck = .not. btest( pickmask(i1,j1,k1), iend) ! btest is true when bit flag is set, so the .not. makes a zero bit OK for propagation
       ELSE
          pickmaskcheck = .true.
       ENDIF


        IF ( pickmaskcheck .eqv. .false. ) THEN
          i1 = i2
          j1 = j2
          k1 = k2
!          npick = -1
        ENDIF

!
! check whether second point goes into another tile
!
         IF ( i1 .eq. 0 .or. j1 .eq. 0 .or.       &
     &      i1 .eq. NLx+1 .or. j1 .eq. NLy+1 ) THEN ! {
         ishift = 0
         jshift = 0
         IF ( i1 .eq. 0 )     ishift = -1
         IF ( j1 .eq. 0 )     jshift = -1
         IF ( i1 .eq. NLx+1 ) ishift = 1
         IF ( j1 .eq. NLy+1 ) jshift = 1
         
         ichoosenextrank = procmap(myproci+ishift, myprocj+jshift,myprock)
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
         write(0,*) 'Next lightning point is in another tile! '
         write(0,*) 'ichooserank,ichoosenextrank = ',ichooserank,ichoosenextrank
         ENDIF
#endif
         ENDIF ! } 
         
         ENDIF ! }
         ichoose_send_buffer(1) = npick
         ichoose_send_buffer(2) = bpick
         ichoose_send_buffer(3) = i1
         ichoose_send_buffer(4) = j1
         ichoose_send_buffer(5) = k1
         ichoose_send_buffer(6) = i2
         ichoose_send_buffer(7) = j2
         ichoose_send_buffer(8) = k2
         ichoose_send_buffer(9) = igslg0
         ichoose_send_buffer(10) = jgslg0
         ichoose_send_buffer(11) = kgslg0
         ichoose_send_buffer(12) = NLx
         ichoose_send_buffer(13) = NLy
         ichoose_send_buffer(14) = NLz
         ichoose_send_buffer(15) = ichoosenextrank
         IF ( pickmaskcheck ) THEN
           ichoose_send_buffer(16) = 1
         ELSE
           ichoose_send_buffer(16) = 0
         ENDIF

       ELSE
  !       tmp = ran0mpi(iseed) 
       ENDIF ! my_rank == ichooserank !}
#if DEBUG > 0
       IF ( ndebug .ge. 1 ) write(iunit,*) 'my_rank,tmp (ran0) = ',tmp
#endif
#ifdef MPI
! send out info on new channel point
       CALL MPI_Bcast(ichoose_send_buffer, 16, MPI_INTEGER, ichooserank,       &
     &                my_comm, mpi_error_code)
#endif
       
       npick = ichoose_send_buffer(1)
       bpick = ichoose_send_buffer(2)

       i1 = ichoose_send_buffer(3)
       j1 = ichoose_send_buffer(4)
       k1 = ichoose_send_buffer(5)
       i2 = ichoose_send_buffer(6)
       j2 = ichoose_send_buffer(7)
       k2 = ichoose_send_buffer(8)
       
       idif = i1-i2
       jdif = j1-j2
       kdif = k1-k2

       ioffset = ichoose_send_buffer(9)
       joffset = ichoose_send_buffer(10)
       koffset = ichoose_send_buffer(11)
       NLxchoose = ichoose_send_buffer(12)
       NLychoose = ichoose_send_buffer(13)
       NLzchoose = ichoose_send_buffer(14)
       
       ichoosenextrank = ichoose_send_buffer(15)
       
       pickmaskcheck = ( ichoose_send_buffer(16) == 1 )
       
       IF ( ( pickmaskcheck .eqv. .false. ) .and. pickloop > 1) THEN
       
        numextrapicks(iend) = numextrapicks(iend) - 1
!        addpick = addpick - 1
        CYCLE
       
       ELSE
       
       ENDIF


      IF ( ichooserank .ne. ichoosenextrank ) THEN
      
       IF ( my_rank .eq. ichoosenextrank ) THEN
          ichoose_send_buffer(1) = i1 + ioffset - igslg0
          ichoose_send_buffer(2) = j1 + joffset - jgslg0
          ichoose_send_buffer(3) = k1
          ichoose_send_buffer(4) = igslg0
          ichoose_send_buffer(5) = jgslg0
          ichoose_send_buffer(6) = kgslg0
          ichoose_send_buffer(7) = NLx
          ichoose_send_buffer(8) = NLy
          ichoose_send_buffer(9) = NLz
      
       ENDIF ! my_rank == ichoosenextrank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichoosenextrank,       &
     &                my_comm, mpi_error_code)
#endif

          i1 = ichoose_send_buffer(1)
          j1 = ichoose_send_buffer(2)
          k1 = ichoose_send_buffer(3)
          ioffset = ichoose_send_buffer(4)
          joffset = ichoose_send_buffer(5)
          koffset = ichoose_send_buffer(6)

       
       ENDIF
!
!  Now we have a new bond to add to the flash!
!
#if DEBUG > 0
       IF (ndebug .ge. 3) write(iunit,'(a,2i5)')       &
     &      'picked link,bond:',npick,bpick
#endif
       IF (npick .eq. 0 .or. bpick .eq. 0) THEN
         write(iunit,*) 'npick or bpick is zero .... stopping'
         write(iunit,*) 'pick, probtot = ',pick,probtot
         write(iunit,*) 'npick,bpick,numchan(1,2) = ',npick,bpick,      &
     &      numchan(1),numchan(2)
         GOTO 201
       END IF
       
       numchan(iend) = numchan(iend) + 1
       linksrc(numchan(iend),iend) = npick
       linkstep(numchan(iend),iend) = nstep + addpick
       addpick = addpick + 1
       linkrank(numchan(iend),iend) = ichoosenextrank

       linkoffset(numchan(iend),1,iend) = ioffset
       linkoffset(numchan(iend),2,iend) = joffset

! ibonds is the source point for the new point
      ibonds(numchan(iend),1,1,iend) = i2
      ibonds(numchan(iend),2,1,iend) = j2
      ibonds(numchan(iend),3,1,iend) = k2
      
     gdstr(:,:,:,:) = 1.0
      IF ( hstretch ) THEN
        gdstr(1,0,0,1) = gdx(i2,1)
        gdstr(0,1,0,1) = gdy(j2,1)
        gdstr(1,1,0,1) = gdxy(i2,j2,1)
        gdstr(1,0,1,1) = gdxz(i2,k2,1)
        gdstr(0,1,1,1) = gdyz(j2,k2,1)
        gdstr(1,1,1,1) = gdxyz(i2,j2,k2,1)
        ! gdstr(Abs(i),Abs(j),Abs(j),2)
      ENDIF

      chandist(numchan(iend),iend) = gd(idif,jdif,kdif,k2,1) +      &
     &                       chandist(npick,iend)   
!     :                       chandist(numchan(iend)-1,iend) 
       phichan(numchan(iend),iend) = phichan(npick,iend) -      &
     &            gd(idif,jdif,kdif,k2,1)*signchan(iend)*eintz(k1)

       IF (       &
     &      (i1 + ioffset .ge. igslg0 .and. i1 + ioffset .le. igslg0 + NLx+1 .and.      &
     &       j1 + joffset .ge. jgslg0 .and. j1 + joffset .le. jgslg0 + NLy+1 ) ) THEN
#if DEBUG > 0
        IF ( ndebug .ge. 1 ) THEN
         write(iunit,*) 'j1 + joffset - jgslg0 = ',j1 + joffset - jgslg0
         write(iunit,*) 'i1 + ioffset - igslg0 = ',i1 + ioffset - igslg0
         write(iunit,*) 'k1 = ',k1
        ENDIF
#endif
        mask(i1 + ioffset - igslg0, j1 + joffset - jgslg0, k1) = phichan(numchan(iend),iend)
        phip(i1 + ioffset - igslg0, j1 + joffset - jgslg0, k1) = phichan(numchan(iend),iend)

      ENDIF
      
       DO j = j1 + joffset - pickmaskwidth, j1 + joffset + pickmaskwidth
       DO i = i1 + ioffset - pickmaskwidth, i1 + ioffset + pickmaskwidth
       IF (       &
     &      (i .ge. igslg0 .and. i .le. igslg0 + NLx+1 .and.      &
     &       j .ge. jgslg0 .and. j .le. jgslg0 + NLy+1 ) ) THEN
        
        ! need to make a loop here to set points around the new bond point. 
        ! Actually need to make an outer loop to see if point is within range of points within this tile
        pickmask(i - igslg0, j - jgslg0, k1) = IBSET( pickmask(i - igslg0, j - jgslg0, k1), iend)

      ENDIF
      ENDDO
      ENDDO
      
      
      IF ( my_rank == ichoosenextrank ) THEN
          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          
          rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdx(i1,1)*gdy(j1,1)
          
         ENDIF

       
       
       
#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,'(a,2(2x,e12.5))')       &
     &      'phip(chan),phip(picked)',phichan(npick,iend),      &
     &       phip(i1,j1,k1)
#endif
       
#ifndef MPI
       IF ( bcy .eq. 2 .and. ( j1 .eq. 0 .or. j1 .eq. nny ) ) THEN
         mask(i1,j1,k1) = phichan(numchan(iend),iend)
         phip(i1,j1,k1) = phichan(numchan(iend),iend)
         IF ( j1 .eq. 0   ) j1 = nny - 1
         IF ( j1 .eq. nny ) j1 = 1
       ENDIF

       IF ( bcx .eq. 2 .and. ( i1 .eq. 0 .or. i1 .eq. nnx ) ) THEN
         mask(i1,j1,k1) = phichan(numchan(iend),iend)
         phip(i1,j1,k1) = phichan(numchan(iend),iend)
         IF ( i1 .eq. 0   ) i1 = nnx - 1
         IF ( i1 .eq. nnx ) i1 = 1
       ENDIF
#endif
       chanijk(numchan(iend),1,iend) = i1
       chanijk(numchan(iend),2,iend) = j1
       chanijk(numchan(iend),3,iend) = k1
      ibonds(numchan(iend),1,2,iend) = i1
      ibonds(numchan(iend),2,2,iend) = j1
      ibonds(numchan(iend),3,2,iend) = k1
       
       
#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,'(a,3i5,a,3i5)')       &
     &      'added bond from ', i2,j2,k2,      &
     &    ' to ',i1,j1,k1
       IF (ndebug .ge. 2) write(iunit,'(a,e13.5)')       &
     &      'new point has phichan = ',phichan(numchan(iend),iend)
#endif

       IF (k1 .ge. nzs-3) THEN
         IF ( my_rank == 0 ) THEN
         write(iunit,*) 'reached the top! stopping now..... ibal= ',ibal,      &
     &      ' icldtop,k1,zcldtop = ',icldtop,k1,zcldtop
         ENDIF
#if DEBUG > 0
         IF ( ndebug .ge. 1 ) THEN
          write(iunit,*) 'goto 201, my_rank = ',my_rank
         ENDIF
#endif
         GO TO 201
! new test version 04.13.2004:
       ELSEIF ( ibal .ge. 2 ) THEN
!        IF ( ndebug .ge. 1 ) write(iunit,*) 'top check, k1 = ',k1
        IF ( icldtop .gt. 0 .and. hgt(k1) .gt. zcldtop + 3000.   &
!     :     .and.  numchan(iend) .gt. numchan(3-iend)      &
     &     ) THEN
         IF ( my_rank == 0 ) THEN
         write(iunit,*) 'reached the top! stopping now..... ibal= ',ibal,      &
     &      ' icldtop,k1,zcldtop = ',icldtop,k1,zcldtop
         ENDIF
         GOTO 201

!       ELSEIF ( ibal .eq. 2 ) THEN
!        IF ( icldtop .gt. 0 .and. k1 .gt. icldtop + 4 .and. 
!     :       numchan(iend) .gt. numchan(3-iend) ) THEN
!         write(iunit,*) 'reached the top! stopping now..... ibal=2 ',
!     :      ' icldtop,k1 = ',icldtop,k1
!         GOTO 201
        ENDIF
       END IF
       
       IF ( iextra .ge. 1 ) THEN
!         IF ( ibal .eq. 2 .and. iextra .ge. 5 ) THEN !  numchan(iend)/3 ) THEN
! 5/4/08 erm test following:
         IF ( ibal .eq. 2 .and. iextra .ge. maxextra ) THEN !  numchan(iend)/3 ) THEN
!         IF ( ibal .eq. 2 .and. iextra .ge. Min(2*ngroundstep,maxextra) ) THEN !  numchan(iend)/3 ) THEN
           IF ( my_rank == 0 ) THEN
             write(iunit,*) 'iextra = ',iextra
             write(iunit,*) 'Max amount of extra growth reached. Stop now.'
           ENDIF
           GOTO 201
         ENDIF
         IF ( ibal .ne. 2  .and. iextra .gt. 50 ) THEN
           IF ( my_rank == 0 ) THEN
             write(iunit,*) 'Max amount of extra growth reached. Stop now.'
           ENDIF
           GOTO 201
         ENDIF
         iextra = iextra + 1
       ENDIF
!       IF ( 0.5*(z1d2lgt(k1,1)+z1d2lgt(Max(1,k1-1),1)) .lt. 0.5*cgthres
!     :     .and. istop(3-iend) .eq. 2 )  THEN
!        IF ( ndebug .ge. 1 ) write(iunit,*) 'CG check, k1 = ',k1
       IF ( 0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. 0.5*cgthres      &
     &     .and. istop(3-iend) .eq. 2 )  THEN
           IF ( my_rank == 0 ) THEN
             write(iunit,*) 'Something screwy! Both ends at ground! ',      &
     &        'Stay with the first and stop now'
           ENDIF
          GOTO 201
       ENDIF
       
!       IF ( ( 0.5*(z1d2lgt(k1,1)+z1d2lgt(Max(1,k1-1),1)) .lt. cgthres  .and. 
!     :      isa .ne. 1 .and. 
!     :        z1d2lgt(kzst,1) .gt. 2001.0 .and. istop(3-iend).le.1) .or.
!     :      ( k1 .le. 1  ) ) THEN
!        IF ( ibal .lt. 4 .or. (ibal .ge. 4 .and. 
!     :       Sign(1.0,phichan(numchan(iend),iend)) .eq. signchan(iend))
!     :      .or. ( k1 .le. 1  )  ) THEN

       IF ( ( 0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and.       &
     &      isa .lt. 2 .and.       &
     &        hgt(kzst) .gt. 2001.0 .and. istop(3-iend).le.1) .or.      &
     &      ( k1 .le. 1  ) ) THEN
        IF ( ibal .eq. 1 .or. ibal.eq.3 .or.       &
     &      ( (ibal .ge. 4 .or. ibal .eq.2) .and.       &
     &       Sign(1.0,phichan(numchan(iend),iend)-      &
     &            hgt(k1)*signchan(iend)*eintz(k1) )       &
     &        .eq. signchan(iend))      &
     &      .or. ( k1 .le. 1  )  ) THEN

         IF ( my_rank == 0 ) THEN
          write(iunit,*) 'reached the ground! stopping now this channel',iend
          write(iunit,*) 'nstep = ',nstep
          ngroundstep = nstep
          write(iunit,'(a,f8.3)') 'fecritn = ',fecritn(iend)
          write(iunit,*) 'touchdown at i,j,k= ',i1,j1,k1
          write(iunit,'(a,2(1x,f8.3))')       &
     &        'fecrit(1,2) = ', fecritn(1),fecritn(2)
         ENDIF
         IF ( istop(3-iend) .eq. 2 ) THEN
           IF ( my_rank == 0 ) THEN
           write(iunit,*) 'Something screwy! Both ends at ground! ',      &
     &      'Stay with the first and stop now'
           ENDIF
          GOTO 201
         ENDIF
         IF ( iextra .eq. 0 ) iextra = 1
         icgyn = 1 ! Nint(signchan(iend))
         icgpn = Int( 1.1*signchan(iend) )
         istop(iend) = 2
         ihold = 0
         iendgnd = iend
         fecritn(3-iend) = Min( fecrit, fecritn(3-iend) )
! draw channel down to ground
          
         phi1 = phichan(numchan(iend),iend)
          
         DO n=k1-1,1,-1
           numchan(iend) = numchan(iend) + 1

           linksrc(numchan(iend),iend) = numchan(iend)-1
           linkstep(numchan(iend),iend) = nstep

           ibonds(numchan(iend),1,1,iend) = i1
           ibonds(numchan(iend),2,1,iend) = j1
           ibonds(numchan(iend),3,1,iend) =  n + 1
           ibonds(numchan(iend),1,2,iend) = i1
           ibonds(numchan(iend),2,2,iend) = j1
           ibonds(numchan(iend),3,2,iend) =  n
           chanijk(numchan(iend),1,iend) = i1
           chanijk(numchan(iend),2,iend) = j1
           chanijk(numchan(iend),3,iend) = n
           chandist(numchan(iend),iend) = gd(0,0,1,n,1) +      &
     &                       chandist(numchan(iend)-1,iend) 

           linkrank(numchan(iend),iend) = ichoosenextrank

           linkoffset(numchan(iend),1,iend) = ioffset
           linkoffset(numchan(iend),2,iend) = joffset
           
           IF ( ibal .eq. 2 ) THEN
            phichan(numchan(iend),iend) =       &
     &            phichan(numchan(iend)-1,iend) -      &
     &            gd(0,0,-1,n,1)*signchan(iend)*eintz(n)
           
           ELSE
           phichan(numchan(iend),iend) =      &
     &        phichan(numchan(iend)-(k1-n),iend)*Real(n)/Real(k1)
           ENDIF

           IF ( my_rank == ichoosenextrank ) THEN
#if DEBUG > 0
           IF ( ndebug .ge. 1 ) write(iunit,*) 'i1,j1 = ',i1,j1
#endif
           mask(i1,j1,n) = phichan(numchan(iend),iend)
           phip(i1,j1,n) = phichan(numchan(iend),iend)

!           IF ( usebox ) THEN
            rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,n)+cofx(i1,2)*phip(i1+1,j1,n) +      &
     &      cofy(j1,1)*phip(i1,j1-1,n)+cofy(j1,2)*phip(i1,j1+1,n) +      &
     &      cofz(n,1)*phip(i1,j1,n-1)+cofz(n,2)*phip(i1,j1,n+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(n,3))*phip(i1,j1,n)  )
          
             rhs(i1,j1,n) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(n,3))*gdx(i1,1)*gdy(j1,1)
          
!         ENDIF


           ELSE ! check if point is in ghost zone
             IF (       &
     &        (i1 + ioffset .ge. igslg0 .and. i1 + ioffset .le. igslg0 + NLx+1 .and.      &
     &         j1 + joffset .ge. jgslg0 .and. j1 + joffset .le. jgslg0 + NLy+1 ) ) THEN
#if DEBUG > 0
          IF ( ndebug .ge. 1 ) write(iunit,*) 'i1 + ioffset - igslg0, j1 + joffset - jgslg0,',      &
     &     i1 + ioffset - igslg0, j1 + joffset - jgslg0
#endif
              mask(i1 + ioffset - igslg0, j1 + joffset - jgslg0, n) = phichan(numchan(iend),iend)
              phip(i1 + ioffset - igslg0, j1 + joffset - jgslg0, n) = phichan(numchan(iend),iend)

! do not need to set pickmask here because it is a CG and will not grow anymore, anyway
!              pickmask(i1 + ioffset - igslg0, j1 + joffset - jgslg0, n) = IBSET( pickmask(i1 + ioffset - igslg0, j1 + joffset - jgslg0, n), iend)
              
        
             ENDIF
           ENDIF
           
           chanlen = chandist(numchan(iend),iend)
           
         END DO
         kmn = 1
         kmn1 = 1

!
! for ibal=2, reset the channel potential values to zero and relax
!
         IF ( ibal .eq. 2 ) THEN 
         
         IF ( cgfr .eq. 1.0 ) THEN ! set whole flash structure ground potential (zero)
         DO iend0 =1,nend
          DO n=1,numchan(iend0)
          
           IF ( linkrank(n,iend0) == my_rank ) THEN

           i1 = ibonds(n,1,2,iend0)
           j1 = ibonds(n,2,2,iend0)
           k1 = ibonds(n,3,2,iend0)
!           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend0) = 0.0
            phip(i1,j1,k1) = phichan(n,iend0)
            
           ENDIF
           
!           ELSE
!            phichan(n,iend0) = -eint*signchan(iend0)*
!     :                    (chanlen+chandist(numchan(iend0),iend0))
!            phip(i1,j1,k1) = phichan(n,iend0)
!           ENDIF
          END DO ! n
         END DO ! iend0
         
         ELSEIF ( cgfr .gt. 1.0 ) THEN ! test code to set channel potential with an internal field

         iextra = 5 ! prevent extra growth
         delta = phichan(numchan(iend),iend)
         DO iend0 =1,nend
          DO n=1,numchan(iend0)
           i1 = ibonds(n,1,2,iend0)
           j1 = ibonds(n,2,2,iend0)
           k1 = ibonds(n,3,2,iend0)
!           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend0) = phichan(n,iend0) - delta
            phip(i1,j1,k1) = phichan(n,iend0)
!           ELSE
!            phichan(n,iend0) = -eint*signchan(iend0)*
!     :                    (chanlen+chandist(numchan(iend0),iend0))
!            phip(i1,j1,k1) = phichan(n,iend0)
!           ENDIF
          END DO ! n
         END DO ! iend0
 
          IF ( .false. ) THEN
! first set ground channel values to a "missing" value
!         DO iend0 =1,nend
          DO n=1,numchan(iend)
!           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend) = 1.e30
          END DO ! n
!         END DO ! iend0
! now go back along channel to initiation point
!         DO iend0 =1,nend
          n = numchan(iend)
          phichan(n,iend) = 0.0
          
          DO m=1,numchan(iend)
            i = linksrc(n,iend)
!            IF ( i .ge. 1 ) THEN
            x = chandist(n,iend) - chandist(i,iend)
            phichan(i,iend) = phichan(n,iend)
!            ELSEIF ( i .eq. 0 ) THEN
!              EXIT  ! then we have reached the init point
!            ENDIF
            IF ( i .eq. 1 ) EXIT
            n = i
          END DO ! m
!         END DO ! iend0
         
         ENDIF ! .false.
         
         ENDIF ! cgfr
         
!         IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
!         itype = 1
!         IF (ndebug .ge. 1) write(iunit,*) 'strzsor 2'
         IF ( ny .gt. 2 ) THEN

#ifdef MPI
         CALL MPI_Barrier(my_comm, mpi_error_code)


         dt1 = MPI_Wtime()
!         CALL BMG3_SymStd_SOLVE_boxmg(                                   
         IF ( .true. ) THEN

        BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 4
         nboxmgiter = nboxmgiter + 4


         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :             my_comm                                              &
     &             my_comm,mask                                              &
     &             )

        BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 1
!         dt2 = MPI_Wtime()
         dt2 = dt2 + ( MPI_Wtime() - dt1 )

!            WRITE(iunit,*) ' time for boxmg lightning = ', dT2 - dT1

         ELSE
         imn2 = Min(imn,Max(1,imn-2*iadd))
         jmn2 = Min(jmn,Max(1,jmn-2*iadd))
         imx2 = Max(imx,Min(nxs+nbw,imx+2*iadd))
         jmx2 = Max(jmx,Min(nys+nbs,jmx+2*iadd))
         kmx2 = nzs-1
         
!         IF ( ndebug >= 1 ) write(iunit,*) 'call strzsormpi, my_rank = ',my_rank
         
         CALL STRZSORmpi(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,      &
     &      cofx,cofy,cofz, bcx, bcy, ibcn,ibcs,ibce,ibcw,      &
     &      imn2,jmn2,kmn,imx2,jmx2,kmx2,itype,rjac,omega,ibc,soreps0,0,numits)

         nSTRZSORmpi = nSTRZSORmpi + 1
         dtSTRZSORmpi = dtSTRZSORmpi + ( MPI_Wtime() - dt1 )
         numitstot = numitstot + numits
         
         ENDIF
         
         CALL MPI_Barrier(my_comm, mpi_error_code)
#endif
         ELSE
          write(iunit,*) '2D lightning not supported in MPI!!'
          call commasmpi_abort()
         ENDIF
!         GOTO 201
! turn this off for MPI (2015/02/23)
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF
         ENDIF ! ( ibal .eq. 2 )
       
       END IF ! ibal checks
       ENDIF ! cgthres
       
       IF ( istop(iend) == 2 ) EXIT
       
       ENDDO pickloop1

       
       

#ifdef MPI       
! this is my current stopping point
!         lgtstp = 1
!         return

!
! Notes to self:
!
! need to check again if new channel point is in a ghost zone to set mask array
!
         
#endif
       
       IF ( my_rank == ichoosenextrank ) THEN
       
!       IF (j1 .le. nbs+1 .or. ( ny .gt. 2 .and. j1 .ge. nys+nbs+1 ) .or.
!     :     i1 .le. nbw+1 .or. i1 .ge. nxs+nbw+1 ) THEN
!         write(iunit,*) 'reached the side! stopping now..... '
!         GO TO 201
!       END IF
       
       ENDIF

!       mask(i1,j1,k1) = 0.0
!       phip(i1,j1,k1) = phichan(numchan(iend),iend)
      
      i0 = Min(i1,i2)
      j0 = Min(j1,j2)
      k0 = Min(k1,k2)

      imn1 = Min(imn1,i1+ioffset)
      jmn1 = Min(jmn1,j1+joffset)
      kmn1 = Min(kmn1,k1)
      imx1 = Max(imx1,i1+ioffset)
      jmx1 = Max(jmx1,j1+joffset)
      kmx1 = Max(kmx1,k1)

!      iadd = 6
      iadd = NInt(4000./ds)
      iaddz = NInt(3000./dz)

      imn = Min(imn,Max(1,imn1-iadd))
      jmn = Min(jmn,Max(1,jmn1-iadd))
      kmn = Min(kmn,Max(1,kmn1-iaddz))
      imx = Max(imx,Min(nxs+nbw,imx1+iadd))
      jmx = Max(jmx,Min(nys+nbs,jmx1+iadd))
      kmx = Max(kmx,Min(nzs,kmx1+iaddz))
      
      IF ( ny .eq. 2 ) THEN ; jmn = 1 ; jmx = 1 ; ENDIF

#ifndef MPI       
       IF ( bcx .eq. 2 ) THEN
         IF ( imn .eq. 1   ) imx = nxs
         IF ( imx .eq. nxs ) imn = 1
       ENDIF

       IF ( bcy .eq. 2 ) THEN
         IF ( jmn .eq. 1   ) jmx = nys
         IF ( jmx .eq. nys ) jmn = 1
       ENDIF
#endif
      
      iloc0(iend) = i0
      jloc0(iend) = j0
      kloc0(iend) = k0
      
      

      END IF ! (istop(iend) .eq. 0) 
      END IF ! (istop(iend) .le. 1)
      END DO ! iend 

!  relax around each end in random order
! MPI: relax full grid

!      i1 = 1 + Int( 0.5 + ran0mpi(iseed) )
!      j1 = 3 - i1
!      k1 = j1 - i1
      
!      DO iend=i1,j1,k1
!       IF (istop(iend) .eq. 0) THEN 
!        IF (ndebug .ge. 2) write(iunit,*) 'call sor' 
!        itype = 0     
!        IF (ndebug .ge. 1) THEN
!         write(iunit,*) 'strzsor 3, iend=',iend,'i,j,k =',
!     :     iloc0(iend),jloc0(iend),kloc0(iend),rjac,omega
!        ENDIF


        ! done with picking new bonds, so reset mask to 0
        pickmask(:,:,:) = 0
        
        IF ( ny .gt. 2 ) THEN
#ifdef MPI
         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()
         
!         IF ( usebox ) THEN

!      IF ((nstep/intsormpi)*intsormpi .eq. nstep) THEN
!        soreps0 = soreps1
!      ELSE
!        soreps0 = soreps
!      ENDIF


!      IF ( ndebug >= 1 ) write(iunit,*) 'nstep,intfgmpi,my_rank = ',nstep,intfgmpi,my_rank
      IF ((nstep/intfgmpi)*intfgmpi .eq. nstep) THEN

!         IF ( ndebug >= 1 ) write(iunit,*) 'call boxmg 1, my_rank = ',my_rank

       i = BMG_iPARMSlg(id_BMG3_MAX_ITERSdg)
       BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = 1 ! Max( i, 2 )
         nboxmgiter = nboxmgiter + 1

!         CALL BMG3_SymStd_SOLVE_boxmg(                                   
         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :             my_comm                                              &
     &             my_comm,mask                                              &
     &             )


         BMG_iPARMSlg(id_BMG3_MAX_ITERSdg) = i

!         dt2 = MPI_Wtime()
         dt2 = dt2 + ( MPI_Wtime() - dt1 )
         
         ELSE
         itype = 1
!         IF ( ndebug >= 1 ) write(iunit,*) 'nstep,intsormpi,my_rank = ',nstep,intsormpi,my_rank
         IF ( (nstep/intsormpi)*intsormpi == nstep ) THEN
!         IF ( ndebug >= 1 ) write(iunit,*) 'call strzsormpi 1, my_rank = ',my_rank
        CALL STRZSORmpi(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,      &
     &       i0,j0,k0,mdel,      &
     &       cofx,cofy,cofz, bcx, bcy,ibcn,ibcs,ibce,ibcw,      &
     &       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps0,ndebug,numits)

         nSTRZSORmpi = nSTRZSORmpi + 1
         dtSTRZSORmpi = dtSTRZSORmpi + ( MPI_Wtime() - dt1 )
         numitstot = numitstot + numits
          ENDIF
         ENDIF
         
         CALL MPI_Barrier(my_comm, mpi_error_code)

!            WRITE(iunit,*) ' time for boxmg lightning = ', dT2 - dT1
#endif
        ELSE
          write(iunit,*) '2d lightning not supported'
        ENDIF
! turn this off for MPI (2015/02/23)
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF
#if DEBUG > 0
        IF (ndebug .ge. 2) write(iunit,*) 'done sor'
#endif
!       END IF
!      END DO

!      IF ((nstep/intfg)*intfg .eq. nstep .or.
!     :       nstep .eq. nstop )
!     :   THEN
!       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
!      itype = 1
!      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 4'
!      IF ( ny .gt. 2 ) THEN
!      CALL strzsor(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
!     :       cofx,cofy,cofz, bcx, bcy,
!     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
!      ELSE
!      CALL strzsor(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
!     :       cofx,cofy,cofz, bcx, bcy,
!     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
!      ENDIF

!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF
       
!      END IF


!  every so often check the channel charges for imbalance.  If one end has
!  > 50% more charge than the other end, then we adjust the propagation
!  thresholds.

      IF ( ibal .eq. 1 .or. ibal .ge. 3) THEN
      IF ( (( nstep/3*3 .eq. nstep .or. nstep .le. 30 )  .and.      &
     &     ( istop(1) .le. 1 .and. istop(2) .le. 1 ))       &
     &       .or. ihold .ne. 0 ) THEN
      
         IF ( ihold .eq. 1 .and. istop(2) .eq. 1 ) THEN ! holding end 1, but end 2 can't grow
           IF ( fecritn(2) .gt. fecmin ) THEN
            IF ( fecritn(1) .gt. 4.0*fecrit ) THEN
              write(iunit,*) 'stoppping: fecritn(1) = ', fecritn(1)
              GOTO 201
            ENDIF
             IF ( ibal .eq. 1 ) THEN
               fecritn(1) = fecritn(1) + fecinc
             ELSEIF (ibal .eq. 5 ) THEN
               fecritn(1) = Min(fecrit, fecritn(1) + fecinc)
             ENDIF
            fecritn(2) = Max ( fecritn(2) - fecinc, fecmin )
            IF (istop(2) .eq. 1 ) istop(2) = 0
            GOTO 200
           ELSE  ! cant ramp down any further, so drop out
             IF ( (numchan(1) + numchan(2) ) .le. 4 .and.       &
     &             numtry .lt. maxtry  .and. ninittot .ge. 1 ) THEN
               GOTO 199
             ELSE
               IF ( ibal .eq. 1 .or. ibal .eq. 5 ) THEN
                 write(iunit,'(a)') 'Cant force growth. Stop now.'
                 GOTO 201
               ELSEIF ( ibal .ge. 3 ) THEN
!                write(iunit,*) 'nstep,istop2,ihold',nstep,istop(2),ihold
                IF (istop(2) .eq. 1 ) istop(2) = 0
                ihold = 0
               ENDIF
             ENDIF
           ENDIF
         ELSEIF ( ihold .eq. 2 .and. istop(1) .eq. 1 ) THEN
           IF ( fecritn(1) .gt. fecmin ) THEN
            IF ( fecritn(2) .gt. 4.0*fecrit ) THEN
              write(iunit,*) 'stoppping: fecritn(2) = ', fecritn(2)
              GOTO 201
            ENDIF
            IF ( ibal .eq. 1 ) THEN
               fecritn(2) = fecritn(2) + fecinc
             ELSEIF (ibal .eq. 5 ) THEN
               fecritn(2) = Min(fecrit, fecritn(2) + fecinc)
            ENDIF 
            fecritn(1) = Max ( fecritn(1) - fecinc, fecmin )
            IF (istop(1) .eq. 1 ) istop(1) = 0
            GOTO 200
           ELSE ! cant ramp down any further, so drop out
             IF ( (numchan(1) + numchan(2) ) .le. 4 .and.       &
     &             numtry .lt. maxtry  .and. ninittot .ge. 1 ) THEN
               GOTO 199
             ELSE
               IF ( ibal .eq. 1 .or. ibal .eq. 5 ) THEN
                 write(iunit,'(a)') 'Cant force growth. Stop now.'
                 GOTO 201
               ELSEIF ( ibal .ge. 3 ) THEN
!                 write(iunit,*) 'nstep,istop1,ihold',nstep,istop(1),ihold
                 IF (istop(1) .eq. 1 ) istop(1) = 0
                 ihold = 0
               ENDIF
             ENDIF
           ENDIF
         
         ENDIF

      DO iend =1,nend
        qtot(iend) = 0.0
      END DO ! iend
      DO iend =1,nend
        DO n=1,numchan(iend)
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
     

        chg = rhotemp - rho(i1,j1,k1)
       
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ENDIF

        END DO ! n

      END DO ! iend
       

!       IF ( ibal .eq. 5 ) THEN
!         IF ( Abs(qtot(1)) .gt. 1.1*Abs(qtot(2)) )  THEN
!              ihold = 1
!         ELSEIF ( Abs(qtot(2)) .gt. 1.1*Abs(qtot(1)) ) THEN
!              ihold = 2
!         ELSE
!             ihold = 0
!         END IF
!       ENDIF ! ibal .eq. 5
       
       IF ( ibal .eq. 1 .and. Abs(qtot(1)) .gt. 1.05*Abs(qtot(2))) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            ihold = 1
       ELSEIF ( ibal .ge. 5 .and.       &
     &            Abs(qtot(1)) .gt. 1.1*Abs(qtot(2))) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            ihold = 1
       ELSEIF ( ibal .ge. 3 .and. ibal .lt. 5 .and.      &
     &          Abs(qtot(1)) .gt. 1.05*Abs(qtot(2))) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            fecritn(2) = Max ( fecritn(2) - fecinc, fecmin )
            fecritn(1) = Min ( fecritn(1) + fecinc, fecrit)
            ihold = 0
       ELSEIF ( ibal .eq. 1 .and.       &
     &          Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            ihold = 2
       ELSEIF ( ibal .ge. 5 .and.       &
     &          Abs(qtot(2)) .gt. 1.1*Abs(qtot(1)) ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            ihold = 2
       ELSEIF ( ibal .ge. 3 .and. ibal .lt. 5 .and.      &
     &           Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
#if DEBUG > 0
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
#endif
            fecritn(1) = Max ( fecritn(1) - fecinc, fecmin )
            fecritn(2) = Min ( fecritn(2) + fecinc, fecrit)
            ihold = 0
       ELSEIF ( fecritn(1) .gt. fecritn(2) ) THEN
            fecritn(1) = fecritn(1) - fecinc
            fecritn(1) = Max ( fecritn(1), fecrit )  
            fecritn(2) = fecritn(2) + fecinc
            fecritn(2) = Min ( fecritn(2), fecrit ) 
            ihold = 0
       ELSEIF ( fecritn(2) .gt. fecritn(1) ) THEN
            fecritn(1) = fecritn(1) + fecinc
            fecritn(1) = Min ( fecritn(1), fecrit )  
            fecritn(2) = fecritn(2) - fecinc
            fecritn(2) = Max ( fecritn(2), fecrit )  
            ihold = 0
       END IF

       
       IF ( nstep .ge. 40 .and. ibal .eq. 1 .and.      &
     &      ( Abs(qtot(1)) .gt. 2.0*Abs(qtot(2)) .or.      &
     &        Abs(qtot(2)) .gt. 2.0*Abs(qtot(1)) ) ) THEN
         IF ( my_rank == 0 ) THEN
           write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         ENDIF
         GOTO 201
       END IF
       IF ( nstep .ge. 40 .and. ibal .eq. 3 .and.      &
     &      ( Abs(qtot(1)) .gt. 3.0*Abs(qtot(2)) .or.      &
     &        Abs(qtot(2)) .gt. 3.0*Abs(qtot(1)) ) ) THEN
         IF ( my_rank == 0 ) THEN
           write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         ENDIF
         GOTO 201
       END IF
       IF ( nstep .ge. 40 .and. ibal .ge. 4 .and.       &
     &   ( Abs(qtot(1)) .gt. 10. .or. Abs(qtot(2)) .gt. 10. ) .and.      &
     &   (idown .ne. 0 .and. fecritn(Max(1,idown)) .lt. fecrit ) .and.      &
     &      ( Abs(qtot(1)) .gt. 1.4*Abs(qtot(2)) .or.      &
     &        Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) ) ) THEN
         IF ( my_rank == 0 ) THEN
           write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         ENDIF
         GOTO 201
       END IF
 
 200  CONTINUE
       
      END IF ! (nstep/5*5 .eq. nstep )
      
      ELSEIF ( ibal .eq. 2 ) THEN

!
!  Adjust channel reference potential to achieve approximate charge neutrality
!
      IF ( (( (nstep/3)*3 .eq. nstep .or. nstep .le. 40 )  .and.      &
     &     ( istop(1) .le. 1 .and. istop(2) .le. 1 ))       &
     &       .or. ihold .ne. 0 .or. nohold .ge. 1 ) THEN
!
! nohold:  >= 1 : adjust ref. potential after every step
!           <= 0 : adjust ref. potential after every 3rd step or when ihold /= 0
!    (if 0, imposes adjustment limit to be within the ambient values of the first two points)
!
!
#if DEBUG > 0
      IF ( ndebug .ge. 1 ) write(iunit,*) 'start of balance, rank =',my_rank
#endif
      
      IF (  nohold .le. 0 .and.      &
     &     .not. ( ihold .eq. 1 .and. istop(2) .eq. 1 ) .and.      &
     &     .not. ( ihold .eq. 2 .and. istop(1) .eq. 1 )  ) THEN
!      

#if DEBUG > 0
      IF ( ndebug .ge. 1 ) write(iunit,*) 'start of balance0, rank = ',my_rank
#endif
        DO iend =1,nend
          qtot(iend) = 0.0
        END DO ! iend
        DO iend =1,nend
         DO n=1,numchan(iend)
         
         IF ( linkrank(n,iend) == my_rank ) THEN
         
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

!         IF ( usebox ) THEN
          
          rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdx(i1,1)*gdy(j1,1)
          
!         ENDIF

        chg = rhotemp - rho(i1,j1,k1)
       
!         IF ( chg  .lt. 0.0 ) THEN
!           qtot(1) = qtot(1) + chg*zmap(k1)
!         ELSE 
!           qtot(2) = qtot(2) + chg*zmap(k1)
!         ENDIF
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ENDIF
         
        ENDIF ! linkrank

         END DO ! n

       END DO ! iend
       
! allreduce       
#ifdef MPI
       mpitotindp(1)  = qtot(1)
       mpitotindp(2)  = qtot(2)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       qtot(1) = mpitotoutdp(1)  
       qtot(2) = mpitotoutdp(2)  
#endif
       
       IF ( Abs(qtot(1)) .gt. 1.05*Abs(qtot(2)) )  THEN
            ihold = 1
       ELSEIF ( Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
            ihold = 2
       ELSE
           ihold = 0
       ENDIF
       
#ifdef MPI
#if DEBUG > 0
        IF (ndebug .gt. 0 ) write(iunit,*) 'balanc0: myrank,qtot1,2,ihold = ',my_rank,qtot(1),qtot(2),ihold
        IF (ndebug .gt. 0 ) write(iunit,*) 'balanc0: myrank,istop1,2 = ',my_rank,istop(1),istop(2)
#endif
#endif

!       write(iunit,*) 'balck01: nstep = ',nstep,ihold,
!     :   qtot(1)*dls3,qtot(2)*dls3

            IF ( ((Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) .or.      &
     &           Abs(qtot(1)) .gt. 1.4*Abs(qtot(2))  ) .and.      &
     &           nstep .ge. 40*maxpicks) .or.      &
     &   (Abs(qtot(2)) .gt. 120.0 .and. Abs(qtot(1)) .gt. 120.0      &
     &          .and. (numchan(1) .gt. 2*numchan(2)*maxpicks .or.      &
     &          numchan(2) .gt. 2*numchan(1)*maxpicks ) ) ) THEN
              IF ( my_rank == 0 ) THEN
               write(iunit,*) 'Flash is too unbalanced, stopping 1'
               write(iunit,*) 'q1,q2= ',qtot(1),qtot(2) !,qtot(1),qtot(2)
               write(iunit,*) 'numchan1,2 = ',numchan(1),numchan(2)
              ENDIF
              GOTO 201
            ENDIF
            

       
!       
       ELSE ! (one end needs to grow but cant, so adjust potential to balance charge)
            ! OR nohold >= 1 and always adjust reference to balance
!       
!       maxiter = Max(1,nohold)
!       DO iter = 1,maxiter
      k = 0  ! counts how many cycles of adjusting phiref and strzsor with qtot1 > qtot2
      j = 0  ! counts how many cycles of adjusting phiref and strzsor with qtot2 > qtot1
      ilast = 0
      iflip = 0
      istop1 = 0
      bfac = 0.02
      qtot1 = 0.0
      qtot2 = 0.0
         
 21   CONTINUE
 
      IF ( iflip .ge. 5 ) THEN
        IF ( my_rank == 0 ) THEN
          write(iunit,*) 'Too many flips'
        ENDIF
        GOTO 31
      ENDIF
       
      DO iend =1,nend
      qtot(iend) = 0.0
      ENDDO

      DO iend =1,nend
        DO n=1,numchan(iend)
         
         IF ( linkrank(n,iend) == my_rank ) THEN
         
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

!         IF ( usebox ) THEN
          
          rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdx(i1,1)*gdy(j1,1)
          
!         ENDIF

        chg = rhotemp - rho(i1,j1,k1)
       

         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ENDIF

        ENDIF ! linkrank
        END DO ! n
      END DO ! iend

!       write(iunit,*) 'BALCK02: nstep = ',nstep,ihold,
!     :  qtot(1)*dls3,qtot(2)*dls3
       
#ifdef MPI
       mpitotindp(1)  = qtot(1)
       mpitotindp(2)  = qtot(2)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       qtot(1) = mpitotoutdp(1)  
       qtot(2) = mpitotoutdp(2)  
#endif

            IF ( ((Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) .or.      &
     &           Abs(qtot(1)) .gt. 1.4*Abs(qtot(2))  ) .and.      &
     &           nstep .ge. 40*maxpicks) .or.      &
     &   (Abs(qtot(2)) .gt. 60.0 .and. Abs(qtot(1)) .gt. 60.0      &
     &          .and. (numchan(1) .gt. 2*numchan(2)*maxpicks .or.      &
     &          numchan(2) .gt. 2*numchan(1)*maxpicks ) ) ) THEN
              IF ( my_rank == 0 ) THEN
               write(iunit,*) 'Flash is too unbalanced, stopping 2'
               write(iunit,*) 'q1,q2= ',qtot(1),qtot(2)
              ENDIF
              GOTO 201
            ENDIF

#if DEBUG > 0
        IF (ndebug .gt. 0 ) THEN
         write(iunit,*) 'balance1: myrank,ihold = ',my_rank,ihold,k,j
         write(iunit,*) 'balance1: qtot(1,2) = ',qtot(1),qtot(2)
         write(iunit,*) 'balance1: phistart = ',phistart
         write(iunit,*) 'balance1: number of iterations (j,k) = ',j, k
        ENDIF
#endif
!       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
       IF ( istop1 .eq. 1 ) THEN 
         IF ( my_rank == 0 ) THEN
          write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
          write(iunit,*) 'phistart = ',phistart
          write(iunit,*) 'number of iterations (j,k) = ',j, k
         ENDIF
         GOTO 201
       ELSEIF ( istop1 .eq. -1 ) THEN
!         write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
!         write(iunit,*) 'phistart = ',phistart
!         write(iunit,*) 'number of iterations (j,k) = ',j, k
!         write(iunit,*) 'goto 30'
          ihold = 0
          GOTO 31
       ENDIF

            IF ( k .gt. 45 .or. j .gt. 45 ) THEN
             IF ( my_rank == 0 ) THEN
              write(iunit,*) 'Too many iterations (j,k) = ',j,k
              write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
              write(iunit,*) 'phistart,bfac,iflip = ',phistart,bfac,iflip
              write(iunit,*)  'phisthi,phistlow = ',phisthi,phistlow
             ENDIF
            ENDIF

            IF ( k .gt. 50 .or. j .gt. 50 ) THEN
             IF ( my_rank == 0 ) THEN
              write(iunit,*) 'Too many iterations (j,k) = ',j,k
             ENDIF
             GOTO 201
!             call commasmpi_abort()
            ENDIF
!c
       IF ( ( Abs(qtot(1)) .gt. 1.02*Abs(qtot(2)) .and.       &
     &              j.eq.0 .and. k.eq.0) .or.      &
     &        Abs(qtot(1)) .gt. 1.01*Abs(qtot(2) ) ) THEN ! {
!c
#if DEBUG > 0
        IF (ndebug .gt. 0 ) THEN
         write(iunit,*) 'balance2a: myrank,qtot1,2,ihold = ',my_rank,qtot(1),qtot(2),ihold,k,j
        ENDIF
#endif
         IF ( qtot1 .eq. 0.0 )  qtot1 = qtot(1)
         
         IF ( qtot(1) .lt. qtot1 ) THEN  ! {
!c       write(iunit,*) 'qtot(1,2) = ',qtot1,qtot2
            qtot1 = qtot(1)
!c       write(iunit,*) 'nstep = ',nstep
!            write(iunit,*) 'Going the wrong way! nstep =',nstep
!c       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
!       write(iunit,*) 'number of iterations (j,k) = ',j, k
!            istop1 = 1
            istop1 = -1
!            GOTO 21
            
            bfac = -Abs(bfac)
         ELSE
            qtot1 = qtot(1)
         ENDIF  ! }

         IF ( ilast .eq. 2 ) THEN
           iflip = iflip + 1
           bfac = 0.5*bfac
         ENDIF
         ilast = 1
         
! phi(i1,j1,k1),phi(i2,j2,k2),phi(i3,j3,k3)
         IF ( j .gt. 0 .and. k .eq. 0 .and. iflip .gt. 0 ) THEN
!           bfac = 0.5*bfac
         ELSEIF ( iflip .gt. 0 ) THEN
!           bfac = 0.5*bfac
         ENDIF

         
         IF ( istop1 .le. 0 ) THEN
         phistart = phistart + bfac*(phidiff)
         adjlast = bfac*phidiff
         ELSE
         phistart = phistart - adjlast
         IF ( my_rank == 0 ) THEN
           write(iunit,*) 'phistart = ',phistart
         ENDIF
         ENDIF
         
! 5.6.2004 ERM [
         IF (  nohold .eq. 0 .and.       &
     &         ( phistart .gt. phisthi .or. phistart .lt. phistlow) ) THEN
!           phistart =  Min( phisthi, Max( phistlow, phistart ) )
           IF ( my_rank == 0 ) THEN
           write(iunit,*) 'Reached limit of phistart adjustment:',      &
     &               ' phistartn =',phistart
           ENDIF
           GOTO 201
         ENDIF
! 5.6.2004 ERM ]
         
         DO iend =1,nend
          DO n=1,numchan(iend)
           i1 = ibonds(n,1,2,iend)
           j1 = ibonds(n,2,2,iend)
           k1 = ibonds(n,3,2,iend)
            x = phichan(n,iend)
            phichan(n,iend) = phistart -       &
     &         ( chandist(n,iend) - chandist(1,iend))*      &
     &           signchan(iend)*eintz(k1)
            IF ( linkrank(n,iend) == my_rank ) THEN
              phip(i1,j1,k1) = phichan(n,iend)
              mask(i1,j1,k1) = phichan(n,iend)
!             IF ( usebox ) THEN
              rhotemp = -eperao*(      &
     &        cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &        cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &        cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &       (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          
             rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdx(i1,1)*gdy(j1,1)
          
!             ENDIF

            ENDIF ! linkrank
!            write(iunit,*) 'nst,iend,n,old,new ',
!     :       nstep,iend,n,x,phichan(n,iend)

          END DO ! n
         END DO ! iend
            k = k + 1

#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
#endif
      itype = 1
#if DEBUG > 0
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 5'
#endif
      IF ( ny .gt. 2 ) THEN
#ifdef MPI

         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()
          
!          IF ( usebox ) THEN
          IF ( .false. ) THEN
          
          
!         CALL BMG3_SymStd_SOLVE_boxmg(                                   
         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :             my_comm                                              &
     &             my_comm,mask                                              &
     &             )

!         dt2 = MPI_Wtime()
         dt2 = dt2 + ( MPI_Wtime() - dt1 )
         
         nboxmgiter = nboxmgiter + 1

        ELSE
         itype = 1
        
!         IF ( ndebug >= 1 ) write(iunit,*) 'call strzsormpi 2, my_rank = ',my_rank
         CALL STRZSORmpi(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,      &
     &       cofx,cofy,cofz, bcx, bcy,ibcn,ibcs,ibce,ibcw,      &
     &       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps0,0,numits)

         nSTRZSORmpi = nSTRZSORmpi + 1
         dtSTRZSORmpi = dtSTRZSORmpi + ( MPI_Wtime() - dt1 )
         numitstot = numitstot + numits

        ENDIF
        
        CALL MPI_Barrier(my_comm, mpi_error_code)

!            WRITE(iunit,*) ' time for boxmg lightning = ', dT2 - dT1

!        ENDIF
#endif
      ELSE
       write(iunit,*) '2D lightning not supported in MPI'
      ENDIF
! turn this off for MPI (2015/02/23)
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF
      

            GOTO 21
!c
       ELSEIF ( ( Abs(qtot(2)) .gt. 1.02*Abs(qtot(1)) .and.       &
     &              j.eq.0 .and. k.eq.0) .or.      &
     &        Abs(qtot(2)) .gt. 1.01*Abs(qtot(1) ) ) THEN
!c
#if DEBUG > 0
        IF (ndebug .gt. 0 ) THEN
         write(iunit,*) 'balance2b: myrank,qtot1,2,ihold = ',my_rank,qtot(1),qtot(2),ihold,k,j
        ENDIF
#endif
         
         IF ( qtot2 .eq. 0.0 )  qtot2 = qtot(2)
         
         IF ( qtot(2) .gt. qtot2 ) THEN
!c       write(iunit,*) 'qtot(1,2) = ',qtot1,qtot2
            qtot2 = qtot(2)
!c       write(iunit,*) 'nstep = ',nstep
!            write(iunit,*) 'Going the wrong way! nstep =',nstep
!c       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
!       write(iunit,*) 'number of iterations (j,k) = ',j, k
!            istop1 = 1
            istop1 = -1
!            GOTO 201
            bfac = -Abs(bfac)
         ELSE
            qtot2 = qtot(2)
         ENDIF

         IF ( ilast .eq. 1 ) THEN
            iflip = iflip + 1
            bfac = 0.5*bfac
         ENDIF
         
         ilast = 2

         IF ( k .gt. 0 .and. j .eq. 0 .and. iflip .gt. 0 ) THEN
!           bfac = 0.5*bfac
         ELSEIF ( iflip .gt. 0 ) THEN
!           bfac = 0.5*bfac
         ENDIF


         IF ( istop1 .le. 0 ) THEN
         phistart = phistart - bfac*(phidiff)
         adjlast = -bfac*phidiff
         ELSE
         phistart = phistart - adjlast
         IF ( my_rank == 0 ) THEN
           write(iunit,*) 'phistart = ',phistart
         ENDIF
         ENDIF
         
         DO iend =1,nend
          DO n=1,numchan(iend)
           i1 = ibonds(n,1,2,iend)
           j1 = ibonds(n,2,2,iend)
           k1 = ibonds(n,3,2,iend)
            x = phichan(n,iend)
            phichan(n,iend) = phistart -       &
     &         ( chandist(n,iend) - chandist(1,iend))*      &
     &           signchan(iend)*eintz(k1)
            IF ( linkrank(n,iend) == my_rank ) THEN
              phip(i1,j1,k1) = phichan(n,iend)
              mask(i1,j1,k1) = phichan(n,iend)
!             IF ( usebox ) THEN
              rhotemp = -eperao*(      &
     &        cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &        cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &        cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &       (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          
             rhs(i1,j1,k1) = rhotemp/eperao*dx*dy*(dlzlgt/z1d2lgt(k1,3))*gdx(i1,1)*gdy(j1,1)
          
!             ENDIF
            ENDIF ! linkrank
!            write(iunit,*) 'iend,n,old,new ',iend,n,x,phichan(n,iend)
          END DO ! n
         END DO ! iend
            j = j + 1

#if DEBUG > 0
       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
#endif
      itype = 1
#if DEBUG > 0
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 6'
#endif
      IF ( ny .gt. 2 ) THEN
#ifdef MPI


         CALL MPI_Barrier(my_comm, mpi_error_code)

         dt1 = MPI_Wtime()
         
!         IF ( usebox ) THEN
          IF ( .false. ) THEN

!         CALL BMG3_SymStd_SOLVE_boxmg(                                   
         CALL commas_SymStd_SOLVE_boxmg(                                         &
     &             NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg, iGslg, jGslg, kGslg,      &
     &             BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,                         &
     &             Qlg, QFlg, BMG_rWORKlg(BMG_pWORKlg(ip_RESdg)), NFlg, NCbmglg,       &
     &             SOlg, NSOlg,                                                    &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_SORdg)), NSORlg,                         &
     &             BMG_rWORKlg(BMG_pWORKlg(ip_CIdg)), NCIlg,                           &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_iGdg)), NOGlg, NOGclg,                     &
     &             BMG_iWORK_PLlg, NBMG_iWORK_PLlg,                                &
     &             BMG_rWORK_PLlg, NBMG_rWORK_PLlg,                                &
     &             BMG_iWORK_CSlg, NBMG_iWORK_CSlg,                                &
     &             BMG_rWORK_CSlg, NBMG_rWORK_CSlg,                                &
     &             BMG_iWORKlg(BMG_pWORKlg(ip_MSGdg)), NMSGilg,                        &
     &             pMSGlg, pMSGSOlg,                                               &
     &             BMG_MSG_iGRIDlg, NBMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,               &
     &             number_of_processes, BMG_rWORKlg(BMG_pWORKlg(ip_MSG_BUFdg)), NMSGrlg, &
!     :             my_comm                                              &
     &             my_comm,mask                                              &
     &             )

         dt2 = dt2 + ( MPI_Wtime() - dt1 )

         nboxmgiter = nboxmgiter + 1

        ELSE
         itype = 1
        
!         IF ( ndebug >= 1 ) write(iunit,*) 'call strzsormpi 3, my_rank = ',my_rank
         CALL STRZSORmpi(mask,rhs0sor,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,      &
     &       cofx,cofy,cofz, bcx, bcy,ibcn,ibcs,ibce,ibcw,      &
     &       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps0,0,numits)

         nSTRZSORmpi = nSTRZSORmpi + 1
         dtSTRZSORmpi = dtSTRZSORmpi + ( MPI_Wtime() - dt1 )
         numitstot = numitstot + numits

        ENDIF

         CALL MPI_Barrier(my_comm, mpi_error_code)
!         dt2 = MPI_Wtime()

!            WRITE(iunit,*) ' time for boxmg lightning = ', dT2 - dT1
!       ENDIF
#endif
      ELSE
       write(iunit,*) '2D lightning not supported in MPI'
      ENDIF

! turn this off for MPI (2015/02/23)
!      IF ( bcy .eq. 2 ) THEN
!       DO n = 1,2
!         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
!         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
!       ENDDO
!      ENDIF

            GOTO 21
!c
       ELSE
!c
#if DEBUG > 0
        IF ( ndebug .ge. 1 ) THEN
         write(iunit,*) 'End of balance qtot(1,2) = ',qtot(1),qtot(2)
         write(iunit,*) 'phistart = ',phistart
         write(iunit,*) 'nstep, number of iterations (j,k) = ',nstep,j, k
        ENDIF
#endif
         ihold = 0
!c
       END IF
!c       
 31   CONTINUE
!       
       ENDIF  ! ( .not. ( ihold ... ) )
!       

      END IF ! (nstep/5*5 .eq. nstep )
      
      ENDIF ! (ibal)
      
      
#if DEBUG > 0
      IF (ndebug .ge. 2) write(iunit,'(a,i5)') 'nstep = ',nstep
#endif

      nstep = nstep + 1 + Maxval( numextrapicks )
      ! numextrapicks

      ENDDO  LightLoop ! nstep 
      
      IF ( my_rank == 0 ) THEN
        write(iunit,'(a)') 'Did not reach boundary'
      ENDIF
 201  CONTINUE
 
      IF ( my_rank == 0 ) THEN
        write(iunit,'(a,i5)') 'Number of lightning steps: nstep = ',nstep
      ENDIF

      tim2 = commas_dtime(timarr)
!      IF ( lsgi ) tim2 = commas_dtime(timarr)
!      IF ( .not. lsgi ) tim2 = second()
      IF ( my_rank == 0 ) THEN
       write(iunit,*) 'time for propagation = ',tim2,timarr(1)
       write(iunit,*) 'propagation time for boxmg = ', dt2, nboxmgiter
       write(iunit,'(a,f13.5,2(2x,i6))') 'propagation time for STRZSORmpi = ', dtSTRZSORmpi,nSTRZSORmpi,numitstot
!      IF ( .not. lsgi ) write(iunit,*) 'time for propagation = ',tim2-tim1
!      write(0,'(a,2(1x,1pe12.5))') ' time for prop = ',tim2,timarr(1)
      write(iunit,'(a,a,3(i4,1x),4x,3(i4,1x))') 'Flash bounding box ',      &
     &    'imn,jmn,kmn,imx,jmx,kmx: ',imn1,jmn1,kmn1,imx1,jmx1,kmx1
      write(iunit,'(a)') 'End of lightning propagation'
      write(iunit,'(a,2(1x,i5))') 'numchan(1),numchan(2)',      &
     &         numchan(1),numchan(2)
      ENDIF

!
!  Check energy
!
      eold = 0.0
      enew = 0.0
      eoldct = 0.0
      enewct = 0.0
      IF ( icgyn .ge. 1 ) THEN
       kmn = 1
       kmx = nzs-1
       iadd = 10
       imn = Min(imn,Max(1,imn-iadd))
       jmn = Min(jmn,Max(1,jmn-iadd))
       imx = Max(imx,Min(nxs+nbw,imx+iadd))
       jmx = Max(jmx,Min(nys+nbs,jmx+iadd))
      ENDIF

!      jmn2 = 2
!      jmx2 = 1
      

      ist = Max(1,nbw+1)
      jst = Max(1,nbs+1)
      jend = Max(1, nny-nbn-1)

      ixe = NLx
      IF ( myproci == nproci ) ixe = nnx - nbe - 1
      
      jye = NLy
      IF ( myprocj == nprocj ) jye = Max(1, nny - nbn - 1)


      DO k1=1,nzs-1
        vol = dls3*zmap(k1)
       DO j1 = jst,jye ! jmn,jmx
        DO i1 = ist,ixe ! imn,imx
         x = 0.5*phi2(i1,j1,k1)*rho(i1,j1,k1)*vol
         eold = eold + x
        
        IF ( mask(i1,j1,k1) .gt. 0.9e30 ) THEN ! non-channel point
          y = 0.5*phip(i1,j1,k1)*rho(i1,j1,k1)*vol
!          enew = enew + 0.5*phip(i1,j1,k1)*rho(i1,j1,k1)*vol
        ELSE ! channel point; use Gauss Law to get charge density -- no, actually it is the Poisson equation here.
           rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

          y = 0.5*phip(i1,j1,k1)*rhotemp*vol*gdx(i1,1)*gdy(j1,1)
!          enew = enew + 0.5*phip(i1,j1,k1)*rhotemp*vol
        ENDIF
        
          enew = enew + y

!         IF ( ibal.eq.2 .and. icldtop .ne. 0 .and. 
!     :        k1 .gt. icldtop + 4 ) THEN
!           eoldct = eoldct + x
!           enewct = enewct + y
!         ENDIF
        
        
        ENDDO
       ENDDO
      ENDDO

! allreduce
#ifdef MPI
       mpitotindp(1)  = eold
       mpitotindp(2)  = enew

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       eold = mpitotoutdp(1)  
       enew = mpitotoutdp(2)  
#endif
      
      IF ( my_rank == 0 .and. ibal .eq. 2 ) THEN
       write(iunit,'(a,3(1x,1pe13.5))') 'phiref, phistart, delta-phi = ',      &
     &   phistart,phistart0,phistart-phistart0
      
      ENDIF
       
       IF ( my_rank == 0 ) THEN
       write(iunit,'(a,3(1x,1pe13.5))')       &
     &    'Energy check (old,new,diff): ',eold,enew,enew-eold
       ENDIF

       lct = .true.

       IF ( my_rank == 0 .and. ( enewct .ne. 0.0 .or. eoldct .ne. 0.0 ) ) THEN
         write(iunit,'(a,3(1x,1pe13.5))')       &
     &  'CT Energy check (old,new,diff): ',eoldct,enewct,enewct-eoldct
        lct = .true.
!        IF ( (eoldct .gt. 0.0 .and. enewct .gt. 0.0 .and.
!     :        enewct .gt. 1.05*eoldct)  .or.
!     :    ( Sign(1.0,enewct) .ne. Sign(1.0,eoldct) ) .or.
!     :    ( Abs(enewct) .gt. 1.1*Abs(eoldct) .and.
!     :      Sign(1.0,enewct) .eq. Sign(1.0,eoldct)  ) ) THEN
!          lct = .false.
!        ENDIF
        
        
       ENDIF
       IF ( .not. (enew .gt. 0.0 .and. enew .lt. eold .and. lct )      &
     &      .and. ibal .eq. 2       &
     &         ) THEN
!     :      .and. bcy .ne. 2 .and. .false. 
!     :      .and. .not. ( cgfr .lt. 1.0 .and. icgyn .ne. 0 ) ) THEN
!
! 2005.10.03: let it pass if cgfr .lt. 1 and is a CG.
!
         ldumpsp = .true.
         IF ( my_rank == 0 ) THEN
           write(iunit,*) 'Failed energy reduction! Try again'
         ENDIF
      DO iend =1,nend
        qtot(iend) = 0.0
        rhotot(iend) = 0.0
        rhotot1(iend) = 0.0
      END DO ! iend
      
      DO iend =1,nend
        DO n=1,numchan(iend)
         
         IF ( linkrank(n,iend) == my_rank ) THEN
         
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)
#ifdef MPI
         i2 = ibonds(n,1,2,iend)-nbw
         j2 = ibonds(n,2,2,iend)-nbs
#else
         i2 = ibonds(n,1,2,iend)-nbw-1
         j2 = ibonds(n,2,2,iend)-nbs-1
#endif
         k2 = ibonds(n,3,2,iend)

         
         rhotemp = 0.0
         chg = 0.0
         IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
        
         chg = rhotemp - rho(i1,j1,k1)
         
         ELSE
          chg = 0.0
          rhotemp = 0.0
         ENDIF

        chg0 = rho(i1,j1,k1)
        
         IF ( (ibal .ge. 5 .or. ibal .eq. 2) .and.       &
     &        0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and.      &
     &        Sign(1.0d0,chg) .ne. signchan(iend) ) THEN
           chg = 0.0
           rhotemp = rho(i1,j1,k1)
         ENDIF
          
!         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)
!         IF ( chg  .lt. 0.0 ) THEN
!           qtot(1) = qtot(1) + chg*zmap(k2)
!           rhotot(1) = rhotot(1) + chg0*zmap(k2)
!         ELSE 
!           qtot(2) = qtot(2) + chg*zmap(k2)
!           rhotot(2) = rhotot(2) + chg0*zmap(k2)
!         ENDIF

!         write(0,*) 'i,j,k,zmap: ',i1,j1,i1+nbw,j1+nbs,k2
!         write(0,*) 'zmap etc. =', zmap(k2),dlzlgt,x1d2lgt(i1),y1d2lgt(j1)
         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
           rhotot(1) = rhotot(1) + chg0*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
           rhotot(2) = rhotot(2) + chg0*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ENDIF
         
         ENDIF ! linkrank
         
        END DO ! n
      END DO ! iend
      
#ifdef MPI
       mpitotindp(1)  = qtot(1)
       mpitotindp(2)  = qtot(2)
       mpitotindp(3)  = rhotot(1)
       mpitotindp(4)  = rhotot(2)
       mpitotindp(5)  = rhotot1(1)
       mpitotindp(6)  = rhotot1(2)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 6, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       qtot(1) = mpitotoutdp(1)
       qtot(2) = mpitotoutdp(2)
       rhotot(1) = mpitotoutdp(3)
       rhotot(2) = mpitotoutdp(4)
       rhotot1(1) = mpitotoutdp(5)
       rhotot1(2) = mpitotoutdp(6)

#endif
      
      DO iend=1,nend
      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,i1,a,1pe13.5)') 'initial charge on channel qtot(',      &
     &             iend, ')= ', qtot(iend) !*dls3
      write(iunit,'(a,1pe13.5)') 'rhotot  = ',rhotot(iend) ! *dls3
      write(iunit,'(a,1pe13.5)') 'rhotot1 = ',rhotot1(iend) ! *dls3
      ENDIF
       chgdis0(iend) = 0.0d0
      END DO ! iend

        IF ( numtry .eq. 1 ) THEN
          IF ( my_rank == 0 ) THEN
            write(iunit,'(a)') 'Try 1'
          ENDIF
        ENDIF
         
         IF ( isa .eq. 1 ) THEN
           call commasmpi_abort()
         ENDIF

         lgtstp = 2
         GOTO 199
       ELSE
         lgtstp = 0
       ENDIF  ! energy check
!
! write text output file
!      
      times = (nstep1-1)*(dtp)
      
      
      IF ( isa .lt. 1 ) THEN
      write(nmliter,915) iliter
 915  format(i2.2)
      write(rstime,'(i6.6)') time
! 912  format(i5.5)
!      open(unit=16,file=hdrcfs(1:lhdrcfs)//'flashes/'//
!     :  filehead(1:lfilehead)
!     :  //'.chnl.'//rstime(1:ihttima)//'.'//nmliter, form='formatted',
!     :  status = 'unknown')
!       rewind(16)
       
      IF ( my_rank == 0 ) THEN

       write(16,*)
       write(16,'(a,f9.2,a,i3)') 'FLASH at time = ', time_real, ', iliter = ',iliter

      ENDIF ! my_rank == 0

      ENDIF ! ( isa .lt. 1 )
      
      IF ( isa .ge. 1 ) THEN
        write(16,'(a,i5)') 'Flash number ',iliter
      ENDIF ! ( isa .eq. 1 )


      IF ( my_rank == 0 ) THEN

      write(16,'(5(1x,i5),3(1x,f8.1),2(1x,i3))')       &
     &   numchan(1),numchan(2),NGxlg,NGylg,nzs,dx,dy,dz,nbw,nbs

      ENDIF ! my_rank == 0

      pchnl = numchan(1)
      nchnl = numchan(2)
!      volbi = 1.0/(dx1*dy1*dz1)
      DO iend =1,nend
      qtot(iend) = 0.0d0
      rhotot(iend) = 0.0
      rhotot1(iend) = 0.0
      chgdis0(iend) = 0.0d0
      chgdis(iend) = 0.0d0
      numchancheck(iend) = 0
      chanlentot(iend) = 0.0
      chanareatot(iend) = 0.0
      footprint(:,:,iend) = 0
!      icnt1(iend) = 0
!      icnt2(iend) = 0
!      icnt3(iend) = 0
!      icnt4(iend) = 0
      ENDDO
       chanareatot(nend+1) = 0
      
       
          chg = 0.0
          rhotemp = 0.0
          rhox = 0.0
          phix = 0.0
      
      
      DO iend =1,nend
      bondinfoin(:,:) = -1.e36 ! initialize values to large negative value so we can use an mpi Max function

        DO n=1,numchan(iend)

         IF ( my_rank == linkrank(n,iend) ) THEN !{
         i1 = ibonds(n,1,2,iend) ! + linkoffset(numchan(iend),1,iend)
         j1 = ibonds(n,2,2,iend) ! + linkoffset(numchan(iend),2,iend)
         k1 = ibonds(n,3,2,iend)
#ifdef MPI
         i2 = ibonds(n,1,2,iend)-nbw
         j2 = ibonds(n,2,2,iend)-nbs
#else
         i2 = ibonds(n,1,2,iend)-nbw-1
         j2 = ibonds(n,2,2,iend)-nbs-1
#endif
         k2 = ibonds(n,3,2,iend)
         numchancheck(iend) = numchancheck(iend) + 1
         
         IF ( footprint(i1,j1,iend) == 0 ) THEN
           footprint(i1,j1,iend) = 1
           chanareatot(iend) = chanareatot(iend) + x1d2lgt(i1)*y1d2lgt(j1) ! + 1
           IF ( iend == 1 ) THEN ! add 1 to total for all on first end
             chanareatot(nend+1) = chanareatot(nend+1) +  x1d2lgt(i1)*y1d2lgt(j1) ! + 1
           ELSE ! add 1 to total if area is covered by the second end but not the first
             IF ( footprint(i1,j1,1) == 0 ) chanareatot(nend+1) = chanareatot(nend+1) +  x1d2lgt(i1)*y1d2lgt(j1) !+ 1
           ENDIF
         ENDIF
         
         IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

         chg = rhotemp - rho(i1,j1,k1)
         ELSE
          rhotemp = 0.0
          chg = 0.0
         ENDIF
         
         rhox = rho(i1,j1,k1)
         phix = phi2(i1,j1,k1) 

#ifdef MPI
          i3 = i1-nbw
          j3 = j1-nbs
#else
          i3 = i1-nbw-1
          j3 = j1-nbs-1
#endif

          x = ds*Real(i3 - 1)
          y = ds*Real(j3 - 1)
          z = hgt(k1)

          IF ( igrid .eq. 0 ) THEN
          i = NInt( (x + eps1)/dx1 ) + 1
          j = NInt( (y + eps1)/dy1 ) + 1
          
          ELSEIF ( igrid .eq. 1 ) THEN
          i = (i3 - 1)/ifac + 1
          j = (j3 - 1)/ifac + 1
          ENDIF

          IF ( istretch .or. dz .eq. dz1 ) THEN
            k   = k1
          ELSE
            k   = Int( ((k1-1)*dz + eps1)/dzl ) + 1
          ENDIF


          energy = t5(i,j,k)

         bondinfoin(1,n) = phix    ! initial potential
         bondinfoin(2,n) = rhox    ! initial charge
         bondinfoin(3,n) = rhotemp ! final charge
         bondinfoin(4,n) = chg     ! channel charge
         bondinfoin(5,n) = energy  ! cell energy
         
         
         ENDIF ! } linkrank
        
        ENDDO

#ifdef MPI
!         k = 5
!         CALL MPI_AllReduce(bondinfoin(1,n), bondinfoout(1,n), k, MPI_real, MPI_Max, my_comm, mpi_error_code)

         k = 5*numchan(iend)
         CALL MPI_AllReduce(bondinfoin(1,1), bondinfoout, k, MPI_real, MPI_Max, my_comm, mpi_error_code)
#else
         bondinfoout(:,:) = bondinfoin(:,:)
#endif



        DO n=1,numchan(iend)
         
         rhotemp = 0.0
         chg = 0.0
         IF ( my_rank == linkrank(n,iend) ) THEN !{
         i1 = ibonds(n,1,2,iend) ! + linkoffset(numchan(iend),1,iend)
         j1 = ibonds(n,2,2,iend) ! + linkoffset(numchan(iend),2,iend)
         k1 = ibonds(n,3,2,iend)
#ifdef MPI
         i2 = ibonds(n,1,2,iend)-nbw
         j2 = ibonds(n,2,2,iend)-nbs
#else
         i2 = ibonds(n,1,2,iend)-nbw-1
         j2 = ibonds(n,2,2,iend)-nbs-1
#endif
         k2 = ibonds(n,3,2,iend)
!         numchancheck(iend) = numchancheck(iend) + 1
         
         
         
         IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN

          rhotemp = -eperao*(      &
     &      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +      &
     &      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +      &
     &      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +      &
     &      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

        
         chg = rhotemp - rho(i1,j1,k1)
         rhox = rho(i1,j1,k1)
         phix = phi2(i1,j1,k1) 
         
         ELSE
          chg = 0.0
          rhotemp = 0.0
          rhox = 0.0
          phix = 0.0
          energy = 0
         ENDIF


        chg0 = rho(i1,j1,k1)
        
!         IF ( (ibal .ge. 5 .or. ibal .eq. 2) .and. 
!     :        0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and.
!     :        Sign(1.0,chg) .ne. signchan(iend) ) THEN
!           chg = 0.0
!           rhotemp = rho(i1,j1,k1)
!         ENDIF
          
!         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)
!!         rhotot1(iend) = rhotot1(iend) + rhow(i1,j1,k1)*zmap(k2)
!         IF ( chg  .lt. 0.0 ) THEN
!           qtot(1) = qtot(1) + chg*zmap(k2)
!           rhotot(1) = rhotot(1) + chg0*zmap(k2)
!         ELSE 
!           qtot(2) = qtot(2) + chg*zmap(k2)
!           rhotot(2) = rhotot(2) + chg0*zmap(k2)
!         ENDIF

!         write(0,*) 'i,j,k,zmap2: ',i1,j1,i1+nbw,j1+nbs,k2
!         write(0,*) 'zmap etc. = ',zmap(k2),dlzlgt,x1d2lgt(i1),y1d2lgt(j1)
         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
           rhotot(1) = rhotot(1) + chg0*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
           rhotot(2) = rhotot(2) + chg0*zmap(k2)*dlzlgt*x1d2lgt(i1)*y1d2lgt(j1)
         ENDIF

         
         ELSE 
          chg = 0.0
          rhotemp = 0.0
          rhox = 0.0
          phix = 0.0
         ENDIF ! } linkrank


         
!         ENDDO ! numchan

         
!         DO n=1,numchan(iend)

!         phix = bondinfoout(1,n) 
!         rhox = bondinfoout(2,n) 
!         rhotemp = bondinfoout(3,n) 
!         chg = bondinfoout(4,n) 
!         energy = bondinfoout(5,n) 

         IF ( n > 1 ) THEN
          chanlen =  chandist(n,iend) - chandist( linksrc(n,iend), iend)
         ELSE
          chanlen =  chandist(n,iend)
         ENDIF
         chanlentot(iend) = chanlentot(iend) + chanlen

         
         IF ( my_rank == 0 ) THEN

           ibonds(n,1,2,iend) = ibonds(n,1,2,iend) + linkoffset(n,1,iend)
           ibonds(n,2,2,iend) = ibonds(n,2,2,iend) + linkoffset(n,2,iend)

           ibonds(n,1,1,iend) = ibonds(n,1,1,iend) + linkoffset(Max(1,linksrc(n,iend)),1,iend)
           ibonds(n,2,1,iend) = ibonds(n,2,1,iend) + linkoffset(Max(1,linksrc(n,iend)),2,iend)

         write(16,'(3(1x,i5),3(2x,i6))') (ibonds(n,i,1,iend), i=1,3),      &
     &         n, linksrc(n,iend),linkstep(n,iend)
         IF ( ibal .lt. 5 .and. ibal .ne. 2 ) THEN
         write(16,'(3(1x,i4),5(1x,1pe13.5))')       &
     &      (ibonds(n,i,2,iend), i=1,3), phix,  & ! phi(i1,j1,k1),      
     &      rhox,rhotemp,                       & ! rho(i1,j1,k1),rhotemp,      &
     &       chandist(n,iend), chg
!     :      (ibonds(n,i,2,iend), i=1,3),phi(i1,j1,k1),
!     :      rho(i1,j1,k1),rhotemp,
!     :       chandist(n,iend),chg
         ELSE
!         write(16,'(3(1x,i4),7(1x,1pe13.5))') 
!     :      (ibonds(n,i,2,iend), i=1,3), phix, ! phi(i1,j1,k1),
!     :      rhox,rhotemp, ! rho(i1,j1,k1),rhotemp,
!     :       chandist(n,iend), chg ,phichan(n,iend),energy

!          bondinfoin(1,n) = phix    ! initial potential
!          bondinfoin(2,n) = rhox    ! initial charge
!          bondinfoin(3,n) = rhotemp ! final charge
!          bondinfoin(4,n) = chg     ! channel charge
!          bondinfoin(5,n) = energy  ! cell energy
        ! bondinfo: 1=init potential, 2=initial charge density, 3=final charge, 4=channel charge, 5=energy (MJ?)
       ! output: channel point indices (3 int); 
       !         initial potential, 
       !         initial charge density
       !         final charge, 
       !         channel distance,
       !         channel charge, 
       !         channel potential
       !         energy
        write(16,'(3(1x,i4),7(1x,1pe13.5))')       &
     &      (ibonds(n,i,2,iend), i=1,3), bondinfoout(1,n),       &
     &      bondinfoout(2,n),bondinfoout(3,n),      &
     &      chandist(n,iend), bondinfoout(4,n) ,phichan(n,iend),      &
     &      bondinfoout(5,n)

! flash2 file
!        write(17,'(3(1x,i4),7(1x,1pe13.5))') 
!     :      (ibonds(n,i,2,iend), i=1,3), bondinfoout(1,n), 
!     :      bondinfoout(2,n),bondinfoout(3,n),
!     :      chandist(n,iend), bondinfoout(4,n) ,phichan(n,iend),
!     :      bondinfoout(5,n)


!     :      (ibonds(n,i,2,iend), i=1,3),phi(i1,j1,k1),
!     :      rho(i1,j1,k1),rhotemp,
!     :       chandist(n,iend),chg,phichan(n,iend)
         
         ENDIF
         
         ENDIF ! my_rank == 0
         

         IF ( my_rank == linkrank(n,iend) ) THEN
         

!         IF ( istretch ) THEN
!          write(iunit,*) 'first i1,j1,k1,nbw,nbs = ',i1,j1,k1,nbw,nbs
!         ENDIF

#ifdef MPI
          i1 = i1-nbw
          j1 = j1-nbs
#else
          i1 = i1-nbw-1
          j1 = j1-nbs-1
#endif
!          x = 0.5*ds*Real(i1 + i2 - 2)
!          y = 0.5*ds*Real(j1 + j2 - 2)
          x = ds*Real(i1 - 1)
          y = ds*Real(j1 - 1)
          z = hgt(k1)
!          i = Int( (x + eps1)/dx1 + 0.5) + 1
!          j = Int( (y + eps1)/dy1 + 0.5) + 1
!          k = Int( (z + eps1)/dz1 + 0.5) + 1
          
          IF ( igrid .eq. 0 ) THEN
          i = NInt( (x + eps1)/dx1 ) + 1
          j = NInt( (y + eps1)/dy1 ) + 1
          
          ELSEIF ( igrid .eq. 1 ) THEN
          i = (i1 - 1)/ifac + 1
          j = (j1 - 1)/ifac + 1
          ENDIF
          IF ( istretch .or. dz .eq. dz1 ) THEN
            k   = k1
            IF ( lightintz .gt. 1 ) k = (k1 + 1)/lightintz
            ds3 = ds2*dzz(k)
            ds3 = ds2*dlzlgt/z1d2lgt(k1,3)
          ELSE
            k   = Int( ((k1-1)*dz + eps1)/dzl ) + 1
!            ds3 = ds2*dlzlgt
            ds3 = ds2*dz
          ENDIF
          

      IF ( .not. ( icgyn .eq. 1 .and. iend .eq. iendgnd )       &
     &          .or. ibal .ge. 5  .or. ibal .eq. 2       &
     &   ) THEN
      
        IF ( (dx .eq. dx1 .and. dy .eq. dy1 .and. dz .eq. dz1) .or.      &
     &        igrid .eq. 1 ) THEN !{
          volbi = gxt(i,3)*gyt(j,3)*gzt(k,3)
!         IF ( istretch ) write(iunit,*) 'volbi,ds3,chg ', volbi,ds3,chg,volbi*ds3
!         IF ( istretch ) THEN
!          write(iunit,*) 'volbi,ds3,chg ', volbi,ds3,chg,volbi*ds3,i,j,k
!          write(iunit,*) 'i1,j1,k1 = ',i1,j1,k1
!          write(iunit,*) 'dxx,dyy,dzz = ',dxx(i),dyy(j),dzz(k),dxx(i)*dyy(j)*dzz(k)
!         ENDIF
         
         ! t3(i,j,k) = t3(i,j,k) + (chg)*volbi*ds3
         t3(i,j,k) = t3(i,j,k) + (chg)*volratio

         IF ( chg .lt. 0.0 ) THEN
           chgdis0(1) = chgdis0(1) + (chg)*volratio/volbi ! ds3 ! *volbi*dxl*dyl*dzl/gt(1,1,k1,imapz)
         ELSE
           chgdis0(2) = chgdis0(2) + (chg)*volratio/volbi ! ds3 ! *volbi*dxl*dyl*dzl/gt(1,1,k1,imapz)
         ENDIF
        ELSE !}
!         IF ( .false. ) THEN
!
! Figure out which grid boxes should recieve charge and what fraction
!
!  First find which coarse point is closest to each corner of the fine grid
!  xw,xe,ys,yn define the edges of the fine grid square
         xw = x - 0.5*ds
         xe = x + 0.5*ds
         ys = y - 0.5*ds
         yn = y + 0.5*ds
! note: zu and zd are only different when using a constant vertical spacing
! and dslight < dz.
         zu = z + 0.5*ds
         zd = z - 0.5*ds
!
! iw,ie,js,jn are the points on the coarse grid closest to the fine grid edges
!         
         iw = NInt( xw/dx1 + 1.0 )
         ie = NInt( xe/dx1 + 1.0 )
         js = NInt( ys/dy1 + 1.0 )
         jn = NInt( yn/dy1 + 1.0 )
         ku = Int( zu/dz1 + 1.0 )
         kd = Int( zd/dz1 + 1.0 )
         
         IF ( js .le. 0 ) THEN
#if DEBUG > 0
           IF ( ndebug .ge. 1 ) THEN
             write(iunit,*) 'js = ',js,y,ys,dy1,ys/dy1
             write(iunit,*) 'j1,nbs = ',j1,nbs,j1+nbs+1
           ENDIF
#endif
         ENDIF
         
         IF ( istretch .or. dz .eq. dz1 ) THEN
          ku = k1
          kd = k1
         ENDIF
          k = kd
         
!
!  Now find out how many coarse grid cells are overlapped (1,2, or 4)
!    (xc,yc,zc) is the intersection point of the coarse grid cells within
!    the fine grid box and is used for determining the area/volume overlaps
!
!          chgdis0(iend) = chgdis0(iend) + (chg)*volbi*ds3*dxl*dyl*dzl/
!     :     gt(1,1,k1,imapz)
      
      IF ( iw .lt. 1 ) THEN
        write(0,*) 'iw,xw,i1 = ',iw,xw,i1
      ENDIF
      
      vol = dxx(iw)*dyy(js)*dzz(k) ! dxl*dyl*dzl/gt(iw,js,k1,imapz)
      volbi = gxt(iw,3)*gyt(js,3)*gzt(k,3)
!      print*, 'iw,js,k,volbi,ds3,chg ', iw,js,k,volbi,ds3,chg

      IF ( chg .lt. 0.0 ) THEN
          chgdis0(1) = chgdis0(1) + (chg)*ds3 ! *volbi*dxl*dyl*dzl/gt(1,1,k1,imapz)
      ELSE
          chgdis0(2) = chgdis0(2) + (chg)*ds3 ! *volbi*dxl*dyl*dzl/gt(1,1,k1,imapz)
      ENDIF
      
      IF ( ku - kd .eq. 0 ) THEN ! (1,2, or 4 cells)
        IF ( ie - iw .eq. 0 .and. jn - js .eq. 0 ) THEN ! one  cell
          t3(ie,js,k) = t3(ie,js,k) + (chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 0 .and. jn - js .eq. 1 ) THEN ! two cells (n-s)
          yc = gyt(js,1) - gyt(1,1) + dy1/2 ! ( Real(js) - 0.5 )*dy1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          IF ( frs.lt.0.0 .or. frs.gt.1.0) THEN
           write(iunit,*) 'nuts! frs,frn,yc,ys,dy,js = ',frs,frn,yc,ys,dy,js
           DO j=1,ny
             write(iunit,*) 'j,gyt = ',j,gyt(j,1)
           ENDDO
           CALL commasmpi_abort()
          ENDIF
          t3(ie,js,k) = t3(ie,js,k) + frs*(chg)*volbi*ds3
          t3(ie,jn,k) = t3(ie,jn,k) + frn*(chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 0 ) THEN ! two cells (e-w)
          xc = gxt(iw,1) - gxt(1,1) + dx1/2 ! ( Real(iw) - 0.5 )*dx1
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          IF ( frw.lt.0.0 .or. frw.gt.1.0) THEN
           write(iunit,*) 'nuts! frw,fre = ',frw,fre
          ENDIF
          t3(iw,js,k) = t3(iw,js,k) + frw*(chg)*volbi*ds3
          t3(ie,js,k) = t3(ie,js,k) + fre*(chg)*volbi*ds3
          
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 1 ) THEN  ! four cells
          yc = gyt(js,1) - gyt(1,1) + dy1/2 ! ( Real(js) - 0.5 )*dy1
          xc = gxt(iw,1) - gxt(1,1) + dx1/2 ! ( Real(iw) - 0.5 )*dx1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          IF ( frs.lt.0.0 .or. frs.gt.1.0) THEN
           write(iunit,*) 'nuts4! frs,frn = ',frs,frn
          ENDIF
          IF ( frw.lt.0.0 .or. frw.gt.1.0) THEN
           write(iunit,*) 'nuts4! frw,fre = ',frw,fre
          ENDIF
          t3(ie,js,k) = t3(ie,js,k) + frs*fre*(chg)*volbi*ds3
          t3(ie,jn,k) = t3(ie,jn,k) + frn*fre*(chg)*volbi*ds3
          t3(iw,js,k) = t3(iw,js,k) + frs*frw*(chg)*volbi*ds3
          t3(iw,jn,k) = t3(iw,jn,k) + frn*frw*(chg)*volbi*ds3
        ENDIF
       ELSEIF ( ku - kd .eq. 1 ) THEN ! (2,4 or 8 cells)
        zc = ( Real(kd) - 0.5 )*dz1
        frd = ( zc - zd )/dz
        fru = 1.0 - frd
        IF ( ie - iw .eq. 0 .and. jn - js .eq. 0 ) THEN ! two  cell
          t3(ie,js,k)   = t3(ie,js,k)   + frd*(chg)*volbi*ds3
          t3(ie,js,k+1) = t3(ie,js,k+1) + fru*(chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 0 .and. jn - js .eq. 1 ) THEN ! four cells (n-s)
          yc = ( Real(js) - 0.5 )*dy1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          t3(ie,js,k) = t3(ie,js,k) + frd*frs*(chg)*volbi*ds3
          t3(ie,jn,k) = t3(ie,jn,k) + frd*frn*(chg)*volbi*ds3
          t3(ie,js,k+1) = t3(ie,js,k+1) + fru*frs*(chg)*volbi*ds3
          t3(ie,jn,k+1) = t3(ie,jn,k+1) + fru*frn*(chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 0 ) THEN ! four cells (e-w)
          xc = ( Real(iw) - 0.5 )*dx1
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          t3(iw,js,k) = t3(iw,js,k) + frd*frw*(chg)*volbi*ds3
          t3(ie,js,k) = t3(ie,js,k) + frd*fre*(chg)*volbi*ds3
          t3(iw,js,k+1) = t3(iw,js,k+1) + fru*frw*(chg)*volbi*ds3
          t3(ie,js,k+1) = t3(ie,js,k+1) + fru*fre*(chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 1 ) THEN  ! eight cells
          yc = ( Real(js) - 0.5 )*dy1
          xc = ( Real(iw) - 0.5 )*dx1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          t3(ie,js,k) = t3(ie,js,k) + frd*frs*fre*(chg)*volbi*ds3
          t3(ie,jn,k) = t3(ie,jn,k) + frd*frn*fre*(chg)*volbi*ds3
          t3(iw,js,k) = t3(iw,js,k) + frd*frs*frw*(chg)*volbi*ds3
          t3(iw,jn,k) = t3(iw,jn,k) + frd*frn*frw*(chg)*volbi*ds3
          t3(ie,js,k+1) = t3(ie,js,k+1) + fru*frs*fre*(chg)*volbi*ds3
          t3(ie,jn,k+1) = t3(ie,jn,k+1) + fru*frn*fre*(chg)*volbi*ds3
          t3(iw,js,k+1) = t3(iw,js,k+1) + fru*frs*frw*(chg)*volbi*ds3
          t3(iw,jn,k+1) = t3(iw,jn,k+1) + fru*frn*frw*(chg)*volbi*ds3
        ENDIF
       
       
       
        ELSE
          write(iunit,*) 'Something wrong! iw,ie,js,jn,kd,ku = ',      &
     &               iw,ie,js,jn,kd,ku
          call commasmpi_abort()
        ENDIF

!         t4(i,j,k)  = -Sign(1.0,t3(i,j,k) )
      ENDIF
      ENDIF
         t2(i,j,k) = t2(i,j,k) + signchan(iend)
         t1(i,j,k) = t1(i,j,k) + chanlen


        ENDIF ! linkrank


         IF ( n .eq. 1 .and. iend .eq. 1 .and.       &
     &             .not. ( isa .ge. 1 .and. irand .eq. 1 )  ) THEN 
            IF ( linkrank(n,iend) == my_rank ) THEN
            ixst = i 
            jyst = j 
            kzst = k
            ELSE
            ixst = 0 
            jyst = 0 
            kzst = 0
            
            ENDIF
          END IF
        
        END DO ! n
      END DO ! iend
      
!      write(iunit,*) 'my_rank,qtot1,qtot2,dls3 = ',my_rank,qtot(1),qtot(2),dls3
! allreduce
#ifdef MPI
       mpitotindp(1)  = qtot(1)
       mpitotindp(2)  = qtot(2)
       mpitotindp(3)  = rhotot(1)
       mpitotindp(4)  = rhotot(2)
       mpitotindp(5)  = rhotot1(1)
       mpitotindp(6)  = rhotot1(2)
       mpitotindp(7)  = chgdis0(1)
       mpitotindp(8)  = chgdis0(2)
       mpitotindp(9)  = numchancheck(1)
       mpitotindp(10)  = numchancheck(2)
       mpitotindp(11) = chanareatot(1)
       mpitotindp(12) = chanareatot(nend)
       mpitotindp(13) = chanareatot(nend+1)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 13, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       qtot(1) = mpitotoutdp(1)
       qtot(2) = mpitotoutdp(2)
       rhotot(1) = mpitotoutdp(3)
       rhotot(2) = mpitotoutdp(4)
       rhotot1(1) = mpitotoutdp(5)
       rhotot1(2) = mpitotoutdp(6)
       chgdis0(1) = mpitotoutdp(7)
       chgdis0(2) = mpitotoutdp(8)
       numchancheck(1) = NInt( mpitotoutdp(9) )
       numchancheck(2) = NInt( mpitotoutdp(10) )
       chanareatot(1) = mpitotoutdp(11)
       chanareatot(nend) = mpitotoutdp(12)
       chanareatot(nend+1) = mpitotoutdp(13)
       


      IF ( my_rank == 0 .and. (numchan(1) .ne. numchancheck(1) .or. numchan(2) .ne. numchancheck(2) ) ) THEN
      write(iunit,'(a,2(1x,i5))') 'numchancheck(1),numchancheck(2)',      &
     &         numchancheck(1),numchancheck(2)
      ENDIF
#endif
      
      DO iend=1,nend
      IF ( my_rank == 0 ) THEN
!      write(iunit,*) 'chan chg check: my_rank,qtot,dls3 = ',my_rank,qtot(iend),dls3
      write(iunit,'(a,i1,a,1pe13.5)') 'initial charge on channel qtot(',      &
     &             iend, ')= ', qtot(iend) ! *dls3
      write(iunit,'(a,1pe13.5)') 'rhotot  = ',rhotot(iend) ! *dls3
      write(iunit,'(a,1pe13.5)') 'rhotot1 = ',rhotot1(iend) ! *dls3
      IF ( rhotot(iend) .ne. 0.0 ) THEN
        write(iunit,'(a,1pe13.5)') 'fraction of charge removed = ',      &
     &       qtot(iend)/rhotot(iend)
      ENDIF
      ENDIF
      END DO ! iend

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(2x,1pe13.5))') 'count chgdis0(1,2) = ',chgdis0(1),chgdis0(2),       &
     &      chgdis0(1)+chgdis0(2)
      write(iunit,'(a,3(2x,1pe13.5))') 'channel length(1,2) = ',chanlentot(1),chanlentot(2),       &
     &      chanlentot(1)+chanlentot(2)
!      write(iunit,'(a,3(2x,1pe13.5))') 'channel area(1,2) = ',chanareatot(1)*ds2,chanareatot(2)*ds2,       &
!     &      chanareatot(nend+1)*ds2
      write(iunit,'(a,3(2x,1pe13.5))') 'channel area(1,2) = ',chanareatot(1),chanareatot(2),       &
     &      chanareatot(nend+1)
      ENDIF
      
      DO iend=1,nend
       chgdis0(iend) = 0.0
      END DO ! iend

      DO iend=1,nend
          IF ( Abs(qtot(iend)) .lt. 0.1 .and.       &
     &         .not. ( icgyn .ge. 1 .and. ibal .eq. 2 ) ) THEN
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'Not enough charge on leader',iend,      &
     &          '... go back and try again.'
            ENDIF
#ifdef MPI
!         call commasmpi_abort()
#endif
         
         t3(:,:,:) = 0.0
         t2(:,:,:) = 0.0

! 06.30.2004 [


! 06.30.2004 ]
            GOTO 199
          ENDIF
      END DO ! iend


       IF ( ibal .eq. 1 .or. ibal .ge. 3 ) THEN
         IF ( my_rank == 0 ) THEN
           write(iunit,'(a,2(1x,f8.3))')       &
     &        'fecrit(1,2) = ', fecritn(1),fecritn(2)
          ENDIF
       ENDIF

      IF ( isa .lt. 1 ) THEN
!        close(16)
      ELSE
        IF ( kmn1 .le. 7 ) THEN
          IF ( my_rank == 0 ) THEN
            write(iunit,'(a,i3)') 'Potential CG: kmn1 = ',kmn1
          ENDIF
        ENDIF
        nkzmn(kmn1) = nkzmn(kmn1) + 1
        GOTO 45
      ENDIF

! now get net charge in t3 array before smoothing:
      chgdis0(1) = 0.0
      chgdis0(2) = 0.0
      DO kz=1,nz
      DO jy=1,ny
      DO ix=1,nx
       IF ( t3(ix,jy,kz) .gt. 0.0 ) THEN
         chgdis0(2) = chgdis0(2) + t3(ix,jy,kz)*dxx(Min(ix,nx-1))*dyy(Min(jy,ny-1))*dlzlgt/z1d2lgt(kz,3) ! *dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
       ELSEIF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
         chgdis0(1) = chgdis0(1) + t3(ix,jy,kz)*dxx(Min(ix,nx-1))*dyy(Min(jy,ny-1))*dlzlgt/z1d2lgt(kz,3) ! *dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
       ENDIF
      END DO
      END DO
      END DO

! allreduce
#ifdef MPI
       mpitotindp(1)  = chgdis0(1)
       mpitotindp(2)  = chgdis0(2)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       chgdis0(1) = mpitotoutdp(1)
       chgdis0(2) = mpitotoutdp(2)

#endif


      IF ( my_rank == 0 ) THEN
        write(iunit,'(a,3(2x,1pe13.5))') 'presmooth chgdis0(1,2) = ',chgdis0(1),chgdis0(2),       &
     &      chgdis0(1)+chgdis0(2)
      ENDIF

! very light smoothing of channel charge...
      IF ( dx1/dx .gt. 1.7 ) THEN
        fx = 0.5 ! changed from 1.0 9/12/2006 erm
        fz = 1.0 ! changed from 0.5 9/12/2006 erm
        i = 1
      ELSE 
        fx = 1.0
        fz = 1.0
        i = 1 ! 2 ! changed 9/12/2006 erm
      ENDIF
      IF ( dz .gt. 250. ) i = 2
!      i=0

! #ifndef MPI
      DO j=1,i
       CALL smthlt1mpi (nx,ny,nz,nor,t3,t4,fx,fz) 
      ENDDO
! #endif
      
! now set t4 array with flags of where lightning charge is dumped
      chgdis0(1) = 0.0
      chgdis0(2) = 0.0
      t4(:,:,:) = 0.0
      DO kz=1,nz
      DO jy=1,ny
      DO ix=1,nx
        IF ( Abs( t3(ix,jy,kz) ) .gt. 1.0e-28 ) THEN
          t4(ix,jy,kz) = Sign(1.0,t3(ix,jy,kz) )
        ELSE
          t4(ix,jy,kz) = 0.0
        END IF
       IF ( .not. ( t3(ix,jy,kz) .lt. 1000. .and. t3(ix,jy,kz) .gt. -1000.) ) THEN
         write(0,*) 'PROBLEM with t3! ',ix,jy,kz,t3(ix,jy,kz) 
         call commasmpi_abort()
       ELSEIF ( t3(ix,jy,kz) .gt. 0.0 ) THEN
         chgdis0(2) = chgdis0(2) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dlzlgt/z1d2lgt(kz,3) !*dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
!         print*,' t3 = ',ix,jy,kz,t3(ix,jy,kz)
       ELSEIF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
         chgdis0(1) = chgdis0(1) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dlzlgt/z1d2lgt(kz,3) ! *dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
!         print*,' t3 = ',ix,jy,kz,t3(ix,jy,kz)
       ENDIF
      END DO
      END DO
      END DO

! allreduce
#ifdef MPI
       mpitotindp(1)  = chgdis0(1)
       mpitotindp(2)  = chgdis0(2)

       CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION,       &
     &                   MPI_SUM, my_comm, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN

       chgdis0(1) = mpitotoutdp(1)
       chgdis0(2) = mpitotoutdp(2)

#endif

      
      IF ( my_rank == 0 ) THEN
        write(iunit,'(a,3(2x,1pe13.5))') 'recount chgdis0(1,2) = ',chgdis0(1),chgdis0(2),       &
     &      chgdis0(1)+chgdis0(2)
      ENDIF
!      call commasmpi_abort()

 999  CONTINUE 

!      deallocate ( mask )
!      deallocate ( phip )
   
      deallocate ( hgt )
      deallocate ( zmap )
      deallocate ( cofx )
      deallocate ( cofy )
      deallocate ( cofz )
      deallocate ( cz )
      deallocate ( czz )
      ! deallocate ( gd ) 
      deallocate ( ipick )
      deallocate ( footprint )

      RETURN
      END   ! of mudslpm

#else
      RETURN
      END
#endif

! #####################################################################
!
!  subroutine smthlt1mpi
!
! #####################################################################
!
      subroutine smthlt1mpi (nx,ny,nz,nor,a,b,f,fz)
!
      
      USE COMMASMPI_MODULE
      
      implicit none
!
      integer   ia,na
      integer   ix,jy,kz
      integer   id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer   ism,jsm,ksm
      integer   nx,ny,nz
      integer   nor
      real f,fz
!
      real      a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real      b(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      integer ixb,ixe,jyb,jye

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      
! #################################################################
      ism = 2
      jsm = 2
      ksm = 2
!
! copy array
!

      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
      b(ix,jy,kz) = a(ix,jy,kz)
      end do
      end do
      end do

#ifdef MPI        
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,nor,1,      &
     &        w_proc(my_rank),e_proc(my_rank),westward_tag,b)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,nor,1,      &
     &        w_proc(my_rank),e_proc(my_rank),eastward_tag,b)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,nor,1,      &
     &        n_proc(my_rank),s_proc(my_rank),southward_tag,b)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,nor,1,      &
     &        n_proc(my_rank),s_proc(my_rank),northward_tag,b)

#endif
!
!  smooth
!
      ixb = 1
      ixe = nx
      
      IF ( myproci == 1 ) ixb = 1+id1
      IF ( myproci == nproci ) ixe = nx-ism*id1
      
      jyb = 1
      jye = ny
      
      IF ( myprocj == 1 ) jyb = 1+jd1
      IF ( myprocj == nprocj )  jye = ny-jsm*jd1
      
      do kz = 1+kd1,nz-ksm*kd1
      do jy = jyb,jye ! 1+jd1,ny-jsm*jd1
      do ix = ixb,ixe ! 1+id1,nx-ism*id1
      a(ix,jy,kz) = b(ix,jy,kz)      &
     &  +f*(1.0/16.0)*(b(ix+id1,jy,kz)-2.*b(ix,jy,kz)+b(ix-id1,jy,kz))      &
     &  +f*(1.0/16.0)*(b(ix,jy+jd1,kz)-2.*b(ix,jy,kz)+b(ix,jy-jd1,kz))      &
     &  +fz*(1.0/16.0)*(b(ix,jy,kz+kd1)-2.*b(ix,jy,kz)+b(ix,jy,kz-kd1))
! new values  changed 9/12/2006 erm
!     >  +f*(1.0/30.0)*(b(ix+id1,jy,kz)-2.*b(ix,jy,kz)+b(ix-id1,jy,kz))
!     >  +f*(1.0/30.0)*(b(ix,jy+jd1,kz)-2.*b(ix,jy,kz)+b(ix,jy-jd1,kz))
!     > +fz*(1.0/100.0)*(b(ix,jy,kz+kd1)-2.*b(ix,jy,kz)+b(ix,jy,kz-kd1))
      end do
      end do
      end do
!
!  return
!
      return
      end
!
!
!
!
!-------------------------------------------------------------------------
!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                 SUBROUTINE MAXFIELD                  ######
!     ######                written by E. Mansell                 ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
      SUBROUTINE MXFLDSTRZmpi(nx,ny,nz,      &
     &    phi,gd,gdx,gdy,pstop,nbw,nbs,      &
     &    emax,emagmax,i1,j1,k1,i2,j2,k2,dx,dy,dz,istretch,      &
     &    z1d2lgt,nz1d,ezmax,iunit,bcx,bcy)
      
      USE COMMASMPI_MODULE
      implicit none
      
      integer nx,ny,nz
      integer iunit
      integer bcx,bcy
      real emag
!      real phi(-1:nx+2,-1:ny+2,0:nz+1)
      real phi(0:nx+1,0:ny+1,0:nz+1)
      real emax,emagmax,exmax,eymax,ezmax
      integer i1,j1,k1,i2,j2,k2
      integer jst,jend
      integer ix,jy,kz,i,j,k
      integer nbw,nbs
      real dx,dy,dz
      logical istretch
      real dxinv2,dyinv2,dzinv2
      real ex,ey,ez
      real gd(-1:1,-1:1,-1:1,0:nz,2)  !  1: distance, 2: inverse of distance
      real gdx(0:nx,2), gdy(0:ny,2)

      real etest
      integer pstop
      integer nz1d
      real z1d2lgt(nz1d,4)
      integer, parameter :: ndebug = 0
      
      integer ixe, jye
      integer :: nz1

      integer :: itype,imxmn
      integer :: irank1,irank2
   
      integer :: listsize

      integer :: mxmn_size,max_rank,min_rank
      integer :: loc_send_size,loc_recv_size
      integer :: exyz_send_size,exyz_recv_size

      integer,allocatable,save :: mpiloc(:,:), mpiloc0(:,:,:)
      integer,allocatable,save :: loc_send_buffer(:,:),loc_recv_buffer(:,:)

      real,allocatable,save :: mpiexyz(:,:), mpiexyz0(:,:,:)
      real,allocatable,save :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      real,allocatable,save :: max_send_buffer(:,:),max_recv_buffer(:,:)
      real,allocatable,save :: min_send_buffer(:,:),min_recv_buffer(:,:)

      real,allocatable,save :: mpimxmn(:,:)
#ifdef MPI
      INCLUDE "mpif.h"
#endif

!#######################################################################
!
!     Begin executable code
!
!#######################################################################

#ifdef MPI
      nz1 = 1
      
      IF ( .not. allocated( mpimxmn ) ) THEN
      allocate(mpimxmn(nz1,4))
      allocate(max_send_buffer(2,nz1))
      allocate(max_recv_buffer(2,nz1))
      allocate(min_send_buffer(2,nz1))
      allocate(min_recv_buffer(2,nz1))

      allocate(mpiloc(nz1,3))
      allocate(mpiloc0(nz1,2,0:number_of_processes-1))
      allocate(loc_send_buffer(nz1,3))
      allocate(loc_recv_buffer(nz1*number_of_processes,2))

!      real,allocatable :: mpiexyz(:,:), mpiexyz0(:,:,:)
!      real,allocatable :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      allocate(mpiexyz(nz1,3))
      allocate(mpiexyz0(nz1,3,0:number_of_processes-1))
      allocate(exyz_send_buffer(nz1,3))
      allocate(exyz_recv_buffer(nz1*number_of_processes,3))
      ENDIF

#endif
      
      i1 = 0
      j1 = 0
      k1 = 0
      i2 = 0 
      j2 = 0
      k2 = 0 
      emagmax = 0.0
      dxinv2 = 0.5/dx ! 0.5 because we are doing a centered difference
      dyinv2 = 0.5/dy
      dzinv2 = 0.5/dz
      IF ( ny .gt. 2 ) THEN
       jst = 2
       IF ( myprocj > 1 ) jst = 1
       jend = ny-1
       jye = ny
       IF ( myprocj == nprocj ) jye = ny-1
      ELSE
       jst = 1
       jend = 1
      ENDIF
      
      ey = 0.0
      
#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: calculate E'
#endif
      ixe = nx
      IF ( myproci == nproci ) ixe = nx-1
      
      DO k=1,nz-1
      DO j=jst,jend
      DO i=1,ixe
         ex = (phi(i-1,j,k) - phi(i+1,j,k))*dxinv2 ! *gdx(i,2)
         IF ( ny .gt. 2 ) THEN
           ey = (phi(i,j-1,k) - phi(i,j+1,k))*dyinv2 ! *gdy(j,2)
         ENDIF
         IF ( istretch ) THEN
           ez = (phi(i,j,k-1) - phi(i,j,k+1))*dzinv2*z1d2lgt(k,3)
         ELSE
           ez = (phi(i,j,k-1) - phi(i,j,k+1))*dzinv2
         ENDIF
         emag = Sqrt(ex**2 + ey**2 + ez**2)
         IF ( emag .gt. emagmax ) THEN
           emagmax = emag
           exmax = ex
           eymax = ey
           ezmax = ez
           i1 = i
           j1 = j
           k1 = k
         END IF
      END DO  ! i
      END DO  ! j
      END DO  ! k
      
#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: myrank,ex,ey,ez = ',my_rank,exmax,eymax,ezmax
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: myrank,emag,i,j,k = ',my_rank,emagmax,i1,j1,k1
#endif
#ifdef MPI
      DO kz = 1,nz1 ! note: nz1 = 1

       max_send_buffer(1,kz) = emagmax
       max_send_buffer(2,kz) = real(my_rank)

      ENDDO
    
      mxmn_size = nz1

! find out the max emag at each level, and which processor it is on
      CALL MPI_Allreduce(max_send_buffer, max_recv_buffer, mxmn_size, MPI_2REAL,        &
     &                  MPI_MAXLOC, my_comm, mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Allreduce attempt"

        CALL MPI_Abort(my_comm, mpi_error_code,mpi_abort_error_code)
      ENDIF

! put global maximum emax into mpimxmn
      emax = 0.0
      DO kz = 1,nz1 ! note nz1 = 1

       emagmax = max_recv_buffer(1,kz)
!       emaxkz(kz) = max_recv_buffer(1,kz)

      ENDDO

      irank1 = int(max_recv_buffer(2,1))

! send e-field components of emag max
      DO kz = 1,nz1
       exyz_send_buffer(kz,1) = exmax
       exyz_send_buffer(kz,2) = eymax
       exyz_send_buffer(kz,3) = ezmax
      ENDDO

      exyz_send_size = nz1*3
      exyz_recv_size = nz1*3

      CALL MPI_Gather(exyz_send_buffer,exyz_send_size,MPI_REAL,      &
     &               exyz_recv_buffer,exyz_recv_size,MPI_REAL,       &
     &               0,my_comm,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE (0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"
        CALL MPI_Abort(my_comm, mpi_error_code,mpi_abort_error_code)
      ENDIF

      IF ( my_rank == 0 ) THEN
        mpiexyz0 = reshape(exyz_recv_buffer,(/nz1,3,number_of_processes/))

           DO kz = 1,nz1
       
            irank1 = int(max_recv_buffer(2,kz))

            exmax = mpiexyz0(kz, 1, irank1)
            eymax = mpiexyz0(kz, 2, irank1)
            ezmax = mpiexyz0(kz, 3, irank1)

          ENDDO
        ENDIF

#else
      irank1 = my_rank
#endif

#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ.2: myrank,ex,ey,ez = ',my_rank,exmax,eymax,ezmax
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ.2: myrank,emag,i,j,k = ',my_rank,emagmax,i1,j1,k1

      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: check for max on boundary'
#endif
      IF ( irank1 == my_rank .and. ((.not. xperiodic .and.       &
     &    ( (i1 .lt. nbw+2 .and. myproci == 1)  .or.       &
     &      (i1 .gt. nx-nbw-2.and. myproci == nproci))) .or.      &
     &     (.not. yperiodic .and. ny .gt. 2  .and.       &
     &      ( (j1 .lt. nbs+2 .and. myprocj == 1) .or.       &
     &        (j1 .gt. ny-nbs-2 .and. myprocj == nprocj))) .or.      &
     &      k1 .lt. 1 .or. k1 .gt. nz-2 ) .and.      &
     &   emagmax .gt. 0.0   ) THEN
         write(iunit,*) 'Harumph! Emagmax along an edge! I dont like it!'
         write(0,*) 'Harumph! Emagmax along an edge! I dont like it!'
         write(0,'(a,3(1x,i5),1x,e12.5)') 'i1,j1,k,emagmax',      &
     &         i1,j1,k1,emagmax,myproci,my_rank
         pstop = 1
      END IF
      
!  find direction of largest 1st order electric field
      IF ( ny .gt. 2 ) THEN ! 3D
        jst = -1
        jend = 1
      ELSE        ! 2D
        jst =  0
        jend = 0
      ENDIF
 
#if DEBUG > 0
       IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: find direction of E'
#endif
      emax = 0.0
      DO k= -1,1
      DO j= jst,jend
      DO i= -1,1
         etest = Abs(phi(i1+i,j1+j,k1+k) - phi(i1,j1,k1))*gd(i,j,k,k1,2)
         IF (etest .gt. emax) THEN
           emax = etest
           i2 = i1+i
           j2 = j1+j
           k2 = k1+k
         END IF
         
      END DO ! i
      END DO ! j
      END DO ! k

#if DEBUG > 0
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: rank,emax = ',my_rank,emax

      IF (ndebug .ge. 2) write(0,*) 'MXFLDSTRZ: return'
#endif

      RETURN
      END



