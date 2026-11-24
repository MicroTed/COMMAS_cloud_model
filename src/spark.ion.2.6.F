CVD$F SKIP
c#######################################################################
c  SUBROUTINE SLPM
c
c  written by ERM      
cc
c
c
c 'phi2' for initial solution, copied into 'phi' which is used for 
c   propagation (values changed by lightning channels)
c 'pot2' was same as 'phi2' but extended one point below ground.
c 'phi2w' and 'pot2w' were set up so that (1,1,1) corresponds (1,1,1) on the
c  storm grid.
c
c  so now:
c    pot2 -> phi
c    pot2w -> phi0
c    f -> rhs
c    fw -> rhs0
c    phi -> phip  (phi-propagate)
c
c   can get rid of phi2 and phi2w
c
c phi, phi0, and phi2 are all the same array expressed differently (pot2)
c       they are all the 'background' (i.e. preflash) potential
c phip (phi-propagate) is a copy of phi that can be updated with 
c       channel potentials
c
c#######################################################################
c##################################################################
c

      subroutine mudslpm
     >  (time_real,time,nstep1,iliter,lgtstp,pchnl,nchnl,ldumpsp
     >  ,nor,nx,ny,nz,pot,exdg
     >  ,eydg,ezdg,emagdg
     >  ,icldtop,zcldtop
     >  ,icgyn,icgpn,ebrkds
     >  ,dtp,dx1,dy1,dz1,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
     >  ,dzmax, istretch
     >  ,gxt,gyt,gzt,dxx,dyy,dzz
     >  ,ixst,jyst,kzst,outname,nkzmn,icheckdims,
     :   nnx,nny,nnz,nnzt,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,nbz2,
     :   llwork,
     :   nxs,nys,nzs,nb,nxb,nyb,nzbt,nzb,nxyz,nwslv,
     :   phi,phi0,phi2,rhs,rhs0,rho,rhow,work,
     :   numic,numcgn,numcgp,
     :   iunit, bcx, bcy)
     
!       USE GRID_MODULE
       USE ELEC_MODULE

c#######################################################################
            
c##################################################################
c  PURPOSE: To generate a lightning channel using the Dielectric
c       Breakdown Model of L. Niemeyer, L. Pietronero, and 
c       H.J. Wiesmann (Phys. Rev. Let., 1984, v. 52, p. 1033)
c
c
c  TO DO list:  vertical refinement on stretched grid
c
c##################################################################
c
c 2006.9:  version 2.6:  Allow for 2-D grid and use staggered E for emag
c
c 2005.10.3: Added input variable cgfr to control whether ibal=2 resets
c            the channel potential along the channel or not.  (A value 
c            of cgfr < 1.0 turns this off an cgfrac = cgfr for the fraction
c            of charge retained on the lower channels.)
c
c 2005.05.25  CHANGED the grid refinement so that the smaller boxes fit inside
c    the larger boxes.  This will preserve the total charge, but grid 
c    points will not be co-located with any coarse grid points.  This
c    requires that dxl = n*ds (integer multiple).
c      In other words, now do
c
c     |   *   |   *   |   *   |
c to
c     | *   * | *   * | *   * |  (igrid=1)
c
c   instead of 
c 
c     |   *   |   *   |   *   |
c to
c     |   * | * | * | * | *   | (igrid=0)
c
c  2007.03.22: erm; changed back to igrid = 0.
c
c  This also makes putting charge back on the coarse grid much simpler!
c
c 2005.01.14: Try normalizing phistart hi and lo to a length of 1km by adding
c             or subtracting the quantity  1.e3*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
c             to phistart: 
c         phistlow = phistart - 5.e2*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
c         phisthi  = phistart + 5.e2*gd(i2-i1,j2-j1,k2-k1,k2,2)*phidiff
c
c         Also change count of one-sided growth to be normalized to a dslight=500m:
c
c            IF ( ibal .eq. 2 .and. istopcnt(iend) .gt. 50*500./ds ) THEN
c          
c         so that a flash can go about the same total length without 
c         growth on the other end
c
c
c 8.2004: some tests with setting channel potential to zero upon becoming a CG
c         if ibal = 2
c
c 6.30.2004: Fixed a bug in the final calculations of the channel charge. 
c            buggy version had
c                IF ( k1 .gt. 1 .and. k1 .lt. kz-2 ) THEN
c         changed to
c                IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN
c         Didn't see an obvious problem with this until trying to rerun
c         the WK Us=50 case.  It could have been a culprit in causing 
c         a sizable fraction of a channel not to be counted, sometimes resulting
c         in a retry of the flash.  It would also affect the actual flash 
c         charge.
c           The problem appeared after a change was made to the z-limit of
c         the initiation point list, which was change from kz=1,nzs-4 to 
c         kz=1,k2 where k2 is the point just below cloud top.  This left kz
c         defined as k2+1, so the channels above k2-2 were being cut off, which
c         is a big problem for flashes that initiate close to cloud top, as 
c         they apparently do in the Us=50 WK supercell case.
c
c         Also removed 'the top' condition of numchan(iend).gt.numchan(3-iend) 
c         for ibal=2.  Not much point to that, really.
c
c 5.06.2004: For ibal=2, added a constraint that the channel reference potential
c            (phistart) cannot be shifted beyond the values of the first 
c            two channel points.
c
c 4.17.2004: Found strange logic in the ibal=5 scheme.  It wasn't really
c            doing ihold like it should and was not stopping the flash
c            if it couldn't force growth on the less-charged end.  The
c            changes cause more total steps in the process, but those are
c            pretty cheap compared to the poisson solver.  Also causes a higher
c            flash rate for the same fecrit
c
c
c 4.8.2004: increased the number of times that charge balance is checked
c           for ibal=5.
c      
c            was:
c      IF ( (( nstep/5*5 .eq. nstep .or. nstep .le. 3 )  .and.
c
c           now:
c      IF ( (( nstep/4*4 .eq. nstep .or. nstep .le. 30 )  .and.
c
c
c 2.14.2004: Two bug fixes: cofz was not set for k=0 and k=1, causing 
c            possible problems with strzsor.  Other bug was incorrect
c            setting of chandist when a new channel point is added.
c
c 1.27.2004: fixed a bug in the calculation of Ez in the initiation point
c            list generation code.  It thought that dz was always dzmax
c            (wrong subscript in z1d2(kz,3))
c
c There was also a bug in the calculation of qtot for the charge balance
c checks (though not at the very end where charge is put into the 'coarse' grid).
c The charge density was being scaled by the stretch factor, but the wrong 
c index was being used (k2 instead of k1).  Not a huge problem, but good to
c have it fixed.
c
c 01.05.2004 Update/bug fixes for ibal=2 option, which seems to be behaving
c            better since other bugs were fixed.  Apparently need to use a 
c            relatively low internal field (1kV/m or less).  Added full flash SOR
c            after adjusting channel potential.  Had a bug that caused internal
c            field to be applied only to the first points, with zero drop after that.
c            Applied the channel potential criterion for ground flashes, but
c            it seems questionable to do so in this case, since the channel is
c            artificially widened, thus diluting the value of the channel potential.
c            Must think about this some more!
c
c            Also leave charge on CG channels as for ibal=5.
c
c
c
c 07.25.2003 sam.spark.ion.2.5.f: update to allow unstretched grid with 
c            dslight < dzl
c
c 01.31.2003 update to documentation: ibal=5 acts like ibal=4 and also
c            feeds back charge from the ground leaders (zeroing points that
c            have the wrong sign of charge, which can happen on the part of the 
c            channel that is extended to ground when a threshold is used.)
c
c            Both ibal=4 and ibal=5 add a second criterion for deciding ground
c            flashes: the potential of the downward leader tip must still have
c            the same sign it started with (i.e. a negative leader tip must have
c            a negative channel potential.)
c
c            Also removed some old commented code to clean things up a bit.
c
c 10.22.2002 version 2.4c:
c            Added ibal = 4 option that is like ibal=3 but does not stop
c            an unbalanced flash if the the downward channels have more 
c            charge.  
c            planned to add: preserve charge on downward CG channels
c
c              In sll, ibal=4 will limit charge on ground channels to 
c            -0.25*(charge on upper channels) (1/31/2003 Now only ibal=5 does this)
c
c 9.16.2002 increased search range for initiation points to get closer
c           to the edge of the domain:
c
c      DO jy = nbs+3,nny-nbn-3
c      DO ix = nbw+3,nnx-nbe-3
c  becomes
c
c      DO jy = nbs+2,nny-nbn-1
c      DO ix = nbw+2,nnx-nbe-1
c
c  also
c     :         i1+i .gt. nbw+1 .and. i1+i .lt. nxs+nbw-1 .and.
c     :         j1+j .gt. nbs+1 .and. j1+j .lt. nys+nbs-1 .and.
c  to
c     :         i1+i .gt. nbw+1 .and. i1+i .le. nxs+nbw .and.
c     :         j1+j .gt. nbs+1 .and. j1+j .le. nys+nbs .and.
c
c 8.26.2002 Had a problem with the growth on the upward channels of a CG
c           also getting to ground.  Suppressed as follows: 
c            1. the flash is stopped if the upper end gets within 0.5*cgthres
c               of the ground
c            2. the upper branches are allowed to grow a maximum of 50 extra 
c               channel segments
c           Also (mostly) fixed an occasional problem where the propagation
c           thresholds of both ends could go below fecrit.  Now, when one end
c           lowered, the other end is raised back toward fecrit (to a maximum
c           of fecrit).
c
c 8.20.2002 Originally left all charge on channel for CGs when ibal=3 but
c           now changed back to the original scheme of setting charge to 
c           zero on downward leaders.  Otherwise can get weird stuff like
c           net negative charge to the storm by a negative CG!
c
c 8.6.2002  Increased CG threshold to 2500m (was 2000m).  
c           Also changed ibal=3 scheme to adjust propagation if imbalance
c           exceeds 5% (was 10%)
c
c
c 5.22.2002  New version (sparkmud.2.4b.f) for vertically stretched grid.
c            For the moment must have dslight = dzmax (i.e. lightning uses the
c            same vertical levels as on the dynamics grid and refinement allowed
c            in horizontal only)  There is no longer a choice about whether to set
c            phi=0 at kz=1 or kz = 0.  Now hard-wired for phi=0 at kz = 0 (i.e.
c            boundary is below ground).
c
c            Also fixed a bug in the code that puts lightning channel charge
c            back into the coarse (dynamics) grid.
c
c 4.8.2002  Fixed bug where still had mseed used for ichoose instead of iseed
c
c 1/17/2002  Now have more variables that are from inlightdbm,
c             including 
c                   ieint (selects constant or fractional internal field)
c                   eint1 (input value for constant internal field)
c                   ibal  (selects channel neutrality balancing method)
c                   fefac (individual propagation settings for pos/neg ends)
c              
c              Added code for the flag 'ieint' to set internal field values
c              
c              Made some changes to ibal=2 method:
c                If an initial potential shift causes imbalance to 'go the wrong way,'
c                then dont try to adjust it on the current step and instead reset 
c                ihold = 0 and let the channel grow a little more.  This seems to be
c                working well.
c                
c                 Also now check channel balance on every step up to 40 
c                 (instead of 4 previously) and after that only every 3rd step
c                 
c
c 1/7/2002  spark.2.4b.f   Try using a second-order difference to calculate e-field
c           next to channel when possible (otherwise revert to 1st order difference)
c
c 1/4/2002  spark.2.4a.f   Local arrays are now automatic or passed in from lightning
c            driver.  This allows changing dimensions without resetting the lightning
c            grid dimensions.
c
c 12/3/2001  Minor fix to reject small flashes (< 4 segments total)
c            when growth cant be forced.
c
c 11/26/2001 version sam.spark.2.4.f:
c            Combine regular and stand-alone versions into one.  Stand-alone
c            behavior is turned on by isa = 1
c
c 10/3/2001 2.3f.f
c    Fix up to allow arbitrary ds (i.e. anything other than the storm dz).

c 9/26/2001
c   increased eint to 5kV/m (from 500 V/m)
c
c 9/25/2001  spark.2.3e.f
c    Removed test for channel outside cloud and put in a difference in prop. threshold
c    for neg. vs. pos. leaders.  Start out with neg. leaders needing 20% higher field.
c
c  9/12/2001 sam.spark.2.3d.f:
c      Increased domain height to improve solution.
c  
c  9/10/2001  sam.spark.2.3c.f:
c        Correct the transfer of lightning charge to the coarse grid.   
c  9/8/2001 sam.spark.2.3b.f
c
c        Changed array t4 to have its logical meaning (i.e. t4 is positive if
c          t3 is positive)
c  8/22/2001
c       put in stronger smoothing (horizontal only) when the lightning grid dx is
c       close to the dynamics dx (dx1) (i.e. when dx1/dx .lt. 1.5 ).  If dx1/dx .gt. 1.5 
c       then the injection into the coarser grid has enough horizontal smoothing already.
c       (The reason for doing this is that a branch from one flash went up above the
c        cloud, and the charge it deposited up there, as 'ions' then instigated another
c        flash.  The second flash initiated _above_ the storm -- not very realistic!)
c
c  version sam.spark.2.3a.f:
c     Still not quite working, so grafted in old code for initial channel...
c     The problem turned out to be accessing the wrong array for ambient charge
c       when calculating the channel charge (should have used rho instead of rhow,
c       resulting in an offset of (nb,nb) gridpoints horizontally.  Oops.)
c
c  version sam.spark.2.3.f:
c     Returned to 2-point initial channel, which seemed to help.
c     
c  7/31/2001 version sam.spark.2.2.f
c            Now trying the 'proper' way to balance charge: by adjusting
c            the reference potential.
c
c            Start using a 3-point initial channel instead of 2-point.  The
c            positive and negative leaders now 'share' the middle point and
c            half of that charge.
c
c            Also fixed up 'signchan' to actually be the sign of charge
c            carried by that end.  (Was reversed before for arcane reasons)
c
c            Also fixed up accounting of the charge on the channel by
c            subtracting the ambient charge density.  This should actually
c            enhance slightly the charge on the channel.
c
c  7/31/2001 version sam.spark.2.1.f
c            Upon suggestion by D. MacGorman, try a new charge balancing 
c            procedure by which the lopsided end is kept from growing 
c            at all until the other end catches up.  (Works pretty well, too.)
c
c
c  Expanded version to feed potential back into SOR  
c  spark2r uses single precision arrays instead of double
c
c added in a call to the FACR routine in FISHPACK to compare with SOR
c   The FACR solution to the initial charge field is fed to the SOR as
c   a first guess.  Then SOR can work on the short wavelengths.
c
c 4/15/99  added diagonal propagation possibility
c
c 4/15/99  Now finally try it with a dipole distribution!
c
c 4/16/99  Now read in a space charge field from the cloud model
c          interpolates from a 1500x1500x500m  grid to 500x500x500m
c
c 4/18/99  Go down to 250m resolution for lightning propagation
c
c 4/21/99  Altered version of spark6b to simulate the Williams (1985)
c          experiments with space charge embedded in PMMA
c 4/25     change to allow for reversal of charge regions
c
c 5/1/99  version sparkw2 for simulating a rectangular slab of charge
c         with the discharge started at one end (Williams et al. (1985),
c         figure 1.)
c
c 5/2/99  Add calculation of charge density at each point in the
c         channel and just add up all the adjustments from the 
c         initial values.
c
c   5/5/99  rename from sparkw3 to sparkc1 and try again with a 
c           cloud model field
c
c  5/9/99  sparkc2
c          Try doing 'local' relaxation: just around the added bond
c
c  5/10/99 sparkc3: 500m grid
c
c  5/25/99 sparkc6: add break-even field option
c
c  6/15/99 sparkc7: testing speed-up on 500m grid
c
c  6/26/99 convert to a subroutine to be called from SAM lightning
c          subroutine
c
c  7/5/99  put a buffer zone extension to the lateral boundaries
c          so there is a charge-free zone next to the boundaries
c
c  8/8/99 allow propagation thresholds to adjust to promote growth 
c         on an end with substantially less charge and retard growth
c         of the end with more charge.  Idea is to check the channel 
c         charges every so many steps and alter the propagation  
c         threshold fraction fecrit for each end.  This should have
c         a similar effect as adjusting the channel potential.
c
c  8/11/99 increase minimum fecrit from 0.5 to 0.6
c  8/12/99 output flash bounding box, put in logical switch 'lsgi'
c          to easily switch between sgi and cray timing routines
c
c          If a CG occurs, don't record the channel charge of the 
c          ground-reaching end in the t3 array (but do calculate the
c          original charge on that end, anyway, and write to output
c          files)  Also, if fecritn of the upward-propagating end
c          back was greater than fecrit, then reset it to fecrit.  
c
c 8/22/99  When calculating probabilities, now check that the local
c          charge density of the bond point is > 1.0e-16.  If it is
c          not, then set Ecrit = Min ( 150kV/m, 2*Ecrit(z) ) to
c          simulate the increased breakdown field for clear air vs.
c          cloudy/hydrometeor-filled air.
c
c 9/12/99  Now after 3 flashes in the same time step we stop using
c          the offset and only start lightning where E exceeds 
c          inititation threshold 
c
c          Also don't relax after the initial starting point is
c          chosen (hope this will help avoid 1-segment channels)
c          (changed back again 9/17)
c
c 9/17/99  Now pick a new starting point if propagation stops with
c          4 or fewer total channel points (do this up to maxtry times-- it's 
c          a lot cheaper than doing the potential all over for a trivial
c          flash)  
c
c          Also reinstate the initial relaxation around the 
c          starting points.
c
c          Also treat the initial 2-point breakdown as a perfect 
c          conductor.  i.e. no internal field for the first points.
c
c 10/2/99  Change to set zero potential below ground level as in
c          sam.08.driver.sgi.f
c
c 10/4/99  set mask array values back to 1 when trying a new point after
c          a small flash
c
c 10/26/99 Stop IC development after 40 steps if charge imbalance
c          is greater than a facter of 2. (CGs are left alone)
c
c 2/12/2000 now stop flashes if level nz-2 is reached (used to be nz-1)
c           Also, changed
c              ecomp = Min( 150.0e3, 2.0*ecritz(k1) )
c           to
c              ecomp = Min( 150.0e3, Max(fecritn(iend),2.0)*ecritz(k1) )
c           just in case fecritn(iend) happens to be greater than 2.
c
c 3/15/2000 Changed lightning file naming convention.  Just reordered
c           the model time in front of the iteration number so that a normal 
c           file listing will put the files in chronological order.
c
c TO DO:
c
c        1. put in check that grid dimensions are OK
c
c#######################################################################
c
c  Variable declarations
c
c#######################################################################
      
      implicit none

!      include 'sam.files.h'

      integer nx,ny,nz,nb,nor  ! nb is number of boundary buffer zone points

      integer, parameter    :: ng1 = 1
      
      integer               :: iunit
      integer               :: bcx, bcy

      real dzz(nz),dxx(nx),dyy(ny)         ! dz(k),dx(i),dy(j)
c      real    :: gxt(nx,4), gyt(ny,4), gzt(nz,4)
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      real ex,ey,ez,emag,dxinv2,dyinv2,dzinv2
      real ex1,ey1,ez1,emag1, dxinv, dyinv, dzinv
      
      double precision tmp
      integer ntry
      
!      real cgfr
      
      integer idown ! which end is going downward
      real phi1
      integer numic,numcgn,numcgp
      
      integer icnt1(2),icnt2(2),icnt3(2),icnt4(2)
      real chgdis(2),chgdis0(2)
      
!      real cgthres       ! altitude threshold to declare a CG
      real chanlen       !  channel length from initiation point to ground
      integer iextra     ! counter for max allowed extra segments after CG
      integer imapz, mzdist, dzmax
      logical istretch  ! whether grid is stretched or not
      logical useinterp ! whether to use interp grid (true) or coarse grid (false) for finding init points
      logical ldumpsp   ! flag to make special dump when something goes screwy
      logical lct       ! above-cloud energy test
      integer iixp,jjyq,kkzr,iiex,jjey,kkez
      integer*8 llwork
      integer nnx,nny,nnz
      integer nxs,nys,nzs   ! indices of 'interior' of slm grid
      integer nbw,nbe,nbs,nbn,nbz2
      
      real eold,enew,eoldct,enewct
      real a,b,c,d
      real etst
c
c note that nnzt = nnz-1
c
      integer nnzt
      real phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      real phi0(-nbw:nxs+nbe,-nbs:nys+nbn,0:nnzt)
      real rhs0(-nbw:nxs+nbe,-nbs:nys+nbn,0:nnzt)
      real phi2(nnx,nny,nnzt)
      
!      real potfair(nnz),rhofair(nnz)

      real rho(nnx,nny,nnzt)
      real rhow(-nbw:nxs+nbe,-nbs:nys+nbn,nnzt)

      integer nzmax
!      parameter (nzmax1 = 400)
!      real z1d2(nzmax1,4)
      real dlzz,dlz2
!      common/strz/z1d2,dlz,nz1d

      real eztop, pottop
      common/bndcz/eztop, pottop
      
      real,allocatable :: hgt(:),zmap(:), hgtdg(:)
      real,allocatable :: cofx(:,:),cofy(:,:),cofz(:,:)
      real,allocatable :: cz(:),czz(:)

c      real hgt(nnz),zmap(nnz)
c      real cofx(nnx,3),cofy(nny,3),cofz(0:nnzt,3)
c      real cz(0:nnz-1),czz(0:nnz-1)

c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling
c
      integer*8 iprm(22)
      integer mgopt(4)
      real fprm(8)


      integer*8 intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
      common/itmsp3/intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc0,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc0,yd,ze,zf,tolmax,relmax
      real dlx,dly ! ,dlz2,x,y,z,cxx,cyy,czz,cx,cy,cz,cex,cey,cez,ce,errm

      equivalence(intl,iprm)
      equivalence(xa,fprm)

      integer*8 iprm2d(17)
      real fprm2d(6)

      integer*8 intl2d,nxa2d,nxb12d,nyc2d,nyd2d,ixp2d,jyq2d,iex2d,jey2d,
     +              nxm2d,nym2d,
     +              iguess2d,maxcy2d,method2d,nwork2d,lwrkqd2d,itero2d
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
       common/itmsp2d/ intl2d,nxa2d,nxb12d,nyc2d,nyd2d,ixp2d,jyq2d,iex2d,jey2d,
     +              nxm2d,nym2d,
     +              iguess2d,maxcy2d,method2d,nwork2d,lwrkqd2d,itero2d
      real xa2d,xb2d,yc2d,yd2d,tolmax2d,relmax2d
      common/ftmsp2d/xa2d,xb2d,yc2d,yd2d,tolmax2d,relmax2d

      equivalence(intl2d,iprm2d)
      equivalence(xa2d,fprm2d)
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,cfz,cfzlgt,cfzu,bndc,bndc2d


      double precision e,e1     ! temp value of efield next to channel
!      integer l2nde  !  =2 to use 2nd order difference for e
                     !  =1 for 1st order
c      parameter (l2nde = 1)
      integer icheckdims
      integer icldtop
      real    zcldtop
      integer ihold
      real fx,fz
      real scintmx, scintmn 
      real phidiff
      real bfac

!      integer ibal  ! flag for which balancing procedure to use

      real qtot1,qtot2
      integer iadd,iaddz,istop1
      real chg,chg0
      integer ilast, iflip
      real adjlast
      
      integer iw, ie, js, jn, ku, kd
      real    xw, xe, yn, zu, zd
      real   xc, yc, zc
      real frs, frn, fre, frw, fru, frd


      
c input stuff
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


                 
c nsavephi:  =2 or 3 read facr calculated phi from file,1=write facr phi to file
c   =3 write phi after sor, =4 read sor-relaxed phi (skips first relaxation)
c  nsavephi =1 1 causes call to fishpack
      integer ndebug
      parameter (ndebug=0)
      integer interphi
      parameter (interphi = 0)  ! 0 = interpolate rho, solve phi
                                ! 1 = solve phi, interpolate phi
      integer ireverse
c      parameter (ireverse = 0)    ! =1 to swap density regions
      integer ntest
      parameter (ntest=0)  ! =2 stops code somewhere....
      integer idphi
      parameter (idphi=0)  ! =1 to dump phi from fishpack and quit

      integer ibg ! flag for ground potential to be set above (+1)
                  ! or below (-1) physical ground
      parameter (ibg=-1)


c
c  vars for stand-alone call
c
      integer ixst,jyst,kzst
      character(len=80) outname
      integer nkzmn(nzs)

      integer nxb,nyb,nzbt
      integer jslab
c use dx = dz = ds 'cause it's easier that way       
      real ds,ds2,ds3,dls3,dsfac
      real eps0,rjac,pi,omega
c      real  mask(nnx,nny,nnzt)
      real,allocatable ::  mask(:,:,:)
      real,allocatable ::  phip(:,:,:)   !   (nnx,nny,0:nnzt)
      integer i,j,k,i1,j1,k1,n,ix,jy,kz,nstop,nstep,time,m
      real    time_real
      integer iter,maxiter
      integer idif, jdif
      integer i2,j2,k2,i3,j3,k3
      real r2
      integer ibc  !  value for SOR bound. conds: 1=Neum. 0=Diri.
      integer itype
      real resid
      real temp
      integer nend    ! 1 or 2 for uni or bidirectional
      parameter(nend=2)
      integer nxyz,nxyz2
      parameter (nxyz2 = 10000/nend)
      integer chanijk(nxyz2,3,nend) ! i,j,k of channel points
      integer ibonds(nxyz2,3,2,nend) ! list of connected points
      integer linksrc(nxyz2,nend) ! list of source points, i.e., the bond number from which this one branched.
      integer linkstep(nxyz2,nend) ! list of source points, i.e., the bond number from which this one branched.
      real    phichan(nxyz2,nend)   ! potential of channel point
      real signchan(2)           ! -1.0 or +1.0 depending on direction
      integer numchan(2)         ! number of points in channel
      integer npbond(nxyz2,nend)    ! number of potential bonds from point nxyz
c
c with diagonals, have 9+8+8=25 possible unbonded neighbors (at least
c  one neighbor is already bonded).  So change 5 to 25
c
      integer bondijk(25,nxyz2,3,nend) ! locations of pot. bond points
c      real  bondphi(25,nxyz2,2)     ! phi of potential bond
      real  bondprob(25,nxyz2,nend)    ! probability of bond
      double precision  probtot, probcum
      real  phistart,phistlow,phisthi,phistart0
      real  chandist(nxyz2,nend)
      real  cgchandist(nxyz2,nend)
      integer istop(2)  ! flag for no propagation on one end or the other
      integer istopcnt(2)
c
c  need arrays for distance between neighboring points
c
      real,allocatable :: gd(:,:,:,:,:)
c      real gd(-1:1,-1:1,-1:1,0:nnzt,2)  !  1: distance, 2: inverse of distance
      
      real ran0   ! random number function
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
      
c      real tem1(nnx,nny,nnzt)
       real emax,emax2,emagmax,ezmax,ehmax
       real dx,dy,dz,dzi
       real x,y,z

      real dzl,dyl,dxl
      real tdx,tdy,tdz
!      integer neelec
      real pot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
      real exdg(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
      real eydg(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
      real ezdg(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
      real emagdg(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) !(-nor:nx+nor,-nor:ny+nor,-nor:nz+nor,neelec)
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

c      character*40 fname   ! input filename with charge density field
c      character*40 outfile  ! name for binary output
c      character*40 tfile    ! name for text file with link points
      character*6  rstime
      character*2  nmliter
c      integer idump   ! input flag for binary output
      real eps1
      parameter (eps1 = 1.e-2)
      real facx,facy,facz
      real scint,potint
      integer ifile,istat
c
c  electrical permitivity of air C / (N m**2) -  check the units
c
      real eperao
      parameter (eperao  = 8.8592e-12 )
c
c  arrays used in FISHPAK for solver for electric potential field
c
      integer id1,jd1,kd1
      parameter (id1=0,jd1=0,kd1=0)

      integer nwslv

      integer nzb

c other stuff for poisson solver
      integer ierror
      real xs, xf, ys, yf, zs ! , zf  ! domain limits in meters
      real pertrb,elmbda

      integer ist,jst,kst,ist1  ! endpoint of needle from groundplane
      integer istride,jstride,kstride
      integer jend
      integer isp,jsp
c
c stuff for charge removed
c

      real qtot(2)
      real rhotemp
      real rhotot(2),rhotot1(2)
      logical ltest
      real arhopos, arhoneg, rhopos1, rhopos2, rhoneg1, rhoneg2  
      
c vars for localsor
      integer i0,j0,k0
      integer iloc0(nend),jloc0(nend),kloc0(nend)
      integer  mdelf
      real omegal
      integer imn,jmn,kmn,imx,jmx,kmx
      integer imn1,jmn1,kmn1,imx1,jmx1,kmx1
      integer imn2,jmn2,kmn2,imx2,jmx2,kmx2
      integer i1mx, k1mx, i0mx, k0mx
      
c
c vars for allowing for breakeven field threshold
c
!      integer ibrkd  ! = 1 for constant breakdown threshold
                     ! = 2 for breakeven field
      real ebrkd(nzs)  !  this array is passed in -- don't change it!
      real ebrkds(nzs)
      real ecritz(nzs)
      real zlev(nzs)
      
       integer ninit  ! number of possible initiation points
       integer npick1
       integer,allocatable :: ipick(:,:)
c       integer ipick(nxs*nys*nzs/10,3)
       integer ichoose
       integer idum, mseed
       data mseed /678/
       real etest
!       real fecrit
       real fecritn(nend)  ! 
!       real fefac(nend)
       real fecinc   ! increment for adjusting fecrit
       real fecmin   ! minimum value for fecrit
c
c timing vars
c
      real tim1,tim2,timarr(2)
      real second,commas_dtime
      logical lsgi
      parameter( lsgi = .true. )
      
      real times
      real volbi,vol
c
c Misc. vars
c
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
            
c #################################################      
c Begin Executable Code      
c #################################################      


      IF ( lightintz .gt. 1 ) igrid = 1
      
      
      allocate ( hgt(nnz) )
      allocate ( hgtdg(nz) )
      allocate ( zmap(nnz) )
      allocate ( cofx(nnx,3) )
      IF ( ny .gt. 2 ) THEN
      allocate ( cofy(nny,3) )
      ELSE
      allocate ( cofy(-2:ny+2,3) )
      ENDIF
      allocate ( cofz(0:nnzt,3) )
      allocate ( cz(0:nnz-1) )
      allocate ( czz(0:nnz-1) )
      allocate ( gd(-1:1,-1:1,-1:1,0:nnzt,2) )

       npick1 = nxs*nys*nzs/10
       
       IF ( ny .eq. 2 ) npick1 = nxs*nzs

      allocate ( ipick(npick1,3) )

       IF ( ibal .eq. 3 .or. ibal .ge. 4) THEN
          fecmin = fecrit/2.0
       ELSE
          fecmin = 0.5
       ENDIF
       
       allocate ( mask(-1:nnx+2,-1:nny+2,nnzt) )
       allocate ( phip(-1:nnx+2,-1:nny+2,0:nnzt+1) )
       

c
c value of ds is now passed in!
c
      ds  = dslight
      dsfac = 1.0
      IF ( ny .eq. 2 ) dsfac = lightintx
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

      dzi = 1.0/dz
      dxinv2 = 0.5/dx
      dyinv2 = 0.5/dy
      dzinv2 = 0.5/dz

      dxinv = 1./dx
      dyinv = 1./dy
      dzinv = 1./dz
      
      DO kz = 1,nnz
        IF ( istretch ) THEN
          hgt(kz) = z1d2lgt(kz,1) 
          zmap(kz) = 1.0/z1d2lgt(kz,3)
        ELSE
          hgt(kz) = 0.5*dz*(dz1/dz) + (kz - 1)*dz
          zmap(kz) = 1.0
        ENDIF
      ENDDO
      
      DO kz = 1,nz
        hgtdg(kz) = z1d2(kz,1)
      ENDDO
      
       DO i=1,nnx
         cofx(i,1) = 1.0/ds2
         cofx(i,2) = 1.0/ds2
         cofx(i,3) = -(cofx(i,1)+cofx(i,2))
       ENDDO

       DO j=1,nny
         cofy(j,1) = 1.0/(dsfac*ds2)
         cofy(j,2) = 1.0/(dsfac*ds2)
         cofy(j,3) = -(cofy(j,1)+cofy(j,2))
       ENDDO
       IF ( ny .le. 2 ) THEN
         cofy(:,:) = 0.0
       ENDIF
       
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
      
      IF ( ny == 2 ) rjac = 0.0
c      rjac = 0.745356
      omega = 2.00/(1.00 + Sqrt(1.00 - rjac**2))
      IF (ndebug .ge. 2) write(iunit,*) 'pi,rjac,r2',pi,rjac,r2,omega
      

       mstop = -1
       icheck = 1

      IF (ndebug .ge. 1) THEN
        write(iunit,*) 'eta = ',ltgeta
        write(iunit,*) 'intfg = ',intfg
        write(iunit,*) 'ds = ',ds
      END IF
      
      IF (mstop .ge. 0) THEN
        write(iunit,'(a)') 'mstop .ge. 0, so setting nstop'
        nstop = Min(mstop,nstop)
      ENDIF   !  
      IF (icheck .ne. 0) icheck = 1

      

       fecritn(1) = fecrit
       fecritn(2) = fecrit

c
c      IF (ibrkd .le. 2) THEN
c       DO kz = 1,nzs
c
cc      eintz(kz) = feint*ebrkds(kz)
cc       eintz(kz) = 5000.0
cc      ecritz(kz) = ebrkds(kz)
c       fecritn(1) = fecrit
c       fecritn(2) = fecrit
c       ENDDO
c      ELSE
c        DO kz=1,nzs
cc        ebrkds(kz) = 150.0e3
cc        fecrit = ecrit/ebrkds(1)
c        fecritn(1) = fecrit
c        fecritn(2) = fecrit
cc        ecritz(kz) = ebrkds(kz)
cc        eintz(kz) = feint*ecrit
c        END DO
c      END IF  !  (ibrkd .le. 2)


      IF ( istretch ) THEN
      
      DO kz=0,nnzt
!      IF (ndebug .ge. 1) write(0,*) 'istretch: kz = ',kz,z1d2lgt(1,3)
      DO k=-1,1
      DO j=-1,1
      DO i=-1,1
        IF (i.eq.0 .and. j.eq.0 .and. k.eq.0) THEN
         gd(i,j,k,kz,1) = 0.0
         gd(i,j,k,kz,2) = 0.0
         cycle
        ENDIF

        gd(i,j,k,kz,1) = 
     :     Sqrt(dx**2*Abs(i) + dy**2*Abs(j) + 
     :           (dz/z1d2lgt(Max(1,kz+Max(0,k)),3))**2*Abs(k))

          gd(i,j,k,kz,2) = 1.0/(gd(i,j,k,kz,1))

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
        gd(i,j,k,kz,1) = 
     :     Sqrt(dx**2*Abs(i) + dy**2*Abs(j) +  dz**2*Abs(k) )
        IF (.not. (i.eq.0 .and. j.eq.0 .and. k.eq.0)) THEN
          gd(i,j,k,kz,2) = 1.0/gd(i,j,k,kz,1)
        ELSE 
          gd(i,j,k,kz,2) = 0.0
        END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      ENDIF ! istretch
      
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
      
c set up charge distribution

c  45  CONTINUE
c IF THIS IS NOT THE FIRST TIME (stand-alone), THEN 
c   SKIP OVER THE CHARGE/POTENTIAL SETUP, but still reset the mask and tem1

      
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

       write(iunit,*) '*********************************************'
       write(iunit,*) '****** iliter = ',iliter,'*************'
       write(iunit,*) '*********************************************'

      ENDIF ! }

      IF ( iliter .gt. 2 .and. isa .ge. 1 ) THEN
c           tem1(i,j,k) = 0.0
       DO k = 1,nnzt
         DO j = -1,nny+2
          DO i = -1,nnx+2
           mask(i,j,k) = 1.0 
          END DO
         END DO
        END DO
       
       IF ( istop(1) .eq. 2 ) THEN
         numcgn = numcgn + 1
       ELSEIF ( istop(2) .eq. 2 ) THEN
         numcgp = numcgp + 1
       ELSE
         numic = numic + 1
       ENDIF
       
       
       GOTO 55

      ENDIF ! ( iliter .gt. 2 .and. isa .ge. 1 )

c set up charge distribution
      IF (ndebug .ge. 1) write(iunit,*) 'setup charge'
c C$DOACROSS LOCAL(i,j,k)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)

        DO k = 1,nnzt
         DO j = -1,nny+2
          DO i = -1,nnx+2
           mask(i,j,k) = 1.0 
          END DO
         END DO
        END DO

      DO k = 2,nnzt
      DO j = 1,nny
      DO i = 1,nnx
        rho(i,j,k) = rhofairlgt(k-1) ! 0.0
        rhs(i,j,k) = -rhofairlgt(k-1)/eperao ! 0.0
        phi(i,j,k) = potfairlgt(k-1)
      END DO
      END DO
      END DO


      DO k = 1,1
      DO j = 1,nny
      DO i = 1,nnx
        rho(i,j,1) = rhofairlgt(1) ! 0.0
        rhs(i,j,1) =  0.0
        phi(i,j,1) = 0.0

        phi(i,j,nnz-1:nnz) = potfairlgt(nnz-1) ! pottop ! potfair(nztop)
        rhs(i,j,nnz) = -rhofairlgt(nnz-1)/eperao
      END DO
      END DO
      END DO

       
      

c     subroutine mlint2
c    > (ndebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,
c    >  nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
c
      scintmx = 0.0
      scintmn = 0.0
      ifac = Nint(dxl/ds)
      jend = Max(1, nys-1)
      DO k=1,nzs-1
        IF ( istretch .or. dz .eq. dz1 ) THEN
         IF ( lightintz .eq. 1 ) THEN
          k1 = k
          facz = 0.0
         ELSEIF  ( lightintz .eq. 2 ) THEN
          k1 = (k + 1)/lightintz
         ENDIF
        ELSE
          k1 = Int( ((k-1)*dz + eps1)/dzl ) + 1
          facz = ( (k-1)*dz - (k1-1)*dzl )/dzl
        ENDIF
      DO j=1,jend
        j1 = Int( ((j-1)*dy + eps1)/dyl ) + 1
        facy = ( (j-1)*dy - (j1-1)*dyl )/dyl
        IF ( ny .eq. 2 ) THEN ; facy = 0.0 ; j1 = 1 ; ENDIF
      DO i=1,nxs-1
        i1 = Int( ((i-1)*dx + eps1)/dxl ) + 1
        facx = ( (i-1)*dx - (i1-1)*dxl )/dxl
       IF ( igrid .eq. 0 ) THEN
       IF ( i1 .ge. 1   .and. j1 .ge. 1   .and. k1 .ge.  1 .and.
     :     i1 .le. nx-1 .and. j1 .le. ny-1 .and. k1 .le. nz-1 ) THEN
        call mlint2
     > (0, t0,  -nor+1, nx+nor, -nor+1, ny+nor, -nor+1, nz+nor, 
     >  1, 1, nx, ny, nz, scint, facx, facy, facz,i1,j1,k1, 1)
c        call mlint2
c     > (0, pot,    1, nx,     1, ny,  1, nz, 
c     >  1, 1, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,1)
        call mlint2
     > (0, pot,   -nor+1, nx+nor, -nor+1, ny+nor, -nor+1, nz+nor,
     >  1, 1, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,1)
c       scint = t0(i1,j1,k1)
c       potint = pot(i1,j1,k1)
       rhow(i,j,k) = scint
       phi0(i,j,k) = potint
       scintmx = Max(scintmx,scint)
       scintmn = Min(scintmn,scint)
       ENDIF
       ELSEIF ( igrid .eq. 1 ) THEN
       IF ( i1 .ge. 1   .and. j1 .ge. 1   .and. k1 .ge.  1 .and.
     :     i1 .le. nx-1 .and. j1 .le. ny-1 .and. k1 .le. nz-1 ) THEN
c        call mlint2
c     > (0, t0,  -nor, nx+nor, -nor, ny+nor, -nor, nz+nor, 
c     >  1, 1, nx, ny, nz, scint, facx, facy, facz,i1,j1,k1, 1)
c        call mlint2
c     > (0, elec,  -nor, nx+nor, -nor, ny+nor, -nor, nz+nor, 
c     >  1, neelec, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,ipot)
       scint = t0(i1,j1,k1)
        call mlint2
     > (0, pot,   -nor+1, nx+nor, -nor+1, ny+nor, -nor+1, nz+nor,
     >  1, 1, nx, ny, nz, potint, facx, facy, facz,i1,j1,k1,1)

c       potint = pot(i1,j1,k1)
       rhow(i,j,k) = scint
       phi0(i,j,k) = potint
       IF ( ny == 2 .and. scint .gt. scintmx ) THEN
         i1mx = i1
         k1mx = k1
         i0mx = i
         k0mx = k
       ENDIF
       
       scintmx = Max(scintmx,scint)
       scintmn = Min(scintmn,scint)
       ENDIF
       ENDIF
      END DO ! i
      END DO ! j
        
        ecritz(k) = ebrkds(k1) + facz/dz*(ebrkds(k1+1) - ebrkds(k1))
        IF ( ieint .eq. 1 ) THEN
          eintz(k) = eint
        ELSE
          eintz(k) = feint*ecritz(k)
        ENDIF
        
      END DO ! k
      
c      IF ( isa .lt. 1 ) THEN
        write(iunit,*) 'scintmx,scintmn = ', scintmx, scintmn 
c      ENDIF
            
      IF (ndebug .ge. 1) write(iunit,*) 'setup forcing'
c C$DOACROSS LOCAL(i,j,k)
      jend = Max(1, nys-1)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      DO k = 1,nzs-1
      DO j = 1,jend
      DO i = 1,nxs-1 
        rhs0(i,j,k)  = -rhow(i,j,k)/eperao
      END DO
      END DO
      END DO
      
      IF ( bcx .eq. 2 ) rhs0(nxs,1:jend,1:nzs-1) = 
     :                          rhs0(1,1:jend,1:nzs-1)
      IF ( bcy .eq. 2 ) rhs0(1:nxs-1,nys,1:nzs-1) = 
     :                          rhs0(1:nxs-1,1,1:nzs-1)
      
!      tim1 = commas_dtime(timarr)  ! for SGI
!      IF ( .not. lsgi ) tim1 = second()  !  for cray

c do fft solve on work grid    
      IF (ndebug .ge. 1) write(iunit,*) 'call mudpack, ny = ', ny 

      IF ( ny .gt. 2 ) THEN ! 3D case
      
      intl = 0
c
c     set boundary condition flags
c
      nxa = 2
      nxb1 = 2

      IF ( bcx .eq. 2 ) THEN
       nxa = 0
       nxb1 = 0
      ENDIF

      nyc = 2
      nyd = 2

      IF ( bcy .eq. 2 ) THEN
       nyc = 0
       nyd = 0
      ENDIF

      nze = 1
      nzf = 1 ! changed 9/5/2006 to be consistent with settings in mud1

c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      kzr = kkzr
      iex = iiex
      jey = jjey
      kez = kkez
      nxm = nnx
      nym = nny
      nzm = nnz
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set method of relaxation--point gauss/seidel only for mud3sp
c
      method = 0
c
c     set default multigrid options
c
      mgopt(1) = 2  ! kcycle (default 2)
!      IF ( igrid .eq. 1 ) THEN
!       mgopt(1) = 1
!      ENDIF
      mgopt(2) = 1  ! iprer  (default 2)
      mgopt(3) = 1  ! ipost  (default 1)
      mgopt(4) = 1  ! intpol (default 3)
c
c     set full multigrid cycling by flagging no initial guess
c
      iguess = 0
c
c     set a limit of three w(2,1) cycles from the finest level
c
      maxcy = 2
!      IF ( igrid .eq. 1 ) THEN
!       maxcy = 1
!      ENDIF
c
c     set mesh increments
c
c     set end points of solution cube in (x,y,z) space
c
      xa = 0.
      xb = (nnx-1)*dx
      yc0 = 0.
      yd = (nny-1)*dy
      IF ( istretch ) THEN
        ze = -z1d2lgt(1,1)
        zf = ze + (nnz-1)*dlzlgt
      ELSEIF ( dz1 .eq. dz ) THEN
        ze = -dz/2.0
        zf = ze + (nnz-1)*dz
      ELSE
        ze = 0.0
        zf = ze + (nnz-1)*dz
      ENDIF
c
c     set for no error control (this allows phi to be equivalenced with work)
c
      tolmax = 0.0
c
c
c     print input arguments
c
      IF ( ndebug .ge. 2 ) THEN
      write(iunit,100)intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +       nxm,nym,nzm,iguess,maxcy,method,nwork,xa,xb,yc0,yd,ze,zf,
     +       tolmax,mgopt(1),mgopt(2),mgopt(3),mgopt(4)
  100 format(/' input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' nze = ',i2, ' nzf = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' kzr = ',i2,
     +/' iex = ',i2, ' jey = ',i2, ' kez = ',i2,
     +/' nxm = ',i4,' nym = ',i4,' nzm = ',i4, ' iguess = ',i2,
     +/' maxcy = ',i2,
     +/' method = ',i2, ' work space length input = ',i12,
     +/' xa = ',f7.2,' xb = ',f12.2,
     +/' yc = ',f7.2,' yd = ',f12.2,
     +/' ze = ',f7.2,' zf = ',f12.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' kcycle = ',i2
     +/' iprer = ',i2
     +/' ipost = ',i2
     +/' intpol = ',i2 )
     
      write(iunit,*) '(zf-ze)/float(nzm-1) = ',(zf-ze)/float(nzm-1)
      

        write(iunit,104) intl
      ENDIF

      tim1 = commas_dtime(timarr) 

  104 format(/' discretization call to mud3sp', ' intl = ', i2)
      IF ( istretch ) THEN
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzlgt,bndc,rhs,phi,mgopt,ierror)
      ELSE
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzu,bndc,rhs,phi,mgopt,ierror)
      ENDIF
      IF (ndebug .ge. 1) THEN
        write (iunit,105) ierror,iprm(21)
        write (iunit,*) 'iprm(21) = ',iprm(21), nwork
      ENDIF

  105 format(' ierror = ',i2, ' minimum work space = ',i12)
      
      IF (ierror.gt.0) THEN
       call exit(0)
      ENDIF
c
c     approximate pde
c
      intl = 1
      IF (ndebug .ge. 1) THEN
        write(iunit,106) intl,method,iguess,maxcy
      ENDIF
  106 format(/' approximation call to mud3sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      IF ( istretch ) THEN
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzlgt,bndc,rhs,phi,mgopt,ierror)
      ELSE
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzu,bndc,rhs,phi,mgopt,ierror)
      ENDIF

      IF (ndebug .ge. 1) THEN
        write (iunit,107) ierror
      ENDIF
  107 format(' ierror = ',i2)

      if (ierror.gt.0) call exit(0)


      
c      write(iunit,*) 'got here'
      
      tim2 = 0.0
      tim2 = commas_dtime(timarr)
!      IF ( lsgi ) tim2 = commas_dtime(timarr)
!      IF ( .not. lsgi ) tim2 = second()
      
      write(iunit,*) 'time for mudpack = ',tim2

!      IF ( lsgi ) write(iunit,*) 'time for mudpack = ',tim2,timarr(1)
!      IF ( .not. lsgi ) write(iunit,*) 'time for mudpack = ',tim2 - tim1
      
      
      ELSE ! 2D code

c
c     set input integer parameters
c
      intl2d = 0
c
c     set boundary condition flags
c
      IF (ndebug .ge. 1) write(0,*) 'start of mud2d1'
      IF (ndebug .ge. 1) write(iunit,*) 'start of mud2d1'
      

      IF ( ndebug .ge. 2 ) THEN
      j = Max(1, nys/2)
      write(iunit,*) 'interpolated potential'
      write(iunit,*) 'k,z,phi,rho through nx/2,ny/2, i,j=',nxs/2,j
      DO k=nnzt-1,2,-1 
          
        write(iunit,*) k,z1d2lgt(k,1),phi0(nxs/2,j,k),
     :   -0.5/dlzlgt*(phi0(nxs/2,j,k+1)-phi0(nxs/2,j,k-1))*
     :     z1d2lgt(k,3),
     :   rhow(nxs/2,j,k)
      ENDDO
      ENDIF

      nxa2d = 2
      nxb12d = 2
      IF ( bcx .eq. 2 ) THEN
       nxa2d = 0
       nxb12d = 0
      ENDIF
      nyc2d = 1
      nyd2d = 1
c      nze = 1 
c      nzf = 1  
      IF ( kbdcnd .eq. 2 ) THEN
      nyc2d = 1
      nyd2d = 1
      ELSEIF ( kbdcnd .eq. 1 ) THEN
      nyc2d = 1
      nyd2d = 1  ! changed from 2 to 1 on 12/31/2004 by erm 
               ! solution seems to be more stable as a result, especially
               ! when there is an upper level grid stretch.
               ! Also seems to preserve net charge better than nzf = 2
      ELSEIF ( kbdcnd .eq. 3 ) THEN
      nyc2d = 2
      nyd2d = 2
      ENDIF
c
c     set grid sizes from parameter statements
c
      ixp2d = iixp
      jyq2d = kkzr! jjyq
      iex2d = iiex
      jey2d = kkez ! jjey
      nxm2d = nnx
      nym2d = nnz
c      nzm = nnz
c      nz1d = nnz
c
c     set work space length approximation from parameter statement
c
      nwork2d = llwork
c
c     set method of relaxation--point gauss/seidel only for mud3sp
c
      method2d = 0
c
c     set default multigrid options
c
      mgopt(1) = 2  ! kcycle (default 2)
      mgopt(2) = 2  ! iprer  (default 2)
      mgopt(3) = 1  ! ipost  (default 1)
      mgopt(4) = 3  ! intpol (default 3)
c
c     set full multigrid cycling by flagging no initial guess
c
      iguess2d = 0
c
c     set a limit of three w(2,1) cycles from the finest level
c
      maxcy2d = 4
c
c     set end points of solution cube in (x,y,z) space
c
      xa2d = 0.
      xb2d = (nxm2d-1)*dx
c      yc = 0.
c      yd = (nym-1)*dy
      yc2d = -dz*0.5
      yd2d = yc2d + (nym2d-1)*dz
      
      IF ( ndebug .ge. 1 ) THEN
        write(iunit,*) 'spark2d: xa2d,xb2d,yc2d,yd2d = ',xa2d,xb2d,yc2d,yd2d
        write(iunit,*) 'dx,dz,dlzlgt,nxm2d,nym2d = ',dx,dz,dlzlgt,nxm2d,nym2d
      ENDIF

c      STOP
c
c     set for no error control (this allows phi to be equivalenced with work)
c
      tolmax2d = 0.0
      
      tim1 = commas_dtime(timarr) 


      IF ( .true. ) THEN

      ierror = 0
      IF (ndebug .ge. 1) THEN
        write(*,204) intl2d
        DO i = 1,16
          write(0,*) 'i, iprm2d(i) = ',i,iprm2d(i)
        ENDDO
        write(0,*) intl2d,nxa2d,nxb12d,nyc2d,nyd2d,ixp2d,jyq2d,iex2d,jey2d,
     +              nxm2d,nym2d,
     +              iguess2d,maxcy2d,method2d,nwork2d,lwrkqd2d,itero2d
        write(0,*) 'fprm2d = ', xa2d,xb2d,yc2d,yd2d,tolmax2d,relmax2d
      ENDIF
  204 format(/' discretization call to mud2sp', ' intl = ', i2)
      IF ( istretch ) THEN
      call mud2sp(iprm2d,fprm2d,work,cfx,cfzlgt,bndc2d,rhs,phi,mgopt,ierror)
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
      ELSE
c      call mud2sp(iprm2d,fprm2d,work,cfx,cfzu,bndc2d,rhs,phi,mgopt,ierror)
      ENDIF
      IF (ndebug .ge. 1) THEN
        write (*,205) ierror,iprm2d(16)
        write (0,*) 'ierror,iprm2d(16),llwork = ',ierror,iprm2d(16),llwork
        write(0,*) 'nwork2d,iprm2d(15) = ',nwork2d,iprm2d(15)
        DO i = 1,16
          write(0,*) 'i, iprm2d(i) = ',i,iprm2d(i)
        ENDDO
      ENDIF

  205 format(' ierror = ',i2, ' minimum work space = ',i7)
      
      IF (ierror.gt.0) THEN
       call exit(0)
      ENDIF
c
c     approximate pde
c
      intl2d = 1
      IF (ndebug .ge. 1) THEN
        write(*,206) intl2d,method2d,iguess2d,maxcy2d
      ENDIF
  206 format(/' approximation call to mud2sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      IF ( istretch ) THEN
      call mud2sp(iprm2d,fprm2d,work,cfx,cfzlgt,bndc2d,rhs,phi,mgopt,ierror)
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
      ELSE
c      call mud2sp(iprm2d,fprm2d,work,cfx,cfzu,bndc2d,rhs,phi,mgopt,ierror)
      ENDIF

      IF (ndebug .ge. 1) THEN
        write (*,207) ierror
      ENDIF
  207 format(' ierror = ',i2)

      if (ierror.gt.0) call exit(0)

         itype = 1
         i0 = nnx/2
         k0 = nnz/2
         j0 = 1

C do not call here!!! makes no difference: phip gets overwritten         
c        DO k = 1,lightintz
c         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,nnx/2,
c     :      cofx,cofy,cofz, bcx, bcy,
c     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,0.2*soreps,ndebug)
c        ENDDO
      
      tim2 = 0.0
      tim2 = commas_dtime(timarr)
      
      write(iunit,*) 'time for mudpack2D = ',tim2

      ENDIF ! false
      ENDIF ! ny .gt. 2
      
  55  CONTINUE


      IF (ndebug .ge. 1) write(0,*) 'copy phi2 into phip'
      
c      IF ( nny .ge. 2 ) THEN
      DO k = 1,nnzt
      DO j = 1,nny
      DO i = 1,nnx
        phip(i,j,k) =  phi2(i,j,k)
      END DO
      END DO
      END DO
      
c      ELSE
      
c       phip(:,:,:) = 0.0
      
c      ENDIF

      IF (ndebug .ge. 1) write(0,*) 'set phip = 0'

      DO j = 1,nny
      DO i = 1,nnx
        phip(i,j,0) =  0.0
      END DO
      END DO
      
      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF

       IF ( ny .eq. 2 ) THEN
         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,nnx/2,
     :      cofx,cofy,cofz, bcx, bcy,
     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,0.1*soreps,ndebug)
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
        ENDIF
       pstop = 0

      IF ( ndebug .ge. 2 ) THEN
      j = Max(1, nys/2)
      write(iunit,*) 'k,z,phi,rho through nx/2,ny/2, i,j=',nxs/2,j
      DO k=nnzt-1,2,-1
         
        write(iunit,*) k,z1d2lgt(k,1),phi0(nxs/2,j,k),
     :   -0.5/dlzlgt*(phi0(nxs/2,j,k+1)-phi0(nxs/2,j,k-1))*
     :     z1d2lgt(k,3),
     :   rhow(nxs/2,j,k)
      ENDDO

      write(iunit,*)
      
      IF ( ny .gt. 2 ) THEN
      ix = nxs/2
      jy = Min(ny,54)
      i1 = ix
      j1 = jy
      write(iunit,*) 'k,z,phi,rho through ix,jy= ',ix,jy
      DO k=nzs-1,2,-1
        k1 = k
        rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
         k1 = k
         i = 0
         j = 0
         kz = 1
            e = gd(i,j,kz,k1+kz,2)/2.0*( 
     :                - 1.0*phip(i1+2*i,j1+2*j,k1+2*kz) 
     :                + 4.0*phip(i1+i,j1+j,k1+kz) 
     :                - 3.0*phip(i1,j1,k1) )
           e1 = gd(i,j,kz,k1+kz,2)*
     :            (phip(i1+i,j1+j,k1+kz) -  phip(i1,j1,k1))


        write(iunit,*) k,z1d2lgt(k,1),phip(ix,jy,k),
     :   -0.5/dlzlgt*(phip(ix,jy,k+1)-phip(ix,jy,k-1))*
     :     z1d2lgt(k,3),e,e1,
     :   rho(ix,jy,k),rhotemp
      ENDDO
      ENDIF
      
      ENDIF ! ndebug

c       write(iunit,*) 'calling mxfldstrz, dx,dy,dlzlgt=',dx,dy,dlzlgt
       
c      SUBROUTINE MXFLDSTRZ(nx,ny,nz,
c     :    phi,gd,pstop,nbw,nbs,
c     :    emax,emagmax,i1,j1,k1,i2,j2,k2,dx,dy,dz,istretch,
c     :    z1d2,nz1d)

      IF (ndebug .ge. 1) write(0,*) 'call msfldstrz'

       CALL mxfldstrz(nnx,nny,nnzt,
     :    phip,gd,pstop,nbw,nbs,
     :    emax,emagmax,i1,j1,k1,i2,j2,k2,ds,ds,dz,istretch,
     :    z1d2lgt,nzmax1lgt,ezmax,iunit,bcx,bcy)
      
      write(iunit,'(a,1pe12.5,a,3(1x,i4))')
     :   'emagmax = ',emagmax,' at ijk ',i1,j1,k1
      write(iunit,'(a,1pe12.4)')
     :   'ezmax = ',ezmax
      write(iunit,'(a,1pe12.5,a,3(1x,i4))')
     :   'emax = ',emax,' at ijk ',i2,j2,k2  
      
      IF ( emagmax .gt. 4.0e6 .or. pstop .eq. 1 ) THEN
        IF ( pstop .eq. 0 ) THEN
           write(iunit,*) 'Efield too large!  Stopping now...'
        ENDIF
        
        write(iunit,*) 'x-line through large field'
        DO ix=1,nnx
          write(iunit,*) ix,k1,phi(ix,j1,k1),rho(ix,j1,k1),
     :     rhs(ix,j1,k1)
        ENDDO
        
        write(iunit,*) 'z-line through large field'
        DO kz=nzs-1,1,-1
          write(iunit,*) kz,phi(i1,j1,kz),rho(i1,j1,kz),
     :     rhs(i1,j1,kz)
        ENDDO
        i1=Max(1,nnx/2)
        j1=Max(1,nny/2)
        DO kz=nzs-1,1,-1
          write(iunit,*) kz,phi(i1,j1,kz),rho(i1,j1,kz),
     :     rhs(i1,j1,kz)
        ENDDO
       
        STOP
      ENDIF
      

c make list of possible initiation points
      ninit = 0
      
      useinterp = .true.
c      IF ( ny .eq. 2 ) useinterp = .false.
      
      IF ( useinterp ) THEN !  if useinterp = true, then find init points from the FINE grid
      
      IF ( icldtop .ne. 0 ) THEN
        DO kz=1,nzs
          IF ( hgt(kz) .lt. zcldtop ) THEN
            k2 = kz
            icldtop = kz
          ENDIF
        ENDDO
        k2 = Min( k2, nzs-4 )
c        k2 = Min( icldtop, nzs-4 )
      ELSE
        k2 = nzs-4
      ENDIF
      
      IF ( iliter .lt. 40 .or. isa .ge. 1 ) THEN
      
      ist = Max(1,nbw+1)
      jst = Max(1,nbs+1)
      jend = Max(1, nny-nbn-1)
      
      istride = Min(2,lightintx)
      jstride = Min(ny-1,istride)
      kstride = Min(2,lightintz)

      DO kz = 3, k2, kstride ! nzs-4
       emagmax = 0
       ezmax = 0
       ehmax = 0
       x = 0
       y = 0
      DO jy = jst,jend, jstride
      DO ix = ist,nnx-nbe-1, istride
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
!          ex1 = (phip(ix-1,jy  ,kz) - phip(ix,jy,kz))*dxinv
!          ey1 = (phip(ix  ,jy-1,kz) - phip(ix,jy,kz))*dyinv
         IF ( istretch ) THEN
           ez  = (phip(ix,jy,kz-1) - phip(ix,jy,kz))*dzinv*z1d2lgt(kz,3)
!           ez1 = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2*z1d2lgt(kz,3)
!           ez1 = ez
         ELSE
           ez = (phip(ix,jy,kz-1) - phip(ix,jy,kz+1))*dzinv2
!           ez1 = (phip(ix,jy,kz-1 ) - phip(ix,jy,kz))*dzinv
         ENDIF
         emag  =  Sqrt(ex**2 + ey**2 + ez**2) 
!         emag1 = Sqrt(ex1**2 + ey1**2 + ez1**2)
!         emag  = Max( emag1, emag )
c         IF ( emag .gt. emagmax ) THEN
c           emagmax = emag
c           ezmax = ez
c           ehmax = Sqrt(ex**2 + ey**2)
c           x = (ds*Real(ix-nbw-1 - 1))/dx1 + 1
c           y =  (ds*Real(jy-nbs-1 - 1))/dy1 + 1
c         ENDIF

       IF (emag .gt. overvolt*elgtfdel*ecritz(kz) ) THEN
!        write(iunit,*) 'kz, emag, elgtfdel,ectritz(kz) = ',kz,emag, elgtfdel,ecritz(kz) 
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(iunit,*) 'Ee-gad! Over-ran the ipick array! ABORT!',npick1,ninit
        STOP
      END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
c       IF ( iliter .eq. 1 ) THEN
c        write(iunit,'(2(f4.1,1x),i4,2x,3(1x,1pe12.4))') 
c     :   x,y,kz,emagmax,ezmax,ehmax
c       ENDIF
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

       IF (emag .gt. overvolt*ecritz(kz) ) THEN
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(iunit,'(a)') 'Ee-gad! Over-ran the ipick array! ABORT!'
        STOP
      END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
      ENDDO  ! k loop
      END IF
 
      ELSE ! if useinterp = false, then find init points from the COARSE grid

     
      IF ( icldtop .ne. 0 ) THEN
        DO kz=1,nz
          IF ( hgtdg(kz) .lt. zcldtop ) THEN
            k2 = kz
            icldtop = kz
          ENDIF
!          write(iunit,*) 'kz,hgtdg,zcldtop,k2 = ',kz,hgtdg(kz),zcldtop,k2
        ENDDO
        k2 = Min( k2, nz-4 )
c        k2 = Min( icldtop, nzs-4 )
      ELSE
        k2 = nz-4
      ENDIF
      
      IF ( .true. .or. iliter .lt. 40 .or. isa .ge. 1 ) THEN
      
      ist = 2 ! Max(1,nbw+1)
      jst = Max(1,nbs+1)
      jend = Max(1, ny-1)

      DO kz = 3, k2 ! nzs-4
       emagmax = 0
       ezmax = 0
       ehmax = 0
       x = 0
       y = 0
      DO jy = jst,jend
      DO ix = ist,nx-2
      
         ex = exdg(ix,jy,kz)
         ey = eydg(ix,jy,kz)
         ez = ezdg(ix,jy,kz)
         
         emag  = emagdg(ix,jy,kz) !  Sqrt(ex**2 + ey**2 + ez**2) 

!        write(iunit,*) 'kz, emag, elgtfdel,ectritz(kz) = ',kz,emag, elgtfdel,ebrkds(kz) 
       IF (emag .gt. overvolt*elgtfdel*ebrkds(kz) ) THEN
!        write(iunit,*) 'kz, emag, elgtfdel,ectritz(kz) = ',kz,emag, elgtfdel,ebrkds(kz) 
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(iunit,*) 'Ee-gad! Over-ran the ipick array! ABORT!',npick1,ninit
        STOP
      END IF
        ipick(ninit,1) = (ix-1)*lightintx + 1 + nbw
        ipick(ninit,2) = Max(1, (jy-1)*lightintx + 1 + nbs )
        IF ( ny .eq. 2 ) ipick(ninit,2) = 1
        ipick(ninit,3) = (kz-1)*lightintz + 1
!        write(iunit,*) 'added init point ',ix,jy,kz,ipick(ninit,1),ipick(ninit,2),ipick(ninit,3)
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
c       IF ( iliter .eq. 1 ) THEN
c        write(iunit,'(2(f4.1,1x),i4,2x,3(1x,1pe12.4))') 
c     :   x,y,kz,emagmax,ezmax,ehmax
c       ENDIF
      ENDDO  ! k loop
      
      IF ( ninit .eq. 0 ) STOP
      
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

       IF (emag .gt. overvolt*ebrkds(kz) ) THEN
        ninit = ninit + 1
      IF (ninit .gt. npick1) THEN
        write(iunit,'(a)') 'Ee-gad! Over-ran the ipick array! ABORT!'
        STOP
      END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
      ENDDO  ! k loop
      END IF

      ENDIF ! useinterp
      
      write(iunit,'(1x,i6,a)') ninit, ' possible starting points'
       
       numtry = 0
 199  CONTINUE
       numtry = numtry + 1
       IF ( numtry .gt. 1 ) THEN
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
c clear out mask array 
           DO iend = 1, nend
             DO n=1,numchan(iend)
               i1 = chanijk(n,1,iend)
               j1 = chanijk(n,2,iend)
               k1 = chanijk(n,3,iend)
               mask(i1,j1,k1 ) = 1.0
               phip(i1,j1,k1) =  phi2(i1,j1,k1) 
               IF ( numchan(iend) .le. 5 ) THEN
                 write(iunit,'(5i5)') iend,n, i1,j1,k1
               ENDIF
             END DO
           END DO
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

            mask(:,:,: ) = 1.0

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
           write(iunit,'(a,i4)') 'Try number ',numtry
       END IF
       IF ( ninit .le. 0 ) THEN
         write(iunit,'(a)') 'No more points to choose!'
         lgtstp = 1
         GOTO 999
       END IF
c      
c here either choose random point (irand = 0) or use
c a specified starting point (irand = 1)
c
      IF ( irand .ne. 1 .or. isa .lt. 1 ) THEN
      ichoose = int(real(ninit)*ran0(iseed)) + 1
        IF ( irand .eq. -1 ) ichoose = Min(iliter,ninit)
      write (iunit,'(2(a,i6))') 'ichoose = ',ichoose
      write(iunit,*) 'iseed = ',iseed
      i1 = ipick(ichoose,1)
      j1 = ipick(ichoose,2)
      k1 = ipick(ichoose,3)
      ixst = i1
      jyst = j1
      kzst = k1
      write(iunit,'(a,3(1x,i4))') 'Picked starting point ',i1,j1,k1
c update the points array (move up the points lower on the list)     
      DO idum=ichoose,ninit-1
        ipick(idum,1) = ipick(idum+1,1)
        ipick(idum,2) = ipick(idum+1,2)
        ipick(idum,3) = ipick(idum+1,3)
      END DO
      ELSE
      write(iunit,*) 'Starting at predetermined point: ',ixst,jyst,kzst
        i1 = ixst
        j1 = jyst
        k1 = kzst
      ENDIF  ! ( irand .ne. 1 .or. isa .ne. 1 )
      
c one fewer init point now      
      ninit = ninit - 1
c find neighboring point with largest E and make it the first bond.
      emax = 0.0
      ist = -1 ; isp = 1
      jst = -1 ; jsp = 1
       IF ( ny .eq. 2 ) THEN
        jst = 0
        jsp = 0
       ENDIF
c      IF ( bcy .eq. 2 .and. ( j1 .eq. 1     ) ) jst = 0
c      IF ( bcy .eq. 2 .and. ( j1 .eq. nny-1 ) ) jsp = 0
c      IF ( bcx .eq. 2 .and. ( i1 .eq. 1     ) ) ist = 0
c      IF ( bcx .eq. 2 .and. ( i1 .eq. nnx-1 ) ) isp = 0
      DO k= -1,1
      DO j=  jst,jsp
      DO i=  ist,isp
         etest = Abs(phip(i1+i,j1+j,k1+k) - phip(i1,j1,k1))*
     :              gd(i,j,k,k1,2)
c         IF ( k .ne. 0 .and. istretch ) THEN
c          etest = etest*z1d2lgt(k1,1)
c         ENDIF
         IF (etest .gt. emax) THEN
           emax = etest
           i2 = i1+i
           j2 = j1+j
           k2 = k1+k
         END IF
      END DO ! i
      END DO ! j
      END DO ! k

       idif = i2 - i1
       jdif = j2 - j1
       mask(i2,j2,k2) = 0.0

!       IF ( bcy .eq. 2 .and. ( j2 .le. 1 .or. j2 .ge. nny-1) ) 
!     :    print*,'0 j1,j2,nny = ',j1,j2,nny
       IF ( bcy .eq. 2 .and. ( j2 .le. 1 .or. j2 .ge. nny-1 ) ) THEN
 !       IF ( j2 .le. 1 .or. j2 .ge. nny-1 ) print*,'1 j2,nny = ',j2,nny
        IF ( j2 .eq. 0   ) j2 = nny - 1
         IF ( j2 .eq. nny ) j2 = 1
 !       IF ( j2 .le. 1 .or. j2 .ge. nny-1 ) print*,'2 j2,nny = ',j2,nny
       ENDIF

       IF ( bcx .eq. 2 .and. ( i2 .eq. 0 .or. i2 .eq. nnx ) ) THEN
         IF ( i2 .eq. 0   ) i2 = nnx - 1
         IF ( i2 .eq. nnx ) i2 = 1
       ENDIF

      IF ( phip(i2,j2,k2) .gt.  phip(i1,j1,k1) ) THEN
      write(iunit,'(a,3(1x,i4),a)') 'First bond point      ',
     :       i2,j2,k2,' (neg end)'
        IF ( k2 .gt. k1 ) THEN
         idown = 2 ! pos end going down
        ELSE 
         idown = 1 ! neg end going down
        ENDIF
      ELSE
      write(iunit,'(a,3(1x,i4),a)') 'First bond point      ',
     :       i2,j2,k2,' (pos end)'
        IF ( k2 .lt. k1 ) THEN
         idown = 2 ! pos end going down
        ELSE 
         idown = 1 ! neg end going down
        ENDIF
      ENDIF
      
      IF ( k1 .eq. k2 ) idown = 0
      
      IF ( isa .ge. 1 ) THEN
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
        
c start the leader at the location of max efield
      IF (nend .eq. 2) THEN
       phistart = ( phip(i1,j1,k1) + phip(i2,j2,k2) )/2.0
       phidiff = Abs( phip(i1,j1,k1) - phip(i2,j2,k2) )
       write(iunit,'(a,3(1x,1pe13.5))') 'phistart,phi1,phi2 = ',
     :       phistart,phip(i1,j1,k1),phip(i2,j2,k2)
       phistart0 = phistart
c         phistlow = Min( phip(i1,j1,k1),phip(i2,j2,k2) )
         phistlow = phistart - 5.e2*gd(idif,jdif,k2-k1,k2,2)*phidiff
c         phisthi  = Max( phip(i1,j1,k1),phip(i2,j2,k2) )
         phisthi  = phistart + 5.e2*gd(idif,jdif,k2-k1,k2,2)*phidiff
         signchan(1) = -1.0
         signchan(2) =  1.0
       IF (phip(i1,j1,k1) .gt. phistart ) THEN
          m = 1
          n = 2
       ELSE
         m = 2
         n = 1 
       END IF
       numchan(m) = 1
       chandist(1,m) = gd(idif,jdif,k2-k1,k2,1)*0.5
       phichan(1,m) = phistart ! - chandist(1,m)*signchan(m)*eintz(k1)
       phip(i1,j1,k1) = phichan(1,m)

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
       
       numchan(n) = 1
       chandist(1,n) = gd(idif,jdif,k2-k1,k2,1)*0.5
       phichan(1,n) = phistart ! - chandist(1,n)*signchan(n)*eintz(k2)
       phip(i2,j2,k2) = phichan(1,n)
       
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

       mask(i1,j1,k1) = 0.0
       mask(i2,j2,k2) = 0.0
       
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

       imn1 = Min(i1,i2)
       jmn1 = Min(j1,j2)
       kmn1 = Min(k1,k2)
       imx1 = Max(i1,i2)
       jmx1 = Max(j1,j2)
       kmx1 = Max(k1,k2)
      
      END IF ! (nend .eq. 2)
       
      omega=1.45

      i0 = Min(i1,i2)
      j0 = Min(j1,j2)
      k0 = Min(k1,k2) 
      itype = 0
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 1'
      IF ( ny .gt. 2 ) THEN
      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,2*mdel,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ELSE
      CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,2*mdel*lightintz,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ENDIF
      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF

      IF (ntest .eq. 2) STOP
      

c start loop for adding points to channel, relaxing, doing it again
c nstop indicates max number of times to repeat
      IF (ndebug .ge. 1) write(iunit,*) 'start propagation'
      
      DO iend = 1,nend
        istop(iend) = 0
        istopcnt(iend) = 0
      END DO
      
      ihold = 0
      iextra = 0
      
!      DO nstep = 1,nstop
      
      nstep = 0
      
      DO while ( numchan(1) .lt. nxyz2 .and. numchan(2) .lt. nxyz2 )
      
      nstep = nstep + 1
      
      IF (icheck .eq. 2) THEN
        IF (nstep .ge. nstop2) THEN
         write(iunit,*) 'stopping second end now....'
         GOTO 201
        END IF
      END IF
            
      DO iend = 1,nend
c look for possible bonds....
c At 'positive' end (end that is going toward higher potential, which is
c really the negative leader since it carries negative charge), the 
c potential at a neighboring point should be _higher_ in order to be considered
c and vice versa for the 'negative' end (seeks neg. charge and lower potential).
c  So, look for   -signchan(iend)*(bondphi-chanphi) .gt. 0.0, and actually want
c  -signchan(iend)*(bondphi-chanphi)/ds .gt. ecrit.  
c  This should prevent travel in an unfavored direction.
c
c       IF ( bcy .eq. 2 ) THEN
c         mask(imn1:imx1,0,kmn1:kmx1) = mask(imn1:imx1,nys-1,kmn1:kmx1)
c       ENDIF
       IF ( ( istop(iend) .le. 1 ) .and. ihold .ne. iend) THEN
       probtot = 0.0
       DO n=1,numchan(iend)
         npbond(n,iend) = 0
         i1 = chanijk(n,1,iend)
         j1 = chanijk(n,2,iend)
         k1 = chanijk(n,3,iend)
       IF (ndebug .ge. 2) write(iunit,*)
     :     'i1,j1,k1: ',i1,j1,k1
         
         jst = -1 ; jend = 1
         IF ( ny .eq. 2 ) THEN ; jst = 0 ; jend = 0 ; ENDIF
         
         DO k= -1,1
         DO j= jst,jend
         DO i= -1,1
           IF (mask(i1+i,j1+j,k1+k) .gt. 0.5 .and.
     :        (bcx .eq. 2 .or. ( i1+i .gt. nbw+1 .and. i1+i .le. nxs+nbw)) .and.
     :        ( (bcy .eq. 2 .or. ny .eq. 2) .or. ( j1+j .gt. nbs+1 .and. j1+j .le. nys+nbs)) .and.
     :         k1+k .gt. 0 .and. k1+k .lt. nzs-1 )    THEN
c            IF ( ( Abs(rho(i1+i,j1+j,k1+k)) .gt. 1.0e-16 ) .or.
c     :           k1 .lt. 10  ) THEN
              ecomp = fefac(iend)*fecritn(iend)*ecritz(k1)
c            ELSE
c              ecomp = Max(fecritn(iend)*ecritz(k1),
c     :                     Min(125.0e3, 2.0*ecritz(k1) ))
c            ENDIF
            
           IF ( l2nde .ge. 2 .and. 
     :            mask(i1+2*i,j1+2*j,k1+2*k) .gt. 0.5 ) THEN
            e = gd(i,j,k,k1+k,2)/2.0*( 
     :                - 1.0*phip(i1+2*i,j1+2*j,k1+2*k) 
     :                + 4.0*phip(i1+i,j1+j,k1+k) 
     :                - 3.0*phichan(n,iend)  )
            IF ( l2nde .eq. 3 ) THEN
              e1 = gd(i,j,k,k1+k,2)*
     :            (phip(i1+i,j1+j,k1+k) -  phichan(n,iend))
c              e = 0.5*( e + e1 )
            ENDIF
           ELSE
            e = gd(i,j,k,k1+k,2)*
     :            (phip(i1+i,j1+j,k1+k) -  phichan(n,iend))
           ENDIF
            
            temp = Max(0.0d0, 
     :      (  - signchan(iend)*e  - ecomp )  )
     
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
       IF (ndebug .ge. 2) write(iunit,*)
     :     'n,npbond: ',n,npbond(n,iend),phip(i1+i,j1+j,k1+k)
           END IF
         END DO ! i
         END DO ! j
         END DO ! k
         
         IF ( npbond(n,iend) .eq. 0 .and. n .eq. numchan(iend) ) THEN
          npbond(n,iend) = 1
          bondprob(npbond(n,iend),n,iend) = 0.0
          IF ( .false. ) THEN
          write(iunit,*) 'Lightning debug: numchan,n,iend',
     :           numchan(iend),n,iend
          write(iunit,*) 'i1,j1,k1',i1,j1,k1

          DO k=1,numchan(iend)
           write(iunit,'(3(1x,i4))') 
     :      (ibonds(k,i,1,iend), i=1,3)
           write(iunit,'(3(1x,i4),1(1x,1pe13.5))') 
     :      (ibonds(k,i,2,iend), i=1,3),phi(i1,j1,k1)
          ENDDO
          ENDIF
          
         ENDIF
     
       END DO ! n


c calculate probabilities
c       DO n=1,numchan(iend)
c         DO m=1,npbond(n,iend)
c           temp = Max(0.0, 
c     :           Abs(bondphi(m,n,iend) - phichan(n,iend))/ds - ecrit)
c           bondprob(m,n,iend) = (temp)**eta
c           probtot = probtot + bondprob(m,n,iend)
c         END DO ! m
c       END DO ! n
       
       IF (ndebug .ge. 2) write(iunit,'(a,e12.5)') 'probtot = ',probtot
       
       IF (probtot .lt. 1.0e-10) THEN
         IF (ndebug .ge. 2) write(iunit,*) 
     :       'No propagation at end ',iend
         IF ( istop(iend) .ne. 2 ) THEN
           istop(iend) = 1
           istopcnt(iend) = istopcnt(iend) + 1
         ENDIF
           IF ( ibal .eq. 2 .and. istopcnt(iend) .gt. 50*500./ds ) THEN
            write(iunit,"(a,i2,a,i4)") 
     :          "istopcnt(",iend," ) = ",istopcnt(iend)
            GOTO 201
           ENDIF
         IF (icheck .eq. 0) THEN
           icheck = 2
           nstop2 = Min(nstop, nstep + nstep/3)
           write(iunit,*) 'New nstop is ',nstop2
         END IF
          IF (nend .eq. 1) GOTO 201
       ELSE
         IF ( istop(iend) .ne. 2 ) THEN
           istop(iend) = 0
           istopcnt(iend) = 0
         ENDIF
       END IF
       IF (nend .gt. 1) THEN
       IF (istop(1) .ge. 1 .and. istop(2) .ge. 1) THEN
c       IF (istop(1) .ge. 1 .and. istop(2) .ge. 1 .and. 
c     :       ihold .eq. 0) THEN
c         IF (ndebug .ge. 1) 
         write(iunit,*)  'No propagation at either end...stopping loop.',
     :       ' nstep=',nstep
         IF ( ( ((numchan(1) + numchan(2) ) .le. 4 .and. ibal.ne.2)
     :       .or. (numchan(1) + numchan(2) ) .le. 2 ) .and. 
     :       numtry .lt. maxtry  .and. ninit .ge. 1 ) THEN
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
c           bondprob(m,n,iend) = bondprob(m,n,iend)/probtot
           probcum = probcum + bondprob(m,n,iend)
           bondprob(m,n,iend) = probcum
       IF (ndebug .ge. 3) write(iunit,'(a,2i5,2x,e12.5)') 
     :     'n,bond,probcum = ',n,m,probcum
         END DO ! m
       END DO ! n

       IF (ndebug .ge. 2) write(iunit,'(a,e12.5)') 'probcum = ',probcum
     
       IF (probcum .lt. 0.98*probtot .or. 
     :     probcum .gt. 1.1*probtot) THEN
         write(iunit,*) 'Somethings gone wrong! probcum out of range: ',
     :       probcum,probtot
         GOTO 201
       END IF

c 'roll the die' and get a bond
       n = numchan(iend)
       m = npbond(n,iend)
       ntry = 0
 99    CONTINUE
       tmp = ran0(iseed) 
       pick = probtot*tmp
       IF ( tmp .ge. 1.0 .or. pick .ge. bondprob(m,n,iend) ) THEN
         ntry = ntry + 1
         IF ( ntry .lt. 100 ) THEN
           GOTO 99
         ELSE
           tmp = 0.0  
           pick = 0.0
         ENDIF
       ENDIF  
       IF ( ntry .gt. 0 ) THEN
        write(iunit,*) 'Lightning debug: ntry,tmp,pick,probtot = ',
     :     ntry,tmp,pick,probtot
       ENDIF
       IF ( pick .gt. bondprob(m,n,iend) .and. probtot .gt. 0.0 .and.
     :   bondprob(m,n,iend) .gt. 0.0   ) THEN
        write(iunit,*) 
     :   'pick,probtot,probcum,bondprob(m,n,iend),m,n,iend = ',
     :   pick,probtot,probcum,bondprob(m,n,iend),m,n,iend
c        GOTO 99
       ELSEIF ( bondprob(m,n,iend) .eq. 0.0 ) THEN
       write(iunit,*) 'm,n,iend,bondprob(m,n,iend),probcum,probtot = ',
     :   m,n,iend,bondprob(m,n,iend),probcum,probtot
       write(iunit,*) 'numchan(iend) =', numchan(iend)
c       DO n=1,numchan(iend)
c        write(iunit,*) 'n,npbond = ', n, npbond(n,iend)
c         DO m=1,npbond(n,iend)
c         write(iunit,*) m,n,bondprob(m,n,iend)
c         END DO ! m
c       END DO ! n
       ENDIF
       
       npick = 0
       bpick = 0
       IF (ndebug .ge. 2) write(iunit,'(a,e12.5)') 'pick = ',pick
       DO n=1,numchan(iend)
         DO m=1,npbond(n,iend)
           IF (pick .lt. bondprob(m,n,iend)) THEN
             npick = n
             bpick = m
             GOTO 101
           END IF
         END DO ! m
       END DO ! n
  
 101  CONTINUE      
c
c  Now we have a new bond to add to the flash!
c
       IF (ndebug .ge. 3) write(iunit,'(a,2i5)') 
     :      'picked link,bond:',npick,bpick
       IF (npick .eq. 0 .or. bpick .eq. 0) THEN
         write(iunit,*) 'npick or bpick is zero .... stopping'
         write(iunit,*) 'pick, probtot = ',pick,probtot
         write(iunit,*) 'npick,bpick,numchan(1,2) = ',npick,bpick,
     :      numchan(1),numchan(2)
         GOTO 201
       END IF
       
       numchan(iend) = numchan(iend) + 1
       linksrc(numchan(iend),iend) = npick
       linkstep(numchan(iend),iend) = nstep
      i2 = chanijk(npick,1,iend)
      j2 = chanijk(npick,2,iend)
      k2 = chanijk(npick,3,iend)
      ibonds(numchan(iend),1,1,iend) = chanijk(npick,1,iend)
      ibonds(numchan(iend),2,1,iend) = chanijk(npick,2,iend)
      ibonds(numchan(iend),3,1,iend) = chanijk(npick,3,iend)
      
       i1 = bondijk(bpick,npick,1,iend)
       j1 = bondijk(bpick,npick,2,iend)
       k1 = bondijk(bpick,npick,3,iend)

      chandist(numchan(iend),iend) = gd(i1-i2,j1-j2,k1-k2,k2,1) +
     :                       chandist(npick,iend) 
c     :                       chandist(numchan(iend)-1,iend) 
       phichan(numchan(iend),iend) = phichan(npick,iend) -
     :            gd(i1-i2,j1-j2,k1-k2,k2,1)*signchan(iend)*eintz(k1)
       
       IF (ndebug .ge. 2) write(iunit,'(a,2(2x,e12.5))') 
     :      'phip(chan),phip(picked)',phichan(npick,iend),
     :       phip(i1,j1,k1)
       
       IF ( bcy .eq. 2 .and. ( j1 .eq. 0 .or. j1 .eq. nny ) ) THEN
         mask(i1,j1,k1) = 0.0
         phip(i1,j1,k1) = phichan(numchan(iend),iend)
         IF ( j1 .eq. 0   ) j1 = nny - 1
         IF ( j1 .eq. nny ) j1 = 1
         mask(i1,j1,k1) = 0.0
         phip(i1,j1,k1) = phichan(numchan(iend),iend)
       ENDIF

       IF ( bcx .eq. 2 .and. ( i1 .eq. 0 .or. i1 .eq. nnx ) ) THEN
         mask(i1,j1,k1) = 0.0
         phip(i1,j1,k1) = phichan(numchan(iend),iend)
         IF ( i1 .eq. 0   ) i1 = nnx - 1
         IF ( i1 .eq. nnx ) i1 = 1
       ENDIF
       chanijk(numchan(iend),1,iend) = i1
       chanijk(numchan(iend),2,iend) = j1
       chanijk(numchan(iend),3,iend) = k1
      ibonds(numchan(iend),1,2,iend) = i1
      ibonds(numchan(iend),2,2,iend) = j1
      ibonds(numchan(iend),3,2,iend) = k1
       
       
       IF (ndebug .ge. 2) write(iunit,'(a,3i5,a,3i5)') 
     :      'added bond from ', i2,j2,k2,
     :    ' to ',i1,j1,k1
       IF (ndebug .ge. 2) write(iunit,'(a,e13.5)') 
     :      'new point has phichan = ',phichan(numchan(iend),iend)
       
       IF (k1 .ge. nzs-3) THEN
         write(iunit,*) 'reached the top! stopping now..... ibal= ',ibal,
     :      ' icldtop,k1,zcldtop = ',icldtop,k1,zcldtop
         GO TO 201
c new test version 04.13.2004:
       ELSEIF ( ibal .ge. 2 ) THEN
        IF ( icldtop .gt. 0 .and. hgt(k1) .gt. zcldtop + 3000. 
c     :     .and.  numchan(iend) .gt. numchan(3-iend)
     :     ) THEN
         write(iunit,*) 'reached the top! stopping now..... ibal= ',ibal,
     :      ' icldtop,k1,zcldtop = ',icldtop,k1,zcldtop
         GOTO 201

c       ELSEIF ( ibal .eq. 2 ) THEN
c        IF ( icldtop .gt. 0 .and. k1 .gt. icldtop + 4 .and. 
c     :       numchan(iend) .gt. numchan(3-iend) ) THEN
c         write(iunit,*) 'reached the top! stopping now..... ibal=2 ',
c     :      ' icldtop,k1 = ',icldtop,k1
c         GOTO 201
        ENDIF
       END IF
       
       IF ( iextra .ge. 1 ) THEN
!         IF ( ibal .eq. 2 .and. iextra .ge. 5 ) THEN !  numchan(iend)/3 ) THEN
! 5/4/08 erm test following:
         IF ( ibal .eq. 2 .and. iextra .ge. 500 ) THEN !  numchan(iend)/3 ) THEN
           write(iunit,*) 'iextra = ',iextra
           write(iunit,*) 'Max amount of extra growth reached. Stop now.'
           GOTO 201
         ENDIF
         IF ( ibal .ne. 2 .and. iextra .gt. 50 ) THEN
           write(iunit,*) 'Max amount of extra growth reached. Stop now.'
           GOTO 201
         ENDIF
         iextra = iextra + 1
       ENDIF
c       IF ( 0.5*(z1d2lgt(k1,1)+z1d2lgt(Max(1,k1-1),1)) .lt. 0.5*cgthres
c     :     .and. istop(3-iend) .eq. 2 )  THEN
       IF ( 0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. 0.5*cgthres
     :     .and. istop(3-iend) .eq. 2 )  THEN
           write(iunit,*) 'Something screwy! Both ends at ground! ',
     :      'Stay with the first and stop now'
          GOTO 201
       ENDIF
       
c       IF ( ( 0.5*(z1d2lgt(k1,1)+z1d2lgt(Max(1,k1-1),1)) .lt. cgthres  .and. 
c     :      isa .ne. 1 .and. 
c     :        z1d2lgt(kzst,1) .gt. 2001.0 .and. istop(3-iend).le.1) .or.
c     :      ( k1 .le. 1  ) ) THEN
c        IF ( ibal .lt. 4 .or. (ibal .ge. 4 .and. 
c     :       Sign(1.0,phichan(numchan(iend),iend)) .eq. signchan(iend))
c     :      .or. ( k1 .le. 1  )  ) THEN

       IF ( ( 0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and. 
     :      isa .lt. 2 .and. 
     :        hgt(kzst) .gt. 2001.0 .and. istop(3-iend).le.1) .or.
     :      ( k1 .le. 1 .or. (k1 .eq. 2 .and. l2nde .ge. 2) ) ) THEN
        IF ( ibal .eq. 1 .or. ibal.eq.3 .or. 
     :      ( (ibal .ge. 4 .or. ibal .eq.2) .and. 
     :       Sign(1.0,phichan(numchan(iend),iend)-
     :            hgt(k1)*signchan(iend)*eintz(k1) ) 
     :        .eq. signchan(iend))
     :      .or. ( k1 .le. 1  )  ) THEN

         write(iunit,*) 'reached the ground! stopping now this channel',iend
         write(iunit,*) 'nstep = ',nstep
         write(iunit,'(a,f8.3)') 'fecritn = ',fecritn(iend)
         write(iunit,*) 'touchdown at i,j,k= ',i1,j1,k1
         write(iunit,'(a,2(1x,f8.3))') 
     :        'fecrit(1,2) = ', fecritn(1),fecritn(2)
         IF ( istop(3-iend) .eq. 2 ) THEN
           write(iunit,*) 'Something screwy! Both ends at ground! ',
     :      'Stay with the first and stop now'
          GOTO 201
         ENDIF
         IF ( iextra .eq. 0 ) iextra = 1
         icgyn = 1 ! Nint(signchan(iend))
         icgpn = Int( 1.1*signchan(iend) )
         istop(iend) = 2
         ihold = 0
         iendgnd = iend
         fecritn(3-iend) = Min( fecrit, fecritn(3-iend) )
c draw channel down to ground
          
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
           chandist(numchan(iend),iend) = gd(0,0,1,n,1) +
     :                       chandist(numchan(iend)-1,iend) 

           IF ( ibal .eq. 2 ) THEN
            phichan(numchan(iend),iend) = 
     :            phichan(numchan(iend)-1,iend) -
     :            gd(0,0,-1,n,1)*signchan(iend)*eintz(n)
           
           ELSE
           phichan(numchan(iend),iend) =
     :        phichan(numchan(iend)-(k1-n),iend)*Real(n)/Real(k1)
           ENDIF

           mask(i1,j1,n) = 0.0
           phip(i1,j1,n) = phichan(numchan(iend),iend)
           
           chanlen = chandist(numchan(iend),iend)
           
         END DO
         kmn = 1
         kmn1 = 1
c
c for ibal=2, reset the channel potential values to zero and relax
c
         IF ( ibal .eq. 2 ) THEN 
         
         IF ( cgfr .eq. 1.0 ) THEN ! set whole flash structure ground potential (zero)
         DO iend0 =1,nend
          DO n=1,numchan(iend0)
           i1 = ibonds(n,1,2,iend0)
           j1 = ibonds(n,2,2,iend0)
           k1 = ibonds(n,3,2,iend0)
c           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend0) = 0.0
            phip(i1,j1,k1) = phichan(n,iend0)
c           ELSE
c            phichan(n,iend0) = -eint*signchan(iend0)*
c     :                    (chanlen+chandist(numchan(iend0),iend0))
c            phip(i1,j1,k1) = phichan(n,iend0)
c           ENDIF
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
c           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend0) = phichan(n,iend0) - delta
            phip(i1,j1,k1) = phichan(n,iend0)
c           ELSE
c            phichan(n,iend0) = -eint*signchan(iend0)*
c     :                    (chanlen+chandist(numchan(iend0),iend0))
c            phip(i1,j1,k1) = phichan(n,iend0)
c           ENDIF
          END DO ! n
         END DO ! iend0
 
          IF ( .false. ) THEN
c first set ground channel values to a "missing" value
c         DO iend0 =1,nend
          DO n=1,numchan(iend)
c           IF ( istop(iend0) .eq. 2 ) THEN
            phichan(n,iend) = 1.e30
          END DO ! n
c         END DO ! iend0
c now go back along channel to initiation point
c         DO iend0 =1,nend
          n = numchan(iend)
          phichan(n,iend) = 0.0
          
          DO m=1,numchan(iend)
            i = linksrc(n,iend)
c            IF ( i .ge. 1 ) THEN
            x = chandist(n,iend) - chandist(i,iend)
            phichan(i,iend) = phichan(n,iend)
c            ELSEIF ( i .eq. 0 ) THEN
c              EXIT  ! then we have reached the init point
c            ENDIF
            IF ( i .eq. 1 ) EXIT
            n = i
          END DO ! m
c         END DO ! iend0
         
         ENDIF ! .false.
         
         ENDIF ! cgfr
         
         imn2 = Min(imn,Max(1,imn-2*iadd))
         jmn2 = Min(jmn,Max(1,jmn-2*iadd))
         imx2 = Max(imx,Min(nxs+nbw,imx+2*iadd))
         jmx2 = Max(jmx,Min(nys+nbs,jmx+2*iadd))
         kmx2 = nzs-1
         IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
         itype = 1
         IF (ndebug .ge. 1) write(iunit,*) 'strzsor 2'
         IF ( ny .gt. 2 ) THEN
         CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
     :      cofx,cofy,cofz, bcx, bcy,
     :      imn2,jmn2,kmn,imx2,jmx2,kmx2,itype,rjac,omega,ibc,soreps,0)
         ELSE
c         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
c     :      cofx,cofy,cofz, bcx, bcy,
c     :      imn2,jmn2,kmn,imx2,jmx2,kmx2,itype,rjac,omega,ibc,soreps,0)
         itype = 1
         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,nnx/2,1, nnz/2,nnx/2,
     :      cofx,cofy,cofz, bcx, bcy,
     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,soreps,0)
         ENDIF
c         GOTO 201
      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
         ENDIF ! ( ibal .eq. 2 )
       
       END IF
       ENDIF
       
       IF (j1 .le. nbs+1 .or. ( ny .gt. 2 .and. j1 .ge. nys+nbs+1 ) .or.
     :     i1 .le. nbw+1 .or. i1 .ge. nxs+nbw+1 ) THEN
         write(iunit,*) 'reached the side! stopping now..... '
         GO TO 201
       END IF

       mask(i1,j1,k1) = 0.0
       phip(i1,j1,k1) = phichan(numchan(iend),iend)
      
      i0 = Min(i1,i2)
      j0 = Min(j1,j2)
      k0 = Min(k1,k2)

      imn1 = Min(imn1,i1)
      jmn1 = Min(jmn1,j1)
      kmn1 = Min(kmn1,k1)
      imx1 = Max(imx1,i1)
      jmx1 = Max(jmx1,j1)
      kmx1 = Max(kmx1,k1)

c      iadd = 6
      iadd = NInt(4000./ds)
      iaddz = NInt(3000./dz)

      imn = Min(imn,Max(1,imn1-iadd))
      jmn = Min(jmn,Max(1,jmn1-iadd))
      kmn = Min(kmn,Max(1,kmn1-iaddz))
      imx = Max(imx,Min(nxs+nbw,imx1+iadd))
      jmx = Max(jmx,Min(nys+nbs,jmx1+iadd))
      kmx = Max(kmx,Min(nzs,kmx1+iaddz))
      
      IF ( ny .eq. 2 ) THEN ; jmn = 1 ; jmx = 1 ; ENDIF

       IF ( bcx .eq. 2 ) THEN
         IF ( imn .eq. 1   ) imx = nxs
         IF ( imx .eq. nxs ) imn = 1
       ENDIF

       IF ( bcy .eq. 2 ) THEN
         IF ( jmn .eq. 1   ) jmx = nys
         IF ( jmx .eq. nys ) jmn = 1
       ENDIF

      
      iloc0(iend) = i0
      jloc0(iend) = j0
      kloc0(iend) = k0
      
      

      END IF ! (istop(iend) .eq. 0) 
      END IF ! (istop(iend) .le. 1)
      END DO ! iend 

c  relax around each end in random order
      i1 = 1 + Int( 0.5 + ran0(iseed) )
      j1 = 3 - i1
      k1 = j1 - i1
      
      IF ( (nstep/intsormpi)*intsormpi == nstep ) THEN !{

      DO iend=i1,j1,k1
       IF (istop(iend) .eq. 0) THEN 
        IF (ndebug .ge. 2) write(iunit,*) 'call sor' 
        itype = 0     
        IF (ndebug .ge. 1) THEN
         write(iunit,*) 'strzsor 3, iend=',iend,'i,j,k =',
     :     iloc0(iend),jloc0(iend),kloc0(iend),rjac,omega
        ENDIF
        
        IF ( ny .gt. 2 ) THEN
        CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,
     :       iloc0(iend),jloc0(iend),kloc0(iend),mdel,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,ndebug)
        ELSE
        CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,
     :       iloc0(iend),jloc0(iend),kloc0(iend),mdel*lightintz,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,ndebug)
        ENDIF
      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
        IF (ndebug .ge. 2) write(iunit,*) 'done sor'
       END IF
      END DO


      IF ((nstep/intfg)*intfg .eq. nstep .or.
     :       nstep .eq. nstop )
     :   THEN
       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
      itype = 1
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 4'
      IF ( ny .gt. 2 ) THEN
      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ELSE
c        CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,
c     :       i0,j0,k0,mdel*lightintz,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,ndebug)
         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,nnx/2,1, nnz/2,nnx/2,
     :      cofx,cofy,cofz, bcx, bcy,
     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,soreps,0)
c      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ENDIF

      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
       
      END IF
      
      ENDIF !} intsormpi


c  every so often check the channel charges for imbalance.  If one end has
c  > 50% more charge than the other end, then we adjust the propagation
c  thresholds.

      IF ( ibal .eq. 1 .or. ibal .ge. 3) THEN
      IF ( (( nstep/3*3 .eq. nstep .or. nstep .le. 30 )  .and.
     :     ( istop(1) .le. 1 .and. istop(2) .le. 1 )) 
     :       .or. ihold .ne. 0 ) THEN
      
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
             IF ( (numchan(1) + numchan(2) ) .le. 4 .and. 
     :             numtry .lt. maxtry  .and. ninit .ge. 1 ) THEN
               GOTO 199
             ELSE
               IF ( ibal .eq. 1 .or. ibal .eq. 5 ) THEN
                 write(iunit,'(a)') 'Cant force growth. Stop now.'
                 GOTO 201
               ELSEIF ( ibal .ge. 3 ) THEN
c                write(iunit,*) 'nstep,istop2,ihold',nstep,istop(2),ihold
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
             IF ( (numchan(1) + numchan(2) ) .le. 4 .and. 
     :             numtry .lt. maxtry  .and. ninit .ge. 1 ) THEN
               GOTO 199
             ELSE
               IF ( ibal .eq. 1 .or. ibal .eq. 5 ) THEN
                 write(iunit,'(a)') 'Cant force growth. Stop now.'
                 GOTO 201
               ELSEIF ( ibal .ge. 3 ) THEN
c                 write(iunit,*) 'nstep,istop1,ihold',nstep,istop(1),ihold
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

          rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
     

        chg = rhotemp - rho(i1,j1,k1)
       
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)
         ENDIF

        END DO ! n

      END DO ! iend
       

c       IF ( ibal .eq. 5 ) THEN
c         IF ( Abs(qtot(1)) .gt. 1.1*Abs(qtot(2)) )  THEN
c              ihold = 1
c         ELSEIF ( Abs(qtot(2)) .gt. 1.1*Abs(qtot(1)) ) THEN
c              ihold = 2
c         ELSE
c             ihold = 0
c         END IF
c       ENDIF ! ibal .eq. 5
       
       IF ( ibal .eq. 1 .and. Abs(qtot(1)) .gt. 1.05*Abs(qtot(2))) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
            ihold = 1
       ELSEIF ( ibal .ge. 5 .and. 
     :            Abs(qtot(1)) .gt. 1.1*Abs(qtot(2))) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
            ihold = 1
       ELSEIF ( ibal .ge. 3 .and. ibal .lt. 5 .and.
     :          Abs(qtot(1)) .gt. 1.05*Abs(qtot(2))) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
            fecritn(2) = Max ( fecritn(2) - fecinc, fecmin )
            fecritn(1) = Min ( fecritn(1) + fecinc, fecrit)
            ihold = 0
       ELSEIF ( ibal .eq. 1 .and. 
     :          Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
            ihold = 2
       ELSEIF ( ibal .ge. 5 .and. 
     :          Abs(qtot(2)) .gt. 1.1*Abs(qtot(1)) ) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
            ihold = 2
       ELSEIF ( ibal .ge. 3 .and. ibal .lt. 5 .and.
     :           Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
         IF ( ndebug .ge. 2 ) THEN
          write(iunit,*) 'nstep,qtot1,qtot2= ',nstep,qtot(1),qtot(2)
         ENDIF
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

       
       IF ( nstep .ge. 40 .and. ibal .eq. 1 .and.
     :      ( Abs(qtot(1)) .gt. 2.0*Abs(qtot(2)) .or.
     :        Abs(qtot(2)) .gt. 2.0*Abs(qtot(1)) ) ) THEN
         write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         GOTO 201
       END IF
       IF ( nstep .ge. 40 .and. ibal .eq. 3 .and.
     :      ( Abs(qtot(1)) .gt. 3.0*Abs(qtot(2)) .or.
     :        Abs(qtot(2)) .gt. 3.0*Abs(qtot(1)) ) ) THEN
         write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         GOTO 201
       END IF
       IF ( nstep .ge. 40 .and. ibal .ge. 4 .and. 
     :   ( Abs(qtot(1)) .gt. 10. .or. Abs(qtot(2)) .gt. 10. ) .and.
     :   (idown .ne. 0 .and. fecritn(Max(1,idown)) .lt. fecrit ) .and.
     :      ( Abs(qtot(1)) .gt. 1.4*Abs(qtot(2)) .or.
     :        Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) ) ) THEN
         write(iunit,'(a)') 'Charge imbalance is too large. Stop now.'
         GOTO 201
       END IF
 
 200  CONTINUE
       
      END IF ! (nstep/5*5 .eq. nstep )
      
      ELSEIF ( ibal .eq. 2 ) THEN

c
c  Adjust channel reference potential to achieve approximate charge neutrality
c
      IF ( (( (nstep/3)*3 .eq. nstep .or. nstep .le. 40*lightintz .or. ny .eq. 2 ) 
     :      .or. ihold .ne. 0 .or. nohold .ge. 1 ) .and.
     :     ( istop(1) .le. 1 .and. istop(2) .le. 1 )) 
     :      THEN
!
! nohold:  >= 1 : adjust ref. potential after every step
!           <= 0 : adjust ref. potential after every 3rd step or when ihold /= 0
!    (if 0, imposes adjustment limit to be within the ambient values of the first two points)
!
c
      IF (  nohold .le. 0 .and.
     :     .not. ( ihold .eq. 1 .and. istop(2) .eq. 1 ) .and.
     :     .not. ( ihold .eq. 2 .and. istop(1) .eq. 1 )  ) THEN
c      
        DO iend =1,nend
          qtot(iend) = 0.0
        END DO ! iend
        DO iend =1,nend
         DO n=1,numchan(iend)
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)

!         IF ( ny > 2 ) THEN
          rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
!         ELSE
!          rhotemp = -eperao*(
!     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
!     +      cofy(j1,1)*phip(i1,j1,k1)+cofy(j1,2)*phip(i1,j1,k1) +
!     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
!     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
!         ENDIF


        chg = rhotemp - rho(i1,j1,k1)
       
       IF ( .not. ( rhotemp > -100. and. rhotemp < 100. ) ) THEN
         write(0,*) 'spark: problem with rhotemp = ',rhotemp
         write(0,*) phip(i1-1,j1,k1), phip(i1+1,j1,k1)
         write(0,*) phip(i1,j1-1,k1), phip(i1,j1+1,k1)
         write(0,*) phip(i1,j1,k1-1), phip(i1,j1,k1+1)
         write(0,*) phip(i1,j1,k1)  
         write(0,*) 'i1,j1,k1 = ',i1,j1,k1
         write(0,*) 'cofx,cofy,cofz = ',cofx(i1,1),cofy(j1,1),cofz(k1,1)
         write(0,*) 'iend,n = ',iend,n
         STOP
       ENDIF
       
       
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)
         ENDIF

         END DO ! n

       END DO ! iend
       
       IF ( Abs(qtot(1)) .gt. 1.05*Abs(qtot(2)) )  THEN
            ihold = 1
       ELSEIF ( Abs(qtot(2)) .gt. 1.05*Abs(qtot(1)) ) THEN
            ihold = 2
       ELSE
           ihold = 0
       ENDIF

c       write(iunit,*) 'balck01: nstep = ',nstep,ihold,
c     :   qtot(1)*dls3,qtot(2)*dls3

            IF ( ( ((Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) .and. Abs(qtot(1))*dls3*dsfac .gt. 2.0) .or.
     :              (Abs(qtot(1)) .gt. 1.4*Abs(qtot(2)) .and. Abs(qtot(2))*dls3*dsfac .gt. 2.0 )) .and.
     :           nstep .ge. 40*lightintx) .or.
     :   (Abs(qtot(2))*dls3*dsfac .gt. 60.0 .and. Abs(qtot(1))*dls3*dsfac .gt. 60.0
     :          .and. (numchan(1) .gt. 2*numchan(2) .or.
     :          numchan(2) .gt. 2*numchan(1) ) ) ) THEN
             write(iunit,*) 'Flash is too unbalanced, stopping 1'
          write(iunit,*) 'q1,q2= ',qtot(1)*dls3*dsfac,qtot(2)*dls3*dsfac
              GOTO 201
            ENDIF
            

       
c       
       ELSE ! (one end needs to grow but cant, so adjust potential to balance charge)
            ! OR nohold >= 1 and always adjust reference to balance
c       
!       maxiter = Max(1,nohold)
!       DO iter = 1,maxiter
      k = 0  ! counts how many cycles of adjusting phiref and strzsor with qtot1 > qtot2
      j = 0  ! counts how many cycles of adjusting phiref and strzsor with qtot2 > qtot1
      ilast = 0
      iflip = 0
      istop1 = 0
      bfac = 0.02
      IF ( ny .le. 2 ) bfac = 0.075
      qtot1 = 0.0
      qtot2 = 0.0
         
 21   CONTINUE
 
      IF ( iflip .ge. 5 ) THEN
        write(iunit,*) 'Too many flips'
        GOTO 31
      ENDIF
       
      DO iend =1,nend
      qtot(iend) = 0.0
      ENDDO

      DO iend =1,nend
        DO n=1,numchan(iend)
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)

          rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
          


        chg = rhotemp - rho(i1,j1,k1)
       

         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k1)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k1)
         ENDIF

        END DO ! n
      END DO ! iend

c       write(iunit,*) 'BALCK02: nstep = ',nstep,ihold,
c     :  qtot(1)*dls3,qtot(2)*dls3
       

            IF ( ( ((Abs(qtot(2)) .gt. 1.4*Abs(qtot(1)) .and. Abs(qtot(1))*dls3*dsfac .gt. 2.0) .or.
     :              (Abs(qtot(1)) .gt. 1.4*Abs(qtot(2)) .and. Abs(qtot(2))*dls3*dsfac .gt. 2.0 )) .and.
     :           nstep .ge. 40*lightintx) .or.
     :   (Abs(qtot(2))*dls3*dsfac .gt. 60.0 .and. Abs(qtot(1))*dls3*dsfac .gt. 60.0
     :          .and. (numchan(1) .gt. 2*numchan(2) .or.
     :          numchan(2) .gt. 2*numchan(1) ) ) ) THEN
             write(iunit,*) 'Flash is too unbalanced, stopping 2'
          write(iunit,*) 'q1,q2= ',qtot(1)*dls3*dsfac,qtot(2)*dls3*dsfac
              GOTO 201
            ENDIF

c       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
       IF ( istop1 .eq. 1 ) THEN 
         write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
         write(iunit,*) 'phistart = ',phistart
         write(iunit,*) 'number of iterations (j,k) = ',j, k
         GOTO 201
       ELSEIF ( istop1 .eq. -1 ) THEN
c         write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
c         write(iunit,*) 'phistart = ',phistart
c         write(iunit,*) 'number of iterations (j,k) = ',j, k
          ihold = 0
          GOTO 31
       ENDIF
 
            IF ( ny .eq. 2 .and. ( k .gt. 5 .or. j .gt. 5 ) ) THEN
             bfac = 1.1*bfac
            ENDIF

            IF ( k .gt. 45 .or. j .gt. 45 ) THEN
             write(iunit,*) 'Too many iterations (j,k) = ',j,k
             write(iunit,*) 'qtot(1,2) = ',qtot(1)*dls3*dsfac,qtot(2)*dls3*dsfac
             write(iunit,*) 'phistart,bfac,iflip = ',phistart,bfac,iflip
             write(iunit,*)  'phisthi,phistlow = ',phisthi,phistlow
             bfac = 2.0*bfac
            ENDIF

            IF ( k .gt. 50 .or. j .gt. 50 ) THEN
             write(iunit,*) 'Too many iterations (j,k) = ',j,k
             GOTO 201
c             STOP
            ENDIF
cc
       IF ( ( Abs(qtot(1)) .gt. 1.02*Abs(qtot(2)) .and. 
     :              j.eq.0 .and. k.eq.0) .or.
     :        Abs(qtot(1)) .gt. 1.01*Abs(qtot(2) ) ) THEN ! {
cc
         IF ( qtot1 .eq. 0.0 )  qtot1 = qtot(1)
         
         IF ( qtot(1) .lt. qtot1 ) THEN  ! {
cc       write(iunit,*) 'qtot(1,2) = ',qtot1,qtot2
            qtot1 = qtot(1)
cc       write(iunit,*) 'nstep = ',nstep
c            write(iunit,*) 'Going the wrong way! nstep =',nstep
cc       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
c       write(iunit,*) 'number of iterations (j,k) = ',j, k
c            istop1 = 1
            istop1 = -1
c            GOTO 21
            
            bfac = -Abs(bfac)
         ELSE
            qtot1 = qtot(1)
         ENDIF  ! }

         IF ( ilast .eq. 2 ) THEN
           iflip = iflip + 1
           bfac = 0.5*bfac
         ENDIF
         ilast = 1
         
c phi(i1,j1,k1),phi(i2,j2,k2),phi(i3,j3,k3)
         IF ( j .gt. 0 .and. k .eq. 0 .and. iflip .gt. 0 ) THEN
c           bfac = 0.5*bfac
         ELSEIF ( iflip .gt. 0 ) THEN
c           bfac = 0.5*bfac
         ENDIF

         
         IF ( istop1 .le. 0 ) THEN
         phistart = phistart + bfac*(phidiff)
         adjlast = bfac*phidiff
         ELSE
         phistart = phistart - adjlast
         write(iunit,*) 'phistart = ',phistart
         ENDIF
         
c 5.6.2004 ERM [
         IF (  nohold .eq. 0 .and. 
     :         ( phistart .gt. phisthi .or. phistart .lt. phistlow) ) THEN
c           phistart =  Min( phisthi, Max( phistlow, phistart ) )
           write(iunit,*) 'Reached limit of phistart adjustment:',
     :               ' phistartn =',phistart
           GOTO 201
         ENDIF
c 5.6.2004 ERM ]
         
         DO iend =1,nend
          DO n=1,numchan(iend)
           i1 = ibonds(n,1,2,iend)
           j1 = ibonds(n,2,2,iend)
           k1 = ibonds(n,3,2,iend)
            x = phichan(n,iend)
            phichan(n,iend) = phistart - 
     :         ( chandist(n,iend) - chandist(1,iend))*
     :           signchan(iend)*eintz(k1)
            phip(i1,j1,k1) = phichan(n,iend)
c            write(iunit,*) 'nst,iend,n,old,new ',
c     :       nstep,iend,n,x,phichan(n,iend)

          END DO ! n
         END DO ! iend
            k = k + 1

       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
      itype = 1
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 5'
      IF ( ny .gt. 2 ) THEN
      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ELSE
c        CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,
c     :       i0,j0,k0,mdel*lightintz,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,ndebug)
         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,nnx/2,1, nnz/2,nnx/2,
     :      cofx,cofy,cofz, bcx, bcy,
     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,soreps,0)
c      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ENDIF
      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF
      

            GOTO 21
cc
       ELSEIF ( ( Abs(qtot(2)) .gt. 1.02*Abs(qtot(1)) .and. 
     :              j.eq.0 .and. k.eq.0) .or.
     :        Abs(qtot(2)) .gt. 1.01*Abs(qtot(1) ) ) THEN
cc
         
         IF ( qtot2 .eq. 0.0 )  qtot2 = qtot(2)
         
         IF ( qtot(2) .gt. qtot2 ) THEN
cc       write(iunit,*) 'qtot(1,2) = ',qtot1,qtot2
            qtot2 = qtot(2)
cc       write(iunit,*) 'nstep = ',nstep
c            write(iunit,*) 'Going the wrong way! nstep =',nstep
cc       write(iunit,*) 'qtot(1,2) = ',qtot(1),qtot(2)
c       write(iunit,*) 'number of iterations (j,k) = ',j, k
c            istop1 = 1
            istop1 = -1
c            GOTO 201
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
c           bfac = 0.5*bfac
         ELSEIF ( iflip .gt. 0 ) THEN
c           bfac = 0.5*bfac
         ENDIF


         IF ( istop1 .le. 0 ) THEN
         phistart = phistart - bfac*(phidiff)
         adjlast = -bfac*phidiff
         ELSE
         phistart = phistart - adjlast
         write(iunit,*) 'phistart = ',phistart
         ENDIF
         
         DO iend =1,nend
          DO n=1,numchan(iend)
           i1 = ibonds(n,1,2,iend)
           j1 = ibonds(n,2,2,iend)
           k1 = ibonds(n,3,2,iend)
            x = phichan(n,iend)
            phichan(n,iend) = phistart - 
     :         ( chandist(n,iend) - chandist(1,iend))*
     :           signchan(iend)*eintz(k1)
            phip(i1,j1,k1) = phichan(n,iend)
c            write(iunit,*) 'iend,n,old,new ',iend,n,x,phichan(n,iend)
          END DO ! n
         END DO ! iend
            j = j + 1

       IF (ndebug .ge. 2) write(iunit,*) 'full flash sor'      
      itype = 1
      IF (ndebug .ge. 1) write(iunit,*) 'strzsor 6'
      IF ( ny .gt. 2 ) THEN
      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
     :       cofx,cofy,cofz, bcx, bcy,
     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ELSE
c        CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,
c     :       i0,j0,k0,mdel*lightintz,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,ndebug)

         CALL strzsor2d(mask,rhs0,phip,nnx,nny,nnzt,nzs,nnx/2,1, nnz/2,nnx/2,
     :      cofx,cofy,cofz, bcx, bcy,
     :      1,1,1,nnx-1,1,nnzt-1,itype,rjac,omega,ibc,soreps,0)
c      CALL strzsor(mask,rhs0,phip,nnx,nny,nnzt,nzs,i0,j0,k0,mdelf,
c     :       cofx,cofy,cofz, bcx, bcy,
c     :       imn,jmn,kmn,imx,jmx,kmx,itype,rjac,omega,ibc,soreps,0)
      ENDIF

      IF ( bcy .eq. 2 ) THEN
       DO n = 1,2
         phip(-1:nnx+2,1-n    ,0:nnzt+1) = phip(-1:nnx+2,nny-n,0:nnzt+1)
         phip(-1:nnx+2,nny-1+n,0:nnzt+1) = phip(-1:nnx+2,n    ,0:nnzt+1)
       ENDDO
      ENDIF
      IF ( ny .eq. 2 ) THEN
         phip(-1:nnx+2,0,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
         phip(-1:nnx+2,2,0:nnzt+1) = phip(-1:nnx+2,1,0:nnzt+1)
      ENDIF

            GOTO 21
cc
       ELSE
cc
c         write(iunit,*) 'End of balance qtot(1,2) = ',qtot(1),qtot(2)
c         write(iunit,*) 'phistart = ',phistart
c         write(iunit,*) 'nstep, number of iterations (j,k) = ',nstep,j, k
         ihold = 0
cc
       END IF
cc       
 31   CONTINUE
c       
       ENDIF  ! ( .not. ( ihold ... ) )
c       

      END IF ! (nstep/5*5 .eq. nstep )
      
      ENDIF ! (ibal)
      
      
      IF (ndebug .ge. 2) write(iunit,'(a,i4)') 'nstep = ',nstep
      
      END DO ! nstep 
      
      write(iunit,'(a)') 'Did not reach boundary'           
      
 201  CONTINUE
 
       write(iunit,'(a,i4)') 'Number of lightning steps: nstep = ',nstep

      tim2 = commas_dtime(timarr)
!      IF ( lsgi ) tim2 = commas_dtime(timarr)
!      IF ( .not. lsgi ) tim2 = second()
      write(iunit,*) 'time for propagation = ',tim2,timarr(1)
c      IF ( .not. lsgi ) write(iunit,*) 'time for propagation = ',tim2-tim1
c      write(0,'(a,2(1x,1pe12.5))') ' time for prop = ',tim2,timarr(1)
      write(iunit,'(a,a,3(i4,1x),4x,3(i4,1x))') 'Flash bounding box ',
     :    'imn,jmn,kmn,imx,jmx,kmx: ',imn1,jmn1,kmn1,imx1,jmx1,kmx1
      write(iunit,'(a)') 'End of lightning propagation'
      write(iunit,'(a,2(1x,i4))') 'numchan(1),numchan(2)',
     :         numchan(1),numchan(2)

c
c  Check energy
c
      eold = 0.0
      enew = 0.0
      etst = 0.0
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

c      jmn2 = 2
c      jmx2 = 1

      ist = Max(1,nbw+2)
      jst = Max(1,nbs+2)
      jend = Max(1, nny-nbn-1)
      DO k1 = 1, nzs-2
        vol = dls3*dsfac*zmap(k1)
      DO j1 = jst,jend
      DO i1 = ist,nnx-nbe-1

c      DO k1=kmn,kmx
c        vol = dls3*dsfac*zmap(k1)
c       DO j1=jmn,jmx
c        DO i1=imn,imx
         x = 0.5*phi2(i1,j1,k1)*rho(i1,j1,k1)*vol
         eold = eold + x

         z = 0.5*phip(i1,j1,k1)*rho(i1,j1,k1)*vol
         etst = etst + z
        
        IF ( mask(i1,j1,k1) .eq. 1 ) THEN
          y = 0.5*phip(i1,j1,k1)*rho(i1,j1,k1)*vol
c          enew = enew + 0.5*phip(i1,j1,k1)*rho(i1,j1,k1)*vol
        ELSE
           rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

          y = 0.5*phip(i1,j1,k1)*rhotemp*vol
c          enew = enew + 0.5*phip(i1,j1,k1)*rhotemp*vol
        ENDIF
        
          enew = enew + y

         IF ( ibal.eq.2 .and. icldtop .ne. 0 .and. 
     :        k1 .gt. icldtop + 4 ) THEN
           eoldct = eoldct + x
           enewct = enewct + y
         ENDIF
        
        
        ENDDO
       ENDDO
      ENDDO
      
      IF ( ibal .eq. 2 ) THEN
       write(iunit,'(a,3(1x,1pe13.5))') 'phiref, phistart, delta-phi = ',
     :   phistart,phistart0,phistart-phistart0
      
      ENDIF
       write(iunit,'(a,3(1x,1pe13.5))') 
     :    'Energy check (old,new,diff): ',eold,enew,enew-eold

       write(iunit,'(a,1(1x,1pe13.5))') 
     :    'Energy test : ',etst

       lct = .true.

       IF ( enewct .ne. 0.0 .or. eoldct .ne. 0.0 ) THEN
         write(iunit,'(a,3(1x,1pe13.5))') 
     :  'CT Energy check (old,new,diff): ',eoldct,enewct,enewct-eoldct
        lct = .true.
c        IF ( (eoldct .gt. 0.0 .and. enewct .gt. 0.0 .and.
c     :        enewct .gt. 1.05*eoldct)  .or.
c     :    ( Sign(1.0,enewct) .ne. Sign(1.0,eoldct) ) .or.
c     :    ( Abs(enewct) .gt. 1.1*Abs(eoldct) .and.
c     :      Sign(1.0,enewct) .eq. Sign(1.0,eoldct)  ) ) THEN
c          lct = .false.
c        ENDIF
        
        
       ENDIF
       IF ( .not. (enew .gt. 0.0 .and. enew .lt. eold .and. lct )
     :      .and. ibal .eq. 2 
     :      .and. bcy .ne. 2 .and. .false. 
     :      .and. .not. ( cgfr .lt. 1.0 .and. icgyn .ne. 0 ) ) THEN
c
c 2005.10.03: let it pass if cgfr .lt. 1 and is a CG.
c
         ldumpsp = .true.
         write(iunit,*) 'Failed energy reduction! Try again'

      DO iend =1,nend
        qtot(iend) = 0.0
        rhotot(iend) = 0.0
        rhotot1(iend) = 0.0
      END DO ! iend
      
      DO iend =1,nend
        DO n=1,numchan(iend)
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)
         i2 = ibonds(n,1,2,iend)-nbw-1
         j2 = ibonds(n,2,2,iend)-nbs-1
         k2 = ibonds(n,3,2,iend)

         rhotemp = 0.0
         chg = 0.0
         IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN

          rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )
        
         chg = rhotemp - rho(i1,j1,k1)
         
         ELSE
          chg = 0.0
          rhotemp = 0.0
         ENDIF

        chg0 = rho(i1,j1,k1)
        
         IF ( (ibal .ge. 5 .or. ibal .eq. 2) .and. 
     :        0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and.
     :        Sign(1.0,chg) .ne. signchan(iend) ) THEN
           chg = 0.0
           rhotemp = rho(i1,j1,k1)
         ENDIF
          
         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k2)
           rhotot(1) = rhotot(1) + chg0*zmap(k2)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k2)
           rhotot(2) = rhotot(2) + chg0*zmap(k2)
         ENDIF
        END DO ! n
      END DO ! iend
      
      
      DO iend=1,nend
      write(iunit,'(a,i1,a,1pe13.5)') 'initial charge on channel qtot(',
     :             iend, ')= ', qtot(iend)*dls3*dsfac
      write(iunit,'(a,1pe13.5)') 'rhotot  = ',rhotot(iend)*dls3*dsfac
      write(iunit,'(a,1pe13.5)') 'rhotot1 = ',rhotot1(iend)*dls3*dsfac
       chgdis0(iend) = 0.0
      END DO ! iend

        IF ( numtry .eq. 1 ) THEN
          write(iunit,'(a)') 'Try 1'
        ENDIF
         
         IF ( isa .eq. 1 ) THEN
           STOP
         ENDIF

         lgtstp = 2
         GOTO 199
       ELSE
         lgtstp = 0
       ENDIF  ! energy check
c
c write text output file
c      
      times = (nstep1-1)*(dtp)
      
      IF ( isa .lt. 1 ) THEN
      write(nmliter,915) iliter
 915  format(i2.2)
      write(rstime,'(i6.6)') time
c 912  format(i5.5)
!      open(unit=16,file=hdrcfs(1:lhdrcfs)//'flashes/'//
!     :  filehead(1:lfilehead)
!     :  //'.chnl.'//rstime(1:ihttima)//'.'//nmliter, form='formatted',
!     :  status = 'unknown')
!       rewind(16)
       
       write(16,*)
       write(16,'(a,i7,a,i4)') 'FLASH at time = ', time, ', iliter = ',iliter

      ENDIF ! ( isa .lt. 1 )
      
      IF ( isa .ge. 1 ) THEN
        write(16,'(a,i5)') 'Flash number ',iliter
      ENDIF ! ( isa .eq. 1 )


      write(16,'(5(1x,i4),3(1x,f8.1),2(1x,i4))') 
     :   numchan(1),numchan(2),nxb,nyb,nzs,dx,dy,dz,nbw,nbs
      pchnl = numchan(1)
      nchnl = numchan(2)
!      volbi = 1.0/(dx1*dy1*dz1)
      DO iend =1,nend
      qtot(iend) = 0.0
      rhotot(iend) = 0.0
      rhotot1(iend) = 0.0
      chgdis0(iend) = 0.0
      chgdis(iend) = 0.0
c      icnt1(iend) = 0
c      icnt2(iend) = 0
c      icnt3(iend) = 0
c      icnt4(iend) = 0
      ENDDO
      DO iend =1,nend
        DO n=1,numchan(iend)
         i1 = ibonds(n,1,2,iend)
         j1 = ibonds(n,2,2,iend)
         k1 = ibonds(n,3,2,iend)
         i2 = ibonds(n,1,2,iend)-nbw-1
         j2 = ibonds(n,2,2,iend)-nbs-1
         k2 = ibonds(n,3,2,iend)

         IF ( n > 1 ) THEN
          chanlen =  chandist(n,iend) - chandist( linksrc(n,iend), iend)
         ELSE
          chanlen =  chandist(n,iend)
         ENDIF

         rhotemp = 0.0
         chg = 0.0
         IF ( k1 .gt. 1 .and. k1 .lt. nzs-2 ) THEN

          rhotemp = -eperao*(
     +      cofx(i1,1)*phip(i1-1,j1,k1)+cofx(i1,2)*phip(i1+1,j1,k1) +
     +      cofy(j1,1)*phip(i1,j1-1,k1)+cofy(j1,2)*phip(i1,j1+1,k1) +
     +      cofz(k1,1)*phip(i1,j1,k1-1)+cofz(k1,2)*phip(i1,j1,k1+1) +
     +      (cofx(i1,3)+cofy(j1,3)+cofz(k1,3))*phip(i1,j1,k1)  )

        
         chg = rhotemp - rho(i1,j1,k1)
         
         ELSE
          chg = 0.0
          rhotemp = 0.0
         ENDIF


        chg0 = rho(i1,j1,k1)
        
         IF ( (ibal .ge. 5 .or. ibal .eq. 2) .and. 
     :        0.5*(hgt(k1)+hgt(Max(1,k1-1))) .lt. cgthres  .and.
     :        Sign(1.0,chg) .ne. signchan(iend) ) THEN
           chg = 0.0
           rhotemp = rho(i1,j1,k1)
         ENDIF
          
         rhotot1(iend) = rhotot1(iend) + rhow(i2,j2,k2)*zmap(k2)
         IF ( chg  .lt. 0.0 ) THEN
           qtot(1) = qtot(1) + chg*zmap(k2)
           rhotot(1) = rhotot(1) + chg0*zmap(k2)
         ELSE 
           qtot(2) = qtot(2) + chg*zmap(k2)
           rhotot(2) = rhotot(2) + chg0*zmap(k2)
         ENDIF

         write(16,'(3(1x,i4),3(2x,i6))') (ibonds(n,i,1,iend), i=1,3),
     :         n, linksrc(n,iend),linkstep(n,iend)
         IF ( ibal .lt. 5 .and. ibal .ne. 2 ) THEN
         write(16,'(3(1x,i4),5(1x,1pe13.5))') 
     :      (ibonds(n,i,2,iend), i=1,3),phi(i1,j1,k1),
     :      rho(i1,j1,k1),rhotemp,
     :       chandist(n,iend),chg
         ELSE
         write(16,'(3(1x,i4),6(1x,1pe13.5))') 
     :      (ibonds(n,i,2,iend), i=1,3),phi(i1,j1,k1),
     :      rho(i1,j1,k1),rhotemp,
     :       chandist(n,iend),chg,phichan(n,iend)
         
         ENDIF

         i1 = i1-nbw-1
         j1 = j1-nbs-1
c          x = 0.5*ds*Real(i1 + i2 - 2)
c          y = 0.5*ds*Real(j1 + j2 - 2)
          x = ds*Real(i1 - 1)
          y = ds*Real(j1 - 1)
          z = hgt(k1)
c          i = Int( (x + eps1)/dx1 + 0.5) + 1
c          j = Int( (y + eps1)/dy1 + 0.5) + 1
c          k = Int( (z + eps1)/dz1 + 0.5) + 1
          
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
          

      IF ( .not. ( icgyn .eq. 1 .and. iend .eq. iendgnd ) 
     :          .or. ibal .ge. 5  .or. ibal .eq. 2 
     :   ) THEN
      
        IF ( (dx .eq. dx1 .and. dy .eq. dy1 .and. dz .eq. dz1) .or.
     :        igrid .eq. 1 ) THEN
          volbi = gxt(i,3)*gyt(j,3)*gzt(k,3)
!         IF ( istretch ) print*, 'volbi,ds3,chg ', volbi,ds3,chg,(chg)*volbi*ds3
         t3(i,j,k) = t3(i,j,k) + (chg)*volbi*ds3*dsfac
        ELSE
c         IF ( .false. ) THEN
c
c Figure out which grid boxes should recieve charge and what fraction
c
c  First find which coarse point is closest to each corner of the fine grid
c  xw,xe,ys,yn define the edges of the fine grid square
         xw = x - 0.5*ds
         xe = x + 0.5*ds
         ys = y - 0.5*ds
         yn = y + 0.5*ds
c note: zu and zd are only different when using a constant vertical spacing
c and dslight < dz.
         zu = z + 0.5*ds
         zd = z - 0.5*ds
c
c iw,ie,js,jn are the points on the coarse grid closest to the fine grid edges
c         
         iw = NInt( xw/dx1 + 1.0 )
         ie = NInt( xe/dx1 + 1.0 )
         js = NInt( ys/dy1 + 1.0 )
         jn = NInt( yn/dy1 + 1.0 )
         ku = Int( zu/dz1 + 1.0 )
         kd = Int( zd/dz1 + 1.0 )
         IF ( istretch .or. dz .eq. dz1 ) THEN
          ku = k1
          kd = k1
         ENDIF
          k = kd
         
c
c  Now find out how many coarse grid cells are overlapped (1,2, or 4)
c    (xc,yc,zc) is the intersection point of the coarse grid cells within
c    the fine grid box and is used for determining the area/volume overlaps
c
c          chgdis0(iend) = chgdis0(iend) + (chg)*volbi*ds3*dxl*dyl*dzl/
c     :     gt(1,1,k1,imapz)
      
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
          yc = gyt(js,1) ! ( Real(js) - 0.5 )*dy1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          IF ( frs.lt.0.0 .or. frs.gt.1.0) THEN
           write(iunit,*) 'nuts! frs,frn = ',frs,frn
          ENDIF
          t3(ie,js,k) = t3(ie,js,k) + frs*(chg)*volbi*ds3
          t3(ie,jn,k) = t3(ie,jn,k) + frn*(chg)*volbi*ds3
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 0 ) THEN ! two cells (e-w)
          xc = gxt(iw,1) ! ( Real(iw) - 0.5 )*dx1
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          IF ( frw.lt.0.0 .or. frw.gt.1.0) THEN
           write(iunit,*) 'nuts! frw,fre = ',frw,fre
          ENDIF
          t3(iw,js,k) = t3(iw,js,k) + frw*(chg)*volbi*ds3
          t3(ie,js,k) = t3(ie,js,k) + fre*(chg)*volbi*ds3
          
        ELSEIF ( ie - iw .eq. 1 .and. jn - js .eq. 1 ) THEN  ! four cells
          yc = gyt(js,1) ! ( Real(js) - 0.5 )*dy1
          xc = gxt(iw,1) ! ( Real(iw) - 0.5 )*dx1
          frs = (yc - ys)/dy
          frn = 1.0 - frs
          frw = (xc - xw)/dx
          fre = 1.0 - frw
          IF ( frs.lt.0.0 .or. frs.gt.1.0) THEN
           write(iunit,*) 'nuts4! frs,frn,dx,dy = ',frs,frn,dx,dy
           write(iunit,*) 'nuts4: yc,ys,yc-ys: ',yc,ys,yc-ys
           frs = Max(0.0, Min( 1.0, frs ) )
           frn = Max(0.0, Min( 1.0, frn ) )
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
          write(iunit,*) 'Something wrong! iw,ie,js,jn,kd,ku = ',
     :               iw,ie,js,jn,kd,ku
          STOP
        ENDIF

c         t4(i,j,k)  = -Sign(1.0,t3(i,j,k) )
      ENDIF
      ENDIF
         t2(i,j,k) = t2(i,j,k) + signchan(iend)
         t1(i,j,k) = t1(i,j,k) + chanlen

         IF ( n .eq. 1 .and. iend .eq. 1 .and. 
     :             .not. ( isa .ge. 1 .and. irand .eq. 1 )  ) THEN 
            ixst = i 
            jyst = j 
            kzst = k
          END IF

        END DO ! n
      END DO ! iend
      
      
      DO iend=1,nend
      write(iunit,'(a,i1,a,1pe13.5)') 'initial charge on channel qtot(',
     :             iend, ')= ', qtot(iend)*dls3*dsfac
      write(iunit,'(a,1pe13.5)') 'rhotot  = ',rhotot(iend)*dls3*dsfac
      write(iunit,'(a,1pe13.5)') 'rhotot1 = ',rhotot1(iend)*dls3*dsfac
      IF ( rhotot(iend) .ne. 0.0 ) THEN
        write(iunit,'(a,1pe13.5)') 'fraction of charge removed = ',
     :       qtot(iend)/rhotot(iend)
      ENDIF
      END DO ! iend

      write(iunit,*) 'count chgdis0(1,2) = ',chgdis0(1),chgdis0(2), 
     :      chgdis0(1)+chgdis0(2)
      DO iend=1,nend
       chgdis0(iend) = 0.0
      END DO ! iend

      DO iend=1,nend
          IF ( ( ( Abs(qtot(iend))*dls3*dsfac .lt. 0.1 .and. ny .gt. 2 ) .or.
     :      (Abs(qtot(iend))*dls3*dsfac .lt. 0.04 .and. ny .eq. 2) ) .and. 
     :         .not. ( icgyn .ge. 1 .and. ibal .eq. 2 ) ) THEN
            write(iunit,*) 'Not enough charge on leader',iend,
     :          '... go back and try again.'
       
         t3(:,:,:) = 0.0
         t2(:,:,:) = 0.0

c 06.30.2004 [


c 06.30.2004 ]
            GOTO 199
          ENDIF
      END DO ! iend


       IF ( ibal .eq. 1 .or. ibal .ge. 3 ) THEN
         write(iunit,'(a,2(1x,f8.3))') 
     :        'fecrit(1,2) = ', fecritn(1),fecritn(2)
       ENDIF

      IF ( isa .lt. 1 ) THEN
!        close(16)
      ELSE
        IF ( kmn1 .le. 7 ) THEN
          write(iunit,'(a,i4)') 'Potential CG: kmn1 = ',kmn1
        ENDIF
        nkzmn(kmn1) = nkzmn(kmn1) + 1
        GOTO 45
      ENDIF

c now get net charge in t3 array before smoothing:
      chgdis0(1) = 0.0
      chgdis0(2) = 0.0
      DO kz=1,nz
      DO jy=1,ny
      DO ix=1,nx
       IF ( .not. ( t3(ix,jy,kz) .lt. 1000. .and. t3(ix,jy,kz) .gt. -1000.) ) THEN
         write(0,*) 'PROBLEM with t3! (presmooth) ',ix,jy,kz,t3(ix,jy,kz) 
         STOP
       ENDIF
       IF ( t3(ix,jy,kz) .gt. 0.0 ) THEN
         chgdis0(2) = chgdis0(2) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
       ELSEIF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
         chgdis0(1) = chgdis0(1) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
       ENDIF
      END DO
      END DO
      END DO
      write(iunit,*) 'presmooth chgdis0(1,2) = ',chgdis0(1),chgdis0(2), 
     :      chgdis0(1)+chgdis0(2)

c very light smoothing of channel charge...
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
c      i=0
      DO j=1,i
       CALL smthlt1 (nx,ny,nz,nor,t3,t4,fx,fz) 
      ENDDO     
      
c now set t4 array with flags of where lightning charge is dumped
      chgdis0(1) = 0.0
      chgdis0(2) = 0.0
      DO kz=1,nz
      DO jy=1,ny
      DO ix=1,nx
        IF ( Abs( t3(ix,jy,kz) ) .gt. 1.0e-18 ) THEN
          t4(ix,jy,kz) = Sign(1.0,t3(ix,jy,kz) )
        ELSE
          t4(ix,jy,kz) = 0.0
        END IF
       IF ( .not. ( t3(ix,jy,kz) .lt. 1000. .and. t3(ix,jy,kz) .gt. -1000.) ) THEN
         write(0,*) 'PROBLEM with t3! (post smooth) ',ix,jy,kz,t3(ix,jy,kz) 
         STOP
       ELSEIF ( t3(ix,jy,kz) .gt. 0.0 ) THEN
         chgdis0(2) = chgdis0(2) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
!         print*,' t3 = ',ix,jy,kz,t3(ix,jy,kz)
       ELSEIF ( t3(ix,jy,kz) .lt. 0.0 ) THEN
         chgdis0(1) = chgdis0(1) + t3(ix,jy,kz)*dxx(ix)*dyy(jy)*dzz(kz) ! *dxl*dyl*dzl/gt(1,1,kz,imapz)
!         print*,' t3 = ',ix,jy,kz,t3(ix,jy,kz)
       ENDIF
      END DO
      END DO
      END DO
      
      write(iunit,*) 'recount chgdis0(1,2) = ',chgdis0(1),chgdis0(2), 
     :      chgdis0(1)+chgdis0(2)
!      STOP

 999  CONTINUE 

      deallocate ( mask )
      deallocate ( phip )
   
      deallocate ( hgt )
      deallocate ( hgtdg )
      deallocate ( zmap )
      deallocate ( cofx )
      deallocate ( cofy )
      deallocate ( cofz )
      deallocate ( cz )
      deallocate ( czz )
      deallocate ( gd ) 
      deallocate ( ipick )

      RETURN
      END   ! of mudslpm

c #####################################################################
c
c  subroutine smthlt1
c
c #####################################################################
c
      subroutine smthlt1 (nx,ny,nz,nor,a,b,f,fz)
c
      implicit none
c
      integer   ia,na
      integer   ix,jy,kz
      integer   id1,jd1,kd1
      parameter (id1=1,jd1=1,kd1=1)
      integer   ism,jsm,ksm
      integer   nx,ny,nz
      integer   nor
      real f,fz
c
      real      a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real      b(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
c #################################################################
      ism = 2
      jsm = 2
      ksm = 2
c
c copy array
c
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
      b(ix,jy,kz) = a(ix,jy,kz)
      end do
      end do
      end do
c
c  smooth
c
      do kz = 1+kd1,nz-ksm*kd1
      do jy = 1+jd1,ny-jsm*jd1
      do ix = 1+id1,nx-ism*id1
      a(ix,jy,kz) = b(ix,jy,kz)
     >  +f*(1.0/16.0)*(b(ix+id1,jy,kz)-2.*b(ix,jy,kz)+b(ix-id1,jy,kz))
     >  +f*(1.0/16.0)*(b(ix,jy+jd1,kz)-2.*b(ix,jy,kz)+b(ix,jy-jd1,kz))
     >  +fz*(1.0/16.0)*(b(ix,jy,kz+kd1)-2.*b(ix,jy,kz)+b(ix,jy,kz-kd1))
! new values  changed 9/12/2006 erm
!     >  +f*(1.0/30.0)*(b(ix+id1,jy,kz)-2.*b(ix,jy,kz)+b(ix-id1,jy,kz))
!     >  +f*(1.0/30.0)*(b(ix,jy+jd1,kz)-2.*b(ix,jy,kz)+b(ix,jy-jd1,kz))
!     > +fz*(1.0/100.0)*(b(ix,jy,kz+kd1)-2.*b(ix,jy,kz)+b(ix,jy,kz-kd1))
      end do
      end do
      end do
c
c  return
c
      return
      end
c
c
c
c
c-------------------------------------------------------------------------
c
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                 SUBROUTINE MAXFIELD                  ######
c     ######                written by E. Mansell                 ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE MXFLDSTRZ(nx,ny,nz,
     :    phi,gd,pstop,nbw,nbs,
     :    emax,emagmax,i1,j1,k1,i2,j2,k2,dx,dy,dz,istretch,
     :    z1d2lgt,nz1d,ezmax,iunit,bcx,bcy)
      
      implicit none
      
      integer nx,ny,nz
      integer iunit
      integer bcx,bcy
      real emag
      real phi(-1:nx+2,-1:ny+2,0:nz+1)
      real emax,emagmax,ezmax
      integer i1,j1,k1,i2,j2,k2
      integer jst,jend
      integer ix,jy,kz,i,j,k
      integer nbw,nbs
      real dx,dy,dz
      logical istretch
      real dxinv2,dyinv2,dzinv2
      real ex,ey,ez
      real gd(-1:1,-1:1,-1:1,0:nz,2)  !  1: distance, 2: inverse of distance
      real etest
      integer pstop
      integer nz1d
      real z1d2lgt(nz1d,4)
      integer, parameter :: ndebug = 0

c#######################################################################
c
c     Begin executable code
c
c#######################################################################
      
      i1 = 0
      j1 = 0
      k1 = 0
      i2 = 0 
      j2 = 0
      k2 = 0 
      emagmax = 0.0
      dxinv2 = 0.5/dx
      dyinv2 = 0.5/dy
      dzinv2 = 0.5/dz
      IF ( ny .gt. 2 ) THEN
       jst = 2
       jend = ny-1
      ELSE
       jst = 1
       jend = 1
      ENDIF
      
      ey = 0.0
      
      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: calculate E'

      DO k=1,nz-1
      DO j=jst,jend
      DO i=2,nx-1
         ex = (phi(i-1,j,k) - phi(i+1,j,k))*dxinv2
         IF ( ny .gt. 2 ) THEN
           ey = (phi(i,j-1,k) - phi(i,j+1,k))*dyinv2
         ENDIF
         IF ( istretch ) THEN
           ez = (phi(i,j,k-1) - phi(i,j,k+1))*dzinv2*z1d2lgt(k,3)
         ELSE
           ez = (phi(i,j,k-1) - phi(i,j,k+1))*dzinv2
         ENDIF
         emag = Sqrt(ex**2 + ey**2 + ez**2)
         IF ( emag .gt. emagmax ) THEN
           emagmax = emag
           ezmax = ez
           i1 = i
           j1 = j
           k1 = k
         END IF
      END DO  ! i
      END DO  ! j
      END DO  ! k
      

      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: check for max on boundary'

      IF (((bcx .ne. 2 .and. (i1 .lt. nbw+2 .or. i1 .gt. nx-nbw-2)) .or.
     :     (bcy .ne. 2 .and. ny .gt. 2 .and. (j1 .lt. nbs+2 .or. j1 .gt. ny-nbs-2)) .or.
     :      k1 .lt. 1 .or. k1 .gt. nz-2 ) .and.
     :   emagmax .gt. 60000.0   ) THEN
         write(iunit,*) 'Harumph! Emagmax along an edge! I dont like it!'
         write(0,*) 'Harumph! Emagmax along an edge! I dont like it!'
         write(iunit,'(a,3(1x,i4),1x,e12.5)') 'i1,j1,k,emagmax',
     :         i1,j1,k1,emagmax
         pstop = 1
      END IF
      
c  find direction of largest 1st order electric field
      IF ( ny .gt. 2 ) THEN ! 3D
        jst = -1
        jend = 1
      ELSE        ! 2D
        jst =  0
        jend = 0
      ENDIF
 
       IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: find direction of E'

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

      IF (ndebug .ge. 1) write(0,*) 'MXFLDSTRZ: return'
      
      RETURN
      END



