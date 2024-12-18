C#include "sam.def.h"
c#define ICE10
c#define SWM
c#define SAM
c
c Things to do:
c
c  Fix Rain evaporation for gamma function (ipconc >= 3)
c
c  convert cloud ice to snow as in Ferrier 1994 (change only mass in cloud ice),
c    then can try turning off direct conversion from cloud ice to graupel and rimed ice
c
c  look at an iterative check on overdepletion;  need to be careful with two-moment
c
c  check ice supersaturation in two-moment.  Getting enough deposition, or need 
c      to do sat adj. when cloud droplets are all gone?
c
c  
c
c new comment
c
c Fix use of gt for SWM IN FALLOUT ROUTINES
c
c  How to remove hl for ipconc=5?  Need to preprocess?
c
c   When the charging rates are moved to a subroutine, need to move the
c   call to be after the wet growth calculations -- or at least the 
c   splashing stuff.  Think about this....
c
c  Think about what to do with cracif
c
c    FIX BIGSAVE! (done)
c
c    Replace qv0 with qx(mgs,lv)? No. qv0 is base val
c
c Need to look at limiting supersaturation to 1 or so by nucleation/condensation
c
c  put in temperature-dependent function for homogeneous freezing
c
cc--------------------------------------------------------------------------
c
c
c--------------------------------------------------------------------------
c
#if !defined( ICE10 ) | defined( MPI ) | 1
      subroutine ice10q01
     >  (ntmul,nx,ny,nz,na,nba,nv,nstep)
c     >  ,elec,neelec)
#ifdef MPI
       write(0,*) 'STOP: 10-ICE is not supported in MPI yet'
#else
       write(0,*) 'STOP: must compile with -DICE10 to run this code and remove the hardwired 1 in the #if'
#endif
       STOP
       
       RETURN
       END
       

#else
      subroutine ice10q01
     >  (ntmul,nx,ny,nz,na,nba,nv,nstep
     >  ,nor,istag,jstag,kstag,itopo
     >  ,iwrite
     :  ,cwmasn,cwmasx
     >  ,dtp1,dxx,dyy,dzz,gx,gy,gz,dx,dy,dz,z1d
     >  ,xfall
     >  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,tq1
     >  ,ab,an,db,dn,p2,pinit
     <  ,pb,pn,u,v,w,iunit
     :  ,ssat,ssfilt,t00,t77,ssati
     >  ,elec,bcx,bcy)

       USE GRID_MODULE
       USE MICRO_MODULE
       USE ELEC_MODULE
       USE INDEX_MODULE
#ifdef MPI
       USE COMMASMPI_MODULE
#endif
c
c--------------------------------------------------------------------------
c                                
c     version sam.08.dd02a.f
c     detailed parameterized microphysics with 
c     1)  cloud water
c     2)  rain
c     3)  column ice 
c     4)  plate ice 
c     5)  rime ice 
c     6)  snow
c     7)  low-density graupel
c     8)  medium-density graupel
c     9)  high-density graupel
c     10) frozen drop
c     11) hail
c     12) large hail 
c
c     modifications for three graupel classes and ice crystal habits
c     note:   shedding by low density ice is off
c     note:   mzdist = 4
c     note:   only use ihole=1
c     note:   schaci and schacw are set to zero
c
c 2007.03.14: (APS) added array for mp process time series
c
c 2007.02.7: added check for frozen drop wet growth to turn off hallett-mossop
c
c 2006.05.17: Fixed bug when rescaling vapor deposition (e.g., qidpv )
c
c 2006.04.20: More changes to riming and rain interactions conversions:
c             1. Disallow riming conversion to frozen drops/hail except for
c                high dens graupel, and then only when wet growth is indicated
c             2. Keep snow-rain interaction the same as before (allow transfer
c                into FD but not hail).  For graupels, only convert to FD/H from
c                GH, and then only when wet growth occurs.
c
c 2006.04.17: Started converting inmicro params to namelists.  Need to initialize
c             CCN and ion densities in initsubs.f90.
c
c 2006.02.28: Added a subroutine (in sam.setvt.F) to set slope intercepts for 
c             single moment schemeto allow changes without recompiling 10ICE.
c             ALSO: turn off conversion of FD to hail and redirect conversion of 
c             GX -> FD to be GX -> H (i.e., riming/rain interaction conversion to
c             FD now goes to hail instead)
c
c 2005.12.7:  turned graupel/hail coll. of rain back on (exr) and redirected 
c             melting gl to gm and melting gm to gh.  Adjust qglmlr upward to get
c             density of gm.  (Likewise for qgmmlr to gh.)  Might be better to have
c             prediction of particle density.....
c
c 2005.11.30: Fixed setting chcngh to zero (was setting qhcngh twice).  
c             Now allow qhcngh if cwacghf > 1. (also for gm)
c             Now allow qhcnf if T < 253 and rimf > 850.
c
c 2005.11.11: (version "e") Fix qxacr/qracx with correct collection equations for the gamma
c             function used for rain.  Also change exr to be the multiple of the
c             fractions of x and r that are greater than 40 microns.
c
c 2005.11.04: Applied mass-weighted number sedimentation correction to number-weighted
c             sedimentation of concentration.  Also instituted reflectivity correction (rz)
c             as in Ferrier for freezing rain and melting graupel to correct for 
c             converting particles to a category with a different shape parameter 
c             or distribution.
c
c 2005.10.11:  Updates: created a subroutine (rainconv) for rain-ice density 
c              interaction analagous to rimconv.  Reduced rimconv to 21 mass bins
c              from 41.  Capacitance term on column ice was set to zero for 
c              diameters less than 10 microns!  Changed to check instead for
c              minimum mixing ratio and nonzero diameter.
c
c 2005.09.27: version nci.03a.f -- start trying to replace some of the GS arrays
c             with 2D arrays like qxw(ngs,lv:lhab)
c
c 2005.08.31: removed the test function theta2temp because it was not
c             accurate enough for determining supersaturation.
c             Also added a subroutine (qvexcess) that performs a saturation
c             adjustment to determine the 
c             amount of vapor that can be condensed to acheive a given
c             supersaturation.  It is used for initial development of cloud 
c             droplets (instead of the KW method) and for limiting the 
c             supersaturation to ssmx (currently set at 1.01 for 1% 
c             supersaturation)
c
c 2005.06.28: Fixed a couple small bugs:
c              cidia and cipdia were not defined for qci .le. qimin
c              could have a divide by zero with xvc = 0.0 in the Hallett-Mossop
c                 loop with ipconc .ge. 2 (fixed)
c              Now, only calculate ventillation coefficients if the mass is
c                above threshold.  Otherwise set to zero
c              Took out the recalculation of the cloud droplet radius in
c                 the riming conversions (just use 0.5*cwdia)
c
c 2005.05.06: turn of HM process for frozen drops, so only graupel does it now
c
c 2005.03:    Test switching from lookup tables to subroutine call for 
c             graupel rime conversions.
c
c 2005.02.23: Test turning off Bigg freezing for droplets smaller than
c             18 micron diameter when ipconc .ge. 2
c
c 2005.02.10: Added input variable eii1 for setting collection efficiencies
c             of graupel/hail collecting ice crystals/snow
c
c 2005.01.18: Fixed out-of-bounds bug for ipconc < 2 (for na < lss) in
c             the down-scaling code section that reduces ssmax out of cloud
c
c 2005.01.14: For ipconc < 2, set ccw to be cwccn*dnz(kz)/dnz(1) 
c             instead of just cwccn.
c
c 2004.12.31: Small fix for ipconc=0 where ccw was not set to cwccn.
c             (also introduced a bug for ipconc .ge. 2, fixed on 2005.01.17)
c
c 2004.11.22: version nci.02c:
c             Add explicit integration of droplet condensation in case of
c             an undershoot in the rate.
c
c 2004.10.25: try turning off ice crystal conc. limit on qiint/ciint
c             Also add variable qvimxd = 0.7*(qwv - qis)/dtp to use for
c             ice deposition instead of qvmxd, which uses (qwv - qvs)
c
c 2004.10.24: noticed that ssat can get large at T < 233 K in strong updrafts 
c             (perhaps ice crystal growth is not keeping up with temperature drop?)
c             So limit ice crystal nucleation (qiint) by setting the maximum ssi 
c             as t8/t9 (water sat/ice sat) for now.  Might need a sat. adj. scheme for this?
c
c             Also put in the concentration gain for rain from shedding (these had
c             been subtracted from the ice categories instead of being added to the rain)
c
c 2004.10.20: version sam.03.dd02a.nci.02a.f: move droplet condensation/nucleation
c             to a second gather-scatter loop in order to use supersaturation calculated
c             on the new time level.  (Instead of using ad to calculate SS.)
c
c 2004.10.05: Start adding cloud droplet concentration
c             Need to do concentration transfer for riming crystals to rimed ice!
c
c
c 2004.08.31:  version sam.08.dd02a.nci.01.f for adding ice crystal concentration
c              (for CI, IP, and IR)
c              Also set erw = 0 if rwdia .le. 50.0e-6
c
c 2004.03.30: new version (032) with new formula for snow-droplet collision 
c             efficiency (esw), based numerical data from Wang and Ji (1992)
c             in Pruppacher and Klett 1997, p. 601.
c             
c  2004.08.08: dd02a.nci.01.f -- switch lrtr,lrth,lrthw to number concentrations
c              for cloud droplets, ice plates, and ice columns.  Add 43rd variable
c              for rimed ice concentration.
c
c  2004.07.14:  Add option IF ( alfarim .eq. 0.0 ) THEN if riming rate is below
c               threshold set Tau to zero.  (Otherwise alfarim has to be a negative number)
c
c
c  2004.01.27: change option itimrimw=1 to use itim corresponding to rimtim
c              instead of just itim = 1 (rimtim of 20s)
c
c
c 10.17.2003  Fix potential bug in saturation adjustment scheme.  Only adjust
c             ice if qitmp(mgs) .gt. qimin instead of just adding 1.e-20 to qitmp(mgs).
c
c change from 030.f: replace amccw with qcw in the bulk rime density equation
c
c 05.13.2003 Clean up cloud ice and water by turning them into vapor
c            if they go below minimum.  Turn off qimlr and qipmlr (melt
c            into cloud water instead).  Threshold qirmlr to go to cloud
c            water if qir .le. 1.0e-6.  Increased qirmin from 1e-9 to 1e-7.
c
c 04.04.2003 Add option to compute lookup tables at runtime using
c            input value of cwccn (not done yet)
c
c 03.19.2003 Took out line that had
c                 qwvp(mgs) = max(qwvp(mgs), 0.0) 
c            which caused erroneous zeroing of the perturbation water vapor
c            Bug was found by Adam Houston in the 3-ICE code
c
c 03.14.2003 Small fix to do zeroing of precip variables when 
c            ipelec .lt. 2 (there was an IF statement of ipelec.ge.2
c            around the code for redistribution of minor amounts of 
c            precipitate).  The only real problem was minor differences
c            between runs with electriciy turned on vs. turned off.
c
c 11.12.2002 Add option nonigrd=-10 to turn off all charge separation
c            (both inductive and noninductive) but allow mass/charge 
c            transfer and ions with ipelec=4
c
c 10.16.2002 version dd02a.ion.2b1.f:
c            set charging to zero for ix,jy .le. 3 and ix,jy .ge. nx-3
c            Hope this will avoid large fields and flash rates when active
c            cells exit lateral boundaries.
c
c 8.30.2002 Fixed vtwbar (had an extraneous factor of rho0).  Thanks to M. Gilmore
c
c 8.28.2002 Added optional input for value of 'ewfac' along with
c           the input of 'nonigrd, rgard'.  'ewfac' has a default value
c           of 1.0 and multiplies the collection efficiency of cloud water
c           by graupel and hail.  This allows adjusting the actual rime
c           accretion rates for testing the effects on noninductive charging.
c
c 8.05.2002 H-M process: increased cimas0 to 2.5e-10 and switched to 
c           Ferrier (1994) parameterization of the temperature function.
c         
c           Also added new input variable for autoconversion threshold
c
c
c 7.22.2002 changed egms to have a factor of 0.01 instead of 0.1 to be consistent
c           with everything else
c
c 6.09.2002 Test increasing minimum cloud droplet mass from that of a 
c           1 micron radius to that of 5 micron radius (cwdia=10e-6).
c           Also recalculate ccw from cwmas in case so that the concentration
c           and mass are always consistent (as before changed for cloud ice)
c
c 5.24.2002  Added charge per particle limits on cloud ice and droplets
c
c 03/12/2001  Added input variable rgard to 'inmicro' which replaces the gf8*gf4i in 
c             the gardiner scheme.
c
c 09/05/2001  Added variable 'rgard1' as a multiplier for the gardiner
c             scheme for cloud ice (columns and plates)
c
c ----------------------------------------
c  3/11/99 took out ihabdo stuff from big jy loop
c  3/11/99 commented a bunch of unused arrays 
c  6/13/99 converted 2D arrays to be passed in so that no arrays are left
c          that will require the subroutine to be recompiled if array sizes
c          are changed.  Arrays all defined and passed in from sam.08.driver.sgi.f
c  6/23/99 add new version of Saunder's noninductive charging based on 
c          rime accretion rate.
c  7/6/99  add option for more liberal inductive mechanism (ipelec = 4)
c          which sets erbnd = 0.015 and costhe = 0.5 (as in Helsdon & Farley, 1987)
c          UPDATE: now ipelec=4 uses values read from 'inmicro'
c  7/9/99  completed addition of Takahashi (JAS, 1984) charging scheme (nonigrd=-1)
c  7/14/99 setting nsave = -1 will synchronize nsave with nhelec (dump intervals)
c  7/17/99 change inductive graupel/hail--cloud droplet charging to match
c          the formulation in Ziegler et al. (1991) equation (7). (Also changed
c          some ()**0.5 to Sqrt()
c  7/31/99 For gardiner charging, now limit low-temp charging rates by
c          a parabolic function with a value of zero at -30 C and goes to
c          zero at -43 C.
c
c 8/13/99  Start working on a subroutine to save lots of rates and temporary
c          array values.
c
c 9/7/99  test turning off inductive charging for wet growth (Already do
c         this for non-inductive charging) 
c
c 4/20/2000 Changed rain-snow multicomponent transform so that snow production
c           is put into low density graupel instead.
c
c           Also fixed charge transfer between cloud water and cloud ice for
c           melting and homogenous freezing in the saturation adjustment
c           section near the very end (it wasn't done before....)
c
c 8/15/2000 Fixed charge sink: when minor amounts of precip are transfered to
c           smaller categories, now transfer the charge, too. 
c
c 8/25/2000 made implicit none and commented out more unused temporary arrays
c
c 8/26/2000 reduced scalar array to remove the unused "mean cloud amount" arrays
c           (kept the one that is used, of course..)
c
c 8/29/2000 fixed a bug encountered when running on the Alpha box: if a mixing
c           ratio turned out to be exactly equal to the 'qxxmin' value, then
c           the mass would end up as zero (due to underflow because of very small
c           characteristic diameter) and then a divide-by-zero would occur.
c
c
c 2/13/2000 REMINDER OF STUFF TO DO:
c
c           put in code to turn off _inductive_ charging for T > 273 since
c           droplets can't even partially freeze to surface of graupel.
c          Also should check collection equations for consistency.  Maybe make
c          parameter statements for gamma function values with integer arguments
c           (or even all the commonly-used gamma function values).
c
c  **DONE (except evap.)
c        Also: track charge lost to evaporation and rain-out. (need to track loss 
c               from 'Wilson ion capture' in driver.)
c         
c      **DONE**   someday should bite the bullet and go to 'implicit none'
c
c         need to update for vertical stretching possibility (fallout, especially)
c
c
c   **DONE**  For 'mean cloud amount' only amccw is used.  So should be able to collapse 
c         the remaining arrays and save the memory and advection calculations.
c
c         Track charge that is zeroed-out when the hydrometeor mixing ratio is zero.
c
c
c  7/22/2001  After putting in some accounting checks for total charge, found
c             out that the precip fallout terms for charge were incorrectly 
c             formulated.  Should not multiply and divide by air density when
c             doing fall terms for charge!  Correcting this cured the charge 
c             non-conservation that resulted from the incorrect formulation.
c
c  8/16/2001  version sam.08.dd02a.lowmem.1b.f
c             Put in limits for inductive charging (to prevent runaway feedback)
c
c  1/8/2002  version sam.08.dd02a.lm.omp.2.f  
c            Changes to work with driver version sam.08.driver.lm.omp.2.f
c
c  3.29.2002  Removed the automatic rezeroing of hail-ice noninductive charge 
c             separation.  Should control this with the flags instead.
c
c  3.29.2002  Fixed bug in which Hallett-Mossop rates were calculated but not 
c             applied.  Now can turn types I and II on/off with flags itype1 and
c             itype2.  H-M is allowed for low/med/high graupel but not for 
c             hail/large hail/frozen drops.
c
c--------------------------------------------------------------------------
c
c  general declarations
c
c--------------------------------------------------------------------------
c
c
c
      implicit none
      
      
      TYPE(VARIABLE)     :: elec(neelec)
      
      integer            :: bcx, bcy
c
      integer ng1
      integer iunit,iunit0

      parameter(ng1 = 1)
!     parameter( iunit = 6 )

c
c      real timetd1,timetd2
c      real tarray(2),etime

      real,allocatable :: db1(:,:),db0(:,:)
c     real theta2temp
      real qvex
      
      character(LEN=100) line
      integer istat1
      
      logical, parameter :: lraincv = .true.
      
      integer iraincv, icgxconv
      parameter ( iraincv = 1, icgxconv = 1)
      real ffrz

      real qiptmp,qirtmp,qcitmp,cirdiatmp
      real ccwtmp,ciptmp,cirtmp,ccitmp
      real cpqc,cpip,cpir,cpci
      real cpqc0,cpip0,cpir0,cpci0
      real scfac,cpip1
      
      double precision dp1
      
!      real delqnw, delqxw
!      real tindmn,tindmx   ! min and max temperatures where inductive charging is allowed
      real frac
      
      real dtp1
      integer ntmul
      
c      real rar  ! rime accretion rate as calculated from qxacw
c extra variables for rimconv subroutine
      integer mr
      parameter(mr=5)
      real rfr(mr),  rfrh(mr),  rfrr(mr)
      real crfr(mr), crfrh(mr), crfrr(mr)


c a few vars for time-split fallout      
      real vtmax
      integer n,ndfall
      
      logical lsavetime  !  flag that it is time to save stuff (open a
                         !  file and call the save subroutine )
      character(LEN=80) savename
      integer   isaveunit,isaveunit2
      
      double precision chgneg,chgpos,sctot
      
      real temgtmp
      integer nstart, nstop, nstep
      integer nx,ny,nz,na,nba,nv
c      integer ng
      integer nor,nht,ngt,imapz,mzdist,igsr,istag,jstag,kstag,itopo
      integer iwrite
      real dtp,dx,dy,dz

      real dzc(nz)                         ! 1/dz(k)
      real dzz(nz),dxx(nx),dyy(ny)         ! dz(k),dx(i),dy(j)
      real gx(-nor+1:nx+nor)
      real gy(-nor+1:ny+nor)
      real gz(-nor+1:nz+nor)
      real z1d(-nor+1:nz+nor,4)

      real dv

      real dtptmp
c      integer nxl,nyl,nzl

      integer itest,nidx,id1,jd1,kd1
      parameter (itest=1)
      parameter (nidx=10)
      parameter (id1=1,jd1=1,kd1=1)
      integer ierr
      integer iend,imake,iread
      save imake
      integer ix,jy,kz, il, ic
      integer i, j, k
      integer i1,j1,k1
      real u1,v1,w1,w2
c      integer nxm,nym,nzm
      real eps,eps2
      parameter (eps=1.e-20,eps2=1.e-5)
c
c  electrical permitivity of air C / (N m**2) -  check the units
c
      real eperao
      parameter (eperao  = 8.8592e-12 )
      
      real ec,eci  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      parameter (eci = 1.0/ec)
c
c Takahashi lookup table
c
      integer nlwc,ntem
      parameter (nlwc=30, ntem=31)
      real takalu(0:ntem,0:nlwc) 
c
c  Other elec. vars
c
      real  temele
c      parameter (trever=-15.)  ! read it in instead

c
c   mixing ratio
c
c      integer lv,lc,lr,li,lir,ls
c      integer lgl,lgm,lgh,lf,lh
c      integer lip,lhl,lhab
c      parameter (lv=3,lc=4,lr=5,li=6,lir=7,ls=8)
c      parameter (lgl=9,lgm=10,lgh=11,lf=12,lh=13)
c      parameter (lip=14,lhl=15,lhab=15)
c
c  space charge (17 to 28)
c
c      integer lscw,lscr,lsci,lscir,lscs
c      integer lscgl,lscgm,lscgh,lscf,lsch
c      integer lscip,lschl
c      integer lscwi
c      integer lscpi,lscni
c      parameter ( lscwi = 17 )
c      parameter (lscw=0+lscwi,lscr=1+lscwi,lsci=2+lscwi)
c      parameter (lscir=3+lscwi,lscs=4+lscwi,lscgl=5+lscwi)
c      parameter (lscgm=6+lscwi,lscgh=7+lscwi,lscf=8+lscwi)
c      parameter (lsch=9+lscwi,lscip=10+lscwi,lschl=11+lscwi)
c      parameter (lscpi=12+lscwi,lscni=13+lscwi)

c      parameter (lscw=28,lscr=29,lsci=30,lscir=31,lscs=32)
c      parameter (lscgl=33,lscgm=34,lscgh=35,lscf=36,lsch=37)
c      parameter (lscip=38,lschl=39)
c
c  number concentration (was 16 to 27)
c
c      integer lnc,lnr,lni,lnir,lns
c      integer lngl,lngm,lngh,lnf
c      integer lnh,lnip,lnhl
c      integer lnshi
c      parameter (lnshi = 43)
c      parameter (lnc=0+lnshi,lnr=1+lnshi,lni=2+lnshi,lnir=3+lnshi)
c      parameter (lns=4+lnshi,lngl=5+lnshi,lngm=6+lnshi,lngh=7+lnshi)
c      parameter (lnf=8+lnshi,lnh=9+lnshi,lnip=10+lnshi,lnhl=11+lnshi)
c      parameter (lnc=32,lnr=44,lni=40,lnir=43)
c      parameter (lns=4+lnshi,lngl=5+lnshi,lngm=6+lnshi,lngh=7+lnshi)
c      parameter (lnf=8+lnshi,lnh=9+lnshi,lnip=42,lnhl=10+lnshi)
      
c
c Maximum supersaturation for a parcel
c
c      integer lss
c      parameter (lss = 45 )
c
c      integer lccn
c      parameter (lccn = 46 )
c
c  mean cloud content exposure (was 52 to 63)
c
c      integer lamcw ! ,lamcr,lamci,lamcir,lamcs
c      integer lamcgl,lamcgm,lamcgh,lamcf,lamch
c      integer lamcip,lamchl
c      integer lashi
c      parameter (lamcw=0+lashi) ! ,lamcr=1+lashi,lamci=2+lashi)
c      parameter (lamcir=3+lashi,lamcs=4+lashi,lamcgl=5+lashi)
c      parameter (lamcgm=6+lashi,lamcgh=7+lashi,lamcf=8+lashi)
c      parameter (lamch=9+lashi,lamcip=10+lashi,lamchl=11+lashi)
c
c  tau  (was 64 to 75)
c
      
      
!      integer lrtc,lrtr,lrti,lrtir,lrts
!      integer lrtgl,lrtgm,lrtgh,lrtf,lrth
!      integer lrtip,lrthl
!      integer lrshi

!      parameter (lrtc=0+lrshi) !,lrtr=1+lrshi,
!      parameter (lrti=1+lrshi,lrtip=2+lrshi)
!      parameter (lrtir=3+lrshi,lrts=4+lrshi,lrtgl=5+lrshi)
!      parameter (lrtgm=6+lrshi,lrtgh=7+lrshi,lrtf=8+lrshi)
c      parameter (lrth=10+lrshi,lrthl=11+lrshi)

      integer lrt(lc:lhab)
      integer lsc(lc:lhab)
      integer ln(lc:lhab)
      integer ipc(lc:lhab)
      
c
c  temperature (was 40 to 51)
c
      integer ltpc,ltpr,ltpi,ltpir,ltps
      integer ltpgl,ltpgm,ltpgh,ltpf,ltph
      integer ltpip,ltphl
      integer ltshi
      parameter (ltshi = 52)
      parameter (ltpc=0+ltshi,ltpr=1+ltshi,ltpi=2+ltshi)
      parameter (ltpir=3+ltshi,ltps=4+ltshi,ltpgl=5+ltshi)
      parameter (ltpgm=6+ltshi,ltpgh=7+ltshi,ltpf=8+ltshi)
      parameter (ltph=9+ltshi,ltpip=10+ltshi,ltphl=11+ltshi)
c
c habits
c
      integer nhab,nhabc,ihabdo
      parameter (nhab=10,nhabc=8,ihabdo=0)
      
c
c temporary arrays-self contained-sizes
c
!      integer iex,iey,iez,iemag,ipot
!      parameter (iex=1,iey=2,iez=3,iemag=4,ipot=5)
!      integer neelec
c
c  wind indicies
c
      integer mu,mv,mw
      parameter (mu=1,mv=2,mw=3)
c
c  conversion parameters
c
      integer mqcw,mqxw,mtem,mrho,mtim
      parameter (mqcw=21,mqxw=21,mtem=21,mrho=5,mtim=6)

      real xftim,xftimi,yftim, xftem,yftem, xfqcw,yfqcw, xfqxw,yfqxw
      parameter (xftim=0.05,xftimi = 1./xftim,yftim=1.)
      parameter (xftem=0.5,yftem=1.)
      parameter (xfqcw=2000.,yfqcw=1.)
      parameter (xfqxw=2000.,yfqxw=1.)
      
c
c  charge fallout arrays
c
      real xfall(nx,ny,na), xfalltot(nx,ny,na)
      real xfall0(nx,ny)
c      real gt0(-nor+ng1:nx+nor,-nor+ng1:1+nor,-nor+ng1:nz+nor,ngt)

c
c params read in from inmicro
c
      integer iptemp
      parameter ( iptemp = 0 )
      integer iptemp0

      real dtfac
      parameter ( dtfac = 1.0 )
      integer infdo
      
      save infdo

      

      real dtrim

      


      real cckm,ccne,ccnefac,cnexp
      save cckm,ccne,ccnefac,cnexp
c
      integer nsave
c
      integer ido(lc:lqmx)
      save ido
      
c      integer idocw, idorw, idoci, idoir, idoip, idosw
c
c      integer idogl, idogm, idogh, idofw, idohw, idohl
c
      integer iexy(lc:lqmx,lc:lqmx)
      save iexy
      
c      integer ieswi, ieswir, ieswip, ieswc, ieswr

c      integer ieglsw, iegli, ieglir, ieglip, ieglc, ieglr
c
c      integer iegmsw, iegmi, iegmir, iegmip, iegmc, iegmr
c
c      integer ieghsw, ieghi, ieghir, ieghip, ieghc, ieghr
c
c      integer iefwsw, iefwi, iefwir, iefwip, iefwc, iefwr
c
c      integer iehwsw, iehwi, iehwir, iehwip, iehwc, iehwr
c
c      integer iehlsw, iehli, iehlir, iehlip, iehlc, iehlr
c
!      real delqnsa, delqxsa, delqnsb, delqxsb, delqnia, delqxia
!      real delqnra, delqxra
      
!      real scippmx,scwppmx

c
c  arrays for transfer percents for ice-cloud & ice-rain 
c  multicomponent interactions
c
      real qfrwasw(mqcw,mqxw,mrho)
      real qfrwagl(mqcw,mqxw,mrho)
      real qfrwagm(mqcw,mqxw,mrho)
      real qfrwagh(mqcw,mqxw,mrho)
c
      real qfswarw(mqcw,mqxw,mrho)
      real qfglarw(mqcw,mqxw,mrho)
      real qfgmarw(mqcw,mqxw,mrho)
      real qfgharw(mqcw,mqxw,mrho)
c
c
c external temporary arrays
c
      real t00(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t77(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

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

      real tq1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)  !qproc temporary


      real pinit(-nor+1:nz+nor)  ! base state Pi
      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)  ! perturbation Pi
      real pb(-nor+1:nz+nor)
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real db(-nor+1:nz+nor)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real ab(-nor+1:nz+nor,na)
c      real an(nx,ny,nz,na)
c      real vn(-nor+1:ny+nor,-nor+1:nz+nor,-nor+1:nx+nor,nv)
      real u(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      real ssati(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      

      real ssat(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real ssfilt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
c 
c  declarations microphyscs and for gather/scatter
c
      integer nxmpb,nzmpb,nxz
      integer jgs,mgs,ngs,numgs
      parameter (ngs=500)
      integer ntt
      parameter (ntt=300)

      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs)
      integer nsvcnt
c      integer isave(ntt)
      integer ncuse
      parameter (ncuse=0)
      integer il0(ngs),il5(ngs)
c      integer il1m(ngs),il2m(ngs),il3m(ngs),il4m(ngs),il5m(ngs)
c
      integer jdx(nidx)
c
c
      real tdtol,temsav,tfrcbw,tfrcbi,thnuc
      
      real tfr, tfrh
      parameter ( tfr = 273.15, tfrh = 233.15 )
      
      real cp, rd
      parameter ( cp = 1004.0, rd = 287.04 )
      
      real cpi
      parameter ( cpi = 1./cp )
      
      real poo,cap
      parameter ( cap = rd/cp, poo = 1.0e+05 )

c      real tmxs(ntt),xmxs(ntt),xmns(ntt)
c
c  Ice Multiplication Arrays.
c
      real  fimt1(ngs),fimta(ngs),fimt2(ngs) !,qmul1(ngs),qmul2(ngs)
      real xcwmas
c
c  gamma function
c
c      integer ngm0,ngm1,ngm2
c      parameter (ngm0=1000,ngm1=500,ngm2=500)
c      real gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)
c
c Variables for Ziegler warm rain microphysics
c      

      real ccnc(ngs)
      real sscb  ! 'cloud base' SS threshold
      parameter ( sscb = 2.0 )
      integer idecss  ! flag to turn on (=1) decay of ssmax when no cloud or ice crystals
      parameter ( idecss = 1 )
      integer iba ! flag to do condensation/nucleation in 1st or 2nd loop
                  ! =0 to use ad to calculate SS
                  ! =1 to use an at end of main jy loop to calculate SS
      parameter (iba = 1)
      integer ifilt   ! =1 to filter ssat, =0 to set ssfilt=ssat
      parameter ( ifilt = 0 ) 
      real temp1,temp2,ssold
      real ssmax(ngs)       ! maximum SS experienced by a parcel
      real ssmx
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real bfnu
      real ventr, ventc
      save ventr, ventc
      parameter ( bfnu = (rnu + 2.0)/(rnu + 1.0)  )
      real volb,rb(ngs), t2s, aa1, aa2, xl2p(ngs)
      parameter ( aa1 = 9.44e15, aa2 = 5.78e3 ) ! a1 in Ziegler
      real rvt
      parameter ( rvt = 0.104 )
      real cautn(ngs), rh(ngs), nh(ngs)
      real ec0, ex1, ft, rhoinv(ngs)
      
      integer kbound
      real ac1,bc, taus, sstdy, super, c1,d1,e1,f1,p380,tmp,tmp1b,tmp1e,tmp2
      real tmp2b,tmp2e,qtmp1b,qtmp1e
      real d1r, d1i, d1s, e1i
      real c1sw   ! integration factor for snow melting with snu = -0.8
      save c1sw
      real, parameter :: vr1mm = 5.23599e-10 ! volume of 1mm diameter sphere (m**3)
      real  :: vmlt
      real rhosm
      parameter ( rhosm = 500. )
      integer nc ! condensation step
      real dtcon,dtcon1,dtcon2 ! condensation time step (dtcon*nc = dtp)
      real delta
      integer ltemq1,ltemq1m ! ,ltemq1m2
      real dqv,qv1,ss1,ss2,qvs1,dqvs,dtemp,dt1   ! temporaries for condensation
c      real  dtemp2,ss1m2
      real ssi1, ssi2, dqvi, dqvis, dqvii,qis1
      real dqvr, dqc, dqr, dqi, dqs
      real qv1m,qvs1m,ss1m,ssi1m,qis1m
      real cwmastmp 
      real as, bs, dcloud,dcloud2
      real cn(ngs) 
c      real xvc(ngs), xvr(ngs)
c      real xvs(ngs),xvgl(ngs),xvgm(ngs),xvgh(ngs),xvf(ngs)
c      real xvh(ngs),xvhl(ngs)
      real,allocatable :: xv(:,:)
      real mwfac
c      parameter ( mwfac = 6.0**(1./3.) ) 
c      ! factor for mass-weighted rain volume diameter
c      real wijk ! wvel
      real  es(ngs) ! ss(ngs),
      real ssbar(ngs) !, sqsat(ngs)
      real ssf(ngs),ssfkp1(ngs),ssfkm1(ngs),ssat0(ngs)
      real ssfjp1(ngs),ssfjm1(ngs)
      real ssfip1(ngs),ssfim1(ngs)
      real supcb, supmx
      parameter (supcb=0.5,supmx=238.0)
      real r2dxm, r2dym, r2dzm
      real dssdz, dssdy, dssdx
c      real cck
c      real xvcmn, xvcmx  ! min, max droplet volumes
c      real xvrmn, xvrmx  ! min, max rain volumes
c      real xvsmn, xvsmx  ! min, max snow volumes
c      real xvfmn, xvfmx  ! min, max frozen drop volumes
c      real xvgmn, xvgmx  ! min, max graupel volumes
c      real xvhmn, xvhmx  ! min, max hail volumes
c      real xvhlmn, xvhlmx  ! min, max lg hail volumes
      real xvmn(lc:lhab), xvmx(lc:lhab)
c      PARAMETER(XVCMN=4.188E-12,XVCMX=6.54E-8)   ! CGS
c      PARAMETER(XVRMN=2.8866E-7,XVRMX=4.1887E-3) ! CGS
c      parameter( xvcmn=4.188e-18 )   ! mks  min volume = 1 micron radius
c      parameter( xvcmx=6.54e-14 )    ! mks, 25 micron max radius
c      parameter( xvcmx=2.89e-13 )    ! mks, 41 micron max radius

c      data xvcmx/2.89e-13/
c      save xvcmx
c      parameter( xvrmn=2.8866e-13, xvrmx=0.523599*(3.e-3)**3 ) !( was 4.1887e-9 )  ! mks
c      parameter( xvsmn=2.8866e-13, xvsmx=0.523599*(6.e-3)**3 ) !( was 4.1887e-9 )  ! mks
c      parameter( xvfmn=0.523599*(0.1e-3)**3, xvfmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvgmn=0.523599*(0.1e-3)**3, xvgmx=0.523599*(6.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvhmn=0.523599*(1.e-3)**3, xvhmx=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
c      parameter( xvhlmn=0.523599*(1.e-3)**3, 
c     :           xvhlmx=0.523599*(40.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      real rwmasn,rwmasx
     
      real vgra,vfrz
      parameter ( vgra = 0.523599*(1.0e-3)**3 ) ! , vfrz = 0.523599*(0.5e-3)**3 )

      real epsi,d
      parameter (epsi = 0.622, d = 0.266)
      real r1,qevap ! ,slv
      
c      real svc(ngs)  !  droplet volume
c
c  misc
c
      real ni,nr,d0
      
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
      real delqci(ngs),delqip(ngs)
      real temp(ngs),tempc(ngs),tempbz(ngs),temp0(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) ! ,tembzg(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real elv(ngs),elf(ngs),els(ngs),elvs(ngs),elss(ngs)
      real gamw(ngs),gams(ngs)   !   qciavl(ngs),
      real tsqr(ngs),ssi(ngs),ssw(ngs)
      real cc3(ngs),cqv1(ngs),cqv2(ngs)
      real f5, prod0, qvs0  ! Kessler condensation factor
      real cwc1,cwc0
c      real qcwdif(ngs) ! ,dcwnc
      real cwnccn(ngs)
      real cwnc(ngs),qcwtmp(ngs),qtmp,ctmp,sctmp,qtot(ngs)
      real cwmasn,cwmasx
      real cwmasn5
      real cwradn
      real cimcnd(ngs) ! ,qcidif(ngs)  ! ,dcinc
      real cincmx(ngs),cincmn(ngs)
c      real cinccn(nz)
c      real cinc(ngs)    !  ,qcitmp(ngs)
      real cimasn,cimasx,ccimx
      real pi,pid4
      real ar,br,cs,ds,gf7,gf6,gf5,gf4,gf3,gf2,gf1
      real aradcw,bradcw,cradcw,dradcw,cwrad,rwrad
      
      real cionp(ngs),cionn(ngs)
c
c  other arrays
c
      
      real rargl(ngs),rargm(ngs),rargh(ngs)
      real rarhw(ngs),rarhl(ngs),rarfw(ngs)
      
      real fwet1(ngs),fwet2(ngs)   !   ,fwet3(ngs)
      real fmlt1(ngs),fmlt2(ngs)   !   ,fmlt3(ngs)
      real fvds(ngs),fvce(ngs),fiinit(ngs) ! ,fcinit(ngs)
      real fvent(ngs),fraci(ngs),fracl(ngs)
c
      real fai(ngs),fav(ngs),fbi(ngs),fbv(ngs)
      real felv(ngs),fels(ngs),felf(ngs)
      real felvs(ngs),felss(ngs)      !   ,felfs(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs)
      real fschm(ngs),fpndl(ngs)
      real fgamw(ngs),fgams(ngs)
      real fcqv1(ngs),fcqv2(ngs),fcc3(ngs)  
c
      real fcci(ngs), fcip(ngs)
c
       real qitmp(ngs)
      
      real qxmin(lc:lqmx)
      save qxmin
      
      real,allocatable :: qx(:,:) ! qx(ngs,lv:lhab)
      real,allocatable :: cx(:,:) ! cx(ngs,lc:lhab)
      real,allocatable :: scx(:,:) ! scx(ngs,lc:lhab)
c
      real swvent(ngs),fwvent(ngs),hwvent(ngs),rwvent(ngs)
      real glvent(ngs),gmvent(ngs),ghvent(ngs)  !  agvent(ngs),
      real civent(ngs),cirvent(ngs),cipvent(ngs)  !  cwvent(ngs),
      real hlvent(ngs)
c
c      real xivent(ngs,nhab),xix(ngs,nhab),xirey(ngs,nhab)
c      real xpvent(ngs)
c
      real rtcw(ngs),rtci(ngs),rtsw(ngs) ! ,rthw(ngs),rtrw(ngs)
      real rtgl(ngs),rtgh(ngs),rtfw(ngs)
      real rtgm(ngs),rtir(ngs) !,rthl(ngs)
      real rtip(ngs)
c
c      real sccw(ngs),scci(ngs),scrw(ngs),scsw(ngs),schw(ngs)
c      real scgl(ngs),scgh(ngs),scfw(ngs)
c      real scgm(ngs),schl(ngs),scir(ngs)
c      real scip(ngs)
c
      real tpcw(ngs),tpci(ngs),tprw(ngs),tpsw(ngs),tphw(ngs)
      real tpgl(ngs),tpgh(ngs),tpfw(ngs)
      real tpgm(ngs),tphl(ngs),tpir(ngs)
      real tpip(ngs)
c
c      real amccw(ngs) !,amcci(ngs),amcrw(ngs),amcsw(ngs),amchw(ngs)

      real rimf(ngs), rimgh(ngs)
c      real rimgl(ngs)
      real rimgm(ngs) 
c      real rimag(ngs), rimhl(ngs)
c
      real rhoglr(ngs) ! ,rhoglw(ngs)
      real rhogmw(ngs), rhogmr(ngs)
      real rhoghw(ngs), rhoghr(ngs)
c
      real fsczz(ngs),fscir(ngs),fscip(ngs)
      real fschw(ngs),fscfw(ngs),fscsw(ngs)
      real fsccw(ngs),fscci(ngs),fscrw(ngs)
      real fscgl(ngs),fscgm(ngs),fscgh(ngs)
      real fschl(ngs)
c
      
      real, allocatable :: vtxbar(:,:,:)
      real, allocatable :: xmas(:,:)
      real, allocatable :: xdn(:,:)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real, allocatable :: xdia(:,:,:)
c      real vtxbar(ngs,nhab),ximas(ngs,nhab),xidn(ngs,nhab)
c
c
c
c
c
c
c
c
c
c
c
      real cilen(ngs) ! ,ciplen(ngs)
c
c
      real rwcap(ngs),swcap(ngs),fwcap(ngs)  !  ,cwcap(ngs)
      real glcap(ngs),gmcap(ngs),ghcap(ngs),hwcap(ngs)
      real cicap(ngs),circap(ngs),cipcap(ngs),hlcap(ngs)

      real rimimas(ngs),rimirmas(ngs),rimipmas(ngs)
      real rimifrc(ngs),rimirfrc(ngs),rimipfrc(ngs)
c
      real qirmxd(ngs),qipmxd(ngs),qvimxd(ngs)
      real qimxd(ngs),qcmxd(ngs),qrmxd(ngs),qsmxd(ngs),qhmxd(ngs)
      real qglmxd(ngs),qghmxd(ngs),qfmxd(ngs)  !  ,qxmxd(ngs,nhab)
      real qgmmxd(ngs),qhlmxd(ngs)  !  ,qhxmxd(ngs),qagmxd(ngs)
c
      real cirmxd(ngs),cipmxd(ngs)  !  ,cvmxd(ngs)
      real cimxd(ngs),ccmxd(ngs),crmxd(ngs),csmxd(ngs),chmxd(ngs)
      real cglmxd(ngs),cghmxd(ngs),cfmxd(ngs)  !  ,cxmxd(ngs,nhab)
      real cgmmxd(ngs)  !  ,chxmxd(ngs),cagmxd(ngs),chlmxd(ngs)
c
      real pamccw(ngs) ! ,pamcci(ngs),pamcrw(ngs),pamcsw(ngs),pamchw(ngs)

c
      real prtcw(ngs),prtci(ngs),prtsw(ngs) ! ,prthw(ngs) ! ,prtrw(ngs)
      real prtgl(ngs),prtgh(ngs),prtfw(ngs)  !  ,prtxi(ngs,nhab)
      real prtgm(ngs),prtir(ngs)  ! ,prthl(ngs) ,prtag(ngs),prthx(ngs)
      real prtip(ngs)
c
c
c  Da electricity array work space...
c
c
      real sctem(ngs)
      real dellwc(ngs)
c
c  space charge arrays.....
c
      real dezcomp(ngs)
c
c
      real scsacw(ngs) 
      

      real scsacr(ngs)  
      real scsaci(ngs)  
      real sxsaci(ngs)  
      real scsacir(ngs),  scsacip(ngs)
      real sxsacir(ngs),  sxsacip(ngs)
c
      real scglacw(ngs) 

      real scglacr(ngs)  !  ,   scracgl(ngs)

      real scglaci(ngs),   scglacs(ngs)
      real sxglaci(ngs),   sxglacs(ngs)
      real scglacir(ngs),  scglacip(ngs)
      real sxglacir(ngs),  sxglacip(ngs)
c
      real scgmacw(ngs)  !  ,   scwacgm(ngs)
c      real sxgmacw(ngs)  !  ,   sxwacgm(ngs)

      real scgmacr(ngs)  !  ,   scracgm(ngs)
c      real sxgmacr(ngs)  !  ,   sxracgm(ngs)

      real scgmaci(ngs),   scgmacs(ngs)
      real sxgmaci(ngs),   sxgmacs(ngs)
      real scgmacir(ngs),  scgmacip(ngs)
      real sxgmacir(ngs),  sxgmacip(ngs)
c
      real scghacw(ngs)  !  ,   scwacgh(ngs)
c      real sxghacw(ngs)  !  ,   sxwacgh(ngs)

      real scghacr(ngs)  !  ,   scracgh(ngs)
c      real sxghacr(ngs)  !  ,   sxracgh(ngs)

      real scghaci(ngs),   scghacs(ngs)
      real sxghaci(ngs),   sxghacs(ngs)
      real scghacir(ngs),  scghacip(ngs)
      real sxghacir(ngs),  sxghacip(ngs)
c
      real scfacw(ngs)  !  ,   scwacf(ngs)
c      real sxfacw(ngs)  !  ,   sxwacf(ngs)

      real scfacr(ngs)  !  ,   scracf(ngs)
c      real sxfacr(ngs)  !  ,   sxracf(ngs)

      real scfaci(ngs),   scfacs(ngs)
      real sxfaci(ngs),   sxfacs(ngs)
      real scfacir(ngs),  scfacip(ngs)
      real sxfacir(ngs),  sxfacip(ngs)
c
      real schacw(ngs), schacr(ngs)
c      real sxhacw(ngs) !, sxhacr(ngs)
      real schaci(ngs), schacs(ngs)
      real sxhaci(ngs), sxhacs(ngs)
      real schacir(ngs), schacip(ngs)
      real sxhacir(ngs), sxhacip(ngs)
c
      real schlacw(ngs), schlacr(ngs)
      real schlaci(ngs), schlacs(ngs)
      real sxhlaci(ngs), sxhlacs(ngs)
      real schlacir(ngs), schlacip(ngs)
      real sxhlacir(ngs), sxhlacip(ngs)
c
c
      real chmul1(ngs),chlmul1(ngs),cfmul1(ngs),csmul1(ngs)
      real cghmul1(ngs),cgmmul1(ngs),cglmul1(ngs)
      real qhmul1(ngs),qhlmul1(ngs),qfmul1(ngs),qsmul1(ngs)
      real qghmul1(ngs),qgmmul1(ngs),qglmul1(ngs)
c
c
c  concentration arrays...
c
c
      real cwacigl(ngs) !, cwacigm(ngs), cwacigh(ngs)
      real cwacipgl(ngs) ! , cwacipgm(ngs), cwacipgh(ngs)

      real chcngh(ngs), chcngm(ngs), chcngl(ngs), chcnf(ngs)
      real chlcnh(ngs)
      real cracif(ngs), ciacrf(ngs)
c      real cracipf(ngs), cipacrf(ngs)
c      real  cipacip(ngs) !  
      real cracr(ngs)
c      real cwacw(ngs), cfacf(ngs)
c      real cglacgl(ngs), cgmacgm(ngs)
      real  ciacip(ngs), cipacir(ngs)  !  ,cghacgh(ngs)
c
      real ciint(ngs), crfrz(ngs), crfrzf(ngs)
      real cicint(ngs)  !  , ciracir(ngs), ciaci(ngs)
      real ciacir(ngs), ciraci(ngs)
      real cipint(ngs) !, cipacwi(ngs)
      real ciracip(ngs)  !  , ciihrp(ngs)
c
      real ciacw(ngs), cwaci(ngs), cwacii(ngs)
      real ciacr(ngs), craci(ngs)
      real ciracw(ngs), cwacir(ngs)
c      real ciracr(ngs)
      real cracir(ngs)
      real cipacw(ngs), cwacip(ngs)
c      real cipacr(ngs), 
      real cracip(ngs)
      real cipaci(ngs)
c
      real csacw(ngs),   cwacs(ngs)
c      real csacwgl(ngs)
      real cwacsgl(ngs)
c      real csacwgm(ngs)
      real cwacsgm(ngs)
c      real csacwgh(ngs)
      real cwacsgh(ngs)
c      real csacwf(ngs)
      real cwacsf(ngs)
      real csacr(ngs),   cracs(ngs)
c      real csacrgl(ngs)
      real cracsgl(ngs)
c      real csacrgm(ngs)
      real cracsgm(ngs)
c      real csacrgh(ngs)
      real cracsgh(ngs)
      real csacrf(ngs),  cracsf(ngs)
      real csaci(ngs),   csacs(ngs)
      real csacir(ngs),  ciracs(ngs)
      real csacip(ngs),  cipacs(ngs)
      real ciacs(ngs)
c
      real cglacw(ngs),   cwacgl(ngs)
c      real cglacwgl(ngs)  !  , cwacglgl(ngs)
c      real cglacwgm(ngs)
      real cwacglgm(ngs)
c      real cglacwgh(ngs)
      real cwacglgh(ngs)
c      real cglacwf(ngs)
      real cwacglf(ngs)
      real cglacr(ngs),   cracgl(ngs)
c      real cglacrgl(ngs)  !  , cracglgl(ngs)
c      real cglacrgm(ngs)
      real cracglgm(ngs)
c      real cglacrgh(ngs)
      real cracglgh(ngs)
c      real cglacrf(ngs)
      real cracglf(ngs)
      real cglaci(ngs),   cglacs(ngs)
      real cglacir(ngs)  !  ,  ciracgl(ngs)
      real cglacip(ngs)  !  ,  cipacgl(ngs)
c
      real cgmacw(ngs),   cwacgm(ngs)
c      real cgmacwgl(ngs)
      real cwacgmgl(ngs)
c      real cgmacwgm(ngs)  !  , cwacgmgm(ngs)
c      real cgmacwgh(ngs)
      real cwacgmgh(ngs)
c      real cgmacwf(ngs)
      real cwacgmf(ngs)
      real cgmacr(ngs),   cracgm(ngs)
c      real cgmacrgl(ngs)
      real cracgmgl(ngs)
c      real cgmacrgm(ngs)  !  , cracgmgm(ngs)
c      real cgmacrgh(ngs)
      real cracgmgh(ngs)
c      real cgmacrf(ngs)
      real cracgmf(ngs)
      real cgmaci(ngs),   cgmacs(ngs)
      real cgmacir(ngs)  !  ,  ciracgm(ngs)
      real cgmacip(ngs)  !  ,  cipacgm(ngs)
c
      real cghacw(ngs),   cwacgh(ngs)
c      real cghacwgl(ngs)
      real cwacghgl(ngs)
c      real cghacwgm(ngs)
      real cwacghgm(ngs)
c      real cghacwgh(ngs)  !  , cwacghgh(ngs)
c      real cghacwf(ngs)
      real cwacghf(ngs)
      real cghacr(ngs),   cracgh(ngs)
c      real cghacrgl(ngs)
      real cracghgl(ngs)
c      real cghacrgm(ngs)
      real cracghgm(ngs)
c      real cghacrgh(ngs)  !  , cracghgh(ngs)
c      real cghacrf(ngs)
      real cracghf(ngs)
      real cghaci(ngs),   cghacs(ngs)
      real cghacir(ngs)  !  ,  ciracgh(ngs)
      real cghacip(ngs)  !  ,  cipacgh(ngs)
c
      real cfacw(ngs),   cwacf(ngs)
c      real cfacwgl(ngs) 
      real cwacfgl(ngs)
c      real cfacwgm(ngs)
      real cwacfgm(ngs)
c      real cfacwgh(ngs)
      real cwacfgh(ngs)
c      real cfacwf(ngs)  !  ,  cwacff(ngs)
      real cfacr(ngs),   cracf(ngs)
c      real cfacrgl(ngs), cracfgl(ngs)
c      real cfacrgm(ngs), cracfgm(ngs)
c      real cfacrgh(ngs), cracfgh(ngs)
c      real cfacrf(ngs)  !  ,  cracff(ngs)
      real cfaci(ngs),   cfacs(ngs)
      real cfacir(ngs)  !  ,  ciracf(ngs)
      real cfacip(ngs)  !  ,  cipacf(ngs)
c
      real cwacr(ngs), cracw(ngs)
      real chacw(ngs), chacr(ngs)
      real chaci(ngs), chacs(ngs)
      real chacir(ngs),chacip(ngs)
c
      real chlacw(ngs), chlacr(ngs)
      real chlaci(ngs), chlacs(ngs)
      real chlacir(ngs),chlacip(ngs)
c
c      real csacx(ngs,nhab)
c
c      real ciacwir(ngs)
      real cwaciir(ngs)
      real ciacwgl(ngs)  !  , cwacrgl(ngs)
      real ciracwgl(ngs), cwacirgl(ngs)  !  , csacipgh(ngs)
      real cwacirir(ngs) !,ciracwir(ngs) ! ,cirirhr(ngs)
      real cwacipir(ngs),cipacwir(ngs)
      real cwacipip(ngs),cipacwip(ngs)
      real crcnw(ngs) ! ,ciacwi(ngs)
      real cidpv(ngs),cisbv(ngs)
      real cimlr(ngs)

      real cfdpv(ngs),cfsbv(ngs)
      real cfmlr(ngs)
      real cfshr(ngs)

      real chlsbv(ngs)
      real chlmlr(ngs) ! ,chlcev(ngs),chldsv(ngs)
      real chlshr(ngs) ! chlwet(ngs),chldry(ngs),

      real chdpv(ngs),chsbv(ngs) ! chcnv(ngs),chevv(ngs),
      real chmlr(ngs) ! ,chcev(ngs),chdsv(ngs)
      real chshr(ngs) ! chwet(ngs),chdry(ngs),
      
      real chmlrr(ngs), cgmmlrr(ngs), cghmlrr(ngs)
      real cglmlrr(ngs), cfmlrr(ngs), chlmlrr(ngs)
      
      real csdpv(ngs),cssbv(ngs) ! cscnv(ngs),csevv(ngs),
      real csmlr(ngs) ! ,cscev(ngs),csdsv(ngs)
      real csshr(ngs) ! cswet(ngs),csdry(ngs),
c      real csshrp(ngs)
      real cipdpv(ngs),cipsbv(ngs) ! cipcnv(ngs),cipevv(ngs),
      real cipmlr(ngs) ! ,cipcev(ngs),cipdsv(ngs)
      real cipshr(ngs) ! cipwet(ngs),cipdry(ngs),

      real cirdpv(ngs),cirsbv(ngs) ! circnv(ngs),cirevv(ngs),
      real cirmlr(ngs),cirmlw(ngs) ! ,circev(ngs),cirdsv(ngs)
      real cirshr(ngs) ! cirwet(ngs),cirdry(ngs),

      real cgldpv(ngs),cglsbv(ngs) ! cglcnv(ngs),cglevv(ngs),
      real cglmlr(ngs) ! ,cglcev(ngs),cgldsv(ngs)
      real cglshr(ngs) ! cglwet(ngs),cgldry(ngs),

      real cgmdpv(ngs),cgmsbv(ngs) ! cgmcnv(ngs),cgmevv(ngs),
      real cgmmlr(ngs) ! ,cgmcev(ngs),cgmdsv(ngs)
      real cgmshr(ngs) ! cgmwet(ngs),cgmdry(ngs),
c      real cgmshrp(ngs)
      real cghdpv(ngs),cghsbv(ngs) ! cghcnv(ngs),cghevv(ngs),
      real cghmlr(ngs) ! ,cghcev(ngs),cghdsv(ngs)
      real cghshr(ngs) ! cghwet(ngs),cghdry(ngs),
c      real cghshrp(ngs)
c      real crztot(ngs),crzmax(ngs),crzfac(ngs)
      real crcev(ngs) ! ,crmlr(ngs)
      real crshr(ngs)
c
c
c arrays for w-ac-x ;  x-ac-w
c
c
      real qsacl(ngs),qglacl(ngs),qgmacl(ngs),qghacl(ngs)
      real qiacl(ngs),qfacl(ngs),qhacl(ngs)  ! qaacl(ngs),
      real qhlacl(ngs)  !  ,xacl(ngs,nhab),qiracl(ngs),qipacl(ngs)
c
      real qrcnw(ngs)
c      real qwacirc(ngs) ! qiracwc(ngs),
c      real qiaci(ngs), qsacs(ngs), 
      real qiracir(ngs) ! , qwacw(ngs)
      real qipacwir(ngs)
      real qiracwir(ngs),qwacirir(ngs)
      real qwacirgl(ngs),qiracwgl(ngs)
c      real qwacirgm(ngs),qiracwgm(ngs)
c      real qwacirgh(ngs),qiracwgh(ngs)
      real qwacipgl(ngs),qipacwgl(ngs)
c      real qwacipgm(ngs),qipacwgm(ngs)
c      real qwacipgh(ngs),qipacwgh(ngs)
c
      real qiacwi(ngs),qwacii(ngs)
      real qipacwip(ngs),qwacipip(ngs)
      real qiacwgl(ngs),qwacigl(ngs)
c      real qiacwgm(ngs),qwacigm(ngs)
c      real qiacwgh(ngs),qwacigh(ngs)
      real qiacwir(ngs),qwaciir(ngs)
c      real qiacwip(ngs),qipacip(ngs)
      real qwacipir(ngs)
c      real tqiacw(ngs),tqwaci(ngs)
c      real tqipacw(ngs),tqwacip(ngs)
      real qracw(ngs) ! qwacr(ngs),
      real qwaci(ngs),qiacw(ngs)
      real qwacip(ngs),qipacw(ngs)
      real qwacir(ngs),qiracw(ngs)
c      real qwaca(ngs),qaacw(ngs)
      real qwacs(ngs),qsacw(ngs)
      real qhacw(ngs) ! qwach(ngs),
      real qwacf(ngs),qfacw(ngs)
c      real qwachx(ngs),qhxacw(ngs)
      real qhlacw(ngs) ! qwachl(ngs),
      real qwacgh(ngs),qghacw(ngs)
      real qwacgm(ngs),qgmacw(ngs)
      real qwacgl(ngs),qglacw(ngs)
c      real qwacx(ngs,nhab),qxacw(ngs,nhab)
c      real qxacwx(ngs,nhab)
c
      real qwacss(ngs),qsacws(ngs)
c      real qwacsa(ngs),qsacwa(ngs)
      real qwacsf(ngs),qsacwf(ngs)
c      real qwacsh(ngs),qsacwh(ngs)
      real qwacsgl(ngs),qsacwgl(ngs)
      real qwacsgm(ngs),qsacwgm(ngs)
      real qwacsgh(ngs),qsacwgh(ngs)
c      real qwacshx(ngs),qsacwhx(ngs)
c      real qwacshl(ngs),qsacwhl(ngs)
c
      real qwacgmf(ngs),qgmacwf(ngs)
c      real qwacgmh(ngs),qgmacwh(ngs)
      real qwacgmgl(ngs),qgmacwgl(ngs)
      real qwacgmgm(ngs),qgmacwgm(ngs)
      real qwacgmgh(ngs),qgmacwgh(ngs)      
c
      real qwacglf(ngs),qglacwf(ngs)
c      real qwacglh(ngs),qglacwh(ngs)
      real qwacglgl(ngs),qglacwgl(ngs)
      real qwacglgm(ngs),qglacwgm(ngs)
      real qwacglgh(ngs),qglacwgh(ngs)      
c
      real qwacghf(ngs),qghacwf(ngs)
c      real qwacghh(ngs),qghacwh(ngs)
      real qwacghgl(ngs),qghacwgl(ngs)
      real qwacghgm(ngs),qghacwgm(ngs)
      real qwacghgh(ngs),qghacwgh(ngs)      
c
      real qwacff(ngs),qfacwf(ngs)
      real qwacfh(ngs),qfacwh(ngs)
      real qwacfgl(ngs),qfacwgl(ngs)
      real qwacfgm(ngs),qfacwgm(ngs)
      real qwacfgh(ngs),qfacwgh(ngs)      

c
c  arrays for x-ac-r and r-ac-x; 
c
      real qfacr(ngs) ! ,qracf(ngs)
c      real qaacr(ngs),qraca(ngs)
      real qsacr(ngs),qracs(ngs)
      real qhacr(ngs) ! ,qrach(ngs)
      real qiacr(ngs),qraci(ngs)
      real qipacr(ngs),qracip(ngs)
      real qiracr(ngs),qracir(ngs)
      real qhlacr(ngs) ! ,qrachl(ngs)
c      real qhxacr(ngs),qrachx(ngs)
      real qglacr(ngs),qracgl(ngs)
      real qgmacr(ngs),qracgm(ngs)
      real qghacr(ngs),qracgh(ngs)
c
      real qracif(ngs),qiacrf(ngs)

c
      real qracss(ngs),qsacrs(ngs)
c      real qracsh(ngs),qsacrh(ngs)
      real qracsf(ngs),qsacrf(ngs)
      real qracsgl(ngs),qsacrgl(ngs)
      real qracsgm(ngs),qsacrgm(ngs)
      real qracsgh(ngs),qsacrgh(ngs)      
c
      real qracglf(ngs),qglacrf(ngs)
c      real qracglh(ngs),qglacrh(ngs)
      real qracglgl(ngs),qglacrgl(ngs)
      real qracglgm(ngs),qglacrgm(ngs)
      real qracglgh(ngs),qglacrgh(ngs)      
c
      real qracgmf(ngs),qgmacrf(ngs)
c      real qracgmh(ngs),qgmacrh(ngs)
      real qracgmgl(ngs),qgmacrgl(ngs)
      real qracgmgm(ngs),qgmacrgm(ngs)
      real qracgmgh(ngs),qgmacrgh(ngs)      
c
      real qracghf(ngs),qghacrf(ngs)
c      real qracghh(ngs),qghacrh(ngs)
      real qracghgl(ngs),qghacrgl(ngs)
      real qracghgm(ngs),qghacrgm(ngs)
      real qracghgh(ngs),qghacrgh(ngs)      
c
      real qfacrf(ngs) ! qracff(ngs),
c
c  ice - ice interactions
c
      real qsaci(ngs),qsacir(ngs)
      real qiracs(ngs)  !  ,qsacx(ngs,nhab)
c
c      real qhacx(ngs,nhab)
      real qhaci(ngs),qhlaci(ngs)  !  ,qhxaci(ngs)
      real qhacip(ngs),qhlacip(ngs)  !  ,qhxacip(ngs)
      real qhacs(ngs),qhlacs(ngs)  !  ,qhxacs(ngs)
c
      real qfaci(ngs) ! ,qiacf(ngs)
      real qfacip(ngs) ! ,qipacf(ngs)
      real qfacs(ngs) ! ,qsacf(ngs)
c
      real qgmaci(ngs) ! ,qiacgm(ngs)
      real qgmacip(ngs) ! ,qipacgm(ngs)
      real qgmacs(ngs) ! ,qsacgm(ngs)
c
      real qglaci(ngs) ! ,qiacgl(ngs)
      real qglacip(ngs) ! ,qipacgl(ngs)
      real qglacs(ngs) ! ,qsacgl(ngs)
c
      real qghaci(ngs) ! ,qiacgh(ngs)
      real qghacip(ngs) ! ,qipacgh(ngs)
      real qghacs(ngs) ! ,qsacgh(ngs)
c
      real qipaci(ngs),qiacip(ngs)
      real qiraci(ngs),qiacir(ngs)
      real qiracip(ngs),qipacir(ngs)
      real qhacir(ngs) ! qirach(ngs),
      real qhlacir(ngs) ! qirachl(ngs),
      real qfacir(ngs) ! qiracf(ngs),
      real qghacir(ngs) ! qiracgh(ngs),
      real qgmacir(ngs) ! qiracgm(ngs),
      real qglacir(ngs) ! qiracgl(ngs),
c
c  conversions
c
      real qrfrz(ngs), qrfrzf(ngs)
      real qwfrz(ngs), qwctfz(ngs)
      real cwfrz(ngs), cwctfz(ngs)
      real qwfrzc(ngs), qwctfzc(ngs)
      real cwfrzc(ngs), cwctfzc(ngs)
      real qwfrzp(ngs), qwctfzp(ngs)
      real cwfrzp(ngs), cwctfzp(ngs)
      real xcolmn(ngs), xplate(ngs)
      real ciihr(ngs), qiihr(ngs)
      real cicichr(ngs), qicichr(ngs)
      real cipiphr(ngs), qipiphr(ngs)
      real qscni(ngs), cscni(ngs)
      real qscnvi(ngs), cscnvi(ngs)
      real qscnir(ngs),cscnir(ngs)
      real qscnip(ngs),cscnip(ngs)
c      real qscnx(ngs,nhab)
c      real qhcns(ngs)
      real qiint(ngs),qipipnt(ngs),qicicnt(ngs)
c      real qsfw(ngs),qsfi(ngs)
      real timflg(ngs)
c      real ssifac(ngs)
      real cninm(ngs),cnina(ngs),cninp(ngs),wvel(ngs),wvelkm1(ngs)
      real uvel(ngs),vvel(ngs)
c
      real qidpv(ngs),qisbv(ngs) ! qicnv(ngs),qievv(ngs),
      real qimlr(ngs),qidsv(ngs),qidsvp(ngs) ! ,qicev(ngs)

c
      real qfdpv(ngs),qfsbv(ngs) ! qfcnv(ngs),qfevv(ngs),
      real qfmlr(ngs),qfdsv(ngs) ! ,qfcev(ngs)
      real qfwet(ngs),qfdry(ngs),qfshr(ngs)
      real qfshrp(ngs)
c
      real qhldpv(ngs),qhlsbv(ngs) ! qhlcnv(ngs),qhlevv(ngs),
      real qhlmlr(ngs),qhldsv(ngs) ! ,qhlcev(ngs)
      real qhlwet(ngs),qhldry(ngs),qhlshr(ngs)
      real qhlshrp(ngs)
c
      real qhdpv(ngs),qhsbv(ngs) ! qhcnv(ngs),qhevv(ngs),
      real qhmlr(ngs),qhdsv(ngs) ! ,qhcev(ngs)
      real qhwet(ngs),qhdry(ngs),qhshr(ngs)
      real qhshrp(ngs)
c
      real qsdpv(ngs),qssbv(ngs) ! qscnv(ngs),qsevv(ngs),
      real qsmlr(ngs),qsdsv(ngs) ! ,qscev(ngs)
      real qswet(ngs),qsdry(ngs),qsshr(ngs)
      real qsshrp(ngs)
c
c
      real qipdpv(ngs),qipsbv(ngs) ! qipcnv(ngs),qipevv(ngs),
      real qipmlr(ngs),qipdsv(ngs) ! ,qipcev(ngs)
c      real qipshr(ngs) ! qipwet(ngs),qipdry(ngs),
c      real qipshrp(ngs)
c
      real qirdpv(ngs),qirsbv(ngs) ! qircnv(ngs),qirevv(ngs),
      real qirmlr(ngs),qirdsv(ngs),qirmlw(ngs)  ! ,qircev(ngs)
c      real qirshr(ngs) ! qirwet(ngs),qirdry(ngs),
c      real qirshrp(ngs)
c
      real qgldpv(ngs),qglsbv(ngs) ! qglcnv(ngs),qglevv(ngs),
      real qglmlr(ngs),qgldsv(ngs) ! ,qglcev(ngs)
      real qglwet(ngs),qgldry(ngs),qglshr(ngs)
      real qglshrp(ngs)
c
      real qgmdpv(ngs),qgmsbv(ngs) ! qgmcnv(ngs),qgmevv(ngs),
      real qgmmlr(ngs),qgmdsv(ngs) ! ,qgmcev(ngs)
      real qgmwet(ngs),qgmdry(ngs),qgmshr(ngs)
      real qgmshrp(ngs)
c
      real qghdpv(ngs),qghsbv(ngs) ! qghcnv(ngs),qghevv(ngs),
      real qghmlr(ngs),qghdsv(ngs) ! ,qghcev(ngs)
      real qghwet(ngs),qghdry(ngs),qghshr(ngs)
      real qghshrp(ngs)
c
      real qrztot(ngs),qrzmax(ngs),qrzfac(ngs)
      real qrcev(ngs)
      real qrshr(ngs)
c
c      real ffglwg(ngs),ffgmwg(ngs),ffghwg(ngs),ffswwg(ngs)
c      real ffhwwg(ngs),ffagwg(ngs),fffwwg(ngs),ffciwg(ngs)
c
      real qhcnf(ngs) ! ,qhcnhl(ngs),qhlcnhx(ngs)
      real qhlcnh(ngs) ! qhxcnhl(ngs),
      real qhcngh(ngs),qhcngm(ngs),qhcngl(ngs)
c
c      real exwidia(nhab),exwwdia(nhab)

      real eiw(ngs),eii(ngs),eiri(ngs),eipir(ngs) ! eww(ngs),
      real erw(ngs),esw(ngs),eglw(ngs),eghw(ngs),efw(ngs)
      real ehxw(ngs),ehlw(ngs),egmw(ngs),ehw(ngs) ! eaw(ngs),
      real err(ngs),esr(ngs),eglr(ngs),eghr(ngs),efr(ngs)
      real ehxr(ngs),ehlr(ngs),egmr(ngs) ! ,eipr(ngs),ear(ngs)
      real eri(ngs),esi(ngs),egli(ngs),eghi(ngs),efi(ngs)
      real ehxi(ngs),ehli(ngs),egmi(ngs),ehi(ngs) ! eai(ngs),
      real ers(ngs),ess(ngs),egls(ngs),eghs(ngs),efs(ngs),ehs(ngs)
      real ehxs(ngs),ehls(ngs),egms(ngs),egmip(ngs) ! eas(ngs),

      real exy(ngs,ls:lhab,lc:ls)
      real scxacy(ngs,ls:lhl,lc:ls)
      real cxacy(ngs,ls:lhl,lc:ls)
      real sxxacy(ngs,ls:lhl,li:ls)

      real ehip(ngs),ehlip(ngs),ehlir(ngs)
      real erir(ngs),esir(ngs),eglir(ngs),egmir(ngs),eghir(ngs)
      real efir(ngs),ehir(ngs),eirw(ngs),eirir(ngs),ehr(ngs)
      real erip(ngs),esip(ngs),eglip(ngs),eghip(ngs)
      real efip(ngs),eipi(ngs),eipw(ngs),eipip(ngs)
c
c  arrays for production terms
c
      double precision ptotal(ngs) ! , pqtot(ngs) 
      logical printrates
c
      real pqcwi(ngs),pqcii(ngs),pqrwi(ngs)
      real pqswi(ngs),pqhwi(ngs),pqwvi(ngs)
      real pqgli(ngs),pqghi(ngs),pqfwi(ngs)
      real pqgmi(ngs),pqhli(ngs) ! ,pqhxi(ngs)
      real pqiri(ngs),pqipi(ngs) ! pqwai(ngs),
c
      real pqcwd(ngs),pqcid(ngs),pqrwd(ngs)
      real pqswd(ngs),pqhwd(ngs),pqwvd(ngs)
      real pqgld(ngs),pqghd(ngs),pqfwd(ngs)
      real pqgmd(ngs),pqhld(ngs) ! ,pqhxd(ngs)
      real pqird(ngs),pqipd(ngs) ! pqwad(ngs),
c
c      real pqxii(ngs,nhab),pqxid(ngs,nhab)
c
      real  pctot(ngs)
      real  pcipi(ngs), pcipd(ngs)
      real  pciri(ngs), pcird(ngs)
      real  pccwi(ngs), pccwd(ngs)
      real  pccii(ngs), pccid(ngs)
      real  pcrwi(ngs), pcrwd(ngs)
      real  pcswi(ngs), pcswd(ngs)
      real  pchwi(ngs), pchwd(ngs)
      real  pchli(ngs), pchld(ngs)
      real  pcfwi(ngs), pcfwd(ngs)
      real  pcgli(ngs), pcgld(ngs)
      real  pcgmi(ngs), pcgmd(ngs)
      real  pcghi(ngs), pcghd(ngs)
c
      real  psctot(ngs) 
      real  psccwi(ngs), psccwd(ngs)
      real  psccii(ngs), psccid(ngs)
      real  pscipi(ngs), pscipd(ngs)
      real  psciri(ngs), pscird(ngs)
      real  pscrwi(ngs), pscrwd(ngs)
      real  pscswi(ngs), pscswd(ngs)
      real  pschwi(ngs), pschwd(ngs)
      real  pschli(ngs), pschld(ngs)
      real  pscfwi(ngs), pscfwd(ngs)
      real  pscgli(ngs), pscgld(ngs)
      real  pscgmi(ngs), pscgmd(ngs)
      real  pscghi(ngs), pscghd(ngs)

      real  psccwmi(ngs), psccwmd(ngs)
      real  psccimi(ngs), psccimd(ngs)
      real  pscipmi(ngs), pscipmd(ngs)
      real  pscirmi(ngs), pscirmd(ngs)
      real  pscrwmi(ngs), pscrwmd(ngs)
      real  pscswmi(ngs), pscswmd(ngs)
      real  pschwmi(ngs), pschwmd(ngs)
      real  pschlmi(ngs), pschlmd(ngs)
      real  pscfwmi(ngs), pscfwmd(ngs)
      real  pscglmi(ngs), pscglmd(ngs)
      real  pscgmmi(ngs), pscgmmd(ngs)
      real  pscghmi(ngs), pscghmd(ngs)
c
c  other arrays
c
c
      real wvdf(ngs),tka(ngs) !,akvisc(ngs),ci(ngs),cw(ngs),thdf(ngs)
      real dqisdt(ngs),advisc(ngs) !dqwsdt(ngs), ,schm(ngs),pndl(ngs)

      real qss0(ngs)

      real advisc0,advisc1,tka0

      real qsacip(ngs)
      real pbz(ngs),pres(ngs),presp(ngs),pres0(ngs)
      real pk(ngs)
      real rho0(ngs),dnz(ngs),pi0(ngs),piz(ngs)
      real rhovt(ngs)
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
c     real pceds(ngs) ! ,ppceds(ngs),pmceds(ngs)
      real times
c     real qwfzi(ngs) ! ,qimlw(ngs)
      real ptwfzi(ngs),ptimlw(ngs)
      real psub(ngs),pvap(ngs),pfrz(ngs),ptem(ngs)
c     real fload(ngs)
c     character(LEN=80) filnam
c     character(LEN=15) rrshcm
c     character(LEN=2)  headr1
c     character(LEN=5)  rstime
      character(LEN=6)  rstime
c     character(LEN=2)  nmliter
c
c  iholef = 1 to do hole filling technique version 1
c  which uses all hydrometerors to do hole filling of all hydrometeors
c  iholef = 2 to do hole filling technique version 2
c  which uses an individual hydrometeror species to do hole 
c  filling of a species of a hydrometeor
c
c  iholen = interval that hole filling is done
c
      integer  iholef
      integer  iholen
      parameter (iholef = 0 ) ! = 1)  ERM changed to 0 because now have holefill in solver
      parameter (iholen = 9999999 )!  1)
      real  cqtotn,cqtotn1
      real  cctotn
      real  citotn
      real  crtotn
      real  cstotn
      real  cvtotn
      real  cftotn
      real  cgltotn
      real  cghtotn
      real  chtotn
      real  cqtotp,cqtotp1
      real  cctotp
      real  citotp
      real  ciptotp
      real  crtotp
      real  cstotp
      real  cvtotp
      real  cftotp
      real  chltotp
c      real  chxtotp
      real  cgltotp
      real  cgmtotp
      real  cghtotp
      real  chtotp
      real  cqfac
      real  ccfac
      real  cifac
      real  cipfac
      real  crfac
      real  csfac
      real  cvfac
      real  cffac
      real  cglfac
      real  cghfac
      real  chfac


      real ssifac, qvapor
c
      
      real,allocatable :: xvt(:,:,:,:) ! (nx,nz,2,lc:lhab) ! 1=mass-weighted, 2=number-weighted
      real,allocatable :: tmpn(:,:,:)
c     
c      real prt2d(nx,ny)
c      real prti2d(nx,ny)
c      real prtl2d(nx,ny)
c      real cwrt2d(nx,ny)
c      real cirt2d(nx,ny)
c      real cirrt2d(nx,ny)
c      real ciprt2d(nx,ny)
c      real rwrt2d(nx,ny)
c      real fwrt2d(nx,ny)
c      real swrt2d(nx,ny)
c      real glrt2d(nx,ny)
c      real gmrt2d(nx,ny)
c      real ghrt2d(nx,ny)
c      real hlrt2d(nx,ny)
c      real hwrt2d(nx,ny)
c     
c      real ptot2d(nx,ny)
c      real ptoti2d(nx,ny)
c      real ptotl2d(nx,ny)
c      real cwtot2d(nx,ny)
c      real citot2d(nx,ny)
c      real cirtot2d(nx,ny)
c      real ciptot2d(nx,ny)
c      real rwtot2d(nx,ny)
c      real fwtot2d(nx,ny)
c      real swtot2d(nx,ny)
c      real gltot2d(nx,ny)
c      real gmtot2d(nx,ny)
c      real ghtot2d(nx,ny)
c      real hltot2d(nx,ny)
c      real hwtot2d(nx,ny)
c
c
c   Miscellaneous variables
c
      real :: alphah0
      integer ireadqf,lrho,lqsw,lqgl,lqgm ,lqgh ! ,ltim,ltem,lqcw,lqfw
      integer lqrw,imkgam 
      real vt
      data imkgam /0/
      integer igam        
      real arg
      real, external :: gamma  ! gamma is a function  
      real gaml02
      real erbnd1, fdgt1, costhe1
      real qeps
      real cp608,cv,bta1,cnit,dragh,dnz00,rho00,pii
      real qccrit,gf4br,gf4ds,gf4p5, gf3ds, gf1ds,gr

c      real cnoi,cnoip,cnoir,cnor,cnos,cnogl,cnogm,cnogh,cnof,cnoh
c      real cnohl,
c      real rwdnmx,cwdnmx,cidnmx,xidnmx,swdnmx,gldnmx,gmdnmx
c      real cirdn0, cwdn0, rwdn0, swdn0, gldn0
c      real gmdn0, ghdn0, fwdn0, hwdn0, hldn0
      
      real xdn0(lc:lhab)

c      real ghdnmx,fwdnmx,hwdnmx,hldnmx,rwdnmn,cwdnmn,xidnmn,cidnmn
c      real swdnmn,gldnmn,gmdnmn,ghdnmn,fwdnmn
      integer l,itermax ,ltemq,inumgs, idelq ! , ib
c      real hwdnmn,hldnmn,
      real c1f3,brz,arz,rw,temq ! ,cmn,cmi40,cmi50
c      real ri50,vti50,bsfw,cm50a,a,cm40b,cm50b
      real cnin20,cnin10,cnin1a
      real cnin2a,cnin2b,ssival,tqvcon
      real cd(5)
      real cdx(lc:lhab)
      real cno(lc:lhab), cnox
      real cval,aval,eval,fval,gval ,qsign,ftelwc,qconkq
      real qconm,qconn,cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5,qdiff,argrcnw
      real c4,bradp,bl2,bt2,dtrh,hrifac, hdia0,hdia1,civenta,civentb
      real civentc,civentd,civente,civentf,civentg,cireyn,xcivent
      real cipventa,cipventb,cipventc,cipventd,cipreyn,cirventa
      real cirventb
      integer igmrwa,igmrwb,igmswa, igmswb,igmfwa,igmfwb,igmhwa,igmhwb
      real rwventa ,rwventb,swventa,swventb,fwventa,fwventb,fwventc
      real hwventa,hwventb 
      integer iptotal ! counts number of times that eqtol is exceeded
      integer ilock1,ilock2,ilock3,ilock4,ilockc1,ilockc2
      double precision ptotalmx,ptotalmn
      real psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p
      integer nscmax,nscmin,iwetg
      
      real    hwventc, hlventa, hlventb,  hlventc
      real  glventa, glventb, glventc 
      real   gmventa, gmventb,  gmventc, ghventa, ghventb, ghventc 

      real  dzfacp,  dzfacm,  cmassin,  cwdiar ! , cwmasr
      real  rimmas, rhobar
      real   argtim, argqcw, argqxw, argtem
      real   frcswsw, frcswgl, frcswgm, frcswgh, frcswfw, frcswsw1
      real   frcglgl, frcglgm, frcglgh,  frcglfw, frcglgl1
      real   frcgmgl, frcgmgm, frcgmgh,  frcgmfw, frcgmgm1
      real   frcghgl, frcghgm, frcghgh,  frcghfw,  frcghgh1
      real   frcfwgl, frcfwgm, frcfwgh, frcfwfw,  frcfwfw1
      real   frcswrsw, frcswrgl,  frcswrgm,  frcswrgh, frcswrfw
      real   frcswrsw1
      real   frcrswsw, frcrswgl, frcrswgm, frcrswgh, frcrswfw
      real  frcrswsw1 
      real  frcglrgl, frcglrgm, frcglrgh,  frcglrfw, frcglrgl1
      real  frcrglgl  
      real  frcrglgm,  frcrglgh, frcrglfw, frcrglgl1  
      real  frcgmrgl, frcgmrgm, frcgmrgh, frcgmrfw,  frcgmrgm1
      real  frcrgmgl, frcrgmgm,  frcrgmgh, frcrgmfw, frcrgmgm1
      real  sum,  qweps,  gf2a, gf4a, dqldt, dqidt, dqdt
      real frcghrgl, frcghrgm, frcghrgh, frcghrfw, frcghrgh1, frcrghgl
      real frcrghgm, frcrghgh,  frcrghfw, frcrghgh1
      real    a1,a2,a3,a4,a5,a6
      real   gamss
      real cdw, cdi, denom1, denom2, delqci1, delqip1 ! , dtz1, dtz2
      real cirtotn,  ciptotn, cgmtotn, chltotn,  cirtotp
      real  cgmfac, chlfac,  cirfac
      integer igmhla, igmhlb, igmgla, igmglb, igmgma,  igmgmb
      integer igmgha, igmghb
      integer idqis, item, itim0 ! ,  itim
      integer  iqgl, iqgm, iqgh, iqrw, iqsw ! ,iqcw, iqfw
      integer  itertd, ia

      save iread
      
#ifndef MPI
C$    volatile iptotal,ptotalmx,ptotalmn,nscmax,nscmin,iwetg,
C$   &   psctotmx,psctotmn,psctot1,psctot2,sctot1n,sctot1p
#endif
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

#ifdef MPI
      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag

      logical :: debug_mpi = .false.
      integer       :: nampi
#endif

c
c ####################################################################
c
c  Start routine
c
c ####################################################################
c

c
c  slope intercepts
c
      IF ( ngs .lt. nz ) THEN
       write(0,*) 'Error in ICE10: Must have ngs .ge. nz!'
       STOP
      ENDIF

#ifdef MPI
      write(0,*) 'ICE10: Not converted to MPI yet, stopping...'
      STOP
#endif

      allocate ( db1(nx,nz) )
      allocate ( db0(nx,nz) )
      allocate ( xvt(nx,nz,2,lc:lhab) )
      allocate ( tmpn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor) )

      lrt(lc) = lrtc
      lrt(lr) = 0 !lrtr
      lrt(li) = lrti
      lrt(lip) = lrtip
      lrt(lir) = lrtir
      lrt(ls) = lrts
      lrt(lgl) = lrtgl
      lrt(lgm) = lrtgm
      lrt(lgh) = lrtgh
      lrt(lf) = lrtf
      lrt(lh) = 0 !lrth
      lrt(lhl) = 0 !lrthl
      
      IF ( lrtc > 1 ) iptime = 1
      IF ( lrtgm > 1 ) iptime = 2
      
      lsc(lc) = lscw
      lsc(lr) = lscr
      lsc(li) = lsci
      lsc(lip) = lscip
      lsc(lir) = lscir
      lsc(ls) = lscs
      lsc(lgl) = lscgl
      lsc(lgm) = lscgm
      lsc(lgh) = lscgh
      lsc(lf) = lscf
      lsc(lh) = lsch
      lsc(lhl) = lschl
      

      ln(lc) = lnc
      ln(lr) = lnr
      ln(li) = lni
      ln(lip) = lnip
      ln(lir) = lnir
      ln(ls) = lns
      ln(lgl) = lngl
      ln(lgm) = lngm
      ln(lgh) = lngh
      ln(lf) = lnf
      ln(lh) = lnh
      ln(lhl) = lnhl

      ipc(lc) = 2
      ipc(lr) = 3
      ipc(li) = 1
      ipc(lip) = 1
      ipc(lir) = 1
      ipc(ls) = 4
      ipc(lgl) = 5
      ipc(lgm) = 5
      ipc(lgh) = 5
      ipc(lf) = 5
      ipc(lh) = 5
      ipc(lhl) = 5
      
      alphah0 = 0.0
       rz =  ((4. + alphah0)*(5. + alphah0)*(6. + alphah0)*(1. + rnu))/ 
     &  ((1. + alphah0)*(2. + alphah0)*(3. + alphah0)*(2. + rnu))

      
C      cno(lc)  = 1.0e+08
C      cno(li)  = 1.0e+08
C      cno(lip) = 1.0e+08 
C      cno(lir) = 1.0e+08 
C      cno(lr)  = 8.0e+06 
C      cno(ls)  = 8.0e+06 
C      cno(lgl) = 4.0e+05 
C      cno(lgm) = 4.0e+05 
C      cno(lgh) = 4.0e+05 
C      cno(lf)  = 4.0e+05
C      cno(lh)  = 4.0e+04
C      cno(lhl) = 1.0e+03

      CALL setcno(cno)
c
c  density maximums and minimums
c
      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(lip) =  917.0
      xdnmx(lir) =  400.0
      xdnmx(ls) =  200.0
      xdnmx(lgl) =  400.0
      xdnmx(lgm) =  600.0
      xdnmx(lgh) =  800.0
      xdnmx(lf) =  900.0
      xdnmx(lh) =  900.0
      xdnmx(lhl) =  900.0
c
      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(lip) =  100.0
      xdnmn(lir) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lgl) =  200.0
      xdnmn(lgm) =  400.0
      xdnmn(lgh) =  600.0
      xdnmn(lf) =  700.0
      xdnmn(lh) =  700.0
      xdnmn(lhl) =  900.0

      xdn0(lir) = 300.0
      xdn0(lc) = 1000.0
      xdn0(lr) = rho_qr_i10
      xdn0(ls) = rho_qs_i10
      xdn0(lgl) = rho_qgl_i10 ! (0.5)*(xdnmn(lgl)+xdnmx(lgl))
      xdn0(lgm) = rho_qgm_i10 ! (0.5)*(xdnmn(lgm)+xdnmx(lgm))
      xdn0(lgh) = rho_qgh_i10 ! (0.5)*(xdnmn(lgh)+xdnmx(lgh))
      xdn0(lf) =  rho_qf_i10  ! (0.5)*(xdnmn(lf)+xdnmx(lf))
      xdn0(lh) =  rho_qh_i10  ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      xdn0(lhl) = rho_qhl_i10 ! (0.5)*(xdnmn(lhl)+xdnmx(lhl))

c
c  Set terminal velocities...
c    also set drag coefficients
c
      cdx(lr) = 0.60
      cdx(lh) = 0.45
      cdx(lhl) = 0.45
      cdx(lf) = 0.45
      cdx(lgh) = 0.60
      cdx(lgm) = 0.80
      cdx(lgl) = 0.80
      cdx(lir) = 2.00
      cdx(ls) = 2.00
      cd(1) = cdx(ls)
      cd(2) = cdx(lgl)
      cd(3) = cdx(lgm)
      cd(4) = cdx(lgh)
      cd(5) = cdx(lf)

      
      dtp = ntmul*dtp1
      write(iunit,*) 'IN TDMPHY'
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix)

       call sett09s(nx,ny,nz,nor,nor,na,lt,
     :             t77,pinit,p2,an,pn,pb,
     :          t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t00 )

c      DO kz=-nor+ng1,nz+nor
c      DO jy=-nor+ng1,ny+nor
c      DO ix=-nor+ng1,nx+nor
c        t77(ix,jy,kz) = ((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
c        t0(ix,jy,kz) = ac(ix,jy,kz,lt)*t77(ix,jy,kz)
c        t1(ix,jy,kz) = 0.0
c        t2(ix,jy,kz) = 0.0
c        t3(ix,jy,kz) = 0.0
c        t4(ix,jy,kz) = 0.0
c        t5(ix,jy,kz) = 0.0
c        t6(ix,jy,kz) = 0.0
c        t7(ix,jy,kz) = 0.0
c        t8(ix,jy,kz) = 0.0
c        t9(ix,jy,kz) = 0.0
c
c      END DO
c      END DO
c      END DO


c
c  open input data file for the microphysics code
c
      if ( iread .eq. 0 ) then
c
      iread = 1
c

      write(iunit,*) 'INPUT:  ndebug (-1=no debug)'
      write(iunit,*) 'ndebug= ',ndebug

      ido(lc) = idocw ; ido(lr) = idorw ; ido(li) = idoci ; ido(lir) = idoir
      ido(lip) = idoip; ido(ls) = idosw

c
      ido(lgl) = idogl ; ido(lgm) = idogm ; ido(lgh) = idogh 
      ido(lf)  = idofw ; ido(lh)  = idohw ; ido(lhl) = idohl

      write(iunit,'(a,6(1x,i2))') 'ido(lc:ls)', ! idocw, idorw, idoci, idoir, idoip, idosw
     : ido(lc), ido(lr), ido(li), ido(lir), ido(lip), ido(ls)

c
      write(iunit,'(a,6(1x,i2))') 'ido(lgl:lhl)', ! idogl, idogm, idogh, idofw, idohw, idohl
     : ido(lgl), ido(lgm), ido(lgh), ido(lf),ido(lh),ido(lhl)

      iexy(:,:) = 0
      
      

      iexy(ls,li) = ieswi ; iexy(ls,lir) = ieswir ; 
      iexy(ls,lip) = ieswip ; iexy(ls,lc) = ieswc ; iexy(ls,lr) = ieswr ;

      iexy(lgl,ls)  = ieglsw ; iexy(lgl,li) = iegli ; iexy(lgl,lir) = ieglir ; 
      iexy(lgl,lip) = ieglip ; iexy(lgl,lc) = ieglc ; iexy(lgl,lr)  = ieglr ;

      iexy(lgm,ls)  = iegmsw ; iexy(lgm,li) = iegmi ; iexy(lgm,lir) = iegmir ; 
      iexy(lgm,lip) = iegmip ; iexy(lgm,lc) = iegmc ; iexy(lgm,lr)  = iegmr ;

      iexy(lgh,ls)  = ieghsw ; iexy(lgh,li) = ieghi ; iexy(lgh,lir) = ieghir ; 
      iexy(lgh,lip) = ieghip ; iexy(lgh,lc) = ieghc ; iexy(lgh,lr)  = ieghr ;

      iexy(lf,ls)  = iefwsw ; iexy(lf,li) = iefwi ; iexy(lf,lir) = iefwir ; 
      iexy(lf,lip) = iefwip ; iexy(lf,lc) = iefwc ; iexy(lf,lr)  = iefwr ;

      iexy(lh,ls)  = iehwsw ; iexy(lh,li) = iehwi ; iexy(lh,lir) = iehwir ; 
      iexy(lh,lip) = iehwip ; iexy(lh,lc) = iehwc ; iexy(lh,lr)  = iehwr ;

      iexy(lhl,ls)  = iehlsw ; iexy(lhl,li) = iehli ; iexy(lhl,lir) = iehlir ; 
      iexy(lhl,lip) = iehlip ; iexy(lhl,lc) = iehlc ; iexy(lhl,lr)  = iehlr ;

        write(iunit,*) 'nonigrd= ',nonigrd, '  rgard= ',rgard,
     :   ' ewfac= ',ewfac, ' eii0 = ',eii0, ' eii1 = ',eii1


      ffrz = fconv
      
      write(iunit,*) 'ircnw, qminrncw=', ircnw, qminrncw

      write(iunit,*) 'cwdiap=', cwdiap

      write(iunit,*) 'cwdisp= ', cwdisp

      write(iunit,*) 'iauttim= ',iauttim

      write(iunit,*) 'auttim= ',auttim

      write(iunit,*) 'qcwmntim= ', qcwmntim

      cwccn = ccn
      alfarim = Min (0.0,alfarim)
      write(iunit,*) 'alfarim = ', alfarim

         pi = 4.0*atan(1.0)
         xvcmx = (4./3.)*pi*xcradmx**3

         write(iunit,*) 'iccwflg, irenuc = ', iccwflg, irenuc
         
         write(iunit,*) 'xcradmx, xvcmx = ', xcradmx, xvcmx
         
         write(iunit,*) 'ciintmx = ',ciintmx
         
      write(iunit,*) 'Intercepts: r,s,gl,gm,gh,lf,lh,lhl ',cno(lr),cno(ls),cno(lgl),
     :              cno(lgm),cno(lgh),cno(lf),cno(lh),cno(lhl)
      
      cck = 0.6
      cckm = cck-1.
c      ccn = cwccn
      ccne = 0.9893*1.e6*(1.e-6*Abs(cwccn))**0.77
      write(iunit,*) 'cwccn, cck, ccne = ',cwccn,cck,ccne


      IF ( cwccn .lt. 0.0 ) THEN
      cwccn = Abs(cwccn)
      ccwmx = cwccn
      ELSE
      ccwmx = cwccn*1.4
      ENDIF

      write(iunit,'(a,i2,2x,i2,1x,1pe12.5,3i3)') 
     :  'itype1, itype2,cimas0, icfn, ihrn, ibfc = ',
     :    itype1, itype2, cimas0, icfn, ihrn, ibfc

      write(iunit,*) 'ipelec= ',ipelec, ', icorona = ',icorona

      write(iunit,*) 'ipotslv= ',ipotslv

      write(iunit,*) 'ilghtest= ',ilghtest
c
      write(iunit,*) 'ilight= ',ilight
c
      write(iunit,*) 'nhelec= ',nhelec
c
      write(iunit,*) 'nliter= ',nliter
c
      write(iunit,*) 'iscreen= ',iscreen
c
      write(iunit,*) 'eqtot= ',eqtot

      IF ( ipelec .ge. 1 ) THEN
      
      write(iunit,*) 'nonigrd = ', nonigrd
      
      IF ( nonigrd .eq. 0 ) THEN
        isaund = Nint(rgard)
        write (iunit,*) 'RAR: isaund = ',isaund
      ELSEIF ( nonigrd .eq. 1 ) THEN
         isaund = Nint(rgard)
         IF ( isaund .ge. 1 ) THEN
          write(iunit,*) 'using qxacw-based RAR calculation'
         ELSE
         nonigrd = 0
         isaund = 0
          write (iunit,*) 'Saunders1991: isaund = ',isaund
         ENDIF
      ELSEIF ( nonigrd .eq. -1 ) THEN
        isaund = Nint(rgard)
        write (iunit,*) 'Takahashi: isaund = ',isaund
      ELSEIF ( nonigrd .eq. 2 ) THEN
        rgard1 = rgard
      ENDIF

      ENDIF

! end of cut for SWM -- input vars moved to a module/namelist

      dtrim = rimtim
      argtim = ((rimtim-0.5/xftim)*xftim + yftim)
      itim0 = min(max(ifix(argtim),1), mtim)
      write(iunit,*) 'itim0 = ',itim0


      IF ( nonigrd .eq. -1 ) THEN
       open(unit=90,file='takahashi.txt',form='formatted')
       DO ix=nlwc,0,-1
        read(90,*) (takalu(jy,ix), jy=0,31)
       END DO
       close(90)
      END IF
      
      ireadqf = 0

      IF ( ipconc .lt. 3 .and. .not. lraincv ) THEN

      if ( iptime .ge. 1 ) then
      if ( ireadqf .eq. 0 ) then
      ireadqf = 1
      write(iunit,*) '*****READ DATA FOR LOOK UP TABLES******'


c  r/x fractions
c
c
c

      
      close(90)
      close(92)
      if ( ireadmic .eq. 0 ) then
      open(unit=90,file='sam.qfswarw.data',form='unformatted')
      open(unit=92,file='sam.qfrwasw.data',form='unformatted')
      read(90) qfswarw
      read(92) qfrwasw
      end if
c
      if ( ireadmic .eq. 1 ) then
      open(unit=90,file='sam.qfswarw.data',form='formatted')
      open(unit=92,file='sam.qfrwasw.data',form='formatted')
      do lrho = 1,mrho
      do lqrw = 1,mqcw
      do lqsw = 1,mqxw
      read(90,*) qfswarw(lqrw,lqsw,lrho)
      read(92,*) qfrwasw(lqrw,lqsw,lrho)
      end do
      end do
      end do
      end if
      close(unit=90)
      close(unit=92)
c
c
c
      if ( ireadmic .eq. 0 ) then
      open(unit=90,file='sam.qfglarw.data',form='unformatted')
      open(unit=92,file='sam.qfrwagl.data',form='unformatted')
      read(90) qfglarw
      read(92) qfrwagl
      end if
c
      if ( ireadmic .eq. 1 ) then
      open(unit=90,file='sam.qfglarw.data',form='formatted')
      open(unit=92,file='sam.qfrwagl.data',form='formatted')
      do lrho = 1,mrho
      do lqrw = 1,mqcw
      do lqgl = 1,mqxw
      read(90,*) qfglarw(lqrw,lqgl,lrho)
      read(92,*) qfrwagl(lqrw,lqgl,lrho)
      end do
      end do
      end do
      end if
      close(unit=90)
      close(unit=92)
c
c
c
      if ( ireadmic .eq. 0 ) then
      open(unit=90,file='sam.qfgmarw.data',form='unformatted')
      open(unit=92,file='sam.qfrwagm.data',form='unformatted')
      read(90) qfgmarw
      read(92) qfrwagm
      end if
c
      if ( ireadmic .eq. 1 ) then
      open(unit=90,file='sam.qfgmarw.data',form='formatted')
      open(unit=92,file='sam.qfrwagm.data',form='formatted')
      do lrho = 1,mrho
      do lqrw = 1,mqcw
      do lqgm = 1,mqxw
      read(90,*) qfgmarw(lqrw,lqgm,lrho)
      read(92,*) qfrwagm(lqrw,lqgm,lrho)
      end do
      end do
      end do
      end if
      close(unit=90)
      close(unit=92)
c
c
c
      if ( ireadmic .eq. 0 ) then
      open(unit=90,file='sam.qfgharw.data',form='unformatted')
      open(unit=92,file='sam.qfrwagh.data',form='unformatted')
      read(90) qfgharw
      read(92) qfrwagh
      end if
c
      if ( ireadmic .eq. 1 ) then
      open(unit=90,file='sam.qfgharw.data',form='formatted')
      open(unit=92,file='sam.qfrwagh.data',form='formatted')
      do lrho = 1,mrho
      do lqrw = 1,mqcw
      do lqgh = 1,mqxw
      read(90,*) qfgharw(lqrw,lqgh,lrho)
      read(92,*) qfrwagh(lqrw,lqgh,lrho)
      end do
      end do
      end do
      end if
      close(unit=90)
      close(unit=92)
c
      end if
      end if
      
      ENDIF ! ( ipconc .lt. 5 )
c

      vfrz = 0.523599*(dfrz)**3 
      vmlt = 0.523599*(dmlt)**3 

      do l = 1,nqsat
      temq = 163.15 + (l-1)*fqsat
      tabqvs(l) = exp(caw*(temq-273.16)/(temq-cbw))
      dtabqvs(l) = ((-caw*(-273.16 + temq))/(temq - cbw)**2 +
     :                 caw/(temq - cbw))*tabqvs(l)
      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
      end do
c
      end if    ! ( iread .eq. 0 ) 
c
c

c  Gamma function
c
c
c     compute fgoi(arg*100)
c
      if ( imkgam .eq. 0 ) then
      imkgam = 1
      call makegmoi()
!      do igam = 1,1000
!      arg = 0.01*igam
!      gmoi(igam) = gamma(arg)
c      write(97,*) igam,gmoi(igam)
!      end do
      end if
c
c
c  electricity constants
c
c  mixing ratio epsilon
c
      qeps  = 1.0e-20

c
c  rebound efficiency (erbnd)
c
c

c
c      times = (nstep-1)*(dtp1*dtfac)
      times = (nstep)*(dtp1*dtfac)
c      headr1 = 'pa'
c
c  constants
c
      cp608 = 0.608
      cv = 717.0
      ar = 841.99666  
      br = 0.8
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
      rho00 = 1.225
c      cs = 4.83607122
c      ds = 0.25
c  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pi = 4.0*atan(1.0)
      pii = 1./pi
      pid4 = pi/4.0 
!      qscrit = 6.0e-04
      gf1 = gamma(1.0)
      gf1p5 = 0.8862269255
      gf2 = gamma(2.0)
      gf3 = gamma(3.0)
      gf3p5 = 3.32335097
      gf4 = gamma(4.0)
      gf5 = gamma(5.0)
      gf6 = gamma(6.0)
      gf7 = gamma(7.0)
      gf4br = gamma(4.0+br)
      gf4ds = gamma(4.0+ds)
      gf4p5 = gamma(4.0+0.5)
      gf3ds = gamma(3.0+ds)
      gf1ds = gamma(1.0+ds)
      gr = 9.8
c
c  constants
c
      c1f3 = 1.0/3.0
c
c  general constants for microphysics
c
      brz = 100.0
      arz = 0.66

      call setqxmin(qxmin)
      
c      qxmin(:) = 1.0e-12
      
c      qxmin(lc) = 1.0e-09
c      qxmin(lr) = 1.0e-07
c      qxmin(li) = 1.0e-09
c      qxmin(lip) = 1.0e-09
c      qxmin(lir) = 1.0e-08
c      qxmin(ls) = 1.0e-9
c      qxmin(lgl) = 1.0e-07
c      qxmin(lgm) = 1.0e-07
c      qxmin(lgh) = 1.0e-07
c      qxmin(lf) = 1.0e-07
c      qxmin(lh) = 1.0e-07
c      qxmin(lhl) = 1.0e-07
      
      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn
      xvmn(lgl) = xvgmn
      xvmn(lgm) = xvgmn
      xvmn(lgh) = xvgmn
      xvmn(lf) = xvfmn
      xvmn(lh) = xvhmn
      xvmn(lhl) = xvhlmn

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx
      xvmx(ls) = xvsmx
      xvmx(lgl) = xvgmx
      xvmx(lgm) = xvgmx
      xvmx(lgh) = xvgmx
      xvmx(lf) = xvfmx
      xvmx(lh) = xvhmx
      xvmx(lhl) = xvhlmx


      tdtol = 1.0e-05
      thnuc = 233.15
      rw = 461.5
      advisc0 = 1.832e-05
      advisc1 = 1.718e-05
      tka0 = 2.43e-02
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
c
c
      if ( imake .eq. 0 ) then
      imake = 1


      ventr   = Gamma(rnu + 4./3.)/(rnu + 1.)**(1./3.)/Gamma(rnu + 1.)
      ventc   = Gamma(cnu + 4./3.)/(cnu + 1.)**(1./3.)/Gamma(cnu + 1.)
c
c  Set up look up tables for supersaturation w.r.t. liq and ice
c
CVD$L SKIP
c      do l = 1,nqsat
c      temq = 163.15 + (l-1)*fqsat
c      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
c      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
c      end do
            
c
      itermax = 0
      DO mgs = 1,ngs
        il0(mgs) = 1
      ENDDO

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag
      
      do kz = kzb,kze
#else
      DO kz = 1,nz-kstag
#endif
      piz(kz) = pinit(kz) ! (pb(kz)/poo)**cap ! pinit(kz)+p2(igs(mgs),jy,kgs(mgs)) 
      dnz(kz) = db(kz)
      pbz(kz) = pb(kz)
      temp(kz) = piz(kz)*ab(kz,lt)

      tempc(kz) = temp(kz) - tfr
      tempbz(kz) = temp(kz) 
      elv(kz) = 2500300. - 2369.276*tempc(kz) 
      elf(kz) = 335717. - 2369.276*tempc(kz) 
      els(kz) = elv(kz) + elf(kz)
      elvs(kz) = elv(kz)*elv(kz)
      elss(kz) = els(kz)*els(kz)
      gamw(kz) = elv(kz)*cpi/piz(kz)
      gams(kz) = els(kz)*cpi/piz(kz)
      cqv1(kz) = 4098.0258*piz(kz)*gamw(kz)
      cqv2(kz) = 5807.6953*piz(kz)*gams(kz)
      cc3(kz) = cpi*elf(kz)/piz(kz)
      
      ENDDO

c
c  cw constants in mks units
c
c      cwmasn = 4.25e-15  ! radius of 1.0e-6
      cwmasn = 5.23e-13  ! radius of 5.0e-6
      cwmasn5 =  5.23e-13
      cwradn = 5.0e-6
      cwmasx = 5.25e-10  ! radius of 50.0e-6
      mwfac = 6.0**(1./3.)
      IF ( ipconc .ge. 2 ) THEN
        cwmasn = xvmn(lc)*1000.
        cwradn = 1.0e-6
        cwmasx = xvmx(lc)*1000.
      ENDIF
        rwmasn = xvmn(lr)*1000.
        rwmasx = xvmx(lr)*1000.

#ifdef MPI
      kzb=1
      kze=ktile
      if (kzend.eq.nzend) kze=kzend-kzbeg+1-kstag
      
      do kz = kzb,kze
#else
      do kz = 1,nz-kstag 
#endif
      cwnccn(kz) = cwccn*dnz(kz)/dnz(1)  
 
      ENDDO
      
      cwc1 = 6.0/(pi*1000.)
      cwc0 = pii ! 6.0*pii
c
c  ci constants in mks units
c
      cimasn = 6.88e-13 
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter
      do 40 kz = 1,nz-kstag 
c      cidn(kz) = 900.0
      cimcnd(kz) = 1.0e-12 
      cincmn(kz) = 1.0e+01
      cincmx(kz) = 1.0e+10
c      cinccn(kz) = 1.0e+09
  40  continue

c 
c  constants for paramerization
c
      if ( ncuse .eq. 0 ) then
      do 50 kz = 1,nz-kstag
      wvdf(kz) = (2.11e-05)*((temp(kz)/tfr)**1.94)*
     >  (101325.0/(pbz(kz)))
      advisc(kz) = advisc0*(416.16/(temp(kz)+120.0))*
     >  (temp(kz)/296.0)**(1.5)
      tka(kz) = tka0*advisc(kz)/advisc1
  50  continue
      end if
c
      if (ndebug .gt. 0 ) print*,'dbg = 1'
c
c
c polar conductivity interpolation
c   (moved to driver subroutine)
c         
c  Set initial temperatures for hydrometeors
c

c
c  end if for constants (imake)
c
      end if
c
c  zero the precip flux arrays (2d)
c
c      if ( nstep .eq. nstart ) then
c      do jy = 1,ny-jstag
c      do ix = 1,nx-istag
c      ciptot2d(ix,jy) = 0.0
c      cirtot2d(ix,jy) = 0.0
c      cwtot2d(ix,jy) = 0.0
c      rwtot2d(ix,jy) = 0.0
c      swtot2d(ix,jy) = 0.0
c      fwtot2d(ix,jy) = 0.0
c      gltot2d(ix,jy) = 0.0
c      gmtot2d(ix,jy) = 0.0
c      ghtot2d(ix,jy) = 0.0
c      hwtot2d(ix,jy) = 0.0
c      hltot2d(ix,jy) = 0.0
c      end do
c      end do
c      if ( ihabdo .lt. 1 ) then
c      do jy = 1,ny-jstag
c      do ix = 1,nx-istag
c      citot2d(ix,jy) = 0.0
c      end do
c      end do
c      end if     
c      end if     


c
c  set 3-d temperature and saturation mixing ratio variables 
c
      ierr = 0
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ltemq)
!! c$omp PARALLEL DO  IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix,ltemq)
      do kz = 1,nz-1
      do jy = 1,ny-jstag
      do ix = 1,nx-istag

c
c  temperature
c
c      t0(ix,jy,kz) = an(ix,jy,kz,lt)
c     >  *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
      
      
c      t0(ix,jy,kz) = 
c     :     theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)))
c
c  look-up index
c 
      ltemq = Int( (t0(ix,jy,kz)-163.15)/fqsat+1.5 )
      IF ( ltemq .ge. nqsat .or. ltemq .le. 0 ) THEN
        ierr = ierr + 1
        write(iunit,*) 'WARNING 1: ltemq .ge. nqsat! ierr = ',ierr
        write(iunit,*) 'ix,jy,kz,t0,an,ltemq,t77 = ',ix,jy,kz,t0(ix,jy,kz),
     :    an(ix,jy,kz,lt), ltemq,t77(ix,jy,kz)
        write(0,*) 'WARNING 1: ltemq .ge. nqsat! ierr = ',ierr
        write(0,*) 'ix,jy,kz,t0,an,ltemq = ',ix,jy,kz,t0(ix,jy,kz),
     :    an(ix,jy,kz,lt), ltemq
        ltemq = Max(1, Min (nqsat,ltemq))
       
      ENDIF
      IF ( ierr .gt. 4 ) THEN
       write(0,*) 'too many warnings! STOP'
       write(iunit,*) 'too many warnings! STOP'
       STOP
      ENDIF
c
c saturation mixing ratio for liquid
c
      a1 = t00(ix,jy,kz)

      t8(ix,jy,kz) = a1*tabqvs(ltemq)

c
c saturation mixing ratio for ice
c
      t9(ix,jy,kz) = a1*tabqis(ltemq)
c
      end do
      end do
      end do
      
c      IF ( ierr .gt. 4 ) THEN
c        call flush(iunit)
c        STOP
c      ENDIF
c
c
c  Stuff for Meyers et al. (1992) and Ferrier (1994) ice nucleation
c
c  constants
c
      cnin20 = 1.0e3
      cnin10 = 5.0e1
      cnin1a = 4.5
      cnin2a = 12.96
      cnin2b = 0.639
c
c  calculate rate of nucleation
c


! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix,ssival)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix,ssival)
      do kz = 1,nz-1
      do jy = 1,ny-jstag
      do ix = 1,nx-istag
c
c filter supersaturation field (if ifilt = 1)
c
c
      ssival = Min(t8(ix,jy,kz),max(an(ix,jy,kz,lv),0.0))/t9(ix,jy,kz)
c
c      t7(ix,jy,kz) = 0.0
c
      if ( ssival .gt. 1.0 ) then
c
      if ( t0(ix,jy,kz).le.268.15 ) then
c     if ( t0(ix,jy,kz).le.273.15 ) then
        
       dp1 = cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,jy,kz) = Min(dp1, 1.0d30)
      end if
      
c
c   have turned off nucleation by Meyer at warm temperatures
c
      IF ( imeyers5 ) THEN
      if ( t0(ix,jy,kz).lt.273.15.and.t0(ix,jy,kz).gt.268.15 ) then
      qvapor = max(an(ix,jy,kz,lv),0.0) 
      ssifac = 0.0
      if ( (qvapor-t9(ix,jy,kz)) .gt. 1.0e-5 ) then
      if ( (t8(ix,jy,kz)-t9(ix,jy,kz)) .gt. 1.0e-5 ) then
      ssifac = (qvapor-t9(ix,jy,kz))
     >             /(t8(ix,jy,kz)-t9(ix,jy,kz))
      ssifac = ssifac**cnin1a   
      end if
      end if
      t7(ix,jy,kz) = cnin10*ssifac
     >  *exp(-(t0(ix,jy,kz)-tfr)*bta1)
      end if
      ENDIF
c
      end if
c
      end do
      end do
      end do
c
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO IF (nx*ny*nz .gt. 40*40*42), 
!! c$omp+ DEFAULT(SHARED),PRIVATE(kz,jy,ix)
      DO kz=-nor+ng1,nz+nor
      DO jy=-nor+ng1,ny+nor
      DO ix=-nor+ng1,nx+nor
        t8(ix,jy,kz) = 0.0
        t9(ix,jy,kz) = 0.0
      END DO
      END DO
      END DO
c
c  set save counter (number of saves):  nsvcnt
c
      nsvcnt = 0
      iend = 0
      iptotal = 0
      nscmax = 0
      nscmin = 0
      ptotalmx = 0.0
      ptotalmn = 0.0
      psctotmx = 0.0
      psctotmn = 0.0     
      psctot1 = 0.0
      psctot2 = 0.0 
      sctot1n = 0.0
      sctot1p = 0.0
      ilock1 = 0
      ilock2 = 0
      ilock3 = 0
      iwetg = 0
      ilock4 = 0
      ilockc1 = 1
      ilockc2 = 2
     

c      timetd1 = etime(tarray)
c      timetd1 = tarray(1)

c 
#ifndef MPI
C$     ndebug = -1
c cmic$  cncall
#endif
c***********************************************************
c  start jy loop
c***********************************************************
c
      db0(:,:) = 1.0
      
      do 9999 jy = 1,ny-jstag
c
c  VERY IMPORTANT:  SET jgs = jy
c
      jgs = jy
      
c
c  Do fallout stuff now if itfall=0 (boxfall) or 2 (crowfall)
c
      IF ( itfall .eq. 0 .or. itfall .eq. 1 ) THEN
      

c
c  zero the precip flux arrays (2d)
c

      xvt(:,:,:,:) = 0.0
      
      do kz = 1,nz-1
      DO ix = 1,nx-1
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo


c set up vt arrays:
c      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
c     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)

      call ddfallnx(nx,ny,nz,nor,na,dtp,dz,jgs,ngs,
     :  xvt,
     :  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx,
     :  rwmasn,rwmasx,cwradn,
     :  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn,
     :  itype1,itype2,infdo)
c     :  rwmasn,rwmasx)
      
c charge

      IF ( ipelec .ge. 2 .and. iscfall .ge. 1 ) THEN
      
      DO il = lc,lhab
      
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :    an,db0,imapz,mzdist,lsc(il),0,xfall)

      ENDDO
            
      ENDIF ! ipelec
c mixing ratio      
      DO il = lc,lhab

      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,il,1,xfall)
       
      ENDDO
      
c time
       IF ( iptime .ge. 1 ) THEN     

      DO il = lc,lhab
      IF ( lrt(il) .gt. 1 ) THEN
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,lrt(il),1,xfall)
      ENDIF
      ENDDO

     
      ENDIF


      
      IF ( ipconc .gt. 0 ) THEN
      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN

      IF ( infall .eq. 2 ) THEN
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!
          DO kz = 1,nz
            DO ix = 1,nx
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
      ENDIF

         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)


         IF ( infall .eq. 2 .and.
     :       ( il .eq. lr .or. (il .ge. ls .and. il .le. lhl) )) THEN
c     :        .or. il .eq. lhl )) THEN

           
           call boxfall(nx,ny,nz,nor,1,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :       tmpn,db0,imapz,mzdist,1,0,xfall0)

           DO kz=1,nz-1
           DO ix=1,nx-1
             
             an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
             
           ENDDO
           ENDDO
           
c          ENDDO

         ENDIF
        ENDIF
      ENDDO
      
      ENDIF


      
      ELSEIF ( itfall .eq. 2 ) THEN
      
c
c  zero the precip flux arrays (2d)
c
      xvt(:,:,:,:) = 0.0
      
      
      do kz = 1,nz-1
      DO ix = 1,nx-1
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

      call ddfallnx(nx,ny,nz,nor,na,dtp,dz,jgs,ngs,
     :  xvt,
     :  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx,
     :  rwmasn,rwmasx,cwradn,
     :  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn,
     :  itype1,itype2,infdo)
c     :  rwmasn,rwmasx)


      vtmax = 0.0
      
      do kz = 1,nz-1
      do ix = 1,nx-1
c
      
! fix gt for SWM
      DO il = lr,lhab
      vtmax = Max(vtmax,xvt(ix,kz,1,il)*gz(kz))
      ENDDO

      ENDDO
      ENDDO


      
      IF ( dtp*vtmax .lt. 0.9 ) THEN
        ndfall = 1
      ELSE
        ndfall = Max(2, Int(dtp*vtmax) + 1)
      ENDIF
      
      IF ( ndfall .gt. 1 ) THEN
        dtptmp = dtp
        dtp = dtp/Real(ndfall)
      ENDIF

      DO n = 1,ndfall

      IF ( n .ge. 2 ) THEN
c
c  zero the precip flux arrays (2d)
c
      xvt(:,:,:,:) = 0.0
      
      
      do kz = 1,nz-1
      DO ix = 1,nx-1
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

      call ddfallnx(nx,ny,nz,nor,na,dtp,dz,jgs,ngs,
     :  xvt,
     :  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx,
     :  rwmasn,rwmasx,cwradn,
     :  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn,
     :  itype1,itype2,infdo)
c     :  rwmasn,rwmasx)
      
      ENDIF ! (n .ge. 2)
      
c mixing ratio      

      DO il = lc,lhab
       call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,il,1,xfall)
      ENDDO
      
c time
       IF ( iptime .ge. 1 ) THEN     

      DO il = lc,lhab
      IF ( lrt(il) .gt. 1 ) THEN
      IF ( .not. ( il .eq. lr .or. il .eq. lh .or. il .eq. lhl) ) THEN
      call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,lrt(il),1,xfall)
      ENDIF
      ENDIF
      ENDDO

     
      ENDIF


c charge
      IF ( ipelec .ge. 2 .and. iscfall .ge. 1 ) THEN
      
c      iscfall = 1
      tmp2b = 0.0
      tmp2e = 0.0

      DO il = lc,lhab

      tmp = 0
       DO k = 1,nz-1
       DO i = 1,nx-1
         dv = dxx(i)*dyy(jy)*dzz(k)
         tmp = tmp + an(i,jy,k,lsc(il))*dv
       ENDDO
       ENDDO
       tmp1b = tmp1b + tmp
       tmp2b = tmp2b + tmp
       
!       print*,'Start total for lsc = ',lsc(il),' is ',tmp


        call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,iscfall,il),
     :    an,db0,imapz,mzdist,lsc(il),0,xfall)
      
      
      tmp2 = 0
       DO k = 1,nz-1
       DO i = 1,nx-1
         tmp2 = tmp2 + an(i,jy,k,lsc(il))*dxx(i)*dyy(jy)*dzz(k)
       ENDDO
       ENDDO
       IF ( tmp .ne. 0.0 ) THEN
       IF ( Abs(tmp2 - tmp)/tmp .gt. 0.0001 ) THEN
!       print*, 'Fallout for lsc(il),il = ',lsc(il),il,'  jy = ',jy
!       print*,'
!       print*,'Begin/End total for lsc = ',lsc(il),il,' is ',tmp,tmp2
       ENDIF
       ENDIF
       tmp1e = tmp1e + tmp2
       tmp2e = tmp2e + tmp2
      
      ENDDO
      print*, 'Charge: jy,tmp2b,tmp2e = ',jy,tmp2b,tmp2e
      
      
      ENDIF ! ipelec
      
c      infall = 2
c      DO il = lc,lhab
c        IF ( ipconc .ge. ipc(il) ) THEN
c         call crowfall(nx,ny,nz,nor,na,gt,ngt,dtp,dz,jgs,xvt(1,1,infall,il),
c     :        an,db0,imapz,mzdist,ln(il),0,xfall)
c        ENDIF
c      ENDDO
      IF ( ipconc .gt. 0 ) THEN
      DO il = lc,lhab

      IF ( ipconc .ge. ipc(il) ) THEN
      IF ( infall .eq. 2 ) THEN
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!
          DO kz = 1,nz
            DO ix = 1,nx
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
            ENDDO
          ENDDO
      ENDIF

         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)
c         call boxfall(nx,ny,nz,nor,na,gt,ngt,dtp,dz,jgs,xvt(1,1,infall,il),
c     :        an,db0,imapz,mzdist,ln(il),0,xfall)

         IF ( infall .eq. 2 .and.
     :       ( il .eq. lr .or. (il .ge. ls .and. il .le. lhl) )) THEN
c     :        .or. il .eq. lhl )) THEN

         call crowfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)
           
c           call boxfall(nx,ny,nz,nor,1,gt,ngt,dtp,dz,jgs,xvt(1,1,1,il),
c     :       tmpn,db0,imapz,mzdist,1,0,xfall0)

           DO kz=1,nz-1
           DO ix=1,nx-1
             
             an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
             
           ENDDO
           ENDDO
           
c          ENDDO

         ENDIF
        ENDIF
      ENDDO
      
      ENDIF


      ENDDO ! n=1,ndfall
      
      IF ( ndfall .gt. 1 ) THEN
        dtp = dtptmp
      ENDIF

      ENDIF ! itfall.eq.2


c
c  zero the precip rate arrays (2d)
c
c      do ix = 1,nx-istag
c      ciprt2d(ix,jy) = 0.0
c      cirrt2d(ix,jy) = 0.0
c      cwrt2d(ix,jy) = 0.0
c      rwrt2d(ix,jy) = 0.0
c      swrt2d(ix,jy) = 0.0
c      fwrt2d(ix,jy) = 0.0
c      glrt2d(ix,jy) = 0.0
c      gmrt2d(ix,jy) = 0.0
c      ghrt2d(ix,jy) = 0.0
c      hwrt2d(ix,jy) = 0.0
c      hlrt2d(ix,jy) = 0.0
c      prt2d(ix,jy) = 0.0
c      prti2d(ix,jy) = 0.0
c      prtl2d(ix,jy) = 0.0
c      end do
c      do ix = 1,nx-istag
c      cirt2d(ix,jy) = 0.0
c      end do
c
c  zero the precip flux arrays (2d)
c
      xvt(:,:,:,:) = 0.0
      
      
      do kz = 1,nz-1
      DO ix = 1,nx-1
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo
      
      vtmax = 0.0
      
      
c      do kz = 1,nz-kstag
c      do ix = 1,nx-istag
c      civt(ix,kz) = 0.0
c      end do
c      end do
c
c..Gather microphysics  
c
      if (ndebug .gt. 0 ) print*,'Gather stage'
      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do 1000 inumgs = 1,numgs 
      ngscnt = 0
      do 1001 kz = nzmpb,nz-kstag-1 
      do 1002 ix = nxmpb,nx-istag

      pqs(kz) = t00(ix,jy,kz)
c      pqs(kz) = 380.0/(pn(ix,jy,kz)+pb(kz))

      theta(kz) = an(ix,jy,kz,lt) 
      temg(kz) = t0(ix,jy,kz)
      temcg(kz) = temg(kz) - tfr
      tqvcon = temg(kz)-cbw
      ltemq = (temg(kz)-163.15)/fqsat+1.5
      qvs(kz) = pqs(kz)*tabqvs(ltemq)
      qis(kz) = pqs(kz)*tabqis(ltemq)

      qss(kz) = qvs(kz)
      if ( temg(kz) .lt. tfr ) then
c      if( qcw(kz) .le. qxmin(lc) .and. qci(kz) .gt. qxmin(li)) 
c     >  qss(kz) = qis(kz)
c      if( qcw(kz) .gt. qxmin(lc) .and. qci(kz) .gt. qxmin(li))
c     >   qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz)) /
c     >   (qcw(kz) + qci(kz))
      qss(kz) = qis(kz)
      end if
c
      if ( an(ix,jy,kz,lv)  .gt. qss(kz) .or.
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
      if ( ndebug .gt. 0 ) print*,'dbg = 5'
      
c      write(0,*) 'allocating qc
      allocate (  qx(ngscnt,lv:lhab) )
      allocate (  cx(ngscnt,lc:lhab) )
      allocate (  scx(ngscnt,lc:lhab) )
      allocate (  xv(ngscnt,lc:lhab) )
      allocate (  vtxbar(ngscnt,lc:lhab,2) )
      allocate (  xmas(ngscnt,lc:lhab) )
      allocate (  xdn(ngscnt,lc:lhab) )
      allocate (  xdia(ngscnt,lc:lhab,3) )
      
      xv(:,:)       = 0.0
      xmas(:,:)     = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:)   = 0.0
      scx(:,:)      = 0.0
c
c  define temporaries for state variables to be used in calculations
c
      do 1010 mgs = 1,ngscnt
      temp0(mgs) = tempbz(kgs(mgs))

      theta0(mgs) = ab(kgs(mgs),lt)
      thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
     > - ab(kgs(mgs),lt)
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      pres0(mgs) = pb(kgs(mgs))
      qv0(mgs) = ab(kgs(mgs),lv)
      qwvp(mgs) = an(igs(mgs),jy,kgs(mgs),lv) 
     > - ab(kgs(mgs),lv)

      presp(mgs) = pn(igs(mgs),jy,kgs(mgs))
      pres(mgs) = presp(mgs) + pres0(mgs)
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhoinv(mgs) = 1.0/rho0(mgs)
      rhovt(mgs) = Sqrt(rho00/rho0(mgs))
      pi0(mgs) = piz(kgs(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      pk(mgs)   = t77(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
      qss0(mgs) = (380.0)/(pres(mgs))
      pqs(mgs) = (380.0)/(pres(mgs))
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
c
      il5(mgs) = 0
      if ( temg(mgs) .lt. tfr ) then 
      il5(mgs) = 1
      end if
 1010 continue

c
c zero arrays that are used but not otherwise set (tm)
c
      do mgs = 1,ngscnt
         qfshr(mgs) = 0.0 
         qghshr(mgs) = 0.0 
         qglshr(mgs) = 0.0 
         qgmshr(mgs) = 0.0 
         qhlshr(mgs) = 0.0 
         qhshr(mgs) = 0.0 
       end do
c
c  set temporaries for microphysics variables
c
      DO il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      ENDDO
      end do

c      qx(mgs,lv) = an(igs(mgs),jy,kgs(mgs),lv)
c      qx(mgs,lip) = max(an(igs(mgs),jy,kgs(mgs),lip), 0.0) 
c      qx(mgs,lir) = max(an(igs(mgs),jy,kgs(mgs),lir), 0.0) 
c      qx(mgs,lc) = max(an(igs(mgs),jy,kgs(mgs),lc), 0.0) 
c      qx(mgs,lr) = max(an(igs(mgs),jy,kgs(mgs),lr), 0.0) 
c      qx(mgs,ls) = max(an(igs(mgs),jy,kgs(mgs),ls), 0.0) 
c      qx(mgs,lgl) = max(an(igs(mgs),jy,kgs(mgs),lgl), 0.0) 
c      qx(mgs,lgm) = max(an(igs(mgs),jy,kgs(mgs),lgm), 0.0) 
c      qx(mgs,lgh) = max(an(igs(mgs),jy,kgs(mgs),lgh), 0.0) 
c      qx(mgs,lf) = max(an(igs(mgs),jy,kgs(mgs),lf), 0.0) 
c      qx(mgs,lh) = max(an(igs(mgs),jy,kgs(mgs),lh), 0.0) 
c      qx(mgs,lhl) = max(an(igs(mgs),jy,kgs(mgs),lhl), 0.0) 
c      qx(mgs,li) = max(an(igs(mgs),jy,kgs(mgs),li), 0.0) 
c
c
c  set charge density variables
c
      if ( ipelec .ge. 2 ) then
      DO il = lc,lhab
      do mgs = 1,ngscnt
      scx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lsc(il))
      ENDDO
      ENDDO


      do mgs = 1,ngscnt
      cionp(mgs) = an(igs(mgs),jy,kgs(mgs),lscpi)
      cionn(mgs) = an(igs(mgs),jy,kgs(mgs),lscni)
      end do
      
C$PAR CRITICAL SECTION
c cmic$  guard  
!! c$omp  critical
      DO il = lc,lhab
      DO mgs = 1,ngscnt
       psctot1 = psctot1 + scx(mgs,il)
      ENDDO
      ENDDO
c      DO mgs = 1,ngscnt
c       psctot1 = psctot1 + 
c     :   scx(mgs,lc) + scx(mgs,li) + scx(mgs,lip) + scx(mgs,lir) +
c     :   scx(mgs,lr) + scx(mgs,ls) + scx(mgs,lgl) + scx(mgs,lgm) + scx(mgs,lgh) +
c     :   scx(mgs,lf) + scx(mgs,lh) + scx(mgs,lhl) + 
c     :   ec*(cionp(mgs) - cionn(mgs))
c      
c      ENDDO
!! c$omp  end critical
c cmic$ end guard 
C$PAR END CRITICAL SECTION
      end if
c
c  set concentrations
c
c      ssmax = 0.0
      
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
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        IF ( lccn .gt. 1 ) THEN
         ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        ELSE
         ccnc(mgs) = 0.0
        ENDIF
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

         IF ( ilimit .ge. ipc(ls) ) THEN
            tmp = (xdn0(ls)*cx(mgs,ls))/(rho0(mgs)*qx(mgs,ls))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,ls)*(tmp2)
         IF ( cnox .gt. 3.0*cno(ls) ) THEN
           cx(mgs,ls) = 3.0*cno(ls)/tmp2
         ENDIF
         ENDIF
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
c test check for intercept set max of 4.e6
            
         IF ( ilimit .ge. ipc(lgl) ) THEN
            tmp = (xdn0(lgl)*cx(mgs,lgl))/(rho0(mgs)*qx(mgs,lgl))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgl)*(tmp2)
c           IF ( qx(mgs,lgl) .gt. 0.5e-3 ) THEN
c           write(6,*) 'nogl1: ',igs(mgs),kgs(mgs),cx(mgs,lgl),rho0(mgs),
c     :      Log10(cnox),cnox*1.e-6,tmp2,qx(mgs,lgl)*1000.
c           ENDIF
         IF ( cnox .gt. 3.0*cno(lgl) ) THEN
           cx(mgs,lgl) = 3.0*cno(lgl)/tmp2
c           IF ( qx(mgs,lgl) .gt. 0.5e-3 ) THEN
c           write(6,*) 'nogl2: ',cx(mgs,lgl)
c           ENDIF
         ENDIF
         ENDIF
        ENDIF
        
        

        cx(mgs,lgm) = Max(an(igs(mgs),jy,kgs(mgs),lngm), 0.0)
        IF ( qx(mgs,lgm) .le. qxmin(lh) ) THEN
          cx(mgs,lgm) = 0.0
        ELSEIF ( cx(mgs,lgm) .eq. 0.0 .and. qx(mgs,lgm) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgm) = 0.0
        ELSE
          cx(mgs,lgm) = Max( 1.e-3, cx(mgs,lgm) )

         IF ( ilimit .ge. ipc(lgm) ) THEN
            tmp = (xdn0(lgm)*cx(mgs,lgm))/(rho0(mgs)*qx(mgs,lgm))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgm)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lgm) ) THEN
           cx(mgs,lgm) = 3.0*cno(lgm)/tmp2
         ENDIF
         ENDIF
        ENDIF

        cx(mgs,lgh) = Max(an(igs(mgs),jy,kgs(mgs),lngh), 0.0)
        IF ( qx(mgs,lgh) .le. qxmin(lh) ) THEN
          cx(mgs,lgh) = 0.0
        ELSEIF ( cx(mgs,lgh) .eq. 0.0 .and. qx(mgs,lgh) .lt. 3.0*qxmin(lh) ) THEN
          qx(mgs,lgh) = 0.0
        ELSE
          cx(mgs,lgh) = Max( 1.e-3, cx(mgs,lgh) )
         IF ( ilimit .ge. ipc(lf) ) THEN
            tmp = (xdn0(lgh)*cx(mgs,lgh))/(rho0(mgs)*qx(mgs,lgh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lgh)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lgh) ) THEN
           cx(mgs,lgh) = 3.0*cno(lgh)/tmp2
         ENDIF
         ENDIF
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
         IF ( ilimit .ge. ipc(lh) ) THEN
            tmp = (xdn0(lh)*cx(mgs,lh))/(rho0(mgs)*qx(mgs,lh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lh)*(tmp2)
         IF ( cnox .gt. 3.0*cno(lh) ) THEN
           cx(mgs,lh) = 3.0*cno(lh)/tmp2
         ENDIF
         ENDIF
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
c   
c  set temperature's
c
      if ( iptemp .ge. 1 ) then
      do mgs = 1,ngscnt
      tpcw(mgs) = an(igs(mgs),jy,kgs(mgs),ltpc) 
      tpci(mgs) = an(igs(mgs),jy,kgs(mgs),ltpi) 
      tpip(mgs) = an(igs(mgs),jy,kgs(mgs),ltpip) 
      tpir(mgs) = an(igs(mgs),jy,kgs(mgs),ltpir) 
      tprw(mgs) = an(igs(mgs),jy,kgs(mgs),ltpr) 
      tpsw(mgs) = an(igs(mgs),jy,kgs(mgs),ltps) 
      tpgl(mgs) = an(igs(mgs),jy,kgs(mgs),ltpgl) 
      tpgm(mgs) = an(igs(mgs),jy,kgs(mgs),ltpgm) 
      tpgh(mgs) = an(igs(mgs),jy,kgs(mgs),ltpgh) 
      tpfw(mgs) = an(igs(mgs),jy,kgs(mgs),ltpf) 
      tphw(mgs) = an(igs(mgs),jy,kgs(mgs),ltph)
      tphl(mgs) = an(igs(mgs),jy,kgs(mgs),ltphl)
      end do 
      end if
c
c  set tau's
c
      if ( iptime .ge. 1 ) then
      do mgs = 1,ngscnt
      rtcw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtc), 0.0) 
      if ( iptime .ge. 2 ) then
      rtci(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrti), 0.0) 
      rtip(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtip), 0.0) 
      rtir(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtir), 0.0) 
c      rtrw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtr), 0.0) 
      rtsw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrts), 0.0) 
      rtgl(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtgl), 0.0) 
      rtgm(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtgm), 0.0) 
      rtgh(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtgh), 0.0) 
      rtfw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrtf), 0.0) 
c      rthw(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrth), 0.0) 
c      rthl(mgs) = max(an(igs(mgs),jy,kgs(mgs),lrthl), 0.0) 
      ELSE
       rtci(mgs) = rimtim
       rtip(mgs) = rimtim
       rtir(mgs) = rimtim
       rtsw(mgs) = rimtim
       rtgl(mgs) = rimtim
       rtgm(mgs) = rimtim
       rtgh(mgs) = rimtim
       rtfw(mgs) = rimtim
      
      ENDIF
      end do
      end if
c
c  set factors
c
      do mgs = 1,ngscnt
c
      ssi(mgs) = qx(mgs,lv)/qis(mgs)
      ssw(mgs) = qx(mgs,lv)/qvs(mgs)
c
      tsqr(mgs) = temg(mgs)**2
c
      temgx(mgs) = min(temg(mgs),313.15)
      temgx(mgs) = max(temgx(mgs),233.15)
      felv(mgs) = 2500837.367 
     >  * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
c
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),223.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      felf(mgs) = 333690.6098 
     >  + (2030.61425)*temcgx(mgs)
     >  - (10.46708312)*temcgx(mgs)**2
c
      fels(mgs) = felv(mgs) + felf(mgs)
c
      felvs(mgs) = felv(mgs)*felv(mgs)
      felss(mgs) = fels(mgs)*fels(mgs)
c
      fgamw(mgs) = felv(mgs)*cpi/pi0(mgs)
      fgams(mgs) = fels(mgs)*cpi/pi0(mgs)
c
      fcqv1(mgs) = 4098.0258*pi0(mgs)*fgamw(mgs)
      fcqv2(mgs) = 5807.6953*pi0(mgs)*fgams(mgs)
      fcc3(mgs) = cpi*felf(mgs)/pi0(mgs)
c
      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)*
     >  (101325.0/(pres(mgs)))
c
      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*
     >  (temg(mgs)/296.0)**(1.5)
c
      fakvisc(mgs) = fadvisc(mgs)/rho0(mgs)
c
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fci(mgs) = (2.118636 + 0.007371*(temcgx(mgs)))*(1.0e+03)
c
      if ( temg(mgs) .lt. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4203.1548  + (1.30572e-2)*((temcgx(mgs)-35.)**2)
     >                 + (1.60056e-5)*((temcgx(mgs)-35.)**4)
      end if
      if ( temg(mgs) .ge. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),308.15)
      temcgx(mgs) = max(temcgx(mgs),273.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4243.1688  + (3.47104e-1)*(temcgx(mgs)**2)
      end if
c
      ftka(mgs) = tka0*fadvisc(mgs)/advisc1
      fthdf(mgs) = ftka(mgs)/(cp*rho0(mgs))
c
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))
      fpndl(mgs) = (fakvisc(mgs)/fthdf(mgs))
c
      fai(mgs) = (fels(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbi(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qis(mgs)))
      fav(mgs) = (felv(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbv(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qvs(mgs)))
c
      end do       
c
c
c   ice habit fractions
c
c
c
c  Set density
c
      if (ndebug .gt. 0 ) print*,'Set density'
c

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
c  set some values for ice nucleation
c
      if (ndebug .gt. 0 ) print*,'Set uvel etc.'
      do mgs = 1,ngscnt
      i1 = igs(mgs)
      j1 = jgs
      k1 = kgs(mgs)
      if (ndebug .gt. 10) write(0,*) 'i1,j1,k1,mgs = ',i1,j1,k1,mgs
!      uvel(mgs) = (0.5)*(u(igs(mgs),jgs,kgs(mgs))
!     >                  +u(igs(mgs)+1,jgs,kgs(mgs)))
!      vvel(mgs) = (0.5)*(v(igs(mgs),jgs,kgs(mgs))
!     >                  +v(igs(mgs),jgs+1,kgs(mgs)))
!      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)+1)
!     >                  +w(igs(mgs),jgs,kgs(mgs)))
!      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs))
!     >                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))
!      u1 = u(i1,j1,k1) ! uvel(mgs)
!      v1 = v(i1,j1,k1) ! vvel(mgs)
!      w1 = w(i1,j1,k1) ! wvel(mgs)
!      w2 = wvelkm1(mgs)
      kgsm(mgs)= max(kgs(mgs)-1,1)
      kgsp(mgs) = min(kgs(mgs)+1,nz-1)
      if (ndebug .gt. 10 ) write(0,*) 'kgsm,kgsp = ',kgsm(mgs),kgsp(mgs)
      cninm(mgs) = t7(igs(mgs),jgs,kgsm(mgs))
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      cninp(mgs) = t7(igs(mgs),jgs,kgsp(mgs))
      if (ndebug .gt. 10 ) write(0,*) 'cninm,cnina,cninp = ',cninm(mgs),cnina(mgs),cninp(mgs) 
      end do

      if (ndebug .gt. 0 ) print*,'call setvt'
c
c  Set a couple of cloud variables...
c

c      SUBROUTINE setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
c     :                 xmas,xdn,xvmn,xvmx,xv,cdx,
c     :                 ipconc,ndebug)

      call setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
     :                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,
     :                 ipconc,ndebug,ngs,nz,kgs,cwnccn,fadvisc,
     :                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,
     :                 itype1,itype2,temcg,infall-1)

      

c
c  Set number concentrations (need xdia from setvt)
c
      if ( ndebug .gt. 0 ) print*,'Set concentration'
      if ( ipconc .lt. 5 ) then
      do mgs = 1,ngscnt

      IF ( ipconc .eq. 0 ) THEN
      cx(mgs,lir) = 0.0
      if ( qx(mgs,lir) .gt. qxmin(li) )  then
      cx(mgs,lir) = cno(lir)*xdia(mgs,lir,1)
      end if
      ENDIF

      IF ( ipconc .lt. 3 ) THEN
      cx(mgs,lr) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) )  then
      cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
      end if
      ENDIF

      IF ( ipconc .lt. 4 ) THEN
      tmp = cx(mgs,ls)
      cx(mgs,ls) = 0.0
      xv(mgs,ls) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) )  then
      cx(mgs,ls) = cno(ls)*xdia(mgs,ls,1)
      xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
      end if
      ENDIF ! ( ipconc .lt. 4 )

      IF ( ipconc .lt. 5 ) THEN 

      cx(mgs,lgl) = 0.0
      xv(mgs,lgl) = 0.0
      if ( qx(mgs,lgl) .gt. qxmin(lh) )  then
      cx(mgs,lgl) = cno(lgl)*xdia(mgs,lgl,1)
      xv(mgs,lgl) = Max(xvmn(lgl), rho0(mgs)*qx(mgs,lgl)/(xdn(mgs,lgl)*cx(mgs,lgl)) )
      end if

      cx(mgs,lgm) = 0.0
      xv(mgs,lgm) = 0.0
      if ( qx(mgs,lgm) .gt. qxmin(lh) )  then
      cx(mgs,lgm) = cno(lgm)*xdia(mgs,lgm,1)
      xv(mgs,lgm) = Max(xvmn(lgm), rho0(mgs)*qx(mgs,lgm)/(xdn(mgs,lgm)*cx(mgs,lgm)) )
      end if

      cx(mgs,lgh) = 0.0
      xv(mgs,lgh) = 0.0
      if ( qx(mgs,lgh) .gt. qxmin(lh) )  then
      cx(mgs,lgh) = cno(lgh)*xdia(mgs,lgh,1)
      xv(mgs,lgh) = Max(xvmn(lgh), rho0(mgs)*qx(mgs,lgh)/(xdn(mgs,lgh)*cx(mgs,lgh)) )
      end if

      cx(mgs,lf) = 0.0
      xv(mgs,lf) = 0.0
      if ( qx(mgs,lf) .gt. qxmin(lh) )  then
      cx(mgs,lf) = cno(lf)*xdia(mgs,lf,1)
      xv(mgs,lf) = Max(xvmn(lf), rho0(mgs)*qx(mgs,lf)/(xdn(mgs,lf)*cx(mgs,lf)) )
      end if

      cx(mgs,lh) = 0.0
      xv(mgs,lh) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) )  then
      cx(mgs,lh) = cno(lh)*xdia(mgs,lh,1)
      xv(mgs,lh) = Max(xvmn(lh), rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh)) )
      end if

      cx(mgs,lhl) = 0.0
      xv(mgs,lhl) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lh) )  then
      cx(mgs,lhl) = cno(lhl)*xdia(mgs,lhl,1)
      xv(mgs,lhl) = Max(xvmn(lhl), rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl)) )
      end if

      ENDIF ! ( ipconc .lt. 5 )

      end do
      end if
      
      IF ( ipconc .ge. 2 ) THEN
      DO mgs = 1,ngscnt
        rb(mgs) = 0.5*xdia(mgs,lc,1)*(Sqrt(1./(1.+cnu)))**(1./3.)
        xl2p(mgs) = Max(0.0, 2.7e-2*xdn(mgs,lc)*cx(mgs,lc)*xv(mgs,lc)*   
     :           ((0.5e20*rb(mgs)**3*xdia(mgs,lc,1))-0.4) )
        IF ( rb(mgs) .gt. 3.51e-6 ) THEN
!          rh(mgs) = Max( 0.5d0*xdia(mgs,lc,1), 6.3d-4/(1.d6*(rb(mgs) - 3.5d-6)) )
          rh(mgs) = Max( 41.e-6, 6.3e-4/(1.e6*(rb(mgs) - 3.5e-6)) )
        ELSE
          rh(mgs) = 41.e-6
        ENDIF
        IF ( xl2p(mgs) .gt. 0.0 ) THEN
          nh(mgs) = 4.2e9*xl2p(mgs)
        ELSE
          nh(mgs) = 1.e30
        ENDIF
      ENDDO
      ENDIF
      
c
c cloud ice variables
c
c     if ( ihabdo .lt. 1 ) then
c     do mgs = 1,ngscnt
c     cx(mgs,lir) = min(cnit*exp(-temcg(mgs)*bta1),1.e+09)
c     cx(mgs,lir) = cnina(mgs)
c     cx(mgs,lir) = max(1.0,cx(mgs,lir))
c     xmas(mgs,lir) =
c    >  min( max(qx(mgs,lir)*rho0(mgs)/cx(mgs,lir),cimasn),cimasx )
c     if ( temcg(mgs) .gt. 0 ) then
c     xdia(mgs,lir,1) = 0.0
c     else
c     xdia(mgs,lir,1) = 16.7*(xmas(mgs,lir)**(0.5))
c     xdia(mgs,lir,1) = max(xdia(mgs,lir,1), 1.e-5)
c     xdia(mgs,lir,1) = min(xdia(mgs,lir,1), 1.e-3)
c     endif
c     xdia(mgs,lir,2) = xdia(mgs,lir,1)**2
c     vtxbar(mgs,lir,1) = (131.6)*(xdia(mgs,lir,1)**(0.824))
c     vtxbar(mgs,lir,1) = vtxbar(mgs,lir,1)*rhovt(mgs)
c     circap(mgs)  = (0.5)*xdia(mgs,lir,1)
c     cirat(mgs) = (0.25)*pi*xdia(mgs,lir,2)*cx(mgs,lir)
c     end do
c     endif
c
c
c
c              
c
c  maximum depletion tendency by any one source
c
c
      if( ndebug .gt. 0 ) THEN
       write(iunit,*) 'Set depletion max/min1'
       call flush(iunit)
      endif
      do mgs = 1,ngscnt
c      qvmxd(mgs)  = 0.10*(qx(mgs,lv)-qvs(mgs))/dtp
      qvimxd(mgs) = 0.70*(qx(mgs,lv)-qis(mgs))/dtp ! depletion by all vap. dep to ice.
c      qvmxd(mgs)  = max(qvmxd(mgs), 0.0)
      qvimxd(mgs) = max(qvimxd(mgs), 0.0)
      qimxd(mgs)  = 0.10*qx(mgs,li)/dtp
      qirmxd(mgs) = 0.10*qx(mgs,lir)/dtp
      qipmxd(mgs) = 0.10*qx(mgs,lip)/dtp
      qcmxd(mgs)  = 0.10*qx(mgs,lc)/dtp
      qrmxd(mgs)  = 0.10*qx(mgs,lr)/dtp
      qsmxd(mgs)  = 0.10*qx(mgs,ls)/dtp
      qglmxd(mgs) = 0.10*qx(mgs,lgl)/dtp
      qgmmxd(mgs) = 0.10*qx(mgs,lgm)/dtp
      qghmxd(mgs) = 0.10*qx(mgs,lgh)/dtp
      qhmxd(mgs)  = 0.10*qx(mgs,lh)/dtp
      qhlmxd(mgs) = 0.10*qx(mgs,lhl)/dtp
      qfmxd(mgs)  = 0.10*qx(mgs,lf)/dtp
      end do
c
      if( ndebug .gt. 0 ) THEN
       write(iunit,*) 'Set depletion max/min2'
       call flush(iunit)
      endif

      do mgs = 1,ngscnt
c  
      if ( qx(mgs,lc) .le. qxmin(lc) ) then
      ccmxd(mgs)  = 0.10*cx(mgs,lc)/dtp
      else
      IF ( ipconc .ge. 2 ) THEN
        ccmxd(mgs)  = 0.10*cx(mgs,lc)/dtp
      ELSE
        ccmxd(mgs)  = 0.10*qx(mgs,lc)/(xmas(mgs,lc)*rho0(mgs)*dtp)
      ENDIF
      end if
c
      if ( qx(mgs,li) .le. qxmin(li) ) then
      cimxd(mgs)  = 0.10*cx(mgs,li)/dtp
      else
      IF ( ipconc .ge. 1 ) THEN
        cimxd(mgs)  = 0.10*cx(mgs,li)/dtp
      ELSE
        cimxd(mgs)  = 0.10*qx(mgs,li)/(xmas(mgs,li)*rho0(mgs)*dtp)
      ENDIF
      end if
c
      if ( qx(mgs,lip) .le. qxmin(lip) ) then
      cipmxd(mgs)  = 0.10*cx(mgs,lip)/dtp
      else
      IF ( ipconc .ge. 1 ) THEN
        cipmxd(mgs)  = 0.10*cx(mgs,lip)/dtp
      ELSE
      cipmxd(mgs)  = 0.10*qx(mgs,lip)/(xmas(mgs,lip)*rho0(mgs)*dtp)
      ENDIF
      end if
c
c      if ( qx(mgs,lir) .le. qxmin(lir) ) then
      cirmxd(mgs)  = 0.10*cx(mgs,lir)/dtp
c      else
c      cirmxd(mgs)  = 0.10*qx(mgs,lir)/(xmas(mgs,lir)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,lr) .le. qxmin(lr) ) then
      crmxd(mgs)  = 0.10*cx(mgs,lr)/dtp
c      else
c      crmxd(mgs)  = 0.10*qx(mgs,lr)/(xmas(mgs,lr)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,ls) .le. qxmin(ls) ) then
      csmxd(mgs)  = 0.10*cx(mgs,ls)/dtp
c      else
c      csmxd(mgs)  = 0.10*qx(mgs,ls)/(xmas(mgs,ls)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,lgl) .le. qxmin(lgl) ) then
      cglmxd(mgs)  = 0.10*cx(mgs,lgl)/dtp
c      else
c      cglmxd(mgs)  = 0.10*qx(mgs,lgl)/(xmas(mgs,lgl)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,lgm) .le. qxmin(lgm) ) then
      cgmmxd(mgs)  = 0.10*cx(mgs,lgm)/dtp
c      else
c      cgmmxd(mgs)  = 0.10*qx(mgs,lgm)/(xmas(mgs,lgm)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,lgh) .le. qxmin(lgh) ) then
      cghmxd(mgs)  = 0.10*cx(mgs,lgh)/dtp
c      else
c      cghmxd(mgs)  = 0.10*qx(mgs,lgh)/(xmas(mgs,lgh)*rho0(mgs)*dtp)
c      end if
c  
c      if ( qx(mgs,lh) .le. qxmin(lh) ) then
      chmxd(mgs)  = 0.10*cx(mgs,lh)/dtp
c      else
c      chmxd(mgs)  = 0.10*qx(mgs,lh)/(xmas(mgs,lh)*rho0(mgs)*dtp)
c      end if
c
c      if ( qx(mgs,lf) .le. qxmin(lf) ) then
      cfmxd(mgs)  = 0.10*cx(mgs,lf)/dtp
c      else
c      cfmxd(mgs)  = 0.10*qx(mgs,lf)/(xmas(mgs,lf)*rho0(mgs)*dtp)
c      end if
      end do
c
      if( ndebug .gt. 0 ) THEN
       write(iunit,*) 'Setup inductive charging...'
       call flush(iunit)
      endif
 
c
c  set-up stuff inductive space charge parameterization:
c     This parameterization is based on Ziegler, Macgorman et al...
c
c

      if ( ipelec .ge. 1 ) then
c
c  vertical component of e
c
      do mgs = 1,ngscnt
      dezcomp(mgs) = Max(-1.00e5,
     :   Min(elec(iez)%flt3d(igs(mgs),jgs,kgs(mgs)),1.00e5 )) 
      end do
c
c  others...
c
      if ( nonigrd .ge. 2 ) then
      do mgs = 1,ngscnt
c
c
c  temperature factor:
c      if trever = -21: reduces to jayratene
c      if trever = -10: reduces to takahashi
c
      temele = max(-25.0,((-21.0)/trever)*temcg(mgs))
c
      sctem(mgs)  =
     >   (-1.7e-05)*(temele**3)
     > + (-3.0e-03)*(temele**2)
     > + (-5.0e-02)*(temele**1)
     > + ( 1.3e-01)
c
c  critical (threshold) liquid water content factor (1.0e-04 kg/m**3)
c
      if ( temcg(mgs) .gt. trever ) then
      dellwc(mgs) =
     >   (1000.0)*(rho0(mgs)*qx(mgs,lc)-(1.0e-04))
      end if
      if ( temcg(mgs) .lt. trever ) then
      dellwc(mgs) =
     >   (1000.0)*(rho0(mgs)*qx(mgs,lc))
      end if
      if ( qx(mgs,lc) .lt. 1.e-6 ) dellwc(mgs) = 0.0
c
c parabolic function to reduce charging at low temperature
c
      IF ( temcg(mgs) .lt. -30 .and. temcg(mgs) .gt. -43.0 ) THEN
        sctem(mgs) = sctem(mgs)*
     :            (1.0 - ((temcg(mgs)+30.0)/(43.0 - 30.0))**2)
      ELSEIF ( temcg(mgs) .le. -43.0 ) THEN
        sctem(mgs) = 0.0
      END IF

      end do
      end if
c
      end if
c
      if( ndebug .gt. 0 ) THEN
       write(iunit,*) 'Done with setup inductive charging'
       call flush(iunit)
      endif

c#endif
c
c
c
c
c  microphysics source terms (1/s) for mixing ratios 
c
c
c
c  Collection efficiencies:
c
c 
      if (ndebug .gt. 0 ) print*,'Set collection efficiencies'
c

      exy(:,:,:) = 0.0

      do mgs = 1,ngscnt
c
      erip(mgs) = 0.0
      esip(mgs) = 0.0
      eipw(mgs) = 0.0
      eipir(mgs) = 0.0
      eglip(mgs) = 0.0
      egmip(mgs) = 0.0
      eghip(mgs) = 0.0
      ehip(mgs) = 0.0
      ehlip(mgs) = 0.0
      efip(mgs) = 0.0
c
      eiri(mgs) = 0.0
      erir(mgs) = 0.0
      esir(mgs) = 0.0
      eirw(mgs) = 0.0
      eirir(mgs) = 0.0
      eglir(mgs) = 0.0
      egmir(mgs) = 0.0
      eghir(mgs) = 0.0
      ehir(mgs) = 0.0
      ehlir(mgs) = 0.0
      efir(mgs) = 0.0
c
      erw(mgs) = 0.0
      esw(mgs) = 0.0
      efw(mgs) = 0.0
      ehw(mgs) = 0.0
      ehlw(mgs) = 0.0
      ehxw(mgs) = 0.0
      eglw(mgs) = 0.0
      egmw(mgs) = 0.0
      eghw(mgs) = 0.0
c
      err(mgs) = 0.0
      esr(mgs) = 0.0
      efr(mgs) = 0.0
      ehr(mgs) = 0.0
      ehlr(mgs) = 0.0
      ehxr(mgs) = 0.0
      eglr(mgs) = 0.0
      egmr(mgs) = 0.0
      eghr(mgs) = 0.0
c
      eri(mgs) = 0.0
      esi(mgs) = 0.0
      efi(mgs) = 0.0
      ehi(mgs) = 0.0
      ehli(mgs) = 0.0
      ehxi(mgs) = 0.0
      egli(mgs) = 0.0
      egmi(mgs) = 0.0
      eghi(mgs) = 0.0
c
      ers(mgs) = 0.0
      ess(mgs) = 0.0
      efs(mgs) = 0.0
      ehs(mgs) = 0.0
      ehls(mgs) = 0.0
      ehxs(mgs) = 0.0
      egls(mgs) = 0.0
      egms(mgs) = 0.0
      eghs(mgs) = 0.0
c
      eipi(mgs) = 0.0
      eipw(mgs) = 0.0
      eipip(mgs) = 0.0
      eiw(mgs) = 0.0
      eii(mgs) = 0.0
c
c
c  Ice-Ice: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,li).gt.qxmin(li) ) then
      eii(mgs)=0.1*exp(0.1*temcg(mgs))
      if ( temg(mgs) .lt. 243.15 .and. qx(mgs,lc) .gt. 1.e-6 ) then
      eii(mgs)=0.1
      end if
      if ( temg(mgs) .gt. 273.15 ) eii(mgs) = 0.0
      end if
c
      if ( qx(mgs,lip).gt.qxmin(li) .and. qx(mgs,li) .gt. qxmin(li) ) then
      eipi(mgs)=0.1*exp(0.1*temcg(mgs))
      if ( qx(mgs,lc) .gt. 1.e-6 ) then
      eipi(mgs)=0.1
      end if
      if ( temg(mgs) .gt. 273.15 ) eipi(mgs) = 0.0
      end if
c
      if ( qx(mgs,lip).gt.qxmin(li) ) then
      eipip(mgs)=0.1*exp(0.1*temcg(mgs))
      if ( qx(mgs,lc) .gt. 1.e-6 ) then
      eipip(mgs)=0.1
      end if
      if ( temg(mgs) .lt. 264.15 .and. temg(mgs) .gt. 251.15 ) then
      if ( qx(mgs,lc) .gt. 1.e-6 ) then
      eipip(mgs) = 1.4
      end if
      end if
      if ( temg(mgs) .gt. 273.15 ) eipip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lip).gt.qxmin(li) .and. qx(mgs,lir) .gt. qxmin(li) ) then
      eipir(mgs)=0.1*exp(0.1*temcg(mgs))
      if ( qx(mgs,lc) .gt. 1.e-6 ) then
      eipir(mgs)=0.1
      end if
      if ( temg(mgs) .gt. 273.15 ) eipir(mgs) = 0.0
      end if
c
      if ( qx(mgs,lir).gt.qxmin(li) ) then
      eirir(mgs)=0.1*exp(0.1*temcg(mgs))
      if ( qx(mgs,lc) .gt. 1.e-6 ) then
      eirir(mgs)=0.1
      end if
       IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
        eiri(mgs) = eirir(mgs)
       ENDIF
      if ( temg(mgs) .gt. 273.15 ) eirir(mgs) = 0.0
      end if
c
c
c  Ice-cloud water: Collection (cxc) efficiencies
c
c
      eiw(mgs) = 0.0
      if ( qx(mgs,li).gt.qxmin(li) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      if (xdia(mgs,lc,1).gt.15.0e-06 .and. xdia(mgs,li,1).gt.30.0e-06) then
      eiw(mgs) = 0.5
      end if
      if ( temg(mgs) .gt. 273.15 ) eiw(mgs) = 0.0
      end if
c
      eipw(mgs) = 0.0
      if ( qx(mgs,lip).gt.qxmin(li) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      if (xdia(mgs,lc,1).gt.15.0e-06 .and. xdia(mgs,lip,1).gt.300.0e-06) then
      eipw(mgs) = 0.5
      end if
      if ( temg(mgs) .gt. 273.15 ) eipw(mgs) = 0.0
      end if
c
      eirw(mgs) = 0.0
      if ( qx(mgs,lir).gt.qxmin(li) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      if (xdia(mgs,lc,1).gt.15.0e-06 .and. xdia(mgs,lir,1).gt.30.0e-06) then
      eirw(mgs) = 0.5
      end if
      if ( temg(mgs) .gt. 273.15 ) eirw(mgs) = 0.0
      end if
c
c
c  Rain: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      erw(mgs) =
     >  min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
       IF ( xdia(mgs,lc,1) .lt. 2.4e-06 .or. xdia(mgs,lr,1) .le. 50.0e-6 ) THEN
          erw(mgs)=0.0
       ENDIF
      end if
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lr).gt.qxmin(lr) ) then
      err(mgs)=1.0
      end if
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      ers(mgs)=1.0
      end if
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lir).gt.qxmin(lir) ) then
      erir(mgs)=1.0
      if ( xdia(mgs,lir,1) .lt. 40.e-6 ) erir(mgs)=0.0
      end if
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,li).gt.qxmin(li) ) then
      eri(mgs)=1.0
      if ( xdia(mgs,li,1) .lt. 40.e-6 ) eri(mgs)=0.0
      end if
c
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lip).gt.qxmin(li) ) then
      erip(mgs)=1.0
      if ( xdia(mgs,lip,1) .lt. 40.e-6 ) erip(mgs)=0.0
      end if
c
c
c  Snow aggregates: Collection (cxc) efficiencies
c
c Modified by ERM with a linear function for small droplets and large
c snow agg. based numerical data from Wang and Ji (1992) in P&K 1997, which
c allows collection of very small droplets, albeit at low efficiency.  But slow 
c fall speeds of snow make up for the efficiency.
c
      esw(mgs) = 0.0
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lc).gt.qxmin(lc) ) then
        if ( xdia(mgs,lc,1) .gt. 15.e-6 .and. xdia(mgs,ls,1) .gt. 100.e-6) then
          esw(mgs) = 0.5
        ELSEIF ( xdia(mgs,ls,1) .ge. 500.e-6 ) THEN
          esw(mgs) = Min(0.5, 0.05 + (0.8-0.05)/(40.e-6)*xdia(mgs,lc,1) )
        ENDIF
      end if
c
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lr).gt.qxmin(lr) 
c     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      esr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,ls,1))
      end if
c
      if ( qx(mgs,ls).gt.qxmin(ls) ) then
      IF ( ipconc .lt. 4 ) THEN
        ess(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
        ess(mgs)=min(0.1,ess(mgs))
      ELSE
        ess(mgs) = Exp(0.05*Min( temcg(mgs), 0.0 ) )
      ENDIF
      end if
c
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,li).gt.qxmin(li) ) then
      esi(mgs) = 0.0
      IF ( ipconc .lt. 4 ) THEN
      esi(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
      esi(mgs)=min(0.1,esi(mgs))
      ELSE
      esi(mgs)=0.5*exp(0.1*min(temcg(mgs),0.0))
      ENDIF
      if ( temg(mgs) .gt. 273.15 ) esi(mgs) = 0.0
      end if
c
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lir).gt.qxmin(li) ) then
      esir(mgs) = 0.0
      esir(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
      esir(mgs)=min(0.1,esir(mgs))
      if ( temg(mgs) .gt. 273.15 ) esir(mgs) = 0.0
      end if
c
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lip).gt.qxmin(li) ) then
      esip(mgs) = 0.0
      esip(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
      esip(mgs)=min(0.1,esip(mgs))
      if ( temg(mgs) .gt. 273.15 ) esip(mgs) = 0.0
      end if
c
c
c  Graupel-low density: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      eglw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) eglw(mgs)=0.0
      end if
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      eglr(mgs)=1.0*Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lgl,1))
      end if
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      egls(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,li).gt.qxmin(li) ) then
      egli(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) egli(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,lip).gt.qxmin(li) ) then
      eglip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) eglip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgl).gt.qxmin(lgl) .and. qx(mgs,lir).gt.qxmin(li) ) then
      eglir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) eglir(mgs) = 0.0
      end if
c
c
c  Graupel-medium density: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      egmw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) egmw(mgs)=0.0
      end if
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      egmr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lgm,1))
      end if
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      egms(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,li).gt.qxmin(li) ) then
      egmi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) egmi(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,lip).gt.qxmin(li) ) then
      egmip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) egmip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgm).gt.qxmin(lgm) .and. qx(mgs,lir).gt.qxmin(li) ) then
      egmir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) egmir(mgs) = 0.0
      end if
c
c
c  Graupel-high density: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      eghw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) eghw(mgs)=0.0
      end if
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      eghr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lgh,1))
      end if
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      eghs(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,li).gt.qxmin(li) ) then
      eghi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) eghi(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,lir).gt.qxmin(li) ) then
      eghir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) eghir(mgs) = 0.0
      end if
c
      if ( qx(mgs,lgh).gt.qxmin(lgh) .and. qx(mgs,lip).gt.qxmin(li) ) then
      eghip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) eghip(mgs) = 0.0
      end if
c
c
c  Frozen drops: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      efw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) efw(mgs)=0.0
      end if
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      efr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lf,1))
      end if
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      efs(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,li).gt.qxmin(li) ) then
      efi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) efi(mgs) = 0.0
      end if
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lip).gt.qxmin(li) ) then
      efip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) efip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lf).gt.qxmin(lf) .and. qx(mgs,lir).gt.qxmin(li) ) then
      efir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) efir(mgs) = 0.0
      end if
c
c
c  Hail: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      ehw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehw(mgs)=0.0
      end if
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      ehr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lh,1))
      end if
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      ehs(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehi(mgs) = 0.0
      end if
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lip).gt.qxmin(li) ) then
      ehip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lir).gt.qxmin(li) ) then
      ehir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehir(mgs) = 0.0
      end if
c
c
c  Hail-large: Collection (cxc) efficiencies
c
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      ehlw(mgs) =
     >  ewfac*min((aradcw + cwrad*(bradcw + cwrad*
     <  (cradcw + cwrad*(dradcw)))), 1.0)
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehlw(mgs)=0.0
      end if
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,lr).gt.qxmin(lr) 
     :     .and. temg(mgs) .lt. tfr 
     :                               ) then
      ehlr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lhl,1))
      end if
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      ehls(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      end if
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehli(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehlir(mgs) = 0.0
      end if
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,lip).gt.qxmin(li) ) then
      ehlip(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehlip(mgs) = 0.0
      end if
c
      if ( qx(mgs,lhl).gt.qxmin(lh) .and. qx(mgs,lir).gt.qxmin(li) ) then
      ehlir(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      if ( temg(mgs) .gt. 273.15 ) ehlir(mgs) = 0.0
      end if
c
      end do
c
c
c
c  Set flags for plates vs. columns
c
c
      do mgs = 1,ngscnt
c
      xplate(mgs) = 0.0 
      xcolmn(mgs) = 0.0
c
      if ( temcg(mgs) .lt. 0. .and. temcg(mgs) .ge. -4. ) then
      xplate(mgs) = 1.0 
      xcolmn(mgs) = 0.0
      end if
c
      if ( temcg(mgs) .lt. -4. .and. temcg(mgs) .ge. -9. ) then
      xplate(mgs) = 0.0 
      xcolmn(mgs) = 1.0
      end if
c
      if ( temcg(mgs) .lt. -9. .and. temcg(mgs) .ge. -22.5 ) then
      xplate(mgs) = 1.0 
      xcolmn(mgs) = 0.0
      end if
c
      if ( temcg(mgs) .lt. -22.5 .and. temcg(mgs) .ge. -90. ) then
      xplate(mgs) = 0.0 
      xcolmn(mgs) = 1.0
      end if
c
      end do
c
c
c
c  Collection growth equations....
c
c
      if (ndebug .gt. 0 ) print*,'Collection: rain collects xxxxx'      
c
      do mgs = 1,ngscnt
      qracw(mgs) =  0.0
      IF ( qx(mgs,lr) .gt. qxmin(lr) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
      IF ( ipconc .lt. 3 ) THEN
       IF ( erw(mgs) .gt. 0.0 ) THEN
       qracw(mgs) = 
     >   (0.25)*pi*erw(mgs)*qx(mgs,lc)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lr,2) 
     >    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,lc,1) 
     >    + gf1*xdia(mgs,lc,2) )  
       ENDIF
      ELSE
       rwrad = 0.5*xdia(mgs,lr,1)
        IF ( rwrad .gt. rh(mgs) ) THEN ! .or. cx(mgs,lr) .gt. nh(mgs) ) THEN
         IF ( rwrad .gt. 50.e-6 ) THEN
           qracw(mgs) = aa2*cx(mgs,lr)*cx(mgs,lc)*xmas(mgs,lc)*
     :        ((cnu + 2)*xv(mgs,lc)/(cnu + 1) + xv(mgs,lr))/rho0(mgs)
         ELSE
           qracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*xdn(mgs,lc)*
     :        ((cnu + 3)*(cnu + 2)*xv(mgs,lc)**3/(cnu + 1)**2 + 
     :         (rnu + 2)*xv(mgs,lc)*xv(mgs,lr)**2/(rnu + 1))/rho0(mgs)
         ENDIF
        ENDIF
       ENDIF
       qracw(mgs) = Min(qracw(mgs), qcmxd(mgs))
       ENDIF
      end do
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qraci(mgs) =  0.0
      IF ( eri(mgs) .gt. 0.0 ) THEN
      qraci(mgs) = 
     >   min(
     >   (0.25)*pi*eri(mgs)*qx(mgs,li)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lr,2) 
     >    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1) 
     >    + gf1*xdia(mgs,li,2) )  
     <  , qimxd(mgs))
      if ( temg(mgs) .gt. 268.15 ) then
      qraci(mgs) = 0.0
      end if
      ENDIF
      end do
      end if
c
      do mgs = 1,ngscnt
      qracip(mgs) =  0.0
      IF ( erip(mgs) .gt. 0.0 ) THEN
      qracip(mgs) = 
     >   min(
     >   (0.25)*pi*erip(mgs)*qx(mgs,lip)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2) 
     >    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,lip,1) 
     >    + gf3*xdia(mgs,lr,2) )  
     <  , qipmxd(mgs))
      ENDIF
      if ( temg(mgs) .gt. 268.15 ) then
      qracip(mgs) = 0.0
      end if
      end do
c
      do mgs = 1,ngscnt
      qracir(mgs) =  0.0
      IF ( erir(mgs) .gt. 0.0 ) THEN
      qracir(mgs) = 
     >   min(
     >   (0.25/gf4)*pi*erir(mgs)*qx(mgs,lir)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2) 
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lir,1) 
     >    + gf4*gf3*xdia(mgs,lr,2) )  
     <  , qirmxd(mgs))
      if ( temg(mgs) .gt. 268.15 ) then
      qracir(mgs) = 0.0
      end if
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qracs(mgs) =  0.0
      IF ( ers(mgs) .gt. 0.0 ) THEN
      qracs(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*ers(mgs)*qx(mgs,ls)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lr,1) 
     >    + gf4*gf3*xdia(mgs,lr,2) ) 
     <  , qsmxd(mgs))
      ENDIF
      end do
c
c      do mgs = 1,ngscnt
c      qracf(mgs) =  0.0
c      IF ( efr(mgs) .gt. 0.0 ) THEN
c      qracf(mgs) = 
c     >   min(
c     >   ((0.25)*pi/gf4)*efr(mgs)*qx(mgs,lf)*cx(mgs,lr)
c     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lf,1))
c     >  *(  gf6*gf1*xdia(mgs,lf,2)
c     >    + 2.0*gf5*gf2*xdia(mgs,lf,1)*xdia(mgs,lr,1) 
c     >    + gf4*gf3*xdia(mgs,lr,2) ) 
c     <  , qfmxd(mgs))
c      ENDIF
c      end do
c
      do mgs = 1,ngscnt
      qracgl(mgs) =  0.0
      IF ( eglr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgl,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgl,1)*vtxbar(mgs,lr,1) )
      qracgl(mgs) = 
     >   ((0.25)*pi)*eglr(mgs)*qx(mgs,lgl)*cx(mgs,lr)*vt*
     :      (0.698740*xdia(mgs,lr,2) + 
     :      4.9600413*xdia(mgs,lgl,1)*xdia(mgs,lr,1)  + 
     :       20.*xdia(mgs,lgl,2))
      qracgl(mgs) = Min( qracgl(mgs), qglmxd(mgs) )
      ELSE
      qracgl(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*eglr(mgs)*qx(mgs,lgl)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lgl,1))
     >  *(  gf6*gf1*xdia(mgs,lgl,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lgl,1)*xdia(mgs,lr,1) 
     >    + gf4*gf3*xdia(mgs,lr,2) ) 
     <  , qglmxd(mgs))
      ENDIF
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qracgm(mgs) =  0.0
      IF ( egmr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgm,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgm,1)*vtxbar(mgs,lr,1) )
      qracgm(mgs) = 
     >   ((0.25)*pi)*egmr(mgs)*qx(mgs,lgm)*cx(mgs,lr)*vt*
     :      (0.698740*xdia(mgs,lr,2) + 
     :      4.9600413*xdia(mgs,lgm,1)*xdia(mgs,lr,1)  + 
     :       20.*xdia(mgs,lgm,2))
      qracgm(mgs) = Min( qracgm(mgs), qgmmxd(mgs) )
      ELSE
      qracgm(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*egmr(mgs)*qx(mgs,lgm)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lgm,1))
     >  *(  gf6*gf1*xdia(mgs,lgm,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lr,1) 
     >    + gf4*gf3*xdia(mgs,lr,2) ) 
     <  , qgmmxd(mgs))
      ENDIF
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qracgh(mgs) =  0.0
      IF ( eghr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgh,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgh,1)*vtxbar(mgs,lr,1) )
      qracgh(mgs) = 
     >   ((0.25)*pi)*eghr(mgs)*qx(mgs,lgh)*cx(mgs,lr)*vt*
     :      (0.698740*xdia(mgs,lr,2) + 
     :      4.9600413*xdia(mgs,lgh,1)*xdia(mgs,lr,1)  + 
     :       20.*xdia(mgs,lgh,2))
      qracgh(mgs) = Min( qracgh(mgs), qghmxd(mgs) )
      ELSE
      qracgh(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*eghr(mgs)*qx(mgs,lgh)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lgh,1))
     >  *(  gf6*gf1*xdia(mgs,lgh,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lgh,1)*xdia(mgs,lr,1) 
     >    + gf4*gf3*xdia(mgs,lr,2) ) 
     <  , qghmxd(mgs))
      ENDIF
      ENDIF
      end do
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: snow collects xxxxx'
c
      do mgs = 1,ngscnt
      qsacw(mgs) =  0.0
      IF ( esw(mgs) .gt. 0.0 ) THEN
      qsacw(mgs) =
     >   min( 
     >   ((0.25)*pi)*esw(mgs)*qx(mgs,lc)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,ls,2) 
     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,lc,1) 
     >    + gf1*xdia(mgs,lc,2) )  
     <  , qcmxd(mgs))
      ENDIF
      end do
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qsaci(mgs) = 0.0
      IF ( esi(mgs) .gt. 0.0 ) THEN
      qsaci(mgs) = 
     >   min(
     >   ((0.25)*pi)*esi(mgs)*qx(mgs,li)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,ls,2) 
     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,li,1) 
     >    + gf1*xdia(mgs,li,2) )  
     <  , qimxd(mgs))
      ENDIF
      end do
      end if
c
c
      do mgs = 1,ngscnt
      qsacip(mgs) = 0.0
      IF ( esip(mgs) .gt. 0.0 ) THEN
      qsacip(mgs) = 
     >   min(
     >   ((0.25)*pi)*esip(mgs)*qx(mgs,lip)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2) 
     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,lip,1) 
     >    + gf3*xdia(mgs,ls,2) )  
     <  , qipmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qsacir(mgs) = 0.0
      IF ( esir(mgs) .gt. 0.0 ) THEN
      qsacir(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*esir(mgs)*qx(mgs,lir)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2) 
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lir,1) 
     >    + gf4*gf3*xdia(mgs,ls,2) )  
     <  , qirmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qsacr(mgs) = 0.0
      csacr(mgs) = 0.0
      IF ( esr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,ls,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,ls,1)*vtxbar(mgs,lr,1) )
       qsacr(mgs) = esr(mgs)*cx(mgs,ls)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,ls,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,ls,2))
        qsacr(mgs) = Min( qsacr(mgs), qrmxd(mgs) )
        csacr(mgs) = qsacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        csacr(mgs) = min(csacr(mgs),crmxd(mgs))
      ELSE
       qsacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*esr(mgs)*qx(mgs,lr)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,ls,1) 
     >    + gf4*gf3*xdia(mgs,ls,2) ) 
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do
c
c
c
      do mgs = 1,ngscnt
      qipaci(mgs) = 0.0
      cipaci(mgs) = 0.0
      IF ( eipi(mgs) .gt. 0.0 ) THEN
       vt = Sqrt((vtxbar(mgs,li,1)-vtxbar(mgs,lip,1))**2 + 
     :            0.04*vtxbar(mgs,li,1)*vtxbar(mgs,lip,1) )
      qipaci(mgs) = 
     >   min(
     >   ((0.25)*pi)*eipi(mgs)*qx(mgs,li)*cx(mgs,lip)
     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,li,1))
     >  *(  xdia(mgs,lip,2) 
     >    + 2.0*xdia(mgs,lip,1)*xdia(mgs,li,1) 
     >    + xdia(mgs,li,2) )  
     <  , qimxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
      cipaci(mgs) = Min(cimxd(mgs),qipaci(mgs)*rho0(mgs)/xmas(mgs,li))
      ENDIF
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qiacip(mgs) = 0.0
      ciacip(mgs) = 0.0
      IF ( eipi(mgs) .gt. 0.0 ) THEN
       vt = Sqrt((vtxbar(mgs,li,1)-vtxbar(mgs,lip,1))**2 + 
     :            0.04*vtxbar(mgs,li,1)*vtxbar(mgs,lip,1) )
      qiacip(mgs) = 
     >   min(
     >   ((0.25)*pi)*eipi(mgs)*qx(mgs,lip)*cx(mgs,li)
     :   *vt
c     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,li,1))
     >  *(  xdia(mgs,lip,2) 
     >    + 2.0*xdia(mgs,lip,1)*xdia(mgs,li,1) 
     >    + xdia(mgs,li,2) )  
     <  , qipmxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
      ciacip(mgs) = Min(cipmxd(mgs),qiacip(mgs)*rho0(mgs)/xmas(mgs,lip))
      ENDIF
      ENDIF
      end do
c
c
c
      do mgs = 1,ngscnt
      qiracip(mgs) = 0.0
      ciracip(mgs) = 0.0
      IF ( eipir(mgs) .gt. 0.0 ) THEN
      qiracip(mgs) = 
     >   min(
     >   ((0.25)*pi)*eipir(mgs)*qx(mgs,lip)*cx(mgs,lir)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2) 
     >    + 2.0*gf2*xdia(mgs,lir,1)*xdia(mgs,lip,1) 
     >    + gf3*xdia(mgs,lir,2) )  
     <  , qipmxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
      ciracip(mgs) = Min(cipmxd(mgs),qiracip(mgs)*rho0(mgs)/xmas(mgs,lip))
      ENDIF
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qipacir(mgs) = 0.0
      cipacir(mgs) = 0.0
      IF ( eipir(mgs) .gt. 0.0 ) THEN
      qipacir(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*eipir(mgs)*qx(mgs,lir)*cx(mgs,lip)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,lip,1))
     >  *(  gf6*xdia(mgs,lir,2) 
     >    + 2.0*gf5*xdia(mgs,lir,1)*xdia(mgs,lip,1) 
     >    + gf4*xdia(mgs,lip,2) )  
     <  , qirmxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
        cipacir(mgs) = 
     :    Min(cirmxd(mgs),qipacir(mgs)*rho0(mgs)/xmas(mgs,lir))
      ENDIF
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qipacw(mgs) = 0.0
      IF ( eipw(mgs) .gt. 0.0 ) THEN
      qipacw(mgs) = 
     >   min(
     >   ((0.25)*pi)*eipw(mgs)*qx(mgs,lc)*cx(mgs,lip)
     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,lc,1))
     >  *(  xdia(mgs,lip,2) 
     >    + 2.0*xdia(mgs,lip,1)*xdia(mgs,lc,1) 
     >    + xdia(mgs,lc,2) )  
     <  , qcmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qipacr(mgs) = 0.0
      IF ( iacr .gt. 0 .and. erip(mgs) .gt. 0.0 
     :     .and. ipconc .lt. 3 ) THEN
      qipacr(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*erip(mgs)*qx(mgs,lr)*cx(mgs,lip)
     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*xdia(mgs,lr,2)
     >    + 2.0*gf5*xdia(mgs,lip,1)*xdia(mgs,lr,1) 
     >    + gf4*xdia(mgs,lip,2) ) 
     <  , qrmxd(mgs))
      ENDIF
      end do
c
c
c      do mgs = 1,ngscnt
c      qiracr(mgs) = 0.0
c      IF ( erir(mgs) .gt. 0.0 ) THEN
c      qiracr(mgs) = 
c     >   min(
c     >   ((0.25)*pi/gf4)*erir(mgs)*qx(mgs,lr)*cx(mgs,lir)
c     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,lr,1))
c     >  *(  gf6*gf1*xdia(mgs,lr,2)
c     >    + 2.0*gf5*gf2*xdia(mgs,lir,1)*xdia(mgs,lr,1) 
c     >    + gf4*gf3*xdia(mgs,lir,2) ) 
c     <  , qrmxd(mgs))
c      ENDIF
c      end do
c
c
      do mgs = 1,ngscnt
      qiracs(mgs) = 0.0
      IF ( esir(mgs) .gt. 0.0 ) THEN
      qiracs(mgs) = 
     >   min(
     >   ((0.25)*pi/gf4)*esir(mgs)*cx(mgs,ls)*qx(mgs,lir)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lir,1)*xdia(mgs,ls,1) 
     >    + gf4*gf3*xdia(mgs,lir,2) ) 
     <  , qirmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qiraci(mgs) = 0.0
      ciraci(mgs) = 0.0
      IF ( eiri(mgs) .gt. 0.0 ) THEN
      qiraci(mgs) = 
     >   min(
     >   ((0.25)*pi)*eiri(mgs)*qx(mgs,li)*cx(mgs,lir)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lir,2) 
     >    + 2.0*gf2*xdia(mgs,lir,1)*xdia(mgs,li,1) 
     >    + gf1*xdia(mgs,li,2) )  
     <  , qimxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
      ciraci(mgs) = Min(cimxd(mgs),qiraci(mgs)*rho0(mgs)/xmas(mgs,li))
      ENDIF
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qiacir(mgs) = 0.0
      ciacir(mgs) = 0.0
      IF ( eiri(mgs) .gt. 0.0 ) THEN
      qiacir(mgs) = 
     >   min(
     >   ((0.25/gf4)*pi)*eiri(mgs)*qx(mgs,lir)*cx(mgs,li)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,li,1))
     >  *(  gf6*xdia(mgs,lir,2) 
     >    + 2.0*gf5*xdia(mgs,lir,1)*xdia(mgs,li,1) 
     >    + gf4*xdia(mgs,li,2) )  
     <  , qirmxd(mgs))
      IF ( ipconc .ge. 1 ) THEN
      ciacir(mgs) = Min(cirmxd(mgs),qiacir(mgs)*rho0(mgs)/xmas(mgs,lir))
      ENDIF
      ENDIF
      end do
c
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: graupel-l collects xxxxx'
c
      do mgs = 1,ngscnt
      qglacw(mgs) = 0.0
      rargl(mgs) = 0.0
      IF ( eglw(mgs) .gt. 0.0 ) THEN
      qglacw(mgs) =
     >   min(
     >   ((0.25)*pi)*eglw(mgs)*qx(mgs,lc)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lgl,2)
     >    + 2.0*gf2*xdia(mgs,lgl,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      IF ( qx(mgs,lgl) .gt. qxmin(lgl) .and. ipelec .gt. 0 ) THEN
        rargl(mgs) =  
     :    qglacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lgl,2)*cx(mgs,lgl))
      ENDIF
      ENDIF  
      end do
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qglacir(mgs) = 0.0
      IF ( eglir(mgs) .gt. 0.0 ) THEN
      qglacir(mgs) =
     >  min(
     >  ((0.25/gf4)*pi)*eglir(mgs)*qx(mgs,lir)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lgl,1)*xdia(mgs,lir,1)
     >    + gf4*gf3*xdia(mgs,lgl,2) )
     <  , qirmxd(mgs))
      ENDIF
      end do
      end if
c
      do mgs = 1,ngscnt
      qglacip(mgs) = 0.0
      IF ( eglip(mgs) .gt. 0.0 ) THEN
      qglacip(mgs) =
     >  min(
     >  ((0.25)*pi)*eglip(mgs)*qx(mgs,lip)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lgl,1)*xdia(mgs,lip,1)
     >    + gf3*xdia(mgs,lgl,2) )
     <  , qipmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qglaci(mgs) = 0.0
      IF ( egli(mgs) .gt. 0.0 ) THEN
      qglaci(mgs) =
     >  min(
     >  ((0.25)*pi)*egli(mgs)*qx(mgs,li)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lgl,2)
     >    + 2.0*gf2*xdia(mgs,lgl,1)*xdia(mgs,li,1)  ! (tm) changed gl2 to gf2
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qglacs(mgs) = 0.0
      IF ( egls(mgs) .gt. 0.0 ) THEN
      qglacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*egls(mgs)*qx(mgs,ls)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lgl,1)
     >    + gf4*gf3*xdia(mgs,lgl,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qglacr(mgs) = 0.0
      cglacr(mgs) = 0.0
      IF ( eglr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgl,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgl,1)*vtxbar(mgs,lr,1) )
       qglacr(mgs) = eglr(mgs)*cx(mgs,lgl)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lgl,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lgl,2))
c     :     ( 0.6366197723675814*xdia(mgs,lr,1)*xdia(mgs,lr,2) + 
c     :       0.889663560367198*xdia(mgs,lgl,1)*xdia(mgs,lr,2) + 
c     :       0.7894150909710455*xdia(mgs,lgl,2)*xdia(mgs,lr,1) )    
        qglacr(mgs) = Min( qglacr(mgs), qrmxd(mgs) )
        cglacr(mgs) = qglacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        cglacr(mgs) = min(cglacr(mgs),crmxd(mgs))
      ELSE
      qglacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*eglr(mgs)*qx(mgs,lr)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lgl,1)
     >    + gf4*gf3*xdia(mgs,lgl,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do
c
c
      if (ndebug .gt. 0 ) print*,'Collection: graupel-m collects xxxxx'
c
      do mgs = 1,ngscnt
      qgmacw(mgs) = 0.0
      rargm(mgs) = 0.0
      IF ( egmw(mgs) .gt. 0.0 ) THEN
      qgmacw(mgs) =
     >   min(
     >   ((0.25)*pi)*egmw(mgs)*qx(mgs,lc)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lgm,2)
     >    + 2.0*gf2*xdia(mgs,lgm,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      IF ( qx(mgs,lgm) .gt. qxmin(lgm) .and. ipelec .gt. 0 ) THEN
        rargm(mgs) =  
     :   qgmacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lgm,2)*cx(mgs,lgm))
      ENDIF
      ENDIF  
      end do
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qgmaci(mgs) = 0.0
      IF ( egmi(mgs) .gt. 0.0 ) THEN
      qgmaci(mgs) =
     >  min(
     >  ((0.25)*pi)*egmi(mgs)*qx(mgs,li)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lgm,2)
     >    + 2.0*gf2*xdia(mgs,lgm,1)*xdia(mgs,li,1)
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do
      end if
c
      do mgs = 1,ngscnt
      qgmacip(mgs) = 0.0
      IF ( egmip(mgs) .gt. 0.0 ) THEN
      qgmacip(mgs) =
     >  min(
     >  ((0.25)*pi)*egmip(mgs)*qx(mgs,lip)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lgm,1)*xdia(mgs,lip,1)
     >    + gf3*xdia(mgs,lgm,2) )
     <  , qipmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qgmacir(mgs) = 0.0
      IF ( egmir(mgs) .gt. 0.0 ) THEN
      qgmacir(mgs) =
     >  min(
     >  ((0.25)*pi/gf4)*egmir(mgs)*qx(mgs,lir)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lgm,1)*xdia(mgs,lir,1)
     >    + gf4*gf3*xdia(mgs,lgm,2) )
     <  , qirmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qgmacs(mgs) = 0.0
      IF ( egms(mgs) .gt. 0.0 ) THEN
      qgmacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*egms(mgs)*qx(mgs,ls)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lgm,1)
     >    + gf4*gf3*xdia(mgs,lgm,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qgmacr(mgs) = 0.0
      cgmacr(mgs) = 0.0
      IF ( egmr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgm,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgm,1)*vtxbar(mgs,lr,1) )
        qgmacr(mgs) = egmr(mgs)*cx(mgs,lgm)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lgm,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lgm,2))
        qgmacr(mgs) = Min( qgmacr(mgs), qrmxd(mgs) )
        cgmacr(mgs) = qgmacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        cgmacr(mgs) = min(cgmacr(mgs),crmxd(mgs))
      ELSE
      qgmacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*egmr(mgs)*qx(mgs,lr)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lgm,1)
     >    + gf4*gf3*xdia(mgs,lgm,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: graupel-h collects xxxxx'
c
      do mgs = 1,ngscnt
      qghacw(mgs) = 0.0
      rargh(mgs) = 0.0
      IF ( eghw(mgs) .gt. 0.0 ) THEN
      qghacw(mgs) =
     >   min(
     >   ((0.25)*pi)*eghw(mgs)*qx(mgs,lc)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lgh,2)
     >    + 2.0*gf2*xdia(mgs,lgh,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      IF ( qx(mgs,lgh) .gt. qxmin(lgh) .and. ipelec .gt. 0 ) THEN
        rargh(mgs) =  
     :   qghacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lgh,2)*cx(mgs,lgh))
      ENDIF
      ENDIF  
      end do
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qghaci(mgs) = 0.0
      IF ( eghi(mgs) .gt. 0.0 ) THEN
      qghaci(mgs) =
     >  min(
     >  ((0.25)*pi)*eghi(mgs)*qx(mgs,li)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lgh,2)
     >    + 2.0*gf2*xdia(mgs,lgh,1)*xdia(mgs,li,1)
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do
      end if
c
      do mgs = 1,ngscnt
      qghacir(mgs) = 0.0
      IF ( eghir(mgs) .gt. 0.0 ) THEN
      qghacir(mgs) =
     >  min(
     >  ((0.25)*pi/gf4)*eghir(mgs)*qx(mgs,lir)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lgh,1)*xdia(mgs,lir,1)
     >    + gf4*gf3*xdia(mgs,lgh,2) )
     <  , qirmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qghacip(mgs) = 0.0
      IF ( eghip(mgs) .gt. 0.0 ) THEN
      qghacip(mgs) =
     >  min(
     >  ((0.25)*pi)*eghip(mgs)*qx(mgs,lip)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lgh,1)*xdia(mgs,lip,1)
     >    + gf3*xdia(mgs,lgh,2) )
     <  , qipmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qghaci(mgs) = 0.0
      IF ( eghi(mgs) .gt. 0.0 ) THEN
      qghaci(mgs) =
     >  min(
     >  ((0.25)*pi)*eghi(mgs)*qx(mgs,li)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lgh,2)
     >    + 2.0*gf2*xdia(mgs,lgh,1)*xdia(mgs,li,1)
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qghacs(mgs) = 0.0
      IF ( eghs(mgs) .gt. 0.0 ) THEN
      qghacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*eghs(mgs)*qx(mgs,ls)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lgh,1)
     >    + gf4*gf3*xdia(mgs,lgh,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qghacr(mgs) = 0.0
      cghacr(mgs) = 0.0
      IF ( eghr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lgh,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lgh,1)*vtxbar(mgs,lr,1) )
        qghacr(mgs) = eghr(mgs)*cx(mgs,lgh)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lgh,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lgh,2))
        qghacr(mgs) = Min( qghacr(mgs), qrmxd(mgs) )
        cghacr(mgs) = qghacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        cghacr(mgs) = min(cghacr(mgs),crmxd(mgs))
      ELSE
      qghacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*eghr(mgs)*qx(mgs,lr)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lgh,1)
     >    + gf4*gf3*xdia(mgs,lgh,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: frozen rw collects xxxxx'
c
      do mgs = 1,ngscnt
      qfacw(mgs) = 0.0
      rarfw(mgs) = 0.0
      IF ( efw(mgs) .gt. 0.0 ) THEN
      qfacw(mgs) =
     >   min(
     >   ((0.25)*pi)*efw(mgs)*qx(mgs,lc)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lf,2)
     >    + 2.0*gf2*xdia(mgs,lf,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      IF ( qx(mgs,lf) .gt. qxmin(lf) .and. ipelec .gt. 0 ) THEN
        rarfw(mgs) =  
     :    qfacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lf,2)*cx(mgs,lf))
      ENDIF
      ENDIF  
      end do   
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qfaci(mgs) = 0.0
      IF ( efi(mgs) .gt. 0.0 ) THEN
      qfaci(mgs) =
     >  min(
     >  ((0.25)*pi)*efi(mgs)*qx(mgs,li)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lf,2)
     >    + 2.0*gf2*xdia(mgs,lf,1)*xdia(mgs,li,1)
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do   
      end if 
c
      do mgs = 1,ngscnt
      qfacip(mgs) = 0.0
      IF ( efip(mgs) .gt. 0.0 ) THEN
      qfacip(mgs) =
     >  min(
     >  ((0.25)*pi/gf4)*efip(mgs)*qx(mgs,lip)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lf,1)*xdia(mgs,lip,1)
     >    + gf3*xdia(mgs,lf,2) )
     <  , qipmxd(mgs))
      ENDIF
      end do   
c
c
      do mgs = 1,ngscnt
      qfacir(mgs) = 0.0
      IF ( efir(mgs) .gt. 0.0 ) THEN
      qfacir(mgs) =
     >  min(
     >  ((0.25/gf4)*pi)*efir(mgs)*qx(mgs,lir)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lf,1)*xdia(mgs,lir,1)
     >    + gf4*gf3*xdia(mgs,lf,2) )
     <  , qirmxd(mgs))
      ENDIF
      end do   
c
      do mgs = 1,ngscnt
      qfacs(mgs) = 0.0
      IF ( efs(mgs) .gt. 0.0 ) THEN
      qfacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*efs(mgs)*qx(mgs,ls)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lf,1)
     >    + gf4*gf3*xdia(mgs,lf,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do   
c
      do mgs = 1,ngscnt
      qfacr(mgs) = 0.0
      cfacr(mgs) = 0.0
      IF ( efr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lf,1)*vtxbar(mgs,lr,1) )
        qfacr(mgs) = efr(mgs)*cx(mgs,lf)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lf,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lf,2))
        qfacr(mgs) = Min( qfacr(mgs), qrmxd(mgs) )
        cfacr(mgs) = qfacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        cfacr(mgs) = min(cfacr(mgs),crmxd(mgs))
      ELSE
      qfacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*efr(mgs)*qx(mgs,lr)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lf,1)
     >    + gf4*gf3*xdia(mgs,lf,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do   
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: hail collects xxxxx'
c
      do mgs = 1,ngscnt
      qhacw(mgs) = 0.0
      rarhw(mgs) = 0.0
      IF ( ehw(mgs) .gt. 0.0 ) THEN
      qhacw(mgs) = 
     >   min(
     >   ((0.25)*pi)*ehw(mgs)*qx(mgs,lc)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lh,2) 
     >    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,lc,1) 
     >    + gf1*xdia(mgs,lc,2) )  
     <  , qcmxd(mgs))
      IF ( qx(mgs,lh) .gt. qxmin(lh) .and. ipelec .gt. 0 ) THEN
        rarhw(mgs) =  
     :   qhacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lh,2)*cx(mgs,lh))
      ENDIF
      ENDIF  
      end do   
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qhaci(mgs) = 0.0
      IF ( ehi(mgs) .gt. 0.0 ) THEN
      qhaci(mgs) = 
     >  min(
     >  ((0.25)*pi)*ehi(mgs)*qx(mgs,li)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lh,2) 
     >    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,li,1) 
     >    + gf1*xdia(mgs,li,2) )  
     <  , qimxd(mgs))
      ENDIF
      end do   
      end if   
c
      do mgs = 1,ngscnt
      qhacir(mgs) = 0.0
      IF ( ehir(mgs) .gt. 0.0 ) THEN
      qhacir(mgs) = 
     >  min(
     >  ((0.25/gf4)*pi)*ehir(mgs)*qx(mgs,lir)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2) 
     >    + 2.0*gf5*gf2*xdia(mgs,lh,1)*xdia(mgs,lir,1) 
     >    + gf4*gf3*xdia(mgs,lh,2) )  
     <  , qirmxd(mgs))
      ENDIF
      end do   
c
      do mgs = 1,ngscnt
      qhacip(mgs) = 0.0
      IF ( ehip(mgs) .gt. 0.0 ) THEN
      qhacip(mgs) = 
     >  min(
     >  ((0.25)*pi)*ehip(mgs)*qx(mgs,lip)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2) 
     >    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,lip,1) 
     >    + gf3*xdia(mgs,lh,2) )  
     <  , qipmxd(mgs))
      ENDIF
      end do   
c
c
      do mgs = 1,ngscnt
      qhacs(mgs) = 0.0
      IF ( ehs(mgs) .gt. 0.0 ) THEN
      qhacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*ehs(mgs)*qx(mgs,ls)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)
     >    + gf4*gf3*xdia(mgs,lh,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do   
c
      do mgs = 1,ngscnt
      qhacr(mgs) = 0.0
      chacr(mgs) = 0.0
      IF ( ehr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lr,1) )
       qhacr(mgs) = ehr(mgs)*cx(mgs,lh)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lh,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lh,2))
        qhacr(mgs) = Min( qhacr(mgs), qrmxd(mgs) )
        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        chacr(mgs) = min(chacr(mgs),crmxd(mgs))
      ELSE
      qhacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*ehr(mgs)*qx(mgs,lr)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lh,1)
     >    + gf4*gf3*xdia(mgs,lh,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do   
c
c
      if (ndebug .gt. 0 ) print*,'Collection: large hail collects xxxxx'
c
      do mgs = 1,ngscnt
      qhlacw(mgs) = 0.0
      rarhl(mgs) = 0.0
      IF ( ehlw(mgs) .gt. 0.0 ) THEN
      qhlacw(mgs) =
     >   min(
     >   ((0.25)*pi)*ehlw(mgs)*qx(mgs,lc)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lhl,2)
     >    + 2.0*gf2*xdia(mgs,lhl,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      IF ( qx(mgs,lhl) .gt. qxmin(lh) .and. ipelec .gt. 0 ) THEN
        rarhl(mgs) =  
     :    qhlacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lhl,2)*cx(mgs,lhl))
      ENDIF
      ENDIF  
      end do
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qhlaci(mgs) = 0.0
      IF ( ehli(mgs) .gt. 0.0 ) THEN
      qhlaci(mgs) =
     >  min(
     >  ((0.25)*pi)*ehli(mgs)*qx(mgs,li)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))
     >  *(  gf3*xdia(mgs,lhl,2)
     >    + 2.0*gf2*xdia(mgs,lhl,1)*xdia(mgs,li,1)
     >    + gf1*xdia(mgs,li,2) )
     <  , qimxd(mgs))
      ENDIF
      end do
      end if
c
      do mgs = 1,ngscnt
      qhlacir(mgs) = 0.0
      IF ( ehlir(mgs) .gt. 0.0 ) THEN
      qhlacir(mgs) =
     >  min(
     >  ((0.25/gf4)*pi)*ehlir(mgs)*qx(mgs,lir)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lhl,1)*xdia(mgs,lir,1)
     >    + gf4*gf3*xdia(mgs,lhl,2) )
     <  , qirmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qhlacip(mgs) = 0.0
      IF ( ehlip(mgs) .gt. 0.0 ) THEN
      qhlacip(mgs) =
     >  min(
     >  ((0.25)*pi)*ehlip(mgs)*qx(mgs,lip)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lhl,1)*xdia(mgs,lip,1)
     >    + gf3*xdia(mgs,lhl,2) )
     <  , qipmxd(mgs))
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qhlacs(mgs) = 0.0
      IF ( ehls(mgs) .gt. 0.0 ) THEN
      qhlacs(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*ehls(mgs)*qx(mgs,ls)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lhl,1)
     >    + gf4*gf3*xdia(mgs,lhl,2) )
     >  , qsmxd(mgs))
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qhlacr(mgs) = 0.0
      chlacr(mgs) = 0.0
      IF ( ehlr(mgs) .gt. 0.0 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))**2 + 
     :            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,lr,1) )
       qhlacr(mgs) = ehlr(mgs)*cx(mgs,lhl)*vt*
     :     qx(mgs,lr)*0.25*pi*
     :      (3.02787*xdia(mgs,lr,2) + 
     :       3.30669*xdia(mgs,lhl,1)*xdia(mgs,lr,1) + 
     :       2.*xdia(mgs,lhl,2))
        qhlacr(mgs) = Min( qhlacr(mgs), qrmxd(mgs) )
        chlacr(mgs) = qhlacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        chlacr(mgs) = min(chlacr(mgs),crmxd(mgs))
      ELSE
      qhlacr(mgs) =
     >   min(
     >   ((0.25)*pi/gf4)*ehlr(mgs)*qx(mgs,lr)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lhl,1)
     >    + gf4*gf3*xdia(mgs,lhl,2) )
     <  , qrmxd(mgs))
      ENDIF
      ENDIF
      end do                              
c
c
      if (ndebug .gt. 0 ) print*,'Collection: Cloud collects xxxxx'
c
c 
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qwaci(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 ) THEN
      qwaci(mgs) =
     >   min(
     >   (0.25)*pi*eiw(mgs)*qx(mgs,li)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,li,1)-vtxbar(mgs,lc,1))
     >  *(  xdia(mgs,li,2)
     >    + 2.0*xdia(mgs,li,1)*xdia(mgs,lc,1)
     >    + xdia(mgs,lc,2) )
     <  , qimxd(mgs))
      qwaci(mgs) = qwaci(mgs)*il5(mgs)
      ENDIF
      end do
      endif
c
      do mgs = 1,ngscnt
      qwacir(mgs) = 0.0
      IF ( eirw(mgs) .gt. 0.0 ) THEN
      qwacir(mgs) =
     >   min(
     >   (0.25/gf4)*pi*eirw(mgs)*qx(mgs,lir)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,lc,1))
     >  *(  gf6*xdia(mgs,lir,2)
     >    + 2.0*gf5*xdia(mgs,lir,1)*xdia(mgs,lc,1)
     >    + gf4*xdia(mgs,lc,2) )
     <  , qirmxd(mgs))
      qwacir(mgs) = qwacir(mgs)*il5(mgs)
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qwacip(mgs) = 0.0
      IF ( eipw(mgs) .gt. 0.0 ) THEN
      qwacip(mgs) =
     >   min(
     >   (0.25)*pi*eipw(mgs)*qx(mgs,lip)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,lc,1))
     >  *(  xdia(mgs,lip,2)
     >    + 2.0*xdia(mgs,lip,1)*xdia(mgs,lc,1)
     >    + xdia(mgs,lc,2) )
     <  , qipmxd(mgs))
      qwacip(mgs) = qwacip(mgs)*il5(mgs)
      ENDIF
      end do
c
c
      do mgs = 1,ngscnt
      qwacs(mgs) =  0.0
      IF ( esw(mgs) .gt. 0.0 ) THEN
      qwacs(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*esw(mgs)*qx(mgs,ls)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,ls,1))
     >  *(  gf6*xdia(mgs,ls,2) 
     >    + 2.0*gf5*xdia(mgs,ls,1)*xdia(mgs,lc,1) 
     >    + gf4*xdia(mgs,lc,2) )  
     <  , fconv*qsmxd(mgs)) 
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qwacgl(mgs) =  0.0
      IF ( eglw(mgs) .gt. 0.0 ) THEN
      qwacgl(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*eglw(mgs)*qx(mgs,lgl)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lgl,1))
     >  *(  gf6*xdia(mgs,lgl,2) 
     >    + 2.0*gf5*xdia(mgs,lgl,1)*xdia(mgs,lc,1) 
     >    + gf4*xdia(mgs,lc,2) )  
     <  , fconv*qglmxd(mgs)) 
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qwacgm(mgs) = 0.0
      IF ( egmw(mgs) .gt. 0.0 ) THEN
      qwacgm(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*egmw(mgs)*qx(mgs,lgm)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lgm,1))
     >  *(  gf6*xdia(mgs,lgm,2) 
     >    + 2.0*gf5*xdia(mgs,lgm,1)*xdia(mgs,lc,1) 
     >    + gf4*xdia(mgs,lc,2) )  
     <  , fconv*qgmmxd(mgs)) 
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qwacgh(mgs) =  0.0
      IF ( eghw(mgs) .gt. 0.0 ) THEN
      qwacgh(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*eghw(mgs)*qx(mgs,lgh)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lgh,1))
     >  *(  gf6*xdia(mgs,lgh,2) 
     >    + 2.0*gf5*xdia(mgs,lgh,1)*xdia(mgs,lc,1) 
     >    + gf4*xdia(mgs,lc,2) )  
     <  , fconv*qghmxd(mgs)) 
      ENDIF
      end do
c
      do mgs = 1,ngscnt
      qwacf(mgs) =  0.0
      IF ( efw(mgs) .gt. 0.0 ) THEN
      qwacf(mgs) = 
     >  min(     
     >    ((0.25/gf4)*pi)*efw(mgs)*qx(mgs,lf)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lf,1))
     >  *(  gf6*xdia(mgs,lf,2) 
     >    + 2.0*gf5*xdia(mgs,lf,1)*xdia(mgs,lc,1) 
     >    + gf4*xdia(mgs,lc,2) )  
     <  , fconv*qfmxd(mgs)) 
      ENDIF
      end do
c
      if (ndebug .gt. 0 ) print*,'Collection: Cloud collects xxxxx'
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx1'      
c
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx2'      
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qiacw(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 ) THEN
      qiacw(mgs) =
     >   min(
     >   (0.25)*pi*eiw(mgs)*qx(mgs,lc)*cx(mgs,li)
     >  *abs(vtxbar(mgs,li,1)-vtxbar(mgs,lc,1))
     >  *(  xdia(mgs,li,2)
     >    + 2.0*xdia(mgs,li,1)*xdia(mgs,lc,1)
     >    + xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      qiacw(mgs) = qiacw(mgs)*il5(mgs)
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx3'      
c
      do mgs = 1,ngscnt
      qiracw(mgs) = 0.0
      IF ( eirw(mgs) .gt. 0.0 ) THEN
      qiracw(mgs) =
     >   min(
     >   (0.25)*pi*eirw(mgs)*qx(mgs,lc)*cx(mgs,lir)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,lc,1))
     >  *(  gf3*xdia(mgs,lir,2)
     >    + gf2*2.0*xdia(mgs,lir,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
     <  , qcmxd(mgs))
      qiracw(mgs) = qiracw(mgs)*il5(mgs)
      ENDIF
      end do
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx4'      
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx7'      
c
      do mgs = 1,ngscnt
      qiracr(mgs) = 0.0
      IF ( iacr .eq. 1 .and. erir(mgs) .gt. 0.0 .and.
     :     ipconc .lt. 3 ) THEN
      qiracr(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*erir(mgs)*cx(mgs,lir)*qx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lir,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2) 
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lir,1) 
     >    + gf4*gf3*xdia(mgs,lir,2) )  
     <  , qrmxd(mgs)) 
      if ( temg(mgs) .gt. 268.15 ) then
      qiracr(mgs) = 0.0
      end if
      ENDIF
      end do
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx8'      
c
      do mgs = 1,ngscnt
      qiacr(mgs) = 0.0
      ciacr(mgs) = 0.0
      ciacrf(mgs) = 0.0
      IF ( iacr .eq. 1 .and. eri(mgs) .gt. 0.0 
     :     .and. temg(mgs) .le. 270.15 ) THEN
      IF ( ipconc .ge. 3 ) THEN
       ni = 0.0
         IF ( xdia(mgs,li,1) .ge. 10.e-6 ) THEN
          ni = ni + cx(mgs,li)*Exp(- (40.e-6/xdia(mgs,li,1))**3 )
         ENDIF
         IF ( xdia(mgs,lip,1) .ge. 10.e-6 ) THEN
          ni = ni + cx(mgs,lip)*Exp(- (40.e-6/xdia(mgs,lip,1))**3 )
         ENDIF
         IF ( xdia(mgs,lir,1) .ge. 10.e-6 ) THEN
          ni = ni + cx(mgs,lir)*exp( -40.e-6/xdia(mgs,lir,1) )
         ENDIF
c   Set nr to the number of drops greater than 40 microns.
         arg = 1000.*xdia(mgs,lr,1)
         nr = cx(mgs,lr)*gaml02( arg )
       IF ( ni .gt. 0.0 .and. nr .gt. 0.0 ) THEN
       d0 = xdia(mgs,lr,1)
       qiacr(mgs) = xdn(mgs,lr)*rhoinv(mgs)*
     :     (0.217239*(0.522295*(d0**5) + 
     :      49711.81*(d0**6) - 
     :      1.673016e7*(d0**7)+ 
     :      2.404471e9*(d0**8) - 
     :      1.22872e11*(d0**9))*ni*nr)
      qiacr(mgs) = Min( qrmxd(mgs), qiacr(mgs) )
      ciacr(mgs) =
     :   (0.217239*(0.2301947*(d0**2) + 
     :      15823.76*(d0**3) - 
     :      4.167685e6*(d0**4) + 
     :      4.920215e8*(d0**5) - 
     :      2.133344e10*(d0**6))*ni*nr)
      ciacr(mgs) = Min( crmxd(mgs), ciacr(mgs) )
c      ciacr(mgs) = qiacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
      ciacrf(mgs) = rz*ciacr(mgs)
       ENDIF
      ELSE
      qiacr(mgs) = 
     >  min(     
     >   ((0.25/gf4)*pi)*eri(mgs)*cx(mgs,li)*qx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))
     >  *(  gf6*gf1*xdia(mgs,lr,2) 
     >    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1) 
     >    + gf4*gf3*xdia(mgs,li,2) )  
     <  , qrmxd(mgs)) 
      ENDIF
c      if ( temg(mgs) .gt. 268.15 ) then
c      qiacr(mgs) = 0.0
c      ciacr(mgs) = 0.0
c      end if
      ENDIF
      end do
c
c
c
c
      if (ndebug .gt. 0 ) print*,'Collection: cloud ice collects xxxx11'      
c
c
c  Concentration rates.....cxacy
c
c
c snow
c
      if (ndebug .gt. 0 ) print*,'conc 1'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      csacw(mgs) = 0.0
      IF ( esw(mgs) .gt. 0.0 ) THEN
      csacw(mgs) =
     >   ((0.25)*pi)*esw(mgs)*cx(mgs,lc)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      csacw(mgs) = min(csacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 2'
      if ( ipconc .ge. 4 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cwacs(mgs) = 0.0
      IF ( esw(mgs) .gt. 0.0 ) THEN
      cwacs(mgs) =
     >   ((0.25)*pi)*esw(mgs)*cx(mgs,lc)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      cwacs(mgs) = min(cwacs(mgs), fconv*csmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 3'
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cipacs(mgs) = 0.0
      cipacs(mgs) =
     >   ((0.25)*pi)*esip(mgs)*cx(mgs,lip)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      cipacs(mgs) = min(cipacs(mgs),csmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 4'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      csacip(mgs) = 0.0
      IF ( esip(mgs) .gt. 0.0 ) THEN
      csacip(mgs) =
     >   ((0.25)*pi)*esip(mgs)*cx(mgs,lip)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      csacip(mgs) = min(csacip(mgs),cipmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 5'
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      ciacs(mgs) = 0.0
      ciacs(mgs) =
     >   ((0.25)*pi)*esi(mgs)*cx(mgs,li)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      ciacs(mgs) = min(ciacs(mgs),csmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 6'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      csaci(mgs) = 0.0
      IF ( esi(mgs) .gt. 0.0 ) THEN
       vt = Sqrt((vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))**2 + 
     :            0.04*vtxbar(mgs,ls,1)*vtxbar(mgs,li,1) )
      csaci(mgs) =
     >   ((0.25)*pi)*esi(mgs)*cx(mgs,li)*cx(mgs,ls)
     :  * vt
c     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,ls,1)
     >    + gf3*xdia(mgs,ls,2) )
      csaci(mgs) = min(csaci(mgs),cimxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 7'
c snow aggregation here
      if ( ipconc .ge. 4 .or. ipelec .ge. 100 ) then ! 
      do mgs = 1,ngscnt
      csacs(mgs) = 0.0
      IF ( ess(mgs) .gt. 0.0 ) THEN
c      csacs(mgs) = -a2*eps*ess(mgs)*cx(mgs,ls)**2*xv(mgs,ls)
      csacs(mgs) = rvt*aa2*ess(mgs)*cx(mgs,ls)**2*xv(mgs,ls)
c      csacs(mgs) =
c     >   ((0.25)*pi)*ess(mgs)*cx(mgs,ls)*cx(mgs,ls)
c     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,ls,1))
c     >  *(  gf3*gf1*xdia(mgs,ls,2)
c     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,ls,1)
c     >    + gf1*gf3*xdia(mgs,ls,2) )
      csacs(mgs) = min(csacs(mgs),csmxd(mgs))
      ENDIF
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 8'
c      if ( ipconc .ge. 3 .or. ipelec .ge. 1 ) then
c      do mgs = 1,ngscnt
c      csacr(mgs) = 0.0
c      IF ( esr(mgs) .gt. 0.0 ) THEN
c      csacr(mgs) =
c     >   ((0.25)*pi)*esr(mgs)*cx(mgs,lr)*cx(mgs,ls)
c     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lr,1))
c     >  *(  gf3*gf1*xdia(mgs,lr,2)
c     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,ls,1)
c     >    + gf1*gf3*xdia(mgs,ls,2) )
c      csacr(mgs) = min(csacr(mgs),crmxd(mgs))
c      ENDIF
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 9'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cracs(mgs) = 0.0
      IF ( esr(mgs) .gt. 0.0 ) THEN
      cracs(mgs) =
     >   ((0.25)*pi)*esr(mgs)*cx(mgs,lr)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,ls,1)
     >    + gf1*gf3*xdia(mgs,ls,2) )
      cracs(mgs) = min(cracs(mgs),fconv*csmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 10'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cwaci(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 ) THEN
      cwaci(mgs) =
     >   ((0.25)*pi)*eiw(mgs)*cx(mgs,li)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,li,1))
     >  *(  xdia(mgs,li,2)
     >    + 2.0*xdia(mgs,li,1)*xdia(mgs,lc,1)
     >    + xdia(mgs,lc,2) )
      cwaci(mgs) = min(cwaci(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 11'
      if ( ipconc .ge. 2 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      ciacw(mgs) = 0.0
      IF ( eiw(mgs) .gt. 0.0 ) THEN

        ciacw(mgs) = qiacw(mgs)*rho0(mgs)/xmas(mgs,lc)

c      ciacw(mgs) =
c     >   ((0.25)*pi)*eiw(mgs)*cx(mgs,lc)*cx(mgs,li)
c     >  *abs(vtxbar(mgs,li,1)-vtxbar(mgs,lc,1))
c     >  *(  xdia(mgs,lc,2)
c     >    + 2.0*xdia(mgs,lc,1)*xdia(mgs,li,1)
c     >    + xdia(mgs,li,2) )
c      ciacw(mgs) = min(ciacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 12'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cwacip(mgs) = 0.0
      IF ( eipw(mgs) .gt. 0.0 ) THEN
      cwacip(mgs) =
     >   ((0.25)*pi)*eipw(mgs)*cx(mgs,lip)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lip,1))
     >  *(  xdia(mgs,lip,2)
     >    + 2.0*xdia(mgs,lip,1)*xdia(mgs,lc,1)
     >    + xdia(mgs,lc,2) )
      cwacip(mgs) = min(cwacip(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 13'
      if ( ipconc .ge. 2 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cipacw(mgs) = 0.0
      IF ( eipw(mgs) .gt. 0.0 ) THEN
      cipacw(mgs) =
     >   ((0.25)*pi)*eipw(mgs)*cx(mgs,lc)*cx(mgs,lip)
     >  *abs(vtxbar(mgs,lip,1)-vtxbar(mgs,lc,1))
     >  *(  xdia(mgs,lc,2)
     >    + 2.0*xdia(mgs,lc,1)*xdia(mgs,lip,1)
     >    + xdia(mgs,lip,2) )
      cipacw(mgs) = min(cipacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
c
      if (ndebug .gt. 0 ) print*,'conc 14'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cracip(mgs) = 0.0
      IF ( erip(mgs) .gt. 0.0 ) THEN
      cracip(mgs) =
     >   ((0.25)*pi)*erip(mgs)*cx(mgs,lip)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lr,1)
     >    + gf3*xdia(mgs,lr,2) )
      cracip(mgs) = min(cracip(mgs),cipmxd(mgs))
      ENDIF
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 15'
c      if ( ipconc .ge. 3 .or. ipelec .ge. 9 ) then
c      do mgs = 1,ngscnt
c      cipacr(mgs) = 0.0
c      IF ( iacr .ge. 1 .and. erip(mgs) .gt. 0.0 ) THEN
c      cipacr(mgs) =
c     >   ((0.25)*pi)*erip(mgs)*cx(mgs,lip)*cx(mgs,lr)
c     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lip,1))
c     >  *(  gf1*xdia(mgs,lip,2)
c     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lr,1)
c     >    + gf3*xdia(mgs,lr,2) )
c      cipacr(mgs) = min(cipacr(mgs),crmxd(mgs))
c      ENDIF
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 16'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      craci(mgs) = 0.0
      IF ( eri(mgs) .gt. 0.0 ) THEN
      craci(mgs) =
     >   ((0.25)*pi)*eri(mgs)*cx(mgs,li)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lr,1)
     >    + gf3*xdia(mgs,lr,2) )
      craci(mgs) = min(craci(mgs),cimxd(mgs))
      ENDIF
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 17'
c      if ( ipconc .ge. 3 .or. ipelec .ge. 9 ) then
c      do mgs = 1,ngscnt
c      ciacr(mgs) = 0.0
c      IF ( iacr .ge. 1 .and. eri(mgs) .gt. 0.0 ) THEN
c      ciacr(mgs) =
c     >   ((0.25)*pi)*eri(mgs)*cx(mgs,li)*cx(mgs,lr)
c     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))
c     >  *(  gf1*xdia(mgs,lc,2)
c     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lr,1)
c     >    + gf3*xdia(mgs,lr,2) )
c      ciacr(mgs) = min(ciacr(mgs),crmxd(mgs))
c      ENDIF
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 18'
      if ( ipconc .ge. 2 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
       cracw(mgs) = 0.0
       cracr(mgs) = 0.0
       IF ( ipconc .lt. 3 ) THEN
        IF ( erw(mgs) .gt. 0.0 ) THEN
        cracw(mgs) =
     >   ((0.25)*pi)*erw(mgs)*cx(mgs,lc)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lr,1)
     >    + gf3*xdia(mgs,lr,2) )
        ENDIF
       ELSEIF ( ipconc .ge. 3 .and. 
     :     ( 0.5*xdia(mgs,lr,1) .gt. rh(mgs)  ) ) THEN ! .or. cx(mgs,lr) .gt. nh(mgs))) THEN
        IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( 0.5*xdia(mgs,lr,1) .gt. 50.e-6 ) THEN
          cracw(mgs) = aa2*cx(mgs,lr)*cx(mgs,lc)*(xv(mgs,lc) + xv(mgs,lr))
        ELSE
          cracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*
     :        ((cnu + 2)*xv(mgs,lc)**2/(cnu + 1) + 
     :         (rnu + 2)*xv(mgs,lr)**2/(rnu + 1))
        ENDIF
        ENDIF
        
c Rain self collection (cracr) and break-up (factor of ec0)
c       
        IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        rwrad = 0.5*xdia(mgs,lr,1)
        IF ( rwrad .gt. 1.0e-3 ) THEN
          ec0 = 0.0
          cracr(mgs) = 0.0
        ELSE
          IF ( rwrad .lt. 3.0e-4 ) THEN
            ec0 = 1.0
          ELSE
            ec0 = Exp(-50.0*(100.0*rwrad - 0.03))
          ENDIF
        
          IF ( rwrad .ge. 50.e-6 ) THEN
            cracr(mgs) = ec0*aa2*cx(mgs,lr)**2*xv(mgs,lr)
          ELSE
            cracr(mgs) = ec0*aa1*cx(mgs,lr)**2*xv(mgs,lr)**2*
     :                   (rnu + 2)/(rnu + 1)
          ENDIF
          cracr(mgs) = Min(cracr(mgs),crmxd(mgs))
        ENDIF
        ENDIF
        
        
        
       ENDIF
      cracw(mgs) = min(cracw(mgs),ccmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 19'
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cwacr(mgs) = 0.0
      cwacr(mgs) =
     >   ((0.25)*pi)*erw(mgs)*cx(mgs,lc)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lr,1)
     >    + gf3*xdia(mgs,lr,2) )
      cwacr(mgs) = min(cwacr(mgs),crmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 20'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cwacir(mgs) = 0.0
      IF ( eirw(mgs) .gt. 0.0 ) THEN
      cwacir(mgs) =
     >   ((0.25)*pi)*eirw(mgs)*cx(mgs,lir)*cx(mgs,lc)
     >  *abs(vtxbar(mgs,lc,1)-vtxbar(mgs,lir,1))
     >  *(  gf3*xdia(mgs,lir,2)
     >    + 2*gf2*xdia(mgs,lir,1)*xdia(mgs,lc,1)
     >    + gf1*xdia(mgs,lc,2) )
      cwacir(mgs) = min(cwacir(mgs),crmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 21'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      ciracw(mgs) = 0.0
      IF ( eirw(mgs) .gt. 0.0 ) THEN
      ciracw(mgs) =
     >   ((0.25)*pi)*eirw(mgs)*cx(mgs,lc)*cx(mgs,lir)
     >  *abs(vtxbar(mgs,lir,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lir,1)
     >    + gf3*xdia(mgs,lir,2) )
      ciracw(mgs) = min(ciracw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22a'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cracir(mgs) = 0.0
      IF ( erir(mgs) .gt. 0.0 ) THEN
      cracir(mgs) =
     >   ((0.25)*pi)*erir(mgs)*cx(mgs,lir)*cx(mgs,lr)
     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,lr,1)
     >    + gf3*gf1*xdia(mgs,lr,2) )
      cracir(mgs) = min(cracir(mgs),cirmxd(mgs))
      ENDIF
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 22b'
c      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
c      do mgs = 1,ngscnt
c      ciracr(mgs) = 0.0
c      IF ( iacr .ge. 1 .and. erir(mgs) .gt. 0.0 ) THEN
c      ciracr(mgs) =
c     >   ((0.25)*pi)*erir(mgs)*cx(mgs,lir)*cx(mgs,lr)
c     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lir,1))
c     >  *(  gf3*gf1*xdia(mgs,lir,2)
c     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,lr,1)
c     >    + gf1*gf3*xdia(mgs,lr,2) )
c      ciracr(mgs) = min(ciracr(mgs),crmxd(mgs))
c      ENDIF
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22c'
      if ( ipconc .ge. 1 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      ciracs(mgs) = 0.0
      IF ( esir(mgs) .gt. 0.0 ) THEN
      ciracs(mgs) =
     >   ((0.25)*pi)*esir(mgs)*cx(mgs,lir)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lir,1))
     >  *(  gf3*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,ls,1)
     >    + gf1*gf3*xdia(mgs,ls,2) )
      ciracs(mgs) = min(ciracs(mgs),csmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22d'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      csacir(mgs) = 0.0
      IF ( esir(mgs) .gt. 0.0 ) THEN
      csacir(mgs) =
     >   ((0.25)*pi)*esir(mgs)*cx(mgs,lir)*cx(mgs,ls)
     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lir,1))
     >  *(  gf3*gf1*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,ls,1)
     >    + gf1*gf3*xdia(mgs,ls,2) )
      csacir(mgs) = min(csacir(mgs),cirmxd(mgs))
      ENDIF
      end do
      end if
c                    
c
c
c  Graupel: low density
c
      if (ndebug .gt. 0 ) print*,'conc 22e'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cglacw(mgs) = 0.0
      IF ( eglw(mgs) .gt. 0.0 ) THEN
      cglacw(mgs) =
     >   ((0.25)*pi)*eglw(mgs)*cx(mgs,lc)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgl,1)
     >    + gf3*xdia(mgs,lgl,2) )
      cglacw(mgs) = min(cglacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22f'
      if ( ipconc .ge. 5 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cwacgl(mgs) = 0.0
      cwacgl(mgs) =
     >   ((0.25)*pi)*eglw(mgs)*cx(mgs,lc)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgl,1)
     >    + gf3*xdia(mgs,lgl,2) )
      cwacgl(mgs) = min(cwacgl(mgs), fconv*cglmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22g'
      if ( ipconc .ge. 13 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cglacr(mgs) = 0.0
      IF ( eglr(mgs) .gt. 0.0 ) THEN
      cglacr(mgs) =
     >   ((0.25)*pi)*eglr(mgs)*cx(mgs,lr)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgl,1)
     >    + gf1*gf3*xdia(mgs,lgl,2) )
      cglacr(mgs) = min(cglacr(mgs),crmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22h'
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cracgl(mgs) = 0.0
      IF ( eglr(mgs) .gt. 0.0 ) THEN
      cracgl(mgs) =
     >   ((0.25)*pi)*eglr(mgs)*cx(mgs,lr)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgl,1)
     >    + gf1*gf3*xdia(mgs,lgl,2) )
      cracgl(mgs) = min(cracgl(mgs),fconv*cglmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22i'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cglaci(mgs) = 0.0
      IF ( egli(mgs) .gt. 0.0 ) THEN
      cglaci(mgs) =
     >   ((0.25)*pi)*egli(mgs)*cx(mgs,li)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lgl,1)
     >    + gf3*xdia(mgs,lgl,2) )
      cglaci(mgs) = min(cglaci(mgs),cimxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22i'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cglacip(mgs) = 0.0
      IF ( eglip(mgs) .gt. 0.0 ) THEN
      cglacip(mgs) =
     >   ((0.25)*pi)*eglip(mgs)*cx(mgs,lip)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lgl,1)
     >    + gf3*xdia(mgs,lgl,2) )
      cglacip(mgs) = min(cglacip(mgs),cipmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22j'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cglacir(mgs) = 0.0
      IF ( eglir(mgs) .gt. 0.0 ) THEN
      cglacir(mgs) =
     >   ((0.25)*pi)*eglir(mgs)*cx(mgs,lir)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,lgl,1)
     >    + gf3*gf1*xdia(mgs,lgl,2) )
      cglacir(mgs) = min(cglacir(mgs),cirmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22k'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cglacs(mgs) = 0.0
      IF ( egls(mgs) .gt. 0.0 ) THEN
      cglacs(mgs) =
     >   ((0.25)*pi)*egls(mgs)*cx(mgs,ls)*cx(mgs,lgl)
     >  *abs(vtxbar(mgs,lgl,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lgl,1)
     >    + gf1*gf3*xdia(mgs,lgl,2) )
      cglacs(mgs) = min(cglacs(mgs),csmxd(mgs))
      ENDIF
      end do
      end if
c                    
c
c
c
c  Graupel:  medium density
c
      if (ndebug .gt. 0 ) print*,'conc 22l'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cgmacw(mgs) = 0.0
      IF ( egmw(mgs) .gt. 0.0 ) THEN
      cgmacw(mgs) =
     >   ((0.25)*pi)*egmw(mgs)*cx(mgs,lc)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgm,1)
     >    + gf3*xdia(mgs,lgm,2) )
      cgmacw(mgs) = min(cgmacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22m'
      if ( ipconc .ge. 5 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cwacgm(mgs) = 0.0
      IF ( egmw(mgs) .gt. 0.0 ) THEN
      cwacgm(mgs) =
     >   ((0.25)*pi)*egmw(mgs)*cx(mgs,lc)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgm,1)
     >    + gf3*xdia(mgs,lgm,2) )
      cwacgm(mgs) = min(cwacgm(mgs), fconv*cgmmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22n'
      if ( ipconc .ge. 13 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cgmacr(mgs) = 0.0
      IF ( egmr(mgs) .gt. 0.0 ) THEN
      cgmacr(mgs) =
     >   ((0.25)*pi)*egmr(mgs)*cx(mgs,lr)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgm,1)
     >    + gf1*gf3*xdia(mgs,lgm,2) )
      cgmacr(mgs) = min(cgmacr(mgs),crmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22o'
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cracgm(mgs) = 0.0
      IF ( egmr(mgs) .gt. 0.0 ) THEN
      cracgm(mgs) =
     >   ((0.25)*pi)*egmr(mgs)*cx(mgs,lr)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgm,1)
     >    + gf1*gf3*xdia(mgs,lgm,2) )
      cracgm(mgs) = min(cracgm(mgs),fconv*cgmmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22p'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cgmaci(mgs) = 0.0
      IF ( egmi(mgs) .gt. 0.0 ) THEN
      cgmaci(mgs) =
     >   ((0.25)*pi)*egmi(mgs)*cx(mgs,li)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lgm,1)
     >    + gf3*xdia(mgs,lgm,2) )
      cgmaci(mgs) = min(cgmaci(mgs),cimxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22q'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cgmacip(mgs) = 0.0
      IF ( egmip(mgs) .gt. 0.0 ) THEN
      cgmacip(mgs) =
     >   ((0.25)*pi)*egmip(mgs)*cx(mgs,lip)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lgm,1)
     >    + gf3*xdia(mgs,lgm,2) )
      cgmacip(mgs) = min(cgmacip(mgs),cipmxd(mgs))
      ENDIF
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22r'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cgmacir(mgs) = 0.0
      cgmacir(mgs) =
     >   ((0.25)*pi)*egmir(mgs)*cx(mgs,lir)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,lgm,1)
     >    + gf3*gf1*xdia(mgs,lgm,2) )
      cgmacir(mgs) = min(cgmacir(mgs),cirmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22s'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cgmacs(mgs) = 0.0
      cgmacs(mgs) =
     >   ((0.25)*pi)*egms(mgs)*cx(mgs,ls)*cx(mgs,lgm)
     >  *abs(vtxbar(mgs,lgm,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lgm,1)
     >    + gf1*gf3*xdia(mgs,lgm,2) )
      cgmacs(mgs) = min(cgmacs(mgs),csmxd(mgs))
      end do
      end if
c                    
c
c
c  Graupel:  High density
c
      if (ndebug .gt. 0 ) print*,'conc 22t'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cghacw(mgs) = 0.0
      cghacw(mgs) =
     >   ((0.25)*pi)*eghw(mgs)*cx(mgs,lc)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgh,1)
     >    + gf3*xdia(mgs,lgh,2) )
      cghacw(mgs) = min(cghacw(mgs),ccmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22u'
      if ( ipconc .ge. 5 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cwacgh(mgs) = 0.0
      cwacgh(mgs) =
     >   ((0.25)*pi)*eghw(mgs)*cx(mgs,lc)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lgh,1)
     >    + gf3*xdia(mgs,lgh,2) )
      cwacgh(mgs) = min(cwacgh(mgs), fconv*cghmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22v'
      if ( ipconc .ge. 3 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      cghacr(mgs) = 0.0
      cghacr(mgs) =
     >   ((0.25)*pi)*eghr(mgs)*cx(mgs,lr)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgh,1)
     >    + gf1*gf3*xdia(mgs,lgh,2) )
      cghacr(mgs) = min(cghacr(mgs),crmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      if (ndebug .gt. 0 ) print*,'conc 22w'
      do mgs = 1,ngscnt
      cracgh(mgs) = 0.0
      cracgh(mgs) =
     >   ((0.25)*pi)*eghr(mgs)*cx(mgs,lr)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lgh,1)
     >    + gf1*gf3*xdia(mgs,lgh,2) )
      cracgh(mgs) = min(cracgh(mgs),fconv*cghmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22x'
      do mgs = 1,ngscnt
      cghaci(mgs) = 0.0
      cghaci(mgs) =
     >   ((0.25)*pi)*eghi(mgs)*cx(mgs,li)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lgh,1)
     >    + gf3*xdia(mgs,lgh,2) )
      cghaci(mgs) = min(cghaci(mgs),cimxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22y'
      do mgs = 1,ngscnt
      cghacip(mgs) = 0.0
      cghacip(mgs) =
     >   ((0.25)*pi)*eghip(mgs)*cx(mgs,lip)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lgh,1)
     >    + gf3*xdia(mgs,lgh,2) )
      cghacip(mgs) = min(cghacip(mgs),cipmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22z'
      do mgs = 1,ngscnt
      cghacir(mgs) = 0.0
      cghacir(mgs) =
     >   ((0.25)*pi)*eghir(mgs)*cx(mgs,lir)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*xdia(mgs,lir,1)*xdia(mgs,lgh,1)
     >    + gf3*gf1*xdia(mgs,lgh,2) )
      cghacir(mgs) = min(cghacir(mgs),cirmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22aa'
      do mgs = 1,ngscnt
      cghacs(mgs) = 0.0
      cghacs(mgs) =
     >   ((0.25)*pi)*eghs(mgs)*cx(mgs,ls)*cx(mgs,lgh)
     >  *abs(vtxbar(mgs,lgh,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lgh,1)
     >    + gf1*gf3*xdia(mgs,lgh,2) )
      cghacs(mgs) = min(cghacs(mgs),csmxd(mgs))
      end do
      end if
c                    
c
c
c  Frozen drops
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22bb'
      do mgs = 1,ngscnt
      cfacw(mgs) = 0.0
      cfacw(mgs) =
     >   ((0.25)*pi)*efw(mgs)*cx(mgs,lc)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lf,1)
     >    + gf3*xdia(mgs,lf,2) )
      cfacw(mgs) = min(cfacw(mgs),ccmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 5 .or. ipelec .ge. 9 ) then
      if (ndebug .gt. 0 ) print*,'conc 22cc'
      do mgs = 1,ngscnt
      cwacf(mgs) = 0.0
      cwacf(mgs) =
     >   ((0.25)*pi)*efw(mgs)*cx(mgs,lc)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lf,1)
     >    + gf3*xdia(mgs,lf,2) )
      cwacf(mgs) = min(cwacf(mgs), fconv*cfmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 13 .or. ipelec .ge. 9 ) then
      if (ndebug .gt. 0 ) print*,'conc 22cc'
      do mgs = 1,ngscnt
      cfacr(mgs) = 0.0
      cfacr(mgs) =
     >   ((0.25)*pi)*efr(mgs)*cx(mgs,lr)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lf,1)
     >    + gf1*gf3*xdia(mgs,lf,2) )
      cfacr(mgs) = min(cfacr(mgs),crmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 9 .or. ipelec .ge. 9 ) then
      if (ndebug .gt. 0 ) print*,'conc 22dd'
      do mgs = 1,ngscnt
      cracf(mgs) = 0.0
      cracf(mgs) =
     >   ((0.25)*pi)*efr(mgs)*cx(mgs,lr)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lr,1))
     >  *(  gf3*gf1*xdia(mgs,lr,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lf,1)
     >    + gf1*gf3*xdia(mgs,lf,2) )
      cracf(mgs) = min(cracf(mgs),fconv*cfmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22ee'
      do mgs = 1,ngscnt
      cfaci(mgs) = 0.0
      cfaci(mgs) =
     >   ((0.25)*pi)*efi(mgs)*cx(mgs,li)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lf,1)
     >    + gf3*xdia(mgs,lf,2) )
      cfaci(mgs) = min(cfaci(mgs),cimxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22ff'
      do mgs = 1,ngscnt
      cfacip(mgs) = 0.0
      cfacip(mgs) =
     >   ((0.25)*pi)*efip(mgs)*cx(mgs,lip)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lf,1)
     >    + gf3*xdia(mgs,lf,2) )
      cfacip(mgs) = min(cfacip(mgs),cipmxd(mgs))
      end do
      end if
c
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      if (ndebug .gt. 0 ) print*,'conc 22gg'
      do mgs = 1,ngscnt
      cfacir(mgs) = 0.0
      cfacir(mgs) =
     >   ((0.25)*pi)*efir(mgs)*cx(mgs,lir)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*gf2*xdia(mgs,lir,1)*xdia(mgs,lf,1)
     >    + gf3*gf1*xdia(mgs,lf,2) )
      cfacir(mgs) = min(cfacir(mgs),cirmxd(mgs))
      end do
      end if
c
c
      if (ndebug .gt. 0 ) print*,'conc 22hh'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      cfacs(mgs) = 0.0
      cfacs(mgs) =
     >   ((0.25)*pi)*efs(mgs)*cx(mgs,ls)*cx(mgs,lf)
     >  *abs(vtxbar(mgs,lf,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lf,1)
     >    + gf1*gf3*xdia(mgs,lf,2) )
      cfacs(mgs) = min(cfacs(mgs),csmxd(mgs))
      end do
      end if
c                    
c
c
c  Hail
c
      if (ndebug .gt. 0 ) print*,'conc 22ii'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chacw(mgs) = 0.0
      chacw(mgs) =
     >   ((0.25)*pi)*ehw(mgs)*cx(mgs,lc)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lh,1)
     >    + gf3*xdia(mgs,lh,2) )
      chacw(mgs) = min(chacw(mgs),ccmxd(mgs))
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 22jj'
c      if ( ipconc .lt. 3 .and. ipelec .ge. 1 ) then
c      do mgs = 1,ngscnt
c      chacr(mgs) = 0.0
c      chacr(mgs) =
c     >   ((0.25)*pi)*ehr(mgs)*cx(mgs,lr)*cx(mgs,lh)
c     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))
c     >  *(  gf3*gf1*xdia(mgs,lr,2)
c     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lh,1)
c     >    + gf1*gf3*xdia(mgs,lh,2) )
c      chacr(mgs) = min(chacr(mgs),crmxd(mgs))
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22kk'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chaci(mgs) = 0.0
      chaci(mgs) =
     >   ((0.25)*pi)*ehi(mgs)*cx(mgs,li)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lh,1)
     >    + gf3*xdia(mgs,lh,2) )
      chaci(mgs) = min(chaci(mgs),cimxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22ll'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chacip(mgs) = 0.0
      chacip(mgs) =
     >   ((0.25)*pi)*ehip(mgs)*cx(mgs,lip)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lh,1)
     >    + gf3*xdia(mgs,lh,2) )
      chacip(mgs) = min(chacip(mgs),cipmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22mm'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chacir(mgs) = 0.0
      chacir(mgs) =
     >   ((0.25)*pi)*ehir(mgs)*cx(mgs,lir)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*xdia(mgs,lir,1)*xdia(mgs,lh,1)
     >    + gf3*gf1*xdia(mgs,lh,2) )
      chacir(mgs) = min(chacir(mgs),cirmxd(mgs))
      end do
      end if
c
c
      if (ndebug .gt. 0 ) print*,'conc 22nn'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chacs(mgs) = 0.0
      chacs(mgs) =
     >   ((0.25)*pi)*ehs(mgs)*cx(mgs,ls)*cx(mgs,lh)
     >  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)
     >    + gf1*gf3*xdia(mgs,lh,2) )
      chacs(mgs) = min(chacs(mgs),csmxd(mgs))
      end do
      end if
c
c  Hail - Large
c
      if (ndebug .gt. 0 ) print*,'conc 22ii'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chlacw(mgs) = 0.0
      IF ( ehlw(mgs) .gt. 0.0 ) THEN
      chlacw(mgs) =
     >   ((0.25)*pi)*ehlw(mgs)*cx(mgs,lc)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))
     >  *(  gf1*xdia(mgs,lc,2)
     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lhl,1)
     >    + gf3*xdia(mgs,lhl,2) )
      chlacw(mgs) = min(chlacw(mgs),ccmxd(mgs))
      ENDIF
      end do
      end if
c
c      if (ndebug .gt. 0 ) print*,'conc 22jj'
c      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
c      do mgs = 1,ngscnt
c      chlacr(mgs) = 0.0
c      chlacr(mgs) =
c     >   ((0.25)*pi)*ehlr(mgs)*cx(mgs,lr)*cx(mgs,lhl)
c     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))
c     >  *(  gf3*gf1*xdia(mgs,lr,2)
c     >    + 2.0*gf2*gf2*xdia(mgs,lr,1)*xdia(mgs,lhl,1)
c     >    + gf1*gf3*xdia(mgs,lhl,2) )
c      chlacr(mgs) = min(chlacr(mgs),crmxd(mgs))
c      end do
c      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22kk'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chlaci(mgs) = 0.0
      chlaci(mgs) =
     >   ((0.25)*pi)*ehli(mgs)*cx(mgs,li)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))
     >  *(  gf1*xdia(mgs,li,2)
     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lhl,1)
     >    + gf3*xdia(mgs,lhl,2) )
      chlaci(mgs) = min(chlaci(mgs),cimxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22ll'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chlacip(mgs) = 0.0
      chlacip(mgs) =
     >   ((0.25)*pi)*ehlip(mgs)*cx(mgs,lip)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lip,1))
     >  *(  gf1*xdia(mgs,lip,2)
     >    + 2.0*gf2*xdia(mgs,lip,1)*xdia(mgs,lhl,1)
     >    + gf3*xdia(mgs,lhl,2) )
      chlacip(mgs) = min(chlacip(mgs),cipmxd(mgs))
      end do
      end if
c
      if (ndebug .gt. 0 ) print*,'conc 22mm'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chlacir(mgs) = 0.0
      chlacir(mgs) =
     >   ((0.25)*pi)*ehlir(mgs)*cx(mgs,lir)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lir,1))
     >  *(  gf1*gf3*xdia(mgs,lir,2)
     >    + 2.0*gf2*xdia(mgs,lir,1)*xdia(mgs,lhl,1)
     >    + gf3*gf1*xdia(mgs,lhl,2) )
      chlacir(mgs) = min(chlacir(mgs),cirmxd(mgs))
      end do
      end if
c
c
      if (ndebug .gt. 0 ) print*,'conc 22nn'
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      chlacs(mgs) = 0.0
      chlacs(mgs) =
     >   ((0.25)*pi)*ehls(mgs)*cx(mgs,ls)*cx(mgs,lhl)
     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))
     >  *(  gf3*gf1*xdia(mgs,ls,2)
     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lhl,1)
     >    + gf1*gf3*xdia(mgs,lhl,2) )
      chlacs(mgs) = min(chlacs(mgs),csmxd(mgs))
      end do
      end if
        
c
c   Include file for electrification
c

      IF ( ipelec .ge. 1 ) THEN

c
c  set rates to zero
c
      scxacy(:,:,:) = 0.0

      
      IF ( nonigrd .ne. -10 .and.
     :    ( (ny .eq. 2 .or. bcy .eq. 2 ) .or. 
     :     ( ny .gt. 2 .and. jy .ge. 3 .and. jy .lt. ny-2 ) )) THEN
      
      DO mgs = 1,ngscnt


      exy(mgs,ls,lc) = esw(mgs)
      exy(mgs,ls,li) = esi(mgs)
      exy(mgs,ls,lip) = esip(mgs)
      exy(mgs,ls,lir) = esir(mgs)

      exy(mgs,lgl,lc) = eglw(mgs)
      exy(mgs,lgl,li) = egli(mgs)
      exy(mgs,lgl,lip) = eglip(mgs)
      exy(mgs,lgl,lir) = eglir(mgs)
      exy(mgs,lgl,ls) = egls(mgs)

      exy(mgs,lgm,lc) = egmw(mgs)
      exy(mgs,lgm,li) = egmi(mgs)
      exy(mgs,lgm,lip) = egmip(mgs)
      exy(mgs,lgm,lir) = egmir(mgs)
      exy(mgs,lgm,ls) = egms(mgs)

      exy(mgs,lgh,lc) = eghw(mgs)
      exy(mgs,lgh,li) = eghi(mgs)
      exy(mgs,lgh,lip) = eghip(mgs)
      exy(mgs,lgh,lir) = eghir(mgs)
      exy(mgs,lgh,ls) = eghs(mgs)

      exy(mgs,lf,lc) = efw(mgs)
      exy(mgs,lf,li) = efi(mgs)
      exy(mgs,lf,lip) = efip(mgs)
      exy(mgs,lf,lir) = efir(mgs)
      exy(mgs,lf,ls) = efs(mgs)

      exy(mgs,lh,lc) = ehw(mgs)
      exy(mgs,lh,li) = ehi(mgs)
      exy(mgs,lh,lip) = ehip(mgs)
      exy(mgs,lh,lir) = ehir(mgs)
      exy(mgs,lh,ls) = ehs(mgs)

      exy(mgs,lhl,lc) = ehlw(mgs)
      exy(mgs,lhl,li) = ehli(mgs)
      exy(mgs,lhl,lip) = ehlip(mgs)
      exy(mgs,lhl,lir) = ehlir(mgs)
      exy(mgs,lhl,ls) = ehls(mgs)



      cxacy(mgs,ls,lc) = csacw(mgs)
      cxacy(mgs,ls,li) = csaci(mgs)
      cxacy(mgs,ls,lip) = csacip(mgs)
      cxacy(mgs,ls,lir) = csacir(mgs)

      cxacy(mgs,lgl,lc) = cglacw(mgs)
      cxacy(mgs,lgl,li) = cglaci(mgs)
      cxacy(mgs,lgl,lip) = cglacip(mgs)
      cxacy(mgs,lgl,lir) = cglacir(mgs)
      cxacy(mgs,lgl,ls) = cglacs(mgs)

      cxacy(mgs,lgm,lc) = cgmacw(mgs)
      cxacy(mgs,lgm,li) = cgmaci(mgs)
      cxacy(mgs,lgm,lip) = cgmacip(mgs)
      cxacy(mgs,lgm,lir) = cgmacir(mgs)
      cxacy(mgs,lgm,ls) = cgmacs(mgs)

      cxacy(mgs,lgh,lc) = cghacw(mgs)
      cxacy(mgs,lgh,li) = cghaci(mgs)
      cxacy(mgs,lgh,lip) = cghacip(mgs)
      cxacy(mgs,lgh,lir) = cghacir(mgs)
      cxacy(mgs,lgh,ls) = cghacs(mgs)

      cxacy(mgs,lf,lc) = cfacw(mgs)
      cxacy(mgs,lf,li) = cfaci(mgs)
      cxacy(mgs,lf,lip) = cfacip(mgs)
      cxacy(mgs,lf,lir) = cfacir(mgs)
      cxacy(mgs,lf,ls) = cfacs(mgs)

      cxacy(mgs,lh,lc) = chacw(mgs)
      cxacy(mgs,lh,li) = chaci(mgs)
      cxacy(mgs,lh,lip) = chacip(mgs)
      cxacy(mgs,lh,lir) = chacir(mgs)
      cxacy(mgs,lh,ls) = chacs(mgs)

      cxacy(mgs,lhl,lc) = chlacw(mgs)
      cxacy(mgs,lhl,li) = chlaci(mgs)
      cxacy(mgs,lhl,lip) = chlacip(mgs)
      cxacy(mgs,lhl,lir) = chlacir(mgs)
      cxacy(mgs,lhl,ls) = chlacs(mgs)
      
      ENDDO
      

!
!  First do loop between large ice and small ice
!

      DO il=ls,lhab
      DO ic=li,ls
c
      if ( iexy(il,ic) .eq. 1 .and. il .ne. ic ) then
c      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      if ( ndebug .ge. 1 ) print*, 'elec 16'
      do mgs = 1,ngscnt
c
      if ( qx(mgs,il) .gt. qxmin(il) .and. qx(mgs,ic) .gt. qxmin(ic)
     :     .and. exy(mgs,il,ic) .gt. 0.0 ) then
      if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
c
c  del cglarge
c
c
      IF (nonigrd .eq. -1 ) THEN

       vt = vtxbar(mgs,il,1) - vtxbar(mgs,ic,1)
       CALL takax(isaund,ntem,nlwc,takalu,temcg(mgs),qx(mgs,lc),vt,
     :     xdia(mgs,ic,1),rho0(mgs),ftelwc,exy(mgs,il,lc),qx(mgs,lr),rarfac)
       scxacy(mgs,il,ic) = ftelwc

        if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
          scxacy(mgs,il,ic) = min(delqxia,scxacy(mgs,il,ic))
        end if
        if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
          scxacy(mgs,il,ic) = max(delqnia,scxacy(mgs,il,ic))
        end if
      
      ELSEIF ( nonigrd .eq. 1 .or. nonigrd .eq. 0 ) THEN

      IF (nonigrd .eq. 0) THEN
      vt = vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
      call saundx(isaund,temcg(mgs),qx(mgs,lc),eglw(mgs),vt,xdia(mgs,ic,1),
     >   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac)
      
      ELSE
        vt = vtxbar(mgs,il,1) - vtxbar(mgs,lc,1)
        
      call saundy(isaund,temcg(mgs),qx(mgs,lc),rargl(mgs),vt,xdia(mgs,ic,1),
     >   rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq)

c      call saund(temcg(mgs),qx(mgs,lc),eglw(mgs),xdia(mgs,ic,1),
c     >     rho0(mgs),qsign,ftelwc,qconkq,qconm,qconn,idelq)
      ENDIF

      cfce15 = 1.0e-15
      scxacy(mgs,il,ic) = qconkq*(xdia(mgs,ic,1)**qconm)
     >  *((abs(vtxbar(mgs,il,1)-vtxbar(mgs,ic,1)))**qconn)
     >  *ftelwc*cfce15
 
      IF ( ndebug.ge.1 .and. il .eq. lgh .and. vt .gt. 0.1 .and. qx(mgs,il).gt. 1.e-3 ) THEN
        print*,'scxacy: ',il,ic,mgs,scxacy(mgs,il,ic), scxacy(mgs,lgh,li)
      ENDIF
      
      
      if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
      scxacy(mgs,il,ic) = min(delqxia,scxacy(mgs,il,ic))
      end if
      if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
      scxacy(mgs,il,ic) = max(delqnia,scxacy(mgs,il,ic))
      end if
      
      ENDIF
c
      if ( nonigrd .ge. 2 ) then
      scxacy(mgs,il,ic) = 
     > ( (7.3)*rgard1
     > * dellwc(mgs)*sctem(mgs)*(xdia(mgs,ic,1)**4)
     > * (abs(vtxbar(mgs,il,1)-vtxbar(mgs,ic,1))**3) )
      if ( scxacy(mgs,il,ic) .gt. 0.0 ) then
      scxacy(mgs,il,ic) = min(delqxia,scxacy(mgs,il,ic))
      end if
      if ( scxacy(mgs,il,ic) .lt. 0.0 ) then
      scxacy(mgs,il,ic) = max(delqnia,scxacy(mgs,il,ic))
      end if
      end if
c
      if ( nonigrd .eq. 3 ) then
      if ( scxacy(mgs,il,ic) .gt. 0 )  scxacy(mgs,il,ic) = delqxia
      if ( scxacy(mgs,il,ic) .lt. 0 )  scxacy(mgs,il,ic) = delqnia
      end if
c
      sxxacy(mgs,il,ic) = scxacy(mgs,il,ic)
c
      scxacy(mgs,il,ic) =
     >  cxacy(mgs,il,ic)*(1.0-exy(mgs,il,ic))*scxacy(mgs,il,ic)
     >  *(1./(exy(mgs,il,ic)+eps))
c
      end if
      end if
c
      end do
      
      ENDIF
      
      ENDDO
      ENDDO

!
!  End of ice-ice loop
!

!
! Loop for graupel-droplet inductive charging
!

      IF ( erbnd .gt. 0.0 .and. fdgt .gt. 0.0 .and. ipelec .ge. 3 ) THEN
      
      DO il=ls,lhab
c
      IF ( iexy(il,lc) .eq. 1 ) THEN
      
      do mgs = 1,ngscnt
c

c      if ( ieglc .eq. 1 ) then
c      if ( ipelec .ge. 3 .and. itest .eq. 1 ) then
      if ( ndebug .ge. 1 ) print*, 'elec 19'
c
      if ( qx(mgs,il) .gt. qxmin(il) .and. qx(mgs,lc) .gt. qxmin(lc) 
     :     .and. exy(mgs,il,lc) .gt. 0.0 ) THEN
      IF ( temg(mgs) .gt. tindmn .and. temg(mgs) .lt. tindmx ) then
c
      scxacy(mgs,il,lc) =
     >  0.125*pi**3*exy(mgs,il,lc)*erbnd*fdgt*(cx(mgs,il))
     >  *cx(mgs,lc)*vtxbar(mgs,il,1)*6.0/gf4p5*(xdia(mgs,lc,1)**2)
     >  *( pi*eperao*costhe*dezcomp(mgs)
     >      *gf3p5*(xdia(mgs,il,1)**2)
     <     - (1.0/3.0)*gf1p5*scx(mgs,il)/cx(mgs,il) )
c
      if ( scxacy(mgs,il,lc) .gt. 0.0 ) then
        scxacy(mgs,il,lc) = min(delqxw,scxacy(mgs,il,lc))
      elseif ( scxacy(mgs,il,lc) .lt. 0.0 ) then
        scxacy(mgs,il,lc) = max(delqnw,scxacy(mgs,il,lc))
      end if
      
      ELSE
      
       scxacy(mgs,il,lc) = 0.0

      end if
      ENDIF
c
c      end if
c      end if

      end do  ! mgs
      
      ENDIF
      
      ENDDO  ! il
      
      ENDIF

!
! End of Loop for graupel-droplet inductive charging
!


      ENDIF ! nonigrd .ne. -10

      DO mgs = 1,ngscnt
      
      scsacw(mgs) = scxacy(mgs,ls,lc)
      scsaci(mgs) = scxacy(mgs,ls,li)
      scsacip(mgs) = scxacy(mgs,ls,lip)
      scsacir(mgs) = scxacy(mgs,ls,lir)
      scsacr(mgs) = scxacy(mgs,ls,lr)

      scglacw(mgs) = scxacy(mgs,lgl,lc)
      scglaci(mgs) = scxacy(mgs,lgl,li)
      scglacip(mgs) = scxacy(mgs,lgl,lip)
      scglacir(mgs) = scxacy(mgs,lgl,lir)
      scglacs(mgs) = scxacy(mgs,lgl,ls)
      scglacr(mgs) = scxacy(mgs,lgl,lr)

      scgmacw(mgs) = scxacy(mgs,lgm,lc)
      scgmaci(mgs) = scxacy(mgs,lgm,li)
      scgmacip(mgs) = scxacy(mgs,lgm,lip)
      scgmacir(mgs) = scxacy(mgs,lgm,lir)
      scgmacs(mgs) = scxacy(mgs,lgm,ls)
      scgmacr(mgs) = scxacy(mgs,lgm,lr)

      scghacw(mgs) = scxacy(mgs,lgh,lc)
      scghaci(mgs) = scxacy(mgs,lgh,li)
      scghacip(mgs) = scxacy(mgs,lgh,lip)
      scghacir(mgs) = scxacy(mgs,lgh,lir)
      scghacs(mgs) = scxacy(mgs,lgh,ls)
      scghacr(mgs) = scxacy(mgs,lgh,lr)

      scfacw(mgs) = scxacy(mgs,lf,lc)
      scfaci(mgs) = scxacy(mgs,lf,li)
      scfacip(mgs) = scxacy(mgs,lf,lip)
      scfacir(mgs) = scxacy(mgs,lf,lir)
      scfacs(mgs) = scxacy(mgs,lf,ls)
      scfacr(mgs) = scxacy(mgs,lf,lr)

      schacw(mgs) = scxacy(mgs,lh,lc)
      schaci(mgs) = scxacy(mgs,lh,li)
      schacip(mgs) = scxacy(mgs,lh,lip)
      schacir(mgs) = scxacy(mgs,lh,lir)
      schacs(mgs) = scxacy(mgs,lh,ls)
      schacr(mgs) = scxacy(mgs,lh,lr)

      schlacw(mgs) = scxacy(mgs,lhl,lc)
      schlaci(mgs) = scxacy(mgs,lhl,li)
      schlacip(mgs) = scxacy(mgs,lhl,lip)
      schlacir(mgs) = scxacy(mgs,lhl,lir)
      schlacs(mgs) = scxacy(mgs,lhl,ls)
      schlacr(mgs) = scxacy(mgs,lhl,lr)

      IF ( ndebug .ge. 1 .and. qx(mgs,lgh) .gt. 0.001 ) print*,'scghaci(mgs) = ',mgs,scghaci(mgs), sxxacy(mgs,lgh,li)
      
      ENDDO


c      end if
c



c 
      DO mgs = 1,ngscnt
      IF ( bcx .ne. 2 .and. ( igs(mgs) .le. 3 .or. igs(mgs) .ge. nx-3 )) THEN
      scsaci(mgs) = 0.0
      scsacip(mgs) = 0.0
      scsacir(mgs) = 0.0
      scsacw(mgs) = 0.0
      scsacr(mgs) = 0.0

      scglacs(mgs) = 0.0
      scglaci(mgs) = 0.0
      scglacip(mgs) = 0.0
      scglacir(mgs) = 0.0
      scglacw(mgs) = 0.0
      scglacr(mgs) = 0.0

      scgmacs(mgs) = 0.0
      scgmaci(mgs) = 0.0
      scgmacip(mgs) = 0.0
      scgmacir(mgs) = 0.0
      scgmacw(mgs) = 0.0
      scgmacr(mgs) = 0.0

      scghacs(mgs) = 0.0
      scghacir(mgs) = 0.0
      scghacip(mgs) = 0.0
      scghaci(mgs) = 0.0
      scghacw(mgs) = 0.0
      scghacr(mgs) = 0.0

      scfacs(mgs) = 0.0
      scfacir(mgs) = 0.0
      scfacip(mgs) = 0.0
      scfaci(mgs) = 0.0
      scfacw(mgs) = 0.0
      scfacr(mgs) = 0.0

      schacs(mgs) = 0.0
      schacir(mgs) = 0.0
      schacip(mgs) = 0.0
      schaci(mgs) = 0.0
      schacw(mgs) = 0.0
      schacr(mgs) = 0.0

      schlacs(mgs) = 0.0
      schlacir(mgs) = 0.0
      schlacip(mgs) = 0.0
      schlaci(mgs) = 0.0
      schlacw(mgs) = 0.0
      schlacr(mgs) = 0.0

      sxsaci(mgs) = 0.0
      sxsacip(mgs) = 0.0
      sxsacir(mgs) = 0.0
      sxglacs(mgs) = 0.0
      sxglaci(mgs) = 0.0
      sxglacip(mgs) = 0.0
      sxglacir(mgs) = 0.0
      sxgmacs(mgs) = 0.0
      sxgmaci(mgs) = 0.0
      sxgmacip(mgs) = 0.0
      sxgmacir(mgs) = 0.0
      sxghacs(mgs) = 0.0
      sxghacir(mgs) = 0.0
      sxghacip(mgs) = 0.0
      sxghaci(mgs) = 0.0
      sxfacs(mgs) = 0.0
      sxfaci(mgs) = 0.0
      sxfacir(mgs) = 0.0
      sxfacip(mgs) = 0.0
      sxhacs(mgs) = 0.0
      sxhaci(mgs) = 0.0
      sxhacir(mgs) = 0.0
      sxhacip(mgs) = 0.0
      sxhacs(mgs) = 0.0
      sxhlaci(mgs) = 0.0
      sxhlacir(mgs) = 0.0
      sxhlacip(mgs) = 0.0
      ENDIF
      END DO

      ENDIF ! ( ipelec .ge. 1 .and. nonigrd .ne. -10)


c
c
c
c
c Ziegler (1985) autoconversion
c
c
      IF ( ipconc .ge. 2 ) THEN
      
      DO mgs = 1,ngscnt
        qrcnw(mgs) = 0.0
        crcnw(mgs) = 0.0
      ENDDO
      
      DO mgs = 1,ngscnt
c      qracw(mgs) = 0.0
c      cracw(mgs) = 0.0
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1000.  .and. temg(mgs) .gt. tfrh+4. ) THEN
         volb = xv(mgs,lc)*SQRT(1./(1.+CNU))
         cautn(mgs) = Min(ccmxd(mgs),
     :      ((CNU+2.)/(CNU+1.))*aa1*cx(mgs,lc)**2*xv(mgs,lc)**2)
         cautn(mgs) = Max( 0.0, cautn(mgs) )
         
         IF ( rb(mgs) .le. 7.51e-6 ) THEN
           t2s = 1.e30
c           cautn(mgs) = 0.0
         ELSE
c         XL2P=2.7E-2*XNC*XVC*((1.E12*RB**3*RC)-0.4)
         
c        T2S=3.72E-3/(((1.E4*RB)-7.5)*XNC*XVC) 
c           t2s = 3.72E-3/(((1.e6*rb)-7.5)*cx(mgs,lc)*xv(mgs,lc))
!           t2s = 3.72/(((1.e6*rb(mgs))-7.5)*rho0(mgs)*qx(mgs,lc))
           t2s = 3.72/(1.e6*(rb(mgs)-7.500e-6)*rho0(mgs)*qx(mgs,lc))

           qrcnw(mgs) = xl2p(mgs)/(t2s*rho0(mgs))
           crcnw(mgs) = Max( 0.0, Min(3.5e9*xl2p(mgs)/t2s,cautn(mgs)) )

c           IF (  crcnw(mgs) .gt. cautn(mgs) .and. crcnw(mgs) .gt. 1.0 )
c     :          THEN
c             print*, 'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
c     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s,qx(mgs,lr)
c             print*, '            ',qx(mgs,lc),cx(mgs,lc),0.5e6*xdia(mgs,lc,1)
c             print*, '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
c     :         1.e6*(( 3/(4.*pi))*rho0(mgs)*qrcnw(mgs)/
c     :       (crcnw(mgs)*xdn(mgs,lr)))**(1./3.),rh(mgs)*1.e6,rwrad(mgs)
c           ELSEIF ( crcnw(mgs) .gt. 1.0 .and. cautn(mgs) .gt. 0.) THEN
c             print*, 'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
c     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s
c             print*, '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
c     :  1.e6*(( 3*pi/4.)*rho0(mgs)*qrcnw(mgs)/
c     :   (crcnw(mgs)*xdn(mgs,lr)))**(1./3.)
c           ENDIF
c           crcnw(mgs) = Min(cautn(mgs),3.5e9*xl2p(mgs)/t2s)
           
c           IF ( qrcnw(mgs) .gt. 0.3e-2 ) THEN
c            print*, 'QRCNW'
c            print*, qrcnw(mgs),crcnw(mgs),cautn(mgs)
c            print*, xl2p,t2s,rho0(mgs),xv(mgs,lc),cx(mgs,lc),qx(mgs,lc)
c            print*, rb,0.5*xdia(mgs,lc,1),mgs,igs(mgs),kgs(mgs)
c           ENDIF
c           qrcnw(mgs) = Min(qrcnw(mgs),qcmxd(mgs))
         ENDIF
         
         
       ENDIF
      ENDDO

      
      
      ELSE
      
c
c  Berry 1968 auto conversion for rain (Orville & Kopp 1977) with delay
c
c
      if ( ircnw .eq. 3 ) then
      do mgs = 1,ngscnt
c      sconvmix(lcw,mgs) = 0.0
      qrcnw(mgs) =  0.0

      timflg(mgs) = 1.0 
      if ( iauttim .eq. 1 ) then
      if ( rtcw(mgs) .lt. auttim ) then
      timflg(mgs) = 0.0
      end if
      end if

      qdiff = max((qx(mgs,lc)-qminrncw),0.0)
      if ( qdiff .gt. 0.0 ) then
      argrcnw =
     >  ((1.2e-4)+(1.596e-12)*(cx(mgs,lc)*1.0e-6)
     >  /(cwdisp*qdiff*1.0e-3*rho0(mgs)))
      qrcnw(mgs) = timflg(mgs)*rho0(mgs)*1.0e-3*(qdiff**2)/argrcnw
c      sconvmix(lcw,mgs) = max(sconvmix(lcw,mgs),0.0)
      qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      end do
      
      ENDIF
c
c
c  Berry 1968 auto conversion for rain (Orville & Kopp 1977)
c
c
      if ( ircnw .eq. 4 ) then
      do mgs = 1,ngscnt
c      sconvmix(lcw,mgs) = 0.0
      qrcnw(mgs) =  0.0
      qdiff = max((qx(mgs,lc)-qminrncw),0.0)
      if ( qdiff .gt. 0.0 .and. xdia(mgs,lc,1) .gt. 20.0e-6 ) then
      argrcnw =
     >  ((1.2e-4)+(1.596e-12)*(cx(mgs,lc)*1.0e-6)
     >  /(cwdisp*qdiff*1.0e-3*rho0(mgs)))
      qrcnw(mgs) = (rho0(mgs)*1e-3)*(qdiff**2)/argrcnw
c      sconvmix(lcw,mgs) = max(sconvmix(lcw,mgs),0.0)
      qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      end do
      
      ENDIF
c
c
c
c  Berry 1968 auto conversion for rain (Ferrier 1994) with time delay
c
c
      if ( ircnw .eq. 5 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      timflg(mgs) = 1.0 
      if ( iauttim .eq. 1 ) then
      if ( rtcw(mgs) .lt. auttim ) then
      timflg(mgs) = 0.0
      end if
      end if
      qrcnw(mgs) =  0.0
      qccrit = (pi/6.)*(cx(mgs,lc)*cwdiap**3)*xdn(mgs,lc)/rho0(mgs)
      qdiff = max((qx(mgs,lc)-qccrit),0.)
      if ( qdiff .gt. 0.0 ) then
      argrcnw = 
c     >  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)/(cwdisp*rho0(mgs)*qdiff))
     >  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)*1.0e-3/(cwdisp*rho0(mgs)*qdiff))
      qrcnw(mgs) =
c     >  timflg(mgs)*rho0(mgs)*(qdiff**2)/argrcnw
     >  timflg(mgs)*1.0e-3*rho0(mgs)*(qdiff**2)/argrcnw
      qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      end do
      end if
      
c
c
c  kessler auto conversion for rain.
c
      if ( ircnw .eq. 2 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      qrcnw(mgs) = (0.001)*max((qx(mgs,lc)-qminrncw),0.0)
      end do
      end if
c
c  c4 = pi/6
c  c1 = 0.12-0.32 for colorado storms...typically 0.3-0.4
c  berry reinhart type conversion (proctor 1988)
c
      if ( ircnw .eq. 1 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      c1 = 0.2
      c4 = pi/(6.0)
      bradp = 
     > (1.e+06) * ((c1/(0.38))**(1./3.)) * (xdia(mgs,lc,1)*(0.5))
      bl2 =
     > (0.027) * ((100.0)*(bradp**3)*(xdia(mgs,lc,1)*(0.5)) - (0.4))
      bt2 = (bradp -7.5) / (3.72)
      qrcnw(mgs) = 0.0
      if ( bl2 .gt. 0.0 .and. bt2 .gt. 0.0 ) then
      qrcnw(mgs) = bl2 * bt2 * rho0(mgs)
     >  * qx(mgs,lc) * qx(mgs,lc)
      end if
      end do
      end if
      
      
      ENDIF  !  ( ipconc .ge. 2 )

c
c  Aggregation of crystals
c
      if (ndebug .gt. 0 ) print*,'conc 29a'
      do mgs = 1,ngscnt 
      qscni(mgs) =  0.0
      cscni(mgs) = 0.0
      if ( qx(mgs,li) .gt. qxmin(li) ) then
      qscni(mgs) = 
     >  pi*rho0(mgs)*((0.25)/(6.0))
     >  *eii(mgs)*(qx(mgs,li)**2)*(xdia(mgs,li,2))
     <  *vtxbar(mgs,li,1)/xmas(mgs,li)
      cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li))
      end if
      end do
c
c
      if (ndebug .gt. 0 ) print*,'conc 28a'
      do mgs = 1,ngscnt
      qscnip(mgs) =  0.0
      cscnip(mgs) = 0.0
      if ( qx(mgs,lip) .gt. qxmin(li) ) then
      qscnip(mgs) =
     >  pi*rho0(mgs)*((0.25)/(6.0))
     >  *eipip(mgs)*(qx(mgs,lip)**2)*(xdia(mgs,lip,2))
     <  *vtxbar(mgs,lip,1)/xmas(mgs,lip)
        IF ( ipconc .ge. 1 ) THEN
          cscnip(mgs) = 
     :           Min(cipmxd(mgs),qscnip(mgs)*rho0(mgs)/xmas(mgs,lip))
        ENDIF
      end if
      end do
c
      if (ndebug .gt. 0 ) print*,'conc 27a'
      do mgs = 1,ngscnt
      qscnir(mgs) =  0.0
      cscnir(mgs) = 0.0
      if ( qx(mgs,lir) .gt. qxmin(li) ) then
      qscnir(mgs) =
     >  pi*rho0(mgs)*((0.25)/(6.0))
     >  *eirir(mgs)*(qx(mgs,lir)**2)*(xdia(mgs,lir,2))
     <  *vtxbar(mgs,lir,1)/xmas(mgs,lir)
      IF ( ipconc .ge. 1 ) THEN
      cscnir(mgs) = Min(cirmxd(mgs),qscnir(mgs)*rho0(mgs)/xmas(mgs,lir))
      ENDIF
      end  if
      end do
c
c  Bigg Freezing of Rain
c
      if (ndebug .gt. 0 ) print*,'conc 27a'
      do mgs = 1,ngscnt 
c      qrfrz(mgs) =  0.0
c      crfrz(mgs) =  0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) .and. temcg(mgs) .lt. 0. ) then
c      brz = 100.0
c      arz = 0.66
       IF ( ipconc .lt. 3 ) THEN
       qrfrz(mgs) = 
     >  min(
     >  (20.0)*(pi**2)*brz*(xdn(mgs,lr)/rho0(mgs))
     >   *cx(mgs,lr)*(xdia(mgs,lr,1)**6)
     >   *(exp(max(-arz*temcg(mgs), 0.0))-1.0)
     <  , qrmxd(mgs))
       ELSEIF ( ipconc .ge. 3 ) THEN
         tmp = brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
         tmp = xv(mgs,lr)*brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
         crfrz(mgs) = tmp

         qrfrz(mgs) = 6.0*xdn(mgs,lr)*xv(mgs,lr)*rhoinv(mgs)*crfrz(mgs)
c         qrfrz(mgs) = 6.0*xdn(mgs,lr)*xv(mgs,lr)**2*tmp*rhoinv(mgs)
         qrfrz(mgs) = Min( qrfrz(mgs), ffrz*qrmxd(mgs) )
         crfrz(mgs) = Min( crmxd(mgs), ffrz*crfrz(mgs))
!         crfrzf(mgs) = rz*crfrz(mgs)
         crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*vfrz*1000.0)*rho0(mgs) ) ! rz*crfrz(mgs)
c         crfrz(mgs) = Min(crmxd(mgs),qrfrz(mgs)*rho0(mgs)/xmas(mgs,lr))
       ENDIF
c      if ( temg(mgs) .gt. 268.15 ) then
      else
      qrfrz(mgs) = 0.0
      qrfrzf(mgs) = 0.0
      crfrz(mgs) = 0.0
      crfrzf(mgs) = 0.0
c      end if
      end if
      end do
c
c  Homogeneous freezing of cloud drops to ice crystals
c  following Bigg (1953) and Ferrier (1994).
c
      if (ndebug .gt. 0 ) print*,'conc 25b'
      do mgs = 1,ngscnt
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      qwfrzc(mgs) = 0.0
      cwfrzc(mgs) = 0.0
      qwfrzp(mgs) = 0.0
      cwfrzp(mgs) = 0.0
      IF ( ibfc .ge. 1 ) THEN
      if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1. .and.
     :     .not. (ipconc .ge. 2 .and. xdia(mgs,lc,1) .lt. 18.e-6) ) then
      qwfrz(mgs) = ((2.0)*(brz)/(xdn(mgs,lc)*cx(mgs,lc)))
     >  *(exp(max(-arz*temcg(mgs), 0.0))-1.0)
     >  *rho0(mgs)*(qx(mgs,lc)**2)
      qwfrz(mgs) = max(qwfrz(mgs), 0.0)
      qwfrz(mgs) = min(qwfrz(mgs),qcmxd(mgs))
       IF ( ipconc .ge. 2 ) THEN
         cwfrz(mgs) = 0.5*qwfrz(mgs)*rho0(mgs)/xmas(mgs,lc)
         cwfrz(mgs) = Min(cwfrz(mgs),ccmxd(mgs))
       ELSE
         cwfrz(mgs) = qwfrz(mgs)*rho0(mgs)/xmas(mgs,li)
       ENDIF
      if ( temg(mgs) .gt. 268.15 ) then
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      end if
      end if
      ENDIF
c
      if ( xplate(mgs) .eq. 1 ) then
      qwfrzp(mgs) = qwfrz(mgs)
      cwfrzp(mgs) = cwfrz(mgs)
      end if
c
      if ( xcolmn(mgs) .eq. 1 ) then
      qwfrzc(mgs) = qwfrz(mgs)
      cwfrzc(mgs) = cwfrz(mgs)
      end if
c
c     qwfrzp(mgs) = 0.0
c     qwfrzc(mgs) = qwfrz(mgs)
c
      end do 
c
c
c  Contact freezing nucleation:  factor is to convert from L-1
c  T < -2C:  via Meyers et al. JAM July, 1992 (31, 708-721)
c
      if (ndebug .gt. 0 ) print*,'conc 25a'
      do mgs = 1,ngscnt
      cwctfz(mgs) = 0.0
      qwctfz(mgs) = 0.0
      cwctfzc(mgs) = 0.0
      qwctfzc(mgs) = 0.0
      cwctfzp(mgs) = 0.0
      qwctfzp(mgs) = 0.0
      IF ( icfn .ge. 1 .and. ipconc .eq. 0 ) THEN
      if ( temg(mgs) .lt. 271.15  .and. qx(mgs,lc) .gt. qxmin(lc)) then
      cwctfz(mgs) = cfnfac*exp( (-2.80) - (0.262)*temcg(mgs) )
      cwctfz(mgs) = Min((1.0e3)*cwctfz(mgs)/dtp, ccmxd(mgs) )
      IF ( ipconc .ge. 2 ) THEN
        cwctfz(mgs) = Min( cwctfz(mgs), ccmxd(mgs) )
        qwctfz(mgs) = xmas(mgs,lc)*cwctfz(mgs)/rho0(mgs)
      ELSE
        qwctfz(mgs) = (cimasn)*cwctfz(mgs)/rho0(mgs)
        qwctfz(mgs) = max(qwctfz(mgs), 0.0)
        qwctfz(mgs) = min(qwctfz(mgs),qcmxd(mgs))
      ENDIF
c
      if ( xplate(mgs) .eq. 1 ) then
      qwctfzp(mgs) = qwctfz(mgs)
      cwctfzp(mgs) = cwctfz(mgs)
      end if
c
      if ( xcolmn(mgs) .eq. 1 ) then
      qwctfzc(mgs) = qwctfz(mgs)
      cwctfzc(mgs) = cwctfz(mgs)
      end if
c
c     qwctfzc(mgs) = qwctfz(mgs)
c     qwctfzp(mgs) = 0.0
c
      end if
      ENDIF ! icfn
      end do
c
c
c
c Hobbs-Rangno ice enhancement (Ferrier, 1994)
c
      if (ndebug .gt. 0 ) print*,'conc 23a'
      dtrh = 300.0
      hrifac = (1.e-3)*((0.044)*(0.01**3))
      do mgs = 1,ngscnt
      ciihr(mgs) = 0.0
      qiihr(mgs) = 0.0
      cicichr(mgs) = 0.0
      qicichr(mgs) = 0.0
      cipiphr(mgs) = 0.0
      qipiphr(mgs) = 0.0
      IF ( ihrn .ge. 1 ) THEN
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
      if ( temg(mgs) .lt. 273.15 ) then
c      write(iunit,'(3(1x,i3),3(1x,1pe12.5))')
c     : igs(mgs),jgs,kgs(mgs),cx(mgs,lc),rho0(mgs),qx(mgs,lc)
c      write(iunit,'(1pe15.6)')
c     :  log(cx(mgs,lc)*(1.e-6)/(3.0)),
c     :  ((1.e-3)*rho0(mgs)*qx(mgs,lc)),
c     :  (cx(mgs,lc)*(1.e-6)),
c     : ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)),
c     : (alog(cx(mgs,lc)*(1.e-6)/(3.0)) *
c     >  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6))) 

      IF ( Log(cx(mgs,lc)*(1.e-6)/(3.0)) .gt. 0.0 ) THEN
      ciihr(mgs) = ((1.69e17)/dtrh)
     > *(log(cx(mgs,lc)*(1.e-6)/(3.0)) *
     >  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)))**(7./3.)
      ciihr(mgs) = ciihr(mgs)*(1.0e6)
      qiihr(mgs) = hrifac*ciihr(mgs)/rho0(mgs)
      qiihr(mgs) = max(qiihr(mgs), 0.0)
      qiihr(mgs) = min(qiihr(mgs),qcmxd(mgs))
      ENDIF
c
      if ( xplate(mgs) .eq. 1 ) then
      qipiphr(mgs) = qiihr(mgs)
      cipiphr(mgs) = ciihr(mgs)
      end if
c
      if ( xcolmn(mgs) .eq. 1 ) then
      qicichr(mgs) = qiihr(mgs)
      cicichr(mgs) = ciihr(mgs)
      end if
c
c     qipiphr(mgs) = 0.0
c     qicichr(mgs) = qiihr(mgs)
c
      end if
      end if
      ENDIF ! ihrn
      end do
c
c
c
c  simple frozen rain to hail conversion.  All of the 
c  frozen rain larger than 5.0e-3 m in diameter are converted
c  to hail.  This is done by considering the equation for
c  frozen rain mixing ratio:
c
c
c  qfw = [ cno(lf) * pi * fwdn / (6 rhoair) ] 
c
c         /inf
c      *  |     fwdia*3 exp(-dia/fwdia) d(dia)
c         /Do
c
c  The amount to be reclassified as hail is the integral above from
c  Do to inf where Do is 5.0e-3 m.
c
c
c  qfauh = [ cno(lf) * pi * fwdn / (6 rhoair) ] 
c
c
      if (ndebug .gt. 0 ) print*,'conc 24'
      hdia0 = 05.0e-3
      do mgs = 1,ngscnt
      qhcnf(mgs) = 0.0
      chcnf(mgs) = 0.0
c      if ( temg(mgs) .lt. tfr-20 ) then
c      if ( qfacw(mgs).gt.1.e-6 .or. qfacr(mgs).gt.1.e-6 ) then
c      qhcnf(mgs) =
c     >  ((pi*xdn(mgs,lf)*cx(mgs,lf)) / (6.0*rho0(mgs)*dtp))
c     > *exp(-hdia0/xdia(mgs,lf,1))
c     > *( (hdia0**3) + 3.0*(hdia0**2)*xdia(mgs,lf,1) 
c     >  + 6.0*(hdia0)*(xdia(mgs,lf,1)**2) + 6.0*(xdia(mgs,lf,1)**3) )
c      qhcnf(mgs) =
c     >  min(qhcnf(mgs),qfmxd(mgs))
c      IF ( ipconc .ge. 5 ) THEN
c        chcnf(mgs) = Min( cfmxd(mgs), cx(mgs,lf)*Exp(-hdia0/xdia(mgs,lf,1)))
c      ENDIF
c      end if
c      end if
      end do
c
      hdia0 = 05.0e-3
      do mgs = 1,ngscnt
      qhcngh(mgs) = 0.0
      chcngh(mgs) = 0.0
      if ( temg(mgs) .lt. tfr ) then
      if ( qghacw(mgs).gt.1.e-6 .or. qghacr(mgs).gt.1.e-6 ) then
      qhcngh(mgs) =
     >  ((pi*xdn(mgs,lgh)*cx(mgs,lgh)) / (6.0*rho0(mgs)*dtp))
     > *exp(-hdia0/xdia(mgs,lgh,1))
     > *( (hdia0**3) + 3.0*(hdia0**2)*xdia(mgs,lgh,1)
     >  + 6.0*(hdia0)*(xdia(mgs,lgh,1)**2) + 6.0*(xdia(mgs,lgh,1)**3) )
      qhcngh(mgs) =
     >  min(qhcngh(mgs),qghmxd(mgs))
      IF ( ipconc .ge. 5 ) THEN
        chcngh(mgs) = Min( cghmxd(mgs), cx(mgs,lgh)*Exp(-hdia0/xdia(mgs,lgh,1)))
      ENDIF
      end if
      end if
      end do
c
      hdia0 = 05.0e-3
      do mgs = 1,ngscnt
      qhcngm(mgs) = 0.0
      qhcngm(mgs) = 0.0
      if ( temg(mgs) .lt. tfr ) then
      if ( qgmacw(mgs).gt.1.e-6 .or. qgmacr(mgs).gt.1.e-6 ) then
      qhcngm(mgs) =
     >  ((pi*xdn(mgs,lgm)*cx(mgs,lgm)) / (6.0*rho0(mgs)*dtp))
     > *exp(-hdia0/xdia(mgs,lgm,1))
     > *( (hdia0**3) + 3.0*(hdia0**2)*xdia(mgs,lgm,1)
     >  + 6.0*(hdia0)*(xdia(mgs,lgm,1)**2) + 6.0*(xdia(mgs,lgm,1)**3) )
      qhcngm(mgs) =
     >  min(qhcngm(mgs),qgmmxd(mgs))
      IF ( ipconc .ge. 5 ) THEN
        chcngm(mgs) = Min( cgmmxd(mgs), cx(mgs,lgm)*Exp(-hdia0/xdia(mgs,lgm,1)))
      ENDIF
      end if
      end if
      end do
c
      hdia0 = 05.0e-3
      do mgs = 1,ngscnt
      qhcngl(mgs) = 0.0
      chcngl(mgs) = 0.0
c      if ( temg(mgs) .lt. tfr ) then
c      if ( qglacw(mgs).gt.1.e-6 .or. qglacr(mgs).gt.1.e-6 ) then
c      qhcngl(mgs) =
c     >  ((pi*xdn(mgs,lgl)*cx(mgs,lgl)) / (6.0*rho0(mgs)*dtp))
c     > *exp(-hdia0/xdia(mgs,lgl,1))
c     > *( (hdia0**3) + 3.0*(hdia0**2)*xdia(mgs,lgl,1)
c     >  + 6.0*(hdia0)*(xdia(mgs,lgl,1)**2) + 6.0*(xdia(mgs,lgl,1)**3) )
c      qhcngl(mgs) =
c     >  min(qhcngl(mgs),qglmxd(mgs))
c      end if
c      end if
      end do
c
c large hail
c
      hdia1 = 20.0e-3
      do mgs = 1,ngscnt
      qhlcnh(mgs) = 0.0
      chlcnh(mgs) = 0.0
      if ( temg(mgs) .lt. tfr ) then
      if ( qhacw(mgs).gt.1.e-6 .or. qhacr(mgs).gt.1.e-6 ) then
      qhlcnh(mgs) =
     >  ((pi*xdn(mgs,lh)*cx(mgs,lh)) / (6.0*rho0(mgs)*dtp))
     > *exp(-hdia1/xdia(mgs,lh,1))
     > *( (hdia1**3) + 3.0*(hdia1**2)*xdia(mgs,lh,1)
     >  + 6.0*(hdia1)*(xdia(mgs,lh,1)**2) + 6.0*(xdia(mgs,lh,1)**3) )
      qhlcnh(mgs) =
     >  min(qhlcnh(mgs),qhmxd(mgs))
      IF ( ipconc .ge. 5 ) THEN
        chlcnh(mgs) = Min( chmxd(mgs), cx(mgs,lh)*Exp(-hdia1/xdia(mgs,lh,1)))
      ENDIF
      end if
      end if
      end do


      hdia0 = 300.0e-6
      do mgs = 1,ngscnt
      qscnvi(mgs) = 0.0
      cscnvi(mgs) = 0.0
!      IF ( .false. ) THEN
      if ( ipconc .ge. 4 .and. temg(mgs) .lt. tfr ) then
      if ( cx(mgs,li) .gt. 10. .and. xdia(mgs,li,1) .gt. 20.e-6 ) then
      cirdiatmp =
     >  (qx(mgs,li)*rho0(mgs)
     > /(pi*xdn(mgs,li)*cx(mgs,li)))**(0.333333)
      IF ( cirdiatmp .gt. 50.e-6 ) THEN
      qscnvi(mgs) =
     >  ((pi*xdn(mgs,li)*cx(mgs,li)) / (6.0*rho0(mgs)*dtp))
     > *exp(-hdia0/cirdiatmp)
     > *( (hdia0**3) + 3.0*(hdia0**2)*cirdiatmp
     >  + 6.0*(hdia0)*(cirdiatmp**2) + 6.0*(cirdiatmp**3) )
      qscnvi(mgs) =
     >  min(qscnvi(mgs),qimxd(mgs))
        IF ( ipconc .ge. 4 ) THEN
          cscnvi(mgs) = Min( cimxd(mgs), cx(mgs,li)*Exp(-hdia0/cirdiatmp))
        ENDIF
      ENDIF
      end if
      end if
!      ENDIF
      end do
c
c  Ventilation coeficients
c
      do mgs = 1,ngscnt
      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))
      end do
c
c
      if ( ndebug .gt. 0 ) print*,'civent'
c
      civenta = 1.258e4
      civentb = 2.331
      civentc = 5.662e4
      civentd = 2.373
      civente = 0.8241
      civentf = -0.042
      civentg = 1.70
      do mgs = 1,ngscnt
      IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
      cireyn =
     >  (civenta*xdia(mgs,li,1)**civentb
     >  +civentc*xdia(mgs,li,1)**civentd)
     >  /
     >  (civente*xdia(mgs,li,1)**civentf+civentg)
      xcivent = (fschm(mgs)**(1./3.))*((cireyn/fakvisc(mgs))**0.5)
      if ( xcivent .lt. 1.0 ) then
      civent(mgs) = 1.0 + 0.14*xcivent**2
      end if
      if ( xcivent .ge. 1.0 ) then
      civent(mgs) = 0.86 + 0.28*xcivent
      end if
      ELSE
       civent(mgs) = 0.0
      ENDIF
      end do
c
      if ( ndebug .gt. 0 ) print*,'cipvent'
c
      cipventa = 58.25
      cipventb = 1.824
      cipventc = 0.1652
      cipventd = 1.2990
      do mgs = 1,ngscnt
      IF ( qx(mgs,lip) .gt. qxmin(lip) ) THEN
      cipreyn =
     >  (cipventa*xdia(mgs,lip,1)**cipventb
     >  +cipventc*xdia(mgs,lip,1)**cipventd)
      cipvent(mgs) =
     > (0.65 + 0.44*(fschm(mgs)**(1./3.))*
     >         Sqrt(cipreyn/fakvisc(mgs)))
      ELSE
       cipvent(mgs) = 0.0
      ENDIF
      end do
c
      if ( ndebug .gt. 0 ) print*,'cirvent'
c
      cirventa = (0.78)
      cirventb = (0.308)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lir) .gt. qxmin(lir) ) THEN
      cirvent(mgs) =
     >  ( cirventa + cirventb*fvent(mgs)
     <    *Sqrt(xdia(mgs,lir,1)*vtxbar(mgs,lir,1)) )
      ELSE
      cirvent(mgs) = 0.0
      ENDIF
      end do
c
      igmrwa = 100.0*2.0
      igmrwb = 100.*((5.0+br)/2.0)
      rwventa = (0.78)*gmoi(igmrwa)  ! 0.78
      rwventb = (0.308)*gmoi(igmrwb) ! 0.562825
      do mgs = 1,ngscnt
      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
        IF ( ipconc .ge. 3 ) THEN
        rwvent(mgs) = ventr*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)
        ELSE
         rwvent(mgs) =
     >  (rwventa + rwventb*fvent(mgs)
     >   *Sqrt((ar*rhovt(mgs)))
     <    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
        ENDIF
      ELSE
       rwvent(mgs) = 0.0
      ENDIF
      end do
c
      igmswa = 100.0*2.0
      igmswb = 100.*((5.0+ds)/2.0)
      swventa = (0.78)*gmoi(igmswa)
      swventb = (0.308)*gmoi(igmswb)
      do mgs = 1,ngscnt
      IF ( qx(mgs,ls) .gt. qxmin(ls) ) THEN
      swvent(mgs) =
     >  (swventa + swventb*fvent(mgs)
     >   *Sqrt((cs*rhovt(mgs)))
     <   *(xdia(mgs,ls,1)**((1.0+ds)/2.0)) )
      ELSE
      swvent(mgs) = 0.0
      ENDIF
      end do
c
      igmfwa = 100.0*2.0
      igmfwb = 100.0*2.75
      fwventa = (0.78)*gmoi(igmfwa)
      fwventb = (0.308)*gmoi(igmfwb)
      fwventc = (4.0*gr/(3.0*cdx(lf)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lf) .gt. qxmin(lf) ) THEN
      fwvent(mgs) =
     >  ( fwventa + fwventb*fwventc*fvent(mgs)
     >    *((xdn(mgs,lf)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lf,1)**(0.75)))
      ELSE
       fwvent(mgs) = 0.0
      ENDIF
      end do
c
      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
      hwventc = (4.0*gr/(3.0*cdx(lh)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lh) .gt. qxmin(lh) ) THEN
      hwvent(mgs) =
     >  ( hwventa + hwventb*hwventc*fvent(mgs)
     >    *((xdn(mgs,lh)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lh,1)**(0.75)))
      ELSE
      hwvent(mgs) = 0.0
      ENDIF
      end do
c
      igmhla = 100.0*2.0
      igmhlb = 100.0*2.75
      hlventa = (0.78)*gmoi(igmhla)
      hlventb = (0.308)*gmoi(igmhlb)
      hlventc = (4.0*gr/(3.0*cdx(lhl)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lhl) .gt. qxmin(lh) ) THEN
      hlvent(mgs) =
     >  ( hlventa + hlventb*hlventc*fvent(mgs)
     >    *((xdn(mgs,lhl)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lhl,1)**(0.75)))
      ELSE
      hlvent(mgs) = 0.0
      ENDIF
      end do
c
      igmgla = 100.0*2.0
      igmglb = 100.0*2.75
      glventa = (0.78)*gmoi(igmgla)
      glventb = (0.308)*gmoi(igmglb)
      glventc = (4.0*gr/(3.0*cdx(lgl)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lgl) .gt. qxmin(lgl) ) THEN
      glvent(mgs) =
     >  ( glventa + glventb*glventc*fvent(mgs)
     >    *((xdn(mgs,lgl)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lgl,1)**(0.75)))
      ELSE
      glvent(mgs) = 0.0
      ENDIF
      end do
c
      igmgma = 100.0*2.0
      igmgmb = 100.0*2.75
      gmventa = (0.78)*gmoi(igmgma)
      gmventb = (0.308)*gmoi(igmgmb)
      gmventc = (4.0*gr/(3.0*cdx(lgm)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lgm) .gt. qxmin(lgm) ) THEN
      gmvent(mgs) =
     >  ( gmventa + gmventb*gmventc*fvent(mgs)
     >    *((xdn(mgs,lgm)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lgm,1)**(0.75)))
      ELSE
      gmvent(mgs) = 0.0
      ENDIF
      end do
c
      igmgha = 100.0*2.0
      igmghb = 100.0*2.75
      ghventa = (0.78)*gmoi(igmgha)
      ghventb = (0.308)*gmoi(igmghb)
      ghventc = (4.0*gr/(3.0*cdx(lgh)))**(0.25)
      do mgs = 1,ngscnt
      IF ( qx(mgs,lgh) .gt. qxmin(lgh) ) THEN
      ghvent(mgs) =
     >  ( ghventa + ghventb*ghventc*fvent(mgs)
     >    *((xdn(mgs,lgh)/rho0(mgs))**(0.25))
     <    *(xdia(mgs,lgh,1)**(0.75)))
      ELSE
      ghvent(mgs) = 0.0
      ENDIF
      end do
c
c  Wet growth constants
c
      do mgs = 1,ngscnt
      fwet1(mgs) =
     > (2.0*pi)*
     > ( felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qss0(mgs)-qx(mgs,lv))
     >  -ftka(mgs)*temcg(mgs) )
     > / ( rho0(mgs)*(felf(mgs)+fcw(mgs)*temcg(mgs)) )
      fwet2(mgs) =
     >  (1.0)-fci(mgs)*temcg(mgs)
     > / ( felf(mgs)+fcw(mgs)*temcg(mgs) )
      end do 
c
c  Melting constants
c
      do mgs = 1,ngscnt
      fmlt1(mgs) = (2.0*pi)*
     >  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))
     >   -ftka(mgs)*temcg(mgs)/rho0(mgs) ) 
     >  / (felf(mgs))
      fmlt2(mgs) = -fcw(mgs)*temcg(mgs)/felf(mgs) 
      end do
c
c  Vapor Deposition constants
c
      do mgs = 1,ngscnt
      fvds(mgs) = 
     >  (4.0*pi/rho0(mgs))*(ssi(mgs)-1.0)*
     >  (1.0/(fai(mgs)+fbi(mgs)))
      end do
      do mgs = 1,ngscnt
      fvce(mgs) = 
     >  (4.0*pi/rho0(mgs))*(ssw(mgs)-1.0)*
     >  (1.0/(fav(mgs)+fbv(mgs)))
      end do
c
c  deposition, sublimation, and melting of snow, graupel and hail
c
      do mgs = 1,ngscnt
c
      qimlr(mgs) = 0.0
c     >   min(
c     >  (fmlt1(mgs)*cx(mgs,li)*civent(mgs)*xdia(mgs,li,1) +
c     >   fmlt2(mgs)*(qiacr(mgs)+qiacw(mgs)) )
c     <   , 0.0 )
      qipmlr(mgs) = 0.0
c     >   min(
c     >  (fmlt1(mgs)*cx(mgs,lip)*cipvent(mgs)*xdia(mgs,lip,1) +
c     >   fmlt2(mgs)*(qipacr(mgs)+qipacw(mgs)) )
c     <   , 0.0 )
      qirmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lir)*cirvent(mgs)*xdia(mgs,lir,1) +
     >   fmlt2(mgs)*(qiracr(mgs)+qiracw(mgs)) )
     <   , 0.0 )
      qirmlw(mgs) = 0.0
      IF ( qx(mgs,lir) .lt. 1.0e-7 ) THEN
       qirmlw(mgs) = qirmlr(mgs)
       qirmlr(mgs) = 0.0
      ENDIF 
      qsmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,ls)*swvent(mgs)*xdia(mgs,ls,1) +
     >   fmlt2(mgs)*(qsacr(mgs)+qsacw(mgs)) )
     <   , 0.0 )
      qhmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1) +
     >   fmlt2(mgs)*(qhacr(mgs)+qhacw(mgs)) )
     <   , 0.0 )
      qhlmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lhl)*hlvent(mgs)*xdia(mgs,lhl,1) +
     >   fmlt2(mgs)*(qhlacr(mgs)+qhlacw(mgs)) )
     <   , 0.0 )
      qfmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lf)*fwvent(mgs)*xdia(mgs,lf,1) +
     >   fmlt2(mgs)*(qfacr(mgs)+qfacw(mgs)) )
     <   , 0.0 )
      qglmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lgl)*glvent(mgs)*xdia(mgs,lgl,1) +
     >   fmlt2(mgs)*(qglacr(mgs)+qglacw(mgs)) )
     <   , 0.0 )
      qgmmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lgm)*gmvent(mgs)*xdia(mgs,lgm,1) +
     >   fmlt2(mgs)*(qgmacr(mgs)+qgmacw(mgs)) )
     <   , 0.0 )
      qghmlr(mgs) =
     >   min(
     >  (fmlt1(mgs)*cx(mgs,lgh)*ghvent(mgs)*xdia(mgs,lgh,1) +
     >   fmlt2(mgs)*(qghacr(mgs)+qghacw(mgs)) )
     <   , 0.0 ) 
c
      qimlr(mgs)  = max( qimlr(mgs), -qimxd(mgs) ) 
      qirmlr(mgs) = max( qirmlr(mgs), -qirmxd(mgs) ) 
      qirmlw(mgs) = max( qirmlw(mgs), -qirmxd(mgs) )
      qipmlr(mgs) = max( qipmlr(mgs), -qipmxd(mgs) ) 
      qsmlr(mgs)  = max( qsmlr(mgs),  -qsmxd(mgs) ) 
      qhmlr(mgs)  = max( qhmlr(mgs),  -qhmxd(mgs) ) 
      qhlmlr(mgs) = max( qhlmlr(mgs), -qhlmxd(mgs) ) 
      qfmlr(mgs)  = max( qfmlr(mgs),  -qfmxd(mgs) ) 
      qglmlr(mgs) = max( qglmlr(mgs), -qglmxd(mgs) ) 
      qgmmlr(mgs) = max( qgmmlr(mgs), -qgmmxd(mgs) ) 
      qghmlr(mgs) = max( qghmlr(mgs), -qghmxd(mgs) ) 
c
      end do
c
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cimlr(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qimlr(mgs)
      csmlr(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qsmlr(mgs)
      cipmlr(mgs) = (cx(mgs,lip)/(qx(mgs,lip)+1.e-20))*qipmlr(mgs)
      cirmlr(mgs) = (cx(mgs,lir)/(qx(mgs,lir)+1.e-20))*qirmlr(mgs)
      cirmlw(mgs) = (cx(mgs,lir)/(qx(mgs,lir)+1.e-20))*qirmlw(mgs)

c
c ERM: what was I thinking here?  Looks like an attempt to calculate soaking of the 
c meltwater into the graupel before shedding excess?
c
c cglmlr is put into GM instead of into rain, and cgmmlr is put into GH instead of rain
c
c      qglmlr(mgs) = qglmlr(mgs)*(1.0 + (xdn(mgs,lr) - xdn(mgs,lgm))/
c     :                 (xdn(mgs,lgm)-xdn(mgs,lgl)))
c      qgmmlr(mgs) = qgmmlr(mgs)*(1.0 + (xdn(mgs,lr) - xdn(mgs,lgh))/
c     :                 (xdn(mgs,lgh)-xdn(mgs,lgm)))

      IF ( ipconc .ge. 3 ) THEN
      
      cfmlrr(mgs)  = 0.0
      cglmlrr(mgs) = 0.0
      cgmmlrr(mgs) = 0.0
      cghmlrr(mgs) = 0.0
      chmlrr(mgs)  = 0.0
      chlmlrr(mgs) = 0.0
      
      
      cglmlr(mgs) = (cx(mgs,lgl)/(qx(mgs,lgl)+1.e-20))*qglmlr(mgs)
      cgmmlr(mgs) = (cx(mgs,lgm)/(qx(mgs,lgm)+1.e-20))*qgmmlr(mgs)
      cghmlr(mgs) = (cx(mgs,lgh)/(qx(mgs,lgh)+1.e-20))*qghmlr(mgs)
      cfmlr(mgs)  = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfmlr(mgs)
      chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
      chlmlr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)

      IF ( ihmlt .eq. 1 ) THEN
        cfmlrr(mgs)  = Max( cfmlr(mgs), rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cglmlrr(mgs) = Max( cglmlr(mgs), rho0(mgs)*qglmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cgmmlrr(mgs) = Max( cgmmlr(mgs), rho0(mgs)*qgmmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cghmlrr(mgs) = Max( cghmlr(mgs), rho0(mgs)*qghmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        chmlrr(mgs)  = Max( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        chlmlrr(mgs) = Max( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      ELSEIF ( ihmlt .eq. 2 ) THEN

        cfmlrr(mgs)  = Max( cfmlr(mgs), rho0(mgs)*qfmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cglmlrr(mgs) = Max( cglmlr(mgs), rho0(mgs)*qglmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cgmmlrr(mgs) = Max( cgmmlr(mgs), rho0(mgs)*qgmmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        cghmlrr(mgs) = Max( cghmlr(mgs), rho0(mgs)*qghmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        chmlrr(mgs)  = Max( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
        chlmlrr(mgs) = Max( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 

        IF ( xv(mgs,lf ) .gt. 0.0 ) cfmlrr(mgs)  = 
     :       Max( cfmlrr(mgs), rho0(mgs)*qfmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lf)*xv(mgs,lf))) ) ! 
        IF ( xv(mgs,lgl) .gt. 0.0 ) cglmlrr(mgs)  =
     :        Max( cglmlrr(mgs), rho0(mgs)*qglmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lgl)*xv(mgs,lgl))) ) ! 
        IF ( xv(mgs,lgm) .gt. 0.0 ) cgmmlrr(mgs)  = 
     :        Max( cgmmlrr(mgs), rho0(mgs)*qgmmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lgm)*xv(mgs,lgm))) ) ! 
        IF ( xv(mgs,lgh) .gt. 0.0 ) cghmlrr(mgs)  = 
     :        Max( cghmlrr(mgs), rho0(mgs)*qghmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lgh)*xv(mgs,lgh))) ) ! 
        IF ( xv(mgs,lh ) .gt. 0.0 ) chmlrr(mgs)  = 
     :        Max( chmlrr(mgs), rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh))) ) ! 
        IF ( xv(mgs,lhl) .gt. 0.0 ) chlmlrr(mgs)  = 
     :        Max( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl))) ) ! 

      ELSEIF ( ihmlt .eq. 0 ) THEN
      cfmlrr(mgs)  = cfmlr(mgs)
      cglmlrr(mgs) = cglmlr(mgs)
      cgmmlrr(mgs) = cgmmlr(mgs)
      cghmlrr(mgs) = cghmlr(mgs)
      chmlrr(mgs)  = chmlr(mgs)
      chlmlrr(mgs) = chlmlr(mgs)
      ENDIF
      
      ENDIF ! ipconc .ge. 3

      end do
      end if  
      
! vr1mm
c
c  deposition/sublimation of ice
c
      DO mgs = 1,ngscnt
      
      circap(mgs) = (0.5)*xdia(mgs,lir,1)
      fwcap(mgs) = (0.5)*xdia(mgs,lf,1)
      rwcap(mgs) = (0.5)*xdia(mgs,lr,1)
      swcap(mgs) = (0.5)*xdia(mgs,ls,1)
      glcap(mgs) = (0.5)*xdia(mgs,lgl,1)
      gmcap(mgs) = (0.5)*xdia(mgs,lgm,1)
      ghcap(mgs) = (0.5)*xdia(mgs,lgh,1)
      hwcap(mgs) = (0.5)*xdia(mgs,lh,1)
      hlcap(mgs) = (0.5)*xdia(mgs,lhl,1)

      if ( qx(mgs,li).gt.qxmin(li) .and. xdia(mgs,li,1) .gt. 0.0 ) then
c
c from Cotton, 1972 (Part II)
c
        cilen(mgs)   = 0.4764*(xdia(mgs,li,1))**(0.958)
        cval = xdia(mgs,li,1)
        aval = cilen(mgs)
        eval = Sqrt(1.0-(aval**2)/(cval**2))
        fval = min(0.99,eval)
        gval = alog( abs( (1.+fval)/(1.-fval) ) )
        cicap(mgs) = cval*fval / gval
      ELSE
       cicap(mgs) = 0.0
      end if
      cipcap(mgs) = xdia(mgs,lip,1)*pii
      ENDDO
c
c
      do mgs = 1,ngscnt
      qidsv(mgs) =
     >   fvds(mgs)*cx(mgs,li)*civent(mgs)*cicap(mgs)
      qirdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lir)*cirvent(mgs)*circap(mgs)
      qipdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lip)*cipvent(mgs)*cipcap(mgs)
      qsdsv(mgs) =
     >   fvds(mgs)*cx(mgs,ls)*swvent(mgs)*swcap(mgs)
      qhdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lh)*hwvent(mgs)*hwcap(mgs)
      qhldsv(mgs) =
     >   fvds(mgs)*cx(mgs,lhl)*hlvent(mgs)*hlcap(mgs)
      qfdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lf)*fwvent(mgs)*fwcap(mgs)
      qgldsv(mgs) =
     >   fvds(mgs)*cx(mgs,lgl)*glvent(mgs)*glcap(mgs)
      qgmdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lgm)*gmvent(mgs)*gmcap(mgs)
      qghdsv(mgs) =
     >   fvds(mgs)*cx(mgs,lgh)*ghvent(mgs)*ghcap(mgs)
c
      end do
c
      do mgs = 1,ngscnt
      qisbv(mgs) = max( min(qidsv(mgs), 0.0), -qimxd(mgs) )
      qirsbv(mgs) = max( min(qirdsv(mgs), 0.0), -qirmxd(mgs) )
      qipsbv(mgs) = max( min(qipdsv(mgs), 0.0), -qipmxd(mgs) )
      qglsbv(mgs) = max( min(qgldsv(mgs), 0.0), -qglmxd(mgs) )
      qgmsbv(mgs) = max( min(qgmdsv(mgs), 0.0), -qgmmxd(mgs) )
      qghsbv(mgs) = max( min(qghdsv(mgs), 0.0), -qghmxd(mgs) )
      qfsbv(mgs) = max( min(qfdsv(mgs), 0.0), -qfmxd(mgs) )
      qhsbv(mgs) = max( min(qhdsv(mgs), 0.0), -qhmxd(mgs) )
      qhlsbv(mgs) = max( min(qhldsv(mgs), 0.0), -qhlmxd(mgs) )
      qssbv(mgs) = max( min(qsdsv(mgs), 0.0), -qsmxd(mgs) )

      qidpv(mgs) = Max(qidsv(mgs), 0.0)
      qsdpv(mgs) = Max(qsdsv(mgs), 0.0)
      qirdpv(mgs) = Max(qirdsv(mgs), 0.0)
      qipdpv(mgs) = Max(qipdsv(mgs), 0.0)
      qgldpv(mgs) = Max(qgldsv(mgs), 0.0)
      qgmdpv(mgs) = Max(qgmdsv(mgs), 0.0)
      qghdpv(mgs) = Max(qghdsv(mgs), 0.0)
      qfdpv(mgs) = Max(qfdsv(mgs), 0.0)
      qhdpv(mgs) = Max(qhdsv(mgs), 0.0)
      qhldpv(mgs) = Max(qhldsv(mgs), 0.0)
      
      temp1 = qidpv(mgs) + qsdpv(mgs) + qirdpv(mgs) +
     :        qipdpv(mgs) + qgldpv(mgs) + qgmdpv(mgs) +
     :        qghdpv(mgs) + qfdpv(mgs) + qhldpv(mgs) + qhdpv(mgs)
      
      IF ( temp1 .gt. qvimxd(mgs) ) THEN
      
      frac = qvimxd(mgs)/temp1
      
      qidpv(mgs) = frac*qidpv(mgs)
      qsdpv(mgs) = frac*qsdpv(mgs)
      qirdpv(mgs) = frac*qirdpv(mgs)
      qipdpv(mgs) = frac*qipdpv(mgs)
      qgldpv(mgs) = frac*qgldpv(mgs)
      qgmdpv(mgs) = frac*qgmdpv(mgs)
      qghdpv(mgs) = frac*qghdpv(mgs)
      qfdpv(mgs) = frac*qfdpv(mgs)
      qhdpv(mgs) = frac*qhdpv(mgs)
      qhldpv(mgs) = frac*qhldpv(mgs)
      
      ENDIF

      end do
c
c
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cssbv(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qssbv(mgs)
      cisbv(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qisbv(mgs)
      cirsbv(mgs) = (cx(mgs,lir)/(qx(mgs,lir)+1.e-20))*qirsbv(mgs)
      cipsbv(mgs) = (cx(mgs,lip)/(qx(mgs,lip)+1.e-20))*qipsbv(mgs)
      cglsbv(mgs) = (cx(mgs,lgl)/(qx(mgs,lgl)+1.e-20))*qglsbv(mgs)
      cgmsbv(mgs) = (cx(mgs,lgm)/(qx(mgs,lgm)+1.e-20))*qgmsbv(mgs)
      cghsbv(mgs) = (cx(mgs,lgh)/(qx(mgs,lgh)+1.e-20))*qghsbv(mgs)
      cfsbv(mgs)  = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfsbv(mgs)
      chsbv(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhsbv(mgs)
      chlsbv(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlsbv(mgs)
      csdpv(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qsdpv(mgs)
      cidpv(mgs) =  0.0 ! (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qidpv(mgs)
      cirdpv(mgs) = 0.0 ! (cx(mgs,lir)/(qx(mgs,lir)+1.e-20))*qirdpv(mgs)
      cipdpv(mgs) = 0.0 ! (cx(mgs,lip)/(qx(mgs,lip)+1.e-20))*qipdpv(mgs)
      cgldpv(mgs) = 0.0 ! (cx(mgs,lgl)/(qx(mgs,lgl)+1.e-20))*qgldpv(mgs)
      cgmdpv(mgs) = 0.0 ! (cx(mgs,lgm)/(qx(mgs,lgm)+1.e-20))*qgmdpv(mgs)
      cghdpv(mgs) = 0.0 ! (cx(mgs,lgh)/(qx(mgs,lgh)+1.e-20))*qghdpv(mgs)
      cfdpv(mgs)  = 0.0 ! (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfdpv(mgs)
      chdpv(mgs)  = 0.0 ! (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhdpv(mgs)
      end do
      end if
c
c
c  compute dry growth rate of snow, graupel, and hail
c
      do mgs = 1,ngscnt
c
      qsdry(mgs)  = qsacr(mgs)    + qsacw(mgs)
     >            + qsaci(mgs)    + qsacir(mgs)
     >            + qsacip(mgs)
c
      qfdry(mgs)  = qfaci(mgs)    + qfacs(mgs)
     >            + qfacr(mgs)    + qfacir(mgs)  
     >            + qfacw(mgs)    + qfacip(mgs)   
c
      qhdry(mgs)  = qhaci(mgs)    + qhacs(mgs)
     >            + qhacr(mgs)    + qhacir(mgs)   
     >            + qhacw(mgs)    + qhacip(mgs)
c
      qhldry(mgs) = qhlaci(mgs)   + qhlacs(mgs)
     >            + qhlacr(mgs)   + qhlacir(mgs)
     >            + qhlacw(mgs)   + qhlacip(mgs)
c
      qgldry(mgs) = qglaci(mgs)   + qglacs(mgs)
     >            + qglacr(mgs)   + qglacir(mgs)
     >            + qglacw(mgs)   + qglacip(mgs)
c
      qgmdry(mgs) = qgmaci(mgs)   + qgmacs(mgs)
     >            + qgmacr(mgs)   + qgmacir(mgs)
     >            + qgmacw(mgs)   + qgmacip(mgs)
c
      qghdry(mgs) = qghaci(mgs)   + qghacs(mgs)
     >            + qghacr(mgs)   + qghacir(mgs)
     >            + qghacw(mgs)   + qghacip(mgs)
c
      end do
c
c  set wet growth and shedding
c
      do mgs = 1,ngscnt
c
      qswet(mgs) =
     >  ( xdia(mgs,ls,1)*swvent(mgs)*cx(mgs,ls)*fwet1(mgs)
     >  + fwet2(mgs)*(qsaci(mgs)+qsacir(mgs)
     >               +qsacip(mgs)) )
      qswet(mgs) = max( 0.0, qswet(mgs))
c
      qfwet(mgs) =
     >  ( xdia(mgs,lf,1)*fwvent(mgs)*cx(mgs,lf)*fwet1(mgs)
     >  + fwet2(mgs)*(qfaci(mgs)+qfacir(mgs)
     >               +qfacip(mgs)+qfacs(mgs)) )
      qfwet(mgs) = max( 0.0, qfwet(mgs))
c
      qhwet(mgs) =
     >  ( xdia(mgs,lh,1)*hwvent(mgs)*cx(mgs,lh)*fwet1(mgs)
     >  + fwet2(mgs)*(qhaci(mgs)+qhacir(mgs)
     >               +qhacip(mgs)+qhacs(mgs)) )
      qhwet(mgs) = max( 0.0, qhwet(mgs))
c
      qhlwet(mgs) =
     >  ( xdia(mgs,lhl,1)*hlvent(mgs)*cx(mgs,lhl)*fwet1(mgs)
     >  + fwet2(mgs)*(qhlaci(mgs)+qhlacir(mgs)
     >               +qhlacip(mgs)+qhlacs(mgs)) )
      qhlwet(mgs) = max( 0.0, qhlwet(mgs))
c
      qglwet(mgs) =
     >  ( xdia(mgs,lgl,1)*glvent(mgs)*cx(mgs,lgl)*fwet1(mgs)
     >  + fwet2(mgs)*(qglaci(mgs)+qglacir(mgs)
     >               +qglacip(mgs)+qglacs(mgs)) )
      qglwet(mgs) = max( 0.0, qglwet(mgs))
c
      qgmwet(mgs) =
     >  ( xdia(mgs,lgm,1)*gmvent(mgs)*cx(mgs,lgm)*fwet1(mgs)
     >  + fwet2(mgs)*(qgmaci(mgs)+qgmacir(mgs)
     >               +qgmacip(mgs)+qgmacs(mgs)) )
      qgmwet(mgs) = max( 0.0, qgmwet(mgs))
c
      qghwet(mgs) =
     >  ( xdia(mgs,lgh,1)*ghvent(mgs)*cx(mgs,lgh)*fwet1(mgs)
     >  + fwet2(mgs)*(qghaci(mgs)+qghacir(mgs)
     >               +qghacip(mgs)+qghacs(mgs)) )
      qghwet(mgs) = max( 0.0, qghwet(mgs)) 
c
      end do
c
c shedding rate
c
      do mgs = 1,ngscnt
c
      qsshr(mgs)  =  0.0
      qfshr(mgs)  =  0.0
      qhshr(mgs)  =  0.0
      qhlshr(mgs)  =  0.0
      qglshr(mgs)  =  0.0
      qgmshr(mgs)  =  0.0
      qghshr(mgs)  =  0.0
c
c      qsshr(mgs)  = 
c     >   qswet(mgs) - qsdry(mgs) 
c
      qfshr(mgs)  = 
     >   Min( 0.0, qfwet(mgs) - qfdry(mgs) )
c
      qhshr(mgs)  = 
     >    Min( 0.0, qhwet(mgs) - qhdry(mgs) )
c
      qhlshr(mgs)  = 
     >    Min( 0.0, qhlwet(mgs) - qhldry(mgs) )
c
c      qglshr(mgs)  = 
c     >   qglwet(mgs) - qgldry(mgs) 
c
c      qgmshr(mgs)  = 
c     >   qgmwet(mgs) - qgmdry(mgs) 
c
      qghshr(mgs)  = 
     >    Min( 0.0, qghwet(mgs) - qghdry(mgs) )
c
c limit wet growth to only higher density particles
c
      qsshr(mgs)  =  0.0
      qglshr(mgs)  =  0.0
      qgmshr(mgs)  =  0.0
c
c
c  if T<243 then set shedding to zero....
c
      if ( temg(mgs) .lt. 243.15 ) then
      qsshr(mgs)  =  0.0
      qfshr(mgs)  =  0.0
      qhshr(mgs)  =  0.0
      qhlshr(mgs)  =  0.0
      qglshr(mgs)  =  0.0
      qgmshr(mgs)  =  0.0
      qghshr(mgs)  =  0.0
      end if
c
c  shed all at temperatures > 273.15
c
      if ( temg(mgs) .gt. tfr ) then
      qsshr(mgs)  = -qsdry(mgs)
      qhshr(mgs)  = -qhdry(mgs)
      qhlshr(mgs) = -qhldry(mgs)
      qfshr(mgs)  = -qfdry(mgs)
      qglshr(mgs) = -qgldry(mgs)
      qgmshr(mgs) = -qgmdry(mgs)
      qghshr(mgs) = -qghdry(mgs)
      end if
c
      end do
c
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      csshr(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*Min(0.0,qsshr(mgs))
      cirshr(mgs) = 0.0 ! (cx(mgs,lir)/(qx(mgs,lir)+1.e-20))*qirshr(mgs)
      cipshr(mgs) = 0.0 ! (cx(mgs,lip)/(qx(mgs,lip)+1.e-20))*qipshr(mgs)
      cglshr(mgs) = 0.0 ! (cx(mgs,lgl)/(qx(mgs,lgl)+1.e-20))*qglshr(mgs)
      cgmshr(mgs) = 0.0 ! (cx(mgs,lgm)/(qx(mgs,lgm)+1.e-20))*qgmshr(mgs)
      cghshr(mgs) = (cx(mgs,lgh)/(qx(mgs,lgh)+1.e-20))*qghshr(mgs)
      cfshr(mgs)  = (cx(mgs,lf)/(qx(mgs,lf)+1.e-20))*qfshr(mgs)
      chshr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhshr(mgs)
      chlshr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlshr(mgs)
      end do
      end if
c
c  liquid accretted before wet grow shedding + eps (1e-30)
c
      do mgs = 1,ngscnt
      qiacl(mgs) = qiacw(mgs) + qiacr(mgs) + 1.e-20
      qsacl(mgs) = qsacw(mgs) + qsacr(mgs) + 1.e-20
      qfacl(mgs) = qfacw(mgs) + qfacr(mgs) + 1.e-20
      qhacl(mgs) = qhacw(mgs) + qhacr(mgs) + 1.e-20
      qhlacl(mgs) = qhlacw(mgs) + qhlacr(mgs) + 1.e-20
      qglacl(mgs) = qglacw(mgs) + qglacr(mgs) + 1.e-20
      qgmacl(mgs) = qgmacw(mgs) + qgmacr(mgs) + 1.e-20
      qghacl(mgs) = qghacw(mgs) + qghacr(mgs) + 1.e-20
      end do


c count up number of times that wet growth occurs:
      DO mgs = 1,ngscnt

      IF ( (qghshr(mgs) .lt. 0.0 .or.
     :      qhshr(mgs) .lt. 0.0 .or.
     :      qfshr(mgs) .lt. 0.0 .or.
     :      qhlshr(mgs) .lt. 0.0 ) .and. 
     :     temg(mgs) .lt. 273.15 ) THEN

c C$PAR CRITICAL SECTION (ilock4)
       iwetg = iwetg + 1
c C$PAR END CRITICAL SECTION

      END IF

      END DO
c
c  final decisions
c
      do mgs = 1,ngscnt
c
c  Snow
c
      if ( qsshr(mgs) .lt. 0.0 ) then
      qsdpv(mgs) = 0.0
      qssbv(mgs) = 0.0
      else
      qsshr(mgs) = 0.0
      end if
c
c     if ( qsdry(mgs) .lt. qswet(mgs) ) then
c     qswet(mgs) = 0.0
c     else
c     qsdry(mgs) = 0.0
c     end if 
c
c  hail
c
      if ( qhshr(mgs) .lt. 0.0 ) then
      qhdpv(mgs) = 0.0
      qhsbv(mgs) = 0.0

      schacs(mgs) = 0.0
      schacir(mgs) = 0.0
      schacip(mgs) = 0.0
      schaci(mgs) = 0.0
      
      schacw(mgs) = 0.0

      IF ( ehi(mgs) .gt. 0.0 ) THEN
        qhaci(mgs) = Min(qimxd(mgs),qhaci(mgs)/ehi(mgs))
      ENDIF
      IF ( ehip(mgs) .gt. 0.0 ) THEN
        qhacip(mgs) = Min(qipmxd(mgs),qhacip(mgs)/ehip(mgs))
      ENDIF
      IF ( ehir(mgs) .gt. 0.0 ) THEN
        qhacir(mgs) = Min(qirmxd(mgs),qhacir(mgs)/ehir(mgs))
      ENDIF
      IF ( ehs(mgs) .gt. 0.0 ) THEN
        qhacs(mgs) = Min(qsmxd(mgs),qhacs(mgs)/ehs(mgs))
      ENDIF
      
      qhwet(mgs) = 1.0


      else
      qhshr(mgs) = 0.0
      qhwet(mgs) = 0.0
      end if
c
c  hail - large
c
      if ( qhlshr(mgs) .lt. 0.0 ) then
      qhldpv(mgs) = 0.0
      qhlsbv(mgs) = 0.0

      schlacs(mgs) = 0.0
      schlacir(mgs) = 0.0
      schlacip(mgs) = 0.0
      schlaci(mgs) = 0.0
      
      schlacw(mgs) = 0.0

      IF ( ehli(mgs) .gt. 0.0 ) THEN
        qhlaci(mgs) = Min(qimxd(mgs),qhlaci(mgs)/ehli(mgs))
      ENDIF
      IF ( ehlip(mgs) .gt. 0.0 ) THEN
        qhlacip(mgs) = Min(qipmxd(mgs),qhlacip(mgs)/ehlip(mgs))
      ENDIF
      IF ( ehlir(mgs) .gt. 0.0 ) THEN
        qhlacir(mgs) = Min(qirmxd(mgs),qhlacir(mgs)/ehlir(mgs))
      ENDIF
      IF ( ehls(mgs) .gt. 0.0 ) THEN
        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs(mgs)/ehls(mgs))
      ENDIF

      qhlwet(mgs) = 1.0
      
      else
      qhlshr(mgs) = 0.0
      qhlwet(mgs) = 0.0
      end if
c
c     if ( qhdry(mgs) .lt. qhwet(mgs) ) then
c     qhwet(mgs) = 0.0
c     else
c     qhdry(mgs) = 0.0
c     end if 
c
c  frozen drops
c
      if ( qfshr(mgs) .lt. 0.0 ) then
      qfdpv(mgs) = 0.0
      qfsbv(mgs) = 0.0

      scfacs(mgs) = 0.0
      scfacir(mgs) = 0.0
      scfacip(mgs) = 0.0
      scfaci(mgs) = 0.0
      
      scfacw(mgs) = 0.0

      IF ( efi(mgs) .gt. 0.0 ) THEN
        qfaci(mgs) = Min(qimxd(mgs),qfaci(mgs)/efi(mgs))
      ENDIF
      IF ( efip(mgs) .gt. 0.0 ) THEN
        qfacip(mgs) = Min(qipmxd(mgs),qfacip(mgs)/efip(mgs))
      ENDIF
      IF ( efir(mgs) .gt. 0.0 ) THEN
        qfacir(mgs) = Min(qirmxd(mgs),qfacir(mgs)/efir(mgs))
      ENDIF
      IF ( efs(mgs) .gt. 0.0 ) THEN
        qfacs(mgs) = Min(qsmxd(mgs),qfacs(mgs)/efs(mgs))
      ENDIF
      
      qfwet(mgs) = 1.0
      
      else
      qfshr(mgs) = 0.0
      qfwet(mgs) = 0.0
      end if
c
c     if ( qfdry(mgs) .lt. qfwet(mgs) ) then
c     qfwet(mgs) = 0.0
c     else
c     qfdry(mgs) = 0.0
c     end if 
c
c  graupel-l
c
      if ( qglshr(mgs) .lt. 0.0 ) then
      qgldpv(mgs) = 0.0
      qglsbv(mgs) = 0.0
      else
      qglshr(mgs) = 0.0
      end if
c
c     if ( qgldry(mgs) .lt. qglwet(mgs) ) then
c     qglwet(mgs) = 0.0
c     else
c     qgldry(mgs) = 0.0
c     end if 
c
c  graupel-m
c
      if ( qgmshr(mgs) .lt. 0.0 ) then
      qgmdpv(mgs) = 0.0
      qgmsbv(mgs) = 0.0
      else
      qgmshr(mgs) = 0.0
      end if
c
c     if ( qgmdry(mgs) .lt. qgmwet(mgs) ) then
c     qgmwet(mgs) = 0.0
c     else
c     qgmdry(mgs) = 0.0
c     end if 
c
c  graupel-h
c
      if ( qghshr(mgs) .lt. 0.0 ) then
      qghdpv(mgs) = 0.0
      qghsbv(mgs) = 0.0

      scghacs(mgs) = 0.0
      scghacir(mgs) = 0.0
      scghacip(mgs) = 0.0
      scghaci(mgs) = 0.0
      
      scghacw(mgs) = 0.0
      
      IF ( eghi(mgs) .gt. 0.0 ) THEN
        qghaci(mgs) = Min(qimxd(mgs),qghaci(mgs)/eghi(mgs))
      ENDIF
      IF ( eghip(mgs) .gt. 0.0 ) THEN
        qghacip(mgs) = Min(qipmxd(mgs),qghacip(mgs)/eghip(mgs))
      ENDIF
      IF ( eghir(mgs) .gt. 0.0 ) THEN
        qghacir(mgs) = Min(qirmxd(mgs),qghacir(mgs)/eghir(mgs))
      ENDIF
      
      qghwet(mgs) = 1.0
      
      else
      qghshr(mgs) = 0.0
      qghwet(mgs) = 0.0
      end if
c
c     if ( qghdry(mgs) .lt. qghwet(mgs) ) then
c     qghwet(mgs) = 0.0
c     else
c     qghdry(mgs) = 0.0
c     end if 
c
      end do
c
c
c
c  heat budget for rain---not all rain that collects ice can freeze
c
c
c
      if ( irwfrz .gt. 0 ) then
c
      do mgs = 1,ngscnt
c
c  compute total rain that freeze when it interacts with cloud ice
c
      qrztot(mgs) = 
     > + qrfrz(mgs) + qiacr(mgs) + qipacr(mgs) + qiracr(mgs)
c
c  compute the maximum amount of rain that can freeze
c  The four is from four source terms
c
      qrzmax(mgs) =
     >  ( xdia(mgs,lr,1)*rwvent(mgs)*cx(mgs,lr)*fwet1(mgs) )
      qrzmax(mgs) = max(qrzmax(mgs), 0.0)
      qrzmax(mgs) = min(4.*qrmxd(mgs), qrzmax(mgs))
c
c  compute the correction factor
c
      IF ( qrztot(mgs) .gt. 0.5*qxmin(lr) ) THEN
        qrzfac(mgs) = qrzmax(mgs)/(qrztot(mgs)+1.0e-20)
      ELSE
        qrzfac(mgs) = 1.0
      ENDIF
      qrzfac(mgs) = min(1.0, qrzfac(mgs))
c
      end do
c
c
c now correct the above sources
c
c
      do mgs = 1,ngscnt
      if ( temg(mgs) .le. 273.15 ) then
      qrfrz(mgs)   = qrzfac(mgs)*qrfrz(mgs)
      qiacr(mgs)   = qrzfac(mgs)*qiacr(mgs)
      qiracr(mgs)  = qrzfac(mgs)*qiracr(mgs)
      qipacr(mgs)  = qrzfac(mgs)*qipacr(mgs)
      end if
      end do
c
c
c      if (  nstep .eq. 360 ) then
c      do mgs = 1,ngscnt
c      write (69,969) 'QRFZ(TQTXF)',temg(mgs),qx(mgs,lr),qrztot(mgs)
c     > ,qrzmax(mgs),qrzfac(mgs)
c 969  format(1x,a4,1x,5f14.9)
c      end do
c      end if
c
      end if
c
c
c
c  evaporation of rain
c
c
c
      do mgs = 1,ngscnt
c
      qrcev(mgs) =
     >  fvce(mgs)*cx(mgs,lr)*rwvent(mgs)*rwcap(mgs)
! this line to allow condensation on rain:
      IF ( rcond .eq. 1 ) THEN
        qrcev(mgs) = min(qrcev(mgs), qrmxd(mgs))
! this line to have evaporation only:
      ELSE
        qrcev(mgs) = min(qrcev(mgs), 0.0)
      ENDIF

      qrcev(mgs) = max(qrcev(mgs), -qrmxd(mgs))
c      crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)+1.e-20))*qrcev(mgs)
      if ( temg(mgs) .lt. 273.15 ) qrcev(mgs) = 0.0
      crcev(mgs) = (rho0(mgs)/(cx(mgs,lr)+1.e-20))*qrcev(mgs)
c      if ( temg(mgs) .lt. 273.15 ) crcev(mgs) = 0.0
c  
      end do
c
c
c
c  ICE MULTIPLICATION: Two modes (rimpa, and rimpb)
c  (following Cotton et al. 1986)
c
 
      do mgs = 1,ngscnt
 
      chmul1(mgs) =  0.0
      chlmul1(mgs) = 0.0
      cfmul1(mgs) =  0.0
      cghmul1(mgs) = 0.0
      cgmmul1(mgs) = 0.0
      cglmul1(mgs) = 0.0
      csmul1(mgs) = 0.0
c
      qhmul1(mgs) =  0.0
      qhlmul1(mgs) = 0.0
      qfmul1(mgs) =  0.0
      qghmul1(mgs) = 0.0
      qgmmul1(mgs) = 0.0
      qglmul1(mgs) = 0.0
      qsmul1(mgs) =  0.0

      IF ( (itype1 .ge. 1 .or. itype2 .ge. 1 )
     :              .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
      if ( temg(mgs) .ge. 265.16 .and. temg(mgs) .le. 271.16 ) then
       IF ( ipconc .ge. 0 ) THEN
        IF ( xv(mgs,lc) .gt. 0.0 .and. 
     :       ( qx(mgs,lgl) .gt. qxmin(lgl) .or. 
     :         qx(mgs,lgm) .gt. qxmin(lgm) .or.
     :         qx(mgs,lgh) .gt. qxmin(lgh) .or.
     :        (qx(mgs,lf) .gt. qxmin(lf) .and. itype2 .ge. 2 ) ) ) THEN
c
c  Ziegler et al. 1986 Hallett-Mossop process.  VSTAR = 7.23e-15 (vol of 12micron radius)
c
         ex1 = (1./250.)*Exp(-7.23e-15/xv(mgs,lc))
       IF ( itype2 .le. 2 ) THEN
         ft = Max(0.0,Min(1.0,-0.11*temcg(mgs)**2 - 1.1*temcg(mgs)-1.7))
       ELSE
        IF ( temg(mgs) .ge. 265.16 .and. temg(mgs) .le. 267.16 ) THEN
          ft = 0.5
        ELSEIF (temg(mgs) .ge. 267.16 .and. temg(mgs) .le. 269.16 ) THEN
          ft = 1.0
        ELSEIF (temg(mgs) .ge. 269.16 .and. temg(mgs) .le. 271.16 ) THEN
          ft = 0.5
        ELSE 
          ft = 0.0
        ENDIF
       ENDIF
c        rhoinv = 1./rho0(mgs)
c        DNSTAR = ex1*cglacw(mgs)
         IF ( qx(mgs,lgl) .gt. qxmin(lgl)  ) THEN
          cglmul1(mgs) = (ft*ex1*cglacw(mgs))
          qglmul1(mgs) = cimas0*cglmul1(mgs)*rhoinv(mgs)
         ENDIF
         IF ( qx(mgs,lgm) .gt. qxmin(lgm) ) THEN
          cgmmul1(mgs) = (ft*ex1*cgmacw(mgs))
          qgmmul1(mgs) = cimas0*cgmmul1(mgs)*rhoinv(mgs)
         ENDIF
         IF ( qx(mgs,lgh) .gt. qxmin(lgh) ) THEN
          cghmul1(mgs) = (ft*ex1*cghacw(mgs))
          qghmul1(mgs) = cimas0*cghmul1(mgs)*rhoinv(mgs)
         ENDIF
          IF ( qx(mgs,lf) .gt. qxmin(lf) .and. itype2 .ge. 2 ) THEN
            cfmul1(mgs) = (ft*ex1*cfacw(mgs))
            qfmul1(mgs) = cimas0*cfmul1(mgs)*rhoinv(mgs)
          ENDIF

        ENDIF ! xv(mgs,lc) .gt. 0.0 .and.

       ELSE ! ipconc .lt. 2
c
c  define the temperature function
c
      fimt1(mgs) = 0.0
c
c Cotton et al. (1986) version
c
      if ( temg(mgs) .ge. 268.16 .and. temg(mgs) .le. 270.16 ) then
        fimt1(mgs) = 1.0 -(temg(mgs)-268.16)/2.0
      elseif (temg(mgs) .le. 268.16 .and. temg(mgs) .ge. 265.16 ) then
        fimt1(mgs) = 1.0 +(temg(mgs)-268.16)/3.0
      ELSE 
        fimt1(mgs) = 0.0
      end if
c
c Ferrier (1994) version
c
      if ( temg(mgs) .ge. 265.16 .and. temg(mgs) .le. 267.16 ) then
        fimt1(mgs) = 0.5
      elseif (temg(mgs) .ge. 267.16 .and. temg(mgs) .le. 269.16 ) then
        fimt1(mgs) = 1.0
      elseif (temg(mgs) .ge. 269.16 .and. temg(mgs) .le. 271.16 ) then
        fimt1(mgs) = 0.5
      ELSE 
        fimt1(mgs) = 0.0
      end if
c
c
c   type I:  350 splinters are formed for every 1e-3 grams of cloud
c            water accreted by graupel/hail (note converted to MKS units)
c            3.5e+8 has units of 1/kg
c
      IF ( itype1 .ge. 1 ) THEN
       fimta(mgs) = (3.5e+08)*rho0(mgs)
      ELSE
       fimta(mgs) = 0.0
      ENDIF

c
c
c   type II:  1 splinter formed for every 250 cloud droplets larger than
c             24 micons in diameter (12 microns in radius) accreted by
c             graupel/hail
c
c
      fimt2(mgs) = 0.0
      xcwmas = xmas(mgs,lc) * 1000.
c
      IF ( itype2 .ge. 1 ) THEN
      if ( xcwmas.lt.1.26e-9 ) then
        fimt2(mgs) = 0.0
      end if
      if ( xcwmas .le. 3.55e-9 .and. xcwmas .ge. 1.26e-9 ) then
        fimt2(mgs) = (2.27)*alog(xcwmas) + 13.39
      end if
      if ( xcwmas .gt. 3.55e-9 ) then
        fimt2(mgs) = 1.0
      end if

      fimt2(mgs) = min(fimt2(mgs),1.0)
      fimt2(mgs) = max(fimt2(mgs),0.0)
      
      ENDIF
c
c     qhmul2 = 0.0
c     qsmul2 = 0.0
c
c     qhmul2 =
c    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qhacw(mgs)
c     qsmul2 =
c    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qsacw(mgs)
c
c      cimas0 = (1.0e-12)
c      cimas0 = 2.5e-10
      IF ( qhwet(mgs) .eq. 0.0 ) THEN
c      chmul1(mgs) =  fimt1(mgs)*(fimta(mgs) +
c     :                           (4.0e-03)*fimt2(mgs))*qhacw(mgs)
      ENDIF
c      chlmul1(mgs) = fimt1(mgs)*(fimta(mgs) +
c     :                           (4.0e-03)*fimt2(mgs))*qhlacw(mgs)
      IF ( qfwet(mgs) .eq. 0.0 .and. itype2 .ge. 2 ) THEN
       cfmul1(mgs) =  fimt1(mgs)*(fimta(mgs) +
     :                           (4.0e-03)*fimt2(mgs))*qfacw(mgs)
      ENDIF
      IF ( qghwet(mgs) .eq. 0.0 ) THEN
      cghmul1(mgs) = fimt1(mgs)*(fimta(mgs) +
     :                           (4.0e-03)*fimt2(mgs))*qghacw(mgs)
      ENDIF
      cgmmul1(mgs) = fimt1(mgs)*(fimta(mgs) +
     :                           (4.0e-03)*fimt2(mgs))*qgmacw(mgs)
      cglmul1(mgs) = fimt1(mgs)*(fimta(mgs) +
     :                           (4.0e-03)*fimt2(mgs))*qglacw(mgs)
c      csmul1(mgs) =  fimt1(mgs)*(fimta(mgs) +
c     :                           (4.0e-03)*fimt2(mgs))*qsacw(mgs)
c
      qhmul1(mgs) =  chmul1(mgs)*(cimas0/rho0(mgs))
      qhlmul1(mgs) = chlmul1(mgs)*(cimas0/rho0(mgs))
      qfmul1(mgs) =  cfmul1(mgs)*(cimas0/rho0(mgs))
      qghmul1(mgs) = cghmul1(mgs)*(cimas0/rho0(mgs))
      qgmmul1(mgs) = cgmmul1(mgs)*(cimas0/rho0(mgs))
      qglmul1(mgs) = cglmul1(mgs)*(cimas0/rho0(mgs))
      qsmul1(mgs) =  csmul1(mgs)*(cimas0/rho0(mgs))
c
      ENDIF ! ( ipconc .ge. 2 )
      
      IF ( qghshr(mgs) .lt. 0.0 ) THEN
       cghmul1(mgs) = 0.0
       qghmul1(mgs) = 0.0
      ENDIF
      IF ( qfshr(mgs) .lt. 0.0 ) THEN
       cfmul1(mgs) = 0.0
       qfmul1(mgs) = 0.0
      ENDIF
      IF ( qhshr(mgs) .lt. 0.0 ) THEN
       chmul1(mgs) = 0.0
       qhmul1(mgs) = 0.0
      ENDIF
      IF ( qhlshr(mgs) .lt. 0.0 ) THEN
       chlmul1(mgs) = 0.0
       qhlmul1(mgs) = 0.0
      ENDIF
      
      end if ! (in temperature range)
      
      ENDIF ! ( itype1 .eq. 1 .or. itype2 .eq. 1)
c
      end do
c
c
c
c     end if
c
c     end do
c
c
c
c 
c  vapor to pristine ice crystals   UP
c
c
c
c  compute the nucleation rate
c
c     do mgs = 1,ngscnt
c     idqis = 0
c     if ( ssi(mgs) .gt. 1.0 ) idqis = 1
c     fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
c     dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/
c    >  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
c     qidsvp(mgs) = dqisdt(mgs)
c     cnnt = min(cnit*exp(-temcg(mgs)*bta1),1.0e+09)
c     qiint(mgs) = 
c    >  il5(mgs)*idqis*(1.0/dtp)
c    <  *min((6.88e-13)*cnnt/rho0(mgs), 0.25*dqisdt(mgs)) 
c     end do
c
c  Meyers et al. (1992; JAS) and Ferrier (1994) primary ice nucleation
c
      cmassin = cimasn  ! 6.88e-13
      do mgs = 1,ngscnt
      qiint(mgs) = 0.0
      ciint(mgs) = 0.0
      qicicnt(mgs) = 0.0
      cicint(mgs) = 0.0
      qipipnt(mgs) = 0.0
      cipint(mgs) = 0.0
      if ( ( temg(mgs) .lt. 268.15 .or. ( imeyers5 .and. temg(mgs) .lt.  273.0) ) .and. 
     :    ciintmx .gt. (cx(mgs,li)+cx(mgs,lip)) ) then
      fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
      dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/
     >  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
      qidsvp(mgs) = dqisdt(mgs)
      idqis = 0
      if ( ssi(mgs) .gt. 1.0 ) idqis = 1 
      dzfacp = max( float(kgsp(mgs)-kgs(mgs)), 0.0 )
      dzfacm = max( float(kgs(mgs)-kgsm(mgs)), 0.0 )
      qiint(mgs) =
     >  idqis*il5(mgs)
     >  *(cmassin/rho0(mgs))
     >  *max(0.0,wvel(mgs))
     >  *max((cninp(mgs)-cninm(mgs)),0.0)*gz(kgs(mgs))
     <  /((dzfacp+dzfacm))
      qiint(mgs) = min(qiint(mgs), max(0.25*dqisdt(mgs),0.0)) 
      ciint(mgs) = qiint(mgs)*rho0(mgs)/cmassin
c
c limit new crystals so it does not increase the current concentration
c  above ciintmx (default value is 1000 per liter (1.e6 per m**3))
c
c      ciintmx = 1.e9
      IF ( ciint(mgs) .gt. (ciintmx - (cx(mgs,li)+cx(mgs,lip)))) THEN
        ciint(mgs) = Max(0.0, ciintmx - (cx(mgs,li)+cx(mgs,lip)) )
        qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
      ENDIF
c
      if ( xplate(mgs) .eq. 1 ) then
      qipipnt(mgs) = qiint(mgs)
      cipint(mgs) = ciint(mgs)
      end if
c
      if ( xcolmn(mgs) .eq. 1 ) then
      qicicnt(mgs) = qiint(mgs)
      cicint(mgs) = ciint(mgs)
      end if
c
c     qipipnt(mgs) = 0.0
c     qicicnt(mgs) = qiint(mgs)
c
      end if
      end do
c
c 

c
c  vapor to cloud droplets   UP
c
c  NOT USED YET???? (qcwdv and qcwev do not appear after this)
c      do mgs = 1,ngscnt
c      qcwdv(mgs) = 0.0
c      qcwev(mgs) = 0.0
c      fcinit(mgs) = (felv(mgs)**2)/(cp*rw)
c      dqwsdt(mgs) = (qx(mgs,lv)-qvs(mgs))/
c     >  (1.0 + fcinit(mgs)*qvs(mgs)/tsqr(mgs))
c      if ( (dqwsdt(mgs)+qx(mgs,lc)) .ge. 0.0 ) then
c      qcwdv(mgs) = dqwsdt(mgs)/dtp 
c      end if
c      if ( (dqwsdt(mgs)+qx(mgs,lc)) .lt. 0.0 ) then
c      qcwev(mgs) = qx(mgs,lc)/dtp 
c      end if
c      end do
c
c
      if ( ndebug .gt. 0 ) print*,'dbg=7.inf:  rime time calc'
c
c set the riming time source term and zero out rt** if no riming is
c occuring
c
c
c
c  cloud age source
c
      do mgs = 1,ngscnt
      prtcw(mgs) = 0.0
      if ( qx(mgs,lc) .ge. qcwmntim ) then
      prtcw(mgs) =  1.0
      end if
      if ( qx(mgs,lc) .lt. qcwmntim ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtcw(mgs) =  0.0
        ELSE
          prtcw(mgs) =  alfarim*rtcw(mgs)
        ENDIF
      end if
      end do
c
c  ice riming time source
c
      do mgs = 1,ngscnt
      prtci(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qiacw(mgs).ge.qrimmnc ) then
      prtci(mgs) =  1.0
      end if
      end if
      if ( qiacw(mgs) .lt. qrimmnc ) then
         IF ( alfarim .eq. 0.0 ) THEN
           rtci(mgs) =  0.0
         ELSE
           prtci(mgs) =  alfarim*rtci(mgs)
         ENDIF
      end if
      end do
c
c  ice riming time source
c
      do mgs = 1,ngscnt
      prtip(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qipacw(mgs).ge.qrimmnc ) then
      prtip(mgs) =  1.0
      end if
      end if
      if ( qipacw(mgs) .lt. qrimmnc ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtip(mgs) =  0.0
        ELSE
          prtip(mgs) =  alfarim*rtip(mgs)
        ENDIF
      end if
      end do
c
c  ice riming time source
c
      do mgs = 1,ngscnt
      prtir(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qiracw(mgs).ge.qrimmnc ) then
      prtir(mgs) =  1.0
      end if
      end if
      if ( qiracw(mgs) .lt. qrimmnc ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtir(mgs) =  0.0
        ELSE
          prtir(mgs) =  alfarim*rtir(mgs)
        ENDIF
      end if
      end do
c
c  snow riming time source
c
      do mgs = 1,ngscnt
      prtsw(mgs) = 0.0
      if ( qsacw(mgs) .ge. qrimmnp ) then
      if ( temg(mgs) .lt. 273.15 ) then
      prtsw(mgs) =  1.0
      end if
      end if
      if ( qsacw(mgs) .lt. qrimmnp ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtsw(mgs) =  0.0
        ELSE
          prtsw(mgs) =  alfarim*rtsw(mgs)
        ENDIF
      end if
      end do
c
c  graupel-low density riming time source
c
      do mgs = 1,ngscnt
      prtgl(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qglacw(mgs) .ge. qrimmnp ) then
      prtgl(mgs) =  1.0
      end if
      end if
      if ( qglacw(mgs) .lt. qrimmnp ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtgl(mgs) = 0.0
        ELSE
          prtgl(mgs) =  alfarim*rtgl(mgs)
        ENDIF
      end if
      end do
c
c  graupel-medium density riming time source
c
      do mgs = 1,ngscnt
      prtgm(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qgmacw(mgs) .ge. qrimmnp ) then
      prtgm(mgs) =  1.0
      end if
      end if
      if ( qgmacw(mgs) .lt. qrimmnp ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtgm(mgs) =  0.0
        ELSE
          prtgm(mgs) =  alfarim*rtgm(mgs)
        ENDIF
      end if
      end do
c
c  graupel-high density riming time source
c
      do mgs = 1,ngscnt
      prtgh(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qghacw(mgs) .ge. qrimmnp ) then
      prtgh(mgs) =  1.0
      end if
      end if
      if ( qghacw(mgs) .lt. qrimmnp ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtgh(mgs) =  0.0
        ELSE
          prtgh(mgs) =  alfarim*rtgh(mgs)
        ENDIF
      end if
      end do
c
c  frozen drop riming time source
c
      do mgs = 1,ngscnt
      prtfw(mgs) = 0.0
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qfacw(mgs) .ge. qrimmnp ) then
      prtfw(mgs) =  1.0
      end if
      end if
      if ( qfacw(mgs) .lt. qrimmnp ) then
        IF ( alfarim .eq. 0.0 ) THEN
          rtfw(mgs) =  0.0
        ELSE
          prtfw(mgs) =  alfarim*rtfw(mgs)
        ENDIF
      end if
      end do
c
c
c
      if (ndebug .gt. 0 ) print*,'dbg = 8'
c
c
      if (ndebug .gt. 0 ) print*,'Collection: set 3-component'
c
c  time for riming....
c
c     rimtim = 240.0
c     dtrim = rimtim
c     xacrtim  = 120.0
c     tranfr = 0.50
c     tranfw = 0.50
c
c  coefficients for riming
c
c     rimc1 = 300.00
c     rimc2 = 0.44
c
c 
c  zero som arrays
c
c
      do mgs = 1,ngscnt
      qrshr(mgs) = 0.0
!      qsshrp(mgs) = 0.0
!      qfshrp(mgs) = 0.0
!      qhshrp(mgs) = 0.0
!      qhlshrp(mgs) = 0.0
!      qghshrp(mgs) = 0.0
!      qgmshrp(mgs) = 0.0
!      qglshrp(mgs) = 0.0
      end do
c
c
c  first sum all of the shed rain
c
c
      do mgs = 1,ngscnt
      qrshr(mgs) = qsshr(mgs) + qfshr(mgs) 
     >  + qglshr(mgs) + qgmshr(mgs) + qghshr(mgs) + qhshr(mgs)
     >  + qhlshr(mgs)
      crshr(mgs) = cfshr(mgs) 
     >  + cghshr(mgs) + chshr(mgs)
     >  + chlshr(mgs)
     >  + cglshr(mgs) + cgmshr(mgs)
      IF ( ipconc .ge. 3 ) THEN
       crshr(mgs) = Max(crshr(mgs), rho0(mgs)*qrshr(mgs)/(xdn(mgs,lr)*vr1mm) )
      ENDIF
      
      end do 
c
c
c
c  cloud water - cloud ice interaction:  Multi-component transformation
c
c
c
      do mgs = 1,ngscnt
      qipacwip(mgs) = 0.0
      qwacipip(mgs) = 0.0
      qipacwir(mgs) = 0.0
      qwacipir(mgs) = 0.0
      qipacwgl(mgs) = 0.0
      qwacipgl(mgs) = 0.0
      cipacwip(mgs) = 0.0
      cwacipip(mgs) = 0.0
      cipacwir(mgs) = 0.0
      cwacipir(mgs) = 0.0
c      cipacwgl(mgs) = 0.0
      cwacipgl(mgs) = 0.0
c      qipacwgm(mgs) = 0.0
c      qwacipgm(mgs) = 0.0
c      qipacwgh(mgs) = 0.0
c      qwacipgh(mgs) = 0.0
      rimipmas(mgs) = 0.0
      rimipfrc(mgs) = 0.0
      end do
c
      do mgs = 1,ngscnt
      qiacwi(mgs) = 0.0
      qwacii(mgs) = 0.0
      qiacwir(mgs) = 0.0
      qwaciir(mgs) = 0.0
      qiacwgl(mgs) = 0.0
      qwacigl(mgs) = 0.0
c      ciacwi(mgs) = 0.0
      cwacii(mgs) = 0.0
c      ciacwir(mgs) = 0.0
      cwaciir(mgs) = 0.0
      ciacwgl(mgs) = 0.0
      cwacigl(mgs) = 0.0
c      qiacwgm(mgs) = 0.0
c      qwacigm(mgs) = 0.0
c      qiacwgh(mgs) = 0.0
c      qwacigh(mgs) = 0.0
      rimimas(mgs) = 0.0
      rimifrc(mgs) = 0.0
      end do
c
      do mgs = 1,ngscnt
c      qiracwgh(mgs) =  0.0
c      qwacirgh(mgs) =  0.0
c      qiracwgm(mgs) =  0.0
c      qwacirgm(mgs) =  0.0
      qiracwgl(mgs) =  0.0
      qwacirgl(mgs) =  0.0
      qiracwir(mgs) =  0.0
      qwacirir(mgs) =  0.0
      ciracwgl(mgs) =  0.0
      cwacirgl(mgs) =  0.0
c      ciracwir(mgs) =  0.0
      cwacirir(mgs) =  0.0
      end do
c
c  qci version
c
      IF ( iglcnvi .ge. 1 ) THEN
      
      do mgs = 1,ngscnt
      if ( qx(mgs,li) .gt. qxmin(li) .and. qiacw(mgs) .gt. 1.e-12 ) then
      rimimas(mgs) = rho0(mgs)*qiacw(mgs)*rtci(mgs)/cx(mgs,li)
      if ( (temg(mgs) .lt. tfr) ) then
       qiacwi(mgs)  = qiacw(mgs)
       if ( (rimimas(mgs) .gt. 0.1*xmas(mgs,li)) .and.
     >     (qiacw(mgs) .gt. 1.0e-12) ) then
        ! at this point, qiacwgl is qiacw0 in the paper, the amount of qiacw that goes to ir/gl
        qiacwgl(mgs) = max((qiacw(mgs)-qidpv(mgs)), 0.0)  ! note that qidpv is always >= 0.
        qiacwi(mgs)  = qiacw(mgs)-qiacwgl(mgs) ! portion of qiacw that stays in qi
        qwacigl(mgs) = qwaci(mgs)*qiacwgl(mgs)/qiacw(mgs) ! portion of qi that goes to qgl and qir
        qwacii(mgs)  = qwaci(mgs) - qwacigl(mgs)   ! qi that stays in qi
        cwacigl(mgs) = cx(mgs,li)*qwacigl(mgs)/qx(mgs,li) ! cwaci0 = conc. leaving qi
        cwacii(mgs) =  cwaci(mgs) - cwacigl(mgs) ! cx(mgs,li)*qwacii(mgs)/qx(mgs,li)  ! cwaci staying in qi
       else
        qwaci(mgs) = 0.0
        cwaci(mgs) = 0.0
       end if
      end if
      else 
      qiacwi(mgs) = 0.0
      qiacw(mgs) = 0.0
      qwacii(mgs) = 0.0
      qwaci(mgs) = 0.0
      cwaci(mgs) = 0.0
      end if
      end do
c
c qip version
c 
      do mgs = 1,ngscnt
      if ( qx(mgs,lip) .gt. qxmin(li) .and. qipacw(mgs) .gt. 1.e-12 ) then
      rimipmas(mgs) = rho0(mgs)*qipacw(mgs)*rtip(mgs)/cx(mgs,lip)
      if ( (temg(mgs) .lt. tfr) ) then
      qipacwip(mgs)  = qipacw(mgs)
      if ( (rimipmas(mgs) .gt. 0.1*xmas(mgs,lip)) .and.
     >     (qipacw(mgs) .gt. 1.0e-12) ) then
        qipacwgl(mgs)  = max((qipacw(mgs)-qipdpv(mgs)), 0.0)
        qipacwip(mgs)  = qipacw(mgs)-qipacwgl(mgs)
        qwacipgl(mgs)  = qwacip(mgs)*qipacwgl(mgs)/qipacw(mgs)
        qwacipip(mgs)  = qwacip(mgs)-qwacipgl(mgs)
        cwacipgl(mgs) = cx(mgs,lip)*qwacipgl(mgs)/qx(mgs,lip) ! cwacip0 = conc. leaving qip
        cwacipip(mgs) =  cwacip(mgs) - cwacipgl(mgs)  ! cwacip staying in qip
      else
        qwacip(mgs) = 0.0
        cwacip(mgs) = 0.0
      end if
      end if
      else 
      qipacwip(mgs) = 0.0
      qipacw(mgs) = 0.0
      qwacipip(mgs) = 0.0
      qwacip(mgs) = 0.0
      end if
      end do
c 
c  qir version
c
      do mgs = 1,ngscnt
      rimirmas(mgs) = 0.0
      rimirfrc(mgs) = 0.0
      if ( qx(mgs,lir) .gt. qxmin(li) .and. qiracw(mgs) .gt. 1.e-12 ) then
      rimirmas(mgs) = rho0(mgs)*qiracw(mgs)*rtir(mgs)/cx(mgs,lir)
      if ( (temg(mgs) .lt. tfr) ) then
      qiracwir(mgs)  = qiracw(mgs)
c
c  If rime mass is greater than 10% of the particle mass (and if the
c   riming rate is greater than 1.0e-12), then 
c   initiate conversion to low density graupel
c
      if ( (rimirmas(mgs) .gt. 0.1*xmas(mgs,lir)) .and.
     >     (qiracw(mgs) .gt. 1.0e-12) ) then
        qiracwgl(mgs)  = max((qiracw(mgs)-qirdpv(mgs)), 0.0)
        qiracwir(mgs)  = qiracw(mgs)-qiracwgl(mgs)
        qwacirgl(mgs)  = qwacir(mgs)*qiracwgl(mgs)/qiracw(mgs)
        qwacirir(mgs)  = qwacir(mgs)-qwacirgl(mgs)
        cwacirgl(mgs) = cx(mgs,lir)*qwacirgl(mgs)/qx(mgs,lir) ! cwacir0 = conc. leaving qir
        cwacirir(mgs) =  cwacir(mgs) - cwacirgl(mgs)  ! cwacir staying in qir
      else
        qwacir(mgs) = 0.0
        cwacir(mgs) = 0.0
      end if
      end if
      else 
      qiracwir(mgs)  = 0.0
      qiracir(mgs)  = 0.0
      qwacirir(mgs) = 0.0
      qwacir(mgs) = 0.0
      end if
      end do
      
      ELSE
       DO mgs = 1,ngscnt
        
        qiacwi(mgs)   = qiacw(mgs)
        qipacwip(mgs) = qipacw(mgs)
        qiracwir(mgs) = qiracw(mgs)
        
        qwaci(mgs)  = 0.0
        qwacip(mgs) = 0.0
        qwacir(mgs) = 0.0

        cwaci(mgs)  = 0.0
        cwacip(mgs) = 0.0
        cwacir(mgs) = 0.0
        
       ENDDO
      ENDIF
c
c  T > 273.15
c
      do mgs = 1,ngscnt
      if ( (temg(mgs) .gt. tfr) ) then
      qiracwir(mgs) = 0.0
      qiracwgl(mgs) = 0.0
      qiracw(mgs) = 0.0
      qwacirgl(mgs) = 0.0
      qwacirir(mgs) = 0.0
      qwacir(mgs) = 0.0
      end if
      end do
c
      do mgs = 1,ngscnt
      if ( (temg(mgs) .gt. tfr) ) then
      qiacwi(mgs)  = 0.0
      qiacwgl(mgs) = 0.0
      qiacw(mgs)   = 0.0
      qwacigl(mgs) = 0.0
      qwaci(mgs)   = 0.0
      qwacii(mgs)   = 0.0
      end if
      end do
c
      do mgs = 1,ngscnt
      if ( (temg(mgs) .gt. tfr) ) then
      qipacwip(mgs)  = 0.0
      qipacwgl(mgs) = 0.0
      qipacw(mgs)   = 0.0
      qwacipgl(mgs) = 0.0
      qwacip(mgs)   = 0.0
      qwacipip(mgs)   = 0.0
      end if
      end do
c
      IF ( iglcnvi .ge. 1 ) THEN

      do mgs = 1,ngscnt
      if ( qipacw(mgs) .gt. 1.e-12 ) then
      if ( rimipmas(mgs) .gt. 0.1*xmas(mgs,lip) ) then
      rimipfrc(mgs) = rimipmas(mgs)/(xmas(mgs,lip)+1.e-30)
      rimipfrc(mgs) = min(1.0,rimipfrc(mgs))
      qipacwir(mgs) = (1.0-rimipfrc(mgs))*qipacwgl(mgs)
      qwacipir(mgs) = (1.0-rimipfrc(mgs))*qwacipgl(mgs)
      qipacwgl(mgs) = rimipfrc(mgs)*qipacwgl(mgs)
      qwacipgl(mgs) = rimipfrc(mgs)*qwacipgl(mgs)
c I think that cipacwir/cipacwgl will double-count for cwacipir, so turned off 
c for now.  cipacwir would represent the number concentration of droplets 
c transfered to ir, but that is not meaningful because qiracw does not affect 
c the ice concentration.
c      cipacwir(mgs) = (1.0-rimipfrc(mgs))*cipacwgl(mgs)
      cwacipir(mgs) = (1.0-rimipfrc(mgs))*cwacipgl(mgs)
c      cipacwgl(mgs) = rimipfrc(mgs)*cipacwgl(mgs)
      cwacipgl(mgs) = rimipfrc(mgs)*cwacipgl(mgs)
      end if
      end if
      end do
c
      do mgs = 1,ngscnt
      if ( qiacw(mgs) .gt. 1.e-12 ) then
      if ( rimimas(mgs) .gt. 0.1*xmas(mgs,li) ) then
      rimifrc(mgs) = rimimas(mgs)/(xmas(mgs,li)+1.e-30)
      rimifrc(mgs) = min(1.0,rimifrc(mgs))
      qiacwir(mgs) = (1.0-rimifrc(mgs))*qiacwgl(mgs)
      qwaciir(mgs) = (1.0-rimifrc(mgs))*qwacigl(mgs)
      cwaciir(mgs) = (1.0-rimifrc(mgs))*cwacigl(mgs)
      qiacwgl(mgs) = rimifrc(mgs)*qiacwgl(mgs)
      qwacigl(mgs) = rimifrc(mgs)*qwacigl(mgs)
      cwacigl(mgs) = rimifrc(mgs)*cwacigl(mgs)
      end if
      end if
      end do
c
      do mgs = 1,ngscnt
      if ( rimirmas(mgs) .gt. 0.1*xmas(mgs,lir) ) then
      qiracwgl(mgs) = 1.00*qiracwgl(mgs)
      qwacirgl(mgs) = 1.00*qwacirgl(mgs)
      end if
      end do
      
      ENDIF
c
c
c  cloud water-snow interaction
c
c
c
      if ( ndebug .eq. 1 ) then
      write(iunit,*) 'cloud water-snow interaction'
      end if
      do mgs = 1,ngscnt
c
      qwacss(mgs)  = 0.0
      qsacws(mgs)  = 0.0
      qwacsgl(mgs) = 0.0
      qsacwgl(mgs) = 0.0
      qwacsgm(mgs) = 0.0
      qsacwgm(mgs) = 0.0
      qwacsgh(mgs) = 0.0
      qsacwgh(mgs) = 0.0
      qwacsf(mgs)  = 0.0
      qsacwf(mgs)  = 0.0

      IF ( ipconc .ge. 4 ) THEN
      cwacsgl(mgs) = 0.0
      cwacsgm(mgs) = 0.0
      cwacsgh(mgs) = 0.0
      cwacsf(mgs)  = 0.0
      ENDIF
      
      IF ( qx(mgs,ls) .gt. qxmin(ls) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
c      dtrim = Min(rimtim,rtsw(mgs))
c      cwmasr = qx(mgs,lc)*rho0(mgs)/max(cx(mgs,lc),1.e-20)
c      cwdiar = (6.*cwmasr/xdn(mgs,lc))**(1./3.)
c      rims(mgs) = 1000.0
c      if ( temg(mgs) .lt. 273.15 ) then
c      rims(mgs) = rimc1*(-((0.5)*(1.e+06)*cwdiar)
c     >                *((0.60)*vtxbar(mgs,ls,1))
c     >                /(temg(mgs)-273.15))**(rimc2)
c      rims(mgs) = min( max(100.0, rims(mgs)), 900.0)
c      end if
c
c      rimmas = qsacw(mgs)*dtrim
c
c      rhobar =
c     >  (qx(mgs,ls)*xdn(mgs,ls)+rimmas*rims(mgs)+eps)
c     > /(qx(mgs,ls)+rimmas+eps)
c      rhobar = min( max(100.0, rhobar), 900.0)
c
c
c
c          
      if ( temg(mgs) .le. tfr ) then
c     
c      if ( itimrimw .eq. 1 ) then
c      itim = itim0
c      end if
c
c      if ( itimrimw .eq. 2 ) then
c      argtim = ((rtsw(mgs)-0.5/xftim)*xftim + yftim)
c      itim = min(max(ifix(argtim),1), mtim)
c      end if
c
c      argqcw = (qx(mgs,lc) + 5.e-4 - qxmin(lc))*xfqcw+yfqcw
c      argqxw = (qx(mgs,ls) + 5.e-4 - qxmin(ls))*xfqxw+yfqxw
c      argtem = -temcg(mgs)*xftem+yftem
c      iqcw = min(max(ifix(argqcw),1), mqcw)
c      iqsw = min(max(ifix(argqxw),1), mqxw)
c      item = min(max(ifix(argtem),1), mtem)
c
c      subroutine rimconv(qcw,qfw,ccw,cfw,rtfw,cdx(lf),cno(lf),
c     :           fwdia,cwmas,cwdia,cwdia2,vtwbar,
c     :           rfr,intyp,temg,rho0)
c Note that there is no cdx(ls), and cdgl will not be used
c because intyp=1
      call rimconv(qx(mgs,lc),qx(mgs,ls),cx(mgs,lc),cx(mgs,ls),rtsw(mgs),
     :           cdx(lgl),cno(ls),xdia(mgs,ls,1),xmas(mgs,lc),xdia(mgs,lc,1),
     :           xdia(mgs,lc,2),vtxbar(mgs,lc,1),
     :           rfr,1,temg(mgs),rho0(mgs),cd,cs,ds,crfr)

      frcswsw = rfr(1) ! qfswacw(iqcw,iqsw,item,1,itim)
      frcswgl = rfr(2) ! qfswacw(iqcw,iqsw,item,2,itim)
      frcswgm = rfr(3) ! qfswacw(iqcw,iqsw,item,3,itim)
      frcswgh = rfr(4) + rfr(5) ! qfswacw(iqcw,iqsw,item,4,itim)
! do not allow direct conversion of snow to hail or frozen drops
      frcswfw = 0.0   !rfr(5) ! qfswacw(iqcw,iqsw,item,5,itim)
c
      frcswsw1 = (1.0-frcswsw)
c
      qsacws(mgs)  =   frcswsw*qsacw(mgs)
      qsacwgl(mgs) =   frcswgl*qsacw(mgs)
      qwacsgl(mgs) =   frcswgl*qwacs(mgs)
      qsacwgm(mgs) =   frcswgm*qsacw(mgs)
      qwacsgm(mgs) =   frcswgm*qwacs(mgs)
      qsacwgh(mgs) =   frcswgh*qsacw(mgs)
      qwacsgh(mgs) =   frcswgh*qwacs(mgs)
      qsacwf(mgs)  =   frcswfw*qsacw(mgs)
      qwacsf(mgs)  =   frcswfw*qwacs(mgs)
      qwacs(mgs)   =  frcswsw1*qwacs(mgs)
      
      IF ( ipconc .ge. 4 ) THEN
      IF ( icgxconv .eq. 0 ) THEN
        cwacsgl(mgs) = Min( csmxd(mgs), qwacsgl(mgs)*rho0(mgs)/xmas(mgs,ls))
        cwacsgm(mgs) = Min( csmxd(mgs), qwacsgm(mgs)*rho0(mgs)/xmas(mgs,ls))
        cwacsgh(mgs) = Min( csmxd(mgs), qwacsgh(mgs)*rho0(mgs)/xmas(mgs,ls))
        cwacsf(mgs)  = Min( csmxd(mgs), qwacsf(mgs)*rho0(mgs)/xmas(mgs,ls))
      ELSE
        cwacsgl(mgs) = crfr(2)*cwacs(mgs)
        cwacsgm(mgs) = crfr(3)*cwacs(mgs)
        cwacsgh(mgs) = (crfr(4) + crfr(5))*cwacs(mgs)
        cwacsf(mgs)  = 0.0
      ENDIF
      ENDIF
c
      end if
c
c
      if ( temg(mgs) .gt. tfr ) then
      qwacs(mgs) = 0.0
      qwacss(mgs) = qwacs(mgs)
      qsacws(mgs) = qsacw(mgs)
      end if
c
      ENDIF ! ( qsw .gt. qxmin(ls) .and. qcw .gt. qxmin(lc) )
      
      end do
c
c
c  cloud water-graupel-l interaction:  update
c
c
      if ( ndebug .eq. 1 ) then
      write(iunit,*)  'cloud water-graupel-l interaction'
      end if
      do mgs = 1,ngscnt
      
      qwacglgl(mgs) = 0.0
      qglacwgl(mgs) = 0.0
      qwacglgm(mgs) = 0.0
      qglacwgm(mgs) = 0.0
      qwacglgh(mgs) = 0.0
      qglacwgh(mgs) = 0.0
      qwacglf(mgs)  = 0.0
      qglacwf(mgs)  = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cwacglgm(mgs) = 0.0
      cwacglgh(mgs) = 0.0
      cwacglf(mgs)  = 0.0
      ENDIF
      
      IF ( qx(mgs,lgl) .gt. qxmin(lgl) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
c
c      dtrim = Min(rimtim,rtgl(mgs))
c      cwmasr = qx(mgs,lc)*rho0(mgs)/max(cx(mgs,lc),1.e-20)
c      cwdiar = (6.*cwmasr/xdn(mgs,lc))**(1./3.)
c      rimgl(mgs) = 1000.0
c      if ( temg(mgs) .lt. 273.15 ) then
c      rimgl(mgs) = rimc1*(-((0.5)*(1.e+06)*cwdiar)
c     >                *((0.60)*vtxbar(mgs,lgl,1))
c     >                /(temg(mgs)-273.15))**(rimc2)
c      rimgl(mgs) = min( max(100.0, rimgl(mgs)), 900.0)
c      end if
c
c      rimmas = qglacw(mgs)*dtrim
c
c      rhoglw(mgs) = 0.0
c      rhobar =
c     >  (qx(mgs,lgl)*xdn(mgs,lgl)+rimmas*rimgl(mgs)+eps)
c     > /(qx(mgs,lgl)+rimmas+eps)
c      rhobar = min( max(100.0, rhobar), 900.0)
c      rhoglw(mgs) = rhobar
c
c
      if ( temg(mgs) .le. tfr ) then
c
c      if ( itimrimw .eq. 1 ) then
c      itim = itim0
c      end if
c
c      if ( itimrimw .eq. 2 ) then
c      argtim = ((rtgl(mgs)-0.5/xftim)*xftim + yftim)
c      itim = min(max(ifix(argtim),1), mtim)
c      end if
c
c      argqcw = (qx(mgs,lc) + 5.e-4 - qxmin(lc))*xfqcw+yfqcw
c      argqxw = (qx(mgs,lgl) + 5.e-4 - qxmin(lgl))*xfqxw+yfqxw
c      argtem = -temcg(mgs)*xftem+yftem
c      iqcw = min(max(ifix(argqcw),1), mqcw)
c      iqgl = min(max(ifix(argqxw),1), mqxw)
c      item = min(max(ifix(argtem),1), mtem)
c
c
c      subroutine rimconv(qcw,qfw,ccw,cfw,rtfw,cdx(lf),cno(lf),
c     :           fwdia,cwmas,cwdia,cwdia2,vtwbar,
c     :           rfr,intyp,temg,rho0)
c Note that there is no cdx(ls), and cdx(lgl) will not be used
c because intyp=1
      call rimconv(qx(mgs,lc),qx(mgs,lgl),cx(mgs,lc),cx(mgs,lgl),rtgl(mgs),
     :           cdx(lgl),cno(lgl),xdia(mgs,lgl,1),xmas(mgs,lc),xdia(mgs,lc,1),
     :           xdia(mgs,lc,2),vtxbar(mgs,lc,1),
     :           rfr,2,temg(mgs),rho0(mgs),cd,cs,ds,crfr)

      frcglgl = rfr(1) ! qfglacw(iqcw,iqgl,item,1,itim)
      frcglgm = rfr(2) ! qfglacw(iqcw,iqgl,item,2,itim)
      frcglgh = rfr(3) + rfr(4) ! qfglacw(iqcw,iqgl,item,3,itim)
! not allowing direct conversion to hail or frozen drops
      frcglfw = 0.0  !rfr(4) ! qfglacw(iqcw,iqgl,item,4,itim)
c
      frcglgl1 = (1.0-frcglgl)
c
      qglacwgl(mgs) =  frcglgl*qglacw(mgs)
      qglacwgm(mgs) =  frcglgm*qglacw(mgs)
      qwacglgm(mgs) =  frcglgm*qwacgl(mgs)
      qglacwgh(mgs) =  frcglgh*qglacw(mgs)
      qwacglgh(mgs) =  frcglgh*qwacgl(mgs)
      qglacwf(mgs)  =  frcglfw*qglacw(mgs)
      qwacglf(mgs)  =  frcglfw*qwacgl(mgs)
      qwacgl(mgs)   = frcglgl1*qwacgl(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
      IF ( icgxconv .eq. 0 ) THEN
c      cwacglgl(mgs) = 
c     :         Min( cglmxd(mgs), qwacglgl(mgs)*rho0(mgs)/xmas(mgs,lgl))
      cwacglgm(mgs) = 
     :         Min( cglmxd(mgs), (qwacglgm(mgs))*rho0(mgs)/xmas(mgs,lgl))
c     :         Min( cglmxd(mgs), (qglacwgm(mgs)+qwacglgm(mgs))*rho0(mgs)/xmas(mgs,lgl))
      cwacglgh(mgs) = 
     :         Min( cglmxd(mgs), (qwacglgh(mgs))*rho0(mgs)/xmas(mgs,lgl))
c     :         Min( cglmxd(mgs), (qglacwgh(mgs)+qwacglgh(mgs))*rho0(mgs)/xmas(mgs,lgl))
      cwacglf(mgs)  = 
     :         Min( cglmxd(mgs), (qwacglf(mgs))*rho0(mgs)/xmas(mgs,lgl))
c     :         Min( cglmxd(mgs), (qglacwf(mgs)+qwacglf(mgs))*rho0(mgs)/xmas(mgs,lgl))

      ELSE
!      cwacglgl(mgs) = crfr(1)*cwacgl(mgs)
      cwacglgm(mgs) = crfr(2)*cwacgl(mgs)
      cwacglgh(mgs) = (crfr(3) + crfr(4))*cwacgl(mgs)
      cwacglf(mgs)  = 0.0
      ENDIF
      ENDIF
      end if
c
      if ( temg(mgs) .gt. tfr ) then
      qwacgl(mgs) = 0.0
      qwacglgl(mgs) = qwacgl(mgs)
      qglacwgl(mgs) = qglacw(mgs)
      end if
c
      ENDIF ! qx(mgs,lgl) .gt. qxmin(lgl) .and. qx(mgs,lc) .gt. qxmin(lc)
      
      end do 
c
c
c
c
c
c
c  cloud water-graupel-m interaction
c
c
      if ( ndebug .eq. 1 ) then
      write(iunit,*) 'cloud water-graupel-m interaction'
      end if
      do mgs = 1,ngscnt
c 
      qwacgmgl(mgs)  = 0.0
      qgmacwgl(mgs)  = 0.0
      qwacgmgm(mgs)  = 0.0
      qgmacwgm(mgs)  = 0.0
      qwacgmgh(mgs)  = 0.0
      qgmacwgh(mgs)  = 0.0
      qwacgmf(mgs)   = 0.0
      qgmacwf(mgs)   = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cwacgmgl(mgs) = 0.0
      cwacgmgh(mgs) = 0.0
      cwacgmf(mgs)  = 0.0
      ENDIF
      
      IF ( qx(mgs,lgm) .gt. qxmin(lgm) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
      dtrim = Min(rimtim,rtgm(mgs))
c      cwmasr = qx(mgs,lc)*rho0(mgs)/max(cx(mgs,lc),1.e-20)
c      cwdiar = (6.*cwmasr/xdn(mgs,lc))**(1./3.)
      cwdiar = 0.5*xdia(mgs,lc,1)
      rimgm(mgs) = 1000.0
      if ( temg(mgs) .lt. 273.15 ) then
      rimgm(mgs) = rimc1*(-((0.5)*(1.e+06)*cwdiar)
     >                *((0.60)*vtxbar(mgs,lgm,1))
     >                /(temg(mgs)-273.15))**(rimc2)
      rimgm(mgs) = min( max(100.0, rimgm(mgs)), 900.0)   
      end if
c
      rimmas = qgmacw(mgs)*dtrim
c
      rhogmw(mgs) = 0.0
      rhobar =
     >  (qx(mgs,lgm)*xdn(mgs,lgm)+rimmas*rimgm(mgs)+eps)
     > /(qx(mgs,lgm)+rimmas+eps)
      rhobar = min( max(100.0, rhobar), 900.0)   
      rhogmw(mgs) = rhobar
c
c
c
      if ( temg(mgs) .le. tfr ) then
c
c      if ( itimrimw .eq. 1 ) then
c      itim = itim0
c      end if
c
c      if ( itimrimw .eq. 2 ) then
c      argtim = ((rtgm(mgs)-0.5/xftim)*xftim + yftim)
c      itim = min(max(ifix(argtim),1), mtim)
c      end if
c
c      argqcw = (qx(mgs,lc) + 5.e-4 - qxmin(lc))*xfqcw+yfqcw
c      argqxw = (qx(mgs,lgm) + 5.e-4 - qxmin(lgm))*xfqxw+yfqxw
c      argtem = -temcg(mgs)*xftem+yftem
c      iqcw = min(max(ifix(argqcw),1), mqcw)
c      iqgm = min(max(ifix(argqxw),1), mqxw)
c      item = min(max(ifix(argtem),1), mtem)
c
      call rimconv(qx(mgs,lc),qx(mgs,lgm),cx(mgs,lc),cx(mgs,lgm),rtgm(mgs),
     :           cdx(lgm),cno(lgm),xdia(mgs,lgm,1),xmas(mgs,lc),xdia(mgs,lc,1),
     :           xdia(mgs,lc,2),vtxbar(mgs,lc,1),
     :           rfr,3,temg(mgs),rho0(mgs),cd,cs,ds,crfr)

      frcgmgl = rfr(1) ! qfgmacw(iqcw,iqgm,item,1,itim)
      frcgmgm = rfr(2) ! qfgmacw(iqcw,iqgm,item,2,itim)
      frcgmgh = rfr(3) + rfr(4) ! qfgmacw(iqcw,iqgm,item,3,itim)
! Not allowing conversion directly to frozen drops/hail
      frcgmfw = 0.0   !  rfr(4) ! qfgmacw(iqcw,iqgm,item,4,itim)
c
      frcgmgm1 = (1.0-frcgmgm)
c
      qgmacwgl(mgs) =  frcgmgl*qgmacw(mgs)
      qwacgmgl(mgs) =  frcgmgl*qwacgm(mgs)
      qgmacwgm(mgs) =  frcgmgm*qgmacw(mgs)
      qgmacwgh(mgs) =  frcgmgh*qgmacw(mgs)
      qwacgmgh(mgs) =  frcgmgh*qwacgm(mgs)
      qgmacwf(mgs)  =  frcgmfw*qgmacw(mgs)
      qwacgmf(mgs)  =  frcgmfw*qwacgm(mgs)
      qwacgm(mgs)   = frcgmgm1*qwacgm(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
      IF ( icgxconv .eq. 0 ) THEN
      cwacgmgl(mgs) = 
     :         Min( cgmmxd(mgs), (qwacgmgl(mgs))*rho0(mgs)/xmas(mgs,lgm))
c     :         Min( cgmmxd(mgs), (qgmacwgl(mgs)+qwacgmgl(mgs))*rho0(mgs)/xmas(mgs,lgm))
c      cwacgmgm(mgs) = 
c     :         Min( cgmmxd(mgs), qwacgmgm(mgs)*rho0(mgs)/xmas(mgs,lgm))
      cwacgmgh(mgs) = 
     :         Min( cgmmxd(mgs), (qwacgmgh(mgs))*rho0(mgs)/xmas(mgs,lgm))
c     :         Min( cgmmxd(mgs), (qgmacwgh(mgs)+qwacgmgh(mgs))*rho0(mgs)/xmas(mgs,lgm))
      cwacgmf(mgs)  = 
     :         Min( cgmmxd(mgs), (qwacgmf(mgs))*rho0(mgs)/xmas(mgs,lgm))
c     :         Min( cgmmxd(mgs), (qgmacwf(mgs)+qwacgmf(mgs))*rho0(mgs)/xmas(mgs,lgm))

      ELSE
        cwacgmgl(mgs) = crfr(1)*cwacgm(mgs)
!        cwacgmgm(mgs) = crfr(2)*cwacgm(mgs)
        cwacgmgh(mgs) = (crfr(3) + crfr(4))*cwacgm(mgs)
        cwacgmf(mgs)  = 0.0
      ENDIF

      ENDIF

      end if
c
c
      if ( temg(mgs) .gt. tfr ) then
      qwacgm(mgs)   = 0.0
      qwacgmgm(mgs) = qwacgm(mgs)
      qgmacwgm(mgs) = qgmacw(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lgm) .gt. qxmin(lgm) .and. qx(mgs,lc) .gt. qxmin(lc) )
      
      end do
c
c
c
c  cloud water-graupel-h interaction
c
c
      if ( ndebug .eq. 1 )  then 
      write(iunit,*) 'cloud water-graupel-h interaction'
      end if
      do mgs = 1,ngscnt
c 
      qwacghgl(mgs)  = 0.0
      qghacwgl(mgs)  = 0.0
      qwacghgm(mgs)  = 0.0
      qghacwgm(mgs)  = 0.0
      qwacghgh(mgs)  = 0.0
      qghacwgh(mgs)  = 0.0
      qwacghf(mgs)   = 0.0
      qghacwf(mgs)   = 0.0
      rimgh(mgs) = 0.0
      rhoghw(mgs) = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cwacghgl(mgs) = 0.0
      cwacghgm(mgs) = 0.0
      cwacghf(mgs)  = 0.0
      ENDIF

      IF ( qx(mgs,lgh) .gt. qxmin(lgh) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN

      dtrim = Min(rimtim,rtgh(mgs))
c      cwmasr = qx(mgs,lc)*rho0(mgs)/max(cx(mgs,lc),1.e-20)
c      cwdiar = (6.*cwmasr/xdn(mgs,lc))**(1./3.)
      cwdiar = 0.5*xdia(mgs,lc,1)
      rimgh(mgs) = 1000.0
      if ( temg(mgs) .lt. 273.15 ) then
      rimgh(mgs) = rimc1*(-((0.5)*(1.e+06)*cwdiar)
     >                *((0.60)*vtxbar(mgs,lgh,1))
     >                /(temg(mgs)-273.15))**(rimc2)
      rimgh(mgs) = min( max(100.0, rimgh(mgs)), 900.0)   
      end if
c
      rimmas = qghacw(mgs)*dtrim
c
      rhobar =
     >  (qx(mgs,lgh)*xdn(mgs,lgh)+rimmas*rimgh(mgs)+eps)
     > /(qx(mgs,lgh)+rimmas+eps)
      rhobar = min( max(100.0, rhobar), 900.0)   
      rhoghw(mgs) = rhobar
c
c
c
      if ( temg(mgs) .le. tfr ) then
c
c      if ( itimrimw .eq. 1 ) then
c      itim = itim0
c      end if
c
c      if ( itimrimw .eq. 2 ) then
c      argtim = ((rtgh(mgs)-0.5/xftim)*xftim + yftim)
c      itim = min(max(ifix(argtim),1), mtim)
c      end if
c
c
c      argqcw = (qx(mgs,lc) + 5.e-4 - qxmin(lc))*xfqcw+yfqcw
c      argqxw = (qx(mgs,lgh) + 5.e-4 - qxmin(lgh))*xfqxw+yfqxw
c      argtem = -temcg(mgs)*xftem+yftem
c      iqcw = min(max(ifix(argqcw),1), mqcw)
c      iqgh = min(max(ifix(argqxw),1), mqxw)
c      item = min(max(ifix(argtem),1), mtem)
c
      call rimconv(qx(mgs,lc),qx(mgs,lgh),cx(mgs,lc),cx(mgs,lgh),rtgh(mgs),
     :           cdx(lgh),cno(lgh),xdia(mgs,lgh,1),xmas(mgs,lc),xdia(mgs,lc,1),
     :           xdia(mgs,lc,2),vtxbar(mgs,lc,1),
     :           rfr,4,temg(mgs),rho0(mgs),cd,cs,ds,crfr)

      frcghgl = rfr(1) ! qfghacw(iqcw,iqgh,item,1,itim)
      frcghgm = rfr(2) ! qfghacw(iqcw,iqgh,item,2,itim)
      frcghgh = rfr(3) + rfr(4) ! qfghacw(iqcw,iqgh,item,3,itim)
! Not allowing conversion directly to frozen drops/hail
      frcghfw = 0.0 ! rfr(4) ! qfghacw(iqcw,iqgh,item,4,itim)
      
c
      frcghgh1 = (1.0-frcghgh)
c
      qghacwgl(mgs) =  frcghgl*qghacw(mgs)
      qwacghgl(mgs) =  frcghgl*qwacgh(mgs)
      qghacwgm(mgs) =  frcghgm*qghacw(mgs)
      qwacghgm(mgs) =  frcghgm*qwacgh(mgs)
      qghacwgh(mgs) =  frcghgh*qghacw(mgs)
! only convert to frozen drops/hail if we have wet growth
!      IF ( qghshr(mgs) .lt. 0.0 ) THEN
!        qghacwf(mgs)  =  frcghfw*qghacw(mgs)
!        qwacghf(mgs)  =  frcghfw*qwacgh(mgs)
!      ELSE
        qghacwf(mgs)  = 0.0
        qwacghf(mgs)  = 0.0
        qghacwgh(mgs) = qghacwgh(mgs) + frcghfw*qghacw(mgs)
!      ENDIF

      qwacgh(mgs)   = frcghgh1*qwacgh(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
      IF ( icgxconv .eq. 0 ) THEN
       cwacghgl(mgs) = 
     :         Min( cghmxd(mgs), (qwacghgl(mgs))*rho0(mgs)/xmas(mgs,lgh))
c     :         Min( cghmxd(mgs), (qghacwgl(mgs)+qwacghgl(mgs))*rho0(mgs)/xmas(mgs,lgh))
       cwacghgm(mgs) = 
     :         Min( cghmxd(mgs), (qwacghgm(mgs))*rho0(mgs)/xmas(mgs,lgh))
c     :         Min( cghmxd(mgs), (qghacwgm(mgs)+qwacghgm(mgs))*rho0(mgs)/xmas(mgs,lgh))
c       cwacghgh(mgs) = 
c     :         Min( cghmxd(mgs), qwacghgh(mgs)*rho0(mgs)/xmas(mgs,lgh))
       cwacghf(mgs)  = 
     :         Min( cghmxd(mgs), (qwacghf(mgs))*rho0(mgs)/xmas(mgs,lgh))
c     :         Min( cghmxd(mgs), (qghacwf(mgs)+qwacghf(mgs))*rho0(mgs)/xmas(mgs,lgh))

      ELSE
      cwacghgl(mgs) = crfr(1)*cwacgh(mgs)
      cwacghgm(mgs) = crfr(2)*cwacgh(mgs)
!       IF ( qghshr(mgs) .lt. 0.0 ) THEN
!        cwacghgh(mgs) = crfr(3)*cwacgh(mgs)
!        cwacghf(mgs)  = crfr(4)*cwacgh(mgs)
!       ELSE
!        cwacghgh(mgs) = (crfr(3) + crfr(4))*cwacgh(mgs)
        cwacghf(mgs)  = 0.0
!       ENDIF
      ENDIF
      
      ENDIF

      end if
c         
c
      if ( temg(mgs) .gt. tfr ) then
      qwacgh(mgs)   = 0.0
      qwacghgh(mgs) = qwacgh(mgs)
      qghacwgh(mgs) = qghacw(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lgh) .gt. qxmin(lgh) .and. qx(mgs,lc) .gt. qxmin(lc) )
      
      end do
c
c
c  cloud water-frozen water interaction
c
c
      if ( ndebug .eq. 1 ) then 
      write(iunit,*)  'cloud water-frozen water interaction '
      end if
      do mgs = 1,ngscnt
c 
      qwacfgl(mgs)  = 0.0
      qfacwgl(mgs)  = 0.0
      qwacfgm(mgs)  = 0.0
      qfacwgm(mgs)  = 0.0
      qwacfgh(mgs)  = 0.0
      qfacwgh(mgs)  = 0.0
      qwacff(mgs)   = 0.0
      qfacwf(mgs)   = 0.0
      qwacfh(mgs)   = 0.0
      qfacwh(mgs)   = 0.0
      rimf(mgs) = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cwacfgl(mgs) = 0.0
      cwacfgm(mgs) = 0.0
      cwacfgh(mgs) = 0.0
      ENDIF

      IF ( qx(mgs,lf) .gt. qxmin(lf) .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
      
      dtrim = Min(rimtim,rtfw(mgs))
c      cwmasr = qx(mgs,lc)*rho0(mgs)/max(cx(mgs,lc),1.e-20)
c      cwdiar = (6.*cwmasr/xdn(mgs,lc))**(1./3.)
      cwdiar = 0.5*xdia(mgs,lc,1)
      rimf(mgs) = 800.0
      if ( temg(mgs) .lt. 273.15 ) then
      rimf(mgs) = rimc1*(-((0.5)*(1.e+06)*cwdiar)
     >                *((0.60)*vtxbar(mgs,lf,1))
     >                /(temg(mgs)-273.15))**(rimc2)
      rimf(mgs) = min( max(100.0, rimf(mgs)), 900.0)   
      end if
c
      rimmas = qfacw(mgs)*dtrim
c
      rhobar =
     >  (qx(mgs,lf)*xdn(mgs,lf)+rimmas*rimf(mgs)+eps)
     > /(qx(mgs,lf)+rimmas+eps)
      rhobar = min( max(100.0, rhobar), 900.0)   
c
c
c          
      if ( temg(mgs) .le. tfr ) then
c
c      if ( itimrimw .eq. 1 ) then
c      itim = itim0
c      end if
c
c      if ( itimrimw .eq. 2 ) then
c      argtim = ((rtfw(mgs)-0.5/xftim)*xftim + yftim)
c      itim = min(max(ifix(argtim),1), mtim)
c      end if
c
c
c      argqcw = (qx(mgs,lc) + 5.e-4 - qxmin(lc))*xfqcw+yfqcw
c      argqxw = (qx(mgs,lf) + 5.e-4 - qxmin(lf))*xfqxw+yfqxw
c      argtem = -temcg(mgs)*xftem+yftem
c      iqcw = min(max(ifix(argqcw),1), mqcw)
c      iqfw = min(max(ifix(argqxw),1), mqxw)
c      item = min(max(ifix(argtem),1), mtem)
c
      call rimconv(qx(mgs,lc),qx(mgs,lf),cx(mgs,lc),cx(mgs,lf),rtfw(mgs),
     :           cdx(lf),cno(lf),xdia(mgs,lf,1),xmas(mgs,lc),xdia(mgs,lc,1),
     :           xdia(mgs,lc,2),vtxbar(mgs,lc,1),
     :           rfr,5,temg(mgs),rho0(mgs),cd,cs,ds,crfr)

      frcfwgl = rfr(1) ! qffwacw(iqcw,iqfw,item,1,itim)
      frcfwgm = rfr(2) ! qffwacw(iqcw,iqfw,item,2,itim)
      frcfwgh = rfr(3) ! qffwacw(iqcw,iqfw,item,3,itim)
      frcfwfw = rfr(4) ! qffwacw(iqcw,iqfw,item,4,itim)
c
      frcfwfw1 = (1.0-frcfwfw)
c
      qfacwgl(mgs) =   frcfwgl*qfacw(mgs)
      qwacfgl(mgs) =   frcfwgl*qwacf(mgs)
      qfacwgm(mgs) =   frcfwgm*qfacw(mgs)
      qwacfgm(mgs) =   frcfwgm*qwacf(mgs)
      qfacwgh(mgs) =   frcfwgh*qfacw(mgs)
      qwacfgh(mgs) =   frcfwgh*qwacf(mgs)
      qfacwf(mgs)  =   frcfwfw*qfacw(mgs)
      qwacf(mgs)   =  frcfwfw1*qwacf(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
      IF ( icgxconv .eq. 0 ) THEN
        cwacfgl(mgs) = 
     :         Min( cfmxd(mgs), (qwacfgl(mgs))*rho0(mgs)/xmas(mgs,lf))
c     :         Min( cfmxd(mgs), (qfacwgl(mgs)+qwacfgl(mgs))*rho0(mgs)/xmas(mgs,lf))
        cwacfgm(mgs) = 
     :         Min( cfmxd(mgs), (qwacfgm(mgs))*rho0(mgs)/xmas(mgs,lf))
c     :         Min( cfmxd(mgs), (qfacwgm(mgs)+qwacfgm(mgs))*rho0(mgs)/xmas(mgs,lf))
        cwacfgh(mgs) = 
     :         Min( cfmxd(mgs), (qwacfgh(mgs))*rho0(mgs)/xmas(mgs,lf))
c     :         Min( cfmxd(mgs), (qfacwgh(mgs)+qwacfgh(mgs))*rho0(mgs)/xmas(mgs,lf))
c        cwacf(mgs)  = 
c     :         Min( cfmxd(mgs), qwacf(mgs)*rho0(mgs)/xmas(mgs,lf))

      ELSE
        cwacfgl(mgs) = crfr(1)*cwacf(mgs)
        cwacfgm(mgs) = crfr(2)*cwacf(mgs)
        cwacfgh(mgs) = crfr(3)*cwacf(mgs)
!        cwacff(mgs)  = crfr(4)*cwacf(mgs)
      ENDIF
      
      ENDIF
      end if
c       
c
      if ( temg(mgs) .gt. tfr ) then
      qwacf(mgs)   = 0.0
      qwacff(mgs) = qwacf(mgs)
      qfacwf(mgs) = qfacw(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lf) .gt. qxmin(lf) .and. qx(mgs,lc) .gt. qxmin(lc) )
      
      end do
c
c
c
c
c  hail-cloud/rain interaction
c
c
c     do mgs = 1,ngscnt
c     qhacrh(mgs) = qhacr(mgs)
c     qhacwh(mgs) = qhacw(mgs)
c     end do
c
c
c  frozen rain-rain interaction....
c
c
      do mgs = 1,ngscnt
      qfacrf(mgs) = qfacr(mgs)
      end do
c
c
c  rain-ice interaction
c
c
      do mgs = 1,ngscnt
      qracif(mgs) = qraci(mgs)
      qiacrf(mgs) = qiacr(mgs)
      cracif(mgs) = craci(mgs)
c      ciacrf(mgs) = ciacr(mgs)
      end do
c
c
c  rain-snow interaction
c
c
c
      do mgs = 1,ngscnt
c
c      qracglgl(mgs) = qracgl(mgs) ! 0.0
c      qglacrgl(mgs) = qglacr(mgs) ! 0.0

      qracsgl(mgs) = qracs(mgs)
      qsacrs(mgs) = qsacr(mgs) ! 0.0
      qracss(mgs) = 0.0
!      qracsgl(mgs) = 0.0
      qsacrgl(mgs) = 0.0
      qracsgm(mgs) = 0.0
      qsacrgm(mgs) = 0.0
      qracsgh(mgs) = 0.0
      qsacrgh(mgs) = 0.0
      qracsf(mgs)  = 0.0
      qsacrf(mgs)  = 0.0
c
      IF ( ipconc .ge. 4 ) THEN
      cracsgl(mgs) = 0.0
      cracsgm(mgs) = 0.0
      cracsgh(mgs) = 0.0
      cracsf(mgs)  = 0.0
      ENDIF
      
!      IF ( qx(mgs,ls) .gt. qxmin(ls) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
      IF ( qx(mgs,ls) .gt. 0.05e-3 .and. qx(mgs,lr) .gt. 0.05e-3 ) THEN
      if ( temg(mgs) .le. tfr ) THEN !.or. (ipconc .ge. 4 .and. qsmlr(mgs) .lt. 0. ) then
c
      IF ( ipconc .ge. 3 .or. lraincv ) THEN
c     subroutine rainconv(ipconc,qcw,qfw,ccw,cfw,
c     :           cdfw,cnof,fwdia,rwmas,rwdia,
c     :           cwdia2,vtrbar,cnor,
c     :           rfr,intyp,temg,rho0,cd,cs,ds)

      call rainconv(ipconc,qx(mgs,lr),qx(mgs,ls),cx(mgs,lr),cx(mgs,ls),
     :           cdx(ls),cno(ls),xdia(mgs,ls,1),xmas(mgs,lr),xdia(mgs,lr,1),
     :           cno(lr),xv(mgs,lr),rhovt(mgs),
     :           rfrh,rfrr,1,temg(mgs),rho0(mgs),cd,cs,ds,crfrh,crfrr)

!   qfrwasw -> rfrh    
!   qfswarw -> rfrr
      frcswrsw = 0.0
! only move into low density graupel
      frcswrgl = rfrr(2) + rfrr(1)
      frcswrgm = rfrr(3)
      frcswrgh = rfrr(4)
      frcswrfw = rfrr(5)
      frcswrsw1 = (1.0-frcswrsw)
c
c changed by (tm)
c      frcrswsw = qfrwasw(iqrw,iqsw,1)
c      frcrswgl = qfrwasw(iqrw,iqsw,2)
      frcrswsw = 0.0
      frcrswgl = rfrh(2) + rfrh(1)
      frcrswgm = rfrh(3)
      frcrswgh = rfrh(4)
      frcrswfw = rfrh(5)
      frcrswsw1 = (1.0-frcrswsw)

      ELSE
      
      argqcw = (qx(mgs,lr) + 5.e-4 - qxmin(lr))*xfqcw+yfqcw
      argqxw = (qx(mgs,ls) + 5.e-4 - qxmin(ls))*xfqxw+yfqxw
      argtem = -temcg(mgs)*xftem+yftem
      iqrw = min(max(ifix(argqcw),1), mqcw)
      iqsw = min(max(ifix(argqxw),1), mqxw)
      item = min(max(ifix(argtem),1), mtem)
c
c changed by (tm)
c      frcswrsw = qfswarw(iqrw,iqsw,1)
c      frcswrgl = qfswarw(iqrw,iqsw,2)
      frcswrsw = 0.0
      frcswrgl = qfswarw(iqrw,iqsw,2) + qfswarw(iqrw,iqsw,1)
      frcswrgm = qfswarw(iqrw,iqsw,3)
      frcswrgh = qfswarw(iqrw,iqsw,4)
      frcswrfw = qfswarw(iqrw,iqsw,5)
      frcswrsw1 = (1.0-frcswrsw)
c
c changed by (tm)
c      frcrswsw = qfrwasw(iqrw,iqsw,1)
c      frcrswgl = qfrwasw(iqrw,iqsw,2)
      frcrswsw = 0.0
      frcrswgl = qfrwasw(iqrw,iqsw,2) + qfrwasw(iqrw,iqsw,1)
      frcrswgm = qfrwasw(iqrw,iqsw,3)
      frcrswgh = qfrwasw(iqrw,iqsw,4)
      frcrswfw = qfrwasw(iqrw,iqsw,5)
      frcrswsw1 = (1.0-frcrswsw)
!      frcrswsw1 = 1.0 - (frcrswgl + frcrswgm + frcrswgh + frcrswfw)
      ENDIF
c
      qsacrgl(mgs) =  frcswrgl*qsacr(mgs)
      qracsgl(mgs) =  frcrswgl*qracs(mgs)
      qsacrgm(mgs) =  frcswrgm*qsacr(mgs)
      qracsgm(mgs) =  frcrswgm*qracs(mgs)
      qsacrgh(mgs) =  frcswrgh*qsacr(mgs)
      qracsgh(mgs) =  frcrswgh*qracs(mgs)
      qsacrf(mgs)  =  frcswrfw*qsacr(mgs)
      qracsf(mgs)  =  frcrswfw*qracs(mgs)
      qsacrs(mgs)  =  frcswrsw*qsacr(mgs)  ! note that frcswrsw=0.0
      qracss(mgs)  =  frcrswsw*qracs(mgs)  ! not used, anyway?
      qracs(mgs)   =  frcrswsw1*qracs(mgs)
c
      IF ( ipconc .ge. 4 ) THEN
!      cracsgl(mgs) = Min( csmxd(mgs), qracsgl(mgs)*rho0(mgs)/xmas(mgs,ls))
!      cracsgm(mgs) = Min( csmxd(mgs), qracsgm(mgs)*rho0(mgs)/xmas(mgs,ls))
!      cracsgh(mgs) = Min( csmxd(mgs), qracsgh(mgs)*rho0(mgs)/xmas(mgs,ls))
!      cracsf(mgs)  = Min( csmxd(mgs), qracsf(mgs)*rho0(mgs)/xmas(mgs,ls))

      cracsgl(mgs) = crfrh(2)*cracs(mgs) 
      cracsgm(mgs) = crfrh(3)*cracs(mgs) 
      cracsgh(mgs) = (crfrh(4) + crfrh(5))*cracs(mgs) 
      cracsf(mgs)  = 0.0 ! crfrh(5)*cracs(mgs) 
      ENDIF

      end if
c
      if ( temg(mgs) .gt. tfr ) then
      qracs(mgs) = 0.0
      qracsgl(mgs) = qracs(mgs)
      qracss(mgs) = qracs(mgs) ! not used, anyway?
      qsacrs(mgs) = qsacr(mgs)
      end if
c
      ENDIF ! ( qx(mgs,ls) .gt. qxmin(ls) .and. qx(mgs,lr) .gt. qxmin(lr) )

      end do
c     
c
c
c
c  rain - graupel-low density
c
c
      do mgs = 1,ngscnt
c
      qracglgl(mgs) = qracgl(mgs) ! 0.0
      qglacrgl(mgs) = qglacr(mgs) ! 0.0
      qracglgm(mgs) = 0.0
      qglacrgm(mgs) = 0.0
      qracglgh(mgs) = 0.0
      qglacrgh(mgs) = 0.0
      qracglf(mgs)  = 0.0
      qglacrf(mgs)  = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cracglgm(mgs) = 0.0
      cracglgh(mgs) = 0.0
      cracglf(mgs)  = 0.0
      ENDIF
c
      rhoglr(mgs) = xdn(mgs,lgl)
      
      IF ( iraincv .ge. 1 ) THEN
!      IF ( qx(mgs,lgl) .gt. qxmin(lgl) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
      IF ( qx(mgs,lgl) .gt. 5.e-5 .and. qx(mgs,lr) .gt. 5.e-5 ) THEN
      rhobar =
     >  (qx(mgs,lgl)*xdn(mgs,lgl)+qglacr(mgs)*dtp*xdn(mgs,lr)+eps)
     > /(qx(mgs,lgl)+qglacr(mgs)*dtp+eps)
      rhobar = min( max(100.0, rhobar), 900.0)
      rhoglr(mgs) = rhobar
c
      if ( temg(mgs) .le. tfr ) then
c
      IF ( ipconc .ge. 3 .or. lraincv ) THEN
c     subroutine rainconv(ipconc,qcw,qfw,ccw,cfw,
c     :           cdfw,cnof,fwdia,rwmas,rwdia,
c     :           cwdia2,vtrbar,cnor,
c     :           rfr,intyp,temg,rho0,cd,cs,ds)

      call rainconv(ipconc,qx(mgs,lr),qx(mgs,lgl),cx(mgs,lr),cx(mgs,lgl),
     :           cdx(lgl),cno(lgl),xdia(mgs,lgl,1),xmas(mgs,lr),xdia(mgs,lr,1),
     :           cno(lr),xv(mgs,lr),rhovt(mgs),
     :           rfrh,rfrr,2,temg(mgs),rho0(mgs),cd,cs,ds,crfrh,crfrr)

!   qfrwacgl -> rfrh    
!   qfglacr -> rfrr

      frcglrgl = rfrr(1)
      frcglrgm = rfrr(2)
      frcglrgh = rfrr(3) + rfrr(4)
! do not allow direct conversion to hail or frozen drops
      frcglrfw = 0.0  ! rfrr(4)
      frcglrgl1 = (1.0-frcglrgl)
c
      frcrglgl = rfrh(1)
      frcrglgm = rfrh(2)
      frcrglgh = rfrh(3) + rfrh(4)
      frcrglfw = 0.0 ! rfrh(4)
      frcrglgl1 = (1.0-frcrglgl)

      ELSE

      argqcw = (qx(mgs,lr) + 5.e-4 - qxmin(lr))*xfqcw+yfqcw
      argqxw = (qx(mgs,lgl) + 5.e-4 - qxmin(lgl))*xfqxw+yfqxw
      argtem = -temcg(mgs)*xftem+yftem
      iqrw = min(max(ifix(argqcw),1), mqcw)
      iqgl = min(max(ifix(argqxw),1), mqxw)
      item = min(max(ifix(argtem),1), mtem)
c
      frcglrgl = qfglarw(iqrw,iqgl,1)
      frcglrgm = qfglarw(iqrw,iqgl,2)
      frcglrgh = qfglarw(iqrw,iqgl,3)
      frcglrfw = qfglarw(iqrw,iqgl,4)
      frcglrgl1 = (1.0-frcglrgl)
c
      frcrglgl = qfrwagl(iqrw,iqgl,1)
      frcrglgm = qfrwagl(iqrw,iqgl,2)
      frcrglgh = qfrwagl(iqrw,iqgl,3)
      frcrglfw = qfrwagl(iqrw,iqgl,4)
      frcrglgl1 = (1.0-frcrglgl)
      
      ENDIF
c
      qglacrgl(mgs) =  frcglrgl*qglacr(mgs)
      qracglgl(mgs) =  frcrglgl*qracgl(mgs)
      
      qglacrgm(mgs) =  frcglrgm*qglacr(mgs)
      qracglgm(mgs) =  frcrglgm*qracgl(mgs)
      
      qglacrgh(mgs) =  frcglrgh*qglacr(mgs)
      qracglgh(mgs) =  frcrglgh*qracgl(mgs)
      
      qglacrf(mgs)  =  frcglrfw*qglacr(mgs)
      qracglf(mgs)  =  frcrglfw*qracgl(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
!c      cracglgl(mgs) = 
!c     :         Min( cglmxd(mgs), qracglgl(mgs)*rho0(mgs)/xmas(mgs,lgl))
!      cracglgm(mgs) = 
!     :         Min( cglmxd(mgs), (qracglgm(mgs))*rho0(mgs)/xmas(mgs,lgl))
!c     :         Min( cglmxd(mgs), (qracglgm(mgs)+qglacrgm(mgs))*rho0(mgs)/xmas(mgs,lgl))
!      cracglgh(mgs) = 
!     :         Min( cglmxd(mgs), (qracglgh(mgs))*rho0(mgs)/xmas(mgs,lgl))
!c     :         Min( cglmxd(mgs), (qracglgh(mgs)+qglacrgh(mgs))*rho0(mgs)/xmas(mgs,lgl))
!      cracglf(mgs)  = 
!     :         Min( cglmxd(mgs), (qracglf(mgs))*rho0(mgs)/xmas(mgs,lgl))
!c     :         Min( cglmxd(mgs), (qracglf(mgs)+qglacrf(mgs))*rho0(mgs)/xmas(mgs,lgl))


      cracglgm(mgs) = crfrh(2)*cracgl(mgs) 
      cracglgh(mgs) = (crfrh(3) + crfrh(4))*cracgl(mgs) 
      cracglf(mgs)  = 0.0 !crfrh(4)*cracgl(mgs) 
      ENDIF
      end if
c
      if ( temg(mgs) .gt. tfr ) then
      qracgl(mgs) = 0.0
      qracglgl(mgs) = qracgl(mgs)
      qglacrgl(mgs) = qglacr(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lgl) .gt. qxmin(lgl) .and. qx(mgs,lr) .gt. qxmin(lr) )
      ENDIF
      
      end do
c
c
c   rain-graupel-medium density
c
c
      do mgs = 1,ngscnt
c
      qracgmgl(mgs) = 0.0
      qgmacrgl(mgs) = 0.0
      qracgmgm(mgs) = 0.0
      qgmacrgm(mgs) = qgmacr(mgs) !  0.0
      qracgmgh(mgs) = 0.0
      qgmacrgh(mgs) = 0.0
      qracgmf(mgs)  = 0.0
      qgmacrf(mgs)  = 0.0

      IF ( ipconc .ge. 5 ) THEN
      cracgmgl(mgs) = 0.0
      cracgmgh(mgs) = 0.0
      cracgmf(mgs)  = 0.0
      ENDIF
c
      rhogmr(mgs) = xdn(mgs,lgm)
      
      IF ( iraincv .ge. 1 ) THEN
!      IF ( qx(mgs,lgm) .gt. qxmin(lgm) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
      IF ( qx(mgs,lgm) .gt. 5.e-5 .and. qx(mgs,lr) .gt. 5.e-5 ) THEN
      rhobar =
     >  (qx(mgs,lgm)*xdn(mgs,lgm)+qgmacr(mgs)*dtp*xdn(mgs,lr)+eps)
     > /(qx(mgs,lgm)+qgmacr(mgs)*dtp+eps)
      rhobar = min( max(100.0, rhobar), 900.0)
      rhogmr(mgs) = rhobar
c
      if ( temg(mgs) .le. tfr ) then
c
      IF ( ipconc .ge. 3 .or. lraincv ) THEN
      
      call rainconv(ipconc,qx(mgs,lr),qx(mgs,lgm),cx(mgs,lr),cx(mgs,lgm),
     :           cdx(lgm),cno(lgm),xdia(mgs,lgm,1),xmas(mgs,lr),xdia(mgs,lr,1),
     :           cno(lr),xv(mgs,lr),rhovt(mgs),
     :           rfrh,rfrr,3,temg(mgs),rho0(mgs),cd,cs,ds,crfrh,crfrr)

!   qfracgm -> rfrh    
!   qfgmacr -> rfrr

      frcgmrgl = rfrr(1)
      frcgmrgm = rfrr(2)
      frcgmrgh = rfrr(3) + rfrr(4)
! do not allow direct conversion to hail or frozen drops
      frcgmrfw = 0.0 ! rfrr(4)
      frcgmrgm1 = (1.0-frcgmrgm)
c
      frcrgmgl = rfrh(1)
      frcrgmgm = rfrh(2)
      frcrgmgh = rfrh(3) + rfrh(4)
      frcrgmfw = 0.0 ! rfrh(4)
      frcrgmgm1 = (1.0-frcrgmgm)
      
      ELSE
      
      argqcw = (qx(mgs,lr) + 5.e-4 - qxmin(lr))*xfqcw+yfqcw
      argqxw = (qx(mgs,lgm) + 5.e-4 - qxmin(lgm))*xfqxw+yfqxw
      argtem = -temcg(mgs)*xftem+yftem
      iqrw = min(max(ifix(argqcw),1), mqcw)
      iqgm = min(max(ifix(argqxw),1), mqxw)
      item = min(max(ifix(argtem),1), mtem)
c
      frcgmrgl = qfgmarw(iqrw,iqgm,1)
      frcgmrgm = qfgmarw(iqrw,iqgm,2)
      frcgmrgh = qfgmarw(iqrw,iqgm,3)
      frcgmrfw = qfgmarw(iqrw,iqgm,4)
      frcgmrgm1 = (1.0-frcgmrgm)
c
      frcrgmgl = qfrwagm(iqrw,iqgm,1)
      frcrgmgm = qfrwagm(iqrw,iqgm,2)
      frcrgmgh = qfrwagm(iqrw,iqgm,3)
      frcrgmfw = qfrwagm(iqrw,iqgm,4)
      frcrgmgm1 = (1.0-frcrgmgm)
      
      ENDIF
c
      qgmacrgl(mgs) =  frcgmrgl*qgmacr(mgs)
      qracgmgl(mgs) =  frcrgmgl*qracgm(mgs)
      qgmacrgm(mgs) =  frcgmrgm*qgmacr(mgs)
      qracgmgm(mgs) =  frcrgmgm*qracgm(mgs)
      qgmacrgh(mgs) =  frcgmrgh*qgmacr(mgs)
      qracgmgh(mgs) =  frcrgmgh*qracgm(mgs)
      qgmacrf(mgs)  =  frcgmrfw*qgmacr(mgs)
      qracgmf(mgs)  =  frcrgmfw*qracgm(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
!      cracgmgl(mgs) = 
!     :         Min( cgmmxd(mgs), (qracgmgl(mgs))*rho0(mgs)/xmas(mgs,lgm))
!c     :         Min( cgmmxd(mgs), (qgmacrgl(mgs)+qracgmgl(mgs))*rho0(mgs)/xmas(mgs,lgm))
!c      cracgmgm(mgs) = 
!c     :         Min( cgmmxd(mgs), qracgmgm(mgs)*rho0(mgs)/xmas(mgs,lgm))
!      cracgmgh(mgs) = 
!     :         Min( cgmmxd(mgs), (qracgmgh(mgs))*rho0(mgs)/xmas(mgs,lgm))
!c     :         Min( cgmmxd(mgs), (qgmacrgh(mgs)+qracgmgh(mgs))*rho0(mgs)/xmas(mgs,lgm))
!      cracgmf(mgs)  = 
!     :         Min( cgmmxd(mgs), (qracgmf(mgs))*rho0(mgs)/xmas(mgs,lgm))
!c     :         Min( cgmmxd(mgs), (qracgmf(mgs)+qgmacrf(mgs))*rho0(mgs)/xmas(mgs,lgm))

      cracgmgl(mgs) = crfrh(1)*cracgm(mgs) 
      cracgmgh(mgs) = (crfrh(3) + crfrh(4))*cracgm(mgs) 
      cracgmf(mgs)  = 0.0 ! crfrh(4)*cracgm(mgs) 

      ENDIF
      end if
c
      if ( temg(mgs) .gt. tfr ) then
      qracgm(mgs) = 0.0
      qracgmgm(mgs) = qracgm(mgs)
      qgmacrgm(mgs) = qgmacr(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lgm) .gt. qxmin(lgm) .and. qx(mgs,lr) .gt. qxmin(lr) )
      ENDIF
      
      end do
c
c
c
c  rain graupel-high density
c
c
      do mgs = 1,ngscnt
c
      qracghgl(mgs) = 0.0
      qghacrgl(mgs) = 0.0
      qracghgm(mgs) = 0.0
      qghacrgm(mgs) = 0.0
      qracghgh(mgs) = 0.0
      qghacrgh(mgs) = qghacr(mgs) ! 0.0
      qracghf(mgs)  = 0.0
      qghacrf(mgs)  = 0.0
c
      IF ( ipconc .ge. 5 ) THEN
      cracghgl(mgs) = 0.0
      cracghgm(mgs) = 0.0
      cracghf(mgs)  = 0.0
      ENDIF

      rhoghr(mgs) = xdn(mgs,lgh)
      
      IF ( iraincv .ge. 1 ) THEN
!      IF ( qx(mgs,lgh) .gt. qxmin(lgh) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
      IF ( qx(mgs,lgh) .gt. 5.e-5 .and. qx(mgs,lr) .gt. 5.e-5 ) THEN
      
      rhobar =
     >  (qx(mgs,lgh)*xdn(mgs,lgh)+qghacr(mgs)*dtp*xdn(mgs,lr)+eps)
     > /(qx(mgs,lgh)+qghacr(mgs)*dtp+eps)
      rhobar = min( max(100.0, rhobar), 900.0)
      rhoghr(mgs) = rhobar
c
      if ( temg(mgs) .le. tfr ) then
c
      IF ( ipconc .ge. 3 .or. lraincv ) THEN
      
      call rainconv(ipconc,qx(mgs,lr),qx(mgs,lgh),cx(mgs,lr),cx(mgs,lgh),
     :           cdx(lgh),cno(lgh),xdia(mgs,lgh,1),xmas(mgs,lr),xdia(mgs,lr,1),
     :           cno(lr),xv(mgs,lr),rhovt(mgs),
     :           rfrh,rfrr,4,temg(mgs),rho0(mgs),cd,cs,ds,crfrh,crfrr)

!   qfrwagm -> rfrh    
!   qfgmarw -> rfrr
      
      frcghrgl = rfrr(1)
      frcghrgm = rfrr(2)
      frcghrgh = rfrr(3)
! only convert to frozen drops/hail if we have wet growth
!      IF ( qghshr(mgs) .lt. 0.0 ) THEN
!        frcghrfw = rfrr(4)
!      ELSE
        frcghrfw = 0.0 ! rfrr(4)
        frcghrgh = rfrr(3) + rfrr(4)
!      ENDIF
      
      frcghrgh1 = (1.0-frcghrgh)
c
      frcrghgl = rfrh(1)
      frcrghgm = rfrh(2)
      frcrghgh = rfrh(3)
! only convert to frozen drops/hail if we have wet growth
!      IF ( qghshr(mgs) .lt. 0.0 ) THEN
!        frcrghfw = rfrh(4)
!      ELSE
        frcrghfw = 0.0
        frcrghgh = rfrh(3) + rfrh(4)
!      ENDIF
      frcrghgh1 = (1.0-frcrghgh)

      ELSE
      
      argqcw = (qx(mgs,lr) + 5.e-4 - qxmin(lr))*xfqcw+yfqcw
      argqxw = (qx(mgs,lgh) + 5.e-4 - qxmin(lgh))*xfqxw+yfqxw
      argtem = -temcg(mgs)*xftem+yftem
      iqrw = min(max(ifix(argqcw),1), mqcw)
      iqgh = min(max(ifix(argqxw),1), mqxw)
      item = min(max(ifix(argtem),1), mtem)
c
      frcghrgl = qfgharw(iqrw,iqgh,1)
      frcghrgm = qfgharw(iqrw,iqgh,2)
      frcghrgh = qfgharw(iqrw,iqgh,3)
      frcghrfw = qfgharw(iqrw,iqgh,4)
      frcghrgh1 = (1.0-frcghrgh)
c
      frcrghgl = qfrwagh(iqrw,iqgh,1)
      frcrghgm = qfrwagh(iqrw,iqgh,2)
      frcrghgh = qfrwagh(iqrw,iqgh,3)
      frcrghfw = qfrwagh(iqrw,iqgh,4)
      frcrghgh1 = (1.0-frcrghgh)
      
      ENDIF
c
      qghacrgl(mgs) =  frcghrgl*qghacr(mgs)
      qracghgl(mgs) =  frcrghgl*qracgh(mgs)
      qghacrgm(mgs) =  frcghrgm*qghacr(mgs)
      qracghgm(mgs) =  frcrghgm*qracgh(mgs)
      qghacrgh(mgs) =  frcghrgh*qghacr(mgs)
      qracghgh(mgs) =  frcrghgh*qracgh(mgs)
      qghacrf(mgs)  =  frcghrfw*qghacr(mgs)
      qracghf(mgs)  =  frcrghfw*qracgh(mgs)
c
      IF ( ipconc .ge. 5 ) THEN
!      cracghgl(mgs) = 
!     :         Min( cghmxd(mgs), (qracghgl(mgs))*rho0(mgs)/xmas(mgs,lgh))
!c     :         Min( cghmxd(mgs), (qghacrgl(mgs)+qracghgl(mgs))*rho0(mgs)/xmas(mgs,lgh))
!      cracghgm(mgs) = 
!     :         Min( cghmxd(mgs), (qracghgm(mgs))*rho0(mgs)/xmas(mgs,lgh))
!c     :         Min( cghmxd(mgs), (qghacrgm(mgs)+qracghgm(mgs))*rho0(mgs)/xmas(mgs,lgh))
!c      cracghgh(mgs) = 
!c     :         Min( cghmxd(mgs), qracghgh(mgs)*rho0(mgs)/xmas(mgs,lgh))
!      cracghf(mgs)  = 
!     :         Min( cghmxd(mgs), (qracghf(mgs))*rho0(mgs)/xmas(mgs,lgh))
!c     :         Min( cghmxd(mgs), (qghacrf(mgs)+qracghf(mgs))*rho0(mgs)/xmas(mgs,lgh))

      cracghgl(mgs) = crfrh(1)*cracgh(mgs) 
      cracghgm(mgs) = crfrh(2)*cracgh(mgs) 
!      IF ( qghshr(mgs) .lt. 0.0 ) THEN
!      cracghf(mgs)  = crfrh(4)*cracgh(mgs) 
!      ELSE
      cracghf(mgs) = 0.0
!      ENDIF

      ENDIF
      end if
c
      if ( temg(mgs) .gt. tfr ) then
      qracgh(mgs) = 0.0
      qracghgh(mgs) = qracgh(mgs)
      qghacrgh(mgs) = qghacr(mgs)
      end if
c
      ENDIF ! ( qx(mgs,lgh) .gt. qxmin(lgh) .and. qx(mgs,lr) .gt. qxmin(lr) )
      ENDIF
      
      end do
c
c
c
      do mgs = 1,ngscnt
c      if ( .not.(rhoglw(mgs).gt.800.0 .or. rhoglr(mgs).gt.800.0 )) then
      qhcngl(mgs) = 0.0
      chcngl(mgs) = 0.0
c      end if
!      if ( .not.(rhogmw(mgs).gt.800.0 .or. rhogmr(mgs).gt.800.0
!     :            .or. cwacgmf(mgs) .gt. 1. )) then
      qhcngm(mgs) = 0.0
      chcngm(mgs) = 0.0
!      end if
      if (  qghshr(mgs) .ge. 0.0 .or. 
     :     .not.(rhoghw(mgs).gt.800.0 .or. rhoghr(mgs).gt.800.0
     :            .or. cwacghf(mgs) .gt. 1. )) then
      qhcngh(mgs) = 0.0
      chcngh(mgs) = 0.0
      end if
c      IF ( rimf(mgs) .le. 850.0 ) THEN
c       qhcnf(mgs) = 0.0
c       chcnf(mgs) = 0.0
c      ENDIF
      end do
c
c  end of basic sources and sinks
c
c
c
c
c
c
c  prognose (diagnose at this stage) temperature of hydrometeors...
c
c
c
c
      if ( iptemp .ge. 1 ) then
c
c
      qweps = 1.e-20
c
      do mgs = 1,ngscnt
      tprw(mgs) = max(temg(mgs),273.15)
      end do
c
      do mgs = 1,ngscnt
      tpcw(mgs) = max(temg(mgs),273.15)
      end do
c
      do mgs = 1,ngscnt
      tpip(mgs) = min(temg(mgs),273.15)
      end do
c
      do mgs = 1,ngscnt
      tpci(mgs) = min(temg(mgs),273.15)
      end do
c
      do mgs = 1,ngscnt
      tpir(mgs) = min(temg(mgs),273.15)
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tpsw(mgs) .gt. 273.15 ) tpsw(mgs) = 273.15
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,ls) .lt. 1.e-7 ) then
      tpsw(mgs) = min(temg(mgs),273.15)
      end if
      if ( qx(mgs,ls) .gt. 1.e-7 ) then
      dqldt = qsacr(mgs)+qsacw(mgs)
      dqidt = qsaci(mgs)
     >      + qsacir(mgs)+qsacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tpsw(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,ls))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,ls)*rho0(mgs))
      a3 = cx(mgs,ls)*2.*pi*xdia(mgs,ls,1)*swvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tpsw(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tpsw(mgs) = tpsw(mgs) + dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tpsw(mgs) .lt. temg(mgs) ) tpsw(mgs) = temg(mgs)
      end if
      end if
      if ( tpsw(mgs) .gt. 273.15 ) tpsw(mgs) = 273.15
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tpgl(mgs) .gt. 273.15 ) tpgl(mgs) = 273.15
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lgl) .lt. 1.e-7 ) then
      tpgl(mgs) = min(temg(mgs),273.15)
      end if
      if ( qx(mgs,lgl) .gt. 1.e-7 ) then
      dqldt = qglacr(mgs)+qglacw(mgs)
      dqidt = qglaci(mgs)+qglacs(mgs)
     >      + qglacir(mgs)+qglacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tpgl(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lgl))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lgl)*rho0(mgs))
      a3 = cx(mgs,lgl)*2.*pi*xdia(mgs,lgl,1)*glvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tpgl(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tpgl(mgs) = tpgl(mgs)+dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tpgl(mgs) .lt. temg(mgs) ) tpgl(mgs) = temg(mgs)
      end if
      end if
      if ( tpgl(mgs) .gt. 273.15 ) tpgl(mgs) = 273.15
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tpgm(mgs) .gt. 273.15 ) tpgm(mgs) = 273.15
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lgm) .lt. 1.e-7 ) then
      tpgm(mgs) = min(temg(mgs),273.15)
      end if
      if ( qx(mgs,lgm) .gt. 1.e-7 ) then
      dqldt = qgmacr(mgs)+qgmacw(mgs)
      dqidt = qgmaci(mgs)+qgmacs(mgs)
     >      + qgmacir(mgs)+qgmacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tpgm(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lgm))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lgm)*rho0(mgs))
      a3 = cx(mgs,lgm)*2.*pi*xdia(mgs,lgm,1)*gmvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tpgm(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tpgm(mgs) = tpgm(mgs) +dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tpgm(mgs) .lt. temg(mgs) ) tpgm(mgs) = temg(mgs)
      end if
      end if
      if ( tpgm(mgs) .gt. 273.15 ) tpgm(mgs) = 273.15
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tpgh(mgs) .gt. 273.15 ) tpgh(mgs) = 273.15
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lgh) .lt. 1.e-7 ) then
      tpgh(mgs) = min(temg(mgs), 273.15)
      end if
      if ( qx(mgs,lgh) .gt. 1.e-7 ) then
      dqldt = qghacr(mgs)+qghacw(mgs)
      dqidt = qghaci(mgs)+qghacs(mgs)
     >      + qghacir(mgs)+qghacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tpgh(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lgh))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lgh)*rho0(mgs))
      a3 = cx(mgs,lgh)*2.*pi*xdia(mgs,lgh,1)*ghvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tpgh(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tpgh(mgs) = tpgh(mgs)+dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tpgh(mgs) .lt. temg(mgs) ) tpgh(mgs) = temg(mgs)
      end if
      end if
      if ( tpgh(mgs) .gt. 273.15 ) tpgh(mgs) = 273.15
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tpfw(mgs) .gt. 273.15 ) tpfw(mgs) = 273.15
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lf) .lt. 1.e-7 ) then
      tpfw(mgs) = min(temg(mgs), 273.15)
      end if
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
      dqldt = qfacr(mgs)+qfacw(mgs)
      dqidt = qfaci(mgs)+qfacs(mgs)
     >      + qfacir(mgs)+qfacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tpfw(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lf))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lf)*rho0(mgs))
      a3 = cx(mgs,lf)*2.*pi*xdia(mgs,lf,1)*fwvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tpfw(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tpfw(mgs) = tpfw(mgs)+dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tpfw(mgs) .lt. temg(mgs) ) tpfw(mgs) = temg(mgs)
      end if
      end if
      if ( tpfw(mgs) .gt. 273.15 ) tpfw(mgs) = 273.15
      end do
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tphw(mgs) .gt. 273.15 ) tphw(mgs) = 273.15 
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lh) .lt. 1.e-7 ) then
      tphw(mgs) = min(temg(mgs), 273.15)
      end if
      if ( qx(mgs,lh) .gt. 1.e-7 ) then
      dqldt = qhacr(mgs)+qhacw(mgs)
      dqidt = qhaci(mgs)+qhacs(mgs)
     >      + qhacir(mgs)+qhacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tphw(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lh))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lh)*rho0(mgs))
      a3 = cx(mgs,lh)*2.*pi*xdia(mgs,lh,1)*hwvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tphw(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tphw(mgs) = tphw(mgs) +dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tphw(mgs) .lt. temg(mgs) ) tphw(mgs) = temg(mgs)
      end if
      end if
      if ( tphw(mgs) .gt. 273.15 ) tphw(mgs) = 273.15 
      end do
c
c
      gf1a = 1.
      gf2a = 1.
      gf4a = 6.
      do mgs = 1,ngscnt
      if ( tphl(mgs) .gt. 273.15 ) tphl(mgs) = 273.15 
      if ( temg(mgs) .lt. 273.15 ) then
      if ( qx(mgs,lhl) .le. qxmin(lh) ) then
      tphl(mgs) = min(temg(mgs), 273.15)
      end if
      if ( qx(mgs,lhl) .gt. qxmin(lh) ) then
      dqldt = qhlacr(mgs)+qhlacw(mgs)
      dqidt = qhlaci(mgs)+qhlacs(mgs)
     >      + qhlacir(mgs)+qhlacip(mgs)
      dqdt = dqldt + dqidt
      a1= -(tphl(mgs)-273.15)*gf1a*dqdt/(gf4a*qx(mgs,lhl))
      a2 = 1./(fci(mgs)*gf4a*qx(mgs,lhl)*rho0(mgs))
      a3 = cx(mgs,lhl)*2.*pi*xdia(mgs,lhl,1)*hlvent(mgs)*gf2a
      a4 = (ftka(mgs)*temcg(mgs)-ftka(mgs)*(tphl(mgs)-273.15)
     >     -fels(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv)))
      a5 = dqldt*(felf(mgs)+fcw(mgs)*(temcg(mgs)))*gf1a
      a6 = dqidt*fci(mgs)*temcg(mgs)*gf1a
      tphl(mgs) = tphl(mgs) +dtp*(a1+a2*(a3*a4+a5+a6))
      if ( tphl(mgs) .lt. temg(mgs) ) tphl(mgs) = temg(mgs)
      end if
      end if
      if ( tphl(mgs) .gt. 273.15 ) tphl(mgs) = 273.15 
      end do
c
c
      end if 
c
c
c
c
c
c
c
c
c
c
c
      if ( nstep  .le. -1000 .and. 
     >     ndebug .ge. 1) then
      do mgs = 1,ngscnt
      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++'
      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++'
      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++'
c
      print*,'   '
      print*,'nstep,mgs',nstep,mgs,igs(mgs),jy,kgs(mgs)
      print*,'   '
      print*,'sxsaci(mgs),csaci(mgs)',sxsaci(mgs)
      print*,'sxsacip(mgs),csacip(mgs)',sxsacip(mgs)
      print*,'sxsacir(mgs),csacir(mgs)',sxsacir(mgs)
c      print*,'sxsacw(mgs),csacw(mgs)',sxsacw(mgs)
c      print*,'sxsacr(mgs),csacr(mgs)',sxsacr(mgs)
c
      print*,'sxglaci(mgs),cglaci(mgs)',sxglaci(mgs)
      print*,'sxglacip(mgs),cglacip(mgs)',sxglacip(mgs)
      print*,'sxglacir(mgs),cglacir(mgs)',sxglacir(mgs)
c      print*,'sxglacw(mgs),cglacw(mgs)',sxglacw(mgs)
c      print*,'sxglacr(mgs),cglacr(mgs)',sxglacr(mgs)
c
      print*,'sxglaci(mgs),cglaci(mgs)',sxglaci(mgs)
      print*,'sxglacip(mgs),cglacip(mgs)',sxglacip(mgs)
      print*,'sxglacir(mgs),cglacir(mgs)',sxglacir(mgs)
c      print*,'sxglacw(mgs),cglacw(mgs)',sxglacw(mgs)
c      print*,'sxglacr(mgs),cglacr(mgs)',sxglacr(mgs)
c
      print*,'sxgmaci(mgs),cgmaci(mgs)',sxgmaci(mgs)
      print*,'sxgmacip(mgs),cgmacip(mgs)',sxgmacip(mgs)
      print*,'sxgmacir(mgs),cgmacir(mgs)',sxgmacir(mgs)
c      print*,'sxgmacw(mgs),cgmacw(mgs)',sxgmacw(mgs)
c      print*,'sxgmacr(mgs),cgmacr(mgs)',sxgmacr(mgs)
c
      print*,'sxghaci(mgs),cghaci(mgs)',sxghaci(mgs)
      print*,'sxghacip(mgs),cghacip(mgs)',sxghacip(mgs)
      print*,'sxghacir(mgs),cghacir(mgs)',sxghacir(mgs)
c      print*,'sxghacw(mgs),cghacw(mgs)',sxghacw(mgs)
c      print*,'sxghacr(mgs),cghacr(mgs)',sxghacr(mgs)
c
      print*,'sxfaci(mgs),cfaci(mgs)',sxfaci(mgs)
      print*,'sxfacip(mgs),cfacip(mgs)',sxfacip(mgs)
      print*,'sxfacir(mgs),cfacir(mgs)',sxfacir(mgs)
c      print*,'sxfacw(mgs),cfacw(mgs)',sxfacw(mgs)
c      print*,'sxfacr(mgs),cfacr(mgs)',sxfacr(mgs)
c
      print*,'sxhaci(mgs),chaci(mgs)',sxhaci(mgs)
      print*,'sxhacip(mgs),chacip(mgs)',sxhacip(mgs)
      print*,'sxhacir(mgs),chacir(mgs)',sxhacir(mgs)
c      print*,'sxhacw(mgs),chacw(mgs)',sxhacw(mgs)
c      print*,'sxhacr(mgs),chacr(mgs)',sxhacr(mgs)
c
      end do
c
      end if
c
c
c
c
c
c
c
c
c
      IF ( ipconc .ge. 1 ) THEN
c
c
c  concentration production terms
c
c  YYY
c
c
       DO mgs = 1,ngscnt
       pccwi(mgs) = 0.0
       pccwd(mgs) = 0.0
       pccii(mgs) = 0.0
       pccid(mgs) = 0.0
       pcipi(mgs) = 0.0
       pcipd(mgs) = 0.0
       pciri(mgs) = 0.0
       pcird(mgs) = 0.0
       pcrwi(mgs) = 0.0
       pcrwd(mgs) = 0.0
       pcswi(mgs) = 0.0
       pcswd(mgs) = 0.0
       pcgli(mgs) = 0.0
       pcgld(mgs) = 0.0
       pcgmi(mgs) = 0.0
       pcgmd(mgs) = 0.0
       pcghi(mgs) = 0.0
       pcghd(mgs) = 0.0
       pcfwi(mgs) = 0.0
       pcfwd(mgs) = 0.0
       pchwi(mgs) = 0.0
       pchwd(mgs) = 0.0
       ENDDO
c
c  Cloud ice
c
c      IF ( ipconc .ge. 1 ) THEN

      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      pccii(mgs) =
     >   il5(mgs)*cicint(mgs)
c     >  +il5(mgs)*cidpv(mgs)
     >  +il5(mgs)*(cwacii(mgs)+cwfrzc(mgs)+cwctfzc(mgs)
     >  +cicichr(mgs))
     >  +chmul1(mgs)+cfmul1(mgs)+cghmul1(mgs)
     >  +cgmmul1(mgs)+cglmul1(mgs)+csmul1(mgs)
      pccid(mgs) =
     >   il5(mgs)*(-cscni(mgs) - cwaci(mgs)
     >  -craci(mgs) -csaci(mgs)
     >  -cglaci(mgs) -cgmaci(mgs) -cghaci(mgs)
     >  -cfaci(mgs)  -chaci(mgs) - chlaci(mgs))
     >  -il5(mgs)*cipaci(mgs)
     >  -il5(mgs)*ciraci(mgs)
     >  
     >  +il5(mgs)*cisbv(mgs)
     >  -(1.-il5(mgs))*cimlr(mgs)
      end do
      end if
c
c  Plate ice
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      pcipi(mgs) =
     >   il5(mgs)*cipint(mgs)
c     >  +il5(mgs)*cipdpv(mgs)
     >  +il5(mgs)*(cwacipip(mgs)+cwfrzp(mgs)+cwctfzp(mgs)  !  cipacwi(mgs)+
     >  +cipiphr(mgs))
c    >  +chmul1(mgs)+cfmul1(mgs)+cghmul1(mgs)
c    >  +cgmmul1(mgs)+cglmul1(mgs)+csmul1(mgs)
      pcipd(mgs) =
     >   il5(mgs)*(-cscnip(mgs) -cwacip(mgs)
     >  -cracip(mgs) -csacip(mgs) -ciacip(mgs)
     >  -cglacip(mgs) -cgmacip(mgs) -cghacip(mgs)
     >  -cfacip(mgs)  -chacip(mgs) -chlacip(mgs))
     >  -il5(mgs)*ciracip(mgs)
     >  +il5(mgs)*cipsbv(mgs)
     >  -(1.-il5(mgs))*cipmlr(mgs)
      end do
      end if
c
c  Rimed ice
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      pciri(mgs) =
c     >  +il5(mgs)*cirdpv(mgs)
     >   il5(mgs)*cwaciir(mgs) ! (ciacwir(mgs)+cwaciir(mgs))
     >  +il5(mgs)*cwacirir(mgs) ! (ciracwir(mgs)+cwacirir(mgs))
     >  +il5(mgs)*(cipacwir(mgs)+cwacipir(mgs))
      pcird(mgs) =
     >   il5(mgs)*(-cscnir(mgs) -cwacir(mgs)
     >  -cracir(mgs) -csacir(mgs)
     >  -cglacir(mgs) -cgmacir(mgs) -cghacir(mgs)
     >  -cfacir(mgs)  -chacir(mgs))
     >  +il5(mgs)*cirsbv(mgs)
     >  -il5(mgs)*ciacir(mgs)
     >  -il5(mgs)*cipacir(mgs)
     >  -il5(mgs)*chlacir(mgs)
     >  +(1-il5(mgs))*(cirmlw(mgs) + cirmlr(mgs))
      end do
      end if
      
c      ENDIF ! ( ipconc .ge. 1 )
c
c  Cloud water
c
      IF ( ipconc .ge. 2 ) THEN
      
      do mgs = 1,ngscnt
      pccwi(mgs) =  (0.0) + (1-il5(mgs))*(-cirmlw(mgs))
      pccwd(mgs) = 
     :  - cautn(mgs) +
     >  il5(mgs)*(-ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)
     >  -cwfrzc(mgs)-cwctfzc(mgs)
c     >  -cicichr(mgs)-cipiphr(mgs) 
     :   )
c     >  -il5(mgs)*(cirirhr(mgs))
     >  -il5(mgs)*(ciracw(mgs))
     >  -il5(mgs)*(cipacw(mgs))
     >  -cracw(mgs) -csacw(mgs)  -chacw(mgs)
     >  -cglacw(mgs) -cgmacw(mgs) -cghacw(mgs) -cfacw(mgs)
     >  -chlacw(mgs)
      end do
      
      ENDIF

c
c  Rain
c
      IF ( ipconc .ge. 3 ) THEN

      do mgs = 1,ngscnt
      pcrwi(mgs) = 0.0
ccc     >   cracw(mgs) + 
     :  + crcnw(mgs)
     >  +(1-il5(mgs))*
     >   (-chmlrr(mgs)-csmlr(mgs) -cglmlrr(mgs)-cgmmlrr(mgs)
     >    -cghmlrr(mgs)-cfmlrr(mgs)
     <    -cirmlr(mgs)-cipmlr(mgs)-chlmlrr(mgs))/rz
     <    -crshr(mgs)
      pcrwd(mgs) = 0.0 
     >   + il5(mgs)*(-ciacr(mgs) - crfrz(mgs) ) ! - cipacr(mgs))
cc     >  -csacr(mgs)
     :  - chacr(mgs) - chlacr(mgs)
     >  - cglacr(mgs) - cgmacr(mgs) - cghacr(mgs)
     >  - cfacr(mgs)
     >  +crcev(mgs)
     :  - cracr(mgs)
cc     >  -il5(mgs)*ciracr(mgs)
      end do
      
      ENDIF
c
c  Snow
c
      IF ( ipconc .ge. 4 ) THEN ! 

      do mgs = 1,ngscnt
      pcswi(mgs) =
     >   il5(mgs)*(cscni(mgs)*0.5 + cscnvi(mgs) )  ! +csdpv(mgs))
c     >   +csacrs(mgs)
c     >  +il5(mgs)*csacir(mgs)
c     >  +il5(mgs)*csacip(mgs)
     >  +il5(mgs)*cscnir(mgs)*0.5
     >  +il5(mgs)*cscnip(mgs)*0.5
     >  +il5(mgs)*(ciacir(mgs) + ciraci(mgs))*0.5
     >  +il5(mgs)*(ciacip(mgs) + cipaci(mgs))*0.5
     >  +il5(mgs)*(ciracip(mgs) + cipacir(mgs))*0.5
      pcswd(mgs) =
c     >  -cwacs(mgs)-
     :  cracs(mgs) - cfacs(mgs) - chlacs(mgs)
     >  -chacs(mgs) - cglacs(mgs) - cgmacs(mgs) - cghacs(mgs)
     >  +(1-il5(mgs))*csmlr(mgs) + csshr(mgs) ! + csshrp(mgs)
     >  +il5(mgs)*(cssbv(mgs))
     :  - csacs(mgs)
     :  - cwacsgl(mgs)
     :  - cwacsgm(mgs)
     :  - cwacsgh(mgs)
     :  - cwacsf(mgs)
c     >  -csmul1(mgs)
c         IF ( ny .eq. 2 .and.
c     : (pcswi(mgs) .gt. 1.0 .or. pcswd(mgs) .ne. 0.0  )) THEN
c          print*, 'CSW: pi/d', igs(mgs),kgs(mgs),
c     :     pcswi(mgs)*1.e-3,pcswd(mgs)*1.e-3,
c     :     il5(mgs)*(cscni(mgs)*0.5),
c     >    il5(mgs)*cscnir(mgs)*0.5,
c     >    il5(mgs)*cscnip(mgs)*0.5,
c     :  qx(mgs,lir)*1.e3,cx(mgs,lir),xmas(mgs,lir)
cc     >    il5(mgs)*(ciacir(mgs) + ciraci(mgs)),
cc     >    il5(mgs)*(ciacip(mgs) + cipacip(mgs)),
cc     >    il5(mgs)*(ciracip(mgs) + cipacir(mgs))
c       ENDIF
      end do
      
      ENDIF
c
c  Graupel-L
c
      IF ( ipconc .ge. 5 ) THEN

      do mgs = 1,ngscnt
      pcgli(mgs) =
c     >  +il5(mgs)*(cgldpv(mgs))
     >  +il5(mgs)*(cwacirgl(mgs))
     >  +il5(mgs)*(cwacipgl(mgs))
     >  +il5(mgs)*(ciacwgl(mgs)+cwacigl(mgs))
c     >  +csacwgl(mgs)
     :  +cwacsgl(mgs)
c     >  +csacrgl(mgs)
     :  +cracsgl(mgs)
c     >  +cglacwgl(mgs)
c     >  +cglacrgl(mgs)
c     >  +cgmacwgl(mgs)
     :  +cwacgmgl(mgs)
c     >  +cgmacrgl(mgs)
     :  +cracgmgl(mgs)
c     >  +cghacwgl(mgs)
     :  +cwacghgl(mgs)
c     >  +cghacrgl(mgs)
     :  +cracghgl(mgs)
     >  +cwacfgl(mgs)
c     >  +cglacs(mgs)+cglaci(mgs)
c     >  +cglacir(mgs)
c     >  +cglacip(mgs)*il5(mgs)
      pcgld(mgs) =
!     >   cglshr(mgs) ! + cglshrp(mgs)
     >  +(1-il5(mgs))*cglmlr(mgs)
     >  +il5(mgs)*cglsbv(mgs)
     >  -chcngl(mgs)
     >  -cracglgh(mgs)-cwacglgh(mgs)
     >  -cracglgm(mgs)-cwacglgm(mgs)
     >  -cracglf(mgs)-cwacglf(mgs)
c     >  -cglmul1(mgs)
      end do
c
c  Graupel-M
c
      do mgs = 1,ngscnt
      pcgmi(mgs) =
c     >  +il5(mgs)*(cgmdpv(mgs))
c     >  +il5(mgs)*(ciracwgm(mgs)+cwacirgm(mgs))
c     >  +il5(mgs)*(cipacwgm(mgs)+cwacipgm(mgs))
c     >  +il5(mgs)*(ciacwgm(mgs)+cwacigm(mgs))
c     >  +csacwgm(mgs)
     :  +cwacsgm(mgs)
c     >  +csacrgm(mgs)
     :  +cracsgm(mgs)
c     >  +cglacwgm(mgs)
     :  +cwacglgm(mgs)
c     >  +cglacrgm(mgs)
     :  +cracglgm(mgs)
c     >  +cgmacwgm(mgs)
c     >  +cgmacrgm(mgs)
c     >  +cghacwgm(mgs)
     :  +cwacghgm(mgs)
c     >  +cghacrgm(mgs)
     :  +cracghgm(mgs)
     >  +cwacfgm(mgs)
c     >  +cgmacs(mgs)+cgmaci(mgs)
c     >  +cgmacir(mgs)
c     >  +cgmacip(mgs)*il5(mgs)
     >  -(1-il5(mgs))*cglmlr(mgs)
      pcgmd(mgs) =
!     >   cgmshr(mgs) ! + cgmshrp(mgs)
     >  +(1-il5(mgs))*cgmmlr(mgs)
     >  +il5(mgs)*cgmsbv(mgs)
     >  -chcngm(mgs)
     >  -cracgmgh(mgs)-cwacgmgh(mgs)
     >  -cracgmgl(mgs)-cwacgmgl(mgs)
     >  -cracgmf(mgs)-cwacgmf(mgs)
c     >  -cgmmul1(mgs)
      end do
c
c  Graupel-H
c
      do mgs = 1,ngscnt
      pcghi(mgs) =
c     >  +il5(mgs)*(cghdpv(mgs))
c     >  +il5(mgs)*(ciracwgh(mgs)+cwacirgh(mgs))
c     >  +il5(mgs)*(cipacwgh(mgs)+cwacipgh(mgs))
c     >  +il5(mgs)*(ciacwgh(mgs)+cwacigh(mgs))
c     >  +csacwgh(mgs)
     :  +cwacsgh(mgs)
c     >  +csacrgh(mgs)
     :  +cracsgh(mgs)
c     >  +cglacwgh(mgs)
     :  +cwacglgh(mgs)
c     >  +cglacrgh(mgs)
     :  +cracglgh(mgs)
c     >  +cgmacwgh(mgs)
     :  +cwacgmgh(mgs)
c     >  +cgmacrgh(mgs)
     :  +cracgmgh(mgs)
c     >  +cghacwgh(mgs)
c     >  +cghacrgh(mgs)
     >  +cwacfgh(mgs)
c     >  +cghacs(mgs)+cghaci(mgs)
c     >  +cghacir(mgs)
c     >  +cghacip(mgs)*il5(mgs)
     >  -(1-il5(mgs))*cgmmlr(mgs)
      pcghd(mgs) =
c     >   cghshr(mgs) + 
c     :   cghshrp(mgs) +
     >  (1-il5(mgs))*cghmlr(mgs)
     >  +il5(mgs)*cghsbv(mgs)-chcngh(mgs)
     >  -cracghgl(mgs)-cwacghgl(mgs)
     >  -cracghgm(mgs)-cwacghgm(mgs)
     >  -cracghf(mgs)-cwacghf(mgs)
c     >  -cghmul1(mgs)
      end do
c
c  frozen drops
c
      do mgs = 1,ngscnt
      pcfwi(mgs) =
     :  +crfrzf(mgs)
c     >  +il5(mgs)*(cracif(mgs)+ciacrf(mgs)) ! +cfdpv(mgs))
c     >  +il5(mgs)*(cracir(mgs)+ciracr(mgs))
c     >  +il5(mgs)*(cracip(mgs)+cipacr(mgs))
     : +il5(mgs)*(ciacrf(mgs) ) ! +ciracr(mgs)+cipacr(mgs))
c     >  +csacwf(mgs)
     :  +cwacsf(mgs)+cracsf(mgs)
c     >  +csacrf(mgs)
c     >  +cglacwf(mgs)
c     :  +cwacglf(mgs)
c     >  +cglacrf(mgs)
c     :  +cracglf(mgs)
c     >  +cgmacwf(mgs)
c     :  +cwacgmf(mgs)
c     >  +cgmacrf(mgs)
c     :  +cracgmf(mgs)
c     >  +cghacwf(mgs)
c     :  +cwacghf(mgs)
c     >  +cracghf(mgs)
c     >  +cfacwf(mgs)
c     >  +cfacs(mgs)+cfaci(mgs)
c     >  +cfacir(mgs)
c     >  +cfacip(mgs)*il5(mgs)
      IF ( iacr .ge. 1 ) THEN
      
      ENDIF
      pcfwd(mgs) =
     >  -chcnf(mgs)
!     >  +cfshr(mgs) ! + cfshrp(mgs)
     >  +(1-il5(mgs))*cfmlr(mgs)
     >  +il5(mgs)*cfsbv(mgs)
     >  -cwacfgl(mgs)-cwacfgm(mgs)-cwacfgh(mgs)
c     >  -cfmul1(mgs)
      end do
c
c  Hail
c
      do mgs = 1,ngscnt
      pchwi(mgs) =
     >   chcnf(mgs)+chcngh(mgs)
     >  +chcngm(mgs)+chcngl(mgs)
c     >  +il5(mgs)*(chdpv(mgs))
c     >  +chacr(mgs)+chacw(mgs)
c     >  +chacs(mgs)+chaci(mgs)
c     >  +chacir(mgs)+chacip(mgs)*il5(mgs)
     :  +cwacglf(mgs)
     :  +cracglf(mgs)
     :  +cwacgmf(mgs)
     :  +cracgmf(mgs)
     >  +cracghf(mgs)
     :  +cwacghf(mgs)
      pchwd(mgs) =
c     >   chshr(mgs) +
     >  (1-il5(mgs))*chmlr(mgs)
     >  + il5(mgs)*chsbv(mgs)
     >  -chlcnh(mgs)
c     >  -chmul1(mgs)
      end do
c
c
c  Hail-Large
c
      do mgs = 1,ngscnt
      pchli(mgs) =
     >   chlcnh(mgs)
c     >  +il5(mgs)*(chldpv(mgs))
c     >  +chlacr(mgs)+chlacw(mgs)
c     >  +chlacs(mgs)+chlaci(mgs)
c     >  +chlacir(mgs)+chlacip(mgs)*il5(mgs)
      pchld(mgs) =
c     >   chlshr(mgs) +
     >  (1-il5(mgs))*chlmlr(mgs)
     >  +il5(mgs)*chlsbv(mgs)
c    >  -chlmul1(mgs)
      end do

      ENDIF ! (ipconc .ge. 5 )
c
c  Balance and checks for continuity.....within machine precision...
c
      do mgs = 1,ngscnt
      pctot(mgs)   = pccwi(mgs) +pccwd(mgs) +
     >               pccii(mgs) +pccid(mgs) +
     >               pcipi(mgs) +pcipd(mgs) +
     >               pciri(mgs) +pcird(mgs) +
     >               pcrwi(mgs) +pcrwd(mgs) +
     >               pcswi(mgs) +pcswd(mgs) +
     >               pcgli(mgs) +pcgld(mgs) +
     >               pcgmi(mgs) +pcgmd(mgs) +
     >               pcghi(mgs) +pcghd(mgs) +
     >               pcfwi(mgs) +pcfwd(mgs) +
     <               pchwi(mgs) +pchwd(mgs)
      end do
c
c
      ENDIF ! ( ipconc .ge. 1 )
c
c
c
c
c
c  GOGO
c  production terms
c
c
c
c  Vapor
c
      do mgs = 1,ngscnt
      pqwvi(mgs) = 
     >  -qrcev(mgs)
     >  +il5(mgs)*(-qhsbv(mgs)  -qssbv(mgs) -qfsbv(mgs)
     >             -qglsbv(mgs) -qgmsbv(mgs) -qghsbv(mgs))
     >  -il5(mgs)*qisbv(mgs)
     >  -il5(mgs)*qipsbv(mgs)
     >  -il5(mgs)*qirsbv(mgs)
     >  -il5(mgs)*qhlsbv(mgs)
      pqwvd(mgs) =  
     >   il5(mgs)*(-qiint(mgs)
     >  -qhdpv(mgs) -qsdpv(mgs) -qfdpv(mgs)
     >  -qgldpv(mgs) -qgmdpv(mgs) -qghdpv(mgs))
     >  -il5(mgs)*qidpv(mgs)
     >  -il5(mgs)*qirdpv(mgs)
     >  -il5(mgs)*qipdpv(mgs)
     >  -il5(mgs)*qhldpv(mgs)
      end do
c
c  Cloud ice
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      pqcii(mgs) =  
     >   il5(mgs)*qicicnt(mgs) 
     >  +il5(mgs)*qidpv(mgs) 
     >  +il5(mgs)*(qiacwi(mgs)+qwacii(mgs))
     >  +il5(mgs)*(qwfrzc(mgs)+qwctfzc(mgs))
     >  +il5(mgs)*(qicichr(mgs))
     >  +qgmmul1(mgs)+qglmul1(mgs)+qghmul1(mgs)
     >  +qhmul1(mgs)+qfmul1(mgs)+qsmul1(mgs)+qhlmul1(mgs)
      pqcid(mgs) =  
     >   il5(mgs)*(-qscni(mgs) -qwaci(mgs) - qscnvi(mgs)
     >  -qraci(mgs) -qsaci(mgs)) 
     >  -qglaci(mgs) -qgmaci(mgs) -qghaci(mgs)
     >  -qfaci(mgs)  -qhaci(mgs)
     >  -il5(mgs)*qiraci(mgs)
     >  -il5(mgs)*qipaci(mgs)
     >  +il5(mgs)*qisbv(mgs) 
     >  -qhlaci(mgs)
     <  +(1.-il5(mgs))*qimlr(mgs)
      end do
      end if
c
c  Plates
c
      do mgs = 1,ngscnt
      pqipi(mgs) =
     >   il5(mgs)*(qipipnt(mgs))
     >  +il5(mgs)*(qipacwip(mgs)+qwacipip(mgs))
     >  +il5(mgs)*(qipdpv(mgs))
     >  +il5(mgs)*(qwfrzp(mgs)+qwctfzp(mgs))
     >  +il5(mgs)*(qipiphr(mgs))
c      IF ( pqipi(mgs) .gt. 0.01 ) THEN
c       write(iunit,*) 
c     :  'pqipi',igs(mgs),jgs,kgs(mgs),qipipnt(mgs),qipacwip(mgs),
c     :  qwacipip(mgs),qipdpv(mgs),qwfrzp(mgs),qwctfzp(mgs),qipiphr(mgs)
c      ENDIF
      pqipd(mgs) =
     >   il5(mgs)*qipsbv(mgs)
     >  -il5(mgs)*qscnip(mgs)
     >  -il5(mgs)*qsacip(mgs)
     >  -il5(mgs)*qfacip(mgs)
     >  -il5(mgs)*qhacip(mgs)
     >  -il5(mgs)*qglacip(mgs)
     >  -il5(mgs)*qgmacip(mgs)
     >  -il5(mgs)*qghacip(mgs)
     >  -il5(mgs)*qiracip(mgs)
     >  -il5(mgs)*qiacip(mgs)
     >  -il5(mgs)*qracip(mgs)
     >  -il5(mgs)*qwacip(mgs)
     >  -il5(mgs)*qhlacip(mgs)
     >  +(1.-il5(mgs))*qipmlr(mgs)
      end do
c
c
c  Rimed ice
c 
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      pqiri(mgs) =  
     >  +il5(mgs)*qirdpv(mgs) 
     >  +il5(mgs)*(qiacwir(mgs)+qwaciir(mgs))
     >  +il5(mgs)*(qiracwir(mgs)+qwacirir(mgs))
     >  +il5(mgs)*(qipacwir(mgs)+qwacipir(mgs))
      pqird(mgs) =  
     >   il5(mgs)*(-qscnir(mgs) -qwacir(mgs)
     >  -qracir(mgs) -qsacir(mgs)) 
     >  -qglacir(mgs) -qgmacir(mgs) -qghacir(mgs) 
     >  -qfacir(mgs)  -qhacir(mgs)
     >  -il5(mgs)*qiacir(mgs)
     >  -il5(mgs)*qipacir(mgs)
     >  +il5(mgs)*qirsbv(mgs) 
     >  +(1.-il5(mgs))*(qirmlw(mgs) + qirmlr(mgs))
     >  -qhlacir(mgs)
      end do
      end if
c
c  Cloud water
c 
      do mgs = 1,ngscnt
      pqcwi(mgs) =  (0.0) - (1.-il5(mgs))*qirmlw(mgs)
      pqcwd(mgs) = 
     >  il5(mgs)*(-qiacw(mgs)-qwfrzc(mgs)-qwctfzc(mgs))
     >  -il5(mgs)*(qicichr(mgs))
     >  -il5(mgs)*(qipiphr(mgs))
     >  -il5(mgs)*(qiracw(mgs))
     >  -il5(mgs)*(qipacw(mgs))
     >  -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs)
     >  -qglacw(mgs) -qgmacw(mgs) -qghacw(mgs) -qfacw(mgs) 
     >  -il5(mgs)*(qwfrzp(mgs)+qwctfzp(mgs))
     >  -qhlacw(mgs)
      end do
c
c  Rain
c
      do mgs = 1,ngscnt
      pqrwi(mgs) =  
     >   qracw(mgs) +qrcnw(mgs) 
     >  +(1-il5(mgs))*
     >   (-qhmlr(mgs)-qsmlr(mgs) -qglmlr(mgs)-qgmmlr(mgs)
     >    -qghmlr(mgs)-qfmlr(mgs)
     <    -qirmlr(mgs)-qipmlr(mgs)-qimlr(mgs)-qhlmlr(mgs))
     <    -qrshr(mgs)
      pqrwd(mgs) =  
     >   il5(mgs)*(-qiacr(mgs)-qrfrz(mgs)) 
     >  -qsacr(mgs)-qhacr(mgs)
     >  -qglacr(mgs)-qgmacr(mgs)-qghacr(mgs)
     >  -qfacr(mgs)
     >  -il5(mgs)*qipacr(mgs)
     >  -il5(mgs)*qiracr(mgs)
     >  +qrcev(mgs)
     >  -qhlacr(mgs)
      
      IF ( -pqrwd(mgs)*dtp .gt. qx(mgs,lr) ) THEN
       write(0,*) 'OUCH! pqrwd(mgs)*dtp .gt. qrw(mgs) ',pqrwd(mgs),qx(mgs,lr)
       STOP
      ENDIF
      
      IF ( .false. .and. ny .eq. 2 .and. qgmmlr(mgs) .lt. 0.0 ) THEN
        write(iunit,*) 'pqrwi = ', pqrwi(mgs), igs(mgs),jy,kgs(mgs)
        write(iunit,*) 'pqrwd = ', pqrwd(mgs)
        write(iunit,*) 'pcrwi = ', pcrwi(mgs)
        write(iunit,*) 'pcrwd = ', pcrwd(mgs)
        write(iunit,*) -qhmlr(mgs),-qsmlr(mgs) ,-qglmlr(mgs),-qgmmlr(mgs)
        write(iunit,*) -qghmlr(mgs),-qfmlr(mgs)
        write(iunit,*) 'qirmlr = ',-qirmlr(mgs),-qipmlr(mgs),-qimlr(mgs),-qhlmlr(mgs)
        write(iunit,*) 'chmlrr = ' ,-chmlrr(mgs),-csmlr(mgs) ,-cglmlrr(mgs),-cgmmlrr(mgs)
        write(iunit,*) -cghmlrr(mgs),-cfmlrr(mgs)
        write(iunit,*)  -cirmlr(mgs),-cipmlr(mgs),-chlmlrr(mgs)
        write(iunit,*) 'cgmmlr = ',-cgmmlr(mgs),-cgmmlrr(mgs)
        tmp = -qhmlr(mgs)-qsmlr(mgs) -qglmlr(mgs)-qgmmlr(mgs)
     >    -qghmlr(mgs)-qfmlr(mgs)
     <    -qirmlr(mgs)-qipmlr(mgs)-qimlr(mgs)-qhlmlr(mgs)
        tmp2 = (-chmlrr(mgs)-csmlr(mgs) -cglmlrr(mgs)-cgmmlrr(mgs)
     >    -cghmlrr(mgs)-cfmlrr(mgs)
     <    -cirmlr(mgs)-cipmlr(mgs)-chlmlrr(mgs))/rz
        write(iunit,*) ' qmltr = ',tmp, tmp2
        IF ( tmp2 .ne. 0. ) THEN
        write(iunit,*) 'new drop size = ', 1.e6*(cwc1*rho0(mgs)*tmp/(tmp2))**(1./3.)
        write(iunit,*) 'from gm: ',1.e6*(cwc1*rho0(mgs)*qgmmlr(mgs)/(cgmmlrr(mgs)))**(1./3.)
        write(iunit,*) 'rho,etc = ', rho0(mgs),qgmmlr(mgs),xdn(mgs,lr),vmlt
        ENDIF
        write(iunit,*) 'qr,cr,qgl,qgm,qgh = ',qx(mgs,lr),cx(mgs,lr),qx(mgs,lgl),qx(mgs,lgm),qx(mgs,lgh)
        write(iunit,*) 'rdia = ',1.e6*xdia(mgs,lr,1)
        write(iunit,*) 'rz, dmlt,vmlt,ihmlt = ',rz,dmlt,vmlt,ihmlt


      write(iunit,*)      'Rain '
c
      write(iunit,*)      'qracw = ',qracw(mgs)
      write(iunit,*)      'qrcnw = ',qrcnw(mgs)
      write(iunit,*)       -(1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qsmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qglmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qgmmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qghmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qfmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qirmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qipmlr(mgs)
      write(iunit,*)       -qrshr(mgs)
      write(iunit,*)       'pqrwi = ', pqrwi(mgs)    
      write(iunit,*)        -il5(mgs)*qiacr(mgs)
      write(iunit,*)        -il5(mgs)*qiracr(mgs)
      write(iunit,*)        -il5(mgs)*qipacr(mgs)
      write(iunit,*)        -il5(mgs)*qrfrz(mgs)
      write(iunit,*)        -qsacr(mgs)
      write(iunit,*)        -qhacr(mgs)
      write(iunit,*)        -qglacr(mgs)
      write(iunit,*)        -qgmacr(mgs)
      write(iunit,*)        -qghacr(mgs)
      write(iunit,*)        -qfacr(mgs)
      write(iunit,*)        qrcev(mgs)
      write(iunit,*)       'pqrwd = ', pqrwd(mgs) 
c
      write(iunit,*)  'Rain concentration'
      write(iunit,*)  pcrwi(mgs) 
      write(iunit,*)    crcnw(mgs)
      write(iunit,*)    1-il5(mgs)
      write(iunit,*)   -chmlr(mgs),-csmlr(mgs),-cglmlr(mgs)
      write(iunit,*)    -cgmmlr(mgs),-cghmlr(mgs),-cfmlr(mgs)
      write(iunit,*)    -cirmlr(mgs),-cipmlr(mgs),-chlmlr(mgs)
      write(iunit,*)     -crshr(mgs)
      write(iunit,*)  pcrwd(mgs) 
      write(iunit,*)    il5(mgs)
      write(iunit,*)   -ciacr(mgs),-crfrz(mgs) ! ,-cipacr(mgs)
      write(iunit,*)   -csacr(mgs),-chacr(mgs),-chlacr(mgs)
      write(iunit,*)   -cglacr(mgs),-cgmacr(mgs),-cghacr(mgs)
      write(iunit,*)   -cfacr(mgs)
      write(iunit,*)   +crcev(mgs)
      write(iunit,*)   cracr(mgs)
      ENDIF
      end do
c
c  Snow
c
      do mgs = 1,ngscnt
      pqswi(mgs) =  
     >   il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs) + qscnvi(mgs))
     >  +qsacws(mgs)+qsacrs(mgs)+qracss(mgs)
     >  +il5(mgs)*qsacir(mgs)
     >  +il5(mgs)*qsacip(mgs)
     >  +il5(mgs)*qscnir(mgs)
     >  +il5(mgs)*qscnip(mgs)
     >  +il5(mgs)*(qiacir(mgs)+qiraci(mgs))
     >  +il5(mgs)*(qiacip(mgs)+qipaci(mgs))
     >  +il5(mgs)*(qiracip(mgs)+qipacir(mgs))
      pqswd(mgs) = 
     >  -qwacs(mgs)-qracs(mgs)-qfacs(mgs)
     >  -qhacs(mgs)-qglacs(mgs)-qgmacs(mgs)-qghacs(mgs)
     >  +(1-il5(mgs))*qsmlr(mgs) + qsshr(mgs) ! + qsshrp(mgs)
     >  +il5(mgs)*(qssbv(mgs))
     >  -qsmul1(mgs)
     >  -qhlacs(mgs)
      end do 
      
c C$PAR CRITICAL SECTION      
c      DO mgs=1,ngscnt
c       IF ( temcg(mgs) .lt. 1.0 .and. temcg(mgs) .gt. -15 .and.
c     :       pqswi(mgs) .gt. 1.0e-9 ) THEN
c        write(iunit,*) igs(mgs),jy,kgs(mgs),temcg(mgs)
c        write(iunit,*) qx(mgs,ls),qx(mgs,lr)
c        write(iunit,*) pqswi(mgs),pqswd(mgs),pqswi(mgs)+pqswd(mgs)
c        write(iunit,*) qscni(mgs),qsaci(mgs),qsdpv(mgs)
c        write(iunit,*) qsacws(mgs),qsacrs(mgs)
c        write(iunit,*) qsacir(mgs),qsacip(mgs)
c        write(iunit,*) qscnir(mgs),qscnip(mgs)
c        write(iunit,*) qiacir(mgs),qiraci(mgs)
c        write(iunit,*) qiacip(mgs),qipaci(mgs)
c        write(iunit,*) qiracip(mgs),qipacir(mgs)       
c        write(iunit,*) -qwacs(mgs),-qracs(mgs),-qfacs(mgs)
c        write(iunit,*) -qhacs(mgs),-qglacs(mgs),-qgmacs(mgs),-qghacs(mgs)
c        write(iunit,*) (1-il5(mgs))*qsmlr(mgs), qsshr(mgs), qsshrp(mgs)
c        write(iunit,*) il5(mgs)*(qssbv(mgs)), -qhlacs(mgs)
c       END IF
c       
c      END DO
c C$PAR END CRITICAL SECTION      

c
c  Graupel-L
c
      do mgs = 1,ngscnt
      pqgli(mgs) =  
     >  +il5(mgs)*(qgldpv(mgs))
     >  +il5(mgs)*(qipacwgl(mgs)+qwacipgl(mgs))
     >  +il5(mgs)*(qiracwgl(mgs)+qwacirgl(mgs))
     >  +il5(mgs)*(qiacwgl(mgs)+qwacigl(mgs))
     >  +qsacwgl(mgs)+qwacsgl(mgs)
     >  +qsacrgl(mgs)+qracsgl(mgs)
     >  +qglacwgl(mgs)
     >  +qglacrgl(mgs)
     >  +qgmacwgl(mgs)+qwacgmgl(mgs)
     >  +qgmacrgl(mgs)+qracgmgl(mgs)
     >  +qghacwgl(mgs)+qwacghgl(mgs)
     >  +qghacrgl(mgs)+qracghgl(mgs)
     >  +qfacwgl(mgs)+qwacfgl(mgs)
     >  +qglacs(mgs)+qglaci(mgs)
     >  +qglacir(mgs)
     >  +il5(mgs)*qglacip(mgs)
      pqgld(mgs) =  
     >   qglshr(mgs) !+ qglshrp(mgs)
     >  +(1-il5(mgs))*qglmlr(mgs) 
     >  +il5(mgs)*qglsbv(mgs)
     >  -qracglgh(mgs)-qwacglgh(mgs)
     >  -qracglgm(mgs)-qwacglgm(mgs)
     >  -qracglf(mgs)-qwacglf(mgs)
     >  -qhcngl(mgs)
     >  -qglmul1(mgs)
      IF ( (-pqgld(mgs) - qglacwgl(mgs) - qglacrgl(mgs) )*dtp .gt. qx(mgs,lgl) .and. 
     :          qx(mgs,lgl) .gt. qxmin(lgl) ) THEN
       write(0,*) 'OUCH! pqgld(mgs)*dtp .gt. qgl(mgs) ',pqgld(mgs),qx(mgs,lgl)
       write(0,*) ' pqgld(mgs)*dt = ',pqgld(mgs)*dtp
       write(0,*) ' pqgli(mgs)*dt = ',pqgli(mgs)*dtp
       write(0,*)    qglshr(mgs),temcg(mgs) 
       write(0,*)   +(1-il5(mgs))*qglmlr(mgs) 
       write(0,*)   +il5(mgs)*qglsbv(mgs)
       write(0,*)   -qracglgh(mgs)
       write(0,*)   -qwacglgh(mgs)
       write(0,*)   -qracglgm(mgs)
       write(0,*)   -qwacglgm(mgs)
       write(0,*)   -qracglf(mgs)
       write(0,*)   -qwacglf(mgs)
       write(0,*)   -qhcngl(mgs)
       write(0,*)   -qglmul1(mgs)
       
      IF ( -(pqgli(mgs) + pqgld(mgs))*dtp .gt. qx(mgs,lgl) ) STOP 
      
       
      ENDIF
      end do
c
c  Graupel-M
c
      do mgs = 1,ngscnt
      pqgmi(mgs) =
!     >  +(1-il5(mgs))*(-qglmlr(mgs))
     >  +il5(mgs)*(qgmdpv(mgs))
c     >  +il5(mgs)*(qipacwgm(mgs)+qwacipgm(mgs))
c     >  +il5(mgs)*(qiracwgm(mgs)+qwacirgm(mgs))
c     >  +il5(mgs)*(qiacwgm(mgs)+qwacigm(mgs))
     >  +qsacwgm(mgs)+qwacsgm(mgs)
     >  +qsacrgm(mgs)+qracsgm(mgs)
     >  +qglacwgm(mgs)+qwacglgm(mgs)
     >  +qglacrgm(mgs)+qracglgm(mgs)
     >  +qgmacwgm(mgs)
     >  +qgmacrgm(mgs)
     >  +qghacwgm(mgs)+qwacghgm(mgs)
     >  +qghacrgm(mgs)+qracghgm(mgs)
     >  +qfacwgm(mgs)+qwacfgm(mgs)
     >  +qgmacs(mgs)+qgmaci(mgs)
     >  +qgmacir(mgs)
     >  +il5(mgs)*qgmacip(mgs)
      pqgmd(mgs) =
     >   qgmshr(mgs) ! + qgmshrp(mgs)
     >  +(1-il5(mgs))*qgmmlr(mgs)
     >  +il5(mgs)*qgmsbv(mgs)
     >  -qracgmgh(mgs)-qwacgmgh(mgs)
     >  -qracgmgl(mgs)-qwacgmgl(mgs)
     >  -qracgmf(mgs)-qwacgmf(mgs)
     >  -qhcngm(mgs)
     >  -qgmmul1(mgs)
      end do
c
c  Graupel-H
c
      do mgs = 1,ngscnt
      pqghi(mgs) =  
!     >  +(1-il5(mgs))*(-qgmmlr(mgs))
     >  +il5(mgs)*(qghdpv(mgs))
c     >  +il5(mgs)*(qipacwgh(mgs)+qwacipgh(mgs))
c     >  +il5(mgs)*(qiracwgh(mgs)+qwacirgh(mgs))
c     >  +il5(mgs)*(qiacwgh(mgs)+qwacigh(mgs))
     >  +qsacwgh(mgs)+qwacsgh(mgs)
     >  +qsacrgh(mgs)+qracsgh(mgs)
     >  +qglacwgh(mgs)+qwacglgh(mgs)
     >  +qglacrgh(mgs)+qracglgh(mgs)
     >  +qgmacwgh(mgs)+qwacgmgh(mgs)
     >  +qgmacrgh(mgs)+qracgmgh(mgs)
     >  +qghacwgh(mgs)
     >  +qghacrgh(mgs)
     >  +qfacwgh(mgs)+qwacfgh(mgs)
     >  +qghacs(mgs)+qghaci(mgs)
     >  +qghacir(mgs)
     >  +il5(mgs)*qghacip(mgs)
      pqghd(mgs) =  
     >   qghshr(mgs) !+ qghshrp(mgs)
     >  +(1-il5(mgs))*qghmlr(mgs) 
     >  +il5(mgs)*qghsbv(mgs)
     >  -qracghgl(mgs)-qwacghgl(mgs)
     >  -qracghgm(mgs)-qwacghgm(mgs)
     >  -qracghf(mgs)-qwacghf(mgs)
     >  -qhcngh(mgs)
     >  -qghmul1(mgs)
      end do
c
c  frozen drops
c
      do mgs = 1,ngscnt
      pqfwi(mgs) =  
     >  +il5(mgs)*(qracif(mgs)+qiacrf(mgs)+qrfrz(mgs)+qfdpv(mgs))
     >  +il5(mgs)*(qracir(mgs)+qiracr(mgs))
     >  +il5(mgs)*(qracip(mgs)+qipacr(mgs))
     >  +qsacwf(mgs)+qwacsf(mgs)
     >  +qsacrf(mgs)+qracsf(mgs)
c     >  +qglacwf(mgs)+qwacglf(mgs)
c     >  +qglacrf(mgs)+qracglf(mgs)
c     >  +qgmacwf(mgs)+qwacgmf(mgs)
c     >  +qgmacrf(mgs)+qracgmf(mgs)
c     >  +qghacwf(mgs)+qwacghf(mgs)
c     >  +qghacrf(mgs)+qracghf(mgs)
     >  +qfacrf(mgs)+qfacwf(mgs)
     >  +qfacs(mgs)+qfaci(mgs)
     >  +qfacir(mgs)
     >  +il5(mgs)*qfacip(mgs)
      pqfwd(mgs) =  
     >  -qhcnf(mgs) 
     >  +qfshr(mgs) ! + qfshrp(mgs)
     >  +(1-il5(mgs))*qfmlr(mgs) 
     >  +il5(mgs)*qfsbv(mgs)
     >  -qwacfgl(mgs)-qwacfgm(mgs)-qwacfgh(mgs)
     >  -qfmul1(mgs)
      end do
c
c  Hail
c
      do mgs = 1,ngscnt
      pqhwi(mgs) = 
     >   qhcnf(mgs)+qhcngh(mgs)
     >  +qhcngm(mgs)+qhcngl(mgs)
     >  +il5(mgs)*(qhdpv(mgs))
     >  +qhacr(mgs)+qhacw(mgs)
     >  +qhacs(mgs)+qhaci(mgs)
     >  +qhacir(mgs)
     >  +il5(mgs)*qhacip(mgs)
     >  +qglacwf(mgs)+qwacglf(mgs)
     >  +qglacrf(mgs)+qracglf(mgs)
     >  +qgmacwf(mgs)+qwacgmf(mgs)
     >  +qgmacrf(mgs)+qracgmf(mgs)
     >  +qghacwf(mgs)+qwacghf(mgs)
     >  +qghacrf(mgs)+qracghf(mgs)
      pqhwd(mgs) =  
     >   qhshr(mgs) 
     >  +(1-il5(mgs))*qhmlr(mgs) 
     >  +il5(mgs)*qhsbv(mgs)
     >  -qhlcnh(mgs)
     >  -qhmul1(mgs)
      end do
c
c  Hail-Large
c
      do mgs = 1,ngscnt
      pqhli(mgs) =
     >   qhlcnh(mgs)
     >  +il5(mgs)*(qhldpv(mgs))
     >  +qhlacr(mgs)+qhlacw(mgs)
     >  +qhlacs(mgs)+qhlaci(mgs)
     >  +qhlacir(mgs)
     >  +il5(mgs)*qhlacip(mgs)
      pqhld(mgs) =
     >   qhlshr(mgs)
     >  +(1-il5(mgs))*qhlmlr(mgs)
     >  +il5(mgs)*qhlsbv(mgs)
     >  -qhlmul1(mgs)
      end do
c
c
c
c
      printrates = .false.
      
c      IF ( ipconc .ge. 1 ) THEN
c      IF ( igs(mgs)==19 .and. jy==18 .and. kgs(mgs)==34 .and. 
c     :   cx(mgs,li) + dtp*(pccii(mgs)+pccid(mgs)) .gt. 1000.   ) THEN
c       printrates = .true.
c      ENDIF
c      ENDIF
      
      if ( nstep/1*1 .eq. nstep ) then
      do mgs = 1,ngscnt
c
      ptotal(mgs) = 0.0d0
      ptotal(mgs) = ptotal(mgs) 
     >  + pqwvi(mgs) + pqwvd(mgs)
     >  + pqcwi(mgs) + pqcwd(mgs)
     >  + pqcii(mgs) + pqcid(mgs)
     >  + pqipi(mgs) + pqipd(mgs)
     >  + pqiri(mgs) + pqird(mgs)
     >  + pqrwi(mgs) + pqrwd(mgs)
     >  + pqswi(mgs) + pqswd(mgs)
     >  + pqgli(mgs) + pqgld(mgs)
     >  + pqgmi(mgs) + pqgmd(mgs)
     >  + pqghi(mgs) + pqghd(mgs)
     >  + pqfwi(mgs) + pqfwd(mgs)
     >  + pqhwi(mgs) + pqhwd(mgs)
     <  + pqhli(mgs) + pqhld(mgs)
c      

      if ( ( ndebug .ge. 0 .and. abs(ptotal(mgs)) .gt. eqtot )
c     :    .or. pqswi(mgs)*dtp .gt. 1.e-3
c     :     .or. dtp*(pqrwi(mgs)+pqrwd(mgs)) .gt. 10.0e-3 
c     :     .or. dtp*(pccii(mgs)+pccid(mgs)) .gt. 1.e7 
c     :     .or. dtp*(pcipi(mgs)+pcipd(mgs)) .gt. 1.e7 
     :  .or.  .not. (ptotal(mgs) .lt. 1.0 .and.
     :            ptotal(mgs) .gt. -1.0)
     :  .or.  printrates ) then
      print*,'YIKES!','ptotal',nstep,mgs,igs(mgs),jgs,
     >       kgs(mgs),ptotal(mgs)
     
      print*, t7(igs(mgs),jgs,kgs(mgs))
      print*, cx(mgs,li),cx(mgs,lip),cx(mgs,lc),cx(mgs,lr),0.5*xdia(mgs,lr,1)
      print*, qx(mgs,lc),qx(mgs,li),qx(mgs,lip),qx(mgs,lir),qx(mgs,lr)
      print*, xmas(mgs,lr),xmas(mgs,lr)*cx(mgs,lr)/rho0(mgs)
      print*, vtxbar(mgs,li,1),vtxbar(mgs,lc,1),eiw(mgs),vtxbar(mgs,lr,1)
      print*, xdia(mgs,li,1),xdia(mgs,lc,1),qcmxd(mgs)
      
      
      print*, 'li: ',xdia(mgs,li,1),xdia(mgs,li,2),xmas(mgs,li),qx(mgs,li),
     :         vtxbar(mgs,li,1)
      
      
      print*, 'rain cx,xv : ',cx(mgs,lr),xv(mgs,lr)
      print*, temg(mgs)
      print*, pqwvi(mgs) ,pqwvd(mgs)
      print*, pqcwi(mgs) ,pqcwd(mgs)
      print*, pqcii(mgs) ,pqcid(mgs)
      print*, pqipi(mgs) ,pqipd(mgs)
      print*, pqiri(mgs) ,pqird(mgs)
      print*, pqrwi(mgs) ,pqrwd(mgs)
      print*, pqswi(mgs) ,pqswd(mgs)
      print*, pqgli(mgs) ,pqgld(mgs)
      print*, pqgmi(mgs) ,pqgmd(mgs)
      print*, pqghi(mgs) ,pqghd(mgs)
      print*, pqfwi(mgs) ,pqfwd(mgs)
      print*, pqhwi(mgs) ,pqhwd(mgs)
      print*, pqhli(mgs) ,pqhld(mgs)
      print*,'END OF OUTPUT OF SOURCE AND SINK'

c
c  print production terms
c
      print*,  'Vapor'
c
      print*,  -qrcev(mgs)
      print*,  -il5(mgs)*qhsbv(mgs) 
      print*,  -il5(mgs)*qssbv(mgs) 
      print*,  -il5(mgs)*qfsbv(mgs)
      print*,  -il5(mgs)*qglsbv(mgs) 
      print*,  -il5(mgs)*qgmsbv(mgs) 
      print*,  -il5(mgs)*qghsbv(mgs)
      print*,  -il5(mgs)*qisbv(mgs)
      print*,  -il5(mgs)*qipsbv(mgs)
      print*,  -il5(mgs)*qirsbv(mgs)
      print*,   pqwvi(mgs) 
      print*,  -il5(mgs)*qiint(mgs)
      print*,  -il5(mgs)*qhdpv(mgs) 
      print*,  -il5(mgs)*qsdpv(mgs) 
      print*,  -il5(mgs)*qfdpv(mgs)
      print*,  -il5(mgs)*qgldpv(mgs) 
      print*,  -il5(mgs)*qgmdpv(mgs) 
      print*,  -il5(mgs)*qghdpv(mgs)
      print*,  -il5(mgs)*qidpv(mgs)
      print*,  -il5(mgs)*qipdpv(mgs)
      print*,  -il5(mgs)*qirdpv(mgs)
      print*,   pqwvd(mgs)
c
      print*,  'Cloud ice'
c
      print*,  il5(mgs)*qicicnt(mgs)
      print*,  il5(mgs)*qidpv(mgs)
      print*,  il5(mgs)*qiacwi(mgs)
      print*,  il5(mgs)*qwfrz(mgs)
      print*,  il5(mgs)*qwctfz(mgs)
      print*,  il5(mgs)*qicichr(mgs)
      print*,  'pqcii = ', pqcii(mgs)  
      print*,  -il5(mgs)*qscni(mgs)
      print*,  -il5(mgs)*qwaci(mgs) 
      print*,  -il5(mgs)*qraci(mgs) 
      print*,  -il5(mgs)*qsaci(mgs)
      print*,  -il5(mgs)*qglaci(mgs) 
      print*,  -il5(mgs)*qgmaci(mgs) 
      print*,  -il5(mgs)*qghaci(mgs) 
      print*,  -il5(mgs)*qfaci(mgs)
      print*,  -il5(mgs)*qhaci(mgs)
      print*,  il5(mgs)*qisbv(mgs)
      print*,  'pqcid = ',pqcid(mgs) 
      print*,  ' Conc:'
      print*,  pccii(mgs),pccid(mgs)
      print*,  il5(mgs),cicint(mgs)
      print*,  cwacii(mgs),cwfrzc(mgs),cwctfzc(mgs)
      print*,  cicichr(mgs)
      print*,  chmul1(mgs),cfmul1(mgs),cghmul1(mgs)
      print*,  cgmmul1(mgs),cglmul1(mgs),csmul1(mgs)
c
      print*, 'Cloud ice plate'
c
      print*,  'pqipi = ',pqipi(mgs)
      print*,  il5(mgs)*qipipnt(mgs)
      print*,  il5(mgs)*qipdpv(mgs)
      print*,  il5(mgs)*qipacwip(mgs)
      print*,  il5(mgs)*qwfrzp(mgs)
      print*,  il5(mgs)*qwctfzp(mgs)
      print*,  il5(mgs)*qipiphr(mgs)
      print*,  'pqipd = ', pqipd(mgs) 
      print*,  -il5(mgs)*qscnip(mgs)
      print*,  -il5(mgs)*qwacip(mgs)
      print*,  -il5(mgs)*qracip(mgs)
      print*,  -il5(mgs)*qsacip(mgs)
      print*,  -il5(mgs)*qglacip(mgs)
      print*,  -il5(mgs)*qgmacip(mgs)
      print*,  -il5(mgs)*qghacip(mgs)
      print*,  -il5(mgs)*qfacip(mgs)
      print*,  -il5(mgs)*qhacip(mgs)
      print*,  -il5(mgs)*qiracip(mgs)
      print*,  -il5(mgs)*qiacip(mgs)
      print*,  +il5(mgs)*qipsbv(mgs)
      print*,  +(1.-il5(mgs))*qipmlr(mgs)
      print*,  ' Conc:'
      print*,  pcipi(mgs),pcipd(mgs)
      print*, il5(mgs),cipint(mgs)
      print*, cipdpv(mgs)
      print*, cwfrzp(mgs),cwctfzp(mgs)
      print*, cipiphr(mgs),cipacwip(mgs),cwacipip(mgs)

c
      print*,  'Rimed ice'
      print*,  'q,c,vt,dia = ',qx(mgs,lir),cx(mgs,lir),vtxbar(mgs,lir,1),xdia(mgs,lir,1)
      print*,  'pqiri = ',pqiri(mgs)
      print*,  il5(mgs)*qirdpv(mgs)
c      print*,  il5(mgs)*qirirhr(mgs)
      print*,  il5(mgs)*qiracwir(mgs)
      print*,  il5(mgs)*qwacirir(mgs)
      print*,  il5(mgs)*qiacwir(mgs)
      print*,  il5(mgs)*qwaciir(mgs)
      print*,  il5(mgs)*qipacwir(mgs)
      print*,  il5(mgs)*qwacipir(mgs)
      print*,  'pqird = ',pqird(mgs)
      print*,  -il5(mgs)*qscnir(mgs)
      print*,  -il5(mgs)*qwacir(mgs), qwacir(mgs), qwacirgl(mgs), qwacirir(mgs), qwacirgl(mgs) + qwacirir(mgs)
      print*,  -il5(mgs)*qracir(mgs)
      print*,  -il5(mgs)*qsacir(mgs)
      print*,  -il5(mgs)*qglacir(mgs)
      print*,  -il5(mgs)*qgmacir(mgs)
      print*,  -il5(mgs)*qghacir(mgs)
      print*,  -il5(mgs)*qfacir(mgs)
      print*,  -il5(mgs)*qhacir(mgs)
      print*,  -il5(mgs)*qhlacir(mgs)
      print*,  -il5(mgs)*qiacir(mgs)
      print*,  -il5(mgs)*qipacir(mgs)
c      print*,  il5(mgs)*qwacirc(mgs)
      print*,  il5(mgs)*qirsbv(mgs)
      print*,  (1-il5(mgs))*qirmlr(mgs)
      print*,  (1-il5(mgs))*qirmlw(mgs)
c
c
      print*,  'Cloud water'
c
      print*,  pqcwi(mgs) 
      print*,  -il5(mgs)*qiacw(mgs)
      print*,  -il5(mgs)*qwfrzc(mgs)
      print*,  -il5(mgs)*qwctfzc(mgs)
      print*,  -il5(mgs)*qwfrzp(mgs)
      print*,  -il5(mgs)*qwctfzp(mgs)
      print*,  -il5(mgs)*qiihr(mgs)
      print*,  -il5(mgs)*qicichr(mgs)
      print*,  -il5(mgs)*qipiphr(mgs)
      print*,  -qracw(mgs)
      print*,  -qsacw(mgs) 
      print*,  -qrcnw(mgs) 
      print*,  -qhacw(mgs)
      print*,  -qglacw(mgs) 
      print*,  -qgmacw(mgs) 
      print*,  -qghacw(mgs) 
      print*,  -qfacw(mgs)
      print*,  pqcwd(mgs) 
c
      print*,     'Rain '
c
      print*,     qracw(mgs)
      print*,     qrcnw(mgs)
      print*,      -(1-il5(mgs))*qhmlr(mgs)
      print*,      -(1-il5(mgs))*qsmlr(mgs)
      print*,      -(1-il5(mgs))*qglmlr(mgs)
      print*,      -(1-il5(mgs))*qgmmlr(mgs)
      print*,      -(1-il5(mgs))*qghmlr(mgs)
      print*,      -(1-il5(mgs))*qfmlr(mgs)
      print*,      -(1-il5(mgs))*qirmlr(mgs)
      print*,      -(1-il5(mgs))*qipmlr(mgs)
      print*,      -qrshr(mgs)
      print*,      'pqrwi = ', pqrwi(mgs)    
      print*,       -il5(mgs)*qiacr(mgs)
      print*,       -il5(mgs)*qiracr(mgs)
      print*,       -il5(mgs)*qipacr(mgs)
      print*,       -il5(mgs)*qrfrz(mgs)
      print*,       -qsacr(mgs)
      print*,       -qhacr(mgs)
      print*,       -qglacr(mgs)
      print*,       -qgmacr(mgs)
      print*,       -qghacr(mgs)
      print*,       -qfacr(mgs)
      print*,       qrcev(mgs)
      print*,      'pqrwd = ', pqrwd(mgs) 
c
      
      print*, 'Rain concentration'
      print*, pcrwi(mgs) 
      print*,   crcnw(mgs)
      print*,   1-il5(mgs)
      print*,  -chmlr(mgs),-csmlr(mgs),-cglmlr(mgs)
      print*,   -cgmmlr(mgs),-cghmlr(mgs),-cfmlr(mgs)
      print*,   -cirmlr(mgs),-cipmlr(mgs),-chlmlr(mgs)
      print*,    -crshr(mgs)
      print*, pcrwd(mgs) 
      print*,   il5(mgs)
      print*,  -ciacr(mgs),-crfrz(mgs) ! ,-cipacr(mgs)
      print*,  -csacr(mgs),-chacr(mgs),-chlacr(mgs)
      print*,  -cglacr(mgs),-cgmacr(mgs),-cghacr(mgs)
      print*,  -cfacr(mgs)
      print*,  +crcev(mgs)
      print*,  cracr(mgs)
c      print*,  -il5(mgs)*ciracr(mgs)


      print*,  'Snow'
c
      print*,       il5(mgs)*qscni(mgs)
      print*,       il5(mgs)*qscnir(mgs)
      print*,       il5(mgs)*qscnip(mgs)
      print*,       il5(mgs)*qsaci(mgs)
      print*,       il5(mgs)*qsacir(mgs)
      print*,       il5(mgs)*qsacip(mgs)
      print*,       il5(mgs)*qsdpv(mgs)
      print*,       il5(mgs)*qiacir(mgs)
      print*,       il5(mgs)*qiacip(mgs)
      print*,       il5(mgs)*qiraci(mgs)
      print*,       il5(mgs)*qipaci(mgs)
      print*,       il5(mgs)*qiracip(mgs)
      print*,       il5(mgs)*qipacir(mgs)
      print*,       qsacws(mgs)
      print*,       qsacrs(mgs)
      print*,       pqswi(mgs)  
      print*,       -qwacs(mgs)
      print*,       -qracs(mgs)
      print*,       -qfacs(mgs)
      print*,       -qhacs(mgs)
      print*,       -qglacs(mgs)
      print*,       -qgmacs(mgs)
      print*,       -qghacs(mgs)
      print*,      (1-il5(mgs))*qsmlr(mgs) 
      print*,      qsshr(mgs) 
!      print*,      qsshrp(mgs)
      print*,      il5(mgs)*(qssbv(mgs))
      print*,      pqswd(mgs)  
c
      print*,  'Graupel-L'
c
      print*,  il5(mgs)*(qgldpv(mgs))
      print*,  il5(mgs)*qiacwgl(mgs)
      print*,  il5(mgs)*qwacigl(mgs)
      print*,  il5(mgs)*qipacwgl(mgs)
      print*,  il5(mgs)*qwacipgl(mgs)
      print*,  il5(mgs)*qiracwgl(mgs)
      print*,  il5(mgs)*qwacirgl(mgs)
      print*,  qsacwgl(mgs),qwacsgl(mgs)
      print*,  qsacrgl(mgs),qracsgl(mgs)
      print*,  qglacwgl(mgs)
      print*,  qglacrgl(mgs)
      print*,  qgmacwgl(mgs),qwacgmgl(mgs)
      print*,  qgmacrgl(mgs),qracgmgl(mgs)
      print*,  qghacwgl(mgs),qwacghgl(mgs)
      print*,  qghacrgl(mgs),qracghgl(mgs)
      print*,  qfacwgl(mgs),qwacfgl(mgs)
      print*,  qglacs(mgs),qglaci(mgs)
      print*,  qglacip(mgs),qglacir(mgs)
      print*,  qglshr(mgs) 
!      print*,  qglshrp(mgs)
      print*,  (1-il5(mgs))*qglmlr(mgs)
      print*,  il5(mgs)*qglsbv(mgs)
      print*,  -qracglgh(mgs), -qwacglgh(mgs)
      print*,  -qracglgm(mgs), -qwacglgm(mgs)
      print*,  -qracglf(mgs), -qwacglf(mgs)
      print*,  pqgld(mgs) 
c
      print*,  'Graupel-M'
c
      print*,  il5(mgs)*(qgmdpv(mgs))
c      print*,  il5(mgs)*qiacwgm(mgs)
c      print*,  il5(mgs)*qwacigm(mgs)
c      print*,  il5(mgs)*qipacwgm(mgs)
c      print*,  il5(mgs)*qwacipgm(mgs)
c      print*,  il5(mgs)*qiracwgm(mgs)
c      print*,  il5(mgs)*qwacirgm(mgs)
      print*,  qsacwgm(mgs),qwacsgm(mgs)
      print*,  qsacrgm(mgs),qracsgm(mgs)
      print*,  qglacwgm(mgs),qwacglgm(mgs)
      print*,  qglacrgm(mgs),qracglgm(mgs)
      print*,  qgmacwgm(mgs)
      print*,  qgmacrgm(mgs)
      print*,  qghacwgm(mgs),qwacghgm(mgs)
      print*,  qghacrgm(mgs),qracghgm(mgs)
      print*,  qfacwgm(mgs),qwacfgm(mgs)
      print*,  qgmacs(mgs),qgmaci(mgs)
      print*,  qgmacir(mgs),qgmacip(mgs)
      print*,  qgmshr(mgs) 
!      print*,  qgmshrp(mgs)
      print*,  (1-il5(mgs))*qgmmlr(mgs)
      print*,  il5(mgs)*qgmsbv(mgs)
      print*,  -qracgmgh(mgs), -qwacgmgh(mgs)
      print*,  -qracgmgl(mgs), -qwacgmgl(mgs)
      print*,  -qracgmf(mgs), -qwacgmf(mgs)
      print*,  pqgmd(mgs) 
c
      print*,  'Graupel-H'
c
      print*,  il5(mgs)*(qghdpv(mgs))
c      print*,  il5(mgs)*qiacwgh(mgs)
c      print*,  il5(mgs)*qwacigh(mgs)
c      print*,  il5(mgs)*qipacwgh(mgs)
c      print*,  il5(mgs)*qwacipgh(mgs)
c      print*,  il5(mgs)*qiracwgh(mgs)
c      print*,  il5(mgs)*qwacirgh(mgs)
      print*,  qsacwgh(mgs),qwacsgh(mgs)
      print*,  qsacrgh(mgs),qracsgh(mgs)
      print*,  qglacwgh(mgs),qwacglgh(mgs)
      print*,  qglacrgh(mgs),qracglgh(mgs)
      print*,  qgmacwgh(mgs),qwacgmgh(mgs)
      print*,  qgmacrgh(mgs),qracgmgh(mgs)
      print*,  qghacwgh(mgs)
      print*,  qghacrgh(mgs)
      print*,  qfacwgh(mgs),qwacfgh(mgs)
      print*,  qghacs(mgs),qghaci(mgs)
      print*,  qghacip(mgs),qghacir(mgs)
      print*,  qghshr(mgs) 
!      print*,  qghshrp(mgs)
      print*,  (1-il5(mgs))*qghmlr(mgs)
      print*,  il5(mgs)*qghsbv(mgs)
      print*,  qhcngh(mgs)
      print*,  -qracghgl(mgs), -qwacghgl(mgs)
      print*,  -qracghgm(mgs), -qwacghgm(mgs)
      print*,  -qracghf(mgs), -qwacghf(mgs)
      print*,  pqghd(mgs) 
c
      print*,  'Frozen drops'
c
      print*,  qx(mgs,lf)
      print*,  pqfwi(mgs)  
      print*,   pqfwd(mgs) 
      print*,  il5(mgs)*qracif(mgs)
      print*,  il5(mgs)*qiacrf(mgs)
      print*,  il5(mgs)*qracir(mgs)
      print*,  il5(mgs)*qiracr(mgs)
      print*,  il5(mgs)*qracip(mgs)
      print*,  il5(mgs)*qipacr(mgs)
      print*,  il5(mgs)*qrfrz(mgs) 
      print*,  il5(mgs)*qfdpv(mgs)
      print*,  qsacwf(mgs),qwacsf(mgs)
      print*,  qsacrf(mgs),qracsf(mgs)
      print*,  qglacwf(mgs),qwacglf(mgs)
      print*,  qglacrf(mgs),qracglf(mgs)
      print*,  qgmacwf(mgs),qwacgmf(mgs)
      print*,  qgmacrf(mgs),qracgmf(mgs)
      print*,  qghacwf(mgs),qwacghf(mgs)
      print*,  qghacrf(mgs),qracghf(mgs)
      print*,  qfacrf(mgs),qfacwf(mgs)
      print*,  qfacs(mgs),qfaci(mgs)
      print*,  qfacip(mgs),qfacir(mgs)
      print*,  -qhcnf(mgs)
      print*,  qfshr(mgs) !+ qfshrp(mgs)
      print*,  (1-il5(mgs))*qfmlr(mgs)
      print*,  il5(mgs)*qfsbv(mgs)
      print*,   -qwacfgl(mgs)
      print*,  -qwacfgm(mgs)
      print*,  -qwacfgh(mgs)
c
      print*,  'Hail'
c
      print*,  qhcnf(mgs)
      print*,  qhcngh(mgs)
      print*,  il5(mgs)*(qhdpv(mgs))
      print*,  qhacr(mgs)
      print*,  qhacw(mgs)
      print*,  qhacs(mgs)
      print*,  qhaci(mgs)
      print*,  qhacir(mgs)
      print*,  qhacip(mgs)
      print*,  pqhwi(mgs)
      print*,  qhshr(mgs)
      print*,  (1-il5(mgs))*qhmlr(mgs)
      print*,  il5(mgs)*qhsbv(mgs)
      print*,  pqhwd(mgs) 
c
c  Balance and checks for continuity.....within machine precision...
c
c
      print*,'END OF OUTPUT OF SOURCE AND SINK'
      print*,'PTOTAL',ptotal(mgs)
      
      IF ( printrates ) STOP
c
      end if
c
      end do
c

      DO mgs = 1,ngscnt
      IF ( .not. (ptotal(mgs) .lt. eqtot .and.
     :            ptotal(mgs) .gt. -eqtot) ) THEN
c cmic$ guard
c C$PAR CRITICAL SECTION    
!! c$omp  critical
        iptotal = iptotal + 1
c        write(iunit,*) igs(mgs),jy,kgs(mgs),ptotal(mgs)
        ptotalmx = Max( ptotalmx, ptotal(mgs) )
        ptotalmn = Min( ptotalmn, ptotal(mgs) )
!! c$omp  end critical
c C$PAR END CRITICAL SECTION   
c cmic$ endguard
      END IF
      END DO

      end if ! ( nstep/12*12 .eq. nstep )
c
c  latent heating from phase changes (except qcw, qci cond, and evap)
c
      do mgs = 1,ngscnt
      pfrz(mgs) = 
     >  (1-il5(mgs))*
     >  (qfmlr(mgs)+qhmlr(mgs)+qhlmlr(mgs)+qsmlr(mgs)
     >  +qghmlr(mgs)+qgmmlr(mgs)+qglmlr(mgs))
     >  +il5(mgs)*(qrshr(mgs)
     >  +qrfrz(mgs)+qwfrz(mgs)
     >  +qwctfz(mgs)+qiihr(mgs)
     >  +qiacw(mgs)+qhacw(mgs)
     >  +qipacw(mgs)+qiracw(mgs)
     >  +qsacw(mgs)+qfacw(mgs)
     >  +qglacw(mgs)+qgmacw(mgs)+qghacw(mgs)
     >  +qiacr(mgs)+qhacr(mgs)
     >  +qipacr(mgs)+qiracr(mgs)
     >  +qsacr(mgs)+qfacr(mgs)+qhlacr(mgs)
     >  +qglacr(mgs)+qgmacr(mgs)+qghacr(mgs))
      psub(mgs) = 
     >   il5(mgs)*(qfdpv(mgs)+qfsbv(mgs)
     >  +qsdpv(mgs)+qssbv(mgs)+qhdpv(mgs)+qhsbv(mgs)
     >  +qhldpv(mgs)+qhlsbv(mgs)
     >  +qghdpv(mgs)+qghsbv(mgs)
     >  +qgmdpv(mgs)+qgmsbv(mgs)
     >  +qgldpv(mgs)+qglsbv(mgs)
     >  +qirdpv(mgs)+qirsbv(mgs) 
     >  +qipdpv(mgs)+qipsbv(mgs)
     >  +qidpv(mgs)+qisbv(mgs))
     >  +il5(mgs)*(qiint(mgs))
      pvap(mgs) = 
     >   qrcev(mgs)
      ptem(mgs) = 
     >  (cpi/pi0(mgs))*
     >  (felf(mgs)*pfrz(mgs)
     >  +fels(mgs)*psub(mgs) 
     >  +felv(mgs)*pvap(mgs))
c      IF ( .not. (ptem(mgs) .lt. 1.0e10 .and.
c     :            ptem(mgs) .gt. -1.0e10 ) ) THEN
c C$PAR CRITICAL SECTION
c       write(iunit,*) 'ptem out of range',ptem(mgs),igs(mgs),jy,kgs(mgs)
c       write(iunit,*) cpi,pi0(mgs),felf(mgs),pfrz(mgs),fels(mgs),psub(mgs)
c       write(iunit,*) felv(mgs),pvap(mgs),qx(mgs,lc)
c       write(iunit,*)
c      write(iunit,*) qfmlr(mgs),qhmlr(mgs),qhlmlr(mgs),qsmlr(mgs)
c      write(iunit,*) qghmlr(mgs),qgmmlr(mgs),qglmlr(mgs)
c      write(iunit,*) qrshr(mgs),qrfrz(mgs),qwfrz(mgs),qwctfz(mgs)
c      write(iunit,*) qiihr(mgs),qiacw(mgs),qhacw(mgs),qipacw(mgs)
c      write(iunit,*) qiracw(mgs),qsacw(mgs),qfacw(mgs)
c      write(iunit,*) qglacw(mgs),qgmacw(mgs),qghacw(mgs)
c      write(iunit,*) qiacr(mgs),qhacr(mgs)
c      write(iunit,*) qipacr(mgs),qiracr(mgs)
c      write(iunit,*) qsacr(mgs),qfacr(mgs),qhlacr(mgs)
c      write(iunit,*) qglacr(mgs),qgmacr(mgs),qghacr(mgs)
c C$PAR END CRITICAL SECTION
c      END IF
      thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)
      end do
c
c
c  production of space charge
c
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fsczz(mgs) = 0.0
      end do
      end if

c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscrw(mgs) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      fscrw(mgs) = scx(mgs,lr)/qx(mgs,lr)
      end if
      end do
      end if

c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fsccw(mgs) = 0.0
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
      fsccw(mgs) = scx(mgs,lc)/qx(mgs,lc)
      end if
      end do
      end if

c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscci(mgs) = 0.0
      if ( qx(mgs,li) .gt. qxmin(lc) ) then
      fscci(mgs) = scx(mgs,li)/qx(mgs,li)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscip(mgs) = 0.0
      if ( qx(mgs,lip) .gt. qxmin(li) ) then
      fscip(mgs) = scx(mgs,lip)/qx(mgs,lip)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscir(mgs) = 0.0
      if ( qx(mgs,lir) .gt. qxmin(li) ) then
      fscir(mgs) = scx(mgs,lir)/qx(mgs,lir)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscsw(mgs) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      fscsw(mgs) = scx(mgs,ls)/qx(mgs,ls)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscfw(mgs) = 0.0
      if ( qx(mgs,lf) .gt. qxmin(lf) ) then
      fscfw(mgs) = scx(mgs,lf)/qx(mgs,lf)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fschw(mgs) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      fschw(mgs) = scx(mgs,lh)/qx(mgs,lh)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fschl(mgs) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lh) ) then
      fschl(mgs) = scx(mgs,lhl)/qx(mgs,lhl)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscgl(mgs) = 0.0
      if ( qx(mgs,lgl) .gt. qxmin(lgl) ) then
      fscgl(mgs) = scx(mgs,lgl)/qx(mgs,lgl)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscgm(mgs) = 0.0
      if ( qx(mgs,lgm) .gt. qxmin(lgm) ) then
      fscgm(mgs) = scx(mgs,lgm)/qx(mgs,lgm)
      end if
      end do
      end if
c
c  space charge production
c  do transfer during conversion
c  accretion, etc
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      fscgh(mgs) = 0.0
      if ( qx(mgs,lgh) .gt. qxmin(lgh) ) then
      fscgh(mgs) = scx(mgs,lgh)/qx(mgs,lgh)
      end if
      end do
      end if
c
c
cSSSS
c
c
c  Cloud ice
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccimi(mgs) = 
     >   il5(mgs)*(fsccw(mgs)*qiacwi(mgs)+fsccw(mgs)*qwfrzc(mgs)
     >  +fsccw(mgs)*qwctfzc(mgs)+fsccw(mgs)*qicichr(mgs)
     >  +fscci(mgs)*qwacii(mgs))
     >  +fscgh(mgs)*qghmul1(mgs)+fscgm(mgs)*qgmmul1(mgs)
     >  +fscgl(mgs)*qglmul1(mgs)
     >  +fschw(mgs)*qhmul1(mgs)+fscfw(mgs)*qfmul1(mgs)
     >  +fschl(mgs)*qhlmul1(mgs)+fscsw(mgs)*qsmul1(mgs)
      psccimd(mgs) = 
     >   il5(mgs)*(-fscci(mgs)*qscni(mgs) -fscci(mgs)*qwaci(mgs)
     >  -fscci(mgs)*qraci(mgs) -fscci(mgs)*qsaci(mgs)) 
     >  -fscci(mgs)*qglaci(mgs) -fscci(mgs)*qgmaci(mgs) 
     >  -fscci(mgs)*qghaci(mgs) -fscci(mgs)*qfaci(mgs) 
     >  -fscci(mgs)*qhaci(mgs)
     >  -fscci(mgs)*qiraci(mgs)*il5(mgs)
     >  -fscci(mgs)*qipaci(mgs)*il5(mgs)
     >  -fscci(mgs)*qhlaci(mgs)
      end do
      end if
c
c  Cloud ice - plate
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscipmi(mgs) = 
     >   il5(mgs)*(fsccw(mgs)*qipacwip(mgs)+fsccw(mgs)*qwfrzp(mgs)
     >  +fsccw(mgs)*qwctfzp(mgs)+fsccw(mgs)*qipiphr(mgs)
     >  +fscip(mgs)*qwacipip(mgs))
c    >  +fschw(mgs)*qhmul1(mgs)+fscfw(mgs)*qfmul1(mgs)
c    >  +fscgh(mgs)*qghmul1(mgs)+fscgm(mgs)*qgmmul1(mgs)
c    >  +fscgl(mgs)*qglmul1(mgs)+fscsw(mgs)*qsmul1(mgs)
      pscipmd(mgs) = 
     >   il5(mgs)*(-fscip(mgs)*qscnip(mgs) -fscip(mgs)*qwacip(mgs)
     >  -fscip(mgs)*qracip(mgs) -fscip(mgs)*qsacip(mgs)
     >  -fscip(mgs)*qglacip(mgs) -fscip(mgs)*qgmacip(mgs) 
     >  -fscip(mgs)*qghacip(mgs) -fscip(mgs)*qfacip(mgs) 
     >  -fscip(mgs)*qhacip(mgs)
     >  -fscip(mgs)*qiacip(mgs)
     >  -fscip(mgs)*qiracip(mgs)
     >  -fscip(mgs)*qhlacip(mgs))
     >  +(1-il5(mgs))*fscip(mgs)*qipmlr(mgs)
      end do
      end if
c
c  rimed ice
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscirmi(mgs) = 
     >   il5(mgs)*fsccw(mgs)*(qipacwir(mgs))
     >  +il5(mgs)*fscip(mgs)*(qwacipir(mgs))
     >  +il5(mgs)*fsccw(mgs)*(qiracwir(mgs))
     >  +il5(mgs)*fscir(mgs)*(qwacirir(mgs))
     >  +il5(mgs)*fsccw(mgs)*(qiacwir(mgs))
     >  +il5(mgs)*fscci(mgs)*(qwaciir(mgs))
      pscirmd(mgs) = 
     >   il5(mgs)*(-fscir(mgs)*qscnir(mgs) -fscir(mgs)*qwacir(mgs)
     >  -fscir(mgs)*qracir(mgs)  -fscir(mgs)*qsacir(mgs)) 
     >  -fscir(mgs)*qglacir(mgs) -fscir(mgs)*qgmacir(mgs) 
     >  -fscir(mgs)*qghacir(mgs) -fscir(mgs)*qfacir(mgs) 
     >  -fscir(mgs)*qhacir(mgs)  -fscir(mgs)*qiacir(mgs)*il5(mgs)
     >  -fscir(mgs)*qipacir(mgs)*il5(mgs)
     >  +(1-il5(mgs))*fscir(mgs)*qirmlr(mgs)
     >  -fscir(mgs)*qhlacir(mgs)
      end do
      end if
c
c  Cloud water
c 
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccwmi(mgs) = 0.0
      psccwmd(mgs) = 
     >  -fsccw(mgs)*qracw(mgs)-fsccw(mgs)*qrcnw(mgs)
     >  +il5(mgs)*(-fsccw(mgs)*qiacw(mgs)
     >  -fsccw(mgs)*qwctfzc(mgs)-fsccw(mgs)*qwctfzp(mgs)
     >  -fsccw(mgs)*qwfrzc(mgs)-fsccw(mgs)*qwfrzp(mgs)
     >  -fsccw(mgs)*qicichr(mgs)-fsccw(mgs)*qipiphr(mgs)
     >  -fsccw(mgs)*qiracw(mgs)-fsccw(mgs)*qipacw(mgs))
     >  -fsccw(mgs)*qsacw(mgs) -fsccw(mgs)*qhacw(mgs)
     >  -fsccw(mgs)*qglacw(mgs) -fsccw(mgs)*qgmacw(mgs) 
     >  -fsccw(mgs)*qghacw(mgs) -fsccw(mgs)*qfacw(mgs) 
     >  -fsccw(mgs)*qhlacw(mgs)
      end do
      end if
c
c  Rain
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscrwmi(mgs) = 
     >    +fsccw(mgs)*qracw(mgs) +fsccw(mgs)*qrcnw(mgs) 
     >  +(1-il5(mgs))*
     >   (-fschw(mgs)*qhmlr(mgs)-fscsw(mgs)*qsmlr(mgs)
     >    -fscgl(mgs)*qglmlr(mgs)-fscgm(mgs)*qgmmlr(mgs)
     >    -fscgh(mgs)*qghmlr(mgs)-fscfw(mgs)*qfmlr(mgs)
     >    -fschl(mgs)*qhlmlr(mgs)-fscir(mgs)*qirmlr(mgs)
     >    -fscip(mgs)*qipmlr(mgs))
     >    -fscrw(mgs)*qhshr(mgs)-fscrw(mgs)*qsshr(mgs)
     >    -fscrw(mgs)*qglshr(mgs)-fscrw(mgs)*qgmshr(mgs)
     >    -fscrw(mgs)*qghshr(mgs)-fscrw(mgs)*qfshr(mgs)
     >    -fscrw(mgs)*qhlshr(mgs)
      pscrwmd(mgs) =  
     >  + il5(mgs)*(-fscrw(mgs)*qiacr(mgs)-fscrw(mgs)*qrfrz(mgs)) 
     >  -fscrw(mgs)*qsacr(mgs)-fscrw(mgs)*qhacr(mgs)
     >  -fscrw(mgs)*qglacr(mgs)-fscrw(mgs)*qgmacr(mgs)
     >  -fscrw(mgs)*qghacr(mgs)-fscrw(mgs)*qfacr(mgs)
     >  -fscrw(mgs)*qiracr(mgs)*il5(mgs)
     >  -fscrw(mgs)*qipacr(mgs)*il5(mgs)
     >  -fscrw(mgs)*qhlacr(mgs)
      end do
      end if
c
c  Snow
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscswmi(mgs) = 
     >  +il5(mgs)*(fscci(mgs)*qscni(mgs)+fscci(mgs)*qsaci(mgs)
     >  +fscir(mgs)*qscnir(mgs)+fscir(mgs)*qsacir(mgs)
     >  +fscip(mgs)*qscnip(mgs)+fscip(mgs)*qsacip(mgs)
     >  +fscci(mgs)*qiraci(mgs)+fscir(mgs)*qiacir(mgs)
     >  +fscci(mgs)*qipaci(mgs)+fscip(mgs)*qiacip(mgs)
     >  +fscip(mgs)*qiracip(mgs)+fscir(mgs)*qipacir(mgs))
     >  +fsccw(mgs)*qsacws(mgs)+fscrw(mgs)*qsacrs(mgs)
      pscswmd(mgs) = 
     >  -fscsw(mgs)*qwacs(mgs)-fscsw(mgs)*qracs(mgs)
     >  -fscsw(mgs)*qfacs(mgs)-fscsw(mgs)*qhacs(mgs)
     >  -fscsw(mgs)*qglacs(mgs)-fscsw(mgs)*qgmacs(mgs)
     >  -fscsw(mgs)*qghacs(mgs)
     >  +(1-il5(mgs))*fscsw(mgs)*qsmlr(mgs) +fscrw(mgs)*qsshr(mgs) 
     >  -fscsw(mgs)*qsmul1(mgs)
     >  -fscsw(mgs)*qhlacs(mgs)
      end do 
      end if
c
c  Graupel-L
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscglmi(mgs) = 
     >  +il5(mgs)*(fsccw(mgs)*qiacwgl(mgs)+fscci(mgs)*qwacigl(mgs))
     >  +il5(mgs)*(fsccw(mgs)*qiracwgl(mgs)+fscir(mgs)*qwacirgl(mgs))
     >  +il5(mgs)*(fsccw(mgs)*qipacwgl(mgs)+fscip(mgs)*qwacipgl(mgs))
     >  +fsccw(mgs)*qsacwgl(mgs)+fscsw(mgs)*qwacsgl(mgs)
     >  +fscrw(mgs)*qsacrgl(mgs)+fscsw(mgs)*qracsgl(mgs)
     >  +fsccw(mgs)*qglacwgl(mgs)
     >  +fscrw(mgs)*qglacrgl(mgs)
     >  +fsccw(mgs)*qgmacwgl(mgs)+fscgm(mgs)*qwacgmgl(mgs)
     >  +fscrw(mgs)*qgmacrgl(mgs)+fscgm(mgs)*qracgmgl(mgs)
     >  +fsccw(mgs)*qghacwgl(mgs)+fscgh(mgs)*qwacghgl(mgs)
     >  +fscrw(mgs)*qghacrgl(mgs)+fscgh(mgs)*qracghgl(mgs)
     >  +fsccw(mgs)*qfacwgl(mgs)+fscfw(mgs)*qwacfgl(mgs)
     >  +fscsw(mgs)*qglacs(mgs)+fscci(mgs)*qglaci(mgs)
     >  +fscir(mgs)*qglacir(mgs)+fscip(mgs)*qglacip(mgs)*il5(mgs)
      pscglmd(mgs) = 
     >  + fscrw(mgs)*qglshr(mgs)-fscgl(mgs)*qhcngl(mgs) 
     >  +(1-il5(mgs))*fscgl(mgs)*qglmlr(mgs) 
     >  -fscgl(mgs)*qracglgh(mgs)-fscgl(mgs)*qwacglgh(mgs)
     >  -fscgl(mgs)*qracglgm(mgs)-fscgl(mgs)*qwacglgm(mgs)
     >  -fscgl(mgs)*qracglf(mgs)-fscgl(mgs)*qwacglf(mgs)
     >  -fscgl(mgs)*qglmul1(mgs)
      end do
      end if
c
c  Graupel-M
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscgmmi(mgs) = 
c     >  +il5(mgs)*(fsccw(mgs)*qiacwgm(mgs)+fscci(mgs)*qwacigm(mgs))
c     >  +il5(mgs)*(fsccw(mgs)*qiracwgm(mgs)+fscir(mgs)*qwacirgm(mgs))
c     >  +il5(mgs)*(fsccw(mgs)*qipacwgm(mgs)+fscip(mgs)*qwacipgm(mgs))
     >  +fsccw(mgs)*qsacwgm(mgs)+fscsw(mgs)*qwacsgm(mgs)
     >  +fscrw(mgs)*qsacrgm(mgs)+fscsw(mgs)*qracsgm(mgs)
     >  +fsccw(mgs)*qglacwgm(mgs)+fscgl(mgs)*qwacglgm(mgs)
     >  +fscrw(mgs)*qglacrgm(mgs)+fscgl(mgs)*qracglgm(mgs)
     >  +fsccw(mgs)*qgmacwgm(mgs)
     >  +fscrw(mgs)*qgmacrgm(mgs)
     >  +fsccw(mgs)*qghacwgm(mgs)+fscgh(mgs)*qwacghgm(mgs)
     >  +fscrw(mgs)*qghacrgm(mgs)+fscgh(mgs)*qracghgm(mgs)
     >  +fsccw(mgs)*qfacwgm(mgs)+fscfw(mgs)*qwacfgm(mgs)
     >  +fscsw(mgs)*qgmacs(mgs)+fscci(mgs)*qgmaci(mgs)
     >  +fscir(mgs)*qgmacir(mgs)+fscip(mgs)*qgmacip(mgs)*il5(mgs)
      pscgmmd(mgs) = 
     >  + fscrw(mgs)*qgmshr(mgs)-fscgm(mgs)*qhcngm(mgs) 
     >  +(1-il5(mgs))*fscgm(mgs)*qgmmlr(mgs)
     >  -fscgm(mgs)*qracgmgh(mgs)-fscgm(mgs)*qwacgmgh(mgs)
     >  -fscgm(mgs)*qracgmgl(mgs)-fscgm(mgs)*qwacgmgl(mgs)
     >  -fscgm(mgs)*qracgmf(mgs)-fscgm(mgs)*qwacgmf(mgs)
     >  -fscgm(mgs)*qgmmul1(mgs)
      end do
      end if
c
c  Graupel-H
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscghmi(mgs) = 
c     >  +il5(mgs)*(fsccw(mgs)*qiacwgh(mgs)+fscci(mgs)*qwacigh(mgs))
c     >  +il5(mgs)*(fsccw(mgs)*qipacwgh(mgs)+fscip(mgs)*qwacipgh(mgs))
c     >  +il5(mgs)*(fsccw(mgs)*qiracwgh(mgs)+fscir(mgs)*qwacirgh(mgs))
     >  +fsccw(mgs)*qsacwgh(mgs)+fscsw(mgs)*qwacsgh(mgs)
     >  +fscrw(mgs)*qsacrgh(mgs)+fscsw(mgs)*qracsgh(mgs)
     >  +fsccw(mgs)*qglacwgh(mgs)+fscgl(mgs)*qwacglgh(mgs)
     >  +fscrw(mgs)*qglacrgh(mgs)+fscgl(mgs)*qracglgh(mgs)
     >  +fsccw(mgs)*qgmacwgh(mgs)+fscgm(mgs)*qwacgmgh(mgs)
     >  +fscrw(mgs)*qgmacrgh(mgs)+fscgm(mgs)*qracgmgh(mgs)
     >  +fsccw(mgs)*qghacwgh(mgs)
     >  +fscrw(mgs)*qghacrgh(mgs)
     >  +fsccw(mgs)*qfacwgh(mgs)+fscfw(mgs)*qwacfgh(mgs)
     >  +fscsw(mgs)*qghacs(mgs)+fscci(mgs)*qghaci(mgs)
     >  +fscir(mgs)*qghacir(mgs) +fscip(mgs)*qghacip(mgs)*il5(mgs)
      pscghmd(mgs) =  
     >  + fscrw(mgs)*qghshr(mgs)-fscgh(mgs)*qhcngh(mgs)
     >  +(1-il5(mgs))*fscgh(mgs)*qghmlr(mgs) 
     >  -fscgh(mgs)*qracghgl(mgs)-fscgh(mgs)*qwacghgl(mgs)
     >  -fscgh(mgs)*qracghgm(mgs)-fscgh(mgs)*qwacghgm(mgs)
     >  -fscgh(mgs)*qracghf(mgs)-fscgh(mgs)*qwacghf(mgs)
     >  -fscgh(mgs)*qghmul1(mgs)
      end do
      end if
c
c  frozen drops
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscfwmi(mgs) = 
     >  +il5(mgs)*(fscci(mgs)*qracif(mgs)+fscrw(mgs)*qiacrf(mgs)
     >  +fscrw(mgs)*qrfrz(mgs))
     >  +fsccw(mgs)*qsacwf(mgs)+fscsw(mgs)*qwacsf(mgs)
     >  +fscrw(mgs)*qsacrf(mgs)+fscsw(mgs)*qracsf(mgs)
     >  +fsccw(mgs)*qglacwf(mgs)+fscgl(mgs)*qwacglf(mgs)
     >  +fscrw(mgs)*qglacrf(mgs)+fscgl(mgs)*qracglf(mgs)
     >  +fsccw(mgs)*qgmacwf(mgs)+fscgm(mgs)*qwacgmf(mgs)
     >  +fscrw(mgs)*qgmacrf(mgs)+fscgm(mgs)*qracgmf(mgs)
     >  +fsccw(mgs)*qghacwf(mgs)+fscgh(mgs)*qwacghf(mgs)
     >  +fscrw(mgs)*qghacrf(mgs)+fscgh(mgs)*qracghf(mgs)
     >  +fscrw(mgs)*qfacrf(mgs)+fsccw(mgs)*qfacwf(mgs)
     >  +fscsw(mgs)*qfacs(mgs)+fscci(mgs)*qfaci(mgs)
     >  +il5(mgs)*fscrw(mgs)*qiracr(mgs)
     >  +il5(mgs)*fscir(mgs)*qracir(mgs)
     >  +il5(mgs)*fscrw(mgs)*qipacr(mgs)
     >  +il5(mgs)*fscip(mgs)*qracip(mgs)
     >  +fscir(mgs)*qfacir(mgs) +fscip(mgs)*qfacip(mgs)*il5(mgs)
      pscfwmd(mgs) =  
     >  -fscfw(mgs)*qhcnf(mgs)+fscrw(mgs)*qfshr(mgs) 
     >  +(1-il5(mgs))*fscfw(mgs)*qfmlr(mgs) 
     >  -fscfw(mgs)*qwacfgl(mgs)-fscfw(mgs)*qwacfgm(mgs)
     >  -fscfw(mgs)*qwacfgh(mgs)
     >  -fscfw(mgs)*qfmul1(mgs)
      end do
      end if
c
c  Hail
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschwmi(mgs) = 
     >  +fscfw(mgs)*qhcnf(mgs)+fscgh(mgs)*qhcngh(mgs) 
     >  +fscgm(mgs)*qhcngm(mgs)+fscgl(mgs)*qhcngl(mgs) 
     >  +fscrw(mgs)*qhacr(mgs)+fsccw(mgs)*qhacw(mgs)
     >  +fscsw(mgs)*qhacs(mgs)+fscci(mgs)*qhaci(mgs)
     >  +fscir(mgs)*qhacir(mgs)+fscip(mgs)*qhacip(mgs)*il5(mgs)
      pschwmd(mgs) = 
     >  -fschw(mgs)*qhlcnh(mgs)
     >  + fscrw(mgs)*qhshr(mgs) 
     >  +(1-il5(mgs))*fschw(mgs)*qhmlr(mgs) 
     >  -fschw(mgs)*qhmul1(mgs)
      end do
      end if
c
c  Hail-Large
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschlmi(mgs) = 
     >  +fschw(mgs)*qhlcnh(mgs)
     >  +fscrw(mgs)*qhlacr(mgs)+fsccw(mgs)*qhlacw(mgs)
     >  +fscsw(mgs)*qhlacs(mgs)+fscci(mgs)*qhlaci(mgs)
     >  +fscir(mgs)*qhlacir(mgs)+fscip(mgs)*qhlacip(mgs)*il5(mgs)
      pschlmd(mgs) = 
     >  + fscrw(mgs)*qhlshr(mgs)
     >  +(1-il5(mgs))*fschl(mgs)*qhlmlr(mgs)
     >  -fschl(mgs)*qhlmul1(mgs)
      end do
      end if
c
c
c Collisional Charging transfers (inductive/noninductive)
c
c
c  rain water
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscrwi(mgs) = 0.0
      pscrwd(mgs) = 0.0
      pscrwi(mgs) =  pscrwi(mgs)
     > + max(-schlacr(mgs),0.0)
     > + max(-schacr(mgs),0.0)
     > + max(-scfacr(mgs),0.0)
     > + max(-scsacr(mgs),0.0)
     > + max(-scglacr(mgs),0.0)
     > + max(-scgmacr(mgs),0.0)
     > + max(-scghacr(mgs),0.0)
      pscrwd(mgs) =  pscrwd(mgs)
     > + min(-schlacr(mgs),0.0)
     > + min(-schacr(mgs),0.0)
     > + min(-scfacr(mgs),0.0)
     > + min(-scsacr(mgs),0.0)
     > + min(-scglacr(mgs),0.0)
     > + min(-scgmacr(mgs),0.0)
     > + min(-scghacr(mgs),0.0)
      end do
      end if
c
c  cloud water
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccwi(mgs) = 0.0
      psccwd(mgs) = 0.0
      psccwi(mgs) =  psccwi(mgs)
     > + max(-schlacw(mgs),0.0)
     > + max(-schacw(mgs),0.0)
     > + max(-scfacw(mgs),0.0)
     > + max(-scsacw(mgs),0.0)
     > + max(-scglacw(mgs),0.0)
     > + max(-scgmacw(mgs),0.0)
     > + max(-scghacw(mgs),0.0)
      psccwd(mgs) =  psccwd(mgs)
     > + min(-schlacw(mgs),0.0)
     > + min(-schacw(mgs),0.0)
     > + min(-scfacw(mgs),0.0)
     > + min(-scsacw(mgs),0.0)
     > + min(-scglacw(mgs),0.0)
     > + min(-scgmacw(mgs),0.0)
     > + min(-scghacw(mgs),0.0)
     
      IF ( .not. ( psccwi(mgs) .lt. 1. .and. psccwi(mgs) .gt. -1. )
     :     .or. 
     :     .not. ( psccwd(mgs) .lt. 1. .and. psccwd(mgs) .gt. -1. )
     :      ) THEN
      
      write (iunit,*) 'psccwi,psccwd = ',psccwi(mgs),psccwd(mgs)
      write(iunit,*) cx(mgs,lc),xdia(mgs,lc,1),dezcomp(mgs)
      write(iunit,*) scsacw(mgs),cx(mgs,ls),vtxbar(mgs,ls,1),esw(mgs),qx(mgs,ls)
      write(iunit,*) scglacw(mgs),cx(mgs,lgl),vtxbar(mgs,lgl,1),eglw(mgs),qx(mgs,lgl)
      write(iunit,*) scgmacw(mgs),cx(mgs,lgm),vtxbar(mgs,lgm,1),egmw(mgs),qx(mgs,lgm)
      write(iunit,*) scghacw(mgs),cx(mgs,lgh),vtxbar(mgs,lgh,1),eghw(mgs),qx(mgs,lgh)
      write(iunit,*) scfacw(mgs),cx(mgs,lf),vtxbar(mgs,lf,1),efw(mgs),qx(mgs,lf)
      write(iunit,*) schacw(mgs),cx(mgs,lh),vtxbar(mgs,lh,1),ehw(mgs),qx(mgs,lh)
      write(iunit,*) schlacw(mgs),cx(mgs,lhl),vtxbar(mgs,lhl,1),ehlw(mgs),qx(mgs,lhl)
       STOP
      ENDIF
      
       qcwtmp(mgs) = Max(0.0,qx(mgs,lc) +
     >   dtp*(pqcwi(mgs)+pqcwd(mgs)) )

       IF ( qcwtmp(mgs) .gt. qxmin(lc) ) then
       
      IF ( ipconc .lt. 2 ) THEN
       ccwtmp = cwccn
       ccwtmp = Max(1.0,qcwtmp(mgs)*rho0(mgs)/xmas(mgs,lc))
      ELSE
 
       ccwtmp = Max(0.001, cx(mgs,lc) + dtp*(pccwi(mgs)+pccwd(mgs)) )
      
      ENDIF ! ( ipconc .lt. 2 )

      
       cpqc = (scx(mgs,lc)
     :   +dtp*(psccwi(mgs)+psccwd(mgs)
     :   +psccwmi(mgs)+psccwmd(mgs)))/ccwtmp
      
         IF ( Abs(cpqc) .gt. scwppmx .and.
     :       ( psccwi(mgs)+psccwd(mgs) ) .ne. 0.0 .and.
     :    Sign(1.0,scx(mgs,lc)) .eq.
     :        Sign(1.0,( psccwi(mgs)+psccwd(mgs) )) ) THEN
          
           cpqc0 = Sign(scwppmx,cpqc)
          
           scfac = (cpqc0*ccwtmp - scx(mgs,lc) -
     :       dtp*(psccwmi(mgs)+psccwmd(mgs)))/
     :              ( dtp*(psccwi(mgs)+psccwd(mgs)) )
c           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: scwfac = ',scfac
             write(0,*) 'OUCH: scwfac = ',scfac
           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)
        
c        IF ( scfac .lt. 1.0 .and. scfac .gt. 0.0 ) THEN
c        write(iunit,*) 'LOOK: sccwfac = ',scfac,cpqc,scx(mgs,lc),ccwtmp,
c     : kgs(mgs),temg(mgs)
c        ENDIF
         
        ELSE
          scfac = 1.0
        ENDIF
      
      ELSE
       scfac = 0.0
      ENDIF ! ( qcwtmp .gt. qxmin(lc) ) 
      
      IF ( scfac .lt. 1.0 .and. qx(mgs,lc) .gt. qxmin(lc) ) THEN
       schlacw(mgs) = scfac*schlacw(mgs)
       schacw(mgs) = scfac*schacw(mgs)
       scfacw(mgs) = scfac*scfacw(mgs)
       scsacw(mgs) = scfac*scsacw(mgs)
       scglacw(mgs) = scfac*scglacw(mgs)
       scgmacw(mgs) = scfac*scgmacw(mgs)
       scghacw(mgs) = scfac*scghacw(mgs)

      psccwi(mgs) = 0.0
      psccwd(mgs) = 0.0
      psccwi(mgs) =  psccwi(mgs)
     > + max(-schlacw(mgs),0.0)
     > + max(-schacw(mgs),0.0)
     > + max(-scfacw(mgs),0.0)
     > + max(-scsacw(mgs),0.0)
     > + max(-scglacw(mgs),0.0)
     > + max(-scgmacw(mgs),0.0)
     > + max(-scghacw(mgs),0.0)
      psccwd(mgs) =  psccwd(mgs)
     > + min(-schlacw(mgs),0.0)
     > + min(-schacw(mgs),0.0)
     > + min(-scfacw(mgs),0.0)
     > + min(-scsacw(mgs),0.0)
     > + min(-scglacw(mgs),0.0)
     > + min(-scgmacw(mgs),0.0)
     > + min(-scghacw(mgs),0.0)
      
      ENDIF
      end do
      
      
      end if
c
c  cloud ice
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psccii(mgs) = 0.0
      psccid(mgs) = 0.0
      psccii(mgs) =  psccii(mgs)
     > + max(-schlaci(mgs),0.0)
     > + max(-schaci(mgs),0.0)
     > + max(-scfaci(mgs),0.0)
     > + max(-scsaci(mgs),0.0)
     > + max(-scglaci(mgs),0.0)
     > + max(-scgmaci(mgs),0.0)
     > + max(-scghaci(mgs),0.0)
      psccid(mgs) =  psccid(mgs)
     > + min(-schlaci(mgs),0.0)
     > + min(-schaci(mgs),0.0)
     > + min(-scfaci(mgs),0.0)
     > + min(-scsaci(mgs),0.0)
     > + min(-scglaci(mgs),0.0)
     > + min(-scgmaci(mgs),0.0)
     > + min(-scghaci(mgs),0.0)

      qcitmp = Max(0.0,qx(mgs,li) +
     >   dtp*(pqcii(mgs)+pqcid(mgs)) )

      IF ( qcitmp .gt. qxmin(li) ) then

      IF ( ipconc .le. 0 ) THEN
      
      ccitmp = cnina(mgs)
      IF ( cimn .gt. 1.0 ) THEN
        cx(mgs,li) = Max(cimn,cx(mgs,li))
      ENDIF
      IF ( cimx .gt. 1.0 ) THEN
        cx(mgs,li) = Min(cimx,cx(mgs,li))
      ENDIF
      IF ( itype1 .ge. 1 .or. itype2 .ge. 1 ) THEN 
        ccitmp = Max(ccitmp,qcitmp*rho0(mgs)/cimasx)
        ccitmp = Min(ccitmp,qcitmp*rho0(mgs)/cimasn)
      ENDIF
      ccitmp = max(1.0,ccitmp)

      ELSE
 
       ccitmp = Max(0.001, cx(mgs,li) + dtp*(pccii(mgs)+pccid(mgs)) )
      
      ENDIF ! ( ipconc .le. 0 )
      
      
      cpci = (scx(mgs,li)
     :   +dtp*( psccii(mgs)+psccid(mgs)
     :         +psccimi(mgs)+psccimd(mgs)))/ccitmp
      
          IF ( Abs(cpci) .gt. scippmx .and.
     :       ( psccii(mgs)+psccid(mgs) ) .ne. 0.0 .and.
     :    Sign(1.0,scx(mgs,li)) .eq.
     :        Sign(1.0,( psccii(mgs)+psccid(mgs) )) ) THEN
          
           cpci0 = Sign(scippmx,cpci)
          
           scfac = (cpci0*ccitmp - scx(mgs,li) - 
     :       dtp*(psccimi(mgs)+psccimd(mgs)))/
     :              ( dtp*(psccii(mgs)+psccid(mgs)) )
c           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: sccifac = ',scfac,
     : dtp*(psccimi(mgs)+psccimd(mgs)),cpci,scx(mgs,li),ccitmp
             write(0,*) 'OUCH: sccifac = ',scfac,
     : dtp*(psccimi(mgs)+psccimd(mgs)),cpci,scx(mgs,li),ccitmp,
     : kgs(mgs),temg(mgs)
           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)

c        IF ( scfac .lt. 1.0 ) THEN
c        write(iunit,*) 'LOOK: sccifac = ',scfac,cpci,scx(mgs,li),ccitmp,
c     : kgs(mgs),temg(mgs)
c        ENDIF
         
        ELSE
          scfac = 1.0
        ENDIF
      ELSE
       scfac = 0.0
      ENDIF ! ( qcitmp .gt. qxmin(li) ) 

      IF ( scfac .lt. 1.0 ) THEN


       schlaci(mgs) = scfac*schlaci(mgs)
       schaci(mgs) = scfac*schaci(mgs)
       scfaci(mgs) = scfac*scfaci(mgs)
       scsaci(mgs) = scfac*scsaci(mgs)
       scglaci(mgs) = scfac*scglaci(mgs)
       scgmaci(mgs) = scfac*scgmaci(mgs)
       scghaci(mgs) = scfac*scghaci(mgs)

      psccii(mgs) = 0.0
      psccid(mgs) = 0.0
      psccii(mgs) =  psccii(mgs)
     > + max(-schlaci(mgs),0.0)
     > + max(-schaci(mgs),0.0)
     > + max(-scfaci(mgs),0.0)
     > + max(-scsaci(mgs),0.0)
     > + max(-scglaci(mgs),0.0)
     > + max(-scgmaci(mgs),0.0)
     > + max(-scghaci(mgs),0.0)
      psccid(mgs) =  psccid(mgs)
     > + min(-schlaci(mgs),0.0)
     > + min(-schaci(mgs),0.0)
     > + min(-scfaci(mgs),0.0)
     > + min(-scsaci(mgs),0.0)
     > + min(-scglaci(mgs),0.0)
     > + min(-scgmaci(mgs),0.0)
     > + min(-scghaci(mgs),0.0)
      
      ENDIF

      end do
      end if
c
c  cloud ice plates
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscipi(mgs) = 0.0
      pscipd(mgs) = 0.0
      pscipi(mgs) =  pscipi(mgs)
     > + max(-schlacip(mgs),0.0)
     > + max(-schacip(mgs),0.0)
     > + max(-scfacip(mgs),0.0)
     > + max(-scsacip(mgs),0.0)
     > + max(-scglacip(mgs),0.0)
     > + max(-scgmacip(mgs),0.0)
     > + max(-scghacip(mgs),0.0)
      pscipd(mgs) =  pscipd(mgs)
     > + min(-schlacip(mgs),0.0)
     > + min(-schacip(mgs),0.0)
     > + min(-scfacip(mgs),0.0)
     > + min(-scsacip(mgs),0.0)
     > + min(-scglacip(mgs),0.0)
     > + min(-scgmacip(mgs),0.0)
     > + min(-scghacip(mgs),0.0)

c
c qiptmp is the soon-to-be new qip
c
      qiptmp = Max(0.0,qx(mgs,lip) +
     >   dtp*(pqipi(mgs)+pqipd(mgs)) )
      
      IF ( qiptmp .gt. qxmin(lip) ) then
      
      IF ( ipconc .le. 0 ) THEN
      
      ciptmp = cnina(mgs)
      IF ( cimn .gt. 1.0 ) THEN
        cx(mgs,lip) = Max(cimn,cx(mgs,lip))
      ENDIF
      IF ( cimx .gt. 1.0 ) THEN
        cx(mgs,lip) = Min(cimx,cx(mgs,lip))
      ENDIF
      IF ( itype1 .ge. 1 .or. itype2 .ge. 1 ) THEN 
       ciptmp = Max(ciptmp,qiptmp*rho0(mgs)/cimasx)
       ciptmp = Min(ciptmp,qiptmp*rho0(mgs)/cimasn)
      ENDIF
       ciptmp = max(1.0,ciptmp)
      
      ELSE
 
       ciptmp = Max(0.001, cx(mgs,lip) + dtp*(pcipi(mgs)+pcipd(mgs)) )
      
      ENDIF ! ( ipconc .le. 0 )
c
c ciptmp is the number density for qiptmp
c
c cpip is the charge per particle for ciptmp and the new charge density
c cpip1 is not used
      cpip =(scx(mgs,lip)
     :   +dtp*( pscipi(mgs)+pscipd(mgs)
     :         +pscipmi(mgs)+pscipmd(mgs)))/ciptmp
      cpip1 = scx(mgs,lip)/cx(mgs,lip)
c
c now check if the magnitude of cpip exceeds the 
c     input max value scippmx
c
         IF ( Abs(cpip) .gt. scippmx .and.
     :       ( pscipi(mgs)+pscipd(mgs) ) .ne. 0.0 .and.
     :    Sign(1.0,scx(mgs,lip)) .eq.
     :        Sign(1.0,( pscipi(mgs)+pscipd(mgs) ))  ) THEN
      
           cpip0 = Sign(scippmx,cpip)
      
           scfac = (cpip0*ciptmp - scx(mgs,lip) - 
     :      dtp*(pscipmi(mgs)+pscipmd(mgs)))/
     :              ( dtp*(pscipi(mgs)+pscipd(mgs)) )
c           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: scipfac = ',scfac
             write(0,*) 'OUCH: scipfac = ',scfac,nstep,
     : dtp*(pscipmi(mgs)+pscipmd(mgs)),dtp*(pscipi(mgs)+pscipd(mgs)),
     :  cpip1,cpip,cpip0,scx(mgs,lip),ciptmp,cx(mgs,lip),
     : kgs(mgs),temg(mgs)

           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)
        
c        IF ( scfac .lt. 1.0 .and. scfac.gt.0.0 ) THEN
c         write(0,*) 'LOOK: scipfac = ',scfac,nstep,
c     : dtp*(pscipmi(mgs)+pscipmd(mgs)),dtp*(pscipi(mgs)+pscipd(mgs)),
c     :  cpip,scx(mgs,lip),ciptmp,
c     : kgs(mgs),temg(mgs)
c        ENDIF
         
        ELSE
          scfac = 1.0
        ENDIF

      ELSE
       scfac = 0.0
      ENDIF ! ( qiptmp .gt. qxmin(lc) ) 

      IF ( scfac .lt. 1.0 .and. qx(mgs,lip) .gt. qxmin(lip) ) THEN

       schlacip(mgs) = scfac*schlacip(mgs)
       schacip(mgs) = scfac*schacip(mgs)
       scfacip(mgs) = scfac*scfacip(mgs)
       scsacip(mgs) = scfac*scsacip(mgs)
       scglacip(mgs) = scfac*scglacip(mgs)
       scgmacip(mgs) = scfac*scgmacip(mgs)
       scghacip(mgs) = scfac*scghacip(mgs)

      pscipi(mgs) = 0.0
      pscipd(mgs) = 0.0
      pscipi(mgs) =  pscipi(mgs)
     > + max(-schlacip(mgs),0.0)
     > + max(-schacip(mgs),0.0)
     > + max(-scfacip(mgs),0.0)
     > + max(-scsacip(mgs),0.0)
     > + max(-scglacip(mgs),0.0)
     > + max(-scgmacip(mgs),0.0)
     > + max(-scghacip(mgs),0.0)
      pscipd(mgs) =  pscipd(mgs)
     > + min(-schlacip(mgs),0.0)
     > + min(-schacip(mgs),0.0)
     > + min(-scfacip(mgs),0.0)
     > + min(-scsacip(mgs),0.0)
     > + min(-scglacip(mgs),0.0)
     > + min(-scgmacip(mgs),0.0)
     > + min(-scghacip(mgs),0.0)
      
      ENDIF
      end do
      end if
c
c  rimed ice
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      psciri(mgs) = 0.0
      pscird(mgs) = 0.0
      psciri(mgs) =  psciri(mgs)
     > + max(-schlacir(mgs),0.0)
     > + max(-schacir(mgs),0.0)
     > + max(-scfacir(mgs),0.0)
     > + max(-scsacir(mgs),0.0)
     > + max(-scglacir(mgs),0.0)
     > + max(-scgmacir(mgs),0.0)
     > + max(-scghacir(mgs),0.0)
      pscird(mgs) =  pscird(mgs)
     > + min(-schlacir(mgs),0.0)
     > + min(-schacir(mgs),0.0)
     > + min(-scfacir(mgs),0.0)
     > + min(-scsacir(mgs),0.0)
     > + min(-scglacir(mgs),0.0)
     > + min(-scgmacir(mgs),0.0)
     > + min(-scghacir(mgs),0.0)
c      end do
c      DO mgs=1,ngscnt

      qirtmp = Max(0.0,qx(mgs,lir) +
     >   dtp*(pqiri(mgs)+pqird(mgs)) )

      IF ( qirtmp .gt. qxmin(li) ) then
      
      if ( ipconc .le. 1 ) then
      cirdiatmp =
     >  (qirtmp*rho0(mgs)/(pi*xdn(mgs,lir)*cno(lir)))**(0.25) 
      cirtmp = cno(lir)*cirdiatmp
      elseif ( ipconc .gt. 1 ) then
      cirtmp = Max(0.001, cx(mgs,lir) + dtp*(pciri(mgs)+pcird(mgs)) )
      xdia(mgs,lir,1) =
     >  (qirtmp*rho0(mgs)
     > /(pi*xdn(mgs,lir)*cirtmp))**(0.333333)
      
      end if

      
      cpir = (scx(mgs,lir)
     :   +dtp*( psciri(mgs)+pscird(mgs)
     :         +pscirmi(mgs)+pscirmd(mgs)))/cirtmp
      
         IF ( Abs(cpir) .gt. scippmx .and.
     :       ( psciri(mgs)+pscird(mgs) ) .ne. 0.0 .and.
     :    Sign(1.0,scx(mgs,lir)) .eq.
     :        Sign(1.0,( psciri(mgs)+pscird(mgs) )) ) THEN
          
           cpir0 = Sign(scippmx,cpir)
          
           scfac = (cpir0*cirtmp - scx(mgs,lir) - 
     :      dtp*(pscirmi(mgs)+pscirmd(mgs)))/
     :              ( dtp*(psciri(mgs)+pscird(mgs)) )
c           IF ( scfac .gt. 1.0 .or. scfac .lt. 0.0 ) THEN
           IF ( scfac .gt. 1.0 ) THEN
             write(iunit,*) 'OUCH: scirfac = ',scfac,cpir,cirtmp
             write(0,*) 'OUCH: scirfac = ',scfac,cpir,cirtmp
           ENDIF
          scfac = Max(scfac,0.0)
          scfac = Min(scfac,1.0)
         
        ELSE
          scfac = 1.0
        ENDIF
      ELSE
       scfac = 0.0
      ENDIF ! ( qirtmp .gt. qxmin(li) ) 
      
      IF ( scfac .lt. 1.0 ) THEN

       schlacir(mgs) = scfac*schlacir(mgs)
       schacir(mgs) = scfac*schacir(mgs)
       scfacir(mgs) = scfac*scfacir(mgs)
       scsacir(mgs) = scfac*scsacir(mgs)
       scglacir(mgs) = scfac*scglacir(mgs)
       scgmacir(mgs) = scfac*scgmacir(mgs)
       scghacir(mgs) = scfac*scghacir(mgs)

      psciri(mgs) = 0.0
      pscird(mgs) = 0.0
      psciri(mgs) =  psciri(mgs)
     > + max(-schlacir(mgs),0.0)
     > + max(-schacir(mgs),0.0)
     > + max(-scfacir(mgs),0.0)
     > + max(-scsacir(mgs),0.0)
     > + max(-scglacir(mgs),0.0)
     > + max(-scgmacir(mgs),0.0)
     > + max(-scghacir(mgs),0.0)
      pscird(mgs) =  pscird(mgs)
     > + min(-schlacir(mgs),0.0)
     > + min(-schacir(mgs),0.0)
     > + min(-scfacir(mgs),0.0)
     > + min(-scsacir(mgs),0.0)
     > + min(-scglacir(mgs),0.0)
     > + min(-scgmacir(mgs),0.0)
     > + min(-scghacir(mgs),0.0)

      ENDIF
      
      ENDDO
      end if
c
c  snow
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscswi(mgs) = 0.0
      pscswd(mgs) = 0.0
      pscswi(mgs) =  pscswi(mgs)
     > + max(scsaci(mgs),0.0)
     > + max(scsacip(mgs),0.0)
     > + max(scsacir(mgs),0.0)
     > + max(scsacw(mgs),0.0)
     > + max(scsacr(mgs),0.0)
     > + max(-schlacs(mgs),0.0)
     > + max(-schacs(mgs),0.0)
     > + max(-scfacs(mgs),0.0)
     > + max(-scglacs(mgs),0.0)
     > + max(-scgmacs(mgs),0.0)
     > + max(-scghacs(mgs),0.0)
      pscswd(mgs) =  pscswd(mgs)
     > + min(scsaci(mgs),0.0)
     > + min(scsacip(mgs),0.0)
     > + min(scsacir(mgs),0.0)
     > + min(scsacw(mgs),0.0)
     > + min(scsacr(mgs),0.0)
     > + min(-schlacs(mgs),0.0)
     > + min(-schacs(mgs),0.0)
     > + min(-scfacs(mgs),0.0)
     > + min(-scglacs(mgs),0.0)
     > + min(-scgmacs(mgs),0.0)
     > + min(-scghacs(mgs),0.0)
      end do
      end if
c
c  Frozen drops
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscfwi(mgs) = 0.0
      pscfwd(mgs) = 0.0
      pscfwi(mgs) =  pscfwi(mgs)
     > + max(scfaci(mgs),0.0)
     > + max(scfacip(mgs),0.0)
     > + max(scfacir(mgs),0.0)
     > + max(scfacs(mgs),0.0)
     > + max(scfacw(mgs),0.0)
     > + max(scfacr(mgs),0.0)
      pscfwd(mgs) =  pscfwd(mgs)
     > + min(scfaci(mgs),0.0)
     > + min(scfacip(mgs),0.0)
     > + min(scfacir(mgs),0.0)
     > + min(scfacs(mgs),0.0)
     > + min(scfacw(mgs),0.0)
     > + min(scfacr(mgs),0.0)
      end do
      end if
c
c  Hail
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschwi(mgs) = 0.0
      pschwd(mgs) = 0.0
      pschwi(mgs) =  pschwi(mgs)
     > + max(schaci(mgs),0.0)
     > + max(schacip(mgs),0.0)
     > + max(schacir(mgs),0.0)
     > + max(schacs(mgs),0.0)
     > + max(schacw(mgs),0.0)
     > + max(schacr(mgs),0.0)
      pschwd(mgs) =  pschwd(mgs)
     > + min(schaci(mgs),0.0)
     > + min(schacip(mgs),0.0)
     > + min(schacir(mgs),0.0)
     > + min(schacs(mgs),0.0)
     > + min(schacw(mgs),0.0)
     > + min(schacr(mgs),0.0)
      end do
      end if
c
c  Hail-Large
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pschli(mgs) = 0.0
      pschld(mgs) = 0.0
      pschli(mgs) =  pschli(mgs)
     > + max(schlaci(mgs),0.0)
     > + max(schlacip(mgs),0.0)
     > + max(schlacir(mgs),0.0)
     > + max(schlacs(mgs),0.0)
     > + max(schlacw(mgs),0.0)
     > + max(schlacr(mgs),0.0)
      pschwd(mgs) =  pschwd(mgs)
     > + min(schlaci(mgs),0.0)
     > + min(schlacip(mgs),0.0)
     > + min(schlacir(mgs),0.0)
     > + min(schlacs(mgs),0.0)
     > + min(schlacw(mgs),0.0)
     > + min(schlacr(mgs),0.0)
      end do
      end if
c   
c
c  Graupel - low density
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscgli(mgs) = 0.0
      pscgld(mgs) = 0.0
      pscgli(mgs) =  pscgli(mgs)
     > + max(scglaci(mgs),0.0)
     > + max(scglacip(mgs),0.0)
     > + max(scglacir(mgs),0.0)
     > + max(scglacs(mgs),0.0)
     > + max(scglacw(mgs),0.0)
     > + max(scglacr(mgs),0.0)
      pscgld(mgs) =  pscgld(mgs)
     > + min(scglaci(mgs),0.0)
     > + min(scglacip(mgs),0.0)
     > + min(scglacir(mgs),0.0)
     > + min(scglacs(mgs),0.0)
     > + min(scglacw(mgs),0.0)
     > + min(scglacr(mgs),0.0)
      end do
      end if
c
c graupel - medium densitgm
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscgmi(mgs) = 0.0
      pscgmd(mgs) = 0.0
      pscgmi(mgs) =  pscgmi(mgs)
     > + max(scgmaci(mgs),0.0)
     > + max(scgmacip(mgs),0.0)
     > + max(scgmacir(mgs),0.0)
     > + max(scgmacs(mgs),0.0)
     > + max(scgmacw(mgs),0.0)
     > + max(scgmacr(mgs),0.0)
      pscgmd(mgs) =  pscgmd(mgs)
     > + min(scgmaci(mgs),0.0)
     > + min(scgmacip(mgs),0.0)
     > + min(scgmacir(mgs),0.0)
     > + min(scgmacs(mgs),0.0)
     > + min(scgmacw(mgs),0.0)
     > + min(scgmacr(mgs),0.0)
      end do
      end if
c
c graupel - high density
c
      if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
      do mgs = 1,ngscnt
      pscghi(mgs) = 0.0
      pscghd(mgs) = 0.0
      pscghi(mgs) =  pscghi(mgs)
     > + max(scghaci(mgs),0.0)
     > + max(scghacip(mgs),0.0)
     > + max(scghacir(mgs),0.0)
     > + max(scghacs(mgs),0.0)
     > + max(scghacw(mgs),0.0)
     > + max(scghacr(mgs),0.0)
      pscghd(mgs) =  pscghd(mgs)
     > + min(scghaci(mgs),0.0)
     > + min(scghacip(mgs),0.0)
     > + min(scghacir(mgs),0.0)
     > + min(scghacs(mgs),0.0)
     > + min(scghacw(mgs),0.0)
     > + min(scghacr(mgs),0.0)
      end do
      end if

c
c
c
c  continuity of space charge
c
      if ( ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      psctot(mgs) = psccwi(mgs) + psccwd(mgs) +       psccii(mgs) + psccid(mgs) +
     >              pscipi(mgs) + pscipd(mgs) +       psciri(mgs) + pscird(mgs) +
     >              pscrwi(mgs) + pscrwd(mgs) +       pscswi(mgs) + pscswd(mgs) +
     >              pscfwi(mgs) + pscfwd(mgs) +       pscgli(mgs) + pscgld(mgs) +
     >              pscgmi(mgs) + pscgmd(mgs) +       pscghi(mgs) + pscghd(mgs) +
     >              pschwi(mgs) + pschwd(mgs) +       pschli(mgs) + pschld(mgs) +
     >
     >              psccwmi(mgs) + psccwmd(mgs) +     psccimi(mgs) + psccimd(mgs) +
     >              pscipmi(mgs) + pscipmd(mgs) +     pscirmi(mgs) + pscirmd(mgs) +
     >              pscrwmi(mgs) + pscrwmd(mgs) +     pscswmi(mgs) + pscswmd(mgs) +
     >              pscfwmi(mgs) + pscfwmd(mgs) +     pscglmi(mgs) + pscglmd(mgs) +
     >              pscgmmi(mgs) + pscgmmd(mgs) +     pscghmi(mgs) + pscghmd(mgs) +
     >              pschwmi(mgs) + pschwmd(mgs) +     pschlmi(mgs) + pschlmd(mgs)
      end do
      
c      IF ( nonigrd .eq. -10 .and. psctot(mgs) .ne. 0.0 ) THEN
c      print*, 'psctot = ',psctot(mgs),mgs, igs(mgs),jy,kgs(mgs)
c      print*,    psccwi(mgs) ,psccwd(mgs) ,              psccii(mgs) ,psccid(mgs)
c      print*,    pscipi(mgs) ,pscipd(mgs) ,              psciri(mgs) ,pscird(mgs) 
c      print*,    pscrwi(mgs) ,pscrwd(mgs) ,              pscswi(mgs) ,pscswd(mgs) 
c      print*,    pscfwi(mgs) ,pscfwd(mgs) ,              pscgli(mgs) ,pscgld(mgs) 
c      print*,    pscgmi(mgs) ,pscgmd(mgs) ,              pscghi(mgs) ,pscghd(mgs) 
c      print*,    pschwi(mgs) ,pschwd(mgs) ,              pschli(mgs) ,pschld(mgs) 
c      print*,    psccwmi(mgs) ,psccwmd(mgs) ,            psccimi(mgs) ,psccimd(mgs) 
c      print*,    pscipmi(mgs) ,pscipmd(mgs) ,            pscirmi(mgs) ,pscirmd(mgs) 
c      print*,    pscrwmi(mgs) ,pscrwmd(mgs) 
c      print*,    pscswmi(mgs) ,pscswmd(mgs) 
c      print*,    pscfwmi(mgs) ,pscfwmd(mgs) 
c      print*,    pscglmi(mgs) ,pscglmd(mgs) 
c      print*,    pscgmmi(mgs) ,pscgmmd(mgs) 
c      print*,    pscghmi(mgs) ,pscghmd(mgs) 
c      print*,    pschwmi(mgs) ,pschwmd(mgs)     
c      print*,    pschlmi(mgs) ,pschlmd(mgs)
c      
c      ENDIF
       
c cmic$ guard
c C$PAR CRITICAL SECTION    
!! c$omp  critical
      do mgs = 1,ngscnt
        psctotmx = Max( psctotmx, psctot(mgs) )
        psctotmn = Min( psctotmn, psctot(mgs) )
      end do
!! c$omp  end critical
c C$PAR END CRITICAL SECTION   
c cmic$ endguard

      
      end if

c
c      if ( ipelec .ge. 1 ) then
c      if ( nstep/60*60 .eq. nstep .or. nstep .eq. 722 .or. 
c     >  nstep .eq. 730 ) then
c     if ( nstep/20*20 .eq. nstep ) then
c      do mgs = 1,ngscnt
c      write(121,*) mgs,temg(mgs),psctot(mgs)
c      end do
c      end if
c      end if
c
c
c
c  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
c
c
      do mgs = 1,ngscnt
      qwvp(mgs) = qwvp(mgs) +     
     >   dtp*(pqwvi(mgs)+pqwvd(mgs)) 
      qx(mgs,lc) = qx(mgs,lc) +
     >   dtp*(pqcwi(mgs)+pqcwd(mgs)) 
      qx(mgs,lip) = qx(mgs,lip) +
     >   dtp*(pqipi(mgs)+pqipd(mgs)) 
      qx(mgs,lir) = qx(mgs,lir) +
     >   dtp*(pqiri(mgs)+pqird(mgs)) 
c      IF ( qx(mgs,lr) .gt. 10.0e-3 )  THEN
c       print*, 'RAIN1a: ',igs(mgs),kgs(mgs),qx(mgs,lr)
c      ENDIF
      qx(mgs,lr) = qx(mgs,lr) +
     >   dtp*(pqrwi(mgs)+pqrwd(mgs)) 
c      IF ( qx(mgs,lr) .gt. 10.0e-3 ) THEN
c        print*, 'RAIN1b: ',igs(mgs),kgs(mgs),qx(mgs,lr)
c        print*, pqrwi(mgs),pqrwd(mgs)
c       ENDIF
      qx(mgs,ls) = qx(mgs,ls) +
     >   dtp*(pqswi(mgs)+pqswd(mgs)) 
      qx(mgs,lgl) = qx(mgs,lgl) +
     >   dtp*(pqgli(mgs)+pqgld(mgs)) 
      qx(mgs,lgm) = qx(mgs,lgm) +
     >   dtp*(pqgmi(mgs)+pqgmd(mgs)) 
      qx(mgs,lgh) = qx(mgs,lgh) +
     >   dtp*(pqghi(mgs)+pqghd(mgs)) 
      qx(mgs,lh) = qx(mgs,lh) + 
     >   dtp*(pqhwi(mgs)+pqhwd(mgs)) 
      qx(mgs,lhl) = qx(mgs,lhl) + 
     >   dtp*(pqhli(mgs)+pqhld(mgs)) 
      qx(mgs,lf) = qx(mgs,lf) + 
     >   dtp*(pqfwi(mgs)+pqfwd(mgs)) 
      end do
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qx(mgs,li) = qx(mgs,li) +
     >   dtp*(pqcii(mgs)+pqcid(mgs)) 
      end do
      end if
c
c concentrations
c
      if ( ipconc .ge. 1  ) then
      do mgs = 1,ngscnt
      cx(mgs,li) = cx(mgs,li) +
     >   dtp*(pccii(mgs)+pccid(mgs)) 
      cx(mgs,lip) = cx(mgs,lip) +
     >   dtp*(pcipi(mgs)+pcipd(mgs)) 
      cx(mgs,lir) = cx(mgs,lir) +
     >   dtp*(pciri(mgs)+pcird(mgs)) 
      IF ( ipconc .ge. 2 ) THEN
      cx(mgs,lc) = cx(mgs,lc) +
     >   dtp*(pccwi(mgs)+pccwd(mgs)) 
      ENDIF
      IF ( ipconc .ge. 3 ) THEN
      cx(mgs,lr) = cx(mgs,lr) +
     >   dtp*(pcrwi(mgs)+pcrwd(mgs)) 
      ENDIF
      IF ( ipconc .ge. 4 ) THEN
      cx(mgs,ls) = cx(mgs,ls) +
     >   dtp*(pcswi(mgs)+pcswd(mgs)) 
      ENDIF
      IF ( ipconc .ge. 5 ) THEN
      cx(mgs,lgl) = cx(mgs,lgl) +
     >   dtp*(pcgli(mgs)+pcgld(mgs)) 
      cx(mgs,lgm) = cx(mgs,lgm) +
     >   dtp*(pcgmi(mgs)+pcgmd(mgs)) 
      cx(mgs,lgh) = cx(mgs,lgh) +
     >   dtp*(pcghi(mgs)+pcghd(mgs)) 
      cx(mgs,lhl) = cx(mgs,lhl) + 
     >   dtp*(pchli(mgs)+pchld(mgs)) 
      cx(mgs,lh) = cx(mgs,lh) + 
     >   dtp*(pchwi(mgs)+pchwd(mgs)) 
      cx(mgs,lf) = cx(mgs,lf) + 
     >   dtp*(pcfwi(mgs)+pcfwd(mgs)) 
      ENDIF
      end do
      end if
c
c
c
c  residence time in cloud
c
c
      if ( iptime .ge. 1 ) then
      do mgs = 1,ngscnt
      rtcw(mgs) = rtcw(mgs) + dtp*prtcw(mgs)
c      rtrw(mgs) = rtrw(mgs) + dtp*prtrw(mgs)
      if ( iptime .ge. 2 ) then
      rtci(mgs) = rtci(mgs) + dtp*prtci(mgs)
      rtip(mgs) = rtip(mgs) + dtp*prtip(mgs)
      rtir(mgs) = rtir(mgs) + dtp*prtir(mgs)
      rtsw(mgs) = rtsw(mgs) + dtp*prtsw(mgs)
      rtgl(mgs) = rtgl(mgs) + dtp*prtgl(mgs)
      rtgm(mgs) = rtgm(mgs) + dtp*prtgm(mgs)
      rtgh(mgs) = rtgh(mgs) + dtp*prtgh(mgs)
      rtfw(mgs) = rtfw(mgs) + dtp*prtfw(mgs)
c      rthw(mgs) = rthw(mgs) + dtp*prthw(mgs)
c      rthl(mgs) = rthl(mgs) + dtp*prthl(mgs)
      ENDIF
      end do
      end if
c
c
c  sum space charge sources and sinks
c
c
      if ( ipelec .ge. 2 ) then
      do mgs = 1,ngscnt
      scx(mgs,lc) = scx(mgs,lc) +
     >   dtp*(psccwi(mgs)+psccwd(mgs)+psccwmi(mgs) +psccwmd(mgs))
      scx(mgs,li) = scx(mgs,li) +
     >   dtp*(psccii(mgs)+psccid(mgs)+psccimi(mgs) +psccimd(mgs))
      scx(mgs,lip) = scx(mgs,lip) +
     >   dtp*(pscipi(mgs)+pscipd(mgs)+pscipmi(mgs) +pscipmd(mgs))
      scx(mgs,lir) = scx(mgs,lir) +
     >   dtp*(psciri(mgs)+pscird(mgs)+pscirmi(mgs) +pscirmd(mgs))
      scx(mgs,lr) = scx(mgs,lr) +
     >   dtp*(pscrwi(mgs)+pscrwd(mgs)+pscrwmi(mgs) +pscrwmd(mgs))
      scx(mgs,ls) = scx(mgs,ls) +
     >   dtp*(pscswi(mgs)+pscswd(mgs)+pscswmi(mgs) +pscswmd(mgs)) 
      scx(mgs,lf) = scx(mgs,lf) +
     >   dtp*(pscfwi(mgs)+pscfwd(mgs)+pscfwmi(mgs) +pscfwmd(mgs))
      scx(mgs,lgl) = scx(mgs,lgl) +
     >   dtp*(pscgli(mgs)+pscgld(mgs)+pscglmi(mgs) +pscglmd(mgs))
      scx(mgs,lgm) = scx(mgs,lgm) +
     >   dtp*(pscgmi(mgs)+pscgmd(mgs)+pscgmmi(mgs) +pscgmmd(mgs))
      scx(mgs,lgh) = scx(mgs,lgh) +
     >   dtp*(pscghi(mgs)+pscghd(mgs)+pscghmi(mgs) +pscghmd(mgs))
      scx(mgs,lh) = scx(mgs,lh) +
     >   dtp*(pschwi(mgs)+pschwd(mgs)+pschwmi(mgs) +pschwmd(mgs))
      scx(mgs,lhl) = scx(mgs,lhl) +
     >   dtp*(pschli(mgs)+pschld(mgs)+pschlmi(mgs) +pschlmd(mgs))

c
c  sum total charge
c
        sctot = 0.0
        DO il = lc,lhab
        sctot = sctot + scx(mgs,il)
        ENDDO
        
c     :   scx(mgs,lc) + scx(mgs,li) +
c     :      scx(mgs,lip) + scx(mgs,lir) + scx(mgs,lr) + scx(mgs,ls) +
c     :      scx(mgs,lh) + scx(mgs,lhl) + scx(mgs,lf) + scx(mgs,lgl) +
c     :      scx(mgs,lgm) + scx(mgs,lgh) 
c cmic$ guard
c C$PAR CRITICAL SECTION    
!! c$omp  critical
       sctot1p = sctot1p + Max(0.0d0,sctot)
       sctot1n = sctot1n + Min(0.0d0,sctot)
!! c$omp  end critical
c C$PAR END CRITICAL SECTION   
c cmic$ endguard


      end do
      
      end if 

c
c start saturation adjustment
c
      if (ndebug .gt. 0 ) print*,'conc 30a'
c      include 'sam.jms.satadj.sgi'
c
c
c
c  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
c
c
c
c  set up temperature and vapor arrays
c
      do mgs = 1,ngscnt
      pqs(mgs) = (380.0)/(pres(mgs))
      theta(mgs) = thetap(mgs) + theta0(mgs)
      qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      end do
c
c  melting of cloud ice
c
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptimlw(mgs) = 0.0
      end do
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      if( temg(mgs) .gt. tfr .and.
     >    qitmp(mgs) .gt. 0.0 ) then
      qx(mgs,lc) = qx(mgs,lc) + qitmp(mgs)
      scx(mgs,lc) = scx(mgs,lc) + scx(mgs,li) + scx(mgs,lip)
      thetap(mgs) = thetap(mgs) -
     >  fcc3(mgs)*qitmp(mgs)
      ptimlw(mgs) = -fcc3(mgs)*qitmp(mgs)/dtp
      cx(mgs,lc) = cx(mgs,lc) + cx(mgs,li) + cx(mgs,lip)
      qx(mgs,li) = 0.0
      qx(mgs,lip) = 0.0
      cx(mgs,li) = 0.0
      cx(mgs,lip) = 0.0
      scx(mgs,li) = 0.0
      scx(mgs,lip) = 0.0
      qitmp(mgs) = 0.0
      end if
      end do
      end if
c
c
c      do mgs = 1,ngscnt
c      qimlw(mgs) = (qcwtmp(mgs)-qx(mgs,lc))/dtp
c      end do
c
c  homogeneous freezing of cloud water
c
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptwfzi(mgs) = 0.0
      end do
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      if( temg(mgs) .lt. thnuc + 2. .and. 
     >  qx(mgs,lc) .gt. 0.0 ) then
      
      frac = Max( 0.25, Min( 1., ((thnuc + 2.) - temg(mgs) )/4.0 ) )
      qtmp = frac*qx(mgs,lc) 

      qx(mgs,li) = qx(mgs,li) + qtmp ! qx(mgs,lc)

      IF ( ipconc .ge. 2 ) THEN
        ctmp = frac*cx(mgs,lc)
!        cx(mgs,li) = cx(mgs,li) + cx(mgs,lc)
      ELSE ! (ipconc .lt. 2 )
        IF ( an(igs(mgs),jgs,kgs(mgs)-1,lc) .gt. qx(mgs,lc) ) THEN
           qtmp = an(igs(mgs),jgs,kgs(mgs)-1,lc)  

!           cx(mgs,lc) = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
           ctmp = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
        ELSE
           cx(mgs,lc) = Max(0.0,wvel(mgs))*dtp*cwccn
     :      *gz(kgs(mgs))
          cx(mgs,lc) = cwnccn(kgs(mgs))
          ctmp = frac*cwnccn(kgs(mgs))
        ENDIF
       
      ENDIF

        cx(mgs,li) = Min(ccimx, cx(mgs,li) + ctmp)

      sctmp = frac*scx(mgs,lc)
!      scx(mgs,li) = scx(mgs,li) + scx(mgs,lc)
      scx(mgs,li) = scx(mgs,li) + sctmp
!      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qx(mgs,lc)
!      ptwfzi(mgs) = fcc3(mgs)*qx(mgs,lc)/dtp
!      qx(mgs,lc) = 0.0
!      cx(mgs,lc) = 0.0
!      scx(mgs,lc) = 0.0
      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qtmp
      ptwfzi(mgs) = fcc3(mgs)*qtmp/dtp
      qx(mgs,lc) = qx(mgs,lc) - qtmp
      cx(mgs,lc) = cx(mgs,lc) - ctmp
      scx(mgs,lc) = scx(mgs,lc) - sctmp
      end if
      end do
      end if


c
c  homogeneous freezing of rain
c
c
      if ( ihabdo .lt. 1 ) then
      do mgs = 1,ngscnt
      if( temg(mgs) .lt. thnuc + 2. .and. 
     >  qx(mgs,lr) .gt. 0.0 ) then
      
      frac = Max( 0.25, Min( 1., ((thnuc + 2.) - temg(mgs) )/4.0 ) )
      qtmp = frac*qx(mgs,lr) 

      qx(mgs,li) = qx(mgs,li) + qtmp ! qx(mgs,lr)

      IF ( ipconc .ge. 3 ) THEN
        ctmp = frac*cx(mgs,lr)
!        cx(mgs,li) = cx(mgs,li) + cx(mgs,lr)
      ELSE ! (ipconc .lt. 3 )
        xdia(mgs,lr,1) =
     >  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
         cx(mgs,lr) = frac*cno(lr)*xdia(mgs,lr,1)
         
      ENDIF

       cx(mgs,li) = Min(ccimx, cx(mgs,li) + ctmp)

      sctmp = frac*scx(mgs,lr)
!      scx(mgs,li) = scx(mgs,li) + scx(mgs,lr)
      scx(mgs,li) = scx(mgs,li) + sctmp
!      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qx(mgs,lr)
!      ptwfzi(mgs) = fcc3(mgs)*qx(mgs,lr)/dtp
!      qx(mgs,lr) = 0.0
!      cx(mgs,lr) = 0.0
!      scx(mgs,lr) = 0.0
      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qtmp
      ptwfzi(mgs) = fcc3(mgs)*qtmp/dtp
      qx(mgs,lr) = qx(mgs,lr) - qtmp
      cx(mgs,lr) = cx(mgs,lr) - ctmp
      scx(mgs,lr) = scx(mgs,lr) - sctmp
      end if
      end do
      end if

c
c      do mgs = 1,ngscnt
c      qwfzi(mgs) = (qcwtmp(mgs)-qx(mgs,lc))/dtp   ! Not used?? (ERM)
c      end do
c
c  reset temporaries for cloud particles and vapor
c
      IF ( ipconc .le. 1 ) THEN
      
      do mgs = 1,ngscnt
      qx(mgs,lv) = max( 0.0, qvap(mgs) )
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qx(mgs,li) = max( 0.0, qx(mgs,li) )
      qx(mgs,lip) = max( 0.0, qx(mgs,lip) )
      qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      end do
c
c
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temgtmp = temg(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      temsav = temg(mgs)
      thsave(mgs) = thetap(mgs)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
c      IF ( ltemq .lt. 1 .or. ltemq .gt. nqsat ) THEN
c C$PAR CRITICAL SECTION
c        write(iunit,*) 'out of range ltemq!',temgtmp,temg(mgs),
c     :      thetap(mgs),theta0(mgs),pres(mgs),theta(mgs),
c     :      ltemq,igs(mgs),jy,kgs(mgs)
c        write(iunit,*) an(igs(mgs),jy,kgs(mgs),lt),
c     :   ab(igs(mgs),jy,kgs(mgs),lt),
c     :   t0(igs(mgs),jy,kgs(mgs))
c        write(iunit,*) fcc3(mgs),qx(mgs,lc),qitmp(mgs),dtp,ptem(mgs)
c        STOP
c C$PAR END CRITICAL SECTION
c      END IF
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
c      qss(kz) = qvs(kz)
c      if ( temg(kz) .lt. tfr ) then
c      if( qcw(kz) .le. qxmin(lc) .and. qci(kz) .gt. qxmin(li)) 
c     >  qss(kz) = qis(kz)
c      if( qcw(kz) .gt. qxmin(lc) .and. qci(kz) .gt. qxmin(li))
c     >   qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz)) /
c     >   (qcw(kz) + qci(kz))
c      qss(kz) = qis(kz)
c      end if
c dont get enough condensation with qcw .le./.gt. qxmin(lc)
c      if ( temg(mgs) .lt. tfr ) then
c      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
c     >  qss(mgs) = qvs(mgs)
c      if( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
c     >  qss(mgs) = qis(mgs)
c      if( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
c     >   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
c     >   (qx(mgs,lc) + qitmp(mgs))
c      else
c      qss(mgs) = qvs(mgs)
c      end if
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
     >  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
     >  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
     >   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
     >   (qx(mgs,lc) + qitmp(mgs))
      end if
      end do
c
c  iterate  adjustment
c
      do itertd = 1,2
c
      do mgs = 1,ngscnt
c
c  calculate super-saturation
c
      qitmp(mgs) = qx(mgs,li)+qx(mgs,lip)
      fcci(mgs) = 0.0
      fcip(mgs) = 0.0
      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qx(mgs,lv) - qss(mgs) )
c
c  evaporation and sublimation adjustment
c
      if( dqwv(mgs) .lt. 0. ) then           !  subsaturated
        if( qx(mgs,lc) .gt. -dqwv(mgs) ) then  ! check if qc can make up all of the deficit
          dqcw(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                 !  otherwise make all qc available for evap
          dqcw(mgs) = -qx(mgs,lc)
          dqwv(mgs) = dqwv(mgs) + qx(mgs,lc)
        end if
c
        if( qitmp(mgs) .gt. -dqwv(mgs) ) then  ! check if qi can make up all the deficit
          dqci(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                  ! otherwise make all ice available for sublimation
          dqci(mgs) = -qitmp(mgs)
          dqwv(mgs) = dqwv(mgs) + qitmp(mgs)
        end if
c
       qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs) + dqci(mgs) )  ! add to perturbation vapor
c
c This next line removed 3/19/2003 thanks to Adam Houston,
c  who found the bug in the 3-ICE code
c      qwvp(mgs) = max(qwvp(mgs), 0.0) 
      qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
        fcip(mgs) = qx(mgs,lip)/(qitmp(mgs))
      ELSE
        fcci(mgs) = 0.0
        fcip(mgs) = 0.0
      ENDIF
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
      qx(mgs,li) = qx(mgs,li) + dqci(mgs) * fcci(mgs)
      qx(mgs,lip) = qx(mgs,lip) + dqci(mgs) * fcip(mgs)
      thetap(mgs) = thetap(mgs) +
     >  cpi/pi0(mgs)*
     >  (felv(mgs)*dqcw(mgs) +fels(mgs)*dqci(mgs))

      end if  ! dqwv(mgs) .lt. 0. (end of evap/sublim)
c
c condensation/deposition
c
      IF ( dqwv(mgs) .ge. 0. ) THEN
c
        qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
        fracl(mgs) = 1.0
        fraci(mgs) = 0.0
        if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
          fracl(mgs) = max(min(1.,(temg(mgs)-233.15)/(20.)),0.0)
          fraci(mgs) = 1.0-fracl(mgs)
        end if
        if ( temg(mgs) .le. thnuc ) then
           fraci(mgs) = 1.0
           fracl(mgs) = 0.0
         end if
        fraci(mgs) = 1.0-fracl(mgs)
c
       gamss = (felv(mgs)*fracl(mgs) + fels(mgs)*fraci(mgs))
     >      / (pi0(mgs)*cp)
c
      IF ( temg(mgs) .lt. tfr ) then
        IF (qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) ) then
         dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/
     >  ((temg(mgs)-cbw)**2))
        END IF
        IF ( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
          dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv2(mgs)*qss(mgs)/
     >  ((temg(mgs)-cbi)**2))
        END IF
        IF ( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
         cdw = caw*pi0(mgs)*tfrcbw/((temg(mgs)-cbw)**2)
         cdi = cai*pi0(mgs)*tfrcbi/((temg(mgs)-cbi)**2)
         denom1 = qx(mgs,lc) + qitmp(mgs)
         denom2 = 1.0 + gamss*
     >    (qx(mgs,lc)*qvs(mgs)*cdw + qitmp(mgs)*qis(mgs)*cdi) / denom1
         dqvcnd(mgs) =  dqwv(mgs) / denom2
        END IF 

      ENDIF  !  temg(mgs) .lt. tfr
c
      if ( temg(mgs) .ge. tfr ) then
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/
     >  ((temg(mgs)-cbw)**2))
      end if
c
      delqci1=qx(mgs,li)
      delqip1=qx(mgs,lip)
c
      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
        fcip(mgs) = qx(mgs,lip)/(qitmp(mgs))
      ELSE
        fcci(mgs) = 0.0
        fcip(mgs) = 0.0
      ENDIF
c
      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs)
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)
c
      thetap(mgs) = thetap(mgs) +
     >   (felv(mgs)*dqcw(mgs) + fels(mgs)*dqci(mgs))
     > / (pi0(mgs)*cp)
      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        qx(mgs,li) = qx(mgs,li) + dqci(mgs)*fcci(mgs)
        qx(mgs,lip) = qx(mgs,lip) + dqci(mgs)*fcip(mgs)
        qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      ENDIF
c
      delqci(mgs) =  dqci(mgs)*fcci(mgs)
      delqip(mgs) =  dqci(mgs)*fcip(mgs)
c
      END IF !  dqwv(mgs) .ge. 0.
      end do
c
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li) + qx(mgs,lip)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      qvap(mgs) = Max((qwvp(mgs) + qv0(mgs)), 0.0)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qitmp(mgs) = max( 0.0, qitmp(mgs) )
      qx(mgs,lv) = max( 0.0, qvap(mgs))
c      if ( temg(mgs) .lt. tfr ) then
c      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
c     >  qss(mgs) = qvs(mgs)
cc      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
c      if( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
c     >  qss(mgs) = qis(mgs)
cc      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
c      if( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
c     >  qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
c     > (qx(mgs,lc) + qitmp(mgs))
c      else
c      qss(mgs) = qvs(mgs)
c      end if
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
     >  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
     >  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
     >   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
     >   (qx(mgs,lc) + qitmp(mgs))
      end if
c      pceds(mgs) = (thetap(mgs) - thsave(mgs))/dtp
      end do
c
c  end the saturation adjustment iteration loop
c
      end do

      ENDIF ! ( ipconc .le. 1 )
      
c
c  spread the growth owing to vapor diffusion onto the
c  ice crystal categories using the
c
c  END OF SATURATION ADJUSTMENT
c            

      if (ndebug .gt. 0 ) print*,'conc 30b'
c
c
c  end of saturation adjustment
c
c
c
c  scatter precipitation fluxes, and thetap, and hydrometeors
c
c
c
c C$PAR CRITICAL SECTION  
c cmic$  guard
!! c$omp  critical
      DO il = lc,lhab
      DO mgs = 1,ngscnt
       psctot2 = psctot2 + scx(mgs,il)
c     :   scx(mgs,lc) + scx(mgs,li) + scx(mgs,lip) + scx(mgs,lir) +
c     :   scx(mgs,lr) + scx(mgs,ls) + scx(mgs,lgl) + scx(mgs,lgm) + scx(mgs,lgh) +
c     :   scx(mgs,lf) + scx(mgs,lh) + scx(mgs,lhl)
      ENDDO
      ENDDO
!! c$omp  end critical
c cmic$ end guard
c C$PAR END CRITICAL SECTION

cCDIR$ IVDEP
c      do mgs = 1,ngscnt
cc
cc      cipvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lip,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lip,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      cirvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lir,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lir,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      cwvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lc,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lc,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      rwvt(igs(mgs),kgs(mgs),1) = vtxbar(mgs,lr,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lr,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      fwvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lf,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lf,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      swvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,ls,1)
c      vtmax = Max(vtmax,vtxbar(mgs,ls,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      glvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lgl,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lgl,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      gmvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lgm,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lgm,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      ghvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lgh,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lgh,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      hwvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lh,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lh,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
cc      hlvt(igs(mgs),kgs(mgs)) = vtxbar(mgs,lhl,1)
c      vtmax = Max(vtmax,vtxbar(mgs,lhl,1)/
c     :       dz*gt(igs(mgs),jgs,kgs(mgs),imapz))
cc
c      end do
c
c
c
cCDIR$ IVDEP
c      do mgs = 1,ngscnt
cc
c      if ( kgs(mgs) .eq. 1 ) then
cc
c      ciptot2d(igs(mgs),jgs) = ciptot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lip)*vtxbar(mgs,lip,1)/xdn(mgs,lip)
c      ciprt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lip)*vtxbar(mgs,lip,1)/xdn(mgs,lip)
cc
c      cirtot2d(igs(mgs),jgs) = cirtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lir)*vtxbar(mgs,lir,1)/xdn(mgs,lir)
c      cirrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lir)*vtxbar(mgs,lir,1)/xdn(mgs,lir)
cc
c      cwtot2d(igs(mgs),jgs) = cwtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lc)*vtxbar(mgs,lc,1)/xdn(mgs,lc)
c      cwrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lc)*vtxbar(mgs,lc,1)/xdn(mgs,lc)
cc
c      rwtot2d(igs(mgs),jgs) = rwtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lr)*vtxbar(mgs,lr,1)/xdn(mgs,lr)
c      rwrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lr)*vtxbar(mgs,lr,1)/xdn(mgs,lr)
cc
c      fwtot2d(igs(mgs),jgs) = fwtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lf)*vtxbar(mgs,lf,1)/xdn(mgs,lf)
c      fwrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lf)*vtxbar(mgs,lf,1)/xdn(mgs,lf)
cc
c      swtot2d(igs(mgs),jgs) = swtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,ls)*vtxbar(mgs,ls,1)/xdn(mgs,ls)
c      swrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,ls)*vtxbar(mgs,ls,1)/xdn(mgs,ls)
cc
c      gltot2d(igs(mgs),jgs) = gltot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lgl)*vtxbar(mgs,lgl,1)/xdn(mgs,lgl)
c      glrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lgl)*vtxbar(mgs,lgl,1)/xdn(mgs,lgl)
cc
c      gmtot2d(igs(mgs),jgs) = gmtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lgm)*vtxbar(mgs,lgm,1)/xdn(mgs,lgm)
c      gmrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lgm)*vtxbar(mgs,lgm,1)/xdn(mgs,lgm)
cc
c      ghtot2d(igs(mgs),jgs) = ghtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lgh)*vtxbar(mgs,lgh,1)/xdn(mgs,lgh)
c      ghrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lgh)*vtxbar(mgs,lgh,1)/xdn(mgs,lgh)
cc
c      hwtot2d(igs(mgs),jgs) = hwtot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lh)*vtxbar(mgs,lh,1)/xdn(mgs,lh)
c      hwrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lh)*vtxbar(mgs,lh,1)/xdn(mgs,lh)
cc
c      hltot2d(igs(mgs),jgs) = hltot2d(igs(mgs),jgs)+1000.0*dtp*
c     > qx(mgs,lhl)*vtxbar(mgs,lhl,1)/xdn(mgs,lhl)
c      hlrt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,lhl)*vtxbar(mgs,lhl,1)/xdn(mgs,lhl)
cc
c      end if
cc
c      end do
c
c
c
c      if ( ihabdo .lt. 1 ) then
cCDIR$ IVDEP
c      do mgs = 1,ngscnt
cc
c      if ( kgs(mgs) .eq. 1 ) then
cc
cc      civt(igs(mgs),kgs(mgs)) = vtxbar(mgs,li,1)
c      citot2d(igs(mgs),jgs) = citot2d(igs(mgs),jgs)+1000.*dtp*
c     > qx(mgs,li)*vtxbar(mgs,li,1)/xdn(mgs,li)
c      cirt2d(igs(mgs),jgs) = 1000.*3600.*
c     > qx(mgs,li)*vtxbar(mgs,li,1)/xdn(mgs,li)
cc
c      end if
cc
c      end do
c      end if
c
c  sum precip rates and totals
c
cCDIR$ IVDEP
c      do mgs = 1,ngscnt
c    
c      if ( kgs(mgs) .eq. 1 ) then
c
c      prt2d(igs(mgs),jgs) = 
c     >cirt2d(igs(mgs),jgs)+cirrt2d(igs(mgs),jgs)+ciprt2d(igs(mgs),jgs)
c     >+cwrt2d(igs(mgs),jgs)+rwrt2d(igs(mgs),jgs)+swrt2d(igs(mgs),jgs)
c     >+glrt2d(igs(mgs),jgs)+gmrt2d(igs(mgs),jgs)+ghrt2d(igs(mgs),jgs)
c     >+fwrt2d(igs(mgs),jgs)+hwrt2d(igs(mgs),jgs)
c     >+hltot2d(igs(mgs),jgs)
c
c      ptot2d(igs(mgs),jgs) = 
c     > citot2d(igs(mgs),jgs)+cirtot2d(igs(mgs),jgs)
c     >+ciptot2d(igs(mgs),jgs)+cwtot2d(igs(mgs),jgs)
c     >+rwtot2d(igs(mgs),jgs)+swtot2d(igs(mgs),jgs)
c     >+gltot2d(igs(mgs),jgs)+gmtot2d(igs(mgs),jgs)
c     >+ghtot2d(igs(mgs),jgs)+fwtot2d(igs(mgs),jgs)
c     >+hwtot2d(igs(mgs),jgs)
c     >+hltot2d(igs(mgs),jgs)
cc
c      prti2d(igs(mgs),jgs) = 
c     >cirt2d(igs(mgs),jgs)+cirrt2d(igs(mgs),jgs)+ciprt2d(igs(mgs),jgs)
c     >+glrt2d(igs(mgs),jgs)+gmrt2d(igs(mgs),jgs)+ghrt2d(igs(mgs),jgs)
c     >+fwrt2d(igs(mgs),jgs)+hwrt2d(igs(mgs),jgs)+swrt2d(igs(mgs),jgs)
c     >+hltot2d(igs(mgs),jgs)
cc
c      ptoti2d(igs(mgs),jgs) = 
c     > citot2d(igs(mgs),jgs)+cirtot2d(igs(mgs),jgs)
c     >+ciptot2d(igs(mgs),jgs)+gltot2d(igs(mgs),jgs)
c     >+gmtot2d(igs(mgs),jgs)+ghtot2d(igs(mgs),jgs)
c     >+fwtot2d(igs(mgs),jgs)+hwtot2d(igs(mgs),jgs)
c     >+swtot2d(igs(mgs),jgs)
c     >+hltot2d(igs(mgs),jgs)
cc
c      prtl2d(igs(mgs),jgs) = 
c     >cwrt2d(igs(mgs),jgs)+rwrt2d(igs(mgs),jgs)
cc
c      ptotl2d(igs(mgs),jgs) = 
c     >cwtot2d(igs(mgs),jgs)+rwtot2d(igs(mgs),jgs)
c
c      end if
c
c      end do
c
c
c  Plot Variables
c
c taken out for now -- doesn't do anything anyway...

c
c
c
cCDIR$ IVDEP
      do mgs = 1,ngscnt
      t0(igs(mgs),jy,kgs(mgs)) =  temg(mgs)
      end do
c
c Load the save arrays
c

! Choose process for qproc max/min/avg calculation in microphysics driver

      DO mgs = 1,ngscnt
       tq1(igs(mgs),jy,kgs(mgs)) = qwctfz(mgs)
      ENDDO

      
      CALL ldtemps (nx,ny,nz,nor,ngs,ngscnt,nstep,nsave,
     :  igs,jy,kgs,
     :  scsaci,scsacip,scsacir,scglacs,
     :  scgmacs,scghacs,
     :  scfacs,schacs,schlacs, scglaci, 
     :  scglacip, scglacir, scgmaci,
     :  scgmacip, scgmacir, scghacir,scghacip,
     :  scghaci, scfacir,
     :  scfacip, scfaci, schacir, schacip, 
     :  schaci, schlacir,
     :  schlacip, schlaci, scsacw, scglacw, 
     :  scgmacw, scghacw,
     :  scfacw, schacw, schlacw, dezcomp, 
     :  vtxbar,
c     :  vthbar, vthlbar,vtglbar,vtgmbar,vtghbar,vtfbar,
     :  sxsaci,sxsacip,sxsacir,sxglacs,sxglaci,sxglacip,sxglacir,
     :  sxgmacs,sxgmaci,sxgmacip,sxgmacir,sxghacs,sxghacir,sxghacip,
     :  sxghaci,sxfacs,sxfaci,sxfacir,sxfacip,sxhacs,sxhaci,sxhacir,
     :  sxhacip,sxhlacs,sxhlaci,sxhlacir,sxhlacip,
     :  t1,t2,t3,t4,t5,t6,t8,t9)
       
c CDIR$ IVDEP
      do mgs = 1,ngscnt
c
      an(igs(mgs),jy,kgs(mgs),lt) = 
     >  ab(kgs(mgs),lt) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) = 
     >  ab(kgs(mgs),lv) + qwvp(mgs) 

c
      
      DO il = lc,lhl
        IF ( ido(il) .eq. 1 ) THEN
         an(igs(mgs),jy,kgs(mgs),il) = qx(mgs,il) +
     >     min( an(igs(mgs),jy,kgs(mgs),il), 0.0 )
         qx(mgs,il) = an(igs(mgs),jy,kgs(mgs),il)
        ENDIF
      ENDDO

c
      end do
c

      if ( ipconc .ge. 1 ) then
cCDIR$ IVDEP
      DO il = lc,lhab

        IF ( ipconc .ge. ipc(il) ) THEN

         IF (  ipconc .ge. 4 .and. ipc(il) .ge. 4 ) THEN
           DO mgs = 1,ngscnt
            IF ( qx(mgs,il) .le. 0.0 ) THEN
              cx(mgs,il) = 0.0
            ELSE
              xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-3,cx(mgs,il)))
              IF ( xv(mgs,il) .lt. xvmn(il) .or. xv(mgs,il) .gt. xvmx(il) ) THEN
                xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
                cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xv(mgs,il)*xdn(mgs,il))
               ENDIF
            ENDIF
          ENDDO
        ENDIF

          DO mgs = 1,ngscnt
            an(igs(mgs),jy,kgs(mgs),ln(il)) = Max(cx(mgs,il), 0.0)
          ENDDO
        ENDIF
      ENDDO

      IF (  ipconc .ge. 2 ) THEN
      do mgs = 1,ngscnt
        an(igs(mgs),jy,kgs(mgs),lss) = ssmax(mgs) 
        an(igs(mgs),jy,kgs(mgs),lccn) = Min(ccwmx,ccnc(mgs))
      end do
      ENDIF

      end if 
c
c
      if ( iptemp .ge. 1 ) then
cCDIR$ IVDEP
      do mgs = 1,ngscnt
      an(igs(mgs),jy,kgs(mgs),ltpc) = tpcw(mgs) 
      an(igs(mgs),jy,kgs(mgs),ltpi) = tpci(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpip) = tpip(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpir) = tpir(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpr) = tprw(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpf) = tpfw(mgs)
      an(igs(mgs),jy,kgs(mgs),ltps) = tpsw(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpgl) = tpgl(mgs) 
      an(igs(mgs),jy,kgs(mgs),ltpgm) = tpgm(mgs)
      an(igs(mgs),jy,kgs(mgs),ltpgh) = tpgh(mgs)
      an(igs(mgs),jy,kgs(mgs),ltph) = tphw(mgs)
      an(igs(mgs),jy,kgs(mgs),ltphl) = tphl(mgs)
      end do
      end if

c
c
      if ( iptime .ge. 1 ) then
cCDIR$ IVDEP
      do mgs = 1,ngscnt
      an(igs(mgs),jy,kgs(mgs),lrtc) = max(rtcw(mgs), 0.0)
c      an(igs(mgs),jy,kgs(mgs),lrtr) = max(rtrw(mgs), 0.0)
      IF ( iptime .ge. 2 ) THEN
      an(igs(mgs),jy,kgs(mgs),lrti) = max(rtci(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrtip) = max(rtip(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrtir) = max(rtir(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrts) = max(rtsw(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrtgl) = max(rtgl(mgs), 0.0) ! need rtmax?
      an(igs(mgs),jy,kgs(mgs),lrtgm) = max(rtgm(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrtgh) = max(rtgh(mgs), 0.0)
      an(igs(mgs),jy,kgs(mgs),lrtf) = max(rtfw(mgs), 0.0)
c      an(igs(mgs),jy,kgs(mgs),lrth) = max(rthw(mgs), 0.0)
c      an(igs(mgs),jy,kgs(mgs),lrthl) = max(rthl(mgs), 0.0)
      ENDIF
      end do
      end if
c
c
c
      if ( ipelec .ge. 2 ) then
cCDIR$ IVDEP
      DO il = lc,lhab
      do mgs = 1,ngscnt
      an(igs(mgs),jy,kgs(mgs),lsc(il)) = scx(mgs,il)

      end do
      ENDDO
      end if
c
c  set ab(sc**) equal to interated area...c
c
c      if ( ipelec .ge. 2 .and. ilight .ge. 1 .and. ion.eq.0 ) then
c
c  Dont need surface area anymore and the ab array is now only for lt,le, and lv.
c
c      end if
      
c      IF ( lsavetime .and. ngscnt .ge. 1) THEN
c
c cmic$ guard
c C$PAR CRITICAL SECTION
c       CALL basave2(ngscnt,isaveunit2,igs(1),jy,kgs(1),
c     :  qcw(1),qrw(1),qci(1),qip(1),qir(1),qsw(1),qgl(1),qgm(1),
c     :  qgh(1),qhw(1),qhl(1),qfw(1),qwv(1),
c     :  sccw(1),scci(1),scir(1),scip(1),scrw(1),scsw(1),scfw(1),
c     :  schw(1),scgh(1),scgm(1),scgl(1),schl(1) )
c C$PAR END CRITICAL SECTION
c cmic$ end guard 
c      
c      END IF
      
c
c
c
      deallocate ( qx )
      deallocate ( cx )
      deallocate ( scx )
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
c
c  end of gather scatter (for this jy slice)
c
c  precipitation fallout contributions (if itfall .eq. 3 ) then
c
c
c  Do fallout stuff now if itfall=3
c
      IF ( itfall .eq. 3 ) THEN

c
c  zero the precip flux arrays (2d)
c
      xvt(:,:,:,:) = 0.0
      
      
      do kz = 1,nz-1
      DO ix = 1,nx-1
      db1(ix,kz) = dn(ix,jy,kz)
      ENDDO
      enddo

c set up vt arrays:
c      subroutine ddfall(nx,ny,nz,nor,na,dtp,dz,jgs,ngs1,
c     :  an,db,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cinccn,cwc1)

      call ddfallnx(nx,ny,nz,nor,na,dtp,dz,jgs,ngs,
     :  xvt,
     :  an,dn,ipconc,t0,t7,cwccn,cwmasn,cwmasx,cimn,cimx,
     :  rwmasn,rwmasx,cwradn,
     :  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,ccwmx,xvmn,xvmx,cwnccn,
     :  itype1,itype2,infdo)
c     :  rwmasn,rwmasx)
      
c mixing ratio      
      DO il = lc,lhab
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,il,1,xfall)
      ENDDO
      
c time
      IF ( iptime .ge. 1 ) THEN     

      DO il = lc,lhab
      IF ( lrt(il) .gt. 1 ) THEN
      call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,1,il),
     :  an,db1,imapz,mzdist,lrt(il),1,xfall)
      ENDIF
      ENDDO

      ENDIF




       IF ( ipelec .ge. 2 ) THEN
      
c      iscfall = 1
      DO il = lc,lhab
        call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,iscfall,il),
     :  an,db0,imapz,mzdist,lsc(il),0,xfall)
      ENDDO
      
      
      ENDIF ! ipelec
      
      IF ( ipconc .ge. 1 ) THEN
c      infall = 2
      DO il = lc,lhab
        IF ( ipconc .ge. ipc(il) ) THEN
         IF ( ln(il) .ne. lnr ) THEN

         call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)

         ELSE
           IF ( itfall .eq. 1 ) THEN ! special test case to use number-weighted regardless of infall

            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,2,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)

           ELSE

            call boxfall(nx,ny,nz,nor,na,gz,z1d,dtp,dz,jgs,xvt(1,1,infall,il),
     :        an,db0,imapz,mzdist,ln(il),0,xfall)

           ENDIF
         ENDIF
        ENDIF
      ENDDO
      ENDIF

      
      ENDIF ! itfall.eq.3
c
c
c
c Maximum courant number is dt*vtmax.  If this is greater than 1 then
c need to split the fallout term into smaller steps.
c
c   
c
c
c  end of jy loop
c
c
c
 9999 continue
c
c
c
c
c
c  Ziegler nucleation 
c
      IF ( ipconc .ge. 2 .and. iba .eq. 1 ) THEN


#ifdef MPI

      IF ( number_of_processes .gt. 1 ) THEN

       CALL cld_cpu('MPICOM-MICRO')

       nampi = 2


       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,ng,nampi,  
     :      w_proc(my_rank),e_proc(my_rank),westward_tag,an(-ng+1,-ng+1,-ng+1,lt) )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,ng,nampi,  
     :      w_proc(my_rank),e_proc(my_rank),eastward_tag,an(-ng+1,-ng+1,-ng+1,lt))

       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,ng,nampi,  
     :      n_proc(my_rank),s_proc(my_rank),southward_tag,an(-ng+1,-ng+1,-ng+1,lt))

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,ng,nampi,  
     :      n_proc(my_rank),s_proc(my_rank),northward_tag,an(-ng+1,-ng+1,-ng+1,lt))


        CALL cld_cpu('MPICOM-MICRO')

       ENDIF


        ssat(:,:,:) = 0.0
      ssfilt(:,:,:) = 0.0

!mpi: need ssfilt(i+1,i-1,j+1...) for ss gradient
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
         

c         t7(ix,jy,kz) = ((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
         temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz)
c     >     *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
        
c        temp1 = 
c     :    theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)) )
        
         t0(ix,jy,kz) = temp1
         
         ltemq = Int( (temp1-163.15)/fqsat+1.5 )
      IF ( ltemq .lt. 1 .or. ltemq .gt. nqsat ) THEN
        write(iunit,*) 'out of range ltemq!',temp1,
     :      an(ix,jy,kz,lt),t77(ix,jy,kz),t00(ix,jy,kz),
c     :       thetap(mgs),theta0(mgs),pres(mgs),theta(mgs),
     :      ltemq,ix,jy,kz
c        write(iunit,*) an(igs(mgs),jy,kgs(mgs),lt),
c     :   ab(igs(mgs),jy,kgs(mgs),lt),
c     :   t0(igs(mgs),jy,kgs(mgs))
c        write(iunit,*) fcc3(mgs),qx(mgs,lc),qitmp(mgs),dtp,ptem(mgs)
        STOP
c C$PAR END CRITICAL SECTION
      END IF

          c1 = t00(ix,jy,kz)*tabqvs(ltemq)

          ssat(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
                    
          ssfilt(ix,jy,kz) = ssat(ix,jy,kz)

         IF ( ny .eq. 2 .and. ssfilt(ix,jy,kz) .gt. 1.0
     :                                .and. ncdebug .ge. 1 ) THEN
          write(iunit,*) 'ssfilt1 = ',ssfilt(ix,jy,kz),
     :     an(ix,jy,kz,lv)*1.e3,c1*1.e3,an(ix,jy,kz,lt),ix,kz,
     :     an(ix,jy,kz,lnc)*1.e-6,an(ix,jy,kz,lc)*1.e3
c
c         temp1 = an(ix,jy,kz,lt)
c     >     *((pn(ix,jy,kz)+pb(ix,jy,kz))/poo)**cap
c        
c         ltemq = Int( (temp1-163.15)/fqsat+1.5 )
c
c          c1 = (380.0/(pn(ix,jy,kz)+pb(ix,jy,kz)))*tabqvs(ltemq)
c          
c          write(iunit,*) 'Alt ssat = ',100.*(an(ix,jy,kz,lv)/c1 - 1.0)


         ENDIF
      
        ENDDO
       ENDDO
      ENDDO
      
      IF ( ifilt .eq. 1 ) THEN
      do kz = 1,nz-1
      do jy = 1,ny-jstag
      do ix = 1,nx-istag
c
c filter supersaturation field
c
        IF ( ipconc .ge. 2 .and. ix .gt. 1 .and. ix .lt. nx-1 
     :       .and. kz .gt. 1 .and. kz .lt. nz-1  ) THEN
         ssfilt(ix,jy,kz)=(0.25*
     :       (ssat(ix+1,jy,kz) + ssat(ix-1,jy,kz) +
     :        ssat(ix,Min(ny-1,jy+1),kz) + ssat(ix,Max(1,jy-1),kz) +
     :        ssat(ix,jy,kz+1) + ssat(ix,jy,kz-1)) 
     :        + 0.5*ssat(ix,jy,kz))/2.
        ENDIF
c
      END DO
      END DO
      END DO
      ENDIF

c
      do 29999 jy = 1,ny-jstag
c
c  VERY IMPORTANT:  SET jgs = jy
c
      jgs = jy
      
c
c..Gather microphysics  
c
      if (ndebug .gt. 0 ) print*,'Gather stage'
      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do 2000 inumgs = 1,numgs 
      ngscnt = 0
      do 2001 kz = nzmpb,nz-kstag-1 
      do 2002 ix = nxmpb,nx-istag

      pqs(kz) = 380.0/(pn(ix,jy,kz)+pb(kz))

c      qwv(kz) = max(an(ix,jy,kz,lv),0.0)
c      qcw(kz) = max(an(ix,jy,kz,lc),0.0)
c      qci(kz) = 0.0
c      qci(kz) = max(an(ix,jy,kz,li),0.0)
c      qci(kz) = 
c     :   Max(an(ix,jy,kz,li)+an(ix,jy,kz,lip)+an(ix,jy,kz,lir),0.0)
      theta(kz) = an(ix,jy,kz,lt) 
      temg(kz) = t0(ix,jy,kz) 
c     temg(kz) = theta(kz)*((pn(ix,jy,kz)+pb(ix,jy,kz)) / poo ) ** cap
c      temg(kz) = 
c     :       theta2temp(an(ix,jy,kz,lt),(pn(ix,jy,kz)+pb(ix,jy,kz)) ) ! t0(ix,jy,kz)
      temcg(kz) = temg(kz) - tfr
      tqvcon = temg(kz)-cbw
      ltemq = (temg(kz)-163.15)/fqsat+1.5
      qvs(kz) = pqs(kz)*tabqvs(ltemq)
      qis(kz) = pqs(kz)*tabqis(ltemq)

      qss(kz) = qvs(kz)
      if ( temg(kz) .lt. tfr ) then
c      if( qcw(kz) .le. qxmin(lc) .and. qci(kz) .gt. qxmin(li)) 
c     >  qss(kz) = qis(kz)
c      if( qcw(kz) .gt. qxmin(lc) .and. qci(kz) .gt. qxmin(li))
c     >   qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz)) /
c     >   (qcw(kz) + qci(kz))
      end if
c
      if ( temg(kz) .gt. tfrh .and. 
     >   ( an(ix,jy,kz,lv)  .gt. qss(kz) .or.
     >     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or. 
     >     an(ix,jy,kz,lr)  .gt. qxmin(lr)  )      ) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
       write(iunit,*) 'add point ',ix,kz,ssat(ix,jy,kz),
     :  an(ix,jy,kz,lv)*1.e3,qss(kz)*1.e3, an(ix,jy,kz,lc)*1.e3,
     :  an(ix,jy,kz,lnc)*1.e-6
      ENDIF
      if ( ngscnt .eq. ngs ) goto 2100
      end if
 2002 continue
      nxmpb = 1
 2001 continue
c      if ( jy .eq. (ny-jstag) ) iend = 1
 2100 continue
      if ( ngscnt .eq. 0 ) go to 29998
      if ( ndebug .gt. 0 ) print*,'dbg = 5'
      
      allocate ( qx(ngscnt,lv:lhab) )
      allocate ( cx(ngscnt,lc:lhab) )
      allocate ( scx(ngscnt,lc:lhab) )
      allocate ( xv(ngscnt,lc:lhab) )
      allocate ( vtxbar(ngscnt,lc:lhab,2) )
      allocate ( xmas(ngscnt,lc:lhab) )
      allocate ( xdn(ngscnt,lc:lhab) )
      allocate ( xdia(ngscnt,lc:lhab,3) )
      
      qx(:,:) = 0.0
      cx(:,:) = 0.0

      xv(:,:) = 0.0
      xmas(:,:) = 0.0
      scx(:,:) = 0.0

c
c  define temporaries for state variables to be used in calculations
c
      DO mgs = 1,ngscnt
      qtot(mgs) = qx(mgs,lv)
      qx(mgs,lv) = an(igs(mgs),jy,kgs(mgs),lv)
       DO il = lc,lhab
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
        qtot(mgs) = qtot(mgs) + qx(mgs,il)
       ENDDO
       qcwtmp(mgs) = qx(mgs,lc)
       

       temp0(mgs) = tempbz(kgs(mgs))

      theta0(mgs) = ab(kgs(mgs),lt)
      thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
     > - ab(kgs(mgs),lt)
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      pres0(mgs) = pb(kgs(mgs))
      qv0(mgs) = ab(kgs(mgs),lv)
      qwvp(mgs) = qx(mgs,lv) - qv0(mgs) ! an(igs(mgs),jy,kgs(mgs),lv) 
!     > - ab(kgs(mgs),lv)

       theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
       presp(mgs) = pn(igs(mgs),jy,kgs(mgs))
       pres(mgs) = presp(mgs) + pres0(mgs)
       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhoinv(mgs) = 1.0/rho0(mgs)
       pi0(mgs) = piz(kgs(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
       pk(mgs) = t77(igs(mgs),jy,kgs(mgs)) ! ( pres(mgs) / poo ) ** cap
c       temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap  ! updated temperature
c       temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
       temcg(mgs) = temg(mgs) - tfr
       qss0(mgs) = (380.0)/(pres(mgs))
       pqs(mgs) = (380.0)/(pres(mgs))
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
c
        qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
        es(mgs) = 6.1078e2*tabqvs(ltemq)
        qss(mgs) = qvs(mgs)

        temgx(mgs) = min(temg(mgs),313.15)
        temgx(mgs) = max(temgx(mgs),233.15)
        felv(mgs) = 2500837.367 
     >   * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
c
        temcgx(mgs) = min(temg(mgs),273.15)
        temcgx(mgs) = max(temcgx(mgs),223.15)
        temcgx(mgs) = temcgx(mgs)-273.15
        felf(mgs) = 333690.6098 
     >  + (2030.61425)*temcgx(mgs)
     >  - (10.46708312)*temcgx(mgs)**2
c
        fels(mgs) = felv(mgs) + felf(mgs)
        fgamw(mgs) = felv(mgs)*cpi/pi0(mgs)
        fcqv1(mgs) = 4098.0258*pi0(mgs)*fgamw(mgs)

       il5(mgs) = 0
       if ( temg(mgs) .lt. tfr ) then 
        il5(mgs) = 1
       end if

      
      ENDDO

c
c load concentrations
c
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        IF ( lccn .gt. 1 ) THEN
          ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        ENDIF
       end do
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
       end do
      end if
      
      DO mgs = 1,ngscnt
      
      uvel(mgs) = (0.5)*(u(igs(mgs),jgs,kgs(mgs))
     >                  +u(igs(mgs)+1,jgs,kgs(mgs)))
      vvel(mgs) = (0.5)*(v(igs(mgs),jgs,kgs(mgs))
     >                  +v(igs(mgs),jgs+1,kgs(mgs)))
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)+1)
     >                  +w(igs(mgs),jgs,kgs(mgs)))
      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs))
     >                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))
      IF ( ssfilt(igs(mgs),jgs,kgs(mgs)) .gt. 0.0 ) THEN
        ssbar(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs))
      ELSE
        ssbar(mgs) = ssat(igs(mgs),jgs,kgs(mgs))
      ENDIF
      ssat0(mgs) = ssat(igs(mgs),jgs,kgs(mgs))
      ssf(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs))
      ssfkp1(mgs) = ssfilt(igs(mgs),jgs,kgs(mgs)+1)
      ssfkm1(mgs) = ssfilt(igs(mgs),jgs,Max(1,kgs(mgs)-1))
      ssfjp1(mgs) = ssfilt(igs(mgs),Min(ny-1,jgs+1),kgs(mgs))
      ssfjm1(mgs) = ssfilt(igs(mgs),Max(1,jgs-1),kgs(mgs))
      ssfip1(mgs) = ssfilt(Min(nx-1,igs(mgs)+1),jgs,kgs(mgs))
      ssfim1(mgs) = ssfilt(Max(1,igs(mgs)-1),jgs,kgs(mgs))
      
      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
       write(iunit,*) 'point1 ',igs(mgs),kgs(mgs),ssat0(mgs),ssbar(mgs),
     :  qx(mgs,lv)*1.e3,qss(mgs)*1.e3, qx(mgs,lc)*1.e3,
     :  cx(mgs,lc)*1.e-6
      ENDIF

      ENDDO

c  Set density
c
      if (ndebug .gt. 0 ) print*,'Set density'
c

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
      end do

c
c  cloud water variables
c
      if ( ndebug .gt. 0 )print*,'Set cloud water variables'
      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      IF ( ipconc .ge. 2 .and. cx(mgs,lc) .gt. 1.0e6 ) THEN
        xmas(mgs,lc) =
     >    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
      ELSE
       IF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 0.01 ) THEN
        xmas(mgs,lc) =
     >     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)),
     :      xdn(mgs,lc)*xvmx(lc) )
        
        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)
        
       ELSEIF ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 0.01 ) THEN
        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        
       ELSE
        xmas(mgs,lc) = cwmasn
       ENDIF
      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3

c        rb(mgs) = 0.5*xdia(mgs,lc,1)*(Sqrt(1./(1.+cnu)))**(1./3.)
c        xl2p(mgs) = 2.7e-2*xdn(mgs,lc)*cx(mgs,lc)*xv(mgs,lc)*
c     :           ((0.5e20*rb(mgs)**3*xdia(mgs,lc,1))-0.4)
c        IF ( rb(mgs) .gt. 3.51e-6 ) THEN
c          rh(mgs) = Max( 0.5*xdia(mgs,lc,1), 6.3e-4/(1.e6*(rb(mgs) - 3.50e-6)) )
c        ELSE
c          rh(mgs) = 1.e9
c        ENDIF
c        nh(mgs) = 4.2e9*xl2p(mgs)

      end do
c
c rain
c
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      
      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-3,cx(mgs,lr)))
c      parameter( xvmn(lr)=2.8866e-13, xvmx(lr)=4.1887e-9 )  ! mks
        IF ( xv(mgs,lr) .gt. xvmx(lr) ) THEN
          xv(mgs,lr) = xvmx(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        ELSEIF ( xv(mgs,lr) .lt. xvmn(lr) ) THEN
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        ENDIF

        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,1) = (xmas(mgs,lr)*cwc1)**(1./3.)
c        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

c Inverse exponential version:
c        xdia(mgs,lr,1) =
c     >  (qx(mgs,lr)*rho0(mgs)
c     > /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      ELSE
        xdia(mgs,lr,1) =
     >  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
c        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if

      end do
c
c  Ziegler nucleation 
c
c
c cloud evaporation, condensation, and nucleation
c  sqsat -> qss(mgs)
      
      DO mgs=1,ngscnt
        dcloud = 0.0
        IF ( temg(mgs) .le. tfrh ) THEN
        
        
         CYCLE
        ENDIF
        
      IF( ssat0(mgs) .GT. 0. .OR. ssf(mgs) .GT. 0. ) GO TO 620
      IF( qvap(mgs) .EQ. qss(mgs) ) GO TO 631
C
C.... EVAPORATION. QV IS LESS THAN qss(mgs).
C.... EVAPORATE CLOUD FIRST
C
      IF ( qx(mgs,lc) .LE. 0. ) GO TO 631
C.... CLOUD EVAPORATION.
c      slv = 3.15017e10 - 2.37681e7*temg(mgs)
c  slv = 1.e4*felv(mgs)
c      R1=1./(1. + A*(273.16-B)*qss(mgs)*SLV/(CP*
c     1 (T-B)**2))
c convert input 'cp' to cgs
      R1=1./(1. + caw*(273.16 - cbw)*qss(mgs)*felv(mgs)/
     :            (cp*(temg(mgs) - cbw)**2))
      QEVAP= Min( qx(mgs,lc), R1*(qss(mgs)-qvap(mgs)) )
            
      IF ( qx(mgs,lc) .LT. QEVAP ) THEN ! GO TO 63
        qwvp(mgs) = qwvp(mgs) + qx(mgs,lc)
        thetap(mgs) = thetap(mgs) - felv(mgs)*qx(mgs,lc)/(cp*pi0(mgs))
        qx(mgs,lc) = 0.
        cx(mgs,lc) = 0.
      ELSE
        qwvp(mgs) = qwvp(mgs) + QEVAP
        qx(mgs,lc) = qx(mgs,lc) - QEVAP
        IF ( qx(mgs,lc) .le. 0. ) cx(mgs,lc) = 0.
        thetap(mgs) = thetap(mgs) - felv(mgs)*QEVAP/(CP*pi0(mgs))
      ENDIF

      GO TO 631

C.... CLOUD CONDENSATION

  620 CONTINUE

        IF ( qx(mgs,lc) .GT. qxmin(lc) .and. cx(mgs,lc) .ge. 1. ) THEN



!       ac1 =  xdn(mgs,lc)*elv(kgs(mgs))**2*epsi/
!     :        (tka(kgs(mgs))*rw*temg(mgs)**2)
! took out xdn factor because it cancels later...
       ac1 =  elv(kgs(mgs))**2*epsi/(tka(kgs(mgs))*rw*temg(mgs)**2)
       

!       bc = xdn(mgs,lc)*rw*temg(mgs)/
!     :       (epsi*wvdf(kgs(mgs))*es(mgs))
! took out xdn factor because it cancels later...
       bc =   rw*temg(mgs)/(epsi*wvdf(kgs(mgs))*es(mgs))

c       bs = rho0(mgs)*((rd*temg(mgs)/(epsi*es(mgs)))+
c     :             (epsi*elv(kgs(mgs))**2/(pres(mgs)*temg(mgs)*cp)))

c       taus = Min(dtp, xdn(mgs,lc)*rho0(mgs)*(ac1+bc)/
c     :        (4*pi*0.89298*BS*0.5*xdia(mgs,lc,1)*cx(mgs,lc)*xdn(mgs,lc)))

c      
      IF ( ssf(mgs) .gt. 0.0 .or. ssat0(mgs) .gt. 0.0 ) THEN
       IF ( ny .le. 2 ) THEN
c        print*, 'undershoot: ',ssf(mgs),
c     :   ( (qx(mgs,lv) - dcloud)/c1 - 1.0)*100.
       ENDIF


       
       IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN

         IF ( xdia(mgs,lc,1) .le. 0.0 ) THEN
          xmas(mgs,lc) = cwmasn
          xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3 
         ENDIF
        d1 = (1./(ac1 + bc))*4.0*pi*ventc
     :        *0.5*xdia(mgs,lc,1)*cx(mgs,lc)*rhoinv(mgs)
       
       ELSE
         d1 = 0.0
       ENDIF

       IF ( rcond .eq. 2 .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
       rwvent(mgs) = ventr*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)

       d1r = (1./(ac1 + bc))*4.0*pi*rwvent(mgs)
     :        *0.5*xdia(mgs,lr,1)*cx(mgs,lr)*rhoinv(mgs)
       ELSE
       d1r = 0.0
       ENDIF
       
       IF ( icond .eq. 2 .and. temg(mgs) .lt. tfr ) THEN

        fvds(mgs) = 
     >    (4.0*pi*rhoinv(mgs))*
     >    (1.0/(fai(mgs)+fbi(mgs)))
        
        IF ( qx(mgs,li) .gt. qxmin(li) ) THEN
         d1i = fvds(mgs)*cx(mgs,li)*civent(mgs)*cicap(mgs)
        ELSE
         d1i = 0.0
        ENDIF

        IF ( qx(mgs,ls) .gt. qxmin(ls) ) THEN
         d1s = fvds(mgs)*cx(mgs,ls)*swvent(mgs)*swcap(mgs)
        ELSE
         d1s = 0.0
        ENDIF

        IF ( d1i .gt. 0.0 .or. d1s .gt. 0.0 ) THEN
         e1i = fels(mgs)/(CP*pi0(mgs))
        ELSE
         e1i = 0.0
        ENDIF

       ELSE
         e1i = 0.0
         d1i = 0.0
         d1s = 0.0
       ENDIF
       
       e1  = felv(mgs)/(CP*pi0(mgs))
       f1 = pk(mgs) ! (pres(mgs)/poo)**cap

c
c  fifth trial to see what happens: 
c
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq1 = ltemq 
       temp1 = temg(mgs)
       p380 = 380.0/pres(mgs)
       
c       taus = Max( 0.05*dtp, Min(taus, 0.25*dtp ) )
c       nc = NInt(dtp/Min(1.0,0.5*taus))
c       dtcon = dtp/float(nc)
       ss1 = qx(mgs,lv)/qvs(mgs)
       ss2 = ss1
       ssi1 = qx(mgs,lv)/qis(mgs)
       ssi2 = ssi1
       temp2 = temp1
       qv1 = qx(mgs,lv)
       qvs1 = qvs(mgs)
       dt1 = 0.0

          
c          dtcon = Max(dtcon,0.2)
c          nc = Nint(dtp/dtcon)

       ltemq1 = ltemq 
c want to start out with a small time step to handle the steep slope
c and fast changes, then can switch to a larger step (dtcon2) for the
c rest of the big time step.
c base the initial time step (dtcon1) on the slope (delta)
       IF ( Abs(ss1 - 1.0) .gt. 1.e-5 ) THEN
         delta = 0.5*(qv1-qvs1)/(d1*(ss1 - 1.0))
       ELSE
         delta = 0.1*dtp
       ENDIF
c delta is the extrapolated time to get halfway from qv1 to qvs1
c want at least 5 time steps to the halfway point, so multiply by 0.2
c for the initial time step
       dtcon1 = Min(0.05,0.2*delta)
       nc = Max(5,2*NInt( (dtp-4.0*dtcon1)/delta))
       dtcon2 = (dtp-4.0*dtcon1)/nc
c       print*,'delt,taus = ',delta,taus,
c     :   delta/taus,dtcon1,dtcon2
c          print*, 'RK2c start',nc,dtcon,qv1*1.e3,qvs1*1.e3,
c     :       ss1*100.,temp1
       
       n = 1
       dt1 = 0.0
       nc = 0
       dqc = 0.0
       dqr = 0.0
       dqi = 0.0
       dqs = 0.0
       
       RK2c: DO WHILE ( dt1 .lt. dtp ) 
          nc = 0
          IF ( n .le. 4 ) THEN
            dtcon = dtcon1
          ELSE
            dtcon = dtcon2
          ENDIF
 609       dqv  = -(ss1 - 1.)*d1*dtcon
           dqvr = -(ss1 - 1.)*d1r*dtcon
           IF ( icond .eq. 2 .and. temg(mgs) .lt. tfr ) THEN

           ELSE
            dtemp = -0.5*e1*f1*(dqv + dqvr)
c          print*,'RK2c dqv1 = ',dqv
c calculate midpoint values:
           ltemq1m = ltemq1 + Nint(dtemp*fqsat + 0.5)
           IF ( ltemq1m .lt. 1 .or. ltemq1m .gt. nqsat ) THEN
             write(0,*) 'STOP in icezvd_dr line 3790 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,f1,dqv,dqvr = ', dtemp,e1,f1,dqv,dqvr
             write(0,*) ' d1,d1r,dtcon,ss1 = ',d1,d1r,dtcon,ss1
             write(0,*) ' dqc, dqr = ',dqc,dqr
             write(0,*) ' qv,qc,qr = ',qx(mgs,lv)*1000.,qx(mgs,lc)*1000.,qx(mgs,lr)*1000.
             write(0,*) ' i, j, k = ',igs(mgs),jy,kgs(mgs)
             write(0,*) ' dtcon1,dtcon2,delta = ',dtcon1,dtcon2,delta
             write(0,*) ' nc,dtp = ',nc,dtp
             write(0,*) ' rwvent,xdia,crw = ', rwvent(mgs),xdia(mgs,lr,1),cx(mgs,lr)
             write(0,*) ' xvr,xmasr,xdnr,cwc1 = ',xv(mgs,lr),xmas(mgs,lr),xdn(mgs,lr),cwc1
           ENDIF
            dqvs = dtemp*p380*dtabqvs(ltemq1m)
            qv1m = qv1 + dqv + dqvr
!          qv1mr = qv1r + dqvr

            qvs1m = qvs1 + dqvs
            ss1m = qv1m/qvs1m
           ENDIF

    ! check for undersaturation when no ice is present, if so, then reduce time step
          IF ( ss1m .lt. 1.  .and. (dqvii + dqvis) .eq. 0.0 ) THEN
!            dtcon = Max(dtcon1,0.5*dtcon)
!            IF ( dtcon .gt. dtcon1 ) GOTO 609
            dtcon = (0.5*dtcon)
            IF ( dtcon .ge. dtcon1 ) THEN
             GOTO 609
            ELSE
 !            print*,'dtcon = ',dtcon,n,dt1,temg(mgs),qx(mgs,lc),qx(mgs,lr),qx(mgs,li),qx(mgs,ls)
 !            print*,'dtcon2 : ',dqc,dqr,dqi,dqs,mgs,igs(mgs),kgs(mgs),ssf(mgs)
             EXIT
            ENDIF
          ENDIF
c calculate full step:
          dqv  = -(ss1m - 1.)*d1*dtcon
          dqvr = -(ss1m - 1.)*d1r*dtcon

          IF ( icond .eq. 2 .and. temg(mgs) .lt. tfr ) THEN

          ELSE

c          print*,'RK2a dqv1m = ',dqv
          dtemp = -e1*f1*(dqv + dqvr)
          ltemq1 = ltemq1 + Nint(dtemp*fqsat + 0.5)
           IF ( ltemq1 .lt. 1 .or. ltemq1 .gt. nqsat ) THEN
             write(0,*) 'STOP in icezvd_dr line 3856 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,dqv,dqvr = ', dtemp,e1,dqv,dqvr
           ENDIF
          dqvs = dtemp*p380*dtabqvs(ltemq1)

          qv1 = qv1 + dqv + dqvr
          ENDIF

          dqc = dqc - dqv
          dqr = dqr - dqvr

          qvs1 = qvs1 + dqvs
          ss1 = qv1/qvs1
          temp1 = temp1 + dtemp
c          print*, 'RK2c ',n,qv1*1000.,qvs1*1000.,ss1*100.,temp1,
c     :     dt1+dtcon,dtcon
          IF ( temp2 .eq. temp1 .or. ss2 .eq. ss1 .or. 
     :           ss1 .eq. 1.00 .or. 
     :      ( n .gt. 10 .and. ss1 .lt. 1.0005 ) ) THEN
c           print*,'RK2c break'
           EXIT
          ELSE
           ss2 = ss1
           ssi2 = ssi1
           temp2 = temp1
           dt1 = dt1 + dtcon
c           IF ( n .ge. 5 ) dtcon = 1.25*dtcon
c           dtcon = Min(dtcon, dtp - dt1)
           n = n + 1
          ENDIF
       ENDDO RK2c
       
        
        dcloud = dqc ! qx(mgs,lv) - qv1
        thetap(mgs) = thetap(mgs) + e1*(DCLOUD + dqr)
        qwvp(mgs) = qwvp(mgs) - (DCLOUD + dqr)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        qx(mgs,lr) = qx(mgs,lr) + dqr

        IF ( icond .eq. 2 .and. temg(mgs) .lt. tfr ) THEN
        thetap(mgs) = thetap(mgs) + e1i*(dqi + dqs)
        qwvp(mgs) = qwvp(mgs) - (dqi + dqs)
        qx(mgs,li) = qx(mgs,li) + dqi
        qx(mgs,ls) = qx(mgs,ls) + dqs
        ENDIF

        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*f1
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        es(mgs) = 6.1078e2*tabqvs(ltemq)
        

      IF ( ny .eq. 2 .and. ncdebug .ge. 1 ) THEN
       write(iunit,*) 'point3 ',igs(mgs),kgs(mgs),dcloud*1.e3,dcloud2*1.e3
       write(iunit,*) 'cwdia,ccw,cwdn,c1 = ',xdia(mgs,lc,1),cx(mgs,lc)*1.e-6,
     :     xdn(mgs,lc),qx(mgs,lc)*1.e3
      ENDIF

c            
      
      ENDIF  ! dcloud .gt. 0.
     
c       write(iunit,*) 'cond:ix,kz,dcloud,qcw,cwdia,ccw,taus = ',
c     :   igs(mgs),kgs(mgs),dcloud,qx(mgs,lc),
c     :   xdia(mgs,lc,1),cx(mgs,lc),taus

      ELSE  ! qc .le. qxmin(lc)

        IF ( ssf(mgs) .gt. 0.0 ) THEN ! .and.  ssmax(mgs) .lt. sscb ) THEN

          IF ( iqcinit == 1 ) THEN

!         temp  = pb(k) * t(i,j,k)
         
         qvs0   = 380.*exp(17.27*(temg(mgs)-273.)/(temg(mgs)- 36.))/pk(mgs)

         dcloud = Max(0.0, (qx(mgs,lv)-qvs0) / (1.+qvs0*f5/(temg(mgs)-36.)**2) )
          
          
          ELSEIF ( iqcinit == 3 ) THEN
              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/ 
     &             (cp*(temg(mgs) - cbw)**2))
            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
          
          ELSEIF ( iqcinit == 2 ) THEN
!              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/
!     :             (cp*(temg(mgs) - cbw)**2))
!            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
         ssmx = ssmxinit


         CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,dcloud,
     :    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felv,ssmx,pk,ngscnt)


         ENDIF
        ELSE
            dcloud = 0.0
        ENDIF

c      PT=PT+SLV*DCLOUD/(CP*PIB(K))
        thetap(mgs) = thetap(mgs) + elv(kgs(mgs))*DCLOUD/(CP*pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD

        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) !( pres(mgs) / poo ) ** cap
c        temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
        es(mgs) = 6.1078e2*tabqvs(ltemq)

        END IF ! qc .gt. 0.

c        ES=EES(PIB(K)*PT)
c        SQSAT=EPSI*ES/(PB(K)*1000.-ES)

C.... CLOUD NUCLEATION
c      T=PIB(K)*PT
c      ES=1.E3*PB(K)*QV/EPSI

      IF ( wvel(mgs) .le. 0. ) GO TO 616
      IF ( cx(mgs,lc) .le. 0. )  GO TO 613                             !TWOMEY (1959) Nucleation
      IF ( kgs(mgs) .GT. 1 .and. qx(mgs,lc) .le. qxmin(lc)) GO TO 613  !TWOMEY (1959) Nucleation
      IF ( kgs(mgs) .eq. 1 .and. wvel(mgs) .gt. 0. ) GO TO 613         !TWOMEY (1959) Nucleation
C.... ATTEMPT ZIEGLER CLOUD NUCLEATION IN CLOUD INTERIOR UNLESS...
  616 IF ( ssf(mgs) .LE. SUPCB .AND. wvel(mgs) .GT. 0. ) GO TO 631 !... weakly saturated updraft
      IF ( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1 .AND. 
     :    (ssfkp1(mgs) .GE. SUPMX .OR.
     :     ssf(mgs)    .GE. SUPMX .OR.
     :     ssfkm1(mgs) .GE. SUPMX)) GO TO 631                      !... too much vapour
      IF (ssf(mgs) .LT. 1.E-10 .OR. ssf(mgs) .GE. SUPMX) GO TO 631 !... at the extremes for ss

c
c get here if ( qc > 0 and ss > supcb) or (w < 0)
c
      DSSDZ=0.
      DSSDY=0.
      DSSDX=0.
      r2dxm=1./(2.*dx)
      r2dym=1./(2.*dy)
      r2dzm=1./(2.)*gz(kgs(mgs))

      IF ( ( ssf(mgs) .gt. ssmax(mgs) .or.  irenuc .eq. 1 )
     :   .and.  ( ( lccn .lt. 1 .and. 
     :            cx(mgs,lc) .lt. cwccn*(Min(1.0,rho0(mgs)))) .or.
     :    ( lccn .gt. 1 .and. ccnc(mgs) .gt. 0. )   )
     :    ) THEN
      IF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1
     :  .and. ssf(mgs) .gt. 0.0
     :  .and. ssfkp1(mgs) .LT. SUPMX .and. ssfkp1(mgs) .ge. 0.0
     1  .AND. ssfkm1(mgs) .LT. SUPMX .AND. ssfkm1(mgs) .ge. 0.0 
     :  .AND. ssfkp1(mgs) .gt. ssfkm1(mgs) 
     :  .and. t0(igs(mgs),jgs,Min(nz-1,kgs(mgs)+3)) .gt. 233.2) THEN
          DSSDZ = (ssfkp1(mgs) - ssfkm1(mgs))*R2DZM
c
c otherwise check for cloud base condition with updraft:
c
        ELSEIF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1
c        IF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1
     :  .and. ssf(mgs) .gt. 0.0  .and. wvel(mgs) .gt. 0.0
     :  .and. ssfkp1(mgs) .gt. 0.0
     1  .AND. ssfkm1(mgs) .le. 0.0 .and. wvelkm1(mgs) .gt. 0.0
     :  .AND. ssf(mgs) .gt. ssfkm1(mgs) 
     :  .and. t0(igs(mgs),jgs,Min(nz-1,kgs(mgs)+1)) .gt. 233.2) THEN
         DSSDZ = 2.*(ssf(mgs) - ssfkm1(mgs))*R2DZM  ! 1-sided difference
        ENDIF
 
      
      IF( jgs .GT. 1 .AND. jgs .LT. NY-1 
     :    .and. ssf(mgs) .gt. 0.0
     :    .and. ssfjm1(mgs) .LT. SUPMX .and. ssfjm1(mgs) .ge. supcb
     1    .AND. ssfjp1(mgs) .LT. SUPMX .and. ssfjp1(mgs) .ge. supcb )
c     :    .AND. ssfjp1(mgs) .GT. ssfjm1(mgs))
     3    DSSDY = (ssfjp1(mgs) - ssfjm1(mgs))*R2DYM
      
      IF( igs(mgs) .GT. 1 .AND. igs(mgs) .LT. NX-1 
     :  .and. ssf(mgs) .gt. 0.0
     :  .and. ssfim1(mgs) .LT. SUPMX .and. ssfim1(mgs) .ge. supcb
     1  .AND. ssfip1(mgs) .LT. SUPMX .AND. ssfip1(mgs) .ge. supcb )
c     :  .AND. ssfip1(mgs) .GT. ssfim1(mgs))
     3    DSSDX = (ssfip1(mgs) - ssfim1(mgs))*R2DXM
     
       ENDIF
C
CCLZ  IF(wijk.LE.0.) CN=CCN*ssfilt(ix,jy,kz)**CCK
c note: CCN -> cwccn, DELT -> dtp
      c1 = Max(0.0, rho0(mgs)*(qx(mgs,lv) - qss(mgs))/
     :        (xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3))
      IF ( lccn .lt. 1 ) THEN
       CN(mgs) = cwccn*CCK*ssf(mgs)**CCKM*dtp*
     : Max(0.0, 
     :         ( uvel(mgs)*DSSDX + 
     :           vvel(mgs)*DSSDY + 
     :           wvel(mgs)*DSSDZ) )  
      ELSE
      CN(mgs) = 
c     :   Min(Min(c1,ccnc(mgs)), cwccn*CCK*ssf(mgs)**CCKM*dtp*
c     :   Min(ccnc(mgs), cwccn*CCK*ssf(mgs)**CCKM*dtp*
     :    ( cwccn*CCK*ssf(mgs)**CCKM*dtp*
     : Max(0.0, 
     :         ( uvel(mgs)*DSSDX + 
     :           vvel(mgs)*DSSDY + 
     :           wvel(mgs)*DSSDZ) )  )
c      IF ( cn(mgs) .gt. 0 ) ccnc(mgs) = ccnc(mgs) - cn(mgs)
      ENDIF
      
      IF ( cn(mgs) .gt. 0 ) THEN
       IF ( ccnc(mgs) .lt. 5.e7 .and. cn(mgs) .ge. 5.e7 ) THEN
          cn(mgs) = 5.e7
          ccnc(mgs) = 0.0
       ELSEIF ( cn(mgs) .gt. ccnc(mgs) ) THEN
         cn(mgs) = ccnc(mgs)
         ccnc(mgs) = 0.0
       ENDIF
      cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      ENDIF

      IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .LE. qxmin(lc)) cx(mgs,lc)=0.
      GO TO 631
C.... NUCLEATION ON CLOUD INFLOW BOUNDARY POINT

  613 CONTINUE
C.... S. TWOMEY (1959)
c Note: get here if there is no previous cloud water and w > 0.
      cn(mgs) = 0.0
      IF ( ssmax(mgs) .lt. sscb .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
       CN(mgs) =   CCNE*wvel(mgs)**0.3465
        IF ( ny .le. 2 .and. cn(mgs) .gt. 0.0 
     :                    .and. ncdebug .ge. 1 ) THEN
          print*, 'CN: ',cn(mgs)*1.e-6, cx(mgs,lc)*1.e-6, qx(mgs,lc)*1.e3,
     :       wvel(mgs), dcloud*1.e3
          IF ( cn(mgs) .gt. 1.0 ) print*, 'cwrad = ',
     :       1.e6*(rho0(mgs)*qx(mgs,lc)/cn(mgs)*cwc1)**c1f3,
     :   igs(mgs),kgs(mgs),temcg(mgs), 
     :   1.e3*an(igs(mgs),jgs,kgs(mgs)-1,lc)
        ENDIF
        IF ( iccwflg .eq. 1 ) THEN
          cn(mgs) = Min(cwccn, Max(cn(mgs),
     :       rho0(mgs)*qx(mgs,lc)/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3)))
        ENDIF
      ELSE
       cn(mgs) = 0.0
      ENDIF

      IF( CN(mgs) .GT. cx(mgs,lc) ) cx(mgs,lc) = CN(mgs)
      IF( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .le. qxmin(lc) ) THEN
        cx(mgs,lc) = 0.
      ELSE
        cx(mgs,lc) = Min(cx(mgs,lc),rho0(mgs)*Max(0.0,qx(mgs,lc))/cwmasn)
      ENDIF

  631  CONTINUE

c
c Check for supersaturation greater than ssmx and adjust down
c
       ssmx = 1.1
       qv1 = qv0(mgs) + qwvp(mgs)
       qvs1 = qvs(mgs)

       IF ( qv1 .gt. (ssmx*qvs1) ) THEN
        
         ss1 = qv1/qvs1

        ssmx = 100.*(ssmx - 1.0)
        CALL QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,qvex,
     :    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felv,ssmx,pk,ngscnt)

       IF ( ny .le. 2 .and. ncdebug .ge. 1 ) THEN
        write(iunit,'(a,3(f8.3,1x),2i4)') 
     :   'Large SS1: ',(ss1-1.0)*100.,
     :   qx(mgs,lv)*1.e3,qvs(mgs)*1.e3,igs(mgs),kgs(mgs)
        write(iunit,*) 'qc,qwvp,qvex',qx(mgs,lc)*1.e3, qwvp(mgs)*1.e3,qvex*1.e3
        write(iunit,*) 'qv1,dcloud ',qv1*1.e3,dcloud*1.e3
        
       ENDIF

        IF ( qvex .gt. 0.0 ) THEN
        thetap(mgs) = thetap(mgs) + elv(kgs(mgs))*qvex/(CP*pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - qvex
        qx(mgs,lc) = qx(mgs,lc) + qvex
        cn(mgs) = Min( ccwmx, qvex/Max( cwmasn5, xmas(mgs,lc) )  )
        ccnc(mgs) = Max( 0.0, ccnc(mgs) - cn(mgs) )
        cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
        
c        write(iunit,*) 'theta = ',theta0(mgs) + thetap(mgs)

c        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap
        
        ENDIF

       IF ( ny .le. 2 .and. ncdebug .ge. 1 ) THEN
        write(iunit,*) 'theta = ',theta0(mgs) + thetap(mgs)
       ENDIF
       
       ENDIF

c
c Calculate droplet volume and check if it is within bounds.
c  Adjust if necessary
c  


      cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
      IF( cx(mgs,lc) .GT. 1.0e7 .AND. qx(mgs,lc) .GT. qxmin(lc)) THEN
c        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc)*xdn(mgs,lc))
        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc))
      ENDIF

c      svc(mgs) = Min( svc(mgs), xvmx(lc) )
c      svc(mgs) = Max( svc(mgs), xvmn(lc) )
      xmas(mgs,lc) = Min( xmas(mgs,lc), cwmasx )
      xmas(mgs,lc) = Max( xmas(mgs,lc), cwmasn )
c      IF( SVC .GT. xvmx(lc) ) SVC=xvmx(lc)
c      IF( SVC .LT. xvmn(lc) ) SVC=xvmn(lc)


      IF( cx(mgs,lc) .GT. 10.e6 .AND. qx(mgs,lc) .GT. qxmin(lc) ) GO TO 681
        ccwtmp = cx(mgs,lc)
        cwmastmp = xmas(mgs,lc)
c       svc(mgs) = xvmn(lc)
       xmas(mgs,lc) = Max(xmas(mgs,lc), cwmasn)
c      IF(qx(mgs,lc) .GT. 0. .AND. cx(mgs,lc) .EQ. 0.) 
       IF(qx(mgs,lc) .GT. qxmin(lc) .AND. cx(mgs,lc) .le. 0.) THEN
          cx(mgs,lc) = Min(0.5*cwccn,rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc))
          xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
       ENDIF
c     :        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(SVC(mgs)*xdn(mgs,lc))
c      IF(cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. 0.) 
      IF(cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. qxmin(lc)) 
     :        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
c     :        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
c      IF(qx(mgs,lc) .GT. 0. .AND. SVC(mgs) .LT. xvmn(lc)) SVC(mgs) = xvmn(lc)
c      IF(qx(mgs,lc) .GT. 0. .AND. SVC(mgs) .GT. xvmx(lc)) SVC(mgs) = xvmx(lc)
c      IF(qx(mgs,lc) .GT. 0.) cx(mgs,lc) = 
c     :     rho0(mgs)*qx(mgs,lc)/(SVC(mgs)*xdn(mgs,lc))
c      IF(qx(mgs,lc) .GT. 0. .AND. xmas(mgs,lc) .LT. cwmasn) 
      IF(qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .LT. cwmasn) 
     :          xmas(mgs,lc) = cwmasn
c      IF(qx(mgs,lc) .GT. 0. .AND. xmas(mgs,lc) .GT. cwmasx) 
      IF(qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .GT. cwmasx) 
     :    xmas(mgs,lc) = cwmasx
c      IF(qx(mgs,lc) .GT. 0.) cx(mgs,lc) = 
      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/Max(cwmasn,xmas(mgs,lc))
c        write(iunit,*) 'ccw check3: i,k,ccw,qcw,cwmas,cwmasn = ',
c     :     igs(mgs),kgs(mgs),1.e-6*cx(mgs,lc),1.e3*qx(mgs,lc),
c     :     xmas(mgs,lc),cwmasn
      ENDIF
        
c      IF ( cx(mgs,lc) .gt. 1.5e9 ) THEN
c         write(iunit,*) 
c     : 'ccw check4: i,k,ccw,ccwtmp,qcw,cwmas,cwmastmp,cwmasn = ',
c     :     igs(mgs),kgs(mgs),1.e-6*cx(mgs,lc),1.e-6*ccwtmp,1.e3*qx(mgs,lc),
c     :     xmas(mgs,lc),cwmastmp,cwmasn
c      ENDIF

 681  CONTINUE
        
      IF ( ipconc .ge. 3 ) THEN

        
        IF(cx(mgs,lr) .GT. 0. .AND. qx(mgs,lr) .GT. qxmin(lr)) 
     :       xv(mgs,lr)=rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
        IF(xv(mgs,lr) .GT. xvmx(lr)) xv(mgs,lr) = xvmx(lr)
        IF(xv(mgs,lr) .LT. xvmn(lr)) xv(mgs,lr) = xvmn(lr)
C        IF ( .not. (cx(mgs,lr) .gt. 0. .and. qx(mgs,lr) .gt. qxmin(lr))) THEN
C          xv(mgs,lr) = xvmn(lr)
C          IF(qx(mgs,lr) .gt. 0. .and. cx(mgs,lr) .eq. 0.) 
C     :       cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*xv(mgs,lr))
C          IF(cx(mgs,lr) .gt. 0. .and. qx(mgs,lr) .gt. 0.) 
C     :       xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
C          IF( qx(mgs,lr) .gt. 0. .and. xv(mgs,lr) .lt. xvmn(lr)) xv(mgs,lr)=xvmn(lr)
C          IF( qx(mgs,lr) .gt. 0. .and. xv(mgs,lr) .gt. xvmx(lr)) xv(mgs,lr)=xvmx(lr)
C          IF( qx(mgs,lr) .gt. 0.) 
C     :        cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*xv(mgs,lr))
C        ENDIF
c      IF ( qx(mgs,lr) .gt. 10.0e-3 )
c     :  print*, 'RAIN2b: ',igs(mgs),kgs(mgs),qx(mgs,lr)

      ENDIF
      


      ENDDO ! mgs



      DO mgs=1,ngscnt
      IF ( ssf(mgs) .gt. ssmax(mgs) 
     :  .and. ( idecss .eq. 0 .or. qx(mgs,lc) .gt. qxmin(lc)) ) THEN
        ssmax(mgs) = ssf(mgs)
      ENDIF
      ENDDO
c

      do mgs = 1,ngscnt

      an(igs(mgs),jy,kgs(mgs),lt) = 
     >  ab(kgs(mgs),lt) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) = 
     >  ab(kgs(mgs),lv) + qwvp(mgs) 

c
       if ( ido(lc) .eq. 1 )  then
        an(igs(mgs),jy,kgs(mgs),lc) = qx(mgs,lc) + 
     >    min( an(igs(mgs),jy,kgs(mgs),lc), 0.0 )  
c        qx(mgs,lc) = an(igs(mgs),jy,kgs(mgs),lc)
       end if
c
       if ( ido(lr) .eq. 1 )  then
        an(igs(mgs),jy,kgs(mgs),lr) = qx(mgs,lr) + 
     >    min( an(igs(mgs),jy,kgs(mgs),lr), 0.0 )  
c        qx(mgs,lr) = an(igs(mgs),jy,kgs(mgs),lr)
       end if
      
       IF (  ipconc .ge. 2 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnc) = Max(cx(mgs,lc) , 0.0)
        an(igs(mgs),jy,kgs(mgs),lss) = Max( 0.0, ssmax(mgs) )
        IF ( lccn .gt. 1 ) THEN
          an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0, Min( ccwmx, ccnc(mgs) ) )
        ENDIF
       ENDIF
       IF (  ipconc .ge. 3 ) THEN
        an(igs(mgs),jy,kgs(mgs),lnr) = Max(cx(mgs,lr) , 0.0)
       ENDIF
      end do

      deallocate ( qx )
      deallocate ( cx )
      deallocate ( scx )
      deallocate ( xv )
      deallocate ( vtxbar )
      deallocate ( xmas )
      deallocate ( xdn )
      deallocate ( xdia )

29998 continue
      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
c
      go to 2200
      else
      nzmpb = kz 
      end if
      if ( ix+1 .gt. nx-1 ) then
      nxmpb = 1
      else
      nxmpb = ix+1
      end if
 2000 continue
 2200 continue
c
c  end of gather scatter (for this jy slice)


c
c
c  end of jy loop
c
c
c
29999  continue
c
c
      ENDIF ! ( ipconc .ge. 2 .and. iba .eq. 1 )



c
c  Intermediate charge totals check
c
c  
c  Count up total net positive (chgpos) and net negative charge (chgneg)
c
      IF ( ipelec .gt. 0 ) THEN
      chgneg = 0.0d0
      chgpos = 0.0d0
! C$DOACROSS LOCAL(kz,jy,ix), REDUCTION(chgpos,chgneg)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : chgpos,chgneg)
      do kz = 1, nz-1*kd1
      do jy = 1, ny-1*jd1
      do ix = 1, nx-1*id1
      sctot = 
     > (  an(ix,jy,kz,lscw) 
     >  + an(ix,jy,kz,lscr) 
     >  + an(ix,jy,kz,lsci) 
     >  + an(ix,jy,kz,lscir) 
     >  + an(ix,jy,kz,lscip) 
     >  + an(ix,jy,kz,lscs)
     >  + an(ix,jy,kz,lscgl)
     >  + an(ix,jy,kz,lscgm)
     >  + an(ix,jy,kz,lscgh)
     >  + an(ix,jy,kz,lscf) 
     >  + an(ix,jy,kz,lschl) 
     >  + an(ix,jy,kz,lsch)
     :  + ec*(an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni)) )
      
      dv = dxx(ix)*dyy(jy)*dzz(kz)
      IF ( sctot .gt. 0 ) THEN
        chgpos = chgpos + sctot*dv
      ELSE
        chgneg = chgneg + sctot*dv
      END IF

      
      end do
      end do
      end do
      

      write(iunit,'(a,3(2x,1pe15.8))') 'post-jy pos/neg/net charge (C):',
     :      chgpos,chgneg,chgpos+chgneg

      
      ENDIF 

c
c      timetd2 = etime(tarray)
c      timetd2 = tarray(1)
c      write(iunit,*) 'tdmphy1 time inside is ',timetd2 - timetd1

      IF ( iptotal .gt. 0 ) THEN
        write(iunit,*) 'YIKES! Violated eqtot ',iptotal,' times!'
        write(iunit,*) 'ptotalmn,ptotalmx = ',ptotalmn,ptotalmx
      END IF
      write(iunit,*) 'Hit graupel-snow max ', nscmax,' times'
      write(iunit,*) 'Hit graupel-snow min ', nscmin,' times'
      write(iunit,*) 'wet growth at ', iwetg, ' points '
      write(iunit,*) 'Max/Min psctot: ',psctotmx,psctotmn
      IF ( ipelec .gt. 0 ) THEN
      write(iunit,'(a,2(2x,1pe12.5))') 
     :    'Charge total begin/end',
     :     psctot1*dx*dy*dz, psctot2*dx*dy*dz
c,(psctot1+psctot2)*dx*dy*dz
      write(iunit,'(a,3(2x,1pe12.5))') 
     : 'Charge pos/neg/net',sctot1p*dx*dy*dz,
     :   sctot1n*dx*dy*dz,(sctot1p+sctot1n)*dx*dy*dz
      ENDIF
      call flush(iunit)
      
c
c
c  set save counter:  nsvcnt
c
c      
      nsvcnt = nsvcnt + 1
c
c/////////////////////////////////////////////////////////////////////
c/////  save microphysical data at interval nsvmp                /////
c/////////////////////////////////////////////////////////////////////
c
c
c  end of nsave
c
      if ( iptime .ge. 1 ) then
c cmic$ do all autoscope
c cmic$1 private(kz,jy,ix),
c cmic$2 shared(an,t0)
! C$DOACROSS LOCAL(kz,jy,ix),
C$&        SHARE(an,t0)      
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),SHARED(an,t0)  
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
c
      if ( an(ix,jy,kz,lc) .le. qxmin(lc) ) then
        an(ix,jy,kz,lrtc) = 0.0
      end if
c
c
      IF ( iptime .ge. 2 ) THEN
      if ( t0(ix,jy,kz) .gt. 273.15 ) then
        an(ix,jy,kz,lrtip) = 0.0
      end if
      if ( an(ix,jy,kz,lip) .le. qxmin(li) ) then
        an(ix,jy,kz,lrtip) = 0.0
      end if
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then
        an(ix,jy,kz,lrtir) = 0.0
      end if
      if ( an(ix,jy,kz,lir) .le. qxmin(li) ) then
        an(ix,jy,kz,lrtir) = 0.0
      end if
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then
        an(ix,jy,kz,lrts) = 0.0
      end if
      if ( an(ix,jy,kz,ls) .le. qxmin(ls) ) then
        an(ix,jy,kz,lrts) = 0.0
      end if
c
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then
        an(ix,jy,kz,lrtf) = 0.0
      end if
      if ( an(ix,jy,kz,lf) .le. qxmin(lh) ) then
        an(ix,jy,kz,lrtf) = 0.0
      end if
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then
        an(ix,jy,kz,lrtgl) = 0.0
      end if
        if ( an(ix,jy,kz,lgl) .le. qxmin(lh) ) then
      an(ix,jy,kz,lrtgl) = 0.0
      end if
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then 
        an(ix,jy,kz,lrtgm) = 0.0
      end if
      if ( an(ix,jy,kz,lgm) .le. qxmin(lh) ) then
        an(ix,jy,kz,lrtgm) = 0.0
      end if
c
      if ( t0(ix,jy,kz) .gt. 273.15 ) then      
        an(ix,jy,kz,lrtgh) = 0.0
      end if
      if ( an(ix,jy,kz,lgh) .le. qxmin(lh) ) then
        an(ix,jy,kz,lrtgh) = 0.0
      end if
      
      ENDIF
c
      end do
      end do
      end do
      end if
c
c
c
c  set temperature arrays each time step to ambient if q<qxmin
c
c
      if ( iptemp .ge. 1 ) then
c cmic$ do all autoscope
c cmic$1 private(kz,jy,ix),
c cmic$2 shared(an)
! C$DOACROSS LOCAL(kz,jy,ix),
C$&        SHARE(an)      
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),SHARED(an)  
      do kz = 1,nz-1
      do jy = 1,ny-1
      do ix = 1,nx-1
c
      if ( an(ix,jy,kz,lc) .lt. qxmin(lc) ) then
      an(ix,jy,kz,ltpc) = an(ix,jy,kz,lt)*piz(kz)
      end if
      if ( an(ix,jy,kz,ltpc) .gt. an(ix,jy,kz,lt)*piz(kz) ) then
      an(ix,jy,kz,ltpc) = an(ix,jy,kz,lt)*piz(kz)
      end if
c
      if ( an(ix,jy,kz,li) .gt. qxmin(li) ) then
      if ( an(ix,jy,kz,ltpi) .gt. 273.15 ) then
      an(ix,jy,kz,ltpi) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,li) .lt. qxmin(li) ) then
      an(ix,jy,kz,ltpi) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lr) .lt. qxmin(lr) ) then
      an(ix,jy,kz,ltpr) = an(ix,jy,kz,lt)*piz(kz)
      end if
      if ( an(ix,jy,kz,ltpr) .gt. an(ix,jy,kz,lt)*piz(kz) ) then
      an(ix,jy,kz,ltpr) = an(ix,jy,kz,lt)*piz(kz)
      end if
c
      if ( an(ix,jy,kz,lip) .gt. qxmin(li) ) then
      if ( an(ix,jy,kz,ltpip) .gt. 273.15 ) then
      an(ix,jy,kz,ltpip) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lip) .lt. qxmin(li) ) then
      an(ix,jy,kz,ltpip) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lir) .gt. qxmin(li) ) then
      if ( an(ix,jy,kz,ltpir) .gt. 273.15 ) then
      an(ix,jy,kz,ltpir) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lir) .lt. qxmin(li) ) then
      an(ix,jy,kz,ltpir) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,ls) .gt. qxmin(ls) ) then
      if ( an(ix,jy,kz,ltps) .gt. 273.15 ) then
      an(ix,jy,kz,ltps) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,ls) .lt. qxmin(ls) ) then
      an(ix,jy,kz,ltps) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lgl) .gt. qxmin(lgl) ) then
      if ( an(ix,jy,kz,ltpgl) .gt. 273.15 ) then
      an(ix,jy,kz,ltpgl) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lgl) .lt. qxmin(lgl) ) then
      an(ix,jy,kz,ltpgl) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lgm) .gt. qxmin(lgm) ) then
      if ( an(ix,jy,kz,ltpgm) .gt. 273.15 ) then
      an(ix,jy,kz,ltpgm) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lgm) .lt. qxmin(lgm) ) then
      an(ix,jy,kz,ltpgm) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lgh) .gt. qxmin(lgh) ) then
      if ( an(ix,jy,kz,ltpgh) .gt. 273.15 ) then
      an(ix,jy,kz,ltpgh) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lgh) .lt. qxmin(lgh) ) then
      an(ix,jy,kz,ltpgh) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lf) .gt. qxmin(lf) ) then
      if ( an(ix,jy,kz,ltpf) .gt. 273.15 ) then
      an(ix,jy,kz,ltpf) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lf) .lt. qxmin(lf) ) then
      an(ix,jy,kz,ltpf) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lhl) .gt. qxmin(lh) ) then
      if ( an(ix,jy,kz,ltphl) .gt. 273.15 ) then
      an(ix,jy,kz,ltphl) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lhl) .lt. qxmin(lh) ) then
      an(ix,jy,kz,ltphl) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      if ( an(ix,jy,kz,lh) .gt. qxmin(lh) ) then
      if ( an(ix,jy,kz,ltph) .gt. 273.15 ) then
      an(ix,jy,kz,ltph) = 273.15
      end if
      end if
      if ( an(ix,jy,kz,lh) .lt. qxmin(lh) ) then
      an(ix,jy,kz,ltph) = min(273.15, an(ix,jy,kz,lt)*piz(kz))
      end if
c
      end do
      end do
      end do
      end if 
c
c
c
c
c  Set limits on riming age:
c
      IF ( iptime .ge. 1 ) THEN
      DO kz =  1, nz-kd1
        DO jy =  1, ny-jd1
          DO ix =  1, nx-id1
          an(ix,jy,kz,lrtc) = Max(an(ix,jy,kz,lrtc),0.0)
c          an(ix,jy,kz,lrtr) = Max(an(ix,jy,kz,lrtr),0.0)
          IF ( iptime .ge. 2 ) THEN
          an(ix,jy,kz,lrti) = Min(Max(an(ix,jy,kz,lrti),0.0),rimtim)
          an(ix,jy,kz,lrtip) = Min(Max(an(ix,jy,kz,lrtip),0.0),rimtim)
          an(ix,jy,kz,lrtir) = Min(Max(an(ix,jy,kz,lrtir),0.0),rimtim)
          an(ix,jy,kz,lrts) = Min(Max(an(ix,jy,kz,lrts),0.0),rimtim)
          an(ix,jy,kz,lrtf) = Min(Max(an(ix,jy,kz,lrtf),0.0),rimtim)
c          an(ix,jy,kz,lrth) = Min(Max(an(ix,jy,kz,lrth),0.0),rimtim)
c          an(ix,jy,kz,lrthl) = Min(Max(an(ix,jy,kz,lrthl),0.0),rimtim)
          an(ix,jy,kz,lrtgl) = Min(Max(an(ix,jy,kz,lrtgl),0.0),rimtim)
          an(ix,jy,kz,lrtgm) = Min(Max(an(ix,jy,kz,lrtgm),0.0),rimtim)
          an(ix,jy,kz,lrtgh) = Min(Max(an(ix,jy,kz,lrtgh),0.0),rimtim)
          ENDIF
         ENDDO
        ENDDO
      ENDDO
      ENDIF
c
c

c
c
c
      if ( nstep/iholen*iholen .eq. nstep ) then
c
c  hole filling routine: done at interval iholen 
c
c  zero values
c
      cqtotn = 0.0
      cctotn = 0.0
      citotn = 0.0
      cirtotn = 0.0
      ciptotn = 0.0
      crtotn = 0.0
      cstotn = 0.0
      cvtotn = 0.0
      cftotn = 0.0
      cgltotn = 0.0
      cgmtotn = 0.0
      cghtotn = 0.0
      chtotn = 0.0
      chltotn = 0.0
      cqtotp = 0.0
      cctotp = 0.0
      citotp = 0.0
      cirtotp = 0.0
      ciptotp = 0.0
      crtotp = 0.0
      cstotp = 0.0
      cvtotp = 0.0
      cftotp = 0.0
      cgltotp = 0.0
      cgmtotp = 0.0
      cghtotp = 0.0
      chtotp = 0.0
      chltotp = 0.0
      cqfac = 0.0
      ccfac = 0.0
      cifac = 0.0
      crfac = 0.0
      csfac = 0.0
      cvfac = 0.0
      cffac = 0.0
      cglfac = 0.0
      cgmfac = 0.0
      cghfac = 0.0
      chfac = 0.0
      chlfac = 0.0
c
c  sum to find totals (negative and positive)
c
      if ( iholef .eq. 1 ) then
! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotn)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : cqtotn)
      do 77611 kz =  1, nz-kd1
      do 77612 jy =  1, ny-jd1
      do 77613 ix =  1, nx-id1

      dv = dxx(ix)*dyy(jy)*dzz(kz)

      cqtotn = cqtotn + dn(ix,jy,kz)*dv*
     >               ( an(ix,jy,kz,lc)
     :                + an(ix,jy,kz,li)
     >                + an(ix,jy,kz,lir)
     >                + an(ix,jy,kz,lip)
     >                + an(ix,jy,kz,lr)
     >                + an(ix,jy,kz,ls)
     >                + an(ix,jy,kz,lv)
     >                + an(ix,jy,kz,lf)
     >                + an(ix,jy,kz,lh)
     >                + an(ix,jy,kz,lgl)
     >                + an(ix,jy,kz,lgm)
     >                + an(ix,jy,kz,lgh)
     >                + an(ix,jy,kz,lhl) )
77613 continue
77612 continue
77611 continue
c      if ( ihabdo .lt. 1 ) then
c      cqtotn1 = 0.0
! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotn1)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : cqtotn1)
c      do kz =  1, nz-kd1
c      do jy =  1, ny-jd1
c      do ix =  1, nx-id1
c      cqtotn1 = cqtotn1 + dn(ix,jy,kz)/gt(ix,jy,kz,imapz)*
c     :                 an(ix,jy,kz,li)
c      end do
c      end do
c      end do
c      
c      cqtotn = cqtotn + cqtotn1
c      
c      end if
c
c  set negatives to zero
c
c cmic$ do all autoscope
c cmic$1 shared (an)
c cmic$2 private (kz,jy,ix)
! C$DOACROSS  LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do 77631 kz =  1, nz-kd1
      do 77632 jy =  1, ny-jd1
      do 77633 ix =  1, nx-id1
      an(ix,jy,kz,lc) = max(an(ix,jy,kz,lc),0.0)
      an(ix,jy,kz,lr) = max(an(ix,jy,kz,lr),0.0)
      an(ix,jy,kz,li) = max(an(ix,jy,kz,li),0.0)
      an(ix,jy,kz,lir) = max(an(ix,jy,kz,lir),0.0)
      an(ix,jy,kz,lip) = max(an(ix,jy,kz,lip),0.0)
      an(ix,jy,kz,ls) = max(an(ix,jy,kz,ls),0.0)
      an(ix,jy,kz,lv) = max(an(ix,jy,kz,lv),0.0)
      an(ix,jy,kz,lf) = max(an(ix,jy,kz,lf),0.0)
      an(ix,jy,kz,lh) = max(an(ix,jy,kz,lh),0.0)
      an(ix,jy,kz,lhl) = max(an(ix,jy,kz,lhl),0.0)
      an(ix,jy,kz,lgl) = max(an(ix,jy,kz,lgl),0.0)
      an(ix,jy,kz,lgm) = max(an(ix,jy,kz,lgm),0.0)
      an(ix,jy,kz,lgh) = max(an(ix,jy,kz,lgh),0.0)
      IF ( ipconc .ge. 1 ) THEN
      an(ix,jy,kz,lni) = max(an(ix,jy,kz,lni),0.0)
      an(ix,jy,kz,lnir) = max(an(ix,jy,kz,lnir),0.0)
      an(ix,jy,kz,lnip) = max(an(ix,jy,kz,lnip),0.0)
      ENDIF
      IF ( ipconc .ge. 2 ) THEN
      an(ix,jy,kz,lnc) = max(an(ix,jy,kz,lnc),0.0)
      ENDIF
      IF ( ipconc .ge. 3 ) THEN
      an(ix,jy,kz,lnr) = max(an(ix,jy,kz,lnr),0.0)
      ENDIF
      IF ( ipconc .ge. 4 ) THEN
      an(ix,jy,kz,lns) = max(an(ix,jy,kz,lns),0.0)
      ENDIF
      IF ( ipconc .ge. 5 ) THEN
      an(ix,jy,kz,lnf) = max(an(ix,jy,kz,lnf),0.0)
      an(ix,jy,kz,lngl) = max(an(ix,jy,kz,lngl),0.0)
      an(ix,jy,kz,lngm) = max(an(ix,jy,kz,lngm),0.0)
      an(ix,jy,kz,lngh) = max(an(ix,jy,kz,lngh),0.0)
      an(ix,jy,kz,lnh) = max(an(ix,jy,kz,lnh),0.0)
      an(ix,jy,kz,lnhl) = max(an(ix,jy,kz,lnhl),0.0)
      ENDIF
77633 continue
77632 continue
77631 continue
c      if ( ihabdo .lt. 1 ) then
c  
! C$DOACROSS IF (nx*ny*nz .gt. 40*40*42), LOCAL(kz,jy,ix)
!! c$omp PARALLEL DO 
!! c$omp+ DEFAULT(SHARED), PRIVATE(kz,jy,ix)
c      do kz =  1, nz-kd1
c      do jy =  1, ny-jd1
c      do ix =  1, nx-id1
c      end do
c      end do
c      end do
c      end if
c
c  resum values
c
! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotp)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
!! c$omp+ REDUCTION (+ : cqtotp)
      do 77651 kz =  1, nz-kd1
      do 77652 jy =  1, ny-jd1
      do 77653 ix =  1, nx-id1

      dv = dxx(ix)*dyy(jy)*dzz(kz)

      cqtotp = cqtotp + dn(ix,jy,kz)*dv*
     >               ( an(ix,jy,kz,lc)
     >                + an(ix,jy,kz,lr)
     >                + an(ix,jy,kz,li)
     >                + an(ix,jy,kz,lir)
     >                + an(ix,jy,kz,lip)
     >                + an(ix,jy,kz,ls)
     >                + an(ix,jy,kz,lv)
     >                + an(ix,jy,kz,lf)
     >                + an(ix,jy,kz,lh)
     >                + an(ix,jy,kz,lhl)
     >                + an(ix,jy,kz,lgl)
     >                + an(ix,jy,kz,lgm)
     >                + an(ix,jy,kz,lgh) )
77653 continue
77652 continue
77651 continue
c      if ( ihabdo .lt. 1 ) then
c      cqtotp1 = 0.0
c! C$DOACROSS  LOCAL(kz,jy,ix), REDUCTION (cqtotp1)
c!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),
c!! c$omp+ REDUCTION (+ : cqtotp1)
c      do kz =  1, nz-kd1
c      do jy =  1, ny-jd1
c      do ix =  1, nx-id1
c      cqtotp1 = cqtotp1 + dn(ix,jy,kz)*an(ix,jy,kz,li)/
c     :                 gt(ix,jy,kz,imapz)
c      end do
c      end do
c      end do
c      
c      cqtotp = cqtotp + cqtotp1
c      
c      end if
c
c  set factors based on ratio of (totn+a)/(totp+a) where a =1.e-20
c
      cqfac =
     > (cqtotn+1.e-20)/(cqtotp+1.e-20)
      write(iunit,*) 'CQFAC = ',cqfac,cqtotn,cqtotp
      
      IF ( .not. ( cqfac .lt. 2.0 .and. cqfac .gt. 0.0 ) ) THEN
        write(iunit,*) 'Blowing up! STOP!'
        write(0,*) 'Blowing up! STOP!'
        
        DO ia = 1,na
        DO kz=1,nz-kd1
        DO jy=1,ny-jd1
        DO ix=1,nx-id1
          IF ( .not. ( an(ix,jy,kz,ia)  .lt. 1.0e10  .and. 
     :              an(ix,jy,kz,ia)  .gt. -1.0e10  )   ) THEN
           write(iunit,'(4i5,1x,1pe13.5)') 
     :      ix,jy,kz,ia,an(ix,jy,kz,ia)
          END IF
        END DO
        END DO
        END DO
        END DO
        
        STOP
      END IF
c
c  now reset all values based on factors
c
c cmic$ do all autoscope
c cmic$1 shared (an,cqfac)
c cmic$2 private (kz,jy,ix)
! C$DOACROSS  LOCAL(kz,jy,ix), share(cqfac,an)
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix)
      do 77671 kz =  1, nz-kd1
      do 77672 jy =  1, ny-jd1
      do 77673 ix =  1, nx-id1
      an(ix,jy,kz,lc) = an(ix,jy,kz,lc)*cqfac
      an(ix,jy,kz,lr) = an(ix,jy,kz,lr)*cqfac
      an(ix,jy,kz,lip) = an(ix,jy,kz,lip)*cqfac
      an(ix,jy,kz,lir) = an(ix,jy,kz,lir)*cqfac
      an(ix,jy,kz,ls) = an(ix,jy,kz,ls)*cqfac
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv)*cqfac
      an(ix,jy,kz,lf) = an(ix,jy,kz,lf)*cqfac
      an(ix,jy,kz,lh) = an(ix,jy,kz,lh)*cqfac
      an(ix,jy,kz,lhl) = an(ix,jy,kz,lhl)*cqfac
      an(ix,jy,kz,lgl) = an(ix,jy,kz,lgl)*cqfac
      an(ix,jy,kz,lgm) = an(ix,jy,kz,lgm)*cqfac
      an(ix,jy,kz,lgh) = an(ix,jy,kz,lgh)*cqfac
      an(ix,jy,kz,li) = an(ix,jy,kz,li)*cqfac
77673 continue
77672 continue
77671 continue
c
      end if
c
c  end of hole filling version 1
c
c
c  sum to find totals (negative and positive)
c
      if ( iholef .eq. 2 ) then
      do 77711 kz =  1, nz-kd1
      do 77712 jy =  1, ny-jd1
      do 77713 ix =  1, nx-id1
      cctotn = cctotn + dn(ix,jy,kz)*an(ix,jy,kz,lc)
      citotn = citotn + dn(ix,jy,kz)*an(ix,jy,kz,li)
      cirtotn = cirtotn + dn(ix,jy,kz)*an(ix,jy,kz,lir)
      ciptotn = ciptotn + dn(ix,jy,kz)*an(ix,jy,kz,lip)
      crtotn = crtotn + dn(ix,jy,kz)*an(ix,jy,kz,lr)
      cstotn = cstotn + dn(ix,jy,kz)*an(ix,jy,kz,ls)
      cvtotn = cvtotn + dn(ix,jy,kz)*an(ix,jy,kz,lv)
      cftotn = cftotn + dn(ix,jy,kz)*an(ix,jy,kz,lf)
      cgltotn = cgltotn + dn(ix,jy,kz)*an(ix,jy,kz,lgl)
      cgmtotn = cgmtotn + dn(ix,jy,kz)*an(ix,jy,kz,lgm)
      cghtotn = cghtotn + dn(ix,jy,kz)*an(ix,jy,kz,lgh)
      chtotn = chtotn + dn(ix,jy,kz)*an(ix,jy,kz,lh)
77713 continue
77712 continue
77711 continue
c
c  set negatives to zero
c
      do 77731 kz =  1, nz-kd1
      do 77732 jy =  1, ny-jd1
      do 77733 ix =  1, nx-id1
      an(ix,jy,kz,lc) = max(an(ix,jy,kz,lc),0.0)
      an(ix,jy,kz,li) = max(an(ix,jy,kz,li),0.0)
      an(ix,jy,kz,lir) = max(an(ix,jy,kz,lir),0.0)
      an(ix,jy,kz,lip) = max(an(ix,jy,kz,lip),0.0)
      an(ix,jy,kz,lr) = max(an(ix,jy,kz,lr),0.0)
      an(ix,jy,kz,ls) = max(an(ix,jy,kz,ls),0.0)
      an(ix,jy,kz,lv) = max(an(ix,jy,kz,lv),0.0)
      an(ix,jy,kz,lf) = max(an(ix,jy,kz,lf),0.0)
      an(ix,jy,kz,lgl) = max(an(ix,jy,kz,lgl),0.0)
      an(ix,jy,kz,lgm) = max(an(ix,jy,kz,lgm),0.0)
      an(ix,jy,kz,lgh) = max(an(ix,jy,kz,lgh),0.0)
      an(ix,jy,kz,lh) = max(an(ix,jy,kz,lh),0.0)
77733 continue
77732 continue
77731 continue
c77740 continue
c
c  resum values
c
      do 77751 kz =  1, nz-kd1
      do 77752 jy =  1, ny-jd1
      do 77753 ix =  1, nx-id1
      cctotp = cctotp + dn(ix,jy,kz)*an(ix,jy,kz,lc)
      citotp = citotp + dn(ix,jy,kz)*an(ix,jy,kz,li)
      cirtotp = cirtotp + dn(ix,jy,kz)*an(ix,jy,kz,lir)
      ciptotp = ciptotp + dn(ix,jy,kz)*an(ix,jy,kz,lip)
      crtotp = crtotp + dn(ix,jy,kz)*an(ix,jy,kz,lr)
      cstotp = cstotp + dn(ix,jy,kz)*an(ix,jy,kz,ls)
      cvtotp = cvtotp + dn(ix,jy,kz)*an(ix,jy,kz,lv)
      cftotp = cftotp + dn(ix,jy,kz)*an(ix,jy,kz,lf)
      cgltotp = cgltotp + dn(ix,jy,kz)*an(ix,jy,kz,lgl)
      cgmtotp = cgmtotp + dn(ix,jy,kz)*an(ix,jy,kz,lgm)
      cghtotp = cghtotp + dn(ix,jy,kz)*an(ix,jy,kz,lgh)
      chtotp = chtotp + dn(ix,jy,kz)*an(ix,jy,kz,lh)
77753 continue
77752 continue
77751 continue
c
c  set factors based on ratio of (totn+a)/(totp+a) where a =1.e-20
c
      csfac =
     > (cstotn+1.e-20)/(cstotp+1.e-20)
      crfac =
     > (crtotn+1.e-20)/(crtotp+1.e-20)
      cirfac =
     > (cirtotn+1.e-20)/(cirtotp+1.e-20)
      cipfac =
     > (ciptotn+1.e-20)/(ciptotp+1.e-20)
      cifac =
     > (citotn+1.e-20)/(citotp+1.e-20)
      ccfac =
     > (cctotn+1.e-20)/(cctotp+1.e-20)
      cvfac =
     > (cvtotn+1.e-20)/(cvtotp+1.e-20)
      cffac =
     > (cftotn+1.e-20)/(cftotp+1.e-20)
      cglfac =
     > (cgltotn+1.e-20)/(cgltotp+1.e-20)
      cgmfac =
     > (cgmtotn+1.e-20)/(cgmtotp+1.e-20)
      cghfac =
     > (cghtotn+1.e-20)/(cghtotp+1.e-20)
      chfac =
     > (chtotn+1.e-20)/(chtotp+1.e-20)
c
c  now reset all values based on factors
c
      do 77771 kz =  1, nz-kd1
      do 77772 jy =  1, ny-jd1
      do 77773 ix =  1, nx-id1
      an(ix,jy,kz,lc) = an(ix,jy,kz,lc)*ccfac
      an(ix,jy,kz,li) = an(ix,jy,kz,li)*cifac
      an(ix,jy,kz,lip) = an(ix,jy,kz,lip)*cipfac
      an(ix,jy,kz,lir) = an(ix,jy,kz,lir)*cirfac
      an(ix,jy,kz,lr) = an(ix,jy,kz,lr)*crfac
c BUG!      an(ix,jy,kz,ls) = an(ix,jy,kz,ls)*crfac
      an(ix,jy,kz,ls) = an(ix,jy,kz,ls)*csfac
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv)*cvfac
      an(ix,jy,kz,lf) = an(ix,jy,kz,lf)*cffac
      an(ix,jy,kz,lgl) = an(ix,jy,kz,lgl)*cglfac
      an(ix,jy,kz,lgm) = an(ix,jy,kz,lgm)*cgmfac
      an(ix,jy,kz,lgh) = an(ix,jy,kz,lgh)*cghfac
      an(ix,jy,kz,lh) = an(ix,jy,kz,lh)*chfac
77773 continue
77772 continue
77771 continue
c
      end if
c
c  end of hole filling version 2
c        
      endif
c
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
         
            t0(ix,jy,kz) = an(ix,jy,kz,lt)*t77(ix,jy,kz)

        ENDDO
       ENDDO
      ENDDO


      IF ( ipconc .ge. 2 ) THEN
      do kz = 1,nz-1
       do jy = 1,ny-jstag
        do ix = 1,nx-istag

         temp1 = t0(ix,jy,kz)

         ltemq = Int( (temp1-163.15)/fqsat+1.5 )

c          c1 = t00(ix,jy,kz)*tabqvs(ltemq)
          c1 = t00(ix,jy,kz)*tabqis(ltemq)

          ssati(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values

          c1 = t00(ix,jy,kz)*tabqvs(ltemq)

          ssfilt(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
          IF ( ny .eq. 2 .and. ssfilt(ix,jy,kz) .gt. 1.0
     :                             .and. ncdebug .ge. 1 ) THEN
          write(iunit,*) 'ssfilt2 = ',ssfilt(ix,jy,kz),
     :     an(ix,jy,kz,lv)*1.e3,c1*1.e3,an(ix,jy,kz,lt),ix,kz,
     :     an(ix,jy,kz,lnc)*1.e-6,an(ix,jy,kz,lc)*1.e3
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      ENDIF
c 
c
c Redistribute inappreciable cloud particles and charge
c
c Redistribution everywhere in the domain...
c
      sctot1n = 0.0
      sctot1p = 0.0
      frac = 1.0 ! 0.2
c cmic$ do all autoscope
c cmic$1 private(kz,jy,ix),
c cmic$2 shared(an)
c C$DOACROSS LOCAL(kz,jy,ix),
c C$&        SHARE(an)      
!! c$omp PARALLEL DO DEFAULT(SHARED), PRIVATE(kz,jy,ix),SHARED(an)  

!      IF ( ipconc .lt. 1 ) THEN
      IF ( .false. ) THEN

      do kz = 1,nz-1
      do jy = 1,ny
      do ix = 1,nx
       

      if ( an(ix,jy,kz,lhl) .lt.  frac*qxmin(lh) 
c     :    .and. (an(ix,jy,kz,lhl) .gt. 0.0 .or.
c     :      an(ix,jy,kz,lschl) .ne. 0.0 ) 
     :       ) then
        
        IF ( an(ix,jy,kz,lhl) .gt. 0 ) THEN
          an(ix,jy,kz,lh) = an(ix,jy,kz,lh) + an(ix,jy,kz,lhl)
          an(ix,jy,kz,lhl) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnhl) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnh) = an(ix,jy,kz,lnh) + an(ix,jy,kz,lnhl)
          an(ix,jy,kz,lnhl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lschl) .ne. 0.0) THEN
          an(ix,jy,kz,lsch) = an(ix,jy,kz,lsch) + an(ix,jy,kz,lschl)
          an(ix,jy,kz,lschl) = 0.0
        ENDIF
      end if

      if ( an(ix,jy,kz,lh) .lt. frac*qxmin(lh) 
c     : .and.
c     :     ( an(ix,jy,kz,lh) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lsch) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lnh) .gt. 0.0 ) 
     :        ) then

        IF ( an(ix,jy,kz,lh) .gt. 0 ) THEN
          an(ix,jy,kz,lf) = an(ix,jy,kz,lf) + an(ix,jy,kz,lh)
          an(ix,jy,kz,lh) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lnh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnf) = an(ix,jy,kz,lnf) + an(ix,jy,kz,lnh)
          an(ix,jy,kz,lnh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2  ) THEN ! .and. an(ix,jy,kz,lsch) .ne. 0.0) THEN
          an(ix,jy,kz,lscf) = an(ix,jy,kz,lscf) + an(ix,jy,kz,lsch)
          an(ix,jy,kz,lsch) = 0.0
        ENDIF
      end if

      if ( an(ix,jy,kz,lf) .lt.  frac*qxmin(lf) !.and.
c     :     (an(ix,jy,kz,lf) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscf) .ne. 0.0 .or. 
c     :      an(ix,jy,kz,lnf) .gt. 0.0 ) 
     :        ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN

        IF ( an(ix,jy,kz,lf) .gt. 0 ) THEN
          an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,lf)
          an(ix,jy,kz,lf) = 0.0
        ENDIF

        IF ( ipconc .ge. 5  ) THEN ! .and. an(ix,jy,kz,lnf) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lnf)
          an(ix,jy,kz,lnf) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscf) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscf)
          an(ix,jy,kz,lscf) = 0.0
        ENDIF
      ELSE

        IF ( an(ix,jy,kz,lf) .gt. 0 ) THEN
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + an(ix,jy,kz,lf)
          an(ix,jy,kz,lf) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lnf) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lnf)
          an(ix,jy,kz,lnf) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscf) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscf)
          an(ix,jy,kz,lscf) = 0.0
        ENDIF
      ENDIF

      end if

      if ( an(ix,jy,kz,lgl) .lt. frac*qxmin(lgl) ! .and.
c     :     (an(ix,jy,kz,lgl) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscgl) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngl) .gt. 0.0) 
     :        ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN

        IF ( an(ix,jy,kz,lgl) .gt. 0 ) THEN
          an(ix,jy,kz,lir) = an(ix,jy,kz,lir) + an(ix,jy,kz,lgl)
          an(ix,jy,kz,lgl) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngl) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngl)
          an(ix,jy,kz,lngl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgl) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgl)
          an(ix,jy,kz,lscgl) = 0.0
        ENDIF
      ELSE
        IF ( an(ix,jy,kz,lgl) .gt. 0 ) THEN
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + an(ix,jy,kz,lgl)
          an(ix,jy,kz,lgl) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngl) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngl)
          an(ix,jy,kz,lngl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgl) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgl)
          an(ix,jy,kz,lscgl) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,lgm) .lt. frac*qxmin(lgm)  !.and.
c     :     (an(ix,jy,kz,lgm) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscgm) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngm) .gt. 0.0) 
     :       ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
        IF ( an(ix,jy,kz,lgm) .gt. 0 ) THEN
          an(ix,jy,kz,lir) = an(ix,jy,kz,lir) + an(ix,jy,kz,lgm)
          an(ix,jy,kz,lgm) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngm) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngm)
          an(ix,jy,kz,lngm) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgm) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgm)
          an(ix,jy,kz,lscgm) = 0.0
        ENDIF
      ELSE
        IF ( an(ix,jy,kz,lgm) .gt. 0 ) THEN
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + an(ix,jy,kz,lgm)
          an(ix,jy,kz,lgm) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngm) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngm)
          an(ix,jy,kz,lngm) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgm) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgm)
          an(ix,jy,kz,lscgm) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,lgh) .lt. frac*qxmin(lgh) !.and.
c     :     ( an(ix,jy,kz,lgh) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscgh) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngh) .gt. 0.0) 
     :       ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
        IF ( an(ix,jy,kz,lgh) .gt. 0 ) THEN
          an(ix,jy,kz,lir) = an(ix,jy,kz,lir) + an(ix,jy,kz,lgh)
          an(ix,jy,kz,lgh) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngh)
          an(ix,jy,kz,lngh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgh) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgh)
          an(ix,jy,kz,lscgh) = 0.0
        ENDIF
      ELSE
        IF ( an(ix,jy,kz,lgh) .gt. 0 ) THEN
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + an(ix,jy,kz,lgh)
          an(ix,jy,kz,lgh) = 0.0
        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lngh) .gt. 0.0 ) THEN
          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngh)
          an(ix,jy,kz,lngh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgh) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgh)
          an(ix,jy,kz,lscgh) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,ls) .lt.  frac*qxmin(ls) ! .and.
c     :     ( an(ix,jy,kz,ls) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lscs) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lns) .gt. 0.0) 
     :         ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
        ENDIF
      
        IF ( ipconc .ge. 4 ) THEN !  .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
          an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF

        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscs) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscs)
          an(ix,jy,kz,lscs) = 0.0
        ENDIF
      ELSE
        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lr) = an(ix,jy,kz,lr) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
        ENDIF

        IF ( ipconc .ge. 4 ) THEN !  .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscs) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscs)
          an(ix,jy,kz,lscs) = 0.0
        ENDIF
      ENDIF
      end if

c      IF ( an(ix,jy,kz,lr)  .gt. 10.0e-3 )
c     :  print*, 'RAIN3: ',ix,kz,an(ix,jy,kz,lr)

      if ( an(ix,jy,kz,lr) .lt. frac*qxmin(lr)
c     :     .and. ( an(ix,jy,kz,lr) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscr) .ne. 0.0 ) 
     :  ) then
        an(ix,jy,kz,lc) = an(ix,jy,kz,lc) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        IF ( ipconc .ge. 3 ) THEN
          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lscw) = an(ix,jy,kz,lscw) + an(ix,jy,kz,lscr)
          an(ix,jy,kz,lscr) = 0.0
        ENDIF
      end if

c for qip        
        
      if ( an(ix,jy,kz,lip) .lt. frac*qxmin(li) 
c     :    .and. ( an(ix,jy,kz,lip) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscip) .ne. 0.0 )
     : ) then
        an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,lip)
        an(ix,jy,kz,lip) = 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lnip) = 0.0
       ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscip)
          an(ix,jy,kz,lscip) = 0.0
        ENDIF
      end if
      
      if ( an(ix,jy,kz,lir) .lt. frac*qxmin(li)
c     :     .and. ( an(ix,jy,kz,lir) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscir) .ne. 0.0 ) 
     :    ) then
        an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,lir)
        an(ix,jy,kz,lir) = 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lnir) = 0.0
       ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscir)
          an(ix,jy,kz,lscir) = 0.0
        ENDIF
      end if
      
c
c  for qci
c
      IF ( an(ix,jy,kz,li) .le. frac*qxmin(li) 
c     :    .and. ( an(ix,jy,kz,li) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lsci) .ne. 0.0 ) 
     :    ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
      an(ix,jy,kz,li)= 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lni) = 0.0
       ENDIF
       IF ( ipelec .ge. 2 ) THEN
        IF ( an(ix,jy,kz,lsci) .gt. 0.0 ) THEN
         an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lsci)*eci
         an(ix,jy,kz,lsci) = 0.0
        ELSE
         an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lsci)*eci
         an(ix,jy,kz,lsci) = 0.0
        ENDIF
       ENDIF
      ENDIF
      
c
c  for qcw
c
      IF ( an(ix,jy,kz,lc) .le. frac*qxmin(lc) 
c     :    .and. ( an(ix,jy,kz,lc) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lscw) .ne. 0.0 ) 
     :       ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0
       IF ( ipconc .ge. 2 ) THEN
        IF ( lccn .gt. 1 ) THEN
         an(ix,jy,kz,lccn) =  
     :      Min( ccwmx, an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc)) )
        ENDIF
         an(ix,jy,kz,lnc) = 0.0
       
       IF ( lss .gt. 1 ) THEN
       IF ( idecss .eq. 1 .and. ( an(ix,jy,kz,lni) +
     :      an(ix,jy,kz,lnip)) .lt. 10.e3 
     :     .and. an(ix,jy,kz,lss) .gt. 0.0 ) THEN
         an(ix,jy,kz,lss) = an(ix,jy,kz,lss)*Exp(-dtp/10.)
         IF ( ssat(ix,jy,kz) .lt. -5.0 ) THEN
          an(ix,jy,kz,lss)  = 0.0 
         ENDIF 
       ENDIF
       ENDIF
       
       ENDIF

       IF ( ipelec .gt. 0 ) THEN
        IF ( an(ix,jy,kz,lscw) .gt. 0.0 ) THEN
         an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscw)*eci
         an(ix,jy,kz,lscw) = 0.0
        ELSE
         an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscw)*eci
         an(ix,jy,kz,lscw) = 0.0
        ENDIF
       ENDIF
      ENDIF

      end do
      end do
      end do
      
      ELSE
!
!  alternate test version for ipconc .ge. 3
!  just vaporize stuff to prevent noise in the number concentrations
!
      do kz = 1,nz-1
      do jy = 1,ny
      do ix = 1,nx

      if ( an(ix,jy,kz,lhl) .lt.  frac*qxmin(lhl) 
c     :    .and. (an(ix,jy,kz,lhl) .gt. 0.0 .or.
c     :      an(ix,jy,kz,lschl) .ne. 0.0 ) 
     :       ) then
        
!        IF ( an(ix,jy,kz,lhl) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lhl)
          an(ix,jy,kz,lhl) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnhl) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnh) = an(ix,jy,kz,lnh) + an(ix,jy,kz,lnhl)
          an(ix,jy,kz,lnhl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lschl) .ne. 0.0) THEN
          an(ix,jy,kz,lsch) = an(ix,jy,kz,lsch) + an(ix,jy,kz,lschl)
          an(ix,jy,kz,lschl) = 0.0
        ENDIF
      end if

      if ( an(ix,jy,kz,lh) .lt. frac*qxmin(lh) 
c     : .and.
c     :     ( an(ix,jy,kz,lh) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lsch) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lnh) .gt. 0.0 ) 
     :        ) then

!        IF ( an(ix,jy,kz,lh) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lh)
          an(ix,jy,kz,lh) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnf) = an(ix,jy,kz,lnf) + an(ix,jy,kz,lnh)
          an(ix,jy,kz,lnh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lsch) .ne. 0.0) THEN
          an(ix,jy,kz,lscf) = an(ix,jy,kz,lscf) + an(ix,jy,kz,lsch)
          an(ix,jy,kz,lsch) = 0.0
        ENDIF
      end if

      if ( an(ix,jy,kz,lf) .lt.  frac*qxmin(lf) !.and.
c     :     (an(ix,jy,kz,lf) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscf) .ne. 0.0 .or. 
c     :      an(ix,jy,kz,lnf) .gt. 0.0 ) 
     :        ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN

!        IF ( an(ix,jy,kz,lf) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lf)
          an(ix,jy,kz,lf) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lnf) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lnf)
          an(ix,jy,kz,lnf) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscf) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscf)
          an(ix,jy,kz,lscf) = 0.0
        ENDIF
      ELSE

!        IF ( an(ix,jy,kz,lf) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lf)
          an(ix,jy,kz,lf) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN !  .and. an(ix,jy,kz,lnf) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lnf)
          an(ix,jy,kz,lnf) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscf) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscf)
          an(ix,jy,kz,lscf) = 0.0
        ENDIF
      ENDIF

      end if

      if ( an(ix,jy,kz,lgl) .lt. frac*qxmin(lgl) ! .and.
c     :     (an(ix,jy,kz,lgl) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscgl) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngl) .gt. 0.0) 
     :        ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN

!        IF ( an(ix,jy,kz,lgl) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgl)
          an(ix,jy,kz,lgl) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngl) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngl)
          an(ix,jy,kz,lngl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgl) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgl)
          an(ix,jy,kz,lscgl) = 0.0
        ENDIF
      ELSE
!        IF ( an(ix,jy,kz,lgl) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgl)
          an(ix,jy,kz,lgl) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngl) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngl)
          an(ix,jy,kz,lngl) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgl) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgl)
          an(ix,jy,kz,lscgl) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,lgm) .lt. frac*qxmin(lgm)  !.and.
c     :     (an(ix,jy,kz,lgm) .ge. 0.0 .or.
c     :      an(ix,jy,kz,lscgm) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngm) .gt. 0.0) 
     :       ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
!        IF ( an(ix,jy,kz,lgm) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgm)
          an(ix,jy,kz,lgm) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngm) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngm)
          an(ix,jy,kz,lngm) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgm) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgm)
          an(ix,jy,kz,lscgm) = 0.0
        ENDIF
      ELSE
!        IF ( an(ix,jy,kz,lgm) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgm)
          an(ix,jy,kz,lgm) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngm) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngm)
          an(ix,jy,kz,lngm) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgm) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgm)
          an(ix,jy,kz,lscgm) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,lgh) .lt. frac*qxmin(lgh) !.and.
c     :     ( an(ix,jy,kz,lgh) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscgh) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lngh) .gt. 0.0) 
     :       ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
!        IF ( an(ix,jy,kz,lgh) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgh)
          an(ix,jy,kz,lgh) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngh) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnir) = an(ix,jy,kz,lnir) + an(ix,jy,kz,lngh)
          an(ix,jy,kz,lngh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgh) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscgh)
          an(ix,jy,kz,lscgh) = 0.0
        ENDIF
      ELSE
!        IF ( an(ix,jy,kz,lgh) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lgh)
          an(ix,jy,kz,lgh) = 0.0
!        ENDIF

        IF ( ipconc .ge. 5 ) THEN ! .and. an(ix,jy,kz,lngh) .gt. 0.0 ) THEN
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lngh)
          an(ix,jy,kz,lngh) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscgh) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscgh)
          an(ix,jy,kz,lscgh) = 0.0
        ENDIF
      ENDIF
      end if

      if ( an(ix,jy,kz,ls) .lt.  frac*qxmin(ls) ! .and.
c     :     ( an(ix,jy,kz,ls) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lscs) .ne. 0.0 .or. 
c     :       an(ix,jy,kz,lns) .gt. 0.0) 
     :         ) then
      IF ( t0(ix,jy,kz) .lt. 273.15 ) THEN
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF
      
        IF ( ipconc .ge. 4 ) THEN ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
!          an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF

        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscs) .ne. 0.0) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscs)
          an(ix,jy,kz,lscs) = 0.0
        ENDIF
      ELSE
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF

        IF ( ipconc .ge. 4 ) THEN ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) THEN ! 
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN !  .and. an(ix,jy,kz,lscs) .ne. 0.0) THEN
          an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + an(ix,jy,kz,lscs)
          an(ix,jy,kz,lscs) = 0.0
        ENDIF
      ENDIF
      end if

c      IF ( an(ix,jy,kz,lr)  .gt. 10.0e-3 )
c     :  print*, 'RAIN3: ',ix,kz,an(ix,jy,kz,lr)

      if ( an(ix,jy,kz,lr) .lt. frac*qxmin(lr)
c     :     .and. ( an(ix,jy,kz,lr) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscr) .ne. 0.0 ) 
     :  ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        IF ( ipconc .ge. 3 ) THEN
!          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
        ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lscw) = an(ix,jy,kz,lscw) + an(ix,jy,kz,lscr)
          an(ix,jy,kz,lscr) = 0.0
        ENDIF
      end if

c for qip        
        
      if ( an(ix,jy,kz,lip) .lt. frac*qxmin(lip) 
c     :    .and. ( an(ix,jy,kz,lip) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscip) .ne. 0.0 )
     : ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lip)
        an(ix,jy,kz,lip) = 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lnip) = 0.0
       ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscip)
          an(ix,jy,kz,lscip) = 0.0
        ENDIF
      end if
      
      if ( an(ix,jy,kz,lir) .lt. frac*qxmin(lir)
c     :     .and. ( an(ix,jy,kz,lir) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lscir) .ne. 0.0 ) 
     :    ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lir)
        an(ix,jy,kz,lir) = 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lnir) = 0.0
       ENDIF
      
        IF ( ipelec .ge. 2 ) THEN
          an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + an(ix,jy,kz,lscir)
          an(ix,jy,kz,lscir) = 0.0
        ENDIF
      end if
      
c
c  for qci
c
      IF ( an(ix,jy,kz,li) .le. frac*qxmin(li) 
c     :    .and. ( an(ix,jy,kz,li) .ge. 0.0 .or.
c     :       an(ix,jy,kz,lsci) .ne. 0.0 ) 
     :    ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
      an(ix,jy,kz,li)= 0.0
       IF ( ipconc .ge. 1 ) THEN
         an(ix,jy,kz,lni) = 0.0
       ENDIF
       IF ( ipelec .ge. 2 ) THEN
        IF ( an(ix,jy,kz,lsci) .gt. 0.0 ) THEN
         an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lsci)*eci
         an(ix,jy,kz,lsci) = 0.0
        ELSE
         an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lsci)*eci
         an(ix,jy,kz,lsci) = 0.0
        ENDIF
       ENDIF
      ENDIF
      
c
c  for qcw
c
      IF ( an(ix,jy,kz,lc) .le. frac*qxmin(lc) 
c     :    .and. ( an(ix,jy,kz,lc) .gt. 0.0 .or.
c     :       an(ix,jy,kz,lscw) .ne. 0.0 ) 
     :       ) THEN
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0
       IF ( ipconc .ge. 2 ) THEN
        IF ( lccn .gt. 1 ) THEN
         an(ix,jy,kz,lccn) =  
     :      Min( ccwmx, an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc)) )
        ENDIF
         an(ix,jy,kz,lnc) = 0.0
       
       IF ( lss .gt. 1 ) THEN
       IF ( idecss .eq. 1 .and. ( an(ix,jy,kz,lni) +
     :      an(ix,jy,kz,lnip)) .lt. 10.e3 
     :     .and. an(ix,jy,kz,lss) .gt. 0.0 ) THEN
         an(ix,jy,kz,lss) = an(ix,jy,kz,lss)*Exp(-dtp/10.)
         IF ( ssat(ix,jy,kz) .lt. -5.0 ) THEN
          an(ix,jy,kz,lss)  = 0.0 
         ENDIF 
       ENDIF
       ENDIF
       
       ENDIF

       IF ( ipelec .gt. 0 ) THEN
        IF ( an(ix,jy,kz,lscw) .gt. 0.0 ) THEN
         an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) + an(ix,jy,kz,lscw)*eci
         an(ix,jy,kz,lscw) = 0.0
        ELSE
         an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) - an(ix,jy,kz,lscw)*eci
         an(ix,jy,kz,lscw) = 0.0
        ENDIF
       ENDIF
      ENDIF

      end do
      end do
      end do
      
      ENDIF

c
c
c
c  return   
c     
C       IF ( ipconc .ge. 1 ) THEN
C      do kz = 2,nz-2
C      do jy = 2,ny-2
C      do ix = 2,nx-2
C        
C        IF ( t0(ix,jy,kz) .lt. tfrh .and. an(ix,jy,kz,lni) .gt. 200.e3 ) THEN
C         IF ( an(ix+1,jy,kz,lni) .lt. 1.e3 .and. an(ix-1,jy,kz,lni) .lt. 1.e3 .and.
C     :        an(ix,jy+1,kz,lni) .lt. 1.e3 .and. an(ix,jy+1,kz,lni) .lt. 1.e3 ) THEN
C          write(iunit,*) 'ICE10: isolated lni at ix,jy,kz ',ix,jy,kz
C          write(iunit,*) 'c,q = ',an(ix,jy,kz,lni),an(ix,jy,kz,li),qxmin(li)
C         ENDIF
C        ELSEIF ( an(ix,jy,kz,lni) .lt. 0.0 ) THEN
C          write(iunit,*) 'ICE10: negative lni at ix,jy,kz ',ix,jy,kz
C          write(iunit,*) 'c,q = ',an(ix,jy,kz,lni),an(ix,jy,kz,li),qxmin(li)
C        ENDIF
C         
C      end do
C      end do
C      end do
C       ENDIF

      deallocate ( db1 )
      deallocate ( db0 )
      deallocate ( xvt )
      deallocate ( tmpn )
      
      return
      end
c
c  end of subroutine
c
c--------------------------------------------------------------------------
c
c  Bottom of code
c
#endif
