module my_tmom_mod
!module my3mom_main_mod  !future name; need to change vkuocon6

! DTD (07/18/2011).  This module has been modified for the interface with ARPS
! The following changes have been made (most originally made by Yunheng):
! 1) Changes to the argument list to pass in a 4D scalar array for the microphysics variables
!    instead of individual 3D arrays for each category.  Internally, pointer variables to the 
!    various portions of the scalar array are declared and used.
! 2) Instead of the surface pressure and sigma vertical coordinate values, the full column
!    of pressure is passed into the subroutine.  Appropriate changes have been made in the code
!    where pressure at a given height is needed (i.e. removed calculation of pressure, since it 
!    is now passed in).  Similarly, instead of geopotential height, physical height is now
!    passed in and dealt with appropriately in the code. Finally, the vertical velocity on height
!    surfaces is passed in instead of omega (on pressure surfaces), since ARPS uses a height-based
!    vertical coordinate system anyway.
! 3) The original code passes back out the microphysics tendency terms to the model.  We don't 
!    need these in ARPS, so instead the adjusted variables (at time tfuture) are passed back out
!    directly.  Also, no approximation at time (t*) is needed, since we pass in the present time
!    level explicitly.
! 4) The three include files (consphy.cdk, dintern.cdk, and fintern.cdk) have been explicitly
!    included in the subroutine, consistent with the WRF interface, replacing their #include 
!    directives in the original code.
! 4) Various other minor changes, such as flags to allow for removal of graupel and hail (originally
!    implemented by Bryan Putnam and Youngsun Jung), and to disallow size-sorting in the double-
!    and triple-moment scheme, among others. EDIT: Jason independently added a hail_ON switch, so
!    I used that instead in this version, but changed the name to hl_ON to be consistent with the
!    previous implementation and to avoid interfering with the ARPS namelist variable hail_ON.
!    Also, the vertical coordinate index in the original code is assumed to increase from upper to
!    lower levels (i.e. toward the ground), which is opposite to the ARPS convention.  Therefore, the
!    ordering has been reversed in several places to reflect this, such as in the new calculation of the
!    top of the sedimentation column in keeping with the original changes in this regard by Yunheng.

  implicit none

  private
  public :: mytmom_main

  contains

!_______________________________________________________________________________________!

! SUBROUTINE MYTMOM_MAIN(Womega,T,Q,QC,QR,QI,QN,QG,QH,NC,NR,NY,NN,NG,NH,ZR,ZI,      &
!    ZN,ZG,ZH,PS,TM,QM,QCM,QRM,QIM,QNM,QGM,QHM,NCM,NRM,NYM,NNM,NGM,NHM,ZRM,ZIM,     &
!    ZNM,ZGM,ZHM,PSM,S,LR,SR,GZ,T_TEND,Q_TEND,QCTEND,QRTEND,QITEND,QNTEND,QGTEND,   &
!    QHTEND,NCTEND,NRTEND,NYTEND,NNTEND,NGTEND,NHTEND,ZRTEND,ZITEND,ZNTEND,ZGTEND,  &
!    ZHTEND,DT_sp,NI,N,NK,J,KOUNT,scheme,SS01,SS02,SS03,SS04,SS05,SS06,SS07,SS08,   &
!    SS09,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20)

  SUBROUTINE MYTMOM_MAIN(WW,T,Q,QNZ,P,TM,QM,QNZM,PM,LR,SR,NCtem,ZP,DT_sp,NI,NK,J,   &
    scheme,nscalar,ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,      &
    alpharain,alphaice,alphasnow,alphagrpl,alphahail,ixe)

  use my2mom_fncs_mod
  use my3mom_fncs_mod
  use my3mom_sedi_mod

  IMPLICIT NONE
  
  !INCLUDE 'globcst.inc'

!CALLING PARAMETERS:
  integer,              intent(in)    :: NI,NK,J,scheme,nscalar
  real, dimension(ni,nk), intent(in)    :: WW,ZP,P,PM
  real,                 intent(in)    :: DT_sp
  real, dimension(ni),   intent(out)   :: LR,SR
  real, dimension(ni,nk,nscalar), intent(inout), target :: qnz,qnzm
  REAL, INTENT(INOUT), TARGET :: NCtem(ni,nk,12)
  real, dimension(ni,nk), intent(inout) :: T,Q,TM,QM
!  real, dimension(:,:), intent(inout) :: T,Q,TM,QM,QC,QCM,NC,NCM,QR,QRM,NR,NRM,     &
!        QI,QIM,NY,NYM,QN,QNM,NN,NNM,QG,QGM,NG,NGM,QH,QHM,NH,NHM,ZR,ZRM,ZI,ZIM,ZN,   &
!        ZNM,ZG,ZGM,ZH,ZHM,T_TEND,QCTEND,QRTEND,QITEND,QNTEND,QGTEND,QHTEND,Q_TEND,  &
!        NCTEND,NRTEND,NYTEND,NNTEND,NGTEND,NHTEND,ZRTEND,ZITEND,ZNTEND,ZGTEND,      &
!        ZHTEND,SS01,SS02,SS03,SS04,SS05,SS06,SS07,SS08,SS09,SS10,SS11,SS12,SS13,    &
!        SS14,SS15,SS16,SS17,SS18,SS19,SS20
!  real, dimension(ni,nk) intent(inout) :: SS01,SS02,SS03,SS04,SS05,SS06,SS07,SS08,     &
!        SS09,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,SS18,SS19,SS20

  ! The following must be passed in explicitly since phycst.inc conflicts with consphy_v43.cdk

  REAL, INTENT(IN) :: ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow
  REAL, INTENT(IN) :: rhogrpl,rhohail,alpharain,alphaice,alphasnow,alphagrpl
  REAL, INTENT(IN) :: alphahail
  integer, intent(IN) :: ixe

!____________________________________________________________________________________________
!                                                                                            !
!                Milbrandt-Yau (2005) Multi-Moment Bulk Microphysics Scheme                  !
!                              -  Triple-Moment Version  -                                   !
!____________________________________________________________________________________________!
!  Package version:   2.20.0 (internal bookkeeping)                                          !
!  Last modified  :   2011-03-02                                                             !
!____________________________________________________________________________________________!

!____________________________________________________________________________________________!
!  Author         :   Jason Milbrandt  (McGill University / RPN-A)                           !
!                                                                                            !
!  Revision:                                                                                 !
!                                                                                            !
!  001  J. Milbrandt  (Dec 2004) - optimization and code clean-up for use on IBM;            !
!        [RPN]                     implemented box-Lagrangian sedimentation                  !
!  002  J. Milbrandt  (Dec 2007) - modifications for interface with GEM;  bug-fixes          !
!                                                                                            !
!                                                                                            !
!  Object:                                                                                   !
!          Computes changes to the temperature, water vapor mixing ratio, and the            !
!          mixing ratios of six hydrometeor species resulting from cloud microphysical       !
!          interactions at saturated grid points. Surface precipitation rates from each      !
!          sedimenting hydrometeor categories are also computed.                             !
!                                                                                            !
!                                                                                            !
!  This code and the associated modules form the multi-moment bulk microphysics scheme       !
!  described in the references below.  The scheme computes the tendencies of the hydrometeor !
!  moments, temperature, and humidity on a NI x NJ slab.                                     !
!                                                                                            !
!  The current version of the code is written such that the user can switch between the      !
!  various versions of the scheme -- single-moment, double-moment (fixed- or diagnosed-      !
!  dispersion parameter), and triple-moment -- with a namelist switch ('my_full_version').   !
!  This enables testing and maintenance of single version of code.  Note, the resulting code !
!  is thus less computationally efficient (particularly for the single-moment version for    !
!  which all concentration (NX)-tendency equations are unnecessarily computed).              !
!                                                                                            !
!  References:   Milbrandt and Yau, (2005a): [Part I ] J.Atmos.Sci., vol.62, 3051-3064       !
!                --------- and ---, (2005b): [Part II] J.Atmos.Sci., vol.62, 3065-3081       !
!                (and references therein)                                                    !
!                                                                                            !
!  Please report bugs to:  jason.milbrandt@ec.gc.ca                                          !
!____________________________________________________________________________________________!
!
! Arguments:         Description:                                         Units:
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Input -
!
! NI                 number of x-dir points (local subdomain)
! NK                 number of vertical levels
! N
! J                  y-dir index (local subdomain)
! KOUNT              current model time step number
! DT_sp              model time step                                      [s]
! Womega             vertical velocity                                    [Pa s-1]
! S                  sigma
! GZ                 geopotential height                                  [m]
! scheme             scheme version
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Input/Output -
!
! T                  air temperature at time (t*)                         [K]
! TM                 air temperature at time (t-dt)                       [K]
! Q                  water vapor mixing ratio at (t*)                     [kg kg-1]
! QM                 water vapor mixing ratio at (t-dt)                   [kg kg-1]
! PS                 surface pressure at time (t*)                        [Pa]
! PSM                surface pressure at time (t-dt)                      [Pa]
!
!  For x = (C,R,I,N,G,H):  C = cloud
!                          R = rain
!                          I = ice (pristine) [except 'NY', not 'NI']
!                          N = snow
!                          G = graupel
!                          H = hail
!
! Q(x)               mixing ratio for hydrometeor x at (t*)               [kg kg-1]
! Q(x)M              mixing ratio for hydrometeor x at (t-dt)             [kg kg-1]
! N(x)               total number concentration for hydrometeor x  (t*)   [m-3]
! N(x)M              total number concentration for hydrometeor x  (t-dt) [m-3]
! Z(x)               reflectivity for hydrometeor x  (t*)                 [m6 m-3]
! Z(x)M              reflectivityfor hydrometeor x  (t-dt)                [m6 m-3]
!
! Note:  The arrays "VM" (e.g. variables TM,QM,QCM etc.) are declared as INTENT(INOUT)
!        such that their values are modified in the code [VM = 0.5*(VM + V)].
!        This is to approxiate the values at time level (t), which are needed by
!        this routine but are unavailable to the PHYSICS.  The new values are discared
!        by the calling routine ('vkuocon6.ftn').  However, care should be taken with
!        interfacing with other modelling systems.  For GEM/MC2, it does not matter if
!        VM is modified since the calling module passes back only the tendencies
!        (VTEND) to the model.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!            - Output -
!
! Q_TEND             tendency for water vapor mixing ratio                [kg kg-1 s-1]
! T_TEND             tendency for air temperature                         [K s-1]
! Q(x)TEND           tendency for mixing ratio for hydrometeor x          [kg kg-1 s-1]
! N(x)TEND           tendency for number concentration for hydrometeor x  [m-3 s-1]
! Z(x)TEND           tendency for reflectivity for hydrometeor x          [m6 m-3 s-1]
! Dm_(x)             mean-mass diameter for hydrometeor x                 [m]
! LR                 precipitation rate (at sfc) of liquid rain (r)       [m+3 m-2 s-1]
! LS                 precipitation rate (at sfc) of total solid (i,s,g,h) [m+3 m-2 s-1]
! SSxx               S/S terms (for testing purposes)
!_______________________________________________________________________________________!


!LOCAL VARIABLES:

  ! Pointer variables to portions of the scalar array for each category

  REAL, POINTER :: QC(:,:), QR(:,:), QI(:,:), QN(:,:), QG(:,:), QH(:,:)
  REAL, POINTER :: NC(:,:), NR(:,:), NY(:,:), NN(:,:), NG(:,:), NH(:,:)
  REAL, POINTER ::          ZR(:,:), ZI(:,:), ZN(:,:), ZG(:,:), ZH(:,:)
  REAL, POINTER :: QCM(:,:),QRM(:,:),QIM(:,:),QNM(:,:),QGM(:,:),QHM(:,:)
  REAL, POINTER :: NCM(:,:),NRM(:,:),NYM(:,:),NNM(:,:),NGM(:,:),NHM(:,:)
  REAL, POINTER ::          ZRM(:,:),ZIM(:,:),ZNM(:,:),ZGM(:,:),ZHM(:,:)

! Parameters/variables to count active grid points:
  logical :: log1,log2,log3,log4,doneK,rainPresent,activePoint(ni,nk)

  integer :: i,k,niter,ll,start!,ktop_sedi
  integer, dimension(ni) :: ktop_sedi
  real,    dimension(ni,nk) :: DE,DP,QSS,QSW,QSI,DZ,RHOQX,iDE
  integer, dimension(nk)    :: FLIM
  real,    dimension(ni)    :: HPS
  real    :: rtmp,idt,TcSP,cmrSP,cmiSP,cmsSP,cmgSP,cmhSP,tmp1,tmp2,tmp3,tmp4

  real*8 :: dt,VDmax,NNUmax,X,D,DEL,QREVP,ES,DEdp,NuDEPSOR,NuCONTA,NuCONTB,NuCONTC, &
       NuCONT,GG,Na,Tcc,F1,F2,Kdiff,PSIa,Kn,source,sink,sour,ratio,qvs0,DELqvs,     &
       ft,esi,Si,Simax,Vq,Vn,Vz,GC1,GC2,GC3,GC4,GC5,GC6,GC7,GC8,GC11,               &
       GC12,GC13,GC14,GC15,GR1,GR2,GR3,GR13,GR14,GR15,GR16,GR17,GR36,GI1,GI3,GI4,   &
       GI5,GI6,GI7,GI8,GI36,GS7,GS8,GS30,ckQr1,ckQr2,ckQi1,ckQi2,ckQs1,ckQs2,ckQg1, &
       ckQg2,ckQh1,ckQh2,cexc9,cexr1,cexr2,cexr3,cexr4,cexr5,cexr6,cexr9,cexi1,     &
       cexi2,cexi9,cexs2,czr,czi,czs,czg,czh,cmr,cmi,cms,cmg,cmh,LAMr,Nor,Noi,Nos,  &
       Nog,iLAMr6,iLAMh2,iABi,ABw,VENTr,VENTs,VENTg,VENTi,Cdiff,Ka,MUdyn,MUkin,DEo, &
       gam,ScTHRD,Tc,mi,ri,ff,Ec,Ntr,Dho,DMrain,Ech,DMice,DMsnow,DMgrpl,DMhail,     &
       ssat,Swmax,dey,Esh,Eii,Eis,Ess,Eig,Eih,FRAC,JJ,Dirg,Dirh,Dsrs,Dsrg,Dsrh,     &
       Dgrg,Dgrh,SIGc,L,TAU,DrAUT,DrINIT,Di,Ds,Dg,Dh,qFact,nFact,Ki,ALFr,Rz,GR37,   &
       GI37,GS37,ckQr3,ckQi3,ckQs3,ckQg3,ckQh3,QCLcs,QCLrs,QCLis,QCLcg,QCLrg,QCLig, &
       QCLch,QCLrh,QCLsh,QMLir,QMLsr,QMLgr,QMLhr,QCLih,QCNig,QVDvg,QVDvh,QSHhr,     &
       QFZci,QNUvi,QCLci,QVDvi,QCNis,QCNis1,QCNis2,QCLir,QCLri,QCNsg,QCLsr,QCNgh,   &
       QCLgr,QHwet,QVDvs,QFZrh,QIMsi,QIMgi,NMLhr,NCNgh,NVDvh,NCLir,NCLri,NCLrh,     &
       NCLch,NCLsr,NCLirg,NCNig,NCLirh,NrFZrh,NhFZrh,NCLsrs,NCLsrg,NCLsrh,NCLgrg,   &
       NCLgrh,NVDvg,NMLgr,NiCNis,NsCNis,NVDvs,NMLsr,NCLsh,NCLss,NNUvi,NFZci,NVDvi,  &
       NCLis,NCLig,NCLih,NMLir,NCLrs,NCNsg,NCLcs,NCLcg,NCLci,NIMsi,NIMgi,NIMii,     &
       NCLgr,NCLrg,NSHhr,RCAUTR,RCACCR,CCACCR,CCSCOC,CCAUTR,CRSCOR,ALFx,GX2,GX5,    &
       LAMx,iLAMx,iLAMxB0,Dx,ffx,iMUc,icmr,icmi,icms,icmg,icmh,                     &
       tmpdp1,tmpdp2,tmpdp3,tmpdp4,tmpdp5,tmpdp6,tmpdp7,tmpdp8,tmpdp9,tmpdp10

  real*8, dimension(ni,nk) :: iLAMr,iLAMr2,iLAMr3,iLAMr4,   &
       iLAMr5,iLAMc,iLAMc2,iLAMc3,iLAMc4,iLAMc5,iLAMc6,iLAMi,iLAMi2,iLAMi3,iLAMi4,  &
       iLAMi5,iLAMiB0,iLAMiB1,iLAMiB2,iLAMs,iLAMs2,iLAMsB0,iLAMsB1,iLAMsB2,iLAMg,   &
       iLAMg2,iLAMgB0,iLAMgB1,iLAMgB2,iLAMh,iLAMhB0,iLAMhB1,iLAMhB2,ALFi,ALFs,ALFg, &
       ALFh,GR5,iGR5,GR31,GR32,GR33,GR34,GR35,GI2,GI31,GI32,GS2,GG2,GH2,Dc,Dr,vr0,  &
       vi0,vs0,vg0,vh0,Noh,VENTh


!Triple-moment variables:  (not needed for double-moment version)
  real   :: delZR,delZI,delZN,delZG,delZH
  real*8 :: ZrCLri,ZrCLrs,ZrCLrg,ZrCLrh,ZFZrh,ZCLyi,ZFZci,ZNUvi,ZiCNis,ZsCNis,      &
       ZiCNig,ZiCLir,ZiCLis,ZiCLig,ZiCLih,ZMLir,ZiIMsi,ZiIMgi,ZCNis,ZsCNsg,ZMLsr,   &
       ZsCLsr,ZsCLsh,ZsCLsrs,ZCLss,ZhCLirh,ZhCLsrh,ZhCLgrh,ZCLyg,ZMLgr,ZCNgh,ZCNig, &
       ZgCLirg,ZgCLsrg,ZgCLgrg,ZgCNsrg,ZCLyh,ZhMLhr,ZIMsi,ZIMgi,ZCLys,ZVDvh,ZgCNig, &
       ZgCNsg,ZrMLhr,GalphaI,GalphaS,GalphaG,GalphaH,GalphaRaut,Galpha0,Galpha1,    &
       Galpha2,Galpha3,Galpha4,Galpha5,DQrDt,DNrDt,DZrDt,DQiDt,DNiDt,DZiDt,DQnDt,   &
       DNnDt,DZnDt,DQgDt,DNgDt,DZgDt,DQhDt,DNhDt,DZhDt,DZrDt1,dZrDt2,DZrDt3

  real*8, dimension(ni,nk) :: GalphaR

!BLG (box-Lagrangian) sedimentation variables:
  integer :: nnn,npassr,npassi,npasss,npassg,npassh,a,counter
  real    :: dtr,dti,dtg,dts,dth,tmp,VqMax,VnMax,VzMax,cr6,ci6,cg6,cs6,ch6,ck7,     &
             CoMax,CoMAXr,CoMAXi,CoMAXs,CoMAXg,CoMAXh
  real, dimension(ni,nk) :: VVQ,VVN,VVZ,gamfact
  real, dimension(ni)    :: rrl,dum
  integer, dimension(ni) :: activeColumn
  logical                :: LOCALLIM,slabHASmass


  !   Cloud/Rain size distribution parameters
  !     Note: The symbols for MU and ALPHA are REVERSED from that of CP2000a,b
  !           Explicit appearance of MUr = 1. has been removed.

  real*8, parameter :: MUc= 3.d0, ALFc= 1.d0
  !------------------------------------!
  ! Symbol convention: (dist. params.) !
  !       MY05    F94       CP00       ! F94:  Ferrier, 1994 (JAS)
  !       ------  --------  ------     ! CP00: Cohard & Pinty, 2000a,b (QJGR)
  !       ALFx    ALPHAx    MUx-1      !
  !       MUx     (1)       ALPHAx     !
  !       ALFx+1  ALPHAx+1  MUx        !
  !------------------------------------!

  ! Fallspeed parameters:
  real*8, parameter :: afr= 4854.000d0,  bfr= 1.0000d0,  ffr= 195.d0 !Ferrier (1994)
  real*8, parameter :: afi=   71.340d0,  bfi= 0.6635d0   !Ferrier (1994)
  real*8, parameter :: afs=   11.720d0,  bfs= 0.4100d0   !Locatelli and Hobbs (1974)
  real*8, parameter :: afg=   19.300d0,  bfg= 0.3700d0   !Ferrier (1994)
  real*8, parameter :: afh=  206.890d0,  bfh= 0.6384d0   !Ferrier (1994)
  !Note:  Implicitly, ffx=0 for all x=i,s,g,h (as in F94)
  !real*8, parameter :: afs=    8.996d0,  bfs= 0.4200d0   ! previous

  real  , parameter :: epsQ  = 1.e-14   !min. allowable mixing ratio
  real  , parameter :: epsN  = 1.e-8    !min. allowable number concentration
  real  , parameter :: epsZ  = 1.e-32   !min. allowable reflectivity

  real  , parameter :: rthres= 0.1      !max. (delZx/ZX) (prevents truncation error for Z-tends)
  real*8, parameter :: epsDQ = 1.d-13   !min. allowable Q-tendency for calc. of Z-tendencies
  real*8, parameter :: iLAMmin1= 1.d-6  !min. iLAMx (prevents underflow in Nox and VENTx calcs)
  real*8, parameter :: iLAMmin2= 1.d-10 !min. iLAMx (prevents underflow in Nox and VENTx calcs)
  real*8, parameter :: eps   = 1.d-32
  real*8, parameter :: EPS9  = 1.d-6
  real*8, parameter :: k1    = 0.001d0
  real*8, parameter :: k2    = 0.0005d0
  real*8, parameter :: k3    = 2.54d0
  real*8, parameter :: CPW   = 4218.d0, CPI=2093.d0

!  real*8, parameter :: deg   =  400.d0, mgo= 1.6d-10
!  real*8, parameter :: deh   =  900.d0
!  real*8, parameter :: dei   =  500.d0, mio=1.d-12, Nti0=1.d3
  real*8, parameter :: dew   = 1000.d0
!  real*8, parameter :: des   =  100.d0, mso= 4.4d-10,  rso= 1.d-4

  ! DTD: bulk ice densities now specified in namelist input
  real*8 :: deg 
  real*8, parameter :: mgo= 1.6d-10
  real*8 :: deh
  real*8 :: dei
  real*8, parameter :: mio=1.d-12, Nti0=1.d3
  real*8 :: des
  real*8, parameter :: mso= 4.4d-10,  rso= 1.d-4
  real*8, parameter :: dmr   =    3.d0, dmi=3.d0, dms=3.d0, dmg=3.d0, dmh=3.d0
  real*8, parameter :: MUi   =    1.d0, MUs= 1.d0, MUg= 1.d0, MUh= 1.d0

! [SM] and [DM] size distribution parmaeters:
!  real*8, parameter :: ALFrfix  =  0.d0  !fixed ALPHA for SM rain
!  real*8, parameter :: ALFifix  =  0.d0  !fixed ALPHA for SM ice
!  real*8, parameter :: ALFsfix  =  0.d0  !fixed ALPHA for SM snow
!  real*8, parameter :: ALFgfix  =  0.d0  !fixed ALPHA for SM graupel
!  real*8, parameter :: ALFhfix  =  0.d0  !fixed ALPHA for SM hail

  ! DTD: These are now passed in from the calling program and are specified in the namelist input
  real*8 :: ALFrfix
  real*8 :: ALFifix
  real*8 :: ALFsfix
  real*8 :: ALFgfix
  real*8 :: ALFhfix

  real*8, parameter :: ALFiMAX  = 10.d0  !max alpha for ICE
  real*8, parameter :: ALFsMAX  = 10.d0  !max alpha for SNOW
  real*8, parameter :: ALFrMIN  =  0.d0  !min alpha for RAIN

!  real  , parameter :: Ncfix    =  1.e8  ![m-3] fixed NC  for SM cloud
!  real*8, parameter :: Norfix   =  1.d6  ![m-4] Fixed Nor for SM rain
!  real*8, parameter :: Nosfix   =  1.d7  ![m-4] Fixed Nos for SM snow
!  real*8, parameter :: Nogfix   =  4.d5  ![m-4] Fixed Nog for SM graupel
!  real*8, parameter :: Nohfix   =  1.d5  ![m-4] Fixed Noh for SM hail

  ! DTD: As above, these are now namelist variables (specified at runtime).
  real   :: Ncfix   ![m-3] fixed NC  for SM cloud
  real*8 :: Norfix  ![m-4] Fixed Nor for SM rain
  real*8 :: Nosfix  ![m-4] Fixed Nos for SM snow
  real*8 :: Nogfix  ![m-4] Fixed Nog for SM graupel
  real*8 :: Nohfix  ![m-4] Fixed Noh for SM hail

  ! NOTE: VxMAX below are the max.allowable mass-weighted fallspeeds for sedimentation.
  !       Thus, they are Vx corresponding to DxMAX times the max. density factor, GAM
  !       [GAMmax=sqrt(DEo/DEmin)=sqrt(1.25/0.4)~2.]  e.g. VrMAX=2.*8.m/s

  real*8, parameter :: DrMAX=  5.d-3;   real, parameter :: VrMAX= 16.
  real*8, parameter :: DiMAX= 10.d-3;   real, parameter :: ViMAX=  4.
  real*8, parameter :: DsMAX= 30.d-3;   real, parameter :: VsMAX=  6.
  real*8, parameter :: DgMAX= 50.d-3;   real, parameter :: VgMAX=  8.
  real*8, parameter :: DhMAX= 80.d-3;   real, parameter :: VhMAX= 40.

  real*8, parameter :: Eci= 1.d0, Ecs= 1.d0, Ecg= 1.d0
  real*8, parameter :: Eri= 1.d0, Ers= 1.d0, Erg= 1.d0, Erh= 1.d0
  real*8, parameter :: Xdisp   = 0.25d0         !dispersion of the fall velocity of ice crystals
  real*8, parameter :: aa11    = 9.44d15, aa22= 5.78d3, Rh= 41.d-6
  real*8, parameter :: Avx     = 0.78d0, Bvx= 0.30d0  !ventilation coefficients [F94 (B.36)]
  real*8, parameter :: Abigg   = 0.66d0, Bbigg= 100.d0!parameters in probabilistic freezing
  real*8, parameter :: fdialec = 0.224d0              !dialectric factor, |K|x^2/|K|w^2
  real*8, parameter :: Drshed  = 0.001d0              !mean diam. of drop shed during wet growth
  real*8, parameter :: SIGcTHRS= 15.d-6               !threshold cld std.dev. before autoconversion
  real*8, parameter :: KK1= 3.03d3, KK2= 2.59d15      !parameters in Long (1974) kernel
  real*8, parameter :: Dhh= 82.d-6                    !diameter that rain hump first appears
  real*8, parameter :: ALFrAUT= 2.d0                  !shape parameter of rain for autoconversion

  real*8, parameter :: Dr_3cmpThrs = 1.0d-3           !size threshold [m] for hail production from 3-comp freezing
  real*8, parameter :: Dg_CNgh     = 2.5d-3           !size threshold [m] for CNgh
  real*8, parameter :: r_CNgh      = 0.05d0           !Dg/Dho ratio threshold for CNgh
  real,   parameter :: w_CNgh      = 3.               !vertical motion  threshold [m s-1] for CNgh
  real,   parameter :: Qr_FZrh     = 1.0e-4           !qr-threshold [kg kg-1] for FZrh
  real,   parameter :: Tc_FZrh     = -10.             !temp-threshold (C) for FZrh
  real*8, parameter :: CNsgThres   = 1.0d0          !threshold for CLcs/VDvs ratio for CNsg -- PREVIOUSLY USED 1.d0
  real*8, parameter :: capFact     = 1.0d0          !capacitace of snow -- PREVIOUS USED 1.d0
  real,   parameter :: qReducFact  = 0.0            !reduction factor for supersaturation water vapor
  real,   parameter :: gzMax_sedi  = 180000.        !GZ value below which sedimentation is computed

! The following include files are explicitly inserted below, as in the WRF interface
!-- For ARPS:
!------------------------------------------------------------------------------!
!#include "consphy.cdk"
  real,   parameter :: CPD      =.100546e+4            !J K-1 kg-1; specific heat of dry air
  real,   parameter :: CPV      =.186946e+4            !J K-1 kg-1; specific heat of water vapour
  real,   parameter :: RGASD    =.28705e+3             !J K-1 kg-1; gas constant for dry air
  real,   parameter :: RGASV    =.46151e+3             !J K-1 kg-1; gas constant for water vapour
  real,   parameter :: TRPL     =.27316e+3             !K; triple point of water
  real,   parameter :: TCDK     =.27315e+3             !; conversion from kelvin to celsius
  real,   parameter :: RAUW     =.1e+4                 !; density of liquid H2O
  real,   parameter :: EPS1     =.62194800221014       !; RGASD/RGASV
  real,   parameter :: EPS2     =.3780199778986        !; 1 - EPS1
  real,   parameter :: DELTA    =.6077686814144        !; 1/EPS1 - 1
  real,   parameter :: CAPPA    =.28549121795          !; RGASD/CPD
  real,   parameter :: TGL      =.27316e+3             !K; ice temperature in the atmosphere
  real,   parameter :: CONSOL   =.1367e+4              !W m-2; solar constant
  real,   parameter :: GRAV     =.980616e+1            !M s-2; gravitational acceleration
  real,   parameter :: RAYT     =.637122e+7            !M; mean radius of the earth
  real,   parameter :: STEFAN   =.566948e-7            !J m-2 s-1 K-4; Stefan-Boltzmann constant
  real,   parameter :: PI       =.314159265359e+1      !; PI constant = ACOS(-1)
  real,   parameter :: OMEGA    =.7292e-4              !s-1; angular speed of rotation of the earth
  real,   parameter :: KNAMS    =.514791               !; conversion from knots to m/s
  real,   parameter :: STLO     =.6628486583943e-3     !K s2 m-2; Schuman-Newell Lapse Rate
  real,   parameter :: KARMAN   =.35                   !; Von Karman constant
  real,   parameter :: RIC      =.2                    !; Critical Richardson number
  real,   parameter :: CHLC     =.2501e+7              !J kg-1; latent heat of condensation
  real,   parameter :: CHLF     =.334e+6               !J kg-1; latent heat of fusion
  ! DTD: Not sure if the following are needed
  real,   parameter :: T1S      =.27316e+3             !K constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: T2S      =.25816e+3             !K constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: AW       =.3135012829948e+4     ! constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: BW       =.2367075766316e+1     ! constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: AI       =.2864887713087e+4     ! constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: BI       =.166093131502         ! constant used to calculate L/Cp in fcn HTVOCP
  real,   parameter :: SLP      =.6666666666667e-1     ! constant used to calculate L/Cp in fcn HTVOCP
!------------------------------------------------------------------------------!
!#include "dintern.cdk"
      REAL   TTT, PRS, QQQ, EEE, TVI, QST, QQH
      REAL   T00, PR0, TF, PF,FFF , DDFF
      REAL   QSM , DLEMX
      REAL*8 FOEW,FOQST !,FODLE,FOEFQ !,FODQS,FOQFE,FOTVT,FOTTV,FOHR
!      REAL*8 FOLV,FOLS,FOPOIT,FOPOIP,FOTTVH,FOTVHT
      REAL*8 FOEWA,FOQSA ! ,FODLA,FODQA,FOHRA
      REAL*8 FESI !,FDLESI,FESMX,FDLESMX,FQSMX,FDQSMX
  
! DTD: In the WRF interface, many of the following functions are declared somewhat differently
! depending on the value of RWORDSIZE and DWORDSIZE using compiler directives.  I didn't include
! these directives at the present time, because I'm unsure how to proceed for ARPS.
!------------------------------------------------------------------------------!
!#include "fintern.cdk"
!
!   DEFINITION DES FONCTIONS THERMODYNAMIQUES DE BASE
!   POUR LES CONSTANTES, UTILISER LE COMMON /CONSPHY/
!     NOTE: TOUTES LES FONCTIONS TRAVAILLENT AVEC LES UNITES S.I.
!           I.E. TTT EN DEG K, PRS EN PA, QQQ EN KG/KG
!          *** N. BRUNET - MAI 90 ***
!          * REVISION 01 - MAI 94 - N. BRUNET
!                          NOUVELLE VERSION POUR FAIBLES PRESSIONS
!          * REVISION 02 - AOUT 2000 - J-P TOVIESSI
!                          CALCUL EN REAL*8
!          * REVISION 03 - SEPT 2000 - N. BRUNET
!                          AJOUT DE NOUVELLES FONCTIONS
!          * REVISION 04 - JANV 2000 - J. MAILHOT
!                          FONCTIONS EN PHASE MIXTE
!          * REVISION 05 - DEC 2001 - G. LEMAY
!                          DOUBLE PRECISION POUR PHASE MIXTE
!          * REVISION 06 - AVR 2002 - A. PLANTE
!                          AJOUT DES NOUVELLES FONCTIONS FOTTVH ET FOTVHT
!
!     FONCTION DE TENSION DE VAPEUR SATURANTE (TETENS) - EW OU EI SELON TT
      FOEW(TTT) = 610.78D0*DEXP( DMIN1(DSIGN(17.269D0,                     &
       DBLE(TTT)-DBLE(TRPL)),DSIGN                                         &
       (21.875D0,DBLE(TTT)-DBLE(TRPL)))*DABS(DBLE(TTT)-DBLE(TRPL))/        &
       (DBLE(TTT)-35.86D0+DMAX1(0.D0,DSIGN                                 &
       (28.2D0,DBLE(TRPL)-DBLE(TTT)))))
!
!     FONCTION CALCULANT LA DERIVEE SELON T DE  LN EW (OU LN EI)
!      FODLE(TTT)=(4097.93D0+DMAX1(0.D0,DSIGN(1709.88D0,                    &
!       DBLE(TRPL)-DBLE(TTT))))                                             &
!       /((DBLE(TTT)-35.86D0+DMAX1(0.D0,DSIGN(28.2D0,                       &
!       DBLE(TRPL)-DBLE(TTT))))*(DBLE(TTT)-35.86D0+DMAX1(0.D0               &
!       ,DSIGN(28.2D0,DBLE(TRPL)-DBLE(TTT)))))
!
!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT)
      FOQST(TTT,PRS) = DBLE(EPS1)/(DMAX1(1.D0,DBLE(PRS)/FOEW(TTT))-        &
       DBLE(EPS2))
!
!     FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
!      FODQS(QST,TTT)=DBLE(QST)*(1.D0+DBLE(DELTA)*DBLE(QST))*FODLE(TTT)
!     QST EST LA SORTIE DE FOQST
!
!     FONCTION CALCULANT TENSION VAP (EEE) FN DE HUM SP (QQQ) ET PRS
!      FOEFQ(QQQ,PRS) = DMIN1(DBLE(PRS),(DBLE(QQQ)*DBLE(PRS)) /             &
!       (DBLE(EPS1) + DBLE(EPS2)*DBLE(QQQ)))
!
!      FONCTION CALCULANT HUM SP (QQQ) DE TENS. VAP (EEE) ET PRES (PRS)
!      FOQFE(EEE,PRS) = DMIN1(1.D0,DBLE(EPS1)*DBLE(EEE)/(DBLE(PRS)-         &
!       DBLE(EPS2)*DBLE(EEE)))
!
!      FONCTION CALCULANT TEMP VIRT. (TVI) DE TEMP (TTT) ET HUM SP (QQQ)
!      FOTVT(TTT,QQQ) = DBLE(TTT) * (1.0D0 + DBLE(DELTA)*DBLE(QQQ))         

!      FONCTION CALCULANT TEMP VIRT. (TVI) DE TEMP (TTT), HUM SP (QQQ) ET
!      MASSE SP DES HYDROMETEORES.
!      FOTVHT(TTT,QQQ,QQH) = DBLE(TTT) *                                    &
!           (1.0D0 + DBLE(DELTA)*DBLE(QQQ) - DBLE(QQH))
!
!      FONCTION CALCULANT TTT DE TEMP VIRT. (TVI) ET HUM SP (QQQ)
!      FOTTV(TVI,QQQ) = DBLE(TVI) / (1.0D0 + DBLE(DELTA)*DBLE(QQQ))

!      FONCTION CALCULANT TTT DE TEMP VIRT. (TVI), HUM SP (QQQ) ET
!      MASSE SP DES HYDROMETEORES (QQH)
!      FOTTVH(TVI,QQQ,QQH) = DBLE(TVI) /                                    &
!           (1.0D0 + DBLE(DELTA)*DBLE(QQQ) - DBLE(QQH))
!
!      FONCTION CALCULANT HUM REL DE HUM SP (QQQ), TEMP (TTT) ET PRES (PRS)
!      HR = E/ESAT
!       FOHR(QQQ,TTT,PRS) = DMIN1(DBLE(PRS),FOEFQ(QQQ,PRS)) / FOEW(TTT)
!
!     FONCTION CALCULANT LA CHALEUR LATENTE DE CONDENSATION
!      FOLV(TTT) =DBLE(CHLC) - 2317.D0*(DBLE(TTT)-DBLE(TRPL))
!
!     FONCTION CALCULANT LA CHALEUR LATENTE DE SUBLIMATION
!      FOLS(TTT) = DBLE(CHLC)+DBLE(CHLF)+(DBLE(CPV)-                        &
!                  (7.24D0*DBLE(TTT)+128.4D0))*(DBLE(TTT)-DBLE(TRPL))
!
!     FONCTION RESOLVANT L'EQN. DE POISSON POUR LA TEMPERATURE
!     NOTE: SI PF=1000*100, "FOPOIT" DONNE LE THETA STANDARD
!      FOPOIT(T00,PR0,PF)=DBLE(T00)*(DBLE(PR0)/DBLE(PF))**                  &
!                       (-DBLE(CAPPA))
!
!     FONCTION RESOLVANT L'EQN. DE POISSON POUR LA PRESSION
!      FOPOIP(T00,TF,PR0)=DBLE(PR0)*DEXP(-(DLOG(DBLE(T00)/DBLE(TF))/        &
!                       DBLE(CAPPA)))
!
!     LES 5 FONCTIONS SUIVANTES SONT VALIDES DANS LE CONTEXTE OU ON
!     NE DESIRE PAS TENIR COMPTE DE LA PHASE GLACE DANS LES CALCULS
!     DE SATURATION.
!   FONCTION DE VAPEUR SATURANTE (TETENS)
      FOEWA(TTT)=610.78D0*DEXP(17.269D0*(DBLE(TTT)-DBLE(TRPL))/            &
       (DBLE(TTT)-35.86D0))
!   FONCTION CALCULANT LA DERIVEE SELON T DE LN EW
!      FODLA(TTT)=17.269D0*(DBLE(TRPL)-35.86D0)/(DBLE(TTT)-35.86D0)**2
!   FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE
      FOQSA(TTT,PRS)=DBLE(EPS1)/(DMAX1(1.D0,DBLE(PRS)/FOEWA(TTT))-         &
       DBLE(EPS2))
!   FONCTION CALCULANT LA DERIVEE DE QSAT SELON T
!      FODQA(QST,TTT)=DBLE(QST)*(1.D0+DBLE(DELTA)*DBLE(QST))*FODLA(TTT)
!   FONCTION CALCULANT L'HUMIDITE RELATIVE
!      FOHRA(QQQ,TTT,PRS)=DMIN1(DBLE(PRS),FOEFQ(QQQ,PRS))/FOEWA(TTT)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!   Definition of basic thermodynamic functions in mixed-phase mode
!     FFF is the fraction of ice and DDFF its derivative w/r to T
!     NOTE: S.I. units are used
!           i.e. TTT in deg K, PRS in Pa
!          *** J. Mailhot - Jan. 2000 ***
!
!     Saturation calculations in presence of liquid phase only
!     Function for saturation vapor pressure (TETENS)
!      FESI(TTT)=610.78D0*DEXP(21.875D0*(DBLE(TTT)-DBLE(TRPL))/             &
!             (DBLE(TTT)-7.66D0)  )
!      FDLESI(TTT)=21.875D0*(DBLE(TRPL)-7.66D0)/(DBLE(TTT)-7.66D0)**2
!      FESMX(TTT,FFF) = (1.D0-DBLE(FFF))*FOEWA(TTT)+DBLE(FFF)*FESI(TTT)
!      FDLESMX(TTT,FFF,DDFF) = ( (1.D0-DBLE(FFF))*FOEWA(TTT)*FODLA(TTT)     &
!                            + DBLE(FFF)*FESI(TTT)*FDLESI(TTT)              &
!                  + DBLE(DDFF)*(FESI(TTT)-FOEWA(TTT)) )/FESMX(TTT,FFF)
!      FQSMX(TTT,PRS,FFF) = DBLE(EPS1)/                                     &
!              (DMAX1(1.D0,DBLE(PRS)/FESMX(TTT,FFF) ) - DBLE(EPS2)  )
!      FDQSMX(QSM,DLEMX) = DBLE(QSM ) *(1.D0 + DBLE(DELTA)* DBLE(QSM ) )    &
!                           * DBLE(DLEMX )
!
! ! !------------------------------------------------------------------------------!
!***** END of Replace 3 #includes (for ARPS) ***

  real              :: LCP, LFP, LSP, ck5, ck6
  real*8, parameter :: thrd  = 1.d0/3.d0
  real*8, parameter :: sixth = 0.5d0*thrd
  real*8            :: PI2, PIov6, CHLS
  real  , parameter :: thrdSP = 1./3.

  ! Constants used for contact ice nucleation:
  real*8, parameter :: LAMa0  = 6.6d-8     !mean free path at T0 and p0 [W95_eqn58]
  real*8, parameter :: T0     = 293.15d0   !ref. temp.
  real*8, parameter :: p0     = 101325.d0  !ref. pres.
  real*8, parameter :: Ra     = 1.d-6      !aerosol (IN) radius         [M92 p.713; W95_eqn60]
  real*8, parameter :: kBoltz = 1.381d-23  !Boltzmann's constant
  real*8, parameter :: KAPa   = 5.39d5     !aerosol thermal conductivity

 !Testing switches:
  integer, parameter :: airtype     =  2        ! 1 = maritime;  2 = continental
  logical, parameter :: icephase_ON = .true.    !.false. to suppress ice-phase (Part I)
  logical, parameter :: iceDep_ON   = .true.    !.false. to suppress depositional growth of ice
  logical, parameter :: snow_ON     = .true.    !.false. to suppress snow initiation
  logical, parameter :: warm_ON     = .true.    !.false. to suppress warm-phase (Part II)
  logical            :: hl_ON       = .true.    !.false. to suppress hail initiation (independently added by Jason M. and Bryan P.)
  logical, parameter :: autoconv_ON = .true.    ! autoconversion ON/OFF
  logical, parameter :: rainAccr_ON = .true.    ! rain accretion and self-collection ON/OFF
  logical, parameter :: sedi_ON     = .true.    !.false. to suppress sedimentation
  ! Additional testing switches added by DTD, YSJ, and others
  logical, parameter :: SS_ON       = .true.    !.false. to suppress size-sorting (i.e. differential 
                                                ! fall speeds of moments).  Added by DTD.
  logical, parameter :: clipCons_ON = .false.   ! .false. to neglect conservation of mass when 
                                                ! clipping is performed (also neglect latent heat
                                                ! changes due to clipping).  Added by DTD.

  logical :: grpl_ON     = .true.               !.false. to suppress graupel initiation. Added by YSJ 
  logical :: firstcall   = .true.
  
!---Testing --- for sedi subroutine
  real, dimension(ni,nk) :: QX,NX,ZX
  real*8  :: DxMAX,afx,bfx,cmx,dmx,ckQx1,ckQx2,ckQx3
  real    :: cmxSP,dmxSP,dtx,cx6
  integer :: npassx
!--------------

!__________________________________________________________________________!
! Scheme version:    (all categories)

! integer, parameter :: scheme = 1  ![SM]  - single-moment
! integer, parameter :: scheme = 2  ![DM1] - double-moment, fixed-ALPHAx
! integer, parameter :: scheme = 3  ![DM2] - double-moment, diagnosed-ALPHAx
! integer, parameter :: scheme = 4  ![TM]  - triple-moment
!__________________________________________________________________________!

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ! Decision for graupel and hail initiation
!  IF (firstcall) THEN
!    IF (graupel_ON == 0) THEN
!      grpl_ON = .false.
!    ENDIF
!    IF (hail_ON == 0) THEN
!      hl_ON = .false.
!    ENDIF
!    firstcall = .false.
!  ENDIF

  ! DTD: Initialize values that were originally parameters
  dei = DBLE(rhoice)
  des = DBLE(rhosnow)
  deg = DBLE(rhogrpl)
  deh = DBLE(rhohail)
  Ncfix = DBLE(ntcloud)
  Norfix = DBLE(n0rain)
  Nosfix = DBLE(n0snow)
  Nogfix = DBLE(n0grpl)
  Nohfix = DBLE(n0hail)
  ALFrfix = DBLE(alpharain)
  ALFifix = DBLE(alphaice)
  ALFsfix = DBLE(alphasnow)
  ALFgfix = DBLE(alphagrpl)
  ALFhfix = DBLE(alphahail)
  
  !==================================================================================!

  !----------------------------------------------------------------------------------!
  !                      PART 1:   Prelimiary Calculations                           !
  !----------------------------------------------------------------------------------!

  if (scheme<1 .or. scheme>4) then
    print*, '**************************************************'
    print*, '*                                                *'
    print*, '*            ABORT in S/R MY_MAIN_TM             *'
    print*, '*                                                *'
    print*, '*   INCORRECT SPECIFICATION OF MY_TMOM_VERSION   *'
    print*, '*             (must be 1, 2, 3, or 4)            *'
    print*, '*                                                *'
    print*, '**************************************************'
    !stop
    !CALL arpsstop('INCORRECT SPECIFICATION OF SCHEME VERSION IN MULTIMOMENT.',1)
  endif

! DTD/WYH: Initialize pointers

  IF (scheme >= 1) THEN
    QC => qnz(:,:,1)
    QR => qnz(:,:,2)
    QI => qnz(:,:,3)
    QN => qnz(:,:,4)
    QG => qnz(:,:,5)
    QH => qnz(:,:,6)

    QCM => qnzm(:,:,1)
    QRM => qnzm(:,:,2)
    QIM => qnzm(:,:,3)
    QNM => qnzm(:,:,4)
    QGM => qnzm(:,:,5)
    QHM => qnzm(:,:,6)
  END IF

  IF (scheme == 1) THEN

   nctem(:,:,:) = 0.0

    NC => nctem(:,:,1)
    NR => nctem(:,:,2)
    NY => nctem(:,:,3)
    NN => nctem(:,:,4)
    NG => nctem(:,:,5)
    NH => nctem(:,:,6)

    NCM => nctem(:,:,7)
    NRM => nctem(:,:,8)
    NYM => nctem(:,:,9)
    NNM => nctem(:,:,10)
    NGM => nctem(:,:,11)
    NHM => nctem(:,:,12)

  ELSE IF (scheme >= 2) THEN
    NC => qnz(:,:,7)
    NR => qnz(:,:,8)
    NY => qnz(:,:,9)
    NN => qnz(:,:,10)
    NG => qnz(:,:,11)
    NH => qnz(:,:,12)

    NCM => qnzm(:,:,7)
    NRM => qnzm(:,:,8)
    NYM => qnzm(:,:,9)
    NNM => qnzm(:,:,10)
    NGM => qnzm(:,:,11)
    NHM => qnzm(:,:,12)
  END IF

  IF (scheme == 4) THEN
    ZR => qnz(:,:,13)
    ZI => qnz(:,:,14)
    ZN => qnz(:,:,15)
    ZG => qnz(:,:,16)
    ZH => qnz(:,:,17)

    ZRM => qnzm(:,:,13)
    ZIM => qnzm(:,:,14)
    ZNM => qnzm(:,:,15)
    ZGM => qnzm(:,:,16)
    ZHM => qnzm(:,:,17)
  END IF

!-------------
!Convert N from #/kg to #/m3:
! DTD: This is now done in model dynamics, so is commented out here.
!  IF (scheme >= 2) THEN
!  do k= 1,nk-1
!    do i= 1,ixe
!      !tmp1= S(i,k)*PSM(i)/(RGASD*TM(i,k))  !air density at time (t-1)
!      !tmp2= S(i,k)*PS(i)/(RGASD*T(i,k))    !air density at time (*)!
!
!      tmp1= PM(i,k)/(RGASD*TM(i,k))  !air density at time (t-1)
!      tmp2= P(i,k)/(RGASD*T(i,k))    !air density at time (*)
!
!      NCM(i,k)= NCM(i,k)*tmp1;   NC(i,k)= NC(i,k)*tmp2
!      NRM(i,k)= NRM(i,k)*tmp1;   NR(i,k)= NR(i,k)*tmp2
!      NYM(i,k)= NYM(i,k)*tmp1;   NY(i,k)= NY(i,k)*tmp2
!      NNM(i,k)= NNM(i,k)*tmp1;   NN(i,k)= NN(i,k)*tmp2
!      NGM(i,k)= NGM(i,k)*tmp1;   NG(i,k)= NG(i,k)*tmp2
!      NHM(i,k)= NHM(i,k)*tmp1;   NH(i,k)= NH(i,k)*tmp2
!      IF (scheme == 4) THEN
!        ZRM(i,k)= ZRM(i,k)*tmp1;   ZR(i,k)= ZR(i,k)*tmp2
!        ZIM(i,k)= ZIM(i,k)*tmp1;   ZI(i,k)= ZI(i,k)*tmp2
!        ZNM(i,k)= ZNM(i,k)*tmp1;   ZN(i,k)= ZN(i,k)*tmp2
!        ZGM(i,k)= ZGM(i,k)*tmp1;   ZG(i,k)= ZG(i,k)*tmp2
!        ZHM(i,k)= ZHM(i,k)*tmp1;   ZH(i,k)= ZH(i,k)*tmp2
!      END IF
!    enddo
!  enddo
!  END IF
!=============

 ! The SSxx arrays are for passed to the volatile bus for output as 3-D diagnostic
 ! output variables, for testing purposes.  For example, to output the
 ! instantanous value of the deposition rate, add 'SS01(i,k) = QVDvi'  in the
 ! appropriate place.  It can then be output as a 3-D physics variable by adding
 ! it to the sortie_p list in 'outcfgs.out'

! DTD: turned off for now.  Can probably reinstate these later to allow for output of
! process terms as in earlier version of the scheme interface.

!  SS01= 0.; SS02= 0.; SS03= 0.; SS04= 0.; SS05= 0.; SS06= 0.; SS07= 0.; SS08= 0.
!  SS09= 0.; SS10= 0.; SS11= 0.; SS12= 0.; SS13= 0.; SS14= 0.; SS15= 0.; SS16= 0.
!  SS17= 0.; SS18= 0.; SS19= 0.; SS20= 0.

  PI2   = PI*2.d0
  PIov6 = PI*sixth
  CHLS  = CHLC+CHLF  !J k-1; latent heat of sublimation
  LCP   = CHLC/CPD
  LFP   = CHLF/CPD
  LSP   = LCP+LFP
  ck5   = 4098.170*LCP
  ck6   = 5806.485*LSP
  dt    = dble(DT_sp) ! DTD: Note! No longer multiplied by 2 here!
  idt   = 1./DT_sp

 !-------------------------------------------------------------------!
 ! Defining constants based on size distribution parameters:

  ! Cloud:
  iMUc= 1.d0/MUc
  GC1=  gammaDP(ALFc+1.0d0)            !i.e. gammaDP(alf + 1)
  GC2=  gammaDP(ALFc+1.d0+3.0d0*iMUc)  !i.e. gammaDP(alf + 4)
  GC3=  gammaDP(ALFc+1.d0+6.0d0*iMUc)  !i.e. gammaDP(alf + 7)
  GC4=  gammaDP(ALFc+1.d0+9.0d0*iMUc)  !i.e. gammaDP(alf + 10)

  GC11= gammaDP(1.0d0*iMUc+1.0d0+ALFc)
  GC12= gammaDP(2.0d0*iMUc+1.0d0+ALFc)
  GC5=  gammaDP(1.0d0+ALFc);             GC13= gammaDP(3.0d0*iMUc+1.0d0+ALFc)
  GC6=  gammaDP(1.0d0+ALFc+1.0d0*iMUc);  GC14= gammaDP(4.0d0*iMUc+1.0d0+ALFc)
  GC7=  gammaDP(1.0d0+ALFc+2.0d0*iMUc);  GC15= gammaDP(5.0d0*iMUc+1.0d0+ALFc)
  GC8=  gammaDP(1.0d0+ALFc+3.0d0*iMUc)
  cexc9= GC2/GC1*PIov6*dew

  ! Mass parameters [ m(D) = cD^d ]
  cmr=  PIov6*dew;  cmrSP= sngl(cmr);  icmr= 1.d0/cmr
  cmi=  440.0d0;    cmiSP= sngl(cmi);  icmi= 1.d0/cmi
  cms=  PIov6*des;  cmsSP= sngl(cms);  icms= 1.d0/cms
  cmg=  PIov6*deg;  cmgSP= sngl(cmg);  icmg= 1.d0/cmg
  cmh=  PIov6*deh;  cmhSP= sngl(cmh);  icmh= 1.d0/cmh

  ! g(alpha) function values for specified integer values of ALPHA:
  ! where g(a)= [(6+a)(5+a)(4+a)]\[(3+a)(2+a)(1+a)]
  GalphaRaut = ((6.d0+ALFrAUT)*(5.d0+ALFrAUT)*(4.d0+ALFrAUT))/  &
               ((3.d0+ALFrAUT)*(2.d0+ALFrAUT)*(1.d0+ALFrAUT))
  Galpha0= 20.000d0;   Galpha2= 5.600d0;   Galpha4= 3.429d0
  Galpha1=  8.750d0;   Galpha3= 4.200d0;   Galpha5= 2.946d0

!=======================================================================================!

! Compute variables for BLG sedimentation:

 !Determine the upper-most level in each column up to which to compute sedimentation:
 !DTD: Though this is still calculated, I changed the sedimentation subroutines that use this
 !variable (ktop_sedi) to just use the top of the model index for now (i.e. sedimentation is done for entire
 !column), but this can be changed back in the future.
 
  ktop_sedi= 0
  tmp = 0
  do i= 1,ixe
     do k=nk-1,1,-1         ! Changed by DTD (reversed order of loop since indices increase as you go
                            ! up in ARPS)
       ktop_sedi(i)= k
       !print*,'i,ktop_sedi = ',i,ktop_sedi(i)
       !if (GZ(i,k)<gzMax_sedi) exit
       ! DTD (07/18/2011).  Use tmp for GZ here
       tmp = ZP(i,k)*GRAV   ! Need to check if this is correct: depends on whether GZ is supposed
                            ! to be geopotential height, or is actually just the geopotential
       if (tmp < gzMax_sedi) exit
     enddo
  enddo

  !Compute thickness of layers for sedimentation calcuation:
!  do k=2,nk
!     DZ(:,k)= (GZ(:,k-1)-GZ(:,k))/GRAV
!  enddo
!  DZ(:,1)= DZ(:,2)

! DTD (originally WYH): No geopotential conversion needed here (also array index increases
! from lower to upper levels, rather than the reverse in the original code).
  DO k=1,nk-1
    DZ(:,k) = ZP(:,k+1)-ZP(:,k)
  END DO
  DZ(:,nk) = DZ(:,nk-1)

!Max. Courant number for BLG sedimentation (used to determine npassx and dtx):
  CoMAXr = 3.;  CoMAXi = 5.;  CoMAXs = 5.;  CoMAXg = 5.;  CoMAXh = 6.
  !Note: DZ(1,nk) is the grid-spacing between the lowest 2 model levels.  The
  !      amount of time-splitting for sedimentation (i.e. the number of times
  !      BLG is called each model time step [npassx]) is estimated from this
  !      min DZ, the maximum possible fall speed (VxMAX), and the specified
  !      max Courant number (which may be greater than 1.0 for BLG scheme).

  !Note: The selection of CoMAXx was deemed appropriate for the IMPROVE2 simulations.
  !      However, these runs had only small hail and little snow in PBL.  Thus,
  !      Thus, the approprite values of CoMAXx may be different for other types
  !      of cases.


!  npassr= max(1, nint( DT*Vrmax/(CoMAXr*DZ(1,nk)) ))
!  npassi= max(1, nint( DT*Vimax/(CoMAXi*DZ(1,nk)) ))
!  npasss= max(1, nint( DT*Vsmax/(CoMAXs*DZ(1,nk)) ))
!  npassg= max(1, nint( DT*Vgmax/(CoMAXg*DZ(1,nk)) ))
!  npassh= max(1, nint( DT*Vhmax/(CoMAXh*DZ(1,nk)) ))

  npassr= max(1, nint( DT*Vrmax/(CoMAXr*DZ(1,2)) ))
  npassi= max(1, nint( DT*Vimax/(CoMAXi*DZ(1,2)) ))
  npasss= max(1, nint( DT*Vsmax/(CoMAXs*DZ(1,2)) ))
  npassg= max(1, nint( DT*Vgmax/(CoMAXg*DZ(1,2)) ))
  npassh= max(1, nint( DT*Vhmax/(CoMAXh*DZ(1,2)) ))

  dtr = DT/float(npassr);   cr6= GRAV*dtr
  dti = DT/float(npassi);   ci6= GRAV*dti
  dts = DT/float(npasss);   cs6= GRAV*dts
  dtg = DT/float(npassg);   cg6= GRAV*dtg
  dth = DT/float(npassh);   ch6= GRAV*dth
  ck7 = 1./(DT*GRAV)

!=======================================================================================!

! DTD/WYH: We don't need these tendencies for the ARPS model; the scheme passes out the 
! adjusted microphysics arrays (at the future time level) directly.

! Temporarily store arrays at time (t*) in order to compute (at the end of subroutine)
! the final VXTEND as VXTEND = ( VX{t+1} - VX{t*} )/dt :
!  T_TEND = T ;  Q_TEND = Q
!  QCTEND = QC;  QRTEND = QR;  QITEND = QI;  QNTEND = QN;  QGTEND = QG;  QHTEND = QH
!  if (scheme > 1) then  ![DM and TM]
!     NCTEND = NC; NRTEND = NR; NYTEND = NY;  NNTEND = NN; NGTEND = NG;  NHTEND = NH
!  endif
!  if (scheme==4) then   ![TM]
!     ZRTEND = ZR;   ZITEND = ZI;  ZNTEND = ZN;  ZGTEND = ZG;  ZHTEND = ZH
!  endif


! Clip all moments to zero if one or more corresponding category moments are less than
!  the minimum allowable value:
! (Note: Clipped mass is added back to water vapor field to conserve total mass and
!        temperature is modified to account for latent heating/cooling from phase change)
! DTD: Added a testing flag to turn off the adding of mass back to vapor and associated
! latent heat change.  Recommended set to false for the ARPS (i.e. do *not* conserve total water or
! adjust latent heat)

  do k= 1,nk-1
     do i= 1,ixe
     
      IF (scheme== 1) THEN   ! prevent negative mass being passed in, by WYH

        if(QC(i,k)<epsQ) QC(i,k)= 0.
        if(QR(i,k)<epsQ) QR(i,k)= 0.
        if(QI(i,k)<epsQ) QI(i,k)= 0.
        if(QN(i,k)<epsQ) QN(i,k)= 0.
        if(QG(i,k)<epsQ) QG(i,k)= 0.
        if(QH(i,k)<epsQ) QH(i,k)= 0.

        if(QCM(i,k)<epsQ) QCM(i,k)= 0.
        if(QRM(i,k)<epsQ) QRM(i,k)= 0.
        if(QIM(i,k)<epsQ) QIM(i,k)= 0.
        if(QNM(i,k)<epsQ) QNM(i,k)= 0.
        if(QGM(i,k)<epsQ) QGM(i,k)= 0.
        if(QHM(i,k)<epsQ) QHM(i,k)= 0.

      ELSE IF (scheme==2 .or. scheme==3) THEN

        if(QC(i,k)<epsQ .or. NC(i,k)<epsN)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QC(i,k)
           T(i,k) = T(i,k) - LCP*QC(i,k)
          endif
           QC(i,k)= 0.;   NC(i,k)= 0.
        endif
        if (QR(i,k)<epsQ .or. NR(i,k)<epsN)   then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QR(i,k)
           T(i,k) = T(i,k) - LCP*QR(i,k)
           endif
           QR(i,k)= 0.;  NR(i,k)= 0.
        endif
        if (QI(i,k)<epsQ .or. NY(i,k)<epsN)   then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QI(i,k)
           T(i,k) = T(i,k) - LSP*QI(i,k)
          endif
           QI(i,k)= 0.;  NY(i,k)= 0.
        endif
        if (QN(i,k)<epsQ .or. NN(i,k)<epsN)   then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QN(i,k)
           T(i,k) = T(i,k) - LSP*QN(i,k)
          endif
           QN(i,k)= 0.;  NN(i,k)= 0.
        endif
        if (QG(i,k)<epsQ .or. NG(i,k)<epsN)   then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QG(i,k)
           T(i,k) = T(i,k) - LSP*QG(i,k)
          endif
           QG(i,k)= 0.;  NG(i,k)= 0.
        endif
        if (QH(i,k)<epsQ .or. NH(i,k)<epsN)   then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QH(i,k)
           T(i,k) = T(i,k) - LSP*QH(i,k)
          endif
           QH(i,k)= 0.;  NH(i,k)= 0.
        endif

        if(QCM(i,k)<epsQ .or. NCM(i,k)<epsN)  then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QCM(i,k)
           TM(i,k) = TM(i,k) - LCP*QCM(i,k)
          endif
           QCM(i,k)= 0.;  NCM(i,k)= 0.
        endif
        if (QRM(i,k)<epsQ .or. NRM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QRM(i,k)
           TM(i,k) = TM(i,k) - LCP*QRM(i,k)
          endif
           QRM(i,k)= 0.;  NRM(i,k)= 0.
        endif
        if (QIM(i,k)<epsQ .or. NYM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QIM(i,k)
           TM(i,k) = TM(i,k) - LSP*QIM(i,k)
          endif
           QIM(i,k)= 0.;  NYM(i,k)= 0.
        endif
        if (QNM(i,k)<epsQ .or. NNM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QNM(i,k)
           TM(i,k) = TM(i,k) - LSP*QNM(i,k)
          endif
           QNM(i,k)= 0.;  NNM(i,k)= 0.
        endif
        if (QGM(i,k)<epsQ .or. NGM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QGM(i,k)
           TM(i,k) = TM(i,k) - LSP*QGM(i,k)
          endif
           QGM(i,k)= 0.;  NGM(i,k)= 0.
        endif
        if (QHM(i,k)<epsQ .or. NHM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QHM(i,k)
           TM(i,k) = TM(i,k) - LSP*QHM(i,k)
          endif
           QHM(i,k)= 0.;  NHM(i,k)= 0.
        endif

      ELSEIF (scheme==4) THEN

        if(QC(i,k)<epsQ .or. NC(i,k)<epsN) then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QC(i,k)
           T(i,k) = T(i,k) - LCP*QC(i,k)
          endif
           QC(i,k)= 0.;   NC(i,k)= 0.
        endif
        if (QR(i,k)<epsQ .or. NR(i,k)<epsN .or. ZR(i,k)<epsZ)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QR(i,k)
           T(i,k) = T(i,k) - LCP*QR(i,k)
          endif
           QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
        endif
        if (QI(i,k)<epsQ .or. NY(i,k)<epsN .or. ZI(i,k)<epsZ)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QI(i,k)
           T(i,k) = T(i,k) - LSP*QI(i,k)
          endif
           QI(i,k)= 0.;  NY(i,k)= 0.;  ZI(i,k)= 0.
        endif
        if (QN(i,k)<epsQ .or. NN(i,k)<epsN .or. ZN(i,k)<epsZ)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QN(i,k)
           T(i,k) = T(i,k) - LSP*QN(i,k)
          endif
           QN(i,k)= 0.;  NN(i,k)= 0.;  ZN(i,k)= 0.
        endif
        if (QG(i,k)<epsQ .or. NG(i,k)<epsN .or. ZG(i,k)<epsZ)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QG(i,k)
           T(i,k) = T(i,k) - LSP*QG(i,k)
          endif
           QG(i,k)= 0.;  NG(i,k)= 0.;  ZG(i,k)= 0.
        endif
        if (QH(i,k)<epsQ .or. NH(i,k)<epsN .or. ZH(i,k)<epsZ)    then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QH(i,k)
           T(i,k) = T(i,k) - LSP*QH(i,k)
          endif
           QH(i,k)= 0.;  NH(i,k)= 0.;  ZH(i,k)= 0.
        endif

        if(QCM(i,k)<epsQ .or. NCM(i,k)<epsN) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QCM(i,k)
           TM(i,k) = TM(i,k) - LCP*QCM(i,k)
          endif
           QCM(i,k)= 0.;  NCM(i,k)= 0.
        endif
        if (QRM(i,k)<epsQ .or. NRM(i,k)<epsN .or. ZRM(i,k)<epsZ) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QRM(i,k)
           TM(i,k) = TM(i,k) - LCP*QRM(i,k)
          endif
           QRM(i,k)= 0.;  NRM(i,k)= 0.;  ZRM(i,k)= 0.
        endif
        if (QIM(i,k)<epsQ .or. NYM(i,k)<epsN .or. ZIM(i,k)<epsZ) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QIM(i,k)
           TM(i,k) = TM(i,k) - LSP*QIM(i,k)
          endif
           QIM(i,k)= 0.;  NYM(i,k)= 0.;  ZIM(i,k)= 0.
        endif
        if (QNM(i,k)<epsQ .or. NNM(i,k)<epsN .or. ZNM(i,k)<epsZ) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QNM(i,k)
           TM(i,k) = TM(i,k) - LSP*QNM(i,k)
          endif
           QNM(i,k)= 0.;  NNM(i,k)= 0.;  ZNM(i,k)= 0.
        endif
        if (QGM(i,k)<epsQ .or. NGM(i,k)<epsN .or. ZGM(i,k)<epsZ) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QGM(i,k)
           TM(i,k) = TM(i,k) - LSP*QGM(i,k)
          endif
           QGM(i,k)= 0.;  NGM(i,k)= 0.;  ZGM(i,k)= 0.
        endif
        if (QHM(i,k)<epsQ .or. NHM(i,k)<epsN .or. ZHM(i,k)<epsZ) then
          if(clipCons_ON) then
           QM(i,k) = QM(i,k) + QHM(i,k)
           TM(i,k) = TM(i,k) - LSP*QHM(i,k)
          endif
           QHM(i,k)= 0.;  NHM(i,k)= 0.;  ZHM(i,k)= 0.
        endif

      ENDIF

    enddo  !i-loop
  enddo    !k-loop;    (clipping)

  QM = max(QM,0.)
  Q  = max(Q ,0.)

! DTD/WYH: Don't need the following approximate values for ARPS

! Approximate values at time {t}:
!  [ ave. of values at {*} (advected, but no physics tendency added) and {t-dt} ]:
!  HPS= 0.5*(PSM+PS);   TM = 0.5*(TM + T);   QM = 0.5*(QM + Q)
!  QCM= 0.5*(QCM+QC);   QRM= 0.5*(QRM+QR);   QIM= 0.5*(QIM+QI)
!  QNM= 0.5*(QNM+QN);   QGM= 0.5*(QGM+QG);   QHM= 0.5*(QHM+QH)
!  if (scheme>1)  then  ![DM] or [TM]
!     NCM= 0.5*(NCM+NC);   NRM= 0.5*(NRM+NR);   NYM= 0.5*(NYM+NY)
!     NNM= 0.5*(NNM+NN);   NGM= 0.5*(NGM+NG);   NHM= 0.5*(NHM+NH)
!  endif
!  if (scheme==4) then  ![TM]
!     ZRM= 0.5*(ZRM+ZR);   ZIM= 0.5*(ZIM+ZI);   ZNM= 0.5*(ZNM+ZN)
!     ZGM= 0.5*(ZGM+ZG);   ZHM= 0.5*(ZHM+ZH)
!  endif

  ! DTD/WYH: column pressures passed in directly, so don't use S array (sigma)

  do k=1,nk-1
     do i=1,ixe
!    ! Saturation mixing ratios:  [used in both Parts I and II]
!        QSW(i,k)= FOQSA(TM(i,k),HPS(i)*S(i,k))      !wrt. liquid water at (t)
!        QSS(i,k)= FOQST( T(i,k), PS(i)*S(i,k))      !wrt. ice surface  at (*)
!        QSI(i,k)= FOQST(TM(i,k),HPS(i)*S(i,k))      !wrt. ice surface  at (t)
!    ! Air density at time (t)
!        DE(i,k)= S(i,k)*HPS(i)/(RGASD*TM(i,k))      !air density at time  (t)
     ! Saturation mixing ratios:  [used in both Parts I and II]
      QSW(i,k)= FOQSA(TM(i,k),PM(i,k))      !wrt. liquid water at (t)
      QSS(i,k)= FOQST( T(i,k), P(i,k))      !wrt. ice surface  at (*)
      QSI(i,k)= FOQST(TM(i,k),PM(i,k))      !wrt. ice surface  at (t)
     ! Air density at time (t)
      DE(i,k)= PM(i,k)/(RGASD*TM(i,k))      !air density at time  (t)
     enddo
  enddo


  do i= 1,ixe
! Air-density factor: (for fall velocity computations)
     !DEo= dble(DE(i,nk))
     DEo= dble(DE(i,2))
     !gamfact(i,:)  = sqrt(DEo/dble(DE(i,:)))
    do k = 1,nk-1    ! Added by WYH.
      gamfact(i,k)  = sqrt(DEo/dble(DE(i,k)))
    enddo
  ! DTD/WYH: Don't need to the following for ARPS; it's already w on height levels
! Convert 'omega' (on thermodynamic levels) to 'w' (on momentum):
!     do k= 2,nk-1
!        WW(i,k)= -0.5/(DE(i,k)*GRAV)*(Womega(i,k-1)+Womega(i,k+1))
!     enddo
!     WW(i,1) = -0.5/(DE(i,1)*GRAV)*Womega(i,1)
!     WW(i,nk)= -0.5/(DE(i,nk)*GRAV)*Womega(i,nk)
  enddo

  !----------------------------------------------------------------------------------!
  !                 End of Preliminary Calculation section (Part 1)                  !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                      PART 2: Cold Microphysics Processes                         !
  !----------------------------------------------------------------------------------!

! Determine the active grid points (i.e. those which scheme should treat):
  activePoint = .false.
  DO k=1,nk-1
     DO i=1,ixe
        log1= ((QIM(i,k)+QGM(i,k)+QNM(i,k)+QHM(i,k))<epsQ)  !no solid  (i,g,s,h)
        log2= ((QCM(i,k)+QRM(i,k))                  <epsQ)  !no liquid (c,r)
        log3= ((TM(i,k)>TRPL) .and. log1)                   !T>0C & no i,g,s,h
        log4= log1.and.log2.and.(QM(i,k)<QSI(i,k))          !no sol. or liq.; subsat(i)
        if (.not.( log3 .or. log4 ) .and. icephase_ON) then
          activePoint(i,k)= .true.
        endif
     ENDDO
  ENDDO

    ! Size distribution parameters:
    !  Note: + 'thrd' should actually be '1/dmx'(but dmx=3 for all categories x)
    !        + If Qx=0, LAMx etc. are never be used in any calculations
    !          (If Qc=0, CLcy etc. will never be calculated. iLAMx is set to 0
    !           to avoid possible compiler problems.)

  DO k= 1,nk-1
    DO i= 1,ixe
      IF (activePoint(i,k)) THEN

    ! Cloud:
       if (QCM(i,k)>epsQ) then
          if (scheme==1) NCM(i,k)= Ncfix
          Dc(i,k)     = (dble(DE(i,k)*QCM(i,k)/NCM(i,k))*icmr)**thrd
          iLAMc(i,k)  = ((dble(DE(i,k)*QCM(i,k)/NCM(i,k)))/cexc9)**thrd
          iLAMc2(i,k) = iLAMc(i,k) *iLAMc(i,k)
          iLAMc3(i,k) = iLAMc2(i,k)*iLAMc(i,k)
          iLAMc4(i,k) = iLAMc2(i,k)*iLAMc2(i,k)
          iLAMc5(i,k) = iLAMc3(i,k)*iLAMc2(i,k)
       else
          Dc(i,k)     = 0.d0;   iLAMc3(i,k)= 0.d0
          iLAMc(i,k)  = 0.d0;   iLAMc4(i,k)= 0.d0
          iLAMc2(i,k) = 0.d0;   iLAMc5(i,k)= 0.d0
       endif

    ! Rain:
       if (QRM(i,k)>epsQ) then
          if      (scheme==1) then
             ALFr= ALFrfix
             tmpdp1   = gammaDP(1.d0+ALFr)
             tmpdp2   = gammaDP(4.d0+ALFr)
             NRM(i,k) = (Norfix*tmpdp1)**(3./(4.+ALFr))*(tmpdp1/tmpdp2*DE(i,k)*   &
                        QRM(i,k)/cmr)**((1.+ALFr)/(4.+ALFr))  !i.e. NRM = f(No,QRM)
          endif
          Dr(i,k)    = (dble(DE(i,k)*QRM(i,k)/NRM(i,k))*icmr)**thrd
          if      (scheme==2) then
             ALFr= ALFrfix
          else if (scheme==3) then
             ALFr= diagAlpha_v33(Dr(i,k),1)
          else if (scheme==4) then
             ALFr= max(ALFrMIN, solveAlpha_v33(QRM(i,k),NRM(i,k),ZRM(i,k),cmrSP,DE(i,k)) )
          endif
          cexr1       = 1.d0+ALFr+dmr+bfr
          cexr2       = 1.d0+ALFr+dmr
          ckQr1       = afr*gammaDP(1.d0+ALFr+dmr+bfr)/gammaDP(1.d0+ALFr+dmr)
          GR5(i,k)    = gammaDP(1.d0+ALFr);   iGR5(i,k)= 1.d0/GR5(i,k)
          GR31(i,k)   = gammaDP(2.d0+ALFr)
          GR32(i,k)   = gammaDP(3.d0+ALFr)
          GR33(i,k)   = gammaDP(4.d0+ALFr)
          GR34(i,k)   = gammaDP(5.d0+ALFr)
          GR35(i,k)   = gammaDP(6.d0+ALFr)
          cexr9       = cmr*GR33(i,k)*iGR5(i,k)
          GalphaR(i,k)= ((6.d0+ALFr)*(5.d0+ALFr)*(4.d0+ALFr))/                           &
                        ((3.d0+ALFr)*(2.d0+ALFr)*(1.d0+ALFr))
          iLAMr(i,k)  = max( (dble(DE(i,k)*QRM(i,k)/NRM(i,k))/cexr9)**thrd, iLAMmin1 )
          iLAMr2(i,k) = iLAMr(i,k) *iLAMr(i,k)
          iLAMr3(i,k) = iLAMr2(i,k)*iLAMr(i,k)
          iLAMr4(i,k) = iLAMr2(i,k)*iLAMr2(i,k)
          iLAMr5(i,k) = iLAMr3(i,k)*iLAMr2(i,k)
          if (Dr(i,k)>40.d-6) then
             vr0(i,k) = dble(gamfact(i,k))*ckQr1*(1.d0/iLAMr(i,k))**cexr2/(1.d0/iLAMr(i,k)  &
                            +ffr)**cexr1
          else
             vr0(i,k) = 0.d0
          endif
       else
          iLAMr(i,k)  = 0.d0;  Dr(i,k)    = 0.d0;  vr0(i,k)   = 0.d0;  GalphaR(i,k)= 0.d0
          iLAMr2(i,k) = 0.d0;  iLAMr3(i,k)= 0.d0;  iLAMr4(i,k)= 0.d0;  iLAMr5(i,k) = 0.d0
          !either initialize GR5(i,k) etc. here or do not initiaize any (test)
       endif

    ! Ice:
       if (QIM(i,k)>epsQ) then
          if      (scheme==1) then
             ALFi(i,k)= ALFifix
             NYM(i,k) = 5.*exp(0.304*(TRPL-max(233.,TM(i,k))))  !Cooper eqn.
          else if (scheme==2) then
             ALFi(i,k)= ALFifix
          else if (scheme==3) then
             Di   = (dble(DE(i,k)*QIM(i,k)/NYM(i,k))*icmi)**thrd
             ALFi(i,k)= diagAlpha_v33(Di,2)
          else if (scheme==4) then
             ALFi(i,k)= min( ALFiMAX,                                                    &
                             solveAlpha_v33(QIM(i,k),NYM(i,k),ZIM(i,k),cmiSP,DE(i,k)) )
          endif
          GI4        = gammaDP(ALFi(i,k)+dmi+bfi)
          GI7        = gammaDP(1.d0+ALFi(i,k)+bfi)
          GI8        = gammaDP(dmi+1.d0+ALFi(i,k))
          ckQi1      = afi*gammaDP(1.d0+ALFi(i,k)+dmi+bfi)/gammaDP(1.d0+ALFi(i,k)+dmi)
          GI2(i,k)    = gammaDP(1.d0+ALFi(i,k))
          GI31(i,k)   = gammaDP(2.d0+ALFi(i,k))
          GI32(i,k)   = gammaDP(3.d0+ALFi(i,k))
          cexi9      = cmi*gammaDP(1.d0+ALFi(i,k)+dmi)/GI2(i,k)
          iLAMi(i,k)  = max( iLAMmin2, (dble(DE(i,k)*QIM(i,k)/NYM(i,k))/cexi9)**thrd )
          iLAMi2(i,k) = iLAMi(i,k) *iLAMi(i,k)
          iLAMi3(i,k) = iLAMi2(i,k)*iLAMi(i,k)
          iLAMi4(i,k) = iLAMi2(i,k)*iLAMi2(i,k)
          iLAMi5(i,k) = iLAMi3(i,k)*iLAMi2(i,k)
        !Note: It may be better to only have iLAM0..iLAM3 and just use iLAM3*iLAM2 in place
        !      of iLAM5, which appears in four places; this would save some arrays..
        !      The only advantage of having iLAMi5 is to reduce 4 multiplication lines and
        !      to make the code slightly more readable (neither of which is of primary importance)
        !      An array "iNYM" = 1/NYM may be more practical, since  /NYM appears often.  It is
        !      better to avoid a division than a couple of multiplications.
          iLAMiB0(i,k)= iLAMi(i,k)**(bfi)
          iLAMiB1(i,k)= iLAMi(i,k)**(bfi+1.d0)
          iLAMiB2(i,k)= iLAMi(i,k)**(bfi+2.d0)
          vi0(i,k)    = dble(gamfact(i,k))*ckQi1*iLAMiB0(i,k)
       else
          ALFi(i,k)   = 0.d0;  iLAMi(i,k)  = 0.d0;  vi0(i,k)    = 0.d0
          iLAMi2(i,k) = 0.d0;  iLAMi3(i,k) = 0.d0;  iLAMi4(i,k) = 0.d0;  iLAMi5(i,k)= 0.d0
          iLAMiB0(i,k)= 0.d0;  iLAMiB1(i,k)= 0.d0;  iLAMiB2(i,k)= 0.d0
          NYM(i,k)    = 0.
       endif

    ! Snow:
       if (QNM(i,k)>epsQ) then
          if      (scheme==1) then
             ALFs(i,k)= ALFsfix
             tmpdp1   = gammaDP(1.d0+ALFs(i,k))
             tmpdp2   = gammaDP(4.d0+ALFs(i,k))
             NNM(i,k) = (Nosfix*tmpdp1)**(3./(4.+ALFs(i,k)))*(tmpdp1/tmpdp2*DE(i,k)*     &
                         QNM(i,k)/cms)**((1.+ALFs(i,k))/(4.+ALFs(i,k)))  !i.e. NXM = f(No,QXM)

          else if (scheme==2) then
             ALFs(i,k)= ALFsfix
          else if (scheme==3) then
             Ds   = (dble(DE(i,k)*QNM(i,k)/NNM(i,k))*icms)**thrd
             ALFs(i,k)= diagAlpha_v33(Ds,3)
          else if (scheme==4) then
             ALFs(i,k)= min( ALFsMAX,                                                    &
                             solveAlpha_v33(QNM(i,k),NNM(i,k),ZNM(i,k),cmsSP,DE(i,k)) )
          endif
          GS2(i,k)    = gammaDP(1.d0+ALFs(i,k))
          GS7         = gammaDP(1.d0+ALFs(i,k)+bfs)
          GS8         = gammaDP(dms+1.d0+ALFs(i,k))
          tmpdp1      = gammaDP(1.d0+ALFs(i,k)+dms)
          tmpdp2      = tmpdp1/GS2(i,k)*cms
          iLAMs(i,k)  = max(iLAMmin2, (dble(DE(i,k)*QNM(i,k)/NNM(i,k))/tmpdp2)**thrd)
          iLAMs2(i,k) = iLAMs(i,k) *iLAMs(i,k)
          iLAMsB0(i,k)= iLAMs(i,k)**(bfs)
          iLAMsB1(i,k)= iLAMs(i,k)**(bfs+1.d0)
          iLAMsB2(i,k)= iLAMs(i,k)**(bfs+2.d0)
          ckQs1       = afs*(gammaDP(1.d0+ALFs(i,k)+dms+bfs)/tmpdp1)
          vs0(i,k)    = dble(gamfact(i,k))*ckQs1*iLAMsB0(i,k)
       else
          ALFs(i,k)   = 0.d0;  iLAMs(i,k)  = 0.d0;  vs0(i,k)    = 0.d0;  iLAMs2(i,k)= 0.d0
          iLAMsB0(i,k)= 0.d0;  iLAMsB1(i,k)= 0.d0;  iLAMsB1(i,k)= 0.d0
       endif

    ! Graupel:
       if (QGM(i,k)>epsQ .and. grpl_ON) then
          if      (scheme==1) then
             ALFg(i,k)= ALFgfix
             tmpdp1   = gammaDP(1.d0+ALFg(i,k))
             tmpdp2   = gammaDP(4.d0+ALFg(i,k))
             NGM(i,k) = (Nogfix*tmpdp1)**(3./(4.+ALFg(i,k)))*(tmpdp1/tmpdp2*DE(i,k)*     &
                         QGM(i,k)/cmg)**((1.+ALFg(i,k))/(4.+ALFg(i,k)))  !i.e. NXM = f(No,QXM)
          else if (scheme==2) then
             ALFg(i,k)= ALFgfix
          else if (scheme==3) then
             Dg       = (dble(DE(i,k)*QGM(i,k)/NGM(i,k))*icmg)**thrd
             ALFg(i,k)= diagAlpha_v33(Dg,4)
          else if (scheme==4) then
             ALFg(i,k)= solveAlpha_v33(QGM(i,k),NGM(i,k),ZGM(i,k),cmgSP,DE(i,k))
          endif
          GG2(i,k)    = gammaDP(1.d0+ALFg(i,k))
          tmpdp1      = gammaDP(1.d0+ALFg(i,k)+dmg)
          tmpdp2      = tmpdp1/GG2(i,k)*cmg
          ckQg1       = afg*(gammaDP(1.d0+ALFg(i,k)+dmg+bfg)/tmpdp1)
          iLAMg(i,k)  = max(iLAMmin1, (dble(DE(i,k)*QGM(i,k)/NGM(i,k))/tmpdp2)**thrd)
          iLAMg2(i,k) = iLAMg(i,k) *iLAMg(i,k)
          iLAMgB0(i,k)= iLAMg(i,k)**(bfg)
          iLAMgB1(i,k)= iLAMg(i,k)**(bfg+1.d0)
          iLAMgB2(i,k)= iLAMg(i,k)**(bfg+2.d0)
          vg0(i,k)    = dble(gamfact(i,k))*ckQg1*iLAMgB0(i,k)
       else
          ALFg(i,k)   = 0.d0;  iLAMg(i,k)  = 0.d0;  vg0(i,k)  = 0.d0
          iLAMg2(i,k) = 0.d0;  iLAMgB0(i,k)= 0.d0;  iLAMgB1(i,k)= 0.d0;  iLAMgB1(i,k)= 0.d0
       endif

    ! Hail:
       if (QHM(i,k)>epsQ .and. hl_ON) then        !BJP OCT 2010 hail switch
          if      (scheme==1) then
             ALFh(i,k)= ALFhfix
             tmpdp1   = gammaDP(1.d0+ALFh(i,k))
             tmpdp2   = gammaDP(4.d0+ALFh(i,k))
             NHM(i,k) = (Nohfix*tmpdp1)**(3./(4.+ALFh(i,k)))*(tmpdp1/tmpdp2*DE(i,k)*     &
                         QHM(i,k)/cmh)**((1.+ALFh(i,k))/(4.+ALFh(i,k)))  !i.e. NXM = f(No,QXM)
          else if (scheme==2) then
             ALFh(i,k)= ALFhfix
          else if (scheme==3) then
             Dh       = (dble(DE(i,k)*QHM(i,k)/NHM(i,k))*icmh)**thrd
             ALFh(i,k)= diagAlpha_v33(Dh,5)
          else if (scheme==4) then
             ALFh(i,k)= solveAlpha_v33(QHM(i,k),NHM(i,k),ZHM(i,k),cmhSP,DE(i,k))
          endif
          GH2(i,k)    = gammaDP(1.d0+ALFh(i,k))
          tmpdp1      = gammaDP(1.d0+ALFh(i,k)+dmh)
          tmpdp2      = tmpdp1/GH2(i,k)*cmh
          ckQh1       = afh*(gammaDP(1.d0+ALFh(i,k)+dmh+bfh)/tmpdp1)
          iLAMh(i,k)  = max(iLAMmin1, (dble(DE(i,k)*QHM(i,k)/NHM(i,k))/tmpdp2)**thrd)
          iLAMhB0(i,k)= iLAMh(i,k)**(bfh)
          iLAMhB1(i,k)= iLAMh(i,k)**(bfh+1.d0)
          iLAMhB2(i,k)= iLAMh(i,k)**(bfh+2.d0)
          vh0(i,k)    = dble(gamfact(i,k))*ckQh1*iLAMhB0(i,k)
       else
          ALFh(i,k)   = 0.d0;  iLAMh(i,k)  = 0.d0;  vh0(i,k)    = 0.d0
          iLAMhB0(i,k)= 0.d0;  iLAMhB1(i,k)= 0.d0;  iLAMhB1(i,k)= 0.d0
       endif
!------

      ENDIF
    ENDDO
  ENDDO


  DO k= 1,nk-1
    DO i= 1,ixe
      IF (activePoint(i,k)) THEN

 !Calculating ice-phase source/sink terms:

 ! Initialize all source terms to zero:
       QNUvi=0.d0;  QVDvi=0.d0;  QVDvs=0.d0;  QVDvg=0.d0;  QVDvh=0.d0;   QCLci=0.d0
       QCLcs=0.d0;  QCLcg=0.d0;  QCLch=0.d0;  QFZci=0.d0;  QCLri=0.d0;   QMLsr=0.d0
       QCLrs=0.d0;  QCLrg=0.d0;  QMLgr=0.d0;  QCLrh=0.d0;  QMLhr=0.d0;   QFZrh=0.d0
       QMLir=0.d0;  QCLci=0.d0;  QCNig=0.d0;  QCLsr=0.d0;  QCLsh=0.d0;   QCLgr=0.d0
       QCNis=0.d0;  QCLir=0.d0;  QCLis=0.d0;  QCLig=0.d0;  QCLih=0.d0;   QCNgh=0.d0
       QIMsi=0.d0;  QIMgi=0.d0;  QCNsg=0.d0;  QHwet=0.d0

       if (scheme>1) then
         NCLci= 0.d0; NCLcs=0.d0;  NCLcg=0.d0;  NCLch=0.d0;  NFZci=0.d0;   NMLhr=0.d0
         NCLri= 0.d0; NCLrs=0.d0;  NCLrg=0.d0;  NCLrh=0.d0;  NMLsr=0.d0;   NMLgr=0.d0
         NMLir= 0.d0; NSHhr=0.d0;  NNUvi=0.d0;  NVDvi=0.d0;  NCNig=0.d0;   NVDvh=0.d0
         NCLir= 0.d0; NCLis=0.d0;  NCLig=0.d0;  NCLih=0.d0;  NIMsi=0.d0;   NIMgi=0.d0
         NiCNis=0.d0; NsCNis=0.d0; NIMii=0.d0;  NVDvs=0.d0;  NCNsg=0.d0;   NCLgr=0.d0
         NCLss= 0.d0; NCLsr=0.d0;  NCLsh=0.d0;  NCLsrs=0.d0; NCLgrg=0.d0;  NCNgh=0.d0
         NVDvg= 0.d0; NCLirg=0.d0; NCLsrg=0.d0; NCLgrh=0.d0; NrFZrh=0.d0;  NhFZrh=0.d0
         NCLirh=0.d0; NCLsrh=0.d0
       endif

       if (scheme==4) then
         ZrCLri=0.d0; ZrCLrs=0.d0; ZrCLrg= 0.d0; ZrCLrh= 0.d0;ZFZrh= 0.d0;  ZgCNsrg=0.d0
         ZCLyi= 0.d0; ZFZci= 0.d0; ZNUvi=  0.d0; ZCNis= 0.d0; ZiCNig=0.d0;  ZiCLir= 0.d0
         ZiCLis=0.d0; ZiCLig=0.d0; ZiCLih= 0.d0; ZMLir= 0.d0; ZiIMsi= 0.d0; ZiIMgi= 0.d0
         ZCLys= 0.d0; ZCNis= 0.d0; ZsCNsg= 0.d0; ZMLsr= 0.d0; ZsCLsr= 0.d0; ZsCLsh= 0.d0
         ZsCLsrs=0.d0;ZCLss= 0.d0; ZhCLirh=0.d0; ZhCLsrh=0.d0;ZhCLgrh=0.d0; ZgCLgrg=0.d0
         ZCLyg= 0.d0; ZMLgr= 0.d0; ZCNgh=  0.d0; ZgCLirg=0.d0;ZgCLsrg=0.d0; ZgCNig =0.d0
         ZCLyh= 0.d0; ZVDvh= 0.d0; ZhMLhr= 0.d0; ZIMsi= 0.d0; ZIMgi= 0.d0;  ZiCNis =0.d0
         ZsCNis=0.d0; ZgCNsg=0.d0; ZrMLhr= 0.d0
       endif

       Dirg=0.d0; Dirh=0.d0; Dsrs= 0.d0; Dsrg= 0.d0; Dsrh= 0.d0; Dgrg=0.d0; Dgrh=0.d0

       TcSP  = TM(i,k)-TRPL
       Tc    = dble(TcSP)
       if (Tc<-120.d0) print*, '***WARNING*** -- In MULTIMOMENT --  Ambient Temp.(C):',Tc
       Cdiff = (2.2157d-5+0.0155d-5*Tc)*1.d5/dble(PM(i,k))
       MUdyn = 1.72d-5*(393.d0/(dble(TM(i,k)+120.)))*dble(TM(i,k)/TRPL)**1.5d0 !RYp.102
       MUkin = MUdyn/dble(DE(i,k))
       ScTHRD= (MUkin/Cdiff)**thrd       ! i.e. Sc^(1/3)
       Ka    = 2.3971d-2 + 0.0078d-2*Tc                                          !therm.cond.(air)
       Kdiff = dble(9.1018e-11*TM(i,k)*TM(i,k)+8.8197e-8*TM(i,k)-(1.0654e-5)) !therm.diff.(air)
       DEdp  = dble(DE(i,k))
       gam   = dble(gamfact(i,k))

        !Collection efficiencies:
       Eis   = min(0.05d0*exp(0.1d0*Tc),1.d0)     !F95 (Table 1)
       Eig   = min(0.01d0*exp(0.1d0*Tc),1.d0)     !dry (Eig=1.0 for wet growth)
       Eii   = 0.1d0*Eis
       Ess   = Eis;   Eih = Eig;   Esh = Eig
       !NOTE:  Eci=Ecs=Ecg=Eri=Ers=Erg=Erh=1. (constant parameters)
       !       Ech is computed in CLch section

       qvs0  = dble(FOQSA(TRPL,PM(i,k)))      !sat.mix.ratio at 0C
       DELqvs= qvs0-dble(QM(i,k))

   !-------------------------------------------------------------------------------------------!

           ! COLLECTION by snow, graupel, hail:
           !  (i.e. wet or dry ice-categories [=> excludes ice crystals])

           ! Collection by SNOW:
       if (QNM(i,k)>epsQ) then
          ! cloud:
          if (QCM(i,k)>epsQ) then
             tmpdp1= gammaDP(1.d0+bfs+ALFs(i,k))
             tmpdp2= gammaDP(2.d0+bfs+ALFs(i,k))
             tmpdp3= gammaDP(3.d0+bfs+ALFs(i,k))

             QCLcs= dt*gam*afs*cmr*Ecs*0.25d0*PI/DEdp*dble(NCM(i,k)*NNM(i,k))/GC5/       &
                    GS2(i,k)*(GC13*tmpdp3*iLAMc3(i,k)*iLAMsB2(i,k)+2.d0*GC14*tmpdp2*     &
                    iLAMc4(i,k)*iLAMsB1(i,k)+GC15*tmpdp1*iLAMc5(i,k)*iLAMsB0(i,k))

             NCLcs= dt*gam*afs*0.25d0*PI*Ecs*dble(NCM(i,k)*NNM(i,k))/GC5/GS2(i,k)*       &
                    (GC5*tmpdp3*iLAMsB2(i,k)+ 2.d0*GC11*tmpdp2*iLAMc(i,k)*iLAMsB1(i,k)+  &
                    GC12*tmpdp1*iLAMc2(i,k)*iLAMsB0(i,k))

             QCLcs= min(QCLcs, dble(QCM(i,k)))
             NCLcs= min(NCLcs, dble(NCM(i,k)))
          else
             QCLcs= 0.;   NCLcs= 0.
          endif

          ! ice:
          if (QIM(i,k)>epsQ) then
             tmpdp1= vs0(i,k)-vi0(i,k)
             tmpdp2= tmpdp1*tmpdp1
             tmpdp3= sqrt(tmpdp2+0.04d0*vs0(i,k)*vi0(i,k))
             tmpdp4= gammaDP(2.d0+ALFs(i,k));  tmpdp5= gammaDP(3.d0+ALFs(i,k))

             QCLis= dt*cmi/DEdp*PI*6.d0*Eis*dble(NYM(i,k)*NNM(i,k))*tmpdp3/GI2(i,k)/     &
                    GS2(i,k)*(0.5d0*iLAMs2(i,k)*iLAMi3(i,k)+2.d0*iLAMs(i,k)*iLAMi4(i,k)+ &
                    5.d0*iLAMi5(i,k))

             NCLis= dt*0.25d0*PI*Eis*dble(NYM(i,k)*NNM(i,k))*GI2(i,k)*GS2(i,k)*tmpdp3*   &
                    (GI32(i,k)*GS2(i,k)*iLAMi2(i,k)+2.d0*GI31(i,k)*tmpdp4*iLAMi(i,k)*    &
                    iLAMs(i,k)+GI2(i,k)*tmpdp5*iLAMs2(i,k))

             QCLis= min(QCLis, dble(QIM(i,k)))
             NCLis= min(QCLis*dble(NYM(i,k)/QIM(i,k)), NCLis)
          else
             QCLis= 0.;   NCLis= 0.
          endif

          ! snow: (i.e. self-collection [aggregation])  2002-04-22
          NCLss= dt*0.91226d0*Ess*(DEdp*dble(QNM(i,k)))**((2.d0+bfs)*thrd)*             &
                 dble(NNM(i,k))**((4.d0-bfs)*thrd)
          NCLss= min(NCLss, 0.5d0*dble(NNM(i,k)))
          !Note: 0.91226 = I(bfs)*afs*PI^((1-bfs)/3)*des^((-2-bfs)/3) where I(bfs=0.41)=1138.

       else
          QCLcs= 0.d0;   NCLcs= 0.d0;   QCLis= 0.d0;   NCLis= 0.d0;  NCLss= 0.d0
       endif

       ! Collection by GRAUPEL:
       if (QGM(i,k)>epsQ .and. grpl_ON) then
          ! cloud:
          if (QCM(i,k)>epsQ) then
             tmpdp1= gammaDP(1.d0+bfg+ALFg(i,k))
             tmpdp2= gammaDP(2.d0+bfg+ALFg(i,k))
             tmpdp3= gammaDP(3.d0+bfg+ALFg(i,k))

             QCLcg= dt*gam*afg*cmr*Ecg*0.25d0*PI/DEdp*dble(NCM(i,k)*NGM(i,k))/GC5/      &
                    GG2(i,k)*(GC13*tmpdp3*iLAMc3(i,k)*iLAMgB2(i,k)+ 2.d0*GC14*tmpdp2*    &
                    iLAMc4(i,k)*iLAMgB1(i,k)+GC15*tmpdp1*iLAMc5(i,k)*iLAMgB0(i,k))

             NCLcg= dt*gam*afg*0.25d0*PI*Ecg*dble(NCM(i,k)*NGM(i,k))/GC5/GG2(i,k)*      &
                    (GC5*tmpdp3*iLAMgB2(i,k)+2.d0*GC11*tmpdp2*iLAMc(i,k)*iLAMgB1(i,k)+   &
                    GC12*tmpdp1*iLAMc2(i,k)*iLAMgB0(i,k))

             QCLcg= min(QCLcg, dble(QCM(i,k)))
             NCLcg= min(NCLcg, dble(NCM(i,k)))
          else
             QCLcg= 0.d0;   NCLcg= 0.d0
          endif

          ! ice:
          if (QIM(i,k)>epsQ) then
             tmpdp1= vg0(i,k)-vi0(i,k)
             tmpdp2= tmpdp1*tmpdp1
             tmpdp3= sqrt(tmpdp2+0.04d0*vg0(i,k)*vi0(i,k))

             QCLig= dt*cmi/DEdp*PI*6.d0*Eig*dble(NYM(i,k)*NGM(i,k))*tmpdp3/GI2(i,k)/     &
                    GG2(i,k)*(0.5d0*iLAMg2(i,k)*iLAMi3(i,k)+2.d0*iLAMg(i,k)*iLAMi4(i,k)+  &
                    5.d0*iLAMi5(i,k))

             NCLig= dt*0.25d0*PI*Eig*dble(NYM(i,k)*NGM(i,k))*GI2(i,k)*GG2(i,k)*tmpdp3*   &
                    (GI32(i,k)*GG2(i,k)*iLAMi2(i,k)+2.d0*GI31(i,k)*gammaDP(2.d0+ALFg(i,k))* &
                    iLAMi(i,k)*iLAMg(i,k)+GI2(i,k)*gammaDP(3.d0+ALFg(i,k))*iLAMg2(i,k))

             QCLig= min(QCLig, dble(QIM(i,k)))
             NCLig= min(QCLig*dble(NYM(i,k)/QIM(i,k)), NCLig)
          else
             QCLig= 0.d0;   NCLig= 0.d0
          endif

       else
          QCLcg= 0.d0;   QCLrg= 0.d0;   QCLig= 0.d0
          NCLcg= 0.d0;   NCLrg= 0.d0;   NCLig= 0.d0
       endif

       ! Collection by HAIL:
       if (QHM(i,k)>epsQ .and. hl_ON) then          !BJP OCT 2010 hail switch

          iLAMh2= iLAMh(i,k)*iLAMh(i,k)
          Noh(i,k)  = dble(NHM(i,k))/gammaDP(1.d0+ALFh(i,k))/iLAMh(i,k)**(1.d0+ALFh(i,k))
          VENTh(i,k)= Avx*gammaDP(2.d0+ALFh(i,k))*iLAMh(i,k)**(2.d0+ALFh(i,k)) + Bvx*ScTHRD*  &
                     sqrt(gam*afh/MUkin)*gammaDP(2.5d0+bfh*0.5d0+ALFh(i,k))*iLAMh(i,k)**      &
                     (2.5d0+0.5d0*bfh+ALFh(i,k))

         ! cloud:
          if (QCM(i,k)>epsQ) then
             Dh   = (dble(DE(i,k)*QHM(i,k)/NHM(i,k))*icmh)**thrd
             Ech  = exp(-8.68d-7*Dc(i,k)**(-1.6d0)*Dh)    !Z85_A24
             tmpdp1= gammaDP(1.d0+bfh+ALFh(i,k))
             tmpdp2= gammaDP(2.d0+bfh+ALFh(i,k))
             tmpdp3= gammaDP(3.d0+bfh+ALFh(i,k))

             QCLch= dt*gam*afh*cmr*Ech*0.25d0*PI/DEdp*dble(NCM(i,k)*NHM(i,k))/GC5/         &
                    GH2(i,k)*(GC13*tmpdp3*iLAMc3(i,k)*iLAMhB2(i,k)+2.d0*GC14*tmpdp2*        &
                    iLAMc4(i,k)*iLAMhB1(i,k)+GC15*tmpdp1*iLAMc5(i,k)*iLAMhB0(i,k))

             NCLch= dt*gam*afh*0.25d0*PI*Ech*dble(NCM(i,k)*NHM(i,k))/GC5/GH2(i,k)*         &
                    (GC5*tmpdp3*iLAMhB2(i,k)+2.d0*GC11*tmpdp2*iLAMc(i,k)*iLAMhB1(i,k)+GC12* &
                    tmpdp1*iLAMc2(i,k)*iLAMhB0(i,k))

             QCLch= min(QCLch, dble(QCM(i,k)))
             NCLch= min(NCLch, dble(NCM(i,k)))
          else
             QCLch= 0.d0;   NCLch= 0.d0
          endif

          ! rain:
          if (QRM(i,k)>epsQ) then
!         if (QRM(i,k)>epsQ .and. TcSP<0.) then  !2006-02-02
             tmpdp1= (vh0(i,k)-vr0(i,k));  tmpdp2= tmpdp1*tmpdp1
             tmpdp3= sqrt(tmpdp2+0.04d0*vh0(i,k)*vr0(i,k))
             tmpdp4= gammaDP(2.d0+ALFh(i,k))
             tmpdp5= gammaDP(3.d0+ALFh(i,k))

             QCLrh= dt*cmr*Erh*0.25d0*PI/DEdp*dble(NHM(i,k)*NRM(i,k))*iGR5(i,k)/           &
                    GH2(i,k)*tmpdp3*(GR35(i,k)*GH2(i,k) *iLAMr5(i,k)+2.d0*GR34(i,k)*tmpdp4* &
                    iLAMr4(i,k)*iLAMh(i,k)+GR33(i,k)*tmpdp5*iLAMr3(i,k)*iLAMh2)

             NCLrh=  dt*0.25d0*PI*Erh*dble(NHM(i,k)*NRM(i,k))*iGR5(i,k)/GH2(i,k)*tmpdp3*   &
                     (GR32(i,k)*GH2(i,k) *iLAMr2(i,k)+2.d0*GR31(i,k)*tmpdp4*iLAMr(i,k)*     &
                     iLAMh(i,k)+GR5(i,k)*tmpdp5*iLAMh2)

             QCLrh= min(QCLrh, dble(QRM(i,k)))
             NCLrh= min(NCLrh, QCLrh*dble(NRM(i,k)/QRM(i,k)))
          else
             QCLrh= 0.d0;   NCLrh= 0.d0
          endif
          ! ice:
          if (QIM(i,k)>epsQ) then
             tmpdp1= (vh0(i,k)-vi0(i,k));  tmpdp2= tmpdp1*tmpdp1
             tmpdp3= sqrt(tmpdp2+0.04d0*vh0(i,k)*vi0(i,k))
             tmpdp4= gammaDP(2.d0+ALFh(i,k))
             tmpdp5= gammaDP(3.d0+ALFh(i,k))

             QCLih= dt*cmi/DEdp*PI*6.d0*Eih*dble(NYM(i,k)*NHM(i,k))*tmpdp3/GI2(i,k)/       &
                    GH2(i,k)*(0.5d0*iLAMh2*iLAMi3(i,k)+2.d0*iLAMh(i,k)*iLAMi4(i,k)+5.d0*    &
                    iLAMi5(i,k))

             NCLih= dt*0.25d0*PI*Eih*dble(NYM(i,k)*NHM(i,k))*GI2(i,k)*GH2(i,k)*tmpdp3*     &
                    (GI32(i,k)*GH2(i,k)*iLAMi2(i,k)+2.d0*GI31(i,k)*tmpdp4*iLAMi(i,k)*       &
                    iLAMh(i,k)+GI2(i,k)*tmpdp5*iLAMh2)

             QCLih= min(QCLih, dble(QIM(i,k)))
             NCLih= min(QCLih*dble(NYM(i,k)/QIM(i,k)), NCLih)
          else
             QCLih= 0.d0;   NCLih= 0.d0
          endif
          ! snow:
          if (QNM(i,k)>epsQ) then
             tmpdp1= (vh0(i,k)-vs0(i,k));  tmpdp2= tmpdp1*tmpdp1
             tmpdp3= sqrt(tmpdp2+0.04d0*vh0(i,k)*vs0(i,k))
             tmpdp4= iLAMs2(i,k)*iLAMs2(i,k)
             tmpdp5= tmpdp4*iLAMs(i,k)
             tmpdp6= gammaDP(2.d0+ALFh(i,k))
             tmpdp7= gammaDP(3.d0+ALFh(i,k))

             QCLsh= dt*cms/DEdp*PI*6.d0*Esh*dble(NNM(i,k)*NHM(i,k))*tmpdp3/                &
                    GS2(i,k)/GH2(i,k)*(0.5d0*iLAMh2*iLAMs2(i,k)*iLAMs(i,k)+2.d0*iLAMh(i,k)* &
                    tmpdp4+5.d0*tmpdp5)

             NCLsh= dt*0.25d0*PI*Esh*dble(NNM(i,k)*NHM(i,k))*GS2(i,k)*GH2(i,k)*tmpdp3*     &
                    (gammaDP(3.d0+ALFs(i,k))*GH2(i,k)*iLAMs2(i,k)+2.d0*gammaDP(2.d0+ALFs(i,k))* &
                    tmpdp6*iLAMs(i,k)*iLAMh(i,k)+GS2(i,k)*tmpdp7*iLAMh2)

             QCLsh= min(QCLsh, dble(QNM(i,k)))
             NCLsh= min(dble(NNM(i,k)/QNM(i,k))*QCLsh, NCLsh, dble(NNM(i,k)))
          else
             QCLsh= 0.d0;   NCLsh= 0.d0
          endif
          ! wet growth:
          QHwet= max(0.d0, dt*PI2*(DEdp*CHLC*Cdiff*DELqvs-Ka*Tc)*Noh(i,k)/DEdp/(CHLF+      &
                 CPW*Tc)*VENTh(i,k)+(QCLih/Eih+QCLsh/Esh)*(1.d0-CPI*Tc/(CHLF+CPW*Tc)) )

       else
          QCLch= 0.d0;   QCLrh= 0.d0;   QCLih= 0.d0;   QCLsh= 0.d0;   QHwet= 0.d0
          NCLch= 0.d0;   NCLrh= 0.d0;   NCLsh= 0.d0;   NCLih= 0.d0

       endif

       IF (TM(i,k)>TRPL .and. warm_ON) THEN
          !**********!
          !  T > To  !
          !**********!

          ! MELTING of frozen particles:
         !  ICE:
          QMLir   = dble(QIM(i,k))
          QIM(i,k)= 0.
          NMLir   = dble(NYM(i,k))

          !  SNOW:
          if (QNM(i,k)>epsQ) then
             Nos  = dble(NNM(i,k))/GS2(i,k)/iLAMs(i,k)**(1.d0+ALFs(i,k))
             VENTs= Avx*gammaDP(2.d0+ALFs(i,k))*iLAMs(i,k)**(2.d0+ALFs(i,k))+Bvx*ScTHRD* &
                    sqrt(gam*afs/MUkin)*gammaDP(2.5d0+bfs*0.5d0+ALFs(i,k))*iLAMs(i,k)**  &
                    (2.5d0+0.5d0*bfs+ALFs(i,k))
             QMLsr= dt*(PI2/DEdp/CHLF*Nos*VENTs*(Ka*Tc-CHLC*Cdiff*DELqvs) + CPW/CHLF*Tc* &
                    idt*dble(QCLcs+QCLrs))
             QMLsr= min(max(QMLsr,0.d0), dble(QNM(i,k)))
             NMLsr= dble(NNM(i,k)/QNM(i,k))*QMLsr
          else
             QMLsr= 0.d0;   NMLsr= 0.d0
          endif

          !  GRAUPEL:
          if (QGM(i,k)>epsQ .and. grpl_ON) then
             Nog  = dble(NGM(i,k))/GG2(i,k)/iLAMg(i,k)**(1.+ALFg(i,k))
             VENTg= Avx*gammaDP(2.d0+ALFg(i,k))*iLAMg(i,k)**(2.d0+ALFg(i,k))+Bvx*ScTHRD* &
                    sqrt(gam*afg/MUkin)*gammaDP(2.5d0+bfg/2.d0+ALFg(i,k))*iLAMg(i,k)**   &
                    (2.5d0+0.5d0*bfg+ALFg(i,k))
             QMLgr= dt*(PI2/DEdp/CHLF*Nog*VENTg*(Ka*Tc-CHLC*Cdiff*DELqvs)+CPW/CHLF*Tc*   &
                    idt*dble(QCLcg+QCLrg))
             QMLgr= min(max(QMLgr,0.d0), dble(QGM(i,k)))
             NMLgr= dble(NGM(i,k)/QGM(i,k))*QMLgr
          else
             QMLgr= 0.d0;   NMLgr= 0.d0
          endif

          !  HAIL:
          if (QHM(i,k)>epsQ.and.Tc>5.d0 .and. hl_ON) then   !BJP OCT 2010 hail switch
             Noh(i,k)  = dble(NHM(i,k))/gammaDP(1.d0+ALFh(i,k))/iLAMh(i,k)**(1.d0+ALFh(i,k))
             VENTh(i,k)= Avx*gammaDP(2.d0+ALFh(i,k))*iLAMh(i,k)**(2.d0+ALFh(i,k))+Bvx*ScTHRD* &
                        sqrt(gam*afh/MUkin)*gammaDP(2.5d0+bfh*0.5d0+ALFh(i,k))*iLAMh(i,k)**   &
                        (2.5d0+0.5d0*bfh+ALFh(i,k))
             QMLhr= dt*(PI2/DEdp/CHLF*Noh(i,k)*VENTh(i,k)*(Ka*Tc-CHLC*Cdiff*DELqvs) +    &
                    CPW/CHLF*Tc*idt*dble(QCLch+QCLrh))
             QMLhr= min(max(QMLhr,0.d0), dble(QHM(i,k)))
             NMLhr= dble(NHM(i,k)/QHM(i,k))*QMLhr
             if(QCLrh>0.) NMLhr= NMLhr*0.1d0 !Prevents problems when hail is ML & CL
          else
             QMLhr= 0.d0;   NMLhr= 0.d0
          endif

         ! Cold (sub-zero) source/sink terms:
          QNUvi= 0.d0;   QFZci= 0.d0;   QVDvi= 0.d0;   QVDvs= 0.d0;   QVDvg= 0.d0;   QVDvh= 0.d0
          QCLci= 0.d0;   QCLis= 0.d0;   QCNig= 0.d0;   QCNis1=0.d0;   QCNis2=0.d0;   QCNsg= 0.d0
          QCNgh= 0.d0;   QIMsi= 0.d0;   QIMgi= 0.d0;   QCLir= 0.d0;   QCLri= 0.d0;   QCLsr= 0.d0
          QCLrs= 0.d0;   QCLgr= 0.d0;   QCLrg= 0.d0;   QCNis= 0.d0

          if (scheme>1) then
            NNUvi= 0.d0;   NFZci= 0.d0;   NCLci= 0.d0;   NCLgr= 0.d0;   NCLrg= 0.d0
            NCLis= 0.d0;   NVDvi= 0.d0;   NVDvs= 0.d0;   NVDvg= 0.d0;   NVDvh= 0.d0
            NCNsg= 0.d0;   NCNgh= 0.d0;   NiCNis=0.d0;   NsCNis=0.d0;   NCLrs= 0.d0
            NIMsi= 0.d0;   NIMgi= 0.d0;   NCLir= 0.d0;   NCLri= 0.d0;   NCLsr= 0.d0
            NCNig= 0.d0;   NIMii= 0.d0
          endif

       ELSE
          !----------!
          !  T < To  !
          !----------!
          Si    = dble(QM(i,k)/QSI(i,k))
          tmpdp1= dble(TM(i,k)*TM(i,k))
          ! DTD: 10/15/2010 bug fix.  Changed CHLF to CHLS in ABi
          iABi  = 1.d0/( CHLS*CHLS/(Ka*RGASV*tmpdp1) + 1.d0/(DEdp*dble(QSI(i,k))*Cdiff) )

          ! Warm-air-only source/sink terms:
          QMLir= 0.d0;   QMLsr= 0.d0;   QMLgr= 0.d0;   QMLhr= 0.d0
          NMLir= 0.d0;   NMLsr= 0.d0;   NMLgr= 0.d0;   NMLhr= 0.d0

          ! Probabilistic freezing (Bigg) of rain:
          if (QRM(i,k)>Qr_FZrh .and. TcSP<Tc_FZrh .and. hl_ON) then   !JAM/BJP (hl_ON switch added)
             NrFZrh= -dt*Bbigg*(exp(Abigg*Tc)-1.d0)*DEdp*dble(QRM(i,k))/dew
           ! The Rz factor serves to conserve reflectivity when a rain distribution
           !  converts to an distribution with a different shape parameter, alpha.
           !  (e.g. when rain freezes to hail)  The factor Rz non-conserves N while
           !  acting to conserve Z for double-moment.  See Ferrier, 1994 App. D)
             Rz= 1.d0  !N and Z (and Q) are conserved for FZrh with triple-moment
! ! !       if  (QHM(i,k)>epsQ .and. QRM(i,k)>epsQ .and. (scheme==2.or.scheme==3) )             &
! ! !           Rz= (gammaDP(7.d0+ALFh(i,k))*GH2(i,k)*GR33(i,k)*GR33(i,k))/(GR36(i,k)*GR5(i,k)*   &
! ! !                gammaDP(4.d0+ALFh(i,k))*gammaDP(4.d0+ALFh(i,k)))
!!! *** GR36 has not been defined  ***
             NhFZrh= Rz*NrFZrh
             QFZrh = NrFZrh*dble(QRM(i,k)/NRM(i,k))
          else
             QFZrh= 0.d0;   NrFZrh= 0.d0;  NhFZrh= 0.d0
          endif

          !--------!
          !  ICE:  !
          !--------!
          ! Homogeneous freezing of cloud to ice:
          if (QCM(i,k)>epsQ) then
             tmp2= TcSP*TcSP; tmp3= tmp2*TcSP; tmp4= tmp2*tmp2
             JJ  = dble(10.**max(-20.,(-606.3952-52.6611*TcSP-1.7439*tmp2-0.0265*tmp3-   &
                   1.536e-4*tmp4)))
             tmpdp1= 1.d6* (DEdp*dble(QCM(i,k)/NCM(i,k))*icmr) !i.e. Dc(i,k)[cm]**3
             FRAC  = 1.d0-exp(-JJ*PIov6*tmpdp1*dt)
!                  Dc(i,k)  = (DEdp*dble(QCM(i,k)/NCM(i,k))*icmr)**thrd
!                  tmpdp1= (100.d0*Dc(i,k))
!                  FRAC= 1.d0-exp(-JJ*PIov6*tmpdp1*tmpdp1*tmpdp1*dt)
             if (TcSP>-30.) FRAC= 0.d0
             if (TcSP<-50.) FRAC= 1.d0
             QFZci=   FRAC*dble(QCM(i,k))
             NFZci=   FRAC*dble(NCM(i,k))
          else
             QFZci= 0.d0;   NFZci= 0.d0
          endif

          ! Primary ice nucleation:
          NuDEPSOR= 0.d0;   NuCONT= 0.d0;   NNUvi= 0.d0;   QNUvi= 0.d0
          Simax   = min(Si, SxFNC_v33(WW(i,k),TcSP,PM(i,k),QSW(i,k),QSI(i,k),airtype,2))
          tmp1    = (T(i,k)-7.66)
          NNUmax  = max(0.d0, DEdp/mio*dble(Q(i,k)-QSS(i,k))/(1.d0+ck6*dble(QSS(i,k)/    &
                    (tmp1*tmp1))))

          ! Deposition/sorption nucleation:
          if (Tc<-5.d0 .and. Si>1.d0) then
            if (scheme==1) then
              NuDEPSOR= 5.*exp(0.304*(TRPL-max(233.,TM(i,k))))                           !Cooper eqn.
            else
              NuDEPSOR= max(0.d0, 1.d3*exp(12.96d0*(Simax-1.d0)-0.639d0)-dble(NYM(i,k))) !Meyers(1992)
            endif
          endif

          ! Contact nucleation:
          if (QCM(i,k)>epsQ .and. TcSP<-2.) then
             GG     =  1.d0/dew/(RGASV*dble(TM(i,k))/dble((QSW(i,k)*PM(i,k))/      &
                       EPS1)/Cdiff+CHLC/Ka/dble(TM(i,k))*(CHLC/RGASV/dble(TM(i,k))-      &
                       1.d0))                                                !CP00a
             Swmax  =  SxFNC_v33(WW(i,k),TcSP,PM(i,k),QSW(i,k),QSI(i,k),airtype,1)
             ssat   =  min(dble(QM(i,k)/QSW(i,k)), Swmax) -1.d0
             Tcc    =  Tc + GG*ssat*CHLC/Kdiff                               !C86_eqn64
             Na     =  exp(4.11d0-0.262d0*Tcc)                               !W95_eqn60/M92_2.6
             Kn     =  LAMa0*dble(TM(i,k))*p0/(T0*dble(PM(i,k))*Ra)    !W95_eqn59
             PSIa   = -kBoltz*Tcc/(6.d0*pi*Ra*MUdyn)*(1.+Kn)                 !W95_eqn58
             ft     =  0.4d0*(1.d0+1.45d0*Kn+0.4d0*Kn*exp(-1.d0/Kn))*(Ka+2.5d0*Kn*KAPa)/ &
                      (1.d0+3.d0*Kn)/(2.d0*Ka+5.d0*KAPa*Kn+KAPa)             !W95_eqn57
             Dc(i,k)     =  (DEdp*dble(QCM(i,k)/NCM(i,k))*icmr)**thrd
             F1     =  PI2*Dc(i,k)*Na*dble(NCM(i,k))                         !W95_eqn55
             F2     =  Ka/dble(PM(i,k))*(Tc-Tcc)                       !W95_eqn56
             NuCONTA= -F1*F2*RGASV*dble(TM(i,k))/CHLC/DEdp                   !Diffusiophoresis
             NuCONTB=  F1*F2*ft/DEdp                                         !Thermeophoresis
             NuCONTC=  F1*PSIa                                               !Brownian diffusion
             NuCONT =  max(0.d0,(NuCONTA+NuCONTB+NuCONTC)*dt)
          endif

          ! Total primary ice nucleation:
          if (icephase_ON) then
             NNUvi= min(NNUmax, NuDEPSOR + NuCONT )
             QNUvi= mio/DEdp*NNUvi
             QNUvi= min(QNUvi,dble(Q(i,k)))
          endif

          IF (QIM(i,k)>epsQ) THEN

             ! Riming (stochastic collection of cloud):
             if (QCM(i,k)>epsQ) then
                tmpdp1= gammaDP(0.d0+bfi+1.d0+ALFi(i,k))
                tmpdp2= gammaDP(1.d0+bfi+1.d0+ALFi(i,k))
                tmpdp3= gammaDP(2.d0+bfi+1.d0+ALFi(i,k))

                QCLci= dt*gam*afi*cmr*Eci*0.25d0*PI/DEdp*dble(NCM(i,k)*NYM(i,k))/GC5/    &
                       GI2(i,k)*(GC13*tmpdp3*iLAMc3(i,k)*iLAMiB2(i,k)+2.d0*GC14*tmpdp2*  &
                       iLAMc4(i,k)*iLAMiB1(i,k)+GC15*tmpdp1*iLAMc5(i,k)*iLAMiB0(i,k))

                NCLci= dt*gam*afi*0.25d0*PI*Eci*dble(NCM(i,k)*NYM(i,k))/GC5/GI2(i,k)*    &
                       (GC5*tmpdp3*iLAMiB2(i,k)+2.d0*GC11*tmpdp2*iLAMc(i,k)*iLAMiB1(i,k)+&
                       GC12*tmpdp1*iLAMc2(i,k)*iLAMiB0(i,k))

                QCLci= min(QCLci, dble(QCM(i,k)))
                NCLci= min(NCLci, dble(NCM(i,k)))
             else
                QCLci= 0.d0;   NCLci= 0.d0
             endif

             ! Deposition/sublimation:
          !***## Combine expressions:
             GI5  = gammaDP(2.d0+ALFi(i,k))
             GI6  = gammaDP(2.5d0+bfi*0.5d0+ALFi(i,k))
             Noi  = dble(NYM(i,k))/GI2(i,k)/iLAMi(i,k)**(1.d0+ALFi(i,k))
             VENTi= Avx*GI5*iLAMi(i,k)**(2.d0+ALFi(i,k)) + Bvx*ScTHRD*sqrt(gam*afi/      &
                    MUkin)*GI6*iLAMi(i,k)**(2.5d0+0.5d0*bfi+ALFi(i,k))
          !***##
             tmpdp1= dble(TM(i,k)*TM(i,k))
             QVDvi= dt*iABi*(1.d0/DEdp)*(PI2*(Si-1.)*Noi*VENTi - idt*QCLci*CHLS*CHLF/(Ka*RGASV*      &
                    TM(i,k)*TM(i,k)))

             ! Prevent overdepletion of vapor:
             tmpdp1= (dble(T(i,k))-7.66d0)
             VDmax = dble(Q(i,k)-QSS(i,k))/(1.d0+ck6*dble(QSS(i,k))/(tmpdp1*tmpdp1))
             if(Si>=1.d0) then
                QVDvi= min(max(QVDvi,0.d0),VDmax)
             else
                if (VDmax<0.d0) QVDvi= max(QVDvi,VDmax)
               !IF prevents subl.(QVDvi<0 at t) changing to dep.(VDmax>0 at t*)  2005-06-28
             endif

             if (.not. iceDep_ON) QVDvi= 0. !test; suppress depositional growth

             NVDvi= min(0.d0, dble(NYM(i,k)/QIM(i,k))*QVDvi) !dNi/dt=0 for deposition

             ! Conversion to graupel:
             QCNig= max(QCLci-max(0.d0,QVDvi), 0.d0)
             if (.not. grpl_ON) QCNig = 0.d0
             NCNig= DE(i,k)/mgo*QCNig
             NCNig= min(NCNig, QCNig*dble(NYM(i,k)/QIM(i,k)))
             QCLci= max(QCLci-QCNig, 0.d0)

             ! Conversion to snow:
             !   +depostion/riming of ice:
             Di= (dble(DE(i,k)*QIM(i,k)/NYM(i,k))*icmi)**thrd
             mi= DEdp*dble(QIM(i,k)/NYM(i,k));   ri= 0.5d0*Di
             if (mi<=0.5d0*mso.and.abs(0.5d0*mso-mi)>1.d-20) then
                QCNis1= (mi/(mso-mi))*(QVDvi+QCLci)
             else
                QCNis1= (QVDvi+QCLci) + (1.d0-0.5d0*mso/mi)*dble(QIM(i,k))/dt
             endif
             QCNis1= max(0.d0, QCNis1)
             !   +aggregation of ice:
             if(ri<0.5d0*rso) then
                Ki    =  PIov6*Di*Di*vi0(i,k)*Eii*Xdisp
                tmpdp1= (log(ri/rso));  tmpdp2= tmpdp1*tmpdp1*tmpdp1
                QCNis2= -dt*0.5d0*dble(QIM(i,k)*NYM(i,k))*Ki/tmpdp2
             else
                Ki= 0.d0;   QCNis2= 0.d0
             endif
             !   +total conversion rate:
             QCNis =  QCNis1 + QCNis2
             NsCNis=  DEdp/mso*QCNis                                !source for snow (Ns)
             tmpdp1= dble(NYM(i,k))
             NiCNis= (DEdp/mso*QCNis1 + 0.5d0*Ki*tmpdp1*tmpdp1) !sink for ice (Ni)
             NiCNis= min(NiCNis, dble(NYM(i,k)*0.1)) !031028 Prevents overdepl. of NY when final QI>0

             if (.not.(snow_ON)) then
                QCNis= 0.d0; NiCNis= 0.d0; NsCNis= 0.d0  !Suppress SNOW initiation
             endif

             ! 3-component freezing (collisions with rain):
             if (QRM(i,k)>epsQ .and. QIM(i,k)>epsQ) then
                tmpdp1= (vr0(i,k)-vi0(i,k))
                tmpdp2= tmpdp1*tmpdp1
                tmpdp3= sqrt(tmpdp2+0.04d0*vr0(i,k)*vi0(i,k))
                tmpdp4= gammaDP(1.d0+ALFi(i,k))
                tmpdp5= gammaDP(4.d0+ALFi(i,k))
                tmpdp6= gammaDP(5.d0+ALFi(i,k))
                tmpdp7= gammaDP(6.d0+ALFi(i,k))
                QCLir= dt*cmi*Eri*0.25d0*PI/DEdp*dble(NRM(i,k)*NYM(i,k))/GI2(i,k)*         &
                       iGR5(i,k)*tmpdp3*(tmpdp7*GR5(i,k)*iLAMi5(i,k)+2.d0*tmpdp6*GR31(i,k)* &
                       iLAMi4(i,k)*iLAMr(i,k)+tmpdp5*GR32(i,k)*iLAMi3(i,k)*iLAMr2(i,k))
                NCLri= dt*0.25d0*PI*Eri*dble(NRM(i,k)*NYM(i,k))/GI2(i,k)*iGR5(i,k)*        &
                       tmpdp3*(GI32(i,k)*GR5(i,k) *iLAMi2(i,k)+2.d0*GI31(i,k)*GR31(i,k)*    &
                       iLAMi(i,k)*iLAMr(i,k)+tmpdp4*GR32(i,k)*iLAMr2(i,k))

                QCLri= dt*cmr*Eri*0.25d0*PI/DEdp*dble(NYM(i,k)*NRM(i,k))*iGR5(i,k)/        &
                       GI2(i,k)*tmpdp3*(GR35(i,k)*GI2(i,k) *iLAMr5(i,k)+2.d0*GR34(i,k)*     &
                       gammaDP(2.d0+ALFi(i,k))*iLAMr4(i,k)*iLAMi(i,k)+GR33(i,k)*gammaDP(3.d0+   &
                       ALFi(i,k))*iLAMr3(i,k)*iLAMi2(i,k))

               !note: For explicit eqns, both NCLri and NCLir are mathematically identical)
                NCLir= min(QCLir*dble(NYM(i,k)/QIM(i,k)), NCLir)

                QCLri= min(QCLri, dble(QRM(i,k)));  QCLir= min(QCLir, dble(QIM(i,k)))
                NCLri= min(NCLri, dble(NRM(i,k)));  NCLir= min(NCLir, dble(NYM(i,k)))

                !Determine destination of 3-comp.freezing:
                tmpdp1= max(Di,Dr(i,k))  !OPTIM NOTE:  Make sure Di gets calculated
                dey= (dei*Di*Di*Di+dew*Dr(i,k)*Dr(i,k)*Dr(i,k))/(tmpdp1*tmpdp1*tmpdp1)

!                 if (dey<0.5d0*(deg+deh)) then
!                    Dirg= 1.d0;  Dirh= 0.d0
!                 else
!                    Dirg= 0.d0;  Dirh= 1.d0
!                 endif
                if (dey>0.5d0*(deg+deh) .and. Dr(i,k)>Dr_3cmpThrs .and. hl_ON) then   !JAM/BJP hail switch
                   Dirg= 0.d0;  Dirh= 1.d0
                else
                  if(grpl_ON) then
                    Dirg= 1.d0
                  else
                    Dirg = 0.d0
                  endif
                   Dirh= 0.d0
                endif

             else
                QCLir= 0.d0;  NCLir= 0.d0;  QCLri= 0.d0;
                NCLri= 0.d0;  Dirh= 0.d0;   Dirg= 0.d0
             endif

             !  Rime-splintering (ice multiplication):
             ff= 0.
             if(Tc>=-8.d0.and.Tc<=-5.d0) ff= 3.5d8*(Tc +8.d0)*thrd
             if(Tc> -5.d0.and.Tc< -3.d0) ff= 3.5d8*(-3.d0-Tc)*0.5d0
             NIMii= DEdp*ff*QCLci
             NIMsi= DEdp*ff*QCLcs
             NIMgi= DEdp*ff*QCLcg
             QIMsi= mio/DEdp*NIMsi
             QIMgi= mio/DEdp*NIMgi

          ELSE

             QCLci= 0.d0;    QVDvi= 0.d0;   QCNig= 0.d0;   QCNis= 0.d0
             QIMsi= 0.d0;    QIMgi= 0.d0;   QCLri= 0.d0;   QCLir= 0.d0
             if (scheme>1) then
               NCLci= 0.d0;  NVDvi= 0.d0;   NCLir= 0.d0;   NIMsi= 0.d0
               NCNig= 0.d0;  NiCNis=0.d0;   NsCNis=0.d0
               NIMgi= 0.d0;  NIMii= 0.d0;   NCLri= 0.d0
             endif

          ENDIF
          !---------!
          !  SNOW:  !
          !---------!
          IF (QNM(i,k)>epsQ) THEN

             ! Deposition/sublimation:
             tmpdp1= dble(TM(i,k))
             if (scheme==1) then
               Nos = Nosfix
             else
               Nos = dble(NNM(i,k))/GS2(i,k)/iLAMs(i,k)**(1.d0+ALFs(i,k))
             endif
             VENTs = Avx*gammaDP(2.d0+ALFs(i,k))*iLAMs(i,k)**(2.d0+ALFs(i,k))+Bvx*       &
                     ScTHRD*sqrt(gam*afs/MUkin)*gammaDP(2.5d0+bfs*0.5d0+ALFs(i,k))*      &
                     iLAMs(i,k)**(2.5d0+0.5d0*bfs+ALFs(i,k))
             QVDvs = dt*capFact*iABi*(1.d0/DEdp)*(PI2*(Si-1.)*Nos*VENTs - CHLS*CHLF/(Ka*RGASV*       &
                     TM(i,k)*TM(i,k))*QCLcs*idt)  ! Bugfix from Ted M.  Added missing factor of 1/DE

             ! Prevent overdepletion of vapor:
             tmpdp1= (dble(T(i,k))-7.66)
             VDmax = dble(Q(i,k)-QSS(i,k))/(1.d0+ck6*dble(QSS(i,k))/(tmpdp1*tmpdp1))  !KY97_A.33

             if(Si>=1.) then
                QVDvs= min(max(QVDvs,0.d0),VDmax)
             else
          !     QVDvs= max(QVDvs,VDmax)
                if (VDmax<0.d0) QVDvs= max(QVDvs,VDmax)
                !IF prevents subl.(QVDvs<0 at t) changing to dep.(VDmax>0 at t*)  2005-06-28
             endif
             NVDvs= -min(0.d0,dble(NNM(i,k)/QNM(i,k))*QVDvs)  !pos. quantity

             ! Conversion to graupel:
             if (QCLcs > CNsgThres*QVDvs .and. grpl_ON) then
                QCNsg= (deg/(deg-des))*QCLcs
             else
                QCNsg= 0.
             endif
             NCNsg= DEdp/mgo*QCNsg
             NCNsg= min(NCNsg, dble(0.5*NNM(i,k)/QNM(i,k))*QCNsg) !Prevents incorrect Ns-depletion

             ! 3-component freezing (collisions with rain):
              if (QRM(i,k)>epsQ .and. QNM(i,k)>epsQ .and. Tc<-5.d0) then
                tmpdp1=(vs0(i,k)-vr0(i,k))
                tmpdp2= sqrt(tmpdp1*tmpdp1+0.04d0*vs0(i,k)*vr0(i,k))
                tmpdp4= gammaDP(2.d0+ALFs(i,k))
                tmpdp5= gammaDP(3.d0+ALFs(i,k))
                tmpdp6= iLAMs2(i,k)*iLAMs2(i,k)*iLAMs(i,k)

                QCLrs= dt*cmr*Ers*0.25d0*PI/DEdp*dble(NNM(i,k)*NRM(i,k))*iGR5(i,k)/       &
                       GS2(i,k)*tmpdp2*(GR35(i,k)*GS2(i,k) *iLAMr5(i,k)+2.d0*GR34(i,k)*    &
                       tmpdp4*iLAMr4(i,k)*iLAMs(i,k)+GR33(i,k)*tmpdp5*iLAMr3(i,k)*         &
                       iLAMs2(i,k))

                NCLrs=  dt*0.25d0*PI*Ers*dble(NNM(i,k)*NRM(i,k))*iGR5(i,k)/GS2(i,k)*      &
                        tmpdp2*(GR32(i,k)*GS2(i,k) *iLAMr2(i,k)+2.d0*GR31(i,k)*tmpdp4*     &
                        iLAMr(i,k)*iLAMs(i,k)+GR5(i,k)*tmpdp5*iLAMs2(i,k))

                QCLsr= dt*cms*Ers*0.25d0*PI/DEdp*dble(NRM(i,k)*NNM(i,k))/GS2(i,k)*        &
                       iGR5(i,k)*tmpdp2*(gammaDP(6.d0+ALFs(i,k))*GR5(i,k)*tmpdp6+2.d0*       &
                       gammaDP(5.d0+ALFs(i,k))*GR31(i,k)*iLAMs2(i,k)*iLAMs2(i,k)*iLAMr(i,k)+ &
                       gammaDP(4.d0+ALFs(i,k))*GR32(i,k)*iLAMs2(i,k)*iLAMs(i,k)*iLAMr2(i,k))

               !note: For explicit eqns, NCLsr = NCLrs
                NCLsr= min(QCLsr*dble(NNM(i,k)/QNM(i,k)), NCLrs) !031028

                QCLrs= min(QCLrs, dble(QRM(i,k)));  QCLsr= min(QCLsr, dble(QIM(i,k)))
                NCLrs= min(NCLrs, dble(NRM(i,k)));  NCLsr= min(NCLsr, dble(NYM(i,k)))

                ! Determine destination of 3-comp.freezing:
                Dsrs= 0.d0;   Dsrg= 0.d0;    Dsrh= 0.d0
                Ds  = (dble(DE(i,k)*QNM(i,k)/NNM(i,k))*icms)**thrd
                tmpdp1= max(Ds,Dr(i,k)); tmpdp2= tmpdp1*tmpdp1*tmpdp1
                dey= (des*Ds*Ds*Ds + dew*Dr(i,k)*Dr(i,k)*Dr(i,k))/tmpdp2
                if (dey<=0.5d0*(des+deg)) Dsrs= 1.d0                          !snow
                if (dey>0.5d0*(des+deg) .and. dey<0.5d0*(deg+deh)) then
                  if(grpl_ON) then
                    Dsrg= 1.d0 !graupel
                  else
                    if(hl_ON) Dsrh = 1.d0
                  endif
                endif
                if (dey>=0.5*(deg+deh)) then
                   Dsrh= 1.d0                                                 !hail
                   if (.not.hl_ON .or. Dr(i,k)<Dr_3cmpThrs) then
                      if(grpl_ON) Dsrg = 1.d0
                      Dsrh= 0.d0                             
                   endif
                endif


             else
                QCLrs= 0.d0;   QCLsr= 0.d0;   NCLrs= 0.d0;   NCLsr= 0.d0
             endif

          ELSE

             QVDvs= 0.d0;  QCLcs= 0.d0;  QCNsg= 0.d0;  QCLsr= 0.d0;  QCLrs= 0.d0
             NVDvs= 0.d0;  NCLcs= 0.d0;  NCLsr= 0.d0;  NCLrs= 0.d0;  NCNsg= 0.d0

          ENDIF
          !------------!
          !  GRAUPEL:  !
          !------------!
          IF (QGM(i,k)>epsQ .and. grpl_ON) THEN
          
             ! Sublimation:
               ! ** Potential problems result e.g. very large temp drops due to
               !    excessively large sublimation rates.  Needs to be corrected.
               ! Nog  = dble(NGM(i,k))/GG2(i,k)/iLAMg(i,k)**(1.+ALFg(i,k))
               ! VENTg= Avx*gammaDP(2.d0+ALFs(i,k))*iLAMg(i,k)**(2.d0+ALFg(i,k)) + Bvx*ScTHRD*sqrt(gam*afg/MUkin)* &
      !SHOULD BE:
      ! VENTg= Avx*gammaDP(2.d0+ALFg(i,k))*iLAMg(i,k)**(2.d0+ALFg(i,k)) + Bvx*ScTHRD*sqrt(gam*afg/MUkin)* &
               !        gammaDP(2.5d0+bfg/2.d0+ALFg(i,k))*iLAMg(i,k)**(2.5d0+0.5d0*bfg+ALFg(i,k))
               ! QVDvg= min(0.d0, dt*PI2*(Si-1.)*Nog*VENTg*iABi)
               ! NVDvg= -dble(NGM(i,k)/QGM(i,k))*QVDvg
                  
             ! DTD: 06/01/07.  Corrected bug in sublimation rate for graupel.  Previously ALFs was 
             ! used in the VENTg equation instead of ALFg.  Could this have been causing the 
             ! excessively large erroneous sublimation rates?
             ! Further comment: 07/18/2011: The above bug in ABi may have also contributed to the 
             ! excessive sublimation rates, if not being the main cause.  Therefore, the sublimation
             ! rates have been tentatively reinstated.

             Nog  = dble(NGM(i,k))/GG2(i,k)/iLAMg(i,k)**(1.+ALFg(i,k))
             VENTg= Avx*gammaDP(2.d0+ALFg(i,k))*iLAMg(i,k)**(2.d0+ALFg(i,k)) + Bvx*ScTHRD*sqrt(gam*afg/MUkin)* &
                     gammaDP(2.5d0+bfg/2.d0+ALFg(i,k))*iLAMg(i,k)**(2.5d0+0.5d0*bfg+ALFg(i,k))
             QVDvg= min(0.d0, dt*PI2*(1.d0/DEdp)*(Si-1.)*Nog*VENTg*iABi) ! Bugfix from Ted M.  Added missing factor of 1/DE
             if(.not. grpl_ON) QVDvg= 0.d0
             NVDvg= -dble(NGM(i,k)/QGM(i,k))*QVDvg
            

             Dg   = (dble(DE(i,k)*QGM(i,k)/NGM(i,k))*icmg)**thrd

           !Conversion to hail:    (Dho given by S-L limit)
             if (Dg>Dg_CNgh .and. WW(i,k)>w_CNgh .and. hl_ON) then
          !  note:  The Dg threshold for CNgh to occur is a surrogate for the
          !         presense of large graupel (that should convert to hail under riming).
          !         However, while this is meaningful for single-moment, it is not appropriate
          !         for double-moment, since a given Dg may have a large or small quantity
          !         of "large" graupel.  For d-m, need to change condition based on estimate
          !         of the amount of large graupel (incomplete gamma distribution).
                Dho  = 0.01d0*(exp(min(20.d0,-Tc/(1.1d4*DEdp*dble(QCM(i,k)+QRM(i,k))-       &
                       1.3d3*DEdp*dble(QIM(i,k))+1.d0)))-1.d0)
                Dho = min(1.d0, max(0.0001d0,Dho))    !smallest Dho=0.1mm; largest=1m
                ratio= Dg/Dho
                if (ratio>r_CNgh) then
                   QCNgh= (0.5d0*ratio)*(QCLcg+QCLrg+QCLig)
                   QCNgh= min(QCNgh,(QGM(i,k))+QCLcg+QCLrg+QCLig)
                   NCNgh= DEdp*QCNgh*icmh/(Dho*Dho*Dho)
                else
                   QCNgh= 0.
                   NCNgh= 0.
                endif
             endif


             ! 3-component freezing (collisions with rain)
             if (QRM(i,k)>epsQ) then
                tmpdp1= vg0(i,k)-vr0(i,k)
                tmpdp2= sqrt(tmpdp1*tmpdp1 + 0.04d0*vg0(i,k)*vr0(i,k))
                tmpdp3= gammaDP(2.d0+ALFg(i,k))
                tmpdp4= gammaDP(3.d0+ALFg(i,k))
                tmpdp5= gammaDP(4.d0+ALFg(i,k))
                tmpdp6= gammaDP(5.d0+ALFg(i,k))
                tmpdp7= gammaDP(6.d0+ALFg(i,k))
                tmpdp8= iLAMg2(i,k)*iLAMg(i,k)   ! iLAMg^3
                tmpdp9= tmpdp8*iLAMg(i,k)        ! iLAMg^4
                tmpdp10=tmpdp9*iLAMg(i,k)        ! iLAMg^5

                QCLrg= dt*cmr*Erg*0.25d0*PI/DEdp*dble(NGM(i,k)*NRM(i,k))*iGR5(i,k)/     &
                       GG2(i,k)*tmpdp2*(GR35(i,k)*GG2(i,k)*iLAMr5(i,k)+2.d0*GR34(i,k)*   &
                       tmpdp3*iLAMr4(i,k)*iLAMg(i,k)+GR33(i,k)*tmpdp4*iLAMr3(i,k)*       &
                       iLAMg2(i,k))

                NCLrg= dt*0.25d0*PI*Erg*dble(NGM(i,k)*NRM(i,k))*iGR5(i,k)/GG2(i,k)*     &
                       tmpdp2*(GR32(i,k)*GG2(i,k) *iLAMr2(i,k)+2.d0*GR31(i,k)*tmpdp3*    &
                       iLAMr(i,k)*iLAMg(i,k)+GR5(i,k) *tmpdp4*iLAMg2(i,k))

                QCLgr= dt*cms*Erg*0.25*PI/DEdp*dble(NRM(i,k)*NGM(i,k))/GG2(i,k)*        &
                       iGR5(i,k)*tmpdp2*(tmpdp7*GR5(i,k) *tmpdp10+2.d0*tmpdp6*GR31(i,k)* &
                       tmpdp9*iLAMr(i,k)+tmpdp5*GR32(i,k)*tmpdp8*iLAMr2(i,k))

               !(note: For explicit eqns, NCLgr= NCLrg)
                NCLgr= min(QCLgr*dble(NGM(i,k)/QGM(i,k)), NCLrg)

                QCLrg= min(QCLrg, dble(QRM(i,k)));  QCLgr= min(QCLgr, dble(QGM(i,k)))
                NCLrg= min(NCLrg, dble(NRM(i,k)));  NCLgr= min(NCLgr, dble(NGM(i,k)))

               ! Determine destination of 3-comp.freezing:
                tmpdp1= max(Dg,Dr(i,k));  tmpdp2= tmpdp1*tmpdp1*tmpdp1
                dey= (deg*Dg*Dg*Dg + dew*Dr(i,k)*Dr(i,k)*Dr(i,k))/tmpdp2
                if (dey<0.5d0*(deg+deh)) then
                   Dgrg= 1.d0;  Dgrh= 0.d0
                else
                   Dgrg= 0.d0
                   if(hl_ON) Dgrh= 1.d0   !BJP OCT 2010 hail switch
                endif
                if (dey>0.5d0*(deg+deh) .and. Dr(i,k)>Dr_3cmpThrs .and. hl_ON) then
                   Dgrg= 0.d0;  Dgrh= 1.d0
                else
                   Dgrg= 1.d0;  Dgrh= 0.d0
                endif

             else
                QCLgr= 0.d0;  QCLrg= 0.d0;  NCLgr= 0.d0;  NCLrg= 0.d0
             endif

          ELSE

             QVDvg= 0.d0;  QCNgh= 0.d0;  QCLgr= 0.d0;  QCLrg= 0.d0
             NVDvg= 0.d0;  NCNgh= 0.d0;  NCLgr= 0.d0;  NCLrg= 0.d0


          ENDIF
          !---------!
          !  HAIL:  !
          !---------!
          IF (QHM(i,k)>epsQ .and. hl_ON) THEN     !BJP OCT 2010 Hail switch

             ! Sublimation:
               !Potential prolems.
               ! Noh(i,k)  = dble(NHM(i,k))/gammaDP(1.d0+ALFh(i,k))/iLAMh(i,k)**(1.d0+ALFh(i,k))
               ! VENTh(i,k)= Avx*gammaDP(2.d0+ALFs(i,k))*iLAMh(i,k)**(2.d0+ALFh(i,k)) + Bvx*ScTHRD*sqrt(gam*afh/MUkin)* &
     !SHOULD BE:
     ! VENTh(i,k)= Avx*gammaDP(2.d0+ALFh(i,k))*iLAMh(i,k)**(2.d0+ALFh(i,k)) + Bvx*ScTHRD*sqrt(gam*afh/MUkin)* &
               !            gammaDP(2.5d0+bfh*0.5d0+ALFh(i,k))*iLAMh(i,k)**(2.5d0+0.5d0*bfh+ALFh(i,k))
               ! QVDvh= min(0.d0, dt*PI2*(si-1.)*Noh(i,k)*VENTh(i,k)*iABi)
               ! NVDvh= min(dble(NYM(i,k)),-dble(NHM(i,k)/QHM(i,k))*QVDvh)  !(positive) F94_B.56
               
             ! DTD: As for graupel, tentatively reinstated sublimation rates for hail due to aforementioned
             ! bugfix to ABi
               
             Noh(i,k)  = dble(NHM(i,k))/gammaDP(1.d0+ALFh(i,k))/iLAMh(i,k)**(1.d0+ALFh(i,k))
             VENTh(i,k)= Avx*gammaDP(2.d0+ALFh(i,k))*iLAMh(i,k)**(2.d0+ALFh(i,k)) + Bvx*ScTHRD*sqrt(gam*afh/MUkin)* &
                         gammaDP(2.5d0+bfh*0.5d0+ALFh(i,k))*iLAMh(i,k)**(2.5d0+0.5d0*bfh+ALFh(i,k))
             
             QVDvh= min(0.d0, dt*PI2*(1.d0/DEdp)*(si-1.)*Noh(i,k)*VENTh(i,k)*iABi)  ! Bugfix from Ted M.  Added missing factor of 1/DE
             NVDvh= min(dble(NYM(i,k)),-dble(NHM(i,k)/QHM(i,k))*QVDvh)  !(positive) F94_B.56

             ! Wet growth:
             if (QHwet<(QCLch+QCLrh+QCLih+QCLsh) .and. Tc>-30.d0) then ! ERM changed -40 to -30
                QCLih= min(QCLih/Eih, dble(QIM(i,k)))  !change Eih to 1. in CLih
                NCLih= min(NCLih/Eih, dble(NYM(i,k)))  !  "    "
                QCLsh= min(QCLsh/Esh, dble(QNM(i,k)))  !change Esh to 1. in CLsh
                NCLsh= min(NCLsh/Esh, dble(NNM(i,k)))  !  "    "
                tmp1 = QCLrh
                QCLrh= QHwet-(QCLch+QCLih+QCLsh)       !actual QCLrh minus QSHhr
                QSHhr= tmp1-QCLrh                      !QSHhr used here only
                NSHhr= DEdp*QSHhr/(cmr*Drshed*Drshed*Drshed)
             else
                NSHhr= 0.d0
             endif

          ELSE
             QVDvh= 0.d0;   NVDvh= 0.d0;   NSHhr= 0.d0
          ENDIF


       ENDIF  ! ( if Tc<0C Block )

           !------------  End of source/sink term calculation  -------------!

!+++ SUPPRESS PROCESSES FOR TESTING:   +++++++!

! QVDvs= 0.d0; NVDvs= 0.d0
! NCLss= 0.d0
! QCNis= 0.d0; NiCNis= 0.d0; NsCNis= 0.d0  !Suppress SNOW
! QCNig= 0.;  NCNig= 0.  !Suppress GRAUPEL (CNig)

! Suppress melting:
! !   QMLir= 0.;  NMLir= 0.
! !   QMLsr= 0.;  NMLsr= 0.
! !   QMLgr= 0.;  NMLgr= 0.
! !   QMLhr= 0.;  NMLhr= 0.



!+++++++++++++++++++++++++++++++++++++++++++++!

           ! Iterative adjustment of sink (and source) terms to prevent overdepletion:

       do niter= 1,2

          ! (1) Vapor:
          source= dble(Q(i,k)) +ddim(-QVDvi,0.0d0)+dim(-QVDvs,0.0d0)+ddim(-QVDvg,0.0d0)+ &
                  ddim(-QVDvh,0.0d0)
          sink  = QNUvi+ddim(QVDvi,0.0d0)+ddim(QVDvs,0.0d0)
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QNUvi= ratio*QNUvi;   NNUvi= ratio*NNUvi
             if(QVDvi>0.) then
               QVDvi= ratio*QVDvi; NVDvi= ratio*NVDvi
             endif
             if(QVDvs>0.d0) then
               QVDvs=ratio*QVDvs;  NVDvs=ratio*NVDvs
             endif
             QVDvg= ratio*QVDvg;   NVDvg= ratio*NVDvg;   QVDvh= ratio*QVDvh
             NVDvh= ratio*NVDvh;   NVDvh= ratio*NVDvh
          endif

          ! (2) Cloud:
          source= dble(QC(i,k))
          sink  = QCLci+QCLcs+QCLcg+QCLch+QFZci
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QFZci= ratio*QFZci;   NFZci= ratio*NFZci;   QCLci= ratio*QCLci
             NCLci= ratio*NCLci;   QCLcs= ratio*QCLcs;   NCLcs= ratio*NCLcs
             QCLcg= ratio*QCLcg;   NCLcg= ratio*NCLcg;   QCLch= ratio*QCLch
             NCLch= ratio*NCLch
          endif

          ! (3) Rain:
          source= dble(QR(i,k))+QMLsr+QMLgr+QMLhr+QMLir
          sink  = QCLri+QCLrs+QCLrg+QCLrh+QFZrh
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QCLrg= ratio*QCLrg;   QCLri= ratio*QCLri;   NCLri= ratio*NCLri
             QCLrs= ratio*QCLrs;   NCLrs= ratio*NCLrs;   QCLrg= ratio*QCLrg
             NCLrg= ratio*NCLrg;   QCLrh= ratio*QCLrh;   NCLrh= ratio*NCLrh
             QFZrh= ratio*QFZrh;   NrFZrh=ratio*NrFZrh;  NhFZrh=ratio*NhFZrh
             if (ratio==0.d0) then
                Dirg= 0.d0; Dirh= 0.d0; Dgrg= 0.0; Dgrh= 0.d0
                Dsrs= 0.d0; Dsrg= 0.d0; Dsrh= 0.
              endif
          endif

          ! (4) Ice:
          source= dble(QI(i,k))+QNUvi+ddim(QVDvi,0.0d0)+QCLci+QFZci
          sink  = QCNig+QCNis+QCLir+ddim(-QVDvi,0.0d0)+QCLis+QCLig+QCLih+QMLir
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QMLir= ratio*QMLir;    NMLir= ratio*NMLir
             if (QVDvi<0.0d0) then
                QVDvi= ratio*QVDvi; NVDvi= ratio*NVDvi
             endif
             QCNig=  ratio*QCNig;   NCNig=  ratio*NCNig
             QCNis=  ratio*QCNis;   NiCNis= ratio*NiCNis;   NsCNis= ratio*NsCNis
             QCLir=  ratio*QCLir;   NCLir=  ratio*NCLir;    QCLig=  ratio*QCLig
             QCLis=  ratio*QCLis;   NCLis=  ratio*NCLis;    NCLig=  ratio*NCLig ! DTD: added 10/22/09 NCLig adjustment
             QCLih=  ratio*QCLih;   NCLih=  ratio*NCLih
             if (ratio==0.d0) then
                Dirg= 0.d0; Dirh= 0.d0
             endif
          endif

          ! (5) Snow:
          source= dble(QN(i,k))+QCNis+ddim(QVDvs,0.0d0)+QCLis+Dsrs*(QCLrs+QCLsr)+QCLcs
          sink  = ddim(-QVDvs,0.0d0)+QCNsg+QMLsr+QCLsr+QCLsh
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             if(QVDvs<=0.0d0) then
                QVDvs= ratio*QVDvs;   NVDvs= ratio*NVDvs
             endif
             QCNsg= ratio*QCNsg;   NCNsg= ratio*NCNsg;   QMLsr= ratio*QMLsr
             NMLsr= ratio*NMLsr;   QCLsr= ratio*QCLsr;   NCLsr= ratio*NCLsr
             QCLsh= ratio*QCLsh;   NCLsh= ratio*NCLsh
             if (ratio==0.d0) then
                Dsrs= 0.d0; Dsrg= 0.d0; Dsrh= 0.d0
             endif
          endif

          !  (6) Graupel:
          source= dble(QG(i,k))+QCNig+QCNsg+ddim(QVDvg,0.0d0)+Dirg*(QCLri+QCLir)+        &
                  Dgrg*(QCLrg+QCLgr)+QCLcg+Dsrg*(QCLrs+QCLsr)+QCLig
          sink  = ddim(-QVDvg,0.0d0)+QMLgr+QCNgh+QCLgr
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QVDvg= ratio*QVDvg;   NVDvg= ratio*NVDvg;   QMLgr= ratio*QMLgr
             NMLgr= ratio*NMLgr;   QCNgh= ratio*QCNgh;   NCNgh= ratio*NCNgh
             QCLgr= ratio*QCLgr;   NCLgr= ratio*NCLgr
             if (ratio==0.d0) then
                Dgrg= 0.d0; Dgrh= 0.d0
             endif
          endif

          !  (7) Hail:
          source= dble(QH(i,k))+ddim(QVDvh,0.0d0)+QCLch+QCLrh+Dirh*(QCLri+QCLir)+QCLih+  &
                  QCLsh+Dsrh*(QCLrs+QCLsr)+QCNgh+Dgrh*(QCLrg+QCLgr)+QFZrh
          sink  = ddim(-QVDvh,0.0d0)+QMLhr
          sour  = max(source,0.d0)
          if(sink>sour) then
             ratio= sour/sink
             QVDvh= ratio*QVDvh;   NVDvh= ratio*NVDvh;   QMLhr= ratio*QMLhr
             NMLhr= ratio*NMLhr
          endif

       enddo
       !---------------  End of iterative adjustment section.  ------------------!

       IF (scheme>1) THEN

        !Compute N-tendencies for destination categories of 3-comp.freezing:
         NCLirg= 0.d0;  NCLirh= 0.d0;  NCLsrs= 0.d0;  NCLsrg= 0.d0
         NCLsrh= 0.d0;  NCLgrg= 0.d0;  NCLgrh= 0.d0

         if (QCLir+QCLri>0.d0) then
            Di    = (dble(DE(i,k)*QIM(i,k)/NYM(i,k))*icmi)**thrd
            tmpdp1= max(Dr(i,k),Di);  tmpdp2= tmpdp1*tmpdp1*tmpdp1*PIov6
            NCLirg= Dirg*DEdp*dble(QCLir+QCLri)/(deg*tmpdp2)
            NCLirh= Dirh*DEdp*dble(QCLir+QCLri)/(deh*tmpdp2)
         endif

         if (QCLsr+QCLrs>0.d0) then
            Ds    = (dble(DE(i,k)*QNM(i,k)/NNM(i,k))*icms)**thrd
            tmpdp1= max(Dr(i,k),Ds);  tmpdp2= tmpdp1*tmpdp1*tmpdp1*PIov6
            NCLsrs= Dsrs*DEdp*dble(QCLsr+QCLrs)/(des*tmpdp2)
            NCLsrg= Dsrg*DEdp*dble(QCLsr+QCLrs)/(deg*tmpdp2)
            NCLsrh= Dsrh*DEdp*dble(QCLsr+QCLrs)/(deh*tmpdp2)
         endif

         if (QCLgr+QCLrg>0.d0) then
            Dg    = (dble(DE(i,k)*QGM(i,k)/NGM(i,k))*icmg)**thrd
            tmpdp1= max(Dr(i,k),Dg);  tmpdp2= tmpdp1*tmpdp1*tmpdp1*PIov6
            NCLgrg= Dgrg*DEdp*dble(QCLgr+QCLrg)/(deg*tmpdp2)
            NCLgrh= Dgrh*DEdp*dble(QCLgr+QCLrg)/(deh*tmpdp2)
         endif

       ENDIF !(if scheme>1)

       IF (scheme==4) THEN
       !-------------------------------------------------------------------------!
       !             Compute Z-tendencies for all categories:                    !
       !-------------------------------------------------------------------------!

       ! NOTE:  QCLrh, etc. are actually dt*(dQ/dt).  When substituted into
       !        Z-tendency equations, the dts factor such that ZCLrh, etc.
       !        are equal to dt*(dZ/dt).
       !
       !        Z-tendencies due to CLyx and VDvx are all incorporated into
       !        ZCLyx, since dNx/dt=0 for those processes, allowing dQx/dt to
       !        be factored out.  (also IMgi for Zg and IMsi for Zs)

       !GalphaX should be uninitialized if QRX=0
       tmpdp1= DEdp*DEdp

       if(QRM(i,k)>epsQ) Czr= GalphaR(i,k)*tmpdp1/(cmr*cmr)

       if(QIM(i,k)>epsQ) Czi= ((6.d0+ALFi(i,k))*(5.d0+ALFi(i,k))*(4.d0+ALFi(i,k)))/      &
                               ((3.d0+ALFi(i,k))*(2.d0+ALFi(i,k))*(1.d0+ALFi(i,k))) *    &
                               tmpdp1/(cmi*cmi)

       if(QNM(i,k)>epsQ) Czs= ((6.d0+ALFs(i,k))*(5.d0+ALFs(i,k))*(4.d0+ALFs(i,k)))/      &
                               ((3.d0+ALFs(i,k))*(2.d0+ALFs(i,k))*(1.d0+ALFs(i,k))) *    &
                               tmpdp1/(cms*cms)

       if(QGM(i,k)>epsQ) Czg= ((6.d0+ALFg(i,k))*(5.d0+ALFg(i,k))*(4.d0+ALFg(i,k)))/      &
                               ((3.d0+ALFg(i,k))*(2.d0+ALFg(i,k))*(1.d0+ALFg(i,k))) *    &
                               tmpdp1/(cmg*cmg)

       if(QHM(i,k)>epsQ) Czh= ((6.d0+ALFh(i,k))*(5.d0+ALFh(i,k))*(4.d0+ALFh(i,k)))/      &
                               ((3.d0+ALFh(i,k))*(2.d0+ALFh(i,k))*(1.d0+ALFh(i,k))) *    &
                               tmpdp1/(cmh*cmh)

       !Rain:
       if (NRM(i,k)>epsN) then
          tmpdp1= dble(QRM(i,k)/NRM(i,k))
          if (QCLri>epsDQ) then
             ZrCLri= Czr*tmpdp1*(2.d0*QCLri-tmpdp1*NCLri)
             if (ZrCLri<0.) ZrCLri= Czr*QCLri*QCLri/NCLri
          endif
          if (QCLrs>epsDQ) then
              ZrCLrs= Czr*tmpdp1*(2.d0*QCLrs-tmpdp1*NCLrs)
             if (ZrCLrs<0.) ZrCLrs= Czr*QCLrs*QCLrs/NCLrs
          endif
          if (QCLrg>epsDQ) then
             ZrCLrg= Czr*tmpdp1*(2.d0*QCLrg-tmpdp1*NCLrg)
             if (ZrCLrg<0.) ZrCLrg= Czr*QCLrg*QCLrg/NCLrg
          endif
          if (QCLrh>epsDQ) then
             ZrCLrh= Czr*tmpdp1*(2.d0*QCLrh-tmpdp1*NCLrh)
             if (ZrCLrh<0.) ZrCLrh= Czr*QCLrh*QCLrh/NCLrh
          endif
          if (QFZrh>epsDQ) then
             ZFZrh= Czr*tmpdp1*(2.d0*QFZrh-tmpdp1*NrFZrh)
             if (ZFZrh<0.) ZFZrh= Czr*QFZrh*QFZrh/NrFZrh
          endif
       endif  !(NRM > epsN)

           !Ice:
       tmpdp1= DEdp*icmi;  tmpdp2= tmpdp1*tmpdp1
       if (QNUvi>epsDQ) ZNUvi= Galpha2* tmpdp2*QNUvi*QNUvi /NNUvi
       if (QFZci>epsDQ) ZFZci= Galpha2* tmpdp2*QFZci*QFZci /NFZci
       if (QIMsi>epsDQ) ZIMsi= Galpha2* tmpdp2*QIMsi*QIMsi /NIMsi
       if (QIMgi>epsDQ) ZIMgi= Galpha2* tmpdp2*QIMgi*QIMgi /NIMgi
       if (NYM(i,k)>epsN) then
          tmpdp1= dble(QIM(i,k)/NYM(i,k))
          if (QCLci+QCLri+QVDvi/=0.d0) ZCLyi= Czi*2.d0*tmpdp1*(QCLci+QCLri+QVDvi)
          if (QCNis>epsDQ) then
           ! ZiCNis= Czi*(2.d0*dble(QIM(i,k)/NYM(i,k))*QCNis-dble(QIM(i,k)/NYM(i,k))**2.d0*NiCNis)
           ! if (ZiCNis<0.) ZiCNis= Czi*QCNis**2.d0/NiCNis
             ZiCNis= Czi*QCNis*QCNis/NiCNis !prevents "zeroing" ALFi(i,k)
             tmpdp2= cmi*icms
             ZsCNis= tmpdp2*tmpdp2 * ZiCNis  !Type-3 eqn
           ! ZsCNis= Galpha2*(DEdp*QCNis*icms)**2.d0/NsCNis !Prescribed ALFs(i,k) for CNis
          endif
          if (QCNig>epsDQ) then
             ZiCNig= Czi*tmpdp1*(2.d0*QCNig-tmpdp1*NCNig)
             if (ZiCNig<0.) ZiCNig= Czi*QCNig*QCNig/NCNig
             tmpdp1= cmi*icmg
             ZgCNig= tmpdp1*tmpdp1 * ZiCNig  !Type-3 eqn
           ! ZgCNig= Galpha1*(DEdp*QCNig*icmg)**2.d0/NCNig  !Prescribed ALFg(i,k) for CNig
          endif
          if (QMLir>epsDQ) then
             ZMLir= Czi*tmpdp1*(2.d0*QMLir-tmpdp1*NMLir)
             if (ZMLir<0.) ZMLir= Czi*QMLir*QMLir/NMLir
          endif
          if (QCLir>1.e-10) then
             ZiCLir= Czi*tmpdp1*(2.d0*QCLir-tmpdp1*NCLir)
             if (ZiCLir<=0.) ZiCLir= Czi*QCLir*QCLir/NCLir
          endif
          if (QCLis>1.e-10) then
             ZiCLis= Czi*tmpdp1*(2.d0*QCLis-tmpdp1*NCLis)
             if (ZiCLis<=0.) ZiCLis= Czi*QCLis*QCLis/NCLis
          endif
          if (QCLig>1.e-10) then
             ZiCLig= Czi*tmpdp1*(2.d0*QCLig-tmpdp1*NCLig)
             if (ZiCLig<=0.) ZiCLig= Czi*QCLig*QCLig/NCLig
          endif
          if (QCLih>1.e-10) then
             ZiCLih= Czi*tmpdp1*(2.d0*QCLih-tmpdp1*NCLih)
             if (ZiCLih<=0.) ZiCLih= Czi*QCLih*QCLih/NCLih
          endif
       endif  !(NYM > epsN)

       !Snow:
       if (NNM(i,k)>epsN) then
          tmpdp1= dble(QNM(i,k)/NNM(i,k))
          if (QCLcs+QCLis+QVDvs-QIMsi>epsDQ)                                             &
            ZCLys= Czs*2.d0*tmpdp1*(QCLcs+QCLis+QVDvg-QIMsi)
          if (QMLsr>epsDQ) then
             ZMLsr= Czs*tmpdp1*(2.d0*QMLsr-tmpdp1*NMLsr)
             if (ZMLsr<0.) ZMLsr= Czs*QMLsr*QMLsr/NMLsr
          endif
          if (QCNsg>epsDQ) then
             ZsCNsg= Czs*tmpdp1*(2.d0*QCNsg-tmpdp1*NCNsg)
             if (ZsCNsg<0.d0) ZsCNsg= Czs*QCNsg*QCNsg/NCNsg
             ZgCNsg= (cms*cms)/(cmg*cmg) * ZsCNsg  !gives ALFg(i,k)~8 if ALFs(i,k)~2  ok
           ! ZgCNsg= Galpha2*(DEdp*QCNsg*icmg)**2.d0/NCNsg   !Prescribed ALFg(i,k)
          endif
          if (QCLsr>epsDQ) then
             ZsCLsr= Czs*tmpdp1*(2.d0*QCLsr-tmpdp1*NCLsr)
             if (ZsCLsr<0.)  ZsCLsr= Czs*QCLsr*QCLsr/NCLsr
          endif
          if (QCLsh>epsDQ) then
             ZsCLsh= Czs*tmpdp1*(2.d0*QCLsh-tmpdp1*NCLsh)
             if (ZsCLsh<0.) ZsCLsh= Czs*QCLsh*QCLsh/NCLsh
          endif
          if (NCLss>0.d0) ZCLss= Czs*tmpdp1*tmpdp1*NCLss
          if (Dsrs==1.d0) then
             tmpdp2= (DEdp*icms)*(QCLsr+QCLrs)
             ZsCLsrs= Galpha2*tmpdp2*tmpdp2/NCLsrs
          endif
       endif  !(NNM > epsN)

       !Graupel:
       if (NGM(i,k)>epsN) then
          tmpdp1= dble(QGM(i,k)/NGM(i,k))
          if (QCLcg+QCLig+QVDvg-QIMgi/=0.d0)                                             &
            ZCLyg= Czg*2.d0*tmpdp1*(QCLcg+QCLig+QVDvg-QIMgi)
          if (QMLgr>epsDQ) then
             ZMLgr= Czg*tmpdp1*(2.d0*QMLgr-tmpdp1*NMLgr)
             if (ZMLgr<0.) ZMLgr= Czg*QMLgr*QMLgr/NMLgr
          endif
          if (QCNgh>epsDQ) then
             ZCNgh= Czg*tmpdp1*(2.d0*QCNgh-tmpdp1*NCNgh)
             if (ZCNgh<0.) ZCNgh= Czg*QCNgh*QCNgh/NCNgh
          endif
       endif  !(NGM > epsN)
       tmpdp1= (DEdp*icmg)
       if (Dirg==1.d0) then
          tmpdp2 = tmpdp1*(QCLir+QCLri)
          ZgCLirg= Galpha2*tmpdp2*tmpdp2/NCLirg
       endif
       if (Dsrg==1.d0) then
          tmpdp2 = tmpdp1*(QCLsr+QCLrs)
          ZgCLsrg= Galpha2*tmpdp2*tmpdp2/NCLsrg
       endif
       if (Dgrg==1.d0) then
          tmpdp2 = tmpdp1*(QCLgr+QCLrg)
          ZgCLgrg= Galpha2*tmpdp2*tmpdp2/NCLgrg
       endif

       !Hail:
       if (NHM(i,k)>epsN) then
          tmpdp1= dble(QHM(i,k)/NHM(i,k))
          if (QCLch+QCLrh+QCLih+QCLsh+QVDvh.ne.0.d0)                                     &
            ZCLyh= Czh*2.d0*tmpdp1*(QCLch+QCLrh+QCLih+QCLsh+QVDvh)
          if (QMLhr>epsDQ) then
             ZhMLhr= Czh*tmpdp1*(2.d0*QMLhr-tmpdp1*NMLhr)
             if (ZhMLhr<0.d0) ZhMLhr= Czh*QMLhr*QMLhr/NMLhr
             ZrMLhr= ((cmh*cmh)/(cmr*cmr))*ZhMLhr
          endif
       endif  !(NHM > epsN)
       tmpdp1= DEdp*icmh
       if (Dirh==1.d0) then
          tmpdp2= tmpdp1*(QCLir+QCLri)
          ZhCLirh= Galpha2*tmpdp2*tmpdp2/NCLirh
       endif
       if (Dsrh==1.d0) then
          tmpdp2= tmpdp1*(QCLsr+QCLrs)
          ZhCLsrh= Galpha2*tmpdp2*tmpdp2/NCLsrh
       endif
       if (Dgrh==1.d0) then
          tmpdp2= tmpdp1*(QCLsr+QCLrs)
          ZhCLgrh= Galpha2*tmpdp2*tmpdp2/NCLgrh
       endif

       !-------------------------------------------------------------------------!
       !             End of computation of Z-tendencies terms                    !
       !-------------------------------------------------------------------------!
       ENDIF  !if (scheme==4)

       !========================================================================!
       !           Add all source/sink terms to all predicted moments:          !
       !========================================================================!

       ! Q-Source/Sink Terms:
       Q(i,k) = Q(i,k)  +sngl( -QNUvi -QVDvi -QVDvs -QVDvg -QVDvh )
       QC(i,k)= QC(i,k) +sngl( -QCLci -QCLcs -QCLcg -QCLch -QFZci )
       QR(i,k)= QR(i,k) +sngl( -QCLri +QMLsr -QCLrs -QCLrg +QMLgr -QCLrh +QMLhr -QFZrh   &
                                 +QMLir )
       QI(i,k)= QI(i,k) +sngl( QNUvi +QVDvi +QCLci -QCNig +QFZci -QCNis -QCLir -QCLis    &
                                -QCLig -QMLir -QCLih +QIMsi +QIMgi )
       QG(i,k)= QG(i,k) +sngl( QCNsg +QVDvg +QCLcg -QCLgr-QMLgr -QCNgh -QIMgi +QCNig     &
                     +QCLig +Dirg*(QCLri+QCLir) +Dgrg*(QCLrg+QCLgr) +Dsrg*(QCLrs+QCLsr) )
       QN(i,k)= QN(i,k) +sngl( QCNis +QVDvs +QCLcs -QCNsg -QMLsr -QIMsi -QCLsr +QCLis    &
                                -QCLsh +Dsrs*(QCLrs+QCLsr) )
       QH(i,k)= QH(i,k) +sngl( Dirh*(QCLri+QCLir) -QMLhr +QVDvh +QCLch                   &
            +Dsrh*(QCLrs+QCLsr) +QCLih +QCLsh +QFZrh +QCLrh +QCNgh +Dgrh*(QCLrg+QCLgr) )
       T(i,k)= T(i,k) + LFP*sngl(QCLci+QCLri+QCLcs+QCLrs+QFZci-QMLsr+QCLcg+QCLrg-QMLir   &
            -QMLgr-QMLhr+QCLch+QCLrh+QFZrh) +LSP*sngl(QNUvi+QVDvi+QVDvs+QVDvg+QVDvh)

      IF (scheme>1) THEN

       ! N-Source/Sink Terms:
       NC(i,k)= NC(i,k) +sngl( -NCLci -NCLcs -NCLcg -NCLch -NFZci )
       NR(i,k)= NR(i,k) +sngl( -NCLri -NCLrs -NCLrg -NCLrh +NMLsr +NMLgr +NMLhr -NrFZrh  &
                                 +NMLir +NSHhr )
       NY(i,k)= NY(i,k) +sngl( NNUvi +NVDvi -NCNig +NFZci -NCLir -NCLis -NCLig  -NCLih   &
                                -NMLir +NIMsi +NIMgi -NiCNis +NIMii )
       NN(i,k)= NN(i,k) +sngl( NsCNis -NVDvs -NCNsg -NMLsr -NCLss -NCLsr -NCLsh +NCLsrs )
       NG(i,k)= NG(i,k) +sngl( NCNig +NCNsg -NCLgr -NVDvg -NMLgr +NCLirg +NCLsrg         &
                                +NCLgrg -NCNgh )
       NH(i,k)= NH(i,k) +sngl( NhFZrh +NCNgh -NMLhr -NVDvh +NCLirh +NCLsrh +NCLgrh )

      ENDIF  ! (if scheme>1)

      IF (scheme==4) THEN

       ! Z-Source/Sink Terms:
       delZR= sngl(-ZrCLri-ZrCLrs-ZrCLrg -ZrCLrh -ZFZrh +ZMLir +ZMLsr +ZMLgr +ZrMLhr )
       delZI= sngl( ZNUvi +ZFZci +ZIMsi +ZIMgi +ZCLyi -ZiCNis -ZiCNig -ZMLir -ZiCLir     &
                   -ZiCLis-ZiCLig -ZiCLih )
       delZN= sngl( ZCLys -ZMLsr +ZsCNis -ZsCNsg -ZsCLsr -ZsCLsh +ZCLss +ZsCLsrs )
       delZG= sngl( ZCLyg -ZMLgr -ZCNgh +ZgCNig +ZgCNsg +ZgCLirg +ZgCLsrg +ZgCLgrg )
       delZH= sngl( ZCLyh -ZhMLhr +ZCNgh +ZFZrh +ZhCLgrh +ZhCLirh +ZhCLsrh )

       ZR(i,k)= ZR(i,k) + delZR
       ZI(i,k)= ZI(i,k) + delZI
       ZN(i,k)= ZN(i,k) + delZN
       ZG(i,k)= ZG(i,k) + delZG
       ZH(i,k)= ZH(i,k) + delZH


!!===
! !  *** Under what conditions do the following situations occur?  (Why is this needed?)
! !  ***

!             if (ZI(i,k)<epsZ .and. QI(i,k)>epsQ .and. NY(i,k)>epsN) then
!               print*, 'HH3: ',QI(i,k),NY(i,k),ZI(i,k)
!
!             endif


!         !Protect against Z<0 while Q,N>0 (i.e. erroneous depletion of Z):


!            if (ZI(i,k)<epsZ .and. QI(i,k)>epsQ .and. NY(i,k)>epsN) then
!               tmpdp1= DEdp*icmi
!               if(QIM(i,k)>epsQ) then
!                  GalphaI = ((6.d0+ALFi(i,k))*(5.d0+ALFi(i,k))*(4.d0+ALFi(i,k)))/ &
!                            ((3.d0+ALFi(i,k))*(2.d0+ALFi(i,k))*(1.d0+ALFi(i,k)))
!                  Czi= GalphaI*tmpdp1*tmpdp1
!               else
!                  Czi= Galpha2*tmpdp1*tmpdp1
!               endif
!               ZI(i,k)= sngl(Czi)*QI(i,k)*QI(i,k)/NY(i,k)
!            endif


!       etc. for all categories....   [copy from a version of code pre- 05-07-08]
!!===


! Protect against large relative change in Zx:   (i.e. truncation error)
!  If O(delZx) ~ O(Zxm), then recompute Zx{t+1} = f( ALPHA{t}, Qx{t+1}, Nx{t+1} )
       if (ZRM(i,k)>0. .and. NR(i,k)>epsN .and. delZR/=ZR(i,k)) then
         if (abs(delZR/(ZR(i,k)-delZR)) > rthres)                                        &
          ZR(i,k)= sngl(Czr)*QR(i,k)*QR(i,k)/NR(i,k)
       endif
       if (ZIM(i,k)>0. .and. NY(i,k)>epsN .and. delZI/=ZI(i,k)) then
         if (abs(delZI/(ZI(i,k)-delZI)) > rthres)                                        &
           ZI(i,k)= sngl(Czi)*QI(i,k)*QI(i,k)/NY(i,k)
       endif
       if (ZNM(i,k)>0. .and. NN(i,k)>epsN .and. delZN/=ZN(i,k)) then
         if (abs(delZN/(ZN(i,k)-delZN)) > rthres)                                        &
          ZN(i,k)= sngl(Czs)*QN(i,k)*QN(i,k)/NN(i,k)
       endif
       if (ZGM(i,k)>0. .and. NG(i,k)>epsN .and. delZG/=ZG(i,k)) then
         if (abs(delZG/(ZG(i,k)-delZG)) > rthres)                                        &
          ZG(i,k)= sngl(Czg)*QG(i,k)*QG(i,k)/NG(i,k)
       endif
       if (ZHM(i,k)>0. .and. NH(i,k)>epsN .and. delZH/=ZH(i,k)) then
         if (abs(delZH/(ZH(i,k)-delZH)) > rthres)                                        &
          ZH(i,k)= sngl(Czh)*QH(i,k)*QH(i,k)/NH(i,k)
       endif

      ENDIF  !if (scheme==4)

! Ensure that all moments for each category are positively correlated:

      IF (scheme==1) THEN

       if(QC(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QC(i,k)
          T(i,k) = T(i,k) - LCP*QC(i,k)
        endif
          QC(i,k)= 0.
       endif
       if(QR(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QR(i,k)
          T(i,k) = T(i,k) - LCP*QR(i,k)
        endif
          QR(i,k)= 0.
       endif
       if(QI(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QI(i,k)
          T(i,k) = T(i,k) - LSP*QI(i,k)
        endif
          QI(i,k)= 0.
       endif
       if(QN(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QN(i,k)
          T(i,k) = T(i,k) - LSP*QN(i,k)
        endif
          QN(i,k)= 0.
       endif
       if(QG(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QG(i,k)
          T(i,k) = T(i,k) - LSP*QG(i,k)
        endif
          QG(i,k)= 0.
       endif
       if(QH(i,k)<epsQ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QH(i,k)
          T(i,k) = T(i,k) - LSP*QH(i,k)
        endif
          QH(i,k)= 0.
       endif

      ELSE IF (scheme==2 .or. scheme==3) THEN

       if(QC(i,k)<epsQ .or. NC(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QC(i,k)
          T(i,k) = T(i,k) - LCP*QC(i,k)
        endif
          QC(i,k)= 0.;  NC(i,k)= 0.
       endif
       if(QR(i,k)<epsQ .or. NR(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QR(i,k)
          T(i,k) = T(i,k) - LCP*QR(i,k)
        endif
          QR(i,k)= 0.;  NR(i,k)= 0.
       endif
       if(QI(i,k)<epsQ .or. NY(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QI(i,k)
          T(i,k) = T(i,k) - LSP*QI(i,k)
        endif
          QI(i,k)= 0.;  NY(i,k)= 0.
       endif
       if(QN(i,k)<epsQ .or. NN(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QN(i,k)
          T(i,k) = T(i,k) - LSP*QN(i,k)
        endif
          QN(i,k)= 0.;  NN(i,k)= 0.
       endif
       if(QG(i,k)<epsQ .or. NG(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QG(i,k)
          T(i,k) = T(i,k) - LSP*QG(i,k)
        endif
          QG(i,k)= 0.;  NG(i,k)= 0.
       endif
       if(QH(i,k)<epsQ .or. NH(i,k)<epsN) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QH(i,k)
          T(i,k) = T(i,k) - LSP*QH(i,k)
        endif
          QH(i,k)= 0.;  NH(i,k)= 0.
       endif

      ELSE IF (scheme==4) THEN

       if(QC(i,k)<epsQ .or. NC(i,k)<epsN)                  then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QC(i,k)
          T(i,k) = T(i,k) - LCP*QC(i,k)
        endif
          QC(i,k)= 0.;  NC(i,k)= 0.
       endif
       if(QR(i,k)<epsQ .or. NR(i,k)<epsN .or. ZR(i,k)<epsZ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QR(i,k)
          T(i,k) = T(i,k) - LCP*QR(i,k)
        endif
          QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
       endif
       if(QI(i,k)<epsQ .or. NY(i,k)<epsN .or. ZI(i,k)<epsZ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QI(i,k)
          T(i,k) = T(i,k) - LSP*QI(i,k)
        endif
          QI(i,k)= 0.;  NY(i,k)= 0.;  ZI(i,k)= 0.
       endif
       if(QN(i,k)<epsQ .or. NN(i,k)<epsN .or. ZN(i,k)<epsZ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QN(i,k)
          T(i,k) = T(i,k) - LSP*QN(i,k)
        endif
          QN(i,k)= 0.;  NN(i,k)= 0.;  ZN(i,k)= 0.
       endif
       if(QG(i,k)<epsQ .or. NG(i,k)<epsN .or. ZG(i,k)<epsZ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QG(i,k)
          T(i,k) = T(i,k) - LSP*QG(i,k)
        endif
          QG(i,k)= 0.;  NG(i,k)= 0.;  ZG(i,k)= 0.
       endif
       if(QH(i,k)<epsQ .or. NH(i,k)<epsN .or. ZH(i,k)<epsZ) then
        if(clipCons_ON) then
          Q(i,k) = Q(i,k) + QH(i,k)
          T(i,k) = T(i,k) - LSP*QH(i,k)
        endif
          QH(i,k)= 0.;  NH(i,k)= 0.;  ZH(i,k)= 0.
       endif

      ENDIF
      Q(i,k)= max(Q(i,k),0.)

      ENDIF  !if (activePoint)
    ENDDO
  ENDDO

  ! DTD: Stopped here 07/18/2011

  !----------------------------------------------------------------------------------!
  !                    End of ice phase microphysics (Part 2)                        !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                       PART 3: Warm Microphysics Processes                        !
  !                                                                                  !
  !  Equations for warm-rain coalescence based on Cohard and Pinty 2000a,b (QJRMS)   !
  !  Condensation/evaportaion equation based on Kong and Yau 1997 (Atmos-Ocean)      !
  !  Equations for rain reflectivity (ZR) based on Milbrandt and Yau 2005b (JAS)     !
  !----------------------------------------------------------------------------------!

  ! Warm-rain Coallescence:

 IF (warm_ON) THEN

  DO k= 1,nk-1
     DO i= 1,ixe

        DEdp  = dble(DE(i,k))
        RCAUTR= 0.d0;  CCACCR= 0.d0;  Dc(i,k)= 0.d0;  iLAMc(i,k)= 0.d0;  L  = 0.d0
        RCACCR= 0.d0;  CCSCOC= 0.d0;  Dr(i,k)= 0.d0;  iLAMr(i,k)= 0.d0;  TAU= 0.d0
        CCAUTR= 0.d0;  CRSCOR= 0.d0;  SIGc   = 0.d0;  DrINIT    = 0.d0
        iLAMc3(i,k)= 0.d0;  iLAMc6(i,k)= 0.d0;  iLAMr3(i,k)= 0.d0;  iLAMr6= 0.d0

        if (scheme==1) then
           NCM(i,k) = Ncfix
           ALFr     = ALFrfix
           tmpdp1   = gammaDP(1.d0+ALFr)
           tmpdp2   = gammaDP(4.d0+ALFr)
           NRM(i,k) = (Norfix*tmpdp1)**(3./(4.+ALFr))*(tmpdp1/tmpdp2*DE(i,k)*   &
                      QRM(i,k)/cmr)**((1.+ALFr)/(4.+ALFr))  !i.e. NRM = f(No,QRM)
           rainPresent= (QRM(i,k)>eps)
        else if (scheme==2.or.scheme==3) then
           rainPresent= (QRM(i,k)>epsQ .and. NRM(i,k)>epsN)
        else if (scheme==4) then
           rainPresent= (QRM(i,k)>epsQ .and. NRM(i,k)>epsN .and. ZRM(i,k)>epsZ)
        endif

        if (QCM(i,k)>epsQ .and. NCM(i,k)>epsN) then
           iLAMc(i,k) = ((DEdp*dble(QCM(i,k)/NCM(i,k)))/cexc9)**thrd
           iLAMc3(i,k)= iLAMc(i,k)*iLAMc(i,k)*iLAMc(i,k)
           iLAMc6(i,k)= iLAMc3(i,k)*iLAMc3(i,k)
           Dc(i,k)    = iLAMc(i,k)*(GC2/GC1)**thrd
           SIGc  = iLAMc(i,k)*( GC3/GC1- (GC2/GC1)*(GC2/GC1) )**sixth
           L     = 0.027d0*DEdp*QCM(i,k)*(6.25d18*SIGc*SIGc*SIGc*Dc(i,k)-0.4d0)
           if (SIGc>SIGcTHRS) TAU= 3.7d0/(DEdp*dble(QCM(i,k))*(0.5d6*SIGc-7.5d0))
        endif

        if (rainPresent) then
           Dr(i,k)   = (DEdp*dble(QRM(i,k)/NRM(i,k))*icmr)**thrd
        !Drop-size limiter [prevents initially large drops from melted hail]
           if (Dr(i,k)>3.d-3) then
              tmpdp1= (Dr(i,k)-3.d-3);  tmpdp2= (Dr(i,k)/DrMAX); tmpdp3= tmpdp2*tmpdp2*tmpdp2
              NRM(i,k)= NRM(i,k)*max((1.d0+2.d4*tmpdp1*tmpdp1),tmpdp3)
              tmp1= (DE(i,k)*QRM(i,k)/cmrSP)
             IF (scheme == 4) ZRM(i,k)= Galpha5*tmp1*tmp1/NRM(i,k)
              !Note: ALFr=5 represents a sufficiently narrow size distribution for large Dr
              Dr(i,k)= (dble(tmp1/NRM(i,k)))**thrd
           endif
           if (scheme==1 .or. scheme==2) then
              ALFr= ALFrfix
           else if (scheme==3) then
              ALFr= diagAlpha_v33(Dr(i,k),1)
           else if (scheme==4) then
              ALFr= max(ALFrMIN, solveAlpha_v33(QRM(i,k),NRM(i,k),ZRM(i,k),cmrSP,DE(i,k)) )
           endif
           GR1   = gammaDP(ALFr+1.d0)
           GR2   = gammaDP(ALFr+4.d0)
           GR3   = gammaDP(ALFr+7.d0)
           cexr9 = GR2/GR1*cmr
           iLAMr(i,k) = ((DEdp*dble(QRM(i,k)/NRM(i,k)))/cexr9)**thrd
           iLAMr3(i,k)= iLAMr(i,k)*iLAMr(i,k)*iLAMr(i,k);  iLAMr6= iLAMr3(i,k)*iLAMr3(i,k)
        endif

        !  Autoconversion:
        if (QCM(i,k)>epsQ .and. SIGc>SIGcTHRS .and. autoconv_ON) then
           RCAUTR= min( max(L/TAU,0.d0), dble(QCM(i,k))/dt )
           DrINIT= max(83.d-6, 12.6d-4/(0.5d6*SIGc-3.5d0))  !initiation regime Dr
           DrAUT = max(DrINIT, Dr(i,k))                     !init. or feeding DrAUT
           CCAUTR= RCAUTR*DEdp/(cmr*DrAUT*DrAUT*DrAUT)

           ! ---------------------------------------------------------------------------- !
           ! NOTE: The formulation for CCAUTR here (dNr/dt|initiation) does NOT follow
           !       eqn (18) in CP2000a, but rather it comes from the F90 code provided
           !       by J-P Pinty (subroutine: 'rain_c2r2.f90').
           !       (See notes: 2001-10-17; 2001-10-22)
           !
           !       Similarly, the condition for the activation of accretion and self-
           !       collection depends on whether or not autoconversion is in the feeding
           !       regime (see notes 2002-01-07).  This is apparent in the F90 code, but
           !       NOT in CP2000a.
           ! ---------------------------------------------------------------------------- !

           ! cloud self-collection: (dNc/dt_autoconversion)   {CP eqn(25)}
           tmpdp1= dble(NCM(i,k))
           CCSCOC= min(KK2*tmpdp1*tmpdp1*GC3/GC1*iLAMc6(i,k),dble(NCM(i,k))/dt) !{CP00a eqn(25)}
        endif

        ! Accretion, rain self-collection, and collisional breakup:
        if ((dble(QRM(i,k))>1.2d0*max(L,0.d0)/DEdp.or.Dr(i,k)>max(5.d-6,DrINIT))               &
             .and. rainAccr_ON .and. rainPresent) then

           !  Accretion:                                                      !{CP00a eqn(22)}
          !if (.false.) then  !suppress accretion (only)
           if (QCM(i,k)>epsQ.and.L>0.) then
              if (Dr(i,k).ge.100.d-6) then
                 CCACCR= KK1*dble(NCM(i,k)*NRM(i,k))*(GC2/GC1*iLAMc3(i,k)+GR2/GR1*iLAMr3(i,k))
                 RCACCR= cmr/DEdp*KK1*dble(NCM(i,k)*NRM(i,k))*iLAMc3(i,k)*(GC3/GC1*iLAMc3(i,k)+ &
                         GC2/GC1*GR2/GR1*iLAMr3(i,k))
              else
                 CCACCR= KK2*dble(NCM(i,k)*NRM(i,k))*(GC3/GC1*iLAMc6(i,k)+GR3/GR1*iLAMr6)
!                  RCACCR= cmr/DEdp*KK2*dble(NCM(i,k)*NRM(i,k))*iLAMc3(i,k)*                   &
!                          (GC4/GR1*iLAMc6(i,k)+GC2/GC1*GR3/GR1*iLAMr6)
!++  The following calculation of RCACCR avoids overflow:
                 tmp1   = cmr/DEdp
                 tmp2   = KK2*dble(NCM(i,k)*NRM(i,k))*iLAMc3(i,k)
                 RCACCR = tmp1 * tmp2
                 tmpdp1 = GC4/GR1
                 tmpdp1 = dble(tmpdp1)*iLAMc6(i,k)
                 tmp2   = GC2/GC1
                 tmp2   = tmp2*GR3/GR1
                 tmpdp2 = dble(tmp2)*iLAMr6
                 RCACCR = RCACCR * (tmpdp1 + tmpdp2)
!++
              endif
              CCACCR= min(CCACCR,dble(NC(i,k))/dt)
              RCACCR= min(RCACCR,dble(QC(i,k))/dt)
            endif

           !  Rain self-collection:
           tmpdp1= dble(NRM(i,k)*NRM(i,k))
           if (Dr(i,k).ge.100.d-6) then
              CRSCOR= KK1*tmpdp1*GR2/GR1*iLAMr3(i,k)                       !{CP00a eqn(24)}
           else
              CRSCOR= KK2*tmpdp1*GR3/GR1*iLAMr6                            !{CP00a eqn(25)}
           endif

           !  Raindrop breakup:                                            !{CP00a eqn(26)}
           Ec= 1.
           if (Dr(i,k) >=  600.d-6) Ec= exp(-2.5d3*(Dr(i,k)-6.d-4))
           if (Dr(i,k) >= 2000.d-6) Ec= 0.d0
           CRSCOR= min(Ec*CRSCOR,dble(0.5*NR(i,k))/dt) !0.5 prevents depletion of NR

        endif  !accretion/self-collection/breakup

        ! Prevent overdepletion of cloud:
         source= MAX(dble(QC(i,k)),0D0)        ! Changed by WYH
         sink  = (RCAUTR+RCACCR)*dt
        if (sink>source) then
           ratio = source/sink
           RCAUTR= ratio*RCAUTR
           RCACCR= ratio*RCACCR
           CCACCR= ratio*CCACCR
        endif

        IF (scheme==4) THEN
        !Zr tendencies:
          DZrDt= 0.d0;  DZrDt1= 0.d0;  DZrDt2= 0.d0;  DZrDt3= 0.d0
          GalphaR(i,k)= ((6.d0+ALFr)*(5.d0+ALFr)*(4.d0+ALFr))/                           &
                        ((3.d0+ALFr)*(2.d0+ALFr)*(1.d0+ALFr))
          tmpdp1 = DEdp*icmr
          Czr    = GalphaR(i,k)*tmpdp1*tmpdp1
          if (RCAUTR>epsDQ .and. CCAUTR>0.d0) then
             tmpdp1= DEdp*RCAUTR*icmr
             DZrDt1= GalphaRaut*tmpdp1*tmpdp1/CCAUTR
          endif
          if (NRM(i,k)>epsN) then
             tmpdp2 = dble(QRM(i,k)/NRM(i,k))
          else
             tmpdp2 = 0.
          endif
          if (RCACCR>epsDQ) DZrDt2= Czr* 2.d0*tmpdp2*RCACCR
          if (CRSCOR>0.d0 ) DZrDt3= Czr*tmpdp2*tmpdp2*CRSCOR
          DZrDt  = DZrDt1 + DZrDt2 + DZrDt3
          ZR(i,k)= max(0., ZR(i,k) + sngl(DZrDt)*DT_sp)
        ENDIF  !(if scheme==4)

        ! Apply tendencies:
        QC(i,k)= max(0., QC(i,k)+sngl(-RCAUTR-RCACCR)*DT_sp )
        NC(i,k)= max(0., NC(i,k)+sngl(-CCACCR-CCSCOC)*DT_sp )
        QR(i,k)= max(0., QR(i,k)+sngl( RCAUTR+RCACCR)*DT_sp )
        NR(i,k)= max(0., NR(i,k)+sngl( CCAUTR-CRSCOR)*DT_sp )
        !(Z-tend applied above)

        if (QR(i,k)>epsQ .and. NR(i,k)>epsN .and. scheme>1) then
           if (scheme==4) then
            !Protect against large relative change in Zr:
            ! If O(delZx) ~ O(Zxm), then recompute Zr{t+1} = f( ALPHA{t},Qr{t+1},Nr{t+1} )
             delZR= sngl(DZrDt)*DT_sp
             if (ZRM(i,k)>0. .and. NR(i,k)>epsN .and. delZR/=ZR(i,k)) then
               if (abs(delZR/(ZR(i,k)-delZR)) > rthres)                                   &
                ZR(i,k)= sngl(Czr)*QR(i,k)*QR(i,k)/NR(i,k)
             endif
           endif
           Dr(i,k) = (DEdp*dble(QR(i,k)/NR(i,k))*icmr)**thrd
           if (Dr(i,k)>3.d-3) then
              tmpdp1= (Dr(i,k)-3.d-3);   tmpdp2= tmpdp1*tmpdp1
              tmpdp3= (Dr(i,k)/DrMAX);   tmpdp4= tmpdp3*tmpdp3*tmpdp3
              NR(i,k)= NR(i,k)*sngl(max((1.d0+2.d4*tmpdp2),tmpdp4))
              if (scheme==4) ZR(i,k)= sngl(Czr)*QR(i,k)*QR(i,k)/NR(i,k) !uses previous ALPHA
           endif
           if (scheme==4) then
             !Protect against Z<0 while Q,N>0 (i.e. erroneous depletion of Z):
              if (ZR(i,k)<epsZ .and. QR(i,k)>epsQ .and. NR(i,k)>epsN)                    &
               ZR(i,k)= sngl(Czr)*QR(i,k)*QR(i,k)/NR(i,k)
              if (QR(i,k)<epsQ.or.NR(i,k)<epsN.or.ZR(i,k)<epsZ) then
                if(clipCons_ON) then
                 Q(i,k) = Q(i,k) + QR(i,k)
                 T(i,k) = T(i,k) - LCP*QR(i,k)
                endif
                 QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
              endif
           endif
        else if (scheme>1) then
           QR(i,k)= 0.;   NR(i,k)= 0.
           if (scheme==4) ZR(i,k)= 0.
        endif  !(Qr,Nr>eps ; scheme>1)

     ENDDO
  ENDDO

  ! Condensation/Evaporation:

  DO k=1,nk-1
     DO i=1,ixe

        DEdp    = dble(DE(i,k))
!         DEo     = dble(DE(i,nk))
         DEo     = dble(DE(i,2))
        gam     = sqrt(DEo/DEdp)
        QSS(i,k)= FOQSA(T(i,k), P(i,k))  ! Re-calculate QS with new T (w.r.t. liquid)
       !----
       !The following removes a fraction of the supersaturation water vapor.  The purpose of this
       ! is to reduce the amount of condensed water, ultimately to reduct the total precipitation.
       ! Note -- there is no physical justification for this within the context of the microphysics.
       ! This adjustment is done entirely to reduce the precipitation bias in the GEM-LAM-2.5.
       ! Q(i,k)= Q(i,k) - max(0., qReducFact*(Q(i,k)-QSS(i,k)))
       !----
        ssat    = dble(Q(i,k)/QSS(i,k)-1.)
        Tc      = dble(T(i,k)-TRPL)
        Cdiff   = max(1.62d-5, (2.2157d-5 + 0.0155d-5*Tc)) *1.d5/dble(P(i,k))
        MUdyn   = max(1.51d-5, (1.7153d-5 + 0.0050d-5*Tc))
        MUkin   = MUdyn/DEdp
        Ka      = max(2.07d-2, (2.3971d-2 + 0.0078d-2*Tc))
        ScTHRD  = (MUkin/Cdiff)**thrd ! i.e. Sc^(1/3)

        !Condensation/evaporation:
        ! Capacity of evap/cond in one time step is determined by saturation
        ! adjustment technique [KY97 App.A].  Equation for rain evaporation rate
        ! comes from CP00a.  Explicit condensation rate is not considered
        ! (as it is in Z85), but rather complete removal of supersaturation
        ! is assumed.

        X= dble(Q(i,k)-QSS(i,k))
        if (scheme==1) then
           ALFr   = ALFrfix
           tmpdp1 = gammaDP(1.d0+ALFr)
           tmpdp2 = gammaDP(4.d0+ALFr)
           NR(i,k)= (Norfix*tmpdp1)**(3./(4.+ALFr))*(tmpdp1/tmpdp2*DE(i,k)*   &
                     QR(i,k)/cmr)**((1.+ALFr)/(4.+ALFr))  !i.e. NR = f(No,QR)
           rainPresent= (QR(i,k)>eps)
        else if (scheme==2 .or. scheme==3) then
           rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN)
        else if (scheme==4) then
           rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN .and. ZR(i,k)>epsZ)
        endif

        IF(X>0.d0 .or. QC(i,k)>epsQ .or. rainPresent) THEN
           tmp1 = (T(i,k)-35.86)
           X    = X/dble(1.+ck5*QSS(i,k)/(tmp1*tmp1))
           ES   =  dble(QSW(i,k)*PM(i,k)/62.2)     !change to 0.622 here and below!  ****

           if (X<dble(-QC(i,k))) then

              D= 0.
              if(rainPresent) then

                 if(QM(i,k)<QSW(i,k)) then
                    MUkin= (1.715d-5+5.d-8*Tc)/DEdp
                    ! Rain evap: (F94)

                    Dr(i,k)= (dble(DE(i,k)*QR(i,k)/NR(i,k))*icmr)**thrd
                    if (scheme==1 .or. scheme==2) then
                       ALFr= ALFrfix
                    else if (scheme==3) then
                       ALFr= diagAlpha_v33(Dr(i,k),1)
                    else if (scheme==4) then
                       ALFr= max(ALFrMIN, solveAlpha_v33(QR(i,k),NR(i,k),ZR(i,k),cmrSP,DE(i,k)) )
                    endif
                    GR1  = gammaDP(ALFr+1.d0)
                    GR2  = gammaDP(ALFr+4.d0)
                    GR3  = gammaDP(ALFr+7.d0)
                    GR16 = gammaDP(ALFr+2.d0)
                    GR17 = gammaDP(2.5d0+ALFr+0.5d0*bfr)
                    cexr4= 1.d0+ALFr
                    cexr5= 2.d0+ALFr
                    cexr6= 2.5d0+ALFr+0.5d0*bfr
                    cexr9= GR2/GR1*cmr
                    LAMr = (cexr9*dble(NR(i,k)/QR(i,k))/DEdp)**thrd
                    if (scheme==1) then
                       Nor = Norfix
                    else  !if scheme=2,3,or 4
                       Nor = dble(NR(i,k))*LAMr**(0.5*cexr4)/GR1*LAMr**(0.5*cexr4)
                       !note: above coding prevents overflow
                    endif
                    VENTr= Avx*GR16/LAMr**cexr5 + Bvx*ScTHRD*sqrt(gam*afr/MUkin)*GR17/   &
                           (LAMr+ffr*0.5)**cexr6
                    tmpdp1= dble(T(i,k)*T(i,k))
                    ! Further bug fix DTD: 10/15/2010 -- changed CHLF to CHLC in ABw
                    ABw   = CHLC*CHLC/(Ka*RGASV*tmpdp1)+1.d0/(DEdp*dble(QSS(i,k))*Cdiff)
                    QREVP = -dt*(1.d0/DEdp)*(PI2*ssat*Nor*VENTr/ABw)  ! Bugfix from Ted M.  Added missing factor of 1/DE
                !!  QREVP= 0.d0  !to suppress evaporation of rain
                    if (dble(QR(i,k))>QREVP) then             !Note: QREVP is [(dQ/dt)*dt]
                       DEL= -QREVP
                    else
                       DEL= -dble(QR(i,k))
                    endif
                    D= max(X+dble(QC(i,k)), DEL)
                 endif  !QM< QSM
              endif   !QR<eps & NR<eps
              X= D - dble(QC(i,k))

              IF (scheme==4) THEN
                 DQrDt  = D/dble(dt)
                 if (rainPresent) then
                    DNrDt  = DQrDt*dble(NR(i,k)/QR(i,k))
                    ALFr   = solveAlpha_v33(QR(i,k),NR(i,k),ZR(i,k),cmrSP,DE(i,k))
                    GalphaR(i,k)= ((6.d0+ALFr)*(5.d0+ALFr)*(4.d0+ALFr))/                 &
                                  ((3.d0+ALFr)*(2.d0+ALFr)*(1.d0+ALFr))
                    tmpdp1= DEdp*icmr
                    tmpdp2= dble(QR(i,k)/NR(i,k))
                    DZrDt  = GalphaR(i,k)*tmpdp1*tmpdp1*tmpdp2*(2.d0*DQrDt - tmpdp2*DNrDt)
                    QR(i,k)= max(0., QR(i,k) + sngl(D)       )  !note: D = dt*DQrDt
                    NR(i,k)= max(0., NR(i,k) + DT_sp*sngl(DNrDt))
                    ZR(i,k)= max(0., ZR(i,k) + DT_sp*sngl(DZrDt))
                 else
                    QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
                 endif
              ELSE !(moments)
                 QR(i,k)= QR(i,k) + sngl(D)
                 if (QR(i,k)>0. .and. scheme>1)                                          &
                   NR(i,k)= max(0.,NR(i,k)+sngl(D)*NR(i,k)/QR(i,k)) !(dNr/dt)|evap
                 ! The above expression of (dNr/dt)|evap is from Ferrier, 1994.
                 ! In CP2000a, Nr is not affected by evap. (except if Qr goes to zero).
              ENDIF

              QC(i,k)= 0.;   NC(i,k)= 0.
              T(i,k) = T(i,k) + LCP*sngl(X)
              Q(i,k) = Q(i,k) - sngl(X)

           else  ![if(X >= -QC)]

              ! Nucleation of cloud droplets:
              if (ssat>0.d0 .and. WW(i,k)>0. .and. scheme>1) NC(i,k)=  &
                   max(NC(i,k),NccnFNC_v33(WW(i,k),TM(i,k),PM(i,k),airtype))

              ! All supersaturation is removed (condensed onto cloud field).
              T(i,k) = T(i,k)  + LCP*sngl(X)
              Q(i,k) = Q(i,k)  - sngl(X)
              QC(i,k)= QC(i,k) + sngl(X)
                        
              if (X<0.d0 .and. scheme>1) then
                  if (QC(i,k)>0.) then
                     NC(i,k)= max(0., NC(i,k) + sngl(X)*NC(i,k)/QC(i,k) ) !(dNc/dt)|evap
                  else
                     NC(i,k)= 0.
                  endif
              endif
              if (QC(i,k)>0..and.NC(i,k)==0.) NC(i,k)= 1.e7 !prevents non-zero_Q & zero_N

              ! Homogeneous freezing of cloud to ice:
              !  Note:  This needs to be calculated here, as well as in the ice-phase
              !         section, in case water condenses at a very low temperature.
              if (QC(i,k)>epsQ .and. T(i,k)<243.15 .and. icephase_ON) then

                IF (SCHEME==1) THEN
                  QFZCI  = DBLE(QC(I,K))
                  QI(I,K)= QI(I,K) + SNGL(QFZCI)
                  QC(I,K)= QC(I,K) - SNGL(QFZCI)
                  T(I,K) = T(I,K)  - SNGL(QFZCI)*LFP
                ELSE
                  TcSP= T(i,k)-TRPL
                  tmp2= TcSP*TcSP; tmp3= tmp2*TcSP; tmp4= tmp2*tmp2
                  JJ  = dble(10.**max(-20.,(-606.3952-52.6611*TcSP-1.7439*tmp2-0.0265*   &
                        tmp3-1.536e-4*tmp4)))
                  tmpdp1= 1.d6* (DEdp*dble(QC(i,k)/NC(i,k))*icmr) !i.e. Dc(i,k)[cm]**3
                  FRAC= 1.d0-exp(-JJ*PIov6*tmpdp1*dt)
                  if (TcSP>-30.) FRAC= 0.d0
                  if (TcSP<-50.) FRAC= 1.d0
                  QFZci=   FRAC*dble(QC(i,k))
                  NFZci=   FRAC*dble(NC(i,k))
                  IF (scheme==4 .and. FRAC>0. ) THEN
                     if (QI(i,k)>epsQ) then
                        ALFi(i,k) = solveAlpha_v33(QI(i,k),NY(i,k),ZI(i,k),cmiSP,DE(i,k))
                        GalphaI = ((6.d0+ALFi(i,k))*(5.d0+ALFi(i,k))*(4.d0+ALFi(i,k)))/  &
                                  ((3.d0+ALFi(i,k))*(2.d0+ALFi(i,k))*(1.d0+ALFi(i,k)))
                     else
                        GalphaI= Galpha2
                     endif
                     QI(i,k)= QI(i,k) + sngl(QFZci)
                     NY(i,k)= NY(i,k) + sngl(NFZci)
                     QC(i,k)= max(0., QC(i,k)-sngl(QFZci))
                     NC(i,k)= max(0., NC(i,k)-sngl(NFZci))
                     T(i,k) = T(i,k)  + sngl(QFZci)*LFP
                     if (QI(i,k)>epsQ .and. NY(i,k)>epsN .and. ZI(i,k)>epsZ) then
                        tmp1= DE(i,k)*QI(i,k)/cmiSP
                        ZI(i,k)= sngl(GalphaI)*tmp1*tmp1/NY(i,k) !**
                        !**Note: The above SHOULD be Type 1 dZ/dt eqn (MY2004b) but it was problematic
                     else
                        tmpdp1 =  DEdp*QFZci*icmi
                        ZFZci  = Galpha1*tmpdp1*tmpdp1/NFZci
                        ZI(i,k)= ZI(i,k) + sngl(ZFZci)
                     endif
                  ELSE
                     QI(i,k)= QI(i,k) + sngl(QFZci)
                     NY(i,k)= NY(i,k) + sngl(NFZci)
                     QC(i,k)= QC(i,k) - sngl(QFZci)
                     NC(i,k)= NC(i,k) - sngl(NFZci)
                     T(i,k) = T(i,k)  + sngl(QFZci)*LFP  ! Bug fix: DTD. Changed sign from - to + for QFZci term
                  ENDIF
                ENDIF  !if (scheme=1)

              endif !homogeneous freezing

           endif

        ENDIF

        !  Protect against negative values due to overdepletion:
        if (scheme==1 .and. QR(i,k)<epsQ) then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QR(i,k)
           T(i,k) = T(i,k) - QR(i,k)*LCP
          endif
           QR(i,k)= 0.;
        else if ((scheme==2.or.scheme==3) .and. (QR(i,k)<epsQ.or.NR(i,k)<epsN))  then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QR(i,k)
           T(i,k) = T(i,k) - QR(i,k)*LCP
          endif
           QR(i,k)= 0.;  NR(i,k)= 0.
        else if (scheme==4 .and. (QR(i,k)<epsQ.or.NR(i,k)<epsN.or.ZR(i,k)<epsZ)) then
          if(clipCons_ON) then
           Q(i,k) = Q(i,k) + QR(i,k)
           T(i,k) = T(i,k) - QR(i,k)*LCP
          endif
           QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
        endif

     ENDDO
  ENDDO    !cond/evap [k-loop]

 ENDIF  !if warm_ON

  !----------------------------------------------------------------------------------!
  !                    End of warm-phase microphysics (Part 3)                       !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  !                            PART 4:  Sedimentation                                !
  !----------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------!
  ! Sedimentation is computed using a modified version of the box-Lagrangian         !
  ! scheme (blg4.ftn).  Sedimentation is only computed for columns containing        !
  ! non-zero hydrometeor quantities (at at least one level).                         !
  !----------------------------------------------------------------------------------!

 IF (sedi_ON) THEN

   LR= 0.;  SR= 0.

!--  RAIN sedimentation:  -------------------------!

! The following computes rain sedimentation.  It is slighlty different from the
! sedimentation if i,s,g,h for a few reasons.  First, there is a third fall
! velocity parmaeter, ffr, which makes the calculations for the bulk velocities
! more complicated.  Second, drop break-up is done slightly differently.
! Also, rain is converted to cloud after sedimentation if Dr is small.
! If a different set of fall velocity parameter for rain is adopted, using
! only two parameter, as with i,g,s,h, a common sedimentation subroutine could
! be coded, with a tailscript to accomodate the differences in breakup calculations.

  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns_v33(QR,ni,nk,epsQ,counter,activeColumn,slabHASmass,ktop_sedi,ixe)

   IF (slabHASmass) THEN

     DO nnn= 1,npassr
        !RHOQX= DE*QR
        ! Changed by DTD
        DO k=1,nk-1
          DO i=1,ixe
            RHOQX(i,k)=DE(i,k)*QR(i,k)
          END DO
        END DO
       VVQ= 0.;  VVN= 0.;  VVZ= 0.;  VqMax= 0.; VnMax= 0.; VzMax= 0.
       do a= 1,counter
         i=activeColumn(a)
         do k= 1,nk-1
           if (scheme==1) then
              rainPresent= (QR(i,k)>eps)
              ALFr   = ALFrfix
              tmpdp1 = gammaDP(1.d0+ALFr)
              tmpdp2 = gammaDP(4.d0+ALFr)
              NR(i,k)= (Norfix*tmpdp1)**(3./(4.+ALFr))*(tmpdp1/tmpdp2*DE(i,k)*   &
                        QR(i,k)/cmr)**((1.+ALFr)/(4.+ALFr))  !i.e. NR = f(No,QR)
           else if (scheme==2) then
              rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN)
              ALFr= ALFrfix
           else if (scheme==3) then
              rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN)
              if (rainPresent) ALFr= diagAlpha_v33(Dr(i,k),1)
           else if (scheme==4) then
              rainPresent= (QR(i,k)>epsQ .and. NR(i,k)>epsN .and. ZR(i,k)>epsZ)
              if (rainPresent)    &
                ALFr= max(ALFrMIN, solveAlpha_v33(QR(i,k),NR(i,k),ZR(i,k),cmrSP,DE(i,k)) )
           endif
           if (rainPresent) then
              cexr1 = 1.d0+dmr+ALFr+bfr
              cexr2 = 1.d0+ALFr+dmr
              cexr3 = 1.d0+bfr+ALFr
              cexr5 = 7.d0+bfr+ALFr
              cexr4 = 1.d0+ALFr
              cexr6 = 7.d0+ALFr
              cexr9 = cmr*gammaDP(4.d0+ALFr)/gammaDP(1.d0+ALFr)
              ckQr1 = afr*gammaDP(1.d0+dmr+ALFr+bfr)/gammaDP(1.d0+dmr+ALFr)
              ckQr2 = afr*gammaDP(1.d0+bfr+ALFr)/gammaDP(1.d0+ALFr)
              ckQr3 = afr*gammaDP(7.d0+ALFr+bfr)/gammaDP(7.d0+ALFr)
              LAMr  = (cexr9*dble(NR(i,k)/(QR(i,k)*DE(i,k))))**thrd
         !The following calculations of VVX avoid over/underflow:
              VVQ(i,k)= -gamfact(i,k)*ckQr1*LAMr**(0.5*cexr2)/(LAMr+ffr)**(0.5*cexr1)    &
                         *LAMr**(0.5*cexr2)/(LAMr+ffr)**(0.5*cexr1)
              VqMax= max(VrMAX,-VVQ(i,k))
              if (scheme>1) then
                if(SS_ON) then
                 VVN(i,k)= -gamfact(i,k)*ckQr2*LAMr**(0.5*cexr4)/(LAMr+ffr)**(0.5*       &
                            cexr3)*LAMr**(0.5*cexr4)/(LAMr+ffr)**(0.5*cexr3)
                else
                  VVN(i,k) = VVQ(i,k)
                endif
                 VnMax= max(VrMAX,-VVN(i,k))
              endif
              if (scheme==4) then
                  if(SS_ON) then
                 VVZ(i,k)= -gamfact(i,k)*ckQr3*LAMr**(0.5*cexr6)/(LAMr+ffr)**(0.5*       &
                            cexr5)*LAMr**(0.5*cexr6)/(LAMr+ffr)**(0.5*cexr5)
                  else
                    VVZ(i,k) = VVQ(i,k)
                  end if
                 VzMax= max(VrMAX,-VVZ(i,k))
              endif
           endif
         enddo  !k-loop
       enddo    !i-loop
       locallim= (nnn==1)

       call blg5sedi(RHOQX,DZ,VVQ,nk,dtr,locallim,VqMax,FLIM,counter,activeColumn,       &
                     ktop_sedi)
       if (scheme >1)  &
          call blg5sedi(NR,DZ,VVN,nk,dtr,locallim,VnMax,FLIM,counter,activeColumn,       &
                        ktop_sedi)
       if (scheme==4)  &
          call blg5sedi(ZR,DZ,VVZ,nk,dtr,locallim,VzMax,FLIM,counter,activeColumn,       &
                        ktop_sedi)

        DO k = 1,nk-1                 ! Changed by WYH
          DO i = 1,ixe
            QR(i,k)= RHOQX(i,k)/DE(i,k)
          END DO
        END DO

    ! Prevent levels with zero N and nonzero Q and size-limiter:
       IF (scheme==2.or.scheme==3) THEN
         do a= 1,counter
           i=activeColumn(a)
            do k= 1,nk-1
              if (QR(i,k)>epsQ .and. NR(i,k)>epsN) then
                 Dx= ( dble(DE(i,k)*QR(i,k)/NR(i,k))*icmr)**thrd
                 ! Convert small raindrops to cloud droplets:
                 if (Dx<0.5d0*Dhh) then
                    QC(i,k)= QC(i,k)+QR(i,k)
                    NC(i,k)= NC(i,k)+NR(i,k)
                    QR(i,k)= 0.;  NR(i,k)= 0.;  Dr(i,k)= 0.
                 endif
                 ! Mean-drop size limiter:
                 if (Dx>3.d-3) then
                    tmp1= sngl(Dx)-3.e-3;  tmp2= tmp1*tmp1
                    tmp3= sngl(Dx/DrMAX);  tmp4= tmp3*tmp3*tmp3
                    NR(i,k)= NR(i,k)*max((1.+2.e4*tmp2),tmp4)
                 endif
              else
              ! Prevent levels with zero N and nonzero Q:
              if(clipCons_ON) then
                 Q(i,k) = Q(i,k) + QR(i,k)
                 T(i,k) = T(i,k) - QR(i,k)*LCP
              endif
                 QR(i,k)= 0.;  NR(i,k)= 0.
              endif
           enddo
         enddo
       ELSE IF (scheme==4) THEN
         do a= 1,counter
           i=activeColumn(a)
            do k= 1,nk-1
              if (QR(i,k)>epsQ .and. NR(i,k)>epsN .and. ZR(i,k)>epsZ) then
                 Dx= ( dble(DE(i,k)*QR(i,k)/NR(i,k))*icmr)**thrd
                 ! Convert small raindrops to cloud droplets:
                 if (Dx<0.5d0*Dhh) then
                    QC(i,k)= QC(i,k)+QR(i,k)
                    NC(i,k)= NC(i,k)+NR(i,k)
                    QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.;  Dr(i,k)= 0.
                 endif
                 ! Mean-drop size limiter:
                 if (Dx>3.d-3) then
                    tmp1= sngl(Dx)-3.e-3;  tmp2= tmp1*tmp1
                    tmp3= sngl(Dx/DrMAX);  tmp4= tmp3*tmp3*tmp3
                    NR(i,k)= NR(i,k)*max((1.+2.e4*tmp2),tmp4)
                 endif
              else
              ! Prevent levels with zero N and nonzero Q:
              if(clipCons_ON) then
                 Q(i,k) = Q(i,k) + QR(i,k)
                 T(i,k) = T(i,k) - QR(i,k)*LCP
              endif
                 QR(i,k)= 0.;  NR(i,k)= 0.;  ZR(i,k)= 0.
              endif
           enddo
         enddo
       ENDIF  !(if scheme>1)
       !LR(:)= LR(:) - cr6*VVQ(:,nk)*DE(:,nk)*QR(:,nk)
       LR(:)= LR(:) - cr6*VVQ(:,2)*DE(:,2)*QR(:,2)

     ENDDO  !nnn-loop

   ENDIF  !slabHASmass

!- - End of rain sedimentation - - - - - - - - - - - - - - - - - - - -  - - - - -!

!--  ICE  sedimentation:
  call SEDI_ISGH_v33(QI,NY,ZI,2,Q,T,DE,gamfact,epsQ,epsN,epsZ,afi,bfi,cmi,dmi,dti,ci6,  &
                 ALFifix,0.d0,LSP,npassi,ni,nk,ViMax,DiMax,DZ,SR,scheme,ktop_sedi,SS_ON,ixe)
!--  SNOW sedimentation:
  call SEDI_ISGH_v33(QN,NN,ZN,3,Q,T,DE,gamfact,epsQ,epsN,epsZ,afs,bfs,cms,dms,dts,cs6,  &
                 ALFsfix,Nosfix,LSP,npasss,ni,nk,VsMax,DsMax,DZ,SR,scheme,ktop_sedi,SS_ON,ixe)
!--  GRAUPEL sedimentation:
  call SEDI_ISGH_v33(QG,NG,ZG,4,Q,T,DE,gamfact,epsQ,epsN,epsZ,afg,bfg,cmg,dmg,dtg,cg6,  &
                 ALFgfix,Nogfix,LSP,npassg,ni,nk,VgMax,DgMax,DZ,SR,scheme,ktop_sedi,SS_ON,ixe)
!--  HAIL sedimentation:
  call SEDI_ISGH_v33(QH,NH,ZH,5,Q,T,DE,gamfact,epsQ,epsN,epsZ,afh,bfh,cmh,dmh,dth,ch6,  &
                 ALFhfix,Nohfix,LSP,npassh,ni,nk,VhMax,DhMax,DZ,SR,scheme,ktop_sedi,SS_ON,ixe)

!---  End of sedimentation for each category --------!

   LR= LR*ck7  !liquid precipitation rate
   SR= SR*ck7  !solid precipitation rate

 ENDIF  ! if (sedi_ON)

 where (Q<0.) Q= 0.

 !-----------------------------------------------------------------------------------!
 !                     End of sedimentation calculations (Part 4)                    !
 !-----------------------------------------------------------------------------------!

 !===================================================================================!
 !                             End of microphysics scheme                            !
 !===================================================================================!

!-------------
!Convert N from #/m3 to #/kg: Changed by DTD.  Previously only Z was being converted back by
!mistake.  Now both Z and N are converted back to "per mass units" properly for the dynamics
!DTD: UPDATE 09/27/2011.  This is now done in the model dynamics
!  IF (scheme >= 2) THEN
!  do k= 1,nk-1
!    do i=1,ixe
      !DE(:,k) = S(:,k)*PS(:)/(RGASD*T(:,k))     !air density at time  (t)
      !iDE(:,k)= 1./DE(:,k)
!      tmp1 = RGASD*T(i,k)/P(i,k)                       
!      tmp2 = RGASD*TM(i,k)/PM(i,k)
!      NCM(i,k)= NCM(i,k)*tmp1;   NC(i,k)= NC(i,k)*tmp2
!      NRM(i,k)= NRM(i,k)*tmp1;   NR(i,k)= NR(i,k)*tmp2
!      NYM(i,k)= NYM(i,k)*tmp1;   NY(i,k)= NY(i,k)*tmp2
!      NNM(i,k)= NNM(i,k)*tmp1;   NN(i,k)= NN(i,k)*tmp2
!      NGM(i,k)= NGM(i,k)*tmp1;   NG(i,k)= NG(i,k)*tmp2
!      NHM(i,k)= NHM(i,k)*tmp1;   NH(i,k)= NH(i,k)*tmp2
!      IF (scheme == 4) THEN
!        ZRM(i,k)= ZRM(i,k)*tmp1;   ZR(i,k)= ZR(i,k)*tmp2
!        ZIM(i,k)= ZIM(i,k)*tmp1;   ZI(i,k)= ZI(i,k)*tmp2
!        ZNM(i,k)= ZNM(i,k)*tmp1;   ZN(i,k)= ZN(i,k)*tmp2
!        ZGM(i,k)= ZGM(i,k)*tmp1;   ZG(i,k)= ZG(i,k)*tmp2
!        ZHM(i,k)= ZHM(i,k)*tmp1;   ZH(i,k)= ZH(i,k)*tmp2
!      END IF
!    enddo
!  enddo
!  END IF
!=============

 !-----------------------------------------------------------------------------------!
 !   Compute the tendencies of  T, Q, QC, etc. (to be passed back to model dynamics) !
 !   and reset the fields to their initial (saved) values at time {*}                !
 !-----------------------------------------------------------------------------------!
 ! DTD/WYH: Don't need the tendencies in ARPS

!      do k= 1,nk
!         do i= 1,ni
!
!            rtmp=T_TEND(i,k);   T_TEND(i,k)=(T(i,k) -T_TEND(i,k))*idt;  T(i,k) = rtmp
!            rtmp=Q_TEND(i,k);   Q_TEND(i,k)=(Q(i,k) -Q_TEND(i,k))*idt;  Q(i,k) = rtmp
!            rtmp=QCTEND(i,k);   QCTEND(i,k)=(QC(i,k)-QCTEND(i,k))*idt;  QC(i,k)= rtmp
!            rtmp=QRTEND(i,k);   QRTEND(i,k)=(QR(i,k)-QRTEND(i,k))*idt;  QR(i,k)= rtmp
!            rtmp=QITEND(i,k);   QITEND(i,k)=(QI(i,k)-QITEND(i,k))*idt;  QI(i,k)= rtmp
!            rtmp=QNTEND(i,k);   QNTEND(i,k)=(QN(i,k)-QNTEND(i,k))*idt;  QN(i,k)= rtmp
!            rtmp=QGTEND(i,k);   QGTEND(i,k)=(QG(i,k)-QGTEND(i,k))*idt;  QG(i,k)= rtmp
!            rtmp=QHTEND(i,k);   QHTEND(i,k)=(QH(i,k)-QHTEND(i,k))*idt;  QH(i,k)= rtmp

!            if (scheme>1) then
!              rtmp=NCTEND(i,k); NCTEND(i,k)=(NC(i,k)-NCTEND(i,k))*idt;  NC(i,k)= rtmp
!              rtmp=NRTEND(i,k); NRTEND(i,k)=(NR(i,k)-NRTEND(i,k))*idt;  NR(i,k)= rtmp
!              rtmp=NYTEND(i,k); NYTEND(i,k)=(NY(i,k)-NYTEND(i,k))*idt;  NY(i,k)= rtmp
!              rtmp=NNTEND(i,k); NNTEND(i,k)=(NN(i,k)-NNTEND(i,k))*idt;  NN(i,k)= rtmp
!              rtmp=NGTEND(i,k); NGTEND(i,k)=(NG(i,k)-NGTEND(i,k))*idt;  NG(i,k)= rtmp
!              rtmp=NHTEND(i,k); NHTEND(i,k)=(NH(i,k)-NHTEND(i,k))*idt;  NH(i,k)= rtmp
!            endif

!            if (scheme==4) then
!              rtmp=ZRTEND(i,k); ZRTEND(i,k)=(ZR(i,k)-ZRTEND(i,k))*idt;  ZR(i,k)= rtmp
!              rtmp=ZITEND(i,k); ZITEND(i,k)=(ZI(i,k)-ZITEND(i,k))*idt;  ZI(i,k)= rtmp
!              rtmp=ZNTEND(i,k); ZNTEND(i,k)=(ZN(i,k)-ZNTEND(i,k))*idt;  ZN(i,k)= rtmp
!              rtmp=ZGTEND(i,k); ZGTEND(i,k)=(ZG(i,k)-ZGTEND(i,k))*idt;  ZG(i,k)= rtmp
!              rtmp=ZHTEND(i,k); ZHTEND(i,k)=(ZH(i,k)-ZHTEND(i,k))*idt;  ZH(i,k)= rtmp
!            endif
!
!         enddo
!      enddo

 END SUBROUTINE MYTMOM_MAIN
!===================================================================================================!

END MODULE my_tmom_mod
