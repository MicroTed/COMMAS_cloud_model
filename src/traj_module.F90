!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\     MODULE TRAJ       ////////////////////
!     
!-----------------------------------------------------------------------------

  MODULE TRAJ_MODULE
  
!-----------------------------------------------------------------------------

  implicit none

  integer :: itraj = 0 ! 0 = no trajectories; 
                       ! 1 = write trajectories to netcdf only
                       ! 2 = write trajectories to both netcdf and text files
                       ! 3 = write trajectories to text files only
  integer :: maxtraj = 0 ! number of trajectories (caculated when itraj >= 1)
  integer :: ntraj   = 0 ! number of trajectories (caculated when itraj >= 1)
  integer :: ntrajtype = 1 ! 1 for parcel only, 2 to add rain trajectory (kessler or ZVD), 
                           ! 3 adds graupel+hail (ZVD), 4 adds rain+graupel + hail (ZVD)
  integer :: ixtrj1 = 1, ixtrj2 = 1 ! initial region for trajectories (global coordinates!)
  integer :: jytrj1 = 1, jytrj2 = 1
  integer :: kztrj1 = 1, kztrj2 = 1
  integer :: dxtraj = 1 ! initial spacing for trajectories (in grid points)
  integer :: dytraj = 1
  integer :: dztraj = 1
  integer :: time_traj1 = 0
  integer :: time_traj2 = 0
  real    :: riserate = 0.0
  integer :: icomtraj = 0 ! flag to tell us if comtraj is being run for a snapshot "frozen storm" sounding
  integer :: iverttraj = 0 ! flag to turn off horizontal motion of the balloon
  character(len=6) :: stimecomtraj
  integer, parameter :: n0 = 1 ! array start
  integer, parameter :: ninfomax = 76 ! = ieoffset + number of elec vars
  integer, parameter :: ieoffset = 48
  integer :: ninfo = ieoffset ! 47
  integer, parameter :: trj_unit = 181
  integer, parameter :: nmicrorates = 5
  integer, parameter :: nelecrates = 14
  real, allocatable :: trjdat(:,:,:)
  real, allocatable :: microrates(:,:,:,:)
  real, allocatable :: elecrates(:,:,:,:)
  integer :: numflashtraj
  integer :: ncidtraj(10)
  character(LEN=160) :: ncfiletraj(10)
  character(LEN=30) :: trajvarnames(ninfomax)
  
  integer :: chgavex = -1 ! number of horizontal points to go out for average charge, i.e., ic-chgavex to ic+chgavex
                          ! will be set to value of ng if left as negative
  integer :: chgavez = 0 ! number of vertical points to go out for average charge, i.e., kc-chgavez to kc+chgavez

 TYPE TRAJVARIABLE

  character(LEN=30)   name
  character(LEN=15)   unit
  character(LEN=255)  descr
  integer             varid
 
 END TYPE TRAJVARIABLE
 
  TYPE(TRAJVARIABLE), private :: trajvars(n0:ninfomax)

      integer, parameter, private :: ltt   = 1
      integer, parameter, private :: let   = 2
      integer, parameter, private :: lqt   = 3
      integer, parameter, private :: lct   = 4
      integer, parameter, private :: lrt   = 5
      integer, parameter, private :: lit   = 6
      integer, parameter, private :: lst   = 7
      integer, parameter, private :: lgt   = 8
      integer, parameter, private :: lut   = 9
      integer, parameter, private :: lvt   = 10
      integer, parameter, private :: lwt   = 11 
      integer, parameter, private :: lp    = 12
      integer, parameter, private :: lpp   = 13
      integer, parameter, private :: lx    = 14
      integer, parameter, private :: ly    = 15
      integer, parameter, private :: lz    = 16
      integer, parameter, private :: lthe  = 17
      integer, parameter, private :: lfrz  = 18
      integer, parameter, private :: lmlt  = 19
      integer, parameter, private :: ldep  = 20
      integer, parameter, private :: lsub  = 21
      integer, parameter, private :: lcnd  = 22
      integer, parameter, private :: levp  = 23
      integer, parameter, private :: lgmlt = 24
      integer, parameter, private :: lrevp = 25
      integer, parameter, private :: lbuoy = 26
      integer, parameter, private :: lpgrd = 27
      integer, parameter, private :: lw2   = 28
      integer, parameter, private :: lgdia = 29
      integer, parameter, private :: lcg   = 30
      integer, parameter, private :: lvtg  = 31
      integer, parameter, private :: load  = 32
      integer, parameter, private :: lres  = 33
      integer, parameter, private :: lht   = 34
      integer, parameter, private :: lvth  = 35 ! hail fall speed
!      integer, parameter, private :: loadh = 35
      integer, parameter, private :: lxn   = 36 !  'n'ext step x,y, and z
      integer, parameter, private :: lyn   = 37
      integer, parameter, private :: lzn   = 38
      integer, parameter, private :: lvtr  = 39 ! rain fall speed
      integer, parameter, private :: lrdia  = 40 ! rain mean diameter
      integer, parameter, private :: lhdia  = 41 ! hail mean diameter
      integer, parameter, private :: lztot = 42
      integer, parameter, private :: lzrt = 43
      integer, parameter, private :: lzht = 44
      integer, parameter, private :: lzhlt = 45
      integer, parameter, private :: lnoxt = 46
      integer, parameter, private :: lcot  = 47
      integer, parameter, private :: ldn   = 48 ! air density
      
      integer, parameter, private :: lex   = ieoffset + 1 ! 48 !  efield components at end of time step
      integer, parameter, private :: ley   = ieoffset + 2 ! 49
      integer, parameter, private :: lez   = ieoffset + 3 ! 50
      integer, parameter, private :: lemag = ieoffset + 4 ! 51 !  'n'ext step x,y, and z
      integer, parameter, private :: lphi      = ieoffset + 5 ! 52
      integer, parameter, private :: lnumflash = ieoffset + 6 ! 53
      integer, parameter, private :: lwtrise = ieoffset + 7 ! 54 ! w plus rise rate
      integer, parameter, private :: ldelex   = ieoffset + 8 ! 55 !  efield components
      integer, parameter, private :: ldeley   = ieoffset + 9 ! 56
      integer, parameter, private :: ldelez   = ieoffset + 10 ! 57
      integer, parameter, private :: ldelphi  = ieoffset + 11 ! 58
      integer, parameter, private :: lnetchg  = ieoffset + 12 ! 59
      integer, parameter, private :: lnetchgave  = ieoffset + 13 ! 60
      integer, parameter, private :: lexpl   = ieoffset + 14 ! 61 !  efield components at end of time step
      integer, parameter, private :: leypl   = ieoffset + 15 ! 62
      integer, parameter, private :: lezpl   = ieoffset + 16 ! 63
      integer, parameter, private :: lemagpl = ieoffset + 17 ! 64 !  'n'ext step x,y, and z
      integer, parameter, private :: lphipl  = ieoffset + 18 ! 65
      integer, parameter, private :: lrarh   = ieoffset + 19 ! 66
      integer, parameter, private :: lcrgis  = ieoffset + 20 ! 67
      integer, parameter, private :: lehw    = ieoffset + 21 ! 68
      integer, parameter, private :: lrarhl  = ieoffset + 22 ! 69
      integer, parameter, private :: lcrhis  = ieoffset + 23 ! 70
      integer, parameter, private :: lehlw   = ieoffset + 24 ! 71
      integer, parameter, private :: lchl    = ieoffset + 25 ! 72 ! hail conc.
      integer, parameter, private :: lcci    = ieoffset + 26 ! 73 ! ice crystal conc.
      integer, parameter, private :: lccw    = ieoffset + 27 ! 74 ! droplet conc.
      integer, parameter, private :: lcsw    = ieoffset + 28 ! 75 ! snow conc.
      
      
  
  integer, private :: timeid, trajstep = 0
  logical, private :: trajrestart
  
  CONTAINS

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! --------------------------------------------------------------------------------
   SUBROUTINE TRAJ_CLOSE(prefix,number)
! figure out how many trajectories and time steps and allocate space
   
   USE COMMASMPI_MODULE, only: my_rank,commasmpi_abort
   USE NETCDF
   implicit none
   
   character(*) :: prefix
   character(LEN=8) number
   integer i,status

   IF ( my_rank == 0 ) THEN
      DO i = 1,ntrajtype
      
      status = nf90_close(ncidtraj(i))
      
      ENDDO ! i : trajtype
   ENDIF

   RETURN
   END SUBROUTINE TRAJ_CLOSE
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE TRAJ_INIT(dt,nxend,nyend,nzend,nx,ny,nz,ng,prefix,number,start_commas)
! figure out how many trajectories and time steps and allocate space
   
   USE COMMASMPI_MODULE, only: my_rank,commasmpi_abort
   USE NETCDF
   USE GRID_MODULE, only: variable
   USE MICRO_MODULE, only: ipelec
   USE PARAM_MODULE, only: tke_type

   implicit none
   
   real :: dt
   integer :: nxend,nyend,nzend,nx,ny,nz,ng,start_commas
   character(*) :: prefix
   character(LEN=8) number
   
   integer :: nsteps ! not used at the moment
   integer :: ix,jy,kz,i,n
   character(LEN=1) trjnum
   character(len=2) s1,s2
   integer cmode, status, ncid, dim_id(0:8),dimvals(0:8)
   integer :: dimids(2), nid, trajnumid,start(1),count(1), record_id, records
   integer, allocatable :: trajects(:)
   
   integer :: it,nt

  real, allocatable :: tarray(:)
  logical, allocatable :: index(:)

      integer :: varid
      integer :: start1d(1),count1d(1)
      integer :: start2d(2),count2d(2)
      real :: data1d(1)
      real, allocatable :: data2d(:,:)
 
   nsteps = (time_traj2 - time_traj1)/dt + 1
   
   trajvars(ldn  )%name = 'RHO0' ; trajvars(ldn  )%unit = 'kg m-3' ; trajvars(ldn  )%descr = 'Air Density'
   trajvars(ltt  )%name = 'THETA' ; trajvars(ltt  )%unit = 'K' ; trajvars(ltt  )%descr = 'Potential Temperature'
   IF ( tke_type == 1 ) THEN
   trajvars(let  )%name = 'KM' ; trajvars(let  )%unit = 'm2 s-1' ; trajvars(let  )%descr = 'Mixing coefficient'
   ELSE
   trajvars(let  )%name = 'TKE' ; trajvars(let  )%unit = 'm s-1' ; trajvars(let  )%descr = 'Sqrt TKE'
   ENDIF
   trajvars(lqt  )%name = 'QV' ; trajvars(lqt  )%unit = 'kg kg-1' ; trajvars(lqt  )%descr = 'Water vapor mixing ratio'
   trajvars(lct  )%name = 'QC' ; trajvars(lct  )%unit = 'kg kg-1' ; trajvars(lct  )%descr = 'Cloud droplet mixing ratio'
   trajvars(lrt  )%name = 'QR' ; trajvars(lrt  )%unit = 'kg kg-1' ; trajvars(lrt  )%descr = 'Rain mixing ratio'
   trajvars(lit  )%name = 'QI' ; trajvars(lit  )%unit = 'kg kg-1' ; trajvars(lit  )%descr = 'Cloud ice mixing ratio'
   trajvars(lst  )%name = 'QS' ; trajvars(lst  )%unit = 'kg kg-1' ; trajvars(lst  )%descr = 'Snow mixing ratio'
   trajvars(lgt  )%name = 'QH' ; trajvars(lgt  )%unit = 'kg kg-1' ; trajvars(lgt  )%descr = 'Graupel mixing ratio'
   trajvars(lut  )%name = 'U' ; trajvars(lut  )%unit = 'm s-1' ; trajvars(lut  )%descr = 'u-wind'
   trajvars(lvt  )%name = 'V' ; trajvars(lvt  )%unit = 'm s-1' ; trajvars(lvt  )%descr = 'v-wind'
   trajvars(lwt  )%name = 'W' ; trajvars(lwt  )%unit = 'm s-1' ; trajvars(lwt  )%descr = 'w-wind'
   trajvars(lp   )%name = 'P0' ; trajvars(lp   )%unit = 'Pa' ; trajvars(lp   )%descr = 'Base state pressure'
   trajvars(lpp  )%name = 'PPERT' ; trajvars(lpp  )%unit = 'Pa' ; trajvars(lpp  )%descr = 'Perturbation pressure'
   trajvars(lx   )%name = 'X' ; trajvars(lx   )%unit = 'm' ; trajvars(lx   )%descr = 'x-position'
   trajvars(ly   )%name = 'Y' ; trajvars(ly   )%unit = 'm' ; trajvars(ly   )%descr = 'y-position'
   trajvars(lz   )%name = 'Z' ; trajvars(lz   )%unit = 'm' ; trajvars(lz   )%descr = 'z-position'
   trajvars(lthe )%name = 'THETAE' ; trajvars(lthe )%unit = '' ; trajvars(lthe )%descr = 'Equivalent potential temperature'
   trajvars(lfrz )%name = 'FRZ' ; trajvars(lfrz )%unit = 'J hr-1' ; trajvars(lfrz )%descr = 'Total freezing'
   trajvars(lmlt )%name = 'MLT' ; trajvars(lmlt )%unit = 'J hr-1' ; trajvars(lmlt )%descr = 'Total melting'
   trajvars(ldep )%name = 'DEP' ; trajvars(ldep )%unit = 'J hr-1' ; trajvars(ldep )%descr = 'Total deposition'
   trajvars(lsub )%name = 'SUB' ; trajvars(lsub )%unit = 'J hr-1' ; trajvars(lsub )%descr = 'Total sublimation'
   trajvars(lcnd )%name = 'COND' ; trajvars(lcnd )%unit = 'J hr-1' ; trajvars(lcnd )%descr = 'Total condensation'
   trajvars(levp )%name = 'EVAP' ; trajvars(levp )%unit = 'J hr-1' ; trajvars(levp )%descr = 'Total evaporation'
   trajvars(lgmlt)%name = 'GMLT' ; trajvars(lgmlt)%unit = 'J hr-1' ; trajvars(lgmlt)%descr = 'Graupel melting'
   trajvars(lrevp)%name = 'REVAP' ; trajvars(lrevp)%unit = 'J hr-1' ; trajvars(lrevp)%descr = 'Rain evaporation'
   trajvars(lbuoy)%name = 'BUOY' ; trajvars(lbuoy)%unit = 'm s-2' ; trajvars(lbuoy)%descr = 'Bouyancy'
   trajvars(lpgrd)%name = 'PGRD' ; trajvars(lpgrd)%unit = 'm s-2' ; trajvars(lpgrd)%descr = 'Vertical pressure gradient'
   trajvars(lw2  )%name = 'WVT'; trajvars(lw2  )%unit = 'm s-1' ; trajvars(lw2  )%descr = 'W - Vt'
   trajvars(lgdia)%name = 'HWDIA' ; trajvars(lgdia)%unit = 'm' ; trajvars(lgdia)%descr = 'graupel mean diameter'
   trajvars(lrdia)%name = 'RWDIA' ; trajvars(lrdia)%unit = 'm' ; trajvars(lrdia)%descr = 'rain mean diameter'
   trajvars(lhdia)%name = 'HLDIA' ; trajvars(lhdia)%unit = 'm' ; trajvars(lhdia)%descr = 'hail mean diameter'
   trajvars(lcg  )%name = 'CHW' ; trajvars(lcg  )%unit = 'm-3' ; trajvars(lcg  )%descr = 'graupel concentration'
   trajvars(lvtg )%name = 'VTH' ; trajvars(lvtg )%unit = 'm s-1' ; trajvars(lvtg )%descr = 'Graupel fall speed'
   trajvars(lvth )%name = 'VTHL' ; trajvars(lvth )%unit = 'm s-1' ; trajvars(lvth )%descr = 'Hail fall speed'
   trajvars(lvtr )%name = 'VTR' ; trajvars(lvtr )%unit = 'm s-1' ; trajvars(lvtr )%descr = 'Rain fall speed'
   trajvars(load )%name = 'LOAD' ; trajvars(load )%unit = 'm s-2' ; trajvars(load )%descr = 'mass loading'
   trajvars(lres )%name = 'T' ; trajvars(lres )%unit = 'K' ; trajvars(lres )%descr = 'Temperature'
   trajvars(lht  )%name = 'QHL' ; trajvars(lht  )%unit = 'kg kg-1' ; trajvars(lht  )%descr = 'Hail mixing ratio'
!   trajvars(loadh)%name = 'LOADH' ; trajvars(loadh)%unit = 'm s-2' ; trajvars(loadh)%descr = 'Hail loading'
   trajvars(lxn  )%name = 'XN' ; trajvars(lxn  )%unit = 'm' ; trajvars(lxn  )%descr = 'Next x-position'
   trajvars(lyn  )%name = 'YN' ; trajvars(lyn  )%unit = 'm' ; trajvars(lyn  )%descr = 'Next y-position'
   trajvars(lzn  )%name = 'ZN' ; trajvars(lzn  )%unit = 'm' ; trajvars(lzn  )%descr = 'Next z-position'
   trajvars(lztot)%name = 'DBZ' ; trajvars(lztot)%unit = 'dBZ' ; trajvars(lztot  )%descr = 'reflectivity'
   trajvars(lzrt)%name = 'ZRW' ; trajvars(lzrt)%unit = 'dBZ' ; trajvars(lzrt  )%descr = 'reflectivity'
   trajvars(lzht)%name = 'ZHW' ; trajvars(lzht)%unit = 'dBZ' ; trajvars(lzht  )%descr = 'reflectivity'
   trajvars(lzhlt)%name = 'ZHL' ; trajvars(lzhlt)%unit = 'dBZ' ; trajvars(lzhlt  )%descr = 'reflectivity'
   trajvars(lnoxt)%name = 'LNOX' ; trajvars(lnoxt)%unit = 'ppbv' ; trajvars(lnoxt  )%descr = 'Ambient NOx concentration'
   trajvars(lcot)%name = 'CO' ; trajvars(lcot)%unit = 'ppbv' ; trajvars(lcot  )%descr = 'CO concentration'
   
   trajvars(lex)%name = 'EX' ; trajvars(lex)%unit = 'V m-1' ; trajvars(lex  )%descr = 'X-comp Electric Field'
   trajvars(ley)%name = 'EY' ; trajvars(ley)%unit = 'V m-1' ; trajvars(ley  )%descr = 'Y-comp Electric Field'
   trajvars(lez)%name = 'EZ' ; trajvars(lez)%unit = 'V m-1' ; trajvars(lez  )%descr = 'Z-comp Electric Field'
   trajvars(lemag)%name = 'EMAG' ; trajvars(lemag)%unit = 'V m-1' ; trajvars(lemag  )%descr = 'Electric Field Magnitude'
   trajvars(lphi)%name = 'POTENTIAL' ; trajvars(lphi)%unit = 'V' ; trajvars(lphi  )%descr = 'Electric Potential'
   trajvars(lnumflash)%name = 'NUMFLASH' ; trajvars(lnumflash)%unit = 'x' ; trajvars(lnumflash  )%descr = 'number of flashes'
   trajvars(lwtrise)%name = 'RISERATE' ; trajvars(lwtrise)%unit = 'm s-1' ; trajvars(lwtrise  )%descr = 'Balloon rise rate'
   trajvars(ldelex)%name = 'Delta-EX' ; trajvars(ldelex)%unit = 'V m-1' ; trajvars(ldelex  )%descr = 'X-comp Electric Field Change'
   trajvars(ldeley)%name = 'Delta-EY' ; trajvars(ldeley)%unit = 'V m-1' ; trajvars(ldeley  )%descr = 'Y-comp Electric Field Change'
   trajvars(ldelez)%name = 'Delta-EZ' ; trajvars(ldelez)%unit = 'V m-1' ; trajvars(ldelez  )%descr = 'Z-comp Electric Field Change'
   trajvars(ldelphi)%name = 'Delta-POTENTIAL' ; trajvars(ldelphi)%unit = 'V'  
      trajvars(ldelphi  )%descr = 'Electric Potential Change'
   trajvars(lnetchg)%name = 'CHGNET' ; trajvars(lnetchg)%unit = 'nC m-3' 
      trajvars(lnetchg  )%descr = 'Net Electric Charge Density'
   trajvars(lnetchgave)%name = 'CHGNETAVE' ; trajvars(lnetchgave)%unit = 'nC m-3' 
      trajvars(lnetchgave  )%descr = 'Avg. Net Electric Charge Density'
   trajvars(lexpl)%name = 'EX-PreLightning' ; trajvars(lexpl)%unit = 'V m-1' 
      trajvars(lexpl  )%descr = 'X-comp Electric Field before lightning'
   trajvars(leypl)%name = 'EY-PreLightning' ; trajvars(leypl)%unit = 'V m-1' 
      trajvars(leypl  )%descr = 'Y-comp Electric Field before lightning'
   trajvars(lezpl)%name = 'EZ-PreLightning' ; trajvars(lezpl)%unit = 'V m-1' 
      trajvars(lezpl  )%descr = 'Z-comp Electric Field before lightning'
   trajvars(lemagpl)%name = 'EMAG-PreLightning' ; trajvars(lemagpl)%unit = 'V m-1' 
      trajvars(lemagpl  )%descr = 'Electric Field Magnitude before lightning'
   trajvars(lphipl)%name = 'POTENTIAL-PreLightning' ; trajvars(lphipl)%unit = 'V' 
      trajvars(lphipl  )%descr = 'Electric Potential before lightning'
   trajvars(lrarh)%name = 'RARH' ; trajvars(lrarh)%unit = 'g m-2 s-1' ; trajvars(lrarh  )%descr = 'Graupel rime acc. rate'
   trajvars(lcrgis)%name = 'RSCGIS' ; trajvars(lcrgis)%unit = 'C s-1' 
      trajvars(lcrgis  )%descr = 'Graupel noninductive charging rate'
   trajvars(lehw)%name = 'EHW' ; trajvars(lehw)%unit = 'X' ; trajvars(lehw  )%descr = 'Graupel-droplet collection efficiency'
   trajvars(lrarhl)%name = 'RARHL' ; trajvars(lrarhl)%unit = 'g m-2 s-1' ; trajvars(lrarhl  )%descr = 'Hail rime acc. rate'
   trajvars(lcrhis)%name = 'RSCHIS' ; trajvars(lcrhis)%unit = 'C s-1' 
      trajvars(lcrhis  )%descr = 'Hail noninductive charging rate'
   trajvars(lehlw)%name = 'EHLW' ; trajvars(lehlw)%unit = 'X' ; trajvars(lehlw  )%descr = 'Hail-droplet collection efficiency'
   trajvars(lchl )%name = 'CHL' ; trajvars(lchl  )%unit = 'm-3' ; trajvars(lchl  )%descr = 'Hail concentration'
   trajvars(lcci )%name = 'CCI' ; trajvars(lcci  )%unit = 'm-3' ; trajvars(lcci  )%descr = 'Ice Xtal concentration'
   trajvars(lccw )%name = 'CCW' ; trajvars(lccw  )%unit = 'm-3' ; trajvars(lccw  )%descr = 'Droplet concentration'
   trajvars(lcsw )%name = 'CSW' ; trajvars(lcsw  )%unit = 'm-3' ; trajvars(lcsw  )%descr = 'Snow concentration'


   ixtrj1 = Max(1,ixtrj1)
   jytrj1 = Max(1,jytrj1)
   kztrj1 = Max(1,kztrj1)
   
   ixtrj2 = Min(nxend-1, ixtrj2)
   jytrj2 = Min(nyend-1, jytrj2)
   kztrj2 = Min(nzend-1, kztrj2)

      maxtraj = 0
      do kz = kztrj1,kztrj2,dztraj
      do jy = jytrj1,jytrj2,dytraj
      do ix = ixtrj1,ixtrj2,dxtraj
        maxtraj = maxtraj+1
      end do
      end do
      end do
      
!      maxtraj = maxtraj+1
      allocate( trajects(maxtraj) )
      DO i = 1,maxtraj
       trajects(i) = i
      ENDDO
   
   IF ( ipelec > 0 ) THEN
     ninfo = ninfomax
   ELSE
     ninfo = ieoffset ! 48
   ENDIF
   
   allocate( trjdat(n0:ninfo,maxtraj,ntrajtype) ) ! maybe store nsteps at some point?
   trjdat(:,:,:) = 0.0
   
   allocate( microrates(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nmicrorates) )
   microrates(:,:,:,:) = 0.0

   IF ( ipelec .gt. 0 ) THEN
     allocate( elecrates(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,nelecrates) )
     elecrates(:,:,:,:) = 0.0
   ENDIF
   
   ! check if restarting and need to read in trajectories or starting a new set
   trajrestart = ( start_commas > time_traj1 .and. start_commas < time_traj2 )
   
   IF ( trajrestart .and.  itraj > 2 ) THEN
      write(0,*) 'Cannot restart trajectories unless itraj < 3! itraj = ',itraj
      write(0,*) 'Will start new trajectories instead'
      trajrestart = .false.
   ENDIF
   
   IF ( .not. trajrestart ) THEN
   
   IF ( my_rank == 0 .and. (itraj == 1 .or. itraj == 2) ) THEN
   cmode = ior(NF90_CLOBBER,NF90_CLASSIC_MODEL)
!   cmode = ior(NF90_CLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
      DO i = 1,ntrajtype
       write(trjnum,'(i1)') i
       IF ( icomtraj == 0 ) THEN
         ncfiletraj(i) = trim(prefix)//trim(number)//'.traj'//trjnum//'.nc'
       ELSEIF ( icomtraj == 1 ) THEN
         IF ( chgavex == ng .and. chgavez == 0 ) THEN
           ncfiletraj(i) = trim(prefix)//stimecomtraj//'.snapshot.traj'//trjnum//'.nc'
         ELSE
           write(s1,'(i2.2)') chgavex
           write(s2,'(i2.2)') chgavez

           ncfiletraj(i) = trim(prefix)//stimecomtraj//'.'//s1//'_'//s2//'.snapshot.traj'//trjnum//'.nc'
         ENDIF
       ELSEIF ( icomtraj == 2 ) THEN
         ncfiletraj(i) = trim(prefix)//stimecomtraj//'.snapshotvert.traj'//trjnum//'.nc'
       ENDIF
!       OPEN(unit=trj_unit+i-1, file=prefix(1:length)//trim(number)//'.traj'//trjnum, status='unknown',form='formatted')
        status = NF90_CREATE(ncfiletraj(i) ,cmode,ncidtraj(i))
        if (status /= nf90_noerr) then
         write(0,*) 'problem with creating netcdf trajectory file'
         call commasmpi_abort()
        endif
        ncid = ncidtraj(i)
        
      ! define dimensions for time (unlimited) and number of trajectories (maxtraj)
        
        status     = NF90_DEF_DIM(ncid, 'TRAJECTORY', maxtraj, dim_id(1))

        IF(status /= NF90_NOERR) print *,'TRAJ_INIT:  Error defining NTRAJ'
        status     = NF90_DEF_DIM(ncid, 'TIME', NF90_UNLIMITED, dim_id(0))
        IF(status /= NF90_NOERR) print *,'TRAJ_INIT:  Error defining OBS dimension'

     ! next define dimension variables.

        
        status = NF90_DEF_VAR(ncid, 'TIME', NF90_FLOAT, dim_id(0), nid)
        timeid = nid
        status = NF90_PUT_ATT(ncid, nid, "long_name", 'model time')
        status = NF90_PUT_ATT(ncid, nid, "units",     's')

        status = NF90_DEF_VAR(ncid, 'TRAJECTORY', NF90_INT, dim_id(1), nid)
        trajnumid = nid
        status = NF90_PUT_ATT(ncid, nid, "long_name", 'Trajectory number')
        status = NF90_PUT_ATT(ncid, nid, "units",     'count')


      dimids(1) = dim_id(1)
      dimids(2) = dim_id(0)
      
    ! define a file variable for each trajectory variable. The dimensions are maxtraj (fixed) and time (unlimited). 
    ! For example, 'X' will have the x position of each trajectory at a given time.
      DO n = n0,ninfo
       trajvars(n)%varid = 0
!       IF ( .not. ( n == lxn .or. n == lyn .or. n == lzn ) ) THEN
        
        status = NF90_DEF_VAR(ncid, trajvars(n)%name, NF90_FLOAT, dimids, trajvars(n)%varid)

  status = NF90_PUT_ATT(ncid, trajvars(n)%varid, "long_name", trajvars(n)%descr)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute long_name: ', trajvars(n)%name

  status = NF90_PUT_ATT(ncid, trajvars(n)%varid, "units",     trajvars(n)%unit)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute units:     ', trajvars(n)%name
        
!       ENDIF
      ENDDO
      
!-----------------------------------------------------------------------------------------------
! End definition section

  status = NF90_ENDDEF(NCID)
  IF(status /= NF90_NOERR) print *,'DEFINE_NETCDF:  Error ending define mode', NF90_STRERROR(status)

      start(1) = 1
      count(1) = maxtraj
      status = nf90_put_var(ncid, trajnumid, trajects, start, count )
      
      
      status = nf90_close(ncid)
      
      ENDDO ! i : trajtype
   
   ENDIF
   
   ELSE ! trajrestart
! read in the trajectory file to continue trajectories  
! need to figure out value of trajstep and read that time into the trajectory arrays

       allocate( data2d(maxtraj,1) )
       
       ntraj = maxtraj
 
       DO it = 1,ntrajtype
      
       write(trjnum,'(i1)') it
       ncfiletraj(it) = trim(prefix)//trim(number)//'.traj'//trjnum//'.nc'
      IF ( my_rank == 0 ) THEN
        write(0,*) 'Trying to open netcdf traj file ',ncfiletraj(it)
      ENDIF
      status = NF90_OPEN(ncfiletraj(it) ,NF90_NOWRITE,ncid)
      IF(status /= NF90_NOERR) write(6,*) 'TRAJ_INIT: ERROR OPENING FILE',ncfiletraj(it),NF90_STRERROR(status)

!-----------------------------------------------------------------------------------------------
! INQUIRE AS TO HOW MANY TIME STEPS HAVE BEEN WRITTEN 

  status = NF90_INQUIRE(ncid, unlimitedDimID = record_id)

  IF(status /= NF90_NOERR) write(6,*) 'TRAJ_INIT: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
  
  status = NF90_INQUIRE_DIMENSION(ncid, record_id, len = records)
  
  IF(status /= NF90_NOERR) write(6,*) 'TRAJ_INIT: ERROR READING TIME DIM ',NF90_STRERROR(status)

    allocate(tarray(records))

    status = NF90_INQ_VARID(ncid,'TIME',timeid)
  IF(status /= NF90_NOERR) write(6,*) 'TRAJ_INIT: ERROR FINDING TIME VARIABLE ',NF90_STRERROR(status)

    status = NF90_GET_VAR(ncid, timeid, tarray)
    IF(status /= NF90_NOERR) print *, 'TRAJ_INIT:  ERROR READING COORDINATE TIME ',NF90_STRERROR(status)


   write(0,*) 'record_id, timeid, records,start_commas = ',record_id,timeid,records,start_commas
   trajstep = -1
   DO i = 1,records
!     write(0,*) 'i,tarray = ',i,tarray(i)
     IF ( start_commas == tarray(i) ) THEN
       trajstep = i
       exit
     ENDIF
   ENDDO
   
   deallocate(tarray)

   
   IF ( trajstep <  0 ) THEN
     write(0,*)  'TRAJ_INIT:  ERROR FINDING  trajstep'
     stop
   ELSE
     write(0,*) 'TRAJ_INIT:   found  trajstep, start_commas = ',trajstep, start_commas
   ENDIF
   
! need to get time id?
      start1d(1) = trajstep
      count1d(1) = 1
      data1d(1)  = start_commas
      status = nf90_get_var(ncid, timeid, data1d, start1d, count1d )
      start2d(1) = trajstep
!      write(0,*) 'TRAJ: trajstep,time_real = ',trajstep,time_real

      DO nt = n0,ninfo
!       write(0,*) 'nt = ',nt,lxn
!      DO nt = 1,1
!       IF ( .not. ( nt == lxn .or. nt == lyn .or. nt == lzn ) ) THEN

      
        
      start2d(1) = 1
      start2d(2) = trajstep
      count2d(1) = ntraj
      count2d(2) = 1

! need to get varid first

        status = NF90_INQ_VARID(ncid,trajvars(nt)%name, trajvars(nt)%varid) 

        status = nf90_get_var(ncid, trajvars(nt)%varid, data2d, start=start2d, count=count2d )
        IF(status /= NF90_NOERR) THEN 
          write(0,*) 'NF90_PUT_VAR:  Error writing variable: ', trajvars(nt)%name
          write(0,*) 'start2d = ',start2d(1),start2d(2)
          write(0,*) 'count2d = ',count2d(1),count2d(2)
          write(0,*) 'start2d = ',start2d(1),start2d(2)
          call commasmpi_abort()
        ENDIF

! changed this to invert the loops (outer loop on ninfo, inner on ntraj) and write all the current values at once
! could try it with the dimensions switched back and use a stride value so that nf90_put_var is still called only 
! once per variable (ninfo times) instead of the original inefficient ninfo*ntraj calls, but writing contiguous
! data is probably a good bit faster. It is also easy enough to transpose the arrays after reading them to get into
! trajectory space
      DO n = 1,ntraj
        trjdat(nt,n,it) = data2d(n,1)
      ENDDO
      
!      IF ( my_rank == 0 .and. nt == lxn ) THEN
        DO n = 1,ntraj
!         write(0,*) trajvars(nt)%name, ' for nt,n,it = ',nt,n,it, trjdat(nt,n,it), data2d(n,1)
        ENDDO
!      ENDIF
      
!      ENDIF
      ENDDO ! nt

      status = nf90_close(ncid)
      IF ( my_rank == 0 ) THEN
        write(0,*) 'Done reading trajectory file ',ncfiletraj(it)
      ENDIF
      
     ENDDO ! it

   ENDIF


   deallocate(trajects)
   
   RETURN
   END SUBROUTINE TRAJ_INIT
   
!-----------------------------------------------------------------------------
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE TRAJ       ////////////////////
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE TRAJ(nx,ny,nz,na,dt,an,ab,u,v,w,t0,t8,t9,pinit,  &
                     p2,km,dbz,elec,         &
                     gxt,gyt,gzt,time,time_real,    &
                     uinit,vinit,ugrid,vgrid,microp,dx,dy)

   USE GRID_MODULE, only: variable
   USE PARAM_MODULE, only: ng, g, cv, cp, rcp
   USE GRIDPARAM_MODULE, only: xdomain,ydomain,zdomain
   USE COMMASMPI_MODULE
   USE INDEX_MODULE !, only: neelec,ntakpd,lt,lv,lvi,lvs,lvh,lvhl,lc,lr,li,ls,lh,lhl, &
                    !       lzi,lzr,lzh,lzhl,lnc,lnr,lni,lns,lnh,lnhl
   USE MICRO_MODULE, only: rho_qs, rho_qh, rho_qhl, itype1,itype2,   &
     &                     ipconc, ccn, cimn, cimx, imurain, cxmin, zxmin, &
     &                     takcxmin,ipelec
   USE NETCDF
   use takcommon,       only: lmax,lfmax,kfmax,iimax,kimax,lr100,ls250,shedsmall,dj, &
     &                        bka,capw,bmu,hlat
   use comm1,           only: rhof,acomm1
   use comm2,           only: sr,sx,srf,vw,bew,szr,szf,szi      &
     &                            ,sri,sfi,di,sxf,sxi,sh,ff,vf,fi,vi   
!-----------------------------------------------------------------------------
   implicit none

   integer :: nx,ny,nz,na
   real    :: dt
   real    :: an(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na)
   real    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: t8(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: t9(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   real    :: uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng)
   real    :: ugrid, vgrid
   real    :: pinit(-ng+1:nz+ng)  ! base state Pi
   real    :: p2(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)  ! perturbation Pi
   real    :: km(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)  ! subgrid mixing
   real    :: dbz(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: ab(-ng+1:nz+ng,na)
   integer :: loop
   integer :: time
   real    :: time_real

   TYPE(VARIABLE)  :: elec(neelec)
   
   integer :: n,nt
   integer :: i,j,k
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   integer :: nnum
   character(len=*) :: microp
   real    :: dx,dy
   
   real    :: xmn,xmx
   integer, save :: initflag = 0
   integer :: in, it, il
   integer :: ix,jy,kz
   integer :: ic,jc,kc
   integer :: k1,k2,x1,x2,y1,y2
   integer :: iu,ju,ku
   integer :: iv,jv,kv
   integer :: iw,jw,kw
   real    :: x,y,z
   real    :: facx,facy,facz
   real    :: facxu,facyu,faczu
   real    :: facxv,facyv,faczv
   real    :: facxw,facyw,faczw
   real    :: uint,vint,wint
   integer :: ntraj_tile
   logical :: lcheck
   real, parameter :: eps = 1.e-08
   integer, parameter :: idebug = 0
   real, allocatable :: dn(:,:,:), pn(:,:,:), chgnet(:,:,:)
   real    :: pb(-ng+1:nz+ng)
   real    :: db(-ng+1:nz+ng)
   real    :: qsum, dumint, exint, eyint, ezint
   real    :: pres, tc, qvc
   real    :: thetae
   real    :: rcgs, vtden, qr, vt
   real    :: chw, z1, rdi, alp
   real, parameter ::  pi = 3.14159265 ! 4.0*atan(1.0)
   

!   integer,allocatable :: loc(3,maxtraj)
   real,allocatable,save :: trjbuf(:,:)
   real,allocatable,save :: trjloc(:,:)


      
! local automatic arrays

      integer, parameter :: ngs = 1
      integer, parameter :: ngscnt = 1
      integer :: mgs
      real :: rho0(ngs), rhovt(ngs), temcg(ngs)
      real :: cwnccn(ngs), fadvisc(ngs)
      real :: cnina(ngs)
      integer :: igs(ngs),kgs(ngs)
      real :: qx(ngscnt,lv:lhab)
      real :: qxw(ngscnt,ls:lhab)
      real :: cx(ngscnt,lc:lhab)
      real :: cxtmp(ngscnt,lc:lhab)
      real :: xv(ngscnt,lc:lhab)
      real :: vtxbar(ngscnt,lc:lhab,3)
      real :: xmas(ngscnt,lc:lhab)
      real :: xdn(ngscnt,lc:lhab)
      real :: cdxgs(ngs,lc:lhab)
      real :: xdia(ngscnt,lc:lhab,3) 
      real :: vx(ngscnt,li:lhab)
      real :: alpha(ngscnt,lc:lhab)
      real :: alphan(ngscnt,lc:lhab)
      real :: zx(ngscnt,lr:lhab)
      real xdnmx(lc:lhab), xdnmn(lc:lhab), xdn0(lc:lhab)
      real qxmin(lc:lhab)

      real cdx(lc:lhab)
      real cno(lc:lhab), cnox
      real xvmn(lc:lhab), xvmx(lc:lhab)

      logical ldovol, ldoliq, ldoz
      integer lvol(lc:lhab)
      integer ln(lc:lhab)
      integer lzx(lc:lhab)
      integer lliq(li:lhab)

      real :: axx(ngs,lh:lhab), bxx(ngs,lh:lhab)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)
      
      real :: denom, numer, qmin 
      real :: a

      integer :: ncid, varid, status
      integer :: start1d(1),count1d(1)
      integer :: start2d(2),count2d(2)
      real :: data1d(1)
      real, allocatable,save :: data2d(:,:)
      
      double precision :: term,term3,term4,term5

      integer nbin,l,kf,ia
      parameter (nbin=ntakpd)  ! number of mass bins for bin model
      real rn(nbin) !,rd(nbin),rm(nbin)
      real rq(nbin),vtr(nbin) !,rdrd(nbin)
      real, save :: rda(nbin,3),rma(nbin,3)=0.0,rva(nbin,3),rvna(nbin,3),rdrda(nbin,3)
      real crbin(nbin),chbin(nbin,2),chlbin(nbin),cibin(ntakid,ntakit)
      double precision :: totn, totvt, totq, totz
      integer, save :: lf75,lf500,lf150


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      integer :: downward_tag, upward_tag

      logical :: debug_mpi = .false.
      integer       :: nampi, nb, proc
#ifdef MPI
      INCLUDE "mpif.h"
      INTEGER :: mpi_status(MPI_Status_size)
#endif
!-----------------------------------------------------------------------------

   work_to_do = .false.
    IF ( time_real .ge. time_traj1 .and. time_real .le. time_traj2 ) THEN
      work_to_do = .true.
      trajstep = trajstep + 1
    ELSE
      RETURN
    ENDIF
   IF ( idebug .ge. 1 ) write(0,*) 'traj: time_real, time_traj1, time_traj2, micro = ', &
      time_real, time_traj1, time_traj2, microp

     IF ( microp(1:1) == 'Z' ) THEN
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

      ln(:) = 0
      ln(lc) = lnc
      ln(lr) = lnr
      IF ( li > 0 ) THEN
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      ENDIF
      IF ( lhl .gt. 1 ) ln(lhl) = lnhl

      lzx(:) = 0
      lzx(lr) = lzr
      IF ( li > 0 ) THEN
      lzx(li) = lzi
      lzx(ls) = lzs
      lzx(lh) = lzh
      ENDIF
      IF ( lhl .gt. 1 .and. lzhl > 1 ) lzx(lhl) = lzhl

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

      CALL setqxminz(qxmin)

      xvmn(lc) = xvcmn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn
      xvmn(lh) = xvhmn0

      xvmx(lc) = xvcmx
      xvmx(lr) = xvrmx0
      xvmx(ls) = xvsmx0
      xvmx(lh) = xvhmx0
      
      IF ( lhl .gt. 1 ) THEN
      xvmn(lhl) = xvhlmn0
      xvmx(lhl) = xvhlmx0
      ENDIF
      CALL setcnoz(cno)
      
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
      xdn0(ls) = rho_qs ! 100.0
      xdn0(lh) = rho_qh ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      IF ( lhl .gt. 1 ) xdn0(lhl) = rho_qhl ! 800.0

!
!  Set terminal velocities...
!    also set drag coefficients
!
      cdx(lr) = 0.60
      cdx(lh) = 0.8 ! 1.0 ! 0.45
      cdx(ls) = 2.00
!      cd(1) = cdx(ls)
      IF ( lhl .gt. 1 ) cdx(lhl) = 0.45

!      cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
!      cwradn = 5.0e-6     ! minimum radius
!      cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
      cnina(1) = 1.0


      ENDIF ! microp = Z
      
      IF ( microp(1:3) == 'TAK' ) THEN
        IF ( sr(1) == 0.d0 ) THEN
        call takinitbin
!        write(0,*) 'BINFORCE call takinitbin, sx(1) = ',sx(1),sxf(1,1),sxf(1,2),lmax,lfmax

        ENDIF
        IF ( rma(1,1) == 0.0 ) THEN
          lf75 = 1
          lf500 = 1
          lf150 = 1
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0075 ) lf75 = l
            IF ( srf(l,2) <= 0.0150 ) lf150 = l
            IF ( srf(l,2) <= 0.0500 ) lf500 = l
          enddo
        DO kf=1,2
        DO l = 1,ntakpd
         rma(l,kf)  = 1.e-3*sxf(l,kf) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
!         rva(l,kf)  = rma(l,kf)/(xden)  ! (m^3) volume (mass/1000.)
         rda(l,kf)  = 2.*1.e-2*srf(l,kf) ! (cm to meters and radius to diameter) (6.*rma(l,kf)/(pi*xdn))**(1./3.)
!         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
!         rvna(l,kf) = 1.
!         rdrda(l,kf) = rda(l,kf)/dj
!         write(0,*) 'TRAJ:l,rma,rda,sr2,rdrda = ',l,rma(l,kf),rda(l,kf),2.*1.e-2*sr(l),rdrda(l,kf)
!         write(0,*) 'TRAJ:l,rma,rda,vf = ',l,rma(l,kf),rda(l,kf),vf(l,kf)
        ENDDO
        ENDDO
        
        DO l = 1,ntakrd
         rma(l,3)  = 1.e-3*sx(l) ! (g to kg)  hmmin*exp(3.0*(l-1)/hjo)
!         rva(l,kf)  = rma(l,kf)/(xden)  ! (m^3) volume (mass/1000.)
         rda(l,3)  = 2.*1.e-2*sr(l) ! (cm to meters and radius to diameter) (6.*rma(l,kf)/(pi*xdn))**(1./3.)
!         rda(l)  = (6.*rma(l)/(pii*xden))**(1./3.)
!         rdamelt(l)  = (6.*rma(l)/(pii*1000.))**(1./3.)
!         rvna(l,3) = 1.
!         rdrda(l,3) = rda(l,3)/dj
!         write(0,*) 'TRAJ:l,rma,rda,sr2,rdrda = ',l,rma(l,kf),rda(l,kf),2.*1.e-2*sr(l),rdrda(l,kf)
!         write(0,*) 'TRAJ:l,rma,rda,vf = ',l,rma(l,kf),rda(l,kf),vf(l,kf)
        ENDDO
        
        ENDIF
      ENDIF
      

   

      if ( initflag == 0 ) then
      initflag = 1
      in = 0
      do kz = kztrj1,kztrj2,dztraj
      do jy = jytrj1,jytrj2,dytraj
      do ix = ixtrj1,ixtrj2,dxtraj
      in = in+1
       IF ( in > maxtraj ) THEN
        IF ( my_rank == 0 ) THEN
          write(0,*) 'TRAJ: Somebody screwed up the number of trajectories!'
          write(0,*) 'maxtraj,in = ',maxtraj,in
        ENDIF
        CALL COMMASMPI_ABORT()
       ENDIF
      IF ( .not. trajrestart ) THEN
      trjdat(lx,in,1:ntrajtype) = (ix-1)*dx + 0.5*dx
      trjdat(ly,in,1:ntrajtype) = (jy-1)*dy + 0.5*dy
      trjdat(lz,in,1:ntrajtype) = gzt(kz,1)
      trjdat(lxn,in,1:ntrajtype) = (ix-1)*dx + 0.5*dx
      trjdat(lyn,in,1:ntrajtype) = (jy-1)*dy + 0.5*dy
      trjdat(lzn,in,1:ntrajtype) = gzt(kz,1)
!       IF ( idebug > 0 ) THEN
!        write(0,*) 'in,x,y,z = ',in,(ix-1)*dx + 0.5*dx,(jy-1)*dy + 0.5*dx,gzt(kz,1)
!       ENDIF

      ELSE
       IF ( my_rank == 0 ) THEN
!        write(0,*) 'in, xn = ',in,trjdat(lxn,in,1) 
!        write(0,*) 'in, yn = ',in,trjdat(lyn,in,1) 
!        write(0,*) 'in, zn = ',in,trjdat(lzn,in,1) 
       ENDIF
      ENDIF
      
      end do
      end do
      end do
      
      ntraj = in

        IF ( my_rank == 0 ) THEN
          write(0,*) 'maxtraj,ntraj = ',maxtraj,ntraj
        ENDIF
        
        IF ( .not. allocated( trjbuf ) ) THEN
          allocate( trjbuf(n0:ninfo+2,ntraj*ntrajtype) )
          allocate( trjloc(3,ntraj*ntrajtype) )
          allocate( data2d(maxtraj,1) )
        ENDIF
      
      ELSE
        ! transfer future values to current values
        DO it = 1,ntrajtype
        DO n = 1,ntraj
          trjdat(lx,n,it) = trjdat(lxn,n,it)
          trjdat(ly,n,it) = trjdat(lyn,n,it)
          trjdat(lz,n,it) = trjdat(lzn,n,it)

          trjdat(lqt:lx-1,n,it) = 0 ! reset other values
          trjdat(lfrz:lxn-1,n,it) = 0 ! reset other values
          IF ( ninfo > lzn ) trjdat(lzn+1:ninfo,n,it) = 0 ! reset other values
          
!          trjdat(1:lx-1,n,it) = 0 ! reset other values
!          trjdat(lz+1:ninfo-3,n,it) = 0 ! reset other values
        ENDDO
        ENDDO
      
      end if

      allocate ( dn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )
      allocate ( pn(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) )

      ixb = 0
      ixe = itile+1
      if (ixbeg .eq. nxbeg) ixb = 1
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      jyb = 0
      jye = jtile+1
      if (jybeg .eq. nybeg) jyb = 1
      if (jyend .eq. nyend) jye = jyend-jybeg

      kzb = -1
      kze = ktile+2
      IF ( kzbeg == nzbeg ) kzb = 1
      if (kzend .eq. nzend) kze = kzend-kzbeg

      DO kz=kzb,kze
       db(kz)= 1.0e5*pinit(kz)**2.509/(287.04*ab(kz,lt))
       pb(kz)= 1.0e5*pinit(kz)**3.509


        DO ix=ixb,ixe
         DO jy=jyb,jye
            dn(ix,jy,kz)= 1.0e5*(pinit(kz)+p2(ix,jy,kz))**2.509/(287.04*an(ix,jy,kz,lt) )
            pn(ix,jy,kz)= 1.0e5*(pinit(kz)+p2(ix,jy,kz))**3.509 - pb(kz)
            IF ( kz >= 2 ) THEN
             t9(ix,jy,kz) = -2.0/(dn(ix,jy,kz)+dn(ix,jy,kz-1)) &
                          *gzt(kz,3)*(pn(ix,jy,kz)-pn(ix,jy,kz-1))
             IF ( microp(1:3) == 'LFO' .or. microp(1:1) == 'Z' .or. microp(1:3) == 'TAK') THEN ! mixing ratios are diagnosed in TAK
              qsum = 0.5*(an(ix,jy,kz,lc)+an(ix,jy,kz,lr)+an(ix,jy,kz,li)      &
                       + an(ix,jy,kz,ls)+an(ix,jy,kz,lh) +                    &
                         an(ix,jy,kz-1,lc)+an(ix,jy,kz-1,lr)+an(ix,jy,kz-1,li) &
                       + an(ix,jy,kz-1,ls)+an(ix,jy,kz-1,lh) )
               IF ( lhl > 1 ) THEN 
                 qsum = qsum + 0.5*(an(ix,jy,kz,lhl)+ an(ix,jy,kz-1,lhl) )
               ENDIF
             ELSEIF ( microp(1:1) == 'K'  ) THEN
             qsum = 0.5*(an(ix,jy,kz,lc)+an(ix,jy,kz,lr) +  &
                         an(ix,jy,kz-1,lc)+an(ix,jy,kz-1,lr) )
             ENDIF
             
             t8(ix,jy,kz) = g*(     &
               0.5*((an(ix,jy,kz,lt)-ab(kz,lt))     &
                    /ab(kz,lt)     &
                   +(an(ix,jy,kz-1,lt)-ab(kz-1,lt))     &
                    /ab(kz-1,lt))     &
             + 0.5*0.608*((an(ix,jy,kz,lv)-ab(kz,lv))     &
                         +(an(ix,jy,kz-1,lv)-ab(kz-1,lv)))     &
             - cv/cp*(0.5*(pn(ix,jy,kz)/pb(kz)     &
                          +pn(ix,jy,kz-1)/pb(kz-1)))     &
             - qsum   )
            ELSE 
             t8(ix,jy,kz) = 0.0
             t9(ix,jy,kz) = 0.0
            ENDIF
            
          ENDDO
        ENDDO
      ENDDO

      IF ( ipelec > 0 ) THEN
        allocate ( chgnet(-chgavex+1:nx+chgavex,-chgavex+1:ny+chgavex,-ng+1:nz+ng) )
        chgnet(:,:,:) = 0.0
        DO kz=kzb,kze
          DO jy=jyb,jye
            DO ix=ixb,ixe
             chgnet(ix,jy,kz) = elec(iscnet)%flt3d(ix,jy,kz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

   !  Fill in values at k = 0 and k = 1
       DO ix=ixb,ixe
         DO jy=jyb,jye
             t8(ix,jy,0:1) = t8(ix,jy,2)
             t9(ix,jy,0:1) = t9(ix,jy,2)
        ENDDO
      ENDDO
      

#ifdef MPI
! communicate microrates by one grid cell.
      IF ( number_of_processes .gt. 1 ) THEN

       nb = 1

       IF ( allocated( microrates ) ) THEN

       nampi = nmicrorates

       IF ( nproci > 1 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,microrates )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,microrates)
       ENDIF
       
       IF ( nprocj > 1 ) THEN
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,microrates)

       northward_tag = 10004
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,microrates)
       ENDIF

       IF ( nprock > 1 ) THEN
       
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,microrates)

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,microrates)

       ENDIF
       ENDIF
      


       IF ( allocated( elecrates ) ) THEN

       nampi = nelecrates

       IF ( nproci > 1 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,elecrates )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,elecrates)
       ENDIF
       
       IF ( nprocj > 1 ) THEN
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,elecrates)

       northward_tag = 10014
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,elecrates)
       ENDIF

       IF ( nprock > 1 ) THEN
       
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,elecrates)

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,elecrates)

       ENDIF
       ENDIF
      
       IF ( ipelec > 0 ) THEN

       nampi = 1

       IF ( nproci > 1 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,elec(iscnet)%flt3d )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,elec(iscnet)%flt3d)
       ENDIF
       
       IF ( nprocj > 1 ) THEN
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,elec(iscnet)%flt3d)

       northward_tag = 10024
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,elec(iscnet)%flt3d)
       ENDIF

       IF ( nprock > 1 ) THEN
       
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,elec(iscnet)%flt3d)

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,elec(iscnet)%flt3d)

       ENDIF

       nampi = 1
       nb = chgavex ! number of ghost zones to communicate


       IF ( nproci > 1 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,chgavex,chgavex,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,chgnet )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,chgavex,chgavex,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,chgnet)
       ENDIF
       
       IF ( nprocj > 1 ) THEN
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,chgavex,chgavex,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,chgnet)

       northward_tag = 10034
       CALL sendrecv_northward(nx,ny,nz,chgavex,chgavex,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,chgnet)
       ENDIF

       ENDIF ! ipelec
      
      
      ENDIF

#endif

    t0(:,:,:) = 0.0
!    loc(:,:) = 0

!      facz = (z - gzt(kc,1))*gzt(kc,3)
! kc is grid level below z
      IF (  ipelec > 0 .and. allocated( elecrates ) ) THEN
      ! for chgnetave, first put the local disk averages into t0 for later interpolation
      term = 0.
      term3 = 0.
      term4 = 0.
      term5 = 0.
      DO k = 1,nz-1
      DO j = 1,ny
      DO i = 1,nx
      
         x1 = Max( -ng, i-chgavex )
         x2 = Min( nx+ng, i+chgavex )
         y1 = Max( -ng, j-chgavex )
         y2 = Min( ny+ng, j+chgavex )
         k1 = k
        ! k2 = k + 1
         kz = 0
         dumint = 0.0
!         term4 = Max(term4,elec(iscnet)%flt3d(i,j,k) )
!         term5 = Min(term5,elec(iscnet)%flt3d(i,j,k) )
         term4 = Max(term4,chgnet(i,j,k) )
         term5 = Min(term5,chgnet(i,j,k) )

!          IF ( k == nz/2 .and. i == nx/2 .and. j == ny/2 ) THEN
!            write(0,*) 'i,j,x1,x2,y1,y2,chgavex = ',i,j,x1,x2,y1,y2,chgavex
!          ENDIF
         DO ix = x1,x2
         DO jy = y1,y2
!          IF ( k == nz/2 .and. i == nx/2 .and. j == ny/2 ) THEN
!            write(0,*) 'ix,jy,diff = ',ix,jy,(jy-j)**2 + (ix-i)**2 - chgavex**2 
!          ENDIF
           IF ( (jy-j)**2 + (ix-i)**2 <= chgavex**2 ) THEN ! approximate a circle
            kz = kz + 1
!            dumint = dumint + elec(iscnet)%flt3d(ix,jy,k) ! +  &
            dumint = dumint + chgnet(ix,jy,k) ! +  &
             !  facz*(elec(iscnet)%flt3d(ix,jy,k2) - elec(iscnet)%flt3d(ix,jy,k1))
           ENDIF
         ENDDO
         ENDDO
         
         IF ( kz > 0 ) THEN
         term = max(term,1.e9*dumint/float(kz) )
         term3 = min(term3,1.e9*dumint/float(kz) )
         t0(i,j,k) = 1.e9*dumint/float(kz)
         ELSE
!         write(
         ENDIF
       
       ENDDO
       ENDDO
       ENDDO

!         write(0,*) 'traj: max/min t0,scnet = ',term,term3,term4,term5
#ifdef MPI
! communicate microrates by one grid cell.
      IF ( number_of_processes .gt. 1 ) THEN

       nb = 1
       nampi = 1

       IF ( nproci > 1 ) THEN
       westward_tag = 10001
       CALL sendrecv_westward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),westward_tag,t0 )

       eastward_tag = 10002
       CALL sendrecv_eastward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      w_proc(my_rank),e_proc(my_rank),eastward_tag,t0)
       ENDIF
       
       IF ( nprocj > 1 ) THEN
       southward_tag = 10003
       CALL sendrecv_southward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),southward_tag,t0)

       northward_tag = 10044
       CALL sendrecv_northward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      n_proc(my_rank),s_proc(my_rank),northward_tag,t0)
       ENDIF

       IF ( nprock > 1 ) THEN
       
       downward_tag = 10005
       CALL sendrecv_downward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),downward_tag,t0)

       upward_tag = 10006
       CALL sendrecv_upward(nx,ny,nz,ng,ng,ng,nb,nampi,   & 
     &      d_proc(my_rank),u_proc(my_rank),upward_tag,t0)

       ENDIF
      
      
      ENDIF
#endif
       ENDIF
     
     ntraj_tile = 0
     
     
        ixe = itile+1
        if (ixend .eq. nxend) ixe = ixend-ixbeg
        jye = jtile+1
        if (jyend .eq. nyend) jye = jyend-jybeg
        kze = ktile+1
        if (kzend .eq. nzend) kze = kzend-kzbeg

!     write(*,*) 'TRAJ: ixe,jye, gxt = ',ixe,jye,gxt(1,1),gxt(ixe,1)
!     write(*,*) 'ixbeg,ixend,nxend = ',ixbeg,ixend,nxend
     
     DO it = 1,ntrajtype
     DO n = 1,ntraj
    ! check if parcel is in this tile
    ! Use the SCALAR points because we do not go 'left' of ic = 1
      lcheck = ( trjdat(lx,n,it) >= gxt(1,1) .and. trjdat(lx,n,it) < gxt(ixe,1) ) .and. &
               ( trjdat(ly,n,it) >= gyt(1,1) .and. trjdat(ly,n,it) < gyt(jye,1) ) .and. &
               ( trjdat(lz,n,it) >= gzt(1,1) .and. trjdat(lz,n,it) < gzt(kze,1) )
      IF ( lcheck ) THEN
       ntraj_tile = ntraj_tile + 1
      ELSE
 !       write(*,*) 'TRAJ: n,x,y,z = ',n, trjdat(lx,n,it), trjdat(ly,n,it), trjdat(lz,n,it)
       CYCLE
      ENDIF
 ! interpolate to parcel location

! reset trajectory values
!          trjdat(1:lx-1,n,it) = 0 ! reset other values
!          trjdat(lz+1:lxn-1,n,it) = 0 ! reset other values
!          IF ( ninfo > lzn ) trjdat(lzn+1:ninfo,n,it) = 0 ! reset other values
 
!      write(*,*) 'TRAJ: time, parcel, my_rank = ',time_real,n,my_rank


!
!  find lower left corner of scalar grid cell that trajectory
!  end point is currently in for interpolation of scalars
!
      x = trjdat(lx,n,it)
      y = trjdat(ly,n,it)
      z = trjdat(lz,n,it)

      ic = ifix((x+eps)/dx + 0.5) - ixbeg + 1
      jc = ifix((y+eps)/dy + 0.5) - jybeg + 1
      
      IF ( (myproci == 1 .and. ic == 0) .or. (myprocj == 1 .and. jc == 0 )) THEN
       ntraj_tile = ntraj_tile - 1
       CYCLE
      ENDIF
      
      kc = 1
      DO kz=1,nz-1
        IF ( gzt(kz,1) .le. z ) THEN
          kc = kz
        ENDIF
      ENDDO
!      loc(1,n,it) = ic
!      loc(2,n,it) = jc
!      loc(3,n,it) = kc
      
! Prevent accessing points with invalid data
      IF ( kc >= nz-1 ) THEN 
         kc = nz - 2
         facz = 0.99
         ntraj_tile = ntraj_tile - 1
         CYCLE
      ENDIF

      facx = (x/dx - Real(ic+ixbeg-1) + 0.5)
      facy = (y/dy - Real(jc+jybeg-1) + 0.5)
      facz = (z - gzt(kc,1))*gzt(kc,3)
      
      IF ( ic < 1 ) THEN
         ic = 1
         facx = 0.0
      ENDIF
      IF ( jc < 1 ) THEN
         jc = 1
         facy = 0.0
      ENDIF
      
      
      IF ( idebug .ge. 1 ) write(6,*) 'ic,jc,kc = ',ic,jc,kc,facx,facy,facz
      IF ( idebug .ge. 1 ) write(6,*) 'facz: ',z,gzt(kc,1),gzt(kc,3)
   

      IF ( idebug .ge. 1 ) THEN
      write(0,*) 'nmicrorates = ',nmicrorates
      write(0,*) 'microrates(nx/2,ny/2,3,1) = ',microrates(nx/2,ny/2,3,1)
      ENDIF

      call mlint2  &
      (idebug, microrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nmicrorates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lfrz,n,it) = Max(0.0,dumint)
       trjdat(lmlt,n,it) = Min(0.0,dumint)
       
       IF (  idebug .ge. 1 .and. Abs(dumint) > 0.0 ) THEN
        write(0,*) 'TRAJ: n,lfrz = ',time_real,n,dumint,trjdat(lfrz,n,it),trjdat(lmlt,n,it)
       ENDIF

      call mlint2  &
      (idebug, microrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nmicrorates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 2)
       
       trjdat(ldep,n,it) = Max(0.0,dumint)
       trjdat(lsub,n,it) = Min(0.0,dumint)


       IF ( idebug .ge. 1 .and. Abs(dumint) > 0.0 ) THEN
        write(0,*) 'TRAJ: n,ldep = ',time_real,n,dumint,trjdat(ldep,n,it),trjdat(lsub,n,it)
       ENDIF

      call mlint2  &
      (idebug, microrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nmicrorates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 3)
       
       trjdat(lcnd,n,it) = Max(0.0,dumint)
       trjdat(levp,n,it) = Min(0.0,dumint)

       IF ( idebug .ge. 1 .and. Abs(dumint) > 0.0 ) THEN
        write(0,*) 'TRAJ: n,lcnd = ',time_real,n,dumint,trjdat(lcnd,n,it),trjdat(levp,n,it)
       ENDIF

      call mlint2  &
      (idebug, microrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nmicrorates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 4)
       
       trjdat(lrevp,n,it) = dumint

      call mlint2  &
      (idebug, microrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nmicrorates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 5)
       
       trjdat(lgmlt,n,it) = dumint

      call mlint2  &
      (idebug, t8,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lbuoy,n,it) = dumint

      call mlint2  &
      (idebug, t9,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lpgrd,n,it) = dumint

      call mlint2  &
      (idebug, pn,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lpp,n,it) = dumint

      call mlint2  &
      (idebug, pb,   1, 1, 1, 1, -ng+1, nz+ng, &
       1, 1, 1, 1, nz, dumint, facx, facy, facz, 1, 1, kc, 1)
       
       trjdat(lp,n,it) = dumint

      call mlint2  &
      (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lt)

       trjdat(ltt,n,it) = dumint

       trjdat(lres,n,it) = trjdat(ltt,n,it)*(1.e-5*(trjdat(lp,n,it) + trjdat(lpp,n,it)))**rcp

      call mlint2  &
      (idebug, km,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(let,n,it) = dumint

      call mlint2  &
      (idebug, dbz,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lztot,n,it) = dumint

      qx(1,:) = 0.0
      DO il = lv,lhab
      call mlint2  &
      (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, il)
       
       qx(1,il) = dumint
       IF ( il == lv) trjdat(lqt,n,it) = dumint
       IF ( lc > 1 .and. il == lc ) trjdat(lct,n,it) = dumint
       IF ( lr > 1 .and. il == lr ) trjdat(lrt,n,it) = dumint
       IF ( li > 1 .and. il == li ) trjdat(lit,n,it) = dumint
       IF ( ls > 1 .and. il == ls ) trjdat(lst,n,it) = dumint
       IF ( lh > 1 .and. il == lh ) trjdat(lgt,n,it) = dumint
       IF ( lhl > 1 .and. il == lhl ) trjdat(lht,n,it) = dumint
!       IF ( lh > 1 ) qx(1,lh) = dumint
!       IF ( lhl > 1 ) qx(1,lhl) = dumint     
      ENDDO

      IF ( lnox > 1 ) THEN
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lnox)
        trjdat(lnoxt,n,it) = dumint
      ENDIF

!      IF ( lco > 1 ) THEN
!        call mlint2  &
!         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
!          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lco)
!        trjdat(lcot,n,it) = dumint
!      ENDIF


      cx(1,:) = 0.0
      cxtmp(1,:) = 0.0
      zx(1,:) = 0.0

! get air density
      call mlint2  &
      (idebug, dn,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       rho0(1) = dumint
       trjdat(ldn,n,it) = dumint
          IF ( .not. ( rho0(1) > 0.00001 .and. rho0(1) < 2. )) THEN
            write(0,*) 'TRAJ: problem with rho0: ',rho0(1),ic,jc,kc,my_rank
            write(0,*) 'dn at corners: ',dn(ic,jc,kc),dn(ic+1,jc,kc),dn(ic,jc+1,kc),dn(ic,jc,kc+1), &
     &      dn(ic+1,jc+1,kc),dn(ic+1,jc,kc+1),dn(ic,jc+1,kc+1),dn(ic+1,jc+1,kc+1)
            call commasmpi_abort()
          ENDIF
       rhovt(1) = Sqrt(1.225/rho0(1))
       temcg(1) = trjdat(ltt,n,it) - 273.15


! fill temp arrays to calculate fall speeds for ZVD scheme, if needed
      IF ( it >= 1 .and. microp(1:1) == 'Z' ) THEN
!       qx(1,lv) = trjdat(lqt,n,it)
!       qx(1,lc) = trjdat(lct,n,it)
!       qx(1,lr) = trjdat(lrt,n,it)
!       qx(1,li) = trjdat(lit,n,it)
!       qx(1,ls) = trjdat(lst,n,it)
!       qx(1,lh) = trjdat(lgt,n,it)
!       IF ( lhl > 1 ) qx(1,lhl) = trjdat(lqt,n,it)
       
       cwnccn(1) = ccn
       kgs(1) = 1
       igs(1) = 1
       fadvisc(1) = 1.832e-05

! get number conc.
      DO il = lc,lhab
        IF ( ln(il) > 1 ) THEN
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, ln(il))
       
        ELSE
          dumint = 0.0
        ENDIF
      !  write(0,*) 'TRAJ: il,ln,dumint,an = ',il,ln(il),dumint,an(ic,jc,kc,ln(il))
        
        cx(1,il) = Max(0.0, dumint)
        cxtmp(1,il) = cx(1,il)
      IF ( cx(1,il) < -1 ) THEN
        write(0,*) 'TRAJ: What the heck? il,cx = ',il,cx
        write(0,*) 'i,j,k,myrank = ',ic,jc,kc,my_rank
        write(0,*) 'q,c = ',qx(1,lh),cx(1,lh)
      ENDIF

       
      ENDDO

! get liquid on ice
      qxw(1,:) = 0.0
      IF ( ldoliq ) THEN
      DO il = ls,lhab
        IF ( lliq(il) > 1 ) THEN
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lliq(il))
       
        ELSE
          dumint = 0.0
        ENDIF
       
        qxw(1,il) = dumint
       
      ENDDO
      ENDIF

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
!        IF ( ls .gt. 1 .and. lvs .eq. 0 ) xdn(mgs,ls) = xdn0(ls)
!        IF ( lh .gt. 1 .and. lvh .eq. 0 ) xdn(mgs,lh) = xdn0(lh)
        IF ( li .gt. 1 )  xdn(mgs,li) = xdn0(li)
        IF ( ls .gt. 1 )  xdn(mgs,ls) = xdn0(ls)
        IF ( lh .gt. 1 )  xdn(mgs,lh) = xdn0(lh)
        IF ( lhl .gt. 1 ) xdn(mgs,lhl) = xdn0(lhl)
      end do
       
! get volumes
      vx(1,:) = 0.0
      IF ( ldovol ) THEN
      DO il = li,lhab
        IF ( lvol(il) > 1 ) THEN
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lvol(il))
       
        ELSE
          dumint = 0.0
        ENDIF
       
        vx(1,il) = dumint
       
      ENDDO
      ENDIF

      DO il = lg,lhab
        alpha(1,il) = dnu(il)
      ENDDO
      
      alpha(1,lr) = xnu(lr)

! get 6th moments
      zx(1,:) = 0.0
      IF ( ipconc >= 6 )  THEN

      IF ( lzx(lr) > 1 .and. imurain == 3 ) THEN
        write(0,*) 'TRAJ: STOP: must have imurain = 1'
        CALL commasmpi_abort()
      ENDIF

      DO il = lr,lhab
       IF ( lzx(il) > 1 ) THEN
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, lzx(il))
       
       
        zx(1,il) = dumint
        mgs = 1
        
        IF ( zx(mgs,il) .gt. zxmin .and. cx(mgs,il) > cxmin .and. qx(mgs,il) > qxmin(il) ) THEN

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z1   = zx(mgs,il)

           rdi = z1*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z1,xv
           DO i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,il))/(Abs(alpha(mgs,il))+1.e-5) .lt. 1. ) EXIT
             IF ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) EXIT
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO

        ENDIF
       ENDIF
        
      ENDDO
      
      IF ( lzr > 1 ) THEN
        trjdat(lzrt,n,it) = Max(-100., 10*Log10(Max(1.e-10, 1.0e18*zx(1,lr) ) ) ) 
      ENDIF

      IF ( lzh > 1 ) THEN
        trjdat(lzht,n,it) = Max(-100., 10*Log10(Max(1.e-10, 1.0e18*0.224*zx(1,lh) ) ) ) 
      ENDIF

      IF ( lzhl > 1 .and. lhl > 1) THEN
        trjdat(lzhlt,n,it) = Max(-100., 10*Log10(Max(1.e-10, 1.0e18*0.224*zx(1,lhl) ) ) ) 
      ENDIF
      
      ENDIF ! ipconc >= 6
      


!      write(0,*) 'traj: rho0,rhovt: ',rho0(1),rhovt(1)

      xdia(:,:,:) = 0.0
      vtxbar(:,:,:) = 0.0
      alphan(:,:) = alpha(:,:)
!       write(0,*) 'call setvt: rho = ',rho0(1)
!       DO il = lc,lhab
!        write(0,*) '1: il,qx,cx,xdia,vt = ',il,qx(1,il),cx(1,il),xdia(1,il,3),vtxbar(1,il,1),xvmn(il),xvmx(il),qxmin(il)
!       ENDDO
      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,  &
     &                 ipconc,0,ngs,1,igs,kgs,cwnccn,fadvisc, &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,    &
     &                 itype1,itype2,temcg,1,alpha,alphan,axx,bxx,0)

!       DO il = lc,lhab
!        write(0,*) '2: il,qx,cx,xdia,vt = ',il,qx(1,il),cx(1,il),xdia(1,il,3),vtxbar(1,il,1)
!       ENDDO
      ENDIF ! microp = Z

      
    ! mass-weighted fall speeds for Takahashi bin
      IF ( microp(1:3) == 'TAK' ) THEN
        xdia(:,:,:) = 0.0
       ! rain vt
          totn = 0.0
          totvt = 0.0
          totq = 0.0
          totz = 0.0
        
        DO l= 1,ntakrd
         ia = lnr + l - 1
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, ia)
          
          crbin(l) = dumint
          
          IF ( l >= 15 .and. dumint > takcxmin ) THEN
          totn = totn + dumint
          totvt = totvt + dumint*rma(l,3)*vw(l)  ! Number * mass * fallspeed
          totq = totq + dumint*rma(l,3)   ! total mass
          totz = totz + szr(l)*dumint
          ENDIF
          
          
        ENDDO
         qx(1,lr) = 1.e6*totq/rho0(1)
         cx(1,lr) = 1.e6*totn
         zx(1,lr) = 1.d6*totz
         trjdat(lzrt,n,it) = Max(-100., 10*Log10(Max(1.e-20, zx(1,lr) ) ) ) 
         trjdat(lrt,n,it) = qx(1,lr)
         IF ( totq > 1.e-20 ) THEN
          vtxbar(1,lr,1) = 0.01*totvt/totq*rhovt(1)
         ELSE
          vtxbar(1,lr,1) = 0.0
         ENDIF

  ! cloud droplets
          totn = 0.0
          totvt = 0.0
          totq = 0.0
          totz = 0.0

        DO l= 1,14
          IF ( crbin(l) > takcxmin ) THEN
          dumint = crbin(l)
          totn = totn + dumint
!          totvt = totvt + dumint*rma(l,3)*vw(l)  ! Number * mass * fallspeed
          totq = totq + dumint*rma(l,3)   ! total mass
          ENDIF
        ENDDO
         qx(1,lc) = 1.e6*totq/rho0(1)
         cx(1,lc) = 1.e6*totn
         trjdat(lct,n,it) = qx(1,lc)

   ! graupel vt
          kf = 1
          totn = 0.0
          totvt = 0.0
          totq = 0.0
          totz = 0.0
          
        DO l= 1,ntakpd
         ia = lnh + l - 1
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, ia)
          
          chbin(l,kf) = dumint
          
          IF ( l >= lf150 .and. dumint > takcxmin ) THEN
          totn = totn + dumint
          totvt = totvt + dumint*rma(l,kf)*vf(l,kf)  ! Number * mass * fallspeed
          totq = totq + dumint*rma(l,kf)   ! total mass
          totz = totz + szf(l,kf)*dumint
          ENDIF
          
          
        ENDDO
         qx(1,lh) = 1.e6*totq/rho0(1)
         cx(1,lh) = 1.e6*totn
         zx(1,lh) = 0.224*1.d6*totz
         trjdat(lzht,n,it) = Max(-100., 10*Log10(Max(1.e-20, zx(1,lh) ) ) ) 
         trjdat(lgt,n,it) = qx(1,lh)
         IF ( totq > 1.e-20 ) THEN
          vtxbar(1,lh,1) = 0.01*totvt/totq*rhovt(1)
         ELSE
          vtxbar(1,lh,1) = 0.0
         ENDIF

  ! frozen drops/hail vt
          kf = 2
          totn = 0.0
          totvt = 0.0
          totq = 0.0
          totz = 0.0
          
        DO l= 1,ntakpd
         ia = lnhl + l - 1
        call mlint2  &
         (idebug, an,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, na, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, ia)
          
          chbin(l,kf) = dumint
          
          IF ( l >= lf150 .and. dumint > takcxmin ) THEN
          totn = totn + dumint
          totvt = totvt + dumint*rma(l,kf)*vf(l,kf)  ! Number * mass * fallspeed
          totq = totq + dumint*rma(l,kf)   ! total mass
          totz = totz + szf(l,kf)*dumint
          IF ( .not. ( totq > -1 .and. totq < 100 )) THEN
            write(0,*) 'TRAJ: problem with TAK hl totq: ',totq,dumint,rma(l,kf),l,kf
            call commasmpi_abort()
          ENDIF
          ENDIF
          
          
        ENDDO
         qx(1,lhl) = 1.e6*totq/rho0(1)
         cx(1,lhl) = 1.e6*totn
          IF ( .not. ( qx(1,lhl) > -1. .and. qx(1,lhl) < 100. )) THEN
            write(0,*) 'TRAJ: problem with TAK hl qx(1,lhl): ',qx(1,lhl),totq,rho0(1)
            write(0,*) 'ic,jc,kc,my_rank = ',ic,jc,kc,my_rank
            write(0,*) 'dn at corners: ',dn(ic,jc,kc),dn(ic+1,jc,kc),dn(ic,jc+1,kc),dn(ic,jc,kc+1), &
     &      dn(ic+1,jc+1,kc),dn(ic+1,jc,kc+1),dn(ic,jc+1,kc+1),dn(ic+1,jc+1,kc+1)
            call commasmpi_abort()
          ENDIF
         trjdat(lht,n,it) = qx(1,lhl)
         zx(1,lhl) = 0.224*1.d6*totz
         trjdat(lzhlt,n,it) = Max(-100., 10*Log10(Max(1.e-20, zx(1,lhl) ) ) ) 
!         write(0,*) 'TRAJ: Tak hl qx,cx,totq = ',qx(1,lhl),cx(1,lhl),totq,rho0(1)
         IF ( totq > 1.e-20 ) THEN
          vtxbar(1,lhl,1) = 0.01*totvt/totq*rhovt(1)
         ELSE
          vtxbar(1,lhl,1) = 0.0
         ENDIF

!           trjdat(lit,n,it) = 0.0
           trjdat(lst,n,it) = 0.0

  ! put small graupel and frozen drops into snow
          totn = 0.0
          totvt = 0.0
          totq = 0.0
          totz = 0.0

        DO kf = 1,2
        DO l= 1,lf150-1
          IF ( chbin(l,kf) > takcxmin ) THEN
          dumint = chbin(l,kf)
          totn = totn + dumint
          totvt = totvt + dumint*rma(l,kf)*vf(l,kf)  ! Number * mass * fallspeed
          totq = totq + dumint*rma(l,kf)   ! total mass
          ENDIF
        ENDDO
        ENDDO
         qx(1,ls) = 1.e6*totq/rho0(1)
         cx(1,ls) = 1.e6*totn
         trjdat(lst,n,it) = qx(1,ls)

      
      ENDIF

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

      call mlint2  &
      (idebug, u, -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, uint, facxu, facyu, faczu, iu,ju,ku, 1)

      
      trjdat(lut,n,it) = uint
      
      IF ( idebug .ge. 1 ) write(6,*) 'iu,ju,ku = ',iu,ju,ku,facxu,facyu,faczu,uint

! Electricity
      IF ( ipelec > 0 .and. allocated( elecrates ) ) THEN

      ! Ex
      call mlint2  &
      (idebug, elec(iex)%flt3d,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facxu, facyu, faczu,iu,ju,ku, 1)
       
       trjdat(lex,n,it) = dumint
       exint = dumint
      

      ! Ex-pre lightning
      IF ( icomtraj == 0 ) THEN
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxu, facyu, faczu,iu,ju,ku, 1)
       
       trjdat(lexpl,n,it) = dumint
      ELSE
       trjdat(lexpl,n,it) = trjdat(lex,n,it)
      ENDIF
      
      
      ! delta-Ex
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxu, facyu, faczu,iu,ju,ku, 5)
       
       trjdat(ldelex,n,it) = dumint
      
      ENDIF
   
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

      call mlint2  &
      (idebug, v,  -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, vint, facxv, facyv, faczv, iv,jv,kv, 1)
      
      trjdat(lvt,n,it) = vint

      IF ( idebug .ge. 1 ) write(6,*) 'iv,jv,kv = ',iv,jv,kv,facxv,facyv,faczv,vint

! Electricity
      IF ( ipelec > 0 .and. allocated( elecrates ) ) THEN

      ! Ey
      call mlint2  &
      (idebug, elec(iey)%flt3d,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facxv, facyv, faczv,iv,jv,kv, 1)
       
       trjdat(ley,n,it) = dumint
       eyint = dumint
      
      ! Ey-prelightning
      IF ( icomtraj == 0 ) THEN
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxv, facyv, faczv,iv,jv,kv, 2)
       
       trjdat(leypl,n,it) = dumint
      ELSE
       trjdat(leypl,n,it) = trjdat(ley,n,it)
      ENDIF


      ! delta-Ey
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxv, facyv, faczv,iv,jv,kv, 6)
       
       trjdat(ldeley,n,it) = dumint
      
      ENDIF
      
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

! Electricity
      IF ( ipelec > 0 .and. allocated( elecrates ) ) THEN

      ! Ez
      call mlint2  &
      (idebug, elec(iez)%flt3d,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facxw, facyw, faczw,iw,jw,kw, 1)
       
       trjdat(lez,n,it) = dumint
       ezint = dumint
       trjdat(lemag,n,it) = ( exint**2 + eyint**2 + ezint**2 )**0.50
      

      ! Ez pre-lightning
      IF ( icomtraj == 0 ) THEN
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxw, facyw, faczw,iw,jw,kw, 3)
       
       trjdat(lezpl,n,it) = dumint
      ELSE
       trjdat(lezpl,n,it) = trjdat(lez,n,it)
      ENDIF

       trjdat(lemagpl,n,it) = Sqrt(dumint**2 + trjdat(leypl,n,it)**2 + trjdat(lexpl,n,it)**2)

      ! delta-Ez
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facxw, facyw, faczw,iw,jw,kw, 7)
       
       trjdat(ldelez,n,it) = dumint
      
      ENDIF
      
      call mlint2  &
      (idebug, w, -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, wint, facxw, facyw, faczw, iw,jw,kw, 1)
      
      vt = 0.0
      
      IF ( it == 2 ) THEN
      IF ( microp(1:1) == 'K' ) THEN ! rain fall speed for Kessler (ntrajtype 2)
!       rcgs(k)   = 1.0e2*pb(k)**2.509/(287.04*tb(k))
!       vtden(k)  = sqrt(0.0011225/rcgs(k))
!        vt(i,j,k)     = 36.34*(Max( 0.0,qr(i,j,k) )*rcgs(k))**0.1364 * vtden(k)
      call mlint2  &
      (idebug, db,   1, 1, 1, 1, -ng+1, nz+ng, &
       1, 1, 1, 1, nz, dumint, facx, facy, facz, 1, 1, kc, lt)

        rcgs =  1.e-3*dumint
        vtden = sqrt(0.0011225/rcgs)
        qr = Max( 0.0, trjdat(lrt,n,it) )
        vt = 36.34*(qr*rcgs)**0.1364 * vtden
        IF ( .not. ( vt >= 0. .and. vt < 100. ) ) THEN
          write(0,*) 'traj: problem with vt! vt = ',vt,qr,vtden,rcgs
          call commasmpi_abort()
        ENDIF
!        wint = wint - vt
      ELSEIF ( microp(1:1) == 'Z' ) THEN ! rain fall speed for ZVD (ntrajtype 2)
      
       IF ( qx(1,lr) > qxmin(lr) ) vt = vtxbar(1,lr,1)
       
!       IF ( qx(1,lr) > 0.1e-3 ) write(0,*) 'rank, qr, vtr = ',my_rank,1.e3*qx(1,lr),vt

      ELSEIF ( microp(1:3) == 'TAK' ) THEN ! rain fall speed for TAK bin micro
        vt = vtxbar(1,lr,1)
      ENDIF
      ENDIF ! it == 2

      IF ( it == 3 ) THEN
      IF ( microp(1:1) == 'Z' ) THEN ! graupel + hail
      
        IF ( lhl > 1 ) THEN
         IF ( qx(1,lh) + qx(1,lhl) > Min(qxmin(lh),qxmin(lhl)) ) THEN
           vt = (qx(1,lh)*vtxbar(1,lh,1) + qx(1,lhl)*vtxbar(1,lhl,1) ) / ( qx(1,lh) + qx(1,lhl) )
         ENDIF
        ELSE
         vt = vtxbar(1,lh,1)
        ENDIF

        numer = qx(1,lh)*vtxbar(1,lh,1)
        denom = qx(1,lh)
        qmin = qxmin(lh)
        
        IF ( lhl > 1 ) THEN
          numer = numer + qx(1,lhl)*vtxbar(1,lhl,1)
          denom = denom + qx(1,lhl)
          qmin = Min( qmin, qxmin(lhl) )
        ENDIF

        IF ( denom > qmin ) vt = numer/denom

      ELSEIF ( microp(1:3) == 'TAK' ) THEN ! rain fall speed for TAK bin micro
         IF ( qx(1,lh) + qx(1,lhl) > Min(qxmin(lh),qxmin(lhl)) ) THEN
           vt = (qx(1,lh)*vtxbar(1,lh,1) + qx(1,lhl)*vtxbar(1,lhl,1) ) / ( qx(1,lh) + qx(1,lhl) )
         ENDIF

      ENDIF
      ENDIF ! it == 3

      IF ( it == 4 ) THEN
      IF ( microp(1:1) == 'Z' ) THEN ! rain + graupel + hail
      
        numer = qx(1,lr)*vtxbar(1,lr,1) + qx(1,lh)*vtxbar(1,lh,1)
        denom = qx(1,lr) + qx(1,lh)
        qmin = Min( qxmin(lr), qxmin(lh) )
        
        IF ( lhl > 1 ) THEN
          numer = numer + qx(1,lhl)*vtxbar(1,lhl,1)
          denom = denom + qx(1,lhl)
          qmin = Min( qmin, qxmin(lhl) )
        ENDIF
        
        IF ( denom > qmin ) vt = numer/denom

      ELSEIF ( microp(1:3) == 'TAK' ) THEN ! rain fall speed for TAK bin micro
         IF ( qx(1,lh) + qx(1,lhl) + qx(1,lr) > Min(qxmin(lh),qxmin(lhl)) ) THEN
           vt = (qx(1,lh)*vtxbar(1,lh,1) + qx(1,lhl)*vtxbar(1,lhl,1) + qx(1,lr)*vtxbar(1,lr,1) )  &
     &          / ( qx(1,lh) + qx(1,lhl) + qx(1,lr) )
         ENDIF
      ENDIF
      ENDIF ! it == 4

      IF ( it == 5 .and.  ( microp(1:1) == 'Z' .or.  microp(1:3) == 'TAK' ) ) THEN ! graupel
      
        vt = vtxbar(1,lh,1)
         
         IF ( trjdat(lres,n,it) > 273.15 .and. qx(1,lh) > qxmin(lh) .and. wint < 3.0 ) THEN ! weighted between rain and graupel for melting in weak updraft
!          IF ( qx(1,lh) + qx(1,lr) > qxmin(lh) ) THEN
           vt = (qx(1,lh)*vtxbar(1,lh,1) + qx(1,lr)*vtxbar(1,lr,1) ) / ( qx(1,lh) + qx(1,lr) )
!          ENDIF
         ENDIF

      ENDIF ! it == 5

      IF ( it == 6 ) THEN
      IF ( microp(1:1) == 'Z' .and. lhl > 1 ) THEN ! hail

        vt = vtxbar(1,lhl,1)
        a = 5.
        IF ( qx(1,lhl) < a*qxmin(lhl) .and. qx(1,lhl) > qxmin(lhl)  ) THEN
         vt = vt*(  qx(1,lhl)/((a-1.)*qxmin(lhl)) - 1./(a-1.) )
        ENDIF
      
      ELSEIF ( microp(1:3) == 'TAK' ) THEN
        vt = vtxbar(1,lhl,1)

      
      ENDIF
      ENDIF ! it == 6

      trjdat(lwt,n,it) = wint

      wint = wint - vt
      
      trjdat(lw2,n,it) = wint

      IF ( ipelec > 0 ) THEN
      
      wint = wint + riserate
      
      trjdat(lwtrise,n,it) = wint
      
      ENDIF

      trjdat(lgdia,n,it) = xdia(1,lh,3)
      trjdat(lrdia,n,it) = xdia(1,lr,3)
      IF ( lhl > 1 ) THEN
        trjdat(lhdia,n,it) = xdia(1,lhl,3)
      ELSE
        trjdat(lhdia,n,it) = 0.0
      ENDIF
!      IF ( xdia(1,lhl,3) < -1 .or. n == 1 ) THEN
!        write(0,*) 'TRAJ: n, it, hwdia = ',n,it,xdia(1,lhl,3)
!        write(0,*) 'i,j,k,myrank = ',ic,jc,kc,my_rank
!        write(0,*) 'x,y,z = ',trjdat(lx,n,it),trjdat(ly,n,it),trjdat(lz,n,it)
!        write(0,*) 'q,c = ',qx(1,lhl),cx(1,lhl)
!      ENDIF
      trjdat(lcg,n,it) = cxtmp(1,lh)
      
      IF ( ninfo >= lchl .and. lhl > 1 ) trjdat(lchl,n,it) = cxtmp(1,lhl)
      IF ( ninfo >= lccw ) trjdat(lccw,n,it) = cxtmp(1,lc)
      IF ( ninfo >= lcci ) trjdat(lcci,n,it) = cxtmp(1,li)
      IF ( ninfo >= lcsw ) trjdat(lcsw,n,it) = cxtmp(1,ls)
      
      trjdat(lvtg,n,it) = vtxbar(1,lh,1)
      IF ( lhl > 1 ) trjdat(lvth,n,it) = vtxbar(1,lhl,1)
      trjdat(lvtr,n,it) = vtxbar(1,lr,1)


      
      IF ( idebug .ge. 1 ) write(6,*) 'iw,jw,kw = ',iw,jw,kw,facxw,facyw,faczw

! Electricity
      IF ( ipelec > 0 ) THEN

      ! potential
      
      call mlint2  &
      (idebug, elec(ipot)%flt3d,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lphi,n,it) = dumint
       trjdat(lphipl,n,it) = trjdat(lphi,n,it) 

      IF (  allocated( elecrates ) ) THEN
      
      ! potential before lightning

      IF ( icomtraj == 0 ) THEN
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 4)
       
       trjdat(lphipl,n,it) = dumint
      ELSE
       trjdat(lphipl,n,it) = trjdat(lphi,n,it) 
      ENDIF

      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 8)
       
       trjdat(ldelphi,n,it) = dumint


      ! RAR
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 9)
       
       trjdat(lrarh,n,it) = dumint

      ! RGHIS
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 10)
       
       trjdat(lcrgis,n,it) = dumint

      ! EHW
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 11)
       
       trjdat(lehw,n,it) = dumint
       

      ! RAR
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 12)
       
       trjdat(lrarhl,n,it) = dumint

      ! RGHIS
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 13)
       
       trjdat(lcrhis,n,it) = dumint

      ! EHW
      call mlint2  &
      (idebug, elecrates,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, nelecrates, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 14)
       
       trjdat(lehlw,n,it) = dumint
       
       
      ENDIF
      
      ! net charge
      
      call mlint2  &
      (idebug, elec(iscnet)%flt3d,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
       1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)
       
       trjdat(lnetchg,n,it) = 1.e9*dumint
      
       dumint = 0.0
       kz = 0
       k1 = Max(1,kc-chgavez)
       k2 = Min(nz-1,kc+chgavez)
       x1 = Max( -ng, ic-chgavex )
       x2 = Min( nx+ng, ic+chgavex )
       y1 = Max( -ng, jc-chgavex )
       y2 = Min( ny+ng, jc+chgavex )
       
       IF ( chgavez > 0 ) THEN ! average vertical layers with equal weight
         DO k = k1,k2
         DO ix = x1,x2
         DO jy = y1,y2
           IF ( (jy-jc)**2 + (ix-ic)**2 <= chgavex**2 ) THEN ! approximate a circle
            kz = kz + 1
            dumint = dumint + elec(iscnet)%flt3d(ix,jy,k)
           ENDIF
         ENDDO
         ENDDO
         ENDDO
        trjdat(lnetchgave,n,it) = 1.e9*dumint/float(kz)

       ELSE ! for chgavez=0, interpolate the charge layer for smoother data

!  The t0 array has the local disk average for each grid cell, so just interpolate in 3D
!   with the subroutine
         call mlint2  &
         (idebug, t0,   -ng+1, nx+ng, -ng+1, ny+ng, -ng+1, nz+ng, &
          1, 1, nx, ny, nz, dumint, facx, facy, facz,ic,jc,kc, 1)

          trjdat(lnetchgave,n,it) = dumint

        ENDIF

!        trjdat(lnetchgave,n,it) = 1.e9*dumint/float(kz)

        trjdat( lnumflash,n,it) = numflashtraj

      
      ENDIF


      IF ( icomtraj < 2 .and. iverttraj == 0 ) THEN ! for icomtraj == 2 or iverttraj > 0, keep the sounding vertical
        trjdat(lxn,n,it) = Max(0.5*dx, Min( xdomain - 0.5*dx, x + dt*uint ) )
        trjdat(lyn,n,it) = Max(0.5*dy, Min( xdomain - 0.5*dy, y + dt*vint ) )
      ELSE
        trjdat(lxn,n,it) = trjdat(lx,n,it) ! x
        trjdat(lyn,n,it) = trjdat(ly,n,it) ! y
      ENDIF
      trjdat(lzn,n,it) = Max(0.0, Min(gzt(nz-1,1),  z + dt*wint) )
      
      if ( trjdat(lp,n,it) .lt. 999999.  ) then
      pres = trjdat(lp,n,it) + trjdat(lpp,n,it)
      tc = trjdat(ltt,n,it)
      qvc = trjdat(lqt,n,it)
      trjdat(lthe,n,it) = thetae(qvc,tc,pres)
      IF ( idebug .ge. 1 ) print*,qvc,tc,pres,trjdat(lthe,n,it)
      end if
!  
      trjdat(load,n,it) =   &
        -g*(trjdat(lct,n,it) &
          +trjdat(lrt,n,it) &
          +trjdat(lit,n,it) &
          +trjdat(lst,n,it) &
          +trjdat(lht,n,it) &
          +trjdat(lgt,n,it))
!      trjdat(loadr,n,it) = g*(trjdat(lrt,n,it))
!      trjdat(loads,n,it) = g*(trjdat(lst,n,it))
!      trjdat(loadg,n,it) = g*(trjdat(lgt,n,it))
!      trjdat(loadh,n,it) = g*(trjdat(lht,n,it))
!      trjdat(loadc,n,it) = g*(trjdat(lct,n,it)+trjdat(lit,n,it))

#ifdef MPI
! load data to send to rank 0 
      trjbuf(n0:ninfo,ntraj_tile) = trjdat(n0:ninfo,n,it)
      trjbuf(ninfo+1,ntraj_tile) = n
      trjbuf(ninfo+2,ntraj_tile) = it
#endif
   
     ENDDO ! n
     ENDDO ! it

! DO MPI COMMS HERE -- borrowed from CM1 (thanks, George!)
#ifdef MPI
      IF ( number_of_processes > 1 ) THEN ! will get stuck if only one process
      IF ( my_rank == 0 ) THEN

        DO proc = 1,number_of_processes-1
          CALL MPI_RECV(ntraj_tile,1,MPI_INTEGER,proc,proc,my_comm,mpi_status,mpi_error_code)
          IF ( ntraj_tile > 0 ) CALL MPI_RECV(trjbuf,(ninfo+2)*ntraj_tile,MPI_REAL,proc,1000+proc,my_comm,mpi_status,mpi_error_code)
          DO n = 1,ntraj_tile
            it = NInt((trjbuf(ninfo+2,n)))
            in = Nint(trjbuf(ninfo+1,n))
            trjdat(n0:ninfo,in,it) = trjbuf(n0:ninfo,n)
          ENDDO
        ENDDO

      ELSE

        CALL MPI_SEND(ntraj_tile,1,MPI_INTEGER,0,my_rank,my_comm,mpi_error_code)
        IF ( ntraj_tile > 0 ) CALL MPI_SEND(trjbuf,(ninfo+2)*ntraj_tile,MPI_REAL,0,1000+my_rank,my_comm,mpi_error_code)

      ENDIF

      IF ( my_rank == 0 ) THEN
       DO it = 1,ntrajtype
        DO n = 1,ntraj
          in = n + ntraj*(it-1)
          trjloc(1,in) = trjdat(lxn,n,it)
          trjloc(2,in) = trjdat(lyn,n,it)
          trjloc(3,in) = trjdat(lzn,n,it)
        ENDDO
       ENDDO
      ENDIF

      CALL MPI_BARRIER (my_comm,mpi_error_code)
      CALL MPI_BCAST(trjloc,3*ntraj*ntrajtype,MPI_REAL,0,my_comm,mpi_error_code)

      IF ( my_rank /= 0 ) THEN
        DO it = 1,ntrajtype
        DO n = 1,ntraj
          in = n + ntraj*(it-1)
          trjdat(lxn,n,it) = trjloc(1,in)
          trjdat(lyn,n,it) = trjloc(2,in)
          trjdat(lzn,n,it) = trjloc(3,in)
        ENDDO
        ENDDO
      ENDIF
      
    ENDIF ! number_of_processes
#endif

   IF ( my_rank == 0 ) THEN
      DO it = 1,ntrajtype
      
      IF ( itraj < 3 ) THEN ! write out netcdf
      
      status = NF90_OPEN(ncfiletraj(it) ,NF90_WRITE,ncid)

      start1d(1) = trajstep
      count1d(1) = 1
      data1d(1)  = time_real
      status = nf90_put_var(ncid, timeid, data1d, start1d, count1d )
      start2d(1) = trajstep
!      write(0,*) 'TRAJ: trajstep,time_real = ',trajstep,time_real

      DO nt = n0,ninfo
!      DO nt = 1,1
 !      IF ( .not. ( nt == lxn .or. nt == lyn .or. nt == lzn ) ) THEN

! changed this to invert the loops (outer loop on ninfo, inner on ntraj) and write all the current values at once
! could try it with the dimensions switched back and use a stride value so that nf90_put_var is still called only 
! once per variable (ninfo times) instead of the original inefficient ninfo*ntraj calls, but writing contiguous
! data is probably a good bit faster. It is also easy enough to transpose the arrays after reading them to get into
! trajectory space
      DO n = 1,ntraj
        data2d(n,1) = trjdat(nt,n,it)
      ENDDO
      
        
      start2d(1) = 1
      start2d(2) = trajstep
      count2d(1) = ntraj
      count2d(2) = 1

        
        status = nf90_put_var(ncid, trajvars(nt)%varid, data2d, start2d, count2d )
        IF(status /= NF90_NOERR) THEN 
          write(0,*) 'NF90_PUT_VAR:  Error writing variable: ', trajvars(nt)%name
          write(0,*) 'start2d = ',start2d(1),start2d(2)
          write(0,*) 'count2d = ',count2d(1),count2d(2)
          write(0,*) 'start2d = ',start2d(1),start2d(2)
          call commasmpi_abort()
        ENDIF
      
!       ENDIF
      ENDDO

      status = nf90_close(ncid)
      
      ENDIF ! itraj < 3

      IF ( itraj >= 2 ) THEN ! write text files
      DO n = 1,ntraj

      ic = ifix(( trjdat(lx,n,it) +eps)/dx + 0.5)
      jc = ifix(( trjdat(ly,n,it)+eps)/dy + 0.5)
      kc = 1
      DO kz=1,nz-1
        IF ( gzt(kz,1) .le. trjdat(lz,n,it) ) THEN
          kc = kz
        ENDIF
      ENDDO
      
      write(trj_unit+it-1,*) '==================================================='
      write(trj_unit+it-1,*) 'TRAJ'
      write(trj_unit+it-1,71) n, 1, time_real, ic, jc, kc
      write(trj_unit+it-1,72) trjdat(lx,n,it),trjdat(ly,n,it),trjdat(lz,n,it)
      write(trj_unit+it-1,73) trjdat(lut,n,it),trjdat(lvt,n,it),trjdat(lwt,n,it),trjdat(lw2,n,it)
      write(trj_unit+it-1,72) trjdat(lp,n,it),trjdat(lpp,n,it),trjdat(ltt,n,it)
      write(trj_unit+it-1,72) trjdat(let,n,it),trjdat(lqt,n,it),trjdat(lct,n,it)
      write(trj_unit+it-1,72) trjdat(lrt,n,it),trjdat(lit,n,it),trjdat(lst,n,it)
      write(trj_unit+it-1,72) trjdat(lgt,n,it),trjdat(lfrz,n,it),trjdat(lmlt,n,it)
      write(trj_unit+it-1,72) trjdat(ldep,n,it),trjdat(lsub,n,it),trjdat(lcnd,n,it)
      write(trj_unit+it-1,72) trjdat(levp,n,it),trjdat(lrevp,n,it),trjdat(lgmlt,n,it)
!      write(trj_unit+it-1,72) trjdat(loadr,n,it),trjdat(loads,n,it),trjdat(loadg,n,it)
      write(trj_unit+it-1,72) trjdat(lht,n,it),trjdat(load,n,it),trjdat(lpgrd,n,it)
      write(trj_unit+it-1,72) trjdat(lbuoy,n,it),trjdat(lthe,n,it),trjdat(lres,n,it)
!      write(trj_unit+it-1,73) trjdat(lht,n,it),trjdat(loadh,n,it)
!      write(trj_unit+it-1,72)  trjdat(lgdia,n,it),trjdat(lcg,n,it),trjdat(lvtg,n,it) 

      
      
      
      
      

#ifndef MPI
     IF ( idebug >= 1 ) THEN
      IF ( lzh > 1 ) THEN
      write(trj_unit+it-1,*) 'qh,ch,zh = ',an(ic,jc,kc,lh),an(ic,jc,kc,lnh),an(ic,jc,kc,lzh)
      write(trj_unit+it-1,*) 'qh,kp1   = ',an(ic,jc,kc+1,lh),an(ic,jc,kc+1,lnh),an(ic,jc,kc+1,lzh)
      ELSE
      write(trj_unit+it-1,*) 'qh,ch = ',an(ic,jc,kc,lh),an(ic,jc,kc,lnh)
      ENDIF
     ENDIF
#endif
     ENDDO ! n
      ENDIF ! itraj >= 2
     ENDDO ! it
   ENDIF 
71    format(1x,2i6,1x,f10.3, 3(1x,i4))
72    format(1x,3(1x,1pe13.5))
73    format(1x,4(1x,1pe13.5))
74    format(1x,3(1x,f13.5))

   
   deallocate( dn )
   deallocate( pn )
   IF ( allocated( chgnet ) ) deallocate( chgnet )

   microrates(:,:,:,:) = 0.0
   IF ( allocated( elecrates ) ) THEN
     elecrates(:,:,:,9:14) = 0.0
   ENDIF
   
   RETURN
   
   END SUBROUTINE TRAJ

   
END MODULE TRAJ_MODULE
   