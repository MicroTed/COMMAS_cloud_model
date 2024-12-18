!-----------------------------------------------------------------------------
!
! COMMASPECTRA: A simple hack to print out size spectra from the Takahashi microphysics
!
!
!-----------------------------------------------------------------------------

 PROGRAM COMSPECTRA

#ifdef NAG
  USE F90_UNIX_PROC
#endif

  USE GRID_MODULE
  USE GRIDIO_MODULE
  USE FILE_MODULE
  USE CLINE_MODULE
  USE CPUTIME_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
!  USE VIS5D_MODULE
  USE INDEX_MODULE

  implicit none

!-----------------------------------------------------------------------------
! GRID DEFINITION
 
  TYPE(GRID) :: gd, ge

!-----------------------------------------------------------------------------
! COMMAND LINE VARIABLES

  character(LEN = 120) :: run_file
  integer              :: status, length, length1,itmp
  real                 :: tmp, chartofloat

!-----------------------------------------------------------------------------
! DATE DIAGNOSTICS

  character( LEN = 80 ) :: command
  character( LEN = 8  ) :: date
  character( LEN = 10 ) :: hhmmss
  character( LEN = 3  ), dimension(12) :: months = (/'JAN','FEB',  &
                                                     'MAR','APR',  &
                                                     'MAY','JUN',  &
                                                     'JUL','AUG',  &
                                                     'SEP','OCT',  &
                                                     'NOV','DEC'/)
!-----------------------------------------------------------------------------
! INTEGER INDEX VARIABLES FOR THE RUN

  integer :: u 
  integer :: v  
  integer :: w   
  integer :: pi   
  integer :: km   
  integer :: s     
  integer :: dbz   
  integer :: vzf
  integer :: wz    
  integer :: elec
  integer :: xtra
  integer :: cion
  integer :: muz
   
  integer :: uinit  
  integer :: vinit  
  integer :: winit  
  integer :: piinit 
  integer :: kminit 
  integer :: sinit  

  integer :: gx     
  integer :: gy     
  integer :: gz     
  integer :: precip 

!-----------------------------------------------------------------------------
! SOLVER SCRATCH MEMORY

  real, allocatable    :: st(:,:,:,:)         ! These are big scratch arrays for solver
  real, allocatable    :: flsh(:,:,:,:,:)
  real, allocatable    :: ft(:,:)
  integer, allocatable :: tarray(:)

!-----------------------------------------------------------------------------
! OTHER MISC VARIABLES 

!  real    :: dx
!  real    :: dy
!  real    :: dz
  integer :: ibeg
  integer :: iend
  integer :: ns
  integer :: n
  integer :: ntime
  integer :: m
  integer :: time
  integer :: coards(6)
  integer :: trst       
  integer :: this       
  integer :: tprt
  integer :: tstt
  integer :: tv5d    
  integer :: tvis5dtmp

  logical :: file_exist = .false.
  logical :: io_flag, lstt
  logical :: hack = .false.

  integer i,j,k, it, nt, count

  character(LEN=100)  :: filenm
  character(LEN=120) :: ncfile

  character(len=6) stime
  character(LEN=4) number1,number2
  character(LEN=8) number
  character(LEN=2) dup

  character(LEN = 255):: v5dfldstmp

  integer  :: istat
  integer  :: visinterval
  integer  :: iskip = 1
  integer  :: recalcdbz = 0  ! =1 sets iusewetgraupel=1 and calls radardbz
  integer  :: iusewethailv5d = 0
  integer  :: iusewetgraupelv5d = 0


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

  logical :: mpisep = .false.

!-----------------------------------------------------------------------------
! Need 3D.RUN file in order to get run parameters - get from the command line

  IF( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN

   write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  NEED 3D.RUN FILENAME'
   write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  EXITING RUN!...'
   call exit(1)

  ELSE

   CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)

!   write(6,*) 'COMMAS:  Input 3D.RUN file is:  ', run_file(1:length)

   IF( status .ne. 0 ) THEN

     write(0,*) 'COMMAS:  COULD NOT RETRIEVE 1ST COMMAND LINE ARG:  EXITING RUN!...', status
     write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ!'
     call exit(1)

   ENDIF

  ENDIF
  
!-----------------------------------------------------------------------------
! Read in 3d.run file...

  INQUIRE(file=run_file(1:length), exist=file_exist)

  IF( .NOT. file_exist ) THEN

    write(0,*) 'COMMAS:  INPUT 3D.RUN FILE:  ', run_file(1:length), ' DOES NOT EXIST!!!'
    write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ, EXITING'
    stop

  ENDIF

! READ IN RUN MODEL NAMELIST

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'COM2V5D:  PROBLEM READING NAMELIST'

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COM2V5D:  PROBLEM READING NAMELIST'

  IF ( historyinput .eq. -1 ) historyinput = historyoutput


! READ IN SPECIAL FLAGS
  

!-----------------------------------------------------------------------------
! CHECK FOR ADDITIONAL COMMAND LINE ARGUMENTS

  IF( COMMAND_ARGUMENT_COUNT() .gt. 1 ) THEN

   write(6,"(1x,80('-'))")
   write(6,*) ''
   write(6,*) 'COMMAS:  USER HAS SUPPLIED ADDITIONAL COMMAND LINE ARGUMENTS, READING...'
   write(6,*) ''
   write(6,"(1x,80('-'))")

   DO n = 2,COMMAND_ARGUMENT_COUNT()

    CALL GET_COMMAND_ARGUMENT(n,command,length1,status)

    IF( status .eq. 0 ) THEN

! Do special check to make sure we are starting at the correct time...

     IF( n .eq. 2 ) THEN
      tmp = chartofloat(command(1:length1))
      member = Int(tmp)
     ENDIF

     IF( n .eq. 3 ) THEN
         tmp     = chartofloat(command(1:length1))
!         start   = Int(tmp)
!      IF( Int(tmp) .ne. start ) THEN
!     write(6,"(1x,80('-'))")
!     write(6,*) ''
!        write(6,*) 'COMMAS:  START TIME ON COMMAND LINE NE START TIME IN FILE, PROBLEM!'
!        write(6,"(1x,'START TIME FROM RESTART FILE: ', f10.2)") start
!        write(6,"(1x,'START TIME FROM COMMAND LINE: ', f10.2)") tmp
!        write(6,*) 'COMMAS IS EXITING RUN'
!     write(6,*) ''
!     write(6,"(1x,80('-'))")
!        call exit(1)
!      ELSE
       start = Int(tmp)
!      ENDIF
     ENDIF

     IF( n .eq. 4 ) THEN
         tmp     = chartofloat(command(1:length1))
         stop    = Int(tmp)
     ENDIF
!     IF( n .eq. 4 ) tprint   = chartofloat(command(1:length1))
!     IF( n .eq. 5 ) thistory = chartofloat(command(1:length1))
!     IF( n .eq. 6 ) trestart = chartofloat(command(1:length1))

    ELSE

     write(6,"(1x,71('-'))")
     write(6,*) ''
     write(6,*) 'COMMAS:  COULD NOT RETRIEVE',n,' COMMAND LINE ARG:  EXITING RUN!...', status
     write(6,*) ''
     write(6,"(1x,71('-'))")

     call exit(1)

    ENDIF

   ENDDO

  ENDIF
  
! READ IN GRID PARAMETERS


  v5dfldstmp = v5dflds


   write(number1, '(a,i3.3)') '.', member
   number = number1


!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   CALL STRING_LIMITS(prefix, ibeg, iend)


!-----------------------------------------------------------------------------
! ELSE, READ IN THE DATA


    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

    IF ( start .lt. 0 ) start = 0
    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)


   IF ( historyinput == 1 ) THEN
     ncfile = prefix(1:length)//trim(number)//'.nc'
   ELSEIF ( historyinput .eq. 3  ) THEN
     write(stime,   '(i6.6)') start
     ncfile = prefix(1:length)//'.'//stime//'.nc'
   ENDIF
  
    CALL GRID_INFO_NETCDF( ncfile, nt, microphys)
    
    allocate( tarray(nt) )
    
    CALL GRID_INFO_NETCDF( ncfile, nt, microphys, tarray)
    
    microp = microphys
    
    IF ( microp(1:3) .ne. 'TAK' ) THEN
      write(0,*) 'Program for TAK microphysics only. Stop.'
      STOP
    ENDIF
    
    CALL GRID_READ_NETCDF( gd, ncfile, tarray(1), nx_or_nxend=.true. )

    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)

    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMMAS:  Problem setting TIME'  
!    IF( .not. SET_VARIABLE(gd,'TRESTART',  trestart)) write(6,*) 'COMMAS:  Problem setting TRESTART'
!    IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMMAS:  Problem setting THISTORY'
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'
    IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COMMAS:  Problem setting TPRINT'
    IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COMMAS:  Problem setting TSTAT'
    IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COMMAS:  Problem setting STOP'
    IF( .not. SET_VARIABLE(gd,'DT',        dt)      ) write(6,*) 'COMMAS:  Problem setting DT'
   
    CALL STRING_LIMITS(v5dfldstmp, i, j)
    IF( .not. SET_ATTRIBUTE(gd,'V5DFIELDS',  v5dfldstmp)) write(6,*) 'COMMAS:  Problem setting V5DFIELDS'

    time = start


!--------------------------------------------------------------------------------------
! Open VIS5D file  

!  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
    
    IF ( start .gt. tarray(nt) ) THEN
      write(0,*) 'Starting time is greater than the last record time! STOP!'
      STOP
    ENDIF
    
    start = Max(start, tarray(1) ) 
    
    
!    v5dtimes(1) = Max( start, tarray(1) )

!    print*, 'setting TIME to ',Max(0,start-tvis5d)
!    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start-tvis5d))   ) write(6,*) 'COMMAS:  Problem setting TIME'  


     gx     = GET_VARIABLE_INDEX(gd, 'XC')
     gy     = GET_VARIABLE_INDEX(gd, 'YC')
     gz     = GET_VARIABLE_INDEX(gd, 'ZC')
   
    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'  
    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMMAS:  Problem setting TIME'  



!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

  CALL GET_VARIABLE(gd, 'BCX',       bcx)
  CALL GET_VARIABLE(gd, 'BCY',       bcy)

  CALL GET_VARIABLE(gd, 'NXEND',    nx)
  CALL GET_VARIABLE(gd, 'NYEND',    ny)
  CALL GET_VARIABLE(gd, 'NZEND',    nz)
!  CALL GET_VARIABLE(gd, 'NX',       nx)
!  CALL GET_VARIABLE(gd, 'NY',       ny)
!  CALL GET_VARIABLE(gd, 'NZ',       nz)
  CALL GET_VARIABLE(gd, 'DX',       dx)
  CALL GET_VARIABLE(gd, 'DY',       dy)
  CALL GET_VARIABLE(gd, 'DZ',       dz)
  CALL GET_VARIABLE(gd, 'DT',       dt)
  CALL GET_VARIABLE(gd, 'NSMALL',   nsmall)
  CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
  CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)
  CALL GET_VARIABLE(gd, 'IPELEC',   ipelec)
  i = ipconc
  CALL GET_VARIABLE(gd, 'IPCONC',   ipconc)
  CALL GET_VARIABLE(gd, 'ISFCPHYS',       isfcphys)
   write(0,*) 'COM2V5D: isfcphys, IPCONC = ',isfcphys,ipconc,i

  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  
  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COMMAS:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COMMAS:  Problem setting NG'

  write(*,*) 'ng = ',ng

! Allocate SOLVER scratch space
  
  m = (nx+2*ng)*(ny+2*ng)*(nz+2*ng)
!  allocate ( st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns) )
  IF ( hack ) THEN 
   allocate ( flsh(nx,ny,nz,3,2) )
   flsh = 0.0
  ENDIF
  
  allocate ( ft(m,4) ) ! temporary arrays for solver

  tprt = start + tprint   - 1
  tstt = start + tstat    - 1
  tv5d = start 

  this = thistory * (1 + start/thistory) - 1
  tv5d = start + tvis5d                  - 1

!-----------------------------------------------------------------------------
! PRINT OUT PARAMS


!-----------------------------------------------------------------------------
! Create run lock file to help facilitate parallel runs


!-----------------------------------------------------------------------------
! GET the variables needed for the integration 

 u      = GET_VARIABLE_INDEX(gd, 'U')
 v      = GET_VARIABLE_INDEX(gd, 'V')
 w      = GET_VARIABLE_INDEX(gd, 'W')
 pi     = GET_VARIABLE_INDEX(gd, 'PI')
 km     = GET_VARIABLE_INDEX(gd, 'KM')
 s      = GET_VARIABLE_INDEX(gd, 'TH')
 dbz    = GET_VARIABLE_INDEX(gd, 'DBZ')
 vzf    = GET_VARIABLE_INDEX(gd, 'VZF')
 wz     = GET_VARIABLE_INDEX(gd, 'WZ')
 elec   = GET_VARIABLE_INDEX(gd, 'EX')
 cion   = GET_VARIABLE_INDEX(gd, 'CPIONINIT')
 muz     = GET_VARIABLE_INDEX(gd, 'MUPOSZC')
 
 ! check for extra diagnostic 3d arrays
  IF ( elec > wz .and.  elec - wz > 1 ) THEN ! elec is on and extra arrays exist
     xtra = wz + 1
     nxtra = elec - wz - 1
  ELSEIF ( u > wz + 1 ) THEN ! no elec and extra arrays exist
     xtra = wz + 1
     nxtra = u - wz - 1
  ELSE ! no extra arrays
     xtra = wz
     nxtra = 0
  ENDIF

 uinit  = GET_VARIABLE_INDEX(gd, 'UINIT')
 vinit  = GET_VARIABLE_INDEX(gd, 'VINIT')
 winit  = GET_VARIABLE_INDEX(gd, 'WINIT')
 piinit = GET_VARIABLE_INDEX(gd, 'PIINIT')
 kminit = GET_VARIABLE_INDEX(gd, 'KMINIT')
 sinit  = GET_VARIABLE_INDEX(gd, 'THINIT')

 gx     = GET_VARIABLE_INDEX(gd, 'XC')
 gy     = GET_VARIABLE_INDEX(gd, 'YC')
 gz     = GET_VARIABLE_INDEX(gd, 'ZC')
 precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')

!-----------------------------------------------------------------------------

   IF ( tvis5d .ne. thistory ) iskip = Max( 1, tvis5d/thistory )
   print*, 'iskip = ',iskip

!-----------------------------------------------------------------------------
! Main time step loop

   write(6,*) 'Welcome to comspectra! Available times are '
   DO it = 1, nt
     write(6,*) ' it, time = ',it,tarray(it)
   ENDDO
   
   write(6,*) 'Which time? (use index it)'
   read(6,*) it
       time     = tarray(it)

 100  CONTINUE

!   CALL GRID_READ_NETCDF( gd, prefix(1:length)//trim(number)//'.nc', tarray(it), 1 )

       IF ( it > 1 ) THEN
        IF ( historyinput == 1 ) THEN
          ncfile = prefix(1:length)//trim(number)//'.nc'
          CALL GRID_READ_NETCDF( gd, ncfile, tarray(it), 1 )
        ELSEIF ( historyinput .eq. 3  ) THEN
          write(stime,   '(i6.6)') time
          ncfile = prefix(1:length)//'.'//stime//'.nc'
          CALL GRID_READ_NETCDF( gd, ncfile, time, 1 )
        ENDIF
       ENDIF

        CALL SPECOUT(gd,                                 &
                    gd%var(u),  gd%var(uinit) ,         &            ! U,  UINIT
                    gd%var(v),  gd%var(vinit) ,         &            ! V,  VINIT
                    gd%var(w),  gd%var(winit) ,         &            ! W,  WINIT
                    gd%var(pi), gd%var(piinit),         &            ! PI, PIINIT
                    gd%var(km), gd%var(kminit),         &            ! KM, KINIT
                    gd%var(s),  gd%var(sinit),          &            ! S,  SINIT
                    gd%var(precip),                     &            ! PRECIP
                    gd%var(gx),                         &            ! XCNTR, XEDGE, DXC, DXE 
                    gd%var(gy),                         &            ! YCNTR, YEDGE, DYC, DYE
                    gd%var(gz),                         &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                                 &            ! DT
                    ugrid, vgrid,                       &            ! GRID MOTION
                    nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS
!                    st,                                 &
                     gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,gd%var(s)%index),  &
                    ft(1,1),ft(1,2),ft(1,3),ft(1,4),    &
                    gd%var(dbz), gd%var(vzf), gd%var(wz), gd%var(elec), gd%var(xtra) )


 
   write(6,*) 'Another time? (use index it) Type 0 to print times'
   read(6,*) it
   IF ( it == 0 ) THEN
     DO it = 1, nt
       write(6,*) ' it, time = ',it,tarray(it)
     ENDDO
     write(6,*) 'Which time? (use index it)'
     read(6,*) it
   ENDIF
   
       time     = tarray(it)
   
   IF ( it > 0 ) GOTO 100

 

   

   
!  DEALLOCATE(st)
  DEALLOCATE(ft)
 

!-----------------------------------------------------------------------------
! Print out date and time


 STOP
 END
 
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
 
 
   SUBROUTINE SPECOUT(gd,                 &
                    u,  uinit ,         &            ! U,  UINIT
                    v,  vinit ,         &            ! V,  VINIT
                    w,  winit ,         &            ! W,  WINIT
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    s,  sinit0,         &            ! S,  SINIT
                    precip,             &            ! PRECIP
                    gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                    gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                    gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                 &            ! DT
                    ugrid, vgrid,       &            ! GRID MOTION
                    nx, ny, nz, na,     &            ! NX,NY,NZ,NS
                    st,                 &            ! scalar array
                    dn, tt7, tt0, pn,   &
                    dbz, vzf, wz, elec, xtra )

   USE GRID_MODULE
   USE PARAM_MODULE,    only: pii,bcx,bcy,ng,isfcphys,luno
   USE MICRO_MODULE
   USE INDEX_MODULE
   USE COMMASMPI_MODULE
!   USE takcommon
   use takcommon,       only: lmax,lfmax,kfmax,iimax,kimax,lr100,ls250,dj
   use comm1,           only: rhof
   use comm2,           only: sr,sx,srf,fw,vw,bew,cdw,szr,szf,szi      &
     &                            ,sri,sfi,di,sxf,sxi,sh,ff,vf,fi,vi       
!-----------------------------------------------------------------------------
! VARIABLE DECLARATIONS
   
   implicit none

   TYPE(GRID), target :: gd

   integer :: nx, ny, nz, na
   integer :: nsmall    
   real    :: dt, ugrid, vgrid     

   TYPE(VARIABLE) :: u, uinit
   TYPE(VARIABLE) :: v, vinit
   TYPE(VARIABLE) :: w, winit
   TYPE(VARIABLE) :: pi, piinit
   TYPE(VARIABLE) :: km, kminit
   TYPE(VARIABLE) :: s(na), sinit0(2)     ! We are going to try and create arrays of the scalar variables needed here
   TYPE(VARIABLE) :: precip(4+neelec2d)
   TYPE(VARIABLE) :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE) :: dbz, vzf, wz
   TYPE(VARIABLE) :: xtra(nxtra)
   TYPE(VARIABLE) :: elec(neelec)
   

!   real :: preciptmp(-ng+1:nx+ng,-ng+1:ny+ng,4)
      
   real :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na)

   real :: dn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: pn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

   real :: tt7 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: tt0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

   real, allocatable :: tem1(:,:,:),tem2(:,:,:)
   real, allocatable ::  sinit(:,:)
   
   double precision :: total
      
!-----------------------------------------------------------------------------
! LOCAL variable declarations


   real    :: pb(nz)
   real    :: db(nz)

   integer :: is
   integer :: js
   integer :: ks
   
   real, parameter :: segmx = 20.0

   real    :: smin, smax


   real :: gtx(nx), gty(ny), gtz(nz)
   real :: gxt(nx,4), gyt(ny,4), gzt(nz,4)

   integer  time
      
   integer nstep1
   real start,tstat
      
   integer i,j,k,l, ix,jy,kz, n, n2, ia
   integer i1,j1,k1
   real    fac
      
   real dx,dy,dz
   real x,tx
   integer nstep, nstop, index


   integer idatime,hh,mm,ss
      

      
    integer len
    TYPE(ATTRIBUTE), pointer    :: microphys
    character(LEN = 15) :: micro
    logical    lice
    
    
    real :: hwdn, tmp, tmpg, tmpn, xdia, cno, tmpmx
    
    real :: pres,qwv,theta1,thetae

    
    integer :: ibc, jbc
    
    
!    real, parameter :: pii = 3.14159265359
    real :: cwch

        real chw,qr,z,alp,alpha1, rdi,pi1,vr,g1,zx,ze,nrx, sum
        integer ii


      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real dv, dax, day, daz, phiw, phie, phis, phin, phid, phiu

      double precision :: ztmp, ztmpr, ztmph, ztmphl, ztmpi, ctmpr, ctmph, ctmphl, ctmpi,ctmpc,xtmpi
      double precision :: ztmphl2,ztmph2,ztmpr2
      real dbzmax,dbzmin
      parameter ( dbzmin = -10.0 )
      integer :: ki,il
      real    :: xrtmp(ntakpd),xhtmp(ntakpd),xhltmp(ntakpd),xitmp(ntakid,ntakit)
      real    :: zitmp(ntakid,ntakit),zrtmp(ntakpd),zhtmp(ntakpd),zhltmp(ntakpd)
      real    :: scrtmp(ntakpd),schtmp(ntakpd),schltmp(ntakpd),scitmp(ntakid,ntakit)
      double precision :: schchg,schlchg,scichg,scrchg
      
      double precision :: massr,massh,masshl
      double precision :: qrmass1,qimass1,qhmass1,qhlmass1,qvmass1,qmasstot1,qmass,vtmass,dmwgt,d0,d0r,d0r2
      double precision :: masscum(45),qmass2,masscum2(45)


   real :: rdamelt(nch)
   
   logical :: doavg
   integer :: icount


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

   logical :: mpisep = .false.

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   logical :: debug_mpi = .false.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc Begin Execute
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    CALL GET_VARIABLE (gd, 'DX', dx)
    CALL GET_VARIABLE (gd, 'DY', dy)
    CALL GET_VARIABLE (gd, 'DZ', dz)


    CALL GET_ATTRIBUTE (gd, 'MICROPHYS', microphys)
    micro(:) = microphys%str(:)

    allocate ( sinit(nz,na) )

    sinit(:,:) = 0.0
      
    DO n = 1,2
     DO k = 1,nz-1
      sinit(k,n) = sinit0(n)%flt1d(k)
     ENDDO
    ENDDO
    
          pi1 = 4.0*atan(1.0)

    ixb = 1
    ixe = nx-1
    jyb = 1
    jye = ny-1
    kzb = 1
    kze = nz-1
    DO kz = 1,nz-1
      db(kz) = 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sinit(kz,lt))
      pb(kz) = 1.0e5*piinit%flt1d(kz)**3.509
      DO jy = 1,ny-1
        DO ix = 1,nx-1
          pn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509 - pb(kz)
          dn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/(287.04*s(lt)%flt3d(ix,jy,kz) )
          tt0(ix,jy,kz) = s(lt)%flt3d(ix,jy,kz)*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))
        ENDDO
      ENDDO
    ENDDO


      lice = .true.
      IF ( li .lt. 1 ) lice = .false.


         IF ( sr(1) == 0.d0 ) THEN
           call takinitbin
           write(0,*) 'call takinitbin, sx(1) = ',sx(1),sxf(1,1),sxf(1,2),lmax,lfmax
         ENDIF

       write(6,*) 'Enter i,j,k location: (negative j for vis5d loc, i=-1 for avg. rain)'
       read(6,*) i,j,k
       i1 = i
       IF ( i < 0 ) THEN
         doavg = .true.
         k1 = k
         ixe = nx-1
         jye = ny-1
       ELSE
         doavg = .false.
         ixe = 1
         jye = 1
       ENDIF
       IF ( j < 0 ) j = nyend + j

 50 CONTINUE


       IF ( doavg ) THEN
         icount = 0
         xrtmp(:) = 0.0

         IF ( i1 == -1 ) THEN
         k = k1
         DO j = 1,jye
           DO i = 1,ixe
           
            IF ( dbz%flt3d(i,j,k) > 15.0 ) THEN
              icount = icount + 1
              DO il = 1,lmax
                xrtmp(il) = xrtmp(il) + st(i,j,k,lnr+il-1)
              ENDDO
            
            
            ENDIF
           
           ENDDO
         ENDDO
         
         ELSEIF ( i1 == -2 ) THEN
         
         DO k = 1,k1
         DO j = 1,jye
           DO i = 1,ixe
           
            IF ( dbz%flt3d(i,j,k) > 15.0 .and. w%flt3d(i,j,k) < 0.0 ) THEN
              icount = icount + 1
              DO il = 1,lmax
                xrtmp(il) = xrtmp(il) + st(i,j,k,lnr+il-1)
              ENDDO
            ENDIF
           
           ENDDO
         ENDDO
         ENDDO
         
         
         ENDIF
         
         DO il = 1,lmax
            xrtmp(il) = xrtmp(il)/icount
         ENDDO
         
       write(6,*) 'il, sr, Dr(mm), nr, Nr(D), m_r(D)'
       DO il = 1,lmax
         write(6,'(i2,5(2x,1pe12.5))'), il, sr(il), 20.*sr(il),xrtmp(il),   &
     &          xrtmp(il)/(2.*sr(il)/dj), &
     &          sx(il)*xrtmp(il)/(2.*sr(il)/dj)
       ENDDO
         
         
       
       ENDIF
       
       IF ( .not. doavg ) THEN
       write(6,*) 'Getting spectra for ',i,j,k
       
       ztmpr = 0.d0
       ztmph = 0.d0
       ztmphl = 0.d0
       ztmpr2 = 0.d0
       ztmph2 = 0.d0
       ztmphl2 = 0.d0

       ctmpc = 0.d0
       ctmpr = 0.d0
       ctmph = 0.d0
       ctmphl = 0.d0
       
       massr = 0
       massh = 0
       masshl = 0
       
       xrtmp(:) = 0.0
       zrtmp(:) = 0.0
       DO il = 1,lmax
         xrtmp(il) = st(i,j,k,lnr+il-1)
         zrtmp(il) = 10.0*Log10(Max(1.d-10, 1.d6*xrtmp(il)*szr(il)))
         IF ( il < lr100 ) THEN
           ctmpc = ctmpc + xrtmp(il)
         ELSE
           ctmpr = ctmpr + xrtmp(il)
           ztmpr = ztmpr + xrtmp(il)*szr(il)
           massr = massr + xrtmp(il)*sx(il)
         ENDIF
       ENDDO
       DO il = 1,lfmax
         xhtmp(il) = st(i,j,k,lnh+il-1)
         xhltmp(il) = st(i,j,k,lnhl+il-1)
         IF ( il > ls250 ) THEN
           ztmph  = ztmph  + xhtmp(il)*szf(il,1)
           ztmphl = ztmphl + xhltmp(il)*szf(il,2)
           massh  = massh  + xhtmp(il)*sxf(il,1)
           masshl = masshl + xhltmp(il)*sxf(il,2)
           ctmph  = ctmph  + xhtmp(il)
           ctmphl = ctmphl + xhltmp(il)
         ENDIF
         zhtmp(il) = 10.0*Log10(Max(1.d-10, 1.d6*0.224*xhtmp(il)*szf(il,1)))
         zhltmp(il) = 10.0*Log10(Max(1.d-10, 1.d6*0.224*xhltmp(il)*szf(il,2)))
       ENDDO

       ctmpi = 0.d0
       xtmpi = 0.d0
       ztmpi = 0.d0
       DO il=1,ntakit
       DO ii=1,ntakid
         ia = lni + ii - 1 + (il-1)*ntakid
         xitmp(ii,il) = st(i,j,k,ia)
         zitmp(ii,il) = 10.0*Log10(Max(1.d-10, 1.0d12*0.224 * xitmp(ii,il)*(((6.0/pii))*sxi(ii,il))**2))
         ztmpi = ztmpi+xitmp(ii,il)*(((6.0/pii))*sxi(ii,il))**2
         xtmpi = xtmpi+xitmp(ii,il)*sxi(ii,il)
         ctmpi = ctmpi + xitmp(ii,il)
       ENDDO
       ENDDO

       IF ( lscr > 1 ) THEN
       DO il = 1,lmax
         scrtmp(il) = st(i,j,k,lscr+il-1)
       ENDDO
        scrtmp(lmax+1:ntakpd) = 0
       ENDIF

       IF ( lsch > 1 .and. lschl > 1 ) THEN
       DO il = 1,lfmax
         schtmp(il) = st(i,j,k,lsch+il-1)
         schltmp(il) = st(i,j,k,lschl+il-1)
       ENDDO
       ENDIF

       IF ( lsci > 1 ) THEN
       DO il=1,ntakit
       DO ii=1,ntakid
         ia = lsci + ii - 1 + (il-1)*ntakid
         scitmp(ii,il) = st(i,j,k,ia)
       ENDDO
       ENDDO
       ENDIF

          qmass = 0.0
          qmass2 = 0.0
          vtmass = 0.0
          dmwgt = 0.0d0
          DO il = 16,ntakrd
            qmass = qmass + sx(il)*xrtmp(il)
            masscum(il) = qmass ! cumulative mass
            dmwgt = dmwgt + (2.d0*sr(il))**3 * sx(il)*xrtmp(il)
          ENDDO

          qmass2 = 0.0
          DO il = 28,ntakrd
            qmass2 = qmass2 + sx(il)*xrtmp(il)
            masscum2(il) = qmass2 ! cumulative mass
          ENDDO
          
          IF ( qmass*1.e-3*1.e6 > 1.e-8*dn(i,j,k) ) THEN

          n = 0
          
            DO il = 16,ntakrd
              IF ( qmass*0.5d0 < masscum(il) ) THEN
                n = il
                EXIT
              ENDIF
            ENDDO

          n2 = 0
          
            DO il = 28,ntakrd
              IF ( qmass2*0.5d0 < masscum2(il) ) THEN
                n2 = il
                EXIT
              ENDIF
            ENDDO
          
            
            IF ( n <= 16 ) THEN
             d0r = 0.0
            ELSEIF ( n == 1 .or. n == 15 ) THEN
             d0r  = 0.0 !  sr(n)
            ELSE
            ! n is index of masscum just above d0
               d0 =  sr(n-1) + (0.5*qmass - masscum(n-1))/(masscum(n) - masscum(n-1))*(sr(n) - sr(n-1)) ! radius
               d0r = 2.*d0*10. ! *2 for diam, *10 for units of mm
          !     IF ( d0 < 0.0 ) THEN
          !       write(0,*) 'd0,n,qmass,masscum = ',d0,n,qmass,masscum(n-1),masscum(n)
          !     ENDIF
            ENDIF
            
            IF ( n2 <= 28 ) THEN
             d0r2 = 0.0
            ELSE
            ! n is index of masscum just above d0
               d0 =  sr(n2-1) + (0.5*qmass2 - masscum2(n2-1))/(masscum2(n2) - masscum2(n2-1))*(sr(n2) - sr(n2-1)) ! radius
               d0r2 = 2.*d0*10. ! *2 for diam, *10 for units of mm
            ENDIF
          
          ELSE
            d0r = 0.0
          ENDIF

       
       ztmp = Max(1.d-10, 1.d6*(1.d6*0.224*ztmpi + ztmpr + 0.224*ztmph + 0.224*ztmphl) )

       ztmpr2  = ztmpr
       ztmph2  = ztmph
       ztmphl2 = ztmphl

       ztmpr  = Max(1.d-10, 1.0d6*ztmpr )
       ztmph  = Max(1.d-10, 1.0d6*ztmph*0.224 )
       ztmphl = Max(1.d-10, 1.0d6*ztmphl*0.224)
!       ztmpi  = Max(1.d-10, 1.0d12*ztmpi*0.224 )
       
       write(6,*) 'dbz,qr,qh,qhl = ',dbz%flt3d(i,j,k),st(i,j,k,lr),st(i,j,k,lh),st(i,j,k,lhl)
       write(6,*) 'massr,h,hl = ',1000.*massr,1000.*massh,1000.*masshl
       tmp = 0.0
       IF ( ctmpi > 0 ) tmp = xtmpi/ctmpi
       write(6,*) 'cimas(g), qi(g/kg), cimasave(mg) = ',xtmpi,xtmpi*1.e6/dn(i,j,k), 1000*tmp
       write(6,*) 'dBZr,dBZh,dBZhl,dbZ = ',10.0*Log10(ztmpr), 10.0*Log10(ztmph), 10.0*Log10(ztmphl),10.0*Log10(ztmp)
       write(6,*) 'Zi,dBZi = ',1.0d12*ztmpi,10.0*Log10(Max(1.d-10, 1.0d12*ztmpi*0.224 ))
       write(6,*) 'Nr,Nh,Nhl,Ni,Nc = ', 1.e6*ctmpr, 1.e6*ctmph, 1.e6*ctmphl,1.e3*ctmpi,ctmpc
       write(6,*) 'Zr,Zh,Zhl = ',ztmpr2,ztmph2,ztmphl2
       write(6,*) 'lr100,ls250 = ',lr100,ls250
       write(6,*) 'median vol. diam d0r,d0r2 = ',d0r,d0r2
       write(6,'(a,5(2x,1pe12.5))') 'T,theta,pp,pres,qv = ', tt0(i,j,k), s(lt)%flt3d(i,j,k), pn(i,j,k), pn(i,j,k)+pb(k), s(lv)%flt3d(i,j,k)
              
       write(6,*) 'il, sr, Dr(mm), nr, nh, nhl, Nr(D), Nh(D), Nhl(D), m_r(D), m_h(D), m_hl(D)'
       DO il = 1,lfmax
         write(6,'(i2,11(2x,1pe12.5))'), il, sr(il), 20.*sr(il),xrtmp(il),xhtmp(il),xhltmp(il),   &
     &          xrtmp(il)/(2.*sr(il)/dj),xhtmp(il)/(2.*sr(il)/dj),xhltmp(il)/(2.*sr(il)/dj), &
     &          sx(il)*xrtmp(il)/(2.*sr(il)/dj),sxf(il,1)*xhtmp(il)/(2.*sr(il)/dj),sxf(il,2)*xhltmp(il)/(2.*sr(il)/dj)
       ENDDO

       write(6,*) 'il, sr, Dr(mm), Zr, Zh, Zhl, dbZr, dbZh, dbZhl'
       DO il = 1,lfmax
         write(6,'(i2,8(2x,1pe12.5))'), il, sr(il), 20.*sr(il) ,xrtmp(il)*szr(il),xhtmp(il)*szf(il,1),xhltmp(il)*szf(il,2),zrtmp(il),zhtmp(il),zhltmp(il) 
       ENDDO

       write(6,*) 'il, sri, xi1'
       DO il = 1,iimax
         write(6,'(i2,6(2x,1pe12.5))'), il, sri(il), (xitmp(il,ii),ii=1,ntakit)
       ENDDO

       write(6,*) 'il, zi1'
       DO il = 1,iimax
         write(6,'(i2,5(2x,1pe12.5))'), il,  (zitmp(il,ii),ii=1,ntakit)
       ENDDO

!       write(6,*) 'il, sri, sxi'
!       DO il = 1,iimax
!         write(6,'(i2,6(2x,1pe12.5))'), il, sri(il), (sxi(il,ii),ii=1,ntakit)
!       ENDDO

       IF ( lscr > 1 ) THEN
       
       
       schchg = 0
       schlchg = 0
       scichg = 0
       scrchg = 0
       
       write(6,*) 'il, sr, scr, sch, schl'
       DO il = 1,lfmax
         write(6,'(i2,7(2x,1pe12.5))'), il, sr(il),scrtmp(il),schtmp(il),schltmp(il)
         schchg = schchg + schtmp(il)
         schlchg = schlchg + schltmp(il)
         scrchg = scrchg + scrtmp(il)
       ENDDO

       write(6,*) 'il, sri, sci1'
       DO il = 1,iimax
         write(6,'(i2,6(2x,1pe12.5))'), il, sri(il), (scitmp(il,ii),ii=1,ntakit)
         scichg = scichg + Sum( scitmp(il,1:ntakit) )
       ENDDO
       
        IF ( lscpi > 1 ) THEN
          write(6,*) 'pos/neg ion concentrations = ',st(i,j,k,lscpi), st(i,j,k,lscni)
        ENDIF
        
        write(6,*) 'scnet = ',elec(iscnet)%flt3d(i,j,k) 
        write(6,*) 'scr, sch, schl, sci = ',scrchg, schchg, schlchg, scichg
       
       ENDIF
       
       ENDIF ! .not. doavg

       write(6,*) 'Enter i,j,k location: (negative j for vis5d loc, i<-10 to return)'
       read(6,*) i,j,k
       IF ( j < 0 ) j = nyend + j
       IF ( i > -10 ) GOTO 50


      RETURN
      END SUBROUTINE SPECOUT

