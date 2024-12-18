!-----------------------------------------------------------------------------
!
! COMMAS-1.0
!
! March 2006
!
!-----------------------------------------------------------------------------

 PROGRAM COMTRAJ
 
! for generating trajectories/soundings through a single history time

#ifdef NAG
  USE F90_UNIX_PROC
#endif

  USE GRID_MODULE
  USE GRIDIO_MODULE
  USE FILE_MODULE
  USE CLINE_MODULE
  USE COMMASMPI_MODULE
  USE CPUTIME_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
  USE VIS5D_MODULE
  USE TRAJ_MODULE
  USE INDEX_MODULE
  USE FOLLOW_MODULE

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

  real    :: time_real

  logical :: file_exist = .false.
  logical :: io_flag, lstt
  logical :: hack = .false.

  integer i,j,k, it, itime, nt, count

  character(LEN=100)  :: filenm
  character(LEN=120) :: ncfile

  character(len=6) stime
  character(LEN=4) number1,number2
  character(LEN=8) number
  character(LEN=2) dup
  character(len=3) column

  character(LEN = 255):: v5dfldstmp

  integer  :: istat
  integer  :: iskip = 1

  character(LEN=1) trjnum
  
  integer :: vertical_only = 0

  NAMELIST /comtraj_param/               &
                      vertical_only
                      

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

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'COMTRAJ:  PROBLEM READING NAMELIST'

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMTRAJ:  PROBLEM READING NAMELIST'

  IF( .not. READ_NAMELIST(run_file(1:length),'GRIDN') ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST GRIDN'
  ENDIF

! Read in 10-ice parameters

  IF( (.not. READ_NAMELIST(run_file(1:length),'TRAJECTORIES'))  ) THEN
    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST TRAJECTORIES'
  ENDIF

  IF( .not. READ_NAMELIST(run_file(1:length),'ICE10_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF

  open(15,file=run_file(1:length),status='old',form='formatted')
  rewind(15)
  read(15,NML=comtraj_param,iostat=istat)
  close(15)
      IF ( istat .ne. 0 ) THEN
        write(0,*) 'READ_NAMELIST: comtraj_param namelist not found -- using default values'
      ENDIF



  IF ( historyinput .eq. -1 ) historyinput = historyoutput

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

  CALL COMMASMPI_BOUNDS(nx,ny,nz)


   write(number1, '(a,i3.3)') '.', member
   number = number1



!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   CALL STRING_LIMITS(prefix, ibeg, iend)
   v5dfilename = prefix(ibeg:iend)//trim(number)//'.v5d'
   
   IF ( nyend > 2 ) THEN
     onedoutput = 0
   ELSE
     write(column,'(i3.3)') onedoutput
     onedfilename = prefix(ibeg:iend)//trim(number)//'.'//column//'.1d'
   ENDIF


!-----------------------------------------------------------------------------
! ELSE, READ IN THE DATA


    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

     IF ( itraj > 0 ) THEN
       start = time_traj1
       stop = time_traj2
     ENDIF


!    IF ( start .lt. 0 ) start = 0
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

    itime = -1

    DO it = 1,nt
      IF ( tarray(it) == start ) THEN
        itime = it
      ENDIF
    ENDDO
    
    IF ( itime == -1 ) THEN
     write(0,*) 'Time not found in history file!'
     write(0,*) 'Was looking for time = ',start,' but available times are'
      DO it = 1,nt
        write(0,*) tarray(it)
      ENDDO
      STOP
    ENDIF
    
    CALL GRID_READ_NETCDF( gd, ncfile, tarray(itime), nx_or_nxend=.true. )

    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)

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
    
    


    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start))   ) write(6,*) 'COMMAS:  Problem setting TIME'  



   count = 0
   
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
  CALL GET_VARIABLE(gd, 'ISFCPHYS', isfcphys)
  CALL GET_VARIABLE(gd, 'IMURAIN',  imurain)
   write(0,*) 'COMTRAJ: isfcphys, IPCONC = ',isfcphys,ipconc,i

!    write(6,*) 'COMTRAJ: dt = ',dt

! for balloon trajectories:
!   IF ( chgavex <= 0 ) chgavex = ng
   IF ( chgavex <= 0 ) chgavex = Max(3, Nint(750./dx) )
    write(6,*) 'COMTRAJ: chgavex = ',chgavex
   ! set up gamma function lookup table for setvtz
   call makegmoi()
  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  
  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COMMAS:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COMMAS:  Problem setting NG'

!  write(*,*) 'ng = ',ng

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


  IF ( itraj > 0 ) THEN
      filenm = prefix
      write(stime,   '(i6.6)') start
      write(stimecomtraj,   '(i6.6)') start
      icomtraj = 1
      IF ( vertical_only > 0 ) icomtraj = 2

    IF ( my_rank == 0 .and. itraj >= 2 ) THEN
     IF ( ntrajtype == 1 ) THEN
     OPEN(unit=trj_unit, file=filenm(1:length)//stime//trim(number)//'.traj', status='new',form='formatted')
     ELSE
      DO i = 1,ntrajtype
       write(trjnum,'(i1)') i
       OPEN(unit=trj_unit+i-1, file=filenm(1:length)//stime//trim(number)//'.traj'//trjnum, status='new',form='formatted')
      ENDDO
     ENDIF
    ENDIF
!    write(0,*) 'call traj_init'
    CALL TRAJ_INIT(dt,nxend,nyend,nzend,nx,ny,nz,ng,filenm,number,start)
!    write(0,*) 'done traj_init'
  ENDIF
#ifdef MPI
   CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif

#ifdef MPI
!   CALL PRINT_MPI( gd )
#else
!   CALL PRINT( gd )
#endif


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


!-----------------------------------------------------------------------------
! SUBTRACT OUT GRID MOTION
    IF ( ugrid .ne. 0.0 ) THEN
       IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(6,*) 'COMMAS:  Problem setting UGRID'  
    ENDIF
    IF ( vgrid .ne. 0.0 ) THEN
       IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(6,*) 'COMMAS:  Problem setting VGRID'
    ENDIF

        CALL FOLLOW(gd, -1.0, nxend)

!-----------------------------------------------------------------------------
! Main time step loop

         nt = (stop - start)/dt + 1
         write(0,*) 'COMTRAJ: nt,dt,start,stop = ',nt,dt,start,stop

   DO it = 1, nt

!----------------------------------------------------------------------
! Time and location of the grid
   time = start + (it-1)*dt
   time_real = start + (it-1)*dt
!       IF ( historyoutput == 3 ) THEN 
!         time = start + (it - 1)*thistory
!       ELSE
!         time = tarray(it)
!       ENDIF

!       write(0,*) 'it loop: it,time1 = ',it,time


!       IF ( time .gt. stop ) EXIT
       
!       IF ( time .ge. start ) THEN

!       write(6,*) 'it loop: it,time = ',it,time

!       IF ( it > 1 ) THEN
!        IF ( historyinput == 1 ) THEN
!          ncfile = prefix(1:length)//trim(number)//'.nc'
!          CALL GRID_READ_NETCDF( gd, ncfile, tarray(it), 1 )
!        ELSEIF ( historyinput .eq. 3  ) THEN
!          write(stime,   '(i6.6)') time
!          ncfile = prefix(1:length)//'.'//stime//'.nc'
!          CALL GRID_READ_NETCDF( gd, ncfile, time, 1 )
!        ENDIF
!         
!         
!       ENDIF


!        CALL GET_VARIABLE(gd, 'TIME',       itmp)
        
!        write(6,*) 'time from netcdf is ',itmp

!    IF( .not. SET_VARIABLE(gd,'TIME',      time)   ) write(6,*) 'COMMAS:  Problem setting TIME'  

        iv5dwritten(:) = 0
        lv5dwrite = .true.
        v5dtimes(itv5d) = time
 

!-----------------------------------------------------------------------
! Call SOLVER to make a time step


!-----------------------------------------------------------------------------
! DTREND the pressure field
!
!
!-----------------------------------------------------------------------------
! PRINTING calls

    
!       CALL PRINT( gd )


!-----------------------------------------------------------------------------
! Vis5D OUTPUT
!

        
     CALL trajdriver(gd,                                 &
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
                 lat, lon,                           &            ! Latitude and Longitude
                 microphys,                          &
                 nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS (mpi = nxt,nyt,nzt)
                 io_flag,                            &            ! IO_FLAG
                 ft(1,1),ft(1,2),ft(1,3),            &
                 gd%var(dbz), gd%var(vzf), gd%var(wz),gd%var(xtra),  &
                 gd%var(elec), gd%var(cion),         &
                 gd%var(muz),                        &
                 x_sw_loc,y_sw_loc,                  &
                 time, time_real, tstat, lstt, stop)
        
       

!        write(6,*) 'done V5DOUT'

!        iv5dwritten(:) = 0

   

!        ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO 
   
!  DEALLOCATE(st)
  DEALLOCATE(ft)
 

!-----------------------------------------------------------------------------
! Print out date and time


 STOP
 END
 
 
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
       subroutine trajdriver(gd,                 &
                    u,  uinit ,         &            ! U,  UINIT
                    v,  vinit ,         &            ! V,  VINIT
                    w,  winit ,         &            ! W,  WINIT
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    s,  sinit,          &            ! S,  SINIT
                    precip,             &            ! PRECIP
                    gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                    gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                    gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                 &            ! DT
                    ugrid, vgrid,       &            ! GRID MOTION
                    lat, lon,           &            ! Latitude and Longitude
                    microphys,          &            ! microphysical scheme
                    nx, ny, nz, ns,     &            ! NX,NY,NZ,NS
                    io_flag,            &            ! IO_FLAG
                    t0, t1, t2,         &
                    dbz, vzf, wz, xtra, &
                    elec, cion,         &
                    muz,                &
                    x_sw_loc,y_sw_loc,  &                                   
                    time, time_real, tstat, lstt, tstop)

   USE GRID_MODULE
   USE CPUTIME_MODULE
   USE MICRO_MODULE
   USE PARAM_MODULE
   USE INDEX_MODULE
   USE FORCE_MODULE
   USE TRAJ_MODULE
   USE INIT_MODULE, only: inhom,bogusvortex,timint ! ,HURRFORCE
   USE COMMASMPI_MODULE

#ifdef MPI
    USE mpi
#endif

!-----------------------------------------------------------------------------
   
   implicit none
 
#ifdef MPI
!  INCLUDE "mpif.h"
#endif

!-----------------------------------------------------------------------------
! GRID DEFINITIONS

   TYPE(GRID)         :: gd 

   character(LEN = *) :: microphys
   integer            :: mscheme,nscalar ! Needed for Milbrandt-Yau scheme
   integer            :: nx, ny, nz, ns
   real               :: dt
   real               :: ugrid, vgrid     
   real               :: lat, lon
   logical            :: io_flag

   real unorm
   real vnorm

   TYPE(VARIABLE)     :: u, uinit
   TYPE(VARIABLE)     :: v, vinit
   TYPE(VARIABLE)     :: w, winit
   TYPE(VARIABLE)     :: pi, piinit
   TYPE(VARIABLE)     :: km, kminit
   TYPE(VARIABLE)     :: s(ns), sinit(2)     ! We are going to try and create arrays of the scalar variables needed here
   TYPE(VARIABLE)     :: precip(nprecip+neelec2d)
   TYPE(VARIABLE)     :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE)     :: dbz, vzf, wz
   TYPE(VARIABLE)     :: xtra(nxtra)
   TYPE(VARIABLE)     :: elec(neelec)
   TYPE(VARIABLE)     :: cion(2)
   TYPE(VARIABLE)     :: muz(4)


   real :: preciptmp(-ng+1:nx+ng,-ng+1:ny+ng,nprecip+neelec2d)

   double precision :: preciptot(4), preciptotall(4)


   real :: t0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: t1 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real :: t2 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real, allocatable, save :: dpdt(:,:,:)
   
   real :: z1d4(nzend,4)
   
   integer :: tstat
   logical :: lstt
   real    :: x_sw_loc,y_sw_loc
      
! Local variables
      
   logical, parameter :: debugsolver = .false.
   logical, parameter :: truetime = .false. ! turns on barriers to get timing of each section for the slowest process
   integer :: loop, ge, i, j, k, n, m, nstep, atype, ntimestep
   real    :: dx, dy, dz
   real    :: dx1, dy1, dz1, dyl, mlen
   real    :: dts, sdt
   real    :: den(-ng+1:nz+ng,2), rrp(-ng+1:nz+ng,2), rrm(-ng+1:nz+ng,2)
   double precision :: deninv(-ng+1:nz+ng,2)
   integer :: nsrk,nsfwd
      
   real    :: gtx(-ng+1:nx+ng), gty(-ng+1:ny+ng), gtz(-ng+1:nz+ng)
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4)
   integer :: nst = 1                                ! dummy timestep; used to flip xy|yx ordering in Crowley advection
   save nst
   integer :: iadiv
   real    :: dt1
   integer :: nsub
   real    :: tmp

   integer :: ix,jy,kz,ia
   integer :: i1,i2, j1,j2, k1,k2
   integer :: is, js, ks
   integer :: im1, jm1, ip1, jp1
   real    :: wmax, wmin, vmax, a
   integer :: time, tstop
   real    :: time_real

   real    :: q(2*lqmx), pii_total, t_total
   real    :: qtodbz

   real, allocatable, save :: sbase(:,:)

   real wzz

   logical, parameter :: filter = .false.
   integer, parameter :: ihole  = 2

! Flag to dump out dpdt stats for noise estimates

   logical, parameter :: print_noise = .false.

! microphysics flags for number of sub-steps and minimum dbz

      
   real, parameter    :: mindbz = 0.0
!   integer, parameter :: nphys = 1
!   iuvwadv is defd in param_module, read in namelist ! 1 = normal; 2 = box scheme; 3 = turns off wind update for pure scalar advection
   

   integer, save  :: ifirst = 0
   
   integer, save :: llen
   
   integer, save :: ib,ie,jb,je,kb,ke,ni,nj,nk
   
   integer :: nrain
   
   real, allocatable, save :: dz3d(:,:,:)
   
   real, allocatable, save :: pb(:), db(:)
   
   integer :: tflag ! flag to store turb mixing tendency (for supersaturation tendency)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   integer :: westward_tag, eastward_tag
   integer :: northward_tag, southward_tag
   integer :: downward_tag, upward_tag

   logical :: debug_mpi = .false.

   integer       :: nampi, nb

! Takahashi microphysics
          real   (kind=8) :: a1,a2,a3,parta,pertrb                            &
     &                       ,term,term3,term4,term5,terma,work,worka          &
     &                       ,cap,farm,tempe,dummy,fbuf,pbuf

      real :: temq
      integer :: l

   CALL GET_VARIABLE (gd, 'DX', dx)
   CALL GET_VARIABLE (gd, 'DY', dy)
   CALL GET_VARIABLE (gd, 'DZ', dz)

   DO i = 1,4
     gxt(:,i) = gx(i)%flt1d(:)
     gyt(:,i) = gy(i)%flt1d(:)
     gzt(:,i) = gz(i)%flt1d(:)
   ENDDO


   IF ( .not. allocated( sbase ) ) THEN

   allocate ( sbase(-ng+1:nz+ng,ns) )
   sbase(:,:) = 0.0

   DO n = 1,2
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
#else
     DO k = 1,nz
#endif
       sbase(k,n) = sinit(n)%flt1d(k)
     ENDDO
   ENDDO

   DO n = 3,ns
#ifdef MPI
     kzb = -ng+1
     kze = ktile+ng
     if (kzend .eq. nzend) kze = kzend-kzbeg+1+ng
     
     do k = kzb,kze
#else
     DO k = 1,nz
#endif
      sbase(k,n) = s(n)%base1d(k)
     ENDDO
   ENDDO

   ENDIF


!
! Do trajectories if needed
!
  IF ( itraj > 0 ) THEN
    CALL cld_cpu('TRAJECTORIES')
!    write(0,*) 'solver1: qh,nh = ',lh,lnh,gd%xyz3d%flt4d(61,50,72,km%index+lh),gd%xyz3d%flt4d(61,50,72,km%index+lnh)
!    write(0,*) 'call traj: microphy = ',microphys,dx,dy,time,time_real
    CALL TRAJ(nx,ny,nz,ns,dt,gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,km%index+1),   &
     &        sbase,u%flt3d,v%flt3d,w%flt3d,t0,t1,t2,piinit%flt1d, &
     &        pi%flt3d, km%flt3d, dbz%flt3d, elec,  &
     &        gxt,gyt,gzt,time,time_real,    &
     &        uinit%flt1d,vinit%flt1d,ugrid,vgrid,microphys,dx,dy)

!    write(0,*) 'solver2: qh,nh = ',st(20,20,30,lh),st(20,20,30,lnh)
!    write(0,*) 'done traj'
    CALL cld_cpu('TRAJECTORIES') 
  ENDIF



    RETURN
    END
 
 !-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

 
