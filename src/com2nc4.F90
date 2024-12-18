!-----------------------------------------------------------------------------
!
! COMMAS-1.0
!
! March 2006
!
!-----------------------------------------------------------------------------

 PROGRAM COM2NC4

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
  USE VIS5D_MODULE
  USE INDEX_MODULE
  USE COMMASMPI_MODULE

  implicit none

#ifdef MPI
  INCLUDE "mpif.h"
#endif

!-----------------------------------------------------------------------------
! GRID DEFINITION
 
  TYPE(GRID) :: gd, ge

!-----------------------------------------------------------------------------
! COMMAND LINE VARIABLES

  character(LEN = 120) :: run_file
  character,DIMENSION(120)  :: run_file_array      !needed for mpi process comm
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
  integer :: wz    
  integer :: elec
  integer :: cion
  integer :: muz
  integer :: ipres,itemp
   
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
  integer :: llen

  logical :: file_exist = .false.
  logical :: io_flag, lstt, lfix
  logical :: hack = .false.

  integer i,j,k, it, nt, count

  character(LEN=4) number1,number2
  character(len=6) stime
  character(LEN=8) number
  character(LEN=2) dup

  character(LEN = 255):: v5dfldstmp
  character(LEN=120) :: ncfile,ncfilein,ncfileout
  character(len=1)   :: sversion
  
  double precision :: rhoair

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

  logical :: mpisep = .false.

#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   integer :: mpi_status(MPI_Status_size)
   integer :: token

!  integer :: mpi_error_code

   logical :: debug_mpi = .TRUE.
#endif


!-----------------------------------------------------------------------------
! MPI STARTUP

#ifdef MPI
  CALL COMMASMPI_STARTUP()
#endif

!-----------------------------------------------------------------------------
! Need 3D.RUN file in order to get run parameters - get from the command line

  DO i = 1, LEN(run_file)
    run_file(i:i) = " "
  END DO

#ifdef MPI
  IF (my_rank == 0) THEN
#endif
    IF( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN

      write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  NEED 3D.RUN FILENAME'
      write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  EXITING RUN!...'
      call exit(1)

    ELSE

     CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)

!     write(6,*) 'COMMAS:  Input 3D.RUN file is:  ', run_file(1:length)

     IF( status .ne. 0 ) THEN

      write(0,*) 'COMMAS:  COULD NOT RETRIEVE 1ST COMMAND LINE ARG:  EXITING RUN!...', status
      write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ!'
      call exit(1)

     ENDIF

    ENDIF

#ifdef MPI
    DO i = 1, length
      run_file_array(i) = run_file(i:i)
    END DO

  END IF !! (my_rank == 0)

  CALL MPI_Bcast(length, 1, MPI_INTEGER, &
 &         0, MPI_COMM_WORLD, mpi_error_code)
  CALL MPI_Bcast(run_file_array, length, MPI_CHARACTER, &
 &         0, MPI_COMM_WORLD, mpi_error_code)

  DO i = 1, length
    run_file(i:i) = run_file_array(i)
  END DO
#endif

  
!-----------------------------------------------------------------------------
! Read in 3d.run file...

  INQUIRE(file=run_file(1:length), exist=file_exist)

  IF( .NOT. file_exist ) THEN

    write(0,*) 'COMMAS:  INPUT 3D.RUN FILE:  ', run_file(1:length), ' DOES NOT EXIST!!!'
    write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ, EXITING'
    stop

  ENDIF

! READ IN RUN MODEL NAMELIST

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

   write(number1, '(a,i3.3)') '.', member
   number = number1

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  netcdfversion = 4

  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
  
  
! READ IN GRID PARAMETERS





!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   CALL STRING_LIMITS(prefix, ibeg, iend)


!-----------------------------------------------------------------------------
! ELSE, READ IN THE DATA


    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

    IF ( start .lt. 0 ) start = 0
     write(stime,   '(i6.6)') start
    IF ( historyinput == 3 ) THEN
      ncfilein = prefix(1:length)//'.'//stime//'.nc'
    ELSE
      ncfilein = prefix(1:length)//trim(number)//'.nc'
    ENDIF
    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',trim(ncfilein)

  

    IF ( historyinput == 3 ) THEN
       
       nt = (stop - start)/thistory + 1
       allocate( tarray(nt) )
       DO it = 1,nt
        tarray(it) = start + (it-1)*thistory
        write(6,*) 'it, time = ',it,tarray(it)
       ENDDO
       
    ELSE
    
      CALL GRID_INFO_NETCDF( ncfilein, nt, microphys)
      allocate( tarray(nt) )
    
      CALL GRID_INFO_NETCDF( ncfilein, nt, microphys, tarray)
    
    ENDIF
    
    CALL GRID_READ_NETCDF( gd, ncfilein, tarray(1), nx_or_nxend=.true.  )

    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION '

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


!
! Create the NC4 file:
!

!  CALL GET_VARIABLE(gd, 'NXEND',    nx)
!  CALL GET_VARIABLE(gd, 'NYEND',    ny)
!  CALL GET_VARIABLE(gd, 'NZEND',    nz)

!  IF( .not. SET_VARIABLE(gd,'NX',        nx)            ) write(6,*) 'COM2NC4:  Problem setting NX'
!  IF( .not. SET_VARIABLE(gd,'NY',        ny)            ) write(6,*) 'COM2NC4:  Problem setting NY'
!  IF( .not. SET_VARIABLE(gd,'NZ',        nz)            ) write(6,*) 'COM2NC4:  Problem setting NZ'
              
  CALL GET_VARIABLE(gd, 'NXEND',    nx)
  CALL GET_VARIABLE(gd, 'NYEND',    ny)
  CALL GET_VARIABLE(gd, 'NZEND',    nz)

  IF( .not. SET_VARIABLE(gd,'NX',       nx)      ) write(6,*) 'COM2NC:  Problem setting NX'
  IF( .not. SET_VARIABLE(gd,'NY',       ny)      ) write(6,*) 'COM2NC:  Problem setting NY'


    IF( .not. SET_ATTRIBUTE(gd,'NX',  nx)) write(6,*) 'COMMAS:  Problem setting V5DFIELDS'
    IF( .not. SET_ATTRIBUTE(gd,'NY',  ny)) write(6,*) 'COMMAS:  Problem setting V5DFIELDS'

  ncxe = nx
  ncye = ny
  ncze = nz
  
          nxbeg = 1
          nxend = nx
          nybeg = 1
          nyend = ny
          nzbeg = 1
          nzend = nz
              ixbeg = 1
              jybeg = 1
              kzbeg = 1
              ixend = nxend
              jyend = nyend
              kzend = nzend
    write(sversion,'(1i1)' ) netcdfversion
    IF ( historyoutput == 3 ) THEN
      ncfileout = prefix(1:length)//'.'//stime//'.nc'//sversion
    ELSE
      ncfileout = prefix(1:length)//trim(number1)//'.nc'//sversion
      CALL GRID_DEFINE_NETCDF(gd, ncfileout)
      CALL GRID_WRITE_NETCDF( gd, time, ncfileout)
    ENDIF

!       IF ( netcdfversion .eq. 4 ) THEN
!      CALL GRID_DEFINE_NETCDF(gd, ncfileout)
!      CALL GRID_WRITE_NETCDF( gd, time, ncfileout)
!       ELSEIF ( netcdfversion .eq. 3 ) THEN
!      CALL GRID_DEFINE_NETCDF(gd, ncfileout)
!      CALL GRID_WRITE_NETCDF( gd, time, ncfileout)
!       ENDIF



!--------------------------------------------------------------------------------------
! Open VIS5D file  

!  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
    
    IF ( start .gt. tarray(nt) ) THEN
      write(0,*) 'Starting time is greater than the last record time! STOP!'
      STOP
    ENDIF
    
!    start = Max(start, tarray(2) ) 
    
!    DO it = nt,1,-1
!      IF ( start .le. tarray(it) ) v5dtimes(1) = tarray(it)
!    ENDDO
    
!    v5dtimes(1) = Max( start, tarray(1) )

!    tvis5dtmp = tvis5d ! v5dtimes(1)
!    print*, 'setting TIME to ',Max(0,start-tvis5d)
!    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start-tvis5d))   ) write(6,*) 'COMMAS:  Problem setting TIME'  
!
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5dtmp)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'  



!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

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
  CALL GET_VARIABLE(gd, 'IPCONC',   ipconc)

  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  
  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COMMAS:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COMMAS:  Problem setting NG'

! Allocate SOLVER scratch space
  
  m = (nx+2*ng)*(ny+2*ng)*(nz+2*ng)
!  allocate ( st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns) )
  
!  allocate ( ft(m,4) ) ! temporary arrays for solver

  tprt = start + tprint   - 1
  tstt = start + tstat    - 1
  tv5d = start 
  ntstart = start + dt

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
 wz     = GET_VARIABLE_INDEX(gd, 'WZ')
 elec   = GET_VARIABLE_INDEX(gd, 'EX')
 cion   = GET_VARIABLE_INDEX(gd, 'CPIONINIT')
 muz     = GET_VARIABLE_INDEX(gd, 'MUPOSZC')
 
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

 llen = index(microphys,' ')-1

  lfix = .false.
  IF ( microphys(llen:llen) == 'W' ) THEN
!    ipres = GET_VARIABLE_INDEX(gd, 'P')
!    itemp = GET_VARIABLE_INDEX(gd, 'T')
!    lfix = .true.
  
  ELSE
  
   ipres = 0
   itemp = 0
  
  ENDIF

!-----------------------------------------------------------------------------
! Main time step loop

   DO it = 1, nt

!----------------------------------------------------------------------
! Time and location of the grid
 
       time     = tarray(it)

       write(6,*) 'it loop: it,time1 = ',it,time



       IF ( time .gt. stop ) EXIT
       
       IF ( time .ge. start ) THEN

       write(6,*) 'it loop: it,time, netcdfversion = ',it,time,netcdfversion


    IF ( historyinput == 3 ) THEN
      write(stime,   '(i6.6)') time
      ncfilein = prefix(1:length)//'.'//stime//'.nc'
    ELSE
      ncfilein = prefix(1:length)//trim(number)//'.nc'
    ENDIF

!       CALL GRID_READ_NETCDF( gd, prefix(1:length)//trim(number)//'.nc', tarray(it), 1 )
       CALL GRID_READ_NETCDF( gd, ncfilein, tarray(it), 1 )

  CALL GET_VARIABLE(gd, 'NXEND',    nx)
  CALL GET_VARIABLE(gd, 'NYEND',    ny)
  CALL GET_VARIABLE(gd, 'NZEND',    nz)

  IF( .not. SET_VARIABLE(gd,'NX',       nx)      ) write(6,*) 'COM2NC:  Problem setting NX'
  IF( .not. SET_VARIABLE(gd,'NY',       ny)      ) write(6,*) 'COM2NC:  Problem setting NY'

!  IF( .not. SET_VARIABLE(gd,'NX',        nx)            ) write(6,*) 'COM2NC4:  Problem setting NX'
!  IF( .not. SET_VARIABLE(gd,'NY',        ny)            ) write(6,*) 'COM2NC4:  Problem setting NY'
!  IF( .not. SET_VARIABLE(gd,'NZ',        nz)            ) write(6,*) 'COM2NC4:  Problem setting NZ'
       
      
      IF ( lfix ) THEN
      
!      write(0,*) 'lqcevap,ld0,wz = ',lqcevap,ld0,wz
      

!$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,rhoair)
      DO k = 1,nzend-1
       DO j = 1,nyend-1
        DO i = 1,nxend-1
        
         IF ( lqcevap > 1 ) gd%var(wz+lqcevap)%flt3d(i,j,k) = -Abs(gd%var(wz+lqcevap)%flt3d(i,j,k))
         
         rhoair = 1.0e5*(gd%var(piinit)%flt1d(k)+gd%var(pi)%flt3d(i,j,k))**2.509/(287.04*gd%var(s)%flt3d(i,j,k))
         
         IF ( lnc > 1 ) gd%var(s+lnc-1)%flt3d(i,j,k) = gd%var(s+lnc-1)%flt3d(i,j,k)/rhoair
         IF ( lnr > 1 ) gd%var(s+lnr-1)%flt3d(i,j,k) = gd%var(s+lnr-1)%flt3d(i,j,k)/rhoair
         IF ( lni > 1 ) gd%var(s+lni-1)%flt3d(i,j,k) = gd%var(s+lni-1)%flt3d(i,j,k)/rhoair
         IF ( lns > 1 ) gd%var(s+lns-1)%flt3d(i,j,k) = gd%var(s+lns-1)%flt3d(i,j,k)/rhoair
         IF ( lnh > 1 ) gd%var(s+lnh-1)%flt3d(i,j,k) = gd%var(s+lnh-1)%flt3d(i,j,k)/rhoair
         IF ( lnhl > 1 ) gd%var(s+lnhl-1)%flt3d(i,j,k) = gd%var(s+lnhl-1)%flt3d(i,j,k)/rhoair
         
         IF ( ld0 > 1 ) gd%var(wz+ld0)%flt3d(i,j,k) = 1000.*gd%var(wz+ld0)%flt3d(i,j,k)
        
        ENDDO
       ENDDO
      ENDDO
      
      
      ENDIF ! lfix
      
      
    IF ( historyoutput == 3 ) THEN
      ncfileout = prefix(1:length)//'.'//stime//'.nc'//sversion
      CALL GRID_DEFINE_NETCDF(gd, ncfileout)
      CALL GRID_WRITE_NETCDF( gd, time, ncfileout)
    ELSE
      ncfileout = prefix(1:length)//trim(number1)//'.nc'//sversion
      CALL GRID_WRITE_NETCDF( gd, time, ncfileout)
    ENDIF

!       IF ( netcdfversion .eq. 4 ) THEN
!       CALL GRID_WRITE_NETCDF( gd, time, prefix(1:length)//trim(number1)//'.nc4')
!       ELSEIF ( netcdfversion .eq. 3 ) THEN
!       CALL GRID_WRITE_NETCDF( gd, time, prefix(1:length)//trim(number1)//'.nc3')
!       ENDIF
 


        ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO 
   
!  DEALLOCATE(st)
!  DEALLOCATE(ft)
 

!-----------------------------------------------------------------------------
! Print out date and time


 STOP
 END
