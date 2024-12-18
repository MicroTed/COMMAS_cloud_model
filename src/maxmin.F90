!-----------------------------------------------------------------------------
!
! COMMAS-1.0
!
! March 2006
!
!-----------------------------------------------------------------------------

 PROGRAM MAXMIN

! This program was set up to test reading of netcdf files and checking for bad
! values returned by the read subroutines. The root cause turned out to be a Lustre
! file system issue, combined with lack of error checking on the netcdf side.

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

  integer i,j,k, it, nt, count, in

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
  
  real :: smax, smin


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

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'MAXMIN:  PROBLEM READING NAMELIST'

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COM2V5D:  PROBLEM READING NAMELIST'

! Read in 10-ice parameters

  IF( .not. READ_NAMELIST(run_file(1:length),'ICE10_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF



  IF ( recalcdbz .ge. 1 ) THEN
    iusewetgraupel = iusewetgraupelv5d
    iusewethail = iusewethailv5d
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
    
    CALL GRID_READ_NETCDF( gd, ncfile, tarray(1), nx_or_nxend=.true. )

    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)

    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMMAS:  Problem setting TIME'  
!    IF( .not. SET_VARIABLE(gd,'TRESTART',  trestart)) write(6,*) 'COMMAS:  Problem setting TRESTART'
!    IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMMAS:  Problem setting THISTORY'
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'
    IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COMMAS:  Problem setting TPRINT'
    IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COMMAS:  Problem setting TSTAT'
    IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COMMAS:  Problem setting STOP'
    IF( .not. SET_VARIABLE(gd,'DT',        dt)      ) write(6,*) 'COMMAS:  Problem setting DT'
   

    time = start


!--------------------------------------------------------------------------------------
! Open VIS5D file  

!  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
    
    IF ( start .gt. tarray(nt) ) THEN
      write(0,*) 'Starting time is greater than the last record time! STOP!'
      STOP
    ENDIF
    
    start = Max(start, tarray(1) ) 
    
    DO it = nt,1,-1
!      IF ( start .le. tarray(it) ) v5dtimes(1) = tarray(it)
    ENDDO
    
!    v5dtimes(1) = Max( start, tarray(1) )

!    tvis5dtmp = tvis5d ! v5dtimes(1)
!    print*, 'setting TIME to ',Max(0,start-tvis5d)
!    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start-tvis5d))   ) write(6,*) 'COMMAS:  Problem setting TIME'  

!    print*, 'setting TIME to ',Max(0,start), ' tvis5d = ',tvis5d
    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start))   ) write(6,*) 'COMMAS:  Problem setting TIME'  

!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5dtmp)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'  


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
!  tv5d = start + tvis5d                  - 1

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
 IF ( tke_type == 1 ) THEN
 km     = GET_VARIABLE_INDEX(gd, 'KM')
 ELSE
! km     = GET_VARIABLE_INDEX(gd, 'KM')
 km     = GET_VARIABLE_INDEX(gd, 'TKE')
 ENDIF
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

!   IF ( tvis5d .ne. thistory ) iskip = Max( 1, tvis5d/thistory )
   print*, 'iskip = ',iskip

!-----------------------------------------------------------------------------
! Main time step loop

       IF ( historyoutput == 3 ) THEN 
         nt = (stop - start)/thistory + 1
         iskip = 1
       ENDIF

   DO it = 1, nt, iskip

!----------------------------------------------------------------------
! Time and location of the grid
 
       IF ( historyoutput == 3 ) THEN 
         time = start + (it - 1)*thistory
       ELSE
         time = tarray(it)
       ENDIF

       write(6,*) 'it loop: it,time1 = ',it,time


!       IF ( time .gt. stop ) EXIT
       
!       IF ( time .ge. start ) THEN

       write(6,*) 'it loop: it,time = ',it,time

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
         
         DO in = 1, 10
         
         smax = maxval( gd%var(s)%flt3d )
         smin = minval( gd%var(s)%flt3d(1:nxend-1,1:nyend-1,1:nzend-1) )

         write(0,*) 'TH max/min = ', smax, smin
         IF ( smin < 10.0 ) THEN
           write(0,*) 'BAD TH VALUE!!!!'
           
           DO k = 1,nzend-1
           count = 0
             DO j = 1,nyend-1
              DO i = 1,nxend-1
               IF ( gd%var(s)%flt3d(i,j,k) < 1.0 ) count = count+1
              ENDDO
             ENDDO
             
             
             IF ( count > 0 ) write(0,*) 'k, count = ',k,count
             
           ENDDO
           
           
         ENDIF
         
          IF (in < 10 ) CALL GRID_READ_NETCDF( gd, ncfile, tarray(it), 1 )
          
         ENDDO


!-----------------------------------------------------------------------------
! SUBTRACT OUT GRID MOTION
    IF ( ugrid .ne. 0.0 ) THEN
       IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(6,*) 'COMMAS:  Problem setting UGRID'  
    ENDIF
    IF ( vgrid .ne. 0.0 ) THEN
       IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(6,*) 'COMMAS:  Problem setting VGRID'
    ENDIF

        CALL FOLLOW(gd, -1.0, nxend)

        CALL GET_VARIABLE(gd, 'TIME',       itmp)
        
        write(6,*) 'time from netcdf is ',itmp

!    IF( .not. SET_VARIABLE(gd,'TIME',      time)   ) write(6,*) 'COMMAS:  Problem setting TIME'  

!        iv5dwritten(:) = 0
!        lv5dwrite = .true.
!        v5dtimes(itv5d) = time
 

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

 
!        DO k = 1,ns
!          st(1:nx,1:ny,1:nz,k) = gd%var(k+s-1)%flt3d(1:nx,1:ny,1:nz)
!        ENDDO

        
        write(0,*) 'call V5DOUT'

        
        
       
   

!        ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO 
   
  DEALLOCATE(ft)
 

!-----------------------------------------------------------------------------
! Print out date and time


 STOP
 END
