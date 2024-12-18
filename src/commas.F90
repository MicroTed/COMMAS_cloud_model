!-----------------------------------------------------------------------------
!
! COMMAS-7
!
! September 2008
!
!-----------------------------------------------------------------------------

 PROGRAM COMMAS

#ifdef NAG
  USE F90_UNIX_PROC
#endif

  USE GRID_MODULE
  USE INIT_MODULE, only: nzmeso, nbble
  USE GRIDIO_MODULE
#ifdef USE_PNETCDF
  USE GRIDPIO_MODULE, only : grid_write_pnetcdf
#endif
  USE FILE_MODULE
  USE CLINE_MODULE
  USE CPUTIME_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
  USE VIS5D_MODULE
  USE INDEX_MODULE
  USE FORCE_MODULE
  USE TRMM_MODULE
  USE TICK_MODULE
  USE TRAJ_MODULE
  USE TROPICAL_CYCLONE
  USE FOLLOW_MODULE
  USE BSS_NML, only: numbss

#ifdef MPI
#undef MPI
!  USE MPI
#define MPI 1
#endif

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

  character(LEN = 120)      :: run_file
  character,DIMENSION(120)  :: run_file_array      !needed for mpi process comm
  integer              :: status, length, length1
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
  integer :: rho
  integer :: km 
  integer :: s     
  integer :: vzf
  integer :: dbz   
  integer :: wz    
  integer :: elec
  integer :: xtra
  integer :: cion
  integer :: muz
   
  integer :: uinit  ! base state for current time
  integer :: vinit  
  integer :: uinit0  ! base state for t=0
  integer :: vinit0  
  integer :: winit  
  integer :: piinit 
  integer :: kminit,kmbaserm 
  integer :: sinit  

  integer :: gx     
  integer :: gy     
  integer :: gz     
  integer :: precip
  
  integer :: tcanp
  integer :: wcanp
  integer :: qav
  integer :: veg
  integer :: stype
  integer :: tsrfc
  integer :: wsfc
  integer :: tsoil
  integer :: wsoil
  integer :: vlai
  integer :: albedo
  integer :: rough
  integer :: qv
  integer :: eflx
  integer :: fflx
  integer :: uflx
  integer :: vflx
  integer :: tflx
  integer :: qflx
  integer :: radsw
  integer :: radlw

  real    :: tsurf, qsurf
  real    :: x_sw_loc_tmp, y_sw_loc_tmp
  integer :: itmp


!-----------------------------------------------------------------------------
! SOLVER SCRATCH MEMORY

  real, allocatable, target :: st(:,:)         ! These are big scratch arrays for solver
  real, allocatable :: ft(:,:)
  real, allocatable :: z1d4(:,:)
  real, allocatable :: precip_old(:,:,:)
!  real, pointer     :: ft(:,:)

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
  real    :: time_real
  integer, pointer :: ensize
  integer :: coards(6)
  integer :: trst       
  integer :: this
  integer :: thislast = -1
  integer :: tprt
  integer :: tstt
  integer :: tv5d
  integer :: i,j,k,ix,jy,kz
  integer :: iadd
  
  integer :: nxetmp, nyetmp
  integer :: tautorestart
  real    :: ugrid1, vgrid1, ugrid0, vgrid0
  
  logical :: file_exist = .false.
  logical :: file_exist_auto = .false.
  logical :: io_flag, lstt 
  logical :: clearflash = .false.
  integer :: rsttype = 1      ! 1 = normal restart; 2 = start dumping output to new files as rstprefix
  character(LEN=100)  :: filenm
  character(LEN=120) :: ncfile

  character(LEN=4) number00,number1,number2
  character(len=6) stime
  character(LEN=1) trjnum
  character(LEN=8) number
  character(LEN=8) number0
  character(LEN=6) timestr

  character(80) :: nclibver
  integer :: libver1,libver2,libver3,libver
  
  real :: uold, vold
  real :: xold = 0.0, yold = 0.0
  real :: xbar, ybar

  real :: tim1,tim2
!  double precision :: dt1, dt2, delt, dtSTRZSORmpi
  real :: commas_dtime2

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

   logical :: mpisep = .false.

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
#ifdef MPI
   integer :: mpi_status(MPI_Status_size)
#endif
   integer :: token

      double precision, allocatable :: mpitotinth(:,:),mpitotoutth(:,:)

!  integer :: mpi_error_code

   logical :: debug_mpi = .false.

!-----------------------------------------------------------------------------
! MPI STARTUP

  CALL COMMASMPI_STARTUP()

   IF ( .false. ) THEN
   write(0,*) 'COMMAS:  STARTUP COMPLETED ON rank ', my_rank
   ENDIF

!-----------------------------------------------------------------------------

  CALL cld_cpu('-0')  

  CALL cld_cpu('MAIN')  

!-----------------------------------------------------------------------------
! Need 3D.RUN file in order to get run parameters - get from the command line

  DO i = 1, LEN(run_file)
    run_file(i:i) = " "
  END DO

  IF (my_rank == 0) THEN
    IF( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN

      write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  NEED 3D.RUN FILENAME'
      write(0,*) 'COMMAS:  INCORRECT ARGUMENTS ON COMMAND LINE:  EXITING RUN!...'
      call commasmpi_abort()
      call exit(1)

    ELSE

     CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)

!     write(6,*) 'COMMAS:  Input 3D.RUN file is:  ', run_file(1:length)

     IF( status .ne. 0 ) THEN

      write(0,*) 'COMMAS:  COULD NOT RETRIEVE 1ST COMMAND LINE ARG:  EXITING RUN!...', status
      write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ!'
      call commasmpi_abort()
      call exit(1)

     ENDIF

    ENDIF

    DO i = 1, length
      run_file_array(i) = run_file(i:i)
    END DO

  END IF !! (my_rank == 0)

#ifdef MPI
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

    IF ( my_rank == 0 ) THEN
      write(0,*) 'COMMAS:  INPUT 3D.RUN FILE:  ', run_file(1:length), ' DOES NOT EXIST!!!'
      write(0,*) 'COMMAS:  DOES NOT HAVE A 3D.RUN FILE TO READ, EXITING'
    ENDIF
    CALL COMMASMPI_ABORT()
    stop

  ENDIF

! READ IN RUN MODEL NAMELIST

 IF( .not. READ_NAMELIST(run_file(1:length),'RUN')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
 ENDIF

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF
! READ IN GRID PARAMETERS

  IF( .not. READ_NAMELIST(run_file(1:length),'GRIDN') ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST GRIDN'
  ENDIF

  x_sw_loc_tmp = x_sw_loc
  y_sw_loc_tmp = y_sw_loc
 
! READ IN MPI PARAMETERS
#ifdef MPI
  IF( .not. READ_NAMELIST(run_file(1:length),'MPI_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF
#endif
! READ IN HOMOGENEOUS PARAMETERS

  IF( .not. READ_NAMELIST(run_file(1:length),'HOMOG_INIT')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF

  IF( .not. READ_NAMELIST(run_file(1:length),'FORCING')  ) THEN
    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST FORCING'
  ENDIF

!  IF( (.not. READ_NAMELIST(run_file(1:length),'BALLOON'))  ) THEN
!    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST BALLOON'
!  ENDIF

  IF( (.not. READ_NAMELIST(run_file(1:length),'TRAJECTORIES'))  ) THEN
    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST TRAJECTORIES'
  ENDIF

  IF( (.not. READ_NAMELIST(run_file(1:length),'OUTPUT_OPTIONS'))  ) THEN
    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST OUTPUT_OPTIONS'
  ENDIF


! READ IN OTHER PARAMETERS
! Read in lfo parameters

  IF( .not. READ_NAMELIST(run_file(1:length),'LFO_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST'
  ENDIF

! Read in 10-ice parameters

  IF( .not. READ_NAMELIST(run_file(1:length),'ICE10_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST ICE10_PARAMS'
  ENDIF
  
! Read in Milbrandt-Yau parameters

  IF ( microphys(1:2) == 'MY' ) THEN
  IF( .not. READ_NAMELIST(run_file(1:length),'MY_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST MY_PARAMS'
  ENDIF
  ENDIF

  microp = microphys
    IF ( microphys(1:4) .eq. 'ZVDM' .or. microphys(1:5) .eq. 'ZVDHM' .or. &
         microphys(1:5) .eq. 'ZIEGM'  .or. microphys(1:6) .eq. 'ZIEGHM' ) mixedphase = .true.

  
  IF ( RKSCHEME == 2 ) THEN
  ! For RK2, if using 6th order then change to 5th order
    IF ( atypes2 == 5 ) atypes2 = 1
  ENDIF

! print*, 'COMMAS: READ PARAMH'

  IF( .not. READ_NAMELIST(run_file(1:length),'PARAMH')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST PARAMH'
  ENDIF
  
  ! If TKE scheme is sqrtE, then force mix_type = 1
  IF ( tke_type > 1 .and. mix_type /= 1 ) THEN
   IF ( my_rank == 0 ) print *, 'COMMAS: WARNING: tke_type > 1 but mix_type /= 1, so setting mix_type=1'
   mix_type = 1
  ENDIF


  IF( .not. READ_NAMELIST(run_file(1:length),'BSS_INIT')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST BSS_INIT'
  ENDIF

! check for 7th-order scheme in MPI
!  IF (( atypes2 == 10 .or. atypem2 == 10 .or. atypes1 == 10 .or. atypem1 == 10 ) .and. number_of_processes > 1) THEN
  IF ( atypes2 == 10 .or. atypem2 == 10 .or. atypes1 == 10 .or. atypem1 == 10 ) THEN
    ng = 4
    ngv = 4
  ENDIF
  IF ( atypes2 == 16 .or. atypem2 == 16 .or. atypes1 == 16 .or. atypem1 == 16 .or.  &
       atypes2 == 17 .or. atypem2 == 17 .or. atypes1 == 17 .or. atypem1 == 17 ) THEN
    ng = 5
    ngv = 5
  ENDIF

   IF ( .not. parallelio ) THEN
     parallelio_type = 1 ! do not set to 2 or else IO gets stuck in a barrier
     parallelio_in = 0
   ELSE
     IF ( parallelio_in == -1 ) THEN
       parallelio_in = 1 ! assume restart is also parallel io if parallelio_in has default value
     ENDIF
     IF ( parallelio_type == 2 .or. parallelio_type == 1 ) THEN
       
       IF ( parallelio_type == 2 ) netcdfversion = 3
#if defined(MPI) & defined(USE_PNETCDF)       
  nclibver = nf90_inq_libvers()
  read(nclibver(1:1),'(i1)') libver1
  read(nclibver(3:3),'(i1)') libver2
  read(nclibver(5:5),'(i1)') libver3
  
  libver = 100*libver1 + 10*libver2 + libver3
  IF ( parallelio_type == 2 ) THEN
  IF ( libver >= 460 ) THEN
     pnetcdf_flag = nf90_64bit_offset
  ELSE
     pnetcdf_flag = NF90_PNETCDF
  ENDIF
  
  IF ( libver > 432 .and. libver < 460 ) THEN
    IF ( my_rank == 0 ) THEN
      write(0,*) 'WARNING! Netcdf version '//nclibver(1:5)//' may not set file type properly for pnetcdf.'
      write(0,*) 'Recommend version 4.3.1.1 or 4.6.0 (or higher)'
      write(0,*) 'Version 4.3.3.1 also works with a fix to line 75 of libsrc5/nc5dispatch.c'
    ENDIF
  ENDIF
  
  ELSEIF ( parallelio_type == 1 ) THEN
  
  IF ( my_rank == 0 ) THEN
    write(0,*) 'Netcdf version '//nclibver(1:5)
    write(0,*) 'parallelio, deflate_level = ', parallelio, deflate_level
  ENDIF
  IF ( libver >= 474 .and. parallelio .and. deflate_level >= 1 ) THEN ! allow compression for 4.7.4
    IF ( my_rank == 0 ) write(0,*) 'Netcdf version '//nclibver(1:5)//' supports compression with parallel'
    parallel_compress = .true.
  ENDIF
  
  ENDIF

#endif
     ENDIF
   ENDIF

   CALL STRING_LIMITS(prefix, ibeg, iend)
   length = iend - ibeg + 1
   IF ( my_rank == 0 ) THEN
    OPEN(unit=105, file=prefix(1:length)//'.namelist1', status='unknown',form='formatted',delim='APOSTROPHE')

    CALL WRITE_NAMELISTS(105)
   
    CLOSE(unit=105)
   ENDIF

  ibsd = 1; iesd = nx-1 ; jbsd = 1 ; jesd = ny-1
  ncxe = nx
  ncye = ny
  ncze = nz

!
!.... CLZ (3/15/13): set nzmeso = nz - 1 here for subsequent use
!

  nzmeso = nz - 1

!  write(luno,*)
!  write (luno,*) 'set nzmeso=',nzmeso,' (nzmeso = nz - 1)'
!  write(luno,*)

!
!.... CLZ (3/26/13): make sure khomog is set
!

      if (khomog .ge. 1 .and. khomog .le. nzmeso) then
!      write(luno,*)'khomog has been input =',khomog

      elseif (khomog .lt. 1 .or. khomog .gt. nzmeso) then
      khomog = nzmeso
!      write(luno,*)'khomog has been set in code =',khomog,' (= nzmeso)'

      endif

  IF ( historyinput .eq. -1 ) historyinput = historyoutput

!-----------------------------------------------------------------------------
! CHECK FOR ADDITIONAL COMMAND LINE ARGUMENTS

#ifndef MPI

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
#endif /* end of non-mpi section */

   write(stime,   '(i6.6)') 0

   write(number1, '(a,i3.3)') '.', member
#ifdef MPI
   write(number2, '(a,i3.3)') '.', my_rank
   write(number00, '(a,i3.3)') '.', 0

   number = number1//number2
   number0 = number1//number00
#else
   number = number1
   number0 = number1
#endif

  rsttype = 1
  IF ( rstprefix(1:1) .ne. ' ' ) THEN
    rsttype = 2
    IF ( rstprefix .eq. prefix ) THEN
      write(0,*) 'ERROR: rstprefix = prefix = ',rstprefix
      write(0,*) 'STOPPING'
      CALL COMMASMPI_ABORT()
      STOP
    ENDIF
  ENDIF

   IF ( rsttype .eq. 1 ) THEN
     CALL STRING_LIMITS(prefix, ibeg, iend)
     v5dfilename = prefix(ibeg:iend)//trim(number)//'.v5d'
     onedfilename = prefix(ibeg:iend)//trim(number)//'.1d'
   ELSEIF ( rsttype .eq. 2 ) THEN
     CALL STRING_LIMITS(rstprefix, ibeg, iend)
     v5dfilename = rstprefix(ibeg:iend)//trim(number)//'.v5d'
     onedfilename = rstprefix(ibeg:iend)//trim(number)//'.1d'
   ENDIF
   

!-----------------------------------------------------------------------------
!  MPI SET UP ( allocate / bmpi / empi ) (non-MPI compile runs alternate subroutines)

  CALL COMMASMPI_BOUNDS(nx,ny,nz)
  CALL COMMASMPI_LINKS(bcx,bcy)
!  CALL COMMASMPI_ALLOCATE(ng)

  IF ( isfcl .eq. -1 ) THEN
!     bcz = -1
  ENDIF

#ifdef MPI
  ibsd = 1; iesd = nxend-1 ; jbsd = 1 ; jesd = nyend-1
!  IF ( ixend == nxend ) iesd = ixend - ixbeg
!  IF ( jyend == nyend ) jesd = jyend - jybeg
#endif
!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   IF ( nyend > 2 ) onedoutput = 0
   IF ( onedoutput > nxend ) THEN
      write(0,*) 'WARNING: onedoutput > nxend! Setting onedoutput = nxend/2'
      onedoutput = nxend/2
   ENDIF

  5 CONTINUE

  IF( start .LT. 0 ) THEN !{

! Compute some needed variables

  CALL cld_cpu('I/O')
  
   CALL STRING_LIMITS(prefix, ibeg, iend)
   length = iend - ibeg + 1

! check for existence of main .out file
   INQUIRE(file=prefix(1:length)//trim(number0)//'.out', exist=file_exist)

#ifdef MPI
! Wait for all processes to get here. Otherwise it can happen that 
! rank 0 creates the .out file before all processes have tested for the file, which can result in their aborting.
   CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
   IF ( file_exist ) THEN
     
! check if a file exists with an automatic restart time. This is used for systems where a 
! job might get preempted and then restarted
     INQUIRE(file=prefix(1:length)//'.autorestart', exist=file_exist_auto)
     
     IF ( file_exist_auto .and. member == 0 ) THEN
!       IF ( my_rank == 0 ) THEN
       open(45,file=prefix(1:length)//'.autorestart',status='old',form='formatted')
       read(45,*) tautorestart
       close(45)
!       ENDIF
#ifdef MPI
! could do a broadcast here if only rank 0 reads the time in the file
#endif
       
       IF ( tautorestart < stop .and. tautorestart > 0 ) THEN
         start = tautorestart
         CALL cld_cpu('I/O')
         GOTO 5
       ENDIF
     ENDIF
     
     IF ( my_rank == 0 .or. verbose_mpi ) THEN
      write(0,*) 'COMMAS: NEW OUTPUT FILE ',prefix(1:length)//trim(number)//'.out',' ALREADY EXISTS!',my_rank
      write(0,*) 'PLEASE DELETE FIRST IF YOU REALLY WANT TO OVERWRITE!'
      write(0,*) 'STOP, rsttype = ',rsttype
     ENDIF
     IF ( rsttype .eq. 2 ) write(0,*) 'REMEMBER TO SET START > 0!'
     CALL COMMASMPI_ABORT()
     STOP
   ELSE
     IF ( my_rank == 0 .or. verbose_mpi ) write(6,*) 'COMMAS:  START < 0, INITIALIZING GRID FOR SIMULATION ', &
              prefix(1:length)//trim(number)
   ENDIF

       
! create checkpoint time file
       IF ( my_rank == 0 .and. member == 0 ) THEN
       open(45,file=prefix(1:length)//'.autorestart',status='unknown',form='formatted')
       write(45,*) 0
       close(45)
       ENDIF

#ifdef MPI
   dx        = xdomain / float(nxend-1)
   dy        = ydomain / float(nyend-1)
   dz        = zdomain / float(nzend-1)
#else
   dx        = xdomain / float(nx-1)
   dy        = ydomain / float(ny-1)
   dz        = zdomain / float(nz-1)
#endif

! Check for 2-D mode and set same grid spacing
   IF ( ny .le. 2 ) THEN 
     dy = dx
     ycntr(:) = 0.5*dy
     ywfcen(:) = 0.5*dy
     yssfcen(:) = 0.5*dy
     topforceY = 0.5*dy
   ENDIF
   IF ( nx .le. 2 ) THEN
     dx = dy
     xcntr(:) = 0.5*dx
     xwfcen(:) = 0.5*dx
   ENDIF

   coards(1) = year
   coards(2) = month
   coards(3) = day
   coards(4) = hour
   coards(5) = minute
   coards(6) = second

   CALL GRID_DEFINE_FROM_LIST(gd, microphys, ipconc, ichaff, inucopt, isfcphys=isfcphys)

   IF ( microphys(1:5) .eq. 'ICE10' .or. microphys(1:1) .eq. 'Z' .or. microphys(1:5) .eq. 'WARMZ' &
        .or. microphys(1:3) .eq. 'HCM' .or. microphys(1:3) .eq. 'TAK' .or. microphys(1:4) .eq. 'MORR' &
        .or. microphys(1:4) .eq. 'THOM' ) THEN

   IF ( my_rank == 0 ) THEN
   OPEN(unit=3, file=prefix(1:length)//trim(number)//'.stat', status='unknown',form='formatted')
!     write(3,'(a,3(i6,1x),i3)') 'nx,ny,nz = ',nx,ny,nz
   
    IF ( onedoutput > 0 .and. nyend <= 2 ) OPEN(unit=33, file=onedfilename, status='unknown',form='formatted')
   
   ENDIF
   
   ENDIF
   IF ( ipelec .ge. 2 .and. my_rank == 0 ) THEN
    OPEN(unit=16, file=prefix(1:length)//trim(number)//'.flash', status='unknown',form='formatted')
!    OPEN(unit=17, file=prefix(1:length)//trim(number)//'.flash2', status='unknown',form='formatted')
   ENDIF

   CALL GRID_SET_ATTRIBUTES(gd,prefix(1:length),ne,member,coards  &
                           ,nx,ny,nz,microphys,v5dflds           &
#ifdef MPI                           
                           ,nxend,nyend,nzend,my_rank)
#else
!( leave for paren balance
                           )
#endif


   CALL GRID_ALLOCATE(gd, nx, ny, nz)

!  CALL GRID_ALLOCATE(gd, runname, member, coards, nx, ny, nz, microphys)

!   lun = FILE_OPEN(gd, 'OUTPUTFILE', 'BEGINNING MODEL RUN')                     ! Open I/O file for model output
   IF ( my_rank == 0 .or. verbose_mpi ) THEN
   
    lun = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'BEGINNING MODEL RUN',.true.)                     ! Open I/O file for model output

    CALL FILE_MESSAGE('COMMAS:  VARIABLE NAMESPACE DEFINED')
    
    luno = lun

   ENDIF

    lune = 0
   
   IF ( historyoutput .ne. 2 .and. historyoutput >= -1 ) THEN
#ifdef MPI
!   IF ( number_of_processes .gt. 1 )  netcdfversion = 3
  IF ( .not. parallelio ) THEN
   IF (  historyoutput .eq. 1 ) THEN
     ncfile = prefix(1:length)//trim(number1)//'.nc'
     IF ( my_rank == 0 ) CALL GRID_DEFINE_NETCDF(gd,ncfile)
   ELSEIF ( historyoutput .eq. 0 .or. historyoutput .eq. -1 ) THEN
     ncfile = prefix(1:length)//trim(number)//'.nc'
     CALL GRID_DEFINE_NETCDF(gd,ncfile)
   ELSEIF ( historyoutput .eq. 3  ) THEN
     ncfile = prefix(1:length)//'.'//stime//'.nc'
     IF ( my_rank == 0 ) CALL GRID_DEFINE_NETCDF(gd,ncfile)
   ENDIF
  ELSE
   IF ( historyoutput /= 3 ) THEN
     ncfile = prefix(1:length)//trim(number1)//'.nc'
   ELSEIF ( historyoutput .eq. 3  ) THEN
     ncfile = prefix(1:length)//'.'//stime//'.nc'
   ENDIF

     CALL GRID_DEFINE_NETCDF(gd,ncfile)
! Wait for all processes to get here.  
    CALL MPI_BARRIER(my_comm, mpi_error_code)


  ENDIF
! #elif defined (MPI)
!   CALL GRID_DEFINE_NETCDF(gd,prefix(1:length)//trim(number)//'.nc')
#else
   IF ( historyoutput .ne. 3 .and. historyoutput >= -1 ) THEN
    CALL GRID_DEFINE_NETCDF(gd)
   ELSE
     ncfile = prefix(1:length)//'.'//stime//'.nc'
     CALL GRID_DEFINE_NETCDF(gd,ncfile)
   ENDIF
#endif

   IF ( my_rank == 0 .or. verbose_mpi )  CALL FILE_MESSAGE('COMMAS:  NETCDF FILE DEFINED')


   ENDIF
   
   IF ( my_rank == 0 .or. verbose_mpi ) lun = FILE_CLOSE()

   IF ( my_rank == 0 ) THEN
    OPEN(unit=105, file=prefix(1:length)//'.namelist', status='unknown',form='formatted',delim='APOSTROPHE')

    CALL WRITE_NAMELISTS(105)
   
    CLOSE(unit=105)
   ENDIF

! 
!
!-----------------------------------------------------------------------------
! INIT GRIDS

   CALL INIT_GRID(gd) 

   CALL PRTINFO( gd )

   CALL INIT_BACKGROUND_1D(gd,.true.) 

! prtbase functions are now done in INIT_BACKGROUND_1D
!   CALL PRTBASE( gd )


!
!-----------------------------------
!



!
!.... initial southwest origin point of model grid in fixed mesoscale domain at surface
!

  CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
  CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)



!
!.... CLZ (11/10/11): added inhom (= 1, 2, 3, 4) options for heterogeneous initialization
!


!
!.... CLZ (3-6-13): mesoscale domain arrays allocated in INIT_*_DEF_*, NOT later!!!
!

   if(inhom .eq. 1) then  ! stationary grid initialized from (stationary) Lagrangian analysis
   CALL INIT_USER_DEF(gd)

   elseif(inhom .eq. 2) then  ! stationary grid initialized from parametric functions
   CALL INIT_IDEAL_DEF(gd)

   elseif(inhom .eq. 3) then  ! moving grid initialized from parametric idealized functions
   CALL INIT_IDEAL_DEF_MOV(gd,x_sw_loc,y_sw_loc)      

   elseif(inhom .eq. 4) then  ! moving grid initialized from sounding-constrained functions
   CALL INIT_DATA3D_DEF_MOV(gd,x_sw_loc,y_sw_loc)      

   elseif(inhom .eq. 5) then ! Stationary grid initialized from 2 soundings to create EW front
!   write(0,*) 'About to call RLM subroutine'
   CALL INIT_DATA3D_DEF_FRONT(gd, x_sw_loc, y_sw_loc)
!   write(0,*) 'Done with RLM subroutine'

   endif  ! endif(inhom .eq. 1) then



   time = 0

   ugrid0 = ugrid
   vgrid0 = vgrid
   IF( .not. SET_VARIABLE(gd,'TIME',      time)    ) write(6,*) 'COMMAS:  Problem setting TIME'
   IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(0,*) 'COMMAS:  Problem setting UGRID'  
   IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(0,*) 'COMMAS:  Problem setting VGRID'
   IF( .not. SET_VARIABLE(gd,'UGRID0',     ugrid)   ) write(0,*) 'COMMAS:  Problem setting UGRID'  
   IF( .not. SET_VARIABLE(gd,'VGRID0',     vgrid)   ) write(0,*) 'COMMAS:  Problem setting VGRID'
   IF( .not. SET_VARIABLE(gd,'TRESTART',trestart)  ) write(6,*) 'COMMAS:  Problem setting TRESTART'
   IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMMAS:  Problem setting THISTORY'
   IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'
   IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COMMAS:  Problem setting TPRINT'
   IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COMMAS:  Problem setting TSTAT'
   IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COMMAS:  Problem setting STOP'
     xbar = 0.0
     ybar = 0.0
   IF( .not. SET_VARIABLE(gd,'XBAR',      xbar)     ) write(6,*) 'COMMAS:  Problem setting XBAR'  
   IF( .not. SET_VARIABLE(gd,'YBAR',      ybar)     ) write(6,*) 'COMMAS:  Problem setting YBAR'  

   IF ( bogusvortex == 1 ) THEN
     CALL HURRFORCE(gd, 0.0, nx, ny, nz, ng)
   ELSEIF ( bogusvortex == 2 .or. bogusvortex == 3 ) THEN
     CALL INIT_TC(gd, 0.0, nx, ny, nz, ng)
   ENDIF
!
!.... CLZ (11/10/11): restored call of INIT_PERT_BBLE as in original COMMAS code
!

   CALL INIT_PERT_BBLE(gd)


!
!.... CLZ (3/16/13): if horizontally inhomogeneous initiation, reimpose hydrostatic balance
!
   
   IF (inhom .gt. 0 ) THEN
     CALL HYDRO_BALANCE(gd)   ! hydrostatically balance heterogeneous initial state
   ENDIF



 IF ( isfcphys >= 1 ) THEN
 

 tcanp   = GET_VARIABLE_INDEX(gd, 'TCANP')
 wcanp   = GET_VARIABLE_INDEX(gd, 'WCANP')
 qav     = GET_VARIABLE_INDEX(gd, 'QAV')
 veg     = GET_VARIABLE_INDEX(gd, 'VEG')
 stype   = GET_VARIABLE_INDEX(gd, 'STYPE')
 tsrfc   = GET_VARIABLE_INDEX(gd, 'TSRFC')
 wsfc    = GET_VARIABLE_INDEX(gd, 'WSFC')
 tsoil   = GET_VARIABLE_INDEX(gd, 'TSOIL')
 wsoil   = GET_VARIABLE_INDEX(gd, 'WSOIL')
 vlai    = GET_VARIABLE_INDEX(gd, 'VLAI')
 albedo  = GET_VARIABLE_INDEX(gd, 'ALBEDO')
 rough   = GET_VARIABLE_INDEX(gd, 'ROUGH')
 qv      = GET_VARIABLE_INDEX(gd, 'QV')

 eflx    = GET_VARIABLE_INDEX(gd, 'EFLX')
 fflx    = GET_VARIABLE_INDEX(gd, 'FFLX')
 uflx    = GET_VARIABLE_INDEX(gd, 'UFLX')
 vflx    = GET_VARIABLE_INDEX(gd, 'VFLX')
 tflx    = GET_VARIABLE_INDEX(gd, 'TFLX')
 qflx    = GET_VARIABLE_INDEX(gd, 'QFLX')  
 radsw   = GET_VARIABLE_INDEX(gd, 'RADSW')  
 radlw   = GET_VARIABLE_INDEX(gd, 'RADLW') 

 pi     = GET_VARIABLE_INDEX(gd, 'PI')
 s      = GET_VARIABLE_INDEX(gd, 'TH')
 piinit = GET_VARIABLE_INDEX(gd, 'PIINIT')
 
 CALL GET_VARIABLE(gd,'TSFC',tsurf) 
 CALL GET_VARIABLE(gd,'QSFC',qsurf) 
  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1

      call set_var_sfcphy(gd%var(tcanp),gd%var(wcanp),  &
           gd%var(qav),gd%var(veg),gd%var(stype),       &
           gd%var(tsrfc),gd%var(wsfc),gd%var(tsoil),    &
           gd%var(wsoil),gd%var(vlai),gd%var(albedo),   &
           gd%var(rough),nx,ny,nz,ns,gd%var(s),         &
           gd%var(eflx),gd%var(fflx),gd%var(uflx),      &
           gd%var(vflx),gd%var(tflx),gd%var(qflx),      &
           gd%var(radsw),gd%var(radlw),gd%var(piinit),  &
           gd%var(pi),tsurf,qsurf)
 
 ENDIF

#ifdef MPI
   CALL PRINT_MPI( gd )
#else
   CALL PRINT( gd )
#endif

   IF ( historyoutput .eq. 1 .or. historyoutput .eq. 3 ) THEN ! writing a single netcdf file, either parallel or round-robin.



#ifdef MPI
   IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
   token = 1

   if (my_rank==0) then
!     CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)    !head node will write first
     IF ( number_of_processes .gt. 1 ) THEN
     CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
     ENDIF
   endif

   CALL MPI_BARRIER(my_comm, mpi_error_code)

   IF ( number_of_processes .gt. 1 ) THEN
   do i = 1,number_of_processes-1
    if (my_rank==i) then
     CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
!     CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')          !opens/writes to nc file
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)          !opens/writes to nc file
     if (my_rank.ne.number_of_processes-1) THEN
      CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
     ENDIF
    endif
   enddo

   CALL MPI_BARRIER(my_comm, mpi_error_code)

   ENDIF
   
   ELSE  ! parallelio
     IF ( parallelio_type == 1 ) THEN
#ifdef NC4
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
     IF ( my_rank==0 ) write(0,*) 'Cannot use parallelio unless netCDF4 support is enabled!'
     call commasmpi_abort()
#endif
     ELSEIF ( parallelio_type == 2 ) THEN

#ifdef USE_PNETCDF
!     write(0,*) 'COMMAS PnetCDF file: rank, ncfile = ',my_rank,ncfile
     CALL GRID_WRITE_PNETCDF( gd, 0, ncfile) 
#else

#if defined( NC4 )
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
! #if !defined( NC4 )
     IF ( my_rank==0 ) write(0,*) 'Cannot use pnetcdf unless netCDF4 or pnetcdf support is enabled!'
     call commasmpi_abort()
#endif

#endif
     
     ENDIF
   ENDIF

  ELSEIF ( historyoutput .eq. 0 .or. historyoutput .eq. -1 ) THEN
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)    ! each node writes to its own file

     CALL MPI_BARRIER(my_comm, mpi_error_code)
   
#else
! Non-mpi
  IF ( historyoutput /= 3 .and. historyoutput >= -1) THEN
    CALL GRID_WRITE_NETCDF( gd, 0 )
  ELSE
    CALL GRID_WRITE_NETCDF( gd, 0, ncfile)
  ENDIF
#endif

    ELSEIF ( historyoutput .eq. 2 ) THEN ! fortran binary output
    
! Write a separate binary file for each tile -- faster than round-robin netcdf
    write(timestr,timfmt) 0
    CALL GRID_WRITE_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )

#ifdef MPI
    CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif

   ENDIF
   
   IF ( my_rank==0 ) lun = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'COMMAS: INIT COMPLETED AND GRID WRITTEN TO A FILE')
  
   start = 0
   
!  Initialize vis5d output if called for

!  IF ( tvis5d .lt. stop .or. ( tvis5dstart .ne. 0 .and. tvis5dstart .lt. stop ) ) THEN
!    n = V5DINIT()
!  ENDIF
     gx     = GET_VARIABLE_INDEX(gd, 'XC')
     gy     = GET_VARIABLE_INDEX(gd, 'YC')
     gz     = GET_VARIABLE_INDEX(gd, 'ZC')
     IF ( stop .ge. tvis5d ) THEN
#ifdef MPI
     mpisep = .true.
#endif
     call Vis5DInit( gd, gd%var(gx), gd%var(gy), gd%var(gz), v5dstridex, v5dstridez, mpisep )

     ENDIF
   IF ( my_rank==0 ) lun = FILE_CLOSE()

  CALL cld_cpu('I/O')


#ifdef MPI
!   call simple_xy_wr(my_rank, number_of_processes)
    IF ( parallelio .and. parallelio_type == 1 .and. historyoutput == 3 ) THEN
!      call GRID_DEFWRITE_NETCDFPAR(gd, 'test.nc')
    ENDIF
#endif


  ELSE  !} model is restarting from a specific time

!-----------------------------------------------------------------------------
! ELSE, READ IN THE DATA !{

  nxetmp = nxend
  nyetmp = nyend
  
  CALL cld_cpu('I/O')

    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

    IF ( my_rank == 0 ) THEN
    write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)
    ENDIF

! grid_info
   itmp = 0
  10 CONTINUE ! can come back to this point to re-read the data if the inhom=4 to reinitialize the inhomogeneous base state

   IF ( historyinput .eq. 1 .or. historyinput .eq. 3) THEN

     IF ( historyinput == 1 ) THEN
       ncfile = prefix(1:length)//trim(number1)//'.nc'
     ELSE
       write(stime,'(i6.6)') start
       ncfile = prefix(1:length)//'.'//stime//'.nc'
       IF ( my_rank == 0 ) write(0,*) 'reading in from ncfile = ',ncfile
     ENDIF

#ifdef MPI
   IF ( .not. parallelio .or. parallelio_in /= 1 ) THEN
!token ring to only allow one process to access nc file at a time.
! BUT don't need to worry about reading -- all processes can read the file simultaneously
   token = 1

   
!   if (my_rank==0) then
     IF ( itmp == 0 ) THEN
     CALL GRID_READ_NETCDF( gd, ncfile, start,microphys_namelist=microp)    !head node will read first
     ELSE
     CALL GRID_READ_NETCDF( gd, ncfile, start, no_alloc=1,microphys_namelist=microp)    !head node will read first, grid is allready allocated
     ENDIF
!     IF ( number_of_processes .gt. 1 ) THEN
!      CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)         !pass token to next proc
!     ENDIF
!   endif

   IF ( .false. .and. number_of_processes .gt. 1 ) THEN
   do i = 1,number_of_processes-1
    if (my_rank==i) then
     CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
     IF ( itmp == 0 ) THEN
     CALL GRID_READ_NETCDF( gd, ncfile, start)       !opens/reads to nc file
     ELSE
     CALL GRID_READ_NETCDF( gd, ncfile, start, no_alloc=1)       !opens/reads to nc file, grid is allready allocated
     ENDIF
     if (my_rank.ne.number_of_processes-1) CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
    endif
   enddo
   ENDIF
   
   CALL MPI_BARRIER(my_comm, mpi_error_code)
   
   ELSE ! parallelio 
     CALL GRID_READ_NETCDF( gd, ncfile, start)
   ENDIF
   
  ELSEIF ( historyinput .eq. 0 ) THEN
     CALL GRID_READ_NETCDF( gd, prefix(1:length)//trim(number)//'.nc', start) 

     CALL MPI_BARRIER(my_comm, mpi_error_code)
   
#else
    IF ( .not. restart_separate .or. ne <= 1 ) THEN
      CALL GRID_READ_NETCDF( gd, ncfile, start,microphys_namelist=microp)
    ELSE
      CALL GRID_READ_NETCDF( gd, prefix(1:length)//'.r'//trim(number)//'.nc', start,microphys_namelist=microp)
    ENDIF
#endif

    ELSEIF ( historyinput .eq. 2 ) THEN ! fortran binary output
    
! Read a separate binary file for each tile -- faster than round-robin netcdf
    CALL GRID_DEFINE_FROM_LIST(gd, microphys,ipconc,ichaff,inucopt, isfcphys=isfcphys)
    CALL GRID_ALLOCATE(gd, nx, ny, nz)
    IF ( ne .gt. 1 ) THEN
      iadd = 1
      write(timestr,timfmt) start+iadd
      INQUIRE(file=prefix(1:length)//trim(number)//'.'//timestr//'.bin', exist=file_exist)
      
      IF ( .not. file_exist ) iadd = 0
      
    ELSE
      iadd = 0
    ENDIF
    write(timestr,timfmt) start+iadd
    write(6,*) 'start, timestr = ',start,timestr
    CALL GRID_READ_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )

   ENDIF
   
   IF (inhom .eq. 4 .and. itmp == 0) THEN  ! moving grid initialized from sounding-constrained functions
   ! this is a hack for restarting when using the inhomogeneous environment. The call to INIT_DATA3D_DEF_MOV
   ! will overwrite some 3d arrays, so have to go back (goto -- I hate them) and read in the data again
    itmp = 1
    CALL INIT_DATA3D_DEF_MOV(gd,x_sw_loc_tmp,y_sw_loc_tmp)
    GOTO 10
   
   ENDIF  ! (inhom .eq. 4)
   
   
#ifdef MPI
  CALL GET_VARIABLE(gd, 'NXEND',    nxend)
  CALL GET_VARIABLE(gd, 'NYEND',    nyend)
    IF ( nxend /= nxetmp .or. nyend /= nyetmp ) THEN
      IF ( my_rank == 0 ) THEN
        write(0,*) 'Problem with nxend or nyend in file does not match input!'
        write(0,*) 'nx,ny,nxend,nyend = ',nxetmp,nyetmp,nxend,nyend
      ENDIF
      CALL commasmpi_abort()
      STOP
    ENDIF
#endif



    IF ( rsttype .eq. 2 ) THEN ! write out a new output file
      CALL STRING_LIMITS(rstprefix, ibeg, iend)
      length = iend - ibeg + 1
      IF ( my_rank == 0 ) THEN
      write(6,*) 'COMMAS:  START = ',start,' CREATING NEW OUTPUT FILE ',rstprefix(1:length)//trim(number)
      ENDIF

   INQUIRE(file=rstprefix(1:length)//trim(number)//'.out', exist=file_exist)
#ifdef MPI
! Wait for all processes to get here. Otherwise it can happen that
! rank 0 creates the .out file before all processes have tested for the file,
! which can result in their aborting.
   CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif

   IF ( file_exist ) THEN
     IF ( my_rank == 0 ) THEN
     write(0,*) 'COMMAS: OUTPUT FILE ',rstprefix(1:length)//trim(number)//'.out',' ALREADY EXISTS!'
     write(0,*) 'PLEASE DELETE FIRST IF YOU REALLY WANT TO OVERWRITE!'
     write(0,*) 'ABORT'
     ENDIF
     CALL COMMASMPI_ABORT()
     STOP
   ENDIF

   IF ( my_rank == 0 ) THEN
     
     OPEN(unit=105, file=rstprefix(1:length)//'.namelist', status='unknown',form='formatted',delim='APOSTROPHE')

     CALL WRITE_NAMELISTS(105)
   
     CLOSE(unit=105)
   ENDIF

     coards(1) = year
     coards(2) = month
     coards(3) = day
     coards(4) = hour
     coards(5) = minute
     coards(6) = second

! set thistory before defining a new netcdf file
   IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMMAS:  Problem setting THISTORY'

#ifdef MPI                              
      CALL GRID_SET_ATTRIBUTES(gd,rstprefix(1:length),1,member,coards  &
                              ,nx,ny,nz,microphys,v5dflds              &
                              ,nxend,nyend,nzend,my_rank)
#else
      CALL GRID_SET_ATTRIBUTES(gd,rstprefix(1:length),1,member,coards  &
                              ,nx,ny,nz,microphys,v5dflds  )
#endif

      prefix(:) = ' '
      prefix = rstprefix
      

! create checkpoint time file
       IF ( my_rank == 0 .and. member == 0 ) THEN
       open(45,file=prefix(1:length)//'.autorestart',status='unknown',form='formatted')
       write(45,*) start
       close(45)
       ENDIF


    IF ( historyoutput .eq. 1 .or. historyoutput == 3 ) THEN ! netcdf output

     IF ( historyoutput == 1 ) THEN
       ncfile = prefix(1:length)//trim(number1)//'.nc'
     ELSE
       write(stime,'(i6.6)') start
       ncfile = prefix(1:length)//'.'//stime//'.nc'
     ENDIF

#ifdef MPI

  IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
   token = 1

   if (my_rank==0) then
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)    !head node will write first
    IF ( number_of_processes .gt. 1 ) THEN
      CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
    ENDIF
   endif

   IF ( number_of_processes .gt. 1 ) THEN
   do i = 1,number_of_processes-1
    if (my_rank==i) then
     CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)          !opens/writes to nc file
     if (my_rank.ne.number_of_processes-1) CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
    endif
   enddo
   ENDIF

   CALL MPI_BARRIER(my_comm, mpi_error_code)
   
   ELSE ! parallel IO
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     IF ( parallelio_type == 1 ) THEN
#ifdef NC4
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
     IF ( my_rank==0 ) write(0,*) 'Cannot use parallelio unless netCDF4 support is enabled!'
     call commasmpi_abort()
#endif
     ELSEIF ( parallelio_type == 2 ) THEN

#ifdef USE_PNETCDF
     CALL GRID_WRITE_PNETCDF( gd, 0, ncfile) 
#else

#if defined( NC4 )
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
! #if !defined( NC4 )
     IF ( my_rank==0 ) write(0,*) 'Cannot use pnetcdf unless netCDF4 or pnetcdf support is enabled!'
     call commasmpi_abort()
#endif

#endif
     
     ENDIF
   ENDIF

   ELSEIF ( historyoutput .eq. 0 .or. historyoutput .eq. -1 ) THEN
     CALL GRID_DEFINE_NETCDF(gd, prefix(1:length)//trim(number)//'.nc')
     CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number)//'.nc') 

#else
      CALL GRID_DEFINE_NETCDF(gd,ncfile)
      CALL GRID_WRITE_NETCDF( gd, 0 ,ncfile)
#endif

     ELSEIF ( historyoutput .eq. 2 ) THEN ! fortran binary output
    
! Write a separate binary file for each tile -- faster than round-robin netcdf
     write(timestr,timfmt) start
     CALL GRID_WRITE_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )

     ENDIF

      
      CALL STRING_LIMITS(prefix, ibeg, iend)
      length = iend - ibeg + 1
    
    ENDIF


#ifdef MPI                           
!    tmpname = prefix(ibeg:iend)//trim(number)//'.out'
   IF( .not. SET_ATTRIBUTE(gd,'OUTPUT_FILE_NAME',prefix(ibeg:iend)//trim(number)//'.out') ) &
     write(0,*) 'COMMAS: error setting OUTPUT_FILE_NAME'
#endif
        
    lun = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'COMMAS:   RESTART')
    write(lun,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)
!    CALL WRITE_NAMELISTS(lun)



    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMMAS:  Problem setting TIME'  
   
   IF ( itrack <= 0 ) THEN

    IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(6,*) 'COMMAS:  Problem setting UGRID'  
    IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(6,*) 'COMMAS:  Problem setting VGRID'

   ELSE

    CALL GET_VARIABLE(gd,'UGRID',ugrid)
    CALL GET_VARIABLE(gd,'VGRID',vgrid)

    CALL GET_VARIABLE(gd,'XBAR',xbar)
    CALL GET_VARIABLE(gd,'YBAR',ybar)
    xold = xbar
    yold = ybar
!    ugrid = ugrid1
!    vgrid = vgrid1

   ENDIF

    CALL GET_VARIABLE(gd,'UGRID0',ugrid0)
    CALL GET_VARIABLE(gd,'VGRID0',vgrid0)
    CALL GET_VARIABLE(gd,'DX_STRETCH',dx_stretch)
    CALL GET_VARIABLE(gd,'DY_STRETCH',dy_stretch)
    
!    IF ( my_rank == 0 ) write(0,*) 'ugrid0,vgrid0 = ',ugrid0,vgrid0

    IF( .not. SET_VARIABLE(gd,'TRESTART',  trestart)) write(6,*) 'COMMAS:  Problem setting TRESTART'
    IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMMAS:  Problem setting THISTORY'
    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMMAS:  Problem setting TVIS5D'
    IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COMMAS:  Problem setting TPRINT'
    IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COMMAS:  Problem setting TSTAT'
    IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COMMAS:  Problem setting STOP'
    IF( .not. SET_VARIABLE(gd,'DT',        dt)      ) write(6,*) 'COMMAS:  Problem setting DT'

#ifdef MPI
!    CALL FOLLOW( gd, -1.0 )
!    CALL MPI_BARRIER(my_comm, mpi_error_code)
    
    IF ( my_rank == 0 ) THEN
      write(6,*) 'Startup: my_rank = ',my_rank
      write(6,*) 'nx,ny,nz,nxend,nyend,nzend = ',nx,ny,nz,nxend,nyend,nzend
    ENDIF
    CALL PRINT_MPI( gd )
    
!    CALL FOLLOW(gd, 1.0)
#else
    CALL PRINT( gd )
#endif

    time = start
    lun = FILE_CLOSE()

   CALL cld_cpu('I/O')

!
! reset flash arrays if this is a normal history time
!
!        print*, 'check for reset flash arrays'
    CALL GET_VARIABLE(gd, 'THISTORY',  this) 
    
    IF ( (start/this)*(this) == start ) THEN
    elec   = GET_VARIABLE_INDEX(gd, 'EX')
    precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')
       IF ( associated(gd%var(elec)%flt3d) ) THEN
!        print*, 'reset flash array,'
          gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( neelec2d > 0 ) THEN
            DO i = precip+nprecip, precip+nprecip-1+neelec2d
              gd%var(i)%flt2d(:,:) = 0.0
            ENDDO
          ENDIF
          
          ! reset net charge tendency arrays
          IF ( ichgtndadv >= 1 ) gd%var(elec+ichgtndadv-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndadvsn >= 1 ) gd%var(elec+ichgtndadvsn-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndmix >= 1 ) gd%var(elec+ichgtndmix-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndnic >= 1 ) gd%var(elec+ichgtndnic-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndsed >= 1 ) gd%var(elec+ichgtndsed-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndion >= 1 ) gd%var(elec+ichgtndion-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndlgt >= 1 ) gd%var(elec+ichgtndlgt-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndave >= 1 ) gd%var(elec+ichgtndave-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          
       ENDIF
       
       IF ( ihwindsfcmx > 0 ) THEN
         gd%var(precip-1+ihwindsfcmx)%flt2d(:,:) = 0.0
       ENDIF
       IF ( iwzsfcmax > 0 ) THEN
         gd%var(precip-1+iwzsfcmax)%flt2d(:,:) = 0.0
       ENDIF
       IF ( iwzsfcmin > 0 ) THEN
         gd%var(precip-1+iwzsfcmin)%flt2d(:,:) = 0.0
       ENDIF
     ENDIF


! iuvwadv check to reset mixing ratios to zero

!      write(0,*) 'iuvwadv = ', iuvwadv
    IF ( .false. .and. iuvwadv >= 3 ) THEN
     ! hack to test bringing in only qh and qr. Set all numbers to zero and let calcnfromq handle it
      DO i = 1,size(gd%var)
        
        IF ( .false. ) THEN
          IF ( gd%var(i)%type(1:4) .eq. 'xyz3' .and.  gd%var(i)%buotype == 3 ) THEN
            IF ( .not. ( trim(gd%var(i)%name) == 'QH' .or. trim(gd%var(i)%name) == 'QR' ) ) THEN
              gd%var(i)%flt3d(:,:,:) = 0.0
            ENDIF
          ENDIF
          
          IF ( gd%var(i)%name .eq. 'CCW' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          IF ( gd%var(i)%name .eq. 'CRW' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          IF ( gd%var(i)%name .eq. 'CCI' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          IF ( gd%var(i)%name .eq. 'CSW' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          IF ( gd%var(i)%name .eq. 'CHW' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          IF ( gd%var(i)%name .eq. 'CHL' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
          
          IF ( gd%var(i)%name .eq. 'PI' ) THEN
            gd%var(i)%flt3d(:,:,:) = 0.0
          ENDIF
  
          IF ( gd%var(i)%name .eq. 'CCCN' ) THEN
  !          gd%var(i)%flt3d(:,:,:) = ccn/1.225
          ENDIF
        
        ENDIF ! T/F
        
        IF ( gd%var(i)%name .eq. 'TH' ) THEN
          sinit  = GET_VARIABLE_INDEX(gd, 'THINIT')
          DO kz = 1,nz
          DO jy = -ng+1,ny+ng
          DO ix = -ng+1,nx+ng
          gd%var(i)%flt3d(ix,jy,kz) = gd%var(sinit)%flt1d(kz)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

        IF ( gd%var(i)%name .eq. 'QV' ) THEN
          sinit  = GET_VARIABLE_INDEX(gd, 'QVINIT')
          DO kz = 1,nz
          DO jy = -ng+1,ny+ng
          DO ix = -ng+1,nx+ng
          gd%var(i)%flt3d(ix,jy,kz) = gd%var(sinit)%flt1d(kz)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

        
      ENDDO
    
    ENDIF

!--------------------------------------------------------------------------------------
! Open VIS5D file  

   INQUIRE(file=v5dfilename, exist=file_exist)
   
   lun = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'COMMAS:')
   
   IF ( stop .gt. tvis5d ) THEN

   IF ( file_exist ) THEN
     write(lun,*) 'COMMAS:  OPEN VIS5D FILE FOR APPENDING: ',prefix(1:length)//trim(number)
      call Vis5DOPEN( gd )
   ELSE
     write(lun,*) 'COMMAS:  CREATE NEW VIS5D FILE: ',prefix(1:length)//trim(number)

     gx     = GET_VARIABLE_INDEX(gd, 'XC')
     gy     = GET_VARIABLE_INDEX(gd, 'YC')
     gz     = GET_VARIABLE_INDEX(gd, 'ZC')
#ifdef MPI

     mpisep = .true.

     call Vis5DInit( gd, gd%var(gx), gd%var(gy), gd%var(gz), v5dstridex, v5dstridez, mpisep )
#else
     call Vis5DInit( gd, gd%var(gx), gd%var(gy), gd%var(gz), v5dstridex, v5dstridez, mpisep )
#endif
   ENDIF
   
   ENDIF
   
   lun = FILE_CLOSE()

   IF ( rsttype .eq. 1 ) THEN
      filenm = prefix
   ELSE
      filenm = rstprefix
   ENDIF
   
    CALL STRING_LIMITS(filenm, ibeg, iend)
    length = iend - ibeg + 1

   IF ( my_rank == 0 ) THEN
   
    INQUIRE(file=filenm(1:length)//trim(number)//'.stat', exist=file_exist)
    IF ( file_exist ) THEN
      OPEN(unit=3, file=filenm(1:length)//trim(number)//'.stat', status='old',form='formatted',position='append')
    ELSE

     IF ( microphys(1:5) .eq. 'ICE10' .or. microphys(1:1) .eq. 'Z'    &
          .or. microphys(1:5) .eq. 'WARMZ' .or.  microphys(1:3) .eq. 'TAK') THEN
      OPEN(unit=3, file=filenm(1:length)//trim(number)//'.stat', status='new',form='formatted')
     ENDIF
    
    ENDIF


    IF ( onedoutput > 0 .and. nyend <= 2 ) THEN
      INQUIRE(file=onedfilename, exist=file_exist)
      IF ( file_exist ) THEN
        OPEN(unit=33, file=onedfilename, status='old',form='formatted',position='append')
      ELSE
        OPEN(unit=33, file=onedfilename, status='new',form='formatted')
      ENDIF
    ENDIF
    

    INQUIRE(file=filenm(1:length)//trim(number)//'.flash', exist=file_exist)
    IF ( file_exist .and. ipelec .ge. 2 .and. my_rank == 0 ) THEN
     OPEN(unit=16, file=filenm(1:length)//trim(number)//'.flash', status='old',form='formatted',position='append')
!     OPEN(unit=17, file=filenm(1:length)//trim(number)//'.flash2', status='old',form='formatted',position='append')
     write(16,'(a,i7)') '  COMMAS:   RESTART AT  T = ',start
    ELSEIF ( ipelec .ge. 2 .and. my_rank == 0 ) THEN
     OPEN(unit=16, file=filenm(1:length)//trim(number)//'.flash', status='new',form='formatted')   
!     OPEN(unit=17, file=filenm(1:length)//trim(number)//'.flash2', status='new',form='formatted')   
    ENDIF
   
   ENDIF ! my_rank == 0
   
   ! check for existence of large ion category
   IF ( lscnli > 1 ) largeion = .true.
   
  ENDIF ! start .lt. 0


!--------------------------------------------------------------------------------------
! Open I/O file for model text output

   lun = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'COMMAS: STARTING SIMULATION...')
   
! Set unit numbers for standard out (luno) and standard error (lune)  / luno and lune are stored in param_module

   luno = lun
   lune = 0


!-----------------------------------------------------------------------------
! Compute some needed variables

  IF( float(stop-start) .lt. dt ) THEN

   write(6,*) 'START TIME:  ', start
   write(6,*) 'STOP  TIME:  ', stop
   write(6,*) 'DT:          ', dt
   write(6,*) 'STOP-START < DT, ERROR...exiting...'
   
   call exit(0)

  ENDIF

!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

#ifdef MPI
  CALL GET_VARIABLE(gd, 'NXEND',    nxend)
  CALL GET_VARIABLE(gd, 'NYEND',    nyend)
  CALL GET_VARIABLE(gd, 'NZEND',    nzend)
#else
  CALL GET_VARIABLE(gd, 'NX',       nx)
  IF ( nx == 0 ) THEN
    nx = ncxe
    IF( .not. SET_VARIABLE(gd, 'NX',       nx)   ) write(6,*) 'COMMAS:  Problem setting NX'  
  ENDIF
  CALL GET_VARIABLE(gd, 'NY',       ny)
  CALL GET_VARIABLE(gd, 'NZ',       nz)
#endif
  CALL GET_VARIABLE(gd, 'DX',       dx)
  CALL GET_VARIABLE(gd, 'DY',       dy)
  CALL GET_VARIABLE(gd, 'DZ',       dz)
  CALL GET_VARIABLE(gd, 'DT',       dt)
  CALL GET_VARIABLE(gd, 'NSMALL',   nsmall)
  CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
  CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)
  CALL GET_ATTRIBUTE(gd, 'SIZE_OF_ENSEMBLE',   ensize)

! Now read some parameters back out of the data structure.  This is important for ensemble
! runs that over-ride the namelist settings.

  CALL GET_VARIABLE(gd,'RHO_QR',     rho_qr)
  CALL GET_VARIABLE(gd,'RHO_QS',     rho_qs)
  CALL GET_VARIABLE(gd,'RHO_QH',     rho_qh)
  CALL GET_VARIABLE(gd,'CNOR',       cnor)
  CALL GET_VARIABLE(gd,'CNOS',       cnos)
  CALL GET_VARIABLE(gd,'CNOH',       cnoh)
  CALL GET_VARIABLE(gd,'CCN',        ccn)
  CALL GET_VARIABLE(gd,'ALPHAH',     alphah)
  CALL GET_VARIABLE(gd,'ALPHAHL',    alphahl)
  CALL GET_VARIABLE(gd,'EHSLFO0',    ehslfo0)

  CALL GET_VARIABLE(gd,'LGTSEED',    iseed)

  
  IF ( microp(1:5) == 'ICE10' ) THEN
  CALL GET_VARIABLE(gd,'RHO_QR_10',     rho_qr_i10)
  CALL GET_VARIABLE(gd,'RHO_QS_10',     rho_qs_i10)
  CALL GET_VARIABLE(gd,'RHO_QF_10',     rho_qf_i10)
  CALL GET_VARIABLE(gd,'RHO_QGL_10',    rho_qgl_i10)
  CALL GET_VARIABLE(gd,'RHO_QGM_10',    rho_qgm_i10)
  CALL GET_VARIABLE(gd,'RHO_QGH_10',    rho_qgh_i10)
  CALL GET_VARIABLE(gd,'RHO_QH_10',     rho_qh_i10)
  CALL GET_VARIABLE(gd,'RHO_QHL_10',    rho_qhl_i10)
  CALL GET_VARIABLE(gd,'CNOR_10',       cnor_i10)
  CALL GET_VARIABLE(gd,'CNOS_10',       cnos_i10)
  CALL GET_VARIABLE(gd,'CNOF_10',       cnof_i10)
  CALL GET_VARIABLE(gd,'CNOGL_10',      cnogl_i10)
  CALL GET_VARIABLE(gd,'CNOGM_10',      cnogm_i10)
  CALL GET_VARIABLE(gd,'CNOGH_10',      cnogh_i10)
  CALL GET_VARIABLE(gd,'CNOH_10',       cnoh_i10)
  CALL GET_VARIABLE(gd,'CNOHL_10',      cnohl_i10)
  
  ENDIF
  
  ! DTD: Read in MY-related parameters
  
  IF(microp(1:2) == 'MY' ) THEN
    CALL GET_VARIABLE(gd,'NTC_MY',        ntc_my)
    CALL GET_VARIABLE(gd,'RHO_QR_MY',     rho_qr_my)
    CALL GET_VARIABLE(gd,'RHO_QI_MY',     rho_qi_my)
    CALL GET_VARIABLE(gd,'RHO_QS_MY',     rho_qs_my)
    CALL GET_VARIABLE(gd,'RHO_QG_MY',     rho_qg_my)
    CALL GET_VARIABLE(gd,'RHO_QH_MY',     rho_qh_my)
    CALL GET_VARIABLE(gd,'CNOR_MY',       cnor_my)
    CALL GET_VARIABLE(gd,'CNOS_MY',       cnos_my)
    CALL GET_VARIABLE(gd,'CNOG_MY',       cnog_my)
    CALL GET_VARIABLE(gd,'CNOH_MY',       cnoh_my)
    CALL GET_VARIABLE(gd,'ALPHAR_MY',     alphar_my)
    CALL GET_VARIABLE(gd,'ALPHAI_MY',     alphai_my)
    CALL GET_VARIABLE(gd,'ALPHAS_MY',     alphas_my)
    CALL GET_VARIABLE(gd,'ALPHAG_MY',     alphag_my)
    CALL GET_VARIABLE(gd,'ALPHAH_MY',     alphah_my)
  END IF

  IF(microp(1:3) == 'TAK' ) THEN
    CALL GET_VARIABLE(gd,'RHO_QG_TAK',     rho_qh_tak)
    CALL GET_VARIABLE(gd,'RHO_QH_TAK',     rho_qhl_tak)
  ENDIF

  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1

  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COMMAS:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COMMAS:  Problem setting NG'

! Allocate SOLVER scratch space
  
  m = (nx+2*ng)*(ny+2*ng)*(nz+2*ng)
  allocate ( st(m,0:ns) )
  allocate ( ft(m,14) ) ! temporary arrays for solver
!  ft => st(1:m,ns+1:ns+11)
  allocate ( precip_old(-ng+1:nx+ng,-ng+1:ny+ng,nprecip) )
  
  precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')
!    write(0,*) 'load precip_old, nprecip, precip, iphaildiam= ', nprecip, precip, iphaildiam
    i = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_DIAM', 1)  )
    IF  ( i > 1 ) THEN ! ( iphaildiam > 1 ) THEN
      gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
    ENDIF
    i = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_DIAM', 1)  )
    IF  ( i > 1 ) THEN ! ( iphaildiam > 1 ) THEN
      gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
    ENDIF
    i = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_DIAM', 1)  )
    IF  ( i > 1 ) THEN ! ( iphaildiam > 1 ) THEN
      gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
    ENDIF

     i = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGEH', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGER', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF

     i = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNH', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNF', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF

     DO i = precip,precip+nprecip - 1
!       write(0,*) 'load precip_old i = ',i, i - precip + 1
        precip_old(-ng+1:nx+ng,-ng+1:ny+ng,i - precip + 1) = gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng)
     ENDDO
!    write(0,*) 'load precip_old done'


  st(:,:) = 0.
  ft(:,:) = 0.

  tprt = tprint
  tstt = tstat 
  tv5d = tvis5d
  this = thistory

!  print *, '++++++', tprint, tstat, tvis5d, thistory
!  print *, '++++++', tprt, tstt, tv5d, this

! tprt = tprt * (1 + start/tprt) - 1
! tstt = tstt * (1 + start/tstt) - 1
! this = thistory * (1 + start/thistory) - 1
! tv5d = tvis5d * (1 + start/tvis5d) - 1

!-----------------------------------------------------------------------------
! PRINT OUT PARAMS

  write(lun,"(1x,80('-'))")
  write(lun,*)
  write(lun,*) 'COMMAS_1.0 INPUT DATA'
  write(lun,*)
  write(lun,"(1x,80('-'))")
  write(lun,*)
  write(lun,*) 'SIMULATION NAME            = ', prefix(1:length)//trim(number)
  write(lun,*) 'SIMULATION START TIME      = ', start
  write(lun,*) 'SIMULATION STOP  TIME      = ', stop
  write(lun,*) 'TIME STEP                  = ', dt
  write(lun,*) 'NUMBER OF TIME STEPS       = ', int((stop-start)/dt)
  write(lun,*) 'NUMBER OF SCALAR VARS      = ', ns
  write(lun,*)
  write(lun,*) 'No. of SMALL TIME STEPS    = ', nsmall
  write(lun,*)
  write(lun,*) 'PRINT INTERVAL             = ', tprint
  write(lun,*) 'HISTORY DUMP INTERVAL      = ', thistory
  write(lun,*) 'VIS5D DUMP INTERVAL        = ', tvis5d
  write(lun,*)
  write(lun,*) 'X-GRID MOTION (M/S)        = ', ugrid
  write(lun,*) 'Y-GRID MOTION (M/S)        = ', vgrid
  write(lun,*)


  IF ( my_rank == 0 .or. verbose_mpi ) THEN
    CALL PRINT_ALL( gd )
  ENDIF

  IF ( itraj > 0 ) THEN
! for balloon trajectories:
   IF ( chgavex <= 0 ) chgavex = Min(ng, Max(3, Nint(750./dx) ) )
   chgavex = Min(nx-1, Min(ny-1,chgavex) )
   write(lun,*) 'chgavex = ',chgavex

   IF ( rsttype .eq. 1 ) THEN
      filenm = prefix
   ELSE
      filenm = rstprefix
   ENDIF
    IF ( my_rank == 0 .and. itraj >= 2 ) THEN
     IF ( ntrajtype == 1 ) THEN
     OPEN(unit=trj_unit, file=filenm(1:length)//trim(number)//'.traj', status='unknown',form='formatted')
     ELSE
      DO i = 1,ntrajtype
       write(trjnum,'(i1)') i
       OPEN(unit=trj_unit+i-1, file=filenm(1:length)//trim(number)//'.traj'//trjnum, status='unknown',form='formatted')
      ENDDO
     ENDIF
    ENDIF
    CALL TRAJ_INIT(dt,nxend,nyend,nzend,nx,ny,nz,ng,filenm,number,start)
  ENDIF
#ifdef MPI
   CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif

!-----------------------------------------------------------------------------
! GET the variables needed for the integration 

 u      = GET_VARIABLE_INDEX(gd, 'U')
 v      = GET_VARIABLE_INDEX(gd, 'V')
 w      = GET_VARIABLE_INDEX(gd, 'W')
 pi     = GET_VARIABLE_INDEX(gd, 'PI')
 rho    = GET_VARIABLE_INDEX(gd, 'RHO')
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
     nxtra = 1
  ENDIF
 
 uinit0 = GET_VARIABLE_INDEX(gd, 'UINIT0')
 vinit0 = GET_VARIABLE_INDEX(gd, 'VINIT0')
 uinit  = GET_VARIABLE_INDEX(gd, 'UINIT')
 vinit  = GET_VARIABLE_INDEX(gd, 'VINIT')
 winit  = GET_VARIABLE_INDEX(gd, 'WINIT')
 piinit = GET_VARIABLE_INDEX(gd, 'PIINIT')
 kminit = GET_VARIABLE_INDEX(gd, 'KMINIT')
 kmbaserm = GET_VARIABLE_INDEX(gd, 'KMBASERM')
 sinit  = GET_VARIABLE_INDEX(gd, 'THINIT')

 gx     = GET_VARIABLE_INDEX(gd, 'XC')
 gy     = GET_VARIABLE_INDEX(gd, 'YC')
 gz     = GET_VARIABLE_INDEX(gd, 'ZC')
 precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')
 
 i = GET_VARIABLE_INDEX(gd, 'CRW',1) 

   IF ( microphys(1:5) .eq. 'ICE10' .or. microphys(1:1) .eq. 'Z' .or. microphys(1:5) .eq. 'WARMZ' ) THEN
     IF ( my_rank == 0 ) THEN
     write(3,'(a,3(i6,1x),i3)') 'nx,ny,nz,ns = ',nxend,nyend,nzend,ns
     write(3,'(a,f5.2,3(1x,f9.2))') 'dt,dx,dy,dz = ',dt,dx,dy,dz
     write(3,'(a,a)') 'micro = ',microphys
     write(3,'(a,2(1x,i3))') 'lqb,lqe= ',lqb,lhab
     write(3,'(a)') ' KZ, ZC, DZC:'
     DO k = 1,nz-1
       write(3,'(1x,i4,f8.0,1x,f10.3)') k, gd%var(gz)%flt1d(k), 1./gd%var(gz+2)%flt1d(k)
     ENDDO
     ENDIF
     
   ENDIF


   allocate( z1d4(nzend,4) )
   
     DO i = 1,4
      DO k = 1,nz
       z1d4(kzbeg-1+k,i) = gd%var(gz+i-1)%flt1d(k)
      ENDDO
     ENDDO
#ifdef MPI
      IF ( nprock > 1 ) THEN
 ! communicate temc array using max from each level
      allocate( mpitotinth((nzend),4))
      allocate(mpitotoutth((nzend),4))
      
      mpitotinth(:,:) = 0.0d0
     
      DO i = 1,4
       DO k = kzbeg,kzend
        mpitotinth(k,i) =  z1d4(k,i)
       ENDDO
      ENDDO

      n = 4*(nzend)

      CALL MPI_AllReduce(mpitotinth, mpitotoutth, n, MPI_DOUBLE_PRECISION, MPI_Max, my_comm, mpi_error_code)
      
      DO i = 1,4
       DO k = 1,nzend
        z1d4(k,i) = mpitotoutth(k, i)
       ENDDO
      ENDDO
       
       deallocate( mpitotinth )
       deallocate( mpitotoutth )
      
      ENDIF
#endif

!-----------------------------------------------------------------------------
! Print out date and time

  CALL DATE_AND_TIME(date,hhmmss)

  read(date(5:6),'(i2)') month

  write(lun,"(1x,80('-'))")
  write(lun,*) 
  write(lun,*) 'WALLCLOCK TIME FOR START OF COMMAS MODEL RUN:  ',months(month),' ',date(7:8),' ',date(1:4)
  write(lun,*) 'WALLCLOCK TIME FOR START OF COMMAS MODEL RUN:  ',hhmmss(1:2),':',hhmmss(3:4),':',hhmmss(5:6)
  write(lun,*) 
  write(lun,"(1x,80('-'))")
  write(lun,*) 
    
! LJW 08/22/08:  Added code for millisecond time counter from TICK_MODULE
!-----------------------------------------------------------------------------

!  print *, 'got here'
  IF( .NOT. TICK_INIT(float(start), float(stop), float(time), dt) ) THEN
   print *, 'TICK_INITIALIZATION FAILED!'
   call exit(1)
  ENDIF

! Main time step loop
!-----------------------------------------------------------------------------


! Stat output

   IF ( (microphys(1:5) .eq. 'ICE10' .or. &
         microphys(1:1) .eq. 'Z' .or.     &
         microphys(1:5) .eq. 'WARMZ' .or. &
         microphys(1:3) .eq. 'TAK' ) ) THEN

     CALL cld_cpu('I/O')

     CALL stats( nx,ny,nz,ns,                                 &
                 time,dx,dy,dz,                               &
                 gd%var(gx),gd%var(gy),gd%var(gz),            &
                 ft(1,1),gd%var(sinit),gd%var(s),             &
                 gd%var(u), gd%var(v), gd%var(w),             &
                 gd%var(uinit) ,gd%var(vinit), gd%var(winit), &
                 gd%var(pi), gd%var(piinit),                  &
                 gd%var(km), gd%var(kminit),                  &
                 z1d4,                                        &
                 1, nx-1, 1, ny-1)

     CALL cld_cpu('I/O')

   ENDIF


MAIN: DO WHILE ( .TRUE. )

! LJW 08/22/08:  Increment time to make sure that we are supposed to continue
!----------------------------------------------------------------------------

   IF( TICK_TOCK() ) EXIT

   time = TICK_TIME()
   time_real = TICK_TIME(.true.)/time_scale

   tim1 = commas_dtime2(tim2,0)

!----------------------------------------------------------------------
! Time and location of the grid
 
   x_sw_loc = x_sw_loc + ugrid*dt
   y_sw_loc = y_sw_loc + vgrid*dt

   IF ( ( microphys(1:5) .eq. 'ICE10' .or.          &
          microphys(1:1) .eq. 'Z'     .or.          &
          microphys(1:3) .eq. 'HCM'   .or.          &
          microphys(1:3) .eq. 'TAK'   .or.          &
          microphys(1:8) .eq. 'WARMZIEG'   )        &
#ifdef MPI
          .and. my_rank == 0                        &
#endif
          ) THEN
     write(luno,'(a)') '============================================='
     write(luno,'(a)') 'NSTEP, NSTART, NSTOP, TIME'
!     write(luno,'(1x,i7,1x,i7,1x,i7,1x,i6,".000")')  int(time/dt), Max(1,int(start/dt)), int((stop)/dt), time
     write(luno,'(1x,i7,1x,i7,1x,i7,1x,f10.3)')  int(time/dt), Max(1,int(start/dt)), int((stop)/dt), time_real
     write(luno,'(a)') '============================================='
   ENDIF

   IF( .not. SET_VARIABLE(gd,'TIME',            time)  ) write(6,*) 'COMMAS:  Problem setting TIME'
   IF( .not. SET_VARIABLE(gd,'XG_POS',      x_sw_loc)  ) write(6,*) 'COMMAS:  Problem setting XG_POS'
   IF( .not. SET_VARIABLE(gd,'YG_POS',      y_sw_loc)  ) write(6,*) 'COMMAS:  Problem setting YG_POS'
   IF( .not. SET_VARIABLE(gd,'UGRID',          ugrid)  ) write(6,*) 'COMMAS:  Problem setting UGRID'  
   IF( .not. SET_VARIABLE(gd,'VGRID',          vgrid)  ) write(6,*) 'COMMAS:  Problem setting VGRID'

   IF( .not. SET_VARIABLE(gd,'TRESTART',        time)  ) write(6,*) 'COMMAS:  Problem setting TRESTART'

!----------------------------------------------------------------------

   IF( ALARM(tv5d) ) THEN
     iv5dwritten(:) = 0
     lv5dwrite = .true.
     v5dtimes(itv5d) = time
   ELSE
     lv5dwrite = .false.
   ENDIF
   
   IF ( ALARM(this) .or. float(stop-time) .lt. dt ) THEN
     lncwrite = .true.
   ELSE
     lncwrite = .false. 
   ENDIF

   IF ( float(stop-time) .lt. dt ) THEN
     lrst = .true.
   ELSE
     lrst = .false. 
   ENDIF

   IF( ALARM(tprt) ) THEN
      lprt = .true.
   ELSE
      lprt = .false.
   ENDIF

   IF( (ALARM(tstt) .or. ALARM(start)) .and. tstat .gt. 0 ) THEN
     lstt = .true.
   ELSE
     lstt = .false.
   ENDIF
  
   io_flag  = .false.

   IF( ALARM(tv5d) .or.  &
       ALARM(this) .or.  &
       ALARM(tprt) .or.  &
!      ALARM(trst) .or.  &
       float(stop-time) .lt. dt    ) io_flag = .true.
       

!-----------------------------------------------------------------------------
! SUBTRACT OUT GRID MOTION -- do this every time step to avoid round-off
!  differences that could arise if only done for history dumps/restarts.


    CALL cld_cpu('TIMESTEP')  
     CALL FOLLOW(gd, -1.0)

     IF ( numbss > 0 ) THEN
       CALL BSS(gd, time_real )
     ENDIF
     
     IF ( clearflash ) THEN

       clearflash = .false.

       IF ( ihwindsfcmx > 0 ) THEN
         gd%var(precip-1+ihwindsfcmx)%flt2d(:,:) = 0.0
       ENDIF
       IF ( iwzsfcmax > 0 ) THEN
         gd%var(precip-1+iwzsfcmax)%flt2d(:,:) = 0.0
       ENDIF
       IF ( iwzsfcmin > 0 ) THEN
         gd%var(precip-1+iwzsfcmin)%flt2d(:,:) = 0.0
       ENDIF
       IF ( ihailmax2d > 0 ) THEN
         gd%var(precip-1+ihailmax2d)%flt2d(:,:) = 0.0
       ENDIF
       IF ( ihailmaxk1 > 0 ) THEN
         gd%var(precip-1+ihailmaxk1)%flt2d(:,:) = 0.0
       ENDIF
     
     IF ( associated(gd%var(elec)%flt3d)  ) THEN
     
!     print*, 'reset flash array,'
      gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
      gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
      gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz) = 0.0
      IF ( neelec2d > 0 ) THEN
        DO i = precip+nprecip, precip+nprecip+neelec2d-1
          gd%var(i)%flt2d(:,:) = 0.0
        ENDDO
      ENDIF

          ! reset net charge tendency arrays
          IF ( ichgtndadv >= 1 ) gd%var(elec+ichgtndadv-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndadvsn >= 1 ) gd%var(elec+ichgtndadvsn-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndmix >= 1 ) gd%var(elec+ichgtndmix-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndnic >= 1 ) gd%var(elec+ichgtndnic-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndsed >= 1 ) gd%var(elec+ichgtndsed-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndion >= 1 ) gd%var(elec+ichgtndion-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndlgt >= 1 ) gd%var(elec+ichgtndlgt-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
          IF ( ichgtndave >= 1 ) gd%var(elec+ichgtndave-1)%flt3d(1:nx,1:ny,1:nz) = 0.0
      

     ENDIF
     
     ENDIF

!      IF ( xtra > wz ) THEN
!       gd%var(xtra)%flt3d(1:nx,1:ny,1:nz) = 0.0
!      ENDIF

!   IF ( bogusvortex >= 1 .and. ( time_real <= timint + dt*2 )  ) THEN
!     CALL HURRFORCE(gd, time_real)
!   ENDIF


    CALL cld_cpu('TIMESTEP')  

#ifdef MPI
!     IF (debug_mpi)  write(0,"(A,i2,A,i4)") "COMMAS: PROC ",my_rank," ABOUT TO CALL SOLVER FOR t=",time
#else
!    write(0,"(A,i6)") "COMMAS: ABOUT TO CALL SOLVER FOR t=",time
#endif

   CALL cld_cpu('SOLVER')

     CALL SOLVER(gd,                                 &
                 gd%var(u),  gd%var(uinit) ,         &            ! U,  UINIT
                 gd%var(v),  gd%var(vinit) ,         &            ! V,  VINIT
                 gd%var(w),  gd%var(winit) ,         &            ! W,  WINIT
                 gd%var(pi), gd%var(piinit),         &            ! PI, PIINIT
                 gd%var(km), gd%var(kminit),         &            ! KM, KINIT
                 gd%var(rho),                        &            ! rho (air density)
                 gd%var(kmbaserm),                   &
                 gd%var(s),  gd%var(sinit),          &            ! S,  SINIT
                 gd%var(uinit0), gd%var(vinit0),     &            ! UINIT0, VINIT0
                 gd%var(precip),                     &            ! PRECIP
                 gd%var(gx),                         &            ! XCNTR, XEDGE, DXC, DXE 
                 gd%var(gy),                         &            ! YCNTR, YEDGE, DYC, DYE
                 gd%var(gz),                         &            ! ZCNTR, ZEDGE, DZC, DZE
                 dt,                                 &            ! DT
                 ugrid, vgrid,                       &            ! GRID MOTION
                 ugrid0, vgrid0,                     &            ! GRID MOTION at t=0
                 lat, lon,                           &            ! Latitude and Longitude
                 nsmall,                             &            ! NSMALL
                 microphys,                          &
                 nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS (mpi = nxt,nyt,nzt)
                 io_flag,                            &            ! IO_FLAG
                 st(1,1),ft(1,1),ft(1,2),ft(1,3),ft(1,4), &
                 ft(1,5),ft(1,6),ft(1,7),ft(1,8),    &
                 ft(1,9),st(1,0),ft(1,10),ft(1,11),  &
                 ft(1,12),ft(1,13),ft(1,14),         &
                 precip_old,                         &
                 gd%var(dbz), gd%var(vzf), gd%var(wz),gd%var(xtra),  &
                 gd%var(elec), gd%var(cion),         &
                 gd%var(muz),                        &
                 x_sw_loc,y_sw_loc,                  &
                 time, time_real, tstat, lstt, stop, z1d4, onedoutput, this, &
                 run_file)

#ifdef MPI
!     IF (debug_mpi)  write(0,"(A,i2,A,i4)") "COMMAS: PROC ",my_rank," DONE CALLING SOLVER FOR t=",time
#else
!    write(0,"(A,i6)") "COMMAS: DONE CALLING SOLVER FOR t=",time
#endif


   CALL cld_cpu('SOLVER')
!-----------------------------------------------------------------------------
! DTREND the pressure field
!
!    CALL cld_cpu('TIMESTEP')
!
!      CALL DTREND( gd )
!
!    CALL cld_cpu('TIMESTEP')
!
! LJW 08/22/08:  Add alarms to code for various dumps
!-----------------------------------------------------------------------------
! PRINTING calls

   IF( ALARM(tprt) ) THEN

     CALL cld_cpu('I/O')

#ifdef MPI
     CALL PRINT_MPI( gd )
#else
     CALL PRINT( gd )
#endif

     CALL cld_cpu('I/O')

   ENDIF

! Stat output

   IF ( lstt .and. (microphys(1:5) .eq. 'ICE10' .or. &
                    microphys(1:1) .eq. 'Z' .or.     &
                    microphys(1:5) .eq. 'WARMZ' .or. &
                    microphys(1:3) .eq. 'TAK' ) ) THEN

     CALL cld_cpu('I/O')

     CALL stats( nx,ny,nz,ns,                                 &
                 time,dx,dy,dz,                               &
                 gd%var(gx),gd%var(gy),gd%var(gz),            &
                 ft(1,1),gd%var(sinit),gd%var(s),             &
                 gd%var(u), gd%var(v), gd%var(w),             &
                 gd%var(uinit) ,gd%var(vinit), gd%var(winit), &
                 gd%var(pi), gd%var(piinit),                  &
                 gd%var(km), gd%var(kminit),                  &
                 z1d4,                                        &
                 1, nx-1, 1, ny-1)

     CALL cld_cpu('I/O')

   ENDIF

!-----------------------------------------------------------------------------
! Vis5D OUTPUT
!

   IF ( lv5dwrite ) THEN
 
     CALL cld_cpu('I/O')

      CALL V5DOUT(gd,                                 &
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
!                 st,                                 &
                  gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,gd%var(s)%index),  &
                  ft(1,1),ft(1,2),ft(1,3),ft(1,4),    &
                  gd%var(dbz), gd%var(vzf), gd%var(wz), gd%var(elec),gd%var(xtra), onedoutput )

      itv5d = itv5d + 1

      iv5dwritten(:) = 0

     CALL cld_cpu('I/O')
   
   ENDIF

!-----------------------------------------------------------------------------
! ADD BACK GRID MOTION 

   CALL cld_cpu('FOLLOW')
    CALL FOLLOW(gd, 1.0)
   CALL cld_cpu('FOLLOW')

! ------------------------
! Meso tracking

! set an alarm for 'track'
! call 'box'

   IF ( itrack > 0 ) THEN
   IF ( ALARM( itrack ) ) THEN
      uold = ugrid
      vold = vgrid
      track_tau = float(itrack)

      CALL BOX(gd, gd%var(u),gd%var(v),gd%var(w),gd%var(gx),gd%var(gy),gd%var(gz), &
                     xold,yold,uold,vold,ugrid,vgrid, gd%var(piinit), gd%var(sinit),  &
                     track_xmid,track_ymid, track_tau, track_alpha, track_zeta0,      &
                     track_xmid_delta, track_ymid_delta,luno)

    xbar = xold
    ybar = yold
    IF( .not. SET_VARIABLE(gd,'UGRID',        ugrid)  ) write(6,*) 'COMMAS:  Problem setting UGRID'  
    IF( .not. SET_VARIABLE(gd,'VGRID',        vgrid)  ) write(6,*) 'COMMAS:  Problem setting VGRID'
    IF( .not. SET_VARIABLE(gd,'XBAR',          xbar)  ) write(6,*) 'COMMAS:  Problem setting XBAR'  
    IF( .not. SET_VARIABLE(gd,'YBAR',          ybar)  ) write(6,*) 'COMMAS:  Problem setting YBAR'  

   ENDIF ! alarm
   ENDIF ! itrack > 0
!-----------------------------------------------------------------------------
! HISTORY DUMPS

   IF( alarm(this) ) THEN

     thislast = time


    ! ichgtndadv has both advection and mixing, so take mixing out (2017/10/2)
    ! Not any more!
    IF ( ichgtndadv >= 1 .and. ichgtndmix >= 1 ) THEN 
!       gd%var(elec+ichgtndadv-1)%flt3d(1:nx,1:ny,1:nz) = gd%var(elec+ichgtndadv-1)%flt3d(1:nx,1:ny,1:nz) - & 
!                                                         gd%var(elec+ichgtndmix-1)%flt3d(1:nx,1:ny,1:nz)
    ENDIF
           
        
     CALL cld_cpu('I/O')

    IF ( historyoutput == 3 ) THEN ! netcdf output separate file for each history time

       write(stime,'(i6.6)') time
       ncfile = prefix(1:length)//'.'//stime//'.nc'

#ifdef MPI

  IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
   token = 1

   if (my_rank==0) then
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)    !head node will write first
    IF ( number_of_processes .gt. 1 ) THEN
      CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
    ENDIF
   endif

   IF ( number_of_processes .gt. 1 ) THEN
   do i = 1,number_of_processes-1
    if (my_rank==i) then
     CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)          !opens/writes to nc file
     if (my_rank.ne.number_of_processes-1) CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
    endif
   enddo
   ENDIF

   CALL MPI_BARRIER(my_comm, mpi_error_code)
   
   ELSE ! parallel IO
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     IF ( parallelio_type == 1 ) THEN
#ifdef NC4
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
     IF ( my_rank==0 ) write(0,*) 'Cannot use parallelio unless netCDF4 support is enabled!'
     call commasmpi_abort()
#endif
     ELSEIF ( parallelio_type == 2 ) THEN

#ifdef USE_PNETCDF
     CALL GRID_WRITE_PNETCDF( gd, 0, ncfile) 
#else

#if defined( NC4 )
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
! #if !defined( NC4 )
     IF ( my_rank==0 ) write(0,*) 'Cannot use pnetcdf unless netCDF4 or pnetcdf support is enabled!'
     call commasmpi_abort()
#endif

#endif
     
     ENDIF
   ENDIF


#else
      CALL GRID_DEFINE_NETCDF(gd,ncfile)
      CALL GRID_WRITE_NETCDF( gd, 0 ,ncfile)
#endif

  ENDIF ! historyoutput == 3


     IF ( historyoutput .eq. 1 ) THEN
     
     ncfile = prefix(1:length)//trim(number1)//'.nc'

#ifdef MPI
      IF ( .not. parallelio ) THEN

!token ring to only allow one process to access nc file at a time.
       token = 1

       IF (my_rank==0) THEN
         CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
         IF ( number_of_processes .gt. 1 ) THEN
           CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
         ENDIF
       ENDIF

       IF ( number_of_processes .gt. 1 ) THEN

        DO i = 1,number_of_processes-1

         IF (my_rank==i) THEN
          CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
          CALL GRID_WRITE_NETCDF( gd, start, prefix(1:length)//trim(number1)//'.nc')      !opens/writes to nc file

          IF (my_rank.ne.number_of_processes-1) THEN
           CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
          ENDIF
         ENDIF
        ENDDO
       
       ENDIF

       CALL MPI_BARRIER(my_comm, mpi_error_code)

      ELSE  ! parallelio
!       CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
     IF ( parallelio_type == 1 ) THEN
#ifdef NC4
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
     IF ( my_rank==0 ) write(0,*) 'Cannot use parallelio unless netCDF4 support is enabled!'
     call commasmpi_abort()
#endif
     ELSEIF ( parallelio_type == 2 ) THEN

#ifdef USE_PNETCDF
     CALL GRID_WRITE_PNETCDF( gd, 0, ncfile) 
#else

#if defined( NC4 )
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
! #if !defined( NC4 )
     IF ( my_rank==0 ) write(0,*) 'Cannot use pnetcdf unless netCDF4 or pnetcdf support is enabled!'
     call commasmpi_abort()
#endif

#endif
     
     ENDIF
      ENDIF
   
     ELSEIF ( historyoutput .eq. 0 .or. historyoutput .eq. -1 ) THEN
      CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number)//'.nc') 

      CALL MPI_BARRIER(my_comm, mpi_error_code)
#else
      CALL GRID_WRITE_NETCDF( gd, start )
#endif

     ELSEIF ( historyoutput .eq. 2 ) THEN ! fortran binary output
    
! Write a separate binary file for each tile -- faster than round-robin netcdf

      write(timestr,timfmt) time
      CALL GRID_WRITE_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )

#ifdef MPI
      CALL MPI_BARRIER(my_comm, mpi_error_code)
#endif
    

     ENDIF

     CALL cld_cpu('I/O')

     CALL cld_cpu('MAIN')  

! Dump out cpu timing
 
     CALL cld_cpu('-1')

     CALL cld_cpu('MAIN')  

!
! set flag to reset flash arrays and max/min 2d arrays.  Dont do it right here because then fields will be zero in 
! the history file
!     IF ( associated(gd%var(elec)%flt3d)  ) THEN
     
       clearflash = .true.

!     ENDIF

! update the time in the checkpoint time file
       IF ( my_rank == 0 .and. member == 0 ) THEN
       open(45,file=prefix(1:length)//'.autorestart',status='unknown',form='formatted')
       write(45,*) thislast
       close(45)
       ENDIF
       
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'HAIL_DIAM', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'GR_DIAM', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'FD_DIAM', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng) = 0.0
     ENDIF
     DO i = precip,precip+nprecip-1
        precip_old(:,:,i - precip + 1) = gd%var(i)%flt2d(:,:)
     ENDDO

     i = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGEH', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'DBZCHANGER', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF

     i = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNH', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF
     i = Max( 0, GET_VARIABLE_INDEX(gd, 'CHLCNF', 1)  )
     IF ( i > 1 ) THEN
       gd%var(i)%flt3d(-ng+1:nx+ng,-ng+1:ny+ng,1:nz) = 0.0
     ENDIF


   ENDIF  ! ... ENDIF THIS DUMP...

!#ifdef MPI
!   dt1 = MPI_Wtime()
!#endif
   tim1 = commas_dtime2(tim2,1)

   IF ( microphys(1:5) .eq. 'ICE10' .or.           &
        microphys(1:1) .eq. 'Z'     .or.           &
!        microphys(1:4) .eq. 'THOM'  .or.           &
        microphys(1:3) .eq. 'HCM'   .or.           &
        microphys(1:3) .eq. 'TAK'   .or.           &
        microphys(1:8) .eq. 'WARMZIEG'   ) THEN
      IF ( my_rank == 0 ) THEN
        write(luno,'(a)') '============================================='
        write(luno,'(a,i7,1x,i6)') 'END OF STEP: nstep, time = ',int(time/dt),time
        write(luno,'(a,f13.5)') 'elapsed time = ',tim1
      ENDIF
   ENDIF

!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO MAIN
   
   DEALLOCATE(st)
   DEALLOCATE(ft)
 
!-----------------------------------------------------------------------------
! If done with processing, dump out restart file and plots, etc.
 
 
!  IF( (stop - start) .le. trestart .and. time .gt. thislast ) THEN

   IF( time .gt. thislast .and. .not. ( ensize .gt. 1 .and. restart_separate ) ) THEN
    
    CALL cld_cpu('I/O')

! update the time in the checkpoint time file
       IF ( my_rank == 0 .and. member == 0 ) THEN
       open(45,file=prefix(1:length)//'.autorestart',status='unknown',form='formatted')
       write(45,*) time
       close(45)
       ENDIF


    IF ( historyoutput .eq. 1 ) THEN

#ifdef MPI
     IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
      token = 1

      IF (my_rank==0) THEN
        CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
        IF ( number_of_processes .gt. 1 ) THEN
          CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
        ENDIF
      ENDIF

      CALL MPI_BARRIER(my_comm, mpi_error_code)

      IF ( number_of_processes .gt. 1 ) THEN
       DO i = 1,number_of_processes-1
        IF (my_rank==i) THEN
         CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
         CALL GRID_WRITE_NETCDF( gd, start, prefix(1:length)//trim(number1)//'.nc')      !opens/writes to nc file
         IF (my_rank.ne.number_of_processes-1) CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
        ENDIF

        CALL MPI_BARRIER(my_comm, mpi_error_code)

       ENDDO
      ENDIF

     ELSE  ! parallelio
      CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
     ENDIF

    ELSEIF ( historyoutput .eq. 0 .or. historyoutput .eq. -1 ) THEN
     CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number)//'.nc') 

#else
     CALL GRID_WRITE_NETCDF( gd, start )
#endif

    ELSEIF ( historyoutput == 3 ) THEN
       write(stime,'(i6.6)') time
       ncfile = prefix(1:length)//'.'//stime//'.nc'

#ifdef MPI

  IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
   token = 1

   if (my_rank==0) then
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)    !head node will write first
    IF ( number_of_processes .gt. 1 ) THEN
      CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
    ENDIF
   endif

   CALL MPI_BARRIER(my_comm, mpi_error_code)

   IF ( number_of_processes .gt. 1 ) THEN
   do i = 1,number_of_processes-1
    if (my_rank==i) then
     CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile)          !opens/writes to nc file
     if (my_rank.ne.number_of_processes-1) CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
    endif
    CALL MPI_BARRIER(my_comm, mpi_error_code)

   enddo
   ENDIF

   CALL MPI_BARRIER(my_comm, mpi_error_code)
   
   ELSE ! parallel IO
     CALL GRID_DEFINE_NETCDF(gd, ncfile)
     IF ( parallelio_type == 1 ) THEN
#ifdef NC4
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
     IF ( my_rank==0 ) write(0,*) 'Cannot use parallelio unless netCDF4 support is enabled!'
     call commasmpi_abort()
#endif
     ELSEIF ( parallelio_type == 2 ) THEN

#ifdef USE_PNETCDF
     CALL GRID_WRITE_PNETCDF( gd, 0, ncfile) 
#else

#if defined( NC4 )
     CALL GRID_WRITE_NETCDF( gd, 0, ncfile) 
#else
! #if !defined( NC4 )
     IF ( my_rank==0 ) write(0,*) 'Cannot use pnetcdf unless netCDF4 or pnetcdf support is enabled!'
     call commasmpi_abort()
#endif

#endif
     
     ENDIF
   ENDIF


#else
      CALL GRID_DEFINE_NETCDF(gd,ncfile)
      CALL GRID_WRITE_NETCDF( gd, 0 ,ncfile)
#endif

    ELSEIF ( historyoutput .eq. 2 ) THEN

! Write a separate binary file for each tile --  faster than round-robin netcdf

     write(timestr,timfmt) time
     CALL GRID_WRITE_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )
    
    ENDIF ! historyoutput 

    CALL cld_cpu('I/O')

   ELSEIF ( ensize .gt. 1 ) THEN ! write extra restart time for EnKF
    
    
    CALL cld_cpu('I/O')
   
    IF ( restart_separate ) THEN
    
      writing_restart = .true.
        CALL GRID_DEFINE_NETCDF(gd, prefix(1:length)//'.r'//trim(number1)//'.nc')
        CALL GRID_WRITE_NETCDF( gd, time, prefix(1:length)//'.r'//trim(number1)//'.nc')
      writing_restart = .false.
    
    ELSE

      IF( .not. SET_VARIABLE(gd,'TRESTART',   time+1) ) write(6,*) 'COMMAS:  Problem setting TRESTART'
      IF( .not. SET_VARIABLE(gd,'TIME',       time+1) ) write(6,*) 'COMMAS:  Problem setting TIME'

#ifdef MPI
    IF ( historyoutput .eq. 1 ) THEN

    IF ( .not. parallelio ) THEN
!token ring to only allow one process to access nc file at a time.
     token = 1

     IF (my_rank==0) THEN
      CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    !head node will write first
       IF ( number_of_processes .gt. 1 ) THEN
        CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)      !pass token to next proc
       ENDIF
     ENDIF

     CALL MPI_BARRIER(my_comm, mpi_error_code)

     IF ( number_of_processes .gt. 1 ) THEN
      DO i = 1,number_of_processes-1
       IF (my_rank==i) THEN
        CALL MPI_Recv(token,1,MPI_INTEGER,my_rank-1,6000,my_comm,mpi_status,mpi_error_code) !next proc recv's token
        CALL GRID_WRITE_NETCDF( gd, start, prefix(1:length)//trim(number1)//'.nc')      !opens/writes to nc file
        IF (my_rank.ne.number_of_processes-1) THEN
         CALL MPI_Send(token,1,MPI_INTEGER,my_rank+1,6000,my_comm,mpi_error_code)            !passes token to next proc
        ENDIF
       ENDIF
      CALL MPI_BARRIER(my_comm, mpi_error_code)
      ENDDO
     ENDIF

     ELSE  ! parallelio
      CALL GRID_WRITE_NETCDF( gd, 0, prefix(1:length)//trim(number1)//'.nc')    ! all write together
     ENDIF

    ELSEIF ( historyoutput .eq. 0 ) THEN
      write(0,*) 'CANNOT USE historyoutput=0 for EnKF!'
      CALL COMMASMPI_ABORT()

    ELSEIF ( historyoutput .eq. 2 ) THEN

! Write a separate binary file for each tile --  faster than round-robin netcdf

     write(timestr,timfmt) time+1
     CALL GRID_WRITE_BINARY(gd, prefix(1:length)//trim(number)//'.'//timestr//'.bin' )
    
    ENDIF ! historyoutput 

#else
     CALL GRID_WRITE_NETCDF( gd, start )
#endif
    
    ENDIF ! restart_separate
    
    CALL cld_cpu('I/O')

   ENDIF
 
  IF ( itraj == 1 .or. itraj == 2 ) THEN
   IF ( rsttype .eq. 1 ) THEN
      filenm = prefix
   ELSE
      filenm = rstprefix
   ENDIF
    IF ( my_rank == 0 ) THEN
!    CALL TRAJ_CLOSE(filenm,number)
    ENDIF
  ENDIF

!-----------------------------------------------------------------------------
! Print out date and time

  CALL DATE_AND_TIME(date,hhmmss)

  read(date(5:6),'(i2)') month

  write(luno,"(1x,80('-'))")
  write(luno,*) 
  write(luno,*) 'WALLCLOCK TIME FOR END OF COMMAS MODEL RUN:  ',months(month),' ',date(7:8),' ',date(1:4)
  write(luno,*) 'WALLCLOCK TIME FOR END OF COMMAS MODEL RUN:  ',hhmmss(1:2),':',hhmmss(3:4),':',hhmmss(5:6)
  write(luno,*) 
  write(luno,"(1x,80('-'))")
  write(luno,*) 

  CALL cld_cpu('MAIN')  

! Dump out cpu timing

  CALL cld_cpu('-1')
  CALL cld_cpu('-3')

  write(luno,*) 'Integration Done'

  IF ( itraj > 0 .and. my_rank == 0 ) THEN
   DO i = 1,ntrajtype
   CLOSE(trj_unit+i-1)
   ENDDO
  ENDIF

!-----------------------------------------------------------------------------
! MPI SHUTDOWN

  CALL COMMASMPI_SHUTDOWN()

  luno = FILE_CLOSE()

 STOP
 END
