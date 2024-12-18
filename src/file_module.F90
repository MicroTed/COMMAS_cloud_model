MODULE FILE_MODULE

 implicit none
 
   logical :: restart_separate = .true. ! whether to write separate restart files; 
                                        ! for ensembles only (must have ne > 1)
   logical :: writing_restart = .false.
   integer :: ienkf           = 0
!   character(LEN = 25) :: runname
!   character(LEN = 60) :: path
   integer             :: lun                     = -1           ! output unit
   
   integer             :: membernumber

   integer             :: historyoutput = 1      ! 0=netcdf for each MPI tile; 1=single netcdf, 2=fortran binary
   integer             :: historyinput  = -1     ! 0=netcdf for each MPI tile; 1=single netcdf, 2=fortran binary
   integer             :: restartformat = 1      ! 0=netcdf for each MPI tile; 1=single netcdf, 2=fortran binary
   logical             :: parallelio = .false.
   logical             :: parallel_compress = .false.
   logical             :: parallel_compress_on = .false.
   integer             :: parallelio_in = -1   ! used for restarting from a netcdf4 file but output will be pnetcdf
   integer             :: netcdfversion = 4
   integer             :: parallelio_type = 2  ! 1=netcdf4/hdf5; 2=pnetcdf
   integer             :: nclibversion = 0
   integer             :: pnetcdf_flag = 0 ! used for setting file mode with pnetcdf depending on netcdf version
   integer             :: shuffle = 1, deflate= 1, deflate_level = 2
   integer             :: chunkalg = 2
   integer             :: onedoutput = 0  ! for outputting a column of data from 2D runs

   
   INTEGER           ::  ihttima = 6
   CHARACTER(len= 6) ::  timfmt  = '(i6.6)'
   character(len=40) ::  filehead          ! base name for history files
   integer           ::  lfilehead         ! number of non-blank characters in filehead
   character(len=60) ::  hdrcfs = CHAR(0)  ! storage directory
   integer           ::  lhdrcfs = 0       ! length of hdrcfs

#ifndef CM1

 CONTAINS

!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
 INTEGER FUNCTION FILE_OPEN(gd, filetype, message, destroy)

  USE GRID_MODULE
  USE COMMASMPI_MODULE

  implicit none

  TYPE(GRID)                 :: gd
  character(LEN=*)           :: filetype
  character(LEN=*), optional :: message
  logical, optional          :: destroy

! Local variables

  TYPE(ATTRIBUTE), pointer :: textfile

  integer ibeg, iend
  logical if_exist
  logical remove
  integer, parameter :: lunout = 91

  IF ( my_rank == 0 .or. verbose_mpi ) THEN
  
  remove = .false.
  IF( PRESENT(destroy) ) THEN
   IF(destroy) remove = .true.
  ENDIF

  CALL GET_ATTRIBUTE(gd, filetype, textfile)

  CALL STRING_LIMITS( textfile%str, ibeg, iend )

  INQUIRE(file=textfile%str(ibeg:iend), exist=if_exist)
  
  ENDIF

  close(lunout)

  IF ( my_rank == 0 .or. verbose_mpi ) THEN
  
  IF( .NOT. if_exist .OR. remove ) THEN
    open(lunout, file=textfile%str(ibeg:iend), status='replace')
  ELSE
    open(lunout, file=textfile%str(ibeg:iend), status='old', position='append')
  ENDIF

  IF( PRESENT(message) ) CALL FILE_MESSAGE(message)
  
  ELSE
    open(lunout, file='/dev/null' )
  ENDIF

  FILE_OPEN = lunout

 RETURN
 END FUNCTION FILE_OPEN


!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
 INTEGER FUNCTION FILE_CLOSE(message)

  implicit none

  character(LEN=*), optional :: message

  integer lunout ; parameter( lunout = 91 )

  IF( PRESENT(message) ) CALL FILE_MESSAGE(message)
  
  close(lunout)
  FILE_CLOSE = 0

 RETURN
 END FUNCTION FILE_CLOSE




!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
 SUBROUTINE FILE_MESSAGE(message)

  character(LEN=*) :: message

  integer lunout ; parameter( lunout = 91 )

  write(lunout,*)
  write(lunout,"(1x,79('*'))")
  write(lunout,*)
  write(lunout,*) ' ', message
  write(lunout,*)
  write(lunout,"(1x,79('*'))")
  write(lunout,*)

  RETURN

 END SUBROUTINE FILE_MESSAGE

#endif

END MODULE FILE_MODULE
