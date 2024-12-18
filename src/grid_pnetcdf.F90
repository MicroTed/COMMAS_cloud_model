
SUBROUTINE GRID_WRITE_PNETCDF(gd, start, file, forcewrite, ncidopen)

  USE COMMASMPI_MODULE
  use mpi, only: MPI_OFFSET_KIND
!  USE PNETCDF
  implicit none
  
  integer, INTENT(IN) ::  start
  character(LEN=*), optional :: file
  integer, optional, INTENT(IN) :: forcewrite  ! flag to not overwrite "odd" times (used for enkf fcst/anal files in particular)
  TYPE(GRID)                 :: gd
  TYPE(ATTRIBUTE), pointer   :: attr
  TYPE(VARIABLE),  pointer   :: var
  integer, optional, INTENT(IN)     :: ncidopen    ! File is already open with ncid=ncidopen, so don't open again
  
  integer time, ibeg, iend, thistory
  real dt
  
  integer i
  real    ugrid, vgrid
  integer cmode
  integer n, status
  integer ncid, nt, record_id
  integer(kind=MPI_OFFSET_KIND) records ! G_NY, myOff, block_start, &
  integer(kind=MPI_OFFSET_KIND) start1d(1), count1d(1)
  integer variable_id
  integer, allocatable :: tarray(:), tharray(:), dtarray(:)
  real,    allocatable :: ugarray(:), vgarray(:)
  logical, allocatable :: index(:)
  real,    allocatable :: positions(:,:)
  character(LEN=10)    :: tname
  integer nx, ny, nz
  character(LEN=120) filename
  integer              :: checkhistory

!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and wtype

  IF ( PRESENT( forcewrite ) ) THEN
    checkhistory = Max(1, forcewrite )
  ELSE
    checkhistory = 0
  ENDIF


  IF ( .not. present( ncidopen ) ) THEN

  IF( PRESENT(file) ) THEN
  
    CALL STRING_LIMITS(file, ibeg, iend)
    filename(ibeg:iend) = file(ibeg:iend)
  
  ELSE
  
    CALL GET_ATTRIBUTE(gd, 'MEMBER_NAME', attr)
    CALL STRING_LIMITS(attr%str, ibeg, iend)
    filename(ibeg:iend) = attr%str(ibeg:iend)
  
  ENDIF

  
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  START = ', start

  

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  checkhistory = ', checkhistory
  
 IF ( parallelio ) THEN
   IF ( parallelio_type /= 2 ) THEN ! HDF5 parallel
!   cmode  = ior(NF90_WRITE,nf90_iotype)
   ELSE ! pnetcdf
    cmode  = ior( NF90_WRITE, nf90_share) ! ior(NF90_WRITE, ior(NF90_PNETCDF,NF90_MPIIO))
   ENDIF

!        status = nf90mpi_open(my_comm, filename, NF90_WRITE, my_info, ncid)
        status = nf90mpi_open(my_comm, filename, cmode, my_info, ncid)

 ENDIF
  
  IF (status /= NF90_NOERR) THEN 
    write(0,*) 'WRITE_PNETCDF:  Problem opening file'
    write(0,*) 'WRITE_PNETCDF:  Filename = ', filename(ibeg:iend)
    write(0,*) 'WRITE_PNETCDF:  MY_COMM = ', my_comm
    write(0,*) 'WRITE_PNETCDF:  MY_INFO = ', my_info
    write(0,*) 'WRITE_PNETCDF:  MY_RANK = ', my_rank
    write(0,*) 'WRITE_PNETCDF:  NCID = ', ncid 
    write(0,*) 'status = ', NF90mpi_STRERROR(status)
  ENDIF


  ELSE
    ncid = ncidopen
  ENDIF ! present(ncidopen)

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF: NCID = ', ncid  , NF90mpi_STRERROR(status)

!-----------------------------------------------------------------------------------------------
! Inquire as to how many time steps have been written to determine where to put the data%%%

  status = nf90mpi_INQUIRE(ncid, unlimitedDimID = record_id)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'GRID_WRITE_PNETCDF: ERROR FINDING TIME DIM '  , NF90mpi_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'GRID_WRITE_PNETCDF: ERROR FINDING TIME DIM '  , NF90mpi_STRERROR(status)
#endif

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  record_id = ', record_id

  status = nf90mpi_INQUIRE_DIMENSION(ncid, record_id, len=records)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'GRID_WRITE_PNETCDF: ERROR READING TIME DIM ' !,NF90_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'GRID_WRITE_PNETCDF: ERROR READING TIME DIM '  , NF90mpi_STRERROR(status)
#endif  

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  RECORDS = ', records

!-----------------------------------------------------------------------------------------------
! Get model time array
  
    CALL GET_VARIABLE(gd, 'TIME', time) ! note that if the "+1" time was read, then TIME got set to that. So will write back to the same time.
    CALL GET_VARIABLE(gd, 'THISTORY', thistory)
    CALL GET_VARIABLE(gd, 'UGRID', ugrid)
    CALL GET_VARIABLE(gd, 'VGRID', vgrid)
    CALL GET_VARIABLE(gd, 'DT', dt)

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  TIME = ', time, thistory, ugrid, vgrid, dt

!-----------------------------------------------------------------------------------------------
! IF RECORDS > 0 %%%% determine NT from input time value, ELSE  NT = 1

  IF( records == 0 .or. (writing_restart .and. restart_separate) ) THEN
  
    nt = 1
  
  ELSE
  
! Allocate memory for time and index information
  
    allocate(tarray(records))
    allocate(tharray(records))
    allocate(ugarray(records))
    allocate(vgarray(records))
    allocate(dtarray(records))
    allocate(index(records))

! Read in time array information

! In parallel IO mode, note that all threads get the same values, so just use the default collective access.

   ! switch to independent mode so that rank 0 can write these by itself
    status = NF90MPI_BEGIN_INDEP_DATA(NCID)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error setting independent mode'


    status = nf90mpi_INQ_VARID(ncid,'TIME',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE COORINDATE TIME '  , NF90mpi_STRERROR(status)
    status = nf90mpi_GET_VAR(ncid, variable_id, tarray)
    IF(status /= NF90_NOERR) THEN
    write(0,*) 'WRITE_PNETCDF:  ERROR READING COORDINATE TIME '  , NF90mpi_STRERROR(status)
    write(0,*) 'WRITE_PNETCDF: ncid, variable_id, tarray',ncid, variable_id, records, tarray(1:records)
    ENDIF

    status = nf90mpi_INQ_VARID(ncid,'THISTORY',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE COORINDATE THISTORY '  , NF90mpi_STRERROR(status)
    status = nf90mpi_GET_VAR(ncid, variable_id, tharray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR READING COORDINATE THISTORY '  , NF90mpi_STRERROR(status)

    status = nf90mpi_INQ_VARID(ncid,'UGRID',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE UGRID '  , NF90mpi_STRERROR(status)
    status = nf90mpi_GET_VAR(ncid, variable_id, ugarray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR READING UGRID '  , NF90mpi_STRERROR(status)

    status = nf90mpi_INQ_VARID(ncid,'VGRID',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE VGRID '  , NF90mpi_STRERROR(status)
    status = nf90mpi_GET_VAR(ncid, variable_id, vgarray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR READING VGRID '  , NF90mpi_STRERROR(status)

    status = nf90mpi_INQ_VARID(ncid,'DT',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE DT '  , NF90mpi_STRERROR(status)
    status = nf90mpi_GET_VAR(ncid, variable_id, dtarray)
   IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR READING DT '  , NF90mpi_STRERROR(status)
  
! Search to find the correct time in the netcdf file
  
    index(:) = any(tarray(:) .eq. time) 
  
    IF( count(index) .EQ. 0 .AND. time .LE. tarray(records) ) THEN

#ifdef MPI
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*)
      write(0,*) 'GRID_WRITE_PNETCDF: ERROR, CANNOT FIND A MATCHING TIME '
      write(0,*)
  
      DO n = 1,records
        write(0,*) 'GRID_WRITE_PNETCDF: INPUT TIME ', time, ' DATA TIMES: ', tarray(n), index(n)
      ENDDO
  
      write(0,*)
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*)
#else
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*)
      write(6,*) 'GRID_WRITE_PNETCDF: ERROR, CANNOT FIND A MATCHING TIME '
      write(6,*)
  
      DO n = 1,records
        WRITE(6,*) 'GRID_WRITE_PNETCDF: INPUT TIME ', time, ' DATA TIMES: ', tarray(n), index(n)
      ENDDO
  
      write(6,*)
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*)
#endif      
      DO nt = 1,records ! records,1,-1 ! 1,records
        IF( tarray(nt) .ge. time .or.   &
           (nt .gt. 1 .and. (tarray(nt) - tarray(nt-1)) .ne. tharray(nt) .and. tarray(nt) .gt. start )) EXIT
      ENDDO
      nt = Min(records,nt)
#ifdef MPI
      write(0,*) 'GRID_WRITE_PNETCDF: Overwriting time ',tarray(nt),index(nt)
#else
      write(6,*) 'GRID_WRITE_PNETCDF: Overwriting time ',tarray(nt),index(nt)
#endif  
    ELSE
  
! IF time is in the array or time > tarray(records), use a do loop to find index (SEE NOTE BELOW).
  
      DO nt = 1,records
        IF( tarray(nt) .eq. time ) EXIT
      ENDDO

! Since the DO LOOP exits at NT+1, dont need to add anything to NT to get correct index append data to end of arrays..

! Check if the last time was an 'odd' time (i.e. not a history dump time), unless forcewrite is present
!
! 11/7/2008 modified by erm to also check whether last dump was a valid history time (using the Mod function) for
!           the case of restarting at an 'odd' time into a new output file and saving the initial 'odd' time.  In 
!           this case the second time kept getting overwritten because the first interval is always .ne. thistory.
! 
!  If so, want to overwrite that time unless grid motion has changed or dt has changed

     IF ( nt .eq. records + 1 .and. records .gt. 1 .and. checkhistory .eq. 0 ) THEN

      IF ( (tarray(records) - tarray(records-1)) .ne. tharray(records) .and.   &
             Mod (tarray(records), tharray(records) ) .ne. 0  ) THEN

       IF ( ugrid .eq. ugarray(records) .and. vgrid .eq. vgarray(records) .and. dt .eq. dtarray(records) ) THEN

        IF( DEBUG_IO ) THEN
         print *,'GRID_WRITE_PNETCDF: ugrid,vgrid,dt = ',ugrid,ugarray(records),vgrid,vgarray(records),dt,dtarray(records)
        ENDIF

        nt = records

       ENDIF
      ENDIF
!     ELSEIF ( nt .eq. records + 1 .and. records .eq. 1 .and. checkhistory .eq. 0 ) THEN
!       IF  ( Mod(tharray(records), tarray(records) 
!        nt = records
!       ENDIF
     ENDIF
  
    ENDIF
    
#ifdef MPI
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  INPUT TIME ', time
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  TIME ARRAY(NT) ', tarray(Min(nt,records))
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  INDEX ARRAY(:) ', index
#else
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF:  INPUT TIME ', time
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF:  TIME ARRAY(NT) ', tarray(Min(nt,records))
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF:  INDEX ARRAY(:) ', index
#endif

    deallocate(tarray)
    deallocate(index)
    deallocate(tharray)
    deallocate(ugarray)
    deallocate(vgarray)
    deallocate(dtarray)

    CALL MPI_BARRIER(my_comm, mpi_error_code)
     status = NF90MPI_END_INDEP_DATA(NCID)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error ending independent mode'
  
  ENDIF

#ifdef MPI
  IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF:  NT = ', nt
#else
  IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF:  NT = ', nt
#endif

!-----------------------------------------------------------------------------------------------
! Add to time array

  IF( DEBUG_IO ) write(6,*) 'WRITE_PNETCDF: writing time information, NT = ', nt
  
  status = nf90mpi_INQ_VARID(ncid,'TIME',variable_id)

#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR INQUIRE COORINDATE TIME '  , NF90mpi_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_PNETCDF:  ERROR INQUIRE COORINDATE TIME '  , NF90mpi_STRERROR(status)
#endif
  
  status = nf90_NOERR
  IF ( my_rank == 0 .or. historyoutput == 0 ) THEN
!  status = nf90mpi_PUT_VAR(ncid, variable_id, time, (/ nt /) )
  ENDIF

#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF:  ERROR WRITING COORDINATE TIME', NF90mpi_STRERROR(status), &
      'my_rank =',my_rank,nt,time,status
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_PNETCDF:  ERROR WRITING COORDINATE TIME'  , NF90mpi_STRERROR(status)
#endif

  IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF:  WROTE TIME INFORMATION'

!-----------------------------------------------------------------------
!  NetCDF:  Write out grid dimensions

 
!-----------------------------------------------------------------------
!  NetCDF:  Define variables - scalar arrays are written out as single xyz3d arrays

  DO n = 1,size(gd%var)
  
#ifdef MPI
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF: WRITING OUT VARIABLE:  ', gd%var(n)%name
#else
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF: WRITING OUT VARIABLE:  ', gd%var(n)%name
#endif

    IF( DEBUG_IO ) write(0,*) my_rank,'GRID_WRITE_PNETCDF: WRITING OUT VARIABLE:  ', gd%var(n)%name
    CALL GET_VARIABLE(gd, gd%var(n)%name, var)
    CALL WRITE_NCDF_VAR_PCDF(ncid,var,nt)
    IF( DEBUG_IO ) write(0,*) my_rank,'GRID_WRITE_PNETCDF: WROTE OUT VARIABLE:  ', gd%var(n)%name

    IF ( parallelio_type == 2 ) CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF

#ifdef MPI
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_PNETCDF: VARIABLE:  ', gd%var(n)%name, & 
      ' WAS WRITTEN TO THE FILE. my_rank,time=',my_rank,time
#else  
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_PNETCDF: VARIABLE:  ', gd%var(n)%name, ' WAS WRITTEN TO THE FILE'
#endif
  
  ENDDO
!    IF ( parallelio_type == 2 ) CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)
  IF ( .not. present( ncidopen ) ) THEN
   status = NF90mpi_CLOSE(ncid)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_PNETCDF: Error Closing File: '  , NF90mpi_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_PNETCDF: Error Closing File: '  , NF90mpi_STRERROR(status)
#endif
  ENDIF

  IF ( DEBUG_IO ) THEN
#ifdef MPI
  write(0,*) ''
  write(0,*) 'GRID_WRITE_PNETCDF:  FILE: ',filename(ibeg:iend),'  SUCCESSFULLY WRITTEN'
  write(0,*) ''
#else
  write(lundbg,*) ''
  write(lundbg,*) 'GRID_WRITE_PNETCDF:  FILE: ',filename(ibeg:iend),'  SUCCESSFULLY WRITTEN'
  write(lundbg,*) ''
#endif
  ENDIF

END SUBROUTINE GRID_WRITE_PNETCDF

!===============================================================================
!
!
!
!
!
!
!
!===============================================================================
SUBROUTINE WRITE_NCDF_VAR_PCDF(ncid,var,nt)

  USE COMMASMPI_MODULE
  use mpi, only: MPI_OFFSET_KIND

  implicit none

  integer ncid, nt
  TYPE(variable) :: var

  integer variable_id
!  integer, allocatable :: dims(:),count(:)
  integer(kind=MPI_OFFSET_KIND), allocatable :: dims(:),start(:),count(:)
  integer(kind=MPI_OFFSET_KIND) :: start1d(1),count1d(1),start0d(1),count0d(1)
  integer(kind=MPI_OFFSET_KIND) :: start2d(1),count2d(2)
  integer ireq, num_reqs,reqs(1),sts(1)
  integer dim, status, n1, n2, n3, stag, istag, jstag, kstag
  integer idat
  logical DEBUG_IO2
  real, allocatable :: buf3d(:,:,:)
  
  status = NF90mpi_INQ_VARID(ncid,var%name,variable_id)
  IF(status /= NF90_NOERR) THEN
    print *, 'WRITE_NETCDF_VAR:  ERROR INQUIRE VARIABLE ',var%name  , NF90mpi_STRERROR(status)
    RETURN
  ENDIF


!-------------------------------------------------------------------------------
! DEBUG STATEMENTS

  IF( DEBUG_IO ) THEN
    print *, var%name, var%type, var%istag, var%jstag, var%kstag 
  ENDIF

  DEBUG_IO2 = DEBUG_IO

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, write out the scalar

  dim = var%dim + var%tdepend

    IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR_PCDF: rank, VARIABLE, dim:  ',my_rank,var%name,var%type,dim

  IF( dim == 0 ) THEN

    ! switch to independent mode so that rank 0 can write these by itself
    status = NF90MPI_BEGIN_INDEP_DATA(NCID)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error setting independent mode'
    
    IF ( my_rank == 0 ) THEN
    start0d = 1
    count0d = 1
    
    IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR_PCDF: WRITING OUT VARIABLE:  ',var%name,var%type

      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR_PCDF: before nf90mpi_put_var VARIABLE:  ',var%name,' = ',var%flt
    IF( var%type(1:5) == 'icnst' ) THEN
      status = nf90mpi_put_var(ncid, variable_id, var%int)
      IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR_PCDF: WRITING OUT VARIABLE:  ',var%name,' = ',var%int
    ELSE
      status = nf90mpi_put_var(ncid, variable_id, var%flt) ! , start0d, count0d)
      IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR_PCDF: WRITING OUT VARIABLE:  ',var%name,' = ',var%flt
    ENDIF
      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR_PCDF: after nf90mpi_put_var VARIABLE:  ',var%name,' = ',var%flt

    IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)

    ENDIF ! my_rank == 0
    
    CALL MPI_BARRIER(my_comm, mpi_error_code)
     status = NF90MPI_END_INDEP_DATA(NCID)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error ending independent mode'


!-------------------------------------------------------------------------------
! Else we are writing out a 1, 2, 3D array 

  ELSE

    allocate(dims(dim))
    allocate(count(dim))

    dims(:) = 1

    IF( var%tdepend /= 0 ) dims(dim) = nt

    IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR_PCDF:  nt,dim,rank = ',nt,dim,my_rank


    IF( var%type(1:5) == 'rcnst' .or. var%type(1:5) == 'icnst'  ) THEN ! {

      ! switch to independent mode so that rank 0 can write these by itself
      status = NF90MPI_BEGIN_INDEP_DATA(NCID)
      IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error setting independent mode'
      
      IF ( my_rank == 0 ) THEN
        IF( var%type(1:5) == 'icnst' ) THEN
          status = nf90mpi_put_var(ncid, variable_id, var%int, start=dims)
          IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR_PCDF: WRITING OUT VARIABLE:  ',var%name,' = ',var%int,my_rank
        ELSE
          status = nf90mpi_put_var(ncid, variable_id, var%flt, start=dims)
          IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR_PCDF: WRITING OUT VARIABLE:  ',var%name,' = ',var%flt,my_rank
        ENDIF
        IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)
      
      ENDIF ! my_rank == 0
    
    CALL MPI_BARRIER(my_comm, mpi_error_code)
    status = NF90MPI_END_INDEP_DATA(NCID)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error ending independent mode'

      IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR_PCDF:  ',var%name, n1
    
    ENDIF ! }



    IF( var%type(1:3) == 'x1d' .or. var%type(1:3) == 'y1d' .or. var%type(1:3) == 'z1d' ) THEN 

      stag = -1
      IF( var%type(1:3) == 'x1d' .and. var%istag == 1) stag = 0
      IF( var%type(1:3) == 'y1d' .and. var%jstag == 1) stag = 0
      IF( var%type(1:3) == 'z1d' .and. var%kstag == 1) stag = 0

#ifdef MPI
    IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 .or. parallelio ) THEN ! single file
      IF( var%type(1:3) == 'x1d' ) dims(1) = ixbeg
      IF( var%type(1:3) == 'y1d' ) dims(1) = jybeg
      IF( var%type(1:3) == 'z1d' ) dims(1) = kzbeg
    ENDIF
#endif

      n1 = size(var%flt1d) + stag - 2*var%ng
#ifdef MPI
      IF( var%type(1:3) == 'x1d' .and. var%istag == 0 .and. ixend /= nxend) n1 = n1 + 1
      IF( var%type(1:3) == 'y1d' .and. var%jstag == 0 .and. jyend /= nyend) n1 = n1 + 1
      IF( var%type(1:3) == 'z1d' .and. var%kstag == 0 .and. kzend /= nzend) n1 = n1 + 1
#endif
      
      
! here just write from a minimum of processors that have needed data. So for z1d only write
! from the southwest corner processor(s). x1d only write from the first row (myprocj==1), etc.
! Should still allow for vertical tiling
      status = NF90MPI_BEGIN_INDEP_DATA(NCID)
      IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error setting independent mode'

      IF ( var%type(1:3) == 'z1d' .and. myproci == 1 .and. myprocj == 1 ) THEN
        status = nf90mpi_put_var(ncid, variable_id, var%flt1d(1:n1),  start=dims)
      ELSEIF ( var%type(1:3) == 'x1d' .and. myprock == 1 .and. myprocj == 1 ) THEN
        status = nf90mpi_put_var(ncid, variable_id, var%flt1d(1:n1),  start=dims)
      ELSEIF ( var%type(1:3) == 'y1d' .and. myprock == 1 .and. myproci == 1 ) THEN
        status = nf90mpi_put_var(ncid, variable_id, var%flt1d(1:n1),  start=dims)
      ENDIF

       CALL MPI_BARRIER(my_comm, mpi_error_code)
       status = NF90MPI_END_INDEP_DATA(NCID)

 
!      IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)
      IF(status /= NF90_NOERR) THEN
        print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)
#ifdef MPI
        print *,'ixbeg,ixend,nxbeg,nxend = ',ixbeg,ixend,nxbeg,nxend
#endif
        print *, 'dims = ',dims, ' stag = ',stag
      ENDIF

      IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR_PCDF:  ',var%name, n1

    ENDIF


    IF( var%type(1:3) == 'x2d' .or. var%type(1:3) == 'y2d' .or. var%type(1:3) == 'z2d' ) THEN 
! NOTE: these are special arrays for whats-it visualization and output NOT TESTED
      stag = -1
      IF( var%type(1:3) == 'x2d' .and. var%istag == 1) stag = 0
      IF( var%type(1:3) == 'y2d' .and. var%jstag == 1) stag = 0
      IF( var%type(1:3) == 'z2d' .and. var%kstag == 1) stag = 0

#ifdef MPI
    IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 .or. parallelio ) THEN ! single file
      IF( var%type(1:3) == 'x2d' ) THEN 
        dims(1) = 1
        dims(2) = ixbeg
      ENDIF
      IF( var%type(1:3) == 'y2d' ) THEN
        dims(1) = 1
        dims(2) = jybeg
      ENDIF
      IF( var%type(1:3) == 'z2d' ) THEN
        dims(1) = 1
        dims(2) = kzbeg
      ENDIF
    ENDIF
#endif

      n1 = 3
      n2 = size(var%flt2d,dim=2) + stag - 2*var%ng

#ifdef MPI
      IF( var%type(1:3) == 'x2d' .and. var%istag == 0 .and. ixend /= nxend) n2 = n2 + 1
      IF( var%type(1:3) == 'y2d' .and. var%jstag == 0 .and. jyend /= nyend) n2 = n2 + 1
      IF( var%type(1:3) == 'z2d' .and. var%kstag == 0 .and. kzend /= nzend) n2 = n2 + 1
#endif

      status = nf90mpi_put_var_all(ncid, variable_id, var%flt2d(1:n1,1:n2), dims)

      IF(status /= NF90_NOERR) THEN
        print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)
        RETURN
      ENDIF
      IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR_PCDF:  ',var%name, n1, n2

    ENDIF


    IF( var%type(1:4) == 'xy2d' .or. var%type(1:4) == 'xz2d' .or. var%type(1:4) == 'yz2d' ) THEN 

      istag = -1
      jstag = -1
      IF( var%type(1:4) == 'xy2d' .and. var%istag == 1) istag = 0
      IF( var%type(1:4) == 'xy2d' .and. var%jstag == 1) jstag = 0
      IF( var%type(1:4) == 'xz2d' .and. var%istag == 1) istag = 0
      IF( var%type(1:4) == 'xz2d' .and. var%kstag == 1) jstag = 0
      IF( var%type(1:4) == 'yz2d' .and. var%jstag == 1) istag = 0
      IF( var%type(1:4) == 'yz2d' .and. var%kstag == 1) jstag = 0

#ifdef MPI
    IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 .or. parallelio ) THEN ! single file
      IF( var%type(1:4) == 'xy2d' ) THEN 
        dims(1) = ixbeg
        dims(2) = jybeg
      ENDIF
      IF( var%type(1:4) == 'xz2d' ) THEN
        dims(1) = ixbeg
        dims(2) = kzbeg
      ENDIF
      IF( var%type(1:4) == 'yz2d' ) THEN
        dims(1) = jybeg
        dims(2) = kzbeg
      ENDIF
    ENDIF
#endif

      n1 = size(var%flt2d,dim=1) + istag - 2*var%ng
      n2 = size(var%flt2d,dim=2) + jstag - 2*var%ng

#ifdef MPI
      IF( var%type(1:4) == 'xy2d' ) THEN 
        IF( var%istag == 0 .and. ixend /= nxend ) THEN
         n1 = n1 + 1
        ENDIF
        IF( var%jstag == 0 .and. jyend /= nyend ) THEN
         n2 = n2 + 1
        ENDIF   
      ENDIF
      IF( var%type(1:4) == 'xz2d' .and. var%jstag == 0) THEN
        IF( var%istag == 0 .and. ixend /= nxend ) THEN
         n1 = n1 + 1
        ENDIF
        IF( var%kstag == 0 .and. kzend /= nzend ) THEN
         n2 = n2 + 1
        ENDIF   
      ENDIF
      IF( var%type(1:4) == 'yz2d' .and. var%kstag == 0) THEN
        IF( var%jstag == 0 .and. jyend /= nyend ) THEN
         n1 = n1 + 1
        ENDIF
        IF( var%kstag == 0 .and. kzend /= nzend ) THEN
         n2 = n2 + 1
        ENDIF   
      ENDIF
#endif

      IF ( .not. ( myprock > 1 .and. var%type(1:4) == 'xy2d' )) THEN
      status = nf90mpi_put_var_all(ncid, variable_id, var%flt2d(1:n1,1:n2), start=dims)

      IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)

      IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR_PCDF:  ',var%name, n1, n2
      
      ENDIF

    ENDIF


    IF( var%type(1:4) == 'xyz3' ) THEN

      istag = -1
      jstag = -1
      kstag = -1
      IF( var%istag == 1) istag = 0
      IF( var%jstag == 1) jstag = 0
      IF( var%kstag == 1) kstag = 0
      n1 = size(var%flt3d,dim=1) + istag - 2*var%ng
      n2 = size(var%flt3d,dim=2) + jstag - 2*var%ng
      n3 = size(var%flt3d,dim=3) + kstag - 2*var%ng
      
!      write(0,*) 'WRITE_NCDF_VAR_PCDF: n1,n2,n3 = ',n1,n2,n3, var%name

#ifdef MPI
    ! dims ("start" array) already initialized as 1,1,1; So only change if this is an MPI job
!mpi! #ifdef NC4
    IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3  .or. parallelio ) THEN ! single file
      dims(1) = ixbeg
      dims(2) = jybeg
      dims(3) = kzbeg
    ENDIF
!mpi! #endif
      IF( var%istag == 0 .and. ixend /= nxend ) THEN
       n1 = n1 + 1
      ENDIF
      IF( var%jstag == 0 .and. jyend /= nyend ) THEN
       n2 = n2 + 1
      ENDIF 
      IF( var%kstag == 0 .and. kzend /= nzend ) THEN
       n3 = n3 + 1
      ENDIF
#endif
      
      allocate(buf3d(n1,n2,n3))
      buf3d(1:n1,1:n2,1:n3) = var%flt3d(1:n1,1:n2,1:n3)
!      status = nf90mpi_put_var_all(ncid, variable_id, var%flt3d(1:n1,1:n2,1:n3), start=dims)
      status = nf90mpi_put_var_all(ncid, variable_id, buf3d, start=dims)

      ireq = 0
!      status = nf90mpi_iput_var(ncid, variable_id, buf3d, ireq, start=dims)
      IF(status /= NF90_NOERR) write(0,*)'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)

 !      num_reqs = 1
 !      reqs(1) = ireq
 !      status = nf90mpi_wait_all(ncid, num_reqs, reqs, sts)
       deallocate( buf3d )
      IF(status /= NF90_NOERR) write(0,*)'WRITE_NCDF_VAR_PCDF:  Error nf90mpi_wait_all: ', var%name  , NF90mpi_STRERROR(status)
    ENDIF

    deallocate(dims)
    deallocate(count)

  ENDIF


  IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR_PCDF:  Error writing variable: ', var%name  , NF90mpi_STRERROR(status)

END SUBROUTINE WRITE_NCDF_VAR_PCDF


! ################################################
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
! ################################################
  subroutine handle_err(errcode)
!    use netcdf
    implicit none
    integer, intent(in) :: errcode
    
    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90mpi_strerror(errcode))
       CALL commasmpi_abort()
       stop "Stopped"
    endif
  end subroutine handle_err

