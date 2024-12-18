!===============================================================================
! 
!
!  
!  
!
!  
!  
!===============================================================================

SUBROUTINE GRID_DEFINE_NETCDF(gd, file)

!  USE GRID_MODULE
#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  character(LEN=*), optional :: file
  TYPE(GRID)                 :: gd
  TYPE(ATTRIBUTE), pointer   :: attr

  integer n, status, ibeg, iend
  integer ncid, dim_id(0:8),dimvals(0:8)
  integer cmode
  character(LEN=120) filename
  
  integer mpistag
  integer variable_id
  
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and wtype

  IF( PRESENT(file) ) THEN
 
   CALL STRING_LIMITS(file, ibeg, iend)
   filename(ibeg:iend) = file(ibeg:iend)
  
  ELSE
  
   CALL GET_ATTRIBUTE(gd, 'MEMBER_NAME', attr)
   CALL STRING_LIMITS(attr%str, ibeg, iend)
   filename(ibeg:iend) = attr%str(ibeg:iend)
  
  ENDIF

#ifdef NC4
! linking with netcdf4
! write(0, *) 'netcdfversion = ',netcdfversion
 IF ( netcdfversion .eq. 4 .or. parallelio ) THEN ! turn on hdf5 format to allow compression, but keep classic model
  IF ( parallelio ) THEN
#ifdef MPI
   IF ( parallelio_type == 1 ) THEN ! HDF5 parallel
     cmode  = ior(NF90_NETCDF4,nf90_iotype)
    cmode = IOR(nf90_netcdf4, nf90_classic_model) 
    cmode = IOR(cmode, nf90_mpiio) 
   ELSEIF ( parallelio_type == 2 ) THEN ! pnetcdf
     cmode  = ior(ior(pnetcdf_flag,NF90_MPIIO),nf90_64bit_offset)
!       cmode = ior(IOR(nf90_CLOBBER, nf90_PNETCDF)
   ELSE
     write(0,*) 'GRID_DEFINE_NETCDF: NON-MPI Compile! Please set the value of parallelio to FALSE'
     CALL commasmpi_abort()
   ENDIF
   IF( parallelio_type == 1 ) THEN
     IF ( DEBUG_IO )write(0,*) 'hdf5 create, rank = ',my_rank
     status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
   ELSE
    ! write(0,*) 'pnetcdf create, rank = ',my_rank
    status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
!    status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
   ENDIF
   IF (status /= NF90_NOERR) THEN
    write(0,*) 'GRID_DEFINE_NETCDF:  NETCDF FILE CREATE FAILED FOR PARALLEL IO', &
    NF90_STRERROR(status),'MY_RANK=',my_rank
    
    ! try again with _par version
!    cmode  = ior(NF90_PNETCDF,NF90_MPIIO)
!    write(0,*) 'pnetcdf create again, rank = ',my_rank
!    status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
    !status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
    
   ENDIF

   if (status /= nf90_noerr) call handle_err(status)
#else
   write(0,*) 'GRID_DEFINE_NETCDF: NON-MPI Compile! Please set the value of parallelio to FALSE'
   STOP
#endif
  ELSE
!#else
   IF ( deflate .ne. 0 .and. deflate_level .ge. 1 ) THEN ! use netcdf4
   cmode = ior(NF90_CLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
!   cmode = ior(NF90_CLOBBER,NF90_NETCDF4)
   ELSE ! if no compression, then use the classic model
   cmode = ior(NF90_CLOBBER,NF90_CLASSIC_MODEL)
   ENDIF
!   write(0,*) 'cmode = ',cmode,NF90_CLOBBER,NF90_CLASSIC_MODEL,NF90_HDF5,NF90_FORMAT_NETCDF4_CLASSIC
   status = NF90_CREATE(filename(ibeg:iend),cmode,ncid)
   if (status /= nf90_noerr) call handle_err(status)
  ENDIF
!#endif
 !  status = NF90_CREATE(filename(ibeg:iend),ior(ior(NF90_CLOBBER,NF90_HDF5),NF90_CLASSIC_MODEL),ncid)
 ELSE ! allow version 3 as default for now until utilities all work with netcdf4/hdf5
  status = NF90_CREATE(filename(ibeg:iend),NF90_CLOBBER,ncid)
 ENDIF
#else
! linking with netcdf3:
  status = NF90_CREATE(filename(ibeg:iend),NF90_CLOBBER,ncid)
#endif
! Use this version if the record size will exceed 2GB:
!
!  status = NF90_CREATE(filename(ibeg:iend),or(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)

#ifdef MPI
  IF (status /= NF90_NOERR) THEN 
   write(0,*) 'GRID_DEFINE_NETCDF:  Problem opening file'
   write(0,*) 'GRID_DEFINE_NETCDF:  Filename = ', filename(ibeg:iend)
   write(0,*) 'GRID_DEFINE_NETCDF:  MY_COMM = ', my_comm
   write(0,*) 'GRID_DEFINE_NETCDF:  MY_INFO = ', my_info
   write(0,*) 'GRID_DEFINE_NETCDF:  MY_RANK = ', my_rank
   write(0,*) 'GRID_DEFINE_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
  ENDIF
#else
  IF (status /= NF90_NOERR) THEN 
   write(6,*) 'GRID_DEFINE_NETCDF:  Problem opening file'
   write(6,*) 'GRID_DEFINE_NETCDF:  Filename = ', filename(ibeg:iend)
   write(6,*) 'GRID_DEFINE_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
  ENDIF
#endif
   if (status /= nf90_noerr) call handle_err(status)
  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFINE_NETCDF: FILE CREATED'

!-----------------------------------------------------------------------
!  NetCDF:  Define dimensions

  status     = NF90_DEF_DIM(ncid, 'TIME', NF90_UNLIMITED, dim_id(0))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error unlimited dimension'

#ifdef MPI

 IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 .or. parallelio ) THEN ! creating one file with round-robin or parallel access
  status     = NF90_DEF_DIM(ncid, 'XC', nxend-1, dim_id(1))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nxend-1'
  status     = NF90_DEF_DIM(ncid, 'XE', nxend,   dim_id(2))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nxend'
  
  status     = NF90_DEF_DIM(ncid, 'YC', nyend-1, dim_id(3))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nyend-1'
  status     = NF90_DEF_DIM(ncid, 'YE', nyend,   dim_id(4))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nyend'
  
  status     = NF90_DEF_DIM(ncid, 'ZC', nzend-1, dim_id(5))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nzend-1'
  status     = NF90_DEF_DIM(ncid, 'ZE', nzend,   dim_id(6))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nzend'

   IF ( .false. .and. parallelio ) THEN
    dimvals(1) = nxend-1
    dimvals(2) = nxend
    dimvals(3) = nyend-1
    dimvals(4) = nyend
    dimvals(5) = nzend-1
    dimvals(6) = nzend

  ELSE

  IF ( nproci > 1 ) THEN
  dimvals(1) = itile
  ELSE
  dimvals(1) = itile-1
  ENDIF
  dimvals(2) = itile

  IF ( nprocj > 1 ) THEN
  dimvals(3) = jtile
  ELSE
  dimvals(3) = jtile-1
  ENDIF
  dimvals(4) = jtile

  IF ( nprock > 1 ) THEN
  dimvals(5) = ktile
  ELSE
  dimvals(5) = ktile-1
  ENDIF
  dimvals(6) = ktile
  
  ENDIF
  
 ELSE ! each tile is writing its own file (better for lots of tiles)
 
  mpistag = 1

  if (ixend /= nxend) mpistag = 0

  CALL GET_ATTRIBUTE(gd, 'NX', attr)
  status     = NF90_DEF_DIM(ncid, 'XC', attr%int-mpistag, dim_id(1))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nx-1'
  status     = NF90_DEF_DIM(ncid, 'XE', attr%int,   dim_id(2))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nx'

  dimvals(1) = attr%int-mpistag
  dimvals(2) = attr%int

  mpistag = 1
  if (jyend /= nyend) mpistag = 0

  
  CALL GET_ATTRIBUTE(gd, 'NY', attr)
  status     = NF90_DEF_DIM(ncid, 'YC', attr%int-mpistag, dim_id(3))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining ny-1'
  status     = NF90_DEF_DIM(ncid, 'YE', attr%int,   dim_id(4))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining ny'
  
  dimvals(3) = attr%int-mpistag
  dimvals(4) = attr%int

  mpistag = 1
  if (kzend /= nzend) mpistag = 0
  
  CALL GET_ATTRIBUTE(gd, 'NZ', attr)
  status     = NF90_DEF_DIM(ncid, 'ZC', attr%int-mpistag, dim_id(5))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nz-1'
  status     = NF90_DEF_DIM(ncid, 'ZE', attr%int,   dim_id(6))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nz'

  dimvals(5) = attr%int-mpistag
  dimvals(6) = attr%int

 ENDIF
 
#else

  mpistag = 1

!mpi! #ifdef MPI
!mpi!   if (ixend /= nxend) mpistag = 0  
!mpi! #endif

  CALL GET_ATTRIBUTE(gd, 'NX', attr)
  status     = NF90_DEF_DIM(ncid, 'XC', attr%int-mpistag, dim_id(1))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nx-1'
  status     = NF90_DEF_DIM(ncid, 'XE', attr%int,   dim_id(2))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nx'

  dimvals(1) = attr%int-mpistag
  dimvals(2) = attr%int

!mpi! #ifdef MPI
!mpi!   mpistag = 1
!mpi!   if (jyend /= nyend) mpistag = 0  
!mpi! #endif
  
  CALL GET_ATTRIBUTE(gd, 'NY', attr)
  status     = NF90_DEF_DIM(ncid, 'YC', attr%int-mpistag, dim_id(3))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining ny-1'
  status     = NF90_DEF_DIM(ncid, 'YE', attr%int,   dim_id(4))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining ny'
  
  dimvals(3) = attr%int-mpistag
  dimvals(4) = attr%int
  
!mpi! #ifdef MPI
!mpi!   mpistag = 1
!mpi!   if (kzend /= nzend) mpistag = 0  
!mpi! #endif  
  
  CALL GET_ATTRIBUTE(gd, 'NZ', attr)
  status     = NF90_DEF_DIM(ncid, 'ZC', attr%int-mpistag, dim_id(5))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nz-1'
  status     = NF90_DEF_DIM(ncid, 'ZE', attr%int,   dim_id(6))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining nz'

  dimvals(5) = attr%int-mpistag
  dimvals(6) = attr%int
#endif

#ifdef OPENDX
  status     = NF90_DEF_DIM(ncid, 'dim2', 2, dim_id(7))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining dim2'
  status     = NF90_DEF_DIM(ncid, 'dim3', 3, dim_id(8))
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error defining dim3'
#endif

#ifdef MPI
    IF( DEBUG_IO ) write(0,*) 'GRID_DEFINE_NETCDF: BEFORE DEFINING VARIABLES '
 !   CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
 !   status = NF90_ENDDEF(NCID)
 ! IF(status /= NF90_NOERR) write(0,*)'DEFINE_NETCDF:  Error ending define mode 2,rank', status, NF90_STRERROR(status), my_rank
 !   status = NF90_REDEF(NCID)
 ! IF(status /= NF90_NOERR) write(0,*)'DEFINE_NETCDF:  Error entering redefine mode 2,rank', status, NF90_STRERROR(status), my_rank
#endif

  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFINE_NETCDF: DIMENSIONS ARE DEFINED'

!-----------------------------------------------------------------------
!  NetCDF:  Define variables - scalar arrays are written out as single xyz3d arrays

!  IF ( my_rank == 0 ) write(0,*) 'GRID_DEFINE_NETCDF: size(gd%var) = ',size(gd%var)
  DO n = 1,size(gd%var)
  
  IF( DEBUG_IO ) write(0,*) 'GRID_DEFINE_NETCDF: DEFINING VARIABLE ', gd%var(n)%name,n,size(gd%var), my_rank
  
#ifdef MPI
!    CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
#endif
  
  CALL DEFINE_NCDF_VAR(ncid,dim_id,gd%var(n),dimvals)
  
#ifdef MPI
    IF( DEBUG_IO .and. .false.) THEN
    write(0,*) 'GRID_DEFINE_NETCDF: AFTER DEFINING VARIABLE ', gd%var(n)%name, my_rank
    CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
  status = NF90_ENDDEF(NCID)
  IF(status /= NF90_NOERR) THEN
   write(0,*) 'DEFINE_NETCDF:  Error ending define mode 1,rank', status, NF90_STRERROR(status), my_rank
   call commasmpi_abort()
  ENDIF
    status = NF90_REDEF(NCID)
  IF(status /= NF90_NOERR) write(0,*) 'DEFINE_NETCDF:  Error entering redefine mode 1,rank', status, NF90_STRERROR(status), my_rank
   ENDIF
#endif

  ENDDO

#ifdef MPI
!    IF( DEBUG_IO ) write(0,*) 'GRID_DEFINE_NETCDF: AFTER DEFINING VARIABLE ', gd%var(n)%name
 IF( DEBUG_IO .and. .false. ) THEN
    CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
  status = NF90_ENDDEF(NCID)
  IF(status /= NF90_NOERR) THEN
   write(0,*) 'DEFINE_NETCDF:  Error ending define mode 1,rank', status, NF90_STRERROR(status), my_rank
   call commasmpi_abort()
  ENDIF
    status = NF90_REDEF(NCID)
  IF(status /= NF90_NOERR) write(0,*) 'DEFINE_NETCDF:  Error entering redefine mode 1,rank', status, NF90_STRERROR(status), my_rank
  ENDIF
#endif
 
  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFINE_NETCDF: VARIABLES DEFINED'

!-----------------------------------------------------------------------
! NetCDF:  Define global attributes
  
  IF ( .true. ) THEN
  IF ( my_rank == 0 .or. historyoutput == 0 .or. parallelio_type /= 3 ) THEN
  DO n = 1,size(gd%attr)
  
!   IF ( parallelio .and. parallelio_type == 2 ) THEN ! skip some values that might not be the same for all threads
   IF ( parallelio ) THEN ! skip some values that might not be the same for all threads
     IF ( gd%attr(n)%name == 'NX' .or. gd%attr(n)%name == 'NY' .or. gd%attr(n)%name == 'TILE_INDEX' .or. &
          gd%attr(n)%name == 'MEMBER_NAME' .or. gd%attr(n)%name == 'OUTPUT_FILE_NAME' ) cycle
   ENDIF
  
  IF( DEBUG_IO ) write(6,*) 'GRID_DEFINE_NETCDF: DEFINING ATTRIBUTE ', gd%attr(n)%name,gd%attr(n)%type
  
  IF( gd%attr(n)%type == 'str') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%str)
  IF( gd%attr(n)%type == 'int') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%int)
  IF( gd%attr(n)%type == 'flt') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%flt)
  
  ENDDO
  
  ENDIF
  ENDIF

  IF( DEBUG_IO ) write(6,*) 'GRID_DEFINE_NETCDF: ATTRIBUTES DEFINED'

!-----------------------------------------------------------------------------------------------
! End definition section

!    CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
  status = NF90_ENDDEF(NCID)
  IF(status /= NF90_NOERR) print *,'DEFINE_NETCDF:  Error ending define mode', status, NF90_STRERROR(status)

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)

  status=NF90_CLOSE(NCID)
  IF(status /= NF90_NOERR) print *, 'DEFINE_NETCDF: Error Closing File: ', NF90_STRERROR(status)
  if (status /= nf90_noerr) call handle_err(status)

  IF ( DEBUG_IO ) THEN
  write(lundbg,*) 
  write(lundbg,*) 'GRID_DEFINE_NETCDF:  FILE: ',filename(ibeg:iend), ' SUCCESSFULLY DEFINED'
  write(lundbg,*) 
  ENDIF

END SUBROUTINE GRID_DEFINE_NETCDF





!===============================================================================
!     
!     
!
!    
!
!   
!
!===============================================================================
SUBROUTINE GRID_INFO_NETCDF(file, records, microphys, tarray)

#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  character(LEN=*)  :: file
  integer           :: records
  integer, optional :: tarray(records)
  character(LEN=*)  :: microphys

  integer ncid, record_id, variable_id
  integer cmode
  integer n, status
  integer ibeg,iend
  character(LEN=120) filename


    CALL STRING_LIMITS(file, ibeg, iend)
    filename(ibeg:iend) = file(ibeg:iend)
  
!  status = NF90_OPEN(filename(ibeg:iend),NF90_NOWRITE,ncid)
!mpi! #if defined (NC4) && defined (MPI)
!mpi!   cmode  = ior(NF90_NOWRITE,nf90_iotype)
!mpi!   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
!mpi! #else
  cmode  = NF90_NOWRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
  if (status /= nf90_noerr) call handle_err(status)
!mpi! #endif

#ifdef MPI
  IF (status /= NF90_NOERR) THEN 
    write(0,*) 'GRID_INFO_NETCDF:  Problem opening file'
    write(0,*) 'GRID_INFO_NETCDF:  Filename = ', filename(ibeg:iend)
    write(0,*) 'GRID_INFO_NETCDF:  MY_COMM = ', my_comm
    write(0,*) 'GRID_INFO_NETCDF:  MY_INFO = ', my_info
    write(0,*) 'GRID_INFO_NETCDF:  MY_RANK = ', my_rank
    write(0,*) 'GRID_INFO_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
  ENDIF
#else
  IF (status /= NF90_NOERR) THEN 
    write(lundbg,*) 'GRID_INFO_NETCDF:  Problem opening file'
    write(lundbg,*) 'GRID_INFO_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'GRID_INFO_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
  ENDIF
#endif

!-----------------------------------------------------------------------------------------------
! Inquire as to how many time steps have been written to determine where to put the data%%%

  status = NF90_INQUIRE(ncid, unlimitedDimID = record_id)
  IF(status /= NF90_NOERR) write(6,*) 'GRID_INFO_NETCDF: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
  
  status = NF90_INQUIRE_DIMENSION(ncid, record_id, len = records)
  IF(status /= NF90_NOERR) write(6,*) 'GRID_INFO_NETCDF: ERROR READING TIME DIM ',NF90_STRERROR(status)

  status = NF90_GET_ATT(ncid, NF90_GLOBAL, 'MICROPHYS', microphys)

  IF (status /= NF90_NOERR) THEN
    write(lundbg,*) 'GRID_INFO_NETCDF:  Problem reading attribute'
    write(lundbg,*) 'GRID_INFO_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'GRID_INFO_NETCDF:  NAME = ', microphys
  ENDIF
   if (status /= nf90_noerr) call handle_err(status)


 
  IF ( present(tarray) ) THEN
    status = NF90_INQ_VARID(ncid,'TIME',variable_id)
    IF(status /= NF90_NOERR) print *, 'GRID_INFO_NETCDF:  ERROR INQUIRE COORINDATE TIME ',NF90_STRERROR(status)
  
    status = NF90_GET_VAR(ncid, variable_id, tarray)
    IF(status /= NF90_NOERR) print *, 'GRID_INFO_NETCDF:  ERROR READING COORDINATE TIME ',NF90_STRERROR(status)
  ENDIF
  
  status=NF90_CLOSE(ncid)
  if (status /= nf90_noerr) call handle_err(status)


END SUBROUTINE GRID_INFO_NETCDF

!===============================================================================
!     
!     
!
!    
!
!   
!
!===============================================================================
SUBROUTINE GRID_INFO2_NETCDF(ncid, records, microphys)

#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  integer           :: records
  character(LEN=*)  :: microphys

  integer ncid, record_id, variable_id
  integer cmode
  integer n, status
  integer ibeg,iend
  character(LEN=120) filename


  
  
!  status = NF90_OPEN(filename(ibeg:iend),NF90_NOWRITE,ncid)
!mpi! #if defined (NC4) && defined (MPI)
!mpi!   cmode  = ior(NF90_NOWRITE,nf90_iotype)
!mpi!   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
!mpi! #else
!mpi! #endif


!-----------------------------------------------------------------------------------------------
! Inquire as to how many time steps have been written to determine where to put the data%%%

  status = NF90_INQUIRE(ncid, unlimitedDimID = record_id)
  IF(status /= NF90_NOERR) write(6,*) 'GRID_INFO_NETCDF: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
  
  status = NF90_INQUIRE_DIMENSION(ncid, record_id, len = records)
  IF(status /= NF90_NOERR) write(6,*) 'GRID_INFO_NETCDF: ERROR READING TIME DIM ',NF90_STRERROR(status)

  status = NF90_GET_ATT(ncid, NF90_GLOBAL, 'MICROPHYS', microphys)

  IF (status /= NF90_NOERR) THEN
    write(lundbg,*) 'GRID_INFO_NETCDF:  Problem reading attribute'
!    write(lundbg,*) 'GRID_INFO_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'GRID_INFO_NETCDF:  NAME = ', microphys
  ENDIF
   if (status /= nf90_noerr) call handle_err(status)


  


END SUBROUTINE GRID_INFO2_NETCDF



SUBROUTINE GRID_OPEN_NETCDF(gd,file,cmode,ncid,nx_or_nxend)

#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  TYPE(GRID)                 :: gd
  character(LEN=*)  :: file

  integer ncid, record_id, variable_id
  integer cmode
  logical, optional, INTENT(IN)     :: nx_or_nxend ! tells whether to use nx or nxend for allocation
  integer :: ipconc, ichaff, inucopt, ienkf, isfcphys
  integer n, status
  integer ibeg,iend
  character(LEN=80) filename
  character(LEN=desc_length) microphys

  integer, save :: ensize = 1
  integer, save ::  nx, ny, nz
  integer XEDimID, YEDimID, ZEDimID

    CALL STRING_LIMITS(file, ibeg, iend)
    filename(ibeg:iend) = file(ibeg:iend)
  
!  status = NF90_OPEN(filename(ibeg:iend),NF90_NOWRITE,ncid)
!mpi! #if defined (NC4) && defined (MPI)
!mpi!   cmode  = ior(NF90_NOWRITE,nf90_iotype)
!mpi!   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
!mpi! #else
!  cmode  = NF90_NOWRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
  if (status /= nf90_noerr) call handle_err(status)
!mpi! #endif

!-----------------------------------------------------------------------
! READ IN THE ATTRIBUTES NEEDED TO DEFINE THE GRID NAMESPACE

  status = NF90_GET_ATT(ncid, NF90_GLOBAL, 'MICROPHYS', microphys)

  IF (status /= NF90_NOERR) THEN
    write(lundbg,*) 'GRID_READ_NETCDF:  Problem reading attribute'
    write(lundbg,*) 'GRID_READ_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'GRID_READ_NETCDF:  NAME = ', microphys
  ENDIF

  CALL STRING_LIMITS(microphys, ibeg, iend)

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  MICRO = ', microphys(ibeg:iend)

!
! Get value of ienkf
!
    status = NF90_INQ_VARID(ncid,'IENKF',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE IENKF ',NF90_STRERROR(status)
     ienkf = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, ienkf)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN IENKF:  ',ienkf
    ENDIF

!
! Get value of isfcphys
!
    status = NF90_INQ_VARID(ncid,'ISFCPHYS',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE isfcphys ',NF90_STRERROR(status)
     isfcphys = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, isfcphys)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN isfcphys:  ',isfcphys
    ENDIF


!
! Get value of ichaff
!
    status = NF90_INQ_VARID(ncid,'ICHAFF',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE ICHAFF ',NF90_STRERROR(status)
     ichaff = -1
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, ichaff)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN ICHAFF:  ',ichaff
    ENDIF

!
! Get value of inucopt
!
    status = NF90_INQ_VARID(ncid,'INUCOPT',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE isfcphys ',NF90_STRERROR(status)
     inucopt = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, inucopt)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN inucopt:  ',inucopt
    ENDIF

!
! Get value of ipconc
!
    status = NF90_INQ_VARID(ncid,'IPCONC',variable_id)
    IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR INQUIRE IPCONC ',NF90_STRERROR(status)

    status = NF90_GET_VAR(ncid, variable_id, ipconc)
    IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN IPCONC:  ',ipconc

  IF ( microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) .eq. 'ZVD' ) THEN

! check for special cases with diagnosed CCW and/or CCI
! ipconc 
    IF ( ipconc .eq. 1  .or. ipconc .eq. 0 ) THEN
    status = NF90_INQ_VARID(ncid,'CCW',variable_id)
      IF(status == NF90_NOERR) THEN
        IF ( ipconc .eq. 1 ) THEN
          ipconc = -3
        ELSE
          ipconc = -2
        ENDIF
      ENDIF
 
      IF ( ipconc .eq. 0 ) THEN
        status = NF90_INQ_VARID(ncid,'CCI',variable_id)
        IF(status == NF90_NOERR) ipconc = -1
      ENDIF
   ENDIF
  
  ENDIF


!-----------------------------------------------------------------------
! DEFINE GRID NAMESPACE

  CALL GRID_DEFINE_FROM_LIST(gd, microphys,ipconc, ichaff, inucopt, ienkf=ienkf,isfcphys=isfcphys)

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  NAMESPACE DEFINED'

! READ IN THE GRID ATTRIBUTES

  DO n = 1,size(gd%attr)

    CALL STRING_LIMITS(gd%attr(n)%name, ibeg, iend)

    IF( gd%attr(n)%type .eq. 'int' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%int)
    IF( gd%attr(n)%type .eq. 'flt' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%flt)
    IF( gd%attr(n)%type .eq. 'str' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%str)

    IF (status /= NF90_NOERR) THEN
      write(lundbg,*) 'GRID_OPEN_NETCDF READ_ATTRIBUTE:  Problem reading attribute',gd%attr(n)%type
      write(lundbg,*) 'GRID_OPEN_NETCDF READ_ATTRIBUTE:  Filename = ', filename
      write(lundbg,*) 'GRID_OPEN_NETCDF READ_ATTRIBUTE:  NAME = ', gd%attr(n)%name
    ENDIF

    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  ATTRIBUTES= ', gd%attr(n)%name(ibeg:iend), &
                              gd%attr(n)%int, gd%attr(n)%flt, gd%attr(n)%str

! GET GRID ATTRIBUTES SPECIFYING DIMENSIONS

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NXEND' ) nxend = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NYEND' ) nyend = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NZEND' ) nzend = gd%attr(n)%int

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NX' ) nx = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NY' ) ny = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NZ' ) nz = gd%attr(n)%int

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'SIZE_OF_ENSEMBLE' ) ensize = gd%attr(n)%int
    
  ENDDO
  
  

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  ATTRIBUTES READ IN'
 
!-----------------------------------------------------------------------------------------------
! ALLOCATE SPACE FOR VARIABLES

!
!  Get values of XE,YE,ZE dimensions.  This, along with NX,NY,NZ can tell us if the file has historyoutput 0 or 1.
!

  status = nf90_inq_dimid(ncid, "XE", XEDimID)
  status = nf90_inq_dimid(ncid, "YE", YEDimID)
  status = nf90_inq_dimid(ncid, "ZE", ZEDimID)

  status = nf90_inquire_dimension(ncid, XEDimID, len = ncxe)
  status = nf90_inquire_dimension(ncid, YEDimID, len = ncye)
  status = nf90_inquire_dimension(ncid, ZEDimID, len = ncze)
  
  
  IF ( present( nx_or_nxend ) ) THEN
   IF ( nx_or_nxend ) THEN
    CALL GRID_ALLOCATE(gd, ncxe, ncye, ncze) ! allocate for WHOLE available domain, not just a tile
! set these for when code is compiled with MPI
    ixend = ncxe
    jyend = ncye
    kzend = ncze
   ELSE
    CALL GRID_ALLOCATE(gd, nx, ny, nz)
   ENDIF
  
  ELSE
  
    CALL GRID_ALLOCATE(gd, nx, ny, nz)

  ENDIF
  
!#else 
!  CALL GRID_ALLOCATE(gd, nx, ny, nz)
!#endif




END SUBROUTINE GRID_OPEN_NETCDF

SUBROUTINE GRID_CLOSE_NETCDF(file,ncid)

#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  character(LEN=*)  :: file

  integer ncid, record_id, variable_id
  integer cmode
  integer n, status
  integer ibeg,iend
  character(LEN=80) filename


!    CALL STRING_LIMITS(file, ibeg, iend)
!    filename(ibeg:iend) = file(ibeg:iend)
  
!  status = NF90_OPEN(filename(ibeg:iend),NF90_NOWRITE,ncid)
!mpi! #if defined (NC4) && defined (MPI)
!mpi!   cmode  = ior(NF90_NOWRITE,nf90_iotype)
!mpi!   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
!mpi! #else
!  cmode  = NF90_NOWRITE
!  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
!  if (status /= nf90_noerr) call handle_err(status)
!mpi! #endif

  
  status=NF90_CLOSE(ncid)
  if (status /= nf90_noerr) call handle_err(status)


END SUBROUTINE GRID_CLOSE_NETCDF

!===============================================================================
!     
!     
!
!    
!
!   
!
!===============================================================================

SUBROUTINE GRID_WRITE_NETCDF(gd, start, file, forcewrite, ncidopen)

#ifdef MPI
  USE COMMASMPI_MODULE
#endif
  
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
  integer ncid, nt, record_id, records
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

  
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  START = ', start

  

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  checkhistory = ', checkhistory
  
#if defined (NC4)
 IF ( parallelio ) THEN
#if defined (MPI)
   IF ( parallelio_type /= 2 ) THEN ! HDF5 parallel
   cmode  = ior( ior(NF90_WRITE,nf90_iotype), NF90_MPIIO )
   ELSE ! pnetcdf
   cmode  = ior(NF90_WRITE, ior(pnetcdf_flag,NF90_MPIIO))
!   status = NF90_OPEN(filename(ibeg:iend),cmode,ncid,comm=my_comm,info=my_info)
   ENDIF

   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
#endif
 ELSE
!mpi! #else
  cmode  = NF90_WRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
 ENDIF
#else
  cmode  = NF90_WRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
#endif
!mpi! #endif
  
  IF (status /= NF90_NOERR) THEN 
#ifdef MPI
    write(0,*) 'WRITE_NETCDF:  Problem opening file'
    write(0,*) 'WRITE_NETCDF:  Filename = ', filename(ibeg:iend)
    write(0,*) 'WRITE_NETCDF:  MY_COMM = ', my_comm
    write(0,*) 'WRITE_NETCDF:  MY_INFO = ', my_info
    write(0,*) 'WRITE_NETCDF:  MY_RANK = ', my_rank
    write(0,*) 'WRITE_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
#else
    write(lundbg,*) 'WRITE_NETCDF:  Problem opening file'
    write(lundbg,*) 'WRITE_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'WRITE_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
#endif
  ENDIF


  ELSE
    ncid = ncidopen
  ENDIF ! present(ncidopen)

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF: NCID = ', ncid, NF90_STRERROR(status)

!-----------------------------------------------------------------------------------------------
! Inquire as to how many time steps have been written to determine where to put the data%%%

  status = NF90_INQUIRE(ncid, unlimitedDimID = record_id)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'GRID_WRITE_NETCDF: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'GRID_WRITE_NETCDF: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
#endif

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  record_id = ', record_id

  status = NF90_INQUIRE_DIMENSION(ncid, record_id, len = records)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'GRID_WRITE_NETCDF: ERROR READING TIME DIM ',NF90_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'GRID_WRITE_NETCDF: ERROR READING TIME DIM ',NF90_STRERROR(status)
#endif  

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  RECORDS = ', records

!-----------------------------------------------------------------------------------------------
! Get model time array
  
    CALL GET_VARIABLE(gd, 'TIME', time) ! note that if the "+1" time was read, then TIME got set to that. So will write back to the same time.
    CALL GET_VARIABLE(gd, 'THISTORY', thistory)
    CALL GET_VARIABLE(gd, 'UGRID', ugrid)
    CALL GET_VARIABLE(gd, 'VGRID', vgrid)
    CALL GET_VARIABLE(gd, 'DT', dt)

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  TIME = ', time, thistory, ugrid, vgrid, dt

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

    status = NF90_INQ_VARID(ncid,'TIME',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE TIME ',NF90_STRERROR(status)
    status = NF90_GET_VAR(ncid, variable_id, tarray)
    IF(status /= NF90_NOERR) THEN
    write(0,*) 'WRITE_NETCDF:  ERROR READING COORDINATE TIME ',NF90_STRERROR(status)
    write(0,*) 'WRITE_NETCDF: ncid, variable_id, tarray',ncid, variable_id, records, tarray(1:records)
    ENDIF

    status = NF90_INQ_VARID(ncid,'THISTORY',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE THISTORY ',NF90_STRERROR(status)
    status = NF90_GET_VAR(ncid, variable_id, tharray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR READING COORDINATE THISTORY ',NF90_STRERROR(status)

    status = NF90_INQ_VARID(ncid,'UGRID',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE UGRID ',NF90_STRERROR(status)
    status = NF90_GET_VAR(ncid, variable_id, ugarray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR READING UGRID ',NF90_STRERROR(status)

    status = NF90_INQ_VARID(ncid,'VGRID',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE VGRID ',NF90_STRERROR(status)
    status = NF90_GET_VAR(ncid, variable_id, vgarray)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR READING VGRID ',NF90_STRERROR(status)

    status = NF90_INQ_VARID(ncid,'DT',variable_id)
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE DT ',NF90_STRERROR(status)
    status = NF90_GET_VAR(ncid, variable_id, dtarray)
   IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR READING DT ',NF90_STRERROR(status)
  
! Search to find the correct time in the netcdf file
  
    index(:) = any(tarray(:) .eq. time) 
  
    IF( count(index) .EQ. 0 .AND. time .LE. tarray(records) ) THEN

#ifdef MPI
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*)
      write(0,*) 'GRID_WRITE_NETCDF: ERROR, CANNOT FIND A MATCHING TIME '
      write(0,*)
  
      DO n = 1,records
        write(0,*) 'GRID_WRITE_NETCDF: INPUT TIME ', time, ' DATA TIMES: ', tarray(n), index(n)
      ENDDO
  
      write(0,*)
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*) '==================== **** ERROR *** ========================'
      write(0,*)
#else
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*)
      write(6,*) 'GRID_WRITE_NETCDF: ERROR, CANNOT FIND A MATCHING TIME '
      write(6,*)
  
      DO n = 1,records
        WRITE(6,*) 'GRID_WRITE_NETCDF: INPUT TIME ', time, ' DATA TIMES: ', tarray(n), index(n)
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
      write(0,*) 'GRID_WRITE_NETCDF: Overwriting time ',tarray(nt),index(nt)
#else
      write(6,*) 'GRID_WRITE_NETCDF: Overwriting time ',tarray(nt),index(nt)
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
         print *,'GRID_WRITE_NETCDF: ugrid,vgrid,dt = ',ugrid,ugarray(records),vgrid,vgarray(records),dt,dtarray(records)
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
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  INPUT TIME ', time
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  TIME ARRAY(NT) ', tarray(Min(nt,records))
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  INDEX ARRAY(:) ', index
#else
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  INPUT TIME ', time
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  TIME ARRAY(NT) ', tarray(Min(nt,records))
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  INDEX ARRAY(:) ', index
#endif

    deallocate(tarray)
    deallocate(index)
    deallocate(tharray)
    deallocate(ugarray)
    deallocate(vgarray)
    deallocate(dtarray)
  
  ENDIF

#ifdef MPI
  IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  NT = ', nt
#else
  IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  NT = ', nt
#endif

!-----------------------------------------------------------------------------------------------
! Add to time array

  IF( DEBUG_IO ) write(6,*) 'WRITE_NETCDF: writing time information, NT = ', nt
  
  status = NF90_INQ_VARID(ncid,'TIME',variable_id)

#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE TIME ',NF90_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE TIME ',NF90_STRERROR(status)
#endif
  
  status = NF90_NOERR
  IF ( my_rank == 0 .or. historyoutput == 0 ) THEN
!  status = NF90_PUT_VAR(ncid, variable_id, time, (/ nt /) )
  ENDIF

#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE TIME',NF90_STRERROR(status), &
      'my_rank =',my_rank,nt,time,status
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE TIME',NF90_STRERROR(status)
#endif

  IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  WROTE TIME INFORMATION'

!-----------------------------------------------------------------------
!  NetCDF:  Write out grid dimensions

  IF( nt == 1 .and. .false.) THEN ! {

    CALL GET_VARIABLE(gd, 'XC', var)
    nx = size(var%flt1d) - 1 + var%istag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE XC ',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE XC ',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
   if (ixend /= nxend) nx = nx + 1
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*)  &
      'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS',NF90_STRERROR(status),'XC ON MY_RANK',my_rank
   ENDIF
#endif
      status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx), (/ ixbeg /))
   ELSE
      status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE XC ',var%name, NF90_STRERROR(status)
    IF( DEBUG_IO ) write(0,*) 'ixbeg, nx = ',ixbeg, nx
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  XC DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE XC ',var%name, NF90_STRERROR(status)
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  XC DIMENSION SCALE WRITTEN'
#endif
    

    CALL GET_VARIABLE(gd, 'XE', var)
    nx = size(var%flt1d) - 1 + var%istag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE XE ',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE XE ',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS', &
      NF90_STRERROR(status),'XE ON MY_RANK',my_rank
   ENDIF
#endif
      status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx), (/ ixbeg /))
   ELSE
      status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nx))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE XE ',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  XE DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE XE ',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  XE DIMENSION SCALE WRITTEN'
#endif    


    CALL GET_VARIABLE(gd, 'YC', var)
    ny = size(var%flt1d) - 1 + var%jstag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE YC ',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE YC ',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
     if (jyend /= nyend) ny = ny + 1
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS', &
        NF90_STRERROR(status),'YC ON MY_RANK',my_rank
   ENDIF
#endif
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny), (/ jybeg /))
   ELSE
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE YC ',var%name, NF90_STRERROR(status)   
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  YC DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE YC ',var%name, NF90_STRERROR(status)   
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  YC DIMENSION SCALE WRITTEN'
#endif    


    CALL GET_VARIABLE(gd, 'YE', var)
    ny = size(var%flt1d) - 1 + var%jstag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE YE',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE YE',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS', &
       NF90_STRERROR(status),'YE ON MY_RANK',my_rank
   ENDIF
#endif
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny), (/ jybeg /))
   ELSE
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:ny))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE YE',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  YE DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE YE',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  YE DIMENSION SCALE WRITTEN'
#endif
    

    CALL GET_VARIABLE(gd, 'ZC', var)
    nz = size(var%flt1d) - 1 + var%kstag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE ZC',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE ZC',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
   if (kzend /= nzend) nz = nz + 1
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS', &
       NF90_STRERROR(status),'ZC ON MY_RANK',my_rank
   ENDIF
#endif
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz), (/ kzbeg /))
   ELSE
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE ZC',var%name, NF90_STRERROR(status)
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  ZC DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE ZC',var%name, NF90_STRERROR(status)
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  ZC DIMENSION SCALE WRITTEN'
#endif
    

    CALL GET_VARIABLE(gd, 'ZE', var)
    nz = size(var%flt1d) - 1 + var%kstag - 2*var%ng !!mpidebug note the removal of ghost zones
    status = NF90_INQ_VARID(ncid,var%name,variable_id)
#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE ZE',var%name, NF90_STRERROR(status)
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR INQUIRE COORINDATE ZE',var%name, NF90_STRERROR(status)
#endif

#ifdef MPI
   IF ( Abs(historyoutput) .eq. 1 .or. historyoutput == 3 ) THEN ! single file
#ifdef NC4
   IF (  parallelio ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR ALLOWING PARALLEL ACCESS', &
       NF90_STRERROR(status),'ZE ON MY_RANK',my_rank
   ENDIF
#endif
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz), (/ kzbeg /))
   ELSE
     status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz), (/ 1 /))
   ENDIF
#else
    status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:nz))
#endif

#ifdef MPI
    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE ZE',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF:  ZE DIMENSION SCALE WRITTEN'
#else
    IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF:  ERROR WRITING COORDINATE ZE',NF90_STRERROR(status)
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF:  ZE DIMENSION SCALE WRITTEN'
#endif



#ifdef MPI
  IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF: GRID DIMENSIONS ARE WRITTEN OUT FOR NT=1'
#else
  IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF: GRID DIMENSIONS ARE WRITTEN OUT FOR NT=1'
#endif
    
  ENDIF ! } nt == 1
 
!-----------------------------------------------------------------------
!  NetCDF:  Define variables - scalar arrays are written out as single xyz3d arrays

  DO n = 1,size(gd%var)
  
#ifdef MPI
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF: WRITING OUT VARIABLE:  ', gd%var(n)%name
#else
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF: WRITING OUT VARIABLE:  ', gd%var(n)%name
#endif

    CALL GET_VARIABLE(gd, gd%var(n)%name, var)
    CALL WRITE_NCDF_VAR(ncid,var,nt)
#ifdef MPI
    IF ( parallelio ) CALL MPI_BARRIER(my_comm, mpi_error_code) ! nf90_PNETCDF
    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_NETCDF: VARIABLE:  ', gd%var(n)%name,  &
       ' WAS WRITTEN TO THE FILE. my_rank,time=',my_rank,time
#else  
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_NETCDF: VARIABLE:  ', gd%var(n)%name, ' WAS WRITTEN TO THE FILE'
#endif
  
  ENDDO

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)
  IF ( .not. present( ncidopen ) ) THEN
   status=NF90_CLOSE(ncid)
#ifdef MPI
  IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF: Error Closing File: ', NF90_STRERROR(status)
#else
  IF(status /= NF90_NOERR) write(6,*) 'WRITE_NETCDF: Error Closing File: ', NF90_STRERROR(status)
#endif
  ENDIF

  IF ( DEBUG_IO ) THEN
#ifdef MPI
  write(0,*) ''
  write(0,*) 'GRID_WRITE_NETCDF:  FILE: ',filename(ibeg:iend),'  SUCCESSFULLY WRITTEN'
  write(0,*) ''
#else
  write(lundbg,*) ''
  write(lundbg,*) 'GRID_WRITE_NETCDF:  FILE: ',filename(ibeg:iend),'  SUCCESSFULLY WRITTEN'
  write(lundbg,*) ''
#endif
  ENDIF

END SUBROUTINE GRID_WRITE_NETCDF






!===============================================================================
! 
!
!  
!  
!
!  
!  
!===============================================================================

SUBROUTINE DEFINE_NCDF_VAR(ncid,dim_id,var,dimvals)

#ifdef MPI
  USE mpi
#endif
  USE RUN_ATT_NML, only: microphys
  USE NETCDF 
  
  implicit none

  integer ncid, dim_id(0:8), dimvals(0:8), sizes(1:4)
  type(variable) :: var

  integer, allocatable :: dimids(:)
  integer dim, nid, status, i
  integer, allocatable :: chunksizes(:) , extend_increments(:)
  logical isinteger
  character(20) :: tmpstring

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, write out the scalar 

  dim = var%dim + var%tdepend

  IF( dim == 0 ) THEN

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR: var%name ',var%name
    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_DEF_VAR(ncid, var%name, NF90_INT, nid)
      isinteger = .true.
    ELSE
      status = NF90_DEF_VAR(ncid, var%name, NF90_FLOAT, nid)
      isinteger = .false.
    ENDIF
    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR: var%name,isinteger ', var%name,isinteger

    IF(status /= NF90_NOERR) write(6,*) 'DEFINE_NCDF_VAR:  Error creating variable: ', var%name
    if (status /= nf90_noerr) call handle_err(status)

!-------------------------------------------------------------------------------
! Else we are writing out a 1, 2, 3D array - Use a trick about the type labeling to cut down
!      on the number of conditionals

  ELSE

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR 2: var%name, dim ',var%name, dim
    allocate(dimids(dim))
    
    dimids(:) = 0

    IF( var%tdepend /= 0 ) dimids(dim) = dim_id(0)

    IF( var%type(1:1) == 'x' ) THEN               ! x1d set
      IF( var%istag == 1) THEN 
        dimids(1) = dim_id(2)
        sizes(1) = dimvals(2)
      ELSE
        dimids(1) = dim_id(1)
        sizes(1) = dimvals(1)
      ENDIF
    ENDIF

    IF( var%type(1:1) == 'y' ) THEN               ! y1d set
      IF( var%jstag == 1) THEN 
        dimids(1) = dim_id(4)
      ELSE
        dimids(1) = dim_id(3)
      ENDIF
    ENDIF

    IF( var%type(1:1) == 'z' ) THEN               ! z1d set
      IF( var%kstag == 1) THEN 
        dimids(1) = dim_id(6)
      ELSE
        dimids(1) = dim_id(5)
      ENDIF
    ENDIF

    IF( var%type(2:2) == 'y' ) THEN               ! xy2d set
      IF( var%jstag == 1) THEN 
        dimids(2) = dim_id(4)
        sizes(2) = dimvals(4)
      ELSE
        dimids(2) = dim_id(3)
        sizes(2) = dimvals(3)
      ENDIF
    ENDIF

    IF( var%type(2:2) == 'z' ) THEN               ! xz2d & yz2d set
      IF( var%kstag == 1) THEN 
        dimids(2) = dim_id(6)
      ELSE
        dimids(2) = dim_id(5)
      ENDIF
    ENDIF

    IF( var%type(2:2) == '2' ) THEN               ! x2d & y2d & z2d set
        dimids(2) = dimids(1)
        dimids(1) = dim_id(8)
    ENDIF

    IF( var%type(3:3) == 'z' ) THEN               ! xyz3d set
      IF( var%kstag == 1) THEN 
        dimids(3) = dim_id(6)
        sizes(3) = dimvals(6)
      ELSE
        dimids(3) = dim_id(5)
        sizes(3) = dimvals(5)
      ENDIF
    ENDIF

! Write out information

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR: var%name, dimids ',var%name, dimids

    IF( var%type(1:5) == 'icnst' ) THEN
      IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR INT 1: var%name ',var%name, ncid
      status = NF90_DEF_VAR(ncid, var%name, NF90_INT, dimids, nid)
      isinteger = .true.
      IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR INT: var%name, dimids, nid ',var%name, dimids, nid
#ifdef NC4
     IF ( netcdfversion .eq. 4 .and. parallelio_type /= 2) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, (/ 300 /) )
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF
#endif
    ELSE
      status = NF90_DEF_VAR(ncid, var%name, NF90_real, dimids, nid)
      isinteger = .false.

#if defined (NC4)
      IF( netcdfversion .eq. 4 .and. var%type(1:5) == 'rcnst' ) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, (/ 300 /) )
        if (status /= nf90_noerr) call handle_err(status)
      ENDIF

     IF ( netcdfversion .eq. 4 .and. parallelio_type /= 2) THEN !{
        IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR NC4 setup: var%name ',var%name

      IF( var%type(3:3) == 'z' .or. var%type(1:2) == 'xy' ) THEN               ! xyz3d set
!      IF( var%type(3:3) == 'z' ) THEN               !{ xyz3d set
        allocate( chunksizes(4) )
!        allocate( extend_increments(4) ) 
        chunksizes(1:3) = sizes(1:3)
        
        IF ( var%type(3:3) == 'z' ) THEN
!        IF ( sizes(2) > 2 .and. nprock == 1 ) THEN ! not 2D and not tiled in vertical
        IF ( sizes(2) > 2 ) THEN ! not 2D and not tiled in vertical
          IF ( 4*(sizes(3)/4) == sizes(3) ) THEN
          chunksizes(3)   = sizes(3)/4
          ELSE
          chunksizes(3)   = sizes(3)/4 + 1
          ENDIF
        ENDIF
        ELSE
          chunksizes(3)   = 1
        ENDIF 
        chunksizes(4)   = 1
!        extend_increments(1:4) = chunksizes(1:4)
#ifdef MPI
        ! send chunk sizes from rank 0 so that all have the same value
         IF ( parallelio ) THEN
         CALL MPI_Bcast(chunksizes, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_error_code)
         ENDIF
         
#endif
!         IF ( var%name == 'TH' ) THEN
!           write(0,*) 'rank, TH chunksizes: ',my_rank,chunksizes
!         ENDIF
         IF ( (.not. parallelio) .or. ( parallel_compress .and. var%type(3:3) == 'z'  ) ) THEN
!         IF ( (.not. parallelio) .or. ( parallel_compress .and. .false. ) ) THEN
!        write(0,*) 'DEFINE_NCDF_VAR rank, chunksizes = ',my_rank,chunksizes(:)
!        IF ( var%type(3:3) == 'z' ) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error setting chunking for: ', var%name
!        ENDIF
         ENDIF
!        deallocate( extend_increments ) 

!        IF ( (.not. parallelio) .or. ( var%name(1:1) == 'Q' .and. parallel_compress .and. var%type(3:3) == 'z') ) THEN   
        IF ( (.not. parallelio) .or. ( parallel_compress .and. var%type(3:3) == 'z') ) THEN   
!        IF ( (.not. parallelio) .or. ( parallel_compress .and. .false.) ) THEN   
         IF ( DEBUG_IO ) write(0,*) 'set deflate for ',var%name,var%type,shuffle, deflate, deflate_level
         status = NF90_DEF_VAR_DEFLATE(ncid, nid, shuffle, deflate, deflate_level)
!         write(0,*) 'deflate, deflate_level, shuffle = ',deflate, deflate_level, shuffle
         IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error setting deflate for: ', var%name
        ENDIF

        deallocate( chunksizes )
        ENDIF!}
        
     ELSEIF ( netcdfversion .eq. 4 .and. parallelio .and. parallelio_type /= 2 ) THEN ! } {
        IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VAR NC4 parallel setup: var%name ',var%name
      IF( var%type(3:3) == 'z' ) THEN               ! xyz3d set
        allocate( chunksizes(4) )
!        allocate( extend_increments(4) ) 
       chunksizes(1:3) = sizes(1:3)
        IF ( sizes(2) > 2 ) THEN ! not 2D
          IF ( 4*(sizes(3)/4) == sizes(3) ) THEN
          chunksizes(3)   = sizes(3)/4
          ELSE
          chunksizes(3)   = sizes(3)/4 + 1
          ENDIF
        ENDIF
!        extend_increments(1:4) = chunksizes(1:4)

#ifdef MPI
        ! send chunk sizes from rank 0 so that all have the same value
         IF ( parallelio ) THEN
         CALL MPI_Bcast(chunksizes, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_error_code)
         ENDIF
         
#endif
        
!        write(0,*) 'DEFINE_NCDF_VAR: sizes parallel = ',var%name,chunksizes(:)
!        write(0,*) 'DEFINE_NCDF_VAR: var%kstag = ',var%kstag,dimvals(5),dimvals(6)
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error setting chunking for: ', var%name
        deallocate( chunksizes )
!        deallocate( extend_increments ) 
      ENDIF

     ENDIF !}
#endif
    ENDIF

    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable: ', var%name
    if (status /= nf90_noerr) call handle_err(status)

    deallocate(dimids)

  ENDIF

! DEFINE ATTRIBUTES

  status = NF90_PUT_ATT(ncid, nid, "type", var%type)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute long_name: ', var%name

  status = NF90_PUT_ATT(ncid, nid, "long_name", var%description)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute long_name: ', var%name

  status = NF90_PUT_ATT(ncid, nid, "units",     var%unit)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute units:     ', var%name

  IF ( .not. ( parallelio .and. parallelio_type == 1 .and. microphys(1:3) == 'TAK' ) ) THEN
!   IF ( .true. ) THEN
  status = NF90_PUT_ATT(ncid, nid, "index",     var%index)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute index:     ', var%name

!  ENDIF
!  IF ( .false. ) THEN

!  ELSE
!   i = var%index
!   write(tmpstring,'(i5)' ) var%index
!  status = NF90_PUT_ATT(ncid, nid, "index",     tmpstring)
!  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute index:     ', var%name

!  ENDIF


  status = NF90_PUT_ATT(ncid, nid, "posdef",    var%pdef)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute posdef:    ', var%name

  status = NF90_PUT_ATT(ncid, nid, "istag",      var%istag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute istag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "jstag",      var%jstag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute jstag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "kstag",      var%kstag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute kstag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "dyntype",   var%dyntype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute dyntype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "phytype",   var%phytype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute phytype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "buotype",   var%buotype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute buotype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "basedim",   var%basedim)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute basedim:   ', var%name

  ENDIF ! true/false


  IF ( var%field .ne. CHAR(0) ) THEN
    status = NF90_PUT_ATT(ncid, nid, "field",   var%field)
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute field:   ', var%name
  ENDIF
  
  IF ( var%positions .ne. CHAR(0) ) THEN
    status = NF90_PUT_ATT(ncid, nid, "positions",   var%positions)
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute positions:   ', var%name
  ENDIF

! for CF conventions and VAPOR import
  IF ( var%name == 'XC' .or. var%name == 'XE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'X')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute axis:   ', var%name
  ENDIF

  IF ( var%name == 'YC' .or. var%name == 'YE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Y')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute axis:   ', var%name
  ENDIF

  IF ( var%name == 'ZC' .or. var%name == 'ZE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Z')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute axis:   ', var%name
  ENDIF

  IF ( historyoutput == -1 .and. var%dim .gt. 0 ) THEN
  IF ( isinteger ) THEN
    status = NF90_PUT_ATT(ncid, nid, "_FillValue",   0)
  ELSE
    status = NF90_PUT_ATT(ncid, nid, "_FillValue",   0.0)
  ENDIF
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error creating variable attribute _FillValue:   ', var%name
  ENDIF
  

END SUBROUTINE DEFINE_NCDF_VAR




!===============================================================================
!
!
!
!
!
!
!
!===============================================================================

SUBROUTINE WRITE_NCDF_VAR(ncid,var,nt)

  USE NETCDF 
#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  integer ncid, nt
  TYPE(variable) :: var

  integer variable_id
  integer, allocatable :: dims(:),count(:)
  integer dim, status, n1, n2, n3, stag, istag, jstag, kstag
  logical DEBUG_IO2
  real, allocatable :: buf3d(:,:,:)
  
  status = NF90_INQ_VARID(ncid,var%name,variable_id)
  IF(status /= NF90_NOERR) THEN
    print *, 'WRITE_NETCDF_VAR:  ERROR INQUIRE VARIABLE ',var%name, NF90_STRERROR(status)
    RETURN
  ENDIF

#if defined (NC4) && defined (MPI)
  IF ( parallelio ) THEN
    status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)

    IF(status /= NF90_NOERR) write(0,*) 'WRITE_NETCDF_VAR:  ERROR ALLOWING PARALLEL ACCESS', &
  NF90_STRERROR(status),'ON MY_RANK',my_rank, 'variable_id = ',variable_id
    IF( DEBUG_IO ) write(0,*) 'WRITE_NETCDF_VAR:  PARALLEL VARIABLE ACCESS ALLOWED'
  ENDIF
#endif


!-------------------------------------------------------------------------------
! DEBUG STATEMENTS

  IF( DEBUG_IO ) THEN
    print *, var%name, var%type, var%istag, var%jstag, var%kstag 
  ENDIF

  DEBUG_IO2 = DEBUG_IO

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, write out the scalar

  dim = var%dim + var%tdepend

  IF( dim == 0 ) THEN

    IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR: WRITING OUT VARIABLE:  ',var%name,var%type

    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_PUT_VAR(ncid, variable_id, var%int)
      IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR: WRITING OUT VARIABLE:  ',var%name,' = ',var%int
    ELSE
      status = NF90_PUT_VAR(ncid, variable_id, var%flt)
      IF( DEBUG_IO2 ) write(6,*) 'WRITE_NCDF_VAR: WRITING OUT VARIABLE:  ',var%name,' = ',var%flt
    ENDIF

    IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)

!-------------------------------------------------------------------------------
! Else we are writing out a 1, 2, 3D array 

  ELSE

    allocate(dims(dim))
    allocate(count(dim))

    dims(:) = 1

    IF( var%tdepend /= 0 ) dims(dim) = nt

    IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR:  nt,dim = ',nt,dim


    IF( var%type(1:5) == 'rcnst' .or. var%type(1:5) == 'icnst'  ) THEN !{
      
    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_PUT_VAR(ncid, variable_id, var%int, dims)
      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR: WRITING OUT VARIABLE:  ',var%name,' = ',var%int
    ELSE
      status = NF90_PUT_VAR(ncid, variable_id, var%flt, dims)
      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR: WRITING OUT VARIABLE:  ',var%name,' = ',var%flt
    ENDIF


      IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)

      IF( DEBUG_IO2 ) print *, 'WRITE_NCDF_VAR:  ',var%name
    
    ENDIF !}


    IF( var%type(1:3) == 'x1d' .or. var%type(1:3) == 'y1d' .or. var%type(1:3) == 'z1d' ) THEN !{

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
      status = NF90_PUT_VAR(ncid, variable_id, var%flt1d(1:n1),  dims)
 
!      IF(status /= NF90_NOERR) print *,'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)
      IF(status /= NF90_NOERR) THEN
        print *,'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)
#ifdef MPI
        print *,'ixbeg,ixend,nxbeg,nxend = ',ixbeg,ixend,nxbeg,nxend
#endif
        print *, 'dims = ',dims, ' stag = ',stag
      ENDIF

      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR:  ',var%name, n1

    ENDIF !}

    IF( var%type(1:3) == 'x2d' .or. var%type(1:3) == 'y2d' .or. var%type(1:3) == 'z2d' ) THEN !{

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

      status = NF90_PUT_VAR(ncid, variable_id, var%flt2d(1:n1,1:n2), dims)

      IF(status /= NF90_NOERR) THEN
        write(0,*) 'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)
        RETURN
      ENDIF
      IF( DEBUG_IO2 ) write(0,*) 'WRITE_NCDF_VAR:  ',var%name, n1, n2

    ENDIF !}

    IF( var%type(1:4) == 'xy2d' .or. var%type(1:4) == 'xz2d' .or. var%type(1:4) == 'yz2d' ) THEN !{

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

      IF( DEBUG_IO2 ) then 
        write(0,*)  'WRITE_NCDF_VAR:  ',var%name, n1, n2
        write(0,*) ' dims = ',dims(1),dims(2),dims(3)
      ENDIF

      IF ( .not. ( myprock > 1 .and. var%type(1:4) == 'xy2d' )) THEN
      status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
      status = NF90_PUT_VAR(ncid, variable_id, var%flt2d(1:n1,1:n2), start=dims, count=(/ n1,n2,1 /) ) !,stride=(/1,1,1/) )

      IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)

      IF( DEBUG_IO2 ) write(0,*)  'WRITE_NCDF_VAR:  ',var%name, n1, n2
      
      ENDIF

    ENDIF !}

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
      
!      write(0,*) 'WRITE_NCDF_VAR: n1,n2,n3 = ',n1,n2,n3, var%name

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
      count = (/ n1,n2,n3,1 /) 

      IF( DEBUG_IO2 ) THEN
        write(0,*) 'WRITE_NCDF_VAR:  ',var%name, n1, n2, n3
        write(0,*) 'dims = ',dims
      ENDIF

!      allocate(buf3d(n1,n2,n3))
      
!      buf3d(1:n1,1:n2,1:n3) = var%flt3d(1:n1,1:n2,1:n3)
!      IF ( parallelio ) TH
      status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
!      status = NF90_PUT_VAR(ncid, variable_id, var%flt3d(1:n1,1:n2,1:n3), dims)
      status = NF90_PUT_VAR(ncid, variable_id, var%flt3d(1:n1,1:n2,1:n3), start=dims, count=count)
!      status = NF90_PUT_VAR(ncid, variable_id, buf3d, start=dims, count=(/ n1,n2,n3,1 /),stride=(/1,1,1,1/) )
!      deallocate( buf3d )
      IF(status /= NF90_NOERR) write(0,*)'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)

    ENDIF

    deallocate(dims)
    deallocate(count)

  ENDIF

  IF(status /= NF90_NOERR) write(0,*) 'WRITE_NCDF_VAR:  Error writing variable: ', var%name, NF90_STRERROR(status)

END SUBROUTINE WRITE_NCDF_VAR



!===============================================================================
!
!
!
!   /////////////////////                BEGIN          \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\   SUBROUTINE READ_ATTRIBUTE   ////////////////////
!
!
!
!
! Created by Louis Wicker, December 10, 1990
!
! Latest update: 03-03-06 --> adapted for PYENCOMMAS
!   
!===============================================================================
SUBROUTINE READ_ATTRIBUTE(filename, name, long, float, string)

  USE NETCDF
#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none
   
! Passed variables
   
  character(LEN = *), intent(in)            :: filename
  character(LEN = *), intent(in)            :: name
  integer,            intent(out), optional :: long
  real,               intent(out), optional :: float
  character(LEN = *), intent(out), optional :: string

! Local variables

  integer ncid
  integer cmode
  integer status
  integer ibeg
  integer iend
    
!-----------------------------------------------------------------------
! Open netcdf file...

  CALL STRING_LIMITS(filename, ibeg, iend)

#if defined (NC4) && defined (MPI)
  IF ( parallelio ) THEN

   IF ( parallelio_type /= 2 ) THEN ! HDF5 parallel
   cmode  = ior( ior(NF90_NOWRITE,nf90_iotype),NF90_MPIIO )
   ELSE ! pnetcdf
   cmode  = ior(NF90_NOWRITE, ior(pnetcdf_flag,NF90_MPIIO))
!   status = NF90_OPEN(filename(ibeg:iend),cmode,ncid,comm=my_comm,info=my_info)
   ENDIF


   status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
 
   write(0,*) 'READ_ATTRIBUTE:  NETCDF FILE OPENED SUCCESSFULLY FOR PARALLEL IO', NF90_STRERROR(status)
   write(0,*) 'READ_ATTRIBUTE:  MY_RANK=',my_rank,'MY_INFO=',my_info,'MY_COMM=',my_comm
  ELSE
  cmode  = NF90_NOWRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
  ENDIF
#else
  cmode  = NF90_NOWRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
#endif

#ifdef MPI
  IF (status /= NF90_NOERR) THEN
    write(0,*) 'READ_ATTRIBUTE:  Problem opening file'
    write(0,*) 'READ_ATTRIBUTE:  Filename = ', filename(ibeg:iend)
    write(0,*) 'READ_ATTRIBUTE:  MY_COMM = ', my_comm
    write(0,*) 'READ_ATTRIBUTE:  MY_INFO = ', my_info
    write(0,*) 'READ_ATTRIBUTE:  MY_RANK = ', my_rank
    write(0,*) 'READ_ATTRIBUTE:  NCID = ', ncid, NF90_STRERROR(status)
    CALL commasmpi_abort()
    stop
  ENDIF
#else
  IF (status /= NF90_NOERR) THEN
    write(lundbg,*) 'READ_ATTRIBUTE:  Problem opening file'
    write(lundbg,*) 'READ_ATTRIBUTE:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'READ_ATTRIBUTE:  NCID = ', ncid, NF90_STRERROR(status)
    stop
  ENDIF
#endif

  IF( PRESENT(long) ) THEN

    status = NF90_GET_ATT(ncid, NF90_GLOBAL, name, long)

    IF (status /= NF90_NOERR) THEN 
      write(lundbg,*) 'READ_ATTRIBUTE:  Problem reading long attribute'
      write(lundbg,*) 'READ_ATTRIBUTE:  Filename = ', filename
      write(lundbg,*) 'READ_ATTRIBUTE:  NAME = ', name
    ENDIF

  ENDIF

  IF( PRESENT(float) ) THEN

    status = NF90_GET_ATT(ncid, NF90_GLOBAL, name, float)

    IF (status /= NF90_NOERR) THEN 
      write(lundbg,*) 'READ_ATTRIBUTE:  Problem reading float attribute'
      write(lundbg,*) 'READ_ATTRIBUTE:  Filename = ', filename
      write(lundbg,*) 'READ_ATTRIBUTE:  NAME = ', name
    ENDIF

  ENDIF

  IF( PRESENT(string) ) THEN

    status = NF90_GET_ATT(ncid, NF90_GLOBAL, name, string)

    IF (status /= NF90_NOERR) THEN 
      write(lundbg,*) 'READ_ATTRIBUTE:  Problem reading string attribute'
      write(lundbg,*) 'READ_ATTRIBUTE:  Filename = ', filename
      write(lundbg,*) 'READ_ATTRIBUTE:  NAME = ', name
    ENDIF

  ENDIF

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)

  status=NF90_CLOSE(ncid)
  IF(status /= NF90_NOERR) print *, 'DEFINE_NETCDF: Error Closing File: ', NF90_STRERROR(status)

RETURN
END SUBROUTINE READ_ATTRIBUTE





!===============================================================================
!     
!     
!
!    
!
!   
!
!===============================================================================

SUBROUTINE GRID_READ_NETCDF(gd, filename, start0, no_alloc, nx_or_nxend,ncidopen,irestart,microphys_namelist)

  USE COMMASMPI_MODULE

  implicit none
  
  character(LEN=*)           :: filename
  TYPE(GRID)                 :: gd
  integer                    :: start0

  TYPE(ATTRIBUTE), pointer   :: attr
  TYPE(VARIABLE),  pointer   :: var
  integer, optional          :: no_alloc   ! if present, do not allocate grid
  logical, optional, INTENT(IN)     :: nx_or_nxend ! tells whether to use nx or nxend for allocation
  integer, optional, INTENT(IN)     :: ncidopen    ! File is already open with ncid=ncidopen, so don't open again
  integer, optional, INTENT(IN)     :: irestart    ! control over whether +1 time is used or not
  character(LEN=*), optional, intent(in) :: microphys_namelist
   
  integer :: ipconc,ichaff,inucopt,ienkf,isfcphys
  integer time, ibeg, iend
  integer, save :: ensize = 1
  integer ninit
  integer cmode
  integer n, status, k
  integer ncid, nt, record_id, records
  integer variable_id
  integer, allocatable :: tarray(:)
  logical, allocatable :: index(:)
  integer, save ::  nx, ny, nz
  character(LEN=desc_length) microphys
  character(LEN=name_length) tname
  integer start, restart
  integer XEDimID, YEDimID, ZEDimID
  
  start = start0

!-----------------------------------------------------------------------
! Open netcdf file...

  CALL STRING_LIMITS(filename, ibeg, iend)
  
  IF ( .not. present( ncidopen ) ) THEN
  
!mpi! #if defined (NC4) && defined (MPI)
#ifdef NC4
  IF ( parallelio .and. parallelio_in == 1 ) THEN
#ifdef MPI

   IF ( parallelio_type /= 2 ) THEN ! HDF5 parallel
   cmode  = ior( ior(NF90_NOWRITE,nf90_iotype),NF90_MPIIO)
   ELSE ! pnetcdf
   cmode  = ior(NF90_NOWRITE, ior(pnetcdf_flag,NF90_MPIIO))
!   status = NF90_OPEN(filename(ibeg:iend),cmode,ncid,comm=my_comm,info=my_info)
   ENDIF

    status = NF90_OPEN_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
#else
    write(0,*) 'GRID_READ_NETCDF: not compiled for MPI! set parallelio to .FALSE.'
#endif
  ELSE
    cmode  = NF90_NOWRITE
    status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
  ENDIF
!mpi! 
!mpi!   write(0,*) 'GRID_READ_NETCDF:  NETCDF FILE OPENED SUCCESSFULLY FOR PARALLEL IO', NF90_STRERROR(status)
!mpi!   write(0,*) 'GRID_READ_NETCDF:  MY_RANK=',my_rank,'MY_INFO=',my_info,'MY_COMM=',my_comm
!mpi! #else
#else
  cmode  = NF90_NOWRITE
  status = NF90_OPEN(filename(ibeg:iend),cmode,ncid)
#endif
!mpi! #endif

#ifdef MPI
  IF (status /= NF90_NOERR) THEN 
   write(0,*) 'GRID_READ_NETCDF:  Problem opening file'
   write(0,*) 'GRID_READ_NETCDF:  Filename = ', filename(ibeg:iend)
   write(0,*) 'GRID_READ_NETCDF:  MY_COMM = ', my_comm
   write(0,*) 'GRID_READ_NETCDF:  MY_INFO = ', my_info
   write(0,*) 'GRID_READ_NETCDF:  MY_RANK = ', my_rank
   write(0,*) 'GRID_READ_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
   CALL commasmpi_abort()
   stop
  ENDIF
#else
  IF (status /= NF90_NOERR) THEN 
   write(lundbg,*) 'GRID_READ_NETCDF:  Problem opening file'
   write(lundbg,*) 'GRID_READ_NETCDF:  Filename = ', filename(ibeg:iend)
   write(lundbg,*) 'GRID_READ_NETCDF:  NCID = ', ncid, NF90_STRERROR(status)
   stop
  ENDIF
#endif

  
  ELSE
    ncid = ncidopen
  ENDIF ! present(ncidopen)

  IF ( .not. (present(no_alloc) .or. allocated( gd%xyz3d%flt4d ) ) ) THEN
!-----------------------------------------------------------------------
! READ IN THE ATTRIBUTES NEEDED TO DEFINE THE GRID NAMESPACE

  status = NF90_GET_ATT(ncid, NF90_GLOBAL, 'MICROPHYS', microphys)
  
  IF ( present( microphys_namelist ) ) THEN
    IF ( trim(microphys) /= trim(microphys_namelist) ) THEN
      write(0,*) 'Warning: namelist microphys does not match restart file!'
      write(0,*) 'microphys_namelist = ',microphys_namelist
      write(0,*) 'microphys = ',microphys
      microphys = microphys_namelist
    ENDIF
  ENDIF

  IF (status /= NF90_NOERR) THEN
    write(lundbg,*) 'GRID_READ_NETCDF:  Problem reading attribute'
    write(lundbg,*) 'GRID_READ_NETCDF:  Filename = ', filename(ibeg:iend)
    write(lundbg,*) 'GRID_READ_NETCDF:  NAME = ', microphys
  ENDIF

  CALL STRING_LIMITS(microphys, ibeg, iend)

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  MICRO = ', microphys(ibeg:iend)

!
! Get value of ienkf
!
    status = NF90_INQ_VARID(ncid,'IENKF',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE IENKF ',NF90_STRERROR(status)
     ienkf = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, ienkf)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN IENKF:  ',ienkf
    ENDIF

!
! Get value of isfcphys
!
    status = NF90_INQ_VARID(ncid,'ISFCPHYS',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE isfcphys ',NF90_STRERROR(status)
     isfcphys = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, isfcphys)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN isfcphys:  ',isfcphys
    ENDIF

!
! Get value of ichaff
!
    status = NF90_INQ_VARID(ncid,'ICHAFF',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE ICHAFF ',NF90_STRERROR(status)
     ichaff = -1
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, ichaff)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN ICHAFF:  ',ichaff
    ENDIF

!
! Get value of inucopt
!
    status = NF90_INQ_VARID(ncid,'INUCOPT',variable_id)
    IF(status /= NF90_NOERR) THEN
!     print *, 'GRID_READ_NETCDF:  ERROR INQUIRE isfcphys ',NF90_STRERROR(status)
     inucopt = 0
    ELSE
     status = NF90_GET_VAR(ncid, variable_id, inucopt)
     IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN inucopt:  ',inucopt
    ENDIF


!
! Get value of ipconc
!
    status = NF90_INQ_VARID(ncid,'IPCONC',variable_id)
    IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR INQUIRE IPCONC ',NF90_STRERROR(status)

    status = NF90_GET_VAR(ncid, variable_id, ipconc)
    IF( DEBUG_IO ) write(6,*) 'READ_NCDF_VAR: READING IN IPCONC:  ',ipconc

  IF ( microphys(1:4) .eq. 'ZIEG' .or. microphys(1:3) .eq. 'ZVD' ) THEN

! check for special cases with diagnosed CCW and/or CCI
! ipconc 
    IF ( ipconc .eq. 1  .or. ipconc .eq. 0 ) THEN
    status = NF90_INQ_VARID(ncid,'CCW',variable_id)
      IF(status == NF90_NOERR) THEN
        IF ( ipconc .eq. 1 ) THEN
          ipconc = -3
        ELSE
          ipconc = -2
        ENDIF
      ENDIF
 
      IF ( ipconc .eq. 0 ) THEN
        status = NF90_INQ_VARID(ncid,'CCI',variable_id)
        IF(status == NF90_NOERR) ipconc = -1
      ENDIF
   ENDIF
  
  ENDIF


!-----------------------------------------------------------------------
! DEFINE GRID NAMESPACE

  CALL GRID_DEFINE_FROM_LIST(gd, microphys,ipconc, ichaff, inucopt, ienkf=ienkf,isfcphys=isfcphys)

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  NAMESPACE DEFINED'

! READ IN THE GRID ATTRIBUTES

  DO n = 1,size(gd%attr)

    CALL STRING_LIMITS(gd%attr(n)%name, ibeg, iend)

    IF ( parallelio .and. parallelio_type == 2 .and. parallelio_in == 1 ) THEN ! skip some values that might not be the same for all threads
      IF ( gd%attr(n)%name == 'NX' .or. gd%attr(n)%name == 'NY' .or. gd%attr(n)%name == 'TILE_INDEX' .or. &
           gd%attr(n)%name == 'MEMBER_NAME' .or. gd%attr(n)%name == 'OUTPUT_FILE_NAME' ) cycle
    ENDIF

    IF( gd%attr(n)%type .eq. 'int' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%int)
    IF( gd%attr(n)%type .eq. 'flt' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%flt)
    IF( gd%attr(n)%type .eq. 'str' ) status = NF90_GET_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name(ibeg:iend), gd%attr(n)%str)

    IF (status /= NF90_NOERR) THEN
      write(lundbg,*) 'GRID_READ_NETCDF READ_ATTRIBUTE:  Problem reading attribute ',gd%attr(n)%type
      write(lundbg,*) 'GRID_READ_NETCDF READ_ATTRIBUTE:  Filename = ', filename
      write(lundbg,*) 'GRID_READ_NETCDF READ_ATTRIBUTE:  NAME = ', gd%attr(n)%name
    ENDIF

    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  ATTRIBUTES= ', gd%attr(n)%name(ibeg:iend), &
                              gd%attr(n)%int, gd%attr(n)%flt, gd%attr(n)%str

! GET GRID ATTRIBUTES SPECIFYING DIMENSIONS

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NXEND' ) nxend = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NYEND' ) nyend = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NZEND' ) nzend = gd%attr(n)%int

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NX' ) nx = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NY' ) ny = gd%attr(n)%int
    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'NZ' ) nz = gd%attr(n)%int

    IF( gd%attr(n)%name(ibeg:iend) .EQ. 'SIZE_OF_ENSEMBLE' ) ensize = gd%attr(n)%int
    
  ENDDO
  
  

  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  ATTRIBUTES READ IN'
 
!-----------------------------------------------------------------------------------------------
! ALLOCATE SPACE FOR VARIABLES

!
!  Get values of XE,YE,ZE dimensions.  This, along with NX,NY,NZ can tell us if the file has historyoutput 0 or 1.
!

  status = nf90_inq_dimid(ncid, "XE", XEDimID)
  status = nf90_inq_dimid(ncid, "YE", YEDimID)
  status = nf90_inq_dimid(ncid, "ZE", ZEDimID)

  status = nf90_inquire_dimension(ncid, XEDimID, len = ncxe)
  status = nf90_inquire_dimension(ncid, YEDimID, len = ncye)
  status = nf90_inquire_dimension(ncid, ZEDimID, len = ncze)
  
! This test allows for the tile size to be changed, but only works with single file option (historyout = 1; parallel or round-robin)
  IF ( ( nprocj > 1 .or. nproci > 1 ) .and. ( historyinput == 1  .or. historyinput .eq. 3) ) THEN
    nx = nxt
    ny = nyt
    nz = nzt
  ENDIF
  
  IF ( present( nx_or_nxend ) ) THEN
   IF ( nx_or_nxend ) THEN
!    write(0,*) 'GRID_READ_NETCDF: ncxe,ncye,ncze = ',ncxe, ncye, ncze
    CALL GRID_ALLOCATE(gd, ncxe, ncye, ncze) ! allocate for WHOLE available domain, not just a tile
! set these for when code is compiled with MPI
    ixend = ncxe
    jyend = ncye
    kzend = ncze
   ELSE
    CALL GRID_ALLOCATE(gd, nx, ny, nz)
   ENDIF
  
  ELSE
  
    CALL GRID_ALLOCATE(gd, nx, ny, nz)

  ENDIF
  
!#else 
!  CALL GRID_ALLOCATE(gd, nx, ny, nz)
!#endif

  ENDIF ! present

!-----------------------------------------------------------------------------------------------
! INQUIRE AS TO HOW MANY TIME STEPS HAVE BEEN WRITTEN 

  status = NF90_INQUIRE(ncid, unlimitedDimID = record_id)
  
  IF(status /= NF90_NOERR) write(6,*) 'GRID_READ_NETCDF: ERROR FINDING TIME DIM ',NF90_STRERROR(status)
  
  status = NF90_INQUIRE_DIMENSION(ncid, record_id, len = records)
  
  IF(status /= NF90_NOERR) write(6,*) 'GRID_READ_NETCDF: ERROR READING TIME DIM ',NF90_STRERROR(status)

#ifdef MPI
  IF ( my_rank == 0 ) THEN
  write(0,*) "GRID_READ_NETCDF: NUMBER OF RECORDS=",records,"MY_RANK=",my_rank   !mpidebug
  ENDIF
#endif

!-----------------------------------------------------------------------------------------------
! Determine NT from input time value

  IF( start == 0 ) THEN
  
    nt = 1
  
  ELSE
  
! Allocate memory for time and index information
  
    allocate(tarray(records))
    allocate(index(records))
  
! Read in time array information
  
    status = NF90_INQ_VARID(ncid,'TIME',variable_id)
    IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR INQUIRE COORINDATE TIME ',NF90_STRERROR(status)
  
    status = NF90_GET_VAR(ncid, variable_id, tarray)
    IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR READING COORDINATE TIME ',NF90_STRERROR(status)
  
! Read in restart information if this is an ensemble run
  
    IF ( ensize .gt. 1 ) THEN
      
      status = NF90_INQ_VARID(ncid,'TRESTART',variable_id)
      IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR INQUIRE COORINDATE TRESTART ',NF90_STRERROR(status)
  
      status = NF90_GET_VAR(ncid, variable_id, restart)
      IF(status /= NF90_NOERR) print *, 'GRID_READ_NETCDF:  ERROR READING COORDINATE TRESTART ',NF90_STRERROR(status)
    
      IF ( start .eq. restart - 1  ) THEN ! checks for an extra time level at TIME+1 
        start = restart
      ENDIF
    
    ENDIF
  
! Search to find the correct time in the netcdf file
  
    index(:) = any(tarray(:) .eq. start) 

    IF( count(index) .EQ. 0 ) THEN !  .AND. start .LE. tarray(records) ) THEN
  
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*)
      write(6,*) 'GRID_READ_NETCDF: ERROR, CANNOT FIND A MATCHING TIME IN ',filename
      write(6,*)
  
      DO n = 1,records
        WRITE(6,*) 'GRID_READ_NETCDF: INPUT TIME ', start, ' DATA TIMES: ', tarray(n), index(n)
      ENDDO
  
      write(6,*)
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*) '==================== **** ERROR *** ========================'
      write(6,*)
      write(6,*) '====-= STOPPING RIGHT NOW BEFORE SOMEONE GETS HURT =========='
      write(6,*)

      CALL commasmpi_abort()
      STOP
  
    ELSE
  
! Find the first true element in the index array
  
      DO nt = 1,records
!        print *, 'GRID_READ_NETCDF:  INDEX = ',index(nt)
        IF( tarray(nt) .eq. start ) EXIT
      ENDDO

!  nt = nt + 1  !!!! don't need to add 1 here!!!
  
    ENDIF
  
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  INPUT TIME ', start
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  TIME ARRAY(NT) ', tarray(nt)
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  INDEX ARRAY(:) ', index
  
    deallocate(tarray)
    deallocate(index)
  
  ENDIF
  
  IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF:  NT = ', nt

!-----------------------------------------------------------------------
!  NetCDF:  Define variables - scalar arrays are written out as single xyz3d arrays

  DO n = 1,size(gd%var)
  
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF: READING IN VARIABLE:  ', gd%var(n)%name

    IF ( gd%var(n)%name == 'XTRA1D' ) CYCLE
    
    CALL READ_NCDF_VAR(ncid,gd%var(n),nt)
  
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_NETCDF: VARIABLE:  ', gd%var(n)%name, ' WAS READ FROM FILE'
  
  ENDDO

  IF ( ( nprocj > 1 .or. nproci > 1 ) .and. (historyinput == 1 .or. historyinput .eq. 3)  ) THEN
    IF( .not. SET_VARIABLE(gd,'NX',   nx)      ) write(6,*) 'GRID_READ_NETCDF:  Problem setting NX'
    IF( .not. SET_VARIABLE(gd,'NY',   ny)      ) write(6,*) 'GRID_READ_NETCDF:  Problem setting NY'
    IF( .not. SET_VARIABLE(gd,'NZ',   nz)      ) write(6,*) 'GRID_READ_NETCDF:  Problem setting NZ'
  ENDIF
  

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)

  IF ( .not. present(ncidopen) ) THEN
    status=NF90_CLOSE(ncid)
    IF(status /= NF90_NOERR) write(6,*) 'GRID_READ_NETCDF: Error Closing File: ', NF90_STRERROR(status)
  ENDIF

! Fill in base states:  Only need to do the first time (when allocating)
  
  IF ( .not. present(no_alloc) ) THEN
  
  DO n = 1,size(gd%var)
  
    tname = gd%var(n)%name
!    write(0,*) 'Variable name = ',tname
!    IF ( .false. ) THEN
    IF ( tname .eq. 'U' ) THEN
!    write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'UINIT')
!     write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
 !    write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO
    
    ELSEIF ( tname .eq. 'V' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'VINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
 !    write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'PI' ) THEN
!     IF ( my_rank == 0) write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'PIINIT')
!     IF ( my_rank == 0)write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     IF ( my_rank == 0) write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'KM' .or. tname .eq. 'TKE' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'KMINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'TH' ) THEN
!    IF ( my_rank == 0) write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'THINIT')
!     IF ( my_rank == 0) write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     IF ( my_rank == 0) write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'QV' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'QVINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'CCCN' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'CCNINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'CCCNUF' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'CCNUFINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO

    ELSEIF ( tname .eq. 'CPION' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'CPIONINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO
    
    ELSEIF ( tname .eq. 'CNION' ) THEN
 !   write(0,*) 'Base state for ',tname
     ninit = GET_VARIABLE_INDEX(gd, 'CNIONINIT')
 !    write(0,*) 'ninit = ',ninit
     DO k = 1,nz
     gd%var(n)%base1d(k) = gd%var(ninit)%flt1d(k)
!     write(0,*) gd%var(n)%base1d(k), gd%var(ninit)%flt1d(k)
     ENDDO
    
    ENDIF
    
  
   ENDDO
  ENDIF
  
!  STOP
  
  IF ( DEBUG_IO ) THEN
  write(lundbg,*) ''
  write(lundbg,*) 'GRID_READ_NETCDF:  FILE: ',filename(ibeg:iend),'  SUCCESSFULLY READ'
  write(lundbg,*) ''
  ENDIF

END SUBROUTINE GRID_READ_NETCDF

!===============================================================================
!
!
!
!
!
!
!
!===============================================================================

SUBROUTINE READ_NCDF_VAR(ncid,var,nt)

  USE NETCDF
#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  integer ncid, nt
  TYPE(variable) :: var

  integer variable_id
  integer, allocatable :: dims(:)
  integer dim, status, status2, n1, n2, n3, stag, istag, jstag, kstag,i
  logical DEBUG_IO2
  real, allocatable :: oned(:), twod(:,:), threed(:,:,:)
  real :: offset,scalefactor
  logical :: iscompressed

  status = NF90_INQ_VARID(ncid,var%name,variable_id)
  IF(status /= NF90_NOERR) THEN
    IF ( var%name == 'ISFCPHYS' ) THEN
      print *, 'READ_NETCDF_VAR:  NON-FATAL ERROR: INQUIRE VARIABLE ',var%name, NF90_STRERROR(status)
      print *, '---> SETTING value to 0, vartype = ',var%type
      var%int = 0
      return
    ELSE
      print *, 'READ_NETCDF_VAR:  FATAL ERROR: INQUIRE VARIABLE ',var%name, NF90_STRERROR(status)
      print *, '---> VARIABLE DOES NOT EXIST, SO SKIPPING IT'
      RETURN
#ifdef MPI
      call handle_err(-1)
#else
      call exit(-1)
#endif
    ENDIF
  ENDIF

#if defined (NC4) && defined (MPI)
   IF ( parallelio .and. parallelio_in == 1 ) THEN
     status = NF90_VAR_PAR_ACCESS(ncid,variable_id,INDCOL)
     IF(status /= NF90_NOERR) write(0,*) 'READ_NETCDF_VAR:  ERROR ALLOWING PARALLEL ACCESS',  &
                              NF90_STRERROR(status),'ON MY_RANK',my_rank
     IF( DEBUG_IO ) write(0,*) 'READ_NETCDF_VAR:  PARALLEL VARIABLE ACCESS ALLOWED'
   ENDIF
#endif

!-------------------------------------------------------------------------------
! DEBUG STATEMENTS

  IF( DEBUG_IO ) print *, var%name, var%type, var%istag, var%jstag, var%kstag

  DEBUG_IO2 = DEBUG_IO
!  DEBUG_IO2 = .true.

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, read in the scalar

  dim = var%dim + var%tdepend

  IF( dim == 0 ) THEN

    IF( DEBUG_IO2 ) write(6,*) 'READ_NCDF_VAR: READING IN VARIABLE:  ',var%name,var%type

    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_GET_VAR(ncid, variable_id, var%int)
      IF( DEBUG_IO2 ) write(6,*) 'READ_NCDF_VAR: READING IN VARIABLE:  ',var%name,' = ',var%int
    ELSE
      status = NF90_GET_VAR(ncid, variable_id, var%flt)
      IF( DEBUG_IO2 ) write(6,*) 'READ_NCDF_VAR: READING IN VARIABLE:  ',var%name,' = ',var%flt
! check whether the data have been compressed to integers
       offset = 0.0
       status2 = nf90_inquire_attribute(ncid, variable_id, "add_offset")
       IF(status2 == NF90_NOERR) THEN
         status2 = nf90_get_att(ncid, variable_id, "add_offset", offset)
         IF(status2 == NF90_NOERR) THEN
           var%flt = var%flt + offset
         ENDIF
       ENDIF
    ENDIF


    IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

!-------------------------------------------------------------------------------
! Else we are reading in a 1, 2, 3D array

  ELSE

    allocate(dims(dim))

    dims(:) = 1

    IF( var%tdepend /= 0 ) dims(dim) = nt

    IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  nt, dim = ',nt, dim


    IF( var%type(1:5) == 'rcnst' .or. var%type(1:5) == 'icnst'  ) THEN 
      
    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_GET_VAR(ncid, variable_id, var%int, dims)
      IF( DEBUG_IO2 ) write(6,*) 'READ_NCDF_VAR: READING VARIABLE:  ',var%name,' = ',var%int
    ELSE
      status = NF90_GET_VAR(ncid, variable_id, var%flt, dims)
      IF( DEBUG_IO2 ) write(6,*) 'READ_NCDF_VAR: READING VARIABLE:  ',var%name,' = ',var%flt
    ENDIF


      IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

      IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  ',var%name
    
    ENDIF

    IF( var%type(1:3) == 'x1d' .or. var%type(1:3) == 'y1d' .or. var%type(1:3) == 'z1d' ) THEN

      stag = -1
      IF( var%type(1:3) == 'x1d' .and. var%istag == 1) stag = 0
      IF( var%type(1:3) == 'y1d' .and. var%jstag == 1) stag = 0
      IF( var%type(1:3) == 'z1d' .and. var%kstag == 1) stag = 0

      n1 = size(var%flt1d) + stag - 2*var%ng !!mpidebug note the removal of ghost zones

#ifdef MPI
      IF( var%type(1:3) == 'x1d' .and. var%istag == 0 .and. ixend /= nxend) n1 = n1 + 1
      IF( var%type(1:3) == 'y1d' .and. var%jstag == 0 .and. jyend /= nyend) n1 = n1 + 1
      IF( var%type(1:3) == 'z1d' .and. var%kstag == 0 .and. kzend /= nzend) n1 = n1 + 1
#endif

      IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  ',var%name, n1

      allocate(oned(n1))

#ifdef MPI
     IF ( historyinput .eq. 1 .or. historyinput .eq. 3 ) THEN
      IF( var%type(1:3) == 'x1d' ) dims(1) = ixbeg
      IF( var%type(1:3) == 'y1d' ) dims(1) = jybeg
      IF( var%type(1:3) == 'z1d' ) dims(1) = kzbeg
     ENDIF
      status = NF90_GET_VAR(ncid, variable_id, oned, dims)
#else
      status = NF90_GET_VAR(ncid, variable_id, oned)
#endif

      IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

! check whether the data have been compressed to integers
       offset = 0.0
       status = nf90_inquire_attribute(ncid, variable_id, "add_offset")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "add_offset", offset)
       ENDIF
       
       scalefactor = 1.0
       iscompressed = .false.
       status = nf90_inquire_attribute(ncid, variable_id, "scale_factor")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "scale_factor", scalefactor)
         iscompressed = .true.
       ENDIF
       

      IF ( offset /= 0.0 .or. scalefactor /= 1.0 ) THEN
        var%flt1d(1:n1) = scalefactor*oned(:) + offset
!         DO i=1,n1
!           write(0,*) i,oned(i)
!         ENDDO

!       IF ( var%name == 'DXC') THEN
!         write(0,*) 'READ DXC: offset,scalefactor=',offset,scalefactor
!         DO i=1,n1
!           write(0,*) i,var%flt1d(i)
!         ENDDO
!       ENDIF
      ELSE
        var%flt1d(1:n1) = oned(:)
      ENDIF


      deallocate(oned)


    ENDIF


    IF( var%type(1:3) == 'x2d' .or. var%type(1:3) == 'y2d' .or. var%type(1:3) == 'z2d' ) THEN 

      stag = -1
      IF( var%type(1:3) == 'x2d' .and. var%istag == 1) stag = 0
      IF( var%type(1:3) == 'y2d' .and. var%jstag == 1) stag = 0
      IF( var%type(1:3) == 'z2d' .and. var%kstag == 1) stag = 0

      n1 = 3
      n2 = size(var%flt2d,dim=2) + stag - 2*var%ng !!mpidebug note the removal of ghost zones

#ifdef MPI
     IF ( historyinput .eq. 1 .or. historyinput .eq. 3 ) THEN
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
      IF( var%type(1:3) == 'x2d' .and. var%istag == 0 .and. ixend /= nxend) n2 = n2 + 1
      IF( var%type(1:3) == 'y2d' .and. var%jstag == 0 .and. jyend /= nyend) n2 = n2 + 1
      IF( var%type(1:3) == 'z2d' .and. var%kstag == 0 .and. kzend /= nzend) n2 = n2 + 1
#endif

      allocate(twod(n1,n2))

#ifdef MPI
      status = NF90_GET_VAR(ncid, variable_id, twod, dims)
#else
      status = NF90_GET_VAR(ncid, variable_id, twod, (/1,1/) )
#endif

! check whether the data have been compressed to integers
       offset = 0.0
       status = nf90_inquire_attribute(ncid, variable_id, "add_offset")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "add_offset", offset)
       ENDIF
       
       scalefactor = 1.0
       iscompressed = .false.
       status = nf90_inquire_attribute(ncid, variable_id, "scale_factor")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "scale_factor", scalefactor)
         iscompressed = .true.
       ENDIF

      IF ( offset /= 0.0 .or. scalefactor /= 1.0 ) THEN
        var%flt2d(1:n1,1:n2) = scalefactor*twod(:,:) + offset
      ELSE
        var%flt2d(1:n1,1:n2) = twod(:,:)
      ENDIF

      IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

      deallocate( twod )
      
      IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  ',var%name, n1, n2

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
    
      n1 = size(var%flt2d,dim=1) + istag - 2*var%ng !!mpidebug note the removal of ghost zones
      n2 = size(var%flt2d,dim=2) + jstag - 2*var%ng !!mpidebug note the removal of ghost zones

#ifdef MPI

     IF ( historyinput .eq. 1 .or. historyinput .eq. 3 ) THEN
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

      allocate(twod(n1,n2))

#ifdef MPI
      status = NF90_GET_VAR(ncid, variable_id, twod, dims)
#else
      status = NF90_GET_VAR(ncid, variable_id, twod, (/1, 1, nt/))
#endif
      IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

! check whether the data have been compressed to integers
       offset = 0.0
       status = nf90_inquire_attribute(ncid, variable_id, "add_offset")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "add_offset", offset)
       ENDIF
       
       scalefactor = 1.0
       iscompressed = .false.
       status = nf90_inquire_attribute(ncid, variable_id, "scale_factor")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "scale_factor", scalefactor)
         iscompressed = .true.
       ENDIF

      IF ( offset /= 0.0 .or. scalefactor /= 1.0 ) THEN
        var%flt2d(1:n1,1:n2) = scalefactor*twod(:,:) + offset
      ELSE
        var%flt2d(1:n1,1:n2) = twod(:,:)
      ENDIF

      deallocate(twod)

      IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  ',var%name, n1, n2

    ENDIF

    IF( var%type(1:4) == 'xyz3' ) THEN

      istag = -1
      jstag = -1
      kstag = -1
      IF( var%istag == 1) istag = 0
      IF( var%jstag == 1) jstag = 0
      IF( var%kstag == 1) kstag = 0
      n1 = size(var%flt3d,dim=1) + istag - 2*var%ng !!mpidebug note the removal of ghost zones
      n2 = size(var%flt3d,dim=2) + jstag - 2*var%ng !!mpidebug note the removal of ghost zones
      n3 = size(var%flt3d,dim=3) + kstag - 2*var%ng !!mpidebug note the removal of ghost zones


#ifdef MPI
     IF ( historyinput .eq. 1 .or. historyinput .eq. 3 ) THEN
      dims(1) = ixbeg
      dims(2) = jybeg
      dims(3) = kzbeg
     ENDIF
     
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

      IF( DEBUG_IO2 ) print *, 'READ_NCDF_VAR:  ',var%name, n1, n2, n3

      allocate(threed(n1,n2,n3))

#ifdef MPI
      status = NF90_GET_VAR(ncid, variable_id, threed, dims)
#else
      status = NF90_GET_VAR(ncid, variable_id, threed, (/ 1, 1, 1, nt /))
#endif
      IF(status /= NF90_NOERR) print *,'READ_NCDF_VAR:  Error reading variable: ', var%name, NF90_STRERROR(status)

! check whether the data have been compressed to integers
       offset = 0.0
       status = nf90_inquire_attribute(ncid, variable_id, "add_offset")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "add_offset", offset)
       ENDIF
       
       scalefactor = 1.0
       iscompressed = .false.
       status = nf90_inquire_attribute(ncid, variable_id, "scale_factor")
       IF(status == NF90_NOERR) THEN
         status = nf90_get_att(ncid, variable_id, "scale_factor", scalefactor)
         iscompressed = .true.
       ENDIF

      IF ( offset /= 0.0 .or. scalefactor /= 1.0 ) THEN
        var%flt3d(1:n1,1:n2,1:n3) = scalefactor*threed(:,:,:)+ offset
      ELSE
        var%flt3d(1:n1,1:n2,1:n3) = threed(:,:,:)
      ENDIF

      deallocate(threed)

    ENDIF

    deallocate(dims)

  ENDIF

END SUBROUTINE READ_NCDF_VAR

! ################################################
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
! ################################################
  subroutine handle_err(errcode)
!    use netcdf
    implicit none
    integer, intent(in) :: errcode
    
    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90_strerror(errcode))
       CALL commasmpi_abort()
       stop "Stopped"
    endif
  end subroutine handle_err

! ################################################
! ################################################

#if defined( MPI ) && defined ( NC4 )

SUBROUTINE GRID_DEFWRITE_NETCDFPAR(gd, file)

!  USE GRID_MODULE
  USE COMMASMPI_MODULE

  implicit none

  character(LEN=*), optional :: file
  TYPE(GRID)                 :: gd
  TYPE(ATTRIBUTE), pointer   :: attr

  integer n, status, ibeg, iend
  integer ncid, dim_id(0:8),dimvals(0:8)
  integer cmode
  character(LEN=120) filename
  
  integer mpistag
  integer variable_id
  
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and wtype

  IF( PRESENT(file) ) THEN
 
   CALL STRING_LIMITS(file, ibeg, iend)
   filename(ibeg:iend) = file(ibeg:iend)
  
  ELSE
  
   CALL GET_ATTRIBUTE(gd, 'MEMBER_NAME', attr)
   CALL STRING_LIMITS(attr%str, ibeg, iend)
   filename(ibeg:iend) = attr%str(ibeg:iend)
  
  ENDIF

! linking with netcdf4
! write(0, *) 'netcdfversion = ',netcdfversion
 IF ( netcdfversion .eq. 4 .or. parallelio ) THEN ! turn on hdf5 format to allow compression, but keep classic model
  IF ( parallelio ) THEN
   IF ( parallelio_type == 1 ) THEN ! HDF5 parallel
     cmode  = ior(NF90_NETCDF4,nf90_iotype)
    cmode = IOR(nf90_netcdf4, nf90_classic_model) 
    cmode = IOR(cmode, nf90_mpiio) 
   ELSEIF ( parallelio_type == 2 ) THEN ! pnetcdf
     cmode  = ior(ior(pnetcdf_flag,NF90_MPIIO),nf90_64bit_offset)
!       cmode = ior(IOR(nf90_CLOBBER, nf90_PNETCDF)
   ELSE
     write(0,*) 'GRID_DEFWRITE_NETCDFPAR: NON-MPI Compile! Please set the value of parallelio to FALSE'
     CALL commasmpi_abort()
   ENDIF
   IF( parallelio_type == 1 ) THEN
     IF ( DEBUG_IO )write(0,*) 'hdf5 create, rank = ',my_rank
     status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
   ELSE
    ! write(0,*) 'pnetcdf create, rank = ',my_rank
    status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
!    status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
   ENDIF
   IF (status /= NF90_NOERR) THEN
    write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  NETCDF FILE CREATE FAILED FOR PARALLEL IO', &
    NF90_STRERROR(status),'MY_RANK=',my_rank
    
    ! try again with _par version
!    cmode  = ior(NF90_PNETCDF,NF90_MPIIO)
!    write(0,*) 'pnetcdf create again, rank = ',my_rank
!    status = nf90_create(filename(ibeg:iend), cmode, ncid, comm = my_comm, info = my_info)
    !status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
    
   ENDIF

   if (status /= nf90_noerr) call handle_err(status)
  ELSE
!#else
   IF ( deflate .ne. 0 .and. deflate_level .ge. 1 ) THEN ! use netcdf4
   cmode = ior(NF90_CLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
!   cmode = ior(NF90_CLOBBER,NF90_NETCDF4)
   ELSE ! if no compression, then use the classic model
   cmode = ior(NF90_CLOBBER,NF90_CLASSIC_MODEL)
   ENDIF
!   write(0,*) 'cmode = ',cmode,NF90_CLOBBER,NF90_CLASSIC_MODEL,NF90_HDF5,NF90_FORMAT_NETCDF4_CLASSIC
   status = NF90_CREATE(filename(ibeg:iend),cmode,ncid)
   if (status /= nf90_noerr) call handle_err(status)
  ENDIF
!#endif
 !  status = NF90_CREATE(filename(ibeg:iend),ior(ior(NF90_CLOBBER,NF90_HDF5),NF90_CLASSIC_MODEL),ncid)
 ELSE ! allow version 3 as default for now until utilities all work with netcdf4/hdf5
  status = NF90_CREATE(filename(ibeg:iend),NF90_CLOBBER,ncid)
 ENDIF
! Use this version if the record size will exceed 2GB:
!
!  status = NF90_CREATE(filename(ibeg:iend),or(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)

  IF (status /= NF90_NOERR) THEN 
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  Problem opening file'
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  Filename = ', filename(ibeg:iend)
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  MY_COMM = ', my_comm
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  MY_INFO = ', my_info
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  MY_RANK = ', my_rank
   write(0,*) 'GRID_DEFWRITE_NETCDFPAR:  NCID = ', ncid, NF90_STRERROR(status)
  ENDIF

   if (status /= nf90_noerr) call handle_err(status)
  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFWRITE_NETCDFPAR: FILE CREATED'

!-----------------------------------------------------------------------
!  NetCDF:  Define dimensions

  status     = NF90_DEF_DIM(ncid, 'TIME', NF90_UNLIMITED, dim_id(0))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error unlimited dimension'


  status     = NF90_DEF_DIM(ncid, 'XC', nxend-1, dim_id(1))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nxend-1'
  status     = NF90_DEF_DIM(ncid, 'XE', nxend,   dim_id(2))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nxend'
  
  status     = NF90_DEF_DIM(ncid, 'YC', nyend-1, dim_id(3))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nyend-1'
  status     = NF90_DEF_DIM(ncid, 'YE', nyend,   dim_id(4))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nyend'
  
  status     = NF90_DEF_DIM(ncid, 'ZC', nzend-1, dim_id(5))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nzend-1'
  status     = NF90_DEF_DIM(ncid, 'ZE', nzend,   dim_id(6))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining nzend'


  IF ( nproci > 1 ) THEN
  dimvals(1) = itile
  ELSE
  dimvals(1) = itile-1
  ENDIF
  dimvals(2) = itile

  IF ( nprocj > 1 ) THEN
  dimvals(3) = jtile
  ELSE
  dimvals(3) = jtile-1
  ENDIF
  dimvals(4) = jtile

  IF ( nprock > 1 ) THEN
  dimvals(5) = ktile
  ELSE
  dimvals(5) = ktile-1
  ENDIF
  dimvals(6) = ktile
  
 

#ifdef OPENDX
  status     = NF90_DEF_DIM(ncid, 'dim2', 2, dim_id(7))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining dim2'
  status     = NF90_DEF_DIM(ncid, 'dim3', 3, dim_id(8))
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error defining dim3'
#endif

  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFWRITE_NETCDFPAR: DIMENSIONS ARE DEFINED'

!-----------------------------------------------------------------------
!  NetCDF:  Define variables - scalar arrays are written out as single xyz3d arrays

!  IF ( my_rank == 0 ) write(0,*) 'GRID_DEFWRITE_NETCDFPAR: size(gd%var) = ',size(gd%var)
  DO n = 1,size(gd%var)
  
  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFWRITE_NETCDFPAR: DEFINING VARIABLE ', gd%var(n)%name,n,size(gd%var)
  
  CALL DEFINE_NCDF_VARPAR(ncid,dim_id,gd%var(n),dimvals)

 ! status = NF90_ENDDEF(NCID)
 ! status = NF90_REDEF(NCID)
  
  ENDDO
  
  IF( DEBUG_IO ) write(lundbg,*) 'GRID_DEFWRITE_NETCDFPAR: VARIABLES DEFINED'

!-----------------------------------------------------------------------
! NetCDF:  Define global attributes
  
 ! IF ( .not. ( parallelio .and. parallelio_type == 1 ) ) THEN
  IF ( my_rank == 0 .or. historyoutput == 0 .or. parallelio_type /= 3 ) THEN
  DO n = 1,size(gd%attr)
  
   IF ( parallelio ) THEN ! skip some values that might not be the same for all threads
     IF ( gd%attr(n)%name == 'NX' .or. gd%attr(n)%name == 'NY' .or. gd%attr(n)%name == 'TILE_INDEX' .or. &
          gd%attr(n)%name == 'MEMBER_NAME' .or. gd%attr(n)%name == 'OUTPUT_FILE_NAME' ) cycle
   ENDIF
  
  IF( DEBUG_IO ) write(6,*) 'GRID_DEFWRITE_NETCDFPAR: DEFINING ATTRIBUTE ', gd%attr(n)%name,gd%attr(n)%type
  
  IF( gd%attr(n)%type == 'str') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%str)
  IF( gd%attr(n)%type == 'int') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%int)
  IF( gd%attr(n)%type == 'flt') status = NF90_PUT_ATT(ncid, NF90_GLOBAL, gd%attr(n)%name, gd%attr(n)%flt)
  
  ENDDO
  
  ENDIF
  !ENDIF

  IF( DEBUG_IO ) write(6,*) 'GRID_DEFWRITE_NETCDFPAR: ATTRIBUTES DEFINED'

!-----------------------------------------------------------------------------------------------
! End definition section

  status = NF90_ENDDEF(NCID)
  IF(status /= NF90_NOERR) print *,'GRID_DEFWRITE_NETCDFPAR:  Error ending define mode', status, NF90_STRERROR(status)

!-----------------------------------------------------------------------------------------------
! Close NETCDF file (very important!)

  status=NF90_CLOSE(NCID)
  IF(status /= NF90_NOERR) print *, 'DEFINE_NETCDF: Error Closing File: ', NF90_STRERROR(status)
  if (status /= nf90_noerr) call handle_err(status)

  IF ( DEBUG_IO ) THEN
  write(lundbg,*) 
  write(lundbg,*) 'GRID_DEFWRITE_NETCDFPAR:  FILE: ',filename(ibeg:iend), ' SUCCESSFULLY DEFINED'
  write(lundbg,*) 
  ENDIF

END SUBROUTINE GRID_DEFWRITE_NETCDFPAR

!  
!===============================================================================

SUBROUTINE DEFINE_NCDF_VARPAR(ncid,dim_id,var,dimvals)

  use mpi
  USE NETCDF 
  
  implicit none

  integer ncid, dim_id(0:8), dimvals(0:8), sizes(1:4)
  type(variable) :: var

  integer, allocatable :: dimids(:)
  integer dim, nid, status
  integer, allocatable :: chunksizes(:) , extend_increments(:)
  logical isinteger

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, write out the scalar 

  dim = var%dim + var%tdepend

  IF( dim == 0 ) THEN

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR: var%name ',var%name
    IF( var%type(1:5) == 'icnst' ) THEN
      status = NF90_DEF_VAR(ncid, var%name, NF90_INT, nid)
      isinteger = .true.
    ELSE
      status = NF90_DEF_VAR(ncid, var%name, NF90_FLOAT, nid)
      isinteger = .false.
    ENDIF
    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR: var%name,isinteger ', var%name,isinteger

    IF(status /= NF90_NOERR) write(6,*) 'DEFINE_NCDF_VARPAR:  Error creating variable: ', var%name
    if (status /= nf90_noerr) call handle_err(status)

!-------------------------------------------------------------------------------
! Else we are writing out a 1, 2, 3D array - Use a trick about the type labeling to cut down
!      on the number of conditionals

  ELSE

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR 2: var%name, dim ',var%name, dim
    allocate(dimids(dim))
    
    dimids(:) = 0

    IF( var%tdepend /= 0 ) dimids(dim) = dim_id(0)

    IF( var%type(1:1) == 'x' ) THEN               ! x1d set
      IF( var%istag == 1) THEN 
        dimids(1) = dim_id(2)
        sizes(1) = dimvals(2)
      ELSE
        dimids(1) = dim_id(1)
        sizes(1) = dimvals(1)
      ENDIF
    ENDIF

    IF( var%type(1:1) == 'y' ) THEN               ! y1d set
      IF( var%jstag == 1) THEN 
        dimids(1) = dim_id(4)
      ELSE
        dimids(1) = dim_id(3)
      ENDIF
    ENDIF

    IF( var%type(1:1) == 'z' ) THEN               ! z1d set
      IF( var%kstag == 1) THEN 
        dimids(1) = dim_id(6)
      ELSE
        dimids(1) = dim_id(5)
      ENDIF
    ENDIF

    IF( var%type(2:2) == 'y' ) THEN               ! xy2d set
      IF( var%jstag == 1) THEN 
        dimids(2) = dim_id(4)
        sizes(2) = dimvals(4)
      ELSE
        dimids(2) = dim_id(3)
        sizes(2) = dimvals(3)
      ENDIF
    ENDIF

    IF( var%type(2:2) == 'z' ) THEN               ! xz2d & yz2d set
      IF( var%kstag == 1) THEN 
        dimids(2) = dim_id(6)
      ELSE
        dimids(2) = dim_id(5)
      ENDIF
    ENDIF

    IF( var%type(2:2) == '2' ) THEN               ! x2d & y2d & z2d set
        dimids(2) = dimids(1)
        dimids(1) = dim_id(8)
    ENDIF

    IF( var%type(3:3) == 'z' ) THEN               ! xyz3d set
      IF( var%kstag == 1) THEN 
        dimids(3) = dim_id(6)
        sizes(3) = dimvals(6)
      ELSE
        dimids(3) = dim_id(5)
        sizes(3) = dimvals(5)
      ENDIF
    ENDIF

! Write out information

    IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR: var%name, dimids ',var%name, dimids

    IF( var%type(1:5) == 'icnst' ) THEN
      IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR INT 1: var%name ',var%name, ncid
      status = NF90_DEF_VAR(ncid, var%name, NF90_INT, dimids, nid)
      isinteger = .true.
      IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR INT: var%name, dimids, nid ',var%name, dimids, nid
#ifdef NC4
     IF ( netcdfversion .eq. 4 .and. parallelio_type /= 2) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, (/ 300 /) )
        if (status /= nf90_noerr) call handle_err(status)
     ENDIF
#endif
    ELSE
      status = NF90_DEF_VAR(ncid, var%name, NF90_real, dimids, nid)
      isinteger = .false.

#if defined (NC4)
      IF( netcdfversion .eq. 4 .and. var%type(1:5) == 'rcnst' ) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, (/ 300 /) )
        if (status /= nf90_noerr) call handle_err(status)
      ENDIF

     IF ( netcdfversion .eq. 4 .and. parallelio_type /= 2) THEN !{
        IF( DEBUG_IO ) write(0,*) 'DEFINE_NCDF_VARPAR NC4 setup: var%name ',var%name

      IF( var%type(3:3) == 'z' .or. var%type(1:2) == 'xy' ) THEN               ! xyz3d set
!      IF( var%type(3:3) == 'z' ) THEN               !{ xyz3d set
        allocate( chunksizes(4) )
!        allocate( extend_increments(4) ) 
        chunksizes(1:3) = sizes(1:3)
        
        IF ( var%type(3:3) == 'z' ) THEN
!        IF ( sizes(2) > 2 .and. nprock == 1 ) THEN ! not 2D and not tiled in vertical
        IF ( sizes(2) > 2 ) THEN ! not 2D and not tiled in vertical
          IF ( 4*(sizes(3)/4) == sizes(3) ) THEN
          chunksizes(3)   = sizes(3)/4
          ELSE
          chunksizes(3)   = sizes(3)/4 + 1
          ENDIF
        ENDIF
        ELSE
          chunksizes(3)   = 1
        ENDIF 
        chunksizes(4)   = 1
!        extend_increments(1:4) = chunksizes(1:4)

#ifdef MPI
        ! send chunk sizes from rank 0 so that all have the same value
         IF ( parallelio ) THEN
         CALL MPI_Bcast(chunksizes, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_error_code)
         ENDIF
#endif
        
!        write(0,*) 'DEFINE_NCDF_VARPAR sizes nonparallel = ',chunksizes(:)
!        IF ( var%type(3:3) == 'z' ) THEN
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error setting chunking for: ', var%name
!        ENDIF
        deallocate( chunksizes )
!        deallocate( extend_increments ) 

        IF ( (.not. parallelio) .or. ( var%name(1:1) == 'Q' .and. parallel_compress .and. var%type(3:3) == 'z') ) THEN   
!        IF ( (.not. parallelio) .or. ( parallel_compress .and. var%type(3:3) == 'z') ) THEN   
        write(0,*) 'set deflate for ',var%name,var%type
        status = NF90_DEF_VAR_DEFLATE(ncid, nid, shuffle, deflate, deflate_level)
!        write(0,*) 'deflate, deflate_level, shuffle = ',deflate, deflate_level, shuffle
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error setting deflate for: ', var%name
        ENDIF 
        ENDIF!}
        
     ELSEIF ( netcdfversion .eq. 4 .and. parallelio .and. parallelio_type /= 2 ) THEN ! } {
      IF( var%type(3:3) == 'z' ) THEN               ! xyz3d set
        allocate( chunksizes(4) )
!        allocate( extend_increments(4) ) 
        IF ( sizes(2) > 2 ) THEN ! not 2D
          IF ( 4*(sizes(3)/4) == sizes(3) ) THEN
          chunksizes(3)   = sizes(3)/4
          ELSE
          chunksizes(3)   = sizes(3)/4 + 1
          ENDIF
        ENDIF
!        extend_increments(1:4) = chunksizes(1:4)
        
!        write(0,*) 'DEFINE_NCDF_VARPAR: sizes parallel = ',var%name,chunksizes(:)
!        write(0,*) 'DEFINE_NCDF_VARPAR: var%kstag = ',var%kstag,dimvals(5),dimvals(6)
        status = NF90_DEF_VAR_CHUNKING(ncid, nid, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error setting chunking for: ', var%name
        deallocate( chunksizes )
!        deallocate( extend_increments ) 
      ENDIF

     ENDIF !}
#endif
    ENDIF

    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable: ', var%name
    if (status /= nf90_noerr) call handle_err(status)

    deallocate(dimids)

  ENDIF

! DEFINE ATTRIBUTES
  IF ( .true. ) THEN
  status = NF90_PUT_ATT(ncid, nid, "type", var%type)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute long_name: ', var%name

  status = NF90_PUT_ATT(ncid, nid, "long_name", var%description)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute long_name: ', var%name

  status = NF90_PUT_ATT(ncid, nid, "units",     var%unit)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute units:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "index",     var%index)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute index:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "posdef",    var%pdef)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute posdef:    ', var%name

  status = NF90_PUT_ATT(ncid, nid, "istag",      var%istag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute istag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "jstag",      var%jstag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute jstag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "kstag",      var%kstag)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute kstag:     ', var%name

  status = NF90_PUT_ATT(ncid, nid, "dyntype",   var%dyntype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute dyntype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "phytype",   var%phytype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute phytype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "buotype",   var%buotype)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute buotype:   ', var%name

  status = NF90_PUT_ATT(ncid, nid, "basedim",   var%basedim)
  IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute basedim:   ', var%name

  IF ( var%field .ne. CHAR(0) ) THEN
    status = NF90_PUT_ATT(ncid, nid, "field",   var%field)
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute field:   ', var%name
  ENDIF
  
  IF ( var%positions .ne. CHAR(0) ) THEN
    status = NF90_PUT_ATT(ncid, nid, "positions",   var%positions)
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute positions:   ', var%name
  ENDIF
  
! for CF conventions and VAPOR import
  IF ( var%name == 'XC' .or. var%name == 'XE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'X')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute axis:   ', var%name
  ENDIF

  IF ( var%name == 'YC' .or. var%name == 'YE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Y')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute axis:   ', var%name
  ENDIF

  IF ( var%name == 'ZC' .or. var%name == 'ZE' ) THEN
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Z')
    IF(status /= NF90_NOERR) print *,'DEFINE_NCDF_VARPAR:  Error creating variable attribute axis:   ', var%name
  ENDIF
  
  ENDIF ! t/f


END SUBROUTINE DEFINE_NCDF_VARPAR



#endif \* mpi *\


