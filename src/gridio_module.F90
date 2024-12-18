MODULE GRIDIO_MODULE

 USE GRID_MODULE
 USE NETCDF
 USE FILE_MODULE
 USE COMMASMPI_MODULE
#ifndef MPI
!  USE HDF5
#endif

 implicit none
 
 logical :: DEBUG_IO       = .false.
! logical DEBUG_IO       ; parameter(DEBUG_IO       = .true.)
! #if defined (NC4) && defined (MPI)
#if defined (NC4)
! integer, parameter :: NF90_INDEPENDENT = 0
! integer, parameter :: NF90_COLLECTIVE = 1
!These next two are now defined in netcdf as of version 4.1.2 (or maybe earlier),
!so comment out if you get a compile error.
! integer, parameter :: NF90_MPIIO = 2*4096, NF90_PNETCDF=0
! integer, parameter :: NF90_MPIPOSIX = 4*4096
 integer, parameter :: nf90_iotype = NF90_MPIIO !  NF90_MPIPOSIX
! integer, parameter :: nf90_iotype = NF90_MPIPOSIX
 integer, parameter :: INDCOL = NF90_COLLECTIVE
! integer, parameter :: INDCOL = NF90_INDEPENDENT
#endif

 
 character(len=155),private :: f(1000)
 integer           ,private :: numlines = 0

CONTAINS

#include "grid_create.F90"
#include "grid_netcdf.F90"
#ifndef MPI
! #include "grid_hdf5.F90"
#endif
#include "grid_binary.F90"
!include "grid_ascii.f90"

END MODULE GRIDIO_MODULE

MODULE GRIDPIO_MODULE

 USE GRID_MODULE
! #ifdef USE_PNETCDF
#if defined(USE_PNETCDF) & defined(MPI)
 USE PNETCDF
! use mpi
 USE FILE_MODULE
 USE COMMASMPI_MODULE

 implicit none
 
! include 'mpif.h'
 
 logical DEBUG_IO       ; parameter(DEBUG_IO       = .false.)
! #if defined (NC4) && defined (MPI)
#if defined (NC4)
! integer, parameter :: NF90_INDEPENDENT = 0
! integer, parameter :: NF90_COLLECTIVE = 1
!These next two are now defined in netcdf as of version 4.1.2 (or maybe earlier),
!so comment out if you get a compile error.
! integer, parameter :: NF90_MPIIO = 2*4096, NF90_PNETCDF=0
! integer, parameter :: NF90_MPIPOSIX = 4*4096
! integer, parameter :: nf90_iotype = NF90_MPIIO !  NF90_MPIPOSIX
! integer, parameter :: nf90_iotype = NF90_MPIPOSIX
! integer, parameter :: INDCOL = NF90_COLLECTIVE
! integer, parameter :: INDCOL = NF90_INDEPENDENT
#endif

 

CONTAINS

#include "grid_pnetcdf.F90"
#endif

END MODULE GRIDPIO_MODULE




!     This is part of the netCDF package.  Copyright 2006 University
!     Corporation for Atmospheric Research/Unidata.  See COPYRIGHT
!     file for conditions of use.

!     This is a very simple example which writes a 2D array of sample
!     data. To handle this in netCDF we create two shared dimensions,
!     "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is
!     part of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
      
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_nc4_wr.f90,v 1.6 2010/04/06 19:32:09 ed Exp $
#ifdef MPI
subroutine simple_xy_wr(my_rank,nproc)
  use netcdf
  use mpi
  implicit none

  character (len = *), parameter :: FILE_NAME = "simple_xyz_nc4.nc"
  integer, parameter :: NDIMS = 3
  integer, parameter :: nprocx = 2, nprocy = 2
  integer :: nproc ! = nprocx*nprocy
  integer, parameter :: nxend = 60, nyend = 60, nzend = 50
  integer, parameter :: NX = nxend/nprocx, NY = nyend/nprocx, NZ = nzend
  integer, parameter :: numvar = 14
  integer :: ncid, nid, varid, dimids(NDIMS+1), axisid(4) ! varid2,varid3,varid4,varid5,varid6
  integer :: varids(numvar)
  character(len=40) :: varnames(numvar)
  character(len=3) :: numstr
  integer :: x_dimid, y_dimid, z_dimid, time_dimid
  real :: data_out(NX, NY, NZ)
  integer :: chunks(ndims+1)
  integer :: deflate_level
  integer :: x, y, i,j,k,n, time
  integer cmode, status
  integer, parameter :: nf90_iotype = NF90_MPIIO !  NF90_MPIPOSIX
  integer, parameter :: pnetcdf_flag = NF90_PNETCDF
  integer, parameter :: INDCOL = NF90_COLLECTIVE
  logical, parameter :: DEBUG_IO = .false.
  integer :: parallelio_type = 1
  INTEGER :: my_comm, my_info
  integer p, my_rank, ierr
  integer start(NDIMS+1), count(NDIMS+1)
  real :: xc(nxend), yc(nyend), zc(nzend)
  real :: dx = 500., dy = 500, dz = 400.
  integer, parameter :: ivardef = 2


!      call MPI_Init(ierr)
!      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
!      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

          my_comm = MPI_COMM_WORLD
          my_info = MPI_INFO_NULL
    p = nprocx*nprocy
!     There must be 4 procs for this test.
      if (p .ne. nproc) then
         print *, 'This test program must be run on ',p,' processors.'
         stop 2
      endif
    
    
  ! Create some pretend data. If this wasn't an example program, we
  ! would have some real data to write, for example, model output.
  do k = 1, nz
  do j = 1, NY
     do i = 1, NX
        data_out(i, j, k) = (my_rank+1)*( (j - 1) * NY + (i - 1) )
     end do
  end do
  end do
  
  DO i = 1,nxend
    xc(i) = (i - 0.5)*dx
  ENDDO
  DO j = 1,nyend
    yc(j) = (j - 0.5)*dy
  ENDDO
  DO k = 1,nzend
    zc(k) = (k - 0.5)*dz
  ENDDO
  

   IF ( my_rank == 0 ) write(0,*) 'create file'
   IF ( parallelio_type == 1 ) THEN ! HDF5 parallel
     cmode  = ior(NF90_NETCDF4,nf90_iotype)
    cmode = IOR(nf90_netcdf4, nf90_classic_model) 
    cmode = IOR(cmode, nf90_mpiio) 
   ELSEIF ( parallelio_type == 2 ) THEN ! pnetcdf
     cmode  = ior(ior(pnetcdf_flag,NF90_MPIIO),nf90_64bit_offset)
   ENDIF
!       cmode = ior(IOR(nf90_CLOBBER, nf90_PNETCDF)

  ! Always check the return code of every netCDF function call. In
  ! this example program, wrapping netCDF calls with "call check()"
  ! makes sure that any return which is not equal to nf90_noerr (0)
  ! will print a netCDF error message and exit.

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
!  call check( nf90_create(FILE_NAME, nf90_netcdf4, ncid) )
   IF( parallelio_type == 1 ) THEN
     IF ( DEBUG_IO )write(0,*) 'hdf5 create, rank = ',my_rank
     status = nf90_create(FILE_NAME, cmode, ncid, comm = my_comm, info = my_info)
   ELSE
    ! write(0,*) 'pnetcdf create, rank = ',my_rank
    status = nf90_create(FILE_NAME, cmode, ncid, comm = my_comm, info = my_info)
!    status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
   ENDIF
   call check(status)

!   IF ( my_rank == 0 ) write(0,*) 'define dims'

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  status     = NF90_DEF_DIM(ncid, 'TIME', NF90_UNLIMITED, time_dimid)
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error unlimited dimension'

  call check( nf90_def_dim(ncid, "XC", NXend, x_dimid) )
  call check( nf90_def_dim(ncid, "YC", NYend, y_dimid) )
  call check( nf90_def_dim(ncid, "ZC", NZend, z_dimid) )


  call check( nf90_def_var(ncid, 'XC', NF90_real, (/x_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'X')
    axisid(1) = nid
  call check( nf90_def_var(ncid, 'YC', NF90_real, (/y_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Y')
    axisid(2) = nid
    
  call check( nf90_def_var(ncid, 'ZC', NF90_real, (/z_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Z')
    axisid(3) = nid

  call check( nf90_def_var(ncid, 'TIME', NF90_real, (/time_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'seconds')
    axisid(4) = nid

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ x_dimid, y_dimid, z_dimid, time_dimid/)

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer). Optional parameters chunking, shuffle,
  ! and deflate_level are used.
  chunks(1) = Nx
  chunks(2) = Ny
  chunks(3) = nz/4
  chunks(4) = 1
  deflate_level = 2
  DO n = 1,numvar
  write(numstr,'(i3.3)') n
  varnames(n) = 'data'//numstr(1:3)
  IF ( ivardef == 1 ) then
  call check( nf90_def_var(ncid, varnames(n), NF90_real, dimids, varids(n), &
       chunksizes = chunks, shuffle = .TRUE., deflate_level = deflate_level) )
  ELSE
  
   call check( nf90_def_var(ncid, varnames(n), NF90_real, dimids, varids(n)) )
   call check( NF90_DEF_VAR_CHUNKING(ncid, varids(n), NF90_CHUNKED, chunks) )
   call check( NF90_DEF_VAR_DEFLATE(ncid, varids(n),1, 1, deflate_level) )
  
  ENDIF
    status = NF90_PUT_ATT(ncid, varids(n), "units",   'count')
   ENDDO
  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )


    write(0,*) 'close file after define ',my_rank
    status=NF90_CLOSE(NCID)
    
    CALL MPI_BARRIER(my_comm, ierr)
!    stop
    
    write(0,*) 'open file after define ',my_rank
     cmode  = ior(NF90_WRITE,nf90_iotype)
    status = NF90_OPEN_PAR(file_name,cmode,my_comm,my_info,ncid)

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.

  time = 0
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(4),INDCOL)
  call check( nf90_put_var(ncid, axisid(4), time, (/1/) )) !, count=(/1/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(1),INDCOL)
  call check( nf90_put_var(ncid, axisid(1), xc, start=(/1/), count=(/nxend/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(2),INDCOL)
  call check( nf90_put_var(ncid, axisid(2), yc, start=(/1/), count=(/nyend/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(3),INDCOL)
  call check( nf90_put_var(ncid, axisid(3), zc, start=(/1/), count=(/nzend/) ) )

!     Determine what part of the variable will be written for this
!     processor. It's a checkerboard decomposition.
      count(1) = NX
      count(2) = NY
      count(3) = NZ
      count(4) = 1
      start(3) = 1
      start(4) = 1
      if (my_rank .eq. 0) then
         start(1) = 1
         start(2) = 1
      else if (my_rank .eq. 1) then
         start(1) = NX + 1
         start(2) = 1
      else if (my_rank .eq. 2) then
         start(1) = 1
         start(2) = NY + 1
      else if (my_rank .eq. 3) then
         start(1) = NX + 1
         start(2) = NY + 1
      endif

  DO n = 1,numvar

  do k = 1, nz
  do j = 1, NY
     do i = 1, NX
        data_out(i, j, k) = (float(n)/float(numvar))*(my_rank+1)*( (j - 1) * NY + (i - 1) )
     end do
  end do
  end do

   status = NF90_VAR_PAR_ACCESS(ncid,varids(n),INDCOL)

  call check( nf90_put_var(ncid, varids(n), data_out, start=start, count=count) )
  ENDDO
  
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  print *, '*** SUCCESS writing example file ', FILE_NAME, '!'

!  call MPI_Finalize(ierr)

 contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check  
 end subroutine simple_xy_wr

! ------------------------------------------

subroutine write_test(my_rank,nproc)
  use netcdf
  use mpi
  implicit none

  character (len = *), parameter :: FILE_NAME = "simple_xyz_nc4.nc"
  integer, parameter :: NDIMS = 3
  integer, parameter :: nprocx = 2, nprocy = 2
  integer :: nproc ! = nprocx*nprocy
  integer, parameter :: nxend = 60, nyend = 60, nzend = 50
  integer, parameter :: NX = nxend/nprocx, NY = nyend/nprocx, NZ = nzend
  integer, parameter :: numvar = 14
  integer :: ncid, nid, varid, dimids(NDIMS+1), axisid(4) ! varid2,varid3,varid4,varid5,varid6
  integer :: varids(numvar)
  character(len=40) :: varnames(numvar)
  character(len=3) :: numstr
  integer :: x_dimid, y_dimid, z_dimid, time_dimid
  real :: data_out(NX, NY, NZ)
  integer :: chunks(ndims+1)
  integer :: deflate_level
  integer :: x, y, i,j,k,n, time
  integer cmode, status
  integer, parameter :: nf90_iotype = NF90_MPIIO !  NF90_MPIPOSIX
  integer, parameter :: pnetcdf_flag = NF90_PNETCDF
  integer, parameter :: INDCOL = NF90_COLLECTIVE
  logical, parameter :: DEBUG_IO = .false.
  integer :: parallelio_type = 1
  INTEGER :: my_comm, my_info
  integer p, my_rank, ierr
  integer start(NDIMS+1), count(NDIMS+1)
  real :: xc(nxend), yc(nyend), zc(nzend)
  real :: dx = 500., dy = 500, dz = 400.
  integer, parameter :: ivardef = 2


!      call MPI_Init(ierr)
!      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
!      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

          my_comm = MPI_COMM_WORLD
          my_info = MPI_INFO_NULL
    p = nprocx*nprocy
!     There must be 4 procs for this test.
      if (p .ne. nproc) then
         print *, 'This test program must be run on ',p,' processors.'
         stop 2
      endif
    
    
  ! Create some pretend data. If this wasn't an example program, we
  ! would have some real data to write, for example, model output.
  do k = 1, nz
  do j = 1, NY
     do i = 1, NX
        data_out(i, j, k) = (my_rank+1)*( (j - 1) * NY + (i - 1) )
     end do
  end do
  end do
  
  DO i = 1,nxend
    xc(i) = (i - 0.5)*dx
  ENDDO
  DO j = 1,nyend
    yc(j) = (j - 0.5)*dy
  ENDDO
  DO k = 1,nzend
    zc(k) = (k - 0.5)*dz
  ENDDO
  

   IF ( my_rank == 0 ) write(0,*) 'create file'
   IF ( parallelio_type == 1 ) THEN ! HDF5 parallel
     cmode  = ior(NF90_NETCDF4,nf90_iotype)
    cmode = IOR(nf90_netcdf4, nf90_classic_model) 
    cmode = IOR(cmode, nf90_mpiio) 
   ELSEIF ( parallelio_type == 2 ) THEN ! pnetcdf
     cmode  = ior(ior(pnetcdf_flag,NF90_MPIIO),nf90_64bit_offset)
   ENDIF
!       cmode = ior(IOR(nf90_CLOBBER, nf90_PNETCDF)

  ! Always check the return code of every netCDF function call. In
  ! this example program, wrapping netCDF calls with "call check()"
  ! makes sure that any return which is not equal to nf90_noerr (0)
  ! will print a netCDF error message and exit.

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
!  call check( nf90_create(FILE_NAME, nf90_netcdf4, ncid) )
   IF( parallelio_type == 1 ) THEN
     IF ( DEBUG_IO )write(0,*) 'hdf5 create, rank = ',my_rank
     status = nf90_create(FILE_NAME, cmode, ncid, comm = my_comm, info = my_info)
   ELSE
    ! write(0,*) 'pnetcdf create, rank = ',my_rank
    status = nf90_create(FILE_NAME, cmode, ncid, comm = my_comm, info = my_info)
!    status = NF90_CREATE_PAR(filename(ibeg:iend),cmode,my_comm,my_info,ncid)
   ENDIF
   call check(status)

!   IF ( my_rank == 0 ) write(0,*) 'define dims'

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  status     = NF90_DEF_DIM(ncid, 'TIME', NF90_UNLIMITED, time_dimid)
  IF(status /= NF90_NOERR) print *,'GRID_DEFINE_NETCDF:  Error unlimited dimension'

  call check( nf90_def_dim(ncid, "XC", NXend, x_dimid) )
  call check( nf90_def_dim(ncid, "YC", NYend, y_dimid) )
  call check( nf90_def_dim(ncid, "ZC", NZend, z_dimid) )


  call check( nf90_def_var(ncid, 'XC', NF90_real, (/x_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'X')
    axisid(1) = nid
  call check( nf90_def_var(ncid, 'YC', NF90_real, (/y_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Y')
    axisid(2) = nid
    
  call check( nf90_def_var(ncid, 'ZC', NF90_real, (/z_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'meters')
    status = NF90_PUT_ATT(ncid, nid, "axis",   'Z')
    axisid(3) = nid

  call check( nf90_def_var(ncid, 'TIME', NF90_real, (/time_dimid/), nid) )
    status = NF90_PUT_ATT(ncid, nid, "units",   'seconds')
    axisid(4) = nid

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimids =  (/ x_dimid, y_dimid, z_dimid, time_dimid/)

  ! Define the variable. The type of the variable in this case is
  ! NF90_INT (4-byte integer). Optional parameters chunking, shuffle,
  ! and deflate_level are used.
  chunks(1) = Nx
  chunks(2) = Ny
  chunks(3) = nz/4
  chunks(4) = 1
  deflate_level = 2
  DO n = 1,numvar
  write(numstr,'(i3.3)') n
  varnames(n) = 'data'//numstr(1:3)
  IF ( ivardef == 1 ) then
  call check( nf90_def_var(ncid, varnames(n), NF90_real, dimids, varids(n), &
       chunksizes = chunks, shuffle = .TRUE., deflate_level = deflate_level) )
  ELSE
  
   call check( nf90_def_var(ncid, varnames(n), NF90_real, dimids, varids(n)) )
   call check( NF90_DEF_VAR_CHUNKING(ncid, varids(n), NF90_CHUNKED, chunks) )
   call check( NF90_DEF_VAR_DEFLATE(ncid, varids(n),1, 1, deflate_level) )
  
  ENDIF
    status = NF90_PUT_ATT(ncid, varids(n), "units",   'count')
   ENDDO
  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )


    write(0,*) 'close file after define ',my_rank
    status=NF90_CLOSE(NCID)
    
    CALL MPI_BARRIER(my_comm, ierr)
!    stop
    
    write(0,*) 'open file after define ',my_rank
     cmode  = ior(NF90_WRITE,nf90_iotype)
    status = NF90_OPEN_PAR(file_name,cmode,my_comm,my_info,ncid)

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.

  time = 0
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(4),INDCOL)
  call check( nf90_put_var(ncid, axisid(4), time, (/1/) )) !, count=(/1/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(1),INDCOL)
  call check( nf90_put_var(ncid, axisid(1), xc, start=(/1/), count=(/nxend/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(2),INDCOL)
  call check( nf90_put_var(ncid, axisid(2), yc, start=(/1/), count=(/nyend/) ) )
  status = NF90_VAR_PAR_ACCESS(ncid,axisid(3),INDCOL)
  call check( nf90_put_var(ncid, axisid(3), zc, start=(/1/), count=(/nzend/) ) )

!     Determine what part of the variable will be written for this
!     processor. It's a checkerboard decomposition.
      count(1) = NX
      count(2) = NY
      count(3) = NZ
      count(4) = 1
      start(3) = 1
      start(4) = 1
      if (my_rank .eq. 0) then
         start(1) = 1
         start(2) = 1
      else if (my_rank .eq. 1) then
         start(1) = NX + 1
         start(2) = 1
      else if (my_rank .eq. 2) then
         start(1) = 1
         start(2) = NY + 1
      else if (my_rank .eq. 3) then
         start(1) = NX + 1
         start(2) = NY + 1
      endif

  DO n = 1,numvar

  do k = 1, nz
  do j = 1, NY
     do i = 1, NX
        data_out(i, j, k) = (float(n)/float(numvar))*(my_rank+1)*( (j - 1) * NY + (i - 1) )
     end do
  end do
  end do

   status = NF90_VAR_PAR_ACCESS(ncid,varids(n),INDCOL)

  call check( nf90_put_var(ncid, varids(n), data_out, start=start, count=count) )
  ENDDO
  
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  print *, '*** SUCCESS writing example file ', FILE_NAME, '!'

!  call MPI_Finalize(ierr)

 contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check  
 end subroutine write_test
#endif
