!===============================================================================
! 
!
!  
!  
!
!  
!  
!===============================================================================

SUBROUTINE GRID_DEFINE_HDF5(gd, file)

!  USE GRID_MODULE
#ifdef MPI
  USE COMMASMPI_MODULE
#endif

  implicit none

  character(LEN=*), optional :: file
  TYPE(GRID)                 :: gd
  TYPE(ATTRIBUTE), pointer   :: attr

  integer n, status, ibeg, iend, error
  integer ncid, dim_id(0:8)
  integer cmode
  character(LEN=120) filename
  INTEGER(HID_T) :: file_id                            ! File identifier
  
  integer mpistag
  
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and wtype

  IF( PRESENT(file) ) THEN
 
   CALL STRING_LIMITS(file, ibeg, iend)
   filename(ibeg:iend) = file(ibeg:iend)
  
  ELSE
  
   CALL GET_ATTRIBUTE(gd, 'MEMBER_NAME', attr)
   CALL STRING_LIMITS(attr%str, ibeg, iend)
   filename(ibeg:iend) = attr%str(ibeg:iend-2)//'h5'
  
  ENDIF
  
!
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)
     !
     ! Create a new file using default properties.
     ! 
     CALL h5fcreate_f(file, H5F_ACC_TRUNC_F, file_id, error)

END SUBROUTINE GRID_DEFINE_HDF5
