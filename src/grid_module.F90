!-------------------------------------------------------------------------------
!
!
!
!
!
!
!
!-------------------------------------------------------------------------------

MODULE GRID_MODULE

 implicit none

 logical DEBUG          ; parameter(DEBUG          = .false.)
 logical FAILsoEXIT     ; parameter(FAILsoEXIT     = .true. )
 integer lun2           ; parameter(lun2           = 81     )
 integer lundbg         ; parameter(lundbg         = 0      )
 logical overwrite_file ; parameter(overwrite_file = .true. )
 logical copy_variable  ; parameter(copy_variable  = .true. )

 integer name_length    ; parameter(name_length    = 10)
 integer type_length    ; parameter(type_length    =  5)
 integer unit_length    ; parameter(unit_length    = 15)
 integer desc_length    ; parameter(desc_length    =255)

 TYPE VARIABLE

  character(LEN=name_length)   name
  character(LEN=type_length)   type
  integer                      tdepend
  integer                      dim
  integer                      ng
  integer                      istag
  integer                      jstag
  integer                      kstag
  integer                      index
  integer                      pdef            ! positive definite variable
  integer                      dyntype         ! dynamics type (flag)
  integer                      phytype         ! mixing type   (flag)
  integer                      buotype         ! whether variable used in buoyancy calculation
  character(LEN=unit_length)   unit
  character(LEN=desc_length)   description
  integer                      basedim
  character(LEN=desc_length)   field
  character(LEN=desc_length)   positions
  integer                      binindex
  integer                      numbins

  integer                   :: int  
  real                      :: flt
  real, allocatable         :: flt1d(:) 
  real, allocatable         :: flt2d(:,:) 
  real, allocatable         :: tmp3d(:,:,:) 
  real, pointer             :: flt3d(:,:,:) 

  real, allocatable         :: base1d(:) 
  real, allocatable         :: base2d(:,:) 
  real, allocatable         :: base3d(:,:,:) 

 END TYPE VARIABLE

 TYPE XYZ4D
  real, allocatable         :: flt4d(:,:,:,:) 
 END TYPE XYZ4D
 
 TYPE ATTRIBUTE
  character(LEN=20)             name
  character(LEN=3)              type
  integer                       int
  real                          flt
  character(LEN=desc_length) :: str 
 END TYPE ATTRIBUTE

 TYPE GRID

  TYPE (ATTRIBUTE), allocatable, dimension(:) :: attr
  TYPE (XYZ4D) :: xyz3d
  TYPE (VARIABLE),  allocatable, dimension(:) :: var

 END TYPE GRID

 TYPE MAXMIN
  character(LEN=name_length) name
  real                       max,  min
  integer                    imax, imin
  integer                    jmax, jmin
  integer                    kmax, kmin
 END TYPE MAXMIN
 
 INTERFACE GET_ATTRIBUTE
    module procedure GET_ATTRIBUTE_, GET_ATTRIBUTE_INT, GET_ATTRIBUTE_FLT, GET_ATTRIBUTE_STR
 END INTERFACE
 
 INTERFACE SET_ATTRIBUTE
    module procedure SET_ATTRIBUTE_, SET_ATTRIBUTE_INT, SET_ATTRIBUTE_FLT, SET_ATTRIBUTE_STR
 END INTERFACE
 
 INTERFACE SET_VARIABLE
    module procedure SET_VARIABLE_,        &
                     SET_VARIABLE_INT,     &
                     SET_VARIABLE_FLT
 END INTERFACE
    
 INTERFACE GET_VARIABLE
    module procedure GET_VARIABLE_,        &
                     GET_VARIABLE_INT,     &
                     GET_VARIABLE_FLT,     &
                     GET_VARIABLE_FLT1D,   &
                     GET_VARIABLE_FLT2D,   &
                     GET_VARIABLE_FLT3D

 END INTERFACE GET_VARIABLE

CONTAINS

#include "grid_stat.F90"
#include "grid_variable.F90"
#include "grid_attribute.F90"

END MODULE GRID_MODULE
