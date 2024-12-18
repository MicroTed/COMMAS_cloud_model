!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_VARIABLE_(gd,name,dummy,copy)

   implicit none
   
   TYPE(GRID), target      :: gd  
   TYPE(VARIABLE), pointer :: dummy
   character(LEN=*)        :: name 
   logical, optional       :: copy
 
   integer n, n1, n2, n3
   integer nng,k
 
   n = GET_VARIABLE_INDEX(gd,name)

   IF( n == -1 ) THEN

     print *, 'GET_VARIABLE_:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE
 
    IF( associated(dummy) ) dummy => null()

    IF( present(copy) ) THEN
  
      dummy%name        = gd%var(n)%name
      dummy%tdepend     = gd%var(n)%tdepend
      dummy%type        = gd%var(n)%type
      dummy%dim         = gd%var(n)%dim
      dummy%ng          = gd%var(n)%ng
      dummy%index       = gd%var(n)%index
      dummy%istag       = gd%var(n)%istag
      dummy%jstag       = gd%var(n)%jstag
      dummy%kstag       = gd%var(n)%kstag
      dummy%pdef        = gd%var(n)%pdef
      dummy%dyntype     = gd%var(n)%dyntype
      dummy%phytype     = gd%var(n)%phytype
      dummy%buotype     = gd%var(n)%buotype
      dummy%unit        = gd%var(n)%unit
      dummy%description = gd%var(n)%description
   
      IF( gd%var(n)%type      == 'icnst' ) dummy%int = gd%var(n)%int

      IF( gd%var(n)%type      == 'rcnst' ) dummy%flt = gd%var(n)%flt

      nng = gd%var(n)%ng

      IF( gd%var(n)%type(2:3) == '1d'   ) THEN
        n1 = size(gd%var(n)%flt1d,dim=1) - 2*nng
        allocate(dummy%flt1d(-nng+1:n1+nng))
        dummy%flt1d = gd%var(n)%flt1d
      ENDIF


      IF( gd%var(n)%type(3:4) == '2d'  ) THEN
        n1 = size(gd%var(n)%flt2d,dim=1) - 2*nng
        n2 = size(gd%var(n)%flt2d,dim=2) - 2*nng
        allocate(dummy%flt2d(-nng+1:n1+nng,-nng+1:n2+nng))
        dummy%flt2d = gd%var(n)%flt2d
      ENDIF

      IF( gd%var(n)%type(4:5) == '3d' ) THEN
        n1 = size(gd%var(n)%flt3d,dim=1) - 2*nng
        n2 = size(gd%var(n)%flt3d,dim=2) - 2*nng
        n3 = size(gd%var(n)%flt3d,dim=3) - 2*nng
        allocate(dummy%tmp3d(-nng+1:n1+nng,-nng+1:n2+nng,-nng+1:n3+nng))
        dummy%tmp3d = gd%var(n)%flt3d
        dummy%flt3d => gd%var(n)%tmp3d
      ENDIF

     ELSE

      dummy => gd%var(n)
 
     ENDIF

   ENDIF  
 
  RETURN
  END SUBROUTINE GET_VARIABLE_
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_VARIABLE_INT(gd,name,dummy)

   TYPE(GRID)         :: gd  
   integer            :: dummy
   character(LEN=*)   :: name 

   integer n

   n = GET_VARIABLE_INDEX(gd,name)

   IF( n == -1 ) THEN

     print *, 'GET_VARIABLE_INT:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE

     dummy = gd%var(n)%int

   ENDIF    

  RETURN
  END SUBROUTINE GET_VARIABLE_INT

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_VARIABLE_FLT(gd,name,dummy)

   implicit none
   TYPE(GRID)         :: gd  
   real               :: dummy
   character(LEN=*)   :: name 

   integer n

   n = GET_VARIABLE_INDEX(gd,name)

   IF( n == -1 ) THEN
 
     print *, 'GET_VARIABLE_FLT:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE

     dummy = gd%var(n)%flt

   ENDIF    

  RETURN
  END SUBROUTINE GET_VARIABLE_FLT

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT1D(gd,name,dummy,copy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:)
   character(LEN=*)   :: name 
   logical, optional  :: copy
   
   integer n, dim, nng
   
   n = GET_VARIABLE_INDEX(gd,name)
   nng = gd%var(n)%ng

   IF( n == -1 ) THEN
 
     print *, 'GET_VARIABLE_FLD1D:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE
   
     IF( associated(dummy) ) dummy => null()

     IF( present(copy) ) THEN

      dim = size(gd%var(n)%flt1d) - 2*nng
      allocate(dummy(-nng+1:dim+nng))
      dummy = gd%var(n)%flt1d

     ELSE
      
      dummy => gd%var(n)%flt1d

     ENDIF

   ENDIF    
   
  RETURN
  END SUBROUTINE GET_VARIABLE_FLT1D

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT2D(gd,name,dummy,copy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:,:)
   character(LEN=*)   :: name 
   logical, optional  :: copy
   
   integer n, dims(2), nng
   
   n = GET_VARIABLE_INDEX(gd,name)
   nng = gd%var(n)%ng

   IF( n == -1 ) THEN

     print *, 'GET_VARIABLE_FLD2D:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE
  
     IF( associated(dummy) ) dummy => null()

     IF( present(copy) ) THEN

      dims    = shape(gd%var(n)%flt2d)
      dims(1) = dims(1) - 2*nng
      dims(2) = dims(2) - 2*nng
      allocate(dummy(-nng+1:dims(1)+nng,-nng+1:dims(2)+nng))
      dummy = gd%var(n)%flt2d

     ELSE
     
      dummy => gd%var(n)%flt2d

     ENDIF

     RETURN

   ENDIF    
   
  RETURN
  END SUBROUTINE GET_VARIABLE_FLT2D

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT3D(gd,name,dummy,copy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:,:,:)
   character(LEN=*)   :: name 
   logical, optional  :: copy
   
   integer n, dims(3), nng

   n = GET_VARIABLE_INDEX(gd,name)
   nng = gd%var(n)%ng

   IF( n == -1 ) THEN

     print *, 'GET_VARIABLE_FLD2D:  PROBLEM, VARIABLE ',name,' NOT FOUND'

   ELSE
  
     IF( associated(dummy) ) dummy => null()
   
     IF( PRESENT(copy) ) THEN

      dims = shape(gd%var(n)%flt3d)
      dims(1) = dims(1) - 2*nng
      dims(2) = dims(2) - 2*nng
      dims(3) = dims(3) - 2*nng
      allocate(dummy(-nng+1:dims(1)+nng,-nng+1:dims(2)+nng,-nng+1:dims(3)+nng))
      dummy = gd%var(n)%flt3d

     ELSE

      dummy => gd%var(n)%flt3d

     ENDIF

   ENDIF    
   
  RETURN
  END SUBROUTINE GET_VARIABLE_FLT3D

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_INT_COPY(gd,name,dummy)
   
   TYPE(GRID)         :: gd  
   integer            :: dummy
   character(LEN=*)   :: name 

   CALL GET_VARIABLE_INT(gd, name, dummy )

  END SUBROUTINE GET_VARIABLE_INT_COPY
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT_COPY(gd,name,dummy)
   
   TYPE(GRID)         :: gd  
   real               :: dummy
   character(LEN=*)   :: name 

   CALL GET_VARIABLE_FLT(gd, name, dummy )

  END SUBROUTINE GET_VARIABLE_FLT_COPY
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT1D_COPY(gd,name,dummy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:)
   character(LEN=*)   :: name 

   CALL GET_VARIABLE_FLT1D(gd, name, dummy, .true. )

  END SUBROUTINE GET_VARIABLE_FLT1D_COPY
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT2D_COPY(gd,name,dummy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:,:)
   character(LEN=*)   :: name 

   CALL GET_VARIABLE_FLT2D(gd, name, dummy, .true. )

  END SUBROUTINE GET_VARIABLE_FLT2D_COPY
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------
   
  SUBROUTINE GET_VARIABLE_FLT3D_COPY(gd,name,dummy)
   
   TYPE(GRID), target :: gd  
   real, pointer      :: dummy(:,:,:)
   character(LEN=*)   :: name 

   CALL GET_VARIABLE_FLT3D(gd, name, dummy, .true. )

  END SUBROUTINE GET_VARIABLE_FLT3D_COPY
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_VARIABLE_COPY(gd,name,dummy,copy)
 
   TYPE(GRID), target      :: gd  
   TYPE(VARIABLE), target  :: dummy
   character(LEN=*)        :: name 
   logical, optional       :: copy
 
   integer l, n, n1, n2, n3, nng
   character(LEN=name_length) l_name
 
   l = len(name)
 
   l_name(1:name_length) = ' '
   l_name(1:l)           = name
 
   DO n = 1,size(gd%var)
 
    IF( gd%var(n)%name == l_name ) THEN

      dummy%name        = gd%var(n)%name
      dummy%tdepend     = gd%var(n)%tdepend
      dummy%type        = gd%var(n)%type
      dummy%dim         = gd%var(n)%dim
      dummy%ng          = gd%var(n)%ng
      dummy%index       = gd%var(n)%index
      dummy%istag       = gd%var(n)%istag
      dummy%jstag       = gd%var(n)%jstag
      dummy%kstag       = gd%var(n)%kstag
      dummy%pdef        = gd%var(n)%pdef
      dummy%dyntype     = gd%var(n)%dyntype
      dummy%phytype     = gd%var(n)%phytype
      dummy%buotype     = gd%var(n)%buotype
      dummy%unit        = gd%var(n)%unit
      dummy%description = gd%var(n)%description
   
      IF( gd%var(n)%type      == 'icnst' ) dummy%int = gd%var(n)%int

      IF( gd%var(n)%type      == 'rcnst' ) dummy%flt = gd%var(n)%flt

      nng = gd%var(n)%ng

      IF( gd%var(n)%type(2:3) == '1d'   ) THEN
        n1 = size(gd%var(n)%flt1d,dim=1) - 2*nng
        IF( allocated(dummy%flt1d) ) deallocate(dummy%flt1d)
        allocate(dummy%flt1d(-nng+1:n1+nng))
        dummy%flt1d = gd%var(n)%flt1d
      ENDIF

      IF( gd%var(n)%type(3:4) == '2d'  ) THEN
        n1 = size(gd%var(n)%flt2d,dim=1) - 2*nng
        n2 = size(gd%var(n)%flt2d,dim=2) - 2*nng
        IF( allocated(dummy%flt2d) ) deallocate(dummy%flt2d)
        allocate(dummy%flt2d(-nng+1:n1+nng,-nng+1:n2+nng))
        dummy%flt2d = gd%var(n)%flt2d
      ENDIF

      IF( gd%var(n)%type(4:5) == '3d' ) THEN
        n1 = size(gd%var(n)%flt3d,dim=1) - 2*nng
        n2 = size(gd%var(n)%flt3d,dim=2) - 2*nng
        n3 = size(gd%var(n)%flt3d,dim=3) - 2*nng
        IF( allocated(dummy%tmp3d) ) deallocate(dummy%tmp3d)
        allocate(dummy%tmp3d(-nng+1:n1+nng,-nng+1:n2+nng,-nng+1:n3+nng))
        dummy%tmp3d(:,:,:) = gd%var(n)%flt3d(:,:,:)
        dummy%flt3d => dummy%tmp3d !(-nng+1:n1+nng,-nng+1:n2+nng,-nng+1:n3+nng)
!        dummy%flt3d => dummy%tmp3d(-nng+1:n1+nng,-nng+1:n2+nng,-nng+1:n3+nng)
      ENDIF

     RETURN  
    ENDIF  
 
   ENDDO
 
   print *, 'GET_VARIABLE:  PROBLEM, VARIABLE ',name,' NOT FOUND, STOPPING'
   stop
 
  RETURN
  END SUBROUTINE GET_VARIABLE_COPY
!-------------------------------------------------------------------------------
!
!
!
!
!
!
!
!-------------------------------------------------------------------------------

  SUBROUTINE GET_VARIABLE_LIST(gd,type,list)

   TYPE(GRID) :: gd
   character(LEN=*) :: type
   character(LEN=name_length), pointer :: list(:)
   logical, allocatable :: truefalse(:)

   integer l, ntrue, nvar

   l = len(type)
   nvar = size(gd%var)

   allocate(truefalse(nvar))

   truefalse(:) = merge(.true., .false., gd%var(:)%type(1:l) == type(1:l) )

   ntrue = count( truefalse(:) .eqv. .true. )

   allocate(list(ntrue))

   IF( associated(list) ) list => null()

   ntrue = size( pack(gd%var(:)%name, truefalse(:)) )

   allocate(list(ntrue))

   list = pack(gd%var(:)%name, truefalse(:))

  RETURN
  END SUBROUTINE GET_VARIABLE_LIST


!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_VARIABLE_(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*)  name 
   character(LEN=name_length) l_name
   TYPE(VARIABLE) :: dummy

   integer l, n

   l = len(name)

   l_name(1:name_length) = ' '
   l_name(1:l)           = name

   DO n = 1,size(gd%var)

    IF( gd%var(n)%name == l_name ) THEN
     gd%var(n) = dummy
     SET_VARIABLE_ = .true.
     RETURN  
    ENDIF    

   ENDDO

   IF( FAILsoEXIT ) THEN
     print *, 'SET_VARIABLE:  PROBLEM, VARIABLE ',name,' NOT FOUND, STOPPING'
     stop
   ELSE
     SET_VARIABLE_ = .false.
   ENDIF

  RETURN
  END FUNCTION SET_VARIABLE_

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_VARIABLE_INT(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*)  name 
   character(LEN=name_length) l_name
   integer :: dummy

   integer l, n

   l = len(name)

   l_name(1:name_length) = ' '
   l_name(1:l)           = name

   DO n = 1,size(gd%var)

    IF( gd%var(n)%name == l_name ) THEN
     gd%var(n)%int = dummy
     SET_VARIABLE_INT = .true.
     RETURN  
    ENDIF    

   ENDDO

   IF( FAILsoEXIT ) THEN
     print *, 'SET_VARIABLE:  PROBLEM, VARIABLE ',name,' NOT FOUND, STOPPING'
     stop
   ELSE
     SET_VARIABLE_INT = .false.
   ENDIF

  RETURN
  END FUNCTION SET_VARIABLE_INT

!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_VARIABLE_FLT(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*)  name 
   character(LEN=name_length) l_name
   real :: dummy

   integer l, n

   l = len(name)

   l_name(1:name_length) = ' '
   l_name(1:l)           = name

   DO n = 1,size(gd%var)

    IF( gd%var(n)%name == l_name ) THEN
     gd%var(n)%flt = dummy
     SET_VARIABLE_FLT = .true.
     RETURN  
    ENDIF    

   ENDDO

   IF( FAILsoEXIT ) THEN
     print *, 'SET_VARIABLE:  PROBLEM, VARIABLE ',name,' NOT FOUND, STOPPING'
     stop
   ELSE
     SET_VARIABLE_FLT = .false.
   ENDIF

  RETURN
  END FUNCTION SET_VARIABLE_FLT
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  INTEGER FUNCTION GET_VARIABLE_INDEX(gd,name,nofail)

   TYPE(GRID) :: gd  
   character(LEN=*)  name 
   integer, optional       :: nofail

   integer l, n
   character(LEN=name_length) l_name

   l = len(name)

   l_name(1:name_length) = ' '
   l_name(1:l)           = name

   DO n = 1,size(gd%var)

    IF( gd%var(n)%name == l_name ) THEN
     GET_VARIABLE_INDEX = n
     RETURN  
    ENDIF

   ENDDO

   IF( FAILsoEXIT .and. .not. PRESENT( nofail )  ) THEN
     print *, 'GET_VARIABLE_INDEX:  PROBLEM, VARIABLE ',name,' NOT FOUND, STOPPING'
     stop
   ELSE
     GET_VARIABLE_INDEX = -1
   ENDIF

  RETURN
  END FUNCTION GET_VARIABLE_INDEX
