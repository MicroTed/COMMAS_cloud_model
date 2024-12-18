!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

   SUBROUTINE GET_ATTRIBUTE_(gd,name,dummy)

   TYPE(GRID), target :: gd  
   character(LEN=*)  name 
   character(LEN=20) l_name

   TYPE(ATTRIBUTE), pointer :: dummy

   integer l, n

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   IF( associated(dummy) ) dummy => null()

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN
     dummy => gd%attr(n)
     RETURN  
    ENDIF    

   ENDDO

   print *, 'GET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'
   stop

  RETURN
  END SUBROUTINE GET_ATTRIBUTE_
  
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_ATTRIBUTE_INT(gd,name,dummy)

   TYPE(GRID), target :: gd  
   character(LEN=*)  name 
   character(LEN=20) l_name

   integer, pointer :: dummy

   integer l, n

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   IF( associated(dummy) ) dummy => null()

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN
     dummy => gd%attr(n)%int
     RETURN  
    ENDIF    

   ENDDO

   print *, 'GET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'
   stop

  RETURN
  END SUBROUTINE GET_ATTRIBUTE_INT
  
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_ATTRIBUTE_FLT(gd,name,dummy)

   TYPE(GRID), target :: gd  
   character(LEN=*)  name 
   character(LEN=20) l_name
   real, pointer :: dummy

   integer l, n

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name
 
   IF( associated(dummy) ) dummy => null()

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN
     dummy => gd%attr(n)%flt
     RETURN  
    ENDIF    

   ENDDO

   print *, 'GET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'
   stop

  RETURN
  END SUBROUTINE GET_ATTRIBUTE_FLT
  
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GET_ATTRIBUTE_STR(gd,name,dummy)

   TYPE(GRID), pointer :: gd  
   character(LEN=*)  name 
   character(LEN=20) l_name
   character(len=desc_length), pointer :: dummy

   integer l, n

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   IF( associated(dummy) ) dummy => null()

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN
     dummy => gd%attr(n)%str
     RETURN  
    ENDIF    

   ENDDO

   print *, 'GET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'
   stop

  RETURN
  END SUBROUTINE GET_ATTRIBUTE_STR
  
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_ATTRIBUTE_(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*) name
   TYPE(ATTRIBUTE) :: dummy

   integer l, n
   character(LEN=20) l_name

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN

     gd%attr(n) = dummy
     SET_ATTRIBUTE_ = .true.
     RETURN

    ENDIF    

   ENDDO

   SET_ATTRIBUTE_ = .false.
   print *, 'SET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'

  RETURN
  END FUNCTION SET_ATTRIBUTE_
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_ATTRIBUTE_INT(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN = *) name
   integer dummy

   integer l, n
   character(LEN=20) l_name

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN

     gd%attr(n)%int = dummy
     SET_ATTRIBUTE_INT = .true.
     RETURN

    ENDIF    

   ENDDO

   SET_ATTRIBUTE_INT = .false.
   print *, 'SET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'

  RETURN
  END FUNCTION SET_ATTRIBUTE_INT
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_ATTRIBUTE_FLT(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*) name
   real dummy

   integer l, n
   character(LEN=20) l_name

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN

     gd%attr(n)%flt = dummy
     SET_ATTRIBUTE_FLT = .true.
     RETURN

    ENDIF    

   ENDDO

   SET_ATTRIBUTE_FLT = .false.
   print *, 'SET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'

  RETURN
  END FUNCTION SET_ATTRIBUTE_FLT
!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  LOGICAL FUNCTION SET_ATTRIBUTE_STR(gd,name,dummy)

   TYPE(GRID) :: gd  
   character(LEN=*) name
   character(LEN=*) dummy

   integer l, n
   character(LEN=20) l_name

   l = len(name)

   l_name(1:20) = ' '
   l_name(1:l)  = name

   DO n = 1,size(gd%attr)

    IF( gd%attr(n)%name == l_name ) THEN

     l = len(dummy)
     gd%attr(n)%str(1:l) = dummy
     SET_ATTRIBUTE_STR = .true.
     RETURN

    ENDIF    

   ENDDO

   SET_ATTRIBUTE_STR = .false.
   print *, 'SET_ATTRIBUTE:  PROBLEM, ATTRIBUTE ',name,' NOT FOUND, STOPPING'

  RETURN
  END FUNCTION SET_ATTRIBUTE_STR
