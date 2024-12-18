!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\        STRING_MODULE       ////////////////////
!
!
!
!
!===========================================================================
MODULE STRING_MODULE

 CONTAINS 
!-------------------------------------------------------------------------------
! 
!  SUBROUTINE UCASE
!    
! 
!  Purpose: 
!    To shift a character string to upper case on any processor, 
!    regardless of collating sequence. 
! 
!  Record of revisions: 
!      Date       Programmer          Description of change 
!      ====       ==========          ===================== 
!    01/09/96    S. J. Chapman        Original code 
!    12/22/05    L. J. Wicker         Adapted for NCOMMAS work 
!
!-------------------------------------------------------------------------------
 FUNCTION UCASE ( string ) 
  IMPLICIT NONE 

! Declare calling parameters: 

  CHARACTER(len=*), INTENT(IN) :: string      ! Input string 
  CHARACTER(len=LEN(string))   :: ucase       ! Function 

! Declare local variables: 
  INTEGER :: n                 ! Loop index 
  INTEGER :: length            ! Length of input string 

! Get length of string 

  length = LEN ( string ) 

! Now shift lower case letters to upper case. 

  DO n = 1, length 

   IF ( LGE(string(n:n),'a') .AND. LLE(string(n:n),'z') ) THEN 
    ucase(n:n) = ACHAR ( IACHAR ( string(n:n) ) - 32 ) 
   ELSE 
    ucase(n:n) = string(n:n) 
   END IF 

  END DO 

 END FUNCTION UCASE 

END MODULE STRING_MODULE

