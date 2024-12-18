! F2KCLI : Fortran 200x Command Line Interface
! copyright Interactive Software Services Ltd. 2002
! For conditions of use see manual.txt
!
! Platform    : Mac OS/X
! Compiler    : Absoft Pro Fortran
! To compile  : f95 -c f2kcli.f90
! Implementer : Lawson B. Wakefield, I.S.S. Ltd.
! Date        : June 2002
!
!!!!!!!!!!!
! CHANGE:  LJW, Jan 2006:  Changed name to CLINE_MODULE cause I can never remember F2KCLI
!!!!!!!!!!!

  MODULE CLINE_MODULE

#ifdef NAG
      USE F90_UNIX_ENV, ONLY: IARGC, GETARG
#endif
!
      CONTAINS
!
      SUBROUTINE GET_COMMAND(COMMAND,LENGTH,STATUS)
!
! Description. Returns the entire command by which the program was
!   invoked.
!
! Class. Subroutine.
!
! Arguments.
! COMMAND (optional) shall be scalar and of type default character.
!   It is an INTENT(OUT) argument. It is assigned the entire command
!   by which the program was invoked. If the command cannot be
!   determined, COMMAND is assigned all blanks.
! LENGTH (optional) shall be scalar and of type default integer. It is
!   an INTENT(OUT) argument. It is assigned the significant length
!   of the command by which the program was invoked. The significant
!   length may include trailing blanks if the processor allows commands
!   with significant trailing blanks. This length does not consider any
!   possible truncation or padding in assigning the command to the
!   COMMAND argument; in fact the COMMAND argument need not even be
!   present. If the command length cannot be determined, a length of
!   0 is assigned.
! STATUS (optional) shall be scalar and of type default integer. It is
!   an INTENT(OUT) argument. It is assigned the value 0 if the
!   command retrieval is sucessful. It is assigned a processor-dependent
!   non-zero value if the command retrieval fails.
!
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: COMMAND
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
!
      INTEGER                   :: IARG,NARG,IPOS
      INTEGER            , SAVE :: LENARG
      CHARACTER(LEN=2000), SAVE :: ARGSTR
      LOGICAL            , SAVE :: GETCMD = .TRUE.
!
! Reconstruct the command line from its constituent parts.
! This may not be the original command line.
!
      IF (GETCMD) THEN
          NARG = command_argument_count()
          IF (NARG > 0) THEN
              IPOS = 1
              DO IARG = 1,NARG
                CALL GETARG(IARG,ARGSTR(IPOS:))
                LENARG = LEN_TRIM(ARGSTR)
                IPOS   = LENARG + 2
                IF (IPOS > LEN(ARGSTR)) EXIT
              END DO
          ELSE
              ARGSTR = ' '
              LENARG = 0
          ENDIF
          GETCMD = .FALSE.
      ENDIF
      IF (PRESENT(COMMAND)) COMMAND = ARGSTR
      IF (PRESENT(LENGTH))  LENGTH  = LENARG
      IF (PRESENT(STATUS))  STATUS  = 0
      RETURN
      END SUBROUTINE GET_COMMAND
!
      INTEGER FUNCTION COMMAND_ARGUMENT_COUNT1()
!
! Description. Returns the number of command arguments.
!
! Class. Inquiry function
!
! Arguments. None.
!
! Result Characteristics. Scalar default integer.
!
! Result Value. The result value is equal to the number of command
!   arguments available. If there are no command arguments available
!   or if the processor does not support command arguments, then
!   the result value is 0. If the processor has a concept of a command
!   name, the command name does not count as one of the command
!   arguments.
!
      COMMAND_ARGUMENT_COUNT1 = command_argument_count()
      RETURN
      END FUNCTION COMMAND_ARGUMENT_COUNT1
!
      SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
!
! Description. Returns a command argument.
!
! Class. Subroutine.
!
! Arguments.
! NUMBER shall be scalar and of type default integer. It is an
!   INTENT(IN) argument. It specifies the number of the command
!   argument that the other arguments give information about. Useful
!   values of NUMBER are those between 0 and the argument count
!   returned by the COMMAND_ARGUMENT_COUNT intrinsic.
!   Other values are allowed, but will result in error status return
!   (see below).  Command argument 0 is defined to be the command
!   name by which the program was invoked if the processor has such
!   a concept. It is allowed to call the GET_COMMAND_ARGUMENT
!   procedure for command argument number 0, even if the processor
!   does not define command names or other command arguments.
!   The remaining command arguments are numbered consecutively from
!   1 to the argument count in an order determined by the processor.
! VALUE (optional) shall be scalar and of type default character.
!   It is an INTENT(OUT) argument. It is assigned the value of the
!   command argument specified by NUMBER. If the command argument value
!   cannot be determined, VALUE is assigned all blanks.
! LENGTH (optional) shall be scalar and of type default integer.
!   It is an INTENT(OUT) argument. It is assigned the significant length
!   of the command argument specified by NUMBER. The significant
!   length may include trailing blanks if the processor allows command
!   arguments with significant trailing blanks. This length does not
!   consider any possible truncation or padding in assigning the
!   command argument value to the VALUE argument; in fact the
!   VALUE argument need not even be present. If the command
!   argument length cannot be determined, a length of 0 is assigned.
! STATUS (optional) shall be scalar and of type default integer.
!   It is an INTENT(OUT) argument. It is assigned the value 0 if
!   the argument retrieval is sucessful. It is assigned a
!   processor-dependent non-zero value if the argument retrieval fails.
!
! NOTE
!   One possible reason for failure is that NUMBER is negative or
!   greater than COMMAND_ARGUMENT_COUNT().
!
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
!
!  A temporary variable for the rare case case where LENGTH is
!  specified but VALUE is not. An arbitrary maximum argument length
!  of 1000 characters should cover virtually all situations.
!
      CHARACTER(LEN=1000) :: TMPVAL
!
! Possible error codes:
! 1 = Argument number is less than minimum
! 2 = Argument number exceeds maximum
!
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      ELSE IF (NUMBER > command_argument_count()) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 2
          RETURN
      END IF
!
! Get the argument if VALUE is present
!
      IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
!
! As under Unix, the LENGTH option is probably fairly pointless here,
! but LEN_TRIM is used to ensure at least some sort of meaningful result.
!
      IF (PRESENT(LENGTH)) THEN
          IF (PRESENT(VALUE)) THEN
              LENGTH = LEN_TRIM(VALUE)
          ELSE
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          END IF
      END IF
!
! Since GETARG does not return a result code, assume success
!
      IF (PRESENT(STATUS)) STATUS = 0
      RETURN
      END SUBROUTINE GET_COMMAND_ARGUMENT
!
 END MODULE CLINE_MODULE
 
       
!-----------------------------------------------------------------------------
!
!  SUBROUTINE GETRUNPARAM
!
!     A "generic" command line parser for the cloud models
!
!     tags  (char array) are the list of command line flags (like -dt or -tsave)
!     desc  (char array) are the list of descriptors for the tags (sec/float)
!     param (char array) are the values of the parameters returned, IN CHAR FORM
!
!     tags(0) holds the name of the program (like 'cloud2d')
!-----------------------------------------------------------------------------

      SUBROUTINE GETRUNPARAM(tags, desc, param, ntags)

      implicit none

      integer ntags
      character( LEN = * ) tags(0:ntags), desc(ntags), param(ntags)

      integer maxarg ; parameter( maxarg = 100 )
      character( LEN = 20 ) args(maxarg), farg1(maxarg), farg2(maxarg), arg1, arg2
      integer numargs, n, nn  
      integer iargc, ifile, count, ios

      external iargc

      data ifile / 0 /
      
      logical debug ; parameter( debug = .false. )
      
!-----------------------------------------------------------------------------
! Value to look for errors in the specification

      character( LEN = 5 ) perror ;  parameter( perror = "-999." )

!-----------------------------------------------------------------------------
! Read command line args

      numargs = command_argument_count()

      IF( numargs .eq. 0 ) THEN

         print *, '-----------------------------------------------------------'
         print *, ''
         print *, ' NO COMMAND LINE ARGUEMENTS SUPPLIED'
         print *, ''

      ELSE

       DO n = 1,numargs

        call getarg(n, args(n))

        IF( args(n) .eq. "-help") THEN

         print *, '-----------------------------------------------------------'
         print *, ''
         print *, 'Usage:  ',tags(0)

         DO nn = 1,ntags
          write(6,101) tags(nn),desc(nn)
         ENDDO
101      FORMAT(9x,a10,2x,a40)

         print *, ''
         print *, '-----------------------------------------------------------'

         call exit(1)

        ENDIF

       ENDDO

      ENDIF

!-----------------------------------------------------------------------------
! Parse the arguments

       DO nn = 1,ntags

         IF( ANY(args(:) == tags(nn)) ) THEN

          DO n = 1,numargs
 
           IF( args(n) .EQ. tags(nn) ) THEN
             param(nn) = args(n+1)
           ENDIF

          ENDDO

         ENDIF

       ENDDO

      RETURN
      END SUBROUTINE GETRUNPARAM

!-----------------------------------------------------------------------------
      REAL FUNCTION CHARTOFLOAT(string)

      IMPLICIT NONE

      integer k, l, m, i, n 
      character( LEN = * ) string
      real sign

      l = len(trim(string))
      m = scan(trim(string),'.')

      sign = 1.0
      IF( string(1:1) .EQ. "-" ) THEN
       sign = -1.0
       string = string(2:l)
       l = l - 1
       m = m - 1
      ENDIF

      CHARTOFLOAT = 0.0

! Parse each character (so we dont have to know the exact format)

      IF( m .eq. 0 ) THEN           ! CASE WHERE THERE IS NO DECIMAL PT

        DO n = 1,l
         read(unit=string(n:n), fmt='(i1)') i
         CHARTOFLOAT = CHARTOFLOAT + sign*float(i)*(10**(l-n))
        ENDDO

      ELSE                          ! CASE WHERE THERE IS A DECIMAL PT

        DO n = 1,m-1
         read(unit=string(n:n), fmt='(i1)') i
         CHARTOFLOAT = CHARTOFLOAT + sign*float(i)*(10.**(m-1-n))
        ENDDO

        DO n = m+1,l
         read(unit=string(n:n), fmt='(i1)') i
         CHARTOFLOAT = CHARTOFLOAT + sign*float(i)*(10.**(m-n))
        ENDDO

      ENDIF

      RETURN
      END FUNCTION CHARTOFLOAT

!-----------------------------------------------------------------------------
      INTEGER FUNCTION CHARTOINT(string)
      
      IMPLICIT NONE

      integer l, i, n , sign
      character( LEN = * ) string

      l = len(trim(string))

      sign = 1
      IF( string(1:1) .EQ. "-" ) THEN
       sign = -1
       string = string(2:l)
       l = l - 1
      ENDIF

      CHARTOINT = 0

! Parse each character (so we dont have to know the exact format)

      DO n = 1,l
       read(unit=string(n:n), fmt='(i1)') i
       CHARTOINT = CHARTOINT + sign*i*(10**(l-n))
      ENDDO

      RETURN
      END FUNCTION CHARTOINT

 
