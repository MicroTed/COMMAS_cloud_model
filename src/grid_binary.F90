!-------------------------------------------------------------------------------
!     
!     
!
!    
!
!   
!
!-------------------------------------------------------------------------------
 SUBROUTINE GRID_WRITE_BINARY(gd, file, iounit)

  
   character(LEN=*), optional :: file
   integer, optional          :: iounit
   TYPE(GRID)                 :: gd
   TYPE(ATTRIBUTE), pointer   :: attr
   TYPE(VARIABLE),  pointer   :: var2
   TYPE(VARIABLE)             :: var
    
   integer l, n, m, status
   integer ncid, nt, record_id, records
   integer variable_id
   integer time
   integer n1, n2, n3, n4, stag, istag, jstag, kstag
!   integer, allocatable :: dims(:)
   character(LEN=100) filename
   character :: member_string*4, time_string*15
   logical exist
   integer ibeg, iend, ll
   integer year, month, day, hour, minute, second
   integer iunit
   
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and file type

   iunit = 44
   IF ( present( iounit ) ) iunit = iounit
   
   IF( PRESENT(file) ) THEN

    filename = file
    write(0,*) 'GRID_WRITE_BINARY: infile,unit ', file,filename,iunit
 
!    CALL STRING_LIMITS( file, ibeg, iend)
!    ll = iend - ibeg + 1
!    filename(1:ll) = file(ibeg:iend)
       
   ELSE
       
    filename = ' '
    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    ll = iend - ibeg + 1
    filename(1:ll) = attr%str(ibeg:iend)

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_BINARY:  ', filename(1:ll)

    CALL GET_ATTRIBUTE(gd, 'MEMBER', attr)
    write(member_string,"(a1,i3.3)") '.', attr%int
    filename(1:ll+4) = filename(1:ll) // member_string
    ll = ll + 4

    IF( DEBUG_IO ) write(0,*) 'GRID_WRITE_BINARY:  ', filename(1:ll), attr%int

    CALL GET_VARIABLE(gd, 'YEAR',   year)
    CALL GET_VARIABLE(gd, 'MONTH',  month)
    CALL GET_VARIABLE(gd, 'DAY',    day)
    CALL GET_VARIABLE(gd, 'HOUR',   hour)
    CALL GET_VARIABLE(gd, 'MINUTE', minute)
    CALL GET_VARIABLE(gd, 'SECOND', second)
    CALL GET_VARIABLE(gd, 'TIME', time)

    write(time_string,"(a1,i4.4,5(i2.2))") '.', year, month, day, hour, minute, second
    filename(1:ll+19) = filename(1:ll) // time_string // '.bin'
    ll = ll + 19
    IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_BINARY:  ', filename(1:ll)

   ENDIF

   IF( .NOT. overwrite_file ) THEN

    INQUIRE(file = filename(1:ll), EXIST = exist)

!    IF( exist ) THEN

!     CALL SYSTEM('mv ' // filename(1:ll) // ' ' // filename(1:ll) // '.0')

!    ENDIF

   ENDIF
   
   
   IF ( .not. present( iounit ) ) THEN
   write(0,*) 'GRID_WRITE_BINARY: file,unit ', filename,iunit
   open(unit = iunit, file = trim(filename), status='unknown', form = 'unformatted')
   ENDIF
   rewind(iunit)

!-----------------------------------------------------------------------------------------------
! Write number of attributes and variables

  write(iunit) size(gd%attr)
  write(iunit) size(gd%var)

!-----------------------------------------------------------------------------------------------
! Write Attributes

  DO n = 1,size(gd%attr)
                                  write(iunit) gd%attr(n)%name, gd%attr(n)%type
   IF( gd%attr(n)%type == 'int' ) write(iunit) gd%attr(n)%int
   IF( gd%attr(n)%type == 'flt' ) write(iunit) gd%attr(n)%flt
   IF( gd%attr(n)%type == 'str' ) write(iunit) gd%attr(n)%str

  ENDDO

!-----------------------------------------------------------------------------------------------
! Write Variables

  DO n = 1,size(gd%var)

   write(iunit) gd%var(n)%name,  gd%var(n)%type,    gd%var(n)%dim,      &
             gd%var(n)%istag, gd%var(n)%jstag,   gd%var(n)%kstag,    &
             gd%var(n)%pdef,  gd%var(n)%dyntype, gd%var(n)%phytype,  &
             gd%var(n)%buotype, gd%var(n)%unit, gd%var(n)%description

   IF( DEBUG_IO ) write(6,*) 'GRID_WRITE_BINARY:  ', gd%var(n)%name, gd%var(n)%type

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, write out the scalar

    IF( gd%var(n)%dim == 0 ) THEN

     IF( gd%var(n)%type(1:5) == 'icnst' ) THEN
       write(iunit) gd%var(n)%int
     ELSE
       write(iunit) gd%var(n)%flt
     ENDIF
!     IF( gd%var(n)%type(1:5) == 'rcnst' ) write(iunit) gd%var(n)%flt

!-------------------------------------------------------------------------------
! Else we are writing out a 1, 2, 3D array 

    ELSE

!     allocate(dims(gd%var(n)%dim))

!     dims(:) = 1
   
    IF( gd%var(n)%type(1:3) == 'x1d' .or. gd%var(n)%type(1:3) == 'y1d' .or. gd%var(n)%type(1:3) == 'z1d' ) THEN 

     n1 = size(gd%var(n)%flt1d)

     write(iunit) n1
     write(iunit) gd%var(n)%flt1d ! (1:n1)

    ENDIF
   
    IF( gd%var(n)%type(1:4) == 'xy2d' .or. gd%var(n)%type(1:4) == 'xz2d' .or. gd%var(n)%type(1:4) == 'yz2d' ) THEN 

     n1 = size(gd%var(n)%flt2d,dim=1)
     n2 = size(gd%var(n)%flt2d,dim=2)
    
     write(iunit) n1, n2
     write(iunit) gd%var(n)%flt2d ! (1:n1,1:n2)

    ENDIF
   
    IF( gd%var(n)%type(1:5) == 'xyz3d' ) THEN

     n1 = size(gd%var(n)%flt3d,dim=1) 
     n2 = size(gd%var(n)%flt3d,dim=2) 
     n3 = size(gd%var(n)%flt3d,dim=3) 

     write(iunit) n1, n2, n3
     write(iunit) gd%var(n)%flt3d ! (1:n1,1:n2,1:n3)

    ENDIF

! Write out information

!   deallocate(dims)

  ENDIF

  ENDDO

  IF ( .not. present( iounit ) ) THEN
    close(iunit)
  ENDIF

 END SUBROUTINE GRID_WRITE_BINARY


!-------------------------------------------------------------------------------
!     
!     
!
!    
!
!   
!
!-------------------------------------------------------------------------------
 SUBROUTINE GRID_READ_BINARY(gd, file, member)

   implicit none
   
   character(LEN=*), optional :: file
   integer, optional          :: member
   TYPE(GRID)                 :: gd
   TYPE(ATTRIBUTE), pointer   :: attr
   TYPE(VARIABLE),  pointer   :: var2
   TYPE(VARIABLE)             :: var
    
   integer l, n, m, status
   integer ncid, nt, record_id, records
   integer variable_id
   integer time
   integer n1, n2, n3, n4, stag, istag, jstag, kstag
!   integer, allocatable :: dims(:)
   character(LEN=120) filename
   character :: member_string*4, time_string*15
   logical exist
   integer ibeg, iend, ll
   integer year, month, day, hour, minute, second
   integer numattr, numvar
   integer iunit
   
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and file type

   iunit = 900
   IF ( present( member ) ) iunit = iunit + member
   
   IF( PRESENT(file) ) THEN
 
    CALL STRING_LIMITS( file, ibeg, iend)
    ll = iend - ibeg + 1
    filename(1:ll) = file(ibeg:iend)
       
   ELSE
       
    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    ll = iend - ibeg + 1
    filename(1:ll) = attr%str(ibeg:iend)

    IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', filename(1:ll)

    CALL GET_ATTRIBUTE(gd, 'MEMBER', attr)
    write(member_string,"(a1,i3.3)") '.', attr%int
    filename(1:ll+4) = filename(1:ll) // member_string
    ll = ll + 4

    IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', filename(1:ll), attr%int

    CALL GET_VARIABLE(gd, 'YEAR',   year)
    CALL GET_VARIABLE(gd, 'MONTH',  month)
    CALL GET_VARIABLE(gd, 'DAY',    day)
    CALL GET_VARIABLE(gd, 'HOUR',   hour)
    CALL GET_VARIABLE(gd, 'MINUTE', minute)
    CALL GET_VARIABLE(gd, 'SECOND', second)
    CALL GET_VARIABLE(gd, 'TIME', time)

    write(time_string,"(a1,i4.4,5(i2.2))") '.', year, month, day, hour, minute, second
    filename(1:ll+19) = filename(1:ll) // time_string // '.bin'
    ll = ll + 19
    IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', filename(1:ll)

   ENDIF

   IF( .NOT. overwrite_file ) THEN

    INQUIRE(file = filename(1:ll), EXIST = exist)

!    IF( exist ) THEN

!     CALL SYSTEM('mv ' // filename(1:ll) // ' ' // filename(1:ll) // '.0')

!    ENDIF

   ENDIF
   
   write(6,*) 'GRID_READ_BINARY:  file,unit ', filename(1:ll),iunit
   open(unit = iunit, file = filename(1:ll), status='old', form = 'unformatted')
   rewind(iunit)

!-----------------------------------------------------------------------------------------------
! Write number of attributes and variables

  READ(iunit) numattr ! size(gd%attr)
  READ(iunit) numvar  ! size(gd%var)
  
  IF ( numattr .ne. size(gd%attr) .or. numvar .ne. size(gd%var) ) THEN
    write(0,*) 'GRID_READ_BINARY:  numattr or numvar is incorrect: ',numattr, size(gd%attr),numvar,size(gd%var)
    STOP
  ENDIF
  
!-----------------------------------------------------------------------------------------------
! Read Attributes

  DO n = 1,size(gd%attr)
                                  read(iunit) gd%attr(n)%name, gd%attr(n)%type
   IF( gd%attr(n)%type == 'int' ) read(iunit) gd%attr(n)%int
   IF( gd%attr(n)%type == 'flt' ) read(iunit) gd%attr(n)%flt
   IF( gd%attr(n)%type == 'str' ) read(iunit) gd%attr(n)%str

  ENDDO

!-----------------------------------------------------------------------------------------------
! Read Variables

  DO n = 1,size(gd%var)

   read(iunit) gd%var(n)%name,  gd%var(n)%type,    gd%var(n)%dim,      &
             gd%var(n)%istag, gd%var(n)%jstag,   gd%var(n)%kstag,    &
             gd%var(n)%pdef,  gd%var(n)%dyntype, gd%var(n)%phytype,  &
             gd%var(n)%buotype, gd%var(n)%unit, gd%var(n)%description

   IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', gd%var(n)%name, gd%var(n)%type

!-------------------------------------------------------------------------------
! Determine dims -> if dim == 0, read in the scalar

    IF( gd%var(n)%dim == 0 ) THEN

     IF( gd%var(n)%type(1:5) == 'icnst' )  THEN
       read(iunit) gd%var(n)%int
     ELSE
       read(iunit) gd%var(n)%flt
     ENDIF

!     IF( gd%var(n)%type(1:5) == 'rcnst' ) read(iunit) gd%var(n)%flt

!-------------------------------------------------------------------------------
! Else we are reading in a 1, 2, 3D array 

    ELSE

!     allocate(dims(gd%var(n)%dim))

!     dims(:) = 1
   
    IF( gd%var(n)%type(1:3) == 'x1d' .or. &
        gd%var(n)%type(1:3) == 'y1d' .or. & 
        gd%var(n)%type(1:3) == 'z1d' ) THEN 

     n1 = size(gd%var(n)%flt1d)

     read(iunit) n1
     IF ( n1 .ne. size(gd%var(n)%flt1d) - 2*gd%var(n)%ng ) THEN
       write(0,*) 'GRID_READ_BINARY: problem with 1d array: ',n1,size(gd%var(n)%flt1d)
       write(0,*) 'GRID_READ_BINARY: type = ',gd%var(n)%type(1:3),', ',gd%var(n)%name
!       STOP
     ENDIF
     read(iunit) gd%var(n)%flt1d(1:n1)

    ENDIF
   
    IF( gd%var(n)%type(1:4) == 'xy2d' .or. &
        gd%var(n)%type(1:4) == 'xz2d' .or. &
        gd%var(n)%type(1:4) == 'yz2d' ) THEN 

    
     read(iunit) n1, n2
     IF ( n1 /= size(gd%var(n)%flt2d,dim=1) - 2*gd%var(n)%ng  .or. n2 /= size(gd%var(n)%flt2d,dim=2) - 2*gd%var(n)%ng  ) THEN
       write(0,*) 'GRID_READ_BINARY: problem with 2d array: ',n1,size(gd%var(n)%flt2d,dim=1), &
           n2, size(gd%var(n)%flt2d,dim=2)
!       STOP
     ENDIF
     read(iunit) gd%var(n)%flt2d(1:n1,1:n2)

    ENDIF
   
    IF( gd%var(n)%type(1:5) == 'xyz3d' ) THEN

     read(iunit) n1, n2, n3
     IF ( n1 /= size(gd%var(n)%flt3d,dim=1) - 2*gd%var(n)%ng  .or. n2 /= size(gd%var(n)%flt3d,dim=2) - 2*gd%var(n)%ng  .or. &
          n3 /= size(gd%var(n)%flt3d,dim=3) - 2*gd%var(n)%ng   ) THEN
       write(0,*) 'GRID_READ_BINARY: problem with 2d array: ',n1,size(gd%var(n)%flt3d,dim=1), &
           n2, size(gd%var(n)%flt3d,dim=2), n3, size(gd%var(n)%flt3d,dim=3)
!       STOP
     ENDIF
     read(iunit) gd%var(n)%flt3d(1:n1,1:n2,1:n3)

    ENDIF

! Write out information

!   deallocate(dims)

  ENDIF

  ENDDO

  close(iunit)

 END SUBROUTINE GRID_READ_BINARY

!-------------------------------------------------------------------------------
!     
!     
!
!    
!
!   
!
!-------------------------------------------------------------------------------
! SUBROUTINE GRID_READ_BINARY(gd, file, member, year, month, day, hour, minute, second)
!
!   USE STRING_MODULE
!
!   character(LEN=*)           :: file
!   integer, optional          :: member, year, month, day, hour, minute, second
!   TYPE(GRID)                 :: gd
!   TYPE(ATTRIBUTE), pointer   :: attr
!   TYPE(VARIABLE),  pointer   :: var
!    
!   integer l, n, m, status
!   integer ncid, nt, record_id, records
!   integer variable_id
!   integer n1, n2, n3, n4, nattr, nvar
!   integer, allocatable :: dims(:)
!   character(LEN=120) filename
!   character :: member_string*4, time_string*15, dummy*5
!   logical exist
!   integer ibeg, iend, ll
!   
!!-----------------------------------------------------------------------
!! Create filename based on prefix name, time, and file type
!
!   CALL STRING_LIMITS( file, ibeg, iend)
!   ll = iend - ibeg + 1
!   filename(1:ll) = file(ibeg:iend)
!
!   IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', filename(1:ll)
!
!   IF( PRESENT(member) .and. &
!       PRESENT(year)   .and. &
!       PRESENT(month)  .and. &
!       PRESENT(day)    .and. &
!       PRESENT(hour)   .and. &
!       PRESENT(minute) .and. & 
!       PRESENT(second)       ) THEN
!       
!    write(member_string,"(a1,i3.3)") '.', member
!    filename(1:ll+4) = filename(1:ll) // member_string
!    ll = ll + 4
!    
!    IF( DEBUG_IO ) write(6,*) 'GRID_READ_ASCII:  ', filename(1:ll)
!    
!    write(time_string,"(a1,i4.4,5(i2.2))") '.', year, month, day, hour, minute, second
!    filename(1:ll+19) = filename(1:ll) // time_string // '.bin'
!    ll = ll + 19 
!
!    IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  ', filename(1:ll)
!
!   ENDIF
!
!   open(unit = lun, file = filename(1:ll), status='unknown', form = 'unformatted')
!   rewind(lun)
!
!!-----------------------------------------------------------------------------------------------
!! Read and allocate for the number of attributes and variables
!
!  read(lun) nattr
!  allocate(gd%attr(nattr))
!
!  IF( DEBUG_IO ) write(6,*) 'READ_GRID_BINARY:  NATTR = ', nattr
!
!  read(lun) nvar
!  allocate(gd%var(nvar))  
!
!  IF( DEBUG_IO ) write(6,*) 'READ_GRID_BINARY:  NVAR  = ', nvar 
!
!!-----------------------------------------------------------------------------------------------
!! Write Attributes
!
!  DO n = 1,nattr
!                                  read(lun) gd.attr(n).name, gd.attr(n)%type
!   IF( gd.attr(n)%type == 'int' ) read(lun) gd.attr(n).int
!   IF( gd.attr(n)%type == 'flt' ) read(lun) gd.attr(n).flt
!   IF( gd.attr(n)%type == 'str' ) read(lun) gd.attr(n).str
!
!   IF( DEBUG_IO ) write(6,*) 'READ_GRID_BINARY:  ATTR NAME = ', gd.attr(n).name 
!
!  ENDDO
!
!!-----------------------------------------------------------------------------------------------
!! Write Attributes
!
!102 format(1x,a,2x,a,2x,i1,2x,5(i3,2x),a,2x,a)
!
!  DO n = 1,nvar
!
!   read(lun) gd.var(n).name, gd.var(n)%type, gd.var(n).tdepend, gd.var(n).dim, gd.var(n).stag, &
!             gd.var(n).pdef, gd.var(n).advtype, gd.var(n).mixtype, gd.var(n).unit, gd.var(n).description
!
!   IF( DEBUG_IO ) write(6,*) 'READ_GRID_BINARY:  VAR NAME = ', gd.var(n).name , gd.var(n)%type
!
!!-------------------------------------------------------------------------------
!! Determine dims -> if dim == 0, read in the scalar
!
!   IF( gd.var(n)%type(1:3) == 'int' ) read(lun) gd.var(n).int
!   IF( gd.var(n)%type(1:3) == 'flt' ) read(lun) gd.var(n).flt
!
!!-------------------------------------------------------------------------------
!! Else we are writing out a 1, 2, 3, or 4D array 
!
!   IF( gd.var(n)%type(1:3) == 'x1d' .or. gd.var(n)%type(1:3) == 'y1d' .or. gd.var(n)%type(1:3) == 'z1d' ) THEN 
!
!     read(lun) n1 
!
!     allocate(gd.var(n).flt1d(n1))
!
!     read(lun) gd.var(n).flt1d(1:n1)
!
!   ENDIF
!  
!   IF( gd.var(n)%type(1:4) == 'xy2d' .or. gd.var(n)%type(1:4) == 'xz2d' .or. gd.var(n)%type(1:4) == 'yz2d' ) THEN 
!
!    read(lun) n1, n2
!
!    allocate(gd.var(n).flt2d(n1,n2))
!    
!    read(lun) gd.var(n).flt2d(1:n1,1:n2)
!
!   ENDIF
!   
!   IF( gd.var(n)%type(1:5) == 'xyz3d' ) THEN
!
!    read(lun) n1, n2, n3
!
!    allocate(gd.var(n).flt3d(1:n1,1:n2,1:n3))
!
!    read(lun) gd.var(n).flt3d(1:n1,1:n2,1:n3)
!
!   ENDIF
!
!   IF( gd.var(n)%type(1:5) == 'xyz4d' ) THEN
!
!    read(lun) n1, n2, n3, n4
!
!    IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY:  XYZ4D  = ', n1, n2, n3, n4
!
!    allocate(gd.var(n).flt4d(n1,n2,n3,n4))
!    allocate(gd%var(n)%names(n4))
!    allocate(gd%var(n)%dims(n4))
!    allocate(gd%var(n)%types(n4))
!    allocate(gd%var(n)%tdepends(n4))
!    allocate(gd%var(n)%stags(n4))
!    allocate(gd%var(n)%pdefs(n4))
!    allocate(gd%var(n)%advtypes(n4))
!    allocate(gd%var(n)%mixtypes(n4))
!    allocate(gd%var(n)%units(n4))
!    allocate(gd%var(n)%descriptions(n4))
!
!    DO m = 1,n4
!
!     IF( DEBUG_IO ) write(6,*) 'GRID_READ_BINARY: M = ',m
!
!     read(lun) gd.var(n).names(m), dummy, gd.var(n).tdepends(m), gd.var(n).dims(m), gd.var(n).stags(m), &
!               gd.var(n).pdefs(m), gd.var(n).advtypes(m), gd.var(n).mixtypes(m), gd.var(n).units(m), gd.var(n).descriptions(m)
!
!     IF( DEBUG_IO ) write(6,102) gd.var(n).names(m), dummy, gd.var(n).tdepends(m), gd.var(n).dims(m), gd.var(n).stags(m), &
!                                 gd.var(n).pdefs(m), gd.var(n).advtypes(m), gd.var(n).mixtypes(m), gd.var(n).units(m),    &
!                                 gd.var(n).descriptions(m)
!
!     read(lun) gd.var(n).flt4d(1:n1,1:n2,1:n3,m)
!
!    ENDDO
!
!   ENDIF
!
!  ENDDO
!
!  close(lun)
!
! END SUBROUTINE GRID_READ_BINARY

!-------------------------------------------------------------------------------
 SUBROUTINE GRID_INFO_BINARY(file, nx, ny, nz, member)

   implicit none
   
   character(LEN=*)           :: file
!   TYPE(GRID)                 :: gd
   integer nx,ny,nz
   integer, optional          :: member
   TYPE(ATTRIBUTE), pointer   :: attr
   TYPE(VARIABLE),  pointer   :: var2
   TYPE(VARIABLE)             :: var
    
   integer l, n, m, status
   integer ncid, nt, record_id, records
   integer variable_id
   integer time
   integer n1, n2, n3, n4, stag, istag, jstag, kstag
   integer, allocatable :: dims(:)
   character(LEN=120) filename
   character :: member_string*4, time_string*15
   logical exist
   integer ibeg, iend, ll
   integer year, month, day, hour, minute, second
   integer numattr, numvar
   integer iunit
   
   character(LEN=name_length)   attname
   character(LEN=type_length)   atttype
   integer                      attint
   real                         attflt
   character(LEN=desc_length) :: attstr
   
!-----------------------------------------------------------------------
! Create filename based on prefix name, time, and file type

 
    iunit = 900
    IF ( present( member ) ) iunit = iunit + member
    
    CALL STRING_LIMITS( file, ibeg, iend)
    ll = iend - ibeg + 1
    filename(1:ll) = file(ibeg:iend)
       

   IF( .NOT. overwrite_file ) THEN

    INQUIRE(file = filename(1:ll), EXIST = exist)

!    IF( exist ) THEN

!     CALL SYSTEM('mv ' // filename(1:ll) // ' ' // filename(1:ll) // '.0')

!    ENDIF

   ENDIF
   
   open(unit = iunit, file = filename(1:ll), status='unknown', form = 'unformatted')
   rewind(iunit)

!-----------------------------------------------------------------------------------------------
! Write number of attributes and variables

  READ(iunit) numattr ! size(gd%attr)
   write(0,*) 'GRID_INFO_BINARY: numattr = ', numattr
  READ(iunit) numvar  ! size(gd%var)
   write(0,*) 'GRID_INFO_BINARY: numvar = ', numvar
  
!  IF ( numattr .ne. size(gd%attr) .or. numvar .ne. size(gd%var) ) THEN
!    write(0,*) 'GRID_READ_BINARY:  numattr or numvar is incorrect: ',numattr, size(gd%attr),numvar,size(gd%var)
!    STOP
!  ENDIF
  
!-----------------------------------------------------------------------------------------------
! Read Attributes

  DO n = 1,numattr
   attint = 0
   attflt = 0
   attstr = ' '
   read(iunit) attname, atttype
   write(0,*) 'GRID_INFO_BINARY: attname,atttype = ',attname, atttype
   
   IF( atttype == 'int' ) read(iunit) attint
   IF( atttype == 'flt' ) read(iunit) attflt
   IF( atttype == 'str' ) read(iunit) attstr
   
   IF ( attname == 'NXEND' ) nx = attint
   IF ( attname == 'NYEND' ) ny = attint
   IF ( attname == 'NZEND' ) nz = attint
   
   write(0,*) 'GRID_INFO_BINARY: attname,atttype,attint,attflt,attstr = ',attname,atttype,attint,attflt,attstr

  ENDDO


  RETURN
  
  END SUBROUTINE GRID_INFO_BINARY
