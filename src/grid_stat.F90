!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE GRID_INFO(gd,prt_all)
 
   implicit none
   
   type(GRID) :: gd
   integer, optional :: prt_all

   integer n

   print *, '=========================================================='
   print *, ''
   print *, 'GRID INFORMATION'
   print *, ''
   print *, '# of attributes        ', size(gd%attr)
   print *, '# of vars              ', size(gd%var)
   print *, ''
   print *, '# of real constants    ', count(gd%var(:)%type ==  'icnst')
   print *, '# of integer constants ', count(gd%var(:)%type ==  'rcnst')
   print *, '# of X1D arrays        ', count(gd%var(:)%type ==    'x1d')   
   print *, '# of Y1D arrays        ', count(gd%var(:)%type ==    'y1d')   
   print *, '# of Z1D arrays        ', count(gd%var(:)%type ==    'z1d')   
   print *, '# of XY2D arrays       ', count(gd%var(:)%type ==   'xy2d')  
   print *, '# of XZ2D arrays       ', count(gd%var(:)%type ==   'xz2d')  
   print *, '# of YZ2D arrays       ', count(gd%var(:)%type ==   'yz2d')  
   print *, '# of 3D arrays         ', count(gd%var(:)%type ==  'xyz3d')
   print *, ''
   print *, '=========================================================='

   IF( PRESENT(prt_all) ) THEN

    print *, ''
    print *, 'GRID ATTRIBUTES'
    print *, ''
    DO n = 1,size(gd%attr)
      IF( gd%attr(n)%type == 'str') &
        print *, 'ATTRIBUTE NAME:  ',gd%attr(n)%name,' VALUE:  ',gd%attr(n)%str
      IF( gd%attr(n)%type == 'int') &
        print *, 'ATTRIBUTE NAME:  ',gd%attr(n)%name,' VALUE:  ',gd%attr(n)%int
      IF( gd%attr(n)%type == 'flt') &
        print *, 'ATTRIBUTE NAME:  ',gd%attr(n)%name,' VALUE:  ',gd%attr(n)%flt
    ENDDO
    print *, ''
    print *, '=========================================================='

    IF( prt_all == 2 ) THEN

     print *, ''
     print *, 'GRID VARIABLES'
     print *, ''
     DO n = 1,size(gd%var)
      print *, 'VARIABLE NAME:  ',gd%var(n)%name,' VARIABLE DIMS:  ', gd%var(n)%dim, gd%var(n)%index
     ENDDO
     print *, ''
     print *, '=========================================================='

    ENDIF

   ENDIF

  RETURN
  END SUBROUTINE GRID_INFO





!-------------------------------------------------------------------------------
! 
!
!  
!  
!
!  
!  
!-------------------------------------------------------------------------------

  SUBROUTINE VAR_MAXMIN(var,mxmn,bounds)

   USE COMMASMPI_MODULE

   implicit none
   
   type(VARIABLE)          :: var
   type(MAXMIN)            :: mxmn
   integer,dimension(6), optional       :: bounds

   integer i, j, k, nx, ny, nz
   integer i0, i1, j0, j1, k0, k1
   real value

   IF( var%type(1:5) == 'xyz3d' ) THEN

#ifdef MPI
      if (ixend==nxend) then
        nx = size(var%flt3d,dim=1) - 2*var%ng - 1
      else
        nx = size(var%flt3d,dim=1) - 2*var%ng
      endif
      
      if (jyend==nyend) then
        ny = size(var%flt3d,dim=2) - 2*var%ng - 1
      else
        ny = size(var%flt3d,dim=2) - 2*var%ng
      endif
      
      if (kzend==nzend) then
        nz = size(var%flt3d,dim=3) - 2*var%ng - 1
      else
        nz = size(var%flt3d,dim=3) - 2*var%ng
      endif
      
#else
    nx = size(var%flt3d,dim=1) - 2*var%ng - 1
    ny = size(var%flt3d,dim=2) - 2*var%ng - 1
    nz = size(var%flt3d,dim=3) - 2*var%ng - 1
#endif

#ifdef MPI
    IF( var%istag == 1 .and. ixend == nxend ) nx = nx + 1
    IF( var%jstag == 1 .and. jyend == nyend ) ny = ny + 1
    IF( var%kstag == 1 .and. kzend == nzend ) nz = nz + 1
#else
    IF( var%istag == 1 ) nx = nx + 1
    IF( var%jstag == 1 ) ny = ny + 1
    IF( var%kstag == 1 ) nz = nz + 1
#endif

    i0 = 1
    i1 = nx

    j0 = 1
    j1 = ny

    k0 = 1
    k1 = nz

    IF( PRESENT(bounds) ) THEN

     IF( bounds(1) .NE. -1 ) i0 = max(bounds(1),1)
     IF( bounds(2) .NE. -1 ) i1 = min(bounds(2),nx)
     IF( bounds(3) .NE. -1 ) j0 = max(bounds(3),1)
     IF( bounds(4) .NE. -1 ) j1 = min(bounds(4),ny)
     IF( bounds(5) .NE. -1 ) k0 = max(bounds(5),1)
     IF( bounds(6) .NE. -1 ) k1 = min(bounds(6),nz)

    ENDIF

    mxmn%name = var%name
    mxmn%max  = var%flt3d(1,1,1)
    mxmn%min  = var%flt3d(1,1,1)
    mxmn%imax = i0
    mxmn%jmax = j0
    mxmn%kmax = k0
    mxmn%imin = i0
    mxmn%jmin = j0
    mxmn%kmin = k0

    DO k = k0,k1
     DO j = j0,j1
      DO i = i0,i1

      value = var%flt3d(i,j,k)

      IF( value .gt. mxmn%max ) THEN
       mxmn%max   = value
       mxmn%imax  = i
       mxmn%jmax  = j
       mxmn%kmax  = k
      ENDIF

      IF( value .lt. mxmn%min ) then
       mxmn%min   = value
       mxmn%imin  = i
       mxmn%jmin  = j
       mxmn%kmin  = k
      ENDIF

      ENDDO
     ENDDO
    ENDDO

   ENDIF

  END SUBROUTINE VAR_MAXMIN
