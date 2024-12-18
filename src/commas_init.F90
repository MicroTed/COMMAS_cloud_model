!-----------------------------------------------------------------------------
!    
! INIT for COMMAS-1.0
!    
! 02/27/06
!    
!-----------------------------------------------------------------------------
 PROGRAM INIT_COMMAS
  
  USE ENS_MODULE
  USE GRIDIO_MODULE
  USE FILE_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
  USE COMMASMPI_MODULE

  implicit none

!-----------------------------------------------------------------------------
! Other local variables

  integer ls2, n, ns, i
!  real dx, dy, dz
  integer time
  integer coards(6)

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

  integer :: ixb, jyb, kzb
  integer :: ixe, jye, kze

!-----------------------------------------------------------------------------
! Derived type definitions

!  TYPE(ENSEMBLE)           :: ens

!-----------------------------------------------------------------------------
! Name of input file

  character(LEN=6) :: input_file = '3d.run'

!-----------------------------------------------------------------------------
! MPI STARTUP
#ifdef MPI
   CALL COMMASMPI_STARTUP()
#endif

!-----------------------------------------------------------------------------
! READ IN RUN ATTRIBUTES

  print*, 'COMMAS_INIT: READ RUN'
  IF( .not. READ_NAMELIST(input_file,'RUN') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! READ IN ATTRIBUTES

  print*, 'COMMAS_INIT: READ ATTRIBUTE'
  IF( .not. READ_NAMELIST(input_file,'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! READ IN GRID PARAMETERS

  print*, 'COMMAS_INIT: READ GRIDN'
  IF( .not. READ_NAMELIST(input_file,'GRIDN') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! READ IN MPI PARAMETERS
#ifdef MPI
  IF( .not. READ_NAMELIST(input_file,'MPI_PARAMS') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
#endif
! READ IN HOMOGENEOUS PARAMETERS

  print*, 'COMMAS_INIT: READ HOMOG_INIT'
  IF( .not. READ_NAMELIST(input_file,'HOMOG_INIT') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! READ IN ENSEMBLE PARAMETERS

  print*, 'COMMAS_INIT: READ ENS_INIT'

  IF( .not. READ_NAMELIST(input_file,'ENS_INIT') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! Read in lfo parameters
    IF( .not. READ_NAMELIST(input_file,'LFO_PARAMS') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

! Read in 10-ice parameters
    IF( .not. READ_NAMELIST(input_file,'ICE10_PARAMS') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

  print*, 'COMMAS_INIT: DONE WITH NAMELISTS'

!-----------------------------------------------------------------------------
! Compute some needed variables

  dx     = xdomain / float(nx-1)
  dy     = ydomain / float(ny-1)
  dz     = zdomain / float(nz-1)
  time   = 0

!  CALL STRING_LIMITS(prefix, i, ls)

  CALL STRING_LIMITS(ens_init_file, i, ls2)

!-----------------------------------------------------------------------------
! ALLOCATE ENSEMBLE, DEFINE GRIDS, AND ALLOCATE GRID MEMORY

  allocate(ens%g(member:member+ne-1))

  ens%prefix  = prefix
  coards(1)   = year
  coards(2)   = month
  coards(3)   = day
  coards(4)   = hour
  coards(5)   = minute
  coards(6)   = second

  DO n = member,member+ne-1

   CALL GRID_DEFINE_FROM_LIST(ens%g(n), microphys, ipconc)

   CALL GRID_SET_ATTRIBUTES(ens%g(n),ens%prefix,ne,n,coards  &
                           ,nx,ny,nz,microphys,v5dflds   &
#ifdef MPI   
                           ,nxend,nyend,nzend,my_rank)
#else
                           )
#endif

   CALL GRID_ALLOCATE(ens%g(n), nx, ny, nz)

  ns = size(ens%g(n)%var(:)) - GET_VARIABLE_INDEX(ens%g(n),'TH') + 1
  
  
  IF( .not. SET_VARIABLE(ens%g(n),'TIME',     time)    ) write(6,*) 'COMMAS:  Problem setting TIME'
  IF( .not. SET_VARIABLE(ens%g(n),'XBAR',       0.0)   ) write(6,*) 'COMMAS:  Problem setting XBAR'  
  IF( .not. SET_VARIABLE(ens%g(n),'YBAR',       0.0)   ) write(6,*) 'COMMAS:  Problem setting YBAR'  
  IF( .not. SET_VARIABLE(ens%g(n),'UGRID',    ugrid)   ) write(6,*) 'COMMAS:  Problem setting UGRID'  
  IF( .not. SET_VARIABLE(ens%g(n),'VGRID',    vgrid)   ) write(6,*) 'COMMAS:  Problem setting VGRID'
  IF( .not. SET_VARIABLE(ens%g(n),'NSCALAR',  ns)      ) write(6,*) 'COMMAS:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(ens%g(n),'NG',       ng)      ) write(6,*) 'COMMAS:  Problem setting NG'

   luno = FILE_OPEN(ens%g(n), 'OUTPUT_FILE_NAME', 'START', .true.)

   CALL GRID_DEFINE_NETCDF(ens%g(n))

!--------------------------------------------------------------------------------------
! Open I/O file for model output
                
   CALL FILE_MESSAGE('INIT COMMAS_INIT:  ENSEMBLE IS DEFINED, NETCDF FILE DEFINED, AND MEMORY IS ALLOCATED')

   luno = FILE_CLOSE()

  ENDDO

!-----------------------------------------------------------------------------
! INIT GRIDS

  DO n = member,member+ne-1

   
   CALL INIT_GRID(ens%g(n)) !,                              &
!                  dt, nsmall, dx, dy, dz,                &
!                  n_middle, nbndlyr,                     &
!                  dx_stretch,  dy_stretch, dz_stretch, dzmax,  &
!                  x_sw_loc, y_sw_loc, lat, lon, hgt, dzmaxtop, rtop, ztopstr)

   
   CALL PRTINFO(ens%g(n))

 
   CALL INIT_BACKGROUND_1D(ens%g(n) ) !,                               &
!                           sndtype, sndfile,                       &
!                           wtype, Us, Uz, psfc, tsfc, qsfc, rhmax, &
!                           shape, dudz0, dudz1, dudz2, z0, z1, z2)
                           
   CALL PRTBASE(ens%g(n))
 
   IF( ens_init_type .eq. -1 ) THEN
 
   CALL INIT_PERT_BBLE(ens%g(n)) !,                            &
!                       tbble, xrad,  yrad,  zrad,           &
!                       xcntr, ycntr, zcntr, bbletype, nbble,&
!                       bblsp, ibbleseed)
   ENDIF

   CALL PRINT( ens%g(n) )

  ENDDO 

!--------------------------------------------------------------------------------------
! ENKF ENSEMBLE INIT...

! IF( ens_init_type .eq. 1 .and. ne .lt. 2 ) THEN
!   
!  write(6,*) 'COMMAS_INIT:  ENS INIT TYPE = 1 AND NE < 2:  PROBLEM!!!!'
!  write(6,*) 'COMMAS_INIT:  ENS INIT TYPE = 1 AND NE < 2:  PROBLEM!!!!'
!  write(6,*) 'COMMAS_INIT:  SET ENS INIT TYPE = -1 OR NE > 2 !!!!'
!  write(6,*) 'COMMAS_INIT:  STOPPING !!'
!  stop
! ENDIF
  
! IF( ens_init_type .eq. 1 .and. ne .ge. 2 ) THEN
!
!  open(unit = 11, file=prefix(1:ls)//'.out', status='replace')
!
!  CALL ENSINIT(ens%g, ne, seed,                                      &
!               mean_u_profile_pert, mean_v_profile_pert,             &
!               mean_t_profile_pert, mean_td_profile_pert,            &
!               member_u_profile_pert, member_v_profile_pert,         &
!               member_t_profile_pert, member_td_profile_pert,        &
!               perttype, rbubh, rbubv, nb,                           &
!               upert, vpert, wpert, tpert, qpert,                    &
!               xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!
!  write(6,*) '  INIT:  ENS INIT TYPE = 1 / GRIDS ARE INITIALIZED FOR ENKF RUN'
!
!  close(11)
!
!  DO n = member,member+ne-1
!
!   luno = FILE_OPEN(ens%g(n), 'OUTPUT_FILE_NAME', 'AFTER ENSINIT')
!
!   CALL PRTBASE(ens%g(n))
!
!   CALL PRINT( ens%g(n) )
!
!  ENDDO
!
! ENDIF

!-----------------------------------------------------------------------------
! DUMP MAX/MIN, WRITE OUT NETCDF FILE..

  DO n = member,member+ne-1

   luno = FILE_OPEN(ens%g(n),'OUTPUT_FILE_NAME')

   CALL GRID_WRITE_NETCDF( ens%g(n), time )

   write(6,*) '  INIT:  GRID # ', n, ' HAS BEEN WRITTEN TO A NETCDF FILE...'
 
!   luno = FILE_CLOSE()

  ENDDO 

  deallocate(ens%g)

 END PROGRAM INIT_COMMAS
