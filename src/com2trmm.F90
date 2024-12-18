!-----------------------------------------------------------------------------
!
! COM2TRMM-1.0
!
! March 2006
!
!-----------------------------------------------------------------------------

 PROGRAM COM2TRMM

#ifdef NAG
  USE F90_UNIX_PROC
#endif

  USE GRID_MODULE
  USE GRIDIO_MODULE
  USE FILE_MODULE
  USE CLINE_MODULE
  USE CPUTIME_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
!  USE VIS5D_MODULE
  USE TRMM_MODULE
  USE INDEX_MODULE

  implicit none

!-----------------------------------------------------------------------------
! GRID DEFINITION
 
  TYPE(GRID) :: gd, ge

!-----------------------------------------------------------------------------
! COMMAND LINE VARIABLES

  character(LEN = 25) :: run_file
  integer             :: status, length, length1
  real                :: tmp, chartofloat

!-----------------------------------------------------------------------------
! DATE DIAGNOSTICS

  character( LEN = 80 ) :: command
  character( LEN = 8  ) :: date
  character( LEN = 10 ) :: hhmmss
  character( LEN = 3  ), dimension(12) :: months = (/'JAN','FEB',  &
                                                     'MAR','APR',  &
                                                     'MAY','JUN',  &
                                                     'JUL','AUG',  &
                                                     'SEP','OCT',  &
                                                     'NOV','DEC'/)
!-----------------------------------------------------------------------------
! INTEGER INDEX VARIABLES FOR THE RUN

  integer :: u 
  integer :: v  
  integer :: w   
  integer :: pi   
  integer :: km   
  integer :: s     
  integer :: dbz   
  integer :: wz    
  integer :: elec
  integer :: cion
  integer :: muz
   
  integer :: uinit  
  integer :: vinit  
  integer :: winit  
  integer :: piinit 
  integer :: kminit 
  integer :: sinit  

  integer :: gx     
  integer :: gy     
  integer :: gz     
  integer :: precip 

!-----------------------------------------------------------------------------
! SOLVER SCRATCH MEMORY

  real, allocatable    :: st(:,:,:,:)         ! These are big scratch arrays for solver
  real, allocatable    :: flsh(:,:,:,:,:)
  real, allocatable    :: ft(:,:)
  integer, allocatable :: tarray(:)

!-----------------------------------------------------------------------------
! OTHER MISC VARIABLES 

!  real    :: dx
!  real    :: dy
!  real    :: dz
  integer :: ibeg
  integer :: iend
  integer :: ns
  integer :: n
  integer :: ntime
  integer :: m
  integer :: time
  integer :: coards(6)
  integer :: trst       
  integer :: this       
  integer :: tprt
  integer :: tstt
  integer :: tv5d    
  integer :: tvis5dtmp

  logical :: file_exist = .false.
  logical :: io_flag, lstt
  logical :: hack = .false.
  
  integer :: i,j,k, it, nt, count
  
  integer :: iflash, icgn, icgp

  character(LEN=4) number
  character(LEN=2) dup

  character(LEN=20) :: filenm
  character(LEN = 255):: v5dfldstmp

!-----------------------------------------------------------------------------
! Need 3D.RUN file in order to get run parameters - get from the command line

  IF( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN

   write(0,*) 'COM2TRMM:  INCORRECT ARGUMENTS ON COMMAND LINE:  NEED 3D.RUN FILENAME'
   write(0,*) 'COM2TRMM:  INCORRECT ARGUMENTS ON COMMAND LINE:  EXITING RUN!...'
   call exit(1)

  ELSE

   CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)

!   write(6,*) 'COM2TRMM:  Input 3D.RUN file is:  ', run_file(1:length)

   IF( status .ne. 0 ) THEN

     write(0,*) 'COM2TRMM:  COULD NOT RETRIEVE 1ST COMMAND LINE ARG:  EXITING RUN!...', status
     write(0,*) 'COM2TRMM:  DOES NOT HAVE A 3D.RUN FILE TO READ!'
     call exit(1)

   ENDIF

  ENDIF
  
!-----------------------------------------------------------------------------
! Read in 3d.run file...

  INQUIRE(file=run_file(1:length), exist=file_exist)

  IF( .NOT. file_exist ) THEN

    write(0,*) 'COM2TRMM:  INPUT 3D.RUN FILE:  ', run_file(1:length), ' DOES NOT EXIST!!!'
    write(0,*) 'COM2TRMM:  DOES NOT HAVE A 3D.RUN FILE TO READ, EXITING'
    stop

  ENDIF

! READ IN RUN MODEL NAMELIST

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'COM2TRMM:  PROBLEM READING NAMELIST'

  write(number, '(a,i3.3)') '.', member

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COM2TRMM:  PROBLEM READING NAMELIST'

  IF( .not. READ_NAMELIST(run_file(1:length),'TRMM') ) print *, 'COM2TRMM:  PROBLEM READING NAMELIST'

  
! READ IN GRID PARAMETERS


  ibsd = 1
  iesd = nx-1
  jbsd = 1
  jesd = ny-1
  lpredict = .false.
  
  v5dfldstmp = v5dflds

   start = gridtimes(1)


!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   CALL STRING_LIMITS(prefix, ibeg, iend)
!   v5dfilename = prefix(ibeg:iend)//number//'.v5d'


!-----------------------------------------------------------------------------
! ELSE, READ IN THE DATA


    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

    write(6,*) 'COM2TRMM:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//number
  
    CALL GRID_INFO_NETCDF( prefix(1:length)//number//'.nc', nt, microphys)
    
    allocate( tarray(nt) )
    
    CALL GRID_INFO_NETCDF( prefix(1:length)//number//'.nc', nt, microphys, tarray)
    
    
    CALL GRID_READ_NETCDF( gd, prefix(1:length)//number//'.nc', tarray(1) )

    write(6,*) 'COM2TRMM:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//number

    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COM2TRMM:  Problem setting TIME'  
    IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(6,*) 'COM2TRMM:  Problem setting UGRID'  
    IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(6,*) 'COM2TRMM:  Problem setting VGRID'
    IF( .not. SET_VARIABLE(gd,'TRESTART',  trestart)) write(6,*) 'COM2TRMM:  Problem setting TRESTART'
    IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COM2TRMM:  Problem setting THISTORY'
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COM2TRMM:  Problem setting TVIS5D'
    IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COM2TRMM:  Problem setting TPRINT'
    IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COM2TRMM:  Problem setting TSTAT'
    IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COM2TRMM:  Problem setting STOP'
    IF( .not. SET_VARIABLE(gd,'DT',        dt)      ) write(6,*) 'COM2TRMM:  Problem setting DT'
   
    CALL STRING_LIMITS(v5dfldstmp, i, j)
    IF( .not. SET_ATTRIBUTE(gd,'V5DFIELDS',  v5dfldstmp)) write(6,*) 'COM2TRMM:  Problem setting V5DFIELDS'

    time = start


!--------------------------------------------------------------------------------------
! Open VIS5D file  

!  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COM2TRMM:  PROBLEM READING NAMELIST'
    
    IF ( start .gt. tarray(nt) ) THEN
      write(0,*) 'Starting time is greater than the last record time! STOP!'
      STOP
    ENDIF
    
    start = Max(start, tarray(2) ) 
    
!    DO it = nt,1,-1
!      IF ( start .le. tarray(it) ) v5dtimes(1) = tarray(it)
!    ENDDO
    
!    v5dtimes(1) = Max( start, tarray(1) )

!    tvis5dtmp = tvis5d ! v5dtimes(1)
!    print*, 'setting TIME to ',Max(0,start-tvis5d)
!    IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start-tvis5d))   ) write(6,*) 'COM2TRMM:  Problem setting TIME'  

!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5dtmp)  ) write(6,*) 'COM2TRMM:  Problem setting TVIS5D'  


   filenm = prefix
   
   IF ( trmmstatfile .eq. 'none' ) THEN
    trmmstatfile = filenm(1:length)//number//'.stat'
   ENDIF
   
   count = 0
5   INQUIRE(file=trmmstatfile, exist=file_exist)

   IF ( file_exist ) THEN
     write(6,*) 'stat file already exists!  I will not overwrite -- changing name!'
     count = count + 1
     write(dup, '(i2.2)')  count
     
     IF ( count .ge. 100 ) STOP

     trmmstatfile = prefix(ibeg:iend)//number//'.'//dup//'.stat'

     GOTO 5
   ENDIF
   
   write(6,*) 'statfile = ',trmmstatfile

     OPEN(unit=3, file=trmmstatfile, status='new',form='formatted')

   IF ( trmmoutfile .eq. 'none' ) THEN
    trmmoutfile = filenm(1:length)//number//'.out'
   ENDIF

   count = 0
   
6   INQUIRE(file=trmmoutfile, exist=file_exist)

   IF ( file_exist ) THEN
     write(6,*) 'out file already exists!  I will not overwrite -- changing name!'
     count = count + 1
     write(dup, '(i2.2)')  count
     
     IF ( count .ge. 100 ) STOP

     trmmoutfile = prefix(ibeg:iend)//number//'.'//dup//'.out'

     GOTO 6
   ENDIF

   write(6,*) 'outfile = ',trmmoutfile

     OPEN(unit=91, file=trmmoutfile, status='new',form='formatted')
     
     luno = 91
     lune = luno
     lun = luno
    
!   count = 0
!5  INQUIRE(file=v5dfilename, exist=file_exist)
   
!   IF ( file_exist ) THEN
!     write(6,*) 'vis5d file already exists!  I will not overwrite -- changing name!'
!     count = count + 1
!     write(dup, '(i2.2)')  count
!     
!     IF ( count .ge. 100 ) STOP

!     v5dfilename = prefix(ibeg:iend)//number//'.'//dup//'.v5d'

!     GOTO 5
!   ELSE
!     write(6,*) 'COM2TRMM:  CREATE NEW VIS5D FILE: ', v5dfilename ! prefix(1:length)//number
!     gx     = GET_VARIABLE_INDEX(gd, 'XC')
!     gy     = GET_VARIABLE_INDEX(gd, 'YC')
!     gz     = GET_VARIABLE_INDEX(gd, 'ZC')
!     call Vis5DInit( gd, gd%var(gx), gd%var(gy), gd%var(gz) , 1)
!   ENDIF
   
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COM2TRMM:  Problem setting TVIS5D'  
    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COM2TRMM:  Problem setting TIME'  



!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

  CALL GET_VARIABLE(gd, 'NX',       nx)
  CALL GET_VARIABLE(gd, 'NY',       ny)
  CALL GET_VARIABLE(gd, 'NZ',       nz)
  CALL GET_VARIABLE(gd, 'DX',       dx)
  CALL GET_VARIABLE(gd, 'DY',       dy)
  CALL GET_VARIABLE(gd, 'DZ',       dz)
  CALL GET_VARIABLE(gd, 'DT',       dt)
  CALL GET_VARIABLE(gd, 'NSMALL',   nsmall)
  CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
  CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)
  CALL GET_VARIABLE(gd, 'IPELEC',   ipelec)
  CALL GET_VARIABLE(gd, 'IPCONC',   ipconc)

  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  
  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COM2TRMM:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COM2TRMM:  Problem setting NG'

! Allocate SOLVER scratch space
  
  m = (nx+2*ng)*(ny+2*ng)*(nz+2*ng)
  allocate ( st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns) )
  IF ( hack ) THEN 
   allocate ( flsh(nx,ny,nz,3,2) )
   flsh = 0.0
  ENDIF
  
  allocate ( ft(m,11) ) ! temporary arrays for solver

  tprt = start + tprint   - 1
  tstt = start + tstat    - 1
  tv5d = start 
!  ntstart = start + dt

  this = thistory * (1 + start/thistory) - 1
  tv5d = start + tvis5d                  - 1
  
!  DO it = 1,100
!    IF ( gridtimes(it+1) .eq. 0 ) THEN
!      stop = gridtimes(it)
!      EXIT
!    ENDIF
!  ENDDO
!-----------------------------------------------------------------------------
! PRINT OUT PARAMS

  write(lun,"(1x,80('-'))")
  write(lun,*)
  write(lun,*) 'COMMAS_1.0 INPUT DATA'
  write(lun,*)
  write(lun,"(1x,80('-'))")
  write(lun,*)
  write(lun,*) 'SIMULATION NAME            = ', prefix(1:length)//number
  write(lun,*) 'SIMULATION START TIME      = ', start
  write(lun,*) 'SIMULATION STOP  TIME      = ', stop
  write(lun,*) 'TIME STEP                  = ', dt
  write(lun,*) 'NUMBER OF TIME STEPS       = ', (stop-start)/dt
  write(lun,*) 'NUMBER OF SCALAR VARS      = ', ns
  write(lun,*)
  write(lun,*) 'No. of SMALL TIME STEPS    = ', nsmall
  write(lun,*)
  write(lun,*) 'PRINT INTERVAL             = ', tprint
  write(lun,*) 'HISTORY DUMP INTERVAL      = ', thistory
  write(lun,*) 'VIS5D DUMP INTERVAL        = ', tvis5d
  write(lun,*)
  write(lun,*) 'X-GRID MOTION (M/S)        = ', ugrid
  write(lun,*) 'Y-GRID MOTION (M/S)        = ', vgrid
  write(lun,*)


!-----------------------------------------------------------------------------
! Create run lock file to help facilitate parallel runs


!-----------------------------------------------------------------------------
! GET the variables needed for the integration 

 u      = GET_VARIABLE_INDEX(gd, 'U')
 v      = GET_VARIABLE_INDEX(gd, 'V')
 w      = GET_VARIABLE_INDEX(gd, 'W')
 pi     = GET_VARIABLE_INDEX(gd, 'PI')
 km     = GET_VARIABLE_INDEX(gd, 'KM')
 s      = GET_VARIABLE_INDEX(gd, 'TH')
 dbz    = GET_VARIABLE_INDEX(gd, 'DBZ')
 wz     = GET_VARIABLE_INDEX(gd, 'WZ')
 elec   = GET_VARIABLE_INDEX(gd, 'EX')
 cion   = GET_VARIABLE_INDEX(gd, 'CPIONINIT')
 muz     = GET_VARIABLE_INDEX(gd, 'MUPOSZC')
 
 uinit  = GET_VARIABLE_INDEX(gd, 'UINIT')
 vinit  = GET_VARIABLE_INDEX(gd, 'VINIT')
 winit  = GET_VARIABLE_INDEX(gd, 'WINIT')
 piinit = GET_VARIABLE_INDEX(gd, 'PIINIT')
 kminit = GET_VARIABLE_INDEX(gd, 'KMINIT')
 sinit  = GET_VARIABLE_INDEX(gd, 'THINIT')

 gx     = GET_VARIABLE_INDEX(gd, 'XC')
 gy     = GET_VARIABLE_INDEX(gd, 'YC')
 gz     = GET_VARIABLE_INDEX(gd, 'ZC')
 precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')


   IF ( microphys(1:5) .eq. 'ICE10' .or. microphys(1:1) .eq. 'Z' .or. microphys(1:5) .eq. 'WARMZ' ) THEN
     write(3,'(a,3(i6,1x),i3)') 'nx,ny,nz,ns = ',nx,ny,nz,ns
     write(3,'(a,i3,3(1x,f9.2))') 'dt,dx,dy,dz = ',dt,dx,dy,dz
     write(3,'(a,a)') 'micro = ',microphys
     write(3,'(a,2(1x,i3))') 'lqb,lqe= ',lqb,lhab
     write(3,'(a)') ' KZ, ZC, DZC:'
     DO k = 1,nz-1
       write(3,'(1x,i4,f8.0,1x,f10.3)') k, gd%var(gz)%flt1d(k), 1./gd%var(gz+2)%flt1d(k)
     ENDDO
   ENDIF

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Main time step loop

  ibsd = 1
  iesd = nx-1
  jbsd = 1
  jesd = ny-1
   
   count = 1
!   start = gridtimes(1)
   lstt = .true.
   io_flag = .true.
   
   ! set ibsd etc. here for start time
   
   DO it = 1, nt

!----------------------------------------------------------------------
! Time and location of the grid
 
       time     = tarray(it)

!       IF ( time .gt. gridtimes(count) .and. gridtimes(count + 1 ) .eq. 0 ) EXIT
       
       IF ( time .gt. stop ) EXIT
       
       IF ( time .ge. start ) THEN
       
       IF ( time .ge. gridtimes(count+1) .and. gridtimes(count+1) .ne. 0 ) THEN
         count = count + 1
         
         ! set ibsd etc. here
         
       ENDIF

       write(6,*) 'it loop: it,time = ',it,time

      IF ( microphys(1:5) .eq. 'ICE10' .or.          &
           microphys(1:1) .eq. 'Z'     .or.          &
           microphys(1:8) .eq. 'WARMZIEG'   ) THEN
          write(luno,'(a)') '============================================='
          write(luno,'(a)') 'NSTEP, NSTART, NSTOP, TIME'
          write(luno,'(1x,i7,1x,i7,1x,i7,1x,i6,".000")')  time/dt, Max(1,start/dt), (stop)/dt, time
          write(luno,'(a)') '============================================='
!         write(luno,'(a,i7,1x,i6)') 'END OF STEP: nstep, time = ',time/dt,time
       ENDIF

       CALL GRID_READ_NETCDF( gd, prefix(1:length)//number//'.nc', tarray(it), 1 )

!        iv5dwritten(:) = 0
!        lv5dwrite = .true.
!        v5dtimes(itv5d) = time
 

!-----------------------------------------------------------------------
! Call SOLVER to make a time step


!-----------------------------------------------------------------------------
! DTREND the pressure field
!
!
!-----------------------------------------------------------------------------
! PRINTING calls

    
!       CALL PRINT( gd )


!-----------------------------------------------------------------------------
! Vis5D OUTPUT
!
       ibsd = 1
       iesd = nx-1
       jbsd = 1
       jesd = ny-1

      IF ( gridx1(count) .ne. 0 ) ibsd = Max(1,    gridx1(count))
      IF ( gridx2(count) .ne. 0 ) iesd = Min(nx-1, gridx2(count))
      IF ( gridy1(count) .ne. 0 ) jbsd = Max(1,    gridy1(count))
      IF ( gridy2(count) .ne. 0 ) jesd = Min(ny-1, gridy2(count))

 
!        DO k = 1,ns
!          st(1:nx,1:ny,1:nz,k) = gd%var(k+s-1)%flt3d(1:nx,1:ny,1:nz)
!        ENDDO

!  Get flash rate from inits array:

        IF ( elec .gt. wz ) THEN
        
        iflash = 0
        icgn   = 0
        icgp   = 0
        
        DO k = 1,nz-1
        DO j = jbsd,jesd
        DO i = ibsd,iesd
           iflash = iflash + gd%var(elec+ieinit -1)%flt3d(i,j,k)
          IF ( k .eq. 1 ) THEN
            icgn   = icgn + gd%var(elec+ieflshn-1)%flt3d(i,j,k)
            icgp   = icgp + gd%var(elec+ieflshp-1)%flt3d(i,j,k)
          ENDIF
        ENDDO
        ENDDO
        ENDDO
        
        
        write(6,*) ' total flashes, CGN, CGP = ', iflash, icgn, icgp
        write(lun,*) 'COM2TRMM: IC flashes, CGN, CGP = ', iflash - Abs(icgn) - icgp, Abs(icgn), icgp
        write(6,*) ' ibsd,iesd,jbsd,jesd,nz ',ibsd,iesd,jbsd,jesd,nz
        
        ENDIF
        
        write(0,*) 'call TRMMDATDRIVE'

       CALL TRMMDATDRIVE(                              &
                   gd,                                 &
                   gd%var(u),  gd%var(uinit) ,         &            ! U,  UINIT
                   gd%var(v),  gd%var(vinit) ,         &            ! V,  VINIT
                   gd%var(w),  gd%var(winit) ,         &            ! W,  WINIT
                   gd%var(pi), gd%var(piinit),         &            ! PI, PIINIT
                   gd%var(km), gd%var(kminit),         &            ! KM, KINIT
                   gd%var(s),  gd%var(sinit),          &            ! S,  SINIT
                   gd%var(precip),                     &            ! PRECIP
                   gd%var(gx),                         &            ! XCNTR, XEDGE, DXC, DXE 
                   gd%var(gy),                         &            ! YCNTR, YEDGE, DYC, DYE
                   gd%var(gz),                         &            ! ZCNTR, ZEDGE, DZC, DZE
                   dt,                                 &            ! DT
                   ugrid, vgrid,                       &            ! GRID MOTION
                   nsmall,                             &            ! NSMALL
                   microphys,                          &
                   nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS
                   io_flag,                            &            ! IO_FLAG
                   st,ft(1,1),ft(1,2),ft(1,3),ft(1,4), &
                   ft(1,5),ft(1,6),ft(1,7),ft(1,8),    &
                   ft(1,9),ft(1,10),ft(1,11),          &
                   gd%var(dbz), gd%var(wz),            &
                   gd%var(elec), gd%var(cion),         &
                   gd%var(muz),                        &
                   time, tstat, lstt)
        

         CALL stats( nx,ny,nz,ns,                                 &
                     time,dx,dy,dz,                               &
                     gd%var(gx),gd%var(gy),gd%var(gz),            &
                     ft(1,1),gd%var(sinit),gd%var(s),             &
                     gd%var(u), gd%var(v), gd%var(w),             &
                     gd%var(uinit) ,gd%var(vinit), gd%var(winit), &
                     gd%var(pi), gd%var(piinit),                  &
                     gd%var(km), gd%var(kminit),                  &
                     ibsd, iesd, jbsd, jesd)

        
!        CALL trmmdatdrive(gd,                                 &
!                    gd%var(u),  gd%var(uinit) ,         &            ! U,  UINIT
!                    gd%var(v),  gd%var(vinit) ,         &            ! V,  VINIT
!                    gd%var(w),  gd%var(winit) ,         &            ! W,  WINIT
!                    gd%var(pi), gd%var(piinit),         &            ! PI, PIINIT
!                    gd%var(km), gd%var(kminit),         &            ! KM, KINIT
!                    gd%var(s),  gd%var(sinit),          &            ! S,  SINIT
!                    gd%var(precip),                     &            ! PRECIP
!                    gd%var(gx),                         &            ! XCNTR, XEDGE, DXC, DXE 
!                    gd%var(gy),                         &            ! YCNTR, YEDGE, DYC, DYE
!                    gd%var(gz),                         &            ! ZCNTR, ZEDGE, DZC, DZE
!                    dt,                                 &            ! DT
!                    ugrid, vgrid,                       &            ! GRID MOTION
!                    nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS
!                    st,                                 &
!                    ft(1,1),ft(1,2),ft(1,3),ft(1,4),    &
!                    gd%var(dbz), gd%var(wz), gd%var(elec) )
        
        write(0,*) 'done TRMMDATDRIVE'



      IF ( microphys(1:5) .eq. 'ICE10' .or.           &
           microphys(1:1) .eq. 'Z'     .or.           &
           microphys(1:8) .eq. 'WARMZIEG'   ) THEN
        write(luno,'(a)') '============================================='
        write(luno,'(a,i7,1x,i6)') 'END OF STEP: nstep, time = ',time/dt,time
      ENDIF

!        iv5dwritten(:) = 0

   

        ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO 
   
  DEALLOCATE(st)
  DEALLOCATE(ft)
 

!-----------------------------------------------------------------------------
! Print out date and time

  write(luno,'(a)') ' Integration Done'


 STOP
 END
 
 

