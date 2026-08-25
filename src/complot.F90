!-----------------------------------------------------------------------------
!
! COMMAS-1.0
!
! March 2006
!
!-----------------------------------------------------------------------------

 PROGRAM COMPLOT

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
  USE NAMELIST_MODULE, lowdbz_tmp => lowdbz
  USE VIS5D_MODULE
  USE INDEX_MODULE
  USE PLOT_MODULE
  use dualpara,         only: model_dsd,wgfac ,tosoak, freezesoak, conserveice
   use takcommon,       only: lfmax
   use comm2,           only: sx,sxf,sr

  implicit none

!-----------------------------------------------------------------------------
! GRID DEFINITION
 
  TYPE(GRID) :: gd, ge

!-----------------------------------------------------------------------------
! COMMAND LINE VARIABLES

  character(LEN = 120) :: run_file
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
  integer :: xtra
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

!  real, allocatable    :: st(:,:,:,:)         ! These are big scratch arrays for solver
  real, allocatable    :: flsh(:,:,:,:,:)
  real, allocatable    :: ft(:,:)
  real, allocatable :: z1d4(:,:)
  integer, allocatable :: tarray(:)
  real, allocatable ::  stinit(:,:)
  real, allocatable :: precip_old(:,:,:)

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
  
  integer i,j,k, it, nt, count
  integer ix,jy,kz

  character(LEN=4) number
  character(LEN=2) dup

  character(LEN = 255):: v5dfldstmp
  
  integer :: iunit = 90

  integer  :: iskip = 1, istart = 1

  real, allocatable :: lgtth(:,:)

      integer, parameter :: chans    = 1
      integer, parameter :: chansp   = 2
      integer, parameter :: chansn   = 3
      integer, parameter :: cgchans  = 4
      integer, parameter :: cgchansp = 5
      integer, parameter :: cgchansn = 6
      integer, parameter :: icchans  = 7
      integer, parameter :: icchansp = 8
      integer, parameter :: icchansn = 9
      integer, parameter :: cinit    = 10
      integer, parameter :: cginit   = 11
      integer, parameter :: icinit   = 12
      integer, parameter :: cgninit  = 13
      integer, parameter :: cgpinit  = 14
      
      integer :: loccur
      real    :: xinit
      real    :: xsegn
      real    :: xsegp
      
      double precision    :: volpos, volneg
      real    :: dv, chgrat

  character(len=6) stime
  character(LEN=120) :: ncfile
  integer :: iframe,istat

  NAMELIST /dualpol/ do_dualpol, &
                     dirscatt,   &
                     MPflg, MFflg, wgfac, &
                     dualpol_fulldomain, &
                     tosoak, freezesoak, conserveice


!-----------------------------------------------------------------------------
! Need 3D.RUN file in order to get run parameters - get from the command line

  IF( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN

   write(0,*) 'COMPLOT:  INCORRECT ARGUMENTS ON COMMAND LINE:  NEED 3D.RUN FILENAME'
   write(0,*) 'COMPLOT:  INCORRECT ARGUMENTS ON COMMAND LINE:  EXITING RUN!...'
   call exit(1)

  ELSE

   CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)

!   write(6,*) 'COMPLOT:  Input 3D.RUN file is:  ', run_file(1:length)

   IF( status .ne. 0 ) THEN

     write(0,*) 'COMPLOT:  COULD NOT RETRIEVE 1ST COMMAND LINE ARG:  EXITING RUN!...', status
     write(0,*) 'COMPLOT:  DOES NOT HAVE A 3D.RUN FILE TO READ!'
     call exit(1)

   ENDIF

  ENDIF
  
!-----------------------------------------------------------------------------
! Read in 3d.run file...

  INQUIRE(file=run_file(1:length), exist=file_exist)

  IF( .NOT. file_exist ) THEN

    write(0,*) 'COMPLOT:  INPUT 3D.RUN FILE:  ', run_file(1:length), ' DOES NOT EXIST!!!'
    write(0,*) 'COMPLOT:  DOES NOT HAVE A 3D.RUN FILE TO READ, EXITING'
    stop

  ENDIF

  open(15,file=run_file(1:length),status='old',form='formatted')
  rewind(15)
  write(0,*) 'read plot nml'
  read(15,NML=plot)
  rewind(15)
  write(0,*) 'read microparams nml'
  read(15,NML=micro_params)
  rewind(15)
  write(0,*) 'read dualpol nml'
  read(15,NML=dualpol,iostat=istat)
      IF ( istat .ne. 0 .and. my_rank == 0 ) THEN
        write(0,*) 'Problem reading dualpol namelist -- not found or bad token'
      ENDIF
      write(0,*) 'dirscatt = ',dirscatt

  close(15)

  IF ( recalcdbz .ge. 1 ) iusewetgraupel = recalcdbz
  IF ( cxmin_plot > 0.0 ) cxmin = cxmin_plot
  IF ( zxmin_plot > 0.0 ) zxmin = zxmin_plot
  
  IF ( do_dualpol ) THEN
   call model_dsd(MPflg,cnor,cnos,cnohl,cnoh,rho_qs,rho_qhl,rho_qh)
  ENDIF

! READ IN RUN MODEL NAMELIST

  IF( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'

  write(number, '(a,i3.3)') '.', member

!-----------------------------------------------------------------------------
!  

! READ IN RUN ATTRIBUTES
  
  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
  


  IF( .not. READ_NAMELIST(run_file(1:length),'ICE10_PARAMS')  ) THEN
    IF ( my_rank == 0 ) print *, 'COMMAS:  PROBLEM READING NAMELIST ICE10_PARAMS'
  ENDIF

  IF( (.not. READ_NAMELIST(run_file(1:length),'OUTPUT_OPTIONS'))  ) THEN
    IF ( my_rank == 0 ) write(0,*) 'COMMAS:  PROBLEM READING NAMELIST OUTPUT_OPTIONS'
  ENDIF

! READ IN GRID PARAMETERS


!  v5dfldstmp = v5dflds



!-----------------------------------------------------------------------------
! IF START TIME < 0, then initialize a single grid...

   CALL STRING_LIMITS(prefix, ibeg, iend)
!   v5dfilename = prefix(ibeg:iend)//number//'.v5d'


!-----------------------------------------------------------------------------
! READ IN THE DATA


    CALL STRING_LIMITS(prefix, ibeg, iend)
    length = iend - ibeg + 1

     write(stime,   '(i6.6)') 0

    IF ( historyoutput == 3 ) THEN ! netcdf output separate file for each history time

       write(stime,'(i6.6)') start
       ncfile = prefix(1:length)//'.'//stime//'.nc'

    ELSE
         ncfile = prefix(1:length)//trim(number)//'.nc'
    ENDIF

    write(6,*) 'COMPLOT:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//number
  
    CALL GRID_INFO_NETCDF( ncfile, nt, microphys)
    microp = microphys
    
    allocate( tarray(nt) )
    
    CALL GRID_INFO_NETCDF( ncfile, nt, microphys, tarray)
    
    
    CALL GRID_READ_NETCDF( gd, ncfile, tarray(1), nx_or_nxend=.true. )

    write(6,*) 'COMPLOT:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//number

    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMPLOT:  Problem setting TIME'  
!    IF( .not. SET_VARIABLE(gd,'UGRID',     ugrid)   ) write(6,*) 'COMPLOT:  Problem setting UGRID'  
!    IF( .not. SET_VARIABLE(gd,'VGRID',     vgrid)   ) write(6,*) 'COMPLOT:  Problem setting VGRID'
    CALL GET_VARIABLE(gd, 'UGRID', ugrid)
    CALL GET_VARIABLE(gd, 'VGRID', vgrid)
    CALL GET_VARIABLE(gd, 'THISTORY', thistory)
!    IF( .not. SET_VARIABLE(gd,'TRESTART',  trestart)) write(6,*) 'COMPLOT:  Problem setting TRESTART'
!    IF( .not. SET_VARIABLE(gd,'THISTORY',  thistory)) write(6,*) 'COMPLOT:  Problem setting THISTORY'
!    IF( .not. SET_VARIABLE(gd,'TVIS5D',    tvis5d)  ) write(6,*) 'COMPLOT:  Problem setting TVIS5D'
    IF( .not. SET_VARIABLE(gd,'TPRINT',    tprint)  ) write(6,*) 'COMPLOT:  Problem setting TPRINT'
    IF( .not. SET_VARIABLE(gd,'TSTAT',     tstat)   ) write(6,*) 'COMPLOT:  Problem setting TSTAT'
    IF( .not. SET_VARIABLE(gd,'TIME_STOP', stop)    ) write(6,*) 'COMPLOT:  Problem setting STOP'
    IF( .not. SET_VARIABLE(gd,'DT',        dt)      ) write(6,*) 'COMPLOT:  Problem setting DT'
   
!    CALL STRING_LIMITS(v5dfldstmp, i, j)
!    IF( .not. SET_ATTRIBUTE(gd,'V5DFIELDS',  v5dfldstmp)) write(6,*) 'COMPLOT:  Problem setting V5DFIELDS'

    time = start


!--------------------------------------------------------------------------------------
! Open VIS5D file  

!  IF( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *, 'COMMAS_INIT:  PROBLEM READING NAMELIST'
    
    IF ( start .gt. tarray(nt) ) THEN
      write(0,*) 'Starting time is greater than the last record time! STOP!'
      STOP
    ENDIF
    
    IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMPLOT:  Problem setting TIME'  



!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

  CALL GET_VARIABLE(gd, 'BCX',       bcx)
  CALL GET_VARIABLE(gd, 'BCY',       bcy)

  CALL GET_VARIABLE(gd, 'NXEND',       nx)
  CALL GET_VARIABLE(gd, 'NYEND',       ny)
  CALL GET_VARIABLE(gd, 'NZEND',       nz)

  IF( .not. SET_VARIABLE(gd,'NX',      nx)      ) write(6,*) 'COMPLOT:  Problem setting NX'
  IF( .not. SET_VARIABLE(gd,'NY',      ny)      ) write(6,*) 'COMPLOT:  Problem setting NX'

  CALL GET_VARIABLE(gd, 'DX',       dx)
  CALL GET_VARIABLE(gd, 'DY',       dy)
  CALL GET_VARIABLE(gd, 'DZ',       dz)
  CALL GET_VARIABLE(gd, 'DX_STRETCH',       DX_STRETCH)
  CALL GET_VARIABLE(gd, 'DY_STRETCH',       DY_STRETCH)
  CALL GET_VARIABLE(gd, 'DT',       dt)
  CALL GET_VARIABLE(gd, 'NSMALL',   nsmall)
  CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
  CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)
  CALL GET_VARIABLE(gd, 'IPELEC',   ipelec)
  CALL GET_VARIABLE(gd, 'IPCONC',   ipconc)

  ns = size(gd%var(:)) - GET_VARIABLE_INDEX(gd,'TH') + 1
  
  IF( .not. SET_VARIABLE(gd,'NSCALAR',  ns)      ) write(6,*) 'COMPLOT:  Problem setting NSCALAR'
  IF( .not. SET_VARIABLE(gd,'NG',       ng)      ) write(6,*) 'COMPLOT:  Problem setting NG'

! Allocate SOLVER scratch space
  
  m = (nx+2*ng)*(ny+2*ng)*(nz+2*ng)
!  allocate ( st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,ns) )
  allocate ( stinit(nz,ns) )
  IF ( hack ) THEN 
   allocate ( flsh(nx,ny,nz,3,2) )
   flsh = 0.0
  ENDIF
  
  allocate ( ft(m,4) ) ! temporary arrays for solver

  allocate ( precip_old(-ng+1:nx+ng,-ng+1:ny+ng,nprecip) )

  tprt = start + tprint   - 1
  tstt = start + tstat    - 1
!  tv5d = start 
  ntstart = start + dt

  this = thistory * (1 + start/thistory) - 1
!  tv5d = start + tvis5d                  - 1

!-----------------------------------------------------------------------------
! PRINT OUT PARAMS


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

 ! check for extra diagnostic 3d arrays
  IF ( elec > wz .and.  elec - wz > 1 ) THEN ! elec is on and extra arrays exist
     xtra = wz + 1
     nxtra = elec - wz - 1
  ELSEIF ( u > wz + 1 ) THEN ! no elec and extra arrays exist
     xtra = wz + 1
     nxtra = u - wz - 1
  ELSE ! no extra arrays
     xtra = wz
     nxtra = 0
  ENDIF

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

     DO i = precip,precip+nprecip - 1
!       write(0,*) 'load precip_old i = ',i, i - precip + 1
        precip_old(-ng+1:nx+ng,-ng+1:ny+ng,i - precip + 1) = gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng)
     ENDDO

!-----------------------------------------------------------------------------
    
    OPEN(unit=iunit,file='inplt',form='formatted',status='old')

    IF ( statout ) THEN
      OPEN(unit=3,file=trim(prefix)//'plot.statout.txt',form='formatted',status='unknown')
      allocate( lgtth(nz,14) )
      lgtth(:,:) = 0.0

      write(3,'(a,3(i6,1x),i3)') 'nx,ny,nz,ns = ',nx,ny,nz,ns
      write(3,'(a,1x,4(1x,f9.2))') 'dt,dx,dy,dz = ',dt,dx,dy,dz
      write(3,'(a,a)') 'micro = ',microphys
      write(3,'(a,2(1x,i3))') 'lqb,lqe= ',lqb,lhab
      write(3,'(a)') ' KZ, ZC, DZC:'
      DO k = 1,nz-1
        write(3,'(1x,i4,f8.0,1x,f10.3)') k, gd%var(gz)%flt1d(k), 1./gd%var(gz+2)%flt1d(k)
      ENDDO

    ENDIF


      IF ( microphys(1:3) .eq. 'TAK' ) THEN
         IF ( sr(1) == 0.d0 ) THEN
           call takinitbin
           write(0,*) 'call takinitbin, sx(1) = ',sx(1),sxf(1,1),sxf(1,2),lfmax
         ENDIF
      ENDIF

   allocate( z1d4(nzend,4) )
   
     DO i = 1,4
      DO k = 1,nz
!       z1d4(kzbeg-1+k,i) = gd%var(gz+i-1)%flt1d(k)
       z1d4(k,i) = gd%var(gz+i-1)%flt1d(k)
      ENDDO
     ENDDO


    IF ( lchgratevolume ) THEN
      OPEN(unit=13,file=trim(prefix)//'plot.chgratevolume.txt',form='formatted',status='unknown')
      write(3,'(a)') 'Time (s), Time (min), volpos, volneg'
    ENDIF
    
    IF ( member .gt. 0 ) THEN
!        inameps = len(outname)
!        CALL strlnth(outname,inameps)
        write(0,*) 'member = ',member
        write(0,*) 'number = ',number
        CALL STRING_LIMITS(outname, ibeg, iend)
        outname = outname(ibeg:iend)//number
     ENDIF
!-----------------------------------------------------------------------------


   IF ( ifile2 .le. -1 ) nt = 1
 
       IF ( historyoutput == 3 ) THEN 
         nt = (stop - start)/thistory + 1
!         iskip = 1
!         IF ( tplot < 999999 ) iskip = Max( 1,tplot/thistory )
         IF ( thistorysp > 0 ) THEN
           nt = (stop - start)/thistorysp + 1
         ENDIF
       ENDIF

   IF ( tplot < 999999 .and. tplot .ne. thistory ) THEN
    iskip = Max( 1,tplot/thistory )
    istart = Max( 1, (start-tarray(1))/thistory + 1 )
    print*, 'iskip ,tplot, thistory= ',iskip, tplot, thistory 
    print*, 'start, istart = ',start, istart,nt
   ENDIF
   
!-----------------------------------------------------------------------------
! Main time step loop

   iframe = 0
   DO it = istart, nt, iskip

!----------------------------------------------------------------------
! Time and location of the grid
 
       IF ( historyoutput == 3 ) THEN 
         time = start + (it - 1)*thistory
         IF ( thistorysp > 0 ) THEN
           time = start + (it - 1)*thistorysp
         ENDIF
       ELSE
         time = tarray(it)
       ENDIF
       
       write(0,*) 'time,start,stop = ',time,start,stop

       IF ( ifile2 .ge. 1 ) THEN
       
         IF ( time .gt. stop ) EXIT
       
         IF ( time .lt. start ) CYCLE
       
       ENDIF

      iframe = iframe + 1
       
!       IF ( time .ge. start ) THEN

       write(0,*) 'it loop: it,time = ',it,time

  !     IF ( statout .and. (it/2)*2 .eq. it ) CYCLE
       
       IF ( it > 1 ) THEN 
        IF ( historyoutput == 3 ) THEN ! netcdf output separate file for each history time
          write(stime,'(i6.6)') time
          ncfile = prefix(1:length)//'.'//stime//'.nc'
            INQUIRE(file=trim(ncfile), exist=file_exist)
            IF( .NOT. file_exist ) THEN
              CYCLE
            ENDIF
         CALL GRID_READ_NETCDF( gd, ncfile, time, 1 )
        ELSE
          IF ( ( it == istart .and. istart > 1 ) .or.  &
               ( it > istart .and. iskip > 1 ) ) THEN
            CALL GRID_READ_NETCDF( gd, ncfile, tarray(it-1), 1 )
            DO i = precip,precip+nprecip - 1
!       write(0,*) 'load precip_old i = ',i, i - precip + 1
             precip_old(-ng+1:nx+ng,-ng+1:ny+ng,i - precip + 1) = gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng)
           ENDDO
          
          ENDIF
          CALL GRID_READ_NETCDF( gd, ncfile, tarray(it), 1 )
        ENDIF
       ENDIF

        lstt = .false.
        IF ( statout .and. Mod(time,tstat) == 0 ) lstt = .true.
       
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

!
! Copy variables in XYZ3D array for convenience (can get rid of this eventually, but requires recoding)
! Fill in scalar values at edge so that we can plot all the way to nx or ny
        DO k = 1,ns
!          st(1:nx,1:ny,1:nz,k) = gd%var(k+s-1)%flt3d(1:nx,1:ny,1:nz)
          DO kz = 1,nz-1
           DO jy = 1,ny-1
            gd%var(k+s-1)%flt3d(nx,jy,kz) = gd%var(k+s-1)%flt3d(nx-1,jy,kz)
           ENDDO
          ENDDO
          DO kz = 1,nz-1
           DO ix = 1,nx
            IF ( bcy .ne. 2 ) THEN
            gd%var(k+s-1)%flt3d(ix,ny,kz) = gd%var(k+s-1)%flt3d(ix,ny-1,kz)
            gd%var(k+s-1)%flt3d(ix,0 ,kz) = gd%var(k+s-1)%flt3d(ix,   1,kz)
            ELSE
            gd%var(k+s-1)%flt3d(ix,ny,kz) = gd%var(k+s-1)%flt3d(ix,1,kz)
            gd%var(k+s-1)%flt3d(ix,0,kz) = gd%var(k+s-1)%flt3d(ix,ny-1,kz)
            ENDIF
           ENDDO
          ENDDO
        ENDDO

! fill in edge values for dbz
          DO kz = 1,nz-1
           DO jy = 1,ny-1
            gd%var(dbz)%flt3d(nx,jy,kz) = gd%var(dbz)%flt3d(nx-1,jy,kz)
           ENDDO
          ENDDO
          DO kz = 1,nz-1
           DO ix = 1,nx
            IF ( bcy .ne. 2 ) THEN
            gd%var(dbz)%flt3d(ix,ny,kz) = gd%var(dbz)%flt3d(ix,ny-1,kz)
            gd%var(dbz)%flt3d(ix,0,kz) = gd%var(dbz)%flt3d(ix,1,kz)
            ELSE
            gd%var(dbz)%flt3d(ix,ny,kz) = gd%var(dbz)%flt3d(ix,1,kz)
            gd%var(dbz)%flt3d(ix,0,kz) = gd%var(dbz)%flt3d(ix,ny-1,kz)
            ENDIF
           ENDDO
          ENDDO



        IF ( hack ) THEN
        
    ! code to handle case where history interval was less than 60s and want to generate new 60-second time-height data
          flsh(1:nx,1:ny,1:nz,1,1) = flsh(1:nx,1:ny,1:nz,1,1) + gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz)
          flsh(1:nx,1:ny,1:nz,2,1) = flsh(1:nx,1:ny,1:nz,2,1) + gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz)
          flsh(1:nx,1:ny,1:nz,3,1) = flsh(1:nx,1:ny,1:nz,3,1) + gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz)

          gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz) = flsh(1:nx,1:ny,1:nz,1,1)
          gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz) = flsh(1:nx,1:ny,1:nz,2,1)
          gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz) = flsh(1:nx,1:ny,1:nz,3,1)

!  old code to handle the bug where the flash arrays didn't get reset to zero after history dumps
!          flsh(1:nx,1:ny,1:nz,1,1) = gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz)
!          flsh(1:nx,1:ny,1:nz,2,1) = gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz)
!          flsh(1:nx,1:ny,1:nz,3,1) = gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz)

!          gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz) = gd%var(elec+ieflshn-1)%flt3d(1:nx,1:ny,1:nz) - flsh(1:nx,1:ny,1:nz,1,2)
!          gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz) = gd%var(elec+ieflshp-1)%flt3d(1:nx,1:ny,1:nz) - flsh(1:nx,1:ny,1:nz,2,2)
!          gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz) = gd%var(elec+ieinit -1)%flt3d(1:nx,1:ny,1:nz) - flsh(1:nx,1:ny,1:nz,3,2)
        
        ENDIF
        
!-----------------------------------------------------------------------------
! Stats OUTPUT
!  Note that this portion of code is customized (hacked) for a particular purpose.
!  Normally stats output will come directly from the code!
!
        IF ( statout .and. lstt ) THEN
        
        IF ( associated(gd%var(elec+iez-1)%flt3d) ) THEN
         loccur = 1
          k = 1
   ! treat all as IC for current hack
   !       DO j=1,ny-1
   !        DO i=1,nx-1
   !          IF ( Abs(gd%var(elec+ieflshn-1)%flt3d(i,j,k)) .ge. 0.5 ) loccur = 2 
   !          IF ( Abs(gd%var(elec+ieflshp-1)%flt3d(i,j,k)) .ge. 0.5 ) loccur = 3 
   !        ENDDO
   !       ENDDO

           DO k=1,nz-1
            DO j=1,ny-1
             DO i=1,nx-1
              xinit = gd%var(elec+ieinit-1)%flt3d(i,j,k)
              xsegn = Abs(gd%var(elec+ieflshn-1)%flt3d(i,j,k))
              xsegp = gd%var(elec+ieflshp-1)%flt3d(i,j,k)
              
 
              lgtth(k,chans) = lgtth(k,chans) + xsegn + xsegp
              lgtth(k,chansp) = lgtth(k,chansp) + xsegp
              lgtth(k,chansn) = lgtth(k,chansn) + xsegn
              
             IF ( loccur .eq. 1 ) THEN
               lgtth(k,icchans)  = lgtth(k,icchans)  + xsegp + xsegn
               lgtth(k,icchansp) = lgtth(k,icchansp) + xsegp
               lgtth(k,icchansn) = lgtth(k,icchansn) + xsegn
             ELSEIF ( loccur .ge. 2 ) THEN
               lgtth(k,cgchans)  = lgtth(k,cgchans) + xsegp + xsegn
               lgtth(k,cgchansp) = lgtth(k,cgchansp) + xsegp
               lgtth(k,cgchansn) = lgtth(k,cgchansn) + xsegn
             ENDIF

             lgtth(k,cinit) = lgtth(k,cinit) + xinit

             IF ( loccur .eq. 1 ) lgtth(k,icinit)  = lgtth(k,icinit) + xinit
             IF ( loccur .eq. 2 ) THEN
               lgtth(k,cgninit) = lgtth(k,cgninit) + xinit
               lgtth(k,cginit) = lgtth(k,cginit) + xinit
             ENDIF
             IF ( loccur .eq. 3 ) THEN
               lgtth(k,cgpinit) = lgtth(k,cgpinit) + xinit
               lgtth(k,cginit) = lgtth(k,cginit) + xinit
             ENDIF
             
             ENDDO
            ENDDO
           ENDDO
        
       write(3,'(a,i6)') 'Lightning sums for time = ',time
       write(3,'(3a)') ' Time Altitude Channels chanp chann ',   &
        'cgchan cgchanp cgchann icchan icchnp icchann ',         &
        'cinit cginit icinit cgninit cgpinit'

        DO k=nz-1,1,-1
            write(3,'(1x,f8.1,1x,f8.4,14(1x,f7.0))') &
             float(time)/60.0,0.001*gd%var(gz)%flt1d(k), (lgtth(k,i), i=1,14)
        ENDDO
        
         lgtth(:,:) = 0.0
        
        ENDIF

        write(0,*) 'call stats'
        write(6,*) 'call stats6'
        
        
         CALL stats( nx,ny,nz,ns,                                 &
                     time,dx,dy,dz,                               &
                     gd%var(gx),gd%var(gy),gd%var(gz),            &
                     ft(1,1),gd%var(sinit),gd%var(s),             &
                     gd%var(u), gd%var(v), gd%var(w),             &
                     gd%var(uinit) ,gd%var(vinit), gd%var(winit), &
                     gd%var(pi), gd%var(piinit),                  &
                     gd%var(km), gd%var(kminit),                  &
                     z1d4,                                        &
                     1, nx-1, 1, ny-1)
          
           flsh(1:nx,1:ny,1:nz,1:3,1) = 0
           
          CYCLE
        ENDIF


        IF ( lchgratevolume .and. associated(gd%var(elec+iez-1)%flt3d)  ) THEN
           
           volpos = 0.0d0
           volneg = 0.0d0
           
           DO k=1,nz-1
            DO j=1,ny-1
             DO i=1,nx-1
               dv = 1.e-9/(gd%var(gx+2)%flt1d(i) * gd%var(gy+2)%flt1d(j) * gd%var(gz+2)%flt1d(k))
               IF ( i==1 .and. j==1 .and. k==1 ) THEN
                write(0,*) 'dv,gx,gy,gz = ',dv,gd%var(gx+2)%flt1d(i),gd%var(gy+2)%flt1d(j), gd%var(gz+2)%flt1d(k)
               ENDIF
               chgrat = gd%var(elec+icghis-1)%flt3d(i,j,k)  
               IF ( chgrat .gt. chgratevolume ) THEN
                 volpos = volpos + dv
               ELSEIF ( chgrat .lt. -chgratevolume ) THEN
                 volneg = volneg + dv
               ENDIF
             ENDDO
            ENDDO
           ENDDO
        
          write(13,'(1x,i6,a,f9.2,a,f9.4,a,f9.4)') time, ',', time/60., ',', volpos, ',', volneg
          CYCLE
        ENDIF


        
        IF ( .not. ( lchgratevolume .or. statout ) ) THEN

        write(0,*) 'call cplotba'
        
        
        CALL cplotba(gd,                                &
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
                    ugrid, vgrid,                       &            ! GRID MOTION
                    nx, ny, nz, ns,                     &            ! NX,NY,NZ,NS
!                    st,                                 &
                     gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,gd%var(s)%index),  &
                    ft(1,1),ft(1,2),ft(1,3),ft(1,4),    &
                    gd%var(dbz), gd%var(wz), gd%var(elec), gd%var(xtra), &
                    stinit,                                &
                    iunit,time,member,dx_stretch,dy_stretch,dz,dt,iframe,precip_old)

        
!        itv5d = itv5d + 1
        
        write(6,*) 'done cplotba'

        ENDIF
       
        IF ( hack ) THEN
        
     !     flsh(1:nx,1:ny,1:nz,1,2) = flsh(1:nx,1:ny,1:nz,1,1)
     !     flsh(1:nx,1:ny,1:nz,2,2) = flsh(1:nx,1:ny,1:nz,2,1)
     !     flsh(1:nx,1:ny,1:nz,3,2) = flsh(1:nx,1:ny,1:nz,3,1)
          
        ENDIF

     DO i = precip,precip+nprecip - 1
!       write(0,*) 'load precip_old i = ',i, i - precip + 1
        precip_old(-ng+1:nx+ng,-ng+1:ny+ng,i - precip + 1) = gd%var(i)%flt2d(-ng+1:nx+ng,-ng+1:ny+ng)
     ENDDO


!        iv5dwritten(:) = 0

   

!        ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
 
END DO 
   
  close(iunit)

    IF ( statout ) THEN
      CLOSE(unit=3)
    ELSEIF ( lchgratevolume ) THEN
      CLOSE(unit=13)
    ELSE
      call clsgks
    ENDIF
  

!  DEALLOCATE(st)
  DEALLOCATE(ft)
  DEALLOCATE(stinit)
 

!-----------------------------------------------------------------------------
! Print out date and time


! STOP
 END
