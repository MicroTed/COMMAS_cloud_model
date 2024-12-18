#ifndef GETMEM
! need to define GETMEM to use the external memory usage routine
  subroutine get_mem_used(mem)
  implicit none
  integer mem
  
  mem = 0
  
  return
  end

#endif
!----------------------------------------------------------------------------------------
! This program combines netCDF data from different processors from
! George Bryan's cm1 numerical model. (Now for COMMAS7)
!
! Provided by:  Daniel Kirshbaum, University of Reading
! 
! Last modified:  2 April 2008
!
! Updated:  Early Sept 2008 by Lou Wicker, National Severe Storms Lab
!           ==> Converted it to F90 freeform style (needs compiling for 132 char lines)
!           ==> Created subroutines for reading in variables
!           ==> Made it more compatible with cm1r12 (vertical staggering of tke/kh/kv)
!           OTHER CHANGES
!           * I added in dumping out all of the coordinates / horiz grid stretching not yet supported
!           * There is no terrain coordinate dump - its easy to add in
!----------------------------------------------------------------------------------------
!
! Spring, 2009: adapted by M. Wandishin and E. Mansell for use with COMMAS.  
!     Code is set up to be able to combine either the entire domain or an arbitrary rectangular
!     region, e.g., to remove edge tiles.
!    
!   Now requires HDF5 fortran interface for faster reading of 2D and 3D data (actually faster open/close of the files)
!   Horizontal stretching is supported.
!   Initial file reading by netcdf4 only opens files for tiles along southern and western edges of the domain, which
!     can save a lot of file open/closing by netcdf4.  Even this could be replaced by HDF5 calls, but that may not happen.
!
! erm:
!  compile on sooner (x86_64 cluster):
! ifort -I/home/mansell/opt/netcdf4/include -DNC4 combine_v3.F90 /home/mansell/opt/netcdf4/lib/libnetcdf.a /home/mansell/opt/hdf5/lib/libhdf5_hl.a /home/mansell/opt/hdf5/lib/libhdf5.a -lz -o combine
!
!!
! On Landru, compile with:
! ifort -o combine ~/combine_v3_alt.F90 -L/home/mansell/opt/netcdf4/lib -lnetcdf -L/cluster/local/hdf5/lib -lhdf5_hl -lhdf5 -lz -I/home/mansell/opt/netcdf4/include -DNC4 -g -C 
!

 PROGRAM COMBINE_V3

  USE NAMELIST_MODULE
  USE CLINE_MODULE
  USE FILE_MODULE
  USE COMMASMPI_MODULE
  USE NETCDF
  USE HDF5 ! This module contains all necessary modules 

  implicit none

!-----------------------------------------
! input and auxiliary variables

  integer :: nfile
  character(len=120), dimension(:), allocatable :: fname
  character(len=120) :: npre
  character(len=3)   :: nens, tnum
  integer, allocatable, dimension(:) :: ncid, file_used, hdfid
  INTEGER(HID_T), allocatable, dimension(:) :: hdf5_fid       ! File identifier 

  integer                 :: ntilef_x, &    ! number of tiles in x-direction for entire domain (automatic)
                             ntilef_y       ! number of tiles in y-direction for entire domain (automatic)
!-----------------------------------------
! namelist variables 
  integer                 :: ibeg=-1,  &    ! starting x-index for combining tiles
                             jbeg=-1,  &    ! starting y-index for combining tiles
                             nllc=-1        ! "lower left corner" index.  Calculated from ibeg,jbeg.

   character(len=120) :: prefix_out = " "  ! used for output file name if set.  Otherwise uses original prefix.
   integer            :: ntiles_x = -1, &  ! number of tiles in x-direction of selected region (defaults to entire domain)
                         ntiles_y = -1     ! number of tiles in y-direction of selected region (defaults to entire domain)
   integer            :: begtime=-1, endtime=-1 ! begin/end times for combining.  Used to set file name.
   logical            :: append   = .false.
   integer            :: tstitch  = -1     ! time interval for stitching, if it is to be different from thistory

  NAMELIST /stitcher/               &
                      nllc,         &
                      ntiles_x,     &
                      ntiles_y,     &
                      ibeg, jbeg,   &
                      begtime,      &
                      endtime,      &
                      prefix_out,   &
                      append,       &
                      tstitch

  character(len=120) :: fname1, dump, run_file
  character(len=19)  :: fname2
  character(len=3)   :: smember
  character(len=6)   :: dum1, dum2
  character(len=2)   :: fid
  integer            :: ntbeg, ntend, ntstart
  
  integer            :: ncid1, ncid2
  integer            :: i, j, k, n,lenf, nVars2, ii, jj, vid
  integer            :: ktime
  integer            :: cnt3d
  integer            :: cmode, error
  logical            :: file_exist
  logical            :: one_file_at_a_time, one_time_at_a_time, openclose, one_time, writeout
  logical            :: usehdf5, writehdf5
      real tim1,tim2,timarr(2), tim1a,tim2a
      real tim1r,tim2r,timarrr(2)
      real tim1w,tim2w,timarrw(2)
      real timeread, timewrite
      real commas_dtime, dtime
      integer        :: iunit = 6

!-----------------------------------------
! global vars for the original file

  integer            :: status, nDims, nVars, nGlobalAtts
  integer, dimension(:), allocatable :: GlAttVali, gtype
  character(len=20), dimension(:), allocatable :: GlAttName
  character(len=80), dimension(:), allocatable :: GlAttValc

!-----------------------------------------
! dimensions vars

  integer            :: DimID, DimIDT, DimIDXC, DimIDXE, DimIDYC, DimIDYE, DimIDZC, DimIDZE, DimIDdim2, DimIDdim3
  integer            :: DimIDTapp
  integer, dimension(nf90_max_var_dims) :: DimLen
  character(nf90_max_name) :: DimName

!-----------------------------------------
! attributes vars

  character(nf90_max_name) :: AttName
  character(len=80) :: pre_name,cr_date,cr_time,mem_name,ofile_name,coards,micro,v5dfields
!  integer :: nx,ny,nz,ng,nxend,nyend,nzend,nmem,size_ens,year,month,day,hour,minute,second,t_ind

  character(len=80) :: atype,l_name,units,field,positions
  integer :: indx,posdef,istag,jstag,kstag,dtype,ptype,btype,bdim
  

!-----------------------------------------
! variables vars

  character(nf90_max_name) :: VarName 
  character(nf90_max_name), allocatable :: VarNameArr(:)
  integer :: VarID, VarID1, nDimsVar, nAttsVar, VarType, varid2
  integer, allocatable :: varid2arr(:)
  integer, dimension(4) :: VarDim
  integer, dimension(nf90_max_var_dims) :: VarDimIDs
  integer, dimension(:), allocatable :: VarIDs, Var3Ds
  integer, allocatable :: VarStag(:,:), VarDimArr(:,:)
 
  character(nf90_max_name) :: TimeDimName
  integer :: TimeDimLength, TimeDimID

  integer :: xsize,ysize,zsize
  integer, dimension(:), allocatable :: mx,mxp1,my,myp1,mxi,mxp1i,myj,myp1j
  integer :: mxt,mxp1t,myt,myp1t,nt,mzt,mzp1t,dim2,dim3,mxadd,myadd,ntapp
  integer :: istart,istop,jstart,jstop
  integer :: length

  real    :: value
!  real, dimension(:), allocatable :: var1D, time, dt
  real, dimension(:), allocatable :: var1D, time, timeapp
  real, dimension(:,:), allocatable :: var2D
  real, dimension(:,:,:), allocatable :: var3D
  real, dimension(40,40,40) :: vardum

  integer istat
  integer memusage,memusage2,memusage3
  
  usehdf5 = .true.
  writehdf5 = .false. ! NOT WORKING, LEAVE FALSE; set true to use HDF5 to write vars to combined file.
!  usehdf5 = .false.
  one_file_at_a_time = .false. !  .true.  ! false for one var at a time; true to read everything from one file, then next file...
  one_time_at_a_time = .true. ! only for one_file_at_a_time.  Does not work well...
  one_time =  .true. ! .false.  ! only for nc4; if .true. then use netcdf4 to read one grid at a time, one variable at a time.
  openclose = .false.
!  openclose = .true.
  writeout = .true.
  
  IF ( usehdf5 ) THEN
    one_time = .false.
  ENDIF
  
  shuffle = 1
  deflate = 1
  deflate_level = 2


!-----------------------------------------------------------------------
! get namelist info

  IF ( COMMAND_ARGUMENT_COUNT() .lt. 1 ) THEN 
    write(0,*) 'Need the namelist filename in command line.'
    STOP
  ELSE

    CALL GET_COMMAND_ARGUMENT(1,run_file,length,status)
  
    IF ( status .ne. 0 ) THEN
      write(0,*) 'COULD NOT RETRIEVE 1ST COMMAND LINE ARG: EXITING RUN!...'
      write(0,*) 'DOES NOT HAVE A 3D.RUN FILE TO READ!'
      call exit(1)
    ENDIF
  ENDIF

!-----------------------------------------------------------------------
! Read in 3d.run file...

  INQUIRE(file=run_file(1:length),exist=file_exist)
  IF ( .not. file_exist ) THEN
    write(0,*) 'INPUT FILE:', run_file(1:length), ' DOES NOT EXIST!!!'
    write(0,*) 'EXITING'
    STOP
  ENDIF

! READ IN RUN MODEL NAMELIST

  IF ( .not. READ_NAMELIST(run_file(1:length),'RUN') ) print *,'PROBLEM READING NAMELIST'

! READ IN RUN ATTRIBUTES

  IF ( .not. READ_NAMELIST(run_file(1:length),'ATTRIBUTE') ) print *,'PROBLEM READING NAMELIST'

! READ IN GRID PARAMETERS

  IF ( .not. READ_NAMELIST(run_file(1:length),'GRIDN') ) print *,'PROBLEM READING NAMELIST'

! READ IN MPI PARAMETERS TO GET TILE SIZES

  IF ( .not. READ_NAMELIST(run_file(1:length),'MPI_PARAMS') ) print *,'PROBLEM READING NAMELIST'

  itile = nxt
  jtile = nyt

  ntilef_x = nx/itile    ! # of tiles in x-direction
  ntilef_y = ny/jtile    ! # of tiles in y-direction
  begtime  = start
  endtime  = stop
  
  open(15,file=run_file(1:length),status='old',form='formatted')
  rewind(15)
  read(15,NML=stitcher,iostat=istat)
  close(15)
      IF ( istat .ne. 0 ) THEN
        write(0,*) 'READ_NAMELIST: stitcher namelist not found -- using default values'
        write(0,*) 'Will therefore stitch ALL tiles and times=(start,stop) from RUN namelist.'
      ENDIF

  IF (begtime .lt. 0) begtime = 0
  IF (ntiles_x .lt. 0) ntiles_x = ntilef_x  ! if = -1, use the full domain
  IF (ntiles_y .lt. 0) ntiles_y = ntilef_y  ! if = -1, use the full domain

! READ IN STITCHER PARAMETERS

!  IF ( .not. READ_NAMELIST(run_file(1:length),'STITCHER') ) print *,'PROBLEM READING NAMELIST'

!-----------------------------------------------------------------------
! For now, will simply read the format of the tiles and keep that for the output file
!
! linking with netcdf4

  IF ( netcdfversion .eq. 4) THEN ! turn on hdf5 format to allow compression, but keep classic model
    IF ( deflate .ne. 0 .and. deflate_level .ge. 1 ) THEN ! use netcdf4
     cmode = ior(NF90_CLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
    ELSE ! if no compression, then use the classic model
     cmode = ior(NF90_64BIT_OFFSET, ior(NF90_CLOBBER,NF90_CLASSIC_MODEL) )
    ENDIF
  ENDIF   

!-----------------------------------------------------------------------
! read in the netCDF file name

!  write(*,"(/,'Input tile number of lower-left corner tile to be processed:')")
!  read(*,*) nllc
!  nllc = 5
!  write(*,"(/,'Input number of tiles in x-direction to be processed :')")
!  read(*,*) ntiles_x
!  ntiles_x = 3
!  write(*,"(/,'Input number of tiles in y-direction to be processed :')")
!  read(*,*) ntiles_y
!  ntiles_y = 2
!  ntilef_x = 4
!  ntilef_y = 3

!  write(*,"(/,'Input number of netCDF files to be processed :')")
!  read(*,*) nfile
!  nfile = 12
!  write(*,"(/,'Input the beginning time (in seconds):')")
!  read(*,*) begtime
!  begtime = 0
!  write(*,"(/,'Input the end time (in seconds):')")
!  read(*,*) endtime
!  endtime = 60

!  write(*,"(/,'Input output file prefix (prefix) :')")
!  read(*,*) prefix
!  prefix = 'nco_test_zieg'

!  write(*,"(/,'Input output file ensemble number(member) :')")
!  read(*,*) member
!  member = '000'

  nfile = ntilef_x * ntilef_y
  IF ( ibeg .le. 0 .and. jbeg .le. 0 ) THEN
    nllc = Max(1,nllc)
    ibeg = mod(nllc,ntilef_x)
    jbeg = (nllc-ibeg)/ntiles_y+1
  ELSEIF ( ibeg .ge. 1 .or. jbeg .ge. 1 ) THEN
    ibeg = Max(1,ibeg)
    jbeg = Max(1,jbeg)
    nllc = (jbeg-1)*ntilef_x + ibeg
  ENDIF
  
  nllc = Max(1,nllc)

  write(0,*) 'ibeg,jbeg,nllc = ',ibeg,jbeg,nllc

  print *, 'NFILE    = ',nfile
  print *, 'NPREFIX  = ',prefix
  print *, 'NENS     = ',member
  print *, 'BEG_TIME = ',begtime
  print *, 'END_TIME = ',endtime

  allocate(mx(nfile))
  allocate(mxp1(nfile))
  allocate(my(nfile))
  allocate(myp1(nfile))
  allocate(ncid(nfile))
  allocate(hdf5_fid(nfile))
  allocate(hdfid(nfile))
  allocate(file_used(nfile))
  allocate(fname(nfile))
  
  allocate(mxi(ntilef_x))
  allocate(mxp1i(ntilef_x))
  allocate(myj(ntilef_y))
  allocate(myp1j(ntilef_y))
  
  ncid(:) = -1
  hdf5_fid(:) = -100
  
  file_used(:) = 1

   CALL nf_set_chunk_cache(10,10,50);

  IF ( usehdf5 ) THEN ! old loop

! set file names

    DO i = 1,nfile
    
    write(tnum,"(i3.3)") i-1
!    fname(i) = trim(prefix)//'.'//trim(nens)//'.'//tnum//'.nc'
!    write(tnum,"(i3.3)") jj-2
    write(smember,"(i3.3)") member
    fname(i) = trim(prefix)//'.'//smember//'.'//tnum//'.nc'
    
    ENDDO

  
! first get the x-stuff
  DO i = 1,ntilef_x

    jj = i
    
    write(6,"(/,'jj,Reading ',i4,', ',A)") jj,fname(jj)

! open the netCDF file to read, added in call to disp_err to write out error message always

    CALL DISP_ERR( NF90_OPEN(fname(jj),nf90_nowrite,ncid(jj)), .true.)

! get dimensions here

    status = NF90_INQ_DIMID(ncid(i),"TIME",DimIDT)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDT,len=nt)
    status = NF90_INQ_DIMID(ncid(jj),"XC",DimIDXC)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDXC,len=mxi(i))
    status = NF90_INQ_DIMID(ncid(jj),"XE",DimIDXE)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDXE,len=mxp1i(i))

    IF ( i == 1 ) THEN
    status = NF90_INQ_DIMID(ncid(jj),"ZC",DimIDZC)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDZC,len=mzt)
    status = NF90_INQ_DIMID(ncid(jj),"ZE",DimIDZE)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDZE,len=mzp1t)
    ENDIF

    write(*,97) mxi(i),mxp1i(i),mzt,mzp1t,nt

    IF ( jj /= 1 ) THEN
      write(0,*) 'closing file jj; ncid = ',jj,ncid(jj)
      CALL DISP_ERR(nf90_close(ncid(jj)),.true.)
      ncid(jj) = -1
    ENDIF

    DO j = 1,ntilef_y
      jj = i + (j-1)*ntilef_x
      mx(jj) = mxi(i)
      mxp1(jj) = mxp1i(i)
      write(*,*) 'setting mx for file number ',jj,mx(jj)
    ENDDO


  ENDDO

! then get the y-stuff
  DO i = 1,ntilef_y

    jj =  1+(i-1)*ntilef_x
    write(6,"(/,'jj, Reading ',i4,', ',A)") jj,fname(jj)

! open the netCDF file to read, added in call to disp_err to write out error message always

    IF ( ncid(jj) .lt. 0 ) THEN
    CALL DISP_ERR( NF90_OPEN(fname(jj),nf90_nowrite,ncid(jj)), .true.)
    ENDIF

! get dimensions here

    status = NF90_INQ_DIMID(ncid(jj),"YC",DimIDYC)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDYC,len=myj(i))
    status = NF90_INQ_DIMID(ncid(jj),"YE",DimIDYE)
    status = NF90_INQUIRE_DIMENSION(ncid(jj),DimIDYE,len=myp1j(i))

    write(*,98) myj(i),myp1j(i)

    IF ( jj /= 1 ) THEN
      CALL DISP_ERR(nf90_close(ncid(jj)),.true.)
      ncid(jj) = -1
    ENDIF

    DO j = 1,ntilef_x
      jj = j + (i-1)*ntilef_x
      my(jj) = myj(i)
      myp1(jj) = myp1j(i)
      write(*,*) 'setting my for file number ',jj,my(jj)
    ENDDO


  ENDDO


  ELSEIF ( .true. ) THEN ! old loop
  DO i = 1,nfile

    write(tnum,"(i3.3)") i-1
!    fname(i) = trim(prefix)//'.'//trim(nens)//'.'//tnum//'.nc'
!    write(tnum,"(i3.3)") jj-2
    write(smember,"(i3.3)") member
    fname(i) = trim(prefix)//'.'//smember//'.'//tnum//'.nc'
    write(6,"(/,'Reading ',A)") fname(i)

! open the netCDF file to read, added in call to disp_err to write out error message always

    CALL DISP_ERR( NF90_OPEN(fname(i),nf90_nowrite,ncid(i)), .true.)

! get dimensions here

    status = NF90_INQ_DIMID(ncid(i),"TIME",DimIDT)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDT,len=nt)
    status = NF90_INQ_DIMID(ncid(i),"XC",DimIDXC)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDXC,len=mx(i))
    status = NF90_INQ_DIMID(ncid(i),"XE",DimIDXE)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDXE,len=mxp1(i))
    status = NF90_INQ_DIMID(ncid(i),"YC",DimIDYC)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDYC,len=my(i))
    status = NF90_INQ_DIMID(ncid(i),"YE",DimIDYE)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDYE,len=myp1(i))
    
    IF ( i == 1 ) THEN
    mxi(1) = mx(1)
    myj(1) = my(1)
    status = NF90_INQ_DIMID(ncid(i),"ZC",DimIDZC)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDZC,len=mzt)
    status = NF90_INQ_DIMID(ncid(i),"ZE",DimIDZE)
    status = NF90_INQUIRE_DIMENSION(ncid(i),DimIDZE,len=mzp1t)
    ENDIF

    write(*,96) mx(i),mxp1(i),my(i),myp1(i),mzt,mzp1t,nt
  ENDDO
  
  ELSE ! new test loop
  
  ii = 0
  DO j = nllc,nllc+(ntiles_y-1)*ntilef_x,ntilef_x
  DO i = 1,ntiles_x

    ii = ii+1
    jj = j+i
    write(tnum,"(i3.3)") jj-2
    write(smember,"(i3.3)") member
    fname(ii) = trim(prefix)//'.'//smember//'.'//tnum//'.nc'
    write(6,"(/,'Reading ',A)") fname(ii)

! open the netCDF file to read, added in call to disp_err to write out error message always

    CALL DISP_ERR( NF90_OPEN(fname(ii),nf90_nowrite,ncid(ii)), .true.)

! get dimensions here

    status = NF90_INQ_DIMID(ncid(ii),"TIME",DimIDT)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDT,len=nt)
    status = NF90_INQ_DIMID(ncid(ii),"XC",DimIDXC)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDXC,len=mx(ii))
    status = NF90_INQ_DIMID(ncid(ii),"XE",DimIDXE)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDXE,len=mxp1(ii))
    status = NF90_INQ_DIMID(ncid(ii),"YC",DimIDYC)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDYC,len=my(ii))
    status = NF90_INQ_DIMID(ncid(ii),"YE",DimIDYE)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDYE,len=myp1(ii))
    status = NF90_INQ_DIMID(ncid(ii),"ZC",DimIDZC)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDZC,len=mzt)
    status = NF90_INQ_DIMID(ncid(ii),"ZE",DimIDZE)
    status = NF90_INQUIRE_DIMENSION(ncid(ii),DimIDZE,len=mzp1t)

    write(*,96) mx(ii),mxp1(ii),my(ii),myp1(ii),mzt,mzp1t,nt
  ENDDO
  ENDDO
  
  ENDIF

96  format('nx = ',I5,5x,'nxp1 = ',I5,/,'ny = ',I5,5x,'nyp1 = ', I5,/,'nz = ',I5,5x,'nzp1 = ',I5,/,'nt = ',I5)

97  format('nx = ',I5,5x,'nxp1 = ',I5,/,'nz = ',I5,5x,'nzp1 = ',I5,/,'nt = ',I5)

98  format('ny = ',I5,5x,'nyp1 = ', I5)

!-----------------------------------------------------------------------
! get file info

  write(0,*) 'getting info from file 1; ncid = ',ncid(1)

  status = NF90_INQUIRE(ncid(1),nDims,nVars,nGlobalAtts)

! get global attributes

  allocate(GlAttName(nGlobalAtts))
  allocate(GlAttVali(nGlobalAtts))
  allocate(GlAttValc(nGlobalAtts))
  allocate(gtype(nGlobalAtts))

  DO n = 1,nGlobalAtts
    status = NF90_INQ_ATTNAME(ncid(1),NF90_GLOBAL,n,GlAttName(n))
    status = NF90_INQUIRE_ATTRIBUTE(ncid(1),NF90_GLOBAL,GlAttName(n),gtype(n))
    IF (gtype(n) == 2) THEN    ! type = 2 is character
      status = NF90_GET_ATT(ncid(1),NF90_GLOBAL,GlAttName(n),GlAttValc(n))
    ELSEIF (gtype(n) == 4) THEN    ! type = 4 is integer
      status = NF90_GET_ATT(ncid(1),NF90_GLOBAL,GlAttName(n),GlAttVali(n))
    ELSE
      write(0,"('Attribute ',a,' not a recogized type (i.e., INT or CHAR).')"),trim(GlAttName(n))
      STOP
    ENDIF
  ENDDO


!-----------------------------------------------------------------------
! get times

! times
   allocate(time(nt))
   time(:)=0.0
   status = NF90_INQ_VARID(ncid(1), 'TIME', varid)
   status = NF90_GET_VAR(ncid(1), varid,time(1:nt))
   IF (begtime > endtime) THEN
     write(0,"(/,'Beginning time must be less than or equal to end time. STOP!')")
     STOP
   ENDIF
! check that input times exist in the file 
   ntbeg = -1
   DO n = 1,nt
     IF (begtime == time(n)) ntbeg = n
   ENDDO
   IF (ntbeg == -1) THEN
     write(0,"(/,'Beginning time does not exist in ',a,': STOP!')") trim(fname(1))
     STOP    
   ENDIF
   ntend = -1
   DO n = 1,nt
     IF (endtime == time(n)) ntend = n
   ENDDO
   IF (ntend == -1) THEN
     write(0,"(/,'End time does not exist in ',a,': STOP!')") trim(fname(1))
     STOP
   ENDIF

   write(6,"(/,'Getting ',i2,' times.')") ntend-ntbeg+1
   write(6,"('ntbeg = ',i3,',   ntend = ',i3)") ntbeg, ntend
!   write(6,"(/,'nt / time(nt)', 2x,i3,2x,f8.1)") nt,time(nt)

!-----------------------------------------------------------------------
! Calculate length of x- and y-dimensions for stitched grid

    mxt = 0
    mxp1t = 0
    DO i = ibeg,ntiles_x+ibeg-1
      mxt = mxt + mx(i)
      mxp1t = mxp1t + mxp1(i)
      write(6,*) 'i,mxt,mx(i) = ',i,mxt,mx(i)
    ENDDO
    myt = 0
    myp1t = 0
    DO j = jbeg,ntiles_y+jbeg-1 ! ,ntiles_y*ntiles_x,ntilef_x
      jj = jbeg + (j-1)*ntilef_x
      myt = myt + my(jj)
      myp1t = myp1t + myp1(jj)
      write(6,*) 'jj,myt,my(jj) = ',jj,myt,my(jj)
    ENDDO

    write(6,"(/,'mxt   = ',i,/,'mxp1t = ',i,/,'myt   = ',i,/,'myp1t = ',i)")mxt,mxp1t,myt,myp1t
    write(6,*) 'ntiles_x, ntiles_y = ', ntiles_x, ntiles_y
    write(6,*) 'jbeg,ntiles_y*ntiles_x,ntilef_x = ',jbeg,ntiles_y*ntiles_x,ntilef_x

!-----------------------------------------------------------------------
! Write output, necessary variables

   
   IF ( prefix_out .eq. ' ' ) THEN
   write(*,"(/,'Input prefix of the new netCDF file (name):')") 
!   read(*,*) fname1
!   dump = 'stitcher_test_zieg'
   write(dum1,"(i6.6)") begtime
   write(dum2,"(i6.6)") endtime
   fname1 = trim(prefix)//'_'//dum1//'_'//dum2//'.000.nc'
   ELSE
     fname1 = trim(prefix_out)//'.nc'

   ENDIF
   print *,'output file is ',trim(fname1)
    
     INQUIRE(file=fname1, exist=file_exist)
     
     IF ( file_exist ) THEN
       append = .true.
       write(6,*) 'named output file exists, assuming we are appending...'
     ELSE
       append = .false.
       write(6,*) 'named output file does NOT exist, assuming we are creating...'
     ENDIF
    
    IF ( .not. append ) THEN
    
    ntstart = 1
#ifdef NC4
! linking with netcdf4
   write(0, *) 'netcdfversion = ',netcdfversion
   IF ( netcdfversion .eq. 4 ) THEN
     cmode = ior(NF90_CLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
   ELSE
     cmode = NF90_CLOBBER
   ENDIF
!   status = NF90_CREATE(fname1,cmode,ncid1)
   call disp_err( NF90_CREATE(fname1,cmode,ncid1), .true. )


#else
! linking with netcdf3:
   call disp_err( NF90_CREATE(fname1,nf90_clobber,ncid1), .true. )
#endif
   
   ELSE
   
     INQUIRE(file=fname1, exist=file_exist)
     
     IF ( file_exist ) THEN
        write(6,*) 'Opening file for appending: ',fname1

        CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)

        status = NF90_INQ_DIMID(ncid1,"TIME",DimIDTapp)
        status = NF90_INQUIRE_DIMENSION(ncid1,DimIDTapp,len=ntapp)

        allocate(timeapp(ntapp))
 
        status = NF90_INQ_VARID(ncid1, 'TIME', varid)
        status = NF90_GET_VAR(ncid1, varid,timeapp(1:ntapp))
        
        ntstart = 1

        DO n = 1,ntapp
         IF (begtime > timeapp(n)) ntstart = n + 1
         IF (begtime == timeapp(n)) ntstart = n ! overwrite this time
        ENDDO
        
        write(6,*) 'ntstart = ', ntstart
         
     
     ELSE
       write(0,*) 'Output file not found! ', fname1
       write(0,*) 'Halting program.'
       STOP
     ENDIF
   
   ENDIF

   IF ( .not. append ) THEN
   status = nf90_def_dim(ncid1,"TIME",nf90_UNLIMITED,DimIDT)
   status = nf90_def_dim(ncid1,"XC",mxt,DimIDXC)
   status = nf90_def_dim(ncid1,"XE",mxp1t,DimIDXE)
   status = nf90_def_dim(ncid1,"YC",myt,DimIDYC)
   status = nf90_def_dim(ncid1,"YE",myp1t,DimIDYE)
   status = nf90_def_dim(ncid1,"ZC",mzt,DimIDZC)
   status = nf90_def_dim(ncid1,"ZE",mzp1t,DimIDZE)
!   status = nf90_def_dim(ncid1,"dim2",dim2,DimIDdim2)
!   status = nf90_def_dim(ncid1,"dim3",dim3,DimIDdim3)

   write(*,"(/,'Writing global attributes')")

   DO n = 1,nGlobalAtts
     IF (GlAttName(n) == 'NX') GlAttVali(n) = mxp1t
     IF (GlAttName(n) == 'NY') GlAttVali(n) = myp1t
     IF (GlAttName(n) == 'NZ') GlAttVali(n) = mzp1t
     IF (GlAttName(n) == 'NXEND') GlAttVali(n) = mxp1t
     IF (GlAttName(n) == 'NYEND') GlAttVali(n) = myp1t
     IF (GlAttName(n) == 'NZEND') GlAttVali(n) = mzp1t
     IF (GlAttName(n) == 'MEMBER') GlAttVali(n) = 0
     IF (GlAttName(n) == 'OUTPUT_FILE_NAME') GlAttValc(n) = fname1
     IF (GlAttName(n) == 'TILE_INDEX') GlAttVali(n) = 0
     IF (gtype(n) == 2) THEN
       status = NF90_PUT_ATT(ncid1,NF90_GLOBAL,GlAttName(n),GlAttValc(n))
     ELSE
       status = NF90_PUT_ATT(ncid1,NF90_GLOBAL,GlAttName(n),GlAttVali(n))
     ENDIF
   ENDDO
   
   ENDIF

!--------------------------------------------------------------------------
! Get and write variables

!   allocate(VarName(nVars))
   allocate(VarIDs(nVars))
   allocate(Var3Ds(nVars))
   allocate(VarStag(3,nVars))
   allocate(varid2arr(nVars))
   allocate( VarNameArr(nVars) )
   Var3Ds(:) = 0
   cnt3d = 0
   status = nf90_inquire(ncid(1), unlimitedDimId = DimIDT) ! find dimid for TIME
   
   IF ( .not. one_time ) THEN
   
   
   IF ( one_file_at_a_time ) tim1a = dtime(timarr)
   
   DO vid = 1,nVars
     
     VarDim(:) = 0
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,VarName,ndims=nDimsVar,dimids=VarDim)
!     print *,trim(VarName), nDimsVar,' VarDim =',VarDim
     write(6,*) 'vid loop: nDimsVar = ',nDimsVar,' VarName = ',trim(VarName)
     Var3Ds(vid) = nDimsVar
     
     IF (nDimsVar == 0 .and. .not. append) CALL READ_WRITE_SCALAR(vid)
     IF (nDimsVar == 1) THEN
        IF (VarDim(1) == 1) THEN
          CALL READ_WRITE_SCALAR(vid)
        ELSE
          IF ( .not. append ) THEN
           CALL READ_WRITE_VAR1D(vid)
          ELSE
            IF ( VarDim(1) == DimIDT ) CALL READ_WRITE_VAR1D(vid)
          ENDIF
        ENDIF
     ENDIF
     IF (nDimsVar == 2) THEN
        IF (VarDim(2) == 1) THEN ! time-dependent 1-d variable?
          CALL READ_WRITE_VAR1D(vid)
        ELSE
          IF ( usehdf5 ) THEN
            IF (VarName(3:4) /= 'DX') CALL DEF_VAR3D(vid)
            write(*,*) 'defined 2d variable'
! skip the XCDX, XEDX, YCDX, YEDX, ZCDX,and ZEDX variables
          ELSE
            IF (VarName(3:4) /= 'DX') CALL READ_WRITE_VAR2D(vid)
          ENDIF
        ENDIF
     ENDIF
     IF (nDimsVar == 3) THEN
!     write(6,*) 'nDimsVar = ',nDimsVar,' VarName = ',VarName
      IF ( one_file_at_a_time .or. (usehdf5 .and. VarDim(3) /= 1) ) THEN
!        Var3D(vid) = 4
        CALL DEF_VAR3D(vid)
      ELSE
        IF (VarDim(3) == 1) THEN
          IF ( usehdf5 ) THEN
            CALL DEF_VAR3D(vid)
            write(*,*) 'defined 2d variable ',vid
!            CALL READ_WRITE_VAR2D_hdf5(vid,VarName)
          ELSE
            CALL READ_WRITE_VAR2D(vid)
          ENDIF
        ELSE
          IF ( usehdf5 ) THEN
!            CALL READ_WRITE_VAR3D_hdf5(vid,VarName)
          ELSE
            CALL READ_WRITE_VAR3D(vid)
          ENDIF
        ENDIF
       ENDIF
!     call sleep(5)
     ENDIF
     IF (nDimsVar == 4) THEN
      cnt3d = cnt3d + 1
      
      IF ( one_file_at_a_time .or. usehdf5) THEN
!        Var3D(vid) = 4
          CALL get_mem_used(memusage)
          memusage3 = memusage
        CALL DEF_VAR3D(vid)
          CALL get_mem_used(memusage)
          write(0,*) 'DEF_VAR3D memory increase = ',memusage-memusage3
      ELSE
       tim1 = dtime(timarr) 
       timeread = 0.0
       timewrite = 0.0
       
        CALL READ_WRITE_VAR3D(vid)
!      IF ( VarName == 'V' )  CALL READ_WRITE_VAR3D(vid)
!       IF ( cnt3d .ge. 4 ) EXIT
       
!      CALL h5garbage_collect_f(status)
!      IF ( status .ne. 0 ) THEN
!        write(0,*) 'Error in garbage_collect, status = ',status
!      ENDIF
      
!      IF ( .false. .and. (cnt3d/10)*10 .eq. cnt3d ) THEN
       IF ( openclose ) THEN
       write(6,*) 'close/open files file/ID '
       DO i = 1,nfile ! MIN(10,nfile)

! close and open the netCDF file to read, added in call to disp_err to write out error message always
!        write(6,*) 'close file/ID ',i,ncid(i)
        IF ( ncid(i) .gt. 0 ) THEN
        CALL DISP_ERR(nf90_close(ncid(i)),.true.)
         ncid(i) = -1
!        CALL DISP_ERR( NF90_OPEN(fname(i),nf90_nowrite,ncid(i)), .true.)
        ENDIF

      ENDDO
       DO i = 1,nfile ! MIN(10,nfile)

! close and open the netCDF file to read, added in call to disp_err to write out error message always
!        write(6,*) 'open file/ID ',i,ncid(i)
!        CALL DISP_ERR(nf90_close(ncid(i)),.true.)
        CALL DISP_ERR( NF90_OPEN(fname(i),nf90_nowrite,ncid(i)), .true.)

       ENDDO
       ENDIF ! openclose
      
            tim2 = dtime(timarr)
      
            write(iunit,*) 'time for variable ',trim(VarName),' = ',tim2
            write(iunit,*) 'TOTALS TIME FOR read/write ',trim(VarName),' = ',timeread,timewrite

      ENDIF
     ENDIF  ! one_file_at_a_time
     
!     call sleep(2)
     
   ENDDO


   IF ( usehdf5 ) THEN
   
        CALL h5open_f(error) 

      IF ( writehdf5 ) THEN
       CALL DISP_ERR( nf90_close(ncid1) , .true. )
!(fname1,nf90_write,ncid1)      
         CALL h5fopen_f (fname1, H5F_ACC_RDWR_F, ncid1, error)
         IF ( error /= 0 ) THEN
           write(0,*) 'error opening hdf5 file, error = ',error
           STOP
         ENDIF

      ENDIF
   
   allocate( VarDimArr(4,nVars) )
   
   DO vid = 1,nVars
     
     VarDimArr(:,vid) = 0
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,VarName,ndims=nDimsVar,dimids=VarDim)
!     print *,trim(VarName), nDimsVar,' VarDim =',VarDim
!     write(6,*) 'vid loop: nDimsVar = ',nDimsVar,' VarName = ',trim(VarName)
     Var3Ds(vid) = nDimsVar
      VarDimArr(1:4,vid) = VarDim(1:4)
      VarNameArr(vid) = VarName

   ENDDO
        DO n = 1,nfile ! MIN(10,nfile)
! close and open the netCDF file to read, added in call to disp_err to write out error message always
        IF ( ncid(n) .gt. 0 ) THEN
        write(6,*) 'close file/ID ',n,ncid(n)
        CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        ncid(n) = -1
!         CALL h5fopen_f (fname(n), H5F_ACC_RDONLY_F, hdfid(n), error)
!         IF ( error /= 0 ) THEN
!           write(0,*) 'error opening hdf5 file, error = ',error
!           STOP
!         ENDIF
        ENDIF
        
       ENDDO

          CALL get_mem_used(memusage)
          memusage2 = memusage

       DO ktime = ntbeg,ntend
!       kk = k-tbeg+1
       
       timeread = 0.0
       timewrite = 0.0
       tim1 = dtime(timarr) 

!       CALL DISP_ERR( nf90_close(ncid1) , .true. )
!       CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)

        DO vid = 1,nVars
         nDimsVar = Var3Ds(vid)
         VarDim(1:4) = VarDimArr(1:4,vid)
        
          IF ( Var3Ds(vid) .eq. 4 ) THEN
            CALL  READ_WRITE_VAR3D_hdf5(vid,VarNameArr(vid))
!            CALL FILEREAD_WRITE_3D_one_time(vid,n,ktime,ktime-ntbeg+ntstart)
!            call sleep(0.2)
!          call sleep(1)
          ELSEIF ( Var3Ds(vid) .ge. 2  .and. VarNameArr(vid)(3:4) /= 'DX') THEN
                CALL READ_WRITE_VAR2D_hdf5(vid,VarNameArr(vid))
          ENDIF

            CALL get_mem_used(memusage)
            
            write(0,*) 'H5 memory usage read/write, delta =',memusage, memusage-memusage2
            
            memusage2 = memusage
        
        ENDDO

            tim2 = dtime(timarr)
      
            write(iunit,*) 'TIME FOR TIME ',ktime,' = ',tim2
            write(iunit,*) 'TIME FOR read/write ',ktime,' = ', timeread, timewrite
            write(iunit,*) 'TIME FOR SUM,DIFF: ', timeread+timewrite, tim2-(timeread + timewrite)
          
         timeread = 0.0
         timewrite = 0.0
         
       ENDDO ! ktime

      IF ( writehdf5 ) THEN
         CALL h5fclose_f (ncid1, error)
         IF ( error /= 0 ) THEN
           write(0,*) 'error closing hdf5 file, error = ',error
           STOP
         ENDIF

      ENDIF

   
   GOTO 99
   
   ENDIF

!   IF ( append ) STOP
       
   IF ( one_file_at_a_time ) THEN

   allocate( VarDimArr(4,nVars) )
   
   DO vid = 1,nVars
     
     VarDimArr(:,vid) = 0
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,VarName,ndims=nDimsVar,dimids=VarDim)
!     print *,trim(VarName), nDimsVar,' VarDim =',VarDim
!     write(6,*) 'vid loop: nDimsVar = ',nDimsVar,' VarName = ',trim(VarName)
     Var3Ds(vid) = nDimsVar
      VarDimArr(1:4,vid) = VarDim(1:4)
      VarNameArr(vid) = VarName

   ENDDO

     tim2a = dtime(timarr)
     write(6,*) 'TIME FOR INITIAL WRITING = ',tim2a

       
       IF ( openclose ) THEN

        DO n = 1,nfile ! MIN(10,nfile)


! close and open the netCDF file to read, added in call to disp_err to write out error message always
        write(6,*) 'close file/ID ',n,ncid(n)
        CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        ncid(n) = -1
!        CALL DISP_ERR( NF90_OPEN(fname(n),nf90_nowrite,ncid(n)), .true.)
       ENDDO
       
       ENDIF


       IF ( one_time_at_a_time ) THEN
       
       write(0,*) 'one_time_at_a_time'

           CALL DISP_ERR( nf90_close(ncid1) , .true. )
!         cmode = ior(NF90_NOCLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
!   status = NF90_CREATE(fname1,cmode,ncid1)
!   call disp_err( NF90_CREATE(fname1,cmode,ncid1), .true. )
          CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)


       DO ktime = ntbeg,ntend
!       kk = k-tbeg+1
       write(0,*) 'one_time_at_a_time, ktime = ',ktime
       
       timeread = 0.0
       timewrite = 0.0
       tim1 = dtime(timarr) 

       DO n = 1,nfile ! MIN(10,nfile)
          
         IF ( file_used(n) == 0 ) CYCLE
         
! close and open the netCDF file to read, added in call to disp_err to write out error message always
!        CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        IF ( openclose ) THEN
          CALL DISP_ERR( NF90_OPEN(fname(n),nf90_nowrite,ncid(n)), .true.)
          write(0,*) 'open file/ID ',n,ncid(n),fname(n)
        ENDIF

        DO vid = 1,nVars
         nDimsVar = Var3Ds(vid)
         VarDim(1:4) = VarDimArr(1:4,vid)
        
!         write(0,*) 'ndimsvar = ',nDimsVar
         IF (nDimsVar == 3) THEN
!       write(6,*) 'nDimsVar = ',nDimsVar,' VarName = ',VarName
          IF (VarDim(3) == 1) THEN
!            write(0,*) 'call READ_WRITE_VAR2D_onetime',vid,n,ktime,ktime-ntbeg+ntstart
            CALL READ_WRITE_VAR2D_onetime(vid,n,ktime,ktime-ntbeg+ntstart)
!          call sleep(1)
          ELSE
             write(0,*) 'what variable is this???'
             STOP
!          CALL READ_WRITE_VAR3D(vid)
          ENDIF
!     call sleep(5)
         ENDIF
          IF ( Var3Ds(vid) .eq. 4 ) THEN
            CALL FILEREAD_WRITE_3D_one_time(vid,n,ktime,ktime-ntbeg+ntstart)
!            call sleep(0.2)
!          call sleep(1)
          ENDIF
        
        
        ENDDO

        IF ( openclose ) THEN
          write(0,*) 'close file/ID ',n,ncid(n)
!         call sleep(5)
          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        ENDIF


       ENDDO ! n
       
!   close/open the output netcdf file

           IF ( writeout ) THEN
           
           CALL DISP_ERR( nf90_close(ncid1) , .true. )
!         cmode = ior(NF90_NOCLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
!   status = NF90_CREATE(fname1,cmode,ncid1)
!   call disp_err( NF90_CREATE(fname1,cmode,ncid1), .true. )
          CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)
          
          ENDIF

            tim2 = dtime(timarr)
      
            write(iunit,*) 'TIME FOR TIME ',ktime,' = ',tim2
            write(iunit,*) 'TIME FOR read/write ',ktime,' = ', timeread, timewrite
            write(iunit,*) 'TIME FOR SUM,DIFF: ', timeread+timewrite, tim2-(timeread + timewrite)
          
         timeread = 0.0
         timewrite = 0.0
         
       ENDDO ! ktime
       
       
       ELSE

       DO n = 1,nfile ! MIN(10,nfile)

! close and open the netCDF file to read, added in call to disp_err to write out error message always
        IF ( openclose ) THEN
          write(0,*) 'open file/ID ',n,ncid(n)
!          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
          CALL DISP_ERR( NF90_OPEN(fname(n),nf90_nowrite,ncid(n)), .true.)
        ENDIF

        DO vid = 1,nVars
          IF ( Var3Ds(vid) .eq. 4 ) THEN
            CALL FILEREAD_WRITE_3D(vid,n)
          ENDIF
        
        ENDDO

        IF ( openclose ) THEN
          write(0,*) 'close file/ID ',n,ncid(n)
!          call sleep(5)
        
          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        ENDIF

       ENDDO ! n
       
       ENDIF
   
   
   ENDIF
   
   
   ELSE ! one_time

   allocate( VarDimArr(4,nVars) )
   
   DO vid = 1,nVars
     
     VarDimArr(:,vid) = 0
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,VarName,ndims=nDimsVar,dimids=VarDim)
!     print *,trim(VarName), nDimsVar,' VarDim =',VarDim
!     write(6,*) 'vid loop: nDimsVar = ',nDimsVar,' VarName = ',trim(VarName)
     Var3Ds(vid) = nDimsVar
      VarDimArr(1:4,vid) = VarDim(1:4)
      VarNameArr(vid) = VarName

   ENDDO
   
          CALL get_mem_used(memusage)
          memusage2 = memusage
          memusage3 = memusage

   DO ktime = ntbeg,ntend
       tim1 = dtime(timarr) 

   DO vid = 1,nVars
     
!     print *,trim(VarName), nDimsVar,' VarDim =',VarDim
!     write(6,*) 'vid loop: nDimsVar = ',nDimsVar,' VarName = ',trim(VarName)
     nDimsVar = Var3Ds(vid)
     VarDim(1:4) = VarDimArr(1:4,vid)
     VarName = VarNameArr(vid)
     
     IF (nDimsVar == 0 .and. .not. append) CALL READ_WRITE_SCALAR(vid)
     IF (nDimsVar == 1) THEN
        IF (VarDim(1) == 1) THEN
          CALL READ_WRITE_SCALAR(vid)
        ELSE
          IF ( .not. append .and. ktime .eq. ntbeg ) THEN
           CALL READ_WRITE_VAR1D(vid)
          ELSE
            IF ( VarDim(1) == DimIDT ) CALL READ_WRITE_VAR1D(vid)
          ENDIF
        ENDIF
     ENDIF
     IF (nDimsVar == 2) THEN
        IF (VarDim(2) == 1) THEN
          CALL READ_WRITE_VAR1D(vid)
        ELSE
! skip the XCDX, XEDX, YCDX, YEDX, ZCDX,and ZEDX variables
          IF (VarName(3:4) /= 'DX') CALL READ_WRITE_VAR2D(vid)
        ENDIF
     ENDIF
     IF (nDimsVar == 3) THEN
!     write(6,*) 'nDimsVar = ',nDimsVar,' VarName = ',VarName
        IF (VarDim(3) == 1) THEN
          CALL READ_WRITE_VAR2D(vid)
        ELSE
          CALL READ_WRITE_VAR3D(vid)
        ENDIF
!     call sleep(5)
     ENDIF
     IF (nDimsVar == 4) THEN
      cnt3d = cnt3d + 1
      
!        Var3D(vid) = 4
        IF ( ktime == ntbeg ) THEN
          CALL DEF_VAR3D(vid)
        ENDIF


       IF ( one_file_at_a_time ) THEN
       DO n = 1,nfile ! MIN(10,nfile)
          
         IF ( file_used(n) == 0 ) CYCLE
         
! close and open the netCDF file to read, added in call to disp_err to write out error message always
!        CALL DISP_ERR(nf90_close(ncid(n)),.true.)

            CALL FILEREAD_WRITE_3D_one_time(vid,n,ktime,ktime-ntbeg+ntstart)


       ENDDO ! n
       
       ELSE
       
            CALL get_mem_used(memusage)
            memusage3 = memusage

            CALL READ_WRITE_VAR3D(vid, 1)
            
            CALL get_mem_used(memusage)
            
            memusage2 = memusage

        IF ( openclose ) THEN
         DO n = 1,nfile ! MIN(10,nfile)
          
          IF ( file_used(n) == 0 ) CYCLE
          write(0,*) 'close file/ID ',n,ncid(n)
!          call sleep(5)
        
          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        
          write(0,*) 'open file/ID ',n,ncid(n)
!          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
          CALL DISP_ERR( NF90_OPEN(fname(n),nf90_nowrite,ncid(n)), .true.)
          
          ENDDO
        ENDIF

            CALL get_mem_used(memusage)

            write(0,*) 'NC4 memory usage, delta =',memusage, memusage-memusage2, memusage2-memusage3
            

       ENDIF
       
!   close/open the output netcdf file

!          CALL DISP_ERR( nf90_close(ncid1) , .true. )
!         cmode = ior(NF90_NOCLOBBER,ior(NF90_CLASSIC_MODEL,NF90_NETCDF4))
!   status = NF90_CREATE(fname1,cmode,ncid1)
!   call disp_err( NF90_CREATE(fname1,cmode,ncid1), .true. )
!          CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)

          
        
       ENDIF 
     
!     call sleep(2)
     
   ENDDO ! vid

       DO n = 1,nfile ! MIN(10,nfile)
          
!         IF ( file_used(n) == 0 ) CYCLE
         
        IF ( openclose ) THEN
          write(0,*) 'close file/ID ',n,ncid(n)
!         call sleep(5)
          CALL DISP_ERR(nf90_close(ncid(n)),.true.)
        ENDIF

! close and open the netCDF file to read, added in call to disp_err to write out error message always
!        CALL DISP_ERR(nf90_close(ncid(n)),.true.)
       
       IF ( ktime .lt. ntend ) THEN
        IF ( openclose ) THEN
          CALL DISP_ERR( NF90_OPEN(fname(n),nf90_nowrite,ncid(n)), .true.)
          write(0,*) 'open file/ID ',n,ncid(n)
        ENDIF
       ENDIF



       ENDDO ! n
   

            tim2 = dtime(timarr)
      
            write(iunit,*) 'TIME FOR TIME ',ktime,' = ',tim2
   
     append = .true.
   ENDDO ! ktime


   
   ENDIF
     

   write(6,"(/,'Completed the 1D and 2D fields',/)")

!-----------------------------------------------------------------------
!  Read/Write 3D variables
! 
! LJW:  I created a subroutine to handle to do the read/write to create
!       code that was a bit more manageable

!   CALL READ_WRITE_VAR3D("DBZ",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("WZ",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("U",(/DimIDXE,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("V",(/DimIDXC,DimIDYE,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("W",(/DimIDXC,DimIDYC,DimIDZE/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("PI",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("KM",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("TH",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("QV",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("QC",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)
!   CALL READ_WRITE_VAR3D("QR",(/DimIDXC,DimIDYC,DimIDZC/),ntbeg,ntend)

!-----------------------------------------------------------------------
! close files
!-----------------------------------------------------------------------
   
   IF ( .not. one_file_at_a_time .and. .not. usehdf5 ) THEN
   
   DO i=1,nfile
     write(0,*) 'end of program, close file ',i
     CALL DISP_ERR(nf90_close(ncid(i)),.true.)
   ENDDO
   
   ENDIF

   99 CONTINUE
   
   CALL DISP_ERR( nf90_close(ncid1) , .true. )

   write(6,"(/,'Closed the files',/)")

     !
     ! Close FORTRAN interface.
     !
   IF ( usehdf5 )  CALL h5close_f(error)
   
   
   deallocate(mx)
   deallocate(mxp1)
   deallocate(my)
   deallocate(myp1)
   deallocate(ncid)
   deallocate(file_used)
   deallocate(fname)
   deallocate(VarIDs)
   deallocate(Var3Ds)
   deallocate(VarStag)
   deallocate(varid2arr)

  CONTAINS

   SUBROUTINE READ_WRITE_SCALAR(vid)

     implicit none
     character(len=40) :: vname, AttName, AttValC
     integer       :: AttValI 
     integer, dimension(4) :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, tbeg1

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)
     write(6,"(/,'Read/Write SCALAR: name = ',a)") trim(vname)
     IF ( .not. append ) THEN
       status = NF90_REDEF(ncid1)
       IF (numDims == 1) THEN
         status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
         status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1),varid2)
       ELSE
         status = NF90_DEF_VAR(ncid1,trim(vname),vtype,varid2)
       ENDIF

     ELSE
       status = nf90_inq_varid(ncid1, trim(vname), varid2)
       write(0,*) 'read_write_scalar: vname, varid2 : ',trim(vname),', ',varid2
     ENDIF

! grab and write variable attributes
     DO k = 1,numAtts

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
       ENDIF

     ENDDO
     IF ( .not. append ) status = NF90_ENDDEF(ncid1)

     IF (numDims == 1) THEN
       tbeg = ntbeg
       tend = ntend
       tbeg1 = ntbeg
       IF ( one_time ) THEN
         tbeg = ktime
         tend = ktime
       ENDIF
     ELSE
       tbeg = 1
       tend = 1
       tbeg1 = 1
     ENDIF
     
       
     DO k = tbeg,tend
       kk = k-tbeg1+ntstart
       status = NF90_GET_VAR(ncid(1),vid,value,(/k/))
       IF (vname == 'NX') THEN
         status = NF90_PUT_VAR(ncid1,varid2,mxp1t,(/1/)) 
       ELSE IF (vname == 'NY') THEN
         status = NF90_PUT_VAR(ncid1,varid2,myp1t,(/1/)) 
       ELSE IF (vname == 'NZ') THEN
         status = NF90_PUT_VAR(ncid1,varid2,mzp1t,(/1/)) 
       ELSEIF (vname == 'NXEND') THEN
         status = NF90_PUT_VAR(ncid1,varid2,mxp1t,(/1/)) 
       ELSE IF (vname == 'NYEND') THEN
         status = NF90_PUT_VAR(ncid1,varid2,myp1t,(/1/)) 
       ELSE IF (vname == 'NZEND') THEN
         status = NF90_PUT_VAR(ncid1,varid2,mzp1t,(/1/)) 
       ELSE  
         IF ( vname == 'TIME' ) THEN
           write(6,*) 'read_write_scalar: TIME: value,kk = ',value,kk
         ENDIF
         status = NF90_PUT_VAR(ncid1,varid2,value,(/kk/)) 
       ENDIF
     ENDDO
       
   RETURN
   END SUBROUTINE 

! ###################################################################################

   SUBROUTINE READ_WRITE_VAR1D(vid)

     implicit none

     integer       :: shape
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
     integer, dimension(4) :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, tbeg1
     real, allocatable :: var1d(:)

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

     write(6,"(/,'Read/Write VAR1D: name = ',a,/)") vname
     IF ( .not. append ) status = NF90_REDEF(ncid1)
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
     
     IF ( .not. append ) THEN
       IF (numDims == 1) THEN
         status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1),varid2)
       ELSE
         status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:2),varid2)
       ENDIF
     ELSE
       status = nf90_inq_varid(ncid1, trim(vname), varid2)
       write(0,*) 'read_write_var1d: vname, varid2 : ',trim(vname),', ',varid2
     ENDIF

! grab and write variable attributes
     DO k = 1,numAtts

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     IF ( .not. append ) status = NF90_ENDDEF(ncid1)

     IF (numDims == 2) THEN
       tbeg = ntbeg
       tend = ntend
       tbeg1 = ntbeg
       IF ( one_time ) THEN
         tbeg = ktime
         tend = ktime
       ENDIF
     ELSE
       tbeg = 1
       tend = 1
       tbeg1 = 1
     ENDIF


     DO k = tbeg,tend
       kk = k-tbeg1+ntstart
       istop=0
       IF (vdim == 'z') THEN
         IF (kstag == 1) zsize = mzp1t
         IF (kstag == 0) zsize = mzt
         istop = zsize
         allocate(var1D(zsize))
         var1D(:) = 0.0
         status = NF90_GET_VAR(ncid(1),vid,var1D,(/k/))
       ENDIF

       IF (vdim == 'x') THEN
         IF (istag == 1) xsize = mxp1t
         IF (istag == 0) xsize = mxt
         allocate(var1D(xsize))
         write(0,*) 'x, var1D: ysize = ',xsize
         var1D(:) = 0.0
!         DO i = 1,ntiles_x
!         DO i = nllc,nllc+ntiles_x-1
         DO i = ibeg,ibeg+ntiles_x-1
           IF (istag == 1) mxadd = mxp1(i)
           IF (istag == 0) mxadd = mx(i)
           istart = istop+1
           istop  = istart+mxadd-1
!           write(*,*) 'var1D: i = ',i,ncid(i)
           IF ( ncid(i) < 0 ) THEN
             CALL DISP_ERR( NF90_OPEN(fname(i),nf90_nowrite,ncid(i)), .true.)
           ENDIF
             CALL DISP_ERR(  NF90_GET_VAR(ncid(i),vid,var1D(istart:istop),(/k/)), .true. )
           IF ( i /= 1 ) THEN
             CALL DISP_ERR( NF90_CLOSE(ncid(i)), .true.)
             ncid(i) = -1
           ENDIF
         ENDDO
         write(0,*) 'x, var1D: istop = ',istop
         istop = xsize
       ENDIF

       IF (vdim == 'y') THEN
         IF (jstag == 1) ysize = myp1t
         IF (jstag == 0) ysize = myt
         allocate(var1D(ysize))
         write(0,*) 'y, var1D: ysize = ',ysize
         var1D(:) = 0.0
!         DO i = nllc,nllc+(ntiles_y-1)*ntilef_x,ntilef_x
         DO i = (jbeg-1)*ntilef_x+1, (jbeg + ntiles_y - 2)*ntilef_x+1 ,ntilef_x
           IF ( i .gt. nfile ) THEN
             write(0,*) 'warning! i > nfiles! ',i,nllc,ntiles_y,ntilef_x,nllc+(ntiles_y-1)*ntilef_x
             STOP
           ENDIF
           IF (jstag == 1) myadd = myp1(i)
           IF (jstag == 0) myadd = my(i)
           istart = istop+1
           istop  = istart+myadd-1
!           write(*,*) 'var1D: i = ',i,ncid(i)
           IF ( ncid(i) < 0 ) THEN
             CALL DISP_ERR( NF90_OPEN(fname(i),nf90_nowrite,ncid(i)), .true.)
           ENDIF
             CALL DISP_ERR( NF90_GET_VAR(ncid(i),vid,var1D(istart:istop),(/k/)), .true.)
           IF ( i /= 1 ) THEN
             CALL DISP_ERR( NF90_CLOSE(ncid(i)), .true.)
             ncid(i) = -1
           ENDIF
         ENDDO
         write(0,*) 'y, var1D: istop = ',istop
         istop = ysize
!         write(0,*) var1D
       ENDIF

       write(0,*) 'VAR1D: istop,kk = ',istop,kk
       IF ( kk == 1 ) THEN
         status = NF90_PUT_VAR(ncid1,varid2,var1D,(/1/)) 
       ELSE
         status = NF90_PUT_VAR(ncid1,varid2,var1D,(/1/),(/istop/),(/kk/)) 
       ENDIF
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     ENDDO
       IF ( allocated(var1D) ) deallocate(var1D)
         
   RETURN
   END SUBROUTINE

! ####################################################################################

   SUBROUTINE READ_WRITE_VAR2D(vid)

     implicit none

     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, tbeg1
!     integer       :: ibeg, jbeg
     real, allocatable :: tem2d(:,:)
     integer ii,jj,ijtile

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)
     IF ( .not. append ) status = NF90_REDEF(ncid1)
     write(6,"(/,'Read/Write VAR2D: name = ',a,/)") vname
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,xtype=vtype,ndims=numDims)
     allocate(DimIDV(numDims))
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
     
     IF ( .not. append ) THEN
       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)
     ELSE
       status = nf90_inq_varid(ncid1, trim(vname), varid2)
     ENDIF

!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 2) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:2),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ENDIF

! grab and write variable attributes
     DO k = 1,numAtts

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     IF ( .not. append ) status = NF90_ENDDEF(ncid1)

     IF (numDims == 3) THEN
       tbeg = ntbeg
       tend = ntend
       tbeg1 = ntbeg
       IF ( one_time ) THEN
         tbeg = ktime
         tend = ktime
       ENDIF
     ELSE
       tbeg = 1
       tend = 1
       tbeg1 = 1
     ENDIF


     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     allocate(var2D(xsize,ysize))

! Find the x,y indices of the starting tile 
!     ibeg = mod(nllc,ntilef_x)
!     jbeg = (nllc-ibeg)/ntiles_y+1
!     jstp = 0
!     DO j = 2,jbeg
!       jj = (j-2)*ntilef_x+1
!       IF (jstag == 1) myadd = myp1(jj)
!       IF (jstag == 0) myadd = my(jj)
!       jstp = jstp + myadd
!     ENDDO
!     istp = 0
!     DO i = 2,ibeg
!       IF (istag == 1) mxadd = mxp1(i)
!       IF (istag == 0) mxadd = mx(i)
!       istp = istp + mxadd
!     ENDDO
     
     DO k = tbeg,tend
       kk = k-tbeg1+ntstart
       var2D(:,:) = 0.0

       jstop=0
       DO j = jbeg,ntiles_y+jbeg-1
         istop=0
         jstart = jstop+1
         DO i = ibeg,ntiles_x+ibeg-1
           ijtile = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(ijtile)
           IF (istag == 0) mxadd = mx(ijtile)
           istart = istop+1
           istop  = istart+mxadd-1
           jj = (j-1)*ntilef_x+i
!            print *,i,j,jj
           IF (jstag == 1) myadd = myp1(ijtile)
           IF (jstag == 0) myadd = my(ijtile)
!           jstart = jstop+1
           jstop  = jstart+myadd-1
!            write(0,*) 'xsize =',xsize,'  ysize =',ysize,'  myadd,jstag = ',myadd,jstag
!            write(0,*) 'istop =',istop,'  jstop =',jstop, ' ijtile = ',ijtile
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)

            allocate( tem2d(istop-istart+1, jstop-jstart+1) )
          IF (numDims == 3) THEN
            status = NF90_GET_VAR(ncid(jj),varid1,tem2d,(/1,1,k/))
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1,k/))
          ELSE
            status = NF90_GET_VAR(ncid(jj),varid1,tem2d,(/1,1/))
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1/))
          ENDIF

           DO jj = jstart,jstop !Min(ysize,jstop) 
             var2D(istart:istop,jj) = tem2d(1:istop-istart+1, jj-jstart+1) 
           ENDDO
           deallocate( tem2d )
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1,k/),(/istop-istart+1,jstop-jstart+1,1/))
         ENDDO
       ENDDO
!          print *,'xsize =',xsize,'  ysize =',ysize
!          print *,'istop =',istop,'  jstop =',jstop
       IF (numDims == 3) THEN
          status = NF90_PUT_VAR(ncid1,varid2,var2D,(/1,1,kk/)) 
       ELSE
          status = NF90_PUT_VAR(ncid1,varid2,var2D,(/1,1/)) 
       ENDIF
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     ENDDO
     deallocate(var2D)
     deallocate(DimIDV)

   RETURN
   END SUBROUTINE

! #####################################################################################

   SUBROUTINE READ_WRITE_VAR2D_hdf5_old(vid)

     implicit none

     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, tbeg1
!     integer       :: ibeg, jbeg
     real, allocatable :: tem2d(:,:)
     integer ii,jj,ijtile

     real, allocatable :: tmp2d(:,:), arr1d(:)
     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id       ! Dataset identifier 
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
     INTEGER(HID_T) :: memspace      ! memspace identifier 
     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(3) :: offset,count
     INTEGER(HSIZE_T), DIMENSION(2) :: offset_out,count_out,dimsm, data_dims_out
                                            !hyperslab offset in the file 
!     INTEGER :: dsetrank = 2 ! Dataset rank ( in file )
     INTEGER :: memrank = 2  ! Dataset rank ( in memory )

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)
     IF ( .not. append ) status = NF90_REDEF(ncid1)
     write(6,"(/,'Read/Write VAR2D: name = ',a,/)") vname
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,xtype=vtype,ndims=numDims)
     allocate(DimIDV(numDims))
     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
     
     IF ( .not. append ) THEN
       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)
     ELSE
       status = nf90_inq_varid(ncid1, trim(vname), varid2)
     ENDIF

!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 2) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:2),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ENDIF

! grab and write variable attributes
     DO k = 1,numAtts

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     IF ( .not. append ) status = NF90_ENDDEF(ncid1)

     IF (numDims == 3) THEN
       tbeg = ntbeg
       tend = ntend
       tbeg1 = ntbeg
       IF ( one_time ) THEN
         tbeg = ktime
         tend = ktime
       ENDIF
     ELSE
       tbeg = 1
       tend = 1
       tbeg1 = 1
     ENDIF


     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     allocate(var2D(xsize,ysize))

! Find the x,y indices of the starting tile 
!     ibeg = mod(nllc,ntilef_x)
!     jbeg = (nllc-ibeg)/ntiles_y+1
!     jstp = 0
!     DO j = 2,jbeg
!       jj = (j-2)*ntilef_x+1
!       IF (jstag == 1) myadd = myp1(jj)
!       IF (jstag == 0) myadd = my(jj)
!       jstp = jstp + myadd
!     ENDDO
!     istp = 0
!     DO i = 2,ibeg
!       IF (istag == 1) mxadd = mxp1(i)
!       IF (istag == 0) mxadd = mx(i)
!       istp = istp + mxadd
!     ENDDO
     
     DO k = tbeg,tend
       kk = k-tbeg1+ntstart
       var2D(:,:) = 0.0

       jstop=0
       DO j = jbeg,ntiles_y+jbeg-1
         istop=0
         jstart = jstop+1
         DO i = ibeg,ntiles_x+ibeg-1
           ijtile = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(ijtile)
           IF (istag == 0) mxadd = mx(ijtile)
           istart = istop+1
           istop  = istart+mxadd-1
           jj = (j-1)*ntilef_x+i
!            print *,i,j,jj
           IF (jstag == 1) myadd = myp1(ijtile)
           IF (jstag == 0) myadd = my(ijtile)
!           jstart = jstop+1
           jstop  = jstart+myadd-1
!            write(0,*) 'xsize =',xsize,'  ysize =',ysize,'  myadd,jstag = ',myadd,jstag
!            write(0,*) 'istop =',istop,'  jstop =',jstop, ' ijtile = ',ijtile
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)

            allocate( tem2d(istop-istart+1, jstop-jstart+1) )
          IF (numDims == 3) THEN
            status = NF90_GET_VAR(ncid(jj),varid1,tem2d,(/1,1,k/))
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1,k/))
          ELSE
            status = NF90_GET_VAR(ncid(jj),varid1,tem2d,(/1,1/))
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1/))
          ENDIF

           DO jj = jstart,jstop !Min(ysize,jstop) 
             var2D(istart:istop,jj) = tem2d(1:istop-istart+1, jj-jstart+1) 
           ENDDO
           deallocate( tem2d )
!            status = NF90_GET_VAR(ncid(jj),varid1,var2D(istart:istop,jstart:jstop),(/1,1,k/),(/istop-istart+1,jstop-jstart+1,1/))
         ENDDO
       ENDDO
!          print *,'xsize =',xsize,'  ysize =',ysize
!          print *,'istop =',istop,'  jstop =',jstop
       IF (numDims == 3) THEN
          status = NF90_PUT_VAR(ncid1,varid2,var2D,(/1,1,kk/)) 
       ELSE
          status = NF90_PUT_VAR(ncid1,varid2,var2D,(/1,1/)) 
       ENDIF
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     ENDDO
     deallocate(var2D)
     deallocate(DimIDV)

   RETURN
   END SUBROUTINE

! #####################################################################################

   SUBROUTINE DEF_VAR3D(vid)

     implicit none
     
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     integer i1,j1,k1
     integer       :: n

     integer, allocatable :: chunksizes(:)

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

    IF ( .not. append ) status = NF90_REDEF(ncid1)
    write(0,"(/,'DEF_VAR3D: name = ',a,/)") vname
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,xtype=vtype,ndims=numDims)
    allocate(DimIDV(numDims))
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
    write(0,*) 'DEF_VAR3D'
    IF ( .not. append ) THEN
      status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)
    


!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 3) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:4),varid2)
!     ENDIF
#if defined (NC4)
       IF ( netcdfversion == 4 .and. numDims >= 2 ) THEN

        allocate( chunksizes(4) )
        chunksizes(1) =  mxi(1) !  sizes(1:3)
        chunksizes(2) =  myj(1) !  sizes(1:3)
        chunksizes(3) = 1
        
        IF ( numDims == 4 ) THEN
        IF ( vname == 'W' .or. vname == 'EZ' ) THEN
          chunksizes(3) =  mzp1t !  sizes(1:3)
        ELSE
          chunksizes(3) =  mzt !  sizes(1:3)
        ENDIF
            IF ( 4*(chunksizes(3)/4) == chunksizes(3) ) THEN
            chunksizes(3)   = chunksizes(3)/4
            ELSE
            chunksizes(3)   = chunksizes(3)/4 + 1
            ENDIF
        
        ENDIF
        chunksizes(4)   = 1
        
        status = NF90_DEF_VAR_DEFLATE(ncid1, varid2, shuffle, deflate, deflate_level)
!        write(0,*) 'deflate, deflate_level, shuffle = ',deflate, deflate_level, shuffle
        IF (status /= NF90_NOERR) write(0,*)'DEFINE_NCDF_VAR:  Error setting deflate for: ',vname

        status = NF90_DEF_VAR_CHUNKING(ncid1, varid2, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) THEN
          print *,'DEFINE_NCDF_VAR:  Error setting chunking for: ', vname
          write(*,*) 'chunksizes = ',chunksizes
        ENDIF
        deallocate( chunksizes )
       ENDIF
#endif
    ELSE
     status = nf90_inq_varid(ncid1, vname, varid2)
!     status = NF90_INQUIRE_VARIABLE(ncid(1),varid2,vname,vtype)
     write(0,*) 'DEF_VAR: vname,varid2: ', trim(varname), ',  ', varid2
    ENDIF

    varid2arr(vid) = varid2

! grab and write variable attributes
     kstag = 0
     DO k = 1,numAtts
     IF ( .not. append ) write(0,*) 'Put attribute ',k

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         IF ( .not. append ) status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     
     VarStag(1,vid) = istag
     VarStag(2,vid) = jstag
     VarStag(3,vid) = kstag

     write(0,*) 'call NF90_ENDDEF'
     IF ( .not. append ) status = NF90_ENDDEF(ncid1)
     write(0,*) 'done NF90_ENDDEF'

     deallocate(DimIDV)

   RETURN
   END SUBROUTINE

! #####################################################################################


   SUBROUTINE READ_WRITE_VAR3D(vid, nodef)

     implicit none
     
     integer, optional :: nodef
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     integer i1,j1,k1
     integer       :: n
     real          :: commas_dtime
     integer, allocatable :: chunksizes(:)

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

    IF ( .not. present( nodef ) ) THEN
    IF ( .not. append ) status = NF90_REDEF(ncid1)
    ENDIF
    write(0,"(/,'Read/Write VAR3D: name = ',a,/)") vname
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,xtype=vtype,ndims=numDims)
    allocate(DimIDV(numDims))
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)

    IF ( .not. present( nodef ) ) THEN
    write(0,*) 'DEF_VAR'

    IF ( .not. append ) THEN
    status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)


!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 3) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:4),varid2)
!     ENDIF
#if defined (NC4)
       IF ( netcdfversion == 4 .and. numDims >= 2 ) THEN

        allocate( chunksizes(4) )
        chunksizes(1) =  mxi(1) !  sizes(1:3)
        chunksizes(2) =  myj(1) !  sizes(1:3)
        chunksizes(3) = 1
        
        IF ( numDims == 4 ) THEN
        IF ( vname == 'W' .or. vname == 'EZ' ) THEN
          chunksizes(3) =  mzp1t !  sizes(1:3)
        ELSE
          chunksizes(3) =  mzt !  sizes(1:3)
        ENDIF

            IF ( 4*(chunksizes(3)/4) == chunksizes(3) ) THEN
            chunksizes(3)   = chunksizes(3)/4
            ELSE
            chunksizes(3)   = chunksizes(3)/4 + 1
            ENDIF

        ENDIF
        chunksizes(4)   = 1
        
        status = NF90_DEF_VAR_DEFLATE(ncid1, varid2, shuffle, deflate, deflate_level)
!        write(0,*) 'deflate, deflate_level, shuffle = ',deflate, deflate_level, shuffle
        IF (status /= NF90_NOERR) write(0,*)'DEFINE_NCDF_VAR:  Error setting deflate for: ',vname

        status = NF90_DEF_VAR_CHUNKING(ncid1, varid2, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) THEN
         print *,'DEFINE_NCDF_VAR:  Error setting chunking for: ', vname
         write(*,*) 'chunksizes = ',chunksizes
        ENDIF
        deallocate( chunksizes )
       ENDIF
#endif
    ELSE
     status = nf90_inq_varid(ncid1, vname, varid2)
!     status = NF90_INQUIRE_VARIABLE(ncid(1),varid2,vname,vtype)
     write(0,*) 'DEF_VAR: vname,varid2: ', trim(varname), ',  ', varid2
    ENDIF
! grab and write variable attributes
     DO k = 1,numAtts
     write(0,*) 'Put attribute ',k

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     write(0,*) 'call NF90_ENDDEF'
     status = NF90_ENDDEF(ncid1)
     write(0,*) 'done NF90_ENDDEF'

    ENDIF ! .not. present

     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     IF (kstag == 1) zsize = mzp1t
     IF (kstag == 0) zsize = mzt
     allocate(var3D(xsize,ysize,zsize))

      write(0,*) 'numDims,xsize,ysize,zsize = ',numDims,xsize,ysize,zsize
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

     IF ( one_time ) THEN
       tbeg = ktime
       tend = ktime
     ENDIF

     DO k = tbeg,tend
       kk = k-tbeg+ntstart
       var3D(:,:,:) = 0.0

!       jstop=0
!       DO j = jbeg,ntiles_y+jbeg-1
!         istop=0
!         jstart = jstop+1
!         DO i = ibeg,ntiles_x+ibeg-1
!           jj = (j-1)*ntilef_x+i
!           IF (istag == 1) mxadd = mxp1(jj)
!           IF (istag == 0) mxadd = mx(jj)
!           istart = istop+1
!           istop  = istart+mxadd-1
!           IF (jstag == 1) myadd = myp1(jj)
!           IF (jstag == 0) myadd = my(jj)
!            jstop  = jstart+myadd-1
           tim1r = commas_dtime(timarrr)

       jstop = 0
       DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO


           status = NF90_INQ_VARID(ncid(jj),vname,varid1)
           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
         IF ( openclose ) THEN
         
         ENDIF
         
         ENDDO
       ENDDO
           tim2r = commas_dtime(timarrr)
           timeread = timeread + tim2r

           tim1w = commas_dtime(timarrw)
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       status = NF90_PUT_VAR(ncid1,varid2,var3D,(/1,1,1,kk/)) 
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     vmax = -1.0e20
     vmin =  1.0e20
     DO k1 = 1,zsize
       DO j1 = 1,ysize
         DO i1 = 1,xsize
          vmax = Max( vmax, var3D(i1,j1,k1) )
          vmin = Min( vmin, var3D(i1,j1,k1) )
         ENDDO
       ENDDO
     ENDDO
           tim2w = commas_dtime(timarrw)
           timewrite = timewrite + tim2w
     
     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
     write(6,*) 'TIME FOR VAR read/write',trim(vname),' = ',tim2r,tim2w
     ENDDO
     
     
     deallocate(var3D)
     deallocate(DimIDV)
         
   RETURN
   END SUBROUTINE

! #####################################################################################


   SUBROUTINE READ_WRITE_VAR2D_hdf5(vid,vname)

     implicit none
     
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     real  vmax2, vmin2
     integer i1,j1,k1
     integer       :: n
     real          :: commas_dtime
     integer       :: error
     real, allocatable :: tmp2d(:,:), arr1d(:)
     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id,dset_id1        ! Dataset identifier 
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
     INTEGER(HID_T) :: memspace      ! memspace identifier 
     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims,size
     INTEGER(HSIZE_T), DIMENSION(3) :: offset,count
     INTEGER(HSIZE_T), DIMENSION(2) :: offset_out,count_out,dimsm, data_dims_out
                                            !hyperslab offset in the file 
!     INTEGER :: dsetrank = 2 ! Dataset rank ( in file )
     INTEGER :: memrank = 2  ! Dataset rank ( in memory )

!     vname = VarNameArr(vid)
     
     IF ( .not. writehdf5 ) THEN
     status = nf90_inq_varid(ncid1, trim(vname), varid2)
     
     status = NF90_INQUIRE_VARIABLE(ncid1,varid2,vname,vtype,ndims=numDims,nAtts=numAtts)
     
     ELSE
     ! open dataset
            CALL h5dopen_f(ncid1, vname, dset_id1, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 dataset. dset_id1, error = ',dset_id,error
             STOP
           ENDIF
           numDims = 4
     ENDIF
!     status = nf90_inq_varid(ncid1, trim(vname), varid2)
     
!     status = NF90_INQUIRE_VARIABLE(ncid1,varid2,vname,vtype,ndims=numDims,nAtts=numAtts)
     write(*,*) 'var2d_hdf5: vid, varid2, vname = ',vid, varid2, vname

!    status = NF90_REDEF(ncid1)
!    write(0,"(/,'Read/Write VAR3D: name = ',a,/)") vname
!    status = NF90_INQUIRE_VARIABLE(ncid1,vid,xtype=vtype,ndims=numDims)
!    allocate(DimIDV(numDims))
!    status = NF90_INQUIRE_VARIABLE(ncid1,vid,dimids=DimIDV)
!    write(0,*) 'DEF_VAR'
!    status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)


!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 3) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:4),varid2)
!     ENDIF


     istag = VarStag(1,vid)
     jstag = VarStag(2,vid)
     kstag = VarStag(3,vid)
     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     allocate(var2D(xsize,ysize))

      write(0,*) 'numDims,xsize,ysize = ',numDims,xsize,ysize
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

     IF ( one_time .or. one_time_at_a_time ) THEN
       tbeg = ktime
       tend = ktime
     ENDIF

     DO k = tbeg,tend
       kk = k-ntbeg+ntstart
       var2D(:,:) = 0.0

!       jstop=0
!       DO j = jbeg,ntiles_y+jbeg-1
!         istop=0
!         jstart = jstop+1
!         DO i = ibeg,ntiles_x+ibeg-1
!           jj = (j-1)*ntilef_x+i
!           IF (istag == 1) mxadd = mxp1(jj)
!           IF (istag == 0) mxadd = mx(jj)
!           istart = istop+1
!           istop  = istart+mxadd-1
!           IF (jstag == 1) myadd = myp1(jj)
!           IF (jstag == 0) myadd = my(jj)
!            jstop  = jstart+myadd-1
           tim1r = commas_dtime(timarrr)

       jstop = 0
       DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO

!         CALL h5fopen_f (fname(jj), H5F_ACC_RDONLY_F, file_id, error)
         CALL h5fopen_f (fname(jj), H5F_ACC_RDONLY_F, hdf5_fid(jj), error)
         IF ( error /= 0 ) THEN
           write(0,*) 'error opening hdf5 file, error = ',error
           STOP
         ENDIF

! test check reading TIME array
!   open dataset
!            CALL h5dopen_f(file_id, 'TIME', dset_id, error)
!           IF ( error /= 0 ) THEN
!             write(0,*) 'error opening hdf5 dataset. dset_id, error = ',dset_id,error
!             STOP
!           ENDIF
!           
!           allocate( arr1d(nt) )
!
!            data_dims(1) = nt
!            CALL h5dread_f(dset_id, H5T_NATIVE_REAL, arr1d, data_dims, error)
!           IF ( error /= 0 ) THEN
!             write(0,*) 'error reading TIME in hdf5 dataset. dset_id, error = ',dset_id,error
!             STOP
!           ELSE
!             write(0,*) arr1d(1:5)
!           ENDIF
!           
!           deallocate ( arr1d )
!     !
!     ! Close the dataset.
!     !
!            CALL h5dclose_f(dset_id, error)
!            

!   open dataset
!            CALL h5dopen_f(file_id, vname, dset_id, error)
            CALL h5dopen_f(hdf5_fid(jj), vname, dset_id, error)
!            status = NF90_INQ_VARID(ncid(jj),vname,varid1)
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/))
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 dataset. dset_id, error = ',dset_id,error
             STOP
           ENDIF

     !
     ! Get dataset's dataspace identifier.
     !
            CALL h5dget_space_f(dset_id, dataspace, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dget_space_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF

     !
     ! Select hyperslab in the dataset.
     !
           offset = (/0,0,k-1/)
           count(1) = istop-istart+1
           count(2) = jstop-jstart+1
           count(3) = 1

           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error) 
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f - data, error, = ',error,data_dims
             STOP
           ENDIF
     !
     ! Create memory dataspace.
     !
            dimsm(1:2) = count(1:2)
            CALL h5screate_simple_f(memrank, dimsm, memspace, error)

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5screate_simple_f, error, = ',error,data_dims
             STOP
           ENDIF

     !
     ! Select hyperslab in memory.
     !
            offset_out = (/ 0,0 /)
           count_out(1) = istop-istart+1
           count_out(2) = jstop-jstart+1
            
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                offset_out, count_out, error) 
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f - memory, error, = ',error,data_dims
             STOP
           ENDIF

           
           allocate( tmp2d(istop-istart+1,jstop-jstart+1) )
           data_dims(1) = istop-istart+1
           data_dims(2) = jstop-jstart+1
!           data_dims(3) = zsize
!           data_dims(4) = nt
!           CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tmp3d, data_dims, error)
           CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tmp2d, data_dims, error, &
                    memspace, dataspace)
           IF ( error /= 0 ) THEN
             write(0,*) 'error reading hdf5 file, error, = ',error,data_dims(1:2),memspace
             STOP
           ENDIF
           
           var2D(istart:istop,jstart:jstop) = tmp2d(1:istop-istart+1,1:jstop-jstart+1)
!           var3D(istart:istop,jstart:jstop,:) = tmp2d(1:istop-istart+1,1:jstop-jstart+1,1:zsize,k)
     !
     ! Close the dataspace for the dataset.
     !
           CALL h5sclose_f(dataspace, error)

     !
     ! Close the memoryspace.
     !
           CALL h5sclose_f(memspace, error)

           
     !
     ! Close the dataset.
     !
            CALL h5dclose_f(dset_id, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 file, error = ',error
             STOP
           ENDIF

     !
     ! Close the file.
     !
!            CALL h5fclose_f(file_id, error)
            CALL h5fclose_f(hdf5_fid(jj), error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 file, error = ',error
             STOP
           ENDIF
            
            deallocate(tmp2d)
           

!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
         ENDDO
       ENDDO
           tim2r = commas_dtime(timarrr)
           timeread = timeread + tim2r

           tim1w = commas_dtime(timarrw)
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop

       IF ( writehdf5 ) THEN
       
         allocate( tmp2d(xsize,ysize) )

     !
     ! Get dataset's dataspace identifier.
     !
            CALL h5dget_space_f(dset_id1, dataspace, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dget_space_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF

           offset = (/0,0,kk-1/)
           count(1) = xsize
           count(2) = ysize
           count(3) = 1
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error) 

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF
     !
     ! Create memory dataspace.
     !
            dimsm(1:2) = count(1:2)
            CALL h5screate_simple_f(memrank, dimsm, memspace, error)

     !
     ! Select hyperslab in memory.
     !
            offset_out = (/ 0,0 /)
           count_out(1) = xsize
           count_out(2) = ysize
            
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                offset_out, count_out, error) 
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f - memory, error, = ',error,data_dims
             STOP
           ENDIF
         
         data_dims(1:2) = count(1:2)
         data_dims(3) = 1

!         CALL h5dread_f(dset_id1, H5T_NATIVE_REAL, tmp2d, data_dims, error, &
!                          memspace, dataspace)
!          write(0,*) 'h5dread_f, error = ',error


          size(1)   = count_out(1)
          size(2)   = count_out(2) 
          size(3)   = kk

          CALL h5dextend_f(dset_id1, size, error)

         CALL h5dwrite_f(dset_id1, H5T_NATIVE_REAL, var2D, data_dims, error, &
                          memspace, dataspace)

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dwrite_f var2D error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF

     !
     ! Close the dataspace for the dataset.
     !
           CALL h5sclose_f(dataspace, error)

     !
     ! Close the memoryspace.
     !
           CALL h5sclose_f(memspace, error)

          deallocate( tmp2d )
       ELSE

       status = NF90_PUT_VAR(ncid1,varid2,var2D,(/1,1,kk/)) 
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF

#ifndef __ia64__
! for some reason this is causing problems on itanium systems.  
!         CALL DISP_ERR( nf90_close(ncid1) , .true. )
!         CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)
#endif
       
      ENDIF

     vmax = -1.0e20
     vmin =  1.0e20
       DO j1 = 1,ysize
         DO i1 = 1,xsize
          vmax = Max( vmax, var2D(i1,j1) )
          vmin = Min( vmin, var2D(i1,j1) )
         ENDDO
       ENDDO
           tim2w = commas_dtime(timarrw)
           timewrite = timewrite + tim2w

     IF ( .false. .and. .not. writehdf5 ) THEN
      var2D(:,:) = 0.0
      status = NF90_GET_VAR(ncid1,varid2,var2D,(/1,1,kk/)) 
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     vmax2 = -1.0e20
     vmin2 =  1.0e20
       DO j1 = 1,ysize
         DO i1 = 1,xsize
          vmax2 = Max( vmax2, var2D(i1,j1) )
          vmin2 = Min( vmin2, var2D(i1,j1) )
         ENDDO
       ENDDO
       
     ENDIF
     
     
     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
!     write(0,*) 'k,kk, variable max2/min2 = ',k,kk,vmax2,vmin2
     write(6,*) 'TIME FOR VAR read/write ',trim(vname),' = ',tim2r,tim2w
     ENDDO
     
     
     deallocate(var2D)
!     deallocate(DimIDV)
         
   RETURN
   END SUBROUTINE


! #####################################################################################


   SUBROUTINE READ_WRITE_VAR3D_hdf5(vid,vname)

     implicit none
     
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     real  vmax2, vmin2
     integer i1,j1,k1
     integer       :: n
     real          :: commas_dtime
     integer       :: error
     real, allocatable :: tmp3d(:,:,:), arr1d(:)
     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id,dset_id1       ! Dataset identifier 
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
     INTEGER(HID_T) :: memspace      ! memspace identifier 
     INTEGER(HSIZE_T), DIMENSION(4) :: data_dims
     INTEGER(HSIZE_T), DIMENSION(4) :: offset,count
     INTEGER(HSIZE_T), DIMENSION(3) :: offset_out,count_out,dimsm, data_dims_out
                                            !hyperslab offset in the file 
!     INTEGER :: dsetrank = 2 ! Dataset rank ( in file )
     INTEGER :: memrank = 3  ! Dataset rank ( in memory )
     
     real tim1ra,tim1rb,tim1rc,tim1rd,tim1re
     real tim2ra,tim2rb,tim2rc,tim2rd,tim2re
     real tim2rat, tim2rbt, tim2rct
     
     integer mem0, mem1, mem2, mem3, mem4, mem5, mem6, mem7, mem8, mem9, mem10, mem11
#ifndef GETMEM
     logical :: checkmem = .false.
#else
     logical :: checkmem = .true.
#endif

     IF ( checkmem ) CALL get_mem_used(mem0)

!     vname = VarNameArr(vid)
     varid2 = 0
     
     IF ( .not. writehdf5 ) THEN
     status = nf90_inq_varid(ncid1, trim(vname), varid2)
     
     status = NF90_INQUIRE_VARIABLE(ncid1,varid2,vname,vtype,ndims=numDims,nAtts=numAtts)
     
     ELSE
            CALL h5dopen_f(ncid1, vname, dset_id1, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 dataset. dset_id1, error = ',dset_id,error
             STOP
           ENDIF
     ENDIF
     write(*,*) 'var3d_hdf5: vid, varid2, vname = ',vid, varid2, vname

!    status = NF90_REDEF(ncid1)
!    write(0,"(/,'Read/Write VAR3D: name = ',a,/)") vname
!    status = NF90_INQUIRE_VARIABLE(ncid1,vid,xtype=vtype,ndims=numDims)
!    allocate(DimIDV(numDims))
!    status = NF90_INQUIRE_VARIABLE(ncid1,vid,dimids=DimIDV)
!    write(0,*) 'DEF_VAR'
!    status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)


!     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
!     IF (numDims == 3) THEN
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:3),varid2)
!     ELSE
!       status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV(1:4),varid2)
!     ENDIF

     istag = VarStag(1,vid)
     jstag = VarStag(2,vid)
     kstag = VarStag(3,vid)

     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     IF (kstag == 1) zsize = mzp1t
     IF (kstag == 0) zsize = mzt
     allocate(var3D(xsize,ysize,zsize))

      write(0,*) 'numDims,xsize,ysize,zsize = ',numDims,xsize,ysize,zsize
      write(0,*) 'istag,jstag,kstag = ',istag,jstag,kstag
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

     IF ( one_time .or. one_time_at_a_time ) THEN
       tbeg = ktime
       tend = ktime
     ENDIF

     DO k = tbeg,tend
       kk = k-ntbeg+ntstart
       var3D(:,:,:) = 0.0

!       jstop=0
!       DO j = jbeg,ntiles_y+jbeg-1
!         istop=0
!         jstart = jstop+1
!         DO i = ibeg,ntiles_x+ibeg-1
!           jj = (j-1)*ntilef_x+i
!           IF (istag == 1) mxadd = mxp1(jj)
!           IF (istag == 0) mxadd = mx(jj)
!           istart = istop+1
!           istop  = istart+mxadd-1
!           IF (jstag == 1) myadd = myp1(jj)
!           IF (jstag == 0) myadd = my(jj)
!            jstop  = jstart+myadd-1
           tim1r = 0.0 ! commas_dtime(timarrr)
           tim2r = 0.0
           tim2rat = 0.0
           tim2rbt = 0.0
           tim2rct = 0.0



       jstop = 0
       DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO

           tim1ra = commas_dtime(timarrr)

          IF ( checkmem ) CALL get_mem_used(mem1)
          
!         CALL h5fopen_f (fname(jj), H5F_ACC_RDONLY_F, file_id, error)
         CALL h5fopen_f (fname(jj), H5F_ACC_RDONLY_F, hdf5_fid(jj), error)
         IF ( error /= 0 ) THEN
           write(0,*) 'error opening hdf5 file, error = ',error
           STOP
         ENDIF

          IF ( checkmem ) CALL get_mem_used(mem2)

! test check reading TIME array
!   open dataset
!            CALL h5dopen_f(file_id, 'TIME', dset_id, error)
!           IF ( error /= 0 ) THEN
!             write(0,*) 'error opening hdf5 dataset. dset_id, error = ',dset_id,error
!             STOP
!           ENDIF
!           
!           allocate( arr1d(nt) )
!
!            data_dims(1) = nt
!            CALL h5dread_f(dset_id, H5T_NATIVE_REAL, arr1d, data_dims, error)
!           IF ( error /= 0 ) THEN
!             write(0,*) 'error reading TIME in hdf5 dataset. dset_id, error = ',dset_id,error
!             STOP
!           ELSE
!             write(0,*) arr1d(1:5)
!           ENDIF
!           
!           deallocate ( arr1d )
!     !
!     ! Close the dataset.
!     !
!            CALL h5dclose_f(dset_id, error)
!            

!   open dataset
!            CALL h5dopen_f(file_id, vname, dset_id, error)
            CALL h5dopen_f(hdf5_fid(jj), vname, dset_id, error)
!            status = NF90_INQ_VARID(ncid(jj),vname,varid1)
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/))
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 dataset. dset_id, error = ',dset_id,error
             STOP
           ENDIF

          IF ( checkmem ) CALL get_mem_used(mem3)

     !
     ! Get dataset's dataspace identifier.
     !
            CALL h5dget_space_f(dset_id, dataspace, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dget_space_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF


     !
     ! Select hyperslab in the dataset.
     !
           offset = (/0,0,0,k-1/)
           count(1) = istop-istart+1
           count(2) = jstop-jstart+1
           count(3) = zsize
           count(4) = 1

           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error) 
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f - data, error, = ',error,data_dims
             STOP
           ENDIF
     !
     ! Create memory dataspace.
     !
            dimsm(1:3) = count(1:3)
            CALL h5screate_simple_f(memrank, dimsm, memspace, error)

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5screate_simple_f, error, = ',error,data_dims
             STOP
           ENDIF

            IF ( checkmem ) CALL get_mem_used(mem4)
   !
     ! Select hyperslab in memory.
     !
            offset_out = (/ 0,0,0 /)
           count_out(1) = istop-istart+1
           count_out(2) = jstop-jstart+1
           count_out(3) = zsize
            
            CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                offset_out, count_out, error) 
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f - memory, error, = ',error,data_dims
             STOP
           ENDIF

          IF ( checkmem ) CALL get_mem_used(mem5)

           tim2ra = commas_dtime(timarrr)

           tim1rb = commas_dtime(timarrr)

           
           allocate( tmp3d(istop-istart+1,jstop-jstart+1,zsize) )
           data_dims(1) = istop-istart+1
           data_dims(2) = jstop-jstart+1
           data_dims(3) = zsize

          IF ( checkmem ) CALL get_mem_used(mem6)

!           data_dims(4) = nt
!           CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tmp3d, data_dims, error)
           CALL h5dread_f(dset_id, H5T_NATIVE_REAL, tmp3d, data_dims, error, &
                    memspace, dataspace)
           IF ( error /= 0 ) THEN
             write(0,*) 'error reading hdf5 file, error, = ',error,data_dims(1:3),memspace
             STOP
           ENDIF
           
           var3D(istart:istop,jstart:jstop,:) = tmp3d(1:istop-istart+1,1:jstop-jstart+1,1:zsize)
!           var3D(istart:istop,jstart:jstop,:) = tmp3d(1:istop-istart+1,1:jstop-jstart+1,1:zsize,k)

           tim2rb = commas_dtime(timarrr)

           tim1rc = commas_dtime(timarrr)

     !
     ! Close the dataspace for the dataset.
     !
           CALL h5sclose_f(dataspace, error)

          IF ( checkmem ) CALL get_mem_used(mem7)
     !
     ! Close the memoryspace.
     !
           CALL h5sclose_f(memspace, error)

          IF ( checkmem ) CALL get_mem_used(mem8)
           
     !
     ! Close the dataset.
     !
            CALL h5dclose_f(dset_id, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 file, error = ',error
             STOP
           ENDIF

          IF ( checkmem ) CALL get_mem_used(mem9)

     !
     ! Close the file.
     !
!            CALL h5fclose_f(file_id, error)
            CALL h5fclose_f(hdf5_fid(jj), error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error opening hdf5 file, error = ',error
             STOP
           ENDIF
            
          IF ( checkmem ) CALL get_mem_used(mem10)

            deallocate(tmp3d)

          IF ( checkmem ) CALL get_mem_used(mem11)
          
          IF ( checkmem ) THEN
            write(6,*) 'hdf5 memory usage: ',mem1,mem2,mem3,mem4,mem5,mem6,mem7,mem8,mem9,mem10,mem11
            write(6,*) 'hdf5 memory delta: ',mem11-mem1,mem2-mem1,mem3-mem2,mem4-mem3,mem5-mem4,mem6-mem5,mem7-mem6,mem8-mem7,  &
              mem9-mem8,mem10-mem9,mem11-mem10
          ENDIF
           
           tim2rc = commas_dtime(timarrr)
           
           tim2r = tim2r + tim2ra + tim2rb + tim2rc
           tim2rat = tim2rat + tim2ra
           tim2rbt = tim2rbt + tim2rb
           tim2rct = tim2rct + tim2rc

!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
         ENDDO
       ENDDO
!           tim2r = commas_dtime(timarrr)
           timeread = timeread + tim2r

          IF ( checkmem ) CALL get_mem_used(mem1)

           tim1w = commas_dtime(timarrw)
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       
       IF ( writehdf5 ) THEN

     !
     ! Get dataset's dataspace identifier.
     !
            CALL h5dget_space_f(dset_id1, dataspace, error)
           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dget_space_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF

           offset = (/0,0,0,kk-1/)
           count(1) = xsize
           count(2) = ysize
           count(3) = zsize
           count(4) = 1
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error) 

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5sselect_hyperslab_f error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF
     !
     ! Create memory dataspace.
     !
            dimsm(1:3) = count(1:3)
            CALL h5screate_simple_f(memrank, dimsm, memspace, error)

           data_dims(1:4) = count(1:4)

         CALL h5dwrite_f(dset_id1, H5T_NATIVE_REAL, var3D, data_dims, error, &
                          mem_space_id=memspace, file_space_id=dataspace)

           IF ( error /= 0 ) THEN
             write(0,*) 'error h5dwrite_f var3D error = ',error,data_dims
             STOP
           ELSE
!             write(*,*) 'dataspace, dset_id, vname = ',dataspace, dset_id, vname
           ENDIF

     !
     ! Close the dataspace for the dataset.
     !
           CALL h5sclose_f(dataspace, error)

     !
     ! Close the memoryspace.
     !
           CALL h5sclose_f(memspace, error)

         
       ELSE
         status = NF90_PUT_VAR(ncid1,varid2,var3D,(/1,1,1,kk/)) 
         IF (status /= nf90_noerr) THEN
           print *, 'PUT_VAR  ',trim(nf90_strerror(status))
         ENDIF

#ifndef __ia64__
!         CALL DISP_ERR( nf90_close(ncid1) , .true. )
!         CALL DISP_ERR( NF90_OPEN(fname1,nf90_write,ncid1), .true.)
#endif
       
       ENDIF

          IF ( checkmem ) CALL get_mem_used(mem2)

     vmax = -1.0e20
     vmin =  1.0e20
     DO k1 = 1,zsize
       DO j1 = 1,ysize
         DO i1 = 1,xsize
          vmax = Max( vmax, var3D(i1,j1,k1) )
          vmin = Min( vmin, var3D(i1,j1,k1) )
         ENDDO
       ENDDO
     ENDDO
           tim2w = commas_dtime(timarrw)
           timewrite = timewrite + tim2w
!           timeread = timeread + tim2r

!      var3D(:,:,:) = 0.0
!      status = NF90_GET_VAR(ncid1,varid2,var3D,(/1,1,1,kk/)) 
!       IF (status /= nf90_noerr) THEN
!         print *, 'GET_VAR  ',trim(nf90_strerror(status))
!       ENDIF
!     vmax2 = -1.0e20
!     vmin2 =  1.0e20
!     DO k1 = 1,zsize
!       DO j1 = 1,ysize
!         DO i1 = 1,xsize
!          vmax2 = Max( vmax2, var3D(i1,j1,k1) )
!          vmin2 = Min( vmin2, var3D(i1,j1,k1) )
!         ENDDO
!       ENDDO
!     ENDDO
     
     
     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
!     write(0,*) 'k,kk, variable max2/min2 = ',k,kk,vmax2,vmin2
     write(6,*) 'Times for open,read,close = ',tim2rat,tim2rbt,tim2rct
     write(6,*) 'TIME FOR VAR read/write ',trim(vname),' = ',tim2r,tim2w
     ENDDO
     
     
     deallocate(var3D)
!     deallocate(DimIDV)

          IF ( checkmem ) CALL get_mem_used(mem3)

          IF ( checkmem ) THEN
            write(6,*) 'hdf5sub memory usage: ',mem0,mem3-mem0
            write(6,*) 'nc4 write memory: ',mem2-mem1
          ENDIF
         
   RETURN
   END SUBROUTINE


! #####################################################################################


   SUBROUTINE FILEREAD_WRITE_3D(vid,jfile)

     implicit none
     
     integer jfile
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     integer i1,j1,k1
     integer       :: n

     status = NF90_INQUIRE_VARIABLE(ncid(jfile),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

     istag = VarStag(1,vid)
     jstag = VarStag(2,vid)
     kstag = VarStag(3,vid)
     
     IF (istag == 1) xsize = mxp1(jfile)
     IF (istag == 0) xsize = mx(jfile)
     IF (jstag == 1) ysize = myp1(jfile)
     IF (jstag == 0) ysize = my(jfile)
     IF (kstag == 1) zsize = mzp1t
     IF (kstag == 0) zsize = mzt

!      write(0,*) 'vname,numDims,xsize,ysize,zsize = ',vname,numDims,xsize,ysize,zsize
      write(0,*) 'vname = ',vname
      
      allocate( var3D(xsize,ysize,zsize) )
      
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

     IF ( one_time .or. one_time_at_a_time ) THEN
       tbeg = ktime
       tend = ktime
     ENDIF

     DO k = tbeg,tend
       kk = k-tbeg+ntstart

!       jstop=0
!       DO j = jbeg,ntiles_y+jbeg-1
!         istop=0
!         jstart = jstop+1
!         DO i = ibeg,ntiles_x+ibeg-1
!           jj = (j-1)*ntilef_x+i
!           IF (istag == 1) mxadd = mxp1(jj)
!           IF (istag == 0) mxadd = mx(jj)
!           istart = istop+1
!           istop  = istart+mxadd-1
!           IF (jstag == 1) myadd = myp1(jj)
!           IF (jstag == 0) myadd = my(jj)
!            jstop  = jstart+myadd-1

       jstop = 0
       jj = -1
       outer: DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO


           IF ( jj .eq. jfile ) THEN
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)
           status = NF90_GET_VAR(ncid(jj),varid1,var3D,(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
            EXIT outer
           ENDIF
         ENDDO
       ENDDO outer
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       IF ( jj .ne. -1 ) THEN
         status = NF90_PUT_VAR(ncid1,varid2arr(vid),var3D,(/istart,jstart,1,kk/)) 
         IF (status /= nf90_noerr) THEN
           print *, 'PUT_VAR  ',trim(nf90_strerror(status))
         ENDIF
!     vmax = -1.0e20
!     vmin =  1.0e20
!     DO k1 = 1,zsize
!       DO j1 = 1,ysize
!         DO i1 = 1,xsize
!          vmax = Max( vmax, var3D(i1,j1,k1) )
!          vmin = Min( vmin, var3D(i1,j1,k1) )
!         ENDDO
!       ENDDO
!     ENDDO

       ENDIF
     
!     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
     ENDDO
     
     deallocate( var3D )
     
         
   RETURN

   END SUBROUTINE



! ####################################################################################

   SUBROUTINE READ_WRITE_VAR2D_onetime(vid,jfile,k,kk)

     implicit none

     integer jfile
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, tbeg1
!     integer       :: ibeg, jbeg
     real, allocatable :: tem2d(:,:)
     integer ii,jj,ijtile, found, n
     real  commas_dtime

!     status = NF90_INQUIRE_VARIABLE(ncid(jfile),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

     vname = VarNameArr(vid)
     
     istag = VarStag(1,vid)
     jstag = VarStag(2,vid)
     
     IF (istag == 1) xsize = mxp1(jfile)
     IF (istag == 0) xsize = mx(jfile)
     IF (jstag == 1) ysize = myp1(jfile)
     IF (jstag == 0) ysize = my(jfile)

      IF ( vname .eq. 'RAIN_RAT' .or. jfile .eq. 1 ) write(0,*) 'readwrite2d: vname,k = ',vname,kk
!       write(0,*) 'readwrite2d: vname,k = ',vname,kk

       jstop = 0
       jj = -1
       found = 0
       outer: DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO


           IF ( jj .eq. jfile ) THEN
           found = 1
           tim1r = commas_dtime(timarrr)
           allocate( tem2d(xsize,ysize) )
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)
           status = NF90_GET_VAR(ncid(jj),varid1,tem2d,(/1,1,k/))
           tim2r = commas_dtime(timarrr)
           timeread = timeread + tim2r
            EXIT outer
           ENDIF
         ENDDO
       ENDDO outer
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       IF ( found == 1 ) THEN
           tim1w = commas_dtime(timarrr)
           IF ( writeout ) status = NF90_PUT_VAR(ncid1,varid2arr(vid),tem2d,(/istart,jstart,kk/)) 
           tim2w = commas_dtime(timarrr)
           timewrite = timewrite + tim2w
         IF (status /= nf90_noerr) THEN
           print *, 'PUT_VAR  ',trim(nf90_strerror(status))
         ENDIF
       ELSE
         file_used(jfile) = 0
       ENDIF

     IF ( allocated( tem2d ) ) deallocate(tem2d)

   RETURN
   END SUBROUTINE

! #####################################################################################


   SUBROUTINE FILEREAD_WRITE_3D_one_time(vid,jfile,k,kk)

     implicit none
     
     integer jfile
     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype, found
     real  vmax, vmin
     integer i1,j1,k1
     integer       :: n
     real          :: commas_dtime

!     status = NF90_INQUIRE_VARIABLE(ncid(jfile),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

     vname = VarNameArr(vid)

     istag = VarStag(1,vid)
     jstag = VarStag(2,vid)
     kstag = VarStag(3,vid)
     
     IF (istag == 1) xsize = mxp1(jfile)
     IF (istag == 0) xsize = mx(jfile)
     IF (jstag == 1) ysize = myp1(jfile)
     IF (jstag == 0) ysize = my(jfile)
     IF (kstag == 1) zsize = mzp1t
     IF (kstag == 0) zsize = mzt

!      write(0,*) 'vname,numDims,xsize,ysize,zsize = ',vname,numDims,xsize,ysize,zsize
      IF ( vname .eq. 'DBZ' .or. jfile .eq. 1 ) write(0,*) 'vname,k = ',vname,kk
      
      
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

!     DO k = tbeg,tend
!       kk = k-tbeg+1

!       jstop=0
!       DO j = jbeg,ntiles_y+jbeg-1
!         istop=0
!         jstart = jstop+1
!         DO i = ibeg,ntiles_x+ibeg-1
!           jj = (j-1)*ntilef_x+i
!           IF (istag == 1) mxadd = mxp1(jj)
!           IF (istag == 0) mxadd = mx(jj)
!           istart = istop+1
!           istop  = istart+mxadd-1
!           IF (jstag == 1) myadd = myp1(jj)
!           IF (jstag == 0) myadd = my(jj)
!            jstop  = jstart+myadd-1

       jstop = 0
       jj = -1
       found = 0
       outer: DO j = jbeg,ntiles_y+jbeg-1
         istop = 0
         jstop = 0
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(jj)
           IF (istag == 0) mxadd = mx(jj)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)

           jstop = 0
           DO n = jbeg,j
            jj = (n-1)*ntilef_x+i
            jstart = jstop+1
            IF (jstag == 1) myadd = myp1(jj)
            IF (jstag == 0) myadd = my(jj)         
            jstop  = jstart+myadd-1
           ENDDO


           IF ( jj .eq. jfile ) THEN
           found = 1
           tim1r = commas_dtime(timarrr)
           allocate( var3D(xsize,ysize,zsize) )
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)
           status = NF90_GET_VAR(ncid(jj),varid1,var3D,(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
           tim2r = commas_dtime(timarrr)
           timeread = timeread + tim2r
            EXIT outer
           ENDIF
         ENDDO
       ENDDO outer
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       IF ( found == 1 ) THEN
           tim1w = commas_dtime(timarrr)
           IF ( writeout ) status = NF90_PUT_VAR(ncid1,varid2arr(vid),var3D,(/istart,jstart,1,kk/)) 
           tim2w = commas_dtime(timarrr)
           timewrite = timewrite + tim2w
         IF (status /= nf90_noerr) THEN
           print *, 'PUT_VAR  ',trim(nf90_strerror(status))
         ENDIF
       ELSE
         file_used(jfile) = 0
       ENDIF
!     vmax = -1.0e20
!     vmin =  1.0e20
!     DO k1 = 1,zsize
!       DO j1 = 1,ysize
!         DO i1 = 1,xsize
!          vmax = Max( vmax, var3D(i1,j1,k1) )
!          vmin = Min( vmin, var3D(i1,j1,k1) )
!         ENDDO
!       ENDDO
!     ENDDO
     
!     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
!     ENDDO
     
     IF ( allocated( var3D ) ) deallocate( var3D )
     
         
   RETURN
   END SUBROUTINE

! #####################################################################################

   SUBROUTINE READ_WRITE_VAR3Dtst(vid)

     character(len=40) :: vname, AttName, AttValC
     character(len=1)  :: vdim
     integer       :: AttValI 
!     integer, dimension(4) :: DimIDV
     integer, dimension(:), allocatable :: DimIDV
     integer       :: vid, tbeg, tend, numDims, numAtts, kk, k, vtype, atype
     real  vmax, vmin
     integer i1,j1,k1
     real, allocatable :: tem3d(:,:,:)
     integer, allocatable :: chunksizes(:)

     status = NF90_INQUIRE_VARIABLE(ncid(1),vid,vname,vtype,ndims=numDims,nAtts=numAtts)

!    status = NF90_REDEF(ncid1)
    write(0,"(/,'Read/Write VAR3D: name = ',a,/)") vname
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,xtype=vtype,ndims=numDims)
    allocate(DimIDV(numDims))
    status = NF90_INQUIRE_VARIABLE(ncid(1),vid,dimids=DimIDV)
    write(0,*) 'DEF_VAR'
!    status = NF90_DEF_VAR(ncid1,trim(vname),vtype,DimIDV,varid2)

     IF ( .false. ) THEN

#if defined (NC4)
       IF ( netcdfversion == 4 ) THEN
        allocate( chunksizes(4) )
        chunksizes(1) =  mxi(1) !  sizes(1:3)
        chunksizes(2) =  myj(1) !  sizes(1:3)
        IF ( vname == 'W' .or. vname == 'EZ' ) THEN
          chunksizes(3) =  mzp1t !  sizes(1:3)
        ELSE
          chunksizes(3) =  mzt !  sizes(1:3)
        ENDIF
        chunksizes(4)   = 1
        
        status = NF90_DEF_VAR_DEFLATE(ncid1, varid2, shuffle, deflate, deflate_level)
        IF (status /= NF90_NOERR) write(0,*)'DEFINE_NCDF_VAR:  Error setting deflate for: ',vname

        status = NF90_DEF_VAR_CHUNKING(ncid1, varid2, NF90_CHUNKED, chunksizes)
        IF (status /= NF90_NOERR) print *,'DEFINE_NCDF_VAR:  Error setting chunking for: ', vname
        deallocate( chunksizes )
       ENDIF
#endif

! grab and write variable attributes
     DO k = 1,numAtts
     write(0,*) 'Put attribute ',k

       status = NF90_INQ_ATTNAME(ncid(1),vid,k,AttName)
       status = NF90_INQUIRE_ATTRIBUTE(ncid(1),vid,trim(AttName),atype)
       IF (atype == 2) THEN		! Attribute value is type character
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValC)
         status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),trim(AttValC))
         IF (trim(AttName) == 'type') vdim = AttValC(1:1)
       ENDIF
       IF (atype == 4) THEN		! Attribute value is an integer
         status = NF90_GET_ATT(ncid(1),vid,trim(AttName),AttValI)
         status = NF90_PUT_ATT(ncid1,varid2,trim(AttName),AttValI)
         IF (trim(AttName) == 'istag') istag = AttValI
         IF (trim(AttName) == 'jstag') jstag = AttValI
         IF (trim(AttName) == 'kstag') kstag = AttValI
       ENDIF

     ENDDO
     write(0,*) 'call NF90_ENDDEF'
     status = NF90_ENDDEF(ncid1)
     write(0,*) 'done NF90_ENDDEF'
     
     ENDIF

     IF (istag == 1) xsize = mxp1t
     IF (istag == 0) xsize = mxt
     IF (jstag == 1) ysize = myp1t
     IF (jstag == 0) ysize = myt
     IF (kstag == 1) zsize = mzp1t
     IF (kstag == 0) zsize = mzt
     allocate(var3D(xsize,ysize,zsize))

      write(0,*) 'numDims,xsize,ysize,zsize = ',numDims,xsize,ysize,zsize
!     IF (numDims == 4) THEN
       tbeg = ntbeg
       tend = ntend
!     ELSE
!       tbeg = 1
!       tend = 1
!     ENDIF

     DO k = tbeg,tend
       kk = k-tbeg+ntstart
       var3D(:,:,:) = 0.0

       jstop=0
       DO j = jbeg,ntiles_y+jbeg-1
         istop=0
         jstart = jstop+1
         DO i = ibeg,ntiles_x+ibeg-1
           jj = (j-1)*ntilef_x+i
           IF (istag == 1) mxadd = mxp1(i)
           IF (istag == 0) mxadd = mx(i)
           istart = istop+1
           istop  = istart+mxadd-1
           IF (jstag == 1) myadd = myp1(jj)
           IF (jstag == 0) myadd = my(jj)
           jstop  = jstart+myadd-1
           status = NF90_INQ_VARID(ncid(jj),vname,varid1)
           
!           IF ( j == jbeg .and. i == ibeg ) THEN
!           allocate( tem3d(mxadd,myadd,zsize) )
!            status = NF90_GET_VAR(ncid(jj),varid1,tem3d,(/1,1,1,k/))
           
!           deallocate (tem3d)
!           ENDIF
           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),start=(/1,1,1,k/))
!           status = NF90_GET_VAR(ncid(jj),varid1,var3D(istart:istop,jstart:jstop,:),(/1,1,1,k/),(/istop-istart+1,jstop-jstart+1,kstop,1/))
         ENDDO
       ENDDO
!       print *,'xsize =',xsize,'  ysize =',ysize
!       print *,'istop =',istop,'  jstop =',jstop
       status = NF90_PUT_VAR(ncid1,varid2,var3D,(/1,1,1,kk/)) 
       IF (status /= nf90_noerr) THEN
         print *, 'PUT_VAR  ',trim(nf90_strerror(status))
       ENDIF
     vmax = 0.0
     vmin = 0.0
     DO k1 = 1,zsize
       DO j1 = 1,ysize
         DO i1 = 1,xsize
          vmax = Max( vmax, var3D(i1,j1,k1) )
          vmin = Min( vmin, var3D(i1,j1,k1) )
         ENDDO
       ENDDO
     ENDDO
     
     write(0,*) 'k,kk, variable max/min = ',k,kk,vmax,vmin
     ENDDO
     
     
     deallocate(var3D)
     deallocate(DimIDV)
         
   RETURN
   END SUBROUTINE

999 END PROGRAM COMBINE_V3

!-----------------------------------------------------------------------
! SUBROUTINE DISP_ERROR
! converts error message to text form and prints it
!-----------------------------------------------------------------------
  SUBROUTINE DISP_ERR(status,die_on_error)

   use netcdf

   implicit none
!   include '/usr/local/netcdf3-32/include/netcdf.inc'

   integer, intent (in) :: status
   logical, intent (in) :: die_on_error
  
   IF (status /= nf90_noerr) THEN
     print *, trim(nf90_strerror(status))
    IF( die_on_error) THEN
     stop "STATUS returned an error, program stopped"
    ENDIF
   ENDIF

  END SUBROUTINE DISP_ERR
