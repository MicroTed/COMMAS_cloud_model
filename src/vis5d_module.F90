MODULE VIS5D_MODULE

   implicit none

   include 'v5df90.h'

   logical             :: lncwrite                = .false.      ! flag for history dump
   logical             :: lrst                    = .false.      ! flag for restart
   logical             :: lprt                    = .true.       ! flag for stats print
   integer             :: v5dmaxtimes = MAXTIMES
!   logical             :: lstt                    = .true.       ! flag for stats print
   integer             :: ntstart                                ! model start time
   logical             :: lchgratn = .false. 
   logical             :: lchgrati = .false.

   character(LEN = 120) :: v5dfilename
   character(LEN = 120) :: onedfilename
   integer             :: lenrunname                             ! length of runname
!   integer             :: tvis5d = 999999
!   integer             :: tvis5dstart             = 0
   logical             :: lv5dwrite               = .false.      ! .true. when it's time to write to vis5d file
   logical             :: lv5dcreate              = .false.      ! .true. when vis5d file has been created
   integer             :: itv5d                   = 1            ! current time level for vis5d file
   character(LEN = 120) :: inv5dfile               = 'invis5d'
   integer, parameter  :: nv5dfields              = MAXVARS

   integer             :: iv5dfields(nv5dfields)  = -1 
   character(LEN=10)   :: v5dfields(nv5dfields)   = '          '
   integer             :: iv5dwritten(nv5dfields) = 0

   integer             :: numvars                 = IMISSING
   integer             :: numtimesv5d             = IMISSING

   integer, parameter  :: ntitle = 1000

   character(LEN=10)   :: title(2,ntitle) = '          '
   character(LEN=20)   :: units(ntitle) = '                    '

   character(LEN=20)   :: dunit

   integer             :: compressmode            = 1
   integer             :: projection              = 0  ! generic units
   integer             :: vertical                = 2
   integer             :: nzdbz                   = 0

   integer             :: itimes(MAXTIMES)        = IMISSING
   integer             :: v5dtimes(MAXTIMES)      = 0
   integer             :: idates(MAXTIMES)        = IMISSING
   integer             :: idx(MAXVARS)            = -1
   integer             :: nl(MAXVARS)             = IMISSING
   
   logical             :: idoradar                = .false.
   
   integer             :: vis5dstridex = 1
   integer             :: vis5dstridez = 1
   
   integer             :: ixb_v5d = -1
   integer             :: ixe_v5d = -1
   integer             :: jyb_v5d = -1
   integer             :: jye_v5d = -1
   integer             :: kzb_v5d = -1
   integer             :: kze_v5d = -1
   
   integer             :: nx_v5d = -1
   integer             :: ny_v5d = -1
   integer             :: nz_v5d = -1
   
   real                :: alphamax_v5d = 15., alphamin_v5d = 0.

   
   real, allocatable   :: flshn(:,:,:), flshp(:,:,:), flshi(:,:,:)

CONTAINS


! ############################################################################

  SUBROUTINE Vis5DInit(gd,                 &
                       gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                       gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                       gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                       v5dstridex, v5dstridez, mpisep, iflag   )                 

   USE GRID_MODULE
   USE PARAM_MODULE
   USE COMMASMPI_MODULE
   implicit none

! Passed Variables:
   TYPE(GRID), target          :: gd
   TYPE(VARIABLE)              :: gx(4), gy(4), gz(4)
   integer, optional           :: iflag
   integer                    :: v5dstridex, v5dstridez

! Local Variables:

   TYPE(ATTRIBUTE), pointer    :: v5dflds
   character(LEN= 10)          :: v5dfield_tmp
   
   integer      ::   v5dcreate
   
   integer      ::    nr0,nc,n,it,i,j
   integer      ::    idatime, index
   integer      ::    nx,ny,nxe,nye
   integer      ::    kz,nz,nze
   integer      ::    varnum = 1

   character*10 ::       varname(MAXVARS) = '          '
!   integer      ::       projection
   real         ::       proj_args(MAXPROJARGS)
   real         ::       vert_args(MAXLEVELS)
   real         ::       hh,mm,ss
   integer      ::       time, tstop
   integer      ::       ibeg,slen
   logical      ::       hstretch
   real         ::       dx, dy
   integer      ::       tvis5d
   real         ::       tmp
   integer      ::       itmp,itmp2
   integer      ::       thistory

! MPI LOCAL VARIABLES

   logical :: mpisep
   integer :: ixb,ixe

! --------------------------------------------------------------------

    CALL GET_VARIABLE (gd, 'NX', nx)
    CALL GET_VARIABLE (gd, 'NY', ny)
    CALL GET_VARIABLE (gd, 'NZ', nz)
    CALL GET_VARIABLE (gd, 'NXEND', nxe)
    CALL GET_VARIABLE (gd, 'NYEND', nye)
    CALL GET_VARIABLE (gd, 'NZEND', nze)

    if ( .not. mpisep ) then
     nx = ncxe
     ny = ncye
     nz = ncze
    endif
    
    nx_v5d = nx
    ny_v5d = ny
    nz_v5d = nz
    
    IF ( ixb_v5d > 0 .and. ixe_v5d > 0 ) THEN
       ixe_v5d = Min( ixe_v5d, nx )
       nx_v5d = ixe_v5d - ixb_v5d + 1
    ELSE
       ixb_v5d = 1
       ixe_v5d = nx_v5d
    ENDIF
    
    IF ( jyb_v5d > 0 .and. jye_v5d > 0 ) THEN
       jye_v5d = Min( jye_v5d, ny )
       ny_v5d = jye_v5d - jyb_v5d + 1
    ELSE
       jyb_v5d = 1
       jye_v5d = ny_v5d
    ENDIF

    IF ( kzb_v5d > 0 .and. kze_v5d > 0 ) THEN
       kze_v5d = Min( kze_v5d, nx )
       nz_v5d = kze_v5d - kzb_v5d + 1
    ELSE
       kzb_v5d = 1
       kze_v5d = nz_v5d
    ENDIF
    

    CALL GET_VARIABLE (gd, 'DX', dx)
    CALL GET_VARIABLE (gd, 'DY', dy)

    CALL GET_VARIABLE(gd, 'TIME',          time)
    CALL GET_VARIABLE(gd, 'TIME_STOP',    tstop)
    CALL GET_VARIABLE(gd, 'TVIS5D',      tvis5d)
    CALL GET_VARIABLE(gd, 'THISTORY',  thistory)
    CALL GET_ATTRIBUTE(gd, 'V5DFIELDS', v5dflds)

    vis5dstridex = v5dstridex
    vis5dstridez = v5dstridez

!    print*, 'Vis5dInit: luno = ',luno
    write(luno,*) 'Vis5DInit: time = ', tvis5d + time
!    write(*,*) 'Vis5DInit: time = ', time

        nc=(nx_v5d-1)/v5dstridex
        nl(1)=(nz_v5d-1)/v5dstridez
!
!  Read the v5dflds string for the Vis5D output variables
!
    numvars = 0
    DO n = 1,desc_length
     IF(( v5dflds%str(n:n) .EQ. '/' .or. v5dflds%str(n:n) .EQ. ',')    &
  &       .and. n .gt. 1) numvars = numvars + 1
    ENDDO

!    print*, 'Vis5DInit: numvars = ',numvars

    varnum = 0
    index  = 0
    v5dfield_tmp = '         '

    DO n = 1,desc_length

      IF ( n .eq. 1 .and. ( v5dflds%str(n:n) .EQ. '/' .or.   &
      v5dflds%str(n:n) .EQ. ',' .or. v5dflds%str(n:n) .EQ. ' ') ) CYCLE

      IF(v5dflds%str(n:n) .NE. '/') THEN
        index = index + 1
        v5dfield_tmp(index:index) = v5dflds%str(n:n)
      ELSE
        index = 0
        varnum = varnum + 1
        v5dfields(varnum) = v5dfield_tmp
        v5dfield_tmp = '         '
      ENDIF

      IF( varnum .eq. numvars ) EXIT

    ENDDO

    


! END OF ATTRIBUTE PARSE
!
!    v5dfields(:) = an array of character(LEN=10)
!    numvars      = number of valid fields in v5dfields
!

! Now we will search through the list of available variables for the associated index
!  to get the units, etc.

      CALL settitles(1000,title,units)


      varnum = 0                      ! nv5dfields
      nzdbz  = 0
         
      DO kz=1,nv5dfields

!      print*, 'V5DOUT: v5dfields(',kz,') = ',v5dfields(kz)
          
        IF ( v5dfields(kz) .eq. '   ' ) EXIT

        DO n=1,ntitle
          IF ( v5dfields(kz) .eq. title(1,n) .or. v5dfields(kz) .eq. title(2,n)  ) THEN
            iv5dfields(kz) = n
            EXIT
          ENDIF
        ENDDO
          
        idx(kz) = iv5dfields(kz)
        IF ( idx(kz) .le. 0 ) EXIT
          
        varnum = varnum + 1
        varname(varnum) = title( 1,idx(kz) ) 

        write(luno,*) 'Vis5DInit: adding variable ',varname(varnum)
!        write(*,*) 'Vis5DInit: adding variable ',varname(varnum)

        IF ( idx(kz) .eq. 51 .or. idx(kz) .eq. 22 ) THEN
          idoradar = .true.
          IF ( idx(kz) .eq. 51 ) nzdbz = Max(nzdbz,1)
          IF ( idx(kz) .eq. 22 ) nzdbz = Max(nzdbz,nz_v5d-1)
        ENDIF
        IF ( idx(kz) .eq. 50 .or. idx(kz) .eq. 51 .or. idx(kz) .eq. 41  .or.         &
             idx(kz) .eq. 28 .or. idx(kz) .eq. 29 .or. idx(kz) .eq. 30 .or. idx(kz) .eq. 31 .or.           &
             ( idx(kz) .ge. 58 .and. idx(kz) .le. 90 ) ) THEN
          nl(kz) = 1
        ELSE 
          nl(kz) = (nz_v5d-1)/v5dstridez

        ENDIF
        
        IF ( thistory .ne. tvis5d .and. .not. present( iflag ) ) THEN
          IF ( idx(kz) .eq. 18 ) THEN  ! FLSHI
           allocate( flshi(nx,ny,nz) )
           flshi(:,:,:) = 0.0
          ENDIF
        
          IF ( idx(kz) .eq. 19 ) THEN  ! FLSHP
           allocate( flshp(nx,ny,nz) )
           flshp(:,:,:) = 0.0
          ENDIF

          IF ( idx(kz) .eq. 20 ) THEN  ! FLSHN
           allocate( flshn(nx,ny,nz) )
           flshn(:,:,:) = 0.0
          ENDIF
        
        ENDIF
        
        IF ( idx(kz) .eq. 856 ) THEN  ! cghis_net
           lchgratn = .true.
        ENDIF
        IF ( idx(kz) .eq. 857 ) THEN  ! CGHW_NET
           lchgrati = .true.
        ENDIF

      ENDDO
        
      IF ( varnum .eq. 0 ) THEN
        write(0,*) 'ERROR in Vis5DInit: varnum = 0'
        write(0,*) 'Either set v5dfields in namelist or turn off vis5d output'
        STOP
      ENDIF

      DO kz=1,(nz_v5d - 1)/v5dstridez
        vert_args(kz) = gz(1)%flt1d((kz-1)*v5dstridez + kzb_v5d )/1000.           !  gz(1,1,kz,mzdist)/1000.0
!        print*,'k,gz(1),gz(2) = ',kz, gz(1)%flt1d(kz)/1000., gz(2)%flt1d(kz)/1000.
      ENDDO
        
! Now check whether the grid is stretched horizontally or not
! kinda backwards: hstretch = true means NOT stretched
      hstretch = .true.
     nr0 = Max(2,(ny_v5d-1)/v5dstridex)

     DO j = 1,ny-1
       hstretch = hstretch .and. gy(3)%flt1d(1) .eq. gy(3)%flt1d(j) 
!       print*, 'j, Y: ',j,proj_args(j)
     ENDDO

     DO i = 1,nx-1
       hstretch = hstretch .and. gx(3)%flt1d(1) .eq. gx(3)%flt1d(i) 
     ENDDO

     
     IF ( hstretch ) THEN 
      projection = 0
!
      proj_args(1) =  (ny/v5dstridex - 1) * dy * v5dstridex  ! ymax
      proj_args(2) =  (nx/v5dstridex - 1) * dx * v5dstridex  ! xmin
      proj_args(3) =  dy*v5dstridex       ! yinc
      proj_args(4) =  dx*v5dstridex       ! xinc

     ELSE
     

      IF ( ( gy(3)%flt1d(1 + jyb_v5d - 1) == gy(3)%flt1d(nr0 + jyb_v5d - 1) ) .and. &
           ( gx(3)%flt1d(1 + ixb_v5d - 1) == gx(3)%flt1d(ixe + ixb_v5d - 1) )) THEN
 
       projection = 0
       proj_args(1) =  (nr0/v5dstridex - 1) * v5dstridex/gy(3)%flt1d(1 + jyb_v5d - 1)  ! ymax
       proj_args(2) =  (ixe/v5dstridex - 1) * v5dstridex/gx(3)%flt1d(1 + ixb_v5d - 1)  ! xmin
       proj_args(3) =  v5dstridex/gy(3)%flt1d(1 + jyb_v5d - 1)    ! yinc
       proj_args(4) =  v5dstridex/gx(3)%flt1d(1 + ixb_v5d - 1)       ! xinc

       
      ELSE
       projection = -1
       write(luno,*) 'Vis5DInit: Set projection args for horizontal stretching'

       IF ( v5dstridex == 1 ) THEN
         DO j = 1,nr0
           proj_args(j) = gy(1)%flt1d(j + jyb_v5d - 1)
!           print*, 'j, Y: ',j,proj_args(j)
         ENDDO

         ixb = 1
         ixe = nc ! nx-1
!         ixe = itile
!         IF ( ixend .eq. nxend ) ixe = nx-1
         DO i = 1,ixe ! nx-1
            proj_args(i+nr0) = gx(1)%flt1d(i + ixb_v5d - 1)
!           print*, 'i, X: ',i,proj_args(i+nr0)
         ENDDO

       ELSE
         DO j = v5dstridex, nr0*v5dstridex, v5dstridex
           proj_args(j/v5dstridex) = gy(1)%flt1d(j + jyb_v5d - 1)
!           print*, 'j,j/v5dstridex, Y: ',j,j/v5dstridex,proj_args(j/v5dstridex)
         ENDDO

         DO i = v5dstridex, nc*v5dstridex, v5dstridex
           proj_args(nr0 + i/v5dstridex ) = gx(1)%flt1d(i + ixb_v5d - 1)
!          print*, 'i,i/v5dstridex, X: ',i,i/v5dstridex,proj_args(nr0 + i/v5dstridex)
         ENDDO
       ENDIF
      ENDIF

     ENDIF
     
     IF ( .not. present( iflag ) .or. time .eq. 0 ) THEN
!        idatime = tvis5d + time
        idatime = time
     ELSE
        idatime = time
     ENDIF
        hh      = idatime/3600
        idatime = idatime - 3600*hh
        mm      = idatime/60
        ss      = idatime - 60*mm
     
        it = itv5d
        
        itimes(it) = 10000*hh+100*mm+ss
        write(luno,*) 'Vis5DInit itimes(it) = ',it, itimes(it)
        idates(it) = 2006001

          write(luno,'(a,a)') 'Vis5DInit: Writing Vis5D file: ',v5dfilename

!       IF ( it .eq. 0 ) THEN
        IF ( .not. lv5dcreate ) THEN

!  will need to check if this is a restart and whether to open the 
!  previous file for further writing.
!
        lv5dcreate = .true.
        numtimesv5d = 1
        n = v5dcreate( v5dfilename, numtimesv5d, varnum, nr0, nc, nl, &
                       varname, itimes, idates, compressmode,            &
                       projection, proj_args, vertical, vert_args )

          write(luno,'(a,i3)') 'Creating the Vis5D file, n= ',n
          IF ( n .eq. 0 ) THEN
           write(0,*) 'Error creating file!'
           write(0,*) 'numtimesv5d =',numtimesv5d
           write(0,*) 'varnum = ',varnum
           write(0,'(a,a)') 'varname(1) = ',varname(1)
          ENDIF
        ENDIF
           
! set the max number of time levels before to keep the file under 2 GB
! tmp is the size (in bytes) per time level
   tmp = nr0*(nc)*(nz_v5d-1)/v5dstridez*compressmode*numvars
! itmp is the max number of levels to stay under 2 GB
   itmp = 2*Int(2**30/tmp) - 1
   itmp2 = Int( (tstop - time)/tvis5d )
   
!   print*, 'tmp,itmp,itmp2 = ',tmp,itmp,itmp2
   
   IF ( itmp2 .gt. itmp .or. itmp2 .gt. v5dmaxtimes ) THEN
     write(0,*) 'WARNING!! TOO MANY TIMELEVELS (OR TOO MUCH DATA) FOR 32-bit VIS5D!'
     write(0,*) 'File will not work correctly with 32-bit addressing in unpatched vis5d'
     write(0,*) 'REQUESTED = ',itmp2
     write(0,*) 'MAXIMUM   = ', Min(itmp, v5dmaxtimes)
     IF ( itmp2 .lt. itmp ) THEN
       write(0,*) 'Could increase MAXTIMES and recompile COMMAS and Vis5D'
     STOP
     ELSE
#ifndef V5D64
       write(0,*) 'FILE WILL EXCEED 2 GB!'
       write(0,*) 'OPTIONS: 1. Increase tvis5d'
       write(0,*) '         2. Decrease number of vis5d variables'
       write(0,*) '         3. Set v5dstridex or v5dstridez > 1 to decimate data'
       write(0,*) '            (v5dstridex controls both x and y)'
       write(0,*) '            (v5dstridex/v5dstridez set in ATTRIBUTES namelist)'
       write(0,*) '         4. Use 64-bit addressing v5d.c (still a 32-bit compile)'
     STOP
#endif
     ENDIF
   ENDIF

   v5dmaxtimes = Min( v5dmaxtimes, itmp )
   
   RETURN
   END SUBROUTINE Vis5DInit


! ############################################################################

  SUBROUTINE Vis5DOPEN(gd)

   USE GRID_MODULE
   USE PARAM_MODULE
   
   implicit none

   TYPE(GRID), target :: gd
   
!   integer v5dopen
   
   integer nr0,nc,n,it
   integer kz,nx,ny,nz
   integer :: varnum = 1

   character*10 :: varname(MAXVARS) = '          '

   integer :: istuff(100) = 0

!   integer    projection
   real       proj_args(MAXPROJARGS)
   integer    vertical
   real       vert_args(MAXLEVELS)
   real       hh,mm,ss
   integer    time
   integer    ibeg,slen
   integer    tvis5d
   integer    thistory

!-------------------------------------------------------------------------------
   
    CALL GET_VARIABLE (gd, 'NX', nx)
    CALL GET_VARIABLE (gd, 'NY', ny)
    CALL GET_VARIABLE (gd, 'NZ', nz)

    CALL GET_VARIABLE(gd, 'TIME',       time)
    CALL GET_VARIABLE(gd, 'TVIS5D',      tvis5d)
    CALL GET_VARIABLE(gd, 'THISTORY',  thistory)

    write(luno,*) 'Vis5DOPEN: time = ', time

! Only call this routine from commas, so ixb_v5d < 0 always, i.e., safe to set limits to size of domain/tile

    nx_v5d = nx
    ny_v5d = ny
    nz_v5d = nz
    
    ixb_v5d = 1
    ixe_v5d = nx_v5d
    
    jyb_v5d = 1
    jye_v5d = ny_v5d

    kzb_v5d = 1
    kze_v5d = nz_v5d
    

   nr0 = 0
   
!   print*, 'v5dfilename = ',v5dfilename
!   print*, 'compressmode (before) = ',compressmode
!   print*, 'varnum (before) = ',varnum
!   print*, 'itimes (before) = ',itimes(1)
    call myv5dopen( v5dfilename, numtimesv5d, numvars, nr0, nc, nl, &
                 varname, itimes, idates, compressmode,            &
                 projection, proj_args, vertical, vert_args, istuff )
   
!   print*, 'compressmode (after) = ',compressmode
!   print*, 'varnum,nr0 (after) = ',numvars,nr0
!   print*, 'itimes (after) = ',itimes(1),itimes(10)
!   print*, 'istuff(1) = ',istuff(1)
   
   lv5dcreate = .true.
   itv5d      = numtimesv5d + 1

!   vis5dstridex = 1
   IF ( nc .ne. nx-1 ) vis5dstridex = NInt(float(nx-1)/nc)
   write(luno,*) 'Vis5DOPEN: vis5dstridex = ',vis5dstridex

   IF ( nl(1) .lt. nz-1 ) vis5dstridez = NInt(float(nz-1)/nl(1))
   write(luno,*) 'Vis5DOPEN: vis5dstridez = ',vis5dstridez

      CALL settitles(1000,title,units)


      varnum = 0                      ! nv5dfields
      nzdbz  = 0
         
      DO kz=1,numvars

 !     print*, 'Vis5DOPEN: varname,nl(',kz,') = ',varname(kz),nl(kz)
          

        
        CALL STRING_LIMITS(varname(kz) , ibeg, slen)
        v5dfields(kz)(ibeg:slen) = varname(kz)(ibeg:slen)
        IF ( slen .lt. 10 ) THEN
          v5dfields(kz)(slen+1:10) = ' '
        ENDIF 
        DO n=1,ntitle
          IF ( v5dfields(kz) .eq. title(1,n) .or. v5dfields(kz) .eq. title(2,n)) THEN
            iv5dfields(kz) = n
!            print*, 'Found variable ', v5dfields(kz), ' as ', title(1,n)
            EXIT
          ENDIF
        ENDDO
        
        IF ( n .ge. ntitle ) THEN
         print*,'Warning, did not find variable ',v5dfields(kz), ', n = ',n,', slen = ',slen
        ENDIF
          
        idx(kz) = iv5dfields(kz)
          
        varnum = varnum + 1
!        varname(varnum) = title( idx(kz) ) 

       write(luno,*) 'Vis5DOPEN: adding variable ',varname(varnum)(ibeg:slen), ' as ', v5dfields(kz)
!       write(*,*) 'Vis5DOPEN: adding variable ',v5dfields(kz),', ',units(idx(kz)),', ',idx(kz)

        IF ( idx(kz) .eq. 51 .or. idx(kz) .eq. 22 ) THEN
          idoradar = .true.
        ENDIF
        IF ( idx(kz) .eq. 856 ) THEN  ! cghis_net
           lchgratn = .true.
        ENDIF
        IF ( idx(kz) .eq. 857 ) THEN  ! CGHW_NET
           lchgrati = .true.
        ENDIF

        IF ( thistory .ne. tvis5d  ) THEN
          IF ( idx(kz) .eq. 18 ) THEN  ! FLSHI
           allocate( flshi(nx,ny,nz) )
           flshi(:,:,:) = 0.0
          ENDIF
        
          IF ( idx(kz) .eq. 19 ) THEN  ! FLSHP
           allocate( flshp(nx,ny,nz) )
           flshp(:,:,:) = 0.0
          ENDIF

          IF ( idx(kz) .eq. 20 ) THEN  ! FLSHN
           allocate( flshn(nx,ny,nz) )
           flshn(:,:,:) = 0.0
          ENDIF
        
        ENDIF

      ENDDO
   
    nzdbz = nz - 1
   
   DO it = 1,numtimesv5d
     hh = itimes(it)/10000
     mm = (itimes(it) - hh*10000)/100
     ss = itimes(it) - hh*10000 - mm*100
     v5dtimes(it) = hh*3600 + mm*60 + ss
     IF ( v5dtimes(it) .ge. time ) THEN
       itv5d = it + 1
       EXIT
     ENDIF
!     print*,'Vis5DOPEN: it,times = ',it,v5dtimes(it)
   ENDDO
!   STOP
   
   RETURN
   END SUBROUTINE Vis5DOPEN

!-----------------------------------------------------------------------------
!
!
! SUBROUTINE V5DOUT puts model data into a vis5d file.
!
!
!-----------------------------------------------------------------------------

  SUBROUTINE V5DOUT(gd,                 &
                    u,  uinit ,         &            ! U,  UINIT
                    v,  vinit ,         &            ! V,  VINIT
                    w,  winit ,         &            ! W,  WINIT
                    pi, piinit,         &            ! PI, PIINIT
                    km, kminit,         &            ! KM, KINIT
                    s,  sinit0,         &            ! S,  SINIT
                    precip,             &            ! PRECIP
                    gx,                 &            ! XCNTR, XEDGE, DXC, DXE
                    gy,                 &            ! YCNTR, YEDGE, DYC, DYE
                    gz,                 &            ! ZCNTR, ZEDGE, DZC, DZE
                    dt,                 &            ! DT
                    ugrid, vgrid,       &            ! GRID MOTION
                    nx, ny, nz, na,     &            ! NX,NY,NZ,NS
                    st,                 &            ! scalar array
                    dn, tt7, tt0, pn,   &
                    dbz, vzf, wz, elec, xtra, onedoutput,  &
                    recalcdbz)

   USE GRID_MODULE
   USE PARAM_MODULE,    only: pii,bcx,bcy,ng,isfcphys,luno
   USE MICRO_MODULE
   USE INDEX_MODULE
   USE COMMASMPI_MODULE
!   USE takcommon
   use takcommon,       only: lmax,lfmax,kfmax,iimax,kimax,lr100,ls250
   use comm1,           only: rhof
   use comm2,           only: sr,sx,srf,fw,vw,bew,cdw,szr,szf,szi      &
     &                            ,sri,sfi,di,sxf,sxi,sh,ff,vf,fi,vi       
!-----------------------------------------------------------------------------
! VARIABLE DECLARATIONS
   
   implicit none

   TYPE(GRID), target :: gd

   integer :: nx, ny, nz, na
   integer :: nsmall    
   real    :: dt, ugrid, vgrid     

   TYPE(VARIABLE) :: u, uinit
   TYPE(VARIABLE) :: v, vinit
   TYPE(VARIABLE) :: w, winit
   TYPE(VARIABLE) :: pi, piinit
   TYPE(VARIABLE) :: km, kminit
   TYPE(VARIABLE) :: s(na), sinit0(2)     ! We are going to try and create arrays of the scalar variables needed here
   TYPE(VARIABLE) :: precip(nprecip+neelec2d)
   TYPE(VARIABLE) :: gx(4), gy(4), gz(4)
   TYPE(VARIABLE) :: dbz, vzf, wz
   TYPE(VARIABLE) :: xtra(nxtra)
   TYPE(VARIABLE) :: elec(neelec)
   
   integer, INTENT(IN) :: onedoutput
   integer, optional :: recalcdbz

!   real :: preciptmp(-ng+1:nx+ng,-ng+1:ny+ng,4)
      
   real :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na)

   real :: dn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: pn (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 

   real :: tt7 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) 
   real :: tt0 (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

   real, allocatable :: tem1(:,:,:),tem2(:,:,:)
   real, allocatable ::  sinit(:,:)
      
!-----------------------------------------------------------------------------
! LOCAL variable declarations


   real    :: pb(nz)
   real    :: db(nz)

   integer :: is
   integer :: js
   integer :: ks
   
   real, parameter :: segmx = 20.0

   real    :: smin, smax

   integer v5dcreate
   integer v5dwrite
   integer v5dwriteappend
   integer v5dmcfile
   integer v5dclose
   integer myv5dupdate

   real :: gtx(nx), gty(ny), gtz(nz)
   real :: gxt(nx,4), gyt(ny,4), gzt(nz,4)

   integer  time
      
   integer nstep1
   real start,tstat
      
   integer i,j,k,l, ix,jy,kz, n, ia
   integer i1,j1,k1
   real    fac, fac1d
      
   real dx,dy,dz
   real x,tx
   integer nstep, nstop, index


   integer idatime,hh,mm,ss
      
   character*100 line
   integer istat1

   logical ifv5d
   integer nzdbz
   data ifv5d/.true./
   integer iv5dinterval,iv5dstart
   data iv5dinterval / 2000 /
      
   save ifv5d,iv5dinterval,nzdbz
   integer nr, nc !, nl(MAXVARS)
   integer nr0, nr1, nr2
   integer varnum, it, iv
   integer itimes0,idates0

   TYPE(ATTRIBUTE), pointer    :: v5dflds
   character(LEN= 10) :: v5dfield_tmp

   real, parameter :: eperao  = 8.8592e-12

! Initialize the variables to missing values

    data nr,nc / IMISSING, IMISSING /
      
    data it / 0 /
    data iv5dstart / 1 /
      
    save varnum,it,iv5dstart

    integer initcond
    save initcond
    data initcond/0/
      
    integer len
    TYPE(ATTRIBUTE), pointer    :: microphys
    character(LEN = 15) :: micro
    logical    lice
    
    real x_sw_loc, y_sw_loc
    
    real, parameter  ::  ec = 1.602e-19 ! fundamental unit of charge

    real, parameter  :: ccimx = 200.0
    
    real :: hwdn, tmp, tmpg, tmpn, xdia, cno, tmpmx, voltot, diam
    real :: tmpmax, qmax
    
    real :: hwdnl(nch),hwvol(nch)
    
    real :: pres,qwv,theta1,thetae
    
    real,allocatable :: data_1d(:,:) ! (MAXVARS,nz)
    integer, save :: out1d_flag = 0

    
    integer :: ibc, jbc
    
    integer :: icldtop
    integer :: mask(nx,ny)
    
!    real, parameter :: pii = 3.14159265359
    real :: cwch

        real chw,qr,z,alp,alpha1, rdi,pi1,vr,g1,zx,ze,nrx, sum, area, areaw
        integer ii

      real qxmin(lc:lhab)

  integer :: tcanp
  integer :: wcanp
  integer :: qav
  integer :: veg
  integer :: stype
  integer :: tsrfc
  integer :: wsfc
  integer :: tsoil
  integer :: wsoil
  integer :: vlai
  integer :: albedo
  integer :: rough
  integer :: qv
  integer :: eflx
  integer :: fflx
  integer :: uflx
  integer :: vflx
  integer :: tflx
  integer :: qflx
  integer :: radsw
  integer :: radlw

      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      real dv, dax, day, daz, phiw, phie, phis, phin, phid, phiu

      double precision :: ztmp, ztmpr, ztmph, ztmphl, ztmpi
      real dbzmax,dbzmin
      parameter ( dbzmin = -10.0 )
      integer :: ki,il


   real :: rdamelt(nch)

      real pqs,temg,qvs,qis,temq, fsw

      integer ltemq
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati
!      parameter (fqsat=0.002,fqsati=1./fqsat)
!      parameter (fqsat=0.01,fqsati=1./fqsat)
!      real, save :: tabqvs(nqsat)=0.0,tabqis(nqsat)=0.0
      integer, save :: itabflag = 0


          integer :: lf75 ! 75 micron radius; use as dividing line between ice crystal and graupel/frozen drop

      real gsnow1, gsnow53, gsnow73, gamma
      
      real alphamaxtmp, alphamintmp


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

   logical :: mpisep = .false.

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

   logical :: debug_mpi = .false.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc Begin Execute
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    CALL GET_VARIABLE (gd, 'DX', dx)
    CALL GET_VARIABLE (gd, 'DY', dy)
    CALL GET_VARIABLE (gd, 'DZ', dz)

    CALL GET_VARIABLE(gd, 'XG_POS',   x_sw_loc)
    CALL GET_VARIABLE(gd, 'YG_POS',   y_sw_loc)

    CALL GET_ATTRIBUTE (gd, 'MICROPHYS', microphys)
    micro(:) = microphys%str(:)

    allocate ( sinit(nz,na) )
    allocate ( tem1(Max(2,ny_v5d-1),nx_v5d-1,nz_v5d-1) )
    allocate ( tem2(Max(2,(ny_v5d-1)/vis5dstridex),(nx_v5d-1)/vis5dstridex,(nz_v5d-1)/vis5dstridez) )

    sinit(:,:) = 0.0
      
    DO n = 1,2
     DO k = 1,nz_v5d
      sinit(k,n) = sinit0(n)%flt1d(k)
     ENDDO
    ENDDO
    
    icldtop = -1
    mask(:,:) = 1
          pi1 = 4.0*atan(1.0)


      lice = .true.
      IF ( li .lt. 1 ) lice = .false.

        gsnow1 = gamma(snu + 1.0)
        gsnow53 = gamma(snu + 5./3.)
        gsnow73 = gamma(snu + 7./3.)


!-------------------------------------------------------------------------------
! Now need to get ATTRIBUTE:  V5DFIELDS FROM DATA BASE AND PARSE IT

    
      varnum = numvars

      IF ( micro(1:5) .eq. 'ICE10' ) THEN
         call setqxmin(qxmin)
       ELSEIF ( micro(1:8) .eq. 'WARMZIEG' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 
         call setqxminz(qxmin)
       ELSEIF ( micro(1:4) .eq. 'ZIEG' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 
         call setqxminz(qxmin)
       ELSEIF ( micro(1:3) .eq. 'ZVD' ) THEN !  na .ge. 14 .and. ipconc .ge. 3 ) THEN 
         call setqxminz(qxmin)
       ELSEIF ( micro(1:3) .eq. 'TAK' ) THEN
         call setqxminz(qxmin)
         IF ( sr(1) == 0.d0 ) THEN
           call takinitbin
           write(0,*) 'call takinitbin, sx(1) = ',sx(1),sxf(1,1),sxf(1,2),srf(2,lfmax/2),lmax,lfmax
         ENDIF
        ! options for setting the size cut-off for treating small graupel/frozen drops as ice crystals for collisional charge separation
          lf75 = 1
          IF ( taknicfd <= 1 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0075 ) lf75 = l
          enddo
          ELSEIF ( taknicfd == 2 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0025 ) lf75 = l
          enddo
          ELSEIF ( taknicfd == 3 ) THEN
          do l=1,lfmax
            IF ( srf(l,2) <= 0.0125 ) lf75 = l
          enddo
          ENDIF

       ELSEIF ( micro(1:3) .eq. 'HCM' ) THEN
         IF ( hm(1) == 0.0 ) THEN
         write(0,*) 'v5dout: set up hm'
          DO l = 1,nch
            hm(l) = hmmin*exp(3.0*(l-1)/hjo) 
!            hmdn(l) = 900.
!            hmvol(l) = hm(l)/hmdn(l)
          ENDDO
         ELSE
          DO l = 1,nch
!            write(0,*) 'l,hm = ',l,hm(l)
          ENDDO
         ENDIF

         DO l = 1,nch
           rdamelt(l)  = (6.*hm(l)/(pi1*1000.))**(1./3.)
!          write(0,*) 'l,rdamelt = ',l,rdamelt(l)*1000.
         ENDDO

       ENDIF
    
#ifdef MPI
    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg

    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

    DO kz = kzb,kze
      db(kz) = 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sinit(kz,lt))
      pb(kz) = 1.0e5*piinit%flt1d(kz)**3.509
      DO jy = jyb,jye
        DO ix = ixb,ixe
#else
    ixb = 1
    ixe = nx-1
    jyb = 1
    jye = ny-1
    kzb = 1
    kze = nz-1
    DO kz = 1,nz-1
      db(kz) = 1.0e5*piinit%flt1d(kz)**2.509/(287.04*sinit(kz,lt))
      pb(kz) = 1.0e5*piinit%flt1d(kz)**3.509
!      DO jy=jyb_v5d,jye_v5d
!        DO ix=ixb_v5d,ixe_v5d
      DO jy = 1,ny-1
        DO ix = 1,nx-1
#endif
          pn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**3.509 - pb(kz)
!          dn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/(rd*sinit(kz,1)*(1.0+0.61*sinit(kz,2)))
          dn (ix,jy,kz) = 1.0e5*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))**2.509/(287.04*s(lt)%flt3d(ix,jy,kz) )
          tt0(ix,jy,kz) = s(lt)%flt3d(ix,jy,kz)*(piinit%flt1d(kz)+pi%flt3d(ix,jy,kz))
        ENDDO
      ENDDO
    ENDDO

        dxx(1:ixe) = 1.0/gx(3)%flt1d(1:ixe)
        dyy(1:jye) = 1.0/gy(3)%flt1d(1:jye)
        dzz(1:kze) = 1.0/gz(3)%flt1d(1:kze)


       IF ( present( recalcdbz ) ) THEN
        IF ( recalcdbz .ge. 1 ) THEN
         IF ( micro(1:3) .eq. 'TAK' ) THEN
          
    ixb = 1
    ixe = itile
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    jyb = 1
    jye = jtile
    if (jyend .eq. nyend) jye = jyend-jybeg

    kzb = 1
    kze = ktile
    if (kzend .eq. nzend) kze = kzend-kzbeg

    DO k = kzb,kze
      DO j = jyb,jye
        DO i = ixb,ixe

         ztmpi = 0.d0
         do ki=1,ntakit
           do ii=1,ntakid
              ia = lni + ii - 1 + (ki-1)*ntakid
             ztmpi = ztmpi+st(i,j,k,ia)*(((6.0/pii))*sxi(ii,ki))**2
            enddo
          enddo
         
         ztmpr = 0.d0
         DO il=1,ntakrd
            ztmpr = ztmpr + szr(il)*st(i,j,k,lnr+il-1)
         ENDDO

         ztmph = 0.d0
         DO il=1,ntakpd
            ztmph = ztmph + szf(il,1)*st(i,j,k,lnh+il-1) ! xf(il,1,i,j,k)
         ENDDO
!         ztmph = ztmph

         ztmphl = 0.d0
         DO il=1,ntakpd
            ztmphl = ztmphl + szf(il,2)*st(i,j,k,lnhl+il-1)  ! xf(il,2,i,j,k)
         ENDDO
         
         ztmp = 1.d6*(1.d6*0.224*ztmpi + ztmpr + 0.224*ztmph + 0.224*ztmphl)  ! convert per cm^3 to per m^3
         
         IF ( ztmp > 0.0d0 ) THEN
           dbz%flt3d(i,j,k) = Max(dbzmin, 10.0*Log10(ztmp) )
           dbzmax = Max( dbzmax, dbz%flt3d(i,j,k) )
         ELSE
           dbz%flt3d(i,j,k) = dbzmin
         ENDIF
         
         ENDDO
         ENDDO
         ENDDO
          
         ELSE
           IF ( associated( vzf%flt3d ) ) THEN
             CALL RADARDD02(nx,ny,nz,ng,na,st    &
             ,tt0,dbz%flt3d,dn, nz-1 ,cnoh,rho_qh, ipconc, 6, micro, 0, 0,vzf%flt3d)
           ELSE
             CALL RADARDD02(nx,ny,nz,ng,na,st    &
             ,tt0,dbz%flt3d,dn, nz-1 ,cnoh,rho_qh, ipconc, 6, micro, 0, 0,tt0)
           ENDIF

         ENDIF
        ENDIF
       ENDIF

       IF ( isfcphys > 0 ) THEN
        tcanp   = GET_VARIABLE_INDEX(gd, 'TCANP')
        wcanp   = GET_VARIABLE_INDEX(gd, 'WCANP')
        qav     = GET_VARIABLE_INDEX(gd, 'QAV')
        veg     = GET_VARIABLE_INDEX(gd, 'VEG')
        stype   = GET_VARIABLE_INDEX(gd, 'STYPE')
        tsrfc   = GET_VARIABLE_INDEX(gd, 'TSRFC')
        wsfc    = GET_VARIABLE_INDEX(gd, 'WSFC')
        tsoil   = GET_VARIABLE_INDEX(gd, 'TSOIL')
        wsoil   = GET_VARIABLE_INDEX(gd, 'WSOIL')
        vlai    = GET_VARIABLE_INDEX(gd, 'VLAI')
        albedo  = GET_VARIABLE_INDEX(gd, 'ALBEDO')
        rough   = GET_VARIABLE_INDEX(gd, 'ROUGH')

        eflx    = GET_VARIABLE_INDEX(gd, 'EFLX')
        fflx    = GET_VARIABLE_INDEX(gd, 'FFLX')
        uflx    = GET_VARIABLE_INDEX(gd, 'UFLX')
        vflx    = GET_VARIABLE_INDEX(gd, 'VFLX')
        tflx    = GET_VARIABLE_INDEX(gd, 'TFLX')
        qflx    = GET_VARIABLE_INDEX(gd, 'QFLX')  
        radsw   = GET_VARIABLE_INDEX(gd, 'RADSW')  
        radlw   = GET_VARIABLE_INDEX(gd, 'RADLW') 
       ENDIF
      
!
! Vis5D output
!
        
        IF ( onedoutput /= 0 ) THEN
           IF ( out1d_flag == 0 ) THEN
             out1d_flag = 1
           ! write titles
             write(33,11) (trim(v5dfields(i)), i=1,numvars)
           !  write(33,11) (v5dfields(i), i=1,numvars)
           ENDIF
           
           allocate( data_1d(MAXVARS,nz) )
           
        ENDIF
 11   format(1x,'  Time   Altitude',200('     ',a10))
 12   format(1x,i6,'  ',f9.2,200('  ',1x,1pe12.5))
        
        nr=ny_v5d
        nc=nx_v5d
        nl(1)=nz_v5d-1

        ibc = 0
        jbc = 0
        IF ( bcx .eq. 2 .or. ny .le. 2 ) ibc = 1
        IF ( bcy .eq. 2 .or. nx .le. 2 ) jbc = 1

        nr0 = Max(2,nr-1)
        nr1 = Min(2-jbc,ny_v5d-1)
        nr2 = Max(1    ,ny_v5d-2+jbc)
        
      
        CALL GET_VARIABLE(gd, 'TIME',       time)
        write(luno,*) 'V5DOUT: time = ', time

        idatime = time
        hh      = idatime/3600
        idatime = idatime - 3600*hh
        mm      = idatime/60
        ss      = idatime - 60*mm
     
        it = itv5d
        
        itimes(it) = 10000*hh+100*mm+ss
        write(luno,*) 'V5DOUT itimes(it) = ',it, itimes(it)
        idates(it) = 2006001
        
          write(luno,'(a,a)') 'V5DOUT: Writing Vis5D file: ',v5dfilename


!        it = it+1
        
        DO iv = 1, numvars
        
         fac = 1./Float(vis5dstridex*vis5dstridex*vis5dstridez)
         
         IF ( ny .le. 2 .or. nx .le. 2 ) fac = 1./Float(vis5dstridex*vis5dstridez)

        IF ( iv5dwritten(iv) .eq. 0 ) THEN
          iv5dwritten(iv) = 1
        ELSE
          CYCLE
        ENDIF
         tem1(:,:,:) = 0.0
         
        dunit = '              '
        write(dunit,'(a)') units(idx(iv))
        call v5dsetunits(iv, dunit )

        IF ( idx(iv) .eq. 1 ) THEN   ! THETAP
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lt) - sinit(k+kzb_v5d-1,lt)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 3 ) THEN

          is = u%istag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.5*(u%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
     &                             u%flt3d(i+is+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 4 ) THEN

          js = v%jstag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.5*(v%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
     &                             v%flt3d(i+ixb_v5d-1,nr-j+js+jyb_v5d-1,k+kzb_v5d-1))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 5 ) THEN

          ks = w%kstag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.5*(w%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) + &
     &                             w%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1+ks))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 91 ) THEN ! theta-e

          is = u%istag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                pres = pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) + pb(k)
                qwv = max(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv) , 0.0 )
                theta1 = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lt)
                tem1(j,i,k) = thetae(qwv,theta1,pres)
              ENDDO
            ENDDO
          ENDDO

      ELSEIF ( idx(iv) .eq. 92 ) THEN  ! divergence/convergence

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=2,nr2
            tem1(j,i,k)  =  1000.0*( (u%flt3d(i+ixb_v5d-1+1,nr-j+jyb_v5d-1,k+kzb_v5d-1) -                            &
                                         u%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))*gx(4)%flt1d(i+ixb_v5d-1)   &
                                   + (v%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1+1,k+kzb_v5d-1) -                            &
                                          v%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))*gy(4)%flt1d(nr-j+jyb_v5d-1) )
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 103 ) THEN

          is = u%istag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = u%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 104 ) THEN

          js = v%jstag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = v%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 105 ) THEN

          ks = w%kstag
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = w%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

      ELSEIF ( idx(iv) .eq. 11 ) THEN  ! vorticity

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=2,nr2
            tem1(j,i,k)  =  1000.0*wz%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 12 ) THEN   ! THETA

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lt)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 13 ) THEN   ! THETAV

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                IF ( lthetav > 0 ) THEN
                  tem1(j,i,k) = xtra(lthetav)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                ELSE
                  tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lt)* &
                       (1. + 0.61*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv))
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 16 ) THEN   !  Temperture (K)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = tt0(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 17 ) THEN   ! Dewpoint Temperture (K)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
            x = Log(Max(1.0e-6,st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv))*  &
             (pb(k)+pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))/380.0)
            tx = (35.86*x - 17.27*273.15)/(x - 17.27)
                tem1(j,i,k) = tx
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 23 ) THEN   !  TKE

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = km%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 28 ) THEN ! X-Position

          DO k=1,nl(iv)
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = x_sw_loc + gx(1)%flt1d(i)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 29 ) THEN ! Y-Position

          DO k=1,nl(iv)
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = y_sw_loc + gy(1)%flt1d(nr-j)
              ENDDO
            ENDDO
          ENDDO
      
      
         ELSEIF ( idx(iv) .eq. 45 .or. idx(iv) == 46  ) THEN ! SSI
      
          IF ( tabqvs(1) == 0.0 .and. itabflag == 0 ) THEN
           itabflag = 1

           do l = 1,nqsat
             temq = 163.15 + (l-1)*fqsat
             tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
             tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
           end do
      
          ENDIF
          IF ( idx(iv) .eq. 45 ) THEN ! SSI
           
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

            pqs = 380.0/(pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)+pb(k))
            qwv = max(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv),0.0)
            temg = tt0(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)

            ltemq = (temg-163.15)/fqsat+1.5
            ltemq = Min( nqsat, Max(1,ltemq) )
            qis = pqs*tabqis(ltemq)

            tem1(j,i,k)  = 100.*qwv/qis - 100.

            ENDDO
           ENDDO
          ENDDO
          
          ELSE  ! SSW
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
            pqs = 380.0/(pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)+pb(k))
            qwv = max(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv),0.0)
            temg = tt0(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)

            ltemq = (temg-163.15)/fqsat+1.5
            ltemq = Min( nqsat, Max(1,ltemq) )
            qvs = pqs*tabqvs(ltemq)

            tem1(j,i,k)  = 100.*qwv/qvs - 100.
            ENDDO
           ENDDO
          ENDDO
          
          ENDIF
      

         ELSEIF ( idx(iv) .eq. 55 ) THEN ! pressure perturbation PP (Pa)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 56 ) THEN ! pressure P (mb)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e-2*( pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) + pb(k) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 57 ) THEN ! perturbation exner PI (mb/mb)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = pi%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 420 .and. micro(1:3) == 'TAK' ) THEN ! ccw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

              tmp = 0.0
               DO l = 1,lr100-1
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
               ENDDO
              
              tem1(j,i,k) =  tmp ! leave in #/cm**3

              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 420 .and. lnc .gt. 1 .and. ipconc .ge. 1 ) THEN ! ccw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e-6*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnc) ! per cm**3
              IF ( tem1(j,i,k) .lt. 1.0e-3 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 424 .and. micro(1:3) == 'TAK' ) THEN ! crw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

              tmp = 0.0
               DO l = lr100,ntakrd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
               ENDDO
              
              tem1(j,i,k) =  1.e3*tmp ! convert to #/L

              IF ( tem1(j,i,k) .lt. 1.0e-6 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( idx(iv) .eq. 424 .and. lnr .gt. 1 ) THEN ! crw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) ! per L**3
              ! IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 425 .and. lns .gt. 1 .and. ipconc .ge. 4 ) THEN ! csw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
!              tem1(j,i,k) = Min( ccimx, 1.e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns) ) ! per L**3
              tem1(j,i,k) = 1.e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)  ! per L**3
              IF ( tem1(j,i,k) .lt. 1.0e-3 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 426 .and. lngl .gt. 1 .and. ipconc .ge. 5) THEN ! cgl

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngl) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 427 .and. lngm .gt. 1 .and. ipconc .ge. 5) THEN ! cgm

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngm) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 428 .and. lngh .gt. 1 .and. ipconc .ge. 5) THEN ! cgh

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngh) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 429 .and. micro(1:3) == 'TAK' ) THEN ! cfw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                 tmp = 0.0
                 IF ( lf75 > 1 ) THEN
                  DO l = 1,lf75
                    ia = lnh + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                    ia = lnhl + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDIF

              tem1(j,i,k) = tmp*1000. ! per liter

            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 429 .and. lnf .gt. 1 .and. ipconc .ge. 5 ) THEN ! cfw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnf) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 430 .and. lnhl .gt. 1 .and. ipconc .ge. 5 ) THEN ! chw+chl

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) + &
     &                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 431 .and. ipconc .ge. 5 ) THEN ! cgtot

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngl) +   &
                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngm) +   &
                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngh) +   &
                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnf)    ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
              ELSEIF ( micro(1:1) .eq. 'Z' ) THEN
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh)
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 433 .and. (ipconc .ge. 1  .or.  micro(1:3) == 'TAK' )) THEN ! cwrad

          IF (micro(1:3) == 'TAK' ) THEN
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
              tmpn = 0.0
!               DO n = ls250+1,lfmax
               DO n = 1,lr100-1
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n)
                tmpn = tmpn + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)
               ENDDO
               qr = 1.e6*1.e-3*tmp ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
              
              chw = 1.e6*tmpn
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 1.e6*0.5*(6.*qr/(1000.*3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
          ELSE
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) .gt. qxmin(lc) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnc) .gt. cxmin ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) * &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc)/ &
     &                        (1000.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnc))
                tem1(j,i,k) = 1.e6*(3*tem1(j,i,k)/(4.0*3.14159))**(1./3.)
                tem1(j,i,k) = Min( tem1(j,i,k) , 1.e6*cwradx )
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 436 .and. (ipconc .ge. 3 .or. lnr > 1) .and.  micro(1:3) .ne. 'TAK'  ) THEN ! rwrad

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. cxmin ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)/         &
     &                        (1000.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr))
                tem1(j,i,k) = 1.e6*(3.0*tem1(j,i,k)/(4.0*3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 436 .and.  micro(1:3) .eq. 'TAK'  ) THEN ! rwrad
         
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
              tmpn = 0.0
               DO n = lr100,ntakrd
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n)
                tmpn = tmpn + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)
               ENDDO
               qr = 1.e6*1.e-3*tmp ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.
              
              chw = 1.e6*tmpn
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 0.5*1.e6*(6.*qr/(1000.*3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 464 .and. microp(1:3) == 'TAK' ) THEN ! cidia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
              chw = 0.0
                 DO ki=1,ntakit
                  DO ii=1,ntakid
                    ia = lni + ii - 1 + (ki-1)*ntakid
                    chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)*sxi(ii,ki)
                  ENDDO
                 ENDDO
                 
                 IF ( tmp > 1.e-15 .and. chw > 1.e-20 ) THEN
                  tem1(j,i,k) =  1.e6*0.277823*((1.e-3*tmp/chw)**(0.359971))  ! factor of 1.e6 to convert to microns, factor of 1.e-3 to convert to kg
                 ELSE
                  tem1(j,i,k) = 0.0
                 ENDIF
              ENDDO
            ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 464 .and. ipconc .ge. 1 ) THEN ! cidia (columns)

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) .gt. qxmin(li) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni) .gt. cxmin ) THEN
                tmp = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) *    &
     &                dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)/        &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni)  ! mass per crystal
              IF ( ixtaltype == 1 ) THEN ! column
                 tem1(j,i,k) =  1.e6*0.1871*(tmp**(0.3429))  ! factor of 1.e6 to convert to microns
               ELSEIF  ( ixtaltype == 2 ) THEN ! disk
                 tem1(j,i,k) =  1.e6*0.277823*(tmp**(0.359971))  ! factor of 1.e6 to convert to microns
               ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 465 .and. ipconc .ge. 4 ) THEN ! swdia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              
              hwdn = 100.
              IF ( lsw > 1 ) THEN
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsw) > qxmin(ls) .and.  &
                     st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. qxmin(ls) ) THEN
                   
                   fsw = Min(1.0, st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsw)/st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls))
                   hwdn = fsw*1000. + (1. - fsw)*100.
                   
                ENDIF
              ENDIF
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. qxmin(ls) .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns) .gt. cxmin ) THEN

                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/         &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)

              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 466 .and. lngl .gt. 1 .and. ipconc .ge. 5 ) THEN ! gldia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl) .gt. 1.e-6 .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngl) .gt. 1.e-3 ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl)/        &
     &                        (300.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngl))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 467 .and. lngm .gt. 1 .and. ipconc .ge. 5 ) THEN ! gmdia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm) .gt. 1.e-6 .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngm) .gt. 1.e-3 ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
                              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm)/        &
                              (500.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngm))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 468 .and. lngh .gt. 1 .and. ipconc .ge. 5 ) THEN ! ghdia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh) .gt. 1.e-6 .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngh) .gt. 1.e-3 ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh)/        &
     &                        (700.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lngh))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 469 .and. lnf .gt. 1 .and. ipconc .ge. 5 ) THEN ! fwdia

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) .gt. 1.e-6 .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnf) .gt. 1.e-3 ) THEN
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf)/         &
     &                        (800.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnf))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
         ELSEIF ( idx(iv) .eq. 470 .and.  micro(1:3) .eq. 'HCM'  ) THEN ! hwdia
          
!          write(0,*) 'V5DOUT: HCM hail diam'
          
          hwdn = 900.
          hwdnl(:) = hwdn
          hwvol(1:nch) = hm(1:nch)/hwdn

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              qr = 0.0
               DO n = 1,nch
                qr = qr + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh+n-1)
               ENDDO

             IF ( micro(1:4) .eq. 'HCMD' ) THEN
               DO n = 1,nch
                 hwdnl(n) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh+n-1)
                 hwvol(n) = hm(n)/hwdnl(n)
               ENDDO
             ELSEIF ( micro(1:4) .eq. 'HCMV' ) THEN
               DO n = 1,nch
                 hwvol(n) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh+n-1)
                 hwdnl(n) = hm(n)/hwvol(n)
               ENDDO
             ENDIF
             
              
             IF ( micro(1:4) .eq. 'HCMD' .or. micro(1:4) .eq. 'HCMV'  ) THEN

              voltot = 0.0
              
               DO n = 1,nch
                 voltot = voltot + hwvol(n)
               ENDDO
              
             ELSE
               voltot = qr/hwdn
             ENDIF

              chw = 0.0
               DO n = 1,nch
               ! count number of particles: rho_air*q/mass
               
                chw = chw + dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*  &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh+n-1)/hm(n)
               ENDDO
              
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 1.e6*(6.*voltot/(3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
          
!          ENDIF

         ELSEIF ( idx(iv) .eq. 470 .and.  micro(1:3) .eq. 'TAK'  ) THEN ! hwdia
          
!          write(0,*) 'V5DOUT: HCM hail diam'
          
          hwdn = 1000.*rhof(1)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
!               DO n = ls250+1,lfmax
               DO n = 1,lfmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*sxf(n,1)
               ENDDO
               qr = 1.e6*1.e-3*tmp ! 1.e6 converts cm^-3 to m^-3; 1e-3 converts gm to kg.

              tmp = 0.0
!               DO n = ls250+1,lfmax
               DO n = 1,lfmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)
               ENDDO
              
              chw = 1.e6*tmp
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 470 .and. lnh .gt. 1 .and. ipconc .ge. 5 ) THEN ! hwdia

          hwdn = 900.
          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
!          cwch =  6.0*pii*gamma( (xnu(lh) + 1.)/xmu(lh) )/gamma( (xnu(lh) + 2.)/xmu(lh) )
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and.   &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. cxmin ) THEN
              IF ( lvh .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*      &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ENDIF
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*           &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh))
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

          
         ELSEIF ( idx(iv) .eq. 472 .and.  micro(1:3) .eq. 'TAK'  ) THEN ! hldia
          
!          write(0,*) 'V5DOUT: TAK HLDIA, ls250 = ',ls250
          
          hwdn = 1000.*rhof(2)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*sxf(n,2)
               ENDDO
               qr = 1.e6*1.e-3*tmp! /dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)

              tmp = 0.0
               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)
               ENDDO
              
              chw = 1.e6*tmp
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
         ELSEIF ( idx(iv) .eq. 436 .and.  micro(1:3) .eq. 'TAK'  ) THEN ! rwrad
          
!          write(0,*) 'V5DOUT: TAK RWRAD, lr100 = ',lr100
          
          hwdn = 1000.
          tmpmx = 0.0
          tmpn = 0.0

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO n = lr100,lmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n)
               ENDDO
               qr = 1.e6*1.e-3*tmp ! /dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)

              tmpmx = Max(tmpmx, qr)
              
              tmp = 0.0
               DO n = lr100,lmax
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)
               ENDDO
              
              chw = 1.e6*tmp ! convert cm^-3 to m^-3

              tmpn = Max(tmpn, chw)
              
              IF ( qr > 1.e-8 .and. chw > 1.e-3 ) THEN
              
                tem1(j,i,k) = 0.5e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
!          write(0,*) 'max qr,chw = ',tmpmx,tmpn

         ELSEIF ( idx(iv) .eq. 470 .and. ipconc .lt. 5 ) THEN ! hwdia

          hwdn = 900.
          cno = cnoh
          IF ( micro(1:1) .eq. 'Z' ) THEN
            hwdn = rho_qh
          ELSEIF ( micro(1:5) == 'ICE10' ) THEN
            hwdn = rho_qh_i10
            cno = cnoh_i10
          ENDIF
!          cwch =  6.0*pii*gamma( (xnu(lh) + 1.)/xmu(lh) )/gamma( (xnu(lh) + 2.)/xmu(lh) )
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) ) THEN
              IF ( lvh .gt. 1 ) THEN
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) .gt. 1.e-9 ) THEN
                  hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*     &
     &                   st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/  &
     &                   st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
                ELSE
                  hwdn = rho_qh
                ENDIF
              ENDIF

                xdia = (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)*  &
     &                  dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)/((3.14159)*hwdn*cno))**(0.25) 
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*   &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/(hwdn*xdia*cno)
                tem1(j,i,k) = 1.e6*(6.*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 471 .and. lnh .gt. 1 .and. ipconc .ge. 5 ) THEN ! hwvdia

          hwdn = 900.
          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. cxmin ) THEN
              IF ( lvh .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*     &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/  &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ENDIF
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*           &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh))
                tem1(j,i,k) = 1.e6*(6.0*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 472 .and. lnhl .gt. 1 .and. ipconc .ge. 5 ) THEN ! hldia

          hwdn = 900.
!          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.   &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. cxmin ) THEN
              IF ( lvhl .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
              ENDIF
                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))
                tem1(j,i,k) = 1.e6*(6.0*tem1(j,i,k)/(3.14159))**(1./3.)
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( ( idx(iv) .eq. 490 .or. idx(iv) .eq. 492 .or. idx(iv) .eq. 498 )  .and. lnh .gt. 1  ) THEN ! maximum mass diameter hmdia or mass-weighted diam HMWD

          IF ( micro(1:1) == 'Z' ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               alp = alphah

              tem1(j,i,k) = 0.0

               hwdn = rho_qh
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) < qxmin(lh) .or. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) < cxmin )  CYCLE
              
              IF ( lvh .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ENDIF
               
              IF ( lzh > 1 ) THEN ! find alpha
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzh) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzh)
                alpha1 = 0

                
                rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
                
               ENDIF
              ENDIF ! lzh
              
              ! have alp, now need mean volume diameter
              
                tmp = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh))
                diam = (6.0*tmp/(3.14159))**(1./3.)
                
                IF ( idx(iv) .eq. 490 ) THEN ! max mass diam
                  tem1(j,i,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ELSEIF ( idx(iv) .eq. 492 ) THEN ! mass-weighted
                  tem1(j,i,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ELSEIF ( idx(iv) .eq. 498 ) THEN ! characteristic
                  tem1(j,i,k) = 1.e6*diam*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ENDIF
              
              ENDDO
            ENDDO
          ENDDO

          ELSEIF ( micro(1:3) == 'TAK' ) THEN

              
          hwdn = 1000.*rhof(1)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

                tmp = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0
               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*sxf(n,1) ! total mass
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)  ! total number 
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*sxf(n,1)*srf(n,1)*2. ! D*N(D)*m(D)
                area = area + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*srf(n,1)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*2.*srf(n,1)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO

              
          !    chw = 1.e6*tmp
              
              IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                IF ( idx(iv) .eq. 492 ) THEN
                  tem1(j,i,k) =  0.01*1.e6*tmp/qr  ! 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

!                   qr = 1.e6*1.e-3*qr ! /dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
!                   chw = 1.e6*chw
!                   tem1(j,i,k) = 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

                ELSEIF ( idx(iv) .eq. 490 ) THEN
                  tem1(j,i,k) =  1.e6*0.01*areaw/area  ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
              
            ENDDO
           ENDDO
          ENDDO


          ENDIF
          
!          hwdn = 900.
!!          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
!          DO k=1,nz_v5d-1
!            DO i=1,nc-1
!              DO j=1,nr-1
!              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.   &
!     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. cxmin ) THEN
!              IF ( lvhl .gt. 1 ) THEN
!               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
!     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
!     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
!              ENDIF
!                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
!     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/        &
!     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))
!                tem1(j,i,k) = 1.e6*(6.0*tem1(j,i,k)/(3.14159))**(1./3.)
!              ELSE
!                tem1(j,i,k) = 0.0
!              ENDIF
!            ENDDO
!           ENDDO
!          ENDDO


         ELSEIF ( ( idx(iv) .eq. 491 .or. idx(iv) .eq. 493 .or. idx(iv) .eq. 499 )  .and. lnhl .gt. 1  ) THEN ! maximum mass diameter hlmdia or mass-weighted diam HLMWD

          IF ( micro(1:1) == 'Z' ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               alp = alphahl

              tem1(j,i,k) = 0.0

               hwdn = rho_qhl
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) < qxmin(lhl) .or. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) <= cxmin )  CYCLE
              
              IF ( lvhl .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
              ENDIF
               
              IF ( lzhl > 1 ) THEN ! find alpha
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl)
                alpha1 = 0

                
                rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
                
               ENDIF
              ENDIF ! lzhl
              
              ! have alp, now need mean volume diameter
              
                tmp = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))
                diam = (6.0*tmp/(3.14159))**(1./3.)
                
                IF ( idx(iv) .eq. 491 ) THEN ! max mass diam
                  tem1(j,i,k) = 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ELSEIF ( idx(iv) .eq. 493 ) THEN ! mass-weighted
                  tem1(j,i,k) = 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ELSEIF ( idx(iv) .eq. 499 ) THEN ! characteristic
                  tem1(j,i,k) = 1.e6*diam*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)
                ENDIF
              
              ENDDO
            ENDDO
          ENDDO

          ELSEIF ( micro(1:3) == 'TAK' ) THEN

              
          hwdn = 1000.*rhof(2)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

                tmp = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0
               DO n = ls250+1,lfmax
!               DO n = 1,lfmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*sxf(n,2) ! total mass
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)  ! total number 
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*sxf(n,2)*srf(n,2)*2. ! D*N(D)*m(D)
                area = area + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*srf(n,2)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*2.*srf(n,2)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO

              
          !    chw = 1.e6*tmp
              
              IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                IF ( idx(iv) .eq. 493 ) THEN
                  tem1(j,i,k) =  0.01*1.e6*tmp/qr  ! 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

!                   qr = 1.e6*1.e-3*qr ! /dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
!                   chw = 1.e6*chw
!                   tem1(j,i,k) = 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

                ELSEIF ( idx(iv) .eq. 491 ) THEN
                  tem1(j,i,k) =  1.e6*0.01*areaw/area  ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
              
            ENDDO
           ENDDO
          ENDDO


          ENDIF
          
!          hwdn = 900.
!!          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
!          DO k=1,nz_v5d-1
!            DO i=1,nc-1
!              DO j=1,nr-1
!              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.   &
!     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. cxmin ) THEN
!              IF ( lvhl .gt. 1 ) THEN
!               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
!     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
!     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
!              ENDIF
!                tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
!     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/        &
!     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))
!                tem1(j,i,k) = 1.e6*(6.0*tem1(j,i,k)/(3.14159))**(1./3.)
!              ELSE
!                tem1(j,i,k) = 0.0
!              ENDIF
!            ENDDO
!           ENDDO
!          ENDDO


         ELSEIF ( ( idx(iv) .eq. 494 .or. idx(iv) .eq. 495 .or. idx(iv) == 496 .or. idx(iv) == 497 )  .and. lnr .gt. 1  ) THEN ! mass-weighted diam RWMWD

          IF ( micro(1:1) == 'Z' ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               alp = alphar

              tem1(j,i,k) = 0.0

               hwdn = rho_qr
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) < qxmin(lr) .or. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) <= 0.0 )  CYCLE
              
               
              IF ( lzr > 1 ) THEN ! find alpha
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr)
                alpha1 = 0

                
                IF ( imurain == 1 ) THEN
                
                rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

                alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!               print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
                DO ii = 1,10
                 IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                  alpha1 = Max( alphamin, Min( alphamax, alp ) )
                  alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!                print*,'i,alp = ',i,alp
                  alp = Max( alphamin, Min( alphamax, alp ) )
                ENDDO
                
                ELSE ! imurain == 3
                  nrx = chw
                  vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr/(1000.*chw)
                  alp = 36.*(alpha1+2.0)*nrx*vr**2/(z*pii**2) - 1.
                  DO ii = 1,20
                   IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
                   alpha1 = Max( rnumin, Min( rnumax, alp ) )
                   alp = 36.*(alpha1+2.0)*nrx*vr**2/(z*pii**2) - 1.
                   alp = Max( rnumin, Min( rnumax, alp ) )
                  ENDDO
                
                ENDIF
                
               ENDIF
              
              ELSE
                alp = alphar
              ENDIF ! lzr
              
              ! have alp, now need mean volume diameter
              
              ! tmp is the mean volume (xv)
                tmp = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)/        &
     &                        (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr))
                diam = (6.0*tmp/(3.14159))**(1./3.) ! mean volume diameter
                
                IF ( idx(iv) .eq. 495 ) THEN ! max mass diam
                  tem1(j,i,k) = Min(9000., 1.e6*diam*(3.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.) )
                ELSEIF ( idx(iv) .eq. 494 ) THEN ! mass-weighted
                  IF ( imurain == 1 ) THEN
                    tem1(j,i,k) = Min(9000., 1.e6*diam*(4.0 + alp)*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.) )
                  ELSE
                    tem1(j,i,k) = Min(9000., 1.e6*(4./3.+alp)**(1./3.)*diam )
                  ENDIF
                 ELSEIF ( idx(iv) .eq. 496 ) THEN ! median volume diameter
                  IF ( imurain == 1 ) THEN
                    tem1(j,i,k) = Min(9000., 1.e6*(3.67+alp)*diam*((3.+alp)*(2.+alp)*(1. + alp) )**(-1./3.)  )
                  ELSE
                    tem1(j,i,k) = Min(9000., 1.e6*(1.678+alp)**(1./3.)*diam )
                  ENDIF
                ELSEIF ( idx(iv) .eq. 497 ) THEN ! mean volume (mass) diameter
                    tem1(j,i,k) = 1.e6*diam
                ENDIF
              
              ENDDO
            ENDDO
          ENDDO

          ELSEIF ( micro(1:3) == 'TAK' ) THEN

              
          hwdn = 1000.

!          IF ( idx(iv) .eq. 495 ) THEN
!           write(0,*) 'Calc. ',idx(iv),' for TAK' ! ' RWMDIA'
!          ENDIF
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

                tmp = 0.0
                qr = 0.0
                chw = 0.0
                area = 0.0
                areaw = 0.0
               DO n = lr100,lmax
!               DO n = 1,lfmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n) ! total mass
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)  ! total number 
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n)*sr(n)*2. ! D*N(D)*m(D) for mass-weighted diameter
                area = area + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sr(n)**2. !  total area for area-weighted diameter (pi cancels out)
                areaw = areaw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*2.*sr(n)**3. ! for area-weighted diameter (D=2*r, A=r**2, Pi cancels out)
               ENDDO

              
          !    chw = 1.e6*tmp
              
              IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                IF ( idx(iv) .eq. 494 ) THEN
                  tem1(j,i,k) =  Min(9000., 0.01*1.e6*tmp/qr ) ! 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

!                   qr = 1.e6*1.e-3*qr ! /dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
!                   chw = 1.e6*chw
!                   tem1(j,i,k) = 1.e6*(6.*qr/(hwdn*3.14159*chw))**(1./3.)

                ELSEIF ( idx(iv) .eq. 495 ) THEN
                  tem1(j,i,k) =  Min(9000., 1.e6*0.01*areaw/area ) ! 0.01 converts cm to m, then 1.e6 converts m to microns
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
              
            ENDDO
           ENDDO
          ENDDO


          ENDIF

         ELSEIF ( idx(iv) .eq. 474 .and. (( lnh .gt. 1 .and. ipconc .ge. 5 ) .or. micro(1:3) == 'TAK'  ) ) THEN ! hwmas

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              
              tem1(j,i,k) = 0.0

              IF ( micro(1:3) == 'TAK' ) THEN

                qr = 0.0
                chw = 0.0
               DO n = ls250+1,lfmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)*sxf(n,1) ! total mass (grams per cm^3)
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh+n-1)  ! total number (per cm^3)
               ENDDO
               IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                   tem1(j,i,k) =  1.e3*qr/chw  ! milligrams
               ENDIF

              
              ELSE
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. cxmin ) THEN
                tem1(j,i,k) = 1.e6*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*           &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/        &
     &                        (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh))
              ENDIF
              
              ENDIF
            ENDDO
           ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 475 .and. lnhl .gt. 1 .and. (ipconc .ge. 5 .or. micro(1:3) == 'TAK') ) THEN ! hlmas

!          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              
              tem1(j,i,k) = 0.0

              IF ( micro(1:3) == 'TAK' ) THEN

                qr = 0.0
                chw = 0.0
               DO n = ls250+1,lfmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)*sxf(n,2) ! total mass (grams per cm^3)
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl+n-1)  ! total number (per cm^3)
               ENDDO
               IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                   tem1(j,i,k) =  1.e3*qr/chw  ! milligrams
               ENDIF

              
              ELSE
              
              
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.   &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. cxmin ) THEN
                tem1(j,i,k) = 1.e6*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/        &
     &                        (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))
              ENDIF
              
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 476 .and. lns .gt. 1 .and. ipconc .ge. 5 ) THEN ! swmas

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. qxmin(ls) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns) .gt. cxmin ) THEN
                tem1(j,i,k) = 1.e6*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*           &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/        &
     &                        (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns))
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 477 .and. lnr .gt. 1 .and. (ipconc .ge. 5 .or. micro(1:3) == 'TAK')) THEN ! rwmas

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

              tem1(j,i,k) = 0.0

              IF ( micro(1:3) == 'TAK' ) THEN

                qr = 0.0
                chw = 0.0
               DO n = lr100,lmax
                qr  = qr  + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)*sx(n) ! total mass
                chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr+n-1)  ! total number 
               ENDDO
               IF ( 1.e6*1.e-3*qr > 1.e-8 .and. 1.e6*chw > 1.e-3 ) THEN
                   tem1(j,i,k) =  1.e3*qr/chw  ! milligrams
               ENDIF

              
              ELSE
              

              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) .and.  &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. cxmin ) THEN
                tem1(j,i,k) = 1.e6*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*           &
     &                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)/        &
     &                        (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr))
              ENDIF
              
              ENDIF
              
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 485 .and. lnh .gt. 1 .and. ipconc .ge. 5 ) THEN ! N0HW

          hwdn = 900.
          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. cxmin ) THEN
              IF ( lvh .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*      &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ENDIF
                tmp = (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh))/     &
     &                (dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*               &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh))
                tmpg = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh)*(tmp*(3.14159))**(1./3.)
                IF ( tmpg .gt. 1. ) THEN
                  tem1(j,i,k) = Log10( tmpg )
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 486 .and. lnhl .gt. 1 .and. ipconc .ge. 5 ) THEN ! N0HL

          hwdn = 800.
          IF ( micro(1:1) .eq. 'Z' ) hwdn = rho_qh
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. cxmin ) THEN
              IF ( lvhl .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
              ENDIF
                tmp = (hwdn*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl))/  &
     &                (dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*             &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl))
                tmpg = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl)*(tmp*(3.14159))**(1./3.)
                IF ( tmpg .gt. 0.0 ) THEN
                  tem1(j,i,k) = Log10( tmpg )
                ELSE
                  tem1(j,i,k) = 0.0
                ENDIF
              ELSE
                tem1(j,i,k) = 0.0
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 242 .and. micro(1:3) == 'TAK' ) THEN ! chw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

              tmp = 0.0
               DO l = 1,ntakpd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh-1+l)
               ENDDO
              
              tem1(j,i,k) =  1.e6*tmp

              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 243 .and. micro(1:3) == 'TAK' ) THEN ! chl

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1

              tmp = 0.0
               DO l = ls250,ntakpd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl-1+l)
               ENDDO
              
              tem1(j,i,k) =  1.e6*tmp

              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 243 .and. lnhl .gt. 1 .and. ipconc .ge. 5 ) THEN ! chl

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 244 .and. lnhlf .gt. 1 .and. ipconc .ge. 5 ) THEN ! chlf

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhlf) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 245 .and. lnhf .gt. 1 .and. ipconc .ge. 5 ) THEN ! chf

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhf) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 242 .and. micro(1:3) == 'HCM' ) THEN ! chw

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO n = 1,nch
                tmp = tmp + dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*  &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh+n-1)/hm(n)
!                tmp = tmp + dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) &
!                *st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh+n)
               ENDDO
              
              tem1(j,i,k) =  (  tmp )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 242 .and. lnh .gt. 1 .and. ipconc .ge. 5 ) THEN ! chw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
          
         ELSEIF ( idx(iv) .eq. 241 .and. lnf .gt. 1 .and. ipconc .ge. 5 ) THEN ! cfw

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnf) ! per m**3
              IF ( tem1(j,i,k) .lt. 1.0e-2 ) tem1(j,i,k) = 0.0
            ENDDO
           ENDDO
          ENDDO
          

         ELSEIF ( idx(iv) .eq. 421 .and. microp(1:3) == 'TAK' ) THEN ! cci

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
                 DO ki=1,ntakit
                  DO ii=1,ntakid
                    ia = lni + ii - 1 + (ki-1)*ntakid
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDDO
                 IF ( lf75 > 1 ) THEN
                  DO l = 1,lf75
                    ia = lnh + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                    ia = lnhl + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDIF
                 
                tem1(j,i,k) = 1.e3*tmp
              ENDDO
            ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 421 .and. lni .gt. 1 .and. ipconc .ge. 1 ) THEN ! cci

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
!              tem1(j,i,k) =  Min( ccimx, 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni) )  ! per liter
              tem1(j,i,k) =   1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni)   ! per liter
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 422 .and. microp(1:3) == 'TAK' ) THEN ! cip
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
                 DO ki=1,ntakit
                  DO ii=1,ntakid
                    ia = lni + ii - 1 + (ki-1)*ntakid
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDDO
                 
!                tem1(j,i,k) = Min( ccimx, 1.e3*tmp )
                tem1(j,i,k) = 1.e3*tmp
              ENDDO
            ENDDO
          ENDDO
         

         ELSEIF ( idx(iv) .eq. 422 .and. lnip .gt. 1 .and. ipconc .ge. 1 ) THEN ! cip

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
!              tem1(j,i,k) =  Min( ccimx, 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnip))  ! per liter
              tem1(j,i,k) = 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnip) ! per liter
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 423 .and. lnir .gt. 1 .and. ipconc .ge. 1 ) THEN ! cir

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
!              tem1(j,i,k) =  Min( ccimx, 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnir))  ! per liter
              tem1(j,i,k) =  1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnir)  ! per liter
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 432 .and. microp(1:3) == 'TAK' ) THEN ! citot

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
                 DO ki=1,ntakit
                  DO ii=1,ntakid
                    ia = lni + ii - 1 + (ki-1)*ntakid
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDDO
                 IF ( lf75 > 1 ) THEN
                  DO l = 1,lf75
                    ia = lnh + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                    ia = lnhl + l - 1
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDIF
                 
                tem1(j,i,k) = Min( ccimx, 1.e3*tmp )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 432 .and. ( ipconc .ge. 1 .or. lni > 1 ) ) THEN ! citot

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:7) .eq. 'ICE10DM' ) THEN    !  all ice crystals
              tem1(j,i,k) =  Min( ccimx,          &
                     1.0e-3*(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni)  +  &
                             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnip) +  &
                             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnir) ) )
              ELSEIF ( micro(1:1) .eq. 'Z' .or. micro(1:4) == 'THOM' ) THEN
!              tem1(j,i,k) =  Min( 5.0e4, 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni))  ! per liter
              tem1(j,i,k) =  Min( ccimx, 1.0e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni))  ! per liter
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          
          

         ELSEIF ( idx(iv) .eq. 434 .and. ipconc .ge. 2 .and. lsat .ge. 1 ) THEN ! SS

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = xtra(lsat)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 435 .and. ipconc .ge. 2 .and. lsati .ge. 1 ) THEN ! SSI

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = xtra(lsati)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 438 .and. lccn .gt. 1 ) THEN ! unactivated ccn

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lccn) .lt. 1.e20 .and. &
     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lccn) .gt. -1.e10 ) THEN
              tem1(j,i,k) = 1.0e-6*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lccn) ! per cm**3
              IF ( micro(1:3) == 'TAK' ) tem1(j,i,k) = 1.e6*tem1(j,i,k)
              
              IF ( tem1(j,i,k) .lt. 1.0e-3 ) tem1(j,i,k) = 0.0
              ELSE
              tem1(j,i,k) = MISSING
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 440 .and. lnchaff .gt. 1 ) THEN ! chaff fibers

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnchaff)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 441 .and. lnox .gt. 1 ) THEN ! CLNOX in ppbv

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnox)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 442 .and. lccna .gt. 1 ) THEN ! activated CCN
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.e-6*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lccna)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 443 .and. lcina .gt. 1 ) THEN ! activated ice nuclei

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 0.001*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lcina)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 502 ) THEN ! qvp

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j = 1,nr-1 ! nr1,nr2
              tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv) - sinit(k+kzb_v5d-1,lv) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 503 ) THEN ! qv

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j = 1,nr-1 ! nr1,nr2
              tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lv) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 504 ) THEN ! qc

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j = 1,nr-1 ! nr1,nr2
              tem1(j,i,k) =   1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 604 ) THEN ! cwc

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) )*  &
     &                        dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 605 ) THEN ! rc rain content

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) )*  &
     &                        dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 607 ) THEN ! swc

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) )*  &
     &                        dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 612 ) THEN ! hwc

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )*  &
     &                        dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 613 ) THEN ! hlc

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) )*  &
     &                        dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 505 ) THEN ! rain

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j = 1,nr-1 ! nr1,nr2
                tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 506 .and. lice  ) THEN ! cloud ice (columns)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 507 .and. lice  ) THEN ! snow

          IF ( ls .gt. 1 ) THEN
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) )
            ENDDO
           ENDDO
          ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 508 .and. lgl .gt. 1 ) THEN ! QGL

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl) )
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( idx(iv) .eq. 509 .and. lgm .gt. 1 ) THEN ! QGM

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm) )
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( idx(iv) .eq. 510 .and. lgh .gt. 1 ) THEN ! QGH

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh) )
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( idx(iv) .eq. 511 .and. lf .gt. 1 ) THEN ! frozen drops QFD

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) )
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( idx(iv) .eq. 512 .and. lice ) THEN ! QH

          IF ( micro(1:3) == 'HCM' ) THEN
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO n = 0,nch-1
                tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh+n)
               ENDDO
              
              tem1(j,i,k) =  1.0e3*(  tmp )
            ENDDO
           ENDDO
          ENDDO

          ELSE
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
            ENDDO
           ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 513 .and. lhl .gt. 1 ) THEN ! QHL

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 514 .and. lir .gt. 1 ) THEN ! QIR

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) = 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lir) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 515 .and. lip .gt. 1 ) THEN ! QIP

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lip) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 516 .and. lsw .gt. 1 ) THEN ! QSW

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsw) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 517 .and. lhw .gt. 1 ) THEN ! QHW

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhw) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 518 .and. lhw .gt. 1 ) THEN ! QHI

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) - &
     &                               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhw) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 519 .and. lhlw .gt. 1 ) THEN ! QHW

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhlw) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 521  ) THEN ! FWS

          IF ( lsw > 1 ) THEN
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) > qxmin(ls) ) THEN
              tem1(j,i,k) = 100*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsw)/  &
     &                             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          ELSE
            tem1(:,:,:) = 0.0
          ENDIF

         ELSEIF ( idx(iv) .eq. 522  ) THEN ! FWH

          IF ( lhw > 1 ) THEN
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) > qxmin(lh) ) THEN
              tem1(j,i,k) = 100*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhw)/   &
     &                             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          ELSE
            tem1(:,:,:) = 0.0
          ENDIF

         ELSEIF ( idx(iv) .eq. 523  ) THEN ! FWHL

          IF ( lhlw > 1 ) THEN
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) > qxmin(lhl) ) THEN
              tem1(j,i,k) = 100*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhlw)/   &
     &                             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO
          ELSE
            tem1(:,:,:) = 0.0
          ENDIF

         ELSEIF ( idx(iv) .eq. 524 ) THEN ! rain, QRLOW

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j = 1,nr-1 ! nr1,nr2
                tem1(j,i,k) = Min( 0.02, 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) ) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 520 .and. lice ) THEN ! QHLLOW

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) > qxmin(lhl) ) THEN
                 tem1(j,i,k) =  Min( 0.02, 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) ))
               ELSE
                 tem1(j,i,k) = 0.
               ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 525 .and. lice ) THEN ! QHLOW

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =  Min( 0.02, 1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) ))
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 526 .and. lnr .gt. 1 .and. ipconc .ge. 3 ) THEN ! crwlow

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = Min(200.0,st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)) ! per m**3
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 527 .and. lnr .gt. 1 .and. ipconc .ge. 3 ) THEN ! ccilow

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = Min(500.0,1.e-3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni)) ! per L
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 528 .and. lns .gt. 1 .and. ipconc .ge. 4 ) THEN ! cswlow

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = Min(200.0,st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)) ! per m**3
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 529 .and. lice) THEN ! qitot 
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(           &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lir) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lip) )
              ELSEIF ( li .gt. 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(          &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) )
              
                IF ( lis > 1 ) THEN
                tem1(j,i,k) =    tem1(j,i,k) + 1.0e3*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lis) )
                ENDIF
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 530 .and. lice ) THEN ! qgtot 

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(           &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf))
              ELSEIF ( lh .gt. 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(           &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 531 .and. lice ) THEN ! qci (total cloud particles)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lir) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lip) )
              ELSEIF ( lc .gt. 1 .and. li .gt. 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 532 .and. lice ) THEN ! qcis (total cloud particles + snow)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lir) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lip) )
              ELSEIF ( lc .gt. 1 .and. li .gt. 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) )
              ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 537 ) THEN ! qtot (total hydro)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tem1(j,i,k) =          &
                     1.0e3*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lc) +  &
     &                       st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) )
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) = tem1(j,i,k) +         &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lip) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lir) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl))
            ELSEIF ( lice ) THEN
              tem1(j,i,k) =  tem1(j,i,k) +             &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) +                 &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
            ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 538 ) THEN ! qpre (total precip)

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgl) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgm) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lgh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl))
            ELSEIF ( lice ) THEN
              IF ( lhl > 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) )
              ELSE
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
              ENDIF
            ELSE
              tem1(j,i,k) =          &
                     1.0e3*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) )
            ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 540 .and. lice ) THEN ! qhtot

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(           &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) +         &
                      st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl))
            ELSEIF ( lh .gt. 1 ) THEN
              tem1(j,i,k) =          &
                     1.0e3*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) )
            ENDIF
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 557 .and. lvs .gt. 1 ) THEN  ! integrated snow volume
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                tem1(j,i,k) =  1.0e6*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvs) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 562 .and. lvh .gt. 1 ) THEN  ! integrated graupel volume
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                tem1(j,i,k) =  1.0e6*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) )
              ENDDO
            ENDDO
          ENDDO
         

         ELSEIF ( idx(iv) .eq. 563 .and. lvhl .gt. 1 ) THEN  ! integrated hail volume
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                tem1(j,i,k) =  1.0e6*(  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 564 .and. ( lzh .gt. 1 .or. micro(1:3) == 'TAK') ) THEN  ! graupel shape parameter ALPH
         
          IF ( micro(1:3) == 'TAK') THEN
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

               IF (  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) ) THEN

          ! get Z
              z = 0.0
               DO l = lf75+1,lfmax
                 z = z + szf(l,1)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh-1+l)
               ENDDO
              z = z*1.0e-12/rhof(1)**2  ! divide by rhof because it is contained in szf
          ! get concentration
              chw = 0.0
               DO l = lf75+1,lfmax
                 chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)
               hwdn = rho_qh_tak
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-9 .and.  z .gt. 0.0 ) THEN

            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           
           alphamintmp = alphamin_v5d ! -0.9
           alphamaxtmp = alphamax_v5d ! 25.
           
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamintmp, Min( alphamaxtmp, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamintmp, Min( alphamaxtmp, alp ) )
           ENDDO
                 
                 tem1(j,i,k) = alp

               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
            ENDDO
           ENDDO
          ENDDO
          
          ELSE
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzh) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzh)
                alpha1 = 0

              IF ( lvh .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*     &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/  &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ELSE
               hwdn = rho_qh
              ENDIF
                
            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
!           DO ii = 1,10
!            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
!             alpha1 = Max( alphamin, Min( alphamax, alp ) )
!             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!!           print*,'i,alp = ',i,alp
!             alp = Max( alphamin, Min( alphamax, alp ) )
!           ENDDO
                 
           alphamaxtmp = alphamax_v5d ! 60 ! alphamax
           alphamintmp = alphamin_v5d ! -0.9

           DO ii = 1,15
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamintmp, Min( alphamaxtmp, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamintmp, Min( alphamaxtmp, alp ) )
           ENDDO
           
                 tem1(j,i,k) = alp
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
              
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 565 .and.  lzhl < 1 .and. micro(1:1) == 'Z' ) THEN  ! hail shape parameter ALPHL
           ! constant value
             tem1(:,:,:) = alphahl

         ELSEIF ( idx(iv) .eq. 565 .and. ( lzhl .gt. 1 .or. micro(1:3) == 'TAK') ) THEN  ! hail shape parameter ALPHL

          IF ( micro(1:3) == 'TAK') THEN
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

               IF (  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) ) THEN

          ! get Z
              z = 0.0
               DO l = lf75+1,lfmax
                 z = z + szf(l,2)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl-1+l)
               ENDDO
              z = z*1.0e-12/rhof(2)**2   ! divide by rhof because it is contained in szf
          ! get concentration
              chw = 0.0
               DO l = lf75+1,lfmax
                 chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)
               hwdn = rho_qhl_tak
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-9 .and.  z .gt. 0.0 ) THEN

            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           alphamaxtmp = alphamax_v5d ! 60 ! alphamax
           alphamintmp = alphamin_v5d ! -0.9

           DO ii = 1,15
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamintmp, Min( alphamaxtmp, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamintmp, Min( alphamaxtmp, alp ) )
           ENDDO
                 
                 tem1(j,i,k) = alp

               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
            ENDDO
           ENDDO
          ENDDO
          
          
          ELSE
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl)
                alpha1 = 0

              IF ( lvhl .gt. 1 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
              ELSE
               hwdn = rho_qhl
              ENDIF
                
            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
!           DO ii = 1,10
!            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
!             alpha1 = Max( alphamin, Min( alphamax, alp ) )
!             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!!           print*,'i,alp = ',i,alp
!             alp = Max( alphamin, Min( alphamax, alp ) )
!           ENDDO
                 
           alphamaxtmp = alphamax_v5d ! 60 ! alphamax
           alphamintmp = alphamin_v5d ! -0.9

           DO ii = 1,15
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamintmp, Min( alphamaxtmp, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamintmp, Min( alphamaxtmp, alp ) )
           ENDDO                 
                 tem1(j,i,k) = alp
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
              
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 566 .and. lzr < 1 .and. micro(1:1) == 'Z' ) THEN  ! rain shape parameter ALPR

           tem1(:,:,:) = alphar
         
         ELSEIF ( idx(iv) .eq. 566 .and. ( lzr .gt. 1 .or. micro(1:3) == 'TAK') ) THEN  ! rain shape parameter ALPR
         
          IF ( microp(1:1) == 'Z' .and. imurain == 3 ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr) .gt. 0.0 ) THEN
                
                nrx = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr)
                alpha1 = 0

               hwdn = rho_qr
                
            vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr/(1000.*nrx)

!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
! determine shape parameter alpha by iteration
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha1+2.0)*nrx*vr**2/(z*pi1**2) - 1.
!           print*,'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z,xv
           DO ii = 1,20
!            IF ( 100.*Abs(alp - alpha1)/Abs(alpha1) .lt. 1. ) EXIT
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha1+2.0)*nrx*vr**2/(z*pi1**2) - 1.
!           print*,'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           ENDDO
                 
                 
                 tem1(j,i,k) = alp
               
               ELSE
                 tem1(j,i,k) = rnumin
               ENDIF
              
              
              ENDDO
            ENDDO
          ENDDO
          
          ELSEIF ( microp(1:3) == 'MY3' .or. ( microp(1:1) == 'Z' .and. imurain == 1 ) ) THEN

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
               
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. 1.e-9 .and.   &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) .and.    &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr) .gt. 0.0 ) THEN
                
                chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)
                qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)
                z = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr)
                alpha1 = 0

               hwdn = 1000.
                
            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           ENDDO
                 
                 
                 tem1(j,i,k) = alp
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
              
              ENDDO
            ENDDO
          ENDDO
          
          ELSEIF ( micro(1:3) == 'TAK' ) THEN ! alpr
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2

               IF (  st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. qxmin(lr) ) THEN

          ! get Z
              z = 0.0
               DO l = 15,ntakrd
!                 z = z + 1.d6*szr(l)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
                 z = z + szr(l)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
               ENDDO
              z = z*1.0e-12
          ! get concentration
              chw = 0.0
               DO l = 15,ntakrd
                 chw = chw + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
               ENDDO
               chw = 1.e6*chw ! convert from cm^-3 to m^-3
               
               qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)
               hwdn = 1000.
                alpha1 = 0
                
                
               IF ( chw .gt. 1.e-9 .and.  z .gt. 0.0 ) THEN

            rdi = z*(pi1/6.*hwdn)**2*chw/((dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2)

           alp = (6.0+alpha1)*(5.0+alpha1)*(4.0+alpha1)/ ((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           alphamaxtmp = alphamax
           DO ii = 1,10
            IF ( Abs(alp - alpha1) .lt. 0.01 ) EXIT
             alpha1 = Max( alphamin, Min( alphamaxtmp, alp ) )
             alp = (6.+alpha1)*(5.0+alpha1)*(4.0+alpha1)/((3.0+alpha1)*(2.0+alpha1)*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamaxtmp, alp ) )
           ENDDO
                 
                 tem1(j,i,k) = alp

               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
               
               ELSE
                 tem1(j,i,k) = -1.
               ENDIF
              
            ENDDO
           ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 570 .and. micro(1:3) == 'TAK' ) THEN ! ZRW for TAK

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO l = 24,ntakrd
                 tmp = tmp + 1.d6*szr(l)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr-1+l)
               ENDDO
              
              tem1(j,i,k) =  Max(-10., 10*Log10(Max(0.001, tmp ) ) )
            ENDDO
           ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 570 .or. idx(iv) == 578  ) THEN  ! rain reflectivity ZRW2
!             write(0,*) 'v5dout: ZRW lzr = ',lzr
         
          IF ( lzr .gt. 1 .and. idx(iv) .eq. 570 ) THEN ! predicted moment
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*(  & 
     &                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzr) ) ) ) )
              ENDDO
            ENDDO
          ENDDO
          
          ELSEIF ( lr .gt. 1 .and. lnr .gt. 1 ) THEN ! diagnosed using rnu (same as rnumin, usually)
          
          IF ( micro(1:1) == 'Z' ) THEN
!             write(0,*) 'v5dout: ZRW imurain = ',imurain
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr) .gt. 1.e-6 ) THEN
                  vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                 st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)/         &
     &                 (1000.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr))
               IF ( imurain == 3 ) THEN
                 z = 3.6*(rnu+2.0)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)*vr**2/(rnu+1.0)
               ELSEIF ( imurain == 1 ) THEN
                 g1 = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/  &
     &            ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
                 z = g1*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2   &
     &            *(6.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr)/(1000.*pi1))**2/  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr)
               
!                 IF ( i == 48 .and. j == 33 ) THEN
!                  write(0,*) 'v5dout: i,j,k, ZRW = ',i,j,k,10*Log10(Max(1.e-3,1.0e18*z))
!                  write(0,*) 'g1,qr,cr,alphar = ',g1,1.e3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr), &
!     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr),alphar
!                 ENDIF
               ENDIF

                tem1(j,i,k) = Max(-10., 10.*Log10(Max(1.e-3,1.0e18*z)) )
                
                ELSE 
                tem1(j,i,k) = -10.
                ENDIF
                
                
              ENDDO
            ENDDO
          ENDDO
          ELSEIF ( micro(1:3) == 'MY2' ) THEN
          ENDIF
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 571  ) THEN  ! snow reflectivity ZSW
         
          
          IF ( micro(1:3) == 'TAK' ) THEN
          
          tmpmax = 0.0
          qmax   = 0.0
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2

              tmp = 0.0
              DO ki=1,ntakit
                DO ii=1,ntakid
                   ia = lni + ii - 1 + (ki-1)*ntakid
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)*(((6.0/pii))*sxi(ii,ki))**2
                ENDDO
              ENDDO
              
              tmpmax = Max( tmpmax, tmp )
              qmax = Max( qmax, st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) )
              
              
              tem1(j,i,k) =  Max(-10., 10.0*Log10(Max(0.001, 1.e12*0.224*tmp ) ) )
          
              ENDDO
            ENDDO
          ENDDO
          
           write(0,*) 'TAK ZSW: tmpmax, qmax = ',tmpmax,qmax,Max(-10., 10.0*Log10(Max(0.001, 1.e12*0.224*tmpmax ) ) )
          
          
          ELSEIF ( micro(1:1) == 'Z' ) THEN
          
          
          IF ( ls .gt. 1 .and. lns .gt. 1 ) THEN ! diagnosed using snu
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns) .gt. 1.e-6 ) THEN
                  vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                 st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/         &
     &                 (100.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns))  ! 100. for snow density
               IF ( imusnow == 3 ) THEN
!                 z = 3.6*0.224*(snu+2.0)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*vr**2/(snu+1.0)*(100./1000.)**2

!                 z = 3.6*(snu+2.)*0.224*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)**2/ &
!     &           (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*(snu+1.)*1000.**2)*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2

!                 z = 1.06214**2*(0.189*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)**2)*  &
!    &                      dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2*gsnow1*gsnow73/    &
!     &                   (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*(917.)**2* gsnow53**2)

                 z = 323.3226*0.106214**2*(0.189*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)**2)*  &
    &                      dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2*gsnow73/    &
     &                   (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*(917.)**2 *gsnow1* (1.0+snu)**(4./3.))


!             tmp = Min(1.0,1.e3*(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls))*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!             z = Max( z, 1.e-18*750.0*(tmp)**1.98)

               ELSEIF ( imusnow == 1 ) THEN
                 g1 = (6.0 + alphas)*(5.0 + alphas)*(4.0 + alphas)/  &
     &            ((3.0 + alphas)*(2.0 + alphas)*(1.0 + alphas))
                 z = 0.224*g1*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2   &
     &            *(6.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/(1000.*pi1))**2/  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)
               
!                 IF ( i == 48 .and. j == 33 ) THEN
!                  write(0,*) 'v5dout: i,j,k, ZRW = ',i,j,k,10*Log10(Max(1.e-3,1.0e18*z))
!                  write(0,*) 'g1,qr,cr,alphas = ',g1,1.e3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr), &
!     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr),alphas
!                 ENDIF
               ENDIF

                tem1(j,i,k) = Max(-10., 10*Log10(Max(1.e-3,1.0e18*z)) )
                
                ELSE 
                tem1(j,i,k) = -10.
                ENDIF
                
                
              ENDDO
            ENDDO
          ENDDO
          ENDIF
          ENDIF

         ELSEIF ( idx(iv) .eq. 572  ) THEN  ! test snow reflectivity ZSW (listed as ZGL)
         
          
          IF ( ls .gt. 1 .and. lns .gt. 1 ) THEN ! diagnosed using snu
          
          IF ( micro(1:1) == 'Z' ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns) .gt. 1.e-6 ) THEN
                  vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*            &
     &                 st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/         &
     &                 (100.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns))  ! 100. for snow density
               IF ( imusnow == 3 ) THEN
!                 z = 3.6*0.224*(snu+2.0)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*vr**2/(snu+1.0)*(100./1000.)**2

                 z = 3.6*(snu+2.)*0.224*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)**2/ &
     &           (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*(snu+1.)*1000.**2)*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2

!                 z = 323.3226*0.106214**2*(0.189*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)**2)*  &
!    &                      dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2*gsnow73/    &
!     &                   (st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)*(917.)**2 *gsnow1* (1.0+snu)**(4./3.))

! 323.3226* 0.106214**2*(ksq*an(ix,jy,kz,ls) + (1.-ksq)*qxw)*an(ix,jy,kz,ls)*db(ix,jy,kz)**2*gsnow73/    &
!     &                   (an(ix,jy,kz,lns)*(917.)**2* gsnow1*(1.0+snu)**(4./3.))

!             tmp = Min(1.0,1.e3*(st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls))*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!             z = Max( z, 1.e-18*750.0*(tmp)**1.98)

               ELSEIF ( imusnow == 1 ) THEN
                 g1 = (6.0 + alphas)*(5.0 + alphas)*(4.0 + alphas)/  &
     &            ((3.0 + alphas)*(2.0 + alphas)*(1.0 + alphas))
                 z = 0.224*g1*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2   &
     &            *(6.*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/(1000.*pi1))**2/  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lns)
               
!                 IF ( i == 48 .and. j == 33 ) THEN
!                  write(0,*) 'v5dout: i,j,k, ZRW = ',i,j,k,10*Log10(Max(1.e-3,1.0e18*z))
!                  write(0,*) 'g1,qr,cr,alphas = ',g1,1.e3*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lr), &
!     &             st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnr),alphas
!                 ENDIF
               ENDIF

                tem1(j,i,k) = Max(-10., 10*Log10(Max(1.e-3,1.0e18*z)) )
                
                ELSE 
                tem1(j,i,k) = -10.
                ENDIF
                
                
              ENDDO
            ENDDO
          ENDDO
          ENDIF
          ENDIF


         ELSEIF ( idx(iv) .eq. 576 .and. micro(1:3) == 'HCM' ) THEN ! ZHW for HCM

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO l = 1,nch
                 tmp = tmp + 0.224*(1000.*rdamelt(l))**6*dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)* &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh-1+l)/hm(l)
               ENDDO
              
              tem1(j,i,k) =  Max(-10., 10*Log10(Max(0.001, tmp ) ) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 576 .and. micro(1:3) == 'TAK' ) THEN ! ZHW for TAK

          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO l = ls250+1,ntakpd
                 tmp = tmp + 0.224*1.d6*szf(l,1)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh-1+l)
               ENDDO
              
              tem1(j,i,k) =  Max(-10., 10*Log10(Max(0.001, tmp ) ) )
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 577 .and. micro(1:3) == 'TAK' ) THEN ! ZHL for TAK

!          write(0,*) 'V5DOUT: ZHL for TAK, nz_v5d,lhl = ',nz_v5d,lhl
          
          DO k=1,nz_v5d-1
            DO i=2-ibc,nc-2+ibc
              DO j=nr1,nr2
              tmp = 0.0
               DO l = ls250+1,ntakpd
                 tmp = tmp + 0.224*1.d6*szf(l,2)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl-1+l)
               ENDDO
              
              tem1(j,i,k) =  Max(-10., 10*Log10(Max(0.001, tmp ) ) )
            ENDDO
           ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 575 .and. micro(1:3) /= 'HCM' ) THEN  ! frozen drops 6th moment
         
          IF ( lzf .gt. 1 ) THEN

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
              IF ( lvf .gt. 1 ) THEN
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) > qxmin(lf) .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvf) > 1.0e-28 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvf)
               ELSE
                hwdn = rho_qh
               ENDIF
              ELSE
               hwdn = rho_qh
              ENDIF
              
                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(  & 
     &                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzf)*(hwdn/1000.)**2  ) ) ) )

              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 576 .and. micro(1:3) /= 'HCM' ) THEN  ! graupel 6th moment
         
          IF ( lzh .gt. 1 ) THEN

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
              IF ( lvh .gt. 1 ) THEN
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) > qxmin(lh) .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) > 1.0e-28 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
               ELSE
                hwdn = rho_qh
               ENDIF
              ELSE
               hwdn = rho_qh
              ENDIF
              
                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(  & 
     &                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzh)*(hwdn/1000.)**2  ) ) ) )

              ENDDO
            ENDDO
          ENDDO
          
          ELSEIF ( lnh > 1 ) THEN

          g1 = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh) .gt. 1.e-15 ) THEN

              IF ( lvh .gt. 1 .and. st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) > 0.0) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
              ELSE
               hwdn = rho_qh
              ENDIF

             chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnh)
             qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) 
             zx = g1*0.224*(dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2/chw
             ze =1.e18*zx*(6./(pii*hwdn))**2

                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, ze ) ) )
                
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 577 .and. lzhl .gt. 1 ) THEN  ! hail 6th moment
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2

              IF ( lvhl .gt. 1 ) THEN
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) > qxmin(lhl) .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl) > 1.0e-28 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
               ELSE
                hwdn = rho_qhl
               ENDIF
              ELSE
               hwdn = rho_qhl
              ENDIF
              
                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(  & 
     &                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl)*(hwdn/1000.)**2  ) ) ) )


!                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, 1.0e18*0.224*(   &
!     &                            st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lzhl)) ) ) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 577 .and. lnhl > 1 .and. lzhl .lt. 1 ) THEN  ! hail 6th moment
         
          g1 = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl) .gt. 1.e-15 ) THEN

              IF ( lvhl .gt. 1 .and. st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl) > 0.0 ) THEN
               hwdn = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
              ELSE
               hwdn = rho_qhl
              ENDIF
             chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lnhl)
             qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) 
             zx = g1*0.224*(dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2/chw
             ze =1.e18*zx*(6./(pii*hwdn))**2

                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, ze ) ) )
                
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          
         ELSEIF ( idx(iv) .eq. 579 .and. micro(1:3) == 'TAK'  ) THEN  ! ice crystal reflectivity

          tmpmax = 0.0
          qmax   = 0.0
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2

              tmp = 0.0
              DO ki = 1,2 !ntakit
                DO ii = 1,ntakid
                   ia = lni + ii - 1 + (ki-1)*ntakid
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)*(((6.0/pii))*sxi(ii,ki))**2
                ENDDO
              ENDDO
              
              tmpmax = Max( tmpmax, tmp )
              qmax = Max( qmax, st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) )
              
              
              tem1(j,i,k) =  Max(-10., 10.0*Log10(Max(0.001, 1.e12*0.224*tmp ) ) )
          
              ENDDO
            ENDDO
          ENDDO
          
           write(0,*) 'TAK ZCI: tmpmax, qmax = ',tmpmax,qmax,Max(-10., 10.0*Log10(Max(0.001, 1.e12*0.224*tmpmax ) ) )


         ELSEIF ( idx(iv) .eq. 579 .and. micro(1:1) == 'Z' .and. lni > 1 ) THEN  ! ice crystal reflectivity
         


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
                IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li) .gt. 1.e-9 .and.   &
     &               st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni) .gt. 1.e-15 ) THEN

               hwdn = 900.

             chw = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lni)
             qr = st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,li)
             zx = g1*0.224*(dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr)**2/chw
             ze =1.e18*zx*(6./(pii*hwdn))**2

                 vr = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*qr/(900.*chw)
                 ze =  0.224*3.6e18*(cinu+2.)*chw*vr**2/(cinu+1.)*(900./1000.)**2


                tem1(j,i,k) = Max(-10., 10*Log10(Max(0.001, ze ) ) )
                
                ELSE
                tem1(j,i,k) = -10.
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          

         ELSEIF ( idx(iv) .eq. 587 .and. lvs .gt. 1 ) THEN  ! average snow density
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvs) .gt. 1.e-20 .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls) .gt. qxmin(ls) ) THEN
                 tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*      &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/   &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvs)
!                 tem1(j,i,k) = Min(300.0, Max(50.0,dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ls)/st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvs) ))
               ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 591 .and. lvf .gt. 1 ) THEN  ! average frozen drops density DNF
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvf) .gt. 1.e-20 .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf) .gt. qxmin(lf) ) THEN
                 tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*      &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lf)/   &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvf)
!                 tem1(j,i,k) = Min(900.0, Max(100.0,dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) ))
               ENDIF
              ENDDO
            ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 592 .and. lvh .gt. 1 ) THEN  ! average graupel density DNH
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) .gt. 1.e-20 .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) .gt. qxmin(lh) ) THEN
                 tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*      &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/   &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh)
!                 tem1(j,i,k) = Min(900.0, Max(100.0,dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh)/st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvh) ))
               ENDIF
              ENDDO
            ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 593 .and. lvhl .gt. 1 ) THEN  ! average hail density DNHL
         
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=nr1,nr2
               IF ( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl) .gt. 1.e-20 .and.  &
     &              st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl) .gt. qxmin(lhl) ) THEN
                 tem1(j,i,k) = dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*       &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lhl)/   &
     &                         st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lvhl)
               ENDIF
              ENDDO
            ENDDO
          ENDDO
         
         
         ELSEIF ( idx(iv) .eq. 50 ) THEN ! composite reflectivity (DBZCMP)

            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,1) = 0.0
                DO k=1,nz_v5d-1
                  tem1(j,i,1) = Max( tem1(j,i,1), dbz%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) )
                ENDDO
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 22 .or. idx(iv) .eq. 51) THEN ! reflectivity

          DO k=1,nz_v5d-1 ! (nl(iv))*v5dstridez
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = dbz%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO
         
! electrifical variables
         ELSEIF ( idx(iv) .eq. 25 .and.  lscpi .gt. 1 ) THEN ! CIONP

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-6*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscpi)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 26  .and.  lscni .gt. 1  ) THEN ! CIONN

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-6*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscni)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 250 .and.  lscpli .gt. 1 ) THEN ! CLIONP

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-6*st(i,nr-j,k,lscpli)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 260  .and.  lscnli .gt. 1  ) THEN ! CLIONN

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-6*st(i,nr-j,k,lscnli)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 251  .and.  lscpi .gt. 1 ) THEN ! SCIONP

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*ec*(  st(i,nr-j,k,lscpi)  )
                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCIONP at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpi,lscni = ',st(i,nr-j,k,lscpi), st(i,nr-j,k,lscni)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 252  .and.  lscpli .gt. 1 ) THEN ! SCLIONP

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*ec*(  st(i,nr-j,k,lscpli)  )
                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCLIONP at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpi,lscni = ',st(i,nr-j,k,lscpi), st(i,nr-j,k,lscni)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 261  .and.  lscni .gt. 1 ) THEN ! SCIONN

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*(-1)*ec*( st(i,nr-j,k,lscni) )
                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCLIONN at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpi,lscni = ',st(i,nr-j,k,lscpi), st(i,nr-j,k,lscni)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 262  .and.  lscnli .gt. 1 ) THEN ! SCLIONN

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*(-1)*ec*( st(i,nr-j,k,lscnli) )
                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCLIONN at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpi,lscni = ',st(i,nr-j,k,lscpi), st(i,nr-j,k,lscni)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 27  .and.  lscpi .gt. 1 ) THEN ! SCION

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                IF ( lscpli .gt. 1 ) THEN
                tem1(j,i,k) = 1.0e9*ec*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscpi)    &
     &                                 - st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscni)    &
     &                                 + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscpli)   &
     &                                 - st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscnli) )
                ELSE
                tem1(j,i,k) = 1.0e9*ec*( st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscpi) -    &
     &                                   st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscni) )
                ENDIF

                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCION at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpi,lscni = ',st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscpi),  &
     &                                        st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscni)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 270  .and.  lscpli .gt. 1 ) THEN ! SCLION

          DO k=1,nz-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*ec*( st(i,nr-j,k,lscpli) - st(i,nr-j,k,lscnli) )
                IF ( Abs( tem1(j,i,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCLION at i,j,k =',i,j,k,tem1(j,i,k)
                  write(0,*) 'lscpli,lscnli = ',st(i,nr-j,k,lscpli), st(i,nr-j,k,lscnli)
                  tem1(j,i,k) = 0.0
                ENDIF
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 30 ) THEN ! GLM_8KM

          fac = 1.0

            IF ( iflshr8km > 0 ) THEN
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(iflshr8km)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
            ELSE
             write(0,*) 'V5DOUT: error for GLM_8KM, iflshr8km = 0'
            ENDIF

         ELSEIF ( idx(iv) .eq. 31 ) THEN ! native FED

          fac = 1.0

            IF ( iflshfed > 0 ) THEN
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(iflshfed)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
            ELSE
             write(0,*) 'V5DOUT: error for FLSHFED, iflshfed = 0'
            ENDIF



         ELSEIF ( idx(iv) .eq. 33 .and. ipelec .gt. 0 ) THEN ! FLSH 

          fac = 1.0

          IF ( allocated(flshp) ) THEN
            DO k=1,nz_v5d-1
              DO i=1,nc-1
                DO j=1,nr-1
                  tem1(j,i,k) = Max(-segmx, flshn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)) +  &
                                Min( segmx, flshp(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO k=1,nz_v5d-1
              DO i=1,nc-1
                DO j=1,nr-1
                  tem1(j,i,k) = Max(-segmx, elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)) +  &
                                Min( segmx, elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                ENDDO
              ENDDO
            ENDDO
          ENDIF
         
         ELSEIF ( idx(iv) .eq. 18 .and. ipelec .gt. 0 ) THEN ! FLSHI 

          fac = 1.0

          IF ( allocated(flshi) ) THEN
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                tem1(j,i,k) = flshi(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
               ENDDO
             ENDDO
           ENDDO
          ELSE
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                 tem1(j,i,k) = elec(ieinit)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 19 .and. ipelec .gt. 0 ) THEN ! FLSHP

          fac = 1.0

          IF ( allocated(flshp) ) THEN
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                tem1(j,i,k) = Min( segmx, flshp(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) )
               ENDDO
             ENDDO
           ENDDO
          ELSE
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                 tem1(j,i,k) = Min( segmx, elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 20 .and. ipelec .gt. 0 ) THEN ! FLSHN

          fac = 1.0
          
          IF ( allocated(flshn) ) THEN
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                tem1(j,i,k) = Max( -segmx, flshn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) )
               ENDDO
             ENDDO
           ENDDO
          ELSE
           DO k=1,nz_v5d-1
             DO i=1,nc-1
               DO j=1,nr-1
                 tem1(j,i,k) = Max(-segmx, elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 34 .and. associated( vzf%flt3d ) ) THEN ! reflectivity

          DO k=1,nz_v5d-1 ! (nl(iv))*v5dstridez
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = vzf%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 35 ) THEN ! horizontal wind speed

          DO k=1,nz_v5d-1 ! (nl(iv))*v5dstridez
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = Sqrt(0.25*(u%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
     &                             u%flt3d(i+is+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))**2 + &
     &                         0.25*(v%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +        &
     &                             v%flt3d(i+ixb_v5d-1,nr-j+js+jyb_v5d-1,k+kzb_v5d-1))**2 )
            ENDDO
          ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 60 ) THEN ! IWP ice water path (graupel) (above -10 C level)

            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,1) = 0.0
                DO k=1,nz_v5d-1
                  IF ( tt0(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) .lt. 263.15 ) THEN
                  tem1(j,i,1) = tem1(j,i,1) + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lh) *  &
     &                          dn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)/gz(3)%flt1d(k)
                  ENDIF
                ENDDO
                
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 58 ) THEN ! HAILRAT 

            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(2)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 59 ) THEN ! HAILACC

            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(4)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 61 ) THEN ! FLASHR FLSHR 

          fac = 1.0

            IF ( iflshr > 0 ) THEN
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(iflshr)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
            ELSE
             write(0,*) 'V5DOUT: error for FLSHR, iflshr = 0'
            ENDIF

         ELSEIF ( idx(iv) .eq. 62 ) THEN ! RAINRAT 

            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(1)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 63 ) THEN ! RAINACC

            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,1) = precip(3)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO

!!!!!! surface physics variables

        ELSEIF ( idx(iv) .eq. 65 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(tcanp)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 66 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(wcanp)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 67 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(qav)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 68 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(veg)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 69 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(stype)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 70 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(tsrfc)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 71 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(wsfc)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 72 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(tsoil)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 73 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(wsoil)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 74 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(vlai)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 75 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(albedo)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 76 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(rough)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 78 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(eflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 79 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(fflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 80 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(uflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 81 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(vflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 82 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(tflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 83 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(qflx)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)*2500000.0
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 84 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(radsw)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 85 .and. isfcphys > 0 ) THEN

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = gd%var(radlw)%flt2d(i+ixb_v5d-1,nr-j+jyb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 86  ) THEN ! surface pressure pert.

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

        ELSEIF ( idx(iv) .eq. 87  ) THEN ! surface pressure

          DO k=1,1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e-2*( pn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) + pb(k) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( ( idx(iv) .eq. 118 .or. idx(iv) .eq. 121 .or. idx(iv) .eq. 122 .or. idx(iv) .eq. 123 ) .and. ipelec .gt. 0 ) THEN ! FLSHDEN

!           tem1(:,:,:) = 0
! project to ground:
             DO i=1,nc-1
               DO j=1,nr-1
               
               sum = 0.0
               
                IF ( idx(iv) .eq. 118 ) THEN
                  DO k=1,nz_v5d-1
                   sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                            Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ENDDO
                ELSEIF ( idx(iv) .eq. 121 ) THEN
                  DO k=1,nz_v5d-1
                   sum = sum +  Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ENDDO
                
                ELSEIF ( idx(iv) .eq. 122 ) THEN
                  DO k=1,nz_v5d-1
                   sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                  ENDDO
                ELSEIF ( idx(iv) .eq. 123 ) THEN
                  DO k=1,nz_v5d-1
                   sum = sum + 0.1*elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                            Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ENDDO
                
                ENDIF
                
                IF ( sum .ge. 1.0 ) THEN  
                 tem1(j,i,1) = Log(sum+1)
                 tem1(j,i,nz_v5d-1) = Log(sum+1)
!                 tem1(j,i,nz_v5d) = Log(sum+1)
                ELSE
                 tem1(j,i,1) = 0.0
                 tem1(j,i,nz_v5d-1) = 0.0
                ENDIF
                sum = 0
              ENDDO
            ENDDO

! project to north/south:
             DO i=1,nc-1
                DO k=1,nz_v5d-1
               
               sum = 0.0

                DO j=1,nr-1
               
                  IF ( idx(iv) .eq. 118 ) THEN
                    sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                               Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ELSEIF ( idx(iv) .eq. 121 ) THEN
                    sum = sum +  Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ELSEIF ( idx(iv) .eq. 122 ) THEN
                    sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) 
                  ELSEIF ( idx(iv) .eq. 123 ) THEN
                    sum = sum + 0.1*elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                               Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ENDIF
                ENDDO
                
                IF ( sum .ge. 1.0 ) THEN  
                 tem1(1,i,k) = Log(sum+1)  ! note that j=1 is the northern plane
                 tem1(nr-1,i,k) = Log(sum+1)  ! note that j=1 is the northern plane
                ELSE
                 tem1(1,i,k) = 0.0
                ENDIF
                sum = 0
              ENDDO
            ENDDO

! project to east/west:
               DO j=1,nr-1
                DO k=1,nz_v5d-1
               
               sum = 0.0
               
                 DO i=1,nc-1
                  IF ( idx(iv) .eq. 118 ) THEN
                    sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                            Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ELSEIF ( idx(iv) .eq. 121 ) THEN
                    sum = sum +  Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ELSEIF ( idx(iv) .eq. 122 ) THEN
                    sum = sum + elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) 
                  ELSEIF ( idx(iv) .eq. 123 ) THEN
                    sum = sum + 0.1*elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
                            Abs(elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
                  ENDIF
                ENDDO
                
                IF ( sum .ge. 1.0 ) THEN  
                 tem1(j,nc-1,k) = Log(sum+1)
                 tem1(j,1,k) = Log(sum+1)
                ELSE
                 tem1(j,nc-1,k) = 0.0
                ENDIF
                sum = 0
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 119 .and. ipelec .gt. 0 ) THEN ! FLSHPT with cloud-top cut-off

          fac = 1.0
          
          IF ( icldtop .eq. -1 ) CALL makeflshmask(icldtop,nx,ny,nz,ng,st,na,mask,micro)
          
          IF ( allocated(flshp) .and. icldtop .gt. 1 ) THEN
             DO i=1,nc-1
               DO j=1,nr-1
                DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                tem1(j,i,k) = Min( segmx, flshp(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) )
               ENDDO
             ENDDO
           ENDDO
          ELSEIF ( icldtop .gt. 1 ) THEN
           DO i=1,nc-1
             DO j=1,nr-1
               DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                 tem1(j,i,k) = Min( segmx, elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 120 .and. ipelec .gt. 0 ) THEN ! FLSHNT with cloud-top cut-off

          fac = 1.0
          
          IF ( icldtop .eq. -1 ) CALL makeflshmask(icldtop,nx,ny,nz,ng,st,na,mask,micro)
          
          IF ( allocated(flshn) .and. icldtop .gt. 1 ) THEN
             DO i=1,nc-1
               DO j=1,nr-1
                DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                tem1(j,i,k) = Max( -segmx, flshn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) )
               ENDDO
             ENDDO
           ENDDO
          ELSEIF ( icldtop .gt. 1 ) THEN
           DO i=1,nc-1
             DO j=1,nr-1
               DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                 tem1(j,i,k) = Max( -segmx, elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 133 .and. ipelec .gt. 0 ) THEN ! FLSHT with cloud-top cut-off

          fac = 1.0
          
          IF ( icldtop .eq. -1 ) CALL makeflshmask(icldtop,nx,ny,nz,ng,st,na,mask,micro)
          
          IF ( allocated(flshn) .and. icldtop .gt. 1 ) THEN
             DO i=1,nc-1
               DO j=1,nr-1
                DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                  tem1(j,i,k) = Max(-segmx, flshn(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)) +  &
                                Min( segmx, flshp(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ELSEIF ( icldtop .gt. 1 ) THEN
           DO i=1,nc-1
             DO j=1,nr-1
               DO k=1,mask(i+ixb_v5d-1,nr-j+jyb_v5d-1)
                  tem1(j,i,k) = Max(-segmx, elec(ieflshn)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)) +  &
                                Min( segmx, elec(ieflshp)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
               ENDDO
             ENDDO
           ENDDO
          ENDIF

         ELSEIF ( idx(iv) .eq. 399 .and. lscw .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCCW


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscw)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 400 .and. lsci .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCCI

          IF ( micro(1:3) /= 'TAK' ) THEN
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsci)
              ENDDO
            ENDDO
          ENDDO

          ELSE

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
                 DO ki=1,ntakit
                  DO ii=1,ntakid
                    ia = lsci + ii - 1 + (ki-1)*ntakid
                    tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,ia)
                  ENDDO
                 ENDDO
                tem1(j,i,k) = 1.e9*tmp
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 403 .and. lscr .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCRW


          IF ( micro(1:3) /= 'TAK' ) THEN
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscr)
              ENDDO
            ENDDO
          ENDDO
          
          ELSE

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
               DO l = 1,ntakrd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscr-1+l)
               ENDDO
                tem1(j,i,k) = 1.e9*tmp
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF



         ELSEIF ( idx(iv) .eq. 404 .and. lscs .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCS


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscs)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 405 .and. lscgl .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCGL


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscgl)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 406 .and. lscgm .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCGM


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscgm)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 407 .and. lscgh .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCGH

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscgh)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 408 .and. lscf .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCFD


          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lscf)
              ENDDO
            ENDDO
          ENDDO


         ELSEIF ( idx(iv) .eq. 409 .and. lsch .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCH


          IF ( micro(1:3) /= 'TAK' ) THEN
          
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsch)
              ENDDO
            ENDDO
          ENDDO
          
          ELSE

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
               DO l = 1,ntakpd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lsch-1+l)
               ENDDO
                tem1(j,i,k) = 1.e9*tmp
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

         ELSEIF ( idx(iv) .eq. 396 .and. lschl .gt. 1 .and. ipelec .gt. 0 ) THEN ! SCHL


          IF ( micro(1:3) /= 'TAK' ) THEN

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e9*st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lschl)
              ENDDO
            ENDDO
          ENDDO

          ELSE

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tmp = 0.0
               DO l = 1,ntakpd
                 tmp = tmp + st(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1,lschl-1+l)
               ENDDO
                tem1(j,i,k) = 1.e9*tmp
              ENDDO
            ENDDO
          ENDDO
          
          ENDIF

        ELSEIF ( (idx(iv) .eq. 410 .or. idx(iv) .eq. 393) .and. ipelec .gt. 0 ) THEN ! scnet or energy


         DO k = 1, nz-1
          DO j = 1, ny-1
           DO i = 1, nx-1
          
             if ( lscpli .gt. 1 ) then
              tt7(i,j,k) = ec*(  st(i,j,k,lscpi)  - st(i,j,k,lscni)  &
     &                         + st(i,j,k,lscpli) - st(i,j,k,lscnli) )  
             else
             tt7(i,j,k) = ec*(st(i,j,k,lscpi) - st(i,j,k,lscni))  
             endif

                IF ( Abs( tt7(i,j,k) ) .gt. 20.0 ) THEN
                  write(0,*) 'problem with SCION (tt7) at i,j,k =',i,j,k,tt7(i,j,k)
                  write(0,*) 'lscpi,lscni = ',st(i,j,k,lscpi), st(i,j,k,lscni)
                  tt7(i,j,k) = 0.0
                ENDIF
      
             DO ia = lscb,lsceq
               tt7(i,j,k) = tt7(i,j,k) + st(i,j,k,ia)
             ENDDO

           ENDDO
          ENDDO
         ENDDO

          IF ( idx(iv) == 410 ) THEN ! net charge
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*tt7(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

          ELSEIF (idx(iv) == 392 ) THEN ! electrical energy  from Emag^2 only

           DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                dv = dxx(i+ixb_v5d-1)*dyy(nr-j+jyb_v5d-1)*dzz(k+kzb_v5d-1)

                tem1(j,i,k) = 1.e-6*0.5*eperao*elec(iemag)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2*dv 
              ENDDO
            ENDDO
          ENDDO
          
          ELSEIF (idx(iv) == 393 ) THEN ! electrical energy
           DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                dv = dxx(i+ixb_v5d-1)*dyy(nr-j+jyb_v5d-1)*dzz(k+kzb_v5d-1)
                dax = dyy(nr-j+jyb_v5d-1)*dzz(k+kzb_v5d-1)
                day = dxx(i+ixb_v5d-1)*dzz(k+kzb_v5d-1)
                daz = dxx(i+ixb_v5d-1)*dyy(nr-j+jyb_v5d-1)

!                phiw = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-2,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!                phie = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d  ,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!                phis = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-2,k+kzb_v5d-1) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!                phin = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d  ,k+kzb_v5d-1))
!                phid = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-2) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
!                phiu = 0.5*(elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
!                            elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d  ))

                phiw =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                phie =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                phis =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                phin =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                phid =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
                phiu =  elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)

               IF ( .false. ) THEN ! old method
               tem1(j,i,k) = 1.e-6*0.5*elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)* &
                tt7(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*dv
               ELSE
               tem1(j,i,k) = 1.e-6*0.5*eperao*(elec(iemag)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)**2*dv + &
                  ( -elec(iex)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*phiw*dax    &
                    +elec(iex)%flt3d(i+ixb_v5d  ,nr-j+jyb_v5d-1,k+kzb_v5d-1)*phie*dax    &
                    -elec(iey)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*phis*day    &
                    +elec(iey)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d  ,k+kzb_v5d-1)*phin*day    &
                    -elec(iez)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)*phid*dax    &
                    +elec(iez)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d  )*phiu*daz    &
                  ) )
               ENDIF
              ENDDO
            ENDDO
          ENDDO
         
          ENDIF


        ELSEIF ( idx(iv) .eq. 418 .and. ipelec .gt. 0 ) THEN ! scnet without ion charge

         tt7(:,:,:) = 0.0

         DO k = 1, nz-1
          DO j = 1, ny-1
           DO i = 1, nx-1
          
      
             DO ia = lscb,lsceq
               tt7(i,j,k) = tt7(i,j,k) + st(i,j,k,ia)
             ENDDO

           ENDDO
          ENDDO
         ENDDO

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.0e9*tt7(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 411 .and. ipelec .gt. 0 ) THEN ! Ex

          is = elec(iex)%istag

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.001*0.5*(elec(iex)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +  &
     &                                   elec(iex)%flt3d(i+is+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1))
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 412 .and. ipelec .gt. 0 ) THEN ! Ey

          js = elec(iey)%jstag

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.001*0.5*(elec(iey)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +   &
     &                                   elec(iey)%flt3d(i+ixb_v5d-1,nr-j+js+jyb_v5d-1,k+kzb_v5d-1))
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 413 .and. ipelec .gt. 0 ) THEN ! Ez

          ks = elec(iez)%kstag

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 0.001*0.5*(elec(iez)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1) +    &
     &                                   elec(iez)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1+ks) )
              ENDDO
            ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 414 .and. ipelec .gt. 0 ) THEN ! emag

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-3*elec(iemag)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO
         
         ELSEIF ( idx(iv) .eq. 415 .and. ipelec .gt. 0 ) THEN ! Potential

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
                tem1(j,i,k) = 1.e-6*elec(ipot)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              ENDDO
            ENDDO
          ENDDO

          ELSEIF ( idx(iv) .eq. 41 .and. ipelec .gt. 0 ) THEN ! ezsfc 

          ks = elec(iez)%kstag

            k = 1
            DO i=1,nc-1
              DO j=1,nr-1
               IF ( ks .eq. 1 ) THEN 
                 tem1(j,i,1) = elec(iez)%flt3d(i,nr-j,1)
               ELSE
                 tem1(j,i,1) =  -0.001*(elec(ipot)%flt3d(i,nr-j,1) )*gz(1)%flt1d(1)
               ENDIF
              ENDDO
            ENDDO

         ELSEIF ( idx(iv) .eq. 856 .and. ipelec .ge. 2 ) THEN ! noninductive charging

          tmp = 0.0
          tmpn = 0.0
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e12*elec(icghis)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
              tmp  = Max(tmp,  tem1(j,i,k) )
              tmpn = Min(tmpn, tem1(j,i,k) )
            ENDDO
           ENDDO
          ENDDO
                    

         ELSEIF ( idx(iv) .eq. 857 .and. ipelec .ge. 2  ) THEN ! inductive charging

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e12*elec(icghw)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 871 .and. ipelec .ge. 2 .and. icgaddl > 1 ) THEN ! "alternative" charging

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e12*elec(icgaddl)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO

         ELSEIF ( idx(iv) .eq. 872 .and. ipelec .ge. 2 .and. icgnoliq > 1 ) THEN ! liq-free noninductive charging

          DO k=1,nz_v5d-1
            DO i=1,nc-1
              DO j=1,nr-1
              tem1(j,i,k) = 1.0e12*elec(icgnoliq)%flt3d(i+ixb_v5d-1,nr-j+jyb_v5d-1,k+kzb_v5d-1)
            ENDDO
           ENDDO
          ENDDO



         ENDIF ! idx check
         
!
! Copy slab to j=2 for 2D simulation because vis5d requires a minimum of 2 rows.
!
         IF ( ny .eq. 2 ) THEN
          DO k=1,nz_v5d-1
            DO i=1,nc-1
              tem1(2,i,k) =  tem1(1,i,k)
            ENDDO
          ENDDO

          IF ( onedoutput /= 0 ) THEN
             fac1d = 1.0
             IF ( units(iv) == 'micron' ) fac1d = 0.001 ! convert microns to mm
             DO k = 1,nz_v5d-1
               data_1d(iv,k) = tem1(1,Abs(onedoutput),k)
             ENDDO
          ENDIF

         ENDIF

           itimes0 = itimes(it)
           idates0 = idates(it)
           
!        n = v5dwriteappend(it,iv,tem1,itimes(it),idates(it))
        IF ( vis5dstridex == 1 .and. vis5dstridez == 1 ) THEN
        n = v5dwriteappend(it,iv,tem1,itimes0,idates0)
        
        ELSE
         tem2(:,:,:) = 0.0
!         fac = 1./Float(vis5dstridex*vis5dstridex*vis5dstridez)
         DO k=1,(nz_v5d-1)/vis5dstridez
            DO i=1,Max(1,(nc-1)/vis5dstridex)
              DO j=1,Max(1, (nr-1)/vis5dstridex )
              
              smax = 0.0
              smin = 0.0
              DO k1 = 1,vis5dstridez
              DO i1 = 1,Min(nx-1,vis5dstridex)
              DO j1 = 1,Min(ny-1,vis5dstridex)
               ix = (i-1)*vis5dstridex + i1
               jy = (j-1)*vis5dstridex + j1
               kz = (k-1)*vis5dstridez + k1
               smin = Min( smin, tem1(jy,ix,kz) )
               smax = Max( smax, tem1(jy,ix,kz) )
               tem2(j,i,k) = tem2(j,i,k) + tem1(jy,ix,kz)
              ENDDO
              ENDDO
              ENDDO
              tem2(j,i,k) = Max( smin, Min( smax, fac * tem2(j,i,k) ) )
            ENDDO
           ENDDO
          ENDDO
        
         IF ( ny .eq. 2 ) THEN
          DO k=1,(nz_v5d-1)/vis5dstridez
            DO i=1,(nc-1)/vis5dstridex
              tem2(2,i,k) =  tem2(1,i,k)
            ENDDO
          ENDDO
         
         ENDIF

        n = v5dwriteappend(it,iv,tem2,itimes0,idates0)
        
!        write(0,*) 'writeappend: n,iv,idx(iv) = ',
!     :      n,iv,idx(iv),itimes(it),idates(it)
        ENDIF
        
        
        ENDDO ! iv


         ! write 1d data to file
          IF ( onedoutput /= 0 ) THEN
             DO k = 1,nz_v5d-1
              write(33,12) time,gz(1)%flt1d(k)/1000., ( data_1d(iv,k), iv=1,numvars )
               
             ENDDO
             
             deallocate( data_1d )
          ENDIF

!        n = v5dwrite(it,iv,tem1) ! ,itimes(1),idates(1))
!          write(0,'(a,i3)') 'Writing the Vis5D file, n= ',n

         n = myv5dupdate()
         
!          write(0,'(a,i3)') 'Updating the Vis5D file, n= ',n
!          write(6,'(a,i3)') 'Updating the Vis5D file, n= ',n
        
      
          IF ( allocated(flshp) ) flshp(1:nx,1:ny,1:nz) = 0.0
          IF ( allocated(flshn) ) flshn(1:nx,1:ny,1:nz) = 0.0
          IF ( allocated(flshi) ) flshi(1:nx,1:ny,1:nz) = 0.0
      
      deallocate ( tem1 )
      deallocate ( tem2 )
      
      RETURN
      END SUBROUTINE V5DOUT
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
  SUBROUTINE makeflshmask(icldtop,nx,ny,nz,ng,st,na,mask,micro)

!   USE GRID_MODULE
!   USE PARAM_MODULE
!   USE MICRO_MODULE
   USE INDEX_MODULE
   
   implicit none
   
   integer :: nx, ny, nz, na, ng
   integer :: icldtop
   integer :: mask(nx,ny)
   real    :: st(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng,na)
   character(LEN = 10) :: micro
   
   integer i,j,k
   
   real :: tem1
   
         icldtop = 0
         DO i=1,nx-1    !2-ibc,nc-2+ibc
          DO j=1,ny-1 ! nr1,nr2
           mask(i,j) = 0
           DO k=1,nz_v5d-1
              IF ( micro(1:5) .eq. 'ICE10' ) THEN
              tem1 =          &
                     1.0e3*(           &
                      st(i,j,k,li) +         &
                      st(i,j,k,lir) +         &
                      st(i,j,k,lip) +         &
                      st(i,j,k,lc) )
              ELSEIF ( li .gt. 1 ) THEN
              tem1 =          &
                     1.0e3*(          &
                      st(i,j,k,li) + st(i,j,k,lc) )
              ENDIF
              IF ( tem1 .gt. 1.0e-2 ) THEN
                mask(i,j) = Max( k, mask(i,j) )
                icldtop = Max( k, icldtop )
              ENDIF
            ENDDO
           ENDDO
          ENDDO


         DO i=1,nx-1    !2-ibc,nc-2+ibc
          DO j=1,ny-1 ! nr1,nr2
            IF ( mask(i,j) .eq. 0 ) mask(i,j) = icldtop/2
          ENDDO
         ENDDO
         
         write(0,*) 'vis5d: icldtop = ',icldtop

      
      RETURN
      END SUBROUTINE makeflshmask

!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
END MODULE VIS5D_MODULE
