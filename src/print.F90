!------------------------------------------------------------------------------
!
!   /////////////////////        BEGIN        \\\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE PRINT   ///////////////////////
!
!------------------------------------------------------------------------------
! Created by Louis Wicker, November, 1990
!
! Modified to work with new database structure, December, 2005
!-------------------------------------------------------------------------------
 SUBROUTINE PRINT(gd)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE

  implicit none

  TYPE(GRID) :: gd

!-------------------------------------------------------------------------------
! Local variables

  TYPE(MAXMIN)             :: mxmn
  TYPE(VARIABLE)           :: var

  integer n, nx, ny, nz, ibeg, iend, i,j,k
  integer time

  character(LEN=name_length), pointer :: list(:)
  character(LEN=name_length)          :: name

  real, pointer :: thinit(:)
  real, pointer :: piinit(:)
  real, pointer :: u3(:,:,:), v3(:,:,:)
  real, pointer :: dxe(:)
  real, pointer :: dye(:)

  real, allocatable, dimension(:) ::  den
!  real, allocatable, dimension(:,:,:) ::  wz1

  integer   :: dbound(6) = -1

!-------------------------------------------------------------------------------
! Get all the variables you need

  CALL GET_VARIABLE (gd, 'NX', nx)
  CALL GET_VARIABLE (gd, 'NY', ny)
  CALL GET_VARIABLE (gd, 'NZ', nz)
  CALL GET_VARIABLE (gd, 'DXE', dxe)
  CALL GET_VARIABLE (gd, 'DYE', dye)
  CALL GET_VARIABLE (gd, 'TIME', time)

  CALL GET_VARIABLE (gd, 'THINIT', thinit,.true.)
  CALL GET_VARIABLE (gd, 'PIINIT', piinit,.true.)

  allocate(den(nz-1))
  den(:) = piinit(1:nz-1)**cvr * p00 / (rd * thinit(1:nz-1))

!--------------------------------------------------------------------------------------
! Open I/O file for model output

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME')

  write(luno,"(1x,71('-'))")
  write(luno,"(1x,'T = ',i6,'.00',4x,a3,29x,a3)") time, 'MAX', 'MIN'
  write(luno,*)

!-------------------------------------------------------------------------------
! First, find all the xyz3d variables

  CALL GET_VARIABLE_LIST(gd, 'xyz3d', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   SELECT CASE( name(ibeg:iend) )

   CASE DEFAULT                                         ! get variable, then mxmn
    
!   CALL VAR_MAXMIN(var, mxmn)

! mixing ratio: scale to g/kg
    IF ( name(ibeg:ibeg) .eq. 'Q' ) THEN
    DO k = 1,nz-1
     var%flt3d(:,:,k) = 1000.*var%flt3d(:,:,k)
    ENDDO
    ENDIF

! space charge: scale to nC/m**3
    IF ( name(ibeg:ibeg+1) .eq. 'SC' ) THEN
    DO k = 1,nz-1
     var%flt3d(:,:,k) = 1.e9*var%flt3d(:,:,k)
    ENDDO
    ENDIF

! particle volume but exclude 'V' wind component
!    IF ( iend .gt. ibeg .and. name(ibeg:ibeg) .eq. 'V' .and. trim(name) /= 'VZF') THEN
    IF ( trim(name) == 'VHW' .or. trim(name) == 'VHL' .or. trim(name) == 'VSW' .or. &
         trim(name) == 'VFW' ) THEN
    DO k = 1,nz-1
     var%flt3d(:,:,k) = 1.e6*var%flt3d(:,:,k)
    ENDDO
    ENDIF

    CALL VAR_MAXMIN(var, mxmn)
   
   CASE ( 'PI' )                                        ! Convert PI -> P(mb)
    DO k = 1,nz-1
     var%flt3d(:,:,k) = var%flt3d(:,:,k)*den(k)*thinit(k)*cp/100. 
    ENDDO
    CALL VAR_MAXMIN(var, mxmn)

   CASE ( 'TH' )                                        ! Theta
    DO k = 1,nz-1
     var%flt3d(:,:,k) = var%flt3d(:,:,k) - thinit(k)
    ENDDO
    CALL VAR_MAXMIN(var, mxmn)

   CASE ( 'WZ' )                                        ! Vorticity
    DO k = 1,nz-1
     var%flt3d(:,:,k) = 1000.*var%flt3d(:,:,k)
    ENDDO
    dbound = (/-1, -1, -1, -1, 1, 5/)
    CALL VAR_MAXMIN(var, mxmn, dbound )

   END SELECT

    write(luno,"(1x,a10,2x,'|',2x,2(g13.6,2x,i4,1x,i4,2x,i3,2x,'|',2x))")    &
          mxmn%name,mxmn%max,mxmn%imax,mxmn%jmax,mxmn%kmax,                &
                    mxmn%min,mxmn%imin,mxmn%jmin,mxmn%kmin

  ENDDO

  IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN
    CALL GET_VARIABLE(gd,'U', u3)
    CALL GET_VARIABLE(gd,'V', v3)
    CALL GET_VARIABLE_COPY(gd, 'WZ', var, copy_variable)                    ! This creates a COPY
    var%flt3d(:,:,1) = 0.0
       DO j = 2,ny-1
        DO i = 2,nx-1
         var%flt3d(i,j,1) = 1000.* (                 &
          (v3(i,j,1)-v3(i-1,j  ,1))*dxe(i)    &
         -(u3(i,j,1)-u3(i,  j-1,1))*dye(j)  ) 
        ENDDO
       ENDDO
     dbound = (/-1, -1, -1, -1, 1, 1/)
     CALL VAR_MAXMIN(var, mxmn, dbound )
     write(luno,"(1x,a10,2x,'|',2x,2(g13.6,2x,i4,1x,i4,2x,i3,2x,'|',2x))")    &
          'Wz-sfc    ',mxmn%max,mxmn%imax,mxmn%jmax,mxmn%kmax,                &
                      mxmn%min,mxmn%imin,mxmn%jmin,mxmn%kmin
    
    ENDIF

  write(luno,"(1x,71('-'))")

!  lun = FILE_CLOSE()
  deallocate( den )
  IF( allocated(var%tmp3d) ) deallocate(var%tmp3d)
  IF( allocated(var%flt2d) ) deallocate(var%flt2d)
  IF( allocated(var%flt1d) ) deallocate(var%flt1d)
!  deallocate( var%flt3d )

 END SUBROUTINE PRINT

!------------------------------------------------------------------------------
!
!   /////////////////////        BEGIN         \\\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE PRTBASE  ///////////////////////
!   
!------------------------------------------------------------------------------
!
! Created by Louis Wicker, November, 1990
! Latest update: 08-24-99
!   
!------------------------------------------------------------------------------

 SUBROUTINE PRTBASE(gd)
      
  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  
  implicit none
      
      
! Passed variables

  TYPE(GRID) :: gd
      
! Local variables

!  integer lun
  integer nz
  real, pointer :: tinit(:), uinit(:), qinit(:), vinit(:), piinit(:), dze(:), zc(:), ze(:), dzc(:)
  real dz
  real pinit(1000)
  real rinit(1000)
  real riinit(1000)
  real rivinit(1000)
  real temp(1000)
  real thetav(1000)
  real dewpt(1000)
  real rhw(1000), rhi(1000)
  
  real bb, bsh, td, ps, bbv
  integer k

  TYPE(ATTRIBUTE), pointer   :: attr
  character(LEN=255)         :: outsoundfile
  integer                    :: ibeg,iend

! Open output file for dump

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'INPUT SOUNDING FROM BASE STATE')

! Find vertical dimension size

  CALL GET_VARIABLE(gd,'NZ',nz)
  CALL GET_VARIABLE(gd,'DZ',dz)

! Get 1D base states
      
  CALL GET_VARIABLE(gd, 'ZC',     zc,     copy_variable)
  CALL GET_VARIABLE(gd, 'DZE',    dze,    copy_variable)
  CALL GET_VARIABLE(gd, 'ZE',     ze,     copy_variable)
  CALL GET_VARIABLE(gd, 'DZC',    dzc,    copy_variable)
  CALL GET_VARIABLE(gd, 'UINIT',  uinit,  copy_variable)
  CALL GET_VARIABLE(gd, 'VINIT',  vinit,  copy_variable)
  CALL GET_VARIABLE(gd, 'THINIT', tinit,  copy_variable)
  CALL GET_VARIABLE(gd, 'QVINIT', qinit,  copy_variable)
  CALL GET_VARIABLE(gd, 'PIINIT', piinit, copy_variable)

    CALL GET_ATTRIBUTE(gd, 'PREFIX_NAME', attr)
    CALL STRING_LIMITS( attr%str, ibeg, iend)
    outsoundfile = attr%str(ibeg:iend)//'.output.sound'

  DO k = 1,nz-1
   temp(k)   = tinit(k)*piinit(k)
   thetav(k) = tinit(k)*(1.0 + 0.61*qinit(k))
   pinit(k)  = 1000.*piinit(k)**3.508
   rinit(k)  = pinit(k)*100./(287.04*tinit(k)*piinit(k))
   rhw(k)    = 3.8*exp(17.27*(tinit(k)*piinit(k)-273.16)/(tinit(k)*piinit(k)-36.))/pinit(k)
   rhw(k)    = qinit(k) / rhw(k)
   rhi(k)    = 3.8*exp(21.87*(tinit(k)*piinit(k)-273.16)/(tinit(k)*piinit(k)-7.66))/pinit(k)
   rhi(k)    = qinit(k) / rhi(k)
       
   ps        = 1.0e3*piinit(k)**3.509
   td        = qinit(k)*ps/(0.622+0.001*qinit(k))                ! vapor pressure
   td        = max(td,0.001)                                     ! avoid problems near zero
   dewpt(k)  = 273.16 + (243.5/( (17.67/alog(td/6.112)) - 1.0))  ! Bolton's approximation
  ENDDO

   thetav(nz) = tinit(nz)*(1.0 + 0.61*qinit(nz))
      
  DO k = 2,nz-1
   bb          = 2.0*g*dze(k)*(tinit(k)-tinit(k-1))/(tinit(k)+tinit(k-1))
   bbv         = 2.0*g*dze(k)*(thetav(k)-thetav(k-1))/(thetav(k)+thetav(k-1))
   bsh         = (uinit(k)-uinit(k-1))**2 + (vinit(k)-vinit(k-1))**2
   bsh         = bsh * dze(k)**2
   riinit(k-1) = min(bb/(bsh+1.0e-10),50.)
   rivinit(k-1) = min(bbv/(bsh+1.0e-10),50.)
  ENDDO

  riinit(nz-1) = riinit(nz-2)

  write(luno,*)
  write(luno,"(1x,4x,'HEIGHT',5x,'PRESS',5x,'DEN',5x,'THETA',5x,'TEMP',4x,'THETAV',3x,'DEWPT',6x, &
 &  'QV',10x,'U',10x,'V',7x,'RHw',6x,'RHi',6x,'RI',6x,'RIV')")

  DO k = nz-1,1,-1
   write(luno,"(1x,f11.2,3x,f7.2,3x,f7.4,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f7.4,3x,2(f8.3,3x),4(f6.2,3x))")  &
            zc(k),pinit(k),rinit(k),tinit(k),temp(k),thetav(k),dewpt(k),                           &
            qinit(k)*1000.,uinit(k),vinit(k), rhw(k)*100., rhi(k)*100.,riinit(k),rivinit(k)
  ENDDO

  write(luno,*)
  write(luno,*)
  write(luno,*) 'z1d: level, zc, ze, dzc, dze'
  DO k = 1,nz
    write(luno,'(i3,2x,5(1pe13.5,2x))') k, zc(k),ze(k),dzc(k),dze(k), dz
  ENDDO

   open(21,file=outsoundfile,form='formatted', status='unknown', position='append')
!   rewind(21)
!   write(21,*) psfc, tsfc, qsfc*1000.
   DO k = 1,nz-1
    write(21,'(1x,2(f8.2,2x),f9.4,2x, 2(f10.3,2x),f9.1)') zc(k), tinit(k), qinit(k)*1000., uinit(k), vinit(k), 100.*pinit(k)
   ENDDO

   close(21)

!  luno = FILE_CLOSE('END INPUT SOUNDING')

 RETURN
 END SUBROUTINE PRTBASE
 
 
 
 
!------------------------------------------------------------------------------
!
!   /////////////////////         BEGIN          \\\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE PRTINFO    ///////////////////////
!
!------------------------------------------------------------------------------
! Created by Louis Wicker, summer 2004
! Dump out various information...
!-------------------------------------------------------------------------------
 SUBROUTINE PRTINFO(gd)
  
  USE GRID_MODULE
  USE FILE_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
    
  implicit none
  
  TYPE(GRID) :: gd
  

  integer pos(10)
  integer i, j
  
  integer nx, ny
  real, pointer :: xc(:), dx(:)
  real, pointer :: yc(:), dy(:)
  real flt, flt2
  real scale
 
!-------------------------------------------------------------------------------
! Open output file for dump

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'COMMAS RUN PARAMETERS')

!-------------------------------------------------------------------------------
! Determine grid size

  CALL GET_VARIABLE(gd, 'NX', nx)
  CALL GET_VARIABLE(gd, 'NY', ny)
  CALL GET_VARIABLE(gd, 'XC', xc)
  CALL GET_VARIABLE(gd, 'YC', yc)
  CALL GET_VARIABLE(gd, 'DXC', dx)
  CALL GET_VARIABLE(gd, 'DYC', dy)


  write(luno,'(1x,80("-"))')
  write(luno,'(1x,"COMMAS SOLVER PARAMETERS")')
  write(luno,*)
  write(luno,'(1x,"SOLVER:  RKSCHEME                   :  ",i10)') RKSCHEME
  write(luno,*)
  write(luno,'(1x,"PRANDTL NUMBER (1/Pr)               :  ",f10.5)') 1.0/Pr
  write(luno,'(1x,"DIVERGENCE DAMPING COEFFICIENT      :  ",f10.5)') kdiv0
  write(luno,'(1x,"PARAMETER FOR VERTICAL IMPLICIT     :  ",f10.5)') alpha
  write(luno,'(1x,"SOUND WAVE SPEED                    :  ",f10.5)') cspd
  write(luno,'(1x,"OUTFLOW BOUNDARY GRAVITY WAVE SPEED :  ",f10.5)') dxt
  write(luno,'(1x,"CORIOLIS PARAMETER                  :  ",f10.5)') coriol
  write(luno,'(1x,"RAYLEIGH DAMPING HEIGHT             :  ",f10.1)') rayd_hgt
  write(luno,'(1x,"RAYLEIGH DAMPING COEFFICIENT        :  ",f10.5)') rayd_mag
  write(luno,'(1x,"HORIZONTAL SPONGE LAYER             :  ",i10)'  ) hsponge
  write(luno,'(1x,"HORIZONTAL SPONGE COEFFICIENT       :  ",f10.5)') hrayd_mag
  write(luno,'(1x,"NON-DIMENSIONAL MIXING COEFF Cm/Ce  :  ",2(1x,f10.5))') Cm, Ce
  IF( vert_adv_scheme .eq. 3) write(luno,'(1x,"VERTICAL ADVECTION IS COMPUTED AT 3rd/4th-ORDER")')
  IF( vert_adv_scheme .ne. 3) write(luno,'(1x,"VERTICAL ADVECTION IS COMPUTED AT 5th/6th-ORDER")')
  write(luno,*)

                         write(luno,'(1x,"AHIGHK                              :  ",f10.5)') ahighk
  IF( mix_type .eq. -1 ) write(luno,'(1x,"CONSTANT MIXING COEFFICIENT         :  ",f10.5)') ahighk
  IF( mix_type .eq.  0 ) write(luno,'(1x,"SMAGORINSKY MIXING                 ")')
  IF( mix_type .eq.  1 ) write(luno,'(1x,"TKE MIXING                 ")')
  IF( len_type .eq.  0 ) write(luno,'(1x,"MIXING LENGTH IS VOLUME**1/3     ")')
  IF( len_type .eq.  1 ) write(luno,'(1x,"MIXING LENGTH IS ~ DZ     ")')
  IF( len_type .eq.  2 ) write(luno,'(1x,"MIXING LENGTH IS BOUNDARY LAYER FORMULATION     ")')
  IF( len_type .eq.  2 ) write(luno,'(1x,"MAXIMUM BL MIXING LENGTH            :  ",f10.5)') Lmax
  IF( len_type .eq.  3 ) write(luno,'(1x,"USING USER SPECIFIED MIXING LENGTH  :  ",f10.5)') Lmax
  IF( rmbasekm .eq.  1 ) write(luno,'(1x,"BASE STATE MIXING IS REMOVED")')
                         write(luno,*)
  IF( bcz .eq.  1 )      write(luno,'(1x,"LOWER BOUNDARY CONDITION IS SEMISLIP")')
  IF( bcz .eq.  1 )      write(luno,'(1x,"DRAG COEFFICIENT                    :  ",f10.5)') drag
  IF( bcz .eq.  1 )      write(luno,*)
  IF( MONOTONIC )        write(luno,'(1x,"WENO      ADVECTION ON MOMENTUM/THETA")')
  IF( MONOTONIC )        write(luno,'(1x,"MONOTONIC ADVECTION ON SCALARS/TKE")')

  write(luno,*)

  write(luno,*) '-----------------------------------------------------------------------'
  write(luno,*) 'LFO MICROPHYSICAL CONSTANTS'

  CALL GET_VARIABLE(gd, 'CNOR', flt)
  CALL GET_VARIABLE(gd, 'RHO_QR', flt2)
    
  write(luno,"(1x,a,6(g12.5,2x))") 'CNOR/RHO_QR:        ',flt, flt2

  CALL GET_VARIABLE(gd, 'CNOS', flt)
  CALL GET_VARIABLE(gd, 'RHO_QS', flt2)
    
  write(luno,"(1x,a,6(g12.5,2x))") 'CNOS/RHO_QS:        ',flt, flt2

  CALL GET_VARIABLE(gd, 'CNOH', flt)
  CALL GET_VARIABLE(gd, 'RHO_QH', flt2)
    
  write(luno,"(1x,a,6(g12.5,2x))") 'CNOH/RHO_QH:        ', flt, flt2

  IF( autoconversion .eq. 0 ) THEN
   write(luno,*) 'Simple qc -> qr autoconversion, QC threshold = ',qcmincwrn
  ELSE
   write(luno,*) 'autoconversion = ', autoconversion
   write(luno,*) 'cwdiap, cwdisp = ',cwdiap, cwdisp
   write(luno,*) 'Berry (1968) AC / CCW = ',ccn
  ENDIF
   
  write(luno,*) '-----------------------------------------------------------------------'

  write(luno,*)
  write(luno,'(1x,80("-"))')

! PRINT OUT SOME GRID DOMAIN INFORMATION
  
  write(luno,'(1x,"GRID PARAMETERS")')
  write(luno,*)
  
  CALL GET_VARIABLE(gd, 'LAT', flt)
  write(luno,'(1x,"STARTING LATITUDE  OF GRID: ",1x,f8.3)') flt
  CALL GET_VARIABLE(gd, 'LON', flt)
  write(luno,'(1x,"STARTING LONGITUDE OF GRID: ",1x,f8.3)') flt
  CALL GET_VARIABLE(gd, 'HGT', flt)
  write(luno,'(1x,"STARTING ALTITUDE  OF GRID: ",1x,f8.3)') flt
  
  pos(:) = 1
  pos(6) = nx/2 + 1

  DO i = 1,nx-1
   IF( xc(i) .le. 10000. .and. xc(i+1) .gt. 10000. ) pos(2) = i
   IF( xc(i) .le. 20000. .and. xc(i+1) .gt. 20000. ) pos(3) = i
   IF( xc(i) .le. 30000. .and. xc(i+1) .gt. 30000. ) pos(4) = i
   IF( xc(i) .le. 40000. .and. xc(i+1) .gt. 40000. ) pos(5) = i
  ENDDO
  
  DO i = 1,nx-1
   IF( xc(i) .le. 100. .and. xc(i+1) .gt. 100. ) pos(2) = i
   IF( xc(i) .le. 200. .and. xc(i+1) .gt. 200. ) pos(3) = i
   IF( xc(i) .le. 300. .and. xc(i+1) .gt. 300. ) pos(4) = i
   IF( xc(i) .le. 400. .and. xc(i+1) .gt. 400. ) pos(5) = i
  ENDDO

  scale = 1.0e-3

  write(luno,'(1x,80("-"))')
  write(luno,*)
  write(luno,*) 'X-GRID DISTRIBUTION (ASSUME SYMMETRY ABOUT CENTER)'

  write(luno,"(1x,'GRID LOCATION(KM)     = ',6(f6.0,2x))")  &
                   scale*xc(pos(1)),						&
                   scale*xc(pos(2)),						&
                   scale*xc(pos(3)),						&
                   scale*xc(pos(4)),						&
                   scale*xc(pos(5)),						&
                   scale*xc(pos(6))

  write(luno,"(1x,'GRID SPACING (METERS) = ',6(f6.0,2x))")  &
                   1.0/dx(pos(1)),							&
                   1.0/dx(pos(2)),							&
                   1.0/dx(pos(3)),							&
                   1.0/dx(pos(4)),							&
                   1.0/dx(pos(5)),							&
                   1.0/dx(pos(6))

  write(luno,*)

  pos(:) = 1
  pos(6) = ny/2 + 1

  DO j = 1,ny-1
   IF( yc(j) .le. 10000. .and. yc(j+1) .gt. 10000. ) pos(2) = j
   IF( yc(j) .le. 20000. .and. yc(j+1) .gt. 20000. ) pos(3) = j
   IF( yc(j) .le. 30000. .and. yc(j+1) .gt. 30000. ) pos(4) = j
   IF( yc(j) .le. 40000. .and. yc(j+1) .gt. 40000. ) pos(5) = j
  ENDDO
 
  DO j = 1,ny-1
   IF( yc(j) .le. 100. .and. yc(j+1) .gt. 100. ) pos(2) = j
   IF( yc(j) .le. 200. .and. yc(j+1) .gt. 200. ) pos(3) = j
   IF( yc(j) .le. 300. .and. yc(j+1) .gt. 300. ) pos(4) = j
   IF( yc(j) .le. 400. .and. yc(j+1) .gt. 400. ) pos(5) = j
  ENDDO

  write(luno,'(1x,80("-"))')
  write(luno,*)
  write(luno,*) 'Y-GRID DISTRIBUTION (ASSUME SYMMETRY ABOUT CENTER)'

  write(luno,"(1x,'GRID LOCATION(KM)     = ',6(f6.0,2x))")  &
                   scale*yc(pos(1)),						&
                   scale*yc(pos(2)),						&
                   scale*yc(pos(3)),						&
                   scale*yc(pos(4)),						&
                   scale*yc(pos(5)),						&
                   scale*yc(pos(6))

  write(luno,"(1x,'GRID SPACING (METERS) = ',6(f6.0,2x))")  &
                   1.0/dy(pos(1)),							&
                   1.0/dy(pos(2)),							&
                   1.0/dy(pos(3)),							&
                   1.0/dy(pos(4)),							&
                   1.0/dy(pos(5)),							&
                   1.0/dy(pos(6))
  write(luno,*)
  write(luno,'(1x,80("-"))')
  
!-------------------------------------------------------------------------------
! Open output file for dump

!  luno = FILE_CLOSE()

END SUBROUTINE PRTINFO

!------------------------------------------------------------------------------
!
!   /////////////////////        BEGIN            \\\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE PRINT_ALL   ///////////////////////
!
!------------------------------------------------------------------------------
! Created by Louis Wicker, March, 2003
!
!-------------------------------------------------------------------------------
 SUBROUTINE PRINT_ALL( gd )

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE

  implicit none

  TYPE(GRID) :: gd

!-------------------------------------------------------------------------------
! Local variables

  TYPE(MAXMIN)             :: mxmn
  TYPE(VARIABLE)           :: var

  integer n, nx, ny, nz, ibeg, iend, k
  integer time

  character(LEN=name_length), pointer :: list(:)
  character(LEN=name_length)          :: name

!-------------------------------------------------------------------------------
! Get all the variables you need

  CALL GET_VARIABLE (gd, 'TIME', time)

!--------------------------------------------------------------------------------------
! Open I/O file for model output

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME', 'PRINT_ALL DIAGNOSTICS')

  write(luno,"(1x,71('-'))")
  write(luno,"(1x,'T = ',i6,'.00')") time
  write(luno,*)

!-------------------------------------------------------------------------------
! Find and list all the iconst variables

  CALL GET_VARIABLE_LIST(gd, 'icnst', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   write(luno,*) 'PRINT_ALL:  ICONST = ',name(ibeg:iend), var%int

  ENDDO

!-------------------------------------------------------------------------------
! Find and list all the rconst variables

  CALL GET_VARIABLE_LIST(gd, 'rcnst', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   write(luno,*) 'PRINT_ALL:  RCONST = ',name(ibeg:iend), var%flt

  ENDDO

!-------------------------------------------------------------------------------
! Find and list all the x1d variables

  CALL GET_VARIABLE_LIST(gd, 'x1d', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   nz = size(var%flt1d) - 2*var%ng - 1 + var%istag

   write(luno,*) 'PRINT_ALL: X1D = ',name(ibeg:iend), var%flt1d(1), var%flt1d(nz/2), var%flt1d(nz), &
                 var%istag, nz

  ENDDO

!-------------------------------------------------------------------------------
! Find and list all the y1d variables

  CALL GET_VARIABLE_LIST(gd, 'y1d', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   nz = size(var%flt1d) - 2*var%ng - 1 + var%jstag

   write(luno,*) 'PRINT_ALL: Y1D = ',name(ibeg:iend), var%flt1d(1), var%flt1d(nz/2), var%flt1d(nz), &
                 var%jstag, nz

  ENDDO

!-------------------------------------------------------------------------------
! Find and list all the z1d variables

  CALL GET_VARIABLE_LIST(gd, 'z1d', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   nz = size(var%flt1d) - 2*var%ng - 1 + var%kstag

   write(luno,*) 'PRINT_ALL: Z1D = ',name(ibeg:iend), var%flt1d(1), var%flt1d(nz/2), var%flt1d(nz), &
                 var%kstag, nz

  ENDDO

!-------------------------------------------------------------------------------
! Find and list all the 3D variables

  CALL GET_VARIABLE_LIST(gd, 'xyz3d', list)

  DO n = 1,size(list)

   name = list(n)

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)                    ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   nx = size(var%flt3d,dim=1) - 2*var%ng - 1 + var%istag
   ny = size(var%flt3d,dim=2) - 2*var%ng - 1 + var%jstag
   nz = size(var%flt3d,dim=3) - 2*var%ng - 1 + var%kstag

   write(luno,*) 'PRINT_ALL: XYZ3D = ',name(ibeg:iend), var%flt3d(1,1,1),var%flt3d(nx/2+1,ny/2+1,nz/2),var%flt3d(nx,ny,nz), &
                 var%istag, var%jstag, var%kstag

  ENDDO

  IF( allocated(var%tmp3d) ) deallocate(var%tmp3d)
  IF( allocated(var%flt2d) ) deallocate(var%flt2d)
  IF( allocated(var%flt1d) ) deallocate(var%flt1d)

  write(luno,"(1x,71('-'))")

 END SUBROUTINE PRINT_ALL

#ifdef MPI
!------------------------------------------------------------------------------
!
!   /////////////////////        BEGIN        \\\\\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\  SUBROUTINE PRINT   ///////////////////////
!
!------------------------------------------------------------------------------
! Created by Louis Wicker, November, 1990
!
! Modified to work with new database structure, December, 2005
! Modified to print global max/mins in MPI by APS, January, 2008
!-------------------------------------------------------------------------------
 SUBROUTINE PRINT_MPI(gd)

  USE GRID_MODULE
  USE FILE_MODULE
  USE PARAM_MODULE
  USE COMMASMPI_MODULE

  implicit none

  INCLUDE "mpif.h"

  TYPE(GRID) :: gd

!-------------------------------------------------------------------------------
! Local variables

  TYPE(MAXMIN)             :: mxmn
  TYPE(VARIABLE)           :: var

  integer n, nx, ny, nz, ibeg, iend, i,j,k
  integer time

  character(LEN=name_length), pointer :: list(:)
  character(LEN=name_length)          :: name

  real, pointer :: thinit(:)
  real, pointer :: piinit(:)
  real, pointer :: u3(:,:,:), v3(:,:,:)
  real, pointer :: dxe(:)
  real, pointer :: dye(:)

  real, allocatable, dimension(:) ::  den
!  real, allocatable, dimension(:,:,:) ::  wz1

  integer :: dbound(6) = -1

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
   integer :: itype,imxmn
   integer :: irank1,irank2
   
   integer :: listsize

   integer :: mxmn_size,max_rank,min_rank
   integer :: loc_send_size,loc_recv_size

   integer,allocatable :: mpiloc(:,:,:), mpiloc0(:,:,:,:)
   integer,allocatable :: loc_send_buffer(:,:,:),loc_recv_buffer(:,:,:)

   real,allocatable :: max_send_buffer(:,:),max_recv_buffer(:,:)
   real,allocatable :: min_send_buffer(:,:),min_recv_buffer(:,:)

   real,allocatable :: mpimxmn(:,:)

   character(LEN=name_length), allocatable :: listnames(:)

   logical :: debug_mpi = .false.
   
!-------------------------------------------------------------------------------
! Get all the variables you need

  CALL GET_VARIABLE (gd, 'NX', nx)
  CALL GET_VARIABLE (gd, 'NY', ny)
  CALL GET_VARIABLE (gd, 'NZ', nz)
  CALL GET_VARIABLE (gd, 'DXE', dxe)
  CALL GET_VARIABLE (gd, 'DYE', dye)
  CALL GET_VARIABLE (gd, 'TIME', time)

  CALL GET_VARIABLE (gd, 'THINIT', thinit,.true.)
  CALL GET_VARIABLE (gd, 'PIINIT', piinit,.true.)

  kze = nz
  IF ( myprock == nprock ) kze = nz-1

  allocate(den(nz))
  den(1:kze) = piinit(1:kze)**cvr * p00 / (rd * thinit(1:kze))

!--------------------------------------------------------------------------------------
! Open I/O file for model output

  luno = FILE_OPEN(gd, 'OUTPUT_FILE_NAME')

  write(luno,"(1x,71('-'))")
  write(luno,"(1x,'T = ',i6,'.00',4x,a3,28x,a3)") time, 'MAX', 'MIN'
  write(luno,*)

!-------------------------------------------------------------------------------
! First, find all the xyz3d variables

  CALL GET_VARIABLE_LIST(gd, 'xyz3d', list)
  listsize = size(list) 
  IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN
    listsize = listsize + 1 ! add 1 for Wz-sfc
  ENDIF
  
  allocate(listnames(listsize))
  allocate(mpimxmn(listsize,2))
  allocate(max_send_buffer(2,listsize))
  allocate(max_recv_buffer(2,listsize))
  allocate(min_send_buffer(2,listsize))
  allocate(min_recv_buffer(2,listsize))

  allocate(mpiloc(listsize,3,2))
  allocate(mpiloc0(listsize,3,2,0:number_of_processes-1))
  allocate(loc_send_buffer(listsize,3,2))
  allocate(loc_recv_buffer(listsize*number_of_processes,3,2))

  
  DO n = 1,size(list)

   name = list(n)
   listnames(n) = name

   CALL GET_VARIABLE_COPY(gd, name, var, copy_variable)        ! This creates a COPY

   CALL STRING_LIMITS( name, ibeg, iend )

   SELECT CASE( name(ibeg:iend) )

   CASE DEFAULT                                         ! get variable, then mxmn
    
!   CALL VAR_MAXMIN(var, mxmn)

! mixing ratio: scale to g/kg
    IF ( name(ibeg:ibeg) .eq. 'Q' ) THEN
    DO k = 1,kze
     var%flt3d(:,:,k) = 1000.*var%flt3d(:,:,k)
    ENDDO
    ENDIF

! space charge: scale to nC/m**3
    IF ( name(ibeg:ibeg+1) .eq. 'SC' ) THEN
    DO k = 1,kze
     var%flt3d(:,:,k) = 1.e9*var%flt3d(:,:,k)
    ENDDO
    ENDIF

! particle volume but exclude 'V' wind component
!    IF ( iend .gt. ibeg .and. name(ibeg:ibeg) .eq. 'V' .and. trim(name) /= 'VZF') THEN
    IF ( trim(name) == 'VHW' .or. trim(name) == 'VHL' .or. trim(name) == 'VSW' ) THEN
    DO k = 1,kze
     var%flt3d(:,:,k) = 1.e6*var%flt3d(:,:,k)
    ENDDO
    ENDIF

    CALL VAR_MAXMIN(var, mxmn)

     mpimxmn(n,1) = mxmn%max
     mpiloc(n,1,1) = mxmn%imax
     mpiloc(n,2,1) = mxmn%jmax
     mpiloc(n,3,1) = mxmn%kmax

     mpimxmn(n,2) = mxmn%min
     mpiloc(n,1,2) = mxmn%imin
     mpiloc(n,2,2) = mxmn%jmin
     mpiloc(n,3,2) = mxmn%kmin

   CASE ( 'PI' )                                        ! Convert PI -> P(mb)
    DO k = 1,kze
     var%flt3d(:,:,k) = var%flt3d(:,:,k)*den(k)*thinit(k)*cp/100. 
    ENDDO
    CALL VAR_MAXMIN(var, mxmn)

     mpimxmn(n,1) = mxmn%max
     mpiloc(n,1,1) = mxmn%imax
     mpiloc(n,2,1) = mxmn%jmax
     mpiloc(n,3,1) = mxmn%kmax

     mpimxmn(n,2) = mxmn%min
     mpiloc(n,1,2) = mxmn%imin
     mpiloc(n,2,2) = mxmn%jmin
     mpiloc(n,3,2) = mxmn%kmin

   CASE ( 'TH' )                                        ! Theta
    DO k = 1,kze
     var%flt3d(:,:,k) = var%flt3d(:,:,k) - thinit(k)
    ENDDO
    CALL VAR_MAXMIN(var, mxmn)

     mpimxmn(n,1) = mxmn%max
     mpiloc(n,1,1) = mxmn%imax
     mpiloc(n,2,1) = mxmn%jmax
     mpiloc(n,3,1) = mxmn%kmax

     mpimxmn(n,2) = mxmn%min
     mpiloc(n,1,2) = mxmn%imin
     mpiloc(n,2,2) = mxmn%jmin
     mpiloc(n,3,2) = mxmn%kmin

   CASE ( 'WZ' )                                        ! Vorticity
    IF ( myprock == 1 ) THEN
    DO k = 1,kze
     var%flt3d(:,:,k) = 1000.*var%flt3d(:,:,k)
    ENDDO
    dbound = (/-1, -1, -1, -1, 1, 5/)
    CALL VAR_MAXMIN(var, mxmn, dbound )

     mpimxmn(n,1) = mxmn%max
     mpiloc(n,1,1) = mxmn%imax
     mpiloc(n,2,1) = mxmn%jmax
     mpiloc(n,3,1) = mxmn%kmax

     mpimxmn(n,2) = mxmn%min
     mpiloc(n,1,2) = mxmn%imin
     mpiloc(n,2,2) = mxmn%jmin
     mpiloc(n,3,2) = mxmn%kmin
     
     ELSE

     mpimxmn(n,1)  = 0.0 ! mxmn%max
     mpiloc(n,1,1) = 1 ! mxmn%imax
     mpiloc(n,2,1) = 1 ! mxmn%jmax
     mpiloc(n,3,1) = 1 ! mxmn%kmax

     mpimxmn(n,2)  = 0.0 ! mxmn%min
     mpiloc(n,1,2) = 1 ! mxmn%imin
     mpiloc(n,2,2) = 1 ! mxmn%jmin
     mpiloc(n,3,2) = 1 ! mxmn%kmin
     
     ENDIF

   END SELECT

  ENDDO

  IF ( nx .gt. 2 .and. ny .gt. 2 ) THEN  !note this is only for my_rank=0!
    n = listsize
    listnames(n) = 'Wz-sfc'
    CALL GET_VARIABLE(gd,'U', u3)
    CALL GET_VARIABLE(gd,'V', v3)
    CALL GET_VARIABLE_COPY(gd, 'WZ', var, copy_variable)                    ! This creates a COPY
    var%flt3d(:,:,1) = 0.0

     IF ( myprock == 1 ) THEN
       ixb = 1
       ixe = itile
       if (ixbeg .eq. nxbeg) ixb = 2
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       jyb = 1
       jye = jtile
       if (jybeg .eq. nybeg) jyb = 2
       if (jyend .eq. nyend) jye = jyend-jybeg

       DO k = 1,1
       DO j = jyb,jye
        DO i = ixb,ixe
         var%flt3d(i,j,k) = 1000.* (          &
          (v3(i,j,k)-v3(i-1,j  ,k))*dxe(i)    &
         -(u3(i,j,k)-u3(i,  j-1,k))*dye(j)  ) 
        ENDDO
       ENDDO
       ENDDO
     ENDIF
     
     dbound = (/-1, -1, -1, -1, 1, 1/)
     CALL VAR_MAXMIN(var, mxmn, dbound )


     mpimxmn(n,1) = mxmn%max
     mpiloc(n,1,1) = mxmn%imax
     mpiloc(n,2,1) = mxmn%jmax
     mpiloc(n,3,1) = mxmn%kmax
     
     mpimxmn(n,2) = mxmn%min
     mpiloc(n,1,2) = mxmn%imin
     mpiloc(n,2,2) = mxmn%jmin
     mpiloc(n,3,2) = mxmn%kmin
    
    ENDIF

! Find global maxmin values and which proc they are on

    DO n = 1,listsize

     max_send_buffer(1,n) = mpimxmn(n,1)
     max_send_buffer(2,n) = real(my_rank)

     min_send_buffer(1,n) = mpimxmn(n,2)
     min_send_buffer(2,n) = real(my_rank)

    ENDDO

  mxmn_size = listsize

  CALL MPI_Allreduce(max_send_buffer, max_recv_buffer, mxmn_size, MPI_2REAL,  &
                       MPI_MAXLOC, my_comm, mpi_error_code)

  CALL MPI_Allreduce(min_send_buffer, min_recv_buffer, mxmn_size, MPI_2REAL,  &
                       MPI_MINLOC, my_comm, mpi_error_code)

  IF (mpi_error_code /= MPI_SUCCESS) THEN

   IF (debug_mpi) THEN
    WRITE (0,*) my_rank, "PRINT_MPI: about to CALL MPI_Abort after Allreduce attempt"
   END IF !! (debug_mpi)

   CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

   IF (debug_mpi) THEN
    WRITE (0,*) my_rank, "PRINT_MPI: done calling MPI_Abort after Allreduce attempt"
   END IF !! (debug_mpi)

  END IF !! (mpi_error_code /= MPI_SUCCESS)

  DO n = 1,listsize

   mpimxmn(n,1) = max_recv_buffer(1,n)
   mpimxmn(n,2) = min_recv_buffer(1,n)

  ENDDO

!Now find out the coords of the maxmin

     DO imxmn = 1,2
      DO n = 1,listsize
       loc_send_buffer(n,1,imxmn) = ixbeg-1+mpiloc(n,1,imxmn)
       loc_send_buffer(n,2,imxmn) = jybeg-1+mpiloc(n,2,imxmn)
       loc_send_buffer(n,3,imxmn) = kzbeg-1+mpiloc(n,3,imxmn)
      ENDDO
     ENDDO

    loc_send_size = listsize*3*2
    loc_recv_size = listsize*3*2

    CALL MPI_Gather(loc_send_buffer,loc_send_size,MPI_INTEGER,   &
                    loc_recv_buffer,loc_recv_size,MPI_INTEGER,   &
                    0,my_comm,mpi_error_code)

    IF (mpi_error_code /= MPI_SUCCESS) THEN

     IF (debug_mpi) THEN
      WRITE (0,*) my_rank, "PRINT_MPI: about to CALL MPI_Abort after Gather attempt"
     END IF !! (debug_mpi)

     CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

     IF (debug_mpi) THEN
      WRITE (0,*) my_rank, "PRINT_MPI: done calling MPI_Abort after Gather attempt"
     END IF !! (debug_mpi)

    END IF !! (mpi_error_code /= MPI_SUCCESS)
!
  IF (my_rank == 0) THEN

    mpiloc0 = reshape(loc_recv_buffer,(/listsize,3,2,number_of_processes/))

      DO itype = 1,3
       DO n = 1,listsize
       
        irank1 = int(max_recv_buffer(2,n))
        irank2 = int(min_recv_buffer(2,n))

        mpiloc(n,itype,1) = mpiloc0(n,itype,1,irank1)
        mpiloc(n,itype,2) = mpiloc0(n,itype,2,irank2)

      ENDDO
     ENDDO

  DO n = 1,listsize
    write(luno,"(1x,a10,2x,'|',2x,2(g13.6,2x,i4,1x,i4,2x,i4,2x,'|',2x))")       &
          listnames(n),mpimxmn(n,1),mpiloc(n,1,1),mpiloc(n,2,1),mpiloc(n,3,1),  &
                       mpimxmn(n,2),mpiloc(n,1,2),mpiloc(n,2,2),mpiloc(n,3,2)
  ENDDO


  ENDIF !! my_rank==0

  write(luno,"(1x,71('-'))")

!  lun = FILE_CLOSE()
  deallocate( den )
  IF( allocated(var%tmp3d) ) deallocate(var%tmp3d)
  IF( allocated(var%flt2d) ) deallocate(var%flt2d)
  IF( allocated(var%flt1d) ) deallocate(var%flt1d)
!  deallocate( var%flt3d )

  deallocate(max_send_buffer)
  deallocate(max_recv_buffer)
  deallocate(min_send_buffer)
  deallocate(min_recv_buffer)
  deallocate(loc_send_buffer)
  deallocate(loc_recv_buffer)
  deallocate(mpimxmn)
  deallocate(mpiloc)
  deallocate(mpiloc0)

!mpidebug
  if (debug_mpi) write(0,*) my_rank, "PRINT_MPI: Exiting subroutine..."

 END SUBROUTINE PRINT_MPI
#endif 
