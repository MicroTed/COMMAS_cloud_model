        PROGRAM READCOMMASNC

  USE GRID_MODULE
  USE GRIDIO_MODULE
  USE FILE_MODULE
  USE CLINE_MODULE
  USE MICRO_MODULE
  USE PARAM_MODULE
  USE NAMELIST_MODULE
  USE INDEX_MODULE

       implicit none
!-----------------------------------------------------------------------------
! GRID DEFINITION
        
         TYPE(GRID) :: gd

!-----------------------------------------------------------------------------
! COMMAND LINE VARIABLES

         integer              :: status, length, length1,itmp

! INTEGER INDEX VARIABLES FOR THE RUN

         integer :: u 
         integer :: v  
         integer :: w   
         integer :: pi   
         integer :: km   
         integer :: s     
         integer :: dbz   
         integer :: vzf
         integer :: wz    
          
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

         integer, allocatable :: tarray(:)

!-----------------------------------------------------------------------------
! OTHER MISC VARIABLES 

         integer :: ibeg
         integer :: iend
         integer :: ns
         integer :: n
         integer :: ntime
         integer :: m
         integer :: time

         integer i,j,k, it, nt, count
         
         integer :: ixb,ixe,jyb,jye,kzb,kze
         integer :: nx1, ny1, nz1

         character(LEN=4) number1,number2
         character(LEN=8) number

         integer  :: iskip = 1
         
         real,allocatable,dimension(:,:,:) :: u3d,v3d,w3d,tmp3d,pn,dn
         real,allocatable,dimension(:) :: pb,db

! READ IN GRID PARAMETERS

          write(number1, '(a,i3.3)') '.', member
          number = number1


! READ IN THE DATA


           CALL STRING_LIMITS(prefix, ibeg, iend)
           length = iend - ibeg + 1

           IF ( start .lt. 0 ) start = 0
           write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)

           CALL GRID_INFO_NETCDF( prefix(1:length)//trim(number)//'.nc', nt, microphys)
           
           allocate( tarray(nt) )
           
           CALL GRID_INFO_NETCDF( prefix(1:length)//trim(number)//'.nc', nt, microphys, tarray)
           
           microp = microphys
           
           CALL GRID_READ_NETCDF( gd, prefix(1:length)//trim(number)//'.nc', tarray(1), nx_or_nxend=.true. )

           write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION ',prefix(1:length)//trim(number)

           IF( .not. SET_VARIABLE(gd,'TIME',      start)   ) write(6,*) 'COMMAS:  Problem setting TIME'  
          
           time = start


           start = Max(start, tarray(1) ) 

           IF( .not. SET_VARIABLE(gd,'TIME',      Max(0,start))   ) write(6,*) 'COMMAS:  Problem setting TIME'  

               
            gx     = GET_VARIABLE_INDEX(gd, 'XC')
            gy     = GET_VARIABLE_INDEX(gd, 'YC')
            gz     = GET_VARIABLE_INDEX(gd, 'ZC')
!-----------------------------------------------------------------------------
! READ SOME VARIABLES FROM DATA STRUCTURE

         CALL GET_VARIABLE(gd, 'DX',       dx)
         CALL GET_VARIABLE(gd, 'DY',       dy)
         CALL GET_VARIABLE(gd, 'DZ',       dz)
         CALL GET_VARIABLE(gd, 'DT',       dt)
         CALL GET_VARIABLE(gd, 'NXEND',    nx)
         CALL GET_VARIABLE(gd, 'NYEND',    ny)
         CALL GET_VARIABLE(gd, 'NZEND',    nz)
        u      = GET_VARIABLE_INDEX(gd, 'U')
        v      = GET_VARIABLE_INDEX(gd, 'V')
        w      = GET_VARIABLE_INDEX(gd, 'W')
        pi     = GET_VARIABLE_INDEX(gd, 'PI')
        km     = GET_VARIABLE_INDEX(gd, 'KM')
        s      = GET_VARIABLE_INDEX(gd, 'TH') - 1
        dbz    = GET_VARIABLE_INDEX(gd, 'DBZ')
        vzf    = GET_VARIABLE_INDEX(gd, 'VZF')
        wz     = GET_VARIABLE_INDEX(gd, 'WZ')
!        elec   = GET_VARIABLE_INDEX(gd, 'EX')
!        cion   = GET_VARIABLE_INDEX(gd, 'CPIONINIT')
!        muz     = GET_VARIABLE_INDEX(gd, 'MUPOSZC')
        
        ! check for extra diagnostic 3d arrays
!         IF ( elec > wz .and.  elec - wz > 1 ) THEN ! elec is on and extra arrays exist
!            xtra = wz + 1
!            nxtra = elec - wz - 1
!         ELSEIF ( u > wz + 1 ) THEN ! no elec and extra arrays exist
!            xtra = wz + 1
!            nxtra = u - wz - 1
!         ELSE ! no extra arrays
!            xtra = wz
!            nxtra = 0
!         ENDIF

        uinit  = GET_VARIABLE_INDEX(gd, 'UINIT')
        vinit  = GET_VARIABLE_INDEX(gd, 'VINIT')
        winit  = GET_VARIABLE_INDEX(gd, 'WINIT')
        piinit = GET_VARIABLE_INDEX(gd, 'PIINIT')
        kminit = GET_VARIABLE_INDEX(gd, 'KMINIT')
        sinit  = GET_VARIABLE_INDEX(gd, 'THINIT') - 1
        gx     = GET_VARIABLE_INDEX(gd, 'XC')
        gy     = GET_VARIABLE_INDEX(gd, 'YC')
        gz     = GET_VARIABLE_INDEX(gd, 'ZC')
        precip = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')

!
! limits of the region you want
!
        ixb = 1
        ixe = nx-1
        jyb = 1
        jye = ny-1
        kzb = 1
        kze = nz-1
        
! can pass the arrays as, for example u3d(nx1,ny1,nz1)
        nx1 = ixe-ixb+1
        ny1 = jye-jyb+1
        nz1 = kze-kzb+1


        allocate( db(kzb:kze) )
        allocate( pb(kzb:kze) )
        
        allocate( u3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( v3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( w3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( tmp3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( pn(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( dn(ixb:ixe,jyb:jye,kzb:kze) )
        
        DO k = kzb,kze
          db(k) = 1.0e5*gd%var(piinit)%flt1d(k)**2.509/(287.04*gd%var(sinit+lt)%flt1d(k))
          pb(k) = 1.0e5*gd%var(piinit)%flt1d(k)**3.509
        ENDDO

        

!-----------------------------------------------------------------------------
! Main time step loop

          DO it = 1, nt, iskip

!----------------------------------------------------------------------
! Time and location of the grid
        
              time     = tarray(it)

              write(6,*) 'it loop: it,time1 = ',it,time


              IF ( time .gt. stop ) EXIT
              
              IF ( time .ge. start ) THEN

              write(6,*) 'it loop: it,time = ',it,time

              IF ( it > 1 ) CALL GRID_READ_NETCDF( gd, prefix(1:length)//trim(number)//'.nc', tarray(it), 1 )

               CALL GET_VARIABLE(gd, 'TIME',       itmp)
               
               write(6,*) 'time from netcdf is ',itmp

        
        u3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(u)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! u-wind
        v3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(v)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! v-wind
        w3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(w)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! w-wind

        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(dbz)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! dBZ

        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(wz)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! vertical vorticity

        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lt)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! theta
        
        DO k = kzb,kze
        tmp3d(ixb:ixe,jyb:jye,k)=gd%var(s+lt)%flt3d(ixb:ixe,jyb:jye,k) - gd%var(sinit+lt)%flt1d(k) ! pert. theta
        ENDDO
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(km)%flt3d(ixb:ixe,jyb:jye,kzb:kze)   ! TKE (KM, actually)
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lv)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qvapor
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lc)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qcloud
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lr)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qrain
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+li)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qice
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+ls)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qsnow
        
        tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lh)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qgraupel
        
        IF ( lhl > 1 ) tmp3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lhl)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qhail
        
        ! qtot
        tmp3d(ixb:ixe,jyb:jye,kzb:kze) = gd%var(s+lc)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lr)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+li)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+ls)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lh)%flt3d(ixb:ixe,jyb:jye,kzb:kze)
        IF ( lhl > 1 ) THEN
        tmp3d(ixb:ixe,jyb:jye,kzb:kze) = tmp3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lhl)%flt3d(ixb:ixe,jyb:jye,kzb:kze)
        ENDIF

        DO k = kzb,kze
           pn(ixb:ixe,jyb:jye,k) = 1.0e5*(gd%var(piinit)%flt1d(k)+gd%var(pi)%flt3d(ixb:ixe,jyb:jye,k))**3.509 - pb(k)  ! pert. pressure
           dn(ixb:ixe,jyb:jye,k) = 1.0e5*(gd%var(piinit)%flt1d(k)+gd%var(pi)%flt3d(ixb:ixe,jyb:jye,k))**2.509/(rd*gd%var(sinit+lt)%flt1d(k)*(1.0+0.61*gd%var(sinit+lv)%flt1d(k))) ! air density
        ENDDO


! microphysics

!        CALL V5DOUT(gd,                                 &
!                    gd%var(u),  gd%var(uinit) ,         &            ! U,  UINIT
!                    gd%var(v),  gd%var(vinit) ,         &            ! V,  VINIT
!!                    gd%var(w),  gd%var(winit) ,         &            ! W,  WINIT
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
!!                    st,                                 &
!                     gd%xyz3d%flt4d(-ng+1,-ng+1,-ng+1,gd%var(s)%index),  &
!                    ft(1,1),ft(1,2),ft(1,3),ft(1,4),    &
!                    gd%var(dbz), gd%var(vzf), gd%var(wz), gd%var(elec), gd%var(xtra) ,  &
!                    recalcdbz)


               ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
        
       END DO 
          
        STOP
        END
!-------------------------------------------------------------------------------

