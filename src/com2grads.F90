!
! Reads a commas netcdf file (currently hardwired variable 'filename') and
! converts all times in the file to a grads file 'commas.dat'
!


        PROGRAM COM2GRADS

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

         integer   :: nvars
         parameter (nvars=27)
! INTEGER INDEX VARIABLES FOR THE RUN

         integer :: u,irec10,irec20,irec30,irec40
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
         integer :: rain_rat,hail_acc,hail_rat,rain_acc 
         real    :: tlcl,qvuse,cpmix,cvmix,rgmix,prrcp

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

         integer i,j,k, it, nt, count,itstart,itend
         integer :: ixb,ixe,jyb,jye,kzb,kze
         integer :: nx1, ny1, nz1

         character(LEN=4) number1
         character(LEN=8) number
         character filename*30 /'may29_3mom_traj_500m.000.nc'/
         integer  :: iskip = 3
         
         real,allocatable,dimension(:,:,:) :: precipr,precipr2,preciph,preciph2 
         real,allocatable,dimension(:,:,:) :: temp,thetae3,thetae,thetae2,RH,press2,u3d,v3d,w3d,dbz3d,pn,dn,tke
         real,allocatable,dimension(:,:,:) :: qv,qr,qi,qs,qg,qh,qc,qtot,vort,theta,thetap,thetae4
         real,allocatable,dimension(:) :: pb,db
         real,allocatable,dimension(:,:,:,:) :: gga0
         real,allocatable,dimension(:,:) :: t1,t2,t3,t4

         real :: ntt,thetaee

! READ IN GRID PARAMETERS

          write(number1, '(a,i3.3)') '.', member
          number = number1


      irec10=1
      irec20=0
      irec30=0
      irec40=0

! READ IN THE DATA


           IF ( start .lt. 0 ) start = 0
           write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION '

           CALL GRID_INFO_NETCDF( trim(filename), nt, microphys)
           
           allocate( tarray(nt) )
           
           CALL GRID_INFO_NETCDF( trim(filename), nt, microphys, tarray)
           
           microp = microphys
           
           CALL GRID_READ_NETCDF( gd, trim(filename), tarray(1), nx_or_nxend=.true. )

           write(6,*) 'COMMAS:  START = ',start,' READING IN GRID FOR SIMULATION '

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
        rain_rat = GET_VARIABLE_INDEX(gd, 'RAIN_RAT')
        rain_acc = GET_VARIABLE_INDEX(gd, 'RAIN_ACC')
        hail_rat = GET_VARIABLE_INDEX(gd, 'HAIL_RAT')
        hail_acc = GET_VARIABLE_INDEX(gd, 'HAIL_ACC')

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

        ntt=(nt/iskip) + 1
!       ntt=15 ! for commas2.ctl 3600 sec / 240 sec = 15

       print*, 'nx=',ixe,'ny=',jye,'nz=',kze
       print*, 'nt=',nt
       print*, 'tarray=',tarray(1:nt)

       open(unit=33,file='commas.ctl',form='formatted')

       write(33,*)   'DSET commas.dat'
       write(33,*)   'TITLE COMMAS'
       write(33,*) 'UNDEF -9.99E33'
       write(33,*) 'XDEF ',ixe,'  LEVELS'
       do i=1,ixe
       write(33,*) '',i*dx*0.001
       enddo
       write(33,*) 'YDEF ',jye,'  LEVELS'
      do j=1,jye
      write(33,*) '',j*dy*0.001
      enddo
      write(33,*) 'ZDEF ',kze,'   LEVELS'
      do k=1,kze
      write(33,*) '',k*dz*0.001
      enddo
       write(33,*)  'TDEF ',int(ntt),'LINEAR 00:00Z02jun1995 1hr '
!      do k=1,ntt
!      write(33,*) '',k*(tarray(2)-tarray(1) )*iskip
!      enddo
       write(33,*)   'VARS ',nvars
       write(33,*)'Z1 ',kze,'   201    u'
       write(33,*)'Z2 ',kze,'   201    v'
       write(33,*)'Z3 ',kze,'   201    w'
       write(33,*)'Z4 ',kze,'   201    theta'
       write(33,*)'Z5 ',kze,'   201    tke'
       write(33,*)'Z6 ',kze,'   201    qv'
       write(33,*)'Z7 ',kze,'   201    qc'
       write(33,*)'Z8 ',kze,'   201    qr'
       write(33,*)'Z9 ',kze,'   201    qi'
       write(33,*)'Z10 ',kze,'   201    qs'
       write(33,*)'Z11 ',kze,'   201    qg'
       write(33,*)'Z12 ',kze,'   201    qh'
       write(33,*)'Z13 ',kze,'   201    dbz'
       write(33,*)'Z14 ',kze,'   201    QTOT'
       write(33,*)'Z15 ',kze,'   201    theta pert'
       write(33,*)'Z16 ',kze,'   201    vort'
       write(33,*)'Z17 ',kze,'   201    p pert'
       write(33,*)'Z18 ',kze,'   201    LWC'
       write(33,*)'Z19 ',kze,'   201    thetaE WIKIPEDIA'
       write(33,*)'Z20 ',kze,'   201    thetaE Bolton'
       write(33,*)'Z21 ',kze,'   201    thetaE Emanuel'
       write(33,*)'Z22 ',kze,'   201    thetaE Emanuel using p from model'
       write(33,*)'Z23 ',kze,'   201    RH'
       write(33,*)'Z24 ',kze,'   201    rain_rat'
       write(33,*)'Z25 ',kze,'   201    rain_acc'
       write(33,*)'Z26 ',kze,'   201    hail_rat'
       write(33,*)'Z27 ',kze,'   201    hail_acc'
       write(33,*)'ENDVARS'


        allocate( gga0(1:nvars,kzb:kze,jyb:jye,ixb:ixe) )
        allocate( db(kzb:kze) )
        allocate( pb(kzb:kze) )
        
        allocate( u3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( v3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( w3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( tke(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( theta(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( thetap(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qr(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( vort(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qi(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qc(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qh(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qg(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qv(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qs(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( qtot(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( dbz3d(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( pn(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( dn(ixb:ixe,jyb:jye,kzb:kze) )

        allocate( thetae(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( temp(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( thetae2(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( thetae3(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( thetae4(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( RH(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( press2(ixb:ixe,jyb:jye,kzb:kze) )
       
        allocate( t1(ixb:ixe,jyb:jye) )
        allocate( t2(ixb:ixe,jyb:jye) )
        allocate( t3(ixb:ixe,jyb:jye) )
        allocate( t4(ixb:ixe,jyb:jye) )

        allocate( precipr(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( precipr2(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( preciph(ixb:ixe,jyb:jye,kzb:kze) )
        allocate( preciph2(ixb:ixe,jyb:jye,kzb:kze) )

        DO k = kzb,kze
          db(k) = 1.0e5*gd%var(piinit)%flt1d(k)**2.509/(287.04*gd%var(sinit+lt)%flt1d(k))
          pb(k) = 1.0e5*gd%var(piinit)%flt1d(k)**3.509
        ENDDO

        

!-----------------------------------------------------------------------------
! Main time step loop

          DO it = 1, nt, iskip

!          itstart=(10800/tarray(2)) + 1 ! because array starts at 0 ..uiiii
!          itend=(14400/tarray(2))  + 1

!          print *, 'start and end indices=', itstart, itend

!          DO it = itstart,itend,1 ! from 10800 to 14400. by setp of 1 (i.e., 240 sec)

!----------------------------------------------------------------------
! Time and location of the grid
        
              time     = tarray(it)

              write(6,*) 'it,time = ',it,time

              
              IF ( time .ge. start ) THEN

!              write(6,*) 'it loop: it,time = ',it,time

              IF ( it .ge. 1 ) CALL GRID_READ_NETCDF( gd, trim(filename), tarray(it), 1 )

        t1(ixb:ixe,jyb:jye)=gd%var(rain_rat)%flt2d(ixb:ixe,jyb:jye) ! rain_rate
        t2(ixb:ixe,jyb:jye)=gd%var(rain_acc)%flt2d(ixb:ixe,jyb:jye) ! rain_acc
        t3(ixb:ixe,jyb:jye)=gd%var(hail_rat)%flt2d(ixb:ixe,jyb:jye) ! hail_rate
        t4(ixb:ixe,jyb:jye)=gd%var(hail_acc)%flt2d(ixb:ixe,jyb:jye) ! hail_acc

        
        u3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(u)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! u-wind
        v3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(v)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! v-wind
        w3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(w)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! w-wind

        dbz3d(ixb:ixe,jyb:jye,kzb:kze)=gd%var(dbz)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! dBZ

        vort(ixb:ixe,jyb:jye,kzb:kze)=gd%var(wz)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! vertical vorticity

        theta(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lt)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! theta
        
        DO k = kzb,kze
        thetap(ixb:ixe,jyb:jye,k)=gd%var(s+lt)%flt3d(ixb:ixe,jyb:jye,k) - gd%var(sinit+lt)%flt1d(k) ! pert. theta
        ENDDO
        
        tke(ixb:ixe,jyb:jye,kzb:kze)=gd%var(km)%flt3d(ixb:ixe,jyb:jye,kzb:kze)   ! TKE (KM, actually)
        qv(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lv)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qvapor
        qc(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lc)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qcloud
        qr(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lr)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qrain
        qi(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+li)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qice
        qs(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+ls)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qsnow
        qg(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lh)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qgraupel
        IF ( lhl > 1 ) qh(ixb:ixe,jyb:jye,kzb:kze)=gd%var(s+lhl)%flt3d(ixb:ixe,jyb:jye,kzb:kze) ! qhail
        ! qtot
        qtot(ixb:ixe,jyb:jye,kzb:kze) = gd%var(s+lc)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lr)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+li)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+ls)%flt3d(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lh)%flt3d(ixb:ixe,jyb:jye,kzb:kze)
        IF ( lhl > 1 ) THEN
        qtot(ixb:ixe,jyb:jye,kzb:kze) = qtot(ixb:ixe,jyb:jye,kzb:kze) + &
                                         gd%var(s+lhl)%flt3d(ixb:ixe,jyb:jye,kzb:kze)
        ENDIF
        DO k = kzb,kze
           pn(ixb:ixe,jyb:jye,k) = 1.0e5*(gd%var(piinit)%flt1d(k)+gd%var(pi)%flt3d(ixb:ixe,jyb:jye,k))**3.509 - pb(k)  ! pert. pressure
           dn(ixb:ixe,jyb:jye,k) = 1.0e5*(gd%var(piinit)%flt1d(k)+gd%var(pi)%flt3d(ixb:ixe,jyb:jye,k))**2.509/(rd*gd%var(sinit+lt)%flt1d(k)*(1.0+0.61*gd%var(sinit+lv)%flt1d(k))) ! air density
        ENDDO

!  COMPUTER THETAE iSIMPLE or BOLTON'sFLRMULA                                                                                                      
      do k=1,kze
      do j=1,jye
      do i=1,ixe
      qvuse=qv(i,j,k)
      cpmix=3.5*287.04*(1.+0.94*qvuse)
      cvmix=717.*(1.+1.07*qvuse)
      rgmix=287.04*(1.+0.61*qvuse)
      prrcp=1.e5**(rgmix/cpmix)
      press2(i,j,k)=(dn(i,j,k)*theta(i,j,k)*rgmix/prrcp)**(cpmix/cvmix) ! substitute T from theta equ into ideal gas law                             
!      temp(i,j,k)=theta(i,j,k)/rho(i,j,k)*(press2(i,j,k)/1.e5)**(rgmix/cpmix)                                                                    
      temp(i,j,k)=press2(i,j,k)/(dn(i,j,k)*287.04) ! RHO IN MOFO g/kg and P PA.                                                                      
      thetae(i,j,k)=(temp(i,j,k)+(2.4*1.e6/1004.)*qvuse )*(1.e5/press2(i,j,k))**(287./1004.)
!      if (k.eq.1) print*,press2(i,j,k),qvuse,dn(i,j,k),temp(i,j,k),thetae(i,j,k)
      RH(i,j,k)=qv(i,j,k)/(0.622*(611.2*exp(17.67*(temp(i,j,k)-273.15)/(temp(i,j,k)-29.65))/press2(i,j,k))) 
      tlcl=55.+ 1./(1./(temp(i,j,k)-55.)-alog(RH(i,j,k))/2840.)
      thetae2(i,j,k)=theta(i,j,k)*exp((3.376/tlcl-0.00254)*1000.*qvuse*(1+0.81*qvuse))
      
      thetae3(i,j,k)=thetaee(qvuse,theta(i,j,k),press2(i,j,k))
      press2(i,j,k)=pn(i,j,k)+pb(k)
      thetae4(i,j,k)=thetaee(qvuse,theta(i,j,k),press2(i,j,k))

      precipr(i,j,k)=t1(i,j)
      precipr2(i,j,k)=t2(i,j)
      preciph(i,j,k)=t3(i,j)
      preciph2(i,j,k)=t4(i,j)

      enddo
      enddo
      enddo

      do k=1,kze
      do j=1,jye
      do i=1,ixe
      gga0(1,k,j,i)=u3d(i,j,k)!/rho(i,j,k)
      gga0(2,k,j,i)=v3d(i,j,k)
      gga0(3,k,j,i)=w3d(i,j,k)
      gga0(4,k,j,i)=theta(i,j,k)
      gga0(5,k,j,i)=tke(i,j,k)
      gga0(6,k,j,i)=qv(i,j,k)
      gga0(7,k,j,i)=qc(i,j,k)
      gga0(8,k,j,i)=qr(i,j,k)
      gga0(9,k,j,i)=qi(i,j,k)
      gga0(10,k,j,i)=qs(i,j,k)
      gga0(11,k,j,i)=qg(i,j,k)
      gga0(12,k,j,i)=qh(i,j,k)
      gga0(13,k,j,i)=dbz3d(i,j,k)
      gga0(14,k,j,i)=QTOT(i,j,k)
      gga0(15,k,j,i)=thetap(i,j,k)
      gga0(16,k,j,i)=vort(i,j,k)
      gga0(17,k,j,i)=pn(i,j,k)
      gga0(18,k,j,i)=dn(i,j,k)*(qc(i,j,k)+qr(i,j,k))
      gga0(19,k,j,i)=thetae(i,j,k)
      gga0(20,k,j,i)=thetae2(i,j,k)
!      gga0(20,k,j,i)=thetaee(qv(i,j,k),theta(i,j,k),press2(i,j,k))
      gga0(21,k,j,i)=thetae3(i,j,k)
      gga0(22,k,j,i)=thetae4(i,j,k)
      gga0(23,k,j,i)=RH(i,j,k)
      gga0(24,k,j,i)=precipr(i,j,k)
      gga0(25,k,j,i)=precipr2(i,j,k)
      gga0(26,k,j,i)=preciph(i,j,k)
      gga0(27,k,j,i)=preciph2(i,j,k)
      enddo
      enddo
      enddo

      call ga0(gga0,irec10,irec20,irec30,irec40,nvars,kze,jye,ixe)


               ENDIF
!-----------------------------------------------------------------------------
! END MAIN TIME STEP LOOP
        
       END DO 

      deallocate(u3d,v3d,theta,thetap,qv,qi,qc,qr,qs,qg,qh,qtot,dbz3d,pn,dn,w3d,vort)
      deallocate(temp,thetae,thetae2,thetae3,thetae4,RH,press2)
      deallocate(precipr,precipr2,preciph,preciph2,t1,t2,t3,t4)
      deallocate(gga0)

        STOP
        END
!-------------------
!  function to compute  equivalent potential temperature
!
!
      real function thetaee(qvc,tc,pp)
      implicit none
      
! input
      real  ::   qvc,tc,pp
      
! local
      real, parameter :: cp = 1004.
      real, parameter :: rd = 287.
      real, parameter :: cap = rd/cp
      real, parameter :: capi = 1./cap
      real, parameter :: ep = 0.622
      real, parameter :: eld = 2369.3
      real, parameter :: ar = 18.016/8314.
      real, parameter :: el = 2500300.
      
      real :: pt,tem, qvcs, rh, tt, cd
      real :: c, cc, tlcl
      integer :: nthe
!      print*,qvc,tc,pp 
      pt = pp
      tem  = tc / ((100000./pt)**cap)
      
      qvcs = 380./pt*exp(17.27*(tem-273.16)/(tem-35.86))
      rh = qvc/qvcs
      tt = tem
      cd = qvc*(pt/100.)/(6.11*(tem**3.5)*(qvc + ep))
      if (  rh .gt. .002 .or. pt .gt. 10000. ) then
       DO nthe = 1,20
        c = el - eld*(tt-273.16)
        cc = exp(ar*c*(1./273.16 - 1./tt))
        tlcl = tt - (cc - cd*tt**3.5)*tt*tt/(c*ar*cc)
        if(abs(tt-tlcl).lt..001) EXIT
        tt = tlcl
        thetaee = tc*exp(c*qvc/(cp*tlcl))
       ENDDO
      
       thetaee = tc*exp(c*qvc/(cp*tlcl))
      else
       thetaee = tc
      end if
      
      return
      end
!
!
!-------------------------------------------------------------------------------
!   ==========SUBROUTINES===========
!    NORMAL

      subroutine ga0(a,irec1,irec2,irec3,irec4,nfld,l,m,n)
      integer l,m,n,nfld
      real*4 a(nfld,l,m,n)
     open(29,file='commas.dat',form='binary',access='direct',recl=4*n*m) 
       do la=1,nfld
      do i=1,l
      itot=irec1+irec2+irec3+irec4
      write(29,rec=itot)((a(la,i,j,k),k=1,n),j=1,m)
       irec4=irec4+1
      enddo
       irec4=irec4-1
      irec3=irec3+1
      enddo
      irec3=irec3-1
      irec2=irec2+1
      return
      end

