      program readball

c
c  Alpha compile:
c
c  f90 -o x.readsound -g -convert big_endian readsound.f
c
c  Reads a binary model sounding file and converts to ascii
c
c Takes the name(s) of the sounding file(s) from the command line and prints
c  the ascii to the same file name with '.txt' appended.
c
      
      implicit none
      
c
c sounding
c
      integer llx,lly,llz,llu,llv,llw
      integer llex,lley,llez,llemag,llphi
      integer llsc  ! charge density
      integer llp   ! pressure
      integer llt   ! temp (celsius)
      integer llpt  ! potential temperature (K)
      integer llqv  ! vapor mixing ratio
      parameter (llx=1,lly=2,llz=3,llu=4,llv=5,llw=6)
      parameter (llex=7,lley=8,llez=9,llemag=10,llphi=11)
      parameter (llsc=12,llp=13,llt=14,llpt=15)
      parameter (llqv=15)
      integer lltim
      parameter (lltim=16)
      integer maxsamp
      parameter (maxsamp =  5000)
      integer numvar
      parameter (numvar = 16) 
      real sound(numvar,maxsamp)
      real eh(maxsamp)
      integer indx,iargc
      real brise
      real emax,emin
      
      integer kz,n

      integer ilev,i,j
      
      character*80 infile,outfile

      real qv,p,x,tx,rh,rh2,qsat,tc
      
      real hum
      
      integer, parameter :: ndebug = 0
      integer  :: istat

      
      indx = iargc( )
      
      DO n=1,indx
      call getarg ( n, infile )
      emax = 0.0
      emin = 0.0
      
      kz = index(infile,' ')
      write(6,*) 'kz,infile = ',kz,infile
      write(outfile,'(a)') infile(1:kz-1)//'.txt'

c      read(5,*) infile

      open(unit=10,file=infile,status='old',
     :     form='unformatted')
      
      open(unit=12,file=outfile,status='unknown',
     :     form='formatted')
      
      write(6,*) 'read ilev'
      read(10,iostat=istat) ilev, brise
        IF ( istat .ne. 0 ) THEN
          write(6,*) 'reading old format without brise'
          rewind(10)
          read(10,iostat=istat) ilev
          brise = 5.0
        ENDIF
      write(6,*) 'ilev, brise = ',ilev,brise
      CALL readsound(numvar,ilev,10,sound)
      close(10)
      
      DO j=1,ilev
       sound(llx,j) = sound(llx,j)/1000. ! convert to km
       sound(lly,j) = sound(lly,j)/1000. ! convert to km
       sound(llz,j) = sound(llz,j)/1000. ! convert to km
       sound(llex,j) = sound(llex,j)/1000.  ! convert to kV/m
       sound(lley,j) = sound(lley,j)/1000.  ! convert to kV/m
       eh(j) = Sqrt( sound(llex,j)**2 + sound(lley,j)**2 )
       sound(llez,j) = sound(llez,j)/1000.  ! convert to kV/m
       sound(llemag,j) = sound(llemag,j)/1000. ! convert to kV/m
       sound(llt,j) = sound(llt,j) - 273.1    ! convert to Celsius
       sound(llsc,j) = sound(llsc,j)*1.0e10   ! scale to 0.1nC/m**3
       sound(llphi,j) = sound(llphi,j)*1.0e-6 ! convert from V to MV
       sound(llqv,j) = sound(llqv,j)*1000.0   ! convert from kg/kg to g/kg
       emax = Max(emax, sound(llez,j) )
       emin = Min(emin, sound(llez,j) )
      END DO
      
      write(6,*) 'emin,emax (kV/m) = ',emin,emax
      
      write(12,'(5a)') 
     : 'x (km), y (km), Altitude (km), u (m/s), v (m/s), w (m/s), ',
     : 'Ex (kV/m), Ey (kV/m), Ez (kV/m), Emag (kV/m), ',
     :  'Potential (MV), Charge Density (10^-10 C), P (Pa), ',
     :  'Temperature (C), qv (g/kg), ',
     :  'Time (s), Eh (kV/m), Td (C), RH (%), Rise Rate (m/s), RH2 (%)'

      DO i = 1,ilev
        qv = sound(llqv,i)/1000.
        p = sound(llp,i)
            x = Log(Max(1.0e-6,qv)* (p/380.0))
            tx = (35.88*x - 17.27*273.15)/(x - 17.27) - 273.15
            tc = real(sound(llt,i))
            tx = Min(tx,tc)
            rh = hum(tc,tx)
c   Different rh calculation:   
      tc = sound(llt,i) + 273.1
      qsat = (380./p)*exp(17.2693882*(tc-273.16)/(tc-35.86))
      rh2 = Min(100., 100.*qv/qsat)
      IF ( ndebug .ge. 1 ) THEN
      write(6,*) 'T, TD, rh1, rh2 = ',sound(llt,i),tx,rh,rh2
c        write(12,'(1x,19(e12.5,1x,','))') 
      ENDIF
        write(12,101) 
     :    (sound(j,i),j=1,16),eh(i),tx,rh,sound(llw,i)+brise,rh2
      END DO
c 101  format (18(e12.5,', '),e12.5)
 101  format (20(1pe12.5,', '),1pe12.5)
      CLOSE(12)
      
      ENDDO
      
      STOP
      END

c #####################################################################      
      SUBROUTINE READSOUND(numvar,nlev,ifile,sound)
c #####################################################################      

      implicit none

      integer nlev,ifile,numvar
      real sound(numvar,nlev)
      
      read(ifile) sound
      
      RETURN
      END

c #####################################################################      
        FUNCTION HUM(T,TD)
c #####################################################################      

C       INCLUDE 'LIB_DEV:[GUDOC]EDFVAXBOX.FOR/LIST'
C       G.S. Stipanuk     1973            Original version.
C       Reference Stipanuk paper entitled:
C            "ALGORITHMS FOR GENERATING A SKEW-T, LOG P
C            DIAGRAM AND COMPUTING SELECTED METEOROLOGICAL
C            QUANTITIES."
C            ATMOSPHERIC SCIENCES LABORATORY
C            U.S. ARMY ELECTRONICS COMMAND
C            WHITE SANDS MISSILE RANGE, NEW MEXICO 88002
C            33 PAGES
C       Baker, Schlatter  17-MAY-1982    

C   THIS FUNCTION RETURNS RELATIVE HUMIDITY (%) GIVEN THE
C   TEMPERATURE T AND DEW POINT TD (CELSIUS).  AS CALCULATED HERE,
C   RELATIVE HUMIDITY IS THE RATIO OF THE ACTUAL VAPOR PRESSURE TO
C   THE SATURATION VAPOR PRESSURE.

        implicit none
        
        real hum,esat,t,td
        
        HUM= 100.*(ESAT(TD)/ESAT(T))
        RETURN
        END

c #####################################################################      
        FUNCTION ESAT(T)
c #####################################################################      

C       INCLUDE 'LIB_DEV:[GUDOC]EDFVAXBOX.FOR/LIST'
C       G.S. Stipanuk     1973            Original version.
C       Reference Stipanuk paper entitled:
C            "ALGORITHMS FOR GENERATING A SKEW-T, LOG P
C            DIAGRAM AND COMPUTING SELECTED METEOROLOGICAL
C            QUANTITIES."
C            ATMOSPHERIC SCIENCES LABORATORY
C            U.S. ARMY ELECTRONICS COMMAND
C            WHITE SANDS MISSILE RANGE, NEW MEXICO 88002
C            33 PAGES
C       Baker, Schlatter  17-MAY-1982    

C   THIS FUNCTION RETURNS THE SATURATION VAPOR PRESSURE OVER
C   WATER (MB) GIVEN THE TEMPERATURE (CELSIUS).
C   THE ALGORITHM IS DUE TO NORDQUIST, W.S.,1973: "NUMERICAL APPROXIMA-
C   TIONS OF SELECTED METEORLOLGICAL PARAMETERS FOR CLOUD PHYSICS PROB-
C   LEMS," ECOM-5475, ATMOSPHERIC SCIENCES LABORATORY, U.S. ARMY
C   ELECTRONICS COMMAND, WHITE SANDS MISSILE RANGE, NEW MEXICO 88002.

        implicit none
        real t,tk,p1,p2,c1,esat
        
        TK = T+273.15
        P1 = 11.344-0.0303998*TK
        P2 = 3.49149-1302.8844/TK
        C1 = 23.832241-5.02808*ALOG10(TK)
        ESAT = 10.**(C1-1.3816E-7*10.**P1+8.1328E-3*10.**P2-2949.076/TK)
        RETURN
        END
      
