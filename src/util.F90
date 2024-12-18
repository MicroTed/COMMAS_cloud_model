 SUBROUTINE wrf_debug( level, message )
   implicit none
   integer :: level
   character(*) :: message
   
   IF ( level < 0 ) THEN
     write(0,*) message
   ENDIF
   
 END SUBROUTINE wrf_debug

!
! #####################################################################
!
 SUBROUTINE wrf_message( message )
   implicit none
   character(*) :: message
   
     write(0,*) message
   
 END SUBROUTINE wrf_message

!
! #####################################################################
!
 SUBROUTINE wrf_error_fatal( message )
    USE COMMASMPI_MODULE, only: commasmpi_abort
   implicit none
   character(*) :: message
   
     write(0,*) message
     call commasmpi_abort()
   
 END SUBROUTINE wrf_error_fatal

!
! #####################################################################
!

 LOGICAL FUNCTION wrf_dm_on_monitor()
   USE COMMASMPI_MODULE, only: my_rank
   wrf_dm_on_monitor = my_rank == 0
 END FUNCTION wrf_dm_on_monitor

!
! #####################################################################
!

 SUBROUTINE STRING_LIMITS( string, ibeg, iend ) 

  IMPLICIT NONE 

! List of dummy arguments

  CHARACTER(len=*),INTENT(IN) :: string ! Input string 
  INTEGER,INTENT(OUT)         :: ibeg   ! First non-blank character 
  INTEGER,INTENT(OUT)         :: iend   ! Last non-blank character 

! List of local variables: 

  INTEGER :: length            ! Length of input string 

! Get the length of the input string. 

  length = LEN ( string ) 

! Look for first character used in string.  Use a WHILE loop to find the first non-blank character. 

  ibeg = 0 

  DO 
   ibeg = ibeg + 1 
   IF ( ibeg > length ) EXIT
   IF ( string(ibeg:ibeg) /= ' ' ) EXIT 
  END DO 

! If ibeg > length, the whole string was blank.  Set 
! ibeg = iend = 1.  Otherwise, find the last non-blank 
! character. 

  IF ( ibeg > length ) THEN 
   ibeg = 1 
   iend = 1 
  ELSE 
! Find last nonblank character. 
   iend = length + 1 
   DO 
      iend = iend - 1 
      IF ( string(iend:iend) /= ' ' .and. ICHAR(string(iend:iend)) /= 0 ) EXIT 
   END DO 
  END IF 

 END ! SUBROUTINE STRING_LIMITS



!
!
! #####################################################################
!


         SUBROUTINE GREGJUL(MGREG,JULIAN)
!...CONVERTS GREGORIAN DATE TO JULIAN DATE
!...MGREG IS IN FORM MMDDYYYY
!...JULIAN IS RETURNED IN FORM YYYYJJJ
!...Y2K-ified by tm on 6/13/2000 (4-digit dates)
      INTEGER NUMDAY(12)
      DATA  NUMDAY /31,28,31,30,31,30,31,31,30,31,30,31/
      JULIAN=0
      MONTH=MGREG/1000000
      MDAY=(MGREG-(MONTH*1000000))/10000
      MYEAR=MGREG-((MONTH*1000000)+(MDAY*10000))
!...CHECK FOR LEAP YEAR
! y2000 is also a leap year (divisible by 400), so don't need to worry until 2100.
      IF(MOD(MYEAR,4).EQ.0) NUMDAY(2)=29 
!...ACCUMALATE DAYS TO DATE
      MONTH=MONTH-1
      IF(MONTH.EQ.0) GO TO 40
      DO 30 LP=1,MONTH
      JULIAN=JULIAN+NUMDAY(LP)
30    CONTINUE
40    JULIAN=JULIAN+MDAY
      MYEAR=MYEAR*1000
      JULIAN=JULIAN+MYEAR
      NUMDAY(2)=28
        RETURN
        END
        
!
! #####################################################################
!
!
! #####################################################################
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  STRAKAS ATMOSPHERIC MODEL  (SAM)
!    Designed by Jerry M. Straka
!
!  RNDNUM(ISEED) (must define iseed=-1 on intialization 
!                 in user routine)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!  start function for random numbers....
!
      function rndnum(iseed)
!
!     implicit none
!
      real rndnum
!
      real    ambig
      parameter (ambig=4000000)
      real    amseed
      parameter (amseed=1618033)
      real    amz,amj
      parameter (amz=0)
!
      real    ama(55)
!
      real       fac
      parameter (fac=1./ambig)
!
      save iff, inext, inextp, ama
!
      data iff /0/
!
      if ( iseed .lt. 0 .or. iff .eq. 0 ) then
      iff = 1
      amj = amseed-iabs(iseed)
      amj = mod(amj,ambig)
      ama(55) = amj
      amk = 1
      do n = 1,54
      nn = mod(21*n,55)
      ama(nn) = amk
      amk = amj-amk
      if (amk.lt.amz) amk=amk+ambig
      amj=ama(nn)
      end do
!
      do m = 1,4
      do n = 1,55
      ama(n)=ama(n)-ama(1+mod(n+30,55))
      if (ama(n).lt.amz) ama(n)=ama(n)+ambig
      end do
      end do
!
      inext = 0
      inextp = 31
      iseed = 1
!
      end if
!
! find random number
!
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      amj=ama(inext)-ama(inextp)
      if (amj.lt.amz) amj=amj+ambig
      ama(inext)=amj
      rndnum=amj*fac
!
      return
      end
!
!
!
!
!  end of random number function
!
! 
!
! #####################################################################
      SUBROUTINE SRAND(ISEED)
!
!  This subroutine sets the integer seed to be used with the
!  companion RAND function to the value of ISEED.  A flag is 
!  set to indicate that the sequence of pseudo-random numbers 
!  for the specified seed should start from the beginning.
!
      implicit none
      integer jseed,ifrst,iseed
      COMMON /SEED/JSEED,IFRST
!
      JSEED = ISEED
      IFRST = 0
!
      RETURN
      END
      
! #####################################################################
      REAL FUNCTION RAND()
!
!  This function returns a pseudo-random number for each invocation.
!  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal 
!  standard number generator whose Pascal code appears in the article:
!
!     Park, Steven K. and Miller, Keith W., "Random Number Generators: 
!     Good Ones are Hard to Find", Communications of the ACM, 
!     October, 1988.
!
      implicit none
      integer MPLIER,MODLUS,MOBYMP,MOMDMP
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,MOMDMP=2836)
!
      integer jseed,ifrst
      COMMON  /SEED/JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN
!
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF
!
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAND = REAL(NEXTN)/REAL(MODLUS)
!
      RETURN
      END
! #####################################################################
      BLOCKDATA RANDBD
      integer jseed,ifrst
      COMMON /SEED/JSEED,IFRST
!
      DATA JSEED,IFRST/123456789,0/
!
      END

! #####################################################################
      real function commas_dtime(timarr)
      implicit none
      integer count, count_rate,countmax
      real lasttime,time1,time
      integer :: values(8)
      real timarr(2)
      data lasttime/0.0/
      save lasttime

      call mysystem_clock(count, count_rate, countmax,values)
      time1 = float(count) / float(count_rate)
      
      IF ( time1 .ge. lasttime )  THEN
         time = time1 - lasttime
      ELSE
         time = float(count+countmax) / float(count_rate) - lasttime
      ENDIF
      lasttime = time1
      timarr(:) = time
      commas_dtime = time
      
      RETURN
      END

! #####################################################################
      real function commas_dtime2(lasttime,iop)
      implicit none
      integer :: iop ! = 0 to set lasttime, = 1 to get difference
      integer count, count_rate,countmax
      real lasttime,time1,time
      integer :: values(8)
!      real timarr(2)
!      data lasttime/0.0/
!      save lasttime

      call mysystem_clock(count, count_rate, countmax,values)
      time1 = float(count) / float(count_rate)
      
      IF ( iop == 1 ) THEN
        IF ( time1 .ge. lasttime )  THEN
           time = time1 - lasttime
        ELSE
           time = float(count+countmax) / float(count_rate) - lasttime
        ENDIF
      ELSE
         time = time1
      ENDIF
      
      lasttime = time1
      commas_dtime2 = time
      
      RETURN
      END

! ################################################################
      SUBROUTINE mysystem_clock(count, count_rate, count_max, outvalues)
      implicit none
      
      integer :: values(8)
      integer :: outvalues(8)
      character(LEN=8)  :: date
      character(LEN=10) :: ctime,zone
      integer count,count1, count_rate, count_max
      
            
      count_max = 86400000 - 1 ! 1000(ms/s)*60(s/min)*60(min/hr)*24(hr)
      count_rate = 1000
      
      Call DATE_AND_TIME(date,ctime,zone,values)
!      print*, 'i,values = ',i,(values(k),k=1,8)
      count = 1000*((60*(values(5)*60 + values(6))) + values(7)) + values(8)
      
        outvalues(:) = values(:)
      
      RETURN
      END 

! ################################################################
! 
!
!  function to compute  equivalent potential temperature
!
!
      real function thetae(qvc,tc,pp)
      implicit none
      
     ! input
      real    :: qvc,tc,pp
      
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
        thetae = tc*exp(c*qvc/(cp*tlcl))
       ENDDO
      
       thetae = tc*exp(c*qvc/(cp*tlcl))
      else
       thetae = tc
      end if
      
      return
      end
!
!
!  end of function theta-e
!
!

