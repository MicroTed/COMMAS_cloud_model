!===========================================================================
!
!
!
!
!
!   /////////////////////            BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\       CPUTIME_MODULE       ////////////////////
!
!
!
!
!
!===========================================================================

MODULE CPUTIME_MODULE

 CONTAINS

!----------------------------------------------------------------------
!
!   /////////////////////        BEGIN        \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\       CLD_CPU       ////////////////////
!
!----------------------------------------------------------------------
! Routine sets up and stores off information about cpu times
! To use, simply call cld_cpu with a 'string' which will be the
! name of the cpu timing buffer.  Calling cld_cpu sets a timer up,
! calling it a second time computes the elapsed time and that amount
! is stored in a saved variable.  Example:
!
!  call cld_cpu('ADVECTION')  ! Initialize timing
!  call advectvar(t)          ! do some computing
!  call cld_cpu('ADVECTION')  ! Compute elapsed time
!
!  call cld_cpu('-1')         ! prints out the timming results 
!                             ! using a user supplied routine below.
!
!  call cld_cpu('-0')         ! if you want to clear out the timing
!                             ! buffers during a run, call this.
!----------------------------------------------------------------------
! Written by Louis Wicker, January 1998
!----------------------------------------------------------------------
      subroutine cld_cpu(string)

      implicit none
      character( LEN = * ) string

      integer maxtimers, maxtimers2
      parameter(maxtimers = 100)
      parameter(maxtimers2 = 2*maxtimers)
      real tmp(2)
      double precision :: timer, delta_t
      real time1
      real etime, second

      integer n_timer, slen, m, n
      integer count1, count_rate, count_max
      integer(8) count,addcount
      integer oldcount
      
      character( LEN = 20 ) timer_name(maxtimers), name
      logical timer_log(maxtimers)
      double precision :: timer_time(maxtimers,2)
      double precision :: the_time
      integer(8) itimer_time(maxtimers,2), ithe_time, idelta_t
      integer iunit

      integer :: values1(8)
      integer :: values2(8)
      
      logical debug
      parameter ( debug = .false. )

      data timer_name / maxtimers*'                    ' /
      data timer_log  / maxtimers*.false. /
      data timer_time / maxtimers2*0.0d0 /
      data itimer_time / maxtimers2*0 /
      data n_timer    / 0 /
      data addcount   / 0 /
      data oldcount   / -1 /

! Save timing info between invocations of subroutine

      save timer_name, timer_log, timer_time, n_timer, itimer_time
      save count, addcount, oldcount, values2

!-----------------------------------------------------------------------
! Set here the timing calls for various machines
!-----------------------------------------------------------------------

      call mysystem_clock(count1, count_rate, count_max,values1)
      IF ( oldcount .lt. 0 ) THEN 
        oldcount = count1
        count = 0
      ENDIF
      IF ( oldcount .le. count1 ) THEN
        count = count + ( count1 - oldcount )
        oldcount = count1
      ELSEIF ( oldcount .gt. count1 ) THEN ! clock has restarted at midnight
        IF ( values1(3) .eq. values2(3) ) THEN ! same day of month still -- something screwy
!        write(*,*) 'cputime: Clock has restarted? oldcount,count1 = ',oldcount,count1,count
!        write(*,*) 'old values: ',values2
!        write(*,*) 'new values: ',values1
        count1 = oldcount
        ELSE
        addcount = addcount + count_max
        count = count + (count_max - oldcount + count1)
        oldcount = count1
        ENDIF
      ENDIF
      
      values2(:) = values1(:)
      timer = float(count) / float(count_rate)


!-----------------------------------------------------------------------
! Check to see if the user wants the buffer flushed
!-----------------------------------------------------------------------
      IF( string(1:2) .eq. '-0' ) THEN
            n_timer = 0
            do n = 1,maxtimers
                  timer_name(n) = '                    '
                  timer_log(n)  = .false.
                  timer_time(n,1) = 0.d0
                  timer_time(n,2) = 0.d0
                  itimer_time(n,1) = 0
                  itimer_time(n,2) = 0
            enddo
            RETURN
      ENDIF
!-----------------------------------------------------------------------
! If not, then get on with it
!-----------------------------------------------------------------------
      m = 0
      slen = len(string)

      the_time = timer
      ithe_time = count

     IF ( debug )  print *, 'cpu_time: ', string, the_time, count_rate, count
            
      if( string(1:2) .ne. '-1' .and. string(1:2) .ne. '-2' .and. string(1:2) .ne. '-3' ) then

! Check to see if we have any timers set and if so, find what bin its in.

      if( n_timer .gt. 0 ) then
            do n = 1,n_timer
             name = timer_name(n)
             if(name .eq. string(1:slen)) m = n
            enddo
      endif

! If nothing matches, create a new name and set things up

      IF( m .eq. 0 ) THEN
      
! Error check, make sure we have enough space...

            if( n_timer .eq. maxtimers ) then
             write(0,*) 'CLD_CPU:  n_timer = maxtimers, returning'
             return
            endif
            n_timer = n_timer+1
            name = '                    '
            name(1:slen)          = string
            timer_name(n_timer)   = name
            timer_log(n_timer)    = .true.
            timer_time(n_timer,1) = the_time
            itimer_time(n_timer,1) = ithe_time
            IF ( debug )  print *, 'set up new timer name = ',string,n_timer,timer_name(n_timer)
            return
            
       ENDIF

! IF M > 0, then go ahead and reset stuff

   if( .not. timer_log(m)) then
            timer_time(m,1) = the_time
            itimer_time(m,1) = ithe_time
            timer_log(m)    = .true.
       else
!            IF( the_time .ge. timer_time(m,1) ) THEN
!              delta_t = the_time - timer_time(m,1)
!            ELSE
!              the_time = float(count+count_max) / float(count_rate)
!              delta_t = the_time - timer_time(m,1)
!              ithe_time = count + count_max
!            ENDIF
              delta_t = the_time - timer_time(m,1)
              idelta_t = ithe_time - itimer_time(m,1)
              delta_t = idelta_t/float(count_rate)
              IF ( delta_t .lt. 0 ) THEN
               write(0,*) 'cpu_time, negative delta_t: ',delta_t, idelta_t,count,count_max
               delta_t = float(count+count_max) / float(count_rate) - timer_time(m,1)
               
              ENDIF
            IF ( itimer_time(m,1) .eq. ithe_time ) delta_t = 0.0
          IF ( debug )  print *, 'cpu_time: ', string, ' delta_t = ', delta_t, timer_time(m,1), the_time
!           IF( delta_t .lt. 0.0 ) delta_t = delta_t + float(count_max)/float(count_rate)
            timer_time(m,2) = timer_time(m,2) + delta_t
            itimer_time(m,2) = itimer_time(m,2) + count - itimer_time(m,1)
            timer_log(m)    = .false.
            timer_time(m,1) = the_time
            itimer_time(m,1) = ithe_time
!          IF ( debug )  print *, 'cpu_time: ', string, ' delta_t = ', delta_t
       endif

      ELSE

! User supplied routine for print out timings
        if ( string(1:2) .eq. '-1' ) iunit = 1
        if ( string(1:2) .eq. '-2' ) iunit = 6
        if ( string(1:2) .eq. '-3' ) iunit = 12

        call cld_cpu_print(timer_name, timer_time(1,2), n_timer, iunit)

!           call cld_cpu_print(timer_name, timer_time(1,2), n_timer)

      ENDIF

 RETURN
 END SUBROUTINE CLD_CPU

!----------------------------------------------------------------------
!
!   /////////////////////        BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    CLD_CPU_PRINT     ////////////////////
!
!----------------------------------------------------------------------
! This routine is called whenever cld_cpu is called with the string
! '-1'.  Here, the user can supply a method for printing out all the
! timing statistics stored in cld_cpu.
!
! The user defined 'timing buffers' in the variable 'names' provides 
! a way for the user to format the output.  Any additional timing
! output buffers (such as something put in for debugging) will be
! automatically dumped out at the end at the end.
!----------------------------------------------------------------------
! Written by Louis Wicker, January 1998
!----------------------------------------------------------------------
      subroutine cld_cpu_print(timer_name, timer_time, n_timer, iunit)
      USE PARAM_MODULE
      implicit none
      integer n_timer, n, m, mm, mm1, maxtimers
      character( LEN = 20 ) timer_name(n_timer), name1, name2
      double precision timer_time(n_timer)
      real percentage, total, interface
      
      integer iunit

! Parameters so we can sort names...

      parameter( mm = 20, mm1 = mm+1, maxtimers=100 )
      character( LEN = 20 ) :: names(mm)
      integer list(mm1), list2(maxtimers)
      data list      / mm1*-1 /
      data list2     / maxtimers*-1 /

! This is a list of timing buffers the model uses -> and we want to format
! the output for.  Other buffers that are added will still be dumped at end.

      data names(1)  / 'SOLVER             ' /
      data names(2)  / 'SMLSTEP            ' /
      data names(3)  / 'ADVECT UVW         ' /
      data names(4)  / 'ADVECT SCALAR      ' /  
      data names(5)  / 'MIX                ' /
      data names(6)  / 'FILTER             ' /
      data names(7)  / 'MICROPHYSICS       ' /
      data names(8)  / 'TIMESTEP           ' /
      data names(9)  / 'TKE                ' /
      data names(10) / 'HOLE-FILLING       ' /
      data names(11) / 'MPICOM-SOLVER      ' /
      data names(12) / 'ELECTRICITY        ' /
      data names(13) / 'INTERPOLATE        ' /
      data names(14) / 'ENKF               ' /
      data names(15) / 'CORRECT            ' /
      data names(16) / 'SFC_PHYSICS        ' /

      IF ( iunit .eq. 1 .or. iunit == 6) THEN !{
      
        IF (iunit == 1 ) iunit = luno
      
      write(iunit,999)
      write(iunit,99) 'COMMAS MODEL CPU STATISTICS'
      write(iunit,999)

! Do a sort to find the order in which we want to print....

      do m = 1,mm
       name1 = names(m)
       do n = 1,n_timer
        name2 = timer_name(n)
        if( name1 .eq. name2 ) then
          list(m)  = n
          list2(n) = n
        endif
       enddo
      enddo

      IF ( list(1) .lt. 1 ) THEN
       OUTER: do m = 1,mm
       name1 = names(m)
       do n = 1,n_timer
        name2 = timer_name(n)
        if( name1 .eq. name2 ) then
          list(1)  = n
          EXIT OUTER
        endif
       enddo
      enddo OUTER
      
      ENDIF
! Now write timing for simulation out

      write(iunit,100) timer_name(list(1))//' is ',timer_time(list(1))
      write(iunit,999)

! Check to see if this is an adaptive run, if so, print out that stuff
!     if( list(4) .ne. -1 ) then
!       interface = timer_time(list(1)) 
!     $           - timer_time(list(2)) 
!     $           - timer_time(list(3))
!       write(6,100) 'INTERFACE            is ',interface
!       write(6,999)
!      DO m = 2,mm
!        n = list(m)
!        write(6,100) timer_name(n)//' is ',timer_time(n)
!      ENDDO
!       write(6,999)
!     endif

! Now finish with the user output and solver cpu stats.
!     write(6,100) timer_name(list(2))//' is ',timer_time(list(2))
!     write(6,999)
!     write(6,101) timer_name(list(3)),timer_time(list(3)), 100.
!     write(6,999)

! Now write out the solver statistics
      total = 0.0
      DO m = 2,mm
       n = list(m)
       IF ( n .ge. 1 ) THEN
       IF( timer_time(n) .ge. 0.1 ) THEN
         percentage = timer_time(n)/timer_time(list(1))*100.
         write(iunit,101) timer_name(n),timer_time(n),percentage
       ELSE
         percentage = 0.0
         write(iunit,101) timer_name(n),timer_time(n),percentage
       ENDIF
       total = total + timer_time(n)
       ENDIF
      ENDDO
      write(iunit,999)

! Now dump out any additional buffers that are left. 
        do n = 1,n_timer
         if( list2(n) .eq. -1 ) write(iunit,100) timer_name(n),timer_time(n)
       enddo

! Do an error check at the end to see if the solver cpu timings add up.
      write(iunit,999)
      write(iunit,102) timer_time(list(1)) - total
      write(iunit,999)
      
      ELSEIF (iunit == 12) THEN !}{
      ! just print out total time

        iunit = luno
        do n = 1,n_timer
         IF ( timer_name(n) == 'MAIN' ) THEN
         write(iunit,'(a,f15.5)') 'TOTAL SIMULATION TIME = ',timer_time(n)
         exit
         ENDIF
        enddo
      
      
      ENDIF

99    format(1x,a)
100   format(1x,'Total CPU for ',a,f15.5)
101   format(1x,'Total CPU for ',a,4x,f15.5,5x,f6.2,' %')
102   format(1x, '(SOLVER CPU) - (SUM OF SOLVER ROUTINES CPU) = ',f15.5)
999   format(71('-'))
      
       
      RETURN
      END SUBROUTINE CLD_CPU_PRINT

      
END MODULE CPUTIME_MODULE
