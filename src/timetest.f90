      program timetest
      
      implicit none
      
      integer count1, count_rate,count_max
      integer count2a, count2b, count_rate2,count_max2
      real lasttime,time1
      integer i,j,k
      integer(8) count,addcount
      integer oldcount
      integer :: newvalues(8), oldvalues(8)
      real x,y
      
       call mysystem_clock(count1, count_rate, count_max, newvalues)
       call system_clock(count2a, count_rate2,count_max2)
       print*, 'system_clock count_rate,max = ',count_rate2,count_max2
       oldvalues(:) = newvalues(:)

      oldcount = -1
      addcount = 0
      DO while ( oldvalues(3) == newvalues(3) )
      
       oldvalues(:) = newvalues(:)
       count2a = count2b
       call mysystem_clock(count1, count_rate, count_max, newvalues)
       call system_clock(count2b, count_rate2,count_max2)
       
      IF ( oldcount .lt. 0 ) THEN 
        oldcount = count1
        count = 0
      ENDIF
      IF ( count2b .lt. count2a ) THEN
        write(*,*) 'system_clock went down count2a,count2b,diff = ',count2a,count2b,count2b-count2a
!        EXIT
      ENDIF
      IF ( oldcount .le. count1 ) THEN
        count = count + ( count1 - oldcount )
        oldcount = count1
      ELSEIF ( oldcount .gt. count1 ) THEN ! clock has restarted at midnight
        write(*,*) 'cputime: Clock has restarted? oldcount,count1 = ',oldcount,count1,count
        write(*,*) 'system_clock gave count2a,count2b,diff = ',count2a,count2b,count2b-count2a
        IF ( oldvalues(3) == newvalues(3) ) EXIT
        addcount = addcount + count_max
        count = count + (count_max - oldcount + count1)
        oldcount = count1
      ENDIF
       
       x = i*(float(i))**0.1
      
      
      ENDDO
      
      
      END


! ################################################################
      SUBROUTINE mysystem_clock(count, count_rate, count_max,newvalues)
      implicit none
      
      integer :: values(8), newvalues(8)
      character(LEN=8)  :: date
      character(LEN=10) :: ctime,zone
      integer count,count1, count_rate, count_max
      integer i,k
      
            
      count_max = 86400000 - 1 ! 1000(ms/s)*60(s/min)*60(min/hr)*24(hr)
      count_rate = 1000
      
      Call DATE_AND_TIME(date,ctime,zone,values)
!      print*, 'values = ',(values(k),k=1,8)
      count = 1000*((60*(values(5)*60 + values(6))) + values(7)) + values(8)
      
      newvalues(:) = values(:)
      
      RETURN
      END 

