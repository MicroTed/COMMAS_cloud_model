MODULE TICK_MODULE

 implicit none

 double precision, parameter :: time_scale = 1000.d0

 TYPE TIME_KEEPER
  integer*8   :: scale 
  integer*8   :: start
  integer*8   :: stop
  integer*8   :: time
  integer*8   :: tick
  real        :: dt
 END TYPE TIME_KEEPER

 type(time_keeper) :: clock

CONTAINS
!
 FUNCTION TICK_INIT(start, stop, time, dt) result(output)

  real start, stop, time
  real dt
  logical output

!  print *, 'tick_init:  ', start, stop, time, dt
  output = .false.

  clock%scale    = int(time_scale)
  clock%dt       = dt
  clock%start    = int(start) * clock%scale
  clock%time     = int(time)  * clock%scale
  clock%stop     = int(stop)  * clock%scale
  clock%tick     = int(time_scale * dt)
 
  output = .true.

 END FUNCTION

 FUNCTION TICK_TIME(raw) result(time)

  integer*8 time
  logical, optional :: raw

  IF( present(raw) ) THEN
   IF( raw ) THEN
    time = clock%time 
   ELSE
    time = clock%time / clock%scale
   ENDIF
  ELSE
    time = clock%time / clock%scale
  ENDIF

 END FUNCTION
!
 FUNCTION TICK_TOCK() result(output)

  logical output

  clock%time = clock%time + clock%tick

  output = .false.
  IF(clock%time .gt. clock%stop) THEN
   output = .true.
  ENDIF

 END FUNCTION
!
 FUNCTION ALARM(interval) result(output)

  logical output
  integer interval

  output = .false.
  IF( interval .le. 0 ) RETURN
  IF( mod(clock%time, Int(float(interval*clock%scale),8) ) .lt. clock%tick ) output = .true.


 END FUNCTION
END MODULE
