! v5df.h

! Include file for using v5d functions from FORTRAN programs


! Function prototypes.  See the README file for details.  These are
! the functions you'll want to use for writing v5d file converters.

!      integer v5dcreate

!      integer v5dcreatesimple

!      integer v5dwrite

!      integer v5dwriteappend

!      integer v5dmcfile

!      integer v5dclose
      
!      integer myv5dupdate


! 5-D grid limits, must match those in v5d.h!!!
      integer MAXVARS, MAXTIMES, MAXROWS, MAXCOLUMNS, MAXLEVELS
      integer MAXPROJARGS

      parameter (MAXVARS=200)
      parameter (MAXTIMES=400)
      parameter (MAXROWS=800)
      parameter (MAXCOLUMNS=4000)
      parameter (MAXLEVELS=400)
      parameter (MAXPROJARGS = MAXROWS+MAXCOLUMNS+1)

! Missing values
      real MISSING
      integer IMISSING

      parameter (MISSING=1.0E35)
      parameter (IMISSING=-987654)

      integer maxnl





